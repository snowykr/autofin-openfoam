from __future__ import annotations

import argparse
import math
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class GeometryParams:
    gap: float
    T: float
    H_wall: float
    W_wall: float
    p: float
    t: float
    L: float


@dataclass(frozen=True)
class MeshingParams:
    N_T: int = 5
    N_gapA: int = 3
    N_over: int = 10
    N_t: int = 5
    N_p: int = 5


def _validate_inputs(geom: GeometryParams, mesh: MeshingParams) -> None:
    if geom.gap <= 0 or geom.T <= 0 or geom.H_wall <= 0 or geom.W_wall <= 0:
        raise ValueError("gap, T, H_wall, W_wall must be > 0")
    if geom.t <= 0:
        raise ValueError("t must be > 0")
    if geom.p <= 0:
        raise ValueError(
            "p must be > 0 for reference-style layered blockMesh generation"
        )
    if not math.isfinite(geom.L):
        raise ValueError("L must be finite")
    if geom.L <= 0:
        raise ValueError("L must be > 0")
    if geom.L > geom.gap:
        raise ValueError("require L <= gap to avoid penetrating opposite wall")

    for name, v in (
        ("N_T", mesh.N_T),
        ("N_gapA", mesh.N_gapA),
        ("N_over", mesh.N_over),
        ("N_t", mesh.N_t),
        ("N_p", mesh.N_p),
    ):
        if v <= 0:
            raise ValueError(f"{name} must be >= 1")


def _grid_counts(geom: GeometryParams) -> tuple[int, int, float, float, float]:
    pitch = geom.t + geom.p
    if pitch <= 0:
        raise ValueError("t + p must be > 0")
    # Deterministic nearest-integer rounding (half-up), avoids bankers rounding.
    ny_units = max(2, int(math.floor(geom.H_wall / pitch + 0.5)))
    nz_units = max(2, int(math.floor(geom.W_wall / pitch + 0.5)))
    h_real = ny_units * pitch
    w_real = nz_units * pitch
    return ny_units, nz_units, h_real, w_real, pitch


def generate_blockmeshdict_text(geom: GeometryParams, mesh: MeshingParams) -> str:
    _validate_inputs(geom, mesh)

    ny_units, nz_units, h_real, w_real, pitch = _grid_counts(geom)

    left_outer = -geom.gap / 2 - geom.T
    left_inner = -geom.gap / 2
    right_inner = geom.gap / 2
    right_outer = geom.gap / 2 + geom.T

    right_fin_tip = geom.gap / 2 - geom.L
    left_fin_tip = -geom.gap / 2 + geom.L

    x2 = min(right_fin_tip, left_fin_tip)
    x3 = max(right_fin_tip, left_fin_tip)
    x_coords = [left_outer, left_inner, x2, x3, right_inner, right_outer]

    for i in range(len(x_coords) - 1):
        if not (x_coords[i] < x_coords[i + 1]):
            raise ValueError(
                "Invalid x-layer ordering (zero/negative width layer). "
                "Check L so that 0 < L < gap and avoid exactly L=gap/2."
            )

    y_coords = [0.0]
    for _ in range(ny_units):
        last = y_coords[-1]
        y_coords.append(last + geom.t)
        y_coords.append(last + pitch)

    z_coords = [0.0]
    for _ in range(nz_units):
        last = z_coords[-1]
        z_coords.append(last + geom.t)
        z_coords.append(last + pitch)

    n_x = len(x_coords)
    n_y = len(y_coords)
    n_z = len(z_coords)

    def v_idx(ix: int, iy: int, iz: int) -> int:
        return ix + (n_x * iy) + (n_x * n_y * iz)

    vertices: list[str] = []
    for iz in range(n_z):
        for iy in range(n_y):
            for ix in range(n_x):
                vertices.append(
                    f"    ({x_coords[ix]:.9f} {y_coords[iy]:.9f} {z_coords[iz]:.9f})"
                )

    nx_layer = [mesh.N_T, mesh.N_gapA, mesh.N_over, mesh.N_gapA, mesh.N_T]
    has_overlap = geom.L > geom.gap / 2

    blocks: list[str] = []
    for iz in range(n_z - 1):
        for iy in range(n_y - 1):
            is_y_fin = iy % 2 == 0
            is_z_fin = iz % 2 == 0

            ny_mesh = mesh.N_t if is_y_fin else mesh.N_p
            nz_mesh = mesh.N_t if is_z_fin else mesh.N_p

            is_fin_site = is_y_fin and is_z_fin
            grid_y = iy // 2
            grid_z = iz // 2

            is_left_fin = is_fin_site and ((grid_y + grid_z) % 2 == 0)
            is_right_fin = is_fin_site and not is_left_fin

            for ix_layer in range(5):
                v0 = v_idx(ix_layer, iy, iz)
                v1 = v_idx(ix_layer + 1, iy, iz)
                v2 = v_idx(ix_layer + 1, iy + 1, iz)
                v3 = v_idx(ix_layer, iy + 1, iz)
                v4 = v_idx(ix_layer, iy, iz + 1)
                v5 = v_idx(ix_layer + 1, iy, iz + 1)
                v6 = v_idx(ix_layer + 1, iy + 1, iz + 1)
                v7 = v_idx(ix_layer, iy + 1, iz + 1)

                if ix_layer == 0:
                    zone = "solid_left"
                elif ix_layer == 4:
                    zone = "solid_right"
                elif has_overlap:
                    if is_left_fin and ix_layer in (1, 2):
                        zone = "solid_left"
                    elif is_right_fin and ix_layer in (2, 3):
                        zone = "solid_right"
                    else:
                        zone = "fluid"
                else:
                    if is_left_fin and ix_layer == 1:
                        zone = "solid_left"
                    elif is_right_fin and ix_layer == 3:
                        zone = "solid_right"
                    else:
                        zone = "fluid"

                blocks.append(
                    "    "
                    + f"hex ({v0} {v1} {v2} {v3} {v4} {v5} {v6} {v7}) "
                    + f"{zone} ({nx_layer[ix_layer]} {ny_mesh} {nz_mesh}) simpleGrading (1 1 1)"
                )

    left_wall_faces: list[str] = []
    right_wall_faces: list[str] = []
    side_wall_faces: list[str] = []

    for iz in range(n_z - 1):
        for iy in range(n_y - 1):
            v0 = v_idx(0, iy, iz)
            v3 = v_idx(0, iy + 1, iz)
            v7 = v_idx(0, iy + 1, iz + 1)
            v4 = v_idx(0, iy, iz + 1)
            left_wall_faces.append(f"            ({v0} {v4} {v7} {v3})")

            v1 = v_idx(n_x - 1, iy, iz)
            v2 = v_idx(n_x - 1, iy + 1, iz)
            v6 = v_idx(n_x - 1, iy + 1, iz + 1)
            v5 = v_idx(n_x - 1, iy, iz + 1)
            right_wall_faces.append(f"            ({v1} {v2} {v6} {v5})")

    for iz in range(n_z - 1):
        for ix in range(n_x - 1):
            v0 = v_idx(ix, 0, iz)
            v1 = v_idx(ix + 1, 0, iz)
            v5 = v_idx(ix + 1, 0, iz + 1)
            v4 = v_idx(ix, 0, iz + 1)
            side_wall_faces.append(f"            ({v0} {v1} {v5} {v4})")

    for iz in range(n_z - 1):
        for ix in range(n_x - 1):
            v3 = v_idx(ix, n_y - 1, iz)
            v2 = v_idx(ix + 1, n_y - 1, iz)
            v6 = v_idx(ix + 1, n_y - 1, iz + 1)
            v7 = v_idx(ix, n_y - 1, iz + 1)
            side_wall_faces.append(f"            ({v3} {v7} {v6} {v2})")

    for iy in range(n_y - 1):
        for ix in range(n_x - 1):
            v0 = v_idx(ix, iy, 0)
            v1 = v_idx(ix + 1, iy, 0)
            v2 = v_idx(ix + 1, iy + 1, 0)
            v3 = v_idx(ix, iy + 1, 0)
            side_wall_faces.append(f"            ({v0} {v3} {v2} {v1})")

    for iy in range(n_y - 1):
        for ix in range(n_x - 1):
            v4 = v_idx(ix, iy, n_z - 1)
            v5 = v_idx(ix + 1, iy, n_z - 1)
            v6 = v_idx(ix + 1, iy + 1, n_z - 1)
            v7 = v_idx(ix, iy + 1, n_z - 1)
            side_wall_faces.append(f"            ({v4} {v5} {v6} {v7})")

    vertices_block = "\n".join(vertices)
    blocks_block = "\n".join(blocks)
    left_faces_block = "\n".join(left_wall_faces)
    right_faces_block = "\n".join(right_wall_faces)
    side_faces_block = "\n".join(side_wall_faces)

    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     |                                                 |
|   \\  /    A nd           |                                                 |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}
// Generated by openfoam_automation.py
// Reference-style structured layered decomposition (blockMesh-only)
// Real span after unit rounding: H={h_real:.9f}, W={w_real:.9f}

scale 1;

vertices
(
{vertices_block}
);

blocks
(
{blocks_block}
);

edges
();

boundary
(
    leftWallOuter
    {{
        type wall;
        faces
        (
{left_faces_block}
        );
    }}

    rightWallOuter
    {{
        type wall;
        faces
        (
{right_faces_block}
        );
    }}

    walls
    {{
        type wall;
        faces
        (
{side_faces_block}
        );
    }}
);

mergePatchPairs
(
);

// ************************************************************************* //
"""


def generate_controldict_text() -> str:
    return """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     |                                                 |
|   \\  /    A nd           |                                                 |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}
// Generated by openfoam_automation.py

#includeIfPresent "regionSolvers"

startFrom       startTime;
startTime       0;

stopAt          endTime;
endTime         500;

deltaT          1;

writeControl    timeStep;
writeInterval   50;

purgeWrite      5;

writeFormat     ascii;
writePrecision  8;
writeCompression off;

timeFormat      general;
timePrecision   6;

runTimeModifiable true;

// ************************************************************************* //
"""


def generate_fvschemes_text() -> str:
    return """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     |                                                 |
|   \\  /    A nd           |                                                 |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}
// Generated by openfoam_automation.py

ddtSchemes
{
}

gradSchemes
{
}

divSchemes
{
}

laplacianSchemes
{
}

interpolationSchemes
{
}

snGradSchemes
{
}

fluxRequired
{
}

// ************************************************************************* //
"""


def generate_fvsolution_text() -> str:
    return """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     |                                                 |
|   \\  /    A nd           |                                                 |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
// Generated by openfoam_automation.py

PIMPLE
{
    nOuterCorrectors 1;
}

// ************************************************************************* //
"""


def generate_materialproperties_text(
    *,
    fluid_material: str = "air",
    solid_left_material: str = "copper",
    solid_right_material: str = "copper",
) -> str:
    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     |                                                 |
|   \\  /    A nd           |                                                 |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      materialProperties;
}}
// Generated by openfoam_automation.py

fluid
{{
    solver          fluid;
    material        {fluid_material};
}}

solid_left
{{
    solver          solid;
    material        {solid_left_material};
}}

solid_right
{{
    solver          solid;
    material        {solid_right_material};
}}

// ************************************************************************* //
"""


def write_openfoam_case_files(
    case_dir: Path, geom: GeometryParams, mesh: MeshingParams
) -> None:
    system_dir = case_dir / "system"
    constant_dir = case_dir / "constant"
    system_dir.mkdir(parents=True, exist_ok=True)
    constant_dir.mkdir(parents=True, exist_ok=True)

    try:
        (system_dir / "regionSolvers").unlink()
    except FileNotFoundError:
        pass

    for stale in [
        "surfaceFeaturesDict",
        "surfaceFeatureExtractDict",
        "snappyHexMeshDict",
    ]:
        try:
            (system_dir / stale).unlink()
        except FileNotFoundError:
            pass

    (system_dir / "controlDict").write_text(
        generate_controldict_text(), encoding="utf-8"
    )
    (system_dir / "fvSchemes").write_text(generate_fvschemes_text(), encoding="utf-8")
    (system_dir / "fvSolution").write_text(generate_fvsolution_text(), encoding="utf-8")
    (system_dir / "blockMeshDict").write_text(
        generate_blockmeshdict_text(geom, mesh), encoding="utf-8"
    )
    (constant_dir / "materialProperties").write_text(
        generate_materialproperties_text(), encoding="utf-8"
    )

    allmesh = (
        "#!/bin/sh\n"
        "set -eu\n"
        "blockMesh > log.blockMesh 2>&1\n"
        "splitMeshRegions -cellZones -overwrite -defaultRegionName fluid > log.splitMeshRegions 2>&1\n"
    )
    allmesh_path = case_dir / "Allmesh"
    allmesh_path.write_text(allmesh, encoding="utf-8")
    try:
        allmesh_path.chmod(0o755)
    except OSError:
        pass

    allrun = (
        "#!/bin/sh\n"
        "set -eu\n"
        "./Allmesh\n"
        "foamSetupCHT > log.foamSetupCHT 2>&1\n"
        "foamMultiRun > log.foamMultiRun 2>&1\n"
        "paraFoam -touchAll > /dev/null 2>&1 || true\n"
    )
    allrun_path = case_dir / "Allrun"
    allrun_path.write_text(allrun, encoding="utf-8")
    try:
        allrun_path.chmod(0o755)
    except OSError:
        pass


def run_openfoam_meshing(case_dir: Path) -> None:
    if shutil.which("blockMesh") is None:
        raise FileNotFoundError(
            "blockMesh not found in PATH (did you source OpenFOAM?)"
        )
    if shutil.which("splitMeshRegions") is None:
        raise FileNotFoundError(
            "splitMeshRegions not found in PATH (did you source OpenFOAM?)"
        )

    subprocess.run(["blockMesh"], cwd=str(case_dir), check=True)
    subprocess.run(
        [
            "splitMeshRegions",
            "-cellZones",
            "-overwrite",
            "-defaultRegionName",
            "fluid",
        ],
        cwd=str(case_dir),
        check=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate reference-style blockMeshDict for checkerboard fin walls and run blockMesh + splitMeshRegions."
    )
    parser.add_argument("--case-dir", default=".")
    parser.add_argument("--write-only", action="store_true")

    parser.add_argument("--gap", type=float, default=0.20)
    parser.add_argument("--T", type=float, default=0.01)
    parser.add_argument("--H-wall", type=float, default=0.10)
    parser.add_argument("--W-wall", type=float, default=0.10)
    parser.add_argument("--p", type=float, default=0.02)
    parser.add_argument("--t", type=float, default=0.005)
    parser.add_argument("--L", type=float, default=0.12)

    parser.add_argument("--N-T", type=int, default=5)
    parser.add_argument("--N-gapA", type=int, default=3)
    parser.add_argument("--N-over", type=int, default=10)
    parser.add_argument("--N-t", type=int, default=5)
    parser.add_argument("--N-p", type=int, default=5)

    args = parser.parse_args()

    case_dir = Path(args.case_dir).resolve()
    geom = GeometryParams(
        gap=args.gap,
        T=args.T,
        H_wall=args.H_wall,
        W_wall=args.W_wall,
        p=args.p,
        t=args.t,
        L=args.L,
    )
    mesh = MeshingParams(
        N_T=args.N_T,
        N_gapA=args.N_gapA,
        N_over=args.N_over,
        N_t=args.N_t,
        N_p=args.N_p,
    )

    write_openfoam_case_files(case_dir, geom, mesh)
    if not args.write_only:
        run_openfoam_meshing(case_dir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
