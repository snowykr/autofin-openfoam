import math
import os
from dataclasses import dataclass
from typing import Any

import numpy as np
from stl import mesh


def create_box_stl(
    xmin: float,
    xmax: float,
    ymin: float,
    ymax: float,
    zmin: float,
    zmax: float,
) -> Any:
    """육면체 STL 메쉬 생성"""
    vertices = np.array(
        [
            [xmin, ymin, zmin],
            [xmax, ymin, zmin],
            [xmax, ymax, zmin],
            [xmin, ymax, zmin],
            [xmin, ymin, zmax],
            [xmax, ymin, zmax],
            [xmax, ymax, zmax],
            [xmin, ymax, zmax],
        ]
    )

    # 육면체를 구성하는 12개 삼각형 (STL 형식)
    faces = np.array(
        [
            [0, 3, 1],
            [1, 3, 2],  # 아래면 (z=min)
            [4, 5, 6],
            [4, 6, 7],  # 위면 (z=max)
            [0, 1, 5],
            [0, 5, 4],  # 앞면 (y=min)
            [2, 3, 7],
            [2, 7, 6],  # 뒷면 (y=max)
            [0, 4, 7],
            [0, 7, 3],  # 왼쪽면 (x=min)
            [1, 2, 6],
            [1, 6, 5],  # 오른쪽면 (x=max)
        ]
    )

    box_mesh: Any = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            box_mesh.vectors[i][j] = vertices[f[j], :]

    return box_mesh


@dataclass(frozen=True)
class FinFootprint:
    ymin: float
    ymax: float
    zmin: float
    zmax: float


def _add_rect(
    tris: list[np.ndarray],
    v0: tuple[float, float, float],
    v1: tuple[float, float, float],
    v2: tuple[float, float, float],
    v3: tuple[float, float, float],
) -> None:
    """Append two triangles for a rectangle (v0->v1->v2->v3 CCW as seen from outside)."""

    tris.append(np.array([v0, v1, v2], dtype=float))
    tris.append(np.array([v0, v2, v3], dtype=float))


def _add_face_x(
    tris: list[np.ndarray],
    x: float,
    y0: float,
    y1: float,
    z0: float,
    z1: float,
    outward: str,
) -> None:
    if outward == "+x":
        _add_rect(
            tris,
            (x, y0, z0),
            (x, y1, z0),
            (x, y1, z1),
            (x, y0, z1),
        )
        return
    if outward == "-x":
        _add_rect(
            tris,
            (x, y0, z0),
            (x, y0, z1),
            (x, y1, z1),
            (x, y1, z0),
        )
        return
    raise ValueError(f"Invalid outward={outward!r} for x-face")


def _add_face_y(
    tris: list[np.ndarray],
    y: float,
    x0: float,
    x1: float,
    z0: float,
    z1: float,
    outward: str,
) -> None:
    if outward == "+y":
        _add_rect(
            tris,
            (x0, y, z0),
            (x0, y, z1),
            (x1, y, z1),
            (x1, y, z0),
        )
        return
    if outward == "-y":
        _add_rect(
            tris,
            (x0, y, z0),
            (x1, y, z0),
            (x1, y, z1),
            (x0, y, z1),
        )
        return
    raise ValueError(f"Invalid outward={outward!r} for y-face")


def _add_face_z(
    tris: list[np.ndarray],
    z: float,
    x0: float,
    x1: float,
    y0: float,
    y1: float,
    outward: str,
) -> None:
    if outward == "+z":
        _add_rect(
            tris,
            (x0, y0, z),
            (x1, y0, z),
            (x1, y1, z),
            (x0, y1, z),
        )
        return
    if outward == "-z":
        _add_rect(
            tris,
            (x0, y0, z),
            (x0, y1, z),
            (x1, y1, z),
            (x1, y0, z),
        )
        return
    raise ValueError(f"Invalid outward={outward!r} for z-face")


def _mesh_from_tris(tris: list[np.ndarray]) -> Any:
    m: Any = mesh.Mesh(np.zeros(len(tris), dtype=mesh.Mesh.dtype))
    for i, tri in enumerate(tris):
        m.vectors[i] = tri
    return m


def _unique_sorted(vals: list[float], eps: float = 1e-12) -> list[float]:
    if not vals:
        return []
    out: list[float] = []
    for v in sorted(vals):
        if not out or abs(v - out[-1]) > eps:
            out.append(v)
        else:
            # Snap near-equals to the first value (stabilize grid splits)
            out[-1] = out[-1]
    return out


def _in_fin(fp: FinFootprint, y: float, z: float, eps: float = 1e-15) -> bool:
    return (fp.ymin - eps) <= y <= (fp.ymax + eps) and (fp.zmin - eps) <= z <= (
        fp.zmax + eps
    )


def _inner_face_patches(
    y_min: float,
    y_max: float,
    z_min: float,
    z_max: float,
    fins: list[FinFootprint],
) -> list[tuple[float, float, float, float]]:
    """Return rectangle patches on the wall inner face, with fin footprints removed.

    Returned rectangles are in (y0, y1, z0, z1) with y0<y1, z0<z1.
    """

    if not fins:
        return [(y_min, y_max, z_min, z_max)]

    y_coords = [y_min, y_max]
    z_coords = [z_min, z_max]
    for fp in fins:
        y_coords.extend([fp.ymin, fp.ymax])
        z_coords.extend([fp.zmin, fp.zmax])

    ys = _unique_sorted(y_coords)
    zs = _unique_sorted(z_coords)
    patches: list[tuple[float, float, float, float]] = []
    for i in range(len(ys) - 1):
        y0, y1 = ys[i], ys[i + 1]
        if not (y0 < y1):
            continue
        y_mid = 0.5 * (y0 + y1)
        for j in range(len(zs) - 1):
            z0, z1 = zs[j], zs[j + 1]
            if not (z0 < z1):
                continue
            z_mid = 0.5 * (z0 + z1)

            if any(_in_fin(fp, y_mid, z_mid) for fp in fins):
                continue
            patches.append((y0, y1, z0, z1))

    return patches


def _yz_grid_from_fins(
    y_min: float,
    y_max: float,
    z_min: float,
    z_max: float,
    fins: list[FinFootprint],
) -> tuple[list[float], list[float]]:
    y_coords = [y_min, y_max]
    z_coords = [z_min, z_max]
    for fp in fins:
        y_coords.extend([fp.ymin, fp.ymax])
        z_coords.extend([fp.zmin, fp.zmax])
    return (_unique_sorted(y_coords), _unique_sorted(z_coords))


def create_solid_wall_fin_stl(
    *,
    side: str,
    gap: float,
    T: float,
    H_wall: float,
    W_wall: float,
    L: float,
    fins: list[FinFootprint],
) -> Any:
    """Create a single watertight triSurface for wall+fins on one side.

    This avoids coplanar duplicate faces at the wall-fin interface by:
    - perforating the wall inner face with fin footprints
    - omitting the fin base faces on the wall inner plane
    """

    if side not in {"left", "right"}:
        raise ValueError("side must be 'left' or 'right'")

    # Wall extents
    if side == "left":
        x_inner = -gap / 2
        x_outer = -gap / 2 - T
        wall_inner_outward = "+x"
        fin_x0 = x_inner
        fin_x1 = x_inner + L
    else:
        x_inner = gap / 2
        x_outer = gap / 2 + T
        wall_inner_outward = "-x"
        fin_x0 = x_inner - L
        fin_x1 = x_inner

    x0, x1 = (x_outer, x_inner) if x_outer < x_inner else (x_inner, x_outer)

    tris: list[np.ndarray] = []

    # Build a shared (y,z) grid so all wall faces share vertices.
    # Without this, the perforated inner face creates T-junctions against
    # the non-subdivided side faces, which shows up as open edges.
    ys, zs = _yz_grid_from_fins(0.0, H_wall, 0.0, W_wall, fins)

    # --- Wall faces ---
    # Outer face (subdivided on the shared grid)
    outer_outward = "-x" if side == "left" else "+x"
    for i in range(len(ys) - 1):
        for j in range(len(zs) - 1):
            _add_face_x(
                tris,
                x_outer,
                ys[i],
                ys[i + 1],
                zs[j],
                zs[j + 1],
                outward=outer_outward,
            )

    # y-min / y-max (subdivide along z)
    for j in range(len(zs) - 1):
        _add_face_y(tris, 0.0, x0, x1, zs[j], zs[j + 1], outward="-y")
        _add_face_y(tris, H_wall, x0, x1, zs[j], zs[j + 1], outward="+y")

    # z-min / z-max (subdivide along y)
    for i in range(len(ys) - 1):
        _add_face_z(tris, 0.0, x0, x1, ys[i], ys[i + 1], outward="-z")
        _add_face_z(tris, W_wall, x0, x1, ys[i], ys[i + 1], outward="+z")

    # Inner face (perforated on the same shared grid)
    for i in range(len(ys) - 1):
        y0p, y1p = ys[i], ys[i + 1]
        y_mid = 0.5 * (y0p + y1p)
        for j in range(len(zs) - 1):
            z0p, z1p = zs[j], zs[j + 1]
            z_mid = 0.5 * (z0p + z1p)
            if any(_in_fin(fp, y_mid, z_mid) for fp in fins):
                continue
            _add_face_x(
                tris,
                x_inner,
                y0p,
                y1p,
                z0p,
                z1p,
                outward=wall_inner_outward,
            )

    # --- Fin faces (no base face on x_inner plane) ---
    # Normalize fin x range for side faces
    fx0 = min(fin_x0, fin_x1)
    fx1 = max(fin_x0, fin_x1)

    for fp in fins:
        # Tip face
        x_tip = fin_x1 if side == "left" else fin_x0
        tip_outward = "+x" if side == "left" else "-x"
        _add_face_x(
            tris, x_tip, fp.ymin, fp.ymax, fp.zmin, fp.zmax, outward=tip_outward
        )

        # y faces
        _add_face_y(tris, fp.ymin, fx0, fx1, fp.zmin, fp.zmax, outward="-y")
        _add_face_y(tris, fp.ymax, fx0, fx1, fp.zmin, fp.zmax, outward="+y")

        # z faces
        _add_face_z(tris, fp.zmin, fx0, fx1, fp.ymin, fp.ymax, outward="-z")
        _add_face_z(tris, fp.zmax, fx0, fx1, fp.ymin, fp.ymax, outward="+z")

    return _mesh_from_tris(tris)


def _max_count_with_gap(total_length: float, thickness: float, gap: float) -> int:
    if total_length <= 0:
        return 0
    if thickness <= 0:
        return 0
    if gap < 0:
        return 0
    return int(math.floor((total_length + gap) / (thickness + gap) + 1e-12))


def generate_fin_wall_geometry(
    gap: float,
    T: float,
    H_wall: float,
    W_wall: float,
    p: float,
    t: float,
    L: float,
    output_dir: str = "constant/triSurface",
) -> None:
    """
    핀-벽 구조 STL 파일 자동 생성

    Parameters:
    -----------
    gap : float
        wall과 wall 사이의 공간
    T : float
        wall의 두께
    H_wall : float
        wall의 세로 길이 (y 방향)
    W_wall : float
        wall의 가로 길이 (z 방향)
    p : float
        fin 사이의 공간 (세로 방향과 레이어 사이 모두 적용)
    t : float
        fin의 두께
    L : float
        fin의 길이
    내부 계산:
    ---------
    핀 개수/레이어 개수는 주어진 (H_wall, W_wall, p, t)로부터 자동으로 최대화하여 계산
    output_dir : str
        STL 파일 저장 디렉토리
    """

    # 출력 디렉토리 생성
    os.makedirs(output_dir, exist_ok=True)

    # 좌표계 정의
    # x: wall 사이 방향
    # y: wall의 세로 방향
    # z: 레이어 방향

    # === 파라미터 검증 + 핀 배치 최대화 계산 ===
    if gap <= 0 or T <= 0 or H_wall <= 0 or W_wall <= 0:
        raise ValueError("gap, T, H_wall, W_wall must be > 0")
    if p < 0 or t <= 0 or L <= 0:
        raise ValueError("p must be >= 0, and t, L must be > 0")

    # 반대쪽 wall 관통 방지 조건
    if L > gap:
        raise ValueError("To avoid intersecting the opposite wall, require L <= gap")

    num_y = _max_count_with_gap(H_wall, t, p)
    num_z = _max_count_with_gap(W_wall, t, p)
    if num_y < 1:
        raise ValueError("Not enough space in H_wall to place any fins")
    if num_z < 1:
        raise ValueError("Not enough space in W_wall to place any fins")

    total_cells = num_y * num_z
    left_count = (total_cells + 1) // 2
    right_count = total_cells // 2
    if left_count < 1 or right_count < 1:
        raise ValueError("Not enough space to place at least one fin on each wall")

    used_y = num_y * t + (num_y - 1) * p
    y_start = (H_wall - used_y) / 2

    used_z = num_z * t + (num_z - 1) * p
    z_start = (W_wall - used_z) / 2

    # === 핀 배치 (y-z 단면에서 checkerboard 형태로 좌/우 wall에 번갈아 배치) ===
    print(f"\n=== 자동 배치 결과 ===")
    print(
        f"y-z 그리드: {num_y} x {num_z} (총 {total_cells}) (좌 {left_count}, 우 {right_count})"
    )
    print(f"핀 단면 크기: {t} x {t}, 핀 간격(클리어런스): {p}\n")

    left_fins: list[FinFootprint] = []
    right_fins: list[FinFootprint] = []

    for iz in range(num_z):
        fin_zmin = z_start + iz * (t + p)
        fin_zmax = fin_zmin + t
        for iy in range(num_y):
            fin_ymin = y_start + iy * (t + p)
            fin_ymax = fin_ymin + t

            fp = FinFootprint(fin_ymin, fin_ymax, fin_zmin, fin_zmax)
            if (iy + iz) % 2 == 0:
                left_fins.append(fp)
            else:
                right_fins.append(fp)

    # === 단일 solid 표면(STL) 생성: wall+fin (좌/우 각각) ===
    solid_left = create_solid_wall_fin_stl(
        side="left",
        gap=gap,
        T=T,
        H_wall=H_wall,
        W_wall=W_wall,
        L=L,
        fins=left_fins,
    )
    solid_left.save(f"{output_dir}/solid_left.stl")
    print("생성됨: solid_left.stl")

    solid_right = create_solid_wall_fin_stl(
        side="right",
        gap=gap,
        T=T,
        H_wall=H_wall,
        W_wall=W_wall,
        L=L,
        fins=right_fins,
    )
    solid_right.save(f"{output_dir}/solid_right.stl")
    print("생성됨: solid_right.stl")

    print("\n총 2개 STL 파일 생성 완료")
    print(f"저장 위치: {output_dir}/")


# ========== 메인 실행 부분 ==========
if __name__ == "__main__":
    # 파라미터 설정
    gap = 0.20  # wall 사이 공간 (m)
    T = 0.01  # wall 두께 (m)
    H_wall = 0.10  # wall 세로 길이 (m)
    W_wall = 0.10  # wall 가로 길이 (m)
    p = 0.02  # fin 사이 간격 (m)
    t = 0.005  # fin 두께 (m)
    L = 0.12  # fin 길이 (m)

    print("=== 파라미터 ===")
    print(f"wall 간격 (gap): {gap} m")
    print(f"wall 두께 (T): {T} m")
    print(f"wall 크기: {H_wall} x {W_wall} m")
    print(f"fin 간격 (p): {p} m")
    print(f"fin 두께 (t): {t} m")
    print(f"fin 길이 (L): {L} m")
    print("핀/레이어 개수는 자동으로 최대화하여 계산됩니다.\n")

    # STL 파일 생성
    generate_fin_wall_geometry(
        gap=gap,
        T=T,
        H_wall=H_wall,
        W_wall=W_wall,
        p=p,
        t=t,
        L=L,
    )
