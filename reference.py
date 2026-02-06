import math

# =============================================================================
# [1] 사용자 파라미터 입력 (User Parameters)
# =============================================================================

output_file = "system/blockMeshDict"

# --- Geometry Parameters (mm 단위) ---
gap = 300  # 두 wall 사이 거리 (Fluid Gap)
T = 25.0  # Wall 두께
H_target = 100.0  # 목표 Wall 높이 (Y)
W_target = 50.0  # 목표 Wall 폭 (Z)

t = 1.0  # Fin 두께
p = 1.0  # Fin 사이 거리 (Clearance)
L = 240.0  # Fin 길이 (x방향 돌출 길이)

# --- Mesh Density (격자 개수) ---
N_T = 5  # Wall 두께 방향 격자 수
N_gapA = 3  # Wall과 Fin Tip 사이 (Gap A) 격자 수
N_over = 10  # Fin 겹침 구간 (Overlap) 격자 수
N_t = 5  # Fin 두께(t) 방향 격자 수
N_p = 5  # Fin 간격(p) 방향 격자 수

# =============================================================================
# [2] 파라미터 검증 및 자동 보정
# =============================================================================

if L <= gap / 2:
    print(f"Error: L ({L}) must be > gap/2 ({gap / 2}) for overlap.")
    exit(1)
if L >= gap:
    print(f"Warning: L ({L}) >= gap ({gap}). Reducing L slightly.")
    L = gap - 1e-5

pitch = t + p
Ny_units = max(2, round(H_target / pitch))
Nz_units = max(2, round(W_target / pitch))

H_real = Ny_units * pitch
W_real = Nz_units * pitch

print(f"--- Geometry Info ---")
print(f"Real Size   : H={H_real}, W={W_real}")
print(f"Units (Y,Z) : {Ny_units} x {Nz_units}")
print(f"---------------------")

# =============================================================================
# [3] 좌표 생성 (Coordinate Generation)
# =============================================================================

x_coords = [
    -T,  # x0
    0.0,  # x1 (Left Wall Inner)
    gap - L,  # x2 (Right Fin Tip)
    L,  # x3 (Left Fin Tip)
    gap,  # x4 (Right Wall Inner)
    gap + T  # x5
]

y_coords = [0.0]
for i in range(Ny_units):
    last = y_coords[-1]
    y_coords.append(last + t)
    y_coords.append(last + t + p)

z_coords = [0.0]
for i in range(Nz_units):
    last = z_coords[-1]
    z_coords.append(last + t)
    z_coords.append(last + t + p)

nX = len(x_coords)
nY = len(y_coords)
nZ = len(z_coords)


# =============================================================================
# [4] blockMeshDict 작성 함수
# =============================================================================

def get_vtx_idx(ix, iy, iz):
    return ix + (nX * iy) + (nX * nY * iz)


with open(output_file, 'w') as f:
    f.write(r"""/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scale 0.001; // mm

""")

    # Vertices
    f.write("vertices\n(\n")
    for iz in range(nZ):
        for iy in range(nY):
            for ix in range(nX):
                f.write(f"    ({x_coords[ix]:.6f} {y_coords[iy]:.6f} {z_coords[iz]:.6f})\n")
    f.write(");\n\n")

    # Blocks
    f.write("blocks\n(\n")

    nx_layer = [N_T, N_gapA, N_over, N_gapA, N_T]

    for iz in range(nZ - 1):
        for iy in range(nY - 1):

            # 현재 구간이 Fin 두께(t) 구간인지 여부
            is_y_fin = (iy % 2 == 0)
            is_z_fin = (iz % 2 == 0)

            ny_mesh = N_t if is_y_fin else N_p
            nz_mesh = N_t if is_z_fin else N_p

            is_fin_site = is_y_fin and is_z_fin

            grid_y = iy // 2
            grid_z = iz // 2

            is_left_fin = False
            is_right_fin = False

            if is_fin_site:
                if (grid_y + grid_z) % 2 == 0:
                    is_left_fin = True
                else:
                    is_right_fin = True

            for ix_layer in range(5):
                v0 = get_vtx_idx(ix_layer, iy, iz)
                v1 = get_vtx_idx(ix_layer + 1, iy, iz)
                v2 = get_vtx_idx(ix_layer + 1, iy + 1, iz)
                v3 = get_vtx_idx(ix_layer, iy + 1, iz)

                v4 = get_vtx_idx(ix_layer, iy, iz + 1)
                v5 = get_vtx_idx(ix_layer + 1, iy, iz + 1)
                v6 = get_vtx_idx(ix_layer + 1, iy + 1, iz + 1)
                v7 = get_vtx_idx(ix_layer, iy + 1, iz + 1)

                region_name = "fluid"

                # 영역 할당 로직 (Expanded for safety)
                if ix_layer == 0:
                    region_name = "leftSolid"
                elif ix_layer == 4:
                    region_name = "rightSolid"
                else:
                    # 중간 3개 레이어 (GapA, Overlap, GapC)
                    if is_left_fin:
                        if ix_layer == 1:
                            region_name = "leftSolid"
                        elif ix_layer == 2:
                            region_name = "leftSolid"
                        else:
                            region_name = "fluid"
                    elif is_right_fin:
                        if ix_layer == 1:
                            region_name = "fluid"
                        elif ix_layer == 2:
                            region_name = "rightSolid"
                        else:
                            region_name = "rightSolid"
                    else:
                        region_name = "fluid"

                f.write(f"    hex ({v0} {v1} {v2} {v3} {v4} {v5} {v6} {v7}) ")
                f.write(f" {region_name} ({nx_layer[ix_layer]} {ny_mesh} {nz_mesh}) simpleGrading (1 1 1) \n")

    f.write(");\n\n")

    # Boundary
    f.write("boundary\n(\n")

    # Left Wall Outer
    f.write("    leftWallOuter\n    {\n        type wall;\n        faces\n        (\n")
    for iz in range(nZ - 1):
        for iy in range(nY - 1):
            v0 = get_vtx_idx(0, iy, iz)
            v3 = get_vtx_idx(0, iy + 1, iz)
            v7 = get_vtx_idx(0, iy + 1, iz + 1)
            v4 = get_vtx_idx(0, iy, iz + 1)
            f.write(f"            ({v0} {v4} {v7} {v3})\n")
    f.write("        );\n    }\n")

    # Right Wall Outer
    f.write("    rightWallOuter\n    {\n        type wall;\n        faces\n        (\n")
    for iz in range(nZ - 1):
        for iy in range(nY - 1):
            v1 = get_vtx_idx(nX - 1, iy, iz)
            v2 = get_vtx_idx(nX - 1, iy + 1, iz)
            v6 = get_vtx_idx(nX - 1, iy + 1, iz + 1)
            v5 = get_vtx_idx(nX - 1, iy, iz + 1)
            f.write(f"            ({v1} {v2} {v6} {v5})\n")
    f.write("        );\n    }\n")

    # Side Walls (Combined)
    f.write("    walls\n    {\n        type wall;\n        faces\n        (\n")

    # Bottom (Y=min)
    for iz in range(nZ - 1):
        for ix in range(nX - 1):
            v0 = get_vtx_idx(ix, 0, iz)
            v1 = get_vtx_idx(ix + 1, 0, iz)
            v5 = get_vtx_idx(ix + 1, 0, iz + 1)
            v4 = get_vtx_idx(ix, 0, iz + 1)
            f.write(f"            ({v0} {v1} {v5} {v4})\n")

    # Top (Y=max)
    for iz in range(nZ - 1):
        for ix in range(nX - 1):
            v3 = get_vtx_idx(ix, nY - 1, iz)
            v2 = get_vtx_idx(ix + 1, nY - 1, iz)
            v6 = get_vtx_idx(ix + 1, nY - 1, iz + 1)
            v7 = get_vtx_idx(ix, nY - 1, iz + 1)
            f.write(f"            ({v3} {v7} {v6} {v2})\n")

    # Front (Z=min)
    for iy in range(nY - 1):
        for ix in range(nX - 1):
            v0 = get_vtx_idx(ix, iy, 0)
            v1 = get_vtx_idx(ix + 1, iy, 0)
            v2 = get_vtx_idx(ix + 1, iy + 1, 0)
            v3 = get_vtx_idx(ix, iy + 1, 0)
            f.write(f"            ({v0} {v3} {v2} {v1})\n")

    # Back (Z=max)
    for iy in range(nY - 1):
        for ix in range(nX - 1):
            v4 = get_vtx_idx(ix, iy, nZ - 1)
            v5 = get_vtx_idx(ix + 1, iy, nZ - 1)
            v6 = get_vtx_idx(ix + 1, iy + 1, nZ - 1)
            v7 = get_vtx_idx(ix, iy + 1, nZ - 1)
            f.write(f"            ({v4} {v5} {v6} {v7})\n")

    f.write("        );\n    }\n")
    f.write(");\n\n")

    f.write("mergePatchPairs\n(\n);\n")

print(f"Successfully generated {output_file}")
