import numpy as np
from stl import mesh
import os
import math
from typing import Any


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

    # 교차(겹침) + 반대쪽 wall 관통 방지 조건
    if not (L > gap / 2):
        raise ValueError("To make fins cross in the middle, require L > gap/2")
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

    # === 1. 왼쪽 Wall 생성 ===
    wall_left_xmin = -gap / 2 - T
    wall_left_xmax = -gap / 2
    wall_left_ymin = 0
    wall_left_ymax = H_wall
    wall_left_zmin = 0
    wall_left_zmax = W_wall

    wall_left = create_box_stl(
        wall_left_xmin,
        wall_left_xmax,
        wall_left_ymin,
        wall_left_ymax,
        wall_left_zmin,
        wall_left_zmax,
    )
    wall_left.save(f"{output_dir}/wall_left.stl")
    print(f"생성됨: wall_left.stl")

    # === 2. 오른쪽 Wall 생성 ===
    wall_right_xmin = gap / 2
    wall_right_xmax = gap / 2 + T

    wall_right = create_box_stl(
        wall_right_xmin,
        wall_right_xmax,
        wall_left_ymin,
        wall_left_ymax,
        wall_left_zmin,
        wall_left_zmax,
    )
    wall_right.save(f"{output_dir}/wall_right.stl")
    print(f"생성됨: wall_right.stl")

    # === 핀 생성 (y-z 단면에서 checkerboard 형태로 좌/우 wall에 번갈아 배치) ===
    print(f"\n=== 자동 배치 결과 ===")
    print(
        f"y-z 그리드: {num_y} x {num_z} (총 {total_cells}) (좌 {left_count}, 우 {right_count})"
    )
    print(f"핀 단면 크기: {t} x {t}, 핀 간격(클리어런스): {p}\n")

    for iz in range(num_z):
        fin_zmin = z_start + iz * (t + p)
        fin_zmax = fin_zmin + t
        for iy in range(num_y):
            fin_ymin = y_start + iy * (t + p)
            fin_ymax = fin_ymin + t

            if (iy + iz) % 2 == 0:
                fin_xmin = -gap / 2
                fin_xmax = -gap / 2 + L
                fin = create_box_stl(
                    fin_xmin, fin_xmax, fin_ymin, fin_ymax, fin_zmin, fin_zmax
                )
                fin.save(f"{output_dir}/fin_left_y{iy}_z{iz}.stl")
                print(f"생성됨: fin_left_y{iy}_z{iz}.stl")
            else:
                fin_xmin = gap / 2 - L
                fin_xmax = gap / 2
                fin = create_box_stl(
                    fin_xmin, fin_xmax, fin_ymin, fin_ymax, fin_zmin, fin_zmax
                )
                fin.save(f"{output_dir}/fin_right_y{iy}_z{iz}.stl")
                print(f"생성됨: fin_right_y{iy}_z{iz}.stl")

    print(f"\n총 {2 + total_cells}개 STL 파일 생성 완료")
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
