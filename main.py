import numpy as np
from stl import mesh
import os

def create_box_stl(xmin, xmax, ymin, ymax, zmin, zmax):
    """육면체 STL 메쉬 생성"""
    vertices = np.array([
        [xmin, ymin, zmin],
        [xmax, ymin, zmin],
        [xmax, ymax, zmin],
        [xmin, ymax, zmin],
        [xmin, ymin, zmax],
        [xmax, ymin, zmax],
        [xmax, ymax, zmax],
        [xmin, ymax, zmax]
    ])
    
    # 육면체를 구성하는 12개 삼각형 (STL 형식)
    faces = np.array([
        [0,3,1], [1,3,2],  # 아래면 (z=min)
        [4,5,6], [4,6,7],  # 위면 (z=max)
        [0,1,5], [0,5,4],  # 앞면 (y=min)
        [2,3,7], [2,7,6],  # 뒷면 (y=max)
        [0,4,7], [0,7,3],  # 왼쪽면 (x=min)
        [1,2,6], [1,6,5]   # 오른쪽면 (x=max)
    ])
    
    box_mesh = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            box_mesh.vectors[i][j] = vertices[f[j],:]
    
    return box_mesh


def generate_fin_wall_geometry(gap, T, H_wall, W_wall, p, t, L, 
                                num_fins_per_side, num_layers, 
                                output_dir="constant/triSurface"):
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
    num_fins_per_side : int
        한쪽 wall당 핀 개수
    num_layers : int
        z 방향 레이어 개수
    output_dir : str
        STL 파일 저장 디렉토리
    """
    
    # 출력 디렉토리 생성
    os.makedirs(output_dir, exist_ok=True)
    
    # 좌표계 정의
    # x: wall 사이 방향
    # y: wall의 세로 방향
    # z: 레이어 방향
    
    # === 1. 왼쪽 Wall 생성 ===
    wall_left_xmin = -gap/2 - T
    wall_left_xmax = -gap/2
    wall_left_ymin = 0
    wall_left_ymax = H_wall
    wall_left_zmin = 0
    wall_left_zmax = W_wall
    
    wall_left = create_box_stl(wall_left_xmin, wall_left_xmax, 
                                wall_left_ymin, wall_left_ymax,
                                wall_left_zmin, wall_left_zmax)
    wall_left.save(f'{output_dir}/wall_left.stl')
    print(f"생성됨: wall_left.stl")
    
    # === 2. 오른쪽 Wall 생성 ===
    wall_right_xmin = gap/2
    wall_right_xmax = gap/2 + T
    
    wall_right = create_box_stl(wall_right_xmin, wall_right_xmax,
                                 wall_left_ymin, wall_left_ymax,
                                 wall_left_zmin, wall_left_zmax)
    wall_right.save(f'{output_dir}/wall_right.stl')
    print(f"생성됨: wall_right.stl")
    
    # === 3. 왼쪽 핀 생성 (오른쪽 방향으로 돌출) ===
    # 핀은 y 방향으로 균등 배치
    total_fin_space = (num_fins_per_side - 1) * p + num_fins_per_side * t
    y_start = (H_wall - total_fin_space) / 2  # 중앙 정렬
    
    for layer in range(num_layers):
        z_center = layer * (W_wall / num_layers) + (W_wall / (2 * num_layers))
        
        for i in range(num_fins_per_side):
            y_pos = y_start + i * (t + p)
            
            fin_xmin = -gap/2
            fin_xmax = -gap/2 + L
            fin_ymin = y_pos
            fin_ymax = y_pos + t
            fin_zmin = z_center - t/2
            fin_zmax = z_center + t/2
            
            fin = create_box_stl(fin_xmin, fin_xmax,
                                fin_ymin, fin_ymax,
                                fin_zmin, fin_zmax)
            fin.save(f'{output_dir}/fin_left_layer{layer}_pin{i}.stl')
            print(f"생성됨: fin_left_layer{layer}_pin{i}.stl")
    
    # === 4. 오른쪽 핀 생성 (왼쪽 방향으로 돌출, 교차 배치) ===
    # 왼쪽 핀과 교차하도록 y 방향으로 offset
    y_start_right = y_start + (t + p) / 2
    
    for layer in range(num_layers):
        z_center = layer * (W_wall / num_layers) + (W_wall / (2 * num_layers))
        
        for i in range(num_fins_per_side):
            y_pos = y_start_right + i * (t + p)
            
            fin_xmin = gap/2 - L
            fin_xmax = gap/2
            fin_ymin = y_pos
            fin_ymax = y_pos + t
            fin_zmin = z_center - t/2
            fin_zmax = z_center + t/2
            
            fin = create_box_stl(fin_xmin, fin_xmax,
                                fin_ymin, fin_ymax,
                                fin_zmin, fin_zmax)
            fin.save(f'{output_dir}/fin_right_layer{layer}_pin{i}.stl')
            print(f"생성됨: fin_right_layer{layer}_pin{i}.stl")
    
    print(f"\n총 {2 + 2*num_fins_per_side*num_layers}개 STL 파일 생성 완료")
    print(f"저장 위치: {output_dir}/")


# ========== 메인 실행 부분 ==========
if __name__ == "__main__":
    # 파라미터 설정
    gap = 0.20          # wall 사이 공간 (m)
    T = 0.01            # wall 두께 (m)
    H_wall = 0.10       # wall 세로 길이 (m)
    W_wall = 0.10       # wall 가로 길이 (m)
    p = 0.02            # fin 사이 간격 (m)
    t = 0.005           # fin 두께 (m)
    L = 0.08            # fin 길이 (m)
    num_fins_per_side = 4   # 한쪽당 핀 개수
    num_layers = 3          # 레이어 개수
    
    print("=== 파라미터 ===")
    print(f"wall 간격 (gap): {gap} m")
    print(f"wall 두께 (T): {T} m")
    print(f"wall 크기: {H_wall} x {W_wall} m")
    print(f"fin 간격 (p): {p} m")
    print(f"fin 두께 (t): {t} m")
    print(f"fin 길이 (L): {L} m")
    print(f"한쪽당 핀 개수: {num_fins_per_side}")
    print(f"레이어 개수: {num_layers}\n")
    
    # STL 파일 생성
    generate_fin_wall_geometry(
        gap=gap,
        T=T,
        H_wall=H_wall,
        W_wall=W_wall,
        p=p,
        t=t,
        L=L,
        num_fins_per_side=num_fins_per_side,
        num_layers=num_layers
    )

