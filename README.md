# autofin-openfoam

두 개의 wall(좌/우) 사이에, wall에서 돌출된 직육면체 fin(pin)들을 배치한 STL을 자동 생성하고,
해당 STL을 입력으로 OpenFOAM의 `blockMesh` + `snappyHexMesh`까지(또는 dict 생성까지만) 자동화하는 프로젝트입니다.

## Requirements

- Python 3.11+
- `uv`
- OpenFOAM (meshing까지 실행하려면 `blockMesh`, `snappyHexMesh`가 PATH에 있어야 함)

Python 의존성은 `pyproject.toml`에 정의되어 있으며, 기본은 `numpy-stl`입니다.

## Files

- `fin_wall_stl.py`
  - wall + fin STL 생성 로직
  - 출력: `constant/triSurface/*.stl`
- `openfoam_automation.py`
  - STL 생성 + `system/blockMeshDict`, `system/snappyHexMeshDict`, `system/controlDict` 생성
  - 옵션으로 `blockMesh`, `snappyHexMesh -overwrite` 실행

## Geometry / Parameters

좌표계:

- x: 두 wall 사이 방향
- y: wall 세로 방향 (height)
- z: wall 가로 방향 (width)

### Geometry parameters

`openfoam_automation.py`와 `fin_wall_stl.py`에서 공통으로 사용하는 파라미터입니다.

- `gap` : 두 wall 사이 거리
- `T` : wall 두께
- `H_wall` : wall 높이 (y)
- `W_wall` : wall 폭 (z)
- `t` : fin 단면 두께 (현재 구현은 y 두께 = z 두께 = `t`)
- `p` : fin 사이 클리어런스(거리). y, z 방향 모두에 동일 적용
  - 인접 fin 표면-표면 간격 = `p`
  - 그리드 피치(중심 간격) = `t + p`
- `L` : fin 길이 (x 방향, wall에서 gap 내부로 돌출)

### Units (important)

이 프로젝트의 길이 관련 파라미터는 OpenFOAM 좌표계 단위와 동일하며, 현재 설정(`blockMeshDict`의 `scale 1;`)에서는 **모두 미터(m)** 입니다.

- Geometry
  - `gap`, `T`, `H_wall`, `W_wall`, `t`, `p`, `L`: m
- Meshing (CLI)
  - `margin`: m
  - `background_cell_size`: m
  - `target_surface_cell_size`: m (`None`이면 자동으로 `max(t/2, 1e-6)` m 사용)
- Meshing (내부 설정값)
  - `max_local_cells`, `max_global_cells`, `min_refinement_cells`, `n_cells_between_levels`: 무차원(셀 개수/레벨 관련 정수)
  - `resolve_feature_angle_deg`: degree (deg)

즉, 예를 들어 `--t 0.005`는 fin 두께 5 mm(=0.005 m), `--background-cell-size 0.005`는 배경 셀 크기 5 mm를 의미합니다.

### Placement rule (checkerboard)

`y-z` 단면에서 `t x t` 사각 셀을 `p` 간격으로 배열하고, 셀마다 fin을 하나씩 배치합니다.
fin은 좌/우 wall에 번갈아 붙는 checkerboard 패턴이며, fin들은 직육면체입니다.

### Constraint

fin이 중앙에서 서로 "교차(겹침)"하도록 하기 위해 아래 조건을 강제합니다.

- `L > gap/2` : 좌/우 fin이 gap 중앙에서 x방향으로 겹치도록
- `L <= gap` : 반대편 wall을 관통하지 않도록

추가 제약(현재 `fin_wall_stl.py` 구현):

- 좌/우 wall에 fin이 최소 1개씩 있어야 합니다. (즉, y-z 그리드 셀이 최소 2개 필요)

## Usage

### 1) Install python deps (uv)

```bash
uv sync
```

### 2) STL + dict 생성만 (OpenFOAM 실행 없음)

현재 디렉토리를 OpenFOAM case 디렉토리로 보고 파일을 생성합니다.

```bash
uv run python openfoam_automation.py --write-only
```

다른 디렉토리에 생성하려면:

```bash
uv run python openfoam_automation.py --case-dir ./case1 --write-only
```

### 3) meshing까지 자동 실행

meshing(`blockMesh`, `snappyHexMesh`) 실행 방법은 OpenFOAM을 어떻게 설치했는지에 따라 달라집니다.

- 로컬에 OpenFOAM이 설치되어 있고, `blockMesh`/`snappyHexMesh`가 PATH에 잡히는 경우
- macOS에서 `openfoam-dev-macos`(Docker)로 OpenFOAM을 쓰는 경우

```bash
# (A) 로컬 OpenFOAM 설치 (주로 Linux)
# OpenFOAM 환경 로드 후, python runner가 blockMesh/snappyHexMesh까지 호출
source /opt/openfoam*/etc/bashrc
uv run python openfoam_automation.py

# (B) macOS: openfoam-dev-macos (Docker)
# 호스트에서 case 파일만 생성하고, meshing은 컨테이너 안에서 실행
uv run python openfoam_automation.py --write-only
openfoam-dev-macos
blockMesh
snappyHexMesh -overwrite
```

참고:
- zsh에서 `source /opt/openfoam*/etc/bashrc`가 `no matches found`로 실패하면, 해당 경로에 OpenFOAM이 설치되어 있지 않은 경우입니다. 이때는 (B)처럼 Docker를 쓰거나, 실제 설치 경로로 `bashrc`를 지정하세요.

### 4) 파라미터 변경 예시

```bash
uv run python openfoam_automation.py \
  --write-only \
  --gap 0.20 --T 0.01 --H-wall 0.10 --W-wall 0.10 \
  --p 0.02 --t 0.005 --L 0.12 \
  --margin 0.02 \
  --background-cell-size 0.005 \
  --target-surface-cell-size 0.0025
```

### 5) 시각화 (host ParaView)

macOS/Homebrew ParaView에서 `paraview --data=...`로 열 때는 `case.foam` 더미 파일을 쓰는 것이 가장 호환성이 좋습니다.

```bash
touch case.foam
paraview --data="$(pwd)/case.foam"
```

참고:
- 위 명령의 `touch case.foam`은 파일을 "생성"만 합니다(창이 자동으로 뜨지 않음). 실제 실행은 `paraview --data=...`입니다.
- 일부 환경에서는 `case.OpenFOAM`이 ParaView 기본 reader와 매칭되지 않을 수 있습니다.
- 대안으로, OpenFOAM 컨테이너에서 `foamToVTK -constant -allPatches` 후 `VTK/*.vtk`를 열어도 됩니다.

## Meshing notes

- `blockMesh` 배경 도메인은 기본적으로 `y=[0, H_wall]`, `z=[0, W_wall]`로 타이트하게 생성합니다.
  - 목적: wall slab의 모서리로 "바깥 영역"이 돌아나가면서 내부(gap)와 연결되는 것을 방지하고,
    `snappyHexMesh`의 region 선택(`locationInMesh`)이 안정적으로 gap 내부만 남기도록 하기 위함입니다.
- `locationInMesh`는 "남길 영역" 안에 반드시 있어야 하며, STL 표면이나 배경 격자 경계에 너무 가까우면 실패/오동작할 수 있습니다.
  - `openfoam_automation.py`가 자동으로 안전한 점을 찾지만, 극단적인 파라미터(`p=0` + 단면 완전 타일링 + `L≈gap`)에서는 유체 영역이 거의/전혀 없어 실패할 수 있습니다.

CLI 전체 옵션은 아래로 확인할 수 있습니다.

```bash
uv run python openfoam_automation.py -h
```

## Outputs

기본적으로 case 디렉토리에 아래가 생성됩니다.

- `constant/triSurface/*.stl`
  - `wall_left.stl`, `wall_right.stl`
  - `fin_left_y{iy}_z{iz}.stl`, `fin_right_y{iy}_z{iz}.stl`
- `system/blockMeshDict`
- `system/snappyHexMeshDict`
- `system/controlDict`

## Git notes

- `constant/triSurface/*.stl`은 자동 생성 산출물이므로 `.gitignore`에 포함되어 있습니다.
- `case.foam`은 ParaView용 로컬 더미 파일이며 `.gitignore`에 포함되어 있습니다.
- `system/blockMeshDict`, `system/snappyHexMeshDict`, `system/controlDict`는 스크립트로 자동 생성되지만, 기본 케이스 설정 공유를 위해 현재 저장소에서 추적(버전 관리)합니다.
  - 파라미터 스윕 등으로 작업 트리를 더럽히고 싶지 않으면 `--case-dir`로 별도 디렉토리에 생성하세요.

## Manual OpenFOAM commands (if not using python runner)

```bash
blockMesh
snappyHexMesh -overwrite
```
