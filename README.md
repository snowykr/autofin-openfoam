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

### Placement rule (checkerboard)

`y-z` 단면에서 `t x t` 사각 셀을 `p` 간격으로 배열하고, 셀마다 fin을 하나씩 배치합니다.
fin은 좌/우 wall에 번갈아 붙는 checkerboard 패턴이며, fin들은 직육면체입니다.

### Constraint

fin이 중앙에서 서로 "교차(겹침)"하도록 하기 위해 아래 조건을 강제합니다.

- `L > gap/2` : 좌/우 fin이 gap 중앙에서 x방향으로 겹치도록
- `L <= gap` : 반대편 wall을 관통하지 않도록

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

OpenFOAM 환경을 먼저 로드한 뒤 실행합니다.

```bash
source /opt/openfoam*/etc/bashrc
uv run python openfoam_automation.py
```

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

## Manual OpenFOAM commands (if not using python runner)

```bash
blockMesh
snappyHexMesh -overwrite
```
