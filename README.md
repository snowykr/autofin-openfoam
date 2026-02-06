# autofin-openfoam

두 개의 wall(좌/우) 사이에, wall에서 돌출된 직육면체 fin(pin)들을 배치한 STL을 자동 생성하고,
해당 STL을 입력으로 OpenFOAM의 `blockMesh` + `surfaceFeatures`(또는 `surfaceFeatureExtract`) + `snappyHexMesh`까지(또는 dict 생성까지만) 자동화하는 프로젝트입니다.

이 프로젝트의 기본 meshing 설정은 fin-wall 접합부 모서리 보존을 위해 **explicit feature snapping이 필수**라고 가정합니다.
따라서 기본 워크플로우에서 `surfaceFeatures`(또는 `surfaceFeatureExtract`)를 `snappyHexMesh` 전에 실행해야 합니다.

## Requirements

- Python 3.11+
- `uv`
- OpenFOAM (meshing까지 실행하려면 `blockMesh`, `snappyHexMesh`가 PATH에 있어야 함)
  - 기본값은 feature snapping ON 이므로, 추가로 `surfaceFeatures`(OpenFOAM.org/Foundation) 또는 `surfaceFeatureExtract`(일부 배포판)가 필요합니다.
  - feature snapping을 끄려면 `--no-feature-snap`.

Python 의존성은 `pyproject.toml`에 정의되어 있으며, 기본은 `numpy-stl`입니다.

## Files

- `fin_wall_stl.py`
  - wall + fin STL 생성 로직
  - 출력: `constant/triSurface/solid_left.stl`, `constant/triSurface/solid_right.stl`
- `openfoam_automation.py`
  - STL 생성 + `system/blockMeshDict`, `system/snappyHexMeshDict`, `system/controlDict` 생성
  - multi-region CHT용으로 `constant/materialProperties`, `system/fvSchemes`, `system/fvSolution`도 생성
  - 기본값으로 `system/surfaceFeaturesDict`(또는 `system/surfaceFeatureExtractDict`)도 생성
  - `blockMesh` -> feature 추출 -> `snappyHexMesh -overwrite` -> `splitMeshRegions` 순서로 실행되어야 함
  - `Allmesh`(meshing) + `Allrun`(meshing + foamSetupCHT + foamMultiRun) 스크립트를 함께 생성

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
  - `background_cell_size`: m
  - `target_surface_cell_size`: m (`None`이면 자동으로 `max(t/2, 1e-6)` m 사용)
  - `feature_included_angle_deg`: degree (deg, feature 추출 임계각)
- Meshing (내부 설정값)
  - `max_local_cells`, `max_global_cells`, `min_refinement_cells`, `n_cells_between_levels`: 무차원(셀 개수/레벨 관련 정수)
  - `resolve_feature_angle_deg`: degree (deg)

즉, 예를 들어 `--t 0.005`는 fin 두께 5 mm(=0.005 m), `--background-cell-size 0.005`는 배경 셀 크기 5 mm를 의미합니다.

### Placement rule (checkerboard)

`y-z` 단면에서 `t x t` 사각 셀을 `p` 간격으로 배열하고, 셀마다 fin을 하나씩 배치합니다.
fin은 좌/우 wall에 번갈아 붙는 checkerboard 패턴이며, fin들은 직육면체입니다.

### Constraint

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

meshing(`blockMesh`, feature 추출, `snappyHexMesh`, `splitMeshRegions`) 실행 방법은 OpenFOAM을 어떻게 설치했는지에 따라 달라집니다.

- 로컬에 OpenFOAM이 설치되어 있고, `blockMesh`/`snappyHexMesh`가 PATH에 잡히는 경우
- macOS에서 `openfoam-dev-macos`(Docker)로 OpenFOAM을 쓰는 경우

```bash
# (A) 로컬 OpenFOAM 설치 (주로 Linux)
# OpenFOAM 환경 로드 후, python runner가 blockMesh/feature 추출/snappyHexMesh까지 호출
source /opt/openfoam*/etc/bashrc
uv run python openfoam_automation.py

# (B) macOS: openfoam-dev-macos (Docker)
# 호스트에서 case 파일만 생성하고, meshing은 컨테이너 안에서 실행
uv run python openfoam_automation.py --write-only
openfoam-dev-macos
./Allmesh

# (선택) multi-region CHT 셋업 + 실행(OpenFOAM.org dev 워크플로우)
./Allrun

# Allmesh 없이 수동으로 하면:
blockMesh
surfaceFeatures
snappyHexMesh -overwrite
splitMeshRegions -cellZones -overwrite -defaultRegionName fluid

# (선택) CHT 셋업 + 실행
foamSetupCHT
foamMultiRun
```

참고:
- zsh에서 `source /opt/openfoam*/etc/bashrc`가 `no matches found`로 실패하면, 해당 경로에 OpenFOAM이 설치되어 있지 않은 경우입니다. 이때는 (B)처럼 Docker를 쓰거나, 실제 설치 경로로 `bashrc`를 지정하세요.

### 4) 파라미터 변경 예시

```bash
uv run python openfoam_automation.py \
  --write-only \
  --gap 0.20 --T 0.01 --H-wall 0.10 --W-wall 0.10 \
  --p 0.02 --t 0.005 --L 0.12 \
  --background-cell-size 0.005 \
  --target-surface-cell-size 0.0025
```

feature snapping 파라미터 변경 예시:

```bash
uv run python openfoam_automation.py \
  --write-only \
  --feature-tool surfaceFeatures \
  --feature-included-angle-deg 150
```

feature snapping을 끄고(권장하지 않음) 예전처럼 `blockMesh` + `snappyHexMesh`만 쓰려면:

```bash
uv run python openfoam_automation.py --write-only --no-feature-snap
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

- `blockMesh` 배경 도메인은 `x=[-gap/2-T, +gap/2+T]`, `y=[0, H_wall]`, `z=[0, W_wall]`로 타이트하게 생성합니다.
  - 목적: 물리적으로 필요한 영역(좌 고체 + gap 유체 + 우 고체)만 포함하기 위함입니다.
- `snappyHexMesh`는 `insidePoints`(복수 점)로 fluid + solid 영역을 동시에 keep 하도록 생성됩니다.
  - `refinementSurfaces`에서 `cellZone solid_left/solid_right; mode inside;`를 사용해 좌/우 solid cellZone을 만들고,
    `splitMeshRegions -cellZones -overwrite -defaultRegionName fluid`로 region을 분리합니다.
- `insidePoints`는 각 keep 대상 영역(유체 1점 + 좌/우 wall 내부 1점씩)에 있어야 하며, STL 표면이나 배경 격자 경계에 너무 가까우면 실패/오동작할 수 있습니다.

CLI 전체 옵션은 아래로 확인할 수 있습니다.

```bash
uv run python openfoam_automation.py -h
```

## Outputs

기본적으로 case 디렉토리에 아래가 생성됩니다.

- `constant/triSurface/*.stl`
  - `solid_left.stl`, `solid_right.stl`
- `constant/materialProperties` (CHT용: foamSetupCHT 입력)
- `system/blockMeshDict`
- `system/snappyHexMeshDict`
- `system/controlDict`
- `system/fvSchemes`, `system/fvSolution`
  - 기본값으로 아래 중 하나가 생성됨
    - `system/surfaceFeaturesDict` (OpenFOAM.org/Foundation)
    - `system/surfaceFeatureExtractDict` (일부 배포판)
  - `Allmesh` (blockMesh -> feature 추출 -> snappyHexMesh -> splitMeshRegions)
  - `Allrun` (Allmesh -> foamSetupCHT -> foamMultiRun)

## Git notes

- `constant/triSurface/*.stl`은 자동 생성 산출물이므로 `.gitignore`에 포함되어 있습니다.
- `case.foam`은 ParaView용 로컬 더미 파일이며 `.gitignore`에 포함되어 있습니다.
- `system/blockMeshDict`, `system/snappyHexMeshDict`, `system/controlDict`는 스크립트로 자동 생성되지만, 기본 케이스 설정 공유를 위해 현재 저장소에서 추적(버전 관리)합니다.
  - 동일하게 `system/fvSchemes`, `system/fvSolution`, `constant/materialProperties`도 기본 템플릿으로 생성/추적합니다.
  - 파라미터 스윕 등으로 작업 트리를 더럽히고 싶지 않으면 `--case-dir`로 별도 디렉토리에 생성하세요.

## Manual OpenFOAM commands (if not using python runner)

```bash
./Allmesh

./Allrun

# 또는 수동 실행
blockMesh
surfaceFeatures
snappyHexMesh -overwrite
splitMeshRegions -cellZones -overwrite -defaultRegionName fluid

# (선택) CHT 셋업 + 실행
foamSetupCHT
foamMultiRun
```
