# autofin-openfoam

좌/우 wall + checkerboard fin + 가운데 gap 유체를 OpenFOAM `blockMesh` 기반으로
직접 생성하는 케이스 자동화 스크립트입니다.

이 저장소의 현재 기본 전략은 **snappyHexMesh를 사용하지 않고**, 구조적(conformal) `blockMesh`
격자에 `cellZone`을 직접 부여한 뒤 `splitMeshRegions`로 region을 분리하는 방식입니다.

## Requirements

- Python 3.11+
- `uv`
- OpenFOAM (최소 `blockMesh`, `splitMeshRegions`)

## What It Generates

- `system/blockMeshDict`
  - 블록별 `cellZone`: `solid_left`, `fluid`, `solid_right`
- `system/controlDict`
- `system/fvSchemes`, `system/fvSolution`
- `constant/materialProperties`
- `Allmesh`
  - `blockMesh`
  - `splitMeshRegions -cellZones -overwrite -defaultRegionName fluid`
- `Allrun` (선택)
  - `./Allmesh`
  - `foamSetupCHT`
  - `foamMultiRun`

## Geometry Parameters

- `gap`: 두 wall 사이 거리 (x)
- `T`: wall 두께
- `H_wall`: wall 높이 (y)
- `W_wall`: wall 폭 (z)
- `t`: fin 단면 두께 (y, z 동일)
- `p`: fin 간 클리어런스 (y, z 동일)
- `L`: fin 길이 (x, wall에서 gap 내부로 돌출)

모든 길이 단위는 m (`blockMeshDict`의 `scale 1;`).

## Constraints (blockMesh-only)

- `L <= gap` 필수
  - 반대 wall 관통 방지를 위한 제약
- `gap > 0`, `T > 0`, `H_wall > 0`, `W_wall > 0`

## Usage

### 1) case 파일만 생성

```bash
uv run python openfoam_automation.py --write-only
```

### 2) 파라미터 지정 생성

```bash
uv run python openfoam_automation.py \
  --write-only \
  --gap 0.30 --T 0.025 --H-wall 0.100 --W-wall 0.050 \
  --p 0.001 --t 0.001 --L 0.120 \
  --background-cell-size 0.002
```

### 3) meshing 실행

```bash
./Allmesh
```

수동 실행:

```bash
blockMesh
splitMeshRegions -cellZones -overwrite -defaultRegionName fluid
```

### 4) CHT 실행 (선택)

```bash
./Allrun
```

주의: `foamSetupCHT`는 케이스 루트에 `templates/` 트리를 요구합니다.

## Notes

- 배경 도메인은 의미론적으로 정확하게 고정됩니다:
  - `x=[-gap/2-T, +gap/2+T]`, `y=[0, H_wall]`, `z=[0, W_wall]`
  - wall 바깥 유체 슬랩을 만들지 않습니다.
- `background_cell_size`는 블록별 셀 수를 결정합니다.
  - 너무 크면 간극 해상도가 낮아집니다.

## Visualization

```bash
touch case.foam
paraview --data="$(pwd)/case.foam"
```
