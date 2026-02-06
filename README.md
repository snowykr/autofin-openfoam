# autofin-openfoam

좌/우 wall + checkerboard fin + 가운데 gap 유체를 OpenFOAM `blockMesh` 기반으로
직접 생성하는 케이스 자동화 스크립트입니다.

이 저장소는 구조적(conformal) `blockMesh` 격자에 `cellZone`을 직접 부여한 뒤
`splitMeshRegions`로 region을 분리하는 방식을 사용합니다.

## Requirements

- Python 3.11+
- OpenFOAM (최소 `blockMesh`, `splitMeshRegions`)

## What It Generates

- `system/blockMeshDict`
  - 블록별 `cellZone`: `solid_left`, `fluid`, `solid_right`
- `system/controlDict`
- `system/fvSchemes`, `system/fvSolution`
- `Allmesh`
  - `blockMesh`
  - `splitMeshRegions -cellZones -overwrite -defaultRegionName fluid`
- `Allrun` (선택)
  - `paraFoam -touchAll`

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

- `0 < L <= gap` 필수
  - 반대 wall 관통 방지를 위한 제약
- `L <= gap/2`인 경우: 좌/우 fin이 x 방향으로 교차(중첩)되지 않음을 경고합니다.
- `L == gap`인 경우: 극한(full-span) 케이스임을 경고합니다.
- `gap > 0`, `T > 0`, `H_wall > 0`, `W_wall > 0`, `t > 0`, `p > 0`

## Usage

### 1) case 파일만 생성

```bash
python3 openfoam_automation.py --write-only
```

생성 전 기존 생성물을 정리하려면:

```bash
python3 openfoam_automation.py --write-only --clean-generated
```

### 2) 파라미터 지정 생성

```bash
python3 openfoam_automation.py \
  --write-only \
  --gap 0.30 --T 0.025 --H-wall 0.100 --W-wall 0.050 \
  --p 0.001 --t 0.001 --L 0.270 \
  --N-T 5 --N-gapA 3 --N-over 10 --N-t 5 --N-p 5
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

### 4) 시각화 준비 실행 (선택)

```bash
./Allrun
```

## Notes

- 배경 도메인은 의미론적으로 정확하게 고정됩니다:
  - `x=[-gap/2-T, +gap/2+T]`
  - `L` 경계값(`L=gap/2`, `L=gap`)에서는 0 두께 x 레이어를 자동 제거해 유효 블록만 생성합니다.
  - `y/z`는 `pitch=(t+p)`의 정수 배가 되도록 최근접(half-up) 반올림된 실폭으로 생성되며, 각 방향 최소 2 유닛을 보장합니다.
  - wall 바깥 유체 슬랩을 만들지 않습니다.
- 메쉬 밀도는 고정 분할 파라미터로 제어합니다.
  - `N_T`, `N_gapA`, `N_over` : x 방향 분할 수(경계값에서는 0 두께 레이어 제거로 활성 레이어 수가 5보다 작아질 수 있음)
  - `N_t`, `N_p` : y/z 방향 fin 두께/간격 구간 분할 수

## Visualization

```bash
touch case.foam
paraview --data="$(pwd)/case.foam"
```
