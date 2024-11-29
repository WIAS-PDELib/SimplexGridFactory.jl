# Changes
## v2.3.0 Nov 28, 2024
- Allow for TetGen v1, v2
- Deprecate simplexgrid(Module...)

## v2.2.1 Nov 2, 2024
- Move repo to WIAS-PDELib org

## v2.2 June 25, 2024
- bregions! without args and kwargs takes all boundaries instead of none.

## v2.1 May 24, 2024
- Use simplexgrid(tetgenio) from ExtendableGrids
- Require Julia 1.9

## v2.0 March 21, 2024
- Breaking: moved BinnedPointList to ExtendableGrids
- Drop Julia 1.6 support, require Julia 1.9

## v1.0 Nov 15, 2023
Should be essentially non-breaking, but it is time to move to real semver
- Update to recent TetGen/Triangulate versions
- TetGen/Triangulate version check via weakdeps and compat (starting with 1.9)
- allow builderplot with Makie
- Formatting
