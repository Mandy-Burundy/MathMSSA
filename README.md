# MathMSSA

One-step MSSA row forecast in the Wolfram Language. `X` holds exogenous series. This is the forecast kernel, not a full SSA library: there is no grouping UI, no reconstruction pipeline, and no price conversion.

`MathMSSA[y, X, L, ranks]` returns the next target value as `ForecastH`, in the same numeric scale as the series you pass in. The package does not turn that number into a price.

The same call in Rssa 1.1 is:

```r
ssa(channels, L = L, kind = "mssa", svd.method = "svd", neig = max(ranks) + 1)
rforecast(s, groups = list(1:r), len = 1, direction = "row", base = "original", only.new = TRUE)
```

`channels` is the target followed by the columns of `X`. On the synthetic comparison in [`tests/reference_run.md`](tests/reference_run.md), run `20261023T090951Z` has stage status **PASS**. Against Rssa 1.1 the trajectory matrix and the known vector match exactly (absolute difference 0). The largest absolute difference is 9.21×10⁻¹² for singular values, 5.55×10⁻¹³ for the subspace projection residual, and 6.48×10⁻¹⁴ for `ForecastH`.

A run is passed when the stage status is `PASS`: no row is `FAIL`, and none of the ordinary series `CONST`, `D1`, `LINEAR`, `HARMONIC`, `TREND_HARMONIC`, `NOISY` is `NUMERICALLY_SENSITIVE` on singular values, `ForecastH`, the analytic target, or the scale and permutation checks. `UNIDENTIFIABLE` does not fail the run. That status is used for `LeftSubspace`, `RightSubspace`, `OperatorA`, and `OperatorC` when the relative singular gap is at most 10⁻⁸ or the rank exceeds the retained numerical rank. The subspace and the operators built from one particular basis are then not identified, so those rows are neither a match nor a failure. `ForecastH` at the same rank is still compared. A sensitive `ForecastH` becomes `FAIL` when its absolute difference exceeds `10^4 * cond * eps * max(scale, 1)`.

This run: 1938 `PASS`, 0 `FAIL`, 32 `NUMERICALLY_SENSITIVE`, 32 `UNIDENTIFIABLE`. The full counts and maxima are in `tests/reference_run.md`. Directories under `tests/runs/` are gitignored.

## Call

`y` has length `N`. `X` has shape `(N + 1) × d` with `d ≥ 1`. The extra row of `X` is the exogenous observation at the forecast time. `L` is an integer, `2 ≤ L ≤ N − 1`. `ranks` is a strictly increasing list of positive integers, and the largest rank is at most `min(L, K)`. `K` is the number of columns of the stacked trajectory matrix: the target block contributes `N − L + 1` columns and each of the `d` exogenous blocks contributes `N − L + 2`, so `K = (N − L + 1) + d (N − L + 2)`.

From the repository root:

```wl
Get["Kernel/MathMSSA.wl"];
y = Sin[2. Pi Range[64]/17.];
x = Transpose[{0.7 Sin[2. Pi Range[65]/17. + 0.4]}];
MathMSSA`MathMSSA[y, x, 20, {2}]["ByRank"][2]["ForecastH"]
```

`Get["Kernel/MathMSSA.wl"]` resolves that relative path from the current directory, so it works when the current directory is the repository root. From another directory, pass the full path of `Kernel/MathMSSA.wl`.

To install the paclet, build it and install the archive. `PacletBuild` writes `build/MathMSSA-1.0.0.paclet`. The `build/` directory is gitignored.

```wl
Needs["PacletTools`"];
built = PacletBuild["/full/path/to/MathMSSA"];
PacletInstall[built["PacletArchive"]];
Needs["MathMSSA`"];
```

Logs stay off unless you set `` MathMSSA`$MathMSSALogging = True `` before the call. `examples/BasicForecast.wl` is the same call as the snippet above.

The result is `<|"Decomposition" -> ..., "ByRank" -> <|r -> ...|>|>`. For each requested rank, `ByRank[r]["ForecastH"]` is the target forecast. A bad input or a failed decomposition returns `Failure[...]`.

Several ranks of one window share one truncated SVD, cut at `Max[ranks]`. The forecast uses the right singular vectors `V`:

- each channel becomes a Hankel block of `L` rows; the target block has `N − L + 1` columns, and each auxiliary block has `N − L + 2`;
- `S` is the rows of `V` at the end of each block, `V̄` is the other rows;
- `A = I − S Sᵀ`, `C = S V̄ᵀ`, `B = LinearSolve[A, C]`;
- `z` is the tail of every channel after its first `L` samples;
- `ForecastH` is the first component of `B.z`.

There is no condition-number cutoff. `LinearSolve[A, C]` is always attempted. `ByRank[r]["RowOperatorConditionNumber"]` is returned in either case, and it is `Infinity` when the smallest singular value of `A` is zero. If `LinearSolve` fails, `ForecastH` is `Missing["NotComputed"]` and `ByRank[r]["TechnicalFailure"]` holds a `Failure` with tag `MathMSSALinearSolveFailed`. The condition number computed before the solve remains on that result.

`Decomposition["MSSASingularValues"]` holds the retained singular values. The factors `U` and `V` are not returned. `Decomposition["MSSARetainedNumericalRank"]` counts values above the cutoff inside that retained spectrum, so it cannot exceed `Max[ranks]`. It is not the rank of the full trajectory. `ByRank[r]` also carries the condition number of `A` and the relative residual of `A.B − C`.

## Parallel windows

A single `MathMSSA` call is one window. Independent windows do not share state, so each one can go to its own subkernel. `MathMSSAParallel[{job, ...}]` does that with `ParallelTable` (`Method -> "FinestGrained"`, `DistributedContexts -> None`). A job is `<|"y" -> ..., "X" -> ..., "L" -> ..., "Ranks" -> ...|>`. Ranks inside a job stay on that subkernel, because they share the decomposition.

The return value is `<|"ExecutionMode" -> "PARALLEL" | "SERIAL", "Results" -> { ... }|>`. `SERIAL` means the jobs ran on the main kernel because no subkernel started or the worker failed to load `MathMSSA`.

From the repository root:

```wl
Get["Kernel/MathMSSAParallel.wl"];
MathMSSAParallel`MathMSSAParallel[{
  <|"y" -> y, "X" -> x, "L" -> 20, "Ranks" -> {1, 2}|>,
  <|"y" -> y, "X" -> x, "L" -> 40, "Ranks" -> {1, 2}|>
}]
```

`examples/ParallelPoints.wl` sends eight windows and prints `ExecutionMode`. The subkernels load `Kernel/MathMSSA.wl` themselves.

## Check against Rssa

The check needs R and two packages you install yourself: Rssa exactly 1.1, and jsonlite. The test does not install them. On the CRAN package page, Rssa 1.1 is the published version, with `Date/Publication` 2024-09-05 00:00:02 UTC. The `Packaged` field is 2024-09-04 20:53:58 UTC, which is why some documentation pages show 2024-09-04. A later upload to CRAN would make `install.packages("Rssa")` install a different version. Pin 1.1:

```r
install.packages(c("jsonlite", "remotes"))
remotes::install_version("Rssa", "1.1")
packageVersion("Rssa")
```

`packageVersion("Rssa")` must print `1.1`. Any other version is a different reference, and the test stops.

`RSCRIPT` may point at `Rscript`. Otherwise the test searches `PATH`, including `Rscript.exe` on Windows, and then `/usr/local/bin`, `/opt/homebrew/bin`, and `/usr/bin`.

```bash
wolframscript -file tests/MSSA_RSSA_STAGE_A.wl
```

Stage A is this synthetic comparison. It builds constants, lines, sines, a noisy copy, a one-auxiliary series (`d = 1`), a near-collinear family, and a few boundary inputs. It writes `tests/reference_run.md` and `tests/runs/<RunID>/` with `report.md`, `comparisons.csv`, and `manifest.json`. The run directories are gitignored. The summary file is not. Tolerances are fixed in the test: singular values `1e-10`, `ForecastH` absolute `1e-10` and relative `1e-8`, subspace residual `1e-7`, solve residual `1e-10`, relative singular gap `1e-8`.

The R script calls two non-exported Rssa 1.1 helpers, `.get.or.create.trajmat.mssa` and `.F`, and the exported `calc.v`. Those names are part of Rssa 1.1, not a public API, so the comparison is tied to that version.

## License

This package is released under the MIT License. See `LICENSE`. Rssa is distributed under the GPL. The comparison runs it as a separate R process and does not link it into this package.

## References

- Golyandina, N., Korobeynikov, A., and Zhigljavsky, A. (2018). *Singular Spectrum Analysis with R*. Springer. <https://doi.org/10.1007/978-3-662-57380-8>
- Golyandina, N. and Stepanov, D. (2005). SSA-based approaches to analysis and forecast of multidimensional time series. In *Proceedings of the 5th St. Petersburg Workshop on Simulation*, St. Petersburg State University, 293–298.
- Golyandina, N., Korobeynikov, A., Shlemov, A., and Usevich, K. (2015). Multivariate and 2D extensions of singular spectrum analysis with the Rssa package. *Journal of Statistical Software*, 67(2), 1–78. <https://doi.org/10.18637/jss.v067.i02>

## Files

| File | Role |
|---|---|
| `Kernel/MathMSSA.wl` | One window |
| `Kernel/MathMSSAParallel.wl` | Independent windows on subkernels |
| `examples/` | Short calls, not part of the paclet archive |
| `tests/MSSA_RSSA_STAGE_A.wl` | Synthetic comparison with Rssa |
| `tests/rssa_forecast.R` | The Rssa side of that comparison |
| `tests/reference_run.md` | Counts and maximum differences from the latest run |

Prepared market series and the code that turns `ForecastH` into a price are not in this folder.
