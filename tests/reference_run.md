# Reference run

RunID: 20261023T090951Z
Rssa: 1.1
Stage status: PASS

A run is passed when the stage status is PASS. That requires no FAIL row, and no NUMERICALLY_SENSITIVE row on SingularValues, ForecastH, AnalyticTarget, ScaleWolfram, ScaleRssa, PermuteWolfram, or PermuteRssa for the ordinary series CONST, D1, LINEAR, HARMONIC, TREND_HARMONIC, and NOISY. UNIDENTIFIABLE does not fail the run.

PASS 1938, FAIL 0, NUMERICALLY_SENSITIVE 32, UNIDENTIFIABLE 32, NOT_RUN 0

Maximum absolute difference on rows other than UNIDENTIFIABLE:

- Trajectory, Wolfram against the explicit Hankel matrix: 0
- KnownVector: 0
- SingularValues, Wolfram against Rssa: 9.208633855450898e-12
- LeftSubspace, projection residual: 5.521022831266834e-13
- RightSubspace, projection residual: 5.550810109533966e-13
- ForecastH, Wolfram against Rssa: 6.483702463810914e-14

UNIDENTIFIABLE marks LeftSubspace, RightSubspace, OperatorA, and OperatorC when the relative singular gap is at most 1e-8 or the rank exceeds the retained numerical rank. The subspace and those operators are then not identified by the gap criterion. ForecastH at the same rank is still compared.