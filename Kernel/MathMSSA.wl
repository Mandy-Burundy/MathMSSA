(*
  MathMSSA.wl

  Pure MSSA row-forecast kernel in the H scale.
  Public contract:
    MathMSSA`MathMSSA[y, X, L, ranks]

  The kernel accepts y, X, L, and ranks. Configuration, dates, models, roles,
  causality, prices, UI, parallelism, and cache stay outside this file.
*)

BeginPackage["MathMSSA`"]

MathMSSA::usage =
  "MathMSSA[y, X, L, ranks] computes a one-step flat-layout row-forecast MSSA result in the prepared H scale for one point and one window. It returns Decomposition and ByRank, or Failure[...] for a common input/trajectory/SVD failure.";

$MathMSSALogging::usage =
  "$MathMSSALogging controls concise Print logs from MathMSSA. It stays False unless a caller sets it to True before MathMSSA.";

If[!ValueQ[$MathMSSALogging], $MathMSSALogging = False];

Begin["`Private`"]

$finiteTolerance = 1.*^-12;
$rankRelativeTolerance = 1.*^-10;

ClearAll[
  makeFailure,
  finiteRealQ, finiteRealVectorQ, finiteRealMatrixQ,
  toPackedRealVector, toPackedRealMatrix,
  rankGridQ, safeNorm2, numericalRankFromSingularValues,
  mssaLog,
  rankDiagnosticsTemplate, makeRankResult, makeRankFailure,
  buildHankelBlock, buildMSSATrajectory, decomposeTrajectory,
  operatorDiagnostics, evaluateRank, MathMSSA
];

mssaLog[stage_String, details_: <||>] :=
  If[
    TrueQ[MathMSSA`$MathMSSALogging],
    Print["[MathMSSA] ", stage, " | ", InputForm[details]]
  ];

makeFailure[tag_String, stage_String, details_: <||>] :=
  Module[{failure},
    mssaLog["FAILURE", <|"Tag" -> tag, "Stage" -> stage, "Details" -> details|>];
    failure = Failure[tag, <|"Stage" -> stage, "Details" -> details|>];
    failure
  ];

finiteRealQ[x_] := Quiet[
  Check[
    With[{value = N[x, MachinePrecision]},
      MachineNumberQ[value] && TrueQ[Im[value] == 0.]
    ],
    False
  ]
];

finiteRealVectorQ[value_] :=
  ListQ[value] && VectorQ[value, finiteRealQ];

finiteRealMatrixQ[value_] :=
  ListQ[value] && MatrixQ[value, finiteRealQ];

toPackedRealVector[value_List] :=
  Developer`ToPackedArray[N[Re /@ value, MachinePrecision]];

toPackedRealMatrix[value_List] :=
  Developer`ToPackedArray[N[Map[Re, value, {2}], MachinePrecision]];

rankGridQ[ranks_] :=
  ListQ[ranks] && ranks =!= {} &&
    VectorQ[ranks, IntegerQ[#] && # >= 1 &] &&
    Length[DeleteDuplicates[ranks]] == Length[ranks] && OrderedQ[ranks];

safeNorm2[value_] := Quiet@Check[N[Norm[value, 2], MachinePrecision], $Failed];

numericalRankFromSingularValues[singularValues_List] := Module[{cut},
  If[singularValues === {}, Return[0]];
  cut = Max[$finiteTolerance, $rankRelativeTolerance First[singularValues]];
  Count[singularValues, value_ /; TrueQ[value > cut]]
];

rankDiagnosticsTemplate[] := <|
  "RowOperatorMinEigenvalue" -> Missing["NotComputed"],
  "RowOperatorConditionNumber" -> Missing["NotComputed"],
  "RowOperatorReciprocalConditionNumber" -> Missing["NotComputed"],
  "RowOperatorSymmetryRelativeError" -> Missing["NotComputed"],
  "RowOperatorNorm" -> Missing["NotComputed"],
  "RowOperatorSolveRelativeResidual" -> Missing["NotComputed"]
|>;

makeRankResult[forecastH_, diagnostics_Association, technicalFailure_] := <|
  "ForecastH" -> forecastH,
  "RowOperatorMinEigenvalue" ->
    Lookup[diagnostics, "RowOperatorMinEigenvalue", Missing["NotComputed"]],
  "RowOperatorConditionNumber" ->
    Lookup[diagnostics, "RowOperatorConditionNumber", Missing["NotComputed"]],
  "RowOperatorReciprocalConditionNumber" ->
    Lookup[diagnostics, "RowOperatorReciprocalConditionNumber", Missing["NotComputed"]],
  "RowOperatorSymmetryRelativeError" ->
    Lookup[diagnostics, "RowOperatorSymmetryRelativeError", Missing["NotComputed"]],
  "RowOperatorNorm" ->
    Lookup[diagnostics, "RowOperatorNorm", Missing["NotComputed"]],
  "RowOperatorSolveRelativeResidual" ->
    Lookup[diagnostics, "RowOperatorSolveRelativeResidual", Missing["NotComputed"]],
  "TechnicalFailure" -> technicalFailure
|>;

makeRankFailure[diagnostics_Association, tag_String, stage_String, details_: <||>] :=
  makeRankResult[
    Missing["NotComputed"],
    diagnostics,
    makeFailure[tag, stage, details]
  ];

buildHankelBlock[series_List, window_Integer] :=
  Developer`ToPackedArray[
    Table[
      series[[row + column - 1]],
      {row, window},
      {column, Length[series] - window + 1}
    ]
  ];

buildMSSATrajectory[target_List, inputs_List, window_Integer] := Module[
  {channels, blocks, widths, endpoints, trajectory, knownRowVector},
  channels = Join[{target}, Transpose[inputs]];
  blocks = buildHankelBlock[#, window] & /@ channels;
  widths = Last[Dimensions[#]] & /@ blocks;
  endpoints = Accumulate[widths];
  trajectory = toPackedRealMatrix[Join[Sequence @@ blocks, 2]];
  knownRowVector = toPackedRealVector[Join @@ (Drop[#, window] & /@ channels)];
  <|
    "ChannelCount" -> Length[channels],
    "BlockWidths" -> widths,
    "EndpointIndices" -> endpoints,
    "TrajectoryColumns" -> Total[widths],
    "TrajectoryMatrix" -> trajectory,
    "KnownRowVector" -> knownRowVector
  |>
];

decomposeTrajectory[trajectory_?MatrixQ, ranks_List] := Module[
  {dimensions, maximumRank, result, u, w, v, singularValues, numericalRank},
  dimensions = Dimensions[trajectory];
  maximumRank = Last[ranks];
  mssaLog[
    "SVD_START",
    <|"Dimensions" -> dimensions, "MaxRank" -> maximumRank|>
  ];
  result = Quiet@Check[
    SingularValueDecomposition[N[trajectory, MachinePrecision], maximumRank],
    $Failed
  ];
  If[result === $Failed || !MatchQ[result, {_, _, _}],
    Return@makeFailure[
      "MathMSSASVDFailed",
      "SVD",
      <|"Dimensions" -> dimensions, "MaxRank" -> maximumRank|>
    ]
  ];
  {u, w, v} = result;
  singularValues = Quiet@Check[N[Diagonal[w], MachinePrecision], $Failed];
  If[
    !finiteRealMatrixQ[u] || !finiteRealMatrixQ[v] ||
      !finiteRealVectorQ[singularValues] ||
      Dimensions[u] =!= {First[dimensions], maximumRank} ||
      Dimensions[v] =!= {Last[dimensions], maximumRank} ||
      Length[singularValues] =!= maximumRank,
    Return@makeFailure[
      "MathMSSASVDNonFinite",
      "SVD",
      <|"Dimensions" -> dimensions, "MaxRank" -> maximumRank|>
    ]
  ];
  numericalRank = numericalRankFromSingularValues[singularValues];
  mssaLog[
    "SVD_DONE",
    <|
      "MaxRank" -> maximumRank,
      "NumericalRank" -> numericalRank,
      "LargestSingularValue" -> First[singularValues],
      "SmallestRetainedSingularValue" -> Last[singularValues]
    |>
  ];
  <|
    "MaxRank" -> maximumRank,
    "NumericalRank" -> numericalRank,
    "U" -> toPackedRealMatrix[u],
    "SingularValues" -> toPackedRealVector[singularValues],
    "V" -> toPackedRealMatrix[v]
  |>
];

operatorDiagnostics[a_?MatrixQ, c_?MatrixQ] := Module[
  {diagnostics, aNorm, symmetry, eigenvalues, singularValues, condition, reciprocal},
  diagnostics = rankDiagnosticsTemplate[];
  If[!finiteRealMatrixQ[a] || !finiteRealMatrixQ[c],
    Return@<|
      "Diagnostics" -> diagnostics,
      "TechnicalFailure" -> makeFailure["MathMSSANonFiniteOperator", "ROW_OPERATOR"]
    |>
  ];
  mssaLog["ROW_OPERATOR_NORM_START"];
  aNorm = safeNorm2[a];
  If[!finiteRealQ[aNorm],
    Return@<|
      "Diagnostics" -> diagnostics,
      "TechnicalFailure" -> makeFailure["MathMSSAOperatorNormFailed", "ROW_OPERATOR"]
    |>
  ];
  mssaLog["ROW_OPERATOR_NORM_DONE", <|"ANorm" -> aNorm|>];
  symmetry = Quiet@Check[
    N[safeNorm2[a - Transpose[a]]/Max[aNorm, $finiteTolerance], MachinePrecision],
    $Failed
  ];
  If[!finiteRealQ[symmetry],
    Return@<|
      "Diagnostics" -> diagnostics,
      "TechnicalFailure" -> makeFailure["MathMSSASymmetryDiagnosticFailed", "ROW_OPERATOR"]
    |>
  ];
  AssociateTo[diagnostics, "RowOperatorSymmetryRelativeError" -> symmetry];
  mssaLog[
    "ROW_OPERATOR_SYMMETRY_DONE",
    <|"SymmetryRelativeError" -> symmetry|>
  ];
  mssaLog["ROW_OPERATOR_EIGENVALUES_START"];
  eigenvalues = Quiet@Check[
    N[Eigenvalues[(a + Transpose[a])/2.0], MachinePrecision],
    $Failed
  ];
  If[!finiteRealVectorQ[eigenvalues] || eigenvalues === {},
    Return@<|
      "Diagnostics" -> diagnostics,
      "TechnicalFailure" -> makeFailure["MathMSSAEigenvalueDiagnosticFailed", "ROW_OPERATOR"]
    |>
  ];
  AssociateTo[diagnostics, "RowOperatorMinEigenvalue" -> Min[eigenvalues]];
  mssaLog[
    "ROW_OPERATOR_EIGENVALUES_DONE",
    <|"MinEigenvalue" -> diagnostics["RowOperatorMinEigenvalue"]|>
  ];
  mssaLog["ROW_OPERATOR_SINGULAR_VALUES_START"];
  singularValues = Quiet@Check[N[SingularValueList[a], MachinePrecision], $Failed];
  If[!finiteRealVectorQ[singularValues] || singularValues === {},
    Return@<|
      "Diagnostics" -> diagnostics,
      "TechnicalFailure" -> makeFailure["MathMSSAConditionDiagnosticFailed", "ROW_OPERATOR"]
    |>
  ];
  mssaLog[
    "ROW_OPERATOR_SINGULAR_VALUES_DONE",
    <|
      "LargestSingularValue" -> First[singularValues],
      "SmallestSingularValue" -> Last[singularValues]
    |>
  ];
  condition = If[
    TrueQ[Last[singularValues] == 0.],
    Infinity,
    N[First[singularValues]/Last[singularValues], MachinePrecision]
  ];
  reciprocal = If[condition === Infinity, 0.0, N[1.0/condition, MachinePrecision]];
  If[!finiteRealQ[condition] && condition =!= Infinity,
    Return@<|
      "Diagnostics" -> diagnostics,
      "TechnicalFailure" -> makeFailure["MathMSSAConditionDiagnosticFailed", "ROW_OPERATOR"]
    |>
  ];
  If[!finiteRealQ[reciprocal],
    Return@<|
      "Diagnostics" -> diagnostics,
      "TechnicalFailure" -> makeFailure["MathMSSAReciprocalConditionFailed", "ROW_OPERATOR"]
    |>
  ];
  diagnostics = Join[
    diagnostics,
    <|
      "RowOperatorConditionNumber" -> condition,
      "RowOperatorReciprocalConditionNumber" -> reciprocal
    |>
  ];
  mssaLog[
    "ROW_OPERATOR_DIAGNOSTICS_DONE",
    <|
      "ConditionNumber" -> condition,
      "ReciprocalConditionNumber" -> reciprocal
    |>
  ];
  <|"Diagnostics" -> diagnostics, "TechnicalFailure" -> None|>
];

evaluateRank[svd_Association, endpoints_List, knownRowVector_List, rank_Integer] := Module[
  {v, trajectoryColumns, channelCount, complement, s, vBar, a, c,
   outcome, diagnostics, b, aNorm, bNorm, cNorm, denominator, residual, eta},
  mssaLog["RANK_START", <|"Rank" -> rank|>];
  v = svd["V"];
  trajectoryColumns = Length[v];
  channelCount = Length[endpoints];
  complement = Complement[Range[trajectoryColumns], endpoints];
  If[Length[knownRowVector] =!= Length[complement],
    Return@makeRankFailure[
      rankDiagnosticsTemplate[],
      "MathMSSAKnownRowLengthMismatch",
      "ROW_FORECAST",
      <|
        "Rank" -> rank,
        "Expected" -> Length[complement],
        "Actual" -> Length[knownRowVector]
      |>
    ]
  ];
  s = Quiet@Check[
    toPackedRealMatrix[v[[endpoints, 1 ;; rank]]],
    $Failed
  ];
  vBar = Quiet@Check[
    toPackedRealMatrix[v[[complement, 1 ;; rank]]],
    $Failed
  ];
  If[s === $Failed || vBar === $Failed,
    Return@makeRankFailure[
      rankDiagnosticsTemplate[],
      "MathMSSASubspaceExtractionFailed",
      "ROW_OPERATOR",
      <|"Rank" -> rank|>
    ]
  ];
  mssaLog[
    "RANK_SUBSPACE_READY",
    <|
      "Rank" -> rank,
      "SDimensions" -> Dimensions[s],
      "VBarDimensions" -> Dimensions[vBar]
    |>
  ];
  a = Quiet@Check[
    toPackedRealMatrix[IdentityMatrix[channelCount] - s . Transpose[s]],
    $Failed
  ];
  c = Quiet@Check[
    toPackedRealMatrix[s . Transpose[vBar]],
    $Failed
  ];
  If[a === $Failed || c === $Failed,
    Return@makeRankFailure[
      rankDiagnosticsTemplate[],
      "MathMSSAOperatorConstructionFailed",
      "ROW_OPERATOR",
      <|"Rank" -> rank|>
    ]
  ];
  mssaLog[
    "RANK_OPERATOR_READY",
    <|
      "Rank" -> rank,
      "ADimensions" -> Dimensions[a],
      "CDimensions" -> Dimensions[c]
    |>
  ];
  mssaLog["RANK_DIAGNOSTICS_START", <|"Rank" -> rank|>];
  outcome = operatorDiagnostics[a, c];
  diagnostics = outcome["Diagnostics"];
  If[outcome["TechnicalFailure"] =!= None,
    mssaLog["RANK_OPERATOR_FAILED", <|"Rank" -> rank|>];
    Return@makeRankResult[
      Missing["NotComputed"],
      diagnostics,
      outcome["TechnicalFailure"]
    ]
  ];
  mssaLog[
    "RANK_OPERATOR",
    <|
      "Rank" -> rank,
      "MinEigenvalue" -> diagnostics["RowOperatorMinEigenvalue"],
      "ConditionNumber" -> diagnostics["RowOperatorConditionNumber"],
      "ReciprocalConditionNumber" -> diagnostics["RowOperatorReciprocalConditionNumber"],
      "SymmetryRelativeError" -> diagnostics["RowOperatorSymmetryRelativeError"]
    |>
  ];
  mssaLog["RANK_LINEAR_SOLVE_START", <|"Rank" -> rank|>];
  b = Quiet@Check[LinearSolve[a, c], $Failed];
  If[b === $Failed || !finiteRealMatrixQ[b],
    Return@makeRankFailure[
      diagnostics,
      "MathMSSALinearSolveFailed",
      "ROW_OPERATOR",
      <|"Rank" -> rank|>
    ]
  ];
  b = toPackedRealMatrix[b];
  mssaLog[
    "RANK_LINEAR_SOLVE_DONE",
    <|"Rank" -> rank, "BDimensions" -> Dimensions[b]|>
  ];
  aNorm = safeNorm2[a];
  bNorm = safeNorm2[b];
  cNorm = safeNorm2[c];
  If[!finiteRealQ[aNorm] || !finiteRealQ[bNorm] || !finiteRealQ[cNorm],
    Return@makeRankFailure[
      diagnostics,
      "MathMSSAResidualNormFailed",
      "ROW_OPERATOR",
      <|"Rank" -> rank|>
    ]
  ];
  AssociateTo[diagnostics, "RowOperatorNorm" -> bNorm];
  denominator = aNorm bNorm + cNorm;
  residual = If[
    TrueQ[denominator == 0.],
    0.0,
    Quiet@Check[N[safeNorm2[a . b - c]/denominator, MachinePrecision], $Failed]
  ];
  If[!finiteRealQ[residual],
    Return@makeRankFailure[
      diagnostics,
      "MathMSSASolveResidualFailed",
      "ROW_OPERATOR",
      <|"Rank" -> rank|>
    ]
  ];
  AssociateTo[diagnostics, "RowOperatorSolveRelativeResidual" -> residual];
  eta = Quiet@Check[N[b . knownRowVector, MachinePrecision], $Failed];
  If[eta === $Failed || !finiteRealVectorQ[eta] || Length[eta] =!= channelCount,
    Return@makeRankFailure[
      diagnostics,
      "MathMSSANonFiniteForecast",
      "ROW_FORECAST",
      <|"Rank" -> rank|>
    ]
  ];
  mssaLog[
    "RANK_DONE",
    <|
      "Rank" -> rank,
      "ForecastH" -> First[eta],
      "RowOperatorNorm" -> bNorm,
      "SolveRelativeResidual" -> residual
    |>
  ];
  makeRankResult[First[eta], diagnostics, None]
];

MathMSSA[y_, x_, window_, ranks_] := Module[
  {
    target, inputs, n, d, xDimensions, trajectoryColumns, maximumRank,
    trajectory, svd, decomposition, byRank
  },
  mssaLog[
    "START",
    <|
      "YDimensions" -> Dimensions[y],
      "XDimensions" -> Dimensions[x],
      "L" -> window,
      "Ranks" -> ranks
    |>
  ];
  If[!finiteRealVectorQ[y],
    Return@makeFailure["MathMSSAInvalidTargetHistory", "INPUT", <|"Argument" -> "y"|>]
  ];
  xDimensions = Dimensions[x];
  If[!finiteRealMatrixQ[x] || Length[xDimensions] =!= 2,
    Return@makeFailure["MathMSSAInvalidInputHistory", "INPUT", <|"Argument" -> "X"|>]
  ];
  If[!IntegerQ[window],
    Return@makeFailure["MathMSSAInvalidWindow", "INPUT", <|"L" -> window|>]
  ];
  If[!rankGridQ[ranks],
    Return@makeFailure["MathMSSAInvalidRanks", "INPUT", <|"Ranks" -> ranks|>]
  ];
  n = Length[y];
  d = Last[xDimensions];
  If[n < 3,
    Return@makeFailure["MathMSSAInsufficientHistory", "INPUT", <|"N" -> n|>]
  ];
  If[d < 1 || xDimensions =!= {n + 1, d},
    Return@makeFailure[
      "MathMSSAInputDimensionMismatch",
      "INPUT",
      <|"ExpectedRows" -> n + 1, "ActualDimensions" -> xDimensions|>
    ]
  ];
  If[window < 2 || window > n - 1,
    Return@makeFailure[
      "MathMSSAInvalidWindow",
      "INPUT",
      <|"L" -> window, "ExpectedRange" -> {2, n - 1}|>
    ]
  ];
  trajectoryColumns = (n - window + 1) + d (n - window + 2);
  maximumRank = Last[ranks];
  If[maximumRank > Min[window, trajectoryColumns],
    Return@makeFailure[
      "MathMSSAInvalidRanks",
      "INPUT",
      <|"MaxRank" -> maximumRank, "AllowedMaximum" -> Min[window, trajectoryColumns]|>
    ]
  ];
  target = toPackedRealVector[y];
  inputs = toPackedRealMatrix[x];
  mssaLog[
    "INPUT_ACCEPTED",
    <|"N" -> n, "InputChannels" -> d, "L" -> window, "Ranks" -> ranks|>
  ];
  trajectory = buildMSSATrajectory[target, inputs, window];
  If[
    !finiteRealMatrixQ[trajectory["TrajectoryMatrix"]] ||
      !finiteRealVectorQ[trajectory["KnownRowVector"]],
    Return@makeFailure["MathMSSATrajectoryFailed", "TRAJECTORY"]
  ];
  mssaLog[
    "TRAJECTORY_READY",
    <|
      "ChannelCount" -> trajectory["ChannelCount"],
      "BlockWidths" -> trajectory["BlockWidths"],
      "TrajectoryDimensions" -> Dimensions[trajectory["TrajectoryMatrix"]],
      "KnownRowVectorLength" -> Length[trajectory["KnownRowVector"]]
    |>
  ];
  svd = decomposeTrajectory[trajectory["TrajectoryMatrix"], ranks];
  If[FailureQ[svd], Return[svd]];
  decomposition = <|
    "MSSAChannelCount" -> trajectory["ChannelCount"],
    "MSSATrajectoryColumns" -> trajectory["TrajectoryColumns"],
    "MSSABlockWidths" -> trajectory["BlockWidths"],
    "MSSASVDMaxRank" -> svd["MaxRank"],
    "MSSARetainedNumericalRank" -> svd["NumericalRank"],
    "MSSASingularValues" -> svd["SingularValues"]
  |>;
  byRank = AssociationThread[
    ranks,
    evaluateRank[
      svd,
      trajectory["EndpointIndices"],
      trajectory["KnownRowVector"],
      #
    ] & /@ ranks
  ];
  mssaLog[
    "DONE",
    <|
      "Ranks" -> ranks,
      "SuccessfulRanks" ->
        Select[
          ranks,
          Function[trajectoryRank, byRank[trajectoryRank]["TechnicalFailure"] === None]
        ]
    |>
  ];
  <|"Decomposition" -> decomposition, "ByRank" -> byRank|>
];

MathMSSA[___] :=
  makeFailure["MathMSSAInvalidArity", "INPUT"];

End[]

EndPackage[]
