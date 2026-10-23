(* Synthetic check of MathMSSA against Rssa 1.1.
   Market series are not read. MathMSSA.wl is not modified.
   Reference ForecastH comes only from MathMSSA`MathMSSA.
   Rssa: kind="mssa", svd.method="svd", direction="row", base="original". *)

$HistoryLength = 0;

$stageASeed = 20261023;
$stageANoiseScale = 0.05;
$stageAScale = 2.5;
$stageAPeriod = 17.;
$stageASigmaAtol = 1.*^-10;
$stageASigmaRtol = 1.*^-10;
$stageAValueAtol = 1.*^-10;
$stageAValueRtol = 1.*^-8;
$stageASubspaceTol = 1.*^-7;
$stageAResidualTol = 1.*^-10;
$stageACondLimit = 1.*^8;
$stageAGapLimit = 1.*^-8;

jsonAssoc[rules : {(_String -> _) ...}] := Association[Map[#[[1]] -> jsonAssoc[#[[2]]] &, rules]];
jsonAssoc[items_List] := Map[jsonAssoc, items];
jsonAssoc[other_] := other;

hexFile[path_String] := IntegerString[FileHash[path, "SHA256"], 16, 64];

findRscript[] := Module[{env, separator, names, dirs, fromPath, extras},
  env = Environment["RSCRIPT"];
  separator = If[$OperatingSystem === "Windows", ";", ":"];
  names = If[$OperatingSystem === "Windows", {"Rscript.exe", "Rscript"}, {"Rscript"}];
  dirs = StringSplit[If[StringQ[Environment["PATH"]], Environment["PATH"], ""], separator];
  fromPath = Flatten[Outer[FileNameJoin[{#1, #2}] &, dirs, names]];
  extras = If[$OperatingSystem === "Windows", {}, {"/usr/local/bin/Rscript", "/opt/homebrew/bin/Rscript", "/usr/bin/Rscript"}];
  SelectFirst[Join[If[StringQ[env] && env =!= "", {env}, {}], fromPath, extras], FileExistsQ, $Failed]
];

f64Hash[values_List] := Module[{path, hash},
  path = CreateFile[];
  BinaryWrite[path, N[Flatten[values]], "Real64"];
  Close[path];
  hash = hexFile[path];
  DeleteFile[path];
  hash
];

readF64[path_String, n_Integer] := If[!FileExistsQ[path] || FileByteCount[path] =!= 8 n, $Failed, BinaryReadList[path, "Real64"]];

readMatrix[path_String, dims : {__Integer}] := Module[{raw},
  raw = readF64[path, Times @@ dims];
  If[raw === $Failed, $Failed, If[Length[dims] == 1, raw, ArrayReshape[raw, dims]]]
];
readMatrix[_, _] := $Failed;

asReal[x_] := Which[x === "Infinity", Infinity, x === "-Infinity", -Infinity, NumericQ[x], N[x], True, $Failed];

asInts[x_List] := Round[x];
asInts[x_?NumericQ] := {Round[x]};

cell[x_] := Which[
  x === Null || x === "" || MissingQ[x] || x === $Failed, "",
  x === Infinity, "Infinity",
  x === -Infinity, "-Infinity",
  StringQ[x], x,
  IntegerQ[x], ToString[x],
  NumericQ[x], StringReplace[ToString[FortranForm[N[x]]], " " -> ""],
  True, ""
];

csvColumns = {"Stage", "CaseID", "ModelID", "Date", "L", "Rank", "CheckName", "Status", "WolframValue", "RssaValue", "ExpectedValue", "MaxAbsDiff", "Scale", "AbsTolerance", "RelTolerance", "ConditionWolfram", "ConditionRssa", "BoundaryGapWolfram", "BoundaryGapRssa", "Reason"};

makeRow[fields_Association] := Join[AssociationThread[csvColumns, ConstantArray["", Length[csvColumns]]], fields];

seriesLinear[a_, b_, len_Integer] := a + b * N[Range[len]];
seriesHarmonic[amp_, phase_, len_Integer] := amp * Sin[2. Pi * N[Range[len]] / $stageAPeriod + phase];

constChannels[n_Integer] := {ConstantArray[2.5, n], ConstantArray[-1.25, n + 1], ConstantArray[4., n + 1]};
d1Channels[n_Integer] := {ConstantArray[2.5, n], ConstantArray[-1.25, n + 1]};
linearChannels[n_Integer] := {seriesLinear[0.3, 0.05, n], seriesLinear[-1., 0.02, n + 1], seriesLinear[1.5, -0.01, n + 1]};
harmonicChannels[n_Integer] := {seriesHarmonic[1., 0., n], seriesHarmonic[0.7, 0.4, n + 1], seriesHarmonic[1.3, -0.8, n + 1]};
trendChannels[n_Integer] := {
  seriesLinear[0.2, 0.01, n] + seriesHarmonic[0.5, 0., n],
  seriesLinear[-0.4, 0.03, n + 1] + seriesHarmonic[0.8, 0.5, n + 1],
  seriesLinear[0.6, -0.02, n + 1] + seriesHarmonic[0.4, -0.3, n + 1]
};
noiseChannels[n_Integer] := Module[{base, lens, noise},
  base = trendChannels[n];
  lens = Length /@ base;
  SeedRandom[$stageASeed, Method -> "ExtendedCA"];
  noise = RandomVariate[NormalDistribution[0., $stageANoiseScale], Total[lens]];
  MapThread[#1 + #2 &, {base, TakeList[noise, lens]}]
];
nearChannels[n_Integer] := Module[{eps, ext, wave},
  eps = {1.*^-8, 1.*^-6, 1.*^-4, 1.*^-3, 1.*^-2, 0.05};
  ext = seriesLinear[0.5, 0.01, n + 1];
  wave = seriesHarmonic[1., 0., n + 1];
  Join[{seriesLinear[0.5, 0.01, n]}, Map[ext + # * wave &, eps]]
];

channelNames[nChannels_Integer] := Join[{"y"}, Table["x" <> ToString[i], {i, nChannels - 1}]];
expectedTarget["CONST", n_Integer] := 2.5;
expectedTarget["D1", n_Integer] := 2.5;
expectedTarget["LINEAR", n_Integer] := 0.3 + 0.05 * (n + 1);
expectedTarget["HARMONIC", n_Integer] := Sin[2. Pi * (n + 1) / $stageAPeriod];
expectedTarget["TREND_HARMONIC", n_Integer] := 0.2 + 0.01 * (n + 1) + 0.5 * Sin[2. Pi * (n + 1) / $stageAPeriod];
expectedTarget[_, _] := Null;
ranksOf["CONST"] := {1};
ranksOf["D1"] := {1};
ranksOf["LINEAR"] := {2};
ranksOf["HARMONIC"] := {2};
ranksOf["TREND_HARMONIC"] := {4};
ranksOf["NOISY"] := Range[5];
ranksOf["NEAR_COLLINEAR"] := Range[5];
channelsOf["CONST"] := constChannels;
channelsOf["D1"] := d1Channels;
channelsOf["LINEAR"] := linearChannels;
channelsOf["HARMONIC"] := harmonicChannels;
channelsOf["TREND_HARMONIC"] := trendChannels;
channelsOf["NOISY"] := noiseChannels;
channelsOf["NEAR_COLLINEAR"] := nearChannels;
roleOf["D1"] := "match";
roleOf["NOISY"] := "match";
roleOf["NEAR_COLLINEAR"] := "stability";
roleOf[_] := "analytic";

pointRecord[case_String, n_Integer, l_Integer, channels_List, ranks_List, role_String, property_, baseId_] := <|
  "id" -> StringJoin[case, "_N", ToString[n], "_L", ToString[l], Switch[property, "scale", "_SCALE", "permute", "_PERM", _, ""]],
  "caseId" -> case, "role" -> role, "n" -> n, "l" -> l, "ranks" -> ranks,
  "names" -> channelNames[Length[channels]], "channels" -> channels,
  "expected" -> If[property === Null, expectedTarget[case, n], Null],
  "baseId" -> baseId, "property" -> property
|>;

ordinaryPoints[] := Flatten[Table[Module[{channels, base},
  channels = channelsOf[case][n];
  base = pointRecord[case, n, l, channels, ranksOf[case], roleOf[case], Null, Null];
  If[roleOf[case] =!= "analytic", {base}, {
    base,
    pointRecord[case, n, l, $stageAScale * channels, ranksOf[case], "property", "scale", base["id"]],
    pointRecord[case, n, l, ReplacePart[channels, {2 -> channels[[3]], 3 -> channels[[2]]}], ranksOf[case], "property", "permute", base["id"]]
  }]
], {case, {"CONST", "D1", "LINEAR", "HARMONIC", "TREND_HARMONIC", "NOISY", "NEAR_COLLINEAR"}}, {n, {128, 256}}, {l, {20, 40, 80}}], 3];

boundaryPoints[] := {
  pointRecord["ZEROS", 128, 20, {ConstantArray[0., 128], ConstantArray[0., 129], ConstantArray[0., 129]}, {1}, "boundary", Null, Null],
  pointRecord["CONST_RANK2", 128, 20, constChannels[128], {2}, "boundary", Null, Null],
  pointRecord["HARMONIC_SPLIT", 128, 40, harmonicChannels[128], {1}, "boundary", Null, Null],
  pointRecord["CONST_L2", 128, 2, constChannels[128], {1}, "boundary", Null, Null],
  pointRecord["CONST_LMAX", 128, 127, constChannels[128], {1}, "boundary", Null, Null]
};

designMatrix[channels_List] := Transpose[Rest[channels]];

failureTag[f_Failure] := If[StringQ[f[[1]]], f[[1]], ToString[f[[1]]]];

invalidRows[] := Module[{yShort, calls, results},
  calls = {
    <|"CaseID" -> "INVALID_SHORT_Y", "Expect" -> "MathMSSAInsufficientHistory", "Call" -> Hold[MathMSSA`MathMSSA[{1., 2.}, {{0.}, {0.}, {0.}}, 2, {1}]]|>,
    <|"CaseID" -> "INVALID_X_ROWS", "Expect" -> "MathMSSAInputDimensionMismatch", "Call" -> Hold[MathMSSA`MathMSSA[ConstantArray[1., 16], ConstantArray[0., {16, 1}], 8, {1}]]|>,
    <|"CaseID" -> "INVALID_L1", "Expect" -> "MathMSSAInvalidWindow", "Call" -> Hold[MathMSSA`MathMSSA[ConstantArray[1., 32], ConstantArray[0., {33, 1}], 1, {1}]]|>,
    <|"CaseID" -> "INVALID_L_EQ_N", "Expect" -> "MathMSSAInvalidWindow", "Call" -> Hold[MathMSSA`MathMSSA[ConstantArray[1., 32], ConstantArray[0., {33, 1}], 32, {1}]]|>,
    <|"CaseID" -> "INVALID_NONNUMERIC", "Expect" -> "MathMSSAInvalidTargetHistory", "Call" -> Hold[MathMSSA`MathMSSA[{1., 2., Infinity}, {{0.}, {0.}, {0.}, {0.}}, 2, {1}]]|>,
    <|"CaseID" -> "INVALID_RANK", "Expect" -> "MathMSSAInvalidRanks", "Call" -> Hold[MathMSSA`MathMSSA[ConstantArray[1., 32], ConstantArray[0., {33, 1}], 20, {21}]]|>
  };
  Map[Function[spec, Module[{got, tag, ok},
    got = Quiet @ ReleaseHold[spec["Call"]];
    tag = If[FailureQ[got], failureTag[got], "NoFailure"];
    ok = tag === spec["Expect"];
    makeRow[<|
      "Stage" -> "A", "CaseID" -> spec["CaseID"], "L" -> "", "Rank" -> "", "CheckName" -> "InvalidInput",
      "Status" -> If[ok, "PASS", "FAIL"], "WolframValue" -> tag, "ExpectedValue" -> spec["Expect"],
      "Reason" -> If[ok, "", "INPUT_MISMATCH"]
    |>]
  ]], calls]
];

numericalRank[sigma_List] := Module[{cut, finite},
  finite = Select[sigma, NumericQ[#] && # >= 0. &];
  If[finite === {} || First[finite] == 0., Return[0]];
  cut = Max[1.*^-12, 1.*^-10 * First[finite]];
  Count[finite, v_ /; v > cut]
];

gapValue[sigma_List, rank_Integer] := If[
  Length[sigma] < rank + 1 || !NumericQ[sigma[[1]]] || sigma[[1]] == 0. || !NumericQ[sigma[[rank]]] || !NumericQ[sigma[[rank + 1]]],
  $Failed,
  (sigma[[rank]] - sigma[[rank + 1]]) / sigma[[1]]
];

sideSensitive[cond_, gap_, rank_Integer, nrank_Integer] := !NumericQ[cond] || cond > $stageACondLimit || !NumericQ[gap] || gap <= $stageAGapLimit || rank > nrank;

gapUnidentifiable[gap_, rank_Integer, nrank_Integer] := !NumericQ[gap] || gap <= $stageAGapLimit || rank > nrank;

infNorm[a_] := If[a === $Failed, $Failed, Max[Abs[N[Flatten[{a}]]]]];

pairDiff[a_, b_] := If[a === $Failed || b === $Failed || Dimensions[{a}] =!= Dimensions[{b}] && Dimensions[Developer`ToPackedArray[N[{a}]]] =!= Dimensions[Developer`ToPackedArray[N[{b}]]] && Length[Flatten[{a}]] =!= Length[Flatten[{b}]], $Failed, Max[Abs[Flatten[N[{a}]] - Flatten[N[{b}]]]]];

matrixDiff[a_, b_] := If[a === $Failed || b === $Failed || Dimensions[a] =!= Dimensions[b], $Failed, Max[Abs[Flatten[a - b]]]];

projResidual[aCols_, bCols_] := Module[{qa, qb, coeff, projected},
  If[aCols === $Failed || bCols === $Failed || Dimensions[aCols] =!= Dimensions[bCols], Return[$Failed]];
  qa = Quiet @ Check[Orthogonalize[Transpose[N[aCols]]], $Failed];
  qb = Quiet @ Check[Orthogonalize[Transpose[N[bCols]]], $Failed];
  If[qa === $Failed || qb === $Failed || qa === {} || qb === {} || !MatrixQ[qa] || !MatrixQ[qb], Return[$Failed]];
  coeff = qa . Transpose[qb];
  projected = Transpose[qa] . coeff;
  Quiet @ Check[Norm[Transpose[qb] - projected, 2], $Failed]
];

subspaceDistance[aCols_, bCols_] := Module[{left, right},
  left = projResidual[aCols, bCols];
  right = projResidual[bCols, aCols];
  If[!NumericQ[left] || !NumericQ[right], $Failed, Max[left, right]]
];

diagonalAverage[block_?MatrixQ] := Module[{l, w, n},
  {l, w} = Dimensions[block];
  n = l + w - 1;
  Table[Mean[block[[#, t - # + 1]] & /@ Range[Max[1, t - w + 1], Min[l, t]]], {t, n}]
];

lowRank[u_, sigma_, v_, rank_Integer] := u[[All, 1 ;; rank]] . DiagonalMatrix[sigma[[1 ;; rank]]] . Transpose[v[[All, 1 ;; rank]]];

harnessForecast[v_, endpoints_List, z_List, rank_Integer] := Module[{complement, s, vbar, a, c, solved, eta},
  complement = Complement[Range[Length[v]], endpoints];
  If[Length[z] =!= Length[complement] || rank > Dimensions[v][[2]], Return[<|"A" -> $Failed, "C" -> $Failed, "Forecast" -> $Failed|>]];
  s = v[[endpoints, 1 ;; rank]];
  vbar = v[[complement, 1 ;; rank]];
  a = IdentityMatrix[Length[endpoints]] - s . Transpose[s];
  c = s . Transpose[vbar];
  solved = Quiet @ Check[LinearSolve[a, c], $Failed];
  eta = If[solved === $Failed, $Failed, solved . z];
  <|"A" -> a, "C" -> c, "Forecast" -> If[ListQ[eta] && eta =!= {}, First[eta], $Failed]|>
];

sensitiveCeiling[cond_, scale_] := Module[{c, s},
  c = Which[NumericQ[cond] && 0 < cond < Infinity, cond, cond === Infinity, $stageACondLimit, True, 1.];
  s = If[NumericQ[scale] && scale > 1., scale, 1.];
  1.*^4 * c * $MachineEpsilon * s
];

statusOf[diff_, scale_, atol_, rtol_, sensitive_, cond_: 1.] := Which[
  TrueQ[sensitive] && NumericQ[diff] && diff > sensitiveCeiling[cond, scale], "FAIL",
  TrueQ[sensitive], "NUMERICALLY_SENSITIVE",
  NumericQ[diff] && diff <= atol + rtol * Max[If[NumericQ[scale], scale, 0.], 0.], "PASS",
  True, "FAIL"
];

reasonOf[status_String, kind_String, diff_] := Which[
  status === "PASS" || status === "NOT_RUN" || status === "NOT_APPLICABLE", "",
  status === "UNIDENTIFIABLE", "GAP_UNIDENTIFIABLE",
  status === "NUMERICALLY_SENSITIVE" && NumericQ[diff] && diff > $stageAValueAtol, kind <> "_DIFF",
  status === "NUMERICALLY_SENSITIVE", "NUMERICALLY_SENSITIVE",
  True, kind
];

runWolfram[point_Association] := Module[{y, x, l, ranks, result, traj, svd, sigmaDiag},
  y = N[point["channels"][[1]]];
  x = N[designMatrix[point["channels"]]];
  l = Round[point["l"]];
  ranks = asInts[point["ranks"]];
  result = Quiet @ MathMSSA`MathMSSA[y, x, l, ranks];
  If[FailureQ[result], Return[<|"Failure" -> result|>]];
  traj = MathMSSA`Private`buildMSSATrajectory[y, x, l];
  svd = MathMSSA`Private`decomposeTrajectory[traj["TrajectoryMatrix"], ranks];
  If[FailureQ[svd], Return[<|"Failure" -> svd|>]];
  sigmaDiag = Quiet @ Check[SingularValueList[traj["TrajectoryMatrix"], Min[Last[ranks] + 1, Min[Dimensions[traj["TrajectoryMatrix"]]]]], $Failed];
  <|"Result" -> result, "Trajectory" -> traj, "SVD" -> svd, "SigmaDiagnostic" -> If[ListQ[sigmaDiag], sigmaDiag, $Failed]|>
];

rankPayload[remote_Association, rank_Integer] := Lookup[Lookup[remote, "ranks", <||>], ToString[rank], <||>];

comparePoint[point_Association, remote_Association, local_Association] := Module[
  {id, case, role, l, ranks, rows, stop, inputHash, remoteHash, widths, endpoints, tW, zW, tR, zR, uW, vW, uR, vR, sigmaW, sigmaR, sigmaDiag, frob},
  id = point["id"]; case = point["caseId"]; role = point["role"]; l = Round[point["l"]]; ranks = asInts[point["ranks"]];
  rows = {};
  stop = False;
  inputHash = f64Hash[point["channels"]];
  If[!TrueQ[Lookup[remote, "ok", False]],
    Return[Append[rows, makeRow[<|"Stage" -> "A", "CaseID" -> case, "L" -> l, "Rank" -> "", "CheckName" -> "RssaCall", "Status" -> "FAIL", "Reason" -> "FORECAST_MISMATCH", "RssaValue" -> ToString[Lookup[remote, "message", ""]]|>]]]
  ];
  remoteHash = "";
  widths = asInts[Lookup[remote, "blockWidths", {}]];
  endpoints = asInts[Lookup[remote, "endpoints", {}]];
  tW = local["Trajectory"]["TrajectoryMatrix"];
  zW = local["Trajectory"]["KnownRowVector"];
  tR = readMatrix[FileNameJoin[{remote["RawDir"], id <> ".traj.f64"}], asInts[remote["dims"]["traj"]]];
  zR = readF64[FileNameJoin[{remote["RawDir"], id <> ".z.f64"}], Round[remote["dims"]["z"]]];
  uW = local["SVD"]["U"]; vW = local["SVD"]["V"]; sigmaW = local["SVD"]["SingularValues"];
  uR = readMatrix[FileNameJoin[{remote["RawDir"], id <> ".U.f64"}], asInts[remote["dims"]["U"]]];
  vR = readMatrix[FileNameJoin[{remote["RawDir"], id <> ".V.f64"}], asInts[remote["dims"]["V"]]];
  sigmaR = N[Lookup[remote, "sigma", {}]];
  sigmaDiag = local["SigmaDiagnostic"];
  frob = Sqrt[Total[Flatten[tW]^2]];
  rows = {
    makeRow[<|"Stage" -> "A", "CaseID" -> case, "L" -> l, "CheckName" -> "BlockWidths", "Status" -> If[local["Trajectory"]["BlockWidths"] === widths, "PASS", "FAIL"], "MaxAbsDiff" -> If[local["Trajectory"]["BlockWidths"] === widths, 0, 1], "Reason" -> If[local["Trajectory"]["BlockWidths"] === widths, "", "EMBEDDING_MISMATCH"], "WolframValue" -> ToString[local["Trajectory"]["BlockWidths"]], "RssaValue" -> ToString[widths]|>],
    makeRow[<|"Stage" -> "A", "CaseID" -> case, "L" -> l, "CheckName" -> "Trajectory", "Status" -> statusOf[matrixDiff[tW, tR], infNorm[tW], 0., 0., False], "MaxAbsDiff" -> matrixDiff[tW, tR], "Scale" -> infNorm[tW], "Reason" -> reasonOf[statusOf[matrixDiff[tW, tR], infNorm[tW], 0., 0., False], "EMBEDDING_MISMATCH", matrixDiff[tW, tR]]|>],
    makeRow[<|"Stage" -> "A", "CaseID" -> case, "L" -> l, "CheckName" -> "KnownVector", "Status" -> statusOf[pairDiff[zW, zR], infNorm[zW], 0., 0., False], "MaxAbsDiff" -> pairDiff[zW, zR], "Scale" -> infNorm[zW], "Reason" -> reasonOf[statusOf[pairDiff[zW, zR], infNorm[zW], 0., 0., False], "EMBEDDING_MISMATCH", pairDiff[zW, zR]]|>]
  };
  rows = Append[rows, makeRow[<|"Stage" -> "A", "CaseID" -> case, "L" -> l, "CheckName" -> "RssaStoredTrajectory", "Status" -> statusOf[asReal[Lookup[remote, "extmatDiff", $Failed]], infNorm[tW], 1.*^-12, 0., False], "MaxAbsDiff" -> asReal[Lookup[remote, "extmatDiff", $Failed]], "Scale" -> infNorm[tW], "AbsTolerance" -> 1.*^-12, "Reason" -> reasonOf[statusOf[asReal[Lookup[remote, "extmatDiff", $Failed]], infNorm[tW], 1.*^-12, 0., False], "EMBEDDING_MISMATCH", asReal[Lookup[remote, "extmatDiff", $Failed]]]|>]];
  If[MemberQ[rows[[All, "Status"]], "FAIL"], stop = True];
  rows = Join[rows, If[sigmaR === {} || frob > 0. && Max[sigmaR] > frob * (1. + 1.*^-6),
    {makeRow[<|"Stage" -> "A", "CaseID" -> case, "L" -> l, "CheckName" -> "SingularValueScale", "Status" -> "FAIL", "Reason" -> "SVD_OR_SUBSPACE_MISMATCH", "Scale" -> frob, "RssaValue" -> If[sigmaR === {}, "", Max[sigmaR]]|>]},
    {makeRow[<|"Stage" -> "A", "CaseID" -> case, "L" -> l, "CheckName" -> "SingularValueScale", "Status" -> "PASS", "Scale" -> frob, "RssaValue" -> Max[sigmaR]|>]}
  ]];
  If[stop, Return[rows]];
  Join[rows, Flatten[Map[Function[rank, compareRank[point, remote, local, rank, tW, zW, uW, vW, uR, vR, sigmaW, sigmaR, sigmaDiag, widths]], ranks], 1]]
];

compareRank[point_, remote_, local_, rank_Integer, tW_, zW_, uW_, vW_, uR_, vR_, sigmaW_, sigmaR_, sigmaDiag_, widths_] := Module[
  {case, l, role, payload, condW, condR, gapW, gapR, nrankW, nrankR, sensitive, unidentifiable, condCeiling, harness, aR, cR, diffSigma, scaleSigma, rows, forecastW, forecastR, public, byRank, lowW, lowR, avgW, avgR, blocks, forecastDiff, forecastScale, expected, analyticDiff, harnessDiff, resW, resR, resDiff, resScale, basisStatus},
  case = point["caseId"]; l = Round[point["l"]]; role = point["role"];
  payload = rankPayload[remote, rank];
  byRank = local["Result"]["ByRank"][rank];
  public = If[AssociationQ[byRank], byRank["ForecastH"], $Failed];
  condW = If[AssociationQ[byRank], asReal[byRank["RowOperatorConditionNumber"]], $Failed];
  condR = asReal[Lookup[payload, "condition", $Failed]];
  gapW = If[ListQ[sigmaDiag], gapValue[sigmaDiag, rank], $Failed];
  gapR = gapValue[sigmaR, rank];
  nrankW = If[ListQ[sigmaDiag], numericalRank[sigmaDiag], 0];
  nrankR = numericalRank[sigmaR];
  sensitive = sideSensitive[condW, gapW, rank, nrankW] || sideSensitive[condR, gapR, rank, nrankR];
  unidentifiable = gapUnidentifiable[gapW, rank, nrankW] || gapUnidentifiable[gapR, rank, nrankR];
  condCeiling = If[MemberQ[{condW, condR}, Infinity], Infinity, Max[1., Sequence @@ Select[{condW, condR}, NumericQ[#] && # > 0 &]]];
  basisStatus[diff_, scale_, atol_, rtol_] := If[TrueQ[unidentifiable], "UNIDENTIFIABLE", statusOf[diff, scale, atol, rtol, sensitive, condCeiling]];
  harness = harnessForecast[vW, local["Trajectory"]["EndpointIndices"], zW, rank];
  aR = readMatrix[FileNameJoin[{remote["RawDir"], StringJoin[point["id"], ".r", ToString[rank], ".A.f64"]}], asInts[Lookup[payload, "aDim", {}]]];
  cR = readMatrix[FileNameJoin[{remote["RawDir"], StringJoin[point["id"], ".r", ToString[rank], ".C.f64"]}], asInts[Lookup[payload, "cDim", {}]]];
  diffSigma = If[Length[sigmaW] >= rank && Length[sigmaR] >= rank, Max[Abs[sigmaW[[1 ;; rank]] - sigmaR[[1 ;; rank]]]], $Failed];
  scaleSigma = If[diffSigma === $Failed, $Failed, Max[infNorm[sigmaW[[1 ;; rank]]], infNorm[sigmaR[[1 ;; rank]]]]];
  forecastW = If[NumericQ[public], public, $Failed];
  forecastR = asReal[Lookup[payload, "forecast", $Failed]];
  forecastDiff = pairDiff[forecastW, forecastR];
  forecastScale = Max[infNorm[forecastW], infNorm[forecastR]];
  harnessDiff = pairDiff[forecastW, harness["Forecast"]];
  lowW = If[uW === $Failed || vW === $Failed || Length[sigmaW] < rank, $Failed, lowRank[uW, sigmaW, vW, rank]];
  lowR = If[uR === $Failed || vR === $Failed || Length[sigmaR] < rank, $Failed, lowRank[uR, sigmaR[[1 ;; Dimensions[uR][[2]]]], vR, rank]];
  blocks = If[lowW === $Failed, {}, TakeList[Transpose[lowW], widths] // Map[Transpose]];
  avgW = If[blocks === {} || role === "stability" || role === "boundary", $Failed, Map[diagonalAverage, blocks]];
  avgR = If[avgW === $Failed, $Failed, Lookup[payload, "recon", $Failed]];
  expected = point["expected"];
  analyticDiff = If[!NumericQ[expected] || !NumericQ[forecastW], $Failed, Abs[forecastW - expected]];
  resW = If[AssociationQ[byRank], asReal[byRank["RowOperatorSolveRelativeResidual"]], $Failed];
  resR = asReal[Lookup[payload, "residual", $Failed]];
  resDiff = If[NumericQ[resW] && NumericQ[resR], Abs[resW - resR], $Failed];
  resScale = Max[infNorm[resW], infNorm[resR]];
  rows = {
    makeRow[<|"Stage" -> "A", "CaseID" -> case, "L" -> l, "Rank" -> rank, "CheckName" -> "SingularValues", "Status" -> statusOf[diffSigma, scaleSigma, $stageASigmaAtol, $stageASigmaRtol, sensitive, condCeiling], "MaxAbsDiff" -> diffSigma, "Scale" -> scaleSigma, "AbsTolerance" -> $stageASigmaAtol, "RelTolerance" -> $stageASigmaRtol, "ConditionWolfram" -> condW, "ConditionRssa" -> condR, "BoundaryGapWolfram" -> gapW, "BoundaryGapRssa" -> gapR, "Reason" -> reasonOf[statusOf[diffSigma, scaleSigma, $stageASigmaAtol, $stageASigmaRtol, sensitive, condCeiling], "SVD_OR_SUBSPACE_MISMATCH", diffSigma]|>],
    makeRow[<|"Stage" -> "A", "CaseID" -> case, "L" -> l, "Rank" -> rank, "CheckName" -> "LeftSubspace", "Status" -> basisStatus[subspaceDistance[uW[[All, 1 ;; rank]], uR[[All, 1 ;; rank]]], 1., $stageASubspaceTol, 0.], "MaxAbsDiff" -> subspaceDistance[uW[[All, 1 ;; rank]], uR[[All, 1 ;; rank]]], "AbsTolerance" -> $stageASubspaceTol, "ConditionWolfram" -> condW, "ConditionRssa" -> condR, "BoundaryGapWolfram" -> gapW, "BoundaryGapRssa" -> gapR, "Reason" -> reasonOf[basisStatus[subspaceDistance[uW[[All, 1 ;; rank]], uR[[All, 1 ;; rank]]], 1., $stageASubspaceTol, 0.], "SVD_OR_SUBSPACE_MISMATCH", subspaceDistance[uW[[All, 1 ;; rank]], uR[[All, 1 ;; rank]]]]|>],
    makeRow[<|"Stage" -> "A", "CaseID" -> case, "L" -> l, "Rank" -> rank, "CheckName" -> "RightSubspace", "Status" -> basisStatus[subspaceDistance[vW[[All, 1 ;; rank]], vR[[All, 1 ;; rank]]], 1., $stageASubspaceTol, 0.], "MaxAbsDiff" -> subspaceDistance[vW[[All, 1 ;; rank]], vR[[All, 1 ;; rank]]], "AbsTolerance" -> $stageASubspaceTol, "ConditionWolfram" -> condW, "ConditionRssa" -> condR, "BoundaryGapWolfram" -> gapW, "BoundaryGapRssa" -> gapR, "Reason" -> reasonOf[basisStatus[subspaceDistance[vW[[All, 1 ;; rank]], vR[[All, 1 ;; rank]]], 1., $stageASubspaceTol, 0.], "SVD_OR_SUBSPACE_MISMATCH", subspaceDistance[vW[[All, 1 ;; rank]], vR[[All, 1 ;; rank]]]]|>],
    makeRow[<|"Stage" -> "A", "CaseID" -> case, "L" -> l, "Rank" -> rank, "CheckName" -> "LowRank", "Status" -> statusOf[matrixDiff[lowW, lowR], infNorm[lowW], $stageAValueAtol, $stageAValueRtol, sensitive, condCeiling], "MaxAbsDiff" -> matrixDiff[lowW, lowR], "Scale" -> infNorm[lowW], "AbsTolerance" -> $stageAValueAtol, "RelTolerance" -> $stageAValueRtol, "ConditionWolfram" -> condW, "ConditionRssa" -> condR, "BoundaryGapWolfram" -> gapW, "BoundaryGapRssa" -> gapR, "Reason" -> reasonOf[statusOf[matrixDiff[lowW, lowR], infNorm[lowW], $stageAValueAtol, $stageAValueRtol, sensitive, condCeiling], "SVD_OR_SUBSPACE_MISMATCH", matrixDiff[lowW, lowR]]|>],
    makeRow[<|"Stage" -> "A", "CaseID" -> case, "L" -> l, "Rank" -> rank, "CheckName" -> "OperatorA", "Status" -> basisStatus[matrixDiff[harness["A"], aR], infNorm[harness["A"]], $stageAValueAtol, $stageAValueRtol], "MaxAbsDiff" -> matrixDiff[harness["A"], aR], "Scale" -> infNorm[harness["A"]], "AbsTolerance" -> $stageAValueAtol, "RelTolerance" -> $stageAValueRtol, "ConditionWolfram" -> condW, "ConditionRssa" -> condR, "BoundaryGapWolfram" -> gapW, "BoundaryGapRssa" -> gapR, "Reason" -> reasonOf[basisStatus[matrixDiff[harness["A"], aR], infNorm[harness["A"]], $stageAValueAtol, $stageAValueRtol], "FORECAST_MISMATCH", matrixDiff[harness["A"], aR]]|>],
    makeRow[<|"Stage" -> "A", "CaseID" -> case, "L" -> l, "Rank" -> rank, "CheckName" -> "OperatorC", "Status" -> basisStatus[matrixDiff[harness["C"], cR], infNorm[harness["C"]], $stageAValueAtol, $stageAValueRtol], "MaxAbsDiff" -> matrixDiff[harness["C"], cR], "Scale" -> infNorm[harness["C"]], "AbsTolerance" -> $stageAValueAtol, "RelTolerance" -> $stageAValueRtol, "ConditionWolfram" -> condW, "ConditionRssa" -> condR, "BoundaryGapWolfram" -> gapW, "BoundaryGapRssa" -> gapR, "Reason" -> reasonOf[basisStatus[matrixDiff[harness["C"], cR], infNorm[harness["C"]], $stageAValueAtol, $stageAValueRtol], "FORECAST_MISMATCH", matrixDiff[harness["C"], cR]]|>],
    makeRow[<|"Stage" -> "A", "CaseID" -> case, "L" -> l, "Rank" -> rank, "CheckName" -> "HarnessForecast", "Status" -> If[NumericQ[harnessDiff] && harnessDiff <= $stageAValueAtol + $stageAValueRtol * Max[forecastScale, 0.], "PASS", "FAIL"], "WolframValue" -> forecastW, "ExpectedValue" -> harness["Forecast"], "MaxAbsDiff" -> harnessDiff, "Scale" -> forecastScale, "AbsTolerance" -> $stageAValueAtol, "RelTolerance" -> $stageAValueRtol, "Reason" -> If[NumericQ[harnessDiff] && harnessDiff <= $stageAValueAtol + $stageAValueRtol * Max[forecastScale, 0.], "", "HARNESS_FORECAST_DIVERGED"]|>],
    makeRow[<|"Stage" -> "A", "CaseID" -> case, "L" -> l, "Rank" -> rank, "CheckName" -> "ForecastH", "Status" -> statusOf[forecastDiff, forecastScale, $stageAValueAtol, $stageAValueRtol, sensitive, condCeiling], "WolframValue" -> forecastW, "RssaValue" -> forecastR, "MaxAbsDiff" -> forecastDiff, "Scale" -> forecastScale, "AbsTolerance" -> $stageAValueAtol, "RelTolerance" -> $stageAValueRtol, "ConditionWolfram" -> condW, "ConditionRssa" -> condR, "BoundaryGapWolfram" -> gapW, "BoundaryGapRssa" -> gapR, "Reason" -> reasonOf[statusOf[forecastDiff, forecastScale, $stageAValueAtol, $stageAValueRtol, sensitive, condCeiling], "FORECAST_MISMATCH", forecastDiff]|>],
    makeRow[<|"Stage" -> "A", "CaseID" -> case, "L" -> l, "Rank" -> rank, "CheckName" -> "SolveResidual", "Status" -> statusOf[resDiff, resScale, $stageAResidualTol, 0., sensitive, condCeiling], "WolframValue" -> resW, "RssaValue" -> resR, "MaxAbsDiff" -> resDiff, "Scale" -> resScale, "AbsTolerance" -> $stageAResidualTol, "ConditionWolfram" -> condW, "ConditionRssa" -> condR, "BoundaryGapWolfram" -> gapW, "BoundaryGapRssa" -> gapR, "Reason" -> reasonOf[statusOf[resDiff, resScale, $stageAResidualTol, 0., sensitive, condCeiling], "FORECAST_MISMATCH", resDiff]|>]
  };
  If[NumericQ[expected],
    AppendTo[rows, makeRow[<|"Stage" -> "A", "CaseID" -> case, "L" -> l, "Rank" -> rank, "CheckName" -> "AnalyticTarget", "Status" -> statusOf[analyticDiff, infNorm[expected], $stageAValueAtol, $stageAValueRtol, sensitive, condCeiling], "WolframValue" -> forecastW, "RssaValue" -> forecastR, "ExpectedValue" -> expected, "MaxAbsDiff" -> analyticDiff, "Scale" -> infNorm[expected], "AbsTolerance" -> $stageAValueAtol, "RelTolerance" -> $stageAValueRtol, "ConditionWolfram" -> condW, "ConditionRssa" -> condR, "BoundaryGapWolfram" -> gapW, "BoundaryGapRssa" -> gapR, "Reason" -> reasonOf[statusOf[analyticDiff, infNorm[expected], $stageAValueAtol, $stageAValueRtol, sensitive, condCeiling], "FORECAST_MISMATCH", analyticDiff]|>]]
  ];
  If[avgW =!= $Failed && ListQ[avgR],
    AppendTo[rows, makeRow[<|"Stage" -> "A", "CaseID" -> case, "L" -> l, "Rank" -> rank, "CheckName" -> "DiagonalAverage", "Status" -> statusOf[matrixDiff[avgW, N[avgR]], infNorm[avgW], $stageAValueAtol, $stageAValueRtol, sensitive, condCeiling], "MaxAbsDiff" -> matrixDiff[avgW, N[avgR]], "Scale" -> infNorm[avgW], "AbsTolerance" -> $stageAValueAtol, "RelTolerance" -> $stageAValueRtol, "Reason" -> reasonOf[statusOf[matrixDiff[avgW, N[avgR]], infNorm[avgW], $stageAValueAtol, $stageAValueRtol, sensitive, condCeiling], "SVD_OR_SUBSPACE_MISMATCH", matrixDiff[avgW, N[avgR]]]|>]]
  ];
  rows
];

invarianceRows[points_List, forecasts_Association] := Flatten[Map[Function[point, Module[{base, factor, gotW, gotR, expW, expR, diffW, diffR, scale},
  If[!MemberQ[{"scale", "permute"}, point["property"]], Return[{}, Module]];
  base = Lookup[forecasts, point["baseId"], <||>];
  factor = If[point["property"] === "scale", $stageAScale, 1.];
  gotW = Lookup[forecasts, point["id"], <||>]["Wolfram"];
  gotR = Lookup[forecasts, point["id"], <||>]["Rssa"];
  expW = factor * Lookup[base, "Wolfram", $Failed];
  expR = factor * Lookup[base, "Rssa", $Failed];
  diffW = pairDiff[gotW, expW]; diffR = pairDiff[gotR, expR];
  scale = Max[infNorm[expW], infNorm[expR]];
  {makeRow[<|"Stage" -> "A", "CaseID" -> point["caseId"], "L" -> Round[point["l"]], "Rank" -> First[asInts[point["ranks"]]], "CheckName" -> If[point["property"] === "scale", "ScaleWolfram", "PermuteWolfram"], "Status" -> statusOf[diffW, scale, $stageAValueAtol, $stageAValueRtol, False], "WolframValue" -> gotW, "ExpectedValue" -> expW, "MaxAbsDiff" -> diffW, "Scale" -> scale, "AbsTolerance" -> $stageAValueAtol, "RelTolerance" -> $stageAValueRtol, "Reason" -> reasonOf[statusOf[diffW, scale, $stageAValueAtol, $stageAValueRtol, False], "FORECAST_MISMATCH", diffW]|>],
   makeRow[<|"Stage" -> "A", "CaseID" -> point["caseId"], "L" -> Round[point["l"]], "Rank" -> First[asInts[point["ranks"]]], "CheckName" -> If[point["property"] === "scale", "ScaleRssa", "PermuteRssa"], "Status" -> statusOf[diffR, scale, $stageAValueAtol, $stageAValueRtol, False], "RssaValue" -> gotR, "ExpectedValue" -> expR, "MaxAbsDiff" -> diffR, "Scale" -> scale, "AbsTolerance" -> $stageAValueAtol, "RelTolerance" -> $stageAValueRtol, "Reason" -> reasonOf[statusOf[diffR, scale, $stageAValueAtol, $stageAValueRtol, False], "FORECAST_MISMATCH", diffR]|>]}
]], points], 1];

stageStatus[rows_List, points_List] := Module[{byId, role, bad, ordinarySensitive, ordinaryFail, otherFail},
  byId = AssociationThread[points[[All, "id"]], points[[All, "role"]]];
  ordinaryFail = MemberQ[rows, row_ /; row["Status"] === "FAIL" && MemberQ[{"analytic", "match", "property"}, Lookup[byId, row["CaseID"] <> "", ""]]];
  (* CaseID is not the point id. Use the row CaseID joined with role lookup by case via a parallel list. *)
  otherFail = MemberQ[Lookup[rows, "Status"], "FAIL"];
  ordinarySensitive = False;
  Which[
    otherFail, "FAILED",
    True, "PASS"
  ]
];

writeReport[dir_String, runID_String, rows_List, status_String] := Module[{counts, firstFail, forecastRows, lines},
  counts = Counts[Lookup[rows, "Status"]];
  firstFail = SelectFirst[rows, #["Status"] === "FAIL" &, <||>];
  forecastRows = Select[rows, #["CheckName"] === "ForecastH" &];
  lines = {
    "# MSSA vs Rssa, stage A",
    "",
    StringJoin["RunID: ", runID],
    StringJoin["Stage status: ", status],
    "",
    StringJoin["PASS ", ToString[Lookup[counts, "PASS", 0]], ", FAIL ", ToString[Lookup[counts, "FAIL", 0]], ", NUMERICALLY_SENSITIVE ", ToString[Lookup[counts, "NUMERICALLY_SENSITIVE", 0]], ", UNIDENTIFIABLE ", ToString[Lookup[counts, "UNIDENTIFIABLE", 0]], ", NOT_RUN ", ToString[Lookup[counts, "NOT_RUN", 0]]],
    "",
    StringJoin["First recorded error: ", If[firstFail === <||>, "none", StringRiffle[Lookup[firstFail, {"CaseID", "L", "Rank", "CheckName", "Reason"}], " "]]],
    "",
    "## Forecasts",
    "",
    "| Case | L | r | Wolfram | Rssa | |diff| | Status |",
    "|---|---:|---:|---:|---:|---:|---|"
  };
  lines = Join[lines, Map[StringJoin["| ", #["CaseID"], " | ", cell[#["L"]], " | ", cell[#["Rank"]], " | ", cell[#["WolframValue"]], " | ", cell[#["RssaValue"]], " | ", cell[#["MaxAbsDiff"]], " | ", #["Status"], " |"] &, forecastRows]];
  Export[FileNameJoin[{dir, "report.md"}], StringRiffle[lines, "\n"], "Text", CharacterEncoding -> "UTF-8"]
];

recordedMax[rows_List, name_String] := Module[{vals, max},
  vals = Select[rows, #["CheckName"] === name && #["Status"] =!= "UNIDENTIFIABLE" && NumericQ[#["MaxAbsDiff"]] &];
  If[vals === {}, "none",
    max = Max[Lookup[vals, "MaxAbsDiff"]];
    If[max == 0., "0", StringReplace[ToString[max, InputForm], "*^" -> "e"]]
  ]
];

writeReference[validationDir_String, runID_String, rows_List, status_String] := Module[{counts, lines},
  counts = Counts[Lookup[rows, "Status"]];
  lines = {
    "# Reference run",
    "",
    StringJoin["RunID: ", runID],
    "Rssa: 1.1",
    StringJoin["Stage status: ", status],
    "",
    "A run is passed when the stage status is PASS. That requires no FAIL row, and no NUMERICALLY_SENSITIVE row on SingularValues, ForecastH, AnalyticTarget, ScaleWolfram, ScaleRssa, PermuteWolfram, or PermuteRssa for the ordinary series CONST, D1, LINEAR, HARMONIC, TREND_HARMONIC, and NOISY. UNIDENTIFIABLE does not fail the run.",
    "",
    StringJoin["PASS ", ToString[Lookup[counts, "PASS", 0]], ", FAIL ", ToString[Lookup[counts, "FAIL", 0]], ", NUMERICALLY_SENSITIVE ", ToString[Lookup[counts, "NUMERICALLY_SENSITIVE", 0]], ", UNIDENTIFIABLE ", ToString[Lookup[counts, "UNIDENTIFIABLE", 0]], ", NOT_RUN ", ToString[Lookup[counts, "NOT_RUN", 0]]],
    "",
    "Maximum absolute difference on rows other than UNIDENTIFIABLE:",
    "",
    StringJoin["- Trajectory, Wolfram against the explicit Hankel matrix: ", recordedMax[rows, "Trajectory"]],
    StringJoin["- KnownVector: ", recordedMax[rows, "KnownVector"]],
    StringJoin["- SingularValues, Wolfram against Rssa: ", recordedMax[rows, "SingularValues"]],
    StringJoin["- LeftSubspace, projection residual: ", recordedMax[rows, "LeftSubspace"]],
    StringJoin["- RightSubspace, projection residual: ", recordedMax[rows, "RightSubspace"]],
    StringJoin["- ForecastH, Wolfram against Rssa: ", recordedMax[rows, "ForecastH"]],
    "",
    "UNIDENTIFIABLE marks LeftSubspace, RightSubspace, OperatorA, and OperatorC when the relative singular gap is at most 1e-8 or the rank exceeds the retained numerical rank. The subspace and those operators are then not identified by the gap criterion. ForecastH at the same rank is still compared."
  };
  Export[FileNameJoin[{validationDir, "reference_run.md"}], StringRiffle[lines, "\n"], "Text", CharacterEncoding -> "UTF-8"]
];

runStageA[validationDir_String] := Module[
  {projectRoot, mathPath, rscript, rfile, stamp, runID, runDir, rawDir, points, fixturePath, outputPath, reloaded, reloadDiff, proc, payload, remoteById, rows, localCache, forecasts, inv, allRows, status, mathHashBefore, seconds},
  projectRoot = DirectoryName[validationDir];
  mathPath = FileNameJoin[{projectRoot, "Kernel", "MathMSSA.wl"}];
  rscript = findRscript[];
  rfile = FileNameJoin[{validationDir, "rssa_forecast.R"}];
  If[!FileExistsQ[mathPath] || !FileExistsQ[rscript] || !FileExistsQ[rfile], Print["STAGE_A_INCOMPLETE missing executable or kernel"]; Return[$Failed]];
  stamp = DateList[TimeZone -> "UTC"];
  runID = StringJoin[IntegerString[Round[stamp[[#]]] & /@ Range[5], 10, {4, 2, 2, 2, 2}], IntegerString[Round[stamp[[6]]], 10, 2], "Z"];
  runID = StringInsert[runID, "T", 9];
  runDir = FileNameJoin[{validationDir, "runs", runID}];
  rawDir = FileNameJoin[{runDir, "raw"}];
  CreateDirectory[rawDir, CreateIntermediateDirectories -> True];
  mathHashBefore = hexFile[mathPath];
  Get[mathPath];
  MathMSSA`$MathMSSALogging = False;
  points = Join[ordinaryPoints[], boundaryPoints[]];
  fixturePath = FileNameJoin[{runDir, "fixtures.json"}];
  Export[fixturePath, <|"contract" -> "1.0", "noiseSeed" -> $stageASeed, "noiseMethod" -> "ExtendedCA", "noiseScale" -> $stageANoiseScale, "points" -> points|>, "JSON"];
  reloaded = jsonAssoc[Import[fixturePath, "JSON"]]["points"];
  reloadDiff = Max[Map[Max[Abs[N[Flatten[#[[1]]["channels"]]] - N[Flatten[#[[2]]["channels"]]]]] &, Transpose[{points, reloaded}]]];
  If[reloadDiff =!= 0., Print["INPUT_MISMATCH json roundtrip ", reloadDiff]; Return[$Failed, Module]];
  points = reloaded;
  outputPath = FileNameJoin[{runDir, "rssa_output.json"}];
  seconds = AbsoluteTiming[
    proc = RunProcess[{rscript, rfile, fixturePath, outputPath, rawDir}];
    If[proc["ExitCode"] =!= 0,
      Export[FileNameJoin[{runDir, "r.log"}], StringJoin[proc["StandardOutput"], "\n", proc["StandardError"]], "Text"];
      Export[FileNameJoin[{runDir, "report.md"}], StringJoin["# MSSA vs Rssa, stage A\n\nStage status: INCOMPLETE\n\nRscript exit ", ToString[proc["ExitCode"]], ". Log: r.log\n"], "Text"];
      Print["Rssa exit ", proc["ExitCode"]];
      Return[$Failed, Module]
    ];
    payload = jsonAssoc[Import[outputPath, "JSON"]];
    If[payload["rssa"] =!= "1.1", Print["Rssa version ", payload["rssa"]]; Return[$Failed, Module]];
    remoteById = AssociationThread[payload["points"][[All, "id"]], Map[Join[#, <|"RawDir" -> rawDir|>] &, payload["points"]]];
    localCache = Association[Map[#["id"] -> runWolfram[#] &, points]];
    rows = Flatten[Map[Function[point, Module[{local},
      local = localCache[point["id"]];
      If[KeyExistsQ[local, "Failure"],
        {makeRow[<|"Stage" -> "A", "CaseID" -> point["caseId"], "L" -> Round[point["l"]], "CheckName" -> "WolframCall", "Status" -> "FAIL", "Reason" -> "FORECAST_MISMATCH", "WolframValue" -> failureTag[local["Failure"]]|>]},
        comparePoint[point, remoteById[point["id"]], local]
      ]
    ]], points], 1];
    forecasts = Association[Map[Function[point, point["id"] -> <|"Wolfram" -> (localCache[point["id"]]["Result"]["ByRank"][First[asInts[point["ranks"]]]]["ForecastH"]), "Rssa" -> asReal[rankPayload[remoteById[point["id"]], First[asInts[point["ranks"]]]]["forecast"]]|>], Select[points, !KeyExistsQ[localCache[#["id"]], "Failure"] &]]];
    inv = invarianceRows[points, forecasts];
    allRows = Join[invalidRows[], rows, inv];
    status = If[MemberQ[Lookup[allRows, "Status"], "FAIL"], "FAILED", If[MemberQ[Select[allRows, MemberQ[{"SingularValues", "ForecastH", "AnalyticTarget", "ScaleWolfram", "ScaleRssa", "PermuteWolfram", "PermuteRssa"}, #["CheckName"]] && MemberQ[{"CONST", "D1", "LINEAR", "HARMONIC", "TREND_HARMONIC", "NOISY"}, #["CaseID"]] &], row_ /; row["Status"] === "NUMERICALLY_SENSITIVE"], "INCOMPLETE", "PASS"]];
    Export[FileNameJoin[{runDir, "comparisons.csv"}], Prepend[Map[cell /@ Lookup[#, csvColumns] &, allRows], csvColumns], "CSV"];
    writeReport[runDir, runID, allRows, status];
    writeReference[validationDir, runID, allRows, status];
    Export[FileNameJoin[{runDir, "manifest.json"}], <|
      "contract" -> "1.0", "runID" -> runID, "stage" -> "A", "status" -> status,
      "wolfram" -> $Version, "os" -> $System, "rssa" -> "1.1", "r" -> payload["rVersion"],
      "mathMSSA" -> "Kernel/MathMSSA.wl", "mathMSSAHashBefore" -> mathHashBefore, "mathMSSAHashAfter" -> hexFile[mathPath],
      "runnerHash" -> hexFile[FileNameJoin[{validationDir, "MSSA_RSSA_STAGE_A.wl"}]],
      "rScriptHash" -> hexFile[rfile], "fixtureHash" -> hexFile[fixturePath],
      "tolerances" -> <|"sigmaAtol" -> $stageASigmaAtol, "sigmaRtol" -> $stageASigmaRtol, "valueAtol" -> $stageAValueAtol, "valueRtol" -> $stageAValueRtol, "subspace" -> $stageASubspaceTol, "residual" -> $stageAResidualTol, "cond" -> $stageACondLimit, "gap" -> $stageAGapLimit|>,
      "noise" -> <|"seed" -> $stageASeed, "method" -> "ExtendedCA", "scale" -> $stageANoiseScale|>,
      "counts" -> Counts[Lookup[allRows, "Status"]]
    |>, "JSON"];
    status
  ];
  Print["RunID ", runID];
  Print["Status ", seconds[[2]]];
  Print["Seconds ", seconds[[1]]];
  Print["Report ", FileNameJoin[{runDir, "report.md"}]];
  seconds[[2]]
];

If[StringQ[$InputFileName], runStageA[DirectoryName[$InputFileName]]];
