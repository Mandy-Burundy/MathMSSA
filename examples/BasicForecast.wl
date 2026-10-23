(* wolframscript -file examples/BasicForecast.wl *)

$HistoryLength = 0;

Get[FileNameJoin[{ParentDirectory[DirectoryName[$InputFileName]], "Kernel", "MathMSSA.wl"}]];
MathMSSA`$MathMSSALogging = False;

y = Sin[2. Pi Range[64] / 17.];
x = Transpose[{0.7 Sin[2. Pi Range[65] / 17. + 0.4]}];
result = MathMSSA`MathMSSA[y, x, 20, {2}];

Print[result["ByRank"][2]["ForecastH"]];
