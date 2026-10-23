(* wolframscript -file examples/ParallelPoints.wl
   One subkernel per window. Ranks of one window stay together. *)

$HistoryLength = 0;

root = ParentDirectory[DirectoryName[$InputFileName]];
Get[FileNameJoin[{root, "Kernel", "MathMSSAParallel.wl"}]];

jobs = Table[
  With[{phase = 0.3 k},
    <|
      "y" -> Sin[2. Pi Range[64] / 17. + phase],
      "X" -> Transpose[{0.7 Sin[2. Pi Range[65] / 17. + phase + 0.4]}],
      "L" -> 20,
      "Ranks" -> {2}
    |>
  ],
  {k, 0, 7}
];

batch = MathMSSAParallel`MathMSSAParallel[jobs];
Print["Mode ", batch["ExecutionMode"]];
Print["Kernels ", $KernelCount];
Print[Map[#["ByRank"][2]["ForecastH"] &, batch["Results"]]];
