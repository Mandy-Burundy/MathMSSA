(* Independent MathMSSA windows on subkernels.
   One job is one window: its ranks share one SVD inside MathMSSA.
   The subkernel does not see dates, prices, or stored results. *)

If[StringQ[$InputFileName], Get[FileNameJoin[{DirectoryName[$InputFileName], "MathMSSA.wl"}]]];

BeginPackage["MathMSSAParallel`", {"MathMSSA`"}]

MathMSSAParallel::usage =
  "MathMSSAParallel[{job1, job2, ...}] evaluates MathMSSA on independent jobs. Each job is an association with keys y, X, L, and Ranks. One subkernel receives one job. The result is an association with ExecutionMode PARALLEL or SERIAL and Results, the list of MathMSSA return values in job order.";

Begin["`Private`"]

$mathPath = FileNameJoin[{DirectoryName[$InputFileName], "MathMSSA.wl"}];

ensurePool[] := Module[{target, before},
  target = $ProcessorCount;
  before = $KernelCount;
  If[before < target, Block[{$ProgressReporting = False}, Quiet @ Check[LaunchKernels[target - before], Null]]];
  $KernelCount
];

loadWorkers[] := Module[{loaded},
  If[!FileExistsQ[$mathPath], Return[$Failed]];
  loaded = Quiet @ Check[
    With[{path = $mathPath}, ParallelEvaluate[Get[path]; MathMSSA`$MathMSSALogging = False; MemberQ[$Packages, "MathMSSA`"]]],
    $Failed
  ];
  If[loaded === $Failed || !AllTrue[Flatten[{loaded}], TrueQ], $Failed, Null]
];

evaluateJob[job_Association] := If[
  AllTrue[{"y", "X", "L", "Ranks"}, KeyExistsQ[job, #] &],
  MathMSSA`MathMSSA[job["y"], job["X"], job["L"], job["Ranks"]],
  Failure["MathMSSAParallelInvalidJob", <|"Stage" -> "INPUT"|>]
];

evaluateJob[_] := Failure["MathMSSAParallelInvalidJob", <|"Stage" -> "INPUT"|>];

execution[mode_String, results_List] := <|"ExecutionMode" -> mode, "Results" -> results|>;

MathMSSAParallel[jobs_List] := Module[{count, loaded, results},
  If[jobs === {}, Return[execution["SERIAL", {}]]];
  Block[{MathMSSA`$MathMSSALogging = False},
    count = ensurePool[];
    If[count < 1, Return[execution["SERIAL", Map[evaluateJob, jobs]]]];
    loaded = loadWorkers[];
    If[loaded === $Failed, Return[execution["SERIAL", Map[evaluateJob, jobs]]]];
    results = ParallelTable[
      If[
        AssociationQ[job] && AllTrue[{"y", "X", "L", "Ranks"}, KeyExistsQ[job, #] &],
        MathMSSA`MathMSSA[job["y"], job["X"], job["L"], job["Ranks"]],
        Failure["MathMSSAParallelInvalidJob", <|"Stage" -> "INPUT"|>]
      ],
      {job, jobs},
      Method -> "FinestGrained",
      DistributedContexts -> None,
      ProgressReporting -> False
    ];
    execution["PARALLEL", results]
  ]
];

MathMSSAParallel[___] := Failure["MathMSSAParallelInvalidArity", <|"Stage" -> "INPUT"|>];

End[]

EndPackage[]
