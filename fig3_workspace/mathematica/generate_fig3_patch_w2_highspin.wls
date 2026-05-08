(* =================================================================== *)
(* Kerr GW lensing: targeted patch for Fig.3 paper-grid modes           *)
(* Default cases: M omega = 0.6, 2; a = 0 and a = 0.99                   *)
(* =================================================================== *)

EnvOrDefault[name_, default_] := Module[{value},
  value = Environment[name];
  If[StringQ[value] && StringLength[StringTrim[value]] > 0, StringTrim[value], default]
];

swshPath = EnvOrDefault["BHP_SW_SH_PATH", "C:/Users/97747/SpinWeightedSpheroidalHarmonics/Kernel"];
kerrGeodesicsPath = EnvOrDefault["BHP_KERR_GEODESICS_PATH", "C:/Users/97747/KerrGeodesics/Kernel"];
teukolskyPath = EnvOrDefault["BHP_TEUKOLSKY_PATH", "C:/Users/97747/Teukolsky/Kernel"];

AppendTo[$Path, swshPath];
AppendTo[$Path, kerrGeodesicsPath];
AppendTo[$Path, teukolskyPath];

Needs["SpinWeightedSpheroidalHarmonics`"];
Needs["KerrGeodesics`"];
Needs["Teukolsky`"];

Print["[*] Initializing Parallel Kernels..."];
kernelRaw = Environment["FIG3_KERNELS"];
kernelTarget = If[
  StringQ[kernelRaw] && StringLength[StringTrim[kernelRaw]] > 0,
  Quiet[Check[ToExpression[StringTrim[kernelRaw]], 0]],
  0
];
If[IntegerQ[kernelTarget] && kernelTarget > 0,
  LaunchKernels[kernelTarget],
  LaunchKernels[]
];
ParallelEvaluate[
  sw = Environment["BHP_SW_SH_PATH"];
  kg = Environment["BHP_KERR_GEODESICS_PATH"];
  tk = Environment["BHP_TEUKOLSKY_PATH"];
  If[!StringQ[sw] || StringLength[StringTrim[sw]] == 0, sw = "C:/Users/97747/SpinWeightedSpheroidalHarmonics/Kernel"];
  If[!StringQ[kg] || StringLength[StringTrim[kg]] == 0, kg = "C:/Users/97747/KerrGeodesics/Kernel"];
  If[!StringQ[tk] || StringLength[StringTrim[tk]] == 0, tk = "C:/Users/97747/Teukolsky/Kernel"];
  AppendTo[$Path, StringTrim[sw]];
  AppendTo[$Path, StringTrim[kg]];
  AppendTo[$Path, StringTrim[tk]];
  Needs["SpinWeightedSpheroidalHarmonics`"];
  Needs["KerrGeodesics`"];
  Needs["Teukolsky`"];
];
Print["[*] Active kernels: ", $KernelCount];

If[
  Length[DownValues[SpinWeightedSpheroidalEigenvalue]] == 0 ||
  Length[DownValues[KerrGeoOrbit]] == 0 ||
  Length[DownValues[TeukolskyRadial]] == 0,
  Print["[!] BHPToolkit package load failed. Check local paths and dependencies."];
  Quit[1];
];

sWeight = -2;
mIndex = 2;
workingPrec = 80;
retryPrec1 = 120;
retryPrec2 = 180;
harmRelTol = 10^-12;
bratioRelTol = 10^-8;
checkBratio = False;

thetaList = N[Range[40, 180, 1/2] * Pi / 180, workingPrec];

IsFiniteNumeric[x_, prec_: workingPrec] := Module[{y},
  y = Quiet[Check[N[x, prec], $Failed]];
  y =!= $Failed && NumberQ[y] && FreeQ[y, ComplexInfinity | DirectedInfinity[_] | Indeterminate | Infinity]
];

SafeToString[x_, prec_: workingPrec] := Module[{y},
  y = Quiet[Check[N[x, prec], $Failed]];
  If[y === $Failed || !IsFiniteNumeric[y, prec], "0", ToString[y, InputForm]]
];

RelDiff[x_, y_, prec_: workingPrec] := Module[{xx, yy, den},
  xx = N[x, prec];
  yy = N[y, prec];
  den = Max[Abs[xx], Abs[yy], 10^-50];
  N[Abs[xx - yy]/den, prec]
];

MaxRelDiffList[listA_, listB_, prec_: workingPrec] := Module[{n},
  n = Min[Length[listA], Length[listB]];
  If[n == 0, Return[Infinity]];
  Max@Table[RelDiff[listA[[k]], listB[[k]], prec], {k, 1, n}]
];

ComputeHarmonicsWithPrec[l_, c_, prec_] := Module[
  {
    numTermsList, nt, sh, s0, sTheta, sPiMinusTheta,
    prevS0 = 0, prevSTheta = {}, prevSPi = {}, hasPrev = False,
    relErr = Infinity, converged = False,
    lastGood
  },
  If[Abs[N[c, prec]] < 10^-40,
    s0 = Quiet[Check[N[SpinWeightedSphericalHarmonicY[sWeight, l, mIndex, 0, 0], prec], $Failed]];
    sTheta = Quiet[Check[Table[N[SpinWeightedSphericalHarmonicY[sWeight, l, mIndex, th, 0], prec], {th, thetaList}], $Failed]];
    sPiMinusTheta = Quiet[Check[Table[N[SpinWeightedSphericalHarmonicY[sWeight, l, mIndex, Pi - th, 0], prec], {th, thetaList}], $Failed]];
    If[
      s0 === $Failed || sTheta === $Failed || sPiMinusTheta === $Failed ||
        !IsFiniteNumeric[s0, prec] ||
        !AllTrue[sTheta, IsFiniteNumeric[#, prec] &] ||
        !AllTrue[sPiMinusTheta, IsFiniteNumeric[#, prec] &],
      Return[{False, 0, {}, {}, 0, Infinity, False}],
      Return[{True, s0, sTheta, sPiMinusTheta, 0, 0, True}]
    ];
  ];

  numTermsList = DeleteDuplicates[{40, 60, Max[80, l + 20]}];
  lastGood = {False, 0, {}, {}, 0, Infinity, False};
  Do[
    sh = Quiet[
      Check[
        SpinWeightedSpheroidalHarmonicS[
          sWeight, l, mIndex, SetPrecision[c, prec],
          Method -> {"SphericalExpansion", "NumTerms" -> nt}
        ],
        $Failed
      ]
    ];
    If[sh === $Failed,
      sh = Quiet[
        Check[
          SpinWeightedSpheroidalHarmonicS[sWeight, l, mIndex, SetPrecision[c, prec]],
          $Failed
        ]
      ];
    ];
    If[sh === $Failed, Continue[]];
    s0 = Quiet[Check[N[sh[0, 0], prec], $Failed]];
    sTheta = Quiet[Check[Table[N[sh[th, 0], prec], {th, thetaList}], $Failed]];
    sPiMinusTheta = Quiet[Check[Table[N[sh[Pi - th, 0], prec], {th, thetaList}], $Failed]];
    If[
      s0 === $Failed || sTheta === $Failed || sPiMinusTheta === $Failed ||
        !IsFiniteNumeric[s0, prec] ||
        !AllTrue[sTheta, IsFiniteNumeric[#, prec] &] ||
        !AllTrue[sPiMinusTheta, IsFiniteNumeric[#, prec] &],
      Continue[]
    ];
    relErr = If[!TrueQ[hasPrev],
      Infinity,
      Max[
        RelDiff[s0, prevS0, prec],
        MaxRelDiffList[sTheta, prevSTheta, prec],
        MaxRelDiffList[sPiMinusTheta, prevSPi, prec]
      ]
    ];
    converged = NumericQ[relErr] && relErr < harmRelTol;
    lastGood = {True, s0, sTheta, sPiMinusTheta, nt, relErr, converged};
    prevS0 = s0;
    prevSTheta = sTheta;
    prevSPi = sPiMinusTheta;
    hasPrev = True;
    If[converged, Break[]];
    ,
    {nt, numTermsList}
  ];
  lastGood
];

ComputeEigenReCWithPrec[l_, c_, omega_, prec_] := Module[
  {eval, reCExpr, reC},
  eval = Quiet[Check[N[SpinWeightedSpheroidalEigenvalue[sWeight, l, mIndex, SetPrecision[c, prec]], prec], $Failed]];
  If[eval === $Failed || !IsFiniteNumeric[eval, prec], Return[{False, 0, 0}]];
  reCExpr = eval^2 * (eval + 2)^2 - 8 * SetPrecision[c, prec]^2 * eval * (eval + 2) + 144 * SetPrecision[omega, prec]^2;
  reC = Quiet[Check[N[Re[Sqrt[reCExpr]], prec], $Failed]];
  If[reC === $Failed || !IsFiniteNumeric[reC, prec], {False, 0, eval}, {True, reC, eval}]
];

ComputeBratioWithOptions[l_, aPrec_, omegaPrec_, prec_, goalFrac_, useGoals_] := Module[
  {radialIn, amps, incAmp, reflAmp, ratio, goal},
  goal = Max[20, Floor[goalFrac * prec]];
  radialIn = Quiet[
    Check[
      If[TrueQ[useGoals],
        TeukolskyRadial[
          sWeight, l, mIndex, SetPrecision[aPrec, prec], SetPrecision[omegaPrec, prec],
          "BoundaryConditions" -> "In",
          WorkingPrecision -> prec,
          AccuracyGoal -> goal,
          PrecisionGoal -> goal
        ],
        TeukolskyRadial[
          sWeight, l, mIndex, SetPrecision[aPrec, prec], SetPrecision[omegaPrec, prec],
          "BoundaryConditions" -> "In",
          WorkingPrecision -> prec
        ]
      ],
      $Failed
    ]
  ];
  If[radialIn === $Failed, Return[{False, 0, 0}]];

  amps = Quiet[Check[radialIn["Amplitudes"], $Failed]];
  If[amps === $Failed, Return[{False, 0, 0}]];

  incAmp = Quiet[Check[N[amps["Incidence"], prec], $Failed]];
  reflAmp = Quiet[Check[N[amps["Reflection"], prec], $Failed]];
  If[
    incAmp === $Failed || reflAmp === $Failed ||
      !IsFiniteNumeric[incAmp, prec] || !IsFiniteNumeric[reflAmp, prec] || incAmp == 0,
    Return[{False, 0, 0}]
  ];

  ratio = reflAmp / incAmp;
  If[!IsFiniteNumeric[ratio, prec], {False, 0, 0}, {True, Re[ratio], Im[ratio]}]
];

GetTeukolskyDataPatch[l_, a_, omega_] := Module[
  {
    aPrec, omegaPrec, c,
    harmResult, eigResult, ratioResult,
    s0 = 0, sTheta = {}, sPiMinusTheta = {},
    reC = 0, evalEig = 0, ratioRe = 0, ratioIm = 0,
    validH = False, validE = False, validR = False,
    validAll,
    harmonicsPrecList, eigenPrecList, radialTryList, pr, gf, ug,
    harmNumTerms = 0, harmRelErr = Infinity, harmConverged = False,
    ratioPrecUsed = 0, ratioGoalFracUsed = 0.0, ratioUseGoalsUsed = False,
    ratioCheckResult, ratioRelCheck = Infinity, ratioConverged = False
  },
  aPrec = SetPrecision[a, workingPrec];
  omegaPrec = SetPrecision[omega, workingPrec];
  c = SetPrecision[aPrec * omegaPrec, workingPrec];

  harmonicsPrecList = {60, 80, 100, 120, 180};
  eigenPrecList = {60, 80, 100, 120, 160, 220, 300};
  radialTryList = {
    {60, 0.50, True}, {80, 0.45, True}, {100, 0.40, True},
    {30, 0.0, False}, {40, 0.0, False}, {50, 0.0, False}, {60, 0.0, False},
    {70, 0.0, False}, {80, 0.0, False}, {100, 0.0, False},
    {120, 0.35, True}, {160, 0.30, True},
    {200, 0.40, True}, {240, 0.45, True}, {300, 0.50, True}
  };

  Do[
    harmResult = Quiet[Check[ComputeHarmonicsWithPrec[l, c, pr], {False, 0, {}, {}, 0, Infinity, False}]];
    If[
      harmResult[[1]],
      validH = True;
      s0 = harmResult[[2]];
      sTheta = harmResult[[3]];
      sPiMinusTheta = harmResult[[4]];
      harmNumTerms = harmResult[[5]];
      harmRelErr = harmResult[[6]];
      harmConverged = harmResult[[7]];
      If[TrueQ[harmConverged], Break[];];
    ],
    {pr, harmonicsPrecList}
  ];

  Do[
    eigResult = Quiet[Check[ComputeEigenReCWithPrec[l, c, omegaPrec, pr], {False, 0, 0}]];
    If[eigResult[[1]], validE = True; reC = eigResult[[2]]; evalEig = eigResult[[3]]; Break[];],
    {pr, eigenPrecList}
  ];

  Do[
    pr = pair[[1]]; gf = pair[[2]]; ug = pair[[3]];
    ratioResult = Quiet[Check[ComputeBratioWithOptions[l, aPrec, omegaPrec, pr, gf, ug], {False, 0, 0}]];
    If[
      ratioResult[[1]],
      validR = True;
      ratioRe = ratioResult[[2]];
      ratioIm = ratioResult[[3]];
      ratioPrecUsed = pr;
      ratioGoalFracUsed = gf;
      ratioUseGoalsUsed = ug;
      Break[];
    ],
    {pair, radialTryList}
  ];

  If[validR && TrueQ[checkBratio],
    ratioCheckResult = Quiet[
      Check[
        ComputeBratioWithOptions[
          l, aPrec, omegaPrec,
          Min[ratioPrecUsed + 40, 200],
          Max[0.30, ratioGoalFracUsed],
          True
        ],
        {False, 0, 0}
      ]
    ];
    If[ratioCheckResult[[1]],
      ratioRelCheck = RelDiff[
        ratioRe + I ratioIm,
        ratioCheckResult[[2]] + I ratioCheckResult[[3]],
        Min[ratioPrecUsed + 40, 200]
      ];
      ratioConverged = ratioRelCheck < bratioRelTol;
      ,
      ratioConverged = False;
      ratioRelCheck = Infinity;
    ];
  ];

  validAll = validH && validE && validR;
  If[!validAll,
    s0 = 0;
    sTheta = ConstantArray[0, Length[thetaList]];
    sPiMinusTheta = ConstantArray[0, Length[thetaList]];
    reC = 0;
    ratioRe = 0;
    ratioIm = 0;
  ];

  <|
    "l" -> l,
    "a" -> N[aPrec, 20],
    "omega" -> N[omegaPrec, 20],
    "Re_C" -> SafeToString[reC],
    "sep_eigenvalue" -> SafeToString[evalEig],
    "B_ratio_re" -> SafeToString[ratioRe],
    "B_ratio_im" -> SafeToString[ratioIm],
    "S_0" -> SafeToString[s0],
    "S_theta" -> (SafeToString /@ sTheta),
    "S_pi_minus_theta" -> (SafeToString /@ sPiMinusTheta),
    "harmonics_num_terms" -> harmNumTerms,
    "harmonics_rel_check" -> SafeToString[harmRelErr],
    "harmonics_converged" -> harmConverged,
    "bratio_prec_used" -> ratioPrecUsed,
    "bratio_goal_frac_used" -> ratioGoalFracUsed,
    "bratio_use_goals_used" -> ratioUseGoalsUsed,
    "bratio_rel_check" -> SafeToString[ratioRelCheck],
    "bratio_converged" -> ratioConverged,
    "valid" -> validAll
  |>
];

DistributeDefinitions[
  sWeight, mIndex, workingPrec, thetaList, harmRelTol, bratioRelTol, checkBratio,
  IsFiniteNumeric, SafeToString, RelDiff, MaxRelDiffList,
  ComputeHarmonicsWithPrec, ComputeEigenReCWithPrec, ComputeBratioWithOptions,
  GetTeukolskyDataPatch
];

Print["[*] Starting high-frequency targeted patch generation..."];

aList = {0, 0.99};
omegaList = {0.6, 2};
lMin = 2;
lMax = 90;
batchSize = 20;
outName = "fig3_patch_w2_a0_a099.json";

If[StringQ[Environment["FIG3_PATCH_A_LIST"]] && StringLength[StringTrim[Environment["FIG3_PATCH_A_LIST"]]] > 0,
  aList = Quiet[Check[ToExpression /@ StringSplit[StringTrim[Environment["FIG3_PATCH_A_LIST"]], ","], aList]];
];
If[StringQ[Environment["FIG3_PATCH_OMEGA_LIST"]] && StringLength[StringTrim[Environment["FIG3_PATCH_OMEGA_LIST"]]] > 0,
  omegaList = Quiet[Check[ToExpression /@ StringSplit[StringTrim[Environment["FIG3_PATCH_OMEGA_LIST"]], ","], omegaList]];
];
If[StringQ[Environment["FIG3_PATCH_LMIN"]] && StringLength[StringTrim[Environment["FIG3_PATCH_LMIN"]]] > 0,
  lMin = Quiet[Check[ToExpression[Environment["FIG3_PATCH_LMIN"]], lMin]];
];
If[StringQ[Environment["FIG3_PATCH_LMAX"]] && StringLength[StringTrim[Environment["FIG3_PATCH_LMAX"]]] > 0,
  lMax = Quiet[Check[ToExpression[Environment["FIG3_PATCH_LMAX"]], lMax]];
];
If[StringQ[Environment["FIG3_PATCH_BATCH"]] && StringLength[StringTrim[Environment["FIG3_PATCH_BATCH"]]] > 0,
  batchSize = Quiet[Check[ToExpression[Environment["FIG3_PATCH_BATCH"]], batchSize]];
];
If[StringQ[Environment["FIG3_PATCH_OUT"]] && StringLength[StringTrim[Environment["FIG3_PATCH_OUT"]]] > 0,
  outName = StringTrim[Environment["FIG3_PATCH_OUT"]];
];
If[StringQ[Environment["FIG3_PATCH_CHECK_BRATIO"]] && StringLength[StringTrim[Environment["FIG3_PATCH_CHECK_BRATIO"]]] > 0,
  checkBratio = MemberQ[
    {"1", "true", "yes", "on"},
    ToLowerCase[StringTrim[Environment["FIG3_PATCH_CHECK_BRATIO"]]]
  ];
];
If[StringQ[Environment["FIG3_PATCH_BRATIO_TOL"]] && StringLength[StringTrim[Environment["FIG3_PATCH_BRATIO_TOL"]]] > 0,
  bratioRelTol = Quiet[Check[ToExpression[StringTrim[Environment["FIG3_PATCH_BRATIO_TOL"]]], bratioRelTol]];
];
DistributeDefinitions[bratioRelTol, checkBratio];

If[!IntegerQ[lMin], lMin = 2];
If[!IntegerQ[lMax], lMax = 90];
If[!IntegerQ[batchSize] || batchSize <= 0, batchSize = 20];
If[lMin < 2, lMin = 2];
If[lMax < lMin, lMax = lMin];

outPath = ExpandFileName[
  FileNameJoin[{DirectoryName[$InputFileName], "..", "data", outName}]
];

KeyA[a_] := Round[N[a, 20], 10^-3];
KeyW[w_] := Round[N[w, 20], 10^-10];
KeyRow[a_, w_, l_] := {KeyA[a], KeyW[w], Round[l]};
CompletePatchRowQ[row_] := Module[{sTheta, sPiMinusTheta},
  If[!AssociationQ[row], Return[False]];
  If[KeyExistsQ[row, "valid"] && !TrueQ[row["valid"]], Return[False]];
  If[!(KeyExistsQ[row, "a"] && KeyExistsQ[row, "omega"] && KeyExistsQ[row, "l"]), Return[False]];
  sTheta = If[KeyExistsQ[row, "S_theta"], row["S_theta"], {}];
  sPiMinusTheta = If[KeyExistsQ[row, "S_pi_minus_theta"], row["S_pi_minus_theta"], {}];
  ListQ[sTheta] && ListQ[sPiMinusTheta] &&
    Length[sTheta] == Length[thetaList] &&
    Length[sPiMinusTheta] == Length[thetaList]
];

existingRows = If[FileExistsQ[outPath], Quiet[Check[Import[outPath, "RawJSON"], {}]], {}];
existingAssocRows = Association /@ Select[existingRows, AssociationQ[#] || MatchQ[#, {(_Rule) ..}] &];
existingKeys = DeleteDuplicates[
  KeyRow[#["a"], #["omega"], #["l"]] & /@
    Select[existingAssocRows, CompletePatchRowQ]
];

allTasks = Flatten[
  Table[{a, w, l}, {a, aList}, {w, omegaList}, {l, lMin, lMax}],
  2
];
pendingTasks = Select[allTasks, !MemberQ[existingKeys, KeyRow[#[[1]], #[[2]], #[[3]]]] &];
batchTasks = Take[pendingTasks, UpTo[batchSize]];

Print["[*] Config: aList=", N[aList, 8], ", omegaList=", N[omegaList, 8], ", l=[", lMin, ",", lMax, "], batch=", batchSize];
Print["[*] Output file: ", outPath];
Print["[*] Existing rows: ", Length[existingAssocRows], ", done keys: ", Length[existingKeys]];
Print["[*] Pending tasks: ", Length[pendingTasks], ", running now: ", Length[batchTasks]];

If[Length[batchTasks] == 0,
  Print["[*] Nothing to run. All tasks in range are complete."];
  Quit[0];
];

Print["[*] Starting parallel generation across ", $KernelCount, " kernels..."];
batchRows = ParallelTable[
  GetTeukolskyDataPatch[task[[3]], task[[1]], task[[2]]],
  {task, batchTasks},
  Method -> "FinestGrained"
];

Print["[*] Batch rows: ", Length[batchRows]];
Print["[*] Batch valid: ", Count[batchRows[[All, "valid"]], True], "/", Length[batchRows]];

mergedAssoc = Association@Table[
  KeyRow[row["a"], row["omega"], row["l"]] -> row,
  {row, Join[existingAssocRows, batchRows]}
];
mergedRows = SortBy[
  Values[mergedAssoc],
  {KeyA[#["a"]] &, KeyW[#["omega"]] &, Round[#["l"]] &}
];

Export[outPath, mergedRows, "JSON"];

doneNow = Length[mergedRows];
pendingAfter = Length[allTasks] - doneNow;

Do[
  Module[{rows, validRows, validLs},
    rows = Select[
      mergedRows,
      (Abs[N[#["a"], 20] - N[a, 20]] < 10^-12 && Abs[N[#["omega"], 20] - N[w, 20]] < 10^-12) &
    ];
    validRows = Select[rows, #["valid"] === True &];
    validLs = validRows[[All, "l"]];
    Print[
      "    a=", N[a, 6], ", Momega=", N[w, 6],
      " -> stored ", Length[rows], ", valid ", Length[validRows],
      ", lmax(valid)=", If[Length[validLs] > 0, Max[validLs], "None"]
    ];
  ],
  {a, aList}, {w, omegaList}
];

Print["[+] Export complete: ", outPath];
Print["[*] Progress: total target tasks=", Length[allTasks], ", completed=", doneNow, ", remaining=", pendingAfter];
