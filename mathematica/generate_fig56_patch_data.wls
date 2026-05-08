(* =================================================================== *)
(* Kerr GW lensing: patch data for Table 1 and Fig.6                   *)
(* Computes only missing combinations relative to existing fig4_data:   *)
(*   1) theta_o = pi/3 for a = 0,0.8,0.99,0.999,-0.99                  *)
(*   2) theta_o = pi/6 for a = -0.99                                    *)
(* =================================================================== *)

AppendTo[$Path, "C:/Users/97747/SpinWeightedSpheroidalHarmonics/Kernel"];
AppendTo[$Path, "C:/Users/97747/KerrGeodesics/Kernel"];
AppendTo[$Path, "C:/Users/97747/Teukolsky/Kernel"];

Needs["SpinWeightedSpheroidalHarmonics`"];
Needs["KerrGeodesics`"];
Needs["Teukolsky`"];

If[
  Length[DownValues[SpinWeightedSpheroidalEigenvalue]] == 0 ||
  Length[DownValues[KerrGeoOrbit]] == 0 ||
  Length[DownValues[TeukolskyRadial]] == 0,
  Print["[!] BHPToolkit package load failed. Check local paths and dependencies."];
  Quit[1];
];

Print["[*] Initializing parallel kernels..."];
LaunchKernels[];
ParallelEvaluate[
  AppendTo[$Path, "C:/Users/97747/SpinWeightedSpheroidalHarmonics/Kernel"];
  AppendTo[$Path, "C:/Users/97747/KerrGeodesics/Kernel"];
  AppendTo[$Path, "C:/Users/97747/Teukolsky/Kernel"];
  Needs["SpinWeightedSpheroidalHarmonics`"];
  Needs["KerrGeodesics`"];
  Needs["Teukolsky`"];
];
Print["[*] Active kernels: ", $KernelCount];

sWeight = -2;
mIndex = 2;
workingPrec = 45;
retryPrec = 70;
gammaVal = N[0, workingPrec];

wMin = N[10^-2, workingPrec];
wMax = N[2, workingPrec];
wN = 56;
omegaList = N[Exp[Subdivide[Log[wMin], Log[wMax], wN - 1]], workingPrec];

lMin = 2;
lMax = 60;

jobs = {
  <|
    "thetaObs" -> N[Pi/3, workingPrec],
    "aList" -> {0, 8/10, 99/100, 999/1000, -99/100}
  |>,
  <|
    "thetaObs" -> N[Pi/6, workingPrec],
    "aList" -> {-99/100}
  |>
};

IsFiniteNumeric[x_] := Module[{y},
  y = Quiet[Check[N[x, workingPrec], $Failed]];
  y =!= $Failed && NumberQ[y] &&
    FreeQ[y, ComplexInfinity | DirectedInfinity[_] | Indeterminate | Infinity]
];

SafeToString[x_] := Module[{y},
  y = Quiet[Check[N[x, workingPrec], $Failed]];
  If[y === $Failed || !IsFiniteNumeric[y], "0", ToString[y, InputForm]]
];

ComputeHarmonics[l_, c_, thetaObs_, prec_] := Module[
  {sh, sGamma, sTheta, sPiMinusTheta},
  sh = Quiet[
    Check[
      SpinWeightedSpheroidalHarmonicS[
        sWeight, l, mIndex, SetPrecision[c, prec],
        Method -> {"SphericalExpansion", "NumTerms" -> 24}
      ],
      $Failed
    ]
  ];
  If[sh === $Failed, Return[{False, 0, 0, 0}]];

  sGamma = Quiet[Check[N[sh[gammaVal, 0], prec], $Failed]];
  sTheta = Quiet[Check[N[sh[thetaObs, 0], prec], $Failed]];
  sPiMinusTheta = Quiet[Check[N[sh[Pi - thetaObs, 0], prec], $Failed]];

  If[
    sGamma === $Failed || sTheta === $Failed || sPiMinusTheta === $Failed ||
      !IsFiniteNumeric[sGamma] || !IsFiniteNumeric[sTheta] || !IsFiniteNumeric[sPiMinusTheta],
    {False, 0, 0, 0},
    {True, sGamma, sTheta, sPiMinusTheta}
  ]
];

ComputeEigenReC[l_, c_, omega_, prec_] := Module[
  {eval, reCExpr, reC},
  eval = Quiet[
    Check[
      N[SpinWeightedSpheroidalEigenvalue[sWeight, l, mIndex, SetPrecision[c, prec]], prec],
      $Failed
    ]
  ];
  If[eval === $Failed || !IsFiniteNumeric[eval], Return[{False, 0}]];

  reCExpr = eval^2 * (eval + 2)^2 - 8 * SetPrecision[c, prec]^2 * eval * (eval + 2) + 144 * SetPrecision[omega, prec]^2;
  reC = Quiet[Check[N[Re[Sqrt[reCExpr]], prec], $Failed]];

  If[reC === $Failed || !IsFiniteNumeric[reC], {False, 0}, {True, reC}]
];

ComputeBratio[l_, aPrec_, omegaPrec_, prec_] := Module[
  {radialIn, amps, incAmp, reflAmp, ratio},
  radialIn = Quiet[
    Check[
      TeukolskyRadial[
        sWeight, l, mIndex, SetPrecision[aPrec, prec], SetPrecision[omegaPrec, prec],
        "BoundaryConditions" -> "In",
        WorkingPrecision -> prec
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
      !IsFiniteNumeric[incAmp] || !IsFiniteNumeric[reflAmp] || incAmp == 0,
    Return[{False, 0, 0}]
  ];

  ratio = reflAmp / incAmp;
  If[!IsFiniteNumeric[ratio], {False, 0, 0}, {True, Re[ratio], Im[ratio]}]
];

GetRow[l_, a_, omega_, thetaObs_] := Module[
  {
    aPrec, omegaPrec, c,
    harmTry1, harmTry2, eigTry1, eigTry2, ratioTry1, ratioTry2,
    harmOk, eigOk, ratioOk,
    sGamma, sTheta, sPiMinusTheta, reC, ratioRe, ratioIm, validAll
  },
  aPrec = SetPrecision[a, workingPrec];
  omegaPrec = SetPrecision[omega, workingPrec];
  c = SetPrecision[aPrec * omegaPrec, workingPrec];

  harmTry1 = Quiet[Check[ComputeHarmonics[l, c, thetaObs, workingPrec], {False, 0, 0, 0}]];
  harmOk = harmTry1[[1]];
  sGamma = harmTry1[[2]];
  sTheta = harmTry1[[3]];
  sPiMinusTheta = harmTry1[[4]];
  If[!harmOk,
    harmTry2 = Quiet[Check[ComputeHarmonics[l, c, thetaObs, retryPrec], {False, 0, 0, 0}]];
    harmOk = harmTry2[[1]];
    sGamma = harmTry2[[2]];
    sTheta = harmTry2[[3]];
    sPiMinusTheta = harmTry2[[4]];
  ];

  eigTry1 = Quiet[Check[ComputeEigenReC[l, c, omegaPrec, workingPrec], {False, 0}]];
  eigOk = eigTry1[[1]];
  reC = eigTry1[[2]];
  If[!eigOk,
    eigTry2 = Quiet[Check[ComputeEigenReC[l, c, omegaPrec, retryPrec], {False, 0}]];
    eigOk = eigTry2[[1]];
    reC = eigTry2[[2]];
  ];

  ratioTry1 = Quiet[Check[ComputeBratio[l, aPrec, omegaPrec, workingPrec], {False, 0, 0}]];
  ratioOk = ratioTry1[[1]];
  ratioRe = ratioTry1[[2]];
  ratioIm = ratioTry1[[3]];
  If[!ratioOk,
    ratioTry2 = Quiet[Check[ComputeBratio[l, aPrec, omegaPrec, retryPrec], {False, 0, 0}]];
    ratioOk = ratioTry2[[1]];
    ratioRe = ratioTry2[[2]];
    ratioIm = ratioTry2[[3]];
  ];

  validAll = harmOk && eigOk && ratioOk;
  If[!validAll,
    sGamma = 0; sTheta = 0; sPiMinusTheta = 0; reC = 0; ratioRe = 0; ratioIm = 0;
  ];

  <|
    "l" -> l,
    "m" -> mIndex,
    "a" -> N[aPrec, 20],
    "omega" -> N[omegaPrec, 20],
    "theta_obs" -> N[thetaObs, 20],
    "gamma" -> N[gammaVal, 20],
    "Re_C" -> SafeToString[reC],
    "B_ratio_re" -> SafeToString[ratioRe],
    "B_ratio_im" -> SafeToString[ratioIm],
    "S_gamma" -> SafeToString[sGamma],
    "S_theta_obs" -> SafeToString[sTheta],
    "S_pi_minus_theta" -> SafeToString[sPiMinusTheta],
    "valid" -> validAll
  |>
];

DistributeDefinitions[
  sWeight, mIndex, workingPrec, retryPrec, gammaVal, omegaList, lMin, lMax, jobs,
  IsFiniteNumeric, SafeToString, ComputeHarmonics, ComputeEigenReC, ComputeBratio, GetRow
];

Print["[*] Starting patch generation for Table1/Fig6..."];
Print["    jobs: theta=pi/3 (all a), theta=pi/6 (a=-0.99)"];
Print["    omega grid: [", N[wMin, 16], ", ", N[wMax, 16], "], n=", wN];
Print["    l range: [", lMin, ", ", lMax, "]"];

dataTable = Flatten[
  ParallelTable[
    With[{thetaObs = job["thetaObs"], aVals = job["aList"]},
      Table[
        GetRow[l, a, w, thetaObs],
        {a, aVals},
        {w, omegaList},
        {l, lMin, lMax}
      ]
    ],
    {job, jobs},
    Method -> "FinestGrained"
  ]
];

Print["[*] Total rows: ", Length[dataTable]];
Print["[*] Valid rows: ", Count[dataTable[[All, "valid"]], True], "/", Length[dataTable]];

Export[
  ExpandFileName[
    FileNameJoin[{DirectoryName[$InputFileName], "..", "data", "fig56_patch_data.json"}]
  ],
  dataTable,
  "JSON"
];

Print["[+] Export complete: data/fig56_patch_data.json"];

