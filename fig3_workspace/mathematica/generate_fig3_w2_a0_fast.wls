(* =================================================================== *)
(* Kerr GW lensing: fast single-case data for Fig.3                    *)
(* Case: a = 0, M omega = 2.0                                          *)
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

sWeight = -2;
mIndex = 2;
workingPrec = 40;
retryPrec = 60;

(* Faster grid: 1 degree spacing *)
thetaList = N[Range[40, 180, 1] * Pi / 180, workingPrec];

IsFiniteNumeric[x_] := Module[{y},
  y = Quiet[Check[N[x, workingPrec], $Failed]];
  y =!= $Failed && NumberQ[y] && FreeQ[y, ComplexInfinity | DirectedInfinity[_] | Indeterminate | Infinity]
];

SafeToString[x_] := Module[{y},
  y = Quiet[Check[N[x, workingPrec], $Failed]];
  If[y === $Failed || !IsFiniteNumeric[y], "0", ToString[y, InputForm]]
];

ComputeHarmonics[l_, c_, prec_] := Module[
  {sh, s0, sTheta, sPiMinusTheta},
  sh = Quiet[
    Check[
      SpinWeightedSpheroidalHarmonicS[
        sWeight, l, mIndex, SetPrecision[c, prec],
        Method -> {"SphericalExpansion", "NumTerms" -> 24}
      ],
      $Failed
    ]
  ];
  If[sh === $Failed, Return[{False, 0, {}, {}}]];

  s0 = Quiet[Check[N[sh[0, 0], prec], $Failed]];
  sTheta = Quiet[Check[Table[N[sh[th, 0], prec], {th, thetaList}], $Failed]];
  sPiMinusTheta = Quiet[Check[Table[N[sh[Pi - th, 0], prec], {th, thetaList}], $Failed]];

  If[
    s0 === $Failed || sTheta === $Failed || sPiMinusTheta === $Failed ||
      !IsFiniteNumeric[s0] || !AllTrue[sTheta, IsFiniteNumeric] || !AllTrue[sPiMinusTheta, IsFiniteNumeric],
    {False, 0, {}, {}},
    {True, s0, sTheta, sPiMinusTheta}
  ]
];

ComputeEigenReC[l_, c_, omega_, prec_] := Module[
  {eval, reCExpr, reC},
  eval = Quiet[Check[N[SpinWeightedSpheroidalEigenvalue[sWeight, l, mIndex, SetPrecision[c, prec]], prec], $Failed]];
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

GetTeukolskyDataFast[l_, a_, omega_] := Module[
  {
    aPrec, omegaPrec, c,
    harmTry1, harmTry2, eigTry1, eigTry2, ratioTry1, ratioTry2,
    harmOk, eigOk, ratioOk, s0, sTheta, sPiMinusTheta, reC, ratioRe, ratioIm, validAll
  },
  aPrec = SetPrecision[a, workingPrec];
  omegaPrec = SetPrecision[omega, workingPrec];
  c = SetPrecision[aPrec * omegaPrec, workingPrec];

  harmTry1 = Quiet[Check[ComputeHarmonics[l, c, workingPrec], {False, 0, {}, {}}]];
  harmOk = harmTry1[[1]];
  s0 = harmTry1[[2]]; sTheta = harmTry1[[3]]; sPiMinusTheta = harmTry1[[4]];
  If[!harmOk,
    harmTry2 = Quiet[Check[ComputeHarmonics[l, c, retryPrec], {False, 0, {}, {}}]];
    harmOk = harmTry2[[1]];
    s0 = harmTry2[[2]]; sTheta = harmTry2[[3]]; sPiMinusTheta = harmTry2[[4]];
  ];

  eigTry1 = Quiet[Check[ComputeEigenReC[l, c, omegaPrec, workingPrec], {False, 0}]];
  eigOk = eigTry1[[1]]; reC = eigTry1[[2]];
  If[!eigOk,
    eigTry2 = Quiet[Check[ComputeEigenReC[l, c, omegaPrec, retryPrec], {False, 0}]];
    eigOk = eigTry2[[1]]; reC = eigTry2[[2]];
  ];

  ratioTry1 = Quiet[Check[ComputeBratio[l, aPrec, omegaPrec, workingPrec], {False, 0, 0}]];
  ratioOk = ratioTry1[[1]]; ratioRe = ratioTry1[[2]]; ratioIm = ratioTry1[[3]];
  If[!ratioOk,
    ratioTry2 = Quiet[Check[ComputeBratio[l, aPrec, omegaPrec, retryPrec], {False, 0, 0}]];
    ratioOk = ratioTry2[[1]]; ratioRe = ratioTry2[[2]]; ratioIm = ratioTry2[[3]];
  ];

  validAll = harmOk && eigOk && ratioOk;
  If[!validAll,
    s0 = 0; sTheta = ConstantArray[0, Length[thetaList]]; sPiMinusTheta = ConstantArray[0, Length[thetaList]];
    reC = 0; ratioRe = 0; ratioIm = 0;
  ];

  <|
    "l" -> l,
    "a" -> N[aPrec, 20],
    "omega" -> N[omegaPrec, 20],
    "Re_C" -> SafeToString[reC],
    "B_ratio_re" -> SafeToString[ratioRe],
    "B_ratio_im" -> SafeToString[ratioIm],
    "S_0" -> SafeToString[s0],
    "S_theta" -> (SafeToString /@ sTheta),
    "S_pi_minus_theta" -> (SafeToString /@ sPiMinusTheta),
    "valid" -> validAll
  |>
];

Print["[*] Starting fast single-case generation (a=0, Momega=2.0)..."];

aList = {0};
omegaList = {2};
lMax = 60;

dataTable = Flatten[
  Table[
    GetTeukolskyDataFast[l, a, w],
    {a, aList},
    {w, omegaList},
    {l, 2, lMax}
  ]
];

Print["[*] Total rows: ", Length[dataTable]];
Print["[*] Valid rows: ", Count[dataTable[[All, "valid"]], True], "/", Length[dataTable]];

Export[
  ExpandFileName[
    FileNameJoin[{DirectoryName[$InputFileName], "..", "data", "fig3_data_w2_a0_fast.json"}]
  ],
  dataTable,
  "JSON"
];

Print["[+] Export complete: data/fig3_data_w2_a0_fast.json"];
