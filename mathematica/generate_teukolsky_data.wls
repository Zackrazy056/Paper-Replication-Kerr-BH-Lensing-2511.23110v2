(* =================================================================== *)
(* Kerr GW lensing: Teukolsky data generation script                   *)
(* =================================================================== *)

(* Local toolkit paths *)
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
workingPrec = 80;
retryPrec = 120;

IsFiniteNumeric[x_] := Module[{y},
  y = Quiet[Check[N[x, workingPrec], $Failed]];
  y =!= $Failed && FreeQ[y, ComplexInfinity | DirectedInfinity[_] | Indeterminate | Infinity]
];

SafeToString[x_] := Module[{y},
  y = Quiet[Check[N[x, workingPrec], $Failed]];
  If[
    y === $Failed || !IsFiniteNumeric[y],
    "0",
    ToString[y, InputForm]
  ]
];

ComputeBratio[l_, m_, aPrec_, omegaPrec_, prec_] := Module[
  {radialIn, amps, incAmp, reflAmp, ratio},
  radialIn = TeukolskyRadial[
    sWeight, l, m, aPrec, omegaPrec,
    "BoundaryConditions" -> "In",
    WorkingPrecision -> prec,
    AccuracyGoal -> Floor[0.6 prec],
    PrecisionGoal -> Floor[0.6 prec]
  ];
  amps = radialIn["Amplitudes"];
  incAmp = N[amps["Incidence"], prec];
  reflAmp = N[amps["Reflection"], prec];

  If[
    !IsFiniteNumeric[incAmp] || !IsFiniteNumeric[reflAmp] || incAmp == 0,
    {False, 0, 0},
    ratio = reflAmp / incAmp;
    If[
      !IsFiniteNumeric[ratio],
      {False, 0, 0},
      {True, Re[ratio], Im[ratio]}
    ]
  ]
];

GetTeukolskyData[l_, m_, a_, omega_] := Module[
  {
    aPrec, omegaPrec, c, eval, sh, S0, Spi, ReCexpr, ReC,
    ratioOk, ratioRe, ratioIm, ratioTry1, ratioTry2
  },
  aPrec = SetPrecision[a, workingPrec];
  omegaPrec = SetPrecision[omega, workingPrec];
  c = SetPrecision[aPrec * omegaPrec, workingPrec];

  sh = SpinWeightedSpheroidalHarmonicS[sWeight, l, m, c];
  S0 = Quiet[Check[N[sh[0, 0], workingPrec], 0]];
  Spi = Quiet[Check[N[sh[Pi, 0], workingPrec], 0]];

  eval = Quiet[Check[N[SpinWeightedSpheroidalEigenvalue[sWeight, l, m, c], workingPrec], 0]];
  ReCexpr =
    ((eval + 2)^2 + 4 * aPrec * m * omegaPrec - 4 * aPrec^2 * omegaPrec^2) *
    (eval^2 + 36 * aPrec * m * omegaPrec - 36 * aPrec^2 * omegaPrec^2) +
    (2 * eval + 3) * (96 * aPrec^2 * omegaPrec^2 - 48 * aPrec * m * omegaPrec) +
    144 * omegaPrec^2 * (1 - aPrec^2);
  ReC = Quiet[Check[N[Re[Sqrt[ReCexpr]], workingPrec], 0]];
  If[!IsFiniteNumeric[ReC], ReC = 0];

  ratioTry1 = Quiet[Check[ComputeBratio[l, m, aPrec, omegaPrec, workingPrec], {False, 0, 0}]];
  ratioOk = ratioTry1[[1]];
  ratioRe = ratioTry1[[2]];
  ratioIm = ratioTry1[[3]];

  If[!ratioOk,
    ratioTry2 = Quiet[Check[ComputeBratio[l, m, aPrec, omegaPrec, retryPrec], {False, 0, 0}]];
    ratioOk = ratioTry2[[1]];
    ratioRe = ratioTry2[[2]];
    ratioIm = ratioTry2[[3]];
  ];

  If[!ratioOk,
    ratioRe = 0;
    ratioIm = 0;
  ];

  <|
    "l" -> l,
    "m" -> m,
    "a" -> N[aPrec, 20],
    "omega" -> N[omegaPrec, 20],
    "Re_C" -> SafeToString[ReC],
    "B_ratio_re" -> SafeToString[ratioRe],
    "B_ratio_im" -> SafeToString[ratioIm],
    "S_0" -> SafeToString[S0],
    "S_pi" -> SafeToString[Spi],
    "valid" -> ratioOk
  |>
];

Print["[*] Start high-precision data generation..."];

(* Exact parameter grids to avoid precision downgrade warnings *)
aList = {0, 4/5, 99/100, 999/1000, -99/100};
mList = {2, -2};
omegaList = Table[1/10 + i*(19/490), {i, 0, 49}];  (* 50 points in [0.1, 2.0] *)
lMax = 25;

dataTable = Flatten[
  Table[
    GetTeukolskyData[l, m, a, w],
    {a, aList}, {m, mList}, {w, omegaList}, {l, 2, lMax}
  ]
];

validRows = Select[dataTable, TrueQ[#["valid"]] &];
Print["[*] Total rows: ", Length[dataTable]];
Print["[*] Valid radial rows: ", Length[validRows], "/", Length[dataTable]];
Print["[*] Valid by m: ", Counts[validRows[[All, "m"]]]];

Export[
  ExpandFileName[
    FileNameJoin[{DirectoryName[$InputFileName], "..", "data", "teukolsky_data_real.json"}]
  ],
  dataTable,
  "JSON"
];

Print["[+] Export complete: data/teukolsky_data_real.json"];

