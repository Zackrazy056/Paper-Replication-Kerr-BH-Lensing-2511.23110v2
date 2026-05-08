(* =================================================================== *)
(* Custom Teukolsky data generation (targeted ranges / configurable)   *)
(* Usage example:
   wolframscript -file generate_teukolsky_data_custom.wls \
     --output ../data/out.json \
     --alist 0.999 \
     --mlist 2 \
     --wmin 0.7 --wmax 1.05 --wn 36 \
     --lmin 2 --lmax 120 --wp 100 --retry 140
*)
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

args = Rest[$ScriptCommandLine];
opts = Association[];
iArg = 1;
While[iArg <= Length[args],
  arg = args[[iArg]];
  If[StringStartsQ[arg, "--"],
    key = StringDrop[arg, 2];
    If[iArg < Length[args] && !StringStartsQ[args[[iArg + 1]], "--"],
      opts[key] = args[[iArg + 1]];
      iArg += 2;
      Continue[];
      ,
      opts[key] = "1";
    ];
  ];
  iArg++;
];

GetOpt[name_, default_] := If[KeyExistsQ[opts, name], opts[name], default];
ToNum[x_] := ToExpression[x];
ParseList[str_] := ToExpression /@ StringSplit[StringReplace[str, " " -> ""], ","];

workingPrec = ToNum[GetOpt["wp", "80"]];
retryPrec = ToNum[GetOpt["retry", "120"]];
lMin = ToNum[GetOpt["lmin", "2"]];
lMax = ToNum[GetOpt["lmax", "25"]];
wMin = ToNum[GetOpt["wmin", "0.1"]];
wMax = ToNum[GetOpt["wmax", "2.0"]];
wN = ToNum[GetOpt["wn", "50"]];
aList = ParseList[GetOpt["alist", "0,0.8,0.99,0.999,-0.99"]];
mList = ParseList[GetOpt["mlist", "2,-2"]];
outputPath = GetOpt["output", "../data/teukolsky_data_custom.json"];

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

  If[!ratioOk, ratioRe = 0; ratioIm = 0;];

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

If[wN <= 1,
  omegaList = {wMin},
  omegaList = Table[wMin + i * (wMax - wMin)/(wN - 1), {i, 0, wN - 1}]
];

Print["[*] Custom generation config:"];
Print["    aList=", aList];
Print["    mList=", mList];
Print["    omega=[", N[wMin, 16], ", ", N[wMax, 16], "] n=", wN];
Print["    l=[", lMin, ", ", lMax, "]"];
Print["    wp=", workingPrec, ", retry=", retryPrec];
Print["    output=", outputPath];

dataTable = Flatten[
  Table[
    GetTeukolskyData[l, m, a, w],
    {a, aList}, {m, mList}, {w, omegaList}, {l, lMin, lMax}
  ]
];

validRows = Select[dataTable, TrueQ[#["valid"]] &];
Print["[*] Total rows: ", Length[dataTable]];
Print["[*] Valid radial rows: ", Length[validRows], "/", Length[dataTable]];
Print["[*] Valid by m: ", Counts[validRows[[All, "m"]]]];

Export[ExpandFileName[FileNameJoin[{DirectoryName[$InputFileName], outputPath}]], dataTable, "JSON"];
Print["[+] Export complete: ", outputPath];
