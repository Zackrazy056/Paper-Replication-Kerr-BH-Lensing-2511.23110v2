(* =================================================================== *)
(* Single-case strict SR(k=2): a=0, theta_o=30deg                      *)
(* Uses spherical-basis reduction directly for Schwarzschild case.      *)
(* =================================================================== *)

AppendTo[$Path, "C:/Users/97747/SpinWeightedSpheroidalHarmonics/Kernel"];
Needs["SpinWeightedSpheroidalHarmonics`"];

sWeight = -2;
mIndex = 2;
workingPrec = 60;
kSR = 2;
lModeMax = 60;
theta30 = N[Pi/6, workingPrec];

repoDir = DirectoryName[$InputFileName] <> "/..";
dataFile = ExpandFileName[FileNameJoin[{repoDir, "data", "fig4_data.json"}]];
outFile = ExpandFileName[FileNameJoin[{repoDir, "data", "fig56_sr_spherical_a0_theta30.json"}]];

CleanNum[x_] := Module[{s, v},
  s = If[StringQ[x], x, ToString[x, InputForm]];
  s = StringReplace[s, {" " -> ""}];
  s = StringReplace[s, RegularExpression["`[-+0-9.]*"] -> ""];
  s = StringReplace[s, "\"" -> ""];
  If[StringLength[s] == 0 || StringContainsQ[s, "Indeterminate"] || StringContainsQ[s, "Infinity"], Return[0]];
  v = Quiet[Check[ToExpression[s], 0]];
  N[v, workingPrec]
];

RoundA[a_] := Round[N[a, 20], 10^-3];
RoundW[w_] := Round[N[w, 20], 10^-10];

ReduceOnce[coeff_Association, m_:2] := Module[
  {lVals, lMax, out, fl, fp, fm, termDiag, aNum, aDen, bNum, bDen, aL, bL},
  lVals = Sort[Keys[coeff]];
  If[Length[lVals] == 0, Return[<||>]];
  lMax = Max[lVals];
  out = <||>;
  Do[
    If[l >= lMax, Continue[]];
    fl = Lookup[coeff, l, 0];
    fp = Lookup[coeff, l + 1, 0];
    fm = Lookup[coeff, l - 1, 0];
    termDiag = 1 - (2 m)/(l (l + 1));
    aNum = (l - 1) (l + 3) (l - m + 1) (l + m + 1);
    aDen = 4 l (l + 2) + 3;
    aL = If[aNum >= 0, Sqrt[N[aNum/aDen, workingPrec]], 0];
    bNum = (l^2 - m^2) (l^2 - 4);
    bDen = 4 l^2 - 1;
    bL = If[bNum >= 0, Sqrt[N[bNum/bDen, workingPrec]], 0];
    out[l] = termDiag fl - (aL/(l + 1)) fp - (bL/l) fm;
    ,
    {l, 2, lModeMax}
  ];
  out
];

If[!FileExistsQ[dataFile],
  Print["[!] Missing input data: ", dataFile];
  Quit[1];
];

raw = Association /@ Import[dataFile, "JSON"];
rows = Select[
  raw,
  (!KeyExistsQ[#, "valid"] || TrueQ[#["valid"]]) &&
  RoundA[CleanNum[#["a"]]] == 0 &&
  Abs[CleanNum[#["theta_obs"]] - N[Pi/6, 20]] < 10^-8 &
];

If[Length[rows] == 0,
  Print["[!] No valid rows found for a=0, theta=pi/6 in fig4_data.json"];
  Quit[1];
];

omegaKeys = DeleteDuplicates[RoundW[CleanNum[#["omega"]]] & /@ rows];
omegaKeys = Sort[omegaKeys];

Print["[*] Building strict SR rows for a=0, theta=30deg..."];
Print["    omega count: ", Length[omegaKeys], ", lMax=", lModeMax, ", k=", kSR];

srRows = Table[
  Module[{w = wk, modeRows, coeff, red, f30},
    modeRows = Select[rows, RoundW[CleanNum[#["omega"]]] == w &];
    coeff = Association@Table[
      Module[{l, reC, bRatio, sG, pref, bracket},
        l = Round[CleanNum[r["l"]]];
        reC = CleanNum[r["Re_C"]];
        bRatio = CleanNum[r["B_ratio_re"]] + I*CleanNum[r["B_ratio_im"]];
        sG = CleanNum[r["S_gamma"]];
        pref = N[(2 Pi)/(I w), workingPrec];
        bracket = (-1)^(l + 1) (reC/(16 w^4)) bRatio - 1;
        l -> N[pref sG bracket, workingPrec]
      ],
      {r, modeRows}
    ];
    red = coeff;
    Do[red = ReduceOnce[red, mIndex], {k, 1, kSR}];
    f30 = N[
      Total@Table[
        Lookup[red, l, 0] * N[SpinWeightedSphericalHarmonicY[sWeight, l, mIndex, theta30, 0], workingPrec],
        {l, 2, lModeMax}
      ]/(1 - Cos[theta30])^kSR,
      workingPrec
    ];
    <|
      "a" -> 0.0,
      "omega" -> N[w, 20],
      "k_sr" -> kSR,
      "l_mode_max" -> lModeMax,
      "f_sr_theta30_re" -> ToString[N[Re[f30], workingPrec], InputForm],
      "f_sr_theta30_im" -> ToString[N[Im[f30], workingPrec], InputForm]
    |>
  ],
  {wk, omegaKeys}
];

Export[outFile, srRows, "JSON"];
Print["[+] Export complete: data/fig56_sr_spherical_a0_theta30.json"];
