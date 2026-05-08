(* =================================================================== *)
(* Fig.7-Fig.9 theta-grid scattering cache                              *)
(*                                                                       *)
(* Produces:                                                             *)
(*   data/fig789_patch_data.json      modal rows with S(theta), S(pi-th) *)
(*   data/fig789_sr_spherical.json    generic f_SR(a, omega, theta) rows *)
(*                                                                       *)
(* This script reuses existing Teukolsky radial data in fig4/fig56 JSON   *)
(* files. It only recomputes angular harmonics and the Appendix-C         *)
(* spherical-basis series-reduction projection for the requested theta    *)
(* grid, avoiding another expensive radial solve.                         *)
(* =================================================================== *)

AppendTo[$Path, "C:/Users/97747/SpinWeightedSpheroidalHarmonics/Kernel"];
Needs["SpinWeightedSpheroidalHarmonics`"];

sWeight = -2;
mIndex = 2;
workingPrec = 60;
kSR = 2;
lModeMax = 60;
lBasisMax = 70;
nThetaFit = 80;

repoDir = DirectoryName[$InputFileName] <> "/..";
dataFiles = {
  ExpandFileName[FileNameJoin[{repoDir, "data", "fig4_data.json"}]],
  ExpandFileName[FileNameJoin[{repoDir, "data", "fig56_patch_data.json"}]],
  ExpandFileName[FileNameJoin[{repoDir, "data", "fig56_patch_a_minus099_theta30_hi.json"}]]
};
modalOutPath = ExpandFileName[FileNameJoin[{repoDir, "data", "fig789_patch_data.json"}]];
srOutPath = ExpandFileName[FileNameJoin[{repoDir, "data", "fig789_sr_spherical.json"}]];

CleanNum[x_] := Module[{s, v},
  s = If[StringQ[x], x, ToString[x, InputForm]];
  s = StringReplace[s, {" " -> ""}];
  s = StringReplace[s, RegularExpression["`[-+0-9.]*"] -> ""];
  s = StringReplace[s, "\"" -> ""];
  s = StringReplace[s, RegularExpression["([0-9.]+)[eE]([+-]?[0-9]+)"] -> "$1*^$2"];
  If[
    StringLength[s] == 0 ||
      StringContainsQ[s, "Indeterminate"] ||
      StringContainsQ[s, "Infinity"] ||
      StringContainsQ[s, "ComplexInfinity"],
    Return[0]
  ];
  v = Quiet[Check[ToExpression[s], 0]];
  If[NumericQ[v], N[v, workingPrec], 0]
];

KeyA[a_] := Round[N[a, 20], 10^-3];
KeyW[w_] := Round[N[w, 20], 10^-10];

IsFiniteNumeric[x_] := Module[{y},
  y = Quiet[Check[N[x, workingPrec], $Failed]];
  y =!= $Failed && NumberQ[y] &&
    FreeQ[y, ComplexInfinity | DirectedInfinity[_] | Indeterminate | Infinity]
];

SafeToString[x_] := Module[{y},
  y = Quiet[Check[N[x, workingPrec], $Failed]];
  If[y === $Failed || !IsFiniteNumeric[y], "0", ToString[y, InputForm]]
];

thetaDegRaw = Environment["FIG789_THETA_DEG_LIST"];
thetaDegList = If[
  StringQ[thetaDegRaw] && StringLength[StringTrim[thetaDegRaw]] > 0,
  Quiet[Check[ToExpression /@ StringSplit[StringTrim[thetaDegRaw], ","], {30, 40, 50, 60, 70, 80, 89}]],
  {30, 40, 50, 60, 70, 80, 89}
];
thetaDegList = N[thetaDegList, 20];
thetaRecords = Table[
  <|"thetaDeg" -> N[deg, 20], "theta" -> N[Pi deg/180, workingPrec]|>,
  {deg, thetaDegList}
];

rawData = Flatten[
  Table[
    If[FileExistsQ[f], Import[f, "RawJSON"], {}],
    {f, dataFiles}
  ],
  1
];

rawAssoc = Association /@ Select[
  rawData,
  AssociationQ[#] || MatchQ[#, {(_Rule) ..}] &
];
validData = Select[
  rawAssoc,
  KeyExistsQ[#, "a"] &&
    KeyExistsQ[#, "omega"] &&
    KeyExistsQ[#, "l"] &&
    KeyExistsQ[#, "Re_C"] &&
    KeyExistsQ[#, "B_ratio_re"] &&
    KeyExistsQ[#, "B_ratio_im"] &&
    KeyExistsQ[#, "S_gamma"] &&
    (!KeyExistsQ[#, "valid"] || TrueQ[#["valid"]]) &
];

modeAssoc = <||>;
Do[
  Module[{a, w, l, k, rec},
    a = CleanNum[row["a"]];
    w = CleanNum[row["omega"]];
    l = Round[CleanNum[row["l"]]];
    If[l < 2 || l > lModeMax, Continue[]];
    k = {KeyA[a], KeyW[w], l};
    rec = <|
      "a" -> a,
      "w" -> w,
      "l" -> l,
      "ReC" -> CleanNum[row["Re_C"]],
      "B" -> (CleanNum[row["B_ratio_re"]] + I*CleanNum[row["B_ratio_im"]]),
      "Sg" -> CleanNum[row["S_gamma"]]
    |>;
    If[!KeyExistsQ[modeAssoc, k], modeAssoc[k] = rec];
  ],
  {row, validData}
];

comboKeysAll = Sort[DeleteDuplicates[(Keys[modeAssoc]) /. {ka_, kw_, _} :> {ka, kw}]];
aFilter = {};
aFilterRaw = Environment["FIG789_A_FILTER"];
If[StringQ[aFilterRaw] && StringLength[StringTrim[aFilterRaw]] > 0,
  aFilter = Quiet[Check[KeyA /@ (ToExpression /@ StringSplit[StringTrim[aFilterRaw], ","]), {}]];
];
comboKeys = comboKeysAll;
If[ListQ[aFilter] && Length[aFilter] > 0,
  comboKeys = Select[comboKeysAll, MemberQ[aFilter, KeyA[#[[1]]]] &];
];

thetaFitGrid = N[Subdivide[0.03, Pi - 0.03, nThetaFit - 1], workingPrec];
yFitMatGlobal = Table[
  N[SpinWeightedSphericalHarmonicY[sWeight, L, mIndex, th, 0], workingPrec],
  {th, thetaFitGrid},
  {L, 2, lBasisMax}
];
yThetaAssocGlobal = Association@Table[
  i -> Association@Table[
    L -> N[SpinWeightedSphericalHarmonicY[sWeight, L, mIndex, thetaRecords[[i, "theta"]], 0], workingPrec],
    {L, 2, lBasisMax}
  ],
  {i, Length[thetaRecords]}
];

ReduceOnce[coeff_Association, m_:2] := Module[
  {lVals, lMin, lMax, out, fl, fp, fm, termDiag, aNum, aDen, bNum, bDen, aL, bL},
  lVals = Sort[Keys[coeff]];
  If[Length[lVals] == 0, Return[<||>]];
  lMin = Min[lVals];
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
    {l, lMin, lMax}
  ];
  out
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
  sGamma = Quiet[Check[N[sh[0, 0], prec], $Failed]];
  sTheta = Quiet[Check[N[sh[thetaObs, 0], prec], $Failed]];
  sPiMinusTheta = Quiet[Check[N[sh[Pi - thetaObs, 0], prec], $Failed]];
  If[
    sGamma === $Failed || sTheta === $Failed || sPiMinusTheta === $Failed ||
      !IsFiniteNumeric[sGamma] || !IsFiniteNumeric[sTheta] || !IsFiniteNumeric[sPiMinusTheta],
    {False, 0, 0, 0},
    {True, sGamma, sTheta, sPiMinusTheta}
  ]
];

ComputeThetaModeRow[thetaIndex_, ka_, kw_, l_] := Module[
  {base, thetaObs, a, w, c, harm, validAll, sGamma, sTheta, sPiMinusTheta},
  If[!KeyExistsQ[modeAssoc, {ka, kw, l}], Return[Nothing]];
  base = modeAssoc[{ka, kw, l}];
  thetaObs = thetaRecords[[thetaIndex, "theta"]];
  a = base["a"];
  w = base["w"];
  c = SetPrecision[a*w, workingPrec];
  harm = Quiet[Check[ComputeHarmonics[l, c, thetaObs, workingPrec], {False, 0, 0, 0}]];
  validAll = harm[[1]];
  sGamma = If[validAll, harm[[2]], 0];
  sTheta = If[validAll, harm[[3]], 0];
  sPiMinusTheta = If[validAll, harm[[4]], 0];
  <|
    "l" -> l,
    "m" -> mIndex,
    "a" -> N[a, 20],
    "omega" -> N[w, 20],
    "theta_obs" -> N[thetaObs, 20],
    "theta_deg" -> N[thetaRecords[[thetaIndex, "thetaDeg"]], 20],
    "gamma" -> 0,
    "Re_C" -> SafeToString[base["ReC"]],
    "B_ratio_re" -> SafeToString[Re[base["B"]]],
    "B_ratio_im" -> SafeToString[Im[base["B"]]],
    "S_gamma" -> SafeToString[sGamma],
    "S_theta_obs" -> SafeToString[sTheta],
    "S_pi_minus_theta" -> SafeToString[sPiMinusTheta],
    "valid" -> validAll
  |>
];

ComputeSRRows[ka_, kw_] := Module[
  {
    a, w, pref, modeRows, coeffA, shAssoc, fFitVec, sol,
    fCoeff, redCoeff, fitResidual, lKeep
  },
  modeRows = Table[
    If[KeyExistsQ[modeAssoc, {ka, kw, l}], modeAssoc[{ka, kw, l}], Missing["NotFound"]],
    {l, 2, lModeMax}
  ];
  modeRows = Select[modeRows, # =!= Missing["NotFound"] &];
  If[Length[modeRows] < 20, Return[{}]];
  a = modeRows[[1, "a"]];
  w = modeRows[[1, "w"]];
  If[w <= 0, Return[{}]];
  pref = N[(2 Pi)/(I w), workingPrec];

  coeffA = Association@Table[
    Module[{l = row["l"], bracket},
      bracket = (-1)^(l + 1) (row["ReC"]/(16 w^4)) row["B"] - 1;
      l -> N[pref row["Sg"] bracket, workingPrec]
    ],
    {row, modeRows}
  ];

  shAssoc = Association@Table[
    l -> SpinWeightedSpheroidalHarmonicS[
      sWeight, l, mIndex, N[a w, workingPrec],
      Method -> {"SphericalExpansion", "NumTerms" -> 28}
    ],
    {l, 2, lModeMax}
  ];

  fFitVec = Table[
    N[
      Total@Table[
        Lookup[coeffA, l, 0] * shAssoc[l][th, 0],
        {l, 2, lModeMax}
      ],
      workingPrec
    ],
    {th, thetaFitGrid}
  ];
  sol = Quiet[Check[LeastSquares[yFitMatGlobal, fFitVec], $Failed]];
  If[sol === $Failed, Return[{}]];
  fCoeff = Association@Table[(L -> N[sol[[L - 1]], workingPrec]), {L, 2, lBasisMax}];
  redCoeff = fCoeff;
  Do[redCoeff = ReduceOnce[redCoeff, mIndex], {k, 1, kSR}];
  lKeep = Keys[redCoeff];
  fitResidual = N[Norm[yFitMatGlobal.sol - fFitVec]/Max[Norm[fFitVec], 10^-30], 20];

  Table[
    Module[{thetaObs, yAssoc, fVal},
      thetaObs = thetaRecords[[thetaIndex, "theta"]];
      yAssoc = yThetaAssocGlobal[thetaIndex];
      fVal = N[
        Total@Table[redCoeff[L] Lookup[yAssoc, L, 0], {L, lKeep}] / (1 - Cos[thetaObs])^kSR,
        workingPrec
      ];
      <|
        "a" -> N[a, 20],
        "omega" -> N[w, 20],
        "theta_obs" -> N[thetaObs, 20],
        "theta_deg" -> N[thetaRecords[[thetaIndex, "thetaDeg"]], 20],
        "k_sr" -> kSR,
        "l_mode_max" -> lModeMax,
        "l_basis_max" -> lBasisMax,
        "f_sr_re" -> ToString[N[Re[fVal], workingPrec], InputForm],
        "f_sr_im" -> ToString[N[Im[fVal], workingPrec], InputForm],
        "fit_residual" -> ToString[fitResidual, InputForm]
      |>
    ],
    {thetaIndex, Length[thetaRecords]}
  ]
];

Print["[*] Fig789 theta-grid generation"];
Print["    theta deg list: ", thetaDegList];
If[ListQ[aFilter] && Length[aFilter] > 0, Print["    a-filter (KeyA): ", aFilter]];
Print["    combos: ", Length[comboKeys], " / ", Length[comboKeysAll]];
Print["    l_mode_max=", lModeMax, ", l_basis_max=", lBasisMax, ", k=", kSR];

Print["[*] Launching parallel kernels..."];
LaunchKernels[];
ParallelEvaluate[
  AppendTo[$Path, "C:/Users/97747/SpinWeightedSpheroidalHarmonics/Kernel"];
  Needs["SpinWeightedSpheroidalHarmonics`"];
];
DistributeDefinitions[
  sWeight, mIndex, workingPrec, kSR, lModeMax, lBasisMax, nThetaFit,
  thetaRecords, thetaFitGrid, yFitMatGlobal, yThetaAssocGlobal,
  modeAssoc, comboKeys, IsFiniteNumeric, SafeToString, ComputeHarmonics,
  ReduceOnce, ComputeThetaModeRow, ComputeSRRows
];
Print["[*] Active kernels: ", $KernelCount];

Print["[*] Computing theta-specific modal rows..."];
modalRaw = ParallelTable[
  ComputeThetaModeRow[thetaIndex, pair[[1]], pair[[2]], l],
  {thetaIndex, Length[thetaRecords]},
  {pair, comboKeys},
  {l, 2, lModeMax},
  Method -> "FinestGrained"
];
modalRows = Cases[modalRaw, _Association, Infinity];
Print["[*] Modal rows: ", Length[modalRows], "; valid: ", Count[modalRows[[All, "valid"]], True]];

Print["[*] Computing generic spherical-SR f rows..."];
srRaw = ParallelTable[
  ComputeSRRows[pair[[1]], pair[[2]]],
  {pair, comboKeys},
  Method -> "FinestGrained"
];
srRows = Cases[srRaw, _Association, Infinity];
Print["[*] SR rows: ", Length[srRows]];

Export[modalOutPath, modalRows, "JSON"];
Export[srOutPath, srRows, "JSON"];

Print["[+] Export complete: data/fig789_patch_data.json"];
Print["[+] Export complete: data/fig789_sr_spherical.json"];
