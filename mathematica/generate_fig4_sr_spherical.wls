(* =================================================================== *)
(* Fig.4 spherical-basis SR(k=2) data generation                       *)
(* Appendix-C aligned workflow: spheroidal modal sum -> spherical fit  *)
(* -> series reduction -> evaluate at theta_o = pi/6                   *)
(* =================================================================== *)

AppendTo[$Path, "C:/Users/97747/SpinWeightedSpheroidalHarmonics/Kernel"];
Needs["SpinWeightedSpheroidalHarmonics`"];

sWeight = -2;
mIndex = 2;
workingPrec = 60;

thetaObs = N[Pi/6, workingPrec];
kSR = 2;
lBasisPad = 16;
nThetaFit = 80;
minModeCount = 20;

lModeMax = 80;
lModeRaw = Environment["FIG4_SR_LMAX"];
If[StringQ[lModeRaw] && StringLength[StringTrim[lModeRaw]] > 0,
  lModeMax = Quiet[Check[ToExpression[lModeRaw], 80]];
];
If[!IntegerQ[lModeMax] || lModeMax < 10, lModeMax = 80];

repoDir = DirectoryName[$InputFileName] <> "/..";
dataFiles = {
  ExpandFileName[FileNameJoin[{repoDir, "data", "fig4_a0_data.json"}]],
  ExpandFileName[FileNameJoin[{repoDir, "data", "fig4_other_a_data.json"}]]
};

defaultOut = "fig4_sr_spherical_l" <> ToString[lModeMax] <> ".json";
outName = Environment["FIG4_SR_OUT"];
If[!(StringQ[outName] && StringLength[StringTrim[outName]] > 0), outName = defaultOut];
outPath = ExpandFileName[FileNameJoin[{repoDir, "data", outName}]];

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
    modeAssoc[k] = rec;
  ],
  {row, validData}
];

allComboKeys = Sort[DeleteDuplicates[(Keys[modeAssoc]) /. {ka_, kw_, _} :> {ka, kw}]];

aFilter = {};
aFilterRaw = Environment["FIG4_SR_A_FILTER"];
If[StringQ[aFilterRaw] && StringLength[StringTrim[aFilterRaw]] > 0,
  aFilter = Quiet[Check[KeyA /@ (ToExpression /@ StringSplit[StringTrim[aFilterRaw], ","]), {}]];
];
If[ListQ[aFilter] && Length[aFilter] > 0,
  allComboKeys = Select[allComboKeys, MemberQ[aFilter, KeyA[#[[1]]]] &];
];

existingRows = {};
If[FileExistsQ[outPath],
  existingRows = Quiet[Check[Association /@ Import[outPath, "RawJSON"], {}]];
];
existingKeys = DeleteDuplicates[
  {KeyA[CleanNum[#["a"]]], KeyW[CleanNum[#["omega"]]]} & /@ existingRows
];
comboKeys = Select[allComboKeys, !MemberQ[existingKeys, #] &];

maxCombosPerRun = Length[comboKeys];
batchEnvRaw = Environment["FIG4_SR_BATCH"];
If[StringQ[batchEnvRaw] && StringLength[StringTrim[batchEnvRaw]] > 0,
  maxCombosPerRun = Quiet[Check[ToExpression[batchEnvRaw], Length[comboKeys]]];
];
If[!IntegerQ[maxCombosPerRun] || maxCombosPerRun <= 0, maxCombosPerRun = Length[comboKeys]];
comboBatch = Take[comboKeys, UpTo[maxCombosPerRun]];

thetaFitGrid = N[
  Subdivide[0.03, Pi - 0.03, nThetaFit - 1],
  workingPrec
];

lBasisMax = lModeMax + lBasisPad;
yFitMatGlobal = Table[
  N[SpinWeightedSphericalHarmonicY[sWeight, L, mIndex, th, 0], workingPrec],
  {th, thetaFitGrid},
  {L, 2, lBasisMax}
];
yObsAssoc = Association@Table[
  L -> N[SpinWeightedSphericalHarmonicY[sWeight, L, mIndex, thetaObs, 0], workingPrec],
  {L, 2, lBasisMax}
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

ComputeSRRow[ka_, kw_] := Module[
  {a, w, pref, modeRows, coeffA, shAssoc, fFitVec, sol, fCoeff, redCoeff, fObs, fitResidual, lKeep},
  modeRows = Table[
    If[KeyExistsQ[modeAssoc, {ka, kw, l}], modeAssoc[{ka, kw, l}], Missing["NotFound"]],
    {l, 2, lModeMax}
  ];
  modeRows = Select[modeRows, # =!= Missing["NotFound"] &];
  If[Length[modeRows] < minModeCount, Return[Nothing]];

  a = modeRows[[1, "a"]];
  w = modeRows[[1, "w"]];
  If[w <= 0, Return[Nothing]];
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
  If[sol === $Failed, Return[Nothing]];

  fCoeff = Association@Table[(L -> N[sol[[L - 1]], workingPrec]), {L, 2, lBasisMax}];
  redCoeff = fCoeff;
  Do[redCoeff = ReduceOnce[redCoeff, mIndex], {kk, 1, kSR}];
  lKeep = Keys[redCoeff];

  fObs = N[
    Total@Table[redCoeff[L] Lookup[yObsAssoc, L, 0], {L, lKeep}] / (1 - Cos[thetaObs])^kSR,
    workingPrec
  ];
  fitResidual = N[
    Norm[yFitMatGlobal.sol - fFitVec]/Max[Norm[fFitVec], 10^-30],
    20
  ];

  <|
    "a" -> N[a, 20],
    "omega" -> N[w, 20],
    "theta_obs" -> N[thetaObs, 20],
    "k_sr" -> kSR,
    "l_mode_max" -> lModeMax,
    "l_basis_max" -> lBasisMax,
    "f_sr_theta_re" -> ToString[N[Re[fObs], workingPrec], InputForm],
    "f_sr_theta_im" -> ToString[N[Im[fObs], workingPrec], InputForm],
    "fit_residual" -> ToString[fitResidual, InputForm]
  |>
];

Print["[*] Building Fig.4 spherical-basis SR rows..."];
Print["    out=", outPath];
Print["    l_mode_max=", lModeMax, ", l_basis_max=", lBasisMax, ", k=", kSR];
If[ListQ[aFilter] && Length[aFilter] > 0, Print["    a-filter (KeyA): ", aFilter]];
Print["    total combos: ", Length[allComboKeys]];
Print["    already done: ", Length[existingKeys]];
Print["    remaining combos: ", Length[comboKeys]];
Print["    batch combos this run: ", Length[comboBatch], " (max=", maxCombosPerRun, ")"];

If[Length[comboBatch] == 0,
  Print["[*] Nothing to do. SR file is already complete."];
  Quit[0];
];

Print["[*] Launching parallel kernels..."];
LaunchKernels[];
ParallelEvaluate[
  AppendTo[$Path, "C:/Users/97747/SpinWeightedSpheroidalHarmonics/Kernel"];
  Needs["SpinWeightedSpheroidalHarmonics`"];
];
DistributeDefinitions[
  sWeight, mIndex, workingPrec, thetaObs, kSR, lModeMax, lBasisMax, nThetaFit,
  modeAssoc, comboBatch, thetaFitGrid, yFitMatGlobal, yObsAssoc,
  ReduceOnce, ComputeSRRow
];
Print["[*] Active kernels: ", $KernelCount];

srRows = DeleteCases[
  ParallelTable[
    ComputeSRRow[pair[[1]], pair[[2]]],
    {pair, comboBatch},
    Method -> "FinestGrained"
  ],
  Nothing
];
Print["[*] SR rows generated: ", Length[srRows]];
If[Length[srRows] == 0,
  Print["[!] No SR rows generated. Aborting export."];
  Quit[1];
];

mergedAssoc = Association@Table[
  Module[{rk},
    rk = {KeyA[CleanNum[row["a"]]], KeyW[CleanNum[row["omega"]]]};
    rk -> row
  ],
  {row, Join[existingRows, srRows]}
];
allRows = SortBy[
  Values[mergedAssoc],
  {KeyA[CleanNum[#["a"]]] &, KeyW[CleanNum[#["omega"]]] &}
];
Export[outPath, allRows, "JSON"];

finalKeys = DeleteDuplicates[
  {KeyA[CleanNum[#["a"]]], KeyW[CleanNum[#["omega"]]]} & /@ allRows
];
remainingAfter = Length[allComboKeys] - Length[finalKeys];

Print["[+] Export complete: ", outPath];
Print["[*] Accumulated SR rows: ", Length[allRows], "; remaining combos after this run: ", remainingAfter];
