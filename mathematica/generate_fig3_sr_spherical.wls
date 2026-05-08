(* =================================================================== *)
(* Kerr GW lensing: Fig.3 spherical-basis SR data generation           *)
(* Appendix-C aligned workflow:                                        *)
(* modal f(theta) on spheroidal basis -> spherical projection -> SR    *)
(* =================================================================== *)

AppendTo[$Path, "C:/Users/97747/SpinWeightedSpheroidalHarmonics/Kernel"];
Needs["SpinWeightedSpheroidalHarmonics`"];

sWeight = -2;
mIndex = 2;
workingPrec = 70;

kSR = 2;
lModeTarget = 88;
lBasisPad = 16;
minModeCount = 20;
nFullThetaSteps = 400;
spheroidalNumTerms = 40;
defaultComboTargets = {
  {0.0, 0.5},
  {0.0, 0.6},
  {0.0, 2.0},
  {0.99, 0.5},
  {0.99, 0.6},
  {0.99, 2.0}
};

workingPrecRaw = Environment["FIG3_SR_WORKING_PREC"];
If[StringQ[workingPrecRaw] && StringLength[StringTrim[workingPrecRaw]] > 0,
  workingPrec = Quiet[Check[ToExpression[workingPrecRaw], workingPrec]];
];
If[!IntegerQ[workingPrec] || workingPrec < 20, workingPrec = 70];

lModeRaw = Environment["FIG3_SR_LMAX"];
If[StringQ[lModeRaw] && StringLength[StringTrim[lModeRaw]] > 0,
  lModeTarget = Quiet[Check[ToExpression[lModeRaw], lModeTarget]];
];
If[!IntegerQ[lModeTarget] || lModeTarget < 10, lModeTarget = 88];

fullStepsRaw = Environment["FIG3_SR_FULL_STEPS"];
If[StringQ[fullStepsRaw] && StringLength[StringTrim[fullStepsRaw]] > 0,
  nFullThetaSteps = Quiet[Check[ToExpression[fullStepsRaw], nFullThetaSteps]];
];
If[!IntegerQ[nFullThetaSteps] || nFullThetaSteps < 80, nFullThetaSteps = 400];

numTermsRaw = Environment["FIG3_SR_NUM_TERMS"];
If[StringQ[numTermsRaw] && StringLength[StringTrim[numTermsRaw]] > 0,
  spheroidalNumTerms = Quiet[Check[ToExpression[numTermsRaw], spheroidalNumTerms]];
];
If[!IntegerQ[spheroidalNumTerms] || spheroidalNumTerms < 8, spheroidalNumTerms = 40];

thetaDegList = N[Range[40, 180, 1/2], workingPrec];
thetaRadList = N[thetaDegList * Pi / 180, workingPrec];
nTheta = Length[thetaRadList];
fullThetaRadList = N[Subdivide[0, Pi, nFullThetaSteps], workingPrec];
nFullTheta = Length[fullThetaRadList];

repoDir = DirectoryName[$InputFileName] <> "/..";
dataFiles = {
  ExpandFileName[FileNameJoin[{repoDir, "data", "fig3_data.json"}]],
  ExpandFileName[FileNameJoin[{repoDir, "data", "fig3_patch_w2_a0_a099.json"}]]
};

CleanNum[x_] := Module[{s, v},
  s = If[StringQ[x], StringTrim[x], ToString[x, InputForm]];
  s = StringReplace[s, {"\"" -> "", " " -> "", "*^" -> "e"}];
  s = StringReplace[s, RegularExpression["`[-+0-9.]*"] -> ""];
  If[
    StringLength[s] == 0 || StringContainsQ[s, "Indeterminate"] || StringContainsQ[s, "Infinity"],
    Return[0]
  ];
  v = Quiet[Check[ToExpression[s], 0]];
  If[NumericQ[v], N[v, workingPrec], 0]
];

SafeNumList[x_] := Module[{vals},
  If[!ListQ[x], Return[{}]];
  vals = CleanNum /@ x;
  vals
];

FmtKey[val_, ndigits_] := Module[{s},
  s = ToString@NumberForm[N[val, 30], {Infinity, ndigits}, NumberPadding -> {"", "0"}, NumberSeparator -> "", NumberPoint -> ".", ExponentFunction -> (Null &)];
  StringReplace[s, " " -> ""]
];
KeyA[a_] := FmtKey[Round[N[a, 20], 10^-3], 3];
KeyW[w_] := FmtKey[Round[N[w, 20], 10^-4], 4];

rawData = Flatten[
  Table[
    If[FileExistsQ[f], Import[f, "JSON"], {}],
    {f, dataFiles}
  ],
  1
];

rawAssoc = Association /@ rawData;
validRows = Select[
  rawAssoc,
  (
    (!KeyExistsQ[#, "valid"] || TrueQ[#["valid"]]) &&
    KeyExistsQ[#, "l"] &&
    KeyExistsQ[#, "a"] &&
    KeyExistsQ[#, "omega"] &&
    KeyExistsQ[#, "S_0"] &&
    KeyExistsQ[#, "Re_C"] &&
    KeyExistsQ[#, "B_ratio_re"] &&
    KeyExistsQ[#, "B_ratio_im"]
  ) &
];

modeAssoc = <||>;
Do[
  Module[{a, w, l, k, rec},
    a = CleanNum[row["a"]];
    w = CleanNum[row["omega"]];
    l = Round[CleanNum[row["l"]]];
    If[l < 2 || l > lModeTarget, Continue[]];

    k = {KeyA[a], KeyW[w], l};
    rec = <|
      "a" -> a,
      "w" -> w,
      "l" -> l,
      "ReC" -> CleanNum[row["Re_C"]],
      "B" -> (CleanNum[row["B_ratio_re"]] + I * CleanNum[row["B_ratio_im"]]),
      "S0" -> CleanNum[row["S_0"]]
    |>;

    (* Later files override earlier rows for the same (a,w,l). *)
    modeAssoc[k] = rec;
  ],
  {row, validRows}
];

comboKeysRaw = Keys[modeAssoc];
comboKeysRaw = Select[
  comboKeysRaw,
  (ListQ[#] && Length[#] == 3 && StringQ[#[[1]]] && StringQ[#[[2]]] && IntegerQ[Round[#[[3]]]]) &
];
comboKeys = DeleteDuplicates[comboKeysRaw[[All, {1, 2}]]];
comboTargets = ({KeyA[#[[1]]], KeyW[#[[2]]]} &) /@ defaultComboTargets;
comboKeys = Select[
  comboKeys,
  MemberQ[comboTargets, #] &
];

comboRaw = Environment["FIG3_SR_COMBOS"];
If[StringQ[comboRaw] && StringLength[StringTrim[comboRaw]] > 0,
  Module[{pairs, parsed},
    pairs = StringSplit[StringTrim[comboRaw], ";"];
    parsed = DeleteCases[
      Quiet@Check[
        Module[{parts, aVal, wVal},
          parts = StringSplit[#, ","];
          If[Length[parts] =!= 2, Nothing,
            aVal = ToExpression[StringTrim[parts[[1]]]];
            wVal = ToExpression[StringTrim[parts[[2]]]];
            {KeyA[aVal], KeyW[wVal]}
          ]
        ] & /@ pairs,
        {}
      ],
      Nothing
    ];
    If[Length[parsed] > 0, comboKeys = Select[comboKeys, MemberQ[parsed, #] &]];
  ];
];

batchRaw = Environment["FIG3_SR_BATCH"];
maxCombosPerRun = Length[comboKeys];
If[StringQ[batchRaw] && StringLength[StringTrim[batchRaw]] > 0,
  maxCombosPerRun = Quiet[Check[ToExpression[batchRaw], maxCombosPerRun]];
];
If[!IntegerQ[maxCombosPerRun] || maxCombosPerRun <= 0, maxCombosPerRun = Length[comboKeys]];
comboKeys = Take[comboKeys, UpTo[maxCombosPerRun]];

ReduceOnce[coeff_Association, m_:2] := Module[
  {lVals, lMin, lMax, out, fl, fp, fm, termDiag, aNum, aDen, bNum, bDen, aL, bL},
  lVals = Sort[Keys[coeff]];
  If[Length[lVals] == 0, Return[<||>]];
  lMin = Min[lVals];
  lMax = Max[lVals];
  out = <||>;
  Do[
    fl = Lookup[coeff, l, 0];
    fp = Lookup[coeff, l + 1, 0];
    fm = Lookup[coeff, l - 1, 0];

    termDiag = 1 - (2 m)/(l (l + 1));
    aNum = (l - 1) (l + 3) (l - m + 1) (l + m + 1);
    aDen = 4 l (l + 2) + 3;
    aL = If[aNum >= 0, N[Sqrt[aNum/aDen], workingPrec], 0];

    bNum = (l^2 - m^2) (l^2 - 4);
    bDen = 4 l^2 - 1;
    bL = If[bNum >= 0, N[Sqrt[bNum/bDen], workingPrec], 0];

    out[l] = termDiag fl - (aL/(l + 1)) fp - (bL/l) fm;
    ,
    {l, lMin, lMax}
  ];
  out
];

ComputeSRRow[ka_, kw_] := Module[
  {
    modeRows, a, w, pref,
    coeffA, lKeys,
    lModeMax, lBasisMax,
    cVal, fullSThetaAssoc,
    fullYBasisAssoc, fullFFitVec,
    projWeightVec, trapFactorVec, integWeightVec,
    fullWVec, fullFFitRecon, resVecW, dataVecW, fitCondition,
    yOutBasisAssoc,
    fCoeff, redCoeff, redKeys,
    fSrVec, fitResidual
  },
  modeRows = Lookup[
    modeAssoc,
    Select[
      Keys[modeAssoc],
      (ListQ[#] && Length[#] == 3 && #[[1]] === ka && #[[2]] === kw) &
    ],
    {}
  ];
  modeRows = Select[
    modeRows,
    (AssociationQ[#] && KeyExistsQ[#, "l"] && KeyExistsQ[#, "a"] && KeyExistsQ[#, "w"] && KeyExistsQ[#, "S0"]) &
  ];
  If[Length[modeRows] < minModeCount,
    Print["    [skip] aKey=", ka, ", wKey=", kw, ", modes=", Length[modeRows]];
    Return[Nothing]
  ];

  a = modeRows[[1, "a"]];
  w = modeRows[[1, "w"]];
  If[w <= 0, Return[Nothing]];

  pref = N[(2 Pi)/(I w), workingPrec];

  coeffA = Association@Table[
    Module[{l = row["l"], bracket},
      bracket = (-1)^(l + 1) (row["ReC"]/(16 w^4)) row["B"] - 1;
      l -> N[pref row["S0"] bracket, workingPrec]
    ],
    {row, modeRows}
  ];

  lKeys = Sort[Keys[coeffA]];
  lModeMax = Max[lKeys];
  lBasisMax = lModeMax + lBasisPad;

  cVal = SetPrecision[a * w, workingPrec];
  fullSThetaAssoc = Association@Table[
    Module[{sFun, vals},
      sFun = Quiet[
        Check[
          SpinWeightedSpheroidalHarmonicS[
            sWeight, l, mIndex, cVal,
            Method -> {"SphericalExpansion", "NumTerms" -> spheroidalNumTerms}
          ],
          $Failed
        ]
      ];
      If[sFun === $Failed,
        l -> {},
        vals = Quiet[Check[(N[sFun[#, 0], workingPrec] &) /@ fullThetaRadList, {}]];
        l -> vals
      ]
    ],
    {l, lKeys}
  ];
  If[AnyTrue[Values[fullSThetaAssoc], (Not[ListQ[#]] || Length[#] =!= nFullTheta) &],
    Print["    [skip] failed spheroidal full-grid build for a=", a, ", w=", w];
    Return[Nothing]
  ];

  fullYBasisAssoc = Association@Table[
    L -> Table[
      N[SpinWeightedSphericalHarmonicY[sWeight, L, mIndex, th, 0], workingPrec],
      {th, fullThetaRadList}
    ],
    {L, 2, lBasisMax}
  ];

  fullFFitVec = Table[
    N[Total@Table[coeffA[l] * fullSThetaAssoc[l][[i]], {l, lKeys}], workingPrec],
    {i, 1, nFullTheta}
  ];

  (* Appendix-C aligned spherical projection via discrete inner product:
     F_Lm = 2 Pi \int f(theta,0) Conjugate[Y_Lm(theta,0)] Sin[theta] dtheta *)
  projWeightVec = N[(Sin[#] &) /@ fullThetaRadList, workingPrec];
  trapFactorVec = Table[If[i == 1 || i == nFullTheta, 1/2, 1], {i, 1, nFullTheta}];
  integWeightVec = N[
    MapThread[#1 * #2 &, {projWeightVec, trapFactorVec}] * (2 Pi * Pi / nFullThetaSteps),
    workingPrec
  ];

  fCoeff = Association@Table[
    L -> N[
      Total@MapThread[
        #1 * Conjugate[#2] * #3 &,
        {fullFFitVec, fullYBasisAssoc[L], integWeightVec}
      ],
      workingPrec
    ],
    {L, 2, lBasisMax}
  ];
  fitCondition = "NA_orth_projection";

  redCoeff = fCoeff;
  Do[redCoeff = ReduceOnce[redCoeff, mIndex], {kk, 1, kSR}];
  redKeys = Sort[Keys[redCoeff]];

  yOutBasisAssoc = Association@Table[
    L -> Table[
      N[SpinWeightedSphericalHarmonicY[sWeight, L, mIndex, th, 0], workingPrec],
      {th, thetaRadList}
    ],
    {L, 2, lBasisMax}
  ];

  fSrVec = Table[
    N[
      Total@Table[
        redCoeff[L] * yOutBasisAssoc[L][[i]],
        {L, redKeys}
      ] / (1 - Cos[thetaRadList[[i]]])^kSR,
      workingPrec
    ],
    {i, 1, nTheta}
  ];

  fullFFitRecon = Table[
    N[
      Total@Table[
        fCoeff[L] * fullYBasisAssoc[L][[i]],
        {L, 2, lBasisMax}
      ],
      workingPrec
    ],
    {i, 1, nFullTheta}
  ];
  fullWVec = N[(Sqrt[Max[Sin[#], 10^-30]] &) /@ fullThetaRadList, workingPrec];
  resVecW = MapThread[#1 * #2 &, {fullFFitRecon - fullFFitVec, fullWVec}];
  dataVecW = MapThread[#1 * #2 &, {fullFFitVec, fullWVec}];
  fitResidual = N[Norm[resVecW]/Max[Norm[dataVecW], 10^-30], 20];

  <|
    "a" -> N[a, 20],
    "omega" -> N[w, 20],
    "k_sr" -> kSR,
    "l_mode_max" -> lModeMax,
    "l_basis_max" -> lBasisMax,
    "fit_residual" -> ToString[fitResidual, InputForm],
    "fit_condition" -> ToString[fitCondition, InputForm],
    "theta_deg_list" -> (ToString[#, InputForm] & /@ thetaDegList),
    "f_sr_re" -> (ToString[N[Re[#], workingPrec], InputForm] & /@ fSrVec),
    "f_sr_im" -> (ToString[N[Im[#], workingPrec], InputForm] & /@ fSrVec)
  |>
];

Print["[*] Building Fig.3 spherical-basis SR rows..."];
Print["    combos: ", Length[comboKeys]];
Print["    l_mode_target=", lModeTarget, ", k=", kSR, ", nTheta=", nTheta];
Print["    full-grid points=", nFullTheta, ", spheroidal NumTerms=", spheroidalNumTerms];
Print["    working precision=", workingPrec];
Print["    projection method=orth_projection"];

srRows = DeleteCases[
  Table[
    ComputeSRRow[pair[[1]], pair[[2]]],
    {pair, comboKeys}
  ],
  Nothing
];

Print["[*] SR rows generated: ", Length[srRows]];

outPath = ExpandFileName[FileNameJoin[{repoDir, "data", "fig3_sr_spherical.json"}]];
Export[outPath, srRows, "JSON"];
Print["[+] Export complete: data/fig3_sr_spherical.json"];
