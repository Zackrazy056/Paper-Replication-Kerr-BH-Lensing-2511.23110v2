(* =================================================================== *)
(* Kerr GW lensing: densify Fig.3 patch angular grids                  *)
(* Recompute only S(theta) and S(pi-theta), preserving radial data.     *)
(* =================================================================== *)

AppendTo[$Path, "C:/Users/97747/SpinWeightedSpheroidalHarmonics/Kernel"];
Needs["SpinWeightedSpheroidalHarmonics`"];

Print["[*] Initializing Parallel Kernels..."];
LaunchKernels[];
ParallelEvaluate[
  AppendTo[$Path, "C:/Users/97747/SpinWeightedSpheroidalHarmonics/Kernel"];
  Needs["SpinWeightedSpheroidalHarmonics`"];
];
Print["[*] Active kernels: ", $KernelCount];

sWeight = -2;
mIndex = 2;
workingPrec = 80;
thetaList = N[Range[40, 180, 1/2] * Pi / 180, workingPrec];

inputName = "fig3_patch_w2_a0_a099.json";
targetAList = {0.99};
targetOmegaList = {2};

If[StringQ[Environment["FIG3_DENSIFY_FILE"]] && StringLength[StringTrim[Environment["FIG3_DENSIFY_FILE"]]] > 0,
  inputName = StringTrim[Environment["FIG3_DENSIFY_FILE"]];
];
If[StringQ[Environment["FIG3_DENSIFY_A_LIST"]] && StringLength[StringTrim[Environment["FIG3_DENSIFY_A_LIST"]]] > 0,
  targetAList = Quiet[Check[ToExpression /@ StringSplit[StringTrim[Environment["FIG3_DENSIFY_A_LIST"]], ","], targetAList]];
];
If[StringQ[Environment["FIG3_DENSIFY_OMEGA_LIST"]] && StringLength[StringTrim[Environment["FIG3_DENSIFY_OMEGA_LIST"]]] > 0,
  targetOmegaList = Quiet[Check[ToExpression /@ StringSplit[StringTrim[Environment["FIG3_DENSIFY_OMEGA_LIST"]], ","], targetOmegaList]];
];

dataPath = ExpandFileName[
  FileNameJoin[{DirectoryName[$InputFileName], "..", "data", inputName}]
];

If[!FileExistsQ[dataPath],
  Print["[!] Missing input file: ", dataPath];
  Quit[1];
];

KeyA[a_] := Round[N[a, 20], 10^-3];
KeyW[w_] := Round[N[w, 20], 10^-10];
targetKeys = Flatten[Table[{KeyA[a], KeyW[w]}, {a, targetAList}, {w, targetOmegaList}], 1];

IsFiniteNumeric[x_, prec_: workingPrec] := Module[{y},
  y = Quiet[Check[N[x, prec], $Failed]];
  y =!= $Failed && NumberQ[y] && FreeQ[y, ComplexInfinity | DirectedInfinity[_] | Indeterminate | Infinity]
];

SafeToString[x_, prec_: workingPrec] := Module[{y},
  y = Quiet[Check[N[x, prec], $Failed]];
  If[y === $Failed || !IsFiniteNumeric[y, prec], "0", ToString[y, InputForm]]
];

NeedsDensifyQ[row_] := Module[{sTheta, sPi, key},
  If[!AssociationQ[row], Return[False]];
  If[KeyExistsQ[row, "valid"] && !TrueQ[row["valid"]], Return[False]];
  If[!(KeyExistsQ[row, "l"] && KeyExistsQ[row, "a"] && KeyExistsQ[row, "omega"]), Return[False]];
  key = {KeyA[row["a"]], KeyW[row["omega"]]};
  If[!MemberQ[targetKeys, key], Return[False]];
  sTheta = If[KeyExistsQ[row, "S_theta"], row["S_theta"], {}];
  sPi = If[KeyExistsQ[row, "S_pi_minus_theta"], row["S_pi_minus_theta"], {}];
  !(ListQ[sTheta] && ListQ[sPi] && Length[sTheta] == Length[thetaList] && Length[sPi] == Length[thetaList])
];

DensifyRow[row_] := Module[
  {l, a, w, c, ntList, nt, sh, s0, sTheta, sPiMinusTheta},
  l = Round[row["l"]];
  a = N[row["a"], workingPrec];
  w = N[row["omega"], workingPrec];
  c = SetPrecision[a * w, workingPrec];
  ntList = DeleteDuplicates[{Max[80, l + 20], 100, 120, 160, 200}];
  sh = $Failed;
  Do[
    sh = Quiet[
      Check[
        SpinWeightedSpheroidalHarmonicS[
          sWeight, l, mIndex, c,
          Method -> {"SphericalExpansion", "NumTerms" -> nt}
        ],
        $Failed
      ]
    ];
    If[sh =!= $Failed, Break[]],
    {nt, ntList}
  ];
  If[sh === $Failed,
    sh = Quiet[Check[SpinWeightedSpheroidalHarmonicS[sWeight, l, mIndex, c], $Failed]];
  ];
  If[sh === $Failed, Return[Join[row, <|"densify_failed" -> True|>]]];

  s0 = Quiet[Check[N[sh[0, 0], workingPrec], $Failed]];
  sTheta = Quiet[Check[Table[N[sh[th, 0], workingPrec], {th, thetaList}], $Failed]];
  sPiMinusTheta = Quiet[Check[Table[N[sh[Pi - th, 0], workingPrec], {th, thetaList}], $Failed]];
  If[
    s0 === $Failed || sTheta === $Failed || sPiMinusTheta === $Failed ||
      !IsFiniteNumeric[s0] ||
      !AllTrue[sTheta, IsFiniteNumeric] ||
      !AllTrue[sPiMinusTheta, IsFiniteNumeric],
    Return[Join[row, <|"densify_failed" -> True|>]]
  ];

  Join[
    row,
    <|
      "S_0" -> SafeToString[s0],
      "S_theta" -> (SafeToString /@ sTheta),
      "S_pi_minus_theta" -> (SafeToString /@ sPiMinusTheta),
      "valid" -> True
    |>
  ]
];

DistributeDefinitions[
  sWeight, mIndex, workingPrec, thetaList, KeyA, KeyW, targetKeys,
  IsFiniteNumeric, SafeToString, NeedsDensifyQ, DensifyRow
];

rows = Association /@ Import[dataPath, "RawJSON"];
indices = Flatten@Position[NeedsDensifyQ /@ rows, True];

Print["[*] Input: ", dataPath];
Print["[*] Target keys: ", targetKeys];
Print["[*] Rows needing densify: ", Length[indices]];

If[Length[indices] == 0,
  Print["[*] Nothing to densify."];
  Quit[0];
];

newRows = ParallelTable[DensifyRow[rows[[idx]]], {idx, indices}, Method -> "FinestGrained"];

Do[
  rows[[indices[[k]]]] = newRows[[k]],
  {k, 1, Length[indices]}
];

Export[dataPath, rows, "JSON"];

Print["[+] Export complete: ", dataPath];
Print["[*] Densified valid rows: ", Count[newRows[[All, "valid"]], True], "/", Length[newRows]];
