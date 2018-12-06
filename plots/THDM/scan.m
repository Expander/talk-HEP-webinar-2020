Get["addons/SplitTHDMTHDMTower/SplitTHDMTHDMTower_librarylink.m"];

invalid;
Mtpole = 173.21;
MZpole = 91.1876;
GFInput = 1.16637*^-5;
mbAtmb = 4.18;
alphaSAtMZ = 0.1181;
alphaEmAtMZ = 1/127.950;

sigmaAlphaS = 0.0006;
sigmaMt = 0.98;

SMParameters = {
    alphaEmMZ -> alphaEmAtMZ, (* SMINPUTS[1] *)
    GF -> GFInput,            (* SMINPUTS[2] *)
    alphaSMZ -> alphaSAtMZ,   (* SMINPUTS[3] *)
    MZ -> MZpole,             (* SMINPUTS[4] *)
    mbmb -> mbAtmb,           (* SMINPUTS[5] *)
    Mt -> Mtpole,             (* SMINPUTS[6] *)
    Mtau -> 1.777,            (* SMINPUTS[7] *)
    Mv3 -> 0,                 (* SMINPUTS[8] *)
    MW -> 80.384,             (* SMINPUTS[9] *)
    Me -> 0.000510998902,     (* SMINPUTS[11] *)
    Mv1 -> 0,                 (* SMINPUTS[12] *)
    Mm -> 0.1056583715,       (* SMINPUTS[13] *)
    Mv2 -> 0,                 (* SMINPUTS[14] *)
    md2GeV -> 0.00475,        (* SMINPUTS[21] *)
    mu2GeV -> 0.0024,         (* SMINPUTS[22] *)
    ms2GeV -> 0.104,          (* SMINPUTS[23] *)
    mcmc -> 1.27,             (* SMINPUTS[24] *)
    CKMTheta12 -> 0,
    CKMTheta13 -> 0,
    CKMTheta23 -> 0,
    CKMDelta -> 0,
    PMNSTheta12 -> 0,
    PMNSTheta13 -> 0,
    PMNSTheta23 -> 0,
    PMNSDelta -> 0,
    PMNSAlpha1 -> 0,
    PMNSAlpha2 -> 0,
    alphaEm0 -> 1/137.035999074,
    Mh -> 125.09
};

LinearRange[start_, stop_, steps_] :=
    Table[start + i/steps (stop - start), {i, 0, steps}];

LogRange[start_, stop_, steps_] :=
    Module[{i, result = {}},
           For[i = 0, i <= steps, i++,
               result = AppendTo[result, Exp[Log[start] + (Log[stop] - Log[start]) i / steps]];
              ];
           result
          ];

RunSplitTHDMTHDMTower[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, MA_?NumericQ,
                      Mu_?NumericQ, M12_?NumericQ, M3_?NumericQ,
                      ytLoops_:2, Qpole_:0, Qin_:0, Qmat_:0, QDR_:0, calcUncerts_:False, eft_:0, yt_:0] :=
    Module[{handle, spectrum, uncerts = {}, Qi = Qin, Qm = Qmat},
           If[Qi === 0, Qi = M12];
           If[Qm === 0, Qm = M12];
           handle = FSSplitTHDMTHDMTowerOpenHandle[
               fsSettings -> {
                   precisionGoal -> 1.*^-5,           (* FlexibleSUSY[0] *)
                   maxIterations -> 100,              (* FlexibleSUSY[1] *)
                   calculateStandardModelMasses -> 0, (* FlexibleSUSY[3] *)
                   poleMassLoopOrder -> 3,            (* FlexibleSUSY[4] *)
                   ewsbLoopOrder -> 3,                (* FlexibleSUSY[5] *)
                   betaFunctionLoopOrder -> 2,        (* FlexibleSUSY[6] *)
                   thresholdCorrectionsLoopOrder -> ytLoops,(* FlexibleSUSY[7] *)
                   higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
                   higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
                   higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
                   higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
                   forceOutput -> 0,                  (* FlexibleSUSY[12] *)
                   topPoleQCDCorrections -> 1,        (* FlexibleSUSY[13] *)
                   betaZeroThreshold -> 1.*^-11,      (* FlexibleSUSY[14] *)
                   forcePositiveMasses -> 0,          (* FlexibleSUSY[16] *)
                   poleMassScale -> Qpole,            (* FlexibleSUSY[17] *)
                   eftPoleMassScale -> 0,             (* FlexibleSUSY[18] *)
                   eftMatchingScale -> 0,             (* FlexibleSUSY[19] *)
                   eftMatchingLoopOrderUp -> 0,       (* FlexibleSUSY[20] *)
                   eftMatchingLoopOrderDown -> 0,     (* FlexibleSUSY[21] *)
                   eftHiggsIndex -> 0,                (* FlexibleSUSY[22] *)
                   calculateBSMMasses -> 1,           (* FlexibleSUSY[23] *)
                   parameterOutputScale -> QDR        (* MODSEL[12] *)
               },
               fsSMParameters -> SMParameters,
               fsModelParameters -> {
                   TanBeta -> TB,
                   MSUSY -> MS,
                   MEWSB -> Mtpole,
                   MuInput -> Mu,
                   M1Input -> M12,
                   M2Input -> M12,
                   M3Input -> M3,
                   MAInput -> MA,
                   Qinput -> Qi,
                   Qmatch -> Qm,
                   LambdaLoopOrder -> 2,
                   DeltaAlphaS -> sigmaAlphaS,
                   DeltaMTopPole -> sigmaMt,
                   DeltaEFT -> eft,
                   DeltaYt -> yt,
                   AeInput -> 0 MS TB IdentityMatrix[3],
                   AdInput -> 0 MS TB IdentityMatrix[3],
                   AuInput -> {
                       {0, 0, 0            },
                       {0, 0, 0            },
                       {0, 0, Mu/TB + Xt MS}
                   },
                   msqInput -> MS {1,1,1},
                   msuInput -> MS {1,1,1},
                   msdInput -> MS {1,1,1},
                   mslInput -> MS {1,1,1},
                   mseInput -> MS {1,1,1}
               }
           ];
           spectrum = FSSplitTHDMTHDMTowerCalculateSpectrum[handle];
           If[spectrum === $Failed,
              FSSplitTHDMTHDMTowerCloseHandle[handle];
              $Failed
              ,
              spectrum = HGTHDMIIMSSMBCFull /. spectrum;
              If[calcUncerts,
                 uncerts = FSSplitTHDMTHDMTowerCalculateUncertainties[handle];
                 If[uncerts === $Failed,
                    FSSplitTHDMTHDMTowerCloseHandle[handle];
                    Return[$Failed];
                   ];
                 uncerts = HGTHDMIIMSSMBCFull /. uncerts;
                ];
              FSSplitTHDMTHDMTowerCloseHandle[handle];
              If[calcUncerts, uncerts, spectrum]
             ]
          ];

GetPar[spec_, par__] :=
    GetPar[spec, #]& /@ {par};

GetPar[spec_, par_] :=
    If[spec =!= $Failed, (par /. spec), invalid];

GetPar[spec_, par_[n__?IntegerQ]] :=
    If[spec =!= $Failed, (par /. spec)[[n]], invalid];

RunSplitTHDMTHDMTowerMh[args__] :=
    Module[{spec = RunSplitTHDMTHDMTower[args]},
           If[spec === $Failed,
              invalid,
              GetPar[spec, Pole[M[hh]][1]]
             ]
          ];

RunSplitTHDMTHDMTowerUncertainties[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, MA_?NumericQ,
                                   Mu_?NumericQ, M12_?NumericQ, M3_?NumericQ,
                                   ytLoops_:2, Qpole_:0, Qin_:0, Qmat_:0, QDR_:0] :=
    Module[{uncerts, MhEFT, MhYt},
           uncerts = RunSplitTHDMTHDMTower[MS, TB, Xt, MA, Mu, M12, M3, ytLoops, Qpole, Qin, Qmat, QDR, True];
           (* with extra terms ~ v^2/MS^2 *)
           MhEFT = GetPar[RunSplitTHDMTHDMTower[MS, TB, Xt, MA, Mu, M12, M3, ytLoops, Qpole, Qin, Qmat, QDR, False, 1, 0], Pole[M[hh]][1]];
           MhYt  = GetPar[RunSplitTHDMTHDMTower[MS, TB, Xt, MA, Mu, M12, M3, ytLoops, Qpole, Qin, Qmat, QDR, False, 0, 1], Pole[M[hh]][1]];
           If[uncerts === $Failed,
              { invalid, invalid, invalid, invalid, invalid, invalid, MhEFT, MhYt },
              {
                  (Pole[M[hh]] /. (MIN /. (SUSYScale /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MAX /. (SUSYScale /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MIN /. (AlphaSInput /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MAX /. (AlphaSInput /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MIN /. (MTopPoleInput /. uncerts)))[[1]],
                  (Pole[M[hh]] /. (MAX /. (MTopPoleInput /. uncerts)))[[1]],
                  MhEFT,
                  MhYt
              }
             ]
          ];

IsValid[ex_]          := FreeQ[ex, invalid];
RemoveInvalid[l_List] := Select[l, IsValid];
MaxDiff[l_List]       := Max[Abs[RemoveInvalid[l]]] - Min[Abs[RemoveInvalid[l]]];
MaxDiff[{}]           := 0;
MinMax[l_List]        := { Min[Abs[RemoveInvalid[l]]], Max[Abs[RemoveInvalid[l]]] };
MinMax[{}]            := { invalid, invalid };

RunTHDMSplitTHDM[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ, MA_?NumericQ,
                 Mu_?NumericQ, M12_?NumericQ, M3_?NumericQ] :=
    Module[{Mhyt1L, Mhyt2L, Mhyt3L, DMh},
           Mhyt1L = RunSplitTHDMTHDMTowerMh[MS, TB, Xt, MA, Mu, M12, M3, 1];
           Mhyt2L = RunSplitTHDMTHDMTowerMh[MS, TB, Xt, MA, Mu, M12, M3, 2];
           Mhyt3L = RunSplitTHDMTHDMTowerMh[MS, TB, Xt, MA, Mu, M12, M3, 3];
           DMh = RunSplitTHDMTHDMTowerUncertainties[MS, TB, Xt, MA, Mu, M12, M3, 2];
           (* Mhyt1L, Mhyt2L, Mhyt3L, min DMh^Qpole, max DMh^Qpole *)
           {Mhyt1L, Mhyt2L, Mhyt3L, Sequence @@ DMh}
          ];

steps = 60;

MSstart = 500;
MSstop  = 1.0 10^16;
Xtstart = -3.5;
Xtstop  = 3.5;

(********** THDM+split -> THDM tower degenerate masses: MA = [100, 500], MS = [1000, 10^16], Xt = ? **********)

ScanSplitTHDMTHDMTowerMSMA[Xt_, TB_, Mui_, MSstart_:1000, MSstop_:1.0 10^5, MAstart_:100, MAstop_:500, steps_:60] :=
    Module[{res, tuples},
           tuples = Tuples[{LinearRange[MAstart, MAstop, steps], LogRange[MSstart, MSstop, steps]}];
           res = {N[#[[2]]], Mui, Mui, Mui, TB, N[#[[1]]], Xt,
                  Sequence @@ RunTHDMSplitTHDM[#[[2]], TB, Xt, #[[1]], Mui, Mui, Mui]}& /@ tuples;
           Export["SplitTHDMTHDMTower_MS_MA_Xt-" <> ToString[Xt] <> "_TB-" <> ToString[TB] <>
                  "_Mu-M12-M3-" <> ToString[Mui] <> ".dat", res, "Table"];
           res
          ];

ParallelMap[ScanSplitTHDMTHDMTowerMSMA[0, #, 2000]&, {2, 10, 20, 50}];
ParallelMap[ScanSplitTHDMTHDMTowerMSMA[N@Sqrt[6], #, 2000]&, {2, 10, 20, 50}];
