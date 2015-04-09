:Evaluate:  Print["MR - Matching and Running Mathmatica interface"];
:Evaluate:  Print["To see available functions use Names[\"mr`*\"]"];
:Evaluate:  Print["Andrey Pikelner <pikelner@theor.jinr.ru>"];

:Evaluate:  BeginPackage["mr`"]

:Evaluate:  RunQCDnf6::usage  = "RunQCDnf6[oscale,asMZ, MZscale, L, Mtpole, mtth] evaluate alphas at specified oscale in n=6 QCD given asMZ at MZscale (nf=5) with L-loop RGE with top (with Mtpole) threshold evaluated at mtth, if last argument is not specified it is assumed that mtth = Mtpole"
:Evaluate:  RunQCD::usage  = "RunQCD[oscale,as0,inscale,L,nf] run as from the inscale to oscale in nf flavour QCD at L-loops given as(inscale) = as0"
:Evaluate:  RunSM::usage  = "RunSM[gp,g,gs,yt,lam,m,iscale,oscale] return running parameters at specified oscale given the values at specified iscale;
			     RunSM[pars_,oscale] returns a list {g1 -> ..., g2 -> ..., ..., scale -> oscale} of running parameters given a list pars";
:Evaluate:  RunSMcouplings::usage  = "RunSMcouplings[gp,g,gs,yt,lam,iscale,oscale,L] return running dimensionless couplings at specified oscale given the values at specified iscale;
			     RunSMcouplings[pars_,oscale,L] returns a list {g1 -> ..., g2 -> ..., ..., scale -> oscale} of running parameters given as a  list pars, L-loop RGES are used";


:Evaluate:  MW::usage  = "MW[gp,g,gs,yb,yt,lam,m,scale] returns pole W-boson mass MW given  MSbar parameters at specified scale at 2-loop level"
:Evaluate:  MZ::usage  = "MZ[gp,g,gs,yb,yt,lam,m,scale] returns pole Z-boson mass MZ given  MSbar parameters at specified scale at 2-loop level"
:Evaluate:  MH::usage  = "MH[gp,g,gs,yb,yt,lam,m,scale] returns pole H-boson mass MH given  MSbar parameters at specified scale at 2-loop level"
:Evaluate:  MT::usage  = "MT[gp,g,gs,yb,yt,lam,m,scale] returns pole t-quark mass MT given  MSbar parameters at specified scale at 2-loop level + 3-loop QCD"
:Evaluate:  GF::usage  = "GF[gp,g,gs,yb,yt,lam,m,scale] returns Fermi constant GF given  MSbar parameters at specified scale at 2-loop level"


:Evaluate:  XMMW::usage  = "XMMW[gp,g,gs,yb,yt,lam,m,scale] returns contributions to the pole W-boson mass MW^2 given  MSbar parameters at specified scale"
:Evaluate:  XMMZ::usage  = "XMMZ[gp,g,gs,yb,yt,lam,m,scale] returns contributions to the pole Z-boson mass MZ^2 given  MSbar parameters at specified scale"
:Evaluate:  XMMH::usage  = "XMMH[gp,g,gs,yb,yt,lam,m,scale] returns contributions to the pole H-boson mass MH^2 given  MSbar parameters at specified scale"
:Evaluate:  XdRbar::usage  = "XdRbar[gp,g,gs,yb,yt,lam,m,scale] returns contributions to the Fermi constant GF - running vev relation given  MSbar parameters at specified scale"

:Evaluate:  XMT::usage  = "XMT[gp,g,gs,yb,yt,lam,m,scale] returns electroweak contributions to the pole top quark mass MT given  MSbar parameters at specified scale"
:Evaluate:  XMTQCD::usage  = "XMT[gp,g,gs,yb,yt,lam,m,scale] returns pure QCD contributions to the pole top quark mass MT given  MSbar parameters at specified scale"


:Evaluate:  Xb::usage  = "Xb[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "
:Evaluate:  XW::usage  = "XW[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "
:Evaluate:  XZ::usage  = "XZ[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "
:Evaluate:  XH::usage  = "XH[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "
:Evaluate:  Xt::usage  = "Xt[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "

:Evaluate:  dROS::usage  = "dROS[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "
:Evaluate:  XbQCD::usage  = "XbQCD[Mb,MW,MZ,MH,Mt,scale,nf=5]  Pure QCD corrections, nf is a number of light quarks, default is 5." 
:Evaluate:  XtQCD::usage  = "XtQCD[Mb,MW,MZ,MH,Mt,scale,nf=5]  Pure QCD corrections, nf is a number of light quarks, default is 5." 
:Evaluate:  mmWMMW::usage  = "mmWMMW[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  mW^2/MW^2, full correction to relation between MS-bar mass mW and pole MW"
:Evaluate:  mmZMMZ::usage  = "mmZMMZ[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  mZ^2/MZ^2, full correction to relation between MS-bar mass mZ and pole MZ"
:Evaluate:  mmHMMH::usage  = "mmHMMH[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  mH^2/MH^2, full correction to relation between MS-bar mass mH and pole MH"
:Evaluate:  mtMt::usage  = "mtMt[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  mt/Mt, full correction to relation between MS-bar mass mt and pole Mt"
:Evaluate:  mbMb::usage  = "mbMb[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  mb/Mb, full correction to relation between MS-bar mass mb and pole Mb"


:Evaluate:  a1::usage  = "a1[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  returns a1=5/3(g1/(4\[Pi]))^2 as a series in running aEW[scale] = e^2/(4 Pi)^2 and aQCD[scale] = gs^2/(4Pi)^2, we use GUT normalization"
:Evaluate:  a2::usage  = "a2[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  returns a2=(g2/(4\[Pi]))^2 as a series in running aEW[scale] = e^2/(4 Pi)^2 and aQCD[scale] = gs^2/(4Pi)^2"
:Evaluate:  at::usage  = "at[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  returns at=(yt/(4\[Pi]))^2 as a series in running aEW[scale] = e^2/(4 Pi)^2 and aQCD[scale] = gs^2/(4Pi)^2"
:Evaluate:  ab::usage  = "ab[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  returns ab=(yb/(4\[Pi]))^2 as a series in running aEW[scale] = e^2/(4 Pi)^2 and aQCD[scale] = gs^2/(4Pi)^2"
:Evaluate:  alam::usage  = "alam[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1] returns alam=\[Lambda]/(4\[Pi])^2 as a series in running aEW[scale] = e^2/(4 Pi)^2 and aQCD[scale] = gs^2/(4Pi)^2"
:Evaluate:  vev::usage  = "vev[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1]  returns running vev[scale] as a series in running aEW[scale] = e^2/(4 Pi)^2 and aQCD[scale] = gs^2/(4Pi)^2"

:Evaluate:  dalphaGF::usage  = "dalphaGF[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1] - calculates loop contributions to the relation between GF and running alpha[scale]"
:Evaluate:  alphaGF::usage  = "alphaGF[Mb,MW,MZ,MH,Mt,scale,nL=2,nH=1] - calculates running alpha[scale] in terms GF and running alpha_s[scale]" 

:Evaluate:  aEW::usage  = "\[Alpha]/(4\[Pi]), running EW constant"
:Evaluate:  aQCD::usage  = "\[Alpha]_S/(4\[Pi]), running QCD constant"

:Evaluate:  Gf::usage  = "G_F Fermi constant"
:Evaluate:  h::usage  = "h - loop counter"

:Evaluate:  xW::usage  = "xW[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between running mW and the pole mass MW"
:Evaluate:  xZ::usage  = "xZ[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between running mZ and the pole mass MZ"
:Evaluate:  xH::usage  = "xH[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between running mH and the pole mass MH"
:Evaluate:  xt::usage  = "xt[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between running mt and the pole mass Mt"
:Evaluate:  xb::usage  = "xb[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between running mb and the pole mass Mb"

:Evaluate:  yW::usage  = "yW[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the running Higgs coupling to W-boson and the pole mass MW"
:Evaluate:  yZ::usage  = "yZ[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the running Higgs coupling to Z-boson and the pole mass MZ"
:Evaluate:  yH::usage  = "yH[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the running Higgs coupling to H-boson and the pole mass MH"
:Evaluate:  yT::usage  = "yT[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the running Higgs coupling to t-quark and the pole mass Mt"
:Evaluate:  yB::usage  = "yB[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the running Higgs coupling to b-quark and the pole mass Mb"


:Evaluate:  xMMW::usage  = "xMW[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the pole mass MW^2 and running parameters"
:Evaluate:  xMMZ::usage  = "xMZ[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the pole mass MZ^2 and running parameters"
:Evaluate:  xMMH::usage  = "xMH[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the pole mass MH^2 and running parameters" 
:Evaluate:  xMT::usage  = "xMT[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the pole mass MT and running parameters"
:Evaluate:  xMTQCD::usage  = "xMTQCD[0,b] represents a coefficient of pure QCD contribution aQCD^b in the relation between the pole mass MT and running parameters"
:Evaluate:  xdRbar::usage  = "xdRbar[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between GF and running parameters (vev, etc)"

:Evaluate:   RunQCDnf6[oscale_?NumericQ,asMZ_?NumericQ,MZscale_?NumericQ,nL_Integer,mtpole_?NumericQ] := RunQCDnf6[oscale, asMZ, MZscale,nL,mtpole,mtpole] /; 1<=nL<=4

:Evaluate:  dr::usage  = "dr[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the running Higgs vev and the Fermi constant GF"

:Evaluate:  daGF::usage  = "daGF[a,b] represents a coefficient of aEW^a * aQCD^b in the relation between the running electromagnetic alpha and Fermi constant GF. Note that aEW should be again expressed in terms of Fermi constant"

:Evaluate:   g1::usage  = "running U(1) coupling"
:Evaluate:   g2::usage  = "running SU(2) coupling"
:Evaluate:   gs::usage  = "running SU(3) strong coupling"
:Evaluate:   yt::usage  = "running top Yukawa coupling"
:Evaluate:   yb::usage  = "running bottom Yukawa coupling"
:Evaluate:   lam::usage  = "running higgs self-coupling"
:Evaluate:   m::usage    = "running higgs mass parameter"
:Evaluate:   vev::usage    = "running higgs vev parameter"
:Evaluate:   scale::usage  = "renormalization scale"
:Evaluate:    QCD::usage  = "Can be used together with \"IncludedCorrections\" option, QCD[l] represent pure l-loop QCD corrections, propotional to aQCD^l"
:Evaluate:     EW::usage  = "Can be used together with \"IncludedCorrections\" option, EW[l] represent pure l-loop EW corrections, propotional to aEW^l"
:Evaluate:  MIXED::usage  = "Can be used together with \"IncludedCorrections\" option, EW[x,y] represent (x+y)-loop QCD corrections, propotional to aEW^x*aQCD^y"

:Evaluate:  Protect[g1,g2,gs,yt,yb,lam,m,scale];

:Evaluate:  Begin["`Private`"]


// Mathematica part
:Evaluate:  EW[x_] = {x,0}; QCD[x_] = {0,x}; MIXED[x_,y_]={x,y};

// Masses
:Evaluate:  Options[mmWMMW] = { "IncludedCorrections" -> { 
							  EW[1],
							  EW[2],
							 MIXED[1,1]
							}
			      };

:Evaluate:  mmWMMW[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1, OptionsPattern[]] := 
		Block[ {lal,las,cr = Map[(#->xW[Sequence @@ # ]) &, OptionValue["IncludedCorrections"]]},
				(1+FromCoefficientRules[cr, {lal,las}]) /.XW[mb,mW,mZ,mH,mt,scale,nL,nH] /. {lal -> aEW[scale], las->aQCD[scale]} 
				]

:Evaluate:  Options[mmZMMZ] = { "IncludedCorrections" -> { 
							  EW[1],
							  EW[2],
							 MIXED[1,1]
							}
			      };

:Evaluate:  mmZMMZ[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1, OptionsPattern[]] := 
		Block[ {lal,las,cr = Map[(#->xZ[Sequence @@ # ]) &, OptionValue["IncludedCorrections"]]},
				(1+FromCoefficientRules[cr, {lal,las}]) /.XZ[mb,mW,mZ,mH,mt,scale,nL,nH] /. {lal -> aEW[scale], las->aQCD[scale]} 
				]

:Evaluate:  Options[mmHMMH] = { "IncludedCorrections" -> { 
							  EW[1],
							  EW[2],
							 MIXED[1,1]
							}
				};

:Evaluate:  mmHMMH[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1, OptionsPattern[]] := 
		Block[ {lal,las,cr = Map[(#->xH[Sequence @@ # ]) &, OptionValue["IncludedCorrections"]]},
				(1+FromCoefficientRules[cr, {lal,las}]) /.XH[mb,mW,mZ,mH,mt,scale,nL,nH] /. {lal -> aEW[scale], las->aQCD[scale]} 
				]
                                                                                                                             
:Evaluate:  Options[mtMt] = { "IncludedCorrections" -> { 
							  EW[1],
							  EW[2],
							 MIXED[1,1]
							}
				};

:Evaluate:  mtMt[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1, OptionsPattern[]] := 
		Block[ {lal,las,cr = Map[(#->xt[Sequence @@ # ]) &, OptionValue["IncludedCorrections"]]},
				(1+FromCoefficientRules[cr, {lal,las}]) /.Xt[mb,mW,mZ,mH,mt,scale,nL,nH] /. XtQCD[mb,mW,mZ,mH,mt,scale,nL,nH] /. {lal -> aEW[scale], las->aQCD[scale]} 
				]

:Evaluate:  Options[mbMb] = { "IncludedCorrections" -> { 
							  EW[1],
							  EW[2],
							 MIXED[1,1]
							}
				};
:Evaluate:  mbMb[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1, OptionsPattern[]] := 
		Block[ {lal,las,cr = Map[(#->xb[Sequence @@ # ]) &, OptionValue["IncludedCorrections"]]},
				(1+FromCoefficientRules[cr, {lal,las}]) /.Xb[mb,mW,mZ,mH,mt,scale,nL,nH] /. XbQCD[mb,mW,mZ,mH,mt,scale,nL,nH] /. {lal -> aEW[scale], las->aQCD[scale]} 
				]
// Couplings

:Evaluate:  Options[a1] = { "IncludedCorrections" -> { 
						          EW[1],  
						         EW[2],  
						     MIXED[1,1]
						     }};

:Evaluate:  a1[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1,OptionsPattern[]] := Block[
				{lal,las,crulesW = Map[(#->yW[Sequence @@ # ]) &, OptionValue["IncludedCorrections"]],
				 crulesZ = Map[(#->yZ[Sequence @@ # ]) &, OptionValue["IncludedCorrections"]]},	
				5/3*2^(5/2)*Gf/(4*Pi)^2*((mZ^2*(1+FromCoefficientRules[crulesZ, {lal,las}]) /.XZ[mb,mW,mZ,mH,mt,scale,nL,nH]) 
						       - (mW^2*(1+FromCoefficientRules[crulesW, {lal,las}]) /.XW[mb,mW,mZ,mH,mt,scale,nL,nH])) /. {lal -> aEW[scale], las->aQCD[scale]} 
				];

:Evaluate:  Options[a2] = { "IncludedCorrections" -> { 
						       EW[1],  
						       EW[2],
						       MIXED[1,1]
						     }};

:Evaluate:  a2[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1,OptionsPattern[]] := Block[
				{lal,las,crulesW = Map[(#->yW[Sequence @@ # ]) &, OptionValue["IncludedCorrections"]] },	
				2^(5/2)*Gf/(4*Pi)^2*((mW^2*(1+FromCoefficientRules[crulesW, {lal,las}]) /.XW[mb,mW,mZ,mH,mt,scale,nL,nH])) /. {lal -> aEW[scale], las->aQCD[scale]}  
				];

:Evaluate:  Options[at] = { "IncludedCorrections" -> { 
						       QCD[1],  
						       QCD[2],  
						       QCD[3],  
						       EW[1],  
						       EW[2],
						       MIXED[1,1]
						     }};

:Evaluate:  at[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1,OptionsPattern[]] := Block[ 
				{lal,las,crulest = Map[(#->yT[Sequence @@ # ]) &, OptionValue["IncludedCorrections"]] },
				2^(3/2)*Gf/(4*Pi)^2*mt^2*((1+FromCoefficientRules[crulest,{lal,las}])/.Xt[mb,mW,mZ,mH,mt,scale,nL,nH] /. XtQCD[mb,mW,mZ,mH,mt,scale,nL,nH])^2 /. {lal -> aEW[scale], las->aQCD[scale]} 
				];

:Evaluate:  Options[ab] = { "IncludedCorrections" -> { 
						       QCD[1],  
						       QCD[2],  
						       QCD[3],  
						       EW[1],  
						       EW[2],
						       MIXED[1,1]
						     }};

:Evaluate:  ab[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1,OptionsPattern[]] := Block[ 
				{lal,las,crulesb = Map[(#->yB[Sequence @@ # ]) &, OptionValue["IncludedCorrections"]] },
				2^(3/2)*Gf/(4*Pi)^2*mb^2*((1+FromCoefficientRules[crulesb,{lal,las}])/.Xb[mb,mW,mZ,mH,mt,scale,nL,nH] /. XbQCD[mb,mW,mZ,mH,mt,scale,nL,nH])^2 /. {lal -> aEW[scale], las->aQCD[scale]}];

:Evaluate:  Options[alam] = { "IncludedCorrections" -> { 
						       EW[1],  
						       EW[2],
						       MIXED[1,1]
						     }};

:Evaluate:  alam[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1,OptionsPattern[]] :=Block[ 
				{lal,las,crulesh = Map[(#->yH[Sequence @@ # ]) &, OptionValue["IncludedCorrections"]] }, 
				2^(-1/2)*Gf/(4*Pi)^2*mH^2*(1+FromCoefficientRules[crulesh,{lal,las}])/.XH[mb,mW,mZ,mH,mt,scale,nL,nH] /. {lal -> aEW[scale], las->aQCD[scale]} 
						];

:Evaluate:  Options[alphaGF] = { "IncludedCorrections" -> { 
						       EW[1],  
						       EW[2],
						       MIXED[1,1]
						     }};

:Evaluate:  alphaGF[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1,OptionsPattern[] ] := Module[{lal,las,
				alGF = 2^(1/2)*Gf*mW^2*(1-mW^2/mZ^2)/Pi, crulesAGF = Map[(#->daGF[Sequence @@ # ]) &, OptionValue["IncludedCorrections"]]},
					alGF*(1 + FromCoefficientRules[crulesAGF,{lal,las}]) /. lal -> alGF/(4 Pi) /. las->aQCD[scale] /. dalphaGF[mb,mW,mZ,mH,mt,scale,nL,nH]
						];

:Evaluate:  Options[vev] = { "IncludedCorrections" -> { 
						       EW[1],  
						       EW[2],
						       MIXED[1,1]
						     }};

:Evaluate:  vev[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1, OptionsPattern[]] := 
				Block[ {lal,las,crulesvev = Map[(#->dr[Sequence @@ # ]) &, OptionValue["IncludedCorrections"]] },
			2^(-1/4)/Sqrt[Gf]*Sqrt[(1+FromCoefficientRules[crulesvev,{lal,las}])] /. {lal -> aEW[scale], las->aQCD[scale]} /. dROS[mb,mW,mZ,mH,mt,scale,nL,nH]
						];

// Pole masses

:Evaluate:  Options[MW] = { "IncludedCorrections" -> { 
						       EW[1],  
						       EW[2],
						       MIXED[1,1]
						     }};

:Evaluate:  MW[G1_?NumericQ,G2_?NumericQ,GS__?NumericQ,YB_?NumericQ,YT_?NumericQ,LAM_?NumericQ,M_?NumericQ,SC_?NumericQ,OptionsPattern[]] := Block[{lal,las,vev,lc,aEW,aQCD,mW,res, 
				cr = Map[(#->xMMW[Sequence @@ # ]) &, OptionValue["IncludedCorrections"]],
				pars = {G1,G2,GS,YB,YT,LAM,M,SC}}, (* loop corrections *)	lc = XMMW[ Sequence @@ pars];
				vev = Sqrt[M^2/LAM/2]; mW = G2 * vev/2.0;
				{aEW, aQCD} = {G1^2*G2^2/(G1^2 + G2^2), GS^2}/(16 Pi^2);
				res = mW * Sqrt[1 + FromCoefficientRules[cr, {lal,las}]] /. lc  /. {lal->aEW, las -> aQCD};
				Return[ res ]
			];	


:Evaluate:  MW[runpars_List, opt:OptionsPattern[]]:= Block[{pars = {g1,g2,gs,yb,yt,lam,m,scale} /. runpars},
			(* check numeric *) If [ And @@ NumericQ /@ pars, 
			Return[ MW[ Sequence @@ pars, opt ] ],
			(* else *) Print[" Not All parameters specified ", pars, " from ", runpars]]];

:Evaluate:  Options[MZ] = { "IncludedCorrections" -> { 
						       EW[1],  
						       EW[2],
						       MIXED[1,1]
						     }};


:Evaluate:  MZ[G1_?NumericQ,G2_?NumericQ,GS__?NumericQ,YB_?NumericQ,YT_?NumericQ,LAM_?NumericQ,M_?NumericQ,SC_?NumericQ, OptionsPattern[]] := Block[{lal,las,vev,lc,aEW,aQCD,mZ, res,
			cr = Map[(#->xMMZ[Sequence @@ # ]) &, OptionValue["IncludedCorrections"]],
			pars = {G1,G2,GS,YB,YT,LAM,M,SC} },
			(* loop corrections *)	lc = XMMZ[ Sequence @@ pars];
			vev = Sqrt[M^2/LAM/2]; mZ = Sqrt[G1^2 + G2^2] * vev/2.0;
			{aEW, aQCD} = {G1^2*G2^2/(G1^2 + G2^2), GS^2}/(16 Pi^2);
			res = mZ * Sqrt[1 + FromCoefficientRules[cr, {lal,las}]] /. lc  /. {lal->aEW, las -> aQCD};
			Return[res]]	

:Evaluate:  MZ[runpars_List, opt:OptionsPattern[]]:= Block[{pars = {g1,g2,gs,yb,yt,lam,m,scale} /. runpars},
			(* check numeric *) If [ And @@ NumericQ /@ pars, 
			Return[ MZ[ Sequence @@ pars, opt ] ],
			(* else *) Print[" Not All parameters specified ", pars, " from ", runpars]]];



:Evaluate:  Options[MH] = { "IncludedCorrections" -> { 
						       EW[1],  
						       EW[2],
						       MIXED[1,1]
						     }};

:Evaluate:  MH[G1_?NumericQ,G2_?NumericQ,GS__?NumericQ,YB_?NumericQ,YT_?NumericQ,LAM_?NumericQ,M_?NumericQ,SC_?NumericQ, OptionsPattern[]] := Block[{lal,las,vev,lc,aEW,aQCD,res, 
			cr = Map[(#->xMMH[Sequence @@ # ]) &, OptionValue["IncludedCorrections"]],
			pars = {G1,G2,GS,YB,YT,LAM,M,SC} },
			(* loop corrections *)	lc = XMMH[ Sequence @@ pars];
			{aEW, aQCD} = {G1^2*G2^2/(G1^2 + G2^2), GS^2}/(16 Pi^2);
			res = M * Sqrt[1 + FromCoefficientRules[cr, {lal,las}]] /. lc  /. {lal->aEW, las -> aQCD};
			Return[res]]	

:Evaluate:  MH[runpars_List, opt:OptionsPattern[]]:= Block[{pars = {g1,g2,gs,yb,yt,lam,m,scale} /. runpars},
			(* check numeric *) If [ And @@ NumericQ /@ pars, 
			Return[ MH[ Sequence @@ pars, opt ] ],
			(* else *) Print[" Not All parameters specified ", pars, " from ", runpars]]];



:Evaluate:  Options[MT] = { "IncludedCorrections" -> { 
						      QCD[1],  
						      QCD[2],  
						      QCD[3],  
						       EW[1],  
						       EW[2],
						       MIXED[1,1]
						     }};

:Evaluate:  MT[G1_?NumericQ,G2_?NumericQ,GS__?NumericQ,YB_?NumericQ,YT_?NumericQ,LAM_?NumericQ,M_?NumericQ,SC_?NumericQ, OptionsPattern[]] := Block[{lal,las,vev,lc,aEW,aQCD,mt,res, 
			cr = Map[(#->xMT[Sequence @@ # ]) &, OptionValue["IncludedCorrections"]],
			pars = {G1,G2,GS,YB,YT,LAM,M,SC} },
			(* loop corrections *)	lc = Join[XMT[ Sequence @@ pars], XMTQCD[ Sequence @@ pars]];
			{aEW, aQCD} = {G1^2*G2^2/(G1^2 + G2^2), GS^2}/(16 Pi^2);
			vev = Sqrt[M^2/LAM/2]; mt = YT * vev / Sqrt[2];
			res = mt * (1 + FromCoefficientRules[cr, {lal,las}]) /. lc  /. {lal->aEW, las -> aQCD};
			Return[ res ];
			];	

:Evaluate:  MT[runpars_List, opt:OptionsPattern[]]:= Block[{pars = {g1,g2,gs,yb,yt,lam,m,scale} /. runpars},
			(* check numeric *) If [ And @@ NumericQ /@ pars, 
			Return[ MT[ Sequence @@ pars, opt ] ],
			(* else *) Print[" Not All parameters specified ", pars, " from ", runpars]]];


:Evaluate:  Options[GF] = { "IncludedCorrections" -> { 
						       EW[1],  
						       EW[2],
						       MIXED[1,1]
						     }};

:Evaluate:  GF[G1_?NumericQ,G2_?NumericQ,GS__?NumericQ,YB_?NumericQ,YT_?NumericQ,LAM_?NumericQ,M_?NumericQ,SC_?NumericQ, OptionsPattern[]] := Block[{res,lal,las,vv,lc,aEW,aQCD,gf, 
			cr = Map[(#->xdRbar[Sequence @@ # ]) &, OptionValue["IncludedCorrections"]], pars = {G1,G2,GS,YB,YT,LAM,M,SC} },
			(* loop corrections *)	lc = XdRbar[ Sequence @@ pars];
			{aEW, aQCD} = {G1^2*G2^2/(G1^2 + G2^2), GS^2}/(16 Pi^2);
			vv = M^2/LAM/2; gf = 1/Sqrt[2]/vv; 
			res = gf * (1 + FromCoefficientRules[cr, {lal,las}])  /. lc /. {lal->aEW, las -> aQCD}; 
			Return[res];
			];
				
:Evaluate:  GF[runpars_List, opt:OptionsPattern[]]:= Block[{pars = {g1,g2,gs,yb,yt,lam,m,scale} /. runpars},
			(* check numeric *) If [ And @@ NumericQ /@ pars, 
			Return[ GF[ Sequence @@ pars, opt ] ],
			(* else *) Print[" Not All parameters specified ", pars, " from ", runpars]]];



:Evaluate: RunSM[runpars_List, oscale_?NumericQ] :=  Block[{pars = {g1,g2,gs,yb,yt,lam,m,scale} /. runpars, res}, 
			(* check numeric *) If [ And @@ NumericQ /@ pars, 
			res = RunSM[ Sequence @@ pars, oscale];
			res = MapThread[ #1 -> #2 & , {{g1,g2,gs,yb,yt,lam,m,scale},res}];
			Return[ res ],
			(* else *) Print["RunSM: Not All parameters specified " , pars, " from ", runpars ]]];
:Evaluate: RunSMcouplings[runpars_List, oscale_?NumericQ, L_Integer] :=  Block[{pars = {g1,g2,gs,yb,yt,lam,scale} /. runpars, res}, 
			(* check numeric *) If [ And @@ NumericQ /@ pars, 
			res = RunSMcouplings[ Sequence @@ pars, oscale, L];
			res = MapThread[ #1 -> #2 & , {{g1,g2,gs,yb,yt,lam,scale},res}];
			Return[ res ],
			(* else *) Print["RunSMcouplings: Not All parameters specified " , pars, " from ", runpars ]]];
		
// C++ part
:Begin:
:Function: XMMW
:Pattern: XMMW[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale]}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128}
:ReturnType: Manual
:End:


:Begin:
:Function: XMMZ
:Pattern: XMMZ[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale]}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128}
:ReturnType: Manual
:End:

:Begin:
:Function: XMMH
:Pattern: XMMH[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale]}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128}
:ReturnType: Manual
:End:

:Begin:
:Function: XdRbar
:Pattern: XdRbar[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale]}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128}
:ReturnType: Manual
:End:


:Begin:
:Function: XMT
:Pattern: XMT[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale]}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128}
:ReturnType: Manual
:End:

:Begin:
:Function: XMTQCD
:Pattern: XMTQCD[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,scale_?NumericQ]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[scale]}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128}
:ReturnType: Manual
:End:








:Begin:
:Function: RunQCDnf6
:Pattern: RunQCDnf6[oscale_?NumericQ,asMZ_?NumericQ,MZscale_?NumericQ,nL_Integer,mtpole_?NumericQ,mtth_?NumericQ] /; 1<=nL<=4
:Arguments: {N[oscale],N[asMZ],N[MZscale],nL,N[mtpole],N[mtth]}
:ArgumentTypes: {Real128,Real128,Real128,Integer,Real128,Real128}
:ReturnType: Manual
:End:

:Begin:
:Function: RunQCD
:Pattern: RunQCD[oscale_?NumericQ,as0_?NumericQ,inscale_?NumericQ,nL_Integer,nF_Integer] /; 1<=nL<=4
:Arguments: {N[oscale],N[as0],N[inscale],nL,nF}
:ArgumentTypes: {Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: RunSM
:Pattern: RunSM[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,m_?NumericQ,iscale_?NumericQ,oscale_?NumericQ]
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[m],N[iscale],N[oscale]}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128}
:ReturnType: Manual
:End:

:Begin:
:Function: RunSMcouplings
:Pattern: RunSMcouplings[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yb_?NumericQ,yt_?NumericQ,lam_?NumericQ,iscale_?NumericQ,oscale_?NumericQ, L_Integer] /; 1<=L<=3
:Arguments: {N[gp],N[g],N[gs],N[yb],N[yt],N[lam],N[iscale],N[oscale],L}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128, Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: XW
:Pattern: XW[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[scale],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: XZ
:Pattern: XZ[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[scale],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: XH
:Pattern: XH[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[scale],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: Xt
:Pattern: Xt[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[scale],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: dROS
:Pattern: dROS[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[scale],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: dalphaGF
:Pattern: dalphaGF[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[scale],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: Xb
:Pattern: Xb[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[scale],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: XtQCD
:Pattern: XtQCD[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[scale],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: XbQCD
:Pattern: XbQCD[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,scale_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[scale],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:


:Evaluate: End[]
:Evaluate: EndPackage[]
