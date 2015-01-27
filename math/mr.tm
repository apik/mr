:Evaluate:  Print["MR - Matching and Running Mathmatica interface"];
:Evaluate:  Print["To see available functions use Names[\"mr`*\"]"];
:Evaluate:  Print["Andrey Pikelner <pikelner@theor.jinr.ru>"];

:Evaluate:  BeginPackage["mr`"]

:Evaluate:  RunQCD::usage  = "RunQCD[oscale,asMZ, MZscale,L=4, mtth] evaluate alphas at scale oscale given asMZ at MZscale (nf=5) with L-loop RGE with top threshold at mtth"
:Evaluate:  MW::usage  = "MW[gp,g,gs,yt,lam,mu0,mu,L=2] returns pole W-boson mass MW given  MSbar parameters at scale mu at L-loop level"
:Evaluate:  MWp::usage  = "MWp[gp,g,gs,yt,lam,mu0,mu,L=2] returns pole W-boson mass MW given  MSbar parameters at scale mu at L-loop level"
:Evaluate:  MZp::usage  = "MZp[gp,g,gs,yt,lam,mu0,mu,L=2] returns pole Z-boson mass MZ given  MSbar parameters at scale mu at L-loop level"
:Evaluate:  MZ::usage  = "MZ[gp,g,gs,yt,lam,mu0,mu,L=2] returns pole Z-boson mass MZ given  MSbar parameters at scale mu at L-loop level"
:Evaluate:  MH::usage  = "MH[gp,g,gs,yt,lam,mu0,mu,L=2] returns pole H-boson mass MH given  MSbar parameters at scale mu at L-loop level"
:Evaluate:  MHp::usage  = "MHp[gp,g,gs,yt,lam,mu0,mu,L=2] returns pole H-boson mass MH given  MSbar parameters at scale mu at L-loop level"
:Evaluate:  MT::usage  = "MT[gp,g,gs,yt,lam,mu0,mu,L=2] returns pole t-quark mass MT given  MSbar parameters at scale mu at L-loop level"
:Evaluate:  MTp::usage  = "MTp[gp,g,gs,yt,lam,mu0,mu,L=2] returns pole t-quark mass MT given  MSbar parameters at scale mu at L-loop level"
:Evaluate:  GF::usage  = "GF[gp,g,gs,yt,lam,mu0,mu,L=2] returns Fermi constant GF given  MSbar parameters at scale mu at L-loop level"

:Evaluate:  RunSM::usage  = "RunSM[gp,g,gs,yt,lam,mu0,iscale,oscale] return running parameters at scale oscale given the values at scale iscale"

:Evaluate:  Xb::usage  = "Xb[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "

:Evaluate:  XW::usage  = "XW[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "

:Evaluate:  XZ::usage  = "XZ[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "

:Evaluate:  XH::usage  = "XH[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "

:Evaluate:  Xt::usage  = "Xt[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "

:Evaluate:  dROS::usage  = "dROS[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "

:Evaluate:  XbQCD::usage  = "XbQCD[Mb,MW,MZ,MH,Mt,mu,nf=5]  Pure QCD corrections, nf is a number of light quarks, default is 5." 

:Evaluate:  XtQCD::usage  = "XtQCD[Mb,MW,MZ,MH,Mt,mu,nf=5]  Pure QCD corrections, nf is a number of light quarks, default is 5." 

:Evaluate:  mmWMMW::usage  = "mmWMMW[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  mW^2/MW^2, full correction to relation between MS-bar mass mW and pole MW"

:Evaluate:  mmZMMZ::usage  = "mmZMMZ[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  mZ^2/MZ^2, full correction to relation between MS-bar mass mZ and pole MZ"

:Evaluate:  mmHMMH::usage  = "mmHMMH[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  mH^2/MH^2, full correction to relation between MS-bar mass mH and pole MH"

:Evaluate:  mtMt::usage  = "mtMt[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  mt/Mt, full correction to relation between MS-bar mass mt and pole Mt"

:Evaluate:  mbMb::usage  = "mbMb[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  mb/Mb, full correction to relation between MS-bar mass mb and pole Mb"

:Evaluate:  a1::usage  = "a1[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  a1=5/3(g1/(4\[Pi]))^2, we use GUT normalization"

:Evaluate:  a2::usage  = "a2[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  a2=(g2/(4\[Pi]))^2"

:Evaluate:  at::usage  = "a1[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  at=(yt/(4\[Pi]))^2"

:Evaluate:  ab::usage  = "ab[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  ab=(yb/(4\[Pi]))^2"

:Evaluate:  alam::usage  = "alam[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  alam=\[Lambda]/(4\[Pi])^2"

:Evaluate:  vev::usage  = "vev[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  vev - running vev at scale mu"

:Evaluate:  aEW::usage  = "\[Alpha]/(4\[Pi]), running EW constant"

:Evaluate:  aQCD::usage  = "\[Alpha]_S/(4\[Pi]), running QCD constant"

:Evaluate:  Gf::usage  = "G_F Fermi constant"
:Evaluate:  h::usage  = "h - loop counter"

:Evaluate:  xW::usage  = ""
:Evaluate:  xZ::usage  = ""
:Evaluate:  xH::usage  = ""
:Evaluate:  xt::usage  = ""
:Evaluate:  xb::usage  = ""

:Evaluate:  yW::usage  = ""
:Evaluate:  yZ::usage  = ""
:Evaluate:  yH::usage  = ""
:Evaluate:  yt::usage  = ""
:Evaluate:  yb::usage  = ""

:Evaluate:  dr::usage  = ""

:Evaluate:  Begin["`Private`"]

// Mathematica part

// Masses
:Evaluate:  mmWMMW[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1] := (1+aEW[mu]*xW[1,0]+aEW[mu]*aQCD[mu]*xW[1,1]+aEW[mu]^2*xW[2,0])/.XW[mb,mW,mZ,mH,mt,mu,nL,nH];

:Evaluate:  mmZMMZ[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1] := (1+aEW[mu]*xZ[1,0]+aEW[mu]*aQCD[mu]*xZ[1,1]+aEW[mu]^2*xZ[2,0])/.XZ[mb,mW,mZ,mH,mt,mu,nL,nH];

:Evaluate:  mmHMMH[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1] := (1+aEW[mu]*xH[1,0]+aEW[mu]*aQCD[mu]*xH[1,1]+aEW[mu]^2*xH[2,0])/.XH[mb,mW,mZ,mH,mt,mu,nL,nH];
                                                                                                                             
:Evaluate:  mtMt[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1] := (1+aEW[mu]*xt[1,0]+aEW[mu]*aQCD[mu]*xt[1,1]+aEW[mu]^2*xt[2,0])/.Xt[mb,mW,mZ,mH,mt,mu,nL,nH];

:Evaluate:  mbMb[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1] := (1+aEW[mu]*xb[1,0]+aEW[mu]*aQCD[mu]*xb[1,1]+aEW[mu]^2*xb[2,0])/.Xb[mb,mW,mZ,mH,mt,mu,nL,nH];

// Couplings

:Evaluate:  a1[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1] := 5/3*2^(5/2)*Gf/(4*Pi)^2*((mZ^2*(1+aEW[mu]*yZ[1,0]+aEW[mu]*aQCD[mu]*yZ[1,1]+aEW[mu]^2*yZ[2,0])/.XZ[mb,mW,mZ,mH,mt,mu,nL,nH]) - (mW^2*(1+aEW[mu]*yW[1,0]+aEW[mu]*aQCD[mu]*yW[1,1]+aEW[mu]^2*yW[2,0])/.XW[mb,mW,mZ,mH,mt,mu,nL,nH]));

:Evaluate:  a2[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1] := 2^(5/2)*Gf/(4*Pi)^2*mW^2*(1+aEW[mu]*yW[1,0]+aEW[mu]*aQCD[mu]*yW[1,1]+aEW[mu]^2*yW[2,0])/.XW[mb,mW,mZ,mH,mt,mu,nL,nH];

:Evaluate:  at[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1] := 2^(3/4)*Sqrt[Gf]/(4*Pi)^2*mt*(1+aQCD[mu]*xt[0,1] + aEW[mu]*yt[1,0]+aEW[mu]*aQCD[mu]*yt[1,1]+aEW[mu]^2*yt[2,0] + aQCD[mu]^2*xt[0,2] + aQCD[mu]^3*xt[0,3])/.Xt[mb,mW,mZ,mH,mt,mu,nL,nH] /. XtQCD[mb,mW,mZ,mH,mt,mu,nL,nH];

:Evaluate:  ab[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1] := 2^(3/4)*Sqrt[Gf]/(4*Pi)^2*mb*(1+aQCD[mu]*xb[0,1]+aEW[mu]*yb[1,0]+aEW[mu]*aQCD[mu]*yb[1,1]+aEW[mu]^2*yb[2,0] + aQCD[mu]^2*xb[0,2] + aQCD[mu]^3*xb[0,3])/.Xb[mb,mW,mZ,mH,mt,mu,nL,nH] /. XbQCD[mb,mW,mZ,mH,mt,mu,nL,nH];

:Evaluate:  alam[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1] := 2^(-1/2)*Gf/(4*Pi)^2*mH^2*(1+aEW[mu]*yH[1,0]+aEW[mu]*aQCD[mu]*yH[1,1]+aEW[mu]^2*yH[2,0])/.XH[mb,mW,mZ,mH,mt,mu,nL,nH];

:Evaluate:  vev[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1] := 2^(-1/4)/Sqrt[Gf]*Sqrt[(1+aEW[mu]*dr[1,0]+aEW[mu]*aQCD[mu]*dr[1,1]+aEW[mu]^2*dr[2,0])] /.dROS[mb,mW,mZ,mH,mt,mu,nL,nH];

// C++ part
:Begin:
:Function: MW
:Pattern: MW[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yt_?NumericQ,lam_?NumericQ,mu0_?NumericQ,mu_?NumericQ,L_Integer:2]
:Arguments: {N[gp],N[g],N[gs],N[yt],N[lam],N[mu0],N[mu],L}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: MWp
:Pattern: MWp[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yt_?NumericQ,lam_?NumericQ,mu0_?NumericQ,mu_?NumericQ,L_Integer:2]
:Arguments: {N[gp],N[g],N[gs],N[yt],N[lam],N[mu0],N[mu],L}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: MZ
:Pattern: MZ[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yt_?NumericQ,lam_?NumericQ,mu0_?NumericQ,mu_?NumericQ,L_Integer:2]
:Arguments: {N[gp],N[g],N[gs],N[yt],N[lam],N[mu0],N[mu],L}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: MZp
:Pattern: MZp[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yt_?NumericQ,lam_?NumericQ,mu0_?NumericQ,mu_?NumericQ,L_Integer:2]
:Arguments: {N[gp],N[g],N[gs],N[yt],N[lam],N[mu0],N[mu],L}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: MHp
:Pattern: MHp[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yt_?NumericQ,lam_?NumericQ,mu0_?NumericQ,mu_?NumericQ,L_Integer:2]
:Arguments: {N[gp],N[g],N[gs],N[yt],N[lam],N[mu0],N[mu],L}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: MH
:Pattern: MH[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yt_?NumericQ,lam_?NumericQ,mu0_?NumericQ,mu_?NumericQ,L_Integer:2]
:Arguments: {N[gp],N[g],N[gs],N[yt],N[lam],N[mu0],N[mu],L}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: MT
:Pattern: MT[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yt_?NumericQ,lam_?NumericQ,mu0_?NumericQ,mu_?NumericQ,L_Integer:2]
:Arguments: {N[gp],N[g],N[gs],N[yt],N[lam],N[mu0],N[mu],L}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: MTp
:Pattern: MTp[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yt_?NumericQ,lam_?NumericQ,mu0_?NumericQ,mu_?NumericQ,L_Integer:2]
:Arguments: {N[gp],N[g],N[gs],N[yt],N[lam],N[mu0],N[mu],L}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: GF
:Pattern: GF[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yt_?NumericQ,lam_?NumericQ,mu0_?NumericQ,mu_?NumericQ,L_Integer:2]
:Arguments: {N[gp],N[g],N[gs],N[yt],N[lam],N[mu0],N[mu],L}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: RunQCD
:Pattern: RunQCD[oscale_?NumericQ,asMZ_?NumericQ,MZscale_?NumericQ,nL_Integer,mtth_?NumericQ]
:Arguments: {N[oscale],N[asMZ],N[MZscale],nL,N[mtth]}
:ArgumentTypes: {Real128,Real128,Real128,Integer,Real128}
:ReturnType: Manual
:End:

:Begin:
:Function: RunSM
:Pattern: RunSM[gp_?NumericQ,g_?NumericQ,gs_?NumericQ,yt_?NumericQ,lam_?NumericQ,mu0_?NumericQ,iscale_?NumericQ,oscale_?NumericQ]
:Arguments: {N[gp],N[g],N[gs],N[yt],N[lam],N[mu0],N[iscale],N[oscale]}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Real128,Real128}
:ReturnType: Manual
:End:

:Begin:
:Function: XW
:Pattern: XW[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[mu],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: XZ
:Pattern: XZ[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[mu],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: XH
:Pattern: XH[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[mu],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: Xt
:Pattern: Xt[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[mu],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: dROS
:Pattern: dROS[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[mu],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: Xb
:Pattern: Xb[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[mu],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: XtQCD
:Pattern: XtQCD[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[mu],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: XbQCD
:Pattern: XbQCD[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[mu],nL,nH}
:ArgumentTypes: {Real128,Real128,Real128,Real128,Real128,Real128,Integer,Integer}
:ReturnType: Manual
:End:


:Evaluate: End[]
:Evaluate: EndPackage[]
