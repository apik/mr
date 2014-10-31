:Evaluate:  Print["MR - MAtching and Running Mathmatica interface"];
:Evaluate:  Print["To see available functions use Names[\"mr`*\"]"];
:Evaluate:  Print["Andrey Pikelner <pikelner@theor.jinr.ru>"];

:Evaluate:  BeginPackage["mr`"]

:Evaluate:  Xb::usage  = "Xb[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "

:Evaluate:  XW::usage  = "XW[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "

:Evaluate:  XZ::usage  = "XZ[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "

:Evaluate:  XH::usage  = "XH[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "

:Evaluate:  Xt::usage  = "Xt[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  Input is in terms of pole masses and matching scale, nL and nH are number of light and heavy quark genrations "

:Evaluate:  XbQCD::usage  = "XbQCD[Mb,MW,MZ,MH,Mt,mu,nf=5]  Pure QCD corrections, nf is a number of light quarks, default is 5." 

:Evaluate:  XtQCD::usage  = "XtQCD[Mb,MW,MZ,MH,Mt,mu,nf=5]  Pure QCD corrections, nf is a number of light quarks, default is 5." 

:Evaluate:  mmWMMW::usage  = "mmWMMW[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  mW^2/MW^2, full correction to relation between MS-bar mass mW and pole MW"

:Evaluate:  mmZMMZ::usage  = "mmZMMZ[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  mZ^2/MZ^2, full correction to relation between MS-bar mass mZ and pole MZ"

:Evaluate:  mmHMMH::usage  = "mmHMMH[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  mH^2/MH^2, full correction to relation between MS-bar mass mH and pole MH"

:Evaluate:  mtMt::usage  = "mtMt[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  mt/Mt, full correction to relation between MS-bar mass mt and pole Mt"

:Evaluate:  mbMb::usage  = "mbMb[Mb,MW,MZ,MH,Mt,mu,nL=2,nH=1]  mb/Mb, full correction to relation between MS-bar mass mb and pole Mb"

:Evaluate:  aEW::usage  = "alpha/4/Pi, running EW constant"

:Evaluate:  aQCD::usage  = "alpha_S/4/Pi, running QCD constant"

:Evaluate:  xW::usage  = ""
:Evaluate:  xZ::usage  = ""
:Evaluate:  xH::usage  = ""
:Evaluate:  xt::usage  = ""
:Evaluate:  xb::usage  = ""

:Evaluate:  Begin["`Private`"]

// Mathematica part

:Evaluate:  mmWMMW[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1] := (1+aEW[mu]*xW[1,0]+aEW[mu]*aQCD[mu]*xW[1,1]+aEW[mu]^2*xW[2,0])/.XW[mb,mW,mZ,mH,mt,mu,nH,nL];

:Evaluate:  mmZMMZ[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1] := (1+aEW[mu]*xZ[1,0]+aEW[mu]*aQCD[mu]*xZ[1,1]+aEW[mu]^2*xZ[2,0])/.XZ[mb,mW,mZ,mH,mt,mu,nH,nL];

:Evaluate:  mmHMMH[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1] := (1+aEW[mu]*xH[1,0]+aEW[mu]*aQCD[mu]*xH[1,1]+aEW[mu]^2*xH[2,0])/.XH[mb,mW,mZ,mH,mt,mu,nH,nL];
                                                                                                                             
:Evaluate:  mtMt[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1] := (1+aEW[mu]*xt[1,0]+aEW[mu]*aQCD[mu]*xt[1,1]+aEW[mu]^2*xt[2,0])/.Xt[mb,mW,mZ,mH,mt,mu,nH,nL];

:Evaluate:  mbMb[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1] := (1+aEW[mu]*xb[1,0]+aEW[mu]*aQCD[mu]*xb[1,1]+aEW[mu]^2*xb[2,0])/.Xb[mb,mW,mZ,mH,mt,mu,nH,nL];
 
// C++ part
:Begin:
:Function: XW
:Pattern: XW[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[mu],nL,nH}
:ArgumentTypes: {Real64,Real64,Real64,Real64,Real64,Real64,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: XZ
:Pattern: XZ[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[mu],nL,nH}
:ArgumentTypes: {Real64,Real64,Real64,Real64,Real64,Real64,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: XH
:Pattern: XH[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[mu],nL,nH}
:ArgumentTypes: {Real64,Real64,Real64,Real64,Real64,Real64,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: Xt
:Pattern: Xt[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[mu],nL,nH}
:ArgumentTypes: {Real64,Real64,Real64,Real64,Real64,Real64,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: Xb
:Pattern: Xb[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nL_Integer:2,nH_Integer:1]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[mu],nL,nH}
:ArgumentTypes: {Real64,Real64,Real64,Real64,Real64,Real64,Integer,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: XtQCD
:Pattern: XtQCD[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nf_Integer:5]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[mu],nf}
:ArgumentTypes: {Real64,Real64,Real64,Real64,Real64,Real64,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: XbQCD
:Pattern: XbQCD[mb_?NumericQ,mW_?NumericQ,mZ_?NumericQ,mH_?NumericQ,mt_?NumericQ,mu_?NumericQ,nf_Integer:5]
:Arguments: {N[mb],N[mW],N[mZ],N[mH],N[mt],N[mu],nf}
:ArgumentTypes: {Real64,Real64,Real64,Real64,Real64,Real64,Integer}
:ReturnType: Manual
:End:


:Evaluate: End[]
:Evaluate: EndPackage[]