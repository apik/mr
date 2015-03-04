link=Install["mr"];

(* PDG 2014 data *)

PDG`MW = 80.385;
PDG`dMW = 0.015;
PDG`MZ = 91.1876;
PDG`dMZ = 0.0021;
PDG`MH = 125.7;
PDG`dMH = 0.4;
PDG`MT = 173.21;
PDG`dMT = 1.22;
PDG`GF =  0.000011663787;
PDG`dGF = 0.000000000006; 
PDG`asQCD = 0.1185;

(* pole mass - two-loop eq from MS *)

PDG`Mb = 4.78;
PDG`dMb = 0.06;

(* 3 sigma range [x - 3 s, x + 3 s] divided into n subintervals *)

ThreeSigmaRange[x_,s_,n_:2] := Module[{step = 6 s /n}, Table[ x - 3 s + step * (i-1), {i,1,n+1}]] 


(* find running parameters given (pseudo)observables, and a renormalization scale *)

sol[scale_][mzp_:PDG`MZ,mwp_:PDG`MW,mtp_:PDG`MT,mhp_:PDG`MH,gfp_:PDG`GF,asp_:PDG`asQCD] := Join[FindRoot[ {
	   MZ[gp, g, gs, yt, lam, mu0, scale,2] == mzp (* 91.1876 *),
	   MW[gp, g, gs, yt, lam, mu0, scale,2] == mwp (* 80.385 *),
	   MT[gp, g, gs, yt, lam, mu0, scale,2] == mtp (* scale *),
	   MH[gp, g, gs, yt, lam, mu0, scale,2] == mhp (* 125.7 *),
	   GF[gp, g, gs, yt, lam, mu0, scale,2] == gfp (* 0.000011663787 *),
	   gs^2/(4*Pi)==RunQCD[scale, asp,PDG`MZ,4, PDG`MT] 
	  },{
		{gp, 0.357561},
		{g , 0.64822},
		{gs , 1.1666},
		{yt,0.93558},
		{lam,0.12711},
		{mu0,132.03}
	  }],{mu -> scale}];

(* evaluate chi2 function for given set of running parameters at the given scale *)

chiSquare[sol_,scale_]:=(   
		(MZ[gp, g, gs, yt, lam, mu0, scale,2] - PDG`MZ)^2/(PDG`dMZ)^2 + 
	 	(MW[gp, g, gs, yt, lam, mu0, scale,2] - PDG`MW)^2/(PDG`dMW)^2 + 
		(MT[gp, g, gs, yt, lam, mu0, scale,2] - PDG`MT)^2/(PDG`dMT)^2 + 
		(MH[gp, g, gs, yt, lam, mu0, scale,2] - PDG`MH)^2/(PDG`dMH)^2 + 
		(GF[gp, g, gs, yt, lam, mu0, scale,2] - PDG`GF)^2/(PDG`dGF)^2 (* check this *)	 + 
		(gs^2/(4 Pi) -  RunQCD[scale, PDG`asQCD,PDG`MZ,4, PDG`MT])^2/(0.0001)^2
			) /. sol ;

(* theoretical error for observables *)
EstimateTheorUncertaintyInMatchingDegrassi[scale_, scalefactor_:10] := Module[ {
	  		ref = sol[scale][], 
			refe,
	  		ppp = sol[scale scalefactor][],
			pppe,
	  		mmm = sol[scale/scalefactor][],
			mmme
	 		 },
			refe = Last /@ Delete[RunParsFromPoleMassesAndAsWithMb[gs^2/(4 Pi) /. ref, scale] ,{{4},{8}}] /. z_ corr[___,x_]:> x ;
			pppe = Last /@ Delete[RunParsFromPoleMassesAndAsWithMb[gs^2/(4 Pi) /. ppp, scale] ,{{4},{8}}] /. z_ corr[___,x_]:> x ;
			mmme = Last /@ Delete[RunParsFromPoleMassesAndAsWithMb[gs^2/(4 Pi) /. mmm, scale] ,{{4},{8}}] /. z_ corr[___,x_]:> x ;
Print["refe=", refe];
Print["pppe=", pppe];
Print["mmme=", mmme];
			ppp = RunSM[ Sequence @@ (Last /@ ppp), scale]; (* scale*scalefactor -> scale *)
			mmm = RunSM[ Sequence @@ (Last /@ mmm), scale]; (* scale/scalefactor -> scale *)
			MapThread[ #1 /. Rule[a_,b_]:> Rule[a, {b,1-#2/b,1-#3/b}] &,{ref,ppp,mmm}]	

			] 

EstimateTheorUncertainty[ gp_, g_, gs_, yt_, lam_, mu0_,scale_, loop_:2, scalefactor_:2] := Module[ { 
			(*
			refMZ = MZ[gp, g, gs, yt, lam, mu0, scale,loop],
			refMW = MW[gp, g, gs, yt, lam, mu0, scale,loop],
			refMT = MT[gp, g, gs, yt, lam, mu0, scale,loop],
			refMH = MH[gp, g, gs, yt, lam, mu0, scale,loop],
			refGF = GF[gp, g, gs, yt, lam, mu0, scale,loop],
			*)
			parsP = Flatten[{RunSM[gp, g, gs, yt, lam, mu0, scale, scalefactor  * scale],  loop}],
			parsM = Flatten[{RunSM[gp, g, gs, yt, lam, mu0, scale, 1/scalefactor * scale],  loop}],
			ppp, mmm, ref
			},
			Print[ parsP];
			ref = ({
					MZ[Sequence @@ #],
					MW[Sequence @@ #],
					MT[Sequence @@ #],
					MH[Sequence @@ #],
					GF[Sequence @@ #]} &[{gp,g,gs,yt,lam,mu0,scale,loop}]);
			ppp = ({
					MZ[Sequence @@ #],
					MW[Sequence @@ #],
					MT[Sequence @@ #],
					MH[Sequence @@ #],
					GF[Sequence @@ #]} &[parsP]);
			mmm = ({
					MZ[Sequence @@ #],
					MW[Sequence @@ #],
					MT[Sequence @@ #],
					MH[Sequence @@ #],
					GF[Sequence @@ #]} &[parsM]);
			Return[MapThread[List[##] &,{ref,ppp/ref-1,mmm/ref-1}]]
			];


ScaleDependence[O_][runpars_,loop_:2] := O[ Sequence @@ (RunSM[ Sequence @@ Evaluate[{gp, g, gs, yt, lam, mu0, mu} /. runpars], #]),loop] & 

RunParsFromPoleMassesAndAs[mz_:PDG`MZ,mw_:PDG`MW,mh_:PDG`MH,mt_:PDG`MT,gf_:PDG`GF,asmu_,smu_] := Module[{aew,seq,aa1,aa2, impliciteq},
		aa1 = (4 Pi) 3/5 * a1[0,mw,mz,mh,mt,smu] /.aQCD[smu]->asmu/(4 Pi) /. Gf->gf /. aEW[smu]->aew/(4 Pi) ;
		aa2 = (4 Pi)       a2[0,mw,mz,mh,mt,smu] /.aQCD[smu]->asmu/(4 Pi) /. Gf->gf /. aEW[smu]->aew/(4 Pi) ;
		impliciteq = (aew == (Simplify[aa1*aa2/(aa1 + aa2)]));
		impliciteqcheck = (aew == Normal[Series[Simplify[aa1*aa2/(aa1 + aa2)], {aew,0,2}]]);
		(* solution for alpha ew at the scale *)
		solaew = FindRoot[ impliciteq, {aew, 1/127.94}]; 
		(* Couplings *)
		{ gp->Sqrt[ (4 Pi)^2 3/5 a1[0, mw,mz,mh,mt,smu]], 
		  g -> Sqrt[ (4 Pi)^2 a2[0, mw,mz,mh,mt,smu]],
		  gs -> Sqrt[ (4 Pi) asmu],
		  yt -> Sqrt[ (4 Pi)^2 at[0, mw,mz,mh,mt,smu]],
		  lam -> (4 Pi)^2 alam[0, mw,mz,mh,mt,smu], 
		  mu -> smu} /. {aQCD[smu]->asmu/(4 Pi), Gf->gf, aEW[smu]->aew/(4 Pi)} /. solaew
		  ];

RunParsFromPoleMassesAndAsWithMb[mz_:PDG`MZ,mw_:PDG`MW,mh_:PDG`MH,mt_:PDG`MT,mb_:PDG`Mb,gf_:PDG`GF,asmu_,smu_] := Module[{aew,seq,aa1,aa2, impliciteq,res,solaew},
		aa1 = (4 Pi) 3/5 * a1[mb,mw,mz,mh,mt,smu] /.aQCD[smu]->asmu/(4 Pi) /. Gf->gf /. aEW[smu]->aew/(4 Pi) ;
		aa2 = (4 Pi)       a2[mb,mw,mz,mh,mt,smu] /.aQCD[smu]->asmu/(4 Pi) /. Gf->gf /. aEW[smu]->aew/(4 Pi) ;
		impliciteq = (aew == (Simplify[aa1*aa2/(aa1 + aa2)]));
		impliciteqcheck = (aew == Normal[Series[Simplify[aa1*aa2/(aa1 + aa2)], {aew,0,2}]]);
		(* solution for alpha ew at the scale *)
		solaew = FindRoot[ impliciteq, {aew, 1/127.94}]; 
		(* Couplings *)
		res = { gp->Sqrt[ (4 Pi)^2 3/5 a1[mb, mw,mz,mh,mt,smu]], 
		  g -> Sqrt[ (4 Pi)^2 a2[mb, mw,mz,mh,mt,smu]],
		  gs -> Sqrt[ (4 Pi) asmu],
		  yb -> Sqrt[ (4 Pi)^2 ab[mb, mw,mz,mh,mt,smu]],
		  yt -> Sqrt[ (4 Pi)^2 at[mb, mw,mz,mh,mt,smu]],
		  lam -> (4 Pi)^2 alam[mb,mw,mz,mh,mt,smu], 
		(* ms(mu) = mh(mu) ? *)
		  ms -> mh * Sqrt[mmHMMH[mb,mw,mz,mh,mt,smu]],
		  vev -> vev[mb,mw,mz,mh,mt,smu],
		  mu -> smu};
		(*		Print[res, ":::",2^(-1/2)*Gf/(4*Pi)^2*mh^2*(1+aEW[mu]*yH[1,0]+aEW[mu]*aQCD[mu]*yH[1,1]+aEW[mu]^2*yH[2,0]) /. mu->smu/. XH[mb,mw,mz,mh,mt,smu,1,2], "::",alam[mb, mw,mz,mh,mt,smu], ":::", XH[mb,mw,mz,mh,mt,smu,2,1]];*)
	       	res = res /. {aQCD[smu]->h as asmu/(4 Pi), Gf->gf, aEW[smu]->h aw aew/(4 Pi)} /. solaew;
		(*Print["1/a=", 1/aew /. solaew];
		Print["as=", asmu,":", res, "->",res /. h->1];*)
		res /. (a_->b_):>(a->(b/. h->0)*corr[Collect[1/(b/.h->0)*(Normal[Series[b,{h,0,2}]]),h,Expand],(b /. h->1 /. aw->1 /. as->1)])
		  ];

alpha[gp_,g_] := g^2 gp^2/(g^2 + gp^2)/(4 Pi)
alpha[runpars_List] := (g^2 gp^2/(g^2 + gp^2)/(4 Pi)) /. runpars

