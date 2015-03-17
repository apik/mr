link=Install["mr"];
(*
MW[0.35761, 0.64822, 1.1666, 0.93558, 0.12711, 132.03, 173.2,2]
MW[0.35761, 0.64822, 1.1666, 0.93558, 0.12711, 132.03, 173.2,0]
MW[0.35761, 0.64822, 1.1666, 0.93558, 0.12711, 132.03, 173.2,1]



MZ[0.35761, 0.64822, 1.1666, 0.93558, 0.12711, 132.03, 173.2,2]
MZ[0.35761, 0.64822, 1.1666, 0.93558, 0.12711, 132.03, 173.2,0]
MZ[0.35761, 0.64822, 1.1666, 0.93558, 0.12711, 132.03, 173.2,1]


MH[0.35761, 0.64822, 1.1666, 0.93558, 0.12711, 132.03, 173.2,2]
MH[0.35761, 0.64822, 1.1666, 0.93558, 0.12711, 132.03, 173.2,0]
MH[0.35761, 0.64822, 1.1666, 0.93558, 0.12711, 132.03, 173.2,1]


MT[0.35761, 0.64822, 1.1666, 0.93558, 0.12711, 132.03, 173.2,2]
MT[0.35761, 0.64822, 1.1666, 0.93558, 0.12711, 132.03, 173.2,0]
MT[0.35761, 0.64822, 1.1666, 0.93558, 0.12711, 132.03, 173.2,1]


GF[0.35761, 0.64822, 1.1666, 0.93558, 0.12711, 132.03, 173.2,2]
GF[0.35761, 0.64822, 1.1666, 0.93558, 0.12711, 132.03, 173.2,0]
GF[0.35761, 0.64822, 1.1666, 0.93558, 0.12711, 132.03, 173.2,1]
*)


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
(*
sol = FindRoot[ {
	   MZ[gp, g, gs, yt, lam, mu0, 173.21,2] == PDG`MZ (* 91.1876 *),
	   MW[gp, g, gs, yt, lam, mu0, 173.21,2] == PDG`MW (* 80.385 *),
	   MT[gp, g, gs, yt, lam, mu0, 173.21,3] == PDG`MT (* 173.21 *),
	   MH[gp, g, gs, yt, lam, mu0, 173.21,2] == PDG`MH (* 125.7 *),
	   GF[gp, g, gs, yt, lam, mu0, 173.21,2] == PDG`GF (* 0.000011663787 *),
	   gs^2/(4*Pi)==0.108057
	  },{
		{gp, 0.357561},
		{g , 0.64822},
		{gs , 1.1666},
		{yt,0.93558},
		{lam,0.12711},
		{mu0,132.03}
	  }
	  ]
*)
sol  = Get["sol_mt,m"];

MW[gp, g, gs, yt, lam, mu0, 173.21,2] /. sol 
tmw = MWp[gp, g, gs, yt, lam, mu0, 173.21,2] /. sol 

tmwEv = tmw /. aew->1 /. as->1

(tmwEv[[1]] - tmwEv[[2]])/PDG`dMW 

tmwExp = Normal[Series[First[tmw] /. aew-> h /. as->h, {h,0,2}]] /. h->1

Print["MZ"]

MZ[gp, g, gs, yt, lam, mu0, 173.21,2] /. sol 
tMZ = MZp[gp, g, gs, yt, lam, mu0, 173.21,2] /. sol 

tMZEv = tMZ /. aew->1 /. as->1

(tMZEv[[1]] - tMZEv[[2]])/PDG`dMZ 

tMZExp = Normal[Series[First[tMZ] /. aew-> h /. as->h, {h,0,2}]] /. h->1

Print["MH"]

MH[gp, g, gs, yt, lam, mu0, 173.21,2] /. sol 
tMH = MHp[gp, g, gs, yt, lam, mu0, 173.21,2] /. sol 

tMHEv = tMH /. aew->1 /. as->1

(tMHEv[[1]] - tMHEv[[2]])/PDG`dMH 

tMHExp = Normal[Series[First[tMH] /. aew-> h /. as->h, {h,0,2}]] /. h->1

Print["MT"]

MT[gp, g, gs, yt, lam, mu0, 173.21,2] /. sol 

tMT = MTp[gp, g, gs, yt, lam, mu0, 173.21,2] /. sol 

tMTEv = tMT /. aew->1 /. as->1

(tMTEv - PDG`MT)/PDG`dMT


(*
{ MZ[gp, g, gs, yt, lam, mu0, 173.21,2], 
  MH[gp, g, gs, yt, lam, mu0, 173.21,2],  
  MT[gp, g, gs, yt, lam, mu0, 173.21,3],  
  GF[gp, g, gs, yt, lam, mu0, 173.21,2],  
  gs^2/(4 Pi)} /. sol

*)
(*
  FindMinimum[(MZ[gp, g, gs, yt, lam, mu0, 173.21,2] - 91.1876)^2/0.0021^2
 +(MW[gp, g, gs, yt, lam, mu0, 173.21,2] - 80.385)^2/0.015^2  
 +(MH[gp, g, gs, yt, lam, mu0, 173.21,2] - 125.7)^2/0.4^2  
 +(MT[gp, g, gs, yt, lam, mu0, 173.21,2] - 173.21)^2/1.22^2  
 +(GF[gp, g, gs, yt, lam, mu0, 173.21,2] - 0.000011663787)^2/(0.000000000006)^2  
 +(gs^2/(4 Pi) - 0.1080)^2/0.0005^2, 
{
		{gp, 0.357561},
		{g , 0.64822},
		{gs , 1.1666},
		{yt,0.93558},
		{lam,0.12711},
		{mu0,132.03}
	  },
Method->"LevenbergMarquardt"
 (* Method->"Newton" *)
	  ]
*)
(*
a1[0,PDG`MW,PDG`MZ,PDG`MH,PDG`MT,PDG`MT] /. aEW[PDG`MT]->ae (gp*g)^2/(gp^2+g^2)/4/Pi /. aQCD[PDG`MT]-> as gs^2/(4*Pi) /. Gf->PDG`GF /. sol /. ae->1
a2[0,PDG`MW,PDG`MZ,PDG`MH,PDG`MT,PDG`MT] /. aEW[PDG`MT]->ae (gp*g)^2/(gp^2+g^2)/4/Pi /. aQCD[PDG`MT]-> as gs^2/(4*Pi) /. Gf->PDG`GF /. sol /. ae->1

*)
(*
e^2/(4pi) = g1^2/(4 pi) g2^2/(4 pi) / (g1^2/(4pi) + g2^2/(4 pi))
*)
aa1 = (4 Pi) 3/5 * a1[0,PDG`MW,PDG`MZ,PDG`MH,PDG`MT,PDG`MT];
aa2 = (4 Pi) a2[0,PDG`MW,PDG`MZ,PDG`MH,PDG`MT,PDG`MT];

impliciteq1 = (aew == (Simplify[aa1*aa2/(aa1 + aa2)] /. aQCD[PDG`MT]->0.108057  /. Gf->PDG`GF /. aEW[PDG`MT]->aew))
impliciteq = (aew == Normal[Series[(Simplify[aa1*aa2/(aa1 + aa2)] /. aQCD[PDG`MT]->0.108057  /. Gf->PDG`GF /. aEW[PDG`MT]->aew), {aew,0,2}]])

{FindRoot[ impliciteq1, {aew, 1/127.94}],1/127.94} 
{FindRoot[ impliciteq, {aew, 1/127.94}],1/127.94} 

