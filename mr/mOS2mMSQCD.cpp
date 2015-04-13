//
// MR - 2-loop matching and 3-loop Running, including full 2-loop EW corrections
// Copyright (C) 2014 Andrey Pikelner <pikelner@theor.jinr.ru>
//
// This file is part of MR.
//
// MR is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MR is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MR.  If not, see <http://www.gnu.org/licenses/>.
//

#include <bb.hpp>
#include <tt.hpp>
#include <stdexcept>

// nh - number of heavy fermions with mass M
// nl - number of massles fermions
long double mOS2mMS(long double MM, long double mu2, size_t nl, size_t nh, size_t loops)
{
  double NL = double(nl), NH = double (nh);
  if(loops == 1)
    {
      long double mOS2mMS1 =   + log(mu2/MM) * (  - 4 );
      
      mOS2mMS1 +=  - 16./3.;

      return mOS2mMS1;
    }
  else if(loops == 2)
    {
      long double mOS2mMS2 =   + log(mu2/MM) * (  - 314./3. + 52./9.*NL + 52./9.*NH );

      mOS2mMS2 +=  + pow(log(mu2/MM),2) * (  - 14 + 4./3.*NL + 4./3.*NH );

      mOS2mMS2 +=  + Zeta2 * (  - 64./3. - 32./3.*log(2.) + 16./3.*NL - 32./3.*NH );

      mOS2mMS2 +=  - 3305./18. + 8./3.*Zeta3 + 71./9.*NL + 143./9.*NH;

      return mOS2mMS2;
    }
  else if(loops == 3)
    {
      long double mOS2mMS3 = + log(mu2/MM) * (  - 42650./9. + 48*Zeta3 + 14164./27.*NL + 448./9.*NL
                          *Zeta3 - 712./81.*pow(NL,2) + 18052./27.*NH + 448./9.*NH*Zeta3 - 
                          2288./81.*NH*NL - 1576./81.*pow(NH,2) );
      
      mOS2mMS3 +=  + pow(log(mu2/MM),2) * (  - 3034./3. + 428./3.*NL - 104./27.*pow(NL,2) 
                                   + 428./3.*NH - 208./27.*NH*NL - 104./27.*pow(NH,2) );
      
      mOS2mMS3 +=  + pow(log(mu2/MM),3) * (  - 84 + 128./9.*NL - 16./27.*pow(NL,2) + 
                                   128./9.*NH - 32./27.*NH*NL - 16./27.*pow(NH,2) );
      
      mOS2mMS3 +=  + Zeta2 * (  - 99980./27. - 4928./3.*log(2.) + 896./9.*pow(log(2.),2) 
                              + 10648./9.*Zeta3 + 15056./27.*NL + 1408./27.*NL*log(2.) - 
                              256./27.*NL*pow(log(2.),2) - 416./27.*pow(NL,2) - 215728./81.*NH
                              + 81920./27.*NH*log(2.) + 128./27.*NH*pow(log(2.),2) + 96*NH*Zeta3 + 
                              416./27.*NH*NL + 512./135.*pow(NH,2) );

      mOS2mMS3 +=  + Zeta2*log(mu2/MM) * (  - 384 - 192*log(2.) + 1120./9.*NL + 128./9.*NL*
                                log(2.) - 64./9.*pow(NL,2) - 1472./9.*NH + 128./9.*NH*log(2.) + 64./9.
                                *NH*NL + 128./9.*pow(NH,2) );

      mOS2mMS3 +=  - 1259285./162. + 4864./9.*a4 + 608./27.*pow(log(2.),4) - 
        13640./27.*Zeta5 + 6820./9.*Zeta4 + 584./9.*Zeta3 + 172318./243.*NL
          - 512./27.*NL*a4 - 64./81.*NL*pow(log(2.),4) - 4880./27.*NL*Zeta4
        + 5656./27.*NL*Zeta3 - 4706./729.*pow(NL,2) - 224./27.*pow(NL,2)*Zeta3 
        + 315526./243.*NH - 512./27.*NH*a4 - 64./81.*NH*pow(log(2.),4) 
        - 80*NH*Zeta5 - 6560./27.*NH*Zeta4 - 6008./27.*NH*Zeta3 - 
        23668./729.*NH*NL + 128./27.*NH*NL*Zeta3 - 18962./729.*pow(NH,2)
        + 352./27.*pow(NH,2)*Zeta3;
      
      return mOS2mMS3;
    }
  else if(loops == 4)
    {
      long double mOS2mMS4 = + log(mu2/MM) * (  - 839677./3. + 141056./9.*a4 + 17632./27.*pow(
         log(2.),4) - 157960./27.*Zeta5 + 197780./9.*Zeta4 - 75032./27.*Zeta3 + 
         11149742./243.*NL - 44032./27.*NL*a4 - 5504./81.*NL*pow(log(2.),4)
          - 27920./27.*NL*Zeta5 - 206200./27.*NL*Zeta4 + 84584./9.*NL*Zeta3
          - 1349140./729.*pow(NL,2) + 1024./27.*pow(NL,2)*a4 + 128./81.
         *pow(NL,2)*pow(log(2.),4) + 11200./27.*pow(NL,2)*Zeta4 - 6736./9.*
         pow(NL,2)*Zeta3 + 10408./729.*pow(NL,3) + 128./9.*pow(NL,3)*Zeta3
          + 15568454./243.*NH - 44032./27.*NH*a4 - 5504./81.*NH*pow(
         log(2.),4) - 90560./27.*NH*Zeta5 - 254920./27.*NH*Zeta4 - 28168./9.*NH
         *Zeta3 - 4105736./729.*NH*NL + 2048./27.*NH*NL*a4 + 256./81.*NH*
         NL*pow(log(2.),4) + 160*NH*NL*Zeta5 + 25760./27.*NH*NL*Zeta4 - 128./9.
         *NH*NL*Zeta3 + 19912./243.*NH*pow(NL,2) - 2756596./729.*pow(
         NH,2) + 1024./27.*pow(NH,2)*a4 + 128./81.*pow(NH,2)*pow(log(2.),4)
          + 160*pow(NH,2)*Zeta5 + 14560./27.*pow(NH,2)*Zeta4 + 6608./9.*
         pow(NH,2)*Zeta3 + 29416./243.*pow(NH,2)*NL - 128./3.*pow(NH,2)*
         NL*Zeta3 );

      mOS2mMS4 +=  + log(mu2/MM) * ( 38920./729.*pow(NH,3) - 256./9.*pow(NH,3)*Zeta3 );

      mOS2mMS4 +=  + pow(log(mu2/MM),2) * (  - 686026./9. + 696*Zeta3 + 387311./27.*NL
          + 5104./9.*NL*Zeta3 - 59479./81.*pow(NL,2) - 448./9.*pow(NL,2)*
         Zeta3 + 712./81.*pow(NL,3) + 443687./27.*NH + 5104./9.*NH*Zeta3 - 
         143150./81.*NH*NL - 896./9.*NH*NL*Zeta3 + 1000./27.*NH*pow(NL,2)
          - 83671./81.*pow(NH,2) - 448./9.*pow(NH,2)*Zeta3 + 1288./27.*
         pow(NH,2)*NL + 1576./81.*pow(NH,3) );

      mOS2mMS4 +=  + pow(log(mu2/MM),3) * (  - 10414 + 59992./27.*NL - 11552./81.*pow(
         NL,2) + 208./81.*pow(NL,3) + 59992./27.*NH - 23104./81.*NH*NL
          + 208./27.*NH*pow(NL,2) - 11552./81.*pow(NH,2) + 208./27.*
         pow(NH,2)*NL + 208./81.*pow(NH,3) );

      mOS2mMS4 +=  + pow(log(mu2/MM),4) * (  - 609 + 1306./9.*NL - 308./27.*pow(NL,2)
          + 8./27.*pow(NL,3) + 1306./9.*NH - 616./27.*NH*NL + 8./9.*NH*
         pow(NL,2) - 308./27.*pow(NH,2) + 8./9.*pow(NH,2)*NL + 8./27.*
         pow(NH,3) );

      mOS2mMS4 +=  + Zeta2*log(mu2/MM) * (  - 2978140./27. - 441856./9.*log(2.) + 25984./9.*
         pow(log(2.),2) + 308792./9.*Zeta3 + 223192./9.*NL + 136192./27.*NL*
         log(2.) - 12800./27.*NL*pow(log(2.),2) - 21296./9.*NL*Zeta3 - 5056./3.*
         pow(NL,2) - 2816./27.*pow(NL,2)*log(2.) + 512./27.*pow(NL,2)*pow(
         log(2.),2) + 832./27.*pow(NL,3) - 5734376./81.*NH + 91520*NH*log(2.)
          - 1664./27.*NH*pow(log(2.),2) + 3760./9.*NH*Zeta3 + 387296./81.*NH*
         NL - 55552./9.*NH*NL*log(2.) + 256./27.*NH*NL*pow(log(2.),2) - 192*NH*
         NL*Zeta3 + 2301664./405.*pow(NH,2) - 163840./27.*pow(NH,2)*log(2.)
          - 256./27.*pow(NH,2)*pow(log(2.),2) - 192*pow(NH,2)*Zeta3 - 192./5.
         *pow(NH,2)*NL - 1024./135.*pow(NH,3) );

      mOS2mMS4 +=  + Zeta2*pow(log(mu2/MM),2) * (  - 5568 - 2784*log(2.) + 19696./9.*NL + 
         3584./9.*NL*log(2.) - 2048./9.*pow(NL,2) - 128./9.*pow(NL,2)*log(2.)
          + 64./9.*pow(NL,3) - 17888./9.*NH + 3584./9.*NH*log(2.) + 1280./9.
         *NH*NL - 256./9.*NH*NL*log(2.) + 3328./9.*pow(NH,2) - 128./9.*pow(
         NH,2)*log(2.) - 64./3.*pow(NH,2)*NL - 128./9.*pow(NH,3) );

      // Constant part 
      // nf^3,nf^2 from arXiv:1301.6481 [hep-ph]
      // nf^1,nf^0 from arXiv:1502.01030 [hep-ph]


      const long double zm4_3 = -1744.8; // +/- 21.5
      const long double zm4_4 = -1267.0; // +/- 21.5
      const long double zm4_5 = -859.96; // +/- 21.5

      if (NH==1)
        if (NL == 3)
          mOS2mMS4 +=  + 256*(zm4_3);
        else if (NL == 4)
          mOS2mMS4 +=  + 256*(zm4_4);
        else if (NL == 5)
          mOS2mMS4 +=  + 256*(zm4_5);
        else throw std::logic_error("ERROR: relation between mOS and mMS available only for NL=3,4,5"); 
      else throw std::logic_error("ERROR: relation between mOS and mMS available only for NH=1"); 
      
      return mOS2mMS4;
    }
  else 
    return 0;
}



// Botom mass in nf=6 scheme with explicit dependence on top quark
// mass
//
//    M - Pole mass of b-quark
//    LmuM=Log(mu^2/M^2)
//
//    x=M_t/M
//
// 
long double mOS2mMSnm(long double MM, long double x, long double mu2, size_t nl_,size_t nm_, size_t nh_, size_t loops)
{
  double nl = double(nl_), nm = double (nm_), nh = double (nh_);
  
  long double LmuM = log(mu2/MM);
  
  
  if(loops == 1)
    {
      long double zm1l =
        (
         - 16./3.
         );
      
      zm1l +=
        + LmuM * (
                  - 4
                  );
      return zm1l;
    }
  else if(loops == 2)
    {
      long double zm2l =
        (
         - 3305./18.
         + 8./3.*Zeta3
         - 64./3.*Zeta2
         - 32./3.*Zeta2*log(2.)
         );
      
      zm2l +=  nm * (
                     + 71./9.
                     + 16./3.*Zeta2
                     - 16*x*Zeta2
                     + 8*pow(x,2)
                     - 16*pow(x,3)*Zeta2
                     + 16./3.*pow(x,4)*Zeta2
                     + 16./3.*log(x)*pow(x,2)
                     - 16./3.*log(x)*log(1 - x)
                     + 16./3.*log(x)*log(1 - x)*x
                     + 16./3.*log(x)*log(1 - x)*pow(x,3)
                     - 16./3.*log(x)*log(1 - x)*pow(x,4)
                     - 16./3.*log(x)*log(1 + x)
                     - 16./3.*log(x)*log(1 + x)*x
                     - 16./3.*log(x)*log(1 + x)*pow(x,3)
                     );
  
      zm2l +=  nm * (
                     - 16./3.*log(x)*log(1 + x)*pow(x,4)
                     + 16./3.*pow(log(x),2)*pow(x,4)
                     - 16./3.*Li2( - x).real()
                     - 16./3.*Li2( - x).real()*x
                     - 16./3.*Li2( - x).real()*pow(x,3)
                     - 16./3.*Li2( - x).real()*pow(x,4)
                     - 16./3.*Li2(x).real()
                     + 16./3.*Li2(x).real()*x
                     + 16./3.*Li2(x).real()*pow(x,3)
                     - 16./3.*Li2(x).real()*pow(x,4)
                     );
  
      zm2l +=  nl * (
                     + 71./9.
                     + 16./3.*Zeta2
                     );
      
      zm2l +=  nh * (
                     + 143./9.
                     - 32./3.*Zeta2
                     );
      
      zm2l +=  LmuM * (
                       - 314./3.
                       );
  
      zm2l +=  LmuM*nm * (
                          + 52./9.
                          );
  
      zm2l +=  LmuM*nl * (
                          + 52./9.
                          );
  
      zm2l +=  LmuM*nh * (
                          + 52./9.
                          );
  
      zm2l +=  pow(LmuM,2) * (
                              - 14
                              );
      
      zm2l +=  pow(LmuM,2)*nm * (
                                 + 4./3.
                                 );
      
      zm2l +=  pow(LmuM,2)*nl * (
                                 + 4./3.
                                 );
      
      zm2l +=  pow(LmuM,2)*nh * (
                                 + 4./3.
                                 );

      return zm2l;
    }
  else
    return 0;
}






// To convert number of light and heavy generations NL,NH to
// number of fermions we use:
//        nl = 2*NL + NH, nh = NH for M=Mt
//        nl = 2*NL , nh = NH for M=Mb 

long double bb<OS>::x02(size_t nL, size_t nH, size_t boson)
{     
  return mOS2mMS(MMb, mu2, 2*nL, nH, 2);
}

long double bb<OS>::x03(size_t nL, size_t nH, size_t boson)
{     
  return mOS2mMS(MMb, mu2, 2*nL, nH, 3);
}

long double bb<OS>::x04(size_t nL, size_t nH, size_t boson)
{     
  return mOS2mMS(MMb, mu2, 2*nL, nH, 4);
}

long double bb<OS>::y02(size_t nL, size_t nH, size_t boson)
{     
  return mOS2mMS(MMb, mu2, 2*nL, nH, 2);
}

long double bb<OS>::y03(size_t nL, size_t nH, size_t boson)
{     
  return mOS2mMS(MMb, mu2, 2*nL, nH, 3);
}

long double bb<OS>::y04(size_t nL, size_t nH, size_t boson)
{     
  return mOS2mMS(MMb, mu2, 2*nL, nH, 4);
}

long double tt<OS>::x02(size_t nL, size_t nH, size_t boson)
{     
  return mOS2mMS(MMt, mu2, 2*nL + nH, nH, 2);
}

long double tt<OS>::x03(size_t nL, size_t nH, size_t boson)
{     
  return mOS2mMS(MMt, mu2, 2*nL + nH, nH, 3);
}

long double tt<OS>::x04(size_t nL, size_t nH, size_t boson)
{     
  return mOS2mMS(MMt, mu2, 2*nL + nH, nH, 4);
}

long double tt<OS>::y02(size_t nL, size_t nH, size_t boson)
{     
  return mOS2mMS(MMt, mu2, 2*nL + nH, nH, 2);
}

long double tt<OS>::y03(size_t nL, size_t nH, size_t boson)
{     
  return mOS2mMS(MMt, mu2, 2*nL + nH, nH, 3);
}

long double tt<OS>::y04(size_t nL, size_t nH, size_t boson)
{     
  return mOS2mMS(MMt, mu2, 2*nL + nH, nH, 4);
}



// Inverted relations
// Pole mass M as function of running m(mu) and scale mu
long double mMS2mOS(long double mm, long double mu2, size_t nl, size_t nh, size_t loops)
{
  double NL = double(nl), NH = double (nh);
  if(loops == 1)
    {
      long double mMS2mOS1 =   + log(mu2/mm) * (  4 );
      
      mMS2mOS1 +=  16./3.;

      return mMS2mOS1;
    }
  else if(loops == 2)
    {
      long double mMS2mOS2 =   + log(mu2/mm) * ( 346./3. - 52./9.*NL - 52./9.*NH );
      
      mMS2mOS2 +=  + pow(log(mu2/mm),2) * ( 30 - 4./3.*NL - 4./3.*NH );
      
      mMS2mOS2 +=  + Zeta2 * ( 64./3. + 32./3.*log(2.) - 16./3.*NL + 32./3.*NH );
      
      mMS2mOS2 +=  + 3049./18. - 8./3.*Zeta3 - 71./9.*NL - 143./9.*NH;
        
        
      return mMS2mOS2;
    }
  else if(loops == 3)
    {
      long double mMS2mOS3 =   + log(mu2/mm) * ( 45854./9. - 208./3.*Zeta3 - 4756./9.*NL - 448./9.*NL
                                       *Zeta3 + 712./81.*pow(NL,2) - 6628./9.*NH - 448./9.*NH*Zeta3 + 
                                       2288./81.*NH*NL + 1576./81.*pow(NH,2) );
      
      mMS2mOS3 +=  + pow(log(mu2/mm),2) * ( 1598 - 1540./9.*NL + 104./27.*pow(NL,2) - 
                                  1540./9.*NH + 208./27.*NH*NL + 104./27.*pow(NH,2) );
      
      mMS2mOS3 +=  + pow(log(mu2/mm),3) * ( 260 - 224./9.*NL + 16./27.*pow(NL,2) - 224.
                                  /9.*NH + 32./27.*NH*NL + 16./27.*pow(NH,2) );
      
      mMS2mOS3 +=  + Zeta2 * ( 101516./27. + 15040./9.*log(2.) - 896./9.*pow(log(2.),2)
                               - 10648./9.*Zeta3 - 15440./27.*NL - 1408./27.*NL*log(2.) + 
                               256./27.*NL*pow(log(2.),2) + 416./27.*pow(NL,2) + 218032./81.*NH
                               - 81920./27.*NH*log(2.) - 128./27.*NH*pow(log(2.),2) - 96*NH*Zeta3 - 
                               416./27.*NH*NL - 512./135.*pow(NH,2) );
      
      mMS2mOS3 +=  + Zeta2*log(mu2/mm) * ( 1664./3. + 832./3.*log(2.) - 1504./9.*NL - 128./9.
                                 *NL*log(2.) + 64./9.*pow(NL,2) + 2240./9.*NH - 128./9.*NH*log(2.) - 64.
                                 /9.*NH*NL - 128./9.*pow(NH,2) );
      
      mMS2mOS3 +=  + 1145453./162. - 4864./9.*a4 - 608./27.*pow(log(2.),4) + 
        13640./27.*Zeta5 - 6820./9.*Zeta4 - 72*Zeta3 - 162454./243.*NL + 512.
        /27.*NL*a4 + 64./81.*NL*pow(log(2.),4) + 4880./27.*NL*Zeta4 - 5656./
        27.*NL*Zeta3 + 4706./729.*pow(NL,2) + 224./27.*pow(NL,2)*Zeta3 - 
        310846./243.*NH + 512./27.*NH*a4 + 64./81.*NH*pow(log(2.),4) + 80*
        NH*Zeta5 + 6560./27.*NH*Zeta4 + 6008./27.*NH*Zeta3 + 23668./729.*NH*
        NL - 128./27.*NH*NL*Zeta3 + 18962./729.*pow(NH,2) - 352./27.*
        pow(NH,2)*Zeta3;
        
      return mMS2mOS3;
    }
  else if(loops == 4)
    {
      long double mMS2mOS4 = + log(mu2/mm) * ( 23003021./81. - 179968./9.*a4 - 22496./27.*pow(
         log(2.),4) + 267080./27.*Zeta5 - 252340./9.*Zeta4 + 59192./27.*Zeta3 - 
         10895302./243.*NL + 48128./27.*NL*a4 + 6016./81.*NL*pow(log(2.),4)
          + 27920./27.*NL*Zeta5 + 245240./27.*NL*Zeta4 - 10792*NL*Zeta3 + 
         1279820./729.*pow(NL,2) - 1024./27.*pow(NL,2)*a4 - 128./81.*
         pow(NL,2)*pow(log(2.),4) - 11200./27.*pow(NL,2)*Zeta4 + 22000./27.*
         pow(NL,2)*Zeta3 - 10408./729.*pow(NL,3) - 128./9.*pow(NL,3)*Zeta3
          - 16508926./243.*NH + 48128./27.*NH*a4 + 6016./81.*NH*pow(
         log(2.),4) + 107840./27.*NH*Zeta5 + 307400./27.*NH*Zeta4 + 5192*NH*Zeta3
          + 4075960./729.*NH*NL - 2048./27.*NH*NL*a4 - 256./81.*NH*NL*
         pow(log(2.),4) - 160*NH*NL*Zeta5 - 25760./27.*NH*NL*Zeta4 - 640./27.*
         NH*NL*Zeta3 - 19912./243.*NH*pow(NL,2) + 2796140./729.*pow(NH,2)
          - 1024./27.*pow(NH,2)*a4 - 128./81.*pow(NH,2)*pow(log(2.),4) - 
         160*pow(NH,2)*Zeta5 - 14560./27.*pow(NH,2)*Zeta4 - 22640./27.*pow(
         NH,2)*Zeta3 - 29416./243.*pow(NH,2)*NL + 128./3.*pow(NH,2)*NL*
         Zeta3 );

      mMS2mOS4 +=  + log(mu2/mm) * (  - 38920./729.*pow(NH,3) + 256./9.*pow(NH,3)*Zeta3
          );

      mMS2mOS4 +=  + pow(log(mu2/mm),2) * ( 295300./3. - 3848./3.*Zeta3 - 439711./27.*NL
          - 8624./9.*NL*Zeta3 + 60143./81.*pow(NL,2) + 448./9.*pow(NL,2)*
         Zeta3 - 712./81.*pow(NL,3) - 543607./27.*NH - 8624./9.*NH*Zeta3 + 
         153118./81.*NH*NL + 896./9.*NH*NL*Zeta3 - 1000./27.*NH*pow(NL,2)
          + 92975./81.*pow(NH,2) + 448./9.*pow(NH,2)*Zeta3 - 1288./27.*
         pow(NH,2)*NL - 1576./81.*pow(NH,3) );

      mMS2mOS4 +=  + pow(log(mu2/mm),3) * ( 20342 - 91064./27.*NL + 13696./81.*pow(
         NL,2) - 208./81.*pow(NL,3) - 91064./27.*NH + 27392./81.*NH*NL
          - 208./27.*NH*pow(NL,2) + 13696./81.*pow(NH,2) - 208./27.*
         pow(NH,2)*NL - 208./81.*pow(NH,3) );

      mMS2mOS4 +=  + pow(log(mu2/mm),4) * ( 2405 - 3242./9.*NL + 484./27.*pow(NL,2) - 
         8./27.*pow(NL,3) - 3242./9.*NH + 968./27.*NH*NL - 8./9.*NH*
         pow(NL,2) + 484./27.*pow(NH,2) - 8./9.*pow(NH,2)*NL - 8./27.*
         pow(NH,3) );

      mMS2mOS4 +=  + Zeta2 * ( 630304./81. + 90016./27.*log(2.) - 7168./27.*pow(
         log(2.),2) - 88256./27.*Zeta3 - 512./9.*Zeta3*log(2.) - 60784./81.*NL - 
         2624./81.*NL*log(2.) + 2048./81.*NL*pow(log(2.),2) + 256./9.*NL*Zeta3 - 
         992./81.*pow(NL,2) + 1440416./243.*NH - 660544./81.*NH*log(2.) - 
         1024./81.*NH*pow(log(2.),2) - 2816./9.*NH*Zeta3 + 7904./81.*NH*NL - 
         30016./405.*pow(NH,2) );

      mMS2mOS4 +=  + Zeta2*log(mu2/mm) * ( 3792572./27. + 187520./3.*log(2.) - 33152./9.*
         pow(log(2.),2) - 393976./9.*Zeta3 - 793160./27.*NL - 147200./27.*NL*
         log(2.) + 14848./27.*NL*pow(log(2.),2) + 21296./9.*NL*Zeta3 + 48704./27.
         *pow(NL,2) + 2816./27.*pow(NL,2)*log(2.) - 512./27.*pow(NL,2)*pow(
         log(2.),2) - 832./27.*pow(NL,3) + 7483624./81.*NH - 1042048./9.*NH
         *log(2.) + 640./27.*NH*pow(log(2.),2) - 10672./9.*NH*Zeta3 - 396896./81.
         *NH*NL + 55552./9.*NH*NL*log(2.) - 256./27.*NH*NL*pow(log(2.),2) + 192
         *NH*NL*Zeta3 - 2310112./405.*pow(NH,2) + 163840./27.*pow(NH,2)*
         log(2.) + 256./27.*pow(NH,2)*pow(log(2.),2) + 192*pow(NH,2)*Zeta3 + 192./
         5.*pow(NH,2)*NL + 1024./135.*pow(NH,3) );

      mMS2mOS4 +=  + Zeta2*pow(log(mu2/mm),2) * ( 30784./3. + 15392./3.*log(2.) - 32816./9.*
         NL - 4864./9.*NL*log(2.) + 896./3.*pow(NL,2) + 128./9.*pow(NL,2)*
         log(2.) - 64./9.*pow(NL,3) + 36448./9.*NH - 4864./9.*NH*log(2.) - 640./
         3.*NH*NL + 256./9.*NH*NL*log(2.) - 512*pow(NH,2) + 128./9.*pow(
         NH,2)*log(2.) + 64./3.*pow(NH,2)*NL + 128./9.*pow(NH,3) );

      mMS2mOS4 +=  + pow(Zeta2,2) * ( 4096./9. + 4096./9.*log(2.) + 1024./9.*pow(
         log(2.),2) - 2048./9.*NL - 1024./9.*NL*log(2.) + 256./9.*pow(NL,2) + 
         4096./9.*NH + 2048./9.*NH*log(2.) - 1024./9.*NH*NL + 1024./9.*pow(
         NH,2) );

      mMS2mOS4 +=  - 39571813./972. - 38912./27.*a4 - 4864./81.*
         pow(log(2.),4) + 109120./81.*Zeta5 - 54560./27.*Zeta4 + 2392./27.*Zeta3
          + 64./9.*pow(Zeta3,2) + 3477073./729.*NL + 4096./81.*NL*a4 + 
         512./243.*NL*pow(log(2.),4) + 39040./81.*NL*Zeta4 - 1328./81.*NL*Zeta3
          - 230669./2187.*pow(NL,2) + 1792./81.*pow(NL,2)*Zeta3 + 2903593.
         /729.*NH + 4096./81.*NH*a4 + 512./243.*NH*pow(log(2.),4) + 640./3.
         *NH*Zeta5 + 52480./81.*NH*Zeta4 + 95440./81.*NH*Zeta3 - 522250./2187.
         *NH*NL - 1024./81.*NH*NL*Zeta3 - 151613./2187.*pow(NH,2) - 2816./
        81.*pow(NH,2)*Zeta3;

      

      const long double cm4_3 = 1691.2; // +/- 21.5
      const long double cm4_4 = 1224.0; // +/- 21.5
      const long double cm4_5 = 827.37; // +/- 21.5

      const long double zm4_3 = -1744.8; // +/- 21.5
      const long double zm4_4 = -1267.0; // +/- 21.5
      const long double zm4_5 = -859.96; // +/- 21.5

      if (NH==1)
        if (NL == 3)
          mMS2mOS4 +=  - 256*(zm4_3);
        else if (NL == 4)
          mMS2mOS4 +=  - 256*(zm4_4);
        else if (NL == 5)
          mMS2mOS4 +=  - 256*(zm4_5);
        else throw std::logic_error("ERROR: relation between mMS and mOS available only for NL=3,4,5"); 
      else throw std::logic_error("ERROR: relation between mMS and mOS available only for NH=1"); 
      
      return mMS2mOS4;
    }
  else 
    return 0;
}

// bottom
long double bb<MS>::x02(size_t nL, size_t nH, size_t boson)
{     
  return mMS2mOS(mmb, mu2, 2*nL, nH, 2);
}

long double bb<MS>::x03(size_t nL, size_t nH, size_t boson)
{     
  return mMS2mOS(mmb, mu2, 2*nL, nH, 3);
}

long double bb<MS>::x04(size_t nL, size_t nH, size_t boson)
{     
  return mMS2mOS(mmb, mu2, 2*nL, nH, 4);
}

// top
long double tt<MS>::x02(size_t nL, size_t nH, size_t boson)
{     
  return mMS2mOS(mmt, mu2, 2*nL + nH, nH, 2);
}

long double tt<MS>::x03(size_t nL, size_t nH, size_t boson)
{     
  return mMS2mOS(mmt, mu2, 2*nL + nH, nH, 3);
}

long double tt<MS>::x04(size_t nL, size_t nH, size_t boson)
{     
  return mMS2mOS(mmt, mu2, 2*nL + nH, nH, 4);
}

