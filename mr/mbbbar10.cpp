#include <bb.hpp>
std::complex<long double>
bb::m10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armbbbar[21], mbbbarret;

    armbbbar[1]=double(boson);
    armbbbar[2]=pow(CW,-1);
    armbbbar[3]=pow(MMH,-1);
    armbbbar[4]=pow(MMZ,-1);
    armbbbar[5]=pow(SW,-1);
    armbbbar[6]=Tsil::A(MMH,mu2);
    armbbbar[7]=Tsil::A(MMZ,mu2);
    armbbbar[8]=Tsil::A(MMW,mu2);
    armbbbar[9]=Tsil::A(MMt,mu2);
    armbbbar[10]=Tsil::A(MMb,mu2);
    armbbbar[11]=pow(MMb,-1);
    armbbbar[12]=1/(MMt - MMW);
   armbbbar[13]=MMZ + armbbbar[8];
   armbbbar[14]=pow(armbbbar[5],2);
   armbbbar[13]=armbbbar[13]*armbbbar[14];
   armbbbar[15]=pow(armbbbar[2],2);
   armbbbar[16]=armbbbar[15] + armbbbar[14];
   armbbbar[17]=armbbbar[7]*armbbbar[16];
   armbbbar[18]=1./2.*armbbbar[15];
   armbbbar[19]= - 1 + armbbbar[18];
   armbbbar[19]=MMZ*armbbbar[19];
   armbbbar[17]=3./4.*armbbbar[17] + 3./2.*armbbbar[13] + armbbbar[19];
   armbbbar[17]=armbbbar[3]*armbbbar[17];
   armbbbar[19]=3*armbbbar[6] + 1./2.*MMt;
   armbbbar[20]= - armbbbar[3]*MMt;
   armbbbar[20]=1./8. + armbbbar[20];
   armbbbar[20]=armbbbar[9]*armbbbar[20];
   armbbbar[19]=1./8.*armbbbar[19] + 3*armbbbar[20];
   armbbbar[16]=armbbbar[16]*armbbbar[19];
   armbbbar[18]= - 1 - armbbbar[18];
   armbbbar[18]=armbbbar[7]*armbbbar[18];
   armbbbar[16]=1./3.*armbbbar[18] + armbbbar[16];
   armbbbar[16]=armbbbar[4]*armbbbar[16];
   armbbbar[18]= - armbbbar[9] + armbbbar[8];
   armbbbar[19]=armbbbar[14] - 1;
   armbbbar[18]=armbbbar[12]*armbbbar[18]*armbbbar[19]*MMZ;
   armbbbar[13]=armbbbar[18] - MMZ + armbbbar[13];
   armbbbar[13]=armbbbar[12]*armbbbar[13];
   armbbbar[18]=armbbbar[10]*armbbbar[11];
   armbbbar[15]= - 31./48.*armbbbar[15] - 1./2. + armbbbar[18];
   armbbbar[13]=3./8.*armbbbar[13] + armbbbar[16] + 1./3.*armbbbar[15]
    - 3./16.*armbbbar[14] + armbbbar[17];

      mbbbarret = armbbbar[13]*armbbbar[1];
      return mbbbarret;
}
