#include <bb.hpp>
std::complex<long double>
bb::m10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armbbbar[20], mbbbarret;

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
   armbbbar[13]=armbbbar[12]*armbbbar[8];
   armbbbar[13]= - 1./2. + armbbbar[13];
   armbbbar[14]=armbbbar[8] + 1./2.*armbbbar[7];
   armbbbar[14]=armbbbar[3]*armbbbar[14];
   armbbbar[15]=armbbbar[8] - armbbbar[9];
   armbbbar[15]=armbbbar[12]*armbbbar[15];
   armbbbar[15]=1 + armbbbar[15];
   armbbbar[15]=armbbbar[12]*armbbbar[15];
   armbbbar[15]=1./4.*armbbbar[15] + armbbbar[3];
   armbbbar[15]=MMZ*armbbbar[15];
   armbbbar[13]=armbbbar[15] + 1./4.*armbbbar[13] + armbbbar[14];
   armbbbar[14]=3*armbbbar[6] + 1./2.*MMt;
   armbbbar[15]=armbbbar[14] + 3*armbbbar[9];
   armbbbar[16]= - 3*armbbbar[3]*armbbbar[9]*MMt;
   armbbbar[15]=1./8.*armbbbar[15] + armbbbar[16];
   armbbbar[15]=armbbbar[4]*armbbbar[15];
   armbbbar[13]=3./2.*armbbbar[13] + armbbbar[15];
   armbbbar[13]=armbbbar[13]*pow(armbbbar[5],2);
   armbbbar[15]=armbbbar[3]*armbbbar[7];
   armbbbar[15]= - 31./36. + 3*armbbbar[15];
   armbbbar[17]=pow(armbbbar[2],2);
   armbbbar[15]=armbbbar[17]*armbbbar[15];
   armbbbar[18]= - armbbbar[8] + armbbbar[9];
   armbbbar[18]=armbbbar[12]*armbbbar[18];
   armbbbar[18]= - 1 + armbbbar[18];
   armbbbar[18]=armbbbar[12]*armbbbar[18];
   armbbbar[19]=armbbbar[17]*armbbbar[3];
   armbbbar[18]=1./2.*armbbbar[19] + 3./8.*armbbbar[18] - armbbbar[3];
   armbbbar[18]=MMZ*armbbbar[18];
   armbbbar[19]= - 1./3.*armbbbar[7];
   armbbbar[14]=3./4.*armbbbar[9] + 1./4.*armbbbar[14] + armbbbar[19];
   armbbbar[14]=1./2.*armbbbar[14] + armbbbar[16];
   armbbbar[14]=armbbbar[17]*armbbbar[14];
   armbbbar[14]=armbbbar[19] + armbbbar[14];
   armbbbar[14]=armbbbar[4]*armbbbar[14];
   armbbbar[16]=armbbbar[10]*armbbbar[11];
   armbbbar[16]= - 1./2. + armbbbar[16];
   armbbbar[13]=armbbbar[13] + armbbbar[14] + armbbbar[18] + 1./3.*
   armbbbar[16] + 1./4.*armbbbar[15];

      mbbbarret = armbbbar[13]*armbbbar[1];
      return mbbbarret;
}
