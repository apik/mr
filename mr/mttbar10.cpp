#include <tt.hpp>
long double tt::x10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armttbar[25], mttbarret;

    armttbar[1]=double(nH);
    armttbar[2]=double(boson);
    armttbar[3]=Tsil::A(MMt,mu2);
    armttbar[4]=pow(CW,-1);
    armttbar[5]=pow(MMH,-1);
    armttbar[6]=pow(MMZ,-1);
    armttbar[7]=pow(SW,-1);
    armttbar[8]=Tsil::A(MMb,mu2);
    armttbar[9]=Tsil::B(MMH,MMt,MMt,mu2);
    armttbar[10]=Tsil::B(MMZ,MMt,MMt,mu2);
    armttbar[11]=pow(MMt,-1);
    armttbar[12]=Tsil::B(MMW,MMb,MMt,mu2);
    armttbar[13]=Tsil::A(MMH,mu2);
    armttbar[14]=Tsil::A(MMZ,mu2);
    armttbar[15]=Tsil::A(MMW,mu2);
   armttbar[16]=1./2.*armttbar[14];
   armttbar[17]=armttbar[10]*MMZ;
   armttbar[18]= - armttbar[12]*MMZ;
   armttbar[18]= - 1./2.*armttbar[17] + armttbar[18] + armttbar[8] - 
   armttbar[16];
   armttbar[19]=1./2.*armttbar[11];
   armttbar[18]=armttbar[18]*armttbar[19];
   armttbar[20]=3*armttbar[5];
   armttbar[16]=MMZ + armttbar[16];
   armttbar[16]=armttbar[16]*armttbar[20];
   armttbar[19]=armttbar[20] - armttbar[19];
   armttbar[19]=armttbar[15]*armttbar[19];
   armttbar[21]=armttbar[3]*armttbar[11];
   armttbar[22]=armttbar[21] + armttbar[10] - 3 + armttbar[12];
   armttbar[23]=1./4.*armttbar[12];
   armttbar[24]=armttbar[23]*MMb*armttbar[11];
   armttbar[16]=armttbar[24] + armttbar[19] + armttbar[18] + 
   armttbar[16] + 1./4.*armttbar[22];
   armttbar[18]=MMb*armttbar[12];
   armttbar[18]=armttbar[8] + armttbar[18] - armttbar[15];
   armttbar[19]=1./8.*armttbar[11];
   armttbar[18]=armttbar[18]*armttbar[19];
   armttbar[19]=armttbar[20]*armttbar[1];
   armttbar[20]=armttbar[19]*armttbar[8];
   armttbar[18]=armttbar[18] - armttbar[20] - armttbar[23];
   armttbar[18]=armttbar[18]*MMb;
   armttbar[20]=MMH*armttbar[9];
   armttbar[20]=armttbar[20] - armttbar[8];
   armttbar[22]=armttbar[23] + armttbar[9];
   armttbar[22]=armttbar[22]*MMt;
   armttbar[20]=armttbar[22] - 1./4.*armttbar[20];
   armttbar[19]=armttbar[19]*MMt;
   armttbar[19]=armttbar[19] - 1./4.;
   armttbar[19]=armttbar[19]*armttbar[3];
   armttbar[18]= - 1./4.*armttbar[13] - armttbar[18] - 1./2.*
   armttbar[20] + armttbar[19] - 1./8.*armttbar[15];
   armttbar[18]=armttbar[6]*armttbar[18];
   armttbar[16]=1./2.*armttbar[16] - armttbar[18];
   armttbar[16]=armttbar[16]*pow(armttbar[7],2);
   armttbar[17]=armttbar[17] + armttbar[14];
   armttbar[19]= - armttbar[11]*armttbar[17];
   armttbar[19]=armttbar[21] + armttbar[19];
   armttbar[20]= - 1 - 7*armttbar[10];
   armttbar[22]=MMZ + 3./2.*armttbar[14];
   armttbar[22]=armttbar[5]*armttbar[22];
   armttbar[19]=1./36.*armttbar[20] + armttbar[22] + 17./36.*
   armttbar[19];
   armttbar[18]=1./2.*armttbar[19] - armttbar[18];
   armttbar[18]=armttbar[18]*pow(armttbar[4],2);
   armttbar[19]=MMZ*armttbar[23];
   armttbar[17]=armttbar[19] + 4./9.*armttbar[17];
   armttbar[17]=armttbar[11]*armttbar[17];
   armttbar[19]=armttbar[21] - 1 + armttbar[10];
   armttbar[20]= - armttbar[5]*MMZ;
   armttbar[16]=armttbar[18] + armttbar[16] + armttbar[17] + 
   armttbar[20] + 8./9.*armttbar[19];

      mttbarret = armttbar[16]*armttbar[2];
      return mttbarret.real();
}
