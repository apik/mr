#include <tt.hpp>
std::complex<long double>
tt::m10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armttbar[24], mttbarret;

    armttbar[1]=double(nH);
    armttbar[2]=Tsil::A(MMt,mu2);
    armttbar[3]=pow(CW,-1);
    armttbar[4]=pow(MMH,-1);
    armttbar[5]=pow(MMZ,-1);
    armttbar[6]=pow(SW,-1);
    armttbar[7]=double(boson);
    armttbar[8]=Tsil::B(MMH,MMt,MMt,mu2);
    armttbar[9]=Tsil::B(MMZ,MMt,MMt,mu2);
    armttbar[10]=pow(MMt,-1);
    armttbar[11]=Tsil::A(MMH,mu2);
    armttbar[12]=Tsil::A(MMZ,mu2);
    armttbar[13]=Tsil::A(MMW,mu2);
    armttbar[14]=std::real(Tsil::B(0,MMW,MMt,mu2));
   armttbar[15]=pow(armttbar[6],2);
   armttbar[16]=pow(armttbar[3],2);
   armttbar[17]=armttbar[15] + armttbar[16];
   armttbar[18]=MMt*armttbar[14];
   armttbar[18]=armttbar[18] + armttbar[13];
   armttbar[18]=armttbar[2] + armttbar[11] + 1./2.*armttbar[18];
   armttbar[19]= - MMt + 1./4.*MMH;
   armttbar[19]=armttbar[19]*armttbar[8];
   armttbar[18]= - armttbar[19] + 1./2.*armttbar[18];
   armttbar[18]=armttbar[18]*armttbar[17]*armttbar[5];
   armttbar[19]=armttbar[14]*MMZ;
   armttbar[20]= - armttbar[19] - armttbar[13] + 1./2.*armttbar[2];
   armttbar[20]=armttbar[20]*armttbar[15];
   armttbar[19]=armttbar[19] + armttbar[20];
   armttbar[20]=17./72.*armttbar[16];
   armttbar[21]=1./8.*armttbar[15];
   armttbar[22]=armttbar[21] + armttbar[20] - 4./9.;
   armttbar[23]= - armttbar[9]*MMZ*armttbar[22];
   armttbar[20]=armttbar[20] + 8./9.;
   armttbar[20]=armttbar[2]*armttbar[20];
   armttbar[19]=armttbar[23] + armttbar[20] + 1./4.*armttbar[19];
   armttbar[19]=armttbar[10]*armttbar[19];
   armttbar[17]=armttbar[17]*armttbar[4];
   armttbar[20]= - armttbar[10]*armttbar[22];
   armttbar[20]=3./4.*armttbar[17] + armttbar[20];
   armttbar[20]=armttbar[12]*armttbar[20];
   armttbar[22]= - 3 + armttbar[14];
   armttbar[22]=armttbar[22]*armttbar[21];
   armttbar[23]=armttbar[13] + MMZ;
   armttbar[15]=armttbar[23]*armttbar[15];
   armttbar[23]=1./2.*armttbar[16] - 1;
   armttbar[23]=MMZ*armttbar[23];
   armttbar[15]=3./2.*armttbar[15] + armttbar[23];
   armttbar[15]=armttbar[4]*armttbar[15];
   armttbar[21]= - 7./72.*armttbar[16] + 8./9. + armttbar[21];
   armttbar[21]=armttbar[9]*armttbar[21];
   armttbar[15]=armttbar[20] + armttbar[19] + armttbar[21] + 
   armttbar[15] + 1./2.*armttbar[18] - 1./72.*armttbar[16] - 8./9. + 
   armttbar[22];
   armttbar[15]=armttbar[7]*armttbar[15];
   armttbar[16]= - armttbar[5]*armttbar[1]*MMt*armttbar[2]*armttbar[17]
   ;

      mttbarret = armttbar[15] + 3*armttbar[16];
      return mttbarret;
}
