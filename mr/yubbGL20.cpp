#include <bb.hpp>
long double bb::ygl20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryubbGL[29], yubbGLret;

    aryubbGL[1]=double(boson);
    aryubbGL[2]=pow(SW,-1);
    aryubbGL[3]=pow(MMH,-1);
    aryubbGL[4]=pow(MMW,-1);
    aryubbGL[5]=Tsil::I2(MMH,MMH,MMH,mu2);
    aryubbGL[6]=Tsil::I2(MMH,MMt,MMt,mu2);
    aryubbGL[7]=pow(MMt,-1);
    aryubbGL[8]=Tsil::I2(0,MMH,MMt,mu2);
    aryubbGL[9]=Tsil::I2(0,0,MMH,mu2);
    aryubbGL[10]=Tsil::I2(0,0,MMt,mu2);
    aryubbGL[11]=Tsil::B(MMH,MMH,MMH,mu2);
    aryubbGL[12]=Tsil::B(MMH,MMt,MMt,mu2);
    aryubbGL[13]=Tsil::A(MMt,mu2);
    aryubbGL[14]=Tsil::B(MMt,MMt,MMH,mu2);
    aryubbGL[15]=Tsil::A(MMH,mu2);
    aryubbGL[16]=std::real(Tsil::B(0,0,MMt,mu2));
    aryubbGL[17]=Tsil::Aeps(MMH,mu2);
    aryubbGL[18]=Tsil::Aeps(MMt,mu2);
    aryubbGL[19]=std::real(Tsil::B(0,0,MMH,mu2));
   aryubbGL[20]=3*aryubbGL[14];
   aryubbGL[21]=pow(Pi,2);
   aryubbGL[22]= - 149./8. - 11*aryubbGL[16];
   aryubbGL[22]= - 35./16.*aryubbGL[21] + 1./8.*aryubbGL[22] - 
   aryubbGL[20];
   aryubbGL[23]=MMt*aryubbGL[3];
   aryubbGL[24]= - 3 - 1./2.*aryubbGL[21];
   aryubbGL[24]=aryubbGL[24]*aryubbGL[23];
   aryubbGL[25]=11*aryubbGL[18];
   aryubbGL[26]= - 5./4.*aryubbGL[13] - aryubbGL[25];
   aryubbGL[26]= - 11./4.*aryubbGL[8] + 3*aryubbGL[26];
   aryubbGL[26]=aryubbGL[3]*aryubbGL[26];
   aryubbGL[27]=aryubbGL[15]*aryubbGL[3];
   aryubbGL[28]=pow(aryubbGL[27],2);
   aryubbGL[22]=29./8.*aryubbGL[24] - 63./8.*aryubbGL[28] - 11./4.*
   aryubbGL[12] + 1./2.*aryubbGL[22] + aryubbGL[26];
   aryubbGL[22]=MMt*aryubbGL[22];
   aryubbGL[24]= - 3*aryubbGL[13] + 11./8.*MMH;
   aryubbGL[24]=aryubbGL[12]*aryubbGL[24];
   aryubbGL[23]=13*aryubbGL[23];
   aryubbGL[26]= - 31./8. + aryubbGL[23];
   aryubbGL[26]=aryubbGL[10]*aryubbGL[26];
   aryubbGL[24]=aryubbGL[24] + aryubbGL[26];
   aryubbGL[20]=9./2.*aryubbGL[21] + 197./8. + aryubbGL[20];
   aryubbGL[20]=MMH*aryubbGL[20];
   aryubbGL[26]=85./4. - 3*aryubbGL[16];
   aryubbGL[26]=aryubbGL[13]*aryubbGL[26];
   aryubbGL[20]=aryubbGL[26] + aryubbGL[20];
   aryubbGL[26]=aryubbGL[3]*aryubbGL[13];
   aryubbGL[26]=45./4.*aryubbGL[27] - 7./4. - 11*aryubbGL[26];
   aryubbGL[26]=aryubbGL[15]*aryubbGL[26];
   aryubbGL[27]=pow(aryubbGL[13],2);
   aryubbGL[28]=aryubbGL[3]*aryubbGL[27];
   aryubbGL[20]=aryubbGL[22] - 17./8.*aryubbGL[6] + 1./4.*aryubbGL[26]
    + 49./16.*aryubbGL[8] - 51./4.*aryubbGL[28] + aryubbGL[25] + 1./8.*
   aryubbGL[20] + 1./2.*aryubbGL[24];
   aryubbGL[20]=MMt*aryubbGL[20];
   aryubbGL[22]= - 515./4. + 9*aryubbGL[11];
   aryubbGL[24]=aryubbGL[7]*aryubbGL[18];
   aryubbGL[22]=3./4.*aryubbGL[19] + 1./4.*aryubbGL[22] + aryubbGL[24];
   aryubbGL[24]= - 3./2.*aryubbGL[8] + aryubbGL[9];
   aryubbGL[24]=aryubbGL[7]*aryubbGL[24];
   aryubbGL[21]=243./4.*S2 + 1./2.*aryubbGL[22] - 1./3.*aryubbGL[21] + 
   aryubbGL[24];
   aryubbGL[21]=MMH*aryubbGL[21];
   aryubbGL[22]=aryubbGL[7]*aryubbGL[27];
   aryubbGL[22]=aryubbGL[22] + aryubbGL[18];
   aryubbGL[22]=aryubbGL[8] - 17*aryubbGL[9] - 15*aryubbGL[5] - 13./2.*
   aryubbGL[13] + 7*aryubbGL[22];
   aryubbGL[24]=aryubbGL[12]*aryubbGL[13];
   aryubbGL[25]=MMH*aryubbGL[7];
   aryubbGL[25]=5./2. + aryubbGL[25];
   aryubbGL[25]=aryubbGL[6]*aryubbGL[25];
   aryubbGL[21]=1./2.*aryubbGL[25] + 3./2.*aryubbGL[24] + 1./4.*
   aryubbGL[22] + aryubbGL[21];
   aryubbGL[21]=MMH*aryubbGL[21];
   aryubbGL[22]= - aryubbGL[7]*aryubbGL[13];
   aryubbGL[22]=61./2. + aryubbGL[22];
   aryubbGL[22]=MMH*aryubbGL[22];
   aryubbGL[22]= - 27./2.*aryubbGL[15] + 21./2.*aryubbGL[13] + 
   aryubbGL[22];
   aryubbGL[22]=aryubbGL[15]*aryubbGL[22];
   aryubbGL[21]=1./2.*aryubbGL[22] + 123./8.*aryubbGL[27] + 
   aryubbGL[21];
   aryubbGL[22]=75./16. - aryubbGL[23];
   aryubbGL[22]=MMt*aryubbGL[22];
   aryubbGL[22]=7./2.*MMH + aryubbGL[22];
   aryubbGL[22]=aryubbGL[17]*aryubbGL[22];
   aryubbGL[20]=aryubbGL[22] + 1./4.*aryubbGL[21] + aryubbGL[20];

      yubbGLret = 1./4.*aryubbGL[20]*pow(aryubbGL[4],2)*pow(
      aryubbGL[2],4)*aryubbGL[1];
      return yubbGLret.real();
}
