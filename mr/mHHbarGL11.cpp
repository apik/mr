#include <HH.hpp>
std::complex<long double>
HH<OS>::mgl11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armHHbarGL[16], mHHbarGLret;

    armHHbarGL[1]=double(boson);
    armHHbarGL[2]=pow(SW,-1);
    armHHbarGL[3]=pow(MMH,-1);
    armHHbarGL[4]=pow(MMW,-1);
    armHHbarGL[5]=Tsil::B(MMt,MMt,MMH,mu2);
    armHHbarGL[6]=Tsil::A(MMt,mu2);
    armHHbarGL[7]=Tsil::Beps(MMt,MMt,MMH,mu2);
    armHHbarGL[8]=Tsil::Aeps(MMt,mu2);
    armHHbarGL[9]=prottttt0->M(0);
    armHHbarGL[10]=prottttt0->Vxzuv(0);
    armHHbarGL[11]=prottttt0->Suxv(0);
   armHHbarGL[12]=4 + 3*armHHbarGL[5];
   armHHbarGL[12]=armHHbarGL[5]*armHHbarGL[12];
   armHHbarGL[13]=armHHbarGL[9]*MMH;
   armHHbarGL[14]=2*armHHbarGL[10];
   armHHbarGL[15]= - armHHbarGL[14] - 3*armHHbarGL[9];
   armHHbarGL[15]=MMt*armHHbarGL[15];
   armHHbarGL[12]=4*armHHbarGL[15] + 2*armHHbarGL[13] + 12 + 
   armHHbarGL[12];
   armHHbarGL[12]=MMt*armHHbarGL[12];
   armHHbarGL[13]= - 10 - 7*armHHbarGL[5];
   armHHbarGL[13]=armHHbarGL[5]*armHHbarGL[13];
   armHHbarGL[14]=armHHbarGL[14] + armHHbarGL[9];
   armHHbarGL[14]=MMt*armHHbarGL[14];
   armHHbarGL[13]=8*armHHbarGL[14] + 4*armHHbarGL[7] - 9 + 
   armHHbarGL[13];
   armHHbarGL[13]=MMt*armHHbarGL[13];
   armHHbarGL[14]=armHHbarGL[6]*armHHbarGL[5];
   armHHbarGL[13]=armHHbarGL[13] - armHHbarGL[11] - 18*armHHbarGL[14];
   armHHbarGL[13]=MMt*armHHbarGL[13];
   armHHbarGL[15]=pow(armHHbarGL[6],2);
   armHHbarGL[13]= - 12*armHHbarGL[15] + armHHbarGL[13];
   armHHbarGL[15]=armHHbarGL[8]*MMt;
   armHHbarGL[13]=12*armHHbarGL[15] + 2*armHHbarGL[13];
   armHHbarGL[13]=armHHbarGL[3]*armHHbarGL[13];
   armHHbarGL[12]=6*armHHbarGL[14] + armHHbarGL[12] + armHHbarGL[13];

      mHHbarGLret = 2*armHHbarGL[12]*armHHbarGL[4]*pow(armHHbarGL[2],2)
      *armHHbarGL[1];
      return mHHbarGLret;
}
