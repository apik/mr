#include <HH.hpp>
namespace mr
{
  long double HH<OS>::xgl11(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<long double> armHHbarGL[18], mHHbarGLret;

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
    armHHbarGL[12]=armHHbarGL[3]*armHHbarGL[9];
    armHHbarGL[13]=2*armHHbarGL[3];
    armHHbarGL[14]=armHHbarGL[10]*armHHbarGL[13];
    armHHbarGL[12]=armHHbarGL[12] + armHHbarGL[14];
    armHHbarGL[12]=MMt*armHHbarGL[12];
    armHHbarGL[14]= - 10 - 7*armHHbarGL[5];
    armHHbarGL[14]=armHHbarGL[5]*armHHbarGL[14];
    armHHbarGL[14]= - 9 + armHHbarGL[14];
    armHHbarGL[14]=armHHbarGL[3]*armHHbarGL[14];
    armHHbarGL[15]=4*armHHbarGL[3];
    armHHbarGL[16]=armHHbarGL[7]*armHHbarGL[15];
    armHHbarGL[12]=8*armHHbarGL[12] + armHHbarGL[16] - 4*armHHbarGL[10]
      - 6*armHHbarGL[9] + armHHbarGL[14];
    armHHbarGL[12]=MMt*armHHbarGL[12];
    armHHbarGL[14]=MMH*armHHbarGL[9];
    armHHbarGL[12]=armHHbarGL[14] + armHHbarGL[12];
    armHHbarGL[14]=armHHbarGL[5]*armHHbarGL[6];
    armHHbarGL[16]=12*armHHbarGL[8] - 36*armHHbarGL[14];
    armHHbarGL[16]=armHHbarGL[3]*armHHbarGL[16];
    armHHbarGL[17]=4 + 3*armHHbarGL[5];
    armHHbarGL[17]=armHHbarGL[5]*armHHbarGL[17];
    armHHbarGL[13]= - armHHbarGL[11]*armHHbarGL[13];
    armHHbarGL[12]=armHHbarGL[13] + 12 + armHHbarGL[17] + 2*
      armHHbarGL[12] + armHHbarGL[16];
    armHHbarGL[12]=MMt*armHHbarGL[12];
    armHHbarGL[13]= - pow(armHHbarGL[6],2)*armHHbarGL[15];
    armHHbarGL[13]=armHHbarGL[14] + armHHbarGL[13];
    armHHbarGL[12]=6*armHHbarGL[13] + armHHbarGL[12];

    mHHbarGLret = 2*armHHbarGL[12]*armHHbarGL[4]*pow(armHHbarGL[2],2)
      *armHHbarGL[1];
    return mHHbarGLret.real();
  }
} // namespace mr
