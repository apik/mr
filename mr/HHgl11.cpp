#include <HH.hpp>
std::complex<long double> HH::mgl11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mHHgl[17];

    mHHgl[1]=pow(SW,-1);
    mHHgl[2]=pow(MMH,-1);
    mHHgl[3]=pow(MMW,-1);
    mHHgl[4]=Tsil::B(MMt,MMt,MMH,mu2);
    mHHgl[5]=log(pow(mu2,-1)*MMt);
    mHHgl[6]=Tsil::Beps(MMt,MMt,MMH,mu2);
    mHHgl[7]=log(MMt);
    mHHgl[8]=prottttt0->M(0);
    mHHgl[9]=prottttt0->Vxzuv(0);
    mHHgl[10]=prottttt0->Suxv(0);
   mHHgl[11]=6*mHHgl[5];
   mHHgl[12]=2 - mHHgl[7];
   mHHgl[12]=5*mHHgl[12] - 6*mHHgl[4];
   mHHgl[12]=mHHgl[12]*mHHgl[11];
   mHHgl[13]=8 - 7*mHHgl[4];
   mHHgl[13]=mHHgl[4]*mHHgl[13];
   mHHgl[13]=4*mHHgl[6] - 27 + mHHgl[13];
   mHHgl[14]=pow(Pi,2);
   mHHgl[15]=2*mHHgl[9];
   mHHgl[16]=mHHgl[8] + mHHgl[15];
   mHHgl[16]=MMt*mHHgl[16];
   mHHgl[12]=2*mHHgl[13] - mHHgl[14] + 16*mHHgl[16] + mHHgl[12];
   mHHgl[12]=mHHgl[2]*mHHgl[12];
   mHHgl[13]= - 3*mHHgl[8] - mHHgl[15];
   mHHgl[12]=4*mHHgl[13] + mHHgl[12];
   mHHgl[12]=MMt*mHHgl[12];
   mHHgl[13]=MMH*mHHgl[8];
   mHHgl[14]=mHHgl[10]*mHHgl[2];
   mHHgl[13]=mHHgl[13] - mHHgl[14];
   mHHgl[11]=mHHgl[11] - 2 + 3*mHHgl[4];
   mHHgl[11]=mHHgl[4]*mHHgl[11];
   mHHgl[11]=mHHgl[12] + 12 + 2*mHHgl[13] + mHHgl[11];

      return 2*MMt*mHHgl[11]*mHHgl[3]*pow(mHHgl[1],2);
}
