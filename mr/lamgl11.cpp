#include <HH.hpp>
std::complex<long double> HH::lamgl11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mlamgl[18];

    mlamgl[1]=pow(SW,-1);
    mlamgl[2]=pow(MMH,-1);
    mlamgl[3]=pow(MMW,-1);
    mlamgl[4]=Tsil::B(MMt,MMt,MMH,mu2);
    mlamgl[5]=log(pow(mu2,-1)*MMt);
    mlamgl[6]=Tsil::Beps(MMt,MMt,MMH,mu2);
    mlamgl[7]=log(MMt);
    mlamgl[8]=prottttt0->M(0);
    mlamgl[9]=prottttt0->Vxzuv(0);
    mlamgl[10]=prottttt0->Suxv(0);
   mlamgl[11]=6*mlamgl[5];
   mlamgl[12]=3*mlamgl[4];
   mlamgl[13]=1 - mlamgl[12];
   mlamgl[13]=2*mlamgl[13] - mlamgl[7];
   mlamgl[13]=mlamgl[13]*mlamgl[11];
   mlamgl[14]=pow(Pi,2);
   mlamgl[15]=8 - 7*mlamgl[4];
   mlamgl[15]=mlamgl[4]*mlamgl[15];
   mlamgl[15]= - 31 + mlamgl[15];
   mlamgl[16]=2*mlamgl[9];
   mlamgl[17]=mlamgl[16] + mlamgl[8];
   mlamgl[17]=MMt*mlamgl[17];
   mlamgl[13]=16*mlamgl[17] + 8*mlamgl[6] + mlamgl[13] + 2*mlamgl[15]
    - mlamgl[14];
   mlamgl[13]=mlamgl[2]*mlamgl[13];
   mlamgl[15]= - mlamgl[16] - 3*mlamgl[8];
   mlamgl[13]=4*mlamgl[15] + mlamgl[13];
   mlamgl[13]=MMt*mlamgl[13];
   mlamgl[15]=mlamgl[8]*MMH;
   mlamgl[16]=mlamgl[2]*mlamgl[10];
   mlamgl[15]=mlamgl[15] - mlamgl[16];
   mlamgl[16]=2*mlamgl[4];
   mlamgl[12]= - 2 + mlamgl[12];
   mlamgl[12]=mlamgl[12]*mlamgl[16];
   mlamgl[16]= - mlamgl[7] + 1 + mlamgl[16];
   mlamgl[11]=mlamgl[16]*mlamgl[11];
   mlamgl[11]=2*mlamgl[13] + mlamgl[11] - 2./3.*mlamgl[14] + 77./2. + 
   mlamgl[12] + 4*mlamgl[15];

      return MMt*mlamgl[11]*mlamgl[3]*pow(mlamgl[1],2);
}
