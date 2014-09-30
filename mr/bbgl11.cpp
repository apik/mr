#include <bb.hpp>
std::complex<long double> bb::mgl11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mbbgl[12];

    mbbgl[1]=pow(SW,-1);
    mbbgl[2]=pow(MMH,-1);
    mbbgl[3]=pow(MMW,-1);
    mbbgl[4]=log(pow(mu2,-1)*MMb);
    mbbgl[5]=log(pow(mu2,-1)*MMt);
    mbbgl[6]=log(pow(mu2,-1)*MMH);
   mbbgl[7]= - mbbgl[5]*mbbgl[2];
   mbbgl[8]=8*mbbgl[2] + 3*mbbgl[7];
   mbbgl[8]=mbbgl[5]*mbbgl[8];
   mbbgl[8]= - mbbgl[2] + mbbgl[8];
   mbbgl[7]=mbbgl[2] + mbbgl[7];
   mbbgl[7]=mbbgl[4]*mbbgl[7];
   mbbgl[7]=2*mbbgl[8] + 3*mbbgl[7];
   mbbgl[7]=MMt*mbbgl[7];
   mbbgl[8]=pow(Pi,2);
   mbbgl[9]= - 13 + 2*mbbgl[8];
   mbbgl[10]= - 2 + 3./2.*mbbgl[5];
   mbbgl[10]=mbbgl[5]*mbbgl[10];
   mbbgl[11]= - 5./2. + 3*mbbgl[5];
   mbbgl[11]=mbbgl[4]*mbbgl[11];
   mbbgl[7]=4*mbbgl[7] + 1./2.*mbbgl[11] + 1./3.*mbbgl[9] + mbbgl[10];
   mbbgl[7]=MMt*mbbgl[7];
   mbbgl[9]= - 1 + mbbgl[6];
   mbbgl[9]=MMH*mbbgl[9];
   mbbgl[10]= - mbbgl[4]*MMb;
   mbbgl[9]=5*mbbgl[10] + 3./2.*mbbgl[9] + 10*MMb;
   mbbgl[9]=mbbgl[4]*mbbgl[9];
   mbbgl[10]=1 - mbbgl[6];
   mbbgl[10]=MMH*mbbgl[10];
   mbbgl[8]= - 3 - 1./6.*mbbgl[8];
   mbbgl[8]=MMb*mbbgl[8];
   mbbgl[7]=mbbgl[7] + mbbgl[9] + 2*mbbgl[10] + 5./2.*mbbgl[8];

      return mbbgl[7]*mbbgl[3]*pow(mbbgl[1],2);
}
