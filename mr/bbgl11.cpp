#include <bb.hpp>
std::complex<long double> bb::mgl11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mbbgl[13];

    mbbgl[1]=pow(SW,-1);
    mbbgl[2]=pow(MMH,-1);
    mbbgl[3]=pow(MMW,-1);
    mbbgl[4]=log(pow(mu2,-1)*MMb);
    mbbgl[5]=log(pow(mu2,-1)*MMt);
    mbbgl[6]=log(pow(mu2,-1)*MMH);
    mbbgl[7]=log(MMt);
    mbbgl[8]=log(MMb);
   mbbgl[9]=pow(Pi,2);
   mbbgl[10]= - 2 + 3./2.*mbbgl[8];
   mbbgl[11]=3./2.*mbbgl[7] + mbbgl[10];
   mbbgl[11]=mbbgl[5]*mbbgl[11];
   mbbgl[12]= - 6*mbbgl[7] + 16 - 3*mbbgl[8];
   mbbgl[12]=mbbgl[5]*mbbgl[12];
   mbbgl[12]=mbbgl[12] - 2 + 3*mbbgl[4];
   mbbgl[12]=MMt*mbbgl[2]*mbbgl[12];
   mbbgl[11]=4*mbbgl[12] + mbbgl[11] + 2./3.*mbbgl[9] - 13./3. - 5./4.*
   mbbgl[4];
   mbbgl[11]=MMt*mbbgl[11];
   mbbgl[12]=2 - mbbgl[8];
   mbbgl[12]=mbbgl[4]*mbbgl[12];
   mbbgl[9]= - 1./12.*mbbgl[9] - 3./2. + mbbgl[12];
   mbbgl[9]=MMb*mbbgl[9];
   mbbgl[10]=mbbgl[6]*mbbgl[10];
   mbbgl[10]=mbbgl[10] + 2 - 3./2.*mbbgl[4];
   mbbgl[10]=MMH*mbbgl[10];
   mbbgl[9]=mbbgl[11] + 5*mbbgl[9] + mbbgl[10];

      return mbbgl[9]*mbbgl[3]*pow(mbbgl[1],2);
}
