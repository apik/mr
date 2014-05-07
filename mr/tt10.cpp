#include <tt.hpp>
std::complex<long double> tt::m10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mtt[21];

    mtt[1]=pow(CW,-1);
    mtt[2]=pow(MMH,-1);
    mtt[3]=pow(SW,-1);
    mtt[4]=Tsil::B(MMH,MMt,MMt,mu2);
    mtt[5]=pow(MMZ,-1);
    mtt[6]=Tsil::B(MMZ,MMt,MMt,mu2);
    mtt[7]=pow(MMt,-1);
    mtt[8]=Tsil::B(0,MMW,MMt,mu2);
    mtt[9]=Tsil::A(MMH,mu2);
    mtt[10]=Tsil::A(MMZ,mu2);
    mtt[11]=Tsil::A(MMW,mu2);
    mtt[12]=Tsil::A(MMt,mu2);
   mtt[13]=pow(mtt[3],2);
   mtt[14]=1./8.*mtt[13];
   mtt[15]=pow(mtt[1],2);
   mtt[16]=mtt[14] + 17./72.*mtt[15];
   mtt[17]=mtt[16] - 4./9.;
   mtt[18]= - mtt[6]*mtt[17];
   mtt[19]=1./4.*mtt[8];
   mtt[20]=1 - mtt[13];
   mtt[20]=mtt[20]*mtt[19];
   mtt[18]=mtt[20] + mtt[18];
   mtt[18]=MMZ*mtt[18];
   mtt[16]=8./9. + mtt[16];
   mtt[16]=mtt[12]*mtt[16];
   mtt[17]= - mtt[10]*mtt[17];
   mtt[20]=1./4.*mtt[11];
   mtt[21]= - mtt[13]*mtt[20];
   mtt[16]=mtt[21] + mtt[17] + mtt[16] + mtt[18];
   mtt[16]=mtt[7]*mtt[16];
   mtt[17]= - mtt[4]*MMH;
   mtt[17]=mtt[12] + 1./2.*mtt[17];
   mtt[18]=mtt[4] + mtt[19];
   mtt[18]=MMt*mtt[18];
   mtt[17]=1./2.*mtt[17] + mtt[18] + 1./2.*mtt[9] + mtt[20];
   mtt[18]=mtt[15] + mtt[13];
   mtt[17]=mtt[5]*mtt[18]*mtt[17];
   mtt[19]=mtt[5]*MMt*mtt[12];
   mtt[19]= - 3*mtt[19] + 3./4.*mtt[10];
   mtt[18]=mtt[18]*mtt[19];
   mtt[19]=3./2.*mtt[13];
   mtt[20]=1./2.*mtt[15] - 1 + mtt[19];
   mtt[20]=MMZ*mtt[20];
   mtt[19]=mtt[11]*mtt[19];
   mtt[18]=mtt[19] + mtt[20] + mtt[18];
   mtt[18]=mtt[2]*mtt[18];
   mtt[19]=1./8.*mtt[8] - 3./8.;
   mtt[13]=mtt[13]*mtt[19];
   mtt[14]= - 7./72.*mtt[15] + mtt[14] + 8./9.;
   mtt[14]=mtt[6]*mtt[14];

      return  - 8./9. + mtt[13] + mtt[14] - 1./72.*mtt[15] + mtt[16] + 
      1./2.*mtt[17] + mtt[18];
}
