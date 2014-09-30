#include <tt.hpp>
std::complex<long double> tt::drgl20(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mdrgl[25];

    mdrgl[1]=pow(SW,-1);
    mdrgl[2]=pow(MMH,-1);
    mdrgl[3]=pow(MMW,-1);
    mdrgl[4]=Tsil::I2(MMH,MMt,MMt,mu2);
    mdrgl[5]=Tsil::I2(0,MMH,MMt,mu2);
    mdrgl[6]=Tsil::B(MMH,MMH,MMH,mu2);
    mdrgl[7]=Tsil::A(MMH,mu2);
    mdrgl[8]=Tsil::A(MMt,mu2);
    mdrgl[9]=Tsil::B(MMH,MMt,MMt,mu2);
    mdrgl[10]=Tsil::B(MMt,MMt,MMH,mu2);
    mdrgl[11]=std::real(Tsil::B(0,0,MMH,mu2));
    mdrgl[12]=pow(MMt,-1);
    mdrgl[13]=std::real(Tsil::B(0,0,MMt,mu2));
   mdrgl[14]=pow(mdrgl[7],2);
   mdrgl[15]=5*mdrgl[7];
   mdrgl[16]=mdrgl[15] + 11*mdrgl[8];
   mdrgl[16]=mdrgl[8]*mdrgl[16];
   mdrgl[16]= - 15./4.*mdrgl[14] + mdrgl[16];
   mdrgl[16]=mdrgl[2]*mdrgl[16];
   mdrgl[17]=mdrgl[13]*mdrgl[8];
   mdrgl[16]=mdrgl[17] + mdrgl[16];
   mdrgl[18]=mdrgl[10]*mdrgl[7];
   mdrgl[19]=mdrgl[18] - mdrgl[4];
   mdrgl[20]=3./4. + mdrgl[11];
   mdrgl[20]=mdrgl[8]*mdrgl[20];
   mdrgl[19]=3./2.*mdrgl[20] - 9./8.*mdrgl[5] + mdrgl[7] + 3./4.*
   mdrgl[19];
   mdrgl[20]=pow(Pi,2);
   mdrgl[21]=5./16.*mdrgl[10] - 1 - 3./32.*mdrgl[20];
   mdrgl[21]=MMH*mdrgl[21];
   mdrgl[22]=9./4.*mdrgl[6];
   mdrgl[23]=mdrgl[8]*mdrgl[22];
   mdrgl[24]=mdrgl[8] - 1./8.*MMH;
   mdrgl[24]=mdrgl[9]*mdrgl[24];
   mdrgl[16]=3./2.*mdrgl[24] + mdrgl[23] + 1./2.*mdrgl[19] + mdrgl[21]
    + 1./8.*mdrgl[16];
   mdrgl[19]=3./2.*mdrgl[5] + 9*mdrgl[4] - mdrgl[15];
   mdrgl[18]=1./2.*mdrgl[19] - 3*mdrgl[18];
   mdrgl[19]= - 19./8. + 3*mdrgl[10];
   mdrgl[19]=mdrgl[8]*mdrgl[19];
   mdrgl[21]=mdrgl[2]*mdrgl[14];
   mdrgl[17]=21./16.*mdrgl[21] - mdrgl[17] + 1./2.*mdrgl[18] + 
   mdrgl[19];
   mdrgl[18]=231 + 29./2.*mdrgl[20];
   mdrgl[18]=1./8.*mdrgl[18] - 3*mdrgl[13];
   mdrgl[19]=mdrgl[2]*mdrgl[8];
   mdrgl[21]=mdrgl[10]*mdrgl[19];
   mdrgl[18]= - 6*mdrgl[9] + 1./2.*mdrgl[18] - 36*mdrgl[21];
   mdrgl[18]=MMt*mdrgl[18];
   mdrgl[17]=mdrgl[18] + 3*mdrgl[17];
   mdrgl[17]=mdrgl[2]*mdrgl[17];
   mdrgl[18]=3*mdrgl[9];
   mdrgl[19]=5./4. - 4*mdrgl[19];
   mdrgl[19]=mdrgl[19]*mdrgl[18];
   mdrgl[21]=177 + 35*mdrgl[20];
   mdrgl[21]=9./4.*mdrgl[13] + 1./16.*mdrgl[21] - 15*mdrgl[10];
   mdrgl[17]=mdrgl[19] + 1./4.*mdrgl[21] + mdrgl[17];
   mdrgl[17]=MMt*mdrgl[17];
   mdrgl[16]=3*mdrgl[16] + mdrgl[17];
   mdrgl[16]=MMt*mdrgl[16];
   mdrgl[17]=mdrgl[5] - mdrgl[4];
   mdrgl[19]=3./2.*mdrgl[11] - 11./2.;
   mdrgl[19]=mdrgl[7]*mdrgl[19];
   mdrgl[21]=mdrgl[12]*pow(mdrgl[8],2);
   mdrgl[17]= - 3./2.*mdrgl[21] + mdrgl[19] + 3*mdrgl[17];
   mdrgl[19]=5./2.*MMH + 3*mdrgl[7];
   mdrgl[19]=mdrgl[19]*mdrgl[22];
   mdrgl[20]=15./8.*mdrgl[11] + 131./16. + 1./3.*mdrgl[20];
   mdrgl[20]=MMH*mdrgl[20];
   mdrgl[18]= - mdrgl[8]*mdrgl[18];
   mdrgl[17]=mdrgl[18] + mdrgl[19] + 3./2.*mdrgl[17] + mdrgl[20];
   mdrgl[17]=MMH*mdrgl[17];
   mdrgl[15]= - mdrgl[15] + mdrgl[8];
   mdrgl[15]=mdrgl[8]*mdrgl[15];
   mdrgl[14]=9./2.*mdrgl[14] + mdrgl[15];
   mdrgl[15]=S2*pow(MMH,2);
   mdrgl[14]= - 243./4.*mdrgl[15] + 3./2.*mdrgl[14] + mdrgl[17];
   mdrgl[14]=1./8.*mdrgl[14] + mdrgl[16];

      return mdrgl[14]*pow(mdrgl[3],2)*pow(mdrgl[1],4);
}
