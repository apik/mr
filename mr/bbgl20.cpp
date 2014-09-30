#include <bb.hpp>
std::complex<long double> bb::mgl20(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mbbgl[21];

    mbbgl[1]=pow(SW,-1);
    mbbgl[2]=pow(MMH,-1);
    mbbgl[3]=pow(MMW,-1);
    mbbgl[4]=pow(MMt,-1);
    mbbgl[5]=Tsil::I2(MMH,MMH,MMH,mu2);
    mbbgl[6]=Tsil::I2(MMH,MMt,MMt,mu2);
    mbbgl[7]=Tsil::I2(0,MMH,MMt,mu2);
    mbbgl[8]=Tsil::B(MMH,MMt,MMt,mu2);
    mbbgl[9]=log(pow(mu2,-1)*MMt);
    mbbgl[10]=Tsil::B(MMt,MMt,MMH,mu2);
    mbbgl[11]=log(pow(mu2,-1)*MMH);
    mbbgl[12]=Tsil::B(0,0,MMH,mu2);
    mbbgl[13]=Tsil::B(0,0,MMt,mu2);
   mbbgl[14]=pow(Pi,2);
   mbbgl[15]= - 7 - 1./2.*mbbgl[14];
   mbbgl[15]=mbbgl[4]*mbbgl[15];
   mbbgl[16]= - mbbgl[11]*mbbgl[4];
   mbbgl[16]=3*mbbgl[4] + mbbgl[16];
   mbbgl[16]=mbbgl[11]*mbbgl[16];
   mbbgl[15]=1./2.*mbbgl[15] + mbbgl[16];
   mbbgl[15]=MMH*mbbgl[15];
   mbbgl[16]= - 3*mbbgl[7] + mbbgl[6];
   mbbgl[16]=mbbgl[4]*mbbgl[16];
   mbbgl[17]= - 85./2.*mbbgl[11] + 175./2. + 9*mbbgl[12];
   mbbgl[17]=mbbgl[11]*mbbgl[17];
   mbbgl[18]= - 1./4.*mbbgl[9] + 1 - 1./2.*mbbgl[11];
   mbbgl[18]=mbbgl[9]*mbbgl[18];
   mbbgl[15]=mbbgl[15] + mbbgl[18] + 1./4.*mbbgl[17] + 1./2.*mbbgl[16]
    - 17 - 7./48.*mbbgl[14];
   mbbgl[15]=MMH*mbbgl[15];
   mbbgl[16]= - 13*mbbgl[6] - 15*mbbgl[5] + 19*mbbgl[7];
   mbbgl[15]=1./4.*mbbgl[16] + mbbgl[15];
   mbbgl[15]=MMH*mbbgl[15];
   mbbgl[16]= - 1./2. + 2*mbbgl[10];
   mbbgl[17]=pow(mbbgl[2],2);
   mbbgl[16]=mbbgl[17]*mbbgl[16];
   mbbgl[18]=1 - 2*mbbgl[10];
   mbbgl[18]=mbbgl[17]*mbbgl[18];
   mbbgl[17]= - mbbgl[9]*mbbgl[17];
   mbbgl[17]=mbbgl[18] + 1./2.*mbbgl[17];
   mbbgl[17]=mbbgl[9]*mbbgl[17];
   mbbgl[16]=mbbgl[16] + mbbgl[17];
   mbbgl[16]=MMt*mbbgl[16];
   mbbgl[17]=9./8.*mbbgl[14] + 47 + 3*mbbgl[13];
   mbbgl[18]= - 9*mbbgl[10];
   mbbgl[17]=1./2.*mbbgl[17] + mbbgl[18];
   mbbgl[17]=1./2.*mbbgl[17] + 3*mbbgl[8];
   mbbgl[17]=mbbgl[2]*mbbgl[17];
   mbbgl[19]=3*mbbgl[10] - 55./8. - mbbgl[13];
   mbbgl[19]=1./2.*mbbgl[19] - 2*mbbgl[8];
   mbbgl[19]=mbbgl[2]*mbbgl[19];
   mbbgl[20]=mbbgl[9]*mbbgl[2];
   mbbgl[19]=3*mbbgl[19] + 19./4.*mbbgl[20];
   mbbgl[19]=mbbgl[9]*mbbgl[19];
   mbbgl[16]=9*mbbgl[16] + mbbgl[17] + mbbgl[19];
   mbbgl[16]=MMt*mbbgl[16];
   mbbgl[17]= - 11./2.*mbbgl[8] + 125./96.*mbbgl[14] + 1./8.*mbbgl[13]
    + 667./64. - 9*mbbgl[12];
   mbbgl[18]=13./2.*mbbgl[11] - 25./2. + mbbgl[18];
   mbbgl[18]=mbbgl[11]*mbbgl[18];
   mbbgl[19]= - mbbgl[7] + 27*mbbgl[6];
   mbbgl[19]=mbbgl[2]*mbbgl[19];
   mbbgl[17]=1./2.*mbbgl[19] + 1./2.*mbbgl[17] + mbbgl[18];
   mbbgl[18]=5*mbbgl[8] + 1./4.*mbbgl[13] + 191./16. + 3*mbbgl[12];
   mbbgl[18]= - 51./128.*mbbgl[9] + 3./8.*mbbgl[18] - 2*mbbgl[11];
   mbbgl[18]=mbbgl[9]*mbbgl[18];
   mbbgl[16]=mbbgl[16] + 1./4.*mbbgl[17] + mbbgl[18];
   mbbgl[16]=MMt*mbbgl[16];
   mbbgl[17]= - 75./8.*mbbgl[11] + 271./8. + 9*mbbgl[10];
   mbbgl[17]=mbbgl[11]*mbbgl[17];
   mbbgl[14]=mbbgl[17] - 1./4.*mbbgl[8] - 47 - 41./24.*mbbgl[14];
   mbbgl[17]= - 3./4.*mbbgl[11] + 1 - 1./8.*mbbgl[8];
   mbbgl[17]=3*mbbgl[17] - 11./32.*mbbgl[9];
   mbbgl[17]=mbbgl[9]*mbbgl[17];
   mbbgl[14]=1./4.*mbbgl[14] + mbbgl[17];
   mbbgl[14]=MMH*mbbgl[14];
   mbbgl[17]= - 1./2.*mbbgl[7] - 7*mbbgl[6];
   mbbgl[14]=5./8.*mbbgl[17] + mbbgl[14];
   mbbgl[14]=1./4.*mbbgl[14] + mbbgl[16];
   mbbgl[14]=MMt*mbbgl[14];
   mbbgl[14]=1./16.*mbbgl[15] + mbbgl[14];

      return mbbgl[14]*pow(mbbgl[3],2)*pow(mbbgl[1],4);
}
