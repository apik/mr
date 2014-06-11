#include <bb.hpp>
std::complex<long double> bb::mgl20(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mbbgl[26];

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
    mbbgl[14]=log(MMH);
    mbbgl[15]=log(MMt);
   mbbgl[16]=pow(Pi,2);
   mbbgl[17]=mbbgl[15]*mbbgl[9];
   mbbgl[18]= - 51*mbbgl[17] + 667./4. + 573*mbbgl[9];
   mbbgl[19]= - 11./2. + 15*mbbgl[9];
   mbbgl[19]=mbbgl[8]*mbbgl[19];
   mbbgl[20]=27*mbbgl[6] - mbbgl[7];
   mbbgl[20]=mbbgl[2]*mbbgl[20];
   mbbgl[21]=3*mbbgl[9];
   mbbgl[22]=mbbgl[21] + 1./2.;
   mbbgl[23]=mbbgl[13]*mbbgl[22];
   mbbgl[24]=mbbgl[9] - 1;
   mbbgl[25]=mbbgl[12]*mbbgl[24];
   mbbgl[18]=125./96.*mbbgl[16] + 9*mbbgl[25] + 1./4.*mbbgl[23] + 
   mbbgl[20] + 1./16.*mbbgl[18] + mbbgl[19];
   mbbgl[19]=mbbgl[24]*mbbgl[10];
   mbbgl[20]=19*mbbgl[17] + 47 - 165./4.*mbbgl[9];
   mbbgl[23]=1 - 2*mbbgl[9];
   mbbgl[23]=mbbgl[8]*mbbgl[23];
   mbbgl[24]=mbbgl[9] - 1./2.;
   mbbgl[25]=mbbgl[13]*mbbgl[24];
   mbbgl[20]=9./2.*mbbgl[19] + 9./32.*mbbgl[16] - 3./2.*mbbgl[25] + 1./
   4.*mbbgl[20] + 3*mbbgl[23];
   mbbgl[20]=mbbgl[2]*mbbgl[20];
   mbbgl[19]= - 2*mbbgl[19] - 1./2.*mbbgl[17] + mbbgl[24];
   mbbgl[19]=MMt*mbbgl[19]*pow(mbbgl[2],2);
   mbbgl[19]=9*mbbgl[19] + mbbgl[20];
   mbbgl[19]=MMt*mbbgl[19];
   mbbgl[18]=1./8.*mbbgl[18] + mbbgl[19];
   mbbgl[18]=MMt*mbbgl[18];
   mbbgl[19]=mbbgl[8]*mbbgl[22];
   mbbgl[19]= - 41./96.*mbbgl[16] - 1./8.*mbbgl[19] - 11./32.*mbbgl[17]
    - 47./4. + mbbgl[21];
   mbbgl[19]=MMH*mbbgl[19];
   mbbgl[20]= - 7*mbbgl[6] - 1./2.*mbbgl[7];
   mbbgl[19]=5./8.*mbbgl[20] + mbbgl[19];
   mbbgl[18]=1./4.*mbbgl[19] + mbbgl[18];
   mbbgl[18]=MMt*mbbgl[18];
   mbbgl[19]=1./2.*mbbgl[4];
   mbbgl[20]=mbbgl[6] - 3*mbbgl[7];
   mbbgl[19]=mbbgl[20]*mbbgl[19];
   mbbgl[20]=MMH*mbbgl[4];
   mbbgl[21]= - 7 - 1./2.*mbbgl[16];
   mbbgl[21]=mbbgl[21]*mbbgl[20];
   mbbgl[16]=1./2.*mbbgl[21] + mbbgl[19] - 7./48.*mbbgl[16] - 1./4.*
   mbbgl[17] - 17 + mbbgl[9];
   mbbgl[16]=MMH*mbbgl[16];
   mbbgl[17]= - 15*mbbgl[5] - 13*mbbgl[6] + 19*mbbgl[7];
   mbbgl[16]=1./4.*mbbgl[17] + mbbgl[16];
   mbbgl[17]=1./16.*MMH;
   mbbgl[16]=mbbgl[16]*mbbgl[17];
   mbbgl[19]=3 - mbbgl[14];
   mbbgl[19]=mbbgl[19]*mbbgl[20];
   mbbgl[20]=35 - 17*mbbgl[14];
   mbbgl[20]=9./2.*mbbgl[12] + 5./4.*mbbgl[20] - mbbgl[15];
   mbbgl[19]=1./2.*mbbgl[20] + mbbgl[19];
   mbbgl[19]=mbbgl[19]*pow(MMH,2);
   mbbgl[20]=mbbgl[10] - mbbgl[15];
   mbbgl[21]=271 - 75*mbbgl[14];
   mbbgl[20]=1./8.*mbbgl[21] + 9*mbbgl[20];
   mbbgl[17]=mbbgl[20]*mbbgl[17];
   mbbgl[20]= - 25 + 13*mbbgl[14];
   mbbgl[20]= - 9./4.*mbbgl[10] + 1./8.*mbbgl[20] - 2*mbbgl[15];
   mbbgl[20]=MMt*mbbgl[20];
   mbbgl[17]=mbbgl[17] + mbbgl[20];
   mbbgl[17]=MMt*mbbgl[17];
   mbbgl[17]=1./16.*mbbgl[19] + mbbgl[17];
   mbbgl[17]=mbbgl[11]*mbbgl[17];
   mbbgl[16]=mbbgl[17] + mbbgl[16] + mbbgl[18];

      return mbbgl[16]*pow(mbbgl[3],2)*pow(mbbgl[1],4);
}
