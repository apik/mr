#include <tt.hpp>
std::complex<long double> tt::m10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mtt[24], mttret;

    mtt[1]=pow(CW,-1);
    mtt[2]=pow(MMH,-1);
    mtt[3]=pow(SW,-1);
    mtt[4]=double(nH);
    mtt[5]=Tsil::A(MMt,mu2);
    mtt[6]=pow(MMZ,-1);
    mtt[7]=Tsil::B(MMH,MMt,MMt,mu2);
    mtt[8]=Tsil::B(MMZ,MMt,MMt,mu2);
    mtt[9]=pow(MMt,-1);
    mtt[10]=Tsil::A(MMH,mu2);
    mtt[11]=Tsil::A(MMZ,mu2);
    mtt[12]=Tsil::A(MMW,mu2);
    mtt[13]=std::real(Tsil::B(0,MMW,MMt,mu2));
   mtt[14]=pow(mtt[1],2);
   mtt[15]=pow(mtt[3],2);
   mtt[16]=1./2.*mtt[14] - 1 + 3./2.*mtt[15];
   mtt[16]=mtt[2]*mtt[16];
   mtt[17]=1./4.*mtt[9];
   mtt[18]=1 - mtt[15];
   mtt[18]=mtt[13]*mtt[18]*mtt[17];
   mtt[16]=mtt[16] + mtt[18];
   mtt[16]=MMZ*mtt[16];
   mtt[18]=1./8.*mtt[15];
   mtt[19]=mtt[18] + 17./72.*mtt[14];
   mtt[20]=mtt[19] - 4./9.;
   mtt[20]=mtt[20]*mtt[9];
   mtt[21]= - MMZ*mtt[20];
   mtt[18]=mtt[21] - 7./72.*mtt[14] + mtt[18] + 8./9.;
   mtt[18]=mtt[8]*mtt[18];
   mtt[21]=1./8.*mtt[13];
   mtt[22]=1./2.*mtt[7];
   mtt[23]=mtt[5]*mtt[2]*mtt[4];
   mtt[23]=mtt[22] + mtt[21] - 3*mtt[23];
   mtt[23]=MMt*mtt[23];
   mtt[22]= - MMH*mtt[22];
   mtt[22]=mtt[22] + mtt[5] + 1./2.*mtt[12] + mtt[10];
   mtt[22]=1./4.*mtt[22] + mtt[23];
   mtt[23]=mtt[14] + mtt[15];
   mtt[22]=mtt[6]*mtt[23]*mtt[22];
   mtt[23]=mtt[2]*mtt[23];
   mtt[20]= - mtt[20] + 3./4.*mtt[23];
   mtt[20]=mtt[11]*mtt[20];
   mtt[17]=3./2.*mtt[2] - mtt[17];
   mtt[17]=mtt[12]*mtt[17];
   mtt[17]=mtt[21] - 3./8. + mtt[17];
   mtt[15]=mtt[15]*mtt[17];
   mtt[17]=8./9. + mtt[19];
   mtt[17]=mtt[5]*mtt[9]*mtt[17];

      mttret =  - 8./9. - 1./72.*mtt[14] + mtt[15] + mtt[16] + mtt[17]
       + mtt[18] + mtt[20] + mtt[22];
      return mttret;
}
