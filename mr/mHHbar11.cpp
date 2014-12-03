#include <HH.hpp>
long double HH<OS>::x11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armHHbar[17], mHHbarret;

    armHHbar[1]=double(nH);
    armHHbar[2]=pow(CW,-1);
    armHHbar[3]=pow(MMH,-1);
    armHHbar[4]=pow(MMZ,-1);
    armHHbar[5]=pow(SW,-1);
    armHHbar[6]=Tsil::B(MMt,MMt,MMH,mu2);
    armHHbar[7]=Tsil::A(MMt,mu2);
    armHHbar[8]=Tsil::Beps(MMt,MMt,MMH,mu2);
    armHHbar[9]=Tsil::Aeps(MMt,mu2);
    armHHbar[10]=prottttt0->M(0);
    armHHbar[11]=prottttt0->Vxzuv(0);
    armHHbar[12]=prottttt0->Suxv(0);
   armHHbar[13]=2*armHHbar[11];
   armHHbar[14]=armHHbar[13] + 3*armHHbar[10];
   armHHbar[15]=4*MMt;
   armHHbar[14]=armHHbar[14]*armHHbar[15];
   armHHbar[15]=4 + 3*armHHbar[6];
   armHHbar[15]=armHHbar[15]*armHHbar[6];
   armHHbar[16]=armHHbar[10]*MMH;
   armHHbar[16]=armHHbar[16] + 6;
   armHHbar[14]= - armHHbar[14] + armHHbar[15] + 2*armHHbar[16];
   armHHbar[14]=armHHbar[14]*MMt;
   armHHbar[13]=armHHbar[13] + armHHbar[10];
   armHHbar[15]=8*MMt;
   armHHbar[13]=armHHbar[13]*armHHbar[15];
   armHHbar[15]=10 + 7*armHHbar[6];
   armHHbar[15]=armHHbar[15]*armHHbar[6];
   armHHbar[13]= - armHHbar[15] + armHHbar[13] - 9 + 4*armHHbar[8];
   armHHbar[13]=armHHbar[13]*MMt;
   armHHbar[15]=armHHbar[6]*armHHbar[7];
   armHHbar[16]= - armHHbar[9] + 3*armHHbar[15];
   armHHbar[13]= - armHHbar[12] + armHHbar[13] - 6*armHHbar[16];
   armHHbar[13]=MMt*armHHbar[13];
   armHHbar[16]=pow(armHHbar[7],2);
   armHHbar[13]= - 12*armHHbar[16] + armHHbar[13];
   armHHbar[13]=armHHbar[3]*armHHbar[13];
   armHHbar[13]=2*armHHbar[13] + armHHbar[14] + 6*armHHbar[15];
   armHHbar[14]=pow(armHHbar[5],2);
   armHHbar[15]=pow(armHHbar[2],2);
   armHHbar[14]=armHHbar[14] + armHHbar[15];

      mHHbarret = 2*armHHbar[14]*armHHbar[13]*armHHbar[4]*armHHbar[1];
      return mHHbarret.real();
}
