#include <HH.hpp>
std::complex<long double>
HH<OS>::m11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armHHbar[20], mHHbarret;

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
   armHHbar[14]=armHHbar[13] + armHHbar[10];
   armHHbar[15]=8*MMt;
   armHHbar[14]=armHHbar[14]*armHHbar[15];
   armHHbar[14]=4*armHHbar[8] + armHHbar[14] - 9;
   armHHbar[15]=pow(MMt,2);
   armHHbar[14]=armHHbar[15]*armHHbar[14];
   armHHbar[16]=6*MMt;
   armHHbar[17]=armHHbar[16]*armHHbar[9];
   armHHbar[18]=armHHbar[12]*MMt;
   armHHbar[19]=pow(armHHbar[7],2);
   armHHbar[14]= - armHHbar[17] - armHHbar[14] + armHHbar[18] + 12*
   armHHbar[19];
   armHHbar[17]=armHHbar[3]*armHHbar[1];
   armHHbar[14]=armHHbar[17]*armHHbar[14];
   armHHbar[13]=armHHbar[13] + 3*armHHbar[10];
   armHHbar[13]=armHHbar[13]*MMt;
   armHHbar[13]=armHHbar[13] - 3;
   armHHbar[18]=MMH*armHHbar[10];
   armHHbar[13]= - armHHbar[18] + 2*armHHbar[13];
   armHHbar[18]=armHHbar[1]*MMt;
   armHHbar[13]=armHHbar[13]*armHHbar[18];
   armHHbar[13]=armHHbar[13] + armHHbar[14];
   armHHbar[14]=armHHbar[15]*armHHbar[17];
   armHHbar[15]= - armHHbar[18] + 5*armHHbar[14];
   armHHbar[16]=armHHbar[16]*armHHbar[17];
   armHHbar[16]=armHHbar[16] - armHHbar[1];
   armHHbar[17]=3*armHHbar[7];
   armHHbar[16]=armHHbar[16]*armHHbar[17];
   armHHbar[15]=armHHbar[16] + 2*armHHbar[15];
   armHHbar[14]= - 3*armHHbar[18] + 14*armHHbar[14];
   armHHbar[14]=armHHbar[14]*armHHbar[6];
   armHHbar[14]=armHHbar[14] + 2*armHHbar[15];
   armHHbar[14]=armHHbar[14]*armHHbar[6];
   armHHbar[13]=armHHbar[14] + 2*armHHbar[13];
   armHHbar[14]= - pow(armHHbar[5],2);
   armHHbar[15]= - pow(armHHbar[2],2);
   armHHbar[14]=armHHbar[14] + armHHbar[15];

      mHHbarret = 2*armHHbar[14]*armHHbar[13]*armHHbar[4];
      return mHHbarret;
}
