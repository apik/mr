#include <WW.hpp>
std::complex<long double>
WW<OS>::m11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armWWbar[26], mWWbarret;

    armWWbar[1]=double(nH);
    armWWbar[2]=pow(CW,-1);
    armWWbar[3]=pow(MMH,-1);
    armWWbar[4]=pow(MMZ,-1);
    armWWbar[5]=pow(SW,-1);
    armWWbar[6]=Tsil::B(0,MMt,MMW,mu2);
    armWWbar[7]=Tsil::A(MMt,mu2);
    armWWbar[8]=pow(MMt,-1);
    armWWbar[9]=Tsil::Aeps(MMt,mu2);
    armWWbar[10]=prot00tt0->M(0);
    armWWbar[11]=prot00tt0->Tuxv(0);
    armWWbar[12]=double(nL);
    armWWbar[13]=std::real(Tsil::B(0,0,MMW,mu2));
    armWWbar[14]=prot00000->M(0);
   armWWbar[15]= - 1 + 1./3.*armWWbar[6];
   armWWbar[16]=2*armWWbar[6];
   armWWbar[15]=armWWbar[15]*armWWbar[16];
   armWWbar[15]=armWWbar[15] - 13./6.;
   armWWbar[16]=MMt*armWWbar[3];
   armWWbar[15]=64*armWWbar[16] + 8./3.*armWWbar[11] + 5*armWWbar[15];
   armWWbar[15]=armWWbar[15]*MMt;
   armWWbar[16]=armWWbar[6] - 1;
   armWWbar[17]=armWWbar[7]*armWWbar[8];
   armWWbar[18]=5./3.*armWWbar[16] + 2*armWWbar[17];
   armWWbar[18]=armWWbar[18]*armWWbar[7];
   armWWbar[19]=pow(armWWbar[7],2);
   armWWbar[20]=armWWbar[19]*armWWbar[3];
   armWWbar[18]=armWWbar[18] - 12*armWWbar[20];
   armWWbar[15]=armWWbar[15] + 4*armWWbar[18];
   armWWbar[18]=pow(armWWbar[2],2);
   armWWbar[20]=armWWbar[15]*armWWbar[18];
   armWWbar[21]=MMt*armWWbar[10];
   armWWbar[22]= - 2./3.*armWWbar[21] + 4./3.*armWWbar[11] + 3*
   armWWbar[6];
   armWWbar[22]=armWWbar[22]*MMt;
   armWWbar[22]=armWWbar[22] + 5./3.*armWWbar[7];
   armWWbar[23]=armWWbar[22]*MMt;
   armWWbar[19]=2./3.*armWWbar[19];
   armWWbar[23]=armWWbar[23] - armWWbar[19];
   armWWbar[23]= - armWWbar[18]*armWWbar[23];
   armWWbar[24]=armWWbar[18]*armWWbar[9];
   armWWbar[25]= - armWWbar[24]*MMt;
   armWWbar[23]=armWWbar[23] + 4./3.*armWWbar[25];
   armWWbar[18]=armWWbar[18] + 1;
   armWWbar[25]=2*armWWbar[4];
   armWWbar[18]=armWWbar[25]*armWWbar[18]*armWWbar[23];
   armWWbar[18]=armWWbar[18] + armWWbar[20] + 8./3.*armWWbar[24];
   armWWbar[18]=armWWbar[4]*armWWbar[18];
   armWWbar[20]= - 4./3.*armWWbar[9] - armWWbar[22];
   armWWbar[20]=MMt*armWWbar[20];
   armWWbar[19]=armWWbar[19] + armWWbar[20];
   armWWbar[19]=armWWbar[19]*armWWbar[25];
   armWWbar[15]=armWWbar[19] + 8./3.*armWWbar[9] + armWWbar[15];
   armWWbar[15]=armWWbar[4]*armWWbar[15];
   armWWbar[16]=armWWbar[16]*armWWbar[17];
   armWWbar[17]=armWWbar[9]*armWWbar[8];
   armWWbar[16]=armWWbar[11] + armWWbar[16] + armWWbar[17];
   armWWbar[17]=1 + 2./3.*armWWbar[6];
   armWWbar[17]=armWWbar[6]*armWWbar[17];
   armWWbar[17]=armWWbar[21] - armWWbar[17];
   armWWbar[19]=8./3.*MMZ;
   armWWbar[20]=armWWbar[19]*armWWbar[10];
   armWWbar[15]=armWWbar[15] + 13 + armWWbar[20] - 4*armWWbar[17] + 16./
   3.*armWWbar[16];
   armWWbar[16]=pow(armWWbar[5],2);
   armWWbar[15]=armWWbar[15]*armWWbar[16];
   armWWbar[15]=armWWbar[15] - armWWbar[20] + armWWbar[18];
   armWWbar[15]=armWWbar[1]*armWWbar[15];
   armWWbar[17]=armWWbar[19]*armWWbar[14];
   armWWbar[18]=armWWbar[17] + 31./3. + 4*armWWbar[13];
   armWWbar[16]=armWWbar[18]*armWWbar[16];
   armWWbar[16]= - armWWbar[17] + armWWbar[16];
   armWWbar[16]=armWWbar[12]*armWWbar[16];

      mWWbarret = armWWbar[15] + armWWbar[16];
      return mWWbarret;
}
