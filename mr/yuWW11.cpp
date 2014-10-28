#include <WW.hpp>
std::complex<long double>
WW<OS>::my11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuWW[26], yuWWret;

    aryuWW[1]=double(nH);
    aryuWW[2]=pow(CW,-1);
    aryuWW[3]=pow(MMZ,-1);
    aryuWW[4]=pow(SW,-1);
    aryuWW[5]=Tsil::I2(0,0,MMt,mu2);
    aryuWW[6]=Tsil::B(0,MMt,MMW,mu2);
    aryuWW[7]=Tsil::A(MMt,mu2);
    aryuWW[8]=pow(MMt,-1);
    aryuWW[9]=Tsil::Aeps(MMt,mu2);
    aryuWW[10]=prot00tt0->M(0);
    aryuWW[11]=prot00tt0->Tuxv(0);
    aryuWW[12]=double(nL);
    aryuWW[13]=std::real(Tsil::B(0,0,MMW,mu2));
    aryuWW[14]=prot00000->M(0);
   aryuWW[15]=pow(aryuWW[2],2);
   aryuWW[16]=aryuWW[15] + 1;
   aryuWW[16]=aryuWW[15]*aryuWW[16];
   aryuWW[17]=pow(aryuWW[4],2);
   aryuWW[16]=aryuWW[16] + aryuWW[17];
   aryuWW[16]=aryuWW[16]*pow(aryuWW[3],2);
   aryuWW[18]=aryuWW[16]*MMt;
   aryuWW[15]=aryuWW[15] + aryuWW[17];
   aryuWW[15]=aryuWW[15]*aryuWW[3];
   aryuWW[19]=5*aryuWW[15];
   aryuWW[20]= - aryuWW[19] - 3*aryuWW[18];
   aryuWW[20]=MMt*aryuWW[20];
   aryuWW[21]=aryuWW[15]*aryuWW[7];
   aryuWW[22]=MMt*aryuWW[19];
   aryuWW[22]=4*aryuWW[17] + aryuWW[22];
   aryuWW[22]=aryuWW[6]*aryuWW[22];
   aryuWW[23]=2*aryuWW[17];
   aryuWW[20]=1./3.*aryuWW[22] + 10./3.*aryuWW[21] + aryuWW[23] + 
   aryuWW[20];
   aryuWW[20]=aryuWW[6]*aryuWW[20];
   aryuWW[22]=4./3.*aryuWW[9];
   aryuWW[24]=aryuWW[6]*aryuWW[7];
   aryuWW[24]=aryuWW[22] + 4./3.*aryuWW[24];
   aryuWW[24]=aryuWW[17]*aryuWW[24];
   aryuWW[21]= - 4./3.*aryuWW[17] + aryuWW[21];
   aryuWW[21]=aryuWW[7]*aryuWW[21];
   aryuWW[21]=aryuWW[21] + aryuWW[24];
   aryuWW[21]=aryuWW[8]*aryuWW[21];
   aryuWW[24]=aryuWW[5]*aryuWW[15];
   aryuWW[21]=aryuWW[24] + aryuWW[21];
   aryuWW[24]=aryuWW[17] - 1;
   aryuWW[25]=8./3.*MMZ;
   aryuWW[24]=aryuWW[24]*aryuWW[25];
   aryuWW[25]=aryuWW[16]*pow(MMt,2);
   aryuWW[25]= - aryuWW[17] + 1./3.*aryuWW[25];
   aryuWW[25]=MMt*aryuWW[25];
   aryuWW[25]=aryuWW[24] + 4*aryuWW[25];
   aryuWW[25]=aryuWW[10]*aryuWW[25];
   aryuWW[19]= - aryuWW[19] - aryuWW[18];
   aryuWW[16]=2*aryuWW[16];
   aryuWW[16]=aryuWW[7]*aryuWW[16];
   aryuWW[16]=5*aryuWW[19] + aryuWW[16];
   aryuWW[16]=aryuWW[7]*aryuWW[16];
   aryuWW[19]= - 2*aryuWW[18] - aryuWW[15];
   aryuWW[19]=aryuWW[19]*aryuWW[22];
   aryuWW[18]= - aryuWW[18] + aryuWW[15];
   aryuWW[18]=MMt*aryuWW[18];
   aryuWW[18]=aryuWW[23] + aryuWW[18];
   aryuWW[18]=aryuWW[11]*aryuWW[18];
   aryuWW[15]=aryuWW[15]*MMt;
   aryuWW[15]=8./3.*aryuWW[18] + aryuWW[19] + 2*aryuWW[20] + 2./3.*
   aryuWW[16] + 13*aryuWW[17] + 23./3.*aryuWW[15] + aryuWW[25] + 4*
   aryuWW[21];
   aryuWW[15]=aryuWW[1]*aryuWW[15];
   aryuWW[16]=31./3. + 4*aryuWW[13];
   aryuWW[16]=aryuWW[17]*aryuWW[16];
   aryuWW[17]=aryuWW[14]*aryuWW[24];
   aryuWW[16]=aryuWW[17] + aryuWW[16];
   aryuWW[16]=aryuWW[12]*aryuWW[16];

      yuWWret = aryuWW[15] + aryuWW[16];
      return yuWWret;
}
