#include <WW.hpp>
std::complex<long double>
WW<OS>::my11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuWW[24], yuWWret;

    aryuWW[1]=double(nH);
    aryuWW[2]=pow(CW,-1);
    aryuWW[3]=pow(MMH,-1);
    aryuWW[4]=pow(MMZ,-1);
    aryuWW[5]=pow(SW,-1);
    aryuWW[6]=Tsil::B(0,MMt,MMW,mu2);
    aryuWW[7]=Tsil::A(MMt,mu2);
    aryuWW[8]=pow(MMt,-1);
    aryuWW[9]=Tsil::Aeps(MMt,mu2);
    aryuWW[10]=prot00tt0->M(0);
    aryuWW[11]=prot00tt0->Tuxv(0);
    aryuWW[12]=double(nL);
    aryuWW[13]=std::real(Tsil::B(0,0,MMW,mu2));
    aryuWW[14]=prot00000->M(0);
   aryuWW[15]=aryuWW[10]*MMZ;
   aryuWW[16]=aryuWW[9]*aryuWW[8];
   aryuWW[17]=1 + 2./3.*aryuWW[6];
   aryuWW[17]=aryuWW[6]*aryuWW[17];
   aryuWW[18]=aryuWW[6]*aryuWW[8];
   aryuWW[18]= - aryuWW[8] + aryuWW[18];
   aryuWW[18]=aryuWW[7]*aryuWW[18];
   aryuWW[19]= - MMt*aryuWW[10];
   aryuWW[15]=4*aryuWW[19] + 16./3.*aryuWW[18] + 4*aryuWW[17] + 16./3.*
   aryuWW[16] + 8./3.*aryuWW[15] + 13 + 16./3.*aryuWW[11];
   aryuWW[16]=pow(aryuWW[5],2);
   aryuWW[15]=aryuWW[16]*aryuWW[15];
   aryuWW[17]= - 4./3.*aryuWW[11] - 3*aryuWW[6];
   aryuWW[18]=pow(aryuWW[2],2);
   aryuWW[19]=aryuWW[18]*aryuWW[17];
   aryuWW[19]=aryuWW[17] + aryuWW[19];
   aryuWW[19]=aryuWW[18]*aryuWW[19];
   aryuWW[20]=aryuWW[18]*aryuWW[10];
   aryuWW[20]=aryuWW[10] + aryuWW[20];
   aryuWW[20]=MMt*aryuWW[18]*aryuWW[20];
   aryuWW[19]=aryuWW[19] + 2./3.*aryuWW[20];
   aryuWW[19]=MMt*aryuWW[19];
   aryuWW[20]= - 4*aryuWW[9] - 5*aryuWW[7];
   aryuWW[21]=aryuWW[18]*aryuWW[20];
   aryuWW[21]=aryuWW[20] + aryuWW[21];
   aryuWW[21]=aryuWW[18]*aryuWW[21];
   aryuWW[19]=1./3.*aryuWW[21] + aryuWW[19];
   aryuWW[19]=MMt*aryuWW[19];
   aryuWW[21]=pow(aryuWW[7],2);
   aryuWW[22]=aryuWW[18]*aryuWW[21];
   aryuWW[22]=aryuWW[21] + aryuWW[22];
   aryuWW[22]=aryuWW[18]*aryuWW[22];
   aryuWW[23]=MMt*aryuWW[10];
   aryuWW[17]=aryuWW[17] + 2./3.*aryuWW[23];
   aryuWW[17]=MMt*aryuWW[17];
   aryuWW[17]=1./3.*aryuWW[20] + aryuWW[17];
   aryuWW[17]=MMt*aryuWW[17];
   aryuWW[17]=2./3.*aryuWW[21] + aryuWW[17];
   aryuWW[17]=aryuWW[16]*aryuWW[17];
   aryuWW[17]=aryuWW[17] + 2./3.*aryuWW[22] + aryuWW[19];
   aryuWW[17]=aryuWW[4]*aryuWW[17];
   aryuWW[19]= - 1 + aryuWW[6];
   aryuWW[20]= - 6*aryuWW[3] + aryuWW[8];
   aryuWW[20]=aryuWW[7]*aryuWW[20];
   aryuWW[19]=5./3.*aryuWW[19] + 2*aryuWW[20];
   aryuWW[19]=aryuWW[7]*aryuWW[19];
   aryuWW[19]=2./3.*aryuWW[9] + aryuWW[19];
   aryuWW[20]= - 65./2. + 8*aryuWW[11];
   aryuWW[21]= - 1 + 1./3.*aryuWW[6];
   aryuWW[21]=aryuWW[6]*aryuWW[21];
   aryuWW[20]=1./3.*aryuWW[20] + 10*aryuWW[21];
   aryuWW[21]=MMt*aryuWW[3];
   aryuWW[21]=aryuWW[20] + 64*aryuWW[21];
   aryuWW[21]=MMt*aryuWW[21];
   aryuWW[21]=4*aryuWW[19] + aryuWW[21];
   aryuWW[21]=aryuWW[16]*aryuWW[21];
   aryuWW[19]=aryuWW[18]*aryuWW[19];
   aryuWW[20]=aryuWW[18]*aryuWW[20];
   aryuWW[18]=MMt*aryuWW[18]*aryuWW[3];
   aryuWW[18]=aryuWW[20] + 64*aryuWW[18];
   aryuWW[18]=MMt*aryuWW[18];
   aryuWW[17]=2*aryuWW[17] + aryuWW[21] + 4*aryuWW[19] + aryuWW[18];
   aryuWW[17]=aryuWW[4]*aryuWW[17];
   aryuWW[18]= - aryuWW[10]*MMZ;
   aryuWW[15]=aryuWW[17] + 8./3.*aryuWW[18] + aryuWW[15];
   aryuWW[15]=aryuWW[1]*aryuWW[15];
   aryuWW[17]=MMZ*aryuWW[14];
   aryuWW[17]=8./3.*aryuWW[17] + 31./3. + 4*aryuWW[13];
   aryuWW[16]=aryuWW[16]*aryuWW[12]*aryuWW[17];
   aryuWW[17]= - aryuWW[12]*MMZ*aryuWW[14];

      yuWWret = aryuWW[15] + aryuWW[16] + 8./3.*aryuWW[17];
      return yuWWret;
}
