#include <tt.hpp>
std::complex<long double>
tt::my10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryutt[22], yuttret;

    aryutt[1]=double(nH);
    aryutt[2]=Tsil::A(MMt,mu2);
    aryutt[3]=pow(CW,-1);
    aryutt[4]=pow(MMH,-1);
    aryutt[5]=pow(MMZ,-1);
    aryutt[6]=pow(SW,-1);
    aryutt[7]=double(boson);
    aryutt[8]=Tsil::B(MMH,MMt,MMt,mu2);
    aryutt[9]=Tsil::B(MMZ,MMt,MMt,mu2);
    aryutt[10]=pow(MMt,-1);
    aryutt[11]=Tsil::A(MMH,mu2);
    aryutt[12]=Tsil::A(MMZ,mu2);
    aryutt[13]=Tsil::A(MMW,mu2);
    aryutt[14]=std::real(Tsil::B(0,MMW,MMt,mu2));
   aryutt[15]= - aryutt[8]*MMH;
   aryutt[15]=1./2.*aryutt[15] + aryutt[11] + 1./2.*aryutt[13];
   aryutt[16]=1./4.*aryutt[14];
   aryutt[17]=aryutt[8] + aryutt[16];
   aryutt[17]=MMt*aryutt[17];
   aryutt[18]=1./2.*aryutt[2];
   aryutt[15]=aryutt[18] + 1./2.*aryutt[15] + aryutt[17];
   aryutt[15]=aryutt[5]*aryutt[15];
   aryutt[17]= - 1 - 7*aryutt[9];
   aryutt[19]=3./2.*aryutt[12] + MMZ;
   aryutt[19]=aryutt[4]*aryutt[19];
   aryutt[20]= - MMZ*aryutt[9];
   aryutt[20]=aryutt[20] - aryutt[12] + aryutt[2];
   aryutt[20]=aryutt[10]*aryutt[20];
   aryutt[17]=aryutt[15] + 17./36.*aryutt[20] + 1./36.*aryutt[17] + 
   aryutt[19];
   aryutt[19]=pow(aryutt[3],2);
   aryutt[17]=aryutt[19]*aryutt[17];
   aryutt[20]= - aryutt[14] - 1./2.*aryutt[9];
   aryutt[20]=MMZ*aryutt[20];
   aryutt[18]=aryutt[20] + aryutt[18] - aryutt[13] - 1./2.*aryutt[12];
   aryutt[18]=aryutt[10]*aryutt[18];
   aryutt[20]=aryutt[9] - 3 + aryutt[14];
   aryutt[21]=MMZ + aryutt[13] + 1./2.*aryutt[12];
   aryutt[21]=aryutt[4]*aryutt[21];
   aryutt[15]=aryutt[15] + 1./2.*aryutt[18] + 1./4.*aryutt[20] + 3*
   aryutt[21];
   aryutt[18]=pow(aryutt[6],2);
   aryutt[15]=aryutt[18]*aryutt[15];
   aryutt[20]=aryutt[12] + 2*aryutt[2];
   aryutt[16]=aryutt[16] + 4./9.*aryutt[9];
   aryutt[16]=MMZ*aryutt[16];
   aryutt[16]=4./9.*aryutt[20] + aryutt[16];
   aryutt[16]=aryutt[10]*aryutt[16];
   aryutt[20]= - 1 + aryutt[9];
   aryutt[21]= - aryutt[4]*MMZ;
   aryutt[15]=1./2.*aryutt[15] + 1./2.*aryutt[17] + aryutt[16] + 8./9.*
   aryutt[20] + aryutt[21];
   aryutt[15]=aryutt[7]*aryutt[15];
   aryutt[16]= - aryutt[5]*aryutt[4]*aryutt[2]*MMt*aryutt[1];
   aryutt[17]=aryutt[19]*aryutt[16];
   aryutt[16]=aryutt[18]*aryutt[16];
   aryutt[16]=aryutt[17] + aryutt[16];

      yuttret = aryutt[15] + 3*aryutt[16];
      return yuttret;
}
