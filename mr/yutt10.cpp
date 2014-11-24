#include <tt.hpp>
std::complex<long double>
tt::y10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryutt[25], yuttret;

    aryutt[1]=double(nH);
    aryutt[2]=double(boson);
    aryutt[3]=pow(CW,-1);
    aryutt[4]=pow(MMZ,-1);
    aryutt[5]=pow(SW,-1);
    aryutt[6]=Tsil::A(MMt,mu2);
    aryutt[7]=Tsil::B(MMH,MMt,MMt,mu2);
    aryutt[8]=Tsil::B(MMZ,MMt,MMt,mu2);
    aryutt[9]=pow(MMt,-1);
    aryutt[10]=Tsil::B(MMW,MMb,MMt,mu2);
    aryutt[11]=Tsil::A(MMH,mu2);
    aryutt[12]=Tsil::A(MMZ,mu2);
    aryutt[13]=Tsil::A(MMW,mu2);
    aryutt[14]=Tsil::A(MMb,mu2);
    aryutt[15]=1/( - MMb + MMt);
    aryutt[16]=1/( - MMW + MMH);
   aryutt[17]=1./2.*MMt;
   aryutt[18]=aryutt[17] + aryutt[6];
   aryutt[19]=1./2.*MMb;
   aryutt[20]=aryutt[19] - aryutt[18];
   aryutt[20]=aryutt[20]*MMb;
   aryutt[21]=aryutt[14]*MMb;
   aryutt[20]=aryutt[20] + aryutt[21];
   aryutt[20]=aryutt[20]*aryutt[15];
   aryutt[18]=aryutt[20] - aryutt[18];
   aryutt[20]=3*aryutt[1];
   aryutt[18]=aryutt[18]*aryutt[20];
   aryutt[20]=aryutt[11] - aryutt[14];
   aryutt[17]=aryutt[17] - MMb;
   aryutt[17]=aryutt[17]*aryutt[10];
   aryutt[22]=1./4.*MMH;
   aryutt[17]=aryutt[18] + aryutt[22] + aryutt[6] + 7./2.*aryutt[13] + 
   aryutt[17] + 3./2.*aryutt[12] - 1./2.*aryutt[20];
   aryutt[18]=pow(aryutt[5],2);
   aryutt[20]=aryutt[13] - aryutt[12];
   aryutt[20]=aryutt[20]*aryutt[18];
   aryutt[20]=3./2.*aryutt[20] + aryutt[17];
   aryutt[23]=aryutt[10]*pow(MMb,2);
   aryutt[24]=aryutt[13]*MMb;
   aryutt[21]= - aryutt[23] + aryutt[24] - aryutt[21];
   aryutt[23]=1./4.*aryutt[9];
   aryutt[21]=aryutt[21]*aryutt[23];
   aryutt[22]=aryutt[22] - MMt;
   aryutt[22]=aryutt[22]*aryutt[7];
   aryutt[21]=aryutt[21] + aryutt[22];
   aryutt[20]=1./2.*aryutt[20] - aryutt[21];
   aryutt[20]=aryutt[4]*aryutt[18]*aryutt[20];
   aryutt[17]=1./2.*aryutt[17] - aryutt[21];
   aryutt[17]=aryutt[4]*aryutt[17];
   aryutt[21]=1./2. - aryutt[8];
   aryutt[22]=aryutt[8]*MMZ;
   aryutt[22]=aryutt[22] + aryutt[12];
   aryutt[23]=aryutt[22] - aryutt[6];
   aryutt[24]= - aryutt[9]*aryutt[23];
   aryutt[21]=7*aryutt[21] + 17*aryutt[24];
   aryutt[17]=1./36.*aryutt[21] + aryutt[17];
   aryutt[17]=aryutt[17]*pow(aryutt[3],2);
   aryutt[17]=aryutt[20] + aryutt[17];
   aryutt[20]=aryutt[13] - aryutt[11];
   aryutt[20]= - 3*aryutt[20];
   aryutt[20]=aryutt[16]*aryutt[20];
   aryutt[20]=aryutt[8] - 3./2. + aryutt[10] + aryutt[20];
   aryutt[20]=aryutt[20]*aryutt[18];
   aryutt[19]= - MMZ + aryutt[19];
   aryutt[19]=aryutt[10]*aryutt[19];
   aryutt[19]= - aryutt[13] + aryutt[14] + aryutt[19] - 1./2.*
   aryutt[23];
   aryutt[18]=aryutt[19]*aryutt[18];
   aryutt[19]=aryutt[10]*MMZ;
   aryutt[18]=aryutt[19] + aryutt[18];
   aryutt[18]=8./9.*aryutt[6] + 4./9.*aryutt[22] + 1./4.*aryutt[18];
   aryutt[18]=aryutt[9]*aryutt[18];
   aryutt[19]= - 1 + aryutt[8];
   aryutt[17]=aryutt[18] + 8./9.*aryutt[19] + 1./8.*aryutt[20] + 1./2.*
   aryutt[17];

      yuttret = aryutt[17]*aryutt[2];
      return yuttret;
}
