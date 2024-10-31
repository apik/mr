#include <tt.hpp>
namespace mr
{
  double tt<OS>::y10(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<double> aryutt[26], yuttret;

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
    aryutt[17]=aryutt[6] + 1./2.*MMt;
    aryutt[18]=1./2.*MMb;
    aryutt[19]=aryutt[14] + aryutt[18] - aryutt[17];
    aryutt[19]=MMb*aryutt[19]*aryutt[15];
    aryutt[17]=aryutt[19] - aryutt[17];
    aryutt[19]=3*aryutt[1];
    aryutt[17]=aryutt[17]*aryutt[19];
    aryutt[19]=aryutt[14] - aryutt[13];
    aryutt[20]=MMb*aryutt[10];
    aryutt[21]=aryutt[19] + aryutt[20];
    aryutt[18]=aryutt[21]*aryutt[18]*aryutt[9];
    aryutt[21]=MMt*aryutt[10];
    aryutt[21]= - aryutt[11] + 3*aryutt[12] + 7*aryutt[13] + aryutt[21]
      + aryutt[14];
    aryutt[17]=aryutt[17] + aryutt[18] - aryutt[20] + aryutt[6] + 1./2.*
      aryutt[21];
    aryutt[18]=aryutt[7]*MMt;
    aryutt[17]=aryutt[18] + 1./2.*aryutt[17];
    aryutt[18]=aryutt[4]*aryutt[17];
    aryutt[21]=aryutt[7] - 1./2.;
    aryutt[22]=MMH*aryutt[4];
    aryutt[22]=1./4.*aryutt[22];
    aryutt[21]=aryutt[21]*aryutt[22];
    aryutt[22]=aryutt[6] - aryutt[12];
    aryutt[23]=aryutt[9]*aryutt[22];
    aryutt[23]=7./2. + 17*aryutt[23];
    aryutt[18]= - aryutt[21] + 1./36.*aryutt[23] + aryutt[18];
    aryutt[23]=pow(aryutt[3],2);
    aryutt[18]=aryutt[18]*aryutt[23];
    aryutt[24]=pow(aryutt[5],2);
    aryutt[25]=aryutt[13] - aryutt[12];
    aryutt[25]=aryutt[25]*aryutt[24];
    aryutt[17]=3./4.*aryutt[25] + aryutt[17];
    aryutt[17]=aryutt[4]*aryutt[17];
    aryutt[20]=aryutt[22] + aryutt[20];
    aryutt[22]=MMZ*aryutt[10];
    aryutt[19]= - aryutt[22] + aryutt[19] + 1./2.*aryutt[20];
    aryutt[19]=aryutt[9]*aryutt[19];
    aryutt[20]= - 3./2. + aryutt[10];
    aryutt[25]= - aryutt[13] + aryutt[11];
    aryutt[25]=aryutt[16]*aryutt[25];
    aryutt[19]=3./2.*aryutt[25] + 1./2.*aryutt[20] + aryutt[19];
    aryutt[17]= - aryutt[21] + 1./2.*aryutt[19] + aryutt[17];
    aryutt[17]=aryutt[17]*aryutt[24];
    aryutt[17]=aryutt[18] + aryutt[17];
    aryutt[18]=aryutt[9]*MMZ;
    aryutt[19]= - 7 - 17*aryutt[18];
    aryutt[19]=aryutt[19]*aryutt[23];
    aryutt[20]=2 + aryutt[18];
    aryutt[19]=4*aryutt[20] + 1./8.*aryutt[19];
    aryutt[18]=1 - aryutt[18];
    aryutt[18]=aryutt[18]*aryutt[24];
    aryutt[18]=1./9.*aryutt[19] + 1./8.*aryutt[18];
    aryutt[18]=aryutt[8]*aryutt[18];
    aryutt[19]=8./9.*aryutt[6] + 4./9.*aryutt[12] + 1./4.*aryutt[22];
    aryutt[19]=aryutt[9]*aryutt[19];
    aryutt[17]=aryutt[18] - 8./9. + aryutt[19] + 1./2.*aryutt[17];

    yuttret = aryutt[17]*aryutt[2];
    return yuttret.real();
  }
} // namespace mr
