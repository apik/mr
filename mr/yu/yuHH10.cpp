#include <HH.hpp>
namespace mr
{
  double HH<OS>::y10(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<double> aryuHH[29], yuHHret;

    aryuHH[1]=double(nH);
    aryuHH[2]=double(boson);
    aryuHH[3]=pow(CW,-1);
    aryuHH[4]=pow(MMZ,-1);
    aryuHH[5]=pow(SW,-1);
    aryuHH[6]=Tsil::B(MMt,MMt,MMH,mu2);
    aryuHH[7]=pow(MMH,-1);
    aryuHH[8]=Tsil::B(MMb,MMb,MMH,mu2);
    aryuHH[9]=Tsil::A(MMt,mu2);
    aryuHH[10]=Tsil::B(MMH,MMH,MMH,mu2);
    aryuHH[11]=Tsil::B(MMZ,MMZ,MMH,mu2);
    aryuHH[12]=Tsil::B(MMW,MMW,MMH,mu2);
    aryuHH[13]=Tsil::A(MMZ,mu2);
    aryuHH[14]=Tsil::A(MMW,mu2);
    aryuHH[15]=1/( - MMb + MMt);
    aryuHH[16]=Tsil::A(MMb,mu2);
    aryuHH[17]=1/( - MMW + MMH);
    aryuHH[18]=Tsil::A(MMH,mu2);
    aryuHH[19]= - MMt + MMb;
    aryuHH[20]=1./2.*aryuHH[1];
    aryuHH[19]=aryuHH[20]*aryuHH[19];
    aryuHH[21]=aryuHH[16]*aryuHH[1];
    aryuHH[22]=aryuHH[9]*aryuHH[1];
    aryuHH[19]=aryuHH[21] - aryuHH[22] + aryuHH[19];
    aryuHH[21]=3./2.*aryuHH[15];
    aryuHH[19]=aryuHH[19]*aryuHH[21];
    aryuHH[21]=2*aryuHH[7];
    aryuHH[23]=aryuHH[21]*MMb*aryuHH[1];
    aryuHH[23]=aryuHH[23] - aryuHH[20];
    aryuHH[23]=aryuHH[8]*aryuHH[23];
    aryuHH[19]=aryuHH[19] - 3*aryuHH[23];
    aryuHH[19]=aryuHH[19]*MMb;
    aryuHH[23]=aryuHH[11] + 1;
    aryuHH[24]=1./8.*MMH;
    aryuHH[23]=aryuHH[24]*aryuHH[23];
    aryuHH[24]=MMt*aryuHH[1];
    aryuHH[21]=aryuHH[21]*aryuHH[24];
    aryuHH[20]=aryuHH[21] - aryuHH[20];
    aryuHH[21]=aryuHH[6]*MMt;
    aryuHH[21]=3*aryuHH[21];
    aryuHH[20]=aryuHH[20]*aryuHH[21];
    aryuHH[21]=9./8.*aryuHH[10] + 1./4.*aryuHH[12];
    aryuHH[21]=aryuHH[21]*MMH;
    aryuHH[19]=aryuHH[19] + aryuHH[21] - 3./2.*aryuHH[22] + aryuHH[13]
      - aryuHH[20] - 3./4.*aryuHH[24] + 2*aryuHH[14] + aryuHH[23];
    aryuHH[19]=aryuHH[19]*aryuHH[4];
    aryuHH[20]=1./2.*aryuHH[11];
    aryuHH[21]=aryuHH[12] + aryuHH[20];
    aryuHH[22]=MMZ*aryuHH[7];
    aryuHH[23]=3*aryuHH[22];
    aryuHH[24]=aryuHH[23] - 1;
    aryuHH[21]=aryuHH[24]*aryuHH[21];
    aryuHH[25]=pow(aryuHH[5],2);
    aryuHH[26]= - aryuHH[13] + aryuHH[14];
    aryuHH[26]=aryuHH[4]*aryuHH[26]*aryuHH[25];
    aryuHH[27]=aryuHH[22] - 1./8.;
    aryuHH[28]=aryuHH[18] - aryuHH[14];
    aryuHH[28]=3./4.*aryuHH[28];
    aryuHH[28]=aryuHH[17]*aryuHH[28];
    aryuHH[21]=3./4.*aryuHH[26] + aryuHH[19] - 3*aryuHH[27] + aryuHH[28]
      + aryuHH[21];
    aryuHH[21]=aryuHH[25]*aryuHH[21];
    aryuHH[20]=aryuHH[24]*aryuHH[20];
    aryuHH[19]=aryuHH[19] + aryuHH[20] - aryuHH[27];
    aryuHH[19]=aryuHH[19]*pow(aryuHH[3],2);
    aryuHH[20]= - aryuHH[12]*aryuHH[23];
    aryuHH[19]=aryuHH[21] + aryuHH[19] + 2*aryuHH[22] + aryuHH[20];

    yuHHret = aryuHH[19]*aryuHH[2];
    return yuHHret.real();
  }
} // namespace mr
