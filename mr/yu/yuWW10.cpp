#include <WW.hpp>
namespace mr
{
  double WW<OS>::y10(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<double> aryuWW[30], yuWWret;

    aryuWW[1]=double(nH);
    aryuWW[2]=double(boson);
    aryuWW[3]=pow(CW,-1);
    aryuWW[4]=pow(MMZ,-1);
    aryuWW[5]=pow(SW,-1);
    aryuWW[6]=Tsil::B(MMt,MMb,MMW,mu2);
    aryuWW[7]=Tsil::B(0,0,MMW,mu2);
    aryuWW[8]=Tsil::A(MMt,mu2);
    aryuWW[9]=Tsil::A(MMb,mu2);
    aryuWW[10]=double(nL + nH);
    aryuWW[11]=Tsil::B(MMW,MMH,MMW,mu2);
    aryuWW[12]=Tsil::B(MMW,MMZ,MMW,mu2);
    aryuWW[13]=Tsil::A(MMH,mu2);
    aryuWW[14]=Tsil::A(MMZ,mu2);
    aryuWW[15]=Tsil::A(MMW,mu2);
    aryuWW[16]=1/( - MMb + MMt);
    aryuWW[17]=1/( - MMW + MMH);
    aryuWW[18]=aryuWW[11]*pow(MMH,2);
    aryuWW[19]=aryuWW[13]*MMH;
    aryuWW[18]=aryuWW[18] + aryuWW[19];
    aryuWW[19]=aryuWW[9] - aryuWW[8];
    aryuWW[20]=MMb*aryuWW[1];
    aryuWW[19]=aryuWW[20]*aryuWW[19];
    aryuWW[21]=pow(MMb,2);
    aryuWW[22]=aryuWW[6]*aryuWW[1];
    aryuWW[23]=aryuWW[21]*aryuWW[22];
    aryuWW[18]= - aryuWW[19] - aryuWW[23] + 1./6.*aryuWW[18];
    aryuWW[23]=aryuWW[9]*aryuWW[1];
    aryuWW[24]=aryuWW[8]*aryuWW[1];
    aryuWW[24]=aryuWW[23] - aryuWW[24];
    aryuWW[25]=1./2.*MMt;
    aryuWW[26]=aryuWW[25]*aryuWW[22];
    aryuWW[27]=aryuWW[20]*aryuWW[6];
    aryuWW[24]= - aryuWW[26] + aryuWW[27] + 1./2.*aryuWW[24];
    aryuWW[24]=aryuWW[24]*MMt;
    aryuWW[26]=aryuWW[15]*MMH;
    aryuWW[18]=1./2.*aryuWW[18] + aryuWW[24] - 1./12.*aryuWW[26];
    aryuWW[18]=aryuWW[18]*aryuWW[4];
    aryuWW[24]=1./2.*aryuWW[1];
    aryuWW[26]=aryuWW[24] - aryuWW[22];
    aryuWW[26]=aryuWW[26]*aryuWW[25];
    aryuWW[28]=aryuWW[11] + 1./8.;
    aryuWW[29]=1./3.*MMH;
    aryuWW[28]=aryuWW[28]*aryuWW[29];
    aryuWW[29]=aryuWW[24]*aryuWW[8];
    aryuWW[23]= - aryuWW[18] - aryuWW[26] - aryuWW[20] + aryuWW[28] + 1./
      2.*aryuWW[27] + 1./4.*aryuWW[13] + aryuWW[29] - aryuWW[23];
    aryuWW[26]=aryuWW[15] - aryuWW[14];
    aryuWW[18]= - 1./12.*aryuWW[26] + aryuWW[18];
    aryuWW[27]=pow(aryuWW[3],2);
    aryuWW[18]=aryuWW[18]*aryuWW[27];
    aryuWW[18]=aryuWW[18] + 53./12.*aryuWW[15] + 3./2.*aryuWW[14] - 
      aryuWW[23];
    aryuWW[18]=aryuWW[4]*aryuWW[18];
    aryuWW[18]= - 1./24. + aryuWW[18];
    aryuWW[18]=aryuWW[18]*aryuWW[27];
    aryuWW[28]=17 + aryuWW[27];
    aryuWW[28]=aryuWW[28]*aryuWW[27];
    aryuWW[28]=4 + 1./12.*aryuWW[28];
    aryuWW[28]=aryuWW[12]*aryuWW[28];
    aryuWW[18]=aryuWW[28] - 4 + aryuWW[18];
    aryuWW[18]=aryuWW[2]*aryuWW[18];
    aryuWW[23]=5./12.*aryuWW[15] - 1./2.*aryuWW[14] - aryuWW[23];
    aryuWW[23]=aryuWW[4]*aryuWW[23];
    aryuWW[28]=aryuWW[15] - aryuWW[13];
    aryuWW[28]= - 3./4.*aryuWW[28];
    aryuWW[28]=aryuWW[17]*aryuWW[28];
    aryuWW[29]= - aryuWW[1] + 4./3.*aryuWW[10];
    aryuWW[29]=aryuWW[7]*aryuWW[29];
    aryuWW[22]= - 33./4.*aryuWW[12] + aryuWW[23] + aryuWW[22] + 
      aryuWW[29] - 4./9.*aryuWW[10] - 209./72. + aryuWW[11] + aryuWW[28];
    aryuWW[22]=aryuWW[2]*aryuWW[22];
    aryuWW[20]=aryuWW[25]*aryuWW[20];
    aryuWW[21]=aryuWW[21]*aryuWW[24];
    aryuWW[19]= - aryuWW[19] + aryuWW[20] - aryuWW[21];
    aryuWW[20]=aryuWW[2]*aryuWW[4];
    aryuWW[21]=aryuWW[20]*aryuWW[16];
    aryuWW[21]=3./2.*aryuWW[21];
    aryuWW[19]=aryuWW[19]*aryuWW[21];
    aryuWW[21]=pow(aryuWW[5],2);
    aryuWW[20]=aryuWW[26]*aryuWW[21]*aryuWW[20];
    aryuWW[20]=3./4.*aryuWW[20] + aryuWW[22] - aryuWW[19];
    aryuWW[20]=aryuWW[20]*aryuWW[21];
    aryuWW[19]= - aryuWW[27]*aryuWW[19];

    yuWWret = aryuWW[18] + aryuWW[19] + aryuWW[20];
    return yuWWret.real();
  }
} // namespace mr
