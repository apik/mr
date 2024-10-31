#include <ZZ.hpp>
namespace mr
{
  double ZZ<OS>::y10(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<double> aryuZZ[32], yuZZret;

    aryuZZ[1]=double(nH);
    aryuZZ[2]=double(boson);
    aryuZZ[3]=pow(CW,-1);
    aryuZZ[4]=pow(MMZ,-1);
    aryuZZ[5]=pow(SW,-1);
    aryuZZ[6]=Tsil::B(MMt,MMt,MMZ,mu2);
    aryuZZ[7]=Tsil::B(MMb,MMb,MMZ,mu2);
    aryuZZ[8]=Tsil::B(0,0,MMZ,mu2);
    aryuZZ[9]=Tsil::A(MMt,mu2);
    aryuZZ[10]=Tsil::A(MMb,mu2);
    aryuZZ[11]=double(nL + nH);
    aryuZZ[12]=Tsil::B(MMZ,MMH,MMZ,mu2);
    aryuZZ[13]=Tsil::B(MMW,MMW,MMZ,mu2);
    aryuZZ[14]=Tsil::A(MMH,mu2);
    aryuZZ[15]=Tsil::A(MMZ,mu2);
    aryuZZ[16]=Tsil::A(MMW,mu2);
    aryuZZ[17]=1/( - MMb + MMt);
    aryuZZ[18]=1/( - MMW + MMH);
    aryuZZ[19]=pow(aryuZZ[3],2);
    aryuZZ[20]=pow(aryuZZ[5],2);
    aryuZZ[21]=5./9.*aryuZZ[19] + aryuZZ[20] - 8./9.;
    aryuZZ[22]=1./2.*aryuZZ[20];
    aryuZZ[23]=aryuZZ[22] + 17./18.*aryuZZ[19];
    aryuZZ[24]= - 8./9. - aryuZZ[23];
    aryuZZ[24]=aryuZZ[7]*aryuZZ[24];
    aryuZZ[25]=aryuZZ[19] + aryuZZ[20];
    aryuZZ[26]=MMt - MMb;
    aryuZZ[26]=aryuZZ[10] - aryuZZ[9] - 1./2.*aryuZZ[26];
    aryuZZ[26]=aryuZZ[17]*aryuZZ[25]*aryuZZ[26];
    aryuZZ[24]=3./2.*aryuZZ[26] + aryuZZ[21] + aryuZZ[24];
    aryuZZ[26]=aryuZZ[2]*aryuZZ[1];
    aryuZZ[24]=aryuZZ[24]*aryuZZ[26]*MMb;
    aryuZZ[27]= - 3*aryuZZ[20] + 11./3.;
    aryuZZ[27]=aryuZZ[15]*aryuZZ[27];
    aryuZZ[27]= - aryuZZ[14] + aryuZZ[27];
    aryuZZ[27]=aryuZZ[27]*aryuZZ[20];
    aryuZZ[28]= - aryuZZ[14] + 11./3.*aryuZZ[15];
    aryuZZ[28]=aryuZZ[28]*aryuZZ[19];
    aryuZZ[27]=aryuZZ[27] + aryuZZ[28];
    aryuZZ[28]= - 1 + 3./4.*aryuZZ[20];
    aryuZZ[28]=aryuZZ[28]*aryuZZ[20];
    aryuZZ[29]=5./3.*aryuZZ[19];
    aryuZZ[28]=aryuZZ[29] + 4 + aryuZZ[28];
    aryuZZ[28]=aryuZZ[16]*aryuZZ[28];
    aryuZZ[30]= - aryuZZ[22] - 32./9. + 7./18.*aryuZZ[19];
    aryuZZ[31]=aryuZZ[9]*aryuZZ[1]*aryuZZ[30];
    aryuZZ[27]=aryuZZ[31] + 1./4.*aryuZZ[27] + aryuZZ[28];
    aryuZZ[27]=aryuZZ[2]*aryuZZ[27];
    aryuZZ[28]=aryuZZ[6]*aryuZZ[30];
    aryuZZ[28]=aryuZZ[28] + 41./36.*aryuZZ[19] - 32./9. + 1./4.*
      aryuZZ[20];
    aryuZZ[28]=MMt*aryuZZ[28];
    aryuZZ[21]=aryuZZ[10]*aryuZZ[21];
    aryuZZ[21]=aryuZZ[21] + aryuZZ[28];
    aryuZZ[21]=aryuZZ[26]*aryuZZ[21];
    aryuZZ[28]=aryuZZ[15] - aryuZZ[14];
    aryuZZ[28]= - aryuZZ[2]*aryuZZ[28]*aryuZZ[25];
    aryuZZ[25]=aryuZZ[25]*MMH*aryuZZ[2];
    aryuZZ[30]=aryuZZ[12]*aryuZZ[25];
    aryuZZ[28]=aryuZZ[28] + aryuZZ[30];
    aryuZZ[28]=aryuZZ[4]*MMH*aryuZZ[28];
    aryuZZ[30]=aryuZZ[12] + 1./8.;
    aryuZZ[25]= - aryuZZ[30]*aryuZZ[25];
    aryuZZ[21]=1./12.*aryuZZ[28] + 1./3.*aryuZZ[25] + aryuZZ[27] + 
      aryuZZ[21] + aryuZZ[24];
    aryuZZ[21]=aryuZZ[4]*aryuZZ[21];
    aryuZZ[24]= - 11./9.*aryuZZ[19] + 20./9. - aryuZZ[20];
    aryuZZ[24]=aryuZZ[8]*aryuZZ[24];
    aryuZZ[23]= - 16./9. + aryuZZ[23];
    aryuZZ[23]=aryuZZ[6]*aryuZZ[23];
    aryuZZ[23]=aryuZZ[23] + aryuZZ[24];
    aryuZZ[23]=aryuZZ[1]*aryuZZ[23];
    aryuZZ[24]=pow(CW,2);
    aryuZZ[25]=1./12.*aryuZZ[19] - 33./4.*aryuZZ[20] + 29./3. + 4*
      aryuZZ[24];
    aryuZZ[25]=aryuZZ[13]*aryuZZ[25];
    aryuZZ[27]=aryuZZ[14] - aryuZZ[16];
    aryuZZ[27]=aryuZZ[18]*aryuZZ[27];
    aryuZZ[27]=3./4.*aryuZZ[27] - 209./72. + aryuZZ[12];
    aryuZZ[27]=aryuZZ[20]*aryuZZ[27];
    aryuZZ[28]=5./72. + aryuZZ[12];
    aryuZZ[28]=aryuZZ[28]*aryuZZ[19];
    aryuZZ[20]=aryuZZ[29] - 8./3. + aryuZZ[20];
    aryuZZ[29]=aryuZZ[8] - 1./3.;
    aryuZZ[20]=aryuZZ[11]*aryuZZ[29]*aryuZZ[20];
    aryuZZ[24]=2./3. + aryuZZ[24];
    aryuZZ[20]=aryuZZ[25] + aryuZZ[23] + 4./3.*aryuZZ[20] + aryuZZ[28]
      + 4*aryuZZ[24] + aryuZZ[27];
    aryuZZ[20]=aryuZZ[2]*aryuZZ[20];
    aryuZZ[19]=5./18.*aryuZZ[19] - 4./9. + aryuZZ[22];
    aryuZZ[19]=aryuZZ[7]*aryuZZ[19]*aryuZZ[26];

    yuZZret = aryuZZ[19] + aryuZZ[20] + aryuZZ[21];
    return yuZZret.real();
  }
} // namespace mr
