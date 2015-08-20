#include <ZZ.hpp>
namespace mr
{
  long double ZZ<OS>::y11(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<long double> aryuZZ[32], yuZZret;

    aryuZZ[1]=double(nH);
    aryuZZ[2]=pow(CW,-1);
    aryuZZ[3]=pow(MMZ,-1);
    aryuZZ[4]=pow(SW,-1);
    aryuZZ[5]=Tsil::I2(0,0,MMt,mu2);
    aryuZZ[6]=Tsil::B(MMt,MMt,MMZ,mu2);
    aryuZZ[7]=Tsil::A(MMt,mu2);
    aryuZZ[8]=Tsil::Beps(MMt,MMt,MMZ,mu2);
    aryuZZ[9]=pow(MMt,-1);
    aryuZZ[10]=Tsil::Aeps(MMt,mu2);
    aryuZZ[11]=std::real(Tsil::B(0,0,MMZ,mu2));
    aryuZZ[12]=prottttt0->M(0);
    aryuZZ[13]=prot00000->M(0);
    aryuZZ[14]=prottttt0->Vzxyv(0);
    aryuZZ[15]=prottttt0->Suxv(0);
    aryuZZ[16]=double(nL);
    aryuZZ[17]=1/(4*MMt - MMZ);
    aryuZZ[18]=pow(aryuZZ[4],2);
    aryuZZ[19]=pow(aryuZZ[2],2);
    aryuZZ[20]=aryuZZ[18] + aryuZZ[19];
    aryuZZ[21]=aryuZZ[20]*pow(aryuZZ[3],2);
    aryuZZ[22]=7./9.*aryuZZ[3];
    aryuZZ[22]=aryuZZ[22]*aryuZZ[19];
    aryuZZ[23]=aryuZZ[18] + 64./9.;
    aryuZZ[23]=aryuZZ[23]*aryuZZ[3];
    aryuZZ[22]=aryuZZ[23] - aryuZZ[22];
    aryuZZ[23]=aryuZZ[12]*aryuZZ[22];
    aryuZZ[23]=2*aryuZZ[23] - aryuZZ[21];
    aryuZZ[24]=4./3.*MMt;
    aryuZZ[23]=aryuZZ[23]*aryuZZ[24];
    aryuZZ[25]=17./9.*aryuZZ[19] + aryuZZ[18] - 32./9.;
    aryuZZ[26]=MMt*aryuZZ[22];
    aryuZZ[26]=aryuZZ[26] - aryuZZ[25];
    aryuZZ[26]=aryuZZ[14]*aryuZZ[26];
    aryuZZ[27]= - 197./18.*aryuZZ[19] - 128./9. - 29./2.*aryuZZ[18];
    aryuZZ[27]=aryuZZ[3]*aryuZZ[27];
    aryuZZ[28]=aryuZZ[12]*aryuZZ[20];
    aryuZZ[23]=16./3.*aryuZZ[26] + aryuZZ[23] + 1./3.*aryuZZ[27] - 4*
      aryuZZ[28];
    aryuZZ[23]=MMt*aryuZZ[23];
    aryuZZ[26]=aryuZZ[25]*aryuZZ[3];
    aryuZZ[27]=aryuZZ[21]*MMt;
    aryuZZ[28]= - 1./3.*aryuZZ[27] - aryuZZ[26];
    aryuZZ[28]=aryuZZ[15]*aryuZZ[28];
    aryuZZ[26]=2*aryuZZ[26] + aryuZZ[27];
    aryuZZ[29]=25./9.*aryuZZ[19] + aryuZZ[18] - 64./9.;
    aryuZZ[30]=aryuZZ[29]*aryuZZ[17];
    aryuZZ[26]=2./3.*aryuZZ[26] - aryuZZ[30];
    aryuZZ[26]=aryuZZ[10]*aryuZZ[26];
    aryuZZ[20]=aryuZZ[20]*aryuZZ[3];
    aryuZZ[31]=aryuZZ[5]*aryuZZ[20];
    aryuZZ[26]=aryuZZ[31] + aryuZZ[28] + aryuZZ[26];
    aryuZZ[28]=2*MMt;
    aryuZZ[21]=aryuZZ[21]*aryuZZ[28];
    aryuZZ[28]=29./3.*aryuZZ[19] - 128./3. - aryuZZ[18];
    aryuZZ[28]=aryuZZ[3]*aryuZZ[28];
    aryuZZ[28]=aryuZZ[21] + aryuZZ[28];
    aryuZZ[28]=aryuZZ[28]*aryuZZ[24];
    aryuZZ[21]=aryuZZ[21] - aryuZZ[22];
    aryuZZ[21]=MMt*aryuZZ[21];
    aryuZZ[21]=2./3.*aryuZZ[21] - aryuZZ[29];
    aryuZZ[21]=aryuZZ[6]*aryuZZ[21];
    aryuZZ[21]=aryuZZ[21] + aryuZZ[28] + 143./9.*aryuZZ[19] - 320./9. + 
      7*aryuZZ[18];
    aryuZZ[21]=aryuZZ[6]*aryuZZ[21];
    aryuZZ[28]=aryuZZ[27] - aryuZZ[22];
    aryuZZ[28]=2./3.*aryuZZ[28] - aryuZZ[30];
    aryuZZ[28]=aryuZZ[6]*aryuZZ[28];
    aryuZZ[28]=aryuZZ[30] + aryuZZ[28];
    aryuZZ[31]=41./9.*aryuZZ[19] - 128./9. + aryuZZ[18];
    aryuZZ[31]=aryuZZ[3]*aryuZZ[31];
    aryuZZ[27]=7*aryuZZ[31] + 8*aryuZZ[27];
    aryuZZ[31]=2*aryuZZ[7];
    aryuZZ[20]= - aryuZZ[31]*aryuZZ[9]*aryuZZ[20];
    aryuZZ[20]=aryuZZ[20] + 1./3.*aryuZZ[27] + 2*aryuZZ[28];
    aryuZZ[20]=aryuZZ[31]*aryuZZ[20];
    aryuZZ[22]=aryuZZ[22]*aryuZZ[24];
    aryuZZ[22]=aryuZZ[22] - aryuZZ[29];
    aryuZZ[22]=aryuZZ[8]*aryuZZ[22];
    aryuZZ[24]=937./18.*aryuZZ[19] - 860./9. + 77./2.*aryuZZ[18];
    aryuZZ[20]=aryuZZ[20] + aryuZZ[21] + aryuZZ[22] + 1./3.*aryuZZ[24]
      + aryuZZ[23] + 4*aryuZZ[26];
    aryuZZ[20]=aryuZZ[1]*aryuZZ[20];
    aryuZZ[21]=5./9.*aryuZZ[19] + aryuZZ[18] - 8./9.;
    aryuZZ[22]=aryuZZ[13]*aryuZZ[21];
    aryuZZ[23]=aryuZZ[12]*aryuZZ[25];
    aryuZZ[22]=aryuZZ[23] + aryuZZ[22];
    aryuZZ[23]= - aryuZZ[6]*aryuZZ[29];
    aryuZZ[23]=aryuZZ[23] + 25./3.*aryuZZ[19] - 64./3. + 3*aryuZZ[18];
    aryuZZ[23]=aryuZZ[6]*aryuZZ[17]*aryuZZ[23];
    aryuZZ[24]= - aryuZZ[8]*aryuZZ[30];
    aryuZZ[22]=aryuZZ[23] + 4./3.*aryuZZ[22] + aryuZZ[24];
    aryuZZ[22]=aryuZZ[1]*aryuZZ[22];
    aryuZZ[18]=11./9.*aryuZZ[19] + aryuZZ[18] - 20./9.;
    aryuZZ[18]=aryuZZ[18]*aryuZZ[16];
    aryuZZ[19]=aryuZZ[13]*aryuZZ[18];
    aryuZZ[19]=8./3.*aryuZZ[19] + aryuZZ[22];
    aryuZZ[19]=MMZ*aryuZZ[19];
    aryuZZ[21]=aryuZZ[1]*aryuZZ[21];
    aryuZZ[21]=2*aryuZZ[18] + aryuZZ[21];
    aryuZZ[21]=aryuZZ[11]*aryuZZ[21];

    yuZZret = 31./3.*aryuZZ[18] + aryuZZ[19] + aryuZZ[20] + 2*
      aryuZZ[21];
    return yuZZret.real();
  }
} // namespace mr
