#include <bb.hpp>
namespace mr
{
  double bb<OS>::y11(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<double> aryubb[43], yubbret;

    aryubb[1]=double(nH);
    aryubb[2]=pow(CW,-1);
    aryubb[3]=pow(MMZ,-1);
    aryubb[4]=pow(SW,-1);
    aryubb[5]=Tsil::I2(0,0,MMt,mu2);
    aryubb[6]=Tsil::A(MMt,mu2);
    aryubb[7]=pow(MMt,-1);
    aryubb[8]=Tsil::A(MMb,mu2);
    aryubb[9]=pow(MMb,-1);
    aryubb[10]=Tsil::Aeps(MMt,mu2);
    aryubb[11]=double(boson);
    aryubb[12]=Tsil::I2(0,MMW,MMt,mu2);
    aryubb[13]=Tsil::A(MMZ,mu2);
    aryubb[14]=Tsil::A(MMW,mu2);
    aryubb[15]=Tsil::Aeps(MMW,mu2);
    aryubb[16]=Tsil::Aeps(MMb,mu2);
    aryubb[17]=prot0bb0b->Suxv(0);
    aryubb[18]=1/(MMt - MMW);
    aryubb[19]=1/( - MMW + MMH);
    aryubb[20]=Tsil::A(MMH,mu2);
    aryubb[21]=pow(aryubb[2],2);
    aryubb[22]=pow(aryubb[4],2);
    aryubb[23]=aryubb[21] + aryubb[22];
    aryubb[24]=aryubb[23]*aryubb[3];
    aryubb[25]=aryubb[24]*MMt;
    aryubb[26]=MMH*aryubb[24];
    aryubb[27]= - 26 - 11./2.*aryubb[21];
    aryubb[28]=aryubb[20]*aryubb[22]*aryubb[19];
    aryubb[27]=1./4.*aryubb[25] + 3./2.*aryubb[28] + 1./9.*aryubb[27] + 
      1./4.*aryubb[26];
    aryubb[27]=aryubb[9]*aryubb[27];
    aryubb[28]=aryubb[8]*pow(aryubb[9],2);
    aryubb[27]=32./9.*aryubb[28] + aryubb[27];
    aryubb[27]=aryubb[8]*aryubb[27];
    aryubb[28]=aryubb[22] - 1;
    aryubb[29]=pow(aryubb[18],2);
    aryubb[30]=aryubb[28]*aryubb[29];
    aryubb[31]=pow(CW,2);
    aryubb[32]=aryubb[29]*aryubb[31];
    aryubb[32]= - aryubb[30] + aryubb[32];
    aryubb[33]=pow(aryubb[18],3);
    aryubb[34]=aryubb[28]*aryubb[33];
    aryubb[35]=aryubb[31]*aryubb[33];
    aryubb[33]=aryubb[33] + aryubb[35];
    aryubb[31]=aryubb[33]*aryubb[31];
    aryubb[31]= - aryubb[34] + aryubb[31];
    aryubb[31]=MMZ*aryubb[31];
    aryubb[31]=17*aryubb[32] + 6*aryubb[31];
    aryubb[31]=MMZ*aryubb[31];
    aryubb[32]=aryubb[8]*aryubb[9];
    aryubb[33]=1./2.*aryubb[32];
    aryubb[36]= - 10 + aryubb[33];
    aryubb[36]=aryubb[18]*aryubb[28]*aryubb[36];
    aryubb[31]=3*aryubb[36] + aryubb[31];
    aryubb[31]=MMZ*aryubb[31];
    aryubb[35]=aryubb[34] - aryubb[35];
    aryubb[36]=MMZ*aryubb[35];
    aryubb[33]=1 + aryubb[33];
    aryubb[33]=aryubb[30]*aryubb[33];
    aryubb[33]=aryubb[36] + aryubb[33];
    aryubb[36]=3*MMZ;
    aryubb[33]=aryubb[33]*aryubb[36];
    aryubb[37]=1./2.*aryubb[22];
    aryubb[38]=aryubb[37] + 1;
    aryubb[38]=aryubb[38]*aryubb[22];
    aryubb[38]=aryubb[38] + aryubb[21];
    aryubb[38]=aryubb[38]*aryubb[3];
    aryubb[39]=aryubb[37]*aryubb[19];
    aryubb[38]=aryubb[38] - aryubb[39];
    aryubb[40]=aryubb[22]*aryubb[18];
    aryubb[41]=1./2.*aryubb[40] + aryubb[38];
    aryubb[41]=aryubb[41]*aryubb[32];
    aryubb[42]=3*aryubb[40];
    aryubb[33]=aryubb[33] + 3*aryubb[41] + aryubb[42] - aryubb[38];
    aryubb[33]=aryubb[14]*aryubb[33];
    aryubb[35]=aryubb[35]*aryubb[36];
    aryubb[29]=7*aryubb[29];
    aryubb[29]=aryubb[29]*aryubb[28];
    aryubb[29]=aryubb[35] + aryubb[29];
    aryubb[29]=aryubb[29]*MMZ;
    aryubb[38]=aryubb[24] + 2*aryubb[40];
    aryubb[29]=aryubb[29] + 4*aryubb[38];
    aryubb[38]= - aryubb[12] + aryubb[15];
    aryubb[38]=aryubb[29]*aryubb[38];
    aryubb[37]= - 2 + aryubb[37];
    aryubb[37]=aryubb[37]*aryubb[22];
    aryubb[28]=aryubb[28]*aryubb[22];
    aryubb[28]=5./6.*aryubb[21] - 4./3. - 3./2.*aryubb[28];
    aryubb[28]=aryubb[28]*aryubb[32];
    aryubb[28]=aryubb[28] + aryubb[37] - 2*aryubb[21];
    aryubb[28]=aryubb[13]*aryubb[3]*aryubb[28];
    aryubb[37]= - aryubb[16] + aryubb[17];
    aryubb[37]=aryubb[9]*aryubb[37];
    aryubb[39]= - aryubb[20]*aryubb[39];
    aryubb[21]=aryubb[28] + aryubb[33] + aryubb[38] + aryubb[31] + 
      aryubb[27] - 103./12.*aryubb[25] + aryubb[39] + 40./9.*aryubb[37] - 
      1./12.*aryubb[26] - 125./216.*aryubb[21] - 109./27. - 169./8.*
      aryubb[22];
    aryubb[21]=aryubb[11]*aryubb[21];
    aryubb[22]=3./2.*aryubb[32];
    aryubb[25]=11 - aryubb[22];
    aryubb[25]=aryubb[30]*aryubb[25];
    aryubb[25]=aryubb[35] + aryubb[25];
    aryubb[25]=MMZ*aryubb[25];
    aryubb[26]= - aryubb[24] - aryubb[40];
    aryubb[26]=aryubb[18]*aryubb[26];
    aryubb[27]= - aryubb[36]*aryubb[34];
    aryubb[26]=4*aryubb[26] + aryubb[27];
    aryubb[26]=aryubb[14]*aryubb[26];
    aryubb[27]=aryubb[22] + 13./2.;
    aryubb[27]=aryubb[24]*aryubb[27];
    aryubb[25]=aryubb[26] + aryubb[25] + 9*aryubb[40] + aryubb[27];
    aryubb[25]=aryubb[11]*aryubb[25];
    aryubb[26]=aryubb[1]*aryubb[3];
    aryubb[27]=2*aryubb[26];
    aryubb[27]=aryubb[23]*aryubb[27];
    aryubb[28]= - aryubb[7]*aryubb[27];
    aryubb[24]= - aryubb[24] - aryubb[42];
    aryubb[24]=aryubb[11]*aryubb[18]*aryubb[24];
    aryubb[24]=aryubb[28] + 1./2.*aryubb[24];
    aryubb[24]=aryubb[6]*aryubb[24];
    aryubb[23]=aryubb[23]*aryubb[26];
    aryubb[26]= - 4 - 3*aryubb[32];
    aryubb[26]=aryubb[23]*aryubb[26];
    aryubb[24]=aryubb[24] + aryubb[26] + aryubb[25];
    aryubb[24]=aryubb[6]*aryubb[24];
    aryubb[25]=aryubb[11]*aryubb[29];
    aryubb[25]= - aryubb[27] + aryubb[25];
    aryubb[25]=aryubb[10]*aryubb[25];
    aryubb[22]= - aryubb[22] + 39./4.;
    aryubb[22]=aryubb[23]*MMt*aryubb[22];
    aryubb[23]=aryubb[5]*aryubb[27];

    yubbret = aryubb[21] + aryubb[22] + aryubb[23] + aryubb[24] + 
      aryubb[25];
    return yubbret.real();
  }
} // namespace mr
