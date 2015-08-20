#include <dr.hpp>
namespace mr
{
  long double dr<MS>::dr10(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<long double> ardrbar[22], drbarret;

    ardrbar[1]=double(nH);
    ardrbar[2]=pow(mmZ,-1);
    ardrbar[3]=pow(s,-1);
    ardrbar[4]=pow(c,-1);
    ardrbar[5]=Tsil::A(mmt,mu2);
    ardrbar[6]=pow(mmH,-1);
    ardrbar[7]=Tsil::A(mmb,mu2);
    ardrbar[8]=double(boson);
    ardrbar[9]=Tsil::A(mmW,mu2);
    ardrbar[10]=Tsil::A(mmZ,mu2);
    ardrbar[11]=Tsil::A(mmH,mu2);
    ardrbar[12]=1/( - mmb + mmt);
    ardrbar[13]=1/(mmH - mmW);
    ardrbar[14]=mmb*ardrbar[7];
    ardrbar[15]=mmt*ardrbar[5];
    ardrbar[14]=ardrbar[14] + ardrbar[15];
    ardrbar[14]= - ardrbar[6]*ardrbar[14];
    ardrbar[15]=ardrbar[5] + 1./2.*mmt;
    ardrbar[16]=ardrbar[7] + 1./2.*mmb - ardrbar[15];
    ardrbar[16]= - ardrbar[12]*mmb*ardrbar[16];
    ardrbar[14]=3./2.*ardrbar[16] + 6*ardrbar[14];
    ardrbar[16]=pow(ardrbar[4],2);
    ardrbar[17]=pow(ardrbar[3],2);
    ardrbar[18]=ardrbar[16] + ardrbar[17];
    ardrbar[18]=ardrbar[18]*ardrbar[1];
    ardrbar[14]=ardrbar[18]*ardrbar[14];
    ardrbar[19]= - ardrbar[9] + 1./2.*ardrbar[11];
    ardrbar[20]=3./2.*ardrbar[10];
    ardrbar[19]= - ardrbar[20] - 1./4.*mmH + 3*ardrbar[19];
    ardrbar[21]= - ardrbar[9] + ardrbar[10];
    ardrbar[21]=ardrbar[21]*ardrbar[17];
    ardrbar[21]=3./2.*ardrbar[21] + ardrbar[19];
    ardrbar[21]=ardrbar[21]*ardrbar[17];
    ardrbar[19]=ardrbar[19]*ardrbar[16];
    ardrbar[19]=ardrbar[19] + ardrbar[21];
    ardrbar[19]=ardrbar[8]*ardrbar[19];
    ardrbar[15]=ardrbar[15]*ardrbar[18];
    ardrbar[15]=3*ardrbar[15] + ardrbar[19];
    ardrbar[14]=1./2.*ardrbar[15] + ardrbar[14];
    ardrbar[14]=ardrbar[2]*ardrbar[14];
    ardrbar[15]=ardrbar[20] + mmZ;
    ardrbar[15]=ardrbar[15]*ardrbar[16];
    ardrbar[17]=3*ardrbar[17];
    ardrbar[18]=mmZ + ardrbar[9] + 1./2.*ardrbar[10];
    ardrbar[18]=ardrbar[18]*ardrbar[17];
    ardrbar[15]=ardrbar[18] - 2*mmZ + ardrbar[15];
    ardrbar[15]=ardrbar[6]*ardrbar[15];
    ardrbar[18]= - ardrbar[11] + ardrbar[9];
    ardrbar[18]=ardrbar[13]*ardrbar[18];
    ardrbar[18]= - 1./2. + ardrbar[18];
    ardrbar[17]=ardrbar[18]*ardrbar[17];
    ardrbar[16]= - 1./2.*ardrbar[16] + ardrbar[17];
    ardrbar[15]=1./4.*ardrbar[16] + ardrbar[15];
    ardrbar[15]=ardrbar[8]*ardrbar[15];

    drbarret = ardrbar[14] + ardrbar[15];
    return drbarret.real();
  }
} // namespace mr
