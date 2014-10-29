#include <dr.hpp>
std::complex<long double>
dr::dr10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> ardr[21], drret;

    ardr[1]=double(nH);
    ardr[2]=pow(CW,-1);
    ardr[3]=pow(MMZ,-1);
    ardr[4]=pow(SW,-1);
    ardr[5]=Tsil::A(MMt,mu2);
    ardr[6]=pow(MMH,-1);
    ardr[7]=double(boson);
    ardr[8]=Tsil::A(MMH,mu2);
    ardr[9]=Tsil::A(MMZ,mu2);
    ardr[10]=Tsil::A(MMW,mu2);
    ardr[11]=1/( - MMW + MMH);
   ardr[12]=1./2.*ardr[9] + ardr[10];
   ardr[12]=ardr[6]*ardr[12];
   ardr[13]=MMZ*ardr[6];
   ardr[14]=ardr[13] - 1./8.;
   ardr[12]=ardr[14] + ardr[12];
   ardr[15]=pow(ardr[4],2);
   ardr[16]= - ardr[10] + ardr[9];
   ardr[16]=ardr[16]*ardr[15];
   ardr[17]=3./2.*ardr[9];
   ardr[18]=ardr[17] + 3*ardr[10] - 3./2.*ardr[8] + 1./4.*MMH;
   ardr[16]=3./2.*ardr[16] - ardr[18];
   ardr[19]=1./2.*ardr[3];
   ardr[16]=ardr[16]*ardr[19];
   ardr[20]=ardr[10] - ardr[8];
   ardr[20]=ardr[11]*ardr[20];
   ardr[12]=3./4.*ardr[20] + ardr[16] + 3*ardr[12];
   ardr[12]=ardr[15]*ardr[12];
   ardr[16]= - ardr[18]*ardr[19];
   ardr[17]=ardr[6]*ardr[17];
   ardr[14]=ardr[16] + ardr[17] + ardr[14];
   ardr[16]=pow(ardr[2],2);
   ardr[14]=ardr[14]*ardr[16];
   ardr[12]=ardr[14] - 2*ardr[13] + ardr[12];
   ardr[12]=ardr[7]*ardr[12];
   ardr[13]= - ardr[15] - ardr[16];
   ardr[14]=2*ardr[6];
   ardr[14]=ardr[14]*ardr[5];
   ardr[14]=ardr[14] - 1./4.;
   ardr[14]=ardr[14]*MMt;
   ardr[14]=ardr[14] - 1./2.*ardr[5];
   ardr[13]=ardr[1]*ardr[3]*ardr[14]*ardr[13];

      drret = ardr[12] + 3*ardr[13];
      return drret;
}
