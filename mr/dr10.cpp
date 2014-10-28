#include <dr.hpp>
std::complex<long double>
dr::dr10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> ardr[18], drret;

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
   ardr[12]= - 1./2.*MMH + 3*ardr[8];
   ardr[12]= - 3./2.*ardr[9] + 1./2.*ardr[12] - 3*ardr[10];
   ardr[13]=pow(ardr[2],2);
   ardr[14]=ardr[13]*ardr[12];
   ardr[15]= - ardr[10] + ardr[9];
   ardr[16]=pow(ardr[4],2);
   ardr[15]=ardr[16]*ardr[15];
   ardr[12]=ardr[12] + 3./2.*ardr[15];
   ardr[12]=ardr[16]*ardr[12];
   ardr[12]=ardr[14] + ardr[12];
   ardr[12]=ardr[3]*ardr[12];
   ardr[14]=MMZ + 3./2.*ardr[9];
   ardr[14]=ardr[6]*ardr[14];
   ardr[14]= - 1./8. + ardr[14];
   ardr[14]=ardr[13]*ardr[14];
   ardr[15]= - ardr[8]*ardr[11];
   ardr[17]=ardr[10]*ardr[11];
   ardr[15]=ardr[17] - 1./2. + ardr[15];
   ardr[17]=1./2.*ardr[9] + MMZ + ardr[10];
   ardr[17]=ardr[6]*ardr[17];
   ardr[15]=1./4.*ardr[15] + ardr[17];
   ardr[15]=ardr[16]*ardr[15];
   ardr[17]= - ardr[6]*MMZ;
   ardr[12]=1./2.*ardr[12] + 3*ardr[15] + 2*ardr[17] + ardr[14];
   ardr[12]=ardr[7]*ardr[12];
   ardr[14]=1./2.*MMt + ardr[5];
   ardr[14]=ardr[1]*ardr[14];
   ardr[15]= - ardr[6]*ardr[1]*ardr[5]*MMt;
   ardr[14]=1./2.*ardr[14] + 2*ardr[15];
   ardr[13]=ardr[13]*ardr[14];
   ardr[14]=ardr[16]*ardr[14];
   ardr[13]=ardr[13] + ardr[14];
   ardr[13]=ardr[3]*ardr[13];

      drret = ardr[12] + 3*ardr[13];
      return drret;
}
