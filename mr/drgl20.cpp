#include <dr.hpp>
long double dr::drgl20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> ardrgl[24], drglret;

    ardrgl[1]=double(boson);
    ardrgl[2]=pow(SW,-1);
    ardrgl[3]=pow(MMH,-1);
    ardrgl[4]=pow(MMW,-1);
    ardrgl[5]=Tsil::I2(MMH,MMt,MMt,mu2);
    ardrgl[6]=Tsil::I2(0,MMH,MMt,mu2);
    ardrgl[7]=Tsil::B(MMH,MMH,MMH,mu2);
    ardrgl[8]=Tsil::A(MMH,mu2);
    ardrgl[9]=Tsil::A(MMt,mu2);
    ardrgl[10]=Tsil::B(MMH,MMt,MMt,mu2);
    ardrgl[11]=Tsil::B(MMt,MMt,MMH,mu2);
    ardrgl[12]=std::real(Tsil::B(0,0,MMH,mu2));
    ardrgl[13]=pow(MMt,-1);
    ardrgl[14]=std::real(Tsil::B(0,0,MMt,mu2));
   ardrgl[15]=5*ardrgl[8];
   ardrgl[16]=3./4.*ardrgl[8];
   ardrgl[17]=ardrgl[9] - ardrgl[16];
   ardrgl[17]=ardrgl[17]*ardrgl[15];
   ardrgl[18]=pow(ardrgl[9],2);
   ardrgl[17]=11*ardrgl[18] + ardrgl[17];
   ardrgl[17]=ardrgl[3]*ardrgl[17];
   ardrgl[19]=ardrgl[14]*ardrgl[9];
   ardrgl[17]=ardrgl[19] + ardrgl[17];
   ardrgl[20]=3./2.*ardrgl[12];
   ardrgl[21]=9./2.*ardrgl[7] + ardrgl[20] + 9./8.;
   ardrgl[21]=ardrgl[9]*ardrgl[21];
   ardrgl[22]=ardrgl[10]*ardrgl[9];
   ardrgl[16]=ardrgl[11]*ardrgl[16];
   ardrgl[16]= - 9./8.*ardrgl[6] + 3*ardrgl[22] - 3./4.*ardrgl[5] + 
   ardrgl[16] + ardrgl[8] + ardrgl[21] + 1./4.*ardrgl[17];
   ardrgl[17]=pow(Pi,2);
   ardrgl[21]= - 3./32.*ardrgl[17] - 3./16.*ardrgl[10] - 1 + 5./16.*
   ardrgl[11];
   ardrgl[21]=MMH*ardrgl[21];
   ardrgl[16]=1./2.*ardrgl[16] + ardrgl[21];
   ardrgl[15]= - 19./2.*ardrgl[9] - ardrgl[15];
   ardrgl[21]=ardrgl[9] - 1./2.*ardrgl[8];
   ardrgl[21]=ardrgl[11]*ardrgl[21];
   ardrgl[23]=ardrgl[3]*pow(ardrgl[8],2);
   ardrgl[15]=3./8.*ardrgl[6] + 21./16.*ardrgl[23] + 9./4.*ardrgl[5] - 
   ardrgl[19] + 1./4.*ardrgl[15] + 3*ardrgl[21];
   ardrgl[15]=ardrgl[3]*ardrgl[15];
   ardrgl[19]=3./4.*ardrgl[14] + 59./16. - 5*ardrgl[11];
   ardrgl[21]=ardrgl[3]*ardrgl[9];
   ardrgl[23]=5./4. - 4*ardrgl[21];
   ardrgl[23]=ardrgl[10]*ardrgl[23];
   ardrgl[15]=ardrgl[23] + 1./4.*ardrgl[19] + ardrgl[15];
   ardrgl[19]=77./8. - ardrgl[14];
   ardrgl[21]=ardrgl[11]*ardrgl[21];
   ardrgl[19]= - 2*ardrgl[10] + 1./2.*ardrgl[19] - 12*ardrgl[21];
   ardrgl[19]=3*ardrgl[19] + 29./32.*ardrgl[17];
   ardrgl[19]=MMt*ardrgl[3]*ardrgl[19];
   ardrgl[15]=ardrgl[19] + 3*ardrgl[15] + 35./64.*ardrgl[17];
   ardrgl[15]=MMt*ardrgl[15];
   ardrgl[15]=3*ardrgl[16] + ardrgl[15];
   ardrgl[16]=pow(ardrgl[4],2);
   ardrgl[15]=MMt*ardrgl[16]*ardrgl[15];
   ardrgl[19]=ardrgl[13]*ardrgl[18];
   ardrgl[19]= - 3*ardrgl[19] - 11*ardrgl[8];
   ardrgl[20]=ardrgl[8]*ardrgl[20];
   ardrgl[19]=ardrgl[20] + 1./2.*ardrgl[19] - 3*ardrgl[5];
   ardrgl[20]=ardrgl[7]*ardrgl[8];
   ardrgl[19]=9./4.*ardrgl[20] + 3./2.*ardrgl[6] + 1./2.*ardrgl[19] - 
   ardrgl[22];
   ardrgl[20]=15./2.*ardrgl[12] + 131./4. - 243*S2;
   ardrgl[17]=45./8.*ardrgl[7] + 1./4.*ardrgl[20] + 1./3.*ardrgl[17];
   ardrgl[17]=MMH*ardrgl[17];
   ardrgl[17]=3*ardrgl[19] + ardrgl[17];
   ardrgl[17]=MMH*ardrgl[17];
   ardrgl[19]= - 5*ardrgl[9] + 9./2.*ardrgl[8];
   ardrgl[19]=ardrgl[8]*ardrgl[19];
   ardrgl[18]=ardrgl[18] + ardrgl[19];
   ardrgl[17]=3./2.*ardrgl[18] + ardrgl[17];
   ardrgl[16]=ardrgl[17]*ardrgl[16];
   ardrgl[15]=1./8.*ardrgl[16] + ardrgl[15];

      drglret = ardrgl[15]*pow(ardrgl[2],4)*ardrgl[1];
      return drglret.real();
}
