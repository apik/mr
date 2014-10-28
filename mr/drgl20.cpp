#include <dr.hpp>
std::complex<long double>
dr::drgl20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> ardrgl[21], drglret;

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
   ardrgl[15]= - 3./2.*ardrgl[6] - ardrgl[5];
   ardrgl[16]=1 + 3./4.*ardrgl[11];
   ardrgl[16]=ardrgl[8]*ardrgl[16];
   ardrgl[15]=3./4.*ardrgl[15] + ardrgl[16];
   ardrgl[16]=5*ardrgl[8] + 11*ardrgl[9];
   ardrgl[16]=ardrgl[9]*ardrgl[16];
   ardrgl[17]=pow(ardrgl[8],2);
   ardrgl[16]= - 15./4.*ardrgl[17] + ardrgl[16];
   ardrgl[16]=ardrgl[3]*ardrgl[16];
   ardrgl[18]=pow(Pi,2);
   ardrgl[19]= - 3./16.*ardrgl[10] + 5./16.*ardrgl[11] - 1 - 3./32.*
   ardrgl[18];
   ardrgl[19]=MMH*ardrgl[19];
   ardrgl[20]=3*ardrgl[7] + 3./4. + ardrgl[12];
   ardrgl[20]=3*ardrgl[20] + 1./2.*ardrgl[14];
   ardrgl[20]=1./2.*ardrgl[20] + 3*ardrgl[10];
   ardrgl[20]=ardrgl[9]*ardrgl[20];
   ardrgl[15]=1./8.*ardrgl[16] + 1./2.*ardrgl[20] + 1./2.*ardrgl[15] + 
   ardrgl[19];
   ardrgl[16]=1./2.*ardrgl[6] + 3*ardrgl[5];
   ardrgl[19]= - 5./2. - 3*ardrgl[11];
   ardrgl[19]=ardrgl[8]*ardrgl[19];
   ardrgl[16]=3./2.*ardrgl[16] + ardrgl[19];
   ardrgl[19]= - 4*ardrgl[10] + 3*ardrgl[11] - 19./8. - ardrgl[14];
   ardrgl[19]=ardrgl[9]*ardrgl[19];
   ardrgl[20]=ardrgl[3]*ardrgl[17];
   ardrgl[16]=21./16.*ardrgl[20] + 1./2.*ardrgl[16] + ardrgl[19];
   ardrgl[16]=ardrgl[3]*ardrgl[16];
   ardrgl[19]=231 + 29./2.*ardrgl[18];
   ardrgl[19]=1./8.*ardrgl[19] - 3*ardrgl[14];
   ardrgl[20]= - ardrgl[3]*ardrgl[9]*ardrgl[11];
   ardrgl[19]=36*ardrgl[20] + 1./2.*ardrgl[19] - 6*ardrgl[10];
   ardrgl[19]=MMt*ardrgl[3]*ardrgl[19];
   ardrgl[20]=177 + 35*ardrgl[18];
   ardrgl[20]=1./4.*ardrgl[20] + 9*ardrgl[14];
   ardrgl[20]=15*ardrgl[10] + 1./4.*ardrgl[20] - 15*ardrgl[11];
   ardrgl[16]=ardrgl[19] + 1./4.*ardrgl[20] + 3*ardrgl[16];
   ardrgl[16]=MMt*ardrgl[16];
   ardrgl[15]=3*ardrgl[15] + ardrgl[16];
   ardrgl[15]=MMt*ardrgl[15];
   ardrgl[16]=ardrgl[6] - ardrgl[5];
   ardrgl[19]=9*ardrgl[7] - 11 + 3*ardrgl[12];
   ardrgl[19]=ardrgl[8]*ardrgl[19];
   ardrgl[16]=3*ardrgl[16] + 1./2.*ardrgl[19];
   ardrgl[19]=45./2.*ardrgl[7] + 15./2.*ardrgl[12] + 131./4. - 243*S2;
   ardrgl[18]=1./4.*ardrgl[19] + 1./3.*ardrgl[18];
   ardrgl[18]=MMH*ardrgl[18];
   ardrgl[16]=3./2.*ardrgl[16] + ardrgl[18];
   ardrgl[16]=MMH*ardrgl[16];
   ardrgl[18]= - MMH*ardrgl[10];
   ardrgl[19]= - MMH*ardrgl[13];
   ardrgl[19]=1 + 3./2.*ardrgl[19];
   ardrgl[19]=ardrgl[9]*ardrgl[19];
   ardrgl[18]=1./2.*ardrgl[19] - 5./2.*ardrgl[8] + ardrgl[18];
   ardrgl[18]=ardrgl[9]*ardrgl[18];
   ardrgl[16]=3*ardrgl[18] + 27./4.*ardrgl[17] + ardrgl[16];
   ardrgl[15]=1./8.*ardrgl[16] + ardrgl[15];

      drglret = ardrgl[15]*pow(ardrgl[4],2)*pow(ardrgl[2],4)*ardrgl[1];
      return drglret;
}
