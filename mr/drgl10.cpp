#include <dr.hpp>
std::complex<long double>
dr::drgl10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> ardrgl[9], drglret;

    ardrgl[1]=double(boson);
    ardrgl[2]=pow(SW,-1);
    ardrgl[3]=pow(MMW,-1);
    ardrgl[4]=Tsil::A(MMH,mu2);
    ardrgl[5]=Tsil::A(MMt,mu2);
    ardrgl[6]=pow(MMH,-1);
   ardrgl[7]=3*MMt - 1./2.*MMH + 3*ardrgl[4];
   ardrgl[8]= - MMt*ardrgl[6];
   ardrgl[8]=1./2. + 2*ardrgl[8];
   ardrgl[8]=ardrgl[5]*ardrgl[8];
   ardrgl[7]=1./4.*ardrgl[7] + 3*ardrgl[8];

      drglret = ardrgl[7]*ardrgl[3]*pow(ardrgl[2],2)*ardrgl[1];
      return drglret;
}
