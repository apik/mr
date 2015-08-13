#include <dr.hpp>
namespace mr
{
  long double dr<OS>::drgl11(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<long double> ardrgl[9], drglret;

    ardrgl[1]=double(boson);
    ardrgl[2]=pow(SW,-1);
    ardrgl[3]=pow(MMH,-1);
    ardrgl[4]=pow(MMW,-1);
    ardrgl[5]=Tsil::A(MMt,mu2);
    ardrgl[6]=pow(MMt,-1);
    ardrgl[7]=MMt*ardrgl[3];
    ardrgl[8]=pow(Pi,2);
    ardrgl[7]=2./3.*ardrgl[8] - 29./2. + 64*ardrgl[7];
    ardrgl[7]=MMt*ardrgl[7];
    ardrgl[8]=ardrgl[6] - 8*ardrgl[3];
    ardrgl[8]=ardrgl[5]*ardrgl[8];
    ardrgl[8]=1 + ardrgl[8];
    ardrgl[8]=ardrgl[5]*ardrgl[8];
    ardrgl[7]=6*ardrgl[8] + ardrgl[7];

    drglret = ardrgl[7]*ardrgl[4]*pow(ardrgl[2],2)*ardrgl[1];
    return drglret.real();
  }
} // namespace mr
