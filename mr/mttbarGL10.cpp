#include <tt.hpp>
std::complex<long double>
tt::mgl10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armttbarGL[11], mttbarGLret;

    armttbarGL[1]=double(boson);
    armttbarGL[2]=Tsil::B(MMH,MMt,MMt,mu2);
    armttbarGL[3]=pow(SW,-1);
    armttbarGL[4]=pow(MMW,-1);
    armttbarGL[5]=Tsil::A(MMH,mu2);
    armttbarGL[6]=Tsil::A(MMt,mu2);
    armttbarGL[7]=pow(MMH,-1);
    armttbarGL[8]=std::real(Tsil::B(0,0,MMt,mu2));
   armttbarGL[9]=1./2.*armttbarGL[2];
   armttbarGL[10]=armttbarGL[6]*armttbarGL[7];
   armttbarGL[10]=1./8.*armttbarGL[8] + armttbarGL[9] - 3*
   armttbarGL[10];
   armttbarGL[10]=MMt*armttbarGL[10];
   armttbarGL[9]= - MMH*armttbarGL[9];
   armttbarGL[9]=armttbarGL[9] + armttbarGL[5] + armttbarGL[6];
   armttbarGL[9]=1./4.*armttbarGL[9] + armttbarGL[10];

      mttbarGLret = armttbarGL[9]*armttbarGL[4]*pow(armttbarGL[3],2)*
      armttbarGL[1];
      return mttbarGLret;
}
