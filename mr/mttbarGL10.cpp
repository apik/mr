#include <tt.hpp>
long double tt::xgl10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armttbarGL[12], mttbarGLret;

    armttbarGL[1]=double(boson);
    armttbarGL[2]=Tsil::B(MMH,MMt,MMt,mu2);
    armttbarGL[3]=pow(SW,-1);
    armttbarGL[4]=pow(MMW,-1);
    armttbarGL[5]=Tsil::A(MMH,mu2);
    armttbarGL[6]=Tsil::A(MMt,mu2);
    armttbarGL[7]=pow(MMH,-1);
    armttbarGL[8]=std::real(Tsil::B(0,0,MMt,mu2));
   armttbarGL[9]=armttbarGL[7]*armttbarGL[6];
   armttbarGL[9]=1./2.*armttbarGL[2] + 1./8.*armttbarGL[8] - 3*
   armttbarGL[9];
   armttbarGL[9]=MMt*armttbarGL[9];
   armttbarGL[10]=armttbarGL[6] + armttbarGL[5];
   armttbarGL[11]=MMH*armttbarGL[2];
   armttbarGL[9]= - 1./8.*armttbarGL[11] + 1./4.*armttbarGL[10] + 
   armttbarGL[9];

      mttbarGLret = armttbarGL[9]*armttbarGL[4]*pow(armttbarGL[3],2)*
      armttbarGL[1];
      return mttbarGLret.real();
}
