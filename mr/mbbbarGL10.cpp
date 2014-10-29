#include <bb.hpp>
std::complex<long double>
bb::mgl10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armbbbarGL[9], mbbbarGLret;

    armbbbarGL[1]=double(boson);
    armbbbarGL[2]=pow(SW,-1);
    armbbbarGL[3]=pow(MMW,-1);
    armbbbarGL[4]=Tsil::A(MMH,mu2);
    armbbbarGL[5]=Tsil::A(MMt,mu2);
    armbbbarGL[6]=pow(MMH,-1);
   armbbbarGL[7]=armbbbarGL[4] + armbbbarGL[5];
   armbbbarGL[8]=armbbbarGL[6]*armbbbarGL[5];
   armbbbarGL[8]=1./16. - 3*armbbbarGL[8];
   armbbbarGL[8]=MMt*armbbbarGL[8];
   armbbbarGL[7]=3./8.*armbbbarGL[7] + armbbbarGL[8];

      mbbbarGLret = armbbbarGL[7]*armbbbarGL[3]*pow(armbbbarGL[2],2)*
      armbbbarGL[1];
      return mbbbarGLret;
}
