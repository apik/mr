#include <bb.hpp>
std::complex<long double>
bb::mgl11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armbbbarGL[18], mbbbarGLret;

    armbbbarGL[1]=double(boson);
    armbbbarGL[2]=pow(SW,-1);
    armbbbarGL[3]=pow(MMH,-1);
    armbbbarGL[4]=pow(MMW,-1);
    armbbbarGL[5]=Tsil::I2(0,0,MMt,mu2);
    armbbbarGL[6]=Tsil::A(MMH,mu2);
    armbbbarGL[7]=Tsil::A(MMb,mu2);
    armbbbarGL[8]=pow(MMb,-1);
    armbbbarGL[9]=Tsil::A(MMt,mu2);
    armbbbarGL[10]=pow(MMt,-1);
    armbbbarGL[11]=Tsil::Aeps(MMt,mu2);
    armbbbarGL[12]=Tsil::Aeps(MMb,mu2);
   armbbbarGL[13]=pow(MMt,2);
   armbbbarGL[14]=armbbbarGL[8]*armbbbarGL[7];
   armbbbarGL[15]=3*armbbbarGL[14];
   armbbbarGL[16]=armbbbarGL[15] - 1;
   armbbbarGL[17]= - MMt*armbbbarGL[16];
   armbbbarGL[17]=armbbbarGL[17] - 6*armbbbarGL[9];
   armbbbarGL[17]=armbbbarGL[9]*armbbbarGL[17];
   armbbbarGL[13]=8*armbbbarGL[13] + armbbbarGL[17];
   armbbbarGL[13]=armbbbarGL[3]*armbbbarGL[13];
   armbbbarGL[13]= - armbbbarGL[5] + armbbbarGL[13] + armbbbarGL[11];
   armbbbarGL[17]= - armbbbarGL[9]*armbbbarGL[10];
   armbbbarGL[15]=armbbbarGL[17] + 13 + armbbbarGL[15];
   armbbbarGL[15]=armbbbarGL[9]*armbbbarGL[15];
   armbbbarGL[16]=armbbbarGL[6]*armbbbarGL[16];
   armbbbarGL[15]=armbbbarGL[15] + armbbbarGL[16];
   armbbbarGL[16]=armbbbarGL[8]*pow(armbbbarGL[7],2);
   armbbbarGL[14]= - 103./3. + armbbbarGL[14];
   armbbbarGL[14]=MMt*armbbbarGL[14];
   armbbbarGL[13]=1./4.*armbbbarGL[14] - 5./2.*armbbbarGL[16] + 5*
   armbbbarGL[12] + 1./2.*armbbbarGL[15] + 4*armbbbarGL[13];

      mbbbarGLret = armbbbarGL[13]*armbbbarGL[4]*pow(armbbbarGL[2],2)*
      armbbbarGL[1];
      return mbbbarGLret;
}
