#include <bb.hpp>
std::complex<long double>
bb::mgl11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armbbbarGL[17], mbbbarGLret;

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
   armbbbarGL[13]=MMt*armbbbarGL[3];
   armbbbarGL[14]= - MMt*armbbbarGL[3];
   armbbbarGL[14]=1./2. + 4*armbbbarGL[14];
   armbbbarGL[14]=armbbbarGL[7]*armbbbarGL[8]*armbbbarGL[14];
   armbbbarGL[15]= - 1./2.*armbbbarGL[10] - 24*armbbbarGL[3];
   armbbbarGL[15]=armbbbarGL[9]*armbbbarGL[15];
   armbbbarGL[14]=armbbbarGL[15] + 3*armbbbarGL[14] + 13./2. + 4*
   armbbbarGL[13];
   armbbbarGL[14]=armbbbarGL[9]*armbbbarGL[14];
   armbbbarGL[15]=3*armbbbarGL[6] + 1./2.*MMt;
   armbbbarGL[15]=armbbbarGL[8]*armbbbarGL[15];
   armbbbarGL[16]= - armbbbarGL[7]*armbbbarGL[8];
   armbbbarGL[15]=armbbbarGL[15] + 5*armbbbarGL[16];
   armbbbarGL[15]=armbbbarGL[7]*armbbbarGL[15];
   armbbbarGL[13]= - 103./12. + 32*armbbbarGL[13];
   armbbbarGL[13]=MMt*armbbbarGL[13];
   armbbbarGL[13]=armbbbarGL[14] + 1./2.*armbbbarGL[15] + 
   armbbbarGL[13] - 1./2.*armbbbarGL[6] - 4*armbbbarGL[5] + 5*
   armbbbarGL[12] + 4*armbbbarGL[11];

      mbbbarGLret = armbbbarGL[13]*armbbbarGL[4]*pow(armbbbarGL[2],2)*
      armbbbarGL[1];
      return mbbbarGLret;
}
