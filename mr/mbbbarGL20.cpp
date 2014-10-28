#include <bb.hpp>
std::complex<long double>
bb::mgl20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armbbbarGL[27], mbbbarGLret;

    armbbbarGL[1]=double(boson);
    armbbbarGL[2]=pow(SW,-1);
    armbbbarGL[3]=pow(MMH,-1);
    armbbbarGL[4]=pow(MMW,-1);
    armbbbarGL[5]=Tsil::I2(MMH,MMH,MMH,mu2);
    armbbbarGL[6]=Tsil::I2(MMH,MMt,MMt,mu2);
    armbbbarGL[7]=pow(MMt,-1);
    armbbbarGL[8]=Tsil::I2(0,MMH,MMt,mu2);
    armbbbarGL[9]=Tsil::I2(0,0,MMH,mu2);
    armbbbarGL[10]=Tsil::I2(0,0,MMt,mu2);
    armbbbarGL[11]=Tsil::B(MMH,MMH,MMH,mu2);
    armbbbarGL[12]=Tsil::A(MMH,mu2);
    armbbbarGL[13]=Tsil::A(MMt,mu2);
    armbbbarGL[14]=Tsil::B(MMH,MMt,MMt,mu2);
    armbbbarGL[15]=Tsil::B(MMt,MMt,MMH,mu2);
    armbbbarGL[16]=std::real(Tsil::B(0,0,MMH,mu2));
    armbbbarGL[17]=std::real(Tsil::B(0,0,MMt,mu2));
    armbbbarGL[18]=Tsil::Aeps(MMH,mu2);
    armbbbarGL[19]=Tsil::Aeps(MMt,mu2);
   armbbbarGL[20]=9*armbbbarGL[16];
   armbbbarGL[21]=27*armbbbarGL[11];
   armbbbarGL[22]=armbbbarGL[3]*armbbbarGL[12];
   armbbbarGL[22]=11*armbbbarGL[22] + 15*armbbbarGL[14] + 3./4.*
   armbbbarGL[17] + armbbbarGL[21] + 61./16. + armbbbarGL[20];
   armbbbarGL[23]=armbbbarGL[13]*armbbbarGL[3];
   armbbbarGL[22]=1./2.*armbbbarGL[22] + 9*armbbbarGL[23];
   armbbbarGL[22]=armbbbarGL[13]*armbbbarGL[22];
   armbbbarGL[23]= - 31*armbbbarGL[10] + 75*armbbbarGL[18];
   armbbbarGL[24]=9*armbbbarGL[15];
   armbbbarGL[25]=49./8. + armbbbarGL[24];
   armbbbarGL[25]=armbbbarGL[12]*armbbbarGL[25];
   armbbbarGL[24]= - 7./4.*armbbbarGL[14] - 85./8. + armbbbarGL[24];
   armbbbarGL[24]=MMH*armbbbarGL[24];
   armbbbarGL[22]=1./4.*armbbbarGL[24] + armbbbarGL[22] + 1./4.*
   armbbbarGL[25] - 35./8.*armbbbarGL[6] - 5./16.*armbbbarGL[8] + 1./16.
   *armbbbarGL[23] + 11*armbbbarGL[19];
   armbbbarGL[23]=1./2.*armbbbarGL[10] - armbbbarGL[18];
   armbbbarGL[24]= - 5./2. - 3*armbbbarGL[15];
   armbbbarGL[24]=armbbbarGL[12]*armbbbarGL[24];
   armbbbarGL[23]=3*armbbbarGL[24] + 27./2.*armbbbarGL[6] - 1./2.*
   armbbbarGL[8] + 13*armbbbarGL[23] - 33*armbbbarGL[19];
   armbbbarGL[23]=armbbbarGL[3]*armbbbarGL[23];
   armbbbarGL[24]=463./8. + 7*armbbbarGL[17];
   armbbbarGL[23]=armbbbarGL[23] + 19./4.*armbbbarGL[14] + 1./16.*
   armbbbarGL[24] - 9*armbbbarGL[15];
   armbbbarGL[24]=3*armbbbarGL[15] - 13./8. - armbbbarGL[17];
   armbbbarGL[24]=1./2.*armbbbarGL[24] - 2*armbbbarGL[14];
   armbbbarGL[24]=armbbbarGL[3]*armbbbarGL[24];
   armbbbarGL[25]=pow(armbbbarGL[3],2);
   armbbbarGL[26]= - armbbbarGL[13]*armbbbarGL[25];
   armbbbarGL[24]=armbbbarGL[24] + 3./2.*armbbbarGL[26];
   armbbbarGL[24]=armbbbarGL[13]*armbbbarGL[24];
   armbbbarGL[26]=3 - 1./2.*armbbbarGL[17];
   armbbbarGL[26]=1./2.*armbbbarGL[26] - armbbbarGL[14];
   armbbbarGL[26]=armbbbarGL[3]*armbbbarGL[26];
   armbbbarGL[25]= - armbbbarGL[13]*armbbbarGL[25]*armbbbarGL[15];
   armbbbarGL[25]=armbbbarGL[26] + 6*armbbbarGL[25];
   armbbbarGL[25]=MMt*armbbbarGL[25];
   armbbbarGL[23]=3*armbbbarGL[25] + 1./4.*armbbbarGL[23] + 3*
   armbbbarGL[24];
   armbbbarGL[23]=MMt*armbbbarGL[23];
   armbbbarGL[22]=1./4.*armbbbarGL[22] + armbbbarGL[23];
   armbbbarGL[22]=MMt*armbbbarGL[22];
   armbbbarGL[20]=armbbbarGL[21] + 31 + armbbbarGL[20];
   armbbbarGL[20]=armbbbarGL[12]*armbbbarGL[20];
   armbbbarGL[21]= - armbbbarGL[12]*armbbbarGL[7];
   armbbbarGL[23]= - armbbbarGL[13]*armbbbarGL[7];
   armbbbarGL[21]=armbbbarGL[23] + armbbbarGL[21] + 1./2. - 3*
   armbbbarGL[14];
   armbbbarGL[21]=armbbbarGL[13]*armbbbarGL[21];
   armbbbarGL[23]= - 15*armbbbarGL[5] - 17*armbbbarGL[9];
   armbbbarGL[20]=1./4.*armbbbarGL[21] + 1./8.*armbbbarGL[20] - 13./8.*
   armbbbarGL[6] + 19./8.*armbbbarGL[8] + 7./8.*armbbbarGL[19] + 1./8.*
   armbbbarGL[23] + 7*armbbbarGL[18];
   armbbbarGL[21]=1./2.*armbbbarGL[6] - 3./2.*armbbbarGL[8] + 
   armbbbarGL[9] + 1./2.*armbbbarGL[19];
   armbbbarGL[21]=armbbbarGL[7]*armbbbarGL[21];
   armbbbarGL[21]=1./8.*armbbbarGL[21] + 27./32.*armbbbarGL[11] - 1 + 9.
   /32.*armbbbarGL[16];
   armbbbarGL[21]=MMH*armbbbarGL[21];
   armbbbarGL[20]=1./4.*armbbbarGL[20] + armbbbarGL[21];
   armbbbarGL[20]=MMH*armbbbarGL[20];
   armbbbarGL[21]= - armbbbarGL[12] + 7./8.*armbbbarGL[13];
   armbbbarGL[21]=armbbbarGL[13]*armbbbarGL[21];
   armbbbarGL[23]=pow(armbbbarGL[12],2);
   armbbbarGL[21]= - 1./8.*armbbbarGL[23] + armbbbarGL[21];
   armbbbarGL[20]=9./8.*armbbbarGL[21] + armbbbarGL[20];
   armbbbarGL[20]=1./2.*armbbbarGL[20] + armbbbarGL[22];

      mbbbarGLret = armbbbarGL[20]*pow(armbbbarGL[4],2)*pow(
      armbbbarGL[2],4)*armbbbarGL[1];
      return mbbbarGLret;
}
