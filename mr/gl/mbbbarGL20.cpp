#include <bb.hpp>
namespace mr
{
  long double bb<OS>::xgl20(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<long double> armbbbarGL[29], mbbbarGLret;

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
    armbbbarGL[20]=3*armbbbarGL[3];
    armbbbarGL[21]=3*armbbbarGL[15];
    armbbbarGL[22]= - 5./2. - armbbbarGL[21];
    armbbbarGL[22]=armbbbarGL[12]*armbbbarGL[22]*armbbbarGL[20];
    armbbbarGL[23]=13./2.*armbbbarGL[10] + 27./2.*armbbbarGL[6] - 13*
      armbbbarGL[18] - 33*armbbbarGL[19];
    armbbbarGL[23]=armbbbarGL[3]*armbbbarGL[23];
    armbbbarGL[24]=9*armbbbarGL[15];
    armbbbarGL[25]=1./2.*armbbbarGL[3];
    armbbbarGL[26]= - armbbbarGL[8]*armbbbarGL[25];
    armbbbarGL[22]=7./16.*armbbbarGL[17] + armbbbarGL[22] + 
      armbbbarGL[26] + 463./128. - armbbbarGL[24] + armbbbarGL[23];
    armbbbarGL[23]= - armbbbarGL[17]*armbbbarGL[25];
    armbbbarGL[20]=armbbbarGL[20] + armbbbarGL[23];
    armbbbarGL[23]=3*MMt;
    armbbbarGL[20]=armbbbarGL[20]*armbbbarGL[23];
    armbbbarGL[20]=1./2.*armbbbarGL[22] + armbbbarGL[20];
    armbbbarGL[20]=MMt*armbbbarGL[20];
    armbbbarGL[22]=1./4.*MMH;
    armbbbarGL[26]= - 85./8. + armbbbarGL[24];
    armbbbarGL[26]=armbbbarGL[26]*armbbbarGL[22];
    armbbbarGL[24]=49./8. + armbbbarGL[24];
    armbbbarGL[24]=armbbbarGL[12]*armbbbarGL[24];
    armbbbarGL[24]= - 31./16.*armbbbarGL[10] + 1./4.*armbbbarGL[24] - 5./
      16.*armbbbarGL[8] + armbbbarGL[26] - 35./8.*armbbbarGL[6] + 75./16.*
      armbbbarGL[18] + 11*armbbbarGL[19];
    armbbbarGL[20]=1./2.*armbbbarGL[24] + armbbbarGL[20];
    armbbbarGL[20]=MMt*armbbbarGL[20];
    armbbbarGL[24]=armbbbarGL[18] + 1./8.*armbbbarGL[19];
    armbbbarGL[24]=7*armbbbarGL[24] - 13./8.*armbbbarGL[6];
    armbbbarGL[26]=27./32.*armbbbarGL[11] - 1 + 9./32.*armbbbarGL[16];
    armbbbarGL[26]=MMH*armbbbarGL[26];
    armbbbarGL[24]= - 15./32.*armbbbarGL[5] + 19./32.*armbbbarGL[8] - 17.
      /32.*armbbbarGL[9] + 1./4.*armbbbarGL[24] + armbbbarGL[26];
    armbbbarGL[24]=MMH*armbbbarGL[24];
    armbbbarGL[26]=27*armbbbarGL[11] + 9*armbbbarGL[16];
    armbbbarGL[27]=31 + armbbbarGL[26];
    armbbbarGL[27]=MMH*armbbbarGL[27];
    armbbbarGL[27]=armbbbarGL[27] - 9./2.*armbbbarGL[12];
    armbbbarGL[27]=armbbbarGL[12]*armbbbarGL[27];
    armbbbarGL[28]=armbbbarGL[19] + armbbbarGL[6];
    armbbbarGL[28]= - 3./2.*armbbbarGL[8] + 1./2.*armbbbarGL[28] + 
      armbbbarGL[9];
    armbbbarGL[28]=armbbbarGL[7]*armbbbarGL[28]*pow(MMH,2);
    armbbbarGL[20]=armbbbarGL[20] + 1./8.*armbbbarGL[28] + 1./32.*
      armbbbarGL[27] + armbbbarGL[24];
    armbbbarGL[20]=armbbbarGL[1]*armbbbarGL[20];
    armbbbarGL[24]=armbbbarGL[12]*armbbbarGL[3];
    armbbbarGL[24]=3./4.*armbbbarGL[17] + 11*armbbbarGL[24] + 61./16. + 
      armbbbarGL[26];
    armbbbarGL[21]= - armbbbarGL[17] - 13./8. + armbbbarGL[21];
    armbbbarGL[21]=armbbbarGL[3]*armbbbarGL[21];
    armbbbarGL[26]=MMt*pow(armbbbarGL[3],2);
    armbbbarGL[27]=armbbbarGL[15]*armbbbarGL[26];
    armbbbarGL[21]=1./2.*armbbbarGL[21] - 6*armbbbarGL[27];
    armbbbarGL[21]=armbbbarGL[21]*armbbbarGL[23];
    armbbbarGL[21]=1./8.*armbbbarGL[24] + armbbbarGL[21];
    armbbbarGL[21]=MMt*armbbbarGL[21];
    armbbbarGL[24]=armbbbarGL[7]*MMH;
    armbbbarGL[27]= - 1./2.*armbbbarGL[24] - 9;
    armbbbarGL[27]=armbbbarGL[12]*armbbbarGL[27];
    armbbbarGL[22]=armbbbarGL[22] + armbbbarGL[27];
    armbbbarGL[25]=armbbbarGL[25] - armbbbarGL[26];
    armbbbarGL[25]=MMt*armbbbarGL[25];
    armbbbarGL[24]=63./4. - armbbbarGL[24];
    armbbbarGL[24]=1./16.*armbbbarGL[24] + 9*armbbbarGL[25];
    armbbbarGL[24]=armbbbarGL[13]*armbbbarGL[24];
    armbbbarGL[21]=1./2.*armbbbarGL[24] + 1./16.*armbbbarGL[22] + 
      armbbbarGL[21];
    armbbbarGL[21]=armbbbarGL[1]*armbbbarGL[21];
    armbbbarGL[22]=armbbbarGL[14]*armbbbarGL[1];
    armbbbarGL[24]=MMt*armbbbarGL[3];
    armbbbarGL[24]=5./8. - 2*armbbbarGL[24];
    armbbbarGL[24]=MMt*armbbbarGL[24];
    armbbbarGL[24]= - 1./32.*MMH + armbbbarGL[24];
    armbbbarGL[24]=armbbbarGL[24]*armbbbarGL[22];
    armbbbarGL[21]=3*armbbbarGL[24] + armbbbarGL[21];
    armbbbarGL[21]=armbbbarGL[13]*armbbbarGL[21];
    armbbbarGL[23]= - armbbbarGL[3]*armbbbarGL[23];
    armbbbarGL[23]=19./16. + armbbbarGL[23];
    armbbbarGL[23]=MMt*armbbbarGL[23];
    armbbbarGL[23]= - 7./64.*MMH + armbbbarGL[23];
    armbbbarGL[22]=MMt*armbbbarGL[23]*armbbbarGL[22];
    armbbbarGL[20]=armbbbarGL[21] + 1./2.*armbbbarGL[20] + 
      armbbbarGL[22];

    mbbbarGLret = armbbbarGL[20]*pow(armbbbarGL[4],2)*pow(
                                                          armbbbarGL[2],4);
    return mbbbarGLret.real();
  }
} // namespace mr
