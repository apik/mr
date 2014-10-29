#include <bb.hpp>
std::complex<long double>
bb::m11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armbbbar[39], mbbbarret;

    armbbbar[1]=double(boson);
    armbbbar[2]=pow(CW,-1);
    armbbbar[3]=pow(MMH,-1);
    armbbbar[4]=pow(MMZ,-1);
    armbbbar[5]=pow(SW,-1);
    armbbbar[6]=Tsil::I2(0,MMW,MMt,mu2);
    armbbbar[7]=Tsil::A(MMH,mu2);
    armbbbar[8]=Tsil::A(MMb,mu2);
    armbbbar[9]=pow(MMb,-1);
    armbbbar[10]=Tsil::A(MMZ,mu2);
    armbbbar[11]=Tsil::A(MMW,mu2);
    armbbbar[12]=Tsil::A(MMt,mu2);
    armbbbar[13]=Tsil::Aeps(MMW,mu2);
    armbbbar[14]=Tsil::Aeps(MMt,mu2);
    armbbbar[15]=Tsil::Aeps(MMb,mu2);
    armbbbar[16]=prot0bb0b->Tvxu(0);
    armbbbar[17]=1/(MMt - MMW);
   armbbbar[18]=armbbbar[7]*armbbbar[4];
   armbbbar[19]=MMt*armbbbar[4];
   armbbbar[20]=1./4.*armbbbar[19] + 3./2.*armbbbar[18];
   armbbbar[21]=pow(armbbbar[2],2);
   armbbbar[22]=pow(armbbbar[5],2);
   armbbbar[23]=armbbbar[21] + armbbbar[22];
   armbbbar[20]=armbbbar[23]*armbbbar[20];
   armbbbar[24]= - armbbbar[12] + armbbbar[11];
   armbbbar[25]=armbbbar[17]*MMZ;
   armbbbar[26]=armbbbar[22] - 1;
   armbbbar[24]=armbbbar[25]*armbbbar[26]*armbbbar[24];
   armbbbar[27]=armbbbar[26]*MMZ;
   armbbbar[28]=armbbbar[22]*armbbbar[11];
   armbbbar[24]=armbbbar[24] + armbbbar[28] + armbbbar[27];
   armbbbar[24]=armbbbar[17]*armbbbar[24];
   armbbbar[29]=armbbbar[23]*armbbbar[12];
   armbbbar[30]=4*armbbbar[4];
   armbbbar[31]=armbbbar[30]*MMt*armbbbar[3];
   armbbbar[32]=1./2.*armbbbar[4] - armbbbar[31];
   armbbbar[32]=armbbbar[32]*armbbbar[29];
   armbbbar[33]=3*armbbbar[22];
   armbbbar[34]=armbbbar[21] - 2 + armbbbar[33];
   armbbbar[35]=2*MMZ;
   armbbbar[36]=armbbbar[35]*armbbbar[3];
   armbbbar[34]=armbbbar[34]*armbbbar[36];
   armbbbar[37]=armbbbar[3]*armbbbar[23];
   armbbbar[38]= - 2 - armbbbar[21];
   armbbbar[38]=armbbbar[4]*armbbbar[38];
   armbbbar[37]=2./3.*armbbbar[38] + 3*armbbbar[37];
   armbbbar[37]=armbbbar[10]*armbbbar[37];
   armbbbar[38]=armbbbar[3]*armbbbar[28];
   armbbbar[20]=3./2.*armbbbar[24] + armbbbar[37] + armbbbar[34] + 6*
   armbbbar[38] + 3*armbbbar[32] - 31./36.*armbbbar[21] + 34./9. - 3./4.
   *armbbbar[22] + armbbbar[20];
   armbbbar[20]=armbbbar[9]*armbbbar[20];
   armbbbar[24]=armbbbar[8]*pow(armbbbar[9],2);
   armbbbar[20]=armbbbar[20] + 32./9.*armbbbar[24];
   armbbbar[20]=armbbbar[8]*armbbbar[20];
   armbbbar[24]= - armbbbar[14] + armbbbar[6] - armbbbar[13];
   armbbbar[32]= - 7*armbbbar[24] + 3*armbbbar[11] + 11*armbbbar[12];
   armbbbar[32]=armbbbar[26]*armbbbar[32];
   armbbbar[34]=pow(CW,2);
   armbbbar[37]=armbbbar[26] - armbbbar[34];
   armbbbar[38]=MMZ*armbbbar[37];
   armbbbar[32]= - 17*armbbbar[38] + armbbbar[32];
   armbbbar[32]=MMZ*armbbbar[32];
   armbbbar[38]=armbbbar[11] + armbbbar[12] - armbbbar[24];
   armbbbar[37]=armbbbar[37]*armbbbar[38];
   armbbbar[38]=armbbbar[34] + 1;
   armbbbar[34]=armbbbar[38]*armbbbar[34];
   armbbbar[34]=armbbbar[34] - armbbbar[26];
   armbbbar[34]=armbbbar[34]*armbbbar[35];
   armbbbar[34]=armbbbar[34] + armbbbar[37];
   armbbbar[34]=MMZ*armbbbar[34];
   armbbbar[26]= - armbbbar[26]*armbbbar[11]*armbbbar[12];
   armbbbar[26]=armbbbar[26] + armbbbar[34];
   armbbbar[25]=armbbbar[26]*armbbbar[25];
   armbbbar[26]=armbbbar[22]*pow(armbbbar[12],2);
   armbbbar[34]=armbbbar[12]*armbbbar[28];
   armbbbar[25]=3*armbbbar[25] + armbbbar[32] - 3./2.*armbbbar[26] - 4*
   armbbbar[34];
   armbbbar[25]=armbbbar[17]*armbbbar[25];
   armbbbar[26]= - armbbbar[22]*armbbbar[24];
   armbbbar[32]=1./2.*armbbbar[23];
   armbbbar[34]=armbbbar[12]*armbbbar[4];
   armbbbar[35]= - armbbbar[32]*armbbbar[34];
   armbbbar[35]=9*armbbbar[22] + armbbbar[35];
   armbbbar[35]=armbbbar[12]*armbbbar[35];
   armbbbar[37]= - armbbbar[30]*armbbbar[29];
   armbbbar[33]=armbbbar[33] + armbbbar[37];
   armbbbar[33]=armbbbar[11]*armbbbar[33];
   armbbbar[25]=armbbbar[25] - 30*armbbbar[27] + armbbbar[33] + 8*
   armbbbar[26] + armbbbar[35];
   armbbbar[25]=armbbbar[17]*armbbbar[25];
   armbbbar[19]=armbbbar[3]*armbbbar[19];
   armbbbar[19]= - 103./12.*armbbbar[4] + 32*armbbbar[19];
   armbbbar[19]=MMt*armbbbar[19];
   armbbbar[26]= - 3./2.*armbbbar[4] - armbbbar[3];
   armbbbar[26]=armbbbar[10]*armbbbar[26];
   armbbbar[19]=armbbbar[26] + armbbbar[19];
   armbbbar[19]=armbbbar[23]*armbbbar[19];
   armbbbar[23]= - armbbbar[24]*armbbbar[23]*armbbbar[30];
   armbbbar[24]=armbbbar[3]*armbbbar[34];
   armbbbar[24]= - 24*armbbbar[24] + 13./2.*armbbbar[4] + armbbbar[31];
   armbbbar[24]=armbbbar[24]*armbbbar[29];
   armbbbar[18]= - armbbbar[18]*armbbbar[32];
   armbbbar[26]=2*armbbbar[3];
   armbbbar[26]= - armbbbar[26]*armbbbar[28];
   armbbbar[27]= - 1./3.*armbbbar[21] + 2./3. - armbbbar[22];
   armbbbar[27]=armbbbar[27]*armbbbar[36];
   armbbbar[28]=armbbbar[9]*armbbbar[15];
   armbbbar[28]=armbbbar[28] + armbbbar[16];
   armbbbar[18]=armbbbar[20] + armbbbar[25] + armbbbar[27] + 
   armbbbar[26] + armbbbar[24] + armbbbar[18] - 107./216.*armbbbar[21]
    - 154./27. - 167./8.*armbbbar[22] + armbbbar[23] - 40./9.*
   armbbbar[28] + armbbbar[19];

      mbbbarret = armbbbar[18]*armbbbar[1];
      return mbbbarret;
}
