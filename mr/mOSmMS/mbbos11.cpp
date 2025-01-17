#include <bb.hpp>
namespace mr
{
  double bb<MS>::x11(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<double> armbbos[46], mbbosret;

    armbbos[1]=double(nH);
    armbbos[2]=pow(mmZ,-1);
    armbbos[3]=pow(mmH,-1);
    armbbos[4]=pow(s,-1);
    armbbos[5]=pow(c,-1);
    armbbos[6]=Tsil::A(mmt,mu2);
    armbbos[7]=Tsil::A(mmb,mu2);
    armbbos[8]=pow(mmb,-1);
    armbbos[9]=double(boson);
    armbbos[10]=Tsil::I2(0,mmW,mmt,mu2);
    armbbos[11]=Tsil::A(mmW,mu2);
    armbbos[12]=Tsil::A(mmZ,mu2);
    armbbos[13]=Tsil::A(mmH,mu2);
    armbbos[14]=Tsil::Aeps(mmW,mu2);
    armbbos[15]=Tsil::Aeps(mmt,mu2);
    armbbos[16]=Tsil::Aeps(mmb,mu2);
    armbbos[17]=prot0bb0b->Suxv(0);
    armbbos[18]=1/(mmt - mmW);
    armbbos[19]=pow(armbbos[5],2);
    armbbos[20]=armbbos[19]*armbbos[3];
    armbbos[21]=pow(armbbos[4],2);
    armbbos[22]= - 2./3. + armbbos[21];
    armbbos[22]=armbbos[3]*armbbos[22];
    armbbos[22]=armbbos[22] + 1./3.*armbbos[20];
    armbbos[23]= - 2 + 3*armbbos[21];
    armbbos[23]=armbbos[3]*armbbos[23];
    armbbos[23]=armbbos[23] + armbbos[20];
    armbbos[24]=armbbos[8]*armbbos[7];
    armbbos[23]=armbbos[23]*armbbos[24];
    armbbos[22]=5*armbbos[22] + armbbos[23];
    armbbos[23]=pow(armbbos[18],2);
    armbbos[25]=pow(armbbos[18],3);
    armbbos[26]=armbbos[10]*armbbos[25];
    armbbos[26]=19*armbbos[23] + 3*armbbos[26];
    armbbos[27]=armbbos[21] - 1;
    armbbos[28]=pow(c,2);
    armbbos[29]=armbbos[27] - armbbos[28];
    armbbos[26]=armbbos[29]*armbbos[26];
    armbbos[30]=armbbos[28] + 1;
    armbbos[28]= - armbbos[30]*armbbos[28];
    armbbos[28]=armbbos[28] + armbbos[27];
    armbbos[25]=armbbos[25]*mmZ;
    armbbos[28]=armbbos[28]*armbbos[25];
    armbbos[26]=6*armbbos[28] + armbbos[26];
    armbbos[26]=mmZ*armbbos[26];
    armbbos[28]=3./2.*armbbos[24];
    armbbos[30]=armbbos[28]*armbbos[27];
    armbbos[31]=34*armbbos[27] + armbbos[30];
    armbbos[31]=armbbos[18]*armbbos[31];
    armbbos[32]=armbbos[27]*armbbos[23];
    armbbos[32]=7*armbbos[32];
    armbbos[33]=armbbos[10]*armbbos[32];
    armbbos[22]=armbbos[26] + armbbos[33] + 2*armbbos[22] + armbbos[31];
    armbbos[22]=mmZ*armbbos[22];
    armbbos[26]=armbbos[23]*armbbos[21];
    armbbos[31]=3*armbbos[25];
    armbbos[33]= - armbbos[27]*armbbos[31];
    armbbos[34]=armbbos[19] + armbbos[21];
    armbbos[35]=armbbos[34]*armbbos[2];
    armbbos[36]=4*armbbos[35];
    armbbos[37]=armbbos[18]*armbbos[36];
    armbbos[33]=armbbos[37] + armbbos[26] + armbbos[33];
    armbbos[33]=armbbos[6]*armbbos[33];
    armbbos[37]=armbbos[29]*armbbos[25];
    armbbos[38]=2*armbbos[27];
    armbbos[39]=armbbos[38] + armbbos[30];
    armbbos[39]=armbbos[39]*armbbos[23];
    armbbos[39]=armbbos[39] - armbbos[37];
    armbbos[39]=mmZ*armbbos[39];
    armbbos[40]=3*armbbos[24];
    armbbos[41]=armbbos[40] + 5;
    armbbos[42]=armbbos[21]*armbbos[3];
    armbbos[43]=armbbos[42]*armbbos[41];
    armbbos[44]=armbbos[21]*armbbos[18];
    armbbos[45]=armbbos[44]*armbbos[28];
    armbbos[33]=armbbos[33] + armbbos[39] + 2*armbbos[43] + armbbos[45];
    armbbos[33]=armbbos[11]*armbbos[33];
    armbbos[39]=armbbos[34]*armbbos[13];
    armbbos[43]=armbbos[34]*mmt;
    armbbos[45]=3*armbbos[39] + 1./2.*armbbos[43];
    armbbos[45]=armbbos[45]*armbbos[24];
    armbbos[39]=armbbos[45] + 5*armbbos[39] + 31./2.*armbbos[43];
    armbbos[43]= - 2 - armbbos[19];
    armbbos[24]=armbbos[43]*armbbos[24];
    armbbos[24]=2./3.*armbbos[24] + 11./18.*armbbos[19] - 16./9. + 3./2.
      *armbbos[21];
    armbbos[24]=armbbos[12]*armbbos[24];
    armbbos[43]=4*armbbos[34];
    armbbos[43]=armbbos[10]*armbbos[43];
    armbbos[24]=armbbos[43] + 1./2.*armbbos[39] + armbbos[24];
    armbbos[24]=armbbos[2]*armbbos[24];
    armbbos[27]= - 20*armbbos[27] - armbbos[30];
    armbbos[23]=armbbos[27]*armbbos[23];
    armbbos[23]=armbbos[23] - 5*armbbos[37];
    armbbos[23]=mmZ*armbbos[23];
    armbbos[25]=armbbos[25]*armbbos[38];
    armbbos[25]= - 1./2.*armbbos[26] + armbbos[25];
    armbbos[26]=armbbos[35]*armbbos[18];
    armbbos[25]=3*armbbos[25] + 7./2.*armbbos[26];
    armbbos[25]=armbbos[6]*armbbos[25];
    armbbos[26]= - 2 + armbbos[28];
    armbbos[26]=armbbos[26]*armbbos[35];
    armbbos[23]=armbbos[25] + armbbos[26] - 9*armbbos[44] + armbbos[23];
    armbbos[23]=armbbos[6]*armbbos[23];
    armbbos[25]=armbbos[31]*armbbos[29];
    armbbos[25]=armbbos[25] + armbbos[32];
    armbbos[25]=armbbos[25]*mmZ;
    armbbos[26]=8*armbbos[44];
    armbbos[25]=armbbos[25] + armbbos[36] + armbbos[26];
    armbbos[27]= - armbbos[14] - armbbos[15];
    armbbos[25]=armbbos[25]*armbbos[27];
    armbbos[27]=2 - 1./4.*armbbos[21];
    armbbos[27]=3*armbbos[27] - 31./36.*armbbos[19];
    armbbos[27]=armbbos[7]*armbbos[27];
    armbbos[28]=armbbos[16] - armbbos[17];
    armbbos[29]=armbbos[8]*pow(armbbos[7],2);
    armbbos[27]=armbbos[27] - 8./9.*armbbos[29] + 40./9.*armbbos[28];
    armbbos[27]=armbbos[8]*armbbos[27];
    armbbos[20]=armbbos[42] + armbbos[20];
    armbbos[20]=armbbos[12]*armbbos[20]*armbbos[41];
    armbbos[26]=armbbos[10]*armbbos[26];
    armbbos[19]=armbbos[33] + armbbos[25] + armbbos[23] + armbbos[24] + 
      armbbos[22] + armbbos[26] + armbbos[20] - 47./72.*armbbos[19] + 61./
      27. + 159./8.*armbbos[21] + armbbos[27];
    armbbos[19]=armbbos[9]*armbbos[19];
    armbbos[20]=armbbos[2]*armbbos[34]*armbbos[1]*armbbos[3];
    armbbos[21]= - 7 - armbbos[40];
    armbbos[21]=armbbos[21]*mmt*armbbos[20];
    armbbos[20]=6*armbbos[20];
    armbbos[22]= - armbbos[6]*armbbos[20];
    armbbos[21]=armbbos[21] + armbbos[22];
    armbbos[21]=armbbos[6]*armbbos[21];
    armbbos[20]= - pow(mmt,2)*armbbos[20];
    armbbos[20]=armbbos[20] + armbbos[21];

    mbbosret = armbbos[19] + 4*armbbos[20];
    return mbbosret.real();
  }
} // namespace mr
