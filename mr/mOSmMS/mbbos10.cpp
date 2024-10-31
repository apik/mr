#include <bb.hpp>
namespace mr
{
  double bb<MS>::x10(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<double> armbbos[22], mbbosret;

    armbbos[1]=double(nH);
    armbbos[2]=Tsil::A(mmt,mu2);
    armbbos[3]=pow(mmZ,-1);
    armbbos[4]=pow(mmH,-1);
    armbbos[5]=pow(s,-1);
    armbbos[6]=pow(c,-1);
    armbbos[7]=double(boson);
    armbbos[8]=Tsil::A(mmW,mu2);
    armbbos[9]=Tsil::A(mmZ,mu2);
    armbbos[10]=Tsil::A(mmH,mu2);
    armbbos[11]=Tsil::A(mmb,mu2);
    armbbos[12]=pow(mmb,-1);
    armbbos[13]=1/(mmt - mmW);
    armbbos[14]=armbbos[2] - armbbos[8];
    armbbos[15]=pow(armbbos[5],2);
    armbbos[16]=armbbos[15] - 1;
    armbbos[14]=armbbos[13]*armbbos[16]*armbbos[14];
    armbbos[14]=armbbos[14] + 1;
    armbbos[14]=mmZ*armbbos[14];
    armbbos[16]=armbbos[8] + mmZ;
    armbbos[16]=armbbos[16]*armbbos[15];
    armbbos[14]= - armbbos[16] + armbbos[14];
    armbbos[14]=armbbos[13]*armbbos[14];
    armbbos[17]=pow(armbbos[6],2);
    armbbos[18]=armbbos[17] + armbbos[15];
    armbbos[19]=armbbos[9]*armbbos[18];
    armbbos[20]=1./2.*armbbos[17];
    armbbos[21]=1 - armbbos[20];
    armbbos[21]=mmZ*armbbos[21];
    armbbos[16]= - 3./4.*armbbos[19] - 3./2.*armbbos[16] + armbbos[21];
    armbbos[16]=armbbos[4]*armbbos[16];
    armbbos[19]=1./2.*mmt + 3*armbbos[2];
    armbbos[19]= - 3./8.*armbbos[10] - 1./8.*armbbos[19];
    armbbos[19]=armbbos[18]*armbbos[19];
    armbbos[20]=1 + armbbos[20];
    armbbos[20]=armbbos[9]*armbbos[20];
    armbbos[19]=1./3.*armbbos[20] + armbbos[19];
    armbbos[19]=armbbos[3]*armbbos[19];
    armbbos[17]=1 + 31./24.*armbbos[17];
    armbbos[20]= - armbbos[12]*armbbos[11];
    armbbos[17]=1./2.*armbbos[17] + armbbos[20];
    armbbos[14]=armbbos[16] + 3./8.*armbbos[14] + 1./3.*armbbos[17] + 3./
      16.*armbbos[15] + armbbos[19];
    armbbos[14]=armbbos[7]*armbbos[14];
    armbbos[15]=armbbos[4]*armbbos[3]*armbbos[2]*armbbos[18]*armbbos[1]*
      mmt;

    mbbosret = armbbos[14] + 3*armbbos[15];
    return mbbosret.real();
  }
} // namespace mr
