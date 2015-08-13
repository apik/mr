#include <WW.hpp>
namespace mr
{
  long double WW<MS>::x11(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<long double> armWWos[22], mWWosret;

    armWWos[1]=double(nH);
    armWWos[2]=pow(mmZ,-1);
    armWWos[3]=pow(mmH,-1);
    armWWos[4]=pow(s,-1);
    armWWos[5]=pow(c,-1);
    armWWos[6]=Tsil::B(0,mmt,mmW,mu2);
    armWWos[7]=Tsil::A(mmt,mu2);
    armWWos[8]=pow(mmt,-1);
    armWWos[9]=Tsil::Aeps(mmt,mu2);
    armWWos[10]=prot00tt0->M(0);
    armWWos[11]=prot00tt0->Tuxv(0);
    armWWos[12]=double(nL);
    armWWos[13]=std::real(Tsil::B(0,0,mmW,mu2));
    armWWos[14]=prot00000->M(0);
    armWWos[15]=pow(mmt,2);
    armWWos[16]=armWWos[15]*armWWos[11];
    armWWos[17]=armWWos[9]*mmt;
    armWWos[16]=armWWos[16] + armWWos[17];
    armWWos[17]=armWWos[10]*pow(mmt,3);
    armWWos[16]= - armWWos[17] + 2*armWWos[16];
    armWWos[17]=5*armWWos[6];
    armWWos[17]=armWWos[17]*armWWos[15];
    armWWos[18]=armWWos[6]*mmt;
    armWWos[19]= - 6*armWWos[18] - 20./3.*armWWos[7] + 11./3.*mmt;
    armWWos[19]=armWWos[19]*armWWos[7];
    armWWos[16]=armWWos[19] + armWWos[17] + 2./3.*armWWos[16];
    armWWos[17]=2*armWWos[2];
    armWWos[17]=armWWos[16]*armWWos[17];
    armWWos[19]=armWWos[7]*armWWos[8];
    armWWos[20]= - 23 + 14*armWWos[6];
    armWWos[20]=2*armWWos[19] + 1./3.*armWWos[20];
    armWWos[21]=4*armWWos[7];
    armWWos[20]=armWWos[20]*armWWos[21];
    armWWos[21]=mmt + 3*armWWos[7];
    armWWos[21]=armWWos[21]*armWWos[7];
    armWWos[15]=armWWos[21] + 3*armWWos[15];
    armWWos[21]=16*armWWos[3];
    armWWos[15]=armWWos[15]*armWWos[21];
    armWWos[18]= - 5./3.*armWWos[18] + 7*mmt;
    armWWos[21]=2*armWWos[6];
    armWWos[18]=armWWos[18]*armWWos[21];
    armWWos[21]= - 17./2. + 8*armWWos[11];
    armWWos[21]=armWWos[21]*mmt;
    armWWos[21]=armWWos[21] + 8*armWWos[9];
    armWWos[15]= - 1./3.*armWWos[21] - armWWos[15] + armWWos[18] + 
      armWWos[17] - armWWos[20];
    armWWos[15]=armWWos[15]*armWWos[2];
    armWWos[17]=pow(armWWos[5],2);
    armWWos[16]=armWWos[16]*pow(armWWos[2],2)*armWWos[17];
    armWWos[16]=armWWos[15] + 2*armWWos[16];
    armWWos[16]=armWWos[16]*armWWos[17];
    armWWos[17]=1 - armWWos[6];
    armWWos[17]=armWWos[17]*armWWos[19];
    armWWos[18]=armWWos[8]*armWWos[9];
    armWWos[17]=armWWos[11] + armWWos[18] - armWWos[17];
    armWWos[18]=armWWos[10]*mmt;
    armWWos[19]= - 1 - 2./3.*armWWos[6];
    armWWos[19]=armWWos[6]*armWWos[19];
    armWWos[18]=armWWos[18] + armWWos[19];
    armWWos[19]=8./3.*mmZ;
    armWWos[20]=armWWos[19]*armWWos[10];
    armWWos[15]=armWWos[15] - armWWos[20] - 13 + 4*armWWos[18] - 16./3.*
      armWWos[17];
    armWWos[17]=pow(armWWos[4],2);
    armWWos[15]=armWWos[15]*armWWos[17];
    armWWos[15]=armWWos[15] + armWWos[20] + armWWos[16];
    armWWos[15]=armWWos[1]*armWWos[15];
    armWWos[16]=armWWos[19]*armWWos[14];
    armWWos[18]= - armWWos[16] - 31./3. - 4*armWWos[13];
    armWWos[17]=armWWos[17]*armWWos[18];
    armWWos[16]=armWWos[16] + armWWos[17];
    armWWos[16]=armWWos[12]*armWWos[16];

    mWWosret = armWWos[15] + armWWos[16];
    return mWWosret.real();
  }
} // namespace mr
