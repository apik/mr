#include <WW.hpp>
std::complex<long double>
WW<MS>::m11(size_t nL, size_t nH, size_t boson)
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
   armWWos[15]=17./2. - 8*armWWos[11];
   armWWos[16]=7 - 5./3.*armWWos[6];
   armWWos[16]=armWWos[6]*armWWos[16];
   armWWos[17]= - armWWos[7]*armWWos[3];
   armWWos[15]=16*armWWos[17] + 1./3.*armWWos[15] + 2*armWWos[16];
   armWWos[16]=pow(armWWos[5],2);
   armWWos[17]=armWWos[16]*armWWos[15];
   armWWos[18]=pow(armWWos[4],2);
   armWWos[15]=armWWos[18]*armWWos[15];
   armWWos[19]= - armWWos[16]*armWWos[3];
   armWWos[20]= - armWWos[18]*armWWos[3];
   armWWos[19]=armWWos[19] + armWWos[20];
   armWWos[19]=mmt*armWWos[19];
   armWWos[15]=48*armWWos[19] + armWWos[17] + armWWos[15];
   armWWos[15]=mmt*armWWos[15];
   armWWos[17]=4./3.*armWWos[11] + 5*armWWos[6];
   armWWos[19]=armWWos[16]*armWWos[17];
   armWWos[19]=armWWos[17] + armWWos[19];
   armWWos[19]=armWWos[16]*armWWos[19];
   armWWos[17]=armWWos[18]*armWWos[17];
   armWWos[20]= - armWWos[16]*armWWos[10];
   armWWos[20]= - armWWos[10] + armWWos[20];
   armWWos[20]=armWWos[16]*armWWos[20];
   armWWos[21]= - armWWos[18]*armWWos[10];
   armWWos[20]=armWWos[20] + armWWos[21];
   armWWos[20]=mmt*armWWos[20];
   armWWos[17]=2./3.*armWWos[20] + armWWos[19] + armWWos[17];
   armWWos[17]=mmt*armWWos[17];
   armWWos[19]=11./3. - 6*armWWos[6];
   armWWos[19]=armWWos[7]*armWWos[19];
   armWWos[19]=4./3.*armWWos[9] + armWWos[19];
   armWWos[20]=armWWos[16]*armWWos[19];
   armWWos[20]=armWWos[19] + armWWos[20];
   armWWos[20]=armWWos[16]*armWWos[20];
   armWWos[19]=armWWos[18]*armWWos[19];
   armWWos[17]=armWWos[17] + armWWos[20] + armWWos[19];
   armWWos[17]=mmt*armWWos[17];
   armWWos[19]=pow(armWWos[7],2);
   armWWos[20]= - armWWos[16]*armWWos[19];
   armWWos[20]= - armWWos[19] + armWWos[20];
   armWWos[20]=armWWos[16]*armWWos[20];
   armWWos[19]= - armWWos[18]*armWWos[19];
   armWWos[19]=armWWos[20] + armWWos[19];
   armWWos[17]=20./3.*armWWos[19] + armWWos[17];
   armWWos[17]=armWWos[2]*armWWos[17];
   armWWos[19]=23 - 14*armWWos[6];
   armWWos[20]= - armWWos[8] - 6*armWWos[3];
   armWWos[20]=armWWos[7]*armWWos[20];
   armWWos[19]=1./3.*armWWos[19] + 2*armWWos[20];
   armWWos[19]=armWWos[7]*armWWos[19];
   armWWos[19]= - 2./3.*armWWos[9] + armWWos[19];
   armWWos[16]=armWWos[16]*armWWos[19];
   armWWos[19]=armWWos[18]*armWWos[19];
   armWWos[16]=armWWos[16] + armWWos[19];
   armWWos[15]=2*armWWos[17] + 4*armWWos[16] + armWWos[15];
   armWWos[15]=armWWos[2]*armWWos[15];
   armWWos[16]= - armWWos[10]*mmZ;
   armWWos[17]= - armWWos[9]*armWWos[8];
   armWWos[19]= - 1 - 2./3.*armWWos[6];
   armWWos[19]=armWWos[6]*armWWos[19];
   armWWos[20]= - armWWos[6]*armWWos[8];
   armWWos[20]=armWWos[8] + armWWos[20];
   armWWos[20]=armWWos[7]*armWWos[20];
   armWWos[16]=16./3.*armWWos[20] + 4*armWWos[19] + 16./3.*armWWos[17]
    + 8./3.*armWWos[16] - 13 - 16./3.*armWWos[11];
   armWWos[16]=armWWos[18]*armWWos[16];
   armWWos[17]=armWWos[10]*mmZ;
   armWWos[19]=mmt*armWWos[18]*armWWos[10];
   armWWos[15]=armWWos[15] + 4*armWWos[19] + 8./3.*armWWos[17] + 
   armWWos[16];
   armWWos[15]=armWWos[1]*armWWos[15];
   armWWos[16]= - mmZ*armWWos[14];
   armWWos[16]=8./3.*armWWos[16] - 31./3. - 4*armWWos[13];
   armWWos[16]=armWWos[18]*armWWos[12]*armWWos[16];
   armWWos[17]=armWWos[12]*mmZ*armWWos[14];

      mWWosret = armWWos[15] + armWWos[16] + 8./3.*armWWos[17];
      return mWWosret;
}
