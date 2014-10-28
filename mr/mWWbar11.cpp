#include <WW.hpp>
std::complex<long double>
WW<OS>::m11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armWWbar[24], mWWbarret;

    armWWbar[1]=double(nH);
    armWWbar[2]=pow(CW,-1);
    armWWbar[3]=pow(MMH,-1);
    armWWbar[4]=pow(MMZ,-1);
    armWWbar[5]=pow(SW,-1);
    armWWbar[6]=Tsil::B(0,MMt,MMW,mu2);
    armWWbar[7]=Tsil::A(MMt,mu2);
    armWWbar[8]=pow(MMt,-1);
    armWWbar[9]=Tsil::Aeps(MMt,mu2);
    armWWbar[10]=prot00tt0->M(0);
    armWWbar[11]=prot00tt0->Tuxv(0);
    armWWbar[12]=double(nL);
    armWWbar[13]=std::real(Tsil::B(0,0,MMW,mu2));
    armWWbar[14]=prot00000->M(0);
   armWWbar[15]=armWWbar[10]*MMZ;
   armWWbar[16]=armWWbar[9]*armWWbar[8];
   armWWbar[17]=1 + 2./3.*armWWbar[6];
   armWWbar[17]=armWWbar[6]*armWWbar[17];
   armWWbar[18]=armWWbar[6]*armWWbar[8];
   armWWbar[18]= - armWWbar[8] + armWWbar[18];
   armWWbar[18]=armWWbar[7]*armWWbar[18];
   armWWbar[19]= - MMt*armWWbar[10];
   armWWbar[15]=4*armWWbar[19] + 16./3.*armWWbar[18] + 4*armWWbar[17]
    + 16./3.*armWWbar[16] + 8./3.*armWWbar[15] + 13 + 16./3.*
   armWWbar[11];
   armWWbar[16]=pow(armWWbar[5],2);
   armWWbar[15]=armWWbar[16]*armWWbar[15];
   armWWbar[17]= - 4./3.*armWWbar[11] - 3*armWWbar[6];
   armWWbar[18]=pow(armWWbar[2],2);
   armWWbar[19]=armWWbar[18]*armWWbar[17];
   armWWbar[19]=armWWbar[17] + armWWbar[19];
   armWWbar[19]=armWWbar[18]*armWWbar[19];
   armWWbar[20]=armWWbar[18]*armWWbar[10];
   armWWbar[20]=armWWbar[10] + armWWbar[20];
   armWWbar[20]=MMt*armWWbar[18]*armWWbar[20];
   armWWbar[19]=armWWbar[19] + 2./3.*armWWbar[20];
   armWWbar[19]=MMt*armWWbar[19];
   armWWbar[20]= - 4*armWWbar[9] - 5*armWWbar[7];
   armWWbar[21]=armWWbar[18]*armWWbar[20];
   armWWbar[21]=armWWbar[20] + armWWbar[21];
   armWWbar[21]=armWWbar[18]*armWWbar[21];
   armWWbar[19]=1./3.*armWWbar[21] + armWWbar[19];
   armWWbar[19]=MMt*armWWbar[19];
   armWWbar[21]=pow(armWWbar[7],2);
   armWWbar[22]=armWWbar[18]*armWWbar[21];
   armWWbar[22]=armWWbar[21] + armWWbar[22];
   armWWbar[22]=armWWbar[18]*armWWbar[22];
   armWWbar[23]=MMt*armWWbar[10];
   armWWbar[17]=armWWbar[17] + 2./3.*armWWbar[23];
   armWWbar[17]=MMt*armWWbar[17];
   armWWbar[17]=1./3.*armWWbar[20] + armWWbar[17];
   armWWbar[17]=MMt*armWWbar[17];
   armWWbar[17]=2./3.*armWWbar[21] + armWWbar[17];
   armWWbar[17]=armWWbar[16]*armWWbar[17];
   armWWbar[17]=armWWbar[17] + 2./3.*armWWbar[22] + armWWbar[19];
   armWWbar[17]=armWWbar[4]*armWWbar[17];
   armWWbar[19]= - 1 + armWWbar[6];
   armWWbar[20]= - 6*armWWbar[3] + armWWbar[8];
   armWWbar[20]=armWWbar[7]*armWWbar[20];
   armWWbar[19]=5./3.*armWWbar[19] + 2*armWWbar[20];
   armWWbar[19]=armWWbar[7]*armWWbar[19];
   armWWbar[19]=2./3.*armWWbar[9] + armWWbar[19];
   armWWbar[20]= - 65./2. + 8*armWWbar[11];
   armWWbar[21]= - 1 + 1./3.*armWWbar[6];
   armWWbar[21]=armWWbar[6]*armWWbar[21];
   armWWbar[20]=1./3.*armWWbar[20] + 10*armWWbar[21];
   armWWbar[21]=MMt*armWWbar[3];
   armWWbar[21]=armWWbar[20] + 64*armWWbar[21];
   armWWbar[21]=MMt*armWWbar[21];
   armWWbar[21]=4*armWWbar[19] + armWWbar[21];
   armWWbar[21]=armWWbar[16]*armWWbar[21];
   armWWbar[19]=armWWbar[18]*armWWbar[19];
   armWWbar[20]=armWWbar[18]*armWWbar[20];
   armWWbar[18]=MMt*armWWbar[18]*armWWbar[3];
   armWWbar[18]=armWWbar[20] + 64*armWWbar[18];
   armWWbar[18]=MMt*armWWbar[18];
   armWWbar[17]=2*armWWbar[17] + armWWbar[21] + 4*armWWbar[19] + 
   armWWbar[18];
   armWWbar[17]=armWWbar[4]*armWWbar[17];
   armWWbar[18]= - armWWbar[10]*MMZ;
   armWWbar[15]=armWWbar[17] + 8./3.*armWWbar[18] + armWWbar[15];
   armWWbar[15]=armWWbar[1]*armWWbar[15];
   armWWbar[17]=MMZ*armWWbar[14];
   armWWbar[17]=8./3.*armWWbar[17] + 31./3. + 4*armWWbar[13];
   armWWbar[16]=armWWbar[16]*armWWbar[12]*armWWbar[17];
   armWWbar[17]= - armWWbar[12]*MMZ*armWWbar[14];

      mWWbarret = armWWbar[15] + armWWbar[16] + 8./3.*armWWbar[17];
      return mWWbarret;
}
