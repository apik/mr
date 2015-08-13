#include <ZZ.hpp>
long double ZZ<OS>::x11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armZZbar[29], mZZbarret;

    armZZbar[1]=double(nH);
    armZZbar[2]=pow(CW,-1);
    armZZbar[3]=pow(MMH,-1);
    armZZbar[4]=pow(MMZ,-1);
    armZZbar[5]=pow(SW,-1);
    armZZbar[6]=Tsil::B(MMt,MMt,MMZ,mu2);
    armZZbar[7]=Tsil::A(MMt,mu2);
    armZZbar[8]=Tsil::Beps(MMt,MMt,MMZ,mu2);
    armZZbar[9]=Tsil::Aeps(MMt,mu2);
    armZZbar[10]=std::real(Tsil::B(0,0,MMZ,mu2));
    armZZbar[11]=prottttt0->M(0);
    armZZbar[12]=prot00000->M(0);
    armZZbar[13]=prottttt0->Vzxyv(0);
    armZZbar[14]=prottttt0->Suxv(0);
    armZZbar[15]=double(nL);
    armZZbar[16]=1/(4*MMt - MMZ);
   armZZbar[17]=pow(armZZbar[5],2);
   armZZbar[18]=pow(armZZbar[2],2);
   armZZbar[19]= - 7./9.*armZZbar[18] + armZZbar[17] + 64./9.;
   armZZbar[20]=2./3.*armZZbar[13] + 1./3.*armZZbar[11];
   armZZbar[20]=armZZbar[19]*armZZbar[20];
   armZZbar[21]=armZZbar[17] + armZZbar[18];
   armZZbar[22]=armZZbar[21]*armZZbar[3];
   armZZbar[20]=8*armZZbar[22] + armZZbar[20];
   armZZbar[23]=4*MMt;
   armZZbar[20]=armZZbar[20]*armZZbar[23];
   armZZbar[24]= - 299./9.*armZZbar[18] - 64./9. - 35*armZZbar[17];
   armZZbar[25]=armZZbar[8]*armZZbar[19];
   armZZbar[20]=2./3.*armZZbar[25] + 1./3.*armZZbar[24] + armZZbar[20];
   armZZbar[20]=MMt*armZZbar[20];
   armZZbar[24]=211./9.*armZZbar[18] - 448./9. + 11*armZZbar[17];
   armZZbar[24]=armZZbar[7]*armZZbar[24];
   armZZbar[25]=7*armZZbar[17];
   armZZbar[26]=95./9.*armZZbar[18] - 128./9. + armZZbar[25];
   armZZbar[26]=armZZbar[9]*armZZbar[26];
   armZZbar[24]=armZZbar[24] + armZZbar[26];
   armZZbar[22]=pow(armZZbar[7],2)*armZZbar[22];
   armZZbar[26]=17./9.*armZZbar[18] + armZZbar[17] - 32./9.;
   armZZbar[27]= - armZZbar[14]*armZZbar[26];
   armZZbar[22]=armZZbar[27] - 12*armZZbar[22] + 1./3.*armZZbar[24];
   armZZbar[24]=2*armZZbar[7];
   armZZbar[27]= - armZZbar[19]*armZZbar[24];
   armZZbar[28]=29./3.*armZZbar[18] - 128./3. - armZZbar[17];
   armZZbar[28]=MMt*armZZbar[28];
   armZZbar[27]=armZZbar[27] + armZZbar[28];
   armZZbar[19]= - armZZbar[6]*MMt*armZZbar[19];
   armZZbar[19]=2*armZZbar[27] + armZZbar[19];
   armZZbar[19]=armZZbar[6]*armZZbar[19];
   armZZbar[27]=armZZbar[7] + MMt;
   armZZbar[27]=MMt*armZZbar[27];
   armZZbar[28]=armZZbar[6]*pow(MMt,2);
   armZZbar[27]=2*armZZbar[27] + armZZbar[28];
   armZZbar[27]=armZZbar[6]*armZZbar[21]*armZZbar[27];
   armZZbar[24]=armZZbar[24] + armZZbar[9];
   armZZbar[24]=2*armZZbar[24] - MMt - armZZbar[14];
   armZZbar[24]=MMt*armZZbar[21]*armZZbar[24];
   armZZbar[24]=armZZbar[24] + armZZbar[27];
   armZZbar[24]=armZZbar[4]*armZZbar[24];
   armZZbar[19]=2./3.*armZZbar[24] + 1./3.*armZZbar[19] + 2*
   armZZbar[22] + armZZbar[20];
   armZZbar[19]=armZZbar[4]*armZZbar[19];
   armZZbar[20]=25./9.*armZZbar[18] + armZZbar[17] - 64./9.;
   armZZbar[22]= - armZZbar[9] + armZZbar[7];
   armZZbar[24]=4*armZZbar[16];
   armZZbar[22]=armZZbar[24]*armZZbar[20]*armZZbar[22];
   armZZbar[21]= - armZZbar[11]*armZZbar[21];
   armZZbar[27]=armZZbar[13]*armZZbar[26];
   armZZbar[21]=armZZbar[21] - 4./3.*armZZbar[27];
   armZZbar[21]=armZZbar[21]*armZZbar[23];
   armZZbar[23]= - armZZbar[7]*armZZbar[24];
   armZZbar[23]= - armZZbar[6] + armZZbar[23];
   armZZbar[23]=armZZbar[20]*armZZbar[23];
   armZZbar[23]=143./9.*armZZbar[18] - 320./9. + armZZbar[25] + 
   armZZbar[23];
   armZZbar[23]=armZZbar[6]*armZZbar[23];
   armZZbar[24]=armZZbar[20]*armZZbar[8];
   armZZbar[25]=937./18.*armZZbar[18] - 860./9. + 77./2.*armZZbar[17];
   armZZbar[27]=5./9.*armZZbar[18] + armZZbar[17] - 8./9.;
   armZZbar[28]=armZZbar[10]*armZZbar[27];
   armZZbar[19]=2*armZZbar[19] + armZZbar[23] - armZZbar[24] + 
   armZZbar[21] + 2*armZZbar[28] + 1./3.*armZZbar[25] + armZZbar[22];
   armZZbar[19]=armZZbar[1]*armZZbar[19];
   armZZbar[20]= - armZZbar[6]*armZZbar[20];
   armZZbar[20]=armZZbar[20] + 25./3.*armZZbar[18] - 64./3. + 3*
   armZZbar[17];
   armZZbar[20]=armZZbar[6]*armZZbar[20];
   armZZbar[20]=armZZbar[20] - armZZbar[24];
   armZZbar[20]=armZZbar[16]*armZZbar[20];
   armZZbar[21]=armZZbar[11]*armZZbar[26];
   armZZbar[22]=armZZbar[12]*armZZbar[27];
   armZZbar[21]=armZZbar[21] + armZZbar[22];
   armZZbar[20]=4./3.*armZZbar[21] + armZZbar[20];
   armZZbar[20]=armZZbar[1]*armZZbar[20];
   armZZbar[17]=11./9.*armZZbar[18] + armZZbar[17] - 20./9.;
   armZZbar[17]=armZZbar[17]*armZZbar[15];
   armZZbar[18]=armZZbar[12]*armZZbar[17];
   armZZbar[18]=8./3.*armZZbar[18] + armZZbar[20];
   armZZbar[18]=MMZ*armZZbar[18];
   armZZbar[20]=31./3. + 4*armZZbar[10];
   armZZbar[17]=armZZbar[20]*armZZbar[17];

      mZZbarret = armZZbar[17] + armZZbar[18] + armZZbar[19];
      return mZZbarret.real();
}
