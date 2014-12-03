#include <ZZ.hpp>
long double ZZ<MS>::x11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armZZos[30], mZZosret;

    armZZos[1]=double(nH);
    armZZos[2]=pow(mmZ,-1);
    armZZos[3]=pow(s,-1);
    armZZos[4]=pow(c,-1);
    armZZos[5]=pow(mmH,-1);
    armZZos[6]=Tsil::B(mmt,mmt,mmZ,mu2);
    armZZos[7]=Tsil::A(mmt,mu2);
    armZZos[8]=Tsil::Beps(mmt,mmt,mmZ,mu2);
    armZZos[9]=Tsil::Aeps(mmt,mu2);
    armZZos[10]=std::real(Tsil::B(0,0,mmZ,mu2));
    armZZos[11]=prottttt0->M(0);
    armZZos[12]=prot00000->M(0);
    armZZos[13]=prottttt0->Vzxyv(0);
    armZZos[14]=prottttt0->Suxv(0);
    armZZos[15]=double(nL);
    armZZos[16]=1/(4*mmt - mmZ);
   armZZos[17]=pow(armZZos[3],2);
   armZZos[18]=pow(armZZos[4],2);
   armZZos[19]= - 7./9.*armZZos[18] + armZZos[17] + 64./9.;
   armZZos[20]= - armZZos[8]*armZZos[19];
   armZZos[20]=armZZos[20] + 85./9.*armZZos[18] + 128./9. + 13*
   armZZos[17];
   armZZos[20]=armZZos[2]*armZZos[20];
   armZZos[21]=armZZos[17] + armZZos[18];
   armZZos[22]=pow(armZZos[2],2);
   armZZos[23]=armZZos[21]*armZZos[22];
   armZZos[24]=armZZos[14]*armZZos[23];
   armZZos[20]=armZZos[20] + armZZos[24];
   armZZos[19]=armZZos[19]*armZZos[2];
   armZZos[24]=armZZos[6]*armZZos[19];
   armZZos[25]= - 13./3.*armZZos[18] + 64./3. + armZZos[17];
   armZZos[25]=armZZos[2]*armZZos[25];
   armZZos[25]=5*armZZos[25] + armZZos[24];
   armZZos[25]=armZZos[6]*armZZos[25];
   armZZos[26]=armZZos[6] + 2;
   armZZos[26]=armZZos[26]*armZZos[23];
   armZZos[27]=armZZos[7]*armZZos[26];
   armZZos[20]= - 4*armZZos[27] + armZZos[25] + 2*armZZos[20];
   armZZos[25]= - armZZos[6]*armZZos[26];
   armZZos[23]=armZZos[23] + armZZos[25];
   armZZos[25]=armZZos[5]*armZZos[2];
   armZZos[26]=armZZos[21]*armZZos[25];
   armZZos[27]= - 4./3.*armZZos[13] - 2./3.*armZZos[11];
   armZZos[19]=armZZos[27]*armZZos[19];
   armZZos[19]=1./3.*armZZos[23] - 12*armZZos[26] + armZZos[19];
   armZZos[23]=2*mmt;
   armZZos[19]=armZZos[19]*armZZos[23];
   armZZos[27]=4./3.*armZZos[9];
   armZZos[22]= - armZZos[22]*armZZos[27];
   armZZos[25]=armZZos[7]*armZZos[25];
   armZZos[22]=armZZos[22] + 2*armZZos[11] - 8*armZZos[25];
   armZZos[21]=armZZos[21]*armZZos[22];
   armZZos[22]=17./9.*armZZos[18] + armZZos[17] - 32./9.;
   armZZos[25]=armZZos[13]*armZZos[22];
   armZZos[19]=armZZos[19] + 8./3.*armZZos[25] + 1./3.*armZZos[20] + 
   armZZos[21];
   armZZos[19]=armZZos[19]*armZZos[23];
   armZZos[20]=1./2.*armZZos[17] - 32./9. + 25./18.*armZZos[18];
   armZZos[21]=25./9.*armZZos[18] + armZZos[17] - 64./9.;
   armZZos[23]=armZZos[21]*armZZos[6];
   armZZos[25]= - 7*armZZos[20] + armZZos[23];
   armZZos[25]=armZZos[6]*armZZos[25];
   armZZos[28]=armZZos[21]*armZZos[8];
   armZZos[20]=armZZos[25] + armZZos[28] + armZZos[20];
   armZZos[20]=mmZ*armZZos[20];
   armZZos[25]=3*armZZos[17] - 64./3. + 25./3.*armZZos[18];
   armZZos[29]=armZZos[7]*armZZos[2];
   armZZos[29]=4*armZZos[29] - 2;
   armZZos[25]=armZZos[25]*armZZos[29];
   armZZos[25]=5*armZZos[23] + armZZos[25];
   armZZos[25]=armZZos[7]*armZZos[25];
   armZZos[21]=armZZos[9]*armZZos[21];
   armZZos[20]=4*armZZos[21] + 2*armZZos[25] + armZZos[20];
   armZZos[20]=armZZos[16]*armZZos[20];
   armZZos[21]=armZZos[14]*armZZos[22];
   armZZos[25]= - 95./9.*armZZos[18] + 128./9. - 7*armZZos[17];
   armZZos[25]=armZZos[25]*armZZos[27];
   armZZos[21]=armZZos[25] + 4*armZZos[21];
   armZZos[21]=armZZos[2]*armZZos[21];
   armZZos[23]=armZZos[23] - 311./18.*armZZos[18] + 352./9. - 15./2.*
   armZZos[17];
   armZZos[23]=armZZos[6]*armZZos[23];
   armZZos[25]= - 55./9.*armZZos[18] + 256./9. + armZZos[17];
   armZZos[25]=armZZos[2]*armZZos[25];
   armZZos[24]=2*armZZos[25] - 5*armZZos[24];
   armZZos[24]=armZZos[7]*armZZos[24];
   armZZos[25]=armZZos[17] - 8./9. + 5./9.*armZZos[18];
   armZZos[27]=armZZos[12]*mmZ;
   armZZos[29]= - 4./3.*armZZos[27] - 2*armZZos[10];
   armZZos[25]=armZZos[25]*armZZos[29];
   armZZos[29]= - 431./9.*armZZos[18] + 764./9. - 37*armZZos[17];
   armZZos[26]=pow(armZZos[7],2)*armZZos[26];
   armZZos[22]=armZZos[11]*mmZ*armZZos[22];
   armZZos[19]=armZZos[19] + armZZos[20] - 4./3.*armZZos[22] - 48*
   armZZos[26] + 2./3.*armZZos[24] + armZZos[23] + 1./3.*armZZos[29] + 
   armZZos[28] + armZZos[21] + armZZos[25];
   armZZos[19]=armZZos[1]*armZZos[19];
   armZZos[20]=11./9.*armZZos[15];
   armZZos[18]=armZZos[20]*armZZos[18];
   armZZos[17]=armZZos[17] - 20./9.;
   armZZos[17]=armZZos[17]*armZZos[15];
   armZZos[17]=armZZos[17] + armZZos[18];
   armZZos[18]= - 8./3.*armZZos[27] - 31./3. - 4*armZZos[10];
   armZZos[17]=armZZos[17]*armZZos[18];

      mZZosret = armZZos[17] + armZZos[19];
      return mZZosret.real();
}
