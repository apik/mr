#include <ZZ.hpp>
std::complex<long double>
ZZ<MS>::m11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armZZos[28], mZZosret;

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
   armZZos[17]= - 3*armZZos[5];
   armZZos[18]=7./27.*armZZos[13] + armZZos[17];
   armZZos[18]=2*armZZos[18] + 7./27.*armZZos[11];
   armZZos[19]=pow(armZZos[4],2);
   armZZos[18]=armZZos[19]*armZZos[18];
   armZZos[17]= - 1./3.*armZZos[13] + armZZos[17];
   armZZos[17]=2*armZZos[17] - 1./3.*armZZos[11];
   armZZos[20]=pow(armZZos[3],2);
   armZZos[17]=armZZos[20]*armZZos[17];
   armZZos[21]= - 2*armZZos[13] - armZZos[11];
   armZZos[17]=armZZos[17] + 64./27.*armZZos[21] + armZZos[18];
   armZZos[17]=mmt*armZZos[17];
   armZZos[18]=2 - armZZos[8];
   armZZos[21]=5 + 1./3.*armZZos[6];
   armZZos[21]=armZZos[6]*armZZos[21];
   armZZos[18]=2./3.*armZZos[18] + armZZos[21];
   armZZos[21]=85 + 7*armZZos[8];
   armZZos[22]= - armZZos[7]*armZZos[5];
   armZZos[23]=4*armZZos[22];
   armZZos[21]=1./27.*armZZos[21] + armZZos[23];
   armZZos[24]= - 65 - 7./3.*armZZos[6];
   armZZos[24]=armZZos[6]*armZZos[24];
   armZZos[21]=2*armZZos[21] + 1./9.*armZZos[24];
   armZZos[21]=armZZos[19]*armZZos[21];
   armZZos[24]=13 - armZZos[8];
   armZZos[23]=1./3.*armZZos[24] + armZZos[23];
   armZZos[24]=5 + armZZos[6];
   armZZos[24]=armZZos[6]*armZZos[24];
   armZZos[23]=2*armZZos[23] + 1./3.*armZZos[24];
   armZZos[23]=armZZos[20]*armZZos[23];
   armZZos[17]=4*armZZos[17] + armZZos[23] + 64./9.*armZZos[18] + 
   armZZos[21];
   armZZos[17]=mmt*armZZos[17];
   armZZos[18]= - 2 - armZZos[6];
   armZZos[18]=armZZos[6]*armZZos[18];
   armZZos[18]=1 + armZZos[18];
   armZZos[21]=armZZos[19]*armZZos[18];
   armZZos[18]=armZZos[20]*armZZos[18];
   armZZos[18]=armZZos[21] + armZZos[18];
   armZZos[18]=mmt*armZZos[18];
   armZZos[21]= - 4*armZZos[7];
   armZZos[23]= - armZZos[6]*armZZos[7];
   armZZos[24]=2*armZZos[23] + armZZos[21] + armZZos[14] - 2*armZZos[9]
   ;
   armZZos[25]=armZZos[19]*armZZos[24];
   armZZos[24]=armZZos[20]*armZZos[24];
   armZZos[18]=armZZos[18] + armZZos[25] + armZZos[24];
   armZZos[18]=armZZos[2]*mmt*armZZos[18];
   armZZos[24]=8./3.*armZZos[7] - armZZos[14] + 4./3.*armZZos[9];
   armZZos[25]=pow(armZZos[7],2);
   armZZos[26]= - armZZos[16]*armZZos[25];
   armZZos[24]=5./9.*armZZos[23] + 1./3.*armZZos[24] + 4*armZZos[26];
   armZZos[26]=17*armZZos[14] - 95./3.*armZZos[9];
   armZZos[22]=12*armZZos[22];
   armZZos[27]= - 55./27. + armZZos[22];
   armZZos[27]=armZZos[7]*armZZos[27];
   armZZos[25]=armZZos[16]*armZZos[25];
   armZZos[26]=50./3.*armZZos[25] + 1./9.*armZZos[26] + armZZos[27];
   armZZos[27]=armZZos[6]*armZZos[7];
   armZZos[26]=2*armZZos[26] + 35./27.*armZZos[27];
   armZZos[26]=armZZos[19]*armZZos[26];
   armZZos[22]=1./3. + armZZos[22];
   armZZos[22]=armZZos[7]*armZZos[22];
   armZZos[22]=6*armZZos[25] + armZZos[22] + armZZos[14] - 7./3.*
   armZZos[9];
   armZZos[22]=2*armZZos[22] + 5./3.*armZZos[23];
   armZZos[22]=armZZos[20]*armZZos[22];
   armZZos[17]=2./3.*armZZos[18] + armZZos[17] + armZZos[22] + 64./3.*
   armZZos[24] + armZZos[26];
   armZZos[17]=armZZos[2]*armZZos[17];
   armZZos[18]=4*armZZos[10];
   armZZos[22]=armZZos[12] + 4*armZZos[11];
   armZZos[22]=mmZ*armZZos[22];
   armZZos[22]=8./3.*armZZos[22] - 16*armZZos[8] + 191./3. + 
   armZZos[18];
   armZZos[23]=7*mmZ - 20*armZZos[7];
   armZZos[23]=armZZos[16]*armZZos[23];
   armZZos[24]= - armZZos[16]*mmZ;
   armZZos[24]= - 1 + armZZos[24];
   armZZos[24]=armZZos[6]*armZZos[24];
   armZZos[23]=2*armZZos[24] + 11 + armZZos[23];
   armZZos[23]=armZZos[6]*armZZos[23];
   armZZos[24]= - 1 - 2*armZZos[8];
   armZZos[24]=mmZ*armZZos[24];
   armZZos[24]= - 8*armZZos[9] + armZZos[24];
   armZZos[24]=1./3.*armZZos[24] + 8*armZZos[7];
   armZZos[24]=armZZos[16]*armZZos[24];
   armZZos[22]=8./3.*armZZos[23] + 1./3.*armZZos[22] + 8*armZZos[24];
   armZZos[23]=1./2. + armZZos[8];
   armZZos[23]=mmZ*armZZos[23];
   armZZos[23]=4*armZZos[9] + armZZos[23];
   armZZos[21]=1./3.*armZZos[23] + armZZos[21];
   armZZos[21]=armZZos[16]*armZZos[21];
   armZZos[24]= - 7./2.*mmZ + 10*armZZos[7];
   armZZos[24]=armZZos[16]*armZZos[24];
   armZZos[25]=armZZos[16]*mmZ;
   armZZos[25]=1 + armZZos[25];
   armZZos[25]=armZZos[6]*armZZos[25];
   armZZos[26]=25*armZZos[25] - 311./2. + 25*armZZos[24];
   armZZos[26]=armZZos[6]*armZZos[26];
   armZZos[27]= - 5*armZZos[12] - 17*armZZos[11];
   armZZos[27]=mmZ*armZZos[27];
   armZZos[27]=4./3.*armZZos[27] + 25*armZZos[8] - 431./3. - 10*
   armZZos[10];
   armZZos[21]=1./3.*armZZos[26] + 1./3.*armZZos[27] + 25*armZZos[21];
   armZZos[21]=armZZos[19]*armZZos[21];
   armZZos[21]=4*armZZos[22] + armZZos[21];
   armZZos[22]=armZZos[25] - 15./2. + armZZos[24];
   armZZos[22]=armZZos[6]*armZZos[22];
   armZZos[24]= - armZZos[12] - armZZos[11];
   armZZos[24]=mmZ*armZZos[24];
   armZZos[23]=armZZos[23] - 12*armZZos[7];
   armZZos[23]=armZZos[16]*armZZos[23];
   armZZos[22]=armZZos[22] + armZZos[23] + 4./3.*armZZos[24] + 
   armZZos[8] - 37./3. - 2*armZZos[10];
   armZZos[22]=armZZos[20]*armZZos[22];
   armZZos[23]=68./27.*armZZos[13] + armZZos[11];
   armZZos[23]=armZZos[19]*armZZos[23];
   armZZos[24]=4./3.*armZZos[13] + armZZos[11];
   armZZos[24]=armZZos[20]*armZZos[24];
   armZZos[23]=armZZos[24] - 128./27.*armZZos[13] + armZZos[23];
   armZZos[23]=mmt*armZZos[23];
   armZZos[17]=2*armZZos[17] + 4*armZZos[23] + 1./3.*armZZos[21] + 
   armZZos[22];
   armZZos[17]=armZZos[1]*armZZos[17];
   armZZos[18]=31./3. + armZZos[18];
   armZZos[18]=armZZos[15]*armZZos[18];
   armZZos[21]=mmZ*armZZos[15]*armZZos[12];
   armZZos[18]=armZZos[18] + 8./3.*armZZos[21];
   armZZos[21]= - 31./3. - 4*armZZos[10];
   armZZos[21]=armZZos[15]*armZZos[21];
   armZZos[22]= - mmZ*armZZos[15]*armZZos[12];
   armZZos[21]=armZZos[21] + 8./3.*armZZos[22];
   armZZos[19]=armZZos[19]*armZZos[21];
   armZZos[18]=20*armZZos[18] + 11*armZZos[19];
   armZZos[19]=armZZos[20]*armZZos[21];

      mZZosret = armZZos[17] + 1./9.*armZZos[18] + armZZos[19];
      return mZZosret;
}
