#include <ZZ.hpp>
std::complex<long double> zz::m11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mZZ2MS[31], mZZ2MSret;

    mZZ2MS[1]=double(nH);
    mZZ2MS[2]=pow(c,-1);
    mZZ2MS[3]=pow(mmH,-1);
    mZZ2MS[4]=pow(mmZ,-1);
    mZZ2MS[5]=pow(s,-1);
    mZZ2MS[6]=Tsil::B(mmt,mmt,mmZ,mu2);
    mZZ2MS[7]=Tsil::A(mmt,mu2);
    mZZ2MS[8]=Tsil::Beps(mmt,mmt,mmZ,mu2);
    mZZ2MS[9]=Tsil::Aeps(mmt,mu2);
    mZZ2MS[10]=std::real(Tsil::B(0,0,mmZ,mu2));
    mZZ2MS[11]=prottttt0->M(0);
    mZZ2MS[12]=prot00000->M(0);
    mZZ2MS[13]=prottttt0->Vzxyv(0);
    mZZ2MS[14]=prottttt0->Suxv(0);
    mZZ2MS[15]=double(nL);
    mZZ2MS[16]=1/(4*mmt - mmZ);
   mZZ2MS[17]=pow(mZZ2MS[5],2);
   mZZ2MS[18]=pow(mZZ2MS[2],2);
   mZZ2MS[19]=mZZ2MS[17] + mZZ2MS[18];
   mZZ2MS[20]=2*mZZ2MS[4];
   mZZ2MS[21]=mZZ2MS[19]*mZZ2MS[20];
   mZZ2MS[22]= - mZZ2MS[9]*mZZ2MS[21];
   mZZ2MS[22]=mZZ2MS[22] + 85./9.*mZZ2MS[18] + 128./9. + 13*mZZ2MS[17];
   mZZ2MS[22]=mZZ2MS[22]*mZZ2MS[20];
   mZZ2MS[23]= - 7./9.*mZZ2MS[18] + mZZ2MS[17] + 64./9.;
   mZZ2MS[23]=mZZ2MS[23]*mZZ2MS[4];
   mZZ2MS[24]=mZZ2MS[6]*mZZ2MS[23];
   mZZ2MS[25]= - 13./3.*mZZ2MS[18] + 64./3. + mZZ2MS[17];
   mZZ2MS[25]=mZZ2MS[4]*mZZ2MS[25];
   mZZ2MS[25]=5*mZZ2MS[25] + mZZ2MS[24];
   mZZ2MS[25]=mZZ2MS[6]*mZZ2MS[25];
   mZZ2MS[26]=17./9.*mZZ2MS[18] + mZZ2MS[17] - 32./9.;
   mZZ2MS[27]=mZZ2MS[13]*mZZ2MS[26];
   mZZ2MS[22]=8*mZZ2MS[27] + mZZ2MS[22] + mZZ2MS[25];
   mZZ2MS[25]=pow(mZZ2MS[4],2);
   mZZ2MS[27]=mZZ2MS[19]*mZZ2MS[25];
   mZZ2MS[28]=mZZ2MS[6] + 2;
   mZZ2MS[28]=mZZ2MS[28]*mZZ2MS[27];
   mZZ2MS[29]= - mZZ2MS[6]*mZZ2MS[28];
   mZZ2MS[27]=mZZ2MS[27] + mZZ2MS[29];
   mZZ2MS[29]=mZZ2MS[4]*mZZ2MS[19]*mZZ2MS[3];
   mZZ2MS[30]= - 2./3.*mZZ2MS[11] - 4./3.*mZZ2MS[13];
   mZZ2MS[30]=mZZ2MS[30]*mZZ2MS[23];
   mZZ2MS[27]=1./3.*mZZ2MS[27] - 12*mZZ2MS[29] + mZZ2MS[30];
   mZZ2MS[30]=2*mmt;
   mZZ2MS[27]=mZZ2MS[27]*mZZ2MS[30];
   mZZ2MS[21]= - mZZ2MS[3]*mZZ2MS[21];
   mZZ2MS[21]= - 1./3.*mZZ2MS[28] + mZZ2MS[21];
   mZZ2MS[28]=4*mZZ2MS[7];
   mZZ2MS[21]=mZZ2MS[21]*mZZ2MS[28];
   mZZ2MS[25]=mZZ2MS[14]*mZZ2MS[25];
   mZZ2MS[25]=2./3.*mZZ2MS[25] + 2*mZZ2MS[11];
   mZZ2MS[19]=mZZ2MS[19]*mZZ2MS[25];
   mZZ2MS[23]=mZZ2MS[8]*mZZ2MS[23];
   mZZ2MS[19]=mZZ2MS[27] - 2./3.*mZZ2MS[23] + mZZ2MS[21] + 1./3.*
   mZZ2MS[22] + mZZ2MS[19];
   mZZ2MS[19]=mZZ2MS[30]*mZZ2MS[19];
   mZZ2MS[21]= - 55./9.*mZZ2MS[18] + 256./9. + mZZ2MS[17];
   mZZ2MS[20]=mZZ2MS[21]*mZZ2MS[20];
   mZZ2MS[20]=mZZ2MS[20] - 5*mZZ2MS[24];
   mZZ2MS[21]=3*mZZ2MS[17] - 64./3. + 25./3.*mZZ2MS[18];
   mZZ2MS[22]=25./9.*mZZ2MS[18] + mZZ2MS[17] - 64./9.;
   mZZ2MS[23]=5*mZZ2MS[6];
   mZZ2MS[23]=mZZ2MS[22]*mZZ2MS[23];
   mZZ2MS[23]= - 2*mZZ2MS[21] + mZZ2MS[23];
   mZZ2MS[23]=mZZ2MS[16]*mZZ2MS[23];
   mZZ2MS[21]=mZZ2MS[16]*mZZ2MS[4]*mZZ2MS[21];
   mZZ2MS[21]= - 6*mZZ2MS[29] + mZZ2MS[21];
   mZZ2MS[21]=mZZ2MS[21]*mZZ2MS[28];
   mZZ2MS[20]=mZZ2MS[21] + 1./3.*mZZ2MS[20] + mZZ2MS[23];
   mZZ2MS[20]=mZZ2MS[7]*mZZ2MS[20];
   mZZ2MS[21]=mZZ2MS[16]*mZZ2MS[9];
   mZZ2MS[23]=mmZ*mZZ2MS[16];
   mZZ2MS[23]=mZZ2MS[23] + 1;
   mZZ2MS[23]=mZZ2MS[8]*mZZ2MS[23];
   mZZ2MS[21]=mZZ2MS[23] + 4*mZZ2MS[21];
   mZZ2MS[21]=mZZ2MS[22]*mZZ2MS[21];
   mZZ2MS[23]= - 95./9.*mZZ2MS[18] + 128./9. - 7*mZZ2MS[17];
   mZZ2MS[24]=4*mZZ2MS[4];
   mZZ2MS[23]=mZZ2MS[24]*mZZ2MS[9]*mZZ2MS[23];
   mZZ2MS[23]=mZZ2MS[23] - 431./9.*mZZ2MS[18] + 764./9. - 37*mZZ2MS[17]
   ;
   mZZ2MS[22]=mZZ2MS[22]*mZZ2MS[6];
   mZZ2MS[25]=mZZ2MS[22] - 311./18.*mZZ2MS[18] + 352./9. - 15./2.*
   mZZ2MS[17];
   mZZ2MS[25]=mZZ2MS[6]*mZZ2MS[25];
   mZZ2MS[27]=1./2.*mZZ2MS[17] - 32./9. + 25./18.*mZZ2MS[18];
   mZZ2MS[22]= - 7*mZZ2MS[27] + mZZ2MS[22];
   mZZ2MS[22]=mZZ2MS[6]*mZZ2MS[22];
   mZZ2MS[22]=mZZ2MS[22] + mZZ2MS[27];
   mZZ2MS[22]=mZZ2MS[16]*mZZ2MS[22];
   mZZ2MS[27]=mZZ2MS[11]*mZZ2MS[26];
   mZZ2MS[22]= - 4./3.*mZZ2MS[27] + mZZ2MS[22];
   mZZ2MS[22]=mmZ*mZZ2MS[22];
   mZZ2MS[27]=mZZ2MS[17] - 8./9. + 5./9.*mZZ2MS[18];
   mZZ2MS[28]=mZZ2MS[12]*mmZ;
   mZZ2MS[29]= - 2*mZZ2MS[10] - 4./3.*mZZ2MS[28];
   mZZ2MS[27]=mZZ2MS[27]*mZZ2MS[29];
   mZZ2MS[24]=mZZ2MS[14]*mZZ2MS[26]*mZZ2MS[24];
   mZZ2MS[19]=mZZ2MS[19] + mZZ2MS[24] + 2*mZZ2MS[20] + mZZ2MS[22] + 1./
   3.*mZZ2MS[23] + mZZ2MS[25] + mZZ2MS[27] + mZZ2MS[21];
   mZZ2MS[19]=mZZ2MS[1]*mZZ2MS[19];
   mZZ2MS[17]=11./9.*mZZ2MS[18] + mZZ2MS[17] - 20./9.;
   mZZ2MS[18]= - 31 - 8*mZZ2MS[28];
   mZZ2MS[18]=1./3.*mZZ2MS[18] - 4*mZZ2MS[10];
   mZZ2MS[17]=mZZ2MS[18]*mZZ2MS[17]*mZZ2MS[15];

      mZZ2MSret = mZZ2MS[17] + mZZ2MS[19];
      return mZZ2MSret;
}
