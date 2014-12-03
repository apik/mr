#include <WW.hpp>
long double WW<MS>::x10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armWWos[26], mWWosret;

    armWWos[1]=double(nL + nH);
    armWWos[2]=pow(s,-1);
    armWWos[3]=std::real(Tsil::B(0,0,mmW,mu2));
    armWWos[4]=double(nH);
    armWWos[5]=pow(mmZ,-1);
    armWWos[6]=pow(c,-1);
    armWWos[7]=Tsil::B(0,mmt,mmW,mu2);
    armWWos[8]=Tsil::A(mmt,mu2);
    armWWos[9]=pow(mmH,-1);
    armWWos[10]=double(nL);
    armWWos[11]=double(boson);
    armWWos[12]=Tsil::B(mmW,mmZ,mmW,mu2);
    armWWos[13]=Tsil::B(mmW,mmH,mmW,mu2);
    armWWos[14]=Tsil::A(mmW,mu2);
    armWWos[15]=Tsil::A(mmZ,mu2);
    armWWos[16]=Tsil::A(mmH,mu2);
   armWWos[17]=armWWos[9]*armWWos[8];
   armWWos[17]=6*armWWos[17] - 1 + 1./2.*armWWos[7];
   armWWos[18]=mmt*armWWos[4];
   armWWos[17]=armWWos[18]*armWWos[17];
   armWWos[19]= - armWWos[4]*armWWos[8];
   armWWos[17]=armWWos[17] + armWWos[19];
   armWWos[19]=pow(armWWos[2],2);
   armWWos[20]=pow(armWWos[6],2);
   armWWos[21]=armWWos[19] + armWWos[20];
   armWWos[17]=armWWos[21]*armWWos[17];
   armWWos[21]=armWWos[13] + 1./2.;
   armWWos[22]=1./3.*mmH;
   armWWos[21]=armWWos[21]*armWWos[22];
   armWWos[22]=1./12.*armWWos[20];
   armWWos[23]= - armWWos[15] + armWWos[14];
   armWWos[23]=armWWos[23]*armWWos[22];
   armWWos[24]=3./2.*armWWos[15];
   armWWos[25]= - armWWos[16] - armWWos[24];
   armWWos[23]=armWWos[23] - 35./12.*armWWos[14] + 1./2.*armWWos[25] + 
   armWWos[21];
   armWWos[23]=armWWos[23]*armWWos[20];
   armWWos[25]= - armWWos[16] + 5./2.*armWWos[15];
   armWWos[21]=13./12.*armWWos[14] + 1./2.*armWWos[25] + armWWos[21];
   armWWos[21]=armWWos[21]*armWWos[19];
   armWWos[21]=armWWos[23] + armWWos[21];
   armWWos[21]=armWWos[11]*armWWos[21];
   armWWos[23]=armWWos[20] + 1;
   armWWos[23]=armWWos[23]*armWWos[20];
   armWWos[23]=armWWos[23] + armWWos[19];
   armWWos[25]=armWWos[4]*armWWos[8]*armWWos[23];
   armWWos[18]=armWWos[18]*armWWos[7]*armWWos[23];
   armWWos[18]=armWWos[25] + armWWos[18];
   armWWos[18]=mmt*armWWos[18];
   armWWos[25]=mmH*armWWos[13];
   armWWos[25]= - armWWos[14] + armWWos[25] + armWWos[16];
   armWWos[23]= - armWWos[11]*armWWos[23]*mmH*armWWos[25];
   armWWos[18]=1./6.*armWWos[23] + armWWos[18];
   armWWos[18]=armWWos[5]*armWWos[18];
   armWWos[17]=1./2.*armWWos[18] + armWWos[21] + armWWos[17];
   armWWos[17]=armWWos[5]*armWWos[17];
   armWWos[18]= - mmZ - armWWos[24];
   armWWos[18]=armWWos[9]*armWWos[18];
   armWWos[18]=1./6. + armWWos[18];
   armWWos[18]=armWWos[18]*armWWos[20];
   armWWos[21]= - armWWos[14] - mmZ - 1./2.*armWWos[15];
   armWWos[23]=3*armWWos[9];
   armWWos[21]=armWWos[23]*armWWos[21];
   armWWos[21]=59./18. - armWWos[13] + armWWos[21];
   armWWos[21]=armWWos[21]*armWWos[19];
   armWWos[20]= - 17 - armWWos[20];
   armWWos[20]=armWWos[20]*armWWos[22];
   armWWos[20]=33./4.*armWWos[19] - 4 + armWWos[20];
   armWWos[20]=armWWos[12]*armWWos[20];
   armWWos[22]=armWWos[9]*mmZ;
   armWWos[22]=2 + armWWos[22];
   armWWos[18]=armWWos[20] + armWWos[21] + 2*armWWos[22] + armWWos[18];
   armWWos[18]=armWWos[11]*armWWos[18];
   armWWos[20]= - armWWos[3] + 1./3.;
   armWWos[21]=armWWos[10] + 1./3.*armWWos[1];
   armWWos[20]=armWWos[21]*armWWos[20];
   armWWos[21]=1./3. - armWWos[7];
   armWWos[21]=armWWos[4]*armWWos[21];
   armWWos[20]=armWWos[21] + armWWos[20];
   armWWos[19]=armWWos[19]*armWWos[20];

      mWWosret = armWWos[17] + armWWos[18] + armWWos[19];
      return mWWosret.real();
}
