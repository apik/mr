#include <WW.hpp>
std::complex<long double>
WW<MS>::m10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armWWos[24], mWWosret;

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
   armWWos[17]= - mmH*armWWos[13];
   armWWos[17]=armWWos[17] - armWWos[16] + armWWos[14];
   armWWos[17]=mmH*armWWos[17];
   armWWos[18]=pow(armWWos[6],2);
   armWWos[19]=armWWos[18]*armWWos[17];
   armWWos[19]=armWWos[17] + armWWos[19];
   armWWos[19]=armWWos[18]*armWWos[19];
   armWWos[20]=pow(armWWos[2],2);
   armWWos[17]=armWWos[20]*armWWos[17];
   armWWos[17]=armWWos[19] + armWWos[17];
   armWWos[17]=armWWos[11]*armWWos[17];
   armWWos[19]=mmt*armWWos[7];
   armWWos[19]=armWWos[8] + armWWos[19];
   armWWos[19]=armWWos[4]*mmt*armWWos[19];
   armWWos[21]=armWWos[18]*armWWos[19];
   armWWos[21]=armWWos[19] + armWWos[21];
   armWWos[21]=armWWos[18]*armWWos[21];
   armWWos[19]=armWWos[20]*armWWos[19];
   armWWos[17]=1./6.*armWWos[17] + armWWos[21] + armWWos[19];
   armWWos[17]=armWWos[5]*armWWos[17];
   armWWos[19]= - armWWos[15] + armWWos[14];
   armWWos[19]=armWWos[18]*armWWos[19];
   armWWos[21]= - 3./2.*armWWos[15];
   armWWos[22]= - 35./6.*armWWos[14] - armWWos[16] + armWWos[21];
   armWWos[23]=1./2. + armWWos[13];
   armWWos[23]=1./3.*mmH*armWWos[23];
   armWWos[19]=1./12.*armWWos[19] + 1./2.*armWWos[22] + armWWos[23];
   armWWos[19]=armWWos[18]*armWWos[19];
   armWWos[22]=13./6.*armWWos[14] - armWWos[16] + 5./2.*armWWos[15];
   armWWos[22]=1./2.*armWWos[22] + armWWos[23];
   armWWos[22]=armWWos[20]*armWWos[22];
   armWWos[19]=armWWos[19] + armWWos[22];
   armWWos[19]=armWWos[11]*armWWos[19];
   armWWos[22]=armWWos[9]*armWWos[8];
   armWWos[22]=6*armWWos[22] - 1 + 1./2.*armWWos[7];
   armWWos[22]=mmt*armWWos[22];
   armWWos[22]= - armWWos[8] + armWWos[22];
   armWWos[22]=armWWos[4]*armWWos[22];
   armWWos[23]=armWWos[18]*armWWos[22];
   armWWos[22]=armWWos[20]*armWWos[22];
   armWWos[17]=1./2.*armWWos[17] + armWWos[19] + armWWos[23] + 
   armWWos[22];
   armWWos[17]=armWWos[5]*armWWos[17];
   armWWos[19]=1 - 17./2.*armWWos[12];
   armWWos[21]= - mmZ + armWWos[21];
   armWWos[21]=armWWos[9]*armWWos[21];
   armWWos[22]= - armWWos[18]*armWWos[12];
   armWWos[19]=1./12.*armWWos[22] + 1./6.*armWWos[19] + armWWos[21];
   armWWos[18]=armWWos[18]*armWWos[19];
   armWWos[19]=1 - armWWos[12];
   armWWos[21]=armWWos[9]*mmZ;
   armWWos[19]=2*armWWos[19] + armWWos[21];
   armWWos[21]=59./9. + 33./2.*armWWos[12];
   armWWos[22]= - armWWos[14] - mmZ - 1./2.*armWWos[15];
   armWWos[22]=armWWos[9]*armWWos[22];
   armWWos[21]=3*armWWos[22] + 1./2.*armWWos[21] - armWWos[13];
   armWWos[21]=armWWos[20]*armWWos[21];
   armWWos[18]=armWWos[21] + 2*armWWos[19] + armWWos[18];
   armWWos[18]=armWWos[11]*armWWos[18];
   armWWos[19]= - armWWos[3]*armWWos[10];
   armWWos[21]=1./3. - armWWos[3];
   armWWos[21]=armWWos[1]*armWWos[21];
   armWWos[22]=1./3. - armWWos[7];
   armWWos[22]=armWWos[4]*armWWos[22];
   armWWos[19]=armWWos[22] + 1./3.*armWWos[21] + 1./3.*armWWos[10] + 
   armWWos[19];
   armWWos[19]=armWWos[20]*armWWos[19];

      mWWosret = armWWos[17] + armWWos[18] + armWWos[19];
      return mWWosret;
}
