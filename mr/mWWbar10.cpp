#include <WW.hpp>
std::complex<long double>
WW<OS>::m10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armWWbar[26], mWWbarret;

    armWWbar[1]=double(nL + nH);
    armWWbar[2]=pow(SW,-1);
    armWWbar[3]=std::real(Tsil::B(0,0,MMW,mu2));
    armWWbar[4]=double(nH);
    armWWbar[5]=pow(CW,-1);
    armWWbar[6]=pow(MMZ,-1);
    armWWbar[7]=Tsil::B(0,MMt,MMW,mu2);
    armWWbar[8]=Tsil::A(MMt,mu2);
    armWWbar[9]=pow(MMH,-1);
    armWWbar[10]=double(nL);
    armWWbar[11]=double(boson);
    armWWbar[12]=Tsil::B(MMW,MMH,MMW,mu2);
    armWWbar[13]=Tsil::B(MMW,MMZ,MMW,mu2);
    armWWbar[14]=Tsil::A(MMH,mu2);
    armWWbar[15]=Tsil::A(MMZ,mu2);
    armWWbar[16]=Tsil::A(MMW,mu2);
   armWWbar[17]=armWWbar[8]*armWWbar[9];
   armWWbar[17]=6*armWWbar[17] - 1 + 1./2.*armWWbar[7];
   armWWbar[17]=armWWbar[17]*MMt;
   armWWbar[17]=armWWbar[17] - armWWbar[8];
   armWWbar[17]=armWWbar[17]*armWWbar[4];
   armWWbar[18]=armWWbar[14]*armWWbar[11];
   armWWbar[17]=armWWbar[17] - 1./2.*armWWbar[18];
   armWWbar[18]=1./3.*armWWbar[11];
   armWWbar[19]=armWWbar[12] + 1./2.;
   armWWbar[19]=armWWbar[19]*MMH;
   armWWbar[20]=35./4.*armWWbar[16] - armWWbar[19];
   armWWbar[20]=armWWbar[20]*armWWbar[18];
   armWWbar[21]=armWWbar[15]*armWWbar[11];
   armWWbar[22]= - armWWbar[11]*armWWbar[16];
   armWWbar[22]=armWWbar[22] + armWWbar[21];
   armWWbar[23]=pow(armWWbar[5],2);
   armWWbar[24]=1./12.*armWWbar[23];
   armWWbar[22]=armWWbar[22]*armWWbar[24];
   armWWbar[20]=armWWbar[22] + 3./4.*armWWbar[21] + armWWbar[20] - 
   armWWbar[17];
   armWWbar[20]=armWWbar[20]*armWWbar[23];
   armWWbar[19]= - 13./4.*armWWbar[16] - armWWbar[19];
   armWWbar[18]=armWWbar[19]*armWWbar[18];
   armWWbar[17]= - 5./4.*armWWbar[21] + armWWbar[18] - armWWbar[17];
   armWWbar[18]=pow(armWWbar[2],2);
   armWWbar[17]=armWWbar[17]*armWWbar[18];
   armWWbar[19]=armWWbar[23] + 1;
   armWWbar[19]=armWWbar[23]*armWWbar[19];
   armWWbar[19]=armWWbar[19] + armWWbar[18];
   armWWbar[22]=MMH*armWWbar[12];
   armWWbar[22]=armWWbar[14] + armWWbar[22] - armWWbar[16];
   armWWbar[25]=armWWbar[11]*MMH;
   armWWbar[25]=1./6.*armWWbar[25];
   armWWbar[22]=armWWbar[22]*armWWbar[25];
   armWWbar[25]=MMt*armWWbar[7];
   armWWbar[25]=armWWbar[25] + armWWbar[8];
   armWWbar[25]=armWWbar[25]*armWWbar[4]*MMt;
   armWWbar[22]=armWWbar[22] - armWWbar[25];
   armWWbar[19]=armWWbar[6]*armWWbar[22]*armWWbar[19];
   armWWbar[17]=1./2.*armWWbar[19] + armWWbar[20] + armWWbar[17];
   armWWbar[17]=armWWbar[6]*armWWbar[17];
   armWWbar[19]=1./3.*armWWbar[1] + armWWbar[10];
   armWWbar[20]=armWWbar[3] - 1./3.;
   armWWbar[19]=armWWbar[20]*armWWbar[19];
   armWWbar[20]=MMZ*armWWbar[9];
   armWWbar[22]=armWWbar[9]*armWWbar[16];
   armWWbar[22]=armWWbar[20] + armWWbar[22];
   armWWbar[22]= - 33./4.*armWWbar[13] - 59./18. + armWWbar[12] + 3*
   armWWbar[22];
   armWWbar[22]=armWWbar[11]*armWWbar[22];
   armWWbar[21]=armWWbar[21]*armWWbar[9];
   armWWbar[21]=3./2.*armWWbar[21];
   armWWbar[25]= - 1./3. + armWWbar[7];
   armWWbar[25]=armWWbar[4]*armWWbar[25];
   armWWbar[19]=armWWbar[21] + armWWbar[22] + armWWbar[25] + 
   armWWbar[19];
   armWWbar[18]=armWWbar[19]*armWWbar[18];
   armWWbar[19]=armWWbar[24] + 17./12.;
   armWWbar[19]=armWWbar[13]*armWWbar[19];
   armWWbar[19]= - 1./6. + armWWbar[20] + armWWbar[19];
   armWWbar[19]=armWWbar[11]*armWWbar[19];
   armWWbar[19]=armWWbar[21] + armWWbar[19];
   armWWbar[19]=armWWbar[19]*armWWbar[23];
   armWWbar[20]=2*armWWbar[13] - 2 - armWWbar[20];
   armWWbar[20]=armWWbar[11]*armWWbar[20];

      mWWbarret = armWWbar[17] + armWWbar[18] + armWWbar[19] + 2*
      armWWbar[20];
      return mWWbarret;
}
