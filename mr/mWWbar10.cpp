#include <WW.hpp>
std::complex<long double>
WW<OS>::m10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armWWbar[24], mWWbarret;

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
   armWWbar[17]=MMH*armWWbar[12];
   armWWbar[17]=armWWbar[17] + armWWbar[14] - armWWbar[16];
   armWWbar[17]=MMH*armWWbar[17];
   armWWbar[18]=pow(armWWbar[5],2);
   armWWbar[19]=armWWbar[18]*armWWbar[17];
   armWWbar[19]=armWWbar[17] + armWWbar[19];
   armWWbar[19]=armWWbar[18]*armWWbar[19];
   armWWbar[20]=pow(armWWbar[2],2);
   armWWbar[17]=armWWbar[20]*armWWbar[17];
   armWWbar[17]=armWWbar[19] + armWWbar[17];
   armWWbar[17]=armWWbar[11]*armWWbar[17];
   armWWbar[19]= - MMt*armWWbar[7];
   armWWbar[19]= - armWWbar[8] + armWWbar[19];
   armWWbar[19]=armWWbar[4]*MMt*armWWbar[19];
   armWWbar[21]=armWWbar[18]*armWWbar[19];
   armWWbar[21]=armWWbar[19] + armWWbar[21];
   armWWbar[21]=armWWbar[18]*armWWbar[21];
   armWWbar[19]=armWWbar[20]*armWWbar[19];
   armWWbar[17]=1./6.*armWWbar[17] + armWWbar[21] + armWWbar[19];
   armWWbar[17]=armWWbar[6]*armWWbar[17];
   armWWbar[19]=armWWbar[15] - armWWbar[16];
   armWWbar[19]=armWWbar[18]*armWWbar[19];
   armWWbar[21]=3./2.*armWWbar[15];
   armWWbar[22]=35./6.*armWWbar[16] + armWWbar[21] + armWWbar[14];
   armWWbar[23]= - 1./2. - armWWbar[12];
   armWWbar[23]=1./3.*MMH*armWWbar[23];
   armWWbar[19]=1./12.*armWWbar[19] + 1./2.*armWWbar[22] + armWWbar[23]
   ;
   armWWbar[19]=armWWbar[18]*armWWbar[19];
   armWWbar[22]= - 13./6.*armWWbar[16] - 5./2.*armWWbar[15] + 
   armWWbar[14];
   armWWbar[22]=1./2.*armWWbar[22] + armWWbar[23];
   armWWbar[22]=armWWbar[20]*armWWbar[22];
   armWWbar[19]=armWWbar[19] + armWWbar[22];
   armWWbar[19]=armWWbar[11]*armWWbar[19];
   armWWbar[22]= - armWWbar[9]*armWWbar[8];
   armWWbar[22]=6*armWWbar[22] + 1 - 1./2.*armWWbar[7];
   armWWbar[22]=MMt*armWWbar[22];
   armWWbar[22]=armWWbar[8] + armWWbar[22];
   armWWbar[22]=armWWbar[4]*armWWbar[22];
   armWWbar[23]=armWWbar[18]*armWWbar[22];
   armWWbar[22]=armWWbar[20]*armWWbar[22];
   armWWbar[17]=1./2.*armWWbar[17] + armWWbar[19] + armWWbar[23] + 
   armWWbar[22];
   armWWbar[17]=armWWbar[6]*armWWbar[17];
   armWWbar[19]= - 1 + 17./2.*armWWbar[13];
   armWWbar[21]=MMZ + armWWbar[21];
   armWWbar[21]=armWWbar[9]*armWWbar[21];
   armWWbar[22]=armWWbar[18]*armWWbar[13];
   armWWbar[19]=1./12.*armWWbar[22] + 1./6.*armWWbar[19] + armWWbar[21]
   ;
   armWWbar[18]=armWWbar[18]*armWWbar[19];
   armWWbar[19]= - 1 + armWWbar[13];
   armWWbar[21]= - armWWbar[9]*MMZ;
   armWWbar[19]=2*armWWbar[19] + armWWbar[21];
   armWWbar[21]= - 59./9. - 33./2.*armWWbar[13];
   armWWbar[22]=armWWbar[16] + MMZ + 1./2.*armWWbar[15];
   armWWbar[22]=armWWbar[9]*armWWbar[22];
   armWWbar[21]=3*armWWbar[22] + 1./2.*armWWbar[21] + armWWbar[12];
   armWWbar[21]=armWWbar[20]*armWWbar[21];
   armWWbar[18]=armWWbar[21] + 2*armWWbar[19] + armWWbar[18];
   armWWbar[18]=armWWbar[11]*armWWbar[18];
   armWWbar[19]=armWWbar[3]*armWWbar[10];
   armWWbar[21]= - 1./3. + armWWbar[3];
   armWWbar[21]=armWWbar[1]*armWWbar[21];
   armWWbar[22]= - 1./3. + armWWbar[7];
   armWWbar[22]=armWWbar[4]*armWWbar[22];
   armWWbar[19]=armWWbar[22] + 1./3.*armWWbar[21] - 1./3.*armWWbar[10]
    + armWWbar[19];
   armWWbar[19]=armWWbar[20]*armWWbar[19];

      mWWbarret = armWWbar[17] + armWWbar[18] + armWWbar[19];
      return mWWbarret;
}
