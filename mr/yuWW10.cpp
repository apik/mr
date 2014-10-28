#include <WW.hpp>
std::complex<long double>
WW<OS>::my10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuWW[24], yuWWret;

    aryuWW[1]=double(nL + nH);
    aryuWW[2]=pow(SW,-1);
    aryuWW[3]=std::real(Tsil::B(0,0,MMW,mu2));
    aryuWW[4]=double(nH);
    aryuWW[5]=pow(CW,-1);
    aryuWW[6]=pow(MMZ,-1);
    aryuWW[7]=Tsil::B(0,MMt,MMW,mu2);
    aryuWW[8]=Tsil::A(MMt,mu2);
    aryuWW[9]=pow(MMH,-1);
    aryuWW[10]=double(nL);
    aryuWW[11]=double(boson);
    aryuWW[12]=Tsil::B(MMW,MMH,MMW,mu2);
    aryuWW[13]=Tsil::B(MMW,MMZ,MMW,mu2);
    aryuWW[14]=Tsil::A(MMH,mu2);
    aryuWW[15]=Tsil::A(MMZ,mu2);
    aryuWW[16]=Tsil::A(MMW,mu2);
   aryuWW[17]=MMH*aryuWW[12];
   aryuWW[17]=aryuWW[17] + aryuWW[14] - aryuWW[16];
   aryuWW[17]=MMH*aryuWW[17];
   aryuWW[18]=pow(aryuWW[5],2);
   aryuWW[19]=aryuWW[18]*aryuWW[17];
   aryuWW[19]=aryuWW[17] + aryuWW[19];
   aryuWW[19]=aryuWW[18]*aryuWW[19];
   aryuWW[20]=pow(aryuWW[2],2);
   aryuWW[17]=aryuWW[20]*aryuWW[17];
   aryuWW[17]=aryuWW[19] + aryuWW[17];
   aryuWW[17]=aryuWW[11]*aryuWW[17];
   aryuWW[19]= - MMt*aryuWW[7];
   aryuWW[19]= - aryuWW[8] + aryuWW[19];
   aryuWW[19]=aryuWW[4]*MMt*aryuWW[19];
   aryuWW[21]=aryuWW[18]*aryuWW[19];
   aryuWW[21]=aryuWW[19] + aryuWW[21];
   aryuWW[21]=aryuWW[18]*aryuWW[21];
   aryuWW[19]=aryuWW[20]*aryuWW[19];
   aryuWW[17]=1./6.*aryuWW[17] + aryuWW[21] + aryuWW[19];
   aryuWW[17]=aryuWW[6]*aryuWW[17];
   aryuWW[19]=aryuWW[15] - aryuWW[16];
   aryuWW[19]=aryuWW[18]*aryuWW[19];
   aryuWW[21]=3./2.*aryuWW[15];
   aryuWW[22]=35./6.*aryuWW[16] + aryuWW[21] + aryuWW[14];
   aryuWW[23]= - 1./2. - aryuWW[12];
   aryuWW[23]=1./3.*MMH*aryuWW[23];
   aryuWW[19]=1./12.*aryuWW[19] + 1./2.*aryuWW[22] + aryuWW[23];
   aryuWW[19]=aryuWW[18]*aryuWW[19];
   aryuWW[22]= - 13./6.*aryuWW[16] - 5./2.*aryuWW[15] + aryuWW[14];
   aryuWW[22]=1./2.*aryuWW[22] + aryuWW[23];
   aryuWW[22]=aryuWW[20]*aryuWW[22];
   aryuWW[19]=aryuWW[19] + aryuWW[22];
   aryuWW[19]=aryuWW[11]*aryuWW[19];
   aryuWW[22]= - aryuWW[9]*aryuWW[8];
   aryuWW[22]=6*aryuWW[22] + 1 - 1./2.*aryuWW[7];
   aryuWW[22]=MMt*aryuWW[22];
   aryuWW[22]=aryuWW[8] + aryuWW[22];
   aryuWW[22]=aryuWW[4]*aryuWW[22];
   aryuWW[23]=aryuWW[18]*aryuWW[22];
   aryuWW[22]=aryuWW[20]*aryuWW[22];
   aryuWW[17]=1./2.*aryuWW[17] + aryuWW[19] + aryuWW[23] + aryuWW[22];
   aryuWW[17]=aryuWW[6]*aryuWW[17];
   aryuWW[19]= - 1 + 17./2.*aryuWW[13];
   aryuWW[21]=MMZ + aryuWW[21];
   aryuWW[21]=aryuWW[9]*aryuWW[21];
   aryuWW[22]=aryuWW[18]*aryuWW[13];
   aryuWW[19]=1./12.*aryuWW[22] + 1./6.*aryuWW[19] + aryuWW[21];
   aryuWW[18]=aryuWW[18]*aryuWW[19];
   aryuWW[19]= - 1 + aryuWW[13];
   aryuWW[21]= - aryuWW[9]*MMZ;
   aryuWW[19]=2*aryuWW[19] + aryuWW[21];
   aryuWW[21]= - 59./9. - 33./2.*aryuWW[13];
   aryuWW[22]=aryuWW[16] + MMZ + 1./2.*aryuWW[15];
   aryuWW[22]=aryuWW[9]*aryuWW[22];
   aryuWW[21]=3*aryuWW[22] + 1./2.*aryuWW[21] + aryuWW[12];
   aryuWW[21]=aryuWW[20]*aryuWW[21];
   aryuWW[18]=aryuWW[21] + 2*aryuWW[19] + aryuWW[18];
   aryuWW[18]=aryuWW[11]*aryuWW[18];
   aryuWW[19]=aryuWW[3]*aryuWW[10];
   aryuWW[21]= - 1./3. + aryuWW[3];
   aryuWW[21]=aryuWW[1]*aryuWW[21];
   aryuWW[22]= - 1./3. + aryuWW[7];
   aryuWW[22]=aryuWW[4]*aryuWW[22];
   aryuWW[19]=aryuWW[22] + 1./3.*aryuWW[21] - 1./3.*aryuWW[10] + 
   aryuWW[19];
   aryuWW[19]=aryuWW[20]*aryuWW[19];

      yuWWret = aryuWW[17] + aryuWW[18] + aryuWW[19];
      return yuWWret;
}
