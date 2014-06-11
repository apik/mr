#include <tt.hpp>
std::complex<long double> tt::mytgl11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mytgl[27];

    mytgl[1]=pow(SW,-1);
    mytgl[2]=pow(MMW,-1);
    mytgl[3]=pow(MMt,-1);
    mytgl[4]=Tsil::I2(MMH,MMt,MMt,mu2);
    mytgl[5]=Tsil::B(MMH,MMt,MMt,mu2);
    mytgl[6]=log(pow(mu2,-1)*MMt);
    mytgl[7]=log(pow(mu2,-1)*MMH);
    mytgl[8]=Tsil::B(0,0,MMt,mu2);
    mytgl[9]=Tsil::Beps(MMH,MMt,MMt,mu2);
    mytgl[10]=log(MMH);
    mytgl[11]=log(MMt);
    mytgl[12]=Tfin(MMH,MMt,0);
    mytgl[13]=Tfin(MMt,0,0);
    mytgl[14]=Vfin(MMt,MMH,0,MMH);
    mytgl[15]=Mfin(MMH,MMt,MMt,0,MMt);
    mytgl[16]=Mfin(0,MMt,MMt,MMH,MMt);
    mytgl[17]=Mfin(0,0,MMt,0,0);
    mytgl[18]=Mfin(0,0,0,MMt,0);
   mytgl[19]= - mytgl[5] + mytgl[7] - mytgl[12];
   mytgl[19]=1 - 1./2.*mytgl[19];
   mytgl[19]=mytgl[3]*mytgl[19];
   mytgl[20]=mytgl[15] + mytgl[16];
   mytgl[19]=4*mytgl[14] - mytgl[20] + mytgl[19];
   mytgl[19]=MMH*mytgl[19];
   mytgl[21]=mytgl[13] - mytgl[9];
   mytgl[22]=pow(Pi,2);
   mytgl[23]= - 1 - 1./3.*mytgl[5];
   mytgl[24]=11./3.*mytgl[11] - 9./2. - 7./3.*mytgl[5];
   mytgl[24]=mytgl[6]*mytgl[24];
   mytgl[25]=mytgl[3]*mytgl[4];
   mytgl[26]=67./2. - 2*mytgl[5];
   mytgl[26]= - 17./6.*mytgl[10] + 1./3.*mytgl[26] - 3./2.*mytgl[11];
   mytgl[26]=mytgl[7]*mytgl[26];
   mytgl[19]= - 1./9.*mytgl[22] + 1./3.*mytgl[19] + mytgl[26] + 1./6.*
   mytgl[25] + 1./2.*mytgl[24] + 23./2.*mytgl[23] - 2*mytgl[12] - 4./3.
   *mytgl[21];
   mytgl[19]=MMH*mytgl[19];
   mytgl[19]= - 11./3.*mytgl[4] + mytgl[19];
   mytgl[21]=mytgl[2]*pow(mytgl[1],2);
   mytgl[19]=mytgl[19]*mytgl[21];
   mytgl[23]= - 1./2.*mytgl[10] + 2 - mytgl[11];
   mytgl[23]=mytgl[7]*mytgl[23];
   mytgl[23]=mytgl[23] + mytgl[12];
   mytgl[24]=7 + mytgl[5];
   mytgl[24]=14./3.*mytgl[24] - 45./4.*mytgl[11];
   mytgl[24]=mytgl[6]*mytgl[24];
   mytgl[25]=5./12.*mytgl[8] - 1./2. + 4./3.*mytgl[6];
   mytgl[25]=mytgl[8]*mytgl[25];
   mytgl[26]= - 8./3.*mytgl[14] + mytgl[20];
   mytgl[26]=MMH*mytgl[26];
   mytgl[20]= - mytgl[17] - 8*mytgl[20] - mytgl[18];
   mytgl[20]=MMt*mytgl[20];
   mytgl[20]=1./3.*mytgl[20] - 19./18.*mytgl[22] + 2*mytgl[26] + 
   mytgl[25] + mytgl[24] + 4*mytgl[13] - 8./3.*mytgl[9] - 23 + 38./3.*
   mytgl[5] + 13./3.*mytgl[23];
   mytgl[20]=MMt*mytgl[21]*mytgl[20];

      return mytgl[19] + mytgl[20];
}
