#include <tt.hpp>
std::complex<long double> tt::mygl11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mytgl[28];

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
    mytgl[12]=prot0ttHt->M(0);
    mytgl[13]=prottH0H->Vxzuv(0);
    mytgl[14]=prot0ttHt->Tuxv(0);
    mytgl[15]=protWt000->Tyzv(0);
    mytgl[16]=protHtt0t->M(0);
    mytgl[17]=prot00t00->M(0);
    mytgl[18]=prot000t0->M(0);
   mytgl[19]=2 - mytgl[11];
   mytgl[19]=mytgl[7]*mytgl[19];
   mytgl[19]=mytgl[19] + mytgl[14];
   mytgl[20]=mytgl[10]*mytgl[7];
   mytgl[21]=pow(Pi,2);
   mytgl[22]=98./3. - 45./4.*mytgl[11];
   mytgl[22]=mytgl[6]*mytgl[22];
   mytgl[23]=19 + 7*mytgl[6];
   mytgl[23]=mytgl[5]*mytgl[23];
   mytgl[24]=5./12.*mytgl[8] - 1./2. + 4./3.*mytgl[6];
   mytgl[24]=mytgl[8]*mytgl[24];
   mytgl[25]=mytgl[12] + mytgl[16];
   mytgl[26]= - 8./3.*mytgl[13] + mytgl[25];
   mytgl[26]=MMH*mytgl[26];
   mytgl[27]= - mytgl[18] - mytgl[17] - 8*mytgl[25];
   mytgl[27]=MMt*mytgl[27];
   mytgl[19]=1./3.*mytgl[27] + 2*mytgl[26] + mytgl[24] + 2./3.*
   mytgl[23] - 13./6.*mytgl[20] - 8./3.*mytgl[9] + mytgl[22] - 19./18.*
   mytgl[21] - 23 + 4*mytgl[15] + 13./3.*mytgl[19];
   mytgl[19]=MMt*mytgl[19];
   mytgl[22]=mytgl[14] + mytgl[5];
   mytgl[23]=1./2.*mytgl[7];
   mytgl[22]=1 - mytgl[23] + 1./2.*mytgl[22];
   mytgl[22]=mytgl[3]*mytgl[22];
   mytgl[22]=4*mytgl[13] - mytgl[25] + mytgl[22];
   mytgl[22]=MMH*mytgl[22];
   mytgl[24]= - 7./2.*mytgl[6] - 23./2. - 2*mytgl[7];
   mytgl[24]=mytgl[5]*mytgl[24];
   mytgl[22]=mytgl[24] + mytgl[22];
   mytgl[24]=67./3. - 3*mytgl[11];
   mytgl[23]=mytgl[24]*mytgl[23];
   mytgl[24]=mytgl[9] - mytgl[15];
   mytgl[25]= - 9./2. + 11./3.*mytgl[11];
   mytgl[25]=mytgl[6]*mytgl[25];
   mytgl[26]=mytgl[3]*mytgl[4];
   mytgl[20]= - 2*mytgl[14] + 1./6.*mytgl[26] - 17./6.*mytgl[20] + 1./2.
   *mytgl[25] + mytgl[23] - 1./9.*mytgl[21] - 23./2. + 1./3.*mytgl[22]
    + 4./3.*mytgl[24];
   mytgl[20]=MMH*mytgl[20];
   mytgl[19]=mytgl[19] - 11./3.*mytgl[4] + mytgl[20];

      return mytgl[19]*mytgl[2]*pow(mytgl[1],2);
}
