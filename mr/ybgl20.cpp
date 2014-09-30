#include <bb.hpp>
std::complex<long double> bb::mygl20(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mybgl[28];

    mybgl[1]=pow(SW,-1);
    mybgl[2]=pow(MMH,-1);
    mybgl[3]=pow(MMW,-1);
    mybgl[4]=Tsil::I2(MMH,MMH,MMH,mu2);
    mybgl[5]=Tsil::I2(MMH,MMt,MMt,mu2);
    mybgl[6]=pow(MMt,-1);
    mybgl[7]=Tsil::I2(0,MMH,MMt,mu2);
    mybgl[8]=Tsil::I2(0,0,MMH,mu2);
    mybgl[9]=Tsil::I2(0,0,MMt,mu2);
    mybgl[10]=Tsil::B(MMH,MMt,MMt,mu2);
    mybgl[11]=log(pow(mu2,-1)*MMt);
    mybgl[12]=Tsil::B(MMt,MMt,MMH,mu2);
    mybgl[13]=Tsil::B(0,0,MMH,mu2);
    mybgl[14]=Tsil::B(0,0,MMt,mu2);
    mybgl[15]=Tsil::Aeps(MMH,mu2);
    mybgl[16]=Tsil::Aeps(MMt,mu2);
    mybgl[17]=log(pow(mu2,-1)*MMH);
   mybgl[18]=1./2.*mybgl[11];
   mybgl[19]=7*mybgl[11] - 83./2. + 21*mybgl[17];
   mybgl[19]=mybgl[19]*mybgl[18];
   mybgl[20]=3*mybgl[11];
   mybgl[21]=mybgl[20] + 5./2.;
   mybgl[22]=mybgl[21]*mybgl[10];
   mybgl[23]=3*mybgl[12];
   mybgl[24]=pow(Pi,2);
   mybgl[25]= - 59 + 45./2.*mybgl[17];
   mybgl[25]=mybgl[17]*mybgl[25];
   mybgl[19]=mybgl[22] + mybgl[23] + mybgl[19] + 9./2.*mybgl[24] + 543./
   8. + mybgl[25];
   mybgl[22]=1./4.*MMt;
   mybgl[19]=mybgl[22]*mybgl[19];
   mybgl[25]=1 - mybgl[17];
   mybgl[25]=mybgl[25]*mybgl[18];
   mybgl[26]=mybgl[16] + mybgl[5];
   mybgl[26]= - 3./2.*mybgl[7] + mybgl[8] + 1./2.*mybgl[26];
   mybgl[26]=mybgl[6]*mybgl[26];
   mybgl[27]=25./2. - 3*mybgl[17];
   mybgl[27]=mybgl[17]*mybgl[27];
   mybgl[27]= - 1181./8. + 9*mybgl[27];
   mybgl[25]=3./8.*mybgl[13] + mybgl[25] + 1./4.*mybgl[27] - 1./3.*
   mybgl[24] + mybgl[26];
   mybgl[26]=1./2.*MMH;
   mybgl[25]=mybgl[25]*mybgl[26];
   mybgl[19]=mybgl[25] + 7./8.*mybgl[16] + 1./8.*mybgl[7] - 17./8.*
   mybgl[8] - 15./8.*mybgl[4] + 5./8.*mybgl[5] + 7*mybgl[15] + 
   mybgl[19];
   mybgl[19]=mybgl[19]*mybgl[26];
   mybgl[25]=123./8.*mybgl[11] - 73./8. - 11*mybgl[17];
   mybgl[18]=mybgl[25]*mybgl[18];
   mybgl[25]=37 - 63./4.*mybgl[17];
   mybgl[25]=mybgl[17]*mybgl[25];
   mybgl[18]= - mybgl[23] + mybgl[18] - 35./16.*mybgl[24] - 1357./64.
    + mybgl[25];
   mybgl[18]=MMt*mybgl[18];
   mybgl[23]=49./2.*mybgl[7] - 17*mybgl[5] + 75./2.*mybgl[15];
   mybgl[18]= - 31./16.*mybgl[9] + 1./2.*mybgl[18] + 1./8.*mybgl[23] + 
   11*mybgl[16];
   mybgl[18]=MMt*mybgl[18];
   mybgl[23]= - 1./8.*mybgl[14] - 1./2.*mybgl[10];
   mybgl[25]=pow(MMt,2);
   mybgl[21]=mybgl[23]*mybgl[21]*mybgl[25];
   mybgl[23]=29 - 17*mybgl[11];
   mybgl[20]=mybgl[23]*mybgl[20];
   mybgl[23]= - 159 - 29./2.*mybgl[24];
   mybgl[20]=1./2.*mybgl[23] + mybgl[20];
   mybgl[20]=mybgl[20]*mybgl[22];
   mybgl[20]=13./2.*mybgl[9] + mybgl[20] - 33*mybgl[16] - 13*mybgl[15]
    - 11./4.*mybgl[7];
   mybgl[20]=mybgl[2]*mybgl[25]*mybgl[20];
   mybgl[22]=S2*pow(MMH,2);
   mybgl[18]=mybgl[20] + 243./16.*mybgl[22] + mybgl[19] + mybgl[18] + 
   mybgl[21];

      return 1./4.*mybgl[18]*pow(mybgl[3],2)*pow(mybgl[1],4);
}
