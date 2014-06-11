#include <bb.hpp>
std::complex<long double> bb::mygl20(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mybgl[25];

    mybgl[1]=pow(SW,-1);
    mybgl[2]=pow(MMH,-1);
    mybgl[3]=pow(MMW,-1);
    mybgl[4]=pow(MMt,-1);
    mybgl[5]=Tsil::I2(MMH,MMH,MMH,mu2);
    mybgl[6]=Tsil::I2(MMH,MMt,MMt,mu2);
    mybgl[7]=Tsil::I2(0,MMH,MMt,mu2);
    mybgl[8]=Tsil::B(MMH,MMt,MMt,mu2);
    mybgl[9]=log(pow(mu2,-1)*MMt);
    mybgl[10]=Tsil::B(MMt,MMt,MMH,mu2);
    mybgl[11]=Tsil::B(0,0,MMH,mu2);
    mybgl[12]=Tsil::B(0,0,MMt,mu2);
    mybgl[13]=log(pow(mu2,-1)*MMH);
    mybgl[14]=log(MMH);
    mybgl[15]=log(MMt);
   mybgl[16]=1./2.*mybgl[13];
   mybgl[17]=21*mybgl[15] - 43 + 15./2.*mybgl[14];
   mybgl[17]=mybgl[17]*mybgl[16];
   mybgl[18]=3*mybgl[9];
   mybgl[19]=mybgl[18] + 5./2.;
   mybgl[20]=mybgl[19]*mybgl[8];
   mybgl[20]=mybgl[20] + 3*mybgl[10];
   mybgl[21]=mybgl[15]*mybgl[9];
   mybgl[22]=pow(Pi,2);
   mybgl[23]=215./2. - 69*mybgl[9];
   mybgl[17]=mybgl[17] + 13./12.*mybgl[22] + 7./4.*mybgl[21] + 1./4.*
   mybgl[23] + mybgl[20];
   mybgl[17]=MMt*mybgl[17];
   mybgl[17]= - 15./2.*mybgl[5] + 5./2.*mybgl[6] + 1./2.*mybgl[7] + 
   mybgl[17];
   mybgl[23]= - mybgl[15] + 235./4. - 19*mybgl[14];
   mybgl[16]=mybgl[23]*mybgl[16];
   mybgl[16]=3./8.*mybgl[11] + 243./4.*S2 + mybgl[16] - 23./48.*
   mybgl[22] - 1./4.*mybgl[21] - 1169./32. + mybgl[9];
   mybgl[16]=MMH*mybgl[16];
   mybgl[16]=1./2.*mybgl[17] + mybgl[16];
   mybgl[17]=pow(mybgl[1],4)*pow(mybgl[3],2);
   mybgl[16]=MMH*mybgl[17]*mybgl[16];
   mybgl[23]= - 7 - 1./2.*mybgl[22];
   mybgl[24]=3 - mybgl[14];
   mybgl[24]=mybgl[13]*mybgl[24];
   mybgl[23]=1./2.*mybgl[23] + mybgl[24];
   mybgl[23]=MMH*mybgl[23];
   mybgl[24]= - 3*mybgl[7] + mybgl[6];
   mybgl[23]=1./2.*mybgl[24] + mybgl[23];
   mybgl[23]=mybgl[4]*pow(MMH,2)*mybgl[23]*mybgl[17];
   mybgl[16]=mybgl[16] + mybgl[23];
   mybgl[23]= - 233./4. + 93*mybgl[9];
   mybgl[20]= - 85./96.*mybgl[22] + 9./16.*mybgl[21] + 1./16.*mybgl[23]
    - mybgl[20];
   mybgl[20]=MMt*mybgl[20];
   mybgl[20]= - 17./4.*mybgl[6] + 49./8.*mybgl[7] + mybgl[20];
   mybgl[20]=MMt*mybgl[20];
   mybgl[18]= - 1./4.*mybgl[22] - mybgl[21] - 7./2. + mybgl[18];
   mybgl[18]=MMt*mybgl[18];
   mybgl[18]= - mybgl[7] + mybgl[18];
   mybgl[18]=mybgl[2]*mybgl[18];
   mybgl[21]= - 1./2.*mybgl[15] + 1 - 1./4.*mybgl[14];
   mybgl[21]=mybgl[13]*mybgl[21];
   mybgl[19]=mybgl[12]*mybgl[19];
   mybgl[18]= - 1./4.*mybgl[19] + 11./2.*mybgl[18] + 11*mybgl[21];
   mybgl[18]=mybgl[18]*pow(MMt,2);
   mybgl[18]=mybgl[20] + mybgl[18];
   mybgl[17]=mybgl[18]*mybgl[17];
   mybgl[16]=mybgl[17] + 1./2.*mybgl[16];

      return 1./8.*mybgl[16];
}
