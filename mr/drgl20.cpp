#include <dr.hpp>
std::complex<long double>
dr::drgl20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> ardrgl[25], drglret;

    ardrgl[1]=double(boson);
    ardrgl[2]=pow(SW,-1);
    ardrgl[3]=pow(MMH,-1);
    ardrgl[4]=pow(MMW,-1);
    ardrgl[5]=Tsil::I2(MMH,MMt,MMt,mu2);
    ardrgl[6]=Tsil::I2(0,MMH,MMt,mu2);
    ardrgl[7]=Tsil::B(MMH,MMH,MMH,mu2);
    ardrgl[8]=Tsil::A(MMH,mu2);
    ardrgl[9]=Tsil::A(MMt,mu2);
    ardrgl[10]=Tsil::B(MMH,MMt,MMt,mu2);
    ardrgl[11]=Tsil::B(MMt,MMt,MMH,mu2);
    ardrgl[12]=std::real(Tsil::B(0,0,MMH,mu2));
    ardrgl[13]=pow(MMt,-1);
    ardrgl[14]=std::real(Tsil::B(0,0,MMt,mu2));
   ardrgl[15]=5*ardrgl[8];
   ardrgl[16]=9*ardrgl[5];
   ardrgl[17]=pow(ardrgl[8],2);
   ardrgl[18]=ardrgl[17]*ardrgl[3];
   ardrgl[19]=ardrgl[16] - ardrgl[15] + 21./4.*ardrgl[18];
   ardrgl[19]=ardrgl[3]*ardrgl[19];
   ardrgl[19]=3./4.*ardrgl[14] + 59./16. + ardrgl[19];
   ardrgl[20]=3*ardrgl[8];
   ardrgl[21]= - ardrgl[3]*ardrgl[20];
   ardrgl[21]= - 5./2. + ardrgl[21];
   ardrgl[21]=ardrgl[11]*ardrgl[21];
   ardrgl[19]=5./2.*ardrgl[10] + 1./2.*ardrgl[19] + ardrgl[21];
   ardrgl[21]=ardrgl[10]*ardrgl[3];
   ardrgl[22]=77./8. - ardrgl[14];
   ardrgl[22]=ardrgl[3]*ardrgl[22];
   ardrgl[23]=ardrgl[9]*ardrgl[11]*pow(ardrgl[3],2);
   ardrgl[22]= - 12*ardrgl[23] + 1./2.*ardrgl[22] - 2*ardrgl[21];
   ardrgl[22]=MMt*ardrgl[22];
   ardrgl[23]=3*ardrgl[11] - 19./8. - ardrgl[14];
   ardrgl[23]=ardrgl[3]*ardrgl[23];
   ardrgl[21]= - 4*ardrgl[21] + ardrgl[23];
   ardrgl[21]=ardrgl[9]*ardrgl[21];
   ardrgl[19]=ardrgl[22] + 1./2.*ardrgl[19] + ardrgl[21];
   ardrgl[19]=MMt*ardrgl[19];
   ardrgl[21]=3./2.*ardrgl[12] + 9./2.*ardrgl[7];
   ardrgl[15]=ardrgl[3]*ardrgl[15];
   ardrgl[15]=ardrgl[14] + 9./2. + ardrgl[15];
   ardrgl[22]=ardrgl[9]*ardrgl[3];
   ardrgl[15]=11./4.*ardrgl[22] + 1./4.*ardrgl[15] + 3*ardrgl[10] + 
   ardrgl[21];
   ardrgl[22]=1./2.*ardrgl[9];
   ardrgl[15]=ardrgl[15]*ardrgl[22];
   ardrgl[23]= - 3./16.*ardrgl[10] - 1;
   ardrgl[23]=MMH*ardrgl[23];
   ardrgl[20]=ardrgl[20] + 5./2.*MMH;
   ardrgl[24]=ardrgl[11]*ardrgl[20];
   ardrgl[15]=ardrgl[19] + ardrgl[15] + 1./8.*ardrgl[24] - 3./8.*
   ardrgl[5] - 15./32.*ardrgl[18] + 1./2.*ardrgl[8] + ardrgl[23];
   ardrgl[15]=MMt*ardrgl[15];
   ardrgl[18]=ardrgl[20]*ardrgl[21];
   ardrgl[16]= - ardrgl[16] + ardrgl[18];
   ardrgl[16]=MMH*ardrgl[16];
   ardrgl[18]= - 33*ardrgl[8] + 131./4.*MMH;
   ardrgl[18]=MMH*ardrgl[18];
   ardrgl[17]=27*ardrgl[17] + ardrgl[18];
   ardrgl[16]=1./2.*ardrgl[17] + ardrgl[16];
   ardrgl[17]= - ardrgl[10]*MMH;
   ardrgl[17]=ardrgl[22] - 5./2.*ardrgl[8] + ardrgl[17];
   ardrgl[17]=ardrgl[9]*ardrgl[17];
   ardrgl[18]=pow(MMH,2);
   ardrgl[19]=MMt*ardrgl[3];
   ardrgl[20]=35./2. + 29*ardrgl[19];
   ardrgl[20]=MMt*ardrgl[20];
   ardrgl[20]= - 9*MMH + ardrgl[20];
   ardrgl[20]=MMt*ardrgl[20];
   ardrgl[20]=1./3.*ardrgl[18] + 1./4.*ardrgl[20];
   ardrgl[20]=ardrgl[20]*pow(Pi,2);
   ardrgl[16]=ardrgl[20] + 1./2.*ardrgl[16] + 3*ardrgl[17];
   ardrgl[17]=ardrgl[13]*MMH*pow(ardrgl[9],2);
   ardrgl[19]= - 3./2. + ardrgl[19];
   ardrgl[19]=MMt*ardrgl[19];
   ardrgl[19]=1./2.*MMH + ardrgl[19];
   ardrgl[19]=ardrgl[6]*ardrgl[19];
   ardrgl[18]=S2*ardrgl[18];
   ardrgl[15]= - 243./32.*ardrgl[18] + 9./8.*ardrgl[19] - 9./32.*
   ardrgl[17] + 3*ardrgl[15] + 1./8.*ardrgl[16];

      drglret = ardrgl[15]*pow(ardrgl[4],2)*pow(ardrgl[2],4)*ardrgl[1];
      return drglret;
}
