#include <dr.hpp>
long double dr<OS>::dr11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> ardr[13], drret;

    ardr[1]=double(nH);
    ardr[2]=pow(CW,-1);
    ardr[3]=pow(MMH,-1);
    ardr[4]=pow(MMZ,-1);
    ardr[5]=pow(SW,-1);
    ardr[6]=Tsil::I2(0,0,MMt,mu2);
    ardr[7]=Tsil::A(MMt,mu2);
    ardr[8]=pow(MMt,-1);
    ardr[9]=Tsil::Aeps(MMt,mu2);
   ardr[10]= - ardr[8] + 12*ardr[3];
   ardr[11]=2*ardr[7];
   ardr[10]=ardr[10]*ardr[11];
   ardr[10]=ardr[10] - 5;
   ardr[10]=ardr[11]*ardr[10];
   ardr[11]=ardr[9] - ardr[6];
   ardr[12]=MMt*ardr[3];
   ardr[12]= - 37./2. + 64*ardr[12];
   ardr[12]=ardr[12]*MMt;
   ardr[10]=ardr[10] - ardr[12] - 4*ardr[11];
   ardr[11]= - pow(ardr[2],2);
   ardr[12]= - pow(ardr[5],2);
   ardr[11]=ardr[11] + ardr[12];

      drret = ardr[11]*ardr[10]*ardr[4]*ardr[1];
      return drret.real();
}
