#include <dr.hpp>
namespace mr
{
  long double dr<MS>::dr11(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<long double> ardrbar[13], drbarret;

    ardrbar[1]=double(nH);
    ardrbar[2]=pow(mmZ,-1);
    ardrbar[3]=pow(mmH,-1);
    ardrbar[4]=pow(s,-1);
    ardrbar[5]=pow(c,-1);
    ardrbar[6]=Tsil::I2(0,0,mmt,mu2);
    ardrbar[7]=Tsil::A(mmt,mu2);
    ardrbar[8]=pow(mmt,-1);
    ardrbar[9]=Tsil::Aeps(mmt,mu2);
    ardrbar[10]=pow(ardrbar[5],2);
    ardrbar[11]=pow(ardrbar[4],2);
    ardrbar[10]=ardrbar[10] + ardrbar[11];
    ardrbar[11]= - ardrbar[8] + 6*ardrbar[3];
    ardrbar[12]=2*ardrbar[7];
    ardrbar[11]=ardrbar[11]*ardrbar[12];
    ardrbar[11]=ardrbar[11] - 1;
    ardrbar[11]=ardrbar[11]*ardrbar[7];
    ardrbar[11]= - ardrbar[6] + ardrbar[11] + ardrbar[9];
    ardrbar[12]=16*ardrbar[7] + 48*mmt;
    ardrbar[12]=ardrbar[3]*ardrbar[12];
    ardrbar[12]= - 25./2. + ardrbar[12];
    ardrbar[12]=ardrbar[12]*mmt;
    ardrbar[11]=ardrbar[12] + 4*ardrbar[11];

    drbarret = ardrbar[11]*ardrbar[10]*ardrbar[2]*ardrbar[1];
    return drbarret.real();
  }
} // namespace mr
