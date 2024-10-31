#include <HH.hpp>
namespace mr
{
  double HH<OS>::y11(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<double> aryuHH[21], yuHHret;

    aryuHH[1]=double(nH);
    aryuHH[2]=pow(CW,-1);
    aryuHH[3]=pow(MMH,-1);
    aryuHH[4]=pow(MMZ,-1);
    aryuHH[5]=pow(SW,-1);
    aryuHH[6]=Tsil::I2(0,0,MMt,mu2);
    aryuHH[7]=Tsil::B(MMt,MMt,MMH,mu2);
    aryuHH[8]=Tsil::A(MMt,mu2);
    aryuHH[9]=Tsil::Beps(MMt,MMt,MMH,mu2);
    aryuHH[10]=pow(MMt,-1);
    aryuHH[11]=Tsil::Aeps(MMt,mu2);
    aryuHH[12]=prottttt0->M(0);
    aryuHH[13]=prottttt0->Vxzuv(0);
    aryuHH[14]=prottttt0->Suxv(0);
    aryuHH[15]=4 + 3*aryuHH[7];
    aryuHH[16]=2*aryuHH[7];
    aryuHH[15]=aryuHH[15]*aryuHH[16];
    aryuHH[16]=2*aryuHH[13];
    aryuHH[17]=aryuHH[16] + 3*aryuHH[12];
    aryuHH[18]=8*MMt;
    aryuHH[17]=aryuHH[17]*aryuHH[18];
    aryuHH[15]= - aryuHH[17] + aryuHH[15] + 85./2.;
    aryuHH[15]=aryuHH[15]*MMt;
    aryuHH[17]=aryuHH[3]*MMt;
    aryuHH[19]=aryuHH[17]*aryuHH[14];
    aryuHH[20]=aryuHH[10]*pow(aryuHH[8],2);
    aryuHH[19]= - aryuHH[6] + aryuHH[19] + aryuHH[20];
    aryuHH[16]=aryuHH[16] + aryuHH[12];
    aryuHH[16]=aryuHH[16]*aryuHH[18];
    aryuHH[18]=10 + 7*aryuHH[7];
    aryuHH[18]=aryuHH[18]*aryuHH[7];
    aryuHH[16]=4*aryuHH[9] - aryuHH[18] + aryuHH[16] - 25;
    aryuHH[18]=4*aryuHH[3];
    aryuHH[16]=aryuHH[18]*aryuHH[16]*pow(MMt,2);
    aryuHH[18]= - 6 + 36*aryuHH[17];
    aryuHH[18]=aryuHH[18]*aryuHH[7];
    aryuHH[18]=aryuHH[18] + 5;
    aryuHH[20]=2*aryuHH[8];
    aryuHH[18]=aryuHH[18]*aryuHH[20];
    aryuHH[17]= - 1 + 6*aryuHH[17];
    aryuHH[20]=4*aryuHH[11];
    aryuHH[17]=aryuHH[17]*aryuHH[20];
    aryuHH[20]=4*MMt;
    aryuHH[20]=aryuHH[20]*MMH*aryuHH[12];
    aryuHH[15]=aryuHH[16] + aryuHH[15] - aryuHH[18] + aryuHH[17] + 
      aryuHH[20] - 4*aryuHH[19];
    aryuHH[16]=pow(aryuHH[5],2);
    aryuHH[17]=pow(aryuHH[2],2);
    aryuHH[16]=aryuHH[16] + aryuHH[17];

    yuHHret = aryuHH[16]*aryuHH[15]*aryuHH[4]*aryuHH[1];
    return yuHHret.real();
  }
} // namespace mr
