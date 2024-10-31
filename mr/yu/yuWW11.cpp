#include <WW.hpp>
namespace mr
{
  double WW<OS>::y11(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<double> aryuWW[22], yuWWret;

    aryuWW[1]=double(nH);
    aryuWW[2]=pow(CW,-1);
    aryuWW[3]=pow(MMZ,-1);
    aryuWW[4]=pow(SW,-1);
    aryuWW[5]=Tsil::I2(0,0,MMt,mu2);
    aryuWW[6]=Tsil::B(0,MMt,MMW,mu2);
    aryuWW[7]=Tsil::A(MMt,mu2);
    aryuWW[8]=pow(MMt,-1);
    aryuWW[9]=Tsil::Aeps(MMt,mu2);
    aryuWW[10]=prot00tt0->M(0);
    aryuWW[11]=prot00tt0->Tuxv(0);
    aryuWW[12]=double(nL);
    aryuWW[13]=std::real(Tsil::B(0,0,MMW,mu2));
    aryuWW[14]=prot00000->M(0);
    aryuWW[15]= - 1 + 1./3.*aryuWW[6];
    aryuWW[16]=10*aryuWW[6];
    aryuWW[15]=aryuWW[15]*aryuWW[16];
    aryuWW[15]=aryuWW[15] + 23./3. + 8./3.*aryuWW[11];
    aryuWW[15]=aryuWW[15]*MMt;
    aryuWW[16]= - 5 + 2*aryuWW[6];
    aryuWW[17]=5./3.*aryuWW[7];
    aryuWW[16]=aryuWW[16]*aryuWW[17];
    aryuWW[17]= - aryuWW[5] + 1./3.*aryuWW[9];
    aryuWW[16]= - aryuWW[16] + 2*aryuWW[17];
    aryuWW[17]=pow(aryuWW[7],2);
    aryuWW[18]=aryuWW[8]*aryuWW[17];
    aryuWW[15]=4*aryuWW[18] + aryuWW[15] - 2*aryuWW[16];
    aryuWW[16]=pow(aryuWW[2],2);
    aryuWW[18]=pow(aryuWW[4],2);
    aryuWW[19]=aryuWW[16] + aryuWW[18];
    aryuWW[15]=aryuWW[19]*aryuWW[15];
    aryuWW[19]= - aryuWW[16] - 1;
    aryuWW[16]=aryuWW[16]*aryuWW[19];
    aryuWW[16]=aryuWW[16] - aryuWW[18];
    aryuWW[19]=MMt*aryuWW[10];
    aryuWW[20]= - 2./3.*aryuWW[19] + 4./3.*aryuWW[11] + 3*aryuWW[6];
    aryuWW[20]=aryuWW[20]*MMt;
    aryuWW[21]=5*aryuWW[7] + 4*aryuWW[9];
    aryuWW[20]=aryuWW[20] + 1./3.*aryuWW[21];
    aryuWW[20]=aryuWW[20]*MMt;
    aryuWW[17]=aryuWW[20] - 2./3.*aryuWW[17];
    aryuWW[16]=aryuWW[3]*aryuWW[17]*aryuWW[16];
    aryuWW[15]=2*aryuWW[16] + aryuWW[15];
    aryuWW[15]=aryuWW[3]*aryuWW[1]*aryuWW[15];
    aryuWW[16]=1 + 2./3.*aryuWW[6];
    aryuWW[16]=aryuWW[6]*aryuWW[16];
    aryuWW[16]=aryuWW[19] - aryuWW[16];
    aryuWW[17]= - 1 + aryuWW[6];
    aryuWW[17]=aryuWW[7]*aryuWW[17];
    aryuWW[17]=aryuWW[9] + aryuWW[17];
    aryuWW[17]=aryuWW[8]*aryuWW[17];
    aryuWW[17]=aryuWW[17] + aryuWW[11];
    aryuWW[16]=13 - 4*aryuWW[16] + 16./3.*aryuWW[17];
    aryuWW[16]=aryuWW[1]*aryuWW[16];
    aryuWW[17]=4*aryuWW[13] + 31./3.;
    aryuWW[17]=aryuWW[12]*aryuWW[17];
    aryuWW[16]=aryuWW[16] + aryuWW[17];
    aryuWW[16]=aryuWW[18]*aryuWW[16];
    aryuWW[17]=aryuWW[14]*aryuWW[12];
    aryuWW[19]=aryuWW[1]*aryuWW[10];
    aryuWW[17]=aryuWW[19] + aryuWW[17];
    aryuWW[18]=aryuWW[18] - 1;
    aryuWW[17]=MMZ*aryuWW[18]*aryuWW[17];

    yuWWret = aryuWW[15] + aryuWW[16] + 8./3.*aryuWW[17];
    return yuWWret.real();
  }
} // namespace mr
