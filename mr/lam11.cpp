#include <HH.hpp>
std::complex<long double> HH::lam11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mlam[19];

    mlam[1]=pow(CW,-1);
    mlam[2]=pow(MMH,-1);
    mlam[3]=pow(MMZ,-1);
    mlam[4]=pow(SW,-1);
    mlam[5]=Tsil::I2(0,0,MMt,mu2);
    mlam[6]=Tsil::B(MMt,MMt,MMH,mu2);
    mlam[7]=Tsil::A(MMt,mu2);
    mlam[8]=Tsil::Beps(MMt,MMt,MMH,mu2);
    mlam[9]=pow(MMt,-1);
    mlam[10]=Tsil::Aeps(MMt,mu2);
    mlam[11]=prottttt0->M(0);
    mlam[12]=prottttt0->Vxzuv(0);
    mlam[13]=prottttt0->Suxv(0);
   mlam[14]=mlam[11] + 2*mlam[12];
   mlam[15]=8*MMt;
   mlam[14]=mlam[15]*mlam[14];
   mlam[15]=10 + 7*mlam[6];
   mlam[15]=mlam[15]*mlam[6];
   mlam[14]=mlam[14] - mlam[15] - 25;
   mlam[14]=mlam[2]*mlam[14];
   mlam[15]=4*mlam[2];
   mlam[15]=mlam[15]*mlam[8];
   mlam[14]=mlam[15] - 6*mlam[11] - 4*mlam[12] + mlam[14];
   mlam[15]=4*MMt;
   mlam[14]=mlam[14]*mlam[15];
   mlam[16]=4 + 3*mlam[6];
   mlam[17]=2*mlam[6];
   mlam[16]=mlam[16]*mlam[17];
   mlam[17]=mlam[7]*mlam[6];
   mlam[17]= - 24*mlam[10] + 72*mlam[17];
   mlam[17]=mlam[17]*mlam[2];
   mlam[14]=mlam[14] - mlam[17] + mlam[16] + 85./2.;
   mlam[14]=mlam[14]*MMt;
   mlam[16]=mlam[10] - mlam[5];
   mlam[17]= - 5 + 6*mlam[6];
   mlam[17]=mlam[17]*mlam[7];
   mlam[16]=mlam[17] - 2*mlam[16];
   mlam[17]=mlam[13]*mlam[2];
   mlam[18]=MMH*mlam[11];
   mlam[17]=mlam[17] - mlam[18];
   mlam[15]=mlam[17]*mlam[15];
   mlam[17]=mlam[9]*pow(mlam[7],2);
   mlam[14]=mlam[14] - mlam[15] - 4*mlam[17] + 2*mlam[16];
   mlam[15]=pow(mlam[4],2);
   mlam[16]=pow(mlam[1],2);
   mlam[15]=mlam[15] + mlam[16];

      return mlam[15]*mlam[14]*mlam[3];
}
