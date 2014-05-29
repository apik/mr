#include <HH.hpp>
std::complex<long double> HH::lam10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mlam[21];

    mlam[1]=pow(CW,-1);
    mlam[2]=pow(MMH,-1);
    mlam[3]=pow(MMZ,-1);
    mlam[4]=pow(SW,-1);
    mlam[5]=Tsil::B(MMH,MMH,MMH,mu2);
    mlam[6]=Tsil::B(MMZ,MMZ,MMH,mu2);
    mlam[7]=Tsil::B(MMW,MMW,MMH,mu2);
    mlam[8]=Tsil::B(MMt,MMt,MMH,mu2);
    mlam[9]=Tsil::A(MMZ,mu2);
    mlam[10]=Tsil::A(MMW,mu2);
    mlam[11]=Tsil::A(MMt,mu2);
    mlam[12]=1/( - MMW + MMH);
    mlam[13]=Tsil::A(MMH,mu2);
   mlam[14]=1 + 9*mlam[5];
   mlam[15]=1./2.*mlam[6];
   mlam[14]=mlam[15] + mlam[7] + 1./2.*mlam[14];
   mlam[16]=1./4.*MMH;
   mlam[14]=mlam[14]*mlam[16];
   mlam[16]=mlam[8]*MMt;
   mlam[16]=mlam[16] - mlam[11];
   mlam[17]=pow(MMt,2)*mlam[2]*mlam[8];
   mlam[14]=mlam[14] - 3./4.*MMt + 2*mlam[10] + mlam[9] - 6*mlam[17] + 
   3./2.*mlam[16];
   mlam[16]=pow(mlam[4],2);
   mlam[17]= - mlam[9] + mlam[10];
   mlam[17]=mlam[17]*mlam[16];
   mlam[17]=3./4.*mlam[17] + mlam[14];
   mlam[17]=mlam[3]*mlam[17];
   mlam[18]=MMZ*mlam[2];
   mlam[19]=3*mlam[18];
   mlam[20]=mlam[19] - 1;
   mlam[15]=mlam[20]*mlam[15];
   mlam[20]= - 1 + mlam[7];
   mlam[19]=mlam[20]*mlam[19];
   mlam[20]=mlam[13] - mlam[10];
   mlam[20]=mlam[12]*mlam[20];
   mlam[17]=3./4.*mlam[20] + mlam[15] + mlam[19] + 3./8. - mlam[7] + 
   mlam[17];
   mlam[16]=mlam[17]*mlam[16];
   mlam[14]=mlam[3]*mlam[14];
   mlam[14]=mlam[14] + mlam[15] + 1./8. - mlam[18];
   mlam[14]=mlam[14]*pow(mlam[1],2);
   mlam[15]=2 - 3*mlam[7];
   mlam[15]=mlam[15]*mlam[18];

      return mlam[14] + mlam[15] + mlam[16];
}
