#include <HH.hpp>
std::complex<long double> HH::my10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mlam[26];

    mlam[1]=pow(CW,-1);
    mlam[2]=pow(MMH,-1);
    mlam[3]=pow(MMZ,-1);
    mlam[4]=pow(SW,-1);
    mlam[5]=double(nH);
    mlam[6]=Tsil::B(MMt,MMt,MMH,mu2);
    mlam[7]=Tsil::B(MMb,MMb,MMH,mu2);
    mlam[8]=Tsil::A(MMt,mu2);
    mlam[9]=Tsil::B(MMH,MMH,MMH,mu2);
    mlam[10]=Tsil::B(MMZ,MMZ,MMH,mu2);
    mlam[11]=Tsil::B(MMW,MMW,MMH,mu2);
    mlam[12]=Tsil::A(MMZ,mu2);
    mlam[13]=Tsil::A(MMW,mu2);
    mlam[14]=1/( - MMb + MMt);
    mlam[15]=Tsil::A(MMb,mu2);
    mlam[16]=1/( - MMW + MMH);
    mlam[17]=Tsil::A(MMH,mu2);
   mlam[18]=MMZ*mlam[2];
   mlam[19]=3*mlam[18];
   mlam[20]=mlam[19] - 1;
   mlam[21]=1./2.*mlam[10];
   mlam[20]=mlam[20]*mlam[21];
   mlam[21]=mlam[11] + 1./2.;
   mlam[21]=1./8.*mlam[10] + 1./4.*mlam[21];
   mlam[21]=mlam[21]*MMH;
   mlam[21]=mlam[21] + mlam[12];
   mlam[21]=mlam[21]*mlam[3];
   mlam[22]=MMH*mlam[9]*mlam[3];
   mlam[23]=mlam[13]*mlam[3];
   mlam[20]=9./8.*mlam[22] + mlam[21] + mlam[20] + 2*mlam[23];
   mlam[21]=1./2.*MMt;
   mlam[22]=mlam[21] + mlam[8];
   mlam[24]= - 1./2.*MMb + mlam[22] - mlam[15];
   mlam[24]=MMb*mlam[24]*mlam[14];
   mlam[22]=mlam[24] + mlam[22];
   mlam[24]=2*mlam[2];
   mlam[25]=mlam[24]*pow(MMt,2);
   mlam[21]=mlam[25] - mlam[21];
   mlam[21]=mlam[21]*mlam[6];
   mlam[24]=mlam[24]*MMb;
   mlam[24]=mlam[24] - 1./2.;
   mlam[24]=mlam[24]*mlam[7]*MMb;
   mlam[21]=mlam[21] + mlam[24] + 1./2.*mlam[22];
   mlam[22]=3*mlam[5];
   mlam[21]=mlam[22]*mlam[21]*mlam[3];
   mlam[22]= - mlam[3]*mlam[12];
   mlam[22]=mlam[22] + mlam[23];
   mlam[23]=pow(mlam[4],2);
   mlam[22]=mlam[22]*mlam[23];
   mlam[24]=mlam[17] - mlam[13];
   mlam[24]=mlam[16]*mlam[24];
   mlam[22]=mlam[24] + mlam[22];
   mlam[24]= - 1 + mlam[11];
   mlam[19]=mlam[24]*mlam[19];
   mlam[19]= - mlam[21] + mlam[19] + 3./8. - mlam[11] + mlam[20] + 3./4.
   *mlam[22];
   mlam[19]=mlam[23]*mlam[19];
   mlam[20]= - mlam[21] + 1./8. - mlam[18] + mlam[20];
   mlam[20]=mlam[20]*pow(mlam[1],2);
   mlam[21]=2 - 3*mlam[11];
   mlam[18]=mlam[21]*mlam[18];

      return mlam[18] + mlam[19] + mlam[20];
}
