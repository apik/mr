#include <alphaGF.hpp>
namespace mr
{
  double alphaGF::a11(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<double> aralphaGF[39], alphaGFret;

    aralphaGF[1]=double(nH);
    aralphaGF[2]=pow(CW,-1);
    aralphaGF[3]=pow(MMZ,-1);
    aralphaGF[4]=pow(SW,-1);
    aralphaGF[5]=Tsil::I2(0,0,MMt,mu2);
    aralphaGF[6]=Tsil::B(MMt,MMt,MMZ,mu2);
    aralphaGF[7]=Tsil::A(MMt,mu2);
    aralphaGF[8]=Tsil::B(0,MMt,MMW,mu2);
    aralphaGF[9]=pow(MMt,-1);
    aralphaGF[10]=Tsil::Beps(MMt,MMt,MMZ,mu2);
    aralphaGF[11]=Tsil::Aeps(MMt,mu2);
    aralphaGF[12]=std::real(Tsil::B(0,0,MMZ,mu2));
    aralphaGF[13]=Wprot00tt0->M(0);
    aralphaGF[14]=Wprot00tt0->Tuxv(0);
    aralphaGF[15]=Zprottttt0->M(0);
    aralphaGF[16]=Zprot00000->M(0);
    aralphaGF[17]=Zprottttt0->Vzxyv(0);
    aralphaGF[18]=Zprottttt0->Suxv(0);
    aralphaGF[19]=double(nL);
    aralphaGF[20]=std::real(Tsil::B(0,0,MMW,mu2));
    aralphaGF[21]=Wprot00000->M(0);
    aralphaGF[22]=1/(4*MMt - MMZ);
    aralphaGF[23]=pow(aralphaGF[4],2);
    aralphaGF[24]=37./9. - 25./2.*aralphaGF[23];
    aralphaGF[24]=aralphaGF[24]*aralphaGF[23];
    aralphaGF[25]=aralphaGF[23] - 1;
    aralphaGF[25]=aralphaGF[25]*aralphaGF[23];
    aralphaGF[26]=pow(aralphaGF[2],2);
    aralphaGF[25]=aralphaGF[25] - aralphaGF[26];
    aralphaGF[27]=10*aralphaGF[8];
    aralphaGF[27]=aralphaGF[25]*aralphaGF[27];
    aralphaGF[28]= - 1./3.*aralphaGF[8] + 1;
    aralphaGF[28]=aralphaGF[28]*aralphaGF[27];
    aralphaGF[25]=aralphaGF[14]*aralphaGF[25];
    aralphaGF[29]=aralphaGF[23] + 16./3.;
    aralphaGF[29]=aralphaGF[29]*aralphaGF[23];
    aralphaGF[29]=aralphaGF[29] - 64./9.;
    aralphaGF[30]=aralphaGF[10]*aralphaGF[29];
    aralphaGF[31]=128./9. + 23*aralphaGF[26];
    aralphaGF[24]=4./3.*aralphaGF[30] - 8./3.*aralphaGF[25] + 
      aralphaGF[28] + 1./3.*aralphaGF[31] + aralphaGF[24];
    aralphaGF[24]=MMt*aralphaGF[24];
    aralphaGF[25]=pow(aralphaGF[2],4);
    aralphaGF[28]=pow(aralphaGF[4],4);
    aralphaGF[30]=aralphaGF[25] - aralphaGF[28];
    aralphaGF[31]=MMt*aralphaGF[13];
    aralphaGF[32]=2./3.*aralphaGF[31] - 4./3.*aralphaGF[14] - 3*
      aralphaGF[8];
    aralphaGF[32]=aralphaGF[30]*aralphaGF[32];
    aralphaGF[33]=2./3.*aralphaGF[28];
    aralphaGF[32]= - aralphaGF[33] + aralphaGF[32];
    aralphaGF[32]=MMt*aralphaGF[32];
    aralphaGF[34]= - aralphaGF[25] + 2*aralphaGF[28];
    aralphaGF[35]=4./3.*aralphaGF[11];
    aralphaGF[34]=aralphaGF[34]*aralphaGF[35];
    aralphaGF[33]= - aralphaGF[18]*aralphaGF[33];
    aralphaGF[32]=aralphaGF[34] + aralphaGF[33] + aralphaGF[32];
    aralphaGF[32]=MMt*aralphaGF[32];
    aralphaGF[25]= - 5*aralphaGF[25] + 13*aralphaGF[28];
    aralphaGF[25]=MMt*aralphaGF[25];
    aralphaGF[33]=2*aralphaGF[7];
    aralphaGF[30]=aralphaGF[30]*aralphaGF[33];
    aralphaGF[25]=aralphaGF[25] + aralphaGF[30];
    aralphaGF[25]=aralphaGF[7]*aralphaGF[25];
    aralphaGF[30]=pow(MMt,2);
    aralphaGF[34]=aralphaGF[7]*MMt;
    aralphaGF[34]=aralphaGF[30] + aralphaGF[34];
    aralphaGF[36]=aralphaGF[6]*aralphaGF[30];
    aralphaGF[34]=2*aralphaGF[34] + aralphaGF[36];
    aralphaGF[36]=2./3.*aralphaGF[6];
    aralphaGF[34]=aralphaGF[36]*aralphaGF[28]*aralphaGF[34];
    aralphaGF[25]=aralphaGF[34] + 1./3.*aralphaGF[25] + aralphaGF[32];
    aralphaGF[25]=aralphaGF[3]*aralphaGF[25];
    aralphaGF[32]=aralphaGF[23] - 8./3.;
    aralphaGF[32]=aralphaGF[32]*aralphaGF[23];
    aralphaGF[32]=aralphaGF[32] + 32./9.;
    aralphaGF[34]=aralphaGF[18]*aralphaGF[32];
    aralphaGF[37]=aralphaGF[26] + aralphaGF[23];
    aralphaGF[37]=aralphaGF[5]*aralphaGF[37];
    aralphaGF[34]=aralphaGF[34] - aralphaGF[37];
    aralphaGF[37]=16./3.*aralphaGF[17];
    aralphaGF[38]=8./3.*aralphaGF[15] + aralphaGF[37];
    aralphaGF[30]=aralphaGF[30]*aralphaGF[29]*aralphaGF[38];
    aralphaGF[38]= - 299./3. + 32*aralphaGF[23];
    aralphaGF[38]=aralphaGF[38]*aralphaGF[23];
    aralphaGF[27]= - aralphaGF[27] + aralphaGF[38] + 896./9. - 25*
      aralphaGF[26];
    aralphaGF[38]=1 - 2*aralphaGF[23];
    aralphaGF[38]=aralphaGF[38]*aralphaGF[23];
    aralphaGF[38]=aralphaGF[26] + aralphaGF[38];
    aralphaGF[38]=aralphaGF[9]*aralphaGF[38]*aralphaGF[33];
    aralphaGF[27]=1./3.*aralphaGF[27] + aralphaGF[38];
    aralphaGF[27]=aralphaGF[27]*aralphaGF[33];
    aralphaGF[38]= - 32 - aralphaGF[23];
    aralphaGF[38]=aralphaGF[38]*aralphaGF[23];
    aralphaGF[38]=128./3. + aralphaGF[38];
    aralphaGF[38]=MMt*aralphaGF[38];
    aralphaGF[33]= - aralphaGF[29]*aralphaGF[33];
    aralphaGF[33]=aralphaGF[38] + aralphaGF[33];
    aralphaGF[29]= - aralphaGF[6]*MMt*aralphaGF[29];
    aralphaGF[29]=2*aralphaGF[33] + aralphaGF[29];
    aralphaGF[29]=aralphaGF[29]*aralphaGF[36];
    aralphaGF[33]= - 7./3. + aralphaGF[23];
    aralphaGF[33]=aralphaGF[33]*aralphaGF[23];
    aralphaGF[26]=5*aralphaGF[33] + 128./9. - aralphaGF[26];
    aralphaGF[26]=aralphaGF[26]*aralphaGF[35];
    aralphaGF[24]=2*aralphaGF[25] + aralphaGF[26] + aralphaGF[29] + 
      aralphaGF[27] + aralphaGF[30] + aralphaGF[24] - 4*aralphaGF[34];
    aralphaGF[24]=aralphaGF[3]*aralphaGF[24];
    aralphaGF[25]= - 80./3. + 7*aralphaGF[23];
    aralphaGF[25]=aralphaGF[25]*aralphaGF[23];
    aralphaGF[26]= - 16 + 3*aralphaGF[23];
    aralphaGF[26]=aralphaGF[26]*aralphaGF[23];
    aralphaGF[26]=64./3. + aralphaGF[26];
    aralphaGF[27]=MMZ*aralphaGF[22];
    aralphaGF[26]=aralphaGF[26]*aralphaGF[27];
    aralphaGF[29]=aralphaGF[23] - 16./3.;
    aralphaGF[29]=aralphaGF[29]*aralphaGF[23];
    aralphaGF[29]=aralphaGF[29] + 64./9.;
    aralphaGF[30]=aralphaGF[29]*aralphaGF[22];
    aralphaGF[33]=4*aralphaGF[7];
    aralphaGF[34]= - aralphaGF[33]*aralphaGF[30];
    aralphaGF[27]=aralphaGF[27] + 1;
    aralphaGF[27]=aralphaGF[27]*aralphaGF[29];
    aralphaGF[29]= - aralphaGF[6]*aralphaGF[27];
    aralphaGF[25]=aralphaGF[29] + aralphaGF[34] + aralphaGF[26] + 320./9.
      + aralphaGF[25];
    aralphaGF[25]=aralphaGF[6]*aralphaGF[25];
    aralphaGF[26]=aralphaGF[23] - 2;
    aralphaGF[26]=aralphaGF[26]*aralphaGF[23];
    aralphaGF[29]=4./3.*aralphaGF[9];
    aralphaGF[29]=aralphaGF[26]*aralphaGF[29];
    aralphaGF[34]= - aralphaGF[30] - aralphaGF[29];
    aralphaGF[34]=aralphaGF[11]*aralphaGF[34];
    aralphaGF[35]= - 1 - 2./3.*aralphaGF[8];
    aralphaGF[35]=aralphaGF[8]*aralphaGF[35]*aralphaGF[26];
    aralphaGF[36]=MMZ*aralphaGF[32];
    aralphaGF[28]= - MMt*aralphaGF[28];
    aralphaGF[28]=1./3.*aralphaGF[36] + aralphaGF[28];
    aralphaGF[28]=aralphaGF[15]*aralphaGF[28];
    aralphaGF[28]=aralphaGF[34] + aralphaGF[35] + aralphaGF[28];
    aralphaGF[34]=1 - aralphaGF[8];
    aralphaGF[29]=aralphaGF[34]*aralphaGF[29];
    aralphaGF[29]=aralphaGF[30] + aralphaGF[29];
    aralphaGF[29]=aralphaGF[29]*aralphaGF[33];
    aralphaGF[30]=4*aralphaGF[31] - 16./3.*aralphaGF[14];
    aralphaGF[30]=aralphaGF[30]*aralphaGF[26];
    aralphaGF[31]= - 4 - 1./2.*aralphaGF[23];
    aralphaGF[31]=aralphaGF[31]*aralphaGF[23];
    aralphaGF[31]=860./9. + aralphaGF[31];
    aralphaGF[33]=aralphaGF[23] - 4./3.;
    aralphaGF[33]=aralphaGF[33]*aralphaGF[23];
    aralphaGF[33]=aralphaGF[33] + 8./9.;
    aralphaGF[34]=aralphaGF[16]*MMZ;
    aralphaGF[35]=aralphaGF[33]*aralphaGF[34];
    aralphaGF[36]= - 1 + 1./3.*aralphaGF[23];
    aralphaGF[23]=aralphaGF[36]*aralphaGF[23];
    aralphaGF[23]=aralphaGF[23] + 2./3.;
    aralphaGF[23]=aralphaGF[23]*MMZ;
    aralphaGF[36]=aralphaGF[13]*aralphaGF[23];
    aralphaGF[27]= - aralphaGF[10]*aralphaGF[27];
    aralphaGF[32]= - MMt*aralphaGF[32]*aralphaGF[37];
    aralphaGF[24]=aralphaGF[24] + aralphaGF[32] + aralphaGF[25] + 
      aralphaGF[27] + aralphaGF[29] - 8*aralphaGF[36] + 4./3.*
      aralphaGF[35] + 1./3.*aralphaGF[31] + aralphaGF[30] + 4*
      aralphaGF[28];
    aralphaGF[24]=aralphaGF[1]*aralphaGF[24];
    aralphaGF[25]= - aralphaGF[20]*aralphaGF[26];
    aralphaGF[26]=aralphaGF[26] + 20./9.;
    aralphaGF[27]=aralphaGF[26]*aralphaGF[34];
    aralphaGF[23]=aralphaGF[21]*aralphaGF[23];
    aralphaGF[23]= - 2*aralphaGF[23] + 2./3.*aralphaGF[27] + 155./27. + 
      aralphaGF[25];
    aralphaGF[23]=aralphaGF[19]*aralphaGF[23];
    aralphaGF[25]=aralphaGF[19]*aralphaGF[26];
    aralphaGF[26]=aralphaGF[1]*aralphaGF[33];
    aralphaGF[25]=2*aralphaGF[25] + aralphaGF[26];
    aralphaGF[25]=aralphaGF[12]*aralphaGF[25];

    alphaGFret = 4*aralphaGF[23] + aralphaGF[24] + 2*aralphaGF[25];
    return alphaGFret.real();
  }
} // namespace mr
