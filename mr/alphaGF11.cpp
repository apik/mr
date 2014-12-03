#include <alphaGF.hpp>
long double alphaGF::a11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aralphaGF[42], alphaGFret;

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
    aralphaGF[13]=prottttt0->M(0);
    aralphaGF[14]=prot00000->M(0);
    aralphaGF[15]=prottttt0->Vzxyv(0);
    aralphaGF[16]=prottttt0->Suxv(0);
    aralphaGF[17]=prot00tt0->M(0);
    aralphaGF[18]=prot00tt0->Tuxv(0);
    aralphaGF[19]=double(nL);
    aralphaGF[20]=std::real(Tsil::B(0,0,MMW,mu2));
    aralphaGF[21]=1/(4*MMt - MMZ);
   aralphaGF[22]=pow(aralphaGF[4],2);
   aralphaGF[23]=aralphaGF[22] + 16./3.;
   aralphaGF[23]=aralphaGF[23]*aralphaGF[22];
   aralphaGF[23]=aralphaGF[23] - 64./9.;
   aralphaGF[24]=aralphaGF[23]*aralphaGF[6];
   aralphaGF[25]= - 32 - aralphaGF[22];
   aralphaGF[25]=aralphaGF[25]*aralphaGF[22];
   aralphaGF[25]=128./3. + aralphaGF[25];
   aralphaGF[25]=2*aralphaGF[25] - aralphaGF[24];
   aralphaGF[26]=2./3.*aralphaGF[6];
   aralphaGF[25]=aralphaGF[25]*aralphaGF[26];
   aralphaGF[27]=4*aralphaGF[7];
   aralphaGF[28]=aralphaGF[6]*aralphaGF[27];
   aralphaGF[28]=aralphaGF[28] - 2*aralphaGF[16];
   aralphaGF[29]=pow(aralphaGF[4],4);
   aralphaGF[28]=aralphaGF[29]*aralphaGF[28];
   aralphaGF[30]=pow(aralphaGF[2],4);
   aralphaGF[31]=13*aralphaGF[29] - 5*aralphaGF[30];
   aralphaGF[31]=aralphaGF[7]*aralphaGF[31];
   aralphaGF[28]=aralphaGF[31] + aralphaGF[28];
   aralphaGF[31]=2./3.*aralphaGF[3];
   aralphaGF[28]=aralphaGF[28]*aralphaGF[31];
   aralphaGF[32]=37./9. - 25./2.*aralphaGF[22];
   aralphaGF[32]=aralphaGF[32]*aralphaGF[22];
   aralphaGF[33]=4./3.*aralphaGF[23];
   aralphaGF[33]=aralphaGF[10]*aralphaGF[33];
   aralphaGF[34]=aralphaGF[22] - 1;
   aralphaGF[34]=aralphaGF[34]*aralphaGF[22];
   aralphaGF[35]=pow(aralphaGF[2],2);
   aralphaGF[34]=aralphaGF[34] - aralphaGF[35];
   aralphaGF[36]=10*aralphaGF[8];
   aralphaGF[36]=aralphaGF[34]*aralphaGF[36];
   aralphaGF[37]= - 1./3.*aralphaGF[8] + 1;
   aralphaGF[37]=aralphaGF[37]*aralphaGF[36];
   aralphaGF[34]=aralphaGF[18]*aralphaGF[34];
   aralphaGF[25]=aralphaGF[28] + aralphaGF[25] - 8./3.*aralphaGF[34] + 
   aralphaGF[37] + 23./3.*aralphaGF[35] + aralphaGF[33] + 128./27. + 
   aralphaGF[32];
   aralphaGF[25]=aralphaGF[3]*aralphaGF[25];
   aralphaGF[28]=2 + aralphaGF[6];
   aralphaGF[26]=aralphaGF[26]*aralphaGF[29]*aralphaGF[28];
   aralphaGF[28]=aralphaGF[30] - aralphaGF[29];
   aralphaGF[32]=MMt*aralphaGF[17]*aralphaGF[28];
   aralphaGF[32]=aralphaGF[29] - aralphaGF[32];
   aralphaGF[33]= - 4./3.*aralphaGF[18] - 3*aralphaGF[8];
   aralphaGF[33]=aralphaGF[28]*aralphaGF[33];
   aralphaGF[26]=aralphaGF[26] + aralphaGF[33] - 2./3.*aralphaGF[32];
   aralphaGF[26]=MMt*aralphaGF[26];
   aralphaGF[30]=2*aralphaGF[29] - aralphaGF[30];
   aralphaGF[30]=aralphaGF[11]*aralphaGF[30];
   aralphaGF[26]=2*aralphaGF[26] + 8./3.*aralphaGF[30];
   aralphaGF[26]=aralphaGF[26]*pow(aralphaGF[3],2);
   aralphaGF[30]=aralphaGF[22] - 2;
   aralphaGF[30]=aralphaGF[30]*aralphaGF[22];
   aralphaGF[32]=aralphaGF[17]*aralphaGF[30];
   aralphaGF[25]=4*aralphaGF[32] + aralphaGF[25] + aralphaGF[26];
   aralphaGF[25]=MMt*aralphaGF[25];
   aralphaGF[26]= - 80./3. + 7*aralphaGF[22];
   aralphaGF[26]=aralphaGF[26]*aralphaGF[22];
   aralphaGF[32]= - 16 + 3*aralphaGF[22];
   aralphaGF[32]=aralphaGF[32]*aralphaGF[22];
   aralphaGF[32]=64./3. + aralphaGF[32];
   aralphaGF[33]=aralphaGF[21]*MMZ;
   aralphaGF[32]=aralphaGF[32]*aralphaGF[33];
   aralphaGF[34]=aralphaGF[22] - 16./3.;
   aralphaGF[34]=aralphaGF[34]*aralphaGF[22];
   aralphaGF[34]=aralphaGF[34] + 64./9.;
   aralphaGF[37]=aralphaGF[34]*aralphaGF[21];
   aralphaGF[38]= - aralphaGF[27]*aralphaGF[37];
   aralphaGF[33]=aralphaGF[33] + 1;
   aralphaGF[33]=aralphaGF[33]*aralphaGF[34];
   aralphaGF[34]= - aralphaGF[6]*aralphaGF[33];
   aralphaGF[26]=aralphaGF[34] + aralphaGF[38] + aralphaGF[32] + 320./9.
    + aralphaGF[26];
   aralphaGF[26]=aralphaGF[6]*aralphaGF[26];
   aralphaGF[32]=aralphaGF[30]*aralphaGF[9];
   aralphaGF[34]= - 7./3. + aralphaGF[22];
   aralphaGF[34]=aralphaGF[34]*aralphaGF[22];
   aralphaGF[34]= - aralphaGF[35] + 128./9. + 5*aralphaGF[34];
   aralphaGF[34]=aralphaGF[3]*aralphaGF[34];
   aralphaGF[34]=1./3.*aralphaGF[34] - 4./3.*aralphaGF[32] - 
   aralphaGF[37];
   aralphaGF[34]=aralphaGF[11]*aralphaGF[34];
   aralphaGF[38]=aralphaGF[22] - 8./3.;
   aralphaGF[39]=aralphaGF[38]*aralphaGF[22];
   aralphaGF[39]=aralphaGF[39] + 32./9.;
   aralphaGF[40]=MMZ*aralphaGF[39];
   aralphaGF[23]=aralphaGF[23]*MMt;
   aralphaGF[41]=aralphaGF[31]*aralphaGF[23];
   aralphaGF[29]= - aralphaGF[29] + aralphaGF[41];
   aralphaGF[29]=MMt*aralphaGF[29];
   aralphaGF[29]=1./3.*aralphaGF[40] + aralphaGF[29];
   aralphaGF[29]=aralphaGF[13]*aralphaGF[29];
   aralphaGF[40]= - 1 - 2./3.*aralphaGF[8];
   aralphaGF[40]=aralphaGF[8]*aralphaGF[40]*aralphaGF[30];
   aralphaGF[29]=aralphaGF[29] + aralphaGF[40] + aralphaGF[34];
   aralphaGF[34]= - 299./3. + 32*aralphaGF[22];
   aralphaGF[34]=aralphaGF[34]*aralphaGF[22];
   aralphaGF[34]= - aralphaGF[36] - 25*aralphaGF[35] + 896./9. + 
   aralphaGF[34];
   aralphaGF[36]=1 - 2*aralphaGF[22];
   aralphaGF[36]=aralphaGF[36]*aralphaGF[22];
   aralphaGF[36]=aralphaGF[36] + aralphaGF[35];
   aralphaGF[36]=aralphaGF[7]*aralphaGF[9]*aralphaGF[36];
   aralphaGF[24]= - 4./3.*aralphaGF[24] + 1./3.*aralphaGF[34] + 2*
   aralphaGF[36];
   aralphaGF[24]=aralphaGF[7]*aralphaGF[24];
   aralphaGF[34]= - aralphaGF[16]*aralphaGF[38];
   aralphaGF[34]=aralphaGF[5] + aralphaGF[34];
   aralphaGF[34]=aralphaGF[22]*aralphaGF[34];
   aralphaGF[35]=aralphaGF[5]*aralphaGF[35];
   aralphaGF[34]=aralphaGF[35] - 32./9.*aralphaGF[16] + aralphaGF[34];
   aralphaGF[28]=aralphaGF[28]*pow(aralphaGF[7],2)*aralphaGF[31];
   aralphaGF[24]=aralphaGF[28] + 2*aralphaGF[34] + aralphaGF[24];
   aralphaGF[24]=aralphaGF[3]*aralphaGF[24];
   aralphaGF[23]=aralphaGF[3]*aralphaGF[23];
   aralphaGF[23]=aralphaGF[23] - aralphaGF[39];
   aralphaGF[23]=aralphaGF[15]*MMt*aralphaGF[23];
   aralphaGF[28]=aralphaGF[18]*aralphaGF[30];
   aralphaGF[23]=aralphaGF[28] - aralphaGF[23];
   aralphaGF[28]=1 - aralphaGF[8];
   aralphaGF[28]=aralphaGF[28]*aralphaGF[32];
   aralphaGF[28]=4./3.*aralphaGF[28] + aralphaGF[37];
   aralphaGF[27]=aralphaGF[28]*aralphaGF[27];
   aralphaGF[28]=aralphaGF[22] - 4./3.;
   aralphaGF[28]=aralphaGF[28]*aralphaGF[22];
   aralphaGF[28]=aralphaGF[28] + 8./9.;
   aralphaGF[31]=aralphaGF[14]*MMZ;
   aralphaGF[32]=4./3.*aralphaGF[31] + 2*aralphaGF[12];
   aralphaGF[28]=aralphaGF[28]*aralphaGF[32];
   aralphaGF[32]= - aralphaGF[10]*aralphaGF[33];
   aralphaGF[33]= - 4 - 1./2.*aralphaGF[22];
   aralphaGF[33]=aralphaGF[33]*aralphaGF[22];
   aralphaGF[33]=860./9. + aralphaGF[33];
   aralphaGF[34]=1 - 1./3.*aralphaGF[22];
   aralphaGF[34]=aralphaGF[34]*aralphaGF[22];
   aralphaGF[34]= - 2./3. + aralphaGF[34];
   aralphaGF[34]=MMZ*aralphaGF[17]*aralphaGF[34];
   aralphaGF[23]=aralphaGF[25] + 2*aralphaGF[24] + aralphaGF[26] + 
   aralphaGF[27] + 8*aralphaGF[34] + 1./3.*aralphaGF[33] + 
   aralphaGF[32] + aralphaGF[28] - 16./3.*aralphaGF[23] + 4*
   aralphaGF[29];
   aralphaGF[23]=aralphaGF[1]*aralphaGF[23];
   aralphaGF[24]=20./9. + aralphaGF[30];
   aralphaGF[24]=aralphaGF[12]*aralphaGF[24];
   aralphaGF[25]= - aralphaGF[20]*aralphaGF[30];
   aralphaGF[22]=2./9. + aralphaGF[22];
   aralphaGF[22]=aralphaGF[22]*aralphaGF[31];
   aralphaGF[22]=2./3.*aralphaGF[22] + aralphaGF[25] + 155./27. + 
   aralphaGF[24];
   aralphaGF[22]=aralphaGF[19]*aralphaGF[22];

      alphaGFret = 4*aralphaGF[22] + aralphaGF[23];
      return alphaGFret.real();
}
