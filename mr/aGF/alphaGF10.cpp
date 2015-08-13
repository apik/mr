#include <alphaGF.hpp>
namespace mr
{
  long double alphaGF::a10(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<long double> aralphaGF[34], alphaGFret;

    aralphaGF[1]=double(nL + nH);
    aralphaGF[2]=std::real(Tsil::B(0,0,MMZ,mu2));
    aralphaGF[3]=pow(SW,-1);
    aralphaGF[4]=std::real(Tsil::B(0,0,MMW,mu2));
    aralphaGF[5]=double(nH);
    aralphaGF[6]=pow(CW,-1);
    aralphaGF[7]=pow(MMZ,-1);
    aralphaGF[8]=Tsil::B(MMt,MMt,MMZ,mu2);
    aralphaGF[9]=Tsil::B(0,MMt,MMW,mu2);
    aralphaGF[10]=Tsil::A(MMt,mu2);
    aralphaGF[11]=double(nL);
    aralphaGF[12]=double(boson);
    aralphaGF[13]=Tsil::B(MMZ,MMH,MMZ,mu2);
    aralphaGF[14]=Tsil::B(MMW,MMH,MMW,mu2);
    aralphaGF[15]=Tsil::B(MMW,MMZ,MMW,mu2);
    aralphaGF[16]=Tsil::B(MMW,MMW,MMZ,mu2);
    aralphaGF[17]=Tsil::A(MMH,mu2);
    aralphaGF[18]=Tsil::A(MMZ,mu2);
    aralphaGF[19]=Tsil::A(MMW,mu2);
    aralphaGF[20]=1/( - MMW + MMH);
    aralphaGF[21]=pow(aralphaGF[6],2);
    aralphaGF[22]=pow(aralphaGF[3],2);
    aralphaGF[23]=aralphaGF[21] + aralphaGF[22];
    aralphaGF[24]=aralphaGF[22] - 1;
    aralphaGF[24]=aralphaGF[24]*aralphaGF[22];
    aralphaGF[24]=aralphaGF[24] - aralphaGF[21];
    aralphaGF[25]=aralphaGF[14]*aralphaGF[24];
    aralphaGF[26]=pow(aralphaGF[3],4);
    aralphaGF[27]=aralphaGF[26]*aralphaGF[13];
    aralphaGF[25]=aralphaGF[25] - 1./8.*aralphaGF[23] - aralphaGF[27];
    aralphaGF[25]=MMH*aralphaGF[25];
    aralphaGF[28]=2*aralphaGF[22];
    aralphaGF[29]= - 31./4. + aralphaGF[28];
    aralphaGF[29]=aralphaGF[29]*aralphaGF[22];
    aralphaGF[30]=1./4.*aralphaGF[21];
    aralphaGF[31]=17 + aralphaGF[21];
    aralphaGF[31]=aralphaGF[31]*aralphaGF[30];
    aralphaGF[29]=aralphaGF[29] + aralphaGF[31];
    aralphaGF[29]=aralphaGF[18]*aralphaGF[29];
    aralphaGF[25]=aralphaGF[25] + aralphaGF[29];
    aralphaGF[28]=19./2. - aralphaGF[28];
    aralphaGF[28]=aralphaGF[28]*aralphaGF[22];
    aralphaGF[29]=1./2.*aralphaGF[21];
    aralphaGF[31]=9 - 1./6.*aralphaGF[21];
    aralphaGF[31]=aralphaGF[31]*aralphaGF[29];
    aralphaGF[28]=aralphaGF[31] - 4 + 1./3.*aralphaGF[28];
    aralphaGF[28]=aralphaGF[19]*aralphaGF[28];
    aralphaGF[31]=pow(aralphaGF[6],4);
    aralphaGF[32]=aralphaGF[31] - aralphaGF[26];
    aralphaGF[33]= - aralphaGF[19]*aralphaGF[32];
    aralphaGF[26]= - aralphaGF[18]*aralphaGF[26];
    aralphaGF[31]=aralphaGF[17]*aralphaGF[31];
    aralphaGF[26]=aralphaGF[31] + aralphaGF[26] + aralphaGF[33];
    aralphaGF[26]=MMH*aralphaGF[26];
    aralphaGF[31]=aralphaGF[14]*aralphaGF[32];
    aralphaGF[31]=aralphaGF[27] + aralphaGF[31];
    aralphaGF[31]=aralphaGF[31]*pow(MMH,2);
    aralphaGF[26]=aralphaGF[26] + aralphaGF[31];
    aralphaGF[26]=aralphaGF[7]*aralphaGF[26];
    aralphaGF[23]=aralphaGF[17]*aralphaGF[23];
    aralphaGF[23]=1./12.*aralphaGF[26] - 1./4.*aralphaGF[23] + 
      aralphaGF[28] + 1./3.*aralphaGF[25];
    aralphaGF[23]=aralphaGF[7]*aralphaGF[23];
    aralphaGF[25]= - 2 + 3./4.*aralphaGF[22];
    aralphaGF[26]=11*aralphaGF[22];
    aralphaGF[25]=aralphaGF[25]*aralphaGF[26];
    aralphaGF[26]=4 + aralphaGF[30];
    aralphaGF[26]=aralphaGF[26]*aralphaGF[21];
    aralphaGF[26]=1./3.*aralphaGF[26] + 8 + aralphaGF[25];
    aralphaGF[26]=aralphaGF[15]*aralphaGF[26];
    aralphaGF[28]=aralphaGF[17] - aralphaGF[19];
    aralphaGF[31]=3./4.*aralphaGF[20];
    aralphaGF[28]=aralphaGF[31]*aralphaGF[28];
    aralphaGF[31]=aralphaGF[22] - 2;
    aralphaGF[33]= - aralphaGF[14]*aralphaGF[31];
    aralphaGF[28]=aralphaGF[33] + 63./8. + aralphaGF[28];
    aralphaGF[28]=aralphaGF[22]*aralphaGF[28];
    aralphaGF[25]= - 41./3. - aralphaGF[25];
    aralphaGF[25]=aralphaGF[16]*aralphaGF[25];
    aralphaGF[33]= - 1 - aralphaGF[16];
    aralphaGF[33]=aralphaGF[33]*pow(CW,2);
    aralphaGF[21]=aralphaGF[23] + 4*aralphaGF[33] + aralphaGF[26] + 
      aralphaGF[25] + aralphaGF[27] - 1./24.*aralphaGF[21] - 44./3. + 
      aralphaGF[28];
    aralphaGF[21]=aralphaGF[12]*aralphaGF[21];
    aralphaGF[23]=1./2.*aralphaGF[22];
    aralphaGF[25]= - 2./3. + aralphaGF[23];
    aralphaGF[25]=aralphaGF[25]*aralphaGF[22];
    aralphaGF[25]=4./9. + aralphaGF[25];
    aralphaGF[25]=aralphaGF[5]*aralphaGF[25];
    aralphaGF[26]=aralphaGF[31]*aralphaGF[22];
    aralphaGF[27]=20./9. + aralphaGF[26];
    aralphaGF[27]=aralphaGF[11]*aralphaGF[27];
    aralphaGF[28]=4 + aralphaGF[26];
    aralphaGF[31]=1./3.*aralphaGF[1];
    aralphaGF[28]=aralphaGF[28]*aralphaGF[31];
    aralphaGF[25]=aralphaGF[28] + aralphaGF[25] + aralphaGF[27];
    aralphaGF[25]=aralphaGF[2]*aralphaGF[25];
    aralphaGF[27]=aralphaGF[10]*aralphaGF[5];
    aralphaGF[28]=aralphaGF[9]*aralphaGF[5];
    aralphaGF[33]= - MMt*aralphaGF[28];
    aralphaGF[33]= - aralphaGF[27] + aralphaGF[33];
    aralphaGF[32]=aralphaGF[7]*aralphaGF[32]*aralphaGF[33];
    aralphaGF[24]=aralphaGF[24]*aralphaGF[28];
    aralphaGF[24]=aralphaGF[24] + aralphaGF[32];
    aralphaGF[32]=32./3. - 29./4.*aralphaGF[22];
    aralphaGF[30]=1./3.*aralphaGF[32] + aralphaGF[30];
    aralphaGF[30]=aralphaGF[5]*aralphaGF[30];
    aralphaGF[32]= - 8./3. - aralphaGF[23];
    aralphaGF[32]=aralphaGF[32]*aralphaGF[22];
    aralphaGF[32]=32./9. + aralphaGF[32];
    aralphaGF[33]=aralphaGF[8]*aralphaGF[5];
    aralphaGF[32]=aralphaGF[32]*aralphaGF[33];
    aralphaGF[24]=aralphaGF[30] + aralphaGF[32] + 1./2.*aralphaGF[24];
    aralphaGF[24]=MMt*aralphaGF[24];
    aralphaGF[30]=32./3. - 19./2.*aralphaGF[22];
    aralphaGF[29]=1./3.*aralphaGF[30] - aralphaGF[29];
    aralphaGF[27]=aralphaGF[29]*aralphaGF[27];
    aralphaGF[24]=aralphaGF[27] + aralphaGF[24];
    aralphaGF[24]=aralphaGF[7]*aralphaGF[24];
    aralphaGF[27]= - aralphaGF[4]*aralphaGF[11];
    aralphaGF[27]= - aralphaGF[28] + aralphaGF[27];
    aralphaGF[27]=aralphaGF[27]*aralphaGF[26];
    aralphaGF[23]= - 4./3. + aralphaGF[23];
    aralphaGF[22]=aralphaGF[23]*aralphaGF[22];
    aralphaGF[22]=16./9. + aralphaGF[22];
    aralphaGF[22]=aralphaGF[22]*aralphaGF[33];
    aralphaGF[23]= - aralphaGF[4]*aralphaGF[26];
    aralphaGF[23]= - 4./3. + aralphaGF[23];
    aralphaGF[23]=aralphaGF[23]*aralphaGF[31];
    aralphaGF[26]= - aralphaGF[5] - aralphaGF[11];

    alphaGFret = aralphaGF[21] + aralphaGF[22] + aralphaGF[23] + 
      aralphaGF[24] + aralphaGF[25] + 20./27.*aralphaGF[26] + 
      aralphaGF[27];
    return alphaGFret.real();
  }
} // namespace mr
