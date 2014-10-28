#include <alphaGF.hpp>
std::complex<long double>
alphaGF::a10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aralphaGF[29], alphaGFret;

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
   aralphaGF[21]=aralphaGF[18] - aralphaGF[19];
   aralphaGF[22]= - aralphaGF[13] + aralphaGF[14];
   aralphaGF[22]=MMH*aralphaGF[22];
   aralphaGF[23]=aralphaGF[13] - aralphaGF[14];
   aralphaGF[23]=MMH*aralphaGF[23];
   aralphaGF[23]=aralphaGF[23] - aralphaGF[18] + aralphaGF[19];
   aralphaGF[23]=aralphaGF[7]*MMH*aralphaGF[23];
   aralphaGF[22]=1./4.*aralphaGF[23] + 2*aralphaGF[21] + aralphaGF[22];
   aralphaGF[22]=aralphaGF[7]*aralphaGF[22];
   aralphaGF[22]=1./3.*aralphaGF[22] - aralphaGF[14] + 33./4.*
   aralphaGF[15] + aralphaGF[13] - 33./4.*aralphaGF[16];
   aralphaGF[22]=aralphaGF[12]*aralphaGF[22];
   aralphaGF[23]= - aralphaGF[8] + aralphaGF[9];
   aralphaGF[23]=aralphaGF[5]*MMt*aralphaGF[23];
   aralphaGF[24]=MMt*aralphaGF[9];
   aralphaGF[24]=aralphaGF[10] + aralphaGF[24];
   aralphaGF[24]=aralphaGF[7]*aralphaGF[5]*MMt*aralphaGF[24];
   aralphaGF[23]=aralphaGF[23] + aralphaGF[24];
   aralphaGF[23]=aralphaGF[7]*aralphaGF[23];
   aralphaGF[24]= - aralphaGF[11]*aralphaGF[4];
   aralphaGF[25]= - aralphaGF[1]*aralphaGF[4];
   aralphaGF[26]=aralphaGF[11] + 1./3.*aralphaGF[1];
   aralphaGF[26]=aralphaGF[2]*aralphaGF[26];
   aralphaGF[27]=1./2.*aralphaGF[2] + 1./2.*aralphaGF[8] - aralphaGF[9]
   ;
   aralphaGF[27]=aralphaGF[5]*aralphaGF[27];
   aralphaGF[22]=aralphaGF[22] + 1./2.*aralphaGF[23] + aralphaGF[27] + 
   aralphaGF[26] + aralphaGF[24] + 1./3.*aralphaGF[25];
   aralphaGF[23]=pow(aralphaGF[3],2);
   aralphaGF[22]=aralphaGF[23]*aralphaGF[22];
   aralphaGF[24]=aralphaGF[11]*aralphaGF[4];
   aralphaGF[25]=aralphaGF[1]*aralphaGF[4];
   aralphaGF[26]= - aralphaGF[11] - 1./3.*aralphaGF[1];
   aralphaGF[26]=aralphaGF[2]*aralphaGF[26];
   aralphaGF[27]= - 1./3.*aralphaGF[2] - 2./3.*aralphaGF[8] + 
   aralphaGF[9];
   aralphaGF[27]=aralphaGF[5]*aralphaGF[27];
   aralphaGF[24]=aralphaGF[27] + aralphaGF[26] + aralphaGF[24] + 1./3.*
   aralphaGF[25];
   aralphaGF[25]= - aralphaGF[17] - 31./3.*aralphaGF[18];
   aralphaGF[25]=1./2.*aralphaGF[25] + 19./3.*aralphaGF[19];
   aralphaGF[26]= - 1./8. - aralphaGF[14];
   aralphaGF[26]=1./3.*MMH*aralphaGF[26];
   aralphaGF[25]=1./2.*aralphaGF[25] + aralphaGF[26];
   aralphaGF[25]=aralphaGF[7]*aralphaGF[25];
   aralphaGF[27]=aralphaGF[17]*aralphaGF[20];
   aralphaGF[27]=21./2. + aralphaGF[27];
   aralphaGF[28]= - aralphaGF[19]*aralphaGF[20];
   aralphaGF[25]=aralphaGF[25] + 3./4.*aralphaGF[28] + 2*aralphaGF[14]
    - 22*aralphaGF[15] + 3./4.*aralphaGF[27] + 22*aralphaGF[16];
   aralphaGF[25]=aralphaGF[12]*aralphaGF[25];
   aralphaGF[27]= - 29./4. - 8*aralphaGF[8];
   aralphaGF[27]=1./3.*aralphaGF[27] - 1./2.*aralphaGF[9];
   aralphaGF[27]=MMt*aralphaGF[27];
   aralphaGF[27]= - 19./6.*aralphaGF[10] + aralphaGF[27];
   aralphaGF[27]=aralphaGF[7]*aralphaGF[5]*aralphaGF[27];
   aralphaGF[22]=aralphaGF[22] + aralphaGF[25] + 2*aralphaGF[24] + 
   aralphaGF[27];
   aralphaGF[22]=aralphaGF[23]*aralphaGF[22];
   aralphaGF[23]=pow(aralphaGF[6],2);
   aralphaGF[21]=aralphaGF[23]*aralphaGF[21];
   aralphaGF[24]= - aralphaGF[17] + 17./3.*aralphaGF[18];
   aralphaGF[24]=1./2.*aralphaGF[24] + 9*aralphaGF[19];
   aralphaGF[21]=1./12.*aralphaGF[21] + 1./2.*aralphaGF[24] + 
   aralphaGF[26];
   aralphaGF[21]=aralphaGF[23]*aralphaGF[21];
   aralphaGF[24]=MMH*aralphaGF[14];
   aralphaGF[24]=aralphaGF[24] + aralphaGF[17] - aralphaGF[19];
   aralphaGF[25]=pow(aralphaGF[6],4);
   aralphaGF[24]=aralphaGF[7]*aralphaGF[25]*MMH*aralphaGF[24];
   aralphaGF[21]=1./12.*aralphaGF[24] - 4*aralphaGF[19] + aralphaGF[21]
   ;
   aralphaGF[21]=aralphaGF[7]*aralphaGF[21];
   aralphaGF[24]=aralphaGF[23]*aralphaGF[15];
   aralphaGF[24]=1./4.*aralphaGF[24] - 1./8. + 4*aralphaGF[15];
   aralphaGF[24]=aralphaGF[23]*aralphaGF[24];
   aralphaGF[26]=pow(CW,2);
   aralphaGF[27]= - 11./3. - aralphaGF[26];
   aralphaGF[26]= - 41./3. - 4*aralphaGF[26];
   aralphaGF[26]=aralphaGF[16]*aralphaGF[26];
   aralphaGF[21]=aralphaGF[21] + 1./3.*aralphaGF[24] + 8*aralphaGF[15]
    + 4*aralphaGF[27] + aralphaGF[26];
   aralphaGF[21]=aralphaGF[12]*aralphaGF[21];
   aralphaGF[24]= - 5./3.*aralphaGF[11] - aralphaGF[1];
   aralphaGF[26]=5./3.*aralphaGF[11] + aralphaGF[1];
   aralphaGF[26]=aralphaGF[2]*aralphaGF[26];
   aralphaGF[27]=aralphaGF[2] - 5./3. + 4*aralphaGF[8];
   aralphaGF[27]=aralphaGF[5]*aralphaGF[27];
   aralphaGF[24]=1./3.*aralphaGF[27] + 1./3.*aralphaGF[24] + 
   aralphaGF[26];
   aralphaGF[26]=1./2. - aralphaGF[9];
   aralphaGF[26]=MMt*aralphaGF[26];
   aralphaGF[26]= - aralphaGF[10] + aralphaGF[26];
   aralphaGF[23]=aralphaGF[23]*aralphaGF[26];
   aralphaGF[26]=1 + aralphaGF[8];
   aralphaGF[26]=MMt*aralphaGF[26];
   aralphaGF[26]=aralphaGF[10] + aralphaGF[26];
   aralphaGF[23]=32./9.*aralphaGF[26] + 1./2.*aralphaGF[23];
   aralphaGF[23]=aralphaGF[5]*aralphaGF[23];
   aralphaGF[26]= - MMt*aralphaGF[9];
   aralphaGF[26]= - aralphaGF[10] + aralphaGF[26];
   aralphaGF[25]=aralphaGF[7]*aralphaGF[5]*aralphaGF[25]*MMt*
   aralphaGF[26];
   aralphaGF[23]=aralphaGF[23] + 1./2.*aralphaGF[25];
   aralphaGF[23]=aralphaGF[7]*aralphaGF[23];

      alphaGFret = aralphaGF[21] + aralphaGF[22] + aralphaGF[23] + 4./3.
      *aralphaGF[24];
      return alphaGFret;
}
