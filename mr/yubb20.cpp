#include <bb.hpp>
std::complex<long double>
bb::my20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryubb[299], yubbret;

    aryubb[1]=double(nL + nH);
    aryubb[2]=pow(CW,-1);
    aryubb[3]=pow(MMH,-1);
    aryubb[4]=pow(MMZ,-1);
    aryubb[5]=pow(SW,-1);
    aryubb[6]=Tsil::I2(0,0,MMZ,mu2);
    aryubb[7]=Tsil::I2(0,0,MMW,mu2);
    aryubb[8]=Tsil::A(MMH,mu2);
    aryubb[9]=std::real(Tsil::B(0,0,MMZ,mu2));
    aryubb[10]=std::real(Tsil::B(0,0,MMW,mu2));
    aryubb[11]=Tsil::A(MMZ,mu2);
    aryubb[12]=Tsil::A(MMW,mu2);
    aryubb[13]=Tsil::A(MMt,mu2);
    aryubb[14]=Tsil::A(MMb,mu2);
    aryubb[15]=pow(MMb,-1);
    aryubb[16]=Tsil::Aeps(MMZ,mu2);
    aryubb[17]=Tsil::Aeps(MMW,mu2);
    aryubb[18]=Tsil::Aeps(MMb,mu2);
    aryubb[19]=prot0bb0b->Tvxu(0);
    aryubb[20]=double(nL);
    aryubb[21]=double(boson);
    aryubb[22]=det(MMW,MMH,MMt);
    aryubb[23]=Tsil::I2(MMW,MMH,MMt,mu2);
    aryubb[24]=Tsil::Aeps(MMH,mu2);
    aryubb[25]=Tsil::Aeps(MMt,mu2);
    aryubb[26]=det(MMW,MMZ,MMt);
    aryubb[27]=Tsil::I2(MMW,MMZ,MMt,mu2);
    aryubb[28]=Tsil::I2(MMH,MMH,MMH,mu2);
    aryubb[29]=Tsil::I2(MMH,MMZ,MMZ,mu2);
    aryubb[30]=Tsil::I2(MMH,MMW,MMW,mu2);
    aryubb[31]=Tsil::I2(MMH,MMt,MMt,mu2);
    aryubb[32]=Tsil::I2(MMZ,MMt,MMt,mu2);
    aryubb[33]=Tsil::I2(MMW,MMW,MMZ,mu2);
    aryubb[34]=Tsil::I2(0,MMZ,MMH,mu2);
    aryubb[35]=Tsil::I2(0,MMW,MMt,mu2);
    aryubb[36]=Tsil::I2(0,0,MMH,mu2);
    aryubb[37]=Tsil::B(MMH,MMH,MMH,mu2);
    aryubb[38]=Tsil::B(MMH,MMt,MMt,mu2);
    aryubb[39]=Tsil::B(MMZ,MMH,MMZ,mu2);
    aryubb[40]=Tsil::B(MMZ,MMZ,MMH,mu2);
    aryubb[41]=Tsil::B(MMZ,MMt,MMt,mu2);
    aryubb[42]=Tsil::B(MMW,MMH,MMW,mu2);
    aryubb[43]=Tsil::B(MMW,MMZ,MMW,mu2);
    aryubb[44]=Tsil::B(MMW,MMW,MMH,mu2);
    aryubb[45]=Tsil::B(MMW,MMW,MMZ,mu2);
    aryubb[46]=Tsil::B(MMt,MMt,MMH,mu2);
    aryubb[47]=Tsil::B(MMt,MMt,MMZ,mu2);
    aryubb[48]=Tsil::B(0,MMt,MMW,mu2);
    aryubb[49]=std::real(Tsil::B(MMt,MMt,MMH,mu2));
    aryubb[50]=std::real(Tsil::B(0,MMW,MMt,mu2));
    aryubb[51]=pow(MMt,-1);
    aryubb[52]=1/(MMt - MMW);
    aryubb[53]=Tsil::I2(0,0,MMt,mu2);
    aryubb[54]=1/(4*MMt - MMZ);
    aryubb[55]=1/( - 4*MMW + MMZ);
    aryubb[56]=1/( - 4*MMW + MMH);
    aryubb[57]=1/( - 4*MMZ + MMH);
   aryubb[58]= - 9*aryubb[37];
   aryubb[59]= - aryubb[57]*aryubb[29];
   aryubb[60]=1./2.*aryubb[59] - 331./32. + aryubb[58];
   aryubb[60]=3*aryubb[60] - 13./2.*aryubb[50];
   aryubb[61]=aryubb[25] + aryubb[17] + aryubb[24] - aryubb[23];
   aryubb[62]=aryubb[22]*aryubb[61];
   aryubb[63]=aryubb[8]*aryubb[22];
   aryubb[64]=aryubb[62] + aryubb[63];
   aryubb[65]= - MMH*aryubb[22];
   aryubb[66]=5./4.*aryubb[64] + aryubb[65];
   aryubb[66]=MMH*aryubb[66];
   aryubb[67]= - aryubb[54]*aryubb[32];
   aryubb[68]=aryubb[24]*aryubb[57];
   aryubb[69]=7./8.*aryubb[54];
   aryubb[70]=3*aryubb[57] + aryubb[69];
   aryubb[71]=aryubb[16]*aryubb[70];
   aryubb[72]=aryubb[25]*aryubb[54];
   aryubb[73]= - aryubb[30] + aryubb[24];
   aryubb[73]=1./2.*aryubb[73] + aryubb[17];
   aryubb[73]=aryubb[56]*aryubb[73];
   aryubb[74]=1./2.*aryubb[57] + aryubb[56];
   aryubb[74]=aryubb[8]*aryubb[74];
   aryubb[70]=aryubb[11]*aryubb[70];
   aryubb[75]=MMH*aryubb[22];
   aryubb[76]=aryubb[56] + 5./4.*aryubb[75];
   aryubb[76]=aryubb[12]*aryubb[76];
   aryubb[77]= - 9./4.*aryubb[48];
   aryubb[78]= - 3*aryubb[38];
   aryubb[79]=3*aryubb[44];
   aryubb[80]=3./2.*aryubb[40];
   aryubb[81]= - 3*aryubb[39];
   aryubb[82]=33./2.*aryubb[45];
   aryubb[83]= - 1./2.*aryubb[9];
   aryubb[60]=3./2.*aryubb[76] + 3./2.*aryubb[66] + 1./4.*aryubb[70] + 
   3./4.*aryubb[74] + 3./2.*aryubb[73] + aryubb[83] - 5./4.*aryubb[42]
    + aryubb[77] + 7./16.*aryubb[72] + aryubb[81] + 9./4.*aryubb[47] - 
   165./16.*aryubb[43] + 1./4.*aryubb[71] + 3./8.*aryubb[68] + 
   aryubb[80] + aryubb[82] + aryubb[78] + 7./32.*aryubb[67] - 7./16.*
   aryubb[41] + 1./4.*aryubb[60] + aryubb[79];
   aryubb[66]= - aryubb[25] - aryubb[17] - aryubb[24] + aryubb[23];
   aryubb[70]=aryubb[56]*aryubb[66];
   aryubb[71]=21./4. + aryubb[70];
   aryubb[71]=aryubb[56]*aryubb[71];
   aryubb[73]=pow(aryubb[56],2);
   aryubb[74]= - aryubb[8]*aryubb[73];
   aryubb[76]= - aryubb[12]*aryubb[73];
   aryubb[84]=3./4.*aryubb[3];
   aryubb[85]= - aryubb[13]*aryubb[73];
   aryubb[85]=aryubb[85] + aryubb[84] + aryubb[76] + aryubb[71] + 
   aryubb[74];
   aryubb[86]=MMZ*aryubb[73];
   aryubb[85]=1./2.*aryubb[85] + 3*aryubb[86];
   aryubb[85]=MMZ*aryubb[85];
   aryubb[86]=aryubb[8]*aryubb[73];
   aryubb[87]= - 1./4.*aryubb[56] + aryubb[86];
   aryubb[88]= - 3./2.*aryubb[3];
   aryubb[76]=aryubb[88] + 1./2.*aryubb[87] + aryubb[76];
   aryubb[76]=aryubb[13]*aryubb[76];
   aryubb[87]= - 3*aryubb[41];
   aryubb[89]=aryubb[87] + 149./2. - 9*aryubb[50];
   aryubb[90]=33./8.*aryubb[43];
   aryubb[91]= - 1./4.*aryubb[25] - 9./4.*aryubb[17] + 1./4.*aryubb[23]
    + aryubb[30] - 5./4.*aryubb[24];
   aryubb[91]=aryubb[56]*aryubb[91];
   aryubb[92]= - 9./4.*aryubb[56] + aryubb[86];
   aryubb[92]=aryubb[12]*aryubb[92];
   aryubb[93]=1./4.*aryubb[11] + aryubb[12];
   aryubb[93]=aryubb[3]*aryubb[93];
   aryubb[94]= - aryubb[8]*aryubb[56];
   aryubb[76]=3*aryubb[85] + 3./2.*aryubb[76] + 9./4.*aryubb[93] + 3./4.
   *aryubb[92] + 15./16.*aryubb[94] + 3./4.*aryubb[91] + 7./8.*
   aryubb[9] + aryubb[77] + 3./4.*aryubb[39] + 21./8.*aryubb[47] + 
   aryubb[90] - 363./16.*aryubb[45] + 1./8.*aryubb[89] + aryubb[78];
   aryubb[76]=MMZ*aryubb[76];
   aryubb[77]= - 3*aryubb[50];
   aryubb[85]= - 113./4. + aryubb[77];
   aryubb[85]=5*aryubb[85] + aryubb[87];
   aryubb[87]=aryubb[8]*aryubb[56];
   aryubb[89]= - aryubb[12]*aryubb[56];
   aryubb[91]=1./2.*aryubb[9];
   aryubb[92]= - 3*aryubb[48];
   aryubb[85]=3./4.*aryubb[89] + 3./4.*aryubb[87] + aryubb[91] + 
   aryubb[92] + 17./4.*aryubb[47] + aryubb[90] - 33./2.*aryubb[45] + 1./
   16.*aryubb[85] + aryubb[78];
   aryubb[85]=aryubb[12]*aryubb[85];
   aryubb[90]=3*aryubb[50];
   aryubb[93]=aryubb[41] - 1./4. + aryubb[90];
   aryubb[95]=3*aryubb[38];
   aryubb[96]= - 55./4.*aryubb[43];
   aryubb[97]=3*aryubb[48];
   aryubb[93]=aryubb[42] + aryubb[97] - aryubb[39] - 7./2.*aryubb[47]
    + aryubb[96] + 121./4.*aryubb[45] + 1./2.*aryubb[93] + aryubb[95];
   aryubb[98]= - 3*aryubb[12];
   aryubb[99]= - 1./2.*aryubb[11];
   aryubb[100]=aryubb[99] + aryubb[98];
   aryubb[100]=aryubb[3]*aryubb[100];
   aryubb[101]=3*aryubb[89];
   aryubb[102]=aryubb[13]*aryubb[3];
   aryubb[93]=9*aryubb[102] + 9./2.*aryubb[100] + aryubb[101] + 3./2.*
   aryubb[87] + 3*aryubb[93] - 7./2.*aryubb[9];
   aryubb[93]=aryubb[13]*aryubb[93];
   aryubb[100]= - aryubb[41] - 49./2. + aryubb[77];
   aryubb[103]= - 121./4.*aryubb[45];
   aryubb[100]=aryubb[42] + aryubb[92] + aryubb[39] + 7./2.*aryubb[47]
    - 11./4.*aryubb[43] + aryubb[103] + 1./2.*aryubb[100] - 5*
   aryubb[38];
   aryubb[104]=7./2.*aryubb[9];
   aryubb[100]=aryubb[101] + 9./2.*aryubb[87] + 3*aryubb[100] + 
   aryubb[104];
   aryubb[100]=aryubb[12]*aryubb[100];
   aryubb[101]=11*aryubb[33];
   aryubb[105]=1./2.*aryubb[30];
   aryubb[106]= - 3./8.*aryubb[31] + aryubb[101] + aryubb[105];
   aryubb[107]=5 + aryubb[95];
   aryubb[107]=1./4.*aryubb[107] - aryubb[42];
   aryubb[107]=MMH*aryubb[107];
   aryubb[108]=1./2.*aryubb[11];
   aryubb[109]=aryubb[108] + aryubb[12];
   aryubb[110]=aryubb[3]*aryubb[12]*aryubb[109];
   aryubb[111]= - 3./4.*aryubb[8];
   aryubb[93]=1./2.*aryubb[93] + 9./4.*aryubb[110] + 1./2.*aryubb[100]
    + aryubb[107] - 33./4.*aryubb[11] + aryubb[111] + 535./16.*
   aryubb[25] - 823./16.*aryubb[17] - 33./4.*aryubb[16] + 3./8.*
   aryubb[23] - 3./4.*aryubb[24] - 287./16.*aryubb[27] + 3*aryubb[106]
    - 109./16.*aryubb[32];
   aryubb[100]= - 3./2.*aryubb[17];
   aryubb[106]= - 1./2.*aryubb[25];
   aryubb[107]=aryubb[106] + aryubb[100] + 1./2.*aryubb[23] + 
   aryubb[105] - aryubb[24];
   aryubb[107]=aryubb[56]*aryubb[107];
   aryubb[112]=aryubb[3]*aryubb[12];
   aryubb[88]= - aryubb[56] + aryubb[88];
   aryubb[88]=aryubb[13]*aryubb[88];
   aryubb[113]= - 33./4.*aryubb[43];
   aryubb[88]=1./2.*aryubb[88] + 3./4.*aryubb[112] + 3./2.*aryubb[89]
    + aryubb[94] + aryubb[107] + aryubb[42] + aryubb[113] + 309./32. - 
   aryubb[38];
   aryubb[89]=MMZ*aryubb[56];
   aryubb[88]=1./2.*aryubb[88] + 3*aryubb[89];
   aryubb[88]=MMZ*aryubb[88];
   aryubb[88]=1./2.*aryubb[93] + 3*aryubb[88];
   aryubb[88]=MMZ*aryubb[88];
   aryubb[89]=pow(MMH,2);
   aryubb[93]= - aryubb[89]*aryubb[42];
   aryubb[114]= - aryubb[13]*MMH;
   aryubb[115]=1./2.*aryubb[114];
   aryubb[116]=aryubb[12]*MMH;
   aryubb[117]=aryubb[115] + aryubb[93] + aryubb[116];
   aryubb[117]=aryubb[13]*aryubb[117];
   aryubb[118]= - aryubb[30] + aryubb[31];
   aryubb[118]= - aryubb[25] + 1./2.*aryubb[118] + aryubb[17];
   aryubb[118]=aryubb[89]*aryubb[118];
   aryubb[119]=aryubb[89]*aryubb[42];
   aryubb[120]= - aryubb[12]*MMH;
   aryubb[119]=aryubb[119] + 1./2.*aryubb[120];
   aryubb[119]=aryubb[12]*aryubb[119];
   aryubb[117]=aryubb[117] + aryubb[118] + aryubb[119];
   aryubb[118]= - 1./2.*aryubb[30];
   aryubb[119]= - 1./2.*aryubb[17];
   aryubb[121]=1./2.*aryubb[25];
   aryubb[122]=aryubb[121] + aryubb[119] + 3./2.*aryubb[23] + 
   aryubb[118] - aryubb[31];
   aryubb[122]=MMH*aryubb[122];
   aryubb[123]=5*aryubb[33] - aryubb[32];
   aryubb[100]=3./2.*aryubb[25] + aryubb[100] + 1./4.*aryubb[123] - 
   aryubb[27];
   aryubb[113]=aryubb[42] + aryubb[113] - 3./16. - aryubb[38];
   aryubb[113]=aryubb[12]*aryubb[113];
   aryubb[123]= - aryubb[42] + 33./4.*aryubb[43] + 3./16. + aryubb[38];
   aryubb[123]=aryubb[13]*aryubb[123];
   aryubb[100]=aryubb[123] + 11./2.*aryubb[100] + aryubb[113];
   aryubb[100]=MMZ*aryubb[100];
   aryubb[113]=3./4.*aryubb[38];
   aryubb[123]=aryubb[113] - aryubb[42];
   aryubb[123]=MMH*aryubb[123];
   aryubb[124]=aryubb[123] + 5./16.*aryubb[12];
   aryubb[124]=aryubb[12]*aryubb[124];
   aryubb[125]= - 3./4.*aryubb[38];
   aryubb[126]=aryubb[125] + aryubb[42];
   aryubb[126]=MMH*aryubb[126];
   aryubb[127]=7./8.*aryubb[13] + aryubb[126] - 19./16.*aryubb[12];
   aryubb[127]=aryubb[13]*aryubb[127];
   aryubb[100]=3*aryubb[100] + aryubb[127] + 1./2.*aryubb[122] + 
   aryubb[124];
   aryubb[100]=MMZ*aryubb[100];
   aryubb[122]=aryubb[12] - 1./2.*aryubb[13];
   aryubb[122]=aryubb[13]*aryubb[122];
   aryubb[124]=pow(aryubb[12],2);
   aryubb[122]= - 1./2.*aryubb[124] + aryubb[122];
   aryubb[127]=pow(MMZ,2);
   aryubb[122]=aryubb[52]*aryubb[127]*aryubb[122];
   aryubb[100]=9./16.*aryubb[122] + 1./4.*aryubb[117] + aryubb[100];
   aryubb[100]=aryubb[52]*aryubb[100];
   aryubb[117]= - 7*aryubb[30];
   aryubb[128]=aryubb[117] - 15./2.*aryubb[31];
   aryubb[129]=1 + aryubb[42];
   aryubb[129]=MMH*aryubb[129];
   aryubb[128]=aryubb[129] - aryubb[8] - 3./4.*aryubb[25] - 5./4.*
   aryubb[17] + 33./4.*aryubb[23] + 1./2.*aryubb[128] - aryubb[24];
   aryubb[128]=MMH*aryubb[128];
   aryubb[129]=3*aryubb[8];
   aryubb[130]=aryubb[129] + 239./4.*aryubb[11];
   aryubb[131]= - 1 + aryubb[95];
   aryubb[131]=3./4.*aryubb[131] - aryubb[39];
   aryubb[131]=1./2.*aryubb[131] - aryubb[42];
   aryubb[131]=MMH*aryubb[131];
   aryubb[130]= - 353./32.*aryubb[12] + 1./4.*aryubb[130] + aryubb[131]
   ;
   aryubb[130]=aryubb[12]*aryubb[130];
   aryubb[131]= - 1 + aryubb[78];
   aryubb[131]=1./4.*aryubb[131] + aryubb[39];
   aryubb[131]=MMH*aryubb[131];
   aryubb[131]= - 51./8.*aryubb[13] + 79./8.*aryubb[12] - 107./8.*
   aryubb[11] + aryubb[131];
   aryubb[131]=aryubb[13]*aryubb[131];
   aryubb[128]=1./2.*aryubb[131] + 1./4.*aryubb[128] + aryubb[130];
   aryubb[88]=1./2.*aryubb[100] + 1./2.*aryubb[128] + aryubb[88];
   aryubb[88]=aryubb[52]*aryubb[88];
   aryubb[100]= - aryubb[41] - 61./4. - aryubb[50];
   aryubb[100]=1./4.*aryubb[100] + aryubb[38];
   aryubb[128]=3./4.*aryubb[94];
   aryubb[130]=aryubb[12]*aryubb[56];
   aryubb[100]=3./4.*aryubb[130] + aryubb[128] + aryubb[97] + 3*
   aryubb[100] - 13./2.*aryubb[47];
   aryubb[131]= - aryubb[3]*aryubb[12];
   aryubb[100]=9./2.*aryubb[102] + 1./2.*aryubb[100] + 9*aryubb[131];
   aryubb[100]=aryubb[13]*aryubb[100];
   aryubb[101]=aryubb[101] - aryubb[35];
   aryubb[101]= - 95./16.*aryubb[16] + 33./16.*aryubb[23] - 15./16.*
   aryubb[24] - 59./16.*aryubb[27] - 11./4.*aryubb[32] - 17./8.*
   aryubb[31] + 9./8.*aryubb[101] + aryubb[30];
   aryubb[132]=23./4. + aryubb[95];
   aryubb[132]= - aryubb[42] + 1./2.*aryubb[132] - aryubb[39];
   aryubb[132]=MMH*aryubb[132];
   aryubb[76]=1./2.*aryubb[88] + 1./2.*aryubb[76] + 1./4.*aryubb[100]
    + 9./16.*aryubb[110] + 1./2.*aryubb[85] + 1./8.*aryubb[132] - 5./32.
   *aryubb[11] - 9./64.*aryubb[8] + 25./8.*aryubb[25] + 1./4.*
   aryubb[101] - 6*aryubb[17];
   aryubb[76]=aryubb[52]*aryubb[76];
   aryubb[85]=aryubb[63] + aryubb[65];
   aryubb[85]=3./4.*MMH*aryubb[85];
   aryubb[88]=165*aryubb[43] - 165*aryubb[45] + 3*aryubb[41] - 319./4.
    + aryubb[90];
   aryubb[88]=1./2.*aryubb[88] + 5*aryubb[47];
   aryubb[100]= - 3./2.*aryubb[48];
   aryubb[101]= - 3*aryubb[42];
   aryubb[110]= - aryubb[11]*aryubb[54];
   aryubb[132]=3./4.*aryubb[12]*aryubb[75];
   aryubb[88]=aryubb[132] + aryubb[85] + 7./4.*aryubb[110] + aryubb[91]
    + aryubb[101] + aryubb[100] + 1./2.*aryubb[88] + aryubb[81];
   aryubb[133]= - 3*aryubb[3];
   aryubb[69]=aryubb[69] + aryubb[133];
   aryubb[69]=aryubb[13]*aryubb[69];
   aryubb[134]=79./4.*aryubb[12] - aryubb[8] + 97./8.*aryubb[11];
   aryubb[134]=aryubb[3]*aryubb[134];
   aryubb[69]=aryubb[69] + 1./2.*aryubb[88] + aryubb[134];
   aryubb[69]=aryubb[13]*aryubb[69];
   aryubb[88]=75./4.*aryubb[24];
   aryubb[134]= - 5./4.*aryubb[23];
   aryubb[135]=9*aryubb[46];
   aryubb[136]=19*aryubb[47] + 167./8. + aryubb[135];
   aryubb[136]=aryubb[8]*aryubb[136];
   aryubb[136]=aryubb[136] + 91./2.*aryubb[25] + 21./2.*aryubb[17] + 33.
   /4.*aryubb[16] + aryubb[134] + aryubb[88] - 33./4.*aryubb[27] - 
   aryubb[35] - 35./2.*aryubb[31];
   aryubb[137]=1./3.*aryubb[11];
   aryubb[138]= - 7./4.*aryubb[38] - 127./12. + aryubb[135];
   aryubb[139]=aryubb[138] + 1./3.*aryubb[39];
   aryubb[139]=MMH*aryubb[139];
   aryubb[140]=1./4. - aryubb[48];
   aryubb[140]=aryubb[12]*aryubb[140];
   aryubb[141]=aryubb[3]*aryubb[12]*aryubb[8];
   aryubb[136]=1./4.*aryubb[141] + 3./4.*aryubb[140] + 1./8.*
   aryubb[139] + 1./8.*aryubb[136] + aryubb[137];
   aryubb[139]=27*aryubb[37];
   aryubb[140]=3./4.*aryubb[50];
   aryubb[142]=aryubb[140] + 667./48. + aryubb[139];
   aryubb[143]=15./2.*aryubb[38];
   aryubb[142]=19./4.*aryubb[47] + aryubb[80] + aryubb[143] + 1./2.*
   aryubb[142] + aryubb[79];
   aryubb[144]=5*aryubb[8];
   aryubb[145]=aryubb[144] - 43*aryubb[11];
   aryubb[145]=1./2.*aryubb[145] + 3*aryubb[12];
   aryubb[145]=aryubb[3]*aryubb[145];
   aryubb[146]= - aryubb[13]*aryubb[3];
   aryubb[142]=29./4.*aryubb[146] + 1./4.*aryubb[145] + 1./4.*
   aryubb[142] - aryubb[39];
   aryubb[142]=aryubb[13]*aryubb[142];
   aryubb[145]= - 9*aryubb[49];
   aryubb[147]=7./16.*aryubb[50];
   aryubb[148]=19./4.*aryubb[38];
   aryubb[149]=3./8.*aryubb[48] + 19./24.*aryubb[47] + aryubb[148] + 
   aryubb[147] + 1477./384. + aryubb[145];
   aryubb[150]=13*aryubb[35] + 27*aryubb[31];
   aryubb[151]= - 1./2.*aryubb[23];
   aryubb[150]=aryubb[151] + 1./2.*aryubb[150] - 13*aryubb[24];
   aryubb[152]= - 3*aryubb[17];
   aryubb[153]= - 5./2. - 3*aryubb[49];
   aryubb[153]=aryubb[8]*aryubb[153];
   aryubb[150]=3./2.*aryubb[153] - 33./2.*aryubb[25] + 1./2.*
   aryubb[150] + aryubb[152];
   aryubb[153]= - 1./2. - aryubb[48];
   aryubb[153]=aryubb[12]*aryubb[153];
   aryubb[154]= - aryubb[11]*aryubb[48];
   aryubb[153]=3*aryubb[153] + aryubb[150] + 3./4.*aryubb[154];
   aryubb[153]=aryubb[3]*aryubb[153];
   aryubb[149]=1./2.*aryubb[149] + aryubb[153];
   aryubb[153]=aryubb[77] - 155./8. + aryubb[135];
   aryubb[153]= - 19./2.*aryubb[47] + 1./2.*aryubb[153] - 6*aryubb[38];
   aryubb[153]=aryubb[3]*aryubb[153];
   aryubb[155]=pow(aryubb[3],2);
   aryubb[156]= - aryubb[13]*aryubb[155];
   aryubb[153]=aryubb[153] + 9./2.*aryubb[156];
   aryubb[153]=aryubb[13]*aryubb[153];
   aryubb[157]=3 - 1./2.*aryubb[50];
   aryubb[158]=1./2.*aryubb[157] - aryubb[38];
   aryubb[158]=aryubb[3]*aryubb[158];
   aryubb[159]= - aryubb[13]*aryubb[155]*aryubb[49];
   aryubb[158]=aryubb[158] + 6*aryubb[159];
   aryubb[158]=3*MMt*aryubb[158];
   aryubb[149]=aryubb[158] + 1./2.*aryubb[149] + aryubb[153];
   aryubb[149]=MMt*aryubb[149];
   aryubb[136]=aryubb[149] + 1./2.*aryubb[136] + aryubb[142];
   aryubb[136]=MMt*aryubb[136];
   aryubb[142]=31 + aryubb[139];
   aryubb[142]=aryubb[80] + 1./2.*aryubb[142] + aryubb[79];
   aryubb[142]=1./4.*aryubb[142] + aryubb[39];
   aryubb[142]=aryubb[8]*aryubb[142];
   aryubb[117]= - 13./2.*aryubb[31] + aryubb[117] - 5./2.*aryubb[29] - 
   15./2.*aryubb[28] + aryubb[34];
   aryubb[149]=aryubb[42] + 3./2. + aryubb[39];
   aryubb[149]=aryubb[11]*aryubb[149];
   aryubb[117]=1./2.*aryubb[149] + aryubb[142] + 7./8.*aryubb[25] + 9./
   8.*aryubb[17] + aryubb[16] + 19./8.*aryubb[23] + 1./4.*aryubb[117]
    + 7*aryubb[24];
   aryubb[142]= - 1./32.*aryubb[42] + 3./32.*aryubb[40] + 3./16.*
   aryubb[44] - 1 + 27./32.*aryubb[37];
   aryubb[142]=MMH*aryubb[142];
   aryubb[117]=1./4.*aryubb[117] + aryubb[142];
   aryubb[117]=MMH*aryubb[117];
   aryubb[142]=pow(aryubb[8],2);
   aryubb[149]=3./8.*aryubb[142];
   aryubb[153]=13*aryubb[8] + aryubb[99];
   aryubb[153]=aryubb[11]*aryubb[153];
   aryubb[153]=aryubb[149] + aryubb[153];
   aryubb[160]=7*aryubb[8];
   aryubb[161]=aryubb[160] + 41*aryubb[11];
   aryubb[162]=13./16. + aryubb[42];
   aryubb[162]=MMH*aryubb[162];
   aryubb[161]= - 1./8.*aryubb[12] + 1./16.*aryubb[161] + aryubb[162];
   aryubb[161]=aryubb[12]*aryubb[161];
   aryubb[162]=23*aryubb[8] + 7./2.*aryubb[11];
   aryubb[163]=1 + aryubb[78];
   aryubb[164]=1./4.*aryubb[163];
   aryubb[165]=aryubb[164] + aryubb[39];
   aryubb[165]=MMH*aryubb[165];
   aryubb[162]=1./4.*aryubb[162] + aryubb[165];
   aryubb[165]= - 5*aryubb[12];
   aryubb[162]=215./32.*aryubb[13] + 1./2.*aryubb[162] + aryubb[165];
   aryubb[162]=aryubb[13]*aryubb[162];
   aryubb[117]=1./2.*aryubb[162] + 1./2.*aryubb[161] + 1./8.*
   aryubb[153] + aryubb[117];
   aryubb[153]= - 3./2.*aryubb[23];
   aryubb[119]=aryubb[121] + aryubb[119] + aryubb[153] + aryubb[30] + 1.
   /2.*aryubb[31];
   aryubb[119]=aryubb[89]*aryubb[119];
   aryubb[161]= - MMH*aryubb[8];
   aryubb[114]=aryubb[114] + aryubb[161] + 3*aryubb[116];
   aryubb[114]=aryubb[13]*aryubb[114];
   aryubb[162]=MMH*aryubb[42];
   aryubb[166]=3./2.*aryubb[8] + aryubb[162];
   aryubb[166]=MMH*aryubb[166];
   aryubb[166]=1./2.*aryubb[166] + aryubb[120];
   aryubb[166]=aryubb[12]*aryubb[166];
   aryubb[166]=1./4.*aryubb[114] + 1./2.*aryubb[119] + aryubb[166];
   aryubb[166]=aryubb[52]*aryubb[166];
   aryubb[167]= - 3*aryubb[8] - 11./2.*aryubb[11];
   aryubb[168]=aryubb[39] - aryubb[42];
   aryubb[169]=MMH*aryubb[168];
   aryubb[170]= - 3./4.*aryubb[13];
   aryubb[167]=aryubb[170] + aryubb[12] + 1./4.*aryubb[167] + 
   aryubb[169];
   aryubb[167]=aryubb[13]*aryubb[167];
   aryubb[171]=3*aryubb[11] - 1./6.*MMH;
   aryubb[171]=aryubb[12]*aryubb[171];
   aryubb[172]= - aryubb[39] + aryubb[42];
   aryubb[173]=MMH*aryubb[172];
   aryubb[174]=aryubb[11] + aryubb[173];
   aryubb[174]=MMH*aryubb[174];
   aryubb[171]=1./6.*aryubb[174] + aryubb[171];
   aryubb[167]=1./8.*aryubb[171] + aryubb[167];
   aryubb[171]=aryubb[16] + aryubb[35] - aryubb[27];
   aryubb[175]= - aryubb[8]*aryubb[48];
   aryubb[171]=1./2.*aryubb[171] + aryubb[175];
   aryubb[176]= - 1./2. + aryubb[92];
   aryubb[176]=1./8.*aryubb[176] + 3*aryubb[102];
   aryubb[176]=aryubb[13]*aryubb[176];
   aryubb[177]=aryubb[3]*aryubb[48];
   aryubb[178]=aryubb[13]*aryubb[177];
   aryubb[178]= - 1./16.*aryubb[48] + 3*aryubb[178];
   aryubb[178]=MMt*aryubb[178];
   aryubb[171]=aryubb[178] + 3./8.*aryubb[171] + aryubb[176];
   aryubb[171]=MMt*aryubb[171];
   aryubb[167]=1./2.*aryubb[167] + aryubb[171];
   aryubb[167]=MMt*aryubb[167];
   aryubb[171]=aryubb[11]*aryubb[8];
   aryubb[179]=MMH*aryubb[8]*aryubb[172];
   aryubb[179]=aryubb[171] + aryubb[179];
   aryubb[179]=MMH*aryubb[179];
   aryubb[180]=aryubb[174] + aryubb[120];
   aryubb[181]=aryubb[13]*aryubb[180];
   aryubb[161]=aryubb[12]*aryubb[161];
   aryubb[161]=aryubb[181] + aryubb[179] + aryubb[161];
   aryubb[161]=1./16.*aryubb[161];
   aryubb[167]=aryubb[161] + aryubb[167];
   aryubb[167]=aryubb[4]*aryubb[167];
   aryubb[117]=1./2.*aryubb[167] + 1./8.*aryubb[166] + 1./2.*
   aryubb[117] + aryubb[136];
   aryubb[117]=aryubb[4]*aryubb[117];
   aryubb[136]= - 3*aryubb[44];
   aryubb[166]= - 3./2.*aryubb[40];
   aryubb[167]=aryubb[91] + aryubb[101] + aryubb[92] + aryubb[81] + 5./
   2.*aryubb[47] + 165./4.*aryubb[43] + aryubb[166] - 165./4.*
   aryubb[45] - 397./16. + aryubb[136];
   aryubb[167]=aryubb[8]*aryubb[167];
   aryubb[179]=1./4. - 3*aryubb[37];
   aryubb[181]= - 1./2.*aryubb[40];
   aryubb[179]=aryubb[181] + 3./2.*aryubb[179] - aryubb[44];
   aryubb[182]= - aryubb[8]*aryubb[57];
   aryubb[183]=aryubb[11]*aryubb[57];
   aryubb[179]=3./8.*aryubb[183] + 3./8.*aryubb[182] + 3./4.*
   aryubb[179] - aryubb[39];
   aryubb[179]=aryubb[11]*aryubb[179];
   aryubb[66]=aryubb[22]*aryubb[66];
   aryubb[184]= - aryubb[8]*aryubb[22];
   aryubb[185]=aryubb[66] + aryubb[184];
   aryubb[185]=1./2.*aryubb[185] + aryubb[75];
   aryubb[185]=MMH*aryubb[185];
   aryubb[186]=1./2.*aryubb[39] + aryubb[181] + 5./8. - aryubb[44];
   aryubb[186]=3./4.*aryubb[185] + 3./2.*aryubb[186] + aryubb[42];
   aryubb[186]=MMH*aryubb[186];
   aryubb[187]= - 3*aryubb[35] + 13*aryubb[29] + aryubb[34] - 3*
   aryubb[6];
   aryubb[187]=1./2.*aryubb[187] + 11*aryubb[30];
   aryubb[187]=15./32.*aryubb[32] + 1./4.*aryubb[187] - aryubb[31];
   aryubb[167]=1./2.*aryubb[186] + 1./2.*aryubb[179] + 1./4.*
   aryubb[167] + 15./16.*aryubb[25] - 75./32.*aryubb[17] - 47./64.*
   aryubb[16] + 25./32.*aryubb[23] - 81./32.*aryubb[24] + 1./2.*
   aryubb[187] - aryubb[27];
   aryubb[179]= - 27*aryubb[37];
   aryubb[186]= - 3./2.*aryubb[50];
   aryubb[187]=aryubb[186] + 259./8. + aryubb[179];
   aryubb[187]=aryubb[166] + aryubb[78] + 1./2.*aryubb[187] + 
   aryubb[136];
   aryubb[188]=aryubb[184] + 1./2.*aryubb[65];
   aryubb[188]=3./16.*MMH*aryubb[188];
   aryubb[189]= - 15./16.*aryubb[48];
   aryubb[190]= - 1./2.*aryubb[42];
   aryubb[191]=3./16.*aryubb[94];
   aryubb[192]=3./16.*aryubb[130];
   aryubb[187]=aryubb[192] + aryubb[188] + aryubb[191] + aryubb[190] + 
   aryubb[189] + 1./8.*aryubb[187] + aryubb[47];
   aryubb[187]=aryubb[12]*aryubb[187];
   aryubb[193]=111./8. - aryubb[50];
   aryubb[96]=aryubb[96] + aryubb[181] + 55./4.*aryubb[45] - 1./2.*
   aryubb[41] + 1./2.*aryubb[193] - aryubb[44];
   aryubb[193]=3*aryubb[39];
   aryubb[109]=aryubb[3]*aryubb[109];
   aryubb[194]=3*aryubb[42];
   aryubb[96]=9./2.*aryubb[109] + aryubb[83] + aryubb[194] + aryubb[97]
    + aryubb[193] + 3*aryubb[96] - 5./2.*aryubb[47];
   aryubb[96]=aryubb[3]*aryubb[96];
   aryubb[109]= - aryubb[11]*aryubb[26];
   aryubb[195]=3./8.*aryubb[75];
   aryubb[196]= - aryubb[26] - 9./8.*aryubb[22];
   aryubb[196]=aryubb[12]*aryubb[196];
   aryubb[196]=aryubb[196] + aryubb[195] + 3./8.*aryubb[184] + 
   aryubb[109];
   aryubb[96]=1./4.*aryubb[196] + aryubb[96];
   aryubb[96]=aryubb[13]*aryubb[96];
   aryubb[64]=1./2.*aryubb[64] + aryubb[65];
   aryubb[64]=3*MMH*aryubb[64];
   aryubb[196]=7*aryubb[41];
   aryubb[197]=55*aryubb[43] - 55*aryubb[45] + aryubb[196] - 179 + 7*
   aryubb[50];
   aryubb[197]=1./2.*aryubb[197] - 43./3.*aryubb[47];
   aryubb[198]=1./6.*aryubb[9];
   aryubb[197]=aryubb[64] + aryubb[198] - aryubb[42] - 7*aryubb[48] + 1.
   /2.*aryubb[197] - aryubb[39];
   aryubb[199]=aryubb[11]*aryubb[26];
   aryubb[200]=3./16.*aryubb[75] + 3./4.*aryubb[63] + aryubb[199];
   aryubb[200]=aryubb[12]*aryubb[200];
   aryubb[197]=1./8.*aryubb[197] + aryubb[200];
   aryubb[200]= - 9*aryubb[46];
   aryubb[201]=59./8. + aryubb[200];
   aryubb[202]=4*aryubb[47];
   aryubb[203]= - 3./4.*aryubb[48];
   aryubb[204]=aryubb[203] + 1./4.*aryubb[201] + aryubb[202];
   aryubb[204]=aryubb[12]*aryubb[204];
   aryubb[205]= - 3*aryubb[16] + aryubb[151] + 1./2.*aryubb[24] + 13./2.
   *aryubb[35] + 3*aryubb[32];
   aryubb[205]=1./2.*aryubb[205] + aryubb[152];
   aryubb[201]=aryubb[201] + 13*aryubb[47];
   aryubb[201]=aryubb[11]*aryubb[201];
   aryubb[201]=aryubb[204] + 1./8.*aryubb[201] + 1./2.*aryubb[205] - 3*
   aryubb[25];
   aryubb[201]=aryubb[3]*aryubb[201];
   aryubb[204]=aryubb[11]*aryubb[49];
   aryubb[205]=aryubb[12]*aryubb[49];
   aryubb[206]=1./2.*aryubb[204] + aryubb[205];
   aryubb[206]=aryubb[3]*aryubb[206];
   aryubb[207]= - aryubb[41] + 7 - aryubb[50];
   aryubb[207]=1./2.*aryubb[207] - aryubb[48];
   aryubb[206]=1./2.*aryubb[207] + 3*aryubb[206];
   aryubb[206]=MMt*aryubb[3]*aryubb[206];
   aryubb[96]=3*aryubb[206] + aryubb[96] + 1./2.*aryubb[197] + 
   aryubb[201];
   aryubb[96]=MMt*aryubb[96];
   aryubb[105]=aryubb[106] + 1./2.*aryubb[17] + aryubb[153] + 
   aryubb[105] + aryubb[31];
   aryubb[89]=1./2.*aryubb[89]*aryubb[105];
   aryubb[153]= - aryubb[8] + aryubb[11];
   aryubb[197]= - aryubb[39] - aryubb[42];
   aryubb[197]=MMH*aryubb[197];
   aryubb[197]=aryubb[153] + aryubb[197];
   aryubb[197]=MMH*aryubb[197];
   aryubb[197]=aryubb[115] + 1./4.*aryubb[197] + aryubb[116];
   aryubb[197]=aryubb[13]*aryubb[197];
   aryubb[201]=aryubb[39] + aryubb[194];
   aryubb[201]=MMH*aryubb[201];
   aryubb[129]=aryubb[201] + aryubb[129] - aryubb[11];
   aryubb[129]=MMH*aryubb[129];
   aryubb[129]=1./4.*aryubb[129] + aryubb[120];
   aryubb[129]=aryubb[12]*aryubb[129];
   aryubb[129]=aryubb[197] + aryubb[89] + aryubb[129];
   aryubb[129]=aryubb[52]*aryubb[129];
   aryubb[197]= - 5*aryubb[31];
   aryubb[201]=9*aryubb[23] - 1./2.*aryubb[24] - 7./2.*aryubb[30] + 
   aryubb[197];
   aryubb[201]=aryubb[121] + 1./2.*aryubb[201] - aryubb[17];
   aryubb[206]=1./4.*aryubb[8];
   aryubb[207]=1 + aryubb[39];
   aryubb[207]=1./2.*aryubb[207] + aryubb[42];
   aryubb[207]=MMH*aryubb[207];
   aryubb[207]=1./2.*aryubb[207] - 1./4.*aryubb[11] + aryubb[201] + 
   aryubb[206];
   aryubb[207]=MMH*aryubb[207];
   aryubb[208]=13./4.*aryubb[8] + 59*aryubb[11];
   aryubb[209]= - 1 + aryubb[38];
   aryubb[209]=3./4.*aryubb[209] - aryubb[42];
   aryubb[209]=MMH*aryubb[209];
   aryubb[208]= - 65./4.*aryubb[12] + 1./2.*aryubb[208] + aryubb[209];
   aryubb[208]=aryubb[12]*aryubb[208];
   aryubb[210]=1./2.*aryubb[8];
   aryubb[211]=7*aryubb[13] - 105./4.*aryubb[12] + aryubb[210] - 17*
   aryubb[11];
   aryubb[211]=aryubb[13]*aryubb[211];
   aryubb[129]=aryubb[129] + 1./2.*aryubb[211] + aryubb[207] + 
   aryubb[208];
   aryubb[129]=aryubb[52]*aryubb[129];
   aryubb[207]= - aryubb[8] + 41./2.*aryubb[11];
   aryubb[207]=aryubb[11]*aryubb[207];
   aryubb[208]=5./2.*aryubb[11];
   aryubb[211]= - 81./2.*aryubb[12] - aryubb[8] + aryubb[208];
   aryubb[211]=aryubb[12]*aryubb[211];
   aryubb[207]=1./2.*aryubb[207] + aryubb[211];
   aryubb[207]=aryubb[3]*aryubb[207];
   aryubb[69]=aryubb[117] + 1./8.*aryubb[129] + aryubb[96] + 1./4.*
   aryubb[69] + 1./8.*aryubb[207] + 1./2.*aryubb[167] + aryubb[187];
   aryubb[69]=aryubb[4]*aryubb[69];
   aryubb[96]=143./8. + aryubb[200];
   aryubb[117]=3./2.*aryubb[41];
   aryubb[96]=aryubb[100] + 23./4.*aryubb[47] + aryubb[117] + 1./2.*
   aryubb[96] + aryubb[90];
   aryubb[96]=aryubb[3]*aryubb[96];
   aryubb[129]= - aryubb[27] + aryubb[16];
   aryubb[167]=aryubb[25] + aryubb[129] + aryubb[17];
   aryubb[167]=aryubb[26]*aryubb[167];
   aryubb[187]=3./8.*aryubb[63];
   aryubb[207]=aryubb[26] + 3./8.*aryubb[22];
   aryubb[211]=aryubb[12]*aryubb[207];
   aryubb[96]=aryubb[96] + aryubb[211] + 21./8.*aryubb[65] + 
   aryubb[199] + aryubb[187] + aryubb[167] + 3./8.*aryubb[62];
   aryubb[211]= - aryubb[22]*aryubb[56];
   aryubb[212]=aryubb[12]*aryubb[211];
   aryubb[213]=aryubb[22]*aryubb[56];
   aryubb[214]=aryubb[8]*aryubb[213];
   aryubb[207]=3./16.*aryubb[212] + aryubb[207] + 3./16.*aryubb[214];
   aryubb[215]=1./2.*aryubb[40];
   aryubb[216]=aryubb[215] - 1./2. + aryubb[44];
   aryubb[216]=aryubb[155]*aryubb[216];
   aryubb[207]=1./2.*aryubb[207] + 9*aryubb[216];
   aryubb[207]=aryubb[13]*aryubb[207];
   aryubb[216]=MMt*aryubb[155]*aryubb[49];
   aryubb[96]=9*aryubb[216] + 1./2.*aryubb[96] + aryubb[207];
   aryubb[96]=MMt*aryubb[96];
   aryubb[207]= - 3*aryubb[47];
   aryubb[217]=aryubb[207] - 99./2.*aryubb[43] - 17 + 99./2.*aryubb[45]
   ;
   aryubb[218]= - 3./2.*aryubb[9];
   aryubb[217]=aryubb[218] + aryubb[194] + aryubb[97] + 1./2.*
   aryubb[217] + aryubb[81];
   aryubb[217]=aryubb[11]*aryubb[217];
   aryubb[219]= - 99*aryubb[43] + 17 + 99*aryubb[45];
   aryubb[219]=1./2.*aryubb[219];
   aryubb[207]=aryubb[219] + aryubb[207];
   aryubb[207]=aryubb[218] + aryubb[194] + aryubb[97] + 1./2.*
   aryubb[207] + aryubb[81];
   aryubb[207]=aryubb[12]*aryubb[207];
   aryubb[207]=1./2.*aryubb[217] + aryubb[207];
   aryubb[207]=aryubb[3]*aryubb[207];
   aryubb[217]=3*aryubb[47];
   aryubb[219]=aryubb[219] + aryubb[217];
   aryubb[218]=aryubb[218] + aryubb[194] - 9./2.*aryubb[48] + 1./2.*
   aryubb[219] + aryubb[81];
   aryubb[218]=aryubb[12]*aryubb[218];
   aryubb[219]=aryubb[45] - aryubb[43];
   aryubb[220]=aryubb[83] + aryubb[42] + 33./4.*aryubb[219] - 
   aryubb[39];
   aryubb[221]=aryubb[12]*aryubb[220];
   aryubb[222]= - aryubb[45] + aryubb[43];
   aryubb[223]=aryubb[91] - aryubb[42] + 33./4.*aryubb[222] + 
   aryubb[39];
   aryubb[223]=aryubb[13]*aryubb[223];
   aryubb[221]=aryubb[221] + aryubb[223];
   aryubb[221]=MMZ*aryubb[221];
   aryubb[223]= - 17./4.*aryubb[11];
   aryubb[224]=aryubb[223] + aryubb[169];
   aryubb[225]=17./4.*aryubb[12];
   aryubb[226]=aryubb[224] + aryubb[225];
   aryubb[226]=aryubb[12]*aryubb[226];
   aryubb[227]=3./2.*aryubb[13] - 23./4.*aryubb[12] + 17./4.*aryubb[11]
    + aryubb[173];
   aryubb[227]=aryubb[13]*aryubb[227];
   aryubb[221]=3*aryubb[221] + aryubb[226] + aryubb[227];
   aryubb[221]=aryubb[52]*aryubb[221];
   aryubb[220]=MMZ*aryubb[220];
   aryubb[226]=aryubb[97] - 1 - aryubb[47];
   aryubb[226]=aryubb[13]*aryubb[226];
   aryubb[218]=aryubb[221] + 3*aryubb[220] + 3./2.*aryubb[226] + 
   aryubb[224] + aryubb[218];
   aryubb[218]=aryubb[52]*aryubb[218];
   aryubb[220]=11./2.*aryubb[222] + aryubb[47];
   aryubb[221]= - 9*aryubb[48];
   aryubb[220]=3./2.*aryubb[9] - 11*aryubb[42] + aryubb[221] + 9./2.*
   aryubb[220] + 11*aryubb[39];
   aryubb[222]=33./2.*aryubb[219] - aryubb[47];
   aryubb[222]=aryubb[83] + aryubb[42] + aryubb[48] + 1./2.*aryubb[222]
    - aryubb[39];
   aryubb[226]=MMZ*aryubb[3]*aryubb[222];
   aryubb[227]=aryubb[47] - aryubb[48];
   aryubb[228]=MMt*aryubb[3]*aryubb[227];
   aryubb[207]=1./4.*aryubb[218] + 3*aryubb[226] + 3./2.*aryubb[228] + 
   1./8.*aryubb[220] + aryubb[207];
   aryubb[168]=aryubb[8]*aryubb[168];
   aryubb[218]=aryubb[42] - 1./4. - aryubb[39];
   aryubb[218]=aryubb[11]*aryubb[218];
   aryubb[168]=1./8.*aryubb[169] + aryubb[168] + 1./2.*aryubb[218];
   aryubb[168]=MMH*aryubb[168];
   aryubb[218]=17./2.*aryubb[8] + aryubb[11];
   aryubb[220]=aryubb[42] + 1./8. - aryubb[39];
   aryubb[220]=MMH*aryubb[220];
   aryubb[218]= - aryubb[12] + 1./2.*aryubb[218] + aryubb[220];
   aryubb[218]=aryubb[12]*aryubb[218];
   aryubb[220]=aryubb[224] + 11./4.*aryubb[12];
   aryubb[220]=aryubb[13]*aryubb[220];
   aryubb[224]= - 17./2.*aryubb[8] + aryubb[11];
   aryubb[224]=aryubb[11]*aryubb[224];
   aryubb[168]=aryubb[220] + aryubb[218] + 1./2.*aryubb[224] + 
   aryubb[168];
   aryubb[218]=aryubb[8]*aryubb[227];
   aryubb[220]=17./12. + aryubb[92];
   aryubb[220]=aryubb[12]*aryubb[220];
   aryubb[218]=aryubb[220] + 1./3.*aryubb[169] + 3*aryubb[218] - 17./12.
   *aryubb[11];
   aryubb[220]=1./2. + aryubb[47];
   aryubb[165]=31./8.*aryubb[11] + aryubb[165];
   aryubb[165]=aryubb[3]*aryubb[165];
   aryubb[165]=aryubb[165] + aryubb[42] - 3./16.*aryubb[48] + 3./16.*
   aryubb[220] - aryubb[39];
   aryubb[165]=aryubb[13]*aryubb[165];
   aryubb[220]= - aryubb[12]*aryubb[48];
   aryubb[220]=1./2.*aryubb[154] + aryubb[220];
   aryubb[220]=aryubb[3]*aryubb[220];
   aryubb[224]=1./2.*aryubb[47] + aryubb[48];
   aryubb[220]=1./4.*aryubb[224] + 3*aryubb[220];
   aryubb[224]= - aryubb[47] + aryubb[48];
   aryubb[224]=aryubb[13]*aryubb[3]*aryubb[224];
   aryubb[220]=1./2.*aryubb[220] + 3*aryubb[224];
   aryubb[220]=MMt*aryubb[220];
   aryubb[165]=1./2.*aryubb[220] + 1./16.*aryubb[218] + aryubb[165];
   aryubb[165]=MMt*aryubb[165];
   aryubb[170]=aryubb[170] + aryubb[12] + aryubb[169] + aryubb[111] - 
   aryubb[11];
   aryubb[170]=aryubb[13]*aryubb[170];
   aryubb[170]=1./48.*aryubb[180] + aryubb[170];
   aryubb[175]=aryubb[178] + 3./8.*aryubb[175] + aryubb[176];
   aryubb[175]=MMt*aryubb[175];
   aryubb[170]=1./2.*aryubb[170] + aryubb[175];
   aryubb[170]=MMt*aryubb[170];
   aryubb[161]=aryubb[161] + aryubb[170];
   aryubb[161]=aryubb[4]*aryubb[161];
   aryubb[170]=aryubb[12]*aryubb[180];
   aryubb[175]=aryubb[52]*aryubb[170];
   aryubb[161]=1./2.*aryubb[161] + 1./32.*aryubb[175] + 1./8.*
   aryubb[168] + aryubb[165];
   aryubb[161]=aryubb[4]*aryubb[161];
   aryubb[165]=aryubb[8]*aryubb[222];
   aryubb[168]= - aryubb[42] + 25./16. + aryubb[39];
   aryubb[168]=aryubb[11]*aryubb[168];
   aryubb[165]=3./4.*aryubb[173] + 3./2.*aryubb[165] + aryubb[168];
   aryubb[82]= - aryubb[47] - 33./2.*aryubb[43] - 1 + aryubb[82];
   aryubb[168]=3./2.*aryubb[48];
   aryubb[82]=aryubb[83] + aryubb[42] + aryubb[168] + 1./2.*aryubb[82]
    - aryubb[39];
   aryubb[82]=aryubb[13]*aryubb[82];
   aryubb[173]=pow(aryubb[11],2);
   aryubb[175]= - 1./2.*aryubb[173];
   aryubb[176]=aryubb[99] + aryubb[12];
   aryubb[176]=aryubb[12]*aryubb[176];
   aryubb[176]=aryubb[175] + aryubb[176];
   aryubb[176]=aryubb[3]*aryubb[176];
   aryubb[178]= - 25./4. + aryubb[217];
   aryubb[178]= - aryubb[42] - 9./8.*aryubb[48] + 1./8.*aryubb[178] + 
   aryubb[39];
   aryubb[178]=aryubb[12]*aryubb[178];
   aryubb[82]=3./4.*aryubb[82] + 17./4.*aryubb[176] + 1./2.*aryubb[165]
    + aryubb[178];
   aryubb[165]= - aryubb[42] - 1./4. + aryubb[39];
   aryubb[165]=MMH*aryubb[165];
   aryubb[165]=aryubb[225] + aryubb[223] + aryubb[165];
   aryubb[165]=aryubb[12]*aryubb[165];
   aryubb[169]= - aryubb[11] + aryubb[169];
   aryubb[169]=MMH*aryubb[169];
   aryubb[169]=aryubb[169] + aryubb[116];
   aryubb[169]=aryubb[13]*aryubb[169];
   aryubb[169]=aryubb[170] + aryubb[169];
   aryubb[169]=aryubb[52]*aryubb[169];
   aryubb[170]= - aryubb[12] + 1./2.*aryubb[13];
   aryubb[170]=aryubb[13]*aryubb[170];
   aryubb[165]=1./4.*aryubb[169] + 3*aryubb[170] + 1./4.*aryubb[174] + 
   aryubb[165];
   aryubb[165]=aryubb[52]*aryubb[165];
   aryubb[169]= - 1./4.*aryubb[9];
   aryubb[174]=aryubb[169] + 1./2.*aryubb[42] - 1./4.*aryubb[48] - 1./2.
   *aryubb[39] + 33./8.*aryubb[219] - aryubb[47];
   aryubb[176]=aryubb[11]*aryubb[227];
   aryubb[178]=aryubb[12]*aryubb[227];
   aryubb[176]=1./2.*aryubb[176] + aryubb[178];
   aryubb[176]=aryubb[3]*aryubb[176];
   aryubb[174]=1./2.*aryubb[174] + 3*aryubb[176];
   aryubb[176]=33*aryubb[43];
   aryubb[178]=aryubb[176] - 1 - 33*aryubb[45];
   aryubb[178]=1./2.*aryubb[178] + aryubb[47];
   aryubb[91]=aryubb[91] - aryubb[42] - aryubb[48] + 1./2.*aryubb[178]
    + aryubb[39];
   aryubb[91]=aryubb[13]*aryubb[3]*aryubb[91];
   aryubb[178]= - aryubb[3]*aryubb[48];
   aryubb[180]=MMt*aryubb[178];
   aryubb[91]=3./4.*aryubb[180] + 1./4.*aryubb[174] + 3*aryubb[91];
   aryubb[91]=MMt*aryubb[91];
   aryubb[82]=aryubb[161] + 1./8.*aryubb[165] + 1./2.*aryubb[82] + 
   aryubb[91];
   aryubb[82]=aryubb[4]*aryubb[82];
   aryubb[82]=1./2.*aryubb[207] + aryubb[82];
   aryubb[91]=pow(aryubb[5],2);
   aryubb[82]=aryubb[91]*aryubb[82];
   aryubb[161]=aryubb[22]*aryubb[73];
   aryubb[165]=aryubb[8]*aryubb[161];
   aryubb[174]=3*aryubb[165];
   aryubb[180]=7./4.*aryubb[211] + aryubb[174];
   aryubb[207]=aryubb[12]*aryubb[180];
   aryubb[217]=3 + 7./4.*aryubb[70];
   aryubb[217]=aryubb[22]*aryubb[217];
   aryubb[218]=aryubb[8]*aryubb[211];
   aryubb[207]=aryubb[207] + aryubb[217] + 7./4.*aryubb[218];
   aryubb[219]=aryubb[181] + 3./4. - aryubb[44];
   aryubb[219]=aryubb[155]*aryubb[219];
   aryubb[220]= - aryubb[22]*aryubb[73];
   aryubb[222]=aryubb[12]*aryubb[220];
   aryubb[180]=1./2.*aryubb[180] + 3*aryubb[222];
   aryubb[180]=aryubb[13]*aryubb[180];
   aryubb[180]=1./2.*aryubb[180] + 1./4.*aryubb[207] + aryubb[219];
   aryubb[61]=aryubb[56]*aryubb[61];
   aryubb[207]= - 9./4. + aryubb[61];
   aryubb[207]=aryubb[22]*aryubb[56]*aryubb[207];
   aryubb[219]=aryubb[12]*aryubb[161];
   aryubb[223]=aryubb[13]*aryubb[161];
   aryubb[224]=aryubb[223] + aryubb[219] + aryubb[207] + aryubb[165];
   aryubb[224]=MMt*aryubb[224];
   aryubb[226]=13./4. + aryubb[70];
   aryubb[226]=aryubb[22]*aryubb[56]*aryubb[226];
   aryubb[227]=aryubb[8]*aryubb[220];
   aryubb[228]=aryubb[13]*aryubb[220];
   aryubb[226]=aryubb[228] + aryubb[222] + aryubb[226] + aryubb[227];
   aryubb[229]=MMt*aryubb[220];
   aryubb[230]=MMZ*aryubb[161];
   aryubb[226]=9*aryubb[230] + 3./2.*aryubb[226] + aryubb[229];
   aryubb[226]=MMZ*aryubb[226];
   aryubb[180]=3*aryubb[226] + 3*aryubb[180] + 1./2.*aryubb[224];
   aryubb[180]=MMZ*aryubb[180];
   aryubb[224]= - 3*aryubb[57] - 7./16.*aryubb[54];
   aryubb[226]= - 1./2.*aryubb[22] + aryubb[214];
   aryubb[226]=aryubb[12]*aryubb[226];
   aryubb[65]=3./2.*aryubb[226] + aryubb[65] + 3./4.*aryubb[184] + 3./4.
   *aryubb[66] + 1./2.*aryubb[224] - 3*aryubb[56];
   aryubb[66]=1./2. - aryubb[44];
   aryubb[224]=aryubb[66] + aryubb[181];
   aryubb[226]=aryubb[11]*aryubb[224];
   aryubb[224]=aryubb[12]*aryubb[224];
   aryubb[224]=1./2.*aryubb[226] + aryubb[224];
   aryubb[224]=aryubb[3]*aryubb[224];
   aryubb[224]=9*aryubb[224] + 17./4.*aryubb[9] + 9./2.*aryubb[39] + 25.
   /4.*aryubb[47] + 33./2.*aryubb[43] + aryubb[80] - 825./8.*aryubb[45]
    + 661./48. + aryubb[79];
   aryubb[224]=aryubb[3]*aryubb[224];
   aryubb[212]=aryubb[212] - aryubb[22] + aryubb[218];
   aryubb[212]=aryubb[13]*aryubb[212];
   aryubb[65]=9./16.*aryubb[212] + 3./4.*aryubb[65] + aryubb[224];
   aryubb[212]=1./4.*aryubb[213] + aryubb[227];
   aryubb[224]=aryubb[12]*aryubb[212];
   aryubb[212]=1./2.*aryubb[212] + aryubb[219];
   aryubb[212]=aryubb[13]*aryubb[212];
   aryubb[226]= - 3 + 1./2.*aryubb[61];
   aryubb[226]=aryubb[22]*aryubb[226];
   aryubb[212]=1./4.*aryubb[212] + 1./8.*aryubb[224] + 1./32.*
   aryubb[214] - aryubb[26] + 1./16.*aryubb[226];
   aryubb[212]=MMt*aryubb[212];
   aryubb[65]=3./2.*aryubb[180] + 1./2.*aryubb[65] + 3*aryubb[212];
   aryubb[65]=MMZ*aryubb[65];
   aryubb[103]=55./4.*aryubb[43] + aryubb[215] + aryubb[103] - 155./16.
    + aryubb[44];
   aryubb[103]=aryubb[104] + aryubb[101] + aryubb[92] + aryubb[193] + 3
   *aryubb[103] + 11./2.*aryubb[47];
   aryubb[103]=aryubb[11]*aryubb[103];
   aryubb[180]=aryubb[34] - 7*aryubb[6];
   aryubb[212]= - 27./4.*aryubb[29];
   aryubb[215]=99*aryubb[33] + aryubb[180] + aryubb[212];
   aryubb[215]= - 147./2.*aryubb[17] - 147./4.*aryubb[16] + aryubb[23]
    + 69./8.*aryubb[24] - 3*aryubb[32] - 27./4.*aryubb[30] + 1./2.*
   aryubb[215] - 13*aryubb[35];
   aryubb[224]=aryubb[80] + 17./4. + aryubb[79];
   aryubb[224]=aryubb[8]*aryubb[224];
   aryubb[103]=1./2.*aryubb[103] + 3./4.*aryubb[224] + 1./2.*
   aryubb[215] + 9*aryubb[25];
   aryubb[215]= - aryubb[11] - aryubb[12];
   aryubb[215]=aryubb[12]*aryubb[215];
   aryubb[224]= - 1./4.*aryubb[173] + aryubb[215];
   aryubb[224]=aryubb[3]*aryubb[224];
   aryubb[226]= - 311./16. + aryubb[44];
   aryubb[226]=11./4.*aryubb[43] + 1./4.*aryubb[40] + 1./2.*aryubb[226]
    - 11*aryubb[45];
   aryubb[226]=aryubb[9] + 3*aryubb[226] + 2*aryubb[47];
   aryubb[226]=aryubb[12]*aryubb[226];
   aryubb[103]=9./8.*aryubb[224] + 1./2.*aryubb[103] + aryubb[226];
   aryubb[103]=aryubb[3]*aryubb[103];
   aryubb[224]=aryubb[12]*aryubb[22];
   aryubb[224]=9*aryubb[224] + 15./2.*aryubb[75] + 7./4.*aryubb[54] + 9
   *aryubb[184];
   aryubb[90]=aryubb[117] + 115./16. + aryubb[90];
   aryubb[90]=aryubb[3]*aryubb[90];
   aryubb[90]=1./16.*aryubb[224] + aryubb[90];
   aryubb[90]=aryubb[13]*aryubb[90];
   aryubb[60]=aryubb[82] + aryubb[69] + aryubb[76] + aryubb[65] + 
   aryubb[96] + aryubb[90] + 1./4.*aryubb[60] + aryubb[103];
   aryubb[60]=aryubb[91]*aryubb[60];
   aryubb[65]=11./2.*aryubb[41];
   aryubb[69]=47./4.*aryubb[45];
   aryubb[76]=13*aryubb[43];
   aryubb[82]= - 17./6.*aryubb[47];
   aryubb[90]= - 5./6.*aryubb[9];
   aryubb[96]=aryubb[11]*aryubb[54];
   aryubb[103]=aryubb[15]*aryubb[14];
   aryubb[117]=aryubb[132] + aryubb[85] + 1./6.*aryubb[96] + aryubb[90]
    + aryubb[101] + aryubb[100] + aryubb[81] + aryubb[82] + aryubb[76]
    + aryubb[69] + aryubb[65] + aryubb[140] + 3895./72. + aryubb[103];
   aryubb[224]= - aryubb[8] - 125./12.*aryubb[11];
   aryubb[226]=aryubb[224] - 211./12.*aryubb[12];
   aryubb[226]=aryubb[3]*aryubb[226];
   aryubb[117]=1./2.*aryubb[117] + aryubb[226];
   aryubb[226]= - 1./24.*aryubb[54];
   aryubb[229]=aryubb[226] - 15*aryubb[3];
   aryubb[229]=aryubb[13]*aryubb[229];
   aryubb[117]=1./2.*aryubb[117] + aryubb[229];
   aryubb[117]=aryubb[13]*aryubb[117];
   aryubb[229]= - aryubb[15]*aryubb[14];
   aryubb[230]= - 3*aryubb[40];
   aryubb[231]= - 13*aryubb[43];
   aryubb[232]=17./6.*aryubb[47];
   aryubb[233]=5./6.*aryubb[9];
   aryubb[234]=aryubb[11] + aryubb[12];
   aryubb[234]=9./2.*aryubb[3]*aryubb[234];
   aryubb[235]= - 11*aryubb[41];
   aryubb[236]=aryubb[234] + aryubb[233] + aryubb[194] + aryubb[97] + 
   aryubb[193] + aryubb[232] + aryubb[231] + aryubb[230] - 47./4.*
   aryubb[45] + aryubb[235] + aryubb[136] + aryubb[186] - 2971./72. + 
   aryubb[229];
   aryubb[236]=aryubb[3]*aryubb[236];
   aryubb[184]=3*aryubb[184];
   aryubb[237]=3*aryubb[75];
   aryubb[238]= - 9*aryubb[22];
   aryubb[239]=35./3.*aryubb[26] + aryubb[238];
   aryubb[239]=aryubb[12]*aryubb[239];
   aryubb[239]=aryubb[239] + aryubb[237] + aryubb[184] + 29./3.*
   aryubb[199];
   aryubb[236]=1./32.*aryubb[239] + aryubb[236];
   aryubb[236]=aryubb[13]*aryubb[236];
   aryubb[239]=6409./72. + aryubb[103];
   aryubb[240]=7./4.*aryubb[50];
   aryubb[241]=77./6.*aryubb[41];
   aryubb[242]=13./3.*aryubb[43];
   aryubb[243]= - 157./18.*aryubb[48];
   aryubb[244]= - 5./18.*aryubb[9];
   aryubb[239]=aryubb[64] + aryubb[244] - aryubb[42] + aryubb[243] - 
   aryubb[39] + 139./9.*aryubb[47] + aryubb[242] + 47./12.*aryubb[45]
    + aryubb[241] + 1./3.*aryubb[239] + aryubb[240];
   aryubb[245]=aryubb[168] - 139./6.*aryubb[47] - 463./24. + 
   aryubb[200];
   aryubb[245]=aryubb[11]*aryubb[245];
   aryubb[246]= - 1./4.*aryubb[23];
   aryubb[152]=1./2.*aryubb[245] - 25*aryubb[25] + aryubb[152] - 11*
   aryubb[16] + aryubb[246] + 1./4.*aryubb[24] + 13./4.*aryubb[35] + 11
   *aryubb[32];
   aryubb[245]= - 719./24. + aryubb[200];
   aryubb[247]= - 16./3.*aryubb[47];
   aryubb[203]=aryubb[203] + 1./4.*aryubb[245] + aryubb[247];
   aryubb[203]=aryubb[12]*aryubb[203];
   aryubb[203]=1./2.*aryubb[152] + aryubb[203];
   aryubb[203]=aryubb[3]*aryubb[203];
   aryubb[75]=3./32.*aryubb[75] + aryubb[187] + 2./3.*aryubb[109];
   aryubb[75]=aryubb[12]*aryubb[75];
   aryubb[187]= - 5./2.*aryubb[48] + aryubb[235] + 133./3. + 
   aryubb[186];
   aryubb[205]=aryubb[204] + aryubb[205];
   aryubb[205]=aryubb[3]*aryubb[205];
   aryubb[187]=1./2.*aryubb[187] + 9*aryubb[205];
   aryubb[187]=MMt*aryubb[3]*aryubb[187];
   aryubb[75]=aryubb[187] + aryubb[236] + aryubb[203] + 1./16.*
   aryubb[239] + aryubb[75];
   aryubb[75]=MMt*aryubb[75];
   aryubb[69]=aryubb[90] + aryubb[101] + aryubb[92] + aryubb[81] + 
   aryubb[82] + aryubb[76] + aryubb[230] + aryubb[69] + aryubb[136] + 
   2245./72. + aryubb[103];
   aryubb[69]=aryubb[8]*aryubb[69];
   aryubb[203]=aryubb[168] - 7./6.*aryubb[47] + 53./24. + aryubb[135];
   aryubb[203]=aryubb[8]*aryubb[203];
   aryubb[205]= - 73./6.*aryubb[35] - 35*aryubb[31];
   aryubb[236]=49./6.*aryubb[17];
   aryubb[205]=aryubb[203] + 335./6.*aryubb[25] + aryubb[236] + 43./6.*
   aryubb[16] + aryubb[134] + aryubb[88] - 5./6.*aryubb[27] + 1./2.*
   aryubb[205] - 19./3.*aryubb[32];
   aryubb[239]= - aryubb[48] - 53./8. + aryubb[47];
   aryubb[239]=aryubb[11]*aryubb[239];
   aryubb[245]=1./6.*aryubb[42] + aryubb[138] + 1./6.*aryubb[39];
   aryubb[245]=1./2.*MMH*aryubb[245];
   aryubb[248]=13./3. + aryubb[221];
   aryubb[248]=1./4.*aryubb[12]*aryubb[248];
   aryubb[205]=aryubb[141] + aryubb[248] + aryubb[245] + 1./2.*
   aryubb[205] + 1./3.*aryubb[239];
   aryubb[239]=aryubb[140] + 893./144. + aryubb[139];
   aryubb[249]= - 7./24.*aryubb[47];
   aryubb[239]=aryubb[249] + aryubb[80] + aryubb[143] + 1./2.*
   aryubb[239] + aryubb[79];
   aryubb[250]=3./16.*aryubb[48];
   aryubb[251]= - 7./2.*aryubb[12] + 5./2.*aryubb[8] - aryubb[11];
   aryubb[251]=1./2.*aryubb[3]*aryubb[251];
   aryubb[252]=25./6.*aryubb[102];
   aryubb[239]=aryubb[252] + aryubb[251] - aryubb[42] + aryubb[250] + 1.
   /2.*aryubb[239] - aryubb[39];
   aryubb[239]=aryubb[13]*aryubb[239];
   aryubb[253]= - 7./144.*aryubb[47];
   aryubb[254]=5./144.*aryubb[48];
   aryubb[255]=aryubb[254] + aryubb[253] + aryubb[148] + aryubb[147] + 
   559./1152. + aryubb[145];
   aryubb[256]= - 1 + aryubb[100];
   aryubb[256]=aryubb[12]*aryubb[256];
   aryubb[256]=aryubb[150] + 3./2.*aryubb[256];
   aryubb[256]=aryubb[3]*aryubb[256];
   aryubb[77]=aryubb[100] + 7./6.*aryubb[47] - 12*aryubb[38] + 
   aryubb[77] - 17./24. + aryubb[135];
   aryubb[77]=aryubb[3]*aryubb[77];
   aryubb[77]=aryubb[77] + 9*aryubb[156];
   aryubb[77]=aryubb[13]*aryubb[77];
   aryubb[135]= - 2*aryubb[38];
   aryubb[157]=aryubb[157] + aryubb[135];
   aryubb[157]=aryubb[3]*aryubb[157];
   aryubb[157]=aryubb[157] + 12*aryubb[159];
   aryubb[157]=3*MMt*aryubb[157];
   aryubb[159]=aryubb[157] + aryubb[77] + 1./2.*aryubb[255] + 
   aryubb[256];
   aryubb[159]=MMt*aryubb[159];
   aryubb[159]=aryubb[159] + 1./4.*aryubb[205] + aryubb[239];
   aryubb[159]=MMt*aryubb[159];
   aryubb[205]=71./2. + aryubb[139];
   aryubb[205]=aryubb[80] + 1./2.*aryubb[205] + aryubb[79];
   aryubb[205]=aryubb[42] + 1./2.*aryubb[205] + aryubb[39];
   aryubb[205]=aryubb[8]*aryubb[205];
   aryubb[239]=27./8.*aryubb[37];
   aryubb[255]=3./4.*aryubb[44];
   aryubb[257]=3./8.*aryubb[40];
   aryubb[258]=1./108.*aryubb[42] - 103./432.*aryubb[39] + aryubb[257]
    + aryubb[255] - 13./3. + aryubb[239];
   aryubb[258]=MMH*aryubb[258];
   aryubb[259]=499./12. + 35*aryubb[39];
   aryubb[259]=1./8.*aryubb[259] - aryubb[42];
   aryubb[259]=aryubb[11]*aryubb[259];
   aryubb[205]=1./2.*aryubb[258] + 1./9.*aryubb[259] + 1./4.*
   aryubb[205] + 7./16.*aryubb[25] + 9./16.*aryubb[17] + 1./4.*
   aryubb[16] + 19./16.*aryubb[23] + 11./3.*aryubb[24] - 13./16.*
   aryubb[31] - 7./8.*aryubb[30] - 11./48.*aryubb[29] + 5./24.*
   aryubb[34] - 1./3.*aryubb[36] - 15./16.*aryubb[28];
   aryubb[205]=1./2.*MMH*aryubb[205];
   aryubb[114]=1./2.*aryubb[114];
   aryubb[258]=aryubb[144] + 3*aryubb[162];
   aryubb[258]=MMH*aryubb[258];
   aryubb[258]=aryubb[258] + 7*aryubb[120];
   aryubb[258]=aryubb[12]*aryubb[258];
   aryubb[258]=aryubb[114] + aryubb[119] + 1./4.*aryubb[258];
   aryubb[259]=aryubb[52]*aryubb[258];
   aryubb[260]=aryubb[16] + aryubb[34] - aryubb[29];
   aryubb[261]= - aryubb[8]*aryubb[39];
   aryubb[172]=1./18.*aryubb[11]*aryubb[172];
   aryubb[261]=aryubb[172] + 1./3.*aryubb[260] + 1./8.*aryubb[261];
   aryubb[261]=MMH*aryubb[261];
   aryubb[262]= - 1./4.*aryubb[142];
   aryubb[263]= - 5./4.*aryubb[8] + 7./3.*aryubb[11];
   aryubb[263]=aryubb[11]*aryubb[263];
   aryubb[263]=aryubb[262] + 1./3.*aryubb[263];
   aryubb[261]=1./2.*aryubb[263] + aryubb[261];
   aryubb[261]=MMH*aryubb[261];
   aryubb[263]= - MMH*aryubb[39];
   aryubb[263]=aryubb[153] + aryubb[263];
   aryubb[263]=MMH*aryubb[263];
   aryubb[264]= - aryubb[12]*aryubb[11];
   aryubb[265]=1./4.*aryubb[263] + 5*aryubb[264];
   aryubb[266]= - aryubb[16] - aryubb[35] + aryubb[27];
   aryubb[267]=5./4.*aryubb[266] + aryubb[154];
   aryubb[267]=MMt*aryubb[267];
   aryubb[268]=MMH*aryubb[39];
   aryubb[268]=aryubb[268] + aryubb[8] - 11./12.*aryubb[11];
   aryubb[268]=aryubb[13]*aryubb[268];
   aryubb[265]=1./3.*aryubb[267] + 1./12.*aryubb[265] + aryubb[268];
   aryubb[265]=MMt*aryubb[265];
   aryubb[263]=aryubb[13]*aryubb[263];
   aryubb[267]= - aryubb[12]*MMH*aryubb[11];
   aryubb[261]=aryubb[265] + 1./8.*aryubb[263] + aryubb[261] + 1./18.*
   aryubb[267];
   aryubb[263]=aryubb[4]*aryubb[261];
   aryubb[99]= - aryubb[8] + aryubb[99];
   aryubb[163]=aryubb[42] + 1./2.*aryubb[163] + aryubb[39];
   aryubb[163]=MMH*aryubb[163];
   aryubb[265]= - 59./6.*aryubb[12];
   aryubb[99]=15./8.*aryubb[13] + aryubb[265] + 43./6.*aryubb[99] + 
   aryubb[163];
   aryubb[99]=aryubb[13]*aryubb[99];
   aryubb[267]=19*aryubb[8];
   aryubb[268]=aryubb[267] + 1403./9.*aryubb[11];
   aryubb[269]=161./54. + aryubb[194];
   aryubb[269]=MMH*aryubb[269];
   aryubb[268]=151./6.*aryubb[12] + 1./4.*aryubb[268] + aryubb[269];
   aryubb[268]=aryubb[12]*aryubb[268];
   aryubb[270]=3./64.*aryubb[142];
   aryubb[271]=aryubb[8] - 7./12.*aryubb[11];
   aryubb[271]=aryubb[11]*aryubb[271];
   aryubb[99]=1./4.*aryubb[263] + 1./8.*aryubb[259] + aryubb[159] + 1./
   8.*aryubb[99] + 1./8.*aryubb[268] + aryubb[205] + aryubb[270] + 1./3.
   *aryubb[271];
   aryubb[99]=aryubb[4]*aryubb[99];
   aryubb[159]=2551./108. + aryubb[179];
   aryubb[159]= - 11*aryubb[43] + aryubb[166] + 11*aryubb[45] + 1./2.*
   aryubb[159] + aryubb[136];
   aryubb[159]= - 11./3.*aryubb[39] + 1./2.*aryubb[159] - 1./3.*
   aryubb[47];
   aryubb[259]=1./3.*aryubb[48];
   aryubb[263]=5./6.*aryubb[42];
   aryubb[268]= - 1./6.*aryubb[9];
   aryubb[271]=11./8.*aryubb[182];
   aryubb[272]=11./8.*aryubb[183];
   aryubb[159]=aryubb[272] + aryubb[271] + aryubb[268] + aryubb[263] + 
   1./2.*aryubb[159] + aryubb[259];
   aryubb[159]=aryubb[11]*aryubb[159];
   aryubb[273]=aryubb[186] - 57761./216. + aryubb[179];
   aryubb[273]=aryubb[166] + aryubb[78] + 1./2.*aryubb[273] + 
   aryubb[136];
   aryubb[273]=aryubb[192] + aryubb[188] + aryubb[191] + aryubb[190] + 
   aryubb[189] + 1./8.*aryubb[273] - 4./3.*aryubb[47];
   aryubb[273]=aryubb[12]*aryubb[273];
   aryubb[93]=aryubb[115] + 1./4.*aryubb[93] + aryubb[116];
   aryubb[93]=aryubb[13]*aryubb[93];
   aryubb[115]=aryubb[8] + 3./2.*aryubb[162];
   aryubb[115]=MMH*aryubb[115];
   aryubb[115]=1./2.*aryubb[115] + aryubb[120];
   aryubb[115]=aryubb[12]*aryubb[115];
   aryubb[89]=aryubb[93] + aryubb[89] + aryubb[115];
   aryubb[93]=aryubb[52]*aryubb[89];
   aryubb[115]=1./2. + aryubb[42];
   aryubb[115]=MMH*aryubb[115];
   aryubb[115]=aryubb[201] + 1./2.*aryubb[115];
   aryubb[115]=MMH*aryubb[115];
   aryubb[201]=25./24.*aryubb[12] + 1./8.*aryubb[209] + 13./64.*
   aryubb[8] + aryubb[11];
   aryubb[201]=aryubb[12]*aryubb[201];
   aryubb[274]= - 145./2.*aryubb[12] + aryubb[8] - 7*aryubb[11];
   aryubb[274]=1./2.*aryubb[274] + 17*aryubb[13];
   aryubb[274]=aryubb[13]*aryubb[274];
   aryubb[93]=1./8.*aryubb[93] + 1./16.*aryubb[274] + 1./8.*aryubb[115]
    + aryubb[201];
   aryubb[93]=aryubb[52]*aryubb[93];
   aryubb[201]=2*aryubb[18];
   aryubb[274]=17./16.*aryubb[34];
   aryubb[275]=23./16.*aryubb[29];
   aryubb[276]=aryubb[275] - 1./16.*aryubb[6] + aryubb[201] + 
   aryubb[274];
   aryubb[185]=3./2.*aryubb[185] + 59./54.*aryubb[42] + 211./54.*
   aryubb[39] + aryubb[230] + 37./24. + aryubb[136];
   aryubb[185]=1./8.*MMH*aryubb[185];
   aryubb[277]= - 7./3.*aryubb[8] + 3./2.*aryubb[11];
   aryubb[277]=aryubb[11]*aryubb[277];
   aryubb[153]=aryubb[153] + 79./2.*aryubb[12];
   aryubb[153]=aryubb[12]*aryubb[153];
   aryubb[153]=aryubb[277] + aryubb[153];
   aryubb[153]=aryubb[3]*aryubb[153];
   aryubb[278]=11./16.*aryubb[30];
   aryubb[279]= - 1./4.*aryubb[31];
   aryubb[280]=pow(aryubb[14],2);
   aryubb[281]= - aryubb[15]*aryubb[280];
   aryubb[282]= - 319./192.*aryubb[24];
   aryubb[283]=25./64.*aryubb[23];
   aryubb[69]=aryubb[99] + aryubb[93] + aryubb[75] + 1./2.*aryubb[117]
    + 1./8.*aryubb[153] + aryubb[273] + aryubb[185] + 1./2.*aryubb[159]
    + 1./8.*aryubb[69] - 435./32.*aryubb[17] - 11./2.*aryubb[16] + 
   aryubb[283] + aryubb[282] - 179./192.*aryubb[27] - 23./64.*
   aryubb[32] + 1./3.*aryubb[281] + aryubb[279] + aryubb[278] + 169./96.
   *aryubb[35] + 1./3.*aryubb[276] + 11./2.*aryubb[33];
   aryubb[69]=aryubb[4]*aryubb[69];
   aryubb[75]=aryubb[121] + 3./2.*aryubb[17] + aryubb[151] + 
   aryubb[118] + aryubb[24];
   aryubb[75]=aryubb[56]*aryubb[75];
   aryubb[93]= - 21*aryubb[41];
   aryubb[99]= - 2225./12. + aryubb[93];
   aryubb[117]=aryubb[56] + aryubb[3];
   aryubb[117]=aryubb[13]*aryubb[117];
   aryubb[118]= - MMZ*aryubb[56];
   aryubb[75]=27*aryubb[118] + 3./2.*aryubb[117] + 3./2.*aryubb[131] + 
   9./2.*aryubb[130] + 3*aryubb[87] + 3*aryubb[75] + aryubb[101] + 
   aryubb[176] + 1./8.*aryubb[99] + aryubb[95];
   aryubb[75]=MMZ*aryubb[75];
   aryubb[99]=aryubb[221] + 69./2.*aryubb[47] + 47./4.*aryubb[43] - 119
   *aryubb[45] - 9*aryubb[38] + 11./4.*aryubb[41] - 9./2.*aryubb[50] - 
   9659./144. + aryubb[229];
   aryubb[117]=3./2.*aryubb[130];
   aryubb[99]=9./2.*aryubb[146] + 27./4.*aryubb[112] + aryubb[117] + 
   aryubb[128] + 1./2.*aryubb[99] + 5./3.*aryubb[9];
   aryubb[99]=aryubb[13]*aryubb[99];
   aryubb[118]=9./2.*aryubb[50];
   aryubb[130]=9*aryubb[48] - 69./2.*aryubb[47] + 283./4.*aryubb[43] + 
   119*aryubb[45] + 15*aryubb[38] - 53./4.*aryubb[41] + aryubb[118] + 
   13469./144. + aryubb[103];
   aryubb[94]=aryubb[117] + 9./4.*aryubb[94] - 5./3.*aryubb[9] + 1./2.*
   aryubb[130] + aryubb[101];
   aryubb[94]=aryubb[12]*aryubb[94];
   aryubb[117]= - 5 + aryubb[78];
   aryubb[117]=1./4.*aryubb[117] + aryubb[42];
   aryubb[117]=MMH*aryubb[117];
   aryubb[130]= - aryubb[3]*aryubb[124];
   aryubb[94]=aryubb[99] + 9./4.*aryubb[130] + aryubb[94] + aryubb[117]
    + 37./4.*aryubb[11] + 3./4.*aryubb[8] - 1373./24.*aryubb[25] + 559./
   8.*aryubb[17] + 37./4.*aryubb[16] - 3./8.*aryubb[23] + 3./4.*
   aryubb[24] + 259./12.*aryubb[27] + 361./24.*aryubb[32] + 9./8.*
   aryubb[31] - 3./2.*aryubb[30] - 367./8.*aryubb[33] + 11./3.*
   aryubb[35];
   aryubb[75]=1./2.*aryubb[94] + aryubb[75];
   aryubb[75]=MMZ*aryubb[75];
   aryubb[93]=17./2. + aryubb[93];
   aryubb[93]=aryubb[101] + aryubb[176] + 1./8.*aryubb[93] + aryubb[95]
   ;
   aryubb[93]=aryubb[12]*aryubb[93];
   aryubb[94]= - 1./2. + 21*aryubb[41];
   aryubb[78]=aryubb[194] - 33*aryubb[43] + 1./8.*aryubb[94] + 
   aryubb[78];
   aryubb[78]=aryubb[13]*aryubb[78];
   aryubb[94]= - 239./4.*aryubb[25] + 247./4.*aryubb[17] + 135./4.*
   aryubb[27] + 27./2.*aryubb[32] - 189./4.*aryubb[33] - aryubb[35];
   aryubb[78]= - MMZ + aryubb[78] + 1./2.*aryubb[94] + aryubb[93];
   aryubb[78]=MMZ*aryubb[78];
   aryubb[93]=MMH*aryubb[105];
   aryubb[94]=aryubb[126] + 419./48.*aryubb[12];
   aryubb[94]=aryubb[12]*aryubb[94];
   aryubb[95]=49./6.*aryubb[13] + aryubb[123] - 859./48.*aryubb[12];
   aryubb[95]=aryubb[13]*aryubb[95];
   aryubb[93]=aryubb[95] + 1./2.*aryubb[93] + aryubb[94];
   aryubb[78]=1./2.*aryubb[93] + aryubb[78];
   aryubb[78]=MMZ*aryubb[78];
   aryubb[93]=1./2.*aryubb[124] + aryubb[170];
   aryubb[93]=aryubb[52]*aryubb[127]*aryubb[93];
   aryubb[78]=aryubb[78] + 9./16.*aryubb[93];
   aryubb[78]=aryubb[52]*aryubb[78];
   aryubb[93]=17*aryubb[12];
   aryubb[94]= - 57./8.*aryubb[11] + aryubb[93];
   aryubb[94]=aryubb[12]*aryubb[94];
   aryubb[95]=235./3.*aryubb[13] + 49*aryubb[11] - 173*aryubb[12];
   aryubb[95]=aryubb[13]*aryubb[95];
   aryubb[94]=aryubb[94] + 1./8.*aryubb[95];
   aryubb[75]=aryubb[78] + 1./2.*aryubb[94] + aryubb[75];
   aryubb[75]=aryubb[52]*aryubb[75];
   aryubb[78]= - 21./4. + aryubb[61];
   aryubb[78]=aryubb[56]*aryubb[78];
   aryubb[94]=aryubb[12]*aryubb[73];
   aryubb[95]=aryubb[13]*aryubb[73];
   aryubb[99]= - MMZ*aryubb[73];
   aryubb[78]=9*aryubb[99] + aryubb[95] - 1./2.*aryubb[3] + aryubb[94]
    + aryubb[78] + aryubb[86];
   aryubb[78]=MMZ*aryubb[78];
   aryubb[86]=aryubb[118] + 8423./144. + aryubb[103];
   aryubb[95]=1./4.*aryubb[56] + aryubb[74];
   aryubb[94]=3./2.*aryubb[3] + 1./2.*aryubb[95] + aryubb[94];
   aryubb[94]=aryubb[13]*aryubb[94];
   aryubb[95]=1./4.*aryubb[25] + 9./4.*aryubb[17] + aryubb[246] - 
   aryubb[30] + 5./4.*aryubb[24];
   aryubb[95]=aryubb[56]*aryubb[95];
   aryubb[99]=9./4.*aryubb[56] + aryubb[74];
   aryubb[99]=aryubb[12]*aryubb[99];
   aryubb[78]=3./2.*aryubb[78] + 3./4.*aryubb[94] + 3./4.*aryubb[131]
    + 3./8.*aryubb[99] + 15./32.*aryubb[87] + 3./8.*aryubb[95] - 5./12.
   *aryubb[9] - 3./8.*aryubb[42] + 9./8.*aryubb[48] - 69./16.*
   aryubb[47] + 59./16.*aryubb[43] + 119./8.*aryubb[45] + 3./2.*
   aryubb[38] + 1./8.*aryubb[86] - aryubb[41];
   aryubb[78]=MMZ*aryubb[78];
   aryubb[86]=53*aryubb[45];
   aryubb[87]=19./4.*aryubb[43];
   aryubb[94]=aryubb[87] + aryubb[86] - 19./4.*aryubb[41] + 2897./144.
    + aryubb[103];
   aryubb[94]=aryubb[268] + 1./8.*aryubb[94] + aryubb[247];
   aryubb[94]=aryubb[12]*aryubb[94];
   aryubb[95]= - 199./6. + aryubb[41];
   aryubb[99]=7./3.*aryubb[47];
   aryubb[95]=1./16.*aryubb[95] + aryubb[99];
   aryubb[95]=aryubb[13]*aryubb[95];
   aryubb[105]= - 1./8.*aryubb[33] + aryubb[35];
   aryubb[105]= - 35./24.*aryubb[11] - 93./8.*aryubb[25] - 63./8.*
   aryubb[17] + 5./4.*aryubb[16] - 3./8.*aryubb[27] + 11*aryubb[105] + 
   1./2.*aryubb[32];
   aryubb[117]=aryubb[12]*aryubb[11];
   aryubb[118]=aryubb[3]*aryubb[117];
   aryubb[75]=1./2.*aryubb[75] + aryubb[78] + aryubb[95] + 9./32.*
   aryubb[118] + 1./4.*aryubb[105] + aryubb[94];
   aryubb[75]=aryubb[52]*aryubb[75];
   aryubb[78]= - 29./6.*aryubb[47];
   aryubb[94]=9*aryubb[39];
   aryubb[95]=7./6.*aryubb[9];
   aryubb[105]=3*aryubb[40];
   aryubb[123]=aryubb[95] + aryubb[101] + aryubb[92] + aryubb[94] + 
   aryubb[78] + aryubb[76] + aryubb[105] + 179./4.*aryubb[45] + 
   aryubb[79] + 4213./72. + aryubb[103];
   aryubb[123]=aryubb[11]*aryubb[123];
   aryubb[86]=aryubb[87] + aryubb[80] + aryubb[86] + 14171./144. + 
   aryubb[103];
   aryubb[126]= - 2./3.*aryubb[9];
   aryubb[86]=aryubb[126] + 1./2.*aryubb[86] - 8./3.*aryubb[47];
   aryubb[86]=aryubb[12]*aryubb[86];
   aryubb[130]=1./3.*aryubb[180] + aryubb[212];
   aryubb[131]=aryubb[130] - 33*aryubb[33];
   aryubb[151]=77./24.*aryubb[24];
   aryubb[131]=33*aryubb[17] + 93./4.*aryubb[16] + aryubb[151] + 1./2.*
   aryubb[131] + aryubb[32];
   aryubb[153]=aryubb[175] + aryubb[264];
   aryubb[153]=9./8.*aryubb[3]*aryubb[153];
   aryubb[159]=17./2. + 9*aryubb[40];
   aryubb[159]=1./8.*aryubb[8]*aryubb[159];
   aryubb[86]=aryubb[153] + aryubb[86] + 1./4.*aryubb[123] + 
   aryubb[159] + 1./2.*aryubb[131] - aryubb[25];
   aryubb[86]=aryubb[3]*aryubb[86];
   aryubb[123]=9./4. + aryubb[70];
   aryubb[123]=aryubb[22]*aryubb[56]*aryubb[123];
   aryubb[123]=aryubb[228] + aryubb[222] + aryubb[123] + aryubb[227];
   aryubb[123]=MMt*aryubb[123];
   aryubb[131]= - 13./4. + aryubb[61];
   aryubb[131]=aryubb[22]*aryubb[56]*aryubb[131];
   aryubb[131]=aryubb[223] + aryubb[219] + aryubb[131] + aryubb[165];
   aryubb[170]=MMt*aryubb[161];
   aryubb[131]=3./2.*aryubb[131] + aryubb[170];
   aryubb[170]=MMZ*aryubb[220];
   aryubb[131]=1./2.*aryubb[131] + 6*aryubb[170];
   aryubb[131]=MMZ*aryubb[131];
   aryubb[170]=7./4.*aryubb[213] + 3*aryubb[227];
   aryubb[175]=aryubb[12]*aryubb[170];
   aryubb[61]= - 3 + 7./4.*aryubb[61];
   aryubb[61]=aryubb[22]*aryubb[61];
   aryubb[61]=aryubb[175] + aryubb[61] + 7./4.*aryubb[214];
   aryubb[175]=aryubb[181] - 3./4. + 2*aryubb[44];
   aryubb[175]=aryubb[155]*aryubb[175];
   aryubb[170]=1./2.*aryubb[170] + 3*aryubb[219];
   aryubb[170]=aryubb[13]*aryubb[170];
   aryubb[61]=9*aryubb[131] + 1./2.*aryubb[123] + 3./2.*aryubb[170] + 3.
   /4.*aryubb[61] + aryubb[175];
   aryubb[61]=MMZ*aryubb[61];
   aryubb[123]= - 33*aryubb[57];
   aryubb[131]=1./2.*aryubb[22] + aryubb[218];
   aryubb[131]=aryubb[12]*aryubb[131];
   aryubb[62]=9./2.*aryubb[131] + aryubb[237] + 9./4.*aryubb[63] + 9./4.
   *aryubb[62] + 9*aryubb[56] + aryubb[123] + 1./16.*aryubb[54];
   aryubb[131]=aryubb[166] - 1./2. + aryubb[79];
   aryubb[131]=aryubb[12]*aryubb[131];
   aryubb[170]=1 + aryubb[230];
   aryubb[175]=aryubb[11]*aryubb[170];
   aryubb[131]=1./2.*aryubb[175] + aryubb[131];
   aryubb[131]=aryubb[3]*aryubb[131];
   aryubb[176]= - 25./3.*aryubb[47] + 151./4.*aryubb[43] + aryubb[40]
    + 205./2.*aryubb[45] - aryubb[44] - 907./24. + aryubb[103];
   aryubb[126]=3./2.*aryubb[131] + aryubb[126] - 2*aryubb[42] - 2*
   aryubb[48] + 1./2.*aryubb[176] + aryubb[193];
   aryubb[126]=aryubb[3]*aryubb[126];
   aryubb[131]=1./4.*aryubb[211] + aryubb[165];
   aryubb[176]=aryubb[12]*aryubb[131];
   aryubb[131]=1./2.*aryubb[131] + aryubb[222];
   aryubb[131]=aryubb[13]*aryubb[131];
   aryubb[70]=3 + 1./2.*aryubb[70];
   aryubb[70]=aryubb[22]*aryubb[70];
   aryubb[70]=3./4.*aryubb[131] + 3./8.*aryubb[176] + 3./32.*
   aryubb[218] + 7*aryubb[26] + 3./16.*aryubb[70];
   aryubb[70]=MMt*aryubb[70];
   aryubb[131]=aryubb[12]*aryubb[213];
   aryubb[131]=aryubb[131] + aryubb[22] + aryubb[214];
   aryubb[131]=aryubb[13]*aryubb[131];
   aryubb[61]=3*aryubb[61] + aryubb[70] + 9./32.*aryubb[131] + 1./8.*
   aryubb[62] + aryubb[126];
   aryubb[61]=MMZ*aryubb[61];
   aryubb[62]=11./2.*aryubb[59];
   aryubb[70]=1./2.*aryubb[229] + aryubb[62] + 10783./288. + aryubb[58]
   ;
   aryubb[126]=7./24.*aryubb[41];
   aryubb[131]=aryubb[54]*aryubb[32];
   aryubb[70]=497./48.*aryubb[45] + 1./48.*aryubb[131] + aryubb[126] + 
   1./2.*aryubb[70] - aryubb[44];
   aryubb[176]=aryubb[27] - aryubb[16];
   aryubb[181]= - aryubb[25] + aryubb[176] - aryubb[17];
   aryubb[181]=aryubb[26]*aryubb[181];
   aryubb[213]= - aryubb[12]*aryubb[26];
   aryubb[214]=aryubb[213] + aryubb[181] + aryubb[109];
   aryubb[219]= - 3*aryubb[46];
   aryubb[221]=1./2.*aryubb[48];
   aryubb[222]=aryubb[221] - 77./3.*aryubb[47] - aryubb[41] - 145./24.
    + aryubb[219];
   aryubb[222]=aryubb[3]*aryubb[222];
   aryubb[214]=11./6.*aryubb[214] + aryubb[222];
   aryubb[222]= - 1 + aryubb[105];
   aryubb[222]=aryubb[155]*aryubb[222];
   aryubb[223]=3*aryubb[222];
   aryubb[227]= - 11./12.*aryubb[26] + aryubb[223];
   aryubb[227]=aryubb[13]*aryubb[227];
   aryubb[228]=6*aryubb[216];
   aryubb[214]=aryubb[228] + 1./2.*aryubb[214] + aryubb[227];
   aryubb[214]=MMt*aryubb[214];
   aryubb[227]=11*aryubb[57];
   aryubb[226]=aryubb[227] + aryubb[226];
   aryubb[246]=aryubb[16]*aryubb[226];
   aryubb[226]=aryubb[11]*aryubb[226];
   aryubb[247]=3*aryubb[213] - 1./12.*aryubb[54] + 3*aryubb[199];
   aryubb[273]= - 203./24. - aryubb[41];
   aryubb[273]=aryubb[3]*aryubb[273];
   aryubb[247]=1./8.*aryubb[247] + aryubb[273];
   aryubb[247]=aryubb[13]*aryubb[247];
   aryubb[273]=11./16.*aryubb[68];
   aryubb[276]= - 247./144.*aryubb[39];
   aryubb[284]= - aryubb[25]*aryubb[54];
   aryubb[285]=31./144.*aryubb[48];
   aryubb[286]=55./144.*aryubb[42];
   aryubb[287]=aryubb[8]*aryubb[57];
   aryubb[288]=11./16.*aryubb[287];
   aryubb[60]=aryubb[60] + aryubb[69] + aryubb[75] + aryubb[61] + 
   aryubb[214] + aryubb[247] + aryubb[86] + 1./8.*aryubb[226] + 
   aryubb[288] - 79./288.*aryubb[9] + aryubb[286] + aryubb[285] + 1./96.
   *aryubb[284] + aryubb[276] - 679./288.*aryubb[47] - 25./24.*
   aryubb[43] + 1./8.*aryubb[246] + aryubb[273] + 1./4.*aryubb[70] + 
   aryubb[40];
   aryubb[60]=aryubb[91]*aryubb[60];
   aryubb[61]=4./3.*aryubb[41];
   aryubb[69]= - 3./4. + aryubb[61];
   aryubb[70]=pow(CW,2);
   aryubb[69]=aryubb[70]*aryubb[69];
   aryubb[75]=337./128. + aryubb[41];
   aryubb[86]=aryubb[32] + aryubb[27];
   aryubb[86]=1./2.*aryubb[86] - aryubb[16];
   aryubb[86]=aryubb[55]*aryubb[86];
   aryubb[214]= - aryubb[17]*aryubb[55];
   aryubb[226]= - aryubb[25]*aryubb[55];
   aryubb[246]= - aryubb[11]*aryubb[55];
   aryubb[247]= - aryubb[55] - 9*aryubb[56];
   aryubb[247]=aryubb[12]*aryubb[247];
   aryubb[289]= - aryubb[3] - aryubb[55] - aryubb[56];
   aryubb[289]=aryubb[13]*aryubb[289];
   aryubb[290]=1 + 1./2.*aryubb[70];
   aryubb[291]=aryubb[56]*aryubb[290];
   aryubb[291]=1./8.*aryubb[55] + 3*aryubb[291];
   aryubb[291]=MMZ*aryubb[291];
   aryubb[292]= - 3*aryubb[70];
   aryubb[293]= - 41./4. + aryubb[292];
   aryubb[293]=aryubb[43]*aryubb[293];
   aryubb[294]=3./4.*aryubb[42];
   aryubb[69]=3*aryubb[291] + 3./8.*aryubb[289] + 3./8.*aryubb[112] + 1.
   /8.*aryubb[247] + 1./4.*aryubb[246] + aryubb[128] + 3./4.*
   aryubb[107] + aryubb[294] + 3./8.*aryubb[226] + 1./8.*aryubb[214] + 
   1./4.*aryubb[86] + aryubb[293] + aryubb[69] + 5./3.*aryubb[75] + 
   aryubb[125];
   aryubb[69]=MMZ*aryubb[69];
   aryubb[75]=3*aryubb[70];
   aryubb[86]=53./4. + aryubb[75];
   aryubb[86]=aryubb[45]*aryubb[86];
   aryubb[107]= - 3 - 4./3.*aryubb[70];
   aryubb[107]=aryubb[47]*aryubb[107];
   aryubb[128]=aryubb[11]*aryubb[55];
   aryubb[247]= - aryubb[12]*aryubb[55];
   aryubb[289]= - aryubb[13]*aryubb[55];
   aryubb[86]=1./2.*aryubb[289] + 1./2.*aryubb[247] + 3./4.*aryubb[128]
    + aryubb[268] + aryubb[107] + 1./2.*aryubb[86] + 1./6.*aryubb[70]
    - 2./3.*aryubb[41] + 25./3. + 1./8.*aryubb[103];
   aryubb[86]=aryubb[13]*aryubb[86];
   aryubb[107]=1./8.*aryubb[229];
   aryubb[268]= - 53./4. + aryubb[292];
   aryubb[268]=1./2.*aryubb[45]*aryubb[268];
   aryubb[291]=3 + 4./3.*aryubb[70];
   aryubb[291]=aryubb[47]*aryubb[291];
   aryubb[292]=1./4.*aryubb[128] + aryubb[198] + aryubb[291] - 3*
   aryubb[43] + aryubb[268] - 1./6.*aryubb[70] + 2*aryubb[41] - 22./3.
    + aryubb[107];
   aryubb[292]=aryubb[12]*aryubb[292];
   aryubb[295]= - 55./96.*aryubb[25] + 55./96.*aryubb[17] + 119./96.*
   aryubb[27] - 29./32.*aryubb[33] - 1./3.*aryubb[32];
   aryubb[296]=17./3.*aryubb[41];
   aryubb[297]=aryubb[296] + 3*aryubb[43];
   aryubb[297]=aryubb[12]*aryubb[297];
   aryubb[298]=aryubb[296] - aryubb[43];
   aryubb[298]=aryubb[13]*aryubb[298];
   aryubb[295]=1./32.*aryubb[298] + 1./3.*aryubb[295] + 1./32.*
   aryubb[297];
   aryubb[297]=pow(aryubb[2],2);
   aryubb[295]=aryubb[297]*aryubb[295];
   aryubb[298]=71./2.*aryubb[33] - 11*aryubb[35];
   aryubb[298]= - 5./2.*aryubb[11] + 349./12.*aryubb[25] - 277./12.*
   aryubb[17] - 5./2.*aryubb[16] - 83./12.*aryubb[27] + 1./2.*
   aryubb[298] - 25./3.*aryubb[32];
   aryubb[69]=aryubb[69] + aryubb[295] + aryubb[86] + 1./6.*aryubb[298]
    + aryubb[292];
   aryubb[69]=MMZ*aryubb[69];
   aryubb[86]= - 1./4. + aryubb[61];
   aryubb[86]=aryubb[70]*aryubb[86];
   aryubb[292]= - 5./64. + 1./3.*aryubb[41];
   aryubb[86]=aryubb[294] + aryubb[293] + aryubb[86] + 5*aryubb[292] + 
   aryubb[125];
   aryubb[86]=aryubb[12]*aryubb[86];
   aryubb[125]= - 1./4. - 4./3.*aryubb[41];
   aryubb[125]=aryubb[70]*aryubb[125];
   aryubb[75]=41./4. + aryubb[75];
   aryubb[75]=aryubb[43]*aryubb[75];
   aryubb[75]= - 3./4.*aryubb[42] + aryubb[75] + aryubb[125] + 
   aryubb[113] - 7./64. - 5./3.*aryubb[41];
   aryubb[75]=aryubb[13]*aryubb[75];
   aryubb[113]=1./4.*aryubb[35];
   aryubb[125]= - 1./3.*aryubb[27] - 2./3.*aryubb[32] + aryubb[33] + 
   aryubb[113];
   aryubb[125]=aryubb[70]*aryubb[125];
   aryubb[292]=1 + aryubb[70];
   aryubb[292]=aryubb[70]*aryubb[292];
   aryubb[292]=1 + aryubb[292];
   aryubb[292]=MMZ*aryubb[292];
   aryubb[293]= - 53 - 23./2.*aryubb[70];
   aryubb[293]=aryubb[17]*aryubb[293];
   aryubb[294]=25 + 17./4.*aryubb[70];
   aryubb[294]=aryubb[25]*aryubb[294];
   aryubb[75]=1./2.*aryubb[292] + aryubb[75] + aryubb[86] + 1./3.*
   aryubb[294] + 1./6.*aryubb[293] + aryubb[125] - 49./12.*aryubb[27]
    - 9./4.*aryubb[32] + 19./3.*aryubb[33] + aryubb[113];
   aryubb[75]=MMZ*aryubb[75];
   aryubb[86]= - 2*aryubb[13];
   aryubb[113]=aryubb[225] + aryubb[86];
   aryubb[113]=aryubb[13]*aryubb[113];
   aryubb[75]=aryubb[75] - 2*aryubb[124] + aryubb[113];
   aryubb[75]=MMZ*aryubb[75];
   aryubb[75]=aryubb[75] + 9./64.*aryubb[122];
   aryubb[75]=aryubb[52]*aryubb[75];
   aryubb[98]=aryubb[137] + aryubb[98];
   aryubb[98]=aryubb[12]*aryubb[98];
   aryubb[113]=23./12.*aryubb[11] + 5*aryubb[12];
   aryubb[113]=aryubb[12]*aryubb[113];
   aryubb[122]= - 19*aryubb[13] + 97*aryubb[11] - 161*aryubb[12];
   aryubb[122]=aryubb[13]*aryubb[122];
   aryubb[113]=aryubb[113] + 1./12.*aryubb[122];
   aryubb[113]=aryubb[297]*aryubb[113];
   aryubb[93]= - aryubb[11] + aryubb[93];
   aryubb[93]=1./2.*aryubb[93] - 4*aryubb[13];
   aryubb[93]=aryubb[13]*aryubb[93];
   aryubb[69]=aryubb[75] + aryubb[69] + 1./24.*aryubb[113] + 1./2.*
   aryubb[98] + 1./3.*aryubb[93];
   aryubb[69]=aryubb[52]*aryubb[69];
   aryubb[75]= - 3*aryubb[73];
   aryubb[93]=pow(aryubb[55],2);
   aryubb[98]= - aryubb[93] + aryubb[75];
   aryubb[98]=aryubb[12]*aryubb[98];
   aryubb[75]= - 5*aryubb[93] + aryubb[75];
   aryubb[75]=aryubb[13]*aryubb[75];
   aryubb[113]= - 3./2.*aryubb[16] + aryubb[32] + 1./2.*aryubb[27];
   aryubb[113]=aryubb[55]*aryubb[113];
   aryubb[113]=13./8. + aryubb[113];
   aryubb[113]=aryubb[55]*aryubb[113];
   aryubb[122]= - aryubb[17]*aryubb[93];
   aryubb[124]= - aryubb[25]*aryubb[93];
   aryubb[125]= - aryubb[11]*aryubb[93];
   aryubb[71]=1./2.*aryubb[75] + aryubb[84] + 1./2.*aryubb[98] + 3./2.*
   aryubb[125] + 3./2.*aryubb[74] + 3./2.*aryubb[71] + 5./2.*
   aryubb[124] + aryubb[113] + 1./2.*aryubb[122];
   aryubb[74]=aryubb[73]*aryubb[290];
   aryubb[74]=1./8.*aryubb[93] + aryubb[74];
   aryubb[74]=MMZ*aryubb[74];
   aryubb[71]=1./2.*aryubb[71] + 9*aryubb[74];
   aryubb[71]=MMZ*aryubb[71];
   aryubb[74]=aryubb[11]*aryubb[93];
   aryubb[75]= - 7./4.*aryubb[55] + 5*aryubb[74];
   aryubb[84]= - aryubb[12]*aryubb[93];
   aryubb[98]= - aryubb[13]*aryubb[93];
   aryubb[75]=2*aryubb[98] + 1./2.*aryubb[75] + aryubb[84];
   aryubb[75]=aryubb[13]*aryubb[75];
   aryubb[84]=3*aryubb[112] + 55./18. + aryubb[43];
   aryubb[84]=aryubb[297]*aryubb[84];
   aryubb[98]= - 13./8.*aryubb[16] + aryubb[32] + 5./8.*aryubb[27];
   aryubb[98]=aryubb[55]*aryubb[98];
   aryubb[74]= - 5./12.*aryubb[55] + aryubb[74];
   aryubb[74]=aryubb[12]*aryubb[74];
   aryubb[61]=aryubb[71] + 1./16.*aryubb[84] + aryubb[75] + 1./2.*
   aryubb[74] + 13./24.*aryubb[246] + aryubb[198] + 7./8.*aryubb[226]
    + 5./24.*aryubb[214] + 1./3.*aryubb[98] + aryubb[291] - 3./2.*
   aryubb[43] + aryubb[268] - 29./12.*aryubb[70] + aryubb[61] - 98./9.
    + aryubb[107];
   aryubb[61]=MMZ*aryubb[61];
   aryubb[71]= - 13./8. - 4*aryubb[47];
   aryubb[71]=4./3.*aryubb[289] + 5./6.*aryubb[247] + 1./3.*aryubb[71]
    + 3./2.*aryubb[128];
   aryubb[71]=aryubb[13]*aryubb[71];
   aryubb[74]=59./6.*aryubb[25] - 59./6.*aryubb[17] + 47./6.*aryubb[27]
    + aryubb[33] - 53./6.*aryubb[32];
   aryubb[75]=aryubb[12]*aryubb[43];
   aryubb[84]= - aryubb[13]*aryubb[41];
   aryubb[74]=17./6.*aryubb[84] + 1./3.*aryubb[74] + 1./2.*aryubb[75];
   aryubb[74]=aryubb[297]*aryubb[74];
   aryubb[75]= - 403./12. + aryubb[196];
   aryubb[75]=1./3.*aryubb[75] + 19*aryubb[43];
   aryubb[75]=aryubb[12]*aryubb[75];
   aryubb[84]=1./9. + aryubb[41];
   aryubb[84]=aryubb[13]*aryubb[84];
   aryubb[74]=1./2.*aryubb[74] + 1./2.*aryubb[84] + 9./4.*aryubb[118]
    + 1./4.*aryubb[75] - 23./9.*aryubb[11] - 49./18.*aryubb[25] - 1./3.
   *aryubb[17] - 55./36.*aryubb[16] - 1./6.*aryubb[27] + 1./4.*
   aryubb[33] + 13./9.*aryubb[32];
   aryubb[74]=aryubb[297]*aryubb[74];
   aryubb[75]= - 3*aryubb[45];
   aryubb[84]=59./12. + aryubb[75];
   aryubb[84]=2./3.*aryubb[128] + 1./2.*aryubb[84] + 8./3.*aryubb[47];
   aryubb[84]=aryubb[12]*aryubb[84];
   aryubb[61]=aryubb[69] + aryubb[61] + 1./8.*aryubb[74] + aryubb[71]
    + aryubb[84] + 1./8.*aryubb[11] - 13./24.*aryubb[25] + 67./24.*
   aryubb[17] + 1./24.*aryubb[16] + 7./24.*aryubb[27] + 2./3.*
   aryubb[32] - aryubb[33] - 13./12.*aryubb[35];
   aryubb[61]=aryubb[52]*aryubb[61];
   aryubb[69]=aryubb[95] + aryubb[101] + aryubb[92] + aryubb[94] + 
   aryubb[78] + aryubb[76] + aryubb[105] + 35./4.*aryubb[45] + 
   aryubb[79] + 2629./72. + aryubb[103];
   aryubb[69]=aryubb[11]*aryubb[69];
   aryubb[71]=aryubb[130] - 9*aryubb[33];
   aryubb[71]=9*aryubb[17] + 45./4.*aryubb[16] + aryubb[151] + 1./2.*
   aryubb[71] + aryubb[32];
   aryubb[74]=19./2.*aryubb[43] + 353./24. + aryubb[105];
   aryubb[74]=aryubb[12]*aryubb[74];
   aryubb[69]=aryubb[153] + 1./4.*aryubb[74] + 1./4.*aryubb[69] + 
   aryubb[159] + 1./2.*aryubb[71] - aryubb[25];
   aryubb[69]=aryubb[3]*aryubb[69];
   aryubb[71]= - 17./2.*aryubb[43];
   aryubb[74]=aryubb[71] + aryubb[105] + 281./24. + aryubb[45];
   aryubb[74]=1./2.*aryubb[74] + 17./3.*aryubb[47];
   aryubb[74]=aryubb[233] + 1./2.*aryubb[74] + aryubb[193];
   aryubb[74]=aryubb[11]*aryubb[74];
   aryubb[78]= - aryubb[33] + 5./9.*aryubb[180] + aryubb[212];
   aryubb[78]=aryubb[17] + 175./12.*aryubb[16] + 223./72.*aryubb[24] + 
   1./2.*aryubb[78] - 17./3.*aryubb[32];
   aryubb[84]= - 17./3. + 1./2.*aryubb[43];
   aryubb[84]=aryubb[12]*aryubb[84];
   aryubb[94]= - aryubb[3]*aryubb[173];
   aryubb[74]=9./16.*aryubb[94] + 1./2.*aryubb[84] + aryubb[74] + 
   aryubb[159] + 1./2.*aryubb[78] + 17./3.*aryubb[25];
   aryubb[74]=aryubb[3]*aryubb[74];
   aryubb[59]=25./6.*aryubb[59] - 33791./7776. + aryubb[58];
   aryubb[59]= - 55./216.*aryubb[45] + 125./144.*aryubb[131] - 119./72.
   *aryubb[41] + 1./2.*aryubb[59] - aryubb[44];
   aryubb[78]=aryubb[57] - 5./24.*aryubb[54];
   aryubb[84]=aryubb[16]*aryubb[78];
   aryubb[78]=aryubb[11]*aryubb[78];
   aryubb[94]=431./24. + 17*aryubb[41];
   aryubb[94]=aryubb[3]*aryubb[94];
   aryubb[94]= - 125./96.*aryubb[54] + aryubb[94];
   aryubb[94]=aryubb[13]*aryubb[94];
   aryubb[95]=1./3.*aryubb[42];
   aryubb[59]=1./3.*aryubb[94] + aryubb[74] + 25./24.*aryubb[78] + 25./
   48.*aryubb[287] - 275./1296.*aryubb[9] + aryubb[95] + 125./288.*
   aryubb[284] - 139./72.*aryubb[39] - 935./1296.*aryubb[47] + 527./864.
   *aryubb[43] + 25./24.*aryubb[84] + 25./48.*aryubb[68] + 1./4.*
   aryubb[59] + aryubb[40];
   aryubb[68]=aryubb[48] + 49./18.*aryubb[47] + aryubb[296] - 275./72.
    + aryubb[219];
   aryubb[68]=aryubb[3]*aryubb[68];
   aryubb[74]=aryubb[13]*aryubb[222];
   aryubb[68]=1./2.*aryubb[68] + 3*aryubb[74];
   aryubb[68]=1./2.*aryubb[68] + 3*aryubb[216];
   aryubb[68]=MMt*aryubb[68];
   aryubb[74]= - 1./3. - 1./2.*aryubb[43];
   aryubb[74]=aryubb[11]*aryubb[74];
   aryubb[74]=aryubb[74] + 1./3.*aryubb[12];
   aryubb[74]=aryubb[3]*aryubb[74];
   aryubb[74]=31./216.*aryubb[43] + aryubb[74];
   aryubb[74]=aryubb[297]*aryubb[74];
   aryubb[59]=1./8.*aryubb[74] + 1./2.*aryubb[59] + aryubb[68];
   aryubb[59]=aryubb[297]*aryubb[59];
   aryubb[58]=31./54.*aryubb[229] + aryubb[62] - 187739./7776. + 
   aryubb[58];
   aryubb[58]= - 983./432.*aryubb[45] + 125./144.*aryubb[67] + 
   aryubb[126] + 1./2.*aryubb[58] - aryubb[44];
   aryubb[62]=aryubb[167] + aryubb[199];
   aryubb[68]=aryubb[12]*aryubb[26];
   aryubb[74]=aryubb[62] + aryubb[68];
   aryubb[78]=aryubb[221] - 119./9.*aryubb[47] - aryubb[41] + 461./72.
    + aryubb[219];
   aryubb[78]=aryubb[3]*aryubb[78];
   aryubb[74]=19./48.*aryubb[74] + aryubb[78];
   aryubb[78]=19./96.*aryubb[26] + aryubb[223];
   aryubb[78]=aryubb[13]*aryubb[78];
   aryubb[74]=aryubb[228] + 1./2.*aryubb[74] + aryubb[78];
   aryubb[74]=MMt*aryubb[74];
   aryubb[78]=125./3.*aryubb[54] + 19*aryubb[199];
   aryubb[78]=1./3.*aryubb[78] + 27*aryubb[213];
   aryubb[84]= - 161./72. - aryubb[41];
   aryubb[84]=aryubb[3]*aryubb[84];
   aryubb[78]=1./32.*aryubb[78] + aryubb[84];
   aryubb[78]=aryubb[13]*aryubb[78];
   aryubb[84]=aryubb[227] + 125./72.*aryubb[54];
   aryubb[94]=aryubb[16]*aryubb[84];
   aryubb[84]=aryubb[11]*aryubb[84];
   aryubb[98]=aryubb[12]*aryubb[199];
   aryubb[58]=aryubb[59] + aryubb[74] + aryubb[78] + aryubb[69] + 31./
   48.*aryubb[98] + 1./8.*aryubb[84] + aryubb[288] - 295./2592.*
   aryubb[9] + aryubb[286] + aryubb[285] + 125./288.*aryubb[72] + 
   aryubb[276] + 449./2592.*aryubb[47] - 67./72.*aryubb[43] + 1./8.*
   aryubb[94] + aryubb[273] + 1./4.*aryubb[58] + aryubb[40];
   aryubb[58]=aryubb[297]*aryubb[58];
   aryubb[59]= - 1./4.*aryubb[45];
   aryubb[65]=aryubb[132] + aryubb[85] + 125./18.*aryubb[110] + 
   aryubb[90] + aryubb[101] + aryubb[100] + aryubb[81] + aryubb[82] + 
   aryubb[76] + aryubb[59] + aryubb[65] + aryubb[140] + 21407./648. + 
   aryubb[103];
   aryubb[69]=aryubb[224] + 15./4.*aryubb[12];
   aryubb[69]=aryubb[3]*aryubb[69];
   aryubb[65]=1./2.*aryubb[65] + aryubb[69];
   aryubb[69]=25./72.*aryubb[54] + aryubb[133];
   aryubb[69]=aryubb[13]*aryubb[69];
   aryubb[65]=1./2.*aryubb[65] + 5*aryubb[69];
   aryubb[65]=aryubb[13]*aryubb[65];
   aryubb[69]=aryubb[234] + aryubb[233] + aryubb[194] + aryubb[97] + 
   aryubb[193] + aryubb[232] + aryubb[231] + aryubb[230] + 1./4.*
   aryubb[45] + aryubb[235] + aryubb[136] + aryubb[186] - 2107./72. + 
   aryubb[229];
   aryubb[69]=aryubb[3]*aryubb[69];
   aryubb[72]=43./3.*aryubb[26] + aryubb[238];
   aryubb[72]=aryubb[12]*aryubb[72];
   aryubb[72]=aryubb[72] + aryubb[237] + aryubb[184] + 19./3.*
   aryubb[109];
   aryubb[69]=1./32.*aryubb[72] + aryubb[69];
   aryubb[69]=aryubb[13]*aryubb[69];
   aryubb[72]=15115./216. + aryubb[103];
   aryubb[64]=aryubb[64] + aryubb[244] - aryubb[42] + aryubb[243] - 
   aryubb[39] + 1763./81.*aryubb[47] + aryubb[242] - 1./12.*aryubb[45]
    + aryubb[241] + 1./3.*aryubb[72] + aryubb[240];
   aryubb[63]=aryubb[195] + 3./2.*aryubb[63] + aryubb[109];
   aryubb[63]=aryubb[12]*aryubb[63];
   aryubb[63]=1./4.*aryubb[64] + aryubb[63];
   aryubb[64]= - aryubb[48] - 23./8. + aryubb[219];
   aryubb[64]=aryubb[12]*aryubb[64];
   aryubb[64]=aryubb[152] + 3./2.*aryubb[64];
   aryubb[64]=aryubb[3]*aryubb[64];
   aryubb[63]=1./2.*aryubb[63] + aryubb[64];
   aryubb[63]=aryubb[187] + 1./2.*aryubb[63] + aryubb[69];
   aryubb[63]=MMt*aryubb[63];
   aryubb[59]=aryubb[90] + aryubb[101] + aryubb[92] + aryubb[81] + 
   aryubb[82] + aryubb[76] + aryubb[230] + aryubb[59] + aryubb[136] + 
   1381./72. + aryubb[103];
   aryubb[59]=aryubb[8]*aryubb[59];
   aryubb[64]=10303./216. + aryubb[179];
   aryubb[64]=aryubb[192] + aryubb[188] + aryubb[191] + aryubb[190] + 
   aryubb[189] - 3./16.*aryubb[40] - 3./8.*aryubb[38] - 3./8.*
   aryubb[44] - 3./32.*aryubb[50] + 1./16.*aryubb[64] + 7./3.*
   aryubb[229];
   aryubb[64]=aryubb[12]*aryubb[64];
   aryubb[69]= - 79./24. + aryubb[200];
   aryubb[69]=aryubb[168] + 1./2.*aryubb[69] + aryubb[99];
   aryubb[69]=aryubb[11]*aryubb[69];
   aryubb[72]= - aryubb[32] + aryubb[16];
   aryubb[72]=1./2.*aryubb[72] + aryubb[25];
   aryubb[69]=7./3.*aryubb[72] + 1./2.*aryubb[69];
   aryubb[69]=aryubb[3]*aryubb[69];
   aryubb[72]=47./8. + aryubb[196];
   aryubb[74]=aryubb[3]*aryubb[11];
   aryubb[72]=9./2.*aryubb[74] + 17./2.*aryubb[43] + 1./3.*aryubb[72]
    + aryubb[230];
   aryubb[72]=aryubb[13]*aryubb[3]*aryubb[72];
   aryubb[76]= - 7*aryubb[41];
   aryubb[78]= - 97./18. + aryubb[76];
   aryubb[81]= - 17*aryubb[43];
   aryubb[78]=7./3.*aryubb[78] + aryubb[81];
   aryubb[78]= - 31./3.*aryubb[48] + 1./2.*aryubb[78] - 385./27.*
   aryubb[47];
   aryubb[82]=11 + aryubb[196];
   aryubb[82]=1./3.*aryubb[82] + aryubb[48];
   aryubb[84]=aryubb[3]*aryubb[204];
   aryubb[82]=1./2.*aryubb[82] + 9*aryubb[84];
   aryubb[82]=MMt*aryubb[3]*aryubb[82];
   aryubb[69]=aryubb[82] + aryubb[72] + 1./48.*aryubb[78] + aryubb[69];
   aryubb[69]=MMt*aryubb[69];
   aryubb[72]= - aryubb[173] + aryubb[117];
   aryubb[72]=aryubb[3]*aryubb[72];
   aryubb[78]= - 1./3. - 1./4.*aryubb[43];
   aryubb[78]=aryubb[8]*aryubb[78];
   aryubb[82]=31./24. + aryubb[43];
   aryubb[82]=aryubb[11]*aryubb[82];
   aryubb[84]= - aryubb[13]*aryubb[43];
   aryubb[85]= - MMH*aryubb[42];
   aryubb[72]=1./4.*aryubb[84] + 1./2.*aryubb[72] + 41./216.*aryubb[12]
    + 1./3.*aryubb[85] + aryubb[78] + 1./9.*aryubb[82];
   aryubb[78]=1 + aryubb[43];
   aryubb[78]=aryubb[13]*aryubb[3]*aryubb[78];
   aryubb[82]=MMt*aryubb[177];
   aryubb[78]=aryubb[82] - 1./48.*aryubb[43] + aryubb[78];
   aryubb[78]=MMt*aryubb[78];
   aryubb[72]=1./2.*aryubb[72] + aryubb[78];
   aryubb[72]=aryubb[297]*aryubb[72];
   aryubb[76]= - 4051./108. + aryubb[76];
   aryubb[74]=227./3.*aryubb[74] + 125./9.*aryubb[96] + 1./3.*
   aryubb[76] + aryubb[81];
   aryubb[76]= - 125./24.*aryubb[54] - 17*aryubb[3];
   aryubb[76]=aryubb[13]*aryubb[76];
   aryubb[74]=1./8.*aryubb[74] + 1./3.*aryubb[76];
   aryubb[74]=aryubb[13]*aryubb[74];
   aryubb[76]=401./324. + aryubb[179];
   aryubb[76]=1./2.*aryubb[76] + aryubb[136];
   aryubb[76]=17./9.*aryubb[43] - 3./4.*aryubb[40] + 1./2.*aryubb[76]
    - 1./9.*aryubb[45];
   aryubb[76]=1./2.*aryubb[76] - 17./27.*aryubb[47];
   aryubb[76]=25./96.*aryubb[183] + 25./96.*aryubb[182] - 5./108.*
   aryubb[9] + 1./4.*aryubb[42] + 1./4.*aryubb[76] - 2./3.*aryubb[39];
   aryubb[76]=aryubb[11]*aryubb[76];
   aryubb[71]=aryubb[71] - 17./8. + aryubb[230];
   aryubb[71]=aryubb[8]*aryubb[71];
   aryubb[78]=85*aryubb[29] + 37*aryubb[34] - 107./9.*aryubb[6];
   aryubb[78]=841./192.*aryubb[32] + 1./16.*aryubb[78] + aryubb[33];
   aryubb[78]= - 389./144.*aryubb[25] - 43./96.*aryubb[17] - 9745./1728.
   *aryubb[16] - 61./24.*aryubb[24] + 1./3.*aryubb[78] - 7./32.*
   aryubb[27];
   aryubb[71]=1./3.*aryubb[78] + 1./8.*aryubb[71];
   aryubb[78]= - 17./9.*aryubb[8] - 27./2.*aryubb[11];
   aryubb[78]=aryubb[11]*aryubb[78];
   aryubb[78]=aryubb[78] + 39*aryubb[264];
   aryubb[78]=aryubb[3]*aryubb[78];
   aryubb[81]= - 49./54.*aryubb[42] + 59./27.*aryubb[39] + 5./9. + 
   aryubb[166];
   aryubb[81]=MMH*aryubb[81];
   aryubb[69]=1./4.*aryubb[72] + 1./2.*aryubb[69] + 1./4.*aryubb[74] + 
   1./16.*aryubb[78] + 349./576.*aryubb[12] + 1./8.*aryubb[81] + 1./2.*
   aryubb[71] + aryubb[76];
   aryubb[69]=aryubb[297]*aryubb[69];
   aryubb[71]= - 1235./36. + aryubb[179];
   aryubb[71]=aryubb[272] + aryubb[271] - 1./54.*aryubb[9] + 
   aryubb[263] + aryubb[259] - 11./6.*aryubb[39] + 23./54.*aryubb[47]
    - 17./12.*aryubb[43] - 3./8.*aryubb[40] - 17./36.*aryubb[45] - 3./4.
   *aryubb[44] + 1./8.*aryubb[71] + 1./9.*aryubb[229];
   aryubb[71]=aryubb[11]*aryubb[71];
   aryubb[72]=137./32.*aryubb[35] + 7./3.*aryubb[33] + aryubb[275] - 
   139./432.*aryubb[6] + 2./9.*aryubb[18] + aryubb[274];
   aryubb[74]= - 17./2.*aryubb[12] - aryubb[8] + 5*aryubb[11];
   aryubb[74]=aryubb[12]*aryubb[74];
   aryubb[74]=aryubb[277] + aryubb[74];
   aryubb[74]=aryubb[3]*aryubb[74];
   aryubb[59]=aryubb[69] + aryubb[63] + 1./2.*aryubb[65] + 1./8.*
   aryubb[74] + aryubb[64] + aryubb[185] + 1./2.*aryubb[71] + 1./8.*
   aryubb[59] - 17./18.*aryubb[25] - 1019./288.*aryubb[17] - 193./162.*
   aryubb[16] + aryubb[283] + aryubb[282] - 697./576.*aryubb[27] + 241./
   576.*aryubb[32] + 1./27.*aryubb[281] + aryubb[279] + 1./3.*
   aryubb[72] + aryubb[278];
   aryubb[59]=aryubb[297]*aryubb[59];
   aryubb[63]= - 5./18.*aryubb[35] - aryubb[31];
   aryubb[63]=7./2.*aryubb[63] - 5./9.*aryubb[32];
   aryubb[63]=aryubb[203] + 877./18.*aryubb[25] + aryubb[236] + 29./6.*
   aryubb[16] + aryubb[134] + aryubb[88] + 5*aryubb[63] - 37./18.*
   aryubb[27];
   aryubb[64]=419./8. + 73*aryubb[47];
   aryubb[64]=1./9.*aryubb[64] - aryubb[48];
   aryubb[64]=aryubb[11]*aryubb[64];
   aryubb[63]=aryubb[141] + aryubb[248] + aryubb[245] + 1./2.*
   aryubb[63] + 1./3.*aryubb[64];
   aryubb[64]=aryubb[140] + 127./48. + aryubb[139];
   aryubb[64]=aryubb[249] + aryubb[80] + aryubb[143] + 1./2.*aryubb[64]
    + aryubb[79];
   aryubb[64]=aryubb[252] + aryubb[251] - aryubb[42] + aryubb[250] + 1./
   2.*aryubb[64] - aryubb[39];
   aryubb[64]=aryubb[13]*aryubb[64];
   aryubb[65]=aryubb[254] + aryubb[253] + aryubb[148] + aryubb[147] + 
   869./384. + aryubb[145];
   aryubb[65]=aryubb[157] + aryubb[77] + 1./2.*aryubb[65] + aryubb[256]
   ;
   aryubb[65]=MMt*aryubb[65];
   aryubb[63]=aryubb[65] + 1./4.*aryubb[63] + aryubb[64];
   aryubb[63]=MMt*aryubb[63];
   aryubb[64]= - 11./9.*aryubb[35] + aryubb[197];
   aryubb[65]=aryubb[48] + 13./8. + 3*aryubb[46];
   aryubb[65]=aryubb[8]*aryubb[65];
   aryubb[64]=3./2.*aryubb[65] + 2441./108.*aryubb[25] + 149./36.*
   aryubb[17] + 401./216.*aryubb[16] - 5./8.*aryubb[23] + 75./8.*
   aryubb[24] - 11./8.*aryubb[27] + 7./4.*aryubb[64] - 13./27.*
   aryubb[32];
   aryubb[65]=aryubb[138] + aryubb[95];
   aryubb[65]=MMH*aryubb[65];
   aryubb[69]= - 11./8. - 7./3.*aryubb[47];
   aryubb[69]=1./3.*aryubb[69] - aryubb[48];
   aryubb[69]=aryubb[11]*aryubb[69];
   aryubb[71]= - 7./6. + aryubb[92];
   aryubb[71]=aryubb[12]*aryubb[71];
   aryubb[64]=1./2.*aryubb[141] + 1./4.*aryubb[71] + 1./4.*aryubb[65]
    + 1./2.*aryubb[64] + 1./3.*aryubb[69];
   aryubb[65]=aryubb[140] + 731./432. + aryubb[139];
   aryubb[65]=3./4.*aryubb[48] + aryubb[80] + aryubb[143] + 1./2.*
   aryubb[65] + aryubb[79];
   aryubb[69]=aryubb[144] + 31*aryubb[11];
   aryubb[69]=1./4.*aryubb[69] + 19*aryubb[12];
   aryubb[69]=aryubb[3]*aryubb[69];
   aryubb[65]=3./4.*aryubb[102] + 1./2.*aryubb[69] + 1./4.*aryubb[65]
    - aryubb[42];
   aryubb[65]=aryubb[13]*aryubb[65];
   aryubb[69]= - 11./36.*aryubb[48] + aryubb[148] + aryubb[147] + 11245.
   /3456. + aryubb[145];
   aryubb[71]= - 1 - aryubb[48];
   aryubb[71]=aryubb[12]*aryubb[71];
   aryubb[72]=aryubb[11]*aryubb[48];
   aryubb[71]=3./2.*aryubb[71] + aryubb[150] + 3./4.*aryubb[72];
   aryubb[71]=aryubb[3]*aryubb[71];
   aryubb[69]=1./2.*aryubb[69] + aryubb[71];
   aryubb[71]= - 3./8. + aryubb[46];
   aryubb[71]=3*aryubb[71] - aryubb[50];
   aryubb[74]= - 1./2.*aryubb[48];
   aryubb[71]=aryubb[74] + 1./2.*aryubb[71] + aryubb[135];
   aryubb[71]=aryubb[3]*aryubb[71];
   aryubb[71]=aryubb[71] + 3./2.*aryubb[156];
   aryubb[71]=aryubb[13]*aryubb[71];
   aryubb[69]=aryubb[158] + 1./2.*aryubb[69] + 3*aryubb[71];
   aryubb[69]=MMt*aryubb[69];
   aryubb[64]=aryubb[69] + 1./4.*aryubb[64] + aryubb[65];
   aryubb[64]=MMt*aryubb[64];
   aryubb[65]=aryubb[42] + aryubb[257] + aryubb[255] + 13./3. + 
   aryubb[239];
   aryubb[65]=aryubb[8]*aryubb[65];
   aryubb[69]= - 17./2.*aryubb[42] + 241./24. + 13*aryubb[39];
   aryubb[69]=aryubb[11]*aryubb[69];
   aryubb[71]=31./432.*aryubb[42] - 55./432.*aryubb[39] + 3./16.*
   aryubb[40] + 3./8.*aryubb[44] - 19./9. + 27./16.*aryubb[37];
   aryubb[71]=MMH*aryubb[71];
   aryubb[65]=1./2.*aryubb[71] + 1./36.*aryubb[69] + 1./4.*aryubb[65]
    + 7./32.*aryubb[25] + 9./32.*aryubb[17] + 1./6.*aryubb[16] + 19./32.
   *aryubb[23] + 65./36.*aryubb[24] - 13./32.*aryubb[31] - 7./16.*
   aryubb[30] - 37./288.*aryubb[29] + 13./144.*aryubb[34] - 1./9.*
   aryubb[36] - 15./32.*aryubb[28];
   aryubb[65]=MMH*aryubb[65];
   aryubb[69]= - 47*aryubb[8] + 161./18.*aryubb[11];
   aryubb[69]=aryubb[11]*aryubb[69];
   aryubb[69]=aryubb[149] + 1./9.*aryubb[69];
   aryubb[71]= - 65*aryubb[8] + 437./9.*aryubb[11];
   aryubb[76]=293./216. + aryubb[42];
   aryubb[76]=MMH*aryubb[76];
   aryubb[71]=419./36.*aryubb[12] + 1./8.*aryubb[71] + aryubb[76];
   aryubb[71]=aryubb[12]*aryubb[71];
   aryubb[76]= - 9*aryubb[8] - 139./6.*aryubb[11];
   aryubb[77]=aryubb[164] + aryubb[42];
   aryubb[77]=MMH*aryubb[77];
   aryubb[76]=2141./432.*aryubb[13] - 257./18.*aryubb[12] + 1./4.*
   aryubb[76] + aryubb[77];
   aryubb[76]=aryubb[13]*aryubb[76];
   aryubb[65]=1./4.*aryubb[76] + 1./4.*aryubb[71] + 1./8.*aryubb[69] + 
   aryubb[65];
   aryubb[69]= - aryubb[11]*aryubb[42];
   aryubb[69]=31./108.*aryubb[162] + 31./108.*aryubb[8] + aryubb[69];
   aryubb[69]=MMH*aryubb[69];
   aryubb[71]=aryubb[111] + 1./9.*aryubb[11];
   aryubb[71]=aryubb[11]*aryubb[71];
   aryubb[76]= - 31./108.*MMH + aryubb[210] + 7./9.*aryubb[11];
   aryubb[76]=aryubb[12]*aryubb[76];
   aryubb[77]= - aryubb[11] + aryubb[12];
   aryubb[78]=aryubb[13]*aryubb[77];
   aryubb[69]=1./4.*aryubb[78] + 1./2.*aryubb[76] + aryubb[71] + 1./2.*
   aryubb[69];
   aryubb[71]=aryubb[208] - aryubb[12];
   aryubb[71]=aryubb[3]*aryubb[71];
   aryubb[71]= - 31./72. + aryubb[71];
   aryubb[71]=aryubb[13]*aryubb[71];
   aryubb[72]=aryubb[3]*aryubb[72];
   aryubb[72]= - 31./36.*aryubb[48] + 3*aryubb[72];
   aryubb[72]=MMt*aryubb[72];
   aryubb[71]=1./2.*aryubb[72] + 1./48.*aryubb[77] + aryubb[71];
   aryubb[71]=MMt*aryubb[71];
   aryubb[69]=1./2.*aryubb[69] + aryubb[71];
   aryubb[69]=aryubb[297]*aryubb[69];
   aryubb[64]=1./4.*aryubb[69] + 1./2.*aryubb[65] + aryubb[64];
   aryubb[64]=aryubb[297]*aryubb[64];
   aryubb[65]= - 43*aryubb[8] + 245./18.*aryubb[11];
   aryubb[65]=391./72.*aryubb[13] + aryubb[265] + 1./6.*aryubb[65] + 
   aryubb[163];
   aryubb[65]=aryubb[13]*aryubb[65];
   aryubb[69]=aryubb[267] + 27*aryubb[11];
   aryubb[69]= - 19./2.*aryubb[12] + 1./4.*aryubb[69] + aryubb[269];
   aryubb[69]=aryubb[12]*aryubb[69];
   aryubb[71]=aryubb[8] - 59./108.*aryubb[11];
   aryubb[71]=aryubb[11]*aryubb[71];
   aryubb[63]=aryubb[64] + aryubb[63] + 1./8.*aryubb[65] + 1./8.*
   aryubb[69] + aryubb[205] + aryubb[270] + 1./3.*aryubb[71];
   aryubb[63]=aryubb[297]*aryubb[63];
   aryubb[64]=aryubb[8] + 1./2.*aryubb[162];
   aryubb[64]=MMH*aryubb[64];
   aryubb[64]=aryubb[64] + 3./2.*aryubb[120];
   aryubb[64]=aryubb[12]*aryubb[64];
   aryubb[64]=aryubb[114] + aryubb[119] + aryubb[64];
   aryubb[64]=aryubb[297]*aryubb[64];
   aryubb[64]=aryubb[258] + 1./2.*aryubb[64];
   aryubb[64]=aryubb[52]*aryubb[297]*aryubb[64];
   aryubb[65]=aryubb[160] - 11./18.*aryubb[11];
   aryubb[69]=3./4.*aryubb[13];
   aryubb[65]=aryubb[69] - aryubb[12] + 1./4.*aryubb[65] + aryubb[162];
   aryubb[65]=aryubb[13]*aryubb[65];
   aryubb[71]= - 13./3.*aryubb[11] + 1./2.*MMH;
   aryubb[71]=aryubb[12]*aryubb[71];
   aryubb[72]= - aryubb[8] + aryubb[85];
   aryubb[72]=MMH*aryubb[72];
   aryubb[71]=1./2.*aryubb[72] + aryubb[71];
   aryubb[65]=1./24.*aryubb[71] + aryubb[65];
   aryubb[71]=aryubb[8]*aryubb[48];
   aryubb[76]=13./18.*aryubb[266] + 3*aryubb[71];
   aryubb[77]=1./3.*aryubb[154];
   aryubb[76]=1./4.*aryubb[76] + aryubb[77];
   aryubb[78]=1./2. + aryubb[97];
   aryubb[78]=1./8.*aryubb[78] + 3*aryubb[146];
   aryubb[78]=aryubb[13]*aryubb[78];
   aryubb[79]=aryubb[13]*aryubb[178];
   aryubb[79]=1./16.*aryubb[48] + 3*aryubb[79];
   aryubb[79]=MMt*aryubb[79];
   aryubb[76]=aryubb[79] + 1./2.*aryubb[76] + aryubb[78];
   aryubb[76]=MMt*aryubb[76];
   aryubb[65]=1./2.*aryubb[65] + aryubb[76];
   aryubb[65]=MMt*aryubb[65];
   aryubb[72]=aryubb[72] + aryubb[116];
   aryubb[69]=aryubb[69] - aryubb[12] + aryubb[162] + 7./4.*aryubb[8]
    - 1./3.*aryubb[11];
   aryubb[69]=aryubb[13]*aryubb[69];
   aryubb[69]=1./48.*aryubb[72] + aryubb[69];
   aryubb[71]=3./4.*aryubb[71] + aryubb[77];
   aryubb[71]=aryubb[79] + 1./2.*aryubb[71] + aryubb[78];
   aryubb[71]=MMt*aryubb[71];
   aryubb[69]=1./2.*aryubb[69] + aryubb[71];
   aryubb[69]=MMt*aryubb[69];
   aryubb[71]= - aryubb[8]*aryubb[42];
   aryubb[76]=aryubb[11]*aryubb[42];
   aryubb[76]=1./4.*aryubb[71] + 1./9.*aryubb[76];
   aryubb[76]=MMH*aryubb[76];
   aryubb[76]=aryubb[76] + aryubb[262] + 1./9.*aryubb[171];
   aryubb[76]=MMH*aryubb[76];
   aryubb[72]=aryubb[13]*aryubb[72];
   aryubb[77]=aryubb[206] - 1./9.*aryubb[11];
   aryubb[77]=aryubb[12]*MMH*aryubb[77];
   aryubb[76]=1./4.*aryubb[72] + aryubb[76] + aryubb[77];
   aryubb[69]=1./4.*aryubb[76] + aryubb[69];
   aryubb[69]=aryubb[297]*aryubb[69];
   aryubb[71]=aryubb[172] + 1./9.*aryubb[260] + 1./8.*aryubb[71];
   aryubb[71]=MMH*aryubb[71];
   aryubb[76]= - 1./3.*aryubb[8] + aryubb[108];
   aryubb[76]=aryubb[11]*aryubb[76];
   aryubb[71]=aryubb[71] - 1./8.*aryubb[142] + 1./3.*aryubb[76];
   aryubb[71]=MMH*aryubb[71];
   aryubb[71]=1./8.*aryubb[72] + aryubb[71] + 1./2.*aryubb[77];
   aryubb[65]=aryubb[69] + 1./2.*aryubb[71] + aryubb[65];
   aryubb[65]=aryubb[297]*aryubb[65];
   aryubb[65]=1./2.*aryubb[261] + aryubb[65];
   aryubb[65]=aryubb[4]*aryubb[297]*aryubb[65];
   aryubb[69]=aryubb[11] + aryubb[13];
   aryubb[69]=aryubb[13]*aryubb[69];
   aryubb[71]=2*MMt + aryubb[86] + aryubb[11] - 2*aryubb[25] + 
   aryubb[32] - aryubb[16];
   aryubb[71]=MMt*aryubb[71];
   aryubb[69]=80./9.*aryubb[71] + 80./9.*aryubb[69] + 1./27.*
   aryubb[173] + 10*aryubb[215];
   aryubb[63]=1./2.*aryubb[65] + 1./8.*aryubb[64] + 1./3.*aryubb[69] + 
   aryubb[63];
   aryubb[63]=aryubb[4]*aryubb[63];
   aryubb[64]=17./3.*aryubb[11] - 5./4.*aryubb[12];
   aryubb[64]=aryubb[12]*aryubb[64];
   aryubb[65]= - 55*aryubb[11] - 47./2.*aryubb[12];
   aryubb[65]=1./4.*aryubb[65] + 13*aryubb[13];
   aryubb[65]=aryubb[13]*aryubb[65];
   aryubb[64]=1./2.*aryubb[64] + 1./3.*aryubb[65];
   aryubb[64]=aryubb[297]*aryubb[64];
   aryubb[65]= - 35./6.*aryubb[12] + 1./2.*aryubb[209] + 13./16.*
   aryubb[8] + 53./9.*aryubb[11];
   aryubb[65]=aryubb[12]*aryubb[65];
   aryubb[69]= - 371./6.*aryubb[12] + aryubb[8] - 7./9.*aryubb[11];
   aryubb[69]=1./2.*aryubb[69] + 67./3.*aryubb[13];
   aryubb[69]=aryubb[13]*aryubb[69];
   aryubb[64]=1./3.*aryubb[64] + 1./4.*aryubb[69] + 1./2.*aryubb[115]
    + aryubb[65];
   aryubb[64]=aryubb[297]*aryubb[64];
   aryubb[65]=aryubb[52]*aryubb[297]*aryubb[89];
   aryubb[64]=aryubb[64] + 1./2.*aryubb[65];
   aryubb[64]=aryubb[52]*aryubb[64];
   aryubb[65]= - 310./9. + aryubb[229];
   aryubb[65]= - 2*aryubb[45] + 1./3.*aryubb[65] - 10*aryubb[70];
   aryubb[65]=aryubb[11]*aryubb[65];
   aryubb[69]= - aryubb[13]*aryubb[54];
   aryubb[69]=2*aryubb[69] - 23./3. + 2*aryubb[96];
   aryubb[69]=aryubb[13]*aryubb[69];
   aryubb[71]= - 32*aryubb[18] + 5./3.*aryubb[6];
   aryubb[72]=aryubb[15]*aryubb[280];
   aryubb[76]= - aryubb[70]*aryubb[33];
   aryubb[77]=1489./108. + 10*aryubb[70];
   aryubb[77]=aryubb[16]*aryubb[77];
   aryubb[78]=3 + 5./3.*aryubb[70];
   aryubb[78]=aryubb[17]*aryubb[78];
   aryubb[79]=13 + 5*aryubb[70];
   aryubb[79]=aryubb[12]*aryubb[79];
   aryubb[80]= - 19 + aryubb[202];
   aryubb[80]=MMt*aryubb[80];
   aryubb[59]=aryubb[63] + 1./4.*aryubb[64] + aryubb[59] + 8./27.*
   aryubb[80] + 8./9.*aryubb[69] + 4./3.*aryubb[79] + 1./3.*aryubb[65]
    - 56./27.*aryubb[25] + 4*aryubb[78] + 1./3.*aryubb[77] + 10./3.*
   aryubb[76] + 11./36.*aryubb[27] + 28./27.*aryubb[32] + 16./27.*
   aryubb[72] - 11./36.*aryubb[35] + 1./27.*aryubb[71] - 6*aryubb[33];
   aryubb[59]=aryubb[4]*aryubb[59];
   aryubb[63]= - 1 - 1./2.*aryubb[70];
   aryubb[64]=aryubb[24]*aryubb[63];
   aryubb[65]=aryubb[23]*aryubb[290];
   aryubb[69]=aryubb[17]*aryubb[63];
   aryubb[71]=aryubb[25]*aryubb[63];
   aryubb[64]=aryubb[71] + aryubb[69] + aryubb[64] + aryubb[65];
   aryubb[64]=aryubb[56]*aryubb[64];
   aryubb[64]=13./4.*aryubb[290] + aryubb[64];
   aryubb[64]=aryubb[22]*aryubb[56]*aryubb[64];
   aryubb[65]=aryubb[55]*aryubb[176];
   aryubb[69]=3./4. + aryubb[65];
   aryubb[69]=aryubb[55]*aryubb[69];
   aryubb[69]=aryubb[124] + aryubb[69] + aryubb[122];
   aryubb[69]=aryubb[26]*aryubb[69];
   aryubb[63]=aryubb[22]*aryubb[73]*aryubb[63];
   aryubb[71]=aryubb[8]*aryubb[63];
   aryubb[76]= - aryubb[26]*aryubb[93];
   aryubb[63]=1./8.*aryubb[76] + 3*aryubb[63];
   aryubb[77]=aryubb[12]*aryubb[63];
   aryubb[78]=aryubb[13]*aryubb[63];
   aryubb[63]=MMt*aryubb[63];
   aryubb[79]=aryubb[70]*aryubb[290];
   aryubb[79]=3./2. + aryubb[79];
   aryubb[73]=aryubb[22]*aryubb[73]*aryubb[79];
   aryubb[79]=aryubb[26]*aryubb[93];
   aryubb[73]=1./32.*aryubb[79] + 3*aryubb[73];
   aryubb[73]=MMZ*aryubb[73];
   aryubb[80]= - aryubb[93]*aryubb[54];
   aryubb[81]=aryubb[11]*aryubb[76];
   aryubb[63]=9*aryubb[73] + aryubb[63] + 3./2.*aryubb[78] + 3./2.*
   aryubb[77] + 3./16.*aryubb[81] + 9./2.*aryubb[71] + 9./2.*aryubb[64]
    + aryubb[80] + 3./16.*aryubb[69];
   aryubb[63]=MMZ*aryubb[63];
   aryubb[64]=aryubb[76] + 3*aryubb[220];
   aryubb[64]=aryubb[12]*aryubb[64];
   aryubb[69]=aryubb[93]*aryubb[54];
   aryubb[71]= - aryubb[26]*aryubb[55];
   aryubb[73]=9./32.*aryubb[71];
   aryubb[76]=aryubb[11]*aryubb[79];
   aryubb[64]=9./4.*aryubb[64] + 9./8.*aryubb[76] + 27./8.*aryubb[165]
    + 63./32.*aryubb[211] + 4*aryubb[69] + aryubb[73];
   aryubb[64]=aryubb[13]*aryubb[64];
   aryubb[77]=aryubb[55]*aryubb[129];
   aryubb[78]= - 5./4. + aryubb[77];
   aryubb[78]=aryubb[55]*aryubb[78];
   aryubb[82]=aryubb[17]*aryubb[93];
   aryubb[84]=aryubb[25]*aryubb[93];
   aryubb[78]=aryubb[84] + aryubb[78] + aryubb[82];
   aryubb[78]=aryubb[26]*aryubb[78];
   aryubb[82]=aryubb[79] + 3*aryubb[161];
   aryubb[84]=aryubb[12]*aryubb[82];
   aryubb[82]=aryubb[13]*aryubb[82];
   aryubb[78]=aryubb[82] + aryubb[84] + aryubb[76] + aryubb[174] + 
   aryubb[78] + 3*aryubb[207];
   aryubb[78]=MMt*aryubb[78];
   aryubb[82]=23./6. + aryubb[70];
   aryubb[82]=aryubb[70]*aryubb[82];
   aryubb[65]=9./8.*aryubb[226] + 9./8.*aryubb[214] + 9./8.*aryubb[65]
    - 83./48. + aryubb[82];
   aryubb[65]=aryubb[26]*aryubb[65];
   aryubb[82]=aryubb[71] + 7*aryubb[211];
   aryubb[76]=aryubb[76] + 1./4.*aryubb[82] + aryubb[174];
   aryubb[76]=aryubb[12]*aryubb[76];
   aryubb[82]=1./2. - aryubb[40];
   aryubb[82]=aryubb[155]*aryubb[82];
   aryubb[84]=aryubb[297]*aryubb[82];
   aryubb[82]=3./2.*aryubb[84] + 19./24.*aryubb[26] + 3*aryubb[82];
   aryubb[82]=aryubb[297]*aryubb[82];
   aryubb[84]=aryubb[16]*aryubb[54];
   aryubb[67]=aryubb[67] + aryubb[84];
   aryubb[67]=aryubb[55]*aryubb[67];
   aryubb[84]= - aryubb[54] + aryubb[67];
   aryubb[84]=aryubb[55]*aryubb[84];
   aryubb[85]=aryubb[25]*aryubb[69];
   aryubb[84]=aryubb[84] + 2*aryubb[85];
   aryubb[73]=2*aryubb[69] + aryubb[73];
   aryubb[73]=aryubb[11]*aryubb[73];
   aryubb[66]=aryubb[155]*aryubb[66];
   aryubb[63]=3*aryubb[63] + 1./2.*aryubb[82] + 1./4.*aryubb[78] + 
   aryubb[64] + 3*aryubb[66] + 9./8.*aryubb[76] + aryubb[73] + 63./32.*
   aryubb[218] + 9./8.*aryubb[217] + 2*aryubb[84] + 1./4.*aryubb[65];
   aryubb[63]=MMZ*aryubb[63];
   aryubb[64]= - 881./72. + aryubb[103];
   aryubb[64]= - 19./9.*aryubb[47] + aryubb[87] + aryubb[40] + 31./6.*
   aryubb[45] + 1./3.*aryubb[64] + aryubb[44];
   aryubb[65]=aryubb[12]*aryubb[170];
   aryubb[65]=aryubb[175] + aryubb[65];
   aryubb[65]=aryubb[3]*aryubb[65];
   aryubb[66]=4./9.*aryubb[9];
   aryubb[64]=3./4.*aryubb[65] + aryubb[66] + aryubb[190] + aryubb[74]
    + 1./2.*aryubb[64] + aryubb[193];
   aryubb[64]=aryubb[3]*aryubb[64];
   aryubb[65]=19./12.*aryubb[213] + 19./12.*aryubb[109] + 19./12.*
   aryubb[181] + aryubb[123] - 125./48.*aryubb[54];
   aryubb[73]= - 3301./18. + 7*aryubb[45];
   aryubb[74]=aryubb[3]*aryubb[175];
   aryubb[73]=3./2.*aryubb[74] + 35./18.*aryubb[9] + 7*aryubb[39] + 119.
   /18.*aryubb[47] - 17./6.*aryubb[43] + 1./12.*aryubb[73] + aryubb[40]
   ;
   aryubb[73]=aryubb[3]*aryubb[73];
   aryubb[74]= - aryubb[57] + 5./48.*aryubb[54];
   aryubb[76]= - aryubb[297]*aryubb[3]*aryubb[43];
   aryubb[73]=1./6.*aryubb[76] + 25./4.*aryubb[74] + aryubb[73];
   aryubb[73]=aryubb[297]*aryubb[73];
   aryubb[74]= - aryubb[13]*aryubb[26];
   aryubb[76]= - MMt*aryubb[26];
   aryubb[64]=1./4.*aryubb[73] + 19./48.*aryubb[76] + 19./96.*
   aryubb[74] + 1./8.*aryubb[65] + aryubb[64];
   aryubb[64]=aryubb[297]*aryubb[64];
   aryubb[65]=aryubb[55]*aryubb[54];
   aryubb[73]= - 1./6. - aryubb[70];
   aryubb[74]=aryubb[26]*aryubb[73];
   aryubb[76]=1./8.*aryubb[74];
   aryubb[78]=aryubb[11]*aryubb[80];
   aryubb[80]=aryubb[12]*aryubb[71];
   aryubb[69]=aryubb[13]*aryubb[69];
   aryubb[69]=8*aryubb[69] + 9./8.*aryubb[80] + 8*aryubb[78] + 8./3.*
   aryubb[65] + aryubb[76];
   aryubb[69]=aryubb[13]*aryubb[69];
   aryubb[78]=aryubb[70]*aryubb[27];
   aryubb[80]=aryubb[16]*aryubb[73];
   aryubb[82]=aryubb[17]*aryubb[73];
   aryubb[73]=aryubb[25]*aryubb[73];
   aryubb[73]=aryubb[73] + aryubb[82] + aryubb[80] + 1./6.*aryubb[27]
    + aryubb[78];
   aryubb[73]=aryubb[26]*aryubb[73];
   aryubb[78]=365./18. + aryubb[229];
   aryubb[80]= - 4*aryubb[70];
   aryubb[82]= - 53./3. + aryubb[80];
   aryubb[82]=aryubb[45]*aryubb[82];
   aryubb[66]=aryubb[66] + 16./9.*aryubb[47] - 6*aryubb[43] + 
   aryubb[82] + 1./3.*aryubb[78] + aryubb[80];
   aryubb[66]=aryubb[3]*aryubb[66];
   aryubb[78]=aryubb[17]*aryubb[55];
   aryubb[80]=aryubb[25]*aryubb[55];
   aryubb[77]=5./6.*aryubb[80] + 5./6.*aryubb[78] + 5./6.*aryubb[77] - 
   85./6. - aryubb[70];
   aryubb[77]=aryubb[26]*aryubb[77];
   aryubb[78]=aryubb[26]*aryubb[55];
   aryubb[80]=aryubb[11]*aryubb[78];
   aryubb[77]=aryubb[77] + 5./6.*aryubb[80];
   aryubb[81]=5./12.*aryubb[78] + aryubb[81];
   aryubb[82]=aryubb[12]*aryubb[81];
   aryubb[77]=1./2.*aryubb[77] + aryubb[82];
   aryubb[79]=aryubb[12]*aryubb[79];
   aryubb[79]=1./2.*aryubb[81] + aryubb[79];
   aryubb[79]=aryubb[13]*aryubb[79];
   aryubb[77]=1./2.*aryubb[77] + aryubb[79];
   aryubb[77]=MMt*aryubb[77];
   aryubb[79]=aryubb[25]*aryubb[65];
   aryubb[67]=4*aryubb[79] + aryubb[54] + 2*aryubb[67];
   aryubb[76]=4./3.*aryubb[65] + aryubb[76];
   aryubb[76]=aryubb[11]*aryubb[76];
   aryubb[74]=aryubb[74] + 9*aryubb[80];
   aryubb[74]=aryubb[12]*aryubb[74];
   aryubb[63]=aryubb[63] + aryubb[64] + aryubb[77] + aryubb[69] + 
   aryubb[66] + 1./8.*aryubb[74] + aryubb[76] + 2./3.*aryubb[67] + 1./8.
   *aryubb[73];
   aryubb[63]=MMZ*aryubb[63];
   aryubb[64]=aryubb[11]*aryubb[71];
   aryubb[66]=aryubb[12]*aryubb[78];
   aryubb[66]=5./3.*aryubb[66] + 1./4.*aryubb[26] + 1./3.*aryubb[64];
   aryubb[66]=aryubb[13]*aryubb[66];
   aryubb[64]=1./8.*aryubb[26] + 2./3.*aryubb[64];
   aryubb[64]=aryubb[12]*aryubb[64];
   aryubb[67]=1 + aryubb[47];
   aryubb[67]=aryubb[3]*aryubb[67];
   aryubb[62]=1./2.*aryubb[66] + 32./9.*aryubb[67] + 1./8.*aryubb[62]
    + aryubb[64];
   aryubb[62]=MMt*aryubb[62];
   aryubb[64]=10*aryubb[72] - 13*aryubb[18] - 101./3.*aryubb[14];
   aryubb[64]=aryubb[15]*aryubb[64];
   aryubb[64]=aryubb[64] + 1849./18. - 13*aryubb[19];
   aryubb[64]=1./3.*aryubb[64] + 4*aryubb[131];
   aryubb[66]= - 1 - 2*aryubb[70];
   aryubb[66]=aryubb[70]*aryubb[66];
   aryubb[64]=1./3.*aryubb[64] + 10*aryubb[66];
   aryubb[66]=aryubb[15]*aryubb[14]*aryubb[51];
   aryubb[66]= - 13./3.*aryubb[51] + 2*aryubb[66];
   aryubb[66]=1./3.*aryubb[66] - aryubb[54];
   aryubb[67]= - aryubb[55]*aryubb[54];
   aryubb[67]=16*aryubb[67] + 5./8.*aryubb[26];
   aryubb[67]=aryubb[11]*aryubb[67];
   aryubb[66]=8./3.*aryubb[66] + aryubb[67];
   aryubb[67]=pow(aryubb[51],2);
   aryubb[65]= - 1./9.*aryubb[67] + 2*aryubb[65];
   aryubb[65]=aryubb[13]*aryubb[65];
   aryubb[65]=8./3.*aryubb[65] + 32./9.*aryubb[3] + 1./3.*aryubb[66] + 
   1./8.*aryubb[68];
   aryubb[65]=aryubb[13]*aryubb[65];
   aryubb[66]= - 5 + aryubb[75];
   aryubb[66]=aryubb[12]*aryubb[66];
   aryubb[66]=aryubb[137] + 2*aryubb[66];
   aryubb[66]=aryubb[3]*aryubb[66];
   aryubb[67]= - 1 - 1./3.*aryubb[70];
   aryubb[67]=aryubb[45]*aryubb[67];
   aryubb[68]= - aryubb[16]*aryubb[54];
   aryubb[69]=aryubb[12]*aryubb[109];
   aryubb[58]=aryubb[60] + aryubb[59] + aryubb[61] + aryubb[63] + 
   aryubb[58] + aryubb[62] + aryubb[65] + aryubb[66] + 1./3.*aryubb[69]
    + 4./9.*aryubb[110] + 4./27.*aryubb[9] + 8./9.*aryubb[284] + 52./27.
   *aryubb[47] + 4./9.*aryubb[68] + 1./3.*aryubb[64] + 4*aryubb[67];
   aryubb[58]=aryubb[21]*aryubb[58];
   aryubb[59]=1./2. + aryubb[9];
   aryubb[60]=aryubb[20]*aryubb[59];
   aryubb[59]=aryubb[1]*aryubb[59];
   aryubb[61]=3*aryubb[60] + aryubb[59];
   aryubb[61]=aryubb[12]*aryubb[61];
   aryubb[62]= - 1 + aryubb[9];
   aryubb[63]=aryubb[20]*aryubb[62];
   aryubb[62]=aryubb[1]*aryubb[62];
   aryubb[62]=3*aryubb[63] + aryubb[62];
   aryubb[62]=MMZ*aryubb[62];
   aryubb[63]= - aryubb[7] + aryubb[17];
   aryubb[64]=aryubb[20]*aryubb[63];
   aryubb[63]=aryubb[1]*aryubb[63];
   aryubb[65]=3*aryubb[20] + aryubb[1];
   aryubb[65]=aryubb[13]*aryubb[65];
   aryubb[62]=3./2.*aryubb[62] + 1./2.*aryubb[65] + aryubb[61] + 3*
   aryubb[64] + aryubb[63];
   aryubb[63]=3*aryubb[9];
   aryubb[64]=1 + aryubb[10];
   aryubb[65]=aryubb[64] + aryubb[63];
   aryubb[66]=aryubb[20]*aryubb[65];
   aryubb[65]=aryubb[1]*aryubb[65];
   aryubb[65]=3*aryubb[66] + aryubb[65];
   aryubb[65]=aryubb[12]*aryubb[65];
   aryubb[66]= - 3*aryubb[9];
   aryubb[64]=aryubb[64] + aryubb[66];
   aryubb[67]=aryubb[20]*aryubb[64];
   aryubb[64]=aryubb[1]*aryubb[64];
   aryubb[64]=3*aryubb[67] + aryubb[64];
   aryubb[64]=aryubb[13]*aryubb[64];
   aryubb[67]= - 1 + aryubb[10];
   aryubb[68]=aryubb[20]*aryubb[67];
   aryubb[67]=aryubb[1]*aryubb[67];
   aryubb[67]=3*aryubb[68] + aryubb[67];
   aryubb[67]=1./2.*MMZ*aryubb[67];
   aryubb[68]=aryubb[106] + aryubb[17] + 1./2.*aryubb[53] - aryubb[7];
   aryubb[69]=aryubb[20]*aryubb[68];
   aryubb[68]=aryubb[1]*aryubb[68];
   aryubb[64]=aryubb[67] + 1./4.*aryubb[64] + 1./4.*aryubb[65] + 3*
   aryubb[69] + aryubb[68];
   aryubb[64]=MMZ*aryubb[64];
   aryubb[65]=aryubb[20]*aryubb[10];
   aryubb[68]=aryubb[1]*aryubb[10];
   aryubb[65]=3*aryubb[65] + aryubb[68];
   aryubb[68]=aryubb[12]*aryubb[65];
   aryubb[69]= - aryubb[20]*aryubb[10];
   aryubb[70]= - aryubb[1]*aryubb[10];
   aryubb[69]=3*aryubb[69] + aryubb[70];
   aryubb[70]=aryubb[13]*aryubb[69];
   aryubb[71]= - aryubb[25] + aryubb[17] + aryubb[53] - aryubb[7];
   aryubb[72]=aryubb[20]*aryubb[71];
   aryubb[71]=aryubb[1]*aryubb[71];
   aryubb[68]=aryubb[70] + aryubb[68] + 3*aryubb[72] + aryubb[71];
   aryubb[68]=1./2.*aryubb[52]*aryubb[127]*aryubb[68];
   aryubb[64]=aryubb[64] + aryubb[68];
   aryubb[64]=aryubb[52]*aryubb[64];
   aryubb[62]=1./2.*aryubb[62] + aryubb[64];
   aryubb[62]=aryubb[52]*aryubb[62];
   aryubb[64]=aryubb[10] - aryubb[9];
   aryubb[70]=aryubb[20]*aryubb[64];
   aryubb[64]=aryubb[1]*aryubb[64];
   aryubb[71]=3*aryubb[70] + aryubb[64];
   aryubb[72]=aryubb[12]*aryubb[71];
   aryubb[73]= - aryubb[10] + aryubb[9];
   aryubb[74]=aryubb[20]*aryubb[73];
   aryubb[73]=aryubb[1]*aryubb[73];
   aryubb[75]=3*aryubb[74] + aryubb[73];
   aryubb[76]=aryubb[13]*aryubb[75];
   aryubb[77]=aryubb[72] + aryubb[76];
   aryubb[77]=aryubb[52]*MMZ*aryubb[77];
   aryubb[78]=MMZ*aryubb[71];
   aryubb[77]=aryubb[77] + aryubb[72] + aryubb[78];
   aryubb[77]=aryubb[52]*aryubb[77];
   aryubb[78]=aryubb[11]*aryubb[71];
   aryubb[72]=1./2.*aryubb[78] + aryubb[72];
   aryubb[72]=aryubb[3]*aryubb[72];
   aryubb[78]=aryubb[3]*aryubb[71];
   aryubb[79]=MMZ*aryubb[78];
   aryubb[72]=1./4.*aryubb[77] + aryubb[79] + 1./8.*aryubb[75] + 
   aryubb[72];
   aryubb[77]=aryubb[8]*aryubb[71];
   aryubb[71]=aryubb[13]*aryubb[71];
   aryubb[71]=aryubb[77] + aryubb[71];
   aryubb[64]=aryubb[70] + 1./3.*aryubb[64];
   aryubb[70]=aryubb[13]*aryubb[3]*aryubb[75];
   aryubb[64]=1./16.*aryubb[64] + aryubb[70];
   aryubb[64]=MMt*aryubb[64];
   aryubb[64]=1./8.*aryubb[71] + aryubb[64];
   aryubb[64]=aryubb[4]*aryubb[64];
   aryubb[64]=1./2.*aryubb[72] + aryubb[64];
   aryubb[64]=aryubb[91]*aryubb[64];
   aryubb[70]=1 - aryubb[10];
   aryubb[71]=aryubb[70] + aryubb[63];
   aryubb[72]=aryubb[20]*aryubb[71];
   aryubb[71]=aryubb[1]*aryubb[71];
   aryubb[71]=3*aryubb[72] + aryubb[71];
   aryubb[71]=aryubb[11]*aryubb[71];
   aryubb[72]=aryubb[17] + 1./2.*aryubb[16] - aryubb[7] - 1./2.*
   aryubb[6];
   aryubb[77]=aryubb[20]*aryubb[72];
   aryubb[72]=aryubb[1]*aryubb[72];
   aryubb[61]=aryubb[61] + 1./4.*aryubb[71] + 3*aryubb[77] + aryubb[72]
   ;
   aryubb[61]=aryubb[3]*aryubb[61];
   aryubb[71]=aryubb[8]*aryubb[75];
   aryubb[71]=aryubb[71] + aryubb[76];
   aryubb[72]=aryubb[74] + 1./3.*aryubb[73];
   aryubb[73]=aryubb[13]*aryubb[78];
   aryubb[72]=1./16.*aryubb[72] + aryubb[73];
   aryubb[72]=MMt*aryubb[72];
   aryubb[71]=1./8.*aryubb[71] + aryubb[72];
   aryubb[71]=aryubb[4]*aryubb[71];
   aryubb[72]=1./4. + aryubb[10];
   aryubb[72]=1./2.*aryubb[72] - aryubb[9];
   aryubb[73]=aryubb[20]*aryubb[72];
   aryubb[72]=aryubb[1]*aryubb[72];
   aryubb[72]=3*aryubb[73] + aryubb[72];
   aryubb[73]= - 11 + 21./2.*aryubb[9];
   aryubb[73]=aryubb[20]*aryubb[73];
   aryubb[74]= - 11./3. + aryubb[104];
   aryubb[75]=aryubb[1]*aryubb[74];
   aryubb[73]=aryubb[73] + aryubb[75];
   aryubb[73]=MMZ*aryubb[3]*aryubb[73];
   aryubb[61]=aryubb[64] + aryubb[71] + 1./2.*aryubb[62] + 1./2.*
   aryubb[73] + 1./8.*aryubb[72] + aryubb[61];
   aryubb[61]=aryubb[91]*aryubb[61];
   aryubb[62]= - 3*aryubb[10];
   aryubb[64]= - 19./3.*aryubb[9];
   aryubb[71]=aryubb[64] - 7./18. + aryubb[62];
   aryubb[71]=aryubb[20]*aryubb[71];
   aryubb[72]=aryubb[66] + 1./6. - aryubb[10];
   aryubb[72]=aryubb[1]*aryubb[72];
   aryubb[71]=aryubb[71] + aryubb[72];
   aryubb[71]=aryubb[12]*aryubb[71];
   aryubb[72]= - 47./6. + 19*aryubb[9];
   aryubb[72]=aryubb[20]*aryubb[72];
   aryubb[73]= - 7./6. + aryubb[63];
   aryubb[73]=aryubb[1]*aryubb[73];
   aryubb[72]=1./3.*aryubb[72] + aryubb[73];
   aryubb[72]=aryubb[13]*aryubb[72];
   aryubb[73]=aryubb[20]*aryubb[70];
   aryubb[70]=aryubb[1]*aryubb[70];
   aryubb[70]=3*aryubb[73] + aryubb[70];
   aryubb[70]=MMZ*aryubb[70];
   aryubb[73]=aryubb[121] - aryubb[17] - 1./2.*aryubb[53] + aryubb[7];
   aryubb[76]=aryubb[20]*aryubb[73];
   aryubb[73]=aryubb[1]*aryubb[73];
   aryubb[70]=aryubb[70] + 1./2.*aryubb[72] + 1./2.*aryubb[71] + 3*
   aryubb[76] + aryubb[73];
   aryubb[70]=MMZ*aryubb[70];
   aryubb[69]=aryubb[12]*aryubb[69];
   aryubb[65]=aryubb[13]*aryubb[65];
   aryubb[71]=aryubb[25] - aryubb[17] - aryubb[53] + aryubb[7];
   aryubb[72]=aryubb[20]*aryubb[71];
   aryubb[71]=aryubb[1]*aryubb[71];
   aryubb[65]=aryubb[65] + aryubb[69] + 3*aryubb[72] + aryubb[71];
   aryubb[65]=aryubb[52]*aryubb[127]*aryubb[65];
   aryubb[65]=aryubb[70] + aryubb[65];
   aryubb[65]=aryubb[52]*aryubb[65];
   aryubb[62]=101./9. + aryubb[62];
   aryubb[62]=1./2.*aryubb[62] + aryubb[64];
   aryubb[62]=aryubb[20]*aryubb[62];
   aryubb[64]=13./3. - aryubb[10];
   aryubb[64]=1./2.*aryubb[64] + aryubb[66];
   aryubb[64]=aryubb[1]*aryubb[64];
   aryubb[62]=aryubb[62] + aryubb[64];
   aryubb[62]=MMZ*aryubb[62];
   aryubb[64]=1./3. - aryubb[9];
   aryubb[66]=aryubb[20]*aryubb[64];
   aryubb[64]=aryubb[1]*aryubb[64];
   aryubb[64]=5./3.*aryubb[66] + aryubb[64];
   aryubb[66]=aryubb[12]*aryubb[64];
   aryubb[62]=aryubb[65] + aryubb[66] + 1./2.*aryubb[62];
   aryubb[62]=aryubb[52]*aryubb[62];
   aryubb[65]= - 11./12.*aryubb[9] + 5./9. - 3./4.*aryubb[10];
   aryubb[65]=aryubb[20]*aryubb[65];
   aryubb[69]= - 1./4.*aryubb[10];
   aryubb[70]=1./3. + aryubb[69];
   aryubb[71]=aryubb[70] - 3./4.*aryubb[9];
   aryubb[71]=aryubb[1]*aryubb[71];
   aryubb[65]=aryubb[65] + aryubb[71];
   aryubb[71]=aryubb[8]*aryubb[65];
   aryubb[72]= - 1./3.*aryubb[9] - 1 + 1./3.*aryubb[10];
   aryubb[73]=aryubb[20]*aryubb[72];
   aryubb[72]=aryubb[1]*aryubb[72];
   aryubb[72]=aryubb[73] + 1./3.*aryubb[72];
   aryubb[72]=aryubb[11]*aryubb[72];
   aryubb[73]=aryubb[13]*aryubb[65];
   aryubb[72]=aryubb[73] + aryubb[71] + aryubb[72];
   aryubb[70]=1./3.*aryubb[70] + aryubb[169];
   aryubb[70]=aryubb[1]*aryubb[70];
   aryubb[69]= - 11./36.*aryubb[9] + 5./27. + aryubb[69];
   aryubb[69]=aryubb[20]*aryubb[69];
   aryubb[69]=aryubb[69] + aryubb[70];
   aryubb[70]=11./3.*aryubb[9] - 20./9. + 3*aryubb[10];
   aryubb[70]=aryubb[20]*aryubb[70];
   aryubb[63]=aryubb[63] - 4./3. + aryubb[10];
   aryubb[63]=aryubb[1]*aryubb[63];
   aryubb[63]=aryubb[70] + aryubb[63];
   aryubb[63]=aryubb[13]*aryubb[3]*aryubb[63];
   aryubb[63]=1./4.*aryubb[69] + aryubb[63];
   aryubb[63]=MMt*aryubb[63];
   aryubb[69]=1./2.*aryubb[72] + aryubb[63];
   aryubb[69]=aryubb[4]*aryubb[69];
   aryubb[65]=aryubb[11]*aryubb[65];
   aryubb[66]=aryubb[65] + 2*aryubb[66];
   aryubb[66]=aryubb[3]*aryubb[66];
   aryubb[70]=31./4.*aryubb[10];
   aryubb[72]= - 43./4.*aryubb[9] + 13 + aryubb[70];
   aryubb[72]=aryubb[20]*aryubb[72];
   aryubb[76]=5./36.*aryubb[9] + 1 + 31./36.*aryubb[10];
   aryubb[76]=aryubb[1]*aryubb[76];
   aryubb[72]=1./3.*aryubb[72] + aryubb[76];
   aryubb[76]= - 2*aryubb[10];
   aryubb[77]= - 29./6.*aryubb[9] + 43./9. + aryubb[76];
   aryubb[77]=aryubb[20]*aryubb[77];
   aryubb[76]=17./3. + aryubb[76];
   aryubb[76]=1./3.*aryubb[76] - 5./2.*aryubb[9];
   aryubb[76]=aryubb[1]*aryubb[76];
   aryubb[76]=aryubb[77] + aryubb[76];
   aryubb[76]=MMZ*aryubb[3]*aryubb[76];
   aryubb[61]=aryubb[61] + aryubb[69] + 1./2.*aryubb[62] + aryubb[76]
    + 1./12.*aryubb[72] + aryubb[66];
   aryubb[61]=aryubb[91]*aryubb[61];
   aryubb[62]= - 1./3. + aryubb[9];
   aryubb[66]=aryubb[20]*aryubb[62];
   aryubb[62]=aryubb[1]*aryubb[62];
   aryubb[62]=5./3.*aryubb[66] + aryubb[62];
   aryubb[66]=aryubb[12]*aryubb[62];
   aryubb[64]=aryubb[13]*aryubb[64];
   aryubb[64]=aryubb[67] + aryubb[66] + aryubb[64];
   aryubb[64]=MMZ*aryubb[64];
   aryubb[64]=aryubb[64] + aryubb[68];
   aryubb[64]=aryubb[52]*aryubb[64];
   aryubb[66]=MMZ*aryubb[62];
   aryubb[64]=aryubb[66] + aryubb[64];
   aryubb[64]=aryubb[52]*aryubb[64];
   aryubb[66]=11./9.*aryubb[9] - 53./27. + aryubb[10];
   aryubb[66]=aryubb[20]*aryubb[66];
   aryubb[67]= - 13./3. + aryubb[10];
   aryubb[67]=1./3.*aryubb[67] + aryubb[9];
   aryubb[67]=aryubb[1]*aryubb[67];
   aryubb[66]=aryubb[66] + aryubb[67];
   aryubb[66]=aryubb[11]*aryubb[66];
   aryubb[66]=aryubb[73] + aryubb[71] + 1./3.*aryubb[66];
   aryubb[67]= - 1 - aryubb[9];
   aryubb[68]=aryubb[20]*aryubb[67];
   aryubb[67]=aryubb[1]*aryubb[67];
   aryubb[67]=11./9.*aryubb[68] + aryubb[67];
   aryubb[67]=aryubb[11]*aryubb[67];
   aryubb[68]=aryubb[6] - aryubb[16];
   aryubb[69]=aryubb[20]*aryubb[68];
   aryubb[68]=aryubb[1]*aryubb[68];
   aryubb[67]=aryubb[67] + 11./9.*aryubb[69] + aryubb[68];
   aryubb[67]=aryubb[297]*aryubb[67];
   aryubb[63]=1./6.*aryubb[67] + 1./2.*aryubb[66] + aryubb[63];
   aryubb[63]=aryubb[297]*aryubb[63];
   aryubb[66]=5./3.*aryubb[20] + aryubb[1];
   aryubb[66]=aryubb[11]*aryubb[66];
   aryubb[66]=2./3.*aryubb[66] + 5./3.*aryubb[69] + aryubb[68];
   aryubb[63]=4./9.*aryubb[66] + aryubb[63];
   aryubb[63]=aryubb[4]*aryubb[63];
   aryubb[59]=11./3.*aryubb[60] + 3*aryubb[59];
   aryubb[59]=aryubb[11]*aryubb[59];
   aryubb[60]= - aryubb[6] + aryubb[16];
   aryubb[66]=aryubb[20]*aryubb[60];
   aryubb[60]=aryubb[1]*aryubb[60];
   aryubb[59]=aryubb[59] + 11./3.*aryubb[66] + 3*aryubb[60];
   aryubb[59]=aryubb[3]*aryubb[59];
   aryubb[60]=631./12. - 55*aryubb[9];
   aryubb[66]=aryubb[20]*aryubb[60];
   aryubb[60]=aryubb[1]*aryubb[60];
   aryubb[60]=11./9.*aryubb[66] + aryubb[60];
   aryubb[59]=1./72.*aryubb[60] + aryubb[59];
   aryubb[59]=aryubb[297]*aryubb[59];
   aryubb[60]=77./36.*aryubb[9] + 43./27. + aryubb[70];
   aryubb[60]=aryubb[20]*aryubb[60];
   aryubb[66]=23./3. + aryubb[70];
   aryubb[66]=1./3.*aryubb[66] + 7./4.*aryubb[9];
   aryubb[66]=aryubb[1]*aryubb[66];
   aryubb[60]=aryubb[60] + aryubb[66];
   aryubb[65]=aryubb[3]*aryubb[65];
   aryubb[59]=1./2.*aryubb[59] + 1./36.*aryubb[60] + aryubb[65];
   aryubb[59]=aryubb[297]*aryubb[59];
   aryubb[60]= - 1./2.*aryubb[10];
   aryubb[65]=2./3. + aryubb[60];
   aryubb[65]=1./3.*aryubb[65] + aryubb[83];
   aryubb[65]=aryubb[1]*aryubb[65];
   aryubb[60]= - 11./18.*aryubb[9] + 10./27. + aryubb[60];
   aryubb[60]=aryubb[20]*aryubb[60];
   aryubb[60]=aryubb[60] + aryubb[65];
   aryubb[60]=aryubb[3]*aryubb[60];
   aryubb[65]=aryubb[20]*aryubb[74];
   aryubb[65]=11./9.*aryubb[65] + aryubb[75];
   aryubb[65]=aryubb[297]*aryubb[3]*aryubb[65];
   aryubb[60]=aryubb[60] + 1./2.*aryubb[65];
   aryubb[60]=aryubb[297]*aryubb[60];
   aryubb[62]=aryubb[3]*aryubb[62];
   aryubb[60]=4./3.*aryubb[62] + aryubb[60];
   aryubb[60]=MMZ*aryubb[60];
   aryubb[62]=aryubb[201] - 13./3.*aryubb[14];
   aryubb[62]=aryubb[15]*aryubb[62];
   aryubb[62]=4*aryubb[9] + 2*aryubb[62] + 19./3. + 4*aryubb[19];
   aryubb[65]=aryubb[20]*aryubb[62];
   aryubb[62]=aryubb[1]*aryubb[62];
   aryubb[62]=5./3.*aryubb[65] + aryubb[62];

      yubbret = aryubb[58] + aryubb[59] + aryubb[60] + aryubb[61] + 1./
      9.*aryubb[62] + aryubb[63] + 1./2.*aryubb[64];
      return yubbret;
}
