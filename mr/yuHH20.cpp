#include <HH.hpp>
std::complex<long double>
HH<OS>::my20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuHH[333], yuHHret;

    aryuHH[1]=double(nL + nH);
    aryuHH[2]=pow(CW,-1);
    aryuHH[3]=pow(MMH,-1);
    aryuHH[4]=pow(SW,-1);
    aryuHH[5]=Tsil::I2(0,0,MMZ,mu2);
    aryuHH[6]=Tsil::I2(0,0,MMW,mu2);
    aryuHH[7]=Tsil::B(MMH,MMH,MMH,mu2);
    aryuHH[8]=pow(MMZ,-1);
    aryuHH[9]=std::real(Tsil::B(0,0,MMZ,mu2));
    aryuHH[10]=std::real(Tsil::B(0,0,MMW,mu2));
    aryuHH[11]=Tsil::B(MMZ,MMZ,MMH,mu2);
    aryuHH[12]=Tsil::B(MMW,MMW,MMH,mu2);
    aryuHH[13]=Tsil::B(0,MMZ,MMH,mu2);
    aryuHH[14]=Tsil::B(0,MMW,MMH,mu2);
    aryuHH[15]=Tsil::Beps(MMZ,MMZ,MMH,mu2);
    aryuHH[16]=Tsil::Beps(MMW,MMW,MMH,mu2);
    aryuHH[17]=Tsil::A(MMH,mu2);
    aryuHH[18]=Tsil::A(MMZ,mu2);
    aryuHH[19]=Tsil::A(MMW,mu2);
    aryuHH[20]=Tsil::Aeps(MMZ,mu2);
    aryuHH[21]=Tsil::Aeps(MMW,mu2);
    aryuHH[22]=protZZ00->Uxzuv(0);
    aryuHH[23]=protWW00->Uxzuv(0);
    aryuHH[24]=protZZ00->Txuv(0);
    aryuHH[25]=protWW00->Txuv(0);
    aryuHH[26]=double(nH);
    aryuHH[27]=Tsil::B(MMt,MMt,MMH,mu2);
    aryuHH[28]=Tsil::A(MMt,mu2);
    aryuHH[29]=Tsil::B(MMt,MMt,MMZ,mu2);
    aryuHH[30]=Tsil::B(0,MMt,MMW,mu2);
    aryuHH[31]=double(nL);
    aryuHH[32]=Tsil::I2(MMH,MMt,MMt,mu2);
    aryuHH[33]=Tsil::I2(MMZ,MMt,MMt,mu2);
    aryuHH[34]=Tsil::I2(0,MMW,MMt,mu2);
    aryuHH[35]=Tsil::B(MMH,MMt,MMt,mu2);
    aryuHH[36]=Tsil::B(MMZ,MMH,MMZ,mu2);
    aryuHH[37]=Tsil::B(MMZ,MMt,MMt,mu2);
    aryuHH[38]=Tsil::B(MMW,MMH,MMW,mu2);
    aryuHH[39]=Tsil::B(MMW,MMZ,MMW,mu2);
    aryuHH[40]=Tsil::B(MMW,MMW,MMZ,mu2);
    aryuHH[41]=std::real(Tsil::B(0,MMW,MMt,mu2));
    aryuHH[42]=Tsil::Beps(MMH,MMH,MMH,mu2);
    aryuHH[43]=Tsil::Beps(MMt,MMt,MMH,mu2);
    aryuHH[44]=Tsil::Aeps(MMH,mu2);
    aryuHH[45]=Tsil::Aeps(MMt,mu2);
    aryuHH[46]=prottttt0->M(0);
    aryuHH[47]=prottttt0->Vxzuv(0);
    aryuHH[48]=prottttt0->Suxv(0);
    aryuHH[49]=protHtHtt->M(0);
    aryuHH[50]=protZtZtt->M(0);
    aryuHH[51]=protWtWt0->M(0);
    aryuHH[52]=protttttH->M(0);
    aryuHH[53]=protttttZ->M(0);
    aryuHH[54]=protHtHtt->Uzxyv(0);
    aryuHH[55]=protZtZtt->Uzxyv(0);
    aryuHH[56]=protWtWt0->Uzxyv(0);
    aryuHH[57]=protHtHtt->Uuyxv(0);
    aryuHH[58]=protZtZtt->Uuyxv(0);
    aryuHH[59]=protWtWt0->Uuyxv(0);
    aryuHH[60]=protZtZtt->Txuv(0);
    aryuHH[61]=protWtWt0->Suxv(0);
    aryuHH[62]=protWtWt0->Txuv(0);
    aryuHH[63]=protZtZtt->Tyzv(0);
    aryuHH[64]=protWtWt0->Tyzv(0);
    aryuHH[65]=protHtHtt->Tyzv(0);
    aryuHH[66]=protHtHtt->Svyz(0);
    aryuHH[67]=protZtZtt->Svyz(0);
    aryuHH[68]=double(boson);
    aryuHH[69]=Tsil::I2(MMH,MMH,MMH,mu2);
    aryuHH[70]=Tsil::I2(MMH,MMZ,MMZ,mu2);
    aryuHH[71]=Tsil::I2(MMH,MMW,MMW,mu2);
    aryuHH[72]=Tsil::I2(MMW,MMW,MMZ,mu2);
    aryuHH[73]=protHHHHH->M(0);
    aryuHH[74]=protHZHZZ->M(0);
    aryuHH[75]=protHWHWW->M(0);
    aryuHH[76]=protZZZZH->M(0);
    aryuHH[77]=protZWZWW->M(0);
    aryuHH[78]=protWWWWH->M(0);
    aryuHH[79]=protWWWWZ->M(0);
    aryuHH[80]=protWWWW0->M(0);
    aryuHH[81]=protWWWW0->Vxzuv(0);
    aryuHH[82]=protHHHHH->Uzxyv(0);
    aryuHH[83]=protHZHZZ->Uzxyv(0);
    aryuHH[84]=protHWHWW->Uzxyv(0);
    aryuHH[85]=protHZHZZ->Uuyxv(0);
    aryuHH[86]=protZWZWW->Uzxyv(0);
    aryuHH[87]=protHWHWW->Uuyxv(0);
    aryuHH[88]=protZWZWW->Uuyxv(0);
    aryuHH[89]=protHZHZZ->Tyzv(0);
    aryuHH[90]=protZWZWW->Txuv(0);
    aryuHH[91]=protZWZWW->Tuxv(0);
    aryuHH[92]=protHWHWW->Tuxv(0);
    aryuHH[93]=protHZHZZ->Svyz(0);
    aryuHH[94]=protHWHWW->Suxv(0);
    aryuHH[95]=protZWZWW->Suxv(0);
    aryuHH[96]=protWWWW0->Suxv(0);
    aryuHH[97]=1/(4*MMt - MMZ);
    aryuHH[98]=1/( - 4*MMW + MMH);
    aryuHH[99]=1/( - 4*MMZ + MMH);
   aryuHH[100]= - 13*aryuHH[39];
   aryuHH[101]=aryuHH[100] - 31 - 13./2.*aryuHH[40];
   aryuHH[102]=aryuHH[99]*aryuHH[17];
   aryuHH[103]=3./4.*aryuHH[102];
   aryuHH[104]= - aryuHH[18]*aryuHH[99];
   aryuHH[105]=3./2.*aryuHH[104];
   aryuHH[106]=3./8.*aryuHH[11];
   aryuHH[107]=27./4.*aryuHH[7];
   aryuHH[108]=3*aryuHH[12];
   aryuHH[109]= - 3*aryuHH[36];
   aryuHH[101]=aryuHH[106] + aryuHH[105] + aryuHH[108] + aryuHH[103] + 
   aryuHH[107] + aryuHH[109] + 1./6.*aryuHH[101] + aryuHH[38];
   aryuHH[101]=aryuHH[11]*aryuHH[101];
   aryuHH[110]=9*aryuHH[86];
   aryuHH[111]=143./2. + aryuHH[110];
   aryuHH[112]= - 27*aryuHH[83];
   aryuHH[111]=1./4.*aryuHH[111] + aryuHH[112];
   aryuHH[113]=9./8.*aryuHH[40];
   aryuHH[114]=3./2.*aryuHH[36];
   aryuHH[111]=27./2.*aryuHH[7] + aryuHH[114] + aryuHH[113] - 33./8.*
   aryuHH[15] - 27./16.*aryuHH[90] + 27./2.*aryuHH[42] - 1./16.*
   aryuHH[91] + 3*aryuHH[85] + 1./2.*aryuHH[111] - 18*aryuHH[89];
   aryuHH[111]=aryuHH[99]*aryuHH[111];
   aryuHH[115]=3*aryuHH[74] + 1./4.*aryuHH[77];
   aryuHH[116]= - 19./2. + aryuHH[88];
   aryuHH[116]=1./2.*aryuHH[39] - 1./2.*aryuHH[16] - 9./2.*aryuHH[90]
    + 1./2.*aryuHH[116] - aryuHH[91];
   aryuHH[116]=aryuHH[98]*aryuHH[116];
   aryuHH[117]= - aryuHH[12]*aryuHH[98]*aryuHH[39];
   aryuHH[118]= - 3./4.*aryuHH[40];
   aryuHH[119]= - aryuHH[36] - 1 + aryuHH[118];
   aryuHH[119]=aryuHH[11]*aryuHH[99]*aryuHH[119];
   aryuHH[120]= - 2*aryuHH[76];
   aryuHH[121]=1./4.*aryuHH[79];
   aryuHH[111]=3./2.*aryuHH[119] + 1./8.*aryuHH[117] + aryuHH[111] + 1./
   4.*aryuHH[116] + aryuHH[120] + 3*aryuHH[115] + aryuHH[121];
   aryuHH[111]=MMZ*aryuHH[111];
   aryuHH[115]= - aryuHH[95] + 1./2.*aryuHH[72];
   aryuHH[116]=3./4.*aryuHH[115] - aryuHH[93];
   aryuHH[117]=7./4.*aryuHH[70];
   aryuHH[119]=5./4.*aryuHH[44];
   aryuHH[122]=5./4.*aryuHH[17];
   aryuHH[116]=aryuHH[122] - 143./4.*aryuHH[20] + aryuHH[119] + 3*
   aryuHH[116] + aryuHH[117];
   aryuHH[116]=aryuHH[99]*aryuHH[116];
   aryuHH[123]=229./6. - 3*aryuHH[86];
   aryuHH[123]= - 17*aryuHH[88] + 1./2.*aryuHH[123] + aryuHH[112];
   aryuHH[124]= - 13./3.*aryuHH[89];
   aryuHH[123]=1./2.*aryuHH[123] + aryuHH[124];
   aryuHH[125]=1 - 5./2.*aryuHH[39];
   aryuHH[126]= - 1./4.*aryuHH[12];
   aryuHH[125]=1./3.*aryuHH[125] + aryuHH[126];
   aryuHH[125]=aryuHH[12]*aryuHH[125];
   aryuHH[127]= - 27*aryuHH[7];
   aryuHH[128]=aryuHH[127] + aryuHH[109] + 43 - 9./4.*aryuHH[40];
   aryuHH[128]=aryuHH[99]*aryuHH[128];
   aryuHH[128]=9./2.*aryuHH[98] + aryuHH[128];
   aryuHH[128]=aryuHH[18]*aryuHH[128];
   aryuHH[129]= - 3*aryuHH[74] - 1./2.*aryuHH[77];
   aryuHH[129]= - aryuHH[76] + 3*aryuHH[129] - 7./4.*aryuHH[79];
   aryuHH[129]=MMH*aryuHH[129];
   aryuHH[130]=aryuHH[115] - aryuHH[21];
   aryuHH[130]=aryuHH[98]*aryuHH[130];
   aryuHH[131]=3 - aryuHH[39];
   aryuHH[131]=aryuHH[98]*aryuHH[131];
   aryuHH[131]=aryuHH[131] + aryuHH[99];
   aryuHH[131]=aryuHH[19]*aryuHH[131];
   aryuHH[132]=3./2.*aryuHH[39];
   aryuHH[133]=27./4.*aryuHH[42];
   aryuHH[134]= - 5./8.*aryuHH[15];
   aryuHH[135]=1./2.*aryuHH[36];
   aryuHH[101]=1./2.*aryuHH[129] + aryuHH[111] + aryuHH[101] + 1./2.*
   aryuHH[128] + 9./4.*aryuHH[131] + aryuHH[125] + 1./2.*aryuHH[116] + 
   9./4.*aryuHH[130] + aryuHH[107] + aryuHH[135] + aryuHH[132] + 3./4.*
   aryuHH[40] + aryuHH[134] + 17./4.*aryuHH[16] - 9./8.*aryuHH[90] + 
   aryuHH[133] - 77./24.*aryuHH[91] + 1./2.*aryuHH[123] + aryuHH[85];
   aryuHH[101]=aryuHH[68]*aryuHH[101];
   aryuHH[111]=1./4.*aryuHH[115] - 3*aryuHH[93];
   aryuHH[116]=aryuHH[122] - 135./4.*aryuHH[20] + aryuHH[119] + 
   aryuHH[111] + aryuHH[117];
   aryuHH[116]=aryuHH[99]*aryuHH[116];
   aryuHH[123]= - 1./4.*aryuHH[40];
   aryuHH[125]=aryuHH[127] + aryuHH[109] + 39 + aryuHH[123];
   aryuHH[125]=aryuHH[99]*aryuHH[125];
   aryuHH[125]=1./2.*aryuHH[98] + aryuHH[125];
   aryuHH[125]=aryuHH[18]*aryuHH[125];
   aryuHH[110]=163./18. + aryuHH[110];
   aryuHH[110]= - 3*aryuHH[88] + 1./2.*aryuHH[110] + aryuHH[112];
   aryuHH[110]=1./2.*aryuHH[110] + aryuHH[124];
   aryuHH[128]= - aryuHH[12]*aryuHH[39];
   aryuHH[110]=1./2.*aryuHH[125] + 1./4.*aryuHH[131] + 1./12.*
   aryuHH[128] + 1./2.*aryuHH[116] + 1./4.*aryuHH[130] + aryuHH[107] + 
   aryuHH[135] + 1./6.*aryuHH[39] + 1./12.*aryuHH[40] - 17./8.*
   aryuHH[15] + 3./4.*aryuHH[16] - 37./12.*aryuHH[90] + aryuHH[133] - 3.
   /4.*aryuHH[91] + 1./2.*aryuHH[110] + aryuHH[85];
   aryuHH[116]=17./6.*aryuHH[39] + 3 + aryuHH[123];
   aryuHH[116]=1./2.*aryuHH[116] + aryuHH[38];
   aryuHH[125]=27./8.*aryuHH[7];
   aryuHH[129]=3./8.*aryuHH[102];
   aryuHH[130]=3./4.*aryuHH[104];
   aryuHH[131]=3./16.*aryuHH[11];
   aryuHH[136]=1./4.*aryuHH[12];
   aryuHH[137]= - 2*aryuHH[36];
   aryuHH[116]=aryuHH[131] + aryuHH[130] + aryuHH[136] + aryuHH[129] + 
   aryuHH[125] + 1./2.*aryuHH[116] + aryuHH[137];
   aryuHH[116]=aryuHH[11]*aryuHH[116];
   aryuHH[138]= - 27./2.*aryuHH[83] + 9 + 1./8.*aryuHH[86];
   aryuHH[139]=1./16.*aryuHH[40];
   aryuHH[140]=3./4.*aryuHH[36];
   aryuHH[138]=aryuHH[107] + aryuHH[140] + aryuHH[139] - 25./16.*
   aryuHH[15] - 3./32.*aryuHH[90] + aryuHH[133] + 3./2.*aryuHH[85] + 1./
   2.*aryuHH[138] - 9*aryuHH[89];
   aryuHH[138]=aryuHH[99]*aryuHH[138];
   aryuHH[141]=9*aryuHH[74];
   aryuHH[142]=1./2.*aryuHH[77];
   aryuHH[143]=aryuHH[141] + aryuHH[142];
   aryuHH[144]=aryuHH[109] - 3 + aryuHH[123];
   aryuHH[144]=aryuHH[11]*aryuHH[99]*aryuHH[144];
   aryuHH[145]= - 1 - aryuHH[90];
   aryuHH[146]=aryuHH[98]*aryuHH[145];
   aryuHH[138]=1./4.*aryuHH[144] + aryuHH[138] + 1./16.*aryuHH[146] + 1.
   /2.*aryuHH[143] - aryuHH[76];
   aryuHH[138]=MMZ*aryuHH[138];
   aryuHH[143]= - aryuHH[76] + aryuHH[121] - 9*aryuHH[74] + aryuHH[142]
   ;
   aryuHH[143]=MMH*aryuHH[143];
   aryuHH[110]=1./4.*aryuHH[143] + aryuHH[138] + 1./2.*aryuHH[110] + 
   aryuHH[116];
   aryuHH[110]=aryuHH[68]*aryuHH[110];
   aryuHH[116]=1./4.*aryuHH[40];
   aryuHH[138]=1 + aryuHH[116];
   aryuHH[143]=15*aryuHH[36];
   aryuHH[144]= - 17./2.*aryuHH[39];
   aryuHH[138]=aryuHH[11] + aryuHH[143] + 5*aryuHH[138] + aryuHH[144];
   aryuHH[138]=aryuHH[11]*aryuHH[138];
   aryuHH[146]= - 43./3. - 5./2.*aryuHH[86];
   aryuHH[147]=MMZ*aryuHH[76];
   aryuHH[148]=4./3.*aryuHH[89];
   aryuHH[149]= - 3./2.*aryuHH[85];
   aryuHH[138]=3./2.*aryuHH[147] + 1./4.*aryuHH[138] + aryuHH[140] + 
   aryuHH[139] + 29./16.*aryuHH[15] + 19./96.*aryuHH[90] + aryuHH[149]
    + 1./8.*aryuHH[146] + aryuHH[148];
   aryuHH[138]=MMZ*aryuHH[138];
   aryuHH[139]=17./3.*aryuHH[93];
   aryuHH[146]= - 37./2.*aryuHH[70];
   aryuHH[147]=77./6.*aryuHH[44];
   aryuHH[150]=41./6.*aryuHH[17];
   aryuHH[151]=aryuHH[150] + 233./6.*aryuHH[20] - 5./3.*aryuHH[21] + 
   aryuHH[147] + aryuHH[146] + aryuHH[139] + 7./3.*aryuHH[95] - 13./4.*
   aryuHH[72];
   aryuHH[116]=1./3. + aryuHH[116];
   aryuHH[116]=aryuHH[143] + 5*aryuHH[116] + aryuHH[144];
   aryuHH[116]=aryuHH[18]*aryuHH[116];
   aryuHH[132]=5./3. + aryuHH[132];
   aryuHH[132]=aryuHH[19]*aryuHH[132];
   aryuHH[143]=13*aryuHH[17];
   aryuHH[152]=aryuHH[18] + aryuHH[143] - 35*aryuHH[19];
   aryuHH[152]=aryuHH[11]*aryuHH[152];
   aryuHH[116]=1./2.*aryuHH[152] + aryuHH[116] + 1./2.*aryuHH[151] + 
   aryuHH[132];
   aryuHH[116]=1./4.*aryuHH[116] + aryuHH[138];
   aryuHH[116]=aryuHH[68]*aryuHH[116];
   aryuHH[132]= - 1./8.*aryuHH[90] - 25./8. - 3*aryuHH[89];
   aryuHH[132]=MMZ*aryuHH[132];
   aryuHH[138]=3./2.*aryuHH[44];
   aryuHH[151]=3./2.*aryuHH[17];
   aryuHH[152]=aryuHH[11]*aryuHH[17];
   aryuHH[153]=1./2.*aryuHH[19];
   aryuHH[111]=aryuHH[132] + 3./2.*aryuHH[152] + 25./4.*aryuHH[18] + 
   aryuHH[153] + aryuHH[151] + aryuHH[138] + aryuHH[111] + 3./2.*
   aryuHH[70];
   aryuHH[111]=aryuHH[68]*MMZ*aryuHH[111];
   aryuHH[132]=aryuHH[13] + 1 + aryuHH[24];
   aryuHH[154]=aryuHH[31]*aryuHH[132];
   aryuHH[132]=aryuHH[1]*aryuHH[132];
   aryuHH[132]=11./3.*aryuHH[154] + 3*aryuHH[132];
   aryuHH[132]=MMZ*aryuHH[132];
   aryuHH[154]=aryuHH[31]*aryuHH[5];
   aryuHH[155]=aryuHH[1]*aryuHH[5];
   aryuHH[156]= - 11./3.*aryuHH[31] - 3*aryuHH[1];
   aryuHH[156]=aryuHH[18]*aryuHH[156];
   aryuHH[132]=aryuHH[132] + aryuHH[156] + 11./3.*aryuHH[154] + 3*
   aryuHH[155];
   aryuHH[132]=MMZ*aryuHH[132];
   aryuHH[111]=1./2.*aryuHH[132] + aryuHH[111];
   aryuHH[111]=aryuHH[3]*aryuHH[111];
   aryuHH[132]= - 1 - 1./3.*aryuHH[24];
   aryuHH[132]=1./2.*aryuHH[132];
   aryuHH[154]=aryuHH[132] - 5./3.*aryuHH[22];
   aryuHH[155]= - 5./9.*aryuHH[13];
   aryuHH[156]=1./6.*aryuHH[9];
   aryuHH[154]=aryuHH[156] + 5./6.*aryuHH[15] + 1./2.*aryuHH[154] + 
   aryuHH[155];
   aryuHH[154]=aryuHH[31]*aryuHH[154];
   aryuHH[157]= - 5*aryuHH[22];
   aryuHH[158]= - 3 - aryuHH[24];
   aryuHH[158]=1./2.*aryuHH[158];
   aryuHH[159]=aryuHH[158] + aryuHH[157];
   aryuHH[159]=3./2.*aryuHH[9] + 15./2.*aryuHH[15] + 3./2.*aryuHH[159]
    - 5*aryuHH[13];
   aryuHH[159]=aryuHH[1]*aryuHH[159];
   aryuHH[154]=11*aryuHH[154] + aryuHH[159];
   aryuHH[159]=5./4.*aryuHH[9];
   aryuHH[160]=1 + aryuHH[159];
   aryuHH[161]=aryuHH[31]*aryuHH[160];
   aryuHH[160]=aryuHH[1]*aryuHH[160];
   aryuHH[160]=11./3.*aryuHH[161] + 3*aryuHH[160];
   aryuHH[160]=aryuHH[11]*aryuHH[160];
   aryuHH[154]=1./2.*aryuHH[154] + aryuHH[160];
   aryuHH[154]=MMZ*aryuHH[154];
   aryuHH[160]= - 11./2.*aryuHH[5] + 5*aryuHH[20];
   aryuHH[161]=aryuHH[31]*aryuHH[160];
   aryuHH[160]=aryuHH[1]*aryuHH[160];
   aryuHH[160]=11./3.*aryuHH[161] + 3*aryuHH[160];
   aryuHH[159]=2./3. + aryuHH[159];
   aryuHH[159]=aryuHH[31]*aryuHH[159];
   aryuHH[161]=2 + 15./4.*aryuHH[9];
   aryuHH[161]=aryuHH[1]*aryuHH[161];
   aryuHH[159]=11./3.*aryuHH[159] + aryuHH[161];
   aryuHH[159]=aryuHH[18]*aryuHH[159];
   aryuHH[111]=aryuHH[111] + aryuHH[116] + aryuHH[154] + 1./4.*
   aryuHH[160] + aryuHH[159];
   aryuHH[111]=aryuHH[3]*aryuHH[111];
   aryuHH[116]=7./6.*aryuHH[13];
   aryuHH[154]= - aryuHH[15] + aryuHH[116] - 55./72. + aryuHH[22];
   aryuHH[159]=1./3.*aryuHH[9];
   aryuHH[154]=1./2.*aryuHH[154] + aryuHH[159];
   aryuHH[154]=aryuHH[31]*aryuHH[154];
   aryuHH[160]=1./2.*aryuHH[5];
   aryuHH[161]=aryuHH[160] - aryuHH[20];
   aryuHH[162]=aryuHH[31]*aryuHH[161];
   aryuHH[161]=aryuHH[1]*aryuHH[161];
   aryuHH[163]=11./3.*aryuHH[162] + 3*aryuHH[161];
   aryuHH[163]=aryuHH[99]*aryuHH[163];
   aryuHH[164]=1./2. - aryuHH[9];
   aryuHH[165]=aryuHH[31]*aryuHH[164];
   aryuHH[164]=aryuHH[1]*aryuHH[164];
   aryuHH[166]=11./3.*aryuHH[165] + 3*aryuHH[164];
   aryuHH[166]=aryuHH[18]*aryuHH[99]*aryuHH[166];
   aryuHH[167]= - 3*aryuHH[15] + 7./2.*aryuHH[13] - 55./24. + 3*
   aryuHH[22];
   aryuHH[167]=1./2.*aryuHH[167] + aryuHH[9];
   aryuHH[167]=aryuHH[1]*aryuHH[167];
   aryuHH[154]=1./2.*aryuHH[166] + 1./2.*aryuHH[163] + 11./3.*
   aryuHH[154] + aryuHH[167];
   aryuHH[163]= - 1./2.*aryuHH[13];
   aryuHH[132]=aryuHH[163] + aryuHH[132] + 1./3.*aryuHH[22];
   aryuHH[166]=aryuHH[159] + aryuHH[132] - 1./3.*aryuHH[15];
   aryuHH[166]=aryuHH[31]*aryuHH[166];
   aryuHH[158]= - 3./2.*aryuHH[13] + aryuHH[158] + aryuHH[22];
   aryuHH[167]=aryuHH[9] + aryuHH[158] - aryuHH[15];
   aryuHH[168]=aryuHH[1]*aryuHH[167];
   aryuHH[166]=11*aryuHH[166] + 3*aryuHH[168];
   aryuHH[166]=aryuHH[99]*aryuHH[166];
   aryuHH[169]= - aryuHH[31]*aryuHH[9];
   aryuHH[170]= - aryuHH[1]*aryuHH[9];
   aryuHH[171]=11./3.*aryuHH[169] + 3*aryuHH[170];
   aryuHH[171]=aryuHH[11]*aryuHH[99]*aryuHH[171];
   aryuHH[166]=aryuHH[166] + aryuHH[171];
   aryuHH[166]=MMZ*aryuHH[166];
   aryuHH[171]= - aryuHH[18]*aryuHH[39];
   aryuHH[172]= - MMZ*aryuHH[11]*aryuHH[39];
   aryuHH[173]=aryuHH[19] - aryuHH[18];
   aryuHH[174]=aryuHH[11]*aryuHH[173];
   aryuHH[171]=aryuHH[172] + aryuHH[171] + aryuHH[174];
   aryuHH[171]=aryuHH[3]*aryuHH[68]*aryuHH[171];
   aryuHH[172]=aryuHH[11]*aryuHH[39];
   aryuHH[145]=aryuHH[145] + 1./3.*aryuHH[172];
   aryuHH[145]=aryuHH[68]*aryuHH[145];
   aryuHH[145]=aryuHH[145] + aryuHH[171];
   aryuHH[171]=pow(aryuHH[2],2);
   aryuHH[145]=aryuHH[171]*aryuHH[145];
   aryuHH[172]= - 1./4.*aryuHH[9];
   aryuHH[175]= - 1./3. + aryuHH[172];
   aryuHH[175]=aryuHH[31]*aryuHH[175];
   aryuHH[176]= - 3./4.*aryuHH[9];
   aryuHH[177]= - 1 + aryuHH[176];
   aryuHH[177]=aryuHH[1]*aryuHH[177];
   aryuHH[175]=11./3.*aryuHH[175] + aryuHH[177];
   aryuHH[175]=aryuHH[11]*aryuHH[175];
   aryuHH[110]=1./8.*aryuHH[145] + aryuHH[111] + aryuHH[110] + 1./4.*
   aryuHH[166] + 1./2.*aryuHH[154] + aryuHH[175];
   aryuHH[110]=aryuHH[171]*aryuHH[110];
   aryuHH[111]=11*aryuHH[40];
   aryuHH[145]=13*aryuHH[39];
   aryuHH[154]= - 3*aryuHH[38];
   aryuHH[166]=aryuHH[154] + aryuHH[145] + 35 + aryuHH[111];
   aryuHH[175]=1./2.*aryuHH[11];
   aryuHH[177]=6*aryuHH[36];
   aryuHH[178]=1./2.*aryuHH[12];
   aryuHH[166]=aryuHH[175] + aryuHH[178] + 1./2.*aryuHH[166] + 
   aryuHH[177];
   aryuHH[166]=aryuHH[11]*aryuHH[166];
   aryuHH[179]=aryuHH[12]*aryuHH[39];
   aryuHH[180]=aryuHH[142] + 3*aryuHH[76];
   aryuHH[180]=MMZ*aryuHH[180];
   aryuHH[181]=8./3.*aryuHH[89];
   aryuHH[182]= - 3*aryuHH[85];
   aryuHH[113]=aryuHH[180] + aryuHH[166] + 3./8.*aryuHH[179] + 
   aryuHH[114] + 1./8.*aryuHH[39] + aryuHH[113] + 73./8.*aryuHH[15] + 1.
   /8.*aryuHH[16] + 65./16.*aryuHH[90] - 7./16.*aryuHH[91] + 
   aryuHH[182] + aryuHH[181] - 1./8.*aryuHH[88] + 2 - 49./8.*aryuHH[86]
   ;
   aryuHH[113]=MMZ*aryuHH[113];
   aryuHH[139]=aryuHH[150] + 119./2.*aryuHH[20] - 19./3.*aryuHH[21] + 
   aryuHH[147] + aryuHH[146] + aryuHH[139] + 91./3.*aryuHH[95] - 149./4.
   *aryuHH[72];
   aryuHH[111]=aryuHH[154] + aryuHH[145] + 43./3. + aryuHH[111];
   aryuHH[146]=3./4.*aryuHH[12];
   aryuHH[111]=aryuHH[146] + 1./2.*aryuHH[111] + aryuHH[177];
   aryuHH[111]=aryuHH[18]*aryuHH[111];
   aryuHH[147]=17*aryuHH[18];
   aryuHH[143]=aryuHH[147] + aryuHH[143] - 19*aryuHH[19];
   aryuHH[143]=1./4.*aryuHH[11]*aryuHH[143];
   aryuHH[150]= - 3./4.*aryuHH[12] + 9./2. + 7*aryuHH[39];
   aryuHH[150]=aryuHH[19]*aryuHH[150];
   aryuHH[111]=aryuHH[113] + aryuHH[143] + aryuHH[111] + 1./4.*
   aryuHH[139] + aryuHH[150];
   aryuHH[111]=aryuHH[68]*aryuHH[111];
   aryuHH[113]= - 6*aryuHH[89];
   aryuHH[139]= - 5./2.*aryuHH[90] - 1./4.*aryuHH[91] - 35./4. + 
   aryuHH[113];
   aryuHH[139]=MMZ*aryuHH[139];
   aryuHH[150]=3*aryuHH[44];
   aryuHH[166]=3*aryuHH[17];
   aryuHH[139]=aryuHH[139] + 3*aryuHH[152] + aryuHH[147] + 10*
   aryuHH[19] + aryuHH[166] + aryuHH[150] + 3*aryuHH[70] + 5*
   aryuHH[115] - 6*aryuHH[93];
   aryuHH[139]=aryuHH[3]*aryuHH[68]*MMZ*aryuHH[139];
   aryuHH[147]= - 3./2.*aryuHH[10];
   aryuHH[179]= - 11./6.*aryuHH[9] + 10./9. + aryuHH[147];
   aryuHH[179]=aryuHH[31]*aryuHH[179];
   aryuHH[180]= - 1./2.*aryuHH[10];
   aryuHH[183]=2./3. + aryuHH[180];
   aryuHH[184]=aryuHH[183] - 3./2.*aryuHH[9];
   aryuHH[184]=aryuHH[1]*aryuHH[184];
   aryuHH[185]=aryuHH[179] + aryuHH[184];
   aryuHH[186]=aryuHH[11]*aryuHH[185];
   aryuHH[187]=MMZ*aryuHH[186];
   aryuHH[185]=aryuHH[18]*aryuHH[185];
   aryuHH[111]=aryuHH[139] + aryuHH[111] + aryuHH[185] + aryuHH[187];
   aryuHH[111]=aryuHH[3]*aryuHH[111];
   aryuHH[139]=1./2.*aryuHH[10];
   aryuHH[187]= - 2./3. + aryuHH[139];
   aryuHH[188]=1./2.*aryuHH[9];
   aryuHH[187]=1./3.*aryuHH[187] + aryuHH[188];
   aryuHH[187]=aryuHH[1]*aryuHH[187];
   aryuHH[189]=11./18.*aryuHH[9] - 10./27. + aryuHH[139];
   aryuHH[189]=aryuHH[31]*aryuHH[189];
   aryuHH[187]=aryuHH[189] + aryuHH[187];
   aryuHH[187]=aryuHH[11]*aryuHH[187];
   aryuHH[101]=aryuHH[110] + aryuHH[111] + aryuHH[187] + aryuHH[101];
   aryuHH[101]=aryuHH[171]*aryuHH[101];
   aryuHH[110]=593./4. - 3*aryuHH[24];
   aryuHH[110]=1./2.*aryuHH[110] - 15*aryuHH[22];
   aryuHH[111]= - 5./2.*aryuHH[13];
   aryuHH[189]= - 1./2.*aryuHH[20];
   aryuHH[190]=5*aryuHH[45] + aryuHH[189];
   aryuHH[190]=aryuHH[97]*aryuHH[190];
   aryuHH[191]=3./4.*aryuHH[9];
   aryuHH[192]=3*aryuHH[30];
   aryuHH[110]=aryuHH[191] + aryuHH[192] + 9./8.*aryuHH[190] + 3./4.*
   aryuHH[29] + 51./8.*aryuHH[15] + 45./32.*aryuHH[63] + 15*aryuHH[16]
    + 53./16.*aryuHH[60] - 21./8.*aryuHH[55] + aryuHH[111] + 9./32.*
   aryuHH[43] - 9./32.*aryuHH[58] - 15*aryuHH[56] + 1./4.*aryuHH[110]
    + 11*aryuHH[62];
   aryuHH[193]=25*aryuHH[29];
   aryuHH[194]=39./2. + aryuHH[193];
   aryuHH[195]=17./4.*aryuHH[9];
   aryuHH[196]=aryuHH[28]*aryuHH[97];
   aryuHH[197]= - 3*aryuHH[30];
   aryuHH[194]=9./2.*aryuHH[196] + aryuHH[195] + 1./4.*aryuHH[194] + 
   aryuHH[197];
   aryuHH[194]=aryuHH[11]*aryuHH[194];
   aryuHH[198]=aryuHH[43] - 3 - aryuHH[58];
   aryuHH[198]= - aryuHH[15] + 5./4.*aryuHH[63] - 1./2.*aryuHH[60] + 1./
   4.*aryuHH[198] + aryuHH[55];
   aryuHH[198]=aryuHH[97]*aryuHH[198];
   aryuHH[199]=aryuHH[27]*aryuHH[97];
   aryuHH[200]= - aryuHH[11]*aryuHH[97];
   aryuHH[198]=aryuHH[200] + aryuHH[198] + 1./4.*aryuHH[199];
   aryuHH[198]=MMZ*aryuHH[198];
   aryuHH[201]= - 3./2.*aryuHH[30];
   aryuHH[202]=7./2.*aryuHH[9] + aryuHH[201] + 6 + 11./2.*aryuHH[29];
   aryuHH[202]=aryuHH[12]*aryuHH[202];
   aryuHH[203]=1./4. + aryuHH[196];
   aryuHH[203]=aryuHH[27]*aryuHH[203];
   aryuHH[204]= - aryuHH[27]*aryuHH[97];
   aryuHH[205]=aryuHH[97] + aryuHH[204];
   aryuHH[205]=aryuHH[18]*aryuHH[205];
   aryuHH[206]= - aryuHH[28]*aryuHH[97];
   aryuHH[110]=9./16.*aryuHH[198] + 1./2.*aryuHH[194] + 9./32.*
   aryuHH[205] + 9./16.*aryuHH[203] + 45./16.*aryuHH[206] + 1./2.*
   aryuHH[110] + aryuHH[202];
   aryuHH[110]=MMZ*aryuHH[110];
   aryuHH[194]= - 1./4.*aryuHH[60];
   aryuHH[202]=aryuHH[194] + 1./4.*aryuHH[13] - aryuHH[62] - 1 + 1./4.*
   aryuHH[24];
   aryuHH[202]=MMZ*aryuHH[202];
   aryuHH[207]=1./4.*aryuHH[33];
   aryuHH[208]= - 1./2.*aryuHH[67];
   aryuHH[209]=3*aryuHH[28];
   aryuHH[202]=aryuHH[202] + 1./4.*aryuHH[18] + aryuHH[209] + 2*
   aryuHH[19] + aryuHH[207] + aryuHH[208] - 2*aryuHH[61] + 1./4.*
   aryuHH[5] + aryuHH[34];
   aryuHH[202]=aryuHH[3]*MMZ*aryuHH[202];
   aryuHH[210]=4*aryuHH[29];
   aryuHH[211]=2*aryuHH[9];
   aryuHH[212]=3./2.*aryuHH[30];
   aryuHH[213]=aryuHH[211] + aryuHH[212] - 13./2. + aryuHH[210];
   aryuHH[213]=aryuHH[19]*aryuHH[213];
   aryuHH[193]= - 47./4. + aryuHH[193];
   aryuHH[193]= - 9./16.*aryuHH[27] + aryuHH[195] + 1./4.*aryuHH[193]
    + aryuHH[197];
   aryuHH[193]=aryuHH[18]*aryuHH[193];
   aryuHH[195]=5*aryuHH[41];
   aryuHH[214]=5./2.*aryuHH[37];
   aryuHH[215]=aryuHH[214] - 41./8. + aryuHH[195];
   aryuHH[215]=3./2.*aryuHH[215] + 5*aryuHH[12];
   aryuHH[215]=aryuHH[28]*aryuHH[215];
   aryuHH[216]= - 11./8.*aryuHH[5] - 7*aryuHH[34];
   aryuHH[216]=111./16.*aryuHH[20] + 15*aryuHH[21] + 111./8.*aryuHH[45]
    - 51./8.*aryuHH[33] + 13./2.*aryuHH[67] + 3*aryuHH[216] + 17*
   aryuHH[61];
   aryuHH[217]=aryuHH[11]*aryuHH[28];
   aryuHH[218]=aryuHH[27]*aryuHH[28];
   aryuHH[110]=3*aryuHH[202] + aryuHH[110] + 19./4.*aryuHH[217] + 1./2.
   *aryuHH[193] + 9./16.*aryuHH[218] + aryuHH[215] + 1./2.*aryuHH[216]
    + aryuHH[213];
   aryuHH[110]=aryuHH[3]*aryuHH[110];
   aryuHH[193]=5./2.*aryuHH[20];
   aryuHH[202]=1./2.*aryuHH[33];
   aryuHH[213]=aryuHH[193] + aryuHH[202] - 15*aryuHH[45];
   aryuHH[213]=aryuHH[97]*aryuHH[213];
   aryuHH[215]=1./2.*aryuHH[22];
   aryuHH[216]=aryuHH[215] - 241./32. - 3*aryuHH[59];
   aryuHH[216]=3*aryuHH[216] - 13*aryuHH[62];
   aryuHH[219]=3*aryuHH[56];
   aryuHH[220]=3./8.*aryuHH[55];
   aryuHH[221]=1./2.*aryuHH[29];
   aryuHH[222]= - 3*aryuHH[16];
   aryuHH[213]=3./16.*aryuHH[213] + aryuHH[221] - 9./8.*aryuHH[15] - 21.
   /32.*aryuHH[63] + aryuHH[222] - 35./16.*aryuHH[60] + aryuHH[220] + 7.
   /8.*aryuHH[13] + 207./32.*aryuHH[43] - 63./32.*aryuHH[58] + 1./2.*
   aryuHH[216] + aryuHH[219];
   aryuHH[216]= - 3./2.*aryuHH[60];
   aryuHH[223]=aryuHH[216] + aryuHH[158] + aryuHH[55];
   aryuHH[223]=aryuHH[188] + aryuHH[221] + 1./2.*aryuHH[223] - 
   aryuHH[15];
   aryuHH[223]=aryuHH[99]*aryuHH[223];
   aryuHH[224]= - 3*aryuHH[43] + 5./2. + 3*aryuHH[58];
   aryuHH[225]=3./2.*aryuHH[60];
   aryuHH[224]=aryuHH[15] - 7./4.*aryuHH[63] + aryuHH[225] + 1./4.*
   aryuHH[224] - aryuHH[55];
   aryuHH[224]=aryuHH[97]*aryuHH[224];
   aryuHH[226]=aryuHH[30] - aryuHH[16] - 3./2.*aryuHH[62] + aryuHH[56];
   aryuHH[226]=aryuHH[98]*aryuHH[226];
   aryuHH[227]= - aryuHH[12]*aryuHH[98]*aryuHH[30];
   aryuHH[228]= - aryuHH[29] - aryuHH[9];
   aryuHH[228]=aryuHH[99]*aryuHH[228];
   aryuHH[228]=1./2.*aryuHH[97] + aryuHH[228];
   aryuHH[228]=aryuHH[11]*aryuHH[228];
   aryuHH[223]=1./4.*aryuHH[228] + 3./32.*aryuHH[204] + aryuHH[227] + 1.
   /2.*aryuHH[223] + 1./8.*aryuHH[224] + aryuHH[226];
   aryuHH[223]=MMZ*aryuHH[223];
   aryuHH[224]= - 9./2. - aryuHH[62];
   aryuHH[226]= - 1./2.*aryuHH[63];
   aryuHH[194]=aryuHH[226] + aryuHH[194] + 1./2.*aryuHH[224] - 
   aryuHH[64];
   aryuHH[194]=MMZ*aryuHH[194];
   aryuHH[224]=1./2.*aryuHH[34] - aryuHH[61];
   aryuHH[227]=1./2.*aryuHH[18];
   aryuHH[228]=2*aryuHH[28];
   aryuHH[194]=aryuHH[194] + aryuHH[227] + aryuHH[228] + aryuHH[19] + 
   aryuHH[207] + aryuHH[224] + aryuHH[208];
   aryuHH[194]=aryuHH[3]*aryuHH[194];
   aryuHH[207]= - 1./2.*aryuHH[30];
   aryuHH[229]= - 3./2.*aryuHH[56];
   aryuHH[230]=3./2.*aryuHH[16];
   aryuHH[231]=1./2.*aryuHH[64];
   aryuHH[232]=aryuHH[207] - 1./4.*aryuHH[29] - 3./4.*aryuHH[15] + 23./
   8.*aryuHH[63] + aryuHH[230] + 9./8.*aryuHH[60] + 3./4.*aryuHH[55] + 
   aryuHH[231] + 15./8.*aryuHH[43] + 1./2.*aryuHH[37] + aryuHH[41] - 7./
   8.*aryuHH[58] + aryuHH[229] + 9./2.*aryuHH[62] + 43./8. - aryuHH[59]
   ;
   aryuHH[195]=aryuHH[214] + 3./8. + aryuHH[195];
   aryuHH[195]=1./2.*aryuHH[195] + 2*aryuHH[12];
   aryuHH[195]=aryuHH[27]*aryuHH[195];
   aryuHH[214]=7 + 23./2.*aryuHH[29];
   aryuHH[233]=3*aryuHH[27];
   aryuHH[214]=1./4.*aryuHH[214] + aryuHH[233];
   aryuHH[214]=aryuHH[11]*aryuHH[214];
   aryuHH[234]= - 3./4.*aryuHH[30];
   aryuHH[235]=aryuHH[234] + 8 + 13./2.*aryuHH[29];
   aryuHH[235]=aryuHH[12]*aryuHH[235];
   aryuHH[236]=2*aryuHH[51] + aryuHH[50];
   aryuHH[236]=MMZ*aryuHH[236];
   aryuHH[194]=3*aryuHH[194] + 3*aryuHH[236] + aryuHH[214] + 3*
   aryuHH[195] + 3./2.*aryuHH[232] + aryuHH[235];
   aryuHH[194]=aryuHH[3]*aryuHH[194];
   aryuHH[195]=3./2.*aryuHH[62];
   aryuHH[214]= - aryuHH[30] + aryuHH[16] - aryuHH[64] - aryuHH[56] - 1
    + aryuHH[195];
   aryuHH[214]=aryuHH[98]*aryuHH[214];
   aryuHH[232]= - aryuHH[29] + aryuHH[15] - aryuHH[63] + aryuHH[225] - 
   1 - aryuHH[55];
   aryuHH[232]=aryuHH[99]*aryuHH[232];
   aryuHH[235]=aryuHH[11]*aryuHH[99]*aryuHH[29];
   aryuHH[236]=aryuHH[12]*aryuHH[98]*aryuHH[30];
   aryuHH[237]=1./4.*aryuHH[236];
   aryuHH[214]=1./8.*aryuHH[235] + aryuHH[237] + 1./8.*aryuHH[232] + 1./
   4.*aryuHH[214] - aryuHH[51] - 1./2.*aryuHH[50];
   aryuHH[232]=1./2.*aryuHH[60];
   aryuHH[238]=aryuHH[226] + aryuHH[232] - 1./2.*aryuHH[64] + 1./2. + 
   aryuHH[62];
   aryuHH[238]=aryuHH[3]*aryuHH[238];
   aryuHH[239]=1./2.*aryuHH[53];
   aryuHH[238]=aryuHH[238] - aryuHH[50] - aryuHH[51] + aryuHH[239];
   aryuHH[238]=MMt*aryuHH[3]*aryuHH[238];
   aryuHH[194]=3*aryuHH[238] + 3*aryuHH[214] + aryuHH[194];
   aryuHH[194]=MMt*aryuHH[194];
   aryuHH[214]= - 13./2. - 19./3.*aryuHH[29];
   aryuHH[238]= - 11./12.*aryuHH[9];
   aryuHH[214]=3./2.*aryuHH[206] + aryuHH[238] + 1./4.*aryuHH[214] + 
   aryuHH[30];
   aryuHH[214]=aryuHH[11]*aryuHH[214];
   aryuHH[240]=aryuHH[202] + aryuHH[160] - aryuHH[67];
   aryuHH[240]=1./2.*aryuHH[240] - aryuHH[20];
   aryuHH[240]=aryuHH[99]*aryuHH[240];
   aryuHH[241]= - 1 - 2./3.*aryuHH[29];
   aryuHH[242]=1./2.*aryuHH[30];
   aryuHH[241]= - 2./3.*aryuHH[9] + 2*aryuHH[241] + aryuHH[242];
   aryuHH[241]=aryuHH[12]*aryuHH[241];
   aryuHH[243]= - aryuHH[9] + 5./2. - aryuHH[29];
   aryuHH[243]=aryuHH[99]*aryuHH[243];
   aryuHH[243]=3./4.*aryuHH[199] - 7./8.*aryuHH[97] + aryuHH[243];
   aryuHH[243]=aryuHH[18]*aryuHH[243];
   aryuHH[244]= - aryuHH[15] + aryuHH[63] - 1 + aryuHH[55];
   aryuHH[244]=aryuHH[97]*aryuHH[244];
   aryuHH[200]=aryuHH[244] + aryuHH[200];
   aryuHH[200]=MMH*aryuHH[200];
   aryuHH[224]=aryuHH[224] - aryuHH[21];
   aryuHH[224]=aryuHH[98]*aryuHH[224];
   aryuHH[244]=1./4.*aryuHH[9];
   aryuHH[245]=1 + aryuHH[207];
   aryuHH[245]=aryuHH[19]*aryuHH[98]*aryuHH[245];
   aryuHH[246]=1./2.*aryuHH[99] + 13./16.*aryuHH[97] + aryuHH[98];
   aryuHH[246]=aryuHH[28]*aryuHH[246];
   aryuHH[247]=3./4.*aryuHH[206] - 1./2.*aryuHH[37] - 3./16. - 
   aryuHH[41];
   aryuHH[247]=aryuHH[27]*aryuHH[247];
   aryuHH[110]=aryuHH[194] + aryuHH[110] + 3./64.*aryuHH[200] + 3./2.*
   aryuHH[223] + 1./2.*aryuHH[214] + 3./8.*aryuHH[243] + 3./4.*
   aryuHH[247] + 3./2.*aryuHH[246] + 3*aryuHH[245] + aryuHH[241] + 3./4.
   *aryuHH[240] + aryuHH[244] + 3./2.*aryuHH[224] + 1./2.*aryuHH[213]
    + aryuHH[30];
   aryuHH[110]=aryuHH[26]*aryuHH[110];
   aryuHH[194]= - 55./2.*aryuHH[39] + 125./3. + 209./4.*aryuHH[40];
   aryuHH[194]=1./2.*aryuHH[194] + aryuHH[38];
   aryuHH[125]=aryuHH[131] + aryuHH[130] + aryuHH[146] + aryuHH[129] + 
   aryuHH[125] + 1./2.*aryuHH[194] - aryuHH[36];
   aryuHH[125]=aryuHH[11]*aryuHH[125];
   aryuHH[129]= - 3*aryuHH[84];
   aryuHH[130]=19./4. + aryuHH[129];
   aryuHH[131]= - 2*aryuHH[92];
   aryuHH[130]= - 11./8.*aryuHH[88] + 1./2.*aryuHH[130] + aryuHH[131];
   aryuHH[194]=9./2.*aryuHH[42];
   aryuHH[213]=1./2.*aryuHH[38];
   aryuHH[214]=9./2.*aryuHH[7];
   aryuHH[130]=aryuHH[214] + aryuHH[213] - 33./8.*aryuHH[39] + 25./8.*
   aryuHH[16] + 33./16.*aryuHH[90] + aryuHH[194] + 33./4.*aryuHH[91] + 
   3*aryuHH[130] + aryuHH[87];
   aryuHH[130]=aryuHH[98]*aryuHH[130];
   aryuHH[223]=19 - 11*aryuHH[86];
   aryuHH[223]=1./4.*aryuHH[223] - 3*aryuHH[83];
   aryuHH[223]=1./4.*aryuHH[223] - aryuHH[89];
   aryuHH[224]=1./4.*aryuHH[36];
   aryuHH[240]=9./4.*aryuHH[7];
   aryuHH[241]=1./2.*aryuHH[85];
   aryuHH[223]=aryuHH[240] + aryuHH[224] - 33./16.*aryuHH[40] + 25./16.
   *aryuHH[15] + 99./32.*aryuHH[90] + 9./4.*aryuHH[42] + 33./16.*
   aryuHH[91] + 3*aryuHH[223] + aryuHH[241];
   aryuHH[223]=aryuHH[99]*aryuHH[223];
   aryuHH[243]=1./2.*aryuHH[74];
   aryuHH[245]=aryuHH[75] + aryuHH[243];
   aryuHH[246]= - aryuHH[38] - 1 + 33./4.*aryuHH[39];
   aryuHH[246]=aryuHH[12]*aryuHH[98]*aryuHH[246];
   aryuHH[247]= - aryuHH[36] - 1 + 33./4.*aryuHH[40];
   aryuHH[247]=aryuHH[11]*aryuHH[99]*aryuHH[247];
   aryuHH[130]=3./4.*aryuHH[247] + 3./2.*aryuHH[246] + 3*aryuHH[223] + 
   3*aryuHH[130] - aryuHH[76] + 39./4.*aryuHH[79] + 39./2.*aryuHH[77]
    + 9*aryuHH[245] - 2*aryuHH[78];
   aryuHH[130]=MMZ*aryuHH[130];
   aryuHH[223]=aryuHH[95] - 1./2.*aryuHH[72];
   aryuHH[223]=33./4.*aryuHH[223] - aryuHH[93];
   aryuHH[223]=aryuHH[122] - 35./4.*aryuHH[20] + aryuHH[119] + 3*
   aryuHH[223] + aryuHH[117];
   aryuHH[223]=aryuHH[99]*aryuHH[223];
   aryuHH[245]=99./4.*aryuHH[95];
   aryuHH[246]= - 99./8.*aryuHH[72];
   aryuHH[247]= - 35./4.*aryuHH[21] + aryuHH[119] + aryuHH[246] + 
   aryuHH[245] - 3*aryuHH[94] + 7./4.*aryuHH[71];
   aryuHH[247]=aryuHH[98]*aryuHH[247];
   aryuHH[248]=aryuHH[17]*aryuHH[98];
   aryuHH[146]=aryuHH[146] + 3./4.*aryuHH[248] + aryuHH[107] - 
   aryuHH[36] - 77./8.*aryuHH[39] + 125./6. + 22*aryuHH[40];
   aryuHH[146]=aryuHH[12]*aryuHH[146];
   aryuHH[249]= - 13 + 9*aryuHH[39];
   aryuHH[249]=aryuHH[127] + 11./4.*aryuHH[249] + aryuHH[154];
   aryuHH[249]=aryuHH[98]*aryuHH[249];
   aryuHH[250]= - aryuHH[12]*aryuHH[98];
   aryuHH[249]=3*aryuHH[250] + aryuHH[249] - 99./4.*aryuHH[99];
   aryuHH[249]=aryuHH[19]*aryuHH[249];
   aryuHH[251]= - 1 + 9./4.*aryuHH[40];
   aryuHH[251]=aryuHH[127] + 11*aryuHH[251] + aryuHH[109];
   aryuHH[251]=aryuHH[99]*aryuHH[251];
   aryuHH[251]= - 99./2.*aryuHH[98] + aryuHH[251];
   aryuHH[251]=aryuHH[18]*aryuHH[251];
   aryuHH[252]= - aryuHH[75] - 1./2.*aryuHH[74];
   aryuHH[253]= - 1./2.*aryuHH[76];
   aryuHH[252]=aryuHH[253] - 11./8.*aryuHH[79] - 11./4.*aryuHH[77] + 9*
   aryuHH[252] - aryuHH[78];
   aryuHH[252]=MMH*aryuHH[252];
   aryuHH[254]=6287./24. - 27*aryuHH[84];
   aryuHH[254]=29./4.*aryuHH[88] - 27./4.*aryuHH[83] + 29./8.*
   aryuHH[86] + 1./2.*aryuHH[254] - 13./3.*aryuHH[92];
   aryuHH[125]=1./2.*aryuHH[252] + aryuHH[130] + aryuHH[125] + 1./4.*
   aryuHH[251] + 1./2.*aryuHH[249] + aryuHH[146] + 1./4.*aryuHH[223] + 
   5./8.*aryuHH[248] + 1./2.*aryuHH[247] + 81./8.*aryuHH[7] + 
   aryuHH[224] + aryuHH[213] - 33./4.*aryuHH[39] - 33./8.*aryuHH[40] - 
   37./16.*aryuHH[15] - 37./8.*aryuHH[16] + 85./4.*aryuHH[90] + 81./8.*
   aryuHH[42] + 85./2.*aryuHH[91] + aryuHH[241] - 13./12.*aryuHH[89] + 
   1./2.*aryuHH[254] + aryuHH[87];
   aryuHH[125]=aryuHH[68]*aryuHH[125];
   aryuHH[130]=99./8.*aryuHH[90] + 99./4.*aryuHH[91] - aryuHH[89] + 273.
   /8. + aryuHH[131];
   aryuHH[130]=MMZ*aryuHH[130];
   aryuHH[131]=1./2.*aryuHH[70];
   aryuHH[146]=aryuHH[12]*aryuHH[17];
   aryuHH[130]=aryuHH[130] + 1./2.*aryuHH[152] - 91./4.*aryuHH[18] - 91.
   /2.*aryuHH[19] + aryuHH[146] + aryuHH[151] + aryuHH[138] + 
   aryuHH[131] - aryuHH[93] + aryuHH[246] + aryuHH[245] - 2*aryuHH[94]
    + aryuHH[71];
   aryuHH[130]=aryuHH[68]*MMZ*aryuHH[130];
   aryuHH[138]=1./2.*aryuHH[24];
   aryuHH[151]=1./2.*aryuHH[13];
   aryuHH[223]=aryuHH[151] + aryuHH[138] + aryuHH[14] + 3./2. + 
   aryuHH[25];
   aryuHH[224]=aryuHH[31]*aryuHH[223];
   aryuHH[223]=aryuHH[1]*aryuHH[223];
   aryuHH[223]=3*aryuHH[224] + aryuHH[223];
   aryuHH[223]=MMZ*aryuHH[223];
   aryuHH[224]=aryuHH[6] + aryuHH[160];
   aryuHH[245]=aryuHH[31]*aryuHH[224];
   aryuHH[224]=aryuHH[1]*aryuHH[224];
   aryuHH[246]= - 3*aryuHH[31] - aryuHH[1];
   aryuHH[247]=aryuHH[19]*aryuHH[246];
   aryuHH[246]=aryuHH[18]*aryuHH[246];
   aryuHH[223]=aryuHH[223] + 1./2.*aryuHH[246] + aryuHH[247] + 3*
   aryuHH[245] + aryuHH[224];
   aryuHH[223]=MMZ*aryuHH[223];
   aryuHH[130]=aryuHH[223] + 3*aryuHH[130];
   aryuHH[130]=aryuHH[3]*aryuHH[130];
   aryuHH[223]=15./4.*aryuHH[15];
   aryuHH[224]=3./2.*aryuHH[10];
   aryuHH[245]= - 1./2.*aryuHH[25];
   aryuHH[246]=aryuHH[245] - 9./4. - 5*aryuHH[23];
   aryuHH[111]=aryuHH[191] + aryuHH[224] + aryuHH[223] + 15./2.*
   aryuHH[16] + aryuHH[111] - 15./4.*aryuHH[22] - 3./8.*aryuHH[24] + 3./
   2.*aryuHH[246] - 5*aryuHH[14];
   aryuHH[111]=aryuHH[31]*aryuHH[111];
   aryuHH[191]=aryuHH[244] + aryuHH[139] + 5./4.*aryuHH[15] + 5./2.*
   aryuHH[16] - 5./6.*aryuHH[13] - 5./4.*aryuHH[22] - 1./8.*aryuHH[24]
    + 1./2.*aryuHH[246] - 5./3.*aryuHH[14];
   aryuHH[191]=aryuHH[1]*aryuHH[191];
   aryuHH[244]=3*aryuHH[9];
   aryuHH[246]=aryuHH[244] + 2 + aryuHH[180];
   aryuHH[247]=aryuHH[31]*aryuHH[246];
   aryuHH[246]=aryuHH[1]*aryuHH[246];
   aryuHH[246]=3*aryuHH[247] + aryuHH[246];
   aryuHH[246]=aryuHH[12]*aryuHH[246];
   aryuHH[247]=7./4.*aryuHH[9];
   aryuHH[249]=1 + aryuHH[180];
   aryuHH[251]=aryuHH[249] + aryuHH[247];
   aryuHH[252]=aryuHH[31]*aryuHH[251];
   aryuHH[251]=aryuHH[1]*aryuHH[251];
   aryuHH[251]=3*aryuHH[252] + aryuHH[251];
   aryuHH[251]=aryuHH[11]*aryuHH[251];
   aryuHH[111]=aryuHH[251] + aryuHH[246] + aryuHH[111] + aryuHH[191];
   aryuHH[111]=MMZ*aryuHH[111];
   aryuHH[191]= - 11./4.*aryuHH[40];
   aryuHH[246]= - 1 + aryuHH[191];
   aryuHH[251]=11./2.*aryuHH[39];
   aryuHH[252]=5*aryuHH[246] + aryuHH[251];
   aryuHH[252]= - 5*aryuHH[12] + aryuHH[114] + 5./2.*aryuHH[252] - 
   aryuHH[38];
   aryuHH[252]=3*aryuHH[252] + aryuHH[175];
   aryuHH[252]=aryuHH[11]*aryuHH[252];
   aryuHH[254]= - aryuHH[38] + 77./4.*aryuHH[39] - 25 - 121./2.*
   aryuHH[40];
   aryuHH[254]=1./2.*aryuHH[254] + aryuHH[36];
   aryuHH[254]=3*aryuHH[254] - 13./4.*aryuHH[12];
   aryuHH[254]=aryuHH[12]*aryuHH[254];
   aryuHH[255]=1./2.*aryuHH[76];
   aryuHH[256]=aryuHH[255] - 33./4.*aryuHH[79] + aryuHH[78] - 33./2.*
   aryuHH[77];
   aryuHH[256]=MMZ*aryuHH[256];
   aryuHH[140]=3*aryuHH[256] + 1./2.*aryuHH[252] + aryuHH[254] + 
   aryuHH[140] + 3./2.*aryuHH[38] - 99./8.*aryuHH[39] - 99./16.*
   aryuHH[40] - 207./16.*aryuHH[15] - 207./8.*aryuHH[16] - 1635./32.*
   aryuHH[90] - 1635./16.*aryuHH[91] + aryuHH[149] + aryuHH[148] - 3*
   aryuHH[87] + 231./8.*aryuHH[88] + 231./16.*aryuHH[86] - 707./8. + 8./
   3.*aryuHH[92];
   aryuHH[140]=MMZ*aryuHH[140];
   aryuHH[148]=3 + aryuHH[191];
   aryuHH[148]=5*aryuHH[148] + aryuHH[251];
   aryuHH[148]= - 13./2.*aryuHH[12] + aryuHH[114] + 5./2.*aryuHH[148]
    - aryuHH[38];
   aryuHH[148]=aryuHH[18]*aryuHH[148];
   aryuHH[149]=13*aryuHH[146] + 41./4.*aryuHH[17] - 2059./12.*
   aryuHH[20] - 2059./6.*aryuHH[21] + 77./4.*aryuHH[44] - 37./4.*
   aryuHH[70] + 17./6.*aryuHH[93] + 1575./8.*aryuHH[72] - 423./2.*
   aryuHH[95] + 17./3.*aryuHH[94] - 37./2.*aryuHH[71];
   aryuHH[191]= - 41./4.*aryuHH[12] + aryuHH[213] + 11./8.*aryuHH[39]
    + 75./2. - 22*aryuHH[40];
   aryuHH[191]=aryuHH[19]*aryuHH[191];
   aryuHH[252]=aryuHH[17] - 15*aryuHH[19];
   aryuHH[252]=13*aryuHH[252] + 33*aryuHH[18];
   aryuHH[252]=aryuHH[11]*aryuHH[252];
   aryuHH[140]=aryuHH[140] + 1./8.*aryuHH[252] + 3./2.*aryuHH[148] + 1./
   4.*aryuHH[149] + 3*aryuHH[191];
   aryuHH[140]=aryuHH[68]*aryuHH[140];
   aryuHH[148]= - aryuHH[6] - 1./2.*aryuHH[5];
   aryuHH[148]=aryuHH[193] + 11./2.*aryuHH[148] + 5*aryuHH[21];
   aryuHH[149]=aryuHH[31]*aryuHH[148];
   aryuHH[148]=aryuHH[1]*aryuHH[148];
   aryuHH[148]=3*aryuHH[149] + aryuHH[148];
   aryuHH[149]=aryuHH[211] + 4./3. + aryuHH[139];
   aryuHH[149]=aryuHH[1]*aryuHH[149];
   aryuHH[191]=6*aryuHH[9] + 4 + aryuHH[224];
   aryuHH[191]=aryuHH[31]*aryuHH[191];
   aryuHH[149]=aryuHH[191] + aryuHH[149];
   aryuHH[149]=aryuHH[19]*aryuHH[149];
   aryuHH[183]=aryuHH[183] + aryuHH[247];
   aryuHH[183]=aryuHH[1]*aryuHH[183];
   aryuHH[191]=21./4.*aryuHH[9] + 2 + aryuHH[147];
   aryuHH[191]=aryuHH[31]*aryuHH[191];
   aryuHH[183]=aryuHH[191] + aryuHH[183];
   aryuHH[183]=aryuHH[18]*aryuHH[183];
   aryuHH[111]=aryuHH[130] + aryuHH[140] + aryuHH[111] + aryuHH[183] + 
   1./2.*aryuHH[148] + aryuHH[149];
   aryuHH[111]=aryuHH[3]*aryuHH[111];
   aryuHH[130]=aryuHH[40] - aryuHH[39];
   aryuHH[140]=33./4.*aryuHH[130];
   aryuHH[148]= - aryuHH[36] + aryuHH[140] + aryuHH[38];
   aryuHH[149]=aryuHH[12]*aryuHH[148];
   aryuHH[183]=aryuHH[11]*aryuHH[148];
   aryuHH[149]=aryuHH[149] + 1./2.*aryuHH[183];
   aryuHH[149]=MMZ*aryuHH[149];
   aryuHH[183]=3*aryuHH[148];
   aryuHH[191]=aryuHH[183] + 17./4.*aryuHH[12];
   aryuHH[191]=aryuHH[19]*aryuHH[191];
   aryuHH[183]=aryuHH[183] - 17./2.*aryuHH[12];
   aryuHH[183]=aryuHH[18]*aryuHH[183];
   aryuHH[149]=3*aryuHH[149] + 17./8.*aryuHH[174] + aryuHH[191] + 1./2.
   *aryuHH[183];
   aryuHH[149]=aryuHH[68]*aryuHH[149];
   aryuHH[183]=aryuHH[10] - aryuHH[9];
   aryuHH[191]=aryuHH[31]*aryuHH[183];
   aryuHH[193]=aryuHH[1]*aryuHH[183];
   aryuHH[211]=3*aryuHH[191] + aryuHH[193];
   aryuHH[224]=aryuHH[12]*aryuHH[211];
   aryuHH[247]=aryuHH[11]*aryuHH[211];
   aryuHH[224]=aryuHH[224] + 1./2.*aryuHH[247];
   aryuHH[224]=MMZ*aryuHH[224];
   aryuHH[247]=aryuHH[19]*aryuHH[211];
   aryuHH[211]=aryuHH[18]*aryuHH[211];
   aryuHH[149]=aryuHH[149] + aryuHH[224] + aryuHH[247] + 1./2.*
   aryuHH[211];
   aryuHH[149]=aryuHH[3]*aryuHH[149];
   aryuHH[211]= - 1./2.*aryuHH[29];
   aryuHH[224]= - 1./2.*aryuHH[9];
   aryuHH[247]=aryuHH[224] + aryuHH[211] + aryuHH[30];
   aryuHH[252]=aryuHH[12]*aryuHH[247];
   aryuHH[254]=1./2.*aryuHH[11]*aryuHH[247];
   aryuHH[256]=aryuHH[252] + aryuHH[254];
   aryuHH[256]=MMZ*aryuHH[256];
   aryuHH[257]=aryuHH[19]*aryuHH[247];
   aryuHH[258]=1./2.*aryuHH[18]*aryuHH[247];
   aryuHH[256]=aryuHH[256] + aryuHH[257] + aryuHH[258];
   aryuHH[256]=aryuHH[3]*aryuHH[256];
   aryuHH[259]=aryuHH[188] + aryuHH[221] - aryuHH[30];
   aryuHH[260]=aryuHH[12]*aryuHH[259];
   aryuHH[259]=aryuHH[11]*aryuHH[259];
   aryuHH[261]=aryuHH[29] - aryuHH[30];
   aryuHH[262]=aryuHH[12]*aryuHH[261];
   aryuHH[263]=aryuHH[11]*aryuHH[261];
   aryuHH[264]=aryuHH[262] + 1./2.*aryuHH[263];
   aryuHH[264]=MMt*aryuHH[3]*aryuHH[264];
   aryuHH[256]=3./2.*aryuHH[264] + 3*aryuHH[256] + aryuHH[260] + 1./2.*
   aryuHH[259];
   aryuHH[256]=aryuHH[26]*aryuHH[256];
   aryuHH[259]= - aryuHH[40] + aryuHH[39];
   aryuHH[260]=2*aryuHH[36];
   aryuHH[264]= - 2*aryuHH[38];
   aryuHH[265]=aryuHH[260] + 33./4.*aryuHH[259] + aryuHH[264];
   aryuHH[265]=aryuHH[12]*aryuHH[265];
   aryuHH[266]=aryuHH[36] + 33./8.*aryuHH[259] - aryuHH[38];
   aryuHH[266]=aryuHH[11]*aryuHH[266];
   aryuHH[265]=aryuHH[265] + aryuHH[266];
   aryuHH[265]=aryuHH[68]*aryuHH[265];
   aryuHH[266]= - aryuHH[10] + aryuHH[9];
   aryuHH[267]=aryuHH[31]*aryuHH[266];
   aryuHH[268]=aryuHH[1]*aryuHH[266];
   aryuHH[269]=aryuHH[267] + 1./3.*aryuHH[268];
   aryuHH[270]=aryuHH[12]*aryuHH[269];
   aryuHH[269]=aryuHH[11]*aryuHH[269];
   aryuHH[149]=aryuHH[256] + aryuHH[149] + aryuHH[265] + aryuHH[270] + 
   1./2.*aryuHH[269];
   aryuHH[256]=pow(aryuHH[4],2);
   aryuHH[149]=aryuHH[256]*aryuHH[149];
   aryuHH[167]=aryuHH[31]*aryuHH[167];
   aryuHH[167]=3*aryuHH[167] + aryuHH[168];
   aryuHH[167]=aryuHH[99]*aryuHH[167];
   aryuHH[168]= - aryuHH[98]*aryuHH[10];
   aryuHH[265]=aryuHH[31]*aryuHH[168];
   aryuHH[168]=aryuHH[1]*aryuHH[168];
   aryuHH[168]=3*aryuHH[265] + aryuHH[168];
   aryuHH[168]=aryuHH[12]*aryuHH[168];
   aryuHH[169]=3*aryuHH[169] + aryuHH[170];
   aryuHH[169]=aryuHH[11]*aryuHH[99]*aryuHH[169];
   aryuHH[170]=aryuHH[10] - aryuHH[16] - 3./2.*aryuHH[14] + aryuHH[245]
    - 3./2. + aryuHH[23];
   aryuHH[170]=aryuHH[98]*aryuHH[170];
   aryuHH[245]=aryuHH[31]*aryuHH[170];
   aryuHH[170]=aryuHH[1]*aryuHH[170];
   aryuHH[167]=1./2.*aryuHH[169] + aryuHH[168] + 1./2.*aryuHH[167] + 3*
   aryuHH[245] + aryuHH[170];
   aryuHH[167]=MMZ*aryuHH[167];
   aryuHH[168]= - 3./2.*aryuHH[15] + aryuHH[222] + 7./4.*aryuHH[13] + 3.
   /2.*aryuHH[22] + 7./2.*aryuHH[14] - 55./16. + 3*aryuHH[23];
   aryuHH[169]=1./2.*aryuHH[6] - aryuHH[21];
   aryuHH[169]=aryuHH[98]*aryuHH[169];
   aryuHH[168]=aryuHH[188] + 3./2.*aryuHH[169] + 1./2.*aryuHH[168] + 
   aryuHH[10];
   aryuHH[168]=aryuHH[31]*aryuHH[168];
   aryuHH[170]= - 1./2.*aryuHH[15];
   aryuHH[215]=aryuHH[170] - aryuHH[16] + 7./12.*aryuHH[13] + 
   aryuHH[215] + 7./6.*aryuHH[14] - 55./48. + aryuHH[23];
   aryuHH[169]=aryuHH[156] + 1./2.*aryuHH[169] + 1./2.*aryuHH[215] + 1./
   3.*aryuHH[10];
   aryuHH[169]=aryuHH[1]*aryuHH[169];
   aryuHH[161]=3*aryuHH[162] + aryuHH[161];
   aryuHH[161]=aryuHH[99]*aryuHH[161];
   aryuHH[162]= - 2*aryuHH[9] - 2 + aryuHH[139];
   aryuHH[215]=aryuHH[31]*aryuHH[162];
   aryuHH[162]=aryuHH[1]*aryuHH[162];
   aryuHH[162]=aryuHH[215] + 1./3.*aryuHH[162];
   aryuHH[162]=aryuHH[12]*aryuHH[162];
   aryuHH[215]=1./2. - aryuHH[10];
   aryuHH[215]=aryuHH[98]*aryuHH[215];
   aryuHH[245]=aryuHH[31]*aryuHH[215];
   aryuHH[215]=aryuHH[1]*aryuHH[215];
   aryuHH[215]=3*aryuHH[245] + aryuHH[215];
   aryuHH[215]=aryuHH[19]*aryuHH[215];
   aryuHH[164]=3*aryuHH[165] + aryuHH[164];
   aryuHH[164]=aryuHH[18]*aryuHH[99]*aryuHH[164];
   aryuHH[139]= - 5./4.*aryuHH[9] - 1 + aryuHH[139];
   aryuHH[165]=aryuHH[31]*aryuHH[139];
   aryuHH[139]=aryuHH[1]*aryuHH[139];
   aryuHH[139]=aryuHH[165] + 1./3.*aryuHH[139];
   aryuHH[139]=aryuHH[11]*aryuHH[139];
   aryuHH[110]=aryuHH[149] + aryuHH[110] + aryuHH[111] + aryuHH[125] + 
   1./2.*aryuHH[167] + aryuHH[139] + 1./4.*aryuHH[164] + 1./2.*
   aryuHH[215] + aryuHH[162] + 1./4.*aryuHH[161] + aryuHH[168] + 
   aryuHH[169];
   aryuHH[110]=aryuHH[256]*aryuHH[110];
   aryuHH[100]=aryuHH[100] - 79 - 73./2.*aryuHH[40];
   aryuHH[100]=aryuHH[106] + aryuHH[105] + aryuHH[108] + aryuHH[103] + 
   aryuHH[107] + aryuHH[109] + 1./6.*aryuHH[100] + aryuHH[38];
   aryuHH[100]=aryuHH[11]*aryuHH[100];
   aryuHH[103]= - 7./2. + 11*aryuHH[86];
   aryuHH[103]=1./4.*aryuHH[103] - 9*aryuHH[83];
   aryuHH[105]= - 19./8.*aryuHH[15];
   aryuHH[103]=aryuHH[214] + aryuHH[135] + 11./8.*aryuHH[40] + 
   aryuHH[105] - 33./16.*aryuHH[90] + aryuHH[194] - 55./16.*aryuHH[91]
    + aryuHH[85] + 1./2.*aryuHH[103] + aryuHH[113];
   aryuHH[103]=aryuHH[99]*aryuHH[103];
   aryuHH[106]= - 125./8. + 9*aryuHH[84];
   aryuHH[111]=55./8.*aryuHH[39];
   aryuHH[113]= - 9./2.*aryuHH[42];
   aryuHH[125]= - 9./2.*aryuHH[7];
   aryuHH[139]= - 1./2.*aryuHH[38];
   aryuHH[106]=aryuHH[125] + aryuHH[139] + aryuHH[111] - 31./8.*
   aryuHH[16] - 11./8.*aryuHH[90] + aryuHH[113] - 55./4.*aryuHH[91] - 
   aryuHH[87] + 55./8.*aryuHH[88] + 1./2.*aryuHH[106] + 6*aryuHH[92];
   aryuHH[106]=aryuHH[98]*aryuHH[106];
   aryuHH[149]= - 1 - 5./4.*aryuHH[39];
   aryuHH[149]=11*aryuHH[149] + aryuHH[38];
   aryuHH[149]=aryuHH[98]*aryuHH[149];
   aryuHH[161]=2*aryuHH[12]*aryuHH[98];
   aryuHH[149]=1./2.*aryuHH[149] + aryuHH[161];
   aryuHH[149]=aryuHH[12]*aryuHH[149];
   aryuHH[162]=aryuHH[246] - aryuHH[36];
   aryuHH[162]=aryuHH[11]*aryuHH[99]*aryuHH[162];
   aryuHH[103]=3./2.*aryuHH[162] + 3*aryuHH[149] + 3*aryuHH[103] + 3*
   aryuHH[106] + aryuHH[120] - 67./4.*aryuHH[79] - 89./4.*aryuHH[77] + 
   2*aryuHH[78] + aryuHH[141] + 10*aryuHH[80] + 8*aryuHH[81] - 9*
   aryuHH[75];
   aryuHH[103]=MMZ*aryuHH[103];
   aryuHH[106]=11./4.*aryuHH[115] - aryuHH[93];
   aryuHH[106]=aryuHH[122] - 167./4.*aryuHH[20] + aryuHH[119] + 3*
   aryuHH[106] + aryuHH[117];
   aryuHH[106]=aryuHH[99]*aryuHH[106];
   aryuHH[117]=5 + aryuHH[118];
   aryuHH[117]=aryuHH[127] + 11*aryuHH[117] + aryuHH[109];
   aryuHH[117]=aryuHH[99]*aryuHH[117];
   aryuHH[117]=33./2.*aryuHH[98] + aryuHH[117];
   aryuHH[117]=aryuHH[18]*aryuHH[117];
   aryuHH[118]= - 785./2. - 27*aryuHH[86];
   aryuHH[112]= - aryuHH[88] + 1./2.*aryuHH[118] + aryuHH[112];
   aryuHH[112]=1./2.*aryuHH[112] + aryuHH[124];
   aryuHH[118]=25 - 11*aryuHH[39];
   aryuHH[118]=aryuHH[98]*aryuHH[118];
   aryuHH[118]=aryuHH[118] + 11*aryuHH[99];
   aryuHH[118]=1./4.*aryuHH[118] + aryuHH[161];
   aryuHH[118]=aryuHH[19]*aryuHH[118];
   aryuHH[119]= - 2*aryuHH[81];
   aryuHH[120]=aryuHH[253] + 17./8.*aryuHH[79] - 3./4.*aryuHH[77] - 9./
   2.*aryuHH[74] + aryuHH[119] - 3*aryuHH[80];
   aryuHH[120]=MMH*aryuHH[120];
   aryuHH[122]=11*aryuHH[115] - 3*aryuHH[21];
   aryuHH[122]=aryuHH[98]*aryuHH[122];
   aryuHH[124]= - 43./4.*aryuHH[12] + 7./2.*aryuHH[39] - 91 - 53*
   aryuHH[40];
   aryuHH[124]=aryuHH[12]*aryuHH[124];
   aryuHH[100]=aryuHH[120] + aryuHH[103] + aryuHH[100] + 1./2.*
   aryuHH[117] + 3*aryuHH[118] + 1./3.*aryuHH[124] + 1./2.*aryuHH[106]
    + 3./4.*aryuHH[122] + aryuHH[107] + aryuHH[135] + aryuHH[251] + 11./
   4.*aryuHH[40] + 19./8.*aryuHH[15] + 9./4.*aryuHH[16] - 89./8.*
   aryuHH[90] + aryuHH[133] - 251./8.*aryuHH[91] + 1./2.*aryuHH[112] + 
   aryuHH[85];
   aryuHH[100]=aryuHH[68]*aryuHH[100];
   aryuHH[103]= - 805./12. - aryuHH[24];
   aryuHH[103]=1./2.*aryuHH[103] + aryuHH[157];
   aryuHH[106]= - 5./3.*aryuHH[13];
   aryuHH[107]= - 5./4.*aryuHH[55];
   aryuHH[112]= - 17./24.*aryuHH[60];
   aryuHH[117]= - 75./16.*aryuHH[63];
   aryuHH[118]= - 5*aryuHH[45];
   aryuHH[120]=aryuHH[118] + 1./2.*aryuHH[20];
   aryuHH[120]=15./4.*aryuHH[97]*aryuHH[120];
   aryuHH[103]=aryuHH[188] + aryuHH[197] + aryuHH[120] + aryuHH[211] + 
   aryuHH[223] + aryuHH[117] - 15*aryuHH[16] + aryuHH[112] + 
   aryuHH[107] + aryuHH[106] - 15./16.*aryuHH[43] + 15./16.*aryuHH[58]
    + 15*aryuHH[56] + 1./2.*aryuHH[103] - 11*aryuHH[62];
   aryuHH[122]= - 8*aryuHH[29];
   aryuHH[124]=215./24. + aryuHH[122];
   aryuHH[133]=5./6.*aryuHH[9];
   aryuHH[124]=15./2.*aryuHH[206] + aryuHH[133] + 1./3.*aryuHH[124] + 
   aryuHH[201];
   aryuHH[124]=aryuHH[11]*aryuHH[124];
   aryuHH[141]= - aryuHH[43] + 3 + aryuHH[58];
   aryuHH[141]=aryuHH[15] - 5./4.*aryuHH[63] + aryuHH[232] + 1./4.*
   aryuHH[141] - aryuHH[55];
   aryuHH[141]=aryuHH[97]*aryuHH[141];
   aryuHH[149]=aryuHH[11]*aryuHH[97];
   aryuHH[141]=aryuHH[149] + aryuHH[141] + 1./4.*aryuHH[204];
   aryuHH[141]=15./8.*MMZ*aryuHH[141];
   aryuHH[161]= - 17./3. - 14*aryuHH[29];
   aryuHH[161]= - 10./3.*aryuHH[9] + 2./3.*aryuHH[161] + aryuHH[201];
   aryuHH[161]=aryuHH[12]*aryuHH[161];
   aryuHH[162]= - 1./4. + aryuHH[206];
   aryuHH[162]=aryuHH[27]*aryuHH[162];
   aryuHH[164]= - aryuHH[97] + aryuHH[199];
   aryuHH[164]=aryuHH[18]*aryuHH[164];
   aryuHH[103]=aryuHH[141] + aryuHH[124] + 15./16.*aryuHH[164] + 15./8.
   *aryuHH[162] + 75./8.*aryuHH[196] + 1./2.*aryuHH[103] + aryuHH[161];
   aryuHH[103]=MMZ*aryuHH[103];
   aryuHH[161]=aryuHH[232] + aryuHH[151] + 6*aryuHH[62] + 7 + 
   aryuHH[138];
   aryuHH[161]=MMZ*aryuHH[161];
   aryuHH[165]= - 3./2.*aryuHH[18];
   aryuHH[167]= - 6*aryuHH[19];
   aryuHH[168]= - 1./2.*aryuHH[33];
   aryuHH[161]=aryuHH[161] + aryuHH[165] - 8*aryuHH[28] + aryuHH[167]
    + aryuHH[168] + aryuHH[67] + 6*aryuHH[61] + aryuHH[160] - 3*
   aryuHH[34];
   aryuHH[161]=aryuHH[3]*MMZ*aryuHH[161];
   aryuHH[122]=385./48. + aryuHH[122];
   aryuHH[122]=15./16.*aryuHH[27] + aryuHH[133] + 1./3.*aryuHH[122] + 
   aryuHH[201];
   aryuHH[122]=aryuHH[18]*aryuHH[122];
   aryuHH[169]=15./8.*aryuHH[20] - 247./12.*aryuHH[45] + 17./4.*
   aryuHH[33] - 11./4.*aryuHH[5] - 13./3.*aryuHH[67];
   aryuHH[194]= - 4*aryuHH[29];
   aryuHH[215]= - aryuHH[9] + 5./3. + aryuHH[194];
   aryuHH[215]=aryuHH[19]*aryuHH[215];
   aryuHH[245]=55./12. - aryuHH[37];
   aryuHH[246]=5./2.*aryuHH[245] - 56./3.*aryuHH[12];
   aryuHH[246]=aryuHH[28]*aryuHH[246];
   aryuHH[251]= - aryuHH[11]*aryuHH[28];
   aryuHH[253]=71./6.*aryuHH[251];
   aryuHH[265]= - aryuHH[27]*aryuHH[28];
   aryuHH[103]=aryuHH[161] + aryuHH[103] + aryuHH[253] + aryuHH[122] + 
   15./8.*aryuHH[265] + aryuHH[246] + 1./2.*aryuHH[169] + 4./3.*
   aryuHH[215];
   aryuHH[103]=aryuHH[3]*aryuHH[103];
   aryuHH[161]= - 3*aryuHH[58];
   aryuHH[215]=3*aryuHH[43] - 5./2. + aryuHH[161];
   aryuHH[215]= - aryuHH[15] + 7./4.*aryuHH[63] + aryuHH[216] + 1./4.*
   aryuHH[215] + aryuHH[55];
   aryuHH[215]=aryuHH[97]*aryuHH[215];
   aryuHH[158]=aryuHH[9] - aryuHH[29] + aryuHH[225] + aryuHH[158] - 
   aryuHH[55];
   aryuHH[158]=aryuHH[99]*aryuHH[158];
   aryuHH[195]= - aryuHH[30] + aryuHH[16] + aryuHH[195] - aryuHH[56];
   aryuHH[195]=aryuHH[98]*aryuHH[195];
   aryuHH[225]=aryuHH[29] - aryuHH[9];
   aryuHH[225]=aryuHH[99]*aryuHH[225];
   aryuHH[225]= - 5./2.*aryuHH[97] + aryuHH[225];
   aryuHH[225]=aryuHH[11]*aryuHH[225];
   aryuHH[195]=1./2.*aryuHH[225] + 15./16.*aryuHH[199] + 3*aryuHH[236]
    + 1./2.*aryuHH[158] + 5./4.*aryuHH[215] + 3*aryuHH[195];
   aryuHH[195]=MMZ*aryuHH[195];
   aryuHH[246]= - 59./2. - 43*aryuHH[29];
   aryuHH[269]= - 2*aryuHH[27];
   aryuHH[270]=3./4.*aryuHH[30];
   aryuHH[246]=aryuHH[269] + 1./3.*aryuHH[246] + aryuHH[270];
   aryuHH[246]=aryuHH[11]*aryuHH[246];
   aryuHH[271]= - 2*aryuHH[28];
   aryuHH[272]=aryuHH[67] + aryuHH[168];
   aryuHH[273]= - aryuHH[18] + aryuHH[272] + aryuHH[271];
   aryuHH[273]=13*aryuHH[273];
   aryuHH[274]=aryuHH[64] + 4 + 1./2.*aryuHH[62];
   aryuHH[274]=aryuHH[63] + 3*aryuHH[274] + 13./2.*aryuHH[60];
   aryuHH[274]=MMZ*aryuHH[274];
   aryuHH[274]=aryuHH[273] + aryuHH[274];
   aryuHH[274]=aryuHH[3]*aryuHH[274];
   aryuHH[275]= - 307./3. + 19*aryuHH[58];
   aryuHH[275]= - 11./2.*aryuHH[29] - 33./2.*aryuHH[15] - 83./4.*
   aryuHH[63] - 77./4.*aryuHH[60] + 33./2.*aryuHH[55] - 19./4.*
   aryuHH[43] + 1./4.*aryuHH[275] - aryuHH[37];
   aryuHH[276]= - 1 - aryuHH[29];
   aryuHH[276]=aryuHH[12]*aryuHH[276];
   aryuHH[277]= - 3./4. - aryuHH[37];
   aryuHH[277]=aryuHH[27]*aryuHH[277];
   aryuHH[278]= - 3*aryuHH[51];
   aryuHH[279]=aryuHH[278] - aryuHH[50];
   aryuHH[279]=MMZ*aryuHH[279];
   aryuHH[274]=aryuHH[274] + 2*aryuHH[279] + aryuHH[246] + 5./2.*
   aryuHH[277] + 1./2.*aryuHH[275] + 56./3.*aryuHH[276];
   aryuHH[274]=aryuHH[3]*aryuHH[274];
   aryuHH[276]= - 11*aryuHH[29] + 11*aryuHH[15] + aryuHH[63] + 33./2.*
   aryuHH[60] + 1 - 11*aryuHH[55];
   aryuHH[276]=aryuHH[99]*aryuHH[276];
   aryuHH[235]=11./4.*aryuHH[235] + 5*aryuHH[50] + 1./4.*aryuHH[276];
   aryuHH[276]=13*aryuHH[63] + 24 + 11*aryuHH[60];
   aryuHH[276]=aryuHH[3]*aryuHH[276];
   aryuHH[276]=aryuHH[276] - aryuHH[53] - 22*aryuHH[50];
   aryuHH[276]=MMt*aryuHH[3]*aryuHH[276];
   aryuHH[274]=aryuHH[276] + aryuHH[235] + aryuHH[274];
   aryuHH[274]=MMt*aryuHH[274];
   aryuHH[279]=3./2.*aryuHH[55];
   aryuHH[280]= - 5./2.*aryuHH[15];
   aryuHH[116]=aryuHH[280] + 35./8.*aryuHH[63] - 1./12.*aryuHH[60] + 
   aryuHH[279] + aryuHH[116] - 9./8.*aryuHH[43] + 9./8.*aryuHH[58] + 43.
   /48. + aryuHH[22];
   aryuHH[160]=aryuHH[168] + aryuHH[160] + aryuHH[67];
   aryuHH[281]=aryuHH[99]*aryuHH[160];
   aryuHH[282]= - 5./2.*aryuHH[20] + aryuHH[168] + 15*aryuHH[45];
   aryuHH[282]=aryuHH[97]*aryuHH[282];
   aryuHH[116]=1./2.*aryuHH[281] + aryuHH[159] + 5./8.*aryuHH[282] + 1./
   2.*aryuHH[116] - 1./3.*aryuHH[29];
   aryuHH[159]= - aryuHH[9] - 3./2. + aryuHH[29];
   aryuHH[159]=aryuHH[99]*aryuHH[159];
   aryuHH[159]=15./4.*aryuHH[204] + 35./8.*aryuHH[97] + aryuHH[159];
   aryuHH[159]=aryuHH[18]*aryuHH[159];
   aryuHH[281]= - 215./12. + 13*aryuHH[29];
   aryuHH[281]=1./9.*aryuHH[281] + aryuHH[30];
   aryuHH[281]=5./2.*aryuHH[196] + 1./2.*aryuHH[281] - 1./9.*aryuHH[9];
   aryuHH[281]=aryuHH[11]*aryuHH[281];
   aryuHH[282]=aryuHH[15] - aryuHH[63] + 1 - aryuHH[55];
   aryuHH[282]=aryuHH[97]*aryuHH[282];
   aryuHH[149]=aryuHH[282] + aryuHH[149];
   aryuHH[149]=5./32.*MMH*aryuHH[149];
   aryuHH[210]=aryuHH[9] - 5./3. + aryuHH[210];
   aryuHH[210]=aryuHH[12]*aryuHH[210];
   aryuHH[282]= - 65./8.*aryuHH[97] - aryuHH[99];
   aryuHH[282]=aryuHH[28]*aryuHH[282];
   aryuHH[283]=15./2.*aryuHH[196] + 15./8. + aryuHH[37];
   aryuHH[283]=aryuHH[27]*aryuHH[283];
   aryuHH[103]=aryuHH[274] + aryuHH[103] + aryuHH[149] + 1./2.*
   aryuHH[195] + aryuHH[281] + 1./4.*aryuHH[159] + 1./4.*aryuHH[283] + 
   1./2.*aryuHH[282] + 1./2.*aryuHH[116] + 4./9.*aryuHH[210];
   aryuHH[103]=aryuHH[26]*aryuHH[103];
   aryuHH[195]=53*aryuHH[40];
   aryuHH[274]=aryuHH[154] + aryuHH[145] + 59 + aryuHH[195];
   aryuHH[175]=aryuHH[175] + 25./2.*aryuHH[12] + 1./2.*aryuHH[274] + 
   aryuHH[177];
   aryuHH[175]=aryuHH[11]*aryuHH[175];
   aryuHH[274]=713./2. - 8*aryuHH[92];
   aryuHH[284]= - 3./2.*aryuHH[38];
   aryuHH[285]=35./3.*aryuHH[12] + aryuHH[284] + 71./8.*aryuHH[39] + 
   737./6. + 119*aryuHH[40];
   aryuHH[285]=aryuHH[12]*aryuHH[285];
   aryuHH[119]=aryuHH[119] - aryuHH[80];
   aryuHH[119]=2*aryuHH[119] - aryuHH[78];
   aryuHH[119]=aryuHH[76] + 22*aryuHH[79] + 2*aryuHH[119] + 55./2.*
   aryuHH[77];
   aryuHH[119]=MMZ*aryuHH[119];
   aryuHH[114]=3*aryuHH[119] + aryuHH[175] + aryuHH[285] + aryuHH[114]
    + aryuHH[284] + 165./8.*aryuHH[39] + 33./8.*aryuHH[40] + 145./8.*
   aryuHH[15] + 269./8.*aryuHH[16] + 953./16.*aryuHH[90] + 2317./16.*
   aryuHH[91] + aryuHH[182] + aryuHH[181] + 3*aryuHH[87] - 341./8.*
   aryuHH[88] + 1./3.*aryuHH[274] - 121./8.*aryuHH[86];
   aryuHH[114]=MMZ*aryuHH[114];
   aryuHH[119]=aryuHH[154] + aryuHH[145] - 105 + aryuHH[195];
   aryuHH[119]=155./4.*aryuHH[12] + 1./2.*aryuHH[119] + aryuHH[177];
   aryuHH[119]=aryuHH[18]*aryuHH[119];
   aryuHH[175]=439./12.*aryuHH[12] + aryuHH[145] + 55./6. + aryuHH[195]
   ;
   aryuHH[175]=aryuHH[19]*aryuHH[175];
   aryuHH[114]=aryuHH[114] + aryuHH[143] + aryuHH[119] + aryuHH[175] + 
   41./24.*aryuHH[17] + 1589./24.*aryuHH[20] + 445./12.*aryuHH[21] + 77.
   /24.*aryuHH[44] - 37./8.*aryuHH[70] + 17./12.*aryuHH[93] - 525./16.*
   aryuHH[72] - 8./3.*aryuHH[96] + 141./4.*aryuHH[95];
   aryuHH[114]=aryuHH[68]*aryuHH[114];
   aryuHH[119]=1./2.*aryuHH[25];
   aryuHH[143]=aryuHH[119] + 3./2. + 5*aryuHH[23];
   aryuHH[175]=aryuHH[147] - 15./2.*aryuHH[16] + 3./2.*aryuHH[143] + 5*
   aryuHH[14];
   aryuHH[175]=aryuHH[31]*aryuHH[175];
   aryuHH[143]=aryuHH[180] - 5./2.*aryuHH[16] + 1./2.*aryuHH[143] + 5./
   3.*aryuHH[14];
   aryuHH[143]=aryuHH[1]*aryuHH[143];
   aryuHH[147]= - 38./3.*aryuHH[9] - 34./9. + aryuHH[147];
   aryuHH[147]=aryuHH[31]*aryuHH[147];
   aryuHH[181]= - 6*aryuHH[9] - 2./3. + aryuHH[180];
   aryuHH[181]=aryuHH[1]*aryuHH[181];
   aryuHH[147]=aryuHH[147] + aryuHH[181];
   aryuHH[147]=aryuHH[12]*aryuHH[147];
   aryuHH[143]=aryuHH[186] + aryuHH[147] + aryuHH[175] + aryuHH[143];
   aryuHH[143]=MMZ*aryuHH[143];
   aryuHH[147]= - aryuHH[14] - 1 - aryuHH[25];
   aryuHH[175]=aryuHH[31]*aryuHH[147];
   aryuHH[147]=aryuHH[1]*aryuHH[147];
   aryuHH[147]=3*aryuHH[175] + aryuHH[147];
   aryuHH[147]=MMZ*aryuHH[147];
   aryuHH[175]= - aryuHH[31]*aryuHH[6];
   aryuHH[181]= - aryuHH[1]*aryuHH[6];
   aryuHH[182]=3*aryuHH[31] + aryuHH[1];
   aryuHH[182]=aryuHH[19]*aryuHH[182];
   aryuHH[147]=2*aryuHH[147] + aryuHH[182] + 3*aryuHH[175] + 
   aryuHH[181];
   aryuHH[147]=MMZ*aryuHH[147];
   aryuHH[175]= - 33./2.*aryuHH[90] - 231./4.*aryuHH[91] - 2*aryuHH[89]
    - 289./4. + 4*aryuHH[92];
   aryuHH[175]=MMZ*aryuHH[175];
   aryuHH[181]= - aryuHH[12]*aryuHH[17];
   aryuHH[152]=aryuHH[175] + aryuHH[152] + 37*aryuHH[18] + 62*
   aryuHH[19] + aryuHH[181] + aryuHH[70] - 2*aryuHH[93] + 33./2.*
   aryuHH[72] - 33*aryuHH[95] + 2*aryuHH[94] - aryuHH[71];
   aryuHH[152]=aryuHH[68]*MMZ*aryuHH[152];
   aryuHH[147]=aryuHH[147] + 3*aryuHH[152];
   aryuHH[147]=aryuHH[3]*aryuHH[147];
   aryuHH[152]=1./3. - aryuHH[9];
   aryuHH[175]=aryuHH[31]*aryuHH[152];
   aryuHH[152]=aryuHH[1]*aryuHH[152];
   aryuHH[152]=5./3.*aryuHH[175] + aryuHH[152];
   aryuHH[152]=aryuHH[19]*aryuHH[152];
   aryuHH[114]=aryuHH[147] + aryuHH[114] + aryuHH[143] + 4*aryuHH[152]
    + aryuHH[185];
   aryuHH[114]=aryuHH[3]*aryuHH[114];
   aryuHH[143]=aryuHH[98]*aryuHH[10];
   aryuHH[147]=aryuHH[31]*aryuHH[143];
   aryuHH[143]=aryuHH[1]*aryuHH[143];
   aryuHH[143]=3*aryuHH[147] + aryuHH[143];
   aryuHH[143]=aryuHH[12]*aryuHH[143];
   aryuHH[119]= - aryuHH[10] + aryuHH[16] + 3./2.*aryuHH[14] + 
   aryuHH[119] + 3./2. - aryuHH[23];
   aryuHH[119]=aryuHH[98]*aryuHH[119];
   aryuHH[147]=aryuHH[31]*aryuHH[119];
   aryuHH[119]=aryuHH[1]*aryuHH[119];
   aryuHH[119]=aryuHH[143] + 3*aryuHH[147] + aryuHH[119];
   aryuHH[119]=MMZ*aryuHH[119];
   aryuHH[143]= - 1./3. + aryuHH[9];
   aryuHH[147]=aryuHH[31]*aryuHH[143];
   aryuHH[143]=aryuHH[1]*aryuHH[143];
   aryuHH[143]=5./3.*aryuHH[147] + aryuHH[143];
   aryuHH[143]=aryuHH[12]*aryuHH[143];
   aryuHH[100]=aryuHH[110] + aryuHH[103] + aryuHH[114] + aryuHH[100] + 
   1./2.*aryuHH[119] + 4./3.*aryuHH[143] + aryuHH[187];
   aryuHH[100]=aryuHH[256]*aryuHH[100];
   aryuHH[103]=13./12.*aryuHH[39];
   aryuHH[110]= - 7./12.*aryuHH[38];
   aryuHH[114]=7./6.*aryuHH[36];
   aryuHH[119]= - 21./4.*aryuHH[7];
   aryuHH[147]= - 7./4.*aryuHH[12];
   aryuHH[152]= - aryuHH[11] + aryuHH[147] + aryuHH[119] + aryuHH[114]
    + aryuHH[110] + aryuHH[103] + 3 - 1./48.*aryuHH[40];
   aryuHH[152]=aryuHH[11]*aryuHH[152];
   aryuHH[175]=aryuHH[145] + 19 + aryuHH[123];
   aryuHH[175]=1./2.*aryuHH[175] + aryuHH[38];
   aryuHH[182]= - 1./2.*aryuHH[36];
   aryuHH[175]=aryuHH[136] + aryuHH[119] + 1./3.*aryuHH[175] + 
   aryuHH[182];
   aryuHH[175]=aryuHH[12]*aryuHH[175];
   aryuHH[185]= - 15./4.*aryuHH[88] + 21./4.*aryuHH[83] + 3./4.*
   aryuHH[86] - 1 + 21./4.*aryuHH[84];
   aryuHH[186]= - 3./8.*aryuHH[36] - 3./8.*aryuHH[38] + 13./8.*
   aryuHH[39] + 2 - 1./32.*aryuHH[40];
   aryuHH[186]=aryuHH[7]*aryuHH[186];
   aryuHH[187]=2*aryuHH[85];
   aryuHH[195]= - 21./4.*aryuHH[42];
   aryuHH[274]=3./8.*aryuHH[90];
   aryuHH[284]= - 1./8.*aryuHH[16];
   aryuHH[285]=aryuHH[76] - 1./2.*aryuHH[79] + aryuHH[142] + aryuHH[80]
    + aryuHH[78];
   aryuHH[285]=1./2.*MMH*aryuHH[285];
   aryuHH[286]=2*aryuHH[87];
   aryuHH[287]= - 1./12.*aryuHH[38];
   aryuHH[288]= - 1./12.*aryuHH[36];
   aryuHH[105]=aryuHH[285] + 1./2.*aryuHH[152] + 1./2.*aryuHH[175] + 3*
   aryuHH[186] + aryuHH[288] + aryuHH[287] + aryuHH[105] + aryuHH[284]
    + aryuHH[274] + aryuHH[195] - 15./8.*aryuHH[91] + aryuHH[187] + 1./
   2.*aryuHH[185] + aryuHH[286];
   aryuHH[105]=MMH*aryuHH[105];
   aryuHH[152]= - 521 - aryuHH[40];
   aryuHH[152]=1./4.*aryuHH[152] + aryuHH[145];
   aryuHH[175]=1./8.*aryuHH[248];
   aryuHH[185]= - 1./8.*aryuHH[19]*aryuHH[98];
   aryuHH[108]=aryuHH[185] + aryuHH[108] + aryuHH[175] + aryuHH[182] + 
   1./6.*aryuHH[152] + aryuHH[264];
   aryuHH[108]=aryuHH[19]*aryuHH[108];
   aryuHH[123]=aryuHH[109] + aryuHH[154] + aryuHH[145] + 67./3. + 
   aryuHH[123];
   aryuHH[123]=aryuHH[17]*aryuHH[123];
   aryuHH[152]= - 1./2.*aryuHH[40];
   aryuHH[186]= - 19 + aryuHH[152];
   aryuHH[186]=1./2.*aryuHH[186] + aryuHH[145];
   aryuHH[248]= - 9*aryuHH[36];
   aryuHH[289]=1./2.*aryuHH[102];
   aryuHH[186]=aryuHH[289] + aryuHH[248] + 1./3.*aryuHH[186] + 
   aryuHH[38];
   aryuHH[290]=1./4.*aryuHH[104];
   aryuHH[186]=aryuHH[290] + 1./2.*aryuHH[186] - 7*aryuHH[12];
   aryuHH[186]=aryuHH[18]*aryuHH[186];
   aryuHH[291]=19./3.*aryuHH[19];
   aryuHH[292]= - 3./2.*aryuHH[17];
   aryuHH[293]= - 53./6.*aryuHH[18] + aryuHH[292] + aryuHH[291];
   aryuHH[293]=1./4.*aryuHH[11]*aryuHH[293];
   aryuHH[294]=17./8.*aryuHH[71];
   aryuHH[295]=aryuHH[294] + 4./3.*aryuHH[96] - aryuHH[94];
   aryuHH[296]=17./8.*aryuHH[70];
   aryuHH[297]= - 9./4.*aryuHH[44];
   aryuHH[298]=1./2.*aryuHH[181];
   aryuHH[105]=aryuHH[105] + aryuHH[293] + 1./2.*aryuHH[186] + 
   aryuHH[108] + aryuHH[298] + 1./4.*aryuHH[123] - 241./24.*aryuHH[20]
    + 73./24.*aryuHH[21] + aryuHH[297] + aryuHH[296] - aryuHH[93] - 3./
   2.*aryuHH[72] + aryuHH[295] + 6*aryuHH[95];
   aryuHH[105]=aryuHH[68]*aryuHH[105];
   aryuHH[108]=101 + aryuHH[144];
   aryuHH[108]=1./8.*aryuHH[104] - 1./3.*aryuHH[12] + 1./8.*aryuHH[102]
    - 5./2.*aryuHH[36] + 1./12.*aryuHH[108] + aryuHH[38];
   aryuHH[108]=aryuHH[18]*aryuHH[108];
   aryuHH[123]= - 61./3.*aryuHH[95] + 13*aryuHH[72];
   aryuHH[144]= - 17*aryuHH[39];
   aryuHH[186]=7./3. + aryuHH[144];
   aryuHH[186]=aryuHH[17]*aryuHH[186];
   aryuHH[299]=139 + aryuHH[144];
   aryuHH[299]=1./4.*aryuHH[299] + aryuHH[12];
   aryuHH[299]=aryuHH[19]*aryuHH[299];
   aryuHH[291]= - aryuHH[17] + aryuHH[291];
   aryuHH[291]=1./2.*aryuHH[291] - aryuHH[18];
   aryuHH[291]=aryuHH[11]*aryuHH[291];
   aryuHH[108]=aryuHH[291] + aryuHH[108] + 1./3.*aryuHH[299] + 1./8.*
   aryuHH[186] - 17./3.*aryuHH[20] - 11./4.*aryuHH[21] - 9./8.*
   aryuHH[44] + aryuHH[296] + 1./4.*aryuHH[123] - aryuHH[93];
   aryuHH[123]= - 1 - 17./4.*aryuHH[39];
   aryuHH[186]=1./3.*aryuHH[123] + aryuHH[126];
   aryuHH[186]=aryuHH[12]*aryuHH[186];
   aryuHH[291]=1./4. + 3*aryuHH[83];
   aryuHH[291]=7./2.*aryuHH[291] + 15*aryuHH[88];
   aryuHH[123]=aryuHH[7]*aryuHH[123];
   aryuHH[299]= - 7./4.*aryuHH[38] + 1 - 17./16.*aryuHH[39];
   aryuHH[299]= - 1./2.*aryuHH[11] - 1./2.*aryuHH[12] - 21./8.*
   aryuHH[7] + 1./3.*aryuHH[299] + aryuHH[36];
   aryuHH[299]=aryuHH[11]*aryuHH[299];
   aryuHH[300]=1./2.*aryuHH[79] + aryuHH[76];
   aryuHH[300]=MMH*aryuHH[300];
   aryuHH[301]= - 1./24.*aryuHH[36];
   aryuHH[123]=1./4.*aryuHH[300] + 1./2.*aryuHH[299] + 1./4.*
   aryuHH[186] + 3./8.*aryuHH[123] + aryuHH[301] - aryuHH[15] - 15./8.*
   aryuHH[16] + 3./2.*aryuHH[90] - 21./16.*aryuHH[42] - 1./8.*
   aryuHH[91] + 1./8.*aryuHH[291] + aryuHH[85];
   aryuHH[123]=MMH*aryuHH[123];
   aryuHH[108]=1./2.*aryuHH[108] + aryuHH[123];
   aryuHH[108]=aryuHH[68]*aryuHH[108];
   aryuHH[123]= - aryuHH[7]*aryuHH[39];
   aryuHH[186]= - 1./12.*aryuHH[39];
   aryuHH[291]=aryuHH[186] - aryuHH[38];
   aryuHH[291]=aryuHH[11]*aryuHH[291];
   aryuHH[123]=aryuHH[291] + 1./6.*aryuHH[128] + 3./4.*aryuHH[123] - 
   aryuHH[16] + aryuHH[90] + 1./4. + aryuHH[88];
   aryuHH[123]=MMH*aryuHH[123];
   aryuHH[128]= - 1./2.*aryuHH[17];
   aryuHH[291]=1./3.*aryuHH[19];
   aryuHH[299]=1./6.*aryuHH[18] + aryuHH[128] + aryuHH[291];
   aryuHH[299]=aryuHH[11]*aryuHH[299];
   aryuHH[300]= - aryuHH[17]*aryuHH[39];
   aryuHH[302]=5 - 1./3.*aryuHH[39];
   aryuHH[302]=aryuHH[19]*aryuHH[302];
   aryuHH[186]=1 + aryuHH[186];
   aryuHH[186]=aryuHH[18]*aryuHH[186];
   aryuHH[115]=1./2.*aryuHH[123] + aryuHH[299] + aryuHH[186] + 1./2.*
   aryuHH[302] + 1./4.*aryuHH[300] + aryuHH[115] - 1./2.*aryuHH[21];
   aryuHH[115]=aryuHH[68]*aryuHH[115];
   aryuHH[123]=aryuHH[3]*aryuHH[68]*aryuHH[18]*aryuHH[173];
   aryuHH[115]=aryuHH[115] + 1./2.*aryuHH[123];
   aryuHH[115]=aryuHH[171]*aryuHH[115];
   aryuHH[123]= - aryuHH[31]*aryuHH[13];
   aryuHH[173]= - aryuHH[1]*aryuHH[13];
   aryuHH[186]=11./9.*aryuHH[31] + aryuHH[1];
   aryuHH[299]=aryuHH[11]*aryuHH[186];
   aryuHH[123]=aryuHH[299] + 11./9.*aryuHH[123] + aryuHH[173];
   aryuHH[123]=MMH*aryuHH[123];
   aryuHH[173]=aryuHH[18]*aryuHH[186];
   aryuHH[123]=aryuHH[173] + aryuHH[123];
   aryuHH[173]=1./2.*aryuHH[17];
   aryuHH[186]= - 43./2.*aryuHH[18] + aryuHH[173] - 91*aryuHH[19];
   aryuHH[186]=aryuHH[18]*aryuHH[186];
   aryuHH[299]=pow(aryuHH[19],2);
   aryuHH[186]= - aryuHH[299] + 1./8.*aryuHH[186];
   aryuHH[186]=aryuHH[3]*aryuHH[68]*aryuHH[186];
   aryuHH[108]=1./4.*aryuHH[115] + 1./3.*aryuHH[186] + 1./4.*
   aryuHH[123] + aryuHH[108];
   aryuHH[108]=aryuHH[171]*aryuHH[108];
   aryuHH[115]=5./3. - 9./4.*aryuHH[10];
   aryuHH[115]=aryuHH[7]*aryuHH[115];
   aryuHH[123]= - aryuHH[9]*aryuHH[7];
   aryuHH[115]=11./4.*aryuHH[123] - aryuHH[14] + aryuHH[115];
   aryuHH[115]=aryuHH[31]*aryuHH[115];
   aryuHH[186]= - 3./4.*aryuHH[10];
   aryuHH[300]=1 + aryuHH[186];
   aryuHH[300]=aryuHH[7]*aryuHH[300];
   aryuHH[300]=9./4.*aryuHH[123] - 1./3.*aryuHH[14] + aryuHH[300];
   aryuHH[300]=aryuHH[1]*aryuHH[300];
   aryuHH[302]=5./3. + aryuHH[180];
   aryuHH[224]=1./3.*aryuHH[302] + aryuHH[224];
   aryuHH[224]=aryuHH[1]*aryuHH[224];
   aryuHH[180]= - 11./18.*aryuHH[9] + 37./27. + aryuHH[180];
   aryuHH[180]=aryuHH[31]*aryuHH[180];
   aryuHH[180]=aryuHH[180] + aryuHH[224];
   aryuHH[180]=aryuHH[12]*aryuHH[180];
   aryuHH[224]= - 1./4.*aryuHH[10];
   aryuHH[302]=1./3. + aryuHH[224];
   aryuHH[172]=1./3.*aryuHH[302] + aryuHH[172];
   aryuHH[172]=aryuHH[1]*aryuHH[172];
   aryuHH[224]= - 11./36.*aryuHH[9] + 5./27. + aryuHH[224];
   aryuHH[224]=aryuHH[31]*aryuHH[224];
   aryuHH[172]=aryuHH[224] + aryuHH[172];
   aryuHH[224]=aryuHH[11]*aryuHH[172];
   aryuHH[115]=aryuHH[224] + aryuHH[180] + aryuHH[115] + aryuHH[300];
   aryuHH[115]=MMH*aryuHH[115];
   aryuHH[180]=aryuHH[238] + 5./9. + aryuHH[186];
   aryuHH[180]=aryuHH[31]*aryuHH[17]*aryuHH[180];
   aryuHH[176]=aryuHH[302] + aryuHH[176];
   aryuHH[176]=aryuHH[1]*aryuHH[17]*aryuHH[176];
   aryuHH[186]= - 11./9.*aryuHH[9] + 47./27. - aryuHH[10];
   aryuHH[186]=aryuHH[31]*aryuHH[186];
   aryuHH[224]=7./3. - aryuHH[10];
   aryuHH[224]=1./3.*aryuHH[224] - aryuHH[9];
   aryuHH[224]=aryuHH[1]*aryuHH[224];
   aryuHH[186]=aryuHH[186] + aryuHH[224];
   aryuHH[186]=aryuHH[19]*aryuHH[186];
   aryuHH[172]=aryuHH[18]*aryuHH[172];
   aryuHH[115]=1./2.*aryuHH[115] + aryuHH[172] + 1./2.*aryuHH[186] + 
   aryuHH[180] + aryuHH[176];
   aryuHH[172]=aryuHH[17] - 55*aryuHH[19];
   aryuHH[172]=aryuHH[19]*aryuHH[172];
   aryuHH[176]=53./24.*aryuHH[18] + 1./24.*aryuHH[17] + 6*aryuHH[19];
   aryuHH[176]=aryuHH[18]*aryuHH[176];
   aryuHH[172]=1./24.*aryuHH[172] + aryuHH[176];
   aryuHH[172]=aryuHH[3]*aryuHH[68]*aryuHH[172];
   aryuHH[105]=aryuHH[108] + aryuHH[172] + aryuHH[115] + aryuHH[105];
   aryuHH[105]=aryuHH[171]*aryuHH[105];
   aryuHH[103]= - aryuHH[11] + aryuHH[147] + aryuHH[119] + aryuHH[114]
    + aryuHH[110] + aryuHH[103] + 7 + 47./48.*aryuHH[40];
   aryuHH[103]=aryuHH[11]*aryuHH[103];
   aryuHH[108]=47./4.*aryuHH[40];
   aryuHH[110]=aryuHH[145] + 31 + aryuHH[108];
   aryuHH[110]=1./2.*aryuHH[110] + aryuHH[38];
   aryuHH[110]=aryuHH[136] + aryuHH[119] + 1./3.*aryuHH[110] + 
   aryuHH[182];
   aryuHH[110]=aryuHH[12]*aryuHH[110];
   aryuHH[114]=1 + 3./2.*aryuHH[84];
   aryuHH[172]=21./2.*aryuHH[83];
   aryuHH[114]= - 15./2.*aryuHH[88] + aryuHH[172] + 7*aryuHH[114] + 15./
   2.*aryuHH[86];
   aryuHH[180]= - 3./4.*aryuHH[36] - 3./4.*aryuHH[38] + 13./4.*
   aryuHH[39] + 7 + 47./16.*aryuHH[40];
   aryuHH[180]=aryuHH[7]*aryuHH[180];
   aryuHH[103]=aryuHH[285] + 1./2.*aryuHH[103] + 1./2.*aryuHH[110] + 3./
   2.*aryuHH[180] + aryuHH[288] + aryuHH[287] - 31./8.*aryuHH[15] + 
   aryuHH[284] + aryuHH[274] + aryuHH[195] - 3./8.*aryuHH[91] + 
   aryuHH[187] + 1./4.*aryuHH[114] + aryuHH[286];
   aryuHH[103]=MMH*aryuHH[103];
   aryuHH[110]= - 113 + 47*aryuHH[40];
   aryuHH[110]=1./4.*aryuHH[110] + aryuHH[145];
   aryuHH[110]=aryuHH[185] - aryuHH[12] + aryuHH[175] + aryuHH[182] + 1.
   /6.*aryuHH[110] + aryuHH[264];
   aryuHH[110]=aryuHH[19]*aryuHH[110];
   aryuHH[108]=aryuHH[109] + aryuHH[154] + aryuHH[145] + 103./3. + 
   aryuHH[108];
   aryuHH[108]=aryuHH[17]*aryuHH[108];
   aryuHH[114]=47./2.*aryuHH[40];
   aryuHH[180]=353 + aryuHH[114];
   aryuHH[180]=1./2.*aryuHH[180] + aryuHH[145];
   aryuHH[180]=aryuHH[289] + aryuHH[248] + 1./3.*aryuHH[180] + 
   aryuHH[38];
   aryuHH[180]=aryuHH[290] + 1./2.*aryuHH[180] - 15*aryuHH[12];
   aryuHH[180]=aryuHH[18]*aryuHH[180];
   aryuHH[103]=aryuHH[103] + aryuHH[293] + 1./2.*aryuHH[180] + 
   aryuHH[110] + aryuHH[298] + 1./4.*aryuHH[108] - 373./24.*aryuHH[20]
    + 145./24.*aryuHH[21] + aryuHH[297] + aryuHH[296] + aryuHH[295] - 
   aryuHH[93];
   aryuHH[103]=aryuHH[68]*aryuHH[103];
   aryuHH[108]= - 5*aryuHH[37] + 51 - 5*aryuHH[41];
   aryuHH[108]= - 55*aryuHH[39] + 1./2.*aryuHH[108] + 55*aryuHH[40];
   aryuHH[108]=2*aryuHH[191] + aryuHH[260] + 1./2.*aryuHH[108] + 2*
   aryuHH[38];
   aryuHH[110]=2*aryuHH[193];
   aryuHH[180]= - 3*aryuHH[12];
   aryuHH[108]= - 3./2.*aryuHH[27] + aryuHH[180] + 3*aryuHH[108] + 
   aryuHH[110];
   aryuHH[108]=aryuHH[27]*aryuHH[108];
   aryuHH[186]=2*aryuHH[61];
   aryuHH[187]= - aryuHH[34] + aryuHH[186];
   aryuHH[195]= - 2*aryuHH[19];
   aryuHH[168]= - aryuHH[18] - 4*aryuHH[28] + aryuHH[195] + aryuHH[168]
    + aryuHH[187] + aryuHH[67];
   aryuHH[168]=aryuHH[3]*aryuHH[168];
   aryuHH[224]= - 11 + 27./2.*aryuHH[59];
   aryuHH[238]=3./4.*aryuHH[56];
   aryuHH[248]= - 3./4.*aryuHH[41];
   aryuHH[274]= - 13./8.*aryuHH[64];
   aryuHH[284]= - 3./4.*aryuHH[16];
   aryuHH[285]= - aryuHH[12]*aryuHH[30];
   aryuHH[289]=9./4.*aryuHH[285];
   aryuHH[290]= - aryuHH[11]*aryuHH[30];
   aryuHH[293]=3./4.*aryuHH[290];
   aryuHH[295]=3*aryuHH[62];
   aryuHH[108]=3*aryuHH[168] + aryuHH[293] + aryuHH[108] + aryuHH[289]
    + aryuHH[234] - 13./4.*aryuHH[63] + aryuHH[284] + aryuHH[216] + 
   aryuHH[274] - 33./4.*aryuHH[43] - 3./4.*aryuHH[37] + aryuHH[248] + 3.
   /2.*aryuHH[58] + aryuHH[238] + 1./2.*aryuHH[224] + aryuHH[295];
   aryuHH[108]=aryuHH[3]*aryuHH[108];
   aryuHH[168]=1 + aryuHH[295];
   aryuHH[168]= - aryuHH[30] + aryuHH[16] + aryuHH[231] + 1./2.*
   aryuHH[168] - aryuHH[56];
   aryuHH[168]=aryuHH[98]*aryuHH[168];
   aryuHH[224]=3*aryuHH[51] + 1./4.*aryuHH[53];
   aryuHH[296]=1 + aryuHH[63];
   aryuHH[296]=aryuHH[99]*aryuHH[296];
   aryuHH[224]=aryuHH[237] + 1./8.*aryuHH[296] + 1./4.*aryuHH[168] + 1./
   2.*aryuHH[224] + aryuHH[50];
   aryuHH[237]=5 + aryuHH[62];
   aryuHH[237]=aryuHH[63] + 1./2.*aryuHH[237] + aryuHH[64];
   aryuHH[237]=aryuHH[3]*aryuHH[237];
   aryuHH[237]=aryuHH[237] - aryuHH[51] - 1./2.*aryuHH[53];
   aryuHH[237]=MMt*aryuHH[3]*aryuHH[237];
   aryuHH[108]=3*aryuHH[237] + 3*aryuHH[224] + aryuHH[108];
   aryuHH[108]=MMt*aryuHH[108];
   aryuHH[224]= - 7*aryuHH[62] - 25./4. + 3*aryuHH[59];
   aryuHH[224]= - 15./4.*aryuHH[43] + 9./4.*aryuHH[58] + 1./2.*
   aryuHH[224] + aryuHH[56];
   aryuHH[237]=aryuHH[272] + aryuHH[20];
   aryuHH[237]=aryuHH[99]*aryuHH[237];
   aryuHH[272]=3./4.*aryuHH[15];
   aryuHH[297]=1./8.*aryuHH[63];
   aryuHH[298]= - 3./2.*aryuHH[16];
   aryuHH[300]= - 1./2.*aryuHH[34];
   aryuHH[302]=aryuHH[300] + aryuHH[61];
   aryuHH[303]=aryuHH[302] + aryuHH[21];
   aryuHH[303]=3./2.*aryuHH[98]*aryuHH[303];
   aryuHH[216]=3./4.*aryuHH[237] + aryuHH[303] - aryuHH[30] + 
   aryuHH[211] + aryuHH[272] + aryuHH[297] + aryuHH[298] + aryuHH[216]
    - 3./4.*aryuHH[55] + 3./2.*aryuHH[224] + aryuHH[64];
   aryuHH[224]= - 33*aryuHH[39] + 33*aryuHH[40] - 3./2.*aryuHH[37] + 19
    - 3./2.*aryuHH[41];
   aryuHH[304]=6*aryuHH[38];
   aryuHH[305]= - 9./2.*aryuHH[12];
   aryuHH[110]=aryuHH[305] + aryuHH[110] + 6*aryuHH[191] + aryuHH[177]
    + 5./2.*aryuHH[224] + aryuHH[304];
   aryuHH[110]=aryuHH[28]*aryuHH[110];
   aryuHH[224]= - 11./4.*aryuHH[61];
   aryuHH[306]=7./2.*aryuHH[21];
   aryuHH[307]= - 9./4.*aryuHH[30];
   aryuHH[308]=aryuHH[307] - 7./2. + 8*aryuHH[29];
   aryuHH[308]=aryuHH[19]*aryuHH[308];
   aryuHH[309]=4*aryuHH[19] - 5./4.*aryuHH[28];
   aryuHH[309]=aryuHH[27]*aryuHH[309];
   aryuHH[310]=39./4.*aryuHH[27] + 1 + 23./4.*aryuHH[29];
   aryuHH[310]=aryuHH[18]*aryuHH[310];
   aryuHH[311]=3./4.*aryuHH[251];
   aryuHH[110]=aryuHH[311] + 1./2.*aryuHH[310] + 3*aryuHH[309] + 
   aryuHH[110] + aryuHH[308] + 5./8.*aryuHH[20] + aryuHH[306] - 41./4.*
   aryuHH[45] + 57./16.*aryuHH[33] - 5./2.*aryuHH[67] + 6*aryuHH[34] + 
   aryuHH[224];
   aryuHH[110]=aryuHH[3]*aryuHH[110];
   aryuHH[308]=55*aryuHH[39] - 55*aryuHH[40] + aryuHH[37] - 51./2. + 
   aryuHH[41];
   aryuHH[308]=aryuHH[267] - aryuHH[36] + 1./4.*aryuHH[308] - 
   aryuHH[38];
   aryuHH[308]=aryuHH[180] + 3*aryuHH[308] + aryuHH[268];
   aryuHH[308]=aryuHH[27]*aryuHH[308];
   aryuHH[194]= - 7 + aryuHH[194];
   aryuHH[309]=1./4.*aryuHH[30];
   aryuHH[194]=2./3.*aryuHH[194] + aryuHH[309];
   aryuHH[194]=aryuHH[12]*aryuHH[194];
   aryuHH[310]= - 19 - 29./2.*aryuHH[29];
   aryuHH[312]= - 3*aryuHH[27];
   aryuHH[310]=1./3.*aryuHH[310] + aryuHH[312];
   aryuHH[310]=aryuHH[11]*aryuHH[310];
   aryuHH[313]= - 1 + aryuHH[242];
   aryuHH[314]=aryuHH[19]*aryuHH[98]*aryuHH[313];
   aryuHH[315]=3./2.*aryuHH[314];
   aryuHH[316]= - aryuHH[98] - aryuHH[99];
   aryuHH[316]=aryuHH[28]*aryuHH[316];
   aryuHH[221]= - 1 + aryuHH[221];
   aryuHH[221]=aryuHH[18]*aryuHH[99]*aryuHH[221];
   aryuHH[108]=aryuHH[108] + aryuHH[110] + 1./4.*aryuHH[310] + 3./4.*
   aryuHH[221] + 1./2.*aryuHH[308] + 3./4.*aryuHH[316] + aryuHH[315] + 
   1./2.*aryuHH[216] + aryuHH[194];
   aryuHH[108]=MMt*aryuHH[108];
   aryuHH[110]=aryuHH[97]*aryuHH[45];
   aryuHH[194]=aryuHH[9]*aryuHH[7];
   aryuHH[216]=aryuHH[7]*aryuHH[29];
   aryuHH[308]= - aryuHH[30]*aryuHH[7];
   aryuHH[220]=3./2.*aryuHH[194] + 9*aryuHH[308] + 15./2.*aryuHH[216]
    + 3./2.*aryuHH[110] - 3./8.*aryuHH[15] + 3./8.*aryuHH[63] + 
   aryuHH[220] + 47./8. - aryuHH[13];
   aryuHH[310]=5./3.*aryuHH[29];
   aryuHH[316]=13./4. + aryuHH[310];
   aryuHH[316]=3./2.*aryuHH[196] + aryuHH[156] + 1./2.*aryuHH[316] - 
   aryuHH[30];
   aryuHH[316]=aryuHH[11]*aryuHH[316];
   aryuHH[317]=1./12.*aryuHH[9] + aryuHH[207] + 1 + 5./12.*aryuHH[29];
   aryuHH[317]=aryuHH[12]*aryuHH[317];
   aryuHH[220]=1./4.*aryuHH[316] + 3./8.*aryuHH[206] + 1./4.*
   aryuHH[220] + aryuHH[317];
   aryuHH[220]=MMH*aryuHH[220];
   aryuHH[316]=91./8. + aryuHH[310];
   aryuHH[316]= - 39./8.*aryuHH[27] + 3./4.*aryuHH[196] + aryuHH[156]
    + 1./2.*aryuHH[316] - aryuHH[30];
   aryuHH[316]=aryuHH[18]*aryuHH[316];
   aryuHH[317]=aryuHH[188] + 5./2.*aryuHH[29] + aryuHH[197];
   aryuHH[318]=aryuHH[27]*aryuHH[317];
   aryuHH[319]=6*aryuHH[30];
   aryuHH[320]= - aryuHH[9] - 5*aryuHH[29] + aryuHH[319];
   aryuHH[321]=aryuHH[3]*aryuHH[28]*aryuHH[320];
   aryuHH[320]=MMt*aryuHH[3]*aryuHH[27]*aryuHH[320];
   aryuHH[318]=aryuHH[320] + 1./2.*aryuHH[318] + aryuHH[321];
   aryuHH[318]=aryuHH[26]*MMt*aryuHH[318];
   aryuHH[317]=aryuHH[17]*aryuHH[317];
   aryuHH[310]=17 + aryuHH[310];
   aryuHH[156]=aryuHH[156] + 1./2.*aryuHH[310] - aryuHH[30];
   aryuHH[156]=aryuHH[19]*aryuHH[156];
   aryuHH[310]= - 5./3.*aryuHH[12];
   aryuHH[320]=3./16.*aryuHH[206] + 27./32. + aryuHH[310];
   aryuHH[320]=aryuHH[28]*aryuHH[320];
   aryuHH[321]=39*aryuHH[19] - 7./2.*aryuHH[28];
   aryuHH[321]=aryuHH[28]*aryuHH[321];
   aryuHH[322]=aryuHH[18]*aryuHH[28];
   aryuHH[321]=aryuHH[321] + 39./2.*aryuHH[322];
   aryuHH[321]=aryuHH[3]*aryuHH[321];
   aryuHH[323]= - 9./4.*aryuHH[21];
   aryuHH[324]= - aryuHH[19] + 5./16.*aryuHH[28];
   aryuHH[324]=aryuHH[27]*aryuHH[324];
   aryuHH[325]= - 1./4.*aryuHH[67];
   aryuHH[108]=aryuHH[318] + aryuHH[108] + 1./4.*aryuHH[321] + 1./2.*
   aryuHH[220] + 19./12.*aryuHH[251] + 1./4.*aryuHH[316] + 3*
   aryuHH[324] + aryuHH[320] + 1./2.*aryuHH[156] + 1./4.*aryuHH[317] - 
   57./64.*aryuHH[20] + aryuHH[323] + 63./32.*aryuHH[45] + 3./64.*
   aryuHH[33] - aryuHH[61] + aryuHH[325];
   aryuHH[108]=aryuHH[26]*aryuHH[108];
   aryuHH[156]= - 55./8.*aryuHH[40];
   aryuHH[220]=aryuHH[111] - 9 + aryuHH[156];
   aryuHH[316]=1./3.*aryuHH[38];
   aryuHH[317]= - 13./8.*aryuHH[12] + aryuHH[119] + aryuHH[220] + 
   aryuHH[316];
   aryuHH[317]=aryuHH[12]*aryuHH[317];
   aryuHH[318]=1./3.*aryuHH[36];
   aryuHH[119]= - aryuHH[11] - 5./2.*aryuHH[12] + aryuHH[119] + 
   aryuHH[220] + aryuHH[318];
   aryuHH[119]=aryuHH[11]*aryuHH[119];
   aryuHH[121]=aryuHH[255] + aryuHH[121] + aryuHH[78] + aryuHH[142];
   aryuHH[121]=MMH*aryuHH[121];
   aryuHH[142]=aryuHH[172] - 13*aryuHH[86] - 467./8. + 21*aryuHH[84];
   aryuHH[142]=1./2.*aryuHH[142] - 13*aryuHH[88];
   aryuHH[172]=5./8.*aryuHH[15];
   aryuHH[220]= - aryuHH[36] - aryuHH[38] + 55./4.*aryuHH[39] - 9 - 55./
   4.*aryuHH[40];
   aryuHH[220]=aryuHH[7]*aryuHH[220];
   aryuHH[119]=1./2.*aryuHH[121] + 1./4.*aryuHH[119] + 1./2.*
   aryuHH[317] + 9./8.*aryuHH[220] + aryuHH[301] + aryuHH[287] + 
   aryuHH[172] + 5./4.*aryuHH[16] - 17./8.*aryuHH[90] - 63./16.*
   aryuHH[42] - 17./4.*aryuHH[91] + aryuHH[85] + 1./4.*aryuHH[142] + 
   aryuHH[286];
   aryuHH[119]=MMH*aryuHH[119];
   aryuHH[111]=aryuHH[111] - 73 + aryuHH[156];
   aryuHH[121]=aryuHH[185] + 49./6.*aryuHH[12] + aryuHH[175] + 
   aryuHH[182] + aryuHH[111] + aryuHH[264];
   aryuHH[121]=aryuHH[19]*aryuHH[121];
   aryuHH[142]=165*aryuHH[39] - 89 - 165*aryuHH[40];
   aryuHH[142]=aryuHH[109] + 1./4.*aryuHH[142] + aryuHH[154];
   aryuHH[142]=aryuHH[17]*aryuHH[142];
   aryuHH[111]=aryuHH[111] + aryuHH[139];
   aryuHH[102]=1./16.*aryuHH[104] + 53./12.*aryuHH[12] + 1./16.*
   aryuHH[102] + 1./2.*aryuHH[111] - aryuHH[36];
   aryuHH[102]=aryuHH[18]*aryuHH[102];
   aryuHH[104]= - 41./24.*aryuHH[18] - 1./8.*aryuHH[17] + 8*aryuHH[19];
   aryuHH[104]=aryuHH[11]*aryuHH[104];
   aryuHH[111]= - 1./2.*aryuHH[93];
   aryuHH[102]=aryuHH[119] + aryuHH[104] + aryuHH[102] + aryuHH[121] + 
   1./4.*aryuHH[181] + 1./4.*aryuHH[142] + 79./24.*aryuHH[20] + 79./12.
   *aryuHH[21] - 27./16.*aryuHH[44] + 17./16.*aryuHH[70] + aryuHH[111]
    - 45./8.*aryuHH[72] + 171./8.*aryuHH[95] - aryuHH[94] + aryuHH[294]
   ;
   aryuHH[102]=aryuHH[68]*aryuHH[102];
   aryuHH[104]= - aryuHH[7]*aryuHH[29];
   aryuHH[119]=aryuHH[30]*aryuHH[7];
   aryuHH[104]=1./2.*aryuHH[123] + 1./2.*aryuHH[104] + aryuHH[119];
   aryuHH[104]=aryuHH[254] + 9./2.*aryuHH[104] + aryuHH[252];
   aryuHH[104]=MMH*aryuHH[104];
   aryuHH[121]=aryuHH[17]*aryuHH[247];
   aryuHH[104]=1./2.*aryuHH[104] + aryuHH[258] + 3./2.*aryuHH[121] + 
   aryuHH[257];
   aryuHH[121]=2*aryuHH[267] + aryuHH[260] + 33./2.*aryuHH[259] + 
   aryuHH[264];
   aryuHH[121]=3*aryuHH[121] + 2*aryuHH[268];
   aryuHH[142]= - 3./2.*aryuHH[12];
   aryuHH[156]=aryuHH[121] + aryuHH[142];
   aryuHH[156]=aryuHH[28]*aryuHH[156];
   aryuHH[175]=aryuHH[18]*aryuHH[261];
   aryuHH[181]=aryuHH[19]*aryuHH[261];
   aryuHH[156]=aryuHH[311] + 3./4.*aryuHH[175] + 3./2.*aryuHH[181] + 
   aryuHH[156];
   aryuHH[156]=aryuHH[3]*aryuHH[156];
   aryuHH[121]=aryuHH[27]*aryuHH[121];
   aryuHH[121]=aryuHH[293] + 3./2.*aryuHH[285] + aryuHH[121];
   aryuHH[121]=MMt*aryuHH[3]*aryuHH[121];
   aryuHH[175]=aryuHH[148] + aryuHH[191];
   aryuHH[175]=3*aryuHH[175] + aryuHH[193];
   aryuHH[175]=aryuHH[27]*aryuHH[175];
   aryuHH[185]= - aryuHH[29] + aryuHH[30];
   aryuHH[220]=aryuHH[12]*aryuHH[185];
   aryuHH[252]=aryuHH[11]*aryuHH[185];
   aryuHH[175]=1./2.*aryuHH[252] + aryuHH[220] + aryuHH[175];
   aryuHH[121]=aryuHH[121] + 1./2.*aryuHH[175] + aryuHH[156];
   aryuHH[121]=MMt*aryuHH[121];
   aryuHH[156]=aryuHH[27]*aryuHH[247];
   aryuHH[175]= - 2*aryuHH[30];
   aryuHH[220]=aryuHH[9] + aryuHH[29] + aryuHH[175];
   aryuHH[247]=aryuHH[3]*aryuHH[28]*aryuHH[220];
   aryuHH[220]=MMt*aryuHH[3]*aryuHH[27]*aryuHH[220];
   aryuHH[156]=aryuHH[220] + 1./2.*aryuHH[156] + aryuHH[247];
   aryuHH[156]=aryuHH[26]*MMt*aryuHH[156];
   aryuHH[104]=3*aryuHH[156] + 1./2.*aryuHH[104] + aryuHH[121];
   aryuHH[104]=aryuHH[26]*aryuHH[104];
   aryuHH[121]=aryuHH[7]*aryuHH[10];
   aryuHH[121]=aryuHH[121] + aryuHH[123];
   aryuHH[156]=aryuHH[31]*aryuHH[121];
   aryuHH[121]=aryuHH[1]*aryuHH[121];
   aryuHH[121]=3*aryuHH[156] + aryuHH[121];
   aryuHH[156]=aryuHH[191] + 1./3.*aryuHH[193];
   aryuHH[191]=aryuHH[12]*aryuHH[156];
   aryuHH[193]=aryuHH[11]*aryuHH[156];
   aryuHH[121]=1./2.*aryuHH[193] + 3./2.*aryuHH[121] + aryuHH[191];
   aryuHH[121]=MMH*aryuHH[121];
   aryuHH[183]=aryuHH[17]*aryuHH[183];
   aryuHH[191]=aryuHH[31]*aryuHH[183];
   aryuHH[183]=aryuHH[1]*aryuHH[183];
   aryuHH[183]=3*aryuHH[191] + aryuHH[183];
   aryuHH[191]=aryuHH[19]*aryuHH[156];
   aryuHH[156]=aryuHH[18]*aryuHH[156];
   aryuHH[121]=1./2.*aryuHH[121] + 1./2.*aryuHH[156] + 1./2.*
   aryuHH[183] + aryuHH[191];
   aryuHH[156]=5./3.*aryuHH[38];
   aryuHH[130]= - 5./3.*aryuHH[36] + 33./8.*aryuHH[130] + aryuHH[156];
   aryuHH[183]=aryuHH[12]*aryuHH[130];
   aryuHH[130]=aryuHH[11]*aryuHH[130];
   aryuHH[191]=aryuHH[7]*aryuHH[148];
   aryuHH[130]=1./2.*aryuHH[130] + 9./4.*aryuHH[191] + aryuHH[183];
   aryuHH[130]=MMH*aryuHH[130];
   aryuHH[140]=aryuHH[36] + aryuHH[140] - aryuHH[38];
   aryuHH[183]=1./2.*aryuHH[140] + aryuHH[310];
   aryuHH[183]=aryuHH[19]*aryuHH[183];
   aryuHH[191]=5./3.*aryuHH[12];
   aryuHH[140]=1./4.*aryuHH[140] + aryuHH[191];
   aryuHH[140]=aryuHH[18]*aryuHH[140];
   aryuHH[148]=aryuHH[17]*aryuHH[148];
   aryuHH[193]= - aryuHH[19] + aryuHH[18];
   aryuHH[193]=aryuHH[11]*aryuHH[193];
   aryuHH[130]=1./2.*aryuHH[130] + 5./6.*aryuHH[193] + aryuHH[140] + 3./
   4.*aryuHH[148] + aryuHH[183];
   aryuHH[130]=aryuHH[68]*aryuHH[130];
   aryuHH[140]= - aryuHH[19] - aryuHH[18];
   aryuHH[140]=aryuHH[18]*aryuHH[140];
   aryuHH[140]=aryuHH[299] + 1./2.*aryuHH[140];
   aryuHH[140]=aryuHH[3]*aryuHH[68]*aryuHH[140];
   aryuHH[104]=aryuHH[104] + 17./4.*aryuHH[140] + 1./2.*aryuHH[121] + 
   aryuHH[130];
   aryuHH[104]=aryuHH[256]*aryuHH[104];
   aryuHH[121]= - aryuHH[14] + aryuHH[163];
   aryuHH[130]= - aryuHH[7]*aryuHH[10];
   aryuHH[140]=9./4.*aryuHH[194] + aryuHH[121] + 9./4.*aryuHH[130];
   aryuHH[140]=aryuHH[31]*aryuHH[140];
   aryuHH[121]=3./4.*aryuHH[194] + 1./3.*aryuHH[121] + 3./4.*
   aryuHH[130];
   aryuHH[121]=aryuHH[1]*aryuHH[121];
   aryuHH[130]=aryuHH[249] + aryuHH[188];
   aryuHH[148]=aryuHH[31]*aryuHH[130];
   aryuHH[130]=aryuHH[1]*aryuHH[130];
   aryuHH[130]=aryuHH[148] + 1./3.*aryuHH[130];
   aryuHH[148]=aryuHH[12]*aryuHH[130];
   aryuHH[130]=aryuHH[11]*aryuHH[130];
   aryuHH[121]=1./2.*aryuHH[130] + aryuHH[148] + aryuHH[140] + 
   aryuHH[121];
   aryuHH[121]=MMH*aryuHH[121];
   aryuHH[130]=aryuHH[17]*aryuHH[266];
   aryuHH[140]=aryuHH[31]*aryuHH[130];
   aryuHH[130]=aryuHH[1]*aryuHH[130];
   aryuHH[130]=3*aryuHH[140] + aryuHH[130];
   aryuHH[140]=aryuHH[9] + 1 - aryuHH[10];
   aryuHH[148]=aryuHH[31]*aryuHH[140];
   aryuHH[140]=aryuHH[1]*aryuHH[140];
   aryuHH[140]=aryuHH[148] + 1./3.*aryuHH[140];
   aryuHH[148]=aryuHH[19]*aryuHH[140];
   aryuHH[140]=aryuHH[18]*aryuHH[140];
   aryuHH[121]=aryuHH[121] + 1./2.*aryuHH[140] + 1./2.*aryuHH[130] + 
   aryuHH[148];
   aryuHH[130]=aryuHH[17] - 199*aryuHH[19];
   aryuHH[130]=aryuHH[19]*aryuHH[130];
   aryuHH[140]=149./6.*aryuHH[18] + 1./6.*aryuHH[17] + 17*aryuHH[19];
   aryuHH[140]=aryuHH[18]*aryuHH[140];
   aryuHH[130]=1./3.*aryuHH[130] + aryuHH[140];
   aryuHH[130]=aryuHH[3]*aryuHH[68]*aryuHH[130];
   aryuHH[102]=aryuHH[104] + aryuHH[108] + 1./8.*aryuHH[130] + 1./2.*
   aryuHH[121] + aryuHH[102];
   aryuHH[102]=aryuHH[256]*aryuHH[102];
   aryuHH[104]= - 745./6. + 27*aryuHH[59];
   aryuHH[108]= - 11*aryuHH[60];
   aryuHH[104]=aryuHH[289] + aryuHH[234] - 155./6.*aryuHH[63] + 
   aryuHH[284] + aryuHH[108] + aryuHH[274] - 149./12.*aryuHH[43] - 11./
   2.*aryuHH[37] + aryuHH[248] + 11*aryuHH[58] + aryuHH[238] + 1./4.*
   aryuHH[104] + aryuHH[295];
   aryuHH[121]= - 15./2.*aryuHH[41];
   aryuHH[130]= - 55*aryuHH[37];
   aryuHH[140]= - 47*aryuHH[40];
   aryuHH[148]=aryuHH[140] + aryuHH[130] - 557./3. + aryuHH[121];
   aryuHH[163]= - 26*aryuHH[39];
   aryuHH[183]=11./3.*aryuHH[9] - 20./9. + 3*aryuHH[10];
   aryuHH[183]=2*aryuHH[31]*aryuHH[183];
   aryuHH[194]=aryuHH[244] - 4./3. + aryuHH[10];
   aryuHH[194]=2*aryuHH[1]*aryuHH[194];
   aryuHH[220]= - 25./3.*aryuHH[27];
   aryuHH[148]=aryuHH[220] + aryuHH[180] + aryuHH[194] + aryuHH[183] + 
   aryuHH[177] + aryuHH[304] + 1./2.*aryuHH[148] + aryuHH[163];
   aryuHH[148]=aryuHH[27]*aryuHH[148];
   aryuHH[167]= - 22*aryuHH[18] - 50*aryuHH[28] + aryuHH[167] - 11*
   aryuHH[33] + 3*aryuHH[187] + 22*aryuHH[67];
   aryuHH[167]=aryuHH[3]*aryuHH[167];
   aryuHH[148]=aryuHH[167] + aryuHH[104] + aryuHH[148];
   aryuHH[148]=aryuHH[3]*aryuHH[148];
   aryuHH[187]= - 2./3.*aryuHH[47] - aryuHH[46];
   aryuHH[168]=3./4.*aryuHH[236] + 11./4.*aryuHH[296] + 3./4.*
   aryuHH[168] + 6*aryuHH[50] + 27./4.*aryuHH[53] + 8*aryuHH[187] + 9./
   2.*aryuHH[51];
   aryuHH[187]=2*aryuHH[47] + aryuHH[46];
   aryuHH[236]=53 + aryuHH[295];
   aryuHH[236]=22*aryuHH[63] + 1./2.*aryuHH[236] + 3*aryuHH[64];
   aryuHH[236]=aryuHH[3]*aryuHH[236];
   aryuHH[187]=aryuHH[236] - 11*aryuHH[53] + 32./3.*aryuHH[187] + 
   aryuHH[278];
   aryuHH[187]=MMt*aryuHH[3]*aryuHH[187];
   aryuHH[148]=aryuHH[187] + aryuHH[168] + aryuHH[148];
   aryuHH[148]=MMt*aryuHH[148];
   aryuHH[236]= - 59./2. - 3*aryuHH[41];
   aryuHH[236]=1./2.*aryuHH[236] - 11*aryuHH[37];
   aryuHH[140]=5*aryuHH[236] + aryuHH[140];
   aryuHH[140]=aryuHH[305] + aryuHH[194] + aryuHH[183] + aryuHH[177] + 
   aryuHH[304] + 1./2.*aryuHH[140] + aryuHH[163];
   aryuHH[140]=aryuHH[28]*aryuHH[140];
   aryuHH[236]=3*aryuHH[34];
   aryuHH[238]= - 2./3.*aryuHH[48] + aryuHH[236];
   aryuHH[224]= - 149./12.*aryuHH[20] + aryuHH[306] - 679./24.*
   aryuHH[45] + 177./8.*aryuHH[33] - 58./3.*aryuHH[67] + 2*aryuHH[238]
    + aryuHH[224];
   aryuHH[238]= - 133./2. - 32*aryuHH[29];
   aryuHH[238]=1./3.*aryuHH[238] + aryuHH[307];
   aryuHH[238]=aryuHH[19]*aryuHH[238];
   aryuHH[244]=12*aryuHH[19] - 55./2.*aryuHH[28];
   aryuHH[244]=aryuHH[27]*aryuHH[244];
   aryuHH[247]= - 1./4.*aryuHH[27] + aryuHH[270] + 8 - 43./3.*
   aryuHH[29];
   aryuHH[247]=aryuHH[18]*aryuHH[247];
   aryuHH[140]=aryuHH[247] + aryuHH[244] + aryuHH[140] + aryuHH[224] + 
   aryuHH[238];
   aryuHH[140]=aryuHH[3]*aryuHH[140];
   aryuHH[238]= - 21*aryuHH[62] + 425./4. + 9*aryuHH[59];
   aryuHH[219]=19./2.*aryuHH[58] + 1./2.*aryuHH[238] + aryuHH[219];
   aryuHH[219]=11./2.*aryuHH[237] + aryuHH[303] - aryuHH[30] - 11./3.*
   aryuHH[29] + aryuHH[280] + 83./12.*aryuHH[63] + aryuHH[298] + 9*
   aryuHH[60] + 5./2.*aryuHH[55] + aryuHH[64] + 1./2.*aryuHH[219] - 7*
   aryuHH[43];
   aryuHH[219]=1./2.*aryuHH[219];
   aryuHH[237]=3*aryuHH[41];
   aryuHH[238]=541./3. + aryuHH[237];
   aryuHH[248]=11*aryuHH[37];
   aryuHH[114]=aryuHH[114] + 1./2.*aryuHH[238] + aryuHH[248];
   aryuHH[114]=aryuHH[109] + aryuHH[154] + 1./2.*aryuHH[114] + 
   aryuHH[145];
   aryuHH[238]=2*aryuHH[27];
   aryuHH[114]=aryuHH[238] + aryuHH[142] + aryuHH[184] + 1./2.*
   aryuHH[114] + aryuHH[179];
   aryuHH[114]=aryuHH[27]*aryuHH[114];
   aryuHH[249]=71 + 53*aryuHH[29];
   aryuHH[249]=aryuHH[312] + 1./9.*aryuHH[249] + aryuHH[207];
   aryuHH[249]=1./2.*aryuHH[11]*aryuHH[249];
   aryuHH[252]=7 + 16*aryuHH[29];
   aryuHH[252]=2./9.*aryuHH[252] + aryuHH[309];
   aryuHH[252]=aryuHH[12]*aryuHH[252];
   aryuHH[254]= - 3./2.*aryuHH[98] - 11*aryuHH[99];
   aryuHH[254]=1./2.*aryuHH[28]*aryuHH[254];
   aryuHH[221]=11./2.*aryuHH[221];
   aryuHH[255]=4./3.*aryuHH[46] - aryuHH[53];
   aryuHH[255]=MMH*aryuHH[255];
   aryuHH[114]=aryuHH[148] + aryuHH[140] + aryuHH[255] + aryuHH[249] + 
   aryuHH[221] + aryuHH[114] + aryuHH[254] + aryuHH[315] + aryuHH[219]
    + aryuHH[252];
   aryuHH[114]=MMt*aryuHH[114];
   aryuHH[140]=65./8. - aryuHH[13];
   aryuHH[148]= - aryuHH[97]*aryuHH[45];
   aryuHH[252]= - 5./8.*aryuHH[63];
   aryuHH[140]=5./2.*aryuHH[148] + aryuHH[172] + aryuHH[252] + 1./3.*
   aryuHH[140] - 5./8.*aryuHH[55];
   aryuHH[148]= - 17*aryuHH[29];
   aryuHH[172]=215./6. + aryuHH[148];
   aryuHH[257]= - 5./18.*aryuHH[9];
   aryuHH[172]=5*aryuHH[206] + aryuHH[257] + 1./18.*aryuHH[172] - 
   aryuHH[30];
   aryuHH[172]=aryuHH[11]*aryuHH[172];
   aryuHH[258]= - 17./4.*aryuHH[29];
   aryuHH[259]=37./3. + aryuHH[258];
   aryuHH[207]= - 5./36.*aryuHH[9] + 1./9.*aryuHH[259] + aryuHH[207];
   aryuHH[207]=aryuHH[12]*aryuHH[207];
   aryuHH[259]=5./3. - 17./8.*aryuHH[29];
   aryuHH[260]=aryuHH[7]*aryuHH[259];
   aryuHH[123]=1./4.*aryuHH[172] + 5./4.*aryuHH[196] + aryuHH[207] + 5./
   8.*aryuHH[123] + 9./4.*aryuHH[308] + 1./2.*aryuHH[140] + aryuHH[260]
   ;
   aryuHH[123]=1./2.*MMH*aryuHH[123];
   aryuHH[140]=10./3. + aryuHH[258];
   aryuHH[140]= - 5./12.*aryuHH[9] + 1./3.*aryuHH[140] + aryuHH[201];
   aryuHH[140]=aryuHH[27]*aryuHH[140];
   aryuHH[172]= - 40./3. + 17*aryuHH[29];
   aryuHH[172]=5./3.*aryuHH[9] + 1./3.*aryuHH[172] + aryuHH[319];
   aryuHH[207]=aryuHH[3]*aryuHH[28]*aryuHH[172];
   aryuHH[172]=MMt*aryuHH[3]*aryuHH[27]*aryuHH[172];
   aryuHH[140]=aryuHH[172] + aryuHH[140] + aryuHH[207];
   aryuHH[140]=MMt*aryuHH[140];
   aryuHH[172]=aryuHH[26]*aryuHH[140];
   aryuHH[207]= - 5./24.*aryuHH[9] + 1./3.*aryuHH[259] + aryuHH[234];
   aryuHH[207]=aryuHH[17]*aryuHH[207];
   aryuHH[234]=499./3. + aryuHH[148];
   aryuHH[234]=aryuHH[257] + 1./18.*aryuHH[234] - aryuHH[30];
   aryuHH[234]=aryuHH[19]*aryuHH[234];
   aryuHH[207]=1./2.*aryuHH[234] + aryuHH[207] - 1./32.*aryuHH[20] + 
   aryuHH[323] + 79./16.*aryuHH[45] - 5./32.*aryuHH[33] - aryuHH[61] + 
   1./6.*aryuHH[67];
   aryuHH[234]=241./12. + aryuHH[148];
   aryuHH[258]=5./2.*aryuHH[206];
   aryuHH[234]=1./4.*aryuHH[27] + aryuHH[258] + aryuHH[257] + 1./18.*
   aryuHH[234] - aryuHH[30];
   aryuHH[234]=1./4.*aryuHH[18]*aryuHH[234];
   aryuHH[257]= - 143./16. + 41./3.*aryuHH[12];
   aryuHH[259]=5./8.*aryuHH[196];
   aryuHH[257]=1./3.*aryuHH[257] + aryuHH[259];
   aryuHH[257]=aryuHH[28]*aryuHH[257];
   aryuHH[260]= - 107*aryuHH[19] - 185*aryuHH[28];
   aryuHH[260]=aryuHH[28]*aryuHH[260];
   aryuHH[266]= - aryuHH[18]*aryuHH[28];
   aryuHH[260]=aryuHH[260] + 71*aryuHH[266];
   aryuHH[260]=aryuHH[3]*aryuHH[260];
   aryuHH[267]= - aryuHH[19] + 13./8.*aryuHH[28];
   aryuHH[267]=3*aryuHH[27]*aryuHH[267];
   aryuHH[268]=71./18.*aryuHH[217];
   aryuHH[114]=aryuHH[172] + aryuHH[114] + 1./12.*aryuHH[260] + 
   aryuHH[123] + aryuHH[268] + aryuHH[234] + aryuHH[267] + aryuHH[207]
    + aryuHH[257];
   aryuHH[114]=aryuHH[26]*aryuHH[114];
   aryuHH[172]=aryuHH[17] + 233*aryuHH[19];
   aryuHH[172]=aryuHH[19]*aryuHH[172];
   aryuHH[172]=1./24.*aryuHH[172] + aryuHH[176];
   aryuHH[172]=aryuHH[3]*aryuHH[68]*aryuHH[172];
   aryuHH[102]=aryuHH[102] + aryuHH[114] + aryuHH[172] + aryuHH[115] + 
   aryuHH[103];
   aryuHH[102]=aryuHH[256]*aryuHH[102];
   aryuHH[103]= - 7./2. + 9*aryuHH[54];
   aryuHH[114]= - 5./4.*aryuHH[59];
   aryuHH[115]=1./2.*aryuHH[35];
   aryuHH[172]=3./2.*aryuHH[56];
   aryuHH[176]=9./4.*aryuHH[43];
   aryuHH[103]=aryuHH[298] + aryuHH[231] + aryuHH[176] - aryuHH[58] + 
   aryuHH[172] + aryuHH[62] + aryuHH[115] + aryuHH[113] + aryuHH[114]
    + 1./2.*aryuHH[103] - 3*aryuHH[65];
   aryuHH[257]= - 5./2.*aryuHH[66] + 11*aryuHH[32];
   aryuHH[236]=3./4.*aryuHH[21] - 33./2.*aryuHH[45] + aryuHH[325] - 7./
   4.*aryuHH[61] - 17./4.*aryuHH[44] + 1./2.*aryuHH[257] + aryuHH[236];
   aryuHH[257]=3*aryuHH[236];
   aryuHH[260]= - 3*aryuHH[17];
   aryuHH[274]=aryuHH[260] - aryuHH[19];
   aryuHH[278]= - 18*aryuHH[28];
   aryuHH[274]=5*aryuHH[274] + aryuHH[278];
   aryuHH[274]=aryuHH[27]*aryuHH[274];
   aryuHH[280]= - 27./4.*aryuHH[17];
   aryuHH[175]= - 9./4. + aryuHH[175];
   aryuHH[175]=aryuHH[19]*aryuHH[175];
   aryuHH[284]= - 5./2.*aryuHH[41];
   aryuHH[286]=aryuHH[284] + 3 - 10*aryuHH[35];
   aryuHH[286]=aryuHH[28]*aryuHH[286];
   aryuHH[289]=47./4. - 18*aryuHH[27];
   aryuHH[289]=aryuHH[18]*aryuHH[289];
   aryuHH[175]=aryuHH[289] + aryuHH[274] + 3*aryuHH[286] + 3*
   aryuHH[175] + aryuHH[280] + aryuHH[257] - 41./4.*aryuHH[20];
   aryuHH[175]=aryuHH[3]*aryuHH[175];
   aryuHH[274]=3./4.*aryuHH[41] + 1 + 21./2.*aryuHH[35];
   aryuHH[127]=aryuHH[180] + aryuHH[127] + aryuHH[137] + aryuHH[274] + 
   aryuHH[264];
   aryuHH[127]=aryuHH[27]*aryuHH[127];
   aryuHH[286]=2*aryuHH[66] - aryuHH[32];
   aryuHH[289]= - aryuHH[27]*aryuHH[17];
   aryuHH[186]=4*aryuHH[289] + aryuHH[278] + aryuHH[195] - 4*aryuHH[17]
    + aryuHH[186] - 4*aryuHH[44] + 4*aryuHH[286] - aryuHH[34];
   aryuHH[186]=aryuHH[3]*aryuHH[186];
   aryuHH[195]= - 5*aryuHH[35];
   aryuHH[278]= - 2 + aryuHH[195];
   aryuHH[269]=aryuHH[269] + 2*aryuHH[278] + aryuHH[284];
   aryuHH[269]=aryuHH[27]*aryuHH[269];
   aryuHH[186]=aryuHH[186] + aryuHH[269] - aryuHH[63] - 5./4.*
   aryuHH[64] - 13./2.*aryuHH[43] - 1./2.*aryuHH[41] - 2*aryuHH[35] + 5.
   /2.*aryuHH[59] - 5*aryuHH[65] + 3 + 4*aryuHH[57];
   aryuHH[186]=aryuHH[3]*aryuHH[186];
   aryuHH[269]=aryuHH[64] + 9 + 8*aryuHH[65];
   aryuHH[269]=aryuHH[3]*aryuHH[269];
   aryuHH[269]= - 4*aryuHH[52] + aryuHH[269];
   aryuHH[269]=MMt*aryuHH[3]*aryuHH[269];
   aryuHH[284]= - 6*aryuHH[49] + aryuHH[52];
   aryuHH[293]=1 + aryuHH[64];
   aryuHH[294]=aryuHH[98]*aryuHH[293];
   aryuHH[186]=aryuHH[269] + aryuHH[186] + 1./4.*aryuHH[294] + 2*
   aryuHH[284] - aryuHH[51];
   aryuHH[186]=MMt*aryuHH[186];
   aryuHH[269]=3*aryuHH[314];
   aryuHH[296]= - 3./2.*aryuHH[28]*aryuHH[98];
   aryuHH[299]= - 27./2.*aryuHH[7];
   aryuHH[301]= - aryuHH[50] + 9*aryuHH[49] + 1./2.*aryuHH[52];
   aryuHH[301]=MMH*aryuHH[301];
   aryuHH[306]=aryuHH[12]*aryuHH[30];
   aryuHH[103]=3*aryuHH[186] + aryuHH[175] + 3./2.*aryuHH[301] + 
   aryuHH[127] + aryuHH[296] + aryuHH[269] + aryuHH[306] + aryuHH[303]
    - aryuHH[30] + aryuHH[299] + 3*aryuHH[103] + 13./4.*aryuHH[63];
   aryuHH[103]=MMt*aryuHH[103];
   aryuHH[127]= - 3*aryuHH[57];
   aryuHH[175]=29./4. + aryuHH[127];
   aryuHH[186]= - 3*aryuHH[54];
   aryuHH[307]=3*aryuHH[42];
   aryuHH[310]= - 1./2.*aryuHH[56];
   aryuHH[311]=3./4.*aryuHH[43];
   aryuHH[175]=aryuHH[311] + aryuHH[310] + aryuHH[307] + 1./4.*
   aryuHH[175] + aryuHH[186];
   aryuHH[314]=1 - 3./2.*aryuHH[35];
   aryuHH[317]=aryuHH[36] + aryuHH[314] + aryuHH[38];
   aryuHH[319]=9*aryuHH[7];
   aryuHH[320]=3./2.*aryuHH[12];
   aryuHH[317]=aryuHH[320] + 1./2.*aryuHH[317] + aryuHH[319];
   aryuHH[317]=aryuHH[27]*aryuHH[317];
   aryuHH[321]= - 7./2.*aryuHH[29];
   aryuHH[323]= - 71 + aryuHH[321];
   aryuHH[323]=1./9.*aryuHH[323] + aryuHH[242];
   aryuHH[323]=1./2.*aryuHH[323] + aryuHH[233];
   aryuHH[323]=aryuHH[11]*aryuHH[323];
   aryuHH[324]=1 - 1./8.*aryuHH[29];
   aryuHH[325]=aryuHH[7]*aryuHH[324];
   aryuHH[324]=7./9.*aryuHH[324] + 1./8.*aryuHH[30];
   aryuHH[324]=aryuHH[12]*aryuHH[324];
   aryuHH[175]=1./4.*aryuHH[323] + aryuHH[317] + aryuHH[324] + 9./16.*
   aryuHH[119] + 7./2.*aryuHH[325] + 11./8.*aryuHH[15] + aryuHH[252] + 
   aryuHH[230] + 3*aryuHH[175] - 11./8.*aryuHH[55];
   aryuHH[175]=MMH*aryuHH[175];
   aryuHH[252]=5./2.*aryuHH[12] + aryuHH[299] + aryuHH[137] + 
   aryuHH[264] + 8 + 15./2.*aryuHH[35];
   aryuHH[252]=aryuHH[28]*aryuHH[252];
   aryuHH[317]=1 + aryuHH[321];
   aryuHH[317]=1./3.*aryuHH[317] + aryuHH[212];
   aryuHH[317]=aryuHH[17]*aryuHH[317];
   aryuHH[323]=9*aryuHH[17];
   aryuHH[324]=9*aryuHH[28];
   aryuHH[325]=aryuHH[324] + aryuHH[323] + 5./2.*aryuHH[19];
   aryuHH[325]=aryuHH[27]*aryuHH[325];
   aryuHH[326]= - 9*aryuHH[17];
   aryuHH[327]=aryuHH[326] - 53*aryuHH[19];
   aryuHH[328]= - 3*aryuHH[28];
   aryuHH[327]=1./4.*aryuHH[327] + aryuHH[328];
   aryuHH[327]=aryuHH[28]*aryuHH[327];
   aryuHH[327]=aryuHH[327] + 31./4.*aryuHH[266];
   aryuHH[327]=aryuHH[3]*aryuHH[327];
   aryuHH[150]=aryuHH[150] - 3*aryuHH[32] - aryuHH[34];
   aryuHH[150]=1./2.*aryuHH[150] - aryuHH[61];
   aryuHH[329]=3*aryuHH[150];
   aryuHH[330]= - 3./4.*aryuHH[21];
   aryuHH[331]=35 - aryuHH[29];
   aryuHH[331]=7./9.*aryuHH[331] + aryuHH[30];
   aryuHH[331]=aryuHH[19]*aryuHH[331];
   aryuHH[332]= - 73 - aryuHH[29];
   aryuHH[332]=7./9.*aryuHH[332] + aryuHH[30];
   aryuHH[332]=1./4.*aryuHH[332] + 9*aryuHH[27];
   aryuHH[332]=aryuHH[18]*aryuHH[332];
   aryuHH[103]=aryuHH[103] + aryuHH[327] + aryuHH[175] + 1./2.*
   aryuHH[332] + 1./2.*aryuHH[325] + aryuHH[252] + 1./4.*aryuHH[331] + 
   1./4.*aryuHH[317] + 27./8.*aryuHH[20] + aryuHH[330] + 7*aryuHH[45]
    - 11./8.*aryuHH[33] + aryuHH[329] + 5./2.*aryuHH[67];
   aryuHH[103]=MMt*aryuHH[103];
   aryuHH[175]= - 1./3.*aryuHH[19];
   aryuHH[128]=aryuHH[128] + aryuHH[175];
   aryuHH[252]=13./3.*aryuHH[128] + 5./8.*aryuHH[28];
   aryuHH[252]=aryuHH[28]*aryuHH[252];
   aryuHH[317]=5./2. - 13*aryuHH[7];
   aryuHH[317]=1./2.*aryuHH[317] - 13./9.*aryuHH[12];
   aryuHH[317]=aryuHH[28]*aryuHH[317];
   aryuHH[317]=71./36.*aryuHH[251] - 5./4.*aryuHH[45] + aryuHH[317];
   aryuHH[317]=MMH*aryuHH[317];
   aryuHH[103]=aryuHH[103] + 1./2.*aryuHH[317] + aryuHH[252] + 97./72.*
   aryuHH[266];
   aryuHH[252]=1./2.*aryuHH[314];
   aryuHH[314]=aryuHH[320] + aryuHH[319] + aryuHH[252] + aryuHH[38];
   aryuHH[314]=aryuHH[27]*aryuHH[314];
   aryuHH[317]=19./3. + aryuHH[30];
   aryuHH[317]=1./2.*aryuHH[317] + aryuHH[233];
   aryuHH[317]=aryuHH[11]*aryuHH[317];
   aryuHH[325]=127./4. - 9*aryuHH[57];
   aryuHH[327]=1 + aryuHH[309];
   aryuHH[327]=aryuHH[12]*aryuHH[327];
   aryuHH[229]=1./4.*aryuHH[317] + aryuHH[314] + aryuHH[327] + 9./8.*
   aryuHH[119] + aryuHH[214] - 7./24.*aryuHH[15] + 25./24.*aryuHH[63]
    + aryuHH[230] + 7./24.*aryuHH[55] + aryuHH[176] + aryuHH[229] + 9*
   aryuHH[42] + 1./4.*aryuHH[325] - 9*aryuHH[54];
   aryuHH[229]=MMH*aryuHH[229];
   aryuHH[314]= - 61./2. + 27*aryuHH[54];
   aryuHH[317]=1./2.*aryuHH[306];
   aryuHH[161]=aryuHH[296] + aryuHH[269] + aryuHH[317] + aryuHH[303] - 
   aryuHH[30] + aryuHH[299] - 41./12.*aryuHH[63] - 9./2.*aryuHH[16] + 3.
   /2.*aryuHH[64] + 27./4.*aryuHH[43] + aryuHH[161] + 9./2.*aryuHH[56]
    + aryuHH[295] + 3./2.*aryuHH[35] - 27./2.*aryuHH[42] - 15./4.*
   aryuHH[59] + 1./2.*aryuHH[314] - 9*aryuHH[65];
   aryuHH[295]= - 3./2. - aryuHH[30];
   aryuHH[295]=aryuHH[19]*aryuHH[295];
   aryuHH[257]=9./2.*aryuHH[295] + aryuHH[280] + aryuHH[257] + 37./12.*
   aryuHH[20];
   aryuHH[280]=7*aryuHH[19];
   aryuHH[295]=aryuHH[260] + aryuHH[280];
   aryuHH[314]= - 9*aryuHH[28];
   aryuHH[295]=5./2.*aryuHH[295] + aryuHH[314];
   aryuHH[295]=aryuHH[27]*aryuHH[295];
   aryuHH[325]= - 19./6. + aryuHH[192];
   aryuHH[327]=17./3.*aryuHH[27];
   aryuHH[325]=1./4.*aryuHH[325] + aryuHH[327];
   aryuHH[325]=aryuHH[18]*aryuHH[325];
   aryuHH[331]= - 5./4.*aryuHH[41];
   aryuHH[195]=aryuHH[331] + 3./2. + aryuHH[195];
   aryuHH[195]=aryuHH[28]*aryuHH[195];
   aryuHH[257]=aryuHH[325] + aryuHH[295] + 1./2.*aryuHH[257] + 3*
   aryuHH[195];
   aryuHH[257]=aryuHH[3]*aryuHH[257];
   aryuHH[274]=1./2.*aryuHH[274];
   aryuHH[295]=aryuHH[142] + aryuHH[299] + aryuHH[274] + aryuHH[264];
   aryuHH[295]=aryuHH[27]*aryuHH[295];
   aryuHH[286]=2*aryuHH[289] + aryuHH[314] - aryuHH[19] - 2*aryuHH[17]
    + aryuHH[61] - 2*aryuHH[44] + 2*aryuHH[286] + aryuHH[300];
   aryuHH[286]=aryuHH[3]*aryuHH[286];
   aryuHH[278]= - aryuHH[27] + aryuHH[278] + aryuHH[331];
   aryuHH[278]=aryuHH[27]*aryuHH[278];
   aryuHH[226]=aryuHH[286] + aryuHH[278] + aryuHH[226] - 5./8.*
   aryuHH[64] - 13./4.*aryuHH[43] - 1./4.*aryuHH[41] - aryuHH[35] + 5./
   4.*aryuHH[59] - 5./2.*aryuHH[65] + 3./2. + 2*aryuHH[57];
   aryuHH[226]=aryuHH[3]*aryuHH[226];
   aryuHH[278]=aryuHH[231] + 9./2. + 4*aryuHH[65];
   aryuHH[278]=aryuHH[3]*aryuHH[278];
   aryuHH[278]= - 2*aryuHH[52] + aryuHH[278];
   aryuHH[278]=MMt*aryuHH[3]*aryuHH[278];
   aryuHH[226]=aryuHH[278] + aryuHH[226] + 1./8.*aryuHH[294] + 
   aryuHH[284] - 1./2.*aryuHH[51];
   aryuHH[226]=3*MMt*aryuHH[226];
   aryuHH[278]=3./4.*aryuHH[301];
   aryuHH[161]=aryuHH[226] + aryuHH[257] + aryuHH[278] + 1./4.*
   aryuHH[290] + 1./2.*aryuHH[161] + aryuHH[295];
   aryuHH[161]=MMt*aryuHH[161];
   aryuHH[257]=1 + aryuHH[30];
   aryuHH[257]=aryuHH[17]*aryuHH[257];
   aryuHH[284]=29./2. + aryuHH[30];
   aryuHH[284]=aryuHH[19]*aryuHH[284];
   aryuHH[257]=1./2.*aryuHH[284] + 3./4.*aryuHH[257] - 13./8.*
   aryuHH[20] + aryuHH[330] + 31./3.*aryuHH[45] + 7./24.*aryuHH[33] + 
   aryuHH[329] - 25./6.*aryuHH[67];
   aryuHH[284]= - 27./4.*aryuHH[7];
   aryuHH[264]=aryuHH[12] + aryuHH[284] + aryuHH[264] + 37./3. + 15./4.
   *aryuHH[35];
   aryuHH[264]=aryuHH[28]*aryuHH[264];
   aryuHH[286]=aryuHH[324] + aryuHH[323] - 35./2.*aryuHH[19];
   aryuHH[286]=aryuHH[27]*aryuHH[286];
   aryuHH[294]=65./2. + aryuHH[30];
   aryuHH[295]= - 17./3.*aryuHH[27];
   aryuHH[294]=1./2.*aryuHH[294] + aryuHH[295];
   aryuHH[294]=aryuHH[18]*aryuHH[294];
   aryuHH[300]=aryuHH[326] + 113*aryuHH[19];
   aryuHH[300]=1./4.*aryuHH[300] + aryuHH[328];
   aryuHH[300]=aryuHH[28]*aryuHH[300];
   aryuHH[300]=aryuHH[300] + 39./4.*aryuHH[322];
   aryuHH[300]=aryuHH[3]*aryuHH[300];
   aryuHH[161]=aryuHH[161] + 1./2.*aryuHH[300] + 1./2.*aryuHH[229] + 1./
   4.*aryuHH[251] + 1./4.*aryuHH[294] + 1./4.*aryuHH[286] + 1./2.*
   aryuHH[257] + aryuHH[264];
   aryuHH[161]=MMt*aryuHH[161];
   aryuHH[229]= - aryuHH[18]*aryuHH[27];
   aryuHH[257]=aryuHH[27]*aryuHH[19];
   aryuHH[264]=aryuHH[257] + aryuHH[229];
   aryuHH[264]=1./2.*aryuHH[264] + aryuHH[251];
   aryuHH[286]= - aryuHH[28]*aryuHH[19];
   aryuHH[294]=aryuHH[286] + 5./2.*aryuHH[322];
   aryuHH[294]=aryuHH[3]*aryuHH[294];
   aryuHH[300]=aryuHH[212] + aryuHH[27];
   aryuHH[300]=aryuHH[18]*aryuHH[300];
   aryuHH[301]= - aryuHH[27]*aryuHH[19];
   aryuHH[300]=aryuHH[301] + aryuHH[300];
   aryuHH[300]=aryuHH[3]*aryuHH[300];
   aryuHH[300]=1./2.*aryuHH[290] + aryuHH[300];
   aryuHH[300]=MMt*aryuHH[300];
   aryuHH[264]=aryuHH[300] + 1./2.*aryuHH[264] + aryuHH[294];
   aryuHH[264]=aryuHH[171]*MMt*aryuHH[264];
   aryuHH[294]= - 25./6. - 9*aryuHH[7];
   aryuHH[294]=1./2.*aryuHH[294] - aryuHH[12];
   aryuHH[294]=aryuHH[28]*aryuHH[294];
   aryuHH[300]=19./12.*aryuHH[217];
   aryuHH[294]=aryuHH[300] + 25./12.*aryuHH[45] + aryuHH[294];
   aryuHH[294]=MMH*aryuHH[294];
   aryuHH[314]= - 25./24.*aryuHH[28] + aryuHH[292] - aryuHH[19];
   aryuHH[314]=aryuHH[28]*aryuHH[314];
   aryuHH[294]=1./2.*aryuHH[294] + aryuHH[314] + 13./24.*aryuHH[322];
   aryuHH[161]=1./2.*aryuHH[264] + 1./2.*aryuHH[294] + aryuHH[161];
   aryuHH[161]=aryuHH[171]*aryuHH[161];
   aryuHH[161]=aryuHH[103] + aryuHH[161];
   aryuHH[161]=aryuHH[171]*aryuHH[161];
   aryuHH[264]= - 13 - 7./4.*aryuHH[29];
   aryuHH[264]=1./3.*aryuHH[264] + aryuHH[270];
   aryuHH[264]=aryuHH[27]*aryuHH[264];
   aryuHH[270]=52 + 7*aryuHH[29];
   aryuHH[197]=1./3.*aryuHH[270] + aryuHH[197];
   aryuHH[270]=aryuHH[28]*aryuHH[197];
   aryuHH[270]=aryuHH[270] + 52./3.*aryuHH[218];
   aryuHH[270]=aryuHH[3]*aryuHH[270];
   aryuHH[197]=MMt*aryuHH[3]*aryuHH[27]*aryuHH[197];
   aryuHH[197]=aryuHH[197] + aryuHH[264] + aryuHH[270];
   aryuHH[197]=MMt*aryuHH[197];
   aryuHH[264]=pow(aryuHH[28],2);
   aryuHH[270]=aryuHH[3]*aryuHH[264];
   aryuHH[294]=aryuHH[265] + 4*aryuHH[270];
   aryuHH[197]=13./3.*aryuHH[294] + aryuHH[197];
   aryuHH[197]=MMt*aryuHH[197];
   aryuHH[294]=2 - aryuHH[30];
   aryuHH[314]=aryuHH[28]*aryuHH[294];
   aryuHH[314]=aryuHH[314] + 2*aryuHH[218];
   aryuHH[314]=aryuHH[3]*aryuHH[314];
   aryuHH[313]=aryuHH[27]*aryuHH[313];
   aryuHH[294]=MMt*aryuHH[3]*aryuHH[27]*aryuHH[294];
   aryuHH[294]=aryuHH[294] + 1./2.*aryuHH[313] + aryuHH[314];
   aryuHH[294]=MMt*aryuHH[294];
   aryuHH[294]=aryuHH[294] + 1./2.*aryuHH[265] + 2*aryuHH[270];
   aryuHH[294]=aryuHH[171]*MMt*aryuHH[294];
   aryuHH[294]=aryuHH[197] + 3*aryuHH[294];
   aryuHH[294]=aryuHH[26]*aryuHH[171]*aryuHH[294];
   aryuHH[161]=aryuHH[161] + aryuHH[294];
   aryuHH[161]=aryuHH[26]*aryuHH[161];
   aryuHH[197]=aryuHH[26]*aryuHH[197];
   aryuHH[103]=aryuHH[103] + aryuHH[197];
   aryuHH[103]=aryuHH[26]*aryuHH[103];
   aryuHH[197]= - 5./2. + 3*aryuHH[54];
   aryuHH[197]=1./2.*aryuHH[197] - aryuHH[65];
   aryuHH[113]=aryuHH[125] - 1./4.*aryuHH[63] + aryuHH[298] + 
   aryuHH[231] + aryuHH[176] - aryuHH[58] + aryuHH[172] + aryuHH[62] + 
   aryuHH[115] + aryuHH[113] + 3*aryuHH[197] + aryuHH[114];
   aryuHH[113]=aryuHH[296] + aryuHH[269] + 3./2.*aryuHH[306] + 
   aryuHH[303] + 3*aryuHH[113] - aryuHH[30];
   aryuHH[114]=aryuHH[142] + aryuHH[299] + aryuHH[274] + aryuHH[137];
   aryuHH[114]=aryuHH[27]*aryuHH[114];
   aryuHH[115]= - 9./2. - 5*aryuHH[30];
   aryuHH[115]=aryuHH[19]*aryuHH[115];
   aryuHH[115]=1./2.*aryuHH[115] - 9./4.*aryuHH[17] + aryuHH[236] - 3./
   4.*aryuHH[20];
   aryuHH[125]= - 5*aryuHH[17] + aryuHH[19];
   aryuHH[125]=1./2.*aryuHH[125] + aryuHH[328];
   aryuHH[125]=aryuHH[27]*aryuHH[125];
   aryuHH[115]=aryuHH[125] + 1./2.*aryuHH[115] + aryuHH[195];
   aryuHH[125]=5./2. - aryuHH[30];
   aryuHH[125]=3./4.*aryuHH[125] - 13*aryuHH[27];
   aryuHH[125]=aryuHH[18]*aryuHH[125];
   aryuHH[115]=3*aryuHH[115] + aryuHH[125];
   aryuHH[115]=aryuHH[3]*aryuHH[115];
   aryuHH[125]=aryuHH[11]*aryuHH[30];
   aryuHH[172]=1./4.*aryuHH[125];
   aryuHH[113]=aryuHH[226] + aryuHH[115] + aryuHH[278] + aryuHH[172] + 
   1./2.*aryuHH[113] + aryuHH[114];
   aryuHH[113]=MMt*aryuHH[113];
   aryuHH[114]=37./4. + aryuHH[127];
   aryuHH[115]=7 + 19./4.*aryuHH[29];
   aryuHH[127]=aryuHH[7]*aryuHH[115];
   aryuHH[114]=1./2.*aryuHH[127] + 1./8.*aryuHH[15] + aryuHH[297] + 1./
   2.*aryuHH[16] - 1./8.*aryuHH[55] + aryuHH[311] + aryuHH[310] + 
   aryuHH[307] + 1./4.*aryuHH[114] + aryuHH[186];
   aryuHH[127]=aryuHH[320] + aryuHH[319] + aryuHH[252] + aryuHH[36];
   aryuHH[127]=aryuHH[27]*aryuHH[127];
   aryuHH[176]=1 + aryuHH[29];
   aryuHH[186]=19./6.*aryuHH[176] + aryuHH[233];
   aryuHH[186]=aryuHH[11]*aryuHH[186];
   aryuHH[115]=aryuHH[12]*aryuHH[115];
   aryuHH[114]=1./4.*aryuHH[186] + aryuHH[127] + 3*aryuHH[114] + 1./3.*
   aryuHH[115];
   aryuHH[114]=MMH*aryuHH[114];
   aryuHH[115]=3 + 5./4.*aryuHH[35];
   aryuHH[115]=aryuHH[320] + aryuHH[284] + 3*aryuHH[115] + aryuHH[137];
   aryuHH[115]=aryuHH[28]*aryuHH[115];
   aryuHH[127]=1./8.*aryuHH[20] - 1./4.*aryuHH[21] + 3*aryuHH[45] - 1./
   8.*aryuHH[33] + aryuHH[150] + aryuHH[208];
   aryuHH[137]=aryuHH[17]*aryuHH[176];
   aryuHH[150]=19*aryuHH[29];
   aryuHH[186]=119./2. + aryuHH[150];
   aryuHH[186]=aryuHH[19]*aryuHH[186];
   aryuHH[127]=1./6.*aryuHH[186] + 3*aryuHH[127] + 19./4.*aryuHH[137];
   aryuHH[137]=aryuHH[209] + aryuHH[166] - 1./2.*aryuHH[19];
   aryuHH[137]=aryuHH[27]*aryuHH[137];
   aryuHH[166]= - aryuHH[17] - 3*aryuHH[19];
   aryuHH[166]=3./4.*aryuHH[166] - aryuHH[28];
   aryuHH[166]=aryuHH[28]*aryuHH[166];
   aryuHH[166]=3*aryuHH[166] + 101./4.*aryuHH[266];
   aryuHH[166]=aryuHH[3]*aryuHH[166];
   aryuHH[150]=83./2. + aryuHH[150];
   aryuHH[150]=1./6.*aryuHH[150] + 13*aryuHH[27];
   aryuHH[150]=aryuHH[18]*aryuHH[150];
   aryuHH[186]=1./4.*aryuHH[217];
   aryuHH[113]=aryuHH[113] + 1./2.*aryuHH[166] + 1./2.*aryuHH[114] + 
   aryuHH[186] + 1./4.*aryuHH[150] + 3./4.*aryuHH[137] + 1./2.*
   aryuHH[127] + aryuHH[115];
   aryuHH[113]=MMt*aryuHH[113];
   aryuHH[114]=aryuHH[173] + aryuHH[291];
   aryuHH[114]=5*aryuHH[114] - 3./8.*aryuHH[28];
   aryuHH[114]=aryuHH[28]*aryuHH[114];
   aryuHH[115]= - 1./2. + 5*aryuHH[7];
   aryuHH[115]=3./2.*aryuHH[115] + aryuHH[191];
   aryuHH[115]=aryuHH[28]*aryuHH[115];
   aryuHH[115]=aryuHH[300] + 3./4.*aryuHH[45] + aryuHH[115];
   aryuHH[115]=MMH*aryuHH[115];
   aryuHH[114]=1./2.*aryuHH[115] + aryuHH[114] + 29./24.*aryuHH[322];
   aryuHH[115]= - 10 - 19*aryuHH[29];
   aryuHH[127]=aryuHH[28]*aryuHH[115];
   aryuHH[127]=aryuHH[127] + 10*aryuHH[265];
   aryuHH[127]=aryuHH[3]*aryuHH[127];
   aryuHH[137]=5 + 19./2.*aryuHH[29];
   aryuHH[137]=aryuHH[27]*aryuHH[137];
   aryuHH[115]=MMt*aryuHH[3]*aryuHH[27]*aryuHH[115];
   aryuHH[115]=aryuHH[115] + 1./2.*aryuHH[137] + aryuHH[127];
   aryuHH[115]=MMt*aryuHH[115];
   aryuHH[127]= - aryuHH[3]*aryuHH[264];
   aryuHH[137]=1./2.*aryuHH[218] + 2*aryuHH[127];
   aryuHH[115]=5*aryuHH[137] + aryuHH[115];
   aryuHH[115]=aryuHH[26]*MMt*aryuHH[115];
   aryuHH[113]=aryuHH[115] + 1./2.*aryuHH[114] + aryuHH[113];
   aryuHH[113]=aryuHH[26]*aryuHH[113];
   aryuHH[114]=aryuHH[201] + 17*aryuHH[27];
   aryuHH[114]=aryuHH[18]*aryuHH[114];
   aryuHH[115]= - aryuHH[19]*aryuHH[30];
   aryuHH[114]=aryuHH[114] + 3*aryuHH[115] + 17*aryuHH[301];
   aryuHH[114]=aryuHH[3]*aryuHH[114];
   aryuHH[137]=aryuHH[38] - aryuHH[36];
   aryuHH[150]=aryuHH[27]*aryuHH[137];
   aryuHH[114]=1./2.*aryuHH[114] + aryuHH[172] + aryuHH[317] + 2*
   aryuHH[150];
   aryuHH[114]=MMt*aryuHH[114];
   aryuHH[166]=aryuHH[216] + aryuHH[308];
   aryuHH[166]=9./2.*aryuHH[166] + aryuHH[262];
   aryuHH[172]= - aryuHH[38] + aryuHH[36];
   aryuHH[191]=aryuHH[27]*aryuHH[172];
   aryuHH[166]=1./8.*aryuHH[263] + 1./4.*aryuHH[166] + aryuHH[191];
   aryuHH[166]=MMH*aryuHH[166];
   aryuHH[195]=aryuHH[17]*aryuHH[261];
   aryuHH[181]=3./2.*aryuHH[195] + aryuHH[181];
   aryuHH[195]=2*aryuHH[137] + aryuHH[178];
   aryuHH[195]=aryuHH[28]*aryuHH[195];
   aryuHH[197]=10*aryuHH[286] + 31./4.*aryuHH[322];
   aryuHH[197]=aryuHH[3]*aryuHH[197];
   aryuHH[201]=aryuHH[261] - 17*aryuHH[27];
   aryuHH[201]=aryuHH[18]*aryuHH[201];
   aryuHH[114]=aryuHH[114] + aryuHH[197] + 1./2.*aryuHH[166] + 
   aryuHH[186] + 1./8.*aryuHH[201] + 17./8.*aryuHH[257] + 1./4.*
   aryuHH[181] + aryuHH[195];
   aryuHH[114]=MMt*aryuHH[114];
   aryuHH[166]=aryuHH[27]*aryuHH[261];
   aryuHH[181]=aryuHH[3]*aryuHH[28]*aryuHH[185];
   aryuHH[185]=MMt*aryuHH[3]*aryuHH[27]*aryuHH[185];
   aryuHH[166]=aryuHH[185] + 1./4.*aryuHH[166] + aryuHH[181];
   aryuHH[181]=pow(MMt,2);
   aryuHH[166]=aryuHH[26]*aryuHH[181]*aryuHH[166];
   aryuHH[114]=aryuHH[114] + 3*aryuHH[166];
   aryuHH[114]=aryuHH[26]*aryuHH[114];
   aryuHH[166]=aryuHH[17]*aryuHH[172];
   aryuHH[185]=1./3.*aryuHH[137];
   aryuHH[186]=7./4.*aryuHH[12] + aryuHH[185] + 51./8.*aryuHH[7];
   aryuHH[186]=aryuHH[19]*aryuHH[186];
   aryuHH[195]= - 7./2.*aryuHH[12] + aryuHH[185] - 51./4.*aryuHH[7];
   aryuHH[195]=aryuHH[18]*aryuHH[195];
   aryuHH[166]=7./8.*aryuHH[174] + 1./2.*aryuHH[195] + aryuHH[166] + 
   aryuHH[186];
   aryuHH[174]=aryuHH[7]*aryuHH[172];
   aryuHH[186]=aryuHH[12]*aryuHH[172];
   aryuHH[195]=aryuHH[11]*aryuHH[172];
   aryuHH[174]=1./6.*aryuHH[195] + 3./4.*aryuHH[174] + 1./3.*
   aryuHH[186];
   aryuHH[174]=MMH*aryuHH[174];
   aryuHH[166]=1./2.*aryuHH[166] + aryuHH[174];
   aryuHH[166]=MMH*aryuHH[166];
   aryuHH[174]=11./3.*aryuHH[19];
   aryuHH[186]=17./2.*aryuHH[17] + aryuHH[174];
   aryuHH[186]=aryuHH[19]*aryuHH[186];
   aryuHH[195]= - 11./3.*aryuHH[18] - 17*aryuHH[17] - 11./3.*aryuHH[19]
   ;
   aryuHH[195]=aryuHH[18]*aryuHH[195];
   aryuHH[186]=aryuHH[186] + 1./2.*aryuHH[195];
   aryuHH[166]=1./4.*aryuHH[186] + aryuHH[166];
   aryuHH[166]=aryuHH[68]*aryuHH[166];
   aryuHH[114]=1./2.*aryuHH[166] + aryuHH[114];
   aryuHH[114]=aryuHH[256]*aryuHH[114];
   aryuHH[166]= - 101./8. - 9*aryuHH[82];
   aryuHH[129]=1./2.*aryuHH[166] + aryuHH[129];
   aryuHH[166]=1./6.*aryuHH[38];
   aryuHH[129]=1./12.*aryuHH[36] + aryuHH[166] + aryuHH[272] + 
   aryuHH[230] + 27./8.*aryuHH[42] - 3./4.*aryuHH[85] - 1./2.*
   aryuHH[89] - 3./2.*aryuHH[87] - 9./16.*aryuHH[83] + 3./8.*
   aryuHH[129] - aryuHH[92];
   aryuHH[186]= - 1 - aryuHH[38];
   aryuHH[195]=aryuHH[186] + aryuHH[135];
   aryuHH[197]= - 1./8.*aryuHH[11];
   aryuHH[195]=aryuHH[197] + aryuHH[178] + 1./3.*aryuHH[195] + 
   aryuHH[214];
   aryuHH[195]=aryuHH[11]*aryuHH[195];
   aryuHH[201]=aryuHH[214] + aryuHH[318] - 1./3. + aryuHH[139];
   aryuHH[201]=aryuHH[12]*aryuHH[201];
   aryuHH[208]=63./8.*aryuHH[7];
   aryuHH[209]=aryuHH[208] + 5 + aryuHH[36];
   aryuHH[209]=aryuHH[7]*aryuHH[209];
   aryuHH[216]=aryuHH[243] + 27./4.*aryuHH[73] + aryuHH[75];
   aryuHH[216]=1./4.*aryuHH[76] + 3*aryuHH[216] + 1./2.*aryuHH[78];
   aryuHH[216]=1./4.*MMH*aryuHH[216];
   aryuHH[195]=aryuHH[216] + 1./4.*aryuHH[195] + 1./2.*aryuHH[201] + 
   aryuHH[129] + 3./4.*aryuHH[209];
   aryuHH[195]=MMH*aryuHH[195];
   aryuHH[201]= - 71./16. + aryuHH[36];
   aryuHH[201]=1./2.*aryuHH[201] + aryuHH[319];
   aryuHH[201]=aryuHH[17]*aryuHH[201];
   aryuHH[209]= - 33./2.*aryuHH[44] - 9./2.*aryuHH[70] + 11./2.*
   aryuHH[93] - 9*aryuHH[71] + 9./2.*aryuHH[69] + 11*aryuHH[94];
   aryuHH[201]=1./4.*aryuHH[146] + aryuHH[201] - 9./2.*aryuHH[20] + 1./
   8.*aryuHH[209] - 9*aryuHH[21];
   aryuHH[226]= - 27./16.*aryuHH[12] - 9./32.*aryuHH[7] + 1./6.*
   aryuHH[36] + 4 + aryuHH[213];
   aryuHH[226]=aryuHH[19]*aryuHH[226];
   aryuHH[230]= - 73./48.*aryuHH[12] + 69./32.*aryuHH[7] + 5./24.*
   aryuHH[36] + 2 + 1./8.*aryuHH[38];
   aryuHH[230]=aryuHH[18]*aryuHH[230];
   aryuHH[231]=37./6.*aryuHH[18] + aryuHH[17] - 191./6.*aryuHH[19];
   aryuHH[231]=aryuHH[11]*aryuHH[231];
   aryuHH[195]=1./2.*aryuHH[195] + 1./16.*aryuHH[231] + aryuHH[230] + 1.
   /2.*aryuHH[201] + aryuHH[226];
   aryuHH[195]=MMH*aryuHH[195];
   aryuHH[201]=11./4.*aryuHH[17] - aryuHH[19];
   aryuHH[201]=aryuHH[19]*aryuHH[201];
   aryuHH[226]=23./6.*aryuHH[18] + 15*aryuHH[17] + 49./6.*aryuHH[19];
   aryuHH[226]=aryuHH[18]*aryuHH[226];
   aryuHH[230]=pow(aryuHH[17],2);
   aryuHH[201]=1./2.*aryuHH[226] + 165./16.*aryuHH[230] + aryuHH[201];
   aryuHH[195]=1./4.*aryuHH[201] + aryuHH[195];
   aryuHH[195]=aryuHH[68]*aryuHH[195];
   aryuHH[113]=aryuHH[114] + aryuHH[195] + aryuHH[113];
   aryuHH[113]=aryuHH[256]*aryuHH[113];
   aryuHH[114]= - 1 + aryuHH[213];
   aryuHH[182]=aryuHH[114] + aryuHH[182];
   aryuHH[178]=aryuHH[197] + aryuHH[178] + 1./3.*aryuHH[182] + 
   aryuHH[214];
   aryuHH[178]=aryuHH[11]*aryuHH[178];
   aryuHH[182]=aryuHH[208] + aryuHH[135] + 5 + aryuHH[213];
   aryuHH[182]=aryuHH[7]*aryuHH[182];
   aryuHH[135]=aryuHH[135] - 1 + aryuHH[139];
   aryuHH[135]=1./3.*aryuHH[135] + aryuHH[214];
   aryuHH[135]=aryuHH[12]*aryuHH[135];
   aryuHH[135]=aryuHH[216] + 1./4.*aryuHH[178] + 1./2.*aryuHH[135] + 
   aryuHH[129] + 3./4.*aryuHH[182];
   aryuHH[135]=MMH*aryuHH[135];
   aryuHH[178]=aryuHH[36] - 71./8. + aryuHH[38];
   aryuHH[178]=1./4.*aryuHH[178] + aryuHH[319];
   aryuHH[178]=aryuHH[17]*aryuHH[178];
   aryuHH[182]=3./2. + aryuHH[316];
   aryuHH[182]= - 11./2.*aryuHH[12] + 15./8.*aryuHH[7] + 11./2.*
   aryuHH[182] + aryuHH[318];
   aryuHH[182]=aryuHH[19]*aryuHH[182];
   aryuHH[195]=aryuHH[20] + aryuHH[209] - 41*aryuHH[21];
   aryuHH[197]=11./3.*aryuHH[12] + 21./4.*aryuHH[7] + 11./6.*aryuHH[36]
    - 5./2. + aryuHH[316];
   aryuHH[197]=aryuHH[18]*aryuHH[197];
   aryuHH[201]=23./8.*aryuHH[18] + aryuHH[17] - 19./8.*aryuHH[19];
   aryuHH[201]=aryuHH[11]*aryuHH[201];
   aryuHH[135]=aryuHH[135] + 1./6.*aryuHH[201] + 1./4.*aryuHH[197] + 1./
   2.*aryuHH[182] + 1./3.*aryuHH[146] + 1./8.*aryuHH[195] + aryuHH[178]
   ;
   aryuHH[135]=MMH*aryuHH[135];
   aryuHH[178]=17*aryuHH[17] - 65./3.*aryuHH[19];
   aryuHH[178]=aryuHH[19]*aryuHH[178];
   aryuHH[178]=165./4.*aryuHH[230] + aryuHH[178];
   aryuHH[182]=13./2.*aryuHH[17];
   aryuHH[195]=aryuHH[182] + 53./3.*aryuHH[19];
   aryuHH[197]=1./3.*aryuHH[18];
   aryuHH[195]=1./2.*aryuHH[195] + aryuHH[197];
   aryuHH[195]=aryuHH[18]*aryuHH[195];
   aryuHH[178]=1./4.*aryuHH[178] + aryuHH[195];
   aryuHH[135]=1./2.*aryuHH[178] + aryuHH[135];
   aryuHH[135]=aryuHH[68]*aryuHH[135];
   aryuHH[103]=aryuHH[113] + aryuHH[135] + aryuHH[103];
   aryuHH[103]=aryuHH[256]*aryuHH[103];
   aryuHH[113]= - aryuHH[7]*aryuHH[38];
   aryuHH[178]= - aryuHH[12]*aryuHH[38];
   aryuHH[195]= - aryuHH[11]*aryuHH[38];
   aryuHH[201]=1./24.*aryuHH[195] + 1./12.*aryuHH[178] + 3./8.*
   aryuHH[113] - 1./4.*aryuHH[15] - aryuHH[16] + 1./4.*aryuHH[85] - 15./
   16. + aryuHH[87];
   aryuHH[201]=MMH*aryuHH[201];
   aryuHH[213]= - 3./2.*aryuHH[7] + 5 - aryuHH[38];
   aryuHH[213]=aryuHH[17]*aryuHH[213];
   aryuHH[226]= - 11./24.*aryuHH[12] + 3./16.*aryuHH[7] + 3 + 
   aryuHH[287];
   aryuHH[226]=aryuHH[19]*aryuHH[226];
   aryuHH[231]=3 - 1./6.*aryuHH[38];
   aryuHH[233]=aryuHH[18]*aryuHH[231];
   aryuHH[236]=5*aryuHH[17] + aryuHH[19];
   aryuHH[236]=1./6.*aryuHH[236] - aryuHH[18];
   aryuHH[236]=aryuHH[11]*aryuHH[236];
   aryuHH[243]= - aryuHH[94] + 1./2.*aryuHH[71];
   aryuHH[201]=1./2.*aryuHH[201] + 1./8.*aryuHH[236] + 1./4.*
   aryuHH[233] + aryuHH[226] + 11./24.*aryuHH[146] + 1./8.*aryuHH[213]
    - 1./4.*aryuHH[20] - aryuHH[21] + 5./8.*aryuHH[44] + 1./8.*
   aryuHH[70] + aryuHH[243] - 1./4.*aryuHH[93];
   aryuHH[201]=MMH*aryuHH[201];
   aryuHH[182]=aryuHH[182] - 5*aryuHH[19];
   aryuHH[182]=aryuHH[19]*aryuHH[182];
   aryuHH[153]=aryuHH[17] + aryuHH[153];
   aryuHH[153]=1./3.*aryuHH[153] - 1./2.*aryuHH[18];
   aryuHH[153]=aryuHH[18]*aryuHH[153];
   aryuHH[153]=aryuHH[153] - 1./2.*aryuHH[230] + 1./3.*aryuHH[182];
   aryuHH[153]=1./4.*aryuHH[153] + aryuHH[201];
   aryuHH[153]=aryuHH[68]*MMH*aryuHH[153];
   aryuHH[113]=1./12.*aryuHH[195] + 1./6.*aryuHH[178] + 3./4.*
   aryuHH[113] - aryuHH[16] - 3./4. + aryuHH[87];
   aryuHH[113]=MMH*aryuHH[113];
   aryuHH[178]= - 3./4.*aryuHH[7];
   aryuHH[139]=aryuHH[178] + 1 + aryuHH[139];
   aryuHH[139]=aryuHH[17]*aryuHH[139];
   aryuHH[182]=3./8.*aryuHH[7];
   aryuHH[195]= - 5./12.*aryuHH[12] + aryuHH[231] + aryuHH[182];
   aryuHH[195]=aryuHH[19]*aryuHH[195];
   aryuHH[201]=5./12.*aryuHH[146];
   aryuHH[213]= - aryuHH[18]*aryuHH[38];
   aryuHH[226]= - aryuHH[17] + aryuHH[19];
   aryuHH[231]=aryuHH[11]*aryuHH[226];
   aryuHH[113]=1./2.*aryuHH[113] + 1./24.*aryuHH[231] + 1./12.*
   aryuHH[213] + aryuHH[195] + aryuHH[201] + 1./2.*aryuHH[139] - 
   aryuHH[21] + aryuHH[243] + 1./2.*aryuHH[44];
   aryuHH[113]=MMH*aryuHH[113];
   aryuHH[139]=7./4.*aryuHH[17] - aryuHH[19];
   aryuHH[139]=aryuHH[19]*aryuHH[139];
   aryuHH[195]=aryuHH[18]*aryuHH[226];
   aryuHH[113]=aryuHH[113] + 1./12.*aryuHH[195] - 1./4.*aryuHH[230] + 1.
   /3.*aryuHH[139];
   aryuHH[113]=aryuHH[171]*aryuHH[68]*MMH*aryuHH[113];
   aryuHH[113]=aryuHH[153] + 1./2.*aryuHH[113];
   aryuHH[113]=aryuHH[171]*aryuHH[113];
   aryuHH[139]= - 1 + aryuHH[87];
   aryuHH[153]= - aryuHH[7]*aryuHH[36];
   aryuHH[195]= - aryuHH[12]*aryuHH[36];
   aryuHH[213]= - aryuHH[11]*aryuHH[36];
   aryuHH[139]=1./12.*aryuHH[213] + 1./6.*aryuHH[195] + 3./4.*
   aryuHH[153] - aryuHH[15] + aryuHH[222] + 3*aryuHH[139] + aryuHH[85];
   aryuHH[139]=MMH*aryuHH[139];
   aryuHH[131]=aryuHH[131] + 3*aryuHH[243] - aryuHH[93];
   aryuHH[153]= - 1./6.*aryuHH[36];
   aryuHH[195]=aryuHH[142] + 9 + aryuHH[153];
   aryuHH[195]=aryuHH[19]*aryuHH[195];
   aryuHH[182]=1./12.*aryuHH[12] + aryuHH[182] + 3 + aryuHH[288];
   aryuHH[182]=aryuHH[18]*aryuHH[182];
   aryuHH[213]= - 3./16.*aryuHH[7] + 1 - 1./8.*aryuHH[36];
   aryuHH[213]=aryuHH[17]*aryuHH[213];
   aryuHH[222]=aryuHH[17] - aryuHH[18];
   aryuHH[222]=aryuHH[11]*aryuHH[222];
   aryuHH[131]=1./4.*aryuHH[139] + 11./48.*aryuHH[222] + 1./2.*
   aryuHH[182] + 1./2.*aryuHH[195] + 17./24.*aryuHH[146] + aryuHH[213]
    + aryuHH[189] - 3./2.*aryuHH[21] + 1./2.*aryuHH[131] + aryuHH[44];
   aryuHH[131]=MMH*aryuHH[131];
   aryuHH[139]=1./3.*aryuHH[17];
   aryuHH[182]=aryuHH[139] - 3./8.*aryuHH[19];
   aryuHH[182]=aryuHH[19]*aryuHH[182];
   aryuHH[195]= - 5./8.*aryuHH[18] + aryuHH[17] + 1./4.*aryuHH[19];
   aryuHH[195]=aryuHH[18]*aryuHH[195];
   aryuHH[131]=1./2.*aryuHH[131] + 1./6.*aryuHH[195] - 1./16.*
   aryuHH[230] + aryuHH[182];
   aryuHH[131]=aryuHH[68]*MMH*aryuHH[131];
   aryuHH[113]=aryuHH[131] + 1./2.*aryuHH[113];
   aryuHH[113]=aryuHH[171]*aryuHH[113];
   aryuHH[182]=aryuHH[7]*aryuHH[137];
   aryuHH[195]=aryuHH[12]*aryuHH[137];
   aryuHH[213]=aryuHH[11]*aryuHH[137];
   aryuHH[170]=1./12.*aryuHH[213] + 1./6.*aryuHH[195] + 3./4.*
   aryuHH[182] + aryuHH[170] - aryuHH[16] + aryuHH[241] - 9./8. + 
   aryuHH[87];
   aryuHH[170]=MMH*aryuHH[170];
   aryuHH[153]=aryuHH[153] + 3 + aryuHH[166];
   aryuHH[166]= - 7./12.*aryuHH[12] + aryuHH[153] - 3./8.*aryuHH[7];
   aryuHH[166]=aryuHH[19]*aryuHH[166];
   aryuHH[222]=3./4.*aryuHH[7];
   aryuHH[153]=1./6.*aryuHH[12] + aryuHH[153] + aryuHH[222];
   aryuHH[153]=aryuHH[18]*aryuHH[153];
   aryuHH[231]=3 + aryuHH[38];
   aryuHH[233]=aryuHH[231] - aryuHH[36];
   aryuHH[233]=aryuHH[17]*aryuHH[233];
   aryuHH[236]= - 5./6.*aryuHH[18] + aryuHH[17] - 1./6.*aryuHH[19];
   aryuHH[236]=aryuHH[11]*aryuHH[236];
   aryuHH[111]=1./2.*aryuHH[170] + 1./4.*aryuHH[236] + 1./2.*
   aryuHH[153] + aryuHH[166] + 1./2.*aryuHH[146] + 1./4.*aryuHH[233] + 
   aryuHH[189] - aryuHH[21] + 3./4.*aryuHH[44] + 1./4.*aryuHH[70] + 
   aryuHH[243] + aryuHH[111];
   aryuHH[111]=MMH*aryuHH[111];
   aryuHH[146]=1./8.*aryuHH[17] + aryuHH[175];
   aryuHH[146]=aryuHH[19]*aryuHH[146];
   aryuHH[153]= - 1./3.*aryuHH[18];
   aryuHH[166]=aryuHH[153] + aryuHH[17] + 1./6.*aryuHH[19];
   aryuHH[166]=aryuHH[18]*aryuHH[166];
   aryuHH[111]=1./2.*aryuHH[111] + aryuHH[146] + 1./4.*aryuHH[166];
   aryuHH[111]=aryuHH[68]*MMH*aryuHH[111];
   aryuHH[146]=aryuHH[16] - aryuHH[64] - 1./4. - aryuHH[56];
   aryuHH[166]=aryuHH[146] + 3./4.*aryuHH[308];
   aryuHH[166]=3*aryuHH[166] + 1./2.*aryuHH[285];
   aryuHH[170]=1./8.*aryuHH[290];
   aryuHH[166]=aryuHH[170] + 1./2.*aryuHH[166] + aryuHH[191];
   aryuHH[166]=MMH*aryuHH[166];
   aryuHH[189]=aryuHH[302] + 1./2.*aryuHH[21];
   aryuHH[233]= - aryuHH[17]*aryuHH[30];
   aryuHH[236]=aryuHH[189] + 1./4.*aryuHH[233];
   aryuHH[241]= - 9 - aryuHH[30];
   aryuHH[241]=aryuHH[19]*aryuHH[241];
   aryuHH[243]=MMt*aryuHH[293];
   aryuHH[252]=3./2.*aryuHH[243];
   aryuHH[261]= - 1./4.*aryuHH[30] - aryuHH[27];
   aryuHH[261]=aryuHH[18]*aryuHH[261];
   aryuHH[166]=aryuHH[252] + aryuHH[166] + aryuHH[261] + aryuHH[257] + 
   aryuHH[328] + 3*aryuHH[236] + 1./2.*aryuHH[241];
   aryuHH[166]=MMt*aryuHH[166];
   aryuHH[236]= - 9./8.*aryuHH[7];
   aryuHH[147]=aryuHH[147] + aryuHH[236] + aryuHH[36] + 3./2. - 
   aryuHH[38];
   aryuHH[147]=aryuHH[28]*aryuHH[147];
   aryuHH[241]= - 3./2.*aryuHH[45];
   aryuHH[262]=1./4.*aryuHH[301];
   aryuHH[263]=1./4.*aryuHH[18]*aryuHH[27];
   aryuHH[269]=1./8.*aryuHH[251];
   aryuHH[150]=1./4.*MMH*aryuHH[150];
   aryuHH[147]=aryuHH[150] + aryuHH[269] + aryuHH[263] + aryuHH[262] + 
   aryuHH[241] + aryuHH[147];
   aryuHH[147]=MMH*aryuHH[147];
   aryuHH[272]= - 3./4.*aryuHH[17] - aryuHH[19];
   aryuHH[272]=aryuHH[28]*aryuHH[272];
   aryuHH[147]=aryuHH[166] + aryuHH[147] + aryuHH[272] + 5./4.*
   aryuHH[266];
   aryuHH[147]=MMt*aryuHH[147];
   aryuHH[166]=aryuHH[28]*aryuHH[30];
   aryuHH[166]=aryuHH[166] + aryuHH[218];
   aryuHH[166]=aryuHH[3]*aryuHH[166];
   aryuHH[272]= - aryuHH[27]*aryuHH[30];
   aryuHH[274]=aryuHH[27]*aryuHH[30];
   aryuHH[278]=MMt*aryuHH[3]*aryuHH[274];
   aryuHH[166]=aryuHH[278] + 1./4.*aryuHH[272] + aryuHH[166];
   aryuHH[166]=MMt*aryuHH[166];
   aryuHH[166]=aryuHH[166] + 1./4.*aryuHH[265] + aryuHH[270];
   aryuHH[166]=3*aryuHH[26]*aryuHH[181]*aryuHH[166];
   aryuHH[147]=1./2.*aryuHH[147] + aryuHH[166];
   aryuHH[147]=aryuHH[26]*aryuHH[147];
   aryuHH[182]=1./6.*aryuHH[213] + 3./2.*aryuHH[182] + 1./3.*
   aryuHH[195];
   aryuHH[182]=MMH*aryuHH[182];
   aryuHH[195]= - 1./6.*aryuHH[12] + aryuHH[185] + aryuHH[178];
   aryuHH[195]=aryuHH[19]*aryuHH[195];
   aryuHH[137]=aryuHH[17]*aryuHH[137];
   aryuHH[185]=1./3.*aryuHH[12] + aryuHH[185] + 3./2.*aryuHH[7];
   aryuHH[185]=aryuHH[18]*aryuHH[185];
   aryuHH[137]=1./2.*aryuHH[182] + 1./12.*aryuHH[193] + 1./2.*
   aryuHH[185] + 1./2.*aryuHH[137] + aryuHH[195];
   aryuHH[137]=MMH*aryuHH[137];
   aryuHH[128]=aryuHH[19]*aryuHH[128];
   aryuHH[182]=aryuHH[197] + aryuHH[17] + aryuHH[291];
   aryuHH[182]=aryuHH[18]*aryuHH[182];
   aryuHH[128]=aryuHH[137] + aryuHH[128] + 1./2.*aryuHH[182];
   aryuHH[128]=aryuHH[68]*MMH*aryuHH[128];
   aryuHH[126]=aryuHH[126] + aryuHH[172] + aryuHH[236];
   aryuHH[126]=aryuHH[28]*aryuHH[126];
   aryuHH[126]=aryuHH[150] + aryuHH[269] + aryuHH[263] + aryuHH[126] + 
   aryuHH[262];
   aryuHH[126]=MMH*aryuHH[126];
   aryuHH[137]=9./2.*aryuHH[308] + aryuHH[285];
   aryuHH[137]=aryuHH[170] + 1./4.*aryuHH[137] + aryuHH[191];
   aryuHH[137]=MMH*aryuHH[137];
   aryuHH[115]=3./2.*aryuHH[233] + aryuHH[115];
   aryuHH[115]=aryuHH[137] + aryuHH[261] + 1./2.*aryuHH[115] + 
   aryuHH[257];
   aryuHH[115]=MMt*aryuHH[115];
   aryuHH[137]=aryuHH[292] + aryuHH[19];
   aryuHH[137]=aryuHH[28]*aryuHH[137];
   aryuHH[137]=aryuHH[137] + 5./2.*aryuHH[266];
   aryuHH[115]=aryuHH[115] + 1./2.*aryuHH[137] + aryuHH[126];
   aryuHH[115]=MMt*aryuHH[115];
   aryuHH[115]=1./2.*aryuHH[115] + aryuHH[166];
   aryuHH[115]=aryuHH[26]*aryuHH[115];
   aryuHH[115]=1./8.*aryuHH[128] + aryuHH[115];
   aryuHH[115]=aryuHH[256]*aryuHH[115];
   aryuHH[111]=aryuHH[115] + 1./2.*aryuHH[111] + aryuHH[147];
   aryuHH[111]=aryuHH[256]*aryuHH[111];
   aryuHH[115]=aryuHH[27]*aryuHH[36];
   aryuHH[115]=9./2.*aryuHH[146] + aryuHH[115];
   aryuHH[115]=MMH*aryuHH[115];
   aryuHH[126]= - aryuHH[28] + aryuHH[189] - 3./2.*aryuHH[19];
   aryuHH[128]=aryuHH[27]*aryuHH[17];
   aryuHH[115]=9./2.*aryuHH[243] + aryuHH[115] + aryuHH[229] + 9*
   aryuHH[126] + aryuHH[128];
   aryuHH[115]=MMt*aryuHH[115];
   aryuHH[126]=aryuHH[305] + 9./2. + aryuHH[36];
   aryuHH[126]=aryuHH[28]*aryuHH[126];
   aryuHH[128]= - MMH*aryuHH[27]*aryuHH[36];
   aryuHH[126]=1./4.*aryuHH[128] + aryuHH[263] + 1./4.*aryuHH[289] - 9./
   2.*aryuHH[45] + aryuHH[126];
   aryuHH[126]=MMH*aryuHH[126];
   aryuHH[128]=aryuHH[17] - 9./2.*aryuHH[19];
   aryuHH[128]=aryuHH[28]*aryuHH[128];
   aryuHH[115]=aryuHH[115] + aryuHH[126] + aryuHH[128] + aryuHH[266];
   aryuHH[115]=MMt*aryuHH[115];
   aryuHH[126]=aryuHH[26]*aryuHH[115];
   aryuHH[111]=aryuHH[111] + aryuHH[131] + 1./2.*aryuHH[126];
   aryuHH[111]=aryuHH[256]*aryuHH[111];
   aryuHH[126]=aryuHH[146] + 3./8.*aryuHH[119];
   aryuHH[128]=aryuHH[27]*aryuHH[38];
   aryuHH[131]=1./8.*aryuHH[125];
   aryuHH[126]=aryuHH[131] + aryuHH[128] + 3*aryuHH[126] + 1./4.*
   aryuHH[306];
   aryuHH[126]=MMH*aryuHH[126];
   aryuHH[137]=aryuHH[17]*aryuHH[30];
   aryuHH[147]=aryuHH[189] + 1./8.*aryuHH[137];
   aryuHH[150]= - 9 + aryuHH[242];
   aryuHH[150]=aryuHH[19]*aryuHH[150];
   aryuHH[166]=aryuHH[17] - aryuHH[19];
   aryuHH[166]=aryuHH[27]*aryuHH[166];
   aryuHH[170]=aryuHH[18]*aryuHH[30];
   aryuHH[126]=aryuHH[252] + 1./2.*aryuHH[126] + 1./8.*aryuHH[170] + 1./
   2.*aryuHH[166] + aryuHH[328] + 3*aryuHH[147] + 1./2.*aryuHH[150];
   aryuHH[126]=MMt*aryuHH[126];
   aryuHH[147]=9./8.*aryuHH[7];
   aryuHH[150]= - 11./4.*aryuHH[12] + aryuHH[231] + aryuHH[147];
   aryuHH[150]=aryuHH[28]*aryuHH[150];
   aryuHH[172]=1./4.*aryuHH[27]*aryuHH[226];
   aryuHH[182]=1./8.*aryuHH[217];
   aryuHH[185]= - 1./4.*MMH*aryuHH[27]*aryuHH[38];
   aryuHH[150]=aryuHH[185] + aryuHH[182] + aryuHH[172] - 3*aryuHH[45]
    + aryuHH[150];
   aryuHH[150]=MMH*aryuHH[150];
   aryuHH[173]=aryuHH[173] - aryuHH[19];
   aryuHH[173]=aryuHH[28]*aryuHH[173];
   aryuHH[173]=7*aryuHH[173] + 1./2.*aryuHH[322];
   aryuHH[150]=1./2.*aryuHH[173] + aryuHH[150];
   aryuHH[126]=1./2.*aryuHH[150] + aryuHH[126];
   aryuHH[126]=MMt*aryuHH[126];
   aryuHH[119]=aryuHH[146] + 3./4.*aryuHH[119];
   aryuHH[119]=3*aryuHH[119] + aryuHH[317];
   aryuHH[119]=aryuHH[131] + 1./2.*aryuHH[119] + aryuHH[128];
   aryuHH[119]=MMH*aryuHH[119];
   aryuHH[128]=aryuHH[189] + 1./4.*aryuHH[137];
   aryuHH[131]= - 9 + aryuHH[30];
   aryuHH[131]=aryuHH[19]*aryuHH[131];
   aryuHH[119]=aryuHH[252] + aryuHH[119] + 1./4.*aryuHH[170] + 
   aryuHH[166] + aryuHH[328] + 3*aryuHH[128] + 1./2.*aryuHH[131];
   aryuHH[119]=MMt*aryuHH[119];
   aryuHH[128]= - 5./4.*aryuHH[12] + aryuHH[147] + 3./2. + aryuHH[38];
   aryuHH[128]=aryuHH[28]*aryuHH[128];
   aryuHH[128]=aryuHH[185] + aryuHH[182] + aryuHH[172] + aryuHH[241] + 
   aryuHH[128];
   aryuHH[128]=MMH*aryuHH[128];
   aryuHH[131]=7./8.*aryuHH[17] - aryuHH[19];
   aryuHH[131]=aryuHH[28]*aryuHH[131];
   aryuHH[119]=1./2.*aryuHH[119] + 1./2.*aryuHH[128] + aryuHH[131] + 1./
   8.*aryuHH[322];
   aryuHH[119]=aryuHH[171]*MMt*aryuHH[119];
   aryuHH[119]=aryuHH[126] + aryuHH[119];
   aryuHH[119]=aryuHH[171]*aryuHH[119];
   aryuHH[115]=1./2.*aryuHH[115] + aryuHH[119];
   aryuHH[115]=aryuHH[171]*aryuHH[115];
   aryuHH[119]= - aryuHH[28]*aryuHH[30];
   aryuHH[119]=aryuHH[119] + aryuHH[265];
   aryuHH[119]=aryuHH[3]*aryuHH[119];
   aryuHH[126]=MMt*aryuHH[3]*aryuHH[272];
   aryuHH[119]=aryuHH[126] + 1./4.*aryuHH[274] + aryuHH[119];
   aryuHH[119]=MMt*aryuHH[119];
   aryuHH[119]=aryuHH[119] + 1./4.*aryuHH[218] + aryuHH[127];
   aryuHH[119]=aryuHH[181]*aryuHH[119];
   aryuHH[126]=aryuHH[171]*aryuHH[119];
   aryuHH[119]=aryuHH[119] + aryuHH[126];
   aryuHH[119]=aryuHH[26]*aryuHH[119]*pow(aryuHH[2],4);
   aryuHH[115]=aryuHH[115] + 3*aryuHH[119];
   aryuHH[115]=aryuHH[26]*aryuHH[115];
   aryuHH[111]=aryuHH[111] + aryuHH[113] + aryuHH[115];
   aryuHH[111]=aryuHH[8]*aryuHH[111];
   aryuHH[113]= - 1./2. + aryuHH[38];
   aryuHH[113]= - 1./16.*aryuHH[11] + aryuHH[136] + aryuHH[240] + 1./3.
   *aryuHH[113] - 1./4.*aryuHH[36];
   aryuHH[113]=aryuHH[11]*aryuHH[113];
   aryuHH[114]=1./3.*aryuHH[114] + aryuHH[214];
   aryuHH[114]=aryuHH[12]*aryuHH[114];
   aryuHH[115]=aryuHH[208] + 5 + aryuHH[38];
   aryuHH[115]=aryuHH[7]*aryuHH[115];
   aryuHH[113]=aryuHH[216] + 1./2.*aryuHH[113] + 1./2.*aryuHH[114] + 
   aryuHH[129] + 3./4.*aryuHH[115];
   aryuHH[113]=MMH*aryuHH[113];
   aryuHH[114]= - 71./16. + aryuHH[38];
   aryuHH[114]=1./2.*aryuHH[114] + aryuHH[319];
   aryuHH[114]=aryuHH[17]*aryuHH[114];
   aryuHH[115]=11*aryuHH[20] + 1./2.*aryuHH[209] - 17*aryuHH[21];
   aryuHH[119]= - 33./4.*aryuHH[12] - 105./8.*aryuHH[7] + 13./2. + 
   aryuHH[156];
   aryuHH[119]=aryuHH[19]*aryuHH[119];
   aryuHH[126]= - 13 - 1./3.*aryuHH[38];
   aryuHH[126]=23./4.*aryuHH[12] - 27./8.*aryuHH[7] + 1./2.*aryuHH[126]
    + aryuHH[36];
   aryuHH[126]=aryuHH[18]*aryuHH[126];
   aryuHH[127]=3./2.*aryuHH[18] + 5./3.*aryuHH[17] - 13./2.*aryuHH[19];
   aryuHH[127]=aryuHH[11]*aryuHH[127];
   aryuHH[113]=aryuHH[113] + 1./8.*aryuHH[127] + 1./2.*aryuHH[126] + 1./
   2.*aryuHH[119] + aryuHH[201] + 1./4.*aryuHH[115] + aryuHH[114];
   aryuHH[113]=MMH*aryuHH[113];
   aryuHH[114]= - 25./2.*aryuHH[17] - 71./3.*aryuHH[19];
   aryuHH[114]=aryuHH[19]*aryuHH[114];
   aryuHH[114]=165./8.*aryuHH[230] + aryuHH[114];
   aryuHH[115]= - 5./4.*aryuHH[18] - aryuHH[17] + 43./12.*aryuHH[19];
   aryuHH[115]=aryuHH[18]*aryuHH[115];
   aryuHH[114]=1./2.*aryuHH[114] + aryuHH[115];
   aryuHH[113]=1./2.*aryuHH[114] + aryuHH[113];
   aryuHH[113]=aryuHH[68]*aryuHH[113];
   aryuHH[114]=5./6.*aryuHH[12] + aryuHH[186] + aryuHH[178];
   aryuHH[114]=aryuHH[18]*aryuHH[114];
   aryuHH[115]= - 5./6.*aryuHH[12] + 1 + aryuHH[222];
   aryuHH[115]=aryuHH[19]*aryuHH[115];
   aryuHH[119]= - 1./12.*aryuHH[18] + aryuHH[139] - 1./4.*aryuHH[19];
   aryuHH[119]=aryuHH[11]*aryuHH[119];
   aryuHH[126]=MMH*aryuHH[11]*aryuHH[38];
   aryuHH[114]=1./3.*aryuHH[126] + aryuHH[119] + aryuHH[114] + 
   aryuHH[115] - aryuHH[21] + aryuHH[20];
   aryuHH[114]=MMH*aryuHH[114];
   aryuHH[115]=aryuHH[153] + aryuHH[260] + aryuHH[174];
   aryuHH[115]=aryuHH[18]*aryuHH[115];
   aryuHH[119]=1./4.*aryuHH[17] + aryuHH[175];
   aryuHH[119]=aryuHH[19]*aryuHH[119];
   aryuHH[114]=1./2.*aryuHH[114] + aryuHH[119] + 1./4.*aryuHH[115];
   aryuHH[114]=aryuHH[171]*aryuHH[68]*aryuHH[114];
   aryuHH[113]=aryuHH[113] + 1./2.*aryuHH[114];
   aryuHH[113]=aryuHH[171]*aryuHH[113];
   aryuHH[113]=aryuHH[135] + 1./2.*aryuHH[113];
   aryuHH[113]=aryuHH[171]*aryuHH[113];
   aryuHH[103]=aryuHH[111] + aryuHH[103] + aryuHH[113] + aryuHH[161];
   aryuHH[103]=aryuHH[8]*aryuHH[103];
   aryuHH[111]=aryuHH[40] + aryuHH[130] - 413./3. + aryuHH[121];
   aryuHH[111]=aryuHH[220] + aryuHH[180] + aryuHH[194] + aryuHH[183] + 
   aryuHH[177] + aryuHH[304] + 1./2.*aryuHH[111] + aryuHH[163];
   aryuHH[111]=aryuHH[27]*aryuHH[111];
   aryuHH[104]=aryuHH[167] + aryuHH[104] + aryuHH[111];
   aryuHH[104]=aryuHH[3]*aryuHH[104];
   aryuHH[104]=aryuHH[187] + aryuHH[168] + aryuHH[104];
   aryuHH[104]=MMt*aryuHH[104];
   aryuHH[111]= - 103./2. - 15*aryuHH[41];
   aryuHH[111]=aryuHH[40] + 1./2.*aryuHH[111] + aryuHH[130];
   aryuHH[111]=aryuHH[305] + aryuHH[194] + aryuHH[183] + aryuHH[177] + 
   aryuHH[304] + 1./2.*aryuHH[111] + aryuHH[163];
   aryuHH[111]=aryuHH[28]*aryuHH[111];
   aryuHH[113]= - 23 - 9./2.*aryuHH[30];
   aryuHH[113]=aryuHH[19]*aryuHH[113];
   aryuHH[111]=aryuHH[247] + aryuHH[244] + aryuHH[111] + aryuHH[224] + 
   1./2.*aryuHH[113];
   aryuHH[111]=aryuHH[3]*aryuHH[111];
   aryuHH[113]=397./3. + aryuHH[237];
   aryuHH[113]=aryuHH[152] + 1./2.*aryuHH[113] + aryuHH[248];
   aryuHH[109]=aryuHH[109] + aryuHH[154] + 1./2.*aryuHH[113] + 
   aryuHH[145];
   aryuHH[109]=aryuHH[238] + aryuHH[142] + aryuHH[184] + 1./2.*
   aryuHH[109] + aryuHH[179];
   aryuHH[109]=aryuHH[27]*aryuHH[109];
   aryuHH[113]= - 2 + aryuHH[309];
   aryuHH[113]=aryuHH[12]*aryuHH[113];
   aryuHH[104]=aryuHH[104] + aryuHH[111] + aryuHH[255] + aryuHH[249] + 
   aryuHH[221] + aryuHH[109] + aryuHH[254] + aryuHH[315] + aryuHH[219]
    + aryuHH[113];
   aryuHH[104]=MMt*aryuHH[104];
   aryuHH[109]=5./2. - 7*aryuHH[37];
   aryuHH[109]=1./3.*aryuHH[109] + aryuHH[144];
   aryuHH[109]=aryuHH[27]*aryuHH[109];
   aryuHH[111]= - 19 + aryuHH[321];
   aryuHH[111]=aryuHH[312] + 1./3.*aryuHH[111] - aryuHH[30];
   aryuHH[111]=aryuHH[11]*aryuHH[111];
   aryuHH[113]= - 25*aryuHH[55] - 97./2.*aryuHH[43] - 73 + 97./2.*
   aryuHH[58];
   aryuHH[114]= - aryuHH[67] + aryuHH[202];
   aryuHH[115]=aryuHH[114] - aryuHH[20];
   aryuHH[115]=aryuHH[99]*aryuHH[115];
   aryuHH[119]=aryuHH[28]*aryuHH[99];
   aryuHH[121]=1 + aryuHH[211];
   aryuHH[121]=aryuHH[18]*aryuHH[99]*aryuHH[121];
   aryuHH[108]=aryuHH[111] + 7./3.*aryuHH[121] + 1./2.*aryuHH[109] + 7./
   3.*aryuHH[119] + 7./6.*aryuHH[115] + 7./9.*aryuHH[29] + 25./6.*
   aryuHH[15] - 151./36.*aryuHH[63] + 1./6.*aryuHH[113] + aryuHH[108];
   aryuHH[109]= - 169./2. + 35*aryuHH[37];
   aryuHH[111]=17*aryuHH[39];
   aryuHH[109]=1./6.*aryuHH[109] + aryuHH[111];
   aryuHH[109]=aryuHH[28]*aryuHH[109];
   aryuHH[113]= - 23./3. + 5./4.*aryuHH[29];
   aryuHH[113]=197./12.*aryuHH[27] + 7./3.*aryuHH[113] + aryuHH[212];
   aryuHH[113]=aryuHH[18]*aryuHH[113];
   aryuHH[115]=3./2.*aryuHH[217];
   aryuHH[109]=aryuHH[115] + aryuHH[113] + 61./6.*aryuHH[265] + 
   aryuHH[109] + 373./36.*aryuHH[20] + 109./36.*aryuHH[45] + 41./9.*
   aryuHH[67] - 23./8.*aryuHH[33];
   aryuHH[109]=aryuHH[3]*aryuHH[109];
   aryuHH[113]= - 1 + 7*aryuHH[37];
   aryuHH[111]=aryuHH[295] + 5./6.*aryuHH[113] + aryuHH[111];
   aryuHH[111]=aryuHH[27]*aryuHH[111];
   aryuHH[113]=115./6.*aryuHH[63] + 7*aryuHH[60] + 7*aryuHH[43] + 7./2.
   *aryuHH[37] + 241./12. - 7*aryuHH[58];
   aryuHH[119]=3./2.*aryuHH[125];
   aryuHH[111]=aryuHH[119] + 1./3.*aryuHH[113] + aryuHH[111];
   aryuHH[113]=aryuHH[18] + aryuHH[114] + aryuHH[228];
   aryuHH[113]=aryuHH[3]*aryuHH[113];
   aryuHH[111]=1./2.*aryuHH[111] + 7./3.*aryuHH[113];
   aryuHH[111]=aryuHH[3]*aryuHH[111];
   aryuHH[113]= - 1 - aryuHH[63];
   aryuHH[121]=aryuHH[3]*aryuHH[113];
   aryuHH[121]=aryuHH[239] + aryuHH[121];
   aryuHH[121]=MMt*aryuHH[3]*aryuHH[121];
   aryuHH[113]=aryuHH[99]*aryuHH[113];
   aryuHH[111]=7./3.*aryuHH[121] + aryuHH[111] + 7./24.*aryuHH[113] - 
   13./8.*aryuHH[53] + 3*aryuHH[50];
   aryuHH[111]=MMt*aryuHH[111];
   aryuHH[113]=MMH*aryuHH[53];
   aryuHH[108]=aryuHH[111] + 1./2.*aryuHH[109] + 1./4.*aryuHH[108] + 1./
   3.*aryuHH[113];
   aryuHH[108]=MMt*aryuHH[108];
   aryuHH[109]=19./8. - aryuHH[13];
   aryuHH[109]=aryuHH[258] + 5./2.*aryuHH[110] + aryuHH[134] + 5./8.*
   aryuHH[63] + 1./3.*aryuHH[109] + 5./8.*aryuHH[55];
   aryuHH[110]=25*aryuHH[196];
   aryuHH[111]=101./12. + aryuHH[110];
   aryuHH[111]=aryuHH[11]*aryuHH[111];
   aryuHH[109]=5*aryuHH[109] + 1./2.*aryuHH[111];
   aryuHH[109]=MMH*aryuHH[109];
   aryuHH[111]= - 119*aryuHH[264] + 215*aryuHH[322];
   aryuHH[111]=aryuHH[3]*aryuHH[111];
   aryuHH[113]=25*aryuHH[206];
   aryuHH[121]=499./6. + aryuHH[113];
   aryuHH[121]=aryuHH[28]*aryuHH[121];
   aryuHH[125]= - 197./2.*aryuHH[27] + 1283./12. + aryuHH[110];
   aryuHH[125]=aryuHH[18]*aryuHH[125];
   aryuHH[109]=1./6.*aryuHH[111] + 1./2.*aryuHH[109] + 19*aryuHH[251]
    + 1./4.*aryuHH[125] + 61./4.*aryuHH[218] + 1./4.*aryuHH[121] - 283./
   16.*aryuHH[20] - 35./8.*aryuHH[45] - 17./3.*aryuHH[67] + 25./16.*
   aryuHH[33];
   aryuHH[111]=aryuHH[28]*aryuHH[39];
   aryuHH[111]=aryuHH[111] + aryuHH[115];
   aryuHH[111]=aryuHH[3]*aryuHH[111];
   aryuHH[115]=aryuHH[27]*aryuHH[39];
   aryuHH[115]=aryuHH[115] + aryuHH[119];
   aryuHH[115]=MMt*aryuHH[3]*aryuHH[115];
   aryuHH[119]= - aryuHH[27]*aryuHH[39];
   aryuHH[111]=aryuHH[115] + 1./4.*aryuHH[119] + aryuHH[111];
   aryuHH[111]=aryuHH[171]*MMt*aryuHH[111];
   aryuHH[108]=1./2.*aryuHH[111] + 1./12.*aryuHH[109] + aryuHH[108];
   aryuHH[108]=aryuHH[171]*aryuHH[108];
   aryuHH[109]=aryuHH[280] - 185./3.*aryuHH[28];
   aryuHH[109]=aryuHH[28]*aryuHH[109];
   aryuHH[109]=aryuHH[109] + 71./3.*aryuHH[266];
   aryuHH[109]=aryuHH[3]*aryuHH[109];
   aryuHH[111]=aryuHH[259] - 143./48. + aryuHH[12];
   aryuHH[111]=aryuHH[28]*aryuHH[111];
   aryuHH[104]=aryuHH[108] + aryuHH[104] + 1./4.*aryuHH[109] + 
   aryuHH[123] + aryuHH[268] + aryuHH[234] + aryuHH[267] + aryuHH[207]
    + aryuHH[111];
   aryuHH[104]=aryuHH[171]*aryuHH[104];
   aryuHH[108]=aryuHH[26]*aryuHH[171]*aryuHH[140];
   aryuHH[104]=aryuHH[104] + aryuHH[108];
   aryuHH[104]=aryuHH[26]*aryuHH[104];
   aryuHH[102]=aryuHH[103] + aryuHH[102] + aryuHH[105] + aryuHH[104];
   aryuHH[102]=aryuHH[8]*aryuHH[102];
   aryuHH[103]= - 61./12. - aryuHH[24];
   aryuHH[103]= - 15./8.*aryuHH[43] + 15./8.*aryuHH[58] + 1./2.*
   aryuHH[103] + aryuHH[157];
   aryuHH[103]=15./8.*aryuHH[164] + 15./4.*aryuHH[162] + 75./4.*
   aryuHH[196] + aryuHH[188] + aryuHH[120] + aryuHH[211] + aryuHH[223]
    + aryuHH[117] + aryuHH[112] + aryuHH[107] + 1./2.*aryuHH[103] + 
   aryuHH[106];
   aryuHH[103]=aryuHH[141] + 1./2.*aryuHH[103] + aryuHH[124];
   aryuHH[103]=MMZ*aryuHH[103];
   aryuHH[104]=aryuHH[232] + aryuHH[151] + 1 + aryuHH[138];
   aryuHH[104]=MMZ*aryuHH[104];
   aryuHH[104]=aryuHH[104] + aryuHH[165] + aryuHH[160] + aryuHH[271];
   aryuHH[104]=aryuHH[3]*MMZ*aryuHH[104];
   aryuHH[105]=aryuHH[28]*aryuHH[245];
   aryuHH[105]=15./4.*aryuHH[265] + aryuHH[169] + 5*aryuHH[105];
   aryuHH[103]=aryuHH[104] + aryuHH[103] + aryuHH[253] + 1./2.*
   aryuHH[105] + aryuHH[122];
   aryuHH[103]=aryuHH[3]*aryuHH[103];
   aryuHH[104]= - 17./2.*aryuHH[60] + 5*aryuHH[132] + 17./3.*aryuHH[55]
   ;
   aryuHH[105]=17./6.*aryuHH[29];
   aryuHH[104]=aryuHH[133] + aryuHH[105] + 1./2.*aryuHH[104] - 11./3.*
   aryuHH[15];
   aryuHH[104]=aryuHH[99]*aryuHH[104];
   aryuHH[106]= - aryuHH[43] + 5./6. + aryuHH[58];
   aryuHH[106]=1./3.*aryuHH[15] - 7./12.*aryuHH[63] + aryuHH[232] + 1./
   4.*aryuHH[106] - 1./3.*aryuHH[55];
   aryuHH[106]=aryuHH[97]*aryuHH[106];
   aryuHH[107]= - 5*aryuHH[9];
   aryuHH[108]=aryuHH[148] + aryuHH[107];
   aryuHH[108]=aryuHH[99]*aryuHH[108];
   aryuHH[108]=25./2.*aryuHH[97] + aryuHH[108];
   aryuHH[108]=aryuHH[11]*aryuHH[108];
   aryuHH[104]=1./6.*aryuHH[108] + 25./16.*aryuHH[204] + 25./4.*
   aryuHH[106] + aryuHH[104];
   aryuHH[104]=MMZ*aryuHH[104];
   aryuHH[106]= - 7313./144. + 5*aryuHH[22];
   aryuHH[106]= - 19./6.*aryuHH[15] - 175./24.*aryuHH[63] - 523./36.*
   aryuHH[60] + aryuHH[279] + 35./18.*aryuHH[13] + 111./8.*aryuHH[43]
    + 1./3.*aryuHH[106] - 111./8.*aryuHH[58];
   aryuHH[108]=aryuHH[107] + 73./2. + aryuHH[148];
   aryuHH[108]=aryuHH[99]*aryuHH[108];
   aryuHH[108]= - 175./8.*aryuHH[97] + aryuHH[108];
   aryuHH[108]=1./3.*aryuHH[108] + 25./4.*aryuHH[199];
   aryuHH[108]=aryuHH[18]*aryuHH[108];
   aryuHH[107]=aryuHH[107] - 101./6. + aryuHH[148];
   aryuHH[107]=1./2.*aryuHH[107] + aryuHH[113];
   aryuHH[107]=aryuHH[11]*aryuHH[107];
   aryuHH[109]=5./6.*aryuHH[20] + 1./6.*aryuHH[33] + aryuHH[118];
   aryuHH[109]=aryuHH[97]*aryuHH[109];
   aryuHH[111]=17./2.*aryuHH[33] + 5./2.*aryuHH[5] - 17*aryuHH[67];
   aryuHH[111]=1./2.*aryuHH[111];
   aryuHH[112]=aryuHH[111] - 11*aryuHH[20];
   aryuHH[112]=aryuHH[99]*aryuHH[112];
   aryuHH[113]=325./8.*aryuHH[97] + 17*aryuHH[99];
   aryuHH[113]=aryuHH[28]*aryuHH[113];
   aryuHH[115]=25./2.*aryuHH[206] - 25./8. - 17./3.*aryuHH[37];
   aryuHH[115]=aryuHH[27]*aryuHH[115];
   aryuHH[104]=25./48.*aryuHH[200] + aryuHH[104] + 1./3.*aryuHH[107] + 
   1./2.*aryuHH[108] + 1./2.*aryuHH[115] + 1./3.*aryuHH[113] + 1./3.*
   aryuHH[112] + 5./9.*aryuHH[9] + 25./8.*aryuHH[109] + 1./2.*
   aryuHH[106] + 17./9.*aryuHH[29];
   aryuHH[106]=373./12. - aryuHH[24];
   aryuHH[106]=1./2.*aryuHH[106] + aryuHH[157];
   aryuHH[106]=5./8.*aryuHH[43] + 1./3.*aryuHH[106] - 5./8.*aryuHH[58];
   aryuHH[106]= - 19./12.*aryuHH[55] + 1./2.*aryuHH[106] + aryuHH[155];
   aryuHH[107]=85*aryuHH[29];
   aryuHH[108]=25*aryuHH[9];
   aryuHH[109]=aryuHH[108] + 101./2. + aryuHH[107];
   aryuHH[109]=1./6.*aryuHH[109] + aryuHH[110];
   aryuHH[109]=aryuHH[11]*aryuHH[109];
   aryuHH[105]=25./4.*aryuHH[198] + aryuHH[109] + 25./8.*aryuHH[205] + 
   25./4.*aryuHH[203] + 125./4.*aryuHH[206] + aryuHH[133] + 25./4.*
   aryuHH[190] + aryuHH[105] + 145./12.*aryuHH[15] + 125./16.*
   aryuHH[63] + 5*aryuHH[106] + 829./72.*aryuHH[60];
   aryuHH[105]=MMZ*aryuHH[105];
   aryuHH[106]=aryuHH[108] - 1111./12. + aryuHH[107];
   aryuHH[106]=1./3.*aryuHH[106] - 25./4.*aryuHH[27];
   aryuHH[106]=aryuHH[18]*aryuHH[106];
   aryuHH[107]=17*aryuHH[37];
   aryuHH[108]= - 395./12. + aryuHH[107];
   aryuHH[108]=aryuHH[28]*aryuHH[108];
   aryuHH[108]=5*aryuHH[108] + 365./8.*aryuHH[20] + 1499./12.*
   aryuHH[45] - 289./4.*aryuHH[33] - 55./4.*aryuHH[5] + 221./3.*
   aryuHH[67];
   aryuHH[105]=aryuHH[105] + 19*aryuHH[217] + 1./2.*aryuHH[106] + 1./3.
   *aryuHH[108] + 25./4.*aryuHH[218];
   aryuHH[106]=29./4.*aryuHH[18] + aryuHH[111] + 17*aryuHH[28];
   aryuHH[108]= - 17./12.*aryuHH[60] + 5./12.*aryuHH[13] - 1 + 5./12.*
   aryuHH[24];
   aryuHH[108]=MMZ*aryuHH[108];
   aryuHH[106]=1./3.*aryuHH[106] + aryuHH[108];
   aryuHH[106]=aryuHH[3]*MMZ*aryuHH[106];
   aryuHH[105]=1./4.*aryuHH[105] + aryuHH[106];
   aryuHH[105]=aryuHH[3]*aryuHH[105];
   aryuHH[106]=aryuHH[227] + 1./2.*aryuHH[114] + aryuHH[28];
   aryuHH[108]= - 25 - 41./3.*aryuHH[60];
   aryuHH[109]= - 17./3.*aryuHH[63];
   aryuHH[108]=1./2.*aryuHH[108] + aryuHH[109];
   aryuHH[108]=MMZ*aryuHH[108];
   aryuHH[106]=41./3.*aryuHH[106] + 1./2.*aryuHH[108];
   aryuHH[106]=aryuHH[3]*aryuHH[106];
   aryuHH[108]=83./3. - 11*aryuHH[58];
   aryuHH[107]=143./4.*aryuHH[43] + 13./4.*aryuHH[108] + aryuHH[107];
   aryuHH[108]=5./4. + 17./3.*aryuHH[37];
   aryuHH[108]=aryuHH[27]*aryuHH[108];
   aryuHH[107]=5*aryuHH[108] + 7./6.*aryuHH[29] + 7./2.*aryuHH[15] + 
   511./12.*aryuHH[63] + 289./12.*aryuHH[60] + 1./3.*aryuHH[107] - 7./2.
   *aryuHH[55];
   aryuHH[108]= - 11 + 35./2.*aryuHH[29];
   aryuHH[108]=1./3.*aryuHH[108] + aryuHH[192];
   aryuHH[108]=1./4.*aryuHH[108] + aryuHH[327];
   aryuHH[108]=aryuHH[11]*aryuHH[108];
   aryuHH[110]=MMZ*aryuHH[50];
   aryuHH[106]=aryuHH[106] + 17./3.*aryuHH[110] + 1./4.*aryuHH[107] + 
   aryuHH[108];
   aryuHH[106]=aryuHH[3]*aryuHH[106];
   aryuHH[107]= - 17 + 7*aryuHH[55];
   aryuHH[107]=7./3.*aryuHH[29] - 7./3.*aryuHH[15] + aryuHH[109] + 1./3.
   *aryuHH[107] - 7./2.*aryuHH[60];
   aryuHH[107]=aryuHH[99]*aryuHH[107];
   aryuHH[108]= - aryuHH[11]*aryuHH[99]*aryuHH[29];
   aryuHH[107]=7./12.*aryuHH[108] - 25./3.*aryuHH[50] + 1./4.*
   aryuHH[107];
   aryuHH[108]=17./2.*aryuHH[53] + 7*aryuHH[50];
   aryuHH[109]= - 41./6.*aryuHH[63] - 8 - 7./6.*aryuHH[60];
   aryuHH[109]=aryuHH[3]*aryuHH[109];
   aryuHH[108]=1./3.*aryuHH[108] + aryuHH[109];
   aryuHH[108]=MMt*aryuHH[3]*aryuHH[108];
   aryuHH[106]=aryuHH[108] + 1./2.*aryuHH[107] + aryuHH[106];
   aryuHH[106]=MMt*aryuHH[106];
   aryuHH[104]=aryuHH[106] + 1./4.*aryuHH[104] + aryuHH[105];
   aryuHH[104]=aryuHH[171]*aryuHH[104];
   aryuHH[105]=1./2.*aryuHH[159] + 1./2.*aryuHH[283] + aryuHH[116] + 
   aryuHH[282];
   aryuHH[106]=aryuHH[225] + 15./8.*aryuHH[199] + 5./2.*aryuHH[215] + 
   aryuHH[158];
   aryuHH[106]=MMZ*aryuHH[106];
   aryuHH[107]=aryuHH[275] + 5*aryuHH[277];
   aryuHH[108]=15 + 13*aryuHH[60];
   aryuHH[108]=1./2.*aryuHH[108] + aryuHH[63];
   aryuHH[108]=MMZ*aryuHH[108];
   aryuHH[108]=aryuHH[273] + aryuHH[108];
   aryuHH[108]=aryuHH[3]*aryuHH[108];
   aryuHH[109]= - MMZ*aryuHH[50];
   aryuHH[107]=aryuHH[108] + 2*aryuHH[109] + 1./2.*aryuHH[107] + 
   aryuHH[246];
   aryuHH[107]=aryuHH[3]*aryuHH[107];
   aryuHH[107]=aryuHH[276] + aryuHH[235] + aryuHH[107];
   aryuHH[107]=MMt*aryuHH[107];
   aryuHH[103]=aryuHH[104] + aryuHH[107] + aryuHH[103] + aryuHH[149] + 
   1./4.*aryuHH[106] + 1./2.*aryuHH[105] + aryuHH[281];
   aryuHH[103]=aryuHH[171]*aryuHH[103];
   aryuHH[104]=aryuHH[28]*aryuHH[12];
   aryuHH[105]=MMZ*aryuHH[210];
   aryuHH[104]=8*aryuHH[104] + aryuHH[105];
   aryuHH[105]= - 1 - aryuHH[62];
   aryuHH[106]=pow(MMZ,2);
   aryuHH[105]=aryuHH[3]*aryuHH[106]*aryuHH[105];
   aryuHH[104]=4./3.*aryuHH[104] + 3*aryuHH[105];
   aryuHH[104]=aryuHH[3]*aryuHH[104];
   aryuHH[105]=MMt*aryuHH[3]*aryuHH[12]*aryuHH[176];
   aryuHH[103]=aryuHH[103] + aryuHH[104] + 32./3.*aryuHH[105];
   aryuHH[103]=aryuHH[26]*aryuHH[103];
   aryuHH[104]=3 + aryuHH[39];
   aryuHH[104]=aryuHH[98]*aryuHH[104];
   aryuHH[104]=aryuHH[104] + aryuHH[250];
   aryuHH[104]=aryuHH[12]*aryuHH[104];
   aryuHH[105]=5*aryuHH[79] - 4*aryuHH[81] - 5*aryuHH[80];
   aryuHH[107]= - 2*aryuHH[39] + 4*aryuHH[91] + 1 - 2*aryuHH[88];
   aryuHH[107]=aryuHH[98]*aryuHH[107];
   aryuHH[108]=1 + aryuHH[91];
   aryuHH[108]=aryuHH[99]*aryuHH[108];
   aryuHH[104]=6*aryuHH[104] + 3*aryuHH[108] + 2*aryuHH[105] + 3*
   aryuHH[107];
   aryuHH[104]=MMZ*aryuHH[104];
   aryuHH[105]=3 + 2*aryuHH[91];
   aryuHH[107]=1 + aryuHH[40];
   aryuHH[107]=aryuHH[12]*aryuHH[107];
   aryuHH[104]=aryuHH[104] + 3*aryuHH[105] + 4*aryuHH[107];
   aryuHH[104]=aryuHH[68]*aryuHH[104];
   aryuHH[105]=pow(CW,2);
   aryuHH[107]= - 64./3. - 3*aryuHH[105];
   aryuHH[108]= - 12*aryuHH[105];
   aryuHH[109]= - 53 + aryuHH[108];
   aryuHH[109]=aryuHH[40]*aryuHH[109];
   aryuHH[110]= - 6*aryuHH[39];
   aryuHH[107]= - 26./3.*aryuHH[12] + aryuHH[110] + 4*aryuHH[107] + 
   aryuHH[109];
   aryuHH[107]=aryuHH[12]*aryuHH[107];
   aryuHH[109]=aryuHH[105]*aryuHH[81];
   aryuHH[109]=aryuHH[81] + aryuHH[109];
   aryuHH[111]=1 + aryuHH[105];
   aryuHH[111]=aryuHH[80]*aryuHH[111];
   aryuHH[109]=2*aryuHH[109] + aryuHH[111];
   aryuHH[109]= - 8*aryuHH[77] + 4*aryuHH[109] + aryuHH[78];
   aryuHH[108]= - 41 + aryuHH[108];
   aryuHH[108]=aryuHH[79]*aryuHH[108];
   aryuHH[108]=3*aryuHH[109] + aryuHH[108];
   aryuHH[108]=MMZ*aryuHH[108];
   aryuHH[109]= - 29./3. + 3*aryuHH[88];
   aryuHH[107]=aryuHH[108] + aryuHH[107] + aryuHH[110] - 6*aryuHH[16]
    - 12*aryuHH[90] + 4*aryuHH[109] - 143./3.*aryuHH[91];
   aryuHH[107]=MMZ*aryuHH[107];
   aryuHH[108]= - aryuHH[12] - 1 - aryuHH[40];
   aryuHH[108]=aryuHH[19]*aryuHH[108];
   aryuHH[109]=1 - aryuHH[12];
   aryuHH[109]=aryuHH[18]*aryuHH[109];
   aryuHH[108]=aryuHH[109] - aryuHH[20] + aryuHH[108];
   aryuHH[107]=12*aryuHH[108] + aryuHH[107];
   aryuHH[107]=aryuHH[68]*aryuHH[107];
   aryuHH[108]=aryuHH[14] + 1 + aryuHH[25];
   aryuHH[109]=aryuHH[31]*aryuHH[108];
   aryuHH[108]=aryuHH[1]*aryuHH[108];
   aryuHH[108]=3*aryuHH[109] + aryuHH[108];
   aryuHH[106]=aryuHH[106]*aryuHH[108];
   aryuHH[108]=25 + 6*aryuHH[105];
   aryuHH[105]=47 + 12*aryuHH[105];
   aryuHH[105]=aryuHH[91]*aryuHH[105];
   aryuHH[105]=6*aryuHH[90] + aryuHH[105] + 2*aryuHH[108] - 3*
   aryuHH[92];
   aryuHH[105]=MMZ*aryuHH[105];
   aryuHH[108]= - 2*aryuHH[18] - 4*aryuHH[19] + 2*aryuHH[95] - 
   aryuHH[72];
   aryuHH[105]=6*aryuHH[108] + aryuHH[105];
   aryuHH[105]=aryuHH[68]*MMZ*aryuHH[105];
   aryuHH[105]=aryuHH[106] + 2*aryuHH[105];
   aryuHH[105]=aryuHH[3]*aryuHH[105];
   aryuHH[106]=MMZ*aryuHH[143];
   aryuHH[105]=aryuHH[105] + 4*aryuHH[106] + aryuHH[107];
   aryuHH[105]=aryuHH[3]*aryuHH[105];

      yuHHret = aryuHH[100] + aryuHH[101] + aryuHH[102] + aryuHH[103]
       + aryuHH[104] + aryuHH[105];
      return yuHHret;
}
