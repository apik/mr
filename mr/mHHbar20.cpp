#include <HH.hpp>
std::complex<long double>
HH<OS>::m20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armHHbar[333], mHHbarret;

    armHHbar[1]=double(nL + nH);
    armHHbar[2]=pow(CW,-1);
    armHHbar[3]=pow(MMH,-1);
    armHHbar[4]=pow(SW,-1);
    armHHbar[5]=Tsil::I2(0,0,MMZ,mu2);
    armHHbar[6]=Tsil::I2(0,0,MMW,mu2);
    armHHbar[7]=Tsil::B(MMH,MMH,MMH,mu2);
    armHHbar[8]=pow(MMZ,-1);
    armHHbar[9]=std::real(Tsil::B(0,0,MMZ,mu2));
    armHHbar[10]=std::real(Tsil::B(0,0,MMW,mu2));
    armHHbar[11]=Tsil::B(MMZ,MMZ,MMH,mu2);
    armHHbar[12]=Tsil::B(MMW,MMW,MMH,mu2);
    armHHbar[13]=Tsil::B(0,MMZ,MMH,mu2);
    armHHbar[14]=Tsil::B(0,MMW,MMH,mu2);
    armHHbar[15]=Tsil::Beps(MMZ,MMZ,MMH,mu2);
    armHHbar[16]=Tsil::Beps(MMW,MMW,MMH,mu2);
    armHHbar[17]=Tsil::A(MMH,mu2);
    armHHbar[18]=Tsil::A(MMZ,mu2);
    armHHbar[19]=Tsil::A(MMW,mu2);
    armHHbar[20]=Tsil::Aeps(MMZ,mu2);
    armHHbar[21]=Tsil::Aeps(MMW,mu2);
    armHHbar[22]=protZZ00->Uxzuv(0);
    armHHbar[23]=protWW00->Uxzuv(0);
    armHHbar[24]=protZZ00->Txuv(0);
    armHHbar[25]=protWW00->Txuv(0);
    armHHbar[26]=double(nH);
    armHHbar[27]=Tsil::B(MMt,MMt,MMH,mu2);
    armHHbar[28]=Tsil::A(MMt,mu2);
    armHHbar[29]=Tsil::B(MMt,MMt,MMZ,mu2);
    armHHbar[30]=Tsil::B(0,MMt,MMW,mu2);
    armHHbar[31]=double(nL);
    armHHbar[32]=Tsil::I2(MMH,MMt,MMt,mu2);
    armHHbar[33]=Tsil::I2(MMZ,MMt,MMt,mu2);
    armHHbar[34]=Tsil::I2(0,MMW,MMt,mu2);
    armHHbar[35]=Tsil::B(MMH,MMt,MMt,mu2);
    armHHbar[36]=Tsil::B(MMZ,MMH,MMZ,mu2);
    armHHbar[37]=Tsil::B(MMZ,MMt,MMt,mu2);
    armHHbar[38]=Tsil::B(MMW,MMH,MMW,mu2);
    armHHbar[39]=Tsil::B(MMW,MMZ,MMW,mu2);
    armHHbar[40]=Tsil::B(MMW,MMW,MMZ,mu2);
    armHHbar[41]=std::real(Tsil::B(0,MMW,MMt,mu2));
    armHHbar[42]=Tsil::Beps(MMH,MMH,MMH,mu2);
    armHHbar[43]=Tsil::Beps(MMt,MMt,MMH,mu2);
    armHHbar[44]=Tsil::Aeps(MMH,mu2);
    armHHbar[45]=Tsil::Aeps(MMt,mu2);
    armHHbar[46]=prottttt0->M(0);
    armHHbar[47]=prottttt0->Vxzuv(0);
    armHHbar[48]=prottttt0->Suxv(0);
    armHHbar[49]=protHtHtt->M(0);
    armHHbar[50]=protZtZtt->M(0);
    armHHbar[51]=protWtWt0->M(0);
    armHHbar[52]=protttttH->M(0);
    armHHbar[53]=protttttZ->M(0);
    armHHbar[54]=protHtHtt->Uzxyv(0);
    armHHbar[55]=protZtZtt->Uzxyv(0);
    armHHbar[56]=protWtWt0->Uzxyv(0);
    armHHbar[57]=protHtHtt->Uuyxv(0);
    armHHbar[58]=protZtZtt->Uuyxv(0);
    armHHbar[59]=protWtWt0->Uuyxv(0);
    armHHbar[60]=protZtZtt->Txuv(0);
    armHHbar[61]=protWtWt0->Suxv(0);
    armHHbar[62]=protWtWt0->Txuv(0);
    armHHbar[63]=protZtZtt->Tyzv(0);
    armHHbar[64]=protWtWt0->Tyzv(0);
    armHHbar[65]=protHtHtt->Tyzv(0);
    armHHbar[66]=protHtHtt->Svyz(0);
    armHHbar[67]=protZtZtt->Svyz(0);
    armHHbar[68]=double(boson);
    armHHbar[69]=Tsil::I2(MMH,MMH,MMH,mu2);
    armHHbar[70]=Tsil::I2(MMH,MMZ,MMZ,mu2);
    armHHbar[71]=Tsil::I2(MMH,MMW,MMW,mu2);
    armHHbar[72]=Tsil::I2(MMW,MMW,MMZ,mu2);
    armHHbar[73]=protHHHHH->M(0);
    armHHbar[74]=protHZHZZ->M(0);
    armHHbar[75]=protHWHWW->M(0);
    armHHbar[76]=protZZZZH->M(0);
    armHHbar[77]=protZWZWW->M(0);
    armHHbar[78]=protWWWWH->M(0);
    armHHbar[79]=protWWWWZ->M(0);
    armHHbar[80]=protWWWW0->M(0);
    armHHbar[81]=protWWWW0->Vxzuv(0);
    armHHbar[82]=protHHHHH->Uzxyv(0);
    armHHbar[83]=protHZHZZ->Uzxyv(0);
    armHHbar[84]=protHWHWW->Uzxyv(0);
    armHHbar[85]=protHZHZZ->Uuyxv(0);
    armHHbar[86]=protZWZWW->Uzxyv(0);
    armHHbar[87]=protHWHWW->Uuyxv(0);
    armHHbar[88]=protZWZWW->Uuyxv(0);
    armHHbar[89]=protHZHZZ->Tyzv(0);
    armHHbar[90]=protZWZWW->Txuv(0);
    armHHbar[91]=protZWZWW->Tuxv(0);
    armHHbar[92]=protHWHWW->Tuxv(0);
    armHHbar[93]=protHZHZZ->Svyz(0);
    armHHbar[94]=protHWHWW->Suxv(0);
    armHHbar[95]=protZWZWW->Suxv(0);
    armHHbar[96]=protWWWW0->Suxv(0);
    armHHbar[97]=1/(4*MMt - MMZ);
    armHHbar[98]=1/( - 4*MMW + MMH);
    armHHbar[99]=1/( - 4*MMZ + MMH);
   armHHbar[100]= - 13*armHHbar[39];
   armHHbar[101]=armHHbar[100] - 31 - 13./2.*armHHbar[40];
   armHHbar[102]=armHHbar[99]*armHHbar[17];
   armHHbar[103]=3./4.*armHHbar[102];
   armHHbar[104]= - armHHbar[18]*armHHbar[99];
   armHHbar[105]=3./2.*armHHbar[104];
   armHHbar[106]=3./8.*armHHbar[11];
   armHHbar[107]=27./4.*armHHbar[7];
   armHHbar[108]=3*armHHbar[12];
   armHHbar[109]= - 3*armHHbar[36];
   armHHbar[101]=armHHbar[106] + armHHbar[105] + armHHbar[108] + 
   armHHbar[103] + armHHbar[107] + armHHbar[109] + 1./6.*armHHbar[101]
    + armHHbar[38];
   armHHbar[101]=armHHbar[11]*armHHbar[101];
   armHHbar[110]=9*armHHbar[86];
   armHHbar[111]=143./2. + armHHbar[110];
   armHHbar[112]= - 27*armHHbar[83];
   armHHbar[111]=1./4.*armHHbar[111] + armHHbar[112];
   armHHbar[113]=9./8.*armHHbar[40];
   armHHbar[114]=3./2.*armHHbar[36];
   armHHbar[111]=27./2.*armHHbar[7] + armHHbar[114] + armHHbar[113] - 
   33./8.*armHHbar[15] - 27./16.*armHHbar[90] + 27./2.*armHHbar[42] - 1.
   /16.*armHHbar[91] + 3*armHHbar[85] + 1./2.*armHHbar[111] - 18*
   armHHbar[89];
   armHHbar[111]=armHHbar[99]*armHHbar[111];
   armHHbar[115]=3*armHHbar[74] + 1./4.*armHHbar[77];
   armHHbar[116]= - 19./2. + armHHbar[88];
   armHHbar[116]=1./2.*armHHbar[39] - 1./2.*armHHbar[16] - 9./2.*
   armHHbar[90] + 1./2.*armHHbar[116] - armHHbar[91];
   armHHbar[116]=armHHbar[98]*armHHbar[116];
   armHHbar[117]= - armHHbar[12]*armHHbar[98]*armHHbar[39];
   armHHbar[118]= - 3./4.*armHHbar[40];
   armHHbar[119]= - armHHbar[36] - 1 + armHHbar[118];
   armHHbar[119]=armHHbar[11]*armHHbar[99]*armHHbar[119];
   armHHbar[120]= - 2*armHHbar[76];
   armHHbar[121]=1./4.*armHHbar[79];
   armHHbar[111]=3./2.*armHHbar[119] + 1./8.*armHHbar[117] + 
   armHHbar[111] + 1./4.*armHHbar[116] + armHHbar[120] + 3*
   armHHbar[115] + armHHbar[121];
   armHHbar[111]=MMZ*armHHbar[111];
   armHHbar[115]= - armHHbar[95] + 1./2.*armHHbar[72];
   armHHbar[116]=3./4.*armHHbar[115] - armHHbar[93];
   armHHbar[117]=7./4.*armHHbar[70];
   armHHbar[119]=5./4.*armHHbar[44];
   armHHbar[122]=5./4.*armHHbar[17];
   armHHbar[116]=armHHbar[122] - 143./4.*armHHbar[20] + armHHbar[119]
    + 3*armHHbar[116] + armHHbar[117];
   armHHbar[116]=armHHbar[99]*armHHbar[116];
   armHHbar[123]=229./6. - 3*armHHbar[86];
   armHHbar[123]= - 17*armHHbar[88] + 1./2.*armHHbar[123] + 
   armHHbar[112];
   armHHbar[124]= - 13./3.*armHHbar[89];
   armHHbar[123]=1./2.*armHHbar[123] + armHHbar[124];
   armHHbar[125]=1 - 5./2.*armHHbar[39];
   armHHbar[126]= - 1./4.*armHHbar[12];
   armHHbar[125]=1./3.*armHHbar[125] + armHHbar[126];
   armHHbar[125]=armHHbar[12]*armHHbar[125];
   armHHbar[127]= - 27*armHHbar[7];
   armHHbar[128]=armHHbar[127] + armHHbar[109] + 43 - 9./4.*
   armHHbar[40];
   armHHbar[128]=armHHbar[99]*armHHbar[128];
   armHHbar[128]=9./2.*armHHbar[98] + armHHbar[128];
   armHHbar[128]=armHHbar[18]*armHHbar[128];
   armHHbar[129]= - 3*armHHbar[74] - 1./2.*armHHbar[77];
   armHHbar[129]= - armHHbar[76] + 3*armHHbar[129] - 7./4.*armHHbar[79]
   ;
   armHHbar[129]=MMH*armHHbar[129];
   armHHbar[130]=armHHbar[115] - armHHbar[21];
   armHHbar[130]=armHHbar[98]*armHHbar[130];
   armHHbar[131]=3 - armHHbar[39];
   armHHbar[131]=armHHbar[98]*armHHbar[131];
   armHHbar[131]=armHHbar[131] + armHHbar[99];
   armHHbar[131]=armHHbar[19]*armHHbar[131];
   armHHbar[132]=3./2.*armHHbar[39];
   armHHbar[133]=27./4.*armHHbar[42];
   armHHbar[134]= - 5./8.*armHHbar[15];
   armHHbar[135]=1./2.*armHHbar[36];
   armHHbar[101]=1./2.*armHHbar[129] + armHHbar[111] + armHHbar[101] + 
   1./2.*armHHbar[128] + 9./4.*armHHbar[131] + armHHbar[125] + 1./2.*
   armHHbar[116] + 9./4.*armHHbar[130] + armHHbar[107] + armHHbar[135]
    + armHHbar[132] + 3./4.*armHHbar[40] + armHHbar[134] + 17./4.*
   armHHbar[16] - 9./8.*armHHbar[90] + armHHbar[133] - 77./24.*
   armHHbar[91] + 1./2.*armHHbar[123] + armHHbar[85];
   armHHbar[101]=armHHbar[68]*armHHbar[101];
   armHHbar[111]=1./4.*armHHbar[115] - 3*armHHbar[93];
   armHHbar[116]=armHHbar[122] - 135./4.*armHHbar[20] + armHHbar[119]
    + armHHbar[111] + armHHbar[117];
   armHHbar[116]=armHHbar[99]*armHHbar[116];
   armHHbar[123]= - 1./4.*armHHbar[40];
   armHHbar[125]=armHHbar[127] + armHHbar[109] + 39 + armHHbar[123];
   armHHbar[125]=armHHbar[99]*armHHbar[125];
   armHHbar[125]=1./2.*armHHbar[98] + armHHbar[125];
   armHHbar[125]=armHHbar[18]*armHHbar[125];
   armHHbar[110]=163./18. + armHHbar[110];
   armHHbar[110]= - 3*armHHbar[88] + 1./2.*armHHbar[110] + 
   armHHbar[112];
   armHHbar[110]=1./2.*armHHbar[110] + armHHbar[124];
   armHHbar[128]= - armHHbar[12]*armHHbar[39];
   armHHbar[110]=1./2.*armHHbar[125] + 1./4.*armHHbar[131] + 1./12.*
   armHHbar[128] + 1./2.*armHHbar[116] + 1./4.*armHHbar[130] + 
   armHHbar[107] + armHHbar[135] + 1./6.*armHHbar[39] + 1./12.*
   armHHbar[40] - 17./8.*armHHbar[15] + 3./4.*armHHbar[16] - 37./12.*
   armHHbar[90] + armHHbar[133] - 3./4.*armHHbar[91] + 1./2.*
   armHHbar[110] + armHHbar[85];
   armHHbar[116]=17./6.*armHHbar[39] + 3 + armHHbar[123];
   armHHbar[116]=1./2.*armHHbar[116] + armHHbar[38];
   armHHbar[125]=27./8.*armHHbar[7];
   armHHbar[129]=3./8.*armHHbar[102];
   armHHbar[130]=3./4.*armHHbar[104];
   armHHbar[131]=3./16.*armHHbar[11];
   armHHbar[136]=1./4.*armHHbar[12];
   armHHbar[137]= - 2*armHHbar[36];
   armHHbar[116]=armHHbar[131] + armHHbar[130] + armHHbar[136] + 
   armHHbar[129] + armHHbar[125] + 1./2.*armHHbar[116] + armHHbar[137];
   armHHbar[116]=armHHbar[11]*armHHbar[116];
   armHHbar[138]= - 27./2.*armHHbar[83] + 9 + 1./8.*armHHbar[86];
   armHHbar[139]=1./16.*armHHbar[40];
   armHHbar[140]=3./4.*armHHbar[36];
   armHHbar[138]=armHHbar[107] + armHHbar[140] + armHHbar[139] - 25./16.
   *armHHbar[15] - 3./32.*armHHbar[90] + armHHbar[133] + 3./2.*
   armHHbar[85] + 1./2.*armHHbar[138] - 9*armHHbar[89];
   armHHbar[138]=armHHbar[99]*armHHbar[138];
   armHHbar[141]=9*armHHbar[74];
   armHHbar[142]=1./2.*armHHbar[77];
   armHHbar[143]=armHHbar[141] + armHHbar[142];
   armHHbar[144]=armHHbar[109] - 3 + armHHbar[123];
   armHHbar[144]=armHHbar[11]*armHHbar[99]*armHHbar[144];
   armHHbar[145]= - 1 - armHHbar[90];
   armHHbar[146]=armHHbar[98]*armHHbar[145];
   armHHbar[138]=1./4.*armHHbar[144] + armHHbar[138] + 1./16.*
   armHHbar[146] + 1./2.*armHHbar[143] - armHHbar[76];
   armHHbar[138]=MMZ*armHHbar[138];
   armHHbar[143]= - armHHbar[76] + armHHbar[121] - 9*armHHbar[74] + 
   armHHbar[142];
   armHHbar[143]=MMH*armHHbar[143];
   armHHbar[110]=1./4.*armHHbar[143] + armHHbar[138] + 1./2.*
   armHHbar[110] + armHHbar[116];
   armHHbar[110]=armHHbar[68]*armHHbar[110];
   armHHbar[116]=1./4.*armHHbar[40];
   armHHbar[138]=1 + armHHbar[116];
   armHHbar[143]=15*armHHbar[36];
   armHHbar[144]= - 17./2.*armHHbar[39];
   armHHbar[138]=armHHbar[11] + armHHbar[143] + 5*armHHbar[138] + 
   armHHbar[144];
   armHHbar[138]=armHHbar[11]*armHHbar[138];
   armHHbar[146]= - 43./3. - 5./2.*armHHbar[86];
   armHHbar[147]=MMZ*armHHbar[76];
   armHHbar[148]=4./3.*armHHbar[89];
   armHHbar[149]= - 3./2.*armHHbar[85];
   armHHbar[138]=3./2.*armHHbar[147] + 1./4.*armHHbar[138] + 
   armHHbar[140] + armHHbar[139] + 29./16.*armHHbar[15] + 19./96.*
   armHHbar[90] + armHHbar[149] + 1./8.*armHHbar[146] + armHHbar[148];
   armHHbar[138]=MMZ*armHHbar[138];
   armHHbar[139]=17./3.*armHHbar[93];
   armHHbar[146]= - 37./2.*armHHbar[70];
   armHHbar[147]=77./6.*armHHbar[44];
   armHHbar[150]=41./6.*armHHbar[17];
   armHHbar[151]=armHHbar[150] + 233./6.*armHHbar[20] - 5./3.*
   armHHbar[21] + armHHbar[147] + armHHbar[146] + armHHbar[139] + 7./3.
   *armHHbar[95] - 13./4.*armHHbar[72];
   armHHbar[116]=1./3. + armHHbar[116];
   armHHbar[116]=armHHbar[143] + 5*armHHbar[116] + armHHbar[144];
   armHHbar[116]=armHHbar[18]*armHHbar[116];
   armHHbar[132]=5./3. + armHHbar[132];
   armHHbar[132]=armHHbar[19]*armHHbar[132];
   armHHbar[143]=13*armHHbar[17];
   armHHbar[152]=armHHbar[18] + armHHbar[143] - 35*armHHbar[19];
   armHHbar[152]=armHHbar[11]*armHHbar[152];
   armHHbar[116]=1./2.*armHHbar[152] + armHHbar[116] + 1./2.*
   armHHbar[151] + armHHbar[132];
   armHHbar[116]=1./4.*armHHbar[116] + armHHbar[138];
   armHHbar[116]=armHHbar[68]*armHHbar[116];
   armHHbar[132]= - 1./8.*armHHbar[90] - 25./8. - 3*armHHbar[89];
   armHHbar[132]=MMZ*armHHbar[132];
   armHHbar[138]=3./2.*armHHbar[44];
   armHHbar[151]=3./2.*armHHbar[17];
   armHHbar[152]=armHHbar[11]*armHHbar[17];
   armHHbar[153]=1./2.*armHHbar[19];
   armHHbar[111]=armHHbar[132] + 3./2.*armHHbar[152] + 25./4.*
   armHHbar[18] + armHHbar[153] + armHHbar[151] + armHHbar[138] + 
   armHHbar[111] + 3./2.*armHHbar[70];
   armHHbar[111]=armHHbar[68]*MMZ*armHHbar[111];
   armHHbar[132]=armHHbar[13] + 1 + armHHbar[24];
   armHHbar[154]=armHHbar[31]*armHHbar[132];
   armHHbar[132]=armHHbar[1]*armHHbar[132];
   armHHbar[132]=11./3.*armHHbar[154] + 3*armHHbar[132];
   armHHbar[132]=MMZ*armHHbar[132];
   armHHbar[154]=armHHbar[31]*armHHbar[5];
   armHHbar[155]=armHHbar[1]*armHHbar[5];
   armHHbar[156]= - 11./3.*armHHbar[31] - 3*armHHbar[1];
   armHHbar[156]=armHHbar[18]*armHHbar[156];
   armHHbar[132]=armHHbar[132] + armHHbar[156] + 11./3.*armHHbar[154]
    + 3*armHHbar[155];
   armHHbar[132]=MMZ*armHHbar[132];
   armHHbar[111]=1./2.*armHHbar[132] + armHHbar[111];
   armHHbar[111]=armHHbar[3]*armHHbar[111];
   armHHbar[132]= - 1 - 1./3.*armHHbar[24];
   armHHbar[132]=1./2.*armHHbar[132];
   armHHbar[154]=armHHbar[132] - 5./3.*armHHbar[22];
   armHHbar[155]= - 5./9.*armHHbar[13];
   armHHbar[156]=1./6.*armHHbar[9];
   armHHbar[154]=armHHbar[156] + 5./6.*armHHbar[15] + 1./2.*
   armHHbar[154] + armHHbar[155];
   armHHbar[154]=armHHbar[31]*armHHbar[154];
   armHHbar[157]= - 5*armHHbar[22];
   armHHbar[158]= - 3 - armHHbar[24];
   armHHbar[158]=1./2.*armHHbar[158];
   armHHbar[159]=armHHbar[158] + armHHbar[157];
   armHHbar[159]=3./2.*armHHbar[9] + 15./2.*armHHbar[15] + 3./2.*
   armHHbar[159] - 5*armHHbar[13];
   armHHbar[159]=armHHbar[1]*armHHbar[159];
   armHHbar[154]=11*armHHbar[154] + armHHbar[159];
   armHHbar[159]=5./4.*armHHbar[9];
   armHHbar[160]=1 + armHHbar[159];
   armHHbar[161]=armHHbar[31]*armHHbar[160];
   armHHbar[160]=armHHbar[1]*armHHbar[160];
   armHHbar[160]=11./3.*armHHbar[161] + 3*armHHbar[160];
   armHHbar[160]=armHHbar[11]*armHHbar[160];
   armHHbar[154]=1./2.*armHHbar[154] + armHHbar[160];
   armHHbar[154]=MMZ*armHHbar[154];
   armHHbar[160]= - 11./2.*armHHbar[5] + 5*armHHbar[20];
   armHHbar[161]=armHHbar[31]*armHHbar[160];
   armHHbar[160]=armHHbar[1]*armHHbar[160];
   armHHbar[160]=11./3.*armHHbar[161] + 3*armHHbar[160];
   armHHbar[159]=2./3. + armHHbar[159];
   armHHbar[159]=armHHbar[31]*armHHbar[159];
   armHHbar[161]=2 + 15./4.*armHHbar[9];
   armHHbar[161]=armHHbar[1]*armHHbar[161];
   armHHbar[159]=11./3.*armHHbar[159] + armHHbar[161];
   armHHbar[159]=armHHbar[18]*armHHbar[159];
   armHHbar[111]=armHHbar[111] + armHHbar[116] + armHHbar[154] + 1./4.*
   armHHbar[160] + armHHbar[159];
   armHHbar[111]=armHHbar[3]*armHHbar[111];
   armHHbar[116]=7./6.*armHHbar[13];
   armHHbar[154]= - armHHbar[15] + armHHbar[116] - 55./72. + 
   armHHbar[22];
   armHHbar[159]=1./3.*armHHbar[9];
   armHHbar[154]=1./2.*armHHbar[154] + armHHbar[159];
   armHHbar[154]=armHHbar[31]*armHHbar[154];
   armHHbar[160]=1./2.*armHHbar[5];
   armHHbar[161]=armHHbar[160] - armHHbar[20];
   armHHbar[162]=armHHbar[31]*armHHbar[161];
   armHHbar[161]=armHHbar[1]*armHHbar[161];
   armHHbar[163]=11./3.*armHHbar[162] + 3*armHHbar[161];
   armHHbar[163]=armHHbar[99]*armHHbar[163];
   armHHbar[164]=1./2. - armHHbar[9];
   armHHbar[165]=armHHbar[31]*armHHbar[164];
   armHHbar[164]=armHHbar[1]*armHHbar[164];
   armHHbar[166]=11./3.*armHHbar[165] + 3*armHHbar[164];
   armHHbar[166]=armHHbar[18]*armHHbar[99]*armHHbar[166];
   armHHbar[167]= - 3*armHHbar[15] + 7./2.*armHHbar[13] - 55./24. + 3*
   armHHbar[22];
   armHHbar[167]=1./2.*armHHbar[167] + armHHbar[9];
   armHHbar[167]=armHHbar[1]*armHHbar[167];
   armHHbar[154]=1./2.*armHHbar[166] + 1./2.*armHHbar[163] + 11./3.*
   armHHbar[154] + armHHbar[167];
   armHHbar[163]= - 1./2.*armHHbar[13];
   armHHbar[132]=armHHbar[163] + armHHbar[132] + 1./3.*armHHbar[22];
   armHHbar[166]=armHHbar[159] + armHHbar[132] - 1./3.*armHHbar[15];
   armHHbar[166]=armHHbar[31]*armHHbar[166];
   armHHbar[158]= - 3./2.*armHHbar[13] + armHHbar[158] + armHHbar[22];
   armHHbar[167]=armHHbar[9] + armHHbar[158] - armHHbar[15];
   armHHbar[168]=armHHbar[1]*armHHbar[167];
   armHHbar[166]=11*armHHbar[166] + 3*armHHbar[168];
   armHHbar[166]=armHHbar[99]*armHHbar[166];
   armHHbar[169]= - armHHbar[31]*armHHbar[9];
   armHHbar[170]= - armHHbar[1]*armHHbar[9];
   armHHbar[171]=11./3.*armHHbar[169] + 3*armHHbar[170];
   armHHbar[171]=armHHbar[11]*armHHbar[99]*armHHbar[171];
   armHHbar[166]=armHHbar[166] + armHHbar[171];
   armHHbar[166]=MMZ*armHHbar[166];
   armHHbar[171]= - armHHbar[18]*armHHbar[39];
   armHHbar[172]= - MMZ*armHHbar[11]*armHHbar[39];
   armHHbar[173]=armHHbar[19] - armHHbar[18];
   armHHbar[174]=armHHbar[11]*armHHbar[173];
   armHHbar[171]=armHHbar[172] + armHHbar[171] + armHHbar[174];
   armHHbar[171]=armHHbar[3]*armHHbar[68]*armHHbar[171];
   armHHbar[172]=armHHbar[11]*armHHbar[39];
   armHHbar[145]=armHHbar[145] + 1./3.*armHHbar[172];
   armHHbar[145]=armHHbar[68]*armHHbar[145];
   armHHbar[145]=armHHbar[145] + armHHbar[171];
   armHHbar[171]=pow(armHHbar[2],2);
   armHHbar[145]=armHHbar[171]*armHHbar[145];
   armHHbar[172]= - 1./4.*armHHbar[9];
   armHHbar[175]= - 1./3. + armHHbar[172];
   armHHbar[175]=armHHbar[31]*armHHbar[175];
   armHHbar[176]= - 3./4.*armHHbar[9];
   armHHbar[177]= - 1 + armHHbar[176];
   armHHbar[177]=armHHbar[1]*armHHbar[177];
   armHHbar[175]=11./3.*armHHbar[175] + armHHbar[177];
   armHHbar[175]=armHHbar[11]*armHHbar[175];
   armHHbar[110]=1./8.*armHHbar[145] + armHHbar[111] + armHHbar[110] + 
   1./4.*armHHbar[166] + 1./2.*armHHbar[154] + armHHbar[175];
   armHHbar[110]=armHHbar[171]*armHHbar[110];
   armHHbar[111]=11*armHHbar[40];
   armHHbar[145]=13*armHHbar[39];
   armHHbar[154]= - 3*armHHbar[38];
   armHHbar[166]=armHHbar[154] + armHHbar[145] + 35 + armHHbar[111];
   armHHbar[175]=1./2.*armHHbar[11];
   armHHbar[177]=6*armHHbar[36];
   armHHbar[178]=1./2.*armHHbar[12];
   armHHbar[166]=armHHbar[175] + armHHbar[178] + 1./2.*armHHbar[166] + 
   armHHbar[177];
   armHHbar[166]=armHHbar[11]*armHHbar[166];
   armHHbar[179]=armHHbar[12]*armHHbar[39];
   armHHbar[180]=armHHbar[142] + 3*armHHbar[76];
   armHHbar[180]=MMZ*armHHbar[180];
   armHHbar[181]=8./3.*armHHbar[89];
   armHHbar[182]= - 3*armHHbar[85];
   armHHbar[113]=armHHbar[180] + armHHbar[166] + 3./8.*armHHbar[179] + 
   armHHbar[114] + 1./8.*armHHbar[39] + armHHbar[113] + 73./8.*
   armHHbar[15] + 1./8.*armHHbar[16] + 65./16.*armHHbar[90] - 7./16.*
   armHHbar[91] + armHHbar[182] + armHHbar[181] - 1./8.*armHHbar[88] + 
   2 - 49./8.*armHHbar[86];
   armHHbar[113]=MMZ*armHHbar[113];
   armHHbar[139]=armHHbar[150] + 119./2.*armHHbar[20] - 19./3.*
   armHHbar[21] + armHHbar[147] + armHHbar[146] + armHHbar[139] + 91./3.
   *armHHbar[95] - 149./4.*armHHbar[72];
   armHHbar[111]=armHHbar[154] + armHHbar[145] + 43./3. + armHHbar[111]
   ;
   armHHbar[146]=3./4.*armHHbar[12];
   armHHbar[111]=armHHbar[146] + 1./2.*armHHbar[111] + armHHbar[177];
   armHHbar[111]=armHHbar[18]*armHHbar[111];
   armHHbar[147]=17*armHHbar[18];
   armHHbar[143]=armHHbar[147] + armHHbar[143] - 19*armHHbar[19];
   armHHbar[143]=1./4.*armHHbar[11]*armHHbar[143];
   armHHbar[150]= - 3./4.*armHHbar[12] + 9./2. + 7*armHHbar[39];
   armHHbar[150]=armHHbar[19]*armHHbar[150];
   armHHbar[111]=armHHbar[113] + armHHbar[143] + armHHbar[111] + 1./4.*
   armHHbar[139] + armHHbar[150];
   armHHbar[111]=armHHbar[68]*armHHbar[111];
   armHHbar[113]= - 6*armHHbar[89];
   armHHbar[139]= - 5./2.*armHHbar[90] - 1./4.*armHHbar[91] - 35./4. + 
   armHHbar[113];
   armHHbar[139]=MMZ*armHHbar[139];
   armHHbar[150]=3*armHHbar[44];
   armHHbar[166]=3*armHHbar[17];
   armHHbar[139]=armHHbar[139] + 3*armHHbar[152] + armHHbar[147] + 10*
   armHHbar[19] + armHHbar[166] + armHHbar[150] + 3*armHHbar[70] + 5*
   armHHbar[115] - 6*armHHbar[93];
   armHHbar[139]=armHHbar[3]*armHHbar[68]*MMZ*armHHbar[139];
   armHHbar[147]= - 3./2.*armHHbar[10];
   armHHbar[179]= - 11./6.*armHHbar[9] + 10./9. + armHHbar[147];
   armHHbar[179]=armHHbar[31]*armHHbar[179];
   armHHbar[180]= - 1./2.*armHHbar[10];
   armHHbar[183]=2./3. + armHHbar[180];
   armHHbar[184]=armHHbar[183] - 3./2.*armHHbar[9];
   armHHbar[184]=armHHbar[1]*armHHbar[184];
   armHHbar[185]=armHHbar[179] + armHHbar[184];
   armHHbar[186]=armHHbar[11]*armHHbar[185];
   armHHbar[187]=MMZ*armHHbar[186];
   armHHbar[185]=armHHbar[18]*armHHbar[185];
   armHHbar[111]=armHHbar[139] + armHHbar[111] + armHHbar[185] + 
   armHHbar[187];
   armHHbar[111]=armHHbar[3]*armHHbar[111];
   armHHbar[139]=1./2.*armHHbar[10];
   armHHbar[187]= - 2./3. + armHHbar[139];
   armHHbar[188]=1./2.*armHHbar[9];
   armHHbar[187]=1./3.*armHHbar[187] + armHHbar[188];
   armHHbar[187]=armHHbar[1]*armHHbar[187];
   armHHbar[189]=11./18.*armHHbar[9] - 10./27. + armHHbar[139];
   armHHbar[189]=armHHbar[31]*armHHbar[189];
   armHHbar[187]=armHHbar[189] + armHHbar[187];
   armHHbar[187]=armHHbar[11]*armHHbar[187];
   armHHbar[101]=armHHbar[110] + armHHbar[111] + armHHbar[187] + 
   armHHbar[101];
   armHHbar[101]=armHHbar[171]*armHHbar[101];
   armHHbar[110]=593./4. - 3*armHHbar[24];
   armHHbar[110]=1./2.*armHHbar[110] - 15*armHHbar[22];
   armHHbar[111]= - 5./2.*armHHbar[13];
   armHHbar[189]= - 1./2.*armHHbar[20];
   armHHbar[190]=5*armHHbar[45] + armHHbar[189];
   armHHbar[190]=armHHbar[97]*armHHbar[190];
   armHHbar[191]=3./4.*armHHbar[9];
   armHHbar[192]=3*armHHbar[30];
   armHHbar[110]=armHHbar[191] + armHHbar[192] + 9./8.*armHHbar[190] + 
   3./4.*armHHbar[29] + 51./8.*armHHbar[15] + 45./32.*armHHbar[63] + 15
   *armHHbar[16] + 53./16.*armHHbar[60] - 21./8.*armHHbar[55] + 
   armHHbar[111] + 9./32.*armHHbar[43] - 9./32.*armHHbar[58] - 15*
   armHHbar[56] + 1./4.*armHHbar[110] + 11*armHHbar[62];
   armHHbar[193]=25*armHHbar[29];
   armHHbar[194]=39./2. + armHHbar[193];
   armHHbar[195]=17./4.*armHHbar[9];
   armHHbar[196]=armHHbar[28]*armHHbar[97];
   armHHbar[197]= - 3*armHHbar[30];
   armHHbar[194]=9./2.*armHHbar[196] + armHHbar[195] + 1./4.*
   armHHbar[194] + armHHbar[197];
   armHHbar[194]=armHHbar[11]*armHHbar[194];
   armHHbar[198]=armHHbar[43] - 3 - armHHbar[58];
   armHHbar[198]= - armHHbar[15] + 5./4.*armHHbar[63] - 1./2.*
   armHHbar[60] + 1./4.*armHHbar[198] + armHHbar[55];
   armHHbar[198]=armHHbar[97]*armHHbar[198];
   armHHbar[199]=armHHbar[27]*armHHbar[97];
   armHHbar[200]= - armHHbar[11]*armHHbar[97];
   armHHbar[198]=armHHbar[200] + armHHbar[198] + 1./4.*armHHbar[199];
   armHHbar[198]=MMZ*armHHbar[198];
   armHHbar[201]= - 3./2.*armHHbar[30];
   armHHbar[202]=7./2.*armHHbar[9] + armHHbar[201] + 6 + 11./2.*
   armHHbar[29];
   armHHbar[202]=armHHbar[12]*armHHbar[202];
   armHHbar[203]=1./4. + armHHbar[196];
   armHHbar[203]=armHHbar[27]*armHHbar[203];
   armHHbar[204]= - armHHbar[27]*armHHbar[97];
   armHHbar[205]=armHHbar[97] + armHHbar[204];
   armHHbar[205]=armHHbar[18]*armHHbar[205];
   armHHbar[206]= - armHHbar[28]*armHHbar[97];
   armHHbar[110]=9./16.*armHHbar[198] + 1./2.*armHHbar[194] + 9./32.*
   armHHbar[205] + 9./16.*armHHbar[203] + 45./16.*armHHbar[206] + 1./2.
   *armHHbar[110] + armHHbar[202];
   armHHbar[110]=MMZ*armHHbar[110];
   armHHbar[194]= - 1./4.*armHHbar[60];
   armHHbar[202]=armHHbar[194] + 1./4.*armHHbar[13] - armHHbar[62] - 1
    + 1./4.*armHHbar[24];
   armHHbar[202]=MMZ*armHHbar[202];
   armHHbar[207]=1./4.*armHHbar[33];
   armHHbar[208]= - 1./2.*armHHbar[67];
   armHHbar[209]=3*armHHbar[28];
   armHHbar[202]=armHHbar[202] + 1./4.*armHHbar[18] + armHHbar[209] + 2
   *armHHbar[19] + armHHbar[207] + armHHbar[208] - 2*armHHbar[61] + 1./
   4.*armHHbar[5] + armHHbar[34];
   armHHbar[202]=armHHbar[3]*MMZ*armHHbar[202];
   armHHbar[210]=4*armHHbar[29];
   armHHbar[211]=2*armHHbar[9];
   armHHbar[212]=3./2.*armHHbar[30];
   armHHbar[213]=armHHbar[211] + armHHbar[212] - 13./2. + armHHbar[210]
   ;
   armHHbar[213]=armHHbar[19]*armHHbar[213];
   armHHbar[193]= - 47./4. + armHHbar[193];
   armHHbar[193]= - 9./16.*armHHbar[27] + armHHbar[195] + 1./4.*
   armHHbar[193] + armHHbar[197];
   armHHbar[193]=armHHbar[18]*armHHbar[193];
   armHHbar[195]=5*armHHbar[41];
   armHHbar[214]=5./2.*armHHbar[37];
   armHHbar[215]=armHHbar[214] - 41./8. + armHHbar[195];
   armHHbar[215]=3./2.*armHHbar[215] + 5*armHHbar[12];
   armHHbar[215]=armHHbar[28]*armHHbar[215];
   armHHbar[216]= - 11./8.*armHHbar[5] - 7*armHHbar[34];
   armHHbar[216]=111./16.*armHHbar[20] + 15*armHHbar[21] + 111./8.*
   armHHbar[45] - 51./8.*armHHbar[33] + 13./2.*armHHbar[67] + 3*
   armHHbar[216] + 17*armHHbar[61];
   armHHbar[217]=armHHbar[11]*armHHbar[28];
   armHHbar[218]=armHHbar[27]*armHHbar[28];
   armHHbar[110]=3*armHHbar[202] + armHHbar[110] + 19./4.*armHHbar[217]
    + 1./2.*armHHbar[193] + 9./16.*armHHbar[218] + armHHbar[215] + 1./2.
   *armHHbar[216] + armHHbar[213];
   armHHbar[110]=armHHbar[3]*armHHbar[110];
   armHHbar[193]=5./2.*armHHbar[20];
   armHHbar[202]=1./2.*armHHbar[33];
   armHHbar[213]=armHHbar[193] + armHHbar[202] - 15*armHHbar[45];
   armHHbar[213]=armHHbar[97]*armHHbar[213];
   armHHbar[215]=1./2.*armHHbar[22];
   armHHbar[216]=armHHbar[215] - 241./32. - 3*armHHbar[59];
   armHHbar[216]=3*armHHbar[216] - 13*armHHbar[62];
   armHHbar[219]=3*armHHbar[56];
   armHHbar[220]=3./8.*armHHbar[55];
   armHHbar[221]=1./2.*armHHbar[29];
   armHHbar[222]= - 3*armHHbar[16];
   armHHbar[213]=3./16.*armHHbar[213] + armHHbar[221] - 9./8.*
   armHHbar[15] - 21./32.*armHHbar[63] + armHHbar[222] - 35./16.*
   armHHbar[60] + armHHbar[220] + 7./8.*armHHbar[13] + 207./32.*
   armHHbar[43] - 63./32.*armHHbar[58] + 1./2.*armHHbar[216] + 
   armHHbar[219];
   armHHbar[216]= - 3./2.*armHHbar[60];
   armHHbar[223]=armHHbar[216] + armHHbar[158] + armHHbar[55];
   armHHbar[223]=armHHbar[188] + armHHbar[221] + 1./2.*armHHbar[223] - 
   armHHbar[15];
   armHHbar[223]=armHHbar[99]*armHHbar[223];
   armHHbar[224]= - 3*armHHbar[43] + 5./2. + 3*armHHbar[58];
   armHHbar[225]=3./2.*armHHbar[60];
   armHHbar[224]=armHHbar[15] - 7./4.*armHHbar[63] + armHHbar[225] + 1./
   4.*armHHbar[224] - armHHbar[55];
   armHHbar[224]=armHHbar[97]*armHHbar[224];
   armHHbar[226]=armHHbar[30] - armHHbar[16] - 3./2.*armHHbar[62] + 
   armHHbar[56];
   armHHbar[226]=armHHbar[98]*armHHbar[226];
   armHHbar[227]= - armHHbar[12]*armHHbar[98]*armHHbar[30];
   armHHbar[228]= - armHHbar[29] - armHHbar[9];
   armHHbar[228]=armHHbar[99]*armHHbar[228];
   armHHbar[228]=1./2.*armHHbar[97] + armHHbar[228];
   armHHbar[228]=armHHbar[11]*armHHbar[228];
   armHHbar[223]=1./4.*armHHbar[228] + 3./32.*armHHbar[204] + 
   armHHbar[227] + 1./2.*armHHbar[223] + 1./8.*armHHbar[224] + 
   armHHbar[226];
   armHHbar[223]=MMZ*armHHbar[223];
   armHHbar[224]= - 9./2. - armHHbar[62];
   armHHbar[226]= - 1./2.*armHHbar[63];
   armHHbar[194]=armHHbar[226] + armHHbar[194] + 1./2.*armHHbar[224] - 
   armHHbar[64];
   armHHbar[194]=MMZ*armHHbar[194];
   armHHbar[224]=1./2.*armHHbar[34] - armHHbar[61];
   armHHbar[227]=1./2.*armHHbar[18];
   armHHbar[228]=2*armHHbar[28];
   armHHbar[194]=armHHbar[194] + armHHbar[227] + armHHbar[228] + 
   armHHbar[19] + armHHbar[207] + armHHbar[224] + armHHbar[208];
   armHHbar[194]=armHHbar[3]*armHHbar[194];
   armHHbar[207]= - 1./2.*armHHbar[30];
   armHHbar[229]= - 3./2.*armHHbar[56];
   armHHbar[230]=3./2.*armHHbar[16];
   armHHbar[231]=1./2.*armHHbar[64];
   armHHbar[232]=armHHbar[207] - 1./4.*armHHbar[29] - 3./4.*
   armHHbar[15] + 23./8.*armHHbar[63] + armHHbar[230] + 9./8.*
   armHHbar[60] + 3./4.*armHHbar[55] + armHHbar[231] + 15./8.*
   armHHbar[43] + 1./2.*armHHbar[37] + armHHbar[41] - 7./8.*
   armHHbar[58] + armHHbar[229] + 9./2.*armHHbar[62] + 43./8. - 
   armHHbar[59];
   armHHbar[195]=armHHbar[214] + 3./8. + armHHbar[195];
   armHHbar[195]=1./2.*armHHbar[195] + 2*armHHbar[12];
   armHHbar[195]=armHHbar[27]*armHHbar[195];
   armHHbar[214]=7 + 23./2.*armHHbar[29];
   armHHbar[233]=3*armHHbar[27];
   armHHbar[214]=1./4.*armHHbar[214] + armHHbar[233];
   armHHbar[214]=armHHbar[11]*armHHbar[214];
   armHHbar[234]= - 3./4.*armHHbar[30];
   armHHbar[235]=armHHbar[234] + 8 + 13./2.*armHHbar[29];
   armHHbar[235]=armHHbar[12]*armHHbar[235];
   armHHbar[236]=2*armHHbar[51] + armHHbar[50];
   armHHbar[236]=MMZ*armHHbar[236];
   armHHbar[194]=3*armHHbar[194] + 3*armHHbar[236] + armHHbar[214] + 3*
   armHHbar[195] + 3./2.*armHHbar[232] + armHHbar[235];
   armHHbar[194]=armHHbar[3]*armHHbar[194];
   armHHbar[195]=3./2.*armHHbar[62];
   armHHbar[214]= - armHHbar[30] + armHHbar[16] - armHHbar[64] - 
   armHHbar[56] - 1 + armHHbar[195];
   armHHbar[214]=armHHbar[98]*armHHbar[214];
   armHHbar[232]= - armHHbar[29] + armHHbar[15] - armHHbar[63] + 
   armHHbar[225] - 1 - armHHbar[55];
   armHHbar[232]=armHHbar[99]*armHHbar[232];
   armHHbar[235]=armHHbar[11]*armHHbar[99]*armHHbar[29];
   armHHbar[236]=armHHbar[12]*armHHbar[98]*armHHbar[30];
   armHHbar[237]=1./4.*armHHbar[236];
   armHHbar[214]=1./8.*armHHbar[235] + armHHbar[237] + 1./8.*
   armHHbar[232] + 1./4.*armHHbar[214] - armHHbar[51] - 1./2.*
   armHHbar[50];
   armHHbar[232]=1./2.*armHHbar[60];
   armHHbar[238]=armHHbar[226] + armHHbar[232] - 1./2.*armHHbar[64] + 1.
   /2. + armHHbar[62];
   armHHbar[238]=armHHbar[3]*armHHbar[238];
   armHHbar[239]=1./2.*armHHbar[53];
   armHHbar[238]=armHHbar[238] - armHHbar[50] - armHHbar[51] + 
   armHHbar[239];
   armHHbar[238]=MMt*armHHbar[3]*armHHbar[238];
   armHHbar[194]=3*armHHbar[238] + 3*armHHbar[214] + armHHbar[194];
   armHHbar[194]=MMt*armHHbar[194];
   armHHbar[214]= - 13./2. - 19./3.*armHHbar[29];
   armHHbar[238]= - 11./12.*armHHbar[9];
   armHHbar[214]=3./2.*armHHbar[206] + armHHbar[238] + 1./4.*
   armHHbar[214] + armHHbar[30];
   armHHbar[214]=armHHbar[11]*armHHbar[214];
   armHHbar[240]=armHHbar[202] + armHHbar[160] - armHHbar[67];
   armHHbar[240]=1./2.*armHHbar[240] - armHHbar[20];
   armHHbar[240]=armHHbar[99]*armHHbar[240];
   armHHbar[241]= - 1 - 2./3.*armHHbar[29];
   armHHbar[242]=1./2.*armHHbar[30];
   armHHbar[241]= - 2./3.*armHHbar[9] + 2*armHHbar[241] + armHHbar[242]
   ;
   armHHbar[241]=armHHbar[12]*armHHbar[241];
   armHHbar[243]= - armHHbar[9] + 5./2. - armHHbar[29];
   armHHbar[243]=armHHbar[99]*armHHbar[243];
   armHHbar[243]=3./4.*armHHbar[199] - 7./8.*armHHbar[97] + 
   armHHbar[243];
   armHHbar[243]=armHHbar[18]*armHHbar[243];
   armHHbar[244]= - armHHbar[15] + armHHbar[63] - 1 + armHHbar[55];
   armHHbar[244]=armHHbar[97]*armHHbar[244];
   armHHbar[200]=armHHbar[244] + armHHbar[200];
   armHHbar[200]=MMH*armHHbar[200];
   armHHbar[224]=armHHbar[224] - armHHbar[21];
   armHHbar[224]=armHHbar[98]*armHHbar[224];
   armHHbar[244]=1./4.*armHHbar[9];
   armHHbar[245]=1 + armHHbar[207];
   armHHbar[245]=armHHbar[19]*armHHbar[98]*armHHbar[245];
   armHHbar[246]=1./2.*armHHbar[99] + 13./16.*armHHbar[97] + 
   armHHbar[98];
   armHHbar[246]=armHHbar[28]*armHHbar[246];
   armHHbar[247]=3./4.*armHHbar[206] - 1./2.*armHHbar[37] - 3./16. - 
   armHHbar[41];
   armHHbar[247]=armHHbar[27]*armHHbar[247];
   armHHbar[110]=armHHbar[194] + armHHbar[110] + 3./64.*armHHbar[200]
    + 3./2.*armHHbar[223] + 1./2.*armHHbar[214] + 3./8.*armHHbar[243]
    + 3./4.*armHHbar[247] + 3./2.*armHHbar[246] + 3*armHHbar[245] + 
   armHHbar[241] + 3./4.*armHHbar[240] + armHHbar[244] + 3./2.*
   armHHbar[224] + 1./2.*armHHbar[213] + armHHbar[30];
   armHHbar[110]=armHHbar[26]*armHHbar[110];
   armHHbar[194]= - 55./2.*armHHbar[39] + 125./3. + 209./4.*
   armHHbar[40];
   armHHbar[194]=1./2.*armHHbar[194] + armHHbar[38];
   armHHbar[125]=armHHbar[131] + armHHbar[130] + armHHbar[146] + 
   armHHbar[129] + armHHbar[125] + 1./2.*armHHbar[194] - armHHbar[36];
   armHHbar[125]=armHHbar[11]*armHHbar[125];
   armHHbar[129]= - 3*armHHbar[84];
   armHHbar[130]=19./4. + armHHbar[129];
   armHHbar[131]= - 2*armHHbar[92];
   armHHbar[130]= - 11./8.*armHHbar[88] + 1./2.*armHHbar[130] + 
   armHHbar[131];
   armHHbar[194]=9./2.*armHHbar[42];
   armHHbar[213]=1./2.*armHHbar[38];
   armHHbar[214]=9./2.*armHHbar[7];
   armHHbar[130]=armHHbar[214] + armHHbar[213] - 33./8.*armHHbar[39] + 
   25./8.*armHHbar[16] + 33./16.*armHHbar[90] + armHHbar[194] + 33./4.*
   armHHbar[91] + 3*armHHbar[130] + armHHbar[87];
   armHHbar[130]=armHHbar[98]*armHHbar[130];
   armHHbar[223]=19 - 11*armHHbar[86];
   armHHbar[223]=1./4.*armHHbar[223] - 3*armHHbar[83];
   armHHbar[223]=1./4.*armHHbar[223] - armHHbar[89];
   armHHbar[224]=1./4.*armHHbar[36];
   armHHbar[240]=9./4.*armHHbar[7];
   armHHbar[241]=1./2.*armHHbar[85];
   armHHbar[223]=armHHbar[240] + armHHbar[224] - 33./16.*armHHbar[40]
    + 25./16.*armHHbar[15] + 99./32.*armHHbar[90] + 9./4.*armHHbar[42]
    + 33./16.*armHHbar[91] + 3*armHHbar[223] + armHHbar[241];
   armHHbar[223]=armHHbar[99]*armHHbar[223];
   armHHbar[243]=1./2.*armHHbar[74];
   armHHbar[245]=armHHbar[75] + armHHbar[243];
   armHHbar[246]= - armHHbar[38] - 1 + 33./4.*armHHbar[39];
   armHHbar[246]=armHHbar[12]*armHHbar[98]*armHHbar[246];
   armHHbar[247]= - armHHbar[36] - 1 + 33./4.*armHHbar[40];
   armHHbar[247]=armHHbar[11]*armHHbar[99]*armHHbar[247];
   armHHbar[130]=3./4.*armHHbar[247] + 3./2.*armHHbar[246] + 3*
   armHHbar[223] + 3*armHHbar[130] - armHHbar[76] + 39./4.*armHHbar[79]
    + 39./2.*armHHbar[77] + 9*armHHbar[245] - 2*armHHbar[78];
   armHHbar[130]=MMZ*armHHbar[130];
   armHHbar[223]=armHHbar[95] - 1./2.*armHHbar[72];
   armHHbar[223]=33./4.*armHHbar[223] - armHHbar[93];
   armHHbar[223]=armHHbar[122] - 35./4.*armHHbar[20] + armHHbar[119] + 
   3*armHHbar[223] + armHHbar[117];
   armHHbar[223]=armHHbar[99]*armHHbar[223];
   armHHbar[245]=99./4.*armHHbar[95];
   armHHbar[246]= - 99./8.*armHHbar[72];
   armHHbar[247]= - 35./4.*armHHbar[21] + armHHbar[119] + armHHbar[246]
    + armHHbar[245] - 3*armHHbar[94] + 7./4.*armHHbar[71];
   armHHbar[247]=armHHbar[98]*armHHbar[247];
   armHHbar[248]=armHHbar[17]*armHHbar[98];
   armHHbar[146]=armHHbar[146] + 3./4.*armHHbar[248] + armHHbar[107] - 
   armHHbar[36] - 77./8.*armHHbar[39] + 125./6. + 22*armHHbar[40];
   armHHbar[146]=armHHbar[12]*armHHbar[146];
   armHHbar[249]= - 13 + 9*armHHbar[39];
   armHHbar[249]=armHHbar[127] + 11./4.*armHHbar[249] + armHHbar[154];
   armHHbar[249]=armHHbar[98]*armHHbar[249];
   armHHbar[250]= - armHHbar[12]*armHHbar[98];
   armHHbar[249]=3*armHHbar[250] + armHHbar[249] - 99./4.*armHHbar[99];
   armHHbar[249]=armHHbar[19]*armHHbar[249];
   armHHbar[251]= - 1 + 9./4.*armHHbar[40];
   armHHbar[251]=armHHbar[127] + 11*armHHbar[251] + armHHbar[109];
   armHHbar[251]=armHHbar[99]*armHHbar[251];
   armHHbar[251]= - 99./2.*armHHbar[98] + armHHbar[251];
   armHHbar[251]=armHHbar[18]*armHHbar[251];
   armHHbar[252]= - armHHbar[75] - 1./2.*armHHbar[74];
   armHHbar[253]= - 1./2.*armHHbar[76];
   armHHbar[252]=armHHbar[253] - 11./8.*armHHbar[79] - 11./4.*
   armHHbar[77] + 9*armHHbar[252] - armHHbar[78];
   armHHbar[252]=MMH*armHHbar[252];
   armHHbar[254]=6287./24. - 27*armHHbar[84];
   armHHbar[254]=29./4.*armHHbar[88] - 27./4.*armHHbar[83] + 29./8.*
   armHHbar[86] + 1./2.*armHHbar[254] - 13./3.*armHHbar[92];
   armHHbar[125]=1./2.*armHHbar[252] + armHHbar[130] + armHHbar[125] + 
   1./4.*armHHbar[251] + 1./2.*armHHbar[249] + armHHbar[146] + 1./4.*
   armHHbar[223] + 5./8.*armHHbar[248] + 1./2.*armHHbar[247] + 81./8.*
   armHHbar[7] + armHHbar[224] + armHHbar[213] - 33./4.*armHHbar[39] - 
   33./8.*armHHbar[40] - 37./16.*armHHbar[15] - 37./8.*armHHbar[16] + 
   85./4.*armHHbar[90] + 81./8.*armHHbar[42] + 85./2.*armHHbar[91] + 
   armHHbar[241] - 13./12.*armHHbar[89] + 1./2.*armHHbar[254] + 
   armHHbar[87];
   armHHbar[125]=armHHbar[68]*armHHbar[125];
   armHHbar[130]=99./8.*armHHbar[90] + 99./4.*armHHbar[91] - 
   armHHbar[89] + 273./8. + armHHbar[131];
   armHHbar[130]=MMZ*armHHbar[130];
   armHHbar[131]=1./2.*armHHbar[70];
   armHHbar[146]=armHHbar[12]*armHHbar[17];
   armHHbar[130]=armHHbar[130] + 1./2.*armHHbar[152] - 91./4.*
   armHHbar[18] - 91./2.*armHHbar[19] + armHHbar[146] + armHHbar[151]
    + armHHbar[138] + armHHbar[131] - armHHbar[93] + armHHbar[246] + 
   armHHbar[245] - 2*armHHbar[94] + armHHbar[71];
   armHHbar[130]=armHHbar[68]*MMZ*armHHbar[130];
   armHHbar[138]=1./2.*armHHbar[24];
   armHHbar[151]=1./2.*armHHbar[13];
   armHHbar[223]=armHHbar[151] + armHHbar[138] + armHHbar[14] + 3./2.
    + armHHbar[25];
   armHHbar[224]=armHHbar[31]*armHHbar[223];
   armHHbar[223]=armHHbar[1]*armHHbar[223];
   armHHbar[223]=3*armHHbar[224] + armHHbar[223];
   armHHbar[223]=MMZ*armHHbar[223];
   armHHbar[224]=armHHbar[6] + armHHbar[160];
   armHHbar[245]=armHHbar[31]*armHHbar[224];
   armHHbar[224]=armHHbar[1]*armHHbar[224];
   armHHbar[246]= - 3*armHHbar[31] - armHHbar[1];
   armHHbar[247]=armHHbar[19]*armHHbar[246];
   armHHbar[246]=armHHbar[18]*armHHbar[246];
   armHHbar[223]=armHHbar[223] + 1./2.*armHHbar[246] + armHHbar[247] + 
   3*armHHbar[245] + armHHbar[224];
   armHHbar[223]=MMZ*armHHbar[223];
   armHHbar[130]=armHHbar[223] + 3*armHHbar[130];
   armHHbar[130]=armHHbar[3]*armHHbar[130];
   armHHbar[223]=15./4.*armHHbar[15];
   armHHbar[224]=3./2.*armHHbar[10];
   armHHbar[245]= - 1./2.*armHHbar[25];
   armHHbar[246]=armHHbar[245] - 9./4. - 5*armHHbar[23];
   armHHbar[111]=armHHbar[191] + armHHbar[224] + armHHbar[223] + 15./2.
   *armHHbar[16] + armHHbar[111] - 15./4.*armHHbar[22] - 3./8.*
   armHHbar[24] + 3./2.*armHHbar[246] - 5*armHHbar[14];
   armHHbar[111]=armHHbar[31]*armHHbar[111];
   armHHbar[191]=armHHbar[244] + armHHbar[139] + 5./4.*armHHbar[15] + 5.
   /2.*armHHbar[16] - 5./6.*armHHbar[13] - 5./4.*armHHbar[22] - 1./8.*
   armHHbar[24] + 1./2.*armHHbar[246] - 5./3.*armHHbar[14];
   armHHbar[191]=armHHbar[1]*armHHbar[191];
   armHHbar[244]=3*armHHbar[9];
   armHHbar[246]=armHHbar[244] + 2 + armHHbar[180];
   armHHbar[247]=armHHbar[31]*armHHbar[246];
   armHHbar[246]=armHHbar[1]*armHHbar[246];
   armHHbar[246]=3*armHHbar[247] + armHHbar[246];
   armHHbar[246]=armHHbar[12]*armHHbar[246];
   armHHbar[247]=7./4.*armHHbar[9];
   armHHbar[249]=1 + armHHbar[180];
   armHHbar[251]=armHHbar[249] + armHHbar[247];
   armHHbar[252]=armHHbar[31]*armHHbar[251];
   armHHbar[251]=armHHbar[1]*armHHbar[251];
   armHHbar[251]=3*armHHbar[252] + armHHbar[251];
   armHHbar[251]=armHHbar[11]*armHHbar[251];
   armHHbar[111]=armHHbar[251] + armHHbar[246] + armHHbar[111] + 
   armHHbar[191];
   armHHbar[111]=MMZ*armHHbar[111];
   armHHbar[191]= - 11./4.*armHHbar[40];
   armHHbar[246]= - 1 + armHHbar[191];
   armHHbar[251]=11./2.*armHHbar[39];
   armHHbar[252]=5*armHHbar[246] + armHHbar[251];
   armHHbar[252]= - 5*armHHbar[12] + armHHbar[114] + 5./2.*
   armHHbar[252] - armHHbar[38];
   armHHbar[252]=3*armHHbar[252] + armHHbar[175];
   armHHbar[252]=armHHbar[11]*armHHbar[252];
   armHHbar[254]= - armHHbar[38] + 77./4.*armHHbar[39] - 25 - 121./2.*
   armHHbar[40];
   armHHbar[254]=1./2.*armHHbar[254] + armHHbar[36];
   armHHbar[254]=3*armHHbar[254] - 13./4.*armHHbar[12];
   armHHbar[254]=armHHbar[12]*armHHbar[254];
   armHHbar[255]=1./2.*armHHbar[76];
   armHHbar[256]=armHHbar[255] - 33./4.*armHHbar[79] + armHHbar[78] - 
   33./2.*armHHbar[77];
   armHHbar[256]=MMZ*armHHbar[256];
   armHHbar[140]=3*armHHbar[256] + 1./2.*armHHbar[252] + armHHbar[254]
    + armHHbar[140] + 3./2.*armHHbar[38] - 99./8.*armHHbar[39] - 99./16.
   *armHHbar[40] - 207./16.*armHHbar[15] - 207./8.*armHHbar[16] - 1635./
   32.*armHHbar[90] - 1635./16.*armHHbar[91] + armHHbar[149] + 
   armHHbar[148] - 3*armHHbar[87] + 231./8.*armHHbar[88] + 231./16.*
   armHHbar[86] - 707./8. + 8./3.*armHHbar[92];
   armHHbar[140]=MMZ*armHHbar[140];
   armHHbar[148]=3 + armHHbar[191];
   armHHbar[148]=5*armHHbar[148] + armHHbar[251];
   armHHbar[148]= - 13./2.*armHHbar[12] + armHHbar[114] + 5./2.*
   armHHbar[148] - armHHbar[38];
   armHHbar[148]=armHHbar[18]*armHHbar[148];
   armHHbar[149]=13*armHHbar[146] + 41./4.*armHHbar[17] - 2059./12.*
   armHHbar[20] - 2059./6.*armHHbar[21] + 77./4.*armHHbar[44] - 37./4.*
   armHHbar[70] + 17./6.*armHHbar[93] + 1575./8.*armHHbar[72] - 423./2.
   *armHHbar[95] + 17./3.*armHHbar[94] - 37./2.*armHHbar[71];
   armHHbar[191]= - 41./4.*armHHbar[12] + armHHbar[213] + 11./8.*
   armHHbar[39] + 75./2. - 22*armHHbar[40];
   armHHbar[191]=armHHbar[19]*armHHbar[191];
   armHHbar[252]=armHHbar[17] - 15*armHHbar[19];
   armHHbar[252]=13*armHHbar[252] + 33*armHHbar[18];
   armHHbar[252]=armHHbar[11]*armHHbar[252];
   armHHbar[140]=armHHbar[140] + 1./8.*armHHbar[252] + 3./2.*
   armHHbar[148] + 1./4.*armHHbar[149] + 3*armHHbar[191];
   armHHbar[140]=armHHbar[68]*armHHbar[140];
   armHHbar[148]= - armHHbar[6] - 1./2.*armHHbar[5];
   armHHbar[148]=armHHbar[193] + 11./2.*armHHbar[148] + 5*armHHbar[21];
   armHHbar[149]=armHHbar[31]*armHHbar[148];
   armHHbar[148]=armHHbar[1]*armHHbar[148];
   armHHbar[148]=3*armHHbar[149] + armHHbar[148];
   armHHbar[149]=armHHbar[211] + 4./3. + armHHbar[139];
   armHHbar[149]=armHHbar[1]*armHHbar[149];
   armHHbar[191]=6*armHHbar[9] + 4 + armHHbar[224];
   armHHbar[191]=armHHbar[31]*armHHbar[191];
   armHHbar[149]=armHHbar[191] + armHHbar[149];
   armHHbar[149]=armHHbar[19]*armHHbar[149];
   armHHbar[183]=armHHbar[183] + armHHbar[247];
   armHHbar[183]=armHHbar[1]*armHHbar[183];
   armHHbar[191]=21./4.*armHHbar[9] + 2 + armHHbar[147];
   armHHbar[191]=armHHbar[31]*armHHbar[191];
   armHHbar[183]=armHHbar[191] + armHHbar[183];
   armHHbar[183]=armHHbar[18]*armHHbar[183];
   armHHbar[111]=armHHbar[130] + armHHbar[140] + armHHbar[111] + 
   armHHbar[183] + 1./2.*armHHbar[148] + armHHbar[149];
   armHHbar[111]=armHHbar[3]*armHHbar[111];
   armHHbar[130]=armHHbar[40] - armHHbar[39];
   armHHbar[140]=33./4.*armHHbar[130];
   armHHbar[148]= - armHHbar[36] + armHHbar[140] + armHHbar[38];
   armHHbar[149]=armHHbar[12]*armHHbar[148];
   armHHbar[183]=armHHbar[11]*armHHbar[148];
   armHHbar[149]=armHHbar[149] + 1./2.*armHHbar[183];
   armHHbar[149]=MMZ*armHHbar[149];
   armHHbar[183]=3*armHHbar[148];
   armHHbar[191]=armHHbar[183] + 17./4.*armHHbar[12];
   armHHbar[191]=armHHbar[19]*armHHbar[191];
   armHHbar[183]=armHHbar[183] - 17./2.*armHHbar[12];
   armHHbar[183]=armHHbar[18]*armHHbar[183];
   armHHbar[149]=3*armHHbar[149] + 17./8.*armHHbar[174] + armHHbar[191]
    + 1./2.*armHHbar[183];
   armHHbar[149]=armHHbar[68]*armHHbar[149];
   armHHbar[183]=armHHbar[10] - armHHbar[9];
   armHHbar[191]=armHHbar[31]*armHHbar[183];
   armHHbar[193]=armHHbar[1]*armHHbar[183];
   armHHbar[211]=3*armHHbar[191] + armHHbar[193];
   armHHbar[224]=armHHbar[12]*armHHbar[211];
   armHHbar[247]=armHHbar[11]*armHHbar[211];
   armHHbar[224]=armHHbar[224] + 1./2.*armHHbar[247];
   armHHbar[224]=MMZ*armHHbar[224];
   armHHbar[247]=armHHbar[19]*armHHbar[211];
   armHHbar[211]=armHHbar[18]*armHHbar[211];
   armHHbar[149]=armHHbar[149] + armHHbar[224] + armHHbar[247] + 1./2.*
   armHHbar[211];
   armHHbar[149]=armHHbar[3]*armHHbar[149];
   armHHbar[211]= - 1./2.*armHHbar[29];
   armHHbar[224]= - 1./2.*armHHbar[9];
   armHHbar[247]=armHHbar[224] + armHHbar[211] + armHHbar[30];
   armHHbar[252]=armHHbar[12]*armHHbar[247];
   armHHbar[254]=1./2.*armHHbar[11]*armHHbar[247];
   armHHbar[256]=armHHbar[252] + armHHbar[254];
   armHHbar[256]=MMZ*armHHbar[256];
   armHHbar[257]=armHHbar[19]*armHHbar[247];
   armHHbar[258]=1./2.*armHHbar[18]*armHHbar[247];
   armHHbar[256]=armHHbar[256] + armHHbar[257] + armHHbar[258];
   armHHbar[256]=armHHbar[3]*armHHbar[256];
   armHHbar[259]=armHHbar[188] + armHHbar[221] - armHHbar[30];
   armHHbar[260]=armHHbar[12]*armHHbar[259];
   armHHbar[259]=armHHbar[11]*armHHbar[259];
   armHHbar[261]=armHHbar[29] - armHHbar[30];
   armHHbar[262]=armHHbar[12]*armHHbar[261];
   armHHbar[263]=armHHbar[11]*armHHbar[261];
   armHHbar[264]=armHHbar[262] + 1./2.*armHHbar[263];
   armHHbar[264]=MMt*armHHbar[3]*armHHbar[264];
   armHHbar[256]=3./2.*armHHbar[264] + 3*armHHbar[256] + armHHbar[260]
    + 1./2.*armHHbar[259];
   armHHbar[256]=armHHbar[26]*armHHbar[256];
   armHHbar[259]= - armHHbar[40] + armHHbar[39];
   armHHbar[260]=2*armHHbar[36];
   armHHbar[264]= - 2*armHHbar[38];
   armHHbar[265]=armHHbar[260] + 33./4.*armHHbar[259] + armHHbar[264];
   armHHbar[265]=armHHbar[12]*armHHbar[265];
   armHHbar[266]=armHHbar[36] + 33./8.*armHHbar[259] - armHHbar[38];
   armHHbar[266]=armHHbar[11]*armHHbar[266];
   armHHbar[265]=armHHbar[265] + armHHbar[266];
   armHHbar[265]=armHHbar[68]*armHHbar[265];
   armHHbar[266]= - armHHbar[10] + armHHbar[9];
   armHHbar[267]=armHHbar[31]*armHHbar[266];
   armHHbar[268]=armHHbar[1]*armHHbar[266];
   armHHbar[269]=armHHbar[267] + 1./3.*armHHbar[268];
   armHHbar[270]=armHHbar[12]*armHHbar[269];
   armHHbar[269]=armHHbar[11]*armHHbar[269];
   armHHbar[149]=armHHbar[256] + armHHbar[149] + armHHbar[265] + 
   armHHbar[270] + 1./2.*armHHbar[269];
   armHHbar[256]=pow(armHHbar[4],2);
   armHHbar[149]=armHHbar[256]*armHHbar[149];
   armHHbar[167]=armHHbar[31]*armHHbar[167];
   armHHbar[167]=3*armHHbar[167] + armHHbar[168];
   armHHbar[167]=armHHbar[99]*armHHbar[167];
   armHHbar[168]= - armHHbar[98]*armHHbar[10];
   armHHbar[265]=armHHbar[31]*armHHbar[168];
   armHHbar[168]=armHHbar[1]*armHHbar[168];
   armHHbar[168]=3*armHHbar[265] + armHHbar[168];
   armHHbar[168]=armHHbar[12]*armHHbar[168];
   armHHbar[169]=3*armHHbar[169] + armHHbar[170];
   armHHbar[169]=armHHbar[11]*armHHbar[99]*armHHbar[169];
   armHHbar[170]=armHHbar[10] - armHHbar[16] - 3./2.*armHHbar[14] + 
   armHHbar[245] - 3./2. + armHHbar[23];
   armHHbar[170]=armHHbar[98]*armHHbar[170];
   armHHbar[245]=armHHbar[31]*armHHbar[170];
   armHHbar[170]=armHHbar[1]*armHHbar[170];
   armHHbar[167]=1./2.*armHHbar[169] + armHHbar[168] + 1./2.*
   armHHbar[167] + 3*armHHbar[245] + armHHbar[170];
   armHHbar[167]=MMZ*armHHbar[167];
   armHHbar[168]= - 3./2.*armHHbar[15] + armHHbar[222] + 7./4.*
   armHHbar[13] + 3./2.*armHHbar[22] + 7./2.*armHHbar[14] - 55./16. + 3
   *armHHbar[23];
   armHHbar[169]=1./2.*armHHbar[6] - armHHbar[21];
   armHHbar[169]=armHHbar[98]*armHHbar[169];
   armHHbar[168]=armHHbar[188] + 3./2.*armHHbar[169] + 1./2.*
   armHHbar[168] + armHHbar[10];
   armHHbar[168]=armHHbar[31]*armHHbar[168];
   armHHbar[170]= - 1./2.*armHHbar[15];
   armHHbar[215]=armHHbar[170] - armHHbar[16] + 7./12.*armHHbar[13] + 
   armHHbar[215] + 7./6.*armHHbar[14] - 55./48. + armHHbar[23];
   armHHbar[169]=armHHbar[156] + 1./2.*armHHbar[169] + 1./2.*
   armHHbar[215] + 1./3.*armHHbar[10];
   armHHbar[169]=armHHbar[1]*armHHbar[169];
   armHHbar[161]=3*armHHbar[162] + armHHbar[161];
   armHHbar[161]=armHHbar[99]*armHHbar[161];
   armHHbar[162]= - 2*armHHbar[9] - 2 + armHHbar[139];
   armHHbar[215]=armHHbar[31]*armHHbar[162];
   armHHbar[162]=armHHbar[1]*armHHbar[162];
   armHHbar[162]=armHHbar[215] + 1./3.*armHHbar[162];
   armHHbar[162]=armHHbar[12]*armHHbar[162];
   armHHbar[215]=1./2. - armHHbar[10];
   armHHbar[215]=armHHbar[98]*armHHbar[215];
   armHHbar[245]=armHHbar[31]*armHHbar[215];
   armHHbar[215]=armHHbar[1]*armHHbar[215];
   armHHbar[215]=3*armHHbar[245] + armHHbar[215];
   armHHbar[215]=armHHbar[19]*armHHbar[215];
   armHHbar[164]=3*armHHbar[165] + armHHbar[164];
   armHHbar[164]=armHHbar[18]*armHHbar[99]*armHHbar[164];
   armHHbar[139]= - 5./4.*armHHbar[9] - 1 + armHHbar[139];
   armHHbar[165]=armHHbar[31]*armHHbar[139];
   armHHbar[139]=armHHbar[1]*armHHbar[139];
   armHHbar[139]=armHHbar[165] + 1./3.*armHHbar[139];
   armHHbar[139]=armHHbar[11]*armHHbar[139];
   armHHbar[110]=armHHbar[149] + armHHbar[110] + armHHbar[111] + 
   armHHbar[125] + 1./2.*armHHbar[167] + armHHbar[139] + 1./4.*
   armHHbar[164] + 1./2.*armHHbar[215] + armHHbar[162] + 1./4.*
   armHHbar[161] + armHHbar[168] + armHHbar[169];
   armHHbar[110]=armHHbar[256]*armHHbar[110];
   armHHbar[100]=armHHbar[100] - 79 - 73./2.*armHHbar[40];
   armHHbar[100]=armHHbar[106] + armHHbar[105] + armHHbar[108] + 
   armHHbar[103] + armHHbar[107] + armHHbar[109] + 1./6.*armHHbar[100]
    + armHHbar[38];
   armHHbar[100]=armHHbar[11]*armHHbar[100];
   armHHbar[103]= - 7./2. + 11*armHHbar[86];
   armHHbar[103]=1./4.*armHHbar[103] - 9*armHHbar[83];
   armHHbar[105]= - 19./8.*armHHbar[15];
   armHHbar[103]=armHHbar[214] + armHHbar[135] + 11./8.*armHHbar[40] + 
   armHHbar[105] - 33./16.*armHHbar[90] + armHHbar[194] - 55./16.*
   armHHbar[91] + armHHbar[85] + 1./2.*armHHbar[103] + armHHbar[113];
   armHHbar[103]=armHHbar[99]*armHHbar[103];
   armHHbar[106]= - 125./8. + 9*armHHbar[84];
   armHHbar[111]=55./8.*armHHbar[39];
   armHHbar[113]= - 9./2.*armHHbar[42];
   armHHbar[125]= - 9./2.*armHHbar[7];
   armHHbar[139]= - 1./2.*armHHbar[38];
   armHHbar[106]=armHHbar[125] + armHHbar[139] + armHHbar[111] - 31./8.
   *armHHbar[16] - 11./8.*armHHbar[90] + armHHbar[113] - 55./4.*
   armHHbar[91] - armHHbar[87] + 55./8.*armHHbar[88] + 1./2.*
   armHHbar[106] + 6*armHHbar[92];
   armHHbar[106]=armHHbar[98]*armHHbar[106];
   armHHbar[149]= - 1 - 5./4.*armHHbar[39];
   armHHbar[149]=11*armHHbar[149] + armHHbar[38];
   armHHbar[149]=armHHbar[98]*armHHbar[149];
   armHHbar[161]=2*armHHbar[12]*armHHbar[98];
   armHHbar[149]=1./2.*armHHbar[149] + armHHbar[161];
   armHHbar[149]=armHHbar[12]*armHHbar[149];
   armHHbar[162]=armHHbar[246] - armHHbar[36];
   armHHbar[162]=armHHbar[11]*armHHbar[99]*armHHbar[162];
   armHHbar[103]=3./2.*armHHbar[162] + 3*armHHbar[149] + 3*
   armHHbar[103] + 3*armHHbar[106] + armHHbar[120] - 67./4.*
   armHHbar[79] - 89./4.*armHHbar[77] + 2*armHHbar[78] + armHHbar[141]
    + 10*armHHbar[80] + 8*armHHbar[81] - 9*armHHbar[75];
   armHHbar[103]=MMZ*armHHbar[103];
   armHHbar[106]=11./4.*armHHbar[115] - armHHbar[93];
   armHHbar[106]=armHHbar[122] - 167./4.*armHHbar[20] + armHHbar[119]
    + 3*armHHbar[106] + armHHbar[117];
   armHHbar[106]=armHHbar[99]*armHHbar[106];
   armHHbar[117]=5 + armHHbar[118];
   armHHbar[117]=armHHbar[127] + 11*armHHbar[117] + armHHbar[109];
   armHHbar[117]=armHHbar[99]*armHHbar[117];
   armHHbar[117]=33./2.*armHHbar[98] + armHHbar[117];
   armHHbar[117]=armHHbar[18]*armHHbar[117];
   armHHbar[118]= - 785./2. - 27*armHHbar[86];
   armHHbar[112]= - armHHbar[88] + 1./2.*armHHbar[118] + armHHbar[112];
   armHHbar[112]=1./2.*armHHbar[112] + armHHbar[124];
   armHHbar[118]=25 - 11*armHHbar[39];
   armHHbar[118]=armHHbar[98]*armHHbar[118];
   armHHbar[118]=armHHbar[118] + 11*armHHbar[99];
   armHHbar[118]=1./4.*armHHbar[118] + armHHbar[161];
   armHHbar[118]=armHHbar[19]*armHHbar[118];
   armHHbar[119]= - 2*armHHbar[81];
   armHHbar[120]=armHHbar[253] + 17./8.*armHHbar[79] - 3./4.*
   armHHbar[77] - 9./2.*armHHbar[74] + armHHbar[119] - 3*armHHbar[80];
   armHHbar[120]=MMH*armHHbar[120];
   armHHbar[122]=11*armHHbar[115] - 3*armHHbar[21];
   armHHbar[122]=armHHbar[98]*armHHbar[122];
   armHHbar[124]= - 43./4.*armHHbar[12] + 7./2.*armHHbar[39] - 91 - 53*
   armHHbar[40];
   armHHbar[124]=armHHbar[12]*armHHbar[124];
   armHHbar[100]=armHHbar[120] + armHHbar[103] + armHHbar[100] + 1./2.*
   armHHbar[117] + 3*armHHbar[118] + 1./3.*armHHbar[124] + 1./2.*
   armHHbar[106] + 3./4.*armHHbar[122] + armHHbar[107] + armHHbar[135]
    + armHHbar[251] + 11./4.*armHHbar[40] + 19./8.*armHHbar[15] + 9./4.
   *armHHbar[16] - 89./8.*armHHbar[90] + armHHbar[133] - 251./8.*
   armHHbar[91] + 1./2.*armHHbar[112] + armHHbar[85];
   armHHbar[100]=armHHbar[68]*armHHbar[100];
   armHHbar[103]= - 805./12. - armHHbar[24];
   armHHbar[103]=1./2.*armHHbar[103] + armHHbar[157];
   armHHbar[106]= - 5./3.*armHHbar[13];
   armHHbar[107]= - 5./4.*armHHbar[55];
   armHHbar[112]= - 17./24.*armHHbar[60];
   armHHbar[117]= - 75./16.*armHHbar[63];
   armHHbar[118]= - 5*armHHbar[45];
   armHHbar[120]=armHHbar[118] + 1./2.*armHHbar[20];
   armHHbar[120]=15./4.*armHHbar[97]*armHHbar[120];
   armHHbar[103]=armHHbar[188] + armHHbar[197] + armHHbar[120] + 
   armHHbar[211] + armHHbar[223] + armHHbar[117] - 15*armHHbar[16] + 
   armHHbar[112] + armHHbar[107] + armHHbar[106] - 15./16.*armHHbar[43]
    + 15./16.*armHHbar[58] + 15*armHHbar[56] + 1./2.*armHHbar[103] - 11
   *armHHbar[62];
   armHHbar[122]= - 8*armHHbar[29];
   armHHbar[124]=215./24. + armHHbar[122];
   armHHbar[133]=5./6.*armHHbar[9];
   armHHbar[124]=15./2.*armHHbar[206] + armHHbar[133] + 1./3.*
   armHHbar[124] + armHHbar[201];
   armHHbar[124]=armHHbar[11]*armHHbar[124];
   armHHbar[141]= - armHHbar[43] + 3 + armHHbar[58];
   armHHbar[141]=armHHbar[15] - 5./4.*armHHbar[63] + armHHbar[232] + 1./
   4.*armHHbar[141] - armHHbar[55];
   armHHbar[141]=armHHbar[97]*armHHbar[141];
   armHHbar[149]=armHHbar[11]*armHHbar[97];
   armHHbar[141]=armHHbar[149] + armHHbar[141] + 1./4.*armHHbar[204];
   armHHbar[141]=15./8.*MMZ*armHHbar[141];
   armHHbar[161]= - 17./3. - 14*armHHbar[29];
   armHHbar[161]= - 10./3.*armHHbar[9] + 2./3.*armHHbar[161] + 
   armHHbar[201];
   armHHbar[161]=armHHbar[12]*armHHbar[161];
   armHHbar[162]= - 1./4. + armHHbar[206];
   armHHbar[162]=armHHbar[27]*armHHbar[162];
   armHHbar[164]= - armHHbar[97] + armHHbar[199];
   armHHbar[164]=armHHbar[18]*armHHbar[164];
   armHHbar[103]=armHHbar[141] + armHHbar[124] + 15./16.*armHHbar[164]
    + 15./8.*armHHbar[162] + 75./8.*armHHbar[196] + 1./2.*armHHbar[103]
    + armHHbar[161];
   armHHbar[103]=MMZ*armHHbar[103];
   armHHbar[161]=armHHbar[232] + armHHbar[151] + 6*armHHbar[62] + 7 + 
   armHHbar[138];
   armHHbar[161]=MMZ*armHHbar[161];
   armHHbar[165]= - 3./2.*armHHbar[18];
   armHHbar[167]= - 6*armHHbar[19];
   armHHbar[168]= - 1./2.*armHHbar[33];
   armHHbar[161]=armHHbar[161] + armHHbar[165] - 8*armHHbar[28] + 
   armHHbar[167] + armHHbar[168] + armHHbar[67] + 6*armHHbar[61] + 
   armHHbar[160] - 3*armHHbar[34];
   armHHbar[161]=armHHbar[3]*MMZ*armHHbar[161];
   armHHbar[122]=385./48. + armHHbar[122];
   armHHbar[122]=15./16.*armHHbar[27] + armHHbar[133] + 1./3.*
   armHHbar[122] + armHHbar[201];
   armHHbar[122]=armHHbar[18]*armHHbar[122];
   armHHbar[169]=15./8.*armHHbar[20] - 247./12.*armHHbar[45] + 17./4.*
   armHHbar[33] - 11./4.*armHHbar[5] - 13./3.*armHHbar[67];
   armHHbar[194]= - 4*armHHbar[29];
   armHHbar[215]= - armHHbar[9] + 5./3. + armHHbar[194];
   armHHbar[215]=armHHbar[19]*armHHbar[215];
   armHHbar[245]=55./12. - armHHbar[37];
   armHHbar[246]=5./2.*armHHbar[245] - 56./3.*armHHbar[12];
   armHHbar[246]=armHHbar[28]*armHHbar[246];
   armHHbar[251]= - armHHbar[11]*armHHbar[28];
   armHHbar[253]=71./6.*armHHbar[251];
   armHHbar[265]= - armHHbar[27]*armHHbar[28];
   armHHbar[103]=armHHbar[161] + armHHbar[103] + armHHbar[253] + 
   armHHbar[122] + 15./8.*armHHbar[265] + armHHbar[246] + 1./2.*
   armHHbar[169] + 4./3.*armHHbar[215];
   armHHbar[103]=armHHbar[3]*armHHbar[103];
   armHHbar[161]= - 3*armHHbar[58];
   armHHbar[215]=3*armHHbar[43] - 5./2. + armHHbar[161];
   armHHbar[215]= - armHHbar[15] + 7./4.*armHHbar[63] + armHHbar[216]
    + 1./4.*armHHbar[215] + armHHbar[55];
   armHHbar[215]=armHHbar[97]*armHHbar[215];
   armHHbar[158]=armHHbar[9] - armHHbar[29] + armHHbar[225] + 
   armHHbar[158] - armHHbar[55];
   armHHbar[158]=armHHbar[99]*armHHbar[158];
   armHHbar[195]= - armHHbar[30] + armHHbar[16] + armHHbar[195] - 
   armHHbar[56];
   armHHbar[195]=armHHbar[98]*armHHbar[195];
   armHHbar[225]=armHHbar[29] - armHHbar[9];
   armHHbar[225]=armHHbar[99]*armHHbar[225];
   armHHbar[225]= - 5./2.*armHHbar[97] + armHHbar[225];
   armHHbar[225]=armHHbar[11]*armHHbar[225];
   armHHbar[195]=1./2.*armHHbar[225] + 15./16.*armHHbar[199] + 3*
   armHHbar[236] + 1./2.*armHHbar[158] + 5./4.*armHHbar[215] + 3*
   armHHbar[195];
   armHHbar[195]=MMZ*armHHbar[195];
   armHHbar[246]= - 59./2. - 43*armHHbar[29];
   armHHbar[269]= - 2*armHHbar[27];
   armHHbar[270]=3./4.*armHHbar[30];
   armHHbar[246]=armHHbar[269] + 1./3.*armHHbar[246] + armHHbar[270];
   armHHbar[246]=armHHbar[11]*armHHbar[246];
   armHHbar[271]= - 2*armHHbar[28];
   armHHbar[272]=armHHbar[67] + armHHbar[168];
   armHHbar[273]= - armHHbar[18] + armHHbar[272] + armHHbar[271];
   armHHbar[273]=13*armHHbar[273];
   armHHbar[274]=armHHbar[64] + 4 + 1./2.*armHHbar[62];
   armHHbar[274]=armHHbar[63] + 3*armHHbar[274] + 13./2.*armHHbar[60];
   armHHbar[274]=MMZ*armHHbar[274];
   armHHbar[274]=armHHbar[273] + armHHbar[274];
   armHHbar[274]=armHHbar[3]*armHHbar[274];
   armHHbar[275]= - 307./3. + 19*armHHbar[58];
   armHHbar[275]= - 11./2.*armHHbar[29] - 33./2.*armHHbar[15] - 83./4.*
   armHHbar[63] - 77./4.*armHHbar[60] + 33./2.*armHHbar[55] - 19./4.*
   armHHbar[43] + 1./4.*armHHbar[275] - armHHbar[37];
   armHHbar[276]= - 1 - armHHbar[29];
   armHHbar[276]=armHHbar[12]*armHHbar[276];
   armHHbar[277]= - 3./4. - armHHbar[37];
   armHHbar[277]=armHHbar[27]*armHHbar[277];
   armHHbar[278]= - 3*armHHbar[51];
   armHHbar[279]=armHHbar[278] - armHHbar[50];
   armHHbar[279]=MMZ*armHHbar[279];
   armHHbar[274]=armHHbar[274] + 2*armHHbar[279] + armHHbar[246] + 5./2.
   *armHHbar[277] + 1./2.*armHHbar[275] + 56./3.*armHHbar[276];
   armHHbar[274]=armHHbar[3]*armHHbar[274];
   armHHbar[276]= - 11*armHHbar[29] + 11*armHHbar[15] + armHHbar[63] + 
   33./2.*armHHbar[60] + 1 - 11*armHHbar[55];
   armHHbar[276]=armHHbar[99]*armHHbar[276];
   armHHbar[235]=11./4.*armHHbar[235] + 5*armHHbar[50] + 1./4.*
   armHHbar[276];
   armHHbar[276]=13*armHHbar[63] + 24 + 11*armHHbar[60];
   armHHbar[276]=armHHbar[3]*armHHbar[276];
   armHHbar[276]=armHHbar[276] - armHHbar[53] - 22*armHHbar[50];
   armHHbar[276]=MMt*armHHbar[3]*armHHbar[276];
   armHHbar[274]=armHHbar[276] + armHHbar[235] + armHHbar[274];
   armHHbar[274]=MMt*armHHbar[274];
   armHHbar[279]=3./2.*armHHbar[55];
   armHHbar[280]= - 5./2.*armHHbar[15];
   armHHbar[116]=armHHbar[280] + 35./8.*armHHbar[63] - 1./12.*
   armHHbar[60] + armHHbar[279] + armHHbar[116] - 9./8.*armHHbar[43] + 
   9./8.*armHHbar[58] + 43./48. + armHHbar[22];
   armHHbar[160]=armHHbar[168] + armHHbar[160] + armHHbar[67];
   armHHbar[281]=armHHbar[99]*armHHbar[160];
   armHHbar[282]= - 5./2.*armHHbar[20] + armHHbar[168] + 15*
   armHHbar[45];
   armHHbar[282]=armHHbar[97]*armHHbar[282];
   armHHbar[116]=1./2.*armHHbar[281] + armHHbar[159] + 5./8.*
   armHHbar[282] + 1./2.*armHHbar[116] - 1./3.*armHHbar[29];
   armHHbar[159]= - armHHbar[9] - 3./2. + armHHbar[29];
   armHHbar[159]=armHHbar[99]*armHHbar[159];
   armHHbar[159]=15./4.*armHHbar[204] + 35./8.*armHHbar[97] + 
   armHHbar[159];
   armHHbar[159]=armHHbar[18]*armHHbar[159];
   armHHbar[281]= - 215./12. + 13*armHHbar[29];
   armHHbar[281]=1./9.*armHHbar[281] + armHHbar[30];
   armHHbar[281]=5./2.*armHHbar[196] + 1./2.*armHHbar[281] - 1./9.*
   armHHbar[9];
   armHHbar[281]=armHHbar[11]*armHHbar[281];
   armHHbar[282]=armHHbar[15] - armHHbar[63] + 1 - armHHbar[55];
   armHHbar[282]=armHHbar[97]*armHHbar[282];
   armHHbar[149]=armHHbar[282] + armHHbar[149];
   armHHbar[149]=5./32.*MMH*armHHbar[149];
   armHHbar[210]=armHHbar[9] - 5./3. + armHHbar[210];
   armHHbar[210]=armHHbar[12]*armHHbar[210];
   armHHbar[282]= - 65./8.*armHHbar[97] - armHHbar[99];
   armHHbar[282]=armHHbar[28]*armHHbar[282];
   armHHbar[283]=15./2.*armHHbar[196] + 15./8. + armHHbar[37];
   armHHbar[283]=armHHbar[27]*armHHbar[283];
   armHHbar[103]=armHHbar[274] + armHHbar[103] + armHHbar[149] + 1./2.*
   armHHbar[195] + armHHbar[281] + 1./4.*armHHbar[159] + 1./4.*
   armHHbar[283] + 1./2.*armHHbar[282] + 1./2.*armHHbar[116] + 4./9.*
   armHHbar[210];
   armHHbar[103]=armHHbar[26]*armHHbar[103];
   armHHbar[195]=53*armHHbar[40];
   armHHbar[274]=armHHbar[154] + armHHbar[145] + 59 + armHHbar[195];
   armHHbar[175]=armHHbar[175] + 25./2.*armHHbar[12] + 1./2.*
   armHHbar[274] + armHHbar[177];
   armHHbar[175]=armHHbar[11]*armHHbar[175];
   armHHbar[274]=713./2. - 8*armHHbar[92];
   armHHbar[284]= - 3./2.*armHHbar[38];
   armHHbar[285]=35./3.*armHHbar[12] + armHHbar[284] + 71./8.*
   armHHbar[39] + 737./6. + 119*armHHbar[40];
   armHHbar[285]=armHHbar[12]*armHHbar[285];
   armHHbar[119]=armHHbar[119] - armHHbar[80];
   armHHbar[119]=2*armHHbar[119] - armHHbar[78];
   armHHbar[119]=armHHbar[76] + 22*armHHbar[79] + 2*armHHbar[119] + 55./
   2.*armHHbar[77];
   armHHbar[119]=MMZ*armHHbar[119];
   armHHbar[114]=3*armHHbar[119] + armHHbar[175] + armHHbar[285] + 
   armHHbar[114] + armHHbar[284] + 165./8.*armHHbar[39] + 33./8.*
   armHHbar[40] + 145./8.*armHHbar[15] + 269./8.*armHHbar[16] + 953./16.
   *armHHbar[90] + 2317./16.*armHHbar[91] + armHHbar[182] + 
   armHHbar[181] + 3*armHHbar[87] - 341./8.*armHHbar[88] + 1./3.*
   armHHbar[274] - 121./8.*armHHbar[86];
   armHHbar[114]=MMZ*armHHbar[114];
   armHHbar[119]=armHHbar[154] + armHHbar[145] - 105 + armHHbar[195];
   armHHbar[119]=155./4.*armHHbar[12] + 1./2.*armHHbar[119] + 
   armHHbar[177];
   armHHbar[119]=armHHbar[18]*armHHbar[119];
   armHHbar[175]=439./12.*armHHbar[12] + armHHbar[145] + 55./6. + 
   armHHbar[195];
   armHHbar[175]=armHHbar[19]*armHHbar[175];
   armHHbar[114]=armHHbar[114] + armHHbar[143] + armHHbar[119] + 
   armHHbar[175] + 41./24.*armHHbar[17] + 1589./24.*armHHbar[20] + 445./
   12.*armHHbar[21] + 77./24.*armHHbar[44] - 37./8.*armHHbar[70] + 17./
   12.*armHHbar[93] - 525./16.*armHHbar[72] - 8./3.*armHHbar[96] + 141./
   4.*armHHbar[95];
   armHHbar[114]=armHHbar[68]*armHHbar[114];
   armHHbar[119]=1./2.*armHHbar[25];
   armHHbar[143]=armHHbar[119] + 3./2. + 5*armHHbar[23];
   armHHbar[175]=armHHbar[147] - 15./2.*armHHbar[16] + 3./2.*
   armHHbar[143] + 5*armHHbar[14];
   armHHbar[175]=armHHbar[31]*armHHbar[175];
   armHHbar[143]=armHHbar[180] - 5./2.*armHHbar[16] + 1./2.*
   armHHbar[143] + 5./3.*armHHbar[14];
   armHHbar[143]=armHHbar[1]*armHHbar[143];
   armHHbar[147]= - 38./3.*armHHbar[9] - 34./9. + armHHbar[147];
   armHHbar[147]=armHHbar[31]*armHHbar[147];
   armHHbar[181]= - 6*armHHbar[9] - 2./3. + armHHbar[180];
   armHHbar[181]=armHHbar[1]*armHHbar[181];
   armHHbar[147]=armHHbar[147] + armHHbar[181];
   armHHbar[147]=armHHbar[12]*armHHbar[147];
   armHHbar[143]=armHHbar[186] + armHHbar[147] + armHHbar[175] + 
   armHHbar[143];
   armHHbar[143]=MMZ*armHHbar[143];
   armHHbar[147]= - armHHbar[14] - 1 - armHHbar[25];
   armHHbar[175]=armHHbar[31]*armHHbar[147];
   armHHbar[147]=armHHbar[1]*armHHbar[147];
   armHHbar[147]=3*armHHbar[175] + armHHbar[147];
   armHHbar[147]=MMZ*armHHbar[147];
   armHHbar[175]= - armHHbar[31]*armHHbar[6];
   armHHbar[181]= - armHHbar[1]*armHHbar[6];
   armHHbar[182]=3*armHHbar[31] + armHHbar[1];
   armHHbar[182]=armHHbar[19]*armHHbar[182];
   armHHbar[147]=2*armHHbar[147] + armHHbar[182] + 3*armHHbar[175] + 
   armHHbar[181];
   armHHbar[147]=MMZ*armHHbar[147];
   armHHbar[175]= - 33./2.*armHHbar[90] - 231./4.*armHHbar[91] - 2*
   armHHbar[89] - 289./4. + 4*armHHbar[92];
   armHHbar[175]=MMZ*armHHbar[175];
   armHHbar[181]= - armHHbar[12]*armHHbar[17];
   armHHbar[152]=armHHbar[175] + armHHbar[152] + 37*armHHbar[18] + 62*
   armHHbar[19] + armHHbar[181] + armHHbar[70] - 2*armHHbar[93] + 33./2.
   *armHHbar[72] - 33*armHHbar[95] + 2*armHHbar[94] - armHHbar[71];
   armHHbar[152]=armHHbar[68]*MMZ*armHHbar[152];
   armHHbar[147]=armHHbar[147] + 3*armHHbar[152];
   armHHbar[147]=armHHbar[3]*armHHbar[147];
   armHHbar[152]=1./3. - armHHbar[9];
   armHHbar[175]=armHHbar[31]*armHHbar[152];
   armHHbar[152]=armHHbar[1]*armHHbar[152];
   armHHbar[152]=5./3.*armHHbar[175] + armHHbar[152];
   armHHbar[152]=armHHbar[19]*armHHbar[152];
   armHHbar[114]=armHHbar[147] + armHHbar[114] + armHHbar[143] + 4*
   armHHbar[152] + armHHbar[185];
   armHHbar[114]=armHHbar[3]*armHHbar[114];
   armHHbar[143]=armHHbar[98]*armHHbar[10];
   armHHbar[147]=armHHbar[31]*armHHbar[143];
   armHHbar[143]=armHHbar[1]*armHHbar[143];
   armHHbar[143]=3*armHHbar[147] + armHHbar[143];
   armHHbar[143]=armHHbar[12]*armHHbar[143];
   armHHbar[119]= - armHHbar[10] + armHHbar[16] + 3./2.*armHHbar[14] + 
   armHHbar[119] + 3./2. - armHHbar[23];
   armHHbar[119]=armHHbar[98]*armHHbar[119];
   armHHbar[147]=armHHbar[31]*armHHbar[119];
   armHHbar[119]=armHHbar[1]*armHHbar[119];
   armHHbar[119]=armHHbar[143] + 3*armHHbar[147] + armHHbar[119];
   armHHbar[119]=MMZ*armHHbar[119];
   armHHbar[143]= - 1./3. + armHHbar[9];
   armHHbar[147]=armHHbar[31]*armHHbar[143];
   armHHbar[143]=armHHbar[1]*armHHbar[143];
   armHHbar[143]=5./3.*armHHbar[147] + armHHbar[143];
   armHHbar[143]=armHHbar[12]*armHHbar[143];
   armHHbar[100]=armHHbar[110] + armHHbar[103] + armHHbar[114] + 
   armHHbar[100] + 1./2.*armHHbar[119] + 4./3.*armHHbar[143] + 
   armHHbar[187];
   armHHbar[100]=armHHbar[256]*armHHbar[100];
   armHHbar[103]=13./12.*armHHbar[39];
   armHHbar[110]= - 7./12.*armHHbar[38];
   armHHbar[114]=7./6.*armHHbar[36];
   armHHbar[119]= - 21./4.*armHHbar[7];
   armHHbar[147]= - 7./4.*armHHbar[12];
   armHHbar[152]= - armHHbar[11] + armHHbar[147] + armHHbar[119] + 
   armHHbar[114] + armHHbar[110] + armHHbar[103] + 3 - 1./48.*
   armHHbar[40];
   armHHbar[152]=armHHbar[11]*armHHbar[152];
   armHHbar[175]=armHHbar[145] + 19 + armHHbar[123];
   armHHbar[175]=1./2.*armHHbar[175] + armHHbar[38];
   armHHbar[182]= - 1./2.*armHHbar[36];
   armHHbar[175]=armHHbar[136] + armHHbar[119] + 1./3.*armHHbar[175] + 
   armHHbar[182];
   armHHbar[175]=armHHbar[12]*armHHbar[175];
   armHHbar[185]= - 15./4.*armHHbar[88] + 21./4.*armHHbar[83] + 3./4.*
   armHHbar[86] - 1 + 21./4.*armHHbar[84];
   armHHbar[186]= - 3./8.*armHHbar[36] - 3./8.*armHHbar[38] + 13./8.*
   armHHbar[39] + 2 - 1./32.*armHHbar[40];
   armHHbar[186]=armHHbar[7]*armHHbar[186];
   armHHbar[187]=2*armHHbar[85];
   armHHbar[195]= - 21./4.*armHHbar[42];
   armHHbar[274]=3./8.*armHHbar[90];
   armHHbar[284]= - 1./8.*armHHbar[16];
   armHHbar[285]=armHHbar[76] - 1./2.*armHHbar[79] + armHHbar[142] + 
   armHHbar[80] + armHHbar[78];
   armHHbar[285]=1./2.*MMH*armHHbar[285];
   armHHbar[286]=2*armHHbar[87];
   armHHbar[287]= - 1./12.*armHHbar[38];
   armHHbar[288]= - 1./12.*armHHbar[36];
   armHHbar[105]=armHHbar[285] + 1./2.*armHHbar[152] + 1./2.*
   armHHbar[175] + 3*armHHbar[186] + armHHbar[288] + armHHbar[287] + 
   armHHbar[105] + armHHbar[284] + armHHbar[274] + armHHbar[195] - 15./
   8.*armHHbar[91] + armHHbar[187] + 1./2.*armHHbar[185] + 
   armHHbar[286];
   armHHbar[105]=MMH*armHHbar[105];
   armHHbar[152]= - 521 - armHHbar[40];
   armHHbar[152]=1./4.*armHHbar[152] + armHHbar[145];
   armHHbar[175]=1./8.*armHHbar[248];
   armHHbar[185]= - 1./8.*armHHbar[19]*armHHbar[98];
   armHHbar[108]=armHHbar[185] + armHHbar[108] + armHHbar[175] + 
   armHHbar[182] + 1./6.*armHHbar[152] + armHHbar[264];
   armHHbar[108]=armHHbar[19]*armHHbar[108];
   armHHbar[123]=armHHbar[109] + armHHbar[154] + armHHbar[145] + 67./3.
    + armHHbar[123];
   armHHbar[123]=armHHbar[17]*armHHbar[123];
   armHHbar[152]= - 1./2.*armHHbar[40];
   armHHbar[186]= - 19 + armHHbar[152];
   armHHbar[186]=1./2.*armHHbar[186] + armHHbar[145];
   armHHbar[248]= - 9*armHHbar[36];
   armHHbar[289]=1./2.*armHHbar[102];
   armHHbar[186]=armHHbar[289] + armHHbar[248] + 1./3.*armHHbar[186] + 
   armHHbar[38];
   armHHbar[290]=1./4.*armHHbar[104];
   armHHbar[186]=armHHbar[290] + 1./2.*armHHbar[186] - 7*armHHbar[12];
   armHHbar[186]=armHHbar[18]*armHHbar[186];
   armHHbar[291]=19./3.*armHHbar[19];
   armHHbar[292]= - 3./2.*armHHbar[17];
   armHHbar[293]= - 53./6.*armHHbar[18] + armHHbar[292] + armHHbar[291]
   ;
   armHHbar[293]=1./4.*armHHbar[11]*armHHbar[293];
   armHHbar[294]=17./8.*armHHbar[71];
   armHHbar[295]=armHHbar[294] + 4./3.*armHHbar[96] - armHHbar[94];
   armHHbar[296]=17./8.*armHHbar[70];
   armHHbar[297]= - 9./4.*armHHbar[44];
   armHHbar[298]=1./2.*armHHbar[181];
   armHHbar[105]=armHHbar[105] + armHHbar[293] + 1./2.*armHHbar[186] + 
   armHHbar[108] + armHHbar[298] + 1./4.*armHHbar[123] - 241./24.*
   armHHbar[20] + 73./24.*armHHbar[21] + armHHbar[297] + armHHbar[296]
    - armHHbar[93] - 3./2.*armHHbar[72] + armHHbar[295] + 6*
   armHHbar[95];
   armHHbar[105]=armHHbar[68]*armHHbar[105];
   armHHbar[108]=101 + armHHbar[144];
   armHHbar[108]=1./8.*armHHbar[104] - 1./3.*armHHbar[12] + 1./8.*
   armHHbar[102] - 5./2.*armHHbar[36] + 1./12.*armHHbar[108] + 
   armHHbar[38];
   armHHbar[108]=armHHbar[18]*armHHbar[108];
   armHHbar[123]= - 61./3.*armHHbar[95] + 13*armHHbar[72];
   armHHbar[144]= - 17*armHHbar[39];
   armHHbar[186]=7./3. + armHHbar[144];
   armHHbar[186]=armHHbar[17]*armHHbar[186];
   armHHbar[299]=139 + armHHbar[144];
   armHHbar[299]=1./4.*armHHbar[299] + armHHbar[12];
   armHHbar[299]=armHHbar[19]*armHHbar[299];
   armHHbar[291]= - armHHbar[17] + armHHbar[291];
   armHHbar[291]=1./2.*armHHbar[291] - armHHbar[18];
   armHHbar[291]=armHHbar[11]*armHHbar[291];
   armHHbar[108]=armHHbar[291] + armHHbar[108] + 1./3.*armHHbar[299] + 
   1./8.*armHHbar[186] - 17./3.*armHHbar[20] - 11./4.*armHHbar[21] - 9./
   8.*armHHbar[44] + armHHbar[296] + 1./4.*armHHbar[123] - armHHbar[93]
   ;
   armHHbar[123]= - 1 - 17./4.*armHHbar[39];
   armHHbar[186]=1./3.*armHHbar[123] + armHHbar[126];
   armHHbar[186]=armHHbar[12]*armHHbar[186];
   armHHbar[291]=1./4. + 3*armHHbar[83];
   armHHbar[291]=7./2.*armHHbar[291] + 15*armHHbar[88];
   armHHbar[123]=armHHbar[7]*armHHbar[123];
   armHHbar[299]= - 7./4.*armHHbar[38] + 1 - 17./16.*armHHbar[39];
   armHHbar[299]= - 1./2.*armHHbar[11] - 1./2.*armHHbar[12] - 21./8.*
   armHHbar[7] + 1./3.*armHHbar[299] + armHHbar[36];
   armHHbar[299]=armHHbar[11]*armHHbar[299];
   armHHbar[300]=1./2.*armHHbar[79] + armHHbar[76];
   armHHbar[300]=MMH*armHHbar[300];
   armHHbar[301]= - 1./24.*armHHbar[36];
   armHHbar[123]=1./4.*armHHbar[300] + 1./2.*armHHbar[299] + 1./4.*
   armHHbar[186] + 3./8.*armHHbar[123] + armHHbar[301] - armHHbar[15]
    - 15./8.*armHHbar[16] + 3./2.*armHHbar[90] - 21./16.*armHHbar[42]
    - 1./8.*armHHbar[91] + 1./8.*armHHbar[291] + armHHbar[85];
   armHHbar[123]=MMH*armHHbar[123];
   armHHbar[108]=1./2.*armHHbar[108] + armHHbar[123];
   armHHbar[108]=armHHbar[68]*armHHbar[108];
   armHHbar[123]= - armHHbar[7]*armHHbar[39];
   armHHbar[186]= - 1./12.*armHHbar[39];
   armHHbar[291]=armHHbar[186] - armHHbar[38];
   armHHbar[291]=armHHbar[11]*armHHbar[291];
   armHHbar[123]=armHHbar[291] + 1./6.*armHHbar[128] + 3./4.*
   armHHbar[123] - armHHbar[16] + armHHbar[90] + 1./4. + armHHbar[88];
   armHHbar[123]=MMH*armHHbar[123];
   armHHbar[128]= - 1./2.*armHHbar[17];
   armHHbar[291]=1./3.*armHHbar[19];
   armHHbar[299]=1./6.*armHHbar[18] + armHHbar[128] + armHHbar[291];
   armHHbar[299]=armHHbar[11]*armHHbar[299];
   armHHbar[300]= - armHHbar[17]*armHHbar[39];
   armHHbar[302]=5 - 1./3.*armHHbar[39];
   armHHbar[302]=armHHbar[19]*armHHbar[302];
   armHHbar[186]=1 + armHHbar[186];
   armHHbar[186]=armHHbar[18]*armHHbar[186];
   armHHbar[115]=1./2.*armHHbar[123] + armHHbar[299] + armHHbar[186] + 
   1./2.*armHHbar[302] + 1./4.*armHHbar[300] + armHHbar[115] - 1./2.*
   armHHbar[21];
   armHHbar[115]=armHHbar[68]*armHHbar[115];
   armHHbar[123]=armHHbar[3]*armHHbar[68]*armHHbar[18]*armHHbar[173];
   armHHbar[115]=armHHbar[115] + 1./2.*armHHbar[123];
   armHHbar[115]=armHHbar[171]*armHHbar[115];
   armHHbar[123]= - armHHbar[31]*armHHbar[13];
   armHHbar[173]= - armHHbar[1]*armHHbar[13];
   armHHbar[186]=11./9.*armHHbar[31] + armHHbar[1];
   armHHbar[299]=armHHbar[11]*armHHbar[186];
   armHHbar[123]=armHHbar[299] + 11./9.*armHHbar[123] + armHHbar[173];
   armHHbar[123]=MMH*armHHbar[123];
   armHHbar[173]=armHHbar[18]*armHHbar[186];
   armHHbar[123]=armHHbar[173] + armHHbar[123];
   armHHbar[173]=1./2.*armHHbar[17];
   armHHbar[186]= - 43./2.*armHHbar[18] + armHHbar[173] - 91*
   armHHbar[19];
   armHHbar[186]=armHHbar[18]*armHHbar[186];
   armHHbar[299]=pow(armHHbar[19],2);
   armHHbar[186]= - armHHbar[299] + 1./8.*armHHbar[186];
   armHHbar[186]=armHHbar[3]*armHHbar[68]*armHHbar[186];
   armHHbar[108]=1./4.*armHHbar[115] + 1./3.*armHHbar[186] + 1./4.*
   armHHbar[123] + armHHbar[108];
   armHHbar[108]=armHHbar[171]*armHHbar[108];
   armHHbar[115]=5./3. - 9./4.*armHHbar[10];
   armHHbar[115]=armHHbar[7]*armHHbar[115];
   armHHbar[123]= - armHHbar[9]*armHHbar[7];
   armHHbar[115]=11./4.*armHHbar[123] - armHHbar[14] + armHHbar[115];
   armHHbar[115]=armHHbar[31]*armHHbar[115];
   armHHbar[186]= - 3./4.*armHHbar[10];
   armHHbar[300]=1 + armHHbar[186];
   armHHbar[300]=armHHbar[7]*armHHbar[300];
   armHHbar[300]=9./4.*armHHbar[123] - 1./3.*armHHbar[14] + 
   armHHbar[300];
   armHHbar[300]=armHHbar[1]*armHHbar[300];
   armHHbar[302]=5./3. + armHHbar[180];
   armHHbar[224]=1./3.*armHHbar[302] + armHHbar[224];
   armHHbar[224]=armHHbar[1]*armHHbar[224];
   armHHbar[180]= - 11./18.*armHHbar[9] + 37./27. + armHHbar[180];
   armHHbar[180]=armHHbar[31]*armHHbar[180];
   armHHbar[180]=armHHbar[180] + armHHbar[224];
   armHHbar[180]=armHHbar[12]*armHHbar[180];
   armHHbar[224]= - 1./4.*armHHbar[10];
   armHHbar[302]=1./3. + armHHbar[224];
   armHHbar[172]=1./3.*armHHbar[302] + armHHbar[172];
   armHHbar[172]=armHHbar[1]*armHHbar[172];
   armHHbar[224]= - 11./36.*armHHbar[9] + 5./27. + armHHbar[224];
   armHHbar[224]=armHHbar[31]*armHHbar[224];
   armHHbar[172]=armHHbar[224] + armHHbar[172];
   armHHbar[224]=armHHbar[11]*armHHbar[172];
   armHHbar[115]=armHHbar[224] + armHHbar[180] + armHHbar[115] + 
   armHHbar[300];
   armHHbar[115]=MMH*armHHbar[115];
   armHHbar[180]=armHHbar[238] + 5./9. + armHHbar[186];
   armHHbar[180]=armHHbar[31]*armHHbar[17]*armHHbar[180];
   armHHbar[176]=armHHbar[302] + armHHbar[176];
   armHHbar[176]=armHHbar[1]*armHHbar[17]*armHHbar[176];
   armHHbar[186]= - 11./9.*armHHbar[9] + 47./27. - armHHbar[10];
   armHHbar[186]=armHHbar[31]*armHHbar[186];
   armHHbar[224]=7./3. - armHHbar[10];
   armHHbar[224]=1./3.*armHHbar[224] - armHHbar[9];
   armHHbar[224]=armHHbar[1]*armHHbar[224];
   armHHbar[186]=armHHbar[186] + armHHbar[224];
   armHHbar[186]=armHHbar[19]*armHHbar[186];
   armHHbar[172]=armHHbar[18]*armHHbar[172];
   armHHbar[115]=1./2.*armHHbar[115] + armHHbar[172] + 1./2.*
   armHHbar[186] + armHHbar[180] + armHHbar[176];
   armHHbar[172]=armHHbar[17] - 55*armHHbar[19];
   armHHbar[172]=armHHbar[19]*armHHbar[172];
   armHHbar[176]=53./24.*armHHbar[18] + 1./24.*armHHbar[17] + 6*
   armHHbar[19];
   armHHbar[176]=armHHbar[18]*armHHbar[176];
   armHHbar[172]=1./24.*armHHbar[172] + armHHbar[176];
   armHHbar[172]=armHHbar[3]*armHHbar[68]*armHHbar[172];
   armHHbar[105]=armHHbar[108] + armHHbar[172] + armHHbar[115] + 
   armHHbar[105];
   armHHbar[105]=armHHbar[171]*armHHbar[105];
   armHHbar[103]= - armHHbar[11] + armHHbar[147] + armHHbar[119] + 
   armHHbar[114] + armHHbar[110] + armHHbar[103] + 7 + 47./48.*
   armHHbar[40];
   armHHbar[103]=armHHbar[11]*armHHbar[103];
   armHHbar[108]=47./4.*armHHbar[40];
   armHHbar[110]=armHHbar[145] + 31 + armHHbar[108];
   armHHbar[110]=1./2.*armHHbar[110] + armHHbar[38];
   armHHbar[110]=armHHbar[136] + armHHbar[119] + 1./3.*armHHbar[110] + 
   armHHbar[182];
   armHHbar[110]=armHHbar[12]*armHHbar[110];
   armHHbar[114]=1 + 3./2.*armHHbar[84];
   armHHbar[172]=21./2.*armHHbar[83];
   armHHbar[114]= - 15./2.*armHHbar[88] + armHHbar[172] + 7*
   armHHbar[114] + 15./2.*armHHbar[86];
   armHHbar[180]= - 3./4.*armHHbar[36] - 3./4.*armHHbar[38] + 13./4.*
   armHHbar[39] + 7 + 47./16.*armHHbar[40];
   armHHbar[180]=armHHbar[7]*armHHbar[180];
   armHHbar[103]=armHHbar[285] + 1./2.*armHHbar[103] + 1./2.*
   armHHbar[110] + 3./2.*armHHbar[180] + armHHbar[288] + armHHbar[287]
    - 31./8.*armHHbar[15] + armHHbar[284] + armHHbar[274] + 
   armHHbar[195] - 3./8.*armHHbar[91] + armHHbar[187] + 1./4.*
   armHHbar[114] + armHHbar[286];
   armHHbar[103]=MMH*armHHbar[103];
   armHHbar[110]= - 113 + 47*armHHbar[40];
   armHHbar[110]=1./4.*armHHbar[110] + armHHbar[145];
   armHHbar[110]=armHHbar[185] - armHHbar[12] + armHHbar[175] + 
   armHHbar[182] + 1./6.*armHHbar[110] + armHHbar[264];
   armHHbar[110]=armHHbar[19]*armHHbar[110];
   armHHbar[108]=armHHbar[109] + armHHbar[154] + armHHbar[145] + 103./3.
    + armHHbar[108];
   armHHbar[108]=armHHbar[17]*armHHbar[108];
   armHHbar[114]=47./2.*armHHbar[40];
   armHHbar[180]=353 + armHHbar[114];
   armHHbar[180]=1./2.*armHHbar[180] + armHHbar[145];
   armHHbar[180]=armHHbar[289] + armHHbar[248] + 1./3.*armHHbar[180] + 
   armHHbar[38];
   armHHbar[180]=armHHbar[290] + 1./2.*armHHbar[180] - 15*armHHbar[12];
   armHHbar[180]=armHHbar[18]*armHHbar[180];
   armHHbar[103]=armHHbar[103] + armHHbar[293] + 1./2.*armHHbar[180] + 
   armHHbar[110] + armHHbar[298] + 1./4.*armHHbar[108] - 373./24.*
   armHHbar[20] + 145./24.*armHHbar[21] + armHHbar[297] + armHHbar[296]
    + armHHbar[295] - armHHbar[93];
   armHHbar[103]=armHHbar[68]*armHHbar[103];
   armHHbar[108]= - 5*armHHbar[37] + 51 - 5*armHHbar[41];
   armHHbar[108]= - 55*armHHbar[39] + 1./2.*armHHbar[108] + 55*
   armHHbar[40];
   armHHbar[108]=2*armHHbar[191] + armHHbar[260] + 1./2.*armHHbar[108]
    + 2*armHHbar[38];
   armHHbar[110]=2*armHHbar[193];
   armHHbar[180]= - 3*armHHbar[12];
   armHHbar[108]= - 3./2.*armHHbar[27] + armHHbar[180] + 3*
   armHHbar[108] + armHHbar[110];
   armHHbar[108]=armHHbar[27]*armHHbar[108];
   armHHbar[186]=2*armHHbar[61];
   armHHbar[187]= - armHHbar[34] + armHHbar[186];
   armHHbar[195]= - 2*armHHbar[19];
   armHHbar[168]= - armHHbar[18] - 4*armHHbar[28] + armHHbar[195] + 
   armHHbar[168] + armHHbar[187] + armHHbar[67];
   armHHbar[168]=armHHbar[3]*armHHbar[168];
   armHHbar[224]= - 11 + 27./2.*armHHbar[59];
   armHHbar[238]=3./4.*armHHbar[56];
   armHHbar[248]= - 3./4.*armHHbar[41];
   armHHbar[274]= - 13./8.*armHHbar[64];
   armHHbar[284]= - 3./4.*armHHbar[16];
   armHHbar[285]= - armHHbar[12]*armHHbar[30];
   armHHbar[289]=9./4.*armHHbar[285];
   armHHbar[290]= - armHHbar[11]*armHHbar[30];
   armHHbar[293]=3./4.*armHHbar[290];
   armHHbar[295]=3*armHHbar[62];
   armHHbar[108]=3*armHHbar[168] + armHHbar[293] + armHHbar[108] + 
   armHHbar[289] + armHHbar[234] - 13./4.*armHHbar[63] + armHHbar[284]
    + armHHbar[216] + armHHbar[274] - 33./4.*armHHbar[43] - 3./4.*
   armHHbar[37] + armHHbar[248] + 3./2.*armHHbar[58] + armHHbar[238] + 
   1./2.*armHHbar[224] + armHHbar[295];
   armHHbar[108]=armHHbar[3]*armHHbar[108];
   armHHbar[168]=1 + armHHbar[295];
   armHHbar[168]= - armHHbar[30] + armHHbar[16] + armHHbar[231] + 1./2.
   *armHHbar[168] - armHHbar[56];
   armHHbar[168]=armHHbar[98]*armHHbar[168];
   armHHbar[224]=3*armHHbar[51] + 1./4.*armHHbar[53];
   armHHbar[296]=1 + armHHbar[63];
   armHHbar[296]=armHHbar[99]*armHHbar[296];
   armHHbar[224]=armHHbar[237] + 1./8.*armHHbar[296] + 1./4.*
   armHHbar[168] + 1./2.*armHHbar[224] + armHHbar[50];
   armHHbar[237]=5 + armHHbar[62];
   armHHbar[237]=armHHbar[63] + 1./2.*armHHbar[237] + armHHbar[64];
   armHHbar[237]=armHHbar[3]*armHHbar[237];
   armHHbar[237]=armHHbar[237] - armHHbar[51] - 1./2.*armHHbar[53];
   armHHbar[237]=MMt*armHHbar[3]*armHHbar[237];
   armHHbar[108]=3*armHHbar[237] + 3*armHHbar[224] + armHHbar[108];
   armHHbar[108]=MMt*armHHbar[108];
   armHHbar[224]= - 7*armHHbar[62] - 25./4. + 3*armHHbar[59];
   armHHbar[224]= - 15./4.*armHHbar[43] + 9./4.*armHHbar[58] + 1./2.*
   armHHbar[224] + armHHbar[56];
   armHHbar[237]=armHHbar[272] + armHHbar[20];
   armHHbar[237]=armHHbar[99]*armHHbar[237];
   armHHbar[272]=3./4.*armHHbar[15];
   armHHbar[297]=1./8.*armHHbar[63];
   armHHbar[298]= - 3./2.*armHHbar[16];
   armHHbar[300]= - 1./2.*armHHbar[34];
   armHHbar[302]=armHHbar[300] + armHHbar[61];
   armHHbar[303]=armHHbar[302] + armHHbar[21];
   armHHbar[303]=3./2.*armHHbar[98]*armHHbar[303];
   armHHbar[216]=3./4.*armHHbar[237] + armHHbar[303] - armHHbar[30] + 
   armHHbar[211] + armHHbar[272] + armHHbar[297] + armHHbar[298] + 
   armHHbar[216] - 3./4.*armHHbar[55] + 3./2.*armHHbar[224] + 
   armHHbar[64];
   armHHbar[224]= - 33*armHHbar[39] + 33*armHHbar[40] - 3./2.*
   armHHbar[37] + 19 - 3./2.*armHHbar[41];
   armHHbar[304]=6*armHHbar[38];
   armHHbar[305]= - 9./2.*armHHbar[12];
   armHHbar[110]=armHHbar[305] + armHHbar[110] + 6*armHHbar[191] + 
   armHHbar[177] + 5./2.*armHHbar[224] + armHHbar[304];
   armHHbar[110]=armHHbar[28]*armHHbar[110];
   armHHbar[224]= - 11./4.*armHHbar[61];
   armHHbar[306]=7./2.*armHHbar[21];
   armHHbar[307]= - 9./4.*armHHbar[30];
   armHHbar[308]=armHHbar[307] - 7./2. + 8*armHHbar[29];
   armHHbar[308]=armHHbar[19]*armHHbar[308];
   armHHbar[309]=4*armHHbar[19] - 5./4.*armHHbar[28];
   armHHbar[309]=armHHbar[27]*armHHbar[309];
   armHHbar[310]=39./4.*armHHbar[27] + 1 + 23./4.*armHHbar[29];
   armHHbar[310]=armHHbar[18]*armHHbar[310];
   armHHbar[311]=3./4.*armHHbar[251];
   armHHbar[110]=armHHbar[311] + 1./2.*armHHbar[310] + 3*armHHbar[309]
    + armHHbar[110] + armHHbar[308] + 5./8.*armHHbar[20] + 
   armHHbar[306] - 41./4.*armHHbar[45] + 57./16.*armHHbar[33] - 5./2.*
   armHHbar[67] + 6*armHHbar[34] + armHHbar[224];
   armHHbar[110]=armHHbar[3]*armHHbar[110];
   armHHbar[308]=55*armHHbar[39] - 55*armHHbar[40] + armHHbar[37] - 51./
   2. + armHHbar[41];
   armHHbar[308]=armHHbar[267] - armHHbar[36] + 1./4.*armHHbar[308] - 
   armHHbar[38];
   armHHbar[308]=armHHbar[180] + 3*armHHbar[308] + armHHbar[268];
   armHHbar[308]=armHHbar[27]*armHHbar[308];
   armHHbar[194]= - 7 + armHHbar[194];
   armHHbar[309]=1./4.*armHHbar[30];
   armHHbar[194]=2./3.*armHHbar[194] + armHHbar[309];
   armHHbar[194]=armHHbar[12]*armHHbar[194];
   armHHbar[310]= - 19 - 29./2.*armHHbar[29];
   armHHbar[312]= - 3*armHHbar[27];
   armHHbar[310]=1./3.*armHHbar[310] + armHHbar[312];
   armHHbar[310]=armHHbar[11]*armHHbar[310];
   armHHbar[313]= - 1 + armHHbar[242];
   armHHbar[314]=armHHbar[19]*armHHbar[98]*armHHbar[313];
   armHHbar[315]=3./2.*armHHbar[314];
   armHHbar[316]= - armHHbar[98] - armHHbar[99];
   armHHbar[316]=armHHbar[28]*armHHbar[316];
   armHHbar[221]= - 1 + armHHbar[221];
   armHHbar[221]=armHHbar[18]*armHHbar[99]*armHHbar[221];
   armHHbar[108]=armHHbar[108] + armHHbar[110] + 1./4.*armHHbar[310] + 
   3./4.*armHHbar[221] + 1./2.*armHHbar[308] + 3./4.*armHHbar[316] + 
   armHHbar[315] + 1./2.*armHHbar[216] + armHHbar[194];
   armHHbar[108]=MMt*armHHbar[108];
   armHHbar[110]=armHHbar[97]*armHHbar[45];
   armHHbar[194]=armHHbar[9]*armHHbar[7];
   armHHbar[216]=armHHbar[7]*armHHbar[29];
   armHHbar[308]= - armHHbar[30]*armHHbar[7];
   armHHbar[220]=3./2.*armHHbar[194] + 9*armHHbar[308] + 15./2.*
   armHHbar[216] + 3./2.*armHHbar[110] - 3./8.*armHHbar[15] + 3./8.*
   armHHbar[63] + armHHbar[220] + 47./8. - armHHbar[13];
   armHHbar[310]=5./3.*armHHbar[29];
   armHHbar[316]=13./4. + armHHbar[310];
   armHHbar[316]=3./2.*armHHbar[196] + armHHbar[156] + 1./2.*
   armHHbar[316] - armHHbar[30];
   armHHbar[316]=armHHbar[11]*armHHbar[316];
   armHHbar[317]=1./12.*armHHbar[9] + armHHbar[207] + 1 + 5./12.*
   armHHbar[29];
   armHHbar[317]=armHHbar[12]*armHHbar[317];
   armHHbar[220]=1./4.*armHHbar[316] + 3./8.*armHHbar[206] + 1./4.*
   armHHbar[220] + armHHbar[317];
   armHHbar[220]=MMH*armHHbar[220];
   armHHbar[316]=91./8. + armHHbar[310];
   armHHbar[316]= - 39./8.*armHHbar[27] + 3./4.*armHHbar[196] + 
   armHHbar[156] + 1./2.*armHHbar[316] - armHHbar[30];
   armHHbar[316]=armHHbar[18]*armHHbar[316];
   armHHbar[317]=armHHbar[188] + 5./2.*armHHbar[29] + armHHbar[197];
   armHHbar[318]=armHHbar[27]*armHHbar[317];
   armHHbar[319]=6*armHHbar[30];
   armHHbar[320]= - armHHbar[9] - 5*armHHbar[29] + armHHbar[319];
   armHHbar[321]=armHHbar[3]*armHHbar[28]*armHHbar[320];
   armHHbar[320]=MMt*armHHbar[3]*armHHbar[27]*armHHbar[320];
   armHHbar[318]=armHHbar[320] + 1./2.*armHHbar[318] + armHHbar[321];
   armHHbar[318]=armHHbar[26]*MMt*armHHbar[318];
   armHHbar[317]=armHHbar[17]*armHHbar[317];
   armHHbar[310]=17 + armHHbar[310];
   armHHbar[156]=armHHbar[156] + 1./2.*armHHbar[310] - armHHbar[30];
   armHHbar[156]=armHHbar[19]*armHHbar[156];
   armHHbar[310]= - 5./3.*armHHbar[12];
   armHHbar[320]=3./16.*armHHbar[206] + 27./32. + armHHbar[310];
   armHHbar[320]=armHHbar[28]*armHHbar[320];
   armHHbar[321]=39*armHHbar[19] - 7./2.*armHHbar[28];
   armHHbar[321]=armHHbar[28]*armHHbar[321];
   armHHbar[322]=armHHbar[18]*armHHbar[28];
   armHHbar[321]=armHHbar[321] + 39./2.*armHHbar[322];
   armHHbar[321]=armHHbar[3]*armHHbar[321];
   armHHbar[323]= - 9./4.*armHHbar[21];
   armHHbar[324]= - armHHbar[19] + 5./16.*armHHbar[28];
   armHHbar[324]=armHHbar[27]*armHHbar[324];
   armHHbar[325]= - 1./4.*armHHbar[67];
   armHHbar[108]=armHHbar[318] + armHHbar[108] + 1./4.*armHHbar[321] + 
   1./2.*armHHbar[220] + 19./12.*armHHbar[251] + 1./4.*armHHbar[316] + 
   3*armHHbar[324] + armHHbar[320] + 1./2.*armHHbar[156] + 1./4.*
   armHHbar[317] - 57./64.*armHHbar[20] + armHHbar[323] + 63./32.*
   armHHbar[45] + 3./64.*armHHbar[33] - armHHbar[61] + armHHbar[325];
   armHHbar[108]=armHHbar[26]*armHHbar[108];
   armHHbar[156]= - 55./8.*armHHbar[40];
   armHHbar[220]=armHHbar[111] - 9 + armHHbar[156];
   armHHbar[316]=1./3.*armHHbar[38];
   armHHbar[317]= - 13./8.*armHHbar[12] + armHHbar[119] + armHHbar[220]
    + armHHbar[316];
   armHHbar[317]=armHHbar[12]*armHHbar[317];
   armHHbar[318]=1./3.*armHHbar[36];
   armHHbar[119]= - armHHbar[11] - 5./2.*armHHbar[12] + armHHbar[119]
    + armHHbar[220] + armHHbar[318];
   armHHbar[119]=armHHbar[11]*armHHbar[119];
   armHHbar[121]=armHHbar[255] + armHHbar[121] + armHHbar[78] + 
   armHHbar[142];
   armHHbar[121]=MMH*armHHbar[121];
   armHHbar[142]=armHHbar[172] - 13*armHHbar[86] - 467./8. + 21*
   armHHbar[84];
   armHHbar[142]=1./2.*armHHbar[142] - 13*armHHbar[88];
   armHHbar[172]=5./8.*armHHbar[15];
   armHHbar[220]= - armHHbar[36] - armHHbar[38] + 55./4.*armHHbar[39]
    - 9 - 55./4.*armHHbar[40];
   armHHbar[220]=armHHbar[7]*armHHbar[220];
   armHHbar[119]=1./2.*armHHbar[121] + 1./4.*armHHbar[119] + 1./2.*
   armHHbar[317] + 9./8.*armHHbar[220] + armHHbar[301] + armHHbar[287]
    + armHHbar[172] + 5./4.*armHHbar[16] - 17./8.*armHHbar[90] - 63./16.
   *armHHbar[42] - 17./4.*armHHbar[91] + armHHbar[85] + 1./4.*
   armHHbar[142] + armHHbar[286];
   armHHbar[119]=MMH*armHHbar[119];
   armHHbar[111]=armHHbar[111] - 73 + armHHbar[156];
   armHHbar[121]=armHHbar[185] + 49./6.*armHHbar[12] + armHHbar[175] + 
   armHHbar[182] + armHHbar[111] + armHHbar[264];
   armHHbar[121]=armHHbar[19]*armHHbar[121];
   armHHbar[142]=165*armHHbar[39] - 89 - 165*armHHbar[40];
   armHHbar[142]=armHHbar[109] + 1./4.*armHHbar[142] + armHHbar[154];
   armHHbar[142]=armHHbar[17]*armHHbar[142];
   armHHbar[111]=armHHbar[111] + armHHbar[139];
   armHHbar[102]=1./16.*armHHbar[104] + 53./12.*armHHbar[12] + 1./16.*
   armHHbar[102] + 1./2.*armHHbar[111] - armHHbar[36];
   armHHbar[102]=armHHbar[18]*armHHbar[102];
   armHHbar[104]= - 41./24.*armHHbar[18] - 1./8.*armHHbar[17] + 8*
   armHHbar[19];
   armHHbar[104]=armHHbar[11]*armHHbar[104];
   armHHbar[111]= - 1./2.*armHHbar[93];
   armHHbar[102]=armHHbar[119] + armHHbar[104] + armHHbar[102] + 
   armHHbar[121] + 1./4.*armHHbar[181] + 1./4.*armHHbar[142] + 79./24.*
   armHHbar[20] + 79./12.*armHHbar[21] - 27./16.*armHHbar[44] + 17./16.
   *armHHbar[70] + armHHbar[111] - 45./8.*armHHbar[72] + 171./8.*
   armHHbar[95] - armHHbar[94] + armHHbar[294];
   armHHbar[102]=armHHbar[68]*armHHbar[102];
   armHHbar[104]= - armHHbar[7]*armHHbar[29];
   armHHbar[119]=armHHbar[30]*armHHbar[7];
   armHHbar[104]=1./2.*armHHbar[123] + 1./2.*armHHbar[104] + 
   armHHbar[119];
   armHHbar[104]=armHHbar[254] + 9./2.*armHHbar[104] + armHHbar[252];
   armHHbar[104]=MMH*armHHbar[104];
   armHHbar[121]=armHHbar[17]*armHHbar[247];
   armHHbar[104]=1./2.*armHHbar[104] + armHHbar[258] + 3./2.*
   armHHbar[121] + armHHbar[257];
   armHHbar[121]=2*armHHbar[267] + armHHbar[260] + 33./2.*armHHbar[259]
    + armHHbar[264];
   armHHbar[121]=3*armHHbar[121] + 2*armHHbar[268];
   armHHbar[142]= - 3./2.*armHHbar[12];
   armHHbar[156]=armHHbar[121] + armHHbar[142];
   armHHbar[156]=armHHbar[28]*armHHbar[156];
   armHHbar[175]=armHHbar[18]*armHHbar[261];
   armHHbar[181]=armHHbar[19]*armHHbar[261];
   armHHbar[156]=armHHbar[311] + 3./4.*armHHbar[175] + 3./2.*
   armHHbar[181] + armHHbar[156];
   armHHbar[156]=armHHbar[3]*armHHbar[156];
   armHHbar[121]=armHHbar[27]*armHHbar[121];
   armHHbar[121]=armHHbar[293] + 3./2.*armHHbar[285] + armHHbar[121];
   armHHbar[121]=MMt*armHHbar[3]*armHHbar[121];
   armHHbar[175]=armHHbar[148] + armHHbar[191];
   armHHbar[175]=3*armHHbar[175] + armHHbar[193];
   armHHbar[175]=armHHbar[27]*armHHbar[175];
   armHHbar[185]= - armHHbar[29] + armHHbar[30];
   armHHbar[220]=armHHbar[12]*armHHbar[185];
   armHHbar[252]=armHHbar[11]*armHHbar[185];
   armHHbar[175]=1./2.*armHHbar[252] + armHHbar[220] + armHHbar[175];
   armHHbar[121]=armHHbar[121] + 1./2.*armHHbar[175] + armHHbar[156];
   armHHbar[121]=MMt*armHHbar[121];
   armHHbar[156]=armHHbar[27]*armHHbar[247];
   armHHbar[175]= - 2*armHHbar[30];
   armHHbar[220]=armHHbar[9] + armHHbar[29] + armHHbar[175];
   armHHbar[247]=armHHbar[3]*armHHbar[28]*armHHbar[220];
   armHHbar[220]=MMt*armHHbar[3]*armHHbar[27]*armHHbar[220];
   armHHbar[156]=armHHbar[220] + 1./2.*armHHbar[156] + armHHbar[247];
   armHHbar[156]=armHHbar[26]*MMt*armHHbar[156];
   armHHbar[104]=3*armHHbar[156] + 1./2.*armHHbar[104] + armHHbar[121];
   armHHbar[104]=armHHbar[26]*armHHbar[104];
   armHHbar[121]=armHHbar[7]*armHHbar[10];
   armHHbar[121]=armHHbar[121] + armHHbar[123];
   armHHbar[156]=armHHbar[31]*armHHbar[121];
   armHHbar[121]=armHHbar[1]*armHHbar[121];
   armHHbar[121]=3*armHHbar[156] + armHHbar[121];
   armHHbar[156]=armHHbar[191] + 1./3.*armHHbar[193];
   armHHbar[191]=armHHbar[12]*armHHbar[156];
   armHHbar[193]=armHHbar[11]*armHHbar[156];
   armHHbar[121]=1./2.*armHHbar[193] + 3./2.*armHHbar[121] + 
   armHHbar[191];
   armHHbar[121]=MMH*armHHbar[121];
   armHHbar[183]=armHHbar[17]*armHHbar[183];
   armHHbar[191]=armHHbar[31]*armHHbar[183];
   armHHbar[183]=armHHbar[1]*armHHbar[183];
   armHHbar[183]=3*armHHbar[191] + armHHbar[183];
   armHHbar[191]=armHHbar[19]*armHHbar[156];
   armHHbar[156]=armHHbar[18]*armHHbar[156];
   armHHbar[121]=1./2.*armHHbar[121] + 1./2.*armHHbar[156] + 1./2.*
   armHHbar[183] + armHHbar[191];
   armHHbar[156]=5./3.*armHHbar[38];
   armHHbar[130]= - 5./3.*armHHbar[36] + 33./8.*armHHbar[130] + 
   armHHbar[156];
   armHHbar[183]=armHHbar[12]*armHHbar[130];
   armHHbar[130]=armHHbar[11]*armHHbar[130];
   armHHbar[191]=armHHbar[7]*armHHbar[148];
   armHHbar[130]=1./2.*armHHbar[130] + 9./4.*armHHbar[191] + 
   armHHbar[183];
   armHHbar[130]=MMH*armHHbar[130];
   armHHbar[140]=armHHbar[36] + armHHbar[140] - armHHbar[38];
   armHHbar[183]=1./2.*armHHbar[140] + armHHbar[310];
   armHHbar[183]=armHHbar[19]*armHHbar[183];
   armHHbar[191]=5./3.*armHHbar[12];
   armHHbar[140]=1./4.*armHHbar[140] + armHHbar[191];
   armHHbar[140]=armHHbar[18]*armHHbar[140];
   armHHbar[148]=armHHbar[17]*armHHbar[148];
   armHHbar[193]= - armHHbar[19] + armHHbar[18];
   armHHbar[193]=armHHbar[11]*armHHbar[193];
   armHHbar[130]=1./2.*armHHbar[130] + 5./6.*armHHbar[193] + 
   armHHbar[140] + 3./4.*armHHbar[148] + armHHbar[183];
   armHHbar[130]=armHHbar[68]*armHHbar[130];
   armHHbar[140]= - armHHbar[19] - armHHbar[18];
   armHHbar[140]=armHHbar[18]*armHHbar[140];
   armHHbar[140]=armHHbar[299] + 1./2.*armHHbar[140];
   armHHbar[140]=armHHbar[3]*armHHbar[68]*armHHbar[140];
   armHHbar[104]=armHHbar[104] + 17./4.*armHHbar[140] + 1./2.*
   armHHbar[121] + armHHbar[130];
   armHHbar[104]=armHHbar[256]*armHHbar[104];
   armHHbar[121]= - armHHbar[14] + armHHbar[163];
   armHHbar[130]= - armHHbar[7]*armHHbar[10];
   armHHbar[140]=9./4.*armHHbar[194] + armHHbar[121] + 9./4.*
   armHHbar[130];
   armHHbar[140]=armHHbar[31]*armHHbar[140];
   armHHbar[121]=3./4.*armHHbar[194] + 1./3.*armHHbar[121] + 3./4.*
   armHHbar[130];
   armHHbar[121]=armHHbar[1]*armHHbar[121];
   armHHbar[130]=armHHbar[249] + armHHbar[188];
   armHHbar[148]=armHHbar[31]*armHHbar[130];
   armHHbar[130]=armHHbar[1]*armHHbar[130];
   armHHbar[130]=armHHbar[148] + 1./3.*armHHbar[130];
   armHHbar[148]=armHHbar[12]*armHHbar[130];
   armHHbar[130]=armHHbar[11]*armHHbar[130];
   armHHbar[121]=1./2.*armHHbar[130] + armHHbar[148] + armHHbar[140] + 
   armHHbar[121];
   armHHbar[121]=MMH*armHHbar[121];
   armHHbar[130]=armHHbar[17]*armHHbar[266];
   armHHbar[140]=armHHbar[31]*armHHbar[130];
   armHHbar[130]=armHHbar[1]*armHHbar[130];
   armHHbar[130]=3*armHHbar[140] + armHHbar[130];
   armHHbar[140]=armHHbar[9] + 1 - armHHbar[10];
   armHHbar[148]=armHHbar[31]*armHHbar[140];
   armHHbar[140]=armHHbar[1]*armHHbar[140];
   armHHbar[140]=armHHbar[148] + 1./3.*armHHbar[140];
   armHHbar[148]=armHHbar[19]*armHHbar[140];
   armHHbar[140]=armHHbar[18]*armHHbar[140];
   armHHbar[121]=armHHbar[121] + 1./2.*armHHbar[140] + 1./2.*
   armHHbar[130] + armHHbar[148];
   armHHbar[130]=armHHbar[17] - 199*armHHbar[19];
   armHHbar[130]=armHHbar[19]*armHHbar[130];
   armHHbar[140]=149./6.*armHHbar[18] + 1./6.*armHHbar[17] + 17*
   armHHbar[19];
   armHHbar[140]=armHHbar[18]*armHHbar[140];
   armHHbar[130]=1./3.*armHHbar[130] + armHHbar[140];
   armHHbar[130]=armHHbar[3]*armHHbar[68]*armHHbar[130];
   armHHbar[102]=armHHbar[104] + armHHbar[108] + 1./8.*armHHbar[130] + 
   1./2.*armHHbar[121] + armHHbar[102];
   armHHbar[102]=armHHbar[256]*armHHbar[102];
   armHHbar[104]= - 745./6. + 27*armHHbar[59];
   armHHbar[108]= - 11*armHHbar[60];
   armHHbar[104]=armHHbar[289] + armHHbar[234] - 155./6.*armHHbar[63]
    + armHHbar[284] + armHHbar[108] + armHHbar[274] - 149./12.*
   armHHbar[43] - 11./2.*armHHbar[37] + armHHbar[248] + 11*armHHbar[58]
    + armHHbar[238] + 1./4.*armHHbar[104] + armHHbar[295];
   armHHbar[121]= - 15./2.*armHHbar[41];
   armHHbar[130]= - 55*armHHbar[37];
   armHHbar[140]= - 47*armHHbar[40];
   armHHbar[148]=armHHbar[140] + armHHbar[130] - 557./3. + 
   armHHbar[121];
   armHHbar[163]= - 26*armHHbar[39];
   armHHbar[183]=11./3.*armHHbar[9] - 20./9. + 3*armHHbar[10];
   armHHbar[183]=2*armHHbar[31]*armHHbar[183];
   armHHbar[194]=armHHbar[244] - 4./3. + armHHbar[10];
   armHHbar[194]=2*armHHbar[1]*armHHbar[194];
   armHHbar[220]= - 25./3.*armHHbar[27];
   armHHbar[148]=armHHbar[220] + armHHbar[180] + armHHbar[194] + 
   armHHbar[183] + armHHbar[177] + armHHbar[304] + 1./2.*armHHbar[148]
    + armHHbar[163];
   armHHbar[148]=armHHbar[27]*armHHbar[148];
   armHHbar[167]= - 22*armHHbar[18] - 50*armHHbar[28] + armHHbar[167]
    - 11*armHHbar[33] + 3*armHHbar[187] + 22*armHHbar[67];
   armHHbar[167]=armHHbar[3]*armHHbar[167];
   armHHbar[148]=armHHbar[167] + armHHbar[104] + armHHbar[148];
   armHHbar[148]=armHHbar[3]*armHHbar[148];
   armHHbar[187]= - 2./3.*armHHbar[47] - armHHbar[46];
   armHHbar[168]=3./4.*armHHbar[236] + 11./4.*armHHbar[296] + 3./4.*
   armHHbar[168] + 6*armHHbar[50] + 27./4.*armHHbar[53] + 8*
   armHHbar[187] + 9./2.*armHHbar[51];
   armHHbar[187]=2*armHHbar[47] + armHHbar[46];
   armHHbar[236]=53 + armHHbar[295];
   armHHbar[236]=22*armHHbar[63] + 1./2.*armHHbar[236] + 3*armHHbar[64]
   ;
   armHHbar[236]=armHHbar[3]*armHHbar[236];
   armHHbar[187]=armHHbar[236] - 11*armHHbar[53] + 32./3.*armHHbar[187]
    + armHHbar[278];
   armHHbar[187]=MMt*armHHbar[3]*armHHbar[187];
   armHHbar[148]=armHHbar[187] + armHHbar[168] + armHHbar[148];
   armHHbar[148]=MMt*armHHbar[148];
   armHHbar[236]= - 59./2. - 3*armHHbar[41];
   armHHbar[236]=1./2.*armHHbar[236] - 11*armHHbar[37];
   armHHbar[140]=5*armHHbar[236] + armHHbar[140];
   armHHbar[140]=armHHbar[305] + armHHbar[194] + armHHbar[183] + 
   armHHbar[177] + armHHbar[304] + 1./2.*armHHbar[140] + armHHbar[163];
   armHHbar[140]=armHHbar[28]*armHHbar[140];
   armHHbar[236]=3*armHHbar[34];
   armHHbar[238]= - 2./3.*armHHbar[48] + armHHbar[236];
   armHHbar[224]= - 149./12.*armHHbar[20] + armHHbar[306] - 679./24.*
   armHHbar[45] + 177./8.*armHHbar[33] - 58./3.*armHHbar[67] + 2*
   armHHbar[238] + armHHbar[224];
   armHHbar[238]= - 133./2. - 32*armHHbar[29];
   armHHbar[238]=1./3.*armHHbar[238] + armHHbar[307];
   armHHbar[238]=armHHbar[19]*armHHbar[238];
   armHHbar[244]=12*armHHbar[19] - 55./2.*armHHbar[28];
   armHHbar[244]=armHHbar[27]*armHHbar[244];
   armHHbar[247]= - 1./4.*armHHbar[27] + armHHbar[270] + 8 - 43./3.*
   armHHbar[29];
   armHHbar[247]=armHHbar[18]*armHHbar[247];
   armHHbar[140]=armHHbar[247] + armHHbar[244] + armHHbar[140] + 
   armHHbar[224] + armHHbar[238];
   armHHbar[140]=armHHbar[3]*armHHbar[140];
   armHHbar[238]= - 21*armHHbar[62] + 425./4. + 9*armHHbar[59];
   armHHbar[219]=19./2.*armHHbar[58] + 1./2.*armHHbar[238] + 
   armHHbar[219];
   armHHbar[219]=11./2.*armHHbar[237] + armHHbar[303] - armHHbar[30] - 
   11./3.*armHHbar[29] + armHHbar[280] + 83./12.*armHHbar[63] + 
   armHHbar[298] + 9*armHHbar[60] + 5./2.*armHHbar[55] + armHHbar[64]
    + 1./2.*armHHbar[219] - 7*armHHbar[43];
   armHHbar[219]=1./2.*armHHbar[219];
   armHHbar[237]=3*armHHbar[41];
   armHHbar[238]=541./3. + armHHbar[237];
   armHHbar[248]=11*armHHbar[37];
   armHHbar[114]=armHHbar[114] + 1./2.*armHHbar[238] + armHHbar[248];
   armHHbar[114]=armHHbar[109] + armHHbar[154] + 1./2.*armHHbar[114] + 
   armHHbar[145];
   armHHbar[238]=2*armHHbar[27];
   armHHbar[114]=armHHbar[238] + armHHbar[142] + armHHbar[184] + 1./2.*
   armHHbar[114] + armHHbar[179];
   armHHbar[114]=armHHbar[27]*armHHbar[114];
   armHHbar[249]=71 + 53*armHHbar[29];
   armHHbar[249]=armHHbar[312] + 1./9.*armHHbar[249] + armHHbar[207];
   armHHbar[249]=1./2.*armHHbar[11]*armHHbar[249];
   armHHbar[252]=7 + 16*armHHbar[29];
   armHHbar[252]=2./9.*armHHbar[252] + armHHbar[309];
   armHHbar[252]=armHHbar[12]*armHHbar[252];
   armHHbar[254]= - 3./2.*armHHbar[98] - 11*armHHbar[99];
   armHHbar[254]=1./2.*armHHbar[28]*armHHbar[254];
   armHHbar[221]=11./2.*armHHbar[221];
   armHHbar[255]=4./3.*armHHbar[46] - armHHbar[53];
   armHHbar[255]=MMH*armHHbar[255];
   armHHbar[114]=armHHbar[148] + armHHbar[140] + armHHbar[255] + 
   armHHbar[249] + armHHbar[221] + armHHbar[114] + armHHbar[254] + 
   armHHbar[315] + armHHbar[219] + armHHbar[252];
   armHHbar[114]=MMt*armHHbar[114];
   armHHbar[140]=65./8. - armHHbar[13];
   armHHbar[148]= - armHHbar[97]*armHHbar[45];
   armHHbar[252]= - 5./8.*armHHbar[63];
   armHHbar[140]=5./2.*armHHbar[148] + armHHbar[172] + armHHbar[252] + 
   1./3.*armHHbar[140] - 5./8.*armHHbar[55];
   armHHbar[148]= - 17*armHHbar[29];
   armHHbar[172]=215./6. + armHHbar[148];
   armHHbar[257]= - 5./18.*armHHbar[9];
   armHHbar[172]=5*armHHbar[206] + armHHbar[257] + 1./18.*armHHbar[172]
    - armHHbar[30];
   armHHbar[172]=armHHbar[11]*armHHbar[172];
   armHHbar[258]= - 17./4.*armHHbar[29];
   armHHbar[259]=37./3. + armHHbar[258];
   armHHbar[207]= - 5./36.*armHHbar[9] + 1./9.*armHHbar[259] + 
   armHHbar[207];
   armHHbar[207]=armHHbar[12]*armHHbar[207];
   armHHbar[259]=5./3. - 17./8.*armHHbar[29];
   armHHbar[260]=armHHbar[7]*armHHbar[259];
   armHHbar[123]=1./4.*armHHbar[172] + 5./4.*armHHbar[196] + 
   armHHbar[207] + 5./8.*armHHbar[123] + 9./4.*armHHbar[308] + 1./2.*
   armHHbar[140] + armHHbar[260];
   armHHbar[123]=1./2.*MMH*armHHbar[123];
   armHHbar[140]=10./3. + armHHbar[258];
   armHHbar[140]= - 5./12.*armHHbar[9] + 1./3.*armHHbar[140] + 
   armHHbar[201];
   armHHbar[140]=armHHbar[27]*armHHbar[140];
   armHHbar[172]= - 40./3. + 17*armHHbar[29];
   armHHbar[172]=5./3.*armHHbar[9] + 1./3.*armHHbar[172] + 
   armHHbar[319];
   armHHbar[207]=armHHbar[3]*armHHbar[28]*armHHbar[172];
   armHHbar[172]=MMt*armHHbar[3]*armHHbar[27]*armHHbar[172];
   armHHbar[140]=armHHbar[172] + armHHbar[140] + armHHbar[207];
   armHHbar[140]=MMt*armHHbar[140];
   armHHbar[172]=armHHbar[26]*armHHbar[140];
   armHHbar[207]= - 5./24.*armHHbar[9] + 1./3.*armHHbar[259] + 
   armHHbar[234];
   armHHbar[207]=armHHbar[17]*armHHbar[207];
   armHHbar[234]=499./3. + armHHbar[148];
   armHHbar[234]=armHHbar[257] + 1./18.*armHHbar[234] - armHHbar[30];
   armHHbar[234]=armHHbar[19]*armHHbar[234];
   armHHbar[207]=1./2.*armHHbar[234] + armHHbar[207] - 1./32.*
   armHHbar[20] + armHHbar[323] + 79./16.*armHHbar[45] - 5./32.*
   armHHbar[33] - armHHbar[61] + 1./6.*armHHbar[67];
   armHHbar[234]=241./12. + armHHbar[148];
   armHHbar[258]=5./2.*armHHbar[206];
   armHHbar[234]=1./4.*armHHbar[27] + armHHbar[258] + armHHbar[257] + 1.
   /18.*armHHbar[234] - armHHbar[30];
   armHHbar[234]=1./4.*armHHbar[18]*armHHbar[234];
   armHHbar[257]= - 143./16. + 41./3.*armHHbar[12];
   armHHbar[259]=5./8.*armHHbar[196];
   armHHbar[257]=1./3.*armHHbar[257] + armHHbar[259];
   armHHbar[257]=armHHbar[28]*armHHbar[257];
   armHHbar[260]= - 107*armHHbar[19] - 185*armHHbar[28];
   armHHbar[260]=armHHbar[28]*armHHbar[260];
   armHHbar[266]= - armHHbar[18]*armHHbar[28];
   armHHbar[260]=armHHbar[260] + 71*armHHbar[266];
   armHHbar[260]=armHHbar[3]*armHHbar[260];
   armHHbar[267]= - armHHbar[19] + 13./8.*armHHbar[28];
   armHHbar[267]=3*armHHbar[27]*armHHbar[267];
   armHHbar[268]=71./18.*armHHbar[217];
   armHHbar[114]=armHHbar[172] + armHHbar[114] + 1./12.*armHHbar[260]
    + armHHbar[123] + armHHbar[268] + armHHbar[234] + armHHbar[267] + 
   armHHbar[207] + armHHbar[257];
   armHHbar[114]=armHHbar[26]*armHHbar[114];
   armHHbar[172]=armHHbar[17] + 233*armHHbar[19];
   armHHbar[172]=armHHbar[19]*armHHbar[172];
   armHHbar[172]=1./24.*armHHbar[172] + armHHbar[176];
   armHHbar[172]=armHHbar[3]*armHHbar[68]*armHHbar[172];
   armHHbar[102]=armHHbar[102] + armHHbar[114] + armHHbar[172] + 
   armHHbar[115] + armHHbar[103];
   armHHbar[102]=armHHbar[256]*armHHbar[102];
   armHHbar[103]= - 7./2. + 9*armHHbar[54];
   armHHbar[114]= - 5./4.*armHHbar[59];
   armHHbar[115]=1./2.*armHHbar[35];
   armHHbar[172]=3./2.*armHHbar[56];
   armHHbar[176]=9./4.*armHHbar[43];
   armHHbar[103]=armHHbar[298] + armHHbar[231] + armHHbar[176] - 
   armHHbar[58] + armHHbar[172] + armHHbar[62] + armHHbar[115] + 
   armHHbar[113] + armHHbar[114] + 1./2.*armHHbar[103] - 3*armHHbar[65]
   ;
   armHHbar[257]= - 5./2.*armHHbar[66] + 11*armHHbar[32];
   armHHbar[236]=3./4.*armHHbar[21] - 33./2.*armHHbar[45] + 
   armHHbar[325] - 7./4.*armHHbar[61] - 17./4.*armHHbar[44] + 1./2.*
   armHHbar[257] + armHHbar[236];
   armHHbar[257]=3*armHHbar[236];
   armHHbar[260]= - 3*armHHbar[17];
   armHHbar[274]=armHHbar[260] - armHHbar[19];
   armHHbar[278]= - 18*armHHbar[28];
   armHHbar[274]=5*armHHbar[274] + armHHbar[278];
   armHHbar[274]=armHHbar[27]*armHHbar[274];
   armHHbar[280]= - 27./4.*armHHbar[17];
   armHHbar[175]= - 9./4. + armHHbar[175];
   armHHbar[175]=armHHbar[19]*armHHbar[175];
   armHHbar[284]= - 5./2.*armHHbar[41];
   armHHbar[286]=armHHbar[284] + 3 - 10*armHHbar[35];
   armHHbar[286]=armHHbar[28]*armHHbar[286];
   armHHbar[289]=47./4. - 18*armHHbar[27];
   armHHbar[289]=armHHbar[18]*armHHbar[289];
   armHHbar[175]=armHHbar[289] + armHHbar[274] + 3*armHHbar[286] + 3*
   armHHbar[175] + armHHbar[280] + armHHbar[257] - 41./4.*armHHbar[20];
   armHHbar[175]=armHHbar[3]*armHHbar[175];
   armHHbar[274]=3./4.*armHHbar[41] + 1 + 21./2.*armHHbar[35];
   armHHbar[127]=armHHbar[180] + armHHbar[127] + armHHbar[137] + 
   armHHbar[274] + armHHbar[264];
   armHHbar[127]=armHHbar[27]*armHHbar[127];
   armHHbar[286]=2*armHHbar[66] - armHHbar[32];
   armHHbar[289]= - armHHbar[27]*armHHbar[17];
   armHHbar[186]=4*armHHbar[289] + armHHbar[278] + armHHbar[195] - 4*
   armHHbar[17] + armHHbar[186] - 4*armHHbar[44] + 4*armHHbar[286] - 
   armHHbar[34];
   armHHbar[186]=armHHbar[3]*armHHbar[186];
   armHHbar[195]= - 5*armHHbar[35];
   armHHbar[278]= - 2 + armHHbar[195];
   armHHbar[269]=armHHbar[269] + 2*armHHbar[278] + armHHbar[284];
   armHHbar[269]=armHHbar[27]*armHHbar[269];
   armHHbar[186]=armHHbar[186] + armHHbar[269] - armHHbar[63] - 5./4.*
   armHHbar[64] - 13./2.*armHHbar[43] - 1./2.*armHHbar[41] - 2*
   armHHbar[35] + 5./2.*armHHbar[59] - 5*armHHbar[65] + 3 + 4*
   armHHbar[57];
   armHHbar[186]=armHHbar[3]*armHHbar[186];
   armHHbar[269]=armHHbar[64] + 9 + 8*armHHbar[65];
   armHHbar[269]=armHHbar[3]*armHHbar[269];
   armHHbar[269]= - 4*armHHbar[52] + armHHbar[269];
   armHHbar[269]=MMt*armHHbar[3]*armHHbar[269];
   armHHbar[284]= - 6*armHHbar[49] + armHHbar[52];
   armHHbar[293]=1 + armHHbar[64];
   armHHbar[294]=armHHbar[98]*armHHbar[293];
   armHHbar[186]=armHHbar[269] + armHHbar[186] + 1./4.*armHHbar[294] + 
   2*armHHbar[284] - armHHbar[51];
   armHHbar[186]=MMt*armHHbar[186];
   armHHbar[269]=3*armHHbar[314];
   armHHbar[296]= - 3./2.*armHHbar[28]*armHHbar[98];
   armHHbar[299]= - 27./2.*armHHbar[7];
   armHHbar[301]= - armHHbar[50] + 9*armHHbar[49] + 1./2.*armHHbar[52];
   armHHbar[301]=MMH*armHHbar[301];
   armHHbar[306]=armHHbar[12]*armHHbar[30];
   armHHbar[103]=3*armHHbar[186] + armHHbar[175] + 3./2.*armHHbar[301]
    + armHHbar[127] + armHHbar[296] + armHHbar[269] + armHHbar[306] + 
   armHHbar[303] - armHHbar[30] + armHHbar[299] + 3*armHHbar[103] + 13./
   4.*armHHbar[63];
   armHHbar[103]=MMt*armHHbar[103];
   armHHbar[127]= - 3*armHHbar[57];
   armHHbar[175]=29./4. + armHHbar[127];
   armHHbar[186]= - 3*armHHbar[54];
   armHHbar[307]=3*armHHbar[42];
   armHHbar[310]= - 1./2.*armHHbar[56];
   armHHbar[311]=3./4.*armHHbar[43];
   armHHbar[175]=armHHbar[311] + armHHbar[310] + armHHbar[307] + 1./4.*
   armHHbar[175] + armHHbar[186];
   armHHbar[314]=1 - 3./2.*armHHbar[35];
   armHHbar[317]=armHHbar[36] + armHHbar[314] + armHHbar[38];
   armHHbar[319]=9*armHHbar[7];
   armHHbar[320]=3./2.*armHHbar[12];
   armHHbar[317]=armHHbar[320] + 1./2.*armHHbar[317] + armHHbar[319];
   armHHbar[317]=armHHbar[27]*armHHbar[317];
   armHHbar[321]= - 7./2.*armHHbar[29];
   armHHbar[323]= - 71 + armHHbar[321];
   armHHbar[323]=1./9.*armHHbar[323] + armHHbar[242];
   armHHbar[323]=1./2.*armHHbar[323] + armHHbar[233];
   armHHbar[323]=armHHbar[11]*armHHbar[323];
   armHHbar[324]=1 - 1./8.*armHHbar[29];
   armHHbar[325]=armHHbar[7]*armHHbar[324];
   armHHbar[324]=7./9.*armHHbar[324] + 1./8.*armHHbar[30];
   armHHbar[324]=armHHbar[12]*armHHbar[324];
   armHHbar[175]=1./4.*armHHbar[323] + armHHbar[317] + armHHbar[324] + 
   9./16.*armHHbar[119] + 7./2.*armHHbar[325] + 11./8.*armHHbar[15] + 
   armHHbar[252] + armHHbar[230] + 3*armHHbar[175] - 11./8.*
   armHHbar[55];
   armHHbar[175]=MMH*armHHbar[175];
   armHHbar[252]=5./2.*armHHbar[12] + armHHbar[299] + armHHbar[137] + 
   armHHbar[264] + 8 + 15./2.*armHHbar[35];
   armHHbar[252]=armHHbar[28]*armHHbar[252];
   armHHbar[317]=1 + armHHbar[321];
   armHHbar[317]=1./3.*armHHbar[317] + armHHbar[212];
   armHHbar[317]=armHHbar[17]*armHHbar[317];
   armHHbar[323]=9*armHHbar[17];
   armHHbar[324]=9*armHHbar[28];
   armHHbar[325]=armHHbar[324] + armHHbar[323] + 5./2.*armHHbar[19];
   armHHbar[325]=armHHbar[27]*armHHbar[325];
   armHHbar[326]= - 9*armHHbar[17];
   armHHbar[327]=armHHbar[326] - 53*armHHbar[19];
   armHHbar[328]= - 3*armHHbar[28];
   armHHbar[327]=1./4.*armHHbar[327] + armHHbar[328];
   armHHbar[327]=armHHbar[28]*armHHbar[327];
   armHHbar[327]=armHHbar[327] + 31./4.*armHHbar[266];
   armHHbar[327]=armHHbar[3]*armHHbar[327];
   armHHbar[150]=armHHbar[150] - 3*armHHbar[32] - armHHbar[34];
   armHHbar[150]=1./2.*armHHbar[150] - armHHbar[61];
   armHHbar[329]=3*armHHbar[150];
   armHHbar[330]= - 3./4.*armHHbar[21];
   armHHbar[331]=35 - armHHbar[29];
   armHHbar[331]=7./9.*armHHbar[331] + armHHbar[30];
   armHHbar[331]=armHHbar[19]*armHHbar[331];
   armHHbar[332]= - 73 - armHHbar[29];
   armHHbar[332]=7./9.*armHHbar[332] + armHHbar[30];
   armHHbar[332]=1./4.*armHHbar[332] + 9*armHHbar[27];
   armHHbar[332]=armHHbar[18]*armHHbar[332];
   armHHbar[103]=armHHbar[103] + armHHbar[327] + armHHbar[175] + 1./2.*
   armHHbar[332] + 1./2.*armHHbar[325] + armHHbar[252] + 1./4.*
   armHHbar[331] + 1./4.*armHHbar[317] + 27./8.*armHHbar[20] + 
   armHHbar[330] + 7*armHHbar[45] - 11./8.*armHHbar[33] + armHHbar[329]
    + 5./2.*armHHbar[67];
   armHHbar[103]=MMt*armHHbar[103];
   armHHbar[175]= - 1./3.*armHHbar[19];
   armHHbar[128]=armHHbar[128] + armHHbar[175];
   armHHbar[252]=13./3.*armHHbar[128] + 5./8.*armHHbar[28];
   armHHbar[252]=armHHbar[28]*armHHbar[252];
   armHHbar[317]=5./2. - 13*armHHbar[7];
   armHHbar[317]=1./2.*armHHbar[317] - 13./9.*armHHbar[12];
   armHHbar[317]=armHHbar[28]*armHHbar[317];
   armHHbar[317]=71./36.*armHHbar[251] - 5./4.*armHHbar[45] + 
   armHHbar[317];
   armHHbar[317]=MMH*armHHbar[317];
   armHHbar[103]=armHHbar[103] + 1./2.*armHHbar[317] + armHHbar[252] + 
   97./72.*armHHbar[266];
   armHHbar[252]=1./2.*armHHbar[314];
   armHHbar[314]=armHHbar[320] + armHHbar[319] + armHHbar[252] + 
   armHHbar[38];
   armHHbar[314]=armHHbar[27]*armHHbar[314];
   armHHbar[317]=19./3. + armHHbar[30];
   armHHbar[317]=1./2.*armHHbar[317] + armHHbar[233];
   armHHbar[317]=armHHbar[11]*armHHbar[317];
   armHHbar[325]=127./4. - 9*armHHbar[57];
   armHHbar[327]=1 + armHHbar[309];
   armHHbar[327]=armHHbar[12]*armHHbar[327];
   armHHbar[229]=1./4.*armHHbar[317] + armHHbar[314] + armHHbar[327] + 
   9./8.*armHHbar[119] + armHHbar[214] - 7./24.*armHHbar[15] + 25./24.*
   armHHbar[63] + armHHbar[230] + 7./24.*armHHbar[55] + armHHbar[176]
    + armHHbar[229] + 9*armHHbar[42] + 1./4.*armHHbar[325] - 9*
   armHHbar[54];
   armHHbar[229]=MMH*armHHbar[229];
   armHHbar[314]= - 61./2. + 27*armHHbar[54];
   armHHbar[317]=1./2.*armHHbar[306];
   armHHbar[161]=armHHbar[296] + armHHbar[269] + armHHbar[317] + 
   armHHbar[303] - armHHbar[30] + armHHbar[299] - 41./12.*armHHbar[63]
    - 9./2.*armHHbar[16] + 3./2.*armHHbar[64] + 27./4.*armHHbar[43] + 
   armHHbar[161] + 9./2.*armHHbar[56] + armHHbar[295] + 3./2.*
   armHHbar[35] - 27./2.*armHHbar[42] - 15./4.*armHHbar[59] + 1./2.*
   armHHbar[314] - 9*armHHbar[65];
   armHHbar[295]= - 3./2. - armHHbar[30];
   armHHbar[295]=armHHbar[19]*armHHbar[295];
   armHHbar[257]=9./2.*armHHbar[295] + armHHbar[280] + armHHbar[257] + 
   37./12.*armHHbar[20];
   armHHbar[280]=7*armHHbar[19];
   armHHbar[295]=armHHbar[260] + armHHbar[280];
   armHHbar[314]= - 9*armHHbar[28];
   armHHbar[295]=5./2.*armHHbar[295] + armHHbar[314];
   armHHbar[295]=armHHbar[27]*armHHbar[295];
   armHHbar[325]= - 19./6. + armHHbar[192];
   armHHbar[327]=17./3.*armHHbar[27];
   armHHbar[325]=1./4.*armHHbar[325] + armHHbar[327];
   armHHbar[325]=armHHbar[18]*armHHbar[325];
   armHHbar[331]= - 5./4.*armHHbar[41];
   armHHbar[195]=armHHbar[331] + 3./2. + armHHbar[195];
   armHHbar[195]=armHHbar[28]*armHHbar[195];
   armHHbar[257]=armHHbar[325] + armHHbar[295] + 1./2.*armHHbar[257] + 
   3*armHHbar[195];
   armHHbar[257]=armHHbar[3]*armHHbar[257];
   armHHbar[274]=1./2.*armHHbar[274];
   armHHbar[295]=armHHbar[142] + armHHbar[299] + armHHbar[274] + 
   armHHbar[264];
   armHHbar[295]=armHHbar[27]*armHHbar[295];
   armHHbar[286]=2*armHHbar[289] + armHHbar[314] - armHHbar[19] - 2*
   armHHbar[17] + armHHbar[61] - 2*armHHbar[44] + 2*armHHbar[286] + 
   armHHbar[300];
   armHHbar[286]=armHHbar[3]*armHHbar[286];
   armHHbar[278]= - armHHbar[27] + armHHbar[278] + armHHbar[331];
   armHHbar[278]=armHHbar[27]*armHHbar[278];
   armHHbar[226]=armHHbar[286] + armHHbar[278] + armHHbar[226] - 5./8.*
   armHHbar[64] - 13./4.*armHHbar[43] - 1./4.*armHHbar[41] - 
   armHHbar[35] + 5./4.*armHHbar[59] - 5./2.*armHHbar[65] + 3./2. + 2*
   armHHbar[57];
   armHHbar[226]=armHHbar[3]*armHHbar[226];
   armHHbar[278]=armHHbar[231] + 9./2. + 4*armHHbar[65];
   armHHbar[278]=armHHbar[3]*armHHbar[278];
   armHHbar[278]= - 2*armHHbar[52] + armHHbar[278];
   armHHbar[278]=MMt*armHHbar[3]*armHHbar[278];
   armHHbar[226]=armHHbar[278] + armHHbar[226] + 1./8.*armHHbar[294] + 
   armHHbar[284] - 1./2.*armHHbar[51];
   armHHbar[226]=3*MMt*armHHbar[226];
   armHHbar[278]=3./4.*armHHbar[301];
   armHHbar[161]=armHHbar[226] + armHHbar[257] + armHHbar[278] + 1./4.*
   armHHbar[290] + 1./2.*armHHbar[161] + armHHbar[295];
   armHHbar[161]=MMt*armHHbar[161];
   armHHbar[257]=1 + armHHbar[30];
   armHHbar[257]=armHHbar[17]*armHHbar[257];
   armHHbar[284]=29./2. + armHHbar[30];
   armHHbar[284]=armHHbar[19]*armHHbar[284];
   armHHbar[257]=1./2.*armHHbar[284] + 3./4.*armHHbar[257] - 13./8.*
   armHHbar[20] + armHHbar[330] + 31./3.*armHHbar[45] + 7./24.*
   armHHbar[33] + armHHbar[329] - 25./6.*armHHbar[67];
   armHHbar[284]= - 27./4.*armHHbar[7];
   armHHbar[264]=armHHbar[12] + armHHbar[284] + armHHbar[264] + 37./3.
    + 15./4.*armHHbar[35];
   armHHbar[264]=armHHbar[28]*armHHbar[264];
   armHHbar[286]=armHHbar[324] + armHHbar[323] - 35./2.*armHHbar[19];
   armHHbar[286]=armHHbar[27]*armHHbar[286];
   armHHbar[294]=65./2. + armHHbar[30];
   armHHbar[295]= - 17./3.*armHHbar[27];
   armHHbar[294]=1./2.*armHHbar[294] + armHHbar[295];
   armHHbar[294]=armHHbar[18]*armHHbar[294];
   armHHbar[300]=armHHbar[326] + 113*armHHbar[19];
   armHHbar[300]=1./4.*armHHbar[300] + armHHbar[328];
   armHHbar[300]=armHHbar[28]*armHHbar[300];
   armHHbar[300]=armHHbar[300] + 39./4.*armHHbar[322];
   armHHbar[300]=armHHbar[3]*armHHbar[300];
   armHHbar[161]=armHHbar[161] + 1./2.*armHHbar[300] + 1./2.*
   armHHbar[229] + 1./4.*armHHbar[251] + 1./4.*armHHbar[294] + 1./4.*
   armHHbar[286] + 1./2.*armHHbar[257] + armHHbar[264];
   armHHbar[161]=MMt*armHHbar[161];
   armHHbar[229]= - armHHbar[18]*armHHbar[27];
   armHHbar[257]=armHHbar[27]*armHHbar[19];
   armHHbar[264]=armHHbar[257] + armHHbar[229];
   armHHbar[264]=1./2.*armHHbar[264] + armHHbar[251];
   armHHbar[286]= - armHHbar[28]*armHHbar[19];
   armHHbar[294]=armHHbar[286] + 5./2.*armHHbar[322];
   armHHbar[294]=armHHbar[3]*armHHbar[294];
   armHHbar[300]=armHHbar[212] + armHHbar[27];
   armHHbar[300]=armHHbar[18]*armHHbar[300];
   armHHbar[301]= - armHHbar[27]*armHHbar[19];
   armHHbar[300]=armHHbar[301] + armHHbar[300];
   armHHbar[300]=armHHbar[3]*armHHbar[300];
   armHHbar[300]=1./2.*armHHbar[290] + armHHbar[300];
   armHHbar[300]=MMt*armHHbar[300];
   armHHbar[264]=armHHbar[300] + 1./2.*armHHbar[264] + armHHbar[294];
   armHHbar[264]=armHHbar[171]*MMt*armHHbar[264];
   armHHbar[294]= - 25./6. - 9*armHHbar[7];
   armHHbar[294]=1./2.*armHHbar[294] - armHHbar[12];
   armHHbar[294]=armHHbar[28]*armHHbar[294];
   armHHbar[300]=19./12.*armHHbar[217];
   armHHbar[294]=armHHbar[300] + 25./12.*armHHbar[45] + armHHbar[294];
   armHHbar[294]=MMH*armHHbar[294];
   armHHbar[314]= - 25./24.*armHHbar[28] + armHHbar[292] - armHHbar[19]
   ;
   armHHbar[314]=armHHbar[28]*armHHbar[314];
   armHHbar[294]=1./2.*armHHbar[294] + armHHbar[314] + 13./24.*
   armHHbar[322];
   armHHbar[161]=1./2.*armHHbar[264] + 1./2.*armHHbar[294] + 
   armHHbar[161];
   armHHbar[161]=armHHbar[171]*armHHbar[161];
   armHHbar[161]=armHHbar[103] + armHHbar[161];
   armHHbar[161]=armHHbar[171]*armHHbar[161];
   armHHbar[264]= - 13 - 7./4.*armHHbar[29];
   armHHbar[264]=1./3.*armHHbar[264] + armHHbar[270];
   armHHbar[264]=armHHbar[27]*armHHbar[264];
   armHHbar[270]=52 + 7*armHHbar[29];
   armHHbar[197]=1./3.*armHHbar[270] + armHHbar[197];
   armHHbar[270]=armHHbar[28]*armHHbar[197];
   armHHbar[270]=armHHbar[270] + 52./3.*armHHbar[218];
   armHHbar[270]=armHHbar[3]*armHHbar[270];
   armHHbar[197]=MMt*armHHbar[3]*armHHbar[27]*armHHbar[197];
   armHHbar[197]=armHHbar[197] + armHHbar[264] + armHHbar[270];
   armHHbar[197]=MMt*armHHbar[197];
   armHHbar[264]=pow(armHHbar[28],2);
   armHHbar[270]=armHHbar[3]*armHHbar[264];
   armHHbar[294]=armHHbar[265] + 4*armHHbar[270];
   armHHbar[197]=13./3.*armHHbar[294] + armHHbar[197];
   armHHbar[197]=MMt*armHHbar[197];
   armHHbar[294]=2 - armHHbar[30];
   armHHbar[314]=armHHbar[28]*armHHbar[294];
   armHHbar[314]=armHHbar[314] + 2*armHHbar[218];
   armHHbar[314]=armHHbar[3]*armHHbar[314];
   armHHbar[313]=armHHbar[27]*armHHbar[313];
   armHHbar[294]=MMt*armHHbar[3]*armHHbar[27]*armHHbar[294];
   armHHbar[294]=armHHbar[294] + 1./2.*armHHbar[313] + armHHbar[314];
   armHHbar[294]=MMt*armHHbar[294];
   armHHbar[294]=armHHbar[294] + 1./2.*armHHbar[265] + 2*armHHbar[270];
   armHHbar[294]=armHHbar[171]*MMt*armHHbar[294];
   armHHbar[294]=armHHbar[197] + 3*armHHbar[294];
   armHHbar[294]=armHHbar[26]*armHHbar[171]*armHHbar[294];
   armHHbar[161]=armHHbar[161] + armHHbar[294];
   armHHbar[161]=armHHbar[26]*armHHbar[161];
   armHHbar[197]=armHHbar[26]*armHHbar[197];
   armHHbar[103]=armHHbar[103] + armHHbar[197];
   armHHbar[103]=armHHbar[26]*armHHbar[103];
   armHHbar[197]= - 5./2. + 3*armHHbar[54];
   armHHbar[197]=1./2.*armHHbar[197] - armHHbar[65];
   armHHbar[113]=armHHbar[125] - 1./4.*armHHbar[63] + armHHbar[298] + 
   armHHbar[231] + armHHbar[176] - armHHbar[58] + armHHbar[172] + 
   armHHbar[62] + armHHbar[115] + armHHbar[113] + 3*armHHbar[197] + 
   armHHbar[114];
   armHHbar[113]=armHHbar[296] + armHHbar[269] + 3./2.*armHHbar[306] + 
   armHHbar[303] + 3*armHHbar[113] - armHHbar[30];
   armHHbar[114]=armHHbar[142] + armHHbar[299] + armHHbar[274] + 
   armHHbar[137];
   armHHbar[114]=armHHbar[27]*armHHbar[114];
   armHHbar[115]= - 9./2. - 5*armHHbar[30];
   armHHbar[115]=armHHbar[19]*armHHbar[115];
   armHHbar[115]=1./2.*armHHbar[115] - 9./4.*armHHbar[17] + 
   armHHbar[236] - 3./4.*armHHbar[20];
   armHHbar[125]= - 5*armHHbar[17] + armHHbar[19];
   armHHbar[125]=1./2.*armHHbar[125] + armHHbar[328];
   armHHbar[125]=armHHbar[27]*armHHbar[125];
   armHHbar[115]=armHHbar[125] + 1./2.*armHHbar[115] + armHHbar[195];
   armHHbar[125]=5./2. - armHHbar[30];
   armHHbar[125]=3./4.*armHHbar[125] - 13*armHHbar[27];
   armHHbar[125]=armHHbar[18]*armHHbar[125];
   armHHbar[115]=3*armHHbar[115] + armHHbar[125];
   armHHbar[115]=armHHbar[3]*armHHbar[115];
   armHHbar[125]=armHHbar[11]*armHHbar[30];
   armHHbar[172]=1./4.*armHHbar[125];
   armHHbar[113]=armHHbar[226] + armHHbar[115] + armHHbar[278] + 
   armHHbar[172] + 1./2.*armHHbar[113] + armHHbar[114];
   armHHbar[113]=MMt*armHHbar[113];
   armHHbar[114]=37./4. + armHHbar[127];
   armHHbar[115]=7 + 19./4.*armHHbar[29];
   armHHbar[127]=armHHbar[7]*armHHbar[115];
   armHHbar[114]=1./2.*armHHbar[127] + 1./8.*armHHbar[15] + 
   armHHbar[297] + 1./2.*armHHbar[16] - 1./8.*armHHbar[55] + 
   armHHbar[311] + armHHbar[310] + armHHbar[307] + 1./4.*armHHbar[114]
    + armHHbar[186];
   armHHbar[127]=armHHbar[320] + armHHbar[319] + armHHbar[252] + 
   armHHbar[36];
   armHHbar[127]=armHHbar[27]*armHHbar[127];
   armHHbar[176]=1 + armHHbar[29];
   armHHbar[186]=19./6.*armHHbar[176] + armHHbar[233];
   armHHbar[186]=armHHbar[11]*armHHbar[186];
   armHHbar[115]=armHHbar[12]*armHHbar[115];
   armHHbar[114]=1./4.*armHHbar[186] + armHHbar[127] + 3*armHHbar[114]
    + 1./3.*armHHbar[115];
   armHHbar[114]=MMH*armHHbar[114];
   armHHbar[115]=3 + 5./4.*armHHbar[35];
   armHHbar[115]=armHHbar[320] + armHHbar[284] + 3*armHHbar[115] + 
   armHHbar[137];
   armHHbar[115]=armHHbar[28]*armHHbar[115];
   armHHbar[127]=1./8.*armHHbar[20] - 1./4.*armHHbar[21] + 3*
   armHHbar[45] - 1./8.*armHHbar[33] + armHHbar[150] + armHHbar[208];
   armHHbar[137]=armHHbar[17]*armHHbar[176];
   armHHbar[150]=19*armHHbar[29];
   armHHbar[186]=119./2. + armHHbar[150];
   armHHbar[186]=armHHbar[19]*armHHbar[186];
   armHHbar[127]=1./6.*armHHbar[186] + 3*armHHbar[127] + 19./4.*
   armHHbar[137];
   armHHbar[137]=armHHbar[209] + armHHbar[166] - 1./2.*armHHbar[19];
   armHHbar[137]=armHHbar[27]*armHHbar[137];
   armHHbar[166]= - armHHbar[17] - 3*armHHbar[19];
   armHHbar[166]=3./4.*armHHbar[166] - armHHbar[28];
   armHHbar[166]=armHHbar[28]*armHHbar[166];
   armHHbar[166]=3*armHHbar[166] + 101./4.*armHHbar[266];
   armHHbar[166]=armHHbar[3]*armHHbar[166];
   armHHbar[150]=83./2. + armHHbar[150];
   armHHbar[150]=1./6.*armHHbar[150] + 13*armHHbar[27];
   armHHbar[150]=armHHbar[18]*armHHbar[150];
   armHHbar[186]=1./4.*armHHbar[217];
   armHHbar[113]=armHHbar[113] + 1./2.*armHHbar[166] + 1./2.*
   armHHbar[114] + armHHbar[186] + 1./4.*armHHbar[150] + 3./4.*
   armHHbar[137] + 1./2.*armHHbar[127] + armHHbar[115];
   armHHbar[113]=MMt*armHHbar[113];
   armHHbar[114]=armHHbar[173] + armHHbar[291];
   armHHbar[114]=5*armHHbar[114] - 3./8.*armHHbar[28];
   armHHbar[114]=armHHbar[28]*armHHbar[114];
   armHHbar[115]= - 1./2. + 5*armHHbar[7];
   armHHbar[115]=3./2.*armHHbar[115] + armHHbar[191];
   armHHbar[115]=armHHbar[28]*armHHbar[115];
   armHHbar[115]=armHHbar[300] + 3./4.*armHHbar[45] + armHHbar[115];
   armHHbar[115]=MMH*armHHbar[115];
   armHHbar[114]=1./2.*armHHbar[115] + armHHbar[114] + 29./24.*
   armHHbar[322];
   armHHbar[115]= - 10 - 19*armHHbar[29];
   armHHbar[127]=armHHbar[28]*armHHbar[115];
   armHHbar[127]=armHHbar[127] + 10*armHHbar[265];
   armHHbar[127]=armHHbar[3]*armHHbar[127];
   armHHbar[137]=5 + 19./2.*armHHbar[29];
   armHHbar[137]=armHHbar[27]*armHHbar[137];
   armHHbar[115]=MMt*armHHbar[3]*armHHbar[27]*armHHbar[115];
   armHHbar[115]=armHHbar[115] + 1./2.*armHHbar[137] + armHHbar[127];
   armHHbar[115]=MMt*armHHbar[115];
   armHHbar[127]= - armHHbar[3]*armHHbar[264];
   armHHbar[137]=1./2.*armHHbar[218] + 2*armHHbar[127];
   armHHbar[115]=5*armHHbar[137] + armHHbar[115];
   armHHbar[115]=armHHbar[26]*MMt*armHHbar[115];
   armHHbar[113]=armHHbar[115] + 1./2.*armHHbar[114] + armHHbar[113];
   armHHbar[113]=armHHbar[26]*armHHbar[113];
   armHHbar[114]=armHHbar[201] + 17*armHHbar[27];
   armHHbar[114]=armHHbar[18]*armHHbar[114];
   armHHbar[115]= - armHHbar[19]*armHHbar[30];
   armHHbar[114]=armHHbar[114] + 3*armHHbar[115] + 17*armHHbar[301];
   armHHbar[114]=armHHbar[3]*armHHbar[114];
   armHHbar[137]=armHHbar[38] - armHHbar[36];
   armHHbar[150]=armHHbar[27]*armHHbar[137];
   armHHbar[114]=1./2.*armHHbar[114] + armHHbar[172] + armHHbar[317] + 
   2*armHHbar[150];
   armHHbar[114]=MMt*armHHbar[114];
   armHHbar[166]=armHHbar[216] + armHHbar[308];
   armHHbar[166]=9./2.*armHHbar[166] + armHHbar[262];
   armHHbar[172]= - armHHbar[38] + armHHbar[36];
   armHHbar[191]=armHHbar[27]*armHHbar[172];
   armHHbar[166]=1./8.*armHHbar[263] + 1./4.*armHHbar[166] + 
   armHHbar[191];
   armHHbar[166]=MMH*armHHbar[166];
   armHHbar[195]=armHHbar[17]*armHHbar[261];
   armHHbar[181]=3./2.*armHHbar[195] + armHHbar[181];
   armHHbar[195]=2*armHHbar[137] + armHHbar[178];
   armHHbar[195]=armHHbar[28]*armHHbar[195];
   armHHbar[197]=10*armHHbar[286] + 31./4.*armHHbar[322];
   armHHbar[197]=armHHbar[3]*armHHbar[197];
   armHHbar[201]=armHHbar[261] - 17*armHHbar[27];
   armHHbar[201]=armHHbar[18]*armHHbar[201];
   armHHbar[114]=armHHbar[114] + armHHbar[197] + 1./2.*armHHbar[166] + 
   armHHbar[186] + 1./8.*armHHbar[201] + 17./8.*armHHbar[257] + 1./4.*
   armHHbar[181] + armHHbar[195];
   armHHbar[114]=MMt*armHHbar[114];
   armHHbar[166]=armHHbar[27]*armHHbar[261];
   armHHbar[181]=armHHbar[3]*armHHbar[28]*armHHbar[185];
   armHHbar[185]=MMt*armHHbar[3]*armHHbar[27]*armHHbar[185];
   armHHbar[166]=armHHbar[185] + 1./4.*armHHbar[166] + armHHbar[181];
   armHHbar[181]=pow(MMt,2);
   armHHbar[166]=armHHbar[26]*armHHbar[181]*armHHbar[166];
   armHHbar[114]=armHHbar[114] + 3*armHHbar[166];
   armHHbar[114]=armHHbar[26]*armHHbar[114];
   armHHbar[166]=armHHbar[17]*armHHbar[172];
   armHHbar[185]=1./3.*armHHbar[137];
   armHHbar[186]=7./4.*armHHbar[12] + armHHbar[185] + 51./8.*
   armHHbar[7];
   armHHbar[186]=armHHbar[19]*armHHbar[186];
   armHHbar[195]= - 7./2.*armHHbar[12] + armHHbar[185] - 51./4.*
   armHHbar[7];
   armHHbar[195]=armHHbar[18]*armHHbar[195];
   armHHbar[166]=7./8.*armHHbar[174] + 1./2.*armHHbar[195] + 
   armHHbar[166] + armHHbar[186];
   armHHbar[174]=armHHbar[7]*armHHbar[172];
   armHHbar[186]=armHHbar[12]*armHHbar[172];
   armHHbar[195]=armHHbar[11]*armHHbar[172];
   armHHbar[174]=1./6.*armHHbar[195] + 3./4.*armHHbar[174] + 1./3.*
   armHHbar[186];
   armHHbar[174]=MMH*armHHbar[174];
   armHHbar[166]=1./2.*armHHbar[166] + armHHbar[174];
   armHHbar[166]=MMH*armHHbar[166];
   armHHbar[174]=11./3.*armHHbar[19];
   armHHbar[186]=17./2.*armHHbar[17] + armHHbar[174];
   armHHbar[186]=armHHbar[19]*armHHbar[186];
   armHHbar[195]= - 11./3.*armHHbar[18] - 17*armHHbar[17] - 11./3.*
   armHHbar[19];
   armHHbar[195]=armHHbar[18]*armHHbar[195];
   armHHbar[186]=armHHbar[186] + 1./2.*armHHbar[195];
   armHHbar[166]=1./4.*armHHbar[186] + armHHbar[166];
   armHHbar[166]=armHHbar[68]*armHHbar[166];
   armHHbar[114]=1./2.*armHHbar[166] + armHHbar[114];
   armHHbar[114]=armHHbar[256]*armHHbar[114];
   armHHbar[166]= - 101./8. - 9*armHHbar[82];
   armHHbar[129]=1./2.*armHHbar[166] + armHHbar[129];
   armHHbar[166]=1./6.*armHHbar[38];
   armHHbar[129]=1./12.*armHHbar[36] + armHHbar[166] + armHHbar[272] + 
   armHHbar[230] + 27./8.*armHHbar[42] - 3./4.*armHHbar[85] - 1./2.*
   armHHbar[89] - 3./2.*armHHbar[87] - 9./16.*armHHbar[83] + 3./8.*
   armHHbar[129] - armHHbar[92];
   armHHbar[186]= - 1 - armHHbar[38];
   armHHbar[195]=armHHbar[186] + armHHbar[135];
   armHHbar[197]= - 1./8.*armHHbar[11];
   armHHbar[195]=armHHbar[197] + armHHbar[178] + 1./3.*armHHbar[195] + 
   armHHbar[214];
   armHHbar[195]=armHHbar[11]*armHHbar[195];
   armHHbar[201]=armHHbar[214] + armHHbar[318] - 1./3. + armHHbar[139];
   armHHbar[201]=armHHbar[12]*armHHbar[201];
   armHHbar[208]=63./8.*armHHbar[7];
   armHHbar[209]=armHHbar[208] + 5 + armHHbar[36];
   armHHbar[209]=armHHbar[7]*armHHbar[209];
   armHHbar[216]=armHHbar[243] + 27./4.*armHHbar[73] + armHHbar[75];
   armHHbar[216]=1./4.*armHHbar[76] + 3*armHHbar[216] + 1./2.*
   armHHbar[78];
   armHHbar[216]=1./4.*MMH*armHHbar[216];
   armHHbar[195]=armHHbar[216] + 1./4.*armHHbar[195] + 1./2.*
   armHHbar[201] + armHHbar[129] + 3./4.*armHHbar[209];
   armHHbar[195]=MMH*armHHbar[195];
   armHHbar[201]= - 71./16. + armHHbar[36];
   armHHbar[201]=1./2.*armHHbar[201] + armHHbar[319];
   armHHbar[201]=armHHbar[17]*armHHbar[201];
   armHHbar[209]= - 33./2.*armHHbar[44] - 9./2.*armHHbar[70] + 11./2.*
   armHHbar[93] - 9*armHHbar[71] + 9./2.*armHHbar[69] + 11*armHHbar[94]
   ;
   armHHbar[201]=1./4.*armHHbar[146] + armHHbar[201] - 9./2.*
   armHHbar[20] + 1./8.*armHHbar[209] - 9*armHHbar[21];
   armHHbar[226]= - 27./16.*armHHbar[12] - 9./32.*armHHbar[7] + 1./6.*
   armHHbar[36] + 4 + armHHbar[213];
   armHHbar[226]=armHHbar[19]*armHHbar[226];
   armHHbar[230]= - 73./48.*armHHbar[12] + 69./32.*armHHbar[7] + 5./24.
   *armHHbar[36] + 2 + 1./8.*armHHbar[38];
   armHHbar[230]=armHHbar[18]*armHHbar[230];
   armHHbar[231]=37./6.*armHHbar[18] + armHHbar[17] - 191./6.*
   armHHbar[19];
   armHHbar[231]=armHHbar[11]*armHHbar[231];
   armHHbar[195]=1./2.*armHHbar[195] + 1./16.*armHHbar[231] + 
   armHHbar[230] + 1./2.*armHHbar[201] + armHHbar[226];
   armHHbar[195]=MMH*armHHbar[195];
   armHHbar[201]=11./4.*armHHbar[17] - armHHbar[19];
   armHHbar[201]=armHHbar[19]*armHHbar[201];
   armHHbar[226]=23./6.*armHHbar[18] + 15*armHHbar[17] + 49./6.*
   armHHbar[19];
   armHHbar[226]=armHHbar[18]*armHHbar[226];
   armHHbar[230]=pow(armHHbar[17],2);
   armHHbar[201]=1./2.*armHHbar[226] + 165./16.*armHHbar[230] + 
   armHHbar[201];
   armHHbar[195]=1./4.*armHHbar[201] + armHHbar[195];
   armHHbar[195]=armHHbar[68]*armHHbar[195];
   armHHbar[113]=armHHbar[114] + armHHbar[195] + armHHbar[113];
   armHHbar[113]=armHHbar[256]*armHHbar[113];
   armHHbar[114]= - 1 + armHHbar[213];
   armHHbar[182]=armHHbar[114] + armHHbar[182];
   armHHbar[178]=armHHbar[197] + armHHbar[178] + 1./3.*armHHbar[182] + 
   armHHbar[214];
   armHHbar[178]=armHHbar[11]*armHHbar[178];
   armHHbar[182]=armHHbar[208] + armHHbar[135] + 5 + armHHbar[213];
   armHHbar[182]=armHHbar[7]*armHHbar[182];
   armHHbar[135]=armHHbar[135] - 1 + armHHbar[139];
   armHHbar[135]=1./3.*armHHbar[135] + armHHbar[214];
   armHHbar[135]=armHHbar[12]*armHHbar[135];
   armHHbar[135]=armHHbar[216] + 1./4.*armHHbar[178] + 1./2.*
   armHHbar[135] + armHHbar[129] + 3./4.*armHHbar[182];
   armHHbar[135]=MMH*armHHbar[135];
   armHHbar[178]=armHHbar[36] - 71./8. + armHHbar[38];
   armHHbar[178]=1./4.*armHHbar[178] + armHHbar[319];
   armHHbar[178]=armHHbar[17]*armHHbar[178];
   armHHbar[182]=3./2. + armHHbar[316];
   armHHbar[182]= - 11./2.*armHHbar[12] + 15./8.*armHHbar[7] + 11./2.*
   armHHbar[182] + armHHbar[318];
   armHHbar[182]=armHHbar[19]*armHHbar[182];
   armHHbar[195]=armHHbar[20] + armHHbar[209] - 41*armHHbar[21];
   armHHbar[197]=11./3.*armHHbar[12] + 21./4.*armHHbar[7] + 11./6.*
   armHHbar[36] - 5./2. + armHHbar[316];
   armHHbar[197]=armHHbar[18]*armHHbar[197];
   armHHbar[201]=23./8.*armHHbar[18] + armHHbar[17] - 19./8.*
   armHHbar[19];
   armHHbar[201]=armHHbar[11]*armHHbar[201];
   armHHbar[135]=armHHbar[135] + 1./6.*armHHbar[201] + 1./4.*
   armHHbar[197] + 1./2.*armHHbar[182] + 1./3.*armHHbar[146] + 1./8.*
   armHHbar[195] + armHHbar[178];
   armHHbar[135]=MMH*armHHbar[135];
   armHHbar[178]=17*armHHbar[17] - 65./3.*armHHbar[19];
   armHHbar[178]=armHHbar[19]*armHHbar[178];
   armHHbar[178]=165./4.*armHHbar[230] + armHHbar[178];
   armHHbar[182]=13./2.*armHHbar[17];
   armHHbar[195]=armHHbar[182] + 53./3.*armHHbar[19];
   armHHbar[197]=1./3.*armHHbar[18];
   armHHbar[195]=1./2.*armHHbar[195] + armHHbar[197];
   armHHbar[195]=armHHbar[18]*armHHbar[195];
   armHHbar[178]=1./4.*armHHbar[178] + armHHbar[195];
   armHHbar[135]=1./2.*armHHbar[178] + armHHbar[135];
   armHHbar[135]=armHHbar[68]*armHHbar[135];
   armHHbar[103]=armHHbar[113] + armHHbar[135] + armHHbar[103];
   armHHbar[103]=armHHbar[256]*armHHbar[103];
   armHHbar[113]= - armHHbar[7]*armHHbar[38];
   armHHbar[178]= - armHHbar[12]*armHHbar[38];
   armHHbar[195]= - armHHbar[11]*armHHbar[38];
   armHHbar[201]=1./24.*armHHbar[195] + 1./12.*armHHbar[178] + 3./8.*
   armHHbar[113] - 1./4.*armHHbar[15] - armHHbar[16] + 1./4.*
   armHHbar[85] - 15./16. + armHHbar[87];
   armHHbar[201]=MMH*armHHbar[201];
   armHHbar[213]= - 3./2.*armHHbar[7] + 5 - armHHbar[38];
   armHHbar[213]=armHHbar[17]*armHHbar[213];
   armHHbar[226]= - 11./24.*armHHbar[12] + 3./16.*armHHbar[7] + 3 + 
   armHHbar[287];
   armHHbar[226]=armHHbar[19]*armHHbar[226];
   armHHbar[231]=3 - 1./6.*armHHbar[38];
   armHHbar[233]=armHHbar[18]*armHHbar[231];
   armHHbar[236]=5*armHHbar[17] + armHHbar[19];
   armHHbar[236]=1./6.*armHHbar[236] - armHHbar[18];
   armHHbar[236]=armHHbar[11]*armHHbar[236];
   armHHbar[243]= - armHHbar[94] + 1./2.*armHHbar[71];
   armHHbar[201]=1./2.*armHHbar[201] + 1./8.*armHHbar[236] + 1./4.*
   armHHbar[233] + armHHbar[226] + 11./24.*armHHbar[146] + 1./8.*
   armHHbar[213] - 1./4.*armHHbar[20] - armHHbar[21] + 5./8.*
   armHHbar[44] + 1./8.*armHHbar[70] + armHHbar[243] - 1./4.*
   armHHbar[93];
   armHHbar[201]=MMH*armHHbar[201];
   armHHbar[182]=armHHbar[182] - 5*armHHbar[19];
   armHHbar[182]=armHHbar[19]*armHHbar[182];
   armHHbar[153]=armHHbar[17] + armHHbar[153];
   armHHbar[153]=1./3.*armHHbar[153] - 1./2.*armHHbar[18];
   armHHbar[153]=armHHbar[18]*armHHbar[153];
   armHHbar[153]=armHHbar[153] - 1./2.*armHHbar[230] + 1./3.*
   armHHbar[182];
   armHHbar[153]=1./4.*armHHbar[153] + armHHbar[201];
   armHHbar[153]=armHHbar[68]*MMH*armHHbar[153];
   armHHbar[113]=1./12.*armHHbar[195] + 1./6.*armHHbar[178] + 3./4.*
   armHHbar[113] - armHHbar[16] - 3./4. + armHHbar[87];
   armHHbar[113]=MMH*armHHbar[113];
   armHHbar[178]= - 3./4.*armHHbar[7];
   armHHbar[139]=armHHbar[178] + 1 + armHHbar[139];
   armHHbar[139]=armHHbar[17]*armHHbar[139];
   armHHbar[182]=3./8.*armHHbar[7];
   armHHbar[195]= - 5./12.*armHHbar[12] + armHHbar[231] + armHHbar[182]
   ;
   armHHbar[195]=armHHbar[19]*armHHbar[195];
   armHHbar[201]=5./12.*armHHbar[146];
   armHHbar[213]= - armHHbar[18]*armHHbar[38];
   armHHbar[226]= - armHHbar[17] + armHHbar[19];
   armHHbar[231]=armHHbar[11]*armHHbar[226];
   armHHbar[113]=1./2.*armHHbar[113] + 1./24.*armHHbar[231] + 1./12.*
   armHHbar[213] + armHHbar[195] + armHHbar[201] + 1./2.*armHHbar[139]
    - armHHbar[21] + armHHbar[243] + 1./2.*armHHbar[44];
   armHHbar[113]=MMH*armHHbar[113];
   armHHbar[139]=7./4.*armHHbar[17] - armHHbar[19];
   armHHbar[139]=armHHbar[19]*armHHbar[139];
   armHHbar[195]=armHHbar[18]*armHHbar[226];
   armHHbar[113]=armHHbar[113] + 1./12.*armHHbar[195] - 1./4.*
   armHHbar[230] + 1./3.*armHHbar[139];
   armHHbar[113]=armHHbar[171]*armHHbar[68]*MMH*armHHbar[113];
   armHHbar[113]=armHHbar[153] + 1./2.*armHHbar[113];
   armHHbar[113]=armHHbar[171]*armHHbar[113];
   armHHbar[139]= - 1 + armHHbar[87];
   armHHbar[153]= - armHHbar[7]*armHHbar[36];
   armHHbar[195]= - armHHbar[12]*armHHbar[36];
   armHHbar[213]= - armHHbar[11]*armHHbar[36];
   armHHbar[139]=1./12.*armHHbar[213] + 1./6.*armHHbar[195] + 3./4.*
   armHHbar[153] - armHHbar[15] + armHHbar[222] + 3*armHHbar[139] + 
   armHHbar[85];
   armHHbar[139]=MMH*armHHbar[139];
   armHHbar[131]=armHHbar[131] + 3*armHHbar[243] - armHHbar[93];
   armHHbar[153]= - 1./6.*armHHbar[36];
   armHHbar[195]=armHHbar[142] + 9 + armHHbar[153];
   armHHbar[195]=armHHbar[19]*armHHbar[195];
   armHHbar[182]=1./12.*armHHbar[12] + armHHbar[182] + 3 + 
   armHHbar[288];
   armHHbar[182]=armHHbar[18]*armHHbar[182];
   armHHbar[213]= - 3./16.*armHHbar[7] + 1 - 1./8.*armHHbar[36];
   armHHbar[213]=armHHbar[17]*armHHbar[213];
   armHHbar[222]=armHHbar[17] - armHHbar[18];
   armHHbar[222]=armHHbar[11]*armHHbar[222];
   armHHbar[131]=1./4.*armHHbar[139] + 11./48.*armHHbar[222] + 1./2.*
   armHHbar[182] + 1./2.*armHHbar[195] + 17./24.*armHHbar[146] + 
   armHHbar[213] + armHHbar[189] - 3./2.*armHHbar[21] + 1./2.*
   armHHbar[131] + armHHbar[44];
   armHHbar[131]=MMH*armHHbar[131];
   armHHbar[139]=1./3.*armHHbar[17];
   armHHbar[182]=armHHbar[139] - 3./8.*armHHbar[19];
   armHHbar[182]=armHHbar[19]*armHHbar[182];
   armHHbar[195]= - 5./8.*armHHbar[18] + armHHbar[17] + 1./4.*
   armHHbar[19];
   armHHbar[195]=armHHbar[18]*armHHbar[195];
   armHHbar[131]=1./2.*armHHbar[131] + 1./6.*armHHbar[195] - 1./16.*
   armHHbar[230] + armHHbar[182];
   armHHbar[131]=armHHbar[68]*MMH*armHHbar[131];
   armHHbar[113]=armHHbar[131] + 1./2.*armHHbar[113];
   armHHbar[113]=armHHbar[171]*armHHbar[113];
   armHHbar[182]=armHHbar[7]*armHHbar[137];
   armHHbar[195]=armHHbar[12]*armHHbar[137];
   armHHbar[213]=armHHbar[11]*armHHbar[137];
   armHHbar[170]=1./12.*armHHbar[213] + 1./6.*armHHbar[195] + 3./4.*
   armHHbar[182] + armHHbar[170] - armHHbar[16] + armHHbar[241] - 9./8.
    + armHHbar[87];
   armHHbar[170]=MMH*armHHbar[170];
   armHHbar[153]=armHHbar[153] + 3 + armHHbar[166];
   armHHbar[166]= - 7./12.*armHHbar[12] + armHHbar[153] - 3./8.*
   armHHbar[7];
   armHHbar[166]=armHHbar[19]*armHHbar[166];
   armHHbar[222]=3./4.*armHHbar[7];
   armHHbar[153]=1./6.*armHHbar[12] + armHHbar[153] + armHHbar[222];
   armHHbar[153]=armHHbar[18]*armHHbar[153];
   armHHbar[231]=3 + armHHbar[38];
   armHHbar[233]=armHHbar[231] - armHHbar[36];
   armHHbar[233]=armHHbar[17]*armHHbar[233];
   armHHbar[236]= - 5./6.*armHHbar[18] + armHHbar[17] - 1./6.*
   armHHbar[19];
   armHHbar[236]=armHHbar[11]*armHHbar[236];
   armHHbar[111]=1./2.*armHHbar[170] + 1./4.*armHHbar[236] + 1./2.*
   armHHbar[153] + armHHbar[166] + 1./2.*armHHbar[146] + 1./4.*
   armHHbar[233] + armHHbar[189] - armHHbar[21] + 3./4.*armHHbar[44] + 
   1./4.*armHHbar[70] + armHHbar[243] + armHHbar[111];
   armHHbar[111]=MMH*armHHbar[111];
   armHHbar[146]=1./8.*armHHbar[17] + armHHbar[175];
   armHHbar[146]=armHHbar[19]*armHHbar[146];
   armHHbar[153]= - 1./3.*armHHbar[18];
   armHHbar[166]=armHHbar[153] + armHHbar[17] + 1./6.*armHHbar[19];
   armHHbar[166]=armHHbar[18]*armHHbar[166];
   armHHbar[111]=1./2.*armHHbar[111] + armHHbar[146] + 1./4.*
   armHHbar[166];
   armHHbar[111]=armHHbar[68]*MMH*armHHbar[111];
   armHHbar[146]=armHHbar[16] - armHHbar[64] - 1./4. - armHHbar[56];
   armHHbar[166]=armHHbar[146] + 3./4.*armHHbar[308];
   armHHbar[166]=3*armHHbar[166] + 1./2.*armHHbar[285];
   armHHbar[170]=1./8.*armHHbar[290];
   armHHbar[166]=armHHbar[170] + 1./2.*armHHbar[166] + armHHbar[191];
   armHHbar[166]=MMH*armHHbar[166];
   armHHbar[189]=armHHbar[302] + 1./2.*armHHbar[21];
   armHHbar[233]= - armHHbar[17]*armHHbar[30];
   armHHbar[236]=armHHbar[189] + 1./4.*armHHbar[233];
   armHHbar[241]= - 9 - armHHbar[30];
   armHHbar[241]=armHHbar[19]*armHHbar[241];
   armHHbar[243]=MMt*armHHbar[293];
   armHHbar[252]=3./2.*armHHbar[243];
   armHHbar[261]= - 1./4.*armHHbar[30] - armHHbar[27];
   armHHbar[261]=armHHbar[18]*armHHbar[261];
   armHHbar[166]=armHHbar[252] + armHHbar[166] + armHHbar[261] + 
   armHHbar[257] + armHHbar[328] + 3*armHHbar[236] + 1./2.*
   armHHbar[241];
   armHHbar[166]=MMt*armHHbar[166];
   armHHbar[236]= - 9./8.*armHHbar[7];
   armHHbar[147]=armHHbar[147] + armHHbar[236] + armHHbar[36] + 3./2.
    - armHHbar[38];
   armHHbar[147]=armHHbar[28]*armHHbar[147];
   armHHbar[241]= - 3./2.*armHHbar[45];
   armHHbar[262]=1./4.*armHHbar[301];
   armHHbar[263]=1./4.*armHHbar[18]*armHHbar[27];
   armHHbar[269]=1./8.*armHHbar[251];
   armHHbar[150]=1./4.*MMH*armHHbar[150];
   armHHbar[147]=armHHbar[150] + armHHbar[269] + armHHbar[263] + 
   armHHbar[262] + armHHbar[241] + armHHbar[147];
   armHHbar[147]=MMH*armHHbar[147];
   armHHbar[272]= - 3./4.*armHHbar[17] - armHHbar[19];
   armHHbar[272]=armHHbar[28]*armHHbar[272];
   armHHbar[147]=armHHbar[166] + armHHbar[147] + armHHbar[272] + 5./4.*
   armHHbar[266];
   armHHbar[147]=MMt*armHHbar[147];
   armHHbar[166]=armHHbar[28]*armHHbar[30];
   armHHbar[166]=armHHbar[166] + armHHbar[218];
   armHHbar[166]=armHHbar[3]*armHHbar[166];
   armHHbar[272]= - armHHbar[27]*armHHbar[30];
   armHHbar[274]=armHHbar[27]*armHHbar[30];
   armHHbar[278]=MMt*armHHbar[3]*armHHbar[274];
   armHHbar[166]=armHHbar[278] + 1./4.*armHHbar[272] + armHHbar[166];
   armHHbar[166]=MMt*armHHbar[166];
   armHHbar[166]=armHHbar[166] + 1./4.*armHHbar[265] + armHHbar[270];
   armHHbar[166]=3*armHHbar[26]*armHHbar[181]*armHHbar[166];
   armHHbar[147]=1./2.*armHHbar[147] + armHHbar[166];
   armHHbar[147]=armHHbar[26]*armHHbar[147];
   armHHbar[182]=1./6.*armHHbar[213] + 3./2.*armHHbar[182] + 1./3.*
   armHHbar[195];
   armHHbar[182]=MMH*armHHbar[182];
   armHHbar[195]= - 1./6.*armHHbar[12] + armHHbar[185] + armHHbar[178];
   armHHbar[195]=armHHbar[19]*armHHbar[195];
   armHHbar[137]=armHHbar[17]*armHHbar[137];
   armHHbar[185]=1./3.*armHHbar[12] + armHHbar[185] + 3./2.*armHHbar[7]
   ;
   armHHbar[185]=armHHbar[18]*armHHbar[185];
   armHHbar[137]=1./2.*armHHbar[182] + 1./12.*armHHbar[193] + 1./2.*
   armHHbar[185] + 1./2.*armHHbar[137] + armHHbar[195];
   armHHbar[137]=MMH*armHHbar[137];
   armHHbar[128]=armHHbar[19]*armHHbar[128];
   armHHbar[182]=armHHbar[197] + armHHbar[17] + armHHbar[291];
   armHHbar[182]=armHHbar[18]*armHHbar[182];
   armHHbar[128]=armHHbar[137] + armHHbar[128] + 1./2.*armHHbar[182];
   armHHbar[128]=armHHbar[68]*MMH*armHHbar[128];
   armHHbar[126]=armHHbar[126] + armHHbar[172] + armHHbar[236];
   armHHbar[126]=armHHbar[28]*armHHbar[126];
   armHHbar[126]=armHHbar[150] + armHHbar[269] + armHHbar[263] + 
   armHHbar[126] + armHHbar[262];
   armHHbar[126]=MMH*armHHbar[126];
   armHHbar[137]=9./2.*armHHbar[308] + armHHbar[285];
   armHHbar[137]=armHHbar[170] + 1./4.*armHHbar[137] + armHHbar[191];
   armHHbar[137]=MMH*armHHbar[137];
   armHHbar[115]=3./2.*armHHbar[233] + armHHbar[115];
   armHHbar[115]=armHHbar[137] + armHHbar[261] + 1./2.*armHHbar[115] + 
   armHHbar[257];
   armHHbar[115]=MMt*armHHbar[115];
   armHHbar[137]=armHHbar[292] + armHHbar[19];
   armHHbar[137]=armHHbar[28]*armHHbar[137];
   armHHbar[137]=armHHbar[137] + 5./2.*armHHbar[266];
   armHHbar[115]=armHHbar[115] + 1./2.*armHHbar[137] + armHHbar[126];
   armHHbar[115]=MMt*armHHbar[115];
   armHHbar[115]=1./2.*armHHbar[115] + armHHbar[166];
   armHHbar[115]=armHHbar[26]*armHHbar[115];
   armHHbar[115]=1./8.*armHHbar[128] + armHHbar[115];
   armHHbar[115]=armHHbar[256]*armHHbar[115];
   armHHbar[111]=armHHbar[115] + 1./2.*armHHbar[111] + armHHbar[147];
   armHHbar[111]=armHHbar[256]*armHHbar[111];
   armHHbar[115]=armHHbar[27]*armHHbar[36];
   armHHbar[115]=9./2.*armHHbar[146] + armHHbar[115];
   armHHbar[115]=MMH*armHHbar[115];
   armHHbar[126]= - armHHbar[28] + armHHbar[189] - 3./2.*armHHbar[19];
   armHHbar[128]=armHHbar[27]*armHHbar[17];
   armHHbar[115]=9./2.*armHHbar[243] + armHHbar[115] + armHHbar[229] + 
   9*armHHbar[126] + armHHbar[128];
   armHHbar[115]=MMt*armHHbar[115];
   armHHbar[126]=armHHbar[305] + 9./2. + armHHbar[36];
   armHHbar[126]=armHHbar[28]*armHHbar[126];
   armHHbar[128]= - MMH*armHHbar[27]*armHHbar[36];
   armHHbar[126]=1./4.*armHHbar[128] + armHHbar[263] + 1./4.*
   armHHbar[289] - 9./2.*armHHbar[45] + armHHbar[126];
   armHHbar[126]=MMH*armHHbar[126];
   armHHbar[128]=armHHbar[17] - 9./2.*armHHbar[19];
   armHHbar[128]=armHHbar[28]*armHHbar[128];
   armHHbar[115]=armHHbar[115] + armHHbar[126] + armHHbar[128] + 
   armHHbar[266];
   armHHbar[115]=MMt*armHHbar[115];
   armHHbar[126]=armHHbar[26]*armHHbar[115];
   armHHbar[111]=armHHbar[111] + armHHbar[131] + 1./2.*armHHbar[126];
   armHHbar[111]=armHHbar[256]*armHHbar[111];
   armHHbar[126]=armHHbar[146] + 3./8.*armHHbar[119];
   armHHbar[128]=armHHbar[27]*armHHbar[38];
   armHHbar[131]=1./8.*armHHbar[125];
   armHHbar[126]=armHHbar[131] + armHHbar[128] + 3*armHHbar[126] + 1./4.
   *armHHbar[306];
   armHHbar[126]=MMH*armHHbar[126];
   armHHbar[137]=armHHbar[17]*armHHbar[30];
   armHHbar[147]=armHHbar[189] + 1./8.*armHHbar[137];
   armHHbar[150]= - 9 + armHHbar[242];
   armHHbar[150]=armHHbar[19]*armHHbar[150];
   armHHbar[166]=armHHbar[17] - armHHbar[19];
   armHHbar[166]=armHHbar[27]*armHHbar[166];
   armHHbar[170]=armHHbar[18]*armHHbar[30];
   armHHbar[126]=armHHbar[252] + 1./2.*armHHbar[126] + 1./8.*
   armHHbar[170] + 1./2.*armHHbar[166] + armHHbar[328] + 3*
   armHHbar[147] + 1./2.*armHHbar[150];
   armHHbar[126]=MMt*armHHbar[126];
   armHHbar[147]=9./8.*armHHbar[7];
   armHHbar[150]= - 11./4.*armHHbar[12] + armHHbar[231] + armHHbar[147]
   ;
   armHHbar[150]=armHHbar[28]*armHHbar[150];
   armHHbar[172]=1./4.*armHHbar[27]*armHHbar[226];
   armHHbar[182]=1./8.*armHHbar[217];
   armHHbar[185]= - 1./4.*MMH*armHHbar[27]*armHHbar[38];
   armHHbar[150]=armHHbar[185] + armHHbar[182] + armHHbar[172] - 3*
   armHHbar[45] + armHHbar[150];
   armHHbar[150]=MMH*armHHbar[150];
   armHHbar[173]=armHHbar[173] - armHHbar[19];
   armHHbar[173]=armHHbar[28]*armHHbar[173];
   armHHbar[173]=7*armHHbar[173] + 1./2.*armHHbar[322];
   armHHbar[150]=1./2.*armHHbar[173] + armHHbar[150];
   armHHbar[126]=1./2.*armHHbar[150] + armHHbar[126];
   armHHbar[126]=MMt*armHHbar[126];
   armHHbar[119]=armHHbar[146] + 3./4.*armHHbar[119];
   armHHbar[119]=3*armHHbar[119] + armHHbar[317];
   armHHbar[119]=armHHbar[131] + 1./2.*armHHbar[119] + armHHbar[128];
   armHHbar[119]=MMH*armHHbar[119];
   armHHbar[128]=armHHbar[189] + 1./4.*armHHbar[137];
   armHHbar[131]= - 9 + armHHbar[30];
   armHHbar[131]=armHHbar[19]*armHHbar[131];
   armHHbar[119]=armHHbar[252] + armHHbar[119] + 1./4.*armHHbar[170] + 
   armHHbar[166] + armHHbar[328] + 3*armHHbar[128] + 1./2.*
   armHHbar[131];
   armHHbar[119]=MMt*armHHbar[119];
   armHHbar[128]= - 5./4.*armHHbar[12] + armHHbar[147] + 3./2. + 
   armHHbar[38];
   armHHbar[128]=armHHbar[28]*armHHbar[128];
   armHHbar[128]=armHHbar[185] + armHHbar[182] + armHHbar[172] + 
   armHHbar[241] + armHHbar[128];
   armHHbar[128]=MMH*armHHbar[128];
   armHHbar[131]=7./8.*armHHbar[17] - armHHbar[19];
   armHHbar[131]=armHHbar[28]*armHHbar[131];
   armHHbar[119]=1./2.*armHHbar[119] + 1./2.*armHHbar[128] + 
   armHHbar[131] + 1./8.*armHHbar[322];
   armHHbar[119]=armHHbar[171]*MMt*armHHbar[119];
   armHHbar[119]=armHHbar[126] + armHHbar[119];
   armHHbar[119]=armHHbar[171]*armHHbar[119];
   armHHbar[115]=1./2.*armHHbar[115] + armHHbar[119];
   armHHbar[115]=armHHbar[171]*armHHbar[115];
   armHHbar[119]= - armHHbar[28]*armHHbar[30];
   armHHbar[119]=armHHbar[119] + armHHbar[265];
   armHHbar[119]=armHHbar[3]*armHHbar[119];
   armHHbar[126]=MMt*armHHbar[3]*armHHbar[272];
   armHHbar[119]=armHHbar[126] + 1./4.*armHHbar[274] + armHHbar[119];
   armHHbar[119]=MMt*armHHbar[119];
   armHHbar[119]=armHHbar[119] + 1./4.*armHHbar[218] + armHHbar[127];
   armHHbar[119]=armHHbar[181]*armHHbar[119];
   armHHbar[126]=armHHbar[171]*armHHbar[119];
   armHHbar[119]=armHHbar[119] + armHHbar[126];
   armHHbar[119]=armHHbar[26]*armHHbar[119]*pow(armHHbar[2],4);
   armHHbar[115]=armHHbar[115] + 3*armHHbar[119];
   armHHbar[115]=armHHbar[26]*armHHbar[115];
   armHHbar[111]=armHHbar[111] + armHHbar[113] + armHHbar[115];
   armHHbar[111]=armHHbar[8]*armHHbar[111];
   armHHbar[113]= - 1./2. + armHHbar[38];
   armHHbar[113]= - 1./16.*armHHbar[11] + armHHbar[136] + armHHbar[240]
    + 1./3.*armHHbar[113] - 1./4.*armHHbar[36];
   armHHbar[113]=armHHbar[11]*armHHbar[113];
   armHHbar[114]=1./3.*armHHbar[114] + armHHbar[214];
   armHHbar[114]=armHHbar[12]*armHHbar[114];
   armHHbar[115]=armHHbar[208] + 5 + armHHbar[38];
   armHHbar[115]=armHHbar[7]*armHHbar[115];
   armHHbar[113]=armHHbar[216] + 1./2.*armHHbar[113] + 1./2.*
   armHHbar[114] + armHHbar[129] + 3./4.*armHHbar[115];
   armHHbar[113]=MMH*armHHbar[113];
   armHHbar[114]= - 71./16. + armHHbar[38];
   armHHbar[114]=1./2.*armHHbar[114] + armHHbar[319];
   armHHbar[114]=armHHbar[17]*armHHbar[114];
   armHHbar[115]=11*armHHbar[20] + 1./2.*armHHbar[209] - 17*
   armHHbar[21];
   armHHbar[119]= - 33./4.*armHHbar[12] - 105./8.*armHHbar[7] + 13./2.
    + armHHbar[156];
   armHHbar[119]=armHHbar[19]*armHHbar[119];
   armHHbar[126]= - 13 - 1./3.*armHHbar[38];
   armHHbar[126]=23./4.*armHHbar[12] - 27./8.*armHHbar[7] + 1./2.*
   armHHbar[126] + armHHbar[36];
   armHHbar[126]=armHHbar[18]*armHHbar[126];
   armHHbar[127]=3./2.*armHHbar[18] + 5./3.*armHHbar[17] - 13./2.*
   armHHbar[19];
   armHHbar[127]=armHHbar[11]*armHHbar[127];
   armHHbar[113]=armHHbar[113] + 1./8.*armHHbar[127] + 1./2.*
   armHHbar[126] + 1./2.*armHHbar[119] + armHHbar[201] + 1./4.*
   armHHbar[115] + armHHbar[114];
   armHHbar[113]=MMH*armHHbar[113];
   armHHbar[114]= - 25./2.*armHHbar[17] - 71./3.*armHHbar[19];
   armHHbar[114]=armHHbar[19]*armHHbar[114];
   armHHbar[114]=165./8.*armHHbar[230] + armHHbar[114];
   armHHbar[115]= - 5./4.*armHHbar[18] - armHHbar[17] + 43./12.*
   armHHbar[19];
   armHHbar[115]=armHHbar[18]*armHHbar[115];
   armHHbar[114]=1./2.*armHHbar[114] + armHHbar[115];
   armHHbar[113]=1./2.*armHHbar[114] + armHHbar[113];
   armHHbar[113]=armHHbar[68]*armHHbar[113];
   armHHbar[114]=5./6.*armHHbar[12] + armHHbar[186] + armHHbar[178];
   armHHbar[114]=armHHbar[18]*armHHbar[114];
   armHHbar[115]= - 5./6.*armHHbar[12] + 1 + armHHbar[222];
   armHHbar[115]=armHHbar[19]*armHHbar[115];
   armHHbar[119]= - 1./12.*armHHbar[18] + armHHbar[139] - 1./4.*
   armHHbar[19];
   armHHbar[119]=armHHbar[11]*armHHbar[119];
   armHHbar[126]=MMH*armHHbar[11]*armHHbar[38];
   armHHbar[114]=1./3.*armHHbar[126] + armHHbar[119] + armHHbar[114] + 
   armHHbar[115] - armHHbar[21] + armHHbar[20];
   armHHbar[114]=MMH*armHHbar[114];
   armHHbar[115]=armHHbar[153] + armHHbar[260] + armHHbar[174];
   armHHbar[115]=armHHbar[18]*armHHbar[115];
   armHHbar[119]=1./4.*armHHbar[17] + armHHbar[175];
   armHHbar[119]=armHHbar[19]*armHHbar[119];
   armHHbar[114]=1./2.*armHHbar[114] + armHHbar[119] + 1./4.*
   armHHbar[115];
   armHHbar[114]=armHHbar[171]*armHHbar[68]*armHHbar[114];
   armHHbar[113]=armHHbar[113] + 1./2.*armHHbar[114];
   armHHbar[113]=armHHbar[171]*armHHbar[113];
   armHHbar[113]=armHHbar[135] + 1./2.*armHHbar[113];
   armHHbar[113]=armHHbar[171]*armHHbar[113];
   armHHbar[103]=armHHbar[111] + armHHbar[103] + armHHbar[113] + 
   armHHbar[161];
   armHHbar[103]=armHHbar[8]*armHHbar[103];
   armHHbar[111]=armHHbar[40] + armHHbar[130] - 413./3. + armHHbar[121]
   ;
   armHHbar[111]=armHHbar[220] + armHHbar[180] + armHHbar[194] + 
   armHHbar[183] + armHHbar[177] + armHHbar[304] + 1./2.*armHHbar[111]
    + armHHbar[163];
   armHHbar[111]=armHHbar[27]*armHHbar[111];
   armHHbar[104]=armHHbar[167] + armHHbar[104] + armHHbar[111];
   armHHbar[104]=armHHbar[3]*armHHbar[104];
   armHHbar[104]=armHHbar[187] + armHHbar[168] + armHHbar[104];
   armHHbar[104]=MMt*armHHbar[104];
   armHHbar[111]= - 103./2. - 15*armHHbar[41];
   armHHbar[111]=armHHbar[40] + 1./2.*armHHbar[111] + armHHbar[130];
   armHHbar[111]=armHHbar[305] + armHHbar[194] + armHHbar[183] + 
   armHHbar[177] + armHHbar[304] + 1./2.*armHHbar[111] + armHHbar[163];
   armHHbar[111]=armHHbar[28]*armHHbar[111];
   armHHbar[113]= - 23 - 9./2.*armHHbar[30];
   armHHbar[113]=armHHbar[19]*armHHbar[113];
   armHHbar[111]=armHHbar[247] + armHHbar[244] + armHHbar[111] + 
   armHHbar[224] + 1./2.*armHHbar[113];
   armHHbar[111]=armHHbar[3]*armHHbar[111];
   armHHbar[113]=397./3. + armHHbar[237];
   armHHbar[113]=armHHbar[152] + 1./2.*armHHbar[113] + armHHbar[248];
   armHHbar[109]=armHHbar[109] + armHHbar[154] + 1./2.*armHHbar[113] + 
   armHHbar[145];
   armHHbar[109]=armHHbar[238] + armHHbar[142] + armHHbar[184] + 1./2.*
   armHHbar[109] + armHHbar[179];
   armHHbar[109]=armHHbar[27]*armHHbar[109];
   armHHbar[113]= - 2 + armHHbar[309];
   armHHbar[113]=armHHbar[12]*armHHbar[113];
   armHHbar[104]=armHHbar[104] + armHHbar[111] + armHHbar[255] + 
   armHHbar[249] + armHHbar[221] + armHHbar[109] + armHHbar[254] + 
   armHHbar[315] + armHHbar[219] + armHHbar[113];
   armHHbar[104]=MMt*armHHbar[104];
   armHHbar[109]=5./2. - 7*armHHbar[37];
   armHHbar[109]=1./3.*armHHbar[109] + armHHbar[144];
   armHHbar[109]=armHHbar[27]*armHHbar[109];
   armHHbar[111]= - 19 + armHHbar[321];
   armHHbar[111]=armHHbar[312] + 1./3.*armHHbar[111] - armHHbar[30];
   armHHbar[111]=armHHbar[11]*armHHbar[111];
   armHHbar[113]= - 25*armHHbar[55] - 97./2.*armHHbar[43] - 73 + 97./2.
   *armHHbar[58];
   armHHbar[114]= - armHHbar[67] + armHHbar[202];
   armHHbar[115]=armHHbar[114] - armHHbar[20];
   armHHbar[115]=armHHbar[99]*armHHbar[115];
   armHHbar[119]=armHHbar[28]*armHHbar[99];
   armHHbar[121]=1 + armHHbar[211];
   armHHbar[121]=armHHbar[18]*armHHbar[99]*armHHbar[121];
   armHHbar[108]=armHHbar[111] + 7./3.*armHHbar[121] + 1./2.*
   armHHbar[109] + 7./3.*armHHbar[119] + 7./6.*armHHbar[115] + 7./9.*
   armHHbar[29] + 25./6.*armHHbar[15] - 151./36.*armHHbar[63] + 1./6.*
   armHHbar[113] + armHHbar[108];
   armHHbar[109]= - 169./2. + 35*armHHbar[37];
   armHHbar[111]=17*armHHbar[39];
   armHHbar[109]=1./6.*armHHbar[109] + armHHbar[111];
   armHHbar[109]=armHHbar[28]*armHHbar[109];
   armHHbar[113]= - 23./3. + 5./4.*armHHbar[29];
   armHHbar[113]=197./12.*armHHbar[27] + 7./3.*armHHbar[113] + 
   armHHbar[212];
   armHHbar[113]=armHHbar[18]*armHHbar[113];
   armHHbar[115]=3./2.*armHHbar[217];
   armHHbar[109]=armHHbar[115] + armHHbar[113] + 61./6.*armHHbar[265]
    + armHHbar[109] + 373./36.*armHHbar[20] + 109./36.*armHHbar[45] + 
   41./9.*armHHbar[67] - 23./8.*armHHbar[33];
   armHHbar[109]=armHHbar[3]*armHHbar[109];
   armHHbar[113]= - 1 + 7*armHHbar[37];
   armHHbar[111]=armHHbar[295] + 5./6.*armHHbar[113] + armHHbar[111];
   armHHbar[111]=armHHbar[27]*armHHbar[111];
   armHHbar[113]=115./6.*armHHbar[63] + 7*armHHbar[60] + 7*armHHbar[43]
    + 7./2.*armHHbar[37] + 241./12. - 7*armHHbar[58];
   armHHbar[119]=3./2.*armHHbar[125];
   armHHbar[111]=armHHbar[119] + 1./3.*armHHbar[113] + armHHbar[111];
   armHHbar[113]=armHHbar[18] + armHHbar[114] + armHHbar[228];
   armHHbar[113]=armHHbar[3]*armHHbar[113];
   armHHbar[111]=1./2.*armHHbar[111] + 7./3.*armHHbar[113];
   armHHbar[111]=armHHbar[3]*armHHbar[111];
   armHHbar[113]= - 1 - armHHbar[63];
   armHHbar[121]=armHHbar[3]*armHHbar[113];
   armHHbar[121]=armHHbar[239] + armHHbar[121];
   armHHbar[121]=MMt*armHHbar[3]*armHHbar[121];
   armHHbar[113]=armHHbar[99]*armHHbar[113];
   armHHbar[111]=7./3.*armHHbar[121] + armHHbar[111] + 7./24.*
   armHHbar[113] - 13./8.*armHHbar[53] + 3*armHHbar[50];
   armHHbar[111]=MMt*armHHbar[111];
   armHHbar[113]=MMH*armHHbar[53];
   armHHbar[108]=armHHbar[111] + 1./2.*armHHbar[109] + 1./4.*
   armHHbar[108] + 1./3.*armHHbar[113];
   armHHbar[108]=MMt*armHHbar[108];
   armHHbar[109]=19./8. - armHHbar[13];
   armHHbar[109]=armHHbar[258] + 5./2.*armHHbar[110] + armHHbar[134] + 
   5./8.*armHHbar[63] + 1./3.*armHHbar[109] + 5./8.*armHHbar[55];
   armHHbar[110]=25*armHHbar[196];
   armHHbar[111]=101./12. + armHHbar[110];
   armHHbar[111]=armHHbar[11]*armHHbar[111];
   armHHbar[109]=5*armHHbar[109] + 1./2.*armHHbar[111];
   armHHbar[109]=MMH*armHHbar[109];
   armHHbar[111]= - 119*armHHbar[264] + 215*armHHbar[322];
   armHHbar[111]=armHHbar[3]*armHHbar[111];
   armHHbar[113]=25*armHHbar[206];
   armHHbar[121]=499./6. + armHHbar[113];
   armHHbar[121]=armHHbar[28]*armHHbar[121];
   armHHbar[125]= - 197./2.*armHHbar[27] + 1283./12. + armHHbar[110];
   armHHbar[125]=armHHbar[18]*armHHbar[125];
   armHHbar[109]=1./6.*armHHbar[111] + 1./2.*armHHbar[109] + 19*
   armHHbar[251] + 1./4.*armHHbar[125] + 61./4.*armHHbar[218] + 1./4.*
   armHHbar[121] - 283./16.*armHHbar[20] - 35./8.*armHHbar[45] - 17./3.
   *armHHbar[67] + 25./16.*armHHbar[33];
   armHHbar[111]=armHHbar[28]*armHHbar[39];
   armHHbar[111]=armHHbar[111] + armHHbar[115];
   armHHbar[111]=armHHbar[3]*armHHbar[111];
   armHHbar[115]=armHHbar[27]*armHHbar[39];
   armHHbar[115]=armHHbar[115] + armHHbar[119];
   armHHbar[115]=MMt*armHHbar[3]*armHHbar[115];
   armHHbar[119]= - armHHbar[27]*armHHbar[39];
   armHHbar[111]=armHHbar[115] + 1./4.*armHHbar[119] + armHHbar[111];
   armHHbar[111]=armHHbar[171]*MMt*armHHbar[111];
   armHHbar[108]=1./2.*armHHbar[111] + 1./12.*armHHbar[109] + 
   armHHbar[108];
   armHHbar[108]=armHHbar[171]*armHHbar[108];
   armHHbar[109]=armHHbar[280] - 185./3.*armHHbar[28];
   armHHbar[109]=armHHbar[28]*armHHbar[109];
   armHHbar[109]=armHHbar[109] + 71./3.*armHHbar[266];
   armHHbar[109]=armHHbar[3]*armHHbar[109];
   armHHbar[111]=armHHbar[259] - 143./48. + armHHbar[12];
   armHHbar[111]=armHHbar[28]*armHHbar[111];
   armHHbar[104]=armHHbar[108] + armHHbar[104] + 1./4.*armHHbar[109] + 
   armHHbar[123] + armHHbar[268] + armHHbar[234] + armHHbar[267] + 
   armHHbar[207] + armHHbar[111];
   armHHbar[104]=armHHbar[171]*armHHbar[104];
   armHHbar[108]=armHHbar[26]*armHHbar[171]*armHHbar[140];
   armHHbar[104]=armHHbar[104] + armHHbar[108];
   armHHbar[104]=armHHbar[26]*armHHbar[104];
   armHHbar[102]=armHHbar[103] + armHHbar[102] + armHHbar[105] + 
   armHHbar[104];
   armHHbar[102]=armHHbar[8]*armHHbar[102];
   armHHbar[103]= - 61./12. - armHHbar[24];
   armHHbar[103]= - 15./8.*armHHbar[43] + 15./8.*armHHbar[58] + 1./2.*
   armHHbar[103] + armHHbar[157];
   armHHbar[103]=15./8.*armHHbar[164] + 15./4.*armHHbar[162] + 75./4.*
   armHHbar[196] + armHHbar[188] + armHHbar[120] + armHHbar[211] + 
   armHHbar[223] + armHHbar[117] + armHHbar[112] + armHHbar[107] + 1./2.
   *armHHbar[103] + armHHbar[106];
   armHHbar[103]=armHHbar[141] + 1./2.*armHHbar[103] + armHHbar[124];
   armHHbar[103]=MMZ*armHHbar[103];
   armHHbar[104]=armHHbar[232] + armHHbar[151] + 1 + armHHbar[138];
   armHHbar[104]=MMZ*armHHbar[104];
   armHHbar[104]=armHHbar[104] + armHHbar[165] + armHHbar[160] + 
   armHHbar[271];
   armHHbar[104]=armHHbar[3]*MMZ*armHHbar[104];
   armHHbar[105]=armHHbar[28]*armHHbar[245];
   armHHbar[105]=15./4.*armHHbar[265] + armHHbar[169] + 5*armHHbar[105]
   ;
   armHHbar[103]=armHHbar[104] + armHHbar[103] + armHHbar[253] + 1./2.*
   armHHbar[105] + armHHbar[122];
   armHHbar[103]=armHHbar[3]*armHHbar[103];
   armHHbar[104]= - 17./2.*armHHbar[60] + 5*armHHbar[132] + 17./3.*
   armHHbar[55];
   armHHbar[105]=17./6.*armHHbar[29];
   armHHbar[104]=armHHbar[133] + armHHbar[105] + 1./2.*armHHbar[104] - 
   11./3.*armHHbar[15];
   armHHbar[104]=armHHbar[99]*armHHbar[104];
   armHHbar[106]= - armHHbar[43] + 5./6. + armHHbar[58];
   armHHbar[106]=1./3.*armHHbar[15] - 7./12.*armHHbar[63] + 
   armHHbar[232] + 1./4.*armHHbar[106] - 1./3.*armHHbar[55];
   armHHbar[106]=armHHbar[97]*armHHbar[106];
   armHHbar[107]= - 5*armHHbar[9];
   armHHbar[108]=armHHbar[148] + armHHbar[107];
   armHHbar[108]=armHHbar[99]*armHHbar[108];
   armHHbar[108]=25./2.*armHHbar[97] + armHHbar[108];
   armHHbar[108]=armHHbar[11]*armHHbar[108];
   armHHbar[104]=1./6.*armHHbar[108] + 25./16.*armHHbar[204] + 25./4.*
   armHHbar[106] + armHHbar[104];
   armHHbar[104]=MMZ*armHHbar[104];
   armHHbar[106]= - 7313./144. + 5*armHHbar[22];
   armHHbar[106]= - 19./6.*armHHbar[15] - 175./24.*armHHbar[63] - 523./
   36.*armHHbar[60] + armHHbar[279] + 35./18.*armHHbar[13] + 111./8.*
   armHHbar[43] + 1./3.*armHHbar[106] - 111./8.*armHHbar[58];
   armHHbar[108]=armHHbar[107] + 73./2. + armHHbar[148];
   armHHbar[108]=armHHbar[99]*armHHbar[108];
   armHHbar[108]= - 175./8.*armHHbar[97] + armHHbar[108];
   armHHbar[108]=1./3.*armHHbar[108] + 25./4.*armHHbar[199];
   armHHbar[108]=armHHbar[18]*armHHbar[108];
   armHHbar[107]=armHHbar[107] - 101./6. + armHHbar[148];
   armHHbar[107]=1./2.*armHHbar[107] + armHHbar[113];
   armHHbar[107]=armHHbar[11]*armHHbar[107];
   armHHbar[109]=5./6.*armHHbar[20] + 1./6.*armHHbar[33] + 
   armHHbar[118];
   armHHbar[109]=armHHbar[97]*armHHbar[109];
   armHHbar[111]=17./2.*armHHbar[33] + 5./2.*armHHbar[5] - 17*
   armHHbar[67];
   armHHbar[111]=1./2.*armHHbar[111];
   armHHbar[112]=armHHbar[111] - 11*armHHbar[20];
   armHHbar[112]=armHHbar[99]*armHHbar[112];
   armHHbar[113]=325./8.*armHHbar[97] + 17*armHHbar[99];
   armHHbar[113]=armHHbar[28]*armHHbar[113];
   armHHbar[115]=25./2.*armHHbar[206] - 25./8. - 17./3.*armHHbar[37];
   armHHbar[115]=armHHbar[27]*armHHbar[115];
   armHHbar[104]=25./48.*armHHbar[200] + armHHbar[104] + 1./3.*
   armHHbar[107] + 1./2.*armHHbar[108] + 1./2.*armHHbar[115] + 1./3.*
   armHHbar[113] + 1./3.*armHHbar[112] + 5./9.*armHHbar[9] + 25./8.*
   armHHbar[109] + 1./2.*armHHbar[106] + 17./9.*armHHbar[29];
   armHHbar[106]=373./12. - armHHbar[24];
   armHHbar[106]=1./2.*armHHbar[106] + armHHbar[157];
   armHHbar[106]=5./8.*armHHbar[43] + 1./3.*armHHbar[106] - 5./8.*
   armHHbar[58];
   armHHbar[106]= - 19./12.*armHHbar[55] + 1./2.*armHHbar[106] + 
   armHHbar[155];
   armHHbar[107]=85*armHHbar[29];
   armHHbar[108]=25*armHHbar[9];
   armHHbar[109]=armHHbar[108] + 101./2. + armHHbar[107];
   armHHbar[109]=1./6.*armHHbar[109] + armHHbar[110];
   armHHbar[109]=armHHbar[11]*armHHbar[109];
   armHHbar[105]=25./4.*armHHbar[198] + armHHbar[109] + 25./8.*
   armHHbar[205] + 25./4.*armHHbar[203] + 125./4.*armHHbar[206] + 
   armHHbar[133] + 25./4.*armHHbar[190] + armHHbar[105] + 145./12.*
   armHHbar[15] + 125./16.*armHHbar[63] + 5*armHHbar[106] + 829./72.*
   armHHbar[60];
   armHHbar[105]=MMZ*armHHbar[105];
   armHHbar[106]=armHHbar[108] - 1111./12. + armHHbar[107];
   armHHbar[106]=1./3.*armHHbar[106] - 25./4.*armHHbar[27];
   armHHbar[106]=armHHbar[18]*armHHbar[106];
   armHHbar[107]=17*armHHbar[37];
   armHHbar[108]= - 395./12. + armHHbar[107];
   armHHbar[108]=armHHbar[28]*armHHbar[108];
   armHHbar[108]=5*armHHbar[108] + 365./8.*armHHbar[20] + 1499./12.*
   armHHbar[45] - 289./4.*armHHbar[33] - 55./4.*armHHbar[5] + 221./3.*
   armHHbar[67];
   armHHbar[105]=armHHbar[105] + 19*armHHbar[217] + 1./2.*armHHbar[106]
    + 1./3.*armHHbar[108] + 25./4.*armHHbar[218];
   armHHbar[106]=29./4.*armHHbar[18] + armHHbar[111] + 17*armHHbar[28];
   armHHbar[108]= - 17./12.*armHHbar[60] + 5./12.*armHHbar[13] - 1 + 5./
   12.*armHHbar[24];
   armHHbar[108]=MMZ*armHHbar[108];
   armHHbar[106]=1./3.*armHHbar[106] + armHHbar[108];
   armHHbar[106]=armHHbar[3]*MMZ*armHHbar[106];
   armHHbar[105]=1./4.*armHHbar[105] + armHHbar[106];
   armHHbar[105]=armHHbar[3]*armHHbar[105];
   armHHbar[106]=armHHbar[227] + 1./2.*armHHbar[114] + armHHbar[28];
   armHHbar[108]= - 25 - 41./3.*armHHbar[60];
   armHHbar[109]= - 17./3.*armHHbar[63];
   armHHbar[108]=1./2.*armHHbar[108] + armHHbar[109];
   armHHbar[108]=MMZ*armHHbar[108];
   armHHbar[106]=41./3.*armHHbar[106] + 1./2.*armHHbar[108];
   armHHbar[106]=armHHbar[3]*armHHbar[106];
   armHHbar[108]=83./3. - 11*armHHbar[58];
   armHHbar[107]=143./4.*armHHbar[43] + 13./4.*armHHbar[108] + 
   armHHbar[107];
   armHHbar[108]=5./4. + 17./3.*armHHbar[37];
   armHHbar[108]=armHHbar[27]*armHHbar[108];
   armHHbar[107]=5*armHHbar[108] + 7./6.*armHHbar[29] + 7./2.*
   armHHbar[15] + 511./12.*armHHbar[63] + 289./12.*armHHbar[60] + 1./3.
   *armHHbar[107] - 7./2.*armHHbar[55];
   armHHbar[108]= - 11 + 35./2.*armHHbar[29];
   armHHbar[108]=1./3.*armHHbar[108] + armHHbar[192];
   armHHbar[108]=1./4.*armHHbar[108] + armHHbar[327];
   armHHbar[108]=armHHbar[11]*armHHbar[108];
   armHHbar[110]=MMZ*armHHbar[50];
   armHHbar[106]=armHHbar[106] + 17./3.*armHHbar[110] + 1./4.*
   armHHbar[107] + armHHbar[108];
   armHHbar[106]=armHHbar[3]*armHHbar[106];
   armHHbar[107]= - 17 + 7*armHHbar[55];
   armHHbar[107]=7./3.*armHHbar[29] - 7./3.*armHHbar[15] + 
   armHHbar[109] + 1./3.*armHHbar[107] - 7./2.*armHHbar[60];
   armHHbar[107]=armHHbar[99]*armHHbar[107];
   armHHbar[108]= - armHHbar[11]*armHHbar[99]*armHHbar[29];
   armHHbar[107]=7./12.*armHHbar[108] - 25./3.*armHHbar[50] + 1./4.*
   armHHbar[107];
   armHHbar[108]=17./2.*armHHbar[53] + 7*armHHbar[50];
   armHHbar[109]= - 41./6.*armHHbar[63] - 8 - 7./6.*armHHbar[60];
   armHHbar[109]=armHHbar[3]*armHHbar[109];
   armHHbar[108]=1./3.*armHHbar[108] + armHHbar[109];
   armHHbar[108]=MMt*armHHbar[3]*armHHbar[108];
   armHHbar[106]=armHHbar[108] + 1./2.*armHHbar[107] + armHHbar[106];
   armHHbar[106]=MMt*armHHbar[106];
   armHHbar[104]=armHHbar[106] + 1./4.*armHHbar[104] + armHHbar[105];
   armHHbar[104]=armHHbar[171]*armHHbar[104];
   armHHbar[105]=1./2.*armHHbar[159] + 1./2.*armHHbar[283] + 
   armHHbar[116] + armHHbar[282];
   armHHbar[106]=armHHbar[225] + 15./8.*armHHbar[199] + 5./2.*
   armHHbar[215] + armHHbar[158];
   armHHbar[106]=MMZ*armHHbar[106];
   armHHbar[107]=armHHbar[275] + 5*armHHbar[277];
   armHHbar[108]=15 + 13*armHHbar[60];
   armHHbar[108]=1./2.*armHHbar[108] + armHHbar[63];
   armHHbar[108]=MMZ*armHHbar[108];
   armHHbar[108]=armHHbar[273] + armHHbar[108];
   armHHbar[108]=armHHbar[3]*armHHbar[108];
   armHHbar[109]= - MMZ*armHHbar[50];
   armHHbar[107]=armHHbar[108] + 2*armHHbar[109] + 1./2.*armHHbar[107]
    + armHHbar[246];
   armHHbar[107]=armHHbar[3]*armHHbar[107];
   armHHbar[107]=armHHbar[276] + armHHbar[235] + armHHbar[107];
   armHHbar[107]=MMt*armHHbar[107];
   armHHbar[103]=armHHbar[104] + armHHbar[107] + armHHbar[103] + 
   armHHbar[149] + 1./4.*armHHbar[106] + 1./2.*armHHbar[105] + 
   armHHbar[281];
   armHHbar[103]=armHHbar[171]*armHHbar[103];
   armHHbar[104]=armHHbar[28]*armHHbar[12];
   armHHbar[105]=MMZ*armHHbar[210];
   armHHbar[104]=8*armHHbar[104] + armHHbar[105];
   armHHbar[105]= - 1 - armHHbar[62];
   armHHbar[106]=pow(MMZ,2);
   armHHbar[105]=armHHbar[3]*armHHbar[106]*armHHbar[105];
   armHHbar[104]=4./3.*armHHbar[104] + 3*armHHbar[105];
   armHHbar[104]=armHHbar[3]*armHHbar[104];
   armHHbar[105]=MMt*armHHbar[3]*armHHbar[12]*armHHbar[176];
   armHHbar[103]=armHHbar[103] + armHHbar[104] + 32./3.*armHHbar[105];
   armHHbar[103]=armHHbar[26]*armHHbar[103];
   armHHbar[104]=3 + armHHbar[39];
   armHHbar[104]=armHHbar[98]*armHHbar[104];
   armHHbar[104]=armHHbar[104] + armHHbar[250];
   armHHbar[104]=armHHbar[12]*armHHbar[104];
   armHHbar[105]=5*armHHbar[79] - 4*armHHbar[81] - 5*armHHbar[80];
   armHHbar[107]= - 2*armHHbar[39] + 4*armHHbar[91] + 1 - 2*
   armHHbar[88];
   armHHbar[107]=armHHbar[98]*armHHbar[107];
   armHHbar[108]=1 + armHHbar[91];
   armHHbar[108]=armHHbar[99]*armHHbar[108];
   armHHbar[104]=6*armHHbar[104] + 3*armHHbar[108] + 2*armHHbar[105] + 
   3*armHHbar[107];
   armHHbar[104]=MMZ*armHHbar[104];
   armHHbar[105]=3 + 2*armHHbar[91];
   armHHbar[107]=1 + armHHbar[40];
   armHHbar[107]=armHHbar[12]*armHHbar[107];
   armHHbar[104]=armHHbar[104] + 3*armHHbar[105] + 4*armHHbar[107];
   armHHbar[104]=armHHbar[68]*armHHbar[104];
   armHHbar[105]=pow(CW,2);
   armHHbar[107]= - 64./3. - 3*armHHbar[105];
   armHHbar[108]= - 12*armHHbar[105];
   armHHbar[109]= - 53 + armHHbar[108];
   armHHbar[109]=armHHbar[40]*armHHbar[109];
   armHHbar[110]= - 6*armHHbar[39];
   armHHbar[107]= - 26./3.*armHHbar[12] + armHHbar[110] + 4*
   armHHbar[107] + armHHbar[109];
   armHHbar[107]=armHHbar[12]*armHHbar[107];
   armHHbar[109]=armHHbar[105]*armHHbar[81];
   armHHbar[109]=armHHbar[81] + armHHbar[109];
   armHHbar[111]=1 + armHHbar[105];
   armHHbar[111]=armHHbar[80]*armHHbar[111];
   armHHbar[109]=2*armHHbar[109] + armHHbar[111];
   armHHbar[109]= - 8*armHHbar[77] + 4*armHHbar[109] + armHHbar[78];
   armHHbar[108]= - 41 + armHHbar[108];
   armHHbar[108]=armHHbar[79]*armHHbar[108];
   armHHbar[108]=3*armHHbar[109] + armHHbar[108];
   armHHbar[108]=MMZ*armHHbar[108];
   armHHbar[109]= - 29./3. + 3*armHHbar[88];
   armHHbar[107]=armHHbar[108] + armHHbar[107] + armHHbar[110] - 6*
   armHHbar[16] - 12*armHHbar[90] + 4*armHHbar[109] - 143./3.*
   armHHbar[91];
   armHHbar[107]=MMZ*armHHbar[107];
   armHHbar[108]= - armHHbar[12] - 1 - armHHbar[40];
   armHHbar[108]=armHHbar[19]*armHHbar[108];
   armHHbar[109]=1 - armHHbar[12];
   armHHbar[109]=armHHbar[18]*armHHbar[109];
   armHHbar[108]=armHHbar[109] - armHHbar[20] + armHHbar[108];
   armHHbar[107]=12*armHHbar[108] + armHHbar[107];
   armHHbar[107]=armHHbar[68]*armHHbar[107];
   armHHbar[108]=armHHbar[14] + 1 + armHHbar[25];
   armHHbar[109]=armHHbar[31]*armHHbar[108];
   armHHbar[108]=armHHbar[1]*armHHbar[108];
   armHHbar[108]=3*armHHbar[109] + armHHbar[108];
   armHHbar[106]=armHHbar[106]*armHHbar[108];
   armHHbar[108]=25 + 6*armHHbar[105];
   armHHbar[105]=47 + 12*armHHbar[105];
   armHHbar[105]=armHHbar[91]*armHHbar[105];
   armHHbar[105]=6*armHHbar[90] + armHHbar[105] + 2*armHHbar[108] - 3*
   armHHbar[92];
   armHHbar[105]=MMZ*armHHbar[105];
   armHHbar[108]= - 2*armHHbar[18] - 4*armHHbar[19] + 2*armHHbar[95] - 
   armHHbar[72];
   armHHbar[105]=6*armHHbar[108] + armHHbar[105];
   armHHbar[105]=armHHbar[68]*MMZ*armHHbar[105];
   armHHbar[105]=armHHbar[106] + 2*armHHbar[105];
   armHHbar[105]=armHHbar[3]*armHHbar[105];
   armHHbar[106]=MMZ*armHHbar[143];
   armHHbar[105]=armHHbar[105] + 4*armHHbar[106] + armHHbar[107];
   armHHbar[105]=armHHbar[3]*armHHbar[105];

      mHHbarret = armHHbar[100] + armHHbar[101] + armHHbar[102] + 
      armHHbar[103] + armHHbar[104] + armHHbar[105];
      return mHHbarret;
}
