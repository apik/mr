#include <HH.hpp>
std::complex<long double> HH::m20(size_t nG)
{     
      
      
    std::complex<long double> mHH[493];

    mHH[1]=pow(MMH,-1);
    mHH[2]=pow(MMZ,-1);
    mHH[3]=pow(MMW,-1);
    mHH[4]=double(nG);
    mHH[5]=Tsil::I2(0,0,MMZ,mu2);
    mHH[6]=Tsil::B(MMH,MMH,MMH,mu2);
    mHH[7]=Tsil::B(MMZ,MMZ,MMH,mu2);
    mHH[8]=Tsil::B(MMW,MMW,MMH,mu2);
    mHH[9]=Tsil::B(MMt,MMt,MMH,mu2);
    mHH[10]=Tsil::B(0,0,MMZ,mu2);
    mHH[11]=Tsil::A(MMH,mu2);
    mHH[12]=Tsil::A(MMZ,mu2);
    mHH[13]=Tsil::A(MMW,mu2);
    mHH[14]=Tsil::A(MMt,mu2);
    mHH[15]=Tsil::B(0,0,MMW,mu2);
    mHH[16]=Tsil::Beps(MMZ,MMZ,MMH,mu2);
    mHH[17]=Tsil::Aeps(MMZ,mu2);
    mHH[18]=protZZ00->Uxzuv(0);
    mHH[19]=protZZ00->Txuv(0);
    mHH[20]=protWW00->Txuv(0);
    mHH[21]=protZZ00->Suxv(0);
    mHH[22]=protWW00->Suxv(0);
    mHH[23]=Tsil::I2(MMH,MMH,MMH,mu2);
    mHH[24]=Tsil::I2(MMH,MMt,MMt,mu2);
    mHH[25]=Tsil::I2(MMZ,MMZ,MMH,mu2);
    mHH[26]=Tsil::I2(MMZ,MMt,MMt,mu2);
    mHH[27]=Tsil::I2(MMW,MMW,MMH,mu2);
    mHH[28]=Tsil::I2(MMW,MMW,MMZ,mu2);
    mHH[29]=Tsil::I2(0,MMW,MMt,mu2);
    mHH[30]=Tsil::B(MMZ,MMH,MMZ,mu2);
    mHH[31]=Tsil::B(MMW,MMH,MMW,mu2);
    mHH[32]=Tsil::B(MMW,MMZ,MMW,mu2);
    mHH[33]=Tsil::B(MMW,MMW,MMZ,mu2);
    mHH[34]=Tsil::B(MMt,MMt,MMZ,mu2);
    mHH[35]=Tsil::B(MMH,MMt,MMt,mu2);
    mHH[36]=Tsil::B(MMZ,MMt,MMt,mu2);
    mHH[37]=Tsil::B(0,MMW,MMt,mu2);
    mHH[38]=Tsil::B(0,MMt,MMW,mu2);
    mHH[39]=Tsil::Beps(MMH,MMH,MMH,mu2);
    mHH[40]=Tsil::Beps(MMW,MMW,MMH,mu2);
    mHH[41]=Tsil::Beps(MMt,MMt,MMH,mu2);
    mHH[42]=Tsil::Aeps(MMH,mu2);
    mHH[43]=Tsil::Aeps(MMW,mu2);
    mHH[44]=Tsil::Aeps(MMt,mu2);
    mHH[45]=prottttt0->M(0);
    mHH[46]=prottttt0->Vxzuv(0);
    mHH[47]=prottttt0->Suxv(0);
    mHH[48]=protHHHHH->M(0);
    mHH[49]=protHZHZZ->M(0);
    mHH[50]=protHWHWW->M(0);
    mHH[51]=protHtHtt->M(0);
    mHH[52]=protZZZZH->M(0);
    mHH[53]=protZWZWW->M(0);
    mHH[54]=protZtZtt->M(0);
    mHH[55]=protWWWWH->M(0);
    mHH[56]=protWWWWZ->M(0);
    mHH[57]=protWWWW0->M(0);
    mHH[58]=protWtWt0->M(0);
    mHH[59]=protttttH->M(0);
    mHH[60]=protttttZ->M(0);
    mHH[61]=protWWWW0->Vxzuv(0);
    mHH[62]=protHHHHH->Uzxyv(0);
    mHH[63]=protHZHZZ->Uzxyv(0);
    mHH[64]=protHWHWW->Uzxyv(0);
    mHH[65]=protHtHtt->Uzxyv(0);
    mHH[66]=protHZHZZ->Uuyxv(0);
    mHH[67]=protZWZWW->Uzxyv(0);
    mHH[68]=protZtZtt->Uzxyv(0);
    mHH[69]=protHWHWW->Uuyxv(0);
    mHH[70]=protZWZWW->Uuyxv(0);
    mHH[71]=protWtWt0->Uzxyv(0);
    mHH[72]=protHtHtt->Uuyxv(0);
    mHH[73]=protZtZtt->Uuyxv(0);
    mHH[74]=protWtWt0->Uuyxv(0);
    mHH[75]=protHZHZZ->Tyzv(0);
    mHH[76]=protZWZWW->Txuv(0);
    mHH[77]=protZtZtt->Txuv(0);
    mHH[78]=protZWZWW->Tuxv(0);
    mHH[79]=protHWHWW->Tuxv(0);
    mHH[80]=protWtWt0->Suxv(0);
    mHH[81]=protWtWt0->Txuv(0);
    mHH[82]=protZtZtt->Tyzv(0);
    mHH[83]=protWtWt0->Tyzv(0);
    mHH[84]=protHtHtt->Tyzv(0);
    mHH[85]=protHtHtt->Svyz(0);
    mHH[86]=protHZHZZ->Svyz(0);
    mHH[87]=protZtZtt->Svyz(0);
    mHH[88]=protHWHWW->Suxv(0);
    mHH[89]=protZWZWW->Suxv(0);
    mHH[90]=protWWWW0->Suxv(0);
    mHH[91]=1/(1 - pow(MMZ,-3)*pow(MMW,3) + 3*pow(MMZ,-2)*pow(MMW,2) - 
   3*pow(MMZ,-1)*MMW);
    mHH[92]=1/(1 + pow(MMZ,-2)*pow(MMW,2) - 2*pow(MMZ,-1)*MMW);
    mHH[93]=Tsil::I2(0,0,MMW,mu2);
    mHH[94]=protWW00->Uxzuv(0);
    mHH[95]=1/(4*MMt - MMZ);
    mHH[96]=1/( - 4*MMW + MMH);
    mHH[97]=1/( - 4*MMZ + MMH);
    mHH[98]=1/(1 - pow(MMZ,-1)*MMW);
   mHH[99]= - 139./2. - 13*mHH[33];
   mHH[100]= - 13*mHH[32];
   mHH[99]=1./2.*mHH[99] + mHH[100];
   mHH[99]=1./3.*mHH[99] - mHH[15];
   mHH[101]= - 13./18.*mHH[10];
   mHH[102]=13./18.*mHH[34];
   mHH[103]=mHH[97]*mHH[11];
   mHH[104]=3./4.*mHH[103];
   mHH[105]=5./3.*mHH[10] - 8./9. + mHH[15];
   mHH[105]=mHH[4]*mHH[105];
   mHH[106]=2./3.*mHH[105];
   mHH[107]=3*mHH[8];
   mHH[108]=mHH[14]*mHH[95];
   mHH[109]=5./2.*mHH[108];
   mHH[110]= - mHH[12]*mHH[97];
   mHH[111]=3./2.*mHH[110];
   mHH[112]=3./8.*mHH[7];
   mHH[113]=27./4.*mHH[6];
   mHH[114]= - 3*mHH[30];
   mHH[115]=1./2.*mHH[38];
   mHH[99]=mHH[112] + mHH[111] + mHH[109] + mHH[107] + mHH[106] + 
   mHH[104] + mHH[115] + mHH[113] + mHH[102] + mHH[114] + mHH[101] + 1./
   2.*mHH[99] + mHH[31];
   mHH[99]=mHH[7]*mHH[99];
   mHH[116]=9*mHH[67];
   mHH[117]=143./2. + mHH[116];
   mHH[118]=1./2.*mHH[18];
   mHH[119]= - 27*mHH[63];
   mHH[117]=mHH[118] + 1./4.*mHH[117] + mHH[119];
   mHH[120]= - 18*mHH[75];
   mHH[121]=3*mHH[66];
   mHH[122]= - 3./8.*mHH[19];
   mHH[123]=27./2.*mHH[39];
   mHH[124]= - 1./4.*mHH[68];
   mHH[125]=3./8.*mHH[77];
   mHH[126]=27./2.*mHH[6];
   mHH[127]=9./8.*mHH[33];
   mHH[128]=3./2.*mHH[30];
   mHH[129]=1./4.*mHH[10];
   mHH[130]= - 1./4.*mHH[34];
   mHH[117]=mHH[126] + mHH[130] + mHH[128] + mHH[129] + mHH[127] - 33./
   8.*mHH[16] + mHH[125] + mHH[124] - 27./16.*mHH[76] + mHH[123] - 1./
   16.*mHH[78] + mHH[122] + mHH[121] + 1./2.*mHH[117] + mHH[120];
   mHH[117]=mHH[97]*mHH[117];
   mHH[131]= - 3./4.*mHH[33];
   mHH[132]= - 1 + mHH[131];
   mHH[133]= - 1./2.*mHH[10];
   mHH[134]=1./2.*mHH[34];
   mHH[132]=mHH[134] + mHH[114] + 3*mHH[132] + mHH[133];
   mHH[132]=mHH[97]*mHH[132];
   mHH[135]= - 5./4.*mHH[95];
   mHH[132]=mHH[135] + mHH[132];
   mHH[132]=mHH[7]*mHH[132];
   mHH[136]= - 3*mHH[73];
   mHH[137]=3*mHH[41] - 5./2. + mHH[136];
   mHH[138]= - 3./2.*mHH[77];
   mHH[137]=7./4.*mHH[82] - mHH[16] + mHH[138] + 1./4.*mHH[137] + 
   mHH[68];
   mHH[137]=5./8.*mHH[95]*mHH[137];
   mHH[139]=3*mHH[49] + 1./4.*mHH[53];
   mHH[140]= - 2*mHH[52];
   mHH[141]= - 19./2. + mHH[70];
   mHH[141]=1./2.*mHH[32] - 1./2.*mHH[40] - 9./2.*mHH[76] + 1./2.*
   mHH[141] - mHH[78];
   mHH[141]=mHH[96]*mHH[141];
   mHH[142]= - mHH[8]*mHH[96]*mHH[32];
   mHH[143]=mHH[9]*mHH[95];
   mHH[144]=15./32.*mHH[143];
   mHH[145]=1./4.*mHH[56];
   mHH[117]=1./2.*mHH[132] + mHH[144] + 1./8.*mHH[142] + mHH[117] + 1./
   4.*mHH[141] + mHH[137] + mHH[140] + 3*mHH[139] + mHH[145];
   mHH[117]=MMZ*mHH[117];
   mHH[132]=9./2.*mHH[28] + mHH[5] - 9*mHH[89];
   mHH[139]= - 1./2.*mHH[21];
   mHH[141]=1./2.*mHH[87];
   mHH[142]= - 1./4.*mHH[26];
   mHH[146]=7./4.*mHH[25];
   mHH[147]=5./4.*mHH[42];
   mHH[148]=5./4.*mHH[11];
   mHH[149]= - 3*mHH[86];
   mHH[132]=mHH[148] - 143./4.*mHH[17] + mHH[142] + mHH[141] + mHH[147]
    + mHH[139] + mHH[146] + 1./4.*mHH[132] + mHH[149];
   mHH[132]=mHH[97]*mHH[132];
   mHH[150]= - 27*mHH[6];
   mHH[151]=mHH[150] + mHH[134] + mHH[114] + mHH[133] + 43 - 9./4.*
   mHH[33];
   mHH[151]=mHH[97]*mHH[151];
   mHH[152]=35./8.*mHH[95];
   mHH[153]=mHH[152] + 9*mHH[96];
   mHH[154]= - mHH[9]*mHH[95];
   mHH[155]=15./8.*mHH[154];
   mHH[151]=mHH[155] + 1./2.*mHH[153] + mHH[151];
   mHH[151]=mHH[12]*mHH[151];
   mHH[153]= - 3*mHH[49] - 1./2.*mHH[53];
   mHH[156]= - mHH[82] + mHH[16] + 1 - mHH[68];
   mHH[156]=mHH[95]*mHH[156];
   mHH[157]=mHH[7]*mHH[95];
   mHH[153]=5./16.*mHH[157] + 5./16.*mHH[156] - mHH[52] + 3*mHH[153] - 
   7./4.*mHH[56];
   mHH[153]=MMH*mHH[153];
   mHH[158]=877./24. - 3*mHH[67];
   mHH[158]=mHH[18] - 17*mHH[70] + 1./2.*mHH[158] + mHH[119];
   mHH[159]= - 13./3.*mHH[75];
   mHH[158]=1./2.*mHH[158] + mHH[159];
   mHH[160]=1 - 5./2.*mHH[32];
   mHH[161]= - 1./4.*mHH[8];
   mHH[160]=1./3.*mHH[160] + mHH[161];
   mHH[160]=mHH[8]*mHH[160];
   mHH[162]= - 1./6.*mHH[19];
   mHH[163]=9./32.*mHH[73];
   mHH[164]= - 9./32.*mHH[41];
   mHH[165]= - 1./48.*mHH[77];
   mHH[166]=35./32.*mHH[82];
   mHH[167]=1./6.*mHH[10];
   mHH[168]= - 1./6.*mHH[34];
   mHH[169]= - 1./2.*mHH[26];
   mHH[170]= - 5./2.*mHH[17] + mHH[169] + 15*mHH[44];
   mHH[170]=5./16.*mHH[95]*mHH[170];
   mHH[171]= - 65./8.*mHH[95] - mHH[97];
   mHH[171]=1./2.*mHH[14]*mHH[171];
   mHH[172]=15./2.*mHH[108] + 15./8. + mHH[36];
   mHH[172]=1./4.*mHH[9]*mHH[172];
   mHH[173]=3./8.*mHH[68];
   mHH[174]=1./2.*mHH[28];
   mHH[175]= - mHH[89] + mHH[174];
   mHH[176]=mHH[175] - mHH[43];
   mHH[176]=mHH[96]*mHH[176];
   mHH[177]=3 - mHH[32];
   mHH[177]=mHH[96]*mHH[177];
   mHH[177]=mHH[177] + mHH[97];
   mHH[177]=mHH[13]*mHH[177];
   mHH[178]=27./4.*mHH[39];
   mHH[179]=3./2.*mHH[32];
   mHH[180]=1./2.*mHH[30];
   mHH[99]=mHH[117] + 1./2.*mHH[153] + mHH[99] + 1./2.*mHH[151] + 
   mHH[172] + mHH[171] + 9./4.*mHH[177] + mHH[160] + 1./2.*mHH[132] + 9.
   /4.*mHH[176] + mHH[113] + mHH[170] + mHH[168] + mHH[180] + mHH[167]
    + mHH[179] + 3./4.*mHH[33] + mHH[166] - 5./4.*mHH[16] + 17./4.*
   mHH[40] + mHH[165] + mHH[173] + mHH[164] + mHH[163] - 9./8.*mHH[76]
    + mHH[178] - 77./24.*mHH[78] + mHH[162] + 1./2.*mHH[158] + mHH[66];
   mHH[99]=MMZ*mHH[99];
   mHH[117]=11*mHH[33];
   mHH[132]=13*mHH[32];
   mHH[151]=3*mHH[15];
   mHH[153]= - 3*mHH[31];
   mHH[158]=mHH[153] + mHH[151] + mHH[132] + 155./4. + mHH[117];
   mHH[160]=8./3.*mHH[10];
   mHH[181]= - 8./3.*mHH[34];
   mHH[182]= - 3./2.*mHH[38];
   mHH[183]= - mHH[14]*mHH[95];
   mHH[184]=15./2.*mHH[183];
   mHH[185]=1./2.*mHH[7];
   mHH[186]= - 5./3.*mHH[10];
   mHH[187]=mHH[186] + 8./9. - mHH[15];
   mHH[188]=mHH[4]*mHH[187];
   mHH[189]=2*mHH[188];
   mHH[190]=6*mHH[30];
   mHH[191]=1./2.*mHH[8];
   mHH[158]=mHH[185] + mHH[184] + mHH[191] + mHH[189] + mHH[182] + 
   mHH[181] + mHH[190] + 1./2.*mHH[158] + mHH[160];
   mHH[158]=mHH[7]*mHH[158];
   mHH[192]= - mHH[41] + 3 + mHH[73];
   mHH[193]=1./2.*mHH[77];
   mHH[192]= - 5./4.*mHH[82] + mHH[16] + mHH[193] + 1./4.*mHH[192] - 
   mHH[68];
   mHH[192]=mHH[95]*mHH[192];
   mHH[194]=3*mHH[52];
   mHH[195]=1./2.*mHH[53];
   mHH[196]=15./8.*mHH[157] + 15./32.*mHH[154] + 15./8.*mHH[192] + 
   mHH[195] + mHH[194];
   mHH[196]=MMZ*mHH[196];
   mHH[197]= - mHH[70] + 85./4. - 49*mHH[67];
   mHH[198]= - 5*mHH[18];
   mHH[197]=1./2.*mHH[197] + mHH[198];
   mHH[199]=8./3.*mHH[75];
   mHH[200]= - 3*mHH[66];
   mHH[201]=13./24.*mHH[19];
   mHH[202]=15./32.*mHH[73];
   mHH[203]= - 15./32.*mHH[41];
   mHH[204]= - 5./8.*mHH[68];
   mHH[205]= - 17./48.*mHH[77];
   mHH[206]= - 75./32.*mHH[82];
   mHH[207]= - 5*mHH[44];
   mHH[208]=mHH[207] + 1./2.*mHH[17];
   mHH[208]=15./8.*mHH[95]*mHH[208];
   mHH[209]=mHH[8]*mHH[32];
   mHH[210]=75./8.*mHH[108];
   mHH[211]= - 1./4. + mHH[183];
   mHH[211]=15./8.*mHH[9]*mHH[211];
   mHH[212]= - mHH[95] + mHH[143];
   mHH[212]=15./16.*mHH[12]*mHH[212];
   mHH[213]=11*mHH[16];
   mHH[127]=mHH[196] + mHH[158] + mHH[212] + mHH[211] + mHH[210] + 3./8.
   *mHH[209] + mHH[208] + mHH[130] + mHH[128] + mHH[129] + 1./8.*
   mHH[32] + mHH[127] + mHH[206] + mHH[213] + 1./8.*mHH[40] + mHH[205]
    + mHH[204] + mHH[203] + mHH[202] + 65./16.*mHH[76] - 7./16.*mHH[78]
    + mHH[201] + mHH[200] + 1./4.*mHH[197] + mHH[199];
   mHH[127]=MMZ*mHH[127];
   mHH[117]=mHH[153] + mHH[151] + mHH[132] + 371./24. + mHH[117];
   mHH[158]=3./4.*mHH[8];
   mHH[196]=15./16.*mHH[9];
   mHH[117]=mHH[196] + mHH[158] + mHH[189] + mHH[182] + mHH[181] + 
   mHH[190] + 1./2.*mHH[117] + mHH[160];
   mHH[117]=mHH[12]*mHH[117];
   mHH[197]=17./3.*mHH[86];
   mHH[209]= - 37./2.*mHH[25];
   mHH[214]= - 11./2.*mHH[5];
   mHH[215]=mHH[209] + mHH[197] - 149./4.*mHH[28] + mHH[214] + 91./3.*
   mHH[89];
   mHH[216]=2./3.*mHH[21];
   mHH[217]=77./24.*mHH[42];
   mHH[218]= - 13./6.*mHH[87];
   mHH[219]=17./8.*mHH[26];
   mHH[220]= - 247./24.*mHH[44];
   mHH[221]=41./24.*mHH[11];
   mHH[222]= - 3./4.*mHH[8] + 9./2. + 7*mHH[32];
   mHH[222]=mHH[13]*mHH[222];
   mHH[223]=55./12. - mHH[36];
   mHH[224]=mHH[14]*mHH[223];
   mHH[225]= - mHH[9]*mHH[14];
   mHH[226]=15./8.*mHH[225];
   mHH[227]=13*mHH[11];
   mHH[228]=mHH[227] - 19*mHH[13];
   mHH[228]=17./2.*mHH[12] + 1./2.*mHH[228] - 71./3.*mHH[14];
   mHH[228]=1./2.*mHH[7]*mHH[228];
   mHH[117]=mHH[127] + mHH[228] + mHH[117] + mHH[226] + 5./2.*mHH[224]
    + mHH[222] + mHH[221] + 253./16.*mHH[17] - 19./12.*mHH[43] + 
   mHH[220] + mHH[219] + mHH[218] + mHH[217] + 1./4.*mHH[215] + 
   mHH[216];
   mHH[117]=MMZ*mHH[117];
   mHH[127]= - 1./2.*mHH[19];
   mHH[215]=mHH[193] - 5./2.*mHH[76] - 1./4.*mHH[78] + mHH[127] - 35./4.
    - 6*mHH[75];
   mHH[215]=MMZ*mHH[215];
   mHH[222]= - 6*mHH[86];
   mHH[224]=3*mHH[25];
   mHH[229]=mHH[7]*mHH[11];
   mHH[230]=3*mHH[229];
   mHH[231]= - 2*mHH[14];
   mHH[232]=1./2.*mHH[5];
   mHH[233]=3*mHH[42];
   mHH[234]=3*mHH[11];
   mHH[215]=mHH[215] + mHH[230] + 17*mHH[12] + mHH[231] + 10*mHH[13] + 
   mHH[234] + mHH[169] + mHH[87] + mHH[233] - mHH[21] + mHH[224] + 
   mHH[222] + 5./2.*mHH[28] + mHH[232] - 5*mHH[89];
   mHH[235]=pow(MMZ,2);
   mHH[215]=mHH[1]*mHH[235]*mHH[215];
   mHH[236]=mHH[11] - 55*mHH[13];
   mHH[236]=mHH[13]*mHH[236];
   mHH[237]=7*mHH[13];
   mHH[238]=mHH[237] - 185./3.*mHH[14];
   mHH[238]=mHH[14]*mHH[238];
   mHH[236]=1./6.*mHH[236] + mHH[238];
   mHH[238]= - 71./12.*mHH[14];
   mHH[239]=53./24.*mHH[12] + mHH[238] + 1./24.*mHH[11] + 6*mHH[13];
   mHH[239]=mHH[12]*mHH[239];
   mHH[117]=mHH[215] + mHH[117] + 1./4.*mHH[236] + mHH[239];
   mHH[117]=mHH[1]*mHH[117];
   mHH[215]= - 1./4.*mHH[33];
   mHH[236]=mHH[132] + 19 + mHH[215];
   mHH[236]=1./3.*mHH[236] + mHH[15];
   mHH[240]=17./36.*mHH[10];
   mHH[241]= - 17./36.*mHH[34];
   mHH[242]= - 21./4.*mHH[6];
   mHH[243]= - 1./2.*mHH[38];
   mHH[244]=1./3.*mHH[31];
   mHH[245]= - 1./2.*mHH[30];
   mHH[236]=mHH[243] + mHH[242] + mHH[241] + mHH[245] + mHH[240] + 1./2.
   *mHH[236] + mHH[244];
   mHH[246]=mHH[186] + 26./9. - mHH[15];
   mHH[246]=1./3.*mHH[4]*mHH[246];
   mHH[247]=1./8.*mHH[8];
   mHH[236]=mHH[247] + 1./2.*mHH[236] + mHH[246];
   mHH[236]=mHH[8]*mHH[236];
   mHH[248]=53 - 1./3.*mHH[33];
   mHH[249]=13./3.*mHH[32];
   mHH[250]= - 7./3.*mHH[31];
   mHH[251]=17./18.*mHH[10];
   mHH[248]=mHH[251] + mHH[250] + mHH[15] + 1./4.*mHH[248] + mHH[249];
   mHH[252]=7./3.*mHH[30];
   mHH[253]= - 21./2.*mHH[6];
   mHH[248]=mHH[243] + mHH[253] + mHH[241] + 1./2.*mHH[248] + mHH[252];
   mHH[254]= - 5./6.*mHH[10];
   mHH[255]= - 1./2.*mHH[15];
   mHH[256]=mHH[254] + 4./9. + mHH[255];
   mHH[256]=1./3.*mHH[4]*mHH[256];
   mHH[257]= - 7./8.*mHH[8];
   mHH[258]=5./8.*mHH[183];
   mHH[259]= - 1./2.*mHH[7];
   mHH[248]=mHH[259] + mHH[258] + mHH[257] + 1./4.*mHH[248] + mHH[256];
   mHH[248]=mHH[7]*mHH[248];
   mHH[260]=21*mHH[63];
   mHH[261]= - 15*mHH[70];
   mHH[262]=21*mHH[64];
   mHH[263]=mHH[261] + mHH[260] + 3*mHH[67] - 11./4. + mHH[262];
   mHH[264]= - 3./2.*mHH[15];
   mHH[265]= - 5./2.*mHH[10];
   mHH[266]=mHH[265] + 4./3. + mHH[264];
   mHH[266]=mHH[6]*mHH[266];
   mHH[266]=5./6. + mHH[266];
   mHH[266]=mHH[4]*mHH[266];
   mHH[267]= - 3./8.*mHH[31] + 3./8.*mHH[15] + 13./8.*mHH[32] + 2 - 1./
   32.*mHH[33];
   mHH[267]= - 17./16.*mHH[34] - 9./8.*mHH[30] + 3*mHH[267] + 17./16.*
   mHH[10];
   mHH[267]=mHH[6]*mHH[267];
   mHH[268]=mHH[52] - 1./2.*mHH[56] + mHH[195] + mHH[57] + mHH[55];
   mHH[269]=MMH*mHH[268];
   mHH[270]=2*mHH[66];
   mHH[271]= - 21./4.*mHH[39];
   mHH[272]=3./8.*mHH[76];
   mHH[273]= - 5./32.*mHH[68];
   mHH[274]= - 1./8.*mHH[40];
   mHH[275]= - 5./32.*mHH[82];
   mHH[276]= - 5./8.*mHH[95]*mHH[44];
   mHH[277]=5./8.*mHH[108];
   mHH[278]=2*mHH[69];
   mHH[279]= - mHH[38]*mHH[6];
   mHH[280]=9./8.*mHH[279];
   mHH[281]= - 1./12.*mHH[31];
   mHH[282]= - 1./12.*mHH[30];
   mHH[236]=1./2.*mHH[269] + mHH[248] + mHH[277] + mHH[236] + mHH[266]
    + mHH[280] + mHH[267] + mHH[276] + mHH[282] + mHH[281] + mHH[275]
    - 71./32.*mHH[16] + mHH[274] + mHH[273] + mHH[272] + mHH[271] - 15./
   8.*mHH[78] + mHH[270] + 1./8.*mHH[263] + mHH[278];
   mHH[236]=MMH*mHH[236];
   mHH[248]= - 1 + mHH[33];
   mHH[263]= - 85./3.*mHH[10];
   mHH[267]= - 17*mHH[32];
   mHH[248]=mHH[263] + 5./2.*mHH[248] + mHH[267];
   mHH[269]=15*mHH[30];
   mHH[283]=85./6.*mHH[34];
   mHH[248]=mHH[283] + 1./2.*mHH[248] + mHH[269];
   mHH[284]=5*mHH[10];
   mHH[285]=4 + mHH[284];
   mHH[285]=mHH[4]*mHH[285];
   mHH[286]=1./4.*mHH[7];
   mHH[248]=mHH[286] + 25./4.*mHH[108] + 1./4.*mHH[248] + 5./3.*
   mHH[285];
   mHH[248]=mHH[7]*mHH[248];
   mHH[285]=mHH[41] - 3 - mHH[73];
   mHH[285]=5./4.*mHH[82] - mHH[16] - 1./2.*mHH[77] + 1./4.*mHH[285] + 
   mHH[68];
   mHH[285]=mHH[95]*mHH[285];
   mHH[287]=25./32.*mHH[143];
   mHH[288]= - mHH[7]*mHH[95];
   mHH[194]=25./8.*mHH[288] + mHH[287] + mHH[194] + 25./8.*mHH[285];
   mHH[194]=MMZ*mHH[194];
   mHH[289]=13./6.*mHH[19];
   mHH[290]=5*mHH[16];
   mHH[291]=mHH[10] + mHH[290] + mHH[289] + 11./3. + mHH[198];
   mHH[291]=mHH[4]*mHH[291];
   mHH[292]= - 161./12. - 5*mHH[67];
   mHH[292]=1./2.*mHH[292] + 85./3.*mHH[18];
   mHH[293]=4./3.*mHH[75];
   mHH[294]= - 3./2.*mHH[66];
   mHH[295]=1./16.*mHH[33];
   mHH[296]= - 17./24.*mHH[10];
   mHH[297]=3./4.*mHH[30];
   mHH[298]=17./24.*mHH[34];
   mHH[299]= - 1./2.*mHH[17];
   mHH[300]=5*mHH[44] + mHH[299];
   mHH[300]=mHH[95]*mHH[300];
   mHH[301]=1./4. + mHH[108];
   mHH[301]=mHH[9]*mHH[301];
   mHH[302]=mHH[95] + mHH[154];
   mHH[302]=mHH[12]*mHH[302];
   mHH[194]=1./2.*mHH[194] + mHH[248] + 25./32.*mHH[302] + 25./16.*
   mHH[301] + 125./16.*mHH[183] + 5./3.*mHH[291] + 25./16.*mHH[300] + 
   mHH[298] + mHH[297] + mHH[296] + mHH[295] + 125./64.*mHH[82] + 1./4.
   *mHH[16] + 829./288.*mHH[77] - 95./48.*mHH[68] + 25./64.*mHH[41] - 
   25./64.*mHH[73] + 19./96.*mHH[76] - 221./144.*mHH[19] + mHH[294] + 1.
   /8.*mHH[292] + mHH[293];
   mHH[194]=MMZ*mHH[194];
   mHH[248]=187./2.*mHH[5] + 7*mHH[89];
   mHH[197]=mHH[209] + mHH[197] + 1./3.*mHH[248] - 13./4.*mHH[28];
   mHH[209]=8./3.*mHH[21];
   mHH[248]=5*mHH[17];
   mHH[291]=mHH[248] + mHH[214] + mHH[209];
   mHH[291]=mHH[4]*mHH[291];
   mHH[292]= - 293./6. + 5*mHH[33];
   mHH[263]=mHH[263] + 1./2.*mHH[292] + mHH[267];
   mHH[263]=mHH[283] + 1./2.*mHH[263] + mHH[269];
   mHH[269]= - 4./3. + mHH[284];
   mHH[269]=mHH[4]*mHH[269];
   mHH[263]= - 25./32.*mHH[9] + 1./4.*mHH[263] + 5./3.*mHH[269];
   mHH[263]=mHH[12]*mHH[263];
   mHH[227]=mHH[227] - 35*mHH[13];
   mHH[269]=1./2.*mHH[12];
   mHH[283]=19*mHH[14];
   mHH[227]=mHH[269] + 1./2.*mHH[227] + mHH[283];
   mHH[227]=mHH[7]*mHH[227];
   mHH[179]=5./3. + mHH[179];
   mHH[179]=mHH[13]*mHH[179];
   mHH[292]=17*mHH[36];
   mHH[303]= - 395./12. + mHH[292];
   mHH[303]=mHH[14]*mHH[303];
   mHH[304]=mHH[9]*mHH[14];
   mHH[179]=mHH[194] + 1./4.*mHH[227] + mHH[263] + 25./16.*mHH[304] + 5.
   /12.*mHH[303] + 1./4.*mHH[179] + 5./3.*mHH[291] + 41./48.*mHH[11] + 
   391./96.*mHH[17] - 5./24.*mHH[43] + 1499./144.*mHH[44] - 289./48.*
   mHH[26] + 221./36.*mHH[87] + 77./48.*mHH[42] + 1./8.*mHH[197] - 17./
   9.*mHH[21];
   mHH[179]=MMZ*mHH[179];
   mHH[194]=1./2.*mHH[11];
   mHH[197]= - 43./2.*mHH[12] + 215./3.*mHH[14] + mHH[194] - 91*mHH[13]
   ;
   mHH[197]=mHH[12]*mHH[197];
   mHH[227]=pow(mHH[13],2);
   mHH[263]=pow(mHH[14],2);
   mHH[197]=1./8.*mHH[197] - mHH[227] - 119./24.*mHH[263];
   mHH[179]=1./3.*mHH[197] + mHH[179];
   mHH[179]=MMZ*mHH[179];
   mHH[174]=mHH[174] - 17./3.*mHH[5] - mHH[89];
   mHH[174]=1./4.*mHH[174] + mHH[149];
   mHH[197]= - 3*mHH[75];
   mHH[291]= - 1 - mHH[19];
   mHH[291]=mHH[4]*mHH[291];
   mHH[291]=10./3.*mHH[291] - 17./12.*mHH[77] - 1./8.*mHH[76] + 17./12.
   *mHH[19] - 25./8. + mHH[197];
   mHH[291]=MMZ*mHH[291];
   mHH[303]=17./6.*mHH[21];
   mHH[305]=3./2.*mHH[42];
   mHH[306]= - 17./6.*mHH[87];
   mHH[307]=17./12.*mHH[26];
   mHH[308]=3./2.*mHH[11];
   mHH[309]= - 2*mHH[21];
   mHH[310]=mHH[5] + mHH[309];
   mHH[310]=mHH[4]*mHH[310];
   mHH[311]=5./4. + 4./3.*mHH[4];
   mHH[311]=mHH[12]*mHH[311];
   mHH[229]=3./2.*mHH[229];
   mHH[312]=1./2.*mHH[13];
   mHH[291]=mHH[291] + mHH[229] + 5*mHH[311] + 17./3.*mHH[14] + 
   mHH[312] + 10./3.*mHH[310] + mHH[308] + mHH[307] + mHH[306] + 
   mHH[305] + mHH[303] + mHH[174] + 3./2.*mHH[25];
   mHH[291]=mHH[1]*mHH[291]*pow(MMZ,3);
   mHH[179]=mHH[179] + mHH[291];
   mHH[179]=mHH[1]*mHH[179];
   mHH[291]=61./3. - mHH[33];
   mHH[291]=1./2.*mHH[291] + 17./3.*mHH[32];
   mHH[291]=17./12.*mHH[10] + 1./4.*mHH[291] + mHH[31];
   mHH[310]=27./8.*mHH[6];
   mHH[311]=3./8.*mHH[103];
   mHH[313]= - 4./3. - mHH[10];
   mHH[313]=mHH[4]*mHH[313];
   mHH[314]=1./4.*mHH[8];
   mHH[315]=3./4.*mHH[110];
   mHH[316]=3./16.*mHH[7];
   mHH[317]= - 2*mHH[30];
   mHH[291]=mHH[316] + mHH[315] + 25./12.*mHH[183] + mHH[314] + 5./3.*
   mHH[313] + mHH[311] + mHH[310] - 17./24.*mHH[34] + 1./2.*mHH[291] + 
   mHH[317];
   mHH[291]=mHH[7]*mHH[291];
   mHH[174]=mHH[148] - 135./4.*mHH[17] + mHH[307] + mHH[306] + mHH[147]
    + mHH[303] + mHH[174] + mHH[146];
   mHH[174]=mHH[97]*mHH[174];
   mHH[116]=415./72. + mHH[116];
   mHH[116]= - 17./3.*mHH[18] - 3*mHH[70] + 1./2.*mHH[116] + mHH[119];
   mHH[116]=1./2.*mHH[116] + mHH[159];
   mHH[119]=5./6.*mHH[17] + 1./6.*mHH[26] + mHH[207];
   mHH[119]=mHH[95]*mHH[119];
   mHH[116]=1./2.*mHH[174] + 1./4.*mHH[176] + mHH[113] + 25./16.*
   mHH[119] + 17./18.*mHH[34] + mHH[180] - 17./18.*mHH[10] + 1./6.*
   mHH[32] + 1./12.*mHH[33] - 175./96.*mHH[82] - 13./12.*mHH[16] + 3./4.
   *mHH[40] - 523./144.*mHH[77] + mHH[173] + 111./32.*mHH[41] - 111./32.
   *mHH[73] - 37./12.*mHH[76] + mHH[178] - 3./4.*mHH[78] + 17./18.*
   mHH[19] + 1./2.*mHH[116] + mHH[66];
   mHH[119]= - 17./12.*mHH[18] - 27./2.*mHH[63] + 9 + 1./8.*mHH[67];
   mHH[119]=mHH[113] + mHH[298] + mHH[297] + mHH[296] + mHH[295] - 25./
   16.*mHH[16] - 17./16.*mHH[77] + 17./24.*mHH[68] - 3./32.*mHH[76] + 
   mHH[178] + 17./16.*mHH[19] + 3./2.*mHH[66] + 1./2.*mHH[119] - 9*
   mHH[75];
   mHH[119]=mHH[97]*mHH[119];
   mHH[174]=17./6.*mHH[10];
   mHH[176]= - 17./6.*mHH[34];
   mHH[295]=mHH[176] + mHH[114] + mHH[174] - 3 + mHH[215];
   mHH[295]=mHH[97]*mHH[295];
   mHH[295]=25./12.*mHH[95] + mHH[295];
   mHH[296]= - mHH[4]*mHH[97]*mHH[10];
   mHH[295]=1./4.*mHH[295] + 5./3.*mHH[296];
   mHH[295]=mHH[7]*mHH[295];
   mHH[298]=9*mHH[49];
   mHH[303]=mHH[298] + mHH[195];
   mHH[306]= - mHH[41] + 5./6. + mHH[73];
   mHH[306]= - 7./12.*mHH[82] + 1./3.*mHH[16] + mHH[193] + 1./4.*
   mHH[306] - 1./3.*mHH[68];
   mHH[306]=mHH[95]*mHH[306];
   mHH[307]= - 1 - mHH[76];
   mHH[313]=mHH[96]*mHH[307];
   mHH[318]=1./3.*mHH[10] - 1./3.*mHH[16] + 1./3.*mHH[18] + mHH[127];
   mHH[318]=mHH[4]*mHH[97]*mHH[318];
   mHH[119]=mHH[295] + 25./64.*mHH[154] + 5*mHH[318] + mHH[119] + 1./16.
   *mHH[313] + 25./16.*mHH[306] + 1./2.*mHH[303] - mHH[52];
   mHH[119]=MMZ*mHH[119];
   mHH[295]=mHH[150] + mHH[176] + mHH[114] + mHH[174] + 39 + mHH[215];
   mHH[295]=mHH[97]*mHH[295];
   mHH[303]= - 175./24.*mHH[95] + mHH[96];
   mHH[295]=1./2.*mHH[303] + mHH[295];
   mHH[303]=2 - mHH[10];
   mHH[303]=mHH[4]*mHH[97]*mHH[303];
   mHH[287]=mHH[287] + 1./4.*mHH[295] + 5./3.*mHH[303];
   mHH[287]=mHH[12]*mHH[287];
   mHH[295]=mHH[82] - mHH[16] - 1 + mHH[68];
   mHH[295]=mHH[95]*mHH[295];
   mHH[306]=25./48.*mHH[288] + 25./48.*mHH[295] - mHH[52] + mHH[145] - 
   9*mHH[49] + mHH[195];
   mHH[306]=MMH*mHH[306];
   mHH[313]= - 2./3.*mHH[19];
   mHH[318]=2./3.*mHH[10];
   mHH[319]= - mHH[17] + mHH[232] - mHH[21];
   mHH[319]=mHH[97]*mHH[319];
   mHH[320]=mHH[319] + mHH[318] - mHH[16] + mHH[313] - 89./36. + 
   mHH[18];
   mHH[320]=mHH[4]*mHH[320];
   mHH[321]=325./8.*mHH[95] + 17*mHH[97];
   mHH[321]=mHH[14]*mHH[321];
   mHH[322]=25./2.*mHH[183] - 25./8. - 17./3.*mHH[36];
   mHH[322]=mHH[9]*mHH[322];
   mHH[323]= - mHH[8]*mHH[32];
   mHH[116]=mHH[119] + 1./4.*mHH[306] + mHH[291] + mHH[287] + 1./8.*
   mHH[322] + 1./12.*mHH[321] + 1./8.*mHH[177] + 1./24.*mHH[323] + 1./2.
   *mHH[116] + 5./3.*mHH[320];
   mHH[116]=MMZ*mHH[116];
   mHH[119]=7./2. + mHH[267];
   mHH[119]=1./4.*mHH[119] - 7*mHH[31];
   mHH[119]= - 21./8.*mHH[6] + 1./12.*mHH[119] + mHH[30];
   mHH[177]=25./48.*mHH[108];
   mHH[287]= - 1./4.*mHH[7];
   mHH[119]=mHH[287] + mHH[177] + mHH[161] + 1./2.*mHH[119] + 5./9.*
   mHH[4];
   mHH[119]=mHH[7]*mHH[119];
   mHH[291]= - 1 - 17./4.*mHH[32];
   mHH[161]=1./3.*mHH[291] + mHH[161];
   mHH[161]=mHH[8]*mHH[161];
   mHH[306]= - 1./3. + mHH[260];
   mHH[306]=1./2.*mHH[306] + 15*mHH[70];
   mHH[291]=mHH[6]*mHH[291];
   mHH[320]=1./2.*mHH[56] + mHH[52];
   mHH[320]=MMH*mHH[320];
   mHH[321]= - 1./24.*mHH[30];
   mHH[322]=mHH[95]*mHH[44];
   mHH[119]=1./4.*mHH[320] + mHH[119] + 25./48.*mHH[183] + 1./4.*
   mHH[161] + 25./36.*mHH[4] + 3./8.*mHH[291] + 25./48.*mHH[322] + 
   mHH[321] + 25./192.*mHH[82] - 217./192.*mHH[16] - 15./8.*mHH[40] + 
   25./192.*mHH[68] + 3./2.*mHH[76] - 21./16.*mHH[39] - 1./8.*mHH[78]
    + 1./8.*mHH[306] + mHH[66];
   mHH[119]=MMH*mHH[119];
   mHH[161]= - 61./3.*mHH[89] + 13*mHH[28];
   mHH[291]=7./3. + mHH[267];
   mHH[291]=mHH[11]*mHH[291];
   mHH[306]=17./8.*mHH[25];
   mHH[161]=1./8.*mHH[291] - 827./96.*mHH[17] - 11./4.*mHH[43] - 35./48.
   *mHH[44] + 25./96.*mHH[26] - 17./18.*mHH[87] - 9./8.*mHH[42] + 17./
   18.*mHH[21] + mHH[306] + 1./4.*mHH[161] - mHH[86];
   mHH[291]=347./4. - 17./3.*mHH[32];
   mHH[291]=1./8.*mHH[103] - 5./2.*mHH[30] + 1./8.*mHH[291] + mHH[31];
   mHH[320]=1./16.*mHH[110];
   mHH[177]=mHH[320] - 197./96.*mHH[9] + mHH[177] - 1./6.*mHH[8] + 1./2.
   *mHH[291] + 20./9.*mHH[4];
   mHH[177]=mHH[12]*mHH[177];
   mHH[291]= - mHH[4]*mHH[21];
   mHH[324]=139 + mHH[267];
   mHH[324]=1./4.*mHH[324] + mHH[8];
   mHH[324]=mHH[13]*mHH[324];
   mHH[325]=499./6. + 25*mHH[183];
   mHH[325]=mHH[14]*mHH[325];
   mHH[326]=19./3.*mHH[13];
   mHH[327]= - 19./3.*mHH[14] - mHH[11] + mHH[326];
   mHH[327]=1./2.*mHH[327] - mHH[12];
   mHH[327]=mHH[7]*mHH[327];
   mHH[116]=mHH[116] + mHH[119] + 1./2.*mHH[327] + mHH[177] + 61./48.*
   mHH[304] + 1./48.*mHH[325] + 1./6.*mHH[324] + 1./2.*mHH[161] + 10./9.
   *mHH[291];
   mHH[116]=MMZ*mHH[116];
   mHH[119]=1 - 3./2.*mHH[35];
   mHH[161]=1./2.*mHH[119];
   mHH[177]=9*mHH[6];
   mHH[291]=3./2.*mHH[8];
   mHH[324]=mHH[291] + mHH[177] + mHH[161] + mHH[31];
   mHH[324]=mHH[9]*mHH[324];
   mHH[325]=19./3. + mHH[38];
   mHH[327]=3*mHH[9];
   mHH[325]=1./2.*mHH[325] + mHH[327];
   mHH[325]=mHH[7]*mHH[325];
   mHH[328]=127./4. - 9*mHH[72];
   mHH[329]= - 3./2.*mHH[71];
   mHH[330]=1./4.*mHH[38];
   mHH[331]=1 + mHH[330];
   mHH[331]=mHH[8]*mHH[331];
   mHH[332]=9./2.*mHH[6];
   mHH[333]=3./2.*mHH[40];
   mHH[334]=mHH[38]*mHH[6];
   mHH[335]=9./4.*mHH[41];
   mHH[324]=1./4.*mHH[325] + mHH[324] + mHH[331] + 9./8.*mHH[334] + 
   mHH[332] + 25./24.*mHH[82] - 7./24.*mHH[16] + mHH[333] + 7./24.*
   mHH[68] + mHH[335] + mHH[329] + 9*mHH[39] + 1./4.*mHH[328] - 9*
   mHH[65];
   mHH[324]=MMH*mHH[324];
   mHH[325]= - 61./2. + 27*mHH[65];
   mHH[328]=mHH[8]*mHH[38];
   mHH[331]=1./2.*mHH[328];
   mHH[336]=3*mHH[81];
   mHH[337]= - 1./2.*mHH[29];
   mHH[338]=mHH[337] + mHH[80];
   mHH[339]=mHH[338] + mHH[43];
   mHH[339]=mHH[96]*mHH[339];
   mHH[340]=3./2.*mHH[339];
   mHH[341]=mHH[38]*mHH[96];
   mHH[342]= - mHH[96] + 1./2.*mHH[341];
   mHH[342]=mHH[13]*mHH[342];
   mHH[343]=3*mHH[342];
   mHH[344]= - 3./2.*mHH[14]*mHH[96];
   mHH[345]= - 27./2.*mHH[6];
   mHH[136]=mHH[344] + mHH[343] + mHH[331] - mHH[38] + mHH[340] + 
   mHH[345] - 41./12.*mHH[82] - 9./2.*mHH[40] + 3./2.*mHH[83] + 27./4.*
   mHH[41] + mHH[136] + 9./2.*mHH[71] + mHH[336] + 3./2.*mHH[35] - 27./
   2.*mHH[39] - 15./4.*mHH[74] + 1./2.*mHH[325] - 9*mHH[84];
   mHH[325]= - 5./2.*mHH[85] + 11*mHH[24];
   mHH[346]=3*mHH[29];
   mHH[347]= - 1./4.*mHH[87];
   mHH[325]=3./4.*mHH[43] - 33./2.*mHH[44] + mHH[347] - 7./4.*mHH[80]
    - 17./4.*mHH[42] + 1./2.*mHH[325] + mHH[346];
   mHH[348]=3*mHH[325];
   mHH[349]= - 3./2. - mHH[38];
   mHH[349]=mHH[13]*mHH[349];
   mHH[350]= - 27./4.*mHH[11];
   mHH[349]=9./2.*mHH[349] + mHH[350] + mHH[348] + 37./12.*mHH[17];
   mHH[351]= - 1 + 7*mHH[36];
   mHH[352]=17*mHH[32];
   mHH[353]= - 17./3.*mHH[9];
   mHH[351]=mHH[353] + 5./6.*mHH[351] + mHH[352];
   mHH[351]=mHH[9]*mHH[351];
   mHH[354]=115./6.*mHH[82] + 7*mHH[77] + 7*mHH[41] + 7./2.*mHH[36] + 
   241./12. - 7*mHH[73];
   mHH[355]=mHH[7]*mHH[38];
   mHH[356]=3./2.*mHH[355];
   mHH[351]=mHH[356] + 1./3.*mHH[354] + mHH[351];
   mHH[354]=17./2.*mHH[60] + 7*mHH[54];
   mHH[354]=MMZ*mHH[354];
   mHH[351]=1./2.*mHH[351] + 1./3.*mHH[354];
   mHH[351]=MMZ*mHH[351];
   mHH[354]= - 3*mHH[11];
   mHH[237]=mHH[354] + mHH[237];
   mHH[357]= - 3*mHH[14];
   mHH[237]=5./2.*mHH[237] + mHH[357];
   mHH[237]=mHH[9]*mHH[237];
   mHH[358]=3*mHH[38];
   mHH[359]= - 19./6. + mHH[358];
   mHH[360]=17./3.*mHH[9];
   mHH[359]=1./4.*mHH[359] + mHH[360];
   mHH[359]=mHH[12]*mHH[359];
   mHH[361]=1./2.*mHH[26];
   mHH[362]= - mHH[87] + mHH[361];
   mHH[363]=2*mHH[14];
   mHH[364]=mHH[12] + mHH[362] + mHH[363];
   mHH[365]= - 41./6.*mHH[82] - 8 - 7./6.*mHH[77];
   mHH[365]=MMZ*mHH[365];
   mHH[364]=7./3.*mHH[364] + mHH[365];
   mHH[364]=mHH[1]*MMZ*mHH[364];
   mHH[365]= - mHH[38] - 5./4.*mHH[37] + 7./2. - 5*mHH[35];
   mHH[365]=mHH[14]*mHH[365];
   mHH[237]=mHH[364] + mHH[351] + mHH[359] + mHH[237] + 1./2.*mHH[349]
    + 3*mHH[365];
   mHH[237]=mHH[1]*mHH[237];
   mHH[349]= - 2*mHH[31];
   mHH[351]=3./4.*mHH[38];
   mHH[359]= - 3./2.*mHH[8];
   mHH[364]=mHH[359] + mHH[351] + mHH[345] + mHH[349] + 3./8.*mHH[37]
    - 1 + 21./4.*mHH[35];
   mHH[364]=mHH[9]*mHH[364];
   mHH[365]=2*mHH[85] - mHH[24];
   mHH[366]= - mHH[9]*mHH[11];
   mHH[337]=2*mHH[366] - 9*mHH[14] - mHH[13] - 2*mHH[11] + mHH[80] - 2*
   mHH[42] + 2*mHH[365] + mHH[337];
   mHH[367]= - 1 - mHH[82];
   mHH[368]=MMZ*mHH[367];
   mHH[368]=3*mHH[337] + 7./3.*mHH[368];
   mHH[368]=mHH[1]*mHH[368];
   mHH[369]= - 1./4.*mHH[37];
   mHH[370]= - 1./2.*mHH[82];
   mHH[371]=mHH[370] - 5./8.*mHH[83] - 13./4.*mHH[41] + mHH[369] - 
   mHH[35] + 5./4.*mHH[74] - 5./2.*mHH[84] + 3./2. + 2*mHH[72];
   mHH[369]= - mHH[35] + mHH[369];
   mHH[369]= - mHH[9] + 5*mHH[369] - mHH[38];
   mHH[369]=mHH[9]*mHH[369];
   mHH[369]=mHH[371] + mHH[369];
   mHH[372]=MMZ*mHH[60];
   mHH[368]=mHH[368] + 3*mHH[369] + 7./6.*mHH[372];
   mHH[368]=mHH[1]*mHH[368];
   mHH[369]= - 6*mHH[51] + mHH[59];
   mHH[372]=1 + mHH[83];
   mHH[373]=mHH[96]*mHH[372];
   mHH[374]=1./8.*mHH[373] + mHH[369] - 1./2.*mHH[58];
   mHH[375]=1./2.*mHH[83];
   mHH[376]=mHH[375] + 9./2. + 4*mHH[84];
   mHH[377]=mHH[1]*mHH[376];
   mHH[377]= - 2*mHH[59] + mHH[377];
   mHH[377]=MMt*mHH[1]*mHH[377];
   mHH[368]=3*mHH[377] + 3*mHH[374] + mHH[368];
   mHH[368]=MMt*mHH[368];
   mHH[367]=mHH[97]*mHH[367];
   mHH[367]=7./24.*mHH[367] - 13./8.*mHH[60] + 3*mHH[54];
   mHH[367]=MMZ*mHH[367];
   mHH[377]= - mHH[54] + 9*mHH[51] + 1./2.*mHH[59];
   mHH[378]=MMH*mHH[377];
   mHH[379]=3./4.*mHH[378];
   mHH[380]= - mHH[7]*mHH[38];
   mHH[136]=mHH[368] + mHH[237] + mHH[367] + mHH[379] + 1./4.*mHH[380]
    + 1./2.*mHH[136] + mHH[364];
   mHH[136]=MMt*mHH[136];
   mHH[237]=5./2. - 7*mHH[36];
   mHH[237]=1./3.*mHH[237] + mHH[267];
   mHH[237]=mHH[9]*mHH[237];
   mHH[267]= - 7./2.*mHH[34];
   mHH[364]= - 19 + mHH[267];
   mHH[367]= - 3*mHH[9];
   mHH[364]=mHH[367] + 1./3.*mHH[364] - mHH[38];
   mHH[364]=mHH[7]*mHH[364];
   mHH[368]= - 25*mHH[68] - 97./2.*mHH[41] - 73 + 97./2.*mHH[73];
   mHH[381]=mHH[362] - mHH[17];
   mHH[381]=mHH[97]*mHH[381];
   mHH[382]=mHH[14]*mHH[97];
   mHH[383]= - 1./2.*mHH[34];
   mHH[384]=1 + mHH[383];
   mHH[384]=mHH[12]*mHH[97]*mHH[384];
   mHH[385]= - 11*mHH[77];
   mHH[237]=mHH[364] + 7./3.*mHH[384] + 1./2.*mHH[237] + 7./3.*mHH[382]
    + 7./6.*mHH[381] + 7./9.*mHH[34] - 151./36.*mHH[82] + 25./6.*
   mHH[16] + 1./6.*mHH[368] + mHH[385];
   mHH[364]= - 17 + 7*mHH[68];
   mHH[368]= - 17./3.*mHH[82];
   mHH[381]=7./3.*mHH[34];
   mHH[364]=mHH[381] + mHH[368] - 7./3.*mHH[16] + 1./3.*mHH[364] - 7./2.
   *mHH[77];
   mHH[364]=mHH[97]*mHH[364];
   mHH[382]= - mHH[7]*mHH[97]*mHH[34];
   mHH[364]=7./12.*mHH[382] - 25./3.*mHH[54] + 1./4.*mHH[364];
   mHH[364]=MMZ*mHH[364];
   mHH[382]=MMH*mHH[60];
   mHH[237]=1./2.*mHH[364] + 1./4.*mHH[237] + 1./3.*mHH[382];
   mHH[237]=MMZ*mHH[237];
   mHH[364]=1 + mHH[38];
   mHH[364]=mHH[11]*mHH[364];
   mHH[382]=29./2. + mHH[38];
   mHH[382]=mHH[13]*mHH[382];
   mHH[233]=mHH[233] - 3*mHH[24] - mHH[29];
   mHH[233]=1./2.*mHH[233] - mHH[80];
   mHH[384]=3*mHH[233];
   mHH[386]= - 3./4.*mHH[43];
   mHH[364]=1./2.*mHH[382] + 3./4.*mHH[364] - 13./8.*mHH[17] + mHH[386]
    + 31./3.*mHH[44] + 7./24.*mHH[26] + mHH[384] - 25./6.*mHH[87];
   mHH[382]= - 169./2. + 35*mHH[36];
   mHH[352]=1./6.*mHH[382] + mHH[352];
   mHH[352]=mHH[14]*mHH[352];
   mHH[382]= - 23./3. + 5./4.*mHH[34];
   mHH[387]=3./2.*mHH[38];
   mHH[382]=197./12.*mHH[9] + 7./3.*mHH[382] + mHH[387];
   mHH[382]=mHH[12]*mHH[382];
   mHH[388]=mHH[7]*mHH[14];
   mHH[389]=3./2.*mHH[388];
   mHH[352]=mHH[389] + mHH[382] + 61./6.*mHH[225] + mHH[352] + 373./36.
   *mHH[17] + 109./36.*mHH[44] + 41./9.*mHH[87] - 23./8.*mHH[26];
   mHH[382]=83./3. - 11*mHH[73];
   mHH[292]=143./4.*mHH[41] + 13./4.*mHH[382] + mHH[292];
   mHH[382]=5./4. + 17./3.*mHH[36];
   mHH[382]=mHH[9]*mHH[382];
   mHH[292]=5*mHH[382] + 7./6.*mHH[34] + 511./12.*mHH[82] + 7./2.*
   mHH[16] + 289./12.*mHH[77] + 1./3.*mHH[292] - 7./2.*mHH[68];
   mHH[382]= - 11 + 35./2.*mHH[34];
   mHH[358]=1./3.*mHH[382] + mHH[358];
   mHH[358]=1./4.*mHH[358] + mHH[360];
   mHH[358]=mHH[7]*mHH[358];
   mHH[360]=MMZ*mHH[54];
   mHH[292]=17./3.*mHH[360] + 1./4.*mHH[292] + mHH[358];
   mHH[292]=MMZ*mHH[292];
   mHH[292]=1./2.*mHH[352] + mHH[292];
   mHH[292]=MMZ*mHH[292];
   mHH[352]=mHH[269] + 1./2.*mHH[362] + mHH[14];
   mHH[358]= - 25 - 41./3.*mHH[77];
   mHH[358]=1./2.*mHH[358] + mHH[368];
   mHH[358]=MMZ*mHH[358];
   mHH[352]=41./3.*mHH[352] + 1./2.*mHH[358];
   mHH[352]=mHH[1]*mHH[235]*mHH[352];
   mHH[358]= - 9*mHH[11];
   mHH[360]=mHH[358] + 113*mHH[13];
   mHH[362]=9*mHH[14];
   mHH[360]=1./4.*mHH[360] + mHH[362];
   mHH[360]=mHH[14]*mHH[360];
   mHH[368]=mHH[12]*mHH[14];
   mHH[360]=mHH[360] + 39./4.*mHH[368];
   mHH[292]=mHH[352] + 1./2.*mHH[360] + mHH[292];
   mHH[292]=mHH[1]*mHH[292];
   mHH[352]= - 27./4.*mHH[6];
   mHH[360]=mHH[8] + mHH[352] + mHH[349] + 37./3. + 15./4.*mHH[35];
   mHH[360]=mHH[14]*mHH[360];
   mHH[382]=65./2. + mHH[38];
   mHH[353]=1./2.*mHH[382] + mHH[353];
   mHH[353]=mHH[12]*mHH[353];
   mHH[382]=9*mHH[11];
   mHH[390]=3*mHH[14] + mHH[382] - 35./2.*mHH[13];
   mHH[390]=mHH[9]*mHH[390];
   mHH[391]= - mHH[7]*mHH[14];
   mHH[136]=mHH[136] + mHH[292] + mHH[237] + 1./2.*mHH[324] + 1./4.*
   mHH[391] + 1./4.*mHH[353] + 1./4.*mHH[390] + 1./2.*mHH[364] + 
   mHH[360];
   mHH[136]=MMt*mHH[136];
   mHH[237]= - 101./8. - 9*mHH[62];
   mHH[237]=1./2.*mHH[237] - 3*mHH[64];
   mHH[292]=3./4.*mHH[16];
   mHH[324]=1./6.*mHH[31];
   mHH[237]=1./12.*mHH[30] + mHH[324] + mHH[292] + mHH[333] + 27./8.*
   mHH[39] - 3./4.*mHH[66] - 1./2.*mHH[75] - 3./2.*mHH[69] - 9./16.*
   mHH[63] + 3./8.*mHH[237] - mHH[79];
   mHH[353]= - 1./2. + mHH[31];
   mHH[360]=9./4.*mHH[6];
   mHH[314]= - 1./16.*mHH[7] + mHH[314] + mHH[360] + 1./3.*mHH[353] - 1.
   /4.*mHH[30];
   mHH[314]=mHH[7]*mHH[314];
   mHH[353]=1./2.*mHH[31];
   mHH[364]= - 1 + mHH[353];
   mHH[390]=1./3.*mHH[364] + mHH[332];
   mHH[390]=mHH[8]*mHH[390];
   mHH[392]=63./8.*mHH[6];
   mHH[393]=mHH[392] + 5 + mHH[31];
   mHH[393]=mHH[6]*mHH[393];
   mHH[394]=1./2.*mHH[49];
   mHH[395]=mHH[394] + 27./4.*mHH[48] + mHH[50];
   mHH[395]=1./4.*mHH[52] + 3*mHH[395] + 1./2.*mHH[55];
   mHH[396]=1./4.*MMH*mHH[395];
   mHH[314]=mHH[396] + 1./2.*mHH[314] + 1./2.*mHH[390] + mHH[237] + 3./
   4.*mHH[393];
   mHH[314]=MMH*mHH[314];
   mHH[390]= - 71./16. + mHH[31];
   mHH[390]=1./2.*mHH[390] + mHH[177];
   mHH[390]=mHH[11]*mHH[390];
   mHH[393]= - 33./2.*mHH[42] - 9./2.*mHH[25] + 11./2.*mHH[86] - 9*
   mHH[27] + 9./2.*mHH[23] + 11*mHH[88];
   mHH[397]=mHH[393] + 25./3.*mHH[44];
   mHH[397]=11*mHH[17] + 1./2.*mHH[397] - 17*mHH[43];
   mHH[398]=5./3.*mHH[31];
   mHH[399]= - 33./4.*mHH[8] - 105./8.*mHH[6] + 13./2. + mHH[398];
   mHH[399]=mHH[13]*mHH[399];
   mHH[400]= - 25./6. - 9*mHH[6];
   mHH[400]=1./2.*mHH[400] - mHH[8];
   mHH[400]=mHH[14]*mHH[400];
   mHH[401]= - 13 - 1./3.*mHH[31];
   mHH[401]=23./4.*mHH[8] - 27./8.*mHH[6] + 1./2.*mHH[401] + mHH[30];
   mHH[401]=mHH[12]*mHH[401];
   mHH[402]=19./3.*mHH[14];
   mHH[403]=3./2.*mHH[12] + mHH[402] + 5./3.*mHH[11] - 13./2.*mHH[13];
   mHH[403]=mHH[7]*mHH[403];
   mHH[404]=mHH[8]*mHH[11];
   mHH[405]=5./12.*mHH[404];
   mHH[314]=mHH[314] + 1./8.*mHH[403] + 1./2.*mHH[401] + 1./2.*mHH[400]
    + 1./2.*mHH[399] + mHH[405] + 1./4.*mHH[397] + mHH[390];
   mHH[314]=MMH*mHH[314];
   mHH[390]= - 25./2.*mHH[11] - 71./3.*mHH[13];
   mHH[390]=mHH[13]*mHH[390];
   mHH[397]=pow(mHH[11],2);
   mHH[390]=165./8.*mHH[397] + mHH[390];
   mHH[399]= - 3./2.*mHH[11];
   mHH[400]= - 25./24.*mHH[14] + mHH[399] - mHH[13];
   mHH[400]=mHH[14]*mHH[400];
   mHH[401]= - 5./4.*mHH[12] + 13./12.*mHH[14] - mHH[11] + 43./12.*
   mHH[13];
   mHH[401]=mHH[12]*mHH[401];
   mHH[314]=mHH[314] + 1./2.*mHH[401] + 1./4.*mHH[390] + mHH[400];
   mHH[390]= - mHH[6]*mHH[31];
   mHH[400]= - mHH[8]*mHH[31];
   mHH[401]= - mHH[7]*mHH[31];
   mHH[403]=1./12.*mHH[401] + 1./6.*mHH[400] + 3./4.*mHH[390] - mHH[40]
    - 3./4. + mHH[69];
   mHH[403]=MMH*mHH[403];
   mHH[406]= - 1./2.*mHH[31];
   mHH[407]= - 3./4.*mHH[6];
   mHH[408]=mHH[407] + 1 + mHH[406];
   mHH[408]=mHH[11]*mHH[408];
   mHH[409]=3 - 1./6.*mHH[31];
   mHH[410]=3./8.*mHH[6];
   mHH[411]= - 5./12.*mHH[8] + mHH[409] + mHH[410];
   mHH[411]=mHH[13]*mHH[411];
   mHH[412]= - mHH[12]*mHH[31];
   mHH[413]= - mHH[11] + mHH[13];
   mHH[414]=mHH[7]*mHH[413];
   mHH[415]= - mHH[88] + 1./2.*mHH[27];
   mHH[403]=1./2.*mHH[403] + 1./24.*mHH[414] + 1./12.*mHH[412] + 
   mHH[411] + mHH[405] + 1./2.*mHH[408] - mHH[43] + mHH[415] + 1./2.*
   mHH[42];
   mHH[403]=MMH*mHH[403];
   mHH[405]=7./4.*mHH[11] - mHH[13];
   mHH[405]=mHH[13]*mHH[405];
   mHH[408]=mHH[12]*mHH[413];
   mHH[403]=mHH[403] + 1./12.*mHH[408] - 1./4.*mHH[397] + 1./3.*
   mHH[405];
   mHH[403]=MMH*mHH[403];
   mHH[405]= - mHH[6]*mHH[32];
   mHH[408]= - 1./12.*mHH[32];
   mHH[411]=mHH[408] - mHH[31];
   mHH[411]=mHH[7]*mHH[411];
   mHH[323]=mHH[411] + 1./6.*mHH[323] + 3./4.*mHH[405] - mHH[40] + 
   mHH[76] + 1./4. + mHH[70];
   mHH[323]=MMH*mHH[323];
   mHH[405]= - 1./2.*mHH[11];
   mHH[411]=1./3.*mHH[13];
   mHH[412]=1./6.*mHH[12] + mHH[405] + mHH[411];
   mHH[412]=mHH[7]*mHH[412];
   mHH[414]=mHH[7]*mHH[32];
   mHH[307]=mHH[307] + 1./3.*mHH[414];
   mHH[307]=MMZ*mHH[307];
   mHH[414]= - mHH[11]*mHH[32];
   mHH[416]=5 - 1./3.*mHH[32];
   mHH[416]=mHH[13]*mHH[416];
   mHH[408]=1 + mHH[408];
   mHH[408]=mHH[12]*mHH[408];
   mHH[307]=1./2.*mHH[307] + 1./2.*mHH[323] + mHH[412] + mHH[408] + 1./
   2.*mHH[416] + 1./4.*mHH[414] + mHH[175] - 1./2.*mHH[43];
   mHH[307]=MMZ*mHH[307];
   mHH[323]= - 1 - mHH[31];
   mHH[408]=5./6.*mHH[8] + mHH[323] + mHH[407];
   mHH[408]=mHH[12]*mHH[408];
   mHH[412]=3./4.*mHH[6];
   mHH[414]= - 5./6.*mHH[8] + 1 + mHH[412];
   mHH[414]=mHH[13]*mHH[414];
   mHH[416]=1./3.*mHH[11];
   mHH[417]= - 1./12.*mHH[12] + mHH[416] - 1./4.*mHH[13];
   mHH[417]=mHH[7]*mHH[417];
   mHH[418]=MMH*mHH[7]*mHH[31];
   mHH[408]=1./3.*mHH[418] + mHH[417] + mHH[408] + mHH[414] - mHH[43]
    + mHH[17];
   mHH[408]=MMH*mHH[408];
   mHH[414]= - 1./3.*mHH[12];
   mHH[417]=mHH[414] + mHH[354] + 11./3.*mHH[13];
   mHH[417]=mHH[12]*mHH[417];
   mHH[418]= - 1./3.*mHH[13];
   mHH[419]=1./4.*mHH[11] + mHH[418];
   mHH[419]=mHH[13]*mHH[419];
   mHH[307]=mHH[307] + 1./2.*mHH[408] + mHH[419] + 1./4.*mHH[417];
   mHH[307]=MMZ*mHH[307];
   mHH[408]= - mHH[12]*mHH[32];
   mHH[417]=mHH[13] - mHH[12];
   mHH[419]=mHH[7]*mHH[417];
   mHH[420]= - MMZ*mHH[7]*mHH[32];
   mHH[408]=mHH[420] + mHH[408] + mHH[419];
   mHH[408]=MMZ*mHH[408];
   mHH[417]=mHH[12]*mHH[417];
   mHH[408]=mHH[417] + mHH[408];
   mHH[235]=mHH[1]*mHH[235]*mHH[408];
   mHH[235]=1./2.*mHH[235] + mHH[403] + mHH[307];
   mHH[307]=mHH[40] - mHH[83] - 1./4. - mHH[71];
   mHH[403]=mHH[307] + 3./4.*mHH[334];
   mHH[331]=3*mHH[403] + mHH[331];
   mHH[403]=mHH[9]*mHH[31];
   mHH[408]=1./8.*mHH[355];
   mHH[331]=mHH[408] + 1./2.*mHH[331] + mHH[403];
   mHH[331]=MMH*mHH[331];
   mHH[417]=mHH[11]*mHH[38];
   mHH[338]=mHH[338] + 1./2.*mHH[43];
   mHH[419]=mHH[338] + 1./4.*mHH[417];
   mHH[420]= - 9 + mHH[38];
   mHH[420]=mHH[13]*mHH[420];
   mHH[421]=MMZ*mHH[380];
   mHH[422]=3./2.*mHH[14] + mHH[11] - mHH[13];
   mHH[422]=mHH[9]*mHH[422];
   mHH[423]=mHH[12]*mHH[38];
   mHH[331]=1./2.*mHH[421] + mHH[331] + 1./4.*mHH[423] + mHH[422] + 
   mHH[357] + 3*mHH[419] + 1./2.*mHH[420];
   mHH[419]=mHH[9]*mHH[32];
   mHH[356]=mHH[419] + mHH[356];
   mHH[356]=MMZ*mHH[356];
   mHH[419]=mHH[387] + mHH[9];
   mHH[419]=mHH[12]*mHH[419];
   mHH[420]= - mHH[9]*mHH[13];
   mHH[356]=mHH[356] + mHH[420] + mHH[419];
   mHH[356]=MMZ*mHH[356];
   mHH[356]= - 3*mHH[263] + 1./2.*mHH[356];
   mHH[356]=mHH[1]*mHH[356];
   mHH[419]=mHH[9]*mHH[38];
   mHH[421]=mHH[372] + mHH[419];
   mHH[424]= - mHH[14]*mHH[38];
   mHH[225]=mHH[424] + mHH[225];
   mHH[225]=mHH[1]*mHH[225];
   mHH[424]= - mHH[9]*mHH[38];
   mHH[425]=MMt*mHH[1]*mHH[424];
   mHH[421]=mHH[425] + 1./4.*mHH[421] + mHH[225];
   mHH[421]=MMt*mHH[421];
   mHH[331]=3*mHH[421] + 1./2.*mHH[331] + mHH[356];
   mHH[331]=MMt*mHH[331];
   mHH[356]=9./8.*mHH[6];
   mHH[421]= - 5./4.*mHH[8] + mHH[356] + 3./2. + mHH[31];
   mHH[421]=mHH[14]*mHH[421];
   mHH[413]=1./4.*mHH[9]*mHH[413];
   mHH[426]=1./8.*mHH[388];
   mHH[427]= - 1./4.*MMH*mHH[9]*mHH[31];
   mHH[428]= - 3./2.*mHH[44];
   mHH[421]=mHH[427] + mHH[426] + mHH[413] + mHH[428] + mHH[421];
   mHH[421]=MMH*mHH[421];
   mHH[429]=mHH[9]*mHH[13];
   mHH[430]= - mHH[12]*mHH[9];
   mHH[429]=mHH[429] + mHH[430];
   mHH[431]= - MMZ*mHH[9]*mHH[32];
   mHH[429]=1./2.*mHH[431] + 1./2.*mHH[429] + mHH[391];
   mHH[429]=MMZ*mHH[429];
   mHH[431]=mHH[14]*mHH[32];
   mHH[389]=mHH[431] + mHH[389];
   mHH[389]=MMZ*mHH[389];
   mHH[431]= - mHH[14]*mHH[13];
   mHH[389]=mHH[389] + mHH[431] + 5./2.*mHH[368];
   mHH[389]=mHH[1]*MMZ*mHH[389];
   mHH[431]=7./8.*mHH[11] - mHH[13];
   mHH[431]=mHH[14]*mHH[431];
   mHH[331]=mHH[331] + 1./2.*mHH[389] + 1./4.*mHH[429] + 1./2.*mHH[421]
    + mHH[431] + 1./8.*mHH[368];
   mHH[331]=MMt*mHH[331];
   mHH[235]=1./4.*mHH[235] + mHH[331];
   mHH[235]=mHH[3]*mHH[235];
   mHH[116]=mHH[235] + mHH[136] + mHH[179] + 1./2.*mHH[314] + mHH[116];
   mHH[116]=mHH[3]*mHH[116];
   mHH[136]= - 55*mHH[36];
   mHH[179]= - 15./2.*mHH[37];
   mHH[235]=mHH[33] + mHH[136] - 413./3. + mHH[179];
   mHH[314]= - 26*mHH[32];
   mHH[331]= - 6*mHH[15];
   mHH[389]=6*mHH[31];
   mHH[421]= - 17./3.*mHH[10];
   mHH[429]=17./3.*mHH[34];
   mHH[105]=8*mHH[105];
   mHH[431]= - 25./3.*mHH[9];
   mHH[432]=6*mHH[38];
   mHH[433]= - 3*mHH[8];
   mHH[235]=mHH[431] + mHH[433] + mHH[105] + mHH[432] + mHH[429] + 
   mHH[190] + mHH[421] + mHH[389] + mHH[331] + 1./2.*mHH[235] + 
   mHH[314];
   mHH[235]=mHH[9]*mHH[235];
   mHH[434]= - 745./6. + 27*mHH[74];
   mHH[435]=3./4.*mHH[71];
   mHH[436]= - 3./4.*mHH[37];
   mHH[437]= - 13./8.*mHH[83];
   mHH[438]= - 3./4.*mHH[40];
   mHH[439]= - mHH[8]*mHH[38];
   mHH[440]=9./4.*mHH[439];
   mHH[441]= - 3./4.*mHH[38];
   mHH[385]=mHH[440] + mHH[441] - 155./6.*mHH[82] + mHH[438] + mHH[385]
    + mHH[437] - 149./12.*mHH[41] - 11./2.*mHH[36] + mHH[436] + 11*
   mHH[73] + mHH[435] + 1./4.*mHH[434] + mHH[336];
   mHH[434]=13*mHH[82] + 24 + 11*mHH[77];
   mHH[442]=MMZ*mHH[434];
   mHH[443]=2*mHH[80];
   mHH[444]= - mHH[29] + mHH[443];
   mHH[445]= - 22*mHH[12] - 50*mHH[14] - 6*mHH[13] - 11*mHH[26] + 3*
   mHH[444] + 22*mHH[87];
   mHH[442]=mHH[445] + mHH[442];
   mHH[442]=mHH[1]*mHH[442];
   mHH[446]= - mHH[60] - 22*mHH[54];
   mHH[447]=MMZ*mHH[446];
   mHH[235]=mHH[442] + mHH[447] + mHH[385] + mHH[235];
   mHH[235]=mHH[1]*mHH[235];
   mHH[442]=1 + mHH[336];
   mHH[442]=mHH[40] + mHH[375] + 1./2.*mHH[442] - mHH[71];
   mHH[442]=mHH[96]*mHH[442];
   mHH[447]= - 2./3.*mHH[46] - mHH[45];
   mHH[448]= - mHH[38]*mHH[96];
   mHH[449]=1 + mHH[82];
   mHH[449]=mHH[97]*mHH[449];
   mHH[450]=mHH[8]*mHH[341];
   mHH[447]=3./4.*mHH[450] + 11./4.*mHH[449] + 3./4.*mHH[448] + 3./4.*
   mHH[442] + 6*mHH[54] + 27./4.*mHH[60] + 8*mHH[447] + 9./2.*mHH[58];
   mHH[451]=2*mHH[46] + mHH[45];
   mHH[452]= - 3*mHH[58];
   mHH[451]= - 11*mHH[60] + 32./3.*mHH[451] + mHH[452];
   mHH[453]=53 + mHH[336];
   mHH[453]=22*mHH[82] + 1./2.*mHH[453] + 3*mHH[83];
   mHH[454]=mHH[1]*mHH[453];
   mHH[454]=mHH[451] + mHH[454];
   mHH[454]=MMt*mHH[1]*mHH[454];
   mHH[235]=mHH[454] + mHH[447] + mHH[235];
   mHH[235]=MMt*mHH[235];
   mHH[454]= - 103./2. - 15*mHH[37];
   mHH[454]=mHH[33] + 1./2.*mHH[454] + mHH[136];
   mHH[455]= - 9./2.*mHH[8];
   mHH[454]=mHH[455] + mHH[105] + mHH[432] + mHH[429] + mHH[190] + 
   mHH[421] + mHH[389] + mHH[331] + 1./2.*mHH[454] + mHH[314];
   mHH[454]=mHH[14]*mHH[454];
   mHH[456]= - 307./3. + 19*mHH[73];
   mHH[456]= - 11./2.*mHH[34] - 83./4.*mHH[82] - 33./2.*mHH[16] - 77./4.
   *mHH[77] + 33./2.*mHH[68] - 19./4.*mHH[41] + 1./4.*mHH[456] - 
   mHH[36];
   mHH[457]= - 3./4. - mHH[36];
   mHH[457]=mHH[9]*mHH[457];
   mHH[458]=mHH[456] + 5*mHH[457];
   mHH[459]= - 59./2. - 43*mHH[34];
   mHH[459]= - 2*mHH[9] + 1./3.*mHH[459] + mHH[351];
   mHH[459]=mHH[7]*mHH[459];
   mHH[460]= - MMZ*mHH[54];
   mHH[458]=2*mHH[460] + 1./2.*mHH[458] + mHH[459];
   mHH[458]=MMZ*mHH[458];
   mHH[346]= - 2./3.*mHH[47] + mHH[346];
   mHH[460]= - 11./4.*mHH[80];
   mHH[461]=7./2.*mHH[43];
   mHH[346]= - 149./12.*mHH[17] + mHH[461] - 679./24.*mHH[44] + 177./8.
   *mHH[26] - 58./3.*mHH[87] + 2*mHH[346] + mHH[460];
   mHH[462]=mHH[87] + mHH[169];
   mHH[231]= - mHH[12] + mHH[462] + mHH[231];
   mHH[463]=15 + 13*mHH[77];
   mHH[463]=1./2.*mHH[463] + mHH[82];
   mHH[463]=MMZ*mHH[463];
   mHH[463]=13*mHH[231] + mHH[463];
   mHH[463]=mHH[1]*MMZ*mHH[463];
   mHH[464]= - 23 - 9./2.*mHH[38];
   mHH[464]=mHH[13]*mHH[464];
   mHH[465]=12*mHH[13] - 55./2.*mHH[14];
   mHH[465]=mHH[9]*mHH[465];
   mHH[466]= - 1./4.*mHH[9] + mHH[351] + 8 - 43./3.*mHH[34];
   mHH[466]=mHH[12]*mHH[466];
   mHH[454]=mHH[463] + mHH[458] + mHH[466] + mHH[465] + mHH[454] + 
   mHH[346] + 1./2.*mHH[464];
   mHH[454]=mHH[1]*mHH[454];
   mHH[458]=3*mHH[37];
   mHH[463]=397./3. + mHH[458];
   mHH[464]=11*mHH[36];
   mHH[463]= - 1./2.*mHH[33] + 1./2.*mHH[463] + mHH[464];
   mHH[467]= - 3*mHH[38];
   mHH[463]=mHH[467] + mHH[176] + mHH[114] + mHH[174] + mHH[153] + 
   mHH[151] + 1./2.*mHH[463] + mHH[132];
   mHH[468]=2*mHH[9];
   mHH[463]=mHH[468] + mHH[359] + 1./2.*mHH[463] + mHH[189];
   mHH[463]=mHH[9]*mHH[463];
   mHH[469]= - 21*mHH[81] + 425./4. + 9*mHH[74];
   mHH[469]=19./2.*mHH[73] + 1./2.*mHH[469] + 3*mHH[71];
   mHH[462]=mHH[462] + mHH[17];
   mHH[462]=mHH[97]*mHH[462];
   mHH[470]= - 3./2.*mHH[40];
   mHH[469]=11./2.*mHH[462] - mHH[38] + mHH[340] - 11./3.*mHH[34] + 83./
   12.*mHH[82] - 5./2.*mHH[16] + mHH[470] + 9*mHH[77] + 5./2.*mHH[68]
    + mHH[83] + 1./2.*mHH[469] - 7*mHH[41];
   mHH[469]=1./2.*mHH[469];
   mHH[471]=71 + 53*mHH[34];
   mHH[471]=mHH[367] + 1./9.*mHH[471] + mHH[243];
   mHH[471]=1./2.*mHH[7]*mHH[471];
   mHH[213]= - 11*mHH[34] + mHH[82] + mHH[213] + 33./2.*mHH[77] + 1 - 
   11*mHH[68];
   mHH[213]=mHH[97]*mHH[213];
   mHH[472]=mHH[7]*mHH[97]*mHH[34];
   mHH[213]=11./4.*mHH[472] + 5*mHH[54] + 1./4.*mHH[213];
   mHH[473]=MMZ*mHH[213];
   mHH[474]= - 2 + mHH[330];
   mHH[474]=mHH[8]*mHH[474];
   mHH[475]=4./3.*mHH[45] - mHH[60];
   mHH[476]=MMH*mHH[475];
   mHH[477]= - 3./2.*mHH[96] - 11*mHH[97];
   mHH[477]=1./2.*mHH[14]*mHH[477];
   mHH[478]= - 1 + mHH[134];
   mHH[478]=mHH[12]*mHH[97]*mHH[478];
   mHH[479]=11./2.*mHH[478];
   mHH[342]=3./2.*mHH[342];
   mHH[235]=mHH[235] + mHH[454] + mHH[473] + mHH[476] + mHH[471] + 
   mHH[479] + mHH[463] + mHH[477] + mHH[342] + mHH[469] + mHH[474];
   mHH[235]=MMt*mHH[235];
   mHH[454]= - 467 - mHH[33];
   mHH[454]=1./4.*mHH[454] + mHH[132];
   mHH[454]=1./3.*mHH[454] + mHH[15];
   mHH[186]=mHH[186] + 44./9. - mHH[15];
   mHH[186]=2./3.*mHH[4]*mHH[186];
   mHH[463]=mHH[11]*mHH[96];
   mHH[473]=1./8.*mHH[463];
   mHH[474]= - 1./8.*mHH[13]*mHH[96];
   mHH[454]=mHH[474] + mHH[107] + mHH[186] + mHH[473] + mHH[243] + 
   mHH[241] + mHH[245] + mHH[240] + 1./2.*mHH[454] + mHH[349];
   mHH[454]=mHH[13]*mHH[454];
   mHH[215]=mHH[467] + mHH[176] + mHH[114] + mHH[174] + mHH[153] + 
   mHH[151] + mHH[132] + 67./3. + mHH[215];
   mHH[215]=mHH[11]*mHH[215];
   mHH[476]= - 43./2. - mHH[33];
   mHH[476]=1./4.*mHH[476] + mHH[132];
   mHH[480]= - 9*mHH[30];
   mHH[481]= - 17./18.*mHH[34];
   mHH[482]=1./2.*mHH[103];
   mHH[476]=mHH[482] - mHH[38] + mHH[481] + mHH[480] + mHH[251] + 
   mHH[31] + 1./3.*mHH[476] + mHH[15];
   mHH[188]=1./3.*mHH[188];
   mHH[483]=1./16.*mHH[9];
   mHH[110]=1./8.*mHH[110];
   mHH[476]=mHH[110] + mHH[483] + mHH[258] - 7./2.*mHH[8] + 1./4.*
   mHH[476] + mHH[188];
   mHH[476]=mHH[12]*mHH[476];
   mHH[326]=mHH[399] + mHH[326];
   mHH[326]= - 53./12.*mHH[12] + 1./2.*mHH[326] + 71./9.*mHH[14];
   mHH[326]=1./2.*mHH[7]*mHH[326];
   mHH[399]=mHH[277] - 143./48. + mHH[8];
   mHH[399]=mHH[14]*mHH[399];
   mHH[484]=17./8.*mHH[27];
   mHH[485]=mHH[484] - mHH[88] + 4./3.*mHH[90] + mHH[22];
   mHH[486]= - 1./6.*mHH[21];
   mHH[487]= - 9./4.*mHH[42];
   mHH[488]=1./6.*mHH[87];
   mHH[489]= - 5./32.*mHH[26];
   mHH[490]=79./16.*mHH[44];
   mHH[187]=mHH[11]*mHH[187];
   mHH[187]= - 4./3.*mHH[22] + mHH[187];
   mHH[187]=mHH[4]*mHH[187];
   mHH[491]= - mHH[8]*mHH[11];
   mHH[492]=1./2.*mHH[491];
   mHH[493]= - mHH[13] + 13./8.*mHH[14];
   mHH[493]=3*mHH[9]*mHH[493];
   mHH[99]=mHH[116] + mHH[235] + mHH[117] + mHH[99] + mHH[236] + 
   mHH[326] + mHH[476] + mHH[493] + mHH[399] + mHH[454] + mHH[492] + 
   mHH[187] + 1./4.*mHH[215] - 967./96.*mHH[17] + 19./24.*mHH[43] + 
   mHH[490] + mHH[489] + mHH[488] - mHH[80] + mHH[487] + mHH[486] + 
   mHH[306] - mHH[86] - 3./2.*mHH[28] + mHH[485] + 6*mHH[89];
   mHH[99]=mHH[3]*mHH[99];
   mHH[116]=53*mHH[33];
   mHH[117]=mHH[153] + mHH[151] + mHH[132] + 251./4. + mHH[116];
   mHH[117]=mHH[185] + mHH[184] + 25./2.*mHH[8] + mHH[189] + mHH[182]
    + mHH[181] + mHH[190] + 1./2.*mHH[117] + mHH[160];
   mHH[117]=mHH[7]*mHH[117];
   mHH[184]=3./2.*mHH[15];
   mHH[185]= - 3./2.*mHH[31];
   mHH[215]= - 28./3.*mHH[10] - 20./9. - mHH[15];
   mHH[215]=mHH[4]*mHH[215];
   mHH[215]=35./3.*mHH[8] + 2*mHH[215] + mHH[182] - 28./3.*mHH[34] + 28.
   /3.*mHH[10] + mHH[185] + mHH[184] + 71./8.*mHH[32] + 737./6. + 119*
   mHH[33];
   mHH[215]=mHH[8]*mHH[215];
   mHH[235]= - 11./3. + 5*mHH[94];
   mHH[236]= - 2*mHH[15];
   mHH[235]=mHH[236] - 10*mHH[40] + 2*mHH[235] - 13./3.*mHH[20];
   mHH[235]=mHH[4]*mHH[235];
   mHH[399]= - 3*mHH[94];
   mHH[454]=2251./48. + mHH[399];
   mHH[476]=33./8.*mHH[33];
   mHH[117]=mHH[117] + mHH[212] + mHH[211] + mHH[210] + mHH[215] + 
   mHH[235] + mHH[182] + mHH[208] + mHH[130] + mHH[128] + mHH[129] + 
   mHH[185] + mHH[184] + 165./8.*mHH[32] + mHH[476] + mHH[206] + 20*
   mHH[16] + 269./8.*mHH[40] + mHH[205] + mHH[204] + mHH[203] + 
   mHH[202] + 15./2.*mHH[71] - 11./2.*mHH[81] + 953./16.*mHH[76] + 2317.
   /16.*mHH[78] + mHH[201] + mHH[200] + mHH[199] - 5./4.*mHH[18] + 3*
   mHH[69] - 341./8.*mHH[70] + 13./4.*mHH[20] - 121./8.*mHH[67] + 5./2.
   *mHH[454] - 8./3.*mHH[79];
   mHH[117]=mHH[98]*mHH[117];
   mHH[184]= - 275*mHH[33];
   mHH[185]= - 103 + mHH[184];
   mHH[199]=55*mHH[32];
   mHH[185]=1./2.*mHH[185] + mHH[199];
   mHH[185]= - mHH[31] + 1./4.*mHH[185] + mHH[15];
   mHH[200]= - 25./4.*mHH[10];
   mHH[201]=9./2.*mHH[30];
   mHH[202]=25./4.*mHH[34];
   mHH[185]=mHH[467] + mHH[202] + mHH[201] + 3*mHH[185] + mHH[200];
   mHH[203]=7*mHH[10];
   mHH[204]=2 - mHH[15];
   mHH[205]=2*mHH[204] + mHH[203];
   mHH[205]=mHH[4]*mHH[205];
   mHH[206]= - 15./2.*mHH[8];
   mHH[185]=mHH[286] + 9./4.*mHH[108] + mHH[206] + 1./2.*mHH[185] + 
   mHH[205];
   mHH[185]=mHH[7]*mHH[185];
   mHH[205]=2*mHH[15];
   mHH[198]=mHH[10] + mHH[205] + mHH[290] + 10*mHH[40] + mHH[289] + 
   mHH[198] + 13./3.*mHH[20] + 11 - 10*mHH[94];
   mHH[198]=mHH[4]*mHH[198];
   mHH[208]= - mHH[31] + mHH[15] + 77./4.*mHH[32] - 25 - 121./2.*
   mHH[33];
   mHH[208]=3*mHH[208] - 11*mHH[10];
   mHH[210]=4 - mHH[15];
   mHH[211]=mHH[210] + 6*mHH[10];
   mHH[211]=mHH[4]*mHH[211];
   mHH[208]= - 13./4.*mHH[8] + 2*mHH[211] + mHH[182] + 11./2.*mHH[34]
    + 1./2.*mHH[208] + 3*mHH[30];
   mHH[208]=mHH[8]*mHH[208];
   mHH[211]= - 5467./32. + 15*mHH[94];
   mHH[185]=mHH[185] + 9./32.*mHH[302] + 9./16.*mHH[301] + 45./16.*
   mHH[183] + mHH[208] + mHH[198] + mHH[387] + 9./16.*mHH[300] + 3./8.*
   mHH[34] + mHH[297] - 3./8.*mHH[10] + 3./2.*mHH[31] + mHH[264] - 99./
   8.*mHH[32] - 99./16.*mHH[33] + 45./64.*mHH[82] - 27./2.*mHH[16] - 
   207./8.*mHH[40] + 53./32.*mHH[77] - 21./16.*mHH[68] + 9./64.*mHH[41]
    - 9./64.*mHH[73] - 15./2.*mHH[71] + 11./2.*mHH[81] - 1635./32.*
   mHH[76] - 1635./16.*mHH[78] - 13./16.*mHH[19] + mHH[294] + mHH[293]
    + 15./8.*mHH[18] - 3*mHH[69] + 231./8.*mHH[70] - 13./4.*mHH[20] + 
   231./16.*mHH[67] + 1./2.*mHH[211] + 8./3.*mHH[79];
   mHH[185]=mHH[92]*mHH[185];
   mHH[198]=mHH[33] - mHH[32];
   mHH[198]=33./4.*mHH[198] - mHH[15];
   mHH[208]=1./2.*mHH[10];
   mHH[211]=mHH[383] - mHH[30] + mHH[208] + mHH[198] + mHH[31];
   mHH[212]=mHH[211] + mHH[38];
   mHH[215]=mHH[15] - mHH[10];
   mHH[235]=mHH[4]*mHH[215];
   mHH[264]=3*mHH[212] + 4*mHH[235];
   mHH[264]=mHH[91]*mHH[264];
   mHH[286]= - 53*mHH[33];
   mHH[289]= - 6*mHH[32];
   mHH[290]= - 16./3.*mHH[10];
   mHH[293]=16./3.*mHH[34];
   mHH[294]= - 1./3. + mHH[10];
   mHH[294]=32./3.*mHH[4]*mHH[294];
   mHH[297]= - 26./3.*mHH[8] + mHH[264] + mHH[294] + mHH[293] + 
   mHH[290] + mHH[289] - 256./3. + mHH[286];
   mHH[297]=mHH[8]*mHH[297];
   mHH[300]= - 2*mHH[61];
   mHH[301]=mHH[300] - mHH[57];
   mHH[301]=2*mHH[301] - mHH[55];
   mHH[192]=5./8.*mHH[157] + 5./32.*mHH[154] + 5./8.*mHH[192] + mHH[52]
    + 22*mHH[56] + 2*mHH[301] + 55./2.*mHH[53];
   mHH[192]=mHH[98]*mHH[192];
   mHH[301]=1./2.*mHH[52];
   mHH[285]=3./16.*mHH[288] + 3./64.*mHH[143] + 3./16.*mHH[285] + 
   mHH[301] - 33./4.*mHH[56] + mHH[55] - 33./2.*mHH[53];
   mHH[285]=mHH[92]*mHH[285];
   mHH[302]=2*mHH[61] + mHH[57];
   mHH[302]= - 8*mHH[53] + 4*mHH[302] + mHH[55];
   mHH[192]=3*mHH[285] + 3*mHH[192] + 3*mHH[302] - 41*mHH[56];
   mHH[192]=MMZ*mHH[192];
   mHH[285]=MMW*mHH[61];
   mHH[302]=mHH[57]*MMW;
   mHH[454]= - mHH[56]*MMW;
   mHH[285]=3*mHH[70] + 3*mHH[454] + 3*mHH[302] - 29./3. + 6*mHH[285];
   mHH[302]=3./2.*mHH[212] + 2*mHH[235];
   mHH[302]=mHH[91]*mHH[302];
   mHH[454]=mHH[7]*mHH[302];
   mHH[117]=mHH[192] + mHH[185] + mHH[117] + mHH[454] + mHH[297] + 
   mHH[289] - 6*mHH[40] - 12*mHH[76] + 4*mHH[285] - 143./3.*mHH[78];
   mHH[117]=MMZ*mHH[117];
   mHH[185]=mHH[153] + mHH[151] + mHH[132] - 831./8. + mHH[116];
   mHH[160]=mHH[196] + 155./4.*mHH[8] + mHH[189] + mHH[182] + mHH[181]
    + mHH[190] + 1./2.*mHH[185] + mHH[160];
   mHH[160]=mHH[12]*mHH[160];
   mHH[181]=1./3. - mHH[10];
   mHH[181]=mHH[4]*mHH[181];
   mHH[116]=439./12.*mHH[8] + 32./3.*mHH[181] - 16./3.*mHH[34] + 16./3.
   *mHH[10] + mHH[132] + 55./6. + mHH[116];
   mHH[116]=mHH[13]*mHH[116];
   mHH[181]=5./2.*mHH[223] - 56./3.*mHH[8];
   mHH[181]=mHH[14]*mHH[181];
   mHH[116]=mHH[228] + mHH[160] + mHH[226] + mHH[181] + mHH[116] + 
   mHH[221] + 3223./48.*mHH[17] + 445./12.*mHH[43] + mHH[220] + 
   mHH[219] + mHH[218] + mHH[217] + mHH[216] - 37./8.*mHH[25] + 17./12.
   *mHH[86] - 525./16.*mHH[28] + 141./4.*mHH[89] - 8./3.*mHH[90] - 11./
   8.*mHH[5];
   mHH[116]=mHH[98]*mHH[116];
   mHH[160]= - 2*mHH[22];
   mHH[181]=mHH[93] + mHH[160];
   mHH[182]= - mHH[27] + mHH[181] + 2*mHH[88];
   mHH[185]= - mHH[93] + 2*mHH[22];
   mHH[192]=mHH[4]*mHH[185];
   mHH[196]=93 - 4*mHH[4];
   mHH[196]=mHH[13]*mHH[196];
   mHH[182]=mHH[230] + 111*mHH[12] - 8*mHH[14] + 2*mHH[196] + 3*
   mHH[491] + 4*mHH[192] + mHH[169] + mHH[87] + 6*mHH[80] - 3*mHH[29]
    - mHH[21] + mHH[224] + mHH[222] + 99./2.*mHH[28] - 99*mHH[89] + 3*
   mHH[182] + mHH[232];
   mHH[182]=mHH[98]*mHH[182];
   mHH[192]=99./4.*mHH[89];
   mHH[196]= - 99./8.*mHH[28];
   mHH[216]=1./4.*mHH[26];
   mHH[217]= - 1./2.*mHH[87];
   mHH[218]=1./2.*mHH[25];
   mHH[185]=mHH[308] + mHH[216] + mHH[217] - 2*mHH[80] + mHH[305] + 
   mHH[29] + 1./2.*mHH[21] + mHH[218] - mHH[86] + mHH[196] + mHH[192]
    - 1./4.*mHH[5] + mHH[27] + mHH[185] - 2*mHH[88];
   mHH[219]=mHH[309] + 2*mHH[181] + mHH[5];
   mHH[219]=mHH[4]*mHH[219];
   mHH[220]= - 273./2. + 8*mHH[4];
   mHH[220]=mHH[13]*mHH[220];
   mHH[221]= - 273./4. + 4*mHH[4];
   mHH[221]=mHH[12]*mHH[221];
   mHH[185]=mHH[229] + mHH[221] + mHH[362] + mHH[220] + 3*mHH[404] + 3*
   mHH[185] + 2*mHH[219];
   mHH[185]=mHH[92]*mHH[185];
   mHH[219]= - 2*mHH[20];
   mHH[220]= - 2*mHH[75] + mHH[219] - 289./4. + 4*mHH[79];
   mHH[221]=1 + mHH[20];
   mHH[221]=mHH[4]*mHH[221];
   mHH[127]=8*mHH[221] + mHH[193] + 6*mHH[81] - 99./2.*mHH[76] - 693./4.
   *mHH[78] + 3*mHH[220] + mHH[127];
   mHH[127]=mHH[98]*mHH[127];
   mHH[220]=1./4.*mHH[19];
   mHH[221]= - 1./4.*mHH[77];
   mHH[222]=mHH[221] - mHH[81] + 99./8.*mHH[76] + 99./4.*mHH[78] + 
   mHH[220] - mHH[75] + mHH[20] + 273./8. - 2*mHH[79];
   mHH[219]= - mHH[19] - 3 + mHH[219];
   mHH[219]=mHH[4]*mHH[219];
   mHH[219]=3*mHH[222] + 2*mHH[219];
   mHH[219]=mHH[92]*mHH[219];
   mHH[222]=50 - 3*mHH[79];
   mHH[223]=3*mHH[20];
   mHH[224]= - 1 - mHH[20];
   mHH[224]=mHH[4]*mHH[224];
   mHH[127]=mHH[219] + mHH[127] + 4*mHH[224] - 3*mHH[81] + 12*mHH[76]
    + 94*mHH[78] + 2*mHH[222] + mHH[223];
   mHH[127]=MMZ*mHH[127];
   mHH[219]=MMW + mHH[89];
   mHH[222]=mHH[78]*MMW;
   mHH[219]= - 2*mHH[12] - 4*mHH[13] + 2*mHH[222] + 2*mHH[219] - 
   mHH[28];
   mHH[127]=mHH[127] + mHH[185] + 12*mHH[219] + mHH[182];
   mHH[127]=mHH[1]*MMZ*mHH[127];
   mHH[182]=579./2. + mHH[184];
   mHH[182]=1./2.*mHH[182] + mHH[199];
   mHH[182]= - mHH[31] + 1./4.*mHH[182] + mHH[15];
   mHH[182]=mHH[467] + mHH[202] + mHH[201] + 3*mHH[182] + mHH[200];
   mHH[184]= - 2./3. - mHH[15];
   mHH[184]=2*mHH[184] + mHH[203];
   mHH[184]=mHH[4]*mHH[184];
   mHH[182]= - 9./32.*mHH[9] - 39./4.*mHH[8] + 1./2.*mHH[182] + 
   mHH[184];
   mHH[182]=mHH[12]*mHH[182];
   mHH[184]= - 11*mHH[33];
   mHH[185]=18 + mHH[184];
   mHH[185]=mHH[353] + mHH[255] + 2*mHH[185] + 11./8.*mHH[32];
   mHH[200]= - 4*mHH[10];
   mHH[201]=4*mHH[10] - 4./3. + mHH[15];
   mHH[201]=mHH[4]*mHH[201];
   mHH[185]= - 123./4.*mHH[8] + 2*mHH[201] + mHH[387] + 4*mHH[34] + 3*
   mHH[185] + mHH[200];
   mHH[185]=mHH[13]*mHH[185];
   mHH[201]=mHH[248] + 10*mHH[43] + mHH[209] + mHH[214] - 11*mHH[93] + 
   16./3.*mHH[22];
   mHH[201]=mHH[4]*mHH[201];
   mHH[202]=5*mHH[37];
   mHH[203]=5./2.*mHH[36];
   mHH[209]=mHH[203] - 41./8. + mHH[202];
   mHH[209]=3./2.*mHH[209] + 5*mHH[8];
   mHH[209]=mHH[14]*mHH[209];
   mHH[214]=mHH[11] - 15*mHH[13];
   mHH[214]=33./2.*mHH[12] + 13./2.*mHH[214] + mHH[283];
   mHH[214]=mHH[7]*mHH[214];
   mHH[182]=1./4.*mHH[214] + mHH[182] + 9./16.*mHH[304] + mHH[209] + 
   mHH[185] + 13./4.*mHH[404] + mHH[201] + 41./16.*mHH[11] - 4145./96.*
   mHH[17] - 2059./24.*mHH[43] + 111./16.*mHH[44] - 51./16.*mHH[26] + 
   13./4.*mHH[87] + 17./2.*mHH[80] + 77./16.*mHH[42] - 21./2.*mHH[29]
    - mHH[21] - 37./16.*mHH[25] + 17./24.*mHH[86] + 1575./32.*mHH[28]
    - 423./8.*mHH[89] + 33./16.*mHH[5] - 37./8.*mHH[27] + 17./12.*
   mHH[88] + 33./4.*mHH[93] - 4*mHH[22];
   mHH[182]=mHH[92]*mHH[182];
   mHH[185]= - 1 - mHH[33];
   mHH[201]= - 12 + 17./4.*mHH[91];
   mHH[201]=mHH[8]*mHH[201];
   mHH[185]=mHH[201] + 12*mHH[185] + mHH[264];
   mHH[185]=mHH[13]*mHH[185];
   mHH[201]= - 12 - 17./4.*mHH[91];
   mHH[201]=mHH[8]*mHH[201];
   mHH[201]=mHH[201] + 12 + mHH[302];
   mHH[201]=mHH[12]*mHH[201];
   mHH[209]= - mHH[12]*mHH[91];
   mHH[214]=mHH[13]*mHH[91];
   mHH[219]=mHH[214] + mHH[209];
   mHH[219]=mHH[7]*mHH[219];
   mHH[222]= - mHH[33]*MMW;
   mHH[222]= - MMW + mHH[222];
   mHH[222]=mHH[8]*mHH[222];
   mHH[222]= - mHH[17] + mHH[222];
   mHH[224]=mHH[14]*mHH[8];
   mHH[116]=mHH[127] + mHH[117] + mHH[182] + mHH[116] + 17./8.*mHH[219]
    + mHH[201] + 32./3.*mHH[224] + 12*mHH[222] + mHH[185];
   mHH[116]=mHH[1]*mHH[116];
   mHH[117]= - 331./2. - 73*mHH[33];
   mHH[100]=1./2.*mHH[117] + mHH[100];
   mHH[100]=1./3.*mHH[100] - mHH[15];
   mHH[100]=mHH[112] + mHH[111] + mHH[109] + mHH[107] + mHH[106] + 
   mHH[104] + mHH[115] + mHH[113] + mHH[102] + mHH[114] + mHH[101] + 1./
   2.*mHH[100] + mHH[31];
   mHH[100]=mHH[7]*mHH[100];
   mHH[101]=33./2.*mHH[28] + mHH[5] - 33*mHH[89];
   mHH[101]=mHH[148] - 167./4.*mHH[17] + mHH[142] + mHH[141] + mHH[147]
    + mHH[139] + mHH[146] + 1./4.*mHH[101] + mHH[149];
   mHH[101]=mHH[97]*mHH[101];
   mHH[102]=5 + mHH[131];
   mHH[102]=mHH[150] + mHH[134] + mHH[114] + 11*mHH[102] + mHH[133];
   mHH[102]=mHH[97]*mHH[102];
   mHH[104]=mHH[152] + 33*mHH[96];
   mHH[102]=mHH[155] + 1./2.*mHH[104] + mHH[102];
   mHH[102]=mHH[12]*mHH[102];
   mHH[104]= - 43./4.*mHH[8] + mHH[294] + mHH[293] + mHH[290] + 7./2.*
   mHH[32] - 91 + mHH[286];
   mHH[104]=mHH[8]*mHH[104];
   mHH[106]= - 1051./8. - 9*mHH[67];
   mHH[107]= - 9*mHH[63];
   mHH[106]=1./2.*mHH[106] + mHH[107];
   mHH[106]=mHH[18] + 3*mHH[106] - mHH[70];
   mHH[106]=1./2.*mHH[106] + mHH[159];
   mHH[109]=25 - 11*mHH[32];
   mHH[109]=mHH[96]*mHH[109];
   mHH[109]=mHH[109] + 11*mHH[97];
   mHH[111]=mHH[8]*mHH[96];
   mHH[109]=1./4.*mHH[109] + 2*mHH[111];
   mHH[109]=mHH[13]*mHH[109];
   mHH[112]=11*mHH[175] - 3*mHH[43];
   mHH[112]=mHH[96]*mHH[112];
   mHH[100]=mHH[100] + 1./2.*mHH[102] + mHH[172] + mHH[171] + 3*
   mHH[109] + 1./3.*mHH[104] + 1./2.*mHH[101] + 3./4.*mHH[112] + 
   mHH[113] + mHH[170] + mHH[168] + mHH[180] + mHH[167] + 11./2.*
   mHH[32] + 11./4.*mHH[33] + mHH[166] + 7./4.*mHH[16] + 9./4.*mHH[40]
    + mHH[165] + mHH[173] + mHH[164] + mHH[163] - 89./8.*mHH[76] + 
   mHH[178] - 251./8.*mHH[78] + mHH[162] + 1./2.*mHH[106] + mHH[66];
   mHH[100]=mHH[98]*mHH[100];
   mHH[101]=509./3. + 209*mHH[33];
   mHH[102]= - 55*mHH[32];
   mHH[101]=1./2.*mHH[101] + mHH[102];
   mHH[101]=19./12.*mHH[10] + mHH[31] + 1./4.*mHH[101] - mHH[15];
   mHH[104]= - 2 + mHH[15];
   mHH[104]=2*mHH[104] - 5*mHH[10];
   mHH[104]=mHH[4]*mHH[104];
   mHH[106]=3./4.*mHH[183];
   mHH[101]=mHH[316] + mHH[315] + mHH[106] + mHH[158] + 1./3.*mHH[104]
    + mHH[311] + mHH[115] + mHH[310] - 19./24.*mHH[34] + 1./2.*mHH[101]
    - mHH[30];
   mHH[101]=mHH[7]*mHH[101];
   mHH[104]=mHH[181] - 2*mHH[43];
   mHH[104]=mHH[96]*mHH[104];
   mHH[109]=2*mHH[94];
   mHH[112]= - 2*mHH[40];
   mHH[104]=mHH[319] + mHH[104] + mHH[318] + 4./3.*mHH[15] - mHH[16] + 
   mHH[112] + mHH[313] + mHH[18] - 4./3.*mHH[20] - 89./12. + mHH[109];
   mHH[104]=mHH[4]*mHH[104];
   mHH[117]=mHH[200] - 4 + mHH[15];
   mHH[117]=mHH[4]*mHH[117];
   mHH[113]=mHH[158] + 2./3.*mHH[117] + 3./4.*mHH[463] + mHH[115] + 
   mHH[113] - 4./3.*mHH[34] - mHH[30] + 4./3.*mHH[10] + mHH[255] - 77./
   8.*mHH[32] + 125./6. + 22*mHH[33];
   mHH[113]=mHH[8]*mHH[113];
   mHH[117]= - mHH[88] - 1./2.*mHH[93] + mHH[22];
   mHH[117]= - 35./4.*mHH[43] - 3*mHH[80] + mHH[147] + 3./2.*mHH[29] + 
   mHH[196] + mHH[192] + 3*mHH[117] + 7./4.*mHH[27];
   mHH[117]=mHH[96]*mHH[117];
   mHH[127]= - 33./2.*mHH[28] - mHH[5] + 33*mHH[89];
   mHH[127]=1./4.*mHH[127] - mHH[86];
   mHH[127]=mHH[148] - 35./4.*mHH[17] + 3./4.*mHH[26] - 3./2.*mHH[87]
    + mHH[147] + 3./2.*mHH[21] + 3*mHH[127] + mHH[146];
   mHH[127]=mHH[97]*mHH[127];
   mHH[131]= - 13 + 9*mHH[32];
   mHH[131]=mHH[150] + mHH[153] + 11./4.*mHH[131] + mHH[151];
   mHH[131]=mHH[96]*mHH[131];
   mHH[131]= - 99./4.*mHH[97] + mHH[131] + 3*mHH[448];
   mHH[139]=mHH[4]*mHH[96]*mHH[204];
   mHH[141]= - mHH[8]*mHH[96];
   mHH[131]=3./2.*mHH[141] + 1./2.*mHH[131] + 2*mHH[139];
   mHH[131]=mHH[13]*mHH[131];
   mHH[139]= - mHH[50] - 1./2.*mHH[49];
   mHH[142]= - 1./2.*mHH[52];
   mHH[139]=3./32.*mHH[288] + 3./32.*mHH[295] + mHH[142] - 11./8.*
   mHH[56] - 11./4.*mHH[53] + 9*mHH[139] - mHH[55];
   mHH[139]=MMH*mHH[139];
   mHH[146]= - 1 + 9./4.*mHH[33];
   mHH[146]=mHH[150] - 3./2.*mHH[34] + mHH[114] + 11*mHH[146] + 3./2.*
   mHH[10];
   mHH[146]=mHH[97]*mHH[146];
   mHH[147]= - 7./8.*mHH[95] - 33*mHH[96];
   mHH[146]=3./2.*mHH[147] + mHH[146];
   mHH[143]=9./32.*mHH[143] + 1./4.*mHH[146] + mHH[303];
   mHH[143]=mHH[12]*mHH[143];
   mHH[146]= - 27./4.*mHH[63] - 9./2.*mHH[74] + 29./8.*mHH[67] - 13./3.
   *mHH[79] - 27./2.*mHH[64] + 24869./192. + mHH[399];
   mHH[147]= - 1./4.*mHH[10];
   mHH[148]=1./4.*mHH[34];
   mHH[149]=5./2.*mHH[17] + mHH[361] - 15*mHH[44];
   mHH[149]=mHH[95]*mHH[149];
   mHH[152]=1./2.*mHH[97] + 13./16.*mHH[95] + mHH[96];
   mHH[152]=mHH[14]*mHH[152];
   mHH[106]=mHH[106] - 1./2.*mHH[36] - 3./16. - mHH[37];
   mHH[106]=mHH[9]*mHH[106];
   mHH[155]=1./4.*mHH[30];
   mHH[158]=3./2.*mHH[71];
   mHH[159]=1./2.*mHH[66];
   mHH[101]=1./2.*mHH[139] + mHH[101] + mHH[143] + 3./4.*mHH[106] + 3./
   2.*mHH[152] + mHH[131] + mHH[113] + mHH[104] + 1./4.*mHH[127] + 5./8.
   *mHH[463] + mHH[38] + 1./2.*mHH[117] + 81./8.*mHH[6] + 3./32.*
   mHH[149] + mHH[148] + mHH[155] + mHH[147] + mHH[353] - mHH[15] - 33./
   4.*mHH[32] - 33./8.*mHH[33] - 21./64.*mHH[82] - 17./8.*mHH[16] - 37./
   8.*mHH[40] - 35./32.*mHH[77] + 3./16.*mHH[68] + 207./64.*mHH[41] - 
   63./64.*mHH[73] + mHH[158] - 13./4.*mHH[81] + 85./4.*mHH[76] + 81./8.
   *mHH[39] + 85./2.*mHH[78] + mHH[220] + mHH[159] - 13./12.*mHH[75] - 
   3./8.*mHH[18] + mHH[69] + 29./8.*mHH[70] + 1./2.*mHH[146] + mHH[20];
   mHH[101]=mHH[92]*mHH[101];
   mHH[104]=47./4.*mHH[33];
   mHH[106]=mHH[132] + 31 + mHH[104];
   mHH[106]=1./3.*mHH[106] + mHH[15];
   mHH[106]=mHH[243] + mHH[242] + mHH[241] + mHH[245] + mHH[240] + 1./2.
   *mHH[106] + mHH[244];
   mHH[106]=mHH[247] + 1./2.*mHH[106] + mHH[246];
   mHH[106]=mHH[8]*mHH[106];
   mHH[113]=117 + 47./3.*mHH[33];
   mHH[113]=mHH[251] + mHH[250] + mHH[15] + 1./4.*mHH[113] + mHH[249];
   mHH[113]=mHH[243] + mHH[253] + mHH[241] + 1./2.*mHH[113] + mHH[252];
   mHH[113]=mHH[259] + mHH[258] + mHH[257] + 1./4.*mHH[113] + mHH[256];
   mHH[113]=mHH[7]*mHH[113];
   mHH[117]=mHH[261] + mHH[260] + 15*mHH[67] + 61./4. + mHH[262];
   mHH[127]= - 3./4.*mHH[31] + 3./4.*mHH[15] + 13./4.*mHH[32] + 7 + 47./
   16.*mHH[33];
   mHH[127]= - 17./8.*mHH[34] - 9./4.*mHH[30] + 3*mHH[127] + 17./8.*
   mHH[10];
   mHH[127]=mHH[6]*mHH[127];
   mHH[106]=mHH[113] + mHH[277] + mHH[106] + mHH[266] + mHH[280] + 1./2.
   *mHH[127] + mHH[276] + mHH[282] + mHH[281] + mHH[275] - 119./32.*
   mHH[16] + mHH[274] + mHH[273] + mHH[272] + mHH[271] - 3./8.*mHH[78]
    + mHH[270] + 1./8.*mHH[117] + mHH[278];
   mHH[106]=mHH[98]*mHH[106];
   mHH[113]=mHH[115] + mHH[130] - 5./3.*mHH[30] + mHH[129] + 1./2.*
   mHH[198] + mHH[398];
   mHH[117]=1./3.*mHH[235];
   mHH[113]=1./2.*mHH[113] + mHH[117];
   mHH[113]=mHH[91]*mHH[113];
   mHH[127]=mHH[8]*mHH[113];
   mHH[113]=mHH[7]*mHH[113];
   mHH[131]=mHH[6]*mHH[211];
   mHH[131]=mHH[131] + mHH[334];
   mHH[139]=mHH[4]*mHH[6]*mHH[215];
   mHH[131]=3./4.*mHH[131] + mHH[139];
   mHH[131]=mHH[91]*mHH[131];
   mHH[139]=MMH*mHH[98]*mHH[268];
   mHH[106]=1./2.*mHH[139] + mHH[106] + 1./2.*mHH[113] + 3./2.*mHH[131]
    + mHH[127];
   mHH[106]=MMH*mHH[106];
   mHH[113]=47*mHH[33];
   mHH[127]= - 59 + mHH[113];
   mHH[127]=1./4.*mHH[127] + mHH[132];
   mHH[127]=1./3.*mHH[127] + mHH[15];
   mHH[127]=mHH[474] - mHH[8] + mHH[186] + mHH[473] + mHH[243] + 
   mHH[241] + mHH[245] + mHH[240] + 1./2.*mHH[127] + mHH[349];
   mHH[127]=mHH[13]*mHH[127];
   mHH[104]=mHH[467] + mHH[176] + mHH[114] + mHH[174] + mHH[153] + 
   mHH[151] + mHH[132] + 103./3. + mHH[104];
   mHH[104]=mHH[11]*mHH[104];
   mHH[113]=1445./2. + mHH[113];
   mHH[113]=1./4.*mHH[113] + mHH[132];
   mHH[113]=mHH[482] - mHH[38] + mHH[481] + mHH[480] + mHH[251] + 
   mHH[31] + 1./3.*mHH[113] + mHH[15];
   mHH[110]=mHH[110] + mHH[483] + mHH[258] + mHH[206] + 1./4.*mHH[113]
    + mHH[188];
   mHH[110]=mHH[12]*mHH[110];
   mHH[113]= - 143./16. + 41./3.*mHH[8];
   mHH[113]=1./3.*mHH[113] + mHH[277];
   mHH[113]=mHH[14]*mHH[113];
   mHH[104]=mHH[326] + mHH[110] + mHH[493] + mHH[113] + mHH[127] + 
   mHH[492] + mHH[187] + 1./4.*mHH[104] - 1495./96.*mHH[17] + 91./24.*
   mHH[43] + mHH[490] + mHH[489] + mHH[488] - mHH[80] + mHH[487] + 
   mHH[486] + mHH[306] + mHH[485] - mHH[86];
   mHH[104]=mHH[98]*mHH[104];
   mHH[110]=55./8.*mHH[32];
   mHH[113]=1./2.*mHH[15];
   mHH[127]= - 5./12.*mHH[10];
   mHH[131]=5./12.*mHH[34];
   mHH[139]=mHH[243] + mHH[242] + mHH[131] + mHH[127] + mHH[244] + 
   mHH[113] + mHH[110] - 9 - 55./8.*mHH[33];
   mHH[143]=mHH[204] + mHH[10];
   mHH[143]=mHH[4]*mHH[143];
   mHH[139]= - 13./16.*mHH[8] + 1./2.*mHH[139] + 1./3.*mHH[143];
   mHH[139]=mHH[8]*mHH[139];
   mHH[143]= - 55*mHH[33];
   mHH[146]=mHH[199] - 147./2. + mHH[143];
   mHH[146]=mHH[254] + 1./4.*mHH[146] + mHH[15];
   mHH[149]=1./3.*mHH[30];
   mHH[146]=mHH[243] + mHH[242] + mHH[131] + 1./2.*mHH[146] + mHH[149];
   mHH[152]=mHH[208] + 1 + mHH[255];
   mHH[152]=mHH[4]*mHH[152];
   mHH[108]=3./16.*mHH[108];
   mHH[146]=mHH[287] + mHH[108] - 5./8.*mHH[8] + 1./4.*mHH[146] + 1./3.
   *mHH[152];
   mHH[146]=mHH[7]*mHH[146];
   mHH[152]= - mHH[31] + mHH[15] + 55./4.*mHH[32] - 9 - 55./4.*mHH[33];
   mHH[162]=5./2.*mHH[34];
   mHH[152]=mHH[162] + mHH[114] + 3*mHH[152] + mHH[265];
   mHH[152]=mHH[6]*mHH[152];
   mHH[145]=mHH[301] + mHH[145] + mHH[55] + mHH[195];
   mHH[145]=MMH*mHH[145];
   mHH[163]=21./2.*mHH[63] - 13*mHH[67] - 235./4. + mHH[262];
   mHH[163]=1./2.*mHH[163] - 13*mHH[70];
   mHH[164]= - mHH[15] + mHH[10];
   mHH[165]=mHH[6]*mHH[164];
   mHH[165]=5./2. + 3*mHH[165];
   mHH[165]=mHH[4]*mHH[165];
   mHH[166]=3./16.*mHH[183];
   mHH[139]=1./2.*mHH[145] + mHH[146] + mHH[166] + mHH[139] + 1./2.*
   mHH[165] + mHH[280] + 3./8.*mHH[152] + 3./16.*mHH[322] + mHH[321] + 
   mHH[281] + 3./64.*mHH[82] + 37./64.*mHH[16] + 5./4.*mHH[40] + 3./64.
   *mHH[68] - 17./8.*mHH[76] - 63./16.*mHH[39] - 17./4.*mHH[78] + 
   mHH[66] + 1./4.*mHH[163] + mHH[278];
   mHH[139]=MMH*mHH[139];
   mHH[145]=55./2.*mHH[32] - 283 - 55./2.*mHH[33];
   mHH[145]=1./2.*mHH[145] + mHH[15];
   mHH[146]=mHH[210] + mHH[10];
   mHH[146]=mHH[4]*mHH[146];
   mHH[127]=mHH[474] + 49./6.*mHH[8] + 2./3.*mHH[146] + mHH[473] + 
   mHH[243] + mHH[131] + mHH[245] + mHH[127] + 1./2.*mHH[145] + 
   mHH[349];
   mHH[127]=mHH[13]*mHH[127];
   mHH[131]=165*mHH[32] - 89 - 165*mHH[33];
   mHH[131]=mHH[467] + mHH[162] + mHH[114] + mHH[265] + mHH[153] + 1./4.
   *mHH[131] + mHH[151];
   mHH[131]=mHH[11]*mHH[131];
   mHH[145]=11*mHH[32] - 457./4. + mHH[184];
   mHH[145]=mHH[254] - mHH[31] + 5./4.*mHH[145] + mHH[15];
   mHH[152]= - 1./4.*mHH[38];
   mHH[103]=mHH[320] - 39./32.*mHH[9] + mHH[108] + 53./12.*mHH[8] + 1./
   3.*mHH[146] + 1./16.*mHH[103] + mHH[152] + 5./24.*mHH[34] + 1./4.*
   mHH[145] - mHH[30];
   mHH[103]=mHH[12]*mHH[103];
   mHH[108]=mHH[160] - mHH[21];
   mHH[145]=mHH[11]*mHH[164];
   mHH[108]=2./3.*mHH[108] + mHH[145];
   mHH[108]=mHH[4]*mHH[108];
   mHH[145]=mHH[166] + 27./32. - 5./3.*mHH[8];
   mHH[145]=mHH[14]*mHH[145];
   mHH[146]= - mHH[13] + 5./16.*mHH[14];
   mHH[146]=mHH[9]*mHH[146];
   mHH[160]= - 41./24.*mHH[12] - 19./12.*mHH[14] - 1./8.*mHH[11] + 8*
   mHH[13];
   mHH[160]=mHH[7]*mHH[160];
   mHH[163]= - 1./2.*mHH[86];
   mHH[103]=mHH[139] + mHH[160] + mHH[103] + 3*mHH[146] + mHH[145] + 
   mHH[127] + 1./4.*mHH[491] + mHH[108] + 1./4.*mHH[131] + 461./192.*
   mHH[17] + 13./3.*mHH[43] + 63./32.*mHH[44] + 3./64.*mHH[26] + 
   mHH[347] - mHH[80] - 27./16.*mHH[42] + 1./4.*mHH[21] + 17./16.*
   mHH[25] + mHH[163] - 45./8.*mHH[28] + 171./8.*mHH[89] + mHH[484] + 
   mHH[22] - mHH[88];
   mHH[103]=mHH[92]*mHH[103];
   mHH[108]= - 47*mHH[33];
   mHH[127]=mHH[108] + mHH[136] - 557./3. + mHH[179];
   mHH[127]=mHH[431] + mHH[433] + mHH[105] + mHH[432] + mHH[429] + 
   mHH[190] + mHH[421] + mHH[389] + mHH[331] + 1./2.*mHH[127] + 
   mHH[314];
   mHH[127]=mHH[9]*mHH[127];
   mHH[127]=mHH[385] + mHH[127];
   mHH[127]=mHH[98]*mHH[127];
   mHH[131]= - 5*mHH[36] + 51 - 5*mHH[37];
   mHH[102]=mHH[102] + 1./2.*mHH[131] + 55*mHH[33];
   mHH[102]=2*mHH[31] + 1./2.*mHH[102] + mHH[236];
   mHH[131]= - 5*mHH[34];
   mHH[136]=8*mHH[235];
   mHH[102]= - 3./2.*mHH[9] + mHH[433] + mHH[136] + mHH[432] + mHH[131]
    + mHH[190] + 3*mHH[102] + mHH[284];
   mHH[102]=mHH[9]*mHH[102];
   mHH[139]= - 11 + 27./2.*mHH[74];
   mHH[102]=3./4.*mHH[380] + mHH[102] + mHH[440] + mHH[441] - 13./4.*
   mHH[82] + mHH[438] + mHH[138] + mHH[437] - 33./4.*mHH[41] - 3./4.*
   mHH[36] + mHH[436] + 3./2.*mHH[73] + mHH[435] + 1./2.*mHH[139] + 
   mHH[336];
   mHH[102]=mHH[92]*mHH[102];
   mHH[139]= - mHH[33] + mHH[32];
   mHH[145]=2*mHH[30];
   mHH[146]= - 2*mHH[38];
   mHH[160]=mHH[146] + mHH[34] + mHH[145] - mHH[10] + mHH[349] + 33./2.
   *mHH[139] + mHH[205];
   mHH[164]=mHH[4]*mHH[164];
   mHH[160]=3*mHH[160] + 8*mHH[164];
   mHH[160]=mHH[91]*mHH[160];
   mHH[165]=mHH[9]*mHH[160];
   mHH[166]= - 2*mHH[13];
   mHH[167]= - mHH[12] - 4*mHH[14] + mHH[166] + mHH[169] + mHH[444] + 
   mHH[87];
   mHH[167]=mHH[92]*mHH[167];
   mHH[168]=mHH[98]*mHH[445];
   mHH[167]=mHH[168] + 3*mHH[167];
   mHH[167]=mHH[1]*mHH[167];
   mHH[168]= - mHH[91]*mHH[38];
   mHH[169]=mHH[8]*mHH[168];
   mHH[170]=mHH[7]*mHH[168];
   mHH[102]=mHH[167] + mHH[102] + mHH[127] + 3./4.*mHH[170] + 3./2.*
   mHH[169] + mHH[165];
   mHH[102]=mHH[1]*mHH[102];
   mHH[127]=mHH[98]*mHH[447];
   mHH[165]=3*mHH[58] + 1./4.*mHH[60];
   mHH[167]=1./4.*mHH[448];
   mHH[171]=1./4.*mHH[450];
   mHH[165]=mHH[171] + 1./8.*mHH[449] + mHH[167] + 1./4.*mHH[442] + 1./
   2.*mHH[165] + mHH[54];
   mHH[165]=mHH[92]*mHH[165];
   mHH[172]=mHH[98]*mHH[451];
   mHH[173]=mHH[98]*mHH[453];
   mHH[175]=5 + mHH[81];
   mHH[175]=mHH[82] + 1./2.*mHH[175] + mHH[83];
   mHH[175]=mHH[92]*mHH[175];
   mHH[173]=mHH[173] + 3*mHH[175];
   mHH[173]=mHH[1]*mHH[173];
   mHH[175]= - mHH[58] - 1./2.*mHH[60];
   mHH[175]=mHH[92]*mHH[175];
   mHH[172]=mHH[173] + mHH[172] + 3*mHH[175];
   mHH[172]=MMt*mHH[1]*mHH[172];
   mHH[102]=mHH[172] + mHH[102] + mHH[127] + 3*mHH[165];
   mHH[102]=MMt*mHH[102];
   mHH[127]= - 59./2. - 3*mHH[37];
   mHH[127]=1./2.*mHH[127] - 11*mHH[36];
   mHH[108]=5*mHH[127] + mHH[108];
   mHH[105]=mHH[455] + mHH[105] + mHH[432] + mHH[429] + mHH[190] + 
   mHH[421] + mHH[389] + mHH[331] + 1./2.*mHH[108] + mHH[314];
   mHH[105]=mHH[14]*mHH[105];
   mHH[108]= - 133./2. - 32*mHH[34];
   mHH[127]= - 9./4.*mHH[38];
   mHH[108]=1./3.*mHH[108] + mHH[127];
   mHH[108]=mHH[13]*mHH[108];
   mHH[105]=mHH[466] + mHH[465] + mHH[105] + mHH[346] + mHH[108];
   mHH[105]=mHH[98]*mHH[105];
   mHH[108]= - 33*mHH[32] + 33*mHH[33] - 3./2.*mHH[36] + 19 - 3./2.*
   mHH[37];
   mHH[108]=mHH[455] + mHH[136] + mHH[432] + mHH[131] + mHH[190] + 
   mHH[284] + mHH[389] + 5./2.*mHH[108] + mHH[331];
   mHH[108]=mHH[14]*mHH[108];
   mHH[127]=mHH[127] - 7./2. + 8*mHH[34];
   mHH[127]=mHH[13]*mHH[127];
   mHH[131]=4*mHH[13] - 5./4.*mHH[14];
   mHH[131]=mHH[9]*mHH[131];
   mHH[136]=39./4.*mHH[9] + 1 + 23./4.*mHH[34];
   mHH[136]=mHH[12]*mHH[136];
   mHH[108]=3./4.*mHH[391] + 1./2.*mHH[136] + 3*mHH[131] + mHH[108] + 
   mHH[127] + 5./8.*mHH[17] + mHH[461] - 41./4.*mHH[44] + 57./16.*
   mHH[26] - 5./2.*mHH[87] + 6*mHH[29] + mHH[460];
   mHH[108]=mHH[92]*mHH[108];
   mHH[127]= - mHH[8]*mHH[91];
   mHH[131]=mHH[160] + 3./2.*mHH[127];
   mHH[131]=mHH[14]*mHH[131];
   mHH[136]=mHH[34] - mHH[38];
   mHH[160]=mHH[91]*mHH[136];
   mHH[165]=mHH[12]*mHH[160];
   mHH[172]=mHH[13]*mHH[160];
   mHH[173]= - mHH[14]*mHH[91];
   mHH[175]=mHH[7]*mHH[173];
   mHH[105]=mHH[108] + mHH[105] + 3./4.*mHH[175] + 3./4.*mHH[165] + 3./
   2.*mHH[172] + mHH[131];
   mHH[105]=mHH[1]*mHH[105];
   mHH[108]=541./3. + mHH[458];
   mHH[108]=47./2.*mHH[33] + 1./2.*mHH[108] + mHH[464];
   mHH[108]=mHH[467] + mHH[176] + mHH[114] + mHH[174] + mHH[153] + 
   mHH[151] + 1./2.*mHH[108] + mHH[132];
   mHH[108]=mHH[468] + mHH[359] + 1./2.*mHH[108] + mHH[189];
   mHH[108]=mHH[9]*mHH[108];
   mHH[131]=7 + 16*mHH[34];
   mHH[131]=2./9.*mHH[131] + mHH[330];
   mHH[131]=mHH[8]*mHH[131];
   mHH[108]=mHH[471] + mHH[479] + mHH[108] + mHH[477] + mHH[342] + 
   mHH[469] + mHH[131];
   mHH[108]=mHH[98]*mHH[108];
   mHH[131]= - 7*mHH[81] - 25./4. + 3*mHH[74];
   mHH[131]= - 15./4.*mHH[41] + 9./4.*mHH[73] + 1./2.*mHH[131] + 
   mHH[71];
   mHH[132]=1./8.*mHH[82];
   mHH[131]=3./4.*mHH[462] - mHH[38] + mHH[340] + mHH[383] + mHH[132]
    + mHH[292] + mHH[470] + mHH[138] - 3./4.*mHH[68] + 3./2.*mHH[131]
    + mHH[83];
   mHH[138]=mHH[199] + mHH[143] + mHH[36] - 51./2. + mHH[37];
   mHH[138]= - mHH[31] + 1./4.*mHH[138] + mHH[15];
   mHH[138]=mHH[467] + mHH[162] + mHH[114] + 3*mHH[138] + mHH[265];
   mHH[138]=mHH[359] + 1./2.*mHH[138] + 2*mHH[164];
   mHH[138]=mHH[9]*mHH[138];
   mHH[143]= - 7 - 4*mHH[34];
   mHH[143]=2./3.*mHH[143] + mHH[330];
   mHH[143]=mHH[8]*mHH[143];
   mHH[151]= - 19 - 29./2.*mHH[34];
   mHH[151]=1./3.*mHH[151] + mHH[367];
   mHH[151]=mHH[7]*mHH[151];
   mHH[153]= - mHH[96] - mHH[97];
   mHH[153]=mHH[14]*mHH[153];
   mHH[131]=1./4.*mHH[151] + 3./4.*mHH[478] + mHH[138] + 3./4.*mHH[153]
    + mHH[342] + 1./2.*mHH[131] + mHH[143];
   mHH[131]=mHH[92]*mHH[131];
   mHH[138]=mHH[9]*mHH[302];
   mHH[143]= - mHH[34] + mHH[38];
   mHH[143]=mHH[91]*mHH[143];
   mHH[151]=mHH[8]*mHH[143];
   mHH[153]=mHH[7]*mHH[143];
   mHH[162]=MMH*mHH[98]*mHH[475];
   mHH[102]=mHH[102] + mHH[105] + mHH[131] + mHH[162] + mHH[108] + 1./4.
   *mHH[153] + 1./2.*mHH[151] + mHH[138];
   mHH[102]=MMt*mHH[102];
   mHH[105]= - 5./2. + 3*mHH[65];
   mHH[105]=1./2.*mHH[105] - mHH[84];
   mHH[108]= - 5./4.*mHH[74];
   mHH[131]= - 9./2.*mHH[39];
   mHH[138]=1./2.*mHH[35];
   mHH[151]= - 9./2.*mHH[6];
   mHH[105]=1./2.*mHH[339] + mHH[151] - 1./4.*mHH[82] + mHH[470] + 
   mHH[375] + mHH[335] - mHH[73] + mHH[158] + mHH[81] + mHH[138] + 
   mHH[131] + 3*mHH[105] + mHH[108];
   mHH[105]=mHH[344] + mHH[343] + 3./2.*mHH[328] + 3*mHH[105] - mHH[38]
   ;
   mHH[153]=1./8.*mHH[37] + 1 + 7./4.*mHH[35];
   mHH[162]=19./4.*mHH[34];
   mHH[153]=mHH[359] + mHH[345] + mHH[162] + 3*mHH[153] + mHH[317];
   mHH[153]=mHH[9]*mHH[153];
   mHH[105]=mHH[379] + 1./4.*mHH[355] + 1./2.*mHH[105] + mHH[153];
   mHH[105]=mHH[92]*mHH[105];
   mHH[153]= - 7./2. + 9*mHH[65];
   mHH[108]=mHH[470] + mHH[375] + mHH[335] - mHH[73] + mHH[158] + 
   mHH[81] + mHH[138] + mHH[131] + mHH[108] + 1./2.*mHH[153] - 3*
   mHH[84];
   mHH[138]=mHH[433] + mHH[351] + mHH[150] - 7./12.*mHH[34] + mHH[317]
    + mHH[349] + 3./4.*mHH[37] - 10./3. + 21./2.*mHH[35];
   mHH[138]=mHH[9]*mHH[138];
   mHH[108]=mHH[138] + mHH[344] + mHH[343] + mHH[328] - mHH[38] + 
   mHH[340] + mHH[345] + 3*mHH[108] + 13./4.*mHH[82];
   mHH[138]=mHH[98]*mHH[108];
   mHH[150]=mHH[467] + mHH[381] + mHH[179] + 79./3. - 30*mHH[35];
   mHH[150]=mHH[14]*mHH[150];
   mHH[146]= - 9./4. + mHH[146];
   mHH[146]=mHH[13]*mHH[146];
   mHH[153]=mHH[354] - mHH[13];
   mHH[153]=5*mHH[153] - 2./3.*mHH[14];
   mHH[153]=mHH[9]*mHH[153];
   mHH[158]=47./4. - 18*mHH[9];
   mHH[158]=mHH[12]*mHH[158];
   mHH[146]=mHH[158] + mHH[153] + mHH[150] + 3*mHH[146] + mHH[350] + 
   mHH[348] - 41./4.*mHH[17];
   mHH[150]=mHH[98]*mHH[146];
   mHH[153]= - 9./2. - 5*mHH[38];
   mHH[153]=mHH[13]*mHH[153];
   mHH[153]=1./2.*mHH[153] - 9./4.*mHH[11] + mHH[325] - 3./4.*mHH[17];
   mHH[158]= - 15*mHH[35];
   mHH[165]= - 15./4.*mHH[37];
   mHH[174]= - 19*mHH[34];
   mHH[176]=mHH[174] + mHH[165] - 11./2. + mHH[158];
   mHH[176]=mHH[14]*mHH[176];
   mHH[178]= - 5*mHH[11] + mHH[13];
   mHH[178]=3./2.*mHH[178] - 19*mHH[14];
   mHH[178]=mHH[9]*mHH[178];
   mHH[181]=5./2. - mHH[38];
   mHH[181]=3./4.*mHH[181] - 13*mHH[9];
   mHH[181]=mHH[12]*mHH[181];
   mHH[153]=mHH[181] + mHH[178] + 3./2.*mHH[153] + mHH[176];
   mHH[153]=mHH[92]*mHH[153];
   mHH[176]=mHH[14]*mHH[143];
   mHH[178]=mHH[13]*mHH[168];
   mHH[176]=1./2.*mHH[178] + mHH[176];
   mHH[181]=mHH[9]*mHH[91];
   mHH[182]=3./2.*mHH[168] + 17*mHH[181];
   mHH[182]=mHH[12]*mHH[182];
   mHH[183]= - mHH[13]*mHH[91];
   mHH[184]=mHH[9]*mHH[183];
   mHH[150]=mHH[153] + mHH[150] + 1./2.*mHH[182] + 3*mHH[176] + 17./2.*
   mHH[184];
   mHH[150]=mHH[1]*mHH[150];
   mHH[153]=8./3. + mHH[158];
   mHH[153]= - 6*mHH[9] + mHH[467] + mHH[381] + 2*mHH[153] + mHH[179];
   mHH[153]=mHH[9]*mHH[153];
   mHH[176]= - mHH[82] - 5./4.*mHH[83] - 13./2.*mHH[41] - 1./2.*mHH[37]
    - 2*mHH[35] + 5./2.*mHH[74] - 5*mHH[84] + 3 + 4*mHH[72];
   mHH[153]=3*mHH[176] + mHH[153];
   mHH[176]=mHH[98]*mHH[153];
   mHH[158]=mHH[367] + mHH[174] + mHH[165] - 16 + mHH[158];
   mHH[158]=mHH[9]*mHH[158];
   mHH[158]=3*mHH[371] + mHH[158];
   mHH[158]=mHH[92]*mHH[158];
   mHH[165]=4*mHH[366] - 18*mHH[14] + mHH[166] - 4*mHH[11] + mHH[443]
    - 4*mHH[42] + 4*mHH[365] - mHH[29];
   mHH[166]=mHH[98]*mHH[165];
   mHH[174]=mHH[92]*mHH[337];
   mHH[166]=mHH[166] + mHH[174];
   mHH[166]=mHH[1]*mHH[166];
   mHH[143]=mHH[9]*mHH[143];
   mHH[143]=3*mHH[166] + mHH[158] + 3*mHH[143] + mHH[176];
   mHH[143]=mHH[1]*mHH[143];
   mHH[158]=1./4.*mHH[373] + 2*mHH[369] - mHH[58];
   mHH[166]=mHH[98]*mHH[158];
   mHH[174]=mHH[92]*mHH[374];
   mHH[166]=mHH[166] + mHH[174];
   mHH[174]= - mHH[98]*mHH[59];
   mHH[176]= - mHH[92]*mHH[59];
   mHH[174]=2*mHH[174] + mHH[176];
   mHH[176]=mHH[83] + 9 + 8*mHH[84];
   mHH[179]=mHH[98]*mHH[176];
   mHH[182]=mHH[92]*mHH[376];
   mHH[179]=mHH[179] + mHH[182];
   mHH[179]=mHH[1]*mHH[179];
   mHH[174]=2*mHH[174] + mHH[179];
   mHH[174]=MMt*mHH[1]*mHH[174];
   mHH[143]=3*mHH[174] + 3*mHH[166] + mHH[143];
   mHH[143]=MMt*mHH[143];
   mHH[166]=mHH[31] - mHH[30];
   mHH[174]=mHH[441] + 2*mHH[166] + 3./4.*mHH[34];
   mHH[174]=mHH[9]*mHH[91]*mHH[174];
   mHH[179]=mHH[91]*mHH[38];
   mHH[182]=mHH[8]*mHH[179];
   mHH[185]=mHH[7]*mHH[179];
   mHH[186]=MMH*mHH[98]*mHH[377];
   mHH[105]=mHH[143] + mHH[150] + mHH[105] + 3./2.*mHH[186] + mHH[138]
    + 1./4.*mHH[185] + 1./2.*mHH[182] + mHH[174];
   mHH[105]=MMt*mHH[105];
   mHH[138]=5./2.*mHH[8] + mHH[345] + mHH[317] + mHH[349] + 8 + 15./2.*
   mHH[35];
   mHH[138]=mHH[14]*mHH[138];
   mHH[143]=1 + mHH[267];
   mHH[143]=1./3.*mHH[143] + mHH[387];
   mHH[143]=mHH[11]*mHH[143];
   mHH[150]=35 - mHH[34];
   mHH[150]=7./9.*mHH[150] + mHH[38];
   mHH[150]=mHH[13]*mHH[150];
   mHH[174]=1./3.*mHH[14] + mHH[382] + 5./2.*mHH[13];
   mHH[174]=mHH[9]*mHH[174];
   mHH[182]= - 73 - mHH[34];
   mHH[182]=7./9.*mHH[182] + mHH[38];
   mHH[182]=1./4.*mHH[182] + 9*mHH[9];
   mHH[182]=mHH[12]*mHH[182];
   mHH[138]=1./2.*mHH[182] + 1./2.*mHH[174] + mHH[138] + 1./4.*mHH[150]
    + 1./4.*mHH[143] + 27./8.*mHH[17] + mHH[386] + 7*mHH[44] - 11./8.*
   mHH[26] + mHH[384] + 5./2.*mHH[87];
   mHH[143]=mHH[98]*mHH[138];
   mHH[150]= - 3*mHH[72];
   mHH[174]=29./4. + mHH[150];
   mHH[182]= - 3*mHH[65];
   mHH[185]=3*mHH[39];
   mHH[186]= - 1./2.*mHH[71];
   mHH[187]=3./4.*mHH[41];
   mHH[174]=mHH[187] + mHH[186] + mHH[185] + 1./4.*mHH[174] + mHH[182];
   mHH[119]=mHH[30] + mHH[119] + mHH[31];
   mHH[119]=mHH[291] + 1./2.*mHH[119] + mHH[177];
   mHH[119]=mHH[9]*mHH[119];
   mHH[188]= - 71 + mHH[267];
   mHH[188]=1./9.*mHH[188] + mHH[115];
   mHH[188]=1./2.*mHH[188] + mHH[327];
   mHH[188]=mHH[7]*mHH[188];
   mHH[189]=1 - 1./8.*mHH[34];
   mHH[190]=mHH[6]*mHH[189];
   mHH[189]=7./9.*mHH[189] + 1./8.*mHH[38];
   mHH[189]=mHH[8]*mHH[189];
   mHH[119]=1./4.*mHH[188] + mHH[119] + mHH[189] + 9./16.*mHH[334] + 7./
   2.*mHH[190] - 5./8.*mHH[82] + 11./8.*mHH[16] + mHH[333] + 3*mHH[174]
    - 11./8.*mHH[68];
   mHH[174]=mHH[98]*mHH[119];
   mHH[188]=mHH[6]*mHH[34];
   mHH[188]=mHH[188] + mHH[279];
   mHH[188]=mHH[91]*mHH[188];
   mHH[189]=mHH[8]*mHH[160];
   mHH[188]=9./2.*mHH[188] + mHH[189];
   mHH[189]=mHH[7]*mHH[160];
   mHH[190]= - mHH[31] + mHH[30];
   mHH[192]=mHH[91]*mHH[190];
   mHH[195]=mHH[9]*mHH[192];
   mHH[188]=1./8.*mHH[189] + 1./4.*mHH[188] + mHH[195];
   mHH[174]=1./2.*mHH[188] + mHH[174];
   mHH[174]=MMH*mHH[174];
   mHH[150]=37./4. + mHH[150];
   mHH[162]=7 + mHH[162];
   mHH[188]=mHH[6]*mHH[162];
   mHH[132]=1./2.*mHH[188] + mHH[132] + 1./8.*mHH[16] + 1./2.*mHH[40]
    - 1./8.*mHH[68] + mHH[187] + mHH[186] + mHH[185] + 1./4.*mHH[150]
    + mHH[182];
   mHH[150]=mHH[291] + mHH[177] + mHH[161] + mHH[30];
   mHH[150]=mHH[9]*mHH[150];
   mHH[161]=1 + mHH[34];
   mHH[182]=19./6.*mHH[161] + mHH[327];
   mHH[182]=mHH[7]*mHH[182];
   mHH[162]=mHH[8]*mHH[162];
   mHH[132]=1./4.*mHH[182] + mHH[150] + 3*mHH[132] + 1./3.*mHH[162];
   mHH[132]=MMH*mHH[132];
   mHH[150]=3 + 5./4.*mHH[35];
   mHH[150]=mHH[291] + mHH[352] + 3*mHH[150] + mHH[317];
   mHH[150]=mHH[14]*mHH[150];
   mHH[162]=3*mHH[44];
   mHH[182]=1./8.*mHH[17] - 1./4.*mHH[43] + mHH[162] - 1./8.*mHH[26] + 
   mHH[233] + mHH[217];
   mHH[185]=mHH[11]*mHH[161];
   mHH[187]=19*mHH[34];
   mHH[188]=119./2. + mHH[187];
   mHH[188]=mHH[13]*mHH[188];
   mHH[182]=1./6.*mHH[188] + 3*mHH[182] + 19./4.*mHH[185];
   mHH[185]=mHH[234] - 1./2.*mHH[13];
   mHH[185]=3*mHH[185] + mHH[283];
   mHH[185]=mHH[9]*mHH[185];
   mHH[187]=83./2. + mHH[187];
   mHH[187]=1./6.*mHH[187] + 13*mHH[9];
   mHH[187]=mHH[12]*mHH[187];
   mHH[132]=1./2.*mHH[132] + 1./4.*mHH[388] + 1./4.*mHH[187] + 1./4.*
   mHH[185] + 1./2.*mHH[182] + mHH[150];
   mHH[132]=mHH[92]*mHH[132];
   mHH[150]=mHH[358] - 53*mHH[13];
   mHH[150]=1./4.*mHH[150] + 43./3.*mHH[14];
   mHH[150]=mHH[14]*mHH[150];
   mHH[182]= - mHH[12]*mHH[14];
   mHH[150]=mHH[150] + 31./4.*mHH[182];
   mHH[185]=mHH[98]*mHH[150];
   mHH[187]= - mHH[11] - 3*mHH[13];
   mHH[187]=9./4.*mHH[187] - 23*mHH[14];
   mHH[187]=mHH[14]*mHH[187];
   mHH[187]=mHH[187] + 101./4.*mHH[182];
   mHH[187]=mHH[92]*mHH[187];
   mHH[188]=mHH[14]*mHH[183];
   mHH[196]=mHH[14]*mHH[91];
   mHH[199]=mHH[12]*mHH[196];
   mHH[185]=1./2.*mHH[187] + mHH[185] + 10*mHH[188] + 31./4.*mHH[199];
   mHH[185]=mHH[1]*mHH[185];
   mHH[136]=mHH[91]*mHH[11]*mHH[136];
   mHH[136]=3./2.*mHH[136] + mHH[172];
   mHH[172]=mHH[8]*mHH[91];
   mHH[187]=mHH[91]*mHH[166];
   mHH[188]=2*mHH[187] + 1./2.*mHH[172];
   mHH[188]=mHH[14]*mHH[188];
   mHH[199]= - mHH[9]*mHH[91];
   mHH[200]=mHH[160] + 17*mHH[199];
   mHH[200]=mHH[12]*mHH[200];
   mHH[201]=mHH[9]*mHH[214];
   mHH[204]=mHH[7]*mHH[196];
   mHH[105]=mHH[105] + mHH[185] + mHH[132] + mHH[174] + mHH[143] + 1./4.
   *mHH[204] + 1./8.*mHH[200] + 17./8.*mHH[201] + 1./4.*mHH[136] + 
   mHH[188];
   mHH[105]=MMt*mHH[105];
   mHH[132]=mHH[6]*mHH[166];
   mHH[136]=mHH[8]*mHH[166];
   mHH[143]=mHH[7]*mHH[166];
   mHH[136]=1./12.*mHH[143] + 1./6.*mHH[136] + 3./4.*mHH[132] - 1./2.*
   mHH[16] - mHH[40] + mHH[159] - 9./8. + mHH[69];
   mHH[136]=MMH*mHH[136];
   mHH[143]= - 1./6.*mHH[30];
   mHH[174]=mHH[143] + 3 + mHH[324];
   mHH[185]= - 7./12.*mHH[8] + mHH[174] - 3./8.*mHH[6];
   mHH[185]=mHH[13]*mHH[185];
   mHH[174]=1./6.*mHH[8] + mHH[174] + mHH[412];
   mHH[174]=mHH[12]*mHH[174];
   mHH[188]=3 + mHH[31];
   mHH[200]=mHH[188] - mHH[30];
   mHH[200]=mHH[11]*mHH[200];
   mHH[201]= - 5./6.*mHH[12] + mHH[11] - 1./6.*mHH[13];
   mHH[201]=mHH[7]*mHH[201];
   mHH[136]=1./2.*mHH[136] + 1./4.*mHH[201] + 1./2.*mHH[174] + mHH[185]
    + 1./2.*mHH[404] + 1./4.*mHH[200] + mHH[299] - mHH[43] + 3./4.*
   mHH[42] + 1./4.*mHH[25] + mHH[415] + mHH[163];
   mHH[136]=MMH*mHH[136];
   mHH[163]=1./8.*mHH[11] + mHH[418];
   mHH[163]=mHH[13]*mHH[163];
   mHH[174]=mHH[414] + mHH[11] + 1./6.*mHH[13];
   mHH[174]=mHH[12]*mHH[174];
   mHH[136]=1./2.*mHH[136] + mHH[163] + 1./4.*mHH[174];
   mHH[136]=mHH[92]*MMH*mHH[136];
   mHH[163]=mHH[455] + 9./2. + mHH[30];
   mHH[163]=mHH[14]*mHH[163];
   mHH[174]=1./4.*mHH[12]*mHH[9];
   mHH[163]=mHH[174] + 1./4.*mHH[366] - 9./2.*mHH[44] + mHH[163];
   mHH[185]=mHH[98]*mHH[163];
   mHH[200]= - 9./8.*mHH[6];
   mHH[201]=mHH[190] + mHH[200];
   mHH[201]=mHH[91]*mHH[201];
   mHH[201]=mHH[201] + 1./4.*mHH[127];
   mHH[201]=mHH[14]*mHH[201];
   mHH[204]=mHH[9]*mHH[187];
   mHH[206]= - mHH[9]*mHH[30];
   mHH[210]=mHH[98]*mHH[206];
   mHH[204]=mHH[204] + mHH[210];
   mHH[204]=MMH*mHH[204];
   mHH[181]=mHH[12]*mHH[181];
   mHH[175]=1./4.*mHH[204] + mHH[185] + 1./8.*mHH[175] + 1./4.*mHH[181]
    + mHH[201] + 1./4.*mHH[184];
   mHH[175]=MMH*mHH[175];
   mHH[181]= - 7./4.*mHH[8] + mHH[200] + mHH[30] + 3./2. - mHH[31];
   mHH[181]=mHH[14]*mHH[181];
   mHH[184]=MMH*mHH[9]*mHH[166];
   mHH[174]=1./4.*mHH[184] + 1./8.*mHH[391] + mHH[174] + 1./4.*mHH[420]
    + mHH[428] + mHH[181];
   mHH[174]=MMH*mHH[174];
   mHH[181]= - 3./4.*mHH[11] - mHH[13];
   mHH[181]=mHH[14]*mHH[181];
   mHH[174]=mHH[174] + mHH[181] + 5./4.*mHH[182];
   mHH[174]=mHH[92]*mHH[174];
   mHH[181]= - mHH[91]*mHH[11];
   mHH[184]=3./2.*mHH[181] + mHH[214];
   mHH[184]=mHH[14]*mHH[184];
   mHH[185]=mHH[12]*mHH[173];
   mHH[184]=mHH[184] + 5./2.*mHH[185];
   mHH[185]=mHH[11] - 9./2.*mHH[13];
   mHH[185]=mHH[14]*mHH[185];
   mHH[182]=mHH[185] + mHH[182];
   mHH[185]=mHH[98]*mHH[182];
   mHH[174]=mHH[174] + mHH[175] + 1./2.*mHH[184] + mHH[185];
   mHH[175]=mHH[307] + 3./4.*mHH[279];
   mHH[175]=3*mHH[175] + 1./2.*mHH[439];
   mHH[184]=mHH[9]*mHH[190];
   mHH[175]=1./8.*mHH[380] + 1./2.*mHH[175] + mHH[184];
   mHH[175]=MMH*mHH[175];
   mHH[184]= - mHH[11]*mHH[38];
   mHH[185]=mHH[338] + 1./4.*mHH[184];
   mHH[200]= - 9 - mHH[38];
   mHH[200]=mHH[13]*mHH[200];
   mHH[201]=mHH[13] - 3./2.*mHH[14];
   mHH[201]=mHH[9]*mHH[201];
   mHH[152]=mHH[152] - mHH[9];
   mHH[152]=mHH[12]*mHH[152];
   mHH[152]=mHH[175] + mHH[152] + mHH[201] + mHH[357] + 3*mHH[185] + 1./
   2.*mHH[200];
   mHH[152]=mHH[92]*mHH[152];
   mHH[175]=mHH[91]*mHH[279];
   mHH[169]=9./2.*mHH[175] + mHH[169];
   mHH[175]=mHH[9]*mHH[30];
   mHH[175]=9./2.*mHH[307] + mHH[175];
   mHH[185]=mHH[98]*mHH[175];
   mHH[169]=mHH[185] + 1./8.*mHH[170] + 1./4.*mHH[169] + mHH[195];
   mHH[169]=MMH*mHH[169];
   mHH[170]= - mHH[14] + mHH[338] - 3./2.*mHH[13];
   mHH[185]=mHH[9]*mHH[11];
   mHH[170]=mHH[430] + 9*mHH[170] + mHH[185];
   mHH[185]=mHH[98]*mHH[170];
   mHH[184]=mHH[91]*mHH[184];
   mHH[178]=3./2.*mHH[184] + mHH[178];
   mHH[173]=mHH[214] + 3./2.*mHH[173];
   mHH[173]=mHH[9]*mHH[173];
   mHH[184]=1./4.*mHH[168] + mHH[199];
   mHH[184]=mHH[12]*mHH[184];
   mHH[152]=mHH[152] + mHH[169] + mHH[185] + mHH[184] + 1./2.*mHH[178]
    + mHH[173];
   mHH[169]=mHH[372] + mHH[424];
   mHH[169]=mHH[92]*mHH[169];
   mHH[168]=mHH[9]*mHH[168];
   mHH[173]=mHH[98]*mHH[372];
   mHH[168]=mHH[169] + mHH[168] + 3*mHH[173];
   mHH[169]=mHH[14]*mHH[38];
   mHH[169]=mHH[169] + mHH[304];
   mHH[169]=mHH[92]*mHH[169];
   mHH[173]=mHH[14]*mHH[179];
   mHH[178]=mHH[9]*mHH[196];
   mHH[169]=mHH[169] + mHH[173] + mHH[178];
   mHH[169]=mHH[1]*mHH[169];
   mHH[173]=mHH[9]*mHH[179];
   mHH[178]=mHH[92]*mHH[419];
   mHH[173]=mHH[173] + mHH[178];
   mHH[173]=MMt*mHH[1]*mHH[173];
   mHH[168]=mHH[173] + 1./4.*mHH[168] + mHH[169];
   mHH[168]=MMt*mHH[168];
   mHH[169]=mHH[263]*mHH[91];
   mHH[173]=mHH[92]*mHH[263];
   mHH[169]=mHH[169] + mHH[173];
   mHH[169]=mHH[1]*mHH[169];
   mHH[152]=3*mHH[168] + 1./2.*mHH[152] + 3*mHH[169];
   mHH[152]=MMt*mHH[152];
   mHH[152]=1./2.*mHH[174] + mHH[152];
   mHH[152]=MMt*mHH[152];
   mHH[168]=mHH[218] + 3*mHH[415] - mHH[86];
   mHH[143]=mHH[359] + 9 + mHH[143];
   mHH[143]=mHH[13]*mHH[143];
   mHH[169]=1./12.*mHH[8] + mHH[410] + 3 + mHH[282];
   mHH[169]=mHH[12]*mHH[169];
   mHH[173]= - 3./16.*mHH[6] + 1 - 1./8.*mHH[30];
   mHH[173]=mHH[11]*mHH[173];
   mHH[174]=mHH[11] - mHH[12];
   mHH[174]=mHH[7]*mHH[174];
   mHH[143]=11./48.*mHH[174] + 1./2.*mHH[169] + 1./2.*mHH[143] + 17./24.
   *mHH[404] + mHH[173] + mHH[299] - 3./2.*mHH[43] + 1./2.*mHH[168] + 
   mHH[42];
   mHH[168]=mHH[98]*mHH[143];
   mHH[169]=1./3.*mHH[166];
   mHH[173]=mHH[169] + mHH[407];
   mHH[173]=mHH[91]*mHH[173];
   mHH[173]=mHH[173] + 1./6.*mHH[127];
   mHH[173]=mHH[13]*mHH[173];
   mHH[174]=mHH[169] + 3./2.*mHH[6];
   mHH[174]=mHH[91]*mHH[174];
   mHH[174]=mHH[174] + 1./3.*mHH[172];
   mHH[174]=mHH[12]*mHH[174];
   mHH[178]=mHH[12]*mHH[91];
   mHH[179]=mHH[183] + mHH[178];
   mHH[179]=mHH[7]*mHH[179];
   mHH[166]=mHH[91]*mHH[11]*mHH[166];
   mHH[166]=1./12.*mHH[179] + 1./2.*mHH[174] + 1./2.*mHH[166] + 
   mHH[173];
   mHH[173]= - 1 + mHH[69];
   mHH[174]= - mHH[6]*mHH[30];
   mHH[184]= - mHH[8]*mHH[30];
   mHH[185]= - mHH[7]*mHH[30];
   mHH[173]=1./12.*mHH[185] + 1./6.*mHH[184] + 3./4.*mHH[174] - mHH[16]
    - 3*mHH[40] + 3*mHH[173] + mHH[66];
   mHH[174]=mHH[98]*mHH[173];
   mHH[132]=mHH[91]*mHH[132];
   mHH[184]=mHH[8]*mHH[187];
   mHH[185]=mHH[7]*mHH[187];
   mHH[132]=1./6.*mHH[185] + 3./2.*mHH[132] + 1./3.*mHH[184];
   mHH[132]=1./2.*mHH[132] + mHH[174];
   mHH[132]=MMH*mHH[132];
   mHH[132]=1./4.*mHH[132] + 1./4.*mHH[166] + mHH[168];
   mHH[132]=MMH*mHH[132];
   mHH[166]=mHH[91]*mHH[11];
   mHH[168]=1./3.*mHH[178] + mHH[166] + 1./3.*mHH[214];
   mHH[168]=mHH[12]*mHH[168];
   mHH[174]=1./2.*mHH[181] + 1./3.*mHH[183];
   mHH[174]=mHH[13]*mHH[174];
   mHH[168]=mHH[174] + 1./2.*mHH[168];
   mHH[174]=mHH[416] - 3./8.*mHH[13];
   mHH[174]=mHH[13]*mHH[174];
   mHH[178]= - 5./8.*mHH[12] + mHH[11] + 1./4.*mHH[13];
   mHH[178]=mHH[12]*mHH[178];
   mHH[174]=1./6.*mHH[178] - 1./16.*mHH[397] + mHH[174];
   mHH[178]=mHH[98]*mHH[174];
   mHH[132]=1./2.*mHH[132] + 1./8.*mHH[168] + mHH[178];
   mHH[132]=MMH*mHH[132];
   mHH[132]=mHH[152] + mHH[132] + 1./2.*mHH[136];
   mHH[132]=mHH[2]*mHH[132];
   mHH[136]=mHH[17] - 41*mHH[43] + mHH[393] + mHH[207];
   mHH[152]=mHH[30] - 71./8. + mHH[31];
   mHH[152]=1./4.*mHH[152] + mHH[177];
   mHH[152]=mHH[11]*mHH[152];
   mHH[168]=3./2. + mHH[244];
   mHH[168]= - 11./2.*mHH[8] + 15./8.*mHH[6] + 11./2.*mHH[168] + 
   mHH[149];
   mHH[168]=mHH[13]*mHH[168];
   mHH[178]=5./2. - 13*mHH[6];
   mHH[178]=1./2.*mHH[178] - 13./9.*mHH[8];
   mHH[178]=mHH[14]*mHH[178];
   mHH[184]=11./3.*mHH[8] + 21./4.*mHH[6] + 11./6.*mHH[30] - 5./2. + 
   mHH[244];
   mHH[184]=mHH[12]*mHH[184];
   mHH[185]=23./8.*mHH[12] + mHH[238] + mHH[11] - 19./8.*mHH[13];
   mHH[185]=mHH[7]*mHH[185];
   mHH[136]=1./6.*mHH[185] + 1./4.*mHH[184] + 1./2.*mHH[178] + 1./2.*
   mHH[168] + 1./3.*mHH[404] + 1./8.*mHH[136] + mHH[152];
   mHH[152]=mHH[98]*mHH[136];
   mHH[168]=mHH[364] + mHH[245];
   mHH[178]= - 1./8.*mHH[7];
   mHH[168]=mHH[178] + mHH[191] + 1./3.*mHH[168] + mHH[332];
   mHH[168]=mHH[7]*mHH[168];
   mHH[184]=mHH[392] + mHH[180] + 5 + mHH[353];
   mHH[184]=mHH[6]*mHH[184];
   mHH[185]=mHH[180] - 1 + mHH[406];
   mHH[185]=1./3.*mHH[185] + mHH[332];
   mHH[185]=mHH[8]*mHH[185];
   mHH[168]=1./4.*mHH[168] + 1./2.*mHH[185] + mHH[237] + 3./4.*mHH[184]
   ;
   mHH[184]=mHH[98]*mHH[168];
   mHH[185]=mHH[91]*mHH[6]*mHH[190];
   mHH[187]=mHH[8]*mHH[192];
   mHH[192]=mHH[7]*mHH[192];
   mHH[185]=1./6.*mHH[192] + 3./4.*mHH[185] + 1./3.*mHH[187];
   mHH[187]=MMH*mHH[98]*mHH[395];
   mHH[184]=1./4.*mHH[187] + 1./2.*mHH[185] + mHH[184];
   mHH[184]=MMH*mHH[184];
   mHH[185]=mHH[169] + 51./8.*mHH[6];
   mHH[185]=mHH[91]*mHH[185];
   mHH[185]=mHH[185] + 7./4.*mHH[172];
   mHH[185]=mHH[13]*mHH[185];
   mHH[169]=mHH[169] - 51./4.*mHH[6];
   mHH[169]=mHH[91]*mHH[169];
   mHH[169]=mHH[169] + 7./2.*mHH[127];
   mHH[169]=mHH[12]*mHH[169];
   mHH[187]=mHH[91]*mHH[11]*mHH[190];
   mHH[169]=7./8.*mHH[219] + 1./2.*mHH[169] + mHH[187] + mHH[185];
   mHH[152]=mHH[184] + 1./4.*mHH[169] + mHH[152];
   mHH[152]=MMH*mHH[152];
   mHH[169]=MMH*mHH[173];
   mHH[143]=mHH[143] + 1./4.*mHH[169];
   mHH[143]=MMH*mHH[143];
   mHH[143]=mHH[174] + 1./2.*mHH[143];
   mHH[143]=MMH*mHH[143];
   mHH[169]=MMH*mHH[206];
   mHH[163]=mHH[163] + 1./4.*mHH[169];
   mHH[163]=MMH*mHH[163];
   mHH[169]=MMH*mHH[175];
   mHH[173]=MMt*mHH[372];
   mHH[169]=9./2.*mHH[173] + mHH[170] + mHH[169];
   mHH[169]=MMt*mHH[169];
   mHH[163]=mHH[169] + mHH[182] + mHH[163];
   mHH[163]=MMt*mHH[163];
   mHH[143]=mHH[143] + 1./2.*mHH[163];
   mHH[143]=mHH[3]*mHH[143];
   mHH[163]=mHH[323] + mHH[180];
   mHH[163]=mHH[178] + mHH[191] + 1./3.*mHH[163] + mHH[332];
   mHH[163]=mHH[7]*mHH[163];
   mHH[149]=mHH[332] + mHH[149] - 1./3. + mHH[406];
   mHH[149]=mHH[8]*mHH[149];
   mHH[169]=mHH[392] + 5 + mHH[30];
   mHH[169]=mHH[6]*mHH[169];
   mHH[149]=mHH[396] + 1./4.*mHH[163] + 1./2.*mHH[149] + mHH[237] + 3./
   4.*mHH[169];
   mHH[149]=MMH*mHH[149];
   mHH[162]=mHH[393] + mHH[162];
   mHH[163]= - 71./16. + mHH[30];
   mHH[163]=1./2.*mHH[163] + mHH[177];
   mHH[163]=mHH[11]*mHH[163];
   mHH[162]=1./4.*mHH[404] + mHH[163] - 9./2.*mHH[17] + 1./8.*mHH[162]
    - 9*mHH[43];
   mHH[163]= - 27./16.*mHH[8] - 9./32.*mHH[6] + 1./6.*mHH[30] + 4 + 
   mHH[353];
   mHH[163]=mHH[13]*mHH[163];
   mHH[169]= - 1./2. + 5*mHH[6];
   mHH[169]=3./2.*mHH[169] + 5./3.*mHH[8];
   mHH[169]=mHH[14]*mHH[169];
   mHH[170]= - 73./48.*mHH[8] + 69./32.*mHH[6] + 5./24.*mHH[30] + 2 + 1.
   /8.*mHH[31];
   mHH[170]=mHH[12]*mHH[170];
   mHH[173]=37./6.*mHH[12] + mHH[402] + mHH[11] - 191./6.*mHH[13];
   mHH[173]=mHH[7]*mHH[173];
   mHH[149]=1./2.*mHH[149] + 1./16.*mHH[173] + mHH[170] + 1./4.*
   mHH[169] + 1./2.*mHH[162] + mHH[163];
   mHH[149]=MMH*mHH[149];
   mHH[162]=11./4.*mHH[11] - mHH[13];
   mHH[162]=mHH[13]*mHH[162];
   mHH[162]=165./16.*mHH[397] + mHH[162];
   mHH[163]=mHH[194] + mHH[411];
   mHH[163]=5*mHH[163] - 3./8.*mHH[14];
   mHH[163]=mHH[14]*mHH[163];
   mHH[169]=23./6.*mHH[12] + 29./6.*mHH[14] + 15*mHH[11] + 49./6.*
   mHH[13];
   mHH[169]=mHH[12]*mHH[169];
   mHH[162]=1./4.*mHH[169] + 1./2.*mHH[162] + mHH[163];
   mHH[149]=1./2.*mHH[162] + mHH[149];
   mHH[149]=mHH[92]*mHH[149];
   mHH[162]=11./3.*mHH[209] + 17*mHH[181] + 11./3.*mHH[183];
   mHH[162]=mHH[12]*mHH[162];
   mHH[163]=17./2.*mHH[166] + 11./3.*mHH[214];
   mHH[163]=mHH[13]*mHH[163];
   mHH[162]=mHH[163] + 1./2.*mHH[162];
   mHH[163]=17*mHH[11] - 65./3.*mHH[13];
   mHH[163]=mHH[13]*mHH[163];
   mHH[163]=165./4.*mHH[397] + mHH[163];
   mHH[166]=mHH[405] + mHH[418];
   mHH[166]=13./3.*mHH[166] + 5./8.*mHH[14];
   mHH[166]=mHH[14]*mHH[166];
   mHH[169]=13./2.*mHH[11];
   mHH[170]= - 97./18.*mHH[14] + mHH[169] + 53./3.*mHH[13];
   mHH[170]=1./2.*mHH[170] + 1./3.*mHH[12];
   mHH[170]=mHH[12]*mHH[170];
   mHH[163]=1./2.*mHH[170] + 1./8.*mHH[163] + mHH[166];
   mHH[166]=mHH[98]*mHH[163];
   mHH[105]=mHH[132] + mHH[143] + mHH[105] + mHH[149] + mHH[152] + 1./8.
   *mHH[162] + mHH[166];
   mHH[105]=mHH[2]*mHH[105];
   mHH[132]=mHH[1]*mHH[146];
   mHH[143]=mHH[1]*mHH[165];
   mHH[143]=mHH[153] + 3*mHH[143];
   mHH[143]=mHH[1]*mHH[143];
   mHH[146]=mHH[1]*mHH[176];
   mHH[146]= - 4*mHH[59] + mHH[146];
   mHH[146]=MMt*mHH[1]*mHH[146];
   mHH[143]=3*mHH[146] + 3*mHH[158] + mHH[143];
   mHH[143]=MMt*mHH[143];
   mHH[108]=mHH[143] + mHH[132] + mHH[108] + 3./2.*mHH[378];
   mHH[108]=MMt*mHH[108];
   mHH[119]=MMH*mHH[119];
   mHH[132]=mHH[1]*mHH[150];
   mHH[108]=mHH[108] + mHH[132] + mHH[138] + mHH[119];
   mHH[108]=MMt*mHH[108];
   mHH[119]=mHH[307] + 3./8.*mHH[334];
   mHH[119]=mHH[408] + mHH[403] + 3*mHH[119] + 1./4.*mHH[328];
   mHH[119]=MMH*mHH[119];
   mHH[132]=mHH[372] + 1./2.*mHH[419];
   mHH[132]=mHH[425] + 1./2.*mHH[132] + mHH[225];
   mHH[132]=MMt*mHH[132];
   mHH[138]=mHH[338] + 1./8.*mHH[417];
   mHH[115]= - 9 + mHH[115];
   mHH[115]=mHH[13]*mHH[115];
   mHH[143]= - mHH[1]*mHH[263];
   mHH[115]=3*mHH[132] + 3*mHH[143] + 1./2.*mHH[119] + 1./8.*mHH[423]
    + 1./2.*mHH[422] + mHH[357] + 3*mHH[138] + 1./2.*mHH[115];
   mHH[115]=MMt*mHH[115];
   mHH[119]= - 11./4.*mHH[8] + mHH[188] + mHH[356];
   mHH[119]=mHH[14]*mHH[119];
   mHH[119]=mHH[427] + mHH[426] + mHH[413] - 3*mHH[44] + mHH[119];
   mHH[119]=MMH*mHH[119];
   mHH[132]=mHH[194] - mHH[13];
   mHH[132]=mHH[14]*mHH[132];
   mHH[132]=7*mHH[132] + 1./2.*mHH[368];
   mHH[119]=1./2.*mHH[132] + mHH[119];
   mHH[115]=1./2.*mHH[119] + mHH[115];
   mHH[115]=MMt*mHH[115];
   mHH[119]=1./24.*mHH[401] + 1./12.*mHH[400] + 3./8.*mHH[390] - 1./4.*
   mHH[16] - mHH[40] + 1./4.*mHH[66] - 15./16. + mHH[69];
   mHH[119]=MMH*mHH[119];
   mHH[132]= - 3./2.*mHH[6] + 5 - mHH[31];
   mHH[132]=mHH[11]*mHH[132];
   mHH[138]= - 11./24.*mHH[8] + 3./16.*mHH[6] + 3 + mHH[281];
   mHH[138]=mHH[13]*mHH[138];
   mHH[143]=mHH[12]*mHH[409];
   mHH[146]=5*mHH[11] + mHH[13];
   mHH[146]=1./6.*mHH[146] - mHH[12];
   mHH[146]=mHH[7]*mHH[146];
   mHH[119]=1./2.*mHH[119] + 1./8.*mHH[146] + 1./4.*mHH[143] + mHH[138]
    + 11./24.*mHH[404] + 1./8.*mHH[132] - 1./4.*mHH[17] - mHH[43] + 5./
   8.*mHH[42] + 1./8.*mHH[25] + mHH[415] - 1./4.*mHH[86];
   mHH[119]=MMH*mHH[119];
   mHH[132]=mHH[169] - 5*mHH[13];
   mHH[132]=mHH[13]*mHH[132];
   mHH[138]=mHH[11] + mHH[312];
   mHH[138]=1./3.*mHH[138] - 1./2.*mHH[12];
   mHH[138]=mHH[12]*mHH[138];
   mHH[132]=mHH[138] - 1./2.*mHH[397] + 1./3.*mHH[132];
   mHH[119]=1./4.*mHH[132] + mHH[119];
   mHH[119]=MMH*mHH[119];
   mHH[115]=1./2.*mHH[119] + mHH[115];
   mHH[115]=mHH[3]*mHH[115];
   mHH[119]=mHH[168] + mHH[396];
   mHH[119]=MMH*mHH[119];
   mHH[119]=mHH[136] + mHH[119];
   mHH[119]=MMH*mHH[119];
   mHH[108]=mHH[115] + mHH[108] + mHH[163] + mHH[119];
   mHH[108]=mHH[3]*mHH[108];
   mHH[115]=mHH[11] - 199*mHH[13];
   mHH[115]=mHH[13]*mHH[115];
   mHH[119]=39*mHH[13] - 7./2.*mHH[14];
   mHH[119]=mHH[14]*mHH[119];
   mHH[132]=149./6.*mHH[12] + 39*mHH[14] + 1./6.*mHH[11] + 17*mHH[13];
   mHH[132]=mHH[12]*mHH[132];
   mHH[115]=1./2.*mHH[132] + 1./6.*mHH[115] + mHH[119];
   mHH[115]=mHH[92]*mHH[115];
   mHH[119]=mHH[183] + mHH[209];
   mHH[119]=mHH[12]*mHH[119];
   mHH[132]=mHH[227]*mHH[91];
   mHH[119]=mHH[132] + 1./2.*mHH[119];
   mHH[132]=mHH[11] + 233*mHH[13];
   mHH[132]=mHH[13]*mHH[132];
   mHH[136]= - 107*mHH[13] - 185*mHH[14];
   mHH[136]=mHH[14]*mHH[136];
   mHH[132]=1./2.*mHH[132] + mHH[136];
   mHH[132]=1./12.*mHH[132] + mHH[239];
   mHH[132]=mHH[98]*mHH[132];
   mHH[115]=1./4.*mHH[115] + 17./4.*mHH[119] + mHH[132];
   mHH[115]=mHH[1]*mHH[115];
   mHH[119]=mHH[11]*mHH[212];
   mHH[132]=mHH[4]*mHH[11]*mHH[215];
   mHH[119]=3./4.*mHH[119] + mHH[132];
   mHH[119]=mHH[91]*mHH[119];
   mHH[132]=mHH[38] + mHH[383] + mHH[30] + mHH[208] + mHH[198] - 
   mHH[31];
   mHH[136]=1./2.*mHH[132] + 2./3.*mHH[235];
   mHH[136]=mHH[91]*mHH[136];
   mHH[127]=mHH[136] + 5./3.*mHH[127];
   mHH[127]=mHH[13]*mHH[127];
   mHH[117]=1./4.*mHH[132] + mHH[117];
   mHH[117]=mHH[91]*mHH[117];
   mHH[117]=mHH[117] + 5./3.*mHH[172];
   mHH[117]=mHH[12]*mHH[117];
   mHH[102]=mHH[105] + mHH[108] + mHH[102] + mHH[115] + mHH[103] + 
   mHH[106] + mHH[104] + 5./6.*mHH[179] + mHH[117] + mHH[119] + 
   mHH[127];
   mHH[102]=mHH[2]*mHH[102];
   mHH[103]= - 7./2. + 11*mHH[67];
   mHH[103]=1./4.*mHH[103] + mHH[107];
   mHH[103]=3*mHH[103] + mHH[118];
   mHH[103]=mHH[126] + mHH[130] + mHH[128] + mHH[129] + mHH[476] - 57./
   8.*mHH[16] + mHH[125] + mHH[124] - 99./16.*mHH[76] + mHH[123] - 165./
   16.*mHH[78] + mHH[122] + mHH[121] + 1./2.*mHH[103] + mHH[120];
   mHH[103]=mHH[97]*mHH[103];
   mHH[104]=9*mHH[64] - 125./8. + mHH[94];
   mHH[104]=mHH[151] + mHH[406] + mHH[113] + mHH[110] - 31./8.*mHH[40]
    + mHH[186] + 3./4.*mHH[81] - 11./8.*mHH[76] + mHH[131] - 55./4.*
   mHH[78] - mHH[69] + 55./8.*mHH[70] - 3./4.*mHH[20] + 1./2.*mHH[104]
    + 6*mHH[79];
   mHH[104]=mHH[96]*mHH[104];
   mHH[105]= - 1 - 11./4.*mHH[33];
   mHH[105]=mHH[134] + mHH[114] + 3*mHH[105] + mHH[133];
   mHH[105]=mHH[97]*mHH[105];
   mHH[105]=mHH[135] + mHH[105];
   mHH[105]=mHH[7]*mHH[105];
   mHH[106]= - 1 - 5./4.*mHH[32];
   mHH[106]=mHH[31] + 11*mHH[106] - mHH[15];
   mHH[106]=mHH[96]*mHH[106];
   mHH[106]=mHH[106] + mHH[341];
   mHH[107]=mHH[4]*mHH[96]*mHH[15];
   mHH[106]=6*mHH[111] + 3./2.*mHH[106] + 2*mHH[107];
   mHH[106]=mHH[8]*mHH[106];
   mHH[107]=mHH[236] + 2*mHH[40] - 2*mHH[94] + mHH[223];
   mHH[107]=mHH[4]*mHH[96]*mHH[107];
   mHH[103]=1./2.*mHH[105] + mHH[144] + mHH[106] + mHH[107] + mHH[103]
    + 3./2.*mHH[448] + 3*mHH[104] + mHH[137] + mHH[140] - 67./4.*
   mHH[56] - 89./4.*mHH[53] + 2*mHH[55] + mHH[298] + 10*mHH[57] + 8*
   mHH[61] - 9*mHH[50];
   mHH[103]=mHH[98]*mHH[103];
   mHH[104]=19 - 11*mHH[67];
   mHH[104]=1./4.*mHH[104] - 3*mHH[63];
   mHH[104]=3*mHH[104] - 1./2.*mHH[18];
   mHH[104]=mHH[360] + 1./8.*mHH[34] + mHH[155] - 1./8.*mHH[10] - 33./
   16.*mHH[33] + 25./16.*mHH[16] - 3./16.*mHH[77] + 1./8.*mHH[68] + 99./
   32.*mHH[76] + 9./4.*mHH[39] + 33./16.*mHH[78] + 3./16.*mHH[19] + 
   mHH[159] + 1./4.*mHH[104] + mHH[197];
   mHH[104]=mHH[97]*mHH[104];
   mHH[105]= - 9*mHH[64] + 57./4. - mHH[94];
   mHH[105]=mHH[332] + mHH[353] + mHH[255] - 33./8.*mHH[32] + 25./8.*
   mHH[40] + 1./2.*mHH[71] - 3./4.*mHH[81] + 33./16.*mHH[76] + 9./2.*
   mHH[39] + 33./4.*mHH[78] + mHH[69] - 33./8.*mHH[70] + 3./4.*mHH[20]
    + 1./2.*mHH[105] - 6*mHH[79];
   mHH[105]=mHH[96]*mHH[105];
   mHH[106]=mHH[205] + mHH[112] + mHH[109] - 3*mHH[20];
   mHH[106]=mHH[96]*mHH[106];
   mHH[107]=mHH[10] - mHH[16] + mHH[18] - 3./2.*mHH[19];
   mHH[107]=mHH[97]*mHH[107];
   mHH[106]=mHH[106] + mHH[107];
   mHH[106]=mHH[4]*mHH[106];
   mHH[107]= - 3*mHH[41] + 5./2. + 3*mHH[73];
   mHH[108]=3./2.*mHH[77];
   mHH[107]= - 7./4.*mHH[82] + mHH[16] + mHH[108] + 1./4.*mHH[107] - 
   mHH[68];
   mHH[107]=mHH[95]*mHH[107];
   mHH[109]= - mHH[31] + mHH[15] - 1 + 33./4.*mHH[32];
   mHH[109]=mHH[96]*mHH[109];
   mHH[109]=mHH[109] + mHH[448];
   mHH[110]= - mHH[4]*mHH[96]*mHH[15];
   mHH[109]=3./2.*mHH[109] + 2*mHH[110];
   mHH[109]=mHH[8]*mHH[109];
   mHH[110]=mHH[383] - mHH[30] + mHH[208] - 1 + 33./4.*mHH[33];
   mHH[110]=mHH[97]*mHH[110];
   mHH[110]=1./4.*mHH[95] + mHH[110];
   mHH[110]=3./4.*mHH[110] + mHH[296];
   mHH[110]=mHH[7]*mHH[110];
   mHH[111]=mHH[50] + mHH[394];
   mHH[104]=mHH[110] + 9./64.*mHH[154] + mHH[109] + mHH[106] + 3*
   mHH[104] + 3./2.*mHH[341] + 3*mHH[105] + 3./16.*mHH[107] - mHH[52]
    + 39./4.*mHH[56] + 39./2.*mHH[53] + 9*mHH[111] - 2*mHH[55];
   mHH[104]=mHH[92]*mHH[104];
   mHH[105]=3 + mHH[32];
   mHH[105]=mHH[96]*mHH[105];
   mHH[105]=mHH[105] + mHH[141];
   mHH[105]=mHH[8]*mHH[105];
   mHH[106]=5*mHH[56] - 4*mHH[61] - 5*mHH[57];
   mHH[107]= - 2*mHH[32] + 4*mHH[78] + 1 - 2*mHH[70];
   mHH[107]=mHH[96]*mHH[107];
   mHH[109]=1 + mHH[78];
   mHH[109]=mHH[97]*mHH[109];
   mHH[103]=mHH[104] + mHH[103] + 6*mHH[105] + 3*mHH[109] + 2*mHH[106]
    + 3*mHH[107];
   mHH[103]=MMZ*mHH[103];
   mHH[104]=33./4.*mHH[139] + mHH[15];
   mHH[105]=4./3.*mHH[164] - mHH[38] + mHH[134] + mHH[145] + mHH[133]
    + mHH[104] + mHH[349];
   mHH[105]=mHH[91]*mHH[105];
   mHH[106]=1 + mHH[33];
   mHH[105]=4*mHH[106] + mHH[105];
   mHH[105]=mHH[8]*mHH[105];
   mHH[106]=mHH[243] + mHH[130] + 23./8.*mHH[82] - 3./4.*mHH[16] + 
   mHH[333] + 9./8.*mHH[77] + 3./4.*mHH[68] + mHH[375] + 15./8.*mHH[41]
    + 1./2.*mHH[36] + mHH[37] - 7./8.*mHH[73] + mHH[329] + 9./2.*
   mHH[81] + 43./8. - mHH[74];
   mHH[107]=mHH[203] + 3./8. + mHH[202];
   mHH[107]=1./2.*mHH[107] + 2*mHH[8];
   mHH[107]=mHH[9]*mHH[107];
   mHH[109]=7 + 23./2.*mHH[34];
   mHH[109]=1./4.*mHH[109] + mHH[327];
   mHH[109]=mHH[7]*mHH[109];
   mHH[110]=mHH[441] + 8 + 13./2.*mHH[34];
   mHH[110]=mHH[8]*mHH[110];
   mHH[106]=mHH[109] + 3*mHH[107] + 3./2.*mHH[106] + mHH[110];
   mHH[106]=mHH[92]*mHH[106];
   mHH[107]= - 1 - mHH[34];
   mHH[107]=mHH[8]*mHH[107];
   mHH[107]=mHH[459] + 5./2.*mHH[457] + 1./2.*mHH[456] + 56./3.*
   mHH[107];
   mHH[107]=mHH[98]*mHH[107];
   mHH[109]=mHH[269] + mHH[363] + mHH[13] + mHH[216] + mHH[217] + 1./2.
   *mHH[29] - mHH[80];
   mHH[109]=mHH[92]*mHH[109];
   mHH[110]= - 9./2. - mHH[81];
   mHH[110]=mHH[370] + mHH[221] + 1./2.*mHH[110] - mHH[83];
   mHH[110]=mHH[92]*mHH[110];
   mHH[111]=mHH[83] + 4 + 1./2.*mHH[81];
   mHH[111]=mHH[82] + 3*mHH[111] + 13./2.*mHH[77];
   mHH[111]=mHH[98]*mHH[111];
   mHH[110]=mHH[111] + 3*mHH[110];
   mHH[110]=MMZ*mHH[110];
   mHH[111]=mHH[98]*mHH[231];
   mHH[109]=mHH[110] + 13*mHH[111] + 3*mHH[109];
   mHH[109]=mHH[1]*mHH[109];
   mHH[110]=32./3.*mHH[161] + 3./2.*mHH[160];
   mHH[110]=mHH[8]*mHH[110];
   mHH[111]=mHH[452] - mHH[54];
   mHH[111]=mHH[98]*mHH[111];
   mHH[112]=2*mHH[58] + mHH[54];
   mHH[112]=mHH[92]*mHH[112];
   mHH[111]=2*mHH[111] + 3*mHH[112];
   mHH[111]=MMZ*mHH[111];
   mHH[106]=mHH[109] + mHH[111] + mHH[106] + mHH[107] + mHH[110] + 3./4.
   *mHH[189];
   mHH[106]=mHH[1]*mHH[106];
   mHH[107]=mHH[40] - mHH[83] - mHH[71] - 1 + 3./2.*mHH[81];
   mHH[107]=mHH[96]*mHH[107];
   mHH[108]= - mHH[34] - mHH[82] + mHH[16] + mHH[108] - 1 - mHH[68];
   mHH[108]=mHH[97]*mHH[108];
   mHH[107]=1./8.*mHH[472] + mHH[171] + 1./8.*mHH[108] + mHH[167] + 1./
   4.*mHH[107] - mHH[58] - 1./2.*mHH[54];
   mHH[107]=mHH[92]*mHH[107];
   mHH[108]=mHH[370] + mHH[193] - 1./2.*mHH[83] + 1./2. + mHH[81];
   mHH[108]=mHH[92]*mHH[108];
   mHH[109]=mHH[98]*mHH[434];
   mHH[108]=mHH[109] + 3*mHH[108];
   mHH[108]=mHH[1]*mHH[108];
   mHH[109]=mHH[98]*mHH[446];
   mHH[110]= - mHH[54] - mHH[58] + 1./2.*mHH[60];
   mHH[110]=mHH[92]*mHH[110];
   mHH[108]=mHH[108] + mHH[109] + 3*mHH[110];
   mHH[108]=MMt*mHH[1]*mHH[108];
   mHH[109]=mHH[98]*mHH[213];
   mHH[106]=mHH[108] + mHH[106] + mHH[109] + 3*mHH[107];
   mHH[106]=MMt*mHH[106];
   mHH[104]=2./3.*mHH[164] + mHH[243] + mHH[148] + mHH[30] + mHH[147]
    + 1./2.*mHH[104] - mHH[31];
   mHH[104]=mHH[7]*mHH[91]*mHH[104];
   mHH[107]=5./32.*mHH[157] + 5./32.*mHH[156] + mHH[142] + 17./8.*
   mHH[56] - 3./4.*mHH[53] - 9./2.*mHH[49] + mHH[300] - 3*mHH[57];
   mHH[107]=MMH*mHH[98]*mHH[107];
   mHH[108]=3 + 2*mHH[78];

      return mHH[99] + mHH[100] + mHH[101] + mHH[102] + mHH[103] + 
      mHH[104] + mHH[105] + mHH[106] + mHH[107] + 3*mHH[108] + mHH[116];
}
