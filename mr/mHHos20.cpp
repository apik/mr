#include <HH.hpp>
std::complex<long double>
HH<MS>::m20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armHHos[208], mHHosret;

    armHHos[1]=double(nL + nH);
    armHHos[2]=pow(s,-1);
    armHHos[3]=pow(c,-1);
    armHHos[4]=pow(mmH,-1);
    armHHos[5]=Tsil::I2(0,0,mmW,mu2);
    armHHos[6]=Tsil::I2(0,0,mmZ,mu2);
    armHHos[7]=Tsil::B(mmW,mmW,mmH,mu2);
    armHHos[8]=pow(mmZ,-1);
    armHHos[9]=Tsil::B(mmZ,mmZ,mmH,mu2);
    armHHos[10]=Tsil::B(0,mmW,mmH,mu2);
    armHHos[11]=Tsil::B(0,mmZ,mmH,mu2);
    armHHos[12]=Tsil::Beps(mmW,mmW,mmH,mu2);
    armHHos[13]=Tsil::Beps(mmZ,mmZ,mmH,mu2);
    armHHos[14]=Tsil::A(mmW,mu2);
    armHHos[15]=Tsil::A(mmZ,mu2);
    armHHos[16]=Tsil::Aeps(mmW,mu2);
    armHHos[17]=Tsil::Aeps(mmZ,mu2);
    armHHos[18]=protZZ00->Uxzuv(0);
    armHHos[19]=protWW00->Uxzuv(0);
    armHHos[20]=protZZ00->Txuv(0);
    armHHos[21]=protWW00->Txuv(0);
    armHHos[22]=double(nH);
    armHHos[23]=Tsil::B(mmt,mmt,mmH,mu2);
    armHHos[24]=Tsil::A(mmt,mu2);
    armHHos[25]=Tsil::I2(mmZ,mmt,mmt,mu2);
    armHHos[26]=Tsil::I2(mmH,mmt,mmt,mu2);
    armHHos[27]=Tsil::I2(0,mmW,mmt,mu2);
    armHHos[28]=Tsil::B(mmH,mmH,mmH,mu2);
    armHHos[29]=Tsil::A(mmH,mu2);
    armHHos[30]=Tsil::Beps(mmH,mmH,mmH,mu2);
    armHHos[31]=Tsil::Beps(mmt,mmt,mmH,mu2);
    armHHos[32]=Tsil::Aeps(mmH,mu2);
    armHHos[33]=Tsil::Aeps(mmt,mu2);
    armHHos[34]=prottttt0->M(0);
    armHHos[35]=prottttt0->Vxzuv(0);
    armHHos[36]=prottttt0->Suxv(0);
    armHHos[37]=protHtHtt->M(0);
    armHHos[38]=protZtZtt->M(0);
    armHHos[39]=protWtWt0->M(0);
    armHHos[40]=protttttH->M(0);
    armHHos[41]=protttttZ->M(0);
    armHHos[42]=protHtHtt->Uzxyv(0);
    armHHos[43]=protZtZtt->Uzxyv(0);
    armHHos[44]=protWtWt0->Uzxyv(0);
    armHHos[45]=protHtHtt->Uuyxv(0);
    armHHos[46]=protZtZtt->Uuyxv(0);
    armHHos[47]=protWtWt0->Uuyxv(0);
    armHHos[48]=protZtZtt->Txuv(0);
    armHHos[49]=protWtWt0->Suxv(0);
    armHHos[50]=protWtWt0->Txuv(0);
    armHHos[51]=protZtZtt->Tyzv(0);
    armHHos[52]=protWtWt0->Tyzv(0);
    armHHos[53]=protHtHtt->Tyzv(0);
    armHHos[54]=protHtHtt->Svyz(0);
    armHHos[55]=protZtZtt->Svyz(0);
    armHHos[56]=double(nL);
    armHHos[57]=double(boson);
    armHHos[58]=Tsil::I2(mmW,mmW,mmZ,mu2);
    armHHos[59]=Tsil::I2(mmH,mmW,mmW,mu2);
    armHHos[60]=Tsil::I2(mmH,mmZ,mmZ,mu2);
    armHHos[61]=Tsil::I2(mmH,mmH,mmH,mu2);
    armHHos[62]=protHHHHH->M(0);
    armHHos[63]=protHZHZZ->M(0);
    armHHos[64]=protHWHWW->M(0);
    armHHos[65]=protZZZZH->M(0);
    armHHos[66]=protZWZWW->M(0);
    armHHos[67]=protWWWWH->M(0);
    armHHos[68]=protWWWWZ->M(0);
    armHHos[69]=protWWWW0->M(0);
    armHHos[70]=protWWWW0->Vxzuv(0);
    armHHos[71]=protHHHHH->Uzxyv(0);
    armHHos[72]=protHZHZZ->Uzxyv(0);
    armHHos[73]=protHWHWW->Uzxyv(0);
    armHHos[74]=protHZHZZ->Uuyxv(0);
    armHHos[75]=protZWZWW->Uzxyv(0);
    armHHos[76]=protHWHWW->Uuyxv(0);
    armHHos[77]=protZWZWW->Uuyxv(0);
    armHHos[78]=protHZHZZ->Tyzv(0);
    armHHos[79]=protZWZWW->Txuv(0);
    armHHos[80]=protZWZWW->Tuxv(0);
    armHHos[81]=protHWHWW->Tuxv(0);
    armHHos[82]=protHZHZZ->Svyz(0);
    armHHos[83]=protHWHWW->Suxv(0);
    armHHos[84]=protZWZWW->Suxv(0);
    armHHos[85]=protWWWW0->Suxv(0);
    armHHos[86]=1/(4*mmt - mmZ);
    armHHos[87]=1/(mmH - 4*mmZ);
    armHHos[88]=1/(mmH - 4*mmW);
   armHHos[89]=3*armHHos[20];
   armHHos[90]=5 + armHHos[89];
   armHHos[91]= - 3*armHHos[18];
   armHHos[92]=9./2.*armHHos[11];
   armHHos[90]=9./2.*armHHos[48] - 3*armHHos[43] + armHHos[92] + 1./2.*
   armHHos[90] + armHHos[91];
   armHHos[93]=3*armHHos[13];
   armHHos[90]=armHHos[9] + 1./2.*armHHos[90] + armHHos[93];
   armHHos[90]=armHHos[87]*armHHos[90];
   armHHos[94]=3*armHHos[31] - 5./2. - 3*armHHos[46];
   armHHos[95]= - 3./2.*armHHos[48];
   armHHos[94]= - armHHos[13] + 7./4.*armHHos[51] + armHHos[95] + 1./4.
   *armHHos[94] + armHHos[43];
   armHHos[94]=armHHos[86]*armHHos[94];
   armHHos[96]=3*armHHos[12];
   armHHos[97]= - 3*armHHos[44];
   armHHos[98]=armHHos[7] + armHHos[96] + armHHos[97] - 1 + 9./2.*
   armHHos[50];
   armHHos[98]=armHHos[88]*armHHos[98];
   armHHos[99]= - armHHos[9]*armHHos[86];
   armHHos[100]=armHHos[23]*armHHos[86];
   armHHos[90]=1./2.*armHHos[90] + 3./8.*armHHos[99] + 9./32.*
   armHHos[100] + 3./8.*armHHos[94] + armHHos[98];
   armHHos[90]=mmZ*armHHos[90];
   armHHos[94]= - 5./2.*armHHos[17];
   armHHos[98]= - 1./2.*armHHos[25];
   armHHos[101]=armHHos[94] + armHHos[98] + 15*armHHos[33];
   armHHos[101]=armHHos[86]*armHHos[101];
   armHHos[102]= - 1./2.*armHHos[6];
   armHHos[103]=armHHos[98] + armHHos[102] + armHHos[55];
   armHHos[103]=1./2.*armHHos[103] + armHHos[17];
   armHHos[104]= - armHHos[9]*armHHos[24];
   armHHos[103]=3*armHHos[104] + 3*armHHos[103] - 11./4.*armHHos[15];
   armHHos[103]=armHHos[87]*armHHos[103];
   armHHos[105]=armHHos[13] - armHHos[51] + 1 - armHHos[43];
   armHHos[105]=armHHos[86]*armHHos[105];
   armHHos[106]=armHHos[9]*armHHos[86];
   armHHos[105]=armHHos[105] + armHHos[106];
   armHHos[105]=mmH*armHHos[105];
   armHHos[107]= - 3./2.*armHHos[18];
   armHHos[108]=13*armHHos[50] + armHHos[107] + 659./32. + 9*
   armHHos[47];
   armHHos[109]= - 1./2.*armHHos[27] + armHHos[49];
   armHHos[110]=armHHos[109] + armHHos[16];
   armHHos[110]=armHHos[88]*armHHos[110];
   armHHos[111]= - armHHos[88]*armHHos[7];
   armHHos[111]= - 13./16.*armHHos[86] + armHHos[111];
   armHHos[111]=armHHos[24]*armHHos[111];
   armHHos[112]=armHHos[24]*armHHos[86];
   armHHos[113]=1./4. + armHHos[112];
   armHHos[113]=armHHos[23]*armHHos[113];
   armHHos[114]= - armHHos[23]*armHHos[86];
   armHHos[115]=7./2.*armHHos[86] + 3*armHHos[114];
   armHHos[115]=armHHos[15]*armHHos[115];
   armHHos[116]=17./4. + 3*armHHos[112];
   armHHos[116]=armHHos[9]*armHHos[116];
   armHHos[117]=5*armHHos[7];
   armHHos[118]= - 3./8.*armHHos[43];
   armHHos[119]= - armHHos[14]*armHHos[88];
   armHHos[90]=armHHos[90] + 3./32.*armHHos[105] + 1./2.*armHHos[103]
    + 1./2.*armHHos[116] + 3./16.*armHHos[115] + 9./8.*armHHos[113] + 3
   *armHHos[111] + 5*armHHos[119] + 3*armHHos[110] + armHHos[117] + 3./
   16.*armHHos[101] + 9./8.*armHHos[13] + 21./32.*armHHos[51] + 
   armHHos[96] + 35./16.*armHHos[48] + armHHos[118] - 7./8.*armHHos[11]
    - 207./32.*armHHos[31] + 63./32.*armHHos[46] + 1./2.*armHHos[108]
    + armHHos[97];
   armHHos[101]= - armHHos[31] + 3 + armHHos[46];
   armHHos[101]=armHHos[13] - 5./4.*armHHos[51] + 1./2.*armHHos[48] + 1.
   /4.*armHHos[101] - armHHos[43];
   armHHos[101]=armHHos[86]*armHHos[101];
   armHHos[101]=armHHos[106] + armHHos[101] + 1./4.*armHHos[114];
   armHHos[101]=mmZ*armHHos[101];
   armHHos[103]= - 641./4. + armHHos[89];
   armHHos[108]=15*armHHos[18];
   armHHos[103]=1./2.*armHHos[103] + armHHos[108];
   armHHos[110]=5./2.*armHHos[11];
   armHHos[111]=1./2.*armHHos[17];
   armHHos[115]= - 5*armHHos[33] + armHHos[111];
   armHHos[115]=armHHos[86]*armHHos[115];
   armHHos[116]= - armHHos[24]*armHHos[86];
   armHHos[120]= - 1./4. + armHHos[116];
   armHHos[120]=armHHos[23]*armHHos[120];
   armHHos[121]= - armHHos[86] + armHHos[100];
   armHHos[121]=armHHos[15]*armHHos[121];
   armHHos[122]= - 59./4. + 9*armHHos[116];
   armHHos[122]=armHHos[9]*armHHos[122];
   armHHos[103]=9./8.*armHHos[101] + 1./2.*armHHos[122] + 9./16.*
   armHHos[121] + 9./8.*armHHos[120] + 45./8.*armHHos[112] - 17*
   armHHos[7] + 9./8.*armHHos[115] - 51./8.*armHHos[13] - 45./32.*
   armHHos[51] - 15*armHHos[12] - 53./16.*armHHos[48] + 21./8.*
   armHHos[43] + armHHos[110] - 9./32.*armHHos[31] + 9./32.*armHHos[46]
    + 15*armHHos[44] + 1./4.*armHHos[103] - 11*armHHos[50];
   armHHos[103]=mmZ*armHHos[103];
   armHHos[122]= - 1./2.*armHHos[48];
   armHHos[123]= - 3./2. + 2*armHHos[7];
   armHHos[123]=armHHos[23]*armHHos[123];
   armHHos[124]=1./2. + armHHos[23];
   armHHos[124]=armHHos[9]*armHHos[124];
   armHHos[125]=3*armHHos[7];
   armHHos[126]=1./2.*armHHos[51];
   armHHos[123]=3*armHHos[124] + 3*armHHos[123] + armHHos[125] + 
   armHHos[126] + armHHos[122] + 1./2.*armHHos[52] - 7./2. - 
   armHHos[50];
   armHHos[123]=mmt*armHHos[123];
   armHHos[127]=9./2. + armHHos[50];
   armHHos[128]=1./4.*armHHos[48];
   armHHos[127]=armHHos[126] + armHHos[128] + 1./2.*armHHos[127] + 
   armHHos[52];
   armHHos[127]=mmZ*armHHos[127];
   armHHos[129]= - 1./4.*armHHos[25];
   armHHos[130]= - 3*armHHos[7];
   armHHos[131]= - 11 + armHHos[130];
   armHHos[131]=armHHos[24]*armHHos[131];
   armHHos[132]=1./2.*armHHos[55];
   armHHos[109]=armHHos[123] + armHHos[127] + 3./2.*armHHos[104] - 1./2.
   *armHHos[15] + armHHos[131] - armHHos[14] + armHHos[129] + 
   armHHos[109] + armHHos[132];
   armHHos[109]=mmt*armHHos[109];
   armHHos[123]=armHHos[128] - 1./4.*armHHos[11] + armHHos[50] + 1 - 1./
   4.*armHHos[20];
   armHHos[123]=mmZ*armHHos[123];
   armHHos[127]= - 2*armHHos[14];
   armHHos[123]=armHHos[123] - 1./4.*armHHos[15] - 3*armHHos[24] + 
   armHHos[127] + armHHos[129] + armHHos[132] + 2*armHHos[49] - 1./4.*
   armHHos[6] - armHHos[27];
   armHHos[123]=mmZ*armHHos[123];
   armHHos[109]=armHHos[123] + armHHos[109];
   armHHos[109]=armHHos[4]*armHHos[109];
   armHHos[123]= - 9./2.*armHHos[50];
   armHHos[128]= - 21./4. - 5*armHHos[7];
   armHHos[128]=armHHos[23]*armHHos[128];
   armHHos[129]=3 - 5./2.*armHHos[23];
   armHHos[129]=armHHos[9]*armHHos[129];
   armHHos[131]=3./2.*armHHos[44];
   armHHos[133]= - 3./2.*armHHos[12];
   armHHos[134]= - 1./2.*armHHos[52];
   armHHos[128]=1./2.*armHHos[129] + 1./2.*armHHos[128] + 3./4.*
   armHHos[13] - 23./8.*armHHos[51] + armHHos[133] - 9./8.*armHHos[48]
    - 3./4.*armHHos[43] + armHHos[134] - 15./8.*armHHos[31] + 7./8.*
   armHHos[46] + armHHos[131] + armHHos[123] - 31./8. + armHHos[47];
   armHHos[129]= - 2*armHHos[39] - armHHos[38];
   armHHos[129]=mmZ*armHHos[129];
   armHHos[135]=armHHos[38] + armHHos[39] - 1./2.*armHHos[41];
   armHHos[135]=mmt*armHHos[135];
   armHHos[128]=armHHos[135] + 1./2.*armHHos[128] + armHHos[129];
   armHHos[128]=mmt*armHHos[128];
   armHHos[129]=11./8.*armHHos[6] + 7*armHHos[27];
   armHHos[129]= - 111./16.*armHHos[17] - 15*armHHos[16] - 111./8.*
   armHHos[33] + 51./8.*armHHos[25] - 13./2.*armHHos[55] + 3*
   armHHos[129] - 17*armHHos[49];
   armHHos[135]=4*armHHos[14];
   armHHos[136]=53./8. + armHHos[117];
   armHHos[136]=armHHos[24]*armHHos[136];
   armHHos[137]=9*armHHos[23];
   armHHos[138]=7 + armHHos[137];
   armHHos[138]=armHHos[15]*armHHos[138];
   armHHos[139]=armHHos[9]*armHHos[24];
   armHHos[140]=3./2.*armHHos[139];
   armHHos[141]= - armHHos[23]*armHHos[24];
   armHHos[103]=3*armHHos[109] + 3*armHHos[128] + 1./2.*armHHos[103] + 
   armHHos[140] + 1./32.*armHHos[138] + 9./16.*armHHos[141] + 3./2.*
   armHHos[136] + 1./2.*armHHos[129] + armHHos[135];
   armHHos[103]=armHHos[4]*armHHos[103];
   armHHos[109]= - 3*armHHos[23];
   armHHos[128]=armHHos[109] - armHHos[13] + armHHos[51] + armHHos[95]
    + 3 + armHHos[43];
   armHHos[129]=3./2.*armHHos[23];
   armHHos[136]= - 1 + armHHos[129];
   armHHos[136]=armHHos[9]*armHHos[136];
   armHHos[128]=1./2.*armHHos[128] + armHHos[136];
   armHHos[128]=armHHos[87]*armHHos[128];
   armHHos[136]= - 1./2.*armHHos[50];
   armHHos[138]=1 + armHHos[136];
   armHHos[138]= - armHHos[12] + armHHos[52] + 3*armHHos[138] + 
   armHHos[44];
   armHHos[138]=1./2.*armHHos[138] - armHHos[7];
   armHHos[138]=armHHos[88]*armHHos[138];
   armHHos[142]= - 1 + armHHos[7];
   armHHos[143]=armHHos[88]*armHHos[142];
   armHHos[144]=armHHos[23]*armHHos[143];
   armHHos[128]=1./4.*armHHos[128] + 3./4.*armHHos[144] + 1./2.*
   armHHos[138] + armHHos[39] + 1./2.*armHHos[38];
   armHHos[128]=mmt*armHHos[128];
   armHHos[90]=armHHos[103] + 1./2.*armHHos[90] + 3*armHHos[128];
   armHHos[103]=pow(armHHos[2],2);
   armHHos[90]=armHHos[103]*armHHos[90];
   armHHos[128]=armHHos[31] - 3 - armHHos[46];
   armHHos[128]= - armHHos[13] + 5./4.*armHHos[51] + armHHos[122] + 1./
   4.*armHHos[128] + armHHos[43];
   armHHos[128]=armHHos[86]*armHHos[128];
   armHHos[138]=1./4.*armHHos[100];
   armHHos[128]=armHHos[99] + armHHos[128] + armHHos[138];
   armHHos[128]=15./4.*mmZ*armHHos[128];
   armHHos[144]=853./12. + armHHos[20];
   armHHos[145]=5*armHHos[18];
   armHHos[144]=1./2.*armHHos[144] + armHHos[145];
   armHHos[146]=5./3.*armHHos[11];
   armHHos[147]=5./4.*armHHos[43];
   armHHos[148]=17./24.*armHHos[48];
   armHHos[149]=75./16.*armHHos[51];
   armHHos[150]= - 15./4.*armHHos[13];
   armHHos[151]=5*armHHos[33];
   armHHos[152]=armHHos[151] - 1./2.*armHHos[17];
   armHHos[152]=15./4.*armHHos[86]*armHHos[152];
   armHHos[153]=17*armHHos[7];
   armHHos[154]=75./4.*armHHos[116];
   armHHos[155]=15./4.*armHHos[113];
   armHHos[156]=armHHos[86] + armHHos[114];
   armHHos[156]=15./8.*armHHos[15]*armHHos[156];
   armHHos[157]= - 1./4. + armHHos[112];
   armHHos[157]=armHHos[9]*armHHos[157];
   armHHos[158]=15*armHHos[157];
   armHHos[144]=armHHos[128] + armHHos[158] + armHHos[156] + 
   armHHos[155] + armHHos[154] + armHHos[153] + armHHos[152] + 
   armHHos[150] + armHHos[149] + 15*armHHos[12] + armHHos[148] + 
   armHHos[147] + armHHos[146] + 15./16.*armHHos[31] - 15./16.*
   armHHos[46] - 15*armHHos[44] + 1./2.*armHHos[144] + 11*armHHos[50];
   armHHos[144]=mmZ*armHHos[144];
   armHHos[159]= - 1./2.*armHHos[20];
   armHHos[160]= - 1./2.*armHHos[11];
   armHHos[161]=armHHos[122] + armHHos[160] - 6*armHHos[50] - 7 + 
   armHHos[159];
   armHHos[161]=mmZ*armHHos[161];
   armHHos[162]=3./2.*armHHos[15];
   armHHos[163]=8*armHHos[24];
   armHHos[164]=6*armHHos[14];
   armHHos[165]=1./2.*armHHos[25];
   armHHos[161]=armHHos[161] + armHHos[162] + armHHos[163] + 
   armHHos[164] + armHHos[165] - armHHos[55] - 6*armHHos[49] + 
   armHHos[102] + 3*armHHos[27];
   armHHos[161]=mmZ*armHHos[161];
   armHHos[166]= - armHHos[55] + armHHos[165];
   armHHos[163]=9*armHHos[104] + 13*armHHos[15] + 13*armHHos[166] + 
   armHHos[163];
   armHHos[167]=1 + 2*armHHos[23];
   armHHos[167]=armHHos[9]*armHHos[167];
   armHHos[168]= - 9*armHHos[23];
   armHHos[167]=9*armHHos[167] + armHHos[168] - 13*armHHos[51] - 30 - 
   11*armHHos[48];
   armHHos[167]=mmt*armHHos[167];
   armHHos[136]= - armHHos[52] - 4 + armHHos[136];
   armHHos[136]= - armHHos[51] + 3*armHHos[136] - 13./2.*armHHos[48];
   armHHos[136]=mmZ*armHHos[136];
   armHHos[136]=armHHos[167] + armHHos[163] + armHHos[136];
   armHHos[136]=mmt*armHHos[136];
   armHHos[136]=armHHos[161] + armHHos[136];
   armHHos[136]=armHHos[4]*armHHos[136];
   armHHos[161]=19*armHHos[31] + 295./3. - 19*armHHos[46];
   armHHos[169]= - 3./2.*armHHos[23];
   armHHos[161]=armHHos[169] + 33*armHHos[13] + 83./2.*armHHos[51] + 77.
   /2.*armHHos[48] + 1./2.*armHHos[161] - 33*armHHos[43];
   armHHos[170]= - 3 + 17./2.*armHHos[23];
   armHHos[170]=armHHos[9]*armHHos[170];
   armHHos[161]=1./2.*armHHos[161] + armHHos[170];
   armHHos[161]=1./2.*armHHos[161];
   armHHos[170]=3*armHHos[39];
   armHHos[171]=armHHos[170] + armHHos[38];
   armHHos[171]=mmZ*armHHos[171];
   armHHos[172]=armHHos[41] + 22*armHHos[38];
   armHHos[172]=mmt*armHHos[172];
   armHHos[171]=armHHos[172] + armHHos[161] + 2*armHHos[171];
   armHHos[171]=mmt*armHHos[171];
   armHHos[173]= - 5 + armHHos[109];
   armHHos[173]=armHHos[15]*armHHos[173];
   armHHos[174]=armHHos[23]*armHHos[24];
   armHHos[173]=5./8.*armHHos[173] + 15./4.*armHHos[174] - 287./12.*
   armHHos[24] - 15./8.*armHHos[17] + 247./12.*armHHos[33] - 17./4.*
   armHHos[25] + 11./4.*armHHos[6] + 13./3.*armHHos[55];
   armHHos[173]=1./2.*armHHos[173] + 5*armHHos[139];
   armHHos[136]=armHHos[136] + armHHos[171] + armHHos[173] + 1./2.*
   armHHos[144];
   armHHos[136]=armHHos[4]*armHHos[136];
   armHHos[144]= - 7./6.*armHHos[11];
   armHHos[171]= - 3./2.*armHHos[43];
   armHHos[175]=5./2.*armHHos[17] + armHHos[165] - 15*armHHos[33];
   armHHos[175]=armHHos[86]*armHHos[175];
   armHHos[100]= - 7./2.*armHHos[86] + 3*armHHos[100];
   armHHos[100]=armHHos[15]*armHHos[100];
   armHHos[176]=5./2.*armHHos[13];
   armHHos[100]=5./4.*armHHos[100] + 15./2.*armHHos[120] + 65./4.*
   armHHos[112] + 5./4.*armHHos[175] + armHHos[176] - 35./8.*
   armHHos[51] + 1./12.*armHHos[48] + armHHos[171] + armHHos[144] + 9./
   8.*armHHos[31] - 9./8.*armHHos[46] - 43./48. - armHHos[18];
   armHHos[175]=armHHos[165] + armHHos[102] - armHHos[55];
   armHHos[177]=armHHos[175] + armHHos[162];
   armHHos[177]=1./2.*armHHos[177] + armHHos[139];
   armHHos[177]=armHHos[87]*armHHos[177];
   armHHos[178]= - armHHos[13] + armHHos[51] - 1 + armHHos[43];
   armHHos[178]=armHHos[86]*armHHos[178];
   armHHos[178]=armHHos[178] + armHHos[99];
   armHHos[178]=mmH*armHHos[178];
   armHHos[179]=1./4. + armHHos[116];
   armHHos[179]=armHHos[9]*armHHos[179];
   armHHos[100]=5./16.*armHHos[178] + armHHos[177] + 1./2.*armHHos[100]
    + 5*armHHos[179];
   armHHos[177]=3 + armHHos[20];
   armHHos[178]=3./2.*armHHos[11];
   armHHos[95]=armHHos[95] + armHHos[43] + armHHos[178] + 1./2.*
   armHHos[177] - armHHos[18];
   armHHos[95]=armHHos[87]*armHHos[95];
   armHHos[177]=3*armHHos[46];
   armHHos[179]= - 3*armHHos[31] + 5./2. + armHHos[177];
   armHHos[180]=3./2.*armHHos[48];
   armHHos[179]=armHHos[13] - 7./4.*armHHos[51] + armHHos[180] + 1./4.*
   armHHos[179] - armHHos[43];
   armHHos[179]=armHHos[86]*armHHos[179];
   armHHos[181]= - 3*armHHos[12];
   armHHos[123]= - armHHos[7] + armHHos[181] + 3*armHHos[44] + 1 + 
   armHHos[123];
   armHHos[123]=armHHos[88]*armHHos[123];
   armHHos[123]=1./2.*armHHos[95] + 5./4.*armHHos[106] + 15./16.*
   armHHos[114] + 5./4.*armHHos[179] + armHHos[123];
   armHHos[123]=mmZ*armHHos[123];
   armHHos[123]=armHHos[100] + armHHos[123];
   armHHos[182]=armHHos[168] - 11*armHHos[13] - armHHos[51] - 33./2.*
   armHHos[48] - 3 + 11*armHHos[43];
   armHHos[183]=9./2.*armHHos[23];
   armHHos[184]=1 + armHHos[183];
   armHHos[184]=armHHos[9]*armHHos[184];
   armHHos[182]=1./2.*armHHos[182] + armHHos[184];
   armHHos[182]=armHHos[87]*armHHos[182];
   armHHos[182]= - 5*armHHos[38] + 1./2.*armHHos[182];
   armHHos[182]=mmt*armHHos[182];
   armHHos[90]=armHHos[90] + armHHos[136] + 1./2.*armHHos[123] + 
   armHHos[182];
   armHHos[90]=armHHos[103]*armHHos[90];
   armHHos[123]=61./12. + armHHos[20];
   armHHos[123]=15./8.*armHHos[31] - 15./8.*armHHos[46] + 1./2.*
   armHHos[123] + armHHos[145];
   armHHos[123]=armHHos[128] + armHHos[158] + armHHos[156] + 
   armHHos[155] + armHHos[154] + armHHos[152] + armHHos[150] + 
   armHHos[149] + armHHos[148] + armHHos[147] + 1./2.*armHHos[123] + 
   armHHos[146];
   armHHos[123]=mmZ*armHHos[123];
   armHHos[128]=armHHos[122] + armHHos[160] - 1 + armHHos[159];
   armHHos[128]=mmZ*armHHos[128];
   armHHos[136]=2*armHHos[24];
   armHHos[128]=armHHos[128] + armHHos[162] + armHHos[175] + 
   armHHos[136];
   armHHos[128]=mmZ*armHHos[128];
   armHHos[147]= - 15 - 13*armHHos[48];
   armHHos[147]=1./2.*armHHos[147] - armHHos[51];
   armHHos[147]=mmZ*armHHos[147];
   armHHos[147]=armHHos[167] + armHHos[163] + armHHos[147];
   armHHos[147]=mmt*armHHos[147];
   armHHos[128]=armHHos[128] + armHHos[147];
   armHHos[128]=armHHos[4]*armHHos[128];
   armHHos[147]=mmZ*armHHos[38];
   armHHos[147]=armHHos[172] + armHHos[161] + 2*armHHos[147];
   armHHos[147]=mmt*armHHos[147];
   armHHos[123]=armHHos[128] + armHHos[147] + armHHos[173] + 1./2.*
   armHHos[123];
   armHHos[123]=armHHos[4]*armHHos[123];
   armHHos[106]=armHHos[106] + armHHos[179] + 3./4.*armHHos[114];
   armHHos[95]=5./2.*armHHos[106] + armHHos[95];
   armHHos[95]=mmZ*armHHos[95];
   armHHos[95]=armHHos[100] + 1./2.*armHHos[95];
   armHHos[100]=5*armHHos[20];
   armHHos[106]= - 2041./12. + armHHos[100];
   armHHos[106]=1./2.*armHHos[106] + 25*armHHos[18];
   armHHos[106]= - 25./8.*armHHos[31] + 1./3.*armHHos[106] + 25./8.*
   armHHos[46];
   armHHos[128]=25*armHHos[116];
   armHHos[147]= - 523./36. + armHHos[128];
   armHHos[147]=armHHos[9]*armHHos[147];
   armHHos[101]=25./4.*armHHos[101] + armHHos[147] + 25./8.*
   armHHos[121] + 25./4.*armHHos[120] + 125./4.*armHHos[112] + 25./4.*
   armHHos[115] - 145./12.*armHHos[13] - 125./16.*armHHos[51] - 829./72.
   *armHHos[48] + 95./12.*armHHos[43] + 1./2.*armHHos[106] + 25./9.*
   armHHos[11];
   armHHos[101]=mmZ*armHHos[101];
   armHHos[106]=2179./12.*armHHos[24] - 365./8.*armHHos[17] - 1499./12.
   *armHHos[33] + 289./4.*armHHos[25] + 55./4.*armHHos[6] - 221./3.*
   armHHos[55];
   armHHos[115]=671./9. + 25*armHHos[23];
   armHHos[115]=armHHos[15]*armHHos[115];
   armHHos[106]=1./8.*armHHos[115] + 1./3.*armHHos[106] + 25./4.*
   armHHos[141];
   armHHos[101]=1./2.*armHHos[101] + 1./2.*armHHos[106] + 5./3.*
   armHHos[139];
   armHHos[106]=25 + 41./3.*armHHos[48];
   armHHos[115]=17./3.*armHHos[51];
   armHHos[106]=1./2.*armHHos[106] + armHHos[115];
   armHHos[106]=mmZ*armHHos[106];
   armHHos[120]= - 9./2.*armHHos[23];
   armHHos[121]=9*armHHos[124] + armHHos[120] + 41./6.*armHHos[51] + 5
    + 7./6.*armHHos[48];
   armHHos[121]=mmt*armHHos[121];
   armHHos[98]=armHHos[55] + armHHos[98];
   armHHos[124]= - 41./2.*armHHos[15] + 41./2.*armHHos[98] - 68*
   armHHos[24];
   armHHos[106]=armHHos[121] + 1./2.*armHHos[106] + 1./3.*armHHos[124]
    + 9./2.*armHHos[104];
   armHHos[106]=mmt*armHHos[106];
   armHHos[121]= - 17./2.*armHHos[25] - 5./2.*armHHos[6] + 17*
   armHHos[55];
   armHHos[121]=1./2.*armHHos[121];
   armHHos[124]= - 29./4.*armHHos[15] + armHHos[121] - 17*armHHos[24];
   armHHos[147]=17./12.*armHHos[48] - 5./12.*armHHos[11] + 1 - 5./12.*
   armHHos[20];
   armHHos[147]=mmZ*armHHos[147];
   armHHos[124]=1./3.*armHHos[124] + armHHos[147];
   armHHos[124]=mmZ*armHHos[124];
   armHHos[106]=armHHos[124] + armHHos[106];
   armHHos[106]=armHHos[4]*armHHos[106];
   armHHos[124]= - 143*armHHos[31] - 875./3. + 143*armHHos[46];
   armHHos[124]= - 43./2.*armHHos[23] - 7*armHHos[13] - 511./6.*
   armHHos[51] - 289./6.*armHHos[48] + 1./6.*armHHos[124] + 7*
   armHHos[43];
   armHHos[147]=17 - 109./6.*armHHos[23];
   armHHos[147]=armHHos[9]*armHHos[147];
   armHHos[124]=1./2.*armHHos[124] + armHHos[147];
   armHHos[147]= - mmZ*armHHos[38];
   armHHos[148]= - 17./2.*armHHos[41] - 7*armHHos[38];
   armHHos[148]=mmt*armHHos[148];
   armHHos[124]=1./3.*armHHos[148] + 1./4.*armHHos[124] + 17./3.*
   armHHos[147];
   armHHos[124]=mmt*armHHos[124];
   armHHos[101]=armHHos[106] + 1./2.*armHHos[101] + armHHos[124];
   armHHos[101]=armHHos[4]*armHHos[101];
   armHHos[106]= - 5*armHHos[18];
   armHHos[124]=2203./48. + armHHos[106];
   armHHos[147]= - 5./6.*armHHos[17] - 1./6.*armHHos[25] + armHHos[151]
   ;
   armHHos[147]=armHHos[86]*armHHos[147];
   armHHos[114]=7./6.*armHHos[86] + armHHos[114];
   armHHos[114]=armHHos[15]*armHHos[114];
   armHHos[113]=25./4.*armHHos[114] + 25./2.*armHHos[113] + 325./12.*
   armHHos[116] + 25./4.*armHHos[147] + 19./6.*armHHos[13] + 175./24.*
   armHHos[51] + 523./36.*armHHos[48] + armHHos[171] - 35./18.*
   armHHos[11] - 111./8.*armHHos[31] + 1./3.*armHHos[124] + 111./8.*
   armHHos[46];
   armHHos[114]=armHHos[31] - 5./6. - armHHos[46];
   armHHos[114]= - 1./3.*armHHos[13] + 7./12.*armHHos[51] + 
   armHHos[122] + 1./4.*armHHos[114] + 1./3.*armHHos[43];
   armHHos[114]=armHHos[86]*armHHos[114];
   armHHos[99]=1./3.*armHHos[99] + armHHos[114] + armHHos[138];
   armHHos[100]=1./3. + armHHos[100];
   armHHos[100]=1./2.*armHHos[100] + armHHos[106];
   armHHos[100]=17./2.*armHHos[48] - 17./3.*armHHos[43] + 1./3.*
   armHHos[100] + armHHos[110];
   armHHos[100]=11./9.*armHHos[9] + 1./2.*armHHos[100] + 11./3.*
   armHHos[13];
   armHHos[100]=armHHos[87]*armHHos[100];
   armHHos[99]=25./4.*armHHos[99] + armHHos[100];
   armHHos[99]=mmZ*armHHos[99];
   armHHos[100]=17*armHHos[104] - 175./12.*armHHos[15] + armHHos[121]
    + 11*armHHos[17];
   armHHos[100]=armHHos[87]*armHHos[100];
   armHHos[106]=5*armHHos[112];
   armHHos[114]=29./12. + armHHos[106];
   armHHos[114]=armHHos[9]*armHHos[114];
   armHHos[99]=armHHos[99] + 25./48.*armHHos[105] + 1./3.*armHHos[100]
    + 1./2.*armHHos[113] + 5./3.*armHHos[114];
   armHHos[100]=armHHos[168] + 7./3.*armHHos[13] + armHHos[115] + 7./2.
   *armHHos[48] + 17 - 7./3.*armHHos[43];
   armHHos[105]= - 17./3. + armHHos[183];
   armHHos[105]=armHHos[9]*armHHos[105];
   armHHos[100]=1./2.*armHHos[100] + armHHos[105];
   armHHos[100]=armHHos[87]*armHHos[100];
   armHHos[100]=25./3.*armHHos[38] + 1./2.*armHHos[100];
   armHHos[100]=mmt*armHHos[100];
   armHHos[99]=1./2.*armHHos[99] + armHHos[100];
   armHHos[99]=1./2.*armHHos[99] + armHHos[101];
   armHHos[100]=pow(armHHos[3],2);
   armHHos[99]=armHHos[100]*armHHos[99];
   armHHos[95]=armHHos[99] + armHHos[123] + 1./2.*armHHos[95] + 
   armHHos[182];
   armHHos[95]=armHHos[100]*armHHos[95];
   armHHos[99]=5./2. - 3*armHHos[42];
   armHHos[99]=1./2.*armHHos[99] + armHHos[53];
   armHHos[101]=5./4.*armHHos[47];
   armHHos[105]=9./2.*armHHos[30];
   armHHos[113]= - 3./2.*armHHos[44];
   armHHos[114]= - 9./4.*armHHos[31];
   armHHos[115]=3./2.*armHHos[12];
   armHHos[121]=27./4.*armHHos[28];
   armHHos[122]=1./2.*armHHos[7];
   armHHos[123]=1./2.*armHHos[27];
   armHHos[124]=armHHos[123] - armHHos[49];
   armHHos[138]=armHHos[124] - armHHos[16];
   armHHos[138]=armHHos[88]*armHHos[138];
   armHHos[147]=armHHos[14]*armHHos[88];
   armHHos[148]=armHHos[24]*armHHos[88];
   armHHos[99]=1./2.*armHHos[148] + armHHos[147] + 1./2.*armHHos[138]
    + armHHos[122] + armHHos[121] + 1./4.*armHHos[51] + armHHos[115] + 
   armHHos[134] + armHHos[114] + armHHos[46] + armHHos[113] - 
   armHHos[50] + armHHos[105] + 3*armHHos[99] + armHHos[101];
   armHHos[149]=9./8.*armHHos[28];
   armHHos[151]= - 1./4.*armHHos[7];
   armHHos[152]=3./4.*armHHos[119] + armHHos[151] + 1 + armHHos[149];
   armHHos[152]=armHHos[23]*armHHos[152];
   armHHos[154]=6*armHHos[37] - armHHos[40];
   armHHos[155]= - 1 - armHHos[52];
   armHHos[156]=armHHos[88]*armHHos[155];
   armHHos[158]=1./8.*armHHos[156] + armHHos[154] + 1./2.*armHHos[39];
   armHHos[158]=mmt*armHHos[158];
   armHHos[109]=1 + armHHos[109];
   armHHos[109]=armHHos[9]*armHHos[109];
   armHHos[161]= - armHHos[87]*armHHos[15]*armHHos[23];
   armHHos[163]=armHHos[38] - 9*armHHos[37] - 1./2.*armHHos[40];
   armHHos[163]=mmH*armHHos[163];
   armHHos[99]=armHHos[158] + 1./4.*armHHos[163] + 3./8.*armHHos[161]
    + 1./8.*armHHos[109] + 1./2.*armHHos[99] + armHHos[152];
   armHHos[99]=mmt*armHHos[99];
   armHHos[167]=3*armHHos[45];
   armHHos[171]= - 37./4. + armHHos[167];
   armHHos[172]=3*armHHos[42];
   armHHos[173]= - 3*armHHos[30];
   armHHos[175]=1./2.*armHHos[44];
   armHHos[179]= - 3./4.*armHHos[31];
   armHHos[182]= - 1./8.*armHHos[51];
   armHHos[183]= - 3./4.*armHHos[7];
   armHHos[184]= - 1 + 1./2.*armHHos[28];
   armHHos[184]=3*armHHos[184] + armHHos[7];
   armHHos[184]=armHHos[23]*armHHos[184];
   armHHos[185]= - 1 + 1./2.*armHHos[23];
   armHHos[185]=armHHos[9]*armHHos[185];
   armHHos[171]=1./4.*armHHos[185] + 1./4.*armHHos[184] + armHHos[183]
    - 27./8.*armHHos[28] - 1./8.*armHHos[13] + armHHos[182] - 1./2.*
   armHHos[12] + 1./8.*armHHos[43] + armHHos[179] + armHHos[175] + 
   armHHos[173] + 1./4.*armHHos[171] + armHHos[172];
   armHHos[171]=mmH*armHHos[171];
   armHHos[185]= - 3*armHHos[32];
   armHHos[186]=armHHos[185] + 3*armHHos[26] + armHHos[27];
   armHHos[186]=1./2.*armHHos[186] + armHHos[49];
   armHHos[187]= - 3./2.*armHHos[29];
   armHHos[132]= - 13./4.*armHHos[14] + armHHos[187] - 1./8.*
   armHHos[17] + 1./4.*armHHos[16] - 3*armHHos[33] + 1./8.*armHHos[25]
    + armHHos[186] + armHHos[132];
   armHHos[188]= - 1./2.*armHHos[7];
   armHHos[149]=1./4.*armHHos[147] + armHHos[188] - 2 + armHHos[149];
   armHHos[149]=armHHos[24]*armHHos[149];
   armHHos[189]=3*armHHos[29] + 7*armHHos[14];
   armHHos[189]=1./4.*armHHos[189] - armHHos[24];
   armHHos[189]=armHHos[23]*armHHos[189];
   armHHos[190]= - 9./4. + armHHos[23];
   armHHos[190]=armHHos[15]*armHHos[190];
   armHHos[99]=armHHos[99] + 1./2.*armHHos[171] + 1./8.*armHHos[104] + 
   1./4.*armHHos[190] + 1./2.*armHHos[189] + 1./2.*armHHos[132] + 
   armHHos[149];
   armHHos[99]=mmt*armHHos[99];
   armHHos[132]=5./2.*armHHos[54] - 11*armHHos[26];
   armHHos[149]= - 3*armHHos[27];
   armHHos[171]=1./4.*armHHos[55];
   armHHos[132]= - 3./4.*armHHos[16] + 33./2.*armHHos[33] + 
   armHHos[171] + 7./4.*armHHos[49] + 17./4.*armHHos[32] + 1./2.*
   armHHos[132] + armHHos[149];
   armHHos[190]=11./4.*armHHos[14] + 11./4.*armHHos[29] + armHHos[132]
    + 3./4.*armHHos[17];
   armHHos[191]=2 + armHHos[23];
   armHHos[191]=armHHos[23]*armHHos[191];
   armHHos[192]=mmt*armHHos[40];
   armHHos[126]=2*armHHos[192] + armHHos[191] + armHHos[126] + 5./8.*
   armHHos[52] + 13./4.*armHHos[31] - 5./4.*armHHos[47] + 5./2.*
   armHHos[53] - 3./2. - 2*armHHos[45];
   armHHos[126]=mmt*armHHos[126];
   armHHos[193]= - 2*armHHos[24];
   armHHos[127]=1./2.*armHHos[24] - 5./4.*armHHos[29] + armHHos[127];
   armHHos[127]=armHHos[23]*armHHos[127];
   armHHos[194]= - 3 + armHHos[23];
   armHHos[194]=armHHos[15]*armHHos[194];
   armHHos[190]=armHHos[126] + 1./8.*armHHos[194] + armHHos[127] + 1./2.
   *armHHos[190] + armHHos[193];
   armHHos[190]=mmt*armHHos[190];
   armHHos[194]= - 11*armHHos[14];
   armHHos[195]= - 15*armHHos[29] + armHHos[194];
   armHHos[196]=1./8.*armHHos[195] + armHHos[193];
   armHHos[196]=armHHos[24]*armHHos[196];
   armHHos[197]= - armHHos[15]*armHHos[24];
   armHHos[196]=armHHos[196] + 1./8.*armHHos[197];
   armHHos[190]=armHHos[196] + armHHos[190];
   armHHos[190]=mmt*armHHos[190];
   armHHos[198]= - 2*armHHos[54] + armHHos[26];
   armHHos[199]=armHHos[23]*armHHos[29];
   armHHos[200]=armHHos[134] - 9./2. - 4*armHHos[53];
   armHHos[200]=mmt*armHHos[200];
   armHHos[123]=armHHos[200] + 2*armHHos[199] + 9*armHHos[24] + 
   armHHos[14] + 2*armHHos[29] - armHHos[49] + 2*armHHos[32] + 2*
   armHHos[198] + armHHos[123];
   armHHos[200]=pow(mmt,3);
   armHHos[123]=armHHos[4]*armHHos[200]*armHHos[123];
   armHHos[190]=armHHos[190] + armHHos[123];
   armHHos[190]=armHHos[4]*armHHos[190];
   armHHos[201]=armHHos[104] - armHHos[33] + armHHos[24];
   armHHos[201]=mmH*armHHos[201];
   armHHos[202]=pow(armHHos[24],2);
   armHHos[201]=armHHos[201] + armHHos[202] + armHHos[197];
   armHHos[99]=armHHos[190] + 1./16.*armHHos[201] + armHHos[99];
   armHHos[99]=armHHos[103]*armHHos[99];
   armHHos[190]=7./2. - 9*armHHos[42];
   armHHos[101]=armHHos[115] + armHHos[134] + armHHos[114] + 
   armHHos[46] + armHHos[113] - armHHos[50] + armHHos[105] + 
   armHHos[101] + 1./2.*armHHos[190] + 3*armHHos[53];
   armHHos[105]=9./4.*armHHos[28];
   armHHos[113]=3./2.*armHHos[119];
   armHHos[190]=armHHos[113] + armHHos[188] + 2 + armHHos[105];
   armHHos[190]=armHHos[23]*armHHos[190];
   armHHos[154]=1./4.*armHHos[156] + 2*armHHos[154] + armHHos[39];
   armHHos[154]=mmt*armHHos[154];
   armHHos[156]=81./4.*armHHos[28];
   armHHos[203]=3./2.*armHHos[7];
   armHHos[204]=3./2.*armHHos[138];
   armHHos[205]=3*armHHos[147];
   armHHos[148]=3./2.*armHHos[148];
   armHHos[101]=3*armHHos[154] + 3./2.*armHHos[163] + 9./4.*
   armHHos[161] + 3./4.*armHHos[109] + 3*armHHos[190] + armHHos[148] + 
   armHHos[205] + armHHos[204] + armHHos[203] + armHHos[156] + 3*
   armHHos[101] - 13./4.*armHHos[51];
   armHHos[101]=mmt*armHHos[101];
   armHHos[154]= - 29./4. + armHHos[167];
   armHHos[154]=armHHos[179] + armHHos[175] + armHHos[173] + 1./4.*
   armHHos[154] + armHHos[172];
   armHHos[167]=5./8.*armHHos[51];
   armHHos[172]= - 81./8.*armHHos[28];
   armHHos[173]= - 9./4.*armHHos[7];
   armHHos[175]=3./4.*armHHos[184];
   armHHos[179]=1 + armHHos[129];
   armHHos[179]=armHHos[9]*armHHos[179];
   armHHos[154]=1./4.*armHHos[179] + armHHos[175] + armHHos[173] + 
   armHHos[172] - 11./8.*armHHos[13] + armHHos[167] + armHHos[133] + 3*
   armHHos[154] + 11./8.*armHHos[43];
   armHHos[154]=mmH*armHHos[154];
   armHHos[179]=3./2.*armHHos[147] + armHHos[130] - 2 + armHHos[121];
   armHHos[179]=armHHos[24]*armHHos[179];
   armHHos[184]=3*armHHos[186];
   armHHos[186]=3./4.*armHHos[16];
   armHHos[190]= - 9./2.*armHHos[29];
   armHHos[206]= - 39./4.*armHHos[14];
   armHHos[207]=45./4. - armHHos[23];
   armHHos[207]=armHHos[15]*armHHos[207];
   armHHos[101]=armHHos[101] + armHHos[154] + 3./4.*armHHos[104] + 1./2.
   *armHHos[207] + 3*armHHos[189] + armHHos[179] + armHHos[206] + 
   armHHos[190] - 27./8.*armHHos[17] + armHHos[186] - 7*armHHos[33] + 
   11./8.*armHHos[25] + armHHos[184] - 5./2.*armHHos[55];
   armHHos[101]=mmt*armHHos[101];
   armHHos[132]=3*armHHos[132];
   armHHos[154]=4*armHHos[192] + 2*armHHos[191] + armHHos[51] + 5./4.*
   armHHos[52] + 13./2.*armHHos[31] - 5./2.*armHHos[47] + 5*armHHos[53]
    - 3 - 4*armHHos[45];
   armHHos[154]=mmt*armHHos[154];
   armHHos[179]=33./4.*armHHos[29];
   armHHos[191]=33./4.*armHHos[14];
   armHHos[192]=armHHos[24] - 5./2.*armHHos[29] - 4*armHHos[14];
   armHHos[192]=armHHos[23]*armHHos[192];
   armHHos[207]= - 41 + 35*armHHos[23];
   armHHos[207]=armHHos[15]*armHHos[207];
   armHHos[154]=3*armHHos[154] + 1./4.*armHHos[207] + 3*armHHos[192] - 
   12*armHHos[24] + armHHos[191] + armHHos[179] + armHHos[132] + 41./4.
   *armHHos[17];
   armHHos[154]=mmt*armHHos[154];
   armHHos[192]=1./4.*armHHos[195] - 4*armHHos[24];
   armHHos[192]=armHHos[24]*armHHos[192];
   armHHos[192]=armHHos[192] + 1./4.*armHHos[197];
   armHHos[154]=3*armHHos[192] + armHHos[154];
   armHHos[154]=mmt*armHHos[154];
   armHHos[192]= - 2*armHHos[49];
   armHHos[195]=2*armHHos[14];
   armHHos[207]= - armHHos[52] - 9 - 8*armHHos[53];
   armHHos[207]=mmt*armHHos[207];
   armHHos[198]=armHHos[207] + 4*armHHos[199] + 18*armHHos[24] + 
   armHHos[195] + 4*armHHos[29] + armHHos[192] + 4*armHHos[32] + 4*
   armHHos[198] + armHHos[27];
   armHHos[198]=armHHos[4]*armHHos[200]*armHHos[198];
   armHHos[154]=armHHos[154] + 3*armHHos[198];
   armHHos[154]=armHHos[4]*armHHos[154];
   armHHos[198]=armHHos[15]*armHHos[24];
   armHHos[199]=armHHos[139] + armHHos[33] - armHHos[24];
   armHHos[199]=mmH*armHHos[199];
   armHHos[199]=armHHos[199] - armHHos[202] + armHHos[198];
   armHHos[101]=armHHos[154] + 5./8.*armHHos[199] + armHHos[101];
   armHHos[99]=armHHos[101] + 3*armHHos[99];
   armHHos[99]=armHHos[103]*armHHos[99];
   armHHos[154]=61./2. - 27*armHHos[42];
   armHHos[199]=27./2.*armHHos[30];
   armHHos[200]= - 3*armHHos[50];
   armHHos[148]=armHHos[148] + armHHos[205] + armHHos[204] + 
   armHHos[203] + armHHos[156] + 41./12.*armHHos[51] + 9./2.*
   armHHos[12] - 3./2.*armHHos[52] - 27./4.*armHHos[31] + armHHos[177]
    - 9./2.*armHHos[44] + armHHos[200] + armHHos[199] + 15./4.*
   armHHos[47] + 1./2.*armHHos[154] + 9*armHHos[53];
   armHHos[109]=3*armHHos[158] + 3./4.*armHHos[163] + 9./8.*
   armHHos[161] + 3./8.*armHHos[109] + 1./2.*armHHos[148] + 3*
   armHHos[152];
   armHHos[109]=mmt*armHHos[109];
   armHHos[148]= - 127./4. + 9*armHHos[45];
   armHHos[152]= - 17./3. + armHHos[129];
   armHHos[152]=armHHos[9]*armHHos[152];
   armHHos[114]=1./4.*armHHos[152] + armHHos[175] + armHHos[173] + 
   armHHos[172] + 7./24.*armHHos[13] - 25./24.*armHHos[51] + 
   armHHos[133] - 7./24.*armHHos[43] + armHHos[114] + armHHos[131] - 9*
   armHHos[30] + 1./4.*armHHos[148] + 9*armHHos[42];
   armHHos[114]=mmH*armHHos[114];
   armHHos[131]=armHHos[206] + armHHos[190] + 13./8.*armHHos[17] + 
   armHHos[186] - 31./3.*armHHos[33] - 7./24.*armHHos[25] + 
   armHHos[184] + 25./6.*armHHos[55];
   armHHos[148]=27./8.*armHHos[28];
   armHHos[152]= - 3./2.*armHHos[7];
   armHHos[154]=3./4.*armHHos[147] + armHHos[152] - 28./3. + 
   armHHos[148];
   armHHos[154]=armHHos[24]*armHHos[154];
   armHHos[158]= - 75./4. + 17./3.*armHHos[23];
   armHHos[158]=armHHos[15]*armHHos[158];
   armHHos[109]=armHHos[109] + 1./2.*armHHos[114] + 3./8.*armHHos[104]
    + 1./4.*armHHos[158] + 3./2.*armHHos[189] + 1./2.*armHHos[131] + 
   armHHos[154];
   armHHos[109]=mmt*armHHos[109];
   armHHos[114]=armHHos[191] + armHHos[179] + armHHos[132] - 37./12.*
   armHHos[17];
   armHHos[131]=37 - 55*armHHos[23];
   armHHos[131]=armHHos[15]*armHHos[131];
   armHHos[114]=3*armHHos[126] + 1./24.*armHHos[131] + 3*armHHos[127]
    + 1./2.*armHHos[114] - 6*armHHos[24];
   armHHos[114]=mmt*armHHos[114];
   armHHos[114]=3*armHHos[196] + armHHos[114];
   armHHos[114]=mmt*armHHos[114];
   armHHos[114]=armHHos[114] + 3*armHHos[123];
   armHHos[114]=armHHos[4]*armHHos[114];
   armHHos[109]=armHHos[114] + 25./48.*armHHos[201] + armHHos[109];
   armHHos[109]=armHHos[100]*armHHos[109];
   armHHos[101]=armHHos[101] + armHHos[109];
   armHHos[101]=armHHos[100]*armHHos[101];
   armHHos[109]=3./2.*armHHos[14];
   armHHos[114]= - armHHos[12] + armHHos[52] + 1./4. + armHHos[44];
   armHHos[114]=mmH*armHHos[114];
   armHHos[123]=mmt*armHHos[155];
   armHHos[114]=1./2.*armHHos[123] + 1./2.*armHHos[114] + armHHos[24]
    + armHHos[109] + armHHos[124] - 1./2.*armHHos[16];
   armHHos[114]=mmt*armHHos[114];
   armHHos[123]=armHHos[24]*armHHos[14];
   armHHos[124]=armHHos[24]*armHHos[142];
   armHHos[124]=armHHos[33] + armHHos[124];
   armHHos[124]=mmH*armHHos[124];
   armHHos[123]=armHHos[123] + armHHos[124];
   armHHos[114]=1./2.*armHHos[123] + armHHos[114];
   armHHos[114]=mmt*armHHos[114];
   armHHos[123]=armHHos[100]*armHHos[114];
   armHHos[123]=armHHos[114] + 1./2.*armHHos[123];
   armHHos[123]=armHHos[100]*armHHos[123];
   armHHos[123]=3./2.*armHHos[114] + armHHos[123];
   armHHos[123]=armHHos[100]*armHHos[123];
   armHHos[124]=armHHos[103]*armHHos[114];
   armHHos[114]=3*armHHos[114] + armHHos[124];
   armHHos[114]=armHHos[103]*armHHos[114];
   armHHos[114]=armHHos[123] + 1./2.*armHHos[114];
   armHHos[114]=armHHos[8]*armHHos[114];
   armHHos[99]=3*armHHos[114] + armHHos[101] + armHHos[99];
   armHHos[99]=armHHos[8]*armHHos[99];
   armHHos[101]= - 2*armHHos[35] - armHHos[34];
   armHHos[101]=11*armHHos[41] + 32./3.*armHHos[101] + armHHos[170];
   armHHos[101]=mmt*armHHos[101];
   armHHos[114]= - 27*armHHos[47];
   armHHos[123]=929./6. + armHHos[114];
   armHHos[124]=11*armHHos[48];
   armHHos[126]=25./3.*armHHos[23] + 158./3. + 21./4.*armHHos[7];
   armHHos[126]=armHHos[23]*armHHos[126];
   armHHos[127]=3./4.*armHHos[23];
   armHHos[131]= - 1 + armHHos[127];
   armHHos[131]=armHHos[9]*armHHos[131];
   armHHos[132]=3*armHHos[131];
   armHHos[154]= - 3./4.*armHHos[44];
   armHHos[155]=13./8.*armHHos[52];
   armHHos[158]=3./4.*armHHos[12];
   armHHos[101]=armHHos[101] + armHHos[132] + armHHos[126] + 
   armHHos[130] + 155./6.*armHHos[51] + armHHos[158] + armHHos[124] + 
   armHHos[155] + 149./12.*armHHos[31] - 11*armHHos[46] + armHHos[154]
    + 1./4.*armHHos[123] + armHHos[200];
   armHHos[101]=mmt*armHHos[101];
   armHHos[123]=2./3.*armHHos[36] + armHHos[149];
   armHHos[126]=21*armHHos[7];
   armHHos[149]= - 457./6. + armHHos[126];
   armHHos[149]=armHHos[24]*armHHos[149];
   armHHos[161]= - 9./4.*armHHos[14] - 10*armHHos[24];
   armHHos[161]=armHHos[23]*armHHos[161];
   armHHos[163]=11./4.*armHHos[49];
   armHHos[170]= - 7./2.*armHHos[16];
   armHHos[173]=13*armHHos[14];
   armHHos[101]=armHHos[101] + 3*armHHos[139] - 119./6.*armHHos[15] + 
   armHHos[161] + 1./4.*armHHos[149] + armHHos[173] + 149./12.*
   armHHos[17] + armHHos[170] + 679./24.*armHHos[33] - 177./8.*
   armHHos[25] + 58./3.*armHHos[55] + 2*armHHos[123] + armHHos[163];
   armHHos[101]=mmt*armHHos[101];
   armHHos[123]=22 + armHHos[168];
   armHHos[123]=armHHos[15]*armHHos[123];
   armHHos[149]= - 53 + armHHos[200];
   armHHos[149]= - 22*armHHos[51] + 1./2.*armHHos[149] - 3*armHHos[52];
   armHHos[149]=mmt*armHHos[149];
   armHHos[161]=armHHos[27] + armHHos[192];
   armHHos[168]= - armHHos[23]*armHHos[14];
   armHHos[123]=armHHos[149] + armHHos[123] + 9*armHHos[168] + 50*
   armHHos[24] + armHHos[164] + 11*armHHos[25] + 3*armHHos[161] - 22*
   armHHos[55];
   armHHos[123]=mmt*armHHos[123];
   armHHos[149]= - armHHos[24]*armHHos[14];
   armHHos[164]=armHHos[149] + armHHos[197];
   armHHos[123]=36*armHHos[164] + armHHos[123];
   armHHos[123]=armHHos[4]*mmt*armHHos[123];
   armHHos[164]=41*armHHos[14];
   armHHos[175]=armHHos[164] - 265./3.*armHHos[24];
   armHHos[175]=armHHos[24]*armHHos[175];
   armHHos[175]=armHHos[175] + 41./3.*armHHos[197];
   armHHos[101]=armHHos[123] + 1./4.*armHHos[175] + armHHos[101];
   armHHos[101]=armHHos[4]*armHHos[101];
   armHHos[123]= - armHHos[9]*armHHos[23];
   armHHos[175]=9*armHHos[123];
   armHHos[177]= - 1 - armHHos[51];
   armHHos[179]=armHHos[175] + 11*armHHos[177] + armHHos[137];
   armHHos[179]=armHHos[87]*armHHos[179];
   armHHos[184]= - 1 + armHHos[200];
   armHHos[134]= - armHHos[12] + armHHos[134] + 1./2.*armHHos[184] + 
   armHHos[44];
   armHHos[134]=armHHos[88]*armHHos[134];
   armHHos[184]=2./3.*armHHos[35] + armHHos[34];
   armHHos[186]=1 - armHHos[7];
   armHHos[189]=armHHos[88]*armHHos[186];
   armHHos[190]=armHHos[23]*armHHos[189];
   armHHos[179]=1./4.*armHHos[179] + 9./4.*armHHos[190] + 3./4.*
   armHHos[134] - 6*armHHos[38] - 27./4.*armHHos[41] + 8*armHHos[184]
    - 9./2.*armHHos[39];
   armHHos[179]=mmt*armHHos[179];
   armHHos[184]=21*armHHos[50] - 1211./12. - 9*armHHos[47];
   armHHos[97]= - 19./2.*armHHos[46] + 1./2.*armHHos[184] + armHHos[97]
   ;
   armHHos[97]=armHHos[176] - 83./12.*armHHos[51] + armHHos[115] - 9*
   armHHos[48] - 5./2.*armHHos[43] - armHHos[52] + 1./2.*armHHos[97] + 
   7*armHHos[31];
   armHHos[127]=2 + armHHos[127];
   armHHos[127]=armHHos[15]*armHHos[127];
   armHHos[166]=armHHos[166] - armHHos[17];
   armHHos[166]=1./2.*armHHos[166] + armHHos[24];
   armHHos[127]=11./2.*armHHos[166] + 3*armHHos[127];
   armHHos[127]=armHHos[87]*armHHos[127];
   armHHos[176]=armHHos[205] - 21./2. - armHHos[7];
   armHHos[176]=3./4.*armHHos[176] - 2*armHHos[23];
   armHHos[176]=armHHos[23]*armHHos[176];
   armHHos[184]= - 1 + armHHos[169];
   armHHos[184]=armHHos[9]*armHHos[184];
   armHHos[191]= - 4./3.*armHHos[34] + armHHos[41];
   armHHos[191]=mmH*armHHos[191];
   armHHos[138]=3./4.*armHHos[138];
   armHHos[192]=3./4.*armHHos[24]*armHHos[88]*armHHos[7];
   armHHos[97]=armHHos[179] + armHHos[191] + armHHos[127] + 1./2.*
   armHHos[184] + armHHos[176] + armHHos[192] + armHHos[138] + 1./2.*
   armHHos[97] + armHHos[125];
   armHHos[97]=mmt*armHHos[97];
   armHHos[127]= - 65./8. + armHHos[11];
   armHHos[176]= - 5./8.*armHHos[13];
   armHHos[179]=armHHos[86]*armHHos[33];
   armHHos[127]=5./2.*armHHos[179] + armHHos[176] + armHHos[167] + 1./3.
   *armHHos[127] + 5./8.*armHHos[43];
   armHHos[167]=5./4.*armHHos[116];
   armHHos[127]=5./4.*armHHos[157] + armHHos[167] + 1./2.*armHHos[127]
    - armHHos[7];
   armHHos[127]=mmH*armHHos[127];
   armHHos[157]=armHHos[167] + 3*armHHos[119] + 175./24. + armHHos[130]
   ;
   armHHos[157]=armHHos[24]*armHHos[157];
   armHHos[106]=armHHos[129] - 3./4. + armHHos[106];
   armHHos[106]=armHHos[15]*armHHos[106];
   armHHos[129]=armHHos[14] - 1./2.*armHHos[24];
   armHHos[129]=armHHos[23]*armHHos[129];
   armHHos[167]=armHHos[87]*armHHos[198];
   armHHos[179]=9./4.*armHHos[16];
   armHHos[184]= - 17./4.*armHHos[14];
   armHHos[97]=armHHos[101] + armHHos[97] + 1./2.*armHHos[127] + 1./2.*
   armHHos[167] + 2*armHHos[104] + 1./8.*armHHos[106] + 9./4.*
   armHHos[129] + 1./2.*armHHos[157] + armHHos[184] + 1./32.*
   armHHos[17] + armHHos[179] - 79./16.*armHHos[33] + 5./32.*
   armHHos[25] + armHHos[49] - 1./6.*armHHos[55];
   armHHos[101]=31 + armHHos[114];
   armHHos[104]=27./2. + 7*armHHos[7];
   armHHos[104]=1./2.*armHHos[104] + armHHos[23];
   armHHos[104]=armHHos[23]*armHHos[104];
   armHHos[106]=armHHos[39] + 1./2.*armHHos[41];
   armHHos[106]=mmt*armHHos[106];
   armHHos[101]=3*armHHos[106] + 3./2.*armHHos[131] + 3./2.*
   armHHos[104] + armHHos[130] + 13./4.*armHHos[51] + armHHos[158] + 
   armHHos[180] + armHHos[155] + 33./4.*armHHos[31] - 3./2.*armHHos[46]
    + armHHos[154] + 1./4.*armHHos[101] + armHHos[200];
   armHHos[101]=mmt*armHHos[101];
   armHHos[104]= - 25 + armHHos[126];
   armHHos[104]=armHHos[24]*armHHos[104];
   armHHos[101]=armHHos[101] + armHHos[140] + 17./4.*armHHos[15] + 9./4.
   *armHHos[168] + 1./4.*armHHos[104] + armHHos[173] - 5./8.*
   armHHos[17] + armHHos[170] + 41./4.*armHHos[33] - 57./16.*
   armHHos[25] + 5./2.*armHHos[55] - 6*armHHos[27] + armHHos[163];
   armHHos[101]=mmt*armHHos[101];
   armHHos[104]=1 + armHHos[169];
   armHHos[104]=armHHos[15]*armHHos[104];
   armHHos[106]= - 5 - armHHos[50];
   armHHos[106]= - armHHos[51] + 1./2.*armHHos[106] - armHHos[52];
   armHHos[106]=mmt*armHHos[106];
   armHHos[104]=armHHos[106] + armHHos[104] + 3*armHHos[168] + 4*
   armHHos[24] + armHHos[195] + armHHos[165] + armHHos[161] - 
   armHHos[55];
   armHHos[104]=mmt*armHHos[104];
   armHHos[106]=2*armHHos[149] + armHHos[197];
   armHHos[104]=6*armHHos[106] + armHHos[104];
   armHHos[104]=armHHos[4]*mmt*armHHos[104];
   armHHos[106]=armHHos[164] - 23./2.*armHHos[24];
   armHHos[106]=armHHos[24]*armHHos[106];
   armHHos[106]=armHHos[106] + 41./2.*armHHos[198];
   armHHos[101]=3*armHHos[104] + 1./4.*armHHos[106] + armHHos[101];
   armHHos[101]=armHHos[4]*armHHos[101];
   armHHos[104]=7*armHHos[50] + 41./4. - 3*armHHos[47];
   armHHos[104]=15./4.*armHHos[31] - 9./4.*armHHos[46] + 1./2.*
   armHHos[104] - armHHos[44];
   armHHos[106]= - 3./4.*armHHos[13];
   armHHos[104]=armHHos[106] + armHHos[182] + armHHos[115] + 
   armHHos[180] + 3./4.*armHHos[43] + 3./2.*armHHos[104] - armHHos[52];
   armHHos[114]=3*armHHos[123] + armHHos[177] + 3*armHHos[23];
   armHHos[114]=armHHos[87]*armHHos[114];
   armHHos[115]= - 3*armHHos[39] - 1./4.*armHHos[41];
   armHHos[114]=1./8.*armHHos[114] + 3./4.*armHHos[190] + 1./4.*
   armHHos[134] + 1./2.*armHHos[115] - armHHos[38];
   armHHos[114]=mmt*armHHos[114];
   armHHos[115]=armHHos[15]*armHHos[23];
   armHHos[115]=armHHos[166] + 3./2.*armHHos[115];
   armHHos[115]=armHHos[87]*armHHos[115];
   armHHos[123]=armHHos[205] - 9./4. - armHHos[7];
   armHHos[123]=armHHos[23]*armHHos[123];
   armHHos[126]=1 - 1./2.*armHHos[23];
   armHHos[126]=armHHos[9]*armHHos[126];
   armHHos[104]=3*armHHos[114] + 3./4.*armHHos[115] + 3./4.*
   armHHos[126] + 3./4.*armHHos[123] + armHHos[192] + armHHos[138] + 1./
   2.*armHHos[104] + armHHos[125];
   armHHos[104]=mmt*armHHos[104];
   armHHos[114]= - armHHos[86]*armHHos[33];
   armHHos[115]=3./2.*armHHos[114] + 3./8.*armHHos[13] - 3./8.*
   armHHos[51] + armHHos[118] - 47./8. + armHHos[11];
   armHHos[116]=3*armHHos[116];
   armHHos[118]= - 13./4. + armHHos[116];
   armHHos[118]=armHHos[9]*armHHos[118];
   armHHos[115]=1./8.*armHHos[118] + 3./8.*armHHos[112] + 1./4.*
   armHHos[115] - armHHos[7];
   armHHos[115]=mmH*armHHos[115];
   armHHos[118]=1./8.*armHHos[112] + armHHos[119] + 7./16. - armHHos[7]
   ;
   armHHos[118]=armHHos[24]*armHHos[118];
   armHHos[123]=armHHos[14] - 1./4.*armHHos[24];
   armHHos[123]=armHHos[23]*armHHos[123];
   armHHos[116]=27./2.*armHHos[23] - 91./4. + armHHos[116];
   armHHos[116]=armHHos[15]*armHHos[116];
   armHHos[126]=armHHos[87]*armHHos[197];
   armHHos[101]=armHHos[101] + armHHos[104] + 1./2.*armHHos[115] + 3./4.
   *armHHos[126] + 1./16.*armHHos[116] + 9./4.*armHHos[123] + 3./2.*
   armHHos[118] + armHHos[184] + 57./64.*armHHos[17] + armHHos[179] - 
   63./32.*armHHos[33] - 3./64.*armHHos[25] + armHHos[49] + 
   armHHos[171];
   armHHos[101]=armHHos[103]*armHHos[101];
   armHHos[101]=armHHos[97] + armHHos[101];
   armHHos[101]=armHHos[103]*armHHos[101];
   armHHos[104]= - 19./8. + armHHos[11];
   armHHos[115]=5./8.*armHHos[13];
   armHHos[104]=5./2.*armHHos[112] + 5./2.*armHHos[114] + armHHos[115]
    - 5./8.*armHHos[51] + 1./3.*armHHos[104] - 5./8.*armHHos[43];
   armHHos[114]= - 101./12. + armHHos[128];
   armHHos[114]=armHHos[9]*armHHos[114];
   armHHos[104]=5*armHHos[104] + 1./2.*armHHos[114];
   armHHos[104]=mmH*armHHos[104];
   armHHos[114]=25*armHHos[43] + 97./2.*armHHos[31] + 287./3. - 97./2.*
   armHHos[46];
   armHHos[116]=17./3. + armHHos[169];
   armHHos[116]=armHHos[9]*armHHos[116];
   armHHos[114]=armHHos[116] + 23./12.*armHHos[23] - 25./6.*armHHos[13]
    + 151./36.*armHHos[51] + 1./6.*armHHos[114] + armHHos[124];
   armHHos[116]=1 + armHHos[51];
   armHHos[118]=armHHos[175] + 7./3.*armHHos[116] + armHHos[137];
   armHHos[118]=armHHos[87]*armHHos[118];
   armHHos[118]=1./8.*armHHos[118] + 13./8.*armHHos[41] - 3*armHHos[38]
   ;
   armHHos[118]=mmt*armHHos[118];
   armHHos[123]=armHHos[98] + armHHos[17];
   armHHos[123]=1./2.*armHHos[123] - armHHos[24];
   armHHos[124]= - 2 + 9./8.*armHHos[23];
   armHHos[124]=armHHos[15]*armHHos[124];
   armHHos[123]=7./12.*armHHos[123] + armHHos[124];
   armHHos[123]=armHHos[87]*armHHos[123];
   armHHos[124]= - mmH*armHHos[41];
   armHHos[114]=armHHos[118] + 1./3.*armHHos[124] + 1./4.*armHHos[114]
    + armHHos[123];
   armHHos[114]=mmt*armHHos[114];
   armHHos[118]= - 67./4. + 17*armHHos[23];
   armHHos[118]=armHHos[23]*armHHos[118];
   armHHos[118]=armHHos[118] - 115./6.*armHHos[51] - 7*armHHos[48] - 7*
   armHHos[31] - 235./12. + 7*armHHos[46];
   armHHos[123]= - mmt*armHHos[41];
   armHHos[118]=7./3.*armHHos[123] + 1./3.*armHHos[118] + armHHos[132];
   armHHos[118]=mmt*armHHos[118];
   armHHos[123]=13./12.*armHHos[24] - 373./36.*armHHos[17] - 109./36.*
   armHHos[33] - 41./9.*armHHos[55] + 23./8.*armHHos[25];
   armHHos[118]=1./2.*armHHos[118] + armHHos[140] + 493./36.*
   armHHos[15] + 1./2.*armHHos[123] + 2*armHHos[141];
   armHHos[118]=mmt*armHHos[118];
   armHHos[98]=armHHos[98] + armHHos[193];
   armHHos[120]= - 7./3. + armHHos[120];
   armHHos[120]=armHHos[15]*armHHos[120];
   armHHos[116]=mmt*armHHos[116];
   armHHos[98]=7./3.*armHHos[116] + 7./3.*armHHos[98] + armHHos[120];
   armHHos[98]=mmt*armHHos[98];
   armHHos[98]=18*armHHos[197] + armHHos[98];
   armHHos[98]=armHHos[4]*mmt*armHHos[98];
   armHHos[116]= - 23*armHHos[202] + 41*armHHos[198];
   armHHos[98]=armHHos[98] + 17./72.*armHHos[116] + armHHos[118];
   armHHos[98]=armHHos[4]*armHHos[98];
   armHHos[112]= - 227./6. + 25*armHHos[112];
   armHHos[112]=armHHos[24]*armHHos[112];
   armHHos[112]=1./4.*armHHos[112] + 283./16.*armHHos[17] + 35./8.*
   armHHos[33] + 17./3.*armHHos[55] - 25./16.*armHHos[25];
   armHHos[116]= - 1283./12. + armHHos[128];
   armHHos[116]=1./3.*armHHos[116] + 43./2.*armHHos[23];
   armHHos[116]=armHHos[15]*armHHos[116];
   armHHos[112]=1./4.*armHHos[116] + 1./3.*armHHos[112] + 9./4.*
   armHHos[141];
   armHHos[98]=armHHos[98] + armHHos[114] + 1./24.*armHHos[104] + 17./
   12.*armHHos[126] + 1./4.*armHHos[112] + 2./3.*armHHos[139];
   armHHos[98]=armHHos[100]*armHHos[98];
   armHHos[97]=armHHos[97] + armHHos[98];
   armHHos[97]=armHHos[100]*armHHos[97];
   armHHos[97]=armHHos[99] + armHHos[97] + armHHos[101];
   armHHos[97]=armHHos[8]*armHHos[97];
   armHHos[98]=3 - armHHos[23];
   armHHos[98]=mmt*armHHos[23]*armHHos[98];
   armHHos[99]=1./2.*armHHos[98] + armHHos[24] + 1./2.*armHHos[141];
   armHHos[101]=pow(mmt,2);
   armHHos[99]=armHHos[101]*armHHos[99];
   armHHos[104]= - 1 - armHHos[23];
   armHHos[104]=mmt*armHHos[23]*armHHos[104];
   armHHos[104]=armHHos[174] + armHHos[104];
   armHHos[104]=mmt*armHHos[104];
   armHHos[104]=2*armHHos[202] + armHHos[104];
   armHHos[104]=armHHos[4]*armHHos[101]*armHHos[104];
   armHHos[99]=armHHos[99] + 2*armHHos[104];
   armHHos[99]=armHHos[4]*armHHos[99];
   armHHos[112]= - 1 + armHHos[23];
   armHHos[112]=armHHos[101]*armHHos[23]*armHHos[112];
   armHHos[99]=1./4.*armHHos[112] + armHHos[99];
   armHHos[114]=armHHos[100]*armHHos[99];
   armHHos[98]=armHHos[98] + armHHos[136] + armHHos[141];
   armHHos[98]=armHHos[101]*armHHos[98];
   armHHos[98]=armHHos[98] + 4*armHHos[104];
   armHHos[98]=armHHos[4]*armHHos[98];
   armHHos[98]=1./2.*armHHos[112] + armHHos[98];
   armHHos[101]=armHHos[98] + armHHos[114];
   armHHos[101]=armHHos[100]*armHHos[101];
   armHHos[99]=armHHos[103]*armHHos[99];
   armHHos[98]=armHHos[98] + armHHos[99];
   armHHos[98]=armHHos[103]*armHHos[98];
   armHHos[98]=armHHos[101] + armHHos[98];
   armHHos[98]=armHHos[22]*armHHos[98]*pow(armHHos[8],2);
   armHHos[99]=1 + armHHos[50];
   armHHos[101]=pow(mmZ,2);
   armHHos[99]=armHHos[101]*armHHos[99]*pow(armHHos[4],2);
   armHHos[90]=9*armHHos[98] + armHHos[97] + armHHos[90] + 3*
   armHHos[99] + armHHos[95];
   armHHos[90]=armHHos[22]*armHHos[90];
   armHHos[95]= - 161./2. - 9*armHHos[75];
   armHHos[97]=27*armHHos[72];
   armHHos[95]=1./4.*armHHos[95] + armHHos[97];
   armHHos[98]=9*armHHos[28];
   armHHos[99]=3./4.*armHHos[9];
   armHHos[104]=armHHos[99] + armHHos[7] + 11./4. + armHHos[98];
   armHHos[104]=armHHos[9]*armHHos[104];
   armHHos[112]= - 27./2.*armHHos[30];
   armHHos[114]= - 81./4.*armHHos[28];
   armHHos[95]=3./4.*armHHos[104] + armHHos[183] + armHHos[114] + 33./8.
   *armHHos[13] + 27./16.*armHHos[79] + armHHos[112] + 1./16.*
   armHHos[80] - 3*armHHos[74] + 1./2.*armHHos[95] + 18*armHHos[78];
   armHHos[95]=armHHos[87]*armHHos[95];
   armHHos[104]=25./2. - armHHos[77];
   armHHos[104]=armHHos[152] + 1./2.*armHHos[12] + 9./2.*armHHos[79] + 
   1./2.*armHHos[104] + armHHos[80];
   armHHos[104]=armHHos[88]*armHHos[104];
   armHHos[116]= - 3*armHHos[63] - 1./4.*armHHos[66];
   armHHos[118]=2*armHHos[65];
   armHHos[120]=9./16.*armHHos[9]*armHHos[143];
   armHHos[123]= - 1./4.*armHHos[68];
   armHHos[95]=armHHos[95] + armHHos[120] + 1./4.*armHHos[104] + 
   armHHos[118] + 3*armHHos[116] + armHHos[123];
   armHHos[95]=mmZ*armHHos[95];
   armHHos[104]=armHHos[84] - 1./2.*armHHos[58];
   armHHos[116]=3./4.*armHHos[104] + armHHos[82];
   armHHos[124]= - 7./4.*armHHos[60];
   armHHos[126]= - 5./4.*armHHos[32];
   armHHos[116]= - 7./2.*armHHos[14] - armHHos[29] + 143./4.*
   armHHos[17] + armHHos[126] + 3*armHHos[116] + armHHos[124];
   armHHos[127]=1./4.*armHHos[7];
   armHHos[128]=armHHos[127] - 7 + armHHos[121];
   armHHos[128]=armHHos[15]*armHHos[128];
   armHHos[129]= - 7./4.*armHHos[29];
   armHHos[131]=43./8.*armHHos[15];
   armHHos[132]=armHHos[131] + armHHos[129] - armHHos[14];
   armHHos[132]=1./2.*armHHos[9]*armHHos[132];
   armHHos[116]=armHHos[132] + 1./2.*armHHos[116] + 3*armHHos[128];
   armHHos[116]=armHHos[87]*armHHos[116];
   armHHos[128]= - 277./6. + 3*armHHos[75];
   armHHos[128]=17*armHHos[77] + 1./2.*armHHos[128] + armHHos[97];
   armHHos[134]=13./3.*armHHos[78];
   armHHos[128]=1./2.*armHHos[128] + armHHos[134];
   armHHos[136]=3./8.*armHHos[147];
   armHHos[137]=5./8.*armHHos[9];
   armHHos[138]=armHHos[137] + armHHos[136] - 15./8. - armHHos[7];
   armHHos[138]=armHHos[9]*armHHos[138];
   armHHos[139]= - 7 + armHHos[7];
   armHHos[139]=armHHos[7]*armHHos[139];
   armHHos[140]= - 59 + armHHos[7];
   armHHos[140]=armHHos[14]*armHHos[88]*armHHos[140];
   armHHos[141]= - 17 - armHHos[7];
   armHHos[141]=armHHos[15]*armHHos[88]*armHHos[141];
   armHHos[143]=3*armHHos[63] + 1./2.*armHHos[66];
   armHHos[143]=armHHos[65] + 3*armHHos[143] + 7./4.*armHHos[68];
   armHHos[143]=mmH*armHHos[143];
   armHHos[149]=armHHos[104] + armHHos[16];
   armHHos[149]=armHHos[88]*armHHos[149];
   armHHos[152]= - 99./8.*armHHos[28];
   armHHos[154]= - 27./4.*armHHos[30];
   armHHos[95]=armHHos[95] + 1./2.*armHHos[143] + armHHos[116] + 
   armHHos[138] + 1./8.*armHHos[141] + 1./8.*armHHos[140] + 9./4.*
   armHHos[149] + 1./4.*armHHos[139] + armHHos[152] + armHHos[115] - 17.
   /4.*armHHos[12] + 9./8.*armHHos[79] + armHHos[154] + 77./24.*
   armHHos[80] + 1./2.*armHHos[128] - armHHos[74];
   armHHos[95]=armHHos[57]*armHHos[95];
   armHHos[115]=81*armHHos[28];
   armHHos[116]= - 233./3. + armHHos[115];
   armHHos[116]=1./2.*armHHos[116] + armHHos[125];
   armHHos[116]=armHHos[15]*armHHos[116];
   armHHos[128]=armHHos[131] + armHHos[129] + armHHos[14];
   armHHos[128]=armHHos[9]*armHHos[128];
   armHHos[138]=1./4.*armHHos[104] + 3*armHHos[82];
   armHHos[116]=armHHos[128] + armHHos[116] - 3./2.*armHHos[14] - 
   armHHos[29] + 135./4.*armHHos[17] + armHHos[126] + armHHos[138] + 
   armHHos[124];
   armHHos[116]=armHHos[87]*armHHos[116];
   armHHos[128]= - 3./2. - armHHos[75];
   armHHos[128]=1./2.*armHHos[128] + 3*armHHos[72];
   armHHos[128]=3*armHHos[128] + armHHos[77];
   armHHos[128]=3./2.*armHHos[128] + armHHos[134];
   armHHos[139]=armHHos[137] + 9./8.*armHHos[147] - 125./24. + 
   armHHos[7];
   armHHos[139]=armHHos[9]*armHHos[139];
   armHHos[140]= - 1./2.*armHHos[66];
   armHHos[141]=armHHos[65] + armHHos[123] + 9*armHHos[63] + 
   armHHos[140];
   armHHos[141]=mmH*armHHos[141];
   armHHos[143]= - armHHos[15]*armHHos[88];
   armHHos[113]=1./2.*armHHos[141] + 1./2.*armHHos[116] + armHHos[139]
    + 1./4.*armHHos[143] + armHHos[113] + 1./4.*armHHos[149] - 5./4.*
   armHHos[7] + armHHos[152] + 17./8.*armHHos[13] - 3./4.*armHHos[12]
    + 37./12.*armHHos[79] + armHHos[154] + 3./4.*armHHos[80] + 1./2.*
   armHHos[128] - armHHos[74];
   armHHos[116]= - 233./3. - armHHos[75];
   armHHos[116]=1./4.*armHHos[116] + armHHos[97];
   armHHos[128]=27*armHHos[28];
   armHHos[139]=79./12. + armHHos[128];
   armHHos[139]=9./8.*armHHos[9] + 1./2.*armHHos[139] + armHHos[125];
   armHHos[139]=armHHos[9]*armHHos[139];
   armHHos[141]=9*armHHos[78];
   armHHos[143]= - 3./2.*armHHos[74];
   armHHos[116]=1./4.*armHHos[139] + armHHos[183] + armHHos[172] + 25./
   16.*armHHos[13] + 3./32.*armHHos[79] + armHHos[154] + armHHos[143]
    + 1./4.*armHHos[116] + armHHos[141];
   armHHos[116]=armHHos[87]*armHHos[116];
   armHHos[139]= - 9*armHHos[63];
   armHHos[149]=armHHos[139] + armHHos[140];
   armHHos[155]=1 + armHHos[79];
   armHHos[157]=armHHos[88]*armHHos[155];
   armHHos[116]=armHHos[116] + 1./16.*armHHos[157] + 1./2.*armHHos[149]
    + armHHos[65];
   armHHos[116]=mmZ*armHHos[116];
   armHHos[113]=1./2.*armHHos[113] + armHHos[116];
   armHHos[113]=armHHos[57]*armHHos[113];
   armHHos[116]=armHHos[29] + armHHos[14];
   armHHos[149]= - 23./2.*armHHos[15];
   armHHos[116]=19*armHHos[116] + armHHos[149];
   armHHos[116]=armHHos[9]*armHHos[116];
   armHHos[157]= - 17./3.*armHHos[82];
   armHHos[158]=37./2.*armHHos[60];
   armHHos[161]= - 77./6.*armHHos[32];
   armHHos[163]= - 43./3.*armHHos[29];
   armHHos[116]=1./2.*armHHos[116] + 5*armHHos[15] - 19./3.*armHHos[14]
    + armHHos[163] - 233./6.*armHHos[17] + 5./3.*armHHos[16] + 
   armHHos[161] + armHHos[158] + armHHos[157] - 7./3.*armHHos[84] + 13./
   4.*armHHos[58];
   armHHos[164]=33 + 5*armHHos[75];
   armHHos[165]= - 23./3. - 17*armHHos[9];
   armHHos[165]=armHHos[9]*armHHos[165];
   armHHos[166]= - mmZ*armHHos[65];
   armHHos[167]= - 4./3.*armHHos[78];
   armHHos[168]=3./2.*armHHos[74];
   armHHos[164]=3./2.*armHHos[166] + 1./32.*armHHos[165] - 29./16.*
   armHHos[13] - 19./96.*armHHos[79] + armHHos[168] + 1./16.*
   armHHos[164] + armHHos[167];
   armHHos[164]=mmZ*armHHos[164];
   armHHos[116]=1./8.*armHHos[116] + armHHos[164];
   armHHos[116]=armHHos[57]*armHHos[116];
   armHHos[164]=3./4.*armHHos[15];
   armHHos[165]= - armHHos[29] + armHHos[164];
   armHHos[165]=armHHos[9]*armHHos[165];
   armHHos[166]=1 - armHHos[9];
   armHHos[166]=armHHos[9]*armHHos[166];
   armHHos[169]=9./8.*armHHos[166] + 1./8.*armHHos[79] + 31./8. + 3*
   armHHos[78];
   armHHos[169]=mmZ*armHHos[169];
   armHHos[170]= - 3./2.*armHHos[32];
   armHHos[138]=armHHos[169] + 3./2.*armHHos[165] - 4*armHHos[15] - 1./
   2.*armHHos[14] + armHHos[187] + armHHos[170] + armHHos[138] - 3./2.*
   armHHos[60];
   armHHos[138]=mmZ*armHHos[138];
   armHHos[165]=pow(armHHos[15],2);
   armHHos[138]=9./4.*armHHos[165] + armHHos[138];
   armHHos[138]=armHHos[57]*armHHos[138];
   armHHos[165]= - armHHos[11] - 1 - armHHos[20];
   armHHos[169]=armHHos[56]*armHHos[165];
   armHHos[165]=armHHos[1]*armHHos[165];
   armHHos[165]=11./3.*armHHos[169] + 3*armHHos[165];
   armHHos[165]=mmZ*armHHos[165];
   armHHos[169]= - armHHos[56]*armHHos[6];
   armHHos[171]= - armHHos[1]*armHHos[6];
   armHHos[173]=11./3.*armHHos[56] + 3*armHHos[1];
   armHHos[173]=armHHos[15]*armHHos[173];
   armHHos[165]=armHHos[165] + armHHos[173] + 11./3.*armHHos[169] + 3*
   armHHos[171];
   armHHos[165]=mmZ*armHHos[165];
   armHHos[138]=1./2.*armHHos[165] + armHHos[138];
   armHHos[138]=armHHos[4]*armHHos[138];
   armHHos[165]=11./2.*armHHos[6] - 5*armHHos[17];
   armHHos[169]=armHHos[56]*armHHos[165];
   armHHos[165]=armHHos[1]*armHHos[165];
   armHHos[171]= - 11./9.*armHHos[56] - armHHos[1];
   armHHos[173]=armHHos[15]*armHHos[171];
   armHHos[165]=13*armHHos[173] + 11./3.*armHHos[169] + 3*armHHos[165];
   armHHos[169]=7./3. + armHHos[20];
   armHHos[169]=1./2.*armHHos[169];
   armHHos[145]=armHHos[169] + armHHos[145];
   armHHos[145]= - 5./2.*armHHos[13] + 1./2.*armHHos[145] + 
   armHHos[146];
   armHHos[145]=armHHos[56]*armHHos[145];
   armHHos[89]=7 + armHHos[89];
   armHHos[89]=1./2.*armHHos[89];
   armHHos[108]=armHHos[89] + armHHos[108];
   armHHos[108]= - 15./2.*armHHos[13] + 1./2.*armHHos[108] + 5*
   armHHos[11];
   armHHos[108]=armHHos[1]*armHHos[108];
   armHHos[146]=armHHos[9]*armHHos[171];
   armHHos[108]=17./2.*armHHos[146] + 11./3.*armHHos[145] + 
   armHHos[108];
   armHHos[108]=mmZ*armHHos[108];
   armHHos[108]=1./2.*armHHos[165] + armHHos[108];
   armHHos[108]=armHHos[138] + 1./2.*armHHos[108] + armHHos[116];
   armHHos[108]=armHHos[4]*armHHos[108];
   armHHos[89]=armHHos[93] + armHHos[92] + armHHos[89] + armHHos[91];
   armHHos[92]=armHHos[1]*armHHos[89];
   armHHos[116]=armHHos[169] - armHHos[18];
   armHHos[138]=1./2.*armHHos[11];
   armHHos[145]=1./3.*armHHos[13] + 1./3.*armHHos[116] + armHHos[138];
   armHHos[145]=armHHos[56]*armHHos[145];
   armHHos[165]=11./9.*armHHos[56] + armHHos[1];
   armHHos[165]=armHHos[9]*armHHos[165];
   armHHos[92]=armHHos[165] + 11*armHHos[145] + armHHos[92];
   armHHos[92]=mmZ*armHHos[87]*armHHos[92];
   armHHos[145]=armHHos[102] + armHHos[17];
   armHHos[169]=armHHos[56]*armHHos[145];
   armHHos[145]=armHHos[1]*armHHos[145];
   armHHos[171]=1./2.*armHHos[173] + 11./3.*armHHos[169] + 3*
   armHHos[145];
   armHHos[171]=armHHos[87]*armHHos[171];
   armHHos[91]=armHHos[93] - 7./2.*armHHos[11] + 13./8. + armHHos[91];
   armHHos[91]=armHHos[1]*armHHos[91];
   armHHos[93]=armHHos[13] + armHHos[144] + 13./24. - armHHos[18];
   armHHos[93]=armHHos[56]*armHHos[93];
   armHHos[91]=armHHos[92] + armHHos[171] + 5*armHHos[165] + 11./3.*
   armHHos[93] + armHHos[91];
   armHHos[92]=armHHos[100]*armHHos[57]*armHHos[155];
   armHHos[91]=1./8.*armHHos[92] + armHHos[108] + 1./4.*armHHos[91] + 
   armHHos[113];
   armHHos[91]=armHHos[100]*armHHos[91];
   armHHos[92]=19*armHHos[29];
   armHHos[93]=armHHos[149] + armHHos[92] + 101*armHHos[14];
   armHHos[93]=armHHos[9]*armHHos[93];
   armHHos[108]= - 23 + armHHos[125];
   armHHos[108]=armHHos[14]*armHHos[108];
   armHHos[113]=137./3. - 9*armHHos[7];
   armHHos[113]=armHHos[15]*armHHos[113];
   armHHos[108]=1./2.*armHHos[93] + 1./2.*armHHos[113] + 1./2.*
   armHHos[108] + armHHos[163] - 119./2.*armHHos[17] + 19./3.*
   armHHos[16] + armHHos[161] + armHHos[158] + armHHos[157] - 91./3.*
   armHHos[84] + 149./4.*armHHos[58];
   armHHos[113]=1./2.*armHHos[77] - 5 + 49./2.*armHHos[75];
   armHHos[144]= - 17./4.*armHHos[9];
   armHHos[155]=armHHos[144] - 23 - 17./4.*armHHos[7];
   armHHos[155]=armHHos[9]*armHHos[155];
   armHHos[157]=armHHos[140] - 3*armHHos[65];
   armHHos[157]=mmZ*armHHos[157];
   armHHos[158]= - 8./3.*armHHos[78];
   armHHos[161]=3*armHHos[74];
   armHHos[113]=armHHos[157] + 1./4.*armHHos[155] + 5./8.*armHHos[7] - 
   73./8.*armHHos[13] - 1./8.*armHHos[12] - 65./16.*armHHos[79] + 7./16.
   *armHHos[80] + armHHos[161] + 1./4.*armHHos[113] + armHHos[158];
   armHHos[113]=mmZ*armHHos[113];
   armHHos[108]=1./4.*armHHos[108] + armHHos[113];
   armHHos[108]=armHHos[57]*armHHos[108];
   armHHos[113]=6*armHHos[78];
   armHHos[155]=9./4.*armHHos[166] + 5./2.*armHHos[79] + 1./4.*
   armHHos[80] + 41./4. + armHHos[113];
   armHHos[155]=mmZ*armHHos[155];
   armHHos[157]=armHHos[164] - armHHos[29] + 3./4.*armHHos[14];
   armHHos[157]=armHHos[9]*armHHos[157];
   armHHos[155]=armHHos[155] + 3*armHHos[157] - 25./2.*armHHos[15] - 11.
   /2.*armHHos[14] - 3*armHHos[29] + armHHos[185] - 3*armHHos[60] + 5*
   armHHos[104] + 6*armHHos[82];
   armHHos[155]=mmZ*armHHos[155];
   armHHos[163]=armHHos[14] + 1./2.*armHHos[15];
   armHHos[163]=armHHos[15]*armHHos[163];
   armHHos[155]=9*armHHos[163] + armHHos[155];
   armHHos[155]=armHHos[4]*armHHos[57]*armHHos[155];
   armHHos[108]=armHHos[108] + armHHos[155];
   armHHos[108]=armHHos[4]*armHHos[108];
   armHHos[91]=armHHos[91] + armHHos[95] + armHHos[108];
   armHHos[91]=armHHos[100]*armHHos[91];
   armHHos[95]=armHHos[99] + armHHos[7] - 5./4. + armHHos[98];
   armHHos[95]=armHHos[9]*armHHos[95];
   armHHos[108]=17./2. - 11*armHHos[75];
   armHHos[108]=1./4.*armHHos[108] + 9*armHHos[72];
   armHHos[155]= - 27./4.*armHHos[28];
   armHHos[165]=19./8.*armHHos[13];
   armHHos[95]=1./4.*armHHos[95] + armHHos[151] + armHHos[155] + 
   armHHos[165] + 33./16.*armHHos[79] - 9./2.*armHHos[30] + 55./16.*
   armHHos[80] - armHHos[74] + 1./2.*armHHos[108] + armHHos[113];
   armHHos[95]=armHHos[87]*armHHos[95];
   armHHos[108]=1091./24. - 27*armHHos[73];
   armHHos[113]= - 57./8.*armHHos[7] + 55./3. + armHHos[155];
   armHHos[113]=armHHos[7]*armHHos[113];
   armHHos[151]=3*armHHos[76];
   armHHos[108]=armHHos[113] + armHHos[156] + 93./8.*armHHos[12] + 33./
   8.*armHHos[79] + armHHos[199] + 165./4.*armHHos[80] + armHHos[151]
    - 165./8.*armHHos[77] + 1./2.*armHHos[108] - 18*armHHos[81];
   armHHos[108]=armHHos[88]*armHHos[108];
   armHHos[113]=armHHos[9]*armHHos[189];
   armHHos[95]=3*armHHos[95] + 3./16.*armHHos[113] + armHHos[108] + 
   armHHos[118] + 67./4.*armHHos[68] + 89./4.*armHHos[66] - 2*
   armHHos[67] + armHHos[139] - 10*armHHos[69] - 8*armHHos[70] + 9*
   armHHos[64];
   armHHos[95]=mmZ*armHHos[95];
   armHHos[108]=11./4.*armHHos[104] + armHHos[82];
   armHHos[108]= - 31./2.*armHHos[14] - armHHos[29] + 167./4.*
   armHHos[17] + armHHos[126] + 3*armHHos[108] + armHHos[124];
   armHHos[113]=armHHos[127] - 10 + armHHos[121];
   armHHos[113]=armHHos[15]*armHHos[113];
   armHHos[108]=armHHos[132] + 1./2.*armHHos[108] + 3*armHHos[113];
   armHHos[108]=armHHos[87]*armHHos[108];
   armHHos[113]=737./2. + 27*armHHos[75];
   armHHos[113]=armHHos[77] + 1./2.*armHHos[113] + armHHos[97];
   armHHos[113]=1./2.*armHHos[113] + armHHos[134];
   armHHos[118]=armHHos[137] + armHHos[136] + 9./8. - armHHos[7];
   armHHos[118]=armHHos[9]*armHHos[118];
   armHHos[121]=2*armHHos[70];
   armHHos[127]=1./2.*armHHos[65];
   armHHos[132]=armHHos[127] - 17./8.*armHHos[68] + 3./4.*armHHos[66]
    + 9./2.*armHHos[63] + armHHos[121] + 3*armHHos[69];
   armHHos[132]=mmH*armHHos[132];
   armHHos[134]=65 + 43./3.*armHHos[7];
   armHHos[134]=armHHos[7]*armHHos[134];
   armHHos[136]=11*armHHos[104] + 3*armHHos[16];
   armHHos[136]=armHHos[88]*armHHos[136];
   armHHos[137]= - 59 - 95*armHHos[7];
   armHHos[137]=armHHos[14]*armHHos[88]*armHHos[137];
   armHHos[139]= - 41 - 25*armHHos[7];
   armHHos[139]=armHHos[15]*armHHos[88]*armHHos[139];
   armHHos[95]=armHHos[95] + armHHos[132] + armHHos[108] + armHHos[118]
    + 1./8.*armHHos[139] + 1./8.*armHHos[137] + 3./4.*armHHos[136] + 1./
   4.*armHHos[134] + armHHos[152] - 19./8.*armHHos[13] - 9./4.*
   armHHos[12] + 89./8.*armHHos[79] + armHHos[154] + 251./8.*
   armHHos[80] + 1./2.*armHHos[113] - armHHos[74];
   armHHos[95]=armHHos[57]*armHHos[95];
   armHHos[108]= - 157./3. + 27*armHHos[73];
   armHHos[113]= - 3*armHHos[76];
   armHHos[118]=9./4.*armHHos[7] + 31./3. + 27./2.*armHHos[28];
   armHHos[118]=armHHos[7]*armHHos[118];
   armHHos[132]= - 99./4.*armHHos[80];
   armHHos[108]=1./2.*armHHos[118] + armHHos[114] - 75./8.*armHHos[12]
    - 99./16.*armHHos[79] + armHHos[112] + armHHos[132] + armHHos[113]
    + 99./8.*armHHos[77] + 1./2.*armHHos[108] + 18*armHHos[81];
   armHHos[108]=armHHos[88]*armHHos[108];
   armHHos[97]=armHHos[97] - 157./3. + 99./4.*armHHos[75];
   armHHos[112]=9./4.*armHHos[9] + 9./2.*armHHos[7] + 275./12. + 
   armHHos[128];
   armHHos[112]=armHHos[9]*armHHos[112];
   armHHos[97]=1./8.*armHHos[112] - 9./16.*armHHos[7] + armHHos[172] - 
   75./16.*armHHos[13] - 297./32.*armHHos[79] + armHHos[154] - 99./16.*
   armHHos[80] + armHHos[143] + 1./4.*armHHos[97] + armHHos[141];
   armHHos[97]=armHHos[87]*armHHos[97];
   armHHos[112]= - 1./2.*armHHos[63];
   armHHos[114]= - armHHos[64] + armHHos[112];
   armHHos[97]=armHHos[97] + armHHos[120] + armHHos[108] + armHHos[65]
    - 39./4.*armHHos[68] - 39./2.*armHHos[66] + 9*armHHos[114] + 2*
   armHHos[67];
   armHHos[97]=mmZ*armHHos[97];
   armHHos[108]= - 99./4.*armHHos[84];
   armHHos[114]=99./8.*armHHos[58];
   armHHos[118]= - armHHos[7]*armHHos[29];
   armHHos[120]=7./4.*armHHos[118] - armHHos[29] + 35./4.*armHHos[16]
    + armHHos[126] + armHHos[114] + armHHos[108] + 3*armHHos[83] - 7./4.
   *armHHos[59];
   armHHos[120]=armHHos[88]*armHHos[120];
   armHHos[134]= - armHHos[84] + 1./2.*armHHos[58];
   armHHos[134]=33./4.*armHHos[134] + armHHos[82];
   armHHos[124]=81./2.*armHHos[14] - armHHos[29] + 35./4.*armHHos[17]
    + armHHos[126] + 3*armHHos[134] + armHHos[124];
   armHHos[126]=9*armHHos[14];
   armHHos[129]=armHHos[131] + armHHos[129] + armHHos[126];
   armHHos[129]=armHHos[9]*armHHos[129];
   armHHos[131]=9./8.*armHHos[7] + 29./3. + armHHos[156];
   armHHos[131]=armHHos[15]*armHHos[131];
   armHHos[124]=1./2.*armHHos[129] + 1./2.*armHHos[124] + armHHos[131];
   armHHos[124]=armHHos[87]*armHHos[124];
   armHHos[129]=5./4.*armHHos[9] + 9./4.*armHHos[147] - 331./12. + 
   armHHos[117];
   armHHos[129]=armHHos[9]*armHHos[129];
   armHHos[131]=armHHos[64] + 1./2.*armHHos[63];
   armHHos[127]=armHHos[127] + 11./8.*armHHos[68] + 11./4.*armHHos[66]
    + 9*armHHos[131] + armHHos[67];
   armHHos[127]=mmH*armHHos[127];
   armHHos[131]=3*armHHos[73];
   armHHos[134]= - 253./8. + armHHos[131];
   armHHos[134]= - 29./4.*armHHos[77] + 27./4.*armHHos[72] - 29./8.*
   armHHos[75] + 9./2.*armHHos[134] + 13./3.*armHHos[81];
   armHHos[117]= - 331./6. + armHHos[117];
   armHHos[117]=armHHos[7]*armHHos[117];
   armHHos[115]=475./6. + armHHos[115];
   armHHos[115]=1./2.*armHHos[115] + 11*armHHos[7];
   armHHos[115]=armHHos[14]*armHHos[88]*armHHos[115];
   armHHos[136]=9./2. + armHHos[7];
   armHHos[136]=armHHos[15]*armHHos[88]*armHHos[136];
   armHHos[137]= - 1./2.*armHHos[74];
   armHHos[97]=armHHos[97] + 1./2.*armHHos[127] + 1./2.*armHHos[124] + 
   1./4.*armHHos[129] + 9./4.*armHHos[136] + 1./2.*armHHos[115] + 1./2.
   *armHHos[120] + 1./4.*armHHos[117] - 297./16.*armHHos[28] + 37./16.*
   armHHos[13] + 37./8.*armHHos[12] - 85./4.*armHHos[79] - 81./8.*
   armHHos[30] - 85./2.*armHHos[80] + armHHos[137] + 13./12.*
   armHHos[78] + 1./2.*armHHos[134] - armHHos[76];
   armHHos[97]=armHHos[57]*armHHos[97];
   armHHos[115]=3./2. - armHHos[7];
   armHHos[115]=armHHos[7]*armHHos[115];
   armHHos[117]= - 1./4.*armHHos[9] + 3./4. - armHHos[7];
   armHHos[117]=armHHos[9]*armHHos[117];
   armHHos[115]=3./2.*armHHos[117] + 3./2.*armHHos[115] - 99./8.*
   armHHos[79] + armHHos[132] + armHHos[78] - 255./8. + 2*armHHos[81];
   armHHos[115]=mmZ*armHHos[115];
   armHHos[109]=armHHos[164] - armHHos[29] + armHHos[109];
   armHHos[109]=armHHos[9]*armHHos[109];
   armHHos[117]=50 + armHHos[203];
   armHHos[117]=armHHos[14]*armHHos[117];
   armHHos[120]=25 + 3./4.*armHHos[7];
   armHHos[120]=armHHos[15]*armHHos[120];
   armHHos[124]= - 1./2.*armHHos[60];
   armHHos[108]=armHHos[115] + 1./2.*armHHos[109] + armHHos[120] + 
   armHHos[117] + armHHos[118] + armHHos[187] + armHHos[170] + 
   armHHos[124] + armHHos[82] + armHHos[114] + armHHos[108] + 2*
   armHHos[83] - armHHos[59];
   armHHos[108]=mmZ*armHHos[108];
   armHHos[109]=armHHos[14] + 1./4.*armHHos[15];
   armHHos[109]=armHHos[15]*armHHos[109];
   armHHos[114]=pow(armHHos[14],2);
   armHHos[109]=armHHos[114] + armHHos[109];
   armHHos[108]=3*armHHos[109] + armHHos[108];
   armHHos[108]=armHHos[57]*armHHos[108];
   armHHos[109]=armHHos[160] + armHHos[159] - armHHos[10] - 3./2. - 
   armHHos[21];
   armHHos[115]=armHHos[56]*armHHos[109];
   armHHos[109]=armHHos[1]*armHHos[109];
   armHHos[109]=3*armHHos[115] + armHHos[109];
   armHHos[109]=mmZ*armHHos[109];
   armHHos[102]= - armHHos[5] + armHHos[102];
   armHHos[115]=armHHos[56]*armHHos[102];
   armHHos[102]=armHHos[1]*armHHos[102];
   armHHos[117]=3*armHHos[56] + armHHos[1];
   armHHos[120]=armHHos[14]*armHHos[117];
   armHHos[117]=armHHos[15]*armHHos[117];
   armHHos[102]=armHHos[109] + 1./2.*armHHos[117] + armHHos[120] + 3*
   armHHos[115] + armHHos[102];
   armHHos[102]=mmZ*armHHos[102];
   armHHos[102]=armHHos[102] + 3*armHHos[108];
   armHHos[102]=armHHos[4]*armHHos[102];
   armHHos[108]= - 181./6. + armHHos[153];
   armHHos[108]=armHHos[7]*armHHos[108];
   armHHos[109]=armHHos[144] - 181./12. + 51*armHHos[7];
   armHHos[109]=armHHos[9]*armHHos[109];
   armHHos[115]= - 1./2.*armHHos[65];
   armHHos[117]=armHHos[115] + 33./4.*armHHos[68] - armHHos[67] + 33./2.
   *armHHos[66];
   armHHos[117]=mmZ*armHHos[117];
   armHHos[108]=3*armHHos[117] + 1./8.*armHHos[109] + 1./8.*
   armHHos[108] + 207./16.*armHHos[13] + 207./8.*armHHos[12] + 1635./32.
   *armHHos[79] + 1635./16.*armHHos[80] + armHHos[168] + armHHos[167]
    + armHHos[151] - 231./8.*armHHos[77] - 231./16.*armHHos[75] + 1341./
   16. - 8./3.*armHHos[81];
   armHHos[108]=mmZ*armHHos[108];
   armHHos[92]=armHHos[149] + armHHos[92] + 231*armHHos[14];
   armHHos[92]=armHHos[9]*armHHos[92];
   armHHos[109]=armHHos[7]*armHHos[29];
   armHHos[117]=19./2.*armHHos[109] - 43./2.*armHHos[29] + 2059./12.*
   armHHos[17] + 2059./6.*armHHos[16] - 77./4.*armHHos[32] + 37./4.*
   armHHos[60] - 17./6.*armHHos[82] - 1575./8.*armHHos[58] + 423./2.*
   armHHos[84] - 17./3.*armHHos[83] + 37./2.*armHHos[59];
   armHHos[120]=13*armHHos[7];
   armHHos[127]= - 947./6. + armHHos[120];
   armHHos[127]=armHHos[14]*armHHos[127];
   armHHos[129]= - 947./3. + 231./4.*armHHos[7];
   armHHos[129]=armHHos[15]*armHHos[129];
   armHHos[92]=armHHos[108] + 1./16.*armHHos[92] + 1./4.*armHHos[129]
    + 1./4.*armHHos[117] + armHHos[127];
   armHHos[92]=armHHos[57]*armHHos[92];
   armHHos[108]=armHHos[5] + 1./2.*armHHos[6];
   armHHos[94]=armHHos[94] + 11./2.*armHHos[108] - 5*armHHos[16];
   armHHos[108]=armHHos[56]*armHHos[94];
   armHHos[94]=armHHos[1]*armHHos[94];
   armHHos[117]= - armHHos[56] - 1./3.*armHHos[1];
   armHHos[127]=armHHos[14]*armHHos[117];
   armHHos[129]=armHHos[15]*armHHos[117];
   armHHos[94]=13./2.*armHHos[129] + 13*armHHos[127] + 3*armHHos[108]
    + armHHos[94];
   armHHos[108]=1./2.*armHHos[21];
   armHHos[132]=armHHos[108] + 7./4. + 5*armHHos[19];
   armHHos[110]=armHHos[150] - 15./2.*armHHos[12] + armHHos[110] + 15./
   4.*armHHos[18] + 3./8.*armHHos[20] + 3./2.*armHHos[132] + 5*
   armHHos[10];
   armHHos[110]=armHHos[56]*armHHos[110];
   armHHos[132]= - 5./4.*armHHos[13] - 5./2.*armHHos[12] + 5./6.*
   armHHos[11] + 5./4.*armHHos[18] + 1./8.*armHHos[20] + 1./2.*
   armHHos[132] + 5./3.*armHHos[10];
   armHHos[132]=armHHos[1]*armHHos[132];
   armHHos[134]=armHHos[7]*armHHos[117];
   armHHos[136]=armHHos[9]*armHHos[117];
   armHHos[110]=17./4.*armHHos[136] + 17./2.*armHHos[134] + 
   armHHos[110] + armHHos[132];
   armHHos[110]=mmZ*armHHos[110];
   armHHos[92]=armHHos[102] + armHHos[92] + 1./2.*armHHos[94] + 
   armHHos[110];
   armHHos[92]=armHHos[4]*armHHos[92];
   armHHos[89]=armHHos[56]*armHHos[89];
   armHHos[94]=armHHos[13] + armHHos[116] + armHHos[178];
   armHHos[94]=armHHos[1]*armHHos[94];
   armHHos[102]=armHHos[56] + 1./3.*armHHos[1];
   armHHos[110]=armHHos[9]*armHHos[102];
   armHHos[89]=armHHos[110] + armHHos[89] + armHHos[94];
   armHHos[89]=armHHos[87]*armHHos[89];
   armHHos[94]=armHHos[96] + 9./2.*armHHos[10] + 3./2.*armHHos[21] + 7./
   2. - 3*armHHos[19];
   armHHos[94]=armHHos[56]*armHHos[94];
   armHHos[108]=armHHos[12] + 3./2.*armHHos[10] + armHHos[108] + 7./6.
    - armHHos[19];
   armHHos[108]=armHHos[1]*armHHos[108];
   armHHos[102]=armHHos[7]*armHHos[102];
   armHHos[94]=armHHos[102] + armHHos[94] + armHHos[108];
   armHHos[94]=armHHos[88]*armHHos[94];
   armHHos[89]=armHHos[94] + 1./2.*armHHos[89];
   armHHos[89]=mmZ*armHHos[89];
   armHHos[94]=13./16. - armHHos[19];
   armHHos[107]=3./2.*armHHos[13] + armHHos[96] - 7./4.*armHHos[11] + 
   armHHos[107] + 3*armHHos[94] - 7./2.*armHHos[10];
   armHHos[107]=armHHos[56]*armHHos[107];
   armHHos[108]=1./2.*armHHos[129];
   armHHos[116]=armHHos[108] + 3*armHHos[169] + armHHos[145];
   armHHos[116]=armHHos[87]*armHHos[116];
   armHHos[129]=1./2.*armHHos[13];
   armHHos[94]=armHHos[129] + armHHos[12] - 7./12.*armHHos[11] - 1./2.*
   armHHos[18] + armHHos[94] - 7./6.*armHHos[10];
   armHHos[94]=armHHos[1]*armHHos[94];
   armHHos[132]= - 1./2.*armHHos[5] + armHHos[16];
   armHHos[139]=armHHos[56]*armHHos[132];
   armHHos[132]=armHHos[1]*armHHos[132];
   armHHos[132]=3*armHHos[139] + armHHos[132];
   armHHos[132]=armHHos[88]*armHHos[132];
   armHHos[117]=armHHos[14]*armHHos[88]*armHHos[117];
   armHHos[89]=armHHos[89] + 1./2.*armHHos[116] + 5./2.*armHHos[110] + 
   1./2.*armHHos[117] + armHHos[132] + 5*armHHos[102] + armHHos[107] + 
   armHHos[94];
   armHHos[89]=armHHos[92] + 1./2.*armHHos[89] + armHHos[97];
   armHHos[89]=armHHos[103]*armHHos[89];
   armHHos[92]= - 353 + 8*armHHos[81];
   armHHos[94]= - 215 - 253./6.*armHHos[7];
   armHHos[94]=armHHos[7]*armHHos[94];
   armHHos[97]=armHHos[144] + 13 - 209./4.*armHHos[7];
   armHHos[97]=armHHos[9]*armHHos[97];
   armHHos[107]=armHHos[121] + armHHos[69];
   armHHos[107]=2*armHHos[107] + armHHos[67];
   armHHos[107]= - armHHos[65] - 22*armHHos[68] + 2*armHHos[107] - 55./
   2.*armHHos[66];
   armHHos[107]=mmZ*armHHos[107];
   armHHos[92]=3*armHHos[107] + 1./4.*armHHos[97] + 1./4.*armHHos[94]
    - 145./8.*armHHos[13] - 269./8.*armHHos[12] - 953./16.*armHHos[79]
    - 2317./16.*armHHos[80] + armHHos[161] + armHHos[158] + 
   armHHos[113] + 341./8.*armHHos[77] + 1./3.*armHHos[92] + 121./8.*
   armHHos[75];
   armHHos[92]=mmZ*armHHos[92];
   armHHos[94]=25*armHHos[7];
   armHHos[97]=587 + armHHos[94];
   armHHos[97]=armHHos[14]*armHHos[97];
   armHHos[107]=715 - 289*armHHos[7];
   armHHos[107]=armHHos[15]*armHHos[107];
   armHHos[92]=armHHos[92] + 1./8.*armHHos[93] + 1./8.*armHHos[107] + 1.
   /24.*armHHos[97] - 43./12.*armHHos[29] - 1589./24.*armHHos[17] - 445.
   /12.*armHHos[16] - 77./24.*armHHos[32] + 37./8.*armHHos[60] - 17./12.
   *armHHos[82] + 525./16.*armHHos[58] + 8./3.*armHHos[85] - 141./4.*
   armHHos[84];
   armHHos[92]=armHHos[57]*armHHos[92];
   armHHos[93]=armHHos[10] + 1 + armHHos[21];
   armHHos[97]=armHHos[56]*armHHos[93];
   armHHos[93]=armHHos[1]*armHHos[93];
   armHHos[93]=3*armHHos[97] + armHHos[93];
   armHHos[93]=mmZ*armHHos[93];
   armHHos[97]=armHHos[56]*armHHos[5];
   armHHos[107]=armHHos[1]*armHHos[5];
   armHHos[110]= - 3*armHHos[56] - armHHos[1];
   armHHos[110]=armHHos[14]*armHHos[110];
   armHHos[93]=2*armHHos[93] + armHHos[110] + 3*armHHos[97] + 
   armHHos[107];
   armHHos[93]=mmZ*armHHos[93];
   armHHos[97]=armHHos[7]*armHHos[142];
   armHHos[97]=3./4.*armHHos[166] + 3*armHHos[97] + 33./2.*armHHos[79]
    + 231./4.*armHHos[80] + 2*armHHos[78] + 283./4. - 4*armHHos[81];
   armHHos[97]=mmZ*armHHos[97];
   armHHos[107]= - 127 + armHHos[130];
   armHHos[107]=armHHos[14]*armHHos[107];
   armHHos[97]=armHHos[97] + armHHos[157] - 71./2.*armHHos[15] + 1./2.*
   armHHos[107] + armHHos[109] - armHHos[60] + 2*armHHos[82] - 33./2.*
   armHHos[58] + 33*armHHos[84] - 2*armHHos[83] + armHHos[59];
   armHHos[97]=mmZ*armHHos[97];
   armHHos[97]=3*armHHos[163] + armHHos[97];
   armHHos[97]=armHHos[57]*armHHos[97];
   armHHos[93]=armHHos[93] + 3*armHHos[97];
   armHHos[93]=armHHos[4]*armHHos[93];
   armHHos[97]= - 3./2.*armHHos[21];
   armHHos[107]=armHHos[97] - 7./2. - 15*armHHos[19];
   armHHos[107]=15./2.*armHHos[12] + 1./2.*armHHos[107] - 5*armHHos[10]
   ;
   armHHos[107]=armHHos[56]*armHHos[107];
   armHHos[110]= - 1./2.*armHHos[21];
   armHHos[113]=armHHos[110] - 7./6. - 5*armHHos[19];
   armHHos[113]=5./2.*armHHos[12] + 1./2.*armHHos[113] - 5./3.*
   armHHos[10];
   armHHos[113]=armHHos[1]*armHHos[113];
   armHHos[102]=17./2.*armHHos[102] + armHHos[107] + armHHos[113];
   armHHos[102]=mmZ*armHHos[102];
   armHHos[92]=armHHos[93] + armHHos[102] + armHHos[92];
   armHHos[92]=armHHos[4]*armHHos[92];
   armHHos[93]=armHHos[181] - 9./2.*armHHos[10] + armHHos[97] - 7./2.
    + 3*armHHos[19];
   armHHos[93]=armHHos[56]*armHHos[93];
   armHHos[97]= - armHHos[12] - 3./2.*armHHos[10] + armHHos[110] - 7./6.
    + armHHos[19];
   armHHos[97]=armHHos[1]*armHHos[97];
   armHHos[93]=armHHos[134] + armHHos[93] + armHHos[97];
   armHHos[93]=mmZ*armHHos[88]*armHHos[93];
   armHHos[89]=armHHos[89] + armHHos[92] + 1./2.*armHHos[93] + 
   armHHos[95];
   armHHos[89]=armHHos[103]*armHHos[89];
   armHHos[92]= - 1./2.*armHHos[28];
   armHHos[93]=1 + armHHos[92];
   armHHos[95]=3./2.*armHHos[9];
   armHHos[93]=armHHos[95] + 3*armHHos[93] + armHHos[188];
   armHHos[93]=armHHos[9]*armHHos[93];
   armHHos[97]=15./4.*armHHos[77] - 21./4.*armHHos[72] - 3./4.*
   armHHos[75] + 1./3. - 21./4.*armHHos[73];
   armHHos[102]= - 2*armHHos[74];
   armHHos[107]=21./4.*armHHos[30];
   armHHos[110]= - 3./8.*armHHos[79];
   armHHos[113]=1./8.*armHHos[12];
   armHHos[116]=armHHos[183] + 4./3. - 3./16.*armHHos[28];
   armHHos[116]=armHHos[7]*armHHos[116];
   armHHos[117]= - armHHos[65] + 1./2.*armHHos[68] + armHHos[140] - 
   armHHos[69] - armHHos[67];
   armHHos[117]=1./2.*mmH*armHHos[117];
   armHHos[121]= - 2*armHHos[76];
   armHHos[93]=armHHos[117] + 1./8.*armHHos[93] + armHHos[116] + 
   armHHos[148] + armHHos[165] + armHHos[113] + armHHos[110] + 
   armHHos[107] + 15./8.*armHHos[80] + armHHos[102] + 1./2.*armHHos[97]
    + armHHos[121];
   armHHos[93]=mmH*armHHos[93];
   armHHos[97]=637./3. + armHHos[128];
   armHHos[130]= - 61./3.*armHHos[7];
   armHHos[132]=23./2.*armHHos[119];
   armHHos[139]= - armHHos[88]*armHHos[29];
   armHHos[97]=armHHos[132] + armHHos[139] + 1./2.*armHHos[97] + 
   armHHos[130];
   armHHos[97]=armHHos[14]*armHHos[97];
   armHHos[141]=7*armHHos[119];
   armHHos[94]=armHHos[141] + armHHos[94] + 7 + armHHos[128];
   armHHos[94]=armHHos[15]*armHHos[94];
   armHHos[143]= - 3*armHHos[14];
   armHHos[144]= - armHHos[29] + armHHos[143];
   armHHos[145]=3*armHHos[15];
   armHHos[144]=11./2.*armHHos[144] + armHHos[145];
   armHHos[144]=1./4.*armHHos[9]*armHHos[144];
   armHHos[149]= - 1./2.*armHHos[29];
   armHHos[150]=armHHos[162] + armHHos[149] - armHHos[14];
   armHHos[150]=1./2.*armHHos[87]*armHHos[15]*armHHos[150];
   armHHos[151]= - 17./8.*armHHos[59];
   armHHos[152]=armHHos[151] - 4./3.*armHHos[85] + armHHos[83];
   armHHos[153]= - 17./8.*armHHos[60];
   armHHos[154]=9./4.*armHHos[32];
   armHHos[155]=13./6.*armHHos[29];
   armHHos[156]=11./8.*armHHos[118];
   armHHos[93]=armHHos[93] + armHHos[150] + armHHos[144] + 1./8.*
   armHHos[94] + 1./4.*armHHos[97] + armHHos[156] + armHHos[155] + 241./
   24.*armHHos[17] - 73./24.*armHHos[16] + armHHos[154] + armHHos[153]
    + armHHos[82] + 3./2.*armHHos[58] + armHHos[152] - 6*armHHos[84];
   armHHos[93]=armHHos[57]*armHHos[93];
   armHHos[94]=1./2.*armHHos[147] - 67./3. + armHHos[188];
   armHHos[94]=armHHos[14]*armHHos[94];
   armHHos[97]= - 29 + armHHos[98];
   armHHos[97]=1./4.*armHHos[119] + 3./8.*armHHos[97] + armHHos[7];
   armHHos[97]=armHHos[15]*armHHos[97];
   armHHos[119]=armHHos[162] - 11./4.*armHHos[29] + armHHos[143];
   armHHos[119]=armHHos[9]*armHHos[119];
   armHHos[157]=armHHos[162] + armHHos[149] + armHHos[14];
   armHHos[157]=armHHos[87]*armHHos[15]*armHHos[157];
   armHHos[158]=61./3.*armHHos[84] - 13*armHHos[58];
   armHHos[94]=1./2.*armHHos[157] + 1./2.*armHHos[119] + armHHos[97] + 
   1./2.*armHHos[94] + 13./12.*armHHos[29] + 17./3.*armHHos[17] + 11./4.
   *armHHos[16] + 9./8.*armHHos[32] + armHHos[153] + 1./4.*armHHos[158]
    + armHHos[82];
   armHHos[97]= - 3*armHHos[28];
   armHHos[119]=29./3. + armHHos[97];
   armHHos[119]=armHHos[95] + 1./2.*armHHos[119] - armHHos[7];
   armHHos[119]=armHHos[9]*armHHos[119];
   armHHos[157]= - 37./12. - 21*armHHos[72];
   armHHos[157]=1./2.*armHHos[157] - 15*armHHos[77];
   armHHos[158]=3 + armHHos[7];
   armHHos[158]=armHHos[7]*armHHos[158];
   armHHos[159]= - 1./2.*armHHos[68] - armHHos[65];
   armHHos[159]=mmH*armHHos[159];
   armHHos[119]=1./4.*armHHos[159] + 1./16.*armHHos[119] + 1./16.*
   armHHos[158] + 27./32.*armHHos[28] + armHHos[13] + 15./8.*
   armHHos[12] - 3./2.*armHHos[79] + 21./16.*armHHos[30] + 1./8.*
   armHHos[80] + 1./8.*armHHos[157] - armHHos[74];
   armHHos[119]=mmH*armHHos[119];
   armHHos[94]=1./2.*armHHos[94] + armHHos[119];
   armHHos[94]=armHHos[57]*armHHos[94];
   armHHos[119]=armHHos[56]*armHHos[11];
   armHHos[157]=armHHos[1]*armHHos[11];
   armHHos[119]=armHHos[146] + 11./9.*armHHos[119] + armHHos[157];
   armHHos[119]=mmH*armHHos[119];
   armHHos[119]=armHHos[173] + armHHos[119];
   armHHos[146]=7./4.*armHHos[15] + 5*armHHos[29] + 19./8.*armHHos[14];
   armHHos[146]=armHHos[15]*armHHos[146];
   armHHos[146]= - 1./8.*armHHos[114] + armHHos[146];
   armHHos[146]=armHHos[4]*armHHos[57]*armHHos[146];
   armHHos[157]=armHHos[12] - armHHos[79] - 1./4. - armHHos[77];
   armHHos[157]=mmH*armHHos[157];
   armHHos[104]=1./2.*armHHos[157] - armHHos[15] - 5./2.*armHHos[14] + 
   armHHos[104] + 1./2.*armHHos[16];
   armHHos[104]=armHHos[100]*armHHos[57]*armHHos[104];
   armHHos[94]=1./4.*armHHos[104] + 1./3.*armHHos[146] + 1./4.*
   armHHos[119] + armHHos[94];
   armHHos[94]=armHHos[100]*armHHos[94];
   armHHos[104]=armHHos[56]*armHHos[10];
   armHHos[119]=armHHos[1]*armHHos[10];
   armHHos[104]=armHHos[134] + armHHos[104] + 1./3.*armHHos[119];
   armHHos[104]=mmH*armHHos[104];
   armHHos[104]=armHHos[127] + armHHos[104];
   armHHos[104]=1./2.*armHHos[104];
   armHHos[119]=10*armHHos[29];
   armHHos[146]=armHHos[119] + 373./8.*armHHos[14];
   armHHos[146]=armHHos[14]*armHHos[146];
   armHHos[157]=7./6.*armHHos[15] + 10./3.*armHHos[29] + 53./8.*
   armHHos[14];
   armHHos[157]=armHHos[15]*armHHos[157];
   armHHos[146]=1./3.*armHHos[146] + armHHos[157];
   armHHos[146]=armHHos[4]*armHHos[57]*armHHos[146];
   armHHos[93]=armHHos[94] + armHHos[146] + armHHos[104] + armHHos[93];
   armHHos[93]=armHHos[100]*armHHos[93];
   armHHos[92]= - 3 + armHHos[92];
   armHHos[94]=armHHos[95] + 3*armHHos[92] + armHHos[188];
   armHHos[94]=armHHos[9]*armHHos[94];
   armHHos[95]= - 21./2.*armHHos[72];
   armHHos[157]=15./2.*armHHos[77] + armHHos[95] - 15./2.*armHHos[75]
    - 25./3. - 21./2.*armHHos[73];
   armHHos[94]=armHHos[117] + 1./8.*armHHos[94] + armHHos[116] + 
   armHHos[148] + 31./8.*armHHos[13] + armHHos[113] + armHHos[110] + 
   armHHos[107] + 3./8.*armHHos[80] + armHHos[102] + 1./4.*armHHos[157]
    + armHHos[121];
   armHHos[94]=mmH*armHHos[94];
   armHHos[102]=277./3. + armHHos[128];
   armHHos[102]=armHHos[132] + armHHos[139] + 1./2.*armHHos[102] + 
   armHHos[130];
   armHHos[102]=armHHos[14]*armHHos[102];
   armHHos[107]=armHHos[141] + 57*armHHos[7] - 109 + armHHos[128];
   armHHos[107]=armHHos[15]*armHHos[107];
   armHHos[94]=armHHos[94] + armHHos[150] + armHHos[144] + 1./8.*
   armHHos[107] + 1./4.*armHHos[102] + armHHos[156] + armHHos[155] + 
   373./24.*armHHos[17] - 145./24.*armHHos[16] + armHHos[154] + 
   armHHos[153] + armHHos[152] + armHHos[82];
   armHHos[94]=armHHos[57]*armHHos[94];
   armHHos[102]=armHHos[115] + armHHos[123] - armHHos[67] + 
   armHHos[140];
   armHHos[102]=mmH*armHHos[102];
   armHHos[97]=175./3. + armHHos[97];
   armHHos[107]=armHHos[97] + armHHos[125];
   armHHos[107]=armHHos[7]*armHHos[107];
   armHHos[95]=armHHos[95] + 13*armHHos[75] + 451./8. - 21*armHHos[73];
   armHHos[95]=1./2.*armHHos[95] + 13*armHHos[77];
   armHHos[97]=armHHos[97] + 3*armHHos[9];
   armHHos[97]=armHHos[9]*armHHos[97];
   armHHos[95]=1./2.*armHHos[102] + 1./32.*armHHos[97] + 1./16.*
   armHHos[107] + 81./32.*armHHos[28] + armHHos[176] - 5./4.*
   armHHos[12] + 17./8.*armHHos[79] + 63./16.*armHHos[30] + 17./4.*
   armHHos[80] - armHHos[74] + 1./4.*armHHos[95] + armHHos[121];
   armHHos[95]=mmH*armHHos[95];
   armHHos[97]=armHHos[162] + armHHos[149] + armHHos[126];
   armHHos[97]=armHHos[87]*armHHos[15]*armHHos[97];
   armHHos[102]=1543./3. + armHHos[128];
   armHHos[107]=armHHos[102] - 45*armHHos[7];
   armHHos[107]=1./2.*armHHos[107] + armHHos[139];
   armHHos[107]=1./4.*armHHos[107] + armHHos[205];
   armHHos[107]=armHHos[14]*armHHos[107];
   armHHos[102]=1./2.*armHHos[102] - 51*armHHos[7];
   armHHos[102]=1./2.*armHHos[102] + 9*armHHos[147];
   armHHos[102]=armHHos[15]*armHHos[102];
   armHHos[110]=armHHos[145] - 11./2.*armHHos[29] - 51*armHHos[14];
   armHHos[110]=armHHos[9]*armHHos[110];
   armHHos[113]=1./2.*armHHos[82];
   armHHos[95]=armHHos[95] + 1./4.*armHHos[97] + 1./8.*armHHos[110] + 1.
   /4.*armHHos[102] + armHHos[107] + armHHos[156] + 13./8.*armHHos[29]
    - 79./24.*armHHos[17] - 79./12.*armHHos[16] + 27./16.*armHHos[32]
    - 17./16.*armHHos[60] + armHHos[113] + 45./8.*armHHos[58] - 171./8.
   *armHHos[84] + armHHos[83] + armHHos[151];
   armHHos[95]=armHHos[57]*armHHos[95];
   armHHos[97]=armHHos[10] + armHHos[138];
   armHHos[102]=armHHos[56]*armHHos[97];
   armHHos[97]=armHHos[1]*armHHos[97];
   armHHos[97]=1./2.*armHHos[136] + armHHos[134] + armHHos[102] + 1./3.
   *armHHos[97];
   armHHos[97]=mmH*armHHos[97];
   armHHos[97]=armHHos[97] + armHHos[127] + armHHos[108];
   armHHos[102]=armHHos[119] - 49./4.*armHHos[14];
   armHHos[102]=armHHos[14]*armHHos[102];
   armHHos[107]=7./12.*armHHos[15] + 5./3.*armHHos[29] - 21./2.*
   armHHos[14];
   armHHos[107]=armHHos[15]*armHHos[107];
   armHHos[102]=1./3.*armHHos[102] + armHHos[107];
   armHHos[102]=armHHos[4]*armHHos[57]*armHHos[102];
   armHHos[95]=armHHos[102] + 1./2.*armHHos[97] + armHHos[95];
   armHHos[95]=armHHos[103]*armHHos[95];
   armHHos[94]=armHHos[95] + armHHos[146] + armHHos[104] + armHHos[94];
   armHHos[94]=armHHos[103]*armHHos[94];
   armHHos[95]=101./8. + 9*armHHos[71];
   armHHos[95]=1./2.*armHHos[95] + armHHos[131];
   armHHos[97]= - 1./2. + armHHos[7];
   armHHos[99]=armHHos[97] + armHHos[99];
   armHHos[99]=armHHos[9]*armHHos[99];
   armHHos[92]=armHHos[28]*armHHos[92];
   armHHos[97]=armHHos[7]*armHHos[97];
   armHHos[102]=armHHos[112] - 27./4.*armHHos[62] - armHHos[64];
   armHHos[102]= - 1./4.*armHHos[65] + 3*armHHos[102] - 1./2.*
   armHHos[67];
   armHHos[102]=mmH*armHHos[102];
   armHHos[92]=1./4.*armHHos[102] + 1./8.*armHHos[99] + 1./4.*
   armHHos[97] + 27./16.*armHHos[92] + armHHos[106] + armHHos[133] - 27.
   /8.*armHHos[30] + 3./4.*armHHos[74] + 1./2.*armHHos[78] + 3./2.*
   armHHos[76] + 9./16.*armHHos[72] + 3./8.*armHHos[95] + armHHos[81];
   armHHos[92]=mmH*armHHos[92];
   armHHos[95]=3*armHHos[109];
   armHHos[97]=33./2.*armHHos[32] + 9./2.*armHHos[60] - 11./2.*
   armHHos[82] + 9*armHHos[59] - 9./2.*armHHos[61] - 11*armHHos[83];
   armHHos[99]=17./4. - 9*armHHos[28];
   armHHos[99]=armHHos[29]*armHHos[99];
   armHHos[102]=armHHos[95] + 3./2.*armHHos[99] - 11*armHHos[17] + 1./2.
   *armHHos[97] + 17*armHHos[16];
   armHHos[104]=15./8.*armHHos[7] - 2 + armHHos[105];
   armHHos[104]=armHHos[14]*armHHos[104];
   armHHos[105]= - 11*armHHos[7] + 23./2. + armHHos[98];
   armHHos[105]=armHHos[15]*armHHos[105];
   armHHos[106]=3./4.*armHHos[29] + armHHos[14];
   armHHos[106]=armHHos[9]*armHHos[106];
   armHHos[102]=1./2.*armHHos[92] + 1./4.*armHHos[106] + 1./8.*
   armHHos[105] + 1./8.*armHHos[102] + armHHos[104];
   armHHos[102]=mmH*armHHos[102];
   armHHos[104]=1./2.*armHHos[29] + armHHos[14];
   armHHos[104]=armHHos[14]*armHHos[104];
   armHHos[105]=7./2.*armHHos[29] + armHHos[194];
   armHHos[105]=armHHos[15]*armHHos[105];
   armHHos[106]=pow(armHHos[29],2);
   armHHos[104]=1./2.*armHHos[105] - 117./16.*armHHos[106] + 7*
   armHHos[104];
   armHHos[102]=1./4.*armHHos[104] + armHHos[102];
   armHHos[102]=armHHos[57]*armHHos[102];
   armHHos[104]=armHHos[14]*armHHos[142];
   armHHos[105]=armHHos[15]*armHHos[186];
   armHHos[104]=armHHos[105] + armHHos[104] + armHHos[16] - armHHos[17]
   ;
   armHHos[104]=mmH*armHHos[104];
   armHHos[105]= - armHHos[15]*armHHos[14];
   armHHos[104]=armHHos[104] + armHHos[114] + armHHos[105];
   armHHos[104]=armHHos[100]*armHHos[57]*armHHos[104];
   armHHos[102]=armHHos[102] + 1./8.*armHHos[104];
   armHHos[102]=armHHos[100]*armHHos[102];
   armHHos[104]=3*armHHos[99] - armHHos[17] + armHHos[97] + 41*
   armHHos[16];
   armHHos[95]=1./2.*armHHos[104] + armHHos[95];
   armHHos[104]=armHHos[188] + 1 + armHHos[98];
   armHHos[104]=armHHos[15]*armHHos[104];
   armHHos[105]=11./4.*armHHos[7] - 13./4. + 3*armHHos[28];
   armHHos[105]=armHHos[14]*armHHos[105];
   armHHos[95]=1./2.*armHHos[104] + 1./2.*armHHos[95] + 3*armHHos[105];
   armHHos[104]=3./8.*armHHos[29] + armHHos[14];
   armHHos[104]=armHHos[9]*armHHos[104];
   armHHos[95]=armHHos[92] + 1./2.*armHHos[95] + armHHos[104];
   armHHos[95]=mmH*armHHos[95];
   armHHos[104]=7*armHHos[29];
   armHHos[105]=armHHos[104] + 27./2.*armHHos[14];
   armHHos[105]=armHHos[14]*armHHos[105];
   armHHos[107]=armHHos[29] + armHHos[143];
   armHHos[107]=armHHos[15]*armHHos[107];
   armHHos[106]= - 117./8.*armHHos[106];
   armHHos[105]=7./2.*armHHos[107] + armHHos[106] + armHHos[105];
   armHHos[95]=1./4.*armHHos[105] + armHHos[95];
   armHHos[95]=armHHos[57]*armHHos[95];
   armHHos[102]=armHHos[95] + armHHos[102];
   armHHos[102]=armHHos[100]*armHHos[102];
   armHHos[98]= - 35./2. + armHHos[98];
   armHHos[105]=9*armHHos[7];
   armHHos[107]=armHHos[98] + armHHos[105];
   armHHos[107]=armHHos[14]*armHHos[107];
   armHHos[98]=1./2.*armHHos[98] + armHHos[105];
   armHHos[98]=armHHos[15]*armHHos[98];
   armHHos[105]=1./4.*armHHos[29] + 3*armHHos[14];
   armHHos[105]=armHHos[9]*armHHos[105];
   armHHos[92]=armHHos[92] + 3./2.*armHHos[105] + 1./2.*armHHos[98] + 1.
   /2.*armHHos[107] + 3./4.*armHHos[109] + 3./8.*armHHos[99] + 9./2.*
   armHHos[17] + 1./8.*armHHos[97] + 9*armHHos[16];
   armHHos[92]=mmH*armHHos[92];
   armHHos[97]=armHHos[104] + armHHos[14];
   armHHos[97]=armHHos[14]*armHHos[97];
   armHHos[97]=armHHos[106] + armHHos[97];
   armHHos[98]=7./4.*armHHos[29] + armHHos[14];
   armHHos[98]=armHHos[15]*armHHos[98];
   armHHos[97]=1./2.*armHHos[97] + armHHos[98];
   armHHos[92]=1./2.*armHHos[97] + armHHos[92];
   armHHos[92]=armHHos[103]*armHHos[57]*armHHos[92];
   armHHos[92]=armHHos[95] + 1./2.*armHHos[92];
   armHHos[92]=armHHos[103]*armHHos[92];
   armHHos[95]=armHHos[129] + armHHos[12] + armHHos[137] + 9./8. - 
   armHHos[76];
   armHHos[95]=mmH*armHHos[95];
   armHHos[97]=armHHos[83] - 1./2.*armHHos[59];
   armHHos[98]=1./2.*armHHos[118];
   armHHos[99]= - 3 + armHHos[122];
   armHHos[99]=armHHos[14]*armHHos[99];
   armHHos[104]= - 3./2.*armHHos[15];
   armHHos[105]= - armHHos[29] + armHHos[15];
   armHHos[106]=armHHos[9]*armHHos[105];
   armHHos[107]=1./4.*armHHos[106];
   armHHos[95]=1./2.*armHHos[95] + armHHos[107] + armHHos[104] + 
   armHHos[99] + armHHos[98] - 3./4.*armHHos[29] + armHHos[111] + 
   armHHos[16] - 3./4.*armHHos[32] - 1./4.*armHHos[60] + armHHos[97] + 
   armHHos[113];
   armHHos[95]=mmH*armHHos[95];
   armHHos[108]= - armHHos[29] + armHHos[14];
   armHHos[108]=armHHos[14]*armHHos[108];
   armHHos[105]=armHHos[15]*armHHos[105];
   armHHos[109]=armHHos[108] + 1./2.*armHHos[105];
   armHHos[95]=1./2.*armHHos[109] + armHHos[95];
   armHHos[95]=armHHos[103]*armHHos[57]*mmH*armHHos[95];
   armHHos[109]=armHHos[124] + 3*armHHos[97] + armHHos[82];
   armHHos[110]=1 - armHHos[76];
   armHHos[96]=armHHos[13] + armHHos[96] + 3*armHHos[110] - armHHos[74]
   ;
   armHHos[96]=mmH*armHHos[96];
   armHHos[96]=1./4.*armHHos[96] + armHHos[107] + armHHos[104] + 3./2.*
   armHHos[99] + 3./4.*armHHos[118] - armHHos[29] + armHHos[111] + 3./2.
   *armHHos[16] + 1./2.*armHHos[109] - armHHos[32];
   armHHos[96]=mmH*armHHos[96];
   armHHos[104]=3*armHHos[108] + armHHos[105];
   armHHos[96]=1./4.*armHHos[104] + armHHos[96];
   armHHos[96]=armHHos[57]*mmH*armHHos[96];
   armHHos[95]=armHHos[96] + 1./2.*armHHos[95];
   armHHos[95]=armHHos[103]*armHHos[95];
   armHHos[103]=1./4.*armHHos[13] + armHHos[12] - 1./4.*armHHos[74] + 
   15./16. - armHHos[76];
   armHHos[103]=mmH*armHHos[103];
   armHHos[103]=1./2.*armHHos[103] + 1./8.*armHHos[106] - 3./4.*
   armHHos[15] + armHHos[99] + armHHos[98] - 5./8.*armHHos[29] + 1./4.*
   armHHos[17] + armHHos[16] - 5./8.*armHHos[32] - 1./8.*armHHos[60] + 
   armHHos[97] + 1./4.*armHHos[82];
   armHHos[103]=mmH*armHHos[103];
   armHHos[104]=armHHos[108] + 1./4.*armHHos[105];
   armHHos[103]=1./2.*armHHos[104] + armHHos[103];
   armHHos[103]=armHHos[57]*mmH*armHHos[103];
   armHHos[104]=armHHos[12] + 3./4. - armHHos[76];
   armHHos[104]=mmH*armHHos[104];
   armHHos[97]=1./2.*armHHos[104] + armHHos[99] + armHHos[98] + 
   armHHos[149] + armHHos[16] + armHHos[97] - 1./2.*armHHos[32];
   armHHos[97]=mmH*armHHos[97];
   armHHos[97]=1./2.*armHHos[108] + armHHos[97];
   armHHos[97]=armHHos[100]*armHHos[57]*mmH*armHHos[97];
   armHHos[97]=armHHos[103] + 1./2.*armHHos[97];
   armHHos[97]=armHHos[100]*armHHos[97];
   armHHos[96]=armHHos[96] + armHHos[97];
   armHHos[96]=armHHos[100]*armHHos[96];
   armHHos[95]=armHHos[96] + armHHos[95];
   armHHos[95]=armHHos[8]*armHHos[95];
   armHHos[92]=1./2.*armHHos[95] + armHHos[102] + armHHos[92];
   armHHos[92]=armHHos[8]*armHHos[92];
   armHHos[92]=armHHos[92] + armHHos[93] + armHHos[94];
   armHHos[92]=armHHos[8]*armHHos[92];
   armHHos[93]=pow(c,2);
   armHHos[94]= - 1 - armHHos[93];
   armHHos[95]=armHHos[70]*armHHos[94];
   armHHos[94]=armHHos[69]*armHHos[94];
   armHHos[94]=2*armHHos[95] + armHHos[94];
   armHHos[94]=8*armHHos[66] + 4*armHHos[94] - armHHos[67];
   armHHos[95]=41 + 12*armHHos[93];
   armHHos[95]=armHHos[68]*armHHos[95];
   armHHos[94]=3*armHHos[94] + armHHos[95];
   armHHos[94]=mmZ*armHHos[94];
   armHHos[95]=67./3. - 6*armHHos[77];
   armHHos[96]=89 + armHHos[120];
   armHHos[96]=armHHos[7]*armHHos[96];
   armHHos[94]=armHHos[94] + 2./3.*armHHos[96] + 6*armHHos[12] + 12*
   armHHos[79] + 2*armHHos[95] + 143./3.*armHHos[80];
   armHHos[94]=mmZ*armHHos[94];
   armHHos[95]=armHHos[15]*armHHos[142];
   armHHos[95]=armHHos[17] + armHHos[95];
   armHHos[94]=12*armHHos[95] + armHHos[94];
   armHHos[94]=armHHos[57]*armHHos[94];
   armHHos[95]= - 47 - 12*armHHos[93];
   armHHos[95]=armHHos[80]*armHHos[95];
   armHHos[96]=armHHos[7]*armHHos[186];
   armHHos[93]=9./2.*armHHos[96] - 12*armHHos[79] + 2*armHHos[95] + 6*
   armHHos[81] - 97 - 24*armHHos[93];
   armHHos[93]=mmZ*armHHos[93];
   armHHos[95]=2*armHHos[15] + armHHos[135] - 2*armHHos[84] + 
   armHHos[58];
   armHHos[93]=12*armHHos[95] + armHHos[93];
   armHHos[93]=armHHos[57]*mmZ*armHHos[93];
   armHHos[95]= - armHHos[10] - 1 - armHHos[21];
   armHHos[96]=armHHos[56]*armHHos[95];
   armHHos[95]=armHHos[1]*armHHos[95];
   armHHos[95]=3*armHHos[96] + armHHos[95];
   armHHos[95]=armHHos[101]*armHHos[95];
   armHHos[93]=armHHos[95] + armHHos[93];
   armHHos[93]=armHHos[4]*armHHos[93];
   armHHos[93]=armHHos[94] + armHHos[93];
   armHHos[93]=armHHos[4]*armHHos[93];
   armHHos[94]= - 5*armHHos[68] + 4*armHHos[70] + 5*armHHos[69];
   armHHos[95]= - 4 + armHHos[7];
   armHHos[95]=armHHos[7]*armHHos[95];
   armHHos[95]=2*armHHos[95] - 4*armHHos[80] + 1 + 2*armHHos[77];
   armHHos[95]=armHHos[88]*armHHos[95];
   armHHos[96]= - 1 - armHHos[80];
   armHHos[96]=armHHos[87]*armHHos[96];
   armHHos[94]=3*armHHos[96] + 2*armHHos[94] + 3*armHHos[95];
   armHHos[94]=mmZ*armHHos[94];
   armHHos[95]= - 3 - 2*armHHos[80];
   armHHos[94]=3*armHHos[95] + armHHos[94];
   armHHos[94]=armHHos[57]*armHHos[94];

      mHHosret = armHHos[89] + armHHos[90] + armHHos[91] + armHHos[92]
       + armHHos[93] + armHHos[94];
      return mHHosret;
}
