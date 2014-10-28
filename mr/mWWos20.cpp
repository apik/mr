#include <WW.hpp>
std::complex<long double>
WW<MS>::m20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armWWos[222], mWWosret;

    armWWos[1]=double(nL + nH);
    armWWos[2]=pow(mmZ,-1);
    armWWos[3]=pow(s,-1);
    armWWos[4]=pow(c,-1);
    armWWos[5]=pow(mmH,-1);
    armWWos[6]=std::real(Tsil::B(0,0,mmW,mu2));
    armWWos[7]=Tsil::I2(0,0,mmW,mu2);
    armWWos[8]=Tsil::I2(0,0,mmZ,mu2);
    armWWos[9]=Tsil::B(mmW,mmZ,mmW,mu2);
    armWWos[10]=Tsil::B(mmW,mmH,mmW,mu2);
    armWWos[11]=Tsil::B(0,mmZ,mmW,mu2);
    armWWos[12]=Tsil::B(0,mmH,mmW,mu2);
    armWWos[13]=Tsil::Beps(mmW,mmZ,mmW,mu2);
    armWWos[14]=Tsil::Beps(mmW,mmH,mmW,mu2);
    armWWos[15]=Tsil::A(mmW,mu2);
    armWWos[16]=Tsil::A(mmZ,mu2);
    armWWos[17]=Tsil::A(mmH,mu2);
    armWWos[18]=Tsil::Aeps(mmW,mu2);
    armWWos[19]=Tsil::Aeps(mmZ,mu2);
    armWWos[20]=Tsil::Aeps(mmH,mu2);
    armWWos[21]=protW0Z00->M(0);
    armWWos[22]=prot0W0Z0->M(0);
    armWWos[23]=prot00W00->M(0);
    armWWos[24]=prot0000Z->M(0);
    armWWos[25]=protHW00->Uxzuv(0);
    armWWos[26]=protW0Z00->Uzxyv(0);
    armWWos[27]=protWtZ00->Uxzuv(0);
    armWWos[28]=protHW00->Txuv(0);
    armWWos[29]=prot0W0Z0->Tuxv(0);
    armWWos[30]=protWtZ00->Txuv(0);
    armWWos[31]=double(nH);
    armWWos[32]=pow(mmt,-1);
    armWWos[33]=Tsil::B(0,mmt,mmW,mu2);
    armWWos[34]=Tsil::A(mmt,mu2);
    armWWos[35]=double(nL);
    armWWos[36]=Tsil::I2(mmZ,mmt,mmt,mu2);
    armWWos[37]=Tsil::I2(mmH,mmt,mmt,mu2);
    armWWos[38]=Tsil::I2(0,mmW,mmt,mu2);
    armWWos[39]=Tsil::Beps(0,mmt,mmW,mu2);
    armWWos[40]=Tsil::Aeps(mmt,mu2);
    armWWos[41]=prot00tt0->M(0);
    armWWos[42]=prot00tt0->Tuxv(0);
    armWWos[43]=protWtZ00->M(0);
    armWWos[44]=protW0Htt->M(0);
    armWWos[45]=protW0Ztt->M(0);
    armWWos[46]=prot0Wt0t->M(0);
    armWWos[47]=prot00Wt0->M(0);
    armWWos[48]=prot00ttZ->M(0);
    armWWos[49]=protW0Htt->Uzxyv(0);
    armWWos[50]=protWtZ00->Uzxyv(0);
    armWWos[51]=protW0Htt->Uxzuv(0);
    armWWos[52]=protW0Ztt->Uxzuv(0);
    armWWos[53]=prot0Wt0t->Uyuzv(0);
    armWWos[54]=protW0Htt->Uyuzv(0);
    armWWos[55]=protW0Ztt->Uyuzv(0);
    armWWos[56]=protWtZ00->Uuyxv(0);
    armWWos[57]=protW0Htt->Tzyv(0);
    armWWos[58]=protWtZ00->Tzyv(0);
    armWWos[59]=protW0Htt->Tvyz(0);
    armWWos[60]=protWtZ00->Tyzv(0);
    armWWos[61]=protW0Htt->Suxv(0);
    armWWos[62]=protW0Htt->Svyz(0);
    armWWos[63]=protWtZ00->Svyz(0);
    armWWos[64]=prot00000->M(0);
    armWWos[65]=double(boson);
    armWWos[66]=Tsil::I2(mmW,mmW,mmZ,mu2);
    armWWos[67]=Tsil::I2(mmH,mmW,mmW,mu2);
    armWWos[68]=Tsil::I2(mmH,mmZ,mmZ,mu2);
    armWWos[69]=Tsil::I2(mmH,mmH,mmH,mu2);
    armWWos[70]=protWHHWW->M(0);
    armWWos[71]=protWHZWW->M(0);
    armWWos[72]=protWZZWW->M(0);
    armWWos[73]=protWWHHH->M(0);
    armWWos[74]=protWWHZZ->M(0);
    armWWos[75]=protW0HWW->M(0);
    armWWos[76]=protW0ZWW->M(0);
    armWWos[77]=prot0WW0W->M(0);
    armWWos[78]=protWH0H->Vxzuv(0);
    armWWos[79]=protWZ0Z->Vxzuv(0);
    armWWos[80]=protWHHWW->Uzxyv(0);
    armWWos[81]=protWHZWW->Uyuzv(0);
    armWWos[82]=protWHZWW->Uzxyv(0);
    armWWos[83]=protWZZWW->Uzxyv(0);
    armWWos[84]=protWWHHH->Uyuzv(0);
    armWWos[85]=protWWHZZ->Uxzuv(0);
    armWWos[86]=protWHHWW->Uxzuv(0);
    armWWos[87]=protWWHZZ->Uyuzv(0);
    armWWos[88]=protWHZWW->Uxzuv(0);
    armWWos[89]=protWHZWW->Tyzv(0);
    armWWos[90]=protW0HWW->Tzyv(0);
    armWWos[91]=protWHZWW->Tzyv(0);
    armWWos[92]=protW0ZWW->Tzyv(0);
    armWWos[93]=protWHHWW->Svyz(0);
    armWWos[94]=protWHZWW->Svyz(0);
    armWWos[95]=protWZZWW->Svyz(0);
    armWWos[96]=protWWZZH->M(0);
    armWWos[97]=1/(mmt - mmW);
    armWWos[98]=1/(4*mmt - mmZ);
    armWWos[99]=1/(mmH - 4*mmZ);
    armWWos[100]=1/(mmH - 4*mmW);
   armWWos[101]= - 659./36. - 85*armWWos[9];
   armWWos[102]=2./3. + armWWos[6];
   armWWos[103]=1./3.*armWWos[1]*armWWos[102];
   armWWos[102]=armWWos[35]*armWWos[102];
   armWWos[104]=armWWos[100]*armWWos[17];
   armWWos[105]=1./2.*armWWos[104];
   armWWos[101]=armWWos[105] + armWWos[10] + armWWos[102] + 1./8.*
   armWWos[101] + armWWos[103];
   armWWos[101]=armWWos[33]*armWWos[101];
   armWWos[106]= - 3*armWWos[37];
   armWWos[107]=8*armWWos[40] + armWWos[106] - 2*armWWos[38];
   armWWos[108]=armWWos[17] - 9./4.*armWWos[16];
   armWWos[108]=armWWos[33]*armWWos[108];
   armWWos[109]=2 - armWWos[33];
   armWWos[109]=armWWos[15]*armWWos[109];
   armWWos[110]=3*armWWos[34];
   armWWos[107]=armWWos[110] + 3*armWWos[109] + armWWos[108] - 
   armWWos[17] + 6*armWWos[18] + 3*armWWos[107] + 10*armWWos[20];
   armWWos[107]=armWWos[5]*armWWos[107];
   armWWos[108]= - 4*armWWos[53];
   armWWos[109]= - 31./3.*armWWos[42];
   armWWos[111]=armWWos[109] - 89./12. + armWWos[108];
   armWWos[112]=5./2.*armWWos[20];
   armWWos[106]= - 5./2.*armWWos[17] + armWWos[18] + armWWos[112] + 6*
   armWWos[40] - armWWos[38] + armWWos[106] + 2*armWWos[62];
   armWWos[106]=armWWos[100]*armWWos[106];
   armWWos[113]=armWWos[47] + 1./3.*armWWos[41];
   armWWos[114]= - 2*armWWos[46];
   armWWos[115]=armWWos[113] + armWWos[114];
   armWWos[116]= - armWWos[48] + 1./2.*armWWos[43] + 2*armWWos[115] + 1.
   /2.*armWWos[45];
   armWWos[117]= - 5 + armWWos[59];
   armWWos[117]=armWWos[100]*armWWos[117];
   armWWos[116]= - 9*armWWos[5] + armWWos[117] + 1./3.*armWWos[116] + 2
   *armWWos[44];
   armWWos[116]=mmt*armWWos[116];
   armWWos[118]= - 41./8.*armWWos[50];
   armWWos[119]= - 7./3.*armWWos[56];
   armWWos[120]= - 17./32.*armWWos[55];
   armWWos[121]=1./6.*armWWos[52];
   armWWos[122]= - 9./4.*armWWos[54];
   armWWos[123]=6*armWWos[51];
   armWWos[124]= - 1./6.*armWWos[27];
   armWWos[125]=35./6.*armWWos[58];
   armWWos[126]= - 311./32.*armWWos[60];
   armWWos[127]= - 17./12.*armWWos[30];
   armWWos[128]=491./96.*armWWos[39];
   armWWos[129]=73./24.*armWWos[59];
   armWWos[130]=5./2.*armWWos[57];
   armWWos[131]=41./8.*armWWos[13];
   armWWos[132]= - 6*armWWos[14];
   armWWos[133]= - armWWos[33]*armWWos[100];
   armWWos[134]= - 6*armWWos[100] + armWWos[133];
   armWWos[134]=armWWos[15]*armWWos[134];
   armWWos[135]=armWWos[34]*armWWos[100];
   armWWos[136]=4*armWWos[135];
   armWWos[137]= - mmH*armWWos[44];
   armWWos[138]=3*armWWos[137];
   armWWos[139]= - 2*armWWos[10];
   armWWos[111]=armWWos[116] + armWWos[138] + armWWos[107] + 
   armWWos[136] + armWWos[134] + armWWos[101] + armWWos[106] + 
   armWWos[139] + armWWos[132] + armWWos[131] + armWWos[130] + 
   armWWos[129] + armWWos[128] + armWWos[127] + armWWos[126] + 
   armWWos[125] + armWWos[124] + armWWos[123] + armWWos[122] + 
   armWWos[121] + armWWos[120] + armWWos[119] + 1./3.*armWWos[111] + 
   armWWos[118];
   armWWos[111]=mmt*armWWos[111];
   armWWos[140]= - 9*armWWos[9];
   armWWos[141]=559./72. + armWWos[140];
   armWWos[142]=5*armWWos[10];
   armWWos[143]= - 451./96.*armWWos[33];
   armWWos[144]= - armWWos[15]*armWWos[100];
   armWWos[145]=8*armWWos[144];
   armWWos[146]= - armWWos[34]*armWWos[100];
   armWWos[147]=3*armWWos[146];
   armWWos[141]=armWWos[147] + armWWos[145] + armWWos[143] + 
   armWWos[104] + armWWos[142] + armWWos[102] + 7./4.*armWWos[141] + 
   armWWos[103];
   armWWos[141]=armWWos[34]*armWWos[141];
   armWWos[148]=43*armWWos[63];
   armWWos[149]=armWWos[148] - 4175./12.*armWWos[36];
   armWWos[150]=13*armWWos[61];
   armWWos[149]=1./8.*armWWos[149] + armWWos[150];
   armWWos[151]=469./12. - 3*armWWos[54];
   armWWos[152]= - 3*armWWos[51];
   armWWos[151]=7./2.*armWWos[14] + 21./4.*armWWos[57] - 5./2.*
   armWWos[59] + 3./4.*armWWos[39] - 1./2.*armWWos[49] + 1./4.*
   armWWos[151] + armWWos[152];
   armWWos[153]= - 5 - 7*armWWos[10];
   armWWos[153]=armWWos[33]*armWWos[153];
   armWWos[154]=1./12.*armWWos[153] + 1./2.*armWWos[151] + 8./3.*
   armWWos[10];
   armWWos[154]=mmH*armWWos[154];
   armWWos[155]=71./12.*armWWos[62];
   armWWos[156]=5./12.*armWWos[8];
   armWWos[157]= - 23./8.*armWWos[38];
   armWWos[158]= - 115./24.*armWWos[20];
   armWWos[159]= - 109./24.*armWWos[17];
   armWWos[160]= - 11*armWWos[17];
   armWWos[161]=armWWos[160] - 101./24.*armWWos[16];
   armWWos[161]=1./8.*armWWos[33]*armWWos[161];
   armWWos[162]= - 2*armWWos[33];
   armWWos[163]=109./24. + armWWos[162];
   armWWos[163]=armWWos[15]*armWWos[163];
   armWWos[164]= - 3*armWWos[34] - 27./2.*armWWos[15] - 12*armWWos[17]
    + 1./4.*armWWos[16];
   armWWos[164]=armWWos[5]*armWWos[34]*armWWos[164];
   armWWos[165]=3./2.*armWWos[37];
   armWWos[111]=armWWos[111] + armWWos[154] + armWWos[164] + 
   armWWos[141] + armWWos[163] + armWWos[161] - 6463./576.*armWWos[16]
    + armWWos[159] + 253./72.*armWWos[18] + armWWos[158] + 6571./576.*
   armWWos[19] + 683./288.*armWWos[40] + armWWos[157] + armWWos[156] + 
   armWWos[155] + 1./3.*armWWos[149] + armWWos[165];
   armWWos[111]=mmt*armWWos[111];
   armWWos[141]= - 59*armWWos[9];
   armWWos[149]= - 775./36. + armWWos[141];
   armWWos[166]=1./2.*armWWos[6];
   armWWos[167]=1./3. + armWWos[166];
   armWWos[168]=armWWos[1]*armWWos[167];
   armWWos[169]=1./3.*armWWos[168];
   armWWos[167]=armWWos[35]*armWWos[167];
   armWWos[170]=1./2.*armWWos[10];
   armWWos[171]=1./4.*armWWos[104];
   armWWos[149]=armWWos[171] + armWWos[170] + armWWos[167] + 1./8.*
   armWWos[149] + armWWos[169];
   armWWos[149]=armWWos[33]*armWWos[149];
   armWWos[172]= - 3./2.*armWWos[37];
   armWWos[173]=4*armWWos[40] + armWWos[172] - armWWos[38];
   armWWos[174]=3*armWWos[18];
   armWWos[175]= - 1./2.*armWWos[17];
   armWWos[173]=armWWos[175] + armWWos[174] + 3*armWWos[173] + 5*
   armWWos[20];
   armWWos[176]= - 3./2.*armWWos[16];
   armWWos[177]=armWWos[17] + armWWos[176];
   armWWos[177]=armWWos[33]*armWWos[177];
   armWWos[178]= - 1./2.*armWWos[33];
   armWWos[179]=1 + armWWos[178];
   armWWos[179]=3*armWWos[15]*armWWos[179];
   armWWos[180]=3./2.*armWWos[34];
   armWWos[177]=armWWos[180] + armWWos[179] + armWWos[173] + 1./2.*
   armWWos[177];
   armWWos[177]=armWWos[5]*armWWos[177];
   armWWos[181]=3*armWWos[40];
   armWWos[182]= - 1./2.*armWWos[38];
   armWWos[183]=1./2.*armWWos[18];
   armWWos[172]= - 5./4.*armWWos[17] + armWWos[183] + 5./4.*armWWos[20]
    + armWWos[181] + armWWos[182] + armWWos[172] + armWWos[62];
   armWWos[172]=armWWos[100]*armWWos[172];
   armWWos[184]= - armWWos[45] + armWWos[43];
   armWWos[117]=1./2.*armWWos[117];
   armWWos[185]= - 9./2.*armWWos[5];
   armWWos[184]=armWWos[185] + armWWos[117] + 1./2.*armWWos[184] + 
   armWWos[44];
   armWWos[184]=mmt*armWWos[184];
   armWWos[186]=5./8.*armWWos[55] - 5*armWWos[56] + 947./24. - 13*
   armWWos[50];
   armWWos[186]=armWWos[122] + 1./2.*armWWos[186] - armWWos[52];
   armWWos[187]=1./4.*armWWos[58];
   armWWos[188]=3*armWWos[51];
   armWWos[189]=73./48.*armWWos[59];
   armWWos[190]=5./4.*armWWos[57];
   armWWos[191]= - 3*armWWos[14];
   armWWos[192]= - 3*armWWos[100] + 1./2.*armWWos[133];
   armWWos[192]=armWWos[15]*armWWos[192];
   armWWos[135]=2*armWWos[135];
   armWWos[137]=3./2.*armWWos[137];
   armWWos[149]=armWWos[184] + armWWos[137] + armWWos[177] + 
   armWWos[135] + armWWos[192] + armWWos[149] + armWWos[172] - 
   armWWos[10] + armWWos[191] + 17./4.*armWWos[13] + armWWos[190] + 
   armWWos[189] + 71./32.*armWWos[39] - 3./8.*armWWos[30] - 195./32.*
   armWWos[60] + armWWos[187] - 1./2.*armWWos[27] + 1./2.*armWWos[186]
    + armWWos[188];
   armWWos[149]=mmt*armWWos[149];
   armWWos[177]= - 1921./36. - 93*armWWos[9];
   armWWos[184]=5./2.*armWWos[10];
   armWWos[186]=4*armWWos[144];
   armWWos[146]=3./2.*armWWos[146];
   armWWos[177]=armWWos[146] + armWWos[186] - 23./32.*armWWos[33] + 
   armWWos[105] + armWWos[184] + armWWos[167] + 1./8.*armWWos[177] + 
   armWWos[169];
   armWWos[177]=armWWos[34]*armWWos[177];
   armWWos[193]= - 2*armWWos[17] + 9./4.*armWWos[16];
   armWWos[194]= - 3./2.*armWWos[34];
   armWWos[193]=armWWos[194] + 3*armWWos[193] + 5*armWWos[15];
   armWWos[193]=armWWos[5]*armWWos[34]*armWWos[193];
   armWWos[151]=1./24.*armWWos[153] + 1./4.*armWWos[151] + 4./3.*
   armWWos[10];
   armWWos[151]=mmH*armWWos[151];
   armWWos[153]=13*armWWos[63] - 45./8.*armWWos[36];
   armWWos[195]=armWWos[160] - 197./4.*armWWos[16];
   armWWos[195]=armWWos[33]*armWWos[195];
   armWWos[196]=35./8. - 8*armWWos[33];
   armWWos[196]=armWWos[15]*armWWos[196];
   armWWos[197]=3./4.*armWWos[37];
   armWWos[198]=71./24.*armWWos[62];
   armWWos[199]= - 115./48.*armWWos[20];
   armWWos[200]= - 109./48.*armWWos[17];
   armWWos[149]=armWWos[149] + armWWos[151] + armWWos[193] + 
   armWWos[177] + 1./3.*armWWos[196] + 1./16.*armWWos[195] - 125./64.*
   armWWos[16] + armWWos[200] - 7./8.*armWWos[18] + armWWos[199] - 159./
   64.*armWWos[19] - 363./32.*armWWos[40] - 17./8.*armWWos[38] + 1./4.*
   armWWos[8] + armWWos[198] + armWWos[197] + 1./4.*armWWos[153] + 2*
   armWWos[61];
   armWWos[149]=mmt*armWWos[149];
   armWWos[153]= - 1 - armWWos[57];
   armWWos[177]= - armWWos[33]*armWWos[10];
   armWWos[153]=armWWos[177] + 13./2.*armWWos[153] - 5./3.*armWWos[10];
   armWWos[153]=mmH*armWWos[153];
   armWWos[177]= - 11./6.*armWWos[40] + armWWos[165] - 11./3.*
   armWWos[62];
   armWWos[193]= - armWWos[33]*armWWos[17];
   armWWos[195]=5./3. + armWWos[33];
   armWWos[196]=armWWos[15]*armWWos[195];
   armWWos[201]=13./2. - armWWos[10];
   armWWos[201]=armWWos[34]*armWWos[201];
   armWWos[153]=1./6.*armWWos[153] + 1./3.*armWWos[201] + 1./6.*
   armWWos[196] + 1./6.*armWWos[193] + 1./2.*armWWos[177] + 17./9.*
   armWWos[17];
   armWWos[153]=mmH*armWWos[153];
   armWWos[177]=23./3.*armWWos[17] - 79./4.*armWWos[16];
   armWWos[177]=15./8.*armWWos[34] + 1./2.*armWWos[177] - 41./3.*
   armWWos[15];
   armWWos[177]=armWWos[34]*armWWos[177];
   armWWos[177]=1./2.*armWWos[177] + armWWos[153];
   armWWos[149]=1./2.*armWWos[177] + armWWos[149];
   armWWos[177]=pow(armWWos[3],2);
   armWWos[149]=armWWos[177]*armWWos[149];
   armWWos[193]=23*armWWos[17];
   armWWos[196]=armWWos[193] - 3617./24.*armWWos[16];
   armWWos[196]= - 3727./96.*armWWos[34] + 1./4.*armWWos[196] - 26./3.*
   armWWos[15];
   armWWos[196]=armWWos[34]*armWWos[196];
   armWWos[111]=armWWos[149] + armWWos[111] + 1./3.*armWWos[196] + 
   armWWos[153];
   armWWos[111]=armWWos[177]*armWWos[111];
   armWWos[149]=armWWos[109] + 935./12. + armWWos[108];
   armWWos[101]=armWWos[116] + armWWos[138] + armWWos[107] + 
   armWWos[136] + armWWos[134] + armWWos[101] + armWWos[106] + 
   armWWos[139] + armWWos[132] + armWWos[131] + armWWos[130] + 
   armWWos[129] + armWWos[128] + armWWos[127] + armWWos[126] + 
   armWWos[125] + armWWos[124] + armWWos[123] + armWWos[122] + 
   armWWos[121] + armWWos[120] + armWWos[119] + 1./3.*armWWos[149] + 
   armWWos[118];
   armWWos[101]=mmt*armWWos[101];
   armWWos[106]= - 4279./72. - 63*armWWos[9];
   armWWos[102]=armWWos[147] + armWWos[145] + armWWos[143] + 
   armWWos[104] + armWWos[142] + armWWos[102] + 1./4.*armWWos[106] + 
   armWWos[103];
   armWWos[102]=armWWos[34]*armWWos[102];
   armWWos[103]=armWWos[148] - 79./12.*armWWos[36];
   armWWos[103]=1./8.*armWWos[103] + armWWos[150];
   armWWos[101]=armWWos[101] + armWWos[154] + armWWos[164] + 
   armWWos[102] + armWWos[163] + armWWos[161] + 1729./576.*armWWos[16]
    + armWWos[159] - 1795./72.*armWWos[18] + armWWos[158] - 1621./576.*
   armWWos[19] - 2503./96.*armWWos[40] + armWWos[157] + armWWos[156] + 
   armWWos[155] + 1./3.*armWWos[103] + armWWos[165];
   armWWos[101]=mmt*armWWos[101];
   armWWos[102]=29./18. - 13*armWWos[9];
   armWWos[102]=armWWos[171] + armWWos[170] + armWWos[167] + 1./4.*
   armWWos[102] + armWWos[169];
   armWWos[102]=armWWos[33]*armWWos[102];
   armWWos[103]= - 3*armWWos[16];
   armWWos[104]=armWWos[17] + armWWos[103];
   armWWos[104]=armWWos[33]*armWWos[104];
   armWWos[104]=armWWos[180] + armWWos[179] + armWWos[173] + 1./2.*
   armWWos[104];
   armWWos[104]=armWWos[5]*armWWos[104];
   armWWos[106]=armWWos[115] + armWWos[45];
   armWWos[106]= - armWWos[48] + 2*armWWos[106] - armWWos[43];
   armWWos[106]=armWWos[185] + armWWos[117] + 1./3.*armWWos[106] + 
   armWWos[44];
   armWWos[106]=mmt*armWWos[106];
   armWWos[107]=armWWos[109] + 6749./288. + armWWos[108];
   armWWos[109]=1./3.*armWWos[27];
   armWWos[102]=armWWos[106] + armWWos[137] + armWWos[104] + 
   armWWos[135] + armWWos[192] + armWWos[102] + armWWos[172] - 
   armWWos[10] + armWWos[191] + 7./8.*armWWos[13] + armWWos[190] + 
   armWWos[189] + 139./48.*armWWos[39] - 25./24.*armWWos[30] - 29./8.*
   armWWos[60] + 67./12.*armWWos[58] + armWWos[109] + armWWos[188] - 9./
   8.*armWWos[54] + 2./3.*armWWos[52] - 11./16.*armWWos[55] - 13./12.*
   armWWos[56] + 1./3.*armWWos[107] - 15./8.*armWWos[50];
   armWWos[102]=mmt*armWWos[102];
   armWWos[104]= - 33*armWWos[9];
   armWWos[106]=47./54. + armWWos[104];
   armWWos[105]=armWWos[146] + armWWos[186] - 191./48.*armWWos[33] + 
   armWWos[105] + armWWos[184] + armWWos[167] + 1./8.*armWWos[106] + 
   armWWos[169];
   armWWos[105]=armWWos[34]*armWWos[105];
   armWWos[106]= - 35*armWWos[63] - 1303./18.*armWWos[36];
   armWWos[106]=1./8.*armWWos[106] + 7*armWWos[61];
   armWWos[107]=armWWos[160] + 271./18.*armWWos[16];
   armWWos[107]=armWWos[33]*armWWos[107];
   armWWos[115]=2*armWWos[33];
   armWWos[116]=37./4. + armWWos[115];
   armWWos[116]=armWWos[15]*armWWos[116];
   armWWos[117]=armWWos[194] - 37./2.*armWWos[15] - 6*armWWos[17] - 13./
   2.*armWWos[16];
   armWWos[117]=armWWos[5]*armWWos[34]*armWWos[117];
   armWWos[102]=armWWos[102] + armWWos[151] + armWWos[117] + 
   armWWos[105] + 1./3.*armWWos[116] + 1./16.*armWWos[107] + 2089./864.
   *armWWos[16] + armWWos[200] - 1331./54.*armWWos[18] + armWWos[199]
    + 1907./864.*armWWos[19] - 2779./432.*armWWos[40] - 3./4.*
   armWWos[38] + 1./6.*armWWos[8] + armWWos[198] + 1./3.*armWWos[106]
    + armWWos[197];
   armWWos[102]=mmt*armWWos[102];
   armWWos[105]=83./6. + 15*armWWos[50];
   armWWos[106]= - 53./36. + 21*armWWos[9];
   armWWos[106]=armWWos[33]*armWWos[106];
   armWWos[107]= - armWWos[33]*armWWos[16];
   armWWos[116]=armWWos[5]*armWWos[107];
   armWWos[105]=3./2.*armWWos[116] + 1./4.*armWWos[106] - 15./4.*
   armWWos[13] - 125./144.*armWWos[39] + 109./48.*armWWos[60] + 67./9.*
   armWWos[58] + armWWos[109] - 1./3.*armWWos[52] + 77./144.*
   armWWos[55] + 1./4.*armWWos[105] + 1./3.*armWWos[56];
   armWWos[106]= - armWWos[45] - armWWos[43];
   armWWos[116]= - 1./3.*armWWos[48];
   armWWos[106]=1./2.*armWWos[106] + armWWos[116];
   armWWos[106]=mmt*armWWos[106];
   armWWos[105]=1./2.*armWWos[105] + 1./3.*armWWos[106];
   armWWos[105]=mmt*armWWos[105];
   armWWos[106]= - 43./9.*armWWos[63] + 3./4.*armWWos[36];
   armWWos[117]= - 1./4. - armWWos[33];
   armWWos[117]=armWWos[15]*armWWos[117];
   armWWos[118]=9*armWWos[9];
   armWWos[119]= - 155./144.*armWWos[33] - 293./144. + armWWos[118];
   armWWos[119]=armWWos[34]*armWWos[119];
   armWWos[120]= - 5./2.*armWWos[16];
   armWWos[121]=armWWos[120] + armWWos[15];
   armWWos[121]=armWWos[5]*armWWos[34]*armWWos[121];
   armWWos[122]=armWWos[33]*armWWos[16];
   armWWos[105]=armWWos[105] + 1./2.*armWWos[121] + 1./2.*armWWos[119]
    + 1./6.*armWWos[117] + 299./576.*armWWos[122] + 523./576.*
   armWWos[16] - 1./8.*armWWos[18] - 13./192.*armWWos[19] + 289./288.*
   armWWos[40] - 1./8.*armWWos[38] + 1./36.*armWWos[8] + 1./8.*
   armWWos[106] + 1./3.*armWWos[61];
   armWWos[105]=mmt*armWWos[105];
   armWWos[106]=1 + armWWos[50];
   armWWos[117]= - 1./2.*armWWos[13];
   armWWos[106]=armWWos[117] + 1./2.*armWWos[60] + 1./2.*armWWos[106]
    + armWWos[58];
   armWWos[119]=armWWos[33]*armWWos[9];
   armWWos[106]=1./2.*armWWos[106] + 1./3.*armWWos[119];
   armWWos[106]=mmt*armWWos[106];
   armWWos[119]=7./3.*armWWos[9];
   armWWos[121]= - 1 + armWWos[119];
   armWWos[121]=armWWos[34]*armWWos[121];
   armWWos[121]=armWWos[40] + armWWos[121];
   armWWos[106]=1./4.*armWWos[121] + armWWos[106];
   armWWos[121]=pow(armWWos[4],2);
   armWWos[106]=armWWos[121]*mmt*armWWos[106];
   armWWos[123]= - 47./16.*armWWos[34] + 205./48.*armWWos[16] - 
   armWWos[15];
   armWWos[123]=armWWos[34]*armWWos[123];
   armWWos[105]=1./2.*armWWos[106] + 1./6.*armWWos[123] + armWWos[105];
   armWWos[105]=armWWos[121]*armWWos[105];
   armWWos[106]=armWWos[193] - 7./18.*armWWos[16];
   armWWos[106]= - 1441./36.*armWWos[34] + 1./2.*armWWos[106] + 3193./9.
   *armWWos[15];
   armWWos[106]=armWWos[34]*armWWos[106];
   armWWos[106]=1./6.*armWWos[106] + armWWos[153];
   armWWos[102]=armWWos[105] + 1./2.*armWWos[106] + armWWos[102];
   armWWos[102]=armWWos[121]*armWWos[102];
   armWWos[105]=armWWos[193] + 479./24.*armWWos[16];
   armWWos[105]=1./4.*armWWos[105] + 230./3.*armWWos[15];
   armWWos[105]=1./3.*armWWos[105] + 41./32.*armWWos[34];
   armWWos[105]=armWWos[34]*armWWos[105];
   armWWos[101]=armWWos[102] + armWWos[101] + armWWos[105] + 
   armWWos[153];
   armWWos[101]=armWWos[121]*armWWos[101];
   armWWos[102]= - 5 + armWWos[33];
   armWWos[102]=armWWos[34]*armWWos[102];
   armWWos[105]=pow(armWWos[34],2);
   armWWos[106]=armWWos[5]*armWWos[105];
   armWWos[102]=armWWos[102] + 12*armWWos[106];
   armWWos[102]=armWWos[5]*armWWos[102];
   armWWos[123]= - 1 + 19./6.*armWWos[33];
   armWWos[124]=pow(armWWos[5],2);
   armWWos[125]=mmt*armWWos[124]*armWWos[34];
   armWWos[126]=36*armWWos[125] + 1./2.*armWWos[123] + 3*armWWos[102];
   armWWos[126]=mmt*armWWos[126];
   armWWos[127]=5./6. + armWWos[33];
   armWWos[127]=armWWos[34]*armWWos[127];
   armWWos[128]= - armWWos[5]*armWWos[105];
   armWWos[126]=armWWos[126] + armWWos[127] + 6*armWWos[128];
   armWWos[126]=mmt*armWWos[126];
   armWWos[126]=armWWos[105] + armWWos[126];
   armWWos[127]=armWWos[177]*armWWos[126];
   armWWos[102]=72*armWWos[125] + armWWos[123] + 6*armWWos[102];
   armWWos[102]=mmt*armWWos[102];
   armWWos[123]=5./3. + armWWos[115];
   armWWos[123]=armWWos[34]*armWWos[123];
   armWWos[102]=armWWos[102] + armWWos[123] + 12*armWWos[128];
   armWWos[102]=mmt*armWWos[102];
   armWWos[102]=2*armWWos[105] + armWWos[102];
   armWWos[123]=armWWos[102] + armWWos[127];
   armWWos[123]=armWWos[177]*armWWos[123];
   armWWos[125]=armWWos[121]*armWWos[126];
   armWWos[102]=armWWos[102] + armWWos[125];
   armWWos[102]=armWWos[121]*armWWos[102];
   armWWos[102]=armWWos[123] + armWWos[102];
   armWWos[102]=armWWos[31]*armWWos[102];
   armWWos[101]=armWWos[102] + armWWos[111] + armWWos[101];
   armWWos[101]=armWWos[31]*armWWos[101];
   armWWos[102]=3*armWWos[16];
   armWWos[111]= - 7./2.*armWWos[17];
   armWWos[123]=armWWos[111] + armWWos[102];
   armWWos[123]=armWWos[33]*armWWos[123];
   armWWos[125]=7./2. + 5*armWWos[54];
   armWWos[125]=armWWos[10] + 9./2.*armWWos[14] + 11./4.*armWWos[57] - 
   9./4.*armWWos[59] - 5./4.*armWWos[39] - 3./2.*armWWos[49] + 1./4.*
   armWWos[125] + armWWos[152];
   armWWos[126]= - 5 - 13*armWWos[10];
   armWWos[126]=armWWos[33]*armWWos[126];
   armWWos[127]=armWWos[125] + 1./6.*armWWos[126];
   armWWos[127]=mmH*armWWos[127];
   armWWos[128]= - 41./6. - 25*armWWos[33];
   armWWos[128]=armWWos[15]*armWWos[128];
   armWWos[129]= - 39./2. - armWWos[33];
   armWWos[129]=armWWos[34]*armWWos[129];
   armWWos[130]= - 15./2.*armWWos[37];
   armWWos[131]=27./4.*armWWos[62];
   armWWos[132]=17*armWWos[40];
   armWWos[134]=7./2.*armWWos[20];
   armWWos[135]= - 9./2.*armWWos[17];
   armWWos[123]=armWWos[127] + 1./2.*armWWos[129] + 1./2.*armWWos[128]
    + 1./2.*armWWos[123] - 121./6.*armWWos[16] + armWWos[135] + 39./4.*
   armWWos[18] + armWWos[134] - 85./12.*armWWos[19] + armWWos[132] - 57.
   /4.*armWWos[38] + armWWos[131] + armWWos[130] + 73./4.*armWWos[63]
    - armWWos[61];
   armWWos[128]=1./2.*armWWos[57];
   armWWos[129]=armWWos[33]*armWWos[17];
   armWWos[136]=armWWos[129] + armWWos[20] - armWWos[17];
   armWWos[136]=armWWos[5]*armWWos[136];
   armWWos[137]=mmH*armWWos[44];
   armWWos[138]=armWWos[137] + 1./2.*armWWos[136] - 5./8.*armWWos[33]
    + armWWos[128] + 15./8.*armWWos[59] + 2*armWWos[39] + 9./2.*
   armWWos[60] - 5./4.*armWWos[54] - 3./4.*armWWos[56] - 9./8. - 
   armWWos[42];
   armWWos[138]=mmt*armWWos[138];
   armWWos[123]=1./2.*armWWos[123] + armWWos[138];
   armWWos[123]=mmt*armWWos[123];
   armWWos[138]=3./2.*armWWos[18] - 3./2.*armWWos[20] - 7./4.*
   armWWos[40] + armWWos[182] - 5./2.*armWWos[62] - armWWos[61] + 9./4.
   *armWWos[37];
   armWWos[142]=armWWos[138] + 5./3.*armWWos[17];
   armWWos[143]=11./4. - armWWos[33];
   armWWos[143]=armWWos[15]*armWWos[143];
   armWWos[145]=1./3.*armWWos[129];
   armWWos[142]=1./3.*armWWos[143] + 1./2.*armWWos[142] + armWWos[145];
   armWWos[146]=1./2.*armWWos[49];
   armWWos[147]= - 3./2.*armWWos[14];
   armWWos[148]=armWWos[147] - 7./4.*armWWos[57] + armWWos[146] - 7./4.
    + armWWos[51];
   armWWos[149]=armWWos[148] - 11./6.*armWWos[10];
   armWWos[150]=armWWos[33]*armWWos[10];
   armWWos[151]=1./3.*armWWos[150];
   armWWos[149]=1./2.*armWWos[149] + armWWos[151];
   armWWos[149]=mmH*armWWos[149];
   armWWos[153]= - 4*armWWos[10];
   armWWos[154]=19./4. + armWWos[153];
   armWWos[154]=armWWos[34]*armWWos[154];
   armWWos[154]=armWWos[149] + armWWos[142] + 1./3.*armWWos[154];
   armWWos[154]=mmH*armWWos[154];
   armWWos[155]= - 5*armWWos[17];
   armWWos[156]= - 105*armWWos[15] + armWWos[155] + 103./3.*armWWos[16]
   ;
   armWWos[157]= - 13*armWWos[34];
   armWWos[156]=1./2.*armWWos[156] + armWWos[157];
   armWWos[156]=armWWos[34]*armWWos[156];
   armWWos[123]=armWWos[123] + 1./4.*armWWos[156] + armWWos[154];
   armWWos[123]=mmt*armWWos[123];
   armWWos[154]=armWWos[111] + 145./9.*armWWos[16];
   armWWos[154]=armWWos[33]*armWWos[154];
   armWWos[156]= - 25./2. - 59*armWWos[33];
   armWWos[156]=armWWos[15]*armWWos[156];
   armWWos[158]= - 25./2. - armWWos[33];
   armWWos[158]=armWWos[34]*armWWos[158];
   armWWos[154]=armWWos[127] + 1./2.*armWWos[158] + 1./6.*armWWos[156]
    + 1./2.*armWWos[154] - 833./36.*armWWos[16] + armWWos[135] + 7*
   armWWos[18] + armWWos[134] - 67./36.*armWWos[19] + armWWos[132] - 31.
   /2.*armWWos[38] + armWWos[131] + armWWos[130] + 67./4.*armWWos[63]
    - armWWos[61];
   armWWos[156]=15./16.*armWWos[59];
   armWWos[158]=1./4.*armWWos[57];
   armWWos[159]=1./4.*armWWos[136];
   armWWos[160]=1./2.*armWWos[137];
   armWWos[161]=armWWos[160] + armWWos[159] - 13./16.*armWWos[33] + 
   armWWos[158] + armWWos[156] + armWWos[39] + 33./16.*armWWos[60] - 5./
   8.*armWWos[54] - 3./8.*armWWos[56] - 5./4. - armWWos[42];
   armWWos[161]=mmt*armWWos[161];
   armWWos[154]=1./4.*armWWos[154] + armWWos[161];
   armWWos[154]=mmt*armWWos[154];
   armWWos[161]=19./8. + armWWos[139];
   armWWos[161]=armWWos[34]*armWWos[161];
   armWWos[142]=1./2.*armWWos[149] + 1./2.*armWWos[142] + 1./3.*
   armWWos[161];
   armWWos[142]=mmH*armWWos[142];
   armWWos[149]=armWWos[155] + 119./3.*armWWos[16];
   armWWos[149]=armWWos[157] + 1./2.*armWWos[149] - 125./3.*armWWos[15]
   ;
   armWWos[149]=armWWos[34]*armWWos[149];
   armWWos[149]=armWWos[154] + 1./8.*armWWos[149] + armWWos[142];
   armWWos[149]=mmt*armWWos[149];
   armWWos[154]=armWWos[176] + armWWos[183] + armWWos[63] + 
   armWWos[182];
   armWWos[157]=1 + armWWos[60];
   armWWos[157]=mmt*armWWos[157];
   armWWos[163]= - armWWos[15]*armWWos[33];
   armWWos[164]=1./3.*armWWos[163];
   armWWos[165]= - 1./2.*armWWos[34];
   armWWos[154]=1./4.*armWWos[157] + armWWos[165] + armWWos[164] + 1./2.
   *armWWos[154] + 1./3.*armWWos[122];
   armWWos[154]=mmt*armWWos[154];
   armWWos[157]= - 7./4.*armWWos[15];
   armWWos[170]=armWWos[16] + armWWos[157];
   armWWos[170]=armWWos[34]*armWWos[170];
   armWWos[154]=1./3.*armWWos[170] + armWWos[154];
   armWWos[154]=armWWos[121]*mmt*armWWos[154];
   armWWos[170]= - armWWos[17] + armWWos[15];
   armWWos[170]=armWWos[34]*armWWos[170];
   armWWos[172]= - mmH*armWWos[34]*armWWos[10];
   armWWos[170]=armWWos[170] + armWWos[172];
   armWWos[170]=mmH*armWWos[170];
   armWWos[172]=1./12.*armWWos[170];
   armWWos[149]=1./2.*armWWos[154] + armWWos[172] + armWWos[149];
   armWWos[149]=armWWos[121]*armWWos[149];
   armWWos[123]=armWWos[149] + 1./6.*armWWos[170] + armWWos[123];
   armWWos[123]=armWWos[121]*armWWos[123];
   armWWos[149]=armWWos[131] + armWWos[130] + 75./4.*armWWos[63] - 
   armWWos[61];
   armWWos[149]=21./2.*armWWos[20] - 395./12.*armWWos[19] + 51*
   armWWos[40] + 3*armWWos[149] - 83./2.*armWWos[38];
   armWWos[125]=3*armWWos[125] + 1./2.*armWWos[126];
   armWWos[125]=mmH*armWWos[125];
   armWWos[126]=3./2.*armWWos[137] + 3./4.*armWWos[136] - 7./16.*
   armWWos[33] + 3./4.*armWWos[57] + 45./16.*armWWos[59] + 3*
   armWWos[39] + 111./16.*armWWos[60] - 15./8.*armWWos[54] - 9./8.*
   armWWos[56] - 1 - armWWos[42];
   armWWos[126]=mmt*armWWos[126];
   armWWos[136]= - 21./2.*armWWos[17] - 17*armWWos[16];
   armWWos[136]=armWWos[33]*armWWos[136];
   armWWos[154]= - 139./2. - 241*armWWos[33];
   armWWos[154]=armWWos[15]*armWWos[154];
   armWWos[173]= - 3*armWWos[33];
   armWWos[179]= - 131./2. + armWWos[173];
   armWWos[179]=armWWos[34]*armWWos[179];
   armWWos[125]=armWWos[126] + 1./4.*armWWos[125] + 1./8.*armWWos[179]
    + 1./24.*armWWos[154] + 1./8.*armWWos[136] - 613./48.*armWWos[16]
    - 27./8.*armWWos[17] + 1./4.*armWWos[149] + 8*armWWos[18];
   armWWos[125]=mmt*armWWos[125];
   armWWos[126]=5*armWWos[17];
   armWWos[136]=3*armWWos[138] + armWWos[126];
   armWWos[136]=armWWos[143] + 1./2.*armWWos[136] + armWWos[129];
   armWWos[138]=3*armWWos[148] - 11./2.*armWWos[10];
   armWWos[138]=1./2.*armWWos[138] + armWWos[150];
   armWWos[138]=mmH*armWWos[138];
   armWWos[136]=1./2.*armWWos[138] + 1./2.*armWWos[136] + armWWos[161];
   armWWos[136]=mmH*armWWos[136];
   armWWos[138]= - 15*armWWos[17] + 293./3.*armWWos[16];
   armWWos[138]= - 39*armWWos[34] + 1./2.*armWWos[138] - 505./3.*
   armWWos[15];
   armWWos[138]=armWWos[34]*armWWos[138];
   armWWos[125]=armWWos[125] + 1./8.*armWWos[138] + armWWos[136];
   armWWos[125]=mmt*armWWos[125];
   armWWos[125]=1./4.*armWWos[170] + armWWos[125];
   armWWos[123]=armWWos[125] + armWWos[123];
   armWWos[123]=armWWos[121]*armWWos[123];
   armWWos[111]=armWWos[111] - 23*armWWos[16];
   armWWos[111]=armWWos[33]*armWWos[111];
   armWWos[136]= - 19./2. - 91./3.*armWWos[33];
   armWWos[136]=armWWos[15]*armWWos[136];
   armWWos[138]= - 53./2. - armWWos[33];
   armWWos[138]=armWWos[34]*armWWos[138];
   armWWos[111]=armWWos[127] + 1./2.*armWWos[138] + 1./2.*armWWos[136]
    + 1./2.*armWWos[111] - 43./4.*armWWos[16] + armWWos[135] + 25./2.*
   armWWos[18] + armWWos[134] - 75./4.*armWWos[19] + armWWos[132] - 13*
   armWWos[38] + armWWos[131] + armWWos[130] + 79./4.*armWWos[63] - 
   armWWos[61];
   armWWos[127]=39./2.*armWWos[60] - 5*armWWos[54] + 1 - 3*armWWos[56];
   armWWos[127]=armWWos[160] + armWWos[159] + 3./16.*armWWos[33] + 
   armWWos[158] + armWWos[156] + 1./8.*armWWos[127] + armWWos[39];
   armWWos[127]=mmt*armWWos[127];
   armWWos[111]=1./4.*armWWos[111] + armWWos[127];
   armWWos[111]=mmt*armWWos[111];
   armWWos[127]=armWWos[155] + 29*armWWos[16];
   armWWos[127]= - 13./2.*armWWos[34] + 1./4.*armWWos[127] - 95./3.*
   armWWos[15];
   armWWos[127]=armWWos[34]*armWWos[127];
   armWWos[111]=armWWos[111] + 1./4.*armWWos[127] + armWWos[142];
   armWWos[111]=mmt*armWWos[111];
   armWWos[111]=armWWos[172] + armWWos[111];
   armWWos[111]=armWWos[177]*armWWos[111];
   armWWos[111]=armWWos[125] + armWWos[111];
   armWWos[111]=armWWos[177]*armWWos[111];
   armWWos[125]=5./4. - armWWos[33];
   armWWos[125]=armWWos[33]*armWWos[125];
   armWWos[127]=armWWos[5]*armWWos[34]*armWWos[33];
   armWWos[125]=armWWos[125] + 3*armWWos[127];
   armWWos[125]=mmt*armWWos[125];
   armWWos[130]=5./2. - armWWos[33];
   armWWos[130]=armWWos[34]*armWWos[130];
   armWWos[125]=armWWos[125] + 1./2.*armWWos[130] + 3*armWWos[106];
   armWWos[125]=mmt*armWWos[125];
   armWWos[125]=1./2.*armWWos[105] + armWWos[125];
   armWWos[125]=mmt*armWWos[125];
   armWWos[131]=armWWos[177]*armWWos[125];
   armWWos[132]=3*armWWos[125];
   armWWos[131]=armWWos[132] + armWWos[131];
   armWWos[131]=armWWos[177]*armWWos[131];
   armWWos[134]=5./2. + armWWos[162];
   armWWos[134]=armWWos[33]*armWWos[134];
   armWWos[127]=armWWos[134] + 6*armWWos[127];
   armWWos[127]=mmt*armWWos[127];
   armWWos[106]=armWWos[127] + armWWos[130] + 6*armWWos[106];
   armWWos[106]=mmt*armWWos[106];
   armWWos[106]=armWWos[105] + armWWos[106];
   armWWos[106]=mmt*armWWos[106];
   armWWos[125]=armWWos[121]*armWWos[125];
   armWWos[106]=armWWos[106] + armWWos[125];
   armWWos[106]=armWWos[121]*armWWos[106];
   armWWos[106]=armWWos[132] + armWWos[106];
   armWWos[106]=armWWos[121]*armWWos[106];
   armWWos[106]=armWWos[131] + armWWos[106];
   armWWos[106]=armWWos[31]*armWWos[106];
   armWWos[106]=armWWos[106] + armWWos[111] + armWWos[123];
   armWWos[106]=armWWos[31]*armWWos[106];
   armWWos[111]= - 3./2.*armWWos[17] + armWWos[183] + armWWos[62] + 
   armWWos[182];
   armWWos[123]=armWWos[59] + 1 + armWWos[49];
   armWWos[125]= - 1./2.*armWWos[14];
   armWWos[123]=armWWos[125] + 1./2.*armWWos[123] + armWWos[57];
   armWWos[127]=1./2.*armWWos[123] + armWWos[151];
   armWWos[127]=mmH*armWWos[127];
   armWWos[127]=armWWos[127] + armWWos[165] + armWWos[164] + 1./2.*
   armWWos[111] + armWWos[145];
   armWWos[127]=mmH*armWWos[127];
   armWWos[130]=1 + armWWos[59];
   armWWos[130]=mmt*mmH*armWWos[130];
   armWWos[127]=armWWos[127] + 1./4.*armWWos[130];
   armWWos[127]=mmt*armWWos[127];
   armWWos[131]=armWWos[17] + armWWos[157];
   armWWos[131]=armWWos[34]*armWWos[131];
   armWWos[132]=7./3.*armWWos[10];
   armWWos[134]= - 1 + armWWos[132];
   armWWos[134]=armWWos[34]*armWWos[134];
   armWWos[134]=armWWos[40] + armWWos[134];
   armWWos[134]=mmH*armWWos[134];
   armWWos[135]=1./3.*armWWos[131] + 1./4.*armWWos[134];
   armWWos[135]=mmH*armWWos[135];
   armWWos[127]=armWWos[135] + armWWos[127];
   armWWos[127]=mmt*armWWos[127];
   armWWos[135]=armWWos[121]*armWWos[127];
   armWWos[135]=armWWos[127] + 1./2.*armWWos[135];
   armWWos[135]=armWWos[121]*armWWos[135];
   armWWos[136]=3./2.*armWWos[123] + armWWos[150];
   armWWos[136]=mmH*armWWos[136];
   armWWos[136]=armWWos[136] + armWWos[194] + armWWos[163] + 3./2.*
   armWWos[111] + armWWos[129];
   armWWos[136]=mmH*armWWos[136];
   armWWos[136]=armWWos[136] + 3./4.*armWWos[130];
   armWWos[136]=mmt*armWWos[136];
   armWWos[138]=7*armWWos[10];
   armWWos[142]= - 3 + armWWos[138];
   armWWos[142]=armWWos[34]*armWWos[142];
   armWWos[142]=armWWos[181] + armWWos[142];
   armWWos[142]=mmH*armWWos[142];
   armWWos[131]=armWWos[131] + 1./4.*armWWos[142];
   armWWos[131]=mmH*armWWos[131];
   armWWos[131]=armWWos[131] + armWWos[136];
   armWWos[131]=mmt*armWWos[131];
   armWWos[131]=1./2.*armWWos[131] + armWWos[135];
   armWWos[131]=armWWos[121]*armWWos[131];
   armWWos[123]=armWWos[123] + 2./3.*armWWos[150];
   armWWos[123]=mmH*armWWos[123];
   armWWos[111]=armWWos[123] - armWWos[34] + 2./3.*armWWos[163] + 
   armWWos[111] + 2./3.*armWWos[129];
   armWWos[111]=mmH*armWWos[111];
   armWWos[111]=armWWos[111] + 1./2.*armWWos[130];
   armWWos[111]=mmt*armWWos[111];
   armWWos[123]=2*armWWos[17];
   armWWos[129]=armWWos[123] - 7./2.*armWWos[15];
   armWWos[129]=armWWos[34]*armWWos[129];
   armWWos[129]=1./3.*armWWos[129] + 1./2.*armWWos[134];
   armWWos[129]=mmH*armWWos[129];
   armWWos[111]=armWWos[129] + armWWos[111];
   armWWos[111]=mmt*armWWos[111];
   armWWos[129]=armWWos[111] + armWWos[131];
   armWWos[129]=armWWos[121]*armWWos[129];
   armWWos[127]=armWWos[177]*armWWos[127];
   armWWos[111]=armWWos[111] + 1./2.*armWWos[127];
   armWWos[111]=armWWos[177]*armWWos[111];
   armWWos[127]= - armWWos[34]*armWWos[33];
   armWWos[130]= - mmt*pow(armWWos[33],2);
   armWWos[131]=armWWos[127] + 1./2.*armWWos[130];
   armWWos[131]=mmt*armWWos[131];
   armWWos[131]= - 1./2.*armWWos[105] + armWWos[131];
   armWWos[134]=pow(mmt,2);
   armWWos[131]=armWWos[134]*armWWos[131];
   armWWos[135]=armWWos[177]*armWWos[131];
   armWWos[127]=2*armWWos[127] + armWWos[130];
   armWWos[127]=mmt*armWWos[127];
   armWWos[105]= - armWWos[105] + armWWos[127];
   armWWos[105]=armWWos[134]*armWWos[105];
   armWWos[127]=2*armWWos[105];
   armWWos[130]=armWWos[127] + armWWos[135];
   armWWos[130]=armWWos[177]*armWWos[130];
   armWWos[134]=armWWos[121]*armWWos[131];
   armWWos[105]=armWWos[105] + armWWos[134];
   armWWos[105]=armWWos[121]*armWWos[105];
   armWWos[105]=3*armWWos[131] + armWWos[105];
   armWWos[105]=armWWos[121]*armWWos[105];
   armWWos[105]=armWWos[127] + armWWos[105];
   armWWos[105]=armWWos[121]*armWWos[105];
   armWWos[105]=armWWos[130] + armWWos[105];
   armWWos[105]=armWWos[31]*armWWos[105];
   armWWos[105]=armWWos[105] + armWWos[111] + armWWos[129];
   armWWos[105]=armWWos[31]*armWWos[105];
   armWWos[111]= - 1./6.*armWWos[80];
   armWWos[127]=1./6.*armWWos[14];
   armWWos[129]=pow(armWWos[10],2);
   armWWos[130]= - 1./9.*armWWos[129] + armWWos[127] + 1 + armWWos[111]
   ;
   armWWos[130]=mmH*armWWos[130];
   armWWos[131]=1./2.*armWWos[67];
   armWWos[134]= - 5./2.*armWWos[20] + armWWos[93] + armWWos[131];
   armWWos[134]=1./2.*armWWos[134] - armWWos[18];
   armWWos[135]= - 1 - 7./18.*armWWos[10];
   armWWos[135]=armWWos[17]*armWWos[135];
   armWWos[136]=1 + 7./9.*armWWos[10];
   armWWos[136]=armWWos[15]*armWWos[136];
   armWWos[130]=1./2.*armWWos[130] + 1./4.*armWWos[136] + 1./3.*
   armWWos[134] + 1./2.*armWWos[135];
   armWWos[130]=mmH*armWWos[130];
   armWWos[135]=19*armWWos[17];
   armWWos[136]=armWWos[135] + 13*armWWos[15];
   armWWos[136]=armWWos[15]*armWWos[136];
   armWWos[142]=pow(armWWos[17],2);
   armWWos[136]=armWWos[142] + 1./4.*armWWos[136];
   armWWos[136]=1./9.*armWWos[136] + armWWos[130];
   armWWos[143]=pow(mmH,2);
   armWWos[136]=armWWos[177]*armWWos[65]*armWWos[143]*armWWos[136];
   armWWos[145]=armWWos[135] + 31*armWWos[15];
   armWWos[145]=armWWos[15]*armWWos[145];
   armWWos[145]=armWWos[142] + 1./4.*armWWos[145];
   armWWos[145]=1./9.*armWWos[145] + armWWos[130];
   armWWos[145]=armWWos[65]*armWWos[143]*armWWos[145];
   armWWos[136]=armWWos[145] + 1./4.*armWWos[136];
   armWWos[136]=armWWos[177]*armWWos[136];
   armWWos[148]=armWWos[135] + 49*armWWos[15];
   armWWos[148]=armWWos[15]*armWWos[148];
   armWWos[148]=armWWos[142] + 1./4.*armWWos[148];
   armWWos[148]=1./9.*armWWos[148] + armWWos[130];
   armWWos[148]=armWWos[65]*armWWos[143]*armWWos[148];
   armWWos[149]=85*armWWos[15];
   armWWos[150]=armWWos[135] + armWWos[149];
   armWWos[150]=armWWos[15]*armWWos[150];
   armWWos[150]=armWWos[142] + 1./4.*armWWos[150];
   armWWos[130]=1./9.*armWWos[150] + armWWos[130];
   armWWos[130]=armWWos[121]*armWWos[65]*armWWos[143]*armWWos[130];
   armWWos[130]=armWWos[148] + 1./2.*armWWos[130];
   armWWos[130]=armWWos[121]*armWWos[130];
   armWWos[148]= - 3 - 7./6.*armWWos[10];
   armWWos[148]=armWWos[17]*armWWos[148];
   armWWos[150]=3 + armWWos[132];
   armWWos[150]=armWWos[15]*armWWos[150];
   armWWos[129]= - 1./3.*armWWos[129] + 1./2.*armWWos[14] + 3 - 1./2.*
   armWWos[80];
   armWWos[129]=mmH*armWWos[129];
   armWWos[129]=1./2.*armWWos[129] + 1./4.*armWWos[150] + armWWos[134]
    + 1./2.*armWWos[148];
   armWWos[129]=mmH*armWWos[129];
   armWWos[134]=armWWos[135] + 37*armWWos[15];
   armWWos[134]=armWWos[15]*armWWos[134];
   armWWos[134]=armWWos[142] + 1./4.*armWWos[134];
   armWWos[129]=1./3.*armWWos[134] + armWWos[129];
   armWWos[129]=armWWos[65]*armWWos[143]*armWWos[129];
   armWWos[129]=1./2.*armWWos[129] + armWWos[130];
   armWWos[129]=armWWos[121]*armWWos[129];
   armWWos[129]=armWWos[145] + 1./2.*armWWos[129];
   armWWos[129]=armWWos[121]*armWWos[129];
   armWWos[105]=armWWos[105] + armWWos[136] + armWWos[129];
   armWWos[105]=armWWos[2]*armWWos[105];
   armWWos[129]= - 133./6. + 3*armWWos[84];
   armWWos[130]=5./3.*armWWos[86];
   armWWos[129]=1./2.*armWWos[129] + armWWos[130];
   armWWos[134]= - 17./6.*armWWos[14];
   armWWos[135]=13*armWWos[10];
   armWWos[136]= - 1 + armWWos[135];
   armWWos[136]=armWWos[10]*armWWos[136];
   armWWos[143]= - armWWos[73] - 1./3.*armWWos[70];
   armWWos[143]=mmH*armWWos[143];
   armWWos[129]=1./2.*armWWos[143] + 1./18.*armWWos[136] + armWWos[134]
    - 17./3.*armWWos[89] + 1./4.*armWWos[85] + 1./2.*armWWos[129] + 
   armWWos[80];
   armWWos[129]=mmH*armWWos[129];
   armWWos[136]=25*armWWos[10];
   armWWos[145]=167 + armWWos[136];
   armWWos[145]=armWWos[17]*armWWos[145];
   armWWos[138]=17./3. + armWWos[138];
   armWWos[138]=armWWos[16]*armWWos[138];
   armWWos[148]=79./3. + 41*armWWos[10];
   armWWos[148]=armWWos[15]*armWWos[148];
   armWWos[150]=5*armWWos[69] - armWWos[93];
   armWWos[129]=1./2.*armWWos[129] + 1./12.*armWWos[148] + 1./2.*
   armWWos[138] + 1./36.*armWWos[145] + 11./4.*armWWos[18] - 23./12.*
   armWWos[20] + 17./6.*armWWos[19] + 37./24.*armWWos[68] - 1./8.*
   armWWos[95] + 19./12.*armWWos[67] + 1./8.*armWWos[150] - 17./3.*
   armWWos[94];
   armWWos[129]=mmH*armWWos[129];
   armWWos[138]=armWWos[17] - 17./8.*armWWos[16];
   armWWos[138]=armWWos[16]*armWWos[138];
   armWWos[145]=5./8.*armWWos[142];
   armWWos[138]=armWWos[145] + armWWos[138];
   armWWos[148]=43./9.*armWWos[17] + 23*armWWos[16];
   armWWos[148]=1./2.*armWWos[148] - 119./9.*armWWos[15];
   armWWos[148]=armWWos[15]*armWWos[148];
   armWWos[129]=1./2.*armWWos[129] + 1./3.*armWWos[138] + 1./4.*
   armWWos[148];
   armWWos[129]=armWWos[177]*armWWos[65]*mmH*armWWos[129];
   armWWos[138]=97 + armWWos[136];
   armWWos[138]=armWWos[17]*armWWos[138];
   armWWos[148]=19 + 131./18.*armWWos[10];
   armWWos[148]=armWWos[16]*armWWos[148];
   armWWos[138]=armWWos[148] + 1./6.*armWWos[138] - 5./6.*armWWos[18]
    - 23./2.*armWWos[20] + 53./6.*armWWos[19] + 55./12.*armWWos[68] - 3.
   /4.*armWWos[95] + 31./3.*armWWos[67] + 3./4.*armWWos[150] - 71./3.*
   armWWos[94];
   armWWos[148]=7 + armWWos[135];
   armWWos[148]=armWWos[10]*armWWos[148];
   armWWos[151]= - 3*armWWos[73] - armWWos[70];
   armWWos[151]=mmH*armWWos[151];
   armWWos[148]=1./8.*armWWos[151] + 1./24.*armWWos[148] - 17./8.*
   armWWos[14] - 71./24.*armWWos[89] + 3./16.*armWWos[85] + 3./4.*
   armWWos[80] + 5./8.*armWWos[86] + 9./16.*armWWos[84] - 81./32. + 1./
   3.*armWWos[90];
   armWWos[148]=mmH*armWWos[148];
   armWWos[151]=1 + 11./3.*armWWos[10];
   armWWos[151]=armWWos[15]*armWWos[151];
   armWWos[138]=1./2.*armWWos[148] + 1./8.*armWWos[138] + 2./3.*
   armWWos[151];
   armWWos[138]=mmH*armWWos[138];
   armWWos[148]=5./4.*armWWos[142];
   armWWos[151]= - 7./3.*armWWos[17];
   armWWos[154]=armWWos[151] - 23./4.*armWWos[16];
   armWWos[154]=armWWos[16]*armWWos[154];
   armWWos[154]=armWWos[148] + 1./3.*armWWos[154];
   armWWos[155]=7*armWWos[17];
   armWWos[156]=armWWos[155] + 23./2.*armWWos[16];
   armWWos[156]=5./2.*armWWos[156] - 58*armWWos[15];
   armWWos[156]=armWWos[15]*armWWos[156];
   armWWos[138]=armWWos[138] + 1./4.*armWWos[154] + 1./9.*armWWos[156];
   armWWos[138]=armWWos[65]*mmH*armWWos[138];
   armWWos[129]=armWWos[138] + 1./2.*armWWos[129];
   armWWos[129]=armWWos[177]*armWWos[129];
   armWWos[154]= - 25./32. + armWWos[90];
   armWWos[156]=23 + armWWos[135];
   armWWos[156]=armWWos[10]*armWWos[156];
   armWWos[154]=1./8.*armWWos[143] + 1./72.*armWWos[156] - 17./24.*
   armWWos[14] - 5./8.*armWWos[89] + 1./16.*armWWos[85] + 1./4.*
   armWWos[80] + 5./24.*armWWos[86] + 1./3.*armWWos[154] + 3./16.*
   armWWos[84];
   armWWos[154]=mmH*armWWos[154];
   armWWos[156]=1./4.*armWWos[150] - 5*armWWos[94];
   armWWos[136]=29 + armWWos[136];
   armWWos[136]=armWWos[17]*armWWos[136];
   armWWos[157]=9 - 73./18.*armWWos[10];
   armWWos[157]=armWWos[16]*armWWos[157];
   armWWos[158]= - 13 + 53*armWWos[10];
   armWWos[158]=armWWos[15]*armWWos[158];
   armWWos[136]=armWWos[154] + 1./36.*armWWos[158] + 1./4.*armWWos[157]
    + 1./72.*armWWos[136] - 47./24.*armWWos[18] - 23./24.*armWWos[20]
    - 11./24.*armWWos[19] + 5./48.*armWWos[68] - 1./16.*armWWos[95] + 1.
   /4.*armWWos[156] + armWWos[67];
   armWWos[136]=mmH*armWWos[136];
   armWWos[154]= - 5./3.*armWWos[17] - 1./8.*armWWos[16];
   armWWos[154]=armWWos[16]*armWWos[154];
   armWWos[156]=61*armWWos[17] + 59*armWWos[16];
   armWWos[156]=1./2.*armWWos[156] - 113*armWWos[15];
   armWWos[156]=armWWos[15]*armWWos[156];
   armWWos[145]=1./6.*armWWos[156] + armWWos[145] + armWWos[154];
   armWWos[136]=1./3.*armWWos[145] + armWWos[136];
   armWWos[136]=armWWos[65]*mmH*armWWos[136];
   armWWos[145]= - 1./2.*armWWos[19];
   armWWos[131]=armWWos[17] + armWWos[145] - armWWos[94] + armWWos[131]
   ;
   armWWos[154]=1 - 7./12.*armWWos[10];
   armWWos[154]=armWWos[16]*armWWos[154];
   armWWos[156]=1 + armWWos[132];
   armWWos[156]=armWWos[15]*armWWos[156];
   armWWos[157]= - 1 - armWWos[89];
   armWWos[157]=mmH*armWWos[157];
   armWWos[131]=1./4.*armWWos[157] + 1./4.*armWWos[156] + 1./2.*
   armWWos[131] + armWWos[154];
   armWWos[131]=mmH*armWWos[131];
   armWWos[154]= - armWWos[16]*armWWos[17];
   armWWos[156]=armWWos[17] + armWWos[16];
   armWWos[156]=7./2.*armWWos[156] - 5*armWWos[15];
   armWWos[156]=armWWos[15]*armWWos[156];
   armWWos[154]=armWWos[154] + 1./2.*armWWos[156];
   armWWos[131]=1./3.*armWWos[154] + armWWos[131];
   armWWos[131]=armWWos[121]*armWWos[65]*mmH*armWWos[131];
   armWWos[131]=armWWos[136] + 1./6.*armWWos[131];
   armWWos[131]=armWWos[121]*armWWos[131];
   armWWos[136]= - 1./2.*armWWos[95];
   armWWos[154]=3./2.*armWWos[68];
   armWWos[150]= - 19./3.*armWWos[18] - 23./3.*armWWos[20] + 19./6.*
   armWWos[19] + armWWos[154] + armWWos[136] + 43./6.*armWWos[67] + 1./
   2.*armWWos[150] - 37./3.*armWWos[94];
   armWWos[156]=31 + 25./2.*armWWos[10];
   armWWos[156]=armWWos[17]*armWWos[156];
   armWWos[150]=1./2.*armWWos[150] + 1./9.*armWWos[156];
   armWWos[156]= - 55./16. + armWWos[90];
   armWWos[135]=11 + armWWos[135];
   armWWos[135]=armWWos[10]*armWWos[135];
   armWWos[135]=1./4.*armWWos[143] + 1./36.*armWWos[135] - 17./12.*
   armWWos[14] - 37./24.*armWWos[89] + 1./8.*armWWos[85] + 1./2.*
   armWWos[80] + 5./12.*armWWos[86] + 1./3.*armWWos[156] + 3./8.*
   armWWos[84];
   armWWos[135]=mmH*armWWos[135];
   armWWos[143]=1 + 1./48.*armWWos[10];
   armWWos[143]=armWWos[16]*armWWos[143];
   armWWos[156]=17 + 229*armWWos[10];
   armWWos[156]=armWWos[15]*armWWos[156];
   armWWos[135]=1./2.*armWWos[135] + 1./144.*armWWos[156] + 1./4.*
   armWWos[150] + 5./3.*armWWos[143];
   armWWos[135]=mmH*armWWos[135];
   armWWos[143]= - 13./9.*armWWos[17] - 1./2.*armWWos[16];
   armWWos[143]=armWWos[16]*armWWos[143];
   armWWos[150]=79*armWWos[17] + 253./3.*armWWos[16];
   armWWos[150]=1./2.*armWWos[150] - 115*armWWos[15];
   armWWos[150]=armWWos[15]*armWWos[150];
   armWWos[143]=1./6.*armWWos[150] + 5./6.*armWWos[142] + armWWos[143];
   armWWos[135]=1./4.*armWWos[143] + armWWos[135];
   armWWos[135]=armWWos[65]*mmH*armWWos[135];
   armWWos[131]=armWWos[135] + 1./2.*armWWos[131];
   armWWos[131]=armWWos[121]*armWWos[131];
   armWWos[131]=armWWos[138] + armWWos[131];
   armWWos[131]=armWWos[121]*armWWos[131];
   armWWos[105]=armWWos[105] + armWWos[106] + armWWos[129] + 
   armWWos[131];
   armWWos[105]=armWWos[2]*armWWos[105];
   armWWos[106]= - 2*armWWos[90];
   armWWos[129]=439./192. + armWWos[106];
   armWWos[131]= - 1./3.*armWWos[80];
   armWWos[135]= - 43./48.*armWWos[91];
   armWWos[138]= - 2./3.*armWWos[30];
   armWWos[143]=pow(armWWos[9],2);
   armWWos[150]=1./24.*armWWos[143];
   armWWos[156]= - 11./3.*armWWos[10];
   armWWos[157]=armWWos[156] - 131./9. + 13./2.*armWWos[9];
   armWWos[157]=armWWos[10]*armWWos[157];
   armWWos[158]=2*armWWos[78];
   armWWos[159]=armWWos[158] - armWWos[75];
   armWWos[160]=1./2.*armWWos[71] - 1./4.*armWWos[74] - 1./8.*
   armWWos[96] + armWWos[159] - 1./2.*armWWos[70];
   armWWos[160]=mmH*armWWos[160];
   armWWos[129]=1./3.*armWWos[160] + 1./12.*armWWos[157] + armWWos[150]
    + 61./48.*armWWos[14] - 5./48.*armWWos[13] + armWWos[138] + 137./48.
   *armWWos[89] - 7./48.*armWWos[85] + armWWos[135] + armWWos[131] - 1./
   8.*armWWos[87] + 5./16.*armWWos[81] + 5./16.*armWWos[82] - 1./6.*
   armWWos[88] + 1./12.*armWWos[83] + 1./3.*armWWos[129] - 7./16.*
   armWWos[84];
   armWWos[129]=mmH*armWWos[129];
   armWWos[157]=3./2.*armWWos[69] + armWWos[93];
   armWWos[160]= - 17./6.*armWWos[68] + armWWos[95] - 19./3.*
   armWWos[67] - 11./12.*armWWos[66] + armWWos[157] + 89./6.*
   armWWos[94];
   armWWos[161]= - 5./9. + 1./8.*armWWos[9];
   armWWos[161]=35*armWWos[161] - 23./4.*armWWos[10];
   armWWos[161]=armWWos[17]*armWWos[161];
   armWWos[162]= - 1./2.*armWWos[9];
   armWWos[163]=89./6.*armWWos[10] - 205./3. + armWWos[162];
   armWWos[163]=armWWos[16]*armWWos[163];
   armWWos[164]= - 11*armWWos[9];
   armWWos[170]= - 457./27. + armWWos[164];
   armWWos[170]=1./8.*armWWos[170] - 43./9.*armWWos[10];
   armWWos[170]=armWWos[15]*armWWos[170];
   armWWos[129]=armWWos[129] + 1./2.*armWWos[170] + 1./12.*armWWos[163]
    + 1./6.*armWWos[161] + 3./4.*armWWos[18] - 3./16.*armWWos[20] + 1./
   4.*armWWos[160] + 2./3.*armWWos[19];
   armWWos[129]=mmH*armWWos[129];
   armWWos[160]= - 1./2.*armWWos[142];
   armWWos[161]=7./3.*armWWos[17] + 31./2.*armWWos[16];
   armWWos[161]=armWWos[16]*armWWos[161];
   armWWos[161]=armWWos[160] + armWWos[161];
   armWWos[163]=1909./8.*armWWos[15] + 13*armWWos[17] + 2969./24.*
   armWWos[16];
   armWWos[163]=armWWos[15]*armWWos[163];
   armWWos[161]=1./2.*armWWos[161] + 1./3.*armWWos[163];
   armWWos[129]=1./2.*armWWos[161] + armWWos[129];
   armWWos[129]=armWWos[65]*armWWos[129];
   armWWos[161]=1./2.*armWWos[86];
   armWWos[163]=armWWos[161] - 233./48. + armWWos[88];
   armWWos[170]= - 1./4.*armWWos[85];
   armWWos[163]=11./12.*armWWos[13] - 77./36.*armWWos[89] + 
   armWWos[170] - 13./4.*armWWos[91] - 5./4.*armWWos[81] + 1./3.*
   armWWos[163] - 5./4.*armWWos[82];
   armWWos[172]=1./3.*armWWos[14];
   armWWos[179]= - 7*armWWos[9];
   armWWos[180]=11./9. + armWWos[179];
   armWWos[180]=armWWos[10]*armWWos[180];
   armWWos[184]= - mmH*armWWos[71];
   armWWos[163]=1./12.*armWWos[184] + 1./16.*armWWos[180] + 1./4.*
   armWWos[163] + armWWos[172];
   armWWos[163]=mmH*armWWos[163];
   armWWos[180]=11./16. - 2./3.*armWWos[9];
   armWWos[180]=armWWos[17]*armWWos[180];
   armWWos[132]=armWWos[132] + 1./12. - armWWos[9];
   armWWos[132]=armWWos[16]*armWWos[132];
   armWWos[184]=3*armWWos[9];
   armWWos[185]= - 5./18.*armWWos[10] + 1./9. + armWWos[184];
   armWWos[185]=armWWos[15]*armWWos[185];
   armWWos[132]=armWWos[163] + 1./4.*armWWos[185] + 1./12.*armWWos[132]
    + armWWos[180] - 4./9.*armWWos[18] - 1./8.*armWWos[20] + 1./18.*
   armWWos[19] - 1./6.*armWWos[68] + 1./12.*armWWos[95] - 1./48.*
   armWWos[67] - 2./9.*armWWos[94] + 5./16.*armWWos[66];
   armWWos[132]=mmH*armWWos[132];
   armWWos[163]=1./9.*armWWos[17] + 5./8.*armWWos[16];
   armWWos[163]=armWWos[16]*armWWos[163];
   armWWos[126]=631*armWWos[15] + armWWos[126] + 217*armWWos[16];
   armWWos[126]=armWWos[15]*armWWos[126];
   armWWos[126]=armWWos[132] + armWWos[163] + 1./72.*armWWos[126];
   armWWos[126]=armWWos[65]*armWWos[126];
   armWWos[132]=1./2.*armWWos[66];
   armWWos[163]= - armWWos[94] + armWWos[132];
   armWWos[180]= - 1./2.*armWWos[20];
   armWWos[185]=armWWos[163] + armWWos[180];
   armWWos[186]=1 - 7./12.*armWWos[9];
   armWWos[186]=armWWos[17]*armWWos[186];
   armWWos[189]=1./2.*armWWos[16];
   armWWos[119]=1 + armWWos[119];
   armWWos[119]=armWWos[15]*armWWos[119];
   armWWos[119]=1./4.*armWWos[119] + armWWos[189] + 1./2.*armWWos[185]
    + armWWos[186];
   armWWos[185]= - 1./2.*armWWos[82];
   armWWos[186]= - 1./2.*armWWos[81];
   armWWos[190]=armWWos[186] - 1 + armWWos[185];
   armWWos[191]=1./6.*armWWos[13];
   armWWos[127]=armWWos[127] + armWWos[191] - 1./2.*armWWos[89] + 1./3.
   *armWWos[190] - 1./2.*armWWos[91];
   armWWos[190]= - armWWos[10]*armWWos[9];
   armWWos[127]=1./2.*armWWos[127] + 1./9.*armWWos[190];
   armWWos[127]=mmH*armWWos[127];
   armWWos[119]=1./3.*armWWos[119] + armWWos[127];
   armWWos[119]=mmH*armWWos[119];
   armWWos[127]=19*armWWos[16] + armWWos[149];
   armWWos[127]=armWWos[15]*armWWos[127];
   armWWos[149]=pow(armWWos[16],2);
   armWWos[127]=armWWos[149] + 1./4.*armWWos[127];
   armWWos[119]=1./9.*armWWos[127] + armWWos[119];
   armWWos[119]=armWWos[121]*armWWos[65]*armWWos[119];
   armWWos[119]=armWWos[126] + 1./4.*armWWos[119];
   armWWos[119]=armWWos[121]*armWWos[119];
   armWWos[126]= - 5./3. - armWWos[6];
   armWWos[127]=armWWos[1]*armWWos[126];
   armWWos[126]=armWWos[35]*armWWos[126];
   armWWos[126]=1./3.*armWWos[127] + armWWos[126];
   armWWos[126]=armWWos[10]*armWWos[126];
   armWWos[127]=armWWos[1]*armWWos[12];
   armWWos[190]=armWWos[35]*armWWos[12];
   armWWos[126]=armWWos[126] + 1./3.*armWWos[127] + armWWos[190];
   armWWos[126]=mmH*armWWos[126];
   armWWos[127]= - 1./2.*armWWos[6];
   armWWos[190]= - 1./3. + armWWos[127];
   armWWos[192]=armWWos[1]*armWWos[190];
   armWWos[190]=armWWos[35]*armWWos[190];
   armWWos[190]=1./3.*armWWos[192] + armWWos[190];
   armWWos[190]=armWWos[17]*armWWos[190];
   armWWos[192]=5./3. + armWWos[6];
   armWWos[193]=armWWos[1]*armWWos[192];
   armWWos[192]=armWWos[35]*armWWos[192];
   armWWos[192]=1./3.*armWWos[193] + armWWos[192];
   armWWos[192]=armWWos[15]*armWWos[192];
   armWWos[126]=1./2.*armWWos[126] + armWWos[190] + 1./2.*armWWos[192];
   armWWos[126]=mmH*armWWos[126];
   armWWos[190]=1./6.*armWWos[126];
   armWWos[119]=armWWos[119] + armWWos[190] + armWWos[129];
   armWWos[119]=armWWos[121]*armWWos[119];
   armWWos[129]= - 1./4.*armWWos[96];
   armWWos[192]= - 1./2.*armWWos[74];
   armWWos[159]=1./4.*armWWos[71] + armWWos[192] + armWWos[129] + 
   armWWos[159] - armWWos[70];
   armWWos[159]=mmH*armWWos[159];
   armWWos[193]=19./9. + 17*armWWos[9];
   armWWos[193]=5./8.*armWWos[193] + armWWos[156];
   armWWos[193]=armWWos[10]*armWWos[193];
   armWWos[106]=3023./192. + armWWos[106];
   armWWos[106]=1./3.*armWWos[159] + 1./6.*armWWos[193] + 1./12.*
   armWWos[143] + 13./8.*armWWos[14] - 41./48.*armWWos[13] + 
   armWWos[138] + 143./16.*armWWos[89] - 11./48.*armWWos[85] + 
   armWWos[135] - 2./3.*armWWos[80] - 1./12.*armWWos[87] + 41./48.*
   armWWos[81] + 41./48.*armWWos[82] - 1./24.*armWWos[86] - 1./12.*
   armWWos[88] + 1./6.*armWWos[83] + 1./3.*armWWos[106] - 7./8.*
   armWWos[84];
   armWWos[106]=mmH*armWWos[106];
   armWWos[135]= - 7./2.*armWWos[68] + 5./6.*armWWos[95] - 119./24.*
   armWWos[67] - 41./24.*armWWos[66] + armWWos[157] + 43./2.*
   armWWos[94];
   armWWos[159]= - 5689./18. + 67*armWWos[9];
   armWWos[159]=1./2.*armWWos[159] - 19*armWWos[10];
   armWWos[159]=armWWos[17]*armWWos[159];
   armWWos[193]= - 1./3.*armWWos[10] - 1765./24. - armWWos[9];
   armWWos[193]=armWWos[16]*armWWos[193];
   armWWos[194]= - 113./3.*armWWos[10] - 1327./27. - 21*armWWos[9];
   armWWos[194]=armWWos[15]*armWWos[194];
   armWWos[106]=armWWos[106] + 1./8.*armWWos[194] + 1./6.*armWWos[193]
    + 1./12.*armWWos[159] - 11./3.*armWWos[18] - 5./6.*armWWos[20] + 1./
   2.*armWWos[135] - 5./3.*armWWos[19];
   armWWos[106]=mmH*armWWos[106];
   armWWos[135]=armWWos[17] - 677./8.*armWWos[16];
   armWWos[135]=armWWos[16]*armWWos[135];
   armWWos[135]=armWWos[148] + armWWos[135];
   armWWos[148]=73./8.*armWWos[17];
   armWWos[159]= - 499./6.*armWWos[15];
   armWWos[193]=armWWos[159] + armWWos[148] - 337./3.*armWWos[16];
   armWWos[193]=armWWos[15]*armWWos[193];
   armWWos[193]=armWWos[135] + armWWos[193];
   armWWos[193]=1./3.*armWWos[193] + armWWos[106];
   armWWos[193]=armWWos[65]*armWWos[193];
   armWWos[126]=1./3.*armWWos[126];
   armWWos[119]=armWWos[119] + armWWos[126] + armWWos[193];
   armWWos[119]=armWWos[121]*armWWos[119];
   armWWos[148]=armWWos[159] + armWWos[148] + 215./3.*armWWos[16];
   armWWos[148]=armWWos[15]*armWWos[148];
   armWWos[135]=armWWos[135] + armWWos[148];
   armWWos[106]=1./3.*armWWos[135] + armWWos[106];
   armWWos[106]=armWWos[65]*armWWos[106];
   armWWos[135]=619./9. + 59*armWWos[9];
   armWWos[135]=1./4.*armWWos[135] + armWWos[156];
   armWWos[135]=armWWos[10]*armWWos[135];
   armWWos[148]= - 1./2.*armWWos[71] + armWWos[192] - armWWos[70] + 
   armWWos[129];
   armWWos[148]=mmH*armWWos[148];
   armWWos[156]=619./9. - 7*armWWos[84];
   armWWos[156]=13./6.*armWWos[87] + 13./6.*armWWos[81] + 13./6.*
   armWWos[82] - 1./6.*armWWos[86] + 1./3.*armWWos[88] + 1./4.*
   armWWos[156] + 1./3.*armWWos[83];
   armWWos[131]=1./6.*armWWos[148] + 1./12.*armWWos[135] + armWWos[150]
    + 17./48.*armWWos[14] - 5./4.*armWWos[13] + 85./12.*armWWos[89] - 1.
   /12.*armWWos[85] + 1./4.*armWWos[156] + armWWos[131];
   armWWos[131]=mmH*armWWos[131];
   armWWos[135]= - 43./12.*armWWos[67] - 5./2.*armWWos[66] + 
   armWWos[157] + 169./6.*armWWos[94];
   armWWos[135]= - 25./12.*armWWos[68] + 1./2.*armWWos[135] + 1./3.*
   armWWos[95];
   armWWos[148]= - 4505./27. + 41*armWWos[9];
   armWWos[148]=1./2.*armWWos[148] - 5*armWWos[10];
   armWWos[148]=armWWos[17]*armWWos[148];
   armWWos[150]= - 5*armWWos[9];
   armWWos[156]= - 21*armWWos[10] - 227./6. + armWWos[150];
   armWWos[156]=armWWos[16]*armWWos[156];
   armWWos[157]= - 31*armWWos[9];
   armWWos[159]= - 2197./27. + armWWos[157];
   armWWos[159]=1./2.*armWWos[159] - 167./9.*armWWos[10];
   armWWos[159]=armWWos[15]*armWWos[159];
   armWWos[131]=armWWos[131] + 1./8.*armWWos[159] + 1./8.*armWWos[156]
    + 1./8.*armWWos[148] - 53./12.*armWWos[18] + 17./48.*armWWos[20] + 
   1./2.*armWWos[135] - 14./3.*armWWos[19];
   armWWos[131]=mmH*armWWos[131];
   armWWos[135]= - armWWos[17] - 219*armWWos[16];
   armWWos[135]=armWWos[16]*armWWos[135];
   armWWos[148]= - 6119./18.*armWWos[15] + armWWos[155] - 6563./6.*
   armWWos[16];
   armWWos[148]=armWWos[15]*armWWos[148];
   armWWos[135]=1./2.*armWWos[148] + 13./6.*armWWos[142] + armWWos[135]
   ;
   armWWos[131]=1./4.*armWWos[135] + armWWos[131];
   armWWos[131]=armWWos[65]*armWWos[131];
   armWWos[131]=armWWos[190] + armWWos[131];
   armWWos[131]=armWWos[177]*armWWos[131];
   armWWos[106]=armWWos[131] + armWWos[126] + armWWos[106];
   armWWos[106]=armWWos[177]*armWWos[106];
   armWWos[101]=armWWos[105] + armWWos[101] + armWWos[106] + 
   armWWos[119];
   armWWos[101]=armWWos[2]*armWWos[101];
   armWWos[105]=33*armWWos[9];
   armWWos[106]=245./3. + armWWos[105];
   armWWos[119]=4./3. - armWWos[6];
   armWWos[119]=2*armWWos[1]*armWWos[119];
   armWWos[126]=4 - 3*armWWos[6];
   armWWos[126]=2*armWWos[35]*armWWos[126];
   armWWos[131]= - 3*armWWos[10];
   armWWos[106]=armWWos[131] + armWWos[126] + 1./2.*armWWos[106] + 
   armWWos[119];
   armWWos[106]=armWWos[34]*armWWos[106];
   armWWos[135]=armWWos[175] + armWWos[103];
   armWWos[135]=armWWos[33]*armWWos[135];
   armWWos[148]= - 1./4.*armWWos[20];
   armWWos[155]=1./4.*armWWos[17];
   armWWos[156]=6 + armWWos[178];
   armWWos[156]=3*armWWos[15]*armWWos[156];
   armWWos[159]= - armWWos[16] - armWWos[15];
   armWWos[175]=armWWos[5]*armWWos[34]*armWWos[159];
   armWWos[106]=36*armWWos[175] + armWWos[106] + armWWos[156] + 1./2.*
   armWWos[135] + 10*armWWos[16] + armWWos[155] + armWWos[174] + 
   armWWos[148] + 11*armWWos[19] + 25*armWWos[40] - 11*armWWos[36] - 3*
   armWWos[38];
   armWWos[106]=armWWos[5]*armWWos[106];
   armWWos[135]= - armWWos[61] + 1./2.*armWWos[37];
   armWWos[175]=7./2.*armWWos[18];
   armWWos[135]= - 9./4.*armWWos[17] + armWWos[175] - 7./4.*armWWos[20]
    + armWWos[181] + armWWos[182] + 3*armWWos[135] + armWWos[62];
   armWWos[135]=armWWos[100]*armWWos[135];
   armWWos[181]=3*armWWos[57];
   armWWos[152]=5./2.*armWWos[14] + armWWos[181] + 1./2.*armWWos[59] + 
   armWWos[146] + 43./8. + armWWos[152];
   armWWos[152]=armWWos[100]*armWWos[152];
   armWWos[182]=3*armWWos[10];
   armWWos[159]=armWWos[5]*armWWos[159];
   armWWos[159]=18*armWWos[159] - 5./2.*armWWos[33] - 32 + armWWos[182]
   ;
   armWWos[159]=armWWos[5]*armWWos[159];
   armWWos[190]= - 4*armWWos[46] + 11*armWWos[45];
   armWWos[193]= - 1 + armWWos[10];
   armWWos[194]=armWWos[100]*armWWos[193];
   armWWos[196]=1./2.*armWWos[33]*armWWos[194];
   armWWos[159]=armWWos[159] + armWWos[196] + armWWos[152] + 
   armWWos[44] + 11./12.*armWWos[48] + 2./3.*armWWos[190] + 3./2.*
   armWWos[43];
   armWWos[159]=mmt*armWWos[159];
   armWWos[190]= - 101./8.*armWWos[50];
   armWWos[197]= - 49./8.*armWWos[56];
   armWWos[198]= - 349./128.*armWWos[55];
   armWWos[199]= - 581./16.*armWWos[52] + armWWos[198] + armWWos[197]
    + armWWos[190] - 13./6.*armWWos[42] + 9259./128. + 68./3.*
   armWWos[53];
   armWWos[200]=11351./96. - 67*armWWos[9];
   armWWos[168]=1./4.*armWWos[200] + armWWos[168];
   armWWos[168]= - 17./18.*armWWos[33] + armWWos[171] + armWWos[10] + 1.
   /3.*armWWos[168] + armWWos[167];
   armWWos[168]=armWWos[33]*armWWos[168];
   armWWos[200]=1./4.*armWWos[54];
   armWWos[201]= - 13./12.*armWWos[27];
   armWWos[202]= - 169./192.*armWWos[58];
   armWWos[203]= - 1651./384.*armWWos[60];
   armWWos[204]= - 19./12.*armWWos[30];
   armWWos[205]=1037./384.*armWWos[39];
   armWWos[206]=13./24.*armWWos[59];
   armWWos[207]= - 7./3. + armWWos[6];
   armWWos[208]=1./6.*armWWos[1]*armWWos[207];
   armWWos[207]=1./2.*armWWos[35]*armWWos[207];
   armWWos[209]= - 3./2.*armWWos[10];
   armWWos[133]=7*armWWos[100] + armWWos[133];
   armWWos[133]=1./2.*armWWos[15]*armWWos[133];
   armWWos[210]=4 + armWWos[10];
   armWWos[210]=armWWos[34]*armWWos[100]*armWWos[210];
   armWWos[199]=armWWos[159] + armWWos[137] + armWWos[106] + 
   armWWos[210] + armWWos[133] + armWWos[168] + armWWos[135] + 
   armWWos[209] + armWWos[207] + armWWos[208] - 905./48.*armWWos[9] + 
   armWWos[147] + 835./48.*armWWos[13] + armWWos[128] + armWWos[206] + 
   armWWos[205] + armWWos[204] + armWWos[146] + armWWos[203] + 
   armWWos[202] + armWWos[201] + armWWos[51] + 1./3.*armWWos[199] + 
   armWWos[200];
   armWWos[199]=mmt*armWWos[199];
   armWWos[119]=armWWos[131] + armWWos[126] + armWWos[119] - 16./3. + 
   99./2.*armWWos[9];
   armWWos[119]=armWWos[34]*armWWos[119];
   armWWos[126]= - 1./2.*armWWos[36];
   armWWos[211]=2*armWWos[40];
   armWWos[212]=1./2.*armWWos[19] + armWWos[211] + armWWos[126] - 
   armWWos[38];
   armWWos[103]= - armWWos[17] + armWWos[103];
   armWWos[103]=armWWos[33]*armWWos[103];
   armWWos[213]= - armWWos[16] - 2*armWWos[15];
   armWWos[214]=armWWos[5]*armWWos[34]*armWWos[213];
   armWWos[103]=18*armWWos[214] + armWWos[119] + armWWos[156] + 1./4.*
   armWWos[103] + 9*armWWos[16] + armWWos[155] + armWWos[174] + 3*
   armWWos[212] + armWWos[148];
   armWWos[103]=armWWos[5]*armWWos[103];
   armWWos[119]= - 2479./144. - 81*armWWos[9];
   armWWos[119]=1./2.*armWWos[33] + armWWos[171] + armWWos[10] + 
   armWWos[167] + 1./8.*armWWos[119] + armWWos[169];
   armWWos[119]=armWWos[33]*armWWos[119];
   armWWos[148]=armWWos[5]*armWWos[213];
   armWWos[148]=3*armWWos[148] + armWWos[178] - 3 + armWWos[10];
   armWWos[148]=armWWos[5]*armWWos[148];
   armWWos[155]=armWWos[45] + armWWos[43];
   armWWos[156]=armWWos[155] + 1./2.*armWWos[48];
   armWWos[148]=3*armWWos[148] + armWWos[196] + armWWos[152] + 3./2.*
   armWWos[156] + armWWos[44];
   armWWos[148]=mmt*armWWos[148];
   armWWos[103]=armWWos[148] + armWWos[137] + armWWos[103] + 
   armWWos[210] + armWWos[133] + armWWos[119] + armWWos[135] + 
   armWWos[209] + armWWos[207] + armWWos[208] - 57./16.*armWWos[9] + 
   armWWos[147] + 93./16.*armWWos[13] + armWWos[128] + armWWos[206] + 
   185./128.*armWWos[39] - armWWos[30] + armWWos[146] - 1663./128.*
   armWWos[60] - 25./64.*armWWos[58] - 7./4.*armWWos[27] + armWWos[51]
    + armWWos[200] + 31./16.*armWWos[52] - 41./128.*armWWos[55] - 11./8.
   *armWWos[56] + 28439./1152. - 6*armWWos[50];
   armWWos[103]=mmt*armWWos[103];
   armWWos[119]=armWWos[32] + 51./2.*armWWos[97];
   armWWos[148]=armWWos[119] + 33./2.*armWWos[98];
   armWWos[152]=2*armWWos[100];
   armWWos[148]=1./2.*armWWos[148] + armWWos[152];
   armWWos[148]=armWWos[15]*armWWos[148];
   armWWos[156]= - 9163./72. - 453*armWWos[9];
   armWWos[167]=9./8.*armWWos[98] + armWWos[32] - 3./2.*armWWos[97];
   armWWos[167]=armWWos[16]*armWWos[167];
   armWWos[169]=123./8.*armWWos[98] - armWWos[32] + 69./2.*armWWos[97];
   armWWos[169]=armWWos[34]*armWWos[169];
   armWWos[171]= - 5./3. + 2*armWWos[6];
   armWWos[174]=armWWos[1]*armWWos[171];
   armWWos[171]=armWWos[35]*armWWos[171];
   armWWos[178]= - armWWos[100]*armWWos[17];
   armWWos[148]=1./4.*armWWos[169] + armWWos[148] - 319./128.*
   armWWos[33] + armWWos[178] + 1./4.*armWWos[167] + armWWos[10] + 
   armWWos[171] + 1./16.*armWWos[156] + 1./3.*armWWos[174];
   armWWos[148]=armWWos[34]*armWWos[148];
   armWWos[156]=9*armWWos[17];
   armWWos[167]=armWWos[156] + 719./32.*armWWos[16];
   armWWos[167]=armWWos[33]*armWWos[167];
   armWWos[169]=3307./12. + 83*armWWos[33];
   armWWos[169]=armWWos[15]*armWWos[169];
   armWWos[212]= - 29./6.*armWWos[17];
   armWWos[167]=1./6.*armWWos[169] + 1./2.*armWWos[167] + 3421./192.*
   armWWos[16] + armWWos[212] - 155./4.*armWWos[18] + armWWos[112] + 
   1373./64.*armWWos[19] - 3777./32.*armWWos[40] + 19./2.*armWWos[38]
    + 7./2.*armWWos[8] - 1./3.*armWWos[62] + 159./4.*armWWos[61] - 23*
   armWWos[63] + 41./32.*armWWos[36];
   armWWos[169]=armWWos[189] + armWWos[15];
   armWWos[165]=3*armWWos[169] + armWWos[165];
   armWWos[165]=armWWos[5]*armWWos[34]*armWWos[165];
   armWWos[213]=161./36. + 3*armWWos[54];
   armWWos[213]=armWWos[14] + 13./6.*armWWos[57] - 3./2.*armWWos[39] + 
   1./2.*armWWos[213] - armWWos[49];
   armWWos[193]=armWWos[33]*armWWos[193];
   armWWos[193]=1./6.*armWWos[193] + 1./2.*armWWos[213] + 13./9.*
   armWWos[10];
   armWWos[193]=1./2.*mmH*armWWos[193];
   armWWos[103]=armWWos[103] + armWWos[193] + 3*armWWos[165] + 1./4.*
   armWWos[167] + armWWos[148];
   armWWos[103]=armWWos[177]*armWWos[103];
   armWWos[148]=10627./24. - 1301*armWWos[9];
   armWWos[148]=1./16.*armWWos[148] + armWWos[174];
   armWWos[165]=armWWos[119] - 121./2.*armWWos[98];
   armWWos[165]=1./2.*armWWos[165] + armWWos[152];
   armWWos[165]=armWWos[15]*armWWos[165];
   armWWos[167]=69./4.*armWWos[97];
   armWWos[213]= - 407./16.*armWWos[98];
   armWWos[214]=armWWos[213] - 7*armWWos[32] + armWWos[167];
   armWWos[214]=armWWos[34]*armWWos[214];
   armWWos[215]= - 1./2.*armWWos[97];
   armWWos[216]= - 1./3.*armWWos[32] + armWWos[215];
   armWWos[217]=armWWos[216] - 77./16.*armWWos[98];
   armWWos[217]=armWWos[16]*armWWos[217];
   armWWos[218]= - 11209./1152.*armWWos[33];
   armWWos[148]=1./2.*armWWos[214] + armWWos[165] + armWWos[218] + 
   armWWos[178] + 1./2.*armWWos[217] + armWWos[10] + 1./3.*armWWos[148]
    + armWWos[171];
   armWWos[148]=armWWos[34]*armWWos[148];
   armWWos[165]= - 29*armWWos[63];
   armWWos[214]=97./6.*armWWos[8] - armWWos[62] + 833./12.*armWWos[61]
    + armWWos[165] - 3205./96.*armWWos[36];
   armWWos[217]=armWWos[156] - 3341./96.*armWWos[16];
   armWWos[217]=armWWos[33]*armWWos[217];
   armWWos[219]=109*armWWos[33];
   armWWos[220]= - 625./4. + armWWos[219];
   armWWos[220]=armWWos[15]*armWWos[220];
   armWWos[221]=59./2.*armWWos[38];
   armWWos[214]=1./6.*armWWos[220] + 1./2.*armWWos[217] + 12113./576.*
   armWWos[16] + armWWos[212] - 403./12.*armWWos[18] + armWWos[112] - 
   1087./192.*armWWos[19] - 17695./288.*armWWos[40] + 1./3.*
   armWWos[214] + armWWos[221];
   armWWos[217]= - 15*armWWos[34] + armWWos[16] + 9*armWWos[15];
   armWWos[217]=armWWos[5]*armWWos[34]*armWWos[217];
   armWWos[103]=armWWos[103] + armWWos[199] + armWWos[193] + 
   armWWos[217] + 1./4.*armWWos[214] + armWWos[148];
   armWWos[103]=armWWos[177]*armWWos[103];
   armWWos[148]=128387./128. - 580*armWWos[53];
   armWWos[148]=1./3.*armWWos[148] - 13./2.*armWWos[42];
   armWWos[148]=8851./144.*armWWos[52] + armWWos[198] + armWWos[197] + 
   1./3.*armWWos[148] + armWWos[190];
   armWWos[106]=armWWos[159] + armWWos[137] + armWWos[106] + 
   armWWos[210] + armWWos[133] + armWWos[168] + armWWos[135] + 
   armWWos[209] + armWWos[207] + armWWos[208] + 5935./432.*armWWos[9]
    + armWWos[147] - 6565./432.*armWWos[13] + armWWos[128] + 
   armWWos[206] + armWWos[205] + armWWos[204] + armWWos[146] + 
   armWWos[203] + armWWos[202] + armWWos[201] + armWWos[51] + 1./3.*
   armWWos[148] + armWWos[200];
   armWWos[106]=mmt*armWWos[106];
   armWWos[128]= - 5./3.*armWWos[32];
   armWWos[133]=armWWos[213] + armWWos[128] + armWWos[167];
   armWWos[133]=armWWos[34]*armWWos[133];
   armWWos[135]= - 197495./8. + 2371*armWWos[9];
   armWWos[135]=1./144.*armWWos[135] + armWWos[174];
   armWWos[119]=armWWos[119] - 107./6.*armWWos[98];
   armWWos[119]=1./2.*armWWos[119] + armWWos[152];
   armWWos[119]=armWWos[15]*armWWos[119];
   armWWos[137]=armWWos[216] + 793./48.*armWWos[98];
   armWWos[137]=armWWos[16]*armWWos[137];
   armWWos[119]=1./2.*armWWos[133] + armWWos[119] + armWWos[218] + 
   armWWos[178] + 1./2.*armWWos[137] + armWWos[10] + 1./3.*armWWos[135]
    + armWWos[171];
   armWWos[119]=armWWos[34]*armWWos[119];
   armWWos[133]=43./2.*armWWos[8] - armWWos[62] + 1345./12.*armWWos[61]
    + armWWos[165] + 1915./96.*armWWos[36];
   armWWos[135]=armWWos[156] - 2599./288.*armWWos[16];
   armWWos[135]=armWWos[33]*armWWos[135];
   armWWos[137]= - 1529./36. + armWWos[219];
   armWWos[137]=armWWos[15]*armWWos[137];
   armWWos[112]=1./6.*armWWos[137] + 1./2.*armWWos[135] + 21235./1728.*
   armWWos[16] + armWWos[212] - 2617./36.*armWWos[18] + armWWos[112] - 
   789./64.*armWWos[19] - 36127./288.*armWWos[40] + 1./3.*armWWos[133]
    + armWWos[221];
   armWWos[133]= - 275./24. + 7*armWWos[50];
   armWWos[135]= - 2*armWWos[56];
   armWWos[137]=1./4.*armWWos[30];
   armWWos[147]= - 31./8. + 25*armWWos[9];
   armWWos[147]=1./4.*armWWos[147] + armWWos[33];
   armWWos[147]=armWWos[33]*armWWos[147];
   armWWos[133]=1./2.*armWWos[147] - 671./72.*armWWos[9] + 205./36.*
   armWWos[13] - 143./64.*armWWos[39] + armWWos[137] + 5./64.*
   armWWos[60] - 605./96.*armWWos[58] - 473./72.*armWWos[52] + 271./64.
   *armWWos[55] + 1./8.*armWWos[133] + armWWos[135];
   armWWos[147]=1./2.*armWWos[36];
   armWWos[145]=5*armWWos[16] + armWWos[145] + armWWos[147] - 
   armWWos[40];
   armWWos[140]=5./6. + armWWos[140];
   armWWos[140]=armWWos[34]*armWWos[140];
   armWWos[148]= - armWWos[5]*armWWos[34]*armWWos[16];
   armWWos[107]=18*armWWos[148] + armWWos[140] + 7./3.*armWWos[145] + 3.
   /4.*armWWos[107];
   armWWos[107]=armWWos[5]*armWWos[107];
   armWWos[140]=19*armWWos[45] + armWWos[48];
   armWWos[145]= - armWWos[5]*armWWos[16];
   armWWos[145]=9*armWWos[145] - 5./3. - armWWos[33];
   armWWos[145]=armWWos[5]*armWWos[145];
   armWWos[140]=1./6.*armWWos[140] + armWWos[145];
   armWWos[140]=mmt*armWWos[140];
   armWWos[107]=armWWos[140] + 1./3.*armWWos[133] + armWWos[107];
   armWWos[107]=mmt*armWWos[107];
   armWWos[133]=2191./8.*armWWos[36] + 31*armWWos[61];
   armWWos[133]=1./2.*armWWos[133] + 115*armWWos[8];
   armWWos[140]=4841./108. + armWWos[33];
   armWWos[140]=armWWos[15]*armWWos[140];
   armWWos[133]=1./3.*armWWos[140] + 329./288.*armWWos[122] - 23803./
   2592.*armWWos[16] + 175./54.*armWWos[18] - 8515./864.*armWWos[19] + 
   1./27.*armWWos[133] - 179./16.*armWWos[40];
   armWWos[140]=11*armWWos[32];
   armWWos[145]=armWWos[140] + armWWos[215];
   armWWos[148]=armWWos[145] - 515./12.*armWWos[98];
   armWWos[148]=armWWos[16]*armWWos[148];
   armWWos[152]= - 851./8. + 301./9.*armWWos[9];
   armWWos[148]= - 233./48.*armWWos[33] + 1./2.*armWWos[152] + 
   armWWos[148];
   armWWos[152]= - 11*armWWos[32];
   armWWos[156]=armWWos[152] + 485./4.*armWWos[98];
   armWWos[156]=armWWos[34]*armWWos[156];
   armWWos[159]=armWWos[15]*armWWos[98];
   armWWos[148]=1./2.*armWWos[156] + 1./2.*armWWos[148] + 235./3.*
   armWWos[159];
   armWWos[148]=armWWos[34]*armWWos[148];
   armWWos[156]=43*armWWos[16] - 17*armWWos[34];
   armWWos[156]=armWWos[5]*armWWos[34]*armWWos[156];
   armWWos[133]=1./3.*armWWos[156] + 1./2.*armWWos[133] + 1./3.*
   armWWos[148];
   armWWos[148]= - 41./18.*armWWos[52] + 331./144.*armWWos[55] - 45./16.
    + armWWos[50];
   armWWos[109]= - 5./36.*armWWos[9] + 11./36.*armWWos[13] - 331./288.*
   armWWos[39] - 137./96.*armWWos[60] - 313./144.*armWWos[58] + 1./2.*
   armWWos[148] + armWWos[109];
   armWWos[148]= - armWWos[43] - 1./6.*armWWos[48];
   armWWos[156]= - armWWos[5]*armWWos[33];
   armWWos[148]=1./3.*armWWos[148] + armWWos[156];
   armWWos[148]=mmt*armWWos[148];
   armWWos[156]= - 25./64. + armWWos[9];
   armWWos[156]=armWWos[33]*armWWos[156];
   armWWos[165]= - 1 - armWWos[9];
   armWWos[165]=armWWos[5]*armWWos[34]*armWWos[165];
   armWWos[109]=armWWos[148] + armWWos[165] + 1./2.*armWWos[109] + 1./3.
   *armWWos[156];
   armWWos[109]=mmt*armWWos[109];
   armWWos[148]=armWWos[16]*armWWos[98];
   armWWos[156]=25./2.*armWWos[148] + 281./8. + 31*armWWos[9];
   armWWos[156]=1./3.*armWWos[156] - 25./8.*armWWos[33];
   armWWos[165]=armWWos[34]*armWWos[98];
   armWWos[156]=25./8.*armWWos[165] + 1./4.*armWWos[156] + 25./3.*
   armWWos[159];
   armWWos[156]=armWWos[34]*armWWos[156];
   armWWos[165]= - 1./9.*armWWos[63] + 13./128.*armWWos[36];
   armWWos[109]=1./2.*armWWos[109] + 1./12.*armWWos[156] + 161./864.*
   armWWos[15] + 25./768.*armWWos[122] + 2057./6912.*armWWos[16] + 19./
   144.*armWWos[18] - 679./2304.*armWWos[19] - 797./1152.*armWWos[40]
    + 7./72.*armWWos[8] + 5*armWWos[165] - 1./16.*armWWos[61];
   armWWos[109]=armWWos[121]*armWWos[109];
   armWWos[107]=armWWos[109] + 1./2.*armWWos[133] + armWWos[107];
   armWWos[107]=armWWos[121]*armWWos[107];
   armWWos[106]=armWWos[107] + armWWos[106] + armWWos[193] + 
   armWWos[217] + 1./4.*armWWos[112] + armWWos[119];
   armWWos[106]=armWWos[121]*armWWos[106];
   armWWos[107]=armWWos[33]*armWWos[195];
   armWWos[107]= - 7./3. + armWWos[107];
   armWWos[109]=4 + armWWos[173];
   armWWos[109]=armWWos[5]*armWWos[34]*armWWos[109];
   armWWos[107]=1./2.*armWWos[107] + 2*armWWos[109];
   armWWos[107]=mmt*armWWos[107];
   armWWos[109]= - 5./3. + armWWos[115];
   armWWos[109]=armWWos[34]*armWWos[109];
   armWWos[107]=armWWos[109] + armWWos[107];
   armWWos[109]=armWWos[177]*armWWos[107];
   armWWos[109]=armWWos[107] + armWWos[109];
   armWWos[109]=armWWos[177]*armWWos[109];
   armWWos[107]=armWWos[121]*armWWos[107];
   armWWos[107]=armWWos[109] + armWWos[107];
   armWWos[107]=armWWos[31]*armWWos[107];
   armWWos[109]= - 1 + armWWos[9];
   armWWos[109]=armWWos[34]*armWWos[109];
   armWWos[112]=armWWos[9] - armWWos[13] + armWWos[52] + 1 - 
   armWWos[53];
   armWWos[112]=mmt*armWWos[112];
   armWWos[109]=armWWos[109] + armWWos[112];
   armWWos[103]=armWWos[107] + armWWos[106] + 256./9.*armWWos[109] + 
   armWWos[103];
   armWWos[103]=armWWos[31]*armWWos[103];
   armWWos[106]= - 5425./18. + 167*armWWos[9];
   armWWos[107]=6*armWWos[10];
   armWWos[109]= - armWWos[16]*armWWos[99];
   armWWos[112]=15./4.*armWWos[109];
   armWWos[115]=25*armWWos[17] + 29./2.*armWWos[16];
   armWWos[115]=1./6.*armWWos[100]*armWWos[115];
   armWWos[119]=33./4.*armWWos[15]*armWWos[100];
   armWWos[106]=armWWos[119] + armWWos[115] + armWWos[112] + 1./4.*
   armWWos[106] + armWWos[107];
   armWWos[106]=armWWos[15]*armWWos[106];
   armWWos[133]=10123./576. + 10*armWWos[90];
   armWWos[156]= - 7./12. - armWWos[9];
   armWWos[156]=armWWos[9]*armWWos[156];
   armWWos[165]=931./9. - 71./2.*armWWos[9];
   armWWos[165]=1./6.*armWWos[165] + armWWos[10];
   armWWos[165]=armWWos[10]*armWWos[165];
   armWWos[167]= - 4./3.*armWWos[78] + armWWos[75];
   armWWos[167]= - 23./24.*armWWos[71] + 17./24.*armWWos[74] + 17./48.*
   armWWos[96] + 1./6.*armWWos[70] + 2*armWWos[167] + 3./4.*armWWos[73]
   ;
   armWWos[167]=mmH*armWWos[167];
   armWWos[133]=armWWos[167] + 1./4.*armWWos[165] + 1./4.*armWWos[156]
    - 215./48.*armWWos[14] - 1./8.*armWWos[13] + 8./3.*armWWos[30] - 
   197./24.*armWWos[89] + 11./12.*armWWos[85] + 89./16.*armWWos[91] + 
   armWWos[111] + 25./48.*armWWos[87] - 9./16.*armWWos[81] - 9./16.*
   armWWos[82] + 11./6.*armWWos[86] + 2./3.*armWWos[88] - 1./2.*
   armWWos[83] + 1./3.*armWWos[133] + 9./8.*armWWos[84];
   armWWos[133]=mmH*armWWos[133];
   armWWos[156]=9*armWWos[69] + 11./3.*armWWos[93];
   armWWos[165]=armWWos[156] - 605./6.*armWWos[94];
   armWWos[167]=armWWos[165] + 149./3.*armWWos[66];
   armWWos[168]=47./3.*armWWos[67];
   armWWos[167]=179./6.*armWWos[95] + 1./2.*armWWos[167] + armWWos[168]
   ;
   armWWos[173]=29./3.*armWWos[68];
   armWWos[190]= - 85./8.*armWWos[20];
   armWWos[193]=3./2.*armWWos[10];
   armWWos[195]=armWWos[193] + 95./9. - 49./8.*armWWos[9];
   armWWos[195]=armWWos[17]*armWWos[195];
   armWWos[167]=armWWos[195] + 6175./72.*armWWos[18] + armWWos[190] - 
   599./8.*armWWos[19] + 1./2.*armWWos[167] + armWWos[173];
   armWWos[197]= - 557./27. + 125*armWWos[9];
   armWWos[198]=armWWos[17]*armWWos[99];
   armWWos[197]=73./24.*armWWos[109] + 7./6.*armWWos[198] + 11./16.*
   armWWos[197] - 23./3.*armWWos[10];
   armWWos[197]=armWWos[16]*armWWos[197];
   armWWos[199]=3*armWWos[17];
   armWWos[200]=armWWos[199] - 7./8.*armWWos[16];
   armWWos[200]=armWWos[16]*armWWos[200];
   armWWos[201]=59./4.*armWWos[15] + 14./3.*armWWos[17] + 69./8.*
   armWWos[16];
   armWWos[201]=armWWos[15]*armWWos[201];
   armWWos[200]=armWWos[201] + 3./2.*armWWos[142] + armWWos[200];
   armWWos[200]=armWWos[5]*armWWos[200];
   armWWos[151]=armWWos[151] + 5./2.*armWWos[16];
   armWWos[151]=armWWos[16]*armWWos[151];
   armWWos[151]= - 1./3.*armWWos[142] + armWWos[151];
   armWWos[151]=1./4.*armWWos[100]*armWWos[151];
   armWWos[106]=armWWos[133] + armWWos[200] + armWWos[106] + 
   armWWos[151] + 1./2.*armWWos[167] + armWWos[197];
   armWWos[106]=armWWos[65]*armWWos[106];
   armWWos[167]=11./24.*armWWos[96];
   armWWos[197]=11./12.*armWWos[74];
   armWWos[201]=11./12.*armWWos[71] + armWWos[197] + armWWos[167] + 3./
   2.*armWWos[73] + 1./3.*armWWos[70];
   armWWos[201]=mmH*armWWos[201];
   armWWos[202]=9*armWWos[84];
   armWWos[203]=97./216. + armWWos[202];
   armWWos[204]= - 1./3.*armWWos[9];
   armWWos[205]= - 7./8. + armWWos[204];
   armWWos[205]=armWWos[9]*armWWos[205];
   armWWos[206]= - 1247./9. - 161./2.*armWWos[9];
   armWWos[206]=1./6.*armWWos[206] + armWWos[10];
   armWWos[206]=armWWos[10]*armWWos[206];
   armWWos[111]=1./2.*armWWos[201] + 1./4.*armWWos[206] + 1./2.*
   armWWos[205] - 137./48.*armWWos[14] + 15./8.*armWWos[13] - 95./8.*
   armWWos[89] + 5./6.*armWWos[85] + armWWos[111] - 29./48.*armWWos[87]
    - 29./48.*armWWos[81] - 29./48.*armWWos[82] + armWWos[130] - 1./3.*
   armWWos[88] + 1./8.*armWWos[203] - 1./3.*armWWos[83];
   armWWos[111]=mmH*armWWos[111];
   armWWos[130]= - 327*armWWos[9];
   armWWos[201]=45421./18. + armWWos[130];
   armWWos[201]=1./2.*armWWos[201] + 29*armWWos[10];
   armWWos[201]=19./4.*armWWos[109] + 1./2.*armWWos[201] + 3*
   armWWos[198];
   armWWos[201]=armWWos[16]*armWWos[201];
   armWWos[189]=armWWos[199] + armWWos[189];
   armWWos[189]=armWWos[16]*armWWos[189];
   armWWos[142]=23./6.*armWWos[142] + armWWos[189];
   armWWos[142]=armWWos[100]*armWWos[142];
   armWWos[156]= - 489./2.*armWWos[66] + armWWos[156] - 197./2.*
   armWWos[94];
   armWWos[189]=armWWos[182] + 517./9. - 147./4.*armWWos[9];
   armWWos[189]=armWWos[17]*armWWos[189];
   armWWos[142]=armWWos[142] + armWWos[201] + 1./2.*armWWos[189] + 2349.
   /8.*armWWos[18] - 151./8.*armWWos[20] + 1571./4.*armWWos[19] + 11./2.
   *armWWos[68] - 159*armWWos[95] + 1./4.*armWWos[156] + 17./3.*
   armWWos[67];
   armWWos[156]= - 22003./27. - 1121*armWWos[9];
   armWWos[156]=7*armWWos[109] + 1./4.*armWWos[156] + 31*armWWos[10];
   armWWos[120]=8*armWWos[17] + armWWos[120];
   armWWos[120]=armWWos[100]*armWWos[120];
   armWWos[120]=7./2.*armWWos[144] + 1./4.*armWWos[156] + armWWos[120];
   armWWos[120]=armWWos[15]*armWWos[120];
   armWWos[156]= - 2*armWWos[16];
   armWWos[189]=3./2.*armWWos[17] + armWWos[156];
   armWWos[189]=armWWos[16]*armWWos[189];
   armWWos[123]= - 35./4.*armWWos[15] + armWWos[123] - 109./8.*
   armWWos[16];
   armWWos[123]=armWWos[15]*armWWos[123];
   armWWos[123]=armWWos[123] + armWWos[160] + armWWos[189];
   armWWos[123]=armWWos[5]*armWWos[123];
   armWWos[111]=armWWos[111] + armWWos[123] + 1./2.*armWWos[142] + 
   armWWos[120];
   armWWos[111]=armWWos[65]*armWWos[111];
   armWWos[120]=1./2.*armWWos[8];
   armWWos[123]=armWWos[7] + armWWos[120];
   armWWos[142]=1./2.*armWWos[20];
   armWWos[123]=35./2.*armWWos[18] + armWWos[142] + 11*armWWos[123] + 
   13*armWWos[19];
   armWWos[160]=armWWos[1]*armWWos[123];
   armWWos[123]=armWWos[35]*armWWos[123];
   armWWos[189]= - 23./18. + armWWos[6];
   armWWos[189]=armWWos[1]*armWWos[189];
   armWWos[199]=3*armWWos[6];
   armWWos[201]= - 23./6. + armWWos[199];
   armWWos[201]=armWWos[35]*armWWos[201];
   armWWos[189]=armWWos[189] + armWWos[201];
   armWWos[189]=1./2.*armWWos[17]*armWWos[189];
   armWWos[123]=armWWos[189] + 1./3.*armWWos[160] + armWWos[123];
   armWWos[160]= - armWWos[149]*armWWos[99];
   armWWos[201]=armWWos[100]*armWWos[149];
   armWWos[160]=armWWos[160] + armWWos[201];
   armWWos[201]=armWWos[100]*armWWos[16];
   armWWos[201]=armWWos[109] + armWWos[201];
   armWWos[201]=armWWos[15]*armWWos[201];
   armWWos[160]=1./2.*armWWos[160] + armWWos[201];
   armWWos[160]=armWWos[177]*armWWos[65]*armWWos[160];
   armWWos[201]=13./3. + armWWos[166];
   armWWos[203]=armWWos[1]*armWWos[201];
   armWWos[201]=armWWos[35]*armWWos[201];
   armWWos[201]=1./3.*armWWos[203] + armWWos[201];
   armWWos[201]=armWWos[10]*armWWos[201];
   armWWos[203]= - 1./12.*armWWos[6] + 1./4.*armWWos[14] - 1./24.*
   armWWos[12] - 1./4.*armWWos[25] + 4./9. + 1./4.*armWWos[28];
   armWWos[205]=armWWos[1]*armWWos[203];
   armWWos[203]=armWWos[35]*armWWos[203];
   armWWos[201]=1./6.*armWWos[201] + 1./3.*armWWos[205] + armWWos[203];
   armWWos[201]=mmH*armWWos[201];
   armWWos[203]= - 163./24. + 4*armWWos[6];
   armWWos[205]=armWWos[1]*armWWos[203];
   armWWos[203]=armWWos[35]*armWWos[203];
   armWWos[203]=1./3.*armWWos[205] + armWWos[203];
   armWWos[203]=armWWos[16]*armWWos[203];
   armWWos[205]= - 989./6. + 79*armWWos[6];
   armWWos[206]=armWWos[1]*armWWos[205];
   armWWos[205]=armWWos[35]*armWWos[205];
   armWWos[205]=1./3.*armWWos[206] + armWWos[205];
   armWWos[205]=armWWos[15]*armWWos[205];
   armWWos[111]=3./4.*armWWos[160] + armWWos[111] + armWWos[201] + 1./
   12.*armWWos[205] + 1./2.*armWWos[123] + armWWos[203];
   armWWos[111]=armWWos[177]*armWWos[111];
   armWWos[123]=7*armWWos[7];
   armWWos[160]=1./6.*armWWos[20];
   armWWos[203]=5./18.*armWWos[18] + armWWos[160] - 19./9.*armWWos[19]
    + armWWos[123] + 17./18.*armWWos[8];
   armWWos[203]=armWWos[1]*armWWos[203];
   armWWos[205]=21*armWWos[7];
   armWWos[206]=5./6.*armWWos[18] + armWWos[142] - 19./3.*armWWos[19]
    + armWWos[205] + 17./6.*armWWos[8];
   armWWos[206]=armWWos[35]*armWWos[206];
   armWWos[207]= - 41./6. - 61*armWWos[6];
   armWWos[208]=armWWos[1]*armWWos[207];
   armWWos[207]=armWWos[35]*armWWos[207];
   armWWos[207]=1./3.*armWWos[208] + armWWos[207];
   armWWos[207]=armWWos[16]*armWWos[207];
   armWWos[203]=1./6.*armWWos[207] + armWWos[189] + armWWos[203] + 
   armWWos[206];
   armWWos[206]= - 235./8. + 23./3.*armWWos[6];
   armWWos[207]=armWWos[1]*armWWos[206];
   armWWos[206]=armWWos[35]*armWWos[206];
   armWWos[206]=1./3.*armWWos[207] + armWWos[206];
   armWWos[206]=armWWos[15]*armWWos[206];
   armWWos[106]=armWWos[111] + armWWos[106] + armWWos[201] + 1./2.*
   armWWos[203] + armWWos[206];
   armWWos[106]=armWWos[177]*armWWos[106];
   armWWos[111]= - 2449./6. - 307*armWWos[9];
   armWWos[107]=armWWos[119] + armWWos[115] + armWWos[112] + 1./12.*
   armWWos[111] + armWWos[107];
   armWWos[107]=armWWos[15]*armWWos[107];
   armWWos[111]=armWWos[165] - 467./3.*armWWos[66];
   armWWos[111]= - 365./6.*armWWos[95] + 1./2.*armWWos[111] + 
   armWWos[168];
   armWWos[111]=armWWos[195] + 127./72.*armWWos[18] + armWWos[190] + 
   2315./24.*armWWos[19] + 1./2.*armWWos[111] + armWWos[173];
   armWWos[112]=35633./9. + 989*armWWos[9];
   armWWos[112]=73./8.*armWWos[109] + 7./2.*armWWos[198] + 1./16.*
   armWWos[112] - 11*armWWos[10];
   armWWos[112]=armWWos[16]*armWWos[112];
   armWWos[107]=armWWos[133] + armWWos[200] + armWWos[107] + 
   armWWos[151] + 1./2.*armWWos[111] + 1./3.*armWWos[112];
   armWWos[107]=armWWos[65]*armWWos[107];
   armWWos[111]=1453./3. + 1429*armWWos[9];
   armWWos[111]=1./36.*armWWos[111] + armWWos[10];
   armWWos[112]=armWWos[17] + 11./2.*armWWos[16];
   armWWos[112]=armWWos[100]*armWWos[112];
   armWWos[111]=1./4.*armWWos[144] + 1./6.*armWWos[112] + 1./4.*
   armWWos[111] + 2*armWWos[109];
   armWWos[111]=armWWos[15]*armWWos[111];
   armWWos[112]=13./16. - armWWos[83];
   armWWos[112]=17./8.*armWWos[81] + 17./8.*armWWos[82] + armWWos[161]
    + 1./2.*armWWos[112] - armWWos[88];
   armWWos[115]= - 1./2. - armWWos[9];
   armWWos[115]=armWWos[9]*armWWos[115];
   armWWos[119]=1./2.*armWWos[9];
   armWWos[133]= - 1./3. + armWWos[119];
   armWWos[133]=armWWos[10]*armWWos[133];
   armWWos[144]=1./2.*armWWos[96] + armWWos[74];
   armWWos[151]=armWWos[144] + 7./3.*armWWos[71];
   armWWos[151]=mmH*armWWos[151];
   armWWos[112]=1./4.*armWWos[151] + 29./12.*armWWos[133] + 1./12.*
   armWWos[115] - 23./24.*armWWos[14] - 1./3.*armWWos[13] - 16./9.*
   armWWos[89] + 1./12.*armWWos[85] + 155./48.*armWWos[91] + 1./3.*
   armWWos[112] + 1./8.*armWWos[87];
   armWWos[112]=mmH*armWWos[112];
   armWWos[115]= - 4271./9. - 313./2.*armWWos[9];
   armWWos[115]=1./6.*armWWos[115] - armWWos[10];
   armWWos[133]= - armWWos[17]*armWWos[99];
   armWWos[151]=2./3.*armWWos[109];
   armWWos[115]=armWWos[151] + 1./4.*armWWos[115] + 1./3.*armWWos[133];
   armWWos[115]=armWWos[16]*armWWos[115];
   armWWos[133]= - 119*armWWos[94] + 1031./6.*armWWos[66];
   armWWos[133]= - 2051./6.*armWWos[18] + 73./6.*armWWos[20] - 815./4.*
   armWWos[19] + 25*armWWos[68] + 151./2.*armWWos[95] + 1./2.*
   armWWos[133] + armWWos[67];
   armWWos[161]=13 + 43./6.*armWWos[9];
   armWWos[161]=armWWos[17]*armWWos[161];
   armWWos[133]=1./3.*armWWos[133] + armWWos[161];
   armWWos[161]=armWWos[17] + 3./4.*armWWos[16];
   armWWos[161]=armWWos[16]*armWWos[161];
   armWWos[165]=17./2.*armWWos[16] - armWWos[15];
   armWWos[165]=armWWos[15]*armWWos[165];
   armWWos[161]=3*armWWos[161] + armWWos[165];
   armWWos[161]=armWWos[5]*armWWos[161];
   armWWos[165]= - 1./3.*armWWos[17] + 3./2.*armWWos[16];
   armWWos[165]=armWWos[100]*armWWos[16]*armWWos[165];
   armWWos[111]=armWWos[112] + 1./2.*armWWos[161] + armWWos[111] + 1./4.
   *armWWos[165] + 1./4.*armWWos[133] + armWWos[115];
   armWWos[111]=armWWos[65]*armWWos[111];
   armWWos[112]= - 3*armWWos[87];
   armWWos[115]= - 1./3.*armWWos[91] + armWWos[112] + armWWos[81] - 61./
   36. + armWWos[82];
   armWWos[133]= - 7 + 13./6.*armWWos[9];
   armWWos[133]=armWWos[10]*armWWos[133];
   armWWos[161]= - 1./2.*armWWos[96] - armWWos[74];
   armWWos[165]=armWWos[161] - armWWos[71];
   armWWos[165]=mmH*armWWos[165];
   armWWos[115]=1./3.*armWWos[165] + 1./3.*armWWos[133] + 5./18.*
   armWWos[9] + armWWos[134] + armWWos[13] - armWWos[89] + 1./2.*
   armWWos[115] + 7./3.*armWWos[85];
   armWWos[115]=mmH*armWWos[115];
   armWWos[133]= - 17 + armWWos[179];
   armWWos[133]=23./36.*armWWos[133] + armWWos[10];
   armWWos[133]=1./2.*armWWos[133] + 1./3.*armWWos[198];
   armWWos[133]=1./2.*armWWos[133] + armWWos[151];
   armWWos[133]=armWWos[16]*armWWos[133];
   armWWos[134]=11./2.*armWWos[9];
   armWWos[151]=43./9. + armWWos[134];
   armWWos[109]=1./4.*armWWos[151] + armWWos[109];
   armWWos[109]=armWWos[15]*armWWos[109];
   armWWos[151]= - armWWos[15]*armWWos[16];
   armWWos[151]=armWWos[149] + armWWos[151];
   armWWos[151]=armWWos[5]*armWWos[151];
   armWWos[165]= - 1./9. + 1./4.*armWWos[9];
   armWWos[165]=armWWos[17]*armWWos[165];
   armWWos[165]=armWWos[165] - 355./36.*armWWos[18] - 43./36.*
   armWWos[20] - 121./12.*armWWos[19] + 11./6.*armWWos[68] + 7./3.*
   armWWos[95] - 41./36.*armWWos[94] + 3*armWWos[66];
   armWWos[109]=1./8.*armWWos[115] + 1./8.*armWWos[151] + armWWos[109]
    + 1./4.*armWWos[165] + armWWos[133];
   armWWos[109]=armWWos[65]*armWWos[109];
   armWWos[115]= - 5./2.*armWWos[19] + armWWos[132] + armWWos[95];
   armWWos[115]=1./2.*armWWos[115] - armWWos[18];
   armWWos[132]= - 1 - 7./18.*armWWos[9];
   armWWos[132]=armWWos[16]*armWWos[132];
   armWWos[133]=1 + 7./9.*armWWos[9];
   armWWos[133]=armWWos[15]*armWWos[133];
   armWWos[151]= - 1 - armWWos[91];
   armWWos[151]=mmH*armWWos[151];
   armWWos[115]=1./12.*armWWos[151] + 1./4.*armWWos[133] + 1./3.*
   armWWos[115] + 1./2.*armWWos[132];
   armWWos[115]=armWWos[121]*armWWos[65]*armWWos[115];
   armWWos[132]=armWWos[18] + 3*armWWos[8] - armWWos[19];
   armWWos[132]=armWWos[1]*armWWos[132];
   armWWos[133]= - 11*armWWos[19];
   armWWos[151]=11*armWWos[18] + 17*armWWos[8] + armWWos[133];
   armWWos[151]=armWWos[35]*armWWos[151];
   armWWos[132]=armWWos[132] + 1./9.*armWWos[151];
   armWWos[151]= - armWWos[1] - 5./9.*armWWos[35];
   armWWos[151]=armWWos[16]*armWWos[151];
   armWWos[165]= - armWWos[1] - 11./9.*armWWos[35];
   armWWos[168]=armWWos[15]*armWWos[165];
   armWWos[109]=1./4.*armWWos[115] + armWWos[109] + 5./24.*armWWos[168]
    + 1./4.*armWWos[132] + 2./3.*armWWos[151];
   armWWos[109]=armWWos[121]*armWWos[109];
   armWWos[115]=5./3.*armWWos[18] + 4./3.*armWWos[8] - armWWos[19];
   armWWos[115]=armWWos[1]*armWWos[115];
   armWWos[132]=55*armWWos[18] + 40*armWWos[8] - 37*armWWos[19];
   armWWos[132]=armWWos[35]*armWWos[132];
   armWWos[115]=armWWos[115] + 1./27.*armWWos[132];
   armWWos[132]= - 71./3. + 11./4.*armWWos[6];
   armWWos[132]=armWWos[1]*armWWos[132];
   armWWos[151]= - 179./9. + 1./4.*armWWos[6];
   armWWos[151]=armWWos[35]*armWWos[151];
   armWWos[132]=armWWos[132] + armWWos[151];
   armWWos[132]=armWWos[16]*armWWos[132];
   armWWos[151]= - 235./3. + armWWos[6];
   armWWos[151]=armWWos[1]*armWWos[151];
   armWWos[168]= - 835./27. + armWWos[6];
   armWWos[168]=armWWos[35]*armWWos[168];
   armWWos[151]=1./3.*armWWos[151] + armWWos[168];
   armWWos[151]=armWWos[15]*armWWos[151];
   armWWos[109]=armWWos[109] + armWWos[111] + 1./12.*armWWos[151] + 2*
   armWWos[115] + 1./9.*armWWos[132];
   armWWos[109]=armWWos[121]*armWWos[109];
   armWWos[111]= - 91./18.*armWWos[18] + armWWos[160] - 13./9.*
   armWWos[19] + armWWos[123] + 65./18.*armWWos[8];
   armWWos[111]=armWWos[1]*armWWos[111];
   armWWos[115]= - 145./18.*armWWos[18] + armWWos[142] - 31./9.*
   armWWos[19] + armWWos[205] + 131./18.*armWWos[8];
   armWWos[115]=armWWos[35]*armWWos[115];
   armWWos[123]=79./6. - armWWos[6];
   armWWos[123]=armWWos[1]*armWWos[123];
   armWWos[132]= - 115./18. - 17*armWWos[6];
   armWWos[132]=armWWos[35]*armWWos[132];
   armWWos[123]=1./3.*armWWos[123] + armWWos[132];
   armWWos[123]=armWWos[16]*armWWos[123];
   armWWos[111]=1./6.*armWWos[123] + armWWos[189] + armWWos[111] + 
   armWWos[115];
   armWWos[115]=23*armWWos[6];
   armWWos[123]=191./8. + armWWos[115];
   armWWos[123]=armWWos[1]*armWWos[123];
   armWWos[115]= - 1865./72. + armWWos[115];
   armWWos[115]=armWWos[35]*armWWos[115];
   armWWos[115]=1./3.*armWWos[123] + armWWos[115];
   armWWos[115]=armWWos[15]*armWWos[115];
   armWWos[107]=armWWos[109] + armWWos[107] + armWWos[201] + 1./2.*
   armWWos[111] + 1./3.*armWWos[115];
   armWWos[107]=armWWos[121]*armWWos[107];
   armWWos[109]= - armWWos[16]*armWWos[9];
   armWWos[111]=1 - armWWos[9];
   armWWos[115]=armWWos[15]*armWWos[111];
   armWWos[109]=4*armWWos[115] + 2*armWWos[109] - 4*armWWos[18] - 2*
   armWWos[66] + 17./3.*armWWos[19];
   armWWos[109]=armWWos[65]*armWWos[109];
   armWWos[101]=armWWos[101] + armWWos[103] + armWWos[107] + 8*
   armWWos[109] + armWWos[106];
   armWWos[101]=armWWos[2]*armWWos[101];
   armWWos[103]= - 7./4. - armWWos[58];
   armWWos[103]=1./2.*armWWos[103] - armWWos[60];
   armWWos[106]=pow(armWWos[97],2);
   armWWos[103]=mmZ*armWWos[106]*armWWos[103];
   armWWos[107]=3*armWWos[155] + armWWos[48];
   armWWos[109]=1./2.*armWWos[55] - 3./8. + armWWos[56];
   armWWos[109]= - 3./4.*armWWos[39] + 1./2.*armWWos[109] + armWWos[30]
   ;
   armWWos[109]=armWWos[32]*armWWos[109];
   armWWos[115]= - armWWos[63] + armWWos[120];
   armWWos[115]=armWWos[97]*armWWos[115];
   armWWos[120]=3./2.*armWWos[115] - 1./4.*armWWos[60] + 103./16. + 
   armWWos[58];
   armWWos[120]=armWWos[97]*armWWos[120];
   armWWos[123]= - 205 - 9*armWWos[55];
   armWWos[123]= - 11*armWWos[13] + 9./8.*armWWos[39] + 9./8.*
   armWWos[60] - 9./4.*armWWos[58] + 1./8.*armWWos[123] + 11*
   armWWos[52];
   armWWos[123]=armWWos[98]*armWWos[123];
   armWWos[132]=2./3.*armWWos[10] + armWWos[14] - 6*armWWos[57] - 53./
   12. - armWWos[49];
   armWWos[132]=armWWos[100]*armWWos[132];
   armWWos[151]=1 - armWWos[10];
   armWWos[151]=armWWos[100]*armWWos[151];
   armWWos[151]=81./512.*armWWos[98] + armWWos[151];
   armWWos[151]=armWWos[33]*armWWos[151];
   armWWos[155]=armWWos[16]*armWWos[106];
   armWWos[160]= - armWWos[9]*armWWos[98];
   armWWos[103]=3./4.*armWWos[103] + armWWos[151] + armWWos[132] + 3./4.
   *armWWos[155] + 99./64.*armWWos[160] + 9./64.*armWWos[123] + 1./2.*
   armWWos[120] + 2*armWWos[107] + armWWos[109];
   armWWos[103]=mmZ*armWWos[103];
   armWWos[107]=armWWos[36] + armWWos[8];
   armWWos[107]= - armWWos[19] - 3*armWWos[40] + 1./2.*armWWos[107] + 2
   *armWWos[38];
   armWWos[109]= - 6*armWWos[18];
   armWWos[102]= - armWWos[17] + armWWos[102];
   armWWos[102]=armWWos[33]*armWWos[102];
   armWWos[120]= - 3 + armWWos[33];
   armWWos[120]=armWWos[15]*armWWos[120];
   armWWos[123]=7./2. + 3*armWWos[33];
   armWWos[123]=mmZ*armWWos[123];
   armWWos[102]=armWWos[110] + armWWos[123] + 3*armWWos[120] + 1./2.*
   armWWos[102] - 9./2.*armWWos[16] + 1./2.*armWWos[17] + armWWos[109]
    + 3*armWWos[107] + armWWos[180];
   armWWos[102]=armWWos[5]*armWWos[102];
   armWWos[107]= - 7./2.*armWWos[14];
   armWWos[110]=armWWos[139] + armWWos[107] + armWWos[181] - 
   armWWos[59] + armWWos[146] + 7./4. + armWWos[188];
   armWWos[110]=armWWos[100]*armWWos[110];
   armWWos[120]=19./2. - armWWos[33];
   armWWos[123]= - armWWos[5]*armWWos[34];
   armWWos[132]=21*armWWos[123];
   armWWos[120]=1./2.*armWWos[120] + armWWos[132];
   armWWos[120]=armWWos[5]*armWWos[120];
   armWWos[139]= - armWWos[45] - 1./2.*armWWos[43];
   armWWos[146]= - 2*armWWos[48];
   armWWos[124]= - mmt*armWWos[124];
   armWWos[110]=18*armWWos[124] + 3*armWWos[120] + armWWos[196] + 
   armWWos[110] - 2*armWWos[44] + 7*armWWos[139] + armWWos[146];
   armWWos[110]=mmt*armWWos[110];
   armWWos[120]= - 18919./1536. + armWWos[174];
   armWWos[139]=armWWos[32] - 81./64.*armWWos[98];
   armWWos[139]=armWWos[16]*armWWos[139];
   armWWos[120]= - 5./4.*armWWos[33] + 1./2.*armWWos[178] + 1./4.*
   armWWos[139] + 1./3.*armWWos[120] + armWWos[171];
   armWWos[120]=armWWos[33]*armWWos[120];
   armWWos[139]=armWWos[9]*armWWos[98];
   armWWos[151]=33./4.*armWWos[139] - 561./32.*armWWos[98] + 
   armWWos[32] - 61*armWWos[97];
   armWWos[168]=2 - armWWos[10];
   armWWos[168]=armWWos[100]*armWWos[168];
   armWWos[171]=81./32.*armWWos[98] - armWWos[32] + 9./2.*armWWos[97];
   armWWos[171]=armWWos[33]*armWWos[171];
   armWWos[173]=mmZ*armWWos[106];
   armWWos[151]=9./8.*armWWos[173] + 1./4.*armWWos[171] + 3./4.*
   armWWos[151] + 2*armWWos[168];
   armWWos[151]=armWWos[34]*armWWos[151];
   armWWos[168]=78289./3456. - 33*armWWos[50];
   armWWos[147]= - armWWos[63] + armWWos[147];
   armWWos[171]=armWWos[147] + armWWos[38];
   armWWos[171]=armWWos[32]*armWWos[171];
   armWWos[174]= - 57./2.*armWWos[40] - 9./8.*armWWos[38] - 1./8.*
   armWWos[8] + 1./4.*armWWos[63] + 15*armWWos[61];
   armWWos[174]=armWWos[97]*armWWos[174];
   armWWos[178]=armWWos[32] + 3./2.*armWWos[97];
   armWWos[178]=armWWos[19]*armWWos[178];
   armWWos[179]=armWWos[32] - 111./4.*armWWos[97];
   armWWos[179]=armWWos[18]*armWWos[179];
   armWWos[181]= - 11*armWWos[18] - 33./16.*armWWos[19] - 155./8.*
   armWWos[40] + 3./8.*armWWos[36] + 11*armWWos[61];
   armWWos[181]=armWWos[98]*armWWos[181];
   armWWos[188]=8./3. - 5*armWWos[6];
   armWWos[189]=armWWos[1]*armWWos[188];
   armWWos[188]=armWWos[35]*armWWos[188];
   armWWos[190]= - 7./3.*armWWos[10];
   armWWos[195]=63./64.*armWWos[98] + armWWos[32] - armWWos[97];
   armWWos[195]=armWWos[16]*armWWos[195];
   armWWos[196]=17./6.*armWWos[17] - armWWos[18] + armWWos[142] - 2*
   armWWos[62] + armWWos[38];
   armWWos[196]=armWWos[100]*armWWos[196];
   armWWos[198]=armWWos[32] - 9./4.*armWWos[97];
   armWWos[198]=1./2.*armWWos[198] + armWWos[100];
   armWWos[198]=armWWos[33]*armWWos[198];
   armWWos[198]=armWWos[198] - 2./3.*armWWos[100] + 33./32.*armWWos[98]
    - armWWos[32] + 39./4.*armWWos[97];
   armWWos[198]=armWWos[15]*armWWos[198];
   armWWos[200]= - armWWos[40]*armWWos[32];
   armWWos[102]=armWWos[110] + armWWos[102] + armWWos[151] + 
   armWWos[103] + armWWos[198] + armWWos[120] + armWWos[196] + 1./4.*
   armWWos[195] + armWWos[190] + 1./3.*armWWos[188] + 1./9.*
   armWWos[189] + 1837./64.*armWWos[9] + 3./16.*armWWos[181] + 1./2.*
   armWWos[179] - armWWos[14] + 1./4.*armWWos[178] + armWWos[174] + 957.
   /64.*armWWos[13] + armWWos[200] + 1./2.*armWWos[171] - 2*armWWos[57]
    + 1681./512.*armWWos[39] + 21./8.*armWWos[30] + armWWos[49] + 145./
   512.*armWWos[60] + 207./256.*armWWos[58] - 33./8.*armWWos[27] - 165./
   64.*armWWos[52] - 657./512.*armWWos[55] + 1./4.*armWWos[168] + 
   armWWos[135];
   armWWos[102]=armWWos[177]*armWWos[102];
   armWWos[103]=armWWos[113] + 2*armWWos[46];
   armWWos[110]= - 5*armWWos[45];
   armWWos[103]= - armWWos[48] - 3*armWWos[43] + 2./3.*armWWos[103] + 
   armWWos[110];
   armWWos[113]=armWWos[60] - 307./4. - armWWos[58];
   armWWos[120]=armWWos[63] - 1./2.*armWWos[8];
   armWWos[135]=armWWos[97]*armWWos[120];
   armWWos[113]=1./6.*armWWos[113] + armWWos[135];
   armWWos[113]=armWWos[97]*armWWos[113];
   armWWos[135]=1 + armWWos[187];
   armWWos[135]=1./2.*armWWos[135] + armWWos[60];
   armWWos[135]=mmZ*armWWos[106]*armWWos[135];
   armWWos[151]= - 1./2.*armWWos[58] - 5./6.*armWWos[55] + 1./8. - 
   armWWos[56];
   armWWos[151]=11./12.*armWWos[39] + 1./2.*armWWos[151] - armWWos[30];
   armWWos[151]=armWWos[32]*armWWos[151];
   armWWos[168]=3067 + 123*armWWos[55];
   armWWos[168]=165*armWWos[13] - 123./8.*armWWos[39] - 123./8.*
   armWWos[60] + 123./4.*armWWos[58] + 1./8.*armWWos[168] - 165*
   armWWos[52];
   armWWos[168]=armWWos[98]*armWWos[168];
   armWWos[171]= - armWWos[16]*armWWos[106];
   armWWos[174]= - 2./3.*armWWos[10] - armWWos[14] + 6*armWWos[57] + 53.
   /12. + armWWos[49];
   armWWos[174]=armWWos[100]*armWWos[174];
   armWWos[178]= - 369./512.*armWWos[98] + armWWos[194];
   armWWos[178]=armWWos[33]*armWWos[178];
   armWWos[103]=armWWos[135] + armWWos[178] + armWWos[174] + 1./4.*
   armWWos[171] + 495./64.*armWWos[139] + 3./64.*armWWos[168] + 1./4.*
   armWWos[113] + 2*armWWos[103] + armWWos[151];
   armWWos[103]=mmZ*armWWos[103];
   armWWos[113]= - 85./9.*armWWos[32] - 369./32.*armWWos[98];
   armWWos[113]=armWWos[33]*armWWos[113];
   armWWos[135]= - mmZ*armWWos[106];
   armWWos[113]=3./4.*armWWos[135] + 1./2.*armWWos[113] + 495./8.*
   armWWos[160] + 7387./64.*armWWos[98] + 73./18.*armWWos[32] - 
   armWWos[97];
   armWWos[113]=armWWos[34]*armWWos[113];
   armWWos[122]=3./2.*armWWos[122] + armWWos[156] + armWWos[211] - 
   armWWos[36] + armWWos[8];
   armWWos[135]= - 11./3. - armWWos[33];
   armWWos[135]=mmZ*armWWos[135];
   armWWos[151]=2*armWWos[34];
   armWWos[135]=armWWos[151] + armWWos[122] + armWWos[135];
   armWWos[135]=armWWos[5]*armWWos[135];
   armWWos[156]= - armWWos[47] - 1./3.*armWWos[41];
   armWWos[168]=armWWos[156] + 2./3.*armWWos[46];
   armWWos[123]=42*armWWos[123] - 13./2. - armWWos[33];
   armWWos[123]=armWWos[5]*armWWos[123];
   armWWos[171]=12*armWWos[124];
   armWWos[168]=armWWos[171] + armWWos[123] + 23./6.*armWWos[43] + 2*
   armWWos[168] - 3*armWWos[45];
   armWWos[168]=mmt*armWWos[168];
   armWWos[174]=armWWos[128] + armWWos[215];
   armWWos[174]=armWWos[19]*armWWos[174];
   armWWos[178]=1./2.*armWWos[32] - armWWos[97];
   armWWos[179]= - 215./128.*armWWos[98];
   armWWos[178]=1./3.*armWWos[178] + armWWos[179];
   armWWos[178]=armWWos[16]*armWWos[178];
   armWWos[181]= - 25./3.*armWWos[42] + 28489./1536. + 8*armWWos[53];
   armWWos[187]= - 7./8.*armWWos[36] - 11*armWWos[61];
   armWWos[187]=121*armWWos[18] + 523./16.*armWWos[19] + 11*
   armWWos[187] + 1721./8.*armWWos[40];
   armWWos[187]=armWWos[98]*armWWos[187];
   armWWos[128]=armWWos[128] + 369./64.*armWWos[98];
   armWWos[128]=armWWos[16]*armWWos[128];
   armWWos[128]= - 113./9.*armWWos[33] + armWWos[128] - 12307./384. + 
   armWWos[9];
   armWWos[128]=armWWos[33]*armWWos[128];
   armWWos[126]=armWWos[63] + armWWos[126];
   armWWos[126]=1./3.*armWWos[32]*armWWos[126];
   armWWos[120]=1./3.*armWWos[120] + 1./2.*armWWos[40];
   armWWos[120]=armWWos[97]*armWWos[120];
   armWWos[188]=1./2.*armWWos[120];
   armWWos[189]= - armWWos[15]*armWWos[98];
   armWWos[102]=armWWos[102] + armWWos[168] + armWWos[135] + 1./2.*
   armWWos[113] + armWWos[103] + 121./32.*armWWos[189] + 1./4.*
   armWWos[128] + 1./2.*armWWos[178] - 1265./64.*armWWos[9] + 1./16.*
   armWWos[187] + 1./4.*armWWos[174] + armWWos[188] - 17./64.*
   armWWos[13] + 79./36.*armWWos[200] + armWWos[126] - 659./1536.*
   armWWos[39] - 19./36.*armWWos[30] + 143./512.*armWWos[60] - 4045./
   768.*armWWos[58] + 5./24.*armWWos[27] - 1093./192.*armWWos[52] + 659.
   /1536.*armWWos[55] + 1./3.*armWWos[181] + 23./4.*armWWos[50];
   armWWos[102]=armWWos[177]*armWWos[102];
   armWWos[103]= - 121./3. + 5*armWWos[55];
   armWWos[103]=5./8.*armWWos[103] + 37*armWWos[52];
   armWWos[113]=1./8.*armWWos[36] + armWWos[61];
   armWWos[113]= - 1./3.*armWWos[18] + 1./48.*armWWos[19] + 1./3.*
   armWWos[113] - 7./8.*armWWos[40];
   armWWos[113]=armWWos[98]*armWWos[113];
   armWWos[128]= - armWWos[16]*armWWos[98];
   armWWos[135]= - 1./2. + armWWos[148];
   armWWos[135]=armWWos[33]*armWWos[135];
   armWWos[103]=25./12.*armWWos[159] + 25./32.*armWWos[135] + 125./96.*
   armWWos[128] - 277./72.*armWWos[9] + 25./2.*armWWos[113] - 151./24.*
   armWWos[13] - 25./64.*armWWos[39] - 25./64.*armWWos[60] - 565./96.*
   armWWos[58] + 1./8.*armWWos[103] + 5./3.*armWWos[27];
   armWWos[113]=3 + 1./3.*armWWos[55];
   armWWos[113]=1./9.*armWWos[13] - 1./24.*armWWos[39] - 1./24.*
   armWWos[60] + 1./12.*armWWos[58] + 1./8.*armWWos[113] - 1./9.*
   armWWos[52];
   armWWos[113]=armWWos[98]*armWWos[113];
   armWWos[128]= - armWWos[33]*armWWos[98];
   armWWos[113]=1./24.*armWWos[128] + armWWos[113] + 1./9.*armWWos[139]
   ;
   armWWos[113]=mmZ*armWWos[113];
   armWWos[135]=1./8.*armWWos[128] - 5./8.*armWWos[98] + 1./3.*
   armWWos[160];
   armWWos[135]=armWWos[34]*armWWos[135];
   armWWos[103]=25./6.*armWWos[135] + 1./3.*armWWos[103] + 25./8.*
   armWWos[113];
   armWWos[103]=armWWos[121]*armWWos[103];
   armWWos[113]= - 4819./9. + 163*armWWos[55];
   armWWos[113]=1./8.*armWWos[113] + 1733./9.*armWWos[52];
   armWWos[135]=armWWos[32]*armWWos[147];
   armWWos[113]=1./4.*armWWos[120] - 2437./288.*armWWos[13] + 11./4.*
   armWWos[200] + 11./2.*armWWos[135] - 163./256.*armWWos[39] + 5./24.*
   armWWos[30] - 2057./768.*armWWos[60] - 151./384.*armWWos[58] + 1./32.
   *armWWos[113] + 22./9.*armWWos[27];
   armWWos[120]= - 1 - armWWos[58];
   armWWos[135]=armWWos[32]*armWWos[120];
   armWWos[147]=1 + armWWos[58];
   armWWos[147]=armWWos[97]*armWWos[147];
   armWWos[135]=1./12.*armWWos[147] + armWWos[116] + 11./2.*
   armWWos[135];
   armWWos[147]=31 + 7./3.*armWWos[55];
   armWWos[147]=37./9.*armWWos[13] - 7./24.*armWWos[39] - 7./24.*
   armWWos[60] + 7./12.*armWWos[58] + 1./8.*armWWos[147] - 37./9.*
   armWWos[52];
   armWWos[147]=armWWos[98]*armWWos[147];
   armWWos[135]=35./384.*armWWos[128] + 185./144.*armWWos[139] + 1./3.*
   armWWos[135] + 5./16.*armWWos[147];
   armWWos[135]=mmZ*armWWos[135];
   armWWos[128]=35./32.*armWWos[128] + 185./12.*armWWos[160] + 
   armWWos[145] - 3025./96.*armWWos[98];
   armWWos[128]=armWWos[34]*armWWos[128];
   armWWos[145]=121./6. - armWWos[33];
   armWWos[132]=1./2.*armWWos[145] + armWWos[132];
   armWWos[132]=armWWos[5]*armWWos[132];
   armWWos[145]=2*armWWos[43] + armWWos[116];
   armWWos[124]=6*armWWos[124] + 1./3.*armWWos[145] + armWWos[132];
   armWWos[124]=mmt*armWWos[124];
   armWWos[132]=35./32.*armWWos[148] - 35./64. - armWWos[9];
   armWWos[145]=1./3.*armWWos[33];
   armWWos[132]=1./2.*armWWos[132] + armWWos[145];
   armWWos[132]=armWWos[33]*armWWos[132];
   armWWos[147]=17*armWWos[36] + 5*armWWos[8];
   armWWos[133]=55./6.*mmZ - 55./6.*armWWos[16] + armWWos[133] + 1./2.*
   armWWos[147] - 17*armWWos[40];
   armWWos[133]=1./3.*armWWos[133] + armWWos[34];
   armWWos[133]=armWWos[5]*armWWos[133];
   armWWos[147]= - 103./8.*armWWos[36] + 47*armWWos[61];
   armWWos[147]= - 47./9.*armWWos[18] + 227./144.*armWWos[19] + 1./9.*
   armWWos[147] - 63./8.*armWWos[40];
   armWWos[147]=armWWos[98]*armWWos[147];
   armWWos[148]=925./192.*armWWos[98] + armWWos[140] - 1./6.*
   armWWos[97];
   armWWos[148]=armWWos[16]*armWWos[148];
   armWWos[103]=1./8.*armWWos[103] + armWWos[124] + armWWos[133] + 1./6.
   *armWWos[128] + 1./2.*armWWos[135] + 235./144.*armWWos[159] + 1./6.*
   armWWos[132] + 1./6.*armWWos[148] - 5911./2592.*armWWos[9] + 1./3.*
   armWWos[113] + 5./8.*armWWos[147];
   armWWos[103]=armWWos[121]*armWWos[103];
   armWWos[113]= - 11*armWWos[39] + 13*armWWos[58] + 19./4. + 11*
   armWWos[55];
   armWWos[113]=armWWos[32]*armWWos[113];
   armWWos[113]=armWWos[116] + 1./4.*armWWos[113];
   armWWos[116]= - 1./9.*armWWos[60] + 5./12. + armWWos[58];
   armWWos[115]=1./2.*armWWos[116] + 1./3.*armWWos[115];
   armWWos[115]=armWWos[97]*armWWos[115];
   armWWos[116]=mmZ*armWWos[106]*armWWos[120];
   armWWos[120]= - 3343 - 941./3.*armWWos[55];
   armWWos[120]= - 2201./9.*armWWos[13] + 941./24.*armWWos[39] + 941./
   24.*armWWos[60] - 941./12.*armWWos[58] + 1./8.*armWWos[120] + 2201./
   9.*armWWos[52];
   armWWos[120]=armWWos[98]*armWWos[120];
   armWWos[124]=armWWos[33]*armWWos[98];
   armWWos[113]=1./24.*armWWos[116] + 941./1536.*armWWos[124] + 1./12.*
   armWWos[155] + 2201./576.*armWWos[160] + 1./64.*armWWos[120] + 1./3.
   *armWWos[113] + 1./4.*armWWos[115];
   armWWos[113]=mmZ*armWWos[113];
   armWWos[115]= - 3*armWWos[9];
   armWWos[116]=armWWos[140] - 941./64.*armWWos[98];
   armWWos[116]=armWWos[16]*armWWos[116];
   armWWos[116]=armWWos[145] + 1./3.*armWWos[116] + 455./1152. + 
   armWWos[115];
   armWWos[116]=armWWos[33]*armWWos[116];
   armWWos[120]=armWWos[152] + 941./32.*armWWos[98];
   armWWos[120]=armWWos[33]*armWWos[120];
   armWWos[120]=1./4.*armWWos[173] + 1./6.*armWWos[120] + 2201./72.*
   armWWos[139] - 17485./576.*armWWos[98] + 7./6.*armWWos[32] - 
   armWWos[97];
   armWWos[120]=armWWos[34]*armWWos[120];
   armWWos[124]= - 4./3. + armWWos[33];
   armWWos[124]=mmZ*armWWos[124];
   armWWos[122]=armWWos[151] + armWWos[122] + armWWos[124];
   armWWos[122]=armWWos[5]*armWWos[122];
   armWWos[110]=armWWos[146] + armWWos[110] + 11./2.*armWWos[43];
   armWWos[110]=armWWos[171] + 1./3.*armWWos[110] + armWWos[123];
   armWWos[110]=mmt*armWWos[110];
   armWWos[123]=71./1152.*armWWos[39] + 103./27.*armWWos[30] + 143./128.
   *armWWos[60] - 4711./576.*armWWos[58] - 179./54.*armWWos[27] + 2963./
   432.*armWWos[52] - 71./1152.*armWWos[55] - 258031./10368. + 3*
   armWWos[50];
   armWWos[124]=armWWos[140] + 1./2.*armWWos[97];
   armWWos[124]=armWWos[19]*armWWos[124];
   armWWos[128]=107*armWWos[18] - 2527./16.*armWWos[19] + 1067./8.*
   armWWos[40] + 793./8.*armWWos[36] - 107*armWWos[61];
   armWWos[128]=armWWos[98]*armWWos[128];
   armWWos[132]=armWWos[179] - 5./2.*armWWos[32] - 1./3.*armWWos[97];
   armWWos[132]=armWWos[16]*armWWos[132];
   armWWos[103]=armWWos[103] + armWWos[110] + armWWos[122] + 1./2.*
   armWWos[120] + armWWos[113] + 107./96.*armWWos[189] + 1./4.*
   armWWos[116] + 1./2.*armWWos[132] + 5183./5184.*armWWos[9] + 1./48.*
   armWWos[128] + 1./12.*armWWos[124] + armWWos[188] - 2827./1728.*
   armWWos[13] + 3./4.*armWWos[200] + 1./4.*armWWos[123] + armWWos[126]
   ;
   armWWos[103]=armWWos[121]*armWWos[103];
   armWWos[110]=1./3.*armWWos[48] + armWWos[43] + 2*armWWos[45] + 
   armWWos[156] + armWWos[114];
   armWWos[113]= - 8./3.*armWWos[13];
   armWWos[114]=armWWos[113] - 7 + 8./3.*armWWos[52];
   armWWos[114]=armWWos[98]*armWWos[114];
   armWWos[110]=8./3.*armWWos[160] + 2./3.*armWWos[110] + armWWos[114];
   armWWos[114]= - 1./8. - 1./3.*armWWos[60];
   armWWos[106]=mmZ*armWWos[106]*armWWos[114];
   armWWos[106]=2*armWWos[110] + armWWos[106];
   armWWos[106]=mmZ*armWWos[106];
   armWWos[108]=14./9.*armWWos[9] - 22./3.*armWWos[13] + armWWos[138]
    + 2./3.*armWWos[27] + 20./3.*armWWos[52] - 53./9. + armWWos[108];
   armWWos[110]= - armWWos[98] + armWWos[139];
   armWWos[110]=armWWos[34]*armWWos[110];
   armWWos[114]= - 5./3. + armWWos[33];
   armWWos[114]=armWWos[33]*armWWos[114];
   armWWos[114]=4./9. + armWWos[114];
   armWWos[114]=armWWos[31]*armWWos[114]*pow(armWWos[3],4);
   armWWos[102]=armWWos[114] + armWWos[103] + armWWos[102] + 64./3.*
   armWWos[110] + 4./3.*armWWos[108] + armWWos[106];
   armWWos[102]=armWWos[31]*armWWos[102];
   armWWos[103]=337./6. + armWWos[164];
   armWWos[103]=1./2.*armWWos[103] + armWWos[153];
   armWWos[103]=armWWos[17]*armWWos[103];
   armWWos[106]=149./2. - 25*armWWos[10];
   armWWos[106]=armWWos[16]*armWWos[106];
   armWWos[108]=3./4.*armWWos[95];
   armWWos[110]=1./4.*armWWos[68];
   armWWos[114]= - 2*armWWos[19];
   armWWos[103]=1./6.*armWWos[106] + armWWos[103] - 107./12.*
   armWWos[18] + 35./3.*armWWos[20] + armWWos[114] + armWWos[110] + 
   armWWos[108] + 11*armWWos[163] - 47./12.*armWWos[67];
   armWWos[103]=armWWos[100]*armWWos[103];
   armWWos[106]= - 33./2.*armWWos[9];
   armWWos[116]=armWWos[182] + 79./3. + armWWos[106];
   armWWos[116]=armWWos[16]*armWWos[116];
   armWWos[120]=armWWos[15]*armWWos[16];
   armWWos[120]=1./2.*armWWos[149] + armWWos[120];
   armWWos[122]=armWWos[16] - armWWos[15];
   armWWos[122]=7./2.*armWWos[122] - 5*mmZ;
   armWWos[122]=mmZ*armWWos[122];
   armWWos[122]=3*armWWos[120] + armWWos[122];
   armWWos[122]=armWWos[5]*armWWos[122];
   armWWos[123]=11*armWWos[66] + armWWos[154];
   armWWos[123]=1./2.*armWWos[123] - 7*armWWos[19];
   armWWos[124]=37 + armWWos[104];
   armWWos[124]=armWWos[17]*armWWos[124];
   armWWos[126]= - 419./6. + armWWos[105];
   armWWos[126]=armWWos[15]*armWWos[126];
   armWWos[128]=99*armWWos[9];
   armWWos[132]=armWWos[10] + 3535./9. + armWWos[128];
   armWWos[132]=mmZ*armWWos[132];
   armWWos[116]=3*armWWos[122] + 1./4.*armWWos[132] + 1./2.*
   armWWos[126] + 1./4.*armWWos[116] + 1./8.*armWWos[124] - 33*
   armWWos[18] + 3*armWWos[123] - 217./24.*armWWos[20];
   armWWos[116]=armWWos[5]*armWWos[116];
   armWWos[122]=11*armWWos[87];
   armWWos[123]= - 11*armWWos[91];
   armWWos[124]=7*armWWos[85] + armWWos[123] - 7 + armWWos[122];
   armWWos[126]=11*armWWos[89];
   armWWos[132]= - 11./2.*armWWos[13];
   armWWos[107]=armWWos[107] + armWWos[132] + 1./2.*armWWos[124] + 
   armWWos[126];
   armWWos[107]=1./2.*armWWos[99]*armWWos[107];
   armWWos[124]= - 28./9. - 11./4.*armWWos[9];
   armWWos[124]=5*armWWos[124] + armWWos[10];
   armWWos[124]=armWWos[10]*armWWos[124];
   armWWos[124]=armWWos[124] + 55./2.*armWWos[9] + 40*armWWos[14] + 55./
   4.*armWWos[13] - 110*armWWos[89] + armWWos[170] - 11./2.*armWWos[91]
    + 2*armWWos[80] - 55./4.*armWWos[81] - 55./4.*armWWos[82] - 3*
   armWWos[86] - 23041./144. - 9*armWWos[84];
   armWWos[124]=armWWos[100]*armWWos[124];
   armWWos[133]= - 2*armWWos[79] - armWWos[76];
   armWWos[133]=11*armWWos[133] + 357./16.*armWWos[72];
   armWWos[135]= - armWWos[9]*armWWos[99];
   armWWos[139]=11./4.*armWWos[135];
   armWWos[140]= - armWWos[10]*armWWos[99];
   armWWos[145]=7./4.*armWWos[140];
   armWWos[124]=armWWos[124] + armWWos[145] + armWWos[139] + 
   armWWos[107] - 55./2.*armWWos[71] - 11*armWWos[74] - 11./2.*
   armWWos[96] + armWWos[70] + 3*armWWos[133] + 8*armWWos[75];
   armWWos[124]=mmZ*armWWos[124];
   armWWos[133]=armWWos[158] - 5./3.*armWWos[75];
   armWWos[133]=14./3.*armWWos[71] + armWWos[197] + 4*armWWos[133] + 
   armWWos[167];
   armWWos[133]=mmH*armWWos[133];
   armWWos[146]= - 47./6.*armWWos[10] + 353./24. + 11*armWWos[9];
   armWWos[146]=armWWos[100]*armWWos[146];
   armWWos[147]= - 15./8.*armWWos[99];
   armWWos[146]=armWWos[147] + armWWos[146];
   armWWos[146]=armWWos[15]*armWWos[146];
   armWWos[148]=11*armWWos[135];
   armWWos[151]=armWWos[10]*armWWos[99];
   armWWos[152]=7*armWWos[151] + 119./3.*armWWos[99] + armWWos[148];
   armWWos[152]=1./4.*armWWos[16]*armWWos[152];
   armWWos[153]= - 488245./864. - 57*armWWos[92];
   armWWos[154]=9491./3. + 511*armWWos[9];
   armWWos[154]=armWWos[9]*armWWos[154];
   armWWos[155]= - 53 - 55./6.*armWWos[9];
   armWWos[155]=armWWos[10]*armWWos[155];
   armWWos[156]=1./6.*armWWos[86];
   armWWos[158]=15./2.*armWWos[18] + 5./12.*armWWos[20] + 29./6.*
   armWWos[19] - 15./2.*armWWos[95] + 7./3.*armWWos[68];
   armWWos[158]=1./2.*armWWos[99]*armWWos[158];
   armWWos[159]=armWWos[9]*armWWos[99];
   armWWos[160]=11*armWWos[159];
   armWWos[164]= - 61./3.*armWWos[99] + armWWos[160];
   armWWos[164]=1./8.*armWWos[17]*armWWos[164];
   armWWos[103]=armWWos[133] + armWWos[116] + armWWos[124] + 
   armWWos[146] + armWWos[103] + armWWos[152] + armWWos[164] + 1./2.*
   armWWos[155] + 1./16.*armWWos[154] + armWWos[158] + 121./12.*
   armWWos[14] - 3709./48.*armWWos[13] + 922./9.*armWWos[30] - 55./4.*
   armWWos[89] - 19./6.*armWWos[85] - 51./4.*armWWos[91] + 11./12.*
   armWWos[87] + 11./12.*armWWos[81] + 11./12.*armWWos[82] + 
   armWWos[156] + 1285./16.*armWWos[88] - 215./8.*armWWos[83] + 1./2.*
   armWWos[153] - 8./3.*armWWos[90];
   armWWos[103]=armWWos[65]*armWWos[103];
   armWWos[116]=5*mmZ;
   armWWos[124]=7*armWWos[169] + armWWos[116];
   armWWos[124]=mmZ*armWWos[124];
   armWWos[133]=armWWos[16] + armWWos[15];
   armWWos[146]=armWWos[15]*armWWos[133];
   armWWos[124]=1./2.*armWWos[124] + 1./4.*armWWos[149] + armWWos[146];
   armWWos[124]=armWWos[5]*armWWos[124];
   armWWos[146]=1./2.*armWWos[68];
   armWWos[153]=armWWos[146] - 11*armWWos[66] + armWWos[67];
   armWWos[153]=1./2.*armWWos[153] + 5*armWWos[19];
   armWWos[153]=5*armWWos[18] + 1./2.*armWWos[153] + armWWos[20];
   armWWos[128]= - 85 + armWWos[128];
   armWWos[128]=1./8.*armWWos[128] - armWWos[10];
   armWWos[128]=armWWos[17]*armWWos[128];
   armWWos[154]= - 1 - 11./2.*armWWos[9];
   armWWos[154]=3*armWWos[154] + armWWos[10];
   armWWos[154]=armWWos[16]*armWWos[154];
   armWWos[155]=armWWos[193] - 4 - 99./4.*armWWos[9];
   armWWos[155]=armWWos[15]*armWWos[155];
   armWWos[131]=armWWos[131] - 1847./6. - 99*armWWos[9];
   armWWos[131]=mmZ*armWWos[131];
   armWWos[124]=9*armWWos[124] + 1./4.*armWWos[131] + armWWos[155] + 3./
   4.*armWWos[154] + 9*armWWos[153] + armWWos[128];
   armWWos[124]=armWWos[5]*armWWos[124];
   armWWos[122]=armWWos[85] + armWWos[123] - 5./2. + armWWos[122];
   armWWos[122]=armWWos[125] + armWWos[132] + 1./2.*armWWos[122] + 
   armWWos[126];
   armWWos[122]=armWWos[99]*armWWos[122];
   armWWos[123]= - armWWos[10] + 62./9. + 33./4.*armWWos[9];
   armWWos[123]=armWWos[10]*armWWos[123];
   armWWos[123]=armWWos[123] + armWWos[106] - 19*armWWos[14] - 33./4.*
   armWWos[13] + 66*armWWos[89] + 3./4.*armWWos[85] + 33./4.*
   armWWos[91] - 2*armWWos[80] + 33./4.*armWWos[81] + 33./4.*
   armWWos[82] + 3*armWWos[86] + 13903./144. + armWWos[202];
   armWWos[123]=armWWos[100]*armWWos[123];
   armWWos[122]=armWWos[123] + 3./4.*armWWos[140] + 33./4.*armWWos[135]
    + 3./2.*armWWos[122] + 33./2.*armWWos[71] + 33./2.*armWWos[74] + 33.
   /4.*armWWos[96] - 459./16.*armWWos[72] - armWWos[70];
   armWWos[122]=mmZ*armWWos[122];
   armWWos[123]=armWWos[151] + 13*armWWos[99] + armWWos[148];
   armWWos[123]=armWWos[16]*armWWos[123];
   armWWos[105]= - 1739./9. + armWWos[105];
   armWWos[105]=1./2.*armWWos[105] + armWWos[190];
   armWWos[105]=armWWos[17]*armWWos[105];
   armWWos[126]= - 13./2. + armWWos[10];
   armWWos[126]=armWWos[16]*armWWos[126];
   armWWos[105]=3*armWWos[126] + 1./2.*armWWos[105] - 145./12.*
   armWWos[18] - 661./24.*armWWos[20] - 1./4.*armWWos[19] + 3./8.*
   armWWos[68] - 1./4.*armWWos[95] + 11./3.*armWWos[67] - 33./4.*
   armWWos[66] + 33./2.*armWWos[94] + 9./2.*armWWos[69] + 5*armWWos[93]
   ;
   armWWos[105]=armWWos[100]*armWWos[105];
   armWWos[126]=armWWos[175] + 21./4.*armWWos[20] - 31./2.*armWWos[19]
    - 7./2.*armWWos[95] + 3*armWWos[68];
   armWWos[126]=armWWos[99]*armWWos[126];
   armWWos[128]= - 1697./3. + armWWos[130];
   armWWos[128]=armWWos[9]*armWWos[128];
   armWWos[130]= - 1./2.*armWWos[10] + 41./9. + armWWos[134];
   armWWos[130]=armWWos[10]*armWWos[130];
   armWWos[131]= - 15*armWWos[99] + armWWos[160];
   armWWos[131]=armWWos[17]*armWWos[131];
   armWWos[104]= - 697./36. + armWWos[104];
   armWWos[104]=1./2.*armWWos[104] + 22./3.*armWWos[10];
   armWWos[104]=armWWos[100]*armWWos[104];
   armWWos[104]= - 7./8.*armWWos[99] + armWWos[104];
   armWWos[104]=armWWos[15]*armWWos[104];
   armWWos[132]= - 13./2.*armWWos[71] - 13./2.*armWWos[74] - 13./4.*
   armWWos[96] - 3./2.*armWWos[73] + 2./3.*armWWos[70];
   armWWos[132]=mmH*armWWos[132];
   armWWos[134]=3./2.*armWWos[85];
   armWWos[104]=armWWos[132] + armWWos[124] + armWWos[122] + 
   armWWos[104] + armWWos[105] + 3./4.*armWWos[123] + 3./8.*
   armWWos[131] + armWWos[130] + 1./16.*armWWos[128] + 1./2.*
   armWWos[126] - 7./6.*armWWos[14] - 391./8.*armWWos[13] + 33*
   armWWos[89] + armWWos[134] - 11./4.*armWWos[87] - 11./4.*armWWos[81]
    - 11./4.*armWWos[82] + armWWos[156] + 435./16.*armWWos[88] + 435./
   16.*armWWos[83] + 5416./81. + 9./4.*armWWos[84];
   armWWos[104]=armWWos[65]*armWWos[104];
   armWWos[105]= - armWWos[19] + 2*armWWos[7] + armWWos[8];
   armWWos[109]=armWWos[109] + 3*armWWos[105] + armWWos[180];
   armWWos[109]=armWWos[35]*armWWos[109];
   armWWos[122]=1 - armWWos[6];
   armWWos[123]=armWWos[1]*armWWos[122];
   armWWos[122]=armWWos[35]*armWWos[122];
   armWWos[122]=1./3.*armWWos[123] + armWWos[122];
   armWWos[122]=armWWos[17]*armWWos[122];
   armWWos[123]= - 3 + armWWos[6];
   armWWos[124]=armWWos[1]*armWWos[123];
   armWWos[123]=armWWos[35]*armWWos[123];
   armWWos[123]=armWWos[124] + 3*armWWos[123];
   armWWos[124]=armWWos[16]*armWWos[123];
   armWWos[123]=armWWos[15]*armWWos[123];
   armWWos[126]=7./6. + armWWos[6];
   armWWos[126]=armWWos[1]*armWWos[126];
   armWWos[128]=7./2. + armWWos[199];
   armWWos[128]=armWWos[35]*armWWos[128];
   armWWos[126]=armWWos[126] + armWWos[128];
   armWWos[126]=mmZ*armWWos[126];
   armWWos[105]= - 2*armWWos[18] + armWWos[105] - 1./6.*armWWos[20];
   armWWos[105]=armWWos[1]*armWWos[105];
   armWWos[105]=armWWos[126] + armWWos[123] + 1./2.*armWWos[124] + 1./2.
   *armWWos[122] + armWWos[105] + armWWos[109];
   armWWos[105]=armWWos[5]*armWWos[105];
   armWWos[109]=armWWos[183] + armWWos[136] + armWWos[19];
   armWWos[109]=armWWos[99]*armWWos[109];
   armWWos[122]=armWWos[16]*armWWos[99];
   armWWos[123]=armWWos[176] - 1./2.*armWWos[18] + 1./2.*armWWos[95] - 
   armWWos[19];
   armWWos[123]=armWWos[100]*armWWos[123];
   armWWos[124]= - armWWos[99] + armWWos[100];
   armWWos[126]=armWWos[15]*armWWos[124];
   armWWos[124]=mmZ*armWWos[124];
   armWWos[109]=3./16.*armWWos[124] + 1./4.*armWWos[126] + armWWos[123]
    + armWWos[109] + 3./2.*armWWos[122];
   armWWos[109]=armWWos[177]*armWWos[65]*armWWos[109];
   armWWos[122]= - 5./3. + armWWos[6];
   armWWos[122]=armWWos[6]*armWWos[122];
   armWWos[123]=pow(Pi,2);
   armWWos[122]=armWWos[122] + 4./9. - armWWos[123];
   armWWos[124]=armWWos[1]*armWWos[122];
   armWWos[126]=5./3.*armWWos[123] + 37./3.*armWWos[29] + 11053./648.
    - 11*armWWos[26];
   armWWos[128]= - 1./3.*armWWos[14];
   armWWos[130]= - 7./9. - 1./4.*armWWos[6];
   armWWos[130]=armWWos[6]*armWWos[130];
   armWWos[126]=1./9.*armWWos[124] + 121./12.*armWWos[9] + 5./3.*
   armWWos[130] + armWWos[128] + 11./2.*armWWos[13] + 19./4.*
   armWWos[30] - 11./4.*armWWos[27] + 23./72.*armWWos[12] + 1./4.*
   armWWos[126] + 1./3.*armWWos[25];
   armWWos[126]=armWWos[1]*armWWos[126];
   armWWos[122]=armWWos[35]*armWWos[122];
   armWWos[131]=5*armWWos[123] + 37*armWWos[29] + 11053./216. - 33*
   armWWos[26];
   armWWos[122]=armWWos[122] + 2./3.*armWWos[124] + 121./4.*armWWos[9]
    + 5*armWWos[130] - armWWos[14] + 33./2.*armWWos[13] + 57./4.*
   armWWos[30] - 33./4.*armWWos[27] + 23./24.*armWWos[12] + 1./4.*
   armWWos[131] + armWWos[25];
   armWWos[122]=armWWos[35]*armWWos[122];
   armWWos[124]= - 1./3. - armWWos[6];
   armWWos[130]=armWWos[1]*armWWos[124];
   armWWos[124]=armWWos[35]*armWWos[124];
   armWWos[124]=1./3.*armWWos[130] + armWWos[124];
   armWWos[124]=armWWos[17]*armWWos[124];
   armWWos[130]= - armWWos[18] + armWWos[7] + armWWos[142];
   armWWos[131]=armWWos[1]*armWWos[130];
   armWWos[130]=armWWos[35]*armWWos[130];
   armWWos[124]=1./2.*armWWos[124] + 1./3.*armWWos[131] + armWWos[130];
   armWWos[124]=armWWos[100]*armWWos[124];
   armWWos[130]= - armWWos[25] + 7./3. + 2*armWWos[28];
   armWWos[131]=1./3.*armWWos[6] + armWWos[172] + 1./3.*armWWos[130] + 
   armWWos[12];
   armWWos[131]=armWWos[1]*armWWos[131];
   armWWos[132]=2./3. - armWWos[6];
   armWWos[136]=armWWos[1]*armWWos[132];
   armWWos[132]=armWWos[35]*armWWos[132];
   armWWos[132]=1./3.*armWWos[136] + armWWos[132];
   armWWos[132]=armWWos[10]*armWWos[132];
   armWWos[130]=armWWos[6] + armWWos[14] + armWWos[130] + 3*armWWos[12]
   ;
   armWWos[130]=armWWos[35]*armWWos[130];
   armWWos[130]=armWWos[132] + armWWos[131] + armWWos[130];
   armWWos[130]=armWWos[100]*armWWos[130];
   armWWos[131]=armWWos[22] + armWWos[21];
   armWWos[132]=armWWos[131] + 1./3.*armWWos[24];
   armWWos[132]=armWWos[1]*armWWos[132];
   armWWos[131]=3*armWWos[131] + armWWos[24];
   armWWos[131]=armWWos[35]*armWWos[131];
   armWWos[131]=armWWos[132] + armWWos[131];
   armWWos[130]=2*armWWos[131] + armWWos[130];
   armWWos[130]=mmZ*armWWos[130];
   armWWos[131]= - 2./3. + armWWos[6];
   armWWos[132]=armWWos[1]*armWWos[131];
   armWWos[131]=armWWos[35]*armWWos[131];
   armWWos[131]=1./3.*armWWos[132] + armWWos[131];
   armWWos[132]=armWWos[15]*armWWos[100]*armWWos[131];
   armWWos[136]= - 1./3.*armWWos[1] - armWWos[35];
   armWWos[136]=armWWos[10]*armWWos[136];
   armWWos[104]=3./2.*armWWos[109] + armWWos[104] + armWWos[105] + 
   armWWos[130] + armWWos[132] + armWWos[124] + 7./3.*armWWos[136] + 
   armWWos[126] + armWWos[122];
   armWWos[104]=armWWos[177]*armWWos[104];
   armWWos[105]=37./2.*armWWos[11];
   armWWos[109]= - 7*armWWos[29];
   armWWos[122]=7*armWWos[123];
   armWWos[124]=armWWos[122] + armWWos[109] + armWWos[105] + 1109./36.
    + 23*armWWos[26];
   armWWos[124]=1./3.*armWWos[124] + 3*armWWos[27];
   armWWos[126]= - 7./3. + armWWos[127];
   armWWos[126]=armWWos[6]*armWWos[126];
   armWWos[130]= - 55 + armWWos[166];
   armWWos[130]=armWWos[9]*armWWos[130];
   armWWos[113]=1./6.*armWWos[130] + 7./6.*armWWos[126] + armWWos[113]
    + 1./4.*armWWos[124] - 10./9.*armWWos[30];
   armWWos[113]=armWWos[1]*armWWos[113];
   armWWos[124]=59./36. + armWWos[26];
   armWWos[105]=9*armWWos[27] + armWWos[122] + armWWos[109] + 23*
   armWWos[124] + armWWos[105];
   armWWos[109]= - 15 - 7./2.*armWWos[6];
   armWWos[109]=armWWos[6]*armWWos[109];
   armWWos[105]=1./2.*armWWos[130] + 1./2.*armWWos[109] - 8*armWWos[13]
    + 1./4.*armWWos[105] - 10./3.*armWWos[30];
   armWWos[105]=armWWos[35]*armWWos[105];
   armWWos[109]= - 2*armWWos[22];
   armWWos[122]= - 2*armWWos[21];
   armWWos[124]=armWWos[122] + armWWos[23] + armWWos[109];
   armWWos[124]=2*armWWos[124] - armWWos[24];
   armWWos[124]=armWWos[1]*armWWos[124];
   armWWos[109]=armWWos[122] + armWWos[109] + 1./9.*armWWos[64] + 
   armWWos[23];
   armWWos[109]=2*armWWos[109] - armWWos[24];
   armWWos[109]=armWWos[35]*armWWos[109];
   armWWos[109]=1./3.*armWWos[124] + armWWos[109];
   armWWos[122]=armWWos[25] - 7./3. - 2*armWWos[28];
   armWWos[124]= - 1./3.*armWWos[6] + armWWos[128] + 1./3.*armWWos[122]
    - armWWos[12];
   armWWos[124]=armWWos[1]*armWWos[124];
   armWWos[126]=armWWos[10]*armWWos[131];
   armWWos[122]= - armWWos[6] - armWWos[14] + armWWos[122] - 3*
   armWWos[12];
   armWWos[122]=armWWos[35]*armWWos[122];
   armWWos[122]=armWWos[126] + armWWos[124] + armWWos[122];
   armWWos[122]=armWWos[100]*armWWos[122];
   armWWos[109]=2*armWWos[109] + armWWos[122];
   armWWos[109]=mmZ*armWWos[109];
   armWWos[122]= - 11./3. - armWWos[6];
   armWWos[124]=armWWos[1]*armWWos[122];
   armWWos[122]=armWWos[35]*armWWos[122];
   armWWos[122]=1./3.*armWWos[124] + armWWos[122];
   armWWos[122]=mmZ*armWWos[122];
   armWWos[124]= - 2./3. + armWWos[166];
   armWWos[124]=armWWos[1]*armWWos[124];
   armWWos[126]= - 2 + 3./2.*armWWos[6];
   armWWos[126]=armWWos[35]*armWWos[126];
   armWWos[124]=armWWos[124] + armWWos[126];
   armWWos[124]=armWWos[16]*armWWos[124];
   armWWos[122]=armWWos[124] + armWWos[122];
   armWWos[122]=armWWos[5]*armWWos[122];
   armWWos[103]=armWWos[104] + armWWos[103] + armWWos[122] + 
   armWWos[109] + armWWos[113] + armWWos[105];
   armWWos[103]=armWWos[177]*armWWos[103];
   armWWos[104]=armWWos[192] + armWWos[129] + 23./8.*armWWos[72] + 2*
   armWWos[79] - armWWos[76];
   armWWos[105]=3*armWWos[91] - 31./4. + armWWos[112];
   armWWos[105]= - armWWos[14] + 3./2.*armWWos[13] - 3*armWWos[89] + 1./
   2.*armWWos[105] + armWWos[85];
   armWWos[105]=armWWos[99]*armWWos[105];
   armWWos[109]=2 - 1./4.*armWWos[91];
   armWWos[109]=armWWos[100]*armWWos[109];
   armWWos[104]=1./3.*armWWos[109] + armWWos[140] + 3./2.*armWWos[159]
    + 1./3.*armWWos[104] + armWWos[105];
   armWWos[104]=mmZ*armWWos[104];
   armWWos[105]=armWWos[66] + 9./2.*armWWos[68];
   armWWos[109]=13./2. + armWWos[150];
   armWWos[109]=armWWos[17]*armWWos[109];
   armWWos[105]=1./3.*armWWos[109] - armWWos[18] - 47./12.*armWWos[20]
    + 1./2.*armWWos[105] - 5*armWWos[19];
   armWWos[109]=21./2.*armWWos[16] + armWWos[116];
   armWWos[109]=mmZ*armWWos[109];
   armWWos[109]=9./2.*armWWos[149] + armWWos[109];
   armWWos[109]=armWWos[5]*armWWos[109];
   armWWos[112]= - 4./3. + 9./4.*armWWos[9];
   armWWos[112]=armWWos[16]*armWWos[112];
   armWWos[113]=19./3. + armWWos[9];
   armWWos[113]=armWWos[15]*armWWos[113];
   armWWos[122]=7./4. + armWWos[184];
   armWWos[122]=mmZ*armWWos[122];
   armWWos[105]=1./2.*armWWos[109] + 1./2.*armWWos[122] + 1./4.*
   armWWos[113] + 1./2.*armWWos[105] + armWWos[112];
   armWWos[105]=armWWos[5]*armWWos[105];
   armWWos[109]=25./3.*armWWos[99] + 3*armWWos[159];
   armWWos[109]=1./2.*armWWos[109] + armWWos[151];
   armWWos[109]=armWWos[16]*armWWos[109];
   armWWos[112]= - armWWos[17]*armWWos[9];
   armWWos[112]= - 7./6.*armWWos[16] + 1./6.*armWWos[112] - 3./4.*
   armWWos[20] - 3./2.*armWWos[19] + 1./3.*armWWos[163] + 3./4.*
   armWWos[68];
   armWWos[112]=armWWos[100]*armWWos[112];
   armWWos[113]=1./36.*armWWos[88] - 33./8.*armWWos[83] + 10069./288.
    + 7*armWWos[92];
   armWWos[122]=2*armWWos[18] - 5./12.*armWWos[20] + 37./6.*armWWos[19]
    - 2*armWWos[95] - 1./3.*armWWos[68];
   armWWos[122]=armWWos[99]*armWWos[122];
   armWWos[126]=233./3. - 121./2.*armWWos[9];
   armWWos[126]=armWWos[9]*armWWos[126];
   armWWos[128]=1 + armWWos[162];
   armWWos[128]=armWWos[10]*armWWos[128];
   armWWos[129]=13./3.*armWWos[99] + 3*armWWos[135];
   armWWos[129]=armWWos[17]*armWWos[129];
   armWWos[130]=4 + armWWos[119];
   armWWos[130]=armWWos[100]*armWWos[130];
   armWWos[130]= - armWWos[99] + 1./3.*armWWos[130];
   armWWos[130]=armWWos[15]*armWWos[130];
   armWWos[131]=1./2.*armWWos[161] - 1./3.*armWWos[71];
   armWWos[131]=mmH*armWWos[131];
   armWWos[104]=1./2.*armWWos[131] + armWWos[105] + armWWos[104] + 
   armWWos[130] + 1./2.*armWWos[112] + armWWos[109] + 1./4.*
   armWWos[129] + 1./6.*armWWos[128] + 1./36.*armWWos[126] + 
   armWWos[122] + 5./12.*armWWos[14] + 151./144.*armWWos[13] + 
   armWWos[138] - 7./4.*armWWos[89] - 5./12.*armWWos[85] - 79./36.*
   armWWos[91] + 1./2.*armWWos[113] + 5./3.*armWWos[87];
   armWWos[104]=armWWos[65]*armWWos[104];
   armWWos[105]= - 5./6.*armWWos[91] + 1./2.*armWWos[87] + 13./8.*
   armWWos[88] - 7./2.*armWWos[83] + 1489./96. + armWWos[92];
   armWWos[109]= - 1./6.*armWWos[89];
   armWWos[105]=armWWos[125] + 11./24.*armWWos[13] + armWWos[109] + 1./
   3.*armWWos[105] + 1./2.*armWWos[85];
   armWWos[112]=armWWos[91] - 23 - armWWos[87];
   armWWos[113]=1./12.*armWWos[13];
   armWWos[109]= - armWWos[14] + armWWos[113] + armWWos[109] + 1./12.*
   armWWos[112] + armWWos[85];
   armWWos[109]=armWWos[99]*armWWos[109];
   armWWos[109]=armWWos[140] + 1./12.*armWWos[159] - 3./16.*armWWos[72]
    + armWWos[109];
   armWWos[109]=mmZ*armWWos[109];
   armWWos[111]=armWWos[17]*armWWos[111];
   armWWos[111]= - armWWos[20] + armWWos[111];
   armWWos[112]=1./3. + armWWos[119];
   armWWos[112]=armWWos[16]*armWWos[112];
   armWWos[119]=mmZ*armWWos[9];
   armWWos[111]=1./3.*armWWos[119] - 1./3.*armWWos[15] + 1./6.*
   armWWos[111] + armWWos[112];
   armWWos[111]=armWWos[5]*armWWos[111];
   armWWos[112]=19*armWWos[99] + armWWos[159];
   armWWos[112]=1./12.*armWWos[112] + armWWos[151];
   armWWos[112]=armWWos[16]*armWWos[112];
   armWWos[119]=armWWos[18] - 5./24.*armWWos[20] + 7./4.*armWWos[19] - 
   armWWos[95] + 1./6.*armWWos[68];
   armWWos[119]=armWWos[99]*armWWos[119];
   armWWos[122]=31 + armWWos[141];
   armWWos[122]=armWWos[9]*armWWos[122];
   armWWos[125]= - armWWos[99] + 1./3.*armWWos[135];
   armWWos[125]=armWWos[17]*armWWos[125];
   armWWos[126]= - armWWos[15]*armWWos[99];
   armWWos[128]=mmH*armWWos[161];
   armWWos[105]=1./12.*armWWos[128] + 1./4.*armWWos[111] + armWWos[109]
    + 1./2.*armWWos[126] + armWWos[112] + 1./8.*armWWos[125] - 1./4.*
   armWWos[10] + 1./144.*armWWos[122] + 1./2.*armWWos[105] + 
   armWWos[119];
   armWWos[105]=armWWos[65]*armWWos[105];
   armWWos[109]=armWWos[29] + 1 + armWWos[11];
   armWWos[111]=armWWos[204] + armWWos[117] + armWWos[109] + 1./2.*
   armWWos[27];
   armWWos[111]=armWWos[1]*armWWos[111];
   armWWos[109]= - 11./9.*armWWos[9] - 11./6.*armWWos[13] + 
   armWWos[109] + 11./6.*armWWos[27];
   armWWos[109]=armWWos[35]*armWWos[109];
   armWWos[109]=armWWos[111] + 1./3.*armWWos[109];
   armWWos[111]= - 1./9.*armWWos[143] + armWWos[191] + 1 - 1./6.*
   armWWos[83];
   armWWos[111]=armWWos[121]*armWWos[65]*armWWos[111];
   armWWos[105]=1./8.*armWWos[111] + 1./2.*armWWos[109] + armWWos[105];
   armWWos[105]=armWWos[121]*armWWos[105];
   armWWos[109]= - 29./8. - 5./9.*armWWos[11];
   armWWos[111]= - 251./3. - armWWos[6];
   armWWos[111]=armWWos[9]*armWWos[111];
   armWWos[112]=pow(armWWos[6],2);
   armWWos[109]=1./36.*armWWos[111] + 1./6.*armWWos[112] - 17./6.*
   armWWos[13] + armWWos[137] + 17./6.*armWWos[27] - 1./6.*armWWos[123]
    + 1./4.*armWWos[109] + 2./3.*armWWos[29];
   armWWos[109]=armWWos[1]*armWWos[109];
   armWWos[111]=armWWos[8] - armWWos[19];
   armWWos[117]=armWWos[1]*armWWos[111];
   armWWos[111]=armWWos[35]*armWWos[111];
   armWWos[119]=armWWos[16]*armWWos[165];
   armWWos[122]=armWWos[1] + 11./9.*armWWos[35];
   armWWos[122]=mmZ*armWWos[122];
   armWWos[111]=5./2.*armWWos[122] + 5./2.*armWWos[119] + 3*
   armWWos[117] + 11./3.*armWWos[111];
   armWWos[111]=armWWos[5]*armWWos[111];
   armWWos[117]= - 701./24. + armWWos[11];
   armWWos[112]=1./2.*armWWos[112] - 187./6.*armWWos[13] + 11./4.*
   armWWos[30] + 187./6.*armWWos[27] - 1./2.*armWWos[123] + 1./4.*
   armWWos[117] + 2*armWWos[29];
   armWWos[117]= - 947./27. - armWWos[6];
   armWWos[117]=armWWos[9]*armWWos[117];
   armWWos[112]=1./3.*armWWos[112] + 1./4.*armWWos[117];
   armWWos[112]=armWWos[35]*armWWos[112];
   armWWos[117]= - armWWos[1]*armWWos[24];
   armWWos[119]= - armWWos[35]*armWWos[24];
   armWWos[117]=armWWos[117] + 1./3.*armWWos[119];
   armWWos[117]=mmZ*armWWos[117];
   armWWos[104]=armWWos[105] + armWWos[104] + armWWos[111] + 1./6.*
   armWWos[117] + armWWos[109] + 1./3.*armWWos[112];
   armWWos[104]=armWWos[121]*armWWos[104];
   armWWos[105]=49./6. + armWWos[115];
   armWWos[105]=armWWos[17]*armWWos[105];
   armWWos[109]=5./2. - armWWos[10];
   armWWos[109]=armWWos[16]*armWWos[109];
   armWWos[105]=1./6.*armWWos[109] + 1./2.*armWWos[105] - 11./12.*
   armWWos[18] - 1./3.*armWWos[20] + armWWos[114] + armWWos[110] + 
   armWWos[108] + 3*armWWos[163] + 1./12.*armWWos[67];
   armWWos[105]=armWWos[100]*armWWos[105];
   armWWos[108]=21./2.*armWWos[133] + armWWos[116];
   armWWos[108]=mmZ*armWWos[108];
   armWWos[108]=9*armWWos[120] + armWWos[108];
   armWWos[108]=armWWos[5]*armWWos[108];
   armWWos[106]=armWWos[182] - 113./3. + armWWos[106];
   armWWos[106]=armWWos[16]*armWWos[106];
   armWWos[109]=armWWos[66] + armWWos[146];
   armWWos[109]=1./2.*armWWos[109] - armWWos[19];
   armWWos[110]=5 - armWWos[9];
   armWWos[110]=armWWos[17]*armWWos[110];
   armWWos[111]= - 83./6. + armWWos[118];
   armWWos[111]=armWWos[15]*armWWos[111];
   armWWos[112]=179./3. + armWWos[157];
   armWWos[112]=1./3.*armWWos[112] - armWWos[10];
   armWWos[112]=mmZ*armWWos[112];
   armWWos[106]=armWWos[108] + 1./4.*armWWos[112] + 1./2.*armWWos[111]
    + 1./4.*armWWos[106] + 1./8.*armWWos[110] - 9*armWWos[18] + 9*
   armWWos[109] - 19./8.*armWWos[20];
   armWWos[106]=armWWos[5]*armWWos[106];
   armWWos[108]= - 1./3.*armWWos[81] - 13./4. - 1./3.*armWWos[82];
   armWWos[108]=armWWos[134] + 1./2.*armWWos[108] - 3*armWWos[91];
   armWWos[109]= - 1 - 1./6.*armWWos[9];
   armWWos[109]=armWWos[10]*armWWos[109];
   armWWos[108]=1./2.*armWWos[109] + 1./6.*armWWos[9] - 2./3.*
   armWWos[14] + armWWos[113] + 1./2.*armWWos[108] - 2./3.*armWWos[89];
   armWWos[108]=armWWos[100]*armWWos[108];
   armWWos[107]=armWWos[108] + armWWos[145] + armWWos[139] + 
   armWWos[107] - 1./6.*armWWos[71] - 3*armWWos[74] - 3./2.*armWWos[96]
    + 31./16.*armWWos[72] + 34./3.*armWWos[79] - 5*armWWos[76];
   armWWos[107]=mmZ*armWWos[107];
   armWWos[108]=11./2.*armWWos[87] + armWWos[186] + armWWos[185] + 
   armWWos[86] + 445./24.*armWWos[88] + 331./4.*armWWos[83] - 167989./
   288. - 59*armWWos[92];
   armWWos[108]= - 7./2.*armWWos[89] - 19./3.*armWWos[85] + 1./3.*
   armWWos[108] - 35./2.*armWWos[91];
   armWWos[109]=1./6.*armWWos[10] + 65./24. + armWWos[184];
   armWWos[109]=armWWos[100]*armWWos[109];
   armWWos[109]=armWWos[147] + armWWos[109];
   armWWos[109]=armWWos[15]*armWWos[109];
   armWWos[110]= - 1277 + 653*armWWos[9];
   armWWos[110]=armWWos[9]*armWWos[110];
   armWWos[111]=5 - 19./2.*armWWos[9];
   armWWos[111]=armWWos[10]*armWWos[111];
   armWWos[112]=11./12.*armWWos[144] - 2*armWWos[71];
   armWWos[112]=mmH*armWWos[112];
   armWWos[105]=armWWos[112] + armWWos[106] + armWWos[107] + 
   armWWos[109] + armWWos[105] + armWWos[152] + armWWos[164] + 1./6.*
   armWWos[111] + 1./48.*armWWos[110] + armWWos[158] + 37./12.*
   armWWos[14] - 727./144.*armWWos[13] + 1./2.*armWWos[108] - 34./3.*
   armWWos[30];
   armWWos[105]=armWWos[65]*armWWos[105];
   armWWos[106]= - 85./9.*armWWos[27] - armWWos[123] + 13./3.*
   armWWos[29] + 55./18.*armWWos[11] + 1645./108. + armWWos[26];
   armWWos[107]=1./9. + armWWos[166];
   armWWos[107]=armWWos[6]*armWWos[107];
   armWWos[108]=293./27. + armWWos[127];
   armWWos[108]=armWWos[9]*armWWos[108];
   armWWos[106]=1./2.*armWWos[108] + 1./2.*armWWos[107] + 19./9.*
   armWWos[13] + 1./4.*armWWos[106] + 10./3.*armWWos[30];
   armWWos[106]=armWWos[1]*armWWos[106];
   armWWos[107]= - 317./27.*armWWos[27] - 1./3.*armWWos[123] + 23./3.*
   armWWos[29] + 229./18.*armWWos[11] + 7981./324. + 3*armWWos[26];
   armWWos[108]= - 1./3. + armWWos[166];
   armWWos[108]=armWWos[6]*armWWos[108];
   armWWos[109]=781./81. - 3./2.*armWWos[6];
   armWWos[109]=armWWos[9]*armWWos[109];
   armWWos[107]=1./2.*armWWos[109] + 1./6.*armWWos[108] + 59./27.*
   armWWos[13] + 1./4.*armWWos[107] + 110./27.*armWWos[30];
   armWWos[107]=armWWos[35]*armWWos[107];
   armWWos[108]= - 4./3. + armWWos[6];
   armWWos[109]=armWWos[1]*armWWos[108];
   armWWos[108]=armWWos[35]*armWWos[108];
   armWWos[108]=1./3.*armWWos[109] + armWWos[108];
   armWWos[108]=mmZ*armWWos[108];
   armWWos[108]=armWWos[124] + armWWos[108];
   armWWos[108]=armWWos[5]*armWWos[108];
   armWWos[104]=armWWos[104] + armWWos[105] + armWWos[108] + 1./3.*
   armWWos[117] + armWWos[106] + armWWos[107];
   armWWos[104]=armWWos[121]*armWWos[104];
   armWWos[105]= - armWWos[72] - armWWos[77] + 2*armWWos[76];
   armWWos[106]=pow(c,2);
   armWWos[105]=armWWos[106]*armWWos[105];
   armWWos[105]=armWWos[71] - armWWos[75] + 2*armWWos[105] - 23./3.*
   armWWos[72] + 4*armWWos[79] + 23./3.*armWWos[76];
   armWWos[107]=2 + armWWos[9];
   armWWos[107]=armWWos[10]*armWWos[107];
   armWWos[107]=4*armWWos[107] - 8*armWWos[9] - 20*armWWos[14] - 4*
   armWWos[13] + 32*armWWos[89] + 4*armWWos[81] + 63 + 4*armWWos[82];
   armWWos[107]=armWWos[100]*armWWos[107];
   armWWos[105]=8*armWWos[105] + armWWos[107];
   armWWos[105]=mmZ*armWWos[105];
   armWWos[107]= - 2 - armWWos[9];
   armWWos[107]=mmZ*armWWos[107];
   armWWos[108]=armWWos[5]*pow(mmZ,2);
   armWWos[107]=4*armWWos[107] + 5*armWWos[108];
   armWWos[107]=armWWos[5]*armWWos[107];
   armWWos[108]=917./2. + 88*armWWos[92];
   armWWos[109]= - 61./3. - 8*armWWos[106];
   armWWos[109]=armWWos[88]*armWWos[109];
   armWWos[110]= - 16*armWWos[106];
   armWWos[111]= - 33 + armWWos[110];
   armWWos[111]=armWWos[30]*armWWos[111];
   armWWos[112]=67./3. + 8*armWWos[106];
   armWWos[112]=armWWos[13]*armWWos[112];
   armWWos[110]= - 2*armWWos[9] - 127./3. + armWWos[110];
   armWWos[110]=armWWos[9]*armWWos[110];
   armWWos[105]=2*armWWos[107] + armWWos[105] + 2*armWWos[110] + 4*
   armWWos[112] + 4*armWWos[111] + 4*armWWos[109] + 8*armWWos[83] + 1./
   3.*armWWos[108] + 32*armWWos[106];
   armWWos[105]=armWWos[65]*armWWos[105];
   armWWos[106]=5./3.*armWWos[9] - armWWos[13] - armWWos[30] - 17./3.
    + armWWos[27];
   armWWos[107]=armWWos[1]*armWWos[106];
   armWWos[106]=armWWos[35]*armWWos[106];
   armWWos[106]=armWWos[107] + 5./3.*armWWos[106];
   armWWos[107]=armWWos[21] - 2*armWWos[23] + armWWos[22];
   armWWos[107]=armWWos[1]*armWWos[107];
   armWWos[108]= - 1./9.*armWWos[64] - armWWos[23];
   armWWos[108]=2./9.*armWWos[24] + armWWos[21] + 2*armWWos[108] + 
   armWWos[22];
   armWWos[108]=armWWos[35]*armWWos[108];
   armWWos[107]=1./3.*armWWos[107] + armWWos[108];
   armWWos[107]=mmZ*armWWos[107];
   armWWos[106]=4./3.*armWWos[106] + armWWos[107];

      mWWosret = armWWos[101] + armWWos[102] + armWWos[103] + 
      armWWos[104] + armWWos[105] + 2*armWWos[106];
      return mWWosret;
}
