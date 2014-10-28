#include <WW.hpp>
std::complex<long double>
WW<OS>::m20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armWWbar[342], mWWbarret;

    armWWbar[1]=double(nL + nH);
    armWWbar[2]=pow(CW,-1);
    armWWbar[3]=pow(MMH,-1);
    armWWbar[4]=pow(MMZ,-1);
    armWWbar[5]=pow(SW,-1);
    armWWbar[6]=std::real(Tsil::B(0,0,MMZ,mu2));
    armWWbar[7]=std::real(Tsil::B(0,0,MMW,mu2));
    armWWbar[8]=Tsil::I2(0,0,MMZ,mu2);
    armWWbar[9]=Tsil::I2(0,0,MMW,mu2);
    armWWbar[10]=Tsil::B(MMZ,MMH,MMZ,mu2);
    armWWbar[11]=Tsil::B(MMW,MMH,MMW,mu2);
    armWWbar[12]=Tsil::B(MMW,MMZ,MMW,mu2);
    armWWbar[13]=Tsil::B(MMW,MMW,MMZ,mu2);
    armWWbar[14]=Tsil::B(0,MMH,MMW,mu2);
    armWWbar[15]=Tsil::B(0,MMZ,MMW,mu2);
    armWWbar[16]=Tsil::Beps(MMW,MMH,MMW,mu2);
    armWWbar[17]=Tsil::Beps(MMW,MMZ,MMW,mu2);
    armWWbar[18]=Tsil::A(MMH,mu2);
    armWWbar[19]=Tsil::A(MMZ,mu2);
    armWWbar[20]=Tsil::A(MMW,mu2);
    armWWbar[21]=Tsil::Aeps(MMH,mu2);
    armWWbar[22]=Tsil::Aeps(MMZ,mu2);
    armWWbar[23]=Tsil::Aeps(MMW,mu2);
    armWWbar[24]=protW0Z00->M(0);
    armWWbar[25]=prot0W0Z0->M(0);
    armWWbar[26]=prot00W00->M(0);
    armWWbar[27]=prot0000Z->M(0);
    armWWbar[28]=protHW00->Uxzuv(0);
    armWWbar[29]=protW0Z00->Uzxyv(0);
    armWWbar[30]=protWtZ00->Uxzuv(0);
    armWWbar[31]=protHW00->Txuv(0);
    armWWbar[32]=prot0W0Z0->Tuxv(0);
    armWWbar[33]=protWtZ00->Txuv(0);
    armWWbar[34]=double(nH);
    armWWbar[35]=pow(MMt,-1);
    armWWbar[36]=Tsil::B(MMt,MMt,MMZ,mu2);
    armWWbar[37]=Tsil::B(0,MMt,MMW,mu2);
    armWWbar[38]=Tsil::A(MMt,mu2);
    armWWbar[39]=Tsil::B(MMt,MMt,MMH,mu2);
    armWWbar[40]=double(nL);
    armWWbar[41]=Tsil::I2(MMH,MMt,MMt,mu2);
    armWWbar[42]=Tsil::I2(MMZ,MMt,MMt,mu2);
    armWWbar[43]=Tsil::I2(0,MMW,MMt,mu2);
    armWWbar[44]=Tsil::B(MMH,MMH,MMH,mu2);
    armWWbar[45]=Tsil::B(MMH,MMt,MMt,mu2);
    armWWbar[46]=Tsil::B(MMZ,MMZ,MMH,mu2);
    armWWbar[47]=Tsil::B(MMZ,MMt,MMt,mu2);
    armWWbar[48]=Tsil::B(MMW,MMW,MMH,mu2);
    armWWbar[49]=std::real(Tsil::B(0,MMW,MMt,mu2));
    armWWbar[50]=Tsil::Beps(0,MMt,MMW,mu2);
    armWWbar[51]=Tsil::Aeps(MMt,mu2);
    armWWbar[52]=prot00tt0->M(0);
    armWWbar[53]=prot00tt0->Tuxv(0);
    armWWbar[54]=protWtZ00->M(0);
    armWWbar[55]=protW0Htt->M(0);
    armWWbar[56]=protW0Ztt->M(0);
    armWWbar[57]=prot0Wt0t->M(0);
    armWWbar[58]=prot00Wt0->M(0);
    armWWbar[59]=prot00ttZ->M(0);
    armWWbar[60]=protW0Htt->Uzxyv(0);
    armWWbar[61]=protWtZ00->Uzxyv(0);
    armWWbar[62]=protW0Htt->Uxzuv(0);
    armWWbar[63]=protW0Ztt->Uxzuv(0);
    armWWbar[64]=prot0Wt0t->Uyuzv(0);
    armWWbar[65]=protW0Htt->Uyuzv(0);
    armWWbar[66]=protW0Ztt->Uyuzv(0);
    armWWbar[67]=protWtZ00->Uuyxv(0);
    armWWbar[68]=protW0Htt->Tzyv(0);
    armWWbar[69]=protWtZ00->Tzyv(0);
    armWWbar[70]=protW0Htt->Tvyz(0);
    armWWbar[71]=protWtZ00->Tyzv(0);
    armWWbar[72]=protW0Htt->Suxv(0);
    armWWbar[73]=protW0Htt->Svyz(0);
    armWWbar[74]=protWtZ00->Svyz(0);
    armWWbar[75]=prot00000->M(0);
    armWWbar[76]=double(boson);
    armWWbar[77]=Tsil::I2(MMH,MMH,MMH,mu2);
    armWWbar[78]=Tsil::I2(MMH,MMZ,MMZ,mu2);
    armWWbar[79]=Tsil::I2(MMH,MMW,MMW,mu2);
    armWWbar[80]=Tsil::I2(MMW,MMW,MMZ,mu2);
    armWWbar[81]=protWHHWW->M(0);
    armWWbar[82]=protWHZWW->M(0);
    armWWbar[83]=protWZZWW->M(0);
    armWWbar[84]=protWWHHH->M(0);
    armWWbar[85]=protWWHZZ->M(0);
    armWWbar[86]=protW0HWW->M(0);
    armWWbar[87]=protW0ZWW->M(0);
    armWWbar[88]=prot0WW0W->M(0);
    armWWbar[89]=protWH0H->Vxzuv(0);
    armWWbar[90]=protWZ0Z->Vxzuv(0);
    armWWbar[91]=protWHHWW->Uzxyv(0);
    armWWbar[92]=protWHZWW->Uyuzv(0);
    armWWbar[93]=protWHZWW->Uzxyv(0);
    armWWbar[94]=protWZZWW->Uzxyv(0);
    armWWbar[95]=protWWHHH->Uyuzv(0);
    armWWbar[96]=protWWHZZ->Uxzuv(0);
    armWWbar[97]=protWHHWW->Uxzuv(0);
    armWWbar[98]=protWWHZZ->Uyuzv(0);
    armWWbar[99]=protWHZWW->Uxzuv(0);
    armWWbar[100]=protWHZWW->Tyzv(0);
    armWWbar[101]=protWHHWW->Tyzv(0);
    armWWbar[102]=protW0HWW->Tzyv(0);
    armWWbar[103]=protWHZWW->Tzyv(0);
    armWWbar[104]=protW0ZWW->Tzyv(0);
    armWWbar[105]=protWHHWW->Svyz(0);
    armWWbar[106]=protWHZWW->Svyz(0);
    armWWbar[107]=protWZZWW->Svyz(0);
    armWWbar[108]=protWtZ00->Suxv(0);
    armWWbar[109]=protWWZZH->M(0);
    armWWbar[110]=1/(MMt - MMW);
    armWWbar[111]=1/(4*MMt - MMZ);
    armWWbar[112]=1/( - 4*MMW + MMH);
    armWWbar[113]=1/( - 4*MMZ + MMH);
   armWWbar[114]= - 1./6. - armWWbar[7];
   armWWbar[115]=1./2.*armWWbar[6];
   armWWbar[116]=1./3.*armWWbar[114] + armWWbar[115];
   armWWbar[116]=armWWbar[1]*armWWbar[116];
   armWWbar[117]= - 25*armWWbar[47];
   armWWbar[118]=137./12. + armWWbar[117];
   armWWbar[119]=19./24.*armWWbar[13];
   armWWbar[120]=11./18.*armWWbar[6];
   armWWbar[121]=armWWbar[120] + 7./54. - armWWbar[7];
   armWWbar[121]=armWWbar[40]*armWWbar[121];
   armWWbar[122]=armWWbar[18]*armWWbar[112];
   armWWbar[123]=1./2.*armWWbar[122];
   armWWbar[124]= - 3./2.*armWWbar[49];
   armWWbar[125]= - 3*armWWbar[45];
   armWWbar[118]=armWWbar[123] + armWWbar[121] + armWWbar[116] + 233./
   24.*armWWbar[12] + 3./2.*armWWbar[10] + armWWbar[119] + 
   armWWbar[124] + 1./8.*armWWbar[118] + armWWbar[125];
   armWWbar[118]=armWWbar[37]*armWWbar[118];
   armWWbar[126]=3*armWWbar[41];
   armWWbar[127]= - 8*armWWbar[51] + armWWbar[126] + 2*armWWbar[43];
   armWWbar[128]= - armWWbar[18] - 9./4.*armWWbar[19];
   armWWbar[128]=armWWbar[37]*armWWbar[128];
   armWWbar[129]= - 5./2.*armWWbar[37] - 3 - 2*armWWbar[39];
   armWWbar[129]=armWWbar[20]*armWWbar[129];
   armWWbar[130]= - 2*armWWbar[49] - 3 - 8*armWWbar[45];
   armWWbar[130]=armWWbar[38]*armWWbar[130];
   armWWbar[131]= - 3*armWWbar[39];
   armWWbar[132]= - 5 + armWWbar[131];
   armWWbar[132]=armWWbar[18]*armWWbar[132];
   armWWbar[127]=3*armWWbar[130] + 3*armWWbar[129] + armWWbar[128] + 
   armWWbar[132] - 6*armWWbar[23] + 3*armWWbar[127] - 10*armWWbar[21];
   armWWbar[127]=armWWbar[3]*armWWbar[127];
   armWWbar[128]= - armWWbar[58] - 1./3.*armWWbar[52];
   armWWbar[129]=2*armWWbar[57];
   armWWbar[130]=armWWbar[128] + armWWbar[129];
   armWWbar[133]=armWWbar[59] - 1./2.*armWWbar[54] + 2*armWWbar[130] - 
   1./2.*armWWbar[56];
   armWWbar[134]= - 4*armWWbar[45];
   armWWbar[135]= - armWWbar[49] + 3 + armWWbar[134];
   armWWbar[135]=armWWbar[3]*armWWbar[135];
   armWWbar[136]=5 - armWWbar[70];
   armWWbar[136]=armWWbar[112]*armWWbar[136];
   armWWbar[133]=3*armWWbar[135] + armWWbar[136] + 1./3.*armWWbar[133]
    - 2*armWWbar[55];
   armWWbar[133]=MMt*armWWbar[133];
   armWWbar[135]=4*armWWbar[64];
   armWWbar[137]=31./3.*armWWbar[53];
   armWWbar[138]=armWWbar[137] - 935./12. + armWWbar[135];
   armWWbar[139]= - 5./2.*armWWbar[21];
   armWWbar[126]= - armWWbar[23] + armWWbar[139] - 6*armWWbar[51] + 
   armWWbar[43] + armWWbar[126] - 2*armWWbar[73];
   armWWbar[126]=armWWbar[112]*armWWbar[126];
   armWWbar[140]=1 + armWWbar[39];
   armWWbar[140]=armWWbar[112]*armWWbar[140];
   armWWbar[141]= - armWWbar[37]*armWWbar[112];
   armWWbar[142]=6*armWWbar[140] + armWWbar[141];
   armWWbar[142]=armWWbar[20]*armWWbar[142];
   armWWbar[143]=41./8.*armWWbar[61];
   armWWbar[144]=7./3.*armWWbar[67];
   armWWbar[145]=17./32.*armWWbar[66];
   armWWbar[146]= - 1./6.*armWWbar[63];
   armWWbar[147]= - 6*armWWbar[62];
   armWWbar[148]=17./12.*armWWbar[33];
   armWWbar[149]=1./6.*armWWbar[30];
   armWWbar[150]= - 35./6.*armWWbar[69];
   armWWbar[151]=311./32.*armWWbar[71];
   armWWbar[152]= - 491./96.*armWWbar[50];
   armWWbar[153]= - 73./24.*armWWbar[70];
   armWWbar[154]= - 5./2.*armWWbar[68];
   armWWbar[155]=9*armWWbar[45];
   armWWbar[156]=3./2.*armWWbar[49];
   armWWbar[157]= - 41./8.*armWWbar[17];
   armWWbar[158]=6*armWWbar[16];
   armWWbar[159]= - 6*armWWbar[39];
   armWWbar[160]=3*armWWbar[39];
   armWWbar[161]=1 + armWWbar[160];
   armWWbar[161]=armWWbar[11]*armWWbar[161];
   armWWbar[162]=2*armWWbar[161];
   armWWbar[163]= - armWWbar[38]*armWWbar[112];
   armWWbar[164]=4*armWWbar[163];
   armWWbar[165]=MMH*armWWbar[55];
   armWWbar[166]=3*armWWbar[165];
   armWWbar[167]=9./4.*armWWbar[65];
   armWWbar[168]=5./2. + armWWbar[131];
   armWWbar[168]=armWWbar[18]*armWWbar[112]*armWWbar[168];
   armWWbar[138]=armWWbar[133] + armWWbar[166] + armWWbar[127] + 
   armWWbar[164] + armWWbar[142] + armWWbar[118] + armWWbar[162] + 
   armWWbar[168] + armWWbar[126] + armWWbar[159] + armWWbar[158] + 
   armWWbar[157] + armWWbar[156] + armWWbar[155] + armWWbar[154] + 
   armWWbar[153] + armWWbar[152] + armWWbar[151] + armWWbar[150] + 
   armWWbar[149] + armWWbar[148] + armWWbar[147] + armWWbar[167] + 
   armWWbar[146] + armWWbar[145] + armWWbar[144] + 1./3.*armWWbar[138]
    + armWWbar[143];
   armWWbar[138]=MMt*armWWbar[138];
   armWWbar[169]=2879./36. + armWWbar[117];
   armWWbar[170]=6*armWWbar[45];
   armWWbar[171]=27./2.*armWWbar[44];
   armWWbar[172]=3*armWWbar[48];
   armWWbar[173]=3./2.*armWWbar[46];
   armWWbar[174]=89./6.*armWWbar[12];
   armWWbar[175]=2*armWWbar[11];
   armWWbar[176]= - 89./96.*armWWbar[37];
   armWWbar[177]=armWWbar[20]*armWWbar[112];
   armWWbar[178]=6*armWWbar[177];
   armWWbar[179]=armWWbar[38]*armWWbar[112];
   armWWbar[180]=3*armWWbar[179];
   armWWbar[181]= - 3./4.*armWWbar[49];
   armWWbar[182]= - 1./2.*armWWbar[10];
   armWWbar[169]=armWWbar[180] + armWWbar[178] + armWWbar[176] + 
   armWWbar[175] + armWWbar[121] + armWWbar[116] + armWWbar[174] + 
   armWWbar[182] + armWWbar[119] + armWWbar[173] + armWWbar[172] + 
   armWWbar[171] + armWWbar[181] + 1./8.*armWWbar[169] + armWWbar[170];
   armWWbar[169]=armWWbar[38]*armWWbar[169];
   armWWbar[183]=3*armWWbar[65];
   armWWbar[184]= - 3805./108. + armWWbar[183];
   armWWbar[185]=1./2.*armWWbar[60];
   armWWbar[186]= - 3./4.*armWWbar[50];
   armWWbar[187]=5./2.*armWWbar[70];
   armWWbar[188]= - 21./4.*armWWbar[68];
   armWWbar[189]= - 7./2.*armWWbar[16];
   armWWbar[190]=7./54.*armWWbar[36];
   armWWbar[191]=3*armWWbar[62];
   armWWbar[184]=armWWbar[190] + armWWbar[160] + armWWbar[189] + 
   armWWbar[125] + armWWbar[188] + armWWbar[187] + armWWbar[186] + 
   armWWbar[185] + 1./4.*armWWbar[184] + armWWbar[191];
   armWWbar[192]= - 3./2.*armWWbar[39];
   armWWbar[193]=armWWbar[190] - 55./27. + armWWbar[192];
   armWWbar[193]=armWWbar[11]*armWWbar[193];
   armWWbar[194]=3*armWWbar[45];
   armWWbar[195]=1./3. + armWWbar[194];
   armWWbar[196]= - 1./3.*armWWbar[10];
   armWWbar[195]=1./2.*armWWbar[195] + armWWbar[196];
   armWWbar[197]=1./6.*armWWbar[11];
   armWWbar[198]=armWWbar[195] + armWWbar[197];
   armWWbar[198]=armWWbar[37]*armWWbar[198];
   armWWbar[199]=1./3.*armWWbar[10];
   armWWbar[184]=1./2.*armWWbar[198] + armWWbar[193] + 1./2.*
   armWWbar[184] + armWWbar[199];
   armWWbar[184]=MMH*armWWbar[184];
   armWWbar[193]= - 43*armWWbar[74];
   armWWbar[198]=armWWbar[193] + 79./12.*armWWbar[42];
   armWWbar[200]= - 13*armWWbar[72];
   armWWbar[198]=1./8.*armWWbar[198] + armWWbar[200];
   armWWbar[201]=965 + 1361*armWWbar[36];
   armWWbar[202]= - 7*armWWbar[37];
   armWWbar[201]=1./27.*armWWbar[201] + armWWbar[202];
   armWWbar[201]=armWWbar[20]*armWWbar[201];
   armWWbar[203]= - 7./18.*armWWbar[36];
   armWWbar[204]=armWWbar[203] + 475./36. + armWWbar[160];
   armWWbar[204]=1./2.*armWWbar[18]*armWWbar[204];
   armWWbar[205]=armWWbar[18] - 5./4.*armWWbar[19];
   armWWbar[206]= - 9*armWWbar[38];
   armWWbar[205]=armWWbar[206] + 3*armWWbar[205] - 23*armWWbar[20];
   armWWbar[205]=armWWbar[3]*armWWbar[38]*armWWbar[205];
   armWWbar[207]= - 5129./16. - 331*armWWbar[36];
   armWWbar[207]=armWWbar[19]*armWWbar[207];
   armWWbar[208]= - 71./12.*armWWbar[73];
   armWWbar[209]= - 5./12.*armWWbar[8];
   armWWbar[210]=23./8.*armWWbar[43];
   armWWbar[211]=115./24.*armWWbar[21];
   armWWbar[212]=5*armWWbar[18];
   armWWbar[213]=armWWbar[212] + 541./24.*armWWbar[19];
   armWWbar[213]=1./8.*armWWbar[37]*armWWbar[213];
   armWWbar[214]= - 3./2.*armWWbar[41];
   armWWbar[138]=armWWbar[138] + armWWbar[184] + armWWbar[205] + 
   armWWbar[169] + 1./4.*armWWbar[201] + armWWbar[213] + 1./36.*
   armWWbar[207] + armWWbar[204] + 1795./72.*armWWbar[23] + 
   armWWbar[211] + 1621./576.*armWWbar[22] + 2503./96.*armWWbar[51] + 
   armWWbar[210] + armWWbar[209] + armWWbar[208] + 1./3.*armWWbar[198]
    + armWWbar[214];
   armWWbar[138]=MMt*armWWbar[138];
   armWWbar[169]= - 11*armWWbar[47];
   armWWbar[198]= - 257./18. + armWWbar[169];
   armWWbar[198]=armWWbar[10] + 1./12.*armWWbar[13] + armWWbar[124] + 1.
   /2.*armWWbar[198] + armWWbar[125];
   armWWbar[201]=7./3.*armWWbar[12];
   armWWbar[207]= - 1./3. + armWWbar[6];
   armWWbar[215]=armWWbar[1]*armWWbar[207];
   armWWbar[216]=1./2.*armWWbar[215];
   armWWbar[207]=armWWbar[40]*armWWbar[207];
   armWWbar[217]=11./18.*armWWbar[207];
   armWWbar[218]=1./4.*armWWbar[122];
   armWWbar[219]=1./2.*armWWbar[11];
   armWWbar[198]=armWWbar[219] + armWWbar[218] + armWWbar[217] + 
   armWWbar[216] + 1./2.*armWWbar[198] + armWWbar[201];
   armWWbar[198]=armWWbar[37]*armWWbar[198];
   armWWbar[220]=3./2.*armWWbar[41];
   armWWbar[221]= - 4*armWWbar[51] + armWWbar[220] + armWWbar[43];
   armWWbar[222]= - 3*armWWbar[23];
   armWWbar[132]=1./2.*armWWbar[132] + armWWbar[222] + 3*armWWbar[221]
    - 5*armWWbar[21];
   armWWbar[221]= - armWWbar[18] - 3./2.*armWWbar[19];
   armWWbar[221]=armWWbar[37]*armWWbar[221];
   armWWbar[223]= - 3./2. - armWWbar[39];
   armWWbar[224]=armWWbar[223] - armWWbar[37];
   armWWbar[224]=armWWbar[20]*armWWbar[224];
   armWWbar[134]= - armWWbar[49] - 3./2. + armWWbar[134];
   armWWbar[134]=3*armWWbar[38]*armWWbar[134];
   armWWbar[221]=armWWbar[134] + 3*armWWbar[224] + armWWbar[132] + 1./2.
   *armWWbar[221];
   armWWbar[221]=armWWbar[3]*armWWbar[221];
   armWWbar[224]= - 3*armWWbar[51];
   armWWbar[225]=1./2.*armWWbar[43];
   armWWbar[226]= - 1./2.*armWWbar[23];
   armWWbar[220]=armWWbar[226] - 5./4.*armWWbar[21] + armWWbar[224] + 
   armWWbar[225] + armWWbar[220] - armWWbar[73];
   armWWbar[220]=armWWbar[112]*armWWbar[220];
   armWWbar[130]=armWWbar[130] - armWWbar[56];
   armWWbar[130]=armWWbar[59] + 2*armWWbar[130] + armWWbar[54];
   armWWbar[136]=1./2.*armWWbar[136];
   armWWbar[227]= - 1./2.*armWWbar[49] + 3./2. - 2*armWWbar[45];
   armWWbar[227]=3*armWWbar[3]*armWWbar[227];
   armWWbar[130]=armWWbar[227] + armWWbar[136] + 1./3.*armWWbar[130] - 
   armWWbar[55];
   armWWbar[130]=MMt*armWWbar[130];
   armWWbar[228]=armWWbar[137] - 6749./288. + armWWbar[135];
   armWWbar[141]=1./2.*armWWbar[141];
   armWWbar[140]=3*armWWbar[140] + armWWbar[141];
   armWWbar[140]=armWWbar[20]*armWWbar[140];
   armWWbar[229]= - 1./3.*armWWbar[30];
   armWWbar[230]= - 3*armWWbar[62];
   armWWbar[231]= - 73./48.*armWWbar[70];
   armWWbar[232]= - 5./4.*armWWbar[68];
   armWWbar[233]=3./4.*armWWbar[49];
   armWWbar[234]=3*armWWbar[16];
   armWWbar[235]=1./2.*armWWbar[168];
   armWWbar[236]=2*armWWbar[163];
   armWWbar[165]=3./2.*armWWbar[165];
   armWWbar[237]=9./2.*armWWbar[45];
   armWWbar[130]=armWWbar[130] + armWWbar[165] + armWWbar[221] + 
   armWWbar[236] + armWWbar[140] + armWWbar[198] + armWWbar[161] + 
   armWWbar[235] + armWWbar[220] + armWWbar[131] + armWWbar[234] - 7./8.
   *armWWbar[17] + armWWbar[233] + armWWbar[237] + armWWbar[232] + 
   armWWbar[231] - 139./48.*armWWbar[50] + 29./8.*armWWbar[71] - 67./12.
   *armWWbar[69] + armWWbar[229] + 25./24.*armWWbar[33] + armWWbar[230]
    + 9./8.*armWWbar[65] - 2./3.*armWWbar[63] + 11./16.*armWWbar[66] + 
   13./12.*armWWbar[67] + 1./3.*armWWbar[228] + 15./8.*armWWbar[61];
   armWWbar[130]=MMt*armWWbar[130];
   armWWbar[198]= - 2279./108. + armWWbar[169];
   armWWbar[221]= - 3./8.*armWWbar[49];
   armWWbar[228]=27./4.*armWWbar[44];
   armWWbar[238]=3./2.*armWWbar[48];
   armWWbar[239]=3./4.*armWWbar[46];
   armWWbar[240]=3*armWWbar[177];
   armWWbar[179]=3./2.*armWWbar[179];
   armWWbar[241]=1./2.*armWWbar[10];
   armWWbar[198]=armWWbar[179] + armWWbar[240] - 25./48.*armWWbar[37]
    + armWWbar[219] + armWWbar[217] + armWWbar[216] + 77./24.*
   armWWbar[12] + armWWbar[241] + 1./24.*armWWbar[13] + armWWbar[239]
    + armWWbar[238] + armWWbar[228] + armWWbar[221] + 1./4.*
   armWWbar[198] + armWWbar[194];
   armWWbar[198]=armWWbar[38]*armWWbar[198];
   armWWbar[242]= - 437./12. + armWWbar[183];
   armWWbar[242]=armWWbar[160] + armWWbar[189] + armWWbar[125] + 
   armWWbar[188] + armWWbar[187] + armWWbar[186] + armWWbar[185] + 1./4.
   *armWWbar[242] + armWWbar[191];
   armWWbar[243]= - 3./4.*armWWbar[39];
   armWWbar[244]= - 1 + armWWbar[243];
   armWWbar[244]=armWWbar[11]*armWWbar[244];
   armWWbar[245]=armWWbar[194] - armWWbar[11];
   armWWbar[245]=armWWbar[37]*armWWbar[245];
   armWWbar[242]=1./8.*armWWbar[245] + 1./4.*armWWbar[242] + 
   armWWbar[244];
   armWWbar[242]=MMH*armWWbar[242];
   armWWbar[244]=3*armWWbar[18];
   armWWbar[245]=1./2.*armWWbar[19];
   armWWbar[206]=armWWbar[206] - 5*armWWbar[20] + armWWbar[244] + 
   armWWbar[245];
   armWWbar[206]=armWWbar[3]*armWWbar[38]*armWWbar[206];
   armWWbar[246]=35*armWWbar[74] + 1303./18.*armWWbar[42];
   armWWbar[246]=1./8.*armWWbar[246] - 7*armWWbar[72];
   armWWbar[247]=169./12. + armWWbar[160];
   armWWbar[247]=armWWbar[18]*armWWbar[247];
   armWWbar[248]= - 505./36. - armWWbar[36];
   armWWbar[248]=armWWbar[19]*armWWbar[248];
   armWWbar[249]=armWWbar[18] + 23./18.*armWWbar[19];
   armWWbar[249]=armWWbar[37]*armWWbar[249];
   armWWbar[250]= - 29 - 25./8.*armWWbar[36];
   armWWbar[250]=1./3.*armWWbar[250] - 17./8.*armWWbar[37];
   armWWbar[250]=armWWbar[20]*armWWbar[250];
   armWWbar[251]= - 3./4.*armWWbar[41];
   armWWbar[252]= - 71./24.*armWWbar[73];
   armWWbar[253]=115./48.*armWWbar[21];
   armWWbar[130]=armWWbar[130] + armWWbar[242] + 1./2.*armWWbar[206] + 
   armWWbar[198] + 1./3.*armWWbar[250] + 7./16.*armWWbar[249] + 1./24.*
   armWWbar[248] + 1./4.*armWWbar[247] + 1331./54.*armWWbar[23] + 
   armWWbar[253] - 1907./864.*armWWbar[22] + 2779./432.*armWWbar[51] + 
   3./4.*armWWbar[43] - 1./6.*armWWbar[8] + armWWbar[252] + 1./3.*
   armWWbar[246] + armWWbar[251];
   armWWbar[130]=MMt*armWWbar[130];
   armWWbar[198]= - 83./6. - 15*armWWbar[61];
   armWWbar[206]=7*armWWbar[47];
   armWWbar[242]= - 7*armWWbar[12] + 5./12. + armWWbar[206];
   armWWbar[242]=armWWbar[37]*armWWbar[242];
   armWWbar[198]=1./12.*armWWbar[242] + 15./4.*armWWbar[17] + 125./144.
   *armWWbar[50] - 109./48.*armWWbar[71] - 67./9.*armWWbar[69] + 
   armWWbar[229] + 1./3.*armWWbar[63] - 77./144.*armWWbar[66] + 1./4.*
   armWWbar[198] - 1./3.*armWWbar[67];
   armWWbar[242]=armWWbar[56] + armWWbar[54];
   armWWbar[246]=1./3.*armWWbar[59];
   armWWbar[242]=1./2.*armWWbar[242] + armWWbar[246];
   armWWbar[242]=MMt*armWWbar[242];
   armWWbar[198]=1./2.*armWWbar[198] + 1./3.*armWWbar[242];
   armWWbar[198]=MMt*armWWbar[198];
   armWWbar[242]=35./12. + armWWbar[47];
   armWWbar[247]= - 13*armWWbar[12];
   armWWbar[242]= - 49./48.*armWWbar[37] + 7./4.*armWWbar[242] + 
   armWWbar[247];
   armWWbar[242]=armWWbar[38]*armWWbar[242];
   armWWbar[248]=43./9.*armWWbar[74] - 3./4.*armWWbar[42];
   armWWbar[249]= - 53./8. + armWWbar[36];
   armWWbar[249]=armWWbar[19]*armWWbar[249];
   armWWbar[250]=armWWbar[37]*armWWbar[19];
   armWWbar[254]= - 2 - 7./8.*armWWbar[36];
   armWWbar[254]=armWWbar[20]*armWWbar[254];
   armWWbar[198]=armWWbar[198] + 1./6.*armWWbar[242] + 1./9.*
   armWWbar[254] + 205./576.*armWWbar[250] + 7./72.*armWWbar[249] + 1./
   8.*armWWbar[23] + 13./192.*armWWbar[22] - 289./288.*armWWbar[51] + 1.
   /8.*armWWbar[43] - 1./36.*armWWbar[8] + 1./8.*armWWbar[248] - 1./3.*
   armWWbar[72];
   armWWbar[198]=MMt*armWWbar[198];
   armWWbar[242]= - 1 - armWWbar[61];
   armWWbar[248]=1./2.*armWWbar[17];
   armWWbar[242]=armWWbar[248] - 1./2.*armWWbar[71] + 1./2.*
   armWWbar[242] - armWWbar[69];
   armWWbar[242]=MMt*armWWbar[242];
   armWWbar[249]=1 - armWWbar[12];
   armWWbar[254]=armWWbar[38]*armWWbar[249];
   armWWbar[255]= - armWWbar[51] + armWWbar[254];
   armWWbar[242]=1./2.*armWWbar[255] + armWWbar[242];
   armWWbar[255]=pow(armWWbar[2],2);
   armWWbar[242]=armWWbar[255]*MMt*armWWbar[242];
   armWWbar[256]=37./32.*armWWbar[19] - 2./3.*armWWbar[20];
   armWWbar[256]=1./3.*armWWbar[256] - 7./32.*armWWbar[38];
   armWWbar[256]=armWWbar[38]*armWWbar[256];
   armWWbar[198]=1./4.*armWWbar[242] + armWWbar[256] + armWWbar[198];
   armWWbar[198]=armWWbar[255]*armWWbar[198];
   armWWbar[242]= - 17./2.*armWWbar[36];
   armWWbar[256]= - 5./2.*armWWbar[6];
   armWWbar[257]=armWWbar[256] - 313./3. + armWWbar[242];
   armWWbar[257]=armWWbar[18]*armWWbar[257];
   armWWbar[258]=11./6.*armWWbar[51] + armWWbar[214] + 11./3.*
   armWWbar[73];
   armWWbar[257]=armWWbar[258] + 1./27.*armWWbar[257];
   armWWbar[259]= - 5./4.*armWWbar[6] + 19./3. - 17./4.*armWWbar[36];
   armWWbar[259]=armWWbar[11]*armWWbar[259];
   armWWbar[260]=1 + armWWbar[68];
   armWWbar[259]=13./4.*armWWbar[260] + 1./9.*armWWbar[259];
   armWWbar[259]=MMH*armWWbar[259];
   armWWbar[261]=5./4.*armWWbar[6];
   armWWbar[262]=armWWbar[261] - 19./3. + 17./4.*armWWbar[36];
   armWWbar[262]=armWWbar[20]*armWWbar[262];
   armWWbar[257]=1./3.*armWWbar[259] - 7./3.*armWWbar[38] + 1./2.*
   armWWbar[257] + 1./27.*armWWbar[262];
   armWWbar[257]=MMH*armWWbar[257];
   armWWbar[259]= - 11*armWWbar[18] + 1231./18.*armWWbar[19];
   armWWbar[259]= - 179./36.*armWWbar[38] + 1./2.*armWWbar[259] - 2749./
   9.*armWWbar[20];
   armWWbar[259]=armWWbar[38]*armWWbar[259];
   armWWbar[257]=1./6.*armWWbar[259] + armWWbar[257];
   armWWbar[130]=armWWbar[198] + 1./2.*armWWbar[257] + armWWbar[130];
   armWWbar[130]=armWWbar[255]*armWWbar[130];
   armWWbar[198]= - 11./2.*armWWbar[6];
   armWWbar[257]= - 11./2.*armWWbar[36];
   armWWbar[259]=armWWbar[198] - 655./3. + armWWbar[257];
   armWWbar[259]=armWWbar[18]*armWWbar[259];
   armWWbar[262]=11./2.*armWWbar[6];
   armWWbar[263]=armWWbar[262] - 47./3. + 11./2.*armWWbar[36];
   armWWbar[263]=1./18.*armWWbar[263] - armWWbar[37];
   armWWbar[263]=armWWbar[20]*armWWbar[263];
   armWWbar[264]=1./3.*armWWbar[37]*armWWbar[18];
   armWWbar[259]=1./3.*armWWbar[263] + armWWbar[264] + armWWbar[258] + 
   1./54.*armWWbar[259];
   armWWbar[263]=armWWbar[198] + 47./3. + armWWbar[257];
   armWWbar[263]=armWWbar[11]*armWWbar[263];
   armWWbar[263]=13*armWWbar[260] + 1./9.*armWWbar[263];
   armWWbar[265]=armWWbar[37]*armWWbar[11];
   armWWbar[263]=1./2.*armWWbar[263] + armWWbar[265];
   armWWbar[263]=MMH*armWWbar[263];
   armWWbar[266]= - 1./9.*armWWbar[11] - 59./9. + armWWbar[10];
   armWWbar[266]=armWWbar[38]*armWWbar[266];
   armWWbar[259]=1./6.*armWWbar[263] + 1./2.*armWWbar[259] + 1./3.*
   armWWbar[266];
   armWWbar[259]=MMH*armWWbar[259];
   armWWbar[263]= - 49*armWWbar[18];
   armWWbar[266]= - 1021./3.*armWWbar[20] + armWWbar[263] - 2323./8.*
   armWWbar[19];
   armWWbar[266]=1./9.*armWWbar[266] - 173./8.*armWWbar[38];
   armWWbar[266]=armWWbar[38]*armWWbar[266];
   armWWbar[130]=armWWbar[130] + armWWbar[138] + 1./4.*armWWbar[266] + 
   armWWbar[259];
   armWWbar[130]=armWWbar[255]*armWWbar[130];
   armWWbar[138]=9*armWWbar[39];
   armWWbar[266]=8./3. + armWWbar[138];
   armWWbar[267]=7./3.*armWWbar[36];
   armWWbar[268]=3*armWWbar[37];
   armWWbar[266]=armWWbar[268] + 2*armWWbar[266] + armWWbar[267];
   armWWbar[266]=armWWbar[3]*armWWbar[38]*armWWbar[266];
   armWWbar[269]= - 3./4.*armWWbar[37];
   armWWbar[270]=1./2.*armWWbar[36];
   armWWbar[271]=armWWbar[269] + 11./36.*armWWbar[6] + 2./27. + 
   armWWbar[270];
   armWWbar[271]=armWWbar[37]*armWWbar[271];
   armWWbar[272]= - 7./2.*armWWbar[36];
   armWWbar[273]= - 8 + armWWbar[272];
   armWWbar[274]=pow(armWWbar[3],2);
   armWWbar[275]= - MMt*armWWbar[274]*armWWbar[38]*armWWbar[39];
   armWWbar[266]=72*armWWbar[275] + armWWbar[266] + 1./9.*armWWbar[273]
    + armWWbar[271];
   armWWbar[266]=MMt*armWWbar[266];
   armWWbar[271]= - 1./2.*armWWbar[36];
   armWWbar[273]= - 19./3.*armWWbar[37] + 11./6.*armWWbar[6] - 89./9.
    + armWWbar[271];
   armWWbar[273]=armWWbar[38]*armWWbar[273];
   armWWbar[276]=pow(armWWbar[38],2);
   armWWbar[277]=armWWbar[3]*armWWbar[276];
   armWWbar[273]=1./2.*armWWbar[273] + 16*armWWbar[277];
   armWWbar[266]=1./3.*armWWbar[273] + armWWbar[266];
   armWWbar[266]=MMt*armWWbar[266];
   armWWbar[266]= - 8./9.*armWWbar[276] + armWWbar[266];
   armWWbar[273]=17./2.*armWWbar[36];
   armWWbar[278]=5./2.*armWWbar[6];
   armWWbar[279]=armWWbar[278] - 11./3. + armWWbar[273];
   armWWbar[280]=armWWbar[37]*armWWbar[279];
   armWWbar[281]=armWWbar[3]*armWWbar[38]*armWWbar[39];
   armWWbar[275]=36*armWWbar[275];
   armWWbar[281]=armWWbar[275] + 1./18.*armWWbar[280] + 9*armWWbar[281]
   ;
   armWWbar[281]=MMt*armWWbar[281];
   armWWbar[279]=armWWbar[38]*armWWbar[279];
   armWWbar[279]=1./18.*armWWbar[279] + armWWbar[281];
   armWWbar[279]=armWWbar[255]*MMt*armWWbar[279];
   armWWbar[279]=armWWbar[266] + armWWbar[279];
   armWWbar[279]=armWWbar[34]*armWWbar[255]*armWWbar[279];
   armWWbar[130]=armWWbar[130] + armWWbar[279];
   armWWbar[130]=armWWbar[34]*armWWbar[130];
   armWWbar[135]=armWWbar[137] + 89./12. + armWWbar[135];
   armWWbar[118]=armWWbar[133] + armWWbar[166] + armWWbar[127] + 
   armWWbar[164] + armWWbar[142] + armWWbar[118] + armWWbar[162] + 
   armWWbar[168] + armWWbar[126] + armWWbar[159] + armWWbar[158] + 
   armWWbar[157] + armWWbar[156] + armWWbar[155] + armWWbar[154] + 
   armWWbar[153] + armWWbar[152] + armWWbar[151] + armWWbar[150] + 
   armWWbar[149] + armWWbar[148] + armWWbar[147] + armWWbar[167] + 
   armWWbar[146] + armWWbar[145] + armWWbar[144] + 1./3.*armWWbar[135]
    + armWWbar[143];
   armWWbar[118]=MMt*armWWbar[118];
   armWWbar[117]= - 1771./12. + armWWbar[117];
   armWWbar[116]=armWWbar[180] + armWWbar[178] + armWWbar[176] + 
   armWWbar[175] + armWWbar[121] + armWWbar[116] + armWWbar[174] + 
   armWWbar[182] + armWWbar[119] + armWWbar[173] + armWWbar[172] + 
   armWWbar[171] + armWWbar[181] + 1./8.*armWWbar[117] + armWWbar[170];
   armWWbar[116]=armWWbar[38]*armWWbar[116];
   armWWbar[117]=armWWbar[193] + 4175./12.*armWWbar[42];
   armWWbar[117]=1./8.*armWWbar[117] + armWWbar[200];
   armWWbar[119]=1733 + 2129*armWWbar[36];
   armWWbar[119]=1./27.*armWWbar[119] + armWWbar[202];
   armWWbar[119]=armWWbar[20]*armWWbar[119];
   armWWbar[121]=5111./16. - 203*armWWbar[36];
   armWWbar[121]=armWWbar[19]*armWWbar[121];
   armWWbar[116]=armWWbar[118] + armWWbar[184] + armWWbar[205] + 
   armWWbar[116] + 1./4.*armWWbar[119] + armWWbar[213] + 1./36.*
   armWWbar[121] + armWWbar[204] - 253./72.*armWWbar[23] + 
   armWWbar[211] - 6571./576.*armWWbar[22] - 683./288.*armWWbar[51] + 
   armWWbar[210] + armWWbar[209] + armWWbar[208] + 1./3.*armWWbar[117]
    + armWWbar[214];
   armWWbar[116]=MMt*armWWbar[116];
   armWWbar[117]=armWWbar[34]*armWWbar[266];
   armWWbar[118]=2539./8.*armWWbar[38] + 2819./3.*armWWbar[20] + 
   armWWbar[263] + 2797./8.*armWWbar[19];
   armWWbar[118]=armWWbar[38]*armWWbar[118];
   armWWbar[116]=armWWbar[117] + armWWbar[116] + 1./36.*armWWbar[118]
    + armWWbar[259];
   armWWbar[116]=armWWbar[34]*armWWbar[116];
   armWWbar[117]= - 3*armWWbar[47];
   armWWbar[118]=2015./36. + armWWbar[117];
   armWWbar[119]=11./2.*armWWbar[13];
   armWWbar[118]=armWWbar[119] + armWWbar[124] + 1./4.*armWWbar[118] + 
   armWWbar[125];
   armWWbar[121]=1./3. - armWWbar[7];
   armWWbar[126]=1./3.*armWWbar[1]*armWWbar[121];
   armWWbar[127]=armWWbar[40]*armWWbar[121];
   armWWbar[133]= - 1./2.*armWWbar[11];
   armWWbar[118]=armWWbar[133] + armWWbar[218] + armWWbar[127] + 
   armWWbar[126] + 59./8.*armWWbar[12] + 1./2.*armWWbar[118] + 
   armWWbar[10];
   armWWbar[118]=armWWbar[37]*armWWbar[118];
   armWWbar[135]= - 3./2.*armWWbar[37];
   armWWbar[137]=armWWbar[223] + armWWbar[135];
   armWWbar[137]=armWWbar[20]*armWWbar[137];
   armWWbar[142]= - 3*armWWbar[19];
   armWWbar[143]= - armWWbar[18] + armWWbar[142];
   armWWbar[143]=armWWbar[37]*armWWbar[143];
   armWWbar[132]=armWWbar[134] + 3*armWWbar[137] + armWWbar[132] + 1./2.
   *armWWbar[143];
   armWWbar[132]=armWWbar[3]*armWWbar[132];
   armWWbar[134]=armWWbar[56] - armWWbar[54];
   armWWbar[134]=armWWbar[227] + armWWbar[136] + 1./2.*armWWbar[134] - 
   armWWbar[55];
   armWWbar[134]=MMt*armWWbar[134];
   armWWbar[136]= - 5./8.*armWWbar[66] + 5*armWWbar[67] - 947./24. + 13
   *armWWbar[61];
   armWWbar[136]=armWWbar[167] + 1./2.*armWWbar[136] + armWWbar[63];
   armWWbar[137]= - 1./4.*armWWbar[69];
   armWWbar[118]=armWWbar[134] + armWWbar[165] + armWWbar[132] + 
   armWWbar[236] + armWWbar[140] + armWWbar[118] + armWWbar[161] + 
   armWWbar[235] + armWWbar[220] + armWWbar[131] + armWWbar[234] - 17./
   4.*armWWbar[17] + armWWbar[233] + armWWbar[237] + armWWbar[232] + 
   armWWbar[231] - 71./32.*armWWbar[50] + 195./32.*armWWbar[71] + 
   armWWbar[137] + 1./2.*armWWbar[30] + 3./8.*armWWbar[33] + 1./2.*
   armWWbar[136] + armWWbar[230];
   armWWbar[118]=MMt*armWWbar[118];
   armWWbar[132]=2585./36. + armWWbar[117];
   armWWbar[134]=11./4.*armWWbar[13];
   armWWbar[126]=armWWbar[179] + armWWbar[240] - 13./32.*armWWbar[37]
    + 3./2.*armWWbar[11] + armWWbar[127] + armWWbar[126] + 93./8.*
   armWWbar[12] - armWWbar[10] + armWWbar[134] + armWWbar[239] + 
   armWWbar[238] + armWWbar[228] + armWWbar[221] + 1./8.*armWWbar[132]
    + armWWbar[194];
   armWWbar[126]=armWWbar[38]*armWWbar[126];
   armWWbar[127]= - 1567./36. + armWWbar[183];
   armWWbar[127]= - 19./9.*armWWbar[36] + armWWbar[160] + armWWbar[189]
    + armWWbar[125] + armWWbar[188] + armWWbar[187] + armWWbar[186] + 
   armWWbar[185] + 1./4.*armWWbar[127] + armWWbar[191];
   armWWbar[132]=1./3. + 3./2.*armWWbar[45];
   armWWbar[132]=5./12.*armWWbar[11] + 1./2.*armWWbar[132] + 
   armWWbar[196];
   armWWbar[132]=armWWbar[37]*armWWbar[132];
   armWWbar[136]= - 19./18.*armWWbar[36] - 20./9. + armWWbar[243];
   armWWbar[136]=armWWbar[11]*armWWbar[136];
   armWWbar[127]=1./2.*armWWbar[132] + armWWbar[136] + 1./4.*
   armWWbar[127] + armWWbar[199];
   armWWbar[127]=MMH*armWWbar[127];
   armWWbar[132]= - 4*armWWbar[19];
   armWWbar[136]=3./2.*armWWbar[18];
   armWWbar[140]= - 9./2.*armWWbar[38] + 7./2.*armWWbar[20] + 
   armWWbar[136] + armWWbar[132];
   armWWbar[140]=armWWbar[3]*armWWbar[38]*armWWbar[140];
   armWWbar[143]= - 13*armWWbar[74] + 45./8.*armWWbar[42];
   armWWbar[144]=19./3.*armWWbar[36] + 233./12. + armWWbar[160];
   armWWbar[144]=armWWbar[18]*armWWbar[144];
   armWWbar[145]= - 257./8. - 113*armWWbar[36];
   armWWbar[145]=armWWbar[19]*armWWbar[145];
   armWWbar[146]=armWWbar[244] + 743./12.*armWWbar[19];
   armWWbar[146]=armWWbar[37]*armWWbar[146];
   armWWbar[147]= - 197 - 49./2.*armWWbar[36];
   armWWbar[147]=1./3.*armWWbar[147] + 23./2.*armWWbar[37];
   armWWbar[147]=armWWbar[20]*armWWbar[147];
   armWWbar[118]=armWWbar[118] + armWWbar[127] + armWWbar[140] + 
   armWWbar[126] + 1./12.*armWWbar[147] + 1./16.*armWWbar[146] + 1./24.
   *armWWbar[145] + 1./4.*armWWbar[144] + 7./8.*armWWbar[23] + 
   armWWbar[253] + 159./64.*armWWbar[22] + 363./32.*armWWbar[51] + 17./
   8.*armWWbar[43] - 1./4.*armWWbar[8] + armWWbar[252] + armWWbar[251]
    + 1./4.*armWWbar[143] - 2*armWWbar[72];
   armWWbar[118]=MMt*armWWbar[118];
   armWWbar[126]= - 1./2.*armWWbar[6];
   armWWbar[127]=armWWbar[126] - 19 + armWWbar[270];
   armWWbar[127]=armWWbar[18]*armWWbar[127];
   armWWbar[140]= - 1 - armWWbar[36];
   armWWbar[143]=armWWbar[140] + armWWbar[6];
   armWWbar[143]=1./6.*armWWbar[143] - armWWbar[37];
   armWWbar[143]=armWWbar[20]*armWWbar[143];
   armWWbar[127]=1./3.*armWWbar[143] + armWWbar[264] + 1./2.*
   armWWbar[258] + 1./9.*armWWbar[127];
   armWWbar[143]=1 + armWWbar[36];
   armWWbar[144]=armWWbar[143] - armWWbar[6];
   armWWbar[144]=armWWbar[11]*armWWbar[144];
   armWWbar[144]=13./2.*armWWbar[260] + 1./3.*armWWbar[144];
   armWWbar[144]=1./2.*armWWbar[144] + armWWbar[265];
   armWWbar[144]=MMH*armWWbar[144];
   armWWbar[145]= - 11./3.*armWWbar[11] - 29./6. + armWWbar[10];
   armWWbar[145]=armWWbar[38]*armWWbar[145];
   armWWbar[127]=1./6.*armWWbar[144] + 1./2.*armWWbar[127] + 1./3.*
   armWWbar[145];
   armWWbar[127]=MMH*armWWbar[127];
   armWWbar[144]=1./3.*armWWbar[6];
   armWWbar[145]=armWWbar[135] + armWWbar[144] - 3 + armWWbar[272];
   armWWbar[145]=armWWbar[37]*armWWbar[145];
   armWWbar[146]= - 19*armWWbar[36];
   armWWbar[147]=armWWbar[268] + armWWbar[146] - 16 + armWWbar[138];
   armWWbar[147]=armWWbar[3]*armWWbar[38]*armWWbar[147];
   armWWbar[148]=8 + 19./2.*armWWbar[36];
   armWWbar[145]=armWWbar[275] + armWWbar[147] + 1./3.*armWWbar[148] + 
   1./2.*armWWbar[145];
   armWWbar[145]=MMt*armWWbar[145];
   armWWbar[147]= - 17./6.*armWWbar[37] + 1./6.*armWWbar[6] + 17./3. + 
   3*armWWbar[36];
   armWWbar[147]=armWWbar[38]*armWWbar[147];
   armWWbar[148]= - armWWbar[3]*armWWbar[276];
   armWWbar[145]=armWWbar[145] + armWWbar[147] + 16*armWWbar[148];
   armWWbar[145]=MMt*armWWbar[145];
   armWWbar[145]=8./3.*armWWbar[276] + armWWbar[145];
   armWWbar[145]=armWWbar[34]*armWWbar[145];
   armWWbar[147]=7*armWWbar[18];
   armWWbar[148]=armWWbar[147] - 223./12.*armWWbar[19];
   armWWbar[148]= - 27./32.*armWWbar[38] + 1./8.*armWWbar[148] - 20./9.
   *armWWbar[20];
   armWWbar[148]=armWWbar[38]*armWWbar[148];
   armWWbar[118]=armWWbar[145] + armWWbar[118] + armWWbar[148] + 
   armWWbar[127];
   armWWbar[118]=armWWbar[34]*armWWbar[118];
   armWWbar[127]=1./4.*armWWbar[109];
   armWWbar[145]=1./2.*armWWbar[85];
   armWWbar[148]=1./2.*armWWbar[82] + armWWbar[145] + armWWbar[81] + 
   armWWbar[127];
   armWWbar[148]=MMH*armWWbar[148];
   armWWbar[149]= - 1025./18. + 7*armWWbar[95];
   armWWbar[149]= - 13./6.*armWWbar[98] - 13./6.*armWWbar[92] - 13./6.*
   armWWbar[93] + 1./6.*armWWbar[97] - 1./3.*armWWbar[99] + 13./9.*
   armWWbar[101] + 1./4.*armWWbar[149] - 1./3.*armWWbar[94];
   armWWbar[150]=1./3.*armWWbar[91];
   armWWbar[151]= - 11./2.*armWWbar[10] - armWWbar[12];
   armWWbar[151]=armWWbar[12]*armWWbar[151];
   armWWbar[152]= - 3./16.*armWWbar[46];
   armWWbar[153]=13./36.*armWWbar[11] - 59./48.*armWWbar[12] - 5./18.*
   armWWbar[10] - 11./24.*armWWbar[13] + armWWbar[152] - 3./8.*
   armWWbar[48] - 53./27. - 9./16.*armWWbar[44];
   armWWbar[153]=armWWbar[11]*armWWbar[153];
   armWWbar[154]=9./16.*armWWbar[44];
   armWWbar[155]=1./8.*armWWbar[48];
   armWWbar[156]=1./16.*armWWbar[46];
   armWWbar[148]=1./6.*armWWbar[148] + armWWbar[153] + 1./24.*
   armWWbar[151] + 5./72.*armWWbar[10] + armWWbar[156] - 17./48.*
   armWWbar[16] + armWWbar[155] + 5./4.*armWWbar[17] + armWWbar[154] - 
   85./12.*armWWbar[100] + 1./12.*armWWbar[96] + 1./4.*armWWbar[149] + 
   armWWbar[150];
   armWWbar[148]=MMH*armWWbar[148];
   armWWbar[149]=9./2.*armWWbar[44];
   armWWbar[151]= - 1./2.*armWWbar[46];
   armWWbar[153]= - 67./3.*armWWbar[12] - 11./3.*armWWbar[13] + 
   armWWbar[151] - armWWbar[48] + 1967./27. + armWWbar[149];
   armWWbar[153]=armWWbar[18]*armWWbar[153];
   armWWbar[157]=1./2.*armWWbar[46];
   armWWbar[158]=31./4.*armWWbar[12] + 59./9.*armWWbar[10] + 11./6.*
   armWWbar[13] + armWWbar[157] + 2765./108. + armWWbar[48];
   armWWbar[159]=43./9.*armWWbar[11];
   armWWbar[158]=1./2.*armWWbar[158] + armWWbar[159];
   armWWbar[158]=armWWbar[20]*armWWbar[158];
   armWWbar[161]=43./12.*armWWbar[79] + 5./2.*armWWbar[80] - 169./6.*
   armWWbar[106] - 3./2.*armWWbar[77] + 7./9.*armWWbar[105];
   armWWbar[161]=25./12.*armWWbar[78] + 1./2.*armWWbar[161] - 1./3.*
   armWWbar[107];
   armWWbar[162]=41./8.*armWWbar[12] + 85./3. - 19./4.*armWWbar[10];
   armWWbar[162]=armWWbar[19]*armWWbar[162];
   armWWbar[164]=1./4.*armWWbar[18] + 34./3.*armWWbar[19];
   armWWbar[164]=armWWbar[11]*armWWbar[164];
   armWWbar[148]=armWWbar[148] + 1./2.*armWWbar[158] + 1./3.*
   armWWbar[164] + 1./6.*armWWbar[162] + 1./8.*armWWbar[153] + 143./36.
   *armWWbar[23] - 127./144.*armWWbar[21] + 1./2.*armWWbar[161] + 14./3.
   *armWWbar[22];
   armWWbar[148]=MMH*armWWbar[148];
   armWWbar[153]=armWWbar[147] + 2759./4.*armWWbar[19];
   armWWbar[153]=armWWbar[19]*armWWbar[153];
   armWWbar[158]= - 5*armWWbar[18];
   armWWbar[161]=3725./2.*armWWbar[20] + armWWbar[158] + 18733./4.*
   armWWbar[19];
   armWWbar[161]=armWWbar[20]*armWWbar[161];
   armWWbar[162]=pow(armWWbar[18],2);
   armWWbar[153]=1./3.*armWWbar[161] - 1./3.*armWWbar[162] + 
   armWWbar[153];
   armWWbar[148]=1./12.*armWWbar[153] + armWWbar[148];
   armWWbar[148]=armWWbar[76]*armWWbar[148];
   armWWbar[153]= - 3*armWWbar[10];
   armWWbar[161]=59./9. - 33./2.*armWWbar[13];
   armWWbar[161]=1./2.*armWWbar[161];
   armWWbar[164]=33./2.*armWWbar[12];
   armWWbar[165]=armWWbar[164] + armWWbar[161] + armWWbar[153];
   armWWbar[166]=1./6. - armWWbar[7];
   armWWbar[167]=armWWbar[166] + armWWbar[115];
   armWWbar[168]=1./3.*armWWbar[1]*armWWbar[167];
   armWWbar[167]=armWWbar[40]*armWWbar[167];
   armWWbar[165]=armWWbar[11] + armWWbar[167] + 1./2.*armWWbar[165] + 
   armWWbar[168];
   armWWbar[165]=armWWbar[38]*armWWbar[165];
   armWWbar[161]=armWWbar[164] + armWWbar[161] + armWWbar[10];
   armWWbar[161]= - armWWbar[11] + armWWbar[167] + 1./2.*armWWbar[161]
    + armWWbar[168];
   armWWbar[161]=armWWbar[37]*armWWbar[161];
   armWWbar[170]= - armWWbar[37]*armWWbar[19];
   armWWbar[171]= - armWWbar[20]*armWWbar[37];
   armWWbar[174]=1./2.*armWWbar[170] + armWWbar[171];
   armWWbar[174]=armWWbar[3]*armWWbar[174];
   armWWbar[174]=armWWbar[161] + 3./2.*armWWbar[174];
   armWWbar[174]=MMt*armWWbar[174];
   armWWbar[176]= - 5./2.*armWWbar[36];
   armWWbar[178]= - 17./3. + armWWbar[176];
   armWWbar[178]=armWWbar[19]*armWWbar[178];
   armWWbar[179]=armWWbar[18]*armWWbar[36];
   armWWbar[178]=armWWbar[179] + armWWbar[178];
   armWWbar[180]= - 1 + armWWbar[271];
   armWWbar[183]=armWWbar[11]*armWWbar[180];
   armWWbar[184]=1./2. - armWWbar[10];
   armWWbar[185]=1./2.*armWWbar[184] + armWWbar[11];
   armWWbar[185]=armWWbar[37]*armWWbar[185];
   armWWbar[183]=armWWbar[185] + armWWbar[183] - 1./4.*armWWbar[36] + 
   armWWbar[10];
   armWWbar[183]=MMH*armWWbar[183];
   armWWbar[186]=17 - 13./2.*armWWbar[36];
   armWWbar[186]=1./2.*armWWbar[186] - armWWbar[37];
   armWWbar[186]=armWWbar[20]*armWWbar[186];
   armWWbar[187]=31./4.*armWWbar[19] - 10*armWWbar[20];
   armWWbar[187]=armWWbar[3]*armWWbar[38]*armWWbar[187];
   armWWbar[188]= - 1./4.*armWWbar[18];
   armWWbar[189]=armWWbar[188] + 4./3.*armWWbar[19];
   armWWbar[189]=armWWbar[37]*armWWbar[189];
   armWWbar[165]=armWWbar[174] + 1./3.*armWWbar[183] + armWWbar[187] + 
   armWWbar[165] + 1./6.*armWWbar[186] + 1./4.*armWWbar[178] + 
   armWWbar[189];
   armWWbar[165]=MMt*armWWbar[165];
   armWWbar[174]=3./2.*armWWbar[36];
   armWWbar[178]= - 3*armWWbar[37];
   armWWbar[183]=armWWbar[178] + armWWbar[115] + 1./3. + armWWbar[174];
   armWWbar[183]=armWWbar[38]*armWWbar[183];
   armWWbar[186]=1./4.*armWWbar[6];
   armWWbar[187]=armWWbar[269] - 1./3. + armWWbar[186];
   armWWbar[187]=armWWbar[37]*armWWbar[187];
   armWWbar[193]= - armWWbar[36] + armWWbar[37];
   armWWbar[193]=armWWbar[3]*armWWbar[38]*armWWbar[193];
   armWWbar[187]=3*armWWbar[193] + armWWbar[270] + armWWbar[187];
   armWWbar[187]=MMt*armWWbar[187];
   armWWbar[183]=1./2.*armWWbar[183] + armWWbar[187];
   armWWbar[183]=armWWbar[34]*MMt*armWWbar[183];
   armWWbar[187]= - armWWbar[36] - armWWbar[6];
   armWWbar[193]=armWWbar[18]*armWWbar[187];
   armWWbar[200]=armWWbar[18] + armWWbar[19];
   armWWbar[202]=armWWbar[37]*armWWbar[200];
   armWWbar[204]= - 1./3.*armWWbar[19];
   armWWbar[202]=armWWbar[202] + 1./2.*armWWbar[193] + armWWbar[204];
   armWWbar[205]=armWWbar[115] + 1./3. + armWWbar[270];
   armWWbar[208]=1./2.*armWWbar[205] - armWWbar[37];
   armWWbar[208]=armWWbar[20]*armWWbar[208];
   armWWbar[202]=1./2.*armWWbar[202] + armWWbar[208];
   armWWbar[208]= - 1./3. + armWWbar[271];
   armWWbar[209]=armWWbar[208] + armWWbar[126];
   armWWbar[210]=armWWbar[11]*armWWbar[209];
   armWWbar[211]=armWWbar[199] + armWWbar[210];
   armWWbar[213]=armWWbar[182] + armWWbar[11];
   armWWbar[214]=armWWbar[37]*armWWbar[213];
   armWWbar[211]=1./2.*armWWbar[211] + armWWbar[214];
   armWWbar[211]=MMH*armWWbar[211];
   armWWbar[214]=armWWbar[10] - armWWbar[11];
   armWWbar[214]=armWWbar[38]*armWWbar[214];
   armWWbar[202]=1./2.*armWWbar[211] + 1./2.*armWWbar[202] + 
   armWWbar[214];
   armWWbar[202]=MMH*armWWbar[202];
   armWWbar[211]= - armWWbar[19] + armWWbar[20];
   armWWbar[211]=armWWbar[38]*armWWbar[211];
   armWWbar[202]=17./4.*armWWbar[211] + armWWbar[202];
   armWWbar[165]=armWWbar[183] + 1./3.*armWWbar[202] + armWWbar[165];
   armWWbar[165]=armWWbar[34]*armWWbar[165];
   armWWbar[183]=armWWbar[7] - armWWbar[6];
   armWWbar[202]=armWWbar[1]*armWWbar[183];
   armWWbar[211]=1./3.*armWWbar[202];
   armWWbar[183]=armWWbar[40]*armWWbar[183];
   armWWbar[214]=armWWbar[211] + armWWbar[183];
   armWWbar[214]=armWWbar[18]*armWWbar[214];
   armWWbar[220]= - 1./3. + armWWbar[7];
   armWWbar[221]=armWWbar[1]*armWWbar[220];
   armWWbar[223]=armWWbar[40]*armWWbar[220];
   armWWbar[221]=1./3.*armWWbar[221] + armWWbar[223];
   armWWbar[223]=armWWbar[19]*armWWbar[221];
   armWWbar[223]=armWWbar[214] + armWWbar[223];
   armWWbar[167]=armWWbar[168] + armWWbar[167];
   armWWbar[167]=armWWbar[20]*armWWbar[167];
   armWWbar[168]=armWWbar[10]*armWWbar[121];
   armWWbar[227]=armWWbar[1]*armWWbar[168];
   armWWbar[228]=armWWbar[40]*armWWbar[168];
   armWWbar[227]=1./3.*armWWbar[227] + armWWbar[228];
   armWWbar[228]=armWWbar[126] - 1./6. + armWWbar[7];
   armWWbar[231]=armWWbar[1]*armWWbar[228];
   armWWbar[228]=armWWbar[40]*armWWbar[228];
   armWWbar[228]=1./3.*armWWbar[231] + armWWbar[228];
   armWWbar[228]=armWWbar[11]*armWWbar[228];
   armWWbar[227]=1./2.*armWWbar[227] + armWWbar[228];
   armWWbar[227]=MMH*armWWbar[227];
   armWWbar[167]=armWWbar[227] + 1./2.*armWWbar[223] + armWWbar[167];
   armWWbar[167]=MMH*armWWbar[167];
   armWWbar[223]= - 11./4.*armWWbar[12];
   armWWbar[134]=armWWbar[223] + armWWbar[134] + armWWbar[199];
   armWWbar[134]=armWWbar[18]*armWWbar[134];
   armWWbar[227]= - 1./9. - 13./8.*armWWbar[10];
   armWWbar[228]= - 11./16.*armWWbar[12];
   armWWbar[227]=1./3.*armWWbar[227] + armWWbar[228];
   armWWbar[227]=armWWbar[19]*armWWbar[227];
   armWWbar[231]= - 47./27. + armWWbar[119];
   armWWbar[223]=5./9.*armWWbar[11] + armWWbar[223] + 1./4.*
   armWWbar[231] - 5./9.*armWWbar[10];
   armWWbar[223]=armWWbar[11]*armWWbar[223];
   armWWbar[231]=armWWbar[12]*armWWbar[10];
   armWWbar[232]=47./27.*armWWbar[10] + 11./2.*armWWbar[231];
   armWWbar[223]=1./4.*armWWbar[232] + armWWbar[223];
   armWWbar[223]=MMH*armWWbar[223];
   armWWbar[232]= - armWWbar[18] + 79./6.*armWWbar[19];
   armWWbar[232]=armWWbar[11]*armWWbar[232];
   armWWbar[234]= - 1./36.*armWWbar[11] + 11./8.*armWWbar[12] - 19./36.
   *armWWbar[10] + 1./27. - 11./16.*armWWbar[13];
   armWWbar[234]=armWWbar[20]*armWWbar[234];
   armWWbar[134]=1./2.*armWWbar[223] + armWWbar[234] + 1./12.*
   armWWbar[232] + 1./4.*armWWbar[134] + armWWbar[227];
   armWWbar[134]=MMH*armWWbar[134];
   armWWbar[223]= - 17*armWWbar[18] + 91./2.*armWWbar[19];
   armWWbar[223]=armWWbar[19]*armWWbar[223];
   armWWbar[227]= - 257./48.*armWWbar[20] + 17./8.*armWWbar[18] + 
   armWWbar[204];
   armWWbar[227]=armWWbar[20]*armWWbar[227];
   armWWbar[223]=1./8.*armWWbar[223] + armWWbar[227];
   armWWbar[134]=1./3.*armWWbar[223] + armWWbar[134];
   armWWbar[134]=armWWbar[76]*armWWbar[134];
   armWWbar[134]=armWWbar[165] + 1./6.*armWWbar[167] + armWWbar[134];
   armWWbar[165]=pow(armWWbar[5],2);
   armWWbar[134]=armWWbar[165]*armWWbar[134];
   armWWbar[167]= - armWWbar[1]*armWWbar[14];
   armWWbar[223]= - armWWbar[40]*armWWbar[14];
   armWWbar[167]=1./3.*armWWbar[167] + armWWbar[223];
   armWWbar[223]=1./2.*armWWbar[167];
   armWWbar[227]=1./6. + armWWbar[7];
   armWWbar[232]=armWWbar[1]*armWWbar[227];
   armWWbar[234]=armWWbar[40]*armWWbar[227];
   armWWbar[232]=1./3.*armWWbar[232] + armWWbar[234];
   armWWbar[232]=armWWbar[11]*armWWbar[232];
   armWWbar[232]=armWWbar[223] + armWWbar[232];
   armWWbar[232]=MMH*armWWbar[232];
   armWWbar[221]=armWWbar[18]*armWWbar[221];
   armWWbar[234]=armWWbar[1]*armWWbar[114];
   armWWbar[114]=armWWbar[40]*armWWbar[114];
   armWWbar[114]=1./3.*armWWbar[234] + armWWbar[114];
   armWWbar[114]=armWWbar[20]*armWWbar[114];
   armWWbar[114]=armWWbar[232] + armWWbar[221] + armWWbar[114];
   armWWbar[114]=MMH*armWWbar[114];
   armWWbar[114]=armWWbar[134] + armWWbar[118] + 1./6.*armWWbar[114] + 
   armWWbar[148];
   armWWbar[114]=armWWbar[165]*armWWbar[114];
   armWWbar[118]= - 2*armWWbar[89];
   armWWbar[134]=armWWbar[118] + armWWbar[86];
   armWWbar[148]= - 1./4.*armWWbar[82] + armWWbar[145] + armWWbar[127]
    + armWWbar[134] + armWWbar[81];
   armWWbar[148]=MMH*armWWbar[148];
   armWWbar[221]= - 9./2.*armWWbar[44];
   armWWbar[232]= - 3*armWWbar[48];
   armWWbar[234]=armWWbar[232] + 25./9. + armWWbar[221];
   armWWbar[234]= - 233./72.*armWWbar[12] - 13./18.*armWWbar[10] - 19./
   72.*armWWbar[13] + 1./2.*armWWbar[234] - armWWbar[46];
   armWWbar[234]=1./2.*armWWbar[234] + 4./9.*armWWbar[11];
   armWWbar[234]=armWWbar[11]*armWWbar[234];
   armWWbar[235]=2*armWWbar[102];
   armWWbar[236]= - 2315./192. + armWWbar[235];
   armWWbar[240]=43./48.*armWWbar[103];
   armWWbar[243]=5./2.*armWWbar[10];
   armWWbar[244]=armWWbar[243] - armWWbar[12];
   armWWbar[244]=armWWbar[12]*armWWbar[244];
   armWWbar[251]=9./8.*armWWbar[44];
   armWWbar[252]=1./4.*armWWbar[48];
   armWWbar[253]=1./8.*armWWbar[46];
   armWWbar[148]=1./3.*armWWbar[148] + armWWbar[234] + 1./12.*
   armWWbar[244] - 13./72.*armWWbar[10] + armWWbar[253] - 13./8.*
   armWWbar[16] + armWWbar[252] + 41./48.*armWWbar[17] + armWWbar[251]
    - 143./16.*armWWbar[100] + 11./48.*armWWbar[96] + armWWbar[240] + 2.
   /3.*armWWbar[91] + 1./12.*armWWbar[98] - 41./48.*armWWbar[92] - 41./
   48.*armWWbar[93] + 1./24.*armWWbar[97] + 1./12.*armWWbar[99] + 13./
   18.*armWWbar[101] - 1./6.*armWWbar[94] + 1./3.*armWWbar[236] + 7./8.
   *armWWbar[95];
   armWWbar[148]=MMH*armWWbar[148];
   armWWbar[234]=9*armWWbar[44];
   armWWbar[236]=877./9. + armWWbar[234];
   armWWbar[236]= - 175./18.*armWWbar[12] + armWWbar[196] - 19./36.*
   armWWbar[13] - armWWbar[46] + 1./2.*armWWbar[236] - armWWbar[48];
   armWWbar[236]=armWWbar[18]*armWWbar[236];
   armWWbar[244]=8./3.*armWWbar[108];
   armWWbar[258]=43 - 5*armWWbar[10];
   armWWbar[258]=7*armWWbar[258] - armWWbar[12];
   armWWbar[258]=armWWbar[19]*armWWbar[258];
   armWWbar[259]=armWWbar[18] + 53./12.*armWWbar[19];
   armWWbar[259]=armWWbar[11]*armWWbar[259];
   armWWbar[260]=31./6.*armWWbar[11] + 89./18.*armWWbar[12] + 107./18.*
   armWWbar[10] + 19./72.*armWWbar[13] + armWWbar[239] + 37./12. + 
   armWWbar[48];
   armWWbar[260]=armWWbar[20]*armWWbar[260];
   armWWbar[148]=armWWbar[148] + 1./2.*armWWbar[260] + 1./2.*
   armWWbar[259] + 1./24.*armWWbar[258] + 1./4.*armWWbar[236] + 25./9.*
   armWWbar[23] - 2./9.*armWWbar[21] + 5./3.*armWWbar[22] + 7./4.*
   armWWbar[78] - 5./12.*armWWbar[107] + 119./48.*armWWbar[79] + 41./48.
   *armWWbar[80] - 43./4.*armWWbar[106] + 7./18.*armWWbar[105] + 
   armWWbar[244] - 3./4.*armWWbar[77];
   armWWbar[148]=MMH*armWWbar[148];
   armWWbar[236]=29*armWWbar[18] + 665*armWWbar[19];
   armWWbar[236]=armWWbar[19]*armWWbar[236];
   armWWbar[236]=11./3.*armWWbar[162] + 1./4.*armWWbar[236];
   armWWbar[236]=1./2.*armWWbar[236];
   armWWbar[258]=199./8.*armWWbar[18];
   armWWbar[259]=1057./8.*armWWbar[20] + armWWbar[258] - 370*
   armWWbar[19];
   armWWbar[259]=armWWbar[20]*armWWbar[259];
   armWWbar[259]=armWWbar[236] + 1./3.*armWWbar[259];
   armWWbar[259]=1./3.*armWWbar[259] + armWWbar[148];
   armWWbar[259]=armWWbar[76]*armWWbar[259];
   armWWbar[260]=7./6. + armWWbar[7];
   armWWbar[260]=1./3.*armWWbar[260] + armWWbar[126];
   armWWbar[260]=armWWbar[1]*armWWbar[260];
   armWWbar[264]= - 11./18.*armWWbar[6];
   armWWbar[265]=armWWbar[264] + 47./54. + armWWbar[7];
   armWWbar[265]=armWWbar[40]*armWWbar[265];
   armWWbar[260]=armWWbar[260] + armWWbar[265];
   armWWbar[260]=armWWbar[11]*armWWbar[260];
   armWWbar[167]=armWWbar[167] + armWWbar[260];
   armWWbar[167]=MMH*armWWbar[167];
   armWWbar[227]=1./3.*armWWbar[227] + armWWbar[126];
   armWWbar[227]=armWWbar[1]*armWWbar[227];
   armWWbar[260]=armWWbar[264] - 7./54. + armWWbar[7];
   armWWbar[260]=armWWbar[40]*armWWbar[260];
   armWWbar[227]=armWWbar[227] + armWWbar[260];
   armWWbar[227]=armWWbar[18]*armWWbar[227];
   armWWbar[260]= - 7./6. - armWWbar[7];
   armWWbar[260]=1./3.*armWWbar[260] + armWWbar[115];
   armWWbar[260]=armWWbar[1]*armWWbar[260];
   armWWbar[120]=armWWbar[120] - 47./54. - armWWbar[7];
   armWWbar[120]=armWWbar[40]*armWWbar[120];
   armWWbar[120]=armWWbar[260] + armWWbar[120];
   armWWbar[120]=armWWbar[20]*armWWbar[120];
   armWWbar[120]=armWWbar[167] + armWWbar[227] + armWWbar[120];
   armWWbar[120]=1./6.*MMH*armWWbar[120];
   armWWbar[114]=armWWbar[114] + armWWbar[116] + armWWbar[120] + 
   armWWbar[259];
   armWWbar[114]=armWWbar[165]*armWWbar[114];
   armWWbar[116]= - 3./2.*armWWbar[48];
   armWWbar[167]= - 9./4.*armWWbar[44];
   armWWbar[227]=armWWbar[196] - 1./36.*armWWbar[13] - 5./4.*
   armWWbar[46] + armWWbar[116] + 251./27. + armWWbar[167];
   armWWbar[197]=armWWbar[197] + 1./2.*armWWbar[227] - 7./9.*
   armWWbar[12];
   armWWbar[197]=armWWbar[11]*armWWbar[197];
   armWWbar[227]= - 157./192. + armWWbar[235];
   armWWbar[235]= - 1./12.*armWWbar[10];
   armWWbar[259]=13./2.*armWWbar[10] - armWWbar[12];
   armWWbar[259]=armWWbar[12]*armWWbar[259];
   armWWbar[134]= - 1./2.*armWWbar[82] + 1./4.*armWWbar[85] + 1./8.*
   armWWbar[109] + armWWbar[134] + 1./2.*armWWbar[81];
   armWWbar[134]=MMH*armWWbar[134];
   armWWbar[134]=1./3.*armWWbar[134] + 1./2.*armWWbar[197] + 1./24.*
   armWWbar[259] + armWWbar[235] + armWWbar[156] - 61./48.*armWWbar[16]
    + armWWbar[155] + 5./48.*armWWbar[17] + armWWbar[154] - 137./48.*
   armWWbar[100] + 7./48.*armWWbar[96] + armWWbar[240] + armWWbar[150]
    + 1./8.*armWWbar[98] - 5./16.*armWWbar[92] - 5./16.*armWWbar[93] + 
   1./6.*armWWbar[99] + 13./36.*armWWbar[101] - 1./12.*armWWbar[94] + 1.
   /3.*armWWbar[227] + 7./16.*armWWbar[95];
   armWWbar[134]=MMH*armWWbar[134];
   armWWbar[150]= - 11./36.*armWWbar[12] + armWWbar[235] - 1./144.*
   armWWbar[13] + armWWbar[152] - 1./8.*armWWbar[48] + 92./27. + 
   armWWbar[154];
   armWWbar[150]=armWWbar[18]*armWWbar[150];
   armWWbar[152]=armWWbar[228] + 419./24. - 2*armWWbar[10];
   armWWbar[152]=armWWbar[19]*armWWbar[152];
   armWWbar[197]=armWWbar[212] - 17./6.*armWWbar[19];
   armWWbar[197]=armWWbar[11]*armWWbar[197];
   armWWbar[212]=1./36.*armWWbar[13] + armWWbar[46] - 1667./108. + 
   armWWbar[48];
   armWWbar[212]=55./36.*armWWbar[11] + 77./144.*armWWbar[12] + 1./4.*
   armWWbar[212] + 4./3.*armWWbar[10];
   armWWbar[212]=armWWbar[20]*armWWbar[212];
   armWWbar[134]=armWWbar[134] + armWWbar[212] + 1./12.*armWWbar[197]
    + 1./3.*armWWbar[152] + armWWbar[150] - 43./36.*armWWbar[23] - 49./
   144.*armWWbar[21] - 2./3.*armWWbar[22] + 17./24.*armWWbar[78] - 1./4.
   *armWWbar[107] + 19./12.*armWWbar[79] + 11./48.*armWWbar[80] - 89./
   24.*armWWbar[106] + 7./36.*armWWbar[105] + armWWbar[244] - 3./8.*
   armWWbar[77];
   armWWbar[134]=MMH*armWWbar[134];
   armWWbar[150]=armWWbar[18] - 49./6.*armWWbar[19];
   armWWbar[150]=armWWbar[19]*armWWbar[150];
   armWWbar[150]=23./9.*armWWbar[162] + 5./2.*armWWbar[150];
   armWWbar[152]=65*armWWbar[18] - 2957./2.*armWWbar[19];
   armWWbar[152]=1./4.*armWWbar[152] - 421*armWWbar[20];
   armWWbar[152]=armWWbar[20]*armWWbar[152];
   armWWbar[150]=1./2.*armWWbar[150] + 1./9.*armWWbar[152];
   armWWbar[134]=1./2.*armWWbar[150] + armWWbar[134];
   armWWbar[134]=armWWbar[76]*armWWbar[134];
   armWWbar[150]= - 1./2.*armWWbar[97];
   armWWbar[152]=armWWbar[150] + 233./48. - armWWbar[99];
   armWWbar[197]=1./4.*armWWbar[96];
   armWWbar[152]= - 11./12.*armWWbar[17] + 77./36.*armWWbar[100] + 
   armWWbar[197] + 13./4.*armWWbar[103] + 5./4.*armWWbar[92] + 1./3.*
   armWWbar[152] + 5./4.*armWWbar[93];
   armWWbar[212]= - 1./3.*armWWbar[16];
   armWWbar[227]=7./18.*armWWbar[12] - 1./3. - armWWbar[46];
   armWWbar[227]=armWWbar[11]*armWWbar[227];
   armWWbar[228]=MMH*armWWbar[82];
   armWWbar[152]=1./12.*armWWbar[228] + 1./8.*armWWbar[227] + 1./48.*
   armWWbar[231] + 1./4.*armWWbar[152] + armWWbar[212];
   armWWbar[152]=MMH*armWWbar[152];
   armWWbar[227]= - 1./6. - armWWbar[10];
   armWWbar[228]=1./4.*armWWbar[12];
   armWWbar[227]=1./3.*armWWbar[227] + armWWbar[228];
   armWWbar[227]=armWWbar[19]*armWWbar[227];
   armWWbar[231]=1./4.*armWWbar[46];
   armWWbar[235]=1./6.*armWWbar[10];
   armWWbar[240]= - 13./18.*armWWbar[12] + armWWbar[235] - 1./9. + 
   armWWbar[231];
   armWWbar[240]=armWWbar[20]*armWWbar[240];
   armWWbar[244]=43./18.*armWWbar[12] - 47./9. - armWWbar[46];
   armWWbar[244]=armWWbar[18]*armWWbar[244];
   armWWbar[152]=armWWbar[152] + 1./2.*armWWbar[240] + 1./4.*
   armWWbar[227] + 1./8.*armWWbar[244] + 4./9.*armWWbar[23] + 1./8.*
   armWWbar[21] - 1./18.*armWWbar[22] + 1./6.*armWWbar[78] - 1./12.*
   armWWbar[107] + 1./48.*armWWbar[79] + 2./9.*armWWbar[106] - 5./16.*
   armWWbar[80];
   armWWbar[152]=MMH*armWWbar[152];
   armWWbar[227]=1./3.*armWWbar[18] - 41./8.*armWWbar[19];
   armWWbar[227]=armWWbar[19]*armWWbar[227];
   armWWbar[240]= - armWWbar[18] - 203./8.*armWWbar[19];
   armWWbar[240]=1./3.*armWWbar[240] - 425./16.*armWWbar[20];
   armWWbar[240]=armWWbar[20]*armWWbar[240];
   armWWbar[227]=1./2.*armWWbar[227] + armWWbar[240];
   armWWbar[152]=1./3.*armWWbar[227] + armWWbar[152];
   armWWbar[152]=armWWbar[76]*armWWbar[152];
   armWWbar[227]= - 1./2.*armWWbar[80];
   armWWbar[240]=armWWbar[106] + armWWbar[227];
   armWWbar[244]=1./2.*armWWbar[21];
   armWWbar[259]=armWWbar[240] + armWWbar[244];
   armWWbar[260]= - 1 + armWWbar[228];
   armWWbar[260]=armWWbar[18]*armWWbar[260];
   armWWbar[264]= - 1 - armWWbar[12];
   armWWbar[264]=armWWbar[20]*armWWbar[264];
   armWWbar[265]= - 1./2.*armWWbar[19];
   armWWbar[259]=1./4.*armWWbar[264] + armWWbar[265] + 1./2.*
   armWWbar[259] + armWWbar[260];
   armWWbar[260]=1./2.*armWWbar[93];
   armWWbar[264]=1./2.*armWWbar[92];
   armWWbar[266]=armWWbar[264] + 1 + armWWbar[260];
   armWWbar[269]= - 1./6.*armWWbar[17];
   armWWbar[266]= - 1./6.*armWWbar[16] + armWWbar[269] + 1./2.*
   armWWbar[100] + 1./3.*armWWbar[266] + 1./2.*armWWbar[103];
   armWWbar[266]=MMH*armWWbar[266];
   armWWbar[259]=1./3.*armWWbar[259] + 1./2.*armWWbar[266];
   armWWbar[259]=MMH*armWWbar[259];
   armWWbar[266]= - 5*armWWbar[19];
   armWWbar[272]= - 29*armWWbar[20];
   armWWbar[275]=armWWbar[266] + armWWbar[272];
   armWWbar[275]=armWWbar[20]*armWWbar[275];
   armWWbar[279]=pow(armWWbar[19],2);
   armWWbar[275]= - armWWbar[279] + 1./2.*armWWbar[275];
   armWWbar[259]=1./6.*armWWbar[275] + armWWbar[259];
   armWWbar[259]=armWWbar[255]*armWWbar[76]*armWWbar[259];
   armWWbar[152]=armWWbar[152] + 1./4.*armWWbar[259];
   armWWbar[152]=armWWbar[255]*armWWbar[152];
   armWWbar[259]=1./3. + armWWbar[126];
   armWWbar[259]=armWWbar[1]*armWWbar[259];
   armWWbar[198]=19./3. + armWWbar[198];
   armWWbar[198]=armWWbar[40]*armWWbar[198];
   armWWbar[198]=armWWbar[259] + 1./9.*armWWbar[198];
   armWWbar[198]=armWWbar[11]*armWWbar[198];
   armWWbar[198]=armWWbar[223] + armWWbar[198];
   armWWbar[198]=MMH*armWWbar[198];
   armWWbar[223]=1./3. - armWWbar[6];
   armWWbar[259]=armWWbar[1]*armWWbar[223];
   armWWbar[223]=armWWbar[40]*armWWbar[223];
   armWWbar[275]=11./9.*armWWbar[223];
   armWWbar[281]=armWWbar[259] + armWWbar[275];
   armWWbar[282]=armWWbar[18]*armWWbar[281];
   armWWbar[283]= - 1./3. + armWWbar[115];
   armWWbar[283]=armWWbar[1]*armWWbar[283];
   armWWbar[262]= - 19./3. + armWWbar[262];
   armWWbar[262]=armWWbar[40]*armWWbar[262];
   armWWbar[262]=armWWbar[283] + 1./9.*armWWbar[262];
   armWWbar[262]=armWWbar[20]*armWWbar[262];
   armWWbar[198]=armWWbar[198] + 1./2.*armWWbar[282] + armWWbar[262];
   armWWbar[198]=MMH*armWWbar[198];
   armWWbar[134]=armWWbar[152] + 1./6.*armWWbar[198] + armWWbar[134];
   armWWbar[134]=armWWbar[255]*armWWbar[134];
   armWWbar[152]=1633./8.*armWWbar[20] + armWWbar[258] + 218*
   armWWbar[19];
   armWWbar[152]=armWWbar[20]*armWWbar[152];
   armWWbar[152]=armWWbar[236] + 1./3.*armWWbar[152];
   armWWbar[148]=1./3.*armWWbar[152] + armWWbar[148];
   armWWbar[148]=armWWbar[76]*armWWbar[148];
   armWWbar[120]=armWWbar[134] + armWWbar[120] + armWWbar[148];
   armWWbar[120]=armWWbar[255]*armWWbar[120];
   armWWbar[134]=15./2.*armWWbar[41];
   armWWbar[148]= - 27./4.*armWWbar[73];
   armWWbar[152]= - 17*armWWbar[51];
   armWWbar[198]= - 7./2.*armWWbar[21];
   armWWbar[236]= - 39./4.*armWWbar[23] + armWWbar[198] + 85./12.*
   armWWbar[22] + armWWbar[152] + 57./4.*armWWbar[43] + armWWbar[148]
    + armWWbar[134] - 73./4.*armWWbar[74] + armWWbar[72];
   armWWbar[258]=armWWbar[181] + 5./8. + armWWbar[125];
   armWWbar[258]=armWWbar[37]*armWWbar[258];
   armWWbar[262]= - 1./2.*armWWbar[68];
   armWWbar[282]= - armWWbar[37]*armWWbar[18];
   armWWbar[283]=armWWbar[282] - armWWbar[21] + armWWbar[18];
   armWWbar[283]=armWWbar[3]*armWWbar[283];
   armWWbar[284]= - MMH*armWWbar[55];
   armWWbar[258]=armWWbar[284] + 1./2.*armWWbar[283] + armWWbar[258] + 
   armWWbar[262] - 15./8.*armWWbar[70] - 2*armWWbar[50] - 9./2.*
   armWWbar[71] + 5./4.*armWWbar[65] + 3./4.*armWWbar[67] + 9./8. + 
   armWWbar[53];
   armWWbar[258]=MMt*armWWbar[258];
   armWWbar[195]=armWWbar[195] + 7./6.*armWWbar[11];
   armWWbar[195]=armWWbar[37]*armWWbar[195];
   armWWbar[285]= - 7./2. - 5*armWWbar[65];
   armWWbar[285]= - 9./2.*armWWbar[16] - 11./4.*armWWbar[68] + 9./4.*
   armWWbar[70] + 5./4.*armWWbar[50] + 3./2.*armWWbar[60] + 1./4.*
   armWWbar[285] + armWWbar[191];
   armWWbar[286]= - 1./2. + armWWbar[131];
   armWWbar[286]=armWWbar[11]*armWWbar[286];
   armWWbar[285]=1./2.*armWWbar[285] + armWWbar[286];
   armWWbar[195]=armWWbar[285] + 1./2.*armWWbar[195];
   armWWbar[195]=MMH*armWWbar[195];
   armWWbar[286]= - 1./4.*armWWbar[49];
   armWWbar[287]=armWWbar[286] + 13./8. - armWWbar[45];
   armWWbar[288]= - 5./4.*armWWbar[37];
   armWWbar[287]=3*armWWbar[287] + armWWbar[288];
   armWWbar[287]=armWWbar[38]*armWWbar[287];
   armWWbar[289]=1./2.*armWWbar[18];
   armWWbar[290]=armWWbar[289] + 17./3.*armWWbar[19];
   armWWbar[290]=armWWbar[37]*armWWbar[290];
   armWWbar[291]=armWWbar[37] + 41./24. + armWWbar[160];
   armWWbar[291]=armWWbar[20]*armWWbar[291];
   armWWbar[292]=3./4. - armWWbar[39];
   armWWbar[292]=armWWbar[18]*armWWbar[292];
   armWWbar[293]=3*armWWbar[292];
   armWWbar[195]=armWWbar[258] + armWWbar[195] + armWWbar[287] + 
   armWWbar[291] + 1./4.*armWWbar[290] + 121./12.*armWWbar[19] + 1./2.*
   armWWbar[236] + armWWbar[293];
   armWWbar[195]=MMt*armWWbar[195];
   armWWbar[236]=1./3.*armWWbar[37];
   armWWbar[190]=armWWbar[236] + armWWbar[190] - 55./27. + 
   armWWbar[131];
   armWWbar[190]=armWWbar[20]*armWWbar[190];
   armWWbar[258]=3./2.*armWWbar[39];
   armWWbar[287]= - 7./108.*armWWbar[36] - 67./27. + armWWbar[258];
   armWWbar[287]=armWWbar[18]*armWWbar[287];
   armWWbar[290]= - 3./2.*armWWbar[23] + 3./2.*armWWbar[21] + 7./4.*
   armWWbar[51] + armWWbar[225] + 5./2.*armWWbar[73] + armWWbar[72] - 9.
   /4.*armWWbar[41];
   armWWbar[282]=1./6.*armWWbar[282];
   armWWbar[190]=1./2.*armWWbar[190] + armWWbar[282] + armWWbar[290] + 
   armWWbar[287];
   armWWbar[287]= - 1./2.*armWWbar[60];
   armWWbar[291]=3./2.*armWWbar[16];
   armWWbar[294]=armWWbar[291] + 7./4.*armWWbar[68] + armWWbar[287] + 7.
   /4. - armWWbar[62];
   armWWbar[295]= - 7./54.*armWWbar[36] + 55./27. + armWWbar[160];
   armWWbar[295]=armWWbar[11]*armWWbar[295];
   armWWbar[296]= - armWWbar[37]*armWWbar[11];
   armWWbar[297]=1./6.*armWWbar[296];
   armWWbar[295]=armWWbar[297] + armWWbar[294] + 1./2.*armWWbar[295];
   armWWbar[295]=MMH*armWWbar[295];
   armWWbar[298]= - 23./3. + armWWbar[194];
   armWWbar[298]=1./2.*armWWbar[298] + armWWbar[196];
   armWWbar[299]= - 2./3.*armWWbar[11];
   armWWbar[298]=1./2.*armWWbar[298] + armWWbar[299];
   armWWbar[298]=armWWbar[38]*armWWbar[298];
   armWWbar[190]=1./2.*armWWbar[295] + 1./2.*armWWbar[190] + 
   armWWbar[298];
   armWWbar[190]=MMH*armWWbar[190];
   armWWbar[295]= - 13*armWWbar[18];
   armWWbar[298]=75*armWWbar[20] + armWWbar[295] - 17*armWWbar[19];
   armWWbar[300]=7*armWWbar[38];
   armWWbar[298]=1./2.*armWWbar[298] + armWWbar[300];
   armWWbar[298]=armWWbar[38]*armWWbar[298];
   armWWbar[190]=armWWbar[195] + 1./4.*armWWbar[298] + armWWbar[190];
   armWWbar[190]=MMt*armWWbar[190];
   armWWbar[195]=armWWbar[181] + 13./8. + armWWbar[125];
   armWWbar[195]=armWWbar[37]*armWWbar[195];
   armWWbar[298]= - 15./16.*armWWbar[70];
   armWWbar[301]= - 1./4.*armWWbar[68];
   armWWbar[302]=1./4.*armWWbar[283];
   armWWbar[303]=1./2.*armWWbar[284];
   armWWbar[195]=armWWbar[303] + armWWbar[302] + 1./2.*armWWbar[195] + 
   armWWbar[301] + armWWbar[298] - armWWbar[50] - 33./16.*armWWbar[71]
    + 5./8.*armWWbar[65] + 3./8.*armWWbar[67] + 5./4. + armWWbar[53];
   armWWbar[195]=MMt*armWWbar[195];
   armWWbar[304]= - 7*armWWbar[23] + armWWbar[198] + 67./36.*
   armWWbar[22] + armWWbar[152] + 31./2.*armWWbar[43] + armWWbar[148]
    + armWWbar[134] - 67./4.*armWWbar[74] + armWWbar[72];
   armWWbar[305]=armWWbar[136] - 19./9.*armWWbar[19];
   armWWbar[305]=armWWbar[37]*armWWbar[305];
   armWWbar[304]=1./4.*armWWbar[305] + 833./72.*armWWbar[19] + 1./2.*
   armWWbar[304] + armWWbar[293];
   armWWbar[305]=armWWbar[288] + armWWbar[181] + 25./8. + armWWbar[125]
   ;
   armWWbar[305]=armWWbar[38]*armWWbar[305];
   armWWbar[194]=armWWbar[194] + armWWbar[11];
   armWWbar[194]=armWWbar[37]*armWWbar[194];
   armWWbar[194]=armWWbar[285] + 1./4.*armWWbar[194];
   armWWbar[194]=MMH*armWWbar[194];
   armWWbar[306]=25./24. + armWWbar[160];
   armWWbar[306]=1./2.*armWWbar[306] + 5./3.*armWWbar[37];
   armWWbar[306]=armWWbar[20]*armWWbar[306];
   armWWbar[194]=armWWbar[195] + 1./2.*armWWbar[194] + 1./2.*
   armWWbar[305] + 1./2.*armWWbar[304] + armWWbar[306];
   armWWbar[194]=MMt*armWWbar[194];
   armWWbar[195]= - 5 + armWWbar[160];
   armWWbar[195]=armWWbar[18]*armWWbar[195];
   armWWbar[304]= - 1 + armWWbar[192];
   armWWbar[304]=armWWbar[20]*armWWbar[304];
   armWWbar[195]=armWWbar[304] + armWWbar[290] + 1./2.*armWWbar[195];
   armWWbar[304]=1 + armWWbar[258];
   armWWbar[304]=armWWbar[11]*armWWbar[304];
   armWWbar[305]=armWWbar[294] + armWWbar[304];
   armWWbar[305]=MMH*armWWbar[305];
   armWWbar[306]= - 1./4.*armWWbar[11];
   armWWbar[307]=armWWbar[306] - 1 + 3./8.*armWWbar[45];
   armWWbar[307]=armWWbar[38]*armWWbar[307];
   armWWbar[195]=1./4.*armWWbar[305] + 1./4.*armWWbar[195] + 
   armWWbar[307];
   armWWbar[195]=MMH*armWWbar[195];
   armWWbar[305]= - 5./3.*armWWbar[19];
   armWWbar[307]= - armWWbar[18] + armWWbar[305];
   armWWbar[307]=7./2.*armWWbar[38] + 7./4.*armWWbar[307] + 59./3.*
   armWWbar[20];
   armWWbar[307]=armWWbar[38]*armWWbar[307];
   armWWbar[194]=armWWbar[194] + 1./4.*armWWbar[307] + armWWbar[195];
   armWWbar[194]=MMt*armWWbar[194];
   armWWbar[195]=3./2.*armWWbar[19];
   armWWbar[307]= - 1 - armWWbar[71];
   armWWbar[307]=MMt*armWWbar[307];
   armWWbar[307]=1./2.*armWWbar[307] + armWWbar[38] + armWWbar[195] + 
   armWWbar[226] - armWWbar[74] + armWWbar[225];
   armWWbar[307]=MMt*armWWbar[307];
   armWWbar[308]=armWWbar[38]*armWWbar[20];
   armWWbar[307]=1./2.*armWWbar[308] + armWWbar[307];
   armWWbar[307]=armWWbar[255]*MMt*armWWbar[307];
   armWWbar[194]=armWWbar[194] + 1./4.*armWWbar[307];
   armWWbar[194]=armWWbar[255]*armWWbar[194];
   armWWbar[307]=armWWbar[18] - armWWbar[20];
   armWWbar[307]=armWWbar[38]*armWWbar[307];
   armWWbar[309]=MMH*armWWbar[38]*armWWbar[11];
   armWWbar[307]=armWWbar[307] + armWWbar[309];
   armWWbar[307]=MMH*armWWbar[307];
   armWWbar[190]=armWWbar[194] + 1./108.*armWWbar[307] + armWWbar[190];
   armWWbar[190]=armWWbar[255]*armWWbar[190];
   armWWbar[194]= - 9*armWWbar[45];
   armWWbar[307]= - 9./4.*armWWbar[49];
   armWWbar[309]=armWWbar[307] + 7./8. + armWWbar[194];
   armWWbar[309]=armWWbar[37]*armWWbar[309];
   armWWbar[283]=3./2.*armWWbar[284] + 3./4.*armWWbar[283] + 1./2.*
   armWWbar[309] - 3./4.*armWWbar[68] - 45./16.*armWWbar[70] - 3*
   armWWbar[50] - 111./16.*armWWbar[71] + 15./8.*armWWbar[65] + 9./8.*
   armWWbar[67] + 1 + armWWbar[53];
   armWWbar[283]=MMt*armWWbar[283];
   armWWbar[309]=armWWbar[148] + armWWbar[134] - 75./4.*armWWbar[74] + 
   armWWbar[72];
   armWWbar[309]= - 21./2.*armWWbar[21] + 395./12.*armWWbar[22] - 51*
   armWWbar[51] + 3*armWWbar[309] + 83./2.*armWWbar[43];
   armWWbar[194]= - 15./4.*armWWbar[37] + armWWbar[307] + 131./8. + 
   armWWbar[194];
   armWWbar[194]=armWWbar[38]*armWWbar[194];
   armWWbar[307]=1 + armWWbar[237];
   armWWbar[307]=11./4.*armWWbar[11] + 1./2.*armWWbar[307] - 
   armWWbar[10];
   armWWbar[307]=armWWbar[37]*armWWbar[307];
   armWWbar[307]=3*armWWbar[285] + armWWbar[307];
   armWWbar[307]=MMH*armWWbar[307];
   armWWbar[310]= - 3./2.*armWWbar[18];
   armWWbar[311]=armWWbar[310] + 145./3.*armWWbar[19];
   armWWbar[311]=armWWbar[37]*armWWbar[311];
   armWWbar[312]=13./3.*armWWbar[37] + 139./24. + armWWbar[138];
   armWWbar[312]=armWWbar[20]*armWWbar[312];
   armWWbar[194]=armWWbar[283] + 1./2.*armWWbar[307] + 1./2.*
   armWWbar[194] + 1./2.*armWWbar[312] + 1./8.*armWWbar[311] + 613./48.
   *armWWbar[19] + 9./2.*armWWbar[292] + 1./4.*armWWbar[309] - 8*
   armWWbar[23];
   armWWbar[194]=MMt*armWWbar[194];
   armWWbar[138]=59./27.*armWWbar[36] - 337./27. + armWWbar[138];
   armWWbar[138]=armWWbar[18]*armWWbar[138];
   armWWbar[283]=1./3.*armWWbar[19];
   armWWbar[292]=armWWbar[158] - armWWbar[19];
   armWWbar[307]=1./6.*armWWbar[37]*armWWbar[292];
   armWWbar[138]=armWWbar[307] + armWWbar[283] + 3*armWWbar[290] + 1./2.
   *armWWbar[138];
   armWWbar[309]=3*armWWbar[294] + armWWbar[196];
   armWWbar[311]=59./216.*armWWbar[36] + 31./27. + 9./8.*armWWbar[39];
   armWWbar[311]=armWWbar[11]*armWWbar[311];
   armWWbar[312]=armWWbar[235] - armWWbar[11];
   armWWbar[312]=armWWbar[37]*armWWbar[312];
   armWWbar[309]=1./4.*armWWbar[312] + 1./4.*armWWbar[309] + 
   armWWbar[311];
   armWWbar[309]=MMH*armWWbar[309];
   armWWbar[311]= - 5*armWWbar[11];
   armWWbar[237]=armWWbar[311] - 11 + armWWbar[237];
   armWWbar[237]=armWWbar[38]*armWWbar[237];
   armWWbar[313]=1./4.*armWWbar[37] - 59./216.*armWWbar[36] - 31./27.
    - 9./8.*armWWbar[39];
   armWWbar[313]=armWWbar[20]*armWWbar[313];
   armWWbar[138]=armWWbar[309] + 1./4.*armWWbar[237] + 1./4.*
   armWWbar[138] + armWWbar[313];
   armWWbar[138]=MMH*armWWbar[138];
   armWWbar[237]=armWWbar[263] - 43*armWWbar[19];
   armWWbar[237]=21./4.*armWWbar[38] + 1./8.*armWWbar[237] + 97./3.*
   armWWbar[20];
   armWWbar[237]=armWWbar[38]*armWWbar[237];
   armWWbar[138]=armWWbar[194] + 1./2.*armWWbar[237] + armWWbar[138];
   armWWbar[138]=MMt*armWWbar[138];
   armWWbar[194]= - 43./18.*armWWbar[20] + 17./9.*armWWbar[18] + 
   armWWbar[245];
   armWWbar[194]=armWWbar[38]*armWWbar[194];
   armWWbar[159]= - armWWbar[10] + armWWbar[159];
   armWWbar[159]=MMH*armWWbar[38]*armWWbar[159];
   armWWbar[159]=armWWbar[194] + 1./2.*armWWbar[159];
   armWWbar[159]=MMH*armWWbar[159];
   armWWbar[138]=1./6.*armWWbar[159] + armWWbar[138];
   armWWbar[159]=armWWbar[138] + armWWbar[190];
   armWWbar[159]=armWWbar[255]*armWWbar[159];
   armWWbar[190]= - 1 + 7./2.*armWWbar[36];
   armWWbar[194]=1./9.*armWWbar[190] + armWWbar[37];
   armWWbar[194]=armWWbar[37]*armWWbar[194];
   armWWbar[237]=armWWbar[38]*armWWbar[37];
   armWWbar[263]=armWWbar[3]*armWWbar[237];
   armWWbar[309]=3*armWWbar[263];
   armWWbar[194]=1./2.*armWWbar[194] + armWWbar[309];
   armWWbar[194]=MMt*armWWbar[194];
   armWWbar[190]=1./2.*armWWbar[190] + 4*armWWbar[37];
   armWWbar[190]=armWWbar[38]*armWWbar[190];
   armWWbar[313]=3*armWWbar[277];
   armWWbar[190]=armWWbar[194] + 1./9.*armWWbar[190] + armWWbar[313];
   armWWbar[190]=MMt*armWWbar[190];
   armWWbar[190]= - 1./18.*armWWbar[276] + armWWbar[190];
   armWWbar[190]=armWWbar[255]*MMt*armWWbar[190];
   armWWbar[194]= - 43 - 59./2.*armWWbar[36];
   armWWbar[314]=1./9.*armWWbar[194] + armWWbar[268];
   armWWbar[314]=armWWbar[37]*armWWbar[314];
   armWWbar[314]=1./2.*armWWbar[314] + 9*armWWbar[263];
   armWWbar[314]=MMt*armWWbar[314];
   armWWbar[194]=1./2.*armWWbar[194] - 8*armWWbar[37];
   armWWbar[194]=armWWbar[38]*armWWbar[194];
   armWWbar[194]=armWWbar[314] + 1./9.*armWWbar[194] + 9*armWWbar[277];
   armWWbar[194]=MMt*armWWbar[194];
   armWWbar[194]= - 43./18.*armWWbar[276] + armWWbar[194];
   armWWbar[194]=MMt*armWWbar[194];
   armWWbar[190]=armWWbar[194] + armWWbar[190];
   armWWbar[190]=armWWbar[34]*armWWbar[255]*armWWbar[190];
   armWWbar[159]=armWWbar[159] + armWWbar[190];
   armWWbar[159]=armWWbar[34]*armWWbar[159];
   armWWbar[190]=armWWbar[34]*armWWbar[194];
   armWWbar[138]=armWWbar[138] + armWWbar[190];
   armWWbar[138]=armWWbar[34]*armWWbar[138];
   armWWbar[190]=1./36.*armWWbar[10];
   armWWbar[154]=armWWbar[190] + armWWbar[156] + armWWbar[155] - 25./9.
    + armWWbar[154];
   armWWbar[154]=armWWbar[18]*armWWbar[154];
   armWWbar[155]= - 3*armWWbar[95];
   armWWbar[156]=139./6. + armWWbar[155];
   armWWbar[194]= - 5./3.*armWWbar[97];
   armWWbar[156]=1./2.*armWWbar[156] + armWWbar[194];
   armWWbar[314]=17./6.*armWWbar[16];
   armWWbar[315]=1./9.*armWWbar[10];
   armWWbar[156]=armWWbar[315] + armWWbar[314] + 17./3.*armWWbar[100]
    - 1./4.*armWWbar[96] + 1./2.*armWWbar[156] - armWWbar[91];
   armWWbar[316]=9./4.*armWWbar[44];
   armWWbar[317]=1./2.*armWWbar[48];
   armWWbar[318]=armWWbar[231] + armWWbar[317] + 1./9. + armWWbar[316];
   armWWbar[318]= - 19./36.*armWWbar[11] + 1./2.*armWWbar[318] + 
   armWWbar[199];
   armWWbar[318]=armWWbar[11]*armWWbar[318];
   armWWbar[319]=armWWbar[84] + 1./3.*armWWbar[81];
   armWWbar[319]=MMH*armWWbar[319];
   armWWbar[320]=1./4.*armWWbar[319];
   armWWbar[156]=armWWbar[320] + 1./2.*armWWbar[156] + armWWbar[318];
   armWWbar[156]=MMH*armWWbar[156];
   armWWbar[318]= - 1./2.*armWWbar[48];
   armWWbar[321]= - 1./4.*armWWbar[46];
   armWWbar[322]= - 5./2.*armWWbar[11] - 61./18.*armWWbar[10] + 
   armWWbar[321] + armWWbar[318] - 37./9. + armWWbar[167];
   armWWbar[322]=armWWbar[20]*armWWbar[322];
   armWWbar[323]= - 5*armWWbar[77] + 11./3.*armWWbar[105];
   armWWbar[323]= - 37./12.*armWWbar[23] + 5./4.*armWWbar[21] - 17./6.*
   armWWbar[22] - 37./24.*armWWbar[78] + 1./8.*armWWbar[107] - 19./12.*
   armWWbar[79] + 1./8.*armWWbar[323] + 17./3.*armWWbar[106];
   armWWbar[324]= - 13./3. + 11./8.*armWWbar[10];
   armWWbar[324]=armWWbar[19]*armWWbar[324];
   armWWbar[325]= - armWWbar[11]*armWWbar[19];
   armWWbar[154]=1./2.*armWWbar[156] + 1./4.*armWWbar[322] + 163./72.*
   armWWbar[325] + 1./3.*armWWbar[324] + 1./2.*armWWbar[323] + 
   armWWbar[154];
   armWWbar[154]=MMH*armWWbar[154];
   armWWbar[156]= - armWWbar[18] + armWWbar[19];
   armWWbar[156]=armWWbar[19]*armWWbar[156];
   armWWbar[295]=armWWbar[295] - 121./3.*armWWbar[19];
   armWWbar[295]=1./2.*armWWbar[295] + 83./3.*armWWbar[20];
   armWWbar[295]=armWWbar[20]*armWWbar[295];
   armWWbar[295]=1./6.*armWWbar[295] + 1./3.*armWWbar[162] + 1./2.*
   armWWbar[156];
   armWWbar[154]=1./2.*armWWbar[295] + armWWbar[154];
   armWWbar[154]=armWWbar[76]*MMH*armWWbar[154];
   armWWbar[295]= - 39./2.*armWWbar[71] + 5*armWWbar[65] - 1 + 3*
   armWWbar[67];
   armWWbar[286]=armWWbar[286] - 1./8. - armWWbar[45];
   armWWbar[286]=armWWbar[37]*armWWbar[286];
   armWWbar[286]=armWWbar[303] + armWWbar[302] + 3./2.*armWWbar[286] + 
   armWWbar[301] + armWWbar[298] + 1./8.*armWWbar[295] - armWWbar[50];
   armWWbar[286]=MMt*armWWbar[286];
   armWWbar[134]= - 25./2.*armWWbar[23] + armWWbar[198] + 75./4.*
   armWWbar[22] + armWWbar[152] + 13*armWWbar[43] + armWWbar[148] + 
   armWWbar[134] - 79./4.*armWWbar[74] + armWWbar[72];
   armWWbar[125]=armWWbar[288] + armWWbar[181] + 53./8. + armWWbar[125]
   ;
   armWWbar[125]=armWWbar[38]*armWWbar[125];
   armWWbar[148]= - 5./2.*armWWbar[18] + 37*armWWbar[19];
   armWWbar[148]=armWWbar[37]*armWWbar[148];
   armWWbar[152]=7./3.*armWWbar[37] + 19./8. + armWWbar[160];
   armWWbar[152]=armWWbar[20]*armWWbar[152];
   armWWbar[125]=armWWbar[125] + armWWbar[152] + 1./4.*armWWbar[148] + 
   43./8.*armWWbar[19] + 1./2.*armWWbar[134] + armWWbar[293];
   armWWbar[134]=3./4.*armWWbar[45];
   armWWbar[148]=1./3. + armWWbar[134];
   armWWbar[148]=19./24.*armWWbar[11] + 1./2.*armWWbar[148] + 
   armWWbar[196];
   armWWbar[148]=armWWbar[37]*armWWbar[148];
   armWWbar[148]=1./2.*armWWbar[285] + armWWbar[148];
   armWWbar[148]=MMH*armWWbar[148];
   armWWbar[125]=armWWbar[286] + 1./2.*armWWbar[125] + armWWbar[148];
   armWWbar[125]=MMt*armWWbar[125];
   armWWbar[148]= - 23./9. + armWWbar[160];
   armWWbar[152]=11./9.*armWWbar[36];
   armWWbar[148]=1./2.*armWWbar[148] + armWWbar[152];
   armWWbar[148]=armWWbar[18]*armWWbar[148];
   armWWbar[181]= - armWWbar[18] + armWWbar[204];
   armWWbar[181]=armWWbar[37]*armWWbar[181];
   armWWbar[148]=1./2.*armWWbar[181] + armWWbar[283] + armWWbar[290] + 
   armWWbar[148];
   armWWbar[152]=armWWbar[152] + 23./9. + armWWbar[258];
   armWWbar[152]=armWWbar[11]*armWWbar[152];
   armWWbar[152]=armWWbar[152] + armWWbar[294] + armWWbar[196];
   armWWbar[181]=1./4.*armWWbar[10];
   armWWbar[198]=armWWbar[181] - armWWbar[11];
   armWWbar[198]=armWWbar[37]*armWWbar[198];
   armWWbar[152]=1./2.*armWWbar[152] + 1./3.*armWWbar[198];
   armWWbar[152]=MMH*armWWbar[152];
   armWWbar[204]= - 11./9.*armWWbar[36] - 23./9. + armWWbar[192];
   armWWbar[204]=1./2.*armWWbar[204] + armWWbar[236];
   armWWbar[204]=armWWbar[20]*armWWbar[204];
   armWWbar[134]= - 7./6.*armWWbar[11] + armWWbar[199] - 5./3. + 
   armWWbar[134];
   armWWbar[134]=armWWbar[38]*armWWbar[134];
   armWWbar[134]=armWWbar[152] + armWWbar[134] + 1./2.*armWWbar[148] + 
   armWWbar[204];
   armWWbar[134]=MMH*armWWbar[134];
   armWWbar[148]= - 23*armWWbar[18] - 9*armWWbar[19];
   armWWbar[148]=armWWbar[300] + 1./2.*armWWbar[148] + 163./3.*
   armWWbar[20];
   armWWbar[148]=armWWbar[38]*armWWbar[148];
   armWWbar[134]=1./4.*armWWbar[148] + armWWbar[134];
   armWWbar[125]=1./2.*armWWbar[134] + armWWbar[125];
   armWWbar[125]=MMt*armWWbar[125];
   armWWbar[134]= - 7 + armWWbar[257];
   armWWbar[148]=1./3.*armWWbar[134] + armWWbar[37];
   armWWbar[148]=armWWbar[37]*armWWbar[148];
   armWWbar[148]=armWWbar[148] + 6*armWWbar[263];
   armWWbar[148]=MMt*armWWbar[148];
   armWWbar[134]=armWWbar[134] - 4*armWWbar[37];
   armWWbar[134]=armWWbar[38]*armWWbar[134];
   armWWbar[134]=armWWbar[148] + 1./3.*armWWbar[134] + 6*armWWbar[277];
   armWWbar[134]=MMt*armWWbar[134];
   armWWbar[134]= - 7./3.*armWWbar[276] + armWWbar[134];
   armWWbar[134]=armWWbar[34]*MMt*armWWbar[134];
   armWWbar[148]=11./3.*armWWbar[18] + armWWbar[19];
   armWWbar[148]=1./2.*armWWbar[148] - 7./3.*armWWbar[20];
   armWWbar[148]=armWWbar[38]*armWWbar[148];
   armWWbar[152]=armWWbar[182] + 7./3.*armWWbar[11];
   armWWbar[152]=MMH*armWWbar[38]*armWWbar[152];
   armWWbar[148]=armWWbar[148] + armWWbar[152];
   armWWbar[148]=MMH*armWWbar[148];
   armWWbar[125]=armWWbar[134] + 1./6.*armWWbar[148] + armWWbar[125];
   armWWbar[125]=armWWbar[34]*armWWbar[125];
   armWWbar[134]=armWWbar[180] + armWWbar[37];
   armWWbar[148]=armWWbar[20]*armWWbar[134];
   armWWbar[152]= - armWWbar[18] - armWWbar[19];
   armWWbar[204]=1./2.*armWWbar[37]*armWWbar[152];
   armWWbar[148]=armWWbar[148] + armWWbar[204] + 1./2.*armWWbar[179] + 
   armWWbar[19];
   armWWbar[179]=1 + armWWbar[270];
   armWWbar[236]=armWWbar[11]*armWWbar[179];
   armWWbar[257]=armWWbar[241] - armWWbar[11];
   armWWbar[263]=armWWbar[37]*armWWbar[257];
   armWWbar[270]=armWWbar[263] - armWWbar[10] + armWWbar[236];
   armWWbar[270]=MMH*armWWbar[270];
   armWWbar[277]=armWWbar[133] + 1./4. + armWWbar[10];
   armWWbar[277]=armWWbar[38]*armWWbar[277];
   armWWbar[148]=1./4.*armWWbar[270] + 1./4.*armWWbar[148] + 
   armWWbar[277];
   armWWbar[148]=MMH*armWWbar[148];
   armWWbar[185]=MMH*armWWbar[185];
   armWWbar[171]=1./3.*armWWbar[185] + armWWbar[189] + 1./6.*
   armWWbar[171];
   armWWbar[171]=MMt*armWWbar[171];
   armWWbar[185]= - 1./2.*armWWbar[18];
   armWWbar[189]=armWWbar[185] + 5./3.*armWWbar[19];
   armWWbar[189]=1./2.*armWWbar[189] + 1./3.*armWWbar[20];
   armWWbar[189]=armWWbar[38]*armWWbar[189];
   armWWbar[148]=armWWbar[171] + armWWbar[189] + 1./3.*armWWbar[148];
   armWWbar[148]=MMt*armWWbar[148];
   armWWbar[134]=armWWbar[37]*armWWbar[134];
   armWWbar[134]=1./2.*armWWbar[134] + armWWbar[309];
   armWWbar[134]=MMt*armWWbar[134];
   armWWbar[171]=armWWbar[38]*armWWbar[180];
   armWWbar[134]=armWWbar[134] + 1./2.*armWWbar[171] + armWWbar[313];
   armWWbar[134]=MMt*armWWbar[134];
   armWWbar[134]= - 1./2.*armWWbar[276] + armWWbar[134];
   armWWbar[134]=armWWbar[34]*MMt*armWWbar[134];
   armWWbar[171]=armWWbar[19] - armWWbar[20];
   armWWbar[171]=armWWbar[38]*armWWbar[171];
   armWWbar[189]= - armWWbar[10] + armWWbar[11];
   armWWbar[270]=MMH*armWWbar[38]*armWWbar[189];
   armWWbar[171]=armWWbar[171] + armWWbar[270];
   armWWbar[171]=MMH*armWWbar[171];
   armWWbar[134]=armWWbar[134] + 1./12.*armWWbar[171] + armWWbar[148];
   armWWbar[134]=armWWbar[34]*armWWbar[134];
   armWWbar[148]= - 1./3. + armWWbar[243];
   armWWbar[148]=armWWbar[19]*armWWbar[148];
   armWWbar[171]= - armWWbar[18]*armWWbar[10];
   armWWbar[148]=1./3.*armWWbar[171] + armWWbar[148];
   armWWbar[243]=1./18.*armWWbar[18] - armWWbar[19];
   armWWbar[243]=armWWbar[11]*armWWbar[243];
   armWWbar[241]=1./9. + armWWbar[241];
   armWWbar[241]=1./2.*armWWbar[241] + 1./3.*armWWbar[11];
   armWWbar[241]=armWWbar[20]*armWWbar[241];
   armWWbar[148]=armWWbar[241] + 1./6.*armWWbar[148] + armWWbar[243];
   armWWbar[241]= - armWWbar[11] - 1./4. + armWWbar[10];
   armWWbar[241]=armWWbar[11]*armWWbar[241];
   armWWbar[241]=armWWbar[181] + armWWbar[241];
   armWWbar[241]=MMH*armWWbar[241];
   armWWbar[148]=1./2.*armWWbar[148] + 1./9.*armWWbar[241];
   armWWbar[148]=MMH*armWWbar[148];
   armWWbar[241]= - 11./3.*armWWbar[18] + armWWbar[266];
   armWWbar[241]=armWWbar[19]*armWWbar[241];
   armWWbar[243]=11*armWWbar[18];
   armWWbar[270]=armWWbar[243] + 19*armWWbar[19];
   armWWbar[270]=1./4.*armWWbar[270] - armWWbar[20];
   armWWbar[270]=armWWbar[20]*armWWbar[270];
   armWWbar[241]=1./4.*armWWbar[241] + 1./3.*armWWbar[270];
   armWWbar[148]=1./6.*armWWbar[241] + armWWbar[148];
   armWWbar[148]=armWWbar[76]*MMH*armWWbar[148];
   armWWbar[134]=1./2.*armWWbar[148] + armWWbar[134];
   armWWbar[134]=armWWbar[165]*armWWbar[134];
   armWWbar[125]=armWWbar[134] + 1./2.*armWWbar[154] + armWWbar[125];
   armWWbar[125]=armWWbar[165]*armWWbar[125];
   armWWbar[134]=27./8.*armWWbar[44];
   armWWbar[148]=3./4.*armWWbar[48];
   armWWbar[154]=3./8.*armWWbar[46];
   armWWbar[241]=armWWbar[235] + armWWbar[154] + armWWbar[148] - 97./9.
    + armWWbar[134];
   armWWbar[241]=armWWbar[18]*armWWbar[241];
   armWWbar[270]= - 27./2.*armWWbar[44];
   armWWbar[277]= - 3./2.*armWWbar[46];
   armWWbar[283]=armWWbar[277] + armWWbar[232] + 19./3. + armWWbar[270]
   ;
   armWWbar[283]= - 71./9.*armWWbar[11] + 1./2.*armWWbar[283] - 55./9.*
   armWWbar[10];
   armWWbar[283]=armWWbar[20]*armWWbar[283];
   armWWbar[134]=armWWbar[154] + armWWbar[148] - 1./3. + armWWbar[134];
   armWWbar[134]= - 11./48.*armWWbar[11] + 1./4.*armWWbar[134] + 
   armWWbar[315];
   armWWbar[134]=armWWbar[11]*armWWbar[134];
   armWWbar[190]=armWWbar[190] + 17./8.*armWWbar[16] + 71./24.*
   armWWbar[100] - 3./16.*armWWbar[96] - 3./4.*armWWbar[91] - 5./8.*
   armWWbar[97] - 9./16.*armWWbar[95] + 87./32. - 1./3.*armWWbar[102];
   armWWbar[285]=3*armWWbar[84] + armWWbar[81];
   armWWbar[285]=MMH*armWWbar[285];
   armWWbar[134]=1./16.*armWWbar[285] + 1./2.*armWWbar[190] + 
   armWWbar[134];
   armWWbar[134]=MMH*armWWbar[134];
   armWWbar[190]= - 2./3.*armWWbar[108];
   armWWbar[285]= - 43./3. + 11./4.*armWWbar[10];
   armWWbar[285]=armWWbar[19]*armWWbar[285];
   armWWbar[286]= - armWWbar[18] - 101*armWWbar[19];
   armWWbar[286]=armWWbar[11]*armWWbar[286];
   armWWbar[134]=armWWbar[134] + 1./8.*armWWbar[283] + 1./72.*
   armWWbar[286] + 1./6.*armWWbar[285] + 1./4.*armWWbar[241] - 7./48.*
   armWWbar[23] + 15./16.*armWWbar[21] - 53./48.*armWWbar[22] - 55./96.
   *armWWbar[78] + 3./32.*armWWbar[107] - 31./24.*armWWbar[79] + 71./24.
   *armWWbar[106] + 11./32.*armWWbar[105] + armWWbar[190] - 15./32.*
   armWWbar[77];
   armWWbar[134]=MMH*armWWbar[134];
   armWWbar[241]=armWWbar[243] + armWWbar[19];
   armWWbar[241]=armWWbar[19]*armWWbar[241];
   armWWbar[241]=5*armWWbar[162] + 1./2.*armWWbar[241];
   armWWbar[243]= - 21*armWWbar[18] - 295./9.*armWWbar[19];
   armWWbar[243]=1./8.*armWWbar[243] + 89./9.*armWWbar[20];
   armWWbar[243]=armWWbar[20]*armWWbar[243];
   armWWbar[241]=1./12.*armWWbar[241] + armWWbar[243];
   armWWbar[134]=1./2.*armWWbar[241] + armWWbar[134];
   armWWbar[134]=armWWbar[76]*MMH*armWWbar[134];
   armWWbar[125]=armWWbar[125] + armWWbar[134] + armWWbar[138];
   armWWbar[125]=armWWbar[165]*armWWbar[125];
   armWWbar[138]=armWWbar[253] + armWWbar[252] - 5./3. + armWWbar[251];
   armWWbar[138]=armWWbar[18]*armWWbar[138];
   armWWbar[241]=armWWbar[196] + armWWbar[151] - armWWbar[48] + 19 + 
   armWWbar[221];
   armWWbar[241]=1./8.*armWWbar[241] - 11./9.*armWWbar[11];
   armWWbar[241]=armWWbar[20]*armWWbar[241];
   armWWbar[243]= - 1./6.*armWWbar[11] + armWWbar[231] + armWWbar[317]
    - 1 + armWWbar[316];
   armWWbar[243]=armWWbar[11]*armWWbar[243];
   armWWbar[252]=31./32. - armWWbar[102];
   armWWbar[243]=1./8.*armWWbar[319] + 1./4.*armWWbar[243] + 17./24.*
   armWWbar[16] + 5./8.*armWWbar[100] - 1./16.*armWWbar[96] - 1./4.*
   armWWbar[91] - 5./24.*armWWbar[97] + 1./3.*armWWbar[252] - 3./16.*
   armWWbar[95];
   armWWbar[243]=MMH*armWWbar[243];
   armWWbar[235]= - 9 + armWWbar[235];
   armWWbar[235]=armWWbar[19]*armWWbar[235];
   armWWbar[252]=armWWbar[11]*armWWbar[19];
   armWWbar[253]= - 1./2.*armWWbar[79];
   armWWbar[138]=1./2.*armWWbar[243] + 1./2.*armWWbar[241] + 5./18.*
   armWWbar[252] + 1./8.*armWWbar[235] + 1./4.*armWWbar[138] + 43./48.*
   armWWbar[23] + 5./16.*armWWbar[21] + 11./48.*armWWbar[22] - 5./96.*
   armWWbar[78] + 1./32.*armWWbar[107] + armWWbar[253] + 5./8.*
   armWWbar[106] + 11./96.*armWWbar[105] + armWWbar[190] - 5./32.*
   armWWbar[77];
   armWWbar[138]=MMH*armWWbar[138];
   armWWbar[235]=armWWbar[19]*armWWbar[18];
   armWWbar[241]=armWWbar[162] + 5./3.*armWWbar[235];
   armWWbar[243]= - 25*armWWbar[18];
   armWWbar[283]=armWWbar[243] - 41./2.*armWWbar[19];
   armWWbar[283]=1./2.*armWWbar[283] + 53*armWWbar[20];
   armWWbar[283]=armWWbar[20]*armWWbar[283];
   armWWbar[241]=1./4.*armWWbar[241] + 1./3.*armWWbar[283];
   armWWbar[138]=1./6.*armWWbar[241] + armWWbar[138];
   armWWbar[138]=armWWbar[76]*MMH*armWWbar[138];
   armWWbar[241]=1./2.*armWWbar[22];
   armWWbar[283]= - armWWbar[18] + armWWbar[241] + armWWbar[106] + 
   armWWbar[253];
   armWWbar[285]=1 + armWWbar[100];
   armWWbar[285]=MMH*armWWbar[285];
   armWWbar[286]= - 1 - armWWbar[11];
   armWWbar[286]=armWWbar[20]*armWWbar[286];
   armWWbar[252]=1./4.*armWWbar[285] + 1./4.*armWWbar[286] + 1./4.*
   armWWbar[252] + 1./2.*armWWbar[283] - armWWbar[19];
   armWWbar[252]=MMH*armWWbar[252];
   armWWbar[283]=1./2.*armWWbar[152] + armWWbar[20];
   armWWbar[285]=armWWbar[20]*armWWbar[283];
   armWWbar[252]=1./2.*armWWbar[285] + armWWbar[252];
   armWWbar[252]=armWWbar[255]*armWWbar[76]*MMH*armWWbar[252];
   armWWbar[138]=armWWbar[138] + 1./12.*armWWbar[252];
   armWWbar[138]=armWWbar[255]*armWWbar[138];
   armWWbar[252]=armWWbar[315] + armWWbar[231] + armWWbar[317] - 47./9.
    + armWWbar[316];
   armWWbar[252]=armWWbar[18]*armWWbar[252];
   armWWbar[285]=armWWbar[157] + armWWbar[48] - 7./9. + armWWbar[149];
   armWWbar[285]= - 7./18.*armWWbar[11] + 1./2.*armWWbar[285] + 
   armWWbar[315];
   armWWbar[285]=armWWbar[11]*armWWbar[285];
   armWWbar[288]=61./16. - armWWbar[102];
   armWWbar[285]=armWWbar[320] + 1./2.*armWWbar[285] + 17./12.*
   armWWbar[16] + 37./24.*armWWbar[100] - 1./8.*armWWbar[96] - 1./2.*
   armWWbar[91] - 5./12.*armWWbar[97] + 1./3.*armWWbar[288] - 3./8.*
   armWWbar[95];
   armWWbar[285]=MMH*armWWbar[285];
   armWWbar[288]= - 9*armWWbar[44];
   armWWbar[290]=131./9. + armWWbar[288];
   armWWbar[290]= - 97./18.*armWWbar[11] - 49./18.*armWWbar[10] + 
   armWWbar[151] + 1./2.*armWWbar[290] - armWWbar[48];
   armWWbar[290]=armWWbar[20]*armWWbar[290];
   armWWbar[293]= - 5 + 11./16.*armWWbar[10];
   armWWbar[293]=armWWbar[19]*armWWbar[293];
   armWWbar[294]= - 1./3.*armWWbar[18] - 13./2.*armWWbar[19];
   armWWbar[294]=armWWbar[11]*armWWbar[294];
   armWWbar[190]=1./2.*armWWbar[285] + 1./8.*armWWbar[290] + 1./24.*
   armWWbar[294] + 1./3.*armWWbar[293] + 1./4.*armWWbar[252] + 5./8.*
   armWWbar[23] + 5./8.*armWWbar[21] - 19./48.*armWWbar[22] - 3./16.*
   armWWbar[78] + 1./16.*armWWbar[107] - 43./48.*armWWbar[79] + 37./24.
   *armWWbar[106] + 11./48.*armWWbar[105] + armWWbar[190] - 5./16.*
   armWWbar[77];
   armWWbar[190]=MMH*armWWbar[190];
   armWWbar[252]=17*armWWbar[18] + armWWbar[266];
   armWWbar[252]=armWWbar[19]*armWWbar[252];
   armWWbar[243]=91*armWWbar[20] + armWWbar[243] - 29*armWWbar[19];
   armWWbar[243]=armWWbar[20]*armWWbar[243];
   armWWbar[243]=1./3.*armWWbar[243] + armWWbar[162] + 1./6.*
   armWWbar[252];
   armWWbar[190]=1./8.*armWWbar[243] + armWWbar[190];
   armWWbar[190]=armWWbar[76]*MMH*armWWbar[190];
   armWWbar[138]=armWWbar[190] + armWWbar[138];
   armWWbar[138]=armWWbar[255]*armWWbar[138];
   armWWbar[134]=armWWbar[134] + armWWbar[138];
   armWWbar[134]=armWWbar[255]*armWWbar[134];
   armWWbar[138]=armWWbar[244] + armWWbar[105] + armWWbar[253];
   armWWbar[190]=1./3.*armWWbar[138];
   armWWbar[243]= - 1./2. + armWWbar[199];
   armWWbar[133]=1./3.*armWWbar[243] + armWWbar[133];
   armWWbar[133]=armWWbar[20]*armWWbar[133];
   armWWbar[243]= - 1 - 1./9.*armWWbar[10];
   armWWbar[243]=armWWbar[18]*armWWbar[243];
   armWWbar[252]=7./2.*armWWbar[18] + armWWbar[19];
   armWWbar[252]=armWWbar[11]*armWWbar[252];
   armWWbar[133]=armWWbar[133] + 1./9.*armWWbar[252] + armWWbar[190] + 
   armWWbar[243];
   armWWbar[219]=armWWbar[196] + armWWbar[219];
   armWWbar[219]=armWWbar[11]*armWWbar[219];
   armWWbar[243]= - 1./4.*armWWbar[16] + 1./4.*armWWbar[91] + 1./2. + 
   armWWbar[101];
   armWWbar[219]=armWWbar[243] + 1./2.*armWWbar[219];
   armWWbar[219]=MMH*armWWbar[219];
   armWWbar[133]=1./2.*armWWbar[133] + 1./3.*armWWbar[219];
   armWWbar[133]=MMH*armWWbar[133];
   armWWbar[219]=1./2.*armWWbar[162] + armWWbar[235];
   armWWbar[252]= - 7./2.*armWWbar[18] - armWWbar[19];
   armWWbar[252]=1./6.*armWWbar[252] - armWWbar[20];
   armWWbar[252]=armWWbar[20]*armWWbar[252];
   armWWbar[219]=1./6.*armWWbar[219] + armWWbar[252];
   armWWbar[133]=1./3.*armWWbar[219] + armWWbar[133];
   armWWbar[219]=pow(MMH,2);
   armWWbar[133]=armWWbar[76]*armWWbar[219]*armWWbar[133];
   armWWbar[136]=armWWbar[136] + armWWbar[226] - armWWbar[73] + 
   armWWbar[225];
   armWWbar[252]= - armWWbar[70] - 1 - armWWbar[60];
   armWWbar[253]=1./2.*armWWbar[16];
   armWWbar[252]=armWWbar[253] + 1./2.*armWWbar[252] - armWWbar[68];
   armWWbar[285]=armWWbar[199] - armWWbar[11];
   armWWbar[285]=armWWbar[37]*armWWbar[285];
   armWWbar[285]=armWWbar[252] + armWWbar[285];
   armWWbar[285]=MMH*armWWbar[285];
   armWWbar[290]= - armWWbar[18] + armWWbar[265];
   armWWbar[293]=armWWbar[37]*armWWbar[290];
   armWWbar[294]=armWWbar[20]*armWWbar[37];
   armWWbar[295]=1./2.*armWWbar[294];
   armWWbar[285]=1./2.*armWWbar[285] + 1./2.*armWWbar[38] + 
   armWWbar[295] + 1./2.*armWWbar[136] + 1./3.*armWWbar[293];
   armWWbar[285]=MMH*armWWbar[285];
   armWWbar[293]= - 1 - armWWbar[70];
   armWWbar[293]=MMt*MMH*armWWbar[293];
   armWWbar[285]=armWWbar[285] + 1./4.*armWWbar[293];
   armWWbar[285]=MMt*armWWbar[285];
   armWWbar[298]= - 1./2.*armWWbar[51];
   armWWbar[300]= - 3./2.*armWWbar[11] + 1./2. + armWWbar[199];
   armWWbar[300]=armWWbar[38]*armWWbar[300];
   armWWbar[300]=armWWbar[298] + armWWbar[300];
   armWWbar[300]=MMH*armWWbar[300];
   armWWbar[290]=1./3.*armWWbar[290] + 3./4.*armWWbar[20];
   armWWbar[290]=armWWbar[38]*armWWbar[290];
   armWWbar[290]=armWWbar[290] + 1./2.*armWWbar[300];
   armWWbar[290]=MMH*armWWbar[290];
   armWWbar[285]=armWWbar[290] + armWWbar[285];
   armWWbar[285]=MMt*armWWbar[285];
   armWWbar[290]=MMt*pow(armWWbar[37],2);
   armWWbar[237]=armWWbar[237] + 1./2.*armWWbar[290];
   armWWbar[237]=MMt*armWWbar[237];
   armWWbar[237]=1./2.*armWWbar[276] + armWWbar[237];
   armWWbar[237]=armWWbar[237]*pow(MMt,2);
   armWWbar[276]=armWWbar[34]*armWWbar[237];
   armWWbar[290]=3*armWWbar[276];
   armWWbar[285]=armWWbar[285] + armWWbar[290];
   armWWbar[285]=armWWbar[34]*armWWbar[285];
   armWWbar[200]=armWWbar[11]*armWWbar[200];
   armWWbar[171]=armWWbar[171] + armWWbar[200];
   armWWbar[200]=armWWbar[20]*armWWbar[257];
   armWWbar[189]=MMH*armWWbar[11]*armWWbar[189];
   armWWbar[171]=1./2.*armWWbar[189] + 1./2.*armWWbar[171] + 
   armWWbar[200];
   armWWbar[171]=MMH*armWWbar[171];
   armWWbar[152]=armWWbar[152] + armWWbar[20];
   armWWbar[152]=armWWbar[20]*armWWbar[152];
   armWWbar[152]=armWWbar[235] + armWWbar[152];
   armWWbar[152]=1./2.*armWWbar[152] + armWWbar[171];
   armWWbar[152]=armWWbar[76]*armWWbar[219]*armWWbar[152];
   armWWbar[171]=MMH*armWWbar[263];
   armWWbar[171]=armWWbar[171] + armWWbar[204] + armWWbar[294];
   armWWbar[171]=MMt*MMH*armWWbar[171];
   armWWbar[189]=armWWbar[38]*armWWbar[283];
   armWWbar[200]=MMH*armWWbar[38]*armWWbar[257];
   armWWbar[189]=armWWbar[189] + armWWbar[200];
   armWWbar[189]=MMH*armWWbar[189];
   armWWbar[171]=armWWbar[189] + armWWbar[171];
   armWWbar[171]=MMt*armWWbar[171];
   armWWbar[171]=1./6.*armWWbar[171] + armWWbar[276];
   armWWbar[171]=armWWbar[34]*armWWbar[171];
   armWWbar[152]=1./36.*armWWbar[152] + armWWbar[171];
   armWWbar[152]=armWWbar[165]*armWWbar[152];
   armWWbar[133]=armWWbar[152] + 1./2.*armWWbar[133] + armWWbar[285];
   armWWbar[133]=armWWbar[165]*armWWbar[133];
   armWWbar[152]= - 1 - 1./24.*armWWbar[10];
   armWWbar[152]=armWWbar[18]*armWWbar[152];
   armWWbar[171]=armWWbar[147] + armWWbar[19];
   armWWbar[171]=armWWbar[11]*armWWbar[171];
   armWWbar[189]= - 1 + armWWbar[181];
   armWWbar[189]=1./2.*armWWbar[189] - armWWbar[11];
   armWWbar[189]=armWWbar[20]*armWWbar[189];
   armWWbar[152]=1./3.*armWWbar[189] + 1./24.*armWWbar[171] + 
   armWWbar[190] + armWWbar[152];
   armWWbar[171]=armWWbar[11]*armWWbar[213];
   armWWbar[171]=armWWbar[243] + 1./8.*armWWbar[171];
   armWWbar[171]=MMH*armWWbar[171];
   armWWbar[152]=1./2.*armWWbar[152] + 1./3.*armWWbar[171];
   armWWbar[152]=MMH*armWWbar[152];
   armWWbar[171]=armWWbar[162] + armWWbar[235];
   armWWbar[189]= - 7*armWWbar[18];
   armWWbar[200]=armWWbar[189] - armWWbar[19];
   armWWbar[204]= - 7*armWWbar[20];
   armWWbar[200]=1./6.*armWWbar[200] + armWWbar[204];
   armWWbar[200]=armWWbar[20]*armWWbar[200];
   armWWbar[171]=1./6.*armWWbar[171] + armWWbar[200];
   armWWbar[152]=1./8.*armWWbar[171] + armWWbar[152];
   armWWbar[152]=armWWbar[76]*armWWbar[219]*armWWbar[152];
   armWWbar[171]=armWWbar[252] + 1./2.*armWWbar[198];
   armWWbar[171]=MMH*armWWbar[171];
   armWWbar[198]= - 3*armWWbar[18];
   armWWbar[200]=armWWbar[198] - armWWbar[19];
   armWWbar[213]=armWWbar[37]*armWWbar[200];
   armWWbar[171]=armWWbar[171] + armWWbar[38] + armWWbar[295] + 
   armWWbar[136] + 1./8.*armWWbar[213];
   armWWbar[171]=MMH*armWWbar[171];
   armWWbar[213]=1./2.*armWWbar[293];
   armWWbar[171]=armWWbar[171] + armWWbar[213];
   armWWbar[171]=MMt*armWWbar[171];
   armWWbar[181]=1 + armWWbar[181];
   armWWbar[181]=1./2.*armWWbar[181] - armWWbar[11];
   armWWbar[181]=armWWbar[38]*armWWbar[181];
   armWWbar[181]=armWWbar[298] + armWWbar[181];
   armWWbar[181]=MMH*armWWbar[181];
   armWWbar[200]=1./8.*armWWbar[200] + armWWbar[20];
   armWWbar[200]=armWWbar[38]*armWWbar[200];
   armWWbar[181]=armWWbar[200] + armWWbar[181];
   armWWbar[181]=MMH*armWWbar[181];
   armWWbar[171]=armWWbar[181] + armWWbar[171];
   armWWbar[171]=MMt*armWWbar[171];
   armWWbar[181]=armWWbar[171] + armWWbar[290];
   armWWbar[181]=armWWbar[34]*armWWbar[181];
   armWWbar[133]=1./2.*armWWbar[133] + armWWbar[152] + armWWbar[181];
   armWWbar[133]=armWWbar[165]*armWWbar[133];
   armWWbar[147]=armWWbar[147] + armWWbar[245];
   armWWbar[147]=armWWbar[11]*armWWbar[147];
   armWWbar[181]= - 3 - 1./18.*armWWbar[10];
   armWWbar[181]=armWWbar[18]*armWWbar[181];
   armWWbar[200]= - 5./3.*armWWbar[11] - 1 + armWWbar[315];
   armWWbar[200]=armWWbar[20]*armWWbar[200];
   armWWbar[138]=1./2.*armWWbar[200] + 1./9.*armWWbar[147] + 
   armWWbar[138] + armWWbar[181];
   armWWbar[147]=armWWbar[196] + armWWbar[11];
   armWWbar[147]=armWWbar[11]*armWWbar[147];
   armWWbar[147]=armWWbar[243] + 1./12.*armWWbar[147];
   armWWbar[147]=MMH*armWWbar[147];
   armWWbar[138]=1./2.*armWWbar[138] + armWWbar[147];
   armWWbar[138]=MMH*armWWbar[138];
   armWWbar[147]=armWWbar[162] + 1./2.*armWWbar[235];
   armWWbar[181]=armWWbar[189] + armWWbar[265];
   armWWbar[181]=1./3.*armWWbar[181] - 19*armWWbar[20];
   armWWbar[181]=armWWbar[20]*armWWbar[181];
   armWWbar[147]=1./3.*armWWbar[147] + armWWbar[181];
   armWWbar[138]=1./6.*armWWbar[147] + armWWbar[138];
   armWWbar[138]=armWWbar[76]*armWWbar[219]*armWWbar[138];
   armWWbar[147]=pow(armWWbar[11],2);
   armWWbar[147]=armWWbar[243] + 1./24.*armWWbar[147];
   armWWbar[147]=MMH*armWWbar[147];
   armWWbar[181]=armWWbar[190] - armWWbar[18];
   armWWbar[189]= - 1./3.*armWWbar[11];
   armWWbar[190]= - 1./4. + armWWbar[189];
   armWWbar[190]=armWWbar[20]*armWWbar[190];
   armWWbar[200]=armWWbar[11]*armWWbar[18];
   armWWbar[147]=1./3.*armWWbar[147] + 1./3.*armWWbar[190] + 1./2.*
   armWWbar[181] + 1./9.*armWWbar[200];
   armWWbar[147]=MMH*armWWbar[147];
   armWWbar[190]= - armWWbar[18] - 101./8.*armWWbar[20];
   armWWbar[190]=armWWbar[20]*armWWbar[190];
   armWWbar[190]=1./8.*armWWbar[162] + armWWbar[190];
   armWWbar[147]=1./9.*armWWbar[190] + armWWbar[147];
   armWWbar[147]=armWWbar[76]*armWWbar[219]*armWWbar[147];
   armWWbar[181]=1./6.*armWWbar[286] + armWWbar[181] + 1./6.*
   armWWbar[200];
   armWWbar[190]=MMH*armWWbar[243];
   armWWbar[181]=1./2.*armWWbar[181] + 1./3.*armWWbar[190];
   armWWbar[181]=MMH*armWWbar[181];
   armWWbar[190]= - armWWbar[18] + armWWbar[272];
   armWWbar[190]=armWWbar[20]*armWWbar[190];
   armWWbar[181]=1./12.*armWWbar[190] + armWWbar[181];
   armWWbar[181]=armWWbar[255]*armWWbar[76]*armWWbar[219]*armWWbar[181]
   ;
   armWWbar[147]=armWWbar[147] + 1./2.*armWWbar[181];
   armWWbar[147]=armWWbar[255]*armWWbar[147];
   armWWbar[138]=1./2.*armWWbar[138] + armWWbar[147];
   armWWbar[138]=armWWbar[255]*armWWbar[138];
   armWWbar[138]=armWWbar[152] + 1./2.*armWWbar[138];
   armWWbar[138]=armWWbar[255]*armWWbar[138];
   armWWbar[147]=3*armWWbar[252] + armWWbar[312];
   armWWbar[147]=MMH*armWWbar[147];
   armWWbar[147]=armWWbar[147] + 3*armWWbar[38] + armWWbar[294] + 3*
   armWWbar[136] + armWWbar[307];
   armWWbar[147]=MMH*armWWbar[147];
   armWWbar[147]=armWWbar[147] + 3./2.*armWWbar[293];
   armWWbar[147]=MMt*armWWbar[147];
   armWWbar[152]=armWWbar[311] + 3 + armWWbar[199];
   armWWbar[152]=armWWbar[38]*armWWbar[152];
   armWWbar[152]=armWWbar[224] + armWWbar[152];
   armWWbar[152]=MMH*armWWbar[152];
   armWWbar[181]=1./3.*armWWbar[292] + 5*armWWbar[20];
   armWWbar[181]=armWWbar[38]*armWWbar[181];
   armWWbar[152]=armWWbar[181] + armWWbar[152];
   armWWbar[152]=MMH*armWWbar[152];
   armWWbar[147]=1./2.*armWWbar[152] + armWWbar[147];
   armWWbar[147]=MMt*armWWbar[147];
   armWWbar[152]=armWWbar[252] + armWWbar[297];
   armWWbar[152]=MMH*armWWbar[152];
   armWWbar[152]=armWWbar[152] + armWWbar[38] + 1./6.*armWWbar[294] + 
   armWWbar[136] + armWWbar[282];
   armWWbar[152]=MMH*armWWbar[152];
   armWWbar[152]=armWWbar[152] + armWWbar[213];
   armWWbar[152]=MMt*armWWbar[152];
   armWWbar[181]=armWWbar[188] + armWWbar[20];
   armWWbar[181]=armWWbar[38]*armWWbar[181];
   armWWbar[189]=1./4. + armWWbar[189];
   armWWbar[189]=armWWbar[38]*armWWbar[189];
   armWWbar[189]= - 1./4.*armWWbar[51] + armWWbar[189];
   armWWbar[189]=MMH*armWWbar[189];
   armWWbar[181]=1./3.*armWWbar[181] + armWWbar[189];
   armWWbar[181]=MMH*armWWbar[181];
   armWWbar[152]=armWWbar[181] + 1./2.*armWWbar[152];
   armWWbar[152]=MMt*armWWbar[152];
   armWWbar[181]=MMH*armWWbar[252];
   armWWbar[136]=armWWbar[181] + armWWbar[136] + armWWbar[38];
   armWWbar[136]=MMH*armWWbar[136];
   armWWbar[136]=armWWbar[136] + armWWbar[213];
   armWWbar[136]=MMt*armWWbar[136];
   armWWbar[181]=1 - armWWbar[11];
   armWWbar[181]=armWWbar[38]*armWWbar[181];
   armWWbar[181]= - armWWbar[51] + armWWbar[181];
   armWWbar[181]=MMH*armWWbar[181];
   armWWbar[181]=armWWbar[308] + armWWbar[181];
   armWWbar[181]=MMH*armWWbar[181];
   armWWbar[136]=1./2.*armWWbar[181] + armWWbar[136];
   armWWbar[136]=armWWbar[255]*MMt*armWWbar[136];
   armWWbar[136]=armWWbar[152] + 1./4.*armWWbar[136];
   armWWbar[136]=armWWbar[255]*armWWbar[136];
   armWWbar[136]=1./4.*armWWbar[147] + armWWbar[136];
   armWWbar[136]=armWWbar[255]*armWWbar[136];
   armWWbar[136]=armWWbar[171] + armWWbar[136];
   armWWbar[136]=armWWbar[255]*armWWbar[136];
   armWWbar[147]=3*armWWbar[237];
   armWWbar[152]=armWWbar[255]*armWWbar[237];
   armWWbar[152]=armWWbar[147] + armWWbar[152];
   armWWbar[152]=armWWbar[255]*armWWbar[152];
   armWWbar[147]=armWWbar[147] + 1./2.*armWWbar[152];
   armWWbar[147]=armWWbar[34]*armWWbar[255]*armWWbar[147];
   armWWbar[136]=armWWbar[136] + armWWbar[147];
   armWWbar[136]=armWWbar[34]*armWWbar[136];
   armWWbar[133]=armWWbar[133] + armWWbar[138] + armWWbar[136];
   armWWbar[133]=armWWbar[4]*armWWbar[133];
   armWWbar[125]=armWWbar[133] + armWWbar[125] + armWWbar[134] + 
   armWWbar[159];
   armWWbar[125]=armWWbar[4]*armWWbar[125];
   armWWbar[114]=armWWbar[125] + armWWbar[114] + armWWbar[120] + 
   armWWbar[130];
   armWWbar[114]=armWWbar[4]*armWWbar[114];
   armWWbar[120]=1./2.*armWWbar[13];
   armWWbar[125]= - 22*armWWbar[47];
   armWWbar[130]= - 6*armWWbar[46];
   armWWbar[133]= - 19./2.*armWWbar[12];
   armWWbar[134]= - 1 + 3*armWWbar[6];
   armWWbar[134]=2*armWWbar[1]*armWWbar[134];
   armWWbar[136]=22./3.*armWWbar[207];
   armWWbar[138]= - 3*armWWbar[49];
   armWWbar[147]= - 6*armWWbar[48];
   armWWbar[152]=6*armWWbar[10];
   armWWbar[159]= - 3*armWWbar[11];
   armWWbar[171]=armWWbar[159] + armWWbar[136] + armWWbar[134] + 
   armWWbar[133] + armWWbar[152] + armWWbar[120] + armWWbar[130] + 
   armWWbar[147] + armWWbar[138] - 37./3. + armWWbar[125];
   armWWbar[171]=armWWbar[38]*armWWbar[171];
   armWWbar[181]= - 9*armWWbar[39];
   armWWbar[189]= - 139./6.*armWWbar[36] - 41./3. + armWWbar[181];
   armWWbar[189]=armWWbar[19]*armWWbar[189];
   armWWbar[190]=1./4.*armWWbar[21];
   armWWbar[142]=armWWbar[18] + armWWbar[142];
   armWWbar[142]=1./4.*armWWbar[37]*armWWbar[142];
   armWWbar[189]=armWWbar[142] + 1./2.*armWWbar[189] + armWWbar[188] + 
   armWWbar[222] + armWWbar[190] - 11*armWWbar[22] - 25*armWWbar[51] + 
   11*armWWbar[42] + 3*armWWbar[43];
   armWWbar[200]= - 1 + armWWbar[131];
   armWWbar[200]=1./2.*armWWbar[200] - armWWbar[37];
   armWWbar[200]=armWWbar[20]*armWWbar[200];
   armWWbar[171]=armWWbar[171] + armWWbar[189] + 3*armWWbar[200];
   armWWbar[171]=armWWbar[3]*armWWbar[171];
   armWWbar[200]=1./3.*armWWbar[13];
   armWWbar[213]= - 21*armWWbar[47];
   armWWbar[219]=3*armWWbar[49];
   armWWbar[235]=armWWbar[200] + armWWbar[219] - 14999./144. + 
   armWWbar[213];
   armWWbar[235]=1./4.*armWWbar[235] + armWWbar[10];
   armWWbar[237]=14./3.*armWWbar[12];
   armWWbar[243]=17./18.*armWWbar[37];
   armWWbar[235]=armWWbar[243] + armWWbar[218] + armWWbar[217] + 
   armWWbar[216] + 1./2.*armWWbar[235] + armWWbar[237];
   armWWbar[235]=armWWbar[37]*armWWbar[235];
   armWWbar[252]=armWWbar[19]*armWWbar[39];
   armWWbar[257]=armWWbar[20]*armWWbar[39];
   armWWbar[263]=armWWbar[252] + armWWbar[257];
   armWWbar[263]=armWWbar[3]*armWWbar[263];
   armWWbar[272]= - 1 - armWWbar[39];
   armWWbar[272]=armWWbar[11]*armWWbar[272];
   armWWbar[169]=18*armWWbar[263] - 5*armWWbar[37] + 3*armWWbar[272] + 
   armWWbar[160] + armWWbar[124] + 142./3. + armWWbar[169];
   armWWbar[169]=armWWbar[3]*armWWbar[169];
   armWWbar[263]= - 3*armWWbar[68];
   armWWbar[191]=armWWbar[160] - 5./2.*armWWbar[16] + armWWbar[263] - 1.
   /2.*armWWbar[70] + armWWbar[287] - 43./8. + armWWbar[191];
   armWWbar[191]=armWWbar[112]*armWWbar[191];
   armWWbar[276]=4*armWWbar[57] - 11*armWWbar[56];
   armWWbar[282]= - 3*armWWbar[11]*armWWbar[112]*armWWbar[39];
   armWWbar[283]=armWWbar[11]*armWWbar[112];
   armWWbar[283]= - armWWbar[112] + armWWbar[283];
   armWWbar[283]=1./2.*armWWbar[37]*armWWbar[283];
   armWWbar[169]=armWWbar[169] + armWWbar[283] + armWWbar[282] + 
   armWWbar[191] - armWWbar[55] - 11./12.*armWWbar[59] + 2./3.*
   armWWbar[276] - 3./2.*armWWbar[54];
   armWWbar[169]=MMt*armWWbar[169];
   armWWbar[276]= - 139523./128. + 580*armWWbar[64];
   armWWbar[285]=13./2.*armWWbar[53];
   armWWbar[276]=1./3.*armWWbar[276] + armWWbar[285];
   armWWbar[286]=101./8.*armWWbar[61];
   armWWbar[290]=49./8.*armWWbar[67];
   armWWbar[292]=349./128.*armWWbar[66];
   armWWbar[276]= - 8851./144.*armWWbar[63] + armWWbar[292] + 
   armWWbar[290] + 1./3.*armWWbar[276] + armWWbar[286];
   armWWbar[293]=armWWbar[72] - 1./2.*armWWbar[41];
   armWWbar[294]= - 7./2.*armWWbar[23];
   armWWbar[224]=armWWbar[294] + 7./4.*armWWbar[21] + armWWbar[224] + 
   armWWbar[225] + 3*armWWbar[293] - armWWbar[73];
   armWWbar[224]=armWWbar[112]*armWWbar[224];
   armWWbar[225]= - 1./2. - armWWbar[39];
   armWWbar[225]=armWWbar[112]*armWWbar[225];
   armWWbar[141]=3*armWWbar[225] + armWWbar[141];
   armWWbar[141]=armWWbar[20]*armWWbar[141];
   armWWbar[225]= - 1./12.*armWWbar[13];
   armWWbar[293]= - 10741./36. - 151*armWWbar[36];
   armWWbar[293]=armWWbar[12]*armWWbar[293];
   armWWbar[295]=19./12.*armWWbar[33];
   armWWbar[297]=13./12.*armWWbar[30];
   armWWbar[300]=169./192.*armWWbar[69];
   armWWbar[301]=1651./384.*armWWbar[71];
   armWWbar[302]= - 1037./384.*armWWbar[50];
   armWWbar[303]=11./2.*armWWbar[47];
   armWWbar[307]= - 1./4.*armWWbar[65];
   armWWbar[308]= - 13./24.*armWWbar[70];
   armWWbar[160]=5./2. + armWWbar[160];
   armWWbar[160]=1./2.*armWWbar[18]*armWWbar[112]*armWWbar[160];
   armWWbar[163]=5*armWWbar[163];
   armWWbar[171]=armWWbar[169] + armWWbar[284] + armWWbar[171] + 
   armWWbar[163] + armWWbar[141] + armWWbar[235] + armWWbar[304] + 
   armWWbar[160] + armWWbar[275] + armWWbar[259] + armWWbar[224] + 1./
   12.*armWWbar[293] - armWWbar[10] + 139./36.*armWWbar[36] + 
   armWWbar[225] + armWWbar[192] + armWWbar[291] + 6565./432.*
   armWWbar[17] + armWWbar[233] + armWWbar[303] + armWWbar[262] + 
   armWWbar[308] + armWWbar[302] + armWWbar[287] + armWWbar[301] + 
   armWWbar[300] + armWWbar[297] + armWWbar[295] - armWWbar[62] + 1./3.
   *armWWbar[276] + armWWbar[307];
   armWWbar[171]=MMt*armWWbar[171];
   armWWbar[235]= - 17*armWWbar[36];
   armWWbar[276]= - 5*armWWbar[6];
   armWWbar[293]=armWWbar[276] + 161./6. + armWWbar[235];
   armWWbar[293]=1./9.*armWWbar[18]*armWWbar[293];
   armWWbar[304]=29*armWWbar[74];
   armWWbar[309]= - 43./2.*armWWbar[8] + armWWbar[73] - 1345./12.*
   armWWbar[72] + armWWbar[304] - 1915./96.*armWWbar[42];
   armWWbar[311]=37*armWWbar[6] - 7505./64. + armWWbar[235];
   armWWbar[311]=armWWbar[19]*armWWbar[311];
   armWWbar[312]=armWWbar[158] - 281./288.*armWWbar[19];
   armWWbar[312]=armWWbar[37]*armWWbar[312];
   armWWbar[313]= - 167*armWWbar[6] - 13525./24. + 19*armWWbar[36];
   armWWbar[313]=1./9.*armWWbar[313] + 133./2.*armWWbar[37];
   armWWbar[313]=armWWbar[20]*armWWbar[313];
   armWWbar[315]= - 59./2.*armWWbar[43];
   armWWbar[309]=1./3.*armWWbar[313] + 1./2.*armWWbar[312] + 1./9.*
   armWWbar[311] + armWWbar[293] + 2617./36.*armWWbar[23] + 
   armWWbar[139] + 789./64.*armWWbar[22] + 36127./288.*armWWbar[51] + 1.
   /3.*armWWbar[309] + armWWbar[315];
   armWWbar[311]=5./3.*armWWbar[35];
   armWWbar[312]= - 69./4.*armWWbar[110];
   armWWbar[313]=407./16.*armWWbar[111];
   armWWbar[319]=armWWbar[313] + armWWbar[311] + armWWbar[312];
   armWWbar[319]=armWWbar[38]*armWWbar[319];
   armWWbar[320]=200183./432. + armWWbar[47];
   armWWbar[320]= - 1./3.*armWWbar[13] + 1./2.*armWWbar[320] + 
   armWWbar[219];
   armWWbar[322]=1./2.*armWWbar[110];
   armWWbar[323]=1./3.*armWWbar[35] + armWWbar[322];
   armWWbar[324]=armWWbar[323] - 793./48.*armWWbar[111];
   armWWbar[324]=armWWbar[19]*armWWbar[324];
   armWWbar[326]= - armWWbar[35] - 51./2.*armWWbar[110];
   armWWbar[327]=armWWbar[326] + 107./6.*armWWbar[111];
   armWWbar[327]=armWWbar[20]*armWWbar[327];
   armWWbar[328]=6889./1152.*armWWbar[37];
   armWWbar[275]=1./2.*armWWbar[319] + 1./2.*armWWbar[327] + 
   armWWbar[328] + 1./2.*armWWbar[324] + armWWbar[275] + armWWbar[259]
    - 8959./432.*armWWbar[12] + 1./4.*armWWbar[320] - armWWbar[10];
   armWWbar[275]=armWWbar[38]*armWWbar[275];
   armWWbar[319]= - 1./2.*armWWbar[42];
   armWWbar[241]=armWWbar[241] + armWWbar[319] + armWWbar[51];
   armWWbar[320]=7./3. + armWWbar[181];
   armWWbar[267]=1./2.*armWWbar[320] + armWWbar[267];
   armWWbar[267]=armWWbar[19]*armWWbar[267];
   armWWbar[320]=2 + armWWbar[206];
   armWWbar[324]= - 3*armWWbar[46];
   armWWbar[320]= - 1./2.*armWWbar[12] + 1./3.*armWWbar[320] + 
   armWWbar[324];
   armWWbar[320]=armWWbar[38]*armWWbar[320];
   armWWbar[241]=armWWbar[320] + 7./3.*armWWbar[241] + 1./2.*
   armWWbar[267];
   armWWbar[241]=armWWbar[3]*armWWbar[241];
   armWWbar[206]=11 + armWWbar[206];
   armWWbar[206]=1./3.*armWWbar[206] - armWWbar[37];
   armWWbar[267]=armWWbar[3]*armWWbar[252];
   armWWbar[206]=1./2.*armWWbar[206] + 9*armWWbar[267];
   armWWbar[206]=armWWbar[3]*armWWbar[206];
   armWWbar[267]= - 19*armWWbar[56] - armWWbar[59];
   armWWbar[206]=1./6.*armWWbar[267] + armWWbar[206];
   armWWbar[206]=MMt*armWWbar[206];
   armWWbar[267]= - 85./8. - 7./3.*armWWbar[61];
   armWWbar[320]=1235./3. + 11*armWWbar[36];
   armWWbar[320]=armWWbar[12]*armWWbar[320];
   armWWbar[327]=31*armWWbar[12] - 1./8. + armWWbar[47];
   armWWbar[327]=1./4.*armWWbar[327] - armWWbar[37];
   armWWbar[327]=armWWbar[37]*armWWbar[327];
   armWWbar[203]=armWWbar[206] + armWWbar[241] + 1./6.*armWWbar[327] + 
   1./72.*armWWbar[320] + armWWbar[203] - 205./108.*armWWbar[17] - 7./
   12.*armWWbar[47] + 143./192.*armWWbar[50] - 5./192.*armWWbar[71] + 
   605./288.*armWWbar[69] - 1./12.*armWWbar[33] + 473./216.*
   armWWbar[63] - 271./192.*armWWbar[66] + 1./8.*armWWbar[267] + 2./3.*
   armWWbar[67];
   armWWbar[203]=MMt*armWWbar[203];
   armWWbar[206]=41./18.*armWWbar[63] - 331./144.*armWWbar[66] + 45./16.
    - armWWbar[61];
   armWWbar[241]=43./2. + 7*armWWbar[36];
   armWWbar[241]=armWWbar[12]*armWWbar[241];
   armWWbar[267]=17*armWWbar[47];
   armWWbar[320]=25./16. + armWWbar[267];
   armWWbar[320]=armWWbar[37]*armWWbar[320];
   armWWbar[206]=1./6.*armWWbar[320] + 1./18.*armWWbar[241] - 11./36.*
   armWWbar[17] + 331./288.*armWWbar[50] + 137./96.*armWWbar[71] + 313./
   144.*armWWbar[69] + 1./2.*armWWbar[206] + armWWbar[229];
   armWWbar[241]=armWWbar[54] + 1./6.*armWWbar[59];
   armWWbar[241]=MMt*armWWbar[241];
   armWWbar[206]=1./2.*armWWbar[206] + 1./3.*armWWbar[241];
   armWWbar[206]=MMt*armWWbar[206];
   armWWbar[241]= - armWWbar[19]*armWWbar[111];
   armWWbar[320]=25./16.*armWWbar[37] + 25./12.*armWWbar[241] + 25./6.*
   armWWbar[12] - 281./48. + armWWbar[267];
   armWWbar[327]= - armWWbar[20]*armWWbar[111];
   armWWbar[329]= - armWWbar[38]*armWWbar[111];
   armWWbar[320]=25./8.*armWWbar[329] + 1./2.*armWWbar[320] + 25./3.*
   armWWbar[327];
   armWWbar[320]=armWWbar[38]*armWWbar[320];
   armWWbar[329]=17*armWWbar[36];
   armWWbar[330]=5*armWWbar[6];
   armWWbar[331]=armWWbar[330] - 2761./96. + armWWbar[329];
   armWWbar[331]=armWWbar[19]*armWWbar[331];
   armWWbar[235]=armWWbar[276] - 73./12. + armWWbar[235];
   armWWbar[235]=armWWbar[20]*armWWbar[235];
   armWWbar[276]=1./9.*armWWbar[74] - 13./128.*armWWbar[42];
   armWWbar[206]=1./2.*armWWbar[206] + 1./12.*armWWbar[320] + 1./72.*
   armWWbar[235] + 25./768.*armWWbar[170] + 1./72.*armWWbar[331] - 19./
   144.*armWWbar[23] + 679./2304.*armWWbar[22] + 797./1152.*
   armWWbar[51] - 7./72.*armWWbar[8] + 5*armWWbar[276] + 1./16.*
   armWWbar[72];
   armWWbar[206]=armWWbar[255]*armWWbar[206];
   armWWbar[235]= - 5*armWWbar[36];
   armWWbar[276]= - 71./18. + armWWbar[235];
   armWWbar[320]= - 71*armWWbar[6];
   armWWbar[276]=43*armWWbar[276] + armWWbar[320];
   armWWbar[276]=armWWbar[20]*armWWbar[276];
   armWWbar[331]= - 2191./8.*armWWbar[42] - 31*armWWbar[72];
   armWWbar[331]=1./2.*armWWbar[331] - 115*armWWbar[8];
   armWWbar[332]=17*armWWbar[6] - 965./432. + 49*armWWbar[36];
   armWWbar[332]=armWWbar[19]*armWWbar[332];
   armWWbar[250]=1./18.*armWWbar[276] + 391./288.*armWWbar[250] + 1./6.
   *armWWbar[332] - 175./54.*armWWbar[23] + 8515./864.*armWWbar[22] + 1.
   /27.*armWWbar[331] + 179./16.*armWWbar[51];
   armWWbar[276]= - 11*armWWbar[35];
   armWWbar[331]=armWWbar[276] + armWWbar[322];
   armWWbar[332]=armWWbar[331] + 515./12.*armWWbar[111];
   armWWbar[332]=armWWbar[19]*armWWbar[332];
   armWWbar[332]= - 175./144.*armWWbar[37] + 1./3.*armWWbar[332] + 749./
   54.*armWWbar[12] + 2233./144. - armWWbar[47];
   armWWbar[333]=11*armWWbar[35];
   armWWbar[334]=armWWbar[333] - 485./4.*armWWbar[111];
   armWWbar[334]=armWWbar[38]*armWWbar[334];
   armWWbar[332]=1./6.*armWWbar[334] + 1./2.*armWWbar[332] + 235./9.*
   armWWbar[327];
   armWWbar[332]=armWWbar[38]*armWWbar[332];
   armWWbar[250]=1./2.*armWWbar[250] + armWWbar[332];
   armWWbar[332]=armWWbar[19] - 1./2.*armWWbar[38];
   armWWbar[332]=armWWbar[3]*armWWbar[38]*armWWbar[332];
   armWWbar[203]=armWWbar[206] + armWWbar[203] + 1./2.*armWWbar[250] + 
   17./3.*armWWbar[332];
   armWWbar[203]=armWWbar[255]*armWWbar[203];
   armWWbar[206]= - 3*armWWbar[65];
   armWWbar[250]= - 1121./324. + armWWbar[206];
   armWWbar[332]=3./2.*armWWbar[50];
   armWWbar[334]= - 13./6.*armWWbar[68];
   armWWbar[250]=5./27.*armWWbar[6] + 17./27.*armWWbar[36] - 
   armWWbar[16] + armWWbar[334] + armWWbar[332] + 1./2.*armWWbar[250]
    + armWWbar[60];
   armWWbar[273]=armWWbar[278] - 38./3. + armWWbar[273];
   armWWbar[273]=armWWbar[11]*armWWbar[273];
   armWWbar[296]=1./4.*armWWbar[296];
   armWWbar[250]=armWWbar[296] + 1./4.*armWWbar[250] + 1./27.*
   armWWbar[273];
   armWWbar[250]=MMH*armWWbar[250];
   armWWbar[273]=armWWbar[305] + 3*armWWbar[20];
   armWWbar[305]= - 15*armWWbar[38];
   armWWbar[273]=2*armWWbar[273] + armWWbar[305];
   armWWbar[273]=armWWbar[3]*armWWbar[38]*armWWbar[273];
   armWWbar[171]=armWWbar[203] + armWWbar[171] + armWWbar[250] + 
   armWWbar[273] + 1./4.*armWWbar[309] + armWWbar[275];
   armWWbar[171]=armWWbar[255]*armWWbar[171];
   armWWbar[203]=armWWbar[256] + 11./3. + armWWbar[242];
   armWWbar[273]=armWWbar[203] + 1./2.*armWWbar[280];
   armWWbar[275]=armWWbar[330] - 22./3. + armWWbar[329];
   armWWbar[275]=armWWbar[3]*armWWbar[38]*armWWbar[275];
   armWWbar[273]=1./3.*armWWbar[273] + armWWbar[275];
   armWWbar[273]=MMt*armWWbar[273];
   armWWbar[203]=armWWbar[38]*armWWbar[203];
   armWWbar[203]=1./3.*armWWbar[203] + armWWbar[273];
   armWWbar[203]=armWWbar[34]*armWWbar[255]*armWWbar[203];
   armWWbar[273]= - 5 - armWWbar[36];
   armWWbar[273]=armWWbar[12]*armWWbar[273];
   armWWbar[273]=1./3.*armWWbar[273] + 1./3.*armWWbar[36] + 4./3.*
   armWWbar[17] - 4./3.*armWWbar[63] - 1 + 4./3.*armWWbar[64];
   armWWbar[273]=MMt*armWWbar[273];
   armWWbar[254]=5./3.*armWWbar[254] + armWWbar[273];
   armWWbar[171]=1./3.*armWWbar[203] + 64./3.*armWWbar[254] + 
   armWWbar[171];
   armWWbar[171]=armWWbar[34]*armWWbar[171];
   armWWbar[125]=armWWbar[159] + armWWbar[136] + armWWbar[134] + 
   armWWbar[133] + armWWbar[152] - 47./2.*armWWbar[13] + armWWbar[130]
    + armWWbar[147] + armWWbar[138] - 109./3. + armWWbar[125];
   armWWbar[125]=armWWbar[38]*armWWbar[125];
   armWWbar[130]= - 73./3. + armWWbar[181];
   armWWbar[130]=armWWbar[178] + 1./2.*armWWbar[130] - 32./3.*
   armWWbar[36];
   armWWbar[130]=armWWbar[20]*armWWbar[130];
   armWWbar[125]=armWWbar[125] + armWWbar[189] + armWWbar[130];
   armWWbar[125]=armWWbar[3]*armWWbar[125];
   armWWbar[130]= - 47./3.*armWWbar[13] + armWWbar[219] - 17303./144.
    + armWWbar[213];
   armWWbar[130]=1./4.*armWWbar[130] + armWWbar[10];
   armWWbar[130]=armWWbar[243] + armWWbar[218] + armWWbar[217] + 
   armWWbar[216] + 1./2.*armWWbar[130] + armWWbar[237];
   armWWbar[130]=armWWbar[37]*armWWbar[130];
   armWWbar[133]=armWWbar[285] - 194825./1152. - 68*armWWbar[64];
   armWWbar[133]=581./16.*armWWbar[63] + armWWbar[292] + armWWbar[290]
    + 1./3.*armWWbar[133] + armWWbar[286];
   armWWbar[134]=32*armWWbar[36];
   armWWbar[136]=43 + armWWbar[134];
   armWWbar[140]=32*armWWbar[7]*armWWbar[140];
   armWWbar[189]= - 11*armWWbar[6];
   armWWbar[136]=armWWbar[189] + 1./3.*armWWbar[136] + armWWbar[140];
   armWWbar[136]=armWWbar[40]*armWWbar[136];
   armWWbar[134]=59 + armWWbar[134];
   armWWbar[134]=1./3.*armWWbar[134] + armWWbar[140];
   armWWbar[134]=1./27.*armWWbar[134] - armWWbar[6];
   armWWbar[134]=armWWbar[1]*armWWbar[134];
   armWWbar[140]=47./12.*armWWbar[13];
   armWWbar[203]=107./4. + armWWbar[329];
   armWWbar[203]=armWWbar[12]*armWWbar[203];
   armWWbar[213]= - 32./9.*armWWbar[36] - 23./9. + armWWbar[258];
   armWWbar[213]=armWWbar[11]*armWWbar[213];
   armWWbar[125]=armWWbar[169] + armWWbar[284] + armWWbar[125] + 
   armWWbar[163] + armWWbar[141] + armWWbar[130] + armWWbar[213] + 
   armWWbar[160] + 1./9.*armWWbar[136] + armWWbar[134] + armWWbar[224]
    + 25./12.*armWWbar[203] - armWWbar[10] + 2723./324.*armWWbar[36] + 
   armWWbar[140] + armWWbar[192] + armWWbar[291] - 835./48.*
   armWWbar[17] + armWWbar[233] + armWWbar[303] + armWWbar[262] + 
   armWWbar[308] + armWWbar[302] + armWWbar[287] + armWWbar[301] + 
   armWWbar[300] + armWWbar[297] + armWWbar[295] - armWWbar[62] + 1./3.
   *armWWbar[133] + armWWbar[307];
   armWWbar[125]=MMt*armWWbar[125];
   armWWbar[130]= - 97./6.*armWWbar[8] + armWWbar[73] - 833./12.*
   armWWbar[72] + armWWbar[304] + 3205./96.*armWWbar[42];
   armWWbar[133]=53*armWWbar[6] - 42739./192. + 47*armWWbar[36];
   armWWbar[133]=armWWbar[19]*armWWbar[133];
   armWWbar[134]=armWWbar[158] + 2381./96.*armWWbar[19];
   armWWbar[134]=armWWbar[37]*armWWbar[134];
   armWWbar[136]=armWWbar[320] - 8533./24. + 403*armWWbar[36];
   armWWbar[136]=1./9.*armWWbar[136] + 229./2.*armWWbar[37];
   armWWbar[136]=armWWbar[20]*armWWbar[136];
   armWWbar[130]=1./3.*armWWbar[136] + 1./2.*armWWbar[134] + 1./9.*
   armWWbar[133] + armWWbar[293] + 403./12.*armWWbar[23] + 
   armWWbar[139] + 1087./192.*armWWbar[22] + 17695./288.*armWWbar[51]
    + 1./3.*armWWbar[130] + armWWbar[315];
   armWWbar[133]=997./1296. + armWWbar[47];
   armWWbar[133]=47./3.*armWWbar[13] + 1./2.*armWWbar[133] + 
   armWWbar[219];
   armWWbar[134]= - 32*armWWbar[7];
   armWWbar[136]=armWWbar[189] + 43./3. + armWWbar[134];
   armWWbar[136]=armWWbar[40]*armWWbar[136];
   armWWbar[169]=armWWbar[313] + 7*armWWbar[35] + armWWbar[312];
   armWWbar[169]=armWWbar[38]*armWWbar[169];
   armWWbar[134]=59./3. + armWWbar[134];
   armWWbar[134]=1./27.*armWWbar[134] - armWWbar[6];
   armWWbar[134]=armWWbar[1]*armWWbar[134];
   armWWbar[189]=armWWbar[323] + 77./16.*armWWbar[111];
   armWWbar[189]=armWWbar[19]*armWWbar[189];
   armWWbar[203]=armWWbar[326] + 121./2.*armWWbar[111];
   armWWbar[203]=armWWbar[20]*armWWbar[203];
   armWWbar[133]=1./2.*armWWbar[169] + 1./2.*armWWbar[203] + 
   armWWbar[328] - 32./9.*armWWbar[11] + 1./2.*armWWbar[189] + 1./9.*
   armWWbar[136] + armWWbar[134] + 2873./48.*armWWbar[12] + 1./4.*
   armWWbar[133] - armWWbar[10];
   armWWbar[133]=armWWbar[38]*armWWbar[133];
   armWWbar[134]=13./2.*armWWbar[36];
   armWWbar[136]=43 + armWWbar[134];
   armWWbar[136]=1./3.*armWWbar[136] + armWWbar[256];
   armWWbar[169]=5./6.*armWWbar[6] - 203./9. - 37./2.*armWWbar[36];
   armWWbar[169]=armWWbar[37]*armWWbar[169];
   armWWbar[136]=armWWbar[275] + 1./3.*armWWbar[136] + 1./2.*
   armWWbar[169];
   armWWbar[136]=MMt*armWWbar[136];
   armWWbar[169]= - 32*armWWbar[37] + armWWbar[256] + 43./3. + 
   armWWbar[242];
   armWWbar[169]=armWWbar[38]*armWWbar[169];
   armWWbar[136]=1./3.*armWWbar[169] + armWWbar[136];
   armWWbar[136]=armWWbar[34]*armWWbar[136];
   armWWbar[169]=armWWbar[266] + armWWbar[204];
   armWWbar[169]=2./3.*armWWbar[169] + armWWbar[305];
   armWWbar[169]=armWWbar[3]*armWWbar[38]*armWWbar[169];
   armWWbar[125]=1./3.*armWWbar[136] + armWWbar[125] + armWWbar[250] + 
   armWWbar[169] + 1./4.*armWWbar[130] + armWWbar[133];
   armWWbar[125]=armWWbar[34]*armWWbar[125];
   armWWbar[130]= - 1369./27. + armWWbar[234];
   armWWbar[133]=7./2.*armWWbar[46];
   armWWbar[130]= - 47./18.*armWWbar[13] + armWWbar[133] + 1./4.*
   armWWbar[130] + armWWbar[172];
   armWWbar[136]=323./144.*armWWbar[12];
   armWWbar[130]=armWWbar[306] + armWWbar[136] + 1./2.*armWWbar[130] + 
   armWWbar[199];
   armWWbar[130]=armWWbar[11]*armWWbar[130];
   armWWbar[169]=1./3. - 5./2.*armWWbar[10];
   armWWbar[169]=1./3.*armWWbar[169] + armWWbar[228];
   armWWbar[169]=armWWbar[12]*armWWbar[169];
   armWWbar[189]= - 10*armWWbar[102];
   armWWbar[203]= - 811./576. + armWWbar[189];
   armWWbar[204]= - 9./8.*armWWbar[95];
   armWWbar[213]=1./2.*armWWbar[94];
   armWWbar[216]= - 2./3.*armWWbar[99];
   armWWbar[217]= - 11./6.*armWWbar[97];
   armWWbar[218]=9./16.*armWWbar[93];
   armWWbar[228]=9./16.*armWWbar[92];
   armWWbar[237]= - 25./48.*armWWbar[98];
   armWWbar[242]= - 89./16.*armWWbar[103];
   armWWbar[243]= - 11./12.*armWWbar[96];
   armWWbar[250]=197./24.*armWWbar[100];
   armWWbar[254]=1./8.*armWWbar[17];
   armWWbar[256]=215./48.*armWWbar[16];
   armWWbar[266]= - 5./8.*armWWbar[46];
   armWWbar[273]=7./6.*armWWbar[10];
   armWWbar[275]=4./3.*armWWbar[89] - armWWbar[86];
   armWWbar[275]=23./24.*armWWbar[82] - 17./24.*armWWbar[85] - 17./48.*
   armWWbar[109] - 1./6.*armWWbar[81] + 2*armWWbar[275] - 3./4.*
   armWWbar[84];
   armWWbar[275]=MMH*armWWbar[275];
   armWWbar[280]=23./9.*armWWbar[101];
   armWWbar[285]=1./6.*armWWbar[91];
   armWWbar[286]= - 9./8.*armWWbar[44];
   armWWbar[290]= - 3./4.*armWWbar[48];
   armWWbar[130]=armWWbar[275] + armWWbar[130] + armWWbar[169] + 
   armWWbar[273] - 47./72.*armWWbar[13] + armWWbar[266] + armWWbar[256]
    + armWWbar[290] + armWWbar[254] + armWWbar[286] + armWWbar[250] + 
   armWWbar[243] + armWWbar[242] + armWWbar[285] + armWWbar[237] + 
   armWWbar[228] + armWWbar[218] + armWWbar[217] + armWWbar[216] + 
   armWWbar[280] + armWWbar[213] + 1./3.*armWWbar[203] + armWWbar[204];
   armWWbar[130]=MMH*armWWbar[130];
   armWWbar[203]= - 45*armWWbar[44];
   armWWbar[292]= - 20033./18. + armWWbar[203];
   armWWbar[293]= - 7./4.*armWWbar[46];
   armWWbar[295]= - 49./6.*armWWbar[10];
   armWWbar[297]= - armWWbar[18]*armWWbar[112];
   armWWbar[300]=2*armWWbar[297];
   armWWbar[301]=3*armWWbar[113] - armWWbar[112];
   armWWbar[301]=5./4.*armWWbar[19]*armWWbar[301];
   armWWbar[302]= - armWWbar[20]*armWWbar[112];
   armWWbar[303]=17./4.*armWWbar[302];
   armWWbar[304]= - 2*armWWbar[48];
   armWWbar[292]=armWWbar[303] + 13./6.*armWWbar[11] + armWWbar[301] + 
   armWWbar[300] - 10529./72.*armWWbar[12] + armWWbar[295] - 3523./72.*
   armWWbar[13] + armWWbar[293] + 1./8.*armWWbar[292] + armWWbar[304];
   armWWbar[292]=armWWbar[20]*armWWbar[292];
   armWWbar[305]=119./24.*armWWbar[12];
   armWWbar[140]=armWWbar[305] - armWWbar[10] + armWWbar[140] + 
   armWWbar[173] + armWWbar[48] - 241./9. + armWWbar[316];
   armWWbar[309]=2*armWWbar[122];
   armWWbar[140]=1./2.*armWWbar[140] + armWWbar[309];
   armWWbar[140]=armWWbar[18]*armWWbar[140];
   armWWbar[312]= - 167./6.*armWWbar[13] + armWWbar[277] + 
   armWWbar[232] + 3325./27. + armWWbar[270];
   armWWbar[313]=4*armWWbar[10];
   armWWbar[315]= - armWWbar[18]*armWWbar[113];
   armWWbar[320]=7./6.*armWWbar[315];
   armWWbar[323]=73./3.*armWWbar[113] - 5*armWWbar[112];
   armWWbar[323]=1./8.*armWWbar[19]*armWWbar[323];
   armWWbar[312]=armWWbar[323] + armWWbar[320] - 4435./48.*armWWbar[12]
    + 1./4.*armWWbar[312] + armWWbar[313];
   armWWbar[312]=armWWbar[19]*armWWbar[312];
   armWWbar[198]=armWWbar[198] + 13./2.*armWWbar[19];
   armWWbar[198]=armWWbar[19]*armWWbar[198];
   armWWbar[162]= - 2*armWWbar[162] + 1./4.*armWWbar[198];
   armWWbar[198]= - 41./3.*armWWbar[18] + 67./2.*armWWbar[19];
   armWWbar[198]=1./4.*armWWbar[198];
   armWWbar[328]=armWWbar[198] + 21*armWWbar[20];
   armWWbar[328]=armWWbar[20]*armWWbar[328];
   armWWbar[328]=armWWbar[162] + armWWbar[328];
   armWWbar[328]=armWWbar[3]*armWWbar[328];
   armWWbar[335]= - 9./8.*armWWbar[77];
   armWWbar[336]=61./24.*armWWbar[105];
   armWWbar[337]=605./48.*armWWbar[106];
   armWWbar[338]= - 47./12.*armWWbar[79];
   armWWbar[339]= - 29./6.*armWWbar[78];
   armWWbar[340]=269./144.*armWWbar[21];
   armWWbar[341]= - armWWbar[18] + 55./2.*armWWbar[19];
   armWWbar[341]=armWWbar[11]*armWWbar[341];
   armWWbar[130]=armWWbar[130] + armWWbar[328] + armWWbar[292] + 1./6.*
   armWWbar[341] + armWWbar[312] + armWWbar[140] - 6607./144.*
   armWWbar[23] + armWWbar[340] + 599./16.*armWWbar[22] + armWWbar[339]
    - 179./24.*armWWbar[107] + armWWbar[338] - 149./24.*armWWbar[80] + 
   armWWbar[337] + armWWbar[336] + 922./9.*armWWbar[108] + 
   armWWbar[335];
   armWWbar[130]=armWWbar[76]*armWWbar[130];
   armWWbar[140]=1 - 3*armWWbar[6];
   armWWbar[292]=armWWbar[40]*armWWbar[140];
   armWWbar[117]=armWWbar[159] + 2*armWWbar[292] + 2*armWWbar[259] - 33
   *armWWbar[12] + armWWbar[152] + 165./2.*armWWbar[13] + armWWbar[324]
    + armWWbar[147] + armWWbar[138] + 185./3. + armWWbar[117];
   armWWbar[117]=armWWbar[38]*armWWbar[117];
   armWWbar[138]=1./2.*armWWbar[42];
   armWWbar[147]= - 2*armWWbar[51];
   armWWbar[159]= - 1./2.*armWWbar[22] + armWWbar[147] + armWWbar[138]
    + armWWbar[43];
   armWWbar[181]=13 + armWWbar[181];
   armWWbar[292]=armWWbar[178] + 1./2.*armWWbar[181] + 8*armWWbar[36];
   armWWbar[292]=armWWbar[20]*armWWbar[292];
   armWWbar[181]=armWWbar[181] + 13*armWWbar[36];
   armWWbar[181]=armWWbar[19]*armWWbar[181];
   armWWbar[117]=armWWbar[117] + armWWbar[292] + armWWbar[142] + 1./4.*
   armWWbar[181] + armWWbar[188] + armWWbar[222] + 3*armWWbar[159] + 
   armWWbar[190];
   armWWbar[117]=armWWbar[3]*armWWbar[117];
   armWWbar[142]=armWWbar[252] + 2*armWWbar[257];
   armWWbar[142]=armWWbar[3]*armWWbar[142];
   armWWbar[159]= - armWWbar[49] + 9 - armWWbar[47];
   armWWbar[135]=3*armWWbar[142] + armWWbar[135] + armWWbar[272] + 1./2.
   *armWWbar[159] + armWWbar[39];
   armWWbar[135]=armWWbar[3]*armWWbar[135];
   armWWbar[142]= - armWWbar[56] - armWWbar[54];
   armWWbar[159]=armWWbar[142] - 1./2.*armWWbar[59];
   armWWbar[135]=3*armWWbar[135] + armWWbar[283] + armWWbar[282] + 
   armWWbar[191] + 3./2.*armWWbar[159] - armWWbar[55];
   armWWbar[135]=MMt*armWWbar[135];
   armWWbar[159]=55*armWWbar[13] + 5551./144. + armWWbar[219];
   armWWbar[123]= - armWWbar[37] + armWWbar[123] + armWWbar[223] + 1./3.
   *armWWbar[259] + 15./4.*armWWbar[12] + 1./4.*armWWbar[159] + 
   armWWbar[10];
   armWWbar[123]=armWWbar[37]*armWWbar[123];
   armWWbar[159]= - 11 - 8*armWWbar[36];
   armWWbar[181]=armWWbar[7]*armWWbar[143];
   armWWbar[159]=1./3.*armWWbar[159] + 8*armWWbar[181];
   armWWbar[159]=1./3.*armWWbar[159] + armWWbar[6];
   armWWbar[181]=armWWbar[1]*armWWbar[159];
   armWWbar[159]=armWWbar[40]*armWWbar[159];
   armWWbar[188]=8./3.*armWWbar[36];
   armWWbar[190]=armWWbar[188] + 11./3. + armWWbar[258];
   armWWbar[190]=armWWbar[11]*armWWbar[190];
   armWWbar[191]= - 91./2. - 55*armWWbar[36];
   armWWbar[191]=armWWbar[12]*armWWbar[191];
   armWWbar[117]=armWWbar[135] + armWWbar[284] + armWWbar[117] + 
   armWWbar[163] + armWWbar[141] + 1./2.*armWWbar[123] + armWWbar[190]
    + armWWbar[160] + armWWbar[159] + 1./3.*armWWbar[181] + 
   armWWbar[224] + 3./8.*armWWbar[191] - armWWbar[10] - 1025./108.*
   armWWbar[36] - 55./4.*armWWbar[13] + armWWbar[192] + armWWbar[291]
    - 93./16.*armWWbar[17] + armWWbar[233] + 3./4.*armWWbar[47] + 
   armWWbar[262] + armWWbar[308] - 185./128.*armWWbar[50] + 
   armWWbar[287] + 1663./128.*armWWbar[71] + 25./64.*armWWbar[69] + 7./
   4.*armWWbar[30] + armWWbar[33] - armWWbar[62] + armWWbar[307] - 31./
   16.*armWWbar[63] + 41./128.*armWWbar[66] + 11./8.*armWWbar[67] - 
   145765./3456. + 6*armWWbar[61];
   armWWbar[117]=MMt*armWWbar[117];
   armWWbar[123]=3*armWWbar[47];
   armWWbar[135]= - 39583./432. + armWWbar[123];
   armWWbar[135]= - 55*armWWbar[13] + 1./2.*armWWbar[135] + 
   armWWbar[219];
   armWWbar[141]= - 11./3. + 8*armWWbar[7];
   armWWbar[141]=1./3.*armWWbar[141] + armWWbar[6];
   armWWbar[159]=armWWbar[1]*armWWbar[141];
   armWWbar[141]=armWWbar[40]*armWWbar[141];
   armWWbar[160]= - 9./8.*armWWbar[111] - armWWbar[35] + 3./2.*
   armWWbar[110];
   armWWbar[160]=armWWbar[19]*armWWbar[160];
   armWWbar[163]=armWWbar[326] - 33./2.*armWWbar[111];
   armWWbar[163]=armWWbar[20]*armWWbar[163];
   armWWbar[181]= - 123./8.*armWWbar[111] + armWWbar[35] - 69./2.*
   armWWbar[110];
   armWWbar[181]=armWWbar[38]*armWWbar[181];
   armWWbar[135]=1./4.*armWWbar[181] + 1./2.*armWWbar[163] + 271./128.*
   armWWbar[37] + 8./3.*armWWbar[11] + 1./4.*armWWbar[160] + 
   armWWbar[141] + 1./3.*armWWbar[159] - 75./16.*armWWbar[12] + 1./4.*
   armWWbar[135] - armWWbar[10];
   armWWbar[135]=armWWbar[38]*armWWbar[135];
   armWWbar[141]=armWWbar[144] + 3./2. + 5./3.*armWWbar[36];
   armWWbar[141]=armWWbar[18]*armWWbar[141];
   armWWbar[144]=13*armWWbar[6] - 7831./96. - 7*armWWbar[36];
   armWWbar[144]=armWWbar[19]*armWWbar[144];
   armWWbar[158]=armWWbar[158] - 3341./96.*armWWbar[19];
   armWWbar[158]=armWWbar[37]*armWWbar[158];
   armWWbar[159]= - 211*armWWbar[6] - 2851./4. - 263*armWWbar[36];
   armWWbar[159]=1./9.*armWWbar[159] - 17*armWWbar[37];
   armWWbar[159]=armWWbar[20]*armWWbar[159];
   armWWbar[139]=1./2.*armWWbar[159] + 1./2.*armWWbar[158] + 1./6.*
   armWWbar[144] + armWWbar[141] + 155./4.*armWWbar[23] + armWWbar[139]
    - 1373./64.*armWWbar[22] + 3777./32.*armWWbar[51] - 19./2.*
   armWWbar[43] - 7./2.*armWWbar[8] + 1./3.*armWWbar[73] - 159./4.*
   armWWbar[72] + 23*armWWbar[74] - 41./32.*armWWbar[42];
   armWWbar[141]= - 89./36. + armWWbar[206];
   armWWbar[141]= - 1./9.*armWWbar[6] - 5./9.*armWWbar[36] - 
   armWWbar[16] + armWWbar[334] + armWWbar[332] + 1./2.*armWWbar[141]
    + armWWbar[60];
   armWWbar[144]=armWWbar[126] - 2 + armWWbar[176];
   armWWbar[144]=armWWbar[11]*armWWbar[144];
   armWWbar[141]=armWWbar[296] + 1./4.*armWWbar[141] + 1./9.*
   armWWbar[144];
   armWWbar[141]=MMH*armWWbar[141];
   armWWbar[144]= - 11 + armWWbar[271];
   armWWbar[144]=1./3.*armWWbar[144] + armWWbar[115];
   armWWbar[158]= - 1./6.*armWWbar[6] + 17./3. + 9./2.*armWWbar[36];
   armWWbar[158]=armWWbar[37]*armWWbar[158];
   armWWbar[159]= - armWWbar[6] + 2 + armWWbar[235];
   armWWbar[159]=armWWbar[3]*armWWbar[38]*armWWbar[159];
   armWWbar[144]=armWWbar[159] + 1./3.*armWWbar[144] + 1./2.*
   armWWbar[158];
   armWWbar[144]=MMt*armWWbar[144];
   armWWbar[158]=5./2.*armWWbar[36];
   armWWbar[159]=8*armWWbar[37] + armWWbar[115] - 11./3. + 
   armWWbar[158];
   armWWbar[159]=armWWbar[38]*armWWbar[159];
   armWWbar[144]=1./3.*armWWbar[159] + armWWbar[144];
   armWWbar[144]=armWWbar[34]*armWWbar[144];
   armWWbar[159]=armWWbar[19] + 2*armWWbar[20];
   armWWbar[159]=7*armWWbar[159] - 3./2.*armWWbar[38];
   armWWbar[159]=armWWbar[3]*armWWbar[38]*armWWbar[159];
   armWWbar[117]=armWWbar[144] + armWWbar[117] + armWWbar[141] + 
   armWWbar[159] + 1./4.*armWWbar[139] + armWWbar[135];
   armWWbar[117]=armWWbar[34]*armWWbar[117];
   armWWbar[135]=1931./27. + armWWbar[234];
   armWWbar[135]=55./6.*armWWbar[13] + armWWbar[173] + 1./4.*
   armWWbar[135] + armWWbar[172];
   armWWbar[135]=armWWbar[306] + 95./48.*armWWbar[12] + 1./2.*
   armWWbar[135] + armWWbar[199];
   armWWbar[135]=armWWbar[11]*armWWbar[135];
   armWWbar[139]= - 11./24.*armWWbar[109];
   armWWbar[141]= - 11./12.*armWWbar[85];
   armWWbar[144]= - 11./12.*armWWbar[82] + armWWbar[141] + 
   armWWbar[139] - 3./2.*armWWbar[84] - 1./3.*armWWbar[81];
   armWWbar[144]=MMH*armWWbar[144];
   armWWbar[159]=10895./216. - 9*armWWbar[95];
   armWWbar[160]=5./2. + 11*armWWbar[10];
   armWWbar[160]=1./2.*armWWbar[160] + armWWbar[12];
   armWWbar[160]=armWWbar[12]*armWWbar[160];
   armWWbar[135]=1./2.*armWWbar[144] + armWWbar[135] + 1./6.*
   armWWbar[160] + 1./24.*armWWbar[10] + 55./24.*armWWbar[13] - 3./8.*
   armWWbar[46] + 137./48.*armWWbar[16] + armWWbar[290] - 15./8.*
   armWWbar[17] + armWWbar[286] + 95./8.*armWWbar[100] - 5./6.*
   armWWbar[96] + armWWbar[285] + 29./48.*armWWbar[98] + 29./48.*
   armWWbar[92] + 29./48.*armWWbar[93] + armWWbar[194] + 1./3.*
   armWWbar[99] + armWWbar[280] + 1./8.*armWWbar[159] + 1./3.*
   armWWbar[94];
   armWWbar[135]=MMH*armWWbar[135];
   armWWbar[144]=15377./27. + armWWbar[203];
   armWWbar[159]=7./2.*armWWbar[113] - armWWbar[112];
   armWWbar[159]=armWWbar[19]*armWWbar[159];
   armWWbar[144]=1./3.*armWWbar[302] - 101./12.*armWWbar[11] + 1./2.*
   armWWbar[159] + 23./12.*armWWbar[297] + 1011./16.*armWWbar[12] - 53./
   12.*armWWbar[10] + 2893./48.*armWWbar[13] - armWWbar[46] + 1./8.*
   armWWbar[144] + armWWbar[304];
   armWWbar[144]=armWWbar[20]*armWWbar[144];
   armWWbar[159]=125./16.*armWWbar[12] + armWWbar[182] - 55./8.*
   armWWbar[13] + armWWbar[231] + armWWbar[317] - 245./9. + 
   armWWbar[251];
   armWWbar[159]=armWWbar[18]*armWWbar[159];
   armWWbar[160]=77./2.*armWWbar[13] + armWWbar[277] + armWWbar[232] - 
   66119./27. + armWWbar[270];
   armWWbar[163]=9*armWWbar[10];
   armWWbar[160]=551./3.*armWWbar[12] + 1./2.*armWWbar[160] + 
   armWWbar[163];
   armWWbar[176]=19./2.*armWWbar[113] - armWWbar[112];
   armWWbar[176]=armWWbar[19]*armWWbar[176];
   armWWbar[160]=1./2.*armWWbar[176] + 1./2.*armWWbar[160] + 3*
   armWWbar[315];
   armWWbar[160]=armWWbar[19]*armWWbar[160];
   armWWbar[176]=armWWbar[310] + 5*armWWbar[19];
   armWWbar[176]=armWWbar[19]*armWWbar[176];
   armWWbar[132]= - 43./4.*armWWbar[20] - 3./4.*armWWbar[18] + 
   armWWbar[132];
   armWWbar[132]=armWWbar[20]*armWWbar[132];
   armWWbar[132]=1./4.*armWWbar[176] + armWWbar[132];
   armWWbar[132]=armWWbar[3]*armWWbar[132];
   armWWbar[176]=489./2.*armWWbar[80] + 197./2.*armWWbar[106] - 9*
   armWWbar[77] + 61./3.*armWWbar[105];
   armWWbar[176]= - 2397./8.*armWWbar[23] + 863./72.*armWWbar[21] - 
   1571./4.*armWWbar[22] - 11./2.*armWWbar[78] + 159*armWWbar[107] + 1./
   4.*armWWbar[176] - 17./3.*armWWbar[79];
   armWWbar[181]=armWWbar[185] - 26*armWWbar[19];
   armWWbar[181]=armWWbar[11]*armWWbar[181];
   armWWbar[132]=armWWbar[135] + armWWbar[132] + armWWbar[144] + 1./3.*
   armWWbar[181] + 1./2.*armWWbar[160] + 1./2.*armWWbar[176] + 
   armWWbar[159];
   armWWbar[132]=armWWbar[76]*armWWbar[132];
   armWWbar[135]= - 1 - 33*armWWbar[13];
   armWWbar[135]=armWWbar[164] + 1./2.*armWWbar[135] + 2*armWWbar[10];
   armWWbar[144]= - armWWbar[7] + armWWbar[6];
   armWWbar[159]=armWWbar[1]*armWWbar[144];
   armWWbar[144]=armWWbar[40]*armWWbar[144];
   armWWbar[135]= - 6*armWWbar[11] + 6*armWWbar[144] + 3*armWWbar[135]
    + 2*armWWbar[159];
   armWWbar[135]=armWWbar[38]*armWWbar[135];
   armWWbar[160]=armWWbar[19]*armWWbar[36];
   armWWbar[160]=armWWbar[160] + armWWbar[170];
   armWWbar[170]=armWWbar[36] - armWWbar[37];
   armWWbar[176]=armWWbar[20]*armWWbar[170];
   armWWbar[160]=1./2.*armWWbar[160] + armWWbar[176];
   armWWbar[135]=3./2.*armWWbar[160] + armWWbar[135];
   armWWbar[135]=armWWbar[3]*armWWbar[135];
   armWWbar[160]=33*armWWbar[13] - 59./9.*armWWbar[36];
   armWWbar[176]=armWWbar[12]*armWWbar[180];
   armWWbar[179]=armWWbar[7]*armWWbar[179];
   armWWbar[179]= - armWWbar[6] - 1./6.*armWWbar[36] + armWWbar[179];
   armWWbar[180]=armWWbar[1]*armWWbar[179];
   armWWbar[179]=armWWbar[40]*armWWbar[179];
   armWWbar[181]= - MMt*armWWbar[3]*armWWbar[37];
   armWWbar[135]=3./2.*armWWbar[181] + armWWbar[135] + armWWbar[161] + 
   armWWbar[236] + armWWbar[179] + 1./3.*armWWbar[180] + 33./4.*
   armWWbar[176] + 1./4.*armWWbar[160] - armWWbar[10];
   armWWbar[135]=MMt*armWWbar[135];
   armWWbar[160]=armWWbar[11]*armWWbar[205];
   armWWbar[161]= - 2*armWWbar[11] - 1./2. + armWWbar[10];
   armWWbar[161]=armWWbar[37]*armWWbar[161];
   armWWbar[160]=armWWbar[161] + armWWbar[160] + armWWbar[186] + 1./4.*
   armWWbar[36] + armWWbar[196];
   armWWbar[160]=MMH*armWWbar[160];
   armWWbar[161]= - 33./4.*armWWbar[12];
   armWWbar[176]=armWWbar[161] + 33./4.*armWWbar[13] - armWWbar[10];
   armWWbar[179]=armWWbar[11] + armWWbar[183] + armWWbar[176] + 
   armWWbar[211];
   armWWbar[179]=armWWbar[38]*armWWbar[179];
   armWWbar[180]=armWWbar[115] + 7./3. + armWWbar[174];
   armWWbar[180]=1./2.*armWWbar[180] - armWWbar[37];
   armWWbar[180]=armWWbar[37]*armWWbar[180];
   armWWbar[181]=armWWbar[36] + armWWbar[6];
   armWWbar[182]=armWWbar[181] - 2*armWWbar[37];
   armWWbar[182]=armWWbar[3]*armWWbar[38]*armWWbar[182];
   armWWbar[126]=3*armWWbar[182] + armWWbar[180] - 2./3.*armWWbar[36]
    + armWWbar[126];
   armWWbar[126]=MMt*armWWbar[126];
   armWWbar[180]=1./2.*armWWbar[187] + armWWbar[37];
   armWWbar[182]=armWWbar[38]*armWWbar[180];
   armWWbar[126]=armWWbar[182] + armWWbar[126];
   armWWbar[126]=armWWbar[34]*armWWbar[126];
   armWWbar[158]=armWWbar[278] + 17./9. + armWWbar[158];
   armWWbar[158]=armWWbar[19]*armWWbar[158];
   armWWbar[158]=armWWbar[193] + armWWbar[158];
   armWWbar[182]=armWWbar[289] - 8./3.*armWWbar[19];
   armWWbar[182]=armWWbar[37]*armWWbar[182];
   armWWbar[134]=13./2.*armWWbar[6] - 17./3. + armWWbar[134];
   armWWbar[134]=1./4.*armWWbar[134] + armWWbar[37];
   armWWbar[134]=armWWbar[20]*armWWbar[134];
   armWWbar[126]=armWWbar[126] + armWWbar[135] + 1./3.*armWWbar[160] + 
   armWWbar[179] + 1./3.*armWWbar[134] + 1./4.*armWWbar[158] + 
   armWWbar[182];
   armWWbar[126]=armWWbar[34]*armWWbar[126];
   armWWbar[134]= - 11*armWWbar[13];
   armWWbar[135]=127./27. + armWWbar[134];
   armWWbar[158]=11./2.*armWWbar[12];
   armWWbar[135]=armWWbar[299] + armWWbar[158] + 1./4.*armWWbar[135] + 
   2./3.*armWWbar[10];
   armWWbar[135]=armWWbar[11]*armWWbar[135];
   armWWbar[160]= - 11./2.*armWWbar[13];
   armWWbar[179]=armWWbar[12]*armWWbar[184];
   armWWbar[179]=11*armWWbar[179] + armWWbar[160] - 127./27.*
   armWWbar[10];
   armWWbar[135]=1./4.*armWWbar[179] + armWWbar[135];
   armWWbar[135]=MMH*armWWbar[135];
   armWWbar[179]=1057./27. - 165./2.*armWWbar[13];
   armWWbar[182]=7*armWWbar[10];
   armWWbar[179]=1./2.*armWWbar[179] + armWWbar[182];
   armWWbar[184]=armWWbar[113] - armWWbar[112];
   armWWbar[186]=armWWbar[19]*armWWbar[184];
   armWWbar[179]=3./8.*armWWbar[186] + 1./4.*armWWbar[179] + 22*
   armWWbar[12];
   armWWbar[179]=armWWbar[19]*armWWbar[179];
   armWWbar[190]= - 1057./27. - 143./2.*armWWbar[13];
   armWWbar[186]=3*armWWbar[186] - 11*armWWbar[12] + 1./2.*
   armWWbar[190] + 25./3.*armWWbar[10];
   armWWbar[186]=1./4.*armWWbar[186] + armWWbar[299];
   armWWbar[186]=armWWbar[20]*armWWbar[186];
   armWWbar[190]=armWWbar[265] + armWWbar[20];
   armWWbar[190]=armWWbar[20]*armWWbar[190];
   armWWbar[190]= - 1./2.*armWWbar[279] + armWWbar[190];
   armWWbar[190]=armWWbar[3]*armWWbar[190];
   armWWbar[191]=armWWbar[18]*armWWbar[176];
   armWWbar[192]=armWWbar[18] - 19./3.*armWWbar[19];
   armWWbar[192]=armWWbar[11]*armWWbar[192];
   armWWbar[135]=armWWbar[135] + 17./4.*armWWbar[190] + armWWbar[186]
    + 1./2.*armWWbar[192] + 1./2.*armWWbar[191] + armWWbar[179];
   armWWbar[135]=armWWbar[76]*armWWbar[135];
   armWWbar[179]=armWWbar[10]*armWWbar[220];
   armWWbar[115]=armWWbar[115] - 1./2.*armWWbar[7] + armWWbar[179];
   armWWbar[179]=armWWbar[1]*armWWbar[115];
   armWWbar[115]=armWWbar[40]*armWWbar[115];
   armWWbar[186]= - 2*armWWbar[7];
   armWWbar[190]=armWWbar[6] + 1./3. + armWWbar[186];
   armWWbar[191]=armWWbar[1]*armWWbar[190];
   armWWbar[190]=armWWbar[40]*armWWbar[190];
   armWWbar[190]=1./3.*armWWbar[191] + armWWbar[190];
   armWWbar[190]=armWWbar[11]*armWWbar[190];
   armWWbar[115]=armWWbar[190] + 1./3.*armWWbar[179] + armWWbar[115];
   armWWbar[115]=MMH*armWWbar[115];
   armWWbar[179]=17./12. - 8*armWWbar[7];
   armWWbar[179]=1./3.*armWWbar[179] + armWWbar[261];
   armWWbar[190]=armWWbar[1]*armWWbar[179];
   armWWbar[179]=armWWbar[40]*armWWbar[179];
   armWWbar[179]=1./3.*armWWbar[190] + armWWbar[179];
   armWWbar[179]=armWWbar[19]*armWWbar[179];
   armWWbar[190]=13./4.*armWWbar[6] - 17./12. + armWWbar[7];
   armWWbar[191]=armWWbar[1]*armWWbar[190];
   armWWbar[190]=armWWbar[40]*armWWbar[190];
   armWWbar[190]=1./3.*armWWbar[191] + armWWbar[190];
   armWWbar[190]=armWWbar[20]*armWWbar[190];
   armWWbar[115]=armWWbar[126] + armWWbar[135] + 1./3.*armWWbar[115] + 
   1./3.*armWWbar[190] + 1./2.*armWWbar[214] + armWWbar[179];
   armWWbar[115]=armWWbar[165]*armWWbar[115];
   armWWbar[126]= - 1./2.*armWWbar[8];
   armWWbar[135]= - armWWbar[9] + armWWbar[126];
   armWWbar[179]= - 1./2.*armWWbar[21];
   armWWbar[135]= - 35./2.*armWWbar[23] + armWWbar[179] + 11*
   armWWbar[135] - 13*armWWbar[22];
   armWWbar[190]=armWWbar[1]*armWWbar[135];
   armWWbar[135]=armWWbar[40]*armWWbar[135];
   armWWbar[191]=1./2. - armWWbar[7];
   armWWbar[191]=1./2.*armWWbar[191] + armWWbar[6];
   armWWbar[192]=armWWbar[1]*armWWbar[191];
   armWWbar[191]=armWWbar[40]*armWWbar[191];
   armWWbar[191]=1./3.*armWWbar[192] + armWWbar[191];
   armWWbar[191]=armWWbar[18]*armWWbar[191];
   armWWbar[192]=619./6. - 71*armWWbar[7];
   armWWbar[192]=1./3.*armWWbar[192] + armWWbar[6];
   armWWbar[193]=armWWbar[1]*armWWbar[192];
   armWWbar[192]=armWWbar[40]*armWWbar[192];
   armWWbar[192]=1./3.*armWWbar[193] + armWWbar[192];
   armWWbar[192]=armWWbar[19]*armWWbar[192];
   armWWbar[135]=1./2.*armWWbar[192] + armWWbar[191] + 1./3.*
   armWWbar[190] + armWWbar[135];
   armWWbar[190]=1./6.*armWWbar[14];
   armWWbar[191]= - armWWbar[16] + armWWbar[190] + armWWbar[28] - 7./9.
    - armWWbar[31];
   armWWbar[192]= - 1./3.*armWWbar[6];
   armWWbar[191]=1./2.*armWWbar[191] + armWWbar[192];
   armWWbar[193]=armWWbar[1]*armWWbar[191];
   armWWbar[191]=armWWbar[40]*armWWbar[191];
   armWWbar[191]=1./3.*armWWbar[193] + armWWbar[191];
   armWWbar[193]= - 1./4.*armWWbar[7];
   armWWbar[192]=armWWbar[192] - 2./9. + armWWbar[193];
   armWWbar[194]=armWWbar[1]*armWWbar[192];
   armWWbar[192]=armWWbar[40]*armWWbar[192];
   armWWbar[192]=1./3.*armWWbar[194] + armWWbar[192];
   armWWbar[192]=armWWbar[11]*armWWbar[192];
   armWWbar[191]=1./2.*armWWbar[191] + armWWbar[192];
   armWWbar[191]=MMH*armWWbar[191];
   armWWbar[186]= - 79./36.*armWWbar[6] + 1211./216. + armWWbar[186];
   armWWbar[186]=armWWbar[1]*armWWbar[186];
   armWWbar[192]= - 79./12.*armWWbar[6] + 1211./72. - 6*armWWbar[7];
   armWWbar[192]=armWWbar[40]*armWWbar[192];
   armWWbar[186]=armWWbar[186] + armWWbar[192];
   armWWbar[186]=armWWbar[20]*armWWbar[186];
   armWWbar[115]=armWWbar[115] + armWWbar[117] + armWWbar[132] + 
   armWWbar[191] + 1./2.*armWWbar[135] + armWWbar[186];
   armWWbar[115]=armWWbar[165]*armWWbar[115];
   armWWbar[117]= - 7*armWWbar[9];
   armWWbar[132]= - 1./6.*armWWbar[21];
   armWWbar[135]= - 5./18.*armWWbar[23] + armWWbar[132] + 19./9.*
   armWWbar[22] + armWWbar[117] - 17./18.*armWWbar[8];
   armWWbar[135]=armWWbar[1]*armWWbar[135];
   armWWbar[186]= - 21*armWWbar[9];
   armWWbar[191]= - 5./6.*armWWbar[23] + armWWbar[179] + 19./3.*
   armWWbar[22] + armWWbar[186] - 17./6.*armWWbar[8];
   armWWbar[191]=armWWbar[40]*armWWbar[191];
   armWWbar[192]=19./6. - armWWbar[7];
   armWWbar[192]=1./6.*armWWbar[192] - armWWbar[6];
   armWWbar[192]=armWWbar[1]*armWWbar[192];
   armWWbar[194]=107./54. - armWWbar[7];
   armWWbar[194]=1./2.*armWWbar[194] - 11./9.*armWWbar[6];
   armWWbar[194]=armWWbar[40]*armWWbar[194];
   armWWbar[192]=armWWbar[192] + armWWbar[194];
   armWWbar[192]=armWWbar[18]*armWWbar[192];
   armWWbar[135]=armWWbar[192] + armWWbar[135] + armWWbar[191];
   armWWbar[191]= - armWWbar[16] + armWWbar[190] + armWWbar[28] - 103./
   81. - armWWbar[31];
   armWWbar[194]=11./27.*armWWbar[6];
   armWWbar[191]=1./2.*armWWbar[191] + armWWbar[194];
   armWWbar[191]=armWWbar[40]*armWWbar[191];
   armWWbar[190]= - armWWbar[16] + armWWbar[190] + armWWbar[28] - 5./3.
    - armWWbar[31];
   armWWbar[190]=1./2.*armWWbar[190] + armWWbar[6];
   armWWbar[190]=armWWbar[1]*armWWbar[190];
   armWWbar[190]=1./3.*armWWbar[190] + armWWbar[191];
   armWWbar[191]=armWWbar[194] - 38./81. + armWWbar[193];
   armWWbar[191]=armWWbar[40]*armWWbar[191];
   armWWbar[193]=armWWbar[6] - 2./3. + armWWbar[193];
   armWWbar[193]=armWWbar[1]*armWWbar[193];
   armWWbar[191]=1./3.*armWWbar[193] + armWWbar[191];
   armWWbar[191]=armWWbar[11]*armWWbar[191];
   armWWbar[190]=1./2.*armWWbar[190] + armWWbar[191];
   armWWbar[190]=MMH*armWWbar[190];
   armWWbar[191]=49*armWWbar[7];
   armWWbar[193]= - 115./6. + armWWbar[191];
   armWWbar[193]=1./12.*armWWbar[193] + armWWbar[330];
   armWWbar[193]=armWWbar[1]*armWWbar[193];
   armWWbar[191]= - 185./18. + armWWbar[191];
   armWWbar[191]=1./4.*armWWbar[191] + 25./3.*armWWbar[6];
   armWWbar[191]=armWWbar[40]*armWWbar[191];
   armWWbar[191]=armWWbar[193] + armWWbar[191];
   armWWbar[191]=armWWbar[19]*armWWbar[191];
   armWWbar[193]=17*armWWbar[7];
   armWWbar[194]=1471./24. + armWWbar[193];
   armWWbar[194]=1./3.*armWWbar[194] + armWWbar[330];
   armWWbar[194]=armWWbar[1]*armWWbar[194];
   armWWbar[193]=83./9.*armWWbar[6] + 13655./216. + armWWbar[193];
   armWWbar[193]=armWWbar[40]*armWWbar[193];
   armWWbar[193]=armWWbar[194] + armWWbar[193];
   armWWbar[193]=armWWbar[20]*armWWbar[193];
   armWWbar[115]=armWWbar[115] + armWWbar[125] + armWWbar[130] + 
   armWWbar[190] + 1./3.*armWWbar[193] + 1./2.*armWWbar[135] + 1./3.*
   armWWbar[191];
   armWWbar[115]=armWWbar[165]*armWWbar[115];
   armWWbar[125]= - 1081./27. + armWWbar[234];
   armWWbar[125]=1./18.*armWWbar[13] + armWWbar[133] + 1./4.*
   armWWbar[125] + armWWbar[172];
   armWWbar[125]=armWWbar[306] + armWWbar[136] + 1./2.*armWWbar[125] + 
   armWWbar[199];
   armWWbar[125]=armWWbar[11]*armWWbar[125];
   armWWbar[130]=341./576. + armWWbar[189];
   armWWbar[125]=armWWbar[275] + armWWbar[125] + armWWbar[169] + 
   armWWbar[273] + 1./72.*armWWbar[13] + armWWbar[266] + armWWbar[256]
    + armWWbar[290] + armWWbar[254] + armWWbar[286] + armWWbar[250] + 
   armWWbar[243] + armWWbar[242] + armWWbar[285] + armWWbar[237] + 
   armWWbar[228] + armWWbar[218] + armWWbar[217] + armWWbar[216] + 
   armWWbar[280] + armWWbar[213] + 1./3.*armWWbar[130] + armWWbar[204];
   armWWbar[125]=MMH*armWWbar[125];
   armWWbar[130]=7487./18. + armWWbar[203];
   armWWbar[130]=armWWbar[303] - 11./6.*armWWbar[11] + armWWbar[301] + 
   armWWbar[300] - 113./72.*armWWbar[12] + armWWbar[295] - 307./72.*
   armWWbar[13] + armWWbar[293] + 1./8.*armWWbar[130] + armWWbar[304];
   armWWbar[130]=armWWbar[20]*armWWbar[130];
   armWWbar[135]=armWWbar[305] - armWWbar[10] + armWWbar[225] + 
   armWWbar[173] + armWWbar[48] - 277./9. + armWWbar[316];
   armWWbar[135]=1./2.*armWWbar[135] + armWWbar[309];
   armWWbar[135]=armWWbar[18]*armWWbar[135];
   armWWbar[136]=65./6.*armWWbar[13] + armWWbar[277] + armWWbar[232] - 
   7691./27. + armWWbar[270];
   armWWbar[136]=armWWbar[323] + armWWbar[320] - 305./16.*armWWbar[12]
    + 1./4.*armWWbar[136] + armWWbar[313];
   armWWbar[136]=armWWbar[19]*armWWbar[136];
   armWWbar[169]=armWWbar[198] + 9*armWWbar[20];
   armWWbar[169]=armWWbar[20]*armWWbar[169];
   armWWbar[162]=armWWbar[162] + armWWbar[169];
   armWWbar[162]=armWWbar[3]*armWWbar[162];
   armWWbar[169]= - armWWbar[18] + 7./2.*armWWbar[19];
   armWWbar[169]=armWWbar[11]*armWWbar[169];
   armWWbar[125]=armWWbar[125] + armWWbar[162] + armWWbar[130] + 1./6.*
   armWWbar[169] + armWWbar[136] + armWWbar[135] - 559./144.*
   armWWbar[23] + armWWbar[340] - 2315./48.*armWWbar[22] + 
   armWWbar[339] + 365./24.*armWWbar[107] + armWWbar[338] + 467./24.*
   armWWbar[80] + armWWbar[337] + armWWbar[336] - 398./9.*armWWbar[108]
    + armWWbar[335];
   armWWbar[125]=armWWbar[76]*armWWbar[125];
   armWWbar[130]=1./48. + armWWbar[94];
   armWWbar[130]= - 17./8.*armWWbar[92] - 17./8.*armWWbar[93] + 
   armWWbar[150] + 1./2.*armWWbar[130] + armWWbar[99];
   armWWbar[135]=armWWbar[12] - 7./6. - 13*armWWbar[10];
   armWWbar[135]=armWWbar[12]*armWWbar[135];
   armWWbar[136]= - 29./72.*armWWbar[12] + 1./2. + armWWbar[46];
   armWWbar[136]=armWWbar[11]*armWWbar[136];
   armWWbar[150]= - 1./2.*armWWbar[109] - armWWbar[85];
   armWWbar[162]=armWWbar[150] - 7./3.*armWWbar[82];
   armWWbar[162]=MMH*armWWbar[162];
   armWWbar[130]=1./4.*armWWbar[162] + armWWbar[136] + 1./12.*
   armWWbar[135] + 5./8.*armWWbar[10] + armWWbar[321] + 23./24.*
   armWWbar[16] + 1./3.*armWWbar[17] + 16./9.*armWWbar[100] - 1./12.*
   armWWbar[96] - 155./48.*armWWbar[103] + 1./3.*armWWbar[130] - 1./8.*
   armWWbar[98];
   armWWbar[130]=MMH*armWWbar[130];
   armWWbar[135]=9./2.*armWWbar[13] + armWWbar[277] + armWWbar[232] + 
   4513./27. + armWWbar[270];
   armWWbar[135]=59./4.*armWWbar[12] + 1./2.*armWWbar[135] + 
   armWWbar[182];
   armWWbar[136]=armWWbar[18]*armWWbar[113];
   armWWbar[162]=2./3.*armWWbar[113] - 3./8.*armWWbar[112];
   armWWbar[162]=armWWbar[19]*armWWbar[162];
   armWWbar[135]=armWWbar[162] + 1./4.*armWWbar[135] + 1./3.*
   armWWbar[136];
   armWWbar[135]=armWWbar[19]*armWWbar[135];
   armWWbar[136]=1./3.*armWWbar[297] - 653./36.*armWWbar[12] - 15*
   armWWbar[10] - 31./12.*armWWbar[13] + 5227./108. + armWWbar[324];
   armWWbar[162]=2*armWWbar[113] - 3./4.*armWWbar[112];
   armWWbar[162]=armWWbar[19]*armWWbar[162];
   armWWbar[136]=1./12.*armWWbar[177] - 1./12.*armWWbar[11] + 1./4.*
   armWWbar[136] + armWWbar[162];
   armWWbar[136]=armWWbar[20]*armWWbar[136];
   armWWbar[162]=armWWbar[195] - armWWbar[20];
   armWWbar[162]=armWWbar[20]*armWWbar[162];
   armWWbar[156]=3./2.*armWWbar[156] + armWWbar[162];
   armWWbar[156]=armWWbar[3]*armWWbar[156];
   armWWbar[162]= - 2*armWWbar[108];
   armWWbar[169]=armWWbar[162] + 7./8.*armWWbar[106];
   armWWbar[169]=2051./24.*armWWbar[23] - 73./24.*armWWbar[21] + 815./
   16.*armWWbar[22] - 25./4.*armWWbar[78] - 151./8.*armWWbar[107] - 1./
   4.*armWWbar[79] + 17*armWWbar[169] - 1031./48.*armWWbar[80];
   armWWbar[177]= - armWWbar[12] - 11./3. + armWWbar[157];
   armWWbar[177]=armWWbar[18]*armWWbar[177];
   armWWbar[130]=armWWbar[130] + 1./4.*armWWbar[156] + armWWbar[136] + 
   5./12.*armWWbar[325] + armWWbar[135] + 1./3.*armWWbar[169] + 
   armWWbar[177];
   armWWbar[130]=armWWbar[76]*armWWbar[130];
   armWWbar[135]=3*armWWbar[98];
   armWWbar[136]=1./3.*armWWbar[103] + armWWbar[135] - armWWbar[92] + 
   61./36. - armWWbar[93];
   armWWbar[136]=armWWbar[314] - armWWbar[17] + armWWbar[100] + 1./2.*
   armWWbar[136] - 7./3.*armWWbar[96];
   armWWbar[156]= - 1./2. - armWWbar[10];
   armWWbar[156]=armWWbar[12]*armWWbar[156];
   armWWbar[169]=3*armWWbar[46];
   armWWbar[177]= - 1./6.*armWWbar[12] + 1 + armWWbar[169];
   armWWbar[177]=armWWbar[11]*armWWbar[177];
   armWWbar[189]=1./2.*armWWbar[109] + armWWbar[85];
   armWWbar[191]=armWWbar[189] + armWWbar[82];
   armWWbar[191]=MMH*armWWbar[191];
   armWWbar[136]=1./6.*armWWbar[191] + 1./2.*armWWbar[177] + 1./2.*
   armWWbar[136] + 1./3.*armWWbar[156];
   armWWbar[136]=MMH*armWWbar[136];
   armWWbar[156]= - 1./12.*armWWbar[12];
   armWWbar[177]=armWWbar[156] - 5./9. + armWWbar[173];
   armWWbar[177]=armWWbar[18]*armWWbar[177];
   armWWbar[120]=197./3. + armWWbar[120];
   armWWbar[120]=25./12.*armWWbar[12] + 1./6.*armWWbar[120] + 
   armWWbar[10];
   armWWbar[120]=1./2.*armWWbar[120] + 1./3.*armWWbar[315];
   armWWbar[191]=armWWbar[19]*armWWbar[113];
   armWWbar[120]=1./2.*armWWbar[120] + 2./3.*armWWbar[191];
   armWWbar[120]=armWWbar[19]*armWWbar[120];
   armWWbar[193]= - 1./6.*armWWbar[13] - 7./18. + armWWbar[324];
   armWWbar[194]= - 11./3.*armWWbar[12];
   armWWbar[193]=armWWbar[194] + 1./2.*armWWbar[193] - armWWbar[10];
   armWWbar[191]=1./4.*armWWbar[193] + armWWbar[191];
   armWWbar[191]=armWWbar[20]*armWWbar[191];
   armWWbar[162]=armWWbar[162] + 41./48.*armWWbar[106];
   armWWbar[120]=1./4.*armWWbar[136] + armWWbar[191] + 1./4.*
   armWWbar[325] + armWWbar[120] + 1./4.*armWWbar[177] + 355./144.*
   armWWbar[23] + 43./144.*armWWbar[21] + 121./48.*armWWbar[22] - 11./
   24.*armWWbar[78] - 7./12.*armWWbar[107] + 1./3.*armWWbar[162] - 3./4.
   *armWWbar[80];
   armWWbar[120]=armWWbar[76]*armWWbar[120];
   armWWbar[136]=7./3. + armWWbar[6];
   armWWbar[136]=armWWbar[1]*armWWbar[136];
   armWWbar[162]=29./3. + 11*armWWbar[6];
   armWWbar[162]=armWWbar[40]*armWWbar[162];
   armWWbar[136]=armWWbar[136] + 1./9.*armWWbar[162];
   armWWbar[136]=armWWbar[19]*armWWbar[136];
   armWWbar[162]=7./6. - armWWbar[6];
   armWWbar[177]=armWWbar[1]*armWWbar[162];
   armWWbar[162]=armWWbar[40]*armWWbar[162];
   armWWbar[162]=armWWbar[177] + 11./9.*armWWbar[162];
   armWWbar[162]=armWWbar[20]*armWWbar[162];
   armWWbar[177]= - armWWbar[23] - 3*armWWbar[8] + armWWbar[22];
   armWWbar[177]=armWWbar[1]*armWWbar[177];
   armWWbar[191]=11*armWWbar[22];
   armWWbar[193]= - 11*armWWbar[23] - 17*armWWbar[8] + armWWbar[191];
   armWWbar[193]=armWWbar[40]*armWWbar[193];
   armWWbar[136]=armWWbar[162] + armWWbar[136] + armWWbar[177] + 1./9.*
   armWWbar[193];
   armWWbar[162]=5./2.*armWWbar[22] + armWWbar[227] - armWWbar[107];
   armWWbar[162]=1./2.*armWWbar[162] + armWWbar[23];
   armWWbar[177]=1./6.*armWWbar[12];
   armWWbar[193]=1 + armWWbar[177];
   armWWbar[193]=armWWbar[19]*armWWbar[193];
   armWWbar[196]= - 1./3.*armWWbar[12];
   armWWbar[198]= - 1 + armWWbar[196];
   armWWbar[198]=armWWbar[20]*armWWbar[198];
   armWWbar[203]=1 + armWWbar[103];
   armWWbar[203]=MMH*armWWbar[203];
   armWWbar[162]=1./12.*armWWbar[203] + 1./4.*armWWbar[198] + 1./3.*
   armWWbar[162] + 1./2.*armWWbar[193];
   armWWbar[162]=armWWbar[255]*armWWbar[76]*armWWbar[162];
   armWWbar[120]=1./4.*armWWbar[162] + 1./4.*armWWbar[136] + 
   armWWbar[120];
   armWWbar[120]=armWWbar[255]*armWWbar[120];
   armWWbar[136]= - 5./3.*armWWbar[23] - 4./3.*armWWbar[8] + 
   armWWbar[22];
   armWWbar[136]=armWWbar[1]*armWWbar[136];
   armWWbar[162]= - 55*armWWbar[23] - 40*armWWbar[8] + 37*armWWbar[22];
   armWWbar[162]=armWWbar[40]*armWWbar[162];
   armWWbar[136]=armWWbar[136] + 1./27.*armWWbar[162];
   armWWbar[162]=67./12. - armWWbar[7];
   armWWbar[162]=1./3.*armWWbar[162] + 9./4.*armWWbar[6];
   armWWbar[162]=armWWbar[1]*armWWbar[162];
   armWWbar[193]=401./36. - armWWbar[7];
   armWWbar[193]=1./9.*armWWbar[193] + 11./4.*armWWbar[6];
   armWWbar[193]=armWWbar[40]*armWWbar[193];
   armWWbar[162]=armWWbar[162] + armWWbar[193];
   armWWbar[162]=armWWbar[19]*armWWbar[162];
   armWWbar[193]=59./9. - 13./2.*armWWbar[6];
   armWWbar[193]=armWWbar[1]*armWWbar[193];
   armWWbar[198]=641./9. - 143./2.*armWWbar[6];
   armWWbar[198]=armWWbar[40]*armWWbar[198];
   armWWbar[193]=armWWbar[193] + 1./9.*armWWbar[198];
   armWWbar[193]=armWWbar[20]*armWWbar[193];
   armWWbar[120]=armWWbar[120] + armWWbar[130] + 1./2.*armWWbar[193] + 
   2*armWWbar[136] + armWWbar[162];
   armWWbar[120]=armWWbar[255]*armWWbar[120];
   armWWbar[117]=91./18.*armWWbar[23] + armWWbar[132] + 13./9.*
   armWWbar[22] + armWWbar[117] - 65./18.*armWWbar[8];
   armWWbar[117]=armWWbar[1]*armWWbar[117];
   armWWbar[130]=145./18.*armWWbar[23] + armWWbar[179] + 31./9.*
   armWWbar[22] + armWWbar[186] - 131./18.*armWWbar[8];
   armWWbar[130]=armWWbar[40]*armWWbar[130];
   armWWbar[117]=armWWbar[192] + armWWbar[117] + armWWbar[130];
   armWWbar[130]= - 139./6. - 11*armWWbar[7];
   armWWbar[130]=1./12.*armWWbar[130] + armWWbar[6];
   armWWbar[130]=armWWbar[1]*armWWbar[130];
   armWWbar[132]=5*armWWbar[7];
   armWWbar[136]= - 11./6. + armWWbar[132];
   armWWbar[136]=1./4.*armWWbar[136] + 5./3.*armWWbar[6];
   armWWbar[136]=armWWbar[40]*armWWbar[136];
   armWWbar[130]=armWWbar[130] + armWWbar[136];
   armWWbar[130]=armWWbar[19]*armWWbar[130];
   armWWbar[136]= - 929./24. + armWWbar[132];
   armWWbar[136]=1./9.*armWWbar[136] - armWWbar[6];
   armWWbar[136]=armWWbar[1]*armWWbar[136];
   armWWbar[132]= - 37./9.*armWWbar[6] + 2039./216. + armWWbar[132];
   armWWbar[132]=armWWbar[40]*armWWbar[132];
   armWWbar[132]=armWWbar[136] + 1./3.*armWWbar[132];
   armWWbar[132]=armWWbar[20]*armWWbar[132];
   armWWbar[117]=armWWbar[120] + armWWbar[125] + armWWbar[190] + 
   armWWbar[132] + 1./2.*armWWbar[117] + 1./3.*armWWbar[130];
   armWWbar[117]=armWWbar[255]*armWWbar[117];
   armWWbar[120]= - 4*armWWbar[108] + armWWbar[80];
   armWWbar[120]=4*armWWbar[23] + 2*armWWbar[120] - 17./3.*armWWbar[22]
   ;
   armWWbar[125]=4*armWWbar[12] + 1 + armWWbar[13];
   armWWbar[125]=armWWbar[19]*armWWbar[125];
   armWWbar[130]=5*armWWbar[12] + 8 + armWWbar[13];
   armWWbar[130]=armWWbar[20]*armWWbar[130];
   armWWbar[120]=2*armWWbar[130] + 2*armWWbar[120] + armWWbar[125];
   armWWbar[120]=armWWbar[76]*armWWbar[120];
   armWWbar[114]=armWWbar[114] + armWWbar[115] + armWWbar[171] + 4*
   armWWbar[120] + armWWbar[117];
   armWWbar[114]=armWWbar[4]*armWWbar[114];
   armWWbar[115]=2*armWWbar[36];
   armWWbar[117]=armWWbar[6] + 1 + armWWbar[115];
   armWWbar[117]=2*armWWbar[117] + armWWbar[268];
   armWWbar[117]=armWWbar[20]*armWWbar[117];
   armWWbar[120]= - armWWbar[42] - armWWbar[8];
   armWWbar[120]=armWWbar[22] + 3*armWWbar[51] + 1./2.*armWWbar[120] - 
   2*armWWbar[43];
   armWWbar[125]=6*armWWbar[23];
   armWWbar[130]=17./4.*armWWbar[6] - 12 + 25./4.*armWWbar[36];
   armWWbar[130]=MMZ*armWWbar[130];
   armWWbar[132]=7./4.*armWWbar[6] + 1 + 11./4.*armWWbar[36];
   armWWbar[132]=armWWbar[19]*armWWbar[132];
   armWWbar[136]=3*MMZ + armWWbar[289];
   armWWbar[136]=armWWbar[37]*armWWbar[136];
   armWWbar[162]=6*armWWbar[49] + 37./2. + armWWbar[123];
   armWWbar[162]=armWWbar[38]*armWWbar[162];
   armWWbar[117]=armWWbar[162] + armWWbar[117] + armWWbar[136] + 
   armWWbar[132] + armWWbar[185] + armWWbar[130] + armWWbar[125] + 3*
   armWWbar[120] + armWWbar[244];
   armWWbar[117]=armWWbar[3]*armWWbar[117];
   armWWbar[120]= - 2833./64. + armWWbar[123];
   armWWbar[120]=1./2.*armWWbar[120] + armWWbar[219];
   armWWbar[130]= - 22*armWWbar[13];
   armWWbar[132]= - 81./512.*armWWbar[111] + armWWbar[112];
   armWWbar[132]=MMZ*armWWbar[132];
   armWWbar[136]= - armWWbar[35] + 81./64.*armWWbar[111];
   armWWbar[136]=armWWbar[19]*armWWbar[136];
   armWWbar[162]= - MMZ*armWWbar[112];
   armWWbar[171]= - 1 + armWWbar[162];
   armWWbar[171]=armWWbar[11]*armWWbar[171];
   armWWbar[120]=5./4.*armWWbar[37] + armWWbar[171] + 1./4.*
   armWWbar[136] + 1./2.*armWWbar[297] + armWWbar[132] + 2*
   armWWbar[144] + 2./3.*armWWbar[159] + 33./4.*armWWbar[12] + 1./4.*
   armWWbar[120] + armWWbar[130];
   armWWbar[120]=armWWbar[37]*armWWbar[120];
   armWWbar[132]=7./4. + armWWbar[69];
   armWWbar[132]=1./2.*armWWbar[132] + armWWbar[71];
   armWWbar[136]=pow(armWWbar[110],2);
   armWWbar[132]=MMZ*armWWbar[136]*armWWbar[132];
   armWWbar[142]=3*armWWbar[142] - armWWbar[59];
   armWWbar[144]= - 1./2.*armWWbar[66] + 3./8. - armWWbar[67];
   armWWbar[144]=3./4.*armWWbar[50] + 1./2.*armWWbar[144] - 
   armWWbar[33];
   armWWbar[144]=armWWbar[35]*armWWbar[144];
   armWWbar[126]=armWWbar[74] + armWWbar[126];
   armWWbar[126]=armWWbar[110]*armWWbar[126];
   armWWbar[159]=3./2.*armWWbar[126] + 1./4.*armWWbar[71] - 103./16. - 
   armWWbar[69];
   armWWbar[159]=armWWbar[110]*armWWbar[159];
   armWWbar[171]=205 + 9*armWWbar[66];
   armWWbar[171]=11*armWWbar[17] - 9./8.*armWWbar[50] - 9./8.*
   armWWbar[71] + 9./4.*armWWbar[69] + 1./8.*armWWbar[171] - 11*
   armWWbar[63];
   armWWbar[171]=armWWbar[111]*armWWbar[171];
   armWWbar[185]=armWWbar[12]*armWWbar[111];
   armWWbar[186]= - armWWbar[16] + 6*armWWbar[68] + 15./4. + 
   armWWbar[60];
   armWWbar[186]=armWWbar[112]*armWWbar[186];
   armWWbar[132]=3./4.*armWWbar[132] + armWWbar[186] + 99./64.*
   armWWbar[185] + 9./64.*armWWbar[171] + 1./2.*armWWbar[159] + 2*
   armWWbar[142] + armWWbar[144];
   armWWbar[132]=MMZ*armWWbar[132];
   armWWbar[142]=7./2.*armWWbar[16];
   armWWbar[144]=armWWbar[131] + armWWbar[142] + armWWbar[263] + 
   armWWbar[70] + armWWbar[287] + 1./4. + armWWbar[230];
   armWWbar[144]=armWWbar[112]*armWWbar[144];
   armWWbar[123]=47./2. + armWWbar[123];
   armWWbar[159]= - 1 + armWWbar[48];
   armWWbar[159]=2*armWWbar[159] + armWWbar[46];
   armWWbar[159]=armWWbar[3]*armWWbar[38]*armWWbar[159];
   armWWbar[123]=9*armWWbar[159] + armWWbar[178] + 23./4.*armWWbar[36]
    - 9./2.*armWWbar[39] + 1./2.*armWWbar[123] + armWWbar[219];
   armWWbar[123]=armWWbar[3]*armWWbar[123];
   armWWbar[159]=armWWbar[56] + 1./2.*armWWbar[54];
   armWWbar[171]=2*armWWbar[59];
   armWWbar[178]=armWWbar[11]*armWWbar[112]*armWWbar[39];
   armWWbar[186]=MMt*armWWbar[274]*armWWbar[39];
   armWWbar[123]=18*armWWbar[186] + armWWbar[123] + armWWbar[283] + 3*
   armWWbar[178] + armWWbar[144] + 2*armWWbar[55] + 7*armWWbar[159] + 
   armWWbar[171];
   armWWbar[123]=MMt*armWWbar[123];
   armWWbar[144]= - 2*armWWbar[36];
   armWWbar[159]=1 + armWWbar[144];
   armWWbar[178]=armWWbar[7]*armWWbar[36];
   armWWbar[190]= - 4./3. + armWWbar[7];
   armWWbar[190]=armWWbar[6]*armWWbar[190];
   armWWbar[159]=armWWbar[190] + 1./3.*armWWbar[159] + 2*armWWbar[178];
   armWWbar[178]=armWWbar[1]*armWWbar[159];
   armWWbar[159]=armWWbar[40]*armWWbar[159];
   armWWbar[190]= - armWWbar[12]*armWWbar[111];
   armWWbar[192]=33./4.*armWWbar[190] + 561./32.*armWWbar[111] - 
   armWWbar[35] + 61*armWWbar[110];
   armWWbar[193]= - MMZ*armWWbar[136];
   armWWbar[198]= - 81./32.*armWWbar[111] + armWWbar[35] - 9./2.*
   armWWbar[110];
   armWWbar[198]=armWWbar[37]*armWWbar[198];
   armWWbar[192]=1./4.*armWWbar[198] + 9./8.*armWWbar[193] + 3./4.*
   armWWbar[192] - 2*armWWbar[112];
   armWWbar[192]=armWWbar[38]*armWWbar[192];
   armWWbar[198]= - 4*armWWbar[36];
   armWWbar[203]= - 2*armWWbar[6] + 1 + armWWbar[198];
   armWWbar[204]=armWWbar[115] + armWWbar[6];
   armWWbar[204]=2./3.*armWWbar[204] - armWWbar[37];
   armWWbar[204]=armWWbar[37]*armWWbar[204];
   armWWbar[203]=1./9.*armWWbar[203] + armWWbar[204];
   armWWbar[203]=armWWbar[34]*armWWbar[203];
   armWWbar[204]=34735./3456. + 33*armWWbar[61];
   armWWbar[205]=armWWbar[74] + armWWbar[319];
   armWWbar[206]=armWWbar[205] - armWWbar[43];
   armWWbar[206]=armWWbar[35]*armWWbar[206];
   armWWbar[211]=armWWbar[51]*armWWbar[35];
   armWWbar[213]=57./2.*armWWbar[51] + 9./8.*armWWbar[43] + 1./8.*
   armWWbar[8] - 1./4.*armWWbar[74] - 15*armWWbar[72];
   armWWbar[213]=armWWbar[110]*armWWbar[213];
   armWWbar[214]= - armWWbar[35] - 3./2.*armWWbar[110];
   armWWbar[214]=armWWbar[22]*armWWbar[214];
   armWWbar[216]= - armWWbar[35] + 111./4.*armWWbar[110];
   armWWbar[216]=armWWbar[23]*armWWbar[216];
   armWWbar[217]=11*armWWbar[23] + 33./16.*armWWbar[22] + 155./8.*
   armWWbar[51] - 3./8.*armWWbar[42] - 11*armWWbar[72];
   armWWbar[217]=armWWbar[111]*armWWbar[217];
   armWWbar[218]=22./3.*armWWbar[13];
   armWWbar[219]= - 293./24. - 9*armWWbar[36];
   armWWbar[219]=armWWbar[12]*armWWbar[219];
   armWWbar[222]= - 11./2.*armWWbar[12];
   armWWbar[224]= - 31./27. + armWWbar[222];
   armWWbar[224]=armWWbar[6]*armWWbar[224];
   armWWbar[227]=armWWbar[23] + armWWbar[179] + 2*armWWbar[73] - 
   armWWbar[43];
   armWWbar[227]=armWWbar[112]*armWWbar[227];
   armWWbar[228]=3*armWWbar[193] - 63./64.*armWWbar[111] - armWWbar[35]
    + armWWbar[110];
   armWWbar[228]=armWWbar[19]*armWWbar[228];
   armWWbar[230]=2./3.*armWWbar[6] + 1 + 4./3.*armWWbar[36];
   armWWbar[230]=armWWbar[11]*armWWbar[230];
   armWWbar[233]= - armWWbar[35] + 9./4.*armWWbar[110];
   armWWbar[233]=1./2.*armWWbar[233] + armWWbar[112];
   armWWbar[233]=armWWbar[37]*armWWbar[233];
   armWWbar[233]=armWWbar[233] - 33./32.*armWWbar[111] + armWWbar[35]
    - 39./4.*armWWbar[110];
   armWWbar[233]=armWWbar[20]*armWWbar[233];
   armWWbar[117]=armWWbar[203] + armWWbar[123] + armWWbar[117] + 
   armWWbar[192] + armWWbar[233] + armWWbar[120] + armWWbar[230] + 1./4.
   *armWWbar[228] + 5./2.*armWWbar[297] + armWWbar[132] + 2./3.*
   armWWbar[159] + 2./9.*armWWbar[178] + armWWbar[227] + 5./4.*
   armWWbar[224] + 11./8.*armWWbar[219] - 391./108.*armWWbar[36] + 
   armWWbar[218] + 3./16.*armWWbar[217] + 1./2.*armWWbar[216] + 
   armWWbar[16] + 1./4.*armWWbar[214] + armWWbar[213] - 957./64.*
   armWWbar[17] + armWWbar[211] + 1./2.*armWWbar[206] + armWWbar[124]
    - 3./4.*armWWbar[47] + 2*armWWbar[68] - 1681./512.*armWWbar[50] - 
   armWWbar[60] - 145./512.*armWWbar[71] - 207./256.*armWWbar[69] + 33./
   8.*armWWbar[30] - 21./8.*armWWbar[33] + 165./64.*armWWbar[63] + 657./
   512.*armWWbar[66] + 1./4.*armWWbar[204] + 2*armWWbar[67];
   armWWbar[117]=armWWbar[34]*armWWbar[117];
   armWWbar[120]= - 33./2.*armWWbar[12];
   armWWbar[123]=armWWbar[120] + armWWbar[163] - 825./4.*armWWbar[13]
    + armWWbar[173] + 133./4. + armWWbar[172];
   armWWbar[123]=MMZ*armWWbar[123];
   armWWbar[124]=3*armWWbar[10];
   armWWbar[132]=armWWbar[164] + armWWbar[124] - 363./4.*armWWbar[13]
    + armWWbar[173] - 631./12. + armWWbar[172];
   armWWbar[132]=armWWbar[19]*armWWbar[132];
   armWWbar[151]=armWWbar[151] + 1 - armWWbar[48];
   armWWbar[159]=pow(MMZ,2);
   armWWbar[178]=armWWbar[159]*armWWbar[151];
   armWWbar[151]=MMZ*armWWbar[151];
   armWWbar[192]=armWWbar[19]*armWWbar[151];
   armWWbar[151]=armWWbar[20]*armWWbar[151];
   armWWbar[151]=armWWbar[151] + armWWbar[178] + 1./2.*armWWbar[192];
   armWWbar[151]=armWWbar[3]*armWWbar[151];
   armWWbar[178]= - 99./2.*armWWbar[12] + armWWbar[173] + 115./2. + 
   armWWbar[172];
   armWWbar[178]=armWWbar[18]*armWWbar[178];
   armWWbar[192]= - 1./2.*armWWbar[78];
   armWWbar[203]=armWWbar[192] + 11*armWWbar[80] - armWWbar[79];
   armWWbar[203]=1./2.*armWWbar[203] - 5*armWWbar[22];
   armWWbar[204]=armWWbar[157] + 7./2. + armWWbar[48];
   armWWbar[204]=MMZ*armWWbar[204];
   armWWbar[204]=3./4.*armWWbar[19] + 3./2.*armWWbar[204] + 
   armWWbar[18];
   armWWbar[204]=armWWbar[11]*armWWbar[204];
   armWWbar[161]=9./2.*armWWbar[11] + armWWbar[161] - 66*armWWbar[13]
    + 9./4.*armWWbar[46] - 226./3. + 9./2.*armWWbar[48];
   armWWbar[161]=armWWbar[20]*armWWbar[161];
   armWWbar[123]=9*armWWbar[151] + armWWbar[161] + armWWbar[204] + 1./2.
   *armWWbar[132] + 1./4.*armWWbar[178] + 1./2.*armWWbar[123] - 45*
   armWWbar[23] + 9./2.*armWWbar[203] - 8*armWWbar[21];
   armWWbar[123]=armWWbar[3]*armWWbar[123];
   armWWbar[132]=11./4.*armWWbar[12];
   armWWbar[151]=armWWbar[132] + armWWbar[231] + armWWbar[317] + 1./2.
    + 3*armWWbar[44];
   armWWbar[151]=armWWbar[112]*armWWbar[151];
   armWWbar[151]=1./4.*armWWbar[113] + armWWbar[151];
   armWWbar[151]=MMZ*armWWbar[151];
   armWWbar[161]= - 145./3. + armWWbar[234];
   armWWbar[178]= - armWWbar[19]*armWWbar[113];
   armWWbar[162]= - 1./2. + armWWbar[162];
   armWWbar[162]=armWWbar[11]*armWWbar[162];
   armWWbar[151]=armWWbar[162] + 3./4.*armWWbar[178] + 3*armWWbar[151]
    + 11*armWWbar[12] + armWWbar[130] - armWWbar[46] + 1./4.*
   armWWbar[161] + armWWbar[304];
   armWWbar[151]=armWWbar[11]*armWWbar[151];
   armWWbar[161]= - 11*armWWbar[98];
   armWWbar[162]=11*armWWbar[103];
   armWWbar[203]= - armWWbar[96] + armWWbar[162] + 5./2. + 
   armWWbar[161];
   armWWbar[204]= - 11*armWWbar[100];
   armWWbar[206]=11./2.*armWWbar[17];
   armWWbar[203]=armWWbar[253] + armWWbar[206] + 1./2.*armWWbar[203] + 
   armWWbar[204];
   armWWbar[203]=armWWbar[113]*armWWbar[203];
   armWWbar[155]= - 11./4.*armWWbar[92] - 11./4.*armWWbar[93] - 
   armWWbar[97] + 16*armWWbar[101] - 5./16. + armWWbar[155];
   armWWbar[213]= - 3./4.*armWWbar[46];
   armWWbar[155]=armWWbar[213] + 19*armWWbar[16] + armWWbar[116] + 33./
   4.*armWWbar[17] + armWWbar[288] - 66*armWWbar[100] - 3./4.*
   armWWbar[96] - 33./4.*armWWbar[103] + 3*armWWbar[155] + 2*
   armWWbar[91];
   armWWbar[155]=armWWbar[112]*armWWbar[155];
   armWWbar[214]=armWWbar[12]*armWWbar[113];
   armWWbar[155]=armWWbar[155] + 33./4.*armWWbar[214] + 3./2.*
   armWWbar[203] - 33./2.*armWWbar[82] - 33./2.*armWWbar[85] - 33./4.*
   armWWbar[109] + 459./16.*armWWbar[83] + armWWbar[81];
   armWWbar[155]=MMZ*armWWbar[155];
   armWWbar[148]=armWWbar[154] + armWWbar[148] + 11./3. + armWWbar[149]
   ;
   armWWbar[148]=armWWbar[112]*armWWbar[148];
   armWWbar[149]= - armWWbar[12]*armWWbar[113];
   armWWbar[154]=11*armWWbar[149];
   armWWbar[203]=15*armWWbar[113] + armWWbar[154];
   armWWbar[148]=3./8.*armWWbar[203] + armWWbar[148];
   armWWbar[148]=armWWbar[18]*armWWbar[148];
   armWWbar[203]=armWWbar[213] + armWWbar[116] + 179./24. + 
   armWWbar[288];
   armWWbar[203]=armWWbar[112]*armWWbar[203];
   armWWbar[216]= - armWWbar[11]*armWWbar[112];
   armWWbar[203]=3./2.*armWWbar[216] + 7./8.*armWWbar[113] + 
   armWWbar[203];
   armWWbar[203]=armWWbar[20]*armWWbar[203];
   armWWbar[216]= - 435./8.*armWWbar[94] - 2413./81. - 9./2.*
   armWWbar[95];
   armWWbar[217]= - 1./6.*armWWbar[97];
   armWWbar[219]= - 3./2.*armWWbar[96];
   armWWbar[224]=5./8.*armWWbar[46];
   armWWbar[227]=armWWbar[294] - 21./4.*armWWbar[21] + 31./2.*
   armWWbar[22] + 7./2.*armWWbar[107] - 3*armWWbar[78];
   armWWbar[227]=armWWbar[113]*armWWbar[227];
   armWWbar[228]=6973./9. + 3267./2.*armWWbar[13];
   armWWbar[228]= - 381./2.*armWWbar[12] + 1./2.*armWWbar[228] - 11*
   armWWbar[10];
   armWWbar[228]=armWWbar[12]*armWWbar[228];
   armWWbar[230]= - 47./12.*armWWbar[23] + 181./24.*armWWbar[21] + 1./4.
   *armWWbar[22] - 3./8.*armWWbar[78] + 1./4.*armWWbar[107] - 11./3.*
   armWWbar[79] + 33./4.*armWWbar[80] - 33./2.*armWWbar[106] - 9./2.*
   armWWbar[77] + 11*armWWbar[105];
   armWWbar[230]=armWWbar[112]*armWWbar[230];
   armWWbar[233]=11*armWWbar[214];
   armWWbar[235]= - 13*armWWbar[113] + armWWbar[233];
   armWWbar[235]=1./2.*armWWbar[235] + 11*armWWbar[112];
   armWWbar[235]=armWWbar[19]*armWWbar[235];
   armWWbar[236]=13./2.*armWWbar[82] + 13./2.*armWWbar[85] + 13./4.*
   armWWbar[109] + 3./2.*armWWbar[84] - 2./3.*armWWbar[81];
   armWWbar[236]=MMH*armWWbar[236];
   armWWbar[123]=armWWbar[236] + armWWbar[123] + armWWbar[203] + 
   armWWbar[151] + 3./2.*armWWbar[235] + armWWbar[148] + armWWbar[155]
    + armWWbar[230] + 1./4.*armWWbar[228] + 4301./72.*armWWbar[13] + 1./
   2.*armWWbar[227] + armWWbar[224] + 7./6.*armWWbar[16] + 5./4.*
   armWWbar[48] + 391./8.*armWWbar[17] - 45./8.*armWWbar[44] - 33*
   armWWbar[100] + armWWbar[219] + 11./4.*armWWbar[98] + 11./4.*
   armWWbar[92] + 11./4.*armWWbar[93] + armWWbar[217] - 435./16.*
   armWWbar[99] + 1./2.*armWWbar[216] + 13*armWWbar[101];
   armWWbar[123]=armWWbar[76]*armWWbar[123];
   armWWbar[148]= - 59./9. + 33./2.*armWWbar[13];
   armWWbar[151]= - armWWbar[6] - 1./3. + 2*armWWbar[7];
   armWWbar[155]=1./3.*armWWbar[1]*armWWbar[151];
   armWWbar[151]=armWWbar[40]*armWWbar[151];
   armWWbar[175]=armWWbar[175] + armWWbar[151] + armWWbar[155] + 
   armWWbar[120] + 1./2.*armWWbar[148] - armWWbar[10];
   armWWbar[175]=armWWbar[37]*armWWbar[175];
   armWWbar[203]=MMZ*armWWbar[187];
   armWWbar[187]=armWWbar[19]*armWWbar[187];
   armWWbar[187]=armWWbar[203] + 1./2.*armWWbar[187];
   armWWbar[203]=MMZ + armWWbar[245];
   armWWbar[216]=armWWbar[37]*armWWbar[203];
   armWWbar[180]=armWWbar[20]*armWWbar[180];
   armWWbar[180]=armWWbar[180] + 1./2.*armWWbar[187] + armWWbar[216];
   armWWbar[180]=armWWbar[3]*armWWbar[180];
   armWWbar[187]=armWWbar[7]*armWWbar[208];
   armWWbar[208]=1 - armWWbar[7];
   armWWbar[216]=armWWbar[6]*armWWbar[208];
   armWWbar[187]=1./2.*armWWbar[216] + 1./6.*armWWbar[36] + 
   armWWbar[187];
   armWWbar[216]=armWWbar[1]*armWWbar[187];
   armWWbar[187]=armWWbar[40]*armWWbar[187];
   armWWbar[209]=armWWbar[209] + armWWbar[37];
   armWWbar[209]=armWWbar[37]*armWWbar[209];
   armWWbar[181]=1./6.*armWWbar[181] + armWWbar[209];
   armWWbar[181]=armWWbar[34]*armWWbar[181];
   armWWbar[134]=armWWbar[134] + 59./9.*armWWbar[36];
   armWWbar[174]=1 + armWWbar[174];
   armWWbar[174]=armWWbar[12]*armWWbar[174];
   armWWbar[209]=59./9. + armWWbar[164];
   armWWbar[209]=armWWbar[6]*armWWbar[209];
   armWWbar[170]=MMt*armWWbar[3]*armWWbar[170];
   armWWbar[134]=armWWbar[181] + 3./2.*armWWbar[170] + 3*armWWbar[180]
    + armWWbar[175] + armWWbar[210] + armWWbar[187] + 1./3.*
   armWWbar[216] + 1./4.*armWWbar[209] + 11./4.*armWWbar[174] + 1./4.*
   armWWbar[134] + armWWbar[199];
   armWWbar[134]=armWWbar[34]*armWWbar[134];
   armWWbar[148]=armWWbar[7]*armWWbar[148];
   armWWbar[148]=armWWbar[160] + armWWbar[148];
   armWWbar[160]=armWWbar[7]*armWWbar[220];
   armWWbar[121]=armWWbar[6]*armWWbar[121];
   armWWbar[160]=armWWbar[160] + armWWbar[121];
   armWWbar[170]=armWWbar[1]*armWWbar[160];
   armWWbar[160]=armWWbar[40]*armWWbar[160];
   armWWbar[174]=1./2. - 3*armWWbar[7];
   armWWbar[174]=armWWbar[12]*armWWbar[174];
   armWWbar[148]=armWWbar[160] + 2./3.*armWWbar[170] + 1./2.*
   armWWbar[209] + 11./2.*armWWbar[174] + 1./2.*armWWbar[148] + 
   armWWbar[168];
   armWWbar[148]=armWWbar[40]*armWWbar[148];
   armWWbar[119]= - 59./27. + armWWbar[119];
   armWWbar[119]=armWWbar[7]*armWWbar[119];
   armWWbar[119]= - 11./6.*armWWbar[13] + armWWbar[119];
   armWWbar[160]=armWWbar[12]*armWWbar[166];
   armWWbar[158]=59./27. + armWWbar[158];
   armWWbar[158]=armWWbar[6]*armWWbar[158];
   armWWbar[119]=1./9.*armWWbar[170] + 1./2.*armWWbar[158] + 11./2.*
   armWWbar[160] + 1./2.*armWWbar[119] + 1./3.*armWWbar[168];
   armWWbar[119]=armWWbar[1]*armWWbar[119];
   armWWbar[158]=59./3. - 99./2.*armWWbar[13];
   armWWbar[158]=99./4.*armWWbar[12] + 1./2.*armWWbar[158] + 
   armWWbar[124];
   armWWbar[158]=armWWbar[12]*armWWbar[158];
   armWWbar[160]=armWWbar[226] + 1./2.*armWWbar[107] - armWWbar[22];
   armWWbar[160]=armWWbar[113]*armWWbar[160];
   armWWbar[166]=1./2.*armWWbar[23] - 1./2.*armWWbar[107] + 
   armWWbar[22];
   armWWbar[166]=armWWbar[112]*armWWbar[166];
   armWWbar[168]=MMZ*armWWbar[184];
   armWWbar[170]= - armWWbar[113] + armWWbar[112];
   armWWbar[170]=armWWbar[19]*armWWbar[170];
   armWWbar[158]=9./2.*armWWbar[170] + 9./16.*armWWbar[168] + 3*
   armWWbar[166] + 11./2.*armWWbar[158] + 77./9.*armWWbar[10] + 3*
   armWWbar[160] - 649./12.*armWWbar[13];
   armWWbar[160]= - 17 + 99./2.*armWWbar[13];
   armWWbar[166]= - 99./4.*armWWbar[12];
   armWWbar[160]=armWWbar[166] + 1./2.*armWWbar[160] + armWWbar[153];
   armWWbar[160]=armWWbar[19]*armWWbar[160];
   armWWbar[168]=17 + 99*armWWbar[13];
   armWWbar[153]=3*armWWbar[11] + armWWbar[166] + 1./4.*armWWbar[168]
    + armWWbar[153];
   armWWbar[153]=armWWbar[20]*armWWbar[153];
   armWWbar[166]=MMZ*armWWbar[176];
   armWWbar[168]=armWWbar[11]*armWWbar[203];
   armWWbar[153]=armWWbar[153] + 3*armWWbar[168] + 3*armWWbar[166] + 1./
   2.*armWWbar[160];
   armWWbar[153]=armWWbar[3]*armWWbar[153];
   armWWbar[160]= - 7./9. + 3./2.*armWWbar[13];
   armWWbar[120]=armWWbar[11] + armWWbar[120] + 11./2.*armWWbar[160] - 
   armWWbar[10];
   armWWbar[120]=armWWbar[11]*armWWbar[120];
   armWWbar[160]=armWWbar[20]*armWWbar[184];
   armWWbar[120]=armWWbar[153] + 3./8.*armWWbar[160] + 1./2.*
   armWWbar[158] + armWWbar[120];
   armWWbar[120]=armWWbar[76]*armWWbar[120];
   armWWbar[153]=armWWbar[202] + 3*armWWbar[183];
   armWWbar[158]=MMZ*armWWbar[153];
   armWWbar[160]=armWWbar[19]*armWWbar[153];
   armWWbar[153]=armWWbar[20]*armWWbar[153];
   armWWbar[153]=armWWbar[153] + armWWbar[158] + 1./2.*armWWbar[160];
   armWWbar[153]=armWWbar[3]*armWWbar[153];
   armWWbar[151]=armWWbar[155] + armWWbar[151];
   armWWbar[151]=armWWbar[11]*armWWbar[151];
   armWWbar[119]=armWWbar[134] + armWWbar[120] + armWWbar[153] + 
   armWWbar[151] + armWWbar[119] + armWWbar[148];
   armWWbar[119]=armWWbar[165]*armWWbar[119];
   armWWbar[120]=armWWbar[6]*armWWbar[220];
   armWWbar[134]=pow(Pi,2);
   armWWbar[148]=pow(armWWbar[7],2);
   armWWbar[120]=2*armWWbar[120] - armWWbar[148] + 1./9. + 
   armWWbar[134];
   armWWbar[151]=armWWbar[1]*armWWbar[120];
   armWWbar[153]= - 5./3.*armWWbar[134] - 37./3.*armWWbar[32] - 3989./
   648. + 11*armWWbar[29];
   armWWbar[155]= - 1./3.*armWWbar[28];
   armWWbar[158]=1./3.*armWWbar[16];
   armWWbar[160]= - 5./9. + 1./4.*armWWbar[7];
   armWWbar[160]=armWWbar[12]*armWWbar[160];
   armWWbar[130]=5./4.*armWWbar[7] - 23./4. + armWWbar[130];
   armWWbar[130]=armWWbar[7]*armWWbar[130];
   armWWbar[166]= - 13./9. + armWWbar[222];
   armWWbar[166]=armWWbar[6]*armWWbar[166];
   armWWbar[168]=armWWbar[23] - armWWbar[9] + armWWbar[179];
   armWWbar[168]=armWWbar[112]*armWWbar[168];
   armWWbar[153]=1./9.*armWWbar[151] + 1./3.*armWWbar[168] + 7./6.*
   armWWbar[166] + 11*armWWbar[160] + 1./3.*armWWbar[130] + 22./9.*
   armWWbar[13] + armWWbar[158] - 11./2.*armWWbar[17] + 11./4.*
   armWWbar[30] - 19./4.*armWWbar[33] - 23./72.*armWWbar[14] + 1./4.*
   armWWbar[153] + armWWbar[155];
   armWWbar[153]=armWWbar[1]*armWWbar[153];
   armWWbar[120]=armWWbar[40]*armWWbar[120];
   armWWbar[160]= - 5*armWWbar[134] - 37*armWWbar[32] - 3989./216. + 33
   *armWWbar[29];
   armWWbar[170]= - 5./3. + 3./4.*armWWbar[7];
   armWWbar[170]=armWWbar[12]*armWWbar[170];
   armWWbar[120]=armWWbar[120] + 2./3.*armWWbar[151] + armWWbar[168] + 
   7./2.*armWWbar[166] + 11*armWWbar[170] + armWWbar[130] + 
   armWWbar[218] + armWWbar[16] - 33./2.*armWWbar[17] + 33./4.*
   armWWbar[30] - 57./4.*armWWbar[33] - 23./24.*armWWbar[14] + 1./4.*
   armWWbar[160] - armWWbar[28];
   armWWbar[120]=armWWbar[40]*armWWbar[120];
   armWWbar[130]=armWWbar[22] - 2*armWWbar[9] - armWWbar[8];
   armWWbar[125]=armWWbar[125] + 3*armWWbar[130] + armWWbar[244];
   armWWbar[125]=armWWbar[40]*armWWbar[125];
   armWWbar[151]=7./2.*armWWbar[6];
   armWWbar[160]=armWWbar[151] - 4 + armWWbar[7];
   armWWbar[166]=armWWbar[1]*armWWbar[160];
   armWWbar[160]=armWWbar[40]*armWWbar[160];
   armWWbar[160]=armWWbar[166] + 3*armWWbar[160];
   armWWbar[160]=MMZ*armWWbar[160];
   armWWbar[166]= - 1 + armWWbar[7];
   armWWbar[168]=armWWbar[1]*armWWbar[166];
   armWWbar[166]=armWWbar[40]*armWWbar[166];
   armWWbar[166]=1./3.*armWWbar[168] + armWWbar[166];
   armWWbar[166]=armWWbar[18]*armWWbar[166];
   armWWbar[168]=1./3. + 3./2.*armWWbar[6];
   armWWbar[168]=armWWbar[1]*armWWbar[168];
   armWWbar[170]=1 + 9./2.*armWWbar[6];
   armWWbar[170]=armWWbar[40]*armWWbar[170];
   armWWbar[168]=armWWbar[168] + armWWbar[170];
   armWWbar[168]=armWWbar[19]*armWWbar[168];
   armWWbar[170]=2*armWWbar[6];
   armWWbar[174]=armWWbar[170] + 2./3. + armWWbar[7];
   armWWbar[174]=armWWbar[1]*armWWbar[174];
   armWWbar[175]=6*armWWbar[6] + 2 + 3*armWWbar[7];
   armWWbar[175]=armWWbar[40]*armWWbar[175];
   armWWbar[174]=armWWbar[174] + armWWbar[175];
   armWWbar[174]=armWWbar[20]*armWWbar[174];
   armWWbar[130]=2*armWWbar[23] + armWWbar[130] + 1./6.*armWWbar[21];
   armWWbar[130]=armWWbar[1]*armWWbar[130];
   armWWbar[125]=armWWbar[174] + armWWbar[168] + 1./2.*armWWbar[166] + 
   armWWbar[160] + armWWbar[130] + armWWbar[125];
   armWWbar[125]=armWWbar[3]*armWWbar[125];
   armWWbar[130]=armWWbar[208] + armWWbar[170];
   armWWbar[160]=armWWbar[1]*armWWbar[130];
   armWWbar[130]=armWWbar[40]*armWWbar[130];
   armWWbar[166]= - armWWbar[112]*armWWbar[7];
   armWWbar[168]=armWWbar[1]*armWWbar[166];
   armWWbar[166]=armWWbar[40]*armWWbar[166];
   armWWbar[166]=1./3.*armWWbar[168] + armWWbar[166];
   armWWbar[166]=MMZ*armWWbar[166];
   armWWbar[130]=armWWbar[166] + 1./3.*armWWbar[160] + armWWbar[130];
   armWWbar[130]=armWWbar[11]*armWWbar[130];
   armWWbar[160]= - armWWbar[25] - armWWbar[24];
   armWWbar[166]=armWWbar[160] - 1./3.*armWWbar[27];
   armWWbar[168]=1./3.*armWWbar[7] + armWWbar[212] - armWWbar[14] + 1./
   3.*armWWbar[28] - 1 - 2./3.*armWWbar[31];
   armWWbar[168]=armWWbar[112]*armWWbar[168];
   armWWbar[166]=2*armWWbar[166] + armWWbar[168];
   armWWbar[166]=armWWbar[1]*armWWbar[166];
   armWWbar[160]=3*armWWbar[160] - armWWbar[27];
   armWWbar[168]=armWWbar[7] - armWWbar[16] - 3*armWWbar[14] + 
   armWWbar[28] - 3 - 2*armWWbar[31];
   armWWbar[168]=armWWbar[112]*armWWbar[168];
   armWWbar[160]=2*armWWbar[160] + armWWbar[168];
   armWWbar[160]=armWWbar[40]*armWWbar[160];
   armWWbar[160]=armWWbar[166] + armWWbar[160];
   armWWbar[160]=MMZ*armWWbar[160];
   armWWbar[166]=armWWbar[112]*armWWbar[208];
   armWWbar[168]=armWWbar[1]*armWWbar[166];
   armWWbar[166]=armWWbar[40]*armWWbar[166];
   armWWbar[166]=1./3.*armWWbar[168] + armWWbar[166];
   armWWbar[166]=armWWbar[18]*armWWbar[166];
   armWWbar[168]=armWWbar[112]*armWWbar[7];
   armWWbar[170]=armWWbar[1]*armWWbar[168];
   armWWbar[168]=armWWbar[40]*armWWbar[168];
   armWWbar[168]=1./3.*armWWbar[170] + armWWbar[168];
   armWWbar[170]=armWWbar[20]*armWWbar[168];
   armWWbar[117]=armWWbar[119] + armWWbar[117] + armWWbar[123] + 
   armWWbar[125] + armWWbar[170] + armWWbar[130] + 1./2.*armWWbar[166]
    + armWWbar[160] + armWWbar[153] + armWWbar[120];
   armWWbar[117]=armWWbar[165]*armWWbar[117];
   armWWbar[119]=armWWbar[128] - 2*armWWbar[57];
   armWWbar[120]=5*armWWbar[56];
   armWWbar[119]=armWWbar[59] + 3*armWWbar[54] + 2./3.*armWWbar[119] + 
   armWWbar[120];
   armWWbar[123]= - armWWbar[71] + 307./4. + armWWbar[69];
   armWWbar[125]= - armWWbar[74] + 1./2.*armWWbar[8];
   armWWbar[128]=armWWbar[110]*armWWbar[125];
   armWWbar[123]=1./6.*armWWbar[123] + armWWbar[128];
   armWWbar[123]=armWWbar[110]*armWWbar[123];
   armWWbar[128]= - 1 + armWWbar[137];
   armWWbar[128]=1./2.*armWWbar[128] - armWWbar[71];
   armWWbar[128]=MMZ*armWWbar[136]*armWWbar[128];
   armWWbar[130]=5./6.*armWWbar[66] - 1./8. + armWWbar[67];
   armWWbar[130]= - 11./12.*armWWbar[50] + 1./4.*armWWbar[69] + 1./2.*
   armWWbar[130] + armWWbar[33];
   armWWbar[130]=armWWbar[35]*armWWbar[130];
   armWWbar[137]= - 3067 - 123*armWWbar[66];
   armWWbar[137]= - 165*armWWbar[17] + 123./8.*armWWbar[50] + 123./8.*
   armWWbar[71] - 123./4.*armWWbar[69] + 1./8.*armWWbar[137] + 165*
   armWWbar[63];
   armWWbar[137]=armWWbar[111]*armWWbar[137];
   armWWbar[153]=armWWbar[16] - 6*armWWbar[68] - 15./4. - armWWbar[60];
   armWWbar[153]=armWWbar[112]*armWWbar[153];
   armWWbar[119]=armWWbar[128] + armWWbar[153] + 495./64.*armWWbar[190]
    + 3./64.*armWWbar[137] + 1./4.*armWWbar[123] + 2*armWWbar[119] + 
   armWWbar[130];
   armWWbar[119]=MMZ*armWWbar[119];
   armWWbar[123]=8497./128. - armWWbar[47];
   armWWbar[128]=53./3.*armWWbar[13];
   armWWbar[130]=369./512.*armWWbar[111] - armWWbar[112];
   armWWbar[130]=MMZ*armWWbar[130];
   armWWbar[137]=armWWbar[311] - 369./64.*armWWbar[111];
   armWWbar[137]=armWWbar[19]*armWWbar[137];
   armWWbar[153]=armWWbar[11]*MMZ*armWWbar[112];
   armWWbar[123]=113./36.*armWWbar[37] + armWWbar[153] + 1./4.*
   armWWbar[137] + armWWbar[130] + 20./9.*armWWbar[223] + 4./3.*
   armWWbar[259] + 19./12.*armWWbar[12] + 1./4.*armWWbar[123] + 
   armWWbar[128];
   armWWbar[123]=armWWbar[37]*armWWbar[123];
   armWWbar[130]=armWWbar[147] + armWWbar[42] - armWWbar[8];
   armWWbar[137]= - 4*armWWbar[6] + 89./3. - 25*armWWbar[36];
   armWWbar[137]=MMZ*armWWbar[137];
   armWWbar[147]=armWWbar[151] + 11./3. - 29./2.*armWWbar[36];
   armWWbar[147]=1./6.*armWWbar[19]*armWWbar[147];
   armWWbar[160]= - armWWbar[37]*MMZ;
   armWWbar[166]=5./3. + armWWbar[198];
   armWWbar[170]=armWWbar[166] - armWWbar[6];
   armWWbar[174]=armWWbar[20]*armWWbar[170];
   armWWbar[175]= - 2*armWWbar[47];
   armWWbar[176]= - 47./3. + armWWbar[175];
   armWWbar[176]=armWWbar[38]*armWWbar[176];
   armWWbar[137]=armWWbar[176] + 4./3.*armWWbar[174] + 5*armWWbar[160]
    + armWWbar[147] + armWWbar[130] + 1./3.*armWWbar[137];
   armWWbar[137]=armWWbar[3]*armWWbar[137];
   armWWbar[160]= - 7./3. + armWWbar[115];
   armWWbar[174]=10./3. - armWWbar[7];
   armWWbar[174]=armWWbar[6]*armWWbar[174];
   armWWbar[166]=armWWbar[7]*armWWbar[166];
   armWWbar[160]=armWWbar[174] + 2./3.*armWWbar[160] + armWWbar[166];
   armWWbar[160]=armWWbar[1]*armWWbar[160];
   armWWbar[115]= - 5./3. + armWWbar[115];
   armWWbar[174]=2 - armWWbar[7];
   armWWbar[174]=armWWbar[6]*armWWbar[174];
   armWWbar[115]=armWWbar[174] + 2./3.*armWWbar[115] + armWWbar[166];
   armWWbar[115]=armWWbar[40]*armWWbar[115];
   armWWbar[166]= - 1./2.*armWWbar[35] + armWWbar[110];
   armWWbar[174]=215./128.*armWWbar[111];
   armWWbar[176]=MMZ*armWWbar[136];
   armWWbar[166]=1./2.*armWWbar[176] + 1./3.*armWWbar[166] + 
   armWWbar[174];
   armWWbar[166]=armWWbar[19]*armWWbar[166];
   armWWbar[179]=85./9.*armWWbar[35] + 369./32.*armWWbar[111];
   armWWbar[179]=armWWbar[37]*armWWbar[179];
   armWWbar[176]=1./2.*armWWbar[179] + 3./4.*armWWbar[176] + 495./8.*
   armWWbar[185] - 7387./64.*armWWbar[111] - 73./18.*armWWbar[35] + 
   armWWbar[110];
   armWWbar[176]=armWWbar[38]*armWWbar[176];
   armWWbar[179]= - 1./2.*armWWbar[37];
   armWWbar[180]= - 2 + armWWbar[169];
   armWWbar[180]=armWWbar[3]*armWWbar[38]*armWWbar[180];
   armWWbar[181]=6*armWWbar[180];
   armWWbar[183]=armWWbar[181] + armWWbar[179] - 77./3.*armWWbar[36] + 
   armWWbar[131] - 25./6. - armWWbar[47];
   armWWbar[183]=armWWbar[3]*armWWbar[183];
   armWWbar[184]=armWWbar[58] + 1./3.*armWWbar[52];
   armWWbar[187]=armWWbar[184] - 2./3.*armWWbar[57];
   armWWbar[198]=12*armWWbar[186];
   armWWbar[183]=armWWbar[198] + armWWbar[183] - 23./6.*armWWbar[54] + 
   2*armWWbar[187] + 3*armWWbar[56];
   armWWbar[183]=MMt*armWWbar[183];
   armWWbar[125]=1./3.*armWWbar[125] + armWWbar[298];
   armWWbar[125]=armWWbar[110]*armWWbar[125];
   armWWbar[187]=1./2.*armWWbar[125];
   armWWbar[199]=armWWbar[311] + armWWbar[322];
   armWWbar[199]=armWWbar[22]*armWWbar[199];
   armWWbar[202]=armWWbar[6] - 5./3. + 4*armWWbar[36];
   armWWbar[203]=armWWbar[37]*armWWbar[170];
   armWWbar[203]=1./3.*armWWbar[202] + armWWbar[203];
   armWWbar[203]=armWWbar[34]*armWWbar[203];
   armWWbar[208]=25./3.*armWWbar[53] - 722227./41472. - 8*armWWbar[64];
   armWWbar[138]= - armWWbar[74] + armWWbar[138];
   armWWbar[138]=1./3.*armWWbar[35]*armWWbar[138];
   armWWbar[209]=7./8.*armWWbar[42] + 11*armWWbar[72];
   armWWbar[209]= - 121*armWWbar[23] - 523./16.*armWWbar[22] + 11*
   armWWbar[209] - 1721./8.*armWWbar[51];
   armWWbar[209]=armWWbar[111]*armWWbar[209];
   armWWbar[210]= - 53./9.*armWWbar[13];
   armWWbar[212]=2809./48. + 235*armWWbar[36];
   armWWbar[212]=armWWbar[12]*armWWbar[212];
   armWWbar[216]= - 95./81. + 27*armWWbar[12];
   armWWbar[216]=armWWbar[6]*armWWbar[216];
   armWWbar[170]=armWWbar[11]*armWWbar[170];
   armWWbar[218]=armWWbar[20]*armWWbar[111];
   armWWbar[115]=4./9.*armWWbar[203] + armWWbar[183] + armWWbar[137] + 
   1./2.*armWWbar[176] + 121./32.*armWWbar[218] + armWWbar[123] + 4./9.
   *armWWbar[170] + 1./2.*armWWbar[166] + armWWbar[119] + 4./9.*
   armWWbar[115] + 4./27.*armWWbar[160] + 1./4.*armWWbar[216] + 1./12.*
   armWWbar[212] + 997./324.*armWWbar[36] + armWWbar[210] + 1./16.*
   armWWbar[209] + 1./4.*armWWbar[199] + armWWbar[187] + 17./64.*
   armWWbar[17] + 79./36.*armWWbar[211] + armWWbar[138] + 1./2.*
   armWWbar[47] + 659./1536.*armWWbar[50] - 143./512.*armWWbar[71] + 
   4045./768.*armWWbar[69] - 5./24.*armWWbar[30] + 19./36.*armWWbar[33]
    + 1093./192.*armWWbar[63] - 659./1536.*armWWbar[66] + 1./3.*
   armWWbar[208] - 23./4.*armWWbar[61];
   armWWbar[115]=armWWbar[34]*armWWbar[115];
   armWWbar[119]=19./4.*armWWbar[12];
   armWWbar[123]=armWWbar[119] + armWWbar[163] + 179./4.*armWWbar[13]
    + armWWbar[169] + 241./6. + armWWbar[172];
   armWWbar[123]=armWWbar[19]*armWWbar[123];
   armWWbar[137]= - 1 + 2*armWWbar[48];
   armWWbar[137]=2*armWWbar[137] - armWWbar[46];
   armWWbar[137]=armWWbar[159]*armWWbar[137];
   armWWbar[160]=armWWbar[277] - 1 + armWWbar[172];
   armWWbar[160]=armWWbar[20]*MMZ*armWWbar[160];
   armWWbar[166]=1 + armWWbar[277];
   armWWbar[170]=MMZ*armWWbar[166];
   armWWbar[176]=armWWbar[19]*armWWbar[170];
   armWWbar[137]=armWWbar[160] + armWWbar[137] + armWWbar[176];
   armWWbar[137]=armWWbar[3]*armWWbar[137];
   armWWbar[157]=205*armWWbar[13] + armWWbar[157] - 1975./18. + 
   armWWbar[48];
   armWWbar[157]=125./2.*armWWbar[12] + 1./2.*armWWbar[157] + 
   armWWbar[152];
   armWWbar[157]=MMZ*armWWbar[157];
   armWWbar[160]=armWWbar[164] - 17./2. + armWWbar[169];
   armWWbar[160]=armWWbar[18]*armWWbar[160];
   armWWbar[164]=armWWbar[173] - 23./2. + armWWbar[232];
   armWWbar[164]=MMZ*armWWbar[164];
   armWWbar[164]=armWWbar[164] + armWWbar[195];
   armWWbar[164]=armWWbar[11]*armWWbar[164];
   armWWbar[183]=53*armWWbar[13];
   armWWbar[195]=85./4.*armWWbar[12] + armWWbar[183] + 1151./12. + 
   armWWbar[169];
   armWWbar[195]=armWWbar[20]*armWWbar[195];
   armWWbar[199]= - 11*armWWbar[80] - 3./2.*armWWbar[78];
   armWWbar[199]=1./2.*armWWbar[199] + 7*armWWbar[22];
   armWWbar[123]=3*armWWbar[137] + armWWbar[195] + 1./2.*armWWbar[164]
    + 1./2.*armWWbar[123] + 1./4.*armWWbar[160] + armWWbar[157] + 33*
   armWWbar[23] + 3*armWWbar[199] + 217./24.*armWWbar[21];
   armWWbar[123]=armWWbar[3]*armWWbar[123];
   armWWbar[116]= - 55./4.*armWWbar[12] + armWWbar[321] + armWWbar[116]
    + 29./2. + armWWbar[288];
   armWWbar[116]=armWWbar[112]*armWWbar[116];
   armWWbar[116]=7./4.*armWWbar[113] + armWWbar[116];
   armWWbar[116]=MMZ*armWWbar[116];
   armWWbar[137]= - 2*armWWbar[46];
   armWWbar[157]=7./4.*armWWbar[178];
   armWWbar[116]=armWWbar[153] + armWWbar[157] + 4*armWWbar[122] + 
   armWWbar[116] + 11./12.*armWWbar[12] + armWWbar[128] + 157./4. + 
   armWWbar[137];
   armWWbar[116]=armWWbar[11]*armWWbar[116];
   armWWbar[122]= - 7*armWWbar[96] + armWWbar[162] + 7 + armWWbar[161];
   armWWbar[122]=armWWbar[142] + armWWbar[206] + 1./2.*armWWbar[122] + 
   armWWbar[204];
   armWWbar[122]=1./2.*armWWbar[113]*armWWbar[122];
   armWWbar[142]=armWWbar[231] - 40*armWWbar[16] + armWWbar[238] - 55./
   4.*armWWbar[17] + armWWbar[234] + 110*armWWbar[100] + armWWbar[197]
    + 11./2.*armWWbar[103] - 2*armWWbar[91] + 55./4.*armWWbar[92] + 55./
   4.*armWWbar[93] + 3*armWWbar[97] - 48*armWWbar[101] + 913./16. + 9*
   armWWbar[95];
   armWWbar[142]=armWWbar[112]*armWWbar[142];
   armWWbar[153]=2*armWWbar[90] + armWWbar[87];
   armWWbar[153]=11*armWWbar[153] - 357./16.*armWWbar[83];
   armWWbar[160]=11./4.*armWWbar[214];
   armWWbar[142]=armWWbar[142] + armWWbar[160] + armWWbar[122] + 55./2.
   *armWWbar[82] + 11*armWWbar[85] + 11./2.*armWWbar[109] - 
   armWWbar[81] + 3*armWWbar[153] - 8*armWWbar[86];
   armWWbar[142]=MMZ*armWWbar[142];
   armWWbar[153]= - 3./4.*armWWbar[107];
   armWWbar[161]= - 1./4.*armWWbar[78];
   armWWbar[162]=2*armWWbar[22];
   armWWbar[164]=107./12.*armWWbar[23] - 35./3.*armWWbar[21] + 
   armWWbar[162] + armWWbar[161] + armWWbar[153] + 11*armWWbar[240] + 
   47./12.*armWWbar[79];
   armWWbar[164]=armWWbar[112]*armWWbar[164];
   armWWbar[118]=armWWbar[118] + 5./3.*armWWbar[86];
   armWWbar[118]= - 14./3.*armWWbar[82] + armWWbar[141] + 4*
   armWWbar[118] + armWWbar[139];
   armWWbar[118]=MMH*armWWbar[118];
   armWWbar[139]=61./3.*armWWbar[113] + armWWbar[154];
   armWWbar[139]=1./2.*armWWbar[139];
   armWWbar[141]= - 293./3. + armWWbar[46];
   armWWbar[141]=armWWbar[112]*armWWbar[141];
   armWWbar[141]=armWWbar[139] + armWWbar[141];
   armWWbar[141]=armWWbar[18]*armWWbar[141];
   armWWbar[154]=15./4.*armWWbar[113];
   armWWbar[195]= - 341./12. - armWWbar[46];
   armWWbar[195]=armWWbar[112]*armWWbar[195];
   armWWbar[195]=armWWbar[154] + armWWbar[195];
   armWWbar[195]=armWWbar[20]*armWWbar[195];
   armWWbar[197]=168871./288. + 57*armWWbar[104];
   armWWbar[199]= - 15./2.*armWWbar[23] - 5./12.*armWWbar[21] - 29./6.*
   armWWbar[22] + 15./2.*armWWbar[107] - 7./3.*armWWbar[78];
   armWWbar[199]=armWWbar[113]*armWWbar[199];
   armWWbar[203]= - 12835./27. - 363*armWWbar[13];
   armWWbar[203]=1./4.*armWWbar[203] + armWWbar[10];
   armWWbar[203]=5*armWWbar[203] - 511./8.*armWWbar[12];
   armWWbar[203]=armWWbar[12]*armWWbar[203];
   armWWbar[204]= - 119./3.*armWWbar[113] + armWWbar[233];
   armWWbar[206]=armWWbar[204] - 33*armWWbar[112];
   armWWbar[206]=armWWbar[19]*armWWbar[206];
   armWWbar[116]=armWWbar[118] + armWWbar[123] + 1./2.*armWWbar[195] + 
   armWWbar[116] + 1./4.*armWWbar[206] + 1./4.*armWWbar[141] + 
   armWWbar[142] + armWWbar[164] + 1./2.*armWWbar[203] - 7./2.*
   armWWbar[10] - 2707./216.*armWWbar[13] + 1./2.*armWWbar[199] + 7./4.
   *armWWbar[46] - 121./12.*armWWbar[16] + armWWbar[318] + 3709./48.*
   armWWbar[17] + armWWbar[167] + 55./4.*armWWbar[100] + 19./6.*
   armWWbar[96] + 51./4.*armWWbar[103] - 11./12.*armWWbar[98] - 11./12.
   *armWWbar[92] - 11./12.*armWWbar[93] + armWWbar[217] - 1285./16.*
   armWWbar[99] + 215./8.*armWWbar[94] + 1./2.*armWWbar[197] + 8./3.*
   armWWbar[102];
   armWWbar[116]=armWWbar[76]*armWWbar[116];
   armWWbar[118]= - 23*armWWbar[29];
   armWWbar[123]= - 37./2.*armWWbar[15];
   armWWbar[141]=7*armWWbar[32];
   armWWbar[142]= - 7*armWWbar[134];
   armWWbar[164]=armWWbar[142] + armWWbar[141] + armWWbar[123] - 35171./
   972. + armWWbar[118];
   armWWbar[167]=193./4. + armWWbar[183];
   armWWbar[183]=7./4.*armWWbar[7];
   armWWbar[167]=1./3.*armWWbar[167] + armWWbar[183];
   armWWbar[167]=armWWbar[7]*armWWbar[167];
   armWWbar[121]=1./3.*armWWbar[220] + armWWbar[121];
   armWWbar[195]=armWWbar[1]*armWWbar[121];
   armWWbar[121]=armWWbar[40]*armWWbar[121];
   armWWbar[197]=10./3.*armWWbar[33];
   armWWbar[203]=19./2.*armWWbar[7];
   armWWbar[206]=227./3. + armWWbar[203];
   armWWbar[206]=armWWbar[12]*armWWbar[206];
   armWWbar[208]=451./54. + 79*armWWbar[12];
   armWWbar[208]=armWWbar[6]*armWWbar[208];
   armWWbar[121]=20./9.*armWWbar[121] + 56./27.*armWWbar[195] + 1./3.*
   armWWbar[208] + 1./6.*armWWbar[206] + armWWbar[167] + armWWbar[210]
    + 8*armWWbar[17] - 9./4.*armWWbar[30] + 1./4.*armWWbar[164] + 
   armWWbar[197];
   armWWbar[121]=armWWbar[40]*armWWbar[121];
   armWWbar[118]=armWWbar[142] + armWWbar[141] + armWWbar[123] - 4003./
   108. + armWWbar[118];
   armWWbar[118]=1./4.*armWWbar[118] + armWWbar[197];
   armWWbar[123]=armWWbar[183] + 67./4. + armWWbar[128];
   armWWbar[123]=armWWbar[7]*armWWbar[123];
   armWWbar[128]=8./3.*armWWbar[17];
   armWWbar[141]=139./3. + armWWbar[203];
   armWWbar[141]=armWWbar[12]*armWWbar[141];
   armWWbar[142]=155./18. + 41*armWWbar[12];
   armWWbar[142]=armWWbar[6]*armWWbar[142];
   armWWbar[118]=4./9.*armWWbar[195] + 1./3.*armWWbar[142] + 1./18.*
   armWWbar[141] + 1./3.*armWWbar[123] - 53./27.*armWWbar[13] + 
   armWWbar[128] + 1./3.*armWWbar[118] - 3./4.*armWWbar[30];
   armWWbar[118]=armWWbar[1]*armWWbar[118];
   armWWbar[123]=7./3. - armWWbar[7];
   armWWbar[123]=1./3.*armWWbar[123] - armWWbar[6];
   armWWbar[123]=armWWbar[1]*armWWbar[123];
   armWWbar[141]= - 29./3.*armWWbar[6] + 89./9. - 5*armWWbar[7];
   armWWbar[141]=armWWbar[40]*armWWbar[141];
   armWWbar[123]=5*armWWbar[123] + armWWbar[141];
   armWWbar[123]=MMZ*armWWbar[123];
   armWWbar[140]=armWWbar[1]*armWWbar[140];
   armWWbar[140]=armWWbar[140] + 11./3.*armWWbar[223];
   armWWbar[140]=1./2.*armWWbar[19]*armWWbar[140];
   armWWbar[141]=armWWbar[259] + 5./3.*armWWbar[223];
   armWWbar[142]=armWWbar[20]*armWWbar[141];
   armWWbar[123]=4*armWWbar[142] + armWWbar[123] + armWWbar[140];
   armWWbar[123]=armWWbar[3]*armWWbar[123];
   armWWbar[142]=2*armWWbar[25];
   armWWbar[164]=2*armWWbar[24];
   armWWbar[167]=armWWbar[164] - armWWbar[26] + armWWbar[142];
   armWWbar[167]=2*armWWbar[167] + armWWbar[27];
   armWWbar[155]= - 1./3.*armWWbar[7] + armWWbar[158] + armWWbar[14] + 
   armWWbar[155] + 1 + 2./3.*armWWbar[31];
   armWWbar[155]=armWWbar[112]*armWWbar[155];
   armWWbar[155]=2./3.*armWWbar[167] + armWWbar[155];
   armWWbar[155]=armWWbar[1]*armWWbar[155];
   armWWbar[142]=armWWbar[164] + armWWbar[142] - 1./9.*armWWbar[75] - 
   armWWbar[26];
   armWWbar[142]=2*armWWbar[142] + armWWbar[27];
   armWWbar[158]= - armWWbar[7] + armWWbar[16] + 3*armWWbar[14] - 
   armWWbar[28] + 3 + 2*armWWbar[31];
   armWWbar[158]=armWWbar[112]*armWWbar[158];
   armWWbar[142]=2*armWWbar[142] + armWWbar[158];
   armWWbar[142]=armWWbar[40]*armWWbar[142];
   armWWbar[142]=armWWbar[155] + armWWbar[142];
   armWWbar[142]=MMZ*armWWbar[142];
   armWWbar[155]=MMZ*armWWbar[168];
   armWWbar[141]=4./3.*armWWbar[141] + armWWbar[155];
   armWWbar[141]=armWWbar[11]*armWWbar[141];
   armWWbar[115]=armWWbar[117] + armWWbar[115] + armWWbar[116] + 
   armWWbar[123] + armWWbar[141] + armWWbar[142] + armWWbar[118] + 
   armWWbar[121];
   armWWbar[115]=armWWbar[165]*armWWbar[115];
   armWWbar[116]=121./3. - 5*armWWbar[66];
   armWWbar[116]=5./8.*armWWbar[116] - 37*armWWbar[63];
   armWWbar[117]= - 1./8.*armWWbar[42] - armWWbar[72];
   armWWbar[117]=1./3.*armWWbar[23] - 1./48.*armWWbar[22] + 1./3.*
   armWWbar[117] + 7./8.*armWWbar[51];
   armWWbar[117]=armWWbar[111]*armWWbar[117];
   armWWbar[118]=101./24. + armWWbar[329];
   armWWbar[118]=armWWbar[12]*armWWbar[118];
   armWWbar[121]=armWWbar[6]*armWWbar[12];
   armWWbar[116]=5./3.*armWWbar[121] + 1./3.*armWWbar[118] + 25./2.*
   armWWbar[117] + 151./24.*armWWbar[17] + 25./64.*armWWbar[50] + 25./
   64.*armWWbar[71] + 565./96.*armWWbar[69] + 1./8.*armWWbar[116] - 5./
   3.*armWWbar[30];
   armWWbar[117]= - 3 - 1./3.*armWWbar[66];
   armWWbar[117]= - 1./9.*armWWbar[17] + 1./24.*armWWbar[50] + 1./24.*
   armWWbar[71] - 1./12.*armWWbar[69] + 1./8.*armWWbar[117] + 1./9.*
   armWWbar[63];
   armWWbar[117]=armWWbar[111]*armWWbar[117];
   armWWbar[117]=armWWbar[117] + 1./9.*armWWbar[190];
   armWWbar[117]=MMZ*armWWbar[117];
   armWWbar[118]=MMZ*armWWbar[111];
   armWWbar[123]=1 + armWWbar[118];
   armWWbar[123]=1./2.*armWWbar[123] + armWWbar[241];
   armWWbar[123]=armWWbar[37]*armWWbar[123];
   armWWbar[141]=armWWbar[37]*armWWbar[111];
   armWWbar[142]=1./8.*armWWbar[141] + 5./8.*armWWbar[111] + 1./3.*
   armWWbar[185];
   armWWbar[142]=armWWbar[38]*armWWbar[142];
   armWWbar[155]=armWWbar[19]*armWWbar[111];
   armWWbar[116]=25./6.*armWWbar[142] + 25./36.*armWWbar[327] + 25./96.
   *armWWbar[123] + 125./288.*armWWbar[155] + 1./3.*armWWbar[116] + 25./
   8.*armWWbar[117];
   armWWbar[116]=armWWbar[255]*armWWbar[116];
   armWWbar[117]=armWWbar[330] + 11 + armWWbar[329];
   armWWbar[117]=armWWbar[19]*armWWbar[117];
   armWWbar[123]= - 17*armWWbar[42] - 5*armWWbar[8];
   armWWbar[142]=35./4.*armWWbar[6] - 121./3. + 119./4.*armWWbar[36];
   armWWbar[142]=MMZ*armWWbar[142];
   armWWbar[155]=7./6. + armWWbar[47];
   armWWbar[155]=armWWbar[38]*armWWbar[155];
   armWWbar[117]=17*armWWbar[155] + 1./2.*armWWbar[117] + 1./3.*
   armWWbar[142] + armWWbar[191] + 1./2.*armWWbar[123] + 17*
   armWWbar[51];
   armWWbar[117]=armWWbar[3]*armWWbar[117];
   armWWbar[123]=2545./3. - 163*armWWbar[66];
   armWWbar[123]=1./8.*armWWbar[123] - 1733./9.*armWWbar[63];
   armWWbar[123]=1./4.*armWWbar[123] - 5./3.*armWWbar[33];
   armWWbar[142]=armWWbar[35]*armWWbar[205];
   armWWbar[123]=1./4.*armWWbar[125] + 2437./288.*armWWbar[17] + 11./4.
   *armWWbar[211] + 11./2.*armWWbar[142] - 17./4.*armWWbar[47] + 163./
   256.*armWWbar[50] + 2057./768.*armWWbar[71] + 151./384.*armWWbar[69]
    + 1./8.*armWWbar[123] - 22./9.*armWWbar[30];
   armWWbar[125]=1 + armWWbar[69];
   armWWbar[142]=armWWbar[35]*armWWbar[125];
   armWWbar[155]= - 1 - armWWbar[69];
   armWWbar[155]=armWWbar[110]*armWWbar[155];
   armWWbar[142]=1./12.*armWWbar[155] + armWWbar[246] + 11./2.*
   armWWbar[142];
   armWWbar[155]= - 31 - 7./3.*armWWbar[66];
   armWWbar[155]= - 37./9.*armWWbar[17] + 7./24.*armWWbar[50] + 7./24.*
   armWWbar[71] - 7./12.*armWWbar[69] + 1./8.*armWWbar[155] + 37./9.*
   armWWbar[63];
   armWWbar[155]=armWWbar[111]*armWWbar[155];
   armWWbar[142]=185./144.*armWWbar[190] + 1./3.*armWWbar[142] + 5./16.
   *armWWbar[155];
   armWWbar[142]=MMZ*armWWbar[142];
   armWWbar[118]=35./16.*armWWbar[241] + 35./32.*armWWbar[118] + 35./32.
    + armWWbar[267];
   armWWbar[155]= - 1./3.*armWWbar[37];
   armWWbar[118]=1./4.*armWWbar[118] + armWWbar[155];
   armWWbar[118]=armWWbar[37]*armWWbar[118];
   armWWbar[141]=35./32.*armWWbar[141] + 185./12.*armWWbar[185] + 
   armWWbar[331] + 3025./96.*armWWbar[111];
   armWWbar[141]=armWWbar[38]*armWWbar[141];
   armWWbar[158]= - 35./6. + armWWbar[267];
   armWWbar[158]=49./18.*armWWbar[36] + 1./3.*armWWbar[158] + 
   armWWbar[131];
   armWWbar[158]=1./2.*armWWbar[158] + 3*armWWbar[180];
   armWWbar[158]=armWWbar[3]*armWWbar[158];
   armWWbar[164]= - 2*armWWbar[54] + armWWbar[246];
   armWWbar[158]=6*armWWbar[186] + 1./3.*armWWbar[164] + armWWbar[158];
   armWWbar[158]=MMt*armWWbar[158];
   armWWbar[164]=103./8.*armWWbar[42] - 47*armWWbar[72];
   armWWbar[164]=47./9.*armWWbar[23] - 227./144.*armWWbar[22] + 1./9.*
   armWWbar[164] + 63./8.*armWWbar[51];
   armWWbar[164]=armWWbar[111]*armWWbar[164];
   armWWbar[165]=2863./36. + 181*armWWbar[36];
   armWWbar[165]=armWWbar[12]*armWWbar[165];
   armWWbar[167]= - 5 + 61./4.*armWWbar[12];
   armWWbar[167]=armWWbar[6]*armWWbar[167];
   armWWbar[168]= - 925./192.*armWWbar[111] + armWWbar[276] + 1./6.*
   armWWbar[110];
   armWWbar[168]=armWWbar[19]*armWWbar[168];
   armWWbar[116]=1./8.*armWWbar[116] + armWWbar[158] + 1./3.*
   armWWbar[117] + 1./6.*armWWbar[141] + 235./144.*armWWbar[327] + 1./6.
   *armWWbar[118] + 1./6.*armWWbar[168] + 1./2.*armWWbar[142] + 1./18.*
   armWWbar[167] + 1./72.*armWWbar[165] - 17./18.*armWWbar[36] + 1./3.*
   armWWbar[123] + 5./8.*armWWbar[164];
   armWWbar[116]=armWWbar[255]*armWWbar[116];
   armWWbar[117]=11*armWWbar[50] - 13*armWWbar[69] - 19./4. - 11*
   armWWbar[66];
   armWWbar[117]=armWWbar[35]*armWWbar[117];
   armWWbar[117]=armWWbar[246] + 1./4.*armWWbar[117];
   armWWbar[118]=1./9.*armWWbar[71] - 5./12. - armWWbar[69];
   armWWbar[118]=1./2.*armWWbar[118] + 1./3.*armWWbar[126];
   armWWbar[118]=armWWbar[110]*armWWbar[118];
   armWWbar[123]=MMZ*armWWbar[136]*armWWbar[125];
   armWWbar[125]=3343 + 941./3.*armWWbar[66];
   armWWbar[125]=2201./9.*armWWbar[17] - 941./24.*armWWbar[50] - 941./
   24.*armWWbar[71] + 941./12.*armWWbar[69] + 1./8.*armWWbar[125] - 
   2201./9.*armWWbar[63];
   armWWbar[125]=armWWbar[111]*armWWbar[125];
   armWWbar[117]=1./24.*armWWbar[123] + 2201./576.*armWWbar[185] + 1./
   64.*armWWbar[125] + 1./3.*armWWbar[117] + 1./4.*armWWbar[118];
   armWWbar[117]=MMZ*armWWbar[117];
   armWWbar[118]= - MMZ*armWWbar[111];
   armWWbar[123]=armWWbar[276] + 941./64.*armWWbar[111];
   armWWbar[123]=armWWbar[19]*armWWbar[123];
   armWWbar[118]=armWWbar[155] + 1./3.*armWWbar[123] + 941./384.*
   armWWbar[118] + armWWbar[201] - 839./1152. - armWWbar[47];
   armWWbar[118]=armWWbar[37]*armWWbar[118];
   armWWbar[123]=8*armWWbar[6] + 11./3. + armWWbar[146];
   armWWbar[123]=MMZ*armWWbar[123];
   armWWbar[125]= - 29./9. + armWWbar[175];
   armWWbar[125]=armWWbar[38]*armWWbar[125];
   armWWbar[123]=armWWbar[125] + armWWbar[147] + armWWbar[130] + 1./9.*
   armWWbar[123];
   armWWbar[123]=armWWbar[3]*armWWbar[123];
   armWWbar[125]=armWWbar[333] - 941./32.*armWWbar[111];
   armWWbar[125]=armWWbar[37]*armWWbar[125];
   armWWbar[125]=1./6.*armWWbar[125] + 1./4.*armWWbar[193] + 2201./72.*
   armWWbar[190] + 17485./576.*armWWbar[111] - 7./6.*armWWbar[35] + 
   armWWbar[110];
   armWWbar[125]=armWWbar[38]*armWWbar[125];
   armWWbar[126]=armWWbar[181] + armWWbar[179] - 119./9.*armWWbar[36]
    + armWWbar[131] + 149./18. - armWWbar[47];
   armWWbar[126]=armWWbar[3]*armWWbar[126];
   armWWbar[120]=armWWbar[171] + armWWbar[120] - 11./2.*armWWbar[54];
   armWWbar[120]=armWWbar[198] + 1./3.*armWWbar[120] + armWWbar[126];
   armWWbar[120]=MMt*armWWbar[120];
   armWWbar[126]=1./6.*armWWbar[193] + armWWbar[174] + 5./2.*
   armWWbar[35] + 1./3.*armWWbar[110];
   armWWbar[126]=armWWbar[19]*armWWbar[126];
   armWWbar[130]= - 71./1152.*armWWbar[50] - 143./128.*armWWbar[71] + 
   4711./576.*armWWbar[69] + 179./54.*armWWbar[30] - 103./27.*
   armWWbar[33] - 2963./432.*armWWbar[63] + 71./1152.*armWWbar[66] + 
   257647./10368. - 3*armWWbar[61];
   armWWbar[130]=1./2.*armWWbar[130] + armWWbar[47];
   armWWbar[131]=armWWbar[276] - 1./2.*armWWbar[110];
   armWWbar[131]=armWWbar[22]*armWWbar[131];
   armWWbar[141]= - 107*armWWbar[23] + 2527./16.*armWWbar[22] - 1067./8.
   *armWWbar[51] - 793./8.*armWWbar[42] + 107*armWWbar[72];
   armWWbar[141]=armWWbar[111]*armWWbar[141];
   armWWbar[142]= - 5183./432. - 53*armWWbar[36];
   armWWbar[142]=armWWbar[12]*armWWbar[142];
   armWWbar[146]= - 7./3. + armWWbar[12];
   armWWbar[146]=armWWbar[6]*armWWbar[146];
   armWWbar[116]=armWWbar[116] + armWWbar[120] + armWWbar[123] + 1./2.*
   armWWbar[125] + 107./96.*armWWbar[218] + 1./4.*armWWbar[118] + 1./2.
   *armWWbar[126] + armWWbar[117] + 1./12.*armWWbar[146] + 1./12.*
   armWWbar[142] + 29./36.*armWWbar[36] + 1./48.*armWWbar[141] + 1./12.
   *armWWbar[131] + armWWbar[187] + 2827./1728.*armWWbar[17] + 3./4.*
   armWWbar[211] + 1./2.*armWWbar[130] + armWWbar[138];
   armWWbar[116]=armWWbar[255]*armWWbar[116];
   armWWbar[117]= - 5./3.*armWWbar[63] + 13./9. + armWWbar[64];
   armWWbar[117]=11./3.*armWWbar[17] + armWWbar[229] + 2*armWWbar[117]
    + 1./3.*armWWbar[33];
   armWWbar[118]= - 1./3. + armWWbar[144];
   armWWbar[118]=armWWbar[12]*armWWbar[118];
   armWWbar[120]=armWWbar[6]*armWWbar[249];
   armWWbar[117]=2./3.*armWWbar[120] + 4./3.*armWWbar[118] + 
   armWWbar[188] + 2*armWWbar[117] + armWWbar[13];
   armWWbar[118]= - 1./3.*armWWbar[59] - armWWbar[54] - 2*armWWbar[56]
    + armWWbar[184] + armWWbar[129];
   armWWbar[123]=armWWbar[128] + 7 - 8./3.*armWWbar[63];
   armWWbar[123]=armWWbar[111]*armWWbar[123];
   armWWbar[118]=8./3.*armWWbar[185] + 2./3.*armWWbar[118] + 
   armWWbar[123];
   armWWbar[123]=1./8. + 1./3.*armWWbar[71];
   armWWbar[123]=MMZ*armWWbar[136]*armWWbar[123];
   armWWbar[118]=2*armWWbar[118] + armWWbar[123];
   armWWbar[118]=MMZ*armWWbar[118];
   armWWbar[123]= - 1 - armWWbar[13];
   armWWbar[125]=armWWbar[37]*armWWbar[123];
   armWWbar[126]=armWWbar[111] + armWWbar[190];
   armWWbar[126]=armWWbar[38]*armWWbar[126];
   armWWbar[128]=MMZ*armWWbar[202];
   armWWbar[128]=armWWbar[128] + 8*armWWbar[38];
   armWWbar[128]=armWWbar[3]*armWWbar[128];
   armWWbar[129]=MMt*armWWbar[3]*armWWbar[143];
   armWWbar[116]=armWWbar[116] + 64./9.*armWWbar[129] + 8./9.*
   armWWbar[128] + 64./3.*armWWbar[126] + 4*armWWbar[125] + 4./3.*
   armWWbar[117] + armWWbar[118];
   armWWbar[116]=armWWbar[34]*armWWbar[116];
   armWWbar[117]=armWWbar[145] + armWWbar[127] - 23./8.*armWWbar[83] - 
   2*armWWbar[90] + armWWbar[87];
   armWWbar[118]= - 3*armWWbar[103] + 31./4. + armWWbar[135];
   armWWbar[118]=armWWbar[16] - 3./2.*armWWbar[17] + 3*armWWbar[100] + 
   1./2.*armWWbar[118] - armWWbar[96];
   armWWbar[118]=armWWbar[113]*armWWbar[118];
   armWWbar[125]= - 2 + 1./4.*armWWbar[103];
   armWWbar[125]=armWWbar[112]*armWWbar[125];
   armWWbar[117]=1./3.*armWWbar[125] + 3./2.*armWWbar[149] + 1./3.*
   armWWbar[117] + armWWbar[118];
   armWWbar[117]=MMZ*armWWbar[117];
   armWWbar[118]=armWWbar[177] + armWWbar[182] + 7./12.*armWWbar[13] - 
   143./36. + armWWbar[46];
   armWWbar[118]=MMZ*armWWbar[118];
   armWWbar[125]= - armWWbar[80] - 9./2.*armWWbar[78];
   armWWbar[126]=5./3.*armWWbar[12] + 1./3. + armWWbar[239];
   armWWbar[126]=armWWbar[18]*armWWbar[126];
   armWWbar[118]=armWWbar[126] + armWWbar[118] + armWWbar[23] + 47./12.
   *armWWbar[21] + 1./2.*armWWbar[125] + 5*armWWbar[22];
   armWWbar[125]=armWWbar[13] + 25./6. + armWWbar[169];
   armWWbar[124]=1./8.*armWWbar[12] + 1./4.*armWWbar[125] + 
   armWWbar[124];
   armWWbar[124]=armWWbar[19]*armWWbar[124];
   armWWbar[125]=armWWbar[159]*armWWbar[166];
   armWWbar[125]=armWWbar[125] + 3./2.*armWWbar[176];
   armWWbar[125]=armWWbar[3]*armWWbar[125];
   armWWbar[126]=1 + armWWbar[46];
   armWWbar[127]=3./2.*armWWbar[126] + armWWbar[12];
   armWWbar[127]=armWWbar[20]*armWWbar[127];
   armWWbar[118]=armWWbar[125] + 1./2.*armWWbar[127] + 1./2.*
   armWWbar[118] + armWWbar[124];
   armWWbar[118]=armWWbar[3]*armWWbar[118];
   armWWbar[124]=armWWbar[324] + armWWbar[196];
   armWWbar[125]=MMZ*armWWbar[113];
   armWWbar[124]=armWWbar[178] + 1./4.*armWWbar[124] + armWWbar[125];
   armWWbar[124]=armWWbar[11]*armWWbar[124];
   armWWbar[127]= - 1 + armWWbar[173];
   armWWbar[127]=armWWbar[112]*armWWbar[127];
   armWWbar[127]=armWWbar[127] - 13./3.*armWWbar[113] + 3*armWWbar[214]
   ;
   armWWbar[127]=armWWbar[18]*armWWbar[127];
   armWWbar[128]= - 1./36.*armWWbar[99] + 33./8.*armWWbar[94] - 10477./
   288. - 7*armWWbar[104];
   armWWbar[129]= - 2*armWWbar[23] + 5./12.*armWWbar[21] - 37./6.*
   armWWbar[22] + 2*armWWbar[107] + 1./3.*armWWbar[78];
   armWWbar[129]=armWWbar[113]*armWWbar[129];
   armWWbar[130]= - 103./9. + 29./4.*armWWbar[13];
   armWWbar[130]= - 7./12.*armWWbar[12] + 1./3.*armWWbar[130] + 13*
   armWWbar[10];
   armWWbar[130]=armWWbar[12]*armWWbar[130];
   armWWbar[131]=3./4.*armWWbar[21] + 3./2.*armWWbar[22] + 1./3.*
   armWWbar[240] - 3./4.*armWWbar[78];
   armWWbar[131]=armWWbar[112]*armWWbar[131];
   armWWbar[135]=7./6.*armWWbar[112] - 25./3.*armWWbar[113] + 3*
   armWWbar[149];
   armWWbar[135]=armWWbar[19]*armWWbar[135];
   armWWbar[136]= - 5./3. + armWWbar[277];
   armWWbar[136]=armWWbar[112]*armWWbar[136];
   armWWbar[136]=armWWbar[113] + 1./2.*armWWbar[136];
   armWWbar[136]=armWWbar[20]*armWWbar[136];
   armWWbar[138]=1./2.*armWWbar[189] + 1./3.*armWWbar[82];
   armWWbar[138]=MMH*armWWbar[138];
   armWWbar[117]=1./2.*armWWbar[138] + armWWbar[118] + armWWbar[136] + 
   armWWbar[124] + 1./2.*armWWbar[135] + 1./4.*armWWbar[127] + 
   armWWbar[117] + 1./2.*armWWbar[131] + 1./4.*armWWbar[130] - 13./6.*
   armWWbar[10] + armWWbar[225] + armWWbar[129] + armWWbar[224] - 5./12.
   *armWWbar[16] - 1./4.*armWWbar[48] - 151./144.*armWWbar[17] + 
   armWWbar[286] + 7./4.*armWWbar[100] + 5./12.*armWWbar[96] + 79./36.*
   armWWbar[103] + 1./2.*armWWbar[128] - 5./3.*armWWbar[98];
   armWWbar[117]=armWWbar[76]*armWWbar[117];
   armWWbar[118]=5./6.*armWWbar[103] - 1./2.*armWWbar[98] - 13./8.*
   armWWbar[99] + 7./2.*armWWbar[94] - 1489./96. - armWWbar[104];
   armWWbar[124]=1./6.*armWWbar[100];
   armWWbar[118]=armWWbar[253] - 11./24.*armWWbar[17] + armWWbar[124]
    + 1./3.*armWWbar[118] - 1./2.*armWWbar[96];
   armWWbar[127]= - armWWbar[103] + 23 + armWWbar[98];
   armWWbar[128]= - 1./12.*armWWbar[17];
   armWWbar[124]=armWWbar[16] + armWWbar[128] + armWWbar[124] + 1./12.*
   armWWbar[127] - armWWbar[96];
   armWWbar[124]=armWWbar[113]*armWWbar[124];
   armWWbar[124]=1./12.*armWWbar[149] + 3./16.*armWWbar[83] + 
   armWWbar[124];
   armWWbar[124]=MMZ*armWWbar[124];
   armWWbar[127]= - 7./3. + 1./4.*armWWbar[13];
   armWWbar[127]=armWWbar[177] + 1./3.*armWWbar[127] + armWWbar[10];
   armWWbar[127]=armWWbar[12]*armWWbar[127];
   armWWbar[125]=armWWbar[178] + 1./4. + armWWbar[125];
   armWWbar[125]=armWWbar[11]*armWWbar[125];
   armWWbar[129]= - armWWbar[23] + 5./24.*armWWbar[21] - 7./4.*
   armWWbar[22] + armWWbar[107] - 1./6.*armWWbar[78];
   armWWbar[129]=armWWbar[113]*armWWbar[129];
   armWWbar[130]=armWWbar[113] + 1./3.*armWWbar[214];
   armWWbar[130]=armWWbar[18]*armWWbar[130];
   armWWbar[131]= - 19*armWWbar[113] + armWWbar[149];
   armWWbar[131]=armWWbar[19]*armWWbar[131];
   armWWbar[135]=armWWbar[20]*armWWbar[113];
   armWWbar[136]= - 1 + armWWbar[12];
   armWWbar[136]=armWWbar[18]*armWWbar[136];
   armWWbar[136]=armWWbar[21] + armWWbar[136];
   armWWbar[136]=armWWbar[3]*armWWbar[136];
   armWWbar[138]=MMH*armWWbar[189];
   armWWbar[118]=1./12.*armWWbar[138] + 1./24.*armWWbar[136] + 1./2.*
   armWWbar[135] + armWWbar[125] + 1./12.*armWWbar[131] + 1./8.*
   armWWbar[130] + armWWbar[124] + 1./4.*armWWbar[127] + 1./2.*
   armWWbar[118] + armWWbar[129];
   armWWbar[118]=armWWbar[76]*armWWbar[118];
   armWWbar[124]= - armWWbar[32] - 1 - armWWbar[15];
   armWWbar[125]=1./2.*armWWbar[121] + armWWbar[177] + armWWbar[248] + 
   armWWbar[124] - 1./2.*armWWbar[30];
   armWWbar[125]=armWWbar[1]*armWWbar[125];
   armWWbar[121]=11./6.*armWWbar[121] + 11./18.*armWWbar[12] + 11./6.*
   armWWbar[17] + armWWbar[124] - 11./6.*armWWbar[30];
   armWWbar[121]=armWWbar[40]*armWWbar[121];
   armWWbar[121]=armWWbar[125] + 1./3.*armWWbar[121];
   armWWbar[124]=armWWbar[269] - 1 + 1./6.*armWWbar[94];
   armWWbar[124]=armWWbar[255]*armWWbar[76]*armWWbar[124];
   armWWbar[118]=1./8.*armWWbar[124] + 1./2.*armWWbar[121] + 
   armWWbar[118];
   armWWbar[118]=armWWbar[255]*armWWbar[118];
   armWWbar[121]=39./8. - 1./9.*armWWbar[15];
   armWWbar[124]= - 1 + armWWbar[132];
   armWWbar[124]=armWWbar[6]*armWWbar[124];
   armWWbar[121]=11./9.*armWWbar[124] + 283./162.*armWWbar[12] - 1./18.
   *armWWbar[148] + 187./54.*armWWbar[17] - 187./54.*armWWbar[30] - 11./
   36.*armWWbar[33] + 1./18.*armWWbar[134] + 1./4.*armWWbar[121] - 2./9.
   *armWWbar[32];
   armWWbar[121]=armWWbar[40]*armWWbar[121];
   armWWbar[125]= - 11./3. + armWWbar[151];
   armWWbar[127]=armWWbar[1]*armWWbar[125];
   armWWbar[125]=armWWbar[40]*armWWbar[125];
   armWWbar[125]=armWWbar[127] + 11./9.*armWWbar[125];
   armWWbar[125]=MMZ*armWWbar[125];
   armWWbar[127]=1./2. + armWWbar[6];
   armWWbar[129]=armWWbar[1]*armWWbar[127];
   armWWbar[127]=armWWbar[40]*armWWbar[127];
   armWWbar[127]=3*armWWbar[129] + 11./3.*armWWbar[127];
   armWWbar[127]=armWWbar[19]*armWWbar[127];
   armWWbar[129]= - armWWbar[8] + armWWbar[22];
   armWWbar[130]=armWWbar[1]*armWWbar[129];
   armWWbar[129]=armWWbar[40]*armWWbar[129];
   armWWbar[125]=armWWbar[127] + armWWbar[125] + 3*armWWbar[130] + 11./
   3.*armWWbar[129];
   armWWbar[125]=armWWbar[3]*armWWbar[125];
   armWWbar[127]=119./8. + 5./3.*armWWbar[15];
   armWWbar[127]=1./2.*armWWbar[134] + 1./4.*armWWbar[127] - 2*
   armWWbar[32];
   armWWbar[124]=armWWbar[124] + 25./18.*armWWbar[12] - 1./6.*
   armWWbar[148] + 17./6.*armWWbar[17] - 17./6.*armWWbar[30] + 1./3.*
   armWWbar[127] - 1./4.*armWWbar[33];
   armWWbar[124]=armWWbar[1]*armWWbar[124];
   armWWbar[127]=armWWbar[1]*armWWbar[27];
   armWWbar[129]=armWWbar[40]*armWWbar[27];
   armWWbar[127]=armWWbar[127] + 1./3.*armWWbar[129];
   armWWbar[127]=MMZ*armWWbar[127];
   armWWbar[117]=armWWbar[118] + armWWbar[117] + armWWbar[125] + 1./6.*
   armWWbar[127] + armWWbar[124] + armWWbar[121];
   armWWbar[117]=armWWbar[255]*armWWbar[117];
   armWWbar[118]=armWWbar[119] + armWWbar[163] + 35./4.*armWWbar[13] + 
   armWWbar[169] + 49./6. + armWWbar[172];
   armWWbar[118]=armWWbar[19]*armWWbar[118];
   armWWbar[119]=2 + armWWbar[324];
   armWWbar[119]=armWWbar[159]*armWWbar[119];
   armWWbar[121]=armWWbar[20]*armWWbar[170];
   armWWbar[119]=3*armWWbar[121] + armWWbar[119] + 3*armWWbar[176];
   armWWbar[119]=armWWbar[3]*armWWbar[119];
   armWWbar[121]=13./6.*armWWbar[12] + armWWbar[152] + 31./6.*
   armWWbar[13] + armWWbar[231] - 367./36. + armWWbar[48];
   armWWbar[121]=MMZ*armWWbar[121];
   armWWbar[124]= - armWWbar[80] + armWWbar[192];
   armWWbar[124]=1./2.*armWWbar[124] + armWWbar[22];
   armWWbar[125]=5./2. + armWWbar[46];
   armWWbar[125]=3*armWWbar[125] + 1./2.*armWWbar[12];
   armWWbar[125]=armWWbar[18]*armWWbar[125];
   armWWbar[126]=MMZ*armWWbar[126];
   armWWbar[126]=armWWbar[126] + armWWbar[19];
   armWWbar[126]=armWWbar[11]*armWWbar[126];
   armWWbar[129]=37./4.*armWWbar[12] + 127./12. + armWWbar[169];
   armWWbar[129]=armWWbar[20]*armWWbar[129];
   armWWbar[118]=armWWbar[119] + armWWbar[129] + 3./4.*armWWbar[126] + 
   1./2.*armWWbar[118] + 1./4.*armWWbar[125] + armWWbar[121] + 9*
   armWWbar[23] + 9*armWWbar[124] + 19./8.*armWWbar[21];
   armWWbar[118]=armWWbar[3]*armWWbar[118];
   armWWbar[119]= - 11./2.*armWWbar[98] + armWWbar[264] + armWWbar[260]
    - armWWbar[97] - 445./24.*armWWbar[99] - 331./4.*armWWbar[94] + 
   159221./288. + 59*armWWbar[104];
   armWWbar[121]=7753./27. + 25*armWWbar[13];
   armWWbar[121]= - 557./24.*armWWbar[12] + 1./4.*armWWbar[121] + 5*
   armWWbar[10];
   armWWbar[121]=armWWbar[12]*armWWbar[121];
   armWWbar[119]=armWWbar[121] - 7*armWWbar[10] - 35./12.*armWWbar[13]
    + armWWbar[199] + armWWbar[133] - 37./6.*armWWbar[16] - 
   armWWbar[48] + 727./72.*armWWbar[17] + armWWbar[221] + 7./2.*
   armWWbar[100] + 19./3.*armWWbar[96] + 1./3.*armWWbar[119] + 35./2.*
   armWWbar[103];
   armWWbar[121]=11./12.*armWWbar[23] + 1./3.*armWWbar[21] + 
   armWWbar[162] + armWWbar[161] + armWWbar[153] + 3*armWWbar[240] - 1./
   12.*armWWbar[79];
   armWWbar[121]=armWWbar[112]*armWWbar[121];
   armWWbar[124]=1./3.*armWWbar[92] + 21./4. + 1./3.*armWWbar[93];
   armWWbar[124]=armWWbar[219] + 1./2.*armWWbar[124] + 3*armWWbar[103];
   armWWbar[124]=armWWbar[213] + 2./3.*armWWbar[16] + armWWbar[128] + 1.
   /2.*armWWbar[124] + 2./3.*armWWbar[100];
   armWWbar[124]=armWWbar[112]*armWWbar[124];
   armWWbar[122]=armWWbar[124] + armWWbar[160] + armWWbar[122] + 1./6.*
   armWWbar[82] + 3*armWWbar[85] + 3./2.*armWWbar[109] - 31./16.*
   armWWbar[83] - 34./3.*armWWbar[90] + 5*armWWbar[87];
   armWWbar[122]=MMZ*armWWbar[122];
   armWWbar[124]=armWWbar[169] + armWWbar[196];
   armWWbar[124]=armWWbar[112]*armWWbar[124];
   armWWbar[124]=7*armWWbar[113] + armWWbar[124];
   armWWbar[124]=MMZ*armWWbar[124];
   armWWbar[124]=armWWbar[157] + 1./4.*armWWbar[124] + armWWbar[156] - 
   3./4. + armWWbar[137];
   armWWbar[124]=armWWbar[11]*armWWbar[124];
   armWWbar[125]= - 53./3. + armWWbar[46];
   armWWbar[125]=armWWbar[112]*armWWbar[125];
   armWWbar[125]=armWWbar[139] + armWWbar[125];
   armWWbar[125]=armWWbar[18]*armWWbar[125];
   armWWbar[126]= - 53./12. - armWWbar[46];
   armWWbar[126]=armWWbar[112]*armWWbar[126];
   armWWbar[126]=armWWbar[154] + armWWbar[126];
   armWWbar[126]=armWWbar[20]*armWWbar[126];
   armWWbar[128]=armWWbar[204] - armWWbar[112];
   armWWbar[128]=armWWbar[19]*armWWbar[128];
   armWWbar[129]=11./12.*armWWbar[150] + 2*armWWbar[82];
   armWWbar[129]=MMH*armWWbar[129];
   armWWbar[118]=armWWbar[129] + armWWbar[118] + 1./2.*armWWbar[126] + 
   armWWbar[124] + 1./4.*armWWbar[128] + 1./4.*armWWbar[125] + 
   armWWbar[122] + 1./2.*armWWbar[119] + armWWbar[121];
   armWWbar[118]=armWWbar[76]*armWWbar[118];
   armWWbar[119]=armWWbar[134] - 13./3.*armWWbar[32] - 55./18.*
   armWWbar[15] - 563./36. - armWWbar[29];
   armWWbar[121]=7./2.*armWWbar[7];
   armWWbar[122]= - 253./3. + armWWbar[121];
   armWWbar[122]=armWWbar[12]*armWWbar[122];
   armWWbar[124]=1./2. + armWWbar[194];
   armWWbar[124]=armWWbar[6]*armWWbar[124];
   armWWbar[125]= - 1./3. - armWWbar[7];
   armWWbar[125]=armWWbar[7]*armWWbar[125];
   armWWbar[119]=armWWbar[124] + 1./18.*armWWbar[122] + 1./4.*
   armWWbar[125] - 19./9.*armWWbar[17] + 85./36.*armWWbar[30] + 1./4.*
   armWWbar[119] - 10./3.*armWWbar[33];
   armWWbar[119]=armWWbar[1]*armWWbar[119];
   armWWbar[122]=1./3.*armWWbar[134] - 23./3.*armWWbar[32] - 229./18.*
   armWWbar[15] - 7993./324. - 3*armWWbar[29];
   armWWbar[121]= - 781./27. + armWWbar[121];
   armWWbar[121]=armWWbar[12]*armWWbar[121];
   armWWbar[124]=11./6. + armWWbar[247];
   armWWbar[124]=armWWbar[6]*armWWbar[124];
   armWWbar[121]=1./3.*armWWbar[124] + 1./6.*armWWbar[121] + 1./12.*
   armWWbar[125] - 59./27.*armWWbar[17] + 317./108.*armWWbar[30] + 1./4.
   *armWWbar[122] - 110./27.*armWWbar[33];
   armWWbar[121]=armWWbar[40]*armWWbar[121];
   armWWbar[122]=MMZ*armWWbar[281];
   armWWbar[122]=armWWbar[122] + armWWbar[140];
   armWWbar[122]=armWWbar[3]*armWWbar[122];
   armWWbar[117]=armWWbar[117] + armWWbar[118] + armWWbar[122] + 1./3.*
   armWWbar[127] + armWWbar[119] + armWWbar[121];
   armWWbar[117]=armWWbar[255]*armWWbar[117];
   armWWbar[118]=pow(CW,2);
   armWWbar[119]=17./3. - 2*armWWbar[118];
   armWWbar[121]= - 4*armWWbar[118];
   armWWbar[122]= - 53./3. + armWWbar[121];
   armWWbar[122]=armWWbar[13]*armWWbar[122];
   armWWbar[119]= - 10*armWWbar[12] + 2*armWWbar[119] + armWWbar[122];
   armWWbar[119]=MMZ*armWWbar[119];
   armWWbar[122]=2 + armWWbar[232];
   armWWbar[122]=armWWbar[3]*armWWbar[159]*armWWbar[122];
   armWWbar[124]= - 5 - 3*armWWbar[13];
   armWWbar[124]=armWWbar[20]*armWWbar[124];
   armWWbar[119]=armWWbar[122] + armWWbar[119] + 2*armWWbar[124];
   armWWbar[119]=armWWbar[3]*armWWbar[119];
   armWWbar[122]=4*armWWbar[118];
   armWWbar[124]=13 + armWWbar[122];
   armWWbar[122]=113./3. + armWWbar[122];
   armWWbar[122]=armWWbar[13]*armWWbar[122];
   armWWbar[122]=6*armWWbar[12] + 5*armWWbar[124] + armWWbar[122];
   armWWbar[122]=armWWbar[12]*armWWbar[122];
   armWWbar[124]=armWWbar[83] + armWWbar[88] - 2*armWWbar[87];
   armWWbar[124]=armWWbar[118]*armWWbar[124];
   armWWbar[124]= - armWWbar[82] + 2*armWWbar[124] + armWWbar[86] + 23./
   3.*armWWbar[83] - 4*armWWbar[90] - 23./3.*armWWbar[87];
   armWWbar[125]=20*armWWbar[16] + 4*armWWbar[17] - 32*armWWbar[100] - 
   4*armWWbar[92] - 55 - 4*armWWbar[93];
   armWWbar[125]=armWWbar[112]*armWWbar[125];
   armWWbar[124]=8*armWWbar[124] + armWWbar[125];
   armWWbar[124]=MMZ*armWWbar[124];
   armWWbar[125]= - 4 + armWWbar[12];
   armWWbar[125]=MMZ*armWWbar[112]*armWWbar[125];
   armWWbar[125]=armWWbar[123] + armWWbar[125];
   armWWbar[125]=armWWbar[11]*armWWbar[125];
   armWWbar[126]= - 3637./6. - 88*armWWbar[104];
   armWWbar[127]=61./3. + 8*armWWbar[118];
   armWWbar[127]=armWWbar[99]*armWWbar[127];
   armWWbar[128]= - 67./3. - 8*armWWbar[118];
   armWWbar[128]=armWWbar[17]*armWWbar[128];
   armWWbar[121]= - 73./9. + armWWbar[121];
   armWWbar[121]=armWWbar[13]*armWWbar[121];
   armWWbar[118]=2*armWWbar[119] + 4*armWWbar[125] + armWWbar[124] + 2*
   armWWbar[122] + 2*armWWbar[121] + 4*armWWbar[128] - 44./3.*
   armWWbar[33] + 4*armWWbar[127] - 8*armWWbar[94] + 1./3.*
   armWWbar[126] - 64*armWWbar[118];
   armWWbar[118]=armWWbar[76]*armWWbar[118];
   armWWbar[119]=armWWbar[7]*armWWbar[123];
   armWWbar[121]=2*armWWbar[120] - 8./3.*armWWbar[12] + armWWbar[119]
    + armWWbar[200] + 2*armWWbar[17] - 2*armWWbar[30] + 11 + 2*
   armWWbar[33];
   armWWbar[121]=armWWbar[1]*armWWbar[121];
   armWWbar[122]=10*armWWbar[17] - 10*armWWbar[30] + 169./3. + 10*
   armWWbar[33];
   armWWbar[122]=1./3.*armWWbar[122] + armWWbar[13];
   armWWbar[119]=10./9.*armWWbar[120] - 40./27.*armWWbar[12] + 1./3.*
   armWWbar[122] + armWWbar[119];
   armWWbar[119]=armWWbar[40]*armWWbar[119];
   armWWbar[119]=1./3.*armWWbar[121] + armWWbar[119];
   armWWbar[120]= - armWWbar[24] + 2*armWWbar[26] - armWWbar[25];
   armWWbar[120]=armWWbar[1]*armWWbar[120];
   armWWbar[121]=1./9.*armWWbar[75] + armWWbar[26];
   armWWbar[121]= - 2./9.*armWWbar[27] - armWWbar[24] + 2*armWWbar[121]
    - armWWbar[25];
   armWWbar[121]=armWWbar[40]*armWWbar[121];
   armWWbar[120]=1./3.*armWWbar[120] + armWWbar[121];
   armWWbar[120]=MMZ*armWWbar[120];
   armWWbar[121]=armWWbar[215] + 5./3.*armWWbar[207];
   armWWbar[121]=armWWbar[3]*MMZ*armWWbar[121];
   armWWbar[119]=4./3.*armWWbar[121] + 2*armWWbar[119] + armWWbar[120];

      mWWbarret = armWWbar[114] + armWWbar[115] + armWWbar[116] + 
      armWWbar[117] + armWWbar[118] + 2*armWWbar[119];
      return mWWbarret;
}
