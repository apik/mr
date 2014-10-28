#include <WW.hpp>
std::complex<long double>
WW<OS>::my20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuWW[342], yuWWret;

    aryuWW[1]=double(nL + nH);
    aryuWW[2]=pow(CW,-1);
    aryuWW[3]=pow(MMH,-1);
    aryuWW[4]=pow(MMZ,-1);
    aryuWW[5]=pow(SW,-1);
    aryuWW[6]=std::real(Tsil::B(0,0,MMZ,mu2));
    aryuWW[7]=std::real(Tsil::B(0,0,MMW,mu2));
    aryuWW[8]=Tsil::I2(0,0,MMZ,mu2);
    aryuWW[9]=Tsil::I2(0,0,MMW,mu2);
    aryuWW[10]=Tsil::B(MMZ,MMH,MMZ,mu2);
    aryuWW[11]=Tsil::B(MMW,MMH,MMW,mu2);
    aryuWW[12]=Tsil::B(MMW,MMZ,MMW,mu2);
    aryuWW[13]=Tsil::B(MMW,MMW,MMZ,mu2);
    aryuWW[14]=Tsil::B(0,MMH,MMW,mu2);
    aryuWW[15]=Tsil::B(0,MMZ,MMW,mu2);
    aryuWW[16]=Tsil::Beps(MMW,MMH,MMW,mu2);
    aryuWW[17]=Tsil::Beps(MMW,MMZ,MMW,mu2);
    aryuWW[18]=Tsil::A(MMH,mu2);
    aryuWW[19]=Tsil::A(MMZ,mu2);
    aryuWW[20]=Tsil::A(MMW,mu2);
    aryuWW[21]=Tsil::Aeps(MMH,mu2);
    aryuWW[22]=Tsil::Aeps(MMZ,mu2);
    aryuWW[23]=Tsil::Aeps(MMW,mu2);
    aryuWW[24]=protW0Z00->M(0);
    aryuWW[25]=prot0W0Z0->M(0);
    aryuWW[26]=prot00W00->M(0);
    aryuWW[27]=prot0000Z->M(0);
    aryuWW[28]=protHW00->Uxzuv(0);
    aryuWW[29]=protW0Z00->Uzxyv(0);
    aryuWW[30]=protWtZ00->Uxzuv(0);
    aryuWW[31]=protHW00->Txuv(0);
    aryuWW[32]=prot0W0Z0->Tuxv(0);
    aryuWW[33]=protWtZ00->Txuv(0);
    aryuWW[34]=double(nH);
    aryuWW[35]=pow(MMt,-1);
    aryuWW[36]=Tsil::B(MMt,MMt,MMZ,mu2);
    aryuWW[37]=Tsil::B(0,MMt,MMW,mu2);
    aryuWW[38]=Tsil::A(MMt,mu2);
    aryuWW[39]=Tsil::B(MMt,MMt,MMH,mu2);
    aryuWW[40]=double(nL);
    aryuWW[41]=Tsil::I2(MMH,MMt,MMt,mu2);
    aryuWW[42]=Tsil::I2(MMZ,MMt,MMt,mu2);
    aryuWW[43]=Tsil::I2(0,MMW,MMt,mu2);
    aryuWW[44]=Tsil::B(MMH,MMH,MMH,mu2);
    aryuWW[45]=Tsil::B(MMH,MMt,MMt,mu2);
    aryuWW[46]=Tsil::B(MMZ,MMZ,MMH,mu2);
    aryuWW[47]=Tsil::B(MMZ,MMt,MMt,mu2);
    aryuWW[48]=Tsil::B(MMW,MMW,MMH,mu2);
    aryuWW[49]=std::real(Tsil::B(0,MMW,MMt,mu2));
    aryuWW[50]=Tsil::Beps(0,MMt,MMW,mu2);
    aryuWW[51]=Tsil::Aeps(MMt,mu2);
    aryuWW[52]=prot00tt0->M(0);
    aryuWW[53]=prot00tt0->Tuxv(0);
    aryuWW[54]=protWtZ00->M(0);
    aryuWW[55]=protW0Htt->M(0);
    aryuWW[56]=protW0Ztt->M(0);
    aryuWW[57]=prot0Wt0t->M(0);
    aryuWW[58]=prot00Wt0->M(0);
    aryuWW[59]=prot00ttZ->M(0);
    aryuWW[60]=protW0Htt->Uzxyv(0);
    aryuWW[61]=protWtZ00->Uzxyv(0);
    aryuWW[62]=protW0Htt->Uxzuv(0);
    aryuWW[63]=protW0Ztt->Uxzuv(0);
    aryuWW[64]=prot0Wt0t->Uyuzv(0);
    aryuWW[65]=protW0Htt->Uyuzv(0);
    aryuWW[66]=protW0Ztt->Uyuzv(0);
    aryuWW[67]=protWtZ00->Uuyxv(0);
    aryuWW[68]=protW0Htt->Tzyv(0);
    aryuWW[69]=protWtZ00->Tzyv(0);
    aryuWW[70]=protW0Htt->Tvyz(0);
    aryuWW[71]=protWtZ00->Tyzv(0);
    aryuWW[72]=protW0Htt->Suxv(0);
    aryuWW[73]=protW0Htt->Svyz(0);
    aryuWW[74]=protWtZ00->Svyz(0);
    aryuWW[75]=prot00000->M(0);
    aryuWW[76]=double(boson);
    aryuWW[77]=Tsil::I2(MMH,MMH,MMH,mu2);
    aryuWW[78]=Tsil::I2(MMH,MMZ,MMZ,mu2);
    aryuWW[79]=Tsil::I2(MMH,MMW,MMW,mu2);
    aryuWW[80]=Tsil::I2(MMW,MMW,MMZ,mu2);
    aryuWW[81]=protWHHWW->M(0);
    aryuWW[82]=protWHZWW->M(0);
    aryuWW[83]=protWZZWW->M(0);
    aryuWW[84]=protWWHHH->M(0);
    aryuWW[85]=protWWHZZ->M(0);
    aryuWW[86]=protW0HWW->M(0);
    aryuWW[87]=protW0ZWW->M(0);
    aryuWW[88]=prot0WW0W->M(0);
    aryuWW[89]=protWH0H->Vxzuv(0);
    aryuWW[90]=protWZ0Z->Vxzuv(0);
    aryuWW[91]=protWHHWW->Uzxyv(0);
    aryuWW[92]=protWHZWW->Uyuzv(0);
    aryuWW[93]=protWHZWW->Uzxyv(0);
    aryuWW[94]=protWZZWW->Uzxyv(0);
    aryuWW[95]=protWWHHH->Uyuzv(0);
    aryuWW[96]=protWWHZZ->Uxzuv(0);
    aryuWW[97]=protWHHWW->Uxzuv(0);
    aryuWW[98]=protWWHZZ->Uyuzv(0);
    aryuWW[99]=protWHZWW->Uxzuv(0);
    aryuWW[100]=protWHZWW->Tyzv(0);
    aryuWW[101]=protWHHWW->Tyzv(0);
    aryuWW[102]=protW0HWW->Tzyv(0);
    aryuWW[103]=protWHZWW->Tzyv(0);
    aryuWW[104]=protW0ZWW->Tzyv(0);
    aryuWW[105]=protWHHWW->Svyz(0);
    aryuWW[106]=protWHZWW->Svyz(0);
    aryuWW[107]=protWZZWW->Svyz(0);
    aryuWW[108]=protWtZ00->Suxv(0);
    aryuWW[109]=protWWZZH->M(0);
    aryuWW[110]=1/(MMt - MMW);
    aryuWW[111]=1/(4*MMt - MMZ);
    aryuWW[112]=1/( - 4*MMW + MMH);
    aryuWW[113]=1/( - 4*MMZ + MMH);
   aryuWW[114]= - 1./6. - aryuWW[7];
   aryuWW[115]=1./2.*aryuWW[6];
   aryuWW[116]=1./3.*aryuWW[114] + aryuWW[115];
   aryuWW[116]=aryuWW[1]*aryuWW[116];
   aryuWW[117]= - 25*aryuWW[47];
   aryuWW[118]=137./12. + aryuWW[117];
   aryuWW[119]=19./24.*aryuWW[13];
   aryuWW[120]=11./18.*aryuWW[6];
   aryuWW[121]=aryuWW[120] + 7./54. - aryuWW[7];
   aryuWW[121]=aryuWW[40]*aryuWW[121];
   aryuWW[122]=aryuWW[18]*aryuWW[112];
   aryuWW[123]=1./2.*aryuWW[122];
   aryuWW[124]= - 3./2.*aryuWW[49];
   aryuWW[125]= - 3*aryuWW[45];
   aryuWW[118]=aryuWW[123] + aryuWW[121] + aryuWW[116] + 233./24.*
   aryuWW[12] + 3./2.*aryuWW[10] + aryuWW[119] + aryuWW[124] + 1./8.*
   aryuWW[118] + aryuWW[125];
   aryuWW[118]=aryuWW[37]*aryuWW[118];
   aryuWW[126]=3*aryuWW[41];
   aryuWW[127]= - 8*aryuWW[51] + aryuWW[126] + 2*aryuWW[43];
   aryuWW[128]= - aryuWW[18] - 9./4.*aryuWW[19];
   aryuWW[128]=aryuWW[37]*aryuWW[128];
   aryuWW[129]= - 5./2.*aryuWW[37] - 3 - 2*aryuWW[39];
   aryuWW[129]=aryuWW[20]*aryuWW[129];
   aryuWW[130]= - 2*aryuWW[49] - 3 - 8*aryuWW[45];
   aryuWW[130]=aryuWW[38]*aryuWW[130];
   aryuWW[131]= - 3*aryuWW[39];
   aryuWW[132]= - 5 + aryuWW[131];
   aryuWW[132]=aryuWW[18]*aryuWW[132];
   aryuWW[127]=3*aryuWW[130] + 3*aryuWW[129] + aryuWW[128] + 
   aryuWW[132] - 6*aryuWW[23] + 3*aryuWW[127] - 10*aryuWW[21];
   aryuWW[127]=aryuWW[3]*aryuWW[127];
   aryuWW[128]= - aryuWW[58] - 1./3.*aryuWW[52];
   aryuWW[129]=2*aryuWW[57];
   aryuWW[130]=aryuWW[128] + aryuWW[129];
   aryuWW[133]=aryuWW[59] - 1./2.*aryuWW[54] + 2*aryuWW[130] - 1./2.*
   aryuWW[56];
   aryuWW[134]= - 4*aryuWW[45];
   aryuWW[135]= - aryuWW[49] + 3 + aryuWW[134];
   aryuWW[135]=aryuWW[3]*aryuWW[135];
   aryuWW[136]=5 - aryuWW[70];
   aryuWW[136]=aryuWW[112]*aryuWW[136];
   aryuWW[133]=3*aryuWW[135] + aryuWW[136] + 1./3.*aryuWW[133] - 2*
   aryuWW[55];
   aryuWW[133]=MMt*aryuWW[133];
   aryuWW[135]=4*aryuWW[64];
   aryuWW[137]=31./3.*aryuWW[53];
   aryuWW[138]=aryuWW[137] - 935./12. + aryuWW[135];
   aryuWW[139]= - 5./2.*aryuWW[21];
   aryuWW[126]= - aryuWW[23] + aryuWW[139] - 6*aryuWW[51] + aryuWW[43]
    + aryuWW[126] - 2*aryuWW[73];
   aryuWW[126]=aryuWW[112]*aryuWW[126];
   aryuWW[140]=1 + aryuWW[39];
   aryuWW[140]=aryuWW[112]*aryuWW[140];
   aryuWW[141]= - aryuWW[37]*aryuWW[112];
   aryuWW[142]=6*aryuWW[140] + aryuWW[141];
   aryuWW[142]=aryuWW[20]*aryuWW[142];
   aryuWW[143]=41./8.*aryuWW[61];
   aryuWW[144]=7./3.*aryuWW[67];
   aryuWW[145]=17./32.*aryuWW[66];
   aryuWW[146]= - 1./6.*aryuWW[63];
   aryuWW[147]= - 6*aryuWW[62];
   aryuWW[148]=17./12.*aryuWW[33];
   aryuWW[149]=1./6.*aryuWW[30];
   aryuWW[150]= - 35./6.*aryuWW[69];
   aryuWW[151]=311./32.*aryuWW[71];
   aryuWW[152]= - 491./96.*aryuWW[50];
   aryuWW[153]= - 73./24.*aryuWW[70];
   aryuWW[154]= - 5./2.*aryuWW[68];
   aryuWW[155]=9*aryuWW[45];
   aryuWW[156]=3./2.*aryuWW[49];
   aryuWW[157]= - 41./8.*aryuWW[17];
   aryuWW[158]=6*aryuWW[16];
   aryuWW[159]= - 6*aryuWW[39];
   aryuWW[160]=3*aryuWW[39];
   aryuWW[161]=1 + aryuWW[160];
   aryuWW[161]=aryuWW[11]*aryuWW[161];
   aryuWW[162]=2*aryuWW[161];
   aryuWW[163]= - aryuWW[38]*aryuWW[112];
   aryuWW[164]=4*aryuWW[163];
   aryuWW[165]=MMH*aryuWW[55];
   aryuWW[166]=3*aryuWW[165];
   aryuWW[167]=9./4.*aryuWW[65];
   aryuWW[168]=5./2. + aryuWW[131];
   aryuWW[168]=aryuWW[18]*aryuWW[112]*aryuWW[168];
   aryuWW[138]=aryuWW[133] + aryuWW[166] + aryuWW[127] + aryuWW[164] + 
   aryuWW[142] + aryuWW[118] + aryuWW[162] + aryuWW[168] + aryuWW[126]
    + aryuWW[159] + aryuWW[158] + aryuWW[157] + aryuWW[156] + 
   aryuWW[155] + aryuWW[154] + aryuWW[153] + aryuWW[152] + aryuWW[151]
    + aryuWW[150] + aryuWW[149] + aryuWW[148] + aryuWW[147] + 
   aryuWW[167] + aryuWW[146] + aryuWW[145] + aryuWW[144] + 1./3.*
   aryuWW[138] + aryuWW[143];
   aryuWW[138]=MMt*aryuWW[138];
   aryuWW[169]=2879./36. + aryuWW[117];
   aryuWW[170]=6*aryuWW[45];
   aryuWW[171]=27./2.*aryuWW[44];
   aryuWW[172]=3*aryuWW[48];
   aryuWW[173]=3./2.*aryuWW[46];
   aryuWW[174]=89./6.*aryuWW[12];
   aryuWW[175]=2*aryuWW[11];
   aryuWW[176]= - 89./96.*aryuWW[37];
   aryuWW[177]=aryuWW[20]*aryuWW[112];
   aryuWW[178]=6*aryuWW[177];
   aryuWW[179]=aryuWW[38]*aryuWW[112];
   aryuWW[180]=3*aryuWW[179];
   aryuWW[181]= - 3./4.*aryuWW[49];
   aryuWW[182]= - 1./2.*aryuWW[10];
   aryuWW[169]=aryuWW[180] + aryuWW[178] + aryuWW[176] + aryuWW[175] + 
   aryuWW[121] + aryuWW[116] + aryuWW[174] + aryuWW[182] + aryuWW[119]
    + aryuWW[173] + aryuWW[172] + aryuWW[171] + aryuWW[181] + 1./8.*
   aryuWW[169] + aryuWW[170];
   aryuWW[169]=aryuWW[38]*aryuWW[169];
   aryuWW[183]=3*aryuWW[65];
   aryuWW[184]= - 3805./108. + aryuWW[183];
   aryuWW[185]=1./2.*aryuWW[60];
   aryuWW[186]= - 3./4.*aryuWW[50];
   aryuWW[187]=5./2.*aryuWW[70];
   aryuWW[188]= - 21./4.*aryuWW[68];
   aryuWW[189]= - 7./2.*aryuWW[16];
   aryuWW[190]=7./54.*aryuWW[36];
   aryuWW[191]=3*aryuWW[62];
   aryuWW[184]=aryuWW[190] + aryuWW[160] + aryuWW[189] + aryuWW[125] + 
   aryuWW[188] + aryuWW[187] + aryuWW[186] + aryuWW[185] + 1./4.*
   aryuWW[184] + aryuWW[191];
   aryuWW[192]= - 3./2.*aryuWW[39];
   aryuWW[193]=aryuWW[190] - 55./27. + aryuWW[192];
   aryuWW[193]=aryuWW[11]*aryuWW[193];
   aryuWW[194]=3*aryuWW[45];
   aryuWW[195]=1./3. + aryuWW[194];
   aryuWW[196]= - 1./3.*aryuWW[10];
   aryuWW[195]=1./2.*aryuWW[195] + aryuWW[196];
   aryuWW[197]=1./6.*aryuWW[11];
   aryuWW[198]=aryuWW[195] + aryuWW[197];
   aryuWW[198]=aryuWW[37]*aryuWW[198];
   aryuWW[199]=1./3.*aryuWW[10];
   aryuWW[184]=1./2.*aryuWW[198] + aryuWW[193] + 1./2.*aryuWW[184] + 
   aryuWW[199];
   aryuWW[184]=MMH*aryuWW[184];
   aryuWW[193]= - 43*aryuWW[74];
   aryuWW[198]=aryuWW[193] + 79./12.*aryuWW[42];
   aryuWW[200]= - 13*aryuWW[72];
   aryuWW[198]=1./8.*aryuWW[198] + aryuWW[200];
   aryuWW[201]=965 + 1361*aryuWW[36];
   aryuWW[202]= - 7*aryuWW[37];
   aryuWW[201]=1./27.*aryuWW[201] + aryuWW[202];
   aryuWW[201]=aryuWW[20]*aryuWW[201];
   aryuWW[203]= - 7./18.*aryuWW[36];
   aryuWW[204]=aryuWW[203] + 475./36. + aryuWW[160];
   aryuWW[204]=1./2.*aryuWW[18]*aryuWW[204];
   aryuWW[205]=aryuWW[18] - 5./4.*aryuWW[19];
   aryuWW[206]= - 9*aryuWW[38];
   aryuWW[205]=aryuWW[206] + 3*aryuWW[205] - 23*aryuWW[20];
   aryuWW[205]=aryuWW[3]*aryuWW[38]*aryuWW[205];
   aryuWW[207]= - 5129./16. - 331*aryuWW[36];
   aryuWW[207]=aryuWW[19]*aryuWW[207];
   aryuWW[208]= - 71./12.*aryuWW[73];
   aryuWW[209]= - 5./12.*aryuWW[8];
   aryuWW[210]=23./8.*aryuWW[43];
   aryuWW[211]=115./24.*aryuWW[21];
   aryuWW[212]=5*aryuWW[18];
   aryuWW[213]=aryuWW[212] + 541./24.*aryuWW[19];
   aryuWW[213]=1./8.*aryuWW[37]*aryuWW[213];
   aryuWW[214]= - 3./2.*aryuWW[41];
   aryuWW[138]=aryuWW[138] + aryuWW[184] + aryuWW[205] + aryuWW[169] + 
   1./4.*aryuWW[201] + aryuWW[213] + 1./36.*aryuWW[207] + aryuWW[204]
    + 1795./72.*aryuWW[23] + aryuWW[211] + 1621./576.*aryuWW[22] + 2503.
   /96.*aryuWW[51] + aryuWW[210] + aryuWW[209] + aryuWW[208] + 1./3.*
   aryuWW[198] + aryuWW[214];
   aryuWW[138]=MMt*aryuWW[138];
   aryuWW[169]= - 11*aryuWW[47];
   aryuWW[198]= - 257./18. + aryuWW[169];
   aryuWW[198]=aryuWW[10] + 1./12.*aryuWW[13] + aryuWW[124] + 1./2.*
   aryuWW[198] + aryuWW[125];
   aryuWW[201]=7./3.*aryuWW[12];
   aryuWW[207]= - 1./3. + aryuWW[6];
   aryuWW[215]=aryuWW[1]*aryuWW[207];
   aryuWW[216]=1./2.*aryuWW[215];
   aryuWW[207]=aryuWW[40]*aryuWW[207];
   aryuWW[217]=11./18.*aryuWW[207];
   aryuWW[218]=1./4.*aryuWW[122];
   aryuWW[219]=1./2.*aryuWW[11];
   aryuWW[198]=aryuWW[219] + aryuWW[218] + aryuWW[217] + aryuWW[216] + 
   1./2.*aryuWW[198] + aryuWW[201];
   aryuWW[198]=aryuWW[37]*aryuWW[198];
   aryuWW[220]=3./2.*aryuWW[41];
   aryuWW[221]= - 4*aryuWW[51] + aryuWW[220] + aryuWW[43];
   aryuWW[222]= - 3*aryuWW[23];
   aryuWW[132]=1./2.*aryuWW[132] + aryuWW[222] + 3*aryuWW[221] - 5*
   aryuWW[21];
   aryuWW[221]= - aryuWW[18] - 3./2.*aryuWW[19];
   aryuWW[221]=aryuWW[37]*aryuWW[221];
   aryuWW[223]= - 3./2. - aryuWW[39];
   aryuWW[224]=aryuWW[223] - aryuWW[37];
   aryuWW[224]=aryuWW[20]*aryuWW[224];
   aryuWW[134]= - aryuWW[49] - 3./2. + aryuWW[134];
   aryuWW[134]=3*aryuWW[38]*aryuWW[134];
   aryuWW[221]=aryuWW[134] + 3*aryuWW[224] + aryuWW[132] + 1./2.*
   aryuWW[221];
   aryuWW[221]=aryuWW[3]*aryuWW[221];
   aryuWW[224]= - 3*aryuWW[51];
   aryuWW[225]=1./2.*aryuWW[43];
   aryuWW[226]= - 1./2.*aryuWW[23];
   aryuWW[220]=aryuWW[226] - 5./4.*aryuWW[21] + aryuWW[224] + 
   aryuWW[225] + aryuWW[220] - aryuWW[73];
   aryuWW[220]=aryuWW[112]*aryuWW[220];
   aryuWW[130]=aryuWW[130] - aryuWW[56];
   aryuWW[130]=aryuWW[59] + 2*aryuWW[130] + aryuWW[54];
   aryuWW[136]=1./2.*aryuWW[136];
   aryuWW[227]= - 1./2.*aryuWW[49] + 3./2. - 2*aryuWW[45];
   aryuWW[227]=3*aryuWW[3]*aryuWW[227];
   aryuWW[130]=aryuWW[227] + aryuWW[136] + 1./3.*aryuWW[130] - 
   aryuWW[55];
   aryuWW[130]=MMt*aryuWW[130];
   aryuWW[228]=aryuWW[137] - 6749./288. + aryuWW[135];
   aryuWW[141]=1./2.*aryuWW[141];
   aryuWW[140]=3*aryuWW[140] + aryuWW[141];
   aryuWW[140]=aryuWW[20]*aryuWW[140];
   aryuWW[229]= - 1./3.*aryuWW[30];
   aryuWW[230]= - 3*aryuWW[62];
   aryuWW[231]= - 73./48.*aryuWW[70];
   aryuWW[232]= - 5./4.*aryuWW[68];
   aryuWW[233]=3./4.*aryuWW[49];
   aryuWW[234]=3*aryuWW[16];
   aryuWW[235]=1./2.*aryuWW[168];
   aryuWW[236]=2*aryuWW[163];
   aryuWW[165]=3./2.*aryuWW[165];
   aryuWW[237]=9./2.*aryuWW[45];
   aryuWW[130]=aryuWW[130] + aryuWW[165] + aryuWW[221] + aryuWW[236] + 
   aryuWW[140] + aryuWW[198] + aryuWW[161] + aryuWW[235] + aryuWW[220]
    + aryuWW[131] + aryuWW[234] - 7./8.*aryuWW[17] + aryuWW[233] + 
   aryuWW[237] + aryuWW[232] + aryuWW[231] - 139./48.*aryuWW[50] + 29./
   8.*aryuWW[71] - 67./12.*aryuWW[69] + aryuWW[229] + 25./24.*
   aryuWW[33] + aryuWW[230] + 9./8.*aryuWW[65] - 2./3.*aryuWW[63] + 11./
   16.*aryuWW[66] + 13./12.*aryuWW[67] + 1./3.*aryuWW[228] + 15./8.*
   aryuWW[61];
   aryuWW[130]=MMt*aryuWW[130];
   aryuWW[198]= - 2279./108. + aryuWW[169];
   aryuWW[221]= - 3./8.*aryuWW[49];
   aryuWW[228]=27./4.*aryuWW[44];
   aryuWW[238]=3./2.*aryuWW[48];
   aryuWW[239]=3./4.*aryuWW[46];
   aryuWW[240]=3*aryuWW[177];
   aryuWW[179]=3./2.*aryuWW[179];
   aryuWW[241]=1./2.*aryuWW[10];
   aryuWW[198]=aryuWW[179] + aryuWW[240] - 25./48.*aryuWW[37] + 
   aryuWW[219] + aryuWW[217] + aryuWW[216] + 77./24.*aryuWW[12] + 
   aryuWW[241] + 1./24.*aryuWW[13] + aryuWW[239] + aryuWW[238] + 
   aryuWW[228] + aryuWW[221] + 1./4.*aryuWW[198] + aryuWW[194];
   aryuWW[198]=aryuWW[38]*aryuWW[198];
   aryuWW[242]= - 437./12. + aryuWW[183];
   aryuWW[242]=aryuWW[160] + aryuWW[189] + aryuWW[125] + aryuWW[188] + 
   aryuWW[187] + aryuWW[186] + aryuWW[185] + 1./4.*aryuWW[242] + 
   aryuWW[191];
   aryuWW[243]= - 3./4.*aryuWW[39];
   aryuWW[244]= - 1 + aryuWW[243];
   aryuWW[244]=aryuWW[11]*aryuWW[244];
   aryuWW[245]=aryuWW[194] - aryuWW[11];
   aryuWW[245]=aryuWW[37]*aryuWW[245];
   aryuWW[242]=1./8.*aryuWW[245] + 1./4.*aryuWW[242] + aryuWW[244];
   aryuWW[242]=MMH*aryuWW[242];
   aryuWW[244]=3*aryuWW[18];
   aryuWW[245]=1./2.*aryuWW[19];
   aryuWW[206]=aryuWW[206] - 5*aryuWW[20] + aryuWW[244] + aryuWW[245];
   aryuWW[206]=aryuWW[3]*aryuWW[38]*aryuWW[206];
   aryuWW[246]=35*aryuWW[74] + 1303./18.*aryuWW[42];
   aryuWW[246]=1./8.*aryuWW[246] - 7*aryuWW[72];
   aryuWW[247]=169./12. + aryuWW[160];
   aryuWW[247]=aryuWW[18]*aryuWW[247];
   aryuWW[248]= - 505./36. - aryuWW[36];
   aryuWW[248]=aryuWW[19]*aryuWW[248];
   aryuWW[249]=aryuWW[18] + 23./18.*aryuWW[19];
   aryuWW[249]=aryuWW[37]*aryuWW[249];
   aryuWW[250]= - 29 - 25./8.*aryuWW[36];
   aryuWW[250]=1./3.*aryuWW[250] - 17./8.*aryuWW[37];
   aryuWW[250]=aryuWW[20]*aryuWW[250];
   aryuWW[251]= - 3./4.*aryuWW[41];
   aryuWW[252]= - 71./24.*aryuWW[73];
   aryuWW[253]=115./48.*aryuWW[21];
   aryuWW[130]=aryuWW[130] + aryuWW[242] + 1./2.*aryuWW[206] + 
   aryuWW[198] + 1./3.*aryuWW[250] + 7./16.*aryuWW[249] + 1./24.*
   aryuWW[248] + 1./4.*aryuWW[247] + 1331./54.*aryuWW[23] + aryuWW[253]
    - 1907./864.*aryuWW[22] + 2779./432.*aryuWW[51] + 3./4.*aryuWW[43]
    - 1./6.*aryuWW[8] + aryuWW[252] + 1./3.*aryuWW[246] + aryuWW[251];
   aryuWW[130]=MMt*aryuWW[130];
   aryuWW[198]= - 83./6. - 15*aryuWW[61];
   aryuWW[206]=7*aryuWW[47];
   aryuWW[242]= - 7*aryuWW[12] + 5./12. + aryuWW[206];
   aryuWW[242]=aryuWW[37]*aryuWW[242];
   aryuWW[198]=1./12.*aryuWW[242] + 15./4.*aryuWW[17] + 125./144.*
   aryuWW[50] - 109./48.*aryuWW[71] - 67./9.*aryuWW[69] + aryuWW[229]
    + 1./3.*aryuWW[63] - 77./144.*aryuWW[66] + 1./4.*aryuWW[198] - 1./3.
   *aryuWW[67];
   aryuWW[242]=aryuWW[56] + aryuWW[54];
   aryuWW[246]=1./3.*aryuWW[59];
   aryuWW[242]=1./2.*aryuWW[242] + aryuWW[246];
   aryuWW[242]=MMt*aryuWW[242];
   aryuWW[198]=1./2.*aryuWW[198] + 1./3.*aryuWW[242];
   aryuWW[198]=MMt*aryuWW[198];
   aryuWW[242]=35./12. + aryuWW[47];
   aryuWW[247]= - 13*aryuWW[12];
   aryuWW[242]= - 49./48.*aryuWW[37] + 7./4.*aryuWW[242] + aryuWW[247];
   aryuWW[242]=aryuWW[38]*aryuWW[242];
   aryuWW[248]=43./9.*aryuWW[74] - 3./4.*aryuWW[42];
   aryuWW[249]= - 53./8. + aryuWW[36];
   aryuWW[249]=aryuWW[19]*aryuWW[249];
   aryuWW[250]=aryuWW[37]*aryuWW[19];
   aryuWW[254]= - 2 - 7./8.*aryuWW[36];
   aryuWW[254]=aryuWW[20]*aryuWW[254];
   aryuWW[198]=aryuWW[198] + 1./6.*aryuWW[242] + 1./9.*aryuWW[254] + 
   205./576.*aryuWW[250] + 7./72.*aryuWW[249] + 1./8.*aryuWW[23] + 13./
   192.*aryuWW[22] - 289./288.*aryuWW[51] + 1./8.*aryuWW[43] - 1./36.*
   aryuWW[8] + 1./8.*aryuWW[248] - 1./3.*aryuWW[72];
   aryuWW[198]=MMt*aryuWW[198];
   aryuWW[242]= - 1 - aryuWW[61];
   aryuWW[248]=1./2.*aryuWW[17];
   aryuWW[242]=aryuWW[248] - 1./2.*aryuWW[71] + 1./2.*aryuWW[242] - 
   aryuWW[69];
   aryuWW[242]=MMt*aryuWW[242];
   aryuWW[249]=1 - aryuWW[12];
   aryuWW[254]=aryuWW[38]*aryuWW[249];
   aryuWW[255]= - aryuWW[51] + aryuWW[254];
   aryuWW[242]=1./2.*aryuWW[255] + aryuWW[242];
   aryuWW[255]=pow(aryuWW[2],2);
   aryuWW[242]=aryuWW[255]*MMt*aryuWW[242];
   aryuWW[256]=37./32.*aryuWW[19] - 2./3.*aryuWW[20];
   aryuWW[256]=1./3.*aryuWW[256] - 7./32.*aryuWW[38];
   aryuWW[256]=aryuWW[38]*aryuWW[256];
   aryuWW[198]=1./4.*aryuWW[242] + aryuWW[256] + aryuWW[198];
   aryuWW[198]=aryuWW[255]*aryuWW[198];
   aryuWW[242]= - 17./2.*aryuWW[36];
   aryuWW[256]= - 5./2.*aryuWW[6];
   aryuWW[257]=aryuWW[256] - 313./3. + aryuWW[242];
   aryuWW[257]=aryuWW[18]*aryuWW[257];
   aryuWW[258]=11./6.*aryuWW[51] + aryuWW[214] + 11./3.*aryuWW[73];
   aryuWW[257]=aryuWW[258] + 1./27.*aryuWW[257];
   aryuWW[259]= - 5./4.*aryuWW[6] + 19./3. - 17./4.*aryuWW[36];
   aryuWW[259]=aryuWW[11]*aryuWW[259];
   aryuWW[260]=1 + aryuWW[68];
   aryuWW[259]=13./4.*aryuWW[260] + 1./9.*aryuWW[259];
   aryuWW[259]=MMH*aryuWW[259];
   aryuWW[261]=5./4.*aryuWW[6];
   aryuWW[262]=aryuWW[261] - 19./3. + 17./4.*aryuWW[36];
   aryuWW[262]=aryuWW[20]*aryuWW[262];
   aryuWW[257]=1./3.*aryuWW[259] - 7./3.*aryuWW[38] + 1./2.*aryuWW[257]
    + 1./27.*aryuWW[262];
   aryuWW[257]=MMH*aryuWW[257];
   aryuWW[259]= - 11*aryuWW[18] + 1231./18.*aryuWW[19];
   aryuWW[259]= - 179./36.*aryuWW[38] + 1./2.*aryuWW[259] - 2749./9.*
   aryuWW[20];
   aryuWW[259]=aryuWW[38]*aryuWW[259];
   aryuWW[257]=1./6.*aryuWW[259] + aryuWW[257];
   aryuWW[130]=aryuWW[198] + 1./2.*aryuWW[257] + aryuWW[130];
   aryuWW[130]=aryuWW[255]*aryuWW[130];
   aryuWW[198]= - 11./2.*aryuWW[6];
   aryuWW[257]= - 11./2.*aryuWW[36];
   aryuWW[259]=aryuWW[198] - 655./3. + aryuWW[257];
   aryuWW[259]=aryuWW[18]*aryuWW[259];
   aryuWW[262]=11./2.*aryuWW[6];
   aryuWW[263]=aryuWW[262] - 47./3. + 11./2.*aryuWW[36];
   aryuWW[263]=1./18.*aryuWW[263] - aryuWW[37];
   aryuWW[263]=aryuWW[20]*aryuWW[263];
   aryuWW[264]=1./3.*aryuWW[37]*aryuWW[18];
   aryuWW[259]=1./3.*aryuWW[263] + aryuWW[264] + aryuWW[258] + 1./54.*
   aryuWW[259];
   aryuWW[263]=aryuWW[198] + 47./3. + aryuWW[257];
   aryuWW[263]=aryuWW[11]*aryuWW[263];
   aryuWW[263]=13*aryuWW[260] + 1./9.*aryuWW[263];
   aryuWW[265]=aryuWW[37]*aryuWW[11];
   aryuWW[263]=1./2.*aryuWW[263] + aryuWW[265];
   aryuWW[263]=MMH*aryuWW[263];
   aryuWW[266]= - 1./9.*aryuWW[11] - 59./9. + aryuWW[10];
   aryuWW[266]=aryuWW[38]*aryuWW[266];
   aryuWW[259]=1./6.*aryuWW[263] + 1./2.*aryuWW[259] + 1./3.*
   aryuWW[266];
   aryuWW[259]=MMH*aryuWW[259];
   aryuWW[263]= - 49*aryuWW[18];
   aryuWW[266]= - 1021./3.*aryuWW[20] + aryuWW[263] - 2323./8.*
   aryuWW[19];
   aryuWW[266]=1./9.*aryuWW[266] - 173./8.*aryuWW[38];
   aryuWW[266]=aryuWW[38]*aryuWW[266];
   aryuWW[130]=aryuWW[130] + aryuWW[138] + 1./4.*aryuWW[266] + 
   aryuWW[259];
   aryuWW[130]=aryuWW[255]*aryuWW[130];
   aryuWW[138]=9*aryuWW[39];
   aryuWW[266]=8./3. + aryuWW[138];
   aryuWW[267]=7./3.*aryuWW[36];
   aryuWW[268]=3*aryuWW[37];
   aryuWW[266]=aryuWW[268] + 2*aryuWW[266] + aryuWW[267];
   aryuWW[266]=aryuWW[3]*aryuWW[38]*aryuWW[266];
   aryuWW[269]= - 3./4.*aryuWW[37];
   aryuWW[270]=1./2.*aryuWW[36];
   aryuWW[271]=aryuWW[269] + 11./36.*aryuWW[6] + 2./27. + aryuWW[270];
   aryuWW[271]=aryuWW[37]*aryuWW[271];
   aryuWW[272]= - 7./2.*aryuWW[36];
   aryuWW[273]= - 8 + aryuWW[272];
   aryuWW[274]=pow(aryuWW[3],2);
   aryuWW[275]= - MMt*aryuWW[274]*aryuWW[38]*aryuWW[39];
   aryuWW[266]=72*aryuWW[275] + aryuWW[266] + 1./9.*aryuWW[273] + 
   aryuWW[271];
   aryuWW[266]=MMt*aryuWW[266];
   aryuWW[271]= - 1./2.*aryuWW[36];
   aryuWW[273]= - 19./3.*aryuWW[37] + 11./6.*aryuWW[6] - 89./9. + 
   aryuWW[271];
   aryuWW[273]=aryuWW[38]*aryuWW[273];
   aryuWW[276]=pow(aryuWW[38],2);
   aryuWW[277]=aryuWW[3]*aryuWW[276];
   aryuWW[273]=1./2.*aryuWW[273] + 16*aryuWW[277];
   aryuWW[266]=1./3.*aryuWW[273] + aryuWW[266];
   aryuWW[266]=MMt*aryuWW[266];
   aryuWW[266]= - 8./9.*aryuWW[276] + aryuWW[266];
   aryuWW[273]=17./2.*aryuWW[36];
   aryuWW[278]=5./2.*aryuWW[6];
   aryuWW[279]=aryuWW[278] - 11./3. + aryuWW[273];
   aryuWW[280]=aryuWW[37]*aryuWW[279];
   aryuWW[281]=aryuWW[3]*aryuWW[38]*aryuWW[39];
   aryuWW[275]=36*aryuWW[275];
   aryuWW[281]=aryuWW[275] + 1./18.*aryuWW[280] + 9*aryuWW[281];
   aryuWW[281]=MMt*aryuWW[281];
   aryuWW[279]=aryuWW[38]*aryuWW[279];
   aryuWW[279]=1./18.*aryuWW[279] + aryuWW[281];
   aryuWW[279]=aryuWW[255]*MMt*aryuWW[279];
   aryuWW[279]=aryuWW[266] + aryuWW[279];
   aryuWW[279]=aryuWW[34]*aryuWW[255]*aryuWW[279];
   aryuWW[130]=aryuWW[130] + aryuWW[279];
   aryuWW[130]=aryuWW[34]*aryuWW[130];
   aryuWW[135]=aryuWW[137] + 89./12. + aryuWW[135];
   aryuWW[118]=aryuWW[133] + aryuWW[166] + aryuWW[127] + aryuWW[164] + 
   aryuWW[142] + aryuWW[118] + aryuWW[162] + aryuWW[168] + aryuWW[126]
    + aryuWW[159] + aryuWW[158] + aryuWW[157] + aryuWW[156] + 
   aryuWW[155] + aryuWW[154] + aryuWW[153] + aryuWW[152] + aryuWW[151]
    + aryuWW[150] + aryuWW[149] + aryuWW[148] + aryuWW[147] + 
   aryuWW[167] + aryuWW[146] + aryuWW[145] + aryuWW[144] + 1./3.*
   aryuWW[135] + aryuWW[143];
   aryuWW[118]=MMt*aryuWW[118];
   aryuWW[117]= - 1771./12. + aryuWW[117];
   aryuWW[116]=aryuWW[180] + aryuWW[178] + aryuWW[176] + aryuWW[175] + 
   aryuWW[121] + aryuWW[116] + aryuWW[174] + aryuWW[182] + aryuWW[119]
    + aryuWW[173] + aryuWW[172] + aryuWW[171] + aryuWW[181] + 1./8.*
   aryuWW[117] + aryuWW[170];
   aryuWW[116]=aryuWW[38]*aryuWW[116];
   aryuWW[117]=aryuWW[193] + 4175./12.*aryuWW[42];
   aryuWW[117]=1./8.*aryuWW[117] + aryuWW[200];
   aryuWW[119]=1733 + 2129*aryuWW[36];
   aryuWW[119]=1./27.*aryuWW[119] + aryuWW[202];
   aryuWW[119]=aryuWW[20]*aryuWW[119];
   aryuWW[121]=5111./16. - 203*aryuWW[36];
   aryuWW[121]=aryuWW[19]*aryuWW[121];
   aryuWW[116]=aryuWW[118] + aryuWW[184] + aryuWW[205] + aryuWW[116] + 
   1./4.*aryuWW[119] + aryuWW[213] + 1./36.*aryuWW[121] + aryuWW[204]
    - 253./72.*aryuWW[23] + aryuWW[211] - 6571./576.*aryuWW[22] - 683./
   288.*aryuWW[51] + aryuWW[210] + aryuWW[209] + aryuWW[208] + 1./3.*
   aryuWW[117] + aryuWW[214];
   aryuWW[116]=MMt*aryuWW[116];
   aryuWW[117]=aryuWW[34]*aryuWW[266];
   aryuWW[118]=2539./8.*aryuWW[38] + 2819./3.*aryuWW[20] + aryuWW[263]
    + 2797./8.*aryuWW[19];
   aryuWW[118]=aryuWW[38]*aryuWW[118];
   aryuWW[116]=aryuWW[117] + aryuWW[116] + 1./36.*aryuWW[118] + 
   aryuWW[259];
   aryuWW[116]=aryuWW[34]*aryuWW[116];
   aryuWW[117]= - 3*aryuWW[47];
   aryuWW[118]=2015./36. + aryuWW[117];
   aryuWW[119]=11./2.*aryuWW[13];
   aryuWW[118]=aryuWW[119] + aryuWW[124] + 1./4.*aryuWW[118] + 
   aryuWW[125];
   aryuWW[121]=1./3. - aryuWW[7];
   aryuWW[126]=1./3.*aryuWW[1]*aryuWW[121];
   aryuWW[127]=aryuWW[40]*aryuWW[121];
   aryuWW[133]= - 1./2.*aryuWW[11];
   aryuWW[118]=aryuWW[133] + aryuWW[218] + aryuWW[127] + aryuWW[126] + 
   59./8.*aryuWW[12] + 1./2.*aryuWW[118] + aryuWW[10];
   aryuWW[118]=aryuWW[37]*aryuWW[118];
   aryuWW[135]= - 3./2.*aryuWW[37];
   aryuWW[137]=aryuWW[223] + aryuWW[135];
   aryuWW[137]=aryuWW[20]*aryuWW[137];
   aryuWW[142]= - 3*aryuWW[19];
   aryuWW[143]= - aryuWW[18] + aryuWW[142];
   aryuWW[143]=aryuWW[37]*aryuWW[143];
   aryuWW[132]=aryuWW[134] + 3*aryuWW[137] + aryuWW[132] + 1./2.*
   aryuWW[143];
   aryuWW[132]=aryuWW[3]*aryuWW[132];
   aryuWW[134]=aryuWW[56] - aryuWW[54];
   aryuWW[134]=aryuWW[227] + aryuWW[136] + 1./2.*aryuWW[134] - 
   aryuWW[55];
   aryuWW[134]=MMt*aryuWW[134];
   aryuWW[136]= - 5./8.*aryuWW[66] + 5*aryuWW[67] - 947./24. + 13*
   aryuWW[61];
   aryuWW[136]=aryuWW[167] + 1./2.*aryuWW[136] + aryuWW[63];
   aryuWW[137]= - 1./4.*aryuWW[69];
   aryuWW[118]=aryuWW[134] + aryuWW[165] + aryuWW[132] + aryuWW[236] + 
   aryuWW[140] + aryuWW[118] + aryuWW[161] + aryuWW[235] + aryuWW[220]
    + aryuWW[131] + aryuWW[234] - 17./4.*aryuWW[17] + aryuWW[233] + 
   aryuWW[237] + aryuWW[232] + aryuWW[231] - 71./32.*aryuWW[50] + 195./
   32.*aryuWW[71] + aryuWW[137] + 1./2.*aryuWW[30] + 3./8.*aryuWW[33]
    + 1./2.*aryuWW[136] + aryuWW[230];
   aryuWW[118]=MMt*aryuWW[118];
   aryuWW[132]=2585./36. + aryuWW[117];
   aryuWW[134]=11./4.*aryuWW[13];
   aryuWW[126]=aryuWW[179] + aryuWW[240] - 13./32.*aryuWW[37] + 3./2.*
   aryuWW[11] + aryuWW[127] + aryuWW[126] + 93./8.*aryuWW[12] - 
   aryuWW[10] + aryuWW[134] + aryuWW[239] + aryuWW[238] + aryuWW[228]
    + aryuWW[221] + 1./8.*aryuWW[132] + aryuWW[194];
   aryuWW[126]=aryuWW[38]*aryuWW[126];
   aryuWW[127]= - 1567./36. + aryuWW[183];
   aryuWW[127]= - 19./9.*aryuWW[36] + aryuWW[160] + aryuWW[189] + 
   aryuWW[125] + aryuWW[188] + aryuWW[187] + aryuWW[186] + aryuWW[185]
    + 1./4.*aryuWW[127] + aryuWW[191];
   aryuWW[132]=1./3. + 3./2.*aryuWW[45];
   aryuWW[132]=5./12.*aryuWW[11] + 1./2.*aryuWW[132] + aryuWW[196];
   aryuWW[132]=aryuWW[37]*aryuWW[132];
   aryuWW[136]= - 19./18.*aryuWW[36] - 20./9. + aryuWW[243];
   aryuWW[136]=aryuWW[11]*aryuWW[136];
   aryuWW[127]=1./2.*aryuWW[132] + aryuWW[136] + 1./4.*aryuWW[127] + 
   aryuWW[199];
   aryuWW[127]=MMH*aryuWW[127];
   aryuWW[132]= - 4*aryuWW[19];
   aryuWW[136]=3./2.*aryuWW[18];
   aryuWW[140]= - 9./2.*aryuWW[38] + 7./2.*aryuWW[20] + aryuWW[136] + 
   aryuWW[132];
   aryuWW[140]=aryuWW[3]*aryuWW[38]*aryuWW[140];
   aryuWW[143]= - 13*aryuWW[74] + 45./8.*aryuWW[42];
   aryuWW[144]=19./3.*aryuWW[36] + 233./12. + aryuWW[160];
   aryuWW[144]=aryuWW[18]*aryuWW[144];
   aryuWW[145]= - 257./8. - 113*aryuWW[36];
   aryuWW[145]=aryuWW[19]*aryuWW[145];
   aryuWW[146]=aryuWW[244] + 743./12.*aryuWW[19];
   aryuWW[146]=aryuWW[37]*aryuWW[146];
   aryuWW[147]= - 197 - 49./2.*aryuWW[36];
   aryuWW[147]=1./3.*aryuWW[147] + 23./2.*aryuWW[37];
   aryuWW[147]=aryuWW[20]*aryuWW[147];
   aryuWW[118]=aryuWW[118] + aryuWW[127] + aryuWW[140] + aryuWW[126] + 
   1./12.*aryuWW[147] + 1./16.*aryuWW[146] + 1./24.*aryuWW[145] + 1./4.
   *aryuWW[144] + 7./8.*aryuWW[23] + aryuWW[253] + 159./64.*aryuWW[22]
    + 363./32.*aryuWW[51] + 17./8.*aryuWW[43] - 1./4.*aryuWW[8] + 
   aryuWW[252] + aryuWW[251] + 1./4.*aryuWW[143] - 2*aryuWW[72];
   aryuWW[118]=MMt*aryuWW[118];
   aryuWW[126]= - 1./2.*aryuWW[6];
   aryuWW[127]=aryuWW[126] - 19 + aryuWW[270];
   aryuWW[127]=aryuWW[18]*aryuWW[127];
   aryuWW[140]= - 1 - aryuWW[36];
   aryuWW[143]=aryuWW[140] + aryuWW[6];
   aryuWW[143]=1./6.*aryuWW[143] - aryuWW[37];
   aryuWW[143]=aryuWW[20]*aryuWW[143];
   aryuWW[127]=1./3.*aryuWW[143] + aryuWW[264] + 1./2.*aryuWW[258] + 1./
   9.*aryuWW[127];
   aryuWW[143]=1 + aryuWW[36];
   aryuWW[144]=aryuWW[143] - aryuWW[6];
   aryuWW[144]=aryuWW[11]*aryuWW[144];
   aryuWW[144]=13./2.*aryuWW[260] + 1./3.*aryuWW[144];
   aryuWW[144]=1./2.*aryuWW[144] + aryuWW[265];
   aryuWW[144]=MMH*aryuWW[144];
   aryuWW[145]= - 11./3.*aryuWW[11] - 29./6. + aryuWW[10];
   aryuWW[145]=aryuWW[38]*aryuWW[145];
   aryuWW[127]=1./6.*aryuWW[144] + 1./2.*aryuWW[127] + 1./3.*
   aryuWW[145];
   aryuWW[127]=MMH*aryuWW[127];
   aryuWW[144]=1./3.*aryuWW[6];
   aryuWW[145]=aryuWW[135] + aryuWW[144] - 3 + aryuWW[272];
   aryuWW[145]=aryuWW[37]*aryuWW[145];
   aryuWW[146]= - 19*aryuWW[36];
   aryuWW[147]=aryuWW[268] + aryuWW[146] - 16 + aryuWW[138];
   aryuWW[147]=aryuWW[3]*aryuWW[38]*aryuWW[147];
   aryuWW[148]=8 + 19./2.*aryuWW[36];
   aryuWW[145]=aryuWW[275] + aryuWW[147] + 1./3.*aryuWW[148] + 1./2.*
   aryuWW[145];
   aryuWW[145]=MMt*aryuWW[145];
   aryuWW[147]= - 17./6.*aryuWW[37] + 1./6.*aryuWW[6] + 17./3. + 3*
   aryuWW[36];
   aryuWW[147]=aryuWW[38]*aryuWW[147];
   aryuWW[148]= - aryuWW[3]*aryuWW[276];
   aryuWW[145]=aryuWW[145] + aryuWW[147] + 16*aryuWW[148];
   aryuWW[145]=MMt*aryuWW[145];
   aryuWW[145]=8./3.*aryuWW[276] + aryuWW[145];
   aryuWW[145]=aryuWW[34]*aryuWW[145];
   aryuWW[147]=7*aryuWW[18];
   aryuWW[148]=aryuWW[147] - 223./12.*aryuWW[19];
   aryuWW[148]= - 27./32.*aryuWW[38] + 1./8.*aryuWW[148] - 20./9.*
   aryuWW[20];
   aryuWW[148]=aryuWW[38]*aryuWW[148];
   aryuWW[118]=aryuWW[145] + aryuWW[118] + aryuWW[148] + aryuWW[127];
   aryuWW[118]=aryuWW[34]*aryuWW[118];
   aryuWW[127]=1./4.*aryuWW[109];
   aryuWW[145]=1./2.*aryuWW[85];
   aryuWW[148]=1./2.*aryuWW[82] + aryuWW[145] + aryuWW[81] + 
   aryuWW[127];
   aryuWW[148]=MMH*aryuWW[148];
   aryuWW[149]= - 1025./18. + 7*aryuWW[95];
   aryuWW[149]= - 13./6.*aryuWW[98] - 13./6.*aryuWW[92] - 13./6.*
   aryuWW[93] + 1./6.*aryuWW[97] - 1./3.*aryuWW[99] + 13./9.*
   aryuWW[101] + 1./4.*aryuWW[149] - 1./3.*aryuWW[94];
   aryuWW[150]=1./3.*aryuWW[91];
   aryuWW[151]= - 11./2.*aryuWW[10] - aryuWW[12];
   aryuWW[151]=aryuWW[12]*aryuWW[151];
   aryuWW[152]= - 3./16.*aryuWW[46];
   aryuWW[153]=13./36.*aryuWW[11] - 59./48.*aryuWW[12] - 5./18.*
   aryuWW[10] - 11./24.*aryuWW[13] + aryuWW[152] - 3./8.*aryuWW[48] - 
   53./27. - 9./16.*aryuWW[44];
   aryuWW[153]=aryuWW[11]*aryuWW[153];
   aryuWW[154]=9./16.*aryuWW[44];
   aryuWW[155]=1./8.*aryuWW[48];
   aryuWW[156]=1./16.*aryuWW[46];
   aryuWW[148]=1./6.*aryuWW[148] + aryuWW[153] + 1./24.*aryuWW[151] + 5.
   /72.*aryuWW[10] + aryuWW[156] - 17./48.*aryuWW[16] + aryuWW[155] + 5.
   /4.*aryuWW[17] + aryuWW[154] - 85./12.*aryuWW[100] + 1./12.*
   aryuWW[96] + 1./4.*aryuWW[149] + aryuWW[150];
   aryuWW[148]=MMH*aryuWW[148];
   aryuWW[149]=9./2.*aryuWW[44];
   aryuWW[151]= - 1./2.*aryuWW[46];
   aryuWW[153]= - 67./3.*aryuWW[12] - 11./3.*aryuWW[13] + aryuWW[151]
    - aryuWW[48] + 1967./27. + aryuWW[149];
   aryuWW[153]=aryuWW[18]*aryuWW[153];
   aryuWW[157]=1./2.*aryuWW[46];
   aryuWW[158]=31./4.*aryuWW[12] + 59./9.*aryuWW[10] + 11./6.*
   aryuWW[13] + aryuWW[157] + 2765./108. + aryuWW[48];
   aryuWW[159]=43./9.*aryuWW[11];
   aryuWW[158]=1./2.*aryuWW[158] + aryuWW[159];
   aryuWW[158]=aryuWW[20]*aryuWW[158];
   aryuWW[161]=43./12.*aryuWW[79] + 5./2.*aryuWW[80] - 169./6.*
   aryuWW[106] - 3./2.*aryuWW[77] + 7./9.*aryuWW[105];
   aryuWW[161]=25./12.*aryuWW[78] + 1./2.*aryuWW[161] - 1./3.*
   aryuWW[107];
   aryuWW[162]=41./8.*aryuWW[12] + 85./3. - 19./4.*aryuWW[10];
   aryuWW[162]=aryuWW[19]*aryuWW[162];
   aryuWW[164]=1./4.*aryuWW[18] + 34./3.*aryuWW[19];
   aryuWW[164]=aryuWW[11]*aryuWW[164];
   aryuWW[148]=aryuWW[148] + 1./2.*aryuWW[158] + 1./3.*aryuWW[164] + 1./
   6.*aryuWW[162] + 1./8.*aryuWW[153] + 143./36.*aryuWW[23] - 127./144.
   *aryuWW[21] + 1./2.*aryuWW[161] + 14./3.*aryuWW[22];
   aryuWW[148]=MMH*aryuWW[148];
   aryuWW[153]=aryuWW[147] + 2759./4.*aryuWW[19];
   aryuWW[153]=aryuWW[19]*aryuWW[153];
   aryuWW[158]= - 5*aryuWW[18];
   aryuWW[161]=3725./2.*aryuWW[20] + aryuWW[158] + 18733./4.*aryuWW[19]
   ;
   aryuWW[161]=aryuWW[20]*aryuWW[161];
   aryuWW[162]=pow(aryuWW[18],2);
   aryuWW[153]=1./3.*aryuWW[161] - 1./3.*aryuWW[162] + aryuWW[153];
   aryuWW[148]=1./12.*aryuWW[153] + aryuWW[148];
   aryuWW[148]=aryuWW[76]*aryuWW[148];
   aryuWW[153]= - 3*aryuWW[10];
   aryuWW[161]=59./9. - 33./2.*aryuWW[13];
   aryuWW[161]=1./2.*aryuWW[161];
   aryuWW[164]=33./2.*aryuWW[12];
   aryuWW[165]=aryuWW[164] + aryuWW[161] + aryuWW[153];
   aryuWW[166]=1./6. - aryuWW[7];
   aryuWW[167]=aryuWW[166] + aryuWW[115];
   aryuWW[168]=1./3.*aryuWW[1]*aryuWW[167];
   aryuWW[167]=aryuWW[40]*aryuWW[167];
   aryuWW[165]=aryuWW[11] + aryuWW[167] + 1./2.*aryuWW[165] + 
   aryuWW[168];
   aryuWW[165]=aryuWW[38]*aryuWW[165];
   aryuWW[161]=aryuWW[164] + aryuWW[161] + aryuWW[10];
   aryuWW[161]= - aryuWW[11] + aryuWW[167] + 1./2.*aryuWW[161] + 
   aryuWW[168];
   aryuWW[161]=aryuWW[37]*aryuWW[161];
   aryuWW[170]= - aryuWW[37]*aryuWW[19];
   aryuWW[171]= - aryuWW[20]*aryuWW[37];
   aryuWW[174]=1./2.*aryuWW[170] + aryuWW[171];
   aryuWW[174]=aryuWW[3]*aryuWW[174];
   aryuWW[174]=aryuWW[161] + 3./2.*aryuWW[174];
   aryuWW[174]=MMt*aryuWW[174];
   aryuWW[176]= - 5./2.*aryuWW[36];
   aryuWW[178]= - 17./3. + aryuWW[176];
   aryuWW[178]=aryuWW[19]*aryuWW[178];
   aryuWW[179]=aryuWW[18]*aryuWW[36];
   aryuWW[178]=aryuWW[179] + aryuWW[178];
   aryuWW[180]= - 1 + aryuWW[271];
   aryuWW[183]=aryuWW[11]*aryuWW[180];
   aryuWW[184]=1./2. - aryuWW[10];
   aryuWW[185]=1./2.*aryuWW[184] + aryuWW[11];
   aryuWW[185]=aryuWW[37]*aryuWW[185];
   aryuWW[183]=aryuWW[185] + aryuWW[183] - 1./4.*aryuWW[36] + 
   aryuWW[10];
   aryuWW[183]=MMH*aryuWW[183];
   aryuWW[186]=17 - 13./2.*aryuWW[36];
   aryuWW[186]=1./2.*aryuWW[186] - aryuWW[37];
   aryuWW[186]=aryuWW[20]*aryuWW[186];
   aryuWW[187]=31./4.*aryuWW[19] - 10*aryuWW[20];
   aryuWW[187]=aryuWW[3]*aryuWW[38]*aryuWW[187];
   aryuWW[188]= - 1./4.*aryuWW[18];
   aryuWW[189]=aryuWW[188] + 4./3.*aryuWW[19];
   aryuWW[189]=aryuWW[37]*aryuWW[189];
   aryuWW[165]=aryuWW[174] + 1./3.*aryuWW[183] + aryuWW[187] + 
   aryuWW[165] + 1./6.*aryuWW[186] + 1./4.*aryuWW[178] + aryuWW[189];
   aryuWW[165]=MMt*aryuWW[165];
   aryuWW[174]=3./2.*aryuWW[36];
   aryuWW[178]= - 3*aryuWW[37];
   aryuWW[183]=aryuWW[178] + aryuWW[115] + 1./3. + aryuWW[174];
   aryuWW[183]=aryuWW[38]*aryuWW[183];
   aryuWW[186]=1./4.*aryuWW[6];
   aryuWW[187]=aryuWW[269] - 1./3. + aryuWW[186];
   aryuWW[187]=aryuWW[37]*aryuWW[187];
   aryuWW[193]= - aryuWW[36] + aryuWW[37];
   aryuWW[193]=aryuWW[3]*aryuWW[38]*aryuWW[193];
   aryuWW[187]=3*aryuWW[193] + aryuWW[270] + aryuWW[187];
   aryuWW[187]=MMt*aryuWW[187];
   aryuWW[183]=1./2.*aryuWW[183] + aryuWW[187];
   aryuWW[183]=aryuWW[34]*MMt*aryuWW[183];
   aryuWW[187]= - aryuWW[36] - aryuWW[6];
   aryuWW[193]=aryuWW[18]*aryuWW[187];
   aryuWW[200]=aryuWW[18] + aryuWW[19];
   aryuWW[202]=aryuWW[37]*aryuWW[200];
   aryuWW[204]= - 1./3.*aryuWW[19];
   aryuWW[202]=aryuWW[202] + 1./2.*aryuWW[193] + aryuWW[204];
   aryuWW[205]=aryuWW[115] + 1./3. + aryuWW[270];
   aryuWW[208]=1./2.*aryuWW[205] - aryuWW[37];
   aryuWW[208]=aryuWW[20]*aryuWW[208];
   aryuWW[202]=1./2.*aryuWW[202] + aryuWW[208];
   aryuWW[208]= - 1./3. + aryuWW[271];
   aryuWW[209]=aryuWW[208] + aryuWW[126];
   aryuWW[210]=aryuWW[11]*aryuWW[209];
   aryuWW[211]=aryuWW[199] + aryuWW[210];
   aryuWW[213]=aryuWW[182] + aryuWW[11];
   aryuWW[214]=aryuWW[37]*aryuWW[213];
   aryuWW[211]=1./2.*aryuWW[211] + aryuWW[214];
   aryuWW[211]=MMH*aryuWW[211];
   aryuWW[214]=aryuWW[10] - aryuWW[11];
   aryuWW[214]=aryuWW[38]*aryuWW[214];
   aryuWW[202]=1./2.*aryuWW[211] + 1./2.*aryuWW[202] + aryuWW[214];
   aryuWW[202]=MMH*aryuWW[202];
   aryuWW[211]= - aryuWW[19] + aryuWW[20];
   aryuWW[211]=aryuWW[38]*aryuWW[211];
   aryuWW[202]=17./4.*aryuWW[211] + aryuWW[202];
   aryuWW[165]=aryuWW[183] + 1./3.*aryuWW[202] + aryuWW[165];
   aryuWW[165]=aryuWW[34]*aryuWW[165];
   aryuWW[183]=aryuWW[7] - aryuWW[6];
   aryuWW[202]=aryuWW[1]*aryuWW[183];
   aryuWW[211]=1./3.*aryuWW[202];
   aryuWW[183]=aryuWW[40]*aryuWW[183];
   aryuWW[214]=aryuWW[211] + aryuWW[183];
   aryuWW[214]=aryuWW[18]*aryuWW[214];
   aryuWW[220]= - 1./3. + aryuWW[7];
   aryuWW[221]=aryuWW[1]*aryuWW[220];
   aryuWW[223]=aryuWW[40]*aryuWW[220];
   aryuWW[221]=1./3.*aryuWW[221] + aryuWW[223];
   aryuWW[223]=aryuWW[19]*aryuWW[221];
   aryuWW[223]=aryuWW[214] + aryuWW[223];
   aryuWW[167]=aryuWW[168] + aryuWW[167];
   aryuWW[167]=aryuWW[20]*aryuWW[167];
   aryuWW[168]=aryuWW[10]*aryuWW[121];
   aryuWW[227]=aryuWW[1]*aryuWW[168];
   aryuWW[228]=aryuWW[40]*aryuWW[168];
   aryuWW[227]=1./3.*aryuWW[227] + aryuWW[228];
   aryuWW[228]=aryuWW[126] - 1./6. + aryuWW[7];
   aryuWW[231]=aryuWW[1]*aryuWW[228];
   aryuWW[228]=aryuWW[40]*aryuWW[228];
   aryuWW[228]=1./3.*aryuWW[231] + aryuWW[228];
   aryuWW[228]=aryuWW[11]*aryuWW[228];
   aryuWW[227]=1./2.*aryuWW[227] + aryuWW[228];
   aryuWW[227]=MMH*aryuWW[227];
   aryuWW[167]=aryuWW[227] + 1./2.*aryuWW[223] + aryuWW[167];
   aryuWW[167]=MMH*aryuWW[167];
   aryuWW[223]= - 11./4.*aryuWW[12];
   aryuWW[134]=aryuWW[223] + aryuWW[134] + aryuWW[199];
   aryuWW[134]=aryuWW[18]*aryuWW[134];
   aryuWW[227]= - 1./9. - 13./8.*aryuWW[10];
   aryuWW[228]= - 11./16.*aryuWW[12];
   aryuWW[227]=1./3.*aryuWW[227] + aryuWW[228];
   aryuWW[227]=aryuWW[19]*aryuWW[227];
   aryuWW[231]= - 47./27. + aryuWW[119];
   aryuWW[223]=5./9.*aryuWW[11] + aryuWW[223] + 1./4.*aryuWW[231] - 5./
   9.*aryuWW[10];
   aryuWW[223]=aryuWW[11]*aryuWW[223];
   aryuWW[231]=aryuWW[12]*aryuWW[10];
   aryuWW[232]=47./27.*aryuWW[10] + 11./2.*aryuWW[231];
   aryuWW[223]=1./4.*aryuWW[232] + aryuWW[223];
   aryuWW[223]=MMH*aryuWW[223];
   aryuWW[232]= - aryuWW[18] + 79./6.*aryuWW[19];
   aryuWW[232]=aryuWW[11]*aryuWW[232];
   aryuWW[234]= - 1./36.*aryuWW[11] + 11./8.*aryuWW[12] - 19./36.*
   aryuWW[10] + 1./27. - 11./16.*aryuWW[13];
   aryuWW[234]=aryuWW[20]*aryuWW[234];
   aryuWW[134]=1./2.*aryuWW[223] + aryuWW[234] + 1./12.*aryuWW[232] + 1.
   /4.*aryuWW[134] + aryuWW[227];
   aryuWW[134]=MMH*aryuWW[134];
   aryuWW[223]= - 17*aryuWW[18] + 91./2.*aryuWW[19];
   aryuWW[223]=aryuWW[19]*aryuWW[223];
   aryuWW[227]= - 257./48.*aryuWW[20] + 17./8.*aryuWW[18] + aryuWW[204]
   ;
   aryuWW[227]=aryuWW[20]*aryuWW[227];
   aryuWW[223]=1./8.*aryuWW[223] + aryuWW[227];
   aryuWW[134]=1./3.*aryuWW[223] + aryuWW[134];
   aryuWW[134]=aryuWW[76]*aryuWW[134];
   aryuWW[134]=aryuWW[165] + 1./6.*aryuWW[167] + aryuWW[134];
   aryuWW[165]=pow(aryuWW[5],2);
   aryuWW[134]=aryuWW[165]*aryuWW[134];
   aryuWW[167]= - aryuWW[1]*aryuWW[14];
   aryuWW[223]= - aryuWW[40]*aryuWW[14];
   aryuWW[167]=1./3.*aryuWW[167] + aryuWW[223];
   aryuWW[223]=1./2.*aryuWW[167];
   aryuWW[227]=1./6. + aryuWW[7];
   aryuWW[232]=aryuWW[1]*aryuWW[227];
   aryuWW[234]=aryuWW[40]*aryuWW[227];
   aryuWW[232]=1./3.*aryuWW[232] + aryuWW[234];
   aryuWW[232]=aryuWW[11]*aryuWW[232];
   aryuWW[232]=aryuWW[223] + aryuWW[232];
   aryuWW[232]=MMH*aryuWW[232];
   aryuWW[221]=aryuWW[18]*aryuWW[221];
   aryuWW[234]=aryuWW[1]*aryuWW[114];
   aryuWW[114]=aryuWW[40]*aryuWW[114];
   aryuWW[114]=1./3.*aryuWW[234] + aryuWW[114];
   aryuWW[114]=aryuWW[20]*aryuWW[114];
   aryuWW[114]=aryuWW[232] + aryuWW[221] + aryuWW[114];
   aryuWW[114]=MMH*aryuWW[114];
   aryuWW[114]=aryuWW[134] + aryuWW[118] + 1./6.*aryuWW[114] + 
   aryuWW[148];
   aryuWW[114]=aryuWW[165]*aryuWW[114];
   aryuWW[118]= - 2*aryuWW[89];
   aryuWW[134]=aryuWW[118] + aryuWW[86];
   aryuWW[148]= - 1./4.*aryuWW[82] + aryuWW[145] + aryuWW[127] + 
   aryuWW[134] + aryuWW[81];
   aryuWW[148]=MMH*aryuWW[148];
   aryuWW[221]= - 9./2.*aryuWW[44];
   aryuWW[232]= - 3*aryuWW[48];
   aryuWW[234]=aryuWW[232] + 25./9. + aryuWW[221];
   aryuWW[234]= - 233./72.*aryuWW[12] - 13./18.*aryuWW[10] - 19./72.*
   aryuWW[13] + 1./2.*aryuWW[234] - aryuWW[46];
   aryuWW[234]=1./2.*aryuWW[234] + 4./9.*aryuWW[11];
   aryuWW[234]=aryuWW[11]*aryuWW[234];
   aryuWW[235]=2*aryuWW[102];
   aryuWW[236]= - 2315./192. + aryuWW[235];
   aryuWW[240]=43./48.*aryuWW[103];
   aryuWW[243]=5./2.*aryuWW[10];
   aryuWW[244]=aryuWW[243] - aryuWW[12];
   aryuWW[244]=aryuWW[12]*aryuWW[244];
   aryuWW[251]=9./8.*aryuWW[44];
   aryuWW[252]=1./4.*aryuWW[48];
   aryuWW[253]=1./8.*aryuWW[46];
   aryuWW[148]=1./3.*aryuWW[148] + aryuWW[234] + 1./12.*aryuWW[244] - 
   13./72.*aryuWW[10] + aryuWW[253] - 13./8.*aryuWW[16] + aryuWW[252]
    + 41./48.*aryuWW[17] + aryuWW[251] - 143./16.*aryuWW[100] + 11./48.
   *aryuWW[96] + aryuWW[240] + 2./3.*aryuWW[91] + 1./12.*aryuWW[98] - 
   41./48.*aryuWW[92] - 41./48.*aryuWW[93] + 1./24.*aryuWW[97] + 1./12.
   *aryuWW[99] + 13./18.*aryuWW[101] - 1./6.*aryuWW[94] + 1./3.*
   aryuWW[236] + 7./8.*aryuWW[95];
   aryuWW[148]=MMH*aryuWW[148];
   aryuWW[234]=9*aryuWW[44];
   aryuWW[236]=877./9. + aryuWW[234];
   aryuWW[236]= - 175./18.*aryuWW[12] + aryuWW[196] - 19./36.*
   aryuWW[13] - aryuWW[46] + 1./2.*aryuWW[236] - aryuWW[48];
   aryuWW[236]=aryuWW[18]*aryuWW[236];
   aryuWW[244]=8./3.*aryuWW[108];
   aryuWW[258]=43 - 5*aryuWW[10];
   aryuWW[258]=7*aryuWW[258] - aryuWW[12];
   aryuWW[258]=aryuWW[19]*aryuWW[258];
   aryuWW[259]=aryuWW[18] + 53./12.*aryuWW[19];
   aryuWW[259]=aryuWW[11]*aryuWW[259];
   aryuWW[260]=31./6.*aryuWW[11] + 89./18.*aryuWW[12] + 107./18.*
   aryuWW[10] + 19./72.*aryuWW[13] + aryuWW[239] + 37./12. + aryuWW[48]
   ;
   aryuWW[260]=aryuWW[20]*aryuWW[260];
   aryuWW[148]=aryuWW[148] + 1./2.*aryuWW[260] + 1./2.*aryuWW[259] + 1./
   24.*aryuWW[258] + 1./4.*aryuWW[236] + 25./9.*aryuWW[23] - 2./9.*
   aryuWW[21] + 5./3.*aryuWW[22] + 7./4.*aryuWW[78] - 5./12.*
   aryuWW[107] + 119./48.*aryuWW[79] + 41./48.*aryuWW[80] - 43./4.*
   aryuWW[106] + 7./18.*aryuWW[105] + aryuWW[244] - 3./4.*aryuWW[77];
   aryuWW[148]=MMH*aryuWW[148];
   aryuWW[236]=29*aryuWW[18] + 665*aryuWW[19];
   aryuWW[236]=aryuWW[19]*aryuWW[236];
   aryuWW[236]=11./3.*aryuWW[162] + 1./4.*aryuWW[236];
   aryuWW[236]=1./2.*aryuWW[236];
   aryuWW[258]=199./8.*aryuWW[18];
   aryuWW[259]=1057./8.*aryuWW[20] + aryuWW[258] - 370*aryuWW[19];
   aryuWW[259]=aryuWW[20]*aryuWW[259];
   aryuWW[259]=aryuWW[236] + 1./3.*aryuWW[259];
   aryuWW[259]=1./3.*aryuWW[259] + aryuWW[148];
   aryuWW[259]=aryuWW[76]*aryuWW[259];
   aryuWW[260]=7./6. + aryuWW[7];
   aryuWW[260]=1./3.*aryuWW[260] + aryuWW[126];
   aryuWW[260]=aryuWW[1]*aryuWW[260];
   aryuWW[264]= - 11./18.*aryuWW[6];
   aryuWW[265]=aryuWW[264] + 47./54. + aryuWW[7];
   aryuWW[265]=aryuWW[40]*aryuWW[265];
   aryuWW[260]=aryuWW[260] + aryuWW[265];
   aryuWW[260]=aryuWW[11]*aryuWW[260];
   aryuWW[167]=aryuWW[167] + aryuWW[260];
   aryuWW[167]=MMH*aryuWW[167];
   aryuWW[227]=1./3.*aryuWW[227] + aryuWW[126];
   aryuWW[227]=aryuWW[1]*aryuWW[227];
   aryuWW[260]=aryuWW[264] - 7./54. + aryuWW[7];
   aryuWW[260]=aryuWW[40]*aryuWW[260];
   aryuWW[227]=aryuWW[227] + aryuWW[260];
   aryuWW[227]=aryuWW[18]*aryuWW[227];
   aryuWW[260]= - 7./6. - aryuWW[7];
   aryuWW[260]=1./3.*aryuWW[260] + aryuWW[115];
   aryuWW[260]=aryuWW[1]*aryuWW[260];
   aryuWW[120]=aryuWW[120] - 47./54. - aryuWW[7];
   aryuWW[120]=aryuWW[40]*aryuWW[120];
   aryuWW[120]=aryuWW[260] + aryuWW[120];
   aryuWW[120]=aryuWW[20]*aryuWW[120];
   aryuWW[120]=aryuWW[167] + aryuWW[227] + aryuWW[120];
   aryuWW[120]=1./6.*MMH*aryuWW[120];
   aryuWW[114]=aryuWW[114] + aryuWW[116] + aryuWW[120] + aryuWW[259];
   aryuWW[114]=aryuWW[165]*aryuWW[114];
   aryuWW[116]= - 3./2.*aryuWW[48];
   aryuWW[167]= - 9./4.*aryuWW[44];
   aryuWW[227]=aryuWW[196] - 1./36.*aryuWW[13] - 5./4.*aryuWW[46] + 
   aryuWW[116] + 251./27. + aryuWW[167];
   aryuWW[197]=aryuWW[197] + 1./2.*aryuWW[227] - 7./9.*aryuWW[12];
   aryuWW[197]=aryuWW[11]*aryuWW[197];
   aryuWW[227]= - 157./192. + aryuWW[235];
   aryuWW[235]= - 1./12.*aryuWW[10];
   aryuWW[259]=13./2.*aryuWW[10] - aryuWW[12];
   aryuWW[259]=aryuWW[12]*aryuWW[259];
   aryuWW[134]= - 1./2.*aryuWW[82] + 1./4.*aryuWW[85] + 1./8.*
   aryuWW[109] + aryuWW[134] + 1./2.*aryuWW[81];
   aryuWW[134]=MMH*aryuWW[134];
   aryuWW[134]=1./3.*aryuWW[134] + 1./2.*aryuWW[197] + 1./24.*
   aryuWW[259] + aryuWW[235] + aryuWW[156] - 61./48.*aryuWW[16] + 
   aryuWW[155] + 5./48.*aryuWW[17] + aryuWW[154] - 137./48.*aryuWW[100]
    + 7./48.*aryuWW[96] + aryuWW[240] + aryuWW[150] + 1./8.*aryuWW[98]
    - 5./16.*aryuWW[92] - 5./16.*aryuWW[93] + 1./6.*aryuWW[99] + 13./36.
   *aryuWW[101] - 1./12.*aryuWW[94] + 1./3.*aryuWW[227] + 7./16.*
   aryuWW[95];
   aryuWW[134]=MMH*aryuWW[134];
   aryuWW[150]= - 11./36.*aryuWW[12] + aryuWW[235] - 1./144.*aryuWW[13]
    + aryuWW[152] - 1./8.*aryuWW[48] + 92./27. + aryuWW[154];
   aryuWW[150]=aryuWW[18]*aryuWW[150];
   aryuWW[152]=aryuWW[228] + 419./24. - 2*aryuWW[10];
   aryuWW[152]=aryuWW[19]*aryuWW[152];
   aryuWW[197]=aryuWW[212] - 17./6.*aryuWW[19];
   aryuWW[197]=aryuWW[11]*aryuWW[197];
   aryuWW[212]=1./36.*aryuWW[13] + aryuWW[46] - 1667./108. + aryuWW[48]
   ;
   aryuWW[212]=55./36.*aryuWW[11] + 77./144.*aryuWW[12] + 1./4.*
   aryuWW[212] + 4./3.*aryuWW[10];
   aryuWW[212]=aryuWW[20]*aryuWW[212];
   aryuWW[134]=aryuWW[134] + aryuWW[212] + 1./12.*aryuWW[197] + 1./3.*
   aryuWW[152] + aryuWW[150] - 43./36.*aryuWW[23] - 49./144.*aryuWW[21]
    - 2./3.*aryuWW[22] + 17./24.*aryuWW[78] - 1./4.*aryuWW[107] + 19./
   12.*aryuWW[79] + 11./48.*aryuWW[80] - 89./24.*aryuWW[106] + 7./36.*
   aryuWW[105] + aryuWW[244] - 3./8.*aryuWW[77];
   aryuWW[134]=MMH*aryuWW[134];
   aryuWW[150]=aryuWW[18] - 49./6.*aryuWW[19];
   aryuWW[150]=aryuWW[19]*aryuWW[150];
   aryuWW[150]=23./9.*aryuWW[162] + 5./2.*aryuWW[150];
   aryuWW[152]=65*aryuWW[18] - 2957./2.*aryuWW[19];
   aryuWW[152]=1./4.*aryuWW[152] - 421*aryuWW[20];
   aryuWW[152]=aryuWW[20]*aryuWW[152];
   aryuWW[150]=1./2.*aryuWW[150] + 1./9.*aryuWW[152];
   aryuWW[134]=1./2.*aryuWW[150] + aryuWW[134];
   aryuWW[134]=aryuWW[76]*aryuWW[134];
   aryuWW[150]= - 1./2.*aryuWW[97];
   aryuWW[152]=aryuWW[150] + 233./48. - aryuWW[99];
   aryuWW[197]=1./4.*aryuWW[96];
   aryuWW[152]= - 11./12.*aryuWW[17] + 77./36.*aryuWW[100] + 
   aryuWW[197] + 13./4.*aryuWW[103] + 5./4.*aryuWW[92] + 1./3.*
   aryuWW[152] + 5./4.*aryuWW[93];
   aryuWW[212]= - 1./3.*aryuWW[16];
   aryuWW[227]=7./18.*aryuWW[12] - 1./3. - aryuWW[46];
   aryuWW[227]=aryuWW[11]*aryuWW[227];
   aryuWW[228]=MMH*aryuWW[82];
   aryuWW[152]=1./12.*aryuWW[228] + 1./8.*aryuWW[227] + 1./48.*
   aryuWW[231] + 1./4.*aryuWW[152] + aryuWW[212];
   aryuWW[152]=MMH*aryuWW[152];
   aryuWW[227]= - 1./6. - aryuWW[10];
   aryuWW[228]=1./4.*aryuWW[12];
   aryuWW[227]=1./3.*aryuWW[227] + aryuWW[228];
   aryuWW[227]=aryuWW[19]*aryuWW[227];
   aryuWW[231]=1./4.*aryuWW[46];
   aryuWW[235]=1./6.*aryuWW[10];
   aryuWW[240]= - 13./18.*aryuWW[12] + aryuWW[235] - 1./9. + 
   aryuWW[231];
   aryuWW[240]=aryuWW[20]*aryuWW[240];
   aryuWW[244]=43./18.*aryuWW[12] - 47./9. - aryuWW[46];
   aryuWW[244]=aryuWW[18]*aryuWW[244];
   aryuWW[152]=aryuWW[152] + 1./2.*aryuWW[240] + 1./4.*aryuWW[227] + 1./
   8.*aryuWW[244] + 4./9.*aryuWW[23] + 1./8.*aryuWW[21] - 1./18.*
   aryuWW[22] + 1./6.*aryuWW[78] - 1./12.*aryuWW[107] + 1./48.*
   aryuWW[79] + 2./9.*aryuWW[106] - 5./16.*aryuWW[80];
   aryuWW[152]=MMH*aryuWW[152];
   aryuWW[227]=1./3.*aryuWW[18] - 41./8.*aryuWW[19];
   aryuWW[227]=aryuWW[19]*aryuWW[227];
   aryuWW[240]= - aryuWW[18] - 203./8.*aryuWW[19];
   aryuWW[240]=1./3.*aryuWW[240] - 425./16.*aryuWW[20];
   aryuWW[240]=aryuWW[20]*aryuWW[240];
   aryuWW[227]=1./2.*aryuWW[227] + aryuWW[240];
   aryuWW[152]=1./3.*aryuWW[227] + aryuWW[152];
   aryuWW[152]=aryuWW[76]*aryuWW[152];
   aryuWW[227]= - 1./2.*aryuWW[80];
   aryuWW[240]=aryuWW[106] + aryuWW[227];
   aryuWW[244]=1./2.*aryuWW[21];
   aryuWW[259]=aryuWW[240] + aryuWW[244];
   aryuWW[260]= - 1 + aryuWW[228];
   aryuWW[260]=aryuWW[18]*aryuWW[260];
   aryuWW[264]= - 1 - aryuWW[12];
   aryuWW[264]=aryuWW[20]*aryuWW[264];
   aryuWW[265]= - 1./2.*aryuWW[19];
   aryuWW[259]=1./4.*aryuWW[264] + aryuWW[265] + 1./2.*aryuWW[259] + 
   aryuWW[260];
   aryuWW[260]=1./2.*aryuWW[93];
   aryuWW[264]=1./2.*aryuWW[92];
   aryuWW[266]=aryuWW[264] + 1 + aryuWW[260];
   aryuWW[269]= - 1./6.*aryuWW[17];
   aryuWW[266]= - 1./6.*aryuWW[16] + aryuWW[269] + 1./2.*aryuWW[100] + 
   1./3.*aryuWW[266] + 1./2.*aryuWW[103];
   aryuWW[266]=MMH*aryuWW[266];
   aryuWW[259]=1./3.*aryuWW[259] + 1./2.*aryuWW[266];
   aryuWW[259]=MMH*aryuWW[259];
   aryuWW[266]= - 5*aryuWW[19];
   aryuWW[272]= - 29*aryuWW[20];
   aryuWW[275]=aryuWW[266] + aryuWW[272];
   aryuWW[275]=aryuWW[20]*aryuWW[275];
   aryuWW[279]=pow(aryuWW[19],2);
   aryuWW[275]= - aryuWW[279] + 1./2.*aryuWW[275];
   aryuWW[259]=1./6.*aryuWW[275] + aryuWW[259];
   aryuWW[259]=aryuWW[255]*aryuWW[76]*aryuWW[259];
   aryuWW[152]=aryuWW[152] + 1./4.*aryuWW[259];
   aryuWW[152]=aryuWW[255]*aryuWW[152];
   aryuWW[259]=1./3. + aryuWW[126];
   aryuWW[259]=aryuWW[1]*aryuWW[259];
   aryuWW[198]=19./3. + aryuWW[198];
   aryuWW[198]=aryuWW[40]*aryuWW[198];
   aryuWW[198]=aryuWW[259] + 1./9.*aryuWW[198];
   aryuWW[198]=aryuWW[11]*aryuWW[198];
   aryuWW[198]=aryuWW[223] + aryuWW[198];
   aryuWW[198]=MMH*aryuWW[198];
   aryuWW[223]=1./3. - aryuWW[6];
   aryuWW[259]=aryuWW[1]*aryuWW[223];
   aryuWW[223]=aryuWW[40]*aryuWW[223];
   aryuWW[275]=11./9.*aryuWW[223];
   aryuWW[281]=aryuWW[259] + aryuWW[275];
   aryuWW[282]=aryuWW[18]*aryuWW[281];
   aryuWW[283]= - 1./3. + aryuWW[115];
   aryuWW[283]=aryuWW[1]*aryuWW[283];
   aryuWW[262]= - 19./3. + aryuWW[262];
   aryuWW[262]=aryuWW[40]*aryuWW[262];
   aryuWW[262]=aryuWW[283] + 1./9.*aryuWW[262];
   aryuWW[262]=aryuWW[20]*aryuWW[262];
   aryuWW[198]=aryuWW[198] + 1./2.*aryuWW[282] + aryuWW[262];
   aryuWW[198]=MMH*aryuWW[198];
   aryuWW[134]=aryuWW[152] + 1./6.*aryuWW[198] + aryuWW[134];
   aryuWW[134]=aryuWW[255]*aryuWW[134];
   aryuWW[152]=1633./8.*aryuWW[20] + aryuWW[258] + 218*aryuWW[19];
   aryuWW[152]=aryuWW[20]*aryuWW[152];
   aryuWW[152]=aryuWW[236] + 1./3.*aryuWW[152];
   aryuWW[148]=1./3.*aryuWW[152] + aryuWW[148];
   aryuWW[148]=aryuWW[76]*aryuWW[148];
   aryuWW[120]=aryuWW[134] + aryuWW[120] + aryuWW[148];
   aryuWW[120]=aryuWW[255]*aryuWW[120];
   aryuWW[134]=15./2.*aryuWW[41];
   aryuWW[148]= - 27./4.*aryuWW[73];
   aryuWW[152]= - 17*aryuWW[51];
   aryuWW[198]= - 7./2.*aryuWW[21];
   aryuWW[236]= - 39./4.*aryuWW[23] + aryuWW[198] + 85./12.*aryuWW[22]
    + aryuWW[152] + 57./4.*aryuWW[43] + aryuWW[148] + aryuWW[134] - 73./
   4.*aryuWW[74] + aryuWW[72];
   aryuWW[258]=aryuWW[181] + 5./8. + aryuWW[125];
   aryuWW[258]=aryuWW[37]*aryuWW[258];
   aryuWW[262]= - 1./2.*aryuWW[68];
   aryuWW[282]= - aryuWW[37]*aryuWW[18];
   aryuWW[283]=aryuWW[282] - aryuWW[21] + aryuWW[18];
   aryuWW[283]=aryuWW[3]*aryuWW[283];
   aryuWW[284]= - MMH*aryuWW[55];
   aryuWW[258]=aryuWW[284] + 1./2.*aryuWW[283] + aryuWW[258] + 
   aryuWW[262] - 15./8.*aryuWW[70] - 2*aryuWW[50] - 9./2.*aryuWW[71] + 
   5./4.*aryuWW[65] + 3./4.*aryuWW[67] + 9./8. + aryuWW[53];
   aryuWW[258]=MMt*aryuWW[258];
   aryuWW[195]=aryuWW[195] + 7./6.*aryuWW[11];
   aryuWW[195]=aryuWW[37]*aryuWW[195];
   aryuWW[285]= - 7./2. - 5*aryuWW[65];
   aryuWW[285]= - 9./2.*aryuWW[16] - 11./4.*aryuWW[68] + 9./4.*
   aryuWW[70] + 5./4.*aryuWW[50] + 3./2.*aryuWW[60] + 1./4.*aryuWW[285]
    + aryuWW[191];
   aryuWW[286]= - 1./2. + aryuWW[131];
   aryuWW[286]=aryuWW[11]*aryuWW[286];
   aryuWW[285]=1./2.*aryuWW[285] + aryuWW[286];
   aryuWW[195]=aryuWW[285] + 1./2.*aryuWW[195];
   aryuWW[195]=MMH*aryuWW[195];
   aryuWW[286]= - 1./4.*aryuWW[49];
   aryuWW[287]=aryuWW[286] + 13./8. - aryuWW[45];
   aryuWW[288]= - 5./4.*aryuWW[37];
   aryuWW[287]=3*aryuWW[287] + aryuWW[288];
   aryuWW[287]=aryuWW[38]*aryuWW[287];
   aryuWW[289]=1./2.*aryuWW[18];
   aryuWW[290]=aryuWW[289] + 17./3.*aryuWW[19];
   aryuWW[290]=aryuWW[37]*aryuWW[290];
   aryuWW[291]=aryuWW[37] + 41./24. + aryuWW[160];
   aryuWW[291]=aryuWW[20]*aryuWW[291];
   aryuWW[292]=3./4. - aryuWW[39];
   aryuWW[292]=aryuWW[18]*aryuWW[292];
   aryuWW[293]=3*aryuWW[292];
   aryuWW[195]=aryuWW[258] + aryuWW[195] + aryuWW[287] + aryuWW[291] + 
   1./4.*aryuWW[290] + 121./12.*aryuWW[19] + 1./2.*aryuWW[236] + 
   aryuWW[293];
   aryuWW[195]=MMt*aryuWW[195];
   aryuWW[236]=1./3.*aryuWW[37];
   aryuWW[190]=aryuWW[236] + aryuWW[190] - 55./27. + aryuWW[131];
   aryuWW[190]=aryuWW[20]*aryuWW[190];
   aryuWW[258]=3./2.*aryuWW[39];
   aryuWW[287]= - 7./108.*aryuWW[36] - 67./27. + aryuWW[258];
   aryuWW[287]=aryuWW[18]*aryuWW[287];
   aryuWW[290]= - 3./2.*aryuWW[23] + 3./2.*aryuWW[21] + 7./4.*
   aryuWW[51] + aryuWW[225] + 5./2.*aryuWW[73] + aryuWW[72] - 9./4.*
   aryuWW[41];
   aryuWW[282]=1./6.*aryuWW[282];
   aryuWW[190]=1./2.*aryuWW[190] + aryuWW[282] + aryuWW[290] + 
   aryuWW[287];
   aryuWW[287]= - 1./2.*aryuWW[60];
   aryuWW[291]=3./2.*aryuWW[16];
   aryuWW[294]=aryuWW[291] + 7./4.*aryuWW[68] + aryuWW[287] + 7./4. - 
   aryuWW[62];
   aryuWW[295]= - 7./54.*aryuWW[36] + 55./27. + aryuWW[160];
   aryuWW[295]=aryuWW[11]*aryuWW[295];
   aryuWW[296]= - aryuWW[37]*aryuWW[11];
   aryuWW[297]=1./6.*aryuWW[296];
   aryuWW[295]=aryuWW[297] + aryuWW[294] + 1./2.*aryuWW[295];
   aryuWW[295]=MMH*aryuWW[295];
   aryuWW[298]= - 23./3. + aryuWW[194];
   aryuWW[298]=1./2.*aryuWW[298] + aryuWW[196];
   aryuWW[299]= - 2./3.*aryuWW[11];
   aryuWW[298]=1./2.*aryuWW[298] + aryuWW[299];
   aryuWW[298]=aryuWW[38]*aryuWW[298];
   aryuWW[190]=1./2.*aryuWW[295] + 1./2.*aryuWW[190] + aryuWW[298];
   aryuWW[190]=MMH*aryuWW[190];
   aryuWW[295]= - 13*aryuWW[18];
   aryuWW[298]=75*aryuWW[20] + aryuWW[295] - 17*aryuWW[19];
   aryuWW[300]=7*aryuWW[38];
   aryuWW[298]=1./2.*aryuWW[298] + aryuWW[300];
   aryuWW[298]=aryuWW[38]*aryuWW[298];
   aryuWW[190]=aryuWW[195] + 1./4.*aryuWW[298] + aryuWW[190];
   aryuWW[190]=MMt*aryuWW[190];
   aryuWW[195]=aryuWW[181] + 13./8. + aryuWW[125];
   aryuWW[195]=aryuWW[37]*aryuWW[195];
   aryuWW[298]= - 15./16.*aryuWW[70];
   aryuWW[301]= - 1./4.*aryuWW[68];
   aryuWW[302]=1./4.*aryuWW[283];
   aryuWW[303]=1./2.*aryuWW[284];
   aryuWW[195]=aryuWW[303] + aryuWW[302] + 1./2.*aryuWW[195] + 
   aryuWW[301] + aryuWW[298] - aryuWW[50] - 33./16.*aryuWW[71] + 5./8.*
   aryuWW[65] + 3./8.*aryuWW[67] + 5./4. + aryuWW[53];
   aryuWW[195]=MMt*aryuWW[195];
   aryuWW[304]= - 7*aryuWW[23] + aryuWW[198] + 67./36.*aryuWW[22] + 
   aryuWW[152] + 31./2.*aryuWW[43] + aryuWW[148] + aryuWW[134] - 67./4.
   *aryuWW[74] + aryuWW[72];
   aryuWW[305]=aryuWW[136] - 19./9.*aryuWW[19];
   aryuWW[305]=aryuWW[37]*aryuWW[305];
   aryuWW[304]=1./4.*aryuWW[305] + 833./72.*aryuWW[19] + 1./2.*
   aryuWW[304] + aryuWW[293];
   aryuWW[305]=aryuWW[288] + aryuWW[181] + 25./8. + aryuWW[125];
   aryuWW[305]=aryuWW[38]*aryuWW[305];
   aryuWW[194]=aryuWW[194] + aryuWW[11];
   aryuWW[194]=aryuWW[37]*aryuWW[194];
   aryuWW[194]=aryuWW[285] + 1./4.*aryuWW[194];
   aryuWW[194]=MMH*aryuWW[194];
   aryuWW[306]=25./24. + aryuWW[160];
   aryuWW[306]=1./2.*aryuWW[306] + 5./3.*aryuWW[37];
   aryuWW[306]=aryuWW[20]*aryuWW[306];
   aryuWW[194]=aryuWW[195] + 1./2.*aryuWW[194] + 1./2.*aryuWW[305] + 1./
   2.*aryuWW[304] + aryuWW[306];
   aryuWW[194]=MMt*aryuWW[194];
   aryuWW[195]= - 5 + aryuWW[160];
   aryuWW[195]=aryuWW[18]*aryuWW[195];
   aryuWW[304]= - 1 + aryuWW[192];
   aryuWW[304]=aryuWW[20]*aryuWW[304];
   aryuWW[195]=aryuWW[304] + aryuWW[290] + 1./2.*aryuWW[195];
   aryuWW[304]=1 + aryuWW[258];
   aryuWW[304]=aryuWW[11]*aryuWW[304];
   aryuWW[305]=aryuWW[294] + aryuWW[304];
   aryuWW[305]=MMH*aryuWW[305];
   aryuWW[306]= - 1./4.*aryuWW[11];
   aryuWW[307]=aryuWW[306] - 1 + 3./8.*aryuWW[45];
   aryuWW[307]=aryuWW[38]*aryuWW[307];
   aryuWW[195]=1./4.*aryuWW[305] + 1./4.*aryuWW[195] + aryuWW[307];
   aryuWW[195]=MMH*aryuWW[195];
   aryuWW[305]= - 5./3.*aryuWW[19];
   aryuWW[307]= - aryuWW[18] + aryuWW[305];
   aryuWW[307]=7./2.*aryuWW[38] + 7./4.*aryuWW[307] + 59./3.*aryuWW[20]
   ;
   aryuWW[307]=aryuWW[38]*aryuWW[307];
   aryuWW[194]=aryuWW[194] + 1./4.*aryuWW[307] + aryuWW[195];
   aryuWW[194]=MMt*aryuWW[194];
   aryuWW[195]=3./2.*aryuWW[19];
   aryuWW[307]= - 1 - aryuWW[71];
   aryuWW[307]=MMt*aryuWW[307];
   aryuWW[307]=1./2.*aryuWW[307] + aryuWW[38] + aryuWW[195] + 
   aryuWW[226] - aryuWW[74] + aryuWW[225];
   aryuWW[307]=MMt*aryuWW[307];
   aryuWW[308]=aryuWW[38]*aryuWW[20];
   aryuWW[307]=1./2.*aryuWW[308] + aryuWW[307];
   aryuWW[307]=aryuWW[255]*MMt*aryuWW[307];
   aryuWW[194]=aryuWW[194] + 1./4.*aryuWW[307];
   aryuWW[194]=aryuWW[255]*aryuWW[194];
   aryuWW[307]=aryuWW[18] - aryuWW[20];
   aryuWW[307]=aryuWW[38]*aryuWW[307];
   aryuWW[309]=MMH*aryuWW[38]*aryuWW[11];
   aryuWW[307]=aryuWW[307] + aryuWW[309];
   aryuWW[307]=MMH*aryuWW[307];
   aryuWW[190]=aryuWW[194] + 1./108.*aryuWW[307] + aryuWW[190];
   aryuWW[190]=aryuWW[255]*aryuWW[190];
   aryuWW[194]= - 9*aryuWW[45];
   aryuWW[307]= - 9./4.*aryuWW[49];
   aryuWW[309]=aryuWW[307] + 7./8. + aryuWW[194];
   aryuWW[309]=aryuWW[37]*aryuWW[309];
   aryuWW[283]=3./2.*aryuWW[284] + 3./4.*aryuWW[283] + 1./2.*
   aryuWW[309] - 3./4.*aryuWW[68] - 45./16.*aryuWW[70] - 3*aryuWW[50]
    - 111./16.*aryuWW[71] + 15./8.*aryuWW[65] + 9./8.*aryuWW[67] + 1 + 
   aryuWW[53];
   aryuWW[283]=MMt*aryuWW[283];
   aryuWW[309]=aryuWW[148] + aryuWW[134] - 75./4.*aryuWW[74] + 
   aryuWW[72];
   aryuWW[309]= - 21./2.*aryuWW[21] + 395./12.*aryuWW[22] - 51*
   aryuWW[51] + 3*aryuWW[309] + 83./2.*aryuWW[43];
   aryuWW[194]= - 15./4.*aryuWW[37] + aryuWW[307] + 131./8. + 
   aryuWW[194];
   aryuWW[194]=aryuWW[38]*aryuWW[194];
   aryuWW[307]=1 + aryuWW[237];
   aryuWW[307]=11./4.*aryuWW[11] + 1./2.*aryuWW[307] - aryuWW[10];
   aryuWW[307]=aryuWW[37]*aryuWW[307];
   aryuWW[307]=3*aryuWW[285] + aryuWW[307];
   aryuWW[307]=MMH*aryuWW[307];
   aryuWW[310]= - 3./2.*aryuWW[18];
   aryuWW[311]=aryuWW[310] + 145./3.*aryuWW[19];
   aryuWW[311]=aryuWW[37]*aryuWW[311];
   aryuWW[312]=13./3.*aryuWW[37] + 139./24. + aryuWW[138];
   aryuWW[312]=aryuWW[20]*aryuWW[312];
   aryuWW[194]=aryuWW[283] + 1./2.*aryuWW[307] + 1./2.*aryuWW[194] + 1./
   2.*aryuWW[312] + 1./8.*aryuWW[311] + 613./48.*aryuWW[19] + 9./2.*
   aryuWW[292] + 1./4.*aryuWW[309] - 8*aryuWW[23];
   aryuWW[194]=MMt*aryuWW[194];
   aryuWW[138]=59./27.*aryuWW[36] - 337./27. + aryuWW[138];
   aryuWW[138]=aryuWW[18]*aryuWW[138];
   aryuWW[283]=1./3.*aryuWW[19];
   aryuWW[292]=aryuWW[158] - aryuWW[19];
   aryuWW[307]=1./6.*aryuWW[37]*aryuWW[292];
   aryuWW[138]=aryuWW[307] + aryuWW[283] + 3*aryuWW[290] + 1./2.*
   aryuWW[138];
   aryuWW[309]=3*aryuWW[294] + aryuWW[196];
   aryuWW[311]=59./216.*aryuWW[36] + 31./27. + 9./8.*aryuWW[39];
   aryuWW[311]=aryuWW[11]*aryuWW[311];
   aryuWW[312]=aryuWW[235] - aryuWW[11];
   aryuWW[312]=aryuWW[37]*aryuWW[312];
   aryuWW[309]=1./4.*aryuWW[312] + 1./4.*aryuWW[309] + aryuWW[311];
   aryuWW[309]=MMH*aryuWW[309];
   aryuWW[311]= - 5*aryuWW[11];
   aryuWW[237]=aryuWW[311] - 11 + aryuWW[237];
   aryuWW[237]=aryuWW[38]*aryuWW[237];
   aryuWW[313]=1./4.*aryuWW[37] - 59./216.*aryuWW[36] - 31./27. - 9./8.
   *aryuWW[39];
   aryuWW[313]=aryuWW[20]*aryuWW[313];
   aryuWW[138]=aryuWW[309] + 1./4.*aryuWW[237] + 1./4.*aryuWW[138] + 
   aryuWW[313];
   aryuWW[138]=MMH*aryuWW[138];
   aryuWW[237]=aryuWW[263] - 43*aryuWW[19];
   aryuWW[237]=21./4.*aryuWW[38] + 1./8.*aryuWW[237] + 97./3.*
   aryuWW[20];
   aryuWW[237]=aryuWW[38]*aryuWW[237];
   aryuWW[138]=aryuWW[194] + 1./2.*aryuWW[237] + aryuWW[138];
   aryuWW[138]=MMt*aryuWW[138];
   aryuWW[194]= - 43./18.*aryuWW[20] + 17./9.*aryuWW[18] + aryuWW[245];
   aryuWW[194]=aryuWW[38]*aryuWW[194];
   aryuWW[159]= - aryuWW[10] + aryuWW[159];
   aryuWW[159]=MMH*aryuWW[38]*aryuWW[159];
   aryuWW[159]=aryuWW[194] + 1./2.*aryuWW[159];
   aryuWW[159]=MMH*aryuWW[159];
   aryuWW[138]=1./6.*aryuWW[159] + aryuWW[138];
   aryuWW[159]=aryuWW[138] + aryuWW[190];
   aryuWW[159]=aryuWW[255]*aryuWW[159];
   aryuWW[190]= - 1 + 7./2.*aryuWW[36];
   aryuWW[194]=1./9.*aryuWW[190] + aryuWW[37];
   aryuWW[194]=aryuWW[37]*aryuWW[194];
   aryuWW[237]=aryuWW[38]*aryuWW[37];
   aryuWW[263]=aryuWW[3]*aryuWW[237];
   aryuWW[309]=3*aryuWW[263];
   aryuWW[194]=1./2.*aryuWW[194] + aryuWW[309];
   aryuWW[194]=MMt*aryuWW[194];
   aryuWW[190]=1./2.*aryuWW[190] + 4*aryuWW[37];
   aryuWW[190]=aryuWW[38]*aryuWW[190];
   aryuWW[313]=3*aryuWW[277];
   aryuWW[190]=aryuWW[194] + 1./9.*aryuWW[190] + aryuWW[313];
   aryuWW[190]=MMt*aryuWW[190];
   aryuWW[190]= - 1./18.*aryuWW[276] + aryuWW[190];
   aryuWW[190]=aryuWW[255]*MMt*aryuWW[190];
   aryuWW[194]= - 43 - 59./2.*aryuWW[36];
   aryuWW[314]=1./9.*aryuWW[194] + aryuWW[268];
   aryuWW[314]=aryuWW[37]*aryuWW[314];
   aryuWW[314]=1./2.*aryuWW[314] + 9*aryuWW[263];
   aryuWW[314]=MMt*aryuWW[314];
   aryuWW[194]=1./2.*aryuWW[194] - 8*aryuWW[37];
   aryuWW[194]=aryuWW[38]*aryuWW[194];
   aryuWW[194]=aryuWW[314] + 1./9.*aryuWW[194] + 9*aryuWW[277];
   aryuWW[194]=MMt*aryuWW[194];
   aryuWW[194]= - 43./18.*aryuWW[276] + aryuWW[194];
   aryuWW[194]=MMt*aryuWW[194];
   aryuWW[190]=aryuWW[194] + aryuWW[190];
   aryuWW[190]=aryuWW[34]*aryuWW[255]*aryuWW[190];
   aryuWW[159]=aryuWW[159] + aryuWW[190];
   aryuWW[159]=aryuWW[34]*aryuWW[159];
   aryuWW[190]=aryuWW[34]*aryuWW[194];
   aryuWW[138]=aryuWW[138] + aryuWW[190];
   aryuWW[138]=aryuWW[34]*aryuWW[138];
   aryuWW[190]=1./36.*aryuWW[10];
   aryuWW[154]=aryuWW[190] + aryuWW[156] + aryuWW[155] - 25./9. + 
   aryuWW[154];
   aryuWW[154]=aryuWW[18]*aryuWW[154];
   aryuWW[155]= - 3*aryuWW[95];
   aryuWW[156]=139./6. + aryuWW[155];
   aryuWW[194]= - 5./3.*aryuWW[97];
   aryuWW[156]=1./2.*aryuWW[156] + aryuWW[194];
   aryuWW[314]=17./6.*aryuWW[16];
   aryuWW[315]=1./9.*aryuWW[10];
   aryuWW[156]=aryuWW[315] + aryuWW[314] + 17./3.*aryuWW[100] - 1./4.*
   aryuWW[96] + 1./2.*aryuWW[156] - aryuWW[91];
   aryuWW[316]=9./4.*aryuWW[44];
   aryuWW[317]=1./2.*aryuWW[48];
   aryuWW[318]=aryuWW[231] + aryuWW[317] + 1./9. + aryuWW[316];
   aryuWW[318]= - 19./36.*aryuWW[11] + 1./2.*aryuWW[318] + aryuWW[199];
   aryuWW[318]=aryuWW[11]*aryuWW[318];
   aryuWW[319]=aryuWW[84] + 1./3.*aryuWW[81];
   aryuWW[319]=MMH*aryuWW[319];
   aryuWW[320]=1./4.*aryuWW[319];
   aryuWW[156]=aryuWW[320] + 1./2.*aryuWW[156] + aryuWW[318];
   aryuWW[156]=MMH*aryuWW[156];
   aryuWW[318]= - 1./2.*aryuWW[48];
   aryuWW[321]= - 1./4.*aryuWW[46];
   aryuWW[322]= - 5./2.*aryuWW[11] - 61./18.*aryuWW[10] + aryuWW[321]
    + aryuWW[318] - 37./9. + aryuWW[167];
   aryuWW[322]=aryuWW[20]*aryuWW[322];
   aryuWW[323]= - 5*aryuWW[77] + 11./3.*aryuWW[105];
   aryuWW[323]= - 37./12.*aryuWW[23] + 5./4.*aryuWW[21] - 17./6.*
   aryuWW[22] - 37./24.*aryuWW[78] + 1./8.*aryuWW[107] - 19./12.*
   aryuWW[79] + 1./8.*aryuWW[323] + 17./3.*aryuWW[106];
   aryuWW[324]= - 13./3. + 11./8.*aryuWW[10];
   aryuWW[324]=aryuWW[19]*aryuWW[324];
   aryuWW[325]= - aryuWW[11]*aryuWW[19];
   aryuWW[154]=1./2.*aryuWW[156] + 1./4.*aryuWW[322] + 163./72.*
   aryuWW[325] + 1./3.*aryuWW[324] + 1./2.*aryuWW[323] + aryuWW[154];
   aryuWW[154]=MMH*aryuWW[154];
   aryuWW[156]= - aryuWW[18] + aryuWW[19];
   aryuWW[156]=aryuWW[19]*aryuWW[156];
   aryuWW[295]=aryuWW[295] - 121./3.*aryuWW[19];
   aryuWW[295]=1./2.*aryuWW[295] + 83./3.*aryuWW[20];
   aryuWW[295]=aryuWW[20]*aryuWW[295];
   aryuWW[295]=1./6.*aryuWW[295] + 1./3.*aryuWW[162] + 1./2.*
   aryuWW[156];
   aryuWW[154]=1./2.*aryuWW[295] + aryuWW[154];
   aryuWW[154]=aryuWW[76]*MMH*aryuWW[154];
   aryuWW[295]= - 39./2.*aryuWW[71] + 5*aryuWW[65] - 1 + 3*aryuWW[67];
   aryuWW[286]=aryuWW[286] - 1./8. - aryuWW[45];
   aryuWW[286]=aryuWW[37]*aryuWW[286];
   aryuWW[286]=aryuWW[303] + aryuWW[302] + 3./2.*aryuWW[286] + 
   aryuWW[301] + aryuWW[298] + 1./8.*aryuWW[295] - aryuWW[50];
   aryuWW[286]=MMt*aryuWW[286];
   aryuWW[134]= - 25./2.*aryuWW[23] + aryuWW[198] + 75./4.*aryuWW[22]
    + aryuWW[152] + 13*aryuWW[43] + aryuWW[148] + aryuWW[134] - 79./4.*
   aryuWW[74] + aryuWW[72];
   aryuWW[125]=aryuWW[288] + aryuWW[181] + 53./8. + aryuWW[125];
   aryuWW[125]=aryuWW[38]*aryuWW[125];
   aryuWW[148]= - 5./2.*aryuWW[18] + 37*aryuWW[19];
   aryuWW[148]=aryuWW[37]*aryuWW[148];
   aryuWW[152]=7./3.*aryuWW[37] + 19./8. + aryuWW[160];
   aryuWW[152]=aryuWW[20]*aryuWW[152];
   aryuWW[125]=aryuWW[125] + aryuWW[152] + 1./4.*aryuWW[148] + 43./8.*
   aryuWW[19] + 1./2.*aryuWW[134] + aryuWW[293];
   aryuWW[134]=3./4.*aryuWW[45];
   aryuWW[148]=1./3. + aryuWW[134];
   aryuWW[148]=19./24.*aryuWW[11] + 1./2.*aryuWW[148] + aryuWW[196];
   aryuWW[148]=aryuWW[37]*aryuWW[148];
   aryuWW[148]=1./2.*aryuWW[285] + aryuWW[148];
   aryuWW[148]=MMH*aryuWW[148];
   aryuWW[125]=aryuWW[286] + 1./2.*aryuWW[125] + aryuWW[148];
   aryuWW[125]=MMt*aryuWW[125];
   aryuWW[148]= - 23./9. + aryuWW[160];
   aryuWW[152]=11./9.*aryuWW[36];
   aryuWW[148]=1./2.*aryuWW[148] + aryuWW[152];
   aryuWW[148]=aryuWW[18]*aryuWW[148];
   aryuWW[181]= - aryuWW[18] + aryuWW[204];
   aryuWW[181]=aryuWW[37]*aryuWW[181];
   aryuWW[148]=1./2.*aryuWW[181] + aryuWW[283] + aryuWW[290] + 
   aryuWW[148];
   aryuWW[152]=aryuWW[152] + 23./9. + aryuWW[258];
   aryuWW[152]=aryuWW[11]*aryuWW[152];
   aryuWW[152]=aryuWW[152] + aryuWW[294] + aryuWW[196];
   aryuWW[181]=1./4.*aryuWW[10];
   aryuWW[198]=aryuWW[181] - aryuWW[11];
   aryuWW[198]=aryuWW[37]*aryuWW[198];
   aryuWW[152]=1./2.*aryuWW[152] + 1./3.*aryuWW[198];
   aryuWW[152]=MMH*aryuWW[152];
   aryuWW[204]= - 11./9.*aryuWW[36] - 23./9. + aryuWW[192];
   aryuWW[204]=1./2.*aryuWW[204] + aryuWW[236];
   aryuWW[204]=aryuWW[20]*aryuWW[204];
   aryuWW[134]= - 7./6.*aryuWW[11] + aryuWW[199] - 5./3. + aryuWW[134];
   aryuWW[134]=aryuWW[38]*aryuWW[134];
   aryuWW[134]=aryuWW[152] + aryuWW[134] + 1./2.*aryuWW[148] + 
   aryuWW[204];
   aryuWW[134]=MMH*aryuWW[134];
   aryuWW[148]= - 23*aryuWW[18] - 9*aryuWW[19];
   aryuWW[148]=aryuWW[300] + 1./2.*aryuWW[148] + 163./3.*aryuWW[20];
   aryuWW[148]=aryuWW[38]*aryuWW[148];
   aryuWW[134]=1./4.*aryuWW[148] + aryuWW[134];
   aryuWW[125]=1./2.*aryuWW[134] + aryuWW[125];
   aryuWW[125]=MMt*aryuWW[125];
   aryuWW[134]= - 7 + aryuWW[257];
   aryuWW[148]=1./3.*aryuWW[134] + aryuWW[37];
   aryuWW[148]=aryuWW[37]*aryuWW[148];
   aryuWW[148]=aryuWW[148] + 6*aryuWW[263];
   aryuWW[148]=MMt*aryuWW[148];
   aryuWW[134]=aryuWW[134] - 4*aryuWW[37];
   aryuWW[134]=aryuWW[38]*aryuWW[134];
   aryuWW[134]=aryuWW[148] + 1./3.*aryuWW[134] + 6*aryuWW[277];
   aryuWW[134]=MMt*aryuWW[134];
   aryuWW[134]= - 7./3.*aryuWW[276] + aryuWW[134];
   aryuWW[134]=aryuWW[34]*MMt*aryuWW[134];
   aryuWW[148]=11./3.*aryuWW[18] + aryuWW[19];
   aryuWW[148]=1./2.*aryuWW[148] - 7./3.*aryuWW[20];
   aryuWW[148]=aryuWW[38]*aryuWW[148];
   aryuWW[152]=aryuWW[182] + 7./3.*aryuWW[11];
   aryuWW[152]=MMH*aryuWW[38]*aryuWW[152];
   aryuWW[148]=aryuWW[148] + aryuWW[152];
   aryuWW[148]=MMH*aryuWW[148];
   aryuWW[125]=aryuWW[134] + 1./6.*aryuWW[148] + aryuWW[125];
   aryuWW[125]=aryuWW[34]*aryuWW[125];
   aryuWW[134]=aryuWW[180] + aryuWW[37];
   aryuWW[148]=aryuWW[20]*aryuWW[134];
   aryuWW[152]= - aryuWW[18] - aryuWW[19];
   aryuWW[204]=1./2.*aryuWW[37]*aryuWW[152];
   aryuWW[148]=aryuWW[148] + aryuWW[204] + 1./2.*aryuWW[179] + 
   aryuWW[19];
   aryuWW[179]=1 + aryuWW[270];
   aryuWW[236]=aryuWW[11]*aryuWW[179];
   aryuWW[257]=aryuWW[241] - aryuWW[11];
   aryuWW[263]=aryuWW[37]*aryuWW[257];
   aryuWW[270]=aryuWW[263] - aryuWW[10] + aryuWW[236];
   aryuWW[270]=MMH*aryuWW[270];
   aryuWW[277]=aryuWW[133] + 1./4. + aryuWW[10];
   aryuWW[277]=aryuWW[38]*aryuWW[277];
   aryuWW[148]=1./4.*aryuWW[270] + 1./4.*aryuWW[148] + aryuWW[277];
   aryuWW[148]=MMH*aryuWW[148];
   aryuWW[185]=MMH*aryuWW[185];
   aryuWW[171]=1./3.*aryuWW[185] + aryuWW[189] + 1./6.*aryuWW[171];
   aryuWW[171]=MMt*aryuWW[171];
   aryuWW[185]= - 1./2.*aryuWW[18];
   aryuWW[189]=aryuWW[185] + 5./3.*aryuWW[19];
   aryuWW[189]=1./2.*aryuWW[189] + 1./3.*aryuWW[20];
   aryuWW[189]=aryuWW[38]*aryuWW[189];
   aryuWW[148]=aryuWW[171] + aryuWW[189] + 1./3.*aryuWW[148];
   aryuWW[148]=MMt*aryuWW[148];
   aryuWW[134]=aryuWW[37]*aryuWW[134];
   aryuWW[134]=1./2.*aryuWW[134] + aryuWW[309];
   aryuWW[134]=MMt*aryuWW[134];
   aryuWW[171]=aryuWW[38]*aryuWW[180];
   aryuWW[134]=aryuWW[134] + 1./2.*aryuWW[171] + aryuWW[313];
   aryuWW[134]=MMt*aryuWW[134];
   aryuWW[134]= - 1./2.*aryuWW[276] + aryuWW[134];
   aryuWW[134]=aryuWW[34]*MMt*aryuWW[134];
   aryuWW[171]=aryuWW[19] - aryuWW[20];
   aryuWW[171]=aryuWW[38]*aryuWW[171];
   aryuWW[189]= - aryuWW[10] + aryuWW[11];
   aryuWW[270]=MMH*aryuWW[38]*aryuWW[189];
   aryuWW[171]=aryuWW[171] + aryuWW[270];
   aryuWW[171]=MMH*aryuWW[171];
   aryuWW[134]=aryuWW[134] + 1./12.*aryuWW[171] + aryuWW[148];
   aryuWW[134]=aryuWW[34]*aryuWW[134];
   aryuWW[148]= - 1./3. + aryuWW[243];
   aryuWW[148]=aryuWW[19]*aryuWW[148];
   aryuWW[171]= - aryuWW[18]*aryuWW[10];
   aryuWW[148]=1./3.*aryuWW[171] + aryuWW[148];
   aryuWW[243]=1./18.*aryuWW[18] - aryuWW[19];
   aryuWW[243]=aryuWW[11]*aryuWW[243];
   aryuWW[241]=1./9. + aryuWW[241];
   aryuWW[241]=1./2.*aryuWW[241] + 1./3.*aryuWW[11];
   aryuWW[241]=aryuWW[20]*aryuWW[241];
   aryuWW[148]=aryuWW[241] + 1./6.*aryuWW[148] + aryuWW[243];
   aryuWW[241]= - aryuWW[11] - 1./4. + aryuWW[10];
   aryuWW[241]=aryuWW[11]*aryuWW[241];
   aryuWW[241]=aryuWW[181] + aryuWW[241];
   aryuWW[241]=MMH*aryuWW[241];
   aryuWW[148]=1./2.*aryuWW[148] + 1./9.*aryuWW[241];
   aryuWW[148]=MMH*aryuWW[148];
   aryuWW[241]= - 11./3.*aryuWW[18] + aryuWW[266];
   aryuWW[241]=aryuWW[19]*aryuWW[241];
   aryuWW[243]=11*aryuWW[18];
   aryuWW[270]=aryuWW[243] + 19*aryuWW[19];
   aryuWW[270]=1./4.*aryuWW[270] - aryuWW[20];
   aryuWW[270]=aryuWW[20]*aryuWW[270];
   aryuWW[241]=1./4.*aryuWW[241] + 1./3.*aryuWW[270];
   aryuWW[148]=1./6.*aryuWW[241] + aryuWW[148];
   aryuWW[148]=aryuWW[76]*MMH*aryuWW[148];
   aryuWW[134]=1./2.*aryuWW[148] + aryuWW[134];
   aryuWW[134]=aryuWW[165]*aryuWW[134];
   aryuWW[125]=aryuWW[134] + 1./2.*aryuWW[154] + aryuWW[125];
   aryuWW[125]=aryuWW[165]*aryuWW[125];
   aryuWW[134]=27./8.*aryuWW[44];
   aryuWW[148]=3./4.*aryuWW[48];
   aryuWW[154]=3./8.*aryuWW[46];
   aryuWW[241]=aryuWW[235] + aryuWW[154] + aryuWW[148] - 97./9. + 
   aryuWW[134];
   aryuWW[241]=aryuWW[18]*aryuWW[241];
   aryuWW[270]= - 27./2.*aryuWW[44];
   aryuWW[277]= - 3./2.*aryuWW[46];
   aryuWW[283]=aryuWW[277] + aryuWW[232] + 19./3. + aryuWW[270];
   aryuWW[283]= - 71./9.*aryuWW[11] + 1./2.*aryuWW[283] - 55./9.*
   aryuWW[10];
   aryuWW[283]=aryuWW[20]*aryuWW[283];
   aryuWW[134]=aryuWW[154] + aryuWW[148] - 1./3. + aryuWW[134];
   aryuWW[134]= - 11./48.*aryuWW[11] + 1./4.*aryuWW[134] + aryuWW[315];
   aryuWW[134]=aryuWW[11]*aryuWW[134];
   aryuWW[190]=aryuWW[190] + 17./8.*aryuWW[16] + 71./24.*aryuWW[100] - 
   3./16.*aryuWW[96] - 3./4.*aryuWW[91] - 5./8.*aryuWW[97] - 9./16.*
   aryuWW[95] + 87./32. - 1./3.*aryuWW[102];
   aryuWW[285]=3*aryuWW[84] + aryuWW[81];
   aryuWW[285]=MMH*aryuWW[285];
   aryuWW[134]=1./16.*aryuWW[285] + 1./2.*aryuWW[190] + aryuWW[134];
   aryuWW[134]=MMH*aryuWW[134];
   aryuWW[190]= - 2./3.*aryuWW[108];
   aryuWW[285]= - 43./3. + 11./4.*aryuWW[10];
   aryuWW[285]=aryuWW[19]*aryuWW[285];
   aryuWW[286]= - aryuWW[18] - 101*aryuWW[19];
   aryuWW[286]=aryuWW[11]*aryuWW[286];
   aryuWW[134]=aryuWW[134] + 1./8.*aryuWW[283] + 1./72.*aryuWW[286] + 1.
   /6.*aryuWW[285] + 1./4.*aryuWW[241] - 7./48.*aryuWW[23] + 15./16.*
   aryuWW[21] - 53./48.*aryuWW[22] - 55./96.*aryuWW[78] + 3./32.*
   aryuWW[107] - 31./24.*aryuWW[79] + 71./24.*aryuWW[106] + 11./32.*
   aryuWW[105] + aryuWW[190] - 15./32.*aryuWW[77];
   aryuWW[134]=MMH*aryuWW[134];
   aryuWW[241]=aryuWW[243] + aryuWW[19];
   aryuWW[241]=aryuWW[19]*aryuWW[241];
   aryuWW[241]=5*aryuWW[162] + 1./2.*aryuWW[241];
   aryuWW[243]= - 21*aryuWW[18] - 295./9.*aryuWW[19];
   aryuWW[243]=1./8.*aryuWW[243] + 89./9.*aryuWW[20];
   aryuWW[243]=aryuWW[20]*aryuWW[243];
   aryuWW[241]=1./12.*aryuWW[241] + aryuWW[243];
   aryuWW[134]=1./2.*aryuWW[241] + aryuWW[134];
   aryuWW[134]=aryuWW[76]*MMH*aryuWW[134];
   aryuWW[125]=aryuWW[125] + aryuWW[134] + aryuWW[138];
   aryuWW[125]=aryuWW[165]*aryuWW[125];
   aryuWW[138]=aryuWW[253] + aryuWW[252] - 5./3. + aryuWW[251];
   aryuWW[138]=aryuWW[18]*aryuWW[138];
   aryuWW[241]=aryuWW[196] + aryuWW[151] - aryuWW[48] + 19 + 
   aryuWW[221];
   aryuWW[241]=1./8.*aryuWW[241] - 11./9.*aryuWW[11];
   aryuWW[241]=aryuWW[20]*aryuWW[241];
   aryuWW[243]= - 1./6.*aryuWW[11] + aryuWW[231] + aryuWW[317] - 1 + 
   aryuWW[316];
   aryuWW[243]=aryuWW[11]*aryuWW[243];
   aryuWW[252]=31./32. - aryuWW[102];
   aryuWW[243]=1./8.*aryuWW[319] + 1./4.*aryuWW[243] + 17./24.*
   aryuWW[16] + 5./8.*aryuWW[100] - 1./16.*aryuWW[96] - 1./4.*
   aryuWW[91] - 5./24.*aryuWW[97] + 1./3.*aryuWW[252] - 3./16.*
   aryuWW[95];
   aryuWW[243]=MMH*aryuWW[243];
   aryuWW[235]= - 9 + aryuWW[235];
   aryuWW[235]=aryuWW[19]*aryuWW[235];
   aryuWW[252]=aryuWW[11]*aryuWW[19];
   aryuWW[253]= - 1./2.*aryuWW[79];
   aryuWW[138]=1./2.*aryuWW[243] + 1./2.*aryuWW[241] + 5./18.*
   aryuWW[252] + 1./8.*aryuWW[235] + 1./4.*aryuWW[138] + 43./48.*
   aryuWW[23] + 5./16.*aryuWW[21] + 11./48.*aryuWW[22] - 5./96.*
   aryuWW[78] + 1./32.*aryuWW[107] + aryuWW[253] + 5./8.*aryuWW[106] + 
   11./96.*aryuWW[105] + aryuWW[190] - 5./32.*aryuWW[77];
   aryuWW[138]=MMH*aryuWW[138];
   aryuWW[235]=aryuWW[19]*aryuWW[18];
   aryuWW[241]=aryuWW[162] + 5./3.*aryuWW[235];
   aryuWW[243]= - 25*aryuWW[18];
   aryuWW[283]=aryuWW[243] - 41./2.*aryuWW[19];
   aryuWW[283]=1./2.*aryuWW[283] + 53*aryuWW[20];
   aryuWW[283]=aryuWW[20]*aryuWW[283];
   aryuWW[241]=1./4.*aryuWW[241] + 1./3.*aryuWW[283];
   aryuWW[138]=1./6.*aryuWW[241] + aryuWW[138];
   aryuWW[138]=aryuWW[76]*MMH*aryuWW[138];
   aryuWW[241]=1./2.*aryuWW[22];
   aryuWW[283]= - aryuWW[18] + aryuWW[241] + aryuWW[106] + aryuWW[253];
   aryuWW[285]=1 + aryuWW[100];
   aryuWW[285]=MMH*aryuWW[285];
   aryuWW[286]= - 1 - aryuWW[11];
   aryuWW[286]=aryuWW[20]*aryuWW[286];
   aryuWW[252]=1./4.*aryuWW[285] + 1./4.*aryuWW[286] + 1./4.*
   aryuWW[252] + 1./2.*aryuWW[283] - aryuWW[19];
   aryuWW[252]=MMH*aryuWW[252];
   aryuWW[283]=1./2.*aryuWW[152] + aryuWW[20];
   aryuWW[285]=aryuWW[20]*aryuWW[283];
   aryuWW[252]=1./2.*aryuWW[285] + aryuWW[252];
   aryuWW[252]=aryuWW[255]*aryuWW[76]*MMH*aryuWW[252];
   aryuWW[138]=aryuWW[138] + 1./12.*aryuWW[252];
   aryuWW[138]=aryuWW[255]*aryuWW[138];
   aryuWW[252]=aryuWW[315] + aryuWW[231] + aryuWW[317] - 47./9. + 
   aryuWW[316];
   aryuWW[252]=aryuWW[18]*aryuWW[252];
   aryuWW[285]=aryuWW[157] + aryuWW[48] - 7./9. + aryuWW[149];
   aryuWW[285]= - 7./18.*aryuWW[11] + 1./2.*aryuWW[285] + aryuWW[315];
   aryuWW[285]=aryuWW[11]*aryuWW[285];
   aryuWW[288]=61./16. - aryuWW[102];
   aryuWW[285]=aryuWW[320] + 1./2.*aryuWW[285] + 17./12.*aryuWW[16] + 
   37./24.*aryuWW[100] - 1./8.*aryuWW[96] - 1./2.*aryuWW[91] - 5./12.*
   aryuWW[97] + 1./3.*aryuWW[288] - 3./8.*aryuWW[95];
   aryuWW[285]=MMH*aryuWW[285];
   aryuWW[288]= - 9*aryuWW[44];
   aryuWW[290]=131./9. + aryuWW[288];
   aryuWW[290]= - 97./18.*aryuWW[11] - 49./18.*aryuWW[10] + aryuWW[151]
    + 1./2.*aryuWW[290] - aryuWW[48];
   aryuWW[290]=aryuWW[20]*aryuWW[290];
   aryuWW[293]= - 5 + 11./16.*aryuWW[10];
   aryuWW[293]=aryuWW[19]*aryuWW[293];
   aryuWW[294]= - 1./3.*aryuWW[18] - 13./2.*aryuWW[19];
   aryuWW[294]=aryuWW[11]*aryuWW[294];
   aryuWW[190]=1./2.*aryuWW[285] + 1./8.*aryuWW[290] + 1./24.*
   aryuWW[294] + 1./3.*aryuWW[293] + 1./4.*aryuWW[252] + 5./8.*
   aryuWW[23] + 5./8.*aryuWW[21] - 19./48.*aryuWW[22] - 3./16.*
   aryuWW[78] + 1./16.*aryuWW[107] - 43./48.*aryuWW[79] + 37./24.*
   aryuWW[106] + 11./48.*aryuWW[105] + aryuWW[190] - 5./16.*aryuWW[77];
   aryuWW[190]=MMH*aryuWW[190];
   aryuWW[252]=17*aryuWW[18] + aryuWW[266];
   aryuWW[252]=aryuWW[19]*aryuWW[252];
   aryuWW[243]=91*aryuWW[20] + aryuWW[243] - 29*aryuWW[19];
   aryuWW[243]=aryuWW[20]*aryuWW[243];
   aryuWW[243]=1./3.*aryuWW[243] + aryuWW[162] + 1./6.*aryuWW[252];
   aryuWW[190]=1./8.*aryuWW[243] + aryuWW[190];
   aryuWW[190]=aryuWW[76]*MMH*aryuWW[190];
   aryuWW[138]=aryuWW[190] + aryuWW[138];
   aryuWW[138]=aryuWW[255]*aryuWW[138];
   aryuWW[134]=aryuWW[134] + aryuWW[138];
   aryuWW[134]=aryuWW[255]*aryuWW[134];
   aryuWW[138]=aryuWW[244] + aryuWW[105] + aryuWW[253];
   aryuWW[190]=1./3.*aryuWW[138];
   aryuWW[243]= - 1./2. + aryuWW[199];
   aryuWW[133]=1./3.*aryuWW[243] + aryuWW[133];
   aryuWW[133]=aryuWW[20]*aryuWW[133];
   aryuWW[243]= - 1 - 1./9.*aryuWW[10];
   aryuWW[243]=aryuWW[18]*aryuWW[243];
   aryuWW[252]=7./2.*aryuWW[18] + aryuWW[19];
   aryuWW[252]=aryuWW[11]*aryuWW[252];
   aryuWW[133]=aryuWW[133] + 1./9.*aryuWW[252] + aryuWW[190] + 
   aryuWW[243];
   aryuWW[219]=aryuWW[196] + aryuWW[219];
   aryuWW[219]=aryuWW[11]*aryuWW[219];
   aryuWW[243]= - 1./4.*aryuWW[16] + 1./4.*aryuWW[91] + 1./2. + 
   aryuWW[101];
   aryuWW[219]=aryuWW[243] + 1./2.*aryuWW[219];
   aryuWW[219]=MMH*aryuWW[219];
   aryuWW[133]=1./2.*aryuWW[133] + 1./3.*aryuWW[219];
   aryuWW[133]=MMH*aryuWW[133];
   aryuWW[219]=1./2.*aryuWW[162] + aryuWW[235];
   aryuWW[252]= - 7./2.*aryuWW[18] - aryuWW[19];
   aryuWW[252]=1./6.*aryuWW[252] - aryuWW[20];
   aryuWW[252]=aryuWW[20]*aryuWW[252];
   aryuWW[219]=1./6.*aryuWW[219] + aryuWW[252];
   aryuWW[133]=1./3.*aryuWW[219] + aryuWW[133];
   aryuWW[219]=pow(MMH,2);
   aryuWW[133]=aryuWW[76]*aryuWW[219]*aryuWW[133];
   aryuWW[136]=aryuWW[136] + aryuWW[226] - aryuWW[73] + aryuWW[225];
   aryuWW[252]= - aryuWW[70] - 1 - aryuWW[60];
   aryuWW[253]=1./2.*aryuWW[16];
   aryuWW[252]=aryuWW[253] + 1./2.*aryuWW[252] - aryuWW[68];
   aryuWW[285]=aryuWW[199] - aryuWW[11];
   aryuWW[285]=aryuWW[37]*aryuWW[285];
   aryuWW[285]=aryuWW[252] + aryuWW[285];
   aryuWW[285]=MMH*aryuWW[285];
   aryuWW[290]= - aryuWW[18] + aryuWW[265];
   aryuWW[293]=aryuWW[37]*aryuWW[290];
   aryuWW[294]=aryuWW[20]*aryuWW[37];
   aryuWW[295]=1./2.*aryuWW[294];
   aryuWW[285]=1./2.*aryuWW[285] + 1./2.*aryuWW[38] + aryuWW[295] + 1./
   2.*aryuWW[136] + 1./3.*aryuWW[293];
   aryuWW[285]=MMH*aryuWW[285];
   aryuWW[293]= - 1 - aryuWW[70];
   aryuWW[293]=MMt*MMH*aryuWW[293];
   aryuWW[285]=aryuWW[285] + 1./4.*aryuWW[293];
   aryuWW[285]=MMt*aryuWW[285];
   aryuWW[298]= - 1./2.*aryuWW[51];
   aryuWW[300]= - 3./2.*aryuWW[11] + 1./2. + aryuWW[199];
   aryuWW[300]=aryuWW[38]*aryuWW[300];
   aryuWW[300]=aryuWW[298] + aryuWW[300];
   aryuWW[300]=MMH*aryuWW[300];
   aryuWW[290]=1./3.*aryuWW[290] + 3./4.*aryuWW[20];
   aryuWW[290]=aryuWW[38]*aryuWW[290];
   aryuWW[290]=aryuWW[290] + 1./2.*aryuWW[300];
   aryuWW[290]=MMH*aryuWW[290];
   aryuWW[285]=aryuWW[290] + aryuWW[285];
   aryuWW[285]=MMt*aryuWW[285];
   aryuWW[290]=MMt*pow(aryuWW[37],2);
   aryuWW[237]=aryuWW[237] + 1./2.*aryuWW[290];
   aryuWW[237]=MMt*aryuWW[237];
   aryuWW[237]=1./2.*aryuWW[276] + aryuWW[237];
   aryuWW[237]=aryuWW[237]*pow(MMt,2);
   aryuWW[276]=aryuWW[34]*aryuWW[237];
   aryuWW[290]=3*aryuWW[276];
   aryuWW[285]=aryuWW[285] + aryuWW[290];
   aryuWW[285]=aryuWW[34]*aryuWW[285];
   aryuWW[200]=aryuWW[11]*aryuWW[200];
   aryuWW[171]=aryuWW[171] + aryuWW[200];
   aryuWW[200]=aryuWW[20]*aryuWW[257];
   aryuWW[189]=MMH*aryuWW[11]*aryuWW[189];
   aryuWW[171]=1./2.*aryuWW[189] + 1./2.*aryuWW[171] + aryuWW[200];
   aryuWW[171]=MMH*aryuWW[171];
   aryuWW[152]=aryuWW[152] + aryuWW[20];
   aryuWW[152]=aryuWW[20]*aryuWW[152];
   aryuWW[152]=aryuWW[235] + aryuWW[152];
   aryuWW[152]=1./2.*aryuWW[152] + aryuWW[171];
   aryuWW[152]=aryuWW[76]*aryuWW[219]*aryuWW[152];
   aryuWW[171]=MMH*aryuWW[263];
   aryuWW[171]=aryuWW[171] + aryuWW[204] + aryuWW[294];
   aryuWW[171]=MMt*MMH*aryuWW[171];
   aryuWW[189]=aryuWW[38]*aryuWW[283];
   aryuWW[200]=MMH*aryuWW[38]*aryuWW[257];
   aryuWW[189]=aryuWW[189] + aryuWW[200];
   aryuWW[189]=MMH*aryuWW[189];
   aryuWW[171]=aryuWW[189] + aryuWW[171];
   aryuWW[171]=MMt*aryuWW[171];
   aryuWW[171]=1./6.*aryuWW[171] + aryuWW[276];
   aryuWW[171]=aryuWW[34]*aryuWW[171];
   aryuWW[152]=1./36.*aryuWW[152] + aryuWW[171];
   aryuWW[152]=aryuWW[165]*aryuWW[152];
   aryuWW[133]=aryuWW[152] + 1./2.*aryuWW[133] + aryuWW[285];
   aryuWW[133]=aryuWW[165]*aryuWW[133];
   aryuWW[152]= - 1 - 1./24.*aryuWW[10];
   aryuWW[152]=aryuWW[18]*aryuWW[152];
   aryuWW[171]=aryuWW[147] + aryuWW[19];
   aryuWW[171]=aryuWW[11]*aryuWW[171];
   aryuWW[189]= - 1 + aryuWW[181];
   aryuWW[189]=1./2.*aryuWW[189] - aryuWW[11];
   aryuWW[189]=aryuWW[20]*aryuWW[189];
   aryuWW[152]=1./3.*aryuWW[189] + 1./24.*aryuWW[171] + aryuWW[190] + 
   aryuWW[152];
   aryuWW[171]=aryuWW[11]*aryuWW[213];
   aryuWW[171]=aryuWW[243] + 1./8.*aryuWW[171];
   aryuWW[171]=MMH*aryuWW[171];
   aryuWW[152]=1./2.*aryuWW[152] + 1./3.*aryuWW[171];
   aryuWW[152]=MMH*aryuWW[152];
   aryuWW[171]=aryuWW[162] + aryuWW[235];
   aryuWW[189]= - 7*aryuWW[18];
   aryuWW[200]=aryuWW[189] - aryuWW[19];
   aryuWW[204]= - 7*aryuWW[20];
   aryuWW[200]=1./6.*aryuWW[200] + aryuWW[204];
   aryuWW[200]=aryuWW[20]*aryuWW[200];
   aryuWW[171]=1./6.*aryuWW[171] + aryuWW[200];
   aryuWW[152]=1./8.*aryuWW[171] + aryuWW[152];
   aryuWW[152]=aryuWW[76]*aryuWW[219]*aryuWW[152];
   aryuWW[171]=aryuWW[252] + 1./2.*aryuWW[198];
   aryuWW[171]=MMH*aryuWW[171];
   aryuWW[198]= - 3*aryuWW[18];
   aryuWW[200]=aryuWW[198] - aryuWW[19];
   aryuWW[213]=aryuWW[37]*aryuWW[200];
   aryuWW[171]=aryuWW[171] + aryuWW[38] + aryuWW[295] + aryuWW[136] + 1.
   /8.*aryuWW[213];
   aryuWW[171]=MMH*aryuWW[171];
   aryuWW[213]=1./2.*aryuWW[293];
   aryuWW[171]=aryuWW[171] + aryuWW[213];
   aryuWW[171]=MMt*aryuWW[171];
   aryuWW[181]=1 + aryuWW[181];
   aryuWW[181]=1./2.*aryuWW[181] - aryuWW[11];
   aryuWW[181]=aryuWW[38]*aryuWW[181];
   aryuWW[181]=aryuWW[298] + aryuWW[181];
   aryuWW[181]=MMH*aryuWW[181];
   aryuWW[200]=1./8.*aryuWW[200] + aryuWW[20];
   aryuWW[200]=aryuWW[38]*aryuWW[200];
   aryuWW[181]=aryuWW[200] + aryuWW[181];
   aryuWW[181]=MMH*aryuWW[181];
   aryuWW[171]=aryuWW[181] + aryuWW[171];
   aryuWW[171]=MMt*aryuWW[171];
   aryuWW[181]=aryuWW[171] + aryuWW[290];
   aryuWW[181]=aryuWW[34]*aryuWW[181];
   aryuWW[133]=1./2.*aryuWW[133] + aryuWW[152] + aryuWW[181];
   aryuWW[133]=aryuWW[165]*aryuWW[133];
   aryuWW[147]=aryuWW[147] + aryuWW[245];
   aryuWW[147]=aryuWW[11]*aryuWW[147];
   aryuWW[181]= - 3 - 1./18.*aryuWW[10];
   aryuWW[181]=aryuWW[18]*aryuWW[181];
   aryuWW[200]= - 5./3.*aryuWW[11] - 1 + aryuWW[315];
   aryuWW[200]=aryuWW[20]*aryuWW[200];
   aryuWW[138]=1./2.*aryuWW[200] + 1./9.*aryuWW[147] + aryuWW[138] + 
   aryuWW[181];
   aryuWW[147]=aryuWW[196] + aryuWW[11];
   aryuWW[147]=aryuWW[11]*aryuWW[147];
   aryuWW[147]=aryuWW[243] + 1./12.*aryuWW[147];
   aryuWW[147]=MMH*aryuWW[147];
   aryuWW[138]=1./2.*aryuWW[138] + aryuWW[147];
   aryuWW[138]=MMH*aryuWW[138];
   aryuWW[147]=aryuWW[162] + 1./2.*aryuWW[235];
   aryuWW[181]=aryuWW[189] + aryuWW[265];
   aryuWW[181]=1./3.*aryuWW[181] - 19*aryuWW[20];
   aryuWW[181]=aryuWW[20]*aryuWW[181];
   aryuWW[147]=1./3.*aryuWW[147] + aryuWW[181];
   aryuWW[138]=1./6.*aryuWW[147] + aryuWW[138];
   aryuWW[138]=aryuWW[76]*aryuWW[219]*aryuWW[138];
   aryuWW[147]=pow(aryuWW[11],2);
   aryuWW[147]=aryuWW[243] + 1./24.*aryuWW[147];
   aryuWW[147]=MMH*aryuWW[147];
   aryuWW[181]=aryuWW[190] - aryuWW[18];
   aryuWW[189]= - 1./3.*aryuWW[11];
   aryuWW[190]= - 1./4. + aryuWW[189];
   aryuWW[190]=aryuWW[20]*aryuWW[190];
   aryuWW[200]=aryuWW[11]*aryuWW[18];
   aryuWW[147]=1./3.*aryuWW[147] + 1./3.*aryuWW[190] + 1./2.*
   aryuWW[181] + 1./9.*aryuWW[200];
   aryuWW[147]=MMH*aryuWW[147];
   aryuWW[190]= - aryuWW[18] - 101./8.*aryuWW[20];
   aryuWW[190]=aryuWW[20]*aryuWW[190];
   aryuWW[190]=1./8.*aryuWW[162] + aryuWW[190];
   aryuWW[147]=1./9.*aryuWW[190] + aryuWW[147];
   aryuWW[147]=aryuWW[76]*aryuWW[219]*aryuWW[147];
   aryuWW[181]=1./6.*aryuWW[286] + aryuWW[181] + 1./6.*aryuWW[200];
   aryuWW[190]=MMH*aryuWW[243];
   aryuWW[181]=1./2.*aryuWW[181] + 1./3.*aryuWW[190];
   aryuWW[181]=MMH*aryuWW[181];
   aryuWW[190]= - aryuWW[18] + aryuWW[272];
   aryuWW[190]=aryuWW[20]*aryuWW[190];
   aryuWW[181]=1./12.*aryuWW[190] + aryuWW[181];
   aryuWW[181]=aryuWW[255]*aryuWW[76]*aryuWW[219]*aryuWW[181];
   aryuWW[147]=aryuWW[147] + 1./2.*aryuWW[181];
   aryuWW[147]=aryuWW[255]*aryuWW[147];
   aryuWW[138]=1./2.*aryuWW[138] + aryuWW[147];
   aryuWW[138]=aryuWW[255]*aryuWW[138];
   aryuWW[138]=aryuWW[152] + 1./2.*aryuWW[138];
   aryuWW[138]=aryuWW[255]*aryuWW[138];
   aryuWW[147]=3*aryuWW[252] + aryuWW[312];
   aryuWW[147]=MMH*aryuWW[147];
   aryuWW[147]=aryuWW[147] + 3*aryuWW[38] + aryuWW[294] + 3*aryuWW[136]
    + aryuWW[307];
   aryuWW[147]=MMH*aryuWW[147];
   aryuWW[147]=aryuWW[147] + 3./2.*aryuWW[293];
   aryuWW[147]=MMt*aryuWW[147];
   aryuWW[152]=aryuWW[311] + 3 + aryuWW[199];
   aryuWW[152]=aryuWW[38]*aryuWW[152];
   aryuWW[152]=aryuWW[224] + aryuWW[152];
   aryuWW[152]=MMH*aryuWW[152];
   aryuWW[181]=1./3.*aryuWW[292] + 5*aryuWW[20];
   aryuWW[181]=aryuWW[38]*aryuWW[181];
   aryuWW[152]=aryuWW[181] + aryuWW[152];
   aryuWW[152]=MMH*aryuWW[152];
   aryuWW[147]=1./2.*aryuWW[152] + aryuWW[147];
   aryuWW[147]=MMt*aryuWW[147];
   aryuWW[152]=aryuWW[252] + aryuWW[297];
   aryuWW[152]=MMH*aryuWW[152];
   aryuWW[152]=aryuWW[152] + aryuWW[38] + 1./6.*aryuWW[294] + 
   aryuWW[136] + aryuWW[282];
   aryuWW[152]=MMH*aryuWW[152];
   aryuWW[152]=aryuWW[152] + aryuWW[213];
   aryuWW[152]=MMt*aryuWW[152];
   aryuWW[181]=aryuWW[188] + aryuWW[20];
   aryuWW[181]=aryuWW[38]*aryuWW[181];
   aryuWW[189]=1./4. + aryuWW[189];
   aryuWW[189]=aryuWW[38]*aryuWW[189];
   aryuWW[189]= - 1./4.*aryuWW[51] + aryuWW[189];
   aryuWW[189]=MMH*aryuWW[189];
   aryuWW[181]=1./3.*aryuWW[181] + aryuWW[189];
   aryuWW[181]=MMH*aryuWW[181];
   aryuWW[152]=aryuWW[181] + 1./2.*aryuWW[152];
   aryuWW[152]=MMt*aryuWW[152];
   aryuWW[181]=MMH*aryuWW[252];
   aryuWW[136]=aryuWW[181] + aryuWW[136] + aryuWW[38];
   aryuWW[136]=MMH*aryuWW[136];
   aryuWW[136]=aryuWW[136] + aryuWW[213];
   aryuWW[136]=MMt*aryuWW[136];
   aryuWW[181]=1 - aryuWW[11];
   aryuWW[181]=aryuWW[38]*aryuWW[181];
   aryuWW[181]= - aryuWW[51] + aryuWW[181];
   aryuWW[181]=MMH*aryuWW[181];
   aryuWW[181]=aryuWW[308] + aryuWW[181];
   aryuWW[181]=MMH*aryuWW[181];
   aryuWW[136]=1./2.*aryuWW[181] + aryuWW[136];
   aryuWW[136]=aryuWW[255]*MMt*aryuWW[136];
   aryuWW[136]=aryuWW[152] + 1./4.*aryuWW[136];
   aryuWW[136]=aryuWW[255]*aryuWW[136];
   aryuWW[136]=1./4.*aryuWW[147] + aryuWW[136];
   aryuWW[136]=aryuWW[255]*aryuWW[136];
   aryuWW[136]=aryuWW[171] + aryuWW[136];
   aryuWW[136]=aryuWW[255]*aryuWW[136];
   aryuWW[147]=3*aryuWW[237];
   aryuWW[152]=aryuWW[255]*aryuWW[237];
   aryuWW[152]=aryuWW[147] + aryuWW[152];
   aryuWW[152]=aryuWW[255]*aryuWW[152];
   aryuWW[147]=aryuWW[147] + 1./2.*aryuWW[152];
   aryuWW[147]=aryuWW[34]*aryuWW[255]*aryuWW[147];
   aryuWW[136]=aryuWW[136] + aryuWW[147];
   aryuWW[136]=aryuWW[34]*aryuWW[136];
   aryuWW[133]=aryuWW[133] + aryuWW[138] + aryuWW[136];
   aryuWW[133]=aryuWW[4]*aryuWW[133];
   aryuWW[125]=aryuWW[133] + aryuWW[125] + aryuWW[134] + aryuWW[159];
   aryuWW[125]=aryuWW[4]*aryuWW[125];
   aryuWW[114]=aryuWW[125] + aryuWW[114] + aryuWW[120] + aryuWW[130];
   aryuWW[114]=aryuWW[4]*aryuWW[114];
   aryuWW[120]=1./2.*aryuWW[13];
   aryuWW[125]= - 22*aryuWW[47];
   aryuWW[130]= - 6*aryuWW[46];
   aryuWW[133]= - 19./2.*aryuWW[12];
   aryuWW[134]= - 1 + 3*aryuWW[6];
   aryuWW[134]=2*aryuWW[1]*aryuWW[134];
   aryuWW[136]=22./3.*aryuWW[207];
   aryuWW[138]= - 3*aryuWW[49];
   aryuWW[147]= - 6*aryuWW[48];
   aryuWW[152]=6*aryuWW[10];
   aryuWW[159]= - 3*aryuWW[11];
   aryuWW[171]=aryuWW[159] + aryuWW[136] + aryuWW[134] + aryuWW[133] + 
   aryuWW[152] + aryuWW[120] + aryuWW[130] + aryuWW[147] + aryuWW[138]
    - 37./3. + aryuWW[125];
   aryuWW[171]=aryuWW[38]*aryuWW[171];
   aryuWW[181]= - 9*aryuWW[39];
   aryuWW[189]= - 139./6.*aryuWW[36] - 41./3. + aryuWW[181];
   aryuWW[189]=aryuWW[19]*aryuWW[189];
   aryuWW[190]=1./4.*aryuWW[21];
   aryuWW[142]=aryuWW[18] + aryuWW[142];
   aryuWW[142]=1./4.*aryuWW[37]*aryuWW[142];
   aryuWW[189]=aryuWW[142] + 1./2.*aryuWW[189] + aryuWW[188] + 
   aryuWW[222] + aryuWW[190] - 11*aryuWW[22] - 25*aryuWW[51] + 11*
   aryuWW[42] + 3*aryuWW[43];
   aryuWW[200]= - 1 + aryuWW[131];
   aryuWW[200]=1./2.*aryuWW[200] - aryuWW[37];
   aryuWW[200]=aryuWW[20]*aryuWW[200];
   aryuWW[171]=aryuWW[171] + aryuWW[189] + 3*aryuWW[200];
   aryuWW[171]=aryuWW[3]*aryuWW[171];
   aryuWW[200]=1./3.*aryuWW[13];
   aryuWW[213]= - 21*aryuWW[47];
   aryuWW[219]=3*aryuWW[49];
   aryuWW[235]=aryuWW[200] + aryuWW[219] - 14999./144. + aryuWW[213];
   aryuWW[235]=1./4.*aryuWW[235] + aryuWW[10];
   aryuWW[237]=14./3.*aryuWW[12];
   aryuWW[243]=17./18.*aryuWW[37];
   aryuWW[235]=aryuWW[243] + aryuWW[218] + aryuWW[217] + aryuWW[216] + 
   1./2.*aryuWW[235] + aryuWW[237];
   aryuWW[235]=aryuWW[37]*aryuWW[235];
   aryuWW[252]=aryuWW[19]*aryuWW[39];
   aryuWW[257]=aryuWW[20]*aryuWW[39];
   aryuWW[263]=aryuWW[252] + aryuWW[257];
   aryuWW[263]=aryuWW[3]*aryuWW[263];
   aryuWW[272]= - 1 - aryuWW[39];
   aryuWW[272]=aryuWW[11]*aryuWW[272];
   aryuWW[169]=18*aryuWW[263] - 5*aryuWW[37] + 3*aryuWW[272] + 
   aryuWW[160] + aryuWW[124] + 142./3. + aryuWW[169];
   aryuWW[169]=aryuWW[3]*aryuWW[169];
   aryuWW[263]= - 3*aryuWW[68];
   aryuWW[191]=aryuWW[160] - 5./2.*aryuWW[16] + aryuWW[263] - 1./2.*
   aryuWW[70] + aryuWW[287] - 43./8. + aryuWW[191];
   aryuWW[191]=aryuWW[112]*aryuWW[191];
   aryuWW[276]=4*aryuWW[57] - 11*aryuWW[56];
   aryuWW[282]= - 3*aryuWW[11]*aryuWW[112]*aryuWW[39];
   aryuWW[283]=aryuWW[11]*aryuWW[112];
   aryuWW[283]= - aryuWW[112] + aryuWW[283];
   aryuWW[283]=1./2.*aryuWW[37]*aryuWW[283];
   aryuWW[169]=aryuWW[169] + aryuWW[283] + aryuWW[282] + aryuWW[191] - 
   aryuWW[55] - 11./12.*aryuWW[59] + 2./3.*aryuWW[276] - 3./2.*
   aryuWW[54];
   aryuWW[169]=MMt*aryuWW[169];
   aryuWW[276]= - 139523./128. + 580*aryuWW[64];
   aryuWW[285]=13./2.*aryuWW[53];
   aryuWW[276]=1./3.*aryuWW[276] + aryuWW[285];
   aryuWW[286]=101./8.*aryuWW[61];
   aryuWW[290]=49./8.*aryuWW[67];
   aryuWW[292]=349./128.*aryuWW[66];
   aryuWW[276]= - 8851./144.*aryuWW[63] + aryuWW[292] + aryuWW[290] + 1.
   /3.*aryuWW[276] + aryuWW[286];
   aryuWW[293]=aryuWW[72] - 1./2.*aryuWW[41];
   aryuWW[294]= - 7./2.*aryuWW[23];
   aryuWW[224]=aryuWW[294] + 7./4.*aryuWW[21] + aryuWW[224] + 
   aryuWW[225] + 3*aryuWW[293] - aryuWW[73];
   aryuWW[224]=aryuWW[112]*aryuWW[224];
   aryuWW[225]= - 1./2. - aryuWW[39];
   aryuWW[225]=aryuWW[112]*aryuWW[225];
   aryuWW[141]=3*aryuWW[225] + aryuWW[141];
   aryuWW[141]=aryuWW[20]*aryuWW[141];
   aryuWW[225]= - 1./12.*aryuWW[13];
   aryuWW[293]= - 10741./36. - 151*aryuWW[36];
   aryuWW[293]=aryuWW[12]*aryuWW[293];
   aryuWW[295]=19./12.*aryuWW[33];
   aryuWW[297]=13./12.*aryuWW[30];
   aryuWW[300]=169./192.*aryuWW[69];
   aryuWW[301]=1651./384.*aryuWW[71];
   aryuWW[302]= - 1037./384.*aryuWW[50];
   aryuWW[303]=11./2.*aryuWW[47];
   aryuWW[307]= - 1./4.*aryuWW[65];
   aryuWW[308]= - 13./24.*aryuWW[70];
   aryuWW[160]=5./2. + aryuWW[160];
   aryuWW[160]=1./2.*aryuWW[18]*aryuWW[112]*aryuWW[160];
   aryuWW[163]=5*aryuWW[163];
   aryuWW[171]=aryuWW[169] + aryuWW[284] + aryuWW[171] + aryuWW[163] + 
   aryuWW[141] + aryuWW[235] + aryuWW[304] + aryuWW[160] + aryuWW[275]
    + aryuWW[259] + aryuWW[224] + 1./12.*aryuWW[293] - aryuWW[10] + 139.
   /36.*aryuWW[36] + aryuWW[225] + aryuWW[192] + aryuWW[291] + 6565./
   432.*aryuWW[17] + aryuWW[233] + aryuWW[303] + aryuWW[262] + 
   aryuWW[308] + aryuWW[302] + aryuWW[287] + aryuWW[301] + aryuWW[300]
    + aryuWW[297] + aryuWW[295] - aryuWW[62] + 1./3.*aryuWW[276] + 
   aryuWW[307];
   aryuWW[171]=MMt*aryuWW[171];
   aryuWW[235]= - 17*aryuWW[36];
   aryuWW[276]= - 5*aryuWW[6];
   aryuWW[293]=aryuWW[276] + 161./6. + aryuWW[235];
   aryuWW[293]=1./9.*aryuWW[18]*aryuWW[293];
   aryuWW[304]=29*aryuWW[74];
   aryuWW[309]= - 43./2.*aryuWW[8] + aryuWW[73] - 1345./12.*aryuWW[72]
    + aryuWW[304] - 1915./96.*aryuWW[42];
   aryuWW[311]=37*aryuWW[6] - 7505./64. + aryuWW[235];
   aryuWW[311]=aryuWW[19]*aryuWW[311];
   aryuWW[312]=aryuWW[158] - 281./288.*aryuWW[19];
   aryuWW[312]=aryuWW[37]*aryuWW[312];
   aryuWW[313]= - 167*aryuWW[6] - 13525./24. + 19*aryuWW[36];
   aryuWW[313]=1./9.*aryuWW[313] + 133./2.*aryuWW[37];
   aryuWW[313]=aryuWW[20]*aryuWW[313];
   aryuWW[315]= - 59./2.*aryuWW[43];
   aryuWW[309]=1./3.*aryuWW[313] + 1./2.*aryuWW[312] + 1./9.*
   aryuWW[311] + aryuWW[293] + 2617./36.*aryuWW[23] + aryuWW[139] + 789.
   /64.*aryuWW[22] + 36127./288.*aryuWW[51] + 1./3.*aryuWW[309] + 
   aryuWW[315];
   aryuWW[311]=5./3.*aryuWW[35];
   aryuWW[312]= - 69./4.*aryuWW[110];
   aryuWW[313]=407./16.*aryuWW[111];
   aryuWW[319]=aryuWW[313] + aryuWW[311] + aryuWW[312];
   aryuWW[319]=aryuWW[38]*aryuWW[319];
   aryuWW[320]=200183./432. + aryuWW[47];
   aryuWW[320]= - 1./3.*aryuWW[13] + 1./2.*aryuWW[320] + aryuWW[219];
   aryuWW[322]=1./2.*aryuWW[110];
   aryuWW[323]=1./3.*aryuWW[35] + aryuWW[322];
   aryuWW[324]=aryuWW[323] - 793./48.*aryuWW[111];
   aryuWW[324]=aryuWW[19]*aryuWW[324];
   aryuWW[326]= - aryuWW[35] - 51./2.*aryuWW[110];
   aryuWW[327]=aryuWW[326] + 107./6.*aryuWW[111];
   aryuWW[327]=aryuWW[20]*aryuWW[327];
   aryuWW[328]=6889./1152.*aryuWW[37];
   aryuWW[275]=1./2.*aryuWW[319] + 1./2.*aryuWW[327] + aryuWW[328] + 1./
   2.*aryuWW[324] + aryuWW[275] + aryuWW[259] - 8959./432.*aryuWW[12]
    + 1./4.*aryuWW[320] - aryuWW[10];
   aryuWW[275]=aryuWW[38]*aryuWW[275];
   aryuWW[319]= - 1./2.*aryuWW[42];
   aryuWW[241]=aryuWW[241] + aryuWW[319] + aryuWW[51];
   aryuWW[320]=7./3. + aryuWW[181];
   aryuWW[267]=1./2.*aryuWW[320] + aryuWW[267];
   aryuWW[267]=aryuWW[19]*aryuWW[267];
   aryuWW[320]=2 + aryuWW[206];
   aryuWW[324]= - 3*aryuWW[46];
   aryuWW[320]= - 1./2.*aryuWW[12] + 1./3.*aryuWW[320] + aryuWW[324];
   aryuWW[320]=aryuWW[38]*aryuWW[320];
   aryuWW[241]=aryuWW[320] + 7./3.*aryuWW[241] + 1./2.*aryuWW[267];
   aryuWW[241]=aryuWW[3]*aryuWW[241];
   aryuWW[206]=11 + aryuWW[206];
   aryuWW[206]=1./3.*aryuWW[206] - aryuWW[37];
   aryuWW[267]=aryuWW[3]*aryuWW[252];
   aryuWW[206]=1./2.*aryuWW[206] + 9*aryuWW[267];
   aryuWW[206]=aryuWW[3]*aryuWW[206];
   aryuWW[267]= - 19*aryuWW[56] - aryuWW[59];
   aryuWW[206]=1./6.*aryuWW[267] + aryuWW[206];
   aryuWW[206]=MMt*aryuWW[206];
   aryuWW[267]= - 85./8. - 7./3.*aryuWW[61];
   aryuWW[320]=1235./3. + 11*aryuWW[36];
   aryuWW[320]=aryuWW[12]*aryuWW[320];
   aryuWW[327]=31*aryuWW[12] - 1./8. + aryuWW[47];
   aryuWW[327]=1./4.*aryuWW[327] - aryuWW[37];
   aryuWW[327]=aryuWW[37]*aryuWW[327];
   aryuWW[203]=aryuWW[206] + aryuWW[241] + 1./6.*aryuWW[327] + 1./72.*
   aryuWW[320] + aryuWW[203] - 205./108.*aryuWW[17] - 7./12.*aryuWW[47]
    + 143./192.*aryuWW[50] - 5./192.*aryuWW[71] + 605./288.*aryuWW[69]
    - 1./12.*aryuWW[33] + 473./216.*aryuWW[63] - 271./192.*aryuWW[66]
    + 1./8.*aryuWW[267] + 2./3.*aryuWW[67];
   aryuWW[203]=MMt*aryuWW[203];
   aryuWW[206]=41./18.*aryuWW[63] - 331./144.*aryuWW[66] + 45./16. - 
   aryuWW[61];
   aryuWW[241]=43./2. + 7*aryuWW[36];
   aryuWW[241]=aryuWW[12]*aryuWW[241];
   aryuWW[267]=17*aryuWW[47];
   aryuWW[320]=25./16. + aryuWW[267];
   aryuWW[320]=aryuWW[37]*aryuWW[320];
   aryuWW[206]=1./6.*aryuWW[320] + 1./18.*aryuWW[241] - 11./36.*
   aryuWW[17] + 331./288.*aryuWW[50] + 137./96.*aryuWW[71] + 313./144.*
   aryuWW[69] + 1./2.*aryuWW[206] + aryuWW[229];
   aryuWW[241]=aryuWW[54] + 1./6.*aryuWW[59];
   aryuWW[241]=MMt*aryuWW[241];
   aryuWW[206]=1./2.*aryuWW[206] + 1./3.*aryuWW[241];
   aryuWW[206]=MMt*aryuWW[206];
   aryuWW[241]= - aryuWW[19]*aryuWW[111];
   aryuWW[320]=25./16.*aryuWW[37] + 25./12.*aryuWW[241] + 25./6.*
   aryuWW[12] - 281./48. + aryuWW[267];
   aryuWW[327]= - aryuWW[20]*aryuWW[111];
   aryuWW[329]= - aryuWW[38]*aryuWW[111];
   aryuWW[320]=25./8.*aryuWW[329] + 1./2.*aryuWW[320] + 25./3.*
   aryuWW[327];
   aryuWW[320]=aryuWW[38]*aryuWW[320];
   aryuWW[329]=17*aryuWW[36];
   aryuWW[330]=5*aryuWW[6];
   aryuWW[331]=aryuWW[330] - 2761./96. + aryuWW[329];
   aryuWW[331]=aryuWW[19]*aryuWW[331];
   aryuWW[235]=aryuWW[276] - 73./12. + aryuWW[235];
   aryuWW[235]=aryuWW[20]*aryuWW[235];
   aryuWW[276]=1./9.*aryuWW[74] - 13./128.*aryuWW[42];
   aryuWW[206]=1./2.*aryuWW[206] + 1./12.*aryuWW[320] + 1./72.*
   aryuWW[235] + 25./768.*aryuWW[170] + 1./72.*aryuWW[331] - 19./144.*
   aryuWW[23] + 679./2304.*aryuWW[22] + 797./1152.*aryuWW[51] - 7./72.*
   aryuWW[8] + 5*aryuWW[276] + 1./16.*aryuWW[72];
   aryuWW[206]=aryuWW[255]*aryuWW[206];
   aryuWW[235]= - 5*aryuWW[36];
   aryuWW[276]= - 71./18. + aryuWW[235];
   aryuWW[320]= - 71*aryuWW[6];
   aryuWW[276]=43*aryuWW[276] + aryuWW[320];
   aryuWW[276]=aryuWW[20]*aryuWW[276];
   aryuWW[331]= - 2191./8.*aryuWW[42] - 31*aryuWW[72];
   aryuWW[331]=1./2.*aryuWW[331] - 115*aryuWW[8];
   aryuWW[332]=17*aryuWW[6] - 965./432. + 49*aryuWW[36];
   aryuWW[332]=aryuWW[19]*aryuWW[332];
   aryuWW[250]=1./18.*aryuWW[276] + 391./288.*aryuWW[250] + 1./6.*
   aryuWW[332] - 175./54.*aryuWW[23] + 8515./864.*aryuWW[22] + 1./27.*
   aryuWW[331] + 179./16.*aryuWW[51];
   aryuWW[276]= - 11*aryuWW[35];
   aryuWW[331]=aryuWW[276] + aryuWW[322];
   aryuWW[332]=aryuWW[331] + 515./12.*aryuWW[111];
   aryuWW[332]=aryuWW[19]*aryuWW[332];
   aryuWW[332]= - 175./144.*aryuWW[37] + 1./3.*aryuWW[332] + 749./54.*
   aryuWW[12] + 2233./144. - aryuWW[47];
   aryuWW[333]=11*aryuWW[35];
   aryuWW[334]=aryuWW[333] - 485./4.*aryuWW[111];
   aryuWW[334]=aryuWW[38]*aryuWW[334];
   aryuWW[332]=1./6.*aryuWW[334] + 1./2.*aryuWW[332] + 235./9.*
   aryuWW[327];
   aryuWW[332]=aryuWW[38]*aryuWW[332];
   aryuWW[250]=1./2.*aryuWW[250] + aryuWW[332];
   aryuWW[332]=aryuWW[19] - 1./2.*aryuWW[38];
   aryuWW[332]=aryuWW[3]*aryuWW[38]*aryuWW[332];
   aryuWW[203]=aryuWW[206] + aryuWW[203] + 1./2.*aryuWW[250] + 17./3.*
   aryuWW[332];
   aryuWW[203]=aryuWW[255]*aryuWW[203];
   aryuWW[206]= - 3*aryuWW[65];
   aryuWW[250]= - 1121./324. + aryuWW[206];
   aryuWW[332]=3./2.*aryuWW[50];
   aryuWW[334]= - 13./6.*aryuWW[68];
   aryuWW[250]=5./27.*aryuWW[6] + 17./27.*aryuWW[36] - aryuWW[16] + 
   aryuWW[334] + aryuWW[332] + 1./2.*aryuWW[250] + aryuWW[60];
   aryuWW[273]=aryuWW[278] - 38./3. + aryuWW[273];
   aryuWW[273]=aryuWW[11]*aryuWW[273];
   aryuWW[296]=1./4.*aryuWW[296];
   aryuWW[250]=aryuWW[296] + 1./4.*aryuWW[250] + 1./27.*aryuWW[273];
   aryuWW[250]=MMH*aryuWW[250];
   aryuWW[273]=aryuWW[305] + 3*aryuWW[20];
   aryuWW[305]= - 15*aryuWW[38];
   aryuWW[273]=2*aryuWW[273] + aryuWW[305];
   aryuWW[273]=aryuWW[3]*aryuWW[38]*aryuWW[273];
   aryuWW[171]=aryuWW[203] + aryuWW[171] + aryuWW[250] + aryuWW[273] + 
   1./4.*aryuWW[309] + aryuWW[275];
   aryuWW[171]=aryuWW[255]*aryuWW[171];
   aryuWW[203]=aryuWW[256] + 11./3. + aryuWW[242];
   aryuWW[273]=aryuWW[203] + 1./2.*aryuWW[280];
   aryuWW[275]=aryuWW[330] - 22./3. + aryuWW[329];
   aryuWW[275]=aryuWW[3]*aryuWW[38]*aryuWW[275];
   aryuWW[273]=1./3.*aryuWW[273] + aryuWW[275];
   aryuWW[273]=MMt*aryuWW[273];
   aryuWW[203]=aryuWW[38]*aryuWW[203];
   aryuWW[203]=1./3.*aryuWW[203] + aryuWW[273];
   aryuWW[203]=aryuWW[34]*aryuWW[255]*aryuWW[203];
   aryuWW[273]= - 5 - aryuWW[36];
   aryuWW[273]=aryuWW[12]*aryuWW[273];
   aryuWW[273]=1./3.*aryuWW[273] + 1./3.*aryuWW[36] + 4./3.*aryuWW[17]
    - 4./3.*aryuWW[63] - 1 + 4./3.*aryuWW[64];
   aryuWW[273]=MMt*aryuWW[273];
   aryuWW[254]=5./3.*aryuWW[254] + aryuWW[273];
   aryuWW[171]=1./3.*aryuWW[203] + 64./3.*aryuWW[254] + aryuWW[171];
   aryuWW[171]=aryuWW[34]*aryuWW[171];
   aryuWW[125]=aryuWW[159] + aryuWW[136] + aryuWW[134] + aryuWW[133] + 
   aryuWW[152] - 47./2.*aryuWW[13] + aryuWW[130] + aryuWW[147] + 
   aryuWW[138] - 109./3. + aryuWW[125];
   aryuWW[125]=aryuWW[38]*aryuWW[125];
   aryuWW[130]= - 73./3. + aryuWW[181];
   aryuWW[130]=aryuWW[178] + 1./2.*aryuWW[130] - 32./3.*aryuWW[36];
   aryuWW[130]=aryuWW[20]*aryuWW[130];
   aryuWW[125]=aryuWW[125] + aryuWW[189] + aryuWW[130];
   aryuWW[125]=aryuWW[3]*aryuWW[125];
   aryuWW[130]= - 47./3.*aryuWW[13] + aryuWW[219] - 17303./144. + 
   aryuWW[213];
   aryuWW[130]=1./4.*aryuWW[130] + aryuWW[10];
   aryuWW[130]=aryuWW[243] + aryuWW[218] + aryuWW[217] + aryuWW[216] + 
   1./2.*aryuWW[130] + aryuWW[237];
   aryuWW[130]=aryuWW[37]*aryuWW[130];
   aryuWW[133]=aryuWW[285] - 194825./1152. - 68*aryuWW[64];
   aryuWW[133]=581./16.*aryuWW[63] + aryuWW[292] + aryuWW[290] + 1./3.*
   aryuWW[133] + aryuWW[286];
   aryuWW[134]=32*aryuWW[36];
   aryuWW[136]=43 + aryuWW[134];
   aryuWW[140]=32*aryuWW[7]*aryuWW[140];
   aryuWW[189]= - 11*aryuWW[6];
   aryuWW[136]=aryuWW[189] + 1./3.*aryuWW[136] + aryuWW[140];
   aryuWW[136]=aryuWW[40]*aryuWW[136];
   aryuWW[134]=59 + aryuWW[134];
   aryuWW[134]=1./3.*aryuWW[134] + aryuWW[140];
   aryuWW[134]=1./27.*aryuWW[134] - aryuWW[6];
   aryuWW[134]=aryuWW[1]*aryuWW[134];
   aryuWW[140]=47./12.*aryuWW[13];
   aryuWW[203]=107./4. + aryuWW[329];
   aryuWW[203]=aryuWW[12]*aryuWW[203];
   aryuWW[213]= - 32./9.*aryuWW[36] - 23./9. + aryuWW[258];
   aryuWW[213]=aryuWW[11]*aryuWW[213];
   aryuWW[125]=aryuWW[169] + aryuWW[284] + aryuWW[125] + aryuWW[163] + 
   aryuWW[141] + aryuWW[130] + aryuWW[213] + aryuWW[160] + 1./9.*
   aryuWW[136] + aryuWW[134] + aryuWW[224] + 25./12.*aryuWW[203] - 
   aryuWW[10] + 2723./324.*aryuWW[36] + aryuWW[140] + aryuWW[192] + 
   aryuWW[291] - 835./48.*aryuWW[17] + aryuWW[233] + aryuWW[303] + 
   aryuWW[262] + aryuWW[308] + aryuWW[302] + aryuWW[287] + aryuWW[301]
    + aryuWW[300] + aryuWW[297] + aryuWW[295] - aryuWW[62] + 1./3.*
   aryuWW[133] + aryuWW[307];
   aryuWW[125]=MMt*aryuWW[125];
   aryuWW[130]= - 97./6.*aryuWW[8] + aryuWW[73] - 833./12.*aryuWW[72]
    + aryuWW[304] + 3205./96.*aryuWW[42];
   aryuWW[133]=53*aryuWW[6] - 42739./192. + 47*aryuWW[36];
   aryuWW[133]=aryuWW[19]*aryuWW[133];
   aryuWW[134]=aryuWW[158] + 2381./96.*aryuWW[19];
   aryuWW[134]=aryuWW[37]*aryuWW[134];
   aryuWW[136]=aryuWW[320] - 8533./24. + 403*aryuWW[36];
   aryuWW[136]=1./9.*aryuWW[136] + 229./2.*aryuWW[37];
   aryuWW[136]=aryuWW[20]*aryuWW[136];
   aryuWW[130]=1./3.*aryuWW[136] + 1./2.*aryuWW[134] + 1./9.*
   aryuWW[133] + aryuWW[293] + 403./12.*aryuWW[23] + aryuWW[139] + 1087.
   /192.*aryuWW[22] + 17695./288.*aryuWW[51] + 1./3.*aryuWW[130] + 
   aryuWW[315];
   aryuWW[133]=997./1296. + aryuWW[47];
   aryuWW[133]=47./3.*aryuWW[13] + 1./2.*aryuWW[133] + aryuWW[219];
   aryuWW[134]= - 32*aryuWW[7];
   aryuWW[136]=aryuWW[189] + 43./3. + aryuWW[134];
   aryuWW[136]=aryuWW[40]*aryuWW[136];
   aryuWW[169]=aryuWW[313] + 7*aryuWW[35] + aryuWW[312];
   aryuWW[169]=aryuWW[38]*aryuWW[169];
   aryuWW[134]=59./3. + aryuWW[134];
   aryuWW[134]=1./27.*aryuWW[134] - aryuWW[6];
   aryuWW[134]=aryuWW[1]*aryuWW[134];
   aryuWW[189]=aryuWW[323] + 77./16.*aryuWW[111];
   aryuWW[189]=aryuWW[19]*aryuWW[189];
   aryuWW[203]=aryuWW[326] + 121./2.*aryuWW[111];
   aryuWW[203]=aryuWW[20]*aryuWW[203];
   aryuWW[133]=1./2.*aryuWW[169] + 1./2.*aryuWW[203] + aryuWW[328] - 32.
   /9.*aryuWW[11] + 1./2.*aryuWW[189] + 1./9.*aryuWW[136] + aryuWW[134]
    + 2873./48.*aryuWW[12] + 1./4.*aryuWW[133] - aryuWW[10];
   aryuWW[133]=aryuWW[38]*aryuWW[133];
   aryuWW[134]=13./2.*aryuWW[36];
   aryuWW[136]=43 + aryuWW[134];
   aryuWW[136]=1./3.*aryuWW[136] + aryuWW[256];
   aryuWW[169]=5./6.*aryuWW[6] - 203./9. - 37./2.*aryuWW[36];
   aryuWW[169]=aryuWW[37]*aryuWW[169];
   aryuWW[136]=aryuWW[275] + 1./3.*aryuWW[136] + 1./2.*aryuWW[169];
   aryuWW[136]=MMt*aryuWW[136];
   aryuWW[169]= - 32*aryuWW[37] + aryuWW[256] + 43./3. + aryuWW[242];
   aryuWW[169]=aryuWW[38]*aryuWW[169];
   aryuWW[136]=1./3.*aryuWW[169] + aryuWW[136];
   aryuWW[136]=aryuWW[34]*aryuWW[136];
   aryuWW[169]=aryuWW[266] + aryuWW[204];
   aryuWW[169]=2./3.*aryuWW[169] + aryuWW[305];
   aryuWW[169]=aryuWW[3]*aryuWW[38]*aryuWW[169];
   aryuWW[125]=1./3.*aryuWW[136] + aryuWW[125] + aryuWW[250] + 
   aryuWW[169] + 1./4.*aryuWW[130] + aryuWW[133];
   aryuWW[125]=aryuWW[34]*aryuWW[125];
   aryuWW[130]= - 1369./27. + aryuWW[234];
   aryuWW[133]=7./2.*aryuWW[46];
   aryuWW[130]= - 47./18.*aryuWW[13] + aryuWW[133] + 1./4.*aryuWW[130]
    + aryuWW[172];
   aryuWW[136]=323./144.*aryuWW[12];
   aryuWW[130]=aryuWW[306] + aryuWW[136] + 1./2.*aryuWW[130] + 
   aryuWW[199];
   aryuWW[130]=aryuWW[11]*aryuWW[130];
   aryuWW[169]=1./3. - 5./2.*aryuWW[10];
   aryuWW[169]=1./3.*aryuWW[169] + aryuWW[228];
   aryuWW[169]=aryuWW[12]*aryuWW[169];
   aryuWW[189]= - 10*aryuWW[102];
   aryuWW[203]= - 811./576. + aryuWW[189];
   aryuWW[204]= - 9./8.*aryuWW[95];
   aryuWW[213]=1./2.*aryuWW[94];
   aryuWW[216]= - 2./3.*aryuWW[99];
   aryuWW[217]= - 11./6.*aryuWW[97];
   aryuWW[218]=9./16.*aryuWW[93];
   aryuWW[228]=9./16.*aryuWW[92];
   aryuWW[237]= - 25./48.*aryuWW[98];
   aryuWW[242]= - 89./16.*aryuWW[103];
   aryuWW[243]= - 11./12.*aryuWW[96];
   aryuWW[250]=197./24.*aryuWW[100];
   aryuWW[254]=1./8.*aryuWW[17];
   aryuWW[256]=215./48.*aryuWW[16];
   aryuWW[266]= - 5./8.*aryuWW[46];
   aryuWW[273]=7./6.*aryuWW[10];
   aryuWW[275]=4./3.*aryuWW[89] - aryuWW[86];
   aryuWW[275]=23./24.*aryuWW[82] - 17./24.*aryuWW[85] - 17./48.*
   aryuWW[109] - 1./6.*aryuWW[81] + 2*aryuWW[275] - 3./4.*aryuWW[84];
   aryuWW[275]=MMH*aryuWW[275];
   aryuWW[280]=23./9.*aryuWW[101];
   aryuWW[285]=1./6.*aryuWW[91];
   aryuWW[286]= - 9./8.*aryuWW[44];
   aryuWW[290]= - 3./4.*aryuWW[48];
   aryuWW[130]=aryuWW[275] + aryuWW[130] + aryuWW[169] + aryuWW[273] - 
   47./72.*aryuWW[13] + aryuWW[266] + aryuWW[256] + aryuWW[290] + 
   aryuWW[254] + aryuWW[286] + aryuWW[250] + aryuWW[243] + aryuWW[242]
    + aryuWW[285] + aryuWW[237] + aryuWW[228] + aryuWW[218] + 
   aryuWW[217] + aryuWW[216] + aryuWW[280] + aryuWW[213] + 1./3.*
   aryuWW[203] + aryuWW[204];
   aryuWW[130]=MMH*aryuWW[130];
   aryuWW[203]= - 45*aryuWW[44];
   aryuWW[292]= - 20033./18. + aryuWW[203];
   aryuWW[293]= - 7./4.*aryuWW[46];
   aryuWW[295]= - 49./6.*aryuWW[10];
   aryuWW[297]= - aryuWW[18]*aryuWW[112];
   aryuWW[300]=2*aryuWW[297];
   aryuWW[301]=3*aryuWW[113] - aryuWW[112];
   aryuWW[301]=5./4.*aryuWW[19]*aryuWW[301];
   aryuWW[302]= - aryuWW[20]*aryuWW[112];
   aryuWW[303]=17./4.*aryuWW[302];
   aryuWW[304]= - 2*aryuWW[48];
   aryuWW[292]=aryuWW[303] + 13./6.*aryuWW[11] + aryuWW[301] + 
   aryuWW[300] - 10529./72.*aryuWW[12] + aryuWW[295] - 3523./72.*
   aryuWW[13] + aryuWW[293] + 1./8.*aryuWW[292] + aryuWW[304];
   aryuWW[292]=aryuWW[20]*aryuWW[292];
   aryuWW[305]=119./24.*aryuWW[12];
   aryuWW[140]=aryuWW[305] - aryuWW[10] + aryuWW[140] + aryuWW[173] + 
   aryuWW[48] - 241./9. + aryuWW[316];
   aryuWW[309]=2*aryuWW[122];
   aryuWW[140]=1./2.*aryuWW[140] + aryuWW[309];
   aryuWW[140]=aryuWW[18]*aryuWW[140];
   aryuWW[312]= - 167./6.*aryuWW[13] + aryuWW[277] + aryuWW[232] + 3325.
   /27. + aryuWW[270];
   aryuWW[313]=4*aryuWW[10];
   aryuWW[315]= - aryuWW[18]*aryuWW[113];
   aryuWW[320]=7./6.*aryuWW[315];
   aryuWW[323]=73./3.*aryuWW[113] - 5*aryuWW[112];
   aryuWW[323]=1./8.*aryuWW[19]*aryuWW[323];
   aryuWW[312]=aryuWW[323] + aryuWW[320] - 4435./48.*aryuWW[12] + 1./4.
   *aryuWW[312] + aryuWW[313];
   aryuWW[312]=aryuWW[19]*aryuWW[312];
   aryuWW[198]=aryuWW[198] + 13./2.*aryuWW[19];
   aryuWW[198]=aryuWW[19]*aryuWW[198];
   aryuWW[162]= - 2*aryuWW[162] + 1./4.*aryuWW[198];
   aryuWW[198]= - 41./3.*aryuWW[18] + 67./2.*aryuWW[19];
   aryuWW[198]=1./4.*aryuWW[198];
   aryuWW[328]=aryuWW[198] + 21*aryuWW[20];
   aryuWW[328]=aryuWW[20]*aryuWW[328];
   aryuWW[328]=aryuWW[162] + aryuWW[328];
   aryuWW[328]=aryuWW[3]*aryuWW[328];
   aryuWW[335]= - 9./8.*aryuWW[77];
   aryuWW[336]=61./24.*aryuWW[105];
   aryuWW[337]=605./48.*aryuWW[106];
   aryuWW[338]= - 47./12.*aryuWW[79];
   aryuWW[339]= - 29./6.*aryuWW[78];
   aryuWW[340]=269./144.*aryuWW[21];
   aryuWW[341]= - aryuWW[18] + 55./2.*aryuWW[19];
   aryuWW[341]=aryuWW[11]*aryuWW[341];
   aryuWW[130]=aryuWW[130] + aryuWW[328] + aryuWW[292] + 1./6.*
   aryuWW[341] + aryuWW[312] + aryuWW[140] - 6607./144.*aryuWW[23] + 
   aryuWW[340] + 599./16.*aryuWW[22] + aryuWW[339] - 179./24.*
   aryuWW[107] + aryuWW[338] - 149./24.*aryuWW[80] + aryuWW[337] + 
   aryuWW[336] + 922./9.*aryuWW[108] + aryuWW[335];
   aryuWW[130]=aryuWW[76]*aryuWW[130];
   aryuWW[140]=1 - 3*aryuWW[6];
   aryuWW[292]=aryuWW[40]*aryuWW[140];
   aryuWW[117]=aryuWW[159] + 2*aryuWW[292] + 2*aryuWW[259] - 33*
   aryuWW[12] + aryuWW[152] + 165./2.*aryuWW[13] + aryuWW[324] + 
   aryuWW[147] + aryuWW[138] + 185./3. + aryuWW[117];
   aryuWW[117]=aryuWW[38]*aryuWW[117];
   aryuWW[138]=1./2.*aryuWW[42];
   aryuWW[147]= - 2*aryuWW[51];
   aryuWW[159]= - 1./2.*aryuWW[22] + aryuWW[147] + aryuWW[138] + 
   aryuWW[43];
   aryuWW[181]=13 + aryuWW[181];
   aryuWW[292]=aryuWW[178] + 1./2.*aryuWW[181] + 8*aryuWW[36];
   aryuWW[292]=aryuWW[20]*aryuWW[292];
   aryuWW[181]=aryuWW[181] + 13*aryuWW[36];
   aryuWW[181]=aryuWW[19]*aryuWW[181];
   aryuWW[117]=aryuWW[117] + aryuWW[292] + aryuWW[142] + 1./4.*
   aryuWW[181] + aryuWW[188] + aryuWW[222] + 3*aryuWW[159] + 
   aryuWW[190];
   aryuWW[117]=aryuWW[3]*aryuWW[117];
   aryuWW[142]=aryuWW[252] + 2*aryuWW[257];
   aryuWW[142]=aryuWW[3]*aryuWW[142];
   aryuWW[159]= - aryuWW[49] + 9 - aryuWW[47];
   aryuWW[135]=3*aryuWW[142] + aryuWW[135] + aryuWW[272] + 1./2.*
   aryuWW[159] + aryuWW[39];
   aryuWW[135]=aryuWW[3]*aryuWW[135];
   aryuWW[142]= - aryuWW[56] - aryuWW[54];
   aryuWW[159]=aryuWW[142] - 1./2.*aryuWW[59];
   aryuWW[135]=3*aryuWW[135] + aryuWW[283] + aryuWW[282] + aryuWW[191]
    + 3./2.*aryuWW[159] - aryuWW[55];
   aryuWW[135]=MMt*aryuWW[135];
   aryuWW[159]=55*aryuWW[13] + 5551./144. + aryuWW[219];
   aryuWW[123]= - aryuWW[37] + aryuWW[123] + aryuWW[223] + 1./3.*
   aryuWW[259] + 15./4.*aryuWW[12] + 1./4.*aryuWW[159] + aryuWW[10];
   aryuWW[123]=aryuWW[37]*aryuWW[123];
   aryuWW[159]= - 11 - 8*aryuWW[36];
   aryuWW[181]=aryuWW[7]*aryuWW[143];
   aryuWW[159]=1./3.*aryuWW[159] + 8*aryuWW[181];
   aryuWW[159]=1./3.*aryuWW[159] + aryuWW[6];
   aryuWW[181]=aryuWW[1]*aryuWW[159];
   aryuWW[159]=aryuWW[40]*aryuWW[159];
   aryuWW[188]=8./3.*aryuWW[36];
   aryuWW[190]=aryuWW[188] + 11./3. + aryuWW[258];
   aryuWW[190]=aryuWW[11]*aryuWW[190];
   aryuWW[191]= - 91./2. - 55*aryuWW[36];
   aryuWW[191]=aryuWW[12]*aryuWW[191];
   aryuWW[117]=aryuWW[135] + aryuWW[284] + aryuWW[117] + aryuWW[163] + 
   aryuWW[141] + 1./2.*aryuWW[123] + aryuWW[190] + aryuWW[160] + 
   aryuWW[159] + 1./3.*aryuWW[181] + aryuWW[224] + 3./8.*aryuWW[191] - 
   aryuWW[10] - 1025./108.*aryuWW[36] - 55./4.*aryuWW[13] + aryuWW[192]
    + aryuWW[291] - 93./16.*aryuWW[17] + aryuWW[233] + 3./4.*aryuWW[47]
    + aryuWW[262] + aryuWW[308] - 185./128.*aryuWW[50] + aryuWW[287] + 
   1663./128.*aryuWW[71] + 25./64.*aryuWW[69] + 7./4.*aryuWW[30] + 
   aryuWW[33] - aryuWW[62] + aryuWW[307] - 31./16.*aryuWW[63] + 41./128.
   *aryuWW[66] + 11./8.*aryuWW[67] - 145765./3456. + 6*aryuWW[61];
   aryuWW[117]=MMt*aryuWW[117];
   aryuWW[123]=3*aryuWW[47];
   aryuWW[135]= - 39583./432. + aryuWW[123];
   aryuWW[135]= - 55*aryuWW[13] + 1./2.*aryuWW[135] + aryuWW[219];
   aryuWW[141]= - 11./3. + 8*aryuWW[7];
   aryuWW[141]=1./3.*aryuWW[141] + aryuWW[6];
   aryuWW[159]=aryuWW[1]*aryuWW[141];
   aryuWW[141]=aryuWW[40]*aryuWW[141];
   aryuWW[160]= - 9./8.*aryuWW[111] - aryuWW[35] + 3./2.*aryuWW[110];
   aryuWW[160]=aryuWW[19]*aryuWW[160];
   aryuWW[163]=aryuWW[326] - 33./2.*aryuWW[111];
   aryuWW[163]=aryuWW[20]*aryuWW[163];
   aryuWW[181]= - 123./8.*aryuWW[111] + aryuWW[35] - 69./2.*aryuWW[110]
   ;
   aryuWW[181]=aryuWW[38]*aryuWW[181];
   aryuWW[135]=1./4.*aryuWW[181] + 1./2.*aryuWW[163] + 271./128.*
   aryuWW[37] + 8./3.*aryuWW[11] + 1./4.*aryuWW[160] + aryuWW[141] + 1./
   3.*aryuWW[159] - 75./16.*aryuWW[12] + 1./4.*aryuWW[135] - aryuWW[10]
   ;
   aryuWW[135]=aryuWW[38]*aryuWW[135];
   aryuWW[141]=aryuWW[144] + 3./2. + 5./3.*aryuWW[36];
   aryuWW[141]=aryuWW[18]*aryuWW[141];
   aryuWW[144]=13*aryuWW[6] - 7831./96. - 7*aryuWW[36];
   aryuWW[144]=aryuWW[19]*aryuWW[144];
   aryuWW[158]=aryuWW[158] - 3341./96.*aryuWW[19];
   aryuWW[158]=aryuWW[37]*aryuWW[158];
   aryuWW[159]= - 211*aryuWW[6] - 2851./4. - 263*aryuWW[36];
   aryuWW[159]=1./9.*aryuWW[159] - 17*aryuWW[37];
   aryuWW[159]=aryuWW[20]*aryuWW[159];
   aryuWW[139]=1./2.*aryuWW[159] + 1./2.*aryuWW[158] + 1./6.*
   aryuWW[144] + aryuWW[141] + 155./4.*aryuWW[23] + aryuWW[139] - 1373./
   64.*aryuWW[22] + 3777./32.*aryuWW[51] - 19./2.*aryuWW[43] - 7./2.*
   aryuWW[8] + 1./3.*aryuWW[73] - 159./4.*aryuWW[72] + 23*aryuWW[74] - 
   41./32.*aryuWW[42];
   aryuWW[141]= - 89./36. + aryuWW[206];
   aryuWW[141]= - 1./9.*aryuWW[6] - 5./9.*aryuWW[36] - aryuWW[16] + 
   aryuWW[334] + aryuWW[332] + 1./2.*aryuWW[141] + aryuWW[60];
   aryuWW[144]=aryuWW[126] - 2 + aryuWW[176];
   aryuWW[144]=aryuWW[11]*aryuWW[144];
   aryuWW[141]=aryuWW[296] + 1./4.*aryuWW[141] + 1./9.*aryuWW[144];
   aryuWW[141]=MMH*aryuWW[141];
   aryuWW[144]= - 11 + aryuWW[271];
   aryuWW[144]=1./3.*aryuWW[144] + aryuWW[115];
   aryuWW[158]= - 1./6.*aryuWW[6] + 17./3. + 9./2.*aryuWW[36];
   aryuWW[158]=aryuWW[37]*aryuWW[158];
   aryuWW[159]= - aryuWW[6] + 2 + aryuWW[235];
   aryuWW[159]=aryuWW[3]*aryuWW[38]*aryuWW[159];
   aryuWW[144]=aryuWW[159] + 1./3.*aryuWW[144] + 1./2.*aryuWW[158];
   aryuWW[144]=MMt*aryuWW[144];
   aryuWW[158]=5./2.*aryuWW[36];
   aryuWW[159]=8*aryuWW[37] + aryuWW[115] - 11./3. + aryuWW[158];
   aryuWW[159]=aryuWW[38]*aryuWW[159];
   aryuWW[144]=1./3.*aryuWW[159] + aryuWW[144];
   aryuWW[144]=aryuWW[34]*aryuWW[144];
   aryuWW[159]=aryuWW[19] + 2*aryuWW[20];
   aryuWW[159]=7*aryuWW[159] - 3./2.*aryuWW[38];
   aryuWW[159]=aryuWW[3]*aryuWW[38]*aryuWW[159];
   aryuWW[117]=aryuWW[144] + aryuWW[117] + aryuWW[141] + aryuWW[159] + 
   1./4.*aryuWW[139] + aryuWW[135];
   aryuWW[117]=aryuWW[34]*aryuWW[117];
   aryuWW[135]=1931./27. + aryuWW[234];
   aryuWW[135]=55./6.*aryuWW[13] + aryuWW[173] + 1./4.*aryuWW[135] + 
   aryuWW[172];
   aryuWW[135]=aryuWW[306] + 95./48.*aryuWW[12] + 1./2.*aryuWW[135] + 
   aryuWW[199];
   aryuWW[135]=aryuWW[11]*aryuWW[135];
   aryuWW[139]= - 11./24.*aryuWW[109];
   aryuWW[141]= - 11./12.*aryuWW[85];
   aryuWW[144]= - 11./12.*aryuWW[82] + aryuWW[141] + aryuWW[139] - 3./2.
   *aryuWW[84] - 1./3.*aryuWW[81];
   aryuWW[144]=MMH*aryuWW[144];
   aryuWW[159]=10895./216. - 9*aryuWW[95];
   aryuWW[160]=5./2. + 11*aryuWW[10];
   aryuWW[160]=1./2.*aryuWW[160] + aryuWW[12];
   aryuWW[160]=aryuWW[12]*aryuWW[160];
   aryuWW[135]=1./2.*aryuWW[144] + aryuWW[135] + 1./6.*aryuWW[160] + 1./
   24.*aryuWW[10] + 55./24.*aryuWW[13] - 3./8.*aryuWW[46] + 137./48.*
   aryuWW[16] + aryuWW[290] - 15./8.*aryuWW[17] + aryuWW[286] + 95./8.*
   aryuWW[100] - 5./6.*aryuWW[96] + aryuWW[285] + 29./48.*aryuWW[98] + 
   29./48.*aryuWW[92] + 29./48.*aryuWW[93] + aryuWW[194] + 1./3.*
   aryuWW[99] + aryuWW[280] + 1./8.*aryuWW[159] + 1./3.*aryuWW[94];
   aryuWW[135]=MMH*aryuWW[135];
   aryuWW[144]=15377./27. + aryuWW[203];
   aryuWW[159]=7./2.*aryuWW[113] - aryuWW[112];
   aryuWW[159]=aryuWW[19]*aryuWW[159];
   aryuWW[144]=1./3.*aryuWW[302] - 101./12.*aryuWW[11] + 1./2.*
   aryuWW[159] + 23./12.*aryuWW[297] + 1011./16.*aryuWW[12] - 53./12.*
   aryuWW[10] + 2893./48.*aryuWW[13] - aryuWW[46] + 1./8.*aryuWW[144]
    + aryuWW[304];
   aryuWW[144]=aryuWW[20]*aryuWW[144];
   aryuWW[159]=125./16.*aryuWW[12] + aryuWW[182] - 55./8.*aryuWW[13] + 
   aryuWW[231] + aryuWW[317] - 245./9. + aryuWW[251];
   aryuWW[159]=aryuWW[18]*aryuWW[159];
   aryuWW[160]=77./2.*aryuWW[13] + aryuWW[277] + aryuWW[232] - 66119./
   27. + aryuWW[270];
   aryuWW[163]=9*aryuWW[10];
   aryuWW[160]=551./3.*aryuWW[12] + 1./2.*aryuWW[160] + aryuWW[163];
   aryuWW[176]=19./2.*aryuWW[113] - aryuWW[112];
   aryuWW[176]=aryuWW[19]*aryuWW[176];
   aryuWW[160]=1./2.*aryuWW[176] + 1./2.*aryuWW[160] + 3*aryuWW[315];
   aryuWW[160]=aryuWW[19]*aryuWW[160];
   aryuWW[176]=aryuWW[310] + 5*aryuWW[19];
   aryuWW[176]=aryuWW[19]*aryuWW[176];
   aryuWW[132]= - 43./4.*aryuWW[20] - 3./4.*aryuWW[18] + aryuWW[132];
   aryuWW[132]=aryuWW[20]*aryuWW[132];
   aryuWW[132]=1./4.*aryuWW[176] + aryuWW[132];
   aryuWW[132]=aryuWW[3]*aryuWW[132];
   aryuWW[176]=489./2.*aryuWW[80] + 197./2.*aryuWW[106] - 9*aryuWW[77]
    + 61./3.*aryuWW[105];
   aryuWW[176]= - 2397./8.*aryuWW[23] + 863./72.*aryuWW[21] - 1571./4.*
   aryuWW[22] - 11./2.*aryuWW[78] + 159*aryuWW[107] + 1./4.*aryuWW[176]
    - 17./3.*aryuWW[79];
   aryuWW[181]=aryuWW[185] - 26*aryuWW[19];
   aryuWW[181]=aryuWW[11]*aryuWW[181];
   aryuWW[132]=aryuWW[135] + aryuWW[132] + aryuWW[144] + 1./3.*
   aryuWW[181] + 1./2.*aryuWW[160] + 1./2.*aryuWW[176] + aryuWW[159];
   aryuWW[132]=aryuWW[76]*aryuWW[132];
   aryuWW[135]= - 1 - 33*aryuWW[13];
   aryuWW[135]=aryuWW[164] + 1./2.*aryuWW[135] + 2*aryuWW[10];
   aryuWW[144]= - aryuWW[7] + aryuWW[6];
   aryuWW[159]=aryuWW[1]*aryuWW[144];
   aryuWW[144]=aryuWW[40]*aryuWW[144];
   aryuWW[135]= - 6*aryuWW[11] + 6*aryuWW[144] + 3*aryuWW[135] + 2*
   aryuWW[159];
   aryuWW[135]=aryuWW[38]*aryuWW[135];
   aryuWW[160]=aryuWW[19]*aryuWW[36];
   aryuWW[160]=aryuWW[160] + aryuWW[170];
   aryuWW[170]=aryuWW[36] - aryuWW[37];
   aryuWW[176]=aryuWW[20]*aryuWW[170];
   aryuWW[160]=1./2.*aryuWW[160] + aryuWW[176];
   aryuWW[135]=3./2.*aryuWW[160] + aryuWW[135];
   aryuWW[135]=aryuWW[3]*aryuWW[135];
   aryuWW[160]=33*aryuWW[13] - 59./9.*aryuWW[36];
   aryuWW[176]=aryuWW[12]*aryuWW[180];
   aryuWW[179]=aryuWW[7]*aryuWW[179];
   aryuWW[179]= - aryuWW[6] - 1./6.*aryuWW[36] + aryuWW[179];
   aryuWW[180]=aryuWW[1]*aryuWW[179];
   aryuWW[179]=aryuWW[40]*aryuWW[179];
   aryuWW[181]= - MMt*aryuWW[3]*aryuWW[37];
   aryuWW[135]=3./2.*aryuWW[181] + aryuWW[135] + aryuWW[161] + 
   aryuWW[236] + aryuWW[179] + 1./3.*aryuWW[180] + 33./4.*aryuWW[176]
    + 1./4.*aryuWW[160] - aryuWW[10];
   aryuWW[135]=MMt*aryuWW[135];
   aryuWW[160]=aryuWW[11]*aryuWW[205];
   aryuWW[161]= - 2*aryuWW[11] - 1./2. + aryuWW[10];
   aryuWW[161]=aryuWW[37]*aryuWW[161];
   aryuWW[160]=aryuWW[161] + aryuWW[160] + aryuWW[186] + 1./4.*
   aryuWW[36] + aryuWW[196];
   aryuWW[160]=MMH*aryuWW[160];
   aryuWW[161]= - 33./4.*aryuWW[12];
   aryuWW[176]=aryuWW[161] + 33./4.*aryuWW[13] - aryuWW[10];
   aryuWW[179]=aryuWW[11] + aryuWW[183] + aryuWW[176] + aryuWW[211];
   aryuWW[179]=aryuWW[38]*aryuWW[179];
   aryuWW[180]=aryuWW[115] + 7./3. + aryuWW[174];
   aryuWW[180]=1./2.*aryuWW[180] - aryuWW[37];
   aryuWW[180]=aryuWW[37]*aryuWW[180];
   aryuWW[181]=aryuWW[36] + aryuWW[6];
   aryuWW[182]=aryuWW[181] - 2*aryuWW[37];
   aryuWW[182]=aryuWW[3]*aryuWW[38]*aryuWW[182];
   aryuWW[126]=3*aryuWW[182] + aryuWW[180] - 2./3.*aryuWW[36] + 
   aryuWW[126];
   aryuWW[126]=MMt*aryuWW[126];
   aryuWW[180]=1./2.*aryuWW[187] + aryuWW[37];
   aryuWW[182]=aryuWW[38]*aryuWW[180];
   aryuWW[126]=aryuWW[182] + aryuWW[126];
   aryuWW[126]=aryuWW[34]*aryuWW[126];
   aryuWW[158]=aryuWW[278] + 17./9. + aryuWW[158];
   aryuWW[158]=aryuWW[19]*aryuWW[158];
   aryuWW[158]=aryuWW[193] + aryuWW[158];
   aryuWW[182]=aryuWW[289] - 8./3.*aryuWW[19];
   aryuWW[182]=aryuWW[37]*aryuWW[182];
   aryuWW[134]=13./2.*aryuWW[6] - 17./3. + aryuWW[134];
   aryuWW[134]=1./4.*aryuWW[134] + aryuWW[37];
   aryuWW[134]=aryuWW[20]*aryuWW[134];
   aryuWW[126]=aryuWW[126] + aryuWW[135] + 1./3.*aryuWW[160] + 
   aryuWW[179] + 1./3.*aryuWW[134] + 1./4.*aryuWW[158] + aryuWW[182];
   aryuWW[126]=aryuWW[34]*aryuWW[126];
   aryuWW[134]= - 11*aryuWW[13];
   aryuWW[135]=127./27. + aryuWW[134];
   aryuWW[158]=11./2.*aryuWW[12];
   aryuWW[135]=aryuWW[299] + aryuWW[158] + 1./4.*aryuWW[135] + 2./3.*
   aryuWW[10];
   aryuWW[135]=aryuWW[11]*aryuWW[135];
   aryuWW[160]= - 11./2.*aryuWW[13];
   aryuWW[179]=aryuWW[12]*aryuWW[184];
   aryuWW[179]=11*aryuWW[179] + aryuWW[160] - 127./27.*aryuWW[10];
   aryuWW[135]=1./4.*aryuWW[179] + aryuWW[135];
   aryuWW[135]=MMH*aryuWW[135];
   aryuWW[179]=1057./27. - 165./2.*aryuWW[13];
   aryuWW[182]=7*aryuWW[10];
   aryuWW[179]=1./2.*aryuWW[179] + aryuWW[182];
   aryuWW[184]=aryuWW[113] - aryuWW[112];
   aryuWW[186]=aryuWW[19]*aryuWW[184];
   aryuWW[179]=3./8.*aryuWW[186] + 1./4.*aryuWW[179] + 22*aryuWW[12];
   aryuWW[179]=aryuWW[19]*aryuWW[179];
   aryuWW[190]= - 1057./27. - 143./2.*aryuWW[13];
   aryuWW[186]=3*aryuWW[186] - 11*aryuWW[12] + 1./2.*aryuWW[190] + 25./
   3.*aryuWW[10];
   aryuWW[186]=1./4.*aryuWW[186] + aryuWW[299];
   aryuWW[186]=aryuWW[20]*aryuWW[186];
   aryuWW[190]=aryuWW[265] + aryuWW[20];
   aryuWW[190]=aryuWW[20]*aryuWW[190];
   aryuWW[190]= - 1./2.*aryuWW[279] + aryuWW[190];
   aryuWW[190]=aryuWW[3]*aryuWW[190];
   aryuWW[191]=aryuWW[18]*aryuWW[176];
   aryuWW[192]=aryuWW[18] - 19./3.*aryuWW[19];
   aryuWW[192]=aryuWW[11]*aryuWW[192];
   aryuWW[135]=aryuWW[135] + 17./4.*aryuWW[190] + aryuWW[186] + 1./2.*
   aryuWW[192] + 1./2.*aryuWW[191] + aryuWW[179];
   aryuWW[135]=aryuWW[76]*aryuWW[135];
   aryuWW[179]=aryuWW[10]*aryuWW[220];
   aryuWW[115]=aryuWW[115] - 1./2.*aryuWW[7] + aryuWW[179];
   aryuWW[179]=aryuWW[1]*aryuWW[115];
   aryuWW[115]=aryuWW[40]*aryuWW[115];
   aryuWW[186]= - 2*aryuWW[7];
   aryuWW[190]=aryuWW[6] + 1./3. + aryuWW[186];
   aryuWW[191]=aryuWW[1]*aryuWW[190];
   aryuWW[190]=aryuWW[40]*aryuWW[190];
   aryuWW[190]=1./3.*aryuWW[191] + aryuWW[190];
   aryuWW[190]=aryuWW[11]*aryuWW[190];
   aryuWW[115]=aryuWW[190] + 1./3.*aryuWW[179] + aryuWW[115];
   aryuWW[115]=MMH*aryuWW[115];
   aryuWW[179]=17./12. - 8*aryuWW[7];
   aryuWW[179]=1./3.*aryuWW[179] + aryuWW[261];
   aryuWW[190]=aryuWW[1]*aryuWW[179];
   aryuWW[179]=aryuWW[40]*aryuWW[179];
   aryuWW[179]=1./3.*aryuWW[190] + aryuWW[179];
   aryuWW[179]=aryuWW[19]*aryuWW[179];
   aryuWW[190]=13./4.*aryuWW[6] - 17./12. + aryuWW[7];
   aryuWW[191]=aryuWW[1]*aryuWW[190];
   aryuWW[190]=aryuWW[40]*aryuWW[190];
   aryuWW[190]=1./3.*aryuWW[191] + aryuWW[190];
   aryuWW[190]=aryuWW[20]*aryuWW[190];
   aryuWW[115]=aryuWW[126] + aryuWW[135] + 1./3.*aryuWW[115] + 1./3.*
   aryuWW[190] + 1./2.*aryuWW[214] + aryuWW[179];
   aryuWW[115]=aryuWW[165]*aryuWW[115];
   aryuWW[126]= - 1./2.*aryuWW[8];
   aryuWW[135]= - aryuWW[9] + aryuWW[126];
   aryuWW[179]= - 1./2.*aryuWW[21];
   aryuWW[135]= - 35./2.*aryuWW[23] + aryuWW[179] + 11*aryuWW[135] - 13
   *aryuWW[22];
   aryuWW[190]=aryuWW[1]*aryuWW[135];
   aryuWW[135]=aryuWW[40]*aryuWW[135];
   aryuWW[191]=1./2. - aryuWW[7];
   aryuWW[191]=1./2.*aryuWW[191] + aryuWW[6];
   aryuWW[192]=aryuWW[1]*aryuWW[191];
   aryuWW[191]=aryuWW[40]*aryuWW[191];
   aryuWW[191]=1./3.*aryuWW[192] + aryuWW[191];
   aryuWW[191]=aryuWW[18]*aryuWW[191];
   aryuWW[192]=619./6. - 71*aryuWW[7];
   aryuWW[192]=1./3.*aryuWW[192] + aryuWW[6];
   aryuWW[193]=aryuWW[1]*aryuWW[192];
   aryuWW[192]=aryuWW[40]*aryuWW[192];
   aryuWW[192]=1./3.*aryuWW[193] + aryuWW[192];
   aryuWW[192]=aryuWW[19]*aryuWW[192];
   aryuWW[135]=1./2.*aryuWW[192] + aryuWW[191] + 1./3.*aryuWW[190] + 
   aryuWW[135];
   aryuWW[190]=1./6.*aryuWW[14];
   aryuWW[191]= - aryuWW[16] + aryuWW[190] + aryuWW[28] - 7./9. - 
   aryuWW[31];
   aryuWW[192]= - 1./3.*aryuWW[6];
   aryuWW[191]=1./2.*aryuWW[191] + aryuWW[192];
   aryuWW[193]=aryuWW[1]*aryuWW[191];
   aryuWW[191]=aryuWW[40]*aryuWW[191];
   aryuWW[191]=1./3.*aryuWW[193] + aryuWW[191];
   aryuWW[193]= - 1./4.*aryuWW[7];
   aryuWW[192]=aryuWW[192] - 2./9. + aryuWW[193];
   aryuWW[194]=aryuWW[1]*aryuWW[192];
   aryuWW[192]=aryuWW[40]*aryuWW[192];
   aryuWW[192]=1./3.*aryuWW[194] + aryuWW[192];
   aryuWW[192]=aryuWW[11]*aryuWW[192];
   aryuWW[191]=1./2.*aryuWW[191] + aryuWW[192];
   aryuWW[191]=MMH*aryuWW[191];
   aryuWW[186]= - 79./36.*aryuWW[6] + 1211./216. + aryuWW[186];
   aryuWW[186]=aryuWW[1]*aryuWW[186];
   aryuWW[192]= - 79./12.*aryuWW[6] + 1211./72. - 6*aryuWW[7];
   aryuWW[192]=aryuWW[40]*aryuWW[192];
   aryuWW[186]=aryuWW[186] + aryuWW[192];
   aryuWW[186]=aryuWW[20]*aryuWW[186];
   aryuWW[115]=aryuWW[115] + aryuWW[117] + aryuWW[132] + aryuWW[191] + 
   1./2.*aryuWW[135] + aryuWW[186];
   aryuWW[115]=aryuWW[165]*aryuWW[115];
   aryuWW[117]= - 7*aryuWW[9];
   aryuWW[132]= - 1./6.*aryuWW[21];
   aryuWW[135]= - 5./18.*aryuWW[23] + aryuWW[132] + 19./9.*aryuWW[22]
    + aryuWW[117] - 17./18.*aryuWW[8];
   aryuWW[135]=aryuWW[1]*aryuWW[135];
   aryuWW[186]= - 21*aryuWW[9];
   aryuWW[191]= - 5./6.*aryuWW[23] + aryuWW[179] + 19./3.*aryuWW[22] + 
   aryuWW[186] - 17./6.*aryuWW[8];
   aryuWW[191]=aryuWW[40]*aryuWW[191];
   aryuWW[192]=19./6. - aryuWW[7];
   aryuWW[192]=1./6.*aryuWW[192] - aryuWW[6];
   aryuWW[192]=aryuWW[1]*aryuWW[192];
   aryuWW[194]=107./54. - aryuWW[7];
   aryuWW[194]=1./2.*aryuWW[194] - 11./9.*aryuWW[6];
   aryuWW[194]=aryuWW[40]*aryuWW[194];
   aryuWW[192]=aryuWW[192] + aryuWW[194];
   aryuWW[192]=aryuWW[18]*aryuWW[192];
   aryuWW[135]=aryuWW[192] + aryuWW[135] + aryuWW[191];
   aryuWW[191]= - aryuWW[16] + aryuWW[190] + aryuWW[28] - 103./81. - 
   aryuWW[31];
   aryuWW[194]=11./27.*aryuWW[6];
   aryuWW[191]=1./2.*aryuWW[191] + aryuWW[194];
   aryuWW[191]=aryuWW[40]*aryuWW[191];
   aryuWW[190]= - aryuWW[16] + aryuWW[190] + aryuWW[28] - 5./3. - 
   aryuWW[31];
   aryuWW[190]=1./2.*aryuWW[190] + aryuWW[6];
   aryuWW[190]=aryuWW[1]*aryuWW[190];
   aryuWW[190]=1./3.*aryuWW[190] + aryuWW[191];
   aryuWW[191]=aryuWW[194] - 38./81. + aryuWW[193];
   aryuWW[191]=aryuWW[40]*aryuWW[191];
   aryuWW[193]=aryuWW[6] - 2./3. + aryuWW[193];
   aryuWW[193]=aryuWW[1]*aryuWW[193];
   aryuWW[191]=1./3.*aryuWW[193] + aryuWW[191];
   aryuWW[191]=aryuWW[11]*aryuWW[191];
   aryuWW[190]=1./2.*aryuWW[190] + aryuWW[191];
   aryuWW[190]=MMH*aryuWW[190];
   aryuWW[191]=49*aryuWW[7];
   aryuWW[193]= - 115./6. + aryuWW[191];
   aryuWW[193]=1./12.*aryuWW[193] + aryuWW[330];
   aryuWW[193]=aryuWW[1]*aryuWW[193];
   aryuWW[191]= - 185./18. + aryuWW[191];
   aryuWW[191]=1./4.*aryuWW[191] + 25./3.*aryuWW[6];
   aryuWW[191]=aryuWW[40]*aryuWW[191];
   aryuWW[191]=aryuWW[193] + aryuWW[191];
   aryuWW[191]=aryuWW[19]*aryuWW[191];
   aryuWW[193]=17*aryuWW[7];
   aryuWW[194]=1471./24. + aryuWW[193];
   aryuWW[194]=1./3.*aryuWW[194] + aryuWW[330];
   aryuWW[194]=aryuWW[1]*aryuWW[194];
   aryuWW[193]=83./9.*aryuWW[6] + 13655./216. + aryuWW[193];
   aryuWW[193]=aryuWW[40]*aryuWW[193];
   aryuWW[193]=aryuWW[194] + aryuWW[193];
   aryuWW[193]=aryuWW[20]*aryuWW[193];
   aryuWW[115]=aryuWW[115] + aryuWW[125] + aryuWW[130] + aryuWW[190] + 
   1./3.*aryuWW[193] + 1./2.*aryuWW[135] + 1./3.*aryuWW[191];
   aryuWW[115]=aryuWW[165]*aryuWW[115];
   aryuWW[125]= - 1081./27. + aryuWW[234];
   aryuWW[125]=1./18.*aryuWW[13] + aryuWW[133] + 1./4.*aryuWW[125] + 
   aryuWW[172];
   aryuWW[125]=aryuWW[306] + aryuWW[136] + 1./2.*aryuWW[125] + 
   aryuWW[199];
   aryuWW[125]=aryuWW[11]*aryuWW[125];
   aryuWW[130]=341./576. + aryuWW[189];
   aryuWW[125]=aryuWW[275] + aryuWW[125] + aryuWW[169] + aryuWW[273] + 
   1./72.*aryuWW[13] + aryuWW[266] + aryuWW[256] + aryuWW[290] + 
   aryuWW[254] + aryuWW[286] + aryuWW[250] + aryuWW[243] + aryuWW[242]
    + aryuWW[285] + aryuWW[237] + aryuWW[228] + aryuWW[218] + 
   aryuWW[217] + aryuWW[216] + aryuWW[280] + aryuWW[213] + 1./3.*
   aryuWW[130] + aryuWW[204];
   aryuWW[125]=MMH*aryuWW[125];
   aryuWW[130]=7487./18. + aryuWW[203];
   aryuWW[130]=aryuWW[303] - 11./6.*aryuWW[11] + aryuWW[301] + 
   aryuWW[300] - 113./72.*aryuWW[12] + aryuWW[295] - 307./72.*
   aryuWW[13] + aryuWW[293] + 1./8.*aryuWW[130] + aryuWW[304];
   aryuWW[130]=aryuWW[20]*aryuWW[130];
   aryuWW[135]=aryuWW[305] - aryuWW[10] + aryuWW[225] + aryuWW[173] + 
   aryuWW[48] - 277./9. + aryuWW[316];
   aryuWW[135]=1./2.*aryuWW[135] + aryuWW[309];
   aryuWW[135]=aryuWW[18]*aryuWW[135];
   aryuWW[136]=65./6.*aryuWW[13] + aryuWW[277] + aryuWW[232] - 7691./27.
    + aryuWW[270];
   aryuWW[136]=aryuWW[323] + aryuWW[320] - 305./16.*aryuWW[12] + 1./4.*
   aryuWW[136] + aryuWW[313];
   aryuWW[136]=aryuWW[19]*aryuWW[136];
   aryuWW[169]=aryuWW[198] + 9*aryuWW[20];
   aryuWW[169]=aryuWW[20]*aryuWW[169];
   aryuWW[162]=aryuWW[162] + aryuWW[169];
   aryuWW[162]=aryuWW[3]*aryuWW[162];
   aryuWW[169]= - aryuWW[18] + 7./2.*aryuWW[19];
   aryuWW[169]=aryuWW[11]*aryuWW[169];
   aryuWW[125]=aryuWW[125] + aryuWW[162] + aryuWW[130] + 1./6.*
   aryuWW[169] + aryuWW[136] + aryuWW[135] - 559./144.*aryuWW[23] + 
   aryuWW[340] - 2315./48.*aryuWW[22] + aryuWW[339] + 365./24.*
   aryuWW[107] + aryuWW[338] + 467./24.*aryuWW[80] + aryuWW[337] + 
   aryuWW[336] - 398./9.*aryuWW[108] + aryuWW[335];
   aryuWW[125]=aryuWW[76]*aryuWW[125];
   aryuWW[130]=1./48. + aryuWW[94];
   aryuWW[130]= - 17./8.*aryuWW[92] - 17./8.*aryuWW[93] + aryuWW[150]
    + 1./2.*aryuWW[130] + aryuWW[99];
   aryuWW[135]=aryuWW[12] - 7./6. - 13*aryuWW[10];
   aryuWW[135]=aryuWW[12]*aryuWW[135];
   aryuWW[136]= - 29./72.*aryuWW[12] + 1./2. + aryuWW[46];
   aryuWW[136]=aryuWW[11]*aryuWW[136];
   aryuWW[150]= - 1./2.*aryuWW[109] - aryuWW[85];
   aryuWW[162]=aryuWW[150] - 7./3.*aryuWW[82];
   aryuWW[162]=MMH*aryuWW[162];
   aryuWW[130]=1./4.*aryuWW[162] + aryuWW[136] + 1./12.*aryuWW[135] + 5.
   /8.*aryuWW[10] + aryuWW[321] + 23./24.*aryuWW[16] + 1./3.*aryuWW[17]
    + 16./9.*aryuWW[100] - 1./12.*aryuWW[96] - 155./48.*aryuWW[103] + 1.
   /3.*aryuWW[130] - 1./8.*aryuWW[98];
   aryuWW[130]=MMH*aryuWW[130];
   aryuWW[135]=9./2.*aryuWW[13] + aryuWW[277] + aryuWW[232] + 4513./27.
    + aryuWW[270];
   aryuWW[135]=59./4.*aryuWW[12] + 1./2.*aryuWW[135] + aryuWW[182];
   aryuWW[136]=aryuWW[18]*aryuWW[113];
   aryuWW[162]=2./3.*aryuWW[113] - 3./8.*aryuWW[112];
   aryuWW[162]=aryuWW[19]*aryuWW[162];
   aryuWW[135]=aryuWW[162] + 1./4.*aryuWW[135] + 1./3.*aryuWW[136];
   aryuWW[135]=aryuWW[19]*aryuWW[135];
   aryuWW[136]=1./3.*aryuWW[297] - 653./36.*aryuWW[12] - 15*aryuWW[10]
    - 31./12.*aryuWW[13] + 5227./108. + aryuWW[324];
   aryuWW[162]=2*aryuWW[113] - 3./4.*aryuWW[112];
   aryuWW[162]=aryuWW[19]*aryuWW[162];
   aryuWW[136]=1./12.*aryuWW[177] - 1./12.*aryuWW[11] + 1./4.*
   aryuWW[136] + aryuWW[162];
   aryuWW[136]=aryuWW[20]*aryuWW[136];
   aryuWW[162]=aryuWW[195] - aryuWW[20];
   aryuWW[162]=aryuWW[20]*aryuWW[162];
   aryuWW[156]=3./2.*aryuWW[156] + aryuWW[162];
   aryuWW[156]=aryuWW[3]*aryuWW[156];
   aryuWW[162]= - 2*aryuWW[108];
   aryuWW[169]=aryuWW[162] + 7./8.*aryuWW[106];
   aryuWW[169]=2051./24.*aryuWW[23] - 73./24.*aryuWW[21] + 815./16.*
   aryuWW[22] - 25./4.*aryuWW[78] - 151./8.*aryuWW[107] - 1./4.*
   aryuWW[79] + 17*aryuWW[169] - 1031./48.*aryuWW[80];
   aryuWW[177]= - aryuWW[12] - 11./3. + aryuWW[157];
   aryuWW[177]=aryuWW[18]*aryuWW[177];
   aryuWW[130]=aryuWW[130] + 1./4.*aryuWW[156] + aryuWW[136] + 5./12.*
   aryuWW[325] + aryuWW[135] + 1./3.*aryuWW[169] + aryuWW[177];
   aryuWW[130]=aryuWW[76]*aryuWW[130];
   aryuWW[135]=3*aryuWW[98];
   aryuWW[136]=1./3.*aryuWW[103] + aryuWW[135] - aryuWW[92] + 61./36.
    - aryuWW[93];
   aryuWW[136]=aryuWW[314] - aryuWW[17] + aryuWW[100] + 1./2.*
   aryuWW[136] - 7./3.*aryuWW[96];
   aryuWW[156]= - 1./2. - aryuWW[10];
   aryuWW[156]=aryuWW[12]*aryuWW[156];
   aryuWW[169]=3*aryuWW[46];
   aryuWW[177]= - 1./6.*aryuWW[12] + 1 + aryuWW[169];
   aryuWW[177]=aryuWW[11]*aryuWW[177];
   aryuWW[189]=1./2.*aryuWW[109] + aryuWW[85];
   aryuWW[191]=aryuWW[189] + aryuWW[82];
   aryuWW[191]=MMH*aryuWW[191];
   aryuWW[136]=1./6.*aryuWW[191] + 1./2.*aryuWW[177] + 1./2.*
   aryuWW[136] + 1./3.*aryuWW[156];
   aryuWW[136]=MMH*aryuWW[136];
   aryuWW[156]= - 1./12.*aryuWW[12];
   aryuWW[177]=aryuWW[156] - 5./9. + aryuWW[173];
   aryuWW[177]=aryuWW[18]*aryuWW[177];
   aryuWW[120]=197./3. + aryuWW[120];
   aryuWW[120]=25./12.*aryuWW[12] + 1./6.*aryuWW[120] + aryuWW[10];
   aryuWW[120]=1./2.*aryuWW[120] + 1./3.*aryuWW[315];
   aryuWW[191]=aryuWW[19]*aryuWW[113];
   aryuWW[120]=1./2.*aryuWW[120] + 2./3.*aryuWW[191];
   aryuWW[120]=aryuWW[19]*aryuWW[120];
   aryuWW[193]= - 1./6.*aryuWW[13] - 7./18. + aryuWW[324];
   aryuWW[194]= - 11./3.*aryuWW[12];
   aryuWW[193]=aryuWW[194] + 1./2.*aryuWW[193] - aryuWW[10];
   aryuWW[191]=1./4.*aryuWW[193] + aryuWW[191];
   aryuWW[191]=aryuWW[20]*aryuWW[191];
   aryuWW[162]=aryuWW[162] + 41./48.*aryuWW[106];
   aryuWW[120]=1./4.*aryuWW[136] + aryuWW[191] + 1./4.*aryuWW[325] + 
   aryuWW[120] + 1./4.*aryuWW[177] + 355./144.*aryuWW[23] + 43./144.*
   aryuWW[21] + 121./48.*aryuWW[22] - 11./24.*aryuWW[78] - 7./12.*
   aryuWW[107] + 1./3.*aryuWW[162] - 3./4.*aryuWW[80];
   aryuWW[120]=aryuWW[76]*aryuWW[120];
   aryuWW[136]=7./3. + aryuWW[6];
   aryuWW[136]=aryuWW[1]*aryuWW[136];
   aryuWW[162]=29./3. + 11*aryuWW[6];
   aryuWW[162]=aryuWW[40]*aryuWW[162];
   aryuWW[136]=aryuWW[136] + 1./9.*aryuWW[162];
   aryuWW[136]=aryuWW[19]*aryuWW[136];
   aryuWW[162]=7./6. - aryuWW[6];
   aryuWW[177]=aryuWW[1]*aryuWW[162];
   aryuWW[162]=aryuWW[40]*aryuWW[162];
   aryuWW[162]=aryuWW[177] + 11./9.*aryuWW[162];
   aryuWW[162]=aryuWW[20]*aryuWW[162];
   aryuWW[177]= - aryuWW[23] - 3*aryuWW[8] + aryuWW[22];
   aryuWW[177]=aryuWW[1]*aryuWW[177];
   aryuWW[191]=11*aryuWW[22];
   aryuWW[193]= - 11*aryuWW[23] - 17*aryuWW[8] + aryuWW[191];
   aryuWW[193]=aryuWW[40]*aryuWW[193];
   aryuWW[136]=aryuWW[162] + aryuWW[136] + aryuWW[177] + 1./9.*
   aryuWW[193];
   aryuWW[162]=5./2.*aryuWW[22] + aryuWW[227] - aryuWW[107];
   aryuWW[162]=1./2.*aryuWW[162] + aryuWW[23];
   aryuWW[177]=1./6.*aryuWW[12];
   aryuWW[193]=1 + aryuWW[177];
   aryuWW[193]=aryuWW[19]*aryuWW[193];
   aryuWW[196]= - 1./3.*aryuWW[12];
   aryuWW[198]= - 1 + aryuWW[196];
   aryuWW[198]=aryuWW[20]*aryuWW[198];
   aryuWW[203]=1 + aryuWW[103];
   aryuWW[203]=MMH*aryuWW[203];
   aryuWW[162]=1./12.*aryuWW[203] + 1./4.*aryuWW[198] + 1./3.*
   aryuWW[162] + 1./2.*aryuWW[193];
   aryuWW[162]=aryuWW[255]*aryuWW[76]*aryuWW[162];
   aryuWW[120]=1./4.*aryuWW[162] + 1./4.*aryuWW[136] + aryuWW[120];
   aryuWW[120]=aryuWW[255]*aryuWW[120];
   aryuWW[136]= - 5./3.*aryuWW[23] - 4./3.*aryuWW[8] + aryuWW[22];
   aryuWW[136]=aryuWW[1]*aryuWW[136];
   aryuWW[162]= - 55*aryuWW[23] - 40*aryuWW[8] + 37*aryuWW[22];
   aryuWW[162]=aryuWW[40]*aryuWW[162];
   aryuWW[136]=aryuWW[136] + 1./27.*aryuWW[162];
   aryuWW[162]=67./12. - aryuWW[7];
   aryuWW[162]=1./3.*aryuWW[162] + 9./4.*aryuWW[6];
   aryuWW[162]=aryuWW[1]*aryuWW[162];
   aryuWW[193]=401./36. - aryuWW[7];
   aryuWW[193]=1./9.*aryuWW[193] + 11./4.*aryuWW[6];
   aryuWW[193]=aryuWW[40]*aryuWW[193];
   aryuWW[162]=aryuWW[162] + aryuWW[193];
   aryuWW[162]=aryuWW[19]*aryuWW[162];
   aryuWW[193]=59./9. - 13./2.*aryuWW[6];
   aryuWW[193]=aryuWW[1]*aryuWW[193];
   aryuWW[198]=641./9. - 143./2.*aryuWW[6];
   aryuWW[198]=aryuWW[40]*aryuWW[198];
   aryuWW[193]=aryuWW[193] + 1./9.*aryuWW[198];
   aryuWW[193]=aryuWW[20]*aryuWW[193];
   aryuWW[120]=aryuWW[120] + aryuWW[130] + 1./2.*aryuWW[193] + 2*
   aryuWW[136] + aryuWW[162];
   aryuWW[120]=aryuWW[255]*aryuWW[120];
   aryuWW[117]=91./18.*aryuWW[23] + aryuWW[132] + 13./9.*aryuWW[22] + 
   aryuWW[117] - 65./18.*aryuWW[8];
   aryuWW[117]=aryuWW[1]*aryuWW[117];
   aryuWW[130]=145./18.*aryuWW[23] + aryuWW[179] + 31./9.*aryuWW[22] + 
   aryuWW[186] - 131./18.*aryuWW[8];
   aryuWW[130]=aryuWW[40]*aryuWW[130];
   aryuWW[117]=aryuWW[192] + aryuWW[117] + aryuWW[130];
   aryuWW[130]= - 139./6. - 11*aryuWW[7];
   aryuWW[130]=1./12.*aryuWW[130] + aryuWW[6];
   aryuWW[130]=aryuWW[1]*aryuWW[130];
   aryuWW[132]=5*aryuWW[7];
   aryuWW[136]= - 11./6. + aryuWW[132];
   aryuWW[136]=1./4.*aryuWW[136] + 5./3.*aryuWW[6];
   aryuWW[136]=aryuWW[40]*aryuWW[136];
   aryuWW[130]=aryuWW[130] + aryuWW[136];
   aryuWW[130]=aryuWW[19]*aryuWW[130];
   aryuWW[136]= - 929./24. + aryuWW[132];
   aryuWW[136]=1./9.*aryuWW[136] - aryuWW[6];
   aryuWW[136]=aryuWW[1]*aryuWW[136];
   aryuWW[132]= - 37./9.*aryuWW[6] + 2039./216. + aryuWW[132];
   aryuWW[132]=aryuWW[40]*aryuWW[132];
   aryuWW[132]=aryuWW[136] + 1./3.*aryuWW[132];
   aryuWW[132]=aryuWW[20]*aryuWW[132];
   aryuWW[117]=aryuWW[120] + aryuWW[125] + aryuWW[190] + aryuWW[132] + 
   1./2.*aryuWW[117] + 1./3.*aryuWW[130];
   aryuWW[117]=aryuWW[255]*aryuWW[117];
   aryuWW[120]= - 4*aryuWW[108] + aryuWW[80];
   aryuWW[120]=4*aryuWW[23] + 2*aryuWW[120] - 17./3.*aryuWW[22];
   aryuWW[125]=4*aryuWW[12] + 1 + aryuWW[13];
   aryuWW[125]=aryuWW[19]*aryuWW[125];
   aryuWW[130]=5*aryuWW[12] + 8 + aryuWW[13];
   aryuWW[130]=aryuWW[20]*aryuWW[130];
   aryuWW[120]=2*aryuWW[130] + 2*aryuWW[120] + aryuWW[125];
   aryuWW[120]=aryuWW[76]*aryuWW[120];
   aryuWW[114]=aryuWW[114] + aryuWW[115] + aryuWW[171] + 4*aryuWW[120]
    + aryuWW[117];
   aryuWW[114]=aryuWW[4]*aryuWW[114];
   aryuWW[115]=2*aryuWW[36];
   aryuWW[117]=aryuWW[6] + 1 + aryuWW[115];
   aryuWW[117]=2*aryuWW[117] + aryuWW[268];
   aryuWW[117]=aryuWW[20]*aryuWW[117];
   aryuWW[120]= - aryuWW[42] - aryuWW[8];
   aryuWW[120]=aryuWW[22] + 3*aryuWW[51] + 1./2.*aryuWW[120] - 2*
   aryuWW[43];
   aryuWW[125]=6*aryuWW[23];
   aryuWW[130]=17./4.*aryuWW[6] - 12 + 25./4.*aryuWW[36];
   aryuWW[130]=MMZ*aryuWW[130];
   aryuWW[132]=7./4.*aryuWW[6] + 1 + 11./4.*aryuWW[36];
   aryuWW[132]=aryuWW[19]*aryuWW[132];
   aryuWW[136]=3*MMZ + aryuWW[289];
   aryuWW[136]=aryuWW[37]*aryuWW[136];
   aryuWW[162]=6*aryuWW[49] + 37./2. + aryuWW[123];
   aryuWW[162]=aryuWW[38]*aryuWW[162];
   aryuWW[117]=aryuWW[162] + aryuWW[117] + aryuWW[136] + aryuWW[132] + 
   aryuWW[185] + aryuWW[130] + aryuWW[125] + 3*aryuWW[120] + 
   aryuWW[244];
   aryuWW[117]=aryuWW[3]*aryuWW[117];
   aryuWW[120]= - 2833./64. + aryuWW[123];
   aryuWW[120]=1./2.*aryuWW[120] + aryuWW[219];
   aryuWW[130]= - 22*aryuWW[13];
   aryuWW[132]= - 81./512.*aryuWW[111] + aryuWW[112];
   aryuWW[132]=MMZ*aryuWW[132];
   aryuWW[136]= - aryuWW[35] + 81./64.*aryuWW[111];
   aryuWW[136]=aryuWW[19]*aryuWW[136];
   aryuWW[162]= - MMZ*aryuWW[112];
   aryuWW[171]= - 1 + aryuWW[162];
   aryuWW[171]=aryuWW[11]*aryuWW[171];
   aryuWW[120]=5./4.*aryuWW[37] + aryuWW[171] + 1./4.*aryuWW[136] + 1./
   2.*aryuWW[297] + aryuWW[132] + 2*aryuWW[144] + 2./3.*aryuWW[159] + 
   33./4.*aryuWW[12] + 1./4.*aryuWW[120] + aryuWW[130];
   aryuWW[120]=aryuWW[37]*aryuWW[120];
   aryuWW[132]=7./4. + aryuWW[69];
   aryuWW[132]=1./2.*aryuWW[132] + aryuWW[71];
   aryuWW[136]=pow(aryuWW[110],2);
   aryuWW[132]=MMZ*aryuWW[136]*aryuWW[132];
   aryuWW[142]=3*aryuWW[142] - aryuWW[59];
   aryuWW[144]= - 1./2.*aryuWW[66] + 3./8. - aryuWW[67];
   aryuWW[144]=3./4.*aryuWW[50] + 1./2.*aryuWW[144] - aryuWW[33];
   aryuWW[144]=aryuWW[35]*aryuWW[144];
   aryuWW[126]=aryuWW[74] + aryuWW[126];
   aryuWW[126]=aryuWW[110]*aryuWW[126];
   aryuWW[159]=3./2.*aryuWW[126] + 1./4.*aryuWW[71] - 103./16. - 
   aryuWW[69];
   aryuWW[159]=aryuWW[110]*aryuWW[159];
   aryuWW[171]=205 + 9*aryuWW[66];
   aryuWW[171]=11*aryuWW[17] - 9./8.*aryuWW[50] - 9./8.*aryuWW[71] + 9./
   4.*aryuWW[69] + 1./8.*aryuWW[171] - 11*aryuWW[63];
   aryuWW[171]=aryuWW[111]*aryuWW[171];
   aryuWW[185]=aryuWW[12]*aryuWW[111];
   aryuWW[186]= - aryuWW[16] + 6*aryuWW[68] + 15./4. + aryuWW[60];
   aryuWW[186]=aryuWW[112]*aryuWW[186];
   aryuWW[132]=3./4.*aryuWW[132] + aryuWW[186] + 99./64.*aryuWW[185] + 
   9./64.*aryuWW[171] + 1./2.*aryuWW[159] + 2*aryuWW[142] + aryuWW[144]
   ;
   aryuWW[132]=MMZ*aryuWW[132];
   aryuWW[142]=7./2.*aryuWW[16];
   aryuWW[144]=aryuWW[131] + aryuWW[142] + aryuWW[263] + aryuWW[70] + 
   aryuWW[287] + 1./4. + aryuWW[230];
   aryuWW[144]=aryuWW[112]*aryuWW[144];
   aryuWW[123]=47./2. + aryuWW[123];
   aryuWW[159]= - 1 + aryuWW[48];
   aryuWW[159]=2*aryuWW[159] + aryuWW[46];
   aryuWW[159]=aryuWW[3]*aryuWW[38]*aryuWW[159];
   aryuWW[123]=9*aryuWW[159] + aryuWW[178] + 23./4.*aryuWW[36] - 9./2.*
   aryuWW[39] + 1./2.*aryuWW[123] + aryuWW[219];
   aryuWW[123]=aryuWW[3]*aryuWW[123];
   aryuWW[159]=aryuWW[56] + 1./2.*aryuWW[54];
   aryuWW[171]=2*aryuWW[59];
   aryuWW[178]=aryuWW[11]*aryuWW[112]*aryuWW[39];
   aryuWW[186]=MMt*aryuWW[274]*aryuWW[39];
   aryuWW[123]=18*aryuWW[186] + aryuWW[123] + aryuWW[283] + 3*
   aryuWW[178] + aryuWW[144] + 2*aryuWW[55] + 7*aryuWW[159] + 
   aryuWW[171];
   aryuWW[123]=MMt*aryuWW[123];
   aryuWW[144]= - 2*aryuWW[36];
   aryuWW[159]=1 + aryuWW[144];
   aryuWW[178]=aryuWW[7]*aryuWW[36];
   aryuWW[190]= - 4./3. + aryuWW[7];
   aryuWW[190]=aryuWW[6]*aryuWW[190];
   aryuWW[159]=aryuWW[190] + 1./3.*aryuWW[159] + 2*aryuWW[178];
   aryuWW[178]=aryuWW[1]*aryuWW[159];
   aryuWW[159]=aryuWW[40]*aryuWW[159];
   aryuWW[190]= - aryuWW[12]*aryuWW[111];
   aryuWW[192]=33./4.*aryuWW[190] + 561./32.*aryuWW[111] - aryuWW[35]
    + 61*aryuWW[110];
   aryuWW[193]= - MMZ*aryuWW[136];
   aryuWW[198]= - 81./32.*aryuWW[111] + aryuWW[35] - 9./2.*aryuWW[110];
   aryuWW[198]=aryuWW[37]*aryuWW[198];
   aryuWW[192]=1./4.*aryuWW[198] + 9./8.*aryuWW[193] + 3./4.*
   aryuWW[192] - 2*aryuWW[112];
   aryuWW[192]=aryuWW[38]*aryuWW[192];
   aryuWW[198]= - 4*aryuWW[36];
   aryuWW[203]= - 2*aryuWW[6] + 1 + aryuWW[198];
   aryuWW[204]=aryuWW[115] + aryuWW[6];
   aryuWW[204]=2./3.*aryuWW[204] - aryuWW[37];
   aryuWW[204]=aryuWW[37]*aryuWW[204];
   aryuWW[203]=1./9.*aryuWW[203] + aryuWW[204];
   aryuWW[203]=aryuWW[34]*aryuWW[203];
   aryuWW[204]=34735./3456. + 33*aryuWW[61];
   aryuWW[205]=aryuWW[74] + aryuWW[319];
   aryuWW[206]=aryuWW[205] - aryuWW[43];
   aryuWW[206]=aryuWW[35]*aryuWW[206];
   aryuWW[211]=aryuWW[51]*aryuWW[35];
   aryuWW[213]=57./2.*aryuWW[51] + 9./8.*aryuWW[43] + 1./8.*aryuWW[8]
    - 1./4.*aryuWW[74] - 15*aryuWW[72];
   aryuWW[213]=aryuWW[110]*aryuWW[213];
   aryuWW[214]= - aryuWW[35] - 3./2.*aryuWW[110];
   aryuWW[214]=aryuWW[22]*aryuWW[214];
   aryuWW[216]= - aryuWW[35] + 111./4.*aryuWW[110];
   aryuWW[216]=aryuWW[23]*aryuWW[216];
   aryuWW[217]=11*aryuWW[23] + 33./16.*aryuWW[22] + 155./8.*aryuWW[51]
    - 3./8.*aryuWW[42] - 11*aryuWW[72];
   aryuWW[217]=aryuWW[111]*aryuWW[217];
   aryuWW[218]=22./3.*aryuWW[13];
   aryuWW[219]= - 293./24. - 9*aryuWW[36];
   aryuWW[219]=aryuWW[12]*aryuWW[219];
   aryuWW[222]= - 11./2.*aryuWW[12];
   aryuWW[224]= - 31./27. + aryuWW[222];
   aryuWW[224]=aryuWW[6]*aryuWW[224];
   aryuWW[227]=aryuWW[23] + aryuWW[179] + 2*aryuWW[73] - aryuWW[43];
   aryuWW[227]=aryuWW[112]*aryuWW[227];
   aryuWW[228]=3*aryuWW[193] - 63./64.*aryuWW[111] - aryuWW[35] + 
   aryuWW[110];
   aryuWW[228]=aryuWW[19]*aryuWW[228];
   aryuWW[230]=2./3.*aryuWW[6] + 1 + 4./3.*aryuWW[36];
   aryuWW[230]=aryuWW[11]*aryuWW[230];
   aryuWW[233]= - aryuWW[35] + 9./4.*aryuWW[110];
   aryuWW[233]=1./2.*aryuWW[233] + aryuWW[112];
   aryuWW[233]=aryuWW[37]*aryuWW[233];
   aryuWW[233]=aryuWW[233] - 33./32.*aryuWW[111] + aryuWW[35] - 39./4.*
   aryuWW[110];
   aryuWW[233]=aryuWW[20]*aryuWW[233];
   aryuWW[117]=aryuWW[203] + aryuWW[123] + aryuWW[117] + aryuWW[192] + 
   aryuWW[233] + aryuWW[120] + aryuWW[230] + 1./4.*aryuWW[228] + 5./2.*
   aryuWW[297] + aryuWW[132] + 2./3.*aryuWW[159] + 2./9.*aryuWW[178] + 
   aryuWW[227] + 5./4.*aryuWW[224] + 11./8.*aryuWW[219] - 391./108.*
   aryuWW[36] + aryuWW[218] + 3./16.*aryuWW[217] + 1./2.*aryuWW[216] + 
   aryuWW[16] + 1./4.*aryuWW[214] + aryuWW[213] - 957./64.*aryuWW[17]
    + aryuWW[211] + 1./2.*aryuWW[206] + aryuWW[124] - 3./4.*aryuWW[47]
    + 2*aryuWW[68] - 1681./512.*aryuWW[50] - aryuWW[60] - 145./512.*
   aryuWW[71] - 207./256.*aryuWW[69] + 33./8.*aryuWW[30] - 21./8.*
   aryuWW[33] + 165./64.*aryuWW[63] + 657./512.*aryuWW[66] + 1./4.*
   aryuWW[204] + 2*aryuWW[67];
   aryuWW[117]=aryuWW[34]*aryuWW[117];
   aryuWW[120]= - 33./2.*aryuWW[12];
   aryuWW[123]=aryuWW[120] + aryuWW[163] - 825./4.*aryuWW[13] + 
   aryuWW[173] + 133./4. + aryuWW[172];
   aryuWW[123]=MMZ*aryuWW[123];
   aryuWW[124]=3*aryuWW[10];
   aryuWW[132]=aryuWW[164] + aryuWW[124] - 363./4.*aryuWW[13] + 
   aryuWW[173] - 631./12. + aryuWW[172];
   aryuWW[132]=aryuWW[19]*aryuWW[132];
   aryuWW[151]=aryuWW[151] + 1 - aryuWW[48];
   aryuWW[159]=pow(MMZ,2);
   aryuWW[178]=aryuWW[159]*aryuWW[151];
   aryuWW[151]=MMZ*aryuWW[151];
   aryuWW[192]=aryuWW[19]*aryuWW[151];
   aryuWW[151]=aryuWW[20]*aryuWW[151];
   aryuWW[151]=aryuWW[151] + aryuWW[178] + 1./2.*aryuWW[192];
   aryuWW[151]=aryuWW[3]*aryuWW[151];
   aryuWW[178]= - 99./2.*aryuWW[12] + aryuWW[173] + 115./2. + 
   aryuWW[172];
   aryuWW[178]=aryuWW[18]*aryuWW[178];
   aryuWW[192]= - 1./2.*aryuWW[78];
   aryuWW[203]=aryuWW[192] + 11*aryuWW[80] - aryuWW[79];
   aryuWW[203]=1./2.*aryuWW[203] - 5*aryuWW[22];
   aryuWW[204]=aryuWW[157] + 7./2. + aryuWW[48];
   aryuWW[204]=MMZ*aryuWW[204];
   aryuWW[204]=3./4.*aryuWW[19] + 3./2.*aryuWW[204] + aryuWW[18];
   aryuWW[204]=aryuWW[11]*aryuWW[204];
   aryuWW[161]=9./2.*aryuWW[11] + aryuWW[161] - 66*aryuWW[13] + 9./4.*
   aryuWW[46] - 226./3. + 9./2.*aryuWW[48];
   aryuWW[161]=aryuWW[20]*aryuWW[161];
   aryuWW[123]=9*aryuWW[151] + aryuWW[161] + aryuWW[204] + 1./2.*
   aryuWW[132] + 1./4.*aryuWW[178] + 1./2.*aryuWW[123] - 45*aryuWW[23]
    + 9./2.*aryuWW[203] - 8*aryuWW[21];
   aryuWW[123]=aryuWW[3]*aryuWW[123];
   aryuWW[132]=11./4.*aryuWW[12];
   aryuWW[151]=aryuWW[132] + aryuWW[231] + aryuWW[317] + 1./2. + 3*
   aryuWW[44];
   aryuWW[151]=aryuWW[112]*aryuWW[151];
   aryuWW[151]=1./4.*aryuWW[113] + aryuWW[151];
   aryuWW[151]=MMZ*aryuWW[151];
   aryuWW[161]= - 145./3. + aryuWW[234];
   aryuWW[178]= - aryuWW[19]*aryuWW[113];
   aryuWW[162]= - 1./2. + aryuWW[162];
   aryuWW[162]=aryuWW[11]*aryuWW[162];
   aryuWW[151]=aryuWW[162] + 3./4.*aryuWW[178] + 3*aryuWW[151] + 11*
   aryuWW[12] + aryuWW[130] - aryuWW[46] + 1./4.*aryuWW[161] + 
   aryuWW[304];
   aryuWW[151]=aryuWW[11]*aryuWW[151];
   aryuWW[161]= - 11*aryuWW[98];
   aryuWW[162]=11*aryuWW[103];
   aryuWW[203]= - aryuWW[96] + aryuWW[162] + 5./2. + aryuWW[161];
   aryuWW[204]= - 11*aryuWW[100];
   aryuWW[206]=11./2.*aryuWW[17];
   aryuWW[203]=aryuWW[253] + aryuWW[206] + 1./2.*aryuWW[203] + 
   aryuWW[204];
   aryuWW[203]=aryuWW[113]*aryuWW[203];
   aryuWW[155]= - 11./4.*aryuWW[92] - 11./4.*aryuWW[93] - aryuWW[97] + 
   16*aryuWW[101] - 5./16. + aryuWW[155];
   aryuWW[213]= - 3./4.*aryuWW[46];
   aryuWW[155]=aryuWW[213] + 19*aryuWW[16] + aryuWW[116] + 33./4.*
   aryuWW[17] + aryuWW[288] - 66*aryuWW[100] - 3./4.*aryuWW[96] - 33./4.
   *aryuWW[103] + 3*aryuWW[155] + 2*aryuWW[91];
   aryuWW[155]=aryuWW[112]*aryuWW[155];
   aryuWW[214]=aryuWW[12]*aryuWW[113];
   aryuWW[155]=aryuWW[155] + 33./4.*aryuWW[214] + 3./2.*aryuWW[203] - 
   33./2.*aryuWW[82] - 33./2.*aryuWW[85] - 33./4.*aryuWW[109] + 459./16.
   *aryuWW[83] + aryuWW[81];
   aryuWW[155]=MMZ*aryuWW[155];
   aryuWW[148]=aryuWW[154] + aryuWW[148] + 11./3. + aryuWW[149];
   aryuWW[148]=aryuWW[112]*aryuWW[148];
   aryuWW[149]= - aryuWW[12]*aryuWW[113];
   aryuWW[154]=11*aryuWW[149];
   aryuWW[203]=15*aryuWW[113] + aryuWW[154];
   aryuWW[148]=3./8.*aryuWW[203] + aryuWW[148];
   aryuWW[148]=aryuWW[18]*aryuWW[148];
   aryuWW[203]=aryuWW[213] + aryuWW[116] + 179./24. + aryuWW[288];
   aryuWW[203]=aryuWW[112]*aryuWW[203];
   aryuWW[216]= - aryuWW[11]*aryuWW[112];
   aryuWW[203]=3./2.*aryuWW[216] + 7./8.*aryuWW[113] + aryuWW[203];
   aryuWW[203]=aryuWW[20]*aryuWW[203];
   aryuWW[216]= - 435./8.*aryuWW[94] - 2413./81. - 9./2.*aryuWW[95];
   aryuWW[217]= - 1./6.*aryuWW[97];
   aryuWW[219]= - 3./2.*aryuWW[96];
   aryuWW[224]=5./8.*aryuWW[46];
   aryuWW[227]=aryuWW[294] - 21./4.*aryuWW[21] + 31./2.*aryuWW[22] + 7./
   2.*aryuWW[107] - 3*aryuWW[78];
   aryuWW[227]=aryuWW[113]*aryuWW[227];
   aryuWW[228]=6973./9. + 3267./2.*aryuWW[13];
   aryuWW[228]= - 381./2.*aryuWW[12] + 1./2.*aryuWW[228] - 11*
   aryuWW[10];
   aryuWW[228]=aryuWW[12]*aryuWW[228];
   aryuWW[230]= - 47./12.*aryuWW[23] + 181./24.*aryuWW[21] + 1./4.*
   aryuWW[22] - 3./8.*aryuWW[78] + 1./4.*aryuWW[107] - 11./3.*
   aryuWW[79] + 33./4.*aryuWW[80] - 33./2.*aryuWW[106] - 9./2.*
   aryuWW[77] + 11*aryuWW[105];
   aryuWW[230]=aryuWW[112]*aryuWW[230];
   aryuWW[233]=11*aryuWW[214];
   aryuWW[235]= - 13*aryuWW[113] + aryuWW[233];
   aryuWW[235]=1./2.*aryuWW[235] + 11*aryuWW[112];
   aryuWW[235]=aryuWW[19]*aryuWW[235];
   aryuWW[236]=13./2.*aryuWW[82] + 13./2.*aryuWW[85] + 13./4.*
   aryuWW[109] + 3./2.*aryuWW[84] - 2./3.*aryuWW[81];
   aryuWW[236]=MMH*aryuWW[236];
   aryuWW[123]=aryuWW[236] + aryuWW[123] + aryuWW[203] + aryuWW[151] + 
   3./2.*aryuWW[235] + aryuWW[148] + aryuWW[155] + aryuWW[230] + 1./4.*
   aryuWW[228] + 4301./72.*aryuWW[13] + 1./2.*aryuWW[227] + aryuWW[224]
    + 7./6.*aryuWW[16] + 5./4.*aryuWW[48] + 391./8.*aryuWW[17] - 45./8.
   *aryuWW[44] - 33*aryuWW[100] + aryuWW[219] + 11./4.*aryuWW[98] + 11./
   4.*aryuWW[92] + 11./4.*aryuWW[93] + aryuWW[217] - 435./16.*
   aryuWW[99] + 1./2.*aryuWW[216] + 13*aryuWW[101];
   aryuWW[123]=aryuWW[76]*aryuWW[123];
   aryuWW[148]= - 59./9. + 33./2.*aryuWW[13];
   aryuWW[151]= - aryuWW[6] - 1./3. + 2*aryuWW[7];
   aryuWW[155]=1./3.*aryuWW[1]*aryuWW[151];
   aryuWW[151]=aryuWW[40]*aryuWW[151];
   aryuWW[175]=aryuWW[175] + aryuWW[151] + aryuWW[155] + aryuWW[120] + 
   1./2.*aryuWW[148] - aryuWW[10];
   aryuWW[175]=aryuWW[37]*aryuWW[175];
   aryuWW[203]=MMZ*aryuWW[187];
   aryuWW[187]=aryuWW[19]*aryuWW[187];
   aryuWW[187]=aryuWW[203] + 1./2.*aryuWW[187];
   aryuWW[203]=MMZ + aryuWW[245];
   aryuWW[216]=aryuWW[37]*aryuWW[203];
   aryuWW[180]=aryuWW[20]*aryuWW[180];
   aryuWW[180]=aryuWW[180] + 1./2.*aryuWW[187] + aryuWW[216];
   aryuWW[180]=aryuWW[3]*aryuWW[180];
   aryuWW[187]=aryuWW[7]*aryuWW[208];
   aryuWW[208]=1 - aryuWW[7];
   aryuWW[216]=aryuWW[6]*aryuWW[208];
   aryuWW[187]=1./2.*aryuWW[216] + 1./6.*aryuWW[36] + aryuWW[187];
   aryuWW[216]=aryuWW[1]*aryuWW[187];
   aryuWW[187]=aryuWW[40]*aryuWW[187];
   aryuWW[209]=aryuWW[209] + aryuWW[37];
   aryuWW[209]=aryuWW[37]*aryuWW[209];
   aryuWW[181]=1./6.*aryuWW[181] + aryuWW[209];
   aryuWW[181]=aryuWW[34]*aryuWW[181];
   aryuWW[134]=aryuWW[134] + 59./9.*aryuWW[36];
   aryuWW[174]=1 + aryuWW[174];
   aryuWW[174]=aryuWW[12]*aryuWW[174];
   aryuWW[209]=59./9. + aryuWW[164];
   aryuWW[209]=aryuWW[6]*aryuWW[209];
   aryuWW[170]=MMt*aryuWW[3]*aryuWW[170];
   aryuWW[134]=aryuWW[181] + 3./2.*aryuWW[170] + 3*aryuWW[180] + 
   aryuWW[175] + aryuWW[210] + aryuWW[187] + 1./3.*aryuWW[216] + 1./4.*
   aryuWW[209] + 11./4.*aryuWW[174] + 1./4.*aryuWW[134] + aryuWW[199];
   aryuWW[134]=aryuWW[34]*aryuWW[134];
   aryuWW[148]=aryuWW[7]*aryuWW[148];
   aryuWW[148]=aryuWW[160] + aryuWW[148];
   aryuWW[160]=aryuWW[7]*aryuWW[220];
   aryuWW[121]=aryuWW[6]*aryuWW[121];
   aryuWW[160]=aryuWW[160] + aryuWW[121];
   aryuWW[170]=aryuWW[1]*aryuWW[160];
   aryuWW[160]=aryuWW[40]*aryuWW[160];
   aryuWW[174]=1./2. - 3*aryuWW[7];
   aryuWW[174]=aryuWW[12]*aryuWW[174];
   aryuWW[148]=aryuWW[160] + 2./3.*aryuWW[170] + 1./2.*aryuWW[209] + 11.
   /2.*aryuWW[174] + 1./2.*aryuWW[148] + aryuWW[168];
   aryuWW[148]=aryuWW[40]*aryuWW[148];
   aryuWW[119]= - 59./27. + aryuWW[119];
   aryuWW[119]=aryuWW[7]*aryuWW[119];
   aryuWW[119]= - 11./6.*aryuWW[13] + aryuWW[119];
   aryuWW[160]=aryuWW[12]*aryuWW[166];
   aryuWW[158]=59./27. + aryuWW[158];
   aryuWW[158]=aryuWW[6]*aryuWW[158];
   aryuWW[119]=1./9.*aryuWW[170] + 1./2.*aryuWW[158] + 11./2.*
   aryuWW[160] + 1./2.*aryuWW[119] + 1./3.*aryuWW[168];
   aryuWW[119]=aryuWW[1]*aryuWW[119];
   aryuWW[158]=59./3. - 99./2.*aryuWW[13];
   aryuWW[158]=99./4.*aryuWW[12] + 1./2.*aryuWW[158] + aryuWW[124];
   aryuWW[158]=aryuWW[12]*aryuWW[158];
   aryuWW[160]=aryuWW[226] + 1./2.*aryuWW[107] - aryuWW[22];
   aryuWW[160]=aryuWW[113]*aryuWW[160];
   aryuWW[166]=1./2.*aryuWW[23] - 1./2.*aryuWW[107] + aryuWW[22];
   aryuWW[166]=aryuWW[112]*aryuWW[166];
   aryuWW[168]=MMZ*aryuWW[184];
   aryuWW[170]= - aryuWW[113] + aryuWW[112];
   aryuWW[170]=aryuWW[19]*aryuWW[170];
   aryuWW[158]=9./2.*aryuWW[170] + 9./16.*aryuWW[168] + 3*aryuWW[166]
    + 11./2.*aryuWW[158] + 77./9.*aryuWW[10] + 3*aryuWW[160] - 649./12.
   *aryuWW[13];
   aryuWW[160]= - 17 + 99./2.*aryuWW[13];
   aryuWW[166]= - 99./4.*aryuWW[12];
   aryuWW[160]=aryuWW[166] + 1./2.*aryuWW[160] + aryuWW[153];
   aryuWW[160]=aryuWW[19]*aryuWW[160];
   aryuWW[168]=17 + 99*aryuWW[13];
   aryuWW[153]=3*aryuWW[11] + aryuWW[166] + 1./4.*aryuWW[168] + 
   aryuWW[153];
   aryuWW[153]=aryuWW[20]*aryuWW[153];
   aryuWW[166]=MMZ*aryuWW[176];
   aryuWW[168]=aryuWW[11]*aryuWW[203];
   aryuWW[153]=aryuWW[153] + 3*aryuWW[168] + 3*aryuWW[166] + 1./2.*
   aryuWW[160];
   aryuWW[153]=aryuWW[3]*aryuWW[153];
   aryuWW[160]= - 7./9. + 3./2.*aryuWW[13];
   aryuWW[120]=aryuWW[11] + aryuWW[120] + 11./2.*aryuWW[160] - 
   aryuWW[10];
   aryuWW[120]=aryuWW[11]*aryuWW[120];
   aryuWW[160]=aryuWW[20]*aryuWW[184];
   aryuWW[120]=aryuWW[153] + 3./8.*aryuWW[160] + 1./2.*aryuWW[158] + 
   aryuWW[120];
   aryuWW[120]=aryuWW[76]*aryuWW[120];
   aryuWW[153]=aryuWW[202] + 3*aryuWW[183];
   aryuWW[158]=MMZ*aryuWW[153];
   aryuWW[160]=aryuWW[19]*aryuWW[153];
   aryuWW[153]=aryuWW[20]*aryuWW[153];
   aryuWW[153]=aryuWW[153] + aryuWW[158] + 1./2.*aryuWW[160];
   aryuWW[153]=aryuWW[3]*aryuWW[153];
   aryuWW[151]=aryuWW[155] + aryuWW[151];
   aryuWW[151]=aryuWW[11]*aryuWW[151];
   aryuWW[119]=aryuWW[134] + aryuWW[120] + aryuWW[153] + aryuWW[151] + 
   aryuWW[119] + aryuWW[148];
   aryuWW[119]=aryuWW[165]*aryuWW[119];
   aryuWW[120]=aryuWW[6]*aryuWW[220];
   aryuWW[134]=pow(Pi,2);
   aryuWW[148]=pow(aryuWW[7],2);
   aryuWW[120]=2*aryuWW[120] - aryuWW[148] + 1./9. + aryuWW[134];
   aryuWW[151]=aryuWW[1]*aryuWW[120];
   aryuWW[153]= - 5./3.*aryuWW[134] - 37./3.*aryuWW[32] - 3989./648. + 
   11*aryuWW[29];
   aryuWW[155]= - 1./3.*aryuWW[28];
   aryuWW[158]=1./3.*aryuWW[16];
   aryuWW[160]= - 5./9. + 1./4.*aryuWW[7];
   aryuWW[160]=aryuWW[12]*aryuWW[160];
   aryuWW[130]=5./4.*aryuWW[7] - 23./4. + aryuWW[130];
   aryuWW[130]=aryuWW[7]*aryuWW[130];
   aryuWW[166]= - 13./9. + aryuWW[222];
   aryuWW[166]=aryuWW[6]*aryuWW[166];
   aryuWW[168]=aryuWW[23] - aryuWW[9] + aryuWW[179];
   aryuWW[168]=aryuWW[112]*aryuWW[168];
   aryuWW[153]=1./9.*aryuWW[151] + 1./3.*aryuWW[168] + 7./6.*
   aryuWW[166] + 11*aryuWW[160] + 1./3.*aryuWW[130] + 22./9.*aryuWW[13]
    + aryuWW[158] - 11./2.*aryuWW[17] + 11./4.*aryuWW[30] - 19./4.*
   aryuWW[33] - 23./72.*aryuWW[14] + 1./4.*aryuWW[153] + aryuWW[155];
   aryuWW[153]=aryuWW[1]*aryuWW[153];
   aryuWW[120]=aryuWW[40]*aryuWW[120];
   aryuWW[160]= - 5*aryuWW[134] - 37*aryuWW[32] - 3989./216. + 33*
   aryuWW[29];
   aryuWW[170]= - 5./3. + 3./4.*aryuWW[7];
   aryuWW[170]=aryuWW[12]*aryuWW[170];
   aryuWW[120]=aryuWW[120] + 2./3.*aryuWW[151] + aryuWW[168] + 7./2.*
   aryuWW[166] + 11*aryuWW[170] + aryuWW[130] + aryuWW[218] + 
   aryuWW[16] - 33./2.*aryuWW[17] + 33./4.*aryuWW[30] - 57./4.*
   aryuWW[33] - 23./24.*aryuWW[14] + 1./4.*aryuWW[160] - aryuWW[28];
   aryuWW[120]=aryuWW[40]*aryuWW[120];
   aryuWW[130]=aryuWW[22] - 2*aryuWW[9] - aryuWW[8];
   aryuWW[125]=aryuWW[125] + 3*aryuWW[130] + aryuWW[244];
   aryuWW[125]=aryuWW[40]*aryuWW[125];
   aryuWW[151]=7./2.*aryuWW[6];
   aryuWW[160]=aryuWW[151] - 4 + aryuWW[7];
   aryuWW[166]=aryuWW[1]*aryuWW[160];
   aryuWW[160]=aryuWW[40]*aryuWW[160];
   aryuWW[160]=aryuWW[166] + 3*aryuWW[160];
   aryuWW[160]=MMZ*aryuWW[160];
   aryuWW[166]= - 1 + aryuWW[7];
   aryuWW[168]=aryuWW[1]*aryuWW[166];
   aryuWW[166]=aryuWW[40]*aryuWW[166];
   aryuWW[166]=1./3.*aryuWW[168] + aryuWW[166];
   aryuWW[166]=aryuWW[18]*aryuWW[166];
   aryuWW[168]=1./3. + 3./2.*aryuWW[6];
   aryuWW[168]=aryuWW[1]*aryuWW[168];
   aryuWW[170]=1 + 9./2.*aryuWW[6];
   aryuWW[170]=aryuWW[40]*aryuWW[170];
   aryuWW[168]=aryuWW[168] + aryuWW[170];
   aryuWW[168]=aryuWW[19]*aryuWW[168];
   aryuWW[170]=2*aryuWW[6];
   aryuWW[174]=aryuWW[170] + 2./3. + aryuWW[7];
   aryuWW[174]=aryuWW[1]*aryuWW[174];
   aryuWW[175]=6*aryuWW[6] + 2 + 3*aryuWW[7];
   aryuWW[175]=aryuWW[40]*aryuWW[175];
   aryuWW[174]=aryuWW[174] + aryuWW[175];
   aryuWW[174]=aryuWW[20]*aryuWW[174];
   aryuWW[130]=2*aryuWW[23] + aryuWW[130] + 1./6.*aryuWW[21];
   aryuWW[130]=aryuWW[1]*aryuWW[130];
   aryuWW[125]=aryuWW[174] + aryuWW[168] + 1./2.*aryuWW[166] + 
   aryuWW[160] + aryuWW[130] + aryuWW[125];
   aryuWW[125]=aryuWW[3]*aryuWW[125];
   aryuWW[130]=aryuWW[208] + aryuWW[170];
   aryuWW[160]=aryuWW[1]*aryuWW[130];
   aryuWW[130]=aryuWW[40]*aryuWW[130];
   aryuWW[166]= - aryuWW[112]*aryuWW[7];
   aryuWW[168]=aryuWW[1]*aryuWW[166];
   aryuWW[166]=aryuWW[40]*aryuWW[166];
   aryuWW[166]=1./3.*aryuWW[168] + aryuWW[166];
   aryuWW[166]=MMZ*aryuWW[166];
   aryuWW[130]=aryuWW[166] + 1./3.*aryuWW[160] + aryuWW[130];
   aryuWW[130]=aryuWW[11]*aryuWW[130];
   aryuWW[160]= - aryuWW[25] - aryuWW[24];
   aryuWW[166]=aryuWW[160] - 1./3.*aryuWW[27];
   aryuWW[168]=1./3.*aryuWW[7] + aryuWW[212] - aryuWW[14] + 1./3.*
   aryuWW[28] - 1 - 2./3.*aryuWW[31];
   aryuWW[168]=aryuWW[112]*aryuWW[168];
   aryuWW[166]=2*aryuWW[166] + aryuWW[168];
   aryuWW[166]=aryuWW[1]*aryuWW[166];
   aryuWW[160]=3*aryuWW[160] - aryuWW[27];
   aryuWW[168]=aryuWW[7] - aryuWW[16] - 3*aryuWW[14] + aryuWW[28] - 3
    - 2*aryuWW[31];
   aryuWW[168]=aryuWW[112]*aryuWW[168];
   aryuWW[160]=2*aryuWW[160] + aryuWW[168];
   aryuWW[160]=aryuWW[40]*aryuWW[160];
   aryuWW[160]=aryuWW[166] + aryuWW[160];
   aryuWW[160]=MMZ*aryuWW[160];
   aryuWW[166]=aryuWW[112]*aryuWW[208];
   aryuWW[168]=aryuWW[1]*aryuWW[166];
   aryuWW[166]=aryuWW[40]*aryuWW[166];
   aryuWW[166]=1./3.*aryuWW[168] + aryuWW[166];
   aryuWW[166]=aryuWW[18]*aryuWW[166];
   aryuWW[168]=aryuWW[112]*aryuWW[7];
   aryuWW[170]=aryuWW[1]*aryuWW[168];
   aryuWW[168]=aryuWW[40]*aryuWW[168];
   aryuWW[168]=1./3.*aryuWW[170] + aryuWW[168];
   aryuWW[170]=aryuWW[20]*aryuWW[168];
   aryuWW[117]=aryuWW[119] + aryuWW[117] + aryuWW[123] + aryuWW[125] + 
   aryuWW[170] + aryuWW[130] + 1./2.*aryuWW[166] + aryuWW[160] + 
   aryuWW[153] + aryuWW[120];
   aryuWW[117]=aryuWW[165]*aryuWW[117];
   aryuWW[119]=aryuWW[128] - 2*aryuWW[57];
   aryuWW[120]=5*aryuWW[56];
   aryuWW[119]=aryuWW[59] + 3*aryuWW[54] + 2./3.*aryuWW[119] + 
   aryuWW[120];
   aryuWW[123]= - aryuWW[71] + 307./4. + aryuWW[69];
   aryuWW[125]= - aryuWW[74] + 1./2.*aryuWW[8];
   aryuWW[128]=aryuWW[110]*aryuWW[125];
   aryuWW[123]=1./6.*aryuWW[123] + aryuWW[128];
   aryuWW[123]=aryuWW[110]*aryuWW[123];
   aryuWW[128]= - 1 + aryuWW[137];
   aryuWW[128]=1./2.*aryuWW[128] - aryuWW[71];
   aryuWW[128]=MMZ*aryuWW[136]*aryuWW[128];
   aryuWW[130]=5./6.*aryuWW[66] - 1./8. + aryuWW[67];
   aryuWW[130]= - 11./12.*aryuWW[50] + 1./4.*aryuWW[69] + 1./2.*
   aryuWW[130] + aryuWW[33];
   aryuWW[130]=aryuWW[35]*aryuWW[130];
   aryuWW[137]= - 3067 - 123*aryuWW[66];
   aryuWW[137]= - 165*aryuWW[17] + 123./8.*aryuWW[50] + 123./8.*
   aryuWW[71] - 123./4.*aryuWW[69] + 1./8.*aryuWW[137] + 165*aryuWW[63]
   ;
   aryuWW[137]=aryuWW[111]*aryuWW[137];
   aryuWW[153]=aryuWW[16] - 6*aryuWW[68] - 15./4. - aryuWW[60];
   aryuWW[153]=aryuWW[112]*aryuWW[153];
   aryuWW[119]=aryuWW[128] + aryuWW[153] + 495./64.*aryuWW[190] + 3./64.
   *aryuWW[137] + 1./4.*aryuWW[123] + 2*aryuWW[119] + aryuWW[130];
   aryuWW[119]=MMZ*aryuWW[119];
   aryuWW[123]=8497./128. - aryuWW[47];
   aryuWW[128]=53./3.*aryuWW[13];
   aryuWW[130]=369./512.*aryuWW[111] - aryuWW[112];
   aryuWW[130]=MMZ*aryuWW[130];
   aryuWW[137]=aryuWW[311] - 369./64.*aryuWW[111];
   aryuWW[137]=aryuWW[19]*aryuWW[137];
   aryuWW[153]=aryuWW[11]*MMZ*aryuWW[112];
   aryuWW[123]=113./36.*aryuWW[37] + aryuWW[153] + 1./4.*aryuWW[137] + 
   aryuWW[130] + 20./9.*aryuWW[223] + 4./3.*aryuWW[259] + 19./12.*
   aryuWW[12] + 1./4.*aryuWW[123] + aryuWW[128];
   aryuWW[123]=aryuWW[37]*aryuWW[123];
   aryuWW[130]=aryuWW[147] + aryuWW[42] - aryuWW[8];
   aryuWW[137]= - 4*aryuWW[6] + 89./3. - 25*aryuWW[36];
   aryuWW[137]=MMZ*aryuWW[137];
   aryuWW[147]=aryuWW[151] + 11./3. - 29./2.*aryuWW[36];
   aryuWW[147]=1./6.*aryuWW[19]*aryuWW[147];
   aryuWW[160]= - aryuWW[37]*MMZ;
   aryuWW[166]=5./3. + aryuWW[198];
   aryuWW[170]=aryuWW[166] - aryuWW[6];
   aryuWW[174]=aryuWW[20]*aryuWW[170];
   aryuWW[175]= - 2*aryuWW[47];
   aryuWW[176]= - 47./3. + aryuWW[175];
   aryuWW[176]=aryuWW[38]*aryuWW[176];
   aryuWW[137]=aryuWW[176] + 4./3.*aryuWW[174] + 5*aryuWW[160] + 
   aryuWW[147] + aryuWW[130] + 1./3.*aryuWW[137];
   aryuWW[137]=aryuWW[3]*aryuWW[137];
   aryuWW[160]= - 7./3. + aryuWW[115];
   aryuWW[174]=10./3. - aryuWW[7];
   aryuWW[174]=aryuWW[6]*aryuWW[174];
   aryuWW[166]=aryuWW[7]*aryuWW[166];
   aryuWW[160]=aryuWW[174] + 2./3.*aryuWW[160] + aryuWW[166];
   aryuWW[160]=aryuWW[1]*aryuWW[160];
   aryuWW[115]= - 5./3. + aryuWW[115];
   aryuWW[174]=2 - aryuWW[7];
   aryuWW[174]=aryuWW[6]*aryuWW[174];
   aryuWW[115]=aryuWW[174] + 2./3.*aryuWW[115] + aryuWW[166];
   aryuWW[115]=aryuWW[40]*aryuWW[115];
   aryuWW[166]= - 1./2.*aryuWW[35] + aryuWW[110];
   aryuWW[174]=215./128.*aryuWW[111];
   aryuWW[176]=MMZ*aryuWW[136];
   aryuWW[166]=1./2.*aryuWW[176] + 1./3.*aryuWW[166] + aryuWW[174];
   aryuWW[166]=aryuWW[19]*aryuWW[166];
   aryuWW[179]=85./9.*aryuWW[35] + 369./32.*aryuWW[111];
   aryuWW[179]=aryuWW[37]*aryuWW[179];
   aryuWW[176]=1./2.*aryuWW[179] + 3./4.*aryuWW[176] + 495./8.*
   aryuWW[185] - 7387./64.*aryuWW[111] - 73./18.*aryuWW[35] + 
   aryuWW[110];
   aryuWW[176]=aryuWW[38]*aryuWW[176];
   aryuWW[179]= - 1./2.*aryuWW[37];
   aryuWW[180]= - 2 + aryuWW[169];
   aryuWW[180]=aryuWW[3]*aryuWW[38]*aryuWW[180];
   aryuWW[181]=6*aryuWW[180];
   aryuWW[183]=aryuWW[181] + aryuWW[179] - 77./3.*aryuWW[36] + 
   aryuWW[131] - 25./6. - aryuWW[47];
   aryuWW[183]=aryuWW[3]*aryuWW[183];
   aryuWW[184]=aryuWW[58] + 1./3.*aryuWW[52];
   aryuWW[187]=aryuWW[184] - 2./3.*aryuWW[57];
   aryuWW[198]=12*aryuWW[186];
   aryuWW[183]=aryuWW[198] + aryuWW[183] - 23./6.*aryuWW[54] + 2*
   aryuWW[187] + 3*aryuWW[56];
   aryuWW[183]=MMt*aryuWW[183];
   aryuWW[125]=1./3.*aryuWW[125] + aryuWW[298];
   aryuWW[125]=aryuWW[110]*aryuWW[125];
   aryuWW[187]=1./2.*aryuWW[125];
   aryuWW[199]=aryuWW[311] + aryuWW[322];
   aryuWW[199]=aryuWW[22]*aryuWW[199];
   aryuWW[202]=aryuWW[6] - 5./3. + 4*aryuWW[36];
   aryuWW[203]=aryuWW[37]*aryuWW[170];
   aryuWW[203]=1./3.*aryuWW[202] + aryuWW[203];
   aryuWW[203]=aryuWW[34]*aryuWW[203];
   aryuWW[208]=25./3.*aryuWW[53] - 722227./41472. - 8*aryuWW[64];
   aryuWW[138]= - aryuWW[74] + aryuWW[138];
   aryuWW[138]=1./3.*aryuWW[35]*aryuWW[138];
   aryuWW[209]=7./8.*aryuWW[42] + 11*aryuWW[72];
   aryuWW[209]= - 121*aryuWW[23] - 523./16.*aryuWW[22] + 11*aryuWW[209]
    - 1721./8.*aryuWW[51];
   aryuWW[209]=aryuWW[111]*aryuWW[209];
   aryuWW[210]= - 53./9.*aryuWW[13];
   aryuWW[212]=2809./48. + 235*aryuWW[36];
   aryuWW[212]=aryuWW[12]*aryuWW[212];
   aryuWW[216]= - 95./81. + 27*aryuWW[12];
   aryuWW[216]=aryuWW[6]*aryuWW[216];
   aryuWW[170]=aryuWW[11]*aryuWW[170];
   aryuWW[218]=aryuWW[20]*aryuWW[111];
   aryuWW[115]=4./9.*aryuWW[203] + aryuWW[183] + aryuWW[137] + 1./2.*
   aryuWW[176] + 121./32.*aryuWW[218] + aryuWW[123] + 4./9.*aryuWW[170]
    + 1./2.*aryuWW[166] + aryuWW[119] + 4./9.*aryuWW[115] + 4./27.*
   aryuWW[160] + 1./4.*aryuWW[216] + 1./12.*aryuWW[212] + 997./324.*
   aryuWW[36] + aryuWW[210] + 1./16.*aryuWW[209] + 1./4.*aryuWW[199] + 
   aryuWW[187] + 17./64.*aryuWW[17] + 79./36.*aryuWW[211] + aryuWW[138]
    + 1./2.*aryuWW[47] + 659./1536.*aryuWW[50] - 143./512.*aryuWW[71]
    + 4045./768.*aryuWW[69] - 5./24.*aryuWW[30] + 19./36.*aryuWW[33] + 
   1093./192.*aryuWW[63] - 659./1536.*aryuWW[66] + 1./3.*aryuWW[208] - 
   23./4.*aryuWW[61];
   aryuWW[115]=aryuWW[34]*aryuWW[115];
   aryuWW[119]=19./4.*aryuWW[12];
   aryuWW[123]=aryuWW[119] + aryuWW[163] + 179./4.*aryuWW[13] + 
   aryuWW[169] + 241./6. + aryuWW[172];
   aryuWW[123]=aryuWW[19]*aryuWW[123];
   aryuWW[137]= - 1 + 2*aryuWW[48];
   aryuWW[137]=2*aryuWW[137] - aryuWW[46];
   aryuWW[137]=aryuWW[159]*aryuWW[137];
   aryuWW[160]=aryuWW[277] - 1 + aryuWW[172];
   aryuWW[160]=aryuWW[20]*MMZ*aryuWW[160];
   aryuWW[166]=1 + aryuWW[277];
   aryuWW[170]=MMZ*aryuWW[166];
   aryuWW[176]=aryuWW[19]*aryuWW[170];
   aryuWW[137]=aryuWW[160] + aryuWW[137] + aryuWW[176];
   aryuWW[137]=aryuWW[3]*aryuWW[137];
   aryuWW[157]=205*aryuWW[13] + aryuWW[157] - 1975./18. + aryuWW[48];
   aryuWW[157]=125./2.*aryuWW[12] + 1./2.*aryuWW[157] + aryuWW[152];
   aryuWW[157]=MMZ*aryuWW[157];
   aryuWW[160]=aryuWW[164] - 17./2. + aryuWW[169];
   aryuWW[160]=aryuWW[18]*aryuWW[160];
   aryuWW[164]=aryuWW[173] - 23./2. + aryuWW[232];
   aryuWW[164]=MMZ*aryuWW[164];
   aryuWW[164]=aryuWW[164] + aryuWW[195];
   aryuWW[164]=aryuWW[11]*aryuWW[164];
   aryuWW[183]=53*aryuWW[13];
   aryuWW[195]=85./4.*aryuWW[12] + aryuWW[183] + 1151./12. + 
   aryuWW[169];
   aryuWW[195]=aryuWW[20]*aryuWW[195];
   aryuWW[199]= - 11*aryuWW[80] - 3./2.*aryuWW[78];
   aryuWW[199]=1./2.*aryuWW[199] + 7*aryuWW[22];
   aryuWW[123]=3*aryuWW[137] + aryuWW[195] + 1./2.*aryuWW[164] + 1./2.*
   aryuWW[123] + 1./4.*aryuWW[160] + aryuWW[157] + 33*aryuWW[23] + 3*
   aryuWW[199] + 217./24.*aryuWW[21];
   aryuWW[123]=aryuWW[3]*aryuWW[123];
   aryuWW[116]= - 55./4.*aryuWW[12] + aryuWW[321] + aryuWW[116] + 29./2.
    + aryuWW[288];
   aryuWW[116]=aryuWW[112]*aryuWW[116];
   aryuWW[116]=7./4.*aryuWW[113] + aryuWW[116];
   aryuWW[116]=MMZ*aryuWW[116];
   aryuWW[137]= - 2*aryuWW[46];
   aryuWW[157]=7./4.*aryuWW[178];
   aryuWW[116]=aryuWW[153] + aryuWW[157] + 4*aryuWW[122] + aryuWW[116]
    + 11./12.*aryuWW[12] + aryuWW[128] + 157./4. + aryuWW[137];
   aryuWW[116]=aryuWW[11]*aryuWW[116];
   aryuWW[122]= - 7*aryuWW[96] + aryuWW[162] + 7 + aryuWW[161];
   aryuWW[122]=aryuWW[142] + aryuWW[206] + 1./2.*aryuWW[122] + 
   aryuWW[204];
   aryuWW[122]=1./2.*aryuWW[113]*aryuWW[122];
   aryuWW[142]=aryuWW[231] - 40*aryuWW[16] + aryuWW[238] - 55./4.*
   aryuWW[17] + aryuWW[234] + 110*aryuWW[100] + aryuWW[197] + 11./2.*
   aryuWW[103] - 2*aryuWW[91] + 55./4.*aryuWW[92] + 55./4.*aryuWW[93]
    + 3*aryuWW[97] - 48*aryuWW[101] + 913./16. + 9*aryuWW[95];
   aryuWW[142]=aryuWW[112]*aryuWW[142];
   aryuWW[153]=2*aryuWW[90] + aryuWW[87];
   aryuWW[153]=11*aryuWW[153] - 357./16.*aryuWW[83];
   aryuWW[160]=11./4.*aryuWW[214];
   aryuWW[142]=aryuWW[142] + aryuWW[160] + aryuWW[122] + 55./2.*
   aryuWW[82] + 11*aryuWW[85] + 11./2.*aryuWW[109] - aryuWW[81] + 3*
   aryuWW[153] - 8*aryuWW[86];
   aryuWW[142]=MMZ*aryuWW[142];
   aryuWW[153]= - 3./4.*aryuWW[107];
   aryuWW[161]= - 1./4.*aryuWW[78];
   aryuWW[162]=2*aryuWW[22];
   aryuWW[164]=107./12.*aryuWW[23] - 35./3.*aryuWW[21] + aryuWW[162] + 
   aryuWW[161] + aryuWW[153] + 11*aryuWW[240] + 47./12.*aryuWW[79];
   aryuWW[164]=aryuWW[112]*aryuWW[164];
   aryuWW[118]=aryuWW[118] + 5./3.*aryuWW[86];
   aryuWW[118]= - 14./3.*aryuWW[82] + aryuWW[141] + 4*aryuWW[118] + 
   aryuWW[139];
   aryuWW[118]=MMH*aryuWW[118];
   aryuWW[139]=61./3.*aryuWW[113] + aryuWW[154];
   aryuWW[139]=1./2.*aryuWW[139];
   aryuWW[141]= - 293./3. + aryuWW[46];
   aryuWW[141]=aryuWW[112]*aryuWW[141];
   aryuWW[141]=aryuWW[139] + aryuWW[141];
   aryuWW[141]=aryuWW[18]*aryuWW[141];
   aryuWW[154]=15./4.*aryuWW[113];
   aryuWW[195]= - 341./12. - aryuWW[46];
   aryuWW[195]=aryuWW[112]*aryuWW[195];
   aryuWW[195]=aryuWW[154] + aryuWW[195];
   aryuWW[195]=aryuWW[20]*aryuWW[195];
   aryuWW[197]=168871./288. + 57*aryuWW[104];
   aryuWW[199]= - 15./2.*aryuWW[23] - 5./12.*aryuWW[21] - 29./6.*
   aryuWW[22] + 15./2.*aryuWW[107] - 7./3.*aryuWW[78];
   aryuWW[199]=aryuWW[113]*aryuWW[199];
   aryuWW[203]= - 12835./27. - 363*aryuWW[13];
   aryuWW[203]=1./4.*aryuWW[203] + aryuWW[10];
   aryuWW[203]=5*aryuWW[203] - 511./8.*aryuWW[12];
   aryuWW[203]=aryuWW[12]*aryuWW[203];
   aryuWW[204]= - 119./3.*aryuWW[113] + aryuWW[233];
   aryuWW[206]=aryuWW[204] - 33*aryuWW[112];
   aryuWW[206]=aryuWW[19]*aryuWW[206];
   aryuWW[116]=aryuWW[118] + aryuWW[123] + 1./2.*aryuWW[195] + 
   aryuWW[116] + 1./4.*aryuWW[206] + 1./4.*aryuWW[141] + aryuWW[142] + 
   aryuWW[164] + 1./2.*aryuWW[203] - 7./2.*aryuWW[10] - 2707./216.*
   aryuWW[13] + 1./2.*aryuWW[199] + 7./4.*aryuWW[46] - 121./12.*
   aryuWW[16] + aryuWW[318] + 3709./48.*aryuWW[17] + aryuWW[167] + 55./
   4.*aryuWW[100] + 19./6.*aryuWW[96] + 51./4.*aryuWW[103] - 11./12.*
   aryuWW[98] - 11./12.*aryuWW[92] - 11./12.*aryuWW[93] + aryuWW[217]
    - 1285./16.*aryuWW[99] + 215./8.*aryuWW[94] + 1./2.*aryuWW[197] + 8.
   /3.*aryuWW[102];
   aryuWW[116]=aryuWW[76]*aryuWW[116];
   aryuWW[118]= - 23*aryuWW[29];
   aryuWW[123]= - 37./2.*aryuWW[15];
   aryuWW[141]=7*aryuWW[32];
   aryuWW[142]= - 7*aryuWW[134];
   aryuWW[164]=aryuWW[142] + aryuWW[141] + aryuWW[123] - 35171./972. + 
   aryuWW[118];
   aryuWW[167]=193./4. + aryuWW[183];
   aryuWW[183]=7./4.*aryuWW[7];
   aryuWW[167]=1./3.*aryuWW[167] + aryuWW[183];
   aryuWW[167]=aryuWW[7]*aryuWW[167];
   aryuWW[121]=1./3.*aryuWW[220] + aryuWW[121];
   aryuWW[195]=aryuWW[1]*aryuWW[121];
   aryuWW[121]=aryuWW[40]*aryuWW[121];
   aryuWW[197]=10./3.*aryuWW[33];
   aryuWW[203]=19./2.*aryuWW[7];
   aryuWW[206]=227./3. + aryuWW[203];
   aryuWW[206]=aryuWW[12]*aryuWW[206];
   aryuWW[208]=451./54. + 79*aryuWW[12];
   aryuWW[208]=aryuWW[6]*aryuWW[208];
   aryuWW[121]=20./9.*aryuWW[121] + 56./27.*aryuWW[195] + 1./3.*
   aryuWW[208] + 1./6.*aryuWW[206] + aryuWW[167] + aryuWW[210] + 8*
   aryuWW[17] - 9./4.*aryuWW[30] + 1./4.*aryuWW[164] + aryuWW[197];
   aryuWW[121]=aryuWW[40]*aryuWW[121];
   aryuWW[118]=aryuWW[142] + aryuWW[141] + aryuWW[123] - 4003./108. + 
   aryuWW[118];
   aryuWW[118]=1./4.*aryuWW[118] + aryuWW[197];
   aryuWW[123]=aryuWW[183] + 67./4. + aryuWW[128];
   aryuWW[123]=aryuWW[7]*aryuWW[123];
   aryuWW[128]=8./3.*aryuWW[17];
   aryuWW[141]=139./3. + aryuWW[203];
   aryuWW[141]=aryuWW[12]*aryuWW[141];
   aryuWW[142]=155./18. + 41*aryuWW[12];
   aryuWW[142]=aryuWW[6]*aryuWW[142];
   aryuWW[118]=4./9.*aryuWW[195] + 1./3.*aryuWW[142] + 1./18.*
   aryuWW[141] + 1./3.*aryuWW[123] - 53./27.*aryuWW[13] + aryuWW[128]
    + 1./3.*aryuWW[118] - 3./4.*aryuWW[30];
   aryuWW[118]=aryuWW[1]*aryuWW[118];
   aryuWW[123]=7./3. - aryuWW[7];
   aryuWW[123]=1./3.*aryuWW[123] - aryuWW[6];
   aryuWW[123]=aryuWW[1]*aryuWW[123];
   aryuWW[141]= - 29./3.*aryuWW[6] + 89./9. - 5*aryuWW[7];
   aryuWW[141]=aryuWW[40]*aryuWW[141];
   aryuWW[123]=5*aryuWW[123] + aryuWW[141];
   aryuWW[123]=MMZ*aryuWW[123];
   aryuWW[140]=aryuWW[1]*aryuWW[140];
   aryuWW[140]=aryuWW[140] + 11./3.*aryuWW[223];
   aryuWW[140]=1./2.*aryuWW[19]*aryuWW[140];
   aryuWW[141]=aryuWW[259] + 5./3.*aryuWW[223];
   aryuWW[142]=aryuWW[20]*aryuWW[141];
   aryuWW[123]=4*aryuWW[142] + aryuWW[123] + aryuWW[140];
   aryuWW[123]=aryuWW[3]*aryuWW[123];
   aryuWW[142]=2*aryuWW[25];
   aryuWW[164]=2*aryuWW[24];
   aryuWW[167]=aryuWW[164] - aryuWW[26] + aryuWW[142];
   aryuWW[167]=2*aryuWW[167] + aryuWW[27];
   aryuWW[155]= - 1./3.*aryuWW[7] + aryuWW[158] + aryuWW[14] + 
   aryuWW[155] + 1 + 2./3.*aryuWW[31];
   aryuWW[155]=aryuWW[112]*aryuWW[155];
   aryuWW[155]=2./3.*aryuWW[167] + aryuWW[155];
   aryuWW[155]=aryuWW[1]*aryuWW[155];
   aryuWW[142]=aryuWW[164] + aryuWW[142] - 1./9.*aryuWW[75] - 
   aryuWW[26];
   aryuWW[142]=2*aryuWW[142] + aryuWW[27];
   aryuWW[158]= - aryuWW[7] + aryuWW[16] + 3*aryuWW[14] - aryuWW[28] + 
   3 + 2*aryuWW[31];
   aryuWW[158]=aryuWW[112]*aryuWW[158];
   aryuWW[142]=2*aryuWW[142] + aryuWW[158];
   aryuWW[142]=aryuWW[40]*aryuWW[142];
   aryuWW[142]=aryuWW[155] + aryuWW[142];
   aryuWW[142]=MMZ*aryuWW[142];
   aryuWW[155]=MMZ*aryuWW[168];
   aryuWW[141]=4./3.*aryuWW[141] + aryuWW[155];
   aryuWW[141]=aryuWW[11]*aryuWW[141];
   aryuWW[115]=aryuWW[117] + aryuWW[115] + aryuWW[116] + aryuWW[123] + 
   aryuWW[141] + aryuWW[142] + aryuWW[118] + aryuWW[121];
   aryuWW[115]=aryuWW[165]*aryuWW[115];
   aryuWW[116]=121./3. - 5*aryuWW[66];
   aryuWW[116]=5./8.*aryuWW[116] - 37*aryuWW[63];
   aryuWW[117]= - 1./8.*aryuWW[42] - aryuWW[72];
   aryuWW[117]=1./3.*aryuWW[23] - 1./48.*aryuWW[22] + 1./3.*aryuWW[117]
    + 7./8.*aryuWW[51];
   aryuWW[117]=aryuWW[111]*aryuWW[117];
   aryuWW[118]=101./24. + aryuWW[329];
   aryuWW[118]=aryuWW[12]*aryuWW[118];
   aryuWW[121]=aryuWW[6]*aryuWW[12];
   aryuWW[116]=5./3.*aryuWW[121] + 1./3.*aryuWW[118] + 25./2.*
   aryuWW[117] + 151./24.*aryuWW[17] + 25./64.*aryuWW[50] + 25./64.*
   aryuWW[71] + 565./96.*aryuWW[69] + 1./8.*aryuWW[116] - 5./3.*
   aryuWW[30];
   aryuWW[117]= - 3 - 1./3.*aryuWW[66];
   aryuWW[117]= - 1./9.*aryuWW[17] + 1./24.*aryuWW[50] + 1./24.*
   aryuWW[71] - 1./12.*aryuWW[69] + 1./8.*aryuWW[117] + 1./9.*
   aryuWW[63];
   aryuWW[117]=aryuWW[111]*aryuWW[117];
   aryuWW[117]=aryuWW[117] + 1./9.*aryuWW[190];
   aryuWW[117]=MMZ*aryuWW[117];
   aryuWW[118]=MMZ*aryuWW[111];
   aryuWW[123]=1 + aryuWW[118];
   aryuWW[123]=1./2.*aryuWW[123] + aryuWW[241];
   aryuWW[123]=aryuWW[37]*aryuWW[123];
   aryuWW[141]=aryuWW[37]*aryuWW[111];
   aryuWW[142]=1./8.*aryuWW[141] + 5./8.*aryuWW[111] + 1./3.*
   aryuWW[185];
   aryuWW[142]=aryuWW[38]*aryuWW[142];
   aryuWW[155]=aryuWW[19]*aryuWW[111];
   aryuWW[116]=25./6.*aryuWW[142] + 25./36.*aryuWW[327] + 25./96.*
   aryuWW[123] + 125./288.*aryuWW[155] + 1./3.*aryuWW[116] + 25./8.*
   aryuWW[117];
   aryuWW[116]=aryuWW[255]*aryuWW[116];
   aryuWW[117]=aryuWW[330] + 11 + aryuWW[329];
   aryuWW[117]=aryuWW[19]*aryuWW[117];
   aryuWW[123]= - 17*aryuWW[42] - 5*aryuWW[8];
   aryuWW[142]=35./4.*aryuWW[6] - 121./3. + 119./4.*aryuWW[36];
   aryuWW[142]=MMZ*aryuWW[142];
   aryuWW[155]=7./6. + aryuWW[47];
   aryuWW[155]=aryuWW[38]*aryuWW[155];
   aryuWW[117]=17*aryuWW[155] + 1./2.*aryuWW[117] + 1./3.*aryuWW[142]
    + aryuWW[191] + 1./2.*aryuWW[123] + 17*aryuWW[51];
   aryuWW[117]=aryuWW[3]*aryuWW[117];
   aryuWW[123]=2545./3. - 163*aryuWW[66];
   aryuWW[123]=1./8.*aryuWW[123] - 1733./9.*aryuWW[63];
   aryuWW[123]=1./4.*aryuWW[123] - 5./3.*aryuWW[33];
   aryuWW[142]=aryuWW[35]*aryuWW[205];
   aryuWW[123]=1./4.*aryuWW[125] + 2437./288.*aryuWW[17] + 11./4.*
   aryuWW[211] + 11./2.*aryuWW[142] - 17./4.*aryuWW[47] + 163./256.*
   aryuWW[50] + 2057./768.*aryuWW[71] + 151./384.*aryuWW[69] + 1./8.*
   aryuWW[123] - 22./9.*aryuWW[30];
   aryuWW[125]=1 + aryuWW[69];
   aryuWW[142]=aryuWW[35]*aryuWW[125];
   aryuWW[155]= - 1 - aryuWW[69];
   aryuWW[155]=aryuWW[110]*aryuWW[155];
   aryuWW[142]=1./12.*aryuWW[155] + aryuWW[246] + 11./2.*aryuWW[142];
   aryuWW[155]= - 31 - 7./3.*aryuWW[66];
   aryuWW[155]= - 37./9.*aryuWW[17] + 7./24.*aryuWW[50] + 7./24.*
   aryuWW[71] - 7./12.*aryuWW[69] + 1./8.*aryuWW[155] + 37./9.*
   aryuWW[63];
   aryuWW[155]=aryuWW[111]*aryuWW[155];
   aryuWW[142]=185./144.*aryuWW[190] + 1./3.*aryuWW[142] + 5./16.*
   aryuWW[155];
   aryuWW[142]=MMZ*aryuWW[142];
   aryuWW[118]=35./16.*aryuWW[241] + 35./32.*aryuWW[118] + 35./32. + 
   aryuWW[267];
   aryuWW[155]= - 1./3.*aryuWW[37];
   aryuWW[118]=1./4.*aryuWW[118] + aryuWW[155];
   aryuWW[118]=aryuWW[37]*aryuWW[118];
   aryuWW[141]=35./32.*aryuWW[141] + 185./12.*aryuWW[185] + aryuWW[331]
    + 3025./96.*aryuWW[111];
   aryuWW[141]=aryuWW[38]*aryuWW[141];
   aryuWW[158]= - 35./6. + aryuWW[267];
   aryuWW[158]=49./18.*aryuWW[36] + 1./3.*aryuWW[158] + aryuWW[131];
   aryuWW[158]=1./2.*aryuWW[158] + 3*aryuWW[180];
   aryuWW[158]=aryuWW[3]*aryuWW[158];
   aryuWW[164]= - 2*aryuWW[54] + aryuWW[246];
   aryuWW[158]=6*aryuWW[186] + 1./3.*aryuWW[164] + aryuWW[158];
   aryuWW[158]=MMt*aryuWW[158];
   aryuWW[164]=103./8.*aryuWW[42] - 47*aryuWW[72];
   aryuWW[164]=47./9.*aryuWW[23] - 227./144.*aryuWW[22] + 1./9.*
   aryuWW[164] + 63./8.*aryuWW[51];
   aryuWW[164]=aryuWW[111]*aryuWW[164];
   aryuWW[165]=2863./36. + 181*aryuWW[36];
   aryuWW[165]=aryuWW[12]*aryuWW[165];
   aryuWW[167]= - 5 + 61./4.*aryuWW[12];
   aryuWW[167]=aryuWW[6]*aryuWW[167];
   aryuWW[168]= - 925./192.*aryuWW[111] + aryuWW[276] + 1./6.*
   aryuWW[110];
   aryuWW[168]=aryuWW[19]*aryuWW[168];
   aryuWW[116]=1./8.*aryuWW[116] + aryuWW[158] + 1./3.*aryuWW[117] + 1./
   6.*aryuWW[141] + 235./144.*aryuWW[327] + 1./6.*aryuWW[118] + 1./6.*
   aryuWW[168] + 1./2.*aryuWW[142] + 1./18.*aryuWW[167] + 1./72.*
   aryuWW[165] - 17./18.*aryuWW[36] + 1./3.*aryuWW[123] + 5./8.*
   aryuWW[164];
   aryuWW[116]=aryuWW[255]*aryuWW[116];
   aryuWW[117]=11*aryuWW[50] - 13*aryuWW[69] - 19./4. - 11*aryuWW[66];
   aryuWW[117]=aryuWW[35]*aryuWW[117];
   aryuWW[117]=aryuWW[246] + 1./4.*aryuWW[117];
   aryuWW[118]=1./9.*aryuWW[71] - 5./12. - aryuWW[69];
   aryuWW[118]=1./2.*aryuWW[118] + 1./3.*aryuWW[126];
   aryuWW[118]=aryuWW[110]*aryuWW[118];
   aryuWW[123]=MMZ*aryuWW[136]*aryuWW[125];
   aryuWW[125]=3343 + 941./3.*aryuWW[66];
   aryuWW[125]=2201./9.*aryuWW[17] - 941./24.*aryuWW[50] - 941./24.*
   aryuWW[71] + 941./12.*aryuWW[69] + 1./8.*aryuWW[125] - 2201./9.*
   aryuWW[63];
   aryuWW[125]=aryuWW[111]*aryuWW[125];
   aryuWW[117]=1./24.*aryuWW[123] + 2201./576.*aryuWW[185] + 1./64.*
   aryuWW[125] + 1./3.*aryuWW[117] + 1./4.*aryuWW[118];
   aryuWW[117]=MMZ*aryuWW[117];
   aryuWW[118]= - MMZ*aryuWW[111];
   aryuWW[123]=aryuWW[276] + 941./64.*aryuWW[111];
   aryuWW[123]=aryuWW[19]*aryuWW[123];
   aryuWW[118]=aryuWW[155] + 1./3.*aryuWW[123] + 941./384.*aryuWW[118]
    + aryuWW[201] - 839./1152. - aryuWW[47];
   aryuWW[118]=aryuWW[37]*aryuWW[118];
   aryuWW[123]=8*aryuWW[6] + 11./3. + aryuWW[146];
   aryuWW[123]=MMZ*aryuWW[123];
   aryuWW[125]= - 29./9. + aryuWW[175];
   aryuWW[125]=aryuWW[38]*aryuWW[125];
   aryuWW[123]=aryuWW[125] + aryuWW[147] + aryuWW[130] + 1./9.*
   aryuWW[123];
   aryuWW[123]=aryuWW[3]*aryuWW[123];
   aryuWW[125]=aryuWW[333] - 941./32.*aryuWW[111];
   aryuWW[125]=aryuWW[37]*aryuWW[125];
   aryuWW[125]=1./6.*aryuWW[125] + 1./4.*aryuWW[193] + 2201./72.*
   aryuWW[190] + 17485./576.*aryuWW[111] - 7./6.*aryuWW[35] + 
   aryuWW[110];
   aryuWW[125]=aryuWW[38]*aryuWW[125];
   aryuWW[126]=aryuWW[181] + aryuWW[179] - 119./9.*aryuWW[36] + 
   aryuWW[131] + 149./18. - aryuWW[47];
   aryuWW[126]=aryuWW[3]*aryuWW[126];
   aryuWW[120]=aryuWW[171] + aryuWW[120] - 11./2.*aryuWW[54];
   aryuWW[120]=aryuWW[198] + 1./3.*aryuWW[120] + aryuWW[126];
   aryuWW[120]=MMt*aryuWW[120];
   aryuWW[126]=1./6.*aryuWW[193] + aryuWW[174] + 5./2.*aryuWW[35] + 1./
   3.*aryuWW[110];
   aryuWW[126]=aryuWW[19]*aryuWW[126];
   aryuWW[130]= - 71./1152.*aryuWW[50] - 143./128.*aryuWW[71] + 4711./
   576.*aryuWW[69] + 179./54.*aryuWW[30] - 103./27.*aryuWW[33] - 2963./
   432.*aryuWW[63] + 71./1152.*aryuWW[66] + 257647./10368. - 3*
   aryuWW[61];
   aryuWW[130]=1./2.*aryuWW[130] + aryuWW[47];
   aryuWW[131]=aryuWW[276] - 1./2.*aryuWW[110];
   aryuWW[131]=aryuWW[22]*aryuWW[131];
   aryuWW[141]= - 107*aryuWW[23] + 2527./16.*aryuWW[22] - 1067./8.*
   aryuWW[51] - 793./8.*aryuWW[42] + 107*aryuWW[72];
   aryuWW[141]=aryuWW[111]*aryuWW[141];
   aryuWW[142]= - 5183./432. - 53*aryuWW[36];
   aryuWW[142]=aryuWW[12]*aryuWW[142];
   aryuWW[146]= - 7./3. + aryuWW[12];
   aryuWW[146]=aryuWW[6]*aryuWW[146];
   aryuWW[116]=aryuWW[116] + aryuWW[120] + aryuWW[123] + 1./2.*
   aryuWW[125] + 107./96.*aryuWW[218] + 1./4.*aryuWW[118] + 1./2.*
   aryuWW[126] + aryuWW[117] + 1./12.*aryuWW[146] + 1./12.*aryuWW[142]
    + 29./36.*aryuWW[36] + 1./48.*aryuWW[141] + 1./12.*aryuWW[131] + 
   aryuWW[187] + 2827./1728.*aryuWW[17] + 3./4.*aryuWW[211] + 1./2.*
   aryuWW[130] + aryuWW[138];
   aryuWW[116]=aryuWW[255]*aryuWW[116];
   aryuWW[117]= - 5./3.*aryuWW[63] + 13./9. + aryuWW[64];
   aryuWW[117]=11./3.*aryuWW[17] + aryuWW[229] + 2*aryuWW[117] + 1./3.*
   aryuWW[33];
   aryuWW[118]= - 1./3. + aryuWW[144];
   aryuWW[118]=aryuWW[12]*aryuWW[118];
   aryuWW[120]=aryuWW[6]*aryuWW[249];
   aryuWW[117]=2./3.*aryuWW[120] + 4./3.*aryuWW[118] + aryuWW[188] + 2*
   aryuWW[117] + aryuWW[13];
   aryuWW[118]= - 1./3.*aryuWW[59] - aryuWW[54] - 2*aryuWW[56] + 
   aryuWW[184] + aryuWW[129];
   aryuWW[123]=aryuWW[128] + 7 - 8./3.*aryuWW[63];
   aryuWW[123]=aryuWW[111]*aryuWW[123];
   aryuWW[118]=8./3.*aryuWW[185] + 2./3.*aryuWW[118] + aryuWW[123];
   aryuWW[123]=1./8. + 1./3.*aryuWW[71];
   aryuWW[123]=MMZ*aryuWW[136]*aryuWW[123];
   aryuWW[118]=2*aryuWW[118] + aryuWW[123];
   aryuWW[118]=MMZ*aryuWW[118];
   aryuWW[123]= - 1 - aryuWW[13];
   aryuWW[125]=aryuWW[37]*aryuWW[123];
   aryuWW[126]=aryuWW[111] + aryuWW[190];
   aryuWW[126]=aryuWW[38]*aryuWW[126];
   aryuWW[128]=MMZ*aryuWW[202];
   aryuWW[128]=aryuWW[128] + 8*aryuWW[38];
   aryuWW[128]=aryuWW[3]*aryuWW[128];
   aryuWW[129]=MMt*aryuWW[3]*aryuWW[143];
   aryuWW[116]=aryuWW[116] + 64./9.*aryuWW[129] + 8./9.*aryuWW[128] + 
   64./3.*aryuWW[126] + 4*aryuWW[125] + 4./3.*aryuWW[117] + aryuWW[118]
   ;
   aryuWW[116]=aryuWW[34]*aryuWW[116];
   aryuWW[117]=aryuWW[145] + aryuWW[127] - 23./8.*aryuWW[83] - 2*
   aryuWW[90] + aryuWW[87];
   aryuWW[118]= - 3*aryuWW[103] + 31./4. + aryuWW[135];
   aryuWW[118]=aryuWW[16] - 3./2.*aryuWW[17] + 3*aryuWW[100] + 1./2.*
   aryuWW[118] - aryuWW[96];
   aryuWW[118]=aryuWW[113]*aryuWW[118];
   aryuWW[125]= - 2 + 1./4.*aryuWW[103];
   aryuWW[125]=aryuWW[112]*aryuWW[125];
   aryuWW[117]=1./3.*aryuWW[125] + 3./2.*aryuWW[149] + 1./3.*
   aryuWW[117] + aryuWW[118];
   aryuWW[117]=MMZ*aryuWW[117];
   aryuWW[118]=aryuWW[177] + aryuWW[182] + 7./12.*aryuWW[13] - 143./36.
    + aryuWW[46];
   aryuWW[118]=MMZ*aryuWW[118];
   aryuWW[125]= - aryuWW[80] - 9./2.*aryuWW[78];
   aryuWW[126]=5./3.*aryuWW[12] + 1./3. + aryuWW[239];
   aryuWW[126]=aryuWW[18]*aryuWW[126];
   aryuWW[118]=aryuWW[126] + aryuWW[118] + aryuWW[23] + 47./12.*
   aryuWW[21] + 1./2.*aryuWW[125] + 5*aryuWW[22];
   aryuWW[125]=aryuWW[13] + 25./6. + aryuWW[169];
   aryuWW[124]=1./8.*aryuWW[12] + 1./4.*aryuWW[125] + aryuWW[124];
   aryuWW[124]=aryuWW[19]*aryuWW[124];
   aryuWW[125]=aryuWW[159]*aryuWW[166];
   aryuWW[125]=aryuWW[125] + 3./2.*aryuWW[176];
   aryuWW[125]=aryuWW[3]*aryuWW[125];
   aryuWW[126]=1 + aryuWW[46];
   aryuWW[127]=3./2.*aryuWW[126] + aryuWW[12];
   aryuWW[127]=aryuWW[20]*aryuWW[127];
   aryuWW[118]=aryuWW[125] + 1./2.*aryuWW[127] + 1./2.*aryuWW[118] + 
   aryuWW[124];
   aryuWW[118]=aryuWW[3]*aryuWW[118];
   aryuWW[124]=aryuWW[324] + aryuWW[196];
   aryuWW[125]=MMZ*aryuWW[113];
   aryuWW[124]=aryuWW[178] + 1./4.*aryuWW[124] + aryuWW[125];
   aryuWW[124]=aryuWW[11]*aryuWW[124];
   aryuWW[127]= - 1 + aryuWW[173];
   aryuWW[127]=aryuWW[112]*aryuWW[127];
   aryuWW[127]=aryuWW[127] - 13./3.*aryuWW[113] + 3*aryuWW[214];
   aryuWW[127]=aryuWW[18]*aryuWW[127];
   aryuWW[128]= - 1./36.*aryuWW[99] + 33./8.*aryuWW[94] - 10477./288.
    - 7*aryuWW[104];
   aryuWW[129]= - 2*aryuWW[23] + 5./12.*aryuWW[21] - 37./6.*aryuWW[22]
    + 2*aryuWW[107] + 1./3.*aryuWW[78];
   aryuWW[129]=aryuWW[113]*aryuWW[129];
   aryuWW[130]= - 103./9. + 29./4.*aryuWW[13];
   aryuWW[130]= - 7./12.*aryuWW[12] + 1./3.*aryuWW[130] + 13*aryuWW[10]
   ;
   aryuWW[130]=aryuWW[12]*aryuWW[130];
   aryuWW[131]=3./4.*aryuWW[21] + 3./2.*aryuWW[22] + 1./3.*aryuWW[240]
    - 3./4.*aryuWW[78];
   aryuWW[131]=aryuWW[112]*aryuWW[131];
   aryuWW[135]=7./6.*aryuWW[112] - 25./3.*aryuWW[113] + 3*aryuWW[149];
   aryuWW[135]=aryuWW[19]*aryuWW[135];
   aryuWW[136]= - 5./3. + aryuWW[277];
   aryuWW[136]=aryuWW[112]*aryuWW[136];
   aryuWW[136]=aryuWW[113] + 1./2.*aryuWW[136];
   aryuWW[136]=aryuWW[20]*aryuWW[136];
   aryuWW[138]=1./2.*aryuWW[189] + 1./3.*aryuWW[82];
   aryuWW[138]=MMH*aryuWW[138];
   aryuWW[117]=1./2.*aryuWW[138] + aryuWW[118] + aryuWW[136] + 
   aryuWW[124] + 1./2.*aryuWW[135] + 1./4.*aryuWW[127] + aryuWW[117] + 
   1./2.*aryuWW[131] + 1./4.*aryuWW[130] - 13./6.*aryuWW[10] + 
   aryuWW[225] + aryuWW[129] + aryuWW[224] - 5./12.*aryuWW[16] - 1./4.*
   aryuWW[48] - 151./144.*aryuWW[17] + aryuWW[286] + 7./4.*aryuWW[100]
    + 5./12.*aryuWW[96] + 79./36.*aryuWW[103] + 1./2.*aryuWW[128] - 5./
   3.*aryuWW[98];
   aryuWW[117]=aryuWW[76]*aryuWW[117];
   aryuWW[118]=5./6.*aryuWW[103] - 1./2.*aryuWW[98] - 13./8.*aryuWW[99]
    + 7./2.*aryuWW[94] - 1489./96. - aryuWW[104];
   aryuWW[124]=1./6.*aryuWW[100];
   aryuWW[118]=aryuWW[253] - 11./24.*aryuWW[17] + aryuWW[124] + 1./3.*
   aryuWW[118] - 1./2.*aryuWW[96];
   aryuWW[127]= - aryuWW[103] + 23 + aryuWW[98];
   aryuWW[128]= - 1./12.*aryuWW[17];
   aryuWW[124]=aryuWW[16] + aryuWW[128] + aryuWW[124] + 1./12.*
   aryuWW[127] - aryuWW[96];
   aryuWW[124]=aryuWW[113]*aryuWW[124];
   aryuWW[124]=1./12.*aryuWW[149] + 3./16.*aryuWW[83] + aryuWW[124];
   aryuWW[124]=MMZ*aryuWW[124];
   aryuWW[127]= - 7./3. + 1./4.*aryuWW[13];
   aryuWW[127]=aryuWW[177] + 1./3.*aryuWW[127] + aryuWW[10];
   aryuWW[127]=aryuWW[12]*aryuWW[127];
   aryuWW[125]=aryuWW[178] + 1./4. + aryuWW[125];
   aryuWW[125]=aryuWW[11]*aryuWW[125];
   aryuWW[129]= - aryuWW[23] + 5./24.*aryuWW[21] - 7./4.*aryuWW[22] + 
   aryuWW[107] - 1./6.*aryuWW[78];
   aryuWW[129]=aryuWW[113]*aryuWW[129];
   aryuWW[130]=aryuWW[113] + 1./3.*aryuWW[214];
   aryuWW[130]=aryuWW[18]*aryuWW[130];
   aryuWW[131]= - 19*aryuWW[113] + aryuWW[149];
   aryuWW[131]=aryuWW[19]*aryuWW[131];
   aryuWW[135]=aryuWW[20]*aryuWW[113];
   aryuWW[136]= - 1 + aryuWW[12];
   aryuWW[136]=aryuWW[18]*aryuWW[136];
   aryuWW[136]=aryuWW[21] + aryuWW[136];
   aryuWW[136]=aryuWW[3]*aryuWW[136];
   aryuWW[138]=MMH*aryuWW[189];
   aryuWW[118]=1./12.*aryuWW[138] + 1./24.*aryuWW[136] + 1./2.*
   aryuWW[135] + aryuWW[125] + 1./12.*aryuWW[131] + 1./8.*aryuWW[130]
    + aryuWW[124] + 1./4.*aryuWW[127] + 1./2.*aryuWW[118] + aryuWW[129]
   ;
   aryuWW[118]=aryuWW[76]*aryuWW[118];
   aryuWW[124]= - aryuWW[32] - 1 - aryuWW[15];
   aryuWW[125]=1./2.*aryuWW[121] + aryuWW[177] + aryuWW[248] + 
   aryuWW[124] - 1./2.*aryuWW[30];
   aryuWW[125]=aryuWW[1]*aryuWW[125];
   aryuWW[121]=11./6.*aryuWW[121] + 11./18.*aryuWW[12] + 11./6.*
   aryuWW[17] + aryuWW[124] - 11./6.*aryuWW[30];
   aryuWW[121]=aryuWW[40]*aryuWW[121];
   aryuWW[121]=aryuWW[125] + 1./3.*aryuWW[121];
   aryuWW[124]=aryuWW[269] - 1 + 1./6.*aryuWW[94];
   aryuWW[124]=aryuWW[255]*aryuWW[76]*aryuWW[124];
   aryuWW[118]=1./8.*aryuWW[124] + 1./2.*aryuWW[121] + aryuWW[118];
   aryuWW[118]=aryuWW[255]*aryuWW[118];
   aryuWW[121]=39./8. - 1./9.*aryuWW[15];
   aryuWW[124]= - 1 + aryuWW[132];
   aryuWW[124]=aryuWW[6]*aryuWW[124];
   aryuWW[121]=11./9.*aryuWW[124] + 283./162.*aryuWW[12] - 1./18.*
   aryuWW[148] + 187./54.*aryuWW[17] - 187./54.*aryuWW[30] - 11./36.*
   aryuWW[33] + 1./18.*aryuWW[134] + 1./4.*aryuWW[121] - 2./9.*
   aryuWW[32];
   aryuWW[121]=aryuWW[40]*aryuWW[121];
   aryuWW[125]= - 11./3. + aryuWW[151];
   aryuWW[127]=aryuWW[1]*aryuWW[125];
   aryuWW[125]=aryuWW[40]*aryuWW[125];
   aryuWW[125]=aryuWW[127] + 11./9.*aryuWW[125];
   aryuWW[125]=MMZ*aryuWW[125];
   aryuWW[127]=1./2. + aryuWW[6];
   aryuWW[129]=aryuWW[1]*aryuWW[127];
   aryuWW[127]=aryuWW[40]*aryuWW[127];
   aryuWW[127]=3*aryuWW[129] + 11./3.*aryuWW[127];
   aryuWW[127]=aryuWW[19]*aryuWW[127];
   aryuWW[129]= - aryuWW[8] + aryuWW[22];
   aryuWW[130]=aryuWW[1]*aryuWW[129];
   aryuWW[129]=aryuWW[40]*aryuWW[129];
   aryuWW[125]=aryuWW[127] + aryuWW[125] + 3*aryuWW[130] + 11./3.*
   aryuWW[129];
   aryuWW[125]=aryuWW[3]*aryuWW[125];
   aryuWW[127]=119./8. + 5./3.*aryuWW[15];
   aryuWW[127]=1./2.*aryuWW[134] + 1./4.*aryuWW[127] - 2*aryuWW[32];
   aryuWW[124]=aryuWW[124] + 25./18.*aryuWW[12] - 1./6.*aryuWW[148] + 
   17./6.*aryuWW[17] - 17./6.*aryuWW[30] + 1./3.*aryuWW[127] - 1./4.*
   aryuWW[33];
   aryuWW[124]=aryuWW[1]*aryuWW[124];
   aryuWW[127]=aryuWW[1]*aryuWW[27];
   aryuWW[129]=aryuWW[40]*aryuWW[27];
   aryuWW[127]=aryuWW[127] + 1./3.*aryuWW[129];
   aryuWW[127]=MMZ*aryuWW[127];
   aryuWW[117]=aryuWW[118] + aryuWW[117] + aryuWW[125] + 1./6.*
   aryuWW[127] + aryuWW[124] + aryuWW[121];
   aryuWW[117]=aryuWW[255]*aryuWW[117];
   aryuWW[118]=aryuWW[119] + aryuWW[163] + 35./4.*aryuWW[13] + 
   aryuWW[169] + 49./6. + aryuWW[172];
   aryuWW[118]=aryuWW[19]*aryuWW[118];
   aryuWW[119]=2 + aryuWW[324];
   aryuWW[119]=aryuWW[159]*aryuWW[119];
   aryuWW[121]=aryuWW[20]*aryuWW[170];
   aryuWW[119]=3*aryuWW[121] + aryuWW[119] + 3*aryuWW[176];
   aryuWW[119]=aryuWW[3]*aryuWW[119];
   aryuWW[121]=13./6.*aryuWW[12] + aryuWW[152] + 31./6.*aryuWW[13] + 
   aryuWW[231] - 367./36. + aryuWW[48];
   aryuWW[121]=MMZ*aryuWW[121];
   aryuWW[124]= - aryuWW[80] + aryuWW[192];
   aryuWW[124]=1./2.*aryuWW[124] + aryuWW[22];
   aryuWW[125]=5./2. + aryuWW[46];
   aryuWW[125]=3*aryuWW[125] + 1./2.*aryuWW[12];
   aryuWW[125]=aryuWW[18]*aryuWW[125];
   aryuWW[126]=MMZ*aryuWW[126];
   aryuWW[126]=aryuWW[126] + aryuWW[19];
   aryuWW[126]=aryuWW[11]*aryuWW[126];
   aryuWW[129]=37./4.*aryuWW[12] + 127./12. + aryuWW[169];
   aryuWW[129]=aryuWW[20]*aryuWW[129];
   aryuWW[118]=aryuWW[119] + aryuWW[129] + 3./4.*aryuWW[126] + 1./2.*
   aryuWW[118] + 1./4.*aryuWW[125] + aryuWW[121] + 9*aryuWW[23] + 9*
   aryuWW[124] + 19./8.*aryuWW[21];
   aryuWW[118]=aryuWW[3]*aryuWW[118];
   aryuWW[119]= - 11./2.*aryuWW[98] + aryuWW[264] + aryuWW[260] - 
   aryuWW[97] - 445./24.*aryuWW[99] - 331./4.*aryuWW[94] + 159221./288.
    + 59*aryuWW[104];
   aryuWW[121]=7753./27. + 25*aryuWW[13];
   aryuWW[121]= - 557./24.*aryuWW[12] + 1./4.*aryuWW[121] + 5*
   aryuWW[10];
   aryuWW[121]=aryuWW[12]*aryuWW[121];
   aryuWW[119]=aryuWW[121] - 7*aryuWW[10] - 35./12.*aryuWW[13] + 
   aryuWW[199] + aryuWW[133] - 37./6.*aryuWW[16] - aryuWW[48] + 727./72.
   *aryuWW[17] + aryuWW[221] + 7./2.*aryuWW[100] + 19./3.*aryuWW[96] + 
   1./3.*aryuWW[119] + 35./2.*aryuWW[103];
   aryuWW[121]=11./12.*aryuWW[23] + 1./3.*aryuWW[21] + aryuWW[162] + 
   aryuWW[161] + aryuWW[153] + 3*aryuWW[240] - 1./12.*aryuWW[79];
   aryuWW[121]=aryuWW[112]*aryuWW[121];
   aryuWW[124]=1./3.*aryuWW[92] + 21./4. + 1./3.*aryuWW[93];
   aryuWW[124]=aryuWW[219] + 1./2.*aryuWW[124] + 3*aryuWW[103];
   aryuWW[124]=aryuWW[213] + 2./3.*aryuWW[16] + aryuWW[128] + 1./2.*
   aryuWW[124] + 2./3.*aryuWW[100];
   aryuWW[124]=aryuWW[112]*aryuWW[124];
   aryuWW[122]=aryuWW[124] + aryuWW[160] + aryuWW[122] + 1./6.*
   aryuWW[82] + 3*aryuWW[85] + 3./2.*aryuWW[109] - 31./16.*aryuWW[83]
    - 34./3.*aryuWW[90] + 5*aryuWW[87];
   aryuWW[122]=MMZ*aryuWW[122];
   aryuWW[124]=aryuWW[169] + aryuWW[196];
   aryuWW[124]=aryuWW[112]*aryuWW[124];
   aryuWW[124]=7*aryuWW[113] + aryuWW[124];
   aryuWW[124]=MMZ*aryuWW[124];
   aryuWW[124]=aryuWW[157] + 1./4.*aryuWW[124] + aryuWW[156] - 3./4. + 
   aryuWW[137];
   aryuWW[124]=aryuWW[11]*aryuWW[124];
   aryuWW[125]= - 53./3. + aryuWW[46];
   aryuWW[125]=aryuWW[112]*aryuWW[125];
   aryuWW[125]=aryuWW[139] + aryuWW[125];
   aryuWW[125]=aryuWW[18]*aryuWW[125];
   aryuWW[126]= - 53./12. - aryuWW[46];
   aryuWW[126]=aryuWW[112]*aryuWW[126];
   aryuWW[126]=aryuWW[154] + aryuWW[126];
   aryuWW[126]=aryuWW[20]*aryuWW[126];
   aryuWW[128]=aryuWW[204] - aryuWW[112];
   aryuWW[128]=aryuWW[19]*aryuWW[128];
   aryuWW[129]=11./12.*aryuWW[150] + 2*aryuWW[82];
   aryuWW[129]=MMH*aryuWW[129];
   aryuWW[118]=aryuWW[129] + aryuWW[118] + 1./2.*aryuWW[126] + 
   aryuWW[124] + 1./4.*aryuWW[128] + 1./4.*aryuWW[125] + aryuWW[122] + 
   1./2.*aryuWW[119] + aryuWW[121];
   aryuWW[118]=aryuWW[76]*aryuWW[118];
   aryuWW[119]=aryuWW[134] - 13./3.*aryuWW[32] - 55./18.*aryuWW[15] - 
   563./36. - aryuWW[29];
   aryuWW[121]=7./2.*aryuWW[7];
   aryuWW[122]= - 253./3. + aryuWW[121];
   aryuWW[122]=aryuWW[12]*aryuWW[122];
   aryuWW[124]=1./2. + aryuWW[194];
   aryuWW[124]=aryuWW[6]*aryuWW[124];
   aryuWW[125]= - 1./3. - aryuWW[7];
   aryuWW[125]=aryuWW[7]*aryuWW[125];
   aryuWW[119]=aryuWW[124] + 1./18.*aryuWW[122] + 1./4.*aryuWW[125] - 
   19./9.*aryuWW[17] + 85./36.*aryuWW[30] + 1./4.*aryuWW[119] - 10./3.*
   aryuWW[33];
   aryuWW[119]=aryuWW[1]*aryuWW[119];
   aryuWW[122]=1./3.*aryuWW[134] - 23./3.*aryuWW[32] - 229./18.*
   aryuWW[15] - 7993./324. - 3*aryuWW[29];
   aryuWW[121]= - 781./27. + aryuWW[121];
   aryuWW[121]=aryuWW[12]*aryuWW[121];
   aryuWW[124]=11./6. + aryuWW[247];
   aryuWW[124]=aryuWW[6]*aryuWW[124];
   aryuWW[121]=1./3.*aryuWW[124] + 1./6.*aryuWW[121] + 1./12.*
   aryuWW[125] - 59./27.*aryuWW[17] + 317./108.*aryuWW[30] + 1./4.*
   aryuWW[122] - 110./27.*aryuWW[33];
   aryuWW[121]=aryuWW[40]*aryuWW[121];
   aryuWW[122]=MMZ*aryuWW[281];
   aryuWW[122]=aryuWW[122] + aryuWW[140];
   aryuWW[122]=aryuWW[3]*aryuWW[122];
   aryuWW[117]=aryuWW[117] + aryuWW[118] + aryuWW[122] + 1./3.*
   aryuWW[127] + aryuWW[119] + aryuWW[121];
   aryuWW[117]=aryuWW[255]*aryuWW[117];
   aryuWW[118]=pow(CW,2);
   aryuWW[119]=17./3. - 2*aryuWW[118];
   aryuWW[121]= - 4*aryuWW[118];
   aryuWW[122]= - 53./3. + aryuWW[121];
   aryuWW[122]=aryuWW[13]*aryuWW[122];
   aryuWW[119]= - 10*aryuWW[12] + 2*aryuWW[119] + aryuWW[122];
   aryuWW[119]=MMZ*aryuWW[119];
   aryuWW[122]=2 + aryuWW[232];
   aryuWW[122]=aryuWW[3]*aryuWW[159]*aryuWW[122];
   aryuWW[124]= - 5 - 3*aryuWW[13];
   aryuWW[124]=aryuWW[20]*aryuWW[124];
   aryuWW[119]=aryuWW[122] + aryuWW[119] + 2*aryuWW[124];
   aryuWW[119]=aryuWW[3]*aryuWW[119];
   aryuWW[122]=4*aryuWW[118];
   aryuWW[124]=13 + aryuWW[122];
   aryuWW[122]=113./3. + aryuWW[122];
   aryuWW[122]=aryuWW[13]*aryuWW[122];
   aryuWW[122]=6*aryuWW[12] + 5*aryuWW[124] + aryuWW[122];
   aryuWW[122]=aryuWW[12]*aryuWW[122];
   aryuWW[124]=aryuWW[83] + aryuWW[88] - 2*aryuWW[87];
   aryuWW[124]=aryuWW[118]*aryuWW[124];
   aryuWW[124]= - aryuWW[82] + 2*aryuWW[124] + aryuWW[86] + 23./3.*
   aryuWW[83] - 4*aryuWW[90] - 23./3.*aryuWW[87];
   aryuWW[125]=20*aryuWW[16] + 4*aryuWW[17] - 32*aryuWW[100] - 4*
   aryuWW[92] - 55 - 4*aryuWW[93];
   aryuWW[125]=aryuWW[112]*aryuWW[125];
   aryuWW[124]=8*aryuWW[124] + aryuWW[125];
   aryuWW[124]=MMZ*aryuWW[124];
   aryuWW[125]= - 4 + aryuWW[12];
   aryuWW[125]=MMZ*aryuWW[112]*aryuWW[125];
   aryuWW[125]=aryuWW[123] + aryuWW[125];
   aryuWW[125]=aryuWW[11]*aryuWW[125];
   aryuWW[126]= - 3637./6. - 88*aryuWW[104];
   aryuWW[127]=61./3. + 8*aryuWW[118];
   aryuWW[127]=aryuWW[99]*aryuWW[127];
   aryuWW[128]= - 67./3. - 8*aryuWW[118];
   aryuWW[128]=aryuWW[17]*aryuWW[128];
   aryuWW[121]= - 73./9. + aryuWW[121];
   aryuWW[121]=aryuWW[13]*aryuWW[121];
   aryuWW[118]=2*aryuWW[119] + 4*aryuWW[125] + aryuWW[124] + 2*
   aryuWW[122] + 2*aryuWW[121] + 4*aryuWW[128] - 44./3.*aryuWW[33] + 4*
   aryuWW[127] - 8*aryuWW[94] + 1./3.*aryuWW[126] - 64*aryuWW[118];
   aryuWW[118]=aryuWW[76]*aryuWW[118];
   aryuWW[119]=aryuWW[7]*aryuWW[123];
   aryuWW[121]=2*aryuWW[120] - 8./3.*aryuWW[12] + aryuWW[119] + 
   aryuWW[200] + 2*aryuWW[17] - 2*aryuWW[30] + 11 + 2*aryuWW[33];
   aryuWW[121]=aryuWW[1]*aryuWW[121];
   aryuWW[122]=10*aryuWW[17] - 10*aryuWW[30] + 169./3. + 10*aryuWW[33];
   aryuWW[122]=1./3.*aryuWW[122] + aryuWW[13];
   aryuWW[119]=10./9.*aryuWW[120] - 40./27.*aryuWW[12] + 1./3.*
   aryuWW[122] + aryuWW[119];
   aryuWW[119]=aryuWW[40]*aryuWW[119];
   aryuWW[119]=1./3.*aryuWW[121] + aryuWW[119];
   aryuWW[120]= - aryuWW[24] + 2*aryuWW[26] - aryuWW[25];
   aryuWW[120]=aryuWW[1]*aryuWW[120];
   aryuWW[121]=1./9.*aryuWW[75] + aryuWW[26];
   aryuWW[121]= - 2./9.*aryuWW[27] - aryuWW[24] + 2*aryuWW[121] - 
   aryuWW[25];
   aryuWW[121]=aryuWW[40]*aryuWW[121];
   aryuWW[120]=1./3.*aryuWW[120] + aryuWW[121];
   aryuWW[120]=MMZ*aryuWW[120];
   aryuWW[121]=aryuWW[215] + 5./3.*aryuWW[207];
   aryuWW[121]=aryuWW[3]*MMZ*aryuWW[121];
   aryuWW[119]=4./3.*aryuWW[121] + 2*aryuWW[119] + aryuWW[120];

      yuWWret = aryuWW[114] + aryuWW[115] + aryuWW[116] + aryuWW[117]
       + aryuWW[118] + 2*aryuWW[119];
      return yuWWret;
}
