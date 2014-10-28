#include <tt.hpp>
std::complex<long double>
tt::m20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armttbar[311], mttbarret;

    armttbar[1]=double(nL + nH);
    armttbar[2]=pow(CW,-1);
    armttbar[3]=pow(MMH,-1);
    armttbar[4]=pow(SW,-1);
    armttbar[5]=Tsil::I2(0,0,MMZ,mu2);
    armttbar[6]=pow(MMt,-1);
    armttbar[7]=Tsil::I2(0,0,MMW,mu2);
    armttbar[8]=Tsil::B(MMH,MMt,MMt,mu2);
    armttbar[9]=pow(MMZ,-1);
    armttbar[10]=std::real(Tsil::B(0,0,MMZ,mu2));
    armttbar[11]=std::real(Tsil::B(0,0,MMW,mu2));
    armttbar[12]=Tsil::B(MMZ,MMt,MMt,mu2);
    armttbar[13]=Tsil::Beps(MMZ,MMt,MMt,mu2);
    armttbar[14]=Tsil::A(MMH,mu2);
    armttbar[15]=Tsil::A(MMZ,mu2);
    armttbar[16]=Tsil::A(MMW,mu2);
    armttbar[17]=Tsil::A(MMt,mu2);
    armttbar[18]=Tsil::Aeps(MMZ,mu2);
    armttbar[19]=Tsil::Aeps(MMW,mu2);
    armttbar[20]=Tsil::Aeps(MMt,mu2);
    armttbar[21]=std::real(Tsil::B(0,MMW,MMt,mu2));
    armttbar[22]=std::real(Tsil::B(0,0,MMt,mu2));
    armttbar[23]=std::real(Tsil::Beps(0,MMW,MMt,mu2));
    armttbar[24]=protZ0tW0->Uzxyv(0);
    armttbar[25]=protWt000->Tyzv(0);
    armttbar[26]=prot0W00->Uxzuv(0);
    armttbar[27]=double(nH);
    armttbar[28]=Tsil::B(MMt,MMt,MMH,mu2);
    armttbar[29]=Tsil::B(MMt,MMt,MMZ,mu2);
    armttbar[30]=Tsil::B(0,MMt,MMW,mu2);
    armttbar[31]=double(nL);
    armttbar[32]=Tsil::I2(MMH,MMt,MMt,mu2);
    armttbar[33]=Tsil::I2(MMZ,MMt,MMt,mu2);
    armttbar[34]=Tsil::I2(0,MMW,MMt,mu2);
    armttbar[35]=Tsil::B(MMH,MMH,MMH,mu2);
    armttbar[36]=Tsil::B(MMZ,MMH,MMZ,mu2);
    armttbar[37]=Tsil::B(MMZ,MMZ,MMH,mu2);
    armttbar[38]=Tsil::B(MMW,MMH,MMW,mu2);
    armttbar[39]=Tsil::B(MMW,MMZ,MMW,mu2);
    armttbar[40]=Tsil::B(MMW,MMW,MMH,mu2);
    armttbar[41]=Tsil::B(MMW,MMW,MMZ,mu2);
    armttbar[42]=Tsil::Aeps(MMH,mu2);
    armttbar[43]=protWt000->Uzxyv(0);
    armttbar[44]=prot0ttHt->Tuxv(0);
    armttbar[45]=prot0ttZt->Tuxv(0);
    armttbar[46]=double(boson);
    armttbar[47]=Tsil::I2(MMH,MMH,MMH,mu2);
    armttbar[48]=Tsil::I2(MMH,MMZ,MMZ,mu2);
    armttbar[49]=Tsil::I2(MMH,MMW,MMW,mu2);
    armttbar[50]=Tsil::I2(MMW,MMW,MMZ,mu2);
    armttbar[51]=Tsil::Beps(MMH,MMt,MMt,mu2);
    armttbar[52]=protWt000->M(0);
    armttbar[53]=prot0ttHt->M(0);
    armttbar[54]=prot0ttZt->M(0);
    armttbar[55]=prot0tt0t->M(0);
    armttbar[56]=prottH0H->Vxzuv(0);
    armttbar[57]=prottZ0Z->Vxzuv(0);
    armttbar[58]=protWt000->Txuv(0);
    armttbar[59]=protHHttH->M(0);
    armttbar[60]=protHZttZ->M(0);
    armttbar[61]=protHWt0W->M(0);
    armttbar[62]=protHttHt->M(0);
    armttbar[63]=protHttZt->M(0);
    armttbar[64]=protZZttH->M(0);
    armttbar[65]=protZWt0W->M(0);
    armttbar[66]=protZttZt->M(0);
    armttbar[67]=protZ0tW0->M(0);
    armttbar[68]=protWW00Z->M(0);
    armttbar[69]=protW00tW->M(0);
    armttbar[70]=prot00WW0->M(0);
    armttbar[71]=protHHttH->Uxzuv(0);
    armttbar[72]=protHZttZ->Uxzuv(0);
    armttbar[73]=protHWt0W->Uxzuv(0);
    armttbar[74]=protHttZt->Uuyxv(0);
    armttbar[75]=protHZttZ->Uyuzv(0);
    armttbar[76]=protZ0tW0->Uxzuv(0);
    armttbar[77]=protHZttZ->Uzxyv(0);
    armttbar[78]=protHWt0W->Uzxyv(0);
    armttbar[79]=protZZttH->Uzxyv(0);
    armttbar[80]=protZWt0W->Uzxyv(0);
    armttbar[81]=protW00tW->Uuyxv(0);
    armttbar[82]=protHWt0W->Uuyxv(0);
    armttbar[83]=protZWt0W->Uuyxv(0);
    armttbar[84]=protHZttZ->Txuv(0);
    armttbar[85]=protHWt0W->Txuv(0);
    armttbar[86]=protHttHt->Txuv(0);
    armttbar[87]=protZWt0W->Txuv(0);
    armttbar[88]=protHttZt->Tuxv(0);
    armttbar[89]=protHWt0W->Tvxu(0);
    armttbar[90]=protZWt0W->Tvxu(0);
    armttbar[91]=protHHttH->Suxv(0);
    armttbar[92]=protHZttZ->Suxv(0);
    armttbar[93]=protHZttZ->Svyz(0);
    armttbar[94]=protHWt0W->Svyz(0);
    armttbar[95]=protHWt0W->Suxv(0);
    armttbar[96]=protZ0tW0->Suxv(0);
    armttbar[97]=protHHttH->Uuyxv(0);
    armttbar[98]=1/(MMt - MMW);
    armttbar[99]=1/(4*MMt - MMZ);
    armttbar[100]=1/( - 4*MMW + MMH);
    armttbar[101]=1/( - 4*MMZ + MMH);
   armttbar[102]=55*armttbar[39];
   armttbar[103]= - 209./2.*armttbar[41] - 367./6. + armttbar[102];
   armttbar[104]=armttbar[14]*armttbar[101];
   armttbar[105]=1./2.*armttbar[104];
   armttbar[106]=5./16.*armttbar[12];
   armttbar[107]=1./2.*armttbar[36];
   armttbar[103]=armttbar[106] - 15./16.*armttbar[8] - 23./64.*
   armttbar[21] + armttbar[105] + armttbar[107] + 1./4.*armttbar[103]
    - armttbar[38];
   armttbar[103]=armttbar[12]*armttbar[103];
   armttbar[108]= - 33*armttbar[41];
   armttbar[109]= - 133./8. + armttbar[108];
   armttbar[110]=27./32.*armttbar[21];
   armttbar[109]=armttbar[110] + 1./4.*armttbar[109] + armttbar[36];
   armttbar[109]=armttbar[99]*armttbar[109];
   armttbar[111]=1./2.*armttbar[98] + armttbar[101];
   armttbar[112]=armttbar[8]*armttbar[101];
   armttbar[113]=3*armttbar[112];
   armttbar[114]= - armttbar[12]*armttbar[101];
   armttbar[109]=3*armttbar[109] + armttbar[114] + 7*armttbar[111] + 
   armttbar[113];
   armttbar[109]=armttbar[15]*armttbar[109];
   armttbar[111]= - 3*armttbar[93];
   armttbar[115]=1./2.*armttbar[42];
   armttbar[116]=3*armttbar[20];
   armttbar[117]=armttbar[116] + 5*armttbar[18] + armttbar[111] + 
   armttbar[115];
   armttbar[117]=armttbar[101]*armttbar[117];
   armttbar[118]= - 11*armttbar[41];
   armttbar[119]=11*armttbar[39];
   armttbar[120]=armttbar[14]*armttbar[100];
   armttbar[121]= - 25./32.*armttbar[21] + 1./4.*armttbar[120] - 
   armttbar[38] + armttbar[118] - 2099./384. + armttbar[119];
   armttbar[121]=armttbar[21]*armttbar[121];
   armttbar[122]=3./2.*armttbar[50] + armttbar[94];
   armttbar[123]=5*armttbar[92];
   armttbar[122]=armttbar[123] + 3*armttbar[93] + 11*armttbar[122] - 27.
   /4.*armttbar[96];
   armttbar[124]= - 5./2.*armttbar[14];
   armttbar[125]= - 1./4.*armttbar[32];
   armttbar[122]= - 237./8.*armttbar[16] + armttbar[124] - 207./16.*
   armttbar[20] - 337./32.*armttbar[18] - 55./2.*armttbar[19] + 9./8.*
   armttbar[33] + 27./16.*armttbar[34] + armttbar[125] + 1./2.*
   armttbar[122] - armttbar[48];
   armttbar[122]=armttbar[99]*armttbar[122];
   armttbar[126]=9./2.*armttbar[19];
   armttbar[127]=1./2.*armttbar[49];
   armttbar[128]=armttbar[126] - 1./4.*armttbar[42] + armttbar[127] - 3
   *armttbar[94];
   armttbar[128]=armttbar[100]*armttbar[128];
   armttbar[129]= - 1./2.*armttbar[37];
   armttbar[130]= - 1./4.*armttbar[21];
   armttbar[131]=armttbar[130] + armttbar[129] - 1 - armttbar[40];
   armttbar[131]=armttbar[8]*armttbar[131];
   armttbar[132]= - 9./4.*armttbar[98] - armttbar[100];
   armttbar[132]=armttbar[21]*armttbar[132];
   armttbar[132]=armttbar[132] - 171./2.*armttbar[98] + 11*
   armttbar[100];
   armttbar[133]=armttbar[8]*armttbar[100];
   armttbar[132]=1./2.*armttbar[132] + 3*armttbar[133];
   armttbar[132]=armttbar[16]*armttbar[132];
   armttbar[133]=1./2.*armttbar[24];
   armttbar[134]=7./8.*armttbar[5] + 9./8.*armttbar[34] + 15*
   armttbar[94] - 7./4.*armttbar[96];
   armttbar[134]=armttbar[98]*armttbar[134];
   armttbar[135]= - armttbar[19]*armttbar[98];
   armttbar[136]=armttbar[18]*armttbar[98];
   armttbar[137]= - 43./8.*armttbar[98] + armttbar[100];
   armttbar[137]=armttbar[20]*armttbar[137];
   armttbar[138]= - 3*armttbar[100] - armttbar[101];
   armttbar[138]=1./8.*armttbar[14]*armttbar[138];
   armttbar[139]=3./2.*armttbar[40];
   armttbar[140]=3./4.*armttbar[37];
   armttbar[141]=pow(Pi,2);
   armttbar[103]=1./4.*armttbar[109] + 3./4.*armttbar[122] + 1./2.*
   armttbar[132] + 1./4.*armttbar[103] + 3./2.*armttbar[131] + 1./2.*
   armttbar[121] + armttbar[138] - 13./8.*armttbar[36] + 1./4.*
   armttbar[38] + 1./4.*armttbar[117] + 561./32.*armttbar[41] + 
   armttbar[140] + 3./2.*armttbar[137] + 3./16.*armttbar[136] - 231./16.
   *armttbar[39] + 1./2.*armttbar[128] - 347./256.*armttbar[13] + 
   armttbar[139] + 123./8.*armttbar[135] + 1./2.*armttbar[134] + 
   armttbar[133] - 1./2.*armttbar[25] + 9./4.*armttbar[23] - 63./64.*
   armttbar[51] - armttbar[45] - 1./8.*armttbar[79] - 59./96.*
   armttbar[88] - 27./8.*armttbar[35] + 365./128.*armttbar[87] + 7./64.
   *armttbar[72] + 7./64.*armttbar[74] + 1./16.*armttbar[75] + 1./4.*
   armttbar[77] + 195./64.*armttbar[90] + 3./32.*armttbar[141] + 1./4.*
   armttbar[82] - 73./96.*armttbar[89] + 31./32.*armttbar[80] + 1./2.*
   armttbar[43] - 41./256.*armttbar[76] + 1./2.*armttbar[78] + 1./8.*
   armttbar[73] + 28655./3072. - 3*armttbar[83];
   armttbar[109]=1./2. - armttbar[40];
   armttbar[117]=armttbar[109] + armttbar[129];
   armttbar[121]=armttbar[16]*armttbar[117];
   armttbar[117]=armttbar[15]*armttbar[117];
   armttbar[117]=armttbar[121] + 1./2.*armttbar[117];
   armttbar[117]=armttbar[3]*armttbar[117];
   armttbar[121]=3*armttbar[40];
   armttbar[122]=3./2.*armttbar[37];
   armttbar[128]=9./2.*armttbar[36];
   armttbar[131]=3./8.*armttbar[12];
   armttbar[117]=9*armttbar[117] + armttbar[131] + 3./8.*armttbar[21]
    + armttbar[128] - 825./8.*armttbar[41] + armttbar[122] + 33./2.*
   armttbar[39] + 331./12. + armttbar[121];
   armttbar[117]=armttbar[3]*armttbar[117];
   armttbar[132]=3*armttbar[78];
   armttbar[134]=1./2.*armttbar[82];
   armttbar[135]= - 1./2.*armttbar[23];
   armttbar[137]= - 1./2.*armttbar[89];
   armttbar[142]=armttbar[135] - 3*armttbar[51] + armttbar[85] + 
   armttbar[134] + armttbar[137] - 13./2. + armttbar[132];
   armttbar[142]=armttbar[100]*armttbar[142];
   armttbar[143]=3./2.*armttbar[77];
   armttbar[144]= - 3./2.*armttbar[51];
   armttbar[145]= - 1./2.*armttbar[13];
   armttbar[146]=1./2.*armttbar[79];
   armttbar[147]=armttbar[145] + armttbar[84] + armttbar[144] + 
   armttbar[146] - 1./2.*armttbar[88] - 5 + armttbar[143];
   armttbar[147]=armttbar[101]*armttbar[147];
   armttbar[148]= - 1./8.*armttbar[8];
   armttbar[110]=armttbar[131] + armttbar[148] + armttbar[110] - 33./4.
   *armttbar[41] + armttbar[36];
   armttbar[110]=armttbar[12]*armttbar[110];
   armttbar[131]=4879./4. - 27*armttbar[76];
   armttbar[131]= - 81./4.*armttbar[87] + 1./2.*armttbar[72] + 1./2.*
   armttbar[74] - 3*armttbar[75] - 27./2.*armttbar[90] + 1./8.*
   armttbar[131] + 33*armttbar[80];
   armttbar[149]= - 1./8.*armttbar[51];
   armttbar[150]= - 27./32.*armttbar[21];
   armttbar[110]=armttbar[110] + armttbar[150] - armttbar[36] + 33./4.*
   armttbar[41] - 185./32.*armttbar[13] + armttbar[149] - armttbar[79]
    + 1./4.*armttbar[131] + 3*armttbar[88];
   armttbar[110]=armttbar[99]*armttbar[110];
   armttbar[131]=1./2.*armttbar[87] + 7./8. + armttbar[90];
   armttbar[151]=pow(armttbar[98],2);
   armttbar[131]=armttbar[151]*armttbar[131];
   armttbar[152]=armttbar[129] + 3./4. - armttbar[40];
   armttbar[153]=pow(armttbar[3],2);
   armttbar[152]=armttbar[153]*armttbar[152];
   armttbar[131]=1./8.*armttbar[131] + 3*armttbar[152];
   armttbar[131]=MMZ*armttbar[131];
   armttbar[152]=25 - 11./2.*armttbar[90];
   armttbar[154]= - 1./2.*armttbar[5];
   armttbar[155]=armttbar[96] + armttbar[154];
   armttbar[155]=armttbar[98]*armttbar[155];
   armttbar[152]=3./2.*armttbar[155] + 1./2.*armttbar[152] - 
   armttbar[87];
   armttbar[152]=armttbar[98]*armttbar[152];
   armttbar[156]= - armttbar[61] - 7./2.*armttbar[65] - 7./8.*
   armttbar[68] + armttbar[66];
   armttbar[157]= - 1./2.*armttbar[101];
   armttbar[158]= - armttbar[100] + armttbar[157];
   armttbar[158]=armttbar[8]*armttbar[158];
   armttbar[159]= - armttbar[16]*armttbar[151];
   armttbar[160]= - armttbar[15]*armttbar[151];
   armttbar[161]= - armttbar[21]*armttbar[100];
   armttbar[110]=3./2.*armttbar[131] + 1./2.*armttbar[117] + 3./16.*
   armttbar[160] + 3./16.*armttbar[110] + 9./32.*armttbar[159] + 1./8.*
   armttbar[114] + 3./4.*armttbar[158] + 1./8.*armttbar[161] + 1./4.*
   armttbar[147] + 1./4.*armttbar[142] + 1./8.*armttbar[152] - 1./4.*
   armttbar[60] - 1./8.*armttbar[64] + 1./2.*armttbar[156] - 
   armttbar[67];
   armttbar[110]=MMZ*armttbar[110];
   armttbar[117]= - armttbar[42] + armttbar[14];
   armttbar[131]= - 1./2.*armttbar[12];
   armttbar[142]=armttbar[131] - 1 - armttbar[21];
   armttbar[147]=armttbar[16]*armttbar[142];
   armttbar[142]=armttbar[15]*armttbar[142];
   armttbar[152]= - armttbar[21]*armttbar[14];
   armttbar[156]= - armttbar[12]*armttbar[14];
   armttbar[142]=3./2.*armttbar[142] + 3*armttbar[147] + 1./2.*
   armttbar[156] + 3./2.*armttbar[117] + armttbar[152];
   armttbar[142]=armttbar[3]*armttbar[142];
   armttbar[147]= - 55*armttbar[39];
   armttbar[158]=457./6. + armttbar[147];
   armttbar[162]=77*armttbar[41];
   armttbar[158]=1./2.*armttbar[158] + armttbar[162];
   armttbar[158]=1./2.*armttbar[158] + armttbar[38];
   armttbar[163]= - armttbar[14]*armttbar[101];
   armttbar[158]=armttbar[131] - 27./16.*armttbar[21] + 1./4.*
   armttbar[163] + 1./2.*armttbar[158] - armttbar[36];
   armttbar[158]=armttbar[12]*armttbar[158];
   armttbar[164]= - 11*armttbar[39];
   armttbar[165]=121./2.*armttbar[41] + 131./6. + armttbar[164];
   armttbar[166]= - armttbar[14]*armttbar[100];
   armttbar[167]=1./2.*armttbar[166];
   armttbar[165]= - 11./8.*armttbar[21] + armttbar[167] + 1./2.*
   armttbar[165] - armttbar[36];
   armttbar[165]=armttbar[21]*armttbar[165];
   armttbar[168]= - 3./2.*armttbar[76] - 229./8. - 11*armttbar[83];
   armttbar[168]= - 13./12.*armttbar[89] - 33./4.*armttbar[80] + 3./2.*
   armttbar[168] - 5*armttbar[43];
   armttbar[169]=3*armttbar[19] - armttbar[49] + armttbar[115];
   armttbar[169]=armttbar[100]*armttbar[169];
   armttbar[170]= - armttbar[48] + armttbar[115];
   armttbar[171]=armttbar[170] + 3*armttbar[18];
   armttbar[171]=armttbar[101]*armttbar[171];
   armttbar[172]=armttbar[100] + 1./2.*armttbar[101];
   armttbar[173]=armttbar[14]*armttbar[172];
   armttbar[174]=armttbar[21]*armttbar[100];
   armttbar[175]=armttbar[100] + armttbar[174];
   armttbar[175]=armttbar[16]*armttbar[175];
   armttbar[176]=armttbar[12]*armttbar[101];
   armttbar[177]=armttbar[101] + armttbar[176];
   armttbar[177]=armttbar[15]*armttbar[177];
   armttbar[178]=armttbar[17]*armttbar[3];
   armttbar[142]=3./4.*armttbar[178] + 1./2.*armttbar[142] + 1./2.*
   armttbar[177] + armttbar[175] + armttbar[158] + armttbar[165] + 3./2.
   *armttbar[173] + 1./2.*armttbar[171] + armttbar[169] + 63./8.*
   armttbar[13] - 5./2.*armttbar[24] + 39./4.*armttbar[23] - 1./2.*
   armttbar[84] + 3*armttbar[45] + armttbar[146] - 13./48.*armttbar[88]
    - armttbar[85] + 61./4.*armttbar[87] - 5./8.*armttbar[75] + 15*
   armttbar[90] - 1./8.*armttbar[141] + 1./2.*armttbar[168] + 
   armttbar[82];
   armttbar[146]=3*armttbar[68];
   armttbar[158]=armttbar[146] - armttbar[66];
   armttbar[165]=3*armttbar[65];
   armttbar[158]=1./2.*armttbar[158] + armttbar[165];
   armttbar[168]= - armttbar[82] - 5 + armttbar[89];
   armttbar[169]=1./2.*armttbar[23];
   armttbar[168]=armttbar[169] + 1./2.*armttbar[168] - armttbar[85];
   armttbar[168]=1./2.*armttbar[100]*armttbar[168];
   armttbar[171]= - armttbar[79] - 5 + armttbar[88];
   armttbar[173]=1./2.*armttbar[13];
   armttbar[171]=armttbar[173] + 1./2.*armttbar[171] - armttbar[84];
   armttbar[171]=armttbar[101]*armttbar[171];
   armttbar[175]=1./4.*armttbar[174];
   armttbar[179]= - armttbar[21] + armttbar[131];
   armttbar[179]=armttbar[3]*armttbar[179];
   armttbar[158]=3./8.*armttbar[179] + 1./8.*armttbar[176] + 
   armttbar[175] + 1./4.*armttbar[171] + armttbar[168] + 1./2.*
   armttbar[158] + armttbar[67];
   armttbar[158]=MMZ*armttbar[158];
   armttbar[142]=1./4.*armttbar[142] + armttbar[158];
   armttbar[142]=MMZ*armttbar[142];
   armttbar[119]=191./2. + armttbar[119];
   armttbar[158]=armttbar[16]*armttbar[100];
   armttbar[179]=11*armttbar[41];
   armttbar[180]= - 1./2.*armttbar[38];
   armttbar[119]=1./2.*armttbar[158] + 61./32.*armttbar[12] + 109./24.*
   armttbar[21] + armttbar[167] + armttbar[180] + 1./8.*armttbar[119]
    + armttbar[179];
   armttbar[119]=armttbar[16]*armttbar[119];
   armttbar[147]= - 233./2. + armttbar[147];
   armttbar[147]=1./2.*armttbar[147] + armttbar[162];
   armttbar[147]=1./2.*armttbar[147] + armttbar[38];
   armttbar[158]= - 1./16.*armttbar[21];
   armttbar[162]=armttbar[15]*armttbar[101];
   armttbar[147]=1./2.*armttbar[162] - 67./48.*armttbar[12] + 
   armttbar[158] + 1./2.*armttbar[163] + 1./2.*armttbar[147] - 
   armttbar[36];
   armttbar[147]=armttbar[15]*armttbar[147];
   armttbar[167]=1./2.*armttbar[48];
   armttbar[181]=armttbar[21]*armttbar[14];
   armttbar[182]=1./8.*armttbar[181];
   armttbar[183]=armttbar[12]*armttbar[14];
   armttbar[184]=1./4.*armttbar[183] + armttbar[182] + 5./8.*
   armttbar[14] + 281./24.*armttbar[20] + 37./2.*armttbar[18] + 413./8.
   *armttbar[19] + 7./16.*armttbar[42] - 3./2.*armttbar[5] - 1./16.*
   armttbar[33] - 17./4.*armttbar[34] + armttbar[167] - 7./12.*
   armttbar[92] - 5./8.*armttbar[93] + 93./4.*armttbar[96] - 57./8.*
   armttbar[94] - 171./8.*armttbar[50] - 7./6.*armttbar[95] + 
   armttbar[49];
   armttbar[185]= - armttbar[16] - 1./4.*armttbar[15];
   armttbar[185]=armttbar[15]*armttbar[185];
   armttbar[186]=pow(armttbar[16],2);
   armttbar[185]= - armttbar[186] + armttbar[185];
   armttbar[185]=armttbar[3]*armttbar[185];
   armttbar[187]= - 17 + 5*armttbar[39];
   armttbar[187]=1./2.*armttbar[187] - 7*armttbar[41];
   armttbar[187]=11./2.*armttbar[187] - armttbar[38];
   armttbar[188]=1./2.*armttbar[15];
   armttbar[189]=armttbar[16] + armttbar[188];
   armttbar[190]=armttbar[3]*armttbar[189];
   armttbar[187]=3./4.*armttbar[190] + 11./16.*armttbar[12] - 45./16.*
   armttbar[21] + 1./2.*armttbar[187] + armttbar[36];
   armttbar[187]=armttbar[17]*armttbar[187];
   armttbar[119]=1./2.*armttbar[187] + 3./4.*armttbar[185] + 1./2.*
   armttbar[147] + 1./2.*armttbar[184] + armttbar[119];
   armttbar[147]= - 77./8.*armttbar[141] - 7./2. + armttbar[76];
   armttbar[184]= - 1./4.*armttbar[12];
   armttbar[187]=armttbar[184] - 5 - armttbar[21];
   armttbar[187]=armttbar[12]*armttbar[187];
   armttbar[191]=1./4.*armttbar[75];
   armttbar[192]=pow(armttbar[21],2);
   armttbar[147]=1./4.*armttbar[187] - 1./4.*armttbar[192] - 3./4.*
   armttbar[13] - 5./4.*armttbar[45] + armttbar[87] + armttbar[191] + 1.
   /2.*armttbar[147] + armttbar[90];
   armttbar[147]=MMZ*armttbar[147];
   armttbar[187]=1./2.*armttbar[19];
   armttbar[192]= - 1./2.*armttbar[34];
   armttbar[193]=5./4.*armttbar[18] + armttbar[187] - 3./2.*
   armttbar[33] + armttbar[192] + armttbar[96] - 1./2.*armttbar[93];
   armttbar[194]=1./2.*armttbar[12];
   armttbar[195]=armttbar[194] - 3 - armttbar[21];
   armttbar[195]=armttbar[16]*armttbar[195];
   armttbar[196]=armttbar[194] + 5 - armttbar[21];
   armttbar[196]=armttbar[15]*armttbar[196];
   armttbar[197]= - 21./2.*armttbar[12] + 7 + armttbar[21];
   armttbar[197]=armttbar[17]*armttbar[197];
   armttbar[147]=1./4.*armttbar[147] + 1./16.*armttbar[197] + 1./16.*
   armttbar[196] + 1./8.*armttbar[195] + 1./4.*armttbar[193] + 
   armttbar[20];
   armttbar[147]=MMZ*armttbar[147];
   armttbar[193]= - armttbar[16] - 9./4.*armttbar[15];
   armttbar[193]=armttbar[15]*armttbar[193];
   armttbar[195]= - 49./2.*armttbar[15];
   armttbar[196]=15./4.*armttbar[17] - armttbar[16] + armttbar[195];
   armttbar[196]=armttbar[17]*armttbar[196];
   armttbar[193]=armttbar[196] - armttbar[186] + armttbar[193];
   armttbar[196]=pow(MMZ,4);
   armttbar[197]= - armttbar[6]*armttbar[196]*armttbar[141];
   armttbar[198]=armttbar[141]*pow(MMZ,3);
   armttbar[197]=21*armttbar[198] + 5*armttbar[197];
   armttbar[197]=armttbar[6]*armttbar[197];
   armttbar[147]=1./16.*armttbar[197] + 1./16.*armttbar[193] + 
   armttbar[147];
   armttbar[147]=armttbar[6]*armttbar[147];
   armttbar[193]= - 1./4.*armttbar[79];
   armttbar[197]= - 1./2.*armttbar[82];
   armttbar[198]= - 1./8.*armttbar[13] + armttbar[169] + 3./8.*
   armttbar[84] - 9./8.*armttbar[51] + armttbar[193] + 11./8.*
   armttbar[88] + 3./8.*armttbar[72] + 3./8.*armttbar[74] + 
   armttbar[197] + 5./4.*armttbar[89] - 1 + 3./4.*armttbar[73];
   armttbar[199]=1./8. + armttbar[36];
   armttbar[200]=1./8.*armttbar[8];
   armttbar[199]=1./3.*armttbar[199] + armttbar[200];
   armttbar[199]=armttbar[12]*armttbar[199];
   armttbar[201]=armttbar[36] + 1./4. + armttbar[38];
   armttbar[201]=armttbar[21]*armttbar[201];
   armttbar[202]=armttbar[8]*armttbar[21];
   armttbar[198]=armttbar[199] + 1./4.*armttbar[202] + 1./2.*
   armttbar[198] + 1./3.*armttbar[201];
   armttbar[198]=MMH*armttbar[198];
   armttbar[119]=1./2.*armttbar[147] + 1./4.*armttbar[198] + 1./2.*
   armttbar[119] + armttbar[142];
   armttbar[119]=armttbar[6]*armttbar[119];
   armttbar[142]=3*armttbar[51];
   armttbar[147]= - 3*armttbar[78];
   armttbar[198]=armttbar[135] + armttbar[142] + armttbar[85] + 
   armttbar[134] + armttbar[137] + 59./8. + armttbar[147];
   armttbar[198]=armttbar[100]*armttbar[198];
   armttbar[199]=1./2.*armttbar[37];
   armttbar[201]=armttbar[199] + 5./2. + armttbar[40];
   armttbar[201]=armttbar[8]*armttbar[201];
   armttbar[202]=1./4.*armttbar[21];
   armttbar[201]=armttbar[201] + armttbar[202] + armttbar[129] - 3./2.
    - armttbar[40];
   armttbar[201]=armttbar[3]*armttbar[201];
   armttbar[146]=armttbar[146] - 5./2.*armttbar[66];
   armttbar[146]=1./2.*armttbar[146] + armttbar[165];
   armttbar[165]=1./2.*armttbar[64];
   armttbar[203]=armttbar[51] + 21./8. - armttbar[77];
   armttbar[203]=armttbar[101]*armttbar[203];
   armttbar[204]=3./2.*armttbar[203];
   armttbar[205]=1./2.*armttbar[161];
   armttbar[172]=armttbar[8]*armttbar[172];
   armttbar[146]=3*armttbar[201] + 3*armttbar[172] + armttbar[205] + 
   armttbar[204] + armttbar[198] + armttbar[60] + armttbar[165] + 3./2.
   *armttbar[67] - armttbar[63] + 1./2.*armttbar[146] + armttbar[61];
   armttbar[146]=MMt*armttbar[146];
   armttbar[172]=7*armttbar[36];
   armttbar[198]= - armttbar[39] + armttbar[41];
   armttbar[198]= - armttbar[36] + 33./4.*armttbar[198] + armttbar[38];
   armttbar[201]=armttbar[12]*armttbar[198];
   armttbar[206]= - 7*armttbar[38];
   armttbar[207]=armttbar[39] - armttbar[41];
   armttbar[208]=armttbar[21]*armttbar[198];
   armttbar[201]=armttbar[201] + armttbar[208] + armttbar[172] + 99./4.
   *armttbar[207] + armttbar[206];
   armttbar[209]= - 17./24.*armttbar[12];
   armttbar[207]=armttbar[36] + 33./4.*armttbar[207] - armttbar[38];
   armttbar[210]=armttbar[209] + armttbar[207] - 17./12.*armttbar[21];
   armttbar[210]=armttbar[16]*armttbar[210];
   armttbar[211]=17./12.*armttbar[12];
   armttbar[212]=armttbar[211] + armttbar[207] + 17./6.*armttbar[21];
   armttbar[212]=armttbar[15]*armttbar[212];
   armttbar[213]=armttbar[21]*armttbar[207];
   armttbar[214]=armttbar[12]*armttbar[207];
   armttbar[213]=armttbar[213] + 1./2.*armttbar[214];
   armttbar[213]=MMZ*armttbar[213];
   armttbar[214]=armttbar[38] - armttbar[36];
   armttbar[215]=armttbar[12]*armttbar[214];
   armttbar[216]=armttbar[21]*armttbar[214];
   armttbar[217]=armttbar[216] + 1./2.*armttbar[215];
   armttbar[217]=MMH*armttbar[217];
   armttbar[218]=armttbar[17]*armttbar[198];
   armttbar[210]=1./3.*armttbar[217] + armttbar[213] + 1./2.*
   armttbar[218] + armttbar[210] + 1./2.*armttbar[212];
   armttbar[210]=armttbar[6]*armttbar[210];
   armttbar[212]=99*armttbar[41] + 17 - 99*armttbar[39];
   armttbar[213]=3*armttbar[38];
   armttbar[217]= - 3*armttbar[36];
   armttbar[212]=armttbar[217] + 1./4.*armttbar[212] + armttbar[213];
   armttbar[212]=armttbar[16]*armttbar[212];
   armttbar[219]=99./2.*armttbar[41] - 17 - 99./2.*armttbar[39];
   armttbar[219]=armttbar[217] + 1./2.*armttbar[219] + armttbar[213];
   armttbar[219]=armttbar[15]*armttbar[219];
   armttbar[212]=armttbar[212] + 1./2.*armttbar[219];
   armttbar[212]=armttbar[3]*armttbar[212];
   armttbar[219]=MMZ*armttbar[3]*armttbar[198];
   armttbar[201]=1./2.*armttbar[210] + 3*armttbar[219] + 1./4.*
   armttbar[201] + armttbar[212];
   armttbar[210]=pow(armttbar[4],2);
   armttbar[201]=armttbar[210]*armttbar[201];
   armttbar[212]=armttbar[140] + 5 + armttbar[139];
   armttbar[212]=armttbar[14]*armttbar[212];
   armttbar[219]= - 1./2.*armttbar[48];
   armttbar[220]=armttbar[219] - armttbar[49] + 11*armttbar[50];
   armttbar[220]=9*armttbar[220] + 25./2.*armttbar[42];
   armttbar[221]=1./8.*armttbar[156];
   armttbar[212]=armttbar[221] + 1./8.*armttbar[152] + 1./2.*
   armttbar[212] - 45./2.*armttbar[18] + 1./4.*armttbar[220] - 45*
   armttbar[19];
   armttbar[220]=armttbar[199] + 11./2.*armttbar[39] - 85./4. + 
   armttbar[40];
   armttbar[222]=1./16.*armttbar[21];
   armttbar[118]=1./16.*armttbar[12] + armttbar[222] + 1./2.*
   armttbar[220] + armttbar[118];
   armttbar[118]=armttbar[16]*armttbar[118];
   armttbar[220]=55./4.*armttbar[39];
   armttbar[223]=1./8.*armttbar[12];
   armttbar[224]=armttbar[223] + 1./8.*armttbar[21] + armttbar[36] - 
   armttbar[38] - 121./4.*armttbar[41] + armttbar[199] + armttbar[220]
    - 93./8. + armttbar[40];
   armttbar[224]=armttbar[15]*armttbar[224];
   armttbar[118]=9./8.*armttbar[185] + 3./4.*armttbar[224] + 1./2.*
   armttbar[212] + 3*armttbar[118];
   armttbar[118]=armttbar[3]*armttbar[118];
   armttbar[185]=7./2.*armttbar[98];
   armttbar[212]=armttbar[21]*armttbar[98];
   armttbar[157]=3./4.*armttbar[212] + armttbar[157] + armttbar[185] - 
   armttbar[100];
   armttbar[212]=13 + armttbar[179];
   armttbar[150]=armttbar[150] + 3./4.*armttbar[212] - armttbar[36];
   armttbar[150]=armttbar[99]*armttbar[150];
   armttbar[150]=1./2.*armttbar[157] + armttbar[150];
   armttbar[157]=1./2.*armttbar[40];
   armttbar[212]=1./4.*armttbar[37];
   armttbar[224]=armttbar[212] + 1 + armttbar[157];
   armttbar[224]=armttbar[3]*armttbar[224];
   armttbar[150]=1./2.*armttbar[150] + armttbar[224];
   armttbar[150]=armttbar[17]*armttbar[150];
   armttbar[224]= - 3*armttbar[72] + 119./4. - 3*armttbar[74];
   armttbar[225]= - armttbar[36] + 3./8.*armttbar[8];
   armttbar[225]=armttbar[12]*armttbar[225];
   armttbar[224]=armttbar[225] + armttbar[36] - 5./8.*armttbar[13] + 15.
   /2.*armttbar[84] + 3./8.*armttbar[51] + armttbar[79] + 1./8.*
   armttbar[224] - 5*armttbar[88];
   armttbar[224]=armttbar[99]*armttbar[224];
   armttbar[225]=1./4.*armttbar[64];
   armttbar[226]=1./2.*armttbar[60];
   armttbar[224]=1./4.*armttbar[224] + armttbar[226] + armttbar[61] + 
   armttbar[225];
   armttbar[224]=MMH*armttbar[224];
   armttbar[103]=1./2.*armttbar[201] + 1./4.*armttbar[146] + 
   armttbar[119] + 1./4.*armttbar[224] + armttbar[110] + 3./2.*
   armttbar[150] + 1./2.*armttbar[103] + armttbar[118];
   armttbar[103]=armttbar[210]*armttbar[103];
   armttbar[110]=armttbar[169] + armttbar[142] - armttbar[85] + 
   armttbar[197] + 1./2.*armttbar[89] + 13./2. + armttbar[147];
   armttbar[110]=armttbar[100]*armttbar[110];
   armttbar[118]= - 3./2.*armttbar[37];
   armttbar[119]=armttbar[118] - 1./2. + armttbar[121];
   armttbar[119]=armttbar[16]*armttbar[119];
   armttbar[142]= - 3*armttbar[37];
   armttbar[146]=1 + armttbar[142];
   armttbar[150]=armttbar[15]*armttbar[146];
   armttbar[119]=armttbar[119] + 1./2.*armttbar[150];
   armttbar[119]=armttbar[3]*armttbar[119];
   armttbar[169]=205./2.*armttbar[41] + armttbar[37] + 151./4.*
   armttbar[39] - 877./18. - armttbar[40];
   armttbar[201]=3*armttbar[36];
   armttbar[119]=3./2.*armttbar[119] + 9./8.*armttbar[12] + 
   armttbar[158] + armttbar[201] + 1./2.*armttbar[169] - 2*armttbar[38]
   ;
   armttbar[119]=armttbar[3]*armttbar[119];
   armttbar[169]= - 11./6.*armttbar[88];
   armttbar[224]= - 11./6.*armttbar[13] + 11./3.*armttbar[84] + 
   armttbar[144] + 11./6.*armttbar[79] + armttbar[169] - 19./3. + 
   armttbar[143];
   armttbar[224]=1./2.*armttbar[101]*armttbar[224];
   armttbar[227]= - 148303./12. + 225*armttbar[76];
   armttbar[227]=675./4.*armttbar[87] - 5./2.*armttbar[72] - 5./2.*
   armttbar[74] + 39*armttbar[75] + 531./4.*armttbar[90] + 1./8.*
   armttbar[227] - 297*armttbar[80];
   armttbar[228]=5*armttbar[36];
   armttbar[229]=225./32.*armttbar[21];
   armttbar[227]=armttbar[229] + armttbar[228] - 297./4.*armttbar[41]
    + 4585./96.*armttbar[13] + 5./8.*armttbar[51] + 5*armttbar[79] + 1./
   4.*armttbar[227] - 15*armttbar[88];
   armttbar[230]= - 5./8.*armttbar[36];
   armttbar[231]= - 39./64.*armttbar[12] + 5./64.*armttbar[8] - 225./
   256.*armttbar[21] + armttbar[230] + 2./3. + 297./32.*armttbar[41];
   armttbar[231]=armttbar[12]*armttbar[231];
   armttbar[227]=1./8.*armttbar[227] + armttbar[231];
   armttbar[227]=armttbar[99]*armttbar[227];
   armttbar[231]= - 1./8.*armttbar[87] - 1./2. - armttbar[90];
   armttbar[231]=armttbar[151]*armttbar[231];
   armttbar[232]=armttbar[129] - 3./4. + 2*armttbar[40];
   armttbar[232]=armttbar[153]*armttbar[232];
   armttbar[231]=1./4.*armttbar[231] + 3*armttbar[232];
   armttbar[231]=MMZ*armttbar[231];
   armttbar[232]= - 17*armttbar[87] - 197 + 11*armttbar[90];
   armttbar[233]=1./2.*armttbar[5];
   armttbar[234]= - armttbar[96] + armttbar[233];
   armttbar[234]=armttbar[98]*armttbar[234];
   armttbar[232]=1./6.*armttbar[232] + armttbar[234];
   armttbar[232]=armttbar[98]*armttbar[232];
   armttbar[235]=armttbar[100] - armttbar[101];
   armttbar[235]=armttbar[8]*armttbar[235];
   armttbar[236]=11./12.*armttbar[114];
   armttbar[237]=armttbar[16]*armttbar[151];
   armttbar[238]=armttbar[15]*armttbar[151];
   armttbar[239]=1./2.*armttbar[61];
   armttbar[110]=armttbar[231] + armttbar[119] + 1./16.*armttbar[238]
    + armttbar[227] + 3./32.*armttbar[237] + armttbar[236] + 3./4.*
   armttbar[235] + 1./8.*armttbar[174] + armttbar[224] + 1./4.*
   armttbar[110] + 1./16.*armttbar[232] + 1./6.*armttbar[60] + 1./12.*
   armttbar[64] + armttbar[67] + armttbar[239] + armttbar[65] - 7./6.*
   armttbar[66] + 11./12.*armttbar[68] + 2./3.*armttbar[54] + 1./3.*
   armttbar[69] - 1./3.*armttbar[52] - 4./9.*armttbar[57] - 1./4.*
   armttbar[70];
   armttbar[110]=MMZ*armttbar[110];
   armttbar[119]=armttbar[82] + 5 - armttbar[89];
   armttbar[119]=armttbar[135] + 1./2.*armttbar[119] + armttbar[85];
   armttbar[119]=armttbar[100]*armttbar[119];
   armttbar[174]=armttbar[79] + 5 - armttbar[88];
   armttbar[174]=armttbar[145] + 1./2.*armttbar[174] + armttbar[84];
   armttbar[174]=armttbar[101]*armttbar[174];
   armttbar[227]= - 2./3.*armttbar[54];
   armttbar[231]=armttbar[227] + 2*armttbar[69] + 2./3.*armttbar[52] + 
   4./3.*armttbar[57] + 1./2.*armttbar[70];
   armttbar[232]=3./4.*armttbar[12];
   armttbar[235]=armttbar[21] + armttbar[232];
   armttbar[235]=armttbar[3]*armttbar[235];
   armttbar[119]=1./2.*armttbar[235] + 1./12.*armttbar[114] + 
   armttbar[205] + 1./6.*armttbar[174] + armttbar[119] - 2*armttbar[67]
    - 4*armttbar[65] + 5./6.*armttbar[66] + 1./3.*armttbar[231] - 3./2.
   *armttbar[68];
   armttbar[119]=MMZ*armttbar[119];
   armttbar[205]=armttbar[42] - armttbar[14];
   armttbar[231]=3*armttbar[21];
   armttbar[235]=5./2.*armttbar[12] + 1 + armttbar[231];
   armttbar[235]=armttbar[16]*armttbar[235];
   armttbar[237]=3 + armttbar[194];
   armttbar[237]=armttbar[15]*armttbar[237];
   armttbar[238]=1./4.*armttbar[181];
   armttbar[235]=1./4.*armttbar[237] + 1./4.*armttbar[235] + 1./12.*
   armttbar[183] + 1./3.*armttbar[205] + armttbar[238];
   armttbar[235]=armttbar[3]*armttbar[235];
   armttbar[237]=17./3.*armttbar[38];
   armttbar[240]= - 5./3.*armttbar[36];
   armttbar[241]=151./8.*armttbar[21] + armttbar[104] + armttbar[240]
    + armttbar[237] - 685./4.*armttbar[41] - 367./2. - 35*armttbar[39];
   armttbar[241]=1./8.*armttbar[241] + armttbar[12];
   armttbar[241]=armttbar[12]*armttbar[241];
   armttbar[242]= - 3*armttbar[19];
   armttbar[243]= - 1./2.*armttbar[42];
   armttbar[244]=armttbar[242] + armttbar[49] + armttbar[243];
   armttbar[244]=armttbar[100]*armttbar[244];
   armttbar[245]= - 119*armttbar[41] - 361./3. - 59./2.*armttbar[39];
   armttbar[246]=1./2.*armttbar[120];
   armttbar[245]=11./8.*armttbar[21] + armttbar[246] + 1./3.*
   armttbar[245] + armttbar[38];
   armttbar[245]=armttbar[21]*armttbar[245];
   armttbar[247]= - 25./6.*armttbar[58] + 2393./64. + 2*armttbar[81];
   armttbar[243]=armttbar[48] + armttbar[243];
   armttbar[243]=1./3.*armttbar[243] - armttbar[18];
   armttbar[243]=1./4.*armttbar[101]*armttbar[243];
   armttbar[248]= - armttbar[100] + armttbar[161];
   armttbar[248]=armttbar[16]*armttbar[248];
   armttbar[249]= - armttbar[101] + armttbar[114];
   armttbar[249]=1./12.*armttbar[15]*armttbar[249];
   armttbar[119]=armttbar[119] + 13./8.*armttbar[178] + 1./2.*
   armttbar[235] + armttbar[249] + 1./4.*armttbar[248] + 1./3.*
   armttbar[241] + 1./4.*armttbar[245] + armttbar[138] + armttbar[243]
    + 1./4.*armttbar[244] - 317./144.*armttbar[13] + 7./8.*armttbar[24]
    - 2./3.*armttbar[25] - 47./12.*armttbar[23] + 1./12.*armttbar[84]
    - 11./6.*armttbar[45] - 1./12.*armttbar[79] + 13./288.*armttbar[88]
    + 1./4.*armttbar[85] - 29./6.*armttbar[87] + 5./48.*armttbar[75] - 
   631./96.*armttbar[90] + 101./96.*armttbar[141] - 1./4.*armttbar[82]
    + 13./96.*armttbar[89] + 37./24.*armttbar[80] + 2./3.*armttbar[43]
    + 5./24.*armttbar[76] + 1./3.*armttbar[247] + 7./2.*armttbar[83];
   armttbar[119]=MMZ*armttbar[119];
   armttbar[138]=armttbar[104] - 5./6.*armttbar[36] + 17./6.*
   armttbar[38] - 685./8.*armttbar[41] - 41./3. - 35./2.*armttbar[39];
   armttbar[235]= - armttbar[15]*armttbar[101];
   armttbar[138]=1./3.*armttbar[235] - 655./216.*armttbar[12] + 1./3.*
   armttbar[138] - 75./16.*armttbar[21];
   armttbar[138]=armttbar[15]*armttbar[138];
   armttbar[241]= - 17./3.*armttbar[38];
   armttbar[244]=5./3.*armttbar[36];
   armttbar[245]=7./12.*armttbar[12] - 365./24.*armttbar[21] + 
   armttbar[244] + armttbar[241] + 685./4.*armttbar[41] + 3593./12. + 
   35*armttbar[39];
   armttbar[247]=9*armttbar[16];
   armttbar[248]=armttbar[247] + 5*armttbar[15];
   armttbar[248]=armttbar[3]*armttbar[248];
   armttbar[245]=1./3.*armttbar[245] + 3./2.*armttbar[248];
   armttbar[245]=armttbar[17]*armttbar[245];
   armttbar[248]=65*armttbar[93] - 239*armttbar[96] + 263*armttbar[50]
    + 197./3.*armttbar[94];
   armttbar[248]=1./2.*armttbar[248] + 7./3.*armttbar[92];
   armttbar[248]=1./2.*armttbar[248] - armttbar[48];
   armttbar[248]=35./48.*armttbar[33] + 1./2.*armttbar[248] + 5./3.*
   armttbar[34];
   armttbar[250]= - 1./4.*armttbar[5];
   armttbar[251]= - 13./144.*armttbar[42];
   armttbar[252]= - 1./36.*armttbar[14];
   armttbar[253]=1./12.*armttbar[156];
   armttbar[254]= - 943./2. - 85*armttbar[39];
   armttbar[254]= - 1303./36.*armttbar[12] - 113./6.*armttbar[21] + 1./
   4.*armttbar[254] - 53*armttbar[41];
   armttbar[254]=armttbar[16]*armttbar[254];
   armttbar[255]=armttbar[3]*armttbar[15]*armttbar[189];
   armttbar[138]=1./4.*armttbar[245] + 1./4.*armttbar[255] + 1./2.*
   armttbar[138] + 1./6.*armttbar[254] + armttbar[253] + armttbar[252]
    - 253./72.*armttbar[20] - 149./8.*armttbar[18] - 3803./144.*
   armttbar[19] + armttbar[251] + 1./3.*armttbar[248] + armttbar[250];
   armttbar[245]=1./9. + armttbar[222];
   armttbar[245]=armttbar[21]*armttbar[245];
   armttbar[106]=armttbar[106] + 59./12. + armttbar[21];
   armttbar[106]=armttbar[12]*armttbar[106];
   armttbar[106]=1./12.*armttbar[106] + armttbar[245] + 13./48.*
   armttbar[13] + 59./144.*armttbar[45] - 1./3.*armttbar[87] - 5./48.*
   armttbar[75] - 19./48.*armttbar[90] + 277./576.*armttbar[141] - 1./6.
   *armttbar[76] + 3./4. + 1./9.*armttbar[58];
   armttbar[106]=MMZ*armttbar[106];
   armttbar[245]= - 109./24.*armttbar[20] - 25./16.*armttbar[18] - 5./
   16.*armttbar[19] + 37./24.*armttbar[33] + 5./6.*armttbar[34] - 
   armttbar[96] + 5./8.*armttbar[93];
   armttbar[248]= - 5./6.*armttbar[12];
   armttbar[254]=armttbar[248] + 47./9. + armttbar[21];
   armttbar[254]=armttbar[16]*armttbar[254];
   armttbar[255]= - 5./8.*armttbar[12] - 77./12. + armttbar[21];
   armttbar[255]=armttbar[15]*armttbar[255];
   armttbar[256]=203./8.*armttbar[12] + 119./4. - 11*armttbar[21];
   armttbar[256]=armttbar[17]*armttbar[256];
   armttbar[106]=armttbar[106] + 1./36.*armttbar[256] + 1./12.*
   armttbar[255] + 1./3.*armttbar[245] + 1./16.*armttbar[254];
   armttbar[106]=MMZ*armttbar[106];
   armttbar[245]=pow(MMZ,2);
   armttbar[254]= - armttbar[245]*armttbar[141];
   armttbar[255]=pow(armttbar[17],2);
   armttbar[254]=1./3.*armttbar[255] + 63./32.*armttbar[254];
   armttbar[254]=MMZ*armttbar[254];
   armttbar[256]=armttbar[6]*armttbar[196]*armttbar[141];
   armttbar[254]=armttbar[254] + 5./8.*armttbar[256];
   armttbar[254]=armttbar[6]*armttbar[254];
   armttbar[256]=1./3.*armttbar[16];
   armttbar[257]=armttbar[256] + 3./2.*armttbar[15];
   armttbar[257]=armttbar[15]*armttbar[257];
   armttbar[258]= - 29./9.*armttbar[16] + 11*armttbar[15];
   armttbar[258]=5*armttbar[258] - 11./6.*armttbar[17];
   armttbar[258]=armttbar[17]*armttbar[258];
   armttbar[257]=5*armttbar[257] + armttbar[258];
   armttbar[106]=armttbar[254] + 1./32.*armttbar[257] + armttbar[106];
   armttbar[106]=armttbar[6]*armttbar[106];
   armttbar[254]= - 1./2.*armttbar[74];
   armttbar[257]= - 1./2.*armttbar[72];
   armttbar[258]=1./3.*armttbar[79];
   armttbar[169]=armttbar[258] + armttbar[169] + armttbar[257] + 13./9.
    + armttbar[254];
   armttbar[241]=armttbar[244] + 11./2. + armttbar[241];
   armttbar[244]= - 1./2.*armttbar[8];
   armttbar[241]=1./3.*armttbar[241] + armttbar[244];
   armttbar[241]=armttbar[12]*armttbar[241];
   armttbar[259]=1./8.*armttbar[51];
   armttbar[169]=1./6.*armttbar[241] + 1./24.*armttbar[13] - 1./8.*
   armttbar[84] + armttbar[259] + 1./4.*armttbar[169] + 1./3.*
   armttbar[45];
   armttbar[169]=1./4.*MMH*armttbar[169];
   armttbar[106]=armttbar[106] + armttbar[169] + 1./2.*armttbar[138] + 
   armttbar[119];
   armttbar[106]=armttbar[6]*armttbar[106];
   armttbar[119]=13./2.*armttbar[39];
   armttbar[138]= - 3./2.*armttbar[38];
   armttbar[241]=3./16.*armttbar[21];
   armttbar[260]=11./8.*armttbar[12];
   armttbar[261]=armttbar[260] + armttbar[241] + armttbar[128] + 
   armttbar[138] + 179./8.*armttbar[41] + armttbar[122] + armttbar[119]
    + 79./3. + armttbar[139];
   armttbar[261]=armttbar[15]*armttbar[261];
   armttbar[262]=19*armttbar[39];
   armttbar[263]=759./2. + armttbar[262];
   armttbar[264]=3*armttbar[37];
   armttbar[263]=1./2.*armttbar[263] + armttbar[264];
   armttbar[265]=53*armttbar[41];
   armttbar[263]=19./8.*armttbar[12] + 1./2.*armttbar[263] + 
   armttbar[265];
   armttbar[263]=armttbar[16]*armttbar[263];
   armttbar[266]= - 11*armttbar[50] - 3./2.*armttbar[48];
   armttbar[267]=5./3.*armttbar[42];
   armttbar[268]=5./3. + armttbar[264];
   armttbar[268]=1./8.*armttbar[14]*armttbar[268];
   armttbar[269]=13./24.*armttbar[183];
   armttbar[270]= - 1./2.*armttbar[15];
   armttbar[271]= - armttbar[16] + armttbar[270];
   armttbar[272]=9./8.*armttbar[3]*armttbar[15]*armttbar[271];
   armttbar[261]=armttbar[272] + 1./2.*armttbar[261] + 1./2.*
   armttbar[263] + armttbar[269] + armttbar[268] + 21./2.*armttbar[18]
    + 33./2.*armttbar[19] + 3./4.*armttbar[266] + armttbar[267];
   armttbar[261]=armttbar[3]*armttbar[261];
   armttbar[263]=1./3.*armttbar[36];
   armttbar[266]=armttbar[263] + armttbar[148];
   armttbar[266]=armttbar[12]*armttbar[266];
   armttbar[273]=armttbar[72] - 119./12. + armttbar[74];
   armttbar[274]= - 1./3.*armttbar[79];
   armttbar[275]= - 1./3.*armttbar[36];
   armttbar[149]=armttbar[266] + armttbar[275] + 5./24.*armttbar[13] - 
   5./2.*armttbar[84] + armttbar[149] + armttbar[274] + 1./8.*
   armttbar[273] + 5./3.*armttbar[88];
   armttbar[149]=armttbar[99]*armttbar[149];
   armttbar[266]= - 1./2.*armttbar[64] - armttbar[60];
   armttbar[149]=1./3.*armttbar[266] + 1./2.*armttbar[149];
   armttbar[149]=5./4.*MMH*armttbar[149];
   armttbar[227]= - armttbar[69] + armttbar[227];
   armttbar[266]=1./4.*armttbar[67];
   armttbar[273]=11./4.*armttbar[64];
   armttbar[276]=11./2.*armttbar[60];
   armttbar[277]=1./2.*armttbar[63];
   armttbar[227]=armttbar[276] + armttbar[273] + armttbar[266] + 
   armttbar[277] + 35./8.*armttbar[65] + 2*armttbar[227] + 1./16.*
   armttbar[66];
   armttbar[278]= - 1 - armttbar[37];
   armttbar[279]=1./2.*armttbar[21];
   armttbar[280]=5 + armttbar[264];
   armttbar[280]=armttbar[8]*armttbar[280];
   armttbar[278]=armttbar[280] + 3*armttbar[278] + armttbar[279];
   armttbar[278]=armttbar[3]*armttbar[278];
   armttbar[280]=1./4.*armttbar[278];
   armttbar[203]=3./4.*armttbar[203];
   armttbar[281]=3./4.*armttbar[112];
   armttbar[227]=armttbar[280] + armttbar[281] + 1./3.*armttbar[227] + 
   armttbar[203];
   armttbar[227]=MMt*armttbar[227];
   armttbar[282]=armttbar[185] + 29*armttbar[101];
   armttbar[282]=11./3.*armttbar[114] + 1./3.*armttbar[282] + 
   armttbar[113];
   armttbar[282]=1./4.*armttbar[282];
   armttbar[283]=4133./24. + 297*armttbar[41];
   armttbar[284]= - 5*armttbar[36];
   armttbar[283]= - 225./32.*armttbar[21] + 1./4.*armttbar[283] + 
   armttbar[284];
   armttbar[283]=1./4.*armttbar[283] + 2./3.*armttbar[12];
   armttbar[283]=armttbar[99]*armttbar[283];
   armttbar[283]=armttbar[282] + armttbar[283];
   armttbar[283]=armttbar[15]*armttbar[283];
   armttbar[285]= - 193 - 297./2.*armttbar[41];
   armttbar[229]=armttbar[229] + 1./2.*armttbar[285] + armttbar[228];
   armttbar[229]=armttbar[99]*armttbar[229];
   armttbar[285]= - 3./4.*armttbar[101];
   armttbar[229]=armttbar[285] + armttbar[229];
   armttbar[286]=1 + armttbar[140];
   armttbar[286]=armttbar[3]*armttbar[286];
   armttbar[229]=1./2.*armttbar[229] + armttbar[286];
   armttbar[229]=armttbar[17]*armttbar[229];
   armttbar[287]=armttbar[116] + 7./3.*armttbar[18] + armttbar[111] + 
   11./6.*armttbar[42];
   armttbar[287]=1./4.*armttbar[101]*armttbar[287];
   armttbar[288]=4913./24. - 113*armttbar[39];
   armttbar[265]=1./4.*armttbar[288] + armttbar[265];
   armttbar[265]=1./4.*armttbar[265] + 13./3.*armttbar[21];
   armttbar[265]=armttbar[21]*armttbar[265];
   armttbar[288]= - 337./32.*armttbar[80] + 1./4.*armttbar[43] - 113./
   256.*armttbar[76] + 43./32.*armttbar[83] - 31./12.*armttbar[58] - 
   13183./1024. + 17./3.*armttbar[81];
   armttbar[289]= - 157./192.*armttbar[51];
   armttbar[290]=7./24.*armttbar[234];
   armttbar[291]=armttbar[19]*armttbar[98];
   armttbar[292]=1./16.*armttbar[291];
   armttbar[293]= - armttbar[18]*armttbar[98];
   armttbar[294]=13./72.*armttbar[38];
   armttbar[295]= - 53./36.*armttbar[36];
   armttbar[296]=11./24.*armttbar[163];
   armttbar[297]= - 1 + armttbar[118];
   armttbar[297]=armttbar[8]*armttbar[297];
   armttbar[298]=1./2.*armttbar[297];
   armttbar[299]=397./64.*armttbar[21] + 11./2.*armttbar[104] + 13./3.*
   armttbar[36] + 7./6.*armttbar[38] + 29./4.*armttbar[41] + 35./3. - 
   31./2.*armttbar[39];
   armttbar[299]=79./48.*armttbar[12] + 1./3.*armttbar[299] + 9./16.*
   armttbar[8];
   armttbar[299]=armttbar[12]*armttbar[299];
   armttbar[300]=armttbar[16]*armttbar[98];
   armttbar[301]=3./8.*armttbar[300];
   armttbar[302]= - 27./2.*armttbar[50] - 7*armttbar[94];
   armttbar[302]= - 25*armttbar[92] - 39*armttbar[93] + 11*
   armttbar[302] + 225./4.*armttbar[96];
   armttbar[302]=1887./8.*armttbar[16] + 25./2.*armttbar[14] + 5279./48.
   *armttbar[20] + 10033./96.*armttbar[18] + 451./2.*armttbar[19] - 69./
   8.*armttbar[33] - 225./16.*armttbar[34] + 5./4.*armttbar[32] + 1./2.
   *armttbar[302] + 5*armttbar[48];
   armttbar[302]=armttbar[99]*armttbar[302];
   armttbar[303]= - 1./6.*armttbar[24];
   armttbar[304]= - 9./8.*armttbar[35];
   armttbar[305]= - 1./4.*armttbar[40];
   armttbar[103]=armttbar[103] + armttbar[227] + armttbar[106] + 
   armttbar[149] + armttbar[110] + armttbar[229] + armttbar[261] + 
   armttbar[283] + 1./4.*armttbar[302] + armttbar[301] + 1./4.*
   armttbar[299] + armttbar[298] + 1./6.*armttbar[265] + armttbar[296]
    + armttbar[295] + armttbar[294] + armttbar[287] - 197./48.*
   armttbar[41] + armttbar[37] + 1./32.*armttbar[293] + 131./48.*
   armttbar[39] + 9055./2304.*armttbar[13] + armttbar[305] + 
   armttbar[292] + armttbar[290] + armttbar[303] + 7./12.*armttbar[25]
    - 17./32.*armttbar[23] + armttbar[84] + armttbar[289] + 7./9.*
   armttbar[45] - 11./24.*armttbar[79] - 145./288.*armttbar[88] + 
   armttbar[304] + 301./384.*armttbar[87] - 19./192.*armttbar[72] - 19./
   192.*armttbar[74] + 11./48.*armttbar[75] + 11./12.*armttbar[77] - 79.
   /128.*armttbar[90] + 1./3.*armttbar[288] + 1./8.*armttbar[141];
   armttbar[103]=armttbar[210]*armttbar[103];
   armttbar[106]=17041./6. + 289*armttbar[39];
   armttbar[110]= - 17*armttbar[41];
   armttbar[106]=1./3.*armttbar[106] + armttbar[110];
   armttbar[227]= - 17*armttbar[36];
   armttbar[106]= - 239./27.*armttbar[12] + 17./6.*armttbar[163] + 1./
   12.*armttbar[106] + armttbar[227];
   armttbar[106]=armttbar[12]*armttbar[106];
   armttbar[229]= - armttbar[3]*armttbar[12];
   armttbar[171]=17./4.*armttbar[229] + 17./2.*armttbar[176] - 257./9.*
   armttbar[66] + 17*armttbar[171];
   armttbar[171]=MMZ*armttbar[171];
   armttbar[117]=armttbar[117] + armttbar[156];
   armttbar[229]= - 1./3. + armttbar[131];
   armttbar[229]=armttbar[15]*armttbar[229];
   armttbar[117]=1./3.*armttbar[117] + armttbar[229];
   armttbar[117]=armttbar[3]*armttbar[117];
   armttbar[229]= - 637./9.*armttbar[75] - 7027./12. + 17*armttbar[80];
   armttbar[229]= - 17./2.*armttbar[84] + 721./9.*armttbar[45] + 17./2.
   *armttbar[79] - 221./48.*armttbar[88] + 1./8.*armttbar[229] + 
   armttbar[87];
   armttbar[229]=1./4.*armttbar[229] - 4./9.*armttbar[13];
   armttbar[170]=1./3.*armttbar[170] + armttbar[18];
   armttbar[170]=armttbar[101]*armttbar[170];
   armttbar[106]=1./12.*armttbar[171] + 17./48.*armttbar[178] + 17./16.
   *armttbar[117] + 17./24.*armttbar[177] + 1./8.*armttbar[106] + 17./
   16.*armttbar[104] + 1./3.*armttbar[229] + 17./8.*armttbar[170];
   armttbar[106]=MMZ*armttbar[106];
   armttbar[110]=armttbar[110] + 2507./2. + 289./3.*armttbar[39];
   armttbar[110]=17./3.*armttbar[162] + 221./72.*armttbar[12] + 17./3.*
   armttbar[163] + 1./12.*armttbar[110] + armttbar[227];
   armttbar[110]=armttbar[15]*armttbar[110];
   armttbar[117]=armttbar[50] - armttbar[94];
   armttbar[117]=17./8.*armttbar[183] + 17./24.*armttbar[14] + 6065./
   144.*armttbar[20] + 977./12.*armttbar[18] + armttbar[19] + 221./96.*
   armttbar[42] + armttbar[154] - 2297./288.*armttbar[33] + 17./4.*
   armttbar[48] - 119./24.*armttbar[92] - 967./48.*armttbar[93] + 17./
   16.*armttbar[117] + armttbar[96];
   armttbar[162]= - 839./3. - 289./8.*armttbar[39];
   armttbar[162]=1./3.*armttbar[162] + 17./8.*armttbar[41];
   armttbar[170]=armttbar[3]*armttbar[15];
   armttbar[162]=17./8.*armttbar[170] + 1475./432.*armttbar[12] + 1./3.
   *armttbar[162] + 17./2.*armttbar[36];
   armttbar[162]=armttbar[17]*armttbar[162];
   armttbar[171]= - 1./2.*armttbar[39];
   armttbar[177]=1253./144.*armttbar[12] - 5./3. + armttbar[171];
   armttbar[177]=armttbar[16]*armttbar[177];
   armttbar[229]=pow(armttbar[15],2);
   armttbar[261]= - armttbar[3]*armttbar[229];
   armttbar[110]=1./2.*armttbar[162] + 17./16.*armttbar[261] + 1./4.*
   armttbar[110] + 1./3.*armttbar[117] + 1./2.*armttbar[177];
   armttbar[117]= - 17./3.*armttbar[79] + 187./6.*armttbar[88] + 17./2.
   *armttbar[72] - 41./9. + 17./2.*armttbar[74];
   armttbar[162]= - 17*armttbar[38];
   armttbar[177]= - 7./4. + armttbar[162];
   armttbar[177]=17./4.*armttbar[8] + 1./3.*armttbar[177] + 17*
   armttbar[36];
   armttbar[177]=armttbar[12]*armttbar[177];
   armttbar[265]= - 1./3.*armttbar[45];
   armttbar[117]=1./6.*armttbar[177] - 17./48.*armttbar[13] + 17./16.*
   armttbar[84] - 17./16.*armttbar[51] + 1./8.*armttbar[117] + 
   armttbar[265];
   armttbar[117]=MMH*armttbar[117];
   armttbar[177]=3211./3. + 139./2.*armttbar[12];
   armttbar[177]=armttbar[15]*armttbar[177];
   armttbar[283]=1607 - 5557./2.*armttbar[12];
   armttbar[283]=armttbar[17]*armttbar[283];
   armttbar[177]=1./6.*armttbar[283] + 1./2.*armttbar[177] + 1927./3.*
   armttbar[20] + 535./2.*armttbar[18] - 107*armttbar[93] - 803./3.*
   armttbar[33];
   armttbar[283]= - 107*armttbar[13] - 1285./3.*armttbar[45] - 3211./3.
    + 107*armttbar[75];
   armttbar[288]= - 257./9. - 5./4.*armttbar[12];
   armttbar[288]=armttbar[12]*armttbar[288];
   armttbar[283]=1./3.*armttbar[283] + 5*armttbar[288];
   armttbar[283]=MMZ*armttbar[283];
   armttbar[177]=1./3.*armttbar[177] + 1./2.*armttbar[283];
   armttbar[177]=MMZ*armttbar[177];
   armttbar[283]= - 13265./3.*armttbar[15] + 1381./2.*armttbar[17];
   armttbar[283]=armttbar[17]*armttbar[283];
   armttbar[283]= - 931./2.*armttbar[229] + armttbar[283];
   armttbar[177]=1./12.*armttbar[283] + armttbar[177];
   armttbar[177]=armttbar[6]*armttbar[177];
   armttbar[106]=1./48.*armttbar[177] + 1./4.*armttbar[117] + 1./2.*
   armttbar[110] + armttbar[106];
   armttbar[106]=armttbar[6]*armttbar[106];
   armttbar[110]= - 17*armttbar[39];
   armttbar[117]=7./2. + armttbar[110];
   armttbar[117]=armttbar[41] + 1./2.*armttbar[117] + armttbar[264];
   armttbar[117]= - 7./48.*armttbar[12] + 1./4.*armttbar[117] + 
   armttbar[201];
   armttbar[117]=armttbar[15]*armttbar[117];
   armttbar[177]= - armttbar[50] - 9./2.*armttbar[48];
   armttbar[283]=95./9. + armttbar[264];
   armttbar[283]=armttbar[14]*armttbar[283];
   armttbar[288]= - 65./6. + armttbar[39];
   armttbar[288]=armttbar[16]*armttbar[288];
   armttbar[117]=9./16.*armttbar[261] + armttbar[117] + 1./4.*
   armttbar[288] + 41./72.*armttbar[156] + 1./8.*armttbar[283] + 5./2.*
   armttbar[18] + armttbar[187] + 1./4.*armttbar[177] + 5./9.*
   armttbar[42];
   armttbar[117]=armttbar[3]*armttbar[117];
   armttbar[177]= - 1./32.*armttbar[21];
   armttbar[148]=25./24.*armttbar[12] + armttbar[148] + armttbar[177]
    + 1./12.*armttbar[41] + armttbar[36];
   armttbar[148]=armttbar[12]*armttbar[148];
   armttbar[187]=993./4. + 1./3.*armttbar[76];
   armttbar[187]=1./4.*armttbar[87] + 1./6.*armttbar[72] + 1./6.*
   armttbar[74] - 25./9.*armttbar[75] + 1./8.*armttbar[187] - 1./9.*
   armttbar[80];
   armttbar[148]=1./3.*armttbar[148] + 1./96.*armttbar[21] + 
   armttbar[275] - 1./36.*armttbar[41] + 289./288.*armttbar[13] - 1./24.
   *armttbar[51] + armttbar[274] + 1./4.*armttbar[187] + armttbar[88];
   armttbar[148]=armttbar[99]*armttbar[148];
   armttbar[187]= - 11 - 17./2.*armttbar[39];
   armttbar[261]=armttbar[3]*armttbar[150];
   armttbar[172]=3./2.*armttbar[261] - 7./36.*armttbar[12] + 
   armttbar[172] + 7./12.*armttbar[41] + 1./3.*armttbar[187] + 
   armttbar[37];
   armttbar[172]=armttbar[3]*armttbar[172];
   armttbar[143]=7./18.*armttbar[13] - 7./9.*armttbar[84] + 
   armttbar[144] - 7./18.*armttbar[79] + 7./18.*armttbar[88] - 37./9.
    + armttbar[143];
   armttbar[143]=armttbar[101]*armttbar[143];
   armttbar[144]=1./2. - armttbar[37];
   armttbar[144]=armttbar[153]*armttbar[144];
   armttbar[187]=MMZ*armttbar[144];
   armttbar[261]= - 1 - armttbar[87];
   armttbar[274]=armttbar[98]*armttbar[261];
   armttbar[274]=7./32.*armttbar[274] - 17./4.*armttbar[60] + 11*
   armttbar[66] - 17./8.*armttbar[64];
   armttbar[283]= - armttbar[8]*armttbar[101];
   armttbar[143]=3./4.*armttbar[187] + 1./4.*armttbar[172] + 25./16.*
   armttbar[148] + 7./72.*armttbar[176] + 3./8.*armttbar[283] + 1./9.*
   armttbar[274] + 1./4.*armttbar[143];
   armttbar[143]=MMZ*armttbar[143];
   armttbar[148]=armttbar[275] + armttbar[200];
   armttbar[148]=armttbar[12]*armttbar[148];
   armttbar[172]= - armttbar[72] + 119./12. - armttbar[74];
   armttbar[148]=armttbar[148] + armttbar[263] - 5./24.*armttbar[13] + 
   5./2.*armttbar[84] + armttbar[259] + armttbar[258] + 1./8.*
   armttbar[172] - 5./3.*armttbar[88];
   armttbar[148]=armttbar[99]*armttbar[148];
   armttbar[165]=armttbar[165] + armttbar[60];
   armttbar[148]=1./3.*armttbar[165] + 1./2.*armttbar[148];
   armttbar[148]=MMH*armttbar[148];
   armttbar[165]= - 5209./24. + armttbar[41];
   armttbar[165]=armttbar[177] + 1./12.*armttbar[165] + armttbar[36];
   armttbar[165]=armttbar[99]*armttbar[165];
   armttbar[172]=armttbar[185] + 47*armttbar[101];
   armttbar[113]=25./3.*armttbar[165] + 7./9.*armttbar[176] + 1./9.*
   armttbar[172] + armttbar[113];
   armttbar[113]=armttbar[15]*armttbar[113];
   armttbar[165]= - armttbar[68] - 437./54.*armttbar[66];
   armttbar[172]= - 1./6.*armttbar[67];
   armttbar[165]= - 7./3.*armttbar[60] - 7./6.*armttbar[64] + 
   armttbar[172] + 1./4.*armttbar[165] - 17./3.*armttbar[63];
   armttbar[112]=1./2.*armttbar[278] + 3./2.*armttbar[112] + 1./3.*
   armttbar[165] + armttbar[204];
   armttbar[112]=MMt*armttbar[112];
   armttbar[165]=armttbar[39] + armttbar[12];
   armttbar[165]=armttbar[15]*armttbar[165];
   armttbar[176]=armttbar[12]*armttbar[39];
   armttbar[185]=MMZ*armttbar[176];
   armttbar[187]= - armttbar[17]*armttbar[39];
   armttbar[200]= - armttbar[16]*armttbar[12];
   armttbar[165]=armttbar[185] + armttbar[187] + armttbar[200] + 
   armttbar[165];
   armttbar[165]=armttbar[6]*armttbar[165];
   armttbar[176]=7./27.*armttbar[176] + armttbar[261] + 1./27.*
   armttbar[39];
   armttbar[185]= - 1./3. + armttbar[171];
   armttbar[185]=armttbar[15]*armttbar[185];
   armttbar[185]=armttbar[256] + armttbar[185];
   armttbar[185]=armttbar[3]*armttbar[185];
   armttbar[204]= - MMZ*armttbar[3]*armttbar[39];
   armttbar[165]=17./108.*armttbar[165] + 1./3.*armttbar[204] + 1./4.*
   armttbar[176] + armttbar[185];
   armttbar[176]=pow(armttbar[2],2);
   armttbar[165]=armttbar[176]*armttbar[165];
   armttbar[185]= - 1./2.*armttbar[50] + armttbar[94];
   armttbar[204]=1./4.*armttbar[96];
   armttbar[123]=armttbar[123] + 25./3.*armttbar[93] + 1./3.*
   armttbar[185] + armttbar[204];
   armttbar[123]= - 11./24.*armttbar[16] + armttbar[124] - 1675./144.*
   armttbar[20] - 2869./288.*armttbar[18] - 1./6.*armttbar[19] + 193./
   72.*armttbar[33] - 1./16.*armttbar[34] + armttbar[125] + 1./2.*
   armttbar[123] - armttbar[48];
   armttbar[123]=armttbar[99]*armttbar[123];
   armttbar[124]= - armttbar[90] - 41./18.*armttbar[80] + 331./144.*
   armttbar[76] - 56965./1728. + armttbar[83];
   armttbar[185]= - 9./2.*armttbar[35];
   armttbar[124]=7./18.*armttbar[79] - 595./216.*armttbar[88] + 
   armttbar[185] - 967./288.*armttbar[87] + 143./144.*armttbar[72] + 
   143./144.*armttbar[74] - 191./324.*armttbar[75] + 1./4.*
   armttbar[124] - 7./9.*armttbar[77];
   armttbar[124]= - 31./576.*armttbar[51] + 1./4.*armttbar[124] - 139./
   81.*armttbar[45];
   armttbar[111]=armttbar[116] + 61./9.*armttbar[18] + armttbar[111] - 
   7./18.*armttbar[42];
   armttbar[111]=armttbar[101]*armttbar[111];
   armttbar[258]= - 3329./2. + 119*armttbar[39];
   armttbar[258]=1./3.*armttbar[258] - 7./2.*armttbar[41];
   armttbar[258]= - 227./96.*armttbar[21] + 7./3.*armttbar[163] + 1./6.
   *armttbar[258] - 7*armttbar[36];
   armttbar[258]= - 235./648.*armttbar[12] + 1./3.*armttbar[258] - 31./
   8.*armttbar[8];
   armttbar[258]=armttbar[12]*armttbar[258];
   armttbar[261]=67./3. - armttbar[41];
   armttbar[261]=1./32.*armttbar[21] + 1./12.*armttbar[261] - 
   armttbar[36];
   armttbar[261]=armttbar[99]*armttbar[261];
   armttbar[261]=armttbar[285] + 25./3.*armttbar[261];
   armttbar[261]=1./2.*armttbar[261] + armttbar[286];
   armttbar[261]=armttbar[17]*armttbar[261];
   armttbar[274]=25./16. - armttbar[39];
   armttbar[274]=armttbar[21]*armttbar[274];
   armttbar[278]=1./6.*armttbar[38];
   armttbar[106]=1./8.*armttbar[165] + 1./4.*armttbar[112] + 1./3.*
   armttbar[106] + 25./24.*armttbar[148] + armttbar[143] + 1./2.*
   armttbar[261] + 1./2.*armttbar[117] + 1./8.*armttbar[113] + 25./24.*
   armttbar[123] + 1./16.*armttbar[300] + 1./16.*armttbar[258] + 1./4.*
   armttbar[297] + 1./96.*armttbar[274] + 7./144.*armttbar[104] - 79./
   144.*armttbar[36] + armttbar[278] + 1./8.*armttbar[111] + 5./1728.*
   armttbar[41] + armttbar[199] + 35./864.*armttbar[39] - 2983./41472.*
   armttbar[13] - 1./8.*armttbar[40] + 1./96.*armttbar[291] + 7./144.*
   armttbar[234] - 1./36.*armttbar[24] - 1./32.*armttbar[23] + 1./2.*
   armttbar[124] - 1./3.*armttbar[84];
   armttbar[106]=armttbar[176]*armttbar[106];
   armttbar[111]=14./27.*armttbar[57] + armttbar[54];
   armttbar[112]=1./8.*armttbar[155] - 13./16.*armttbar[87] - 1 - 11./
   48.*armttbar[90];
   armttbar[112]=armttbar[98]*armttbar[112];
   armttbar[111]=1./2.*armttbar[112] + armttbar[226] + armttbar[225] + 
   armttbar[172] - 7./2.*armttbar[66] + 2*armttbar[111] + 1./4.*
   armttbar[68];
   armttbar[112]=armttbar[16]*armttbar[146];
   armttbar[112]=armttbar[112] + armttbar[150];
   armttbar[112]=armttbar[3]*armttbar[112];
   armttbar[113]= - armttbar[38] + 31./6.*armttbar[41] + armttbar[37]
    + 19./4.*armttbar[39] - 59./18. + armttbar[40];
   armttbar[112]=3./4.*armttbar[112] + 11./24.*armttbar[12] + 
   armttbar[222] + 1./2.*armttbar[113] + armttbar[201];
   armttbar[112]=armttbar[3]*armttbar[112];
   armttbar[113]= - 224603./36. + 7*armttbar[76];
   armttbar[113]=275./3.*armttbar[75] + 5./4.*armttbar[90] + 1./8.*
   armttbar[113] - 37./3.*armttbar[80];
   armttbar[113]=7./4.*armttbar[87] + armttbar[257] + 1./3.*
   armttbar[113] + armttbar[254];
   armttbar[117]=7./96.*armttbar[21];
   armttbar[113]=armttbar[117] + armttbar[36] - 37./36.*armttbar[41] - 
   9091./864.*armttbar[13] + armttbar[259] + armttbar[79] + 1./4.*
   armttbar[113] - 3*armttbar[88];
   armttbar[123]=10./3. + 37./32.*armttbar[41];
   armttbar[123]= - 275./576.*armttbar[12] + 1./64.*armttbar[8] - 7./
   768.*armttbar[21] + 1./9.*armttbar[123] - 1./8.*armttbar[36];
   armttbar[123]=armttbar[12]*armttbar[123];
   armttbar[113]=1./8.*armttbar[113] + armttbar[123];
   armttbar[113]=armttbar[99]*armttbar[113];
   armttbar[123]=1 + armttbar[87];
   armttbar[123]=armttbar[151]*armttbar[123];
   armttbar[123]=1./48.*armttbar[123] + 3*armttbar[144];
   armttbar[123]=MMZ*armttbar[123];
   armttbar[111]=1./2.*armttbar[123] + armttbar[112] + 1./48.*
   armttbar[160] + 5*armttbar[113] + 1./32.*armttbar[159] + 
   armttbar[236] + 3./4.*armttbar[283] + 1./3.*armttbar[111] + 
   armttbar[224];
   armttbar[111]=MMZ*armttbar[111];
   armttbar[112]=46471./216. + armttbar[76];
   armttbar[112]=143./18.*armttbar[80] + 1./2.*armttbar[112] + 1./3.*
   armttbar[43];
   armttbar[113]=1./3.*armttbar[90];
   armttbar[123]=1./4.*armttbar[84];
   armttbar[112]=607./432.*armttbar[13] + 1./8.*armttbar[24] - 68./27.*
   armttbar[25] - 1./12.*armttbar[23] + armttbar[123] - 457./54.*
   armttbar[45] + armttbar[193] + 13./96.*armttbar[88] + 1./12.*
   armttbar[87] - 377./432.*armttbar[75] + 1./4.*armttbar[112] + 
   armttbar[113];
   armttbar[124]=2*armttbar[57] - armttbar[54];
   armttbar[124]=34*armttbar[124] + 263./2.*armttbar[66];
   armttbar[124]=1./9.*armttbar[124] - 1./4.*armttbar[67];
   armttbar[143]=armttbar[3]*armttbar[12];
   armttbar[114]=1./8.*armttbar[143] + 1./4.*armttbar[114] + 1./3.*
   armttbar[124] + 1./2.*armttbar[174];
   armttbar[114]=MMZ*armttbar[114];
   armttbar[124]=armttbar[205] + armttbar[183];
   armttbar[143]= - 1 - 17./6.*armttbar[12];
   armttbar[143]=armttbar[16]*armttbar[143];
   armttbar[144]=1./3. + armttbar[194];
   armttbar[144]=armttbar[15]*armttbar[144];
   armttbar[124]=armttbar[144] + 1./3.*armttbar[124] + armttbar[143];
   armttbar[124]=armttbar[3]*armttbar[124];
   armttbar[143]=17*armttbar[38];
   armttbar[144]=armttbar[284] + armttbar[143] - 821./12.*armttbar[41]
    - 3659./6. - 73*armttbar[39];
   armttbar[144]=1./3.*armttbar[144] + armttbar[104];
   armttbar[144]=1./3.*armttbar[144] - 3./8.*armttbar[21];
   armttbar[144]=1./8.*armttbar[144] + 59./81.*armttbar[12];
   armttbar[144]=armttbar[12]*armttbar[144];
   armttbar[146]= - armttbar[21]*armttbar[39];
   armttbar[112]=1./3.*armttbar[114] + 5./8.*armttbar[178] + 1./8.*
   armttbar[124] + armttbar[249] + armttbar[144] + 1./24.*armttbar[146]
    + 1./8.*armttbar[163] + 1./3.*armttbar[112] + armttbar[243];
   armttbar[112]=MMZ*armttbar[112];
   armttbar[114]= - 5./4.*armttbar[36] + 17./4.*armttbar[38] - 821./48.
   *armttbar[41] - 917./3. - 73./4.*armttbar[39];
   armttbar[105]=1./2.*armttbar[235] - 4013./432.*armttbar[12] - 35./96.
   *armttbar[21] + 1./3.*armttbar[114] + armttbar[105];
   armttbar[105]=armttbar[15]*armttbar[105];
   armttbar[114]=277./12.*armttbar[12] + 19./8.*armttbar[21] + 
   armttbar[228] + armttbar[162] + 821./12.*armttbar[41] + 25201./36.
    + 73*armttbar[39];
   armttbar[124]=17./3.*armttbar[16] + 15*armttbar[15];
   armttbar[124]=armttbar[3]*armttbar[124];
   armttbar[114]=1./9.*armttbar[114] + 1./2.*armttbar[124];
   armttbar[114]=armttbar[17]*armttbar[114];
   armttbar[124]=295*armttbar[50] - 337*armttbar[94];
   armttbar[124]=1993./27.*armttbar[93] + 1./27.*armttbar[124] + 3*
   armttbar[96];
   armttbar[124]=1./2.*armttbar[124] + 7./9.*armttbar[92];
   armttbar[124]=1./2.*armttbar[124] - 1./3.*armttbar[48];
   armttbar[144]=43./6. - 37*armttbar[39];
   armttbar[144]= - 151./18.*armttbar[12] + 1./2.*armttbar[144] + 
   armttbar[21];
   armttbar[144]=armttbar[16]*armttbar[144];
   armttbar[148]= - 13./3.*armttbar[16] + armttbar[188];
   armttbar[148]=armttbar[3]*armttbar[15]*armttbar[148];
   armttbar[105]=1./4.*armttbar[114] + 1./4.*armttbar[148] + 1./3.*
   armttbar[105] + 1./12.*armttbar[144] + armttbar[253] + armttbar[252]
    - 12533./648.*armttbar[20] - 5899./216.*armttbar[18] + 23./48.*
   armttbar[19] + armttbar[251] - 13./36.*armttbar[5] + 547./144.*
   armttbar[33] + 1./2.*armttbar[124] - 1./9.*armttbar[34];
   armttbar[114]=173*armttbar[13] + 1043./3.*armttbar[45] + 4157./3. - 
   173*armttbar[75];
   armttbar[124]=1043./27. + 5./4.*armttbar[12];
   armttbar[124]=armttbar[12]*armttbar[124];
   armttbar[114]=1./9.*armttbar[114] + armttbar[124];
   armttbar[114]=MMZ*armttbar[114];
   armttbar[124]=173*armttbar[93] + 781./3.*armttbar[33];
   armttbar[124]= - 875./9.*armttbar[20] - 865./18.*armttbar[18] + 1./9.
   *armttbar[124] + 11./2.*armttbar[19];
   armttbar[144]= - 11./3. + 9./2.*armttbar[12];
   armttbar[144]=armttbar[16]*armttbar[144];
   armttbar[148]= - 4157./3. - 301./2.*armttbar[12];
   armttbar[148]=armttbar[15]*armttbar[148];
   armttbar[150]= - 529 + 4259./2.*armttbar[12];
   armttbar[150]=armttbar[17]*armttbar[150];
   armttbar[114]=1./6.*armttbar[114] + 1./162.*armttbar[150] + 1./54.*
   armttbar[148] + 1./3.*armttbar[124] + 1./2.*armttbar[144];
   armttbar[114]=MMZ*armttbar[114];
   armttbar[124]=5*armttbar[16];
   armttbar[144]=armttbar[124] + 1429./18.*armttbar[15];
   armttbar[144]=armttbar[15]*armttbar[144];
   armttbar[148]=925./54.*armttbar[17] - 9*armttbar[16] + 11495./81.*
   armttbar[15];
   armttbar[148]=armttbar[17]*armttbar[148];
   armttbar[144]=1./3.*armttbar[144] + armttbar[148];
   armttbar[114]=1./4.*armttbar[144] + armttbar[114];
   armttbar[144]=armttbar[6]*MMZ*armttbar[255];
   armttbar[114]=1./8.*armttbar[114] + 17./27.*armttbar[144];
   armttbar[114]=armttbar[6]*armttbar[114];
   armttbar[105]=armttbar[114] + armttbar[169] + 1./2.*armttbar[105] + 
   armttbar[112];
   armttbar[105]=armttbar[6]*armttbar[105];
   armttbar[112]=armttbar[260] + armttbar[241] + armttbar[128] + 
   armttbar[138] + 35./8.*armttbar[41] + armttbar[122] + armttbar[119]
    + 43./3. + armttbar[139];
   armttbar[112]=armttbar[15]*armttbar[112];
   armttbar[114]=197./6. + armttbar[262];
   armttbar[114]= - 7./12.*armttbar[12] + 1./2.*armttbar[114] + 
   armttbar[264];
   armttbar[114]=armttbar[16]*armttbar[114];
   armttbar[119]= - armttbar[50] + armttbar[219];
   armttbar[112]=armttbar[272] + 1./2.*armttbar[112] + 1./4.*
   armttbar[114] + armttbar[269] + armttbar[268] + 9./2.*armttbar[18]
    + armttbar[126] + 9./4.*armttbar[119] + armttbar[267];
   armttbar[112]=armttbar[3]*armttbar[112];
   armttbar[114]=28*armttbar[54] - 997./16.*armttbar[66];
   armttbar[114]=armttbar[276] + armttbar[273] + armttbar[266] + 
   armttbar[277] + 1./27.*armttbar[114] + 19./8.*armttbar[65];
   armttbar[114]=armttbar[280] + armttbar[281] + 1./3.*armttbar[114] + 
   armttbar[203];
   armttbar[114]=MMt*armttbar[114];
   armttbar[119]=271./8.*armttbar[76] - 253733./864. + 7*armttbar[83];
   armttbar[119]= - 473./72.*armttbar[80] + 1./8.*armttbar[119] + 
   armttbar[43];
   armttbar[119]= - 11./6.*armttbar[79] - 145./72.*armttbar[88] + 
   armttbar[185] + 439./288.*armttbar[87] - 19./48.*armttbar[72] - 19./
   48.*armttbar[74] + 809./324.*armttbar[75] + 11./3.*armttbar[77] - 
   109./96.*armttbar[90] + 1./3.*armttbar[119] + 1./2.*armttbar[141];
   armttbar[126]=5153./8. + 37./3.*armttbar[41];
   armttbar[126]= - 7./96.*armttbar[21] + 1./12.*armttbar[126] - 
   armttbar[36];
   armttbar[126]=1./4.*armttbar[126] + 10./27.*armttbar[12];
   armttbar[126]=armttbar[99]*armttbar[126];
   armttbar[126]=armttbar[282] + 5*armttbar[126];
   armttbar[126]=armttbar[15]*armttbar[126];
   armttbar[128]=13./2.*armttbar[36] + 7./4.*armttbar[38] - 71./24.*
   armttbar[41] - 173./9. - 29./4.*armttbar[39];
   armttbar[104]=269./128.*armttbar[21] + 1./3.*armttbar[128] + 11./4.*
   armttbar[104];
   armttbar[104]=85./2592.*armttbar[12] + 1./3.*armttbar[104] + 9./32.*
   armttbar[8];
   armttbar[104]=armttbar[12]*armttbar[104];
   armttbar[128]= - 37./2.*armttbar[50] + 47*armttbar[94];
   armttbar[128]= - 275./3.*armttbar[93] + 1./3.*armttbar[128] + 7./4.*
   armttbar[96];
   armttbar[128]=1./3.*armttbar[128] - 5*armttbar[92];
   armttbar[138]=1./4.*armttbar[32];
   armttbar[128]= - 437./72.*armttbar[16] + 5./2.*armttbar[14] + 4985./
   144.*armttbar[20] + 29797./864.*armttbar[18] - 19./6.*armttbar[19]
    - 1499./216.*armttbar[33] - 7./48.*armttbar[34] + armttbar[138] + 1.
   /2.*armttbar[128] + armttbar[48];
   armttbar[128]=armttbar[99]*armttbar[128];
   armttbar[144]= - 643./3. - 37./2.*armttbar[41];
   armttbar[117]=armttbar[117] + 1./18.*armttbar[144] + armttbar[36];
   armttbar[117]=armttbar[99]*armttbar[117];
   armttbar[117]=armttbar[285] + 5*armttbar[117];
   armttbar[117]=1./2.*armttbar[117] + armttbar[286];
   armttbar[117]=armttbar[17]*armttbar[117];
   armttbar[144]=49./24. + armttbar[110];
   armttbar[144]=1./8.*armttbar[144] - armttbar[21];
   armttbar[144]=armttbar[21]*armttbar[144];
   armttbar[148]= - 1./18.*armttbar[24];
   armttbar[104]=armttbar[106] + armttbar[114] + armttbar[105] + 
   armttbar[149] + armttbar[111] + armttbar[117] + armttbar[112] + 
   armttbar[126] + 5./4.*armttbar[128] + armttbar[301] + 1./2.*
   armttbar[104] + armttbar[298] + 1./12.*armttbar[144] + armttbar[296]
    + armttbar[295] + armttbar[294] + armttbar[287] + 1./144.*
   armttbar[41] + armttbar[37] + 1./96.*armttbar[136] + 41./144.*
   armttbar[39] + 2461./6912.*armttbar[13] + armttbar[305] + 
   armttbar[292] + armttbar[290] + armttbar[148] - 103./324.*
   armttbar[25] - 5./32.*armttbar[23] + armttbar[84] + armttbar[289] + 
   1./4.*armttbar[119] - 1./9.*armttbar[45];
   armttbar[104]=armttbar[176]*armttbar[104];
   armttbar[105]=3*armttbar[97];
   armttbar[106]=3./4.*armttbar[71];
   armttbar[111]= - 3*armttbar[86];
   armttbar[112]=1./2.*armttbar[77];
   armttbar[114]=1./2.*armttbar[85];
   armttbar[117]=27./4.*armttbar[35];
   armttbar[119]=1./4.*armttbar[79];
   armttbar[126]= - 21./4.*armttbar[51];
   armttbar[128]= - 1./4.*armttbar[13];
   armttbar[136]= - 1./4.*armttbar[84];
   armttbar[144]= - 1./6.*armttbar[38];
   armttbar[149]= - 1./12.*armttbar[36] + armttbar[144] + armttbar[140]
    + armttbar[128] + armttbar[139] + armttbar[135] + armttbar[136] + 
   armttbar[126] + armttbar[119] + armttbar[117] + armttbar[114] + 
   armttbar[112] + armttbar[134] + armttbar[137] + armttbar[78] + 
   armttbar[111] + armttbar[106] - 293./32. + armttbar[105];
   armttbar[150]= - 11./3. - 9*armttbar[35];
   armttbar[150]=armttbar[129] + 1./2.*armttbar[150] - armttbar[40];
   armttbar[155]=1./2.*armttbar[150];
   armttbar[159]=1./4.*armttbar[8];
   armttbar[160]=armttbar[159] + armttbar[155] + armttbar[275];
   armttbar[160]=armttbar[8]*armttbar[160];
   armttbar[162]=1./2.*armttbar[18];
   armttbar[163]= - armttbar[14] - armttbar[20] + armttbar[162] + 
   armttbar[92] + armttbar[219];
   armttbar[163]=armttbar[99]*armttbar[163];
   armttbar[165]=armttbar[38] + armttbar[107];
   armttbar[165]=armttbar[12]*armttbar[165];
   armttbar[169]= - 1 + armttbar[107];
   armttbar[169]=armttbar[15]*armttbar[99]*armttbar[169];
   armttbar[174]= - 1./2. - armttbar[36];
   armttbar[178]=armttbar[17]*armttbar[99]*armttbar[174];
   armttbar[185]=1 + armttbar[84];
   armttbar[193]=MMH*armttbar[99]*armttbar[185];
   armttbar[149]=1./4.*armttbar[193] + 1./2.*armttbar[178] + 1./2.*
   armttbar[169] + 1./2.*armttbar[163] + 1./12.*armttbar[165] + 1./2.*
   armttbar[149] + armttbar[160];
   armttbar[149]=MMH*armttbar[149];
   armttbar[160]= - 1./4.*armttbar[47] - armttbar[91];
   armttbar[165]= - 1./2.*armttbar[32];
   armttbar[160]=armttbar[165] - 3./8.*armttbar[48] + 1./4.*
   armttbar[92] - 1./4.*armttbar[93] - 1./2.*armttbar[94] - 
   armttbar[49] + 3./2.*armttbar[160] + armttbar[95];
   armttbar[203]=25./8.*armttbar[19];
   armttbar[205]=armttbar[116] + 11./8.*armttbar[18] + armttbar[203] + 
   armttbar[160] + 27./8.*armttbar[42];
   armttbar[219]= - 1./8.*armttbar[21];
   armttbar[222]=1./6.*armttbar[36];
   armttbar[224]= - 1./24.*armttbar[12];
   armttbar[225]=armttbar[224] + 7./16.*armttbar[8] + armttbar[219] + 
   armttbar[222] + 1./16. + armttbar[38];
   armttbar[225]=armttbar[16]*armttbar[225];
   armttbar[226]=armttbar[199] + armttbar[40] + 55./3. + 9./2.*
   armttbar[35];
   armttbar[226]=1./8.*armttbar[226] + armttbar[263];
   armttbar[226]=armttbar[14]*armttbar[226];
   armttbar[228]= - armttbar[8]*armttbar[14];
   armttbar[234]=1./4.*armttbar[228];
   armttbar[235]=armttbar[36] + 5./6. + armttbar[38];
   armttbar[235]= - 1./6.*armttbar[12] + 1./2.*armttbar[235] - 5./3.*
   armttbar[8];
   armttbar[235]=armttbar[15]*armttbar[235];
   armttbar[182]=1./2.*armttbar[235] + armttbar[225] + 1./8.*
   armttbar[183] + armttbar[234] + armttbar[182] + 1./2.*armttbar[205]
    + armttbar[226];
   armttbar[205]=1./16.*armttbar[37] + 1./8.*armttbar[40] + 1./3. + 9./
   16.*armttbar[35];
   armttbar[225]= - 7./16.*armttbar[8];
   armttbar[226]=armttbar[225] + armttbar[205] + armttbar[222];
   armttbar[226]=armttbar[17]*armttbar[226];
   armttbar[149]=1./4.*armttbar[149] + 1./2.*armttbar[182] + 
   armttbar[226];
   armttbar[149]=MMH*armttbar[149];
   armttbar[182]= - armttbar[16]*armttbar[38];
   armttbar[180]=armttbar[180] - armttbar[36];
   armttbar[226]=armttbar[15]*armttbar[180];
   armttbar[235]=1./2.*armttbar[38];
   armttbar[236]=armttbar[235] + armttbar[36];
   armttbar[236]=armttbar[17]*armttbar[236];
   armttbar[226]=1./3.*armttbar[236] + armttbar[182] + 1./3.*
   armttbar[226];
   armttbar[226]=MMH*armttbar[226];
   armttbar[236]= - armttbar[14] + armttbar[256];
   armttbar[241]=1./3.*armttbar[15];
   armttbar[236]=1./2.*armttbar[236] + armttbar[241];
   armttbar[236]=armttbar[15]*armttbar[236];
   armttbar[243]= - 1./3.*armttbar[16];
   armttbar[249]=armttbar[14] + armttbar[243];
   armttbar[251]= - 1./3.*armttbar[15];
   armttbar[249]=1./2.*armttbar[249] + armttbar[251];
   armttbar[249]=armttbar[17]*armttbar[249];
   armttbar[252]= - armttbar[14] + armttbar[16];
   armttbar[253]=armttbar[16]*armttbar[252];
   armttbar[226]=armttbar[226] + armttbar[249] + armttbar[253] + 
   armttbar[236];
   armttbar[226]=armttbar[6]*MMH*armttbar[226];
   armttbar[236]=pow(armttbar[14],2);
   armttbar[249]= - 3./2.*armttbar[236];
   armttbar[254]=11*armttbar[14] - 137./2.*armttbar[16];
   armttbar[254]=armttbar[16]*armttbar[254];
   armttbar[257]=107*armttbar[14] + 151*armttbar[16];
   armttbar[257]=1./3.*armttbar[257] - 11./2.*armttbar[15];
   armttbar[257]=armttbar[15]*armttbar[257];
   armttbar[254]=armttbar[257] + armttbar[249] + armttbar[254];
   armttbar[257]=11./3.*armttbar[15] - armttbar[14] - 17./2.*
   armttbar[16];
   armttbar[257]=1./4.*armttbar[257] + armttbar[17];
   armttbar[257]=armttbar[17]*armttbar[257];
   armttbar[254]=1./16.*armttbar[254] + armttbar[257];
   armttbar[149]=1./8.*armttbar[226] + 1./2.*armttbar[254] + 
   armttbar[149];
   armttbar[226]=101./64. - armttbar[71];
   armttbar[226]=3./64.*armttbar[141] - 5./16.*armttbar[73] + 1./2.*
   armttbar[226] + armttbar[86];
   armttbar[254]=15 + 11./4.*armttbar[21];
   armttbar[257]=3*armttbar[8];
   armttbar[254]=1./2.*armttbar[254] + armttbar[257];
   armttbar[254]=armttbar[8]*armttbar[254];
   armttbar[258]=armttbar[181] + 25*armttbar[42] - 9*armttbar[14];
   armttbar[259]=armttbar[8]*armttbar[14];
   armttbar[258]=1./8.*armttbar[258] + armttbar[259];
   armttbar[258]=armttbar[3]*armttbar[258];
   armttbar[260]=3./2.*armttbar[21];
   armttbar[261]=1 + armttbar[260];
   armttbar[261]=armttbar[21]*armttbar[261];
   armttbar[239]=armttbar[239] + 3*armttbar[59] - armttbar[62];
   armttbar[239]=MMH*armttbar[239];
   armttbar[262]=MMt*armttbar[62];
   armttbar[266]=1./2.*armttbar[262] + 1./4.*armttbar[239] + 1./2.*
   armttbar[258] + 1./4.*armttbar[254] + 1./64.*armttbar[261] + 13./32.
   *armttbar[51] + 1./16.*armttbar[85] + 1./2.*armttbar[226] + 
   armttbar[44];
   armttbar[266]=MMt*armttbar[266];
   armttbar[267]= - 5*armttbar[91] - armttbar[95];
   armttbar[268]=1./4.*armttbar[93];
   armttbar[192]=7./4.*armttbar[42] + armttbar[192] - 3*armttbar[32] - 
   3./2.*armttbar[92] + armttbar[268] + armttbar[204] + 3./4.*
   armttbar[267] - armttbar[94];
   armttbar[204]=25./2.*armttbar[20];
   armttbar[267]=35./2.*armttbar[14];
   armttbar[269]=armttbar[238] + armttbar[267] + armttbar[204] - 7*
   armttbar[18] + armttbar[192] - 21./2.*armttbar[19];
   armttbar[272]=11 - 5*armttbar[21];
   armttbar[273]= - 7*armttbar[8];
   armttbar[272]=9./4.*armttbar[272] + armttbar[273];
   armttbar[272]=armttbar[16]*armttbar[272];
   armttbar[269]=1./2.*armttbar[272] + 1./2.*armttbar[269] + 
   armttbar[259];
   armttbar[272]=47./2. - 31./3.*armttbar[21];
   armttbar[272]=1./16.*armttbar[272] + 5./3.*armttbar[8];
   armttbar[272]=armttbar[15]*armttbar[272];
   armttbar[269]=1./4.*armttbar[269] + armttbar[272];
   armttbar[272]= - 31./2. - 9*armttbar[97];
   armttbar[272]=armttbar[72] + armttbar[74] - armttbar[75] - 5./2.*
   armttbar[77] - 3./2.*armttbar[82] + armttbar[137] + armttbar[147] + 
   5./4.*armttbar[73] + 1./2.*armttbar[272] - 5*armttbar[86];
   armttbar[136]= - 1./8.*armttbar[37] + armttbar[305] + 3./8.*
   armttbar[23] + armttbar[136] + 31./16.*armttbar[51] + armttbar[304]
    - 1./4.*armttbar[85] + 1./4.*armttbar[272] - armttbar[44];
   armttbar[272]= - 1./4. + armttbar[36];
   armttbar[272]=armttbar[21]*armttbar[272];
   armttbar[272]=armttbar[136] + 1./6.*armttbar[272];
   armttbar[274]=9*armttbar[35];
   armttbar[276]=11./3. + armttbar[274];
   armttbar[276]=armttbar[199] + 1./2.*armttbar[276] + armttbar[40];
   armttbar[280]=1./8.*armttbar[276];
   armttbar[281]=armttbar[244] + armttbar[280] + armttbar[263];
   armttbar[281]=armttbar[8]*armttbar[281];
   armttbar[282]=armttbar[8] + armttbar[131];
   armttbar[282]=armttbar[12]*armttbar[282];
   armttbar[284]=1./8.*armttbar[282];
   armttbar[285]=armttbar[64] - 9*armttbar[59] - armttbar[62];
   armttbar[285]=1./2.*armttbar[285] + armttbar[60];
   armttbar[285]=MMH*armttbar[285];
   armttbar[286]=1./8.*armttbar[285];
   armttbar[272]=armttbar[286] + armttbar[284] + 1./2.*armttbar[272] + 
   armttbar[281];
   armttbar[272]=MMH*armttbar[272];
   armttbar[257]=armttbar[257] + 7 + armttbar[260];
   armttbar[260]= - armttbar[3]*armttbar[14];
   armttbar[281]=1./8.*armttbar[257] + armttbar[260];
   armttbar[281]=armttbar[17]*armttbar[281];
   armttbar[269]=armttbar[266] + 1./2.*armttbar[272] + 1./2.*
   armttbar[269] + armttbar[281];
   armttbar[269]=MMt*armttbar[269];
   armttbar[272]=1./3.*armttbar[38];
   armttbar[224]=armttbar[224] - 17./24.*armttbar[8] - 1./24.*
   armttbar[21] + armttbar[275] + 1./8. + armttbar[272];
   armttbar[224]=armttbar[16]*armttbar[224];
   armttbar[287]=1./3.*armttbar[216];
   armttbar[288]= - armttbar[38] + armttbar[36];
   armttbar[289]=armttbar[288] + armttbar[287];
   armttbar[290]=armttbar[8]*armttbar[214];
   armttbar[291]=1./3.*armttbar[290];
   armttbar[215]=1./12.*armttbar[215] + 1./4.*armttbar[289] + 
   armttbar[291];
   armttbar[215]=MMH*armttbar[215];
   armttbar[289]=armttbar[14]*armttbar[288];
   armttbar[292]= - 1./2. + armttbar[38];
   armttbar[293]=1./6.*armttbar[12] + 17./6.*armttbar[8] + 1./6.*
   armttbar[21] + armttbar[292] - armttbar[36];
   armttbar[293]=armttbar[15]*armttbar[293];
   armttbar[294]=armttbar[17]*armttbar[288];
   armttbar[215]=1./2.*armttbar[215] + 1./3.*armttbar[294] + 1./4.*
   armttbar[293] + 1./3.*armttbar[289] + armttbar[224];
   armttbar[215]=MMH*armttbar[215];
   armttbar[224]=17*armttbar[14] + 5./2.*armttbar[16];
   armttbar[224]=armttbar[16]*armttbar[224];
   armttbar[289]= - 17*armttbar[14];
   armttbar[293]=armttbar[289] - 11./2.*armttbar[16];
   armttbar[293]=1./3.*armttbar[293] + armttbar[15];
   armttbar[293]=armttbar[15]*armttbar[293];
   armttbar[295]=armttbar[16] - armttbar[15];
   armttbar[296]=armttbar[17]*armttbar[295];
   armttbar[224]=17./3.*armttbar[296] + 1./3.*armttbar[224] + 
   armttbar[293];
   armttbar[293]=armttbar[16]*armttbar[288];
   armttbar[297]=armttbar[15]*armttbar[288];
   armttbar[298]=armttbar[17]*armttbar[214];
   armttbar[293]=1./2.*armttbar[298] + armttbar[293] + 1./2.*
   armttbar[297];
   armttbar[293]=MMH*armttbar[293];
   armttbar[297]= - armttbar[16] - armttbar[15];
   armttbar[297]=armttbar[15]*armttbar[297];
   armttbar[297]=armttbar[186] + 1./2.*armttbar[297];
   armttbar[299]= - armttbar[16] + armttbar[15];
   armttbar[300]=armttbar[17]*armttbar[299];
   armttbar[301]=1./2.*armttbar[300];
   armttbar[293]=armttbar[293] + armttbar[297] + armttbar[301];
   armttbar[293]=armttbar[6]*MMH*armttbar[293];
   armttbar[215]=1./12.*armttbar[293] + 1./4.*armttbar[224] + 
   armttbar[215];
   armttbar[224]=armttbar[202] + armttbar[8];
   armttbar[293]=armttbar[16]*armttbar[224];
   armttbar[302]=armttbar[130] - armttbar[8];
   armttbar[304]=armttbar[15]*armttbar[302];
   armttbar[293]=armttbar[293] + armttbar[304];
   armttbar[304]=armttbar[21]*armttbar[288];
   armttbar[305]=armttbar[8]*armttbar[288];
   armttbar[306]=1./4.*armttbar[304] + armttbar[305];
   armttbar[306]=MMH*armttbar[306];
   armttbar[293]=17./4.*armttbar[293] + armttbar[306];
   armttbar[293]=MMt*armttbar[293];
   armttbar[215]=1./2.*armttbar[215] + 1./3.*armttbar[293];
   armttbar[215]=armttbar[210]*armttbar[215];
   armttbar[149]=1./2.*armttbar[215] + 1./2.*armttbar[149] + 
   armttbar[269];
   armttbar[149]=armttbar[210]*armttbar[149];
   armttbar[215]= - 1./3.*armttbar[38];
   armttbar[150]=1./2.*armttbar[8] + armttbar[275] + armttbar[150] + 
   armttbar[215];
   armttbar[150]=armttbar[8]*armttbar[150];
   armttbar[269]= - armttbar[21]*armttbar[38];
   armttbar[150]=armttbar[150] + 1./12.*armttbar[269] - 17./54.*
   armttbar[36] + 5./54.*armttbar[38] + armttbar[140] - 11./12.*
   armttbar[13] + armttbar[139] + armttbar[135] + 13./12.*armttbar[84]
    + armttbar[126] + 11./12.*armttbar[79] + armttbar[117] + 
   armttbar[114] + armttbar[112] + armttbar[134] + armttbar[137] + 
   armttbar[78] + armttbar[111] + armttbar[106] - 301./32. + 
   armttbar[105];
   armttbar[293]= - 1./2.*armttbar[18];
   armttbar[167]=armttbar[14] + armttbar[20] + armttbar[293] - 
   armttbar[92] + armttbar[167];
   armttbar[167]=armttbar[99]*armttbar[167];
   armttbar[306]=armttbar[38] + 13./8.*armttbar[36];
   armttbar[306]=armttbar[12]*armttbar[306];
   armttbar[307]= - 1./2.*armttbar[36];
   armttbar[308]=1 + armttbar[307];
   armttbar[308]=armttbar[15]*armttbar[99]*armttbar[308];
   armttbar[309]=1./2. + armttbar[36];
   armttbar[309]=armttbar[17]*armttbar[99]*armttbar[309];
   armttbar[310]= - 1 - armttbar[84];
   armttbar[310]=MMH*armttbar[99]*armttbar[310];
   armttbar[150]=5./24.*armttbar[310] + 5./12.*armttbar[309] + 5./12.*
   armttbar[308] + 5./12.*armttbar[167] + 1./4.*armttbar[150] + 1./27.*
   armttbar[306];
   armttbar[150]=MMH*armttbar[150];
   armttbar[167]=307./9. + armttbar[274];
   armttbar[167]=armttbar[199] + 1./2.*armttbar[167] + armttbar[40];
   armttbar[167]=armttbar[263] + 1./4.*armttbar[167] + armttbar[272];
   armttbar[167]=armttbar[14]*armttbar[167];
   armttbar[167]=47./72.*armttbar[183] + 1./2.*armttbar[228] + 5./24.*
   armttbar[181] + armttbar[167] + armttbar[116] - 31./24.*armttbar[18]
    + armttbar[203] + armttbar[160] + 97./24.*armttbar[42];
   armttbar[274]=5./8.*armttbar[8] - 5./8.*armttbar[21] + armttbar[107]
    - 1./72. + 5*armttbar[38];
   armttbar[306]= - 1./9.*armttbar[12];
   armttbar[274]=1./2.*armttbar[274] + armttbar[306];
   armttbar[274]=armttbar[16]*armttbar[274];
   armttbar[308]=331./54. + armttbar[201];
   armttbar[308]= - 125./216.*armttbar[12] + 1./4.*armttbar[308] - 
   armttbar[8];
   armttbar[308]=armttbar[15]*armttbar[308];
   armttbar[167]=1./2.*armttbar[308] + 1./2.*armttbar[167] + 1./3.*
   armttbar[274];
   armttbar[274]=armttbar[225] + 1./12.*armttbar[36] + armttbar[205] + 
   1./12.*armttbar[38];
   armttbar[274]=armttbar[17]*armttbar[274];
   armttbar[150]=1./2.*armttbar[150] + 1./2.*armttbar[167] + 
   armttbar[274];
   armttbar[150]=MMH*armttbar[150];
   armttbar[167]=1./8.*armttbar[85];
   armttbar[226]=armttbar[262] + 1./2.*armttbar[239] + armttbar[258] + 
   1./2.*armttbar[254] + 1./32.*armttbar[261] + 13./16.*armttbar[51] + 
   armttbar[167] + armttbar[226] + 2*armttbar[44];
   armttbar[226]=MMt*armttbar[226];
   armttbar[239]= - armttbar[8] + armttbar[263] + 1./4.*armttbar[276]
    + armttbar[272];
   armttbar[239]=armttbar[8]*armttbar[239];
   armttbar[254]=armttbar[292] + armttbar[36];
   armttbar[254]=armttbar[21]*armttbar[254];
   armttbar[239]=1./4.*armttbar[285] + 1./4.*armttbar[282] + 
   armttbar[239] + armttbar[136] + 1./12.*armttbar[254];
   armttbar[239]=MMH*armttbar[239];
   armttbar[254]=armttbar[238] + armttbar[267] + armttbar[204] + 69./4.
   *armttbar[18] + armttbar[192] - 307./36.*armttbar[19];
   armttbar[258]=205 - 143./2.*armttbar[21];
   armttbar[258]=1./3.*armttbar[258] - 5*armttbar[8];
   armttbar[258]=armttbar[16]*armttbar[258];
   armttbar[254]=1./6.*armttbar[258] + 1./2.*armttbar[254] + 
   armttbar[259];
   armttbar[257]=1./4.*armttbar[257] + 2*armttbar[260];
   armttbar[257]=armttbar[17]*armttbar[257];
   armttbar[258]= - 9 + 37./6.*armttbar[21];
   armttbar[258]=1./16.*armttbar[258] + armttbar[8];
   armttbar[258]=armttbar[15]*armttbar[258];
   armttbar[226]=armttbar[226] + 1./2.*armttbar[239] + armttbar[257] + 
   1./4.*armttbar[254] + armttbar[258];
   armttbar[226]=MMt*armttbar[226];
   armttbar[230]=armttbar[38] + armttbar[230];
   armttbar[239]=armttbar[15]*armttbar[230];
   armttbar[254]= - armttbar[38] + 5./8.*armttbar[36];
   armttbar[254]=armttbar[17]*armttbar[254];
   armttbar[239]=1./9.*armttbar[254] + 5./4.*armttbar[182] + 1./9.*
   armttbar[239];
   armttbar[239]=MMH*armttbar[239];
   armttbar[243]=5./24.*armttbar[15] + 1./8.*armttbar[14] + 
   armttbar[243];
   armttbar[243]=armttbar[15]*armttbar[243];
   armttbar[254]=283./24.*armttbar[15] - 1./8.*armttbar[14] + 
   armttbar[256];
   armttbar[254]=armttbar[17]*armttbar[254];
   armttbar[239]=armttbar[239] + 1./3.*armttbar[254] + 5./4.*
   armttbar[253] + 1./3.*armttbar[243];
   armttbar[239]=armttbar[6]*MMH*armttbar[239];
   armttbar[243]=53*armttbar[14] - 779./6.*armttbar[16];
   armttbar[243]=armttbar[16]*armttbar[243];
   armttbar[243]=armttbar[249] + 1./3.*armttbar[243];
   armttbar[254]=37./2.*armttbar[14] + 61*armttbar[16];
   armttbar[254]=1./3.*armttbar[254] - 15./4.*armttbar[15];
   armttbar[254]=armttbar[15]*armttbar[254];
   armttbar[243]=1./2.*armttbar[243] + armttbar[254];
   armttbar[254]= - 29./3.*armttbar[15] - armttbar[14] - 35./6.*
   armttbar[16];
   armttbar[254]=1./4.*armttbar[254] + armttbar[17];
   armttbar[254]=armttbar[17]*armttbar[254];
   armttbar[243]=1./8.*armttbar[243] + armttbar[254];
   armttbar[150]=armttbar[226] + 1./12.*armttbar[239] + 1./2.*
   armttbar[243] + armttbar[150];
   armttbar[149]=armttbar[150] + armttbar[149];
   armttbar[149]=armttbar[210]*armttbar[149];
   armttbar[226]=1./6.*armttbar[269];
   armttbar[105]=armttbar[226] + 5./108.*armttbar[36] + 19./54.*
   armttbar[38] + armttbar[140] + 7./36.*armttbar[13] + armttbar[139]
    + armttbar[135] - 41./36.*armttbar[84] + armttbar[126] - 7./36.*
   armttbar[79] + armttbar[117] + armttbar[114] + armttbar[112] + 
   armttbar[134] + armttbar[137] + armttbar[78] + armttbar[111] + 
   armttbar[106] - 863./96. + armttbar[105];
   armttbar[106]=armttbar[159] + armttbar[155] + armttbar[215];
   armttbar[106]=armttbar[8]*armttbar[106];
   armttbar[111]=armttbar[272] + armttbar[307];
   armttbar[111]=armttbar[12]*armttbar[111];
   armttbar[105]=25./36.*armttbar[193] + 25./18.*armttbar[178] + 25./18.
   *armttbar[169] + 25./18.*armttbar[163] + 7./36.*armttbar[111] + 1./2.
   *armttbar[105] + armttbar[106];
   armttbar[105]=MMH*armttbar[105];
   armttbar[106]=armttbar[116] + 227./72.*armttbar[18] + armttbar[203]
    + armttbar[160] + 211./72.*armttbar[42];
   armttbar[111]=armttbar[212] + armttbar[157] + 277./27. + 9./4.*
   armttbar[35];
   armttbar[111]=1./4.*armttbar[111] + armttbar[272];
   armttbar[111]=armttbar[14]*armttbar[111];
   armttbar[106]=7./108.*armttbar[156] + armttbar[234] + 1./12.*
   armttbar[181] + 1./2.*armttbar[106] + armttbar[111];
   armttbar[111]=armttbar[225] + armttbar[205] + armttbar[278];
   armttbar[111]=armttbar[17]*armttbar[111];
   armttbar[114]= - 7./144.*armttbar[12];
   armttbar[116]=armttbar[114] + 85./32.*armttbar[8] + armttbar[219] - 
   29./288. + armttbar[38];
   armttbar[116]=armttbar[16]*armttbar[116];
   armttbar[117]= - 17./27. - armttbar[38];
   armttbar[117]=7./36.*armttbar[12] + 13./9.*armttbar[8] + 1./2.*
   armttbar[117] + armttbar[36];
   armttbar[117]=armttbar[15]*armttbar[117];
   armttbar[105]=1./4.*armttbar[105] + armttbar[111] + 1./4.*
   armttbar[117] + 1./2.*armttbar[106] + 1./3.*armttbar[116];
   armttbar[105]=MMH*armttbar[105];
   armttbar[106]=armttbar[272] - armttbar[36];
   armttbar[111]=armttbar[15]*armttbar[106];
   armttbar[116]=armttbar[215] + armttbar[36];
   armttbar[116]=armttbar[17]*armttbar[116];
   armttbar[111]=17./12.*armttbar[116] + armttbar[182] + 17./12.*
   armttbar[111];
   armttbar[111]=MMH*armttbar[111];
   armttbar[116]= - 1./2.*armttbar[16];
   armttbar[117]= - armttbar[14] + armttbar[116];
   armttbar[126]=1./3.*armttbar[117] + armttbar[188];
   armttbar[126]=armttbar[15]*armttbar[126];
   armttbar[134]=armttbar[14] + 1./2.*armttbar[16];
   armttbar[135]=17./3.*armttbar[134] + armttbar[195];
   armttbar[135]=armttbar[17]*armttbar[135];
   armttbar[111]=armttbar[111] + 1./6.*armttbar[135] + armttbar[253] + 
   17./6.*armttbar[126];
   armttbar[111]=armttbar[6]*MMH*armttbar[111];
   armttbar[126]= - 7*armttbar[14] - 53./6.*armttbar[16];
   armttbar[126]=armttbar[16]*armttbar[126];
   armttbar[135]= - 11*armttbar[14];
   armttbar[137]= - 19./2.*armttbar[15] + armttbar[135] + 31*
   armttbar[16];
   armttbar[137]=armttbar[15]*armttbar[137];
   armttbar[126]=armttbar[137] + armttbar[249] + 17./3.*armttbar[126];
   armttbar[137]=1./9.*armttbar[15] - armttbar[14] - 115./6.*
   armttbar[16];
   armttbar[137]=1./4.*armttbar[137] + armttbar[17];
   armttbar[137]=armttbar[17]*armttbar[137];
   armttbar[126]=1./16.*armttbar[126] + armttbar[137];
   armttbar[105]=1./12.*armttbar[111] + 1./2.*armttbar[126] + 
   armttbar[105];
   armttbar[111]=armttbar[238] + armttbar[267] + armttbar[204] - 41./18.
   *armttbar[18] + armttbar[192] - 59./9.*armttbar[19];
   armttbar[126]=107 - 65*armttbar[21];
   armttbar[126]=7./12.*armttbar[126] - 85*armttbar[8];
   armttbar[126]=armttbar[16]*armttbar[126];
   armttbar[111]=1./6.*armttbar[126] + 1./2.*armttbar[111] + 
   armttbar[259];
   armttbar[126]= - 1./16. + 11./3.*armttbar[21];
   armttbar[126]=1./2.*armttbar[126] - 13./3.*armttbar[8];
   armttbar[126]=armttbar[15]*armttbar[126];
   armttbar[111]=1./4.*armttbar[111] + 1./3.*armttbar[126];
   armttbar[126]= - 1./4. + armttbar[38];
   armttbar[126]=armttbar[21]*armttbar[126];
   armttbar[136]=armttbar[136] + 1./6.*armttbar[126];
   armttbar[137]=armttbar[244] + armttbar[280] + armttbar[272];
   armttbar[137]=armttbar[8]*armttbar[137];
   armttbar[136]=armttbar[286] + armttbar[284] + 1./2.*armttbar[136] + 
   armttbar[137];
   armttbar[136]=MMH*armttbar[136];
   armttbar[111]=armttbar[266] + 1./2.*armttbar[136] + 1./2.*
   armttbar[111] + armttbar[281];
   armttbar[111]=MMt*armttbar[111];
   armttbar[136]=armttbar[14] - armttbar[16];
   armttbar[137]=armttbar[16]*armttbar[136];
   armttbar[139]=armttbar[15]*armttbar[252];
   armttbar[137]=1./4.*armttbar[296] + 1./4.*armttbar[137] + 
   armttbar[139];
   armttbar[139]=armttbar[14] + 7*armttbar[183];
   armttbar[140]= - 7./9.*armttbar[12];
   armttbar[155]=armttbar[140] - 1./9. - armttbar[8];
   armttbar[155]=armttbar[16]*armttbar[155];
   armttbar[139]=1./9.*armttbar[139] + armttbar[155];
   armttbar[155]=1./6.*armttbar[8];
   armttbar[160]= - armttbar[38] + armttbar[155];
   armttbar[160]=armttbar[15]*armttbar[160];
   armttbar[163]=armttbar[12]*armttbar[38];
   armttbar[169]=armttbar[38] + 7*armttbar[163];
   armttbar[169]=MMH*armttbar[169];
   armttbar[139]=1./54.*armttbar[169] + 1./6.*armttbar[139] + 
   armttbar[160];
   armttbar[139]=MMH*armttbar[139];
   armttbar[160]=armttbar[15]*armttbar[38];
   armttbar[169]= - armttbar[17]*armttbar[38];
   armttbar[160]=armttbar[160] + armttbar[169];
   armttbar[160]=MMH*armttbar[160];
   armttbar[178]=armttbar[15]*armttbar[136];
   armttbar[182]=armttbar[17]*armttbar[252];
   armttbar[160]=armttbar[160] + armttbar[178] + armttbar[182];
   armttbar[160]=armttbar[6]*MMH*armttbar[160];
   armttbar[178]=1./3.*armttbar[21];
   armttbar[192]= - 1./2. + armttbar[178];
   armttbar[193]= - 1./3.*armttbar[8];
   armttbar[192]=1./2.*armttbar[192] + armttbar[193];
   armttbar[192]=armttbar[15]*armttbar[192];
   armttbar[195]= - 1./3.*armttbar[21];
   armttbar[203]=1./2. + armttbar[195];
   armttbar[204]=1./3.*armttbar[8];
   armttbar[203]=1./2.*armttbar[203] + armttbar[204];
   armttbar[203]=armttbar[16]*armttbar[203];
   armttbar[205]= - armttbar[19] + armttbar[18];
   armttbar[192]=armttbar[192] + 1./4.*armttbar[205] + armttbar[203];
   armttbar[192]=MMt*armttbar[192];
   armttbar[137]=1./2.*armttbar[192] + 17./216.*armttbar[160] + 1./3.*
   armttbar[137] + 1./4.*armttbar[139];
   armttbar[137]=armttbar[176]*armttbar[137];
   armttbar[105]=1./4.*armttbar[137] + 1./2.*armttbar[105] + 
   armttbar[111];
   armttbar[105]=armttbar[176]*armttbar[105];
   armttbar[105]=armttbar[150] + armttbar[105];
   armttbar[105]=armttbar[176]*armttbar[105];
   armttbar[111]=armttbar[8]*armttbar[38];
   armttbar[137]= - 1 - armttbar[85];
   armttbar[139]=armttbar[137] + 1./6.*armttbar[111];
   armttbar[139]=MMH*armttbar[139];
   armttbar[150]= - 1./2.*armttbar[19];
   armttbar[127]=armttbar[150] - armttbar[95] + armttbar[127];
   armttbar[144]=1 + armttbar[144];
   armttbar[160]=armttbar[14]*armttbar[144];
   armttbar[192]= - 1./12.*armttbar[8] + 1 - 1./12.*armttbar[38];
   armttbar[192]=armttbar[16]*armttbar[192];
   armttbar[139]=1./2.*armttbar[139] + 1./6.*armttbar[169] + 
   armttbar[192] + 1./12.*armttbar[259] + armttbar[127] + armttbar[160]
   ;
   armttbar[139]=MMH*armttbar[139];
   armttbar[160]=7*armttbar[14];
   armttbar[192]= - 5*armttbar[16];
   armttbar[205]=armttbar[160] + armttbar[192];
   armttbar[205]=armttbar[16]*armttbar[205];
   armttbar[205]=armttbar[182] - armttbar[236] + 1./2.*armttbar[205];
   armttbar[139]=1./6.*armttbar[205] + armttbar[139];
   armttbar[139]=MMH*armttbar[139];
   armttbar[111]=armttbar[137] + 1./3.*armttbar[111];
   armttbar[111]=MMH*armttbar[111];
   armttbar[205]=1 + armttbar[215];
   armttbar[205]=armttbar[14]*armttbar[205];
   armttbar[144]=armttbar[144] - 1./6.*armttbar[8];
   armttbar[144]=armttbar[16]*armttbar[144];
   armttbar[215]=1./6.*armttbar[259];
   armttbar[111]=1./2.*armttbar[111] + 1./3.*armttbar[169] + 
   armttbar[144] + armttbar[215] + armttbar[127] + armttbar[205];
   armttbar[111]=MMH*armttbar[111];
   armttbar[116]=armttbar[14] + armttbar[116];
   armttbar[116]=armttbar[16]*armttbar[116];
   armttbar[116]=1./2.*armttbar[182] - 1./2.*armttbar[236] + 
   armttbar[116];
   armttbar[111]=1./3.*armttbar[116] + 1./2.*armttbar[111];
   armttbar[111]=MMH*armttbar[111];
   armttbar[116]=1./3.*armttbar[181];
   armttbar[144]= - armttbar[14] + armttbar[42] - armttbar[19];
   armttbar[169]=1./2.*armttbar[144] + armttbar[116];
   armttbar[182]= - armttbar[23] + armttbar[85] + 1./4. + armttbar[82];
   armttbar[205]=armttbar[182] + 1./3.*armttbar[269];
   armttbar[225]= - 1./3.*armttbar[8]*armttbar[38];
   armttbar[205]=1./4.*armttbar[205] + armttbar[225];
   armttbar[205]=MMH*armttbar[205];
   armttbar[234]=1./3.*armttbar[228];
   armttbar[169]=armttbar[205] + armttbar[203] + 1./2.*armttbar[169] + 
   armttbar[234];
   armttbar[169]=MMt*MMH*armttbar[169];
   armttbar[111]=armttbar[111] + armttbar[169];
   armttbar[111]=armttbar[176]*armttbar[111];
   armttbar[169]=armttbar[144] + 5./6.*armttbar[181];
   armttbar[203]=1 - 5./6.*armttbar[21];
   armttbar[203]=1./2.*armttbar[203] + armttbar[204];
   armttbar[203]=armttbar[16]*armttbar[203];
   armttbar[205]=armttbar[182] + armttbar[226];
   armttbar[205]=1./2.*armttbar[205] + armttbar[225];
   armttbar[205]=MMH*armttbar[205];
   armttbar[169]=armttbar[205] + armttbar[203] + 1./2.*armttbar[169] + 
   armttbar[234];
   armttbar[169]=MMt*MMH*armttbar[169];
   armttbar[111]=armttbar[111] + armttbar[139] + armttbar[169];
   armttbar[111]=armttbar[176]*armttbar[111];
   armttbar[139]=armttbar[8]*armttbar[36];
   armttbar[139]=3*armttbar[137] + 1./3.*armttbar[139];
   armttbar[139]=MMH*armttbar[139];
   armttbar[169]=3 + armttbar[275];
   armttbar[169]=armttbar[14]*armttbar[169];
   armttbar[203]= - 1./6.*armttbar[36];
   armttbar[205]=3 + armttbar[203];
   armttbar[205]=armttbar[16]*armttbar[205];
   armttbar[225]= - armttbar[17]*armttbar[36];
   armttbar[226]= - armttbar[15]*armttbar[8];
   armttbar[234]=1./6.*armttbar[226];
   armttbar[139]=1./2.*armttbar[139] + 1./3.*armttbar[225] + 
   armttbar[234] + armttbar[205] + armttbar[215] + 3*armttbar[127] + 
   armttbar[169];
   armttbar[139]=MMH*armttbar[139];
   armttbar[169]=1./3.*armttbar[14] - 3./8.*armttbar[16];
   armttbar[169]=armttbar[16]*armttbar[169];
   armttbar[205]= - armttbar[14] + armttbar[15];
   armttbar[205]=armttbar[17]*armttbar[205];
   armttbar[134]=armttbar[15]*armttbar[134];
   armttbar[139]=1./4.*armttbar[139] + 1./12.*armttbar[205] + 1./12.*
   armttbar[134] - 1./12.*armttbar[236] + armttbar[169];
   armttbar[139]=MMH*armttbar[139];
   armttbar[169]= - armttbar[21]*armttbar[36];
   armttbar[169]=3*armttbar[182] + 1./3.*armttbar[169];
   armttbar[205]= - armttbar[8]*armttbar[36];
   armttbar[169]=1./4.*armttbar[169] + 1./3.*armttbar[205];
   armttbar[169]=MMH*armttbar[169];
   armttbar[205]=1 - armttbar[21];
   armttbar[205]=armttbar[16]*armttbar[205];
   armttbar[215]=armttbar[15]*armttbar[224];
   armttbar[116]=1./2.*armttbar[169] + 1./6.*armttbar[215] + 3./8.*
   armttbar[205] + 1./6.*armttbar[228] + 3./8.*armttbar[144] + 
   armttbar[116];
   armttbar[116]=MMt*MMH*armttbar[116];
   armttbar[116]=armttbar[139] + armttbar[116];
   armttbar[111]=armttbar[116] + 1./2.*armttbar[111];
   armttbar[111]=armttbar[176]*armttbar[111];
   armttbar[139]=armttbar[155] + armttbar[203] + 1 + armttbar[278];
   armttbar[139]=armttbar[16]*armttbar[139];
   armttbar[155]=armttbar[275] + 1 + armttbar[272];
   armttbar[155]=armttbar[14]*armttbar[155];
   armttbar[137]=armttbar[137] + 1./3.*armttbar[305];
   armttbar[137]=MMH*armttbar[137];
   armttbar[127]=1./2.*armttbar[137] + 1./3.*armttbar[298] + 
   armttbar[234] + armttbar[139] + armttbar[127] + armttbar[155];
   armttbar[127]=MMH*armttbar[127];
   armttbar[137]=1./4.*armttbar[14];
   armttbar[139]=armttbar[137] - armttbar[16];
   armttbar[139]=armttbar[16]*armttbar[139];
   armttbar[139]=armttbar[301] + armttbar[139] + 1./2.*armttbar[134];
   armttbar[127]=1./3.*armttbar[139] + 1./2.*armttbar[127];
   armttbar[127]=MMH*armttbar[127];
   armttbar[139]=armttbar[14]*armttbar[214];
   armttbar[155]=armttbar[214] + armttbar[8];
   armttbar[155]=armttbar[16]*armttbar[155];
   armttbar[169]=MMH*armttbar[305];
   armttbar[139]=1./2.*armttbar[169] + armttbar[298] + 1./2.*
   armttbar[226] + armttbar[139] + 1./2.*armttbar[155];
   armttbar[139]=MMH*armttbar[139];
   armttbar[117]=armttbar[16]*armttbar[117];
   armttbar[117]=armttbar[139] + armttbar[300] + armttbar[117] + 
   armttbar[134];
   armttbar[117]=MMH*armttbar[117];
   armttbar[134]=1./4.*armttbar[216] + armttbar[290];
   armttbar[134]=MMH*armttbar[134];
   armttbar[139]=armttbar[16]*armttbar[302];
   armttbar[134]=armttbar[134] + armttbar[139] + armttbar[215];
   armttbar[134]=MMt*MMH*armttbar[134];
   armttbar[117]=1./2.*armttbar[117] + armttbar[134];
   armttbar[117]=armttbar[210]*armttbar[117];
   armttbar[134]=armttbar[144] + armttbar[181];
   armttbar[139]=armttbar[193] + 1./4. + armttbar[195];
   armttbar[139]=armttbar[16]*armttbar[139];
   armttbar[144]=armttbar[182] + armttbar[287];
   armttbar[144]=1./4.*armttbar[144] + armttbar[291];
   armttbar[144]=MMH*armttbar[144];
   armttbar[134]=armttbar[144] + 1./3.*armttbar[215] + 1./4.*
   armttbar[134] + armttbar[139];
   armttbar[134]=MMt*MMH*armttbar[134];
   armttbar[117]=1./3.*armttbar[117] + armttbar[127] + armttbar[134];
   armttbar[117]=armttbar[210]*armttbar[117];
   armttbar[116]=armttbar[116] + 1./2.*armttbar[117];
   armttbar[116]=armttbar[210]*armttbar[116];
   armttbar[111]=armttbar[111] + armttbar[116];
   armttbar[111]=armttbar[9]*armttbar[111];
   armttbar[105]=1./4.*armttbar[111] + armttbar[105] + armttbar[149];
   armttbar[105]=armttbar[9]*armttbar[105];
   armttbar[111]=77./36.*armttbar[76] + 691./144. + 15*armttbar[83];
   armttbar[111]= - 1./4.*armttbar[90] - 1./6.*armttbar[141] + 1./4.*
   armttbar[111] - 1./3.*armttbar[80];
   armttbar[116]=1./2.*armttbar[75];
   armttbar[111]=239./144.*armttbar[87] + 7./36.*armttbar[72] + 7./36.*
   armttbar[74] + armttbar[116] + 1./2.*armttbar[111] + armttbar[77];
   armttbar[111]= - 23./72.*armttbar[88] + 1./2.*armttbar[111] + 17./9.
   *armttbar[44];
   armttbar[117]= - 17./3.*armttbar[39];
   armttbar[127]=7./4. + armttbar[117];
   armttbar[127]=1./2.*armttbar[127] + armttbar[178];
   armttbar[127]=armttbar[21]*armttbar[127];
   armttbar[134]= - 9 - 7./36.*armttbar[21];
   armttbar[139]=1./4.*armttbar[12];
   armttbar[134]=armttbar[139] + 1./8.*armttbar[134] - 7./9.*
   armttbar[8];
   armttbar[134]=armttbar[12]*armttbar[134];
   armttbar[144]=7*armttbar[156] + 61*armttbar[42] + armttbar[160];
   armttbar[149]= - 1 + armttbar[202];
   armttbar[149]=1./2.*armttbar[149] + armttbar[8];
   armttbar[155]=armttbar[15]*armttbar[149];
   armttbar[160]=3*armttbar[155];
   armttbar[144]=1./18.*armttbar[144] + armttbar[160];
   armttbar[144]=armttbar[3]*armttbar[144];
   armttbar[169]= - 1./2.*armttbar[68] - armttbar[65];
   armttbar[169]=1./2.*armttbar[169] - 7./3.*armttbar[63];
   armttbar[178]= - 1./3.*armttbar[67];
   armttbar[169]=1./2.*armttbar[169] + armttbar[178];
   armttbar[169]=MMt*armttbar[169];
   armttbar[182]=77./2. + armttbar[110];
   armttbar[182]=1./3.*armttbar[182] - armttbar[37];
   armttbar[182]=armttbar[8]*armttbar[182];
   armttbar[193]=13./6.*armttbar[63] - armttbar[64];
   armttbar[193]=1./2.*armttbar[193] - armttbar[60];
   armttbar[193]=MMH*armttbar[193];
   armttbar[195]=armttbar[101] + armttbar[283];
   armttbar[195]=armttbar[15]*armttbar[195];
   armttbar[111]=1./6.*armttbar[169] + 1./4.*armttbar[193] + 1./4.*
   armttbar[144] + 3./8.*armttbar[195] + 1./4.*armttbar[134] + 1./8.*
   armttbar[182] + 1./16.*armttbar[127] + 1./8.*armttbar[37] - 55./384.
   *armttbar[13] + armttbar[148] - 1./18.*armttbar[25] - 15./32.*
   armttbar[23] - 7./72.*armttbar[84] - 43./144.*armttbar[51] + 1./2.*
   armttbar[111] - 4./9.*armttbar[45];
   armttbar[111]=MMt*armttbar[111];
   armttbar[127]=MMH*armttbar[163];
   armttbar[127]=armttbar[127] + armttbar[183] + armttbar[200];
   armttbar[127]=MMH*armttbar[127];
   armttbar[134]=armttbar[15]*armttbar[299];
   armttbar[127]=armttbar[127] + armttbar[134] + armttbar[296];
   armttbar[127]=armttbar[6]*armttbar[127];
   armttbar[134]=1./3.*armttbar[146] - armttbar[23] + armttbar[87] + 1./
   4. + armttbar[83];
   armttbar[144]= - armttbar[8]*armttbar[39];
   armttbar[134]=1./4.*armttbar[134] + 1./3.*armttbar[144];
   armttbar[134]=MMt*armttbar[134];
   armttbar[144]=armttbar[150] + 1./2.*armttbar[50] - armttbar[96];
   armttbar[146]= - 1 + armttbar[171];
   armttbar[146]=armttbar[14]*armttbar[146];
   armttbar[140]=armttbar[140] + 89./9. - armttbar[39];
   armttbar[140]=armttbar[16]*armttbar[140];
   armttbar[148]=55 + 7*armttbar[12];
   armttbar[148]=armttbar[15]*armttbar[148];
   armttbar[150]=armttbar[3]*armttbar[15]*armttbar[295];
   armttbar[163]=armttbar[8]*armttbar[39];
   armttbar[163]= - armttbar[38] + 1./4.*armttbar[163];
   armttbar[163]=MMH*armttbar[163];
   armttbar[127]=armttbar[134] + 17./108.*armttbar[127] + 1./3.*
   armttbar[163] + 1./6.*armttbar[187] + 1./2.*armttbar[150] + 1./108.*
   armttbar[148] + 1./12.*armttbar[140] + 1./2.*armttbar[144] + 1./3.*
   armttbar[146];
   armttbar[127]=armttbar[176]*armttbar[127];
   armttbar[134]= - 7./2.*armttbar[92] + armttbar[48];
   armttbar[138]= - 1./3.*armttbar[18] + 1./3.*armttbar[134] + 
   armttbar[138];
   armttbar[138]=7./12.*armttbar[14] + 1./2.*armttbar[138] + 1./3.*
   armttbar[20];
   armttbar[138]=armttbar[99]*armttbar[138];
   armttbar[140]=armttbar[12]*armttbar[36];
   armttbar[144]=1./6.*armttbar[140] + armttbar[203] + 1./6.*
   armttbar[13] - 7./3.*armttbar[84] - 1./6.*armttbar[79] - 31./16. + 1.
   /3.*armttbar[88];
   armttbar[144]=armttbar[99]*armttbar[144];
   armttbar[144]= - 1./3.*armttbar[63] + 25./16.*armttbar[144];
   armttbar[144]=MMH*armttbar[144];
   armttbar[146]= - 3*armttbar[77];
   armttbar[148]= - 97./36.*armttbar[72] - 97./36.*armttbar[74] - 2213./
   432. + armttbar[146];
   armttbar[148]=25./36.*armttbar[79] + 91./72.*armttbar[88] + 1./4.*
   armttbar[148] - 17./9.*armttbar[44];
   armttbar[150]=17*armttbar[39];
   armttbar[163]= - 163./6. + armttbar[150];
   armttbar[163]=1./12.*armttbar[163] + armttbar[37];
   armttbar[163]=armttbar[8]*armttbar[163];
   armttbar[169]=151./4. + armttbar[206];
   armttbar[169]=1./3.*armttbar[169] + 7./2.*armttbar[36];
   armttbar[169]=1./3.*armttbar[169] - 23./16.*armttbar[8];
   armttbar[169]=armttbar[12]*armttbar[169];
   armttbar[171]=17./4. - armttbar[36];
   armttbar[171]=armttbar[15]*armttbar[99]*armttbar[171];
   armttbar[182]=7./8. + armttbar[36];
   armttbar[182]=armttbar[17]*armttbar[99]*armttbar[182];
   armttbar[138]=1./6.*armttbar[144] + 25./36.*armttbar[182] + 25./72.*
   armttbar[171] + 25./12.*armttbar[138] + 1./24.*armttbar[169] + 1./8.
   *armttbar[163] + 29./216.*armttbar[36] - 5./108.*armttbar[38] - 3./
   16.*armttbar[37] - 1./384.*armttbar[13] + 209./1728.*armttbar[84] + 
   205./1152.*armttbar[51] + 1./8.*armttbar[148] + 1./9.*armttbar[45];
   armttbar[138]=MMH*armttbar[138];
   armttbar[144]=47./6. + armttbar[110];
   armttbar[144]=1./3.*armttbar[144] - armttbar[37];
   armttbar[144]=armttbar[14]*armttbar[144];
   armttbar[110]=661./9. + armttbar[110];
   armttbar[110]=403./36.*armttbar[12] + 1./2.*armttbar[110] - 
   armttbar[21];
   armttbar[110]=armttbar[16]*armttbar[110];
   armttbar[148]=3./4.*armttbar[48];
   armttbar[163]=armttbar[99]*armttbar[186];
   armttbar[110]=25./18.*armttbar[163] + 1./12.*armttbar[110] + 11./72.
   *armttbar[156] + 1./4.*armttbar[144] + 383./288.*armttbar[20] - 
   35507./5184.*armttbar[18] - 103./48.*armttbar[19] - 103./432.*
   armttbar[42] - 1./36.*armttbar[5] + 107./144.*armttbar[33] - 155./
   288.*armttbar[34] - 23./18.*armttbar[32] + armttbar[148] + 11./36.*
   armttbar[92] - 7./12.*armttbar[93] - 413./144.*armttbar[96] + 15./8.
   *armttbar[50] + 1./3.*armttbar[94];
   armttbar[144]= - 27*armttbar[35];
   armttbar[169]=29653./648. + armttbar[144];
   armttbar[187]= - 3*armttbar[40];
   armttbar[169]=armttbar[118] + 1./2.*armttbar[169] + armttbar[187];
   armttbar[169]=1./4.*armttbar[169] + armttbar[38];
   armttbar[193]= - 3*armttbar[101];
   armttbar[203]=armttbar[193] + 625./9.*armttbar[99];
   armttbar[203]=armttbar[15]*armttbar[203];
   armttbar[205]= - armttbar[16]*armttbar[98];
   armttbar[169]=1./8.*armttbar[203] + 1./48.*armttbar[205] + 13./1296.
   *armttbar[12] - 197./144.*armttbar[8] - 103./1152.*armttbar[21] + 1./
   2.*armttbar[169] - armttbar[36];
   armttbar[169]=armttbar[15]*armttbar[169];
   armttbar[117]=19./2. + armttbar[117];
   armttbar[117]=79./36.*armttbar[12] + 17./18.*armttbar[8] + 97./144.*
   armttbar[21] + 1./2.*armttbar[117] - armttbar[37];
   armttbar[203]=armttbar[193] + 4825./54.*armttbar[99];
   armttbar[203]=armttbar[15]*armttbar[203];
   armttbar[206]=armttbar[99]*armttbar[16];
   armttbar[117]=1./2.*armttbar[203] + 1./2.*armttbar[117] + 25./9.*
   armttbar[206];
   armttbar[203]=9./8.*armttbar[15];
   armttbar[215]= - 17./9.*armttbar[14] + armttbar[203];
   armttbar[215]=armttbar[3]*armttbar[215];
   armttbar[216]= - armttbar[17]*armttbar[99];
   armttbar[117]=1075./216.*armttbar[216] + 1./2.*armttbar[117] + 
   armttbar[215];
   armttbar[117]=armttbar[17]*armttbar[117];
   armttbar[215]= - armttbar[14] - 25./4.*armttbar[16];
   armttbar[215]=3*armttbar[215] - 7./2.*armttbar[15];
   armttbar[215]=armttbar[3]*armttbar[15]*armttbar[215];
   armttbar[110]=armttbar[117] + 1./4.*armttbar[215] + 1./2.*
   armttbar[110] + armttbar[169];
   armttbar[117]=13./3.*armttbar[92] - armttbar[48];
   armttbar[117]=1./3.*armttbar[117] + armttbar[165];
   armttbar[169]=armttbar[211] + 85./24.*armttbar[8] + 17./3.*
   armttbar[36] - 37./4. - 17./9.*armttbar[38];
   armttbar[169]=armttbar[15]*armttbar[169];
   armttbar[117]=1./2.*armttbar[169] + 17./72.*armttbar[200] + 17./36.*
   armttbar[156] - 187./72.*armttbar[14] - 5./18.*armttbar[20] + 323./
   144.*armttbar[18] + 17./8.*armttbar[117] - 1./3.*armttbar[33];
   armttbar[143]= - 125./8. + armttbar[143];
   armttbar[143]= - 85./8.*armttbar[8] + 1./3.*armttbar[143] + 
   armttbar[227];
   armttbar[143]=1./4.*armttbar[143] - armttbar[12];
   armttbar[143]=armttbar[17]*armttbar[143];
   armttbar[106]=armttbar[12]*armttbar[106];
   armttbar[106]=11./6.*armttbar[185] + armttbar[106];
   armttbar[106]=MMH*armttbar[106];
   armttbar[106]=17./48.*armttbar[106] + 1./2.*armttbar[117] + 1./3.*
   armttbar[143];
   armttbar[106]=MMH*armttbar[106];
   armttbar[117]= - 5423./3.*armttbar[15] + armttbar[289] + 1085*
   armttbar[16];
   armttbar[117]=armttbar[15]*armttbar[117];
   armttbar[117]= - 11*armttbar[186] + 1./12.*armttbar[117];
   armttbar[143]= - 85*armttbar[14] - 1481./2.*armttbar[16];
   armttbar[143]=355./2.*armttbar[17] + 1./4.*armttbar[143] - 1531./3.*
   armttbar[15];
   armttbar[143]=armttbar[17]*armttbar[143];
   armttbar[117]=1./2.*armttbar[117] + 1./3.*armttbar[143];
   armttbar[106]=1./12.*armttbar[117] + armttbar[106];
   armttbar[106]=armttbar[6]*armttbar[106];
   armttbar[106]=1./8.*armttbar[127] + armttbar[111] + 1./6.*
   armttbar[106] + 1./2.*armttbar[110] + armttbar[138];
   armttbar[106]=armttbar[176]*armttbar[106];
   armttbar[110]=11./24.*armttbar[58] + 8185./768. - armttbar[81];
   armttbar[111]=11./12.*armttbar[84];
   armttbar[115]=armttbar[115] - armttbar[19];
   armttbar[115]=armttbar[100]*armttbar[115];
   armttbar[110]=1./16.*armttbar[166] + armttbar[212] + 1./8.*
   armttbar[115] + 73./192.*armttbar[13] + 1./4.*armttbar[40] + 
   armttbar[303] + 25./12.*armttbar[25] + 15./32.*armttbar[23] + 
   armttbar[111] - 403./288.*armttbar[51] + armttbar[45] + 19./24.*
   armttbar[88] + 3./16.*armttbar[85] + 10./9.*armttbar[44] + 13./48.*
   armttbar[87] - 11./24.*armttbar[72] - 11./24.*armttbar[74] + 
   armttbar[191] + armttbar[112] - 33./64.*armttbar[90] - 13./36.*
   armttbar[141] - 15./64.*armttbar[89] + 1./6.*armttbar[80] - 11./64.*
   armttbar[76] + 3./4.*armttbar[78] - 9./32.*armttbar[73] + 1./3.*
   armttbar[110] - 15./32.*armttbar[83];
   armttbar[112]=295./9. - armttbar[40];
   armttbar[117]= - 1./12.*armttbar[41];
   armttbar[127]=13./3.*armttbar[39];
   armttbar[138]=7./16.*armttbar[21];
   armttbar[112]=armttbar[138] - armttbar[36] - armttbar[38] + 
   armttbar[117] + armttbar[129] + 1./2.*armttbar[112] + armttbar[127];
   armttbar[112]=armttbar[8]*armttbar[112];
   armttbar[143]=13*armttbar[39];
   armttbar[169]= - 1./4.*armttbar[41] + 329./12. + armttbar[143];
   armttbar[185]=235./72.*armttbar[21];
   armttbar[169]=armttbar[185] + armttbar[246] - armttbar[36] + 1./3.*
   armttbar[169] - armttbar[38];
   armttbar[169]=armttbar[21]*armttbar[169];
   armttbar[191]=31./3.*armttbar[42] - 5*armttbar[14];
   armttbar[191]=5./3.*armttbar[191] + armttbar[181];
   armttbar[191]=1./2.*armttbar[191] + 11./3.*armttbar[183];
   armttbar[149]=3*armttbar[16]*armttbar[149];
   armttbar[160]=armttbar[160] + 1./2.*armttbar[191] + armttbar[149];
   armttbar[160]=1./2.*armttbar[3]*armttbar[160];
   armttbar[191]=125./3. + 33./4.*armttbar[21];
   armttbar[191]=armttbar[139] + 1./8.*armttbar[191] + 11./3.*
   armttbar[8];
   armttbar[191]=1./2.*armttbar[12]*armttbar[191];
   armttbar[161]=7*armttbar[100] + armttbar[161];
   armttbar[200]= - armttbar[8]*armttbar[100];
   armttbar[161]=1./2.*armttbar[161] + 3*armttbar[200];
   armttbar[161]=armttbar[16]*armttbar[161];
   armttbar[200]=1./4.*armttbar[161];
   armttbar[215]=1./4.*armttbar[70];
   armttbar[224]=1./3.*armttbar[52];
   armttbar[225]=1./2.*armttbar[65] - 1./8.*armttbar[68] - armttbar[69]
    + armttbar[224] - 16./3.*armttbar[53] + armttbar[215];
   armttbar[172]=armttbar[172] + 11./6.*armttbar[63] + 1./3.*
   armttbar[225] + 1./4.*armttbar[61];
   armttbar[172]=MMt*armttbar[172];
   armttbar[225]= - 4./3.*armttbar[56] + armttbar[53];
   armttbar[225]= - 1./2.*armttbar[60] - 1./4.*armttbar[64] - 9./8.*
   armttbar[63] + 4./3.*armttbar[225] - 3./8.*armttbar[61];
   armttbar[225]=MMH*armttbar[225];
   armttbar[195]=3./4.*armttbar[195];
   armttbar[112]=armttbar[172] + armttbar[225] + armttbar[160] + 
   armttbar[195] + armttbar[200] + armttbar[191] + 1./2.*armttbar[112]
    + armttbar[110] + 1./8.*armttbar[169];
   armttbar[112]=MMt*armttbar[112];
   armttbar[169]= - 541./54. + armttbar[144];
   armttbar[226]= - 5*armttbar[38];
   armttbar[227]= - 13./2.*armttbar[8];
   armttbar[169]= - 925./108.*armttbar[12] + armttbar[227] + 45./8.*
   armttbar[21] + armttbar[120] - armttbar[36] + armttbar[226] + 
   armttbar[117] + armttbar[118] + armttbar[127] + 1./2.*armttbar[169]
    + armttbar[187];
   armttbar[234]=51./16.*armttbar[98] - armttbar[100];
   armttbar[234]=armttbar[16]*armttbar[234];
   armttbar[169]=1./4.*armttbar[169] + armttbar[234];
   armttbar[169]=armttbar[16]*armttbar[169];
   armttbar[238]= - 3*armttbar[73];
   armttbar[239]= - 1289./12. + armttbar[238];
   armttbar[243]=9./4.*armttbar[89];
   armttbar[146]= - 19./12.*armttbar[72] - 19./12.*armttbar[74] + 
   armttbar[146] + armttbar[197] + armttbar[243] + 1./4.*armttbar[239]
    + armttbar[147];
   armttbar[146]= - 5./12.*armttbar[79] - 23./24.*armttbar[88] - 5./24.
   *armttbar[85] + 1./4.*armttbar[146] - 7./3.*armttbar[44];
   armttbar[146]=1./24.*armttbar[126] + 11./27.*armttbar[36] + 13./108.
   *armttbar[38] - 3./8.*armttbar[37] + 13./64.*armttbar[13] - 3./8.*
   armttbar[40] - 4./9.*armttbar[25] + 1./32.*armttbar[23] - 277./288.*
   armttbar[84] + 139./144.*armttbar[51] + 1./4.*armttbar[146] + 
   armttbar[265];
   armttbar[239]= - 13./6.*armttbar[39];
   armttbar[244]= - 1./2.*armttbar[21];
   armttbar[249]=armttbar[244] + armttbar[107] + armttbar[235] + 1./24.
   *armttbar[41] + armttbar[37] + armttbar[239] - 365./36. + 
   armttbar[40];
   armttbar[249]=armttbar[8]*armttbar[249];
   armttbar[252]=7./2.*armttbar[92] - armttbar[48];
   armttbar[253]=1./3.*armttbar[18];
   armttbar[252]=armttbar[253] + 1./3.*armttbar[252] + armttbar[125];
   armttbar[252]= - 7./12.*armttbar[14] + 1./2.*armttbar[252] - 1./3.*
   armttbar[20];
   armttbar[252]=5./2.*armttbar[99]*armttbar[252];
   armttbar[254]=2*armttbar[56] - armttbar[53];
   armttbar[254]=2./3.*armttbar[254] + armttbar[277];
   armttbar[256]= - armttbar[12]*armttbar[36];
   armttbar[222]=1./6.*armttbar[256] + armttbar[222] - 1./6.*
   armttbar[13] + 7./3.*armttbar[84] + 1./6.*armttbar[79] + 31./16. - 1.
   /3.*armttbar[88];
   armttbar[222]=armttbar[99]*armttbar[222];
   armttbar[222]=1./3.*armttbar[254] + 5./16.*armttbar[222];
   armttbar[222]=MMH*armttbar[222];
   armttbar[254]= - 83./2. - 7./3.*armttbar[38];
   armttbar[254]=1./2.*armttbar[254] - 13./3.*armttbar[36];
   armttbar[254]=1./9.*armttbar[254] - 21./16.*armttbar[8];
   armttbar[254]=1./4.*armttbar[12]*armttbar[254];
   armttbar[256]= - 17./4. + armttbar[36];
   armttbar[256]=5./12.*armttbar[15]*armttbar[99]*armttbar[256];
   armttbar[257]= - 7./8. - armttbar[36];
   armttbar[257]=5./6.*armttbar[17]*armttbar[99]*armttbar[257];
   armttbar[249]=armttbar[222] + armttbar[257] + armttbar[256] + 
   armttbar[252] + armttbar[254] + armttbar[146] + 1./4.*armttbar[249];
   armttbar[249]=MMH*armttbar[249];
   armttbar[258]=349./144.*armttbar[21];
   armttbar[259]=109./18.*armttbar[8];
   armttbar[260]=13./12.*armttbar[12];
   armttbar[261]=23./4.*armttbar[98] - armttbar[100];
   armttbar[261]=3*armttbar[16]*armttbar[261];
   armttbar[262]=armttbar[261] + armttbar[260] + armttbar[259] + 
   armttbar[258] - armttbar[36] - armttbar[38] + armttbar[117] - 
   armttbar[37] + armttbar[127] + 157./12. - armttbar[40];
   armttbar[265]=armttbar[193] - 9095./54.*armttbar[99];
   armttbar[265]=armttbar[15]*armttbar[265];
   armttbar[262]=1./2.*armttbar[265] + 1./2.*armttbar[262] + 235./9.*
   armttbar[206];
   armttbar[265]=armttbar[203] - 10./9.*armttbar[14] + 9./8.*
   armttbar[16];
   armttbar[265]=armttbar[3]*armttbar[265];
   armttbar[262]=755./216.*armttbar[216] + 1./2.*armttbar[262] + 
   armttbar[265];
   armttbar[262]=armttbar[17]*armttbar[262];
   armttbar[266]= - 5935./648. + armttbar[144];
   armttbar[266]=armttbar[118] + 1./2.*armttbar[266] + armttbar[187];
   armttbar[267]=armttbar[193] - 3325./27.*armttbar[99];
   armttbar[267]=armttbar[15]*armttbar[267];
   armttbar[274]=1./12.*armttbar[8];
   armttbar[275]=1./4.*armttbar[205];
   armttbar[266]=1./2.*armttbar[267] + armttbar[275] - 739./324.*
   armttbar[12] + armttbar[274] + 169./288.*armttbar[21] + 
   armttbar[217] + 1./2.*armttbar[266] + armttbar[38];
   armttbar[266]=armttbar[15]*armttbar[266];
   armttbar[267]=11./16. + armttbar[44];
   armttbar[230]=armttbar[12]*armttbar[230];
   armttbar[167]=1./18.*armttbar[230] + armttbar[204] + 1./4.*
   armttbar[269] - 11./48.*armttbar[84] + 1./3.*armttbar[267] + 
   armttbar[167];
   armttbar[167]=MMH*armttbar[167];
   armttbar[204]=1./3.*armttbar[95] - 1./4.*armttbar[49];
   armttbar[230]=19./24.*armttbar[19];
   armttbar[267]=1./72.*armttbar[183] + 1./6.*armttbar[152] - 1./12.*
   armttbar[14] - 25./9.*armttbar[20] - 19./72.*armttbar[18] + 
   armttbar[230] + 1./3.*armttbar[33] + 25./72.*armttbar[32] + 1./12.*
   armttbar[48] + armttbar[204] - 13./36.*armttbar[92];
   armttbar[269]=1./12.*armttbar[21];
   armttbar[276]=5./16.*armttbar[8];
   armttbar[277]=armttbar[276] + armttbar[269] - 5./16. + armttbar[272]
   ;
   armttbar[278]=armttbar[277] - 1./54.*armttbar[12];
   armttbar[278]=armttbar[16]*armttbar[278];
   armttbar[237]= - 97./4.*armttbar[8] + armttbar[240] - 67./4. + 
   armttbar[237];
   armttbar[237]=1./12.*armttbar[237] + armttbar[12];
   armttbar[237]=armttbar[17]*armttbar[237];
   armttbar[240]=5./108.*armttbar[12] - 5./12.*armttbar[8] + 5./27.*
   armttbar[36] + 1./2. - 17./27.*armttbar[38];
   armttbar[240]=armttbar[15]*armttbar[240];
   armttbar[167]=1./3.*armttbar[167] + 1./3.*armttbar[237] + 1./4.*
   armttbar[240] + 1./2.*armttbar[267] + armttbar[278];
   armttbar[167]=MMH*armttbar[167];
   armttbar[237]= - armttbar[14] + 121./9.*armttbar[16];
   armttbar[237]=armttbar[16]*armttbar[237];
   armttbar[240]=8893./36.*armttbar[15] + armttbar[137] - 209./3.*
   armttbar[16];
   armttbar[240]=armttbar[15]*armttbar[240];
   armttbar[237]=1./4.*armttbar[237] + 1./3.*armttbar[240];
   armttbar[240]=5*armttbar[14];
   armttbar[267]=armttbar[240] - 5971./6.*armttbar[16];
   armttbar[267]=101./6.*armttbar[17] + 1./4.*armttbar[267] + 3335./9.*
   armttbar[15];
   armttbar[267]=armttbar[17]*armttbar[267];
   armttbar[237]=1./2.*armttbar[237] + 1./3.*armttbar[267];
   armttbar[237]=1./6.*armttbar[237] + armttbar[167];
   armttbar[267]=1./3.*armttbar[6]*MMH*armttbar[255];
   armttbar[237]=1./2.*armttbar[237] + armttbar[267];
   armttbar[237]=armttbar[6]*armttbar[237];
   armttbar[278]=745./24. - armttbar[40];
   armttbar[117]= - armttbar[36] - armttbar[38] + armttbar[117] + 
   armttbar[129] + 1./2.*armttbar[278] + armttbar[127];
   armttbar[117]=armttbar[14]*armttbar[117];
   armttbar[278]= - 3*armttbar[14];
   armttbar[280]=armttbar[278] - 5./4.*armttbar[16];
   armttbar[280]=armttbar[16]*armttbar[280];
   armttbar[236]= - 4./9.*armttbar[236];
   armttbar[281]=9./2.*armttbar[15] + armttbar[278] + 17./4.*
   armttbar[16];
   armttbar[281]=1./4.*armttbar[15]*armttbar[281];
   armttbar[280]=armttbar[281] + armttbar[236] + 1./4.*armttbar[280];
   armttbar[280]=armttbar[3]*armttbar[280];
   armttbar[282]= - 1./24.*armttbar[95] + armttbar[49];
   armttbar[148]=armttbar[148] - 13./12.*armttbar[92] + armttbar[268]
    + 139./48.*armttbar[96] + 7./6.*armttbar[94] + armttbar[282] - 11./
   16.*armttbar[50];
   armttbar[148]= - 763./288.*armttbar[19] - 83./144.*armttbar[42] - 1./
   24.*armttbar[5] - 31./96.*armttbar[33] - 341./576.*armttbar[34] + 1./
   2.*armttbar[148] - 8./9.*armttbar[32];
   armttbar[268]=2137./576.*armttbar[20];
   armttbar[283]=3./64.*armttbar[152];
   armttbar[228]=2./9.*armttbar[228];
   armttbar[284]=17./48.*armttbar[156];
   armttbar[106]=armttbar[106] + armttbar[112] + armttbar[237] + 
   armttbar[249] + armttbar[262] + armttbar[280] + 1./4.*armttbar[266]
    + 235./36.*armttbar[163] + 1./2.*armttbar[169] + armttbar[284] + 
   armttbar[228] + armttbar[283] + 1./4.*armttbar[117] + armttbar[268]
    + armttbar[148] - 1649./1152.*armttbar[18];
   armttbar[106]=armttbar[176]*armttbar[106];
   armttbar[112]= - 9./2.*armttbar[73] + 101./96. - 13*armttbar[83];
   armttbar[112]= - 35./16.*armttbar[87] - 1./4.*armttbar[72] - 1./4.*
   armttbar[74] + armttbar[116] + armttbar[77] - 39./16.*armttbar[90]
    + 1./4.*armttbar[141] - 15./16.*armttbar[89] - 1./2.*armttbar[80]
    + 5./32.*armttbar[76] + 1./4.*armttbar[112] + armttbar[132];
   armttbar[102]= - 55*armttbar[41] - 617./12. + armttbar[102];
   armttbar[102]=armttbar[219] + armttbar[246] - armttbar[36] + 1./4.*
   armttbar[102] - armttbar[38];
   armttbar[102]=armttbar[21]*armttbar[102];
   armttbar[116]=55./2.*armttbar[39];
   armttbar[117]= - 55./2.*armttbar[41];
   armttbar[132]=armttbar[117] + armttbar[129] + armttbar[116] - 51./4.
    - armttbar[40];
   armttbar[132]=armttbar[138] - armttbar[36] + 1./2.*armttbar[132] - 
   armttbar[38];
   armttbar[132]=armttbar[8]*armttbar[132];
   armttbar[169]=1./2.*armttbar[183] + 1./2.*armttbar[181] + 3*
   armttbar[42] - armttbar[14];
   armttbar[149]=3./2.*armttbar[155] + 1./2.*armttbar[169] + 
   armttbar[149];
   armttbar[149]=armttbar[3]*armttbar[149];
   armttbar[155]=7 + 17./4.*armttbar[21];
   armttbar[139]=armttbar[139] + 1./8.*armttbar[155] + armttbar[8];
   armttbar[139]=armttbar[12]*armttbar[139];
   armttbar[155]=1./2.*armttbar[25];
   armttbar[169]= - armttbar[64] - 3*armttbar[61] - 1./2.*armttbar[63];
   armttbar[169]=1./2.*armttbar[169] - armttbar[60];
   armttbar[169]=MMH*armttbar[169];
   armttbar[181]=1./2.*armttbar[68];
   armttbar[183]=armttbar[181] - armttbar[65];
   armttbar[183]=armttbar[63] + 1./2.*armttbar[183] + armttbar[61];
   armttbar[183]=MMt*armttbar[183];
   armttbar[102]=1./2.*armttbar[183] + 1./2.*armttbar[169] + 
   armttbar[149] + armttbar[195] + 1./2.*armttbar[161] + 1./2.*
   armttbar[139] + armttbar[132] + 1./4.*armttbar[102] + 1./8.*
   armttbar[166] + armttbar[212] + 1./4.*armttbar[115] + 3./64.*
   armttbar[13] + armttbar[157] + armttbar[155] + 13./8.*armttbar[23]
    + armttbar[123] - 21./16.*armttbar[51] + 1./8.*armttbar[88] + 3./8.
   *armttbar[85] + 1./2.*armttbar[112] + armttbar[44];
   armttbar[102]=MMt*armttbar[102];
   armttbar[112]= - 85./24. + armttbar[238];
   armttbar[112]= - 9./8.*armttbar[72] - 9./8.*armttbar[74] - 3./2.*
   armttbar[77] + armttbar[197] + armttbar[243] + 1./4.*armttbar[112]
    + armttbar[147];
   armttbar[115]=1./4.*armttbar[23];
   armttbar[112]=armttbar[118] + 5./16.*armttbar[13] + armttbar[187] + 
   armttbar[115] - 23./24.*armttbar[84] + 51./16.*armttbar[51] + 
   armttbar[119] + 3./8.*armttbar[88] - 5./12.*armttbar[85] + 1./2.*
   armttbar[112] - armttbar[44];
   armttbar[119]=armttbar[244] + armttbar[107] + armttbar[235] + 55./8.
   *armttbar[41] + armttbar[199] - 55./8.*armttbar[39] + 55./16. + 
   armttbar[40];
   armttbar[119]=armttbar[8]*armttbar[119];
   armttbar[123]=1./2.*armttbar[140] + armttbar[307] + armttbar[173] - 
   7*armttbar[84] - 1./2.*armttbar[79] - 93./16. + armttbar[88];
   armttbar[123]=MMH*armttbar[99]*armttbar[123];
   armttbar[132]= - 13./8.*armttbar[8];
   armttbar[139]=1./3.*armttbar[174] + armttbar[132];
   armttbar[139]=armttbar[12]*armttbar[139];
   armttbar[134]= - armttbar[18] + armttbar[134] + 3./4.*armttbar[32];
   armttbar[134]=7./4.*armttbar[14] + 1./2.*armttbar[134] + 
   armttbar[20];
   armttbar[134]=armttbar[99]*armttbar[134];
   armttbar[112]=1./16.*armttbar[123] + 1./2.*armttbar[182] + 1./4.*
   armttbar[171] + 1./2.*armttbar[134] + 1./8.*armttbar[139] + 1./2.*
   armttbar[119] + 1./12.*armttbar[126] + 5./12.*armttbar[36] + 1./4.*
   armttbar[112] + armttbar[272];
   armttbar[112]=MMH*armttbar[112];
   armttbar[119]= - 199./3. + armttbar[144];
   armttbar[123]= - 55./4.*armttbar[41];
   armttbar[119]= - 71./8.*armttbar[12] + armttbar[227] - 133./24.*
   armttbar[21] + armttbar[120] - armttbar[36] + armttbar[226] + 
   armttbar[123] + armttbar[118] + armttbar[220] + 1./2.*armttbar[119]
    + armttbar[187];
   armttbar[119]=1./4.*armttbar[119] + armttbar[234];
   armttbar[119]=armttbar[16]*armttbar[119];
   armttbar[126]= - armttbar[101] + 9./2.*armttbar[99];
   armttbar[126]=armttbar[15]*armttbar[126];
   armttbar[123]=3./2.*armttbar[126] + 33*armttbar[206] + armttbar[261]
    + 7./8.*armttbar[12] + armttbar[159] + 17./32.*armttbar[21] - 
   armttbar[36] - armttbar[38] + armttbar[123] + armttbar[129] + 
   armttbar[220] - 85./8. - armttbar[40];
   armttbar[126]=armttbar[203] - armttbar[14] + 9./4.*armttbar[16];
   armttbar[126]=armttbar[3]*armttbar[126];
   armttbar[134]=9./8.*armttbar[216];
   armttbar[123]=armttbar[134] + 1./2.*armttbar[123] + armttbar[126];
   armttbar[123]=armttbar[17]*armttbar[123];
   armttbar[126]=armttbar[8]*armttbar[207];
   armttbar[139]=armttbar[12]*armttbar[288];
   armttbar[126]=1./6.*armttbar[139] + 1./2.*armttbar[126] + 
   armttbar[214] + 1./6.*armttbar[304];
   armttbar[126]=MMH*armttbar[126];
   armttbar[140]=armttbar[223] + armttbar[214] + armttbar[202];
   armttbar[140]=armttbar[16]*armttbar[140];
   armttbar[147]=armttbar[184] + armttbar[214] + armttbar[244];
   armttbar[147]=armttbar[15]*armttbar[147];
   armttbar[139]=armttbar[304] + 1./2.*armttbar[139];
   armttbar[139]=MMH*armttbar[139];
   armttbar[139]=1./4.*armttbar[139] + 1./2.*armttbar[294] + 
   armttbar[140] + 1./2.*armttbar[147];
   armttbar[139]=MMH*armttbar[139];
   armttbar[140]=armttbar[16] + armttbar[15];
   armttbar[147]=armttbar[15]*armttbar[140];
   armttbar[147]=1./2.*armttbar[296] - armttbar[186] + 1./2.*
   armttbar[147];
   armttbar[139]=17./4.*armttbar[147] + armttbar[139];
   armttbar[139]=armttbar[6]*armttbar[139];
   armttbar[147]=armttbar[179] - 7 + armttbar[164];
   armttbar[147]=armttbar[36] + 1./4.*armttbar[147] - armttbar[38];
   armttbar[147]=armttbar[211] + 3*armttbar[147] + 17./12.*armttbar[21]
   ;
   armttbar[147]=armttbar[16]*armttbar[147];
   armttbar[149]=armttbar[3]*armttbar[297];
   armttbar[157]=armttbar[14]*armttbar[198];
   armttbar[161]=armttbar[209] - 17./24.*armttbar[21] + armttbar[36] + 
   21./8. - armttbar[38];
   armttbar[161]=armttbar[15]*armttbar[161];
   armttbar[126]=1./3.*armttbar[139] + armttbar[126] + armttbar[218] + 
   17./2.*armttbar[149] + armttbar[161] + armttbar[157] + 1./2.*
   armttbar[147];
   armttbar[139]=armttbar[8]*armttbar[198];
   armttbar[139]=1./4.*armttbar[208] + armttbar[139];
   armttbar[139]=MMt*armttbar[139];
   armttbar[126]=1./2.*armttbar[126] + armttbar[139];
   armttbar[126]=armttbar[210]*armttbar[126];
   armttbar[139]=armttbar[221] + 1./4.*armttbar[152] - 5./8.*
   armttbar[14] - 1./6.*armttbar[20] + 19./48.*armttbar[18] + 
   armttbar[230] - 3./16.*armttbar[32] - 1./8.*armttbar[48] + 
   armttbar[204] + 13./24.*armttbar[92];
   armttbar[147]=armttbar[277] + 1./48.*armttbar[12];
   armttbar[147]=armttbar[16]*armttbar[147];
   armttbar[149]=1./12.*armttbar[12] + armttbar[276] + armttbar[269] - 
   7./8. + armttbar[263];
   armttbar[149]=armttbar[15]*armttbar[149];
   armttbar[111]=armttbar[111] + 5./4. + 1./3.*armttbar[85];
   armttbar[157]= - armttbar[38] + armttbar[307];
   armttbar[157]=armttbar[21]*armttbar[157];
   armttbar[161]=armttbar[12]*armttbar[180];
   armttbar[111]=1./6.*armttbar[161] + 1./2.*armttbar[111] + 1./3.*
   armttbar[157];
   armttbar[111]=MMH*armttbar[111];
   armttbar[157]= - 13./16. - armttbar[36];
   armttbar[157]=1./3.*armttbar[157] - 5./16.*armttbar[8];
   armttbar[157]=armttbar[17]*armttbar[157];
   armttbar[111]=1./4.*armttbar[111] + 1./2.*armttbar[157] + 1./2.*
   armttbar[149] + 1./2.*armttbar[139] + armttbar[147];
   armttbar[111]=MMH*armttbar[111];
   armttbar[139]= - armttbar[14] - 197*armttbar[16];
   armttbar[139]=armttbar[16]*armttbar[139];
   armttbar[147]= - 85*armttbar[15] - armttbar[14] - 13*armttbar[16];
   armttbar[147]=armttbar[15]*armttbar[147];
   armttbar[139]=armttbar[139] + 1./2.*armttbar[147];
   armttbar[147]= - 5./3.*armttbar[14] - 423./2.*armttbar[16];
   armttbar[147]=5./2.*armttbar[17] + 1./4.*armttbar[147] + 
   armttbar[15];
   armttbar[147]=armttbar[17]*armttbar[147];
   armttbar[139]=1./12.*armttbar[139] + armttbar[147];
   armttbar[111]=1./4.*armttbar[139] + armttbar[111];
   armttbar[111]=armttbar[6]*armttbar[111];
   armttbar[139]= - 3113./24. + armttbar[144];
   armttbar[139]=armttbar[118] + 1./2.*armttbar[139] + armttbar[187];
   armttbar[147]= - armttbar[101] + 3*armttbar[99];
   armttbar[147]=armttbar[15]*armttbar[147];
   armttbar[132]=3./4.*armttbar[147] + 3./8.*armttbar[205] + 23./24.*
   armttbar[12] + armttbar[132] - 205./192.*armttbar[21] + 1./4.*
   armttbar[139] - armttbar[36];
   armttbar[132]=armttbar[15]*armttbar[132];
   armttbar[116]=armttbar[117] + armttbar[129] + armttbar[116] - 305./
   24. - armttbar[40];
   armttbar[116]= - armttbar[36] + 1./2.*armttbar[116] - armttbar[38];
   armttbar[116]=armttbar[14]*armttbar[116];
   armttbar[117]= - armttbar[14] - 23./4.*armttbar[16];
   armttbar[117]=armttbar[16]*armttbar[117];
   armttbar[139]=25./2.*armttbar[15] + armttbar[278] + 13./4.*
   armttbar[16];
   armttbar[139]=armttbar[15]*armttbar[139];
   armttbar[117]=3*armttbar[117] + 1./2.*armttbar[139];
   armttbar[117]=armttbar[3]*armttbar[117];
   armttbar[102]=armttbar[126] + armttbar[102] + armttbar[111] + 
   armttbar[112] + armttbar[123] + 1./2.*armttbar[117] + 1./2.*
   armttbar[132] + 33./4.*armttbar[163] + armttbar[119] + 3./16.*
   armttbar[156] + 3./32.*armttbar[152] + 1./2.*armttbar[116] + 27./64.
   *armttbar[20] - 1901./384.*armttbar[18] - 169./32.*armttbar[19] - 
   149./96.*armttbar[42] - 1./8.*armttbar[5] + 3./32.*armttbar[33] - 23.
   /64.*armttbar[34] + armttbar[125] + 3./8.*armttbar[48] - 1./8.*
   armttbar[92] - 1./8.*armttbar[93] + 247./32.*armttbar[96] + 
   armttbar[94] + armttbar[282] - 15./8.*armttbar[50];
   armttbar[102]=armttbar[210]*armttbar[102];
   armttbar[111]=367./9. - armttbar[40];
   armttbar[112]=47./12.*armttbar[41];
   armttbar[111]=armttbar[138] - armttbar[36] - armttbar[38] + 
   armttbar[112] + armttbar[129] + 1./2.*armttbar[111] + armttbar[127];
   armttbar[111]=armttbar[8]*armttbar[111];
   armttbar[116]=47./4.*armttbar[41] + 473./12. + armttbar[143];
   armttbar[116]=armttbar[185] + armttbar[246] - armttbar[36] + 1./3.*
   armttbar[116] - armttbar[38];
   armttbar[116]=armttbar[21]*armttbar[116];
   armttbar[110]=armttbar[172] + armttbar[225] + armttbar[160] + 
   armttbar[195] + armttbar[200] + armttbar[191] + 1./2.*armttbar[111]
    + armttbar[110] + 1./8.*armttbar[116];
   armttbar[110]=MMt*armttbar[110];
   armttbar[111]=3635./54. + armttbar[144];
   armttbar[111]= - 5149./108.*armttbar[12] + armttbar[227] + 77./8.*
   armttbar[21] + armttbar[120] - armttbar[36] + armttbar[226] + 
   armttbar[112] + armttbar[118] + armttbar[127] + 1./2.*armttbar[111]
    + armttbar[187];
   armttbar[111]=1./4.*armttbar[111] + armttbar[234];
   armttbar[111]=armttbar[16]*armttbar[111];
   armttbar[107]=armttbar[244] + armttbar[107] + armttbar[235] - 47./24.
   *armttbar[41] + armttbar[37] + armttbar[239] - 437./36. + 
   armttbar[40];
   armttbar[107]=armttbar[8]*armttbar[107];
   armttbar[107]=armttbar[222] + armttbar[257] + armttbar[256] + 
   armttbar[252] + armttbar[254] + armttbar[146] + 1./4.*armttbar[107];
   armttbar[107]=MMH*armttbar[107];
   armttbar[116]=armttbar[261] + armttbar[260] + armttbar[259] + 
   armttbar[258] - armttbar[36] - armttbar[38] + armttbar[112] - 
   armttbar[37] + armttbar[127] + 205./12. - armttbar[40];
   armttbar[117]= - armttbar[99]*armttbar[16];
   armttbar[119]=armttbar[193] - 271./6.*armttbar[99];
   armttbar[119]=armttbar[15]*armttbar[119];
   armttbar[116]=1./2.*armttbar[119] + 1./2.*armttbar[116] + 77*
   armttbar[117];
   armttbar[116]=armttbar[134] + 1./2.*armttbar[116] + armttbar[265];
   armttbar[116]=armttbar[17]*armttbar[116];
   armttbar[117]= - 6101./216. + armttbar[144];
   armttbar[117]=armttbar[118] + 1./2.*armttbar[117] + armttbar[187];
   armttbar[118]=armttbar[193] - 85./3.*armttbar[99];
   armttbar[118]=armttbar[15]*armttbar[118];
   armttbar[117]=1./2.*armttbar[118] + armttbar[275] + 95./108.*
   armttbar[12] + armttbar[274] + 211./96.*armttbar[21] + armttbar[217]
    + 1./2.*armttbar[117] + armttbar[38];
   armttbar[117]=armttbar[15]*armttbar[117];
   armttbar[118]= - armttbar[14] + 563./3.*armttbar[16];
   armttbar[118]=armttbar[16]*armttbar[118];
   armttbar[119]=575./12.*armttbar[15] + armttbar[137] - 737./3.*
   armttbar[16];
   armttbar[119]=armttbar[15]*armttbar[119];
   armttbar[118]=1./4.*armttbar[118] + 1./3.*armttbar[119];
   armttbar[119]=armttbar[240] - 211./6.*armttbar[16];
   armttbar[119]=1./4.*armttbar[119] + 269./3.*armttbar[15];
   armttbar[120]= - 1./2.*armttbar[17];
   armttbar[119]=1./9.*armttbar[119] + armttbar[120];
   armttbar[119]=armttbar[17]*armttbar[119];
   armttbar[118]=1./6.*armttbar[118] + armttbar[119];
   armttbar[118]=1./2.*armttbar[118] + armttbar[167];
   armttbar[118]=1./2.*armttbar[118] + armttbar[267];
   armttbar[118]=armttbar[6]*armttbar[118];
   armttbar[119]=937./24. - armttbar[40];
   armttbar[112]= - armttbar[36] - armttbar[38] + armttbar[112] + 
   armttbar[129] + 1./2.*armttbar[119] + armttbar[127];
   armttbar[112]=armttbar[14]*armttbar[112];
   armttbar[119]=armttbar[278] + 91./4.*armttbar[16];
   armttbar[119]=armttbar[16]*armttbar[119];
   armttbar[119]=armttbar[281] + armttbar[236] + 1./4.*armttbar[119];
   armttbar[119]=armttbar[3]*armttbar[119];
   armttbar[123]= - armttbar[99]*armttbar[186];
   armttbar[102]=1./2.*armttbar[102] + armttbar[110] + armttbar[118] + 
   armttbar[107] + armttbar[116] + armttbar[119] + 1./4.*armttbar[117]
    + 77./4.*armttbar[123] + 1./2.*armttbar[111] + armttbar[284] + 
   armttbar[228] + armttbar[283] + 1./4.*armttbar[112] + armttbar[268]
    + armttbar[148] - 161./1152.*armttbar[18];
   armttbar[102]=armttbar[210]*armttbar[102];
   armttbar[107]=armttbar[15]*armttbar[99];
   armttbar[110]=armttbar[17]*armttbar[99];
   armttbar[110]=8./9.*armttbar[110] + armttbar[206] + 8./3.*
   armttbar[107];
   armttbar[110]=armttbar[17]*armttbar[110];
   armttbar[111]= - 1 + armttbar[12];
   armttbar[112]=1./3.*armttbar[111] + 4*armttbar[107];
   armttbar[112]=armttbar[15]*armttbar[112];
   armttbar[116]=armttbar[16]*armttbar[111];
   armttbar[116]=2./9.*armttbar[18] + 5*armttbar[116];
   armttbar[110]=4*armttbar[110] + 8./9.*armttbar[112] + 1./3.*
   armttbar[116] + 2*armttbar[163];
   armttbar[112]=armttbar[124] - 16./3.*armttbar[15];
   armttbar[112]=armttbar[15]*armttbar[112];
   armttbar[116]= - 4./3.*armttbar[17] + armttbar[16] - 11./3.*
   armttbar[15];
   armttbar[116]=armttbar[17]*armttbar[116];
   armttbar[112]=4*armttbar[116] + armttbar[186] + armttbar[112];
   armttbar[112]=armttbar[6]*armttbar[112];
   armttbar[110]=2*armttbar[110] + 1./3.*armttbar[112];
   armttbar[102]=armttbar[105] + armttbar[102] + 8./3.*armttbar[110] + 
   armttbar[106];
   armttbar[102]=armttbar[9]*armttbar[102];
   armttbar[105]= - 1./2.*armttbar[70];
   armttbar[106]= - 2./3.*armttbar[52];
   armttbar[110]= - 2*armttbar[69];
   armttbar[112]=2./3.*armttbar[67] + 2*armttbar[65] + armttbar[181] + 
   armttbar[110] + armttbar[105] + armttbar[106];
   armttbar[116]=pow(CW,2);
   armttbar[112]=armttbar[116]*armttbar[112];
   armttbar[105]=64./27.*armttbar[54] + armttbar[110] + armttbar[106]
    - 128./27.*armttbar[57] + armttbar[105];
   armttbar[106]=armttbar[130] - 4./9.*armttbar[12];
   armttbar[106]=armttbar[3]*armttbar[106];
   armttbar[105]=armttbar[106] + armttbar[175] + armttbar[168] + 1./3.*
   armttbar[112] + 17./18.*armttbar[67] + 5./2.*armttbar[65] - 128./81.
   *armttbar[66] + 1./3.*armttbar[105] + 3./4.*armttbar[68];
   armttbar[105]=MMZ*armttbar[105];
   armttbar[106]=20*armttbar[116];
   armttbar[110]=8*armttbar[116];
   armttbar[112]=49./3. + armttbar[110];
   armttbar[112]=armttbar[41]*armttbar[112];
   armttbar[117]= - 2*armttbar[21];
   armttbar[118]= - 32./9.*armttbar[12] + armttbar[117] + 2*
   armttbar[112] + 313./9. + armttbar[106];
   armttbar[118]=armttbar[12]*armttbar[118];
   armttbar[119]=4./3.*armttbar[116] + armttbar[90] + 2./3.*
   armttbar[80] - 38./9. - armttbar[81];
   armttbar[119]=armttbar[116]*armttbar[119];
   armttbar[123]=53./12. + armttbar[116];
   armttbar[123]=armttbar[41]*armttbar[123];
   armttbar[123]=armttbar[123] + armttbar[39] + 22./3. + armttbar[116];
   armttbar[123]=armttbar[21]*armttbar[123];
   armttbar[124]= - 701./27. - armttbar[81];
   armttbar[124]=2*armttbar[124] + 25./6.*armttbar[58];
   armttbar[125]= - 2./9.*armttbar[24];
   armttbar[126]= - 1 - 2./3.*armttbar[116];
   armttbar[126]=armttbar[13]*armttbar[126];
   armttbar[127]= - armttbar[3]*armttbar[15];
   armttbar[132]= - armttbar[17]*armttbar[3];
   armttbar[105]=armttbar[105] + 8./9.*armttbar[132] + 4./9.*
   armttbar[127] + 2./9.*armttbar[118] + armttbar[123] + 4./3.*
   armttbar[126] + 4./3.*armttbar[119] + armttbar[125] + 146./81.*
   armttbar[25] + 5./4.*armttbar[23] + 64./81.*armttbar[45] + 101./36.*
   armttbar[87] + 64./81.*armttbar[75] + 3*armttbar[90] - 29./24.*
   armttbar[141] - 28./27.*armttbar[80] + 2./9.*armttbar[76] + 1./3.*
   armttbar[124] - 5./4.*armttbar[83];
   armttbar[105]=MMZ*armttbar[105];
   armttbar[118]= - 203./6. - armttbar[58];
   armttbar[118]=139./384.*armttbar[141] + 1./3.*armttbar[118] + 
   armttbar[76];
   armttbar[119]=29./4.*armttbar[141] + 2 - armttbar[58];
   armttbar[119]=1./3.*armttbar[119] + armttbar[90];
   armttbar[119]=armttbar[116]*armttbar[119];
   armttbar[123]= - 1 - armttbar[116];
   armttbar[123]=1./9.*armttbar[123] + armttbar[177];
   armttbar[123]=armttbar[21]*armttbar[123];
   armttbar[118]= - 16./27.*armttbar[12] + armttbar[123] - 25./27.*
   armttbar[13] + 1./3.*armttbar[119] - 16./27.*armttbar[45] + 2./3.*
   armttbar[87] + 16./27.*armttbar[75] + 1./3.*armttbar[118] + 1./2.*
   armttbar[90];
   armttbar[118]=MMZ*armttbar[118];
   armttbar[119]=104./27. - armttbar[116];
   armttbar[119]=armttbar[20]*armttbar[119];
   armttbar[119]=2*armttbar[119] + 40./9.*armttbar[18] - 16./9.*
   armttbar[33] - 2./3.*armttbar[34] + armttbar[96] - 16./9.*
   armttbar[93];
   armttbar[123]=85 + 16*armttbar[12];
   armttbar[123]=armttbar[15]*armttbar[123];
   armttbar[124]= - 16./3.*armttbar[12];
   armttbar[126]=armttbar[124] - 26./3. + armttbar[21];
   armttbar[126]=armttbar[17]*armttbar[126];
   armttbar[119]=2./3.*armttbar[126] + 1./9.*armttbar[123] + 2*
   armttbar[119] - 5./3.*armttbar[16];
   armttbar[118]=1./3.*armttbar[119] + armttbar[118];
   armttbar[118]=MMZ*armttbar[118];
   armttbar[119]=armttbar[116]*armttbar[141];
   armttbar[119]=armttbar[141] + 1./2.*armttbar[119];
   armttbar[119]=armttbar[245]*armttbar[119];
   armttbar[119]= - 32./27.*armttbar[255] + 21./16.*armttbar[119];
   armttbar[119]=MMZ*armttbar[119];
   armttbar[123]= - armttbar[116]*armttbar[141];
   armttbar[123]= - armttbar[141] + 1./2.*armttbar[123];
   armttbar[123]=armttbar[116]*armttbar[123];
   armttbar[123]= - 3./2.*armttbar[141] + armttbar[123];
   armttbar[123]=armttbar[6]*armttbar[196]*armttbar[123];
   armttbar[119]=armttbar[119] + 5./16.*armttbar[123];
   armttbar[119]=armttbar[6]*armttbar[119];
   armttbar[123]= - 8./9.*armttbar[17] + armttbar[16] - 80./27.*
   armttbar[15];
   armttbar[123]=armttbar[17]*armttbar[123];
   armttbar[123]= - 8./9.*armttbar[229] + armttbar[123];
   armttbar[118]=armttbar[119] + 4./3.*armttbar[123] + armttbar[118];
   armttbar[118]=armttbar[6]*armttbar[118];
   armttbar[119]= - 23./3.*armttbar[50] - 2*armttbar[94];
   armttbar[123]=armttbar[116]*armttbar[50];
   armttbar[126]= - 4*armttbar[116];
   armttbar[127]=29./3. + armttbar[126];
   armttbar[127]=armttbar[19]*armttbar[127];
   armttbar[119]=2./3.*armttbar[127] + 4./3.*armttbar[123] - 128./27.*
   armttbar[33] - 2*armttbar[34] - 88./9.*armttbar[93] + 1./3.*
   armttbar[119] + 4*armttbar[96];
   armttbar[112]=160./9.*armttbar[12] + 5./4.*armttbar[21] + 4*
   armttbar[112] + 2009./12. + 40*armttbar[116];
   armttbar[112]=armttbar[15]*armttbar[112];
   armttbar[127]= - 16*armttbar[116];
   armttbar[134]= - 49./3. - 8*armttbar[116];
   armttbar[134]=armttbar[41]*armttbar[134];
   armttbar[134]=2*armttbar[134] - 265./3. + armttbar[127];
   armttbar[134]=32./27.*armttbar[12] + 1./3.*armttbar[134] + 2*
   armttbar[21];
   armttbar[134]=armttbar[17]*armttbar[134];
   armttbar[137]=69./4. - 8./9.*armttbar[116];
   armttbar[137]=armttbar[18]*armttbar[137];
   armttbar[138]=217./9. + armttbar[127];
   armttbar[138]=armttbar[20]*armttbar[138];
   armttbar[139]= - 37./3. + armttbar[127];
   armttbar[139]=40./9.*armttbar[12] + armttbar[21] + 1./9.*
   armttbar[139] + armttbar[41];
   armttbar[139]=armttbar[16]*armttbar[139];
   armttbar[105]=armttbar[118] + armttbar[105] + 2./3.*armttbar[134] + 
   1./9.*armttbar[112] + armttbar[139] + 4./9.*armttbar[138] + 2./3.*
   armttbar[119] + armttbar[137];
   armttbar[105]=armttbar[6]*armttbar[105];
   armttbar[112]=4*armttbar[116] - armttbar[90] + 23./3. + armttbar[80]
   ;
   armttbar[112]=armttbar[116]*armttbar[112];
   armttbar[118]= - 64./3.*armttbar[75] - 10*armttbar[90] + 104./3.*
   armttbar[80] + 4211./18. - 4*armttbar[76];
   armttbar[119]=11./9. + armttbar[126];
   armttbar[119]=armttbar[13]*armttbar[119];
   armttbar[134]=13./3. + 2*armttbar[116];
   armttbar[134]=armttbar[41]*armttbar[134];
   armttbar[137]= - 2*armttbar[116];
   armttbar[138]= - 13./3. + armttbar[137];
   armttbar[138]=armttbar[41]*armttbar[138];
   armttbar[139]= - 16./9. + armttbar[138];
   armttbar[139]=8./3.*armttbar[12] + 2*armttbar[139] + armttbar[21];
   armttbar[139]=armttbar[12]*armttbar[139];
   armttbar[112]=4./3.*armttbar[139] - 4./3.*armttbar[21] + 8./3.*
   armttbar[134] + 4./3.*armttbar[119] + 16./3.*armttbar[112] + 1./3.*
   armttbar[118] - 8*armttbar[87];
   armttbar[112]=armttbar[99]*armttbar[112];
   armttbar[113]=1./8. + armttbar[113];
   armttbar[113]=armttbar[151]*armttbar[113];
   armttbar[109]=armttbar[153]*armttbar[109];
   armttbar[109]=1./4.*armttbar[113] + 3*armttbar[109];
   armttbar[109]=MMZ*armttbar[109];
   armttbar[113]= - 3*armttbar[39] + 37./9. + armttbar[137];
   armttbar[118]= - 53./3. + armttbar[126];
   armttbar[118]=armttbar[41]*armttbar[118];
   armttbar[113]= - 8./9.*armttbar[12] + 2*armttbar[113] + 
   armttbar[118];
   armttbar[113]=armttbar[3]*armttbar[113];
   armttbar[109]=armttbar[109] + armttbar[113] + armttbar[112] + 
   armttbar[178] + 1./3.*armttbar[65] - 1./4.*armttbar[68] - 1./3.*
   armttbar[69] + armttbar[224] - 256./81.*armttbar[57] + armttbar[215]
   ;
   armttbar[109]=MMZ*armttbar[109];
   armttbar[112]= - 8./3. - armttbar[116];
   armttbar[112]=armttbar[19]*armttbar[112];
   armttbar[112]=4*armttbar[112] + 2*armttbar[123] + 16./9.*
   armttbar[33] + armttbar[34] + 16./3.*armttbar[93] - 2*armttbar[96]
    + 13./3.*armttbar[50] + armttbar[94];
   armttbar[113]= - 29 + armttbar[126];
   armttbar[113]=armttbar[18]*armttbar[113];
   armttbar[118]= - 29./3. + armttbar[126];
   armttbar[118]=armttbar[16]*armttbar[118];
   armttbar[112]=2*armttbar[118] - 260./9.*armttbar[20] + 2*
   armttbar[112] + armttbar[113];
   armttbar[112]=armttbar[99]*armttbar[112];
   armttbar[110]=armttbar[117] + 4*armttbar[134] + 341./9. + 
   armttbar[110];
   armttbar[110]=armttbar[17]*armttbar[99]*armttbar[110];
   armttbar[113]=4./3.*armttbar[80] + 1 - 4./3.*armttbar[81];
   armttbar[113]=armttbar[116]*armttbar[113];
   armttbar[117]= - 181./9. + armttbar[127];
   armttbar[117]=armttbar[13]*armttbar[117];
   armttbar[113]=2./3.*armttbar[117] + 8*armttbar[113] + 8./27.*
   armttbar[25] + 16./27.*armttbar[45] - 32./27.*armttbar[75] + 110./9.
   *armttbar[80] + 5 - 98./9.*armttbar[81];
   armttbar[117]= - 65./3. + armttbar[127];
   armttbar[117]=armttbar[41]*armttbar[117];
   armttbar[113]=2*armttbar[113] + 1./3.*armttbar[117];
   armttbar[116]=5./3. + armttbar[116];
   armttbar[116]=armttbar[41]*armttbar[116];
   armttbar[106]=16./9.*armttbar[12] + 4*armttbar[116] + 241./9. + 
   armttbar[106];
   armttbar[106]=armttbar[12]*armttbar[106];
   armttbar[116]= - 1 - armttbar[41];
   armttbar[116]=armttbar[21]*armttbar[116];
   armttbar[117]= - 16./9.*armttbar[12] + armttbar[21] - 259./9. + 2*
   armttbar[138];
   armttbar[117]=armttbar[15]*armttbar[99]*armttbar[117];
   armttbar[118]= - 5 - 3*armttbar[41];
   armttbar[118]=armttbar[3]*armttbar[16]*armttbar[118];
   armttbar[119]=armttbar[66] + armttbar[55] - 2*armttbar[54];
   armttbar[119]=MMt*armttbar[119];
   armttbar[102]=armttbar[102] + armttbar[103] + armttbar[104] + 128./
   81.*armttbar[119] + armttbar[105] + armttbar[109] + 8./3.*
   armttbar[110] + 2*armttbar[118] + 8./3.*armttbar[117] + 8./3.*
   armttbar[112] + 4./9.*armttbar[106] + 1./3.*armttbar[113] + 1./2.*
   armttbar[116];
   armttbar[102]=armttbar[46]*armttbar[102];
   armttbar[103]=3./2.*armttbar[20] + armttbar[162] + armttbar[19] + 
   armttbar[250] - armttbar[34] - 1./4.*armttbar[33];
   armttbar[104]=7./2.*armttbar[10];
   armttbar[105]= - 3*armttbar[30];
   armttbar[106]=armttbar[104] + armttbar[105] + 3 + 11./2.*
   armttbar[29];
   armttbar[106]=armttbar[15]*armttbar[106];
   armttbar[109]=armttbar[10] + 3./2. + 2*armttbar[29];
   armttbar[109]=armttbar[16]*armttbar[109];
   armttbar[103]=1./4.*armttbar[106] + 3*armttbar[103] + armttbar[109];
   armttbar[103]=armttbar[3]*armttbar[103];
   armttbar[106]= - 1 + armttbar[29];
   armttbar[109]=3./4.*armttbar[106] + armttbar[10];
   armttbar[109]=armttbar[12]*armttbar[109];
   armttbar[110]= - 3./4.*armttbar[29];
   armttbar[109]=armttbar[109] - armttbar[10] + armttbar[110] + 
   armttbar[13] - armttbar[24] + 3./2. - armttbar[25];
   armttbar[109]=armttbar[99]*armttbar[109];
   armttbar[112]=25./4.*armttbar[29];
   armttbar[113]=17./4.*armttbar[10];
   armttbar[116]=armttbar[113] - 11 + armttbar[112];
   armttbar[116]=armttbar[3]*armttbar[116];
   armttbar[109]=3./16.*armttbar[109] + armttbar[116];
   armttbar[109]=MMZ*armttbar[109];
   armttbar[116]=1 - armttbar[29];
   armttbar[117]=armttbar[232] + 3./4.*armttbar[116] - armttbar[10];
   armttbar[117]=armttbar[99]*armttbar[117];
   armttbar[118]=15./2.*armttbar[12] + 31 + 15*armttbar[21];
   armttbar[118]=armttbar[3]*armttbar[118];
   armttbar[117]=3./2.*armttbar[117] + armttbar[118];
   armttbar[117]=armttbar[17]*armttbar[117];
   armttbar[118]=1./8.*armttbar[13] - 1./8.*armttbar[24] + 9./8.*
   armttbar[25] + 1 - 3./8.*armttbar[45];
   armttbar[118]= - 9./8.*armttbar[10] + 5./2.*armttbar[30] + 3*
   armttbar[118] - 65./32.*armttbar[29];
   armttbar[119]= - 1./8.*armttbar[30];
   armttbar[123]=1./12.*armttbar[10] + armttbar[119] - 2./3. - 3./8.*
   armttbar[29];
   armttbar[123]=armttbar[21]*armttbar[123];
   armttbar[126]= - 239 - 13*armttbar[29];
   armttbar[126]=11./12.*armttbar[10] + 1./48.*armttbar[126] - 
   armttbar[30];
   armttbar[126]=armttbar[12]*armttbar[126];
   armttbar[127]=1./4.*armttbar[33];
   armttbar[134]= - 5./2.*armttbar[20] + 3./4.*armttbar[18] + 
   armttbar[127] - armttbar[5];
   armttbar[134]=armttbar[99]*armttbar[134];
   armttbar[137]=armttbar[10] - 1 + 3./4.*armttbar[29];
   armttbar[137]=armttbar[15]*armttbar[99]*armttbar[137];
   armttbar[103]=1./2.*armttbar[109] + 1./4.*armttbar[117] + 
   armttbar[103] + 3./16.*armttbar[137] + 3./16.*armttbar[134] + 1./8.*
   armttbar[126] + 1./4.*armttbar[118] + armttbar[123];
   armttbar[103]=armttbar[27]*armttbar[103];
   armttbar[109]=armttbar[173] - 1./2.*armttbar[24] + armttbar[155] + 3.
   /8. - armttbar[22];
   armttbar[117]= - 13./2.*armttbar[10] + 3*armttbar[109] + 5*
   armttbar[11];
   armttbar[117]=armttbar[31]*armttbar[117];
   armttbar[109]= - 13./6.*armttbar[10] + armttbar[109] + 5./3.*
   armttbar[11];
   armttbar[109]=armttbar[1]*armttbar[109];
   armttbar[109]=armttbar[117] + armttbar[109];
   armttbar[117]=armttbar[10] + 1./2. - armttbar[11];
   armttbar[118]=armttbar[31]*armttbar[117];
   armttbar[117]=armttbar[1]*armttbar[117];
   armttbar[117]=armttbar[118] + 1./3.*armttbar[117];
   armttbar[117]=armttbar[21]*armttbar[117];
   armttbar[118]=1 - armttbar[11];
   armttbar[123]=armttbar[118] + 5./2.*armttbar[10];
   armttbar[126]=armttbar[31]*armttbar[123];
   armttbar[123]=armttbar[1]*armttbar[123];
   armttbar[123]=armttbar[126] + 1./3.*armttbar[123];
   armttbar[123]=armttbar[12]*armttbar[123];
   armttbar[126]= - armttbar[5] + armttbar[18];
   armttbar[134]=1./2.*armttbar[126] - armttbar[20];
   armttbar[137]=armttbar[31]*armttbar[134];
   armttbar[134]=armttbar[1]*armttbar[134];
   armttbar[138]=3*armttbar[137] + armttbar[134];
   armttbar[138]=armttbar[99]*armttbar[138];
   armttbar[139]=armttbar[31]*armttbar[10];
   armttbar[141]=armttbar[1]*armttbar[10];
   armttbar[143]=3*armttbar[139] + armttbar[141];
   armttbar[144]=armttbar[15]*armttbar[99]*armttbar[143];
   armttbar[109]=1./2.*armttbar[144] + armttbar[138] + 1./2.*
   armttbar[123] + 1./2.*armttbar[109] + armttbar[117];
   armttbar[117]=armttbar[145] + armttbar[133] + armttbar[26] - 
   armttbar[23];
   armttbar[123]=armttbar[31]*armttbar[117];
   armttbar[117]=armttbar[1]*armttbar[117];
   armttbar[133]= - 1 - 3*armttbar[10];
   armttbar[133]=armttbar[31]*armttbar[133];
   armttbar[138]= - 1./3. - armttbar[10];
   armttbar[144]=armttbar[1]*armttbar[138];
   armttbar[133]=armttbar[133] + armttbar[144];
   armttbar[133]=armttbar[21]*armttbar[133];
   armttbar[117]=armttbar[133] + 3*armttbar[123] + armttbar[117];
   armttbar[123]= - 1 + armttbar[11];
   armttbar[123]=1./4.*armttbar[123] - armttbar[10];
   armttbar[133]=armttbar[31]*armttbar[123];
   armttbar[123]=armttbar[1]*armttbar[123];
   armttbar[123]=armttbar[133] + 1./3.*armttbar[123];
   armttbar[133]=armttbar[12]*armttbar[123];
   armttbar[117]=1./2.*armttbar[117] + armttbar[133];
   armttbar[117]=MMZ*armttbar[117];
   armttbar[133]=1./2.*armttbar[20] + armttbar[293] - armttbar[19] + 
   armttbar[7] + armttbar[233];
   armttbar[145]=armttbar[31]*armttbar[133];
   armttbar[133]=armttbar[1]*armttbar[133];
   armttbar[133]=3*armttbar[145] + armttbar[133];
   armttbar[145]= - 1 - armttbar[11];
   armttbar[145]=1./2.*armttbar[145] - armttbar[10];
   armttbar[146]=armttbar[31]*armttbar[145];
   armttbar[145]=armttbar[1]*armttbar[145];
   armttbar[145]=armttbar[146] + 1./3.*armttbar[145];
   armttbar[145]=armttbar[16]*armttbar[145];
   armttbar[123]=armttbar[15]*armttbar[123];
   armttbar[146]= - 7./2. - armttbar[11];
   armttbar[146]=1./4.*armttbar[146] + armttbar[10];
   armttbar[147]=armttbar[31]*armttbar[146];
   armttbar[146]=armttbar[1]*armttbar[146];
   armttbar[146]=armttbar[147] + 1./3.*armttbar[146];
   armttbar[146]=armttbar[17]*armttbar[146];
   armttbar[117]=armttbar[117] + armttbar[146] + armttbar[123] + 1./2.*
   armttbar[133] + armttbar[145];
   armttbar[123]= - 1 + 1./6.*armttbar[29];
   armttbar[133]=5./6.*armttbar[10];
   armttbar[145]= - 1./2.*armttbar[30];
   armttbar[123]=armttbar[248] - 5./3.*armttbar[21] + armttbar[133] + 7
   *armttbar[123] + armttbar[145];
   armttbar[123]=armttbar[17]*armttbar[123];
   armttbar[146]= - 1./3.*armttbar[29];
   armttbar[147]= - 1./6.*armttbar[10];
   armttbar[148]= - 1./4.*armttbar[30];
   armttbar[149]=armttbar[147] + armttbar[148] - 1./4. + armttbar[146];
   armttbar[149]=armttbar[16]*armttbar[149];
   armttbar[146]=1./2. + armttbar[146];
   armttbar[151]= - 5./3.*armttbar[10];
   armttbar[146]=armttbar[151] + 7*armttbar[146] + armttbar[30];
   armttbar[146]=armttbar[12]*armttbar[146];
   armttbar[128]=armttbar[128] + 1./4.*armttbar[24] - armttbar[23] + 1./
   2.*armttbar[45] + 1./4. + armttbar[43];
   armttbar[152]= - 7./6.*armttbar[10] - 1 - 11./6.*armttbar[29];
   armttbar[152]=armttbar[21]*armttbar[152];
   armttbar[128]=1./2.*armttbar[146] + 3*armttbar[128] + armttbar[152];
   armttbar[128]=MMZ*armttbar[128];
   armttbar[127]= - 1./4.*armttbar[20] - armttbar[19] + 1./4.*
   armttbar[5] + armttbar[34] + armttbar[127];
   armttbar[146]=armttbar[151] + armttbar[30] - 1 - 7./3.*armttbar[29];
   armttbar[146]=armttbar[15]*armttbar[146];
   armttbar[123]=1./4.*armttbar[128] + 1./4.*armttbar[123] + 1./8.*
   armttbar[146] + 3./4.*armttbar[127] + armttbar[149];
   armttbar[123]=armttbar[27]*armttbar[123];
   armttbar[127]=1 - armttbar[12];
   armttbar[127]=armttbar[17]*armttbar[127];
   armttbar[128]= - armttbar[12] - 1 - armttbar[45];
   armttbar[128]=MMZ*armttbar[128];
   armttbar[146]=1./2.*armttbar[128] + armttbar[127] + armttbar[188] - 
   1./2.*armttbar[33] + armttbar[20];
   armttbar[146]=MMZ*armttbar[146];
   armttbar[149]= - armttbar[15] + 3./4.*armttbar[17];
   armttbar[149]=armttbar[17]*armttbar[149];
   armttbar[149]=armttbar[149] + 1./2.*armttbar[146];
   armttbar[149]=armttbar[6]*armttbar[27]*armttbar[149];
   armttbar[117]=3./8.*armttbar[149] + 1./2.*armttbar[117] + 
   armttbar[123];
   armttbar[117]=armttbar[6]*armttbar[117];
   armttbar[123]=1./2.*armttbar[29] - armttbar[30];
   armttbar[149]=armttbar[123] + 1./2.*armttbar[10];
   armttbar[152]=armttbar[21]*armttbar[149];
   armttbar[155]=armttbar[12]*armttbar[149];
   armttbar[152]=armttbar[152] + 1./2.*armttbar[155];
   armttbar[152]=MMZ*armttbar[152];
   armttbar[155]=armttbar[16]*armttbar[149];
   armttbar[156]=armttbar[15]*armttbar[149];
   armttbar[157]= - 1./2.*armttbar[10];
   armttbar[160]=armttbar[157] - 1./2.*armttbar[29] + armttbar[30];
   armttbar[161]=armttbar[17]*armttbar[160];
   armttbar[152]=armttbar[152] + 1./2.*armttbar[161] + armttbar[155] + 
   1./2.*armttbar[156];
   armttbar[152]=armttbar[27]*armttbar[152];
   armttbar[155]= - armttbar[11] + armttbar[10];
   armttbar[156]=armttbar[31]*armttbar[155];
   armttbar[155]=armttbar[1]*armttbar[155];
   armttbar[161]=armttbar[156] + 1./3.*armttbar[155];
   armttbar[163]=armttbar[21]*armttbar[161];
   armttbar[164]=armttbar[12]*armttbar[161];
   armttbar[163]=armttbar[163] + 1./2.*armttbar[164];
   armttbar[163]=MMZ*armttbar[163];
   armttbar[164]=armttbar[16]*armttbar[161];
   armttbar[166]=armttbar[15]*armttbar[161];
   armttbar[167]=armttbar[11] - armttbar[10];
   armttbar[168]=armttbar[31]*armttbar[167];
   armttbar[167]=armttbar[1]*armttbar[167];
   armttbar[169]=armttbar[168] + 1./3.*armttbar[167];
   armttbar[171]=armttbar[17]*armttbar[169];
   armttbar[152]=armttbar[152] + armttbar[163] + 1./2.*armttbar[171] + 
   armttbar[164] + 1./2.*armttbar[166];
   armttbar[152]=armttbar[6]*armttbar[152];
   armttbar[163]=armttbar[12]*armttbar[169];
   armttbar[164]=armttbar[21]*armttbar[169];
   armttbar[163]=armttbar[163] + armttbar[164] + 3*armttbar[156] + 
   armttbar[155];
   armttbar[166]=3*armttbar[168] + armttbar[167];
   armttbar[172]=armttbar[16]*armttbar[166];
   armttbar[173]=armttbar[15]*armttbar[166];
   armttbar[172]=armttbar[172] + 1./2.*armttbar[173];
   armttbar[172]=armttbar[3]*armttbar[172];
   armttbar[173]= - 1./4.*armttbar[10];
   armttbar[110]=armttbar[173] + armttbar[110] + armttbar[30];
   armttbar[174]=armttbar[21]*armttbar[110];
   armttbar[175]=3./2.*armttbar[30];
   armttbar[177]=armttbar[157] - armttbar[29] + armttbar[175];
   armttbar[177]=armttbar[12]*armttbar[177];
   armttbar[174]=1./2.*armttbar[177] + 3./2.*armttbar[149] + 
   armttbar[174];
   armttbar[177]=armttbar[16]*armttbar[160];
   armttbar[178]=armttbar[15]*armttbar[160];
   armttbar[177]=armttbar[177] + 1./2.*armttbar[178];
   armttbar[177]=armttbar[3]*armttbar[177];
   armttbar[178]=MMZ*armttbar[3]*armttbar[160];
   armttbar[174]=3*armttbar[178] + 1./2.*armttbar[174] + 3*
   armttbar[177];
   armttbar[174]=armttbar[27]*armttbar[174];
   armttbar[166]=MMZ*armttbar[3]*armttbar[166];
   armttbar[177]=armttbar[29] - armttbar[30];
   armttbar[178]=MMt*armttbar[27]*armttbar[3]*armttbar[177];
   armttbar[152]=3./2.*armttbar[178] + 1./2.*armttbar[152] + 
   armttbar[174] + armttbar[166] + 1./4.*armttbar[163] + armttbar[172];
   armttbar[152]=armttbar[210]*armttbar[152];
   armttbar[163]=3*armttbar[10];
   armttbar[166]=armttbar[118] + armttbar[163];
   armttbar[172]=armttbar[31]*armttbar[166];
   armttbar[166]=armttbar[1]*armttbar[166];
   armttbar[166]=3*armttbar[172] + armttbar[166];
   armttbar[166]=armttbar[15]*armttbar[166];
   armttbar[154]=armttbar[162] + armttbar[19] - armttbar[7] + 
   armttbar[154];
   armttbar[162]=armttbar[31]*armttbar[154];
   armttbar[154]=armttbar[1]*armttbar[154];
   armttbar[172]=1./2. + armttbar[10];
   armttbar[174]=armttbar[31]*armttbar[172];
   armttbar[172]=armttbar[1]*armttbar[172];
   armttbar[178]=3*armttbar[174] + armttbar[172];
   armttbar[178]=armttbar[16]*armttbar[178];
   armttbar[154]=1./4.*armttbar[166] + armttbar[178] + 3*armttbar[162]
    + armttbar[154];
   armttbar[154]=armttbar[3]*armttbar[154];
   armttbar[143]=armttbar[12]*armttbar[143];
   armttbar[162]= - armttbar[10] + armttbar[13] - armttbar[25] - 
   armttbar[24];
   armttbar[166]=armttbar[31]*armttbar[162];
   armttbar[162]=armttbar[1]*armttbar[162];
   armttbar[143]=armttbar[143] + 3*armttbar[166] + armttbar[162];
   armttbar[143]=armttbar[99]*armttbar[143];
   armttbar[178]= - 11 + 21./2.*armttbar[10];
   armttbar[178]=armttbar[31]*armttbar[178];
   armttbar[104]= - 11./3. + armttbar[104];
   armttbar[179]=armttbar[1]*armttbar[104];
   armttbar[178]=armttbar[178] + armttbar[179];
   armttbar[178]=armttbar[3]*armttbar[178];
   armttbar[143]=1./8.*armttbar[143] + armttbar[178];
   armttbar[143]=MMZ*armttbar[143];
   armttbar[178]= - 9*armttbar[28];
   armttbar[105]=armttbar[105] + 23./2.*armttbar[29] + 35./2. + 
   armttbar[178];
   armttbar[105]=3./2.*armttbar[12] + 1./2.*armttbar[105] + 
   armttbar[231];
   armttbar[105]=armttbar[3]*armttbar[105];
   armttbar[180]=armttbar[199] - 1./2. + armttbar[40];
   armttbar[180]=armttbar[17]*armttbar[153]*armttbar[180];
   armttbar[105]=1./2.*armttbar[105] + 9*armttbar[180];
   armttbar[105]=armttbar[27]*armttbar[105];
   armttbar[180]=MMt*armttbar[27]*armttbar[153]*armttbar[28];
   armttbar[105]=armttbar[105] + 9*armttbar[180];
   armttbar[105]=MMt*armttbar[105];
   armttbar[181]=1 - armttbar[10];
   armttbar[182]=armttbar[31]*armttbar[181];
   armttbar[181]=armttbar[1]*armttbar[181];
   armttbar[183]=3*armttbar[182] + armttbar[181];
   armttbar[183]=armttbar[17]*armttbar[99]*armttbar[183];
   armttbar[103]=1./2.*armttbar[152] + armttbar[105] + armttbar[117] + 
   armttbar[103] + 1./2.*armttbar[143] + 1./4.*armttbar[183] + 1./4.*
   armttbar[109] + armttbar[154];
   armttbar[103]=armttbar[210]*armttbar[103];
   armttbar[105]=armttbar[126] - armttbar[20];
   armttbar[109]=armttbar[31]*armttbar[105];
   armttbar[105]=armttbar[1]*armttbar[105];
   armttbar[117]=armttbar[109] + 1./3.*armttbar[105];
   armttbar[143]= - 1./3. + armttbar[10];
   armttbar[152]=armttbar[31]*armttbar[143];
   armttbar[143]=armttbar[1]*armttbar[143];
   armttbar[143]=5./3.*armttbar[152] + armttbar[143];
   armttbar[152]=armttbar[16]*armttbar[143];
   armttbar[154]=17./4.*armttbar[11];
   armttbar[183]=13./3. + armttbar[154];
   armttbar[113]=1./3.*armttbar[183] + armttbar[113];
   armttbar[113]=armttbar[31]*armttbar[113];
   armttbar[183]=67./12.*armttbar[10] + 1 + 17./12.*armttbar[11];
   armttbar[183]=armttbar[1]*armttbar[183];
   armttbar[113]=armttbar[113] + 1./3.*armttbar[183];
   armttbar[183]=armttbar[15]*armttbar[113];
   armttbar[184]= - 17./4.*armttbar[11];
   armttbar[185]=41./3. + armttbar[184];
   armttbar[185]=1./3.*armttbar[185] - 17./4.*armttbar[10];
   armttbar[185]=armttbar[31]*armttbar[185];
   armttbar[186]= - 67./12.*armttbar[10] + 5 - 17./12.*armttbar[11];
   armttbar[186]=armttbar[1]*armttbar[186];
   armttbar[185]=armttbar[185] + 1./3.*armttbar[186];
   armttbar[185]=armttbar[17]*armttbar[185];
   armttbar[117]=1./2.*armttbar[185] + 1./2.*armttbar[183] + 2*
   armttbar[117] + armttbar[152];
   armttbar[152]= - armttbar[26] + armttbar[23];
   armttbar[125]=2./9.*armttbar[13] + 1./4.*armttbar[152] + 
   armttbar[125];
   armttbar[125]=armttbar[1]*armttbar[125];
   armttbar[183]=7./27. + armttbar[11];
   armttbar[183]=1./2.*armttbar[183] + 19./9.*armttbar[10];
   armttbar[183]=armttbar[31]*armttbar[183];
   armttbar[185]= - 1./3. + armttbar[11];
   armttbar[185]=1./6.*armttbar[185] + armttbar[10];
   armttbar[185]=armttbar[1]*armttbar[185];
   armttbar[183]=armttbar[183] + armttbar[185];
   armttbar[183]=armttbar[21]*armttbar[183];
   armttbar[113]=armttbar[12]*armttbar[113];
   armttbar[152]=2./3.*armttbar[13] + 3./4.*armttbar[152] - 2./3.*
   armttbar[24];
   armttbar[152]=armttbar[31]*armttbar[152];
   armttbar[113]=1./6.*armttbar[113] + 1./2.*armttbar[183] + 
   armttbar[152] + armttbar[125];
   armttbar[113]=MMZ*armttbar[113];
   armttbar[125]=1./12. + armttbar[29];
   armttbar[152]=1./4.*armttbar[30];
   armttbar[125]=5./18.*armttbar[10] + 7./9.*armttbar[125] + 
   armttbar[152];
   armttbar[125]=armttbar[21]*armttbar[125];
   armttbar[183]=19./18.*armttbar[13] - 19./18.*armttbar[24] + 3*
   armttbar[23] - armttbar[45] + 11./18. - 3*armttbar[43];
   armttbar[185]=17./3.*armttbar[30];
   armttbar[186]=armttbar[185] - 89./9. + 25./2.*armttbar[29];
   armttbar[186]=1./3.*armttbar[186] + 3./2.*armttbar[10];
   armttbar[186]=armttbar[12]*armttbar[186];
   armttbar[125]=1./8.*armttbar[186] + 1./4.*armttbar[183] + 
   armttbar[125];
   armttbar[125]=MMZ*armttbar[125];
   armttbar[183]=11*armttbar[33];
   armttbar[186]=armttbar[183] - 19*armttbar[5];
   armttbar[188]= - 5*armttbar[18];
   armttbar[186]= - 41./2.*armttbar[20] + 1./2.*armttbar[186] + 
   armttbar[188];
   armttbar[191]=4*armttbar[29];
   armttbar[193]=armttbar[10] - 5./3. + armttbar[191];
   armttbar[195]=armttbar[16]*armttbar[193];
   armttbar[186]=1./4.*armttbar[186] + armttbar[195];
   armttbar[195]= - 25./2.*armttbar[29];
   armttbar[196]= - 17./3.*armttbar[30];
   armttbar[197]=armttbar[196] + 11./9. + armttbar[195];
   armttbar[197]=1./3.*armttbar[197] - 3./2.*armttbar[10];
   armttbar[197]=19./12.*armttbar[12] + 1./8.*armttbar[197] + 14./9.*
   armttbar[21];
   armttbar[197]=armttbar[17]*armttbar[197];
   armttbar[198]=17./12.*armttbar[30] - 17./9. + 25./8.*armttbar[29];
   armttbar[198]=1./3.*armttbar[198] + 3./8.*armttbar[10];
   armttbar[198]=armttbar[15]*armttbar[198];
   armttbar[125]=armttbar[125] + armttbar[197] + 1./9.*armttbar[186] + 
   1./2.*armttbar[198];
   armttbar[125]=armttbar[27]*armttbar[125];
   armttbar[111]=armttbar[17]*armttbar[111];
   armttbar[186]=armttbar[12] + 1 + armttbar[45];
   armttbar[186]=MMZ*armttbar[186];
   armttbar[197]=1./2.*armttbar[33];
   armttbar[111]=1./2.*armttbar[186] + armttbar[111] + armttbar[270] + 
   armttbar[197] - armttbar[20];
   armttbar[111]=MMZ*armttbar[111];
   armttbar[186]=armttbar[241] - 1./4.*armttbar[17];
   armttbar[186]=armttbar[17]*armttbar[186];
   armttbar[111]=armttbar[186] + 1./6.*armttbar[111];
   armttbar[111]=armttbar[6]*armttbar[27]*armttbar[111];
   armttbar[113]=29./12.*armttbar[111] + armttbar[125] + 1./3.*
   armttbar[117] + armttbar[113];
   armttbar[113]=armttbar[6]*armttbar[113];
   armttbar[117]= - 11*armttbar[45];
   armttbar[125]=323./12.*armttbar[29] + 49./3.*armttbar[13] - 49./3.*
   armttbar[24] + 25./3.*armttbar[25] - 31./2. + armttbar[117];
   armttbar[186]=1./3.*armttbar[30];
   armttbar[125]=1./2.*armttbar[125] + armttbar[186];
   armttbar[198]= - 1./3.*armttbar[10];
   armttbar[125]=1./2.*armttbar[125] + armttbar[198];
   armttbar[199]=armttbar[147] + 89./18. + armttbar[191];
   armttbar[199]=armttbar[21]*armttbar[199];
   armttbar[200]=2809./9. + 187*armttbar[29];
   armttbar[200]=1./4.*armttbar[200] - armttbar[30];
   armttbar[200]=1./4.*armttbar[200] + armttbar[10];
   armttbar[200]=armttbar[12]*armttbar[200];
   armttbar[125]=1./4.*armttbar[200] + 1./4.*armttbar[125] + 
   armttbar[199];
   armttbar[199]=5./3. - 29./8.*armttbar[29];
   armttbar[200]= - 3./4.*armttbar[30];
   armttbar[199]=7./24.*armttbar[10] + 1./3.*armttbar[199] + 
   armttbar[200];
   armttbar[199]=armttbar[15]*armttbar[199];
   armttbar[202]=armttbar[33] - armttbar[5];
   armttbar[202]=1./2.*armttbar[202] - armttbar[20];
   armttbar[203]= - 4*armttbar[29];
   armttbar[204]= - armttbar[10] + 5./3. + armttbar[203];
   armttbar[205]=armttbar[16]*armttbar[204];
   armttbar[205]=armttbar[199] + armttbar[202] + 2./3.*armttbar[205];
   armttbar[205]=armttbar[3]*armttbar[205];
   armttbar[206]=5 - 13./4.*armttbar[29];
   armttbar[207]= - 7*armttbar[10];
   armttbar[206]=3*armttbar[206] + armttbar[207];
   armttbar[206]=armttbar[15]*armttbar[99]*armttbar[206];
   armttbar[208]= - 7 + 39*armttbar[29];
   armttbar[209]=7*armttbar[10];
   armttbar[208]= - 39./4.*armttbar[12] + 1./4.*armttbar[208] + 
   armttbar[209];
   armttbar[208]=armttbar[99]*armttbar[208];
   armttbar[211]= - 5./2.*armttbar[12];
   armttbar[212]= - 53./3. + armttbar[211];
   armttbar[212]=armttbar[3]*armttbar[212];
   armttbar[208]=1./2.*armttbar[208] + armttbar[212];
   armttbar[208]=armttbar[17]*armttbar[208];
   armttbar[207]=39./4.*armttbar[116] + armttbar[207];
   armttbar[207]=armttbar[12]*armttbar[207];
   armttbar[207]=armttbar[207] + armttbar[209] + 39./4.*armttbar[29] - 
   7*armttbar[13] + 7*armttbar[24] - 51./2. + 7*armttbar[25];
   armttbar[207]=armttbar[99]*armttbar[207];
   armttbar[195]=43./3. + armttbar[195];
   armttbar[195]= - 2./3.*armttbar[10] + 1./3.*armttbar[195] - 2*
   armttbar[30];
   armttbar[195]=armttbar[3]*armttbar[195];
   armttbar[195]=1./16.*armttbar[207] + armttbar[195];
   armttbar[195]=MMZ*armttbar[195];
   armttbar[207]=7./2.*armttbar[20] - 1./4.*armttbar[18] - 3./4.*
   armttbar[33] + armttbar[5];
   armttbar[207]=armttbar[99]*armttbar[207];
   armttbar[125]=armttbar[195] + 1./2.*armttbar[208] + armttbar[205] + 
   1./8.*armttbar[206] + 1./3.*armttbar[125] + 7./8.*armttbar[207];
   armttbar[125]=armttbar[27]*armttbar[125];
   armttbar[195]=1./3. - armttbar[10];
   armttbar[205]=armttbar[31]*armttbar[195];
   armttbar[195]=armttbar[1]*armttbar[195];
   armttbar[195]=5./3.*armttbar[205] + armttbar[195];
   armttbar[205]=armttbar[21]*armttbar[195];
   armttbar[206]=7./4.*armttbar[11];
   armttbar[207]=41./3. + armttbar[206];
   armttbar[207]=1./3.*armttbar[207] + armttbar[173];
   armttbar[207]=armttbar[31]*armttbar[207];
   armttbar[208]= - 19./12.*armttbar[10] + 5 + 7./12.*armttbar[11];
   armttbar[208]=armttbar[1]*armttbar[208];
   armttbar[207]=armttbar[207] + 1./3.*armttbar[208];
   armttbar[207]=armttbar[12]*armttbar[207];
   armttbar[208]=1./24.*armttbar[11];
   armttbar[209]=35./24.*armttbar[10] + armttbar[208] + armttbar[13] - 
   armttbar[24] + 7./6. + armttbar[25];
   armttbar[209]=armttbar[31]*armttbar[209];
   armttbar[212]=83./24.*armttbar[10] + armttbar[208] + armttbar[13] - 
   armttbar[24] + 1./2. + armttbar[25];
   armttbar[212]=armttbar[1]*armttbar[212];
   armttbar[207]=1./2.*armttbar[207] + 1./2.*armttbar[205] + 
   armttbar[209] + 1./3.*armttbar[212];
   armttbar[209]=armttbar[5] - armttbar[18];
   armttbar[212]=2*armttbar[20];
   armttbar[214]=armttbar[209] + armttbar[212];
   armttbar[215]=armttbar[31]*armttbar[214];
   armttbar[214]=armttbar[1]*armttbar[214];
   armttbar[217]=armttbar[215] + 1./3.*armttbar[214];
   armttbar[217]=armttbar[99]*armttbar[217];
   armttbar[195]=armttbar[16]*armttbar[195];
   armttbar[218]= - 11./12.*armttbar[10] + 5./9. - 3./4.*armttbar[11];
   armttbar[218]=armttbar[31]*armttbar[218];
   armttbar[219]= - 1./4.*armttbar[11];
   armttbar[220]=1./3. + armttbar[219];
   armttbar[221]=armttbar[220] - 3./4.*armttbar[10];
   armttbar[221]=armttbar[1]*armttbar[221];
   armttbar[218]=armttbar[218] + armttbar[221];
   armttbar[218]=armttbar[15]*armttbar[218];
   armttbar[195]=2*armttbar[195] + armttbar[218];
   armttbar[195]=armttbar[3]*armttbar[195];
   armttbar[221]= - armttbar[31]*armttbar[10];
   armttbar[222]= - armttbar[1]*armttbar[10];
   armttbar[223]=armttbar[221] + 1./3.*armttbar[222];
   armttbar[224]=armttbar[12]*armttbar[223];
   armttbar[225]=armttbar[10] - armttbar[13] + armttbar[25] + 
   armttbar[24];
   armttbar[226]=armttbar[31]*armttbar[225];
   armttbar[225]=armttbar[1]*armttbar[225];
   armttbar[224]=armttbar[224] + armttbar[226] + 1./3.*armttbar[225];
   armttbar[224]=armttbar[99]*armttbar[224];
   armttbar[227]= - 2*armttbar[11];
   armttbar[228]=17./3. + armttbar[227];
   armttbar[229]= - 5./2.*armttbar[10];
   armttbar[228]=1./3.*armttbar[228] + armttbar[229];
   armttbar[228]=armttbar[1]*armttbar[228];
   armttbar[227]= - 29./6.*armttbar[10] + 43./9. + armttbar[227];
   armttbar[227]=armttbar[31]*armttbar[227];
   armttbar[227]=armttbar[227] + armttbar[228];
   armttbar[227]=armttbar[3]*armttbar[227];
   armttbar[224]=armttbar[224] + armttbar[227];
   armttbar[224]=MMZ*armttbar[224];
   armttbar[227]= - 3*armttbar[28];
   armttbar[228]=1./2.*armttbar[30];
   armttbar[230]= - armttbar[12] + armttbar[228] - 77./3.*armttbar[29]
    - 37./6. + armttbar[227];
   armttbar[230]=armttbar[3]*armttbar[230];
   armttbar[231]= - 1 + armttbar[264];
   armttbar[231]=3*armttbar[17]*armttbar[153]*armttbar[231];
   armttbar[230]=1./2.*armttbar[230] + armttbar[231];
   armttbar[230]=armttbar[27]*armttbar[230];
   armttbar[232]=6*armttbar[180];
   armttbar[230]=armttbar[230] + armttbar[232];
   armttbar[230]=MMt*armttbar[230];
   armttbar[223]=armttbar[15]*armttbar[99]*armttbar[223];
   armttbar[233]= - 1 + armttbar[10];
   armttbar[234]=armttbar[31]*armttbar[233];
   armttbar[233]=armttbar[1]*armttbar[233];
   armttbar[235]=armttbar[234] + 1./3.*armttbar[233];
   armttbar[235]=armttbar[17]*armttbar[99]*armttbar[235];
   armttbar[103]=armttbar[103] + armttbar[230] + armttbar[113] + 
   armttbar[125] + armttbar[224] + 4*armttbar[235] + armttbar[195] + 2*
   armttbar[223] + 1./3.*armttbar[207] + 2*armttbar[217];
   armttbar[103]=armttbar[210]*armttbar[103];
   armttbar[113]= - 13*armttbar[39];
   armttbar[125]=11./3.*armttbar[10] - 20./9. + 3*armttbar[11];
   armttbar[125]=armttbar[31]*armttbar[125];
   armttbar[163]=armttbar[163] - 4./3. + armttbar[11];
   armttbar[163]=armttbar[1]*armttbar[163];
   armttbar[195]= - 55./4.*armttbar[12];
   armttbar[140]=9./2.*armttbar[3]*armttbar[140];
   armttbar[207]= - 15./8.*armttbar[21];
   armttbar[217]=armttbar[140] + armttbar[195] + armttbar[207] + 
   armttbar[163] + armttbar[125] + armttbar[201] + armttbar[213] + 1./4.
   *armttbar[41] + armttbar[142] + armttbar[113] - 299./12. + 
   armttbar[187];
   armttbar[217]=armttbar[17]*armttbar[3]*armttbar[217];
   armttbar[223]= - 17*armttbar[29];
   armttbar[224]=121./3. + armttbar[223];
   armttbar[230]= - 5./9.*armttbar[10];
   armttbar[235]=3*armttbar[30];
   armttbar[224]=armttbar[230] + 1./9.*armttbar[224] + armttbar[235];
   armttbar[224]=armttbar[21]*armttbar[224];
   armttbar[183]= - 25*armttbar[20] - 11*armttbar[18] + armttbar[242]
    + 3*armttbar[34] + armttbar[183];
   armttbar[236]=armttbar[175] - 139./6.*armttbar[29] - 59./3. + 
   armttbar[178];
   armttbar[236]=armttbar[15]*armttbar[236];
   armttbar[237]= - 1 - armttbar[28];
   armttbar[237]=3*armttbar[237] - armttbar[30];
   armttbar[237]=armttbar[16]*armttbar[237];
   armttbar[237]=1./2.*armttbar[236] + armttbar[183] + 3./2.*
   armttbar[237];
   armttbar[237]=armttbar[3]*armttbar[237];
   armttbar[238]= - 20./3. + 17./2.*armttbar[29];
   armttbar[133]=armttbar[133] + 1./3.*armttbar[238] + armttbar[235];
   armttbar[133]=armttbar[27]*armttbar[17]*armttbar[3]*armttbar[133];
   armttbar[238]=13759./2. + 511*armttbar[29];
   armttbar[239]= - 5./2.*armttbar[30];
   armttbar[238]=1./27.*armttbar[238] + armttbar[239];
   armttbar[238]=armttbar[12]*armttbar[238];
   armttbar[240]=10./3. - 17./4.*armttbar[29];
   armttbar[240]= - 5./36.*armttbar[10] + 1./9.*armttbar[240] + 
   armttbar[145];
   armttbar[240]=armttbar[8]*armttbar[240];
   armttbar[241]= - 3./8.*armttbar[43];
   armttbar[242]= - 49./27. + armttbar[241];
   armttbar[243]=3./16.*armttbar[23];
   armttbar[245]= - 19./144.*armttbar[30];
   armttbar[217]=armttbar[133] + armttbar[217] + 1./2.*armttbar[237] + 
   1./24.*armttbar[238] + armttbar[240] + 1./16.*armttbar[224] + 
   armttbar[245] - 73./216.*armttbar[29] + 400./81.*armttbar[25] + 
   armttbar[243] + 1./2.*armttbar[242] + 400./81.*armttbar[45];
   armttbar[217]=armttbar[27]*armttbar[217];
   armttbar[224]= - 11./9.*armttbar[10] + 47./27. - armttbar[11];
   armttbar[224]=armttbar[31]*armttbar[224];
   armttbar[237]=7./3. - armttbar[11];
   armttbar[237]=1./3.*armttbar[237] - armttbar[10];
   armttbar[237]=armttbar[1]*armttbar[237];
   armttbar[224]=armttbar[224] + armttbar[237];
   armttbar[237]=armttbar[21]*armttbar[224];
   armttbar[238]= - armttbar[31]*armttbar[22];
   armttbar[242]= - armttbar[1]*armttbar[22];
   armttbar[238]=armttbar[238] + 1./3.*armttbar[242];
   armttbar[237]=armttbar[238] + armttbar[237];
   armttbar[242]= - 1./2.*armttbar[11];
   armttbar[246]=2./3. + armttbar[242];
   armttbar[246]=1./3.*armttbar[246] + armttbar[157];
   armttbar[246]=armttbar[1]*armttbar[246];
   armttbar[242]= - 11./18.*armttbar[10] + 10./27. + armttbar[242];
   armttbar[242]=armttbar[31]*armttbar[242];
   armttbar[242]=armttbar[242] + armttbar[246];
   armttbar[246]=armttbar[8]*armttbar[242];
   armttbar[237]=1./8.*armttbar[237] + armttbar[246];
   armttbar[246]=armttbar[16]*armttbar[28];
   armttbar[248]=armttbar[15]*armttbar[28];
   armttbar[249]=armttbar[246] + armttbar[248];
   armttbar[249]=armttbar[3]*armttbar[249];
   armttbar[250]= - 11*armttbar[12] - 3./2.*armttbar[21] + 133./3. + 
   armttbar[239];
   armttbar[249]=1./2.*armttbar[250] + 9*armttbar[249];
   armttbar[249]=MMt*armttbar[27]*armttbar[3]*armttbar[249];
   armttbar[217]=armttbar[249] + armttbar[237] + armttbar[217];
   armttbar[217]=MMt*armttbar[217];
   armttbar[250]=107./4. - 259*armttbar[29];
   armttbar[252]= - 19./4.*armttbar[30];
   armttbar[250]=armttbar[229] + 1./9.*armttbar[250] + armttbar[252];
   armttbar[250]=887./324.*armttbar[12] + 1./9.*armttbar[250] + 
   armttbar[279];
   armttbar[247]=armttbar[247] - 41./3.*armttbar[15];
   armttbar[247]=armttbar[3]*armttbar[247];
   armttbar[254]=2173./54.*armttbar[99] - 9*armttbar[3];
   armttbar[254]=armttbar[17]*armttbar[254];
   armttbar[256]= - armttbar[15]*armttbar[99];
   armttbar[247]=5./2.*armttbar[254] + 1./2.*armttbar[247] + 1./2.*
   armttbar[250] + 1685./27.*armttbar[256];
   armttbar[247]=armttbar[17]*armttbar[247];
   armttbar[250]=5./3. - 17./8.*armttbar[29];
   armttbar[250]=1./9.*armttbar[250];
   armttbar[254]= - 5./72.*armttbar[10];
   armttbar[257]=armttbar[254] + armttbar[250] + armttbar[152];
   armttbar[257]=armttbar[16]*armttbar[257];
   armttbar[250]=armttbar[254] + armttbar[250] + armttbar[148];
   armttbar[250]=armttbar[14]*armttbar[250];
   armttbar[258]=1663 + 4775./6.*armttbar[29];
   armttbar[259]= - 17*armttbar[30];
   armttbar[258]=1./3.*armttbar[258] + armttbar[259];
   armttbar[258]=armttbar[15]*armttbar[258];
   armttbar[253]= - armttbar[33] + armttbar[253];
   armttbar[253]=11./2.*armttbar[253];
   armttbar[260]=armttbar[253] + 7291./81.*armttbar[20];
   armttbar[247]=1./2.*armttbar[247] + 1./144.*armttbar[258] + 1./2.*
   armttbar[257] + 1./16.*armttbar[260] + armttbar[250];
   armttbar[247]=armttbar[27]*armttbar[247];
   armttbar[150]=3./2. + armttbar[150];
   armttbar[150]=9./2.*armttbar[170] + 35./12.*armttbar[12] + 1./2.*
   armttbar[150] + armttbar[142];
   armttbar[150]=armttbar[17]*armttbar[3]*armttbar[150];
   armttbar[257]=103 + 35./2.*armttbar[29];
   armttbar[257]=1./9.*armttbar[257] - armttbar[30];
   armttbar[258]= - 7*armttbar[29];
   armttbar[260]= - 13 + armttbar[258];
   armttbar[260]=7./144.*armttbar[260] - armttbar[30];
   armttbar[260]=armttbar[12]*armttbar[260];
   armttbar[257]=1./24.*armttbar[257] + armttbar[260];
   armttbar[260]= - 11./3. + armttbar[178];
   armttbar[260]=armttbar[175] + 1./2.*armttbar[260] + 7./3.*
   armttbar[29];
   armttbar[260]=armttbar[15]*armttbar[260];
   armttbar[261]= - armttbar[33] + armttbar[18];
   armttbar[261]=1./2.*armttbar[261] + armttbar[20];
   armttbar[260]=7./3.*armttbar[261] + 1./2.*armttbar[260];
   armttbar[260]=armttbar[3]*armttbar[260];
   armttbar[150]=armttbar[150] + 1./3.*armttbar[257] + armttbar[260];
   armttbar[150]=armttbar[27]*armttbar[150];
   armttbar[257]=7./3.*armttbar[12] + 11./3. + armttbar[30];
   armttbar[260]=armttbar[3]*armttbar[248];
   armttbar[257]=1./2.*armttbar[257] + 9*armttbar[260];
   armttbar[257]=MMt*armttbar[27]*armttbar[3]*armttbar[257];
   armttbar[150]=armttbar[150] + armttbar[257];
   armttbar[150]=MMt*armttbar[150];
   armttbar[257]=7*armttbar[29];
   armttbar[260]=323./27. + armttbar[257];
   armttbar[185]= - 193./18.*armttbar[12] + 1./2.*armttbar[260] + 
   armttbar[185];
   armttbar[170]=67./2.*armttbar[170] + 1./4.*armttbar[185] + 725./9.*
   armttbar[107];
   armttbar[185]= - 4775./54.*armttbar[99] - 17*armttbar[3];
   armttbar[185]=armttbar[17]*armttbar[185];
   armttbar[170]=1./3.*armttbar[170] + 1./2.*armttbar[185];
   armttbar[170]=armttbar[17]*armttbar[170];
   armttbar[185]= - 221 - 301./2.*armttbar[29];
   armttbar[185]=1./6.*armttbar[185] + armttbar[259];
   armttbar[185]=armttbar[15]*armttbar[185];
   armttbar[260]=armttbar[33] + 77./9.*armttbar[18];
   armttbar[260]=1./2.*armttbar[260] - armttbar[20];
   armttbar[185]=7./2.*armttbar[260] + 1./3.*armttbar[185];
   armttbar[170]=1./12.*armttbar[185] + armttbar[170];
   armttbar[170]=armttbar[27]*armttbar[170];
   armttbar[185]= - 29./2.*armttbar[15] + 19*armttbar[17];
   armttbar[185]=armttbar[6]*armttbar[27]*armttbar[17]*armttbar[185];
   armttbar[170]=armttbar[170] + 17./27.*armttbar[185];
   armttbar[185]=1 + armttbar[39];
   armttbar[185]=armttbar[17]*armttbar[3]*armttbar[185];
   armttbar[260]= - armttbar[12]*armttbar[30];
   armttbar[185]=17./36.*armttbar[260] + armttbar[185];
   armttbar[185]=armttbar[27]*armttbar[185];
   armttbar[261]=armttbar[3]*armttbar[30];
   armttbar[262]=MMt*armttbar[27]*armttbar[261];
   armttbar[185]=armttbar[185] + armttbar[262];
   armttbar[185]=MMt*armttbar[185];
   armttbar[262]= - armttbar[27]*armttbar[17]*armttbar[12];
   armttbar[185]=17./36.*armttbar[262] + armttbar[185];
   armttbar[185]=armttbar[176]*armttbar[185];
   armttbar[150]=1./2.*armttbar[185] + 1./2.*armttbar[170] + 
   armttbar[150];
   armttbar[150]=armttbar[176]*armttbar[150];
   armttbar[170]=1./3.*armttbar[220] + armttbar[173];
   armttbar[170]=armttbar[1]*armttbar[170];
   armttbar[173]= - 11./36.*armttbar[10] + 5./27. + armttbar[219];
   armttbar[173]=armttbar[31]*armttbar[173];
   armttbar[170]=armttbar[173] + armttbar[170];
   armttbar[173]=armttbar[14]*armttbar[170];
   armttbar[185]=armttbar[16]*armttbar[224];
   armttbar[170]=armttbar[17]*armttbar[170];
   armttbar[170]=armttbar[170] + armttbar[173] + 1./8.*armttbar[185];
   armttbar[173]=1./4.*armttbar[11];
   armttbar[185]=11./36.*armttbar[10] - 5./27. + armttbar[173];
   armttbar[185]=armttbar[31]*armttbar[185];
   armttbar[173]= - 1./3. + armttbar[173];
   armttbar[173]=1./3.*armttbar[173] + 1./4.*armttbar[10];
   armttbar[173]=armttbar[1]*armttbar[173];
   armttbar[173]=armttbar[185] + armttbar[173];
   armttbar[173]=armttbar[8]*armttbar[173];
   armttbar[185]= - 5./3. + 17./8.*armttbar[29];
   armttbar[185]=5./72.*armttbar[10] + 1./9.*armttbar[185] + 
   armttbar[152];
   armttbar[185]=armttbar[27]*armttbar[8]*armttbar[185];
   armttbar[173]=armttbar[173] + armttbar[185];
   armttbar[173]=1./2.*MMH*armttbar[173];
   armttbar[185]=armttbar[16] + 3049./81.*armttbar[15];
   armttbar[185]=1./2.*armttbar[185] - 2165./81.*armttbar[17];
   armttbar[185]=armttbar[6]*armttbar[27]*armttbar[17]*armttbar[185];
   armttbar[150]=1./2.*armttbar[150] + armttbar[217] + 1./2.*
   armttbar[185] + armttbar[173] + armttbar[170] + armttbar[247];
   armttbar[150]=armttbar[176]*armttbar[150];
   armttbar[113]=armttbar[140] + armttbar[195] + armttbar[207] + 
   armttbar[163] + armttbar[125] + armttbar[201] + armttbar[213] - 47./
   4.*armttbar[41] + armttbar[142] + armttbar[113] - 443./12. + 
   armttbar[187];
   armttbar[113]=armttbar[17]*armttbar[3]*armttbar[113];
   armttbar[125]= - 91./3. + armttbar[178];
   armttbar[125]=armttbar[200] + 1./4.*armttbar[125] - 16./3.*
   armttbar[29];
   armttbar[125]=armttbar[16]*armttbar[125];
   armttbar[125]=1./4.*armttbar[236] + 1./2.*armttbar[183] + 
   armttbar[125];
   armttbar[125]=armttbar[3]*armttbar[125];
   armttbar[140]=armttbar[230] + armttbar[235] - 71./27. - 9*
   armttbar[29];
   armttbar[140]=armttbar[21]*armttbar[140];
   armttbar[142]=7./3. + armttbar[241];
   armttbar[163]=armttbar[239] + 391./6. - 19*armttbar[29];
   armttbar[163]=armttbar[12]*armttbar[163];
   armttbar[113]=armttbar[133] + armttbar[113] + armttbar[125] + 1./24.
   *armttbar[163] + armttbar[240] + 1./16.*armttbar[140] + 
   armttbar[245] + 125./72.*armttbar[29] + 16./9.*armttbar[25] + 
   armttbar[243] + 1./2.*armttbar[142] + 16./9.*armttbar[45];
   armttbar[113]=armttbar[27]*armttbar[113];
   armttbar[113]=armttbar[249] + armttbar[237] + armttbar[113];
   armttbar[113]=MMt*armttbar[113];
   armttbar[110]=armttbar[16]*armttbar[110];
   armttbar[125]= - armttbar[29] + armttbar[235];
   armttbar[125]=armttbar[194] + armttbar[21] + 1./2.*armttbar[125] - 
   armttbar[10];
   armttbar[125]=armttbar[17]*armttbar[125];
   armttbar[133]=armttbar[14]*armttbar[160];
   armttbar[140]= - armttbar[29] + armttbar[30];
   armttbar[142]=armttbar[15]*armttbar[140];
   armttbar[110]=1./2.*armttbar[125] + 1./4.*armttbar[142] + 
   armttbar[133] + armttbar[110];
   armttbar[110]=armttbar[27]*armttbar[110];
   armttbar[125]=armttbar[14]*armttbar[169];
   armttbar[133]=armttbar[16]*armttbar[169];
   armttbar[142]=armttbar[8]*armttbar[161];
   armttbar[163]=armttbar[27]*armttbar[8]*armttbar[149];
   armttbar[163]=armttbar[142] + armttbar[163];
   armttbar[163]=MMH*armttbar[163];
   armttbar[110]=1./2.*armttbar[163] + armttbar[110] + armttbar[171] + 
   armttbar[125] + 1./2.*armttbar[133];
   armttbar[125]=armttbar[235] - armttbar[10];
   armttbar[125]=armttbar[21]*armttbar[125];
   armttbar[125]=3*armttbar[140] + armttbar[125];
   armttbar[133]=armttbar[16]*armttbar[177];
   armttbar[163]=armttbar[15]*armttbar[177];
   armttbar[133]=armttbar[133] + 1./2.*armttbar[163];
   armttbar[133]=armttbar[3]*armttbar[133];
   armttbar[160]=armttbar[8]*armttbar[160];
   armttbar[163]=armttbar[12]*armttbar[29];
   armttbar[125]=3./2.*armttbar[133] + 1./8.*armttbar[163] + 1./8.*
   armttbar[125] + armttbar[160];
   armttbar[108]=armttbar[108] - 1 + 33*armttbar[39];
   armttbar[108]=armttbar[156] + armttbar[36] + 1./4.*armttbar[108] - 
   armttbar[38];
   armttbar[108]=3*armttbar[108] + armttbar[155];
   armttbar[108]=armttbar[17]*armttbar[3]*armttbar[108];
   armttbar[133]=armttbar[27]*armttbar[17]*armttbar[3]*armttbar[149];
   armttbar[108]=3*armttbar[133] + 1./2.*armttbar[125] + armttbar[108];
   armttbar[108]=armttbar[27]*armttbar[108];
   armttbar[125]=armttbar[8]*armttbar[169];
   armttbar[133]=1./4.*armttbar[164] + armttbar[125];
   armttbar[149]= - armttbar[3]*armttbar[30];
   armttbar[155]=MMt*armttbar[27]*armttbar[149];
   armttbar[108]=3./4.*armttbar[155] + 1./2.*armttbar[133] + 
   armttbar[108];
   armttbar[108]=MMt*armttbar[108];
   armttbar[108]=1./4.*armttbar[110] + armttbar[108];
   armttbar[108]=armttbar[210]*armttbar[108];
   armttbar[110]= - 55 + 49*armttbar[29];
   armttbar[133]= - 7./2.*armttbar[99] - 3*armttbar[3];
   armttbar[133]=armttbar[17]*armttbar[133];
   armttbar[155]=1./6.*armttbar[10];
   armttbar[110]=3./2.*armttbar[133] + 25*armttbar[190] + 3*
   armttbar[107] + 35./24.*armttbar[12] + 11./6.*armttbar[21] + 
   armttbar[155] + 1./24.*armttbar[110] - armttbar[30];
   armttbar[110]=armttbar[17]*armttbar[110];
   armttbar[133]= - armttbar[33] + armttbar[188];
   armttbar[133]=1./2.*armttbar[133] + armttbar[20];
   armttbar[155]=armttbar[155] + 5./6.*armttbar[29] - armttbar[30];
   armttbar[156]=armttbar[14]*armttbar[155];
   armttbar[133]=3./8.*armttbar[133] + armttbar[156];
   armttbar[156]=1./8.*armttbar[30];
   armttbar[160]=1./48.*armttbar[10] + armttbar[156] - 2./3. - 9./16.*
   armttbar[29];
   armttbar[160]=armttbar[16]*armttbar[160];
   armttbar[163]= - 29 - 49./2.*armttbar[29];
   armttbar[163]=armttbar[15]*armttbar[163];
   armttbar[110]=1./4.*armttbar[110] + 1./96.*armttbar[163] + 1./4.*
   armttbar[133] + armttbar[160];
   armttbar[110]=armttbar[27]*armttbar[110];
   armttbar[129]=armttbar[168] + armttbar[36] + armttbar[38] + 55./4.*
   armttbar[41] + armttbar[129] - 55./4.*armttbar[39] + 57./8. - 
   armttbar[40];
   armttbar[129]=9./2.*armttbar[190] - 15./8.*armttbar[12] + 
   armttbar[207] + 3*armttbar[129] + armttbar[167];
   armttbar[129]=armttbar[17]*armttbar[3]*armttbar[129];
   armttbar[133]=1./3.*armttbar[10];
   armttbar[160]=armttbar[133] + armttbar[235] + 25./3. + armttbar[257]
   ;
   armttbar[160]=armttbar[21]*armttbar[160];
   armttbar[163]=armttbar[23] - 5 - armttbar[43];
   armttbar[163]=3*armttbar[163] - 31./2.*armttbar[29];
   armttbar[160]=1./2.*armttbar[160] + 1./2.*armttbar[163] - 
   armttbar[30];
   armttbar[155]=armttbar[8]*armttbar[155];
   armttbar[163]=1 + armttbar[29];
   armttbar[164]=29./6.*armttbar[163] + armttbar[30];
   armttbar[164]=armttbar[12]*armttbar[164];
   armttbar[155]=1./8.*armttbar[164] + 1./4.*armttbar[160] + 
   armttbar[155];
   armttbar[160]=7 + armttbar[178];
   armttbar[164]=armttbar[200] + 1./4.*armttbar[160] + armttbar[191];
   armttbar[164]=armttbar[16]*armttbar[164];
   armttbar[167]=armttbar[293] - armttbar[19] + armttbar[34] + 
   armttbar[197];
   armttbar[167]=1./2.*armttbar[167] - armttbar[20];
   armttbar[160]=armttbar[160] + 13*armttbar[29];
   armttbar[160]=armttbar[15]*armttbar[160];
   armttbar[160]=1./8.*armttbar[160] + 3*armttbar[167] + armttbar[164];
   armttbar[160]=armttbar[3]*armttbar[160];
   armttbar[164]=armttbar[157] - 5./2.*armttbar[29] + armttbar[235];
   armttbar[164]=armttbar[27]*armttbar[17]*armttbar[3]*armttbar[164];
   armttbar[129]=armttbar[164] + armttbar[129] + 1./2.*armttbar[155] + 
   armttbar[160];
   armttbar[129]=armttbar[27]*armttbar[129];
   armttbar[118]=armttbar[118] + armttbar[10];
   armttbar[155]=armttbar[31]*armttbar[118];
   armttbar[118]=armttbar[1]*armttbar[118];
   armttbar[118]=armttbar[155] + 1./3.*armttbar[118];
   armttbar[155]=armttbar[21]*armttbar[118];
   armttbar[155]=armttbar[238] + armttbar[155];
   armttbar[142]=1./4.*armttbar[155] + armttbar[142];
   armttbar[131]=armttbar[131] + armttbar[244] + 7./2. - armttbar[30];
   armttbar[155]=armttbar[246] + 1./2.*armttbar[248];
   armttbar[155]=armttbar[3]*armttbar[155];
   armttbar[131]=1./2.*armttbar[131] + 3*armttbar[155];
   armttbar[131]=MMt*armttbar[27]*armttbar[3]*armttbar[131];
   armttbar[129]=3*armttbar[131] + 1./2.*armttbar[142] + armttbar[129];
   armttbar[129]=MMt*armttbar[129];
   armttbar[131]=armttbar[14]*armttbar[161];
   armttbar[118]=armttbar[16]*armttbar[118];
   armttbar[142]=armttbar[17]*armttbar[161];
   armttbar[118]=armttbar[142] + armttbar[131] + 1./2.*armttbar[118];
   armttbar[131]=armttbar[147] - 5./6.*armttbar[29] + armttbar[30];
   armttbar[131]=armttbar[27]*armttbar[8]*armttbar[131];
   armttbar[125]=armttbar[125] + armttbar[131];
   armttbar[125]=MMH*armttbar[125];
   armttbar[131]=5*armttbar[271] + 7*armttbar[17];
   armttbar[131]=armttbar[6]*armttbar[27]*armttbar[17]*armttbar[131];
   armttbar[108]=armttbar[108] + armttbar[129] + 1./12.*armttbar[131]
    + 1./8.*armttbar[125] + 1./4.*armttbar[118] + armttbar[110];
   armttbar[108]=armttbar[210]*armttbar[108];
   armttbar[110]= - 23./2.*armttbar[21] + armttbar[229] + armttbar[252]
    + 761./12. - 43*armttbar[29];
   armttbar[110]=1./3.*armttbar[110] - 43./4.*armttbar[12];
   armttbar[118]= - 37*armttbar[16] - 41*armttbar[15];
   armttbar[118]=armttbar[3]*armttbar[118];
   armttbar[125]=11./2.*armttbar[99] - 5*armttbar[3];
   armttbar[125]=armttbar[17]*armttbar[125];
   armttbar[110]=9./2.*armttbar[125] + 1./6.*armttbar[118] + 1./6.*
   armttbar[110] + 15*armttbar[256];
   armttbar[110]=armttbar[17]*armttbar[110];
   armttbar[118]=53./9. + 37./8.*armttbar[29];
   armttbar[118]=armttbar[254] + 1./3.*armttbar[118] + armttbar[152];
   armttbar[118]=armttbar[16]*armttbar[118];
   armttbar[125]=armttbar[253] + 355./9.*armttbar[20];
   armttbar[129]=armttbar[196] + 71 + 101./2.*armttbar[29];
   armttbar[129]=armttbar[15]*armttbar[129];
   armttbar[110]=1./2.*armttbar[110] + 1./48.*armttbar[129] + 1./2.*
   armttbar[118] + 1./16.*armttbar[125] + armttbar[250];
   armttbar[110]=armttbar[27]*armttbar[110];
   armttbar[118]=41*armttbar[16] + 97*armttbar[15];
   armttbar[118]=1./2.*armttbar[118] - 77*armttbar[17];
   armttbar[118]=armttbar[6]*armttbar[27]*armttbar[17]*armttbar[118];
   armttbar[108]=armttbar[108] + armttbar[113] + 1./18.*armttbar[118]
    + armttbar[173] + armttbar[170] + armttbar[110];
   armttbar[108]=armttbar[210]*armttbar[108];
   armttbar[110]=27*armttbar[35];
   armttbar[113]=3*armttbar[28];
   armttbar[118]=armttbar[113] - 365./12. + armttbar[110];
   armttbar[118]=13./36.*armttbar[30] - 7./36.*armttbar[29] + 
   armttbar[122] + 1./2.*armttbar[118] + armttbar[121];
   armttbar[125]= - 2*armttbar[15] + 12*armttbar[14] - 5./2.*
   armttbar[16];
   armttbar[125]=armttbar[3]*armttbar[125];
   armttbar[118]=9*armttbar[132] + armttbar[125] + armttbar[306] + 68./
   9.*armttbar[8] - 97./144.*armttbar[21] - armttbar[36] + 1./2.*
   armttbar[118] - armttbar[38];
   armttbar[118]=armttbar[17]*armttbar[118];
   armttbar[125]=armttbar[19] + 9*armttbar[42] + 3*armttbar[32] - 
   armttbar[34];
   armttbar[129]=3./8.*armttbar[125];
   armttbar[131]= - 3./4.*armttbar[20];
   armttbar[142]=armttbar[156] - 7./72.*armttbar[29] - 20./9. + 3./8.*
   armttbar[28];
   armttbar[142]=armttbar[14]*armttbar[142];
   armttbar[147]=1 - 7./2.*armttbar[29];
   armttbar[147]=1./9.*armttbar[147];
   armttbar[152]=armttbar[147] + 11./2.*armttbar[30];
   armttbar[152]=1./8.*armttbar[16]*armttbar[152];
   armttbar[155]= - armttbar[15]*armttbar[30];
   armttbar[160]=1./18.*armttbar[155];
   armttbar[161]=armttbar[27]*armttbar[255]*armttbar[3];
   armttbar[164]=26./3.*armttbar[161];
   armttbar[167]=armttbar[164] + armttbar[118] + armttbar[160] + 
   armttbar[152] + armttbar[142] + armttbar[131] + armttbar[129] + 400./
   81.*armttbar[18];
   armttbar[167]=armttbar[27]*armttbar[167];
   armttbar[168]= - 4*armttbar[20];
   armttbar[169]= - 3./2.*armttbar[30];
   armttbar[170]= - 1 + armttbar[169];
   armttbar[170]=armttbar[16]*armttbar[170];
   armttbar[171]= - armttbar[19] - 9./2.*armttbar[42] + 3./2.*
   armttbar[32] + armttbar[34];
   armttbar[173]= - 1 - 1./2.*armttbar[28];
   armttbar[173]=armttbar[14]*armttbar[173];
   armttbar[170]=1./2.*armttbar[170] + armttbar[173] + armttbar[171] + 
   armttbar[168];
   armttbar[170]=armttbar[3]*armttbar[170];
   armttbar[178]=armttbar[169] + 7./6.*armttbar[29] + 26./3. + 9*
   armttbar[28];
   armttbar[178]=armttbar[3]*armttbar[178];
   armttbar[183]= - armttbar[17]*armttbar[153];
   armttbar[178]=armttbar[178] + 9*armttbar[183];
   armttbar[178]=armttbar[27]*armttbar[17]*armttbar[178];
   armttbar[185]= - 23./16. - armttbar[43];
   armttbar[115]= - 7./4.*armttbar[28] + 3./8.*armttbar[25] + 
   armttbar[115] + 1./4.*armttbar[185] - 3*armttbar[44];
   armttbar[115]=3*armttbar[115];
   armttbar[147]=armttbar[147] + armttbar[30];
   armttbar[147]=armttbar[21]*armttbar[147];
   armttbar[185]=15*armttbar[28];
   armttbar[187]=armttbar[30] - 7./9.*armttbar[29] - 403./9. + 
   armttbar[185];
   armttbar[187]=armttbar[8]*armttbar[187];
   armttbar[188]= - 5./4.*armttbar[21];
   armttbar[190]=armttbar[188] - 3./2. - armttbar[28];
   armttbar[194]=armttbar[190] - 6*armttbar[8];
   armttbar[194]=armttbar[17]*armttbar[3]*armttbar[194];
   armttbar[147]=armttbar[178] + 3*armttbar[194] + 3*armttbar[170] + 1./
   9.*armttbar[260] + 1./4.*armttbar[187] + 1./8.*armttbar[147] + 
   armttbar[115] - 5./72.*armttbar[30];
   armttbar[147]=armttbar[27]*armttbar[147];
   armttbar[170]= - 3 - armttbar[28];
   armttbar[170]=armttbar[8]*armttbar[170];
   armttbar[170]=armttbar[170] + armttbar[244] + 5./2. + armttbar[28];
   armttbar[170]=armttbar[3]*armttbar[170];
   armttbar[153]= - armttbar[27]*armttbar[17]*armttbar[153]*
   armttbar[28];
   armttbar[178]=armttbar[170] + 12*armttbar[153];
   armttbar[178]=MMt*armttbar[27]*armttbar[178];
   armttbar[147]=armttbar[147] + 3*armttbar[178];
   armttbar[147]=MMt*armttbar[147];
   armttbar[178]= - 3./2.*armttbar[28];
   armttbar[119]=armttbar[119] + 7./72.*armttbar[29] + 101./9. + 
   armttbar[178];
   armttbar[119]=armttbar[8]*armttbar[119];
   armttbar[187]=1./4.*armttbar[28] + 1./2. + armttbar[44];
   armttbar[187]=9*armttbar[187];
   armttbar[119]=armttbar[187] + armttbar[119];
   armttbar[119]=1./2.*MMH*armttbar[27]*armttbar[119];
   armttbar[167]=armttbar[147] + armttbar[167] + armttbar[119];
   armttbar[167]=MMt*armttbar[167];
   armttbar[194]=armttbar[113] - 1091./36. + armttbar[110];
   armttbar[194]=35./36.*armttbar[30] + armttbar[122] + 1./2.*
   armttbar[194] + armttbar[121];
   armttbar[195]=27./8.*armttbar[15] + 6*armttbar[14] + 73./8.*
   armttbar[16];
   armttbar[195]=armttbar[3]*armttbar[195];
   armttbar[132]=9./2.*armttbar[132];
   armttbar[114]=armttbar[132] + armttbar[195] + armttbar[114] + 4*
   armttbar[8] + armttbar[130] + 1./4.*armttbar[194] - armttbar[38];
   armttbar[114]=armttbar[17]*armttbar[114];
   armttbar[130]=3./16.*armttbar[28];
   armttbar[194]=armttbar[156] - 1 + armttbar[130];
   armttbar[194]=armttbar[14]*armttbar[194];
   armttbar[195]=1 + 5*armttbar[30];
   armttbar[195]=armttbar[16]*armttbar[195];
   armttbar[125]=1./2.*armttbar[125] - armttbar[20];
   armttbar[125]=3./8.*armttbar[125];
   armttbar[161]=3*armttbar[161];
   armttbar[114]=armttbar[161] + armttbar[114] + 17./144.*armttbar[155]
    + 1./16.*armttbar[195] + armttbar[125] + armttbar[194];
   armttbar[114]=armttbar[27]*armttbar[114];
   armttbar[194]=1./2. + armttbar[30];
   armttbar[194]=armttbar[21]*armttbar[194];
   armttbar[195]= - 43 + armttbar[185];
   armttbar[195]=1./2.*armttbar[195] + armttbar[30];
   armttbar[195]=armttbar[8]*armttbar[195];
   armttbar[194]=7./72.*armttbar[260] + 1./2.*armttbar[195] + 1./4.*
   armttbar[194] + armttbar[115] - 19./72.*armttbar[30];
   armttbar[195]= - 2*armttbar[20];
   armttbar[171]=1./2.*armttbar[173] + 1./2.*armttbar[171] + 
   armttbar[195];
   armttbar[173]= - 1 - armttbar[30];
   armttbar[173]=armttbar[16]*armttbar[173];
   armttbar[196]=armttbar[15]*armttbar[30];
   armttbar[173]=1./8.*armttbar[196] + armttbar[171] + 1./4.*
   armttbar[173];
   armttbar[173]=armttbar[3]*armttbar[173];
   armttbar[197]=armttbar[145] + 1 + 3./2.*armttbar[28];
   armttbar[197]=armttbar[3]*armttbar[197];
   armttbar[197]=armttbar[197] + 3./2.*armttbar[183];
   armttbar[197]=armttbar[27]*armttbar[17]*armttbar[197];
   armttbar[190]=1./2.*armttbar[190] - 3*armttbar[8];
   armttbar[190]=3*armttbar[17]*armttbar[3]*armttbar[190];
   armttbar[173]=3*armttbar[197] + armttbar[190] + 1./2.*armttbar[194]
    + 3*armttbar[173];
   armttbar[173]=armttbar[27]*armttbar[173];
   armttbar[153]=1./2.*armttbar[170] + 6*armttbar[153];
   armttbar[153]=3*MMt*armttbar[27]*armttbar[153];
   armttbar[170]=armttbar[173] + armttbar[153];
   armttbar[170]=MMt*armttbar[170];
   armttbar[173]=armttbar[148] + 11 + armttbar[178];
   armttbar[173]=armttbar[8]*armttbar[173];
   armttbar[173]=armttbar[187] + armttbar[173];
   armttbar[173]=MMH*armttbar[27]*armttbar[173];
   armttbar[114]=armttbar[170] + armttbar[114] + 1./4.*armttbar[173];
   armttbar[114]=MMt*armttbar[114];
   armttbar[170]=17*armttbar[30];
   armttbar[173]= - 7*armttbar[12] - 1 + armttbar[170];
   armttbar[194]= - armttbar[16] + 5./2.*armttbar[15];
   armttbar[194]=armttbar[3]*armttbar[194];
   armttbar[173]=1./36.*armttbar[173] + armttbar[194];
   armttbar[173]=armttbar[17]*armttbar[173];
   armttbar[173]=17./36.*armttbar[155] + armttbar[173];
   armttbar[173]=armttbar[27]*armttbar[173];
   armttbar[194]= - armttbar[30] + 7*armttbar[260];
   armttbar[197]=armttbar[3]*armttbar[196];
   armttbar[194]=1./18.*armttbar[194] + 3*armttbar[197];
   armttbar[194]=MMt*armttbar[27]*armttbar[194];
   armttbar[173]=armttbar[173] + 1./2.*armttbar[194];
   armttbar[173]=MMt*armttbar[173];
   armttbar[194]= - armttbar[15] + armttbar[17];
   armttbar[194]=armttbar[27]*armttbar[17]*armttbar[194];
   armttbar[173]=17./36.*armttbar[194] + armttbar[173];
   armttbar[173]=armttbar[176]*armttbar[173];
   armttbar[194]=9 + armttbar[273];
   armttbar[194]=armttbar[17]*armttbar[194];
   armttbar[165]=1./2.*armttbar[14] + armttbar[165] + armttbar[20];
   armttbar[165]=9*armttbar[165];
   armttbar[194]=armttbar[165] + armttbar[194];
   armttbar[194]=armttbar[27]*armttbar[194];
   armttbar[197]= - armttbar[8] - 1 - armttbar[44];
   armttbar[197]=9./2.*MMH*armttbar[27]*armttbar[197];
   armttbar[194]=armttbar[194] + armttbar[197];
   armttbar[194]=MMH*armttbar[194];
   armttbar[135]=143./18.*armttbar[17] - 17./18.*armttbar[15] + 
   armttbar[135] + armttbar[16];
   armttbar[135]=armttbar[27]*armttbar[17]*armttbar[135];
   armttbar[135]=armttbar[135] + 1./2.*armttbar[194];
   armttbar[114]=1./4.*armttbar[173] + 1./8.*armttbar[135] + 
   armttbar[114];
   armttbar[114]=armttbar[176]*armttbar[114];
   armttbar[135]= - 107*armttbar[14] + 19./2.*armttbar[16];
   armttbar[173]=1./2.*armttbar[135] - 809./9.*armttbar[15];
   armttbar[194]=19./2.*armttbar[17];
   armttbar[173]=1./3.*armttbar[173] + armttbar[194];
   armttbar[173]=armttbar[27]*armttbar[17]*armttbar[173];
   armttbar[200]=9 - 55./9.*armttbar[8];
   armttbar[200]=armttbar[17]*armttbar[200];
   armttbar[200]=armttbar[165] + armttbar[200];
   armttbar[200]=armttbar[27]*armttbar[200];
   armttbar[200]=armttbar[200] + armttbar[197];
   armttbar[200]=1./4.*MMH*armttbar[200];
   armttbar[173]=1./3.*armttbar[173] + armttbar[200];
   armttbar[114]=armttbar[114] + 1./2.*armttbar[173] + armttbar[167];
   armttbar[114]=armttbar[176]*armttbar[114];
   armttbar[118]=armttbar[164] + armttbar[118] + armttbar[160] + 
   armttbar[152] + armttbar[142] + armttbar[131] + armttbar[129] + 16./
   9.*armttbar[18];
   armttbar[118]=armttbar[27]*armttbar[118];
   armttbar[118]=armttbar[147] + armttbar[118] + armttbar[119];
   armttbar[118]=MMt*armttbar[118];
   armttbar[110]=armttbar[113] - 281./12. + armttbar[110];
   armttbar[110]=armttbar[148] + 19./6.*armttbar[29] + armttbar[122] + 
   1./2.*armttbar[110] + armttbar[121];
   armttbar[113]=2*armttbar[14] + 1./8.*armttbar[16];
   armttbar[113]=3*armttbar[113] - 43./8.*armttbar[15];
   armttbar[113]=armttbar[3]*armttbar[113];
   armttbar[119]= - 1./16.*armttbar[12];
   armttbar[110]=armttbar[132] + armttbar[113] + armttbar[119] + 16./3.
   *armttbar[8] + 1./48.*armttbar[21] + 1./4.*armttbar[110] - 
   armttbar[36];
   armttbar[110]=armttbar[17]*armttbar[110];
   armttbar[113]=19./6.*armttbar[163] + armttbar[235];
   armttbar[113]=armttbar[16]*armttbar[113];
   armttbar[121]=19./24.*armttbar[29] - 1./3. + armttbar[130];
   armttbar[121]=armttbar[14]*armttbar[121];
   armttbar[122]= - armttbar[27]*armttbar[255]*armttbar[3];
   armttbar[110]=5*armttbar[122] + armttbar[110] + 1./16.*armttbar[196]
    + 1./8.*armttbar[113] + armttbar[125] + armttbar[121];
   armttbar[110]=armttbar[27]*armttbar[110];
   armttbar[113]=armttbar[21]*armttbar[163];
   armttbar[121]= - 97./3. + armttbar[185];
   armttbar[121]=1./2.*armttbar[121] + 19./3.*armttbar[29];
   armttbar[121]=armttbar[8]*armttbar[121];
   armttbar[113]=1./8.*armttbar[260] + 1./2.*armttbar[121] + 19./24.*
   armttbar[113] + armttbar[115] + armttbar[156];
   armttbar[115]= - 1./2. - armttbar[30];
   armttbar[115]=armttbar[16]*armttbar[115];
   armttbar[115]=1./8.*armttbar[155] + armttbar[171] + 1./2.*
   armttbar[115];
   armttbar[115]=armttbar[3]*armttbar[115];
   armttbar[121]= - 19./2.*armttbar[29];
   armttbar[125]=armttbar[121] - 5 + 9./2.*armttbar[28];
   armttbar[125]=armttbar[3]*armttbar[125];
   armttbar[125]=armttbar[125] + 9./2.*armttbar[183];
   armttbar[125]=armttbar[27]*armttbar[17]*armttbar[125];
   armttbar[113]=armttbar[125] + armttbar[190] + 1./2.*armttbar[113] + 
   3*armttbar[115];
   armttbar[113]=armttbar[27]*armttbar[113];
   armttbar[113]=armttbar[113] + armttbar[153];
   armttbar[113]=MMt*armttbar[113];
   armttbar[115]= - 19./12.*armttbar[29] + 29./3. + armttbar[178];
   armttbar[115]=armttbar[8]*armttbar[115];
   armttbar[115]=armttbar[187] + armttbar[115];
   armttbar[115]=MMH*armttbar[27]*armttbar[115];
   armttbar[110]=armttbar[113] + armttbar[110] + 1./4.*armttbar[115];
   armttbar[110]=MMt*armttbar[110];
   armttbar[113]=armttbar[169] + 3./2. + armttbar[29];
   armttbar[115]=armttbar[192] + 31./8.*armttbar[15];
   armttbar[115]=armttbar[3]*armttbar[115];
   armttbar[113]=armttbar[115] + armttbar[119] + armttbar[158] - 
   armttbar[36] + 1./8.*armttbar[113] + armttbar[38];
   armttbar[113]=armttbar[17]*armttbar[113];
   armttbar[115]=armttbar[14]*armttbar[177];
   armttbar[119]=armttbar[29] + armttbar[30];
   armttbar[119]=armttbar[16]*armttbar[119];
   armttbar[115]=1./2.*armttbar[196] + armttbar[115] + 1./2.*
   armttbar[119];
   armttbar[113]=1./8.*armttbar[115] + armttbar[113];
   armttbar[113]=armttbar[27]*armttbar[113];
   armttbar[115]=armttbar[21]*armttbar[123];
   armttbar[115]=armttbar[175] + armttbar[115];
   armttbar[119]= - armttbar[16]*armttbar[30];
   armttbar[123]=armttbar[119] + 1./2.*armttbar[155];
   armttbar[123]=armttbar[3]*armttbar[123];
   armttbar[125]=armttbar[8]*armttbar[177];
   armttbar[115]=3*armttbar[123] + 1./4.*armttbar[260] + 1./2.*
   armttbar[115] + armttbar[125];
   armttbar[123]=armttbar[27]*armttbar[17]*armttbar[3]*armttbar[140];
   armttbar[115]=1./2.*armttbar[115] + 3*armttbar[123];
   armttbar[115]=MMt*armttbar[27]*armttbar[115];
   armttbar[123]=MMH*armttbar[27]*armttbar[8]*armttbar[140];
   armttbar[113]=1./2.*armttbar[115] + armttbar[113] + 1./16.*
   armttbar[123];
   armttbar[113]=MMt*armttbar[113];
   armttbar[115]=armttbar[189] + armttbar[120];
   armttbar[115]=armttbar[27]*armttbar[17]*armttbar[115];
   armttbar[113]=1./8.*armttbar[115] + armttbar[113];
   armttbar[113]=armttbar[210]*armttbar[113];
   armttbar[115]=9 - 37./3.*armttbar[8];
   armttbar[115]=armttbar[17]*armttbar[115];
   armttbar[115]=armttbar[165] + armttbar[115];
   armttbar[115]=armttbar[27]*armttbar[115];
   armttbar[115]=armttbar[115] + armttbar[197];
   armttbar[115]=MMH*armttbar[115];
   armttbar[123]= - 17./2.*armttbar[14] + 7*armttbar[16];
   armttbar[123]=71./12.*armttbar[17] + 1./3.*armttbar[123] + 1./4.*
   armttbar[15];
   armttbar[123]=armttbar[27]*armttbar[17]*armttbar[123];
   armttbar[115]=armttbar[123] + 1./4.*armttbar[115];
   armttbar[110]=armttbar[113] + 1./4.*armttbar[115] + armttbar[110];
   armttbar[110]=armttbar[210]*armttbar[110];
   armttbar[113]=armttbar[194] + 1./6.*armttbar[135] - 11*armttbar[15];
   armttbar[113]=armttbar[27]*armttbar[17]*armttbar[113];
   armttbar[113]=1./3.*armttbar[113] + armttbar[200];
   armttbar[110]=armttbar[110] + 1./2.*armttbar[113] + armttbar[118];
   armttbar[110]=armttbar[210]*armttbar[110];
   armttbar[113]=armttbar[14]*armttbar[30];
   armttbar[115]=1./2.*armttbar[16]*armttbar[30];
   armttbar[118]= - armttbar[20] - armttbar[34] + armttbar[19];
   armttbar[123]=armttbar[115] + 3*armttbar[118] + armttbar[113];
   armttbar[125]=armttbar[8] + armttbar[188] + 3 + armttbar[228];
   armttbar[125]=armttbar[17]*armttbar[125];
   armttbar[123]=1./2.*armttbar[123] + armttbar[125];
   armttbar[122]=3*armttbar[122];
   armttbar[123]=1./2.*armttbar[123] + armttbar[122];
   armttbar[123]=armttbar[27]*armttbar[123];
   armttbar[125]=armttbar[21]*armttbar[30];
   armttbar[129]= - armttbar[43] + armttbar[23];
   armttbar[129]=1./2.*armttbar[129] - armttbar[25];
   armttbar[130]=3*armttbar[129];
   armttbar[131]=armttbar[8]*armttbar[30];
   armttbar[132]=armttbar[131] + armttbar[130] + 1./4.*armttbar[125];
   armttbar[135]=3*armttbar[27]*armttbar[17]*armttbar[149];
   armttbar[132]=1./2.*armttbar[132] + armttbar[135];
   armttbar[132]=MMt*armttbar[27]*armttbar[132];
   armttbar[140]= - armttbar[8]*armttbar[30];
   armttbar[142]=1./8.*MMH*armttbar[27]*armttbar[140];
   armttbar[123]=armttbar[132] + armttbar[123] + armttbar[142];
   armttbar[123]=MMt*armttbar[123];
   armttbar[132]=armttbar[14] - 3./2.*armttbar[16];
   armttbar[132]=3*armttbar[132] + armttbar[17];
   armttbar[132]=armttbar[27]*armttbar[17]*armttbar[132];
   armttbar[147]=armttbar[38] - 1./4.*armttbar[8];
   armttbar[147]=MMH*armttbar[27]*armttbar[17]*armttbar[147];
   armttbar[132]=1./2.*armttbar[132] + armttbar[147];
   armttbar[123]=1./2.*armttbar[132] + armttbar[123];
   armttbar[123]=MMt*armttbar[123];
   armttbar[132]=3./2.*armttbar[118];
   armttbar[113]=armttbar[115] + armttbar[132] + armttbar[113];
   armttbar[115]= - armttbar[21] + 3 + armttbar[30];
   armttbar[115]=1./2.*armttbar[115] + armttbar[8];
   armttbar[115]=armttbar[17]*armttbar[115];
   armttbar[113]=1./2.*armttbar[113] + armttbar[115];
   armttbar[113]=1./2.*armttbar[113] + armttbar[122];
   armttbar[113]=armttbar[27]*armttbar[113];
   armttbar[115]=armttbar[130] + 1./2.*armttbar[125];
   armttbar[115]=1./2.*armttbar[115] + armttbar[131];
   armttbar[115]=1./2.*armttbar[115] + armttbar[135];
   armttbar[115]=MMt*armttbar[27]*armttbar[115];
   armttbar[113]=armttbar[115] + armttbar[113] + armttbar[142];
   armttbar[113]=MMt*armttbar[113];
   armttbar[115]=3*armttbar[136] + armttbar[17];
   armttbar[115]=armttbar[27]*armttbar[17]*armttbar[115];
   armttbar[115]=1./2.*armttbar[115] + armttbar[147];
   armttbar[113]=1./2.*armttbar[115] + armttbar[113];
   armttbar[113]=armttbar[176]*MMt*armttbar[113];
   armttbar[113]=armttbar[123] + armttbar[113];
   armttbar[113]=armttbar[176]*armttbar[113];
   armttbar[115]=1 + armttbar[244];
   armttbar[115]=armttbar[17]*armttbar[115];
   armttbar[115]=1./2.*armttbar[118] + armttbar[115];
   armttbar[115]=armttbar[27]*armttbar[115];
   armttbar[118]=MMt*armttbar[27]*armttbar[129];
   armttbar[115]=armttbar[115] + armttbar[118];
   armttbar[115]=MMt*armttbar[115];
   armttbar[118]= - armttbar[15] + armttbar[14] - 9./4.*armttbar[16];
   armttbar[118]=armttbar[27]*armttbar[17]*armttbar[118];
   armttbar[122]=MMH*armttbar[27]*armttbar[17]*armttbar[36];
   armttbar[115]=9./2.*armttbar[115] + armttbar[118] + armttbar[122];
   armttbar[115]=1./2.*MMt*armttbar[115];
   armttbar[113]=armttbar[115] + armttbar[113];
   armttbar[113]=armttbar[176]*armttbar[113];
   armttbar[118]= - armttbar[14] + 3./2.*armttbar[16];
   armttbar[118]=armttbar[120] + 1./2.*armttbar[118] - armttbar[15];
   armttbar[118]=armttbar[27]*armttbar[17]*armttbar[118];
   armttbar[122]=armttbar[288] + armttbar[159];
   armttbar[122]=MMH*armttbar[27]*armttbar[17]*armttbar[122];
   armttbar[118]=armttbar[118] + armttbar[122];
   armttbar[123]= - armttbar[14]*armttbar[30];
   armttbar[119]=1./2.*armttbar[119];
   armttbar[125]=armttbar[123] + armttbar[119];
   armttbar[129]= - armttbar[30] + armttbar[244];
   armttbar[129]=1./2.*armttbar[129] - armttbar[8];
   armttbar[129]=armttbar[17]*armttbar[129];
   armttbar[125]=1./2.*armttbar[125] + armttbar[129];
   armttbar[125]=1./2.*armttbar[125] + armttbar[161];
   armttbar[125]=armttbar[27]*armttbar[125];
   armttbar[129]= - armttbar[21]*armttbar[30];
   armttbar[135]=1./4.*armttbar[129] + armttbar[140];
   armttbar[136]=3*armttbar[27]*armttbar[17]*armttbar[261];
   armttbar[135]=1./2.*armttbar[135] + armttbar[136];
   armttbar[135]=MMt*armttbar[27]*armttbar[135];
   armttbar[131]=1./8.*MMH*armttbar[27]*armttbar[131];
   armttbar[125]=armttbar[135] + armttbar[125] + armttbar[131];
   armttbar[125]=MMt*armttbar[125];
   armttbar[118]=1./2.*armttbar[118] + armttbar[125];
   armttbar[118]=armttbar[210]*MMt*armttbar[118];
   armttbar[119]=armttbar[119] + armttbar[132] + armttbar[123];
   armttbar[123]=3 - armttbar[30];
   armttbar[123]= - armttbar[8] + 1./2.*armttbar[123] - armttbar[21];
   armttbar[123]=armttbar[17]*armttbar[123];
   armttbar[119]=1./2.*armttbar[119] + armttbar[123];
   armttbar[119]=1./2.*armttbar[119] + armttbar[161];
   armttbar[119]=armttbar[27]*armttbar[119];
   armttbar[123]=armttbar[130] + 1./2.*armttbar[129];
   armttbar[123]=1./2.*armttbar[123] + armttbar[140];
   armttbar[123]=1./2.*armttbar[123] + armttbar[136];
   armttbar[123]=MMt*armttbar[27]*armttbar[123];
   armttbar[119]=armttbar[123] + armttbar[119] + armttbar[131];
   armttbar[119]=MMt*armttbar[119];
   armttbar[120]=armttbar[120] - 1./2.*armttbar[14] - armttbar[15];
   armttbar[120]=armttbar[27]*armttbar[17]*armttbar[120];
   armttbar[120]=armttbar[120] + armttbar[122];
   armttbar[119]=1./2.*armttbar[120] + armttbar[119];
   armttbar[119]=MMt*armttbar[119];
   armttbar[118]=armttbar[119] + armttbar[118];
   armttbar[118]=armttbar[210]*armttbar[118];
   armttbar[115]=armttbar[115] + armttbar[118];
   armttbar[115]=armttbar[210]*armttbar[115];
   armttbar[113]=armttbar[113] + armttbar[115];
   armttbar[113]=armttbar[9]*armttbar[113];
   armttbar[115]=armttbar[27]*armttbar[17]*armttbar[15];
   armttbar[118]= - MMt*armttbar[27]*armttbar[18];
   armttbar[115]=armttbar[115] + armttbar[118];
   armttbar[110]=1./2.*armttbar[113] + armttbar[110] + 1024./81.*
   armttbar[115] + armttbar[114];
   armttbar[110]=armttbar[9]*armttbar[110];
   armttbar[113]= - armttbar[29] - armttbar[12];
   armttbar[107]=16*armttbar[216] + 1./3.*armttbar[113] + 10*
   armttbar[107];
   armttbar[107]=armttbar[17]*armttbar[107];
   armttbar[113]= - armttbar[18] - 8./3.*armttbar[20];
   armttbar[114]= - 1 - 1./9.*armttbar[29];
   armttbar[114]=armttbar[15]*armttbar[114];
   armttbar[107]=2./3.*armttbar[107] + 2./3.*armttbar[113] + 
   armttbar[114];
   armttbar[107]=armttbar[27]*armttbar[107];
   armttbar[113]= - 5*armttbar[15] + 11*armttbar[17];
   armttbar[113]=armttbar[6]*armttbar[27]*armttbar[17]*armttbar[113];
   armttbar[114]= - 17 - armttbar[29];
   armttbar[114]=armttbar[12]*armttbar[114];
   armttbar[114]=armttbar[114] + armttbar[29] - 8*armttbar[25] + 1 - 8*
   armttbar[45];
   armttbar[114]=MMt*armttbar[27]*armttbar[114];
   armttbar[107]=2./9.*armttbar[114] + armttbar[107] + 2./9.*
   armttbar[113];
   armttbar[107]=armttbar[110] + armttbar[108] + 64./9.*armttbar[107]
    + armttbar[150];
   armttbar[107]=armttbar[9]*armttbar[107];
   armttbar[108]=2059./3. - 1253*armttbar[29];
   armttbar[108]= - 35./6.*armttbar[10] + 1./24.*armttbar[108] + 
   armttbar[259];
   armttbar[108]=armttbar[12]*armttbar[108];
   armttbar[110]=491./9. + 301./8.*armttbar[45];
   armttbar[113]= - 5*armttbar[5];
   armttbar[114]=41./4.*armttbar[33] + armttbar[113];
   armttbar[114]= - 61./6.*armttbar[20] + 1./3.*armttbar[114] - 7./4.*
   armttbar[18];
   armttbar[114]=armttbar[99]*armttbar[114];
   armttbar[112]=5./3.*armttbar[10] - 29./3. + armttbar[112];
   armttbar[112]=armttbar[15]*armttbar[99]*armttbar[112];
   armttbar[108]=25./12.*armttbar[112] + 25./12.*armttbar[114] + 1./12.
   *armttbar[108] + 25./216.*armttbar[10] - 185./864.*armttbar[29] - 35.
   /72.*armttbar[13] + 35./72.*armttbar[24] + 1./9.*armttbar[110] + 5./
   8.*armttbar[25];
   armttbar[110]= - 17*armttbar[33] + armttbar[113];
   armttbar[112]=5*armttbar[10] + 11 + 17*armttbar[29];
   armttbar[112]=armttbar[15]*armttbar[112];
   armttbar[110]=1./2.*armttbar[112] + 17*armttbar[20] + 1./2.*
   armttbar[110] + 11*armttbar[18];
   armttbar[110]=armttbar[3]*armttbar[110];
   armttbar[112]= - 7 - 25*armttbar[29];
   armttbar[112]=25./4.*armttbar[12] + 1./4.*armttbar[112] + 
   armttbar[151];
   armttbar[112]=armttbar[99]*armttbar[112];
   armttbar[113]=101./3. + 85./2.*armttbar[12];
   armttbar[113]=armttbar[3]*armttbar[113];
   armttbar[112]=25./6.*armttbar[112] + armttbar[113];
   armttbar[112]=armttbar[17]*armttbar[112];
   armttbar[113]=5./4.*armttbar[106] + armttbar[133];
   armttbar[113]=armttbar[12]*armttbar[113];
   armttbar[113]=5./3.*armttbar[113] + armttbar[230] - 25./12.*
   armttbar[29] + 5./9.*armttbar[13] - 5./9.*armttbar[24] + 11./2. - 5./
   9.*armttbar[25];
   armttbar[113]=armttbar[99]*armttbar[113];
   armttbar[114]=35./4.*armttbar[10] - 121./3. + 119./4.*armttbar[29];
   armttbar[114]=armttbar[3]*armttbar[114];
   armttbar[113]=25./16.*armttbar[113] + 1./3.*armttbar[114];
   armttbar[113]=MMZ*armttbar[113];
   armttbar[108]=armttbar[113] + 1./2.*armttbar[112] + 1./2.*
   armttbar[108] + armttbar[110];
   armttbar[108]=armttbar[27]*armttbar[108];
   armttbar[110]=5./27.*armttbar[10] - 7./9.*armttbar[13] + 7./9.*
   armttbar[24] - 23./324. + armttbar[25];
   armttbar[110]=armttbar[31]*armttbar[110];
   armttbar[112]=5./9.*armttbar[10] - 7./3.*armttbar[13] + 7./3.*
   armttbar[24] - 23./108. + 3*armttbar[25];
   armttbar[112]=armttbar[1]*armttbar[112];
   armttbar[110]=11./3.*armttbar[110] + armttbar[112];
   armttbar[112]= - 1./3. + armttbar[157];
   armttbar[113]=armttbar[31]*armttbar[112];
   armttbar[112]=armttbar[1]*armttbar[112];
   armttbar[112]=11./9.*armttbar[113] + armttbar[112];
   armttbar[112]=armttbar[12]*armttbar[112];
   armttbar[110]=1./2.*armttbar[110] + 7./3.*armttbar[112];
   armttbar[112]=11./9.*armttbar[137] + armttbar[134];
   armttbar[112]=armttbar[99]*armttbar[112];
   armttbar[113]=11./9.*armttbar[139] + armttbar[141];
   armttbar[114]=armttbar[15]*armttbar[99]*armttbar[113];
   armttbar[110]=25./6.*armttbar[114] + 1./2.*armttbar[110] + 25./3.*
   armttbar[112];
   armttbar[112]=11./3.*armttbar[174] + 3*armttbar[172];
   armttbar[112]=armttbar[15]*armttbar[112];
   armttbar[114]=armttbar[31]*armttbar[126];
   armttbar[115]=armttbar[1]*armttbar[126];
   armttbar[112]=armttbar[112] + 11./3.*armttbar[114] + 3*armttbar[115]
   ;
   armttbar[112]=armttbar[3]*armttbar[112];
   armttbar[113]=armttbar[12]*armttbar[113];
   armttbar[113]=armttbar[113] + 11./9.*armttbar[166] + armttbar[162];
   armttbar[113]=armttbar[99]*armttbar[113];
   armttbar[104]=armttbar[31]*armttbar[104];
   armttbar[104]=11./9.*armttbar[104] + armttbar[179];
   armttbar[104]=armttbar[3]*armttbar[104];
   armttbar[104]=25./24.*armttbar[113] + armttbar[104];
   armttbar[104]=MMZ*armttbar[104];
   armttbar[113]=11./9.*armttbar[182] + armttbar[181];
   armttbar[113]=armttbar[17]*armttbar[99]*armttbar[113];
   armttbar[104]=1./3.*armttbar[108] + armttbar[104] + 25./6.*
   armttbar[113] + 1./2.*armttbar[110] + armttbar[112];
   armttbar[108]=armttbar[31]*armttbar[138];
   armttbar[110]=11./9.*armttbar[108] + armttbar[144];
   armttbar[112]=armttbar[12]*armttbar[110];
   armttbar[113]=armttbar[24] - armttbar[13];
   armttbar[114]=armttbar[31]*armttbar[113];
   armttbar[113]=armttbar[1]*armttbar[113];
   armttbar[112]=armttbar[112] + 11./9.*armttbar[114] + armttbar[113];
   armttbar[112]=MMZ*armttbar[112];
   armttbar[110]=armttbar[15]*armttbar[110];
   armttbar[115]= - 7./6. + armttbar[10];
   armttbar[118]=armttbar[31]*armttbar[115];
   armttbar[115]=armttbar[1]*armttbar[115];
   armttbar[115]=11./9.*armttbar[118] + armttbar[115];
   armttbar[115]=armttbar[17]*armttbar[115];
   armttbar[118]=armttbar[209] + armttbar[20];
   armttbar[119]=armttbar[31]*armttbar[118];
   armttbar[118]=armttbar[1]*armttbar[118];
   armttbar[110]=armttbar[112] + armttbar[115] + armttbar[110] + 11./9.
   *armttbar[119] + armttbar[118];
   armttbar[112]=armttbar[251] + 1./4.*armttbar[17];
   armttbar[112]=armttbar[17]*armttbar[112];
   armttbar[112]=armttbar[112] + 1./6.*armttbar[146];
   armttbar[112]=armttbar[6]*armttbar[27]*armttbar[112];
   armttbar[115]=armttbar[33] + 5*armttbar[5];
   armttbar[118]=armttbar[229] + 13./3. - 17./2.*armttbar[29];
   armttbar[118]=armttbar[15]*armttbar[118];
   armttbar[115]=1./6.*armttbar[118] + 1./4.*armttbar[20] + 1./12.*
   armttbar[115] + armttbar[18];
   armttbar[118]= - 13./8.*armttbar[12] + 5./16.*armttbar[10] - 2./3.
    + 17./16.*armttbar[29];
   armttbar[118]=armttbar[17]*armttbar[118];
   armttbar[115]=1./4.*armttbar[115] + 1./3.*armttbar[118];
   armttbar[118]= - 5*armttbar[10] + 83./3. + armttbar[223];
   armttbar[118]=armttbar[12]*armttbar[118];
   armttbar[118]=1./18.*armttbar[118] - 5./18.*armttbar[13] + 5./18.*
   armttbar[24] + 1./18. + armttbar[45];
   armttbar[118]=MMZ*armttbar[118];
   armttbar[115]=1./3.*armttbar[115] + 1./8.*armttbar[118];
   armttbar[115]=armttbar[27]*armttbar[115];
   armttbar[110]=17./24.*armttbar[112] + 1./8.*armttbar[110] + 
   armttbar[115];
   armttbar[110]=armttbar[6]*armttbar[110];
   armttbar[112]=17./3.*armttbar[12] + armttbar[30] + 49./18.*
   armttbar[29] - 71./18. + armttbar[227];
   armttbar[112]=armttbar[3]*armttbar[112];
   armttbar[112]=1./2.*armttbar[112] + armttbar[231];
   armttbar[112]=armttbar[27]*armttbar[112];
   armttbar[112]=1./2.*armttbar[112] + 3*armttbar[180];
   armttbar[112]=MMt*armttbar[112];
   armttbar[104]=armttbar[112] + 1./2.*armttbar[104] + 17./3.*
   armttbar[110];
   armttbar[104]=armttbar[176]*armttbar[104];
   armttbar[110]=11./9.*armttbar[215] + armttbar[214];
   armttbar[110]=armttbar[99]*armttbar[110];
   armttbar[112]=517./36.*armttbar[10] + 625./27. + armttbar[206];
   armttbar[112]=armttbar[31]*armttbar[112];
   armttbar[115]=173./3. + armttbar[206];
   armttbar[115]=1./3.*armttbar[115] + 47./4.*armttbar[10];
   armttbar[115]=armttbar[1]*armttbar[115];
   armttbar[112]=armttbar[112] + armttbar[115];
   armttbar[112]=armttbar[12]*armttbar[112];
   armttbar[115]=11./9.*armttbar[221] + armttbar[222];
   armttbar[118]=armttbar[15]*armttbar[99]*armttbar[115];
   armttbar[119]=11*armttbar[13] - 11*armttbar[24] + 131./6. + 11*
   armttbar[25];
   armttbar[119]= - 143./24.*armttbar[10] + 5./9.*armttbar[119] + 1./8.
   *armttbar[11];
   armttbar[119]=armttbar[31]*armttbar[119];
   armttbar[120]=armttbar[208] + 5*armttbar[13] - 5*armttbar[24] + 179./
   18. + 5*armttbar[25];
   armttbar[120]=1./3.*armttbar[120] - 13./8.*armttbar[10];
   armttbar[120]=armttbar[1]*armttbar[120];
   armttbar[110]=10*armttbar[118] + 10*armttbar[110] + 1./6.*
   armttbar[112] + 1./3.*armttbar[119] + armttbar[120];
   armttbar[112]=67./12.*armttbar[29] + 377./27.*armttbar[13] - 377./27.
   *armttbar[24] + 161./27.*armttbar[25] - 11951./162. + armttbar[117];
   armttbar[112]=1./2.*armttbar[112] + armttbar[186];
   armttbar[117]=54977./3. + 16313*armttbar[29];
   armttbar[117]=1./108.*armttbar[117] - armttbar[30];
   armttbar[117]=1./4.*armttbar[117] + 43./27.*armttbar[10];
   armttbar[117]=armttbar[12]*armttbar[117];
   armttbar[118]= - 523./4.*armttbar[33] + 25*armttbar[5];
   armttbar[118]=623./18.*armttbar[20] + 1./9.*armttbar[118] + 47./4.*
   armttbar[18];
   armttbar[118]=armttbar[99]*armttbar[118];
   armttbar[119]= - 25./3.*armttbar[10] + 337./3. - 275./4.*
   armttbar[29];
   armttbar[119]=armttbar[15]*armttbar[99]*armttbar[119];
   armttbar[112]=5./6.*armttbar[119] + 5./2.*armttbar[118] + 
   armttbar[117] + 1./2.*armttbar[112] + armttbar[151];
   armttbar[116]=11./4.*armttbar[116] + armttbar[198];
   armttbar[116]=armttbar[12]*armttbar[116];
   armttbar[117]=25./9.*armttbar[10];
   armttbar[116]=25./3.*armttbar[116] + armttbar[117] + 275./12.*
   armttbar[29] - 25./9.*armttbar[13] + 25./9.*armttbar[24] - 133./2.
    + 25./9.*armttbar[25];
   armttbar[116]=armttbar[99]*armttbar[116];
   armttbar[118]=10./3. + armttbar[121];
   armttbar[118]=4./9.*armttbar[10] + 1./9.*armttbar[118] + 
   armttbar[145];
   armttbar[118]=armttbar[3]*armttbar[118];
   armttbar[116]=5./48.*armttbar[116] + armttbar[118];
   armttbar[116]=MMZ*armttbar[116];
   armttbar[118]=armttbar[202] + armttbar[199];
   armttbar[118]=armttbar[3]*armttbar[118];
   armttbar[119]=47 + 275./3.*armttbar[29];
   armttbar[117]= - 275./12.*armttbar[12] + 1./4.*armttbar[119] + 
   armttbar[117];
   armttbar[117]=armttbar[99]*armttbar[117];
   armttbar[119]= - 47./9. + armttbar[211];
   armttbar[119]=armttbar[3]*armttbar[119];
   armttbar[117]=5./6.*armttbar[117] + armttbar[119];
   armttbar[117]=armttbar[17]*armttbar[117];
   armttbar[112]=armttbar[116] + 1./2.*armttbar[117] + 1./12.*
   armttbar[112] + armttbar[118];
   armttbar[112]=armttbar[27]*armttbar[112];
   armttbar[116]= - armttbar[24] + armttbar[13];
   armttbar[117]=armttbar[31]*armttbar[116];
   armttbar[116]=armttbar[1]*armttbar[116];
   armttbar[116]=11./9.*armttbar[117] + armttbar[116];
   armttbar[117]=1067./36.*armttbar[10] + 245./27. + armttbar[154];
   armttbar[117]=armttbar[31]*armttbar[117];
   armttbar[118]=73./3. + armttbar[154];
   armttbar[118]=1./3.*armttbar[118] + 97./4.*armttbar[10];
   armttbar[118]=armttbar[1]*armttbar[118];
   armttbar[117]=armttbar[117] + armttbar[118];
   armttbar[118]=armttbar[12]*armttbar[117];
   armttbar[116]=10*armttbar[116] + 1./2.*armttbar[118];
   armttbar[116]=MMZ*armttbar[116];
   armttbar[105]=11./9.*armttbar[109] + armttbar[105];
   armttbar[109]=armttbar[15]*armttbar[117];
   armttbar[117]= - 1067./36.*armttbar[10] + 745./27. + armttbar[184];
   armttbar[117]=armttbar[31]*armttbar[117];
   armttbar[118]=197./3. + armttbar[184];
   armttbar[118]=1./3.*armttbar[118] - 97./4.*armttbar[10];
   armttbar[118]=armttbar[1]*armttbar[118];
   armttbar[117]=armttbar[117] + armttbar[118];
   armttbar[117]=armttbar[17]*armttbar[117];
   armttbar[105]=armttbar[116] + 1./2.*armttbar[117] + 10*armttbar[105]
    + 1./2.*armttbar[109];
   armttbar[109]= - 601./3. - 391./2.*armttbar[29];
   armttbar[109]= - 179./18.*armttbar[10] + 5./9.*armttbar[109] + 
   armttbar[259];
   armttbar[109]=1./2.*armttbar[109] + 185*armttbar[12];
   armttbar[109]=armttbar[17]*armttbar[109];
   armttbar[116]=611*armttbar[33] - 107*armttbar[5];
   armttbar[116]= - 443./6.*armttbar[20] + 1./18.*armttbar[116] - 37*
   armttbar[18];
   armttbar[117]= - 793./3. + 1955./8.*armttbar[29];
   armttbar[117]=179./72.*armttbar[10] + 1./9.*armttbar[117] + 17./4.*
   armttbar[30];
   armttbar[117]=armttbar[15]*armttbar[117];
   armttbar[109]=1./2.*armttbar[109] + 1./2.*armttbar[116] + 
   armttbar[117];
   armttbar[116]= - 365./3. + 391./2.*armttbar[29];
   armttbar[116]=179./18.*armttbar[10] + 5./9.*armttbar[116] + 
   armttbar[170];
   armttbar[116]=armttbar[12]*armttbar[116];
   armttbar[116]=1./18.*armttbar[116] + 107./162.*armttbar[13] - 107./
   162.*armttbar[24] + 611./162. - armttbar[45];
   armttbar[116]=MMZ*armttbar[116];
   armttbar[109]=1./9.*armttbar[109] + 1./2.*armttbar[116];
   armttbar[109]=armttbar[27]*armttbar[109];
   armttbar[105]=773./108.*armttbar[111] + 1./9.*armttbar[105] + 1./2.*
   armttbar[109];
   armttbar[105]=armttbar[6]*armttbar[105];
   armttbar[109]=armttbar[12]*armttbar[115];
   armttbar[109]=armttbar[109] + 11./9.*armttbar[226] + armttbar[225];
   armttbar[109]=armttbar[99]*armttbar[109];
   armttbar[111]=armttbar[3]*armttbar[242];
   armttbar[109]=5./3.*armttbar[109] + armttbar[111];
   armttbar[109]=MMZ*armttbar[109];
   armttbar[111]= - armttbar[12] + armttbar[228] - 119./9.*armttbar[29]
    + 113./18. + armttbar[227];
   armttbar[111]=armttbar[3]*armttbar[111];
   armttbar[111]=1./2.*armttbar[111] + armttbar[231];
   armttbar[111]=armttbar[27]*armttbar[111];
   armttbar[111]=armttbar[111] + armttbar[232];
   armttbar[111]=MMt*armttbar[111];
   armttbar[115]=armttbar[3]*armttbar[218];
   armttbar[116]=11./9.*armttbar[234] + armttbar[233];
   armttbar[116]=armttbar[17]*armttbar[99]*armttbar[116];
   armttbar[104]=armttbar[104] + armttbar[111] + armttbar[105] + 
   armttbar[112] + armttbar[109] + 20./3.*armttbar[116] + 1./3.*
   armttbar[110] + armttbar[115];
   armttbar[104]=armttbar[176]*armttbar[104];
   armttbar[105]=61./9. - 4*armttbar[45];
   armttbar[109]= - 2./3.*armttbar[29];
   armttbar[105]=armttbar[133] + armttbar[109] + 1./3.*armttbar[13] - 1.
   /3.*armttbar[24] + 2*armttbar[105] - 23./3.*armttbar[25];
   armttbar[110]= - 5./3. + armttbar[29];
   armttbar[110]=2*armttbar[110] + armttbar[133];
   armttbar[110]=armttbar[15]*armttbar[99]*armttbar[110];
   armttbar[111]= - 1 - armttbar[29];
   armttbar[111]=armttbar[21]*armttbar[111];
   armttbar[112]= - 97./3. + armttbar[258];
   armttbar[112]=2*armttbar[112] - armttbar[10];
   armttbar[112]=armttbar[12]*armttbar[112];
   armttbar[115]=4*armttbar[33] - armttbar[5];
   armttbar[115]= - 10./3.*armttbar[20] + 1./3.*armttbar[115] - 
   armttbar[18];
   armttbar[115]=armttbar[99]*armttbar[115];
   armttbar[116]=2*armttbar[12] + armttbar[198] - 1 - 2*armttbar[29];
   armttbar[116]=armttbar[99]*armttbar[116];
   armttbar[116]=2*armttbar[116] + armttbar[3];
   armttbar[116]=armttbar[17]*armttbar[116];
   armttbar[105]=4*armttbar[116] + 4*armttbar[110] + 4*armttbar[115] + 
   2./9.*armttbar[112] + 2./3.*armttbar[105] + armttbar[111];
   armttbar[106]=2*armttbar[106] + armttbar[133];
   armttbar[106]=armttbar[12]*armttbar[106];
   armttbar[106]=1./3.*armttbar[106] - 1./9.*armttbar[10] + 
   armttbar[109] + 1./9.*armttbar[13] - 1./9.*armttbar[24] + 2 - 1./9.*
   armttbar[25];
   armttbar[106]=armttbar[99]*armttbar[106];
   armttbar[109]=armttbar[3]*armttbar[193];
   armttbar[106]=4*armttbar[106] + 1./3.*armttbar[109];
   armttbar[106]=MMZ*armttbar[106];
   armttbar[105]=2./3.*armttbar[105] + armttbar[106];
   armttbar[105]=armttbar[27]*armttbar[105];
   armttbar[106]= - 4./3. - armttbar[10];
   armttbar[109]=armttbar[31]*armttbar[106];
   armttbar[106]=armttbar[1]*armttbar[106];
   armttbar[106]=5./3.*armttbar[109] + armttbar[106];
   armttbar[106]=armttbar[12]*armttbar[106];
   armttbar[109]=armttbar[10] + armttbar[13] - armttbar[24] + 16./3. + 
   armttbar[25];
   armttbar[110]=armttbar[31]*armttbar[109];
   armttbar[109]=armttbar[1]*armttbar[109];
   armttbar[106]=armttbar[106] + 5./3.*armttbar[110] + armttbar[109];
   armttbar[109]=armttbar[126] + armttbar[195];
   armttbar[110]=armttbar[31]*armttbar[109];
   armttbar[109]=armttbar[1]*armttbar[109];
   armttbar[109]=5./3.*armttbar[110] + armttbar[109];
   armttbar[109]=armttbar[99]*armttbar[109];
   armttbar[110]=5./3.*armttbar[139] + armttbar[141];
   armttbar[111]=armttbar[15]*armttbar[99]*armttbar[110];
   armttbar[112]=5./3.*armttbar[182] + armttbar[181];
   armttbar[112]=armttbar[17]*armttbar[99]*armttbar[112];
   armttbar[106]=4*armttbar[112] + 2*armttbar[111] + 1./3.*
   armttbar[106] + 2*armttbar[109];
   armttbar[109]=armttbar[12]*armttbar[110];
   armttbar[109]=armttbar[109] + 5./3.*armttbar[166] + armttbar[162];
   armttbar[109]=armttbar[99]*armttbar[109];
   armttbar[110]=armttbar[3]*armttbar[143];
   armttbar[109]=4./3.*armttbar[109] + armttbar[110];
   armttbar[109]=MMZ*armttbar[109];
   armttbar[105]=armttbar[105] + 4./3.*armttbar[106] + armttbar[109];
   armttbar[106]=5./3.*armttbar[108] + armttbar[144];
   armttbar[108]=armttbar[15]*armttbar[106];
   armttbar[109]= - 11./3. + armttbar[10];
   armttbar[110]=armttbar[31]*armttbar[109];
   armttbar[109]=armttbar[1]*armttbar[109];
   armttbar[109]=5./3.*armttbar[110] + armttbar[109];
   armttbar[109]=armttbar[17]*armttbar[109];
   armttbar[110]=armttbar[209] + 4*armttbar[20];
   armttbar[111]=armttbar[31]*armttbar[110];
   armttbar[110]=armttbar[1]*armttbar[110];
   armttbar[108]=2*armttbar[109] + 2*armttbar[108] + 5./3.*
   armttbar[111] + armttbar[110];
   armttbar[109]=5./3.*armttbar[114] + armttbar[113];
   armttbar[106]=armttbar[12]*armttbar[106];
   armttbar[106]=32./9.*armttbar[106] + 16./9.*armttbar[109] + 
   armttbar[205];
   armttbar[106]=MMZ*armttbar[106];
   armttbar[109]= - 4*armttbar[33] + armttbar[5];
   armttbar[110]= - armttbar[10] + 7./3. + armttbar[203];
   armttbar[110]=armttbar[15]*armttbar[110];
   armttbar[109]=2./3.*armttbar[110] + armttbar[168] + 1./3.*
   armttbar[109] + armttbar[18];
   armttbar[110]=armttbar[10] + 17./3. + armttbar[191];
   armttbar[110]=armttbar[124] + 4./9.*armttbar[110] - armttbar[21];
   armttbar[110]=armttbar[17]*armttbar[110];
   armttbar[109]=2./3.*armttbar[109] + armttbar[110];
   armttbar[110]= - armttbar[13] - 4 + armttbar[24];
   armttbar[111]=armttbar[21]*armttbar[204];
   armttbar[112]= - armttbar[10] + 1./3. + armttbar[203];
   armttbar[112]=armttbar[12]*armttbar[112];
   armttbar[110]=32./9.*armttbar[112] + 16./9.*armttbar[110] + 
   armttbar[111];
   armttbar[110]=MMZ*armttbar[110];
   armttbar[109]=8*armttbar[109] + armttbar[110];
   armttbar[109]=armttbar[27]*armttbar[109];
   armttbar[110]=armttbar[128] + 2*armttbar[127] + armttbar[15] - 
   armttbar[33] + armttbar[212];
   armttbar[110]=MMZ*armttbar[110];
   armttbar[111]=armttbar[251] + armttbar[17];
   armttbar[111]=armttbar[17]*armttbar[111];
   armttbar[110]=4*armttbar[111] + 1./3.*armttbar[110];
   armttbar[110]=armttbar[6]*armttbar[27]*armttbar[110];
   armttbar[106]=32./9.*armttbar[110] + 1./3.*armttbar[109] + 16./9.*
   armttbar[108] + armttbar[106];
   armttbar[106]=armttbar[6]*armttbar[106];
   armttbar[108]=MMt*armttbar[27]*armttbar[3]*armttbar[163];
   armttbar[105]=32./3.*armttbar[108] + 4*armttbar[105] + armttbar[106]
   ;

      mttbarret = armttbar[102] + armttbar[103] + armttbar[104] + 1./3.
      *armttbar[105] + armttbar[107];
      return mttbarret;
}
