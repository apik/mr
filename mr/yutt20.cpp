#include <tt.hpp>
std::complex<long double>
tt::my20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryutt[311], yuttret;

    aryutt[1]=double(nL + nH);
    aryutt[2]=pow(CW,-1);
    aryutt[3]=pow(MMH,-1);
    aryutt[4]=pow(SW,-1);
    aryutt[5]=Tsil::I2(0,0,MMZ,mu2);
    aryutt[6]=pow(MMt,-1);
    aryutt[7]=Tsil::I2(0,0,MMW,mu2);
    aryutt[8]=Tsil::B(MMH,MMt,MMt,mu2);
    aryutt[9]=pow(MMZ,-1);
    aryutt[10]=std::real(Tsil::B(0,0,MMZ,mu2));
    aryutt[11]=std::real(Tsil::B(0,0,MMW,mu2));
    aryutt[12]=Tsil::B(MMZ,MMt,MMt,mu2);
    aryutt[13]=Tsil::Beps(MMZ,MMt,MMt,mu2);
    aryutt[14]=Tsil::A(MMH,mu2);
    aryutt[15]=Tsil::A(MMZ,mu2);
    aryutt[16]=Tsil::A(MMW,mu2);
    aryutt[17]=Tsil::A(MMt,mu2);
    aryutt[18]=Tsil::Aeps(MMZ,mu2);
    aryutt[19]=Tsil::Aeps(MMW,mu2);
    aryutt[20]=Tsil::Aeps(MMt,mu2);
    aryutt[21]=std::real(Tsil::B(0,MMW,MMt,mu2));
    aryutt[22]=std::real(Tsil::B(0,0,MMt,mu2));
    aryutt[23]=std::real(Tsil::Beps(0,MMW,MMt,mu2));
    aryutt[24]=protZ0tW0->Uzxyv(0);
    aryutt[25]=protWt000->Tyzv(0);
    aryutt[26]=prot0W00->Uxzuv(0);
    aryutt[27]=double(nH);
    aryutt[28]=Tsil::B(MMt,MMt,MMH,mu2);
    aryutt[29]=Tsil::B(MMt,MMt,MMZ,mu2);
    aryutt[30]=Tsil::B(0,MMt,MMW,mu2);
    aryutt[31]=double(nL);
    aryutt[32]=Tsil::I2(MMH,MMt,MMt,mu2);
    aryutt[33]=Tsil::I2(MMZ,MMt,MMt,mu2);
    aryutt[34]=Tsil::I2(0,MMW,MMt,mu2);
    aryutt[35]=Tsil::B(MMH,MMH,MMH,mu2);
    aryutt[36]=Tsil::B(MMZ,MMH,MMZ,mu2);
    aryutt[37]=Tsil::B(MMZ,MMZ,MMH,mu2);
    aryutt[38]=Tsil::B(MMW,MMH,MMW,mu2);
    aryutt[39]=Tsil::B(MMW,MMZ,MMW,mu2);
    aryutt[40]=Tsil::B(MMW,MMW,MMH,mu2);
    aryutt[41]=Tsil::B(MMW,MMW,MMZ,mu2);
    aryutt[42]=Tsil::Aeps(MMH,mu2);
    aryutt[43]=protWt000->Uzxyv(0);
    aryutt[44]=prot0ttHt->Tuxv(0);
    aryutt[45]=prot0ttZt->Tuxv(0);
    aryutt[46]=double(boson);
    aryutt[47]=Tsil::I2(MMH,MMH,MMH,mu2);
    aryutt[48]=Tsil::I2(MMH,MMZ,MMZ,mu2);
    aryutt[49]=Tsil::I2(MMH,MMW,MMW,mu2);
    aryutt[50]=Tsil::I2(MMW,MMW,MMZ,mu2);
    aryutt[51]=Tsil::Beps(MMH,MMt,MMt,mu2);
    aryutt[52]=protWt000->M(0);
    aryutt[53]=prot0ttHt->M(0);
    aryutt[54]=prot0ttZt->M(0);
    aryutt[55]=prot0tt0t->M(0);
    aryutt[56]=prottH0H->Vxzuv(0);
    aryutt[57]=prottZ0Z->Vxzuv(0);
    aryutt[58]=protWt000->Txuv(0);
    aryutt[59]=protHHttH->M(0);
    aryutt[60]=protHZttZ->M(0);
    aryutt[61]=protHWt0W->M(0);
    aryutt[62]=protHttHt->M(0);
    aryutt[63]=protHttZt->M(0);
    aryutt[64]=protZZttH->M(0);
    aryutt[65]=protZWt0W->M(0);
    aryutt[66]=protZttZt->M(0);
    aryutt[67]=protZ0tW0->M(0);
    aryutt[68]=protWW00Z->M(0);
    aryutt[69]=protW00tW->M(0);
    aryutt[70]=prot00WW0->M(0);
    aryutt[71]=protHHttH->Uxzuv(0);
    aryutt[72]=protHZttZ->Uxzuv(0);
    aryutt[73]=protHWt0W->Uxzuv(0);
    aryutt[74]=protHttZt->Uuyxv(0);
    aryutt[75]=protHZttZ->Uyuzv(0);
    aryutt[76]=protZ0tW0->Uxzuv(0);
    aryutt[77]=protHZttZ->Uzxyv(0);
    aryutt[78]=protHWt0W->Uzxyv(0);
    aryutt[79]=protZZttH->Uzxyv(0);
    aryutt[80]=protZWt0W->Uzxyv(0);
    aryutt[81]=protW00tW->Uuyxv(0);
    aryutt[82]=protHWt0W->Uuyxv(0);
    aryutt[83]=protZWt0W->Uuyxv(0);
    aryutt[84]=protHZttZ->Txuv(0);
    aryutt[85]=protHWt0W->Txuv(0);
    aryutt[86]=protHttHt->Txuv(0);
    aryutt[87]=protZWt0W->Txuv(0);
    aryutt[88]=protHttZt->Tuxv(0);
    aryutt[89]=protHWt0W->Tvxu(0);
    aryutt[90]=protZWt0W->Tvxu(0);
    aryutt[91]=protHHttH->Suxv(0);
    aryutt[92]=protHZttZ->Suxv(0);
    aryutt[93]=protHZttZ->Svyz(0);
    aryutt[94]=protHWt0W->Svyz(0);
    aryutt[95]=protHWt0W->Suxv(0);
    aryutt[96]=protZ0tW0->Suxv(0);
    aryutt[97]=protHHttH->Uuyxv(0);
    aryutt[98]=1/(MMt - MMW);
    aryutt[99]=1/(4*MMt - MMZ);
    aryutt[100]=1/( - 4*MMW + MMH);
    aryutt[101]=1/( - 4*MMZ + MMH);
   aryutt[102]=55*aryutt[39];
   aryutt[103]= - 209./2.*aryutt[41] - 367./6. + aryutt[102];
   aryutt[104]=aryutt[14]*aryutt[101];
   aryutt[105]=1./2.*aryutt[104];
   aryutt[106]=5./16.*aryutt[12];
   aryutt[107]=1./2.*aryutt[36];
   aryutt[103]=aryutt[106] - 15./16.*aryutt[8] - 23./64.*aryutt[21] + 
   aryutt[105] + aryutt[107] + 1./4.*aryutt[103] - aryutt[38];
   aryutt[103]=aryutt[12]*aryutt[103];
   aryutt[108]= - 33*aryutt[41];
   aryutt[109]= - 133./8. + aryutt[108];
   aryutt[110]=27./32.*aryutt[21];
   aryutt[109]=aryutt[110] + 1./4.*aryutt[109] + aryutt[36];
   aryutt[109]=aryutt[99]*aryutt[109];
   aryutt[111]=1./2.*aryutt[98] + aryutt[101];
   aryutt[112]=aryutt[8]*aryutt[101];
   aryutt[113]=3*aryutt[112];
   aryutt[114]= - aryutt[12]*aryutt[101];
   aryutt[109]=3*aryutt[109] + aryutt[114] + 7*aryutt[111] + 
   aryutt[113];
   aryutt[109]=aryutt[15]*aryutt[109];
   aryutt[111]= - 3*aryutt[93];
   aryutt[115]=1./2.*aryutt[42];
   aryutt[116]=3*aryutt[20];
   aryutt[117]=aryutt[116] + 5*aryutt[18] + aryutt[111] + aryutt[115];
   aryutt[117]=aryutt[101]*aryutt[117];
   aryutt[118]= - 11*aryutt[41];
   aryutt[119]=11*aryutt[39];
   aryutt[120]=aryutt[14]*aryutt[100];
   aryutt[121]= - 25./32.*aryutt[21] + 1./4.*aryutt[120] - aryutt[38]
    + aryutt[118] - 2099./384. + aryutt[119];
   aryutt[121]=aryutt[21]*aryutt[121];
   aryutt[122]=3./2.*aryutt[50] + aryutt[94];
   aryutt[123]=5*aryutt[92];
   aryutt[122]=aryutt[123] + 3*aryutt[93] + 11*aryutt[122] - 27./4.*
   aryutt[96];
   aryutt[124]= - 5./2.*aryutt[14];
   aryutt[125]= - 1./4.*aryutt[32];
   aryutt[122]= - 237./8.*aryutt[16] + aryutt[124] - 207./16.*
   aryutt[20] - 337./32.*aryutt[18] - 55./2.*aryutt[19] + 9./8.*
   aryutt[33] + 27./16.*aryutt[34] + aryutt[125] + 1./2.*aryutt[122] - 
   aryutt[48];
   aryutt[122]=aryutt[99]*aryutt[122];
   aryutt[126]=9./2.*aryutt[19];
   aryutt[127]=1./2.*aryutt[49];
   aryutt[128]=aryutt[126] - 1./4.*aryutt[42] + aryutt[127] - 3*
   aryutt[94];
   aryutt[128]=aryutt[100]*aryutt[128];
   aryutt[129]= - 1./2.*aryutt[37];
   aryutt[130]= - 1./4.*aryutt[21];
   aryutt[131]=aryutt[130] + aryutt[129] - 1 - aryutt[40];
   aryutt[131]=aryutt[8]*aryutt[131];
   aryutt[132]= - 9./4.*aryutt[98] - aryutt[100];
   aryutt[132]=aryutt[21]*aryutt[132];
   aryutt[132]=aryutt[132] - 171./2.*aryutt[98] + 11*aryutt[100];
   aryutt[133]=aryutt[8]*aryutt[100];
   aryutt[132]=1./2.*aryutt[132] + 3*aryutt[133];
   aryutt[132]=aryutt[16]*aryutt[132];
   aryutt[133]=1./2.*aryutt[24];
   aryutt[134]=7./8.*aryutt[5] + 9./8.*aryutt[34] + 15*aryutt[94] - 7./
   4.*aryutt[96];
   aryutt[134]=aryutt[98]*aryutt[134];
   aryutt[135]= - aryutt[19]*aryutt[98];
   aryutt[136]=aryutt[18]*aryutt[98];
   aryutt[137]= - 43./8.*aryutt[98] + aryutt[100];
   aryutt[137]=aryutt[20]*aryutt[137];
   aryutt[138]= - 3*aryutt[100] - aryutt[101];
   aryutt[138]=1./8.*aryutt[14]*aryutt[138];
   aryutt[139]=3./2.*aryutt[40];
   aryutt[140]=3./4.*aryutt[37];
   aryutt[141]=pow(Pi,2);
   aryutt[103]=1./4.*aryutt[109] + 3./4.*aryutt[122] + 1./2.*
   aryutt[132] + 1./4.*aryutt[103] + 3./2.*aryutt[131] + 1./2.*
   aryutt[121] + aryutt[138] - 13./8.*aryutt[36] + 1./4.*aryutt[38] + 1.
   /4.*aryutt[117] + 561./32.*aryutt[41] + aryutt[140] + 3./2.*
   aryutt[137] + 3./16.*aryutt[136] - 231./16.*aryutt[39] + 1./2.*
   aryutt[128] - 347./256.*aryutt[13] + aryutt[139] + 123./8.*
   aryutt[135] + 1./2.*aryutt[134] + aryutt[133] - 1./2.*aryutt[25] + 9.
   /4.*aryutt[23] - 63./64.*aryutt[51] - aryutt[45] - 1./8.*aryutt[79]
    - 59./96.*aryutt[88] - 27./8.*aryutt[35] + 365./128.*aryutt[87] + 7.
   /64.*aryutt[72] + 7./64.*aryutt[74] + 1./16.*aryutt[75] + 1./4.*
   aryutt[77] + 195./64.*aryutt[90] + 3./32.*aryutt[141] + 1./4.*
   aryutt[82] - 73./96.*aryutt[89] + 31./32.*aryutt[80] + 1./2.*
   aryutt[43] - 41./256.*aryutt[76] + 1./2.*aryutt[78] + 1./8.*
   aryutt[73] + 28655./3072. - 3*aryutt[83];
   aryutt[109]=1./2. - aryutt[40];
   aryutt[117]=aryutt[109] + aryutt[129];
   aryutt[121]=aryutt[16]*aryutt[117];
   aryutt[117]=aryutt[15]*aryutt[117];
   aryutt[117]=aryutt[121] + 1./2.*aryutt[117];
   aryutt[117]=aryutt[3]*aryutt[117];
   aryutt[121]=3*aryutt[40];
   aryutt[122]=3./2.*aryutt[37];
   aryutt[128]=9./2.*aryutt[36];
   aryutt[131]=3./8.*aryutt[12];
   aryutt[117]=9*aryutt[117] + aryutt[131] + 3./8.*aryutt[21] + 
   aryutt[128] - 825./8.*aryutt[41] + aryutt[122] + 33./2.*aryutt[39]
    + 331./12. + aryutt[121];
   aryutt[117]=aryutt[3]*aryutt[117];
   aryutt[132]=3*aryutt[78];
   aryutt[134]=1./2.*aryutt[82];
   aryutt[135]= - 1./2.*aryutt[23];
   aryutt[137]= - 1./2.*aryutt[89];
   aryutt[142]=aryutt[135] - 3*aryutt[51] + aryutt[85] + aryutt[134] + 
   aryutt[137] - 13./2. + aryutt[132];
   aryutt[142]=aryutt[100]*aryutt[142];
   aryutt[143]=3./2.*aryutt[77];
   aryutt[144]= - 3./2.*aryutt[51];
   aryutt[145]= - 1./2.*aryutt[13];
   aryutt[146]=1./2.*aryutt[79];
   aryutt[147]=aryutt[145] + aryutt[84] + aryutt[144] + aryutt[146] - 1.
   /2.*aryutt[88] - 5 + aryutt[143];
   aryutt[147]=aryutt[101]*aryutt[147];
   aryutt[148]= - 1./8.*aryutt[8];
   aryutt[110]=aryutt[131] + aryutt[148] + aryutt[110] - 33./4.*
   aryutt[41] + aryutt[36];
   aryutt[110]=aryutt[12]*aryutt[110];
   aryutt[131]=4879./4. - 27*aryutt[76];
   aryutt[131]= - 81./4.*aryutt[87] + 1./2.*aryutt[72] + 1./2.*
   aryutt[74] - 3*aryutt[75] - 27./2.*aryutt[90] + 1./8.*aryutt[131] + 
   33*aryutt[80];
   aryutt[149]= - 1./8.*aryutt[51];
   aryutt[150]= - 27./32.*aryutt[21];
   aryutt[110]=aryutt[110] + aryutt[150] - aryutt[36] + 33./4.*
   aryutt[41] - 185./32.*aryutt[13] + aryutt[149] - aryutt[79] + 1./4.*
   aryutt[131] + 3*aryutt[88];
   aryutt[110]=aryutt[99]*aryutt[110];
   aryutt[131]=1./2.*aryutt[87] + 7./8. + aryutt[90];
   aryutt[151]=pow(aryutt[98],2);
   aryutt[131]=aryutt[151]*aryutt[131];
   aryutt[152]=aryutt[129] + 3./4. - aryutt[40];
   aryutt[153]=pow(aryutt[3],2);
   aryutt[152]=aryutt[153]*aryutt[152];
   aryutt[131]=1./8.*aryutt[131] + 3*aryutt[152];
   aryutt[131]=MMZ*aryutt[131];
   aryutt[152]=25 - 11./2.*aryutt[90];
   aryutt[154]= - 1./2.*aryutt[5];
   aryutt[155]=aryutt[96] + aryutt[154];
   aryutt[155]=aryutt[98]*aryutt[155];
   aryutt[152]=3./2.*aryutt[155] + 1./2.*aryutt[152] - aryutt[87];
   aryutt[152]=aryutt[98]*aryutt[152];
   aryutt[156]= - aryutt[61] - 7./2.*aryutt[65] - 7./8.*aryutt[68] + 
   aryutt[66];
   aryutt[157]= - 1./2.*aryutt[101];
   aryutt[158]= - aryutt[100] + aryutt[157];
   aryutt[158]=aryutt[8]*aryutt[158];
   aryutt[159]= - aryutt[16]*aryutt[151];
   aryutt[160]= - aryutt[15]*aryutt[151];
   aryutt[161]= - aryutt[21]*aryutt[100];
   aryutt[110]=3./2.*aryutt[131] + 1./2.*aryutt[117] + 3./16.*
   aryutt[160] + 3./16.*aryutt[110] + 9./32.*aryutt[159] + 1./8.*
   aryutt[114] + 3./4.*aryutt[158] + 1./8.*aryutt[161] + 1./4.*
   aryutt[147] + 1./4.*aryutt[142] + 1./8.*aryutt[152] - 1./4.*
   aryutt[60] - 1./8.*aryutt[64] + 1./2.*aryutt[156] - aryutt[67];
   aryutt[110]=MMZ*aryutt[110];
   aryutt[117]= - aryutt[42] + aryutt[14];
   aryutt[131]= - 1./2.*aryutt[12];
   aryutt[142]=aryutt[131] - 1 - aryutt[21];
   aryutt[147]=aryutt[16]*aryutt[142];
   aryutt[142]=aryutt[15]*aryutt[142];
   aryutt[152]= - aryutt[21]*aryutt[14];
   aryutt[156]= - aryutt[12]*aryutt[14];
   aryutt[142]=3./2.*aryutt[142] + 3*aryutt[147] + 1./2.*aryutt[156] + 
   3./2.*aryutt[117] + aryutt[152];
   aryutt[142]=aryutt[3]*aryutt[142];
   aryutt[147]= - 55*aryutt[39];
   aryutt[158]=457./6. + aryutt[147];
   aryutt[162]=77*aryutt[41];
   aryutt[158]=1./2.*aryutt[158] + aryutt[162];
   aryutt[158]=1./2.*aryutt[158] + aryutt[38];
   aryutt[163]= - aryutt[14]*aryutt[101];
   aryutt[158]=aryutt[131] - 27./16.*aryutt[21] + 1./4.*aryutt[163] + 1.
   /2.*aryutt[158] - aryutt[36];
   aryutt[158]=aryutt[12]*aryutt[158];
   aryutt[164]= - 11*aryutt[39];
   aryutt[165]=121./2.*aryutt[41] + 131./6. + aryutt[164];
   aryutt[166]= - aryutt[14]*aryutt[100];
   aryutt[167]=1./2.*aryutt[166];
   aryutt[165]= - 11./8.*aryutt[21] + aryutt[167] + 1./2.*aryutt[165]
    - aryutt[36];
   aryutt[165]=aryutt[21]*aryutt[165];
   aryutt[168]= - 3./2.*aryutt[76] - 229./8. - 11*aryutt[83];
   aryutt[168]= - 13./12.*aryutt[89] - 33./4.*aryutt[80] + 3./2.*
   aryutt[168] - 5*aryutt[43];
   aryutt[169]=3*aryutt[19] - aryutt[49] + aryutt[115];
   aryutt[169]=aryutt[100]*aryutt[169];
   aryutt[170]= - aryutt[48] + aryutt[115];
   aryutt[171]=aryutt[170] + 3*aryutt[18];
   aryutt[171]=aryutt[101]*aryutt[171];
   aryutt[172]=aryutt[100] + 1./2.*aryutt[101];
   aryutt[173]=aryutt[14]*aryutt[172];
   aryutt[174]=aryutt[21]*aryutt[100];
   aryutt[175]=aryutt[100] + aryutt[174];
   aryutt[175]=aryutt[16]*aryutt[175];
   aryutt[176]=aryutt[12]*aryutt[101];
   aryutt[177]=aryutt[101] + aryutt[176];
   aryutt[177]=aryutt[15]*aryutt[177];
   aryutt[178]=aryutt[17]*aryutt[3];
   aryutt[142]=3./4.*aryutt[178] + 1./2.*aryutt[142] + 1./2.*
   aryutt[177] + aryutt[175] + aryutt[158] + aryutt[165] + 3./2.*
   aryutt[173] + 1./2.*aryutt[171] + aryutt[169] + 63./8.*aryutt[13] - 
   5./2.*aryutt[24] + 39./4.*aryutt[23] - 1./2.*aryutt[84] + 3*
   aryutt[45] + aryutt[146] - 13./48.*aryutt[88] - aryutt[85] + 61./4.*
   aryutt[87] - 5./8.*aryutt[75] + 15*aryutt[90] - 1./8.*aryutt[141] + 
   1./2.*aryutt[168] + aryutt[82];
   aryutt[146]=3*aryutt[68];
   aryutt[158]=aryutt[146] - aryutt[66];
   aryutt[165]=3*aryutt[65];
   aryutt[158]=1./2.*aryutt[158] + aryutt[165];
   aryutt[168]= - aryutt[82] - 5 + aryutt[89];
   aryutt[169]=1./2.*aryutt[23];
   aryutt[168]=aryutt[169] + 1./2.*aryutt[168] - aryutt[85];
   aryutt[168]=1./2.*aryutt[100]*aryutt[168];
   aryutt[171]= - aryutt[79] - 5 + aryutt[88];
   aryutt[173]=1./2.*aryutt[13];
   aryutt[171]=aryutt[173] + 1./2.*aryutt[171] - aryutt[84];
   aryutt[171]=aryutt[101]*aryutt[171];
   aryutt[175]=1./4.*aryutt[174];
   aryutt[179]= - aryutt[21] + aryutt[131];
   aryutt[179]=aryutt[3]*aryutt[179];
   aryutt[158]=3./8.*aryutt[179] + 1./8.*aryutt[176] + aryutt[175] + 1./
   4.*aryutt[171] + aryutt[168] + 1./2.*aryutt[158] + aryutt[67];
   aryutt[158]=MMZ*aryutt[158];
   aryutt[142]=1./4.*aryutt[142] + aryutt[158];
   aryutt[142]=MMZ*aryutt[142];
   aryutt[119]=191./2. + aryutt[119];
   aryutt[158]=aryutt[16]*aryutt[100];
   aryutt[179]=11*aryutt[41];
   aryutt[180]= - 1./2.*aryutt[38];
   aryutt[119]=1./2.*aryutt[158] + 61./32.*aryutt[12] + 109./24.*
   aryutt[21] + aryutt[167] + aryutt[180] + 1./8.*aryutt[119] + 
   aryutt[179];
   aryutt[119]=aryutt[16]*aryutt[119];
   aryutt[147]= - 233./2. + aryutt[147];
   aryutt[147]=1./2.*aryutt[147] + aryutt[162];
   aryutt[147]=1./2.*aryutt[147] + aryutt[38];
   aryutt[158]= - 1./16.*aryutt[21];
   aryutt[162]=aryutt[15]*aryutt[101];
   aryutt[147]=1./2.*aryutt[162] - 67./48.*aryutt[12] + aryutt[158] + 1.
   /2.*aryutt[163] + 1./2.*aryutt[147] - aryutt[36];
   aryutt[147]=aryutt[15]*aryutt[147];
   aryutt[167]=1./2.*aryutt[48];
   aryutt[181]=aryutt[21]*aryutt[14];
   aryutt[182]=1./8.*aryutt[181];
   aryutt[183]=aryutt[12]*aryutt[14];
   aryutt[184]=1./4.*aryutt[183] + aryutt[182] + 5./8.*aryutt[14] + 281.
   /24.*aryutt[20] + 37./2.*aryutt[18] + 413./8.*aryutt[19] + 7./16.*
   aryutt[42] - 3./2.*aryutt[5] - 1./16.*aryutt[33] - 17./4.*aryutt[34]
    + aryutt[167] - 7./12.*aryutt[92] - 5./8.*aryutt[93] + 93./4.*
   aryutt[96] - 57./8.*aryutt[94] - 171./8.*aryutt[50] - 7./6.*
   aryutt[95] + aryutt[49];
   aryutt[185]= - aryutt[16] - 1./4.*aryutt[15];
   aryutt[185]=aryutt[15]*aryutt[185];
   aryutt[186]=pow(aryutt[16],2);
   aryutt[185]= - aryutt[186] + aryutt[185];
   aryutt[185]=aryutt[3]*aryutt[185];
   aryutt[187]= - 17 + 5*aryutt[39];
   aryutt[187]=1./2.*aryutt[187] - 7*aryutt[41];
   aryutt[187]=11./2.*aryutt[187] - aryutt[38];
   aryutt[188]=1./2.*aryutt[15];
   aryutt[189]=aryutt[16] + aryutt[188];
   aryutt[190]=aryutt[3]*aryutt[189];
   aryutt[187]=3./4.*aryutt[190] + 11./16.*aryutt[12] - 45./16.*
   aryutt[21] + 1./2.*aryutt[187] + aryutt[36];
   aryutt[187]=aryutt[17]*aryutt[187];
   aryutt[119]=1./2.*aryutt[187] + 3./4.*aryutt[185] + 1./2.*
   aryutt[147] + 1./2.*aryutt[184] + aryutt[119];
   aryutt[147]= - 77./8.*aryutt[141] - 7./2. + aryutt[76];
   aryutt[184]= - 1./4.*aryutt[12];
   aryutt[187]=aryutt[184] - 5 - aryutt[21];
   aryutt[187]=aryutt[12]*aryutt[187];
   aryutt[191]=1./4.*aryutt[75];
   aryutt[192]=pow(aryutt[21],2);
   aryutt[147]=1./4.*aryutt[187] - 1./4.*aryutt[192] - 3./4.*aryutt[13]
    - 5./4.*aryutt[45] + aryutt[87] + aryutt[191] + 1./2.*aryutt[147]
    + aryutt[90];
   aryutt[147]=MMZ*aryutt[147];
   aryutt[187]=1./2.*aryutt[19];
   aryutt[192]= - 1./2.*aryutt[34];
   aryutt[193]=5./4.*aryutt[18] + aryutt[187] - 3./2.*aryutt[33] + 
   aryutt[192] + aryutt[96] - 1./2.*aryutt[93];
   aryutt[194]=1./2.*aryutt[12];
   aryutt[195]=aryutt[194] - 3 - aryutt[21];
   aryutt[195]=aryutt[16]*aryutt[195];
   aryutt[196]=aryutt[194] + 5 - aryutt[21];
   aryutt[196]=aryutt[15]*aryutt[196];
   aryutt[197]= - 21./2.*aryutt[12] + 7 + aryutt[21];
   aryutt[197]=aryutt[17]*aryutt[197];
   aryutt[147]=1./4.*aryutt[147] + 1./16.*aryutt[197] + 1./16.*
   aryutt[196] + 1./8.*aryutt[195] + 1./4.*aryutt[193] + aryutt[20];
   aryutt[147]=MMZ*aryutt[147];
   aryutt[193]= - aryutt[16] - 9./4.*aryutt[15];
   aryutt[193]=aryutt[15]*aryutt[193];
   aryutt[195]= - 49./2.*aryutt[15];
   aryutt[196]=15./4.*aryutt[17] - aryutt[16] + aryutt[195];
   aryutt[196]=aryutt[17]*aryutt[196];
   aryutt[193]=aryutt[196] - aryutt[186] + aryutt[193];
   aryutt[196]=pow(MMZ,4);
   aryutt[197]= - aryutt[6]*aryutt[196]*aryutt[141];
   aryutt[198]=aryutt[141]*pow(MMZ,3);
   aryutt[197]=21*aryutt[198] + 5*aryutt[197];
   aryutt[197]=aryutt[6]*aryutt[197];
   aryutt[147]=1./16.*aryutt[197] + 1./16.*aryutt[193] + aryutt[147];
   aryutt[147]=aryutt[6]*aryutt[147];
   aryutt[193]= - 1./4.*aryutt[79];
   aryutt[197]= - 1./2.*aryutt[82];
   aryutt[198]= - 1./8.*aryutt[13] + aryutt[169] + 3./8.*aryutt[84] - 9.
   /8.*aryutt[51] + aryutt[193] + 11./8.*aryutt[88] + 3./8.*aryutt[72]
    + 3./8.*aryutt[74] + aryutt[197] + 5./4.*aryutt[89] - 1 + 3./4.*
   aryutt[73];
   aryutt[199]=1./8. + aryutt[36];
   aryutt[200]=1./8.*aryutt[8];
   aryutt[199]=1./3.*aryutt[199] + aryutt[200];
   aryutt[199]=aryutt[12]*aryutt[199];
   aryutt[201]=aryutt[36] + 1./4. + aryutt[38];
   aryutt[201]=aryutt[21]*aryutt[201];
   aryutt[202]=aryutt[8]*aryutt[21];
   aryutt[198]=aryutt[199] + 1./4.*aryutt[202] + 1./2.*aryutt[198] + 1./
   3.*aryutt[201];
   aryutt[198]=MMH*aryutt[198];
   aryutt[119]=1./2.*aryutt[147] + 1./4.*aryutt[198] + 1./2.*
   aryutt[119] + aryutt[142];
   aryutt[119]=aryutt[6]*aryutt[119];
   aryutt[142]=3*aryutt[51];
   aryutt[147]= - 3*aryutt[78];
   aryutt[198]=aryutt[135] + aryutt[142] + aryutt[85] + aryutt[134] + 
   aryutt[137] + 59./8. + aryutt[147];
   aryutt[198]=aryutt[100]*aryutt[198];
   aryutt[199]=1./2.*aryutt[37];
   aryutt[201]=aryutt[199] + 5./2. + aryutt[40];
   aryutt[201]=aryutt[8]*aryutt[201];
   aryutt[202]=1./4.*aryutt[21];
   aryutt[201]=aryutt[201] + aryutt[202] + aryutt[129] - 3./2. - 
   aryutt[40];
   aryutt[201]=aryutt[3]*aryutt[201];
   aryutt[146]=aryutt[146] - 5./2.*aryutt[66];
   aryutt[146]=1./2.*aryutt[146] + aryutt[165];
   aryutt[165]=1./2.*aryutt[64];
   aryutt[203]=aryutt[51] + 21./8. - aryutt[77];
   aryutt[203]=aryutt[101]*aryutt[203];
   aryutt[204]=3./2.*aryutt[203];
   aryutt[205]=1./2.*aryutt[161];
   aryutt[172]=aryutt[8]*aryutt[172];
   aryutt[146]=3*aryutt[201] + 3*aryutt[172] + aryutt[205] + 
   aryutt[204] + aryutt[198] + aryutt[60] + aryutt[165] + 3./2.*
   aryutt[67] - aryutt[63] + 1./2.*aryutt[146] + aryutt[61];
   aryutt[146]=MMt*aryutt[146];
   aryutt[172]=7*aryutt[36];
   aryutt[198]= - aryutt[39] + aryutt[41];
   aryutt[198]= - aryutt[36] + 33./4.*aryutt[198] + aryutt[38];
   aryutt[201]=aryutt[12]*aryutt[198];
   aryutt[206]= - 7*aryutt[38];
   aryutt[207]=aryutt[39] - aryutt[41];
   aryutt[208]=aryutt[21]*aryutt[198];
   aryutt[201]=aryutt[201] + aryutt[208] + aryutt[172] + 99./4.*
   aryutt[207] + aryutt[206];
   aryutt[209]= - 17./24.*aryutt[12];
   aryutt[207]=aryutt[36] + 33./4.*aryutt[207] - aryutt[38];
   aryutt[210]=aryutt[209] + aryutt[207] - 17./12.*aryutt[21];
   aryutt[210]=aryutt[16]*aryutt[210];
   aryutt[211]=17./12.*aryutt[12];
   aryutt[212]=aryutt[211] + aryutt[207] + 17./6.*aryutt[21];
   aryutt[212]=aryutt[15]*aryutt[212];
   aryutt[213]=aryutt[21]*aryutt[207];
   aryutt[214]=aryutt[12]*aryutt[207];
   aryutt[213]=aryutt[213] + 1./2.*aryutt[214];
   aryutt[213]=MMZ*aryutt[213];
   aryutt[214]=aryutt[38] - aryutt[36];
   aryutt[215]=aryutt[12]*aryutt[214];
   aryutt[216]=aryutt[21]*aryutt[214];
   aryutt[217]=aryutt[216] + 1./2.*aryutt[215];
   aryutt[217]=MMH*aryutt[217];
   aryutt[218]=aryutt[17]*aryutt[198];
   aryutt[210]=1./3.*aryutt[217] + aryutt[213] + 1./2.*aryutt[218] + 
   aryutt[210] + 1./2.*aryutt[212];
   aryutt[210]=aryutt[6]*aryutt[210];
   aryutt[212]=99*aryutt[41] + 17 - 99*aryutt[39];
   aryutt[213]=3*aryutt[38];
   aryutt[217]= - 3*aryutt[36];
   aryutt[212]=aryutt[217] + 1./4.*aryutt[212] + aryutt[213];
   aryutt[212]=aryutt[16]*aryutt[212];
   aryutt[219]=99./2.*aryutt[41] - 17 - 99./2.*aryutt[39];
   aryutt[219]=aryutt[217] + 1./2.*aryutt[219] + aryutt[213];
   aryutt[219]=aryutt[15]*aryutt[219];
   aryutt[212]=aryutt[212] + 1./2.*aryutt[219];
   aryutt[212]=aryutt[3]*aryutt[212];
   aryutt[219]=MMZ*aryutt[3]*aryutt[198];
   aryutt[201]=1./2.*aryutt[210] + 3*aryutt[219] + 1./4.*aryutt[201] + 
   aryutt[212];
   aryutt[210]=pow(aryutt[4],2);
   aryutt[201]=aryutt[210]*aryutt[201];
   aryutt[212]=aryutt[140] + 5 + aryutt[139];
   aryutt[212]=aryutt[14]*aryutt[212];
   aryutt[219]= - 1./2.*aryutt[48];
   aryutt[220]=aryutt[219] - aryutt[49] + 11*aryutt[50];
   aryutt[220]=9*aryutt[220] + 25./2.*aryutt[42];
   aryutt[221]=1./8.*aryutt[156];
   aryutt[212]=aryutt[221] + 1./8.*aryutt[152] + 1./2.*aryutt[212] - 45.
   /2.*aryutt[18] + 1./4.*aryutt[220] - 45*aryutt[19];
   aryutt[220]=aryutt[199] + 11./2.*aryutt[39] - 85./4. + aryutt[40];
   aryutt[222]=1./16.*aryutt[21];
   aryutt[118]=1./16.*aryutt[12] + aryutt[222] + 1./2.*aryutt[220] + 
   aryutt[118];
   aryutt[118]=aryutt[16]*aryutt[118];
   aryutt[220]=55./4.*aryutt[39];
   aryutt[223]=1./8.*aryutt[12];
   aryutt[224]=aryutt[223] + 1./8.*aryutt[21] + aryutt[36] - aryutt[38]
    - 121./4.*aryutt[41] + aryutt[199] + aryutt[220] - 93./8. + 
   aryutt[40];
   aryutt[224]=aryutt[15]*aryutt[224];
   aryutt[118]=9./8.*aryutt[185] + 3./4.*aryutt[224] + 1./2.*
   aryutt[212] + 3*aryutt[118];
   aryutt[118]=aryutt[3]*aryutt[118];
   aryutt[185]=7./2.*aryutt[98];
   aryutt[212]=aryutt[21]*aryutt[98];
   aryutt[157]=3./4.*aryutt[212] + aryutt[157] + aryutt[185] - 
   aryutt[100];
   aryutt[212]=13 + aryutt[179];
   aryutt[150]=aryutt[150] + 3./4.*aryutt[212] - aryutt[36];
   aryutt[150]=aryutt[99]*aryutt[150];
   aryutt[150]=1./2.*aryutt[157] + aryutt[150];
   aryutt[157]=1./2.*aryutt[40];
   aryutt[212]=1./4.*aryutt[37];
   aryutt[224]=aryutt[212] + 1 + aryutt[157];
   aryutt[224]=aryutt[3]*aryutt[224];
   aryutt[150]=1./2.*aryutt[150] + aryutt[224];
   aryutt[150]=aryutt[17]*aryutt[150];
   aryutt[224]= - 3*aryutt[72] + 119./4. - 3*aryutt[74];
   aryutt[225]= - aryutt[36] + 3./8.*aryutt[8];
   aryutt[225]=aryutt[12]*aryutt[225];
   aryutt[224]=aryutt[225] + aryutt[36] - 5./8.*aryutt[13] + 15./2.*
   aryutt[84] + 3./8.*aryutt[51] + aryutt[79] + 1./8.*aryutt[224] - 5*
   aryutt[88];
   aryutt[224]=aryutt[99]*aryutt[224];
   aryutt[225]=1./4.*aryutt[64];
   aryutt[226]=1./2.*aryutt[60];
   aryutt[224]=1./4.*aryutt[224] + aryutt[226] + aryutt[61] + 
   aryutt[225];
   aryutt[224]=MMH*aryutt[224];
   aryutt[103]=1./2.*aryutt[201] + 1./4.*aryutt[146] + aryutt[119] + 1./
   4.*aryutt[224] + aryutt[110] + 3./2.*aryutt[150] + 1./2.*aryutt[103]
    + aryutt[118];
   aryutt[103]=aryutt[210]*aryutt[103];
   aryutt[110]=aryutt[169] + aryutt[142] - aryutt[85] + aryutt[197] + 1.
   /2.*aryutt[89] + 13./2. + aryutt[147];
   aryutt[110]=aryutt[100]*aryutt[110];
   aryutt[118]= - 3./2.*aryutt[37];
   aryutt[119]=aryutt[118] - 1./2. + aryutt[121];
   aryutt[119]=aryutt[16]*aryutt[119];
   aryutt[142]= - 3*aryutt[37];
   aryutt[146]=1 + aryutt[142];
   aryutt[150]=aryutt[15]*aryutt[146];
   aryutt[119]=aryutt[119] + 1./2.*aryutt[150];
   aryutt[119]=aryutt[3]*aryutt[119];
   aryutt[169]=205./2.*aryutt[41] + aryutt[37] + 151./4.*aryutt[39] - 
   877./18. - aryutt[40];
   aryutt[201]=3*aryutt[36];
   aryutt[119]=3./2.*aryutt[119] + 9./8.*aryutt[12] + aryutt[158] + 
   aryutt[201] + 1./2.*aryutt[169] - 2*aryutt[38];
   aryutt[119]=aryutt[3]*aryutt[119];
   aryutt[169]= - 11./6.*aryutt[88];
   aryutt[224]= - 11./6.*aryutt[13] + 11./3.*aryutt[84] + aryutt[144]
    + 11./6.*aryutt[79] + aryutt[169] - 19./3. + aryutt[143];
   aryutt[224]=1./2.*aryutt[101]*aryutt[224];
   aryutt[227]= - 148303./12. + 225*aryutt[76];
   aryutt[227]=675./4.*aryutt[87] - 5./2.*aryutt[72] - 5./2.*aryutt[74]
    + 39*aryutt[75] + 531./4.*aryutt[90] + 1./8.*aryutt[227] - 297*
   aryutt[80];
   aryutt[228]=5*aryutt[36];
   aryutt[229]=225./32.*aryutt[21];
   aryutt[227]=aryutt[229] + aryutt[228] - 297./4.*aryutt[41] + 4585./
   96.*aryutt[13] + 5./8.*aryutt[51] + 5*aryutt[79] + 1./4.*aryutt[227]
    - 15*aryutt[88];
   aryutt[230]= - 5./8.*aryutt[36];
   aryutt[231]= - 39./64.*aryutt[12] + 5./64.*aryutt[8] - 225./256.*
   aryutt[21] + aryutt[230] + 2./3. + 297./32.*aryutt[41];
   aryutt[231]=aryutt[12]*aryutt[231];
   aryutt[227]=1./8.*aryutt[227] + aryutt[231];
   aryutt[227]=aryutt[99]*aryutt[227];
   aryutt[231]= - 1./8.*aryutt[87] - 1./2. - aryutt[90];
   aryutt[231]=aryutt[151]*aryutt[231];
   aryutt[232]=aryutt[129] - 3./4. + 2*aryutt[40];
   aryutt[232]=aryutt[153]*aryutt[232];
   aryutt[231]=1./4.*aryutt[231] + 3*aryutt[232];
   aryutt[231]=MMZ*aryutt[231];
   aryutt[232]= - 17*aryutt[87] - 197 + 11*aryutt[90];
   aryutt[233]=1./2.*aryutt[5];
   aryutt[234]= - aryutt[96] + aryutt[233];
   aryutt[234]=aryutt[98]*aryutt[234];
   aryutt[232]=1./6.*aryutt[232] + aryutt[234];
   aryutt[232]=aryutt[98]*aryutt[232];
   aryutt[235]=aryutt[100] - aryutt[101];
   aryutt[235]=aryutt[8]*aryutt[235];
   aryutt[236]=11./12.*aryutt[114];
   aryutt[237]=aryutt[16]*aryutt[151];
   aryutt[238]=aryutt[15]*aryutt[151];
   aryutt[239]=1./2.*aryutt[61];
   aryutt[110]=aryutt[231] + aryutt[119] + 1./16.*aryutt[238] + 
   aryutt[227] + 3./32.*aryutt[237] + aryutt[236] + 3./4.*aryutt[235]
    + 1./8.*aryutt[174] + aryutt[224] + 1./4.*aryutt[110] + 1./16.*
   aryutt[232] + 1./6.*aryutt[60] + 1./12.*aryutt[64] + aryutt[67] + 
   aryutt[239] + aryutt[65] - 7./6.*aryutt[66] + 11./12.*aryutt[68] + 2.
   /3.*aryutt[54] + 1./3.*aryutt[69] - 1./3.*aryutt[52] - 4./9.*
   aryutt[57] - 1./4.*aryutt[70];
   aryutt[110]=MMZ*aryutt[110];
   aryutt[119]=aryutt[82] + 5 - aryutt[89];
   aryutt[119]=aryutt[135] + 1./2.*aryutt[119] + aryutt[85];
   aryutt[119]=aryutt[100]*aryutt[119];
   aryutt[174]=aryutt[79] + 5 - aryutt[88];
   aryutt[174]=aryutt[145] + 1./2.*aryutt[174] + aryutt[84];
   aryutt[174]=aryutt[101]*aryutt[174];
   aryutt[227]= - 2./3.*aryutt[54];
   aryutt[231]=aryutt[227] + 2*aryutt[69] + 2./3.*aryutt[52] + 4./3.*
   aryutt[57] + 1./2.*aryutt[70];
   aryutt[232]=3./4.*aryutt[12];
   aryutt[235]=aryutt[21] + aryutt[232];
   aryutt[235]=aryutt[3]*aryutt[235];
   aryutt[119]=1./2.*aryutt[235] + 1./12.*aryutt[114] + aryutt[205] + 1.
   /6.*aryutt[174] + aryutt[119] - 2*aryutt[67] - 4*aryutt[65] + 5./6.*
   aryutt[66] + 1./3.*aryutt[231] - 3./2.*aryutt[68];
   aryutt[119]=MMZ*aryutt[119];
   aryutt[205]=aryutt[42] - aryutt[14];
   aryutt[231]=3*aryutt[21];
   aryutt[235]=5./2.*aryutt[12] + 1 + aryutt[231];
   aryutt[235]=aryutt[16]*aryutt[235];
   aryutt[237]=3 + aryutt[194];
   aryutt[237]=aryutt[15]*aryutt[237];
   aryutt[238]=1./4.*aryutt[181];
   aryutt[235]=1./4.*aryutt[237] + 1./4.*aryutt[235] + 1./12.*
   aryutt[183] + 1./3.*aryutt[205] + aryutt[238];
   aryutt[235]=aryutt[3]*aryutt[235];
   aryutt[237]=17./3.*aryutt[38];
   aryutt[240]= - 5./3.*aryutt[36];
   aryutt[241]=151./8.*aryutt[21] + aryutt[104] + aryutt[240] + 
   aryutt[237] - 685./4.*aryutt[41] - 367./2. - 35*aryutt[39];
   aryutt[241]=1./8.*aryutt[241] + aryutt[12];
   aryutt[241]=aryutt[12]*aryutt[241];
   aryutt[242]= - 3*aryutt[19];
   aryutt[243]= - 1./2.*aryutt[42];
   aryutt[244]=aryutt[242] + aryutt[49] + aryutt[243];
   aryutt[244]=aryutt[100]*aryutt[244];
   aryutt[245]= - 119*aryutt[41] - 361./3. - 59./2.*aryutt[39];
   aryutt[246]=1./2.*aryutt[120];
   aryutt[245]=11./8.*aryutt[21] + aryutt[246] + 1./3.*aryutt[245] + 
   aryutt[38];
   aryutt[245]=aryutt[21]*aryutt[245];
   aryutt[247]= - 25./6.*aryutt[58] + 2393./64. + 2*aryutt[81];
   aryutt[243]=aryutt[48] + aryutt[243];
   aryutt[243]=1./3.*aryutt[243] - aryutt[18];
   aryutt[243]=1./4.*aryutt[101]*aryutt[243];
   aryutt[248]= - aryutt[100] + aryutt[161];
   aryutt[248]=aryutt[16]*aryutt[248];
   aryutt[249]= - aryutt[101] + aryutt[114];
   aryutt[249]=1./12.*aryutt[15]*aryutt[249];
   aryutt[119]=aryutt[119] + 13./8.*aryutt[178] + 1./2.*aryutt[235] + 
   aryutt[249] + 1./4.*aryutt[248] + 1./3.*aryutt[241] + 1./4.*
   aryutt[245] + aryutt[138] + aryutt[243] + 1./4.*aryutt[244] - 317./
   144.*aryutt[13] + 7./8.*aryutt[24] - 2./3.*aryutt[25] - 47./12.*
   aryutt[23] + 1./12.*aryutt[84] - 11./6.*aryutt[45] - 1./12.*
   aryutt[79] + 13./288.*aryutt[88] + 1./4.*aryutt[85] - 29./6.*
   aryutt[87] + 5./48.*aryutt[75] - 631./96.*aryutt[90] + 101./96.*
   aryutt[141] - 1./4.*aryutt[82] + 13./96.*aryutt[89] + 37./24.*
   aryutt[80] + 2./3.*aryutt[43] + 5./24.*aryutt[76] + 1./3.*
   aryutt[247] + 7./2.*aryutt[83];
   aryutt[119]=MMZ*aryutt[119];
   aryutt[138]=aryutt[104] - 5./6.*aryutt[36] + 17./6.*aryutt[38] - 685.
   /8.*aryutt[41] - 41./3. - 35./2.*aryutt[39];
   aryutt[235]= - aryutt[15]*aryutt[101];
   aryutt[138]=1./3.*aryutt[235] - 655./216.*aryutt[12] + 1./3.*
   aryutt[138] - 75./16.*aryutt[21];
   aryutt[138]=aryutt[15]*aryutt[138];
   aryutt[241]= - 17./3.*aryutt[38];
   aryutt[244]=5./3.*aryutt[36];
   aryutt[245]=7./12.*aryutt[12] - 365./24.*aryutt[21] + aryutt[244] + 
   aryutt[241] + 685./4.*aryutt[41] + 3593./12. + 35*aryutt[39];
   aryutt[247]=9*aryutt[16];
   aryutt[248]=aryutt[247] + 5*aryutt[15];
   aryutt[248]=aryutt[3]*aryutt[248];
   aryutt[245]=1./3.*aryutt[245] + 3./2.*aryutt[248];
   aryutt[245]=aryutt[17]*aryutt[245];
   aryutt[248]=65*aryutt[93] - 239*aryutt[96] + 263*aryutt[50] + 197./3.
   *aryutt[94];
   aryutt[248]=1./2.*aryutt[248] + 7./3.*aryutt[92];
   aryutt[248]=1./2.*aryutt[248] - aryutt[48];
   aryutt[248]=35./48.*aryutt[33] + 1./2.*aryutt[248] + 5./3.*
   aryutt[34];
   aryutt[250]= - 1./4.*aryutt[5];
   aryutt[251]= - 13./144.*aryutt[42];
   aryutt[252]= - 1./36.*aryutt[14];
   aryutt[253]=1./12.*aryutt[156];
   aryutt[254]= - 943./2. - 85*aryutt[39];
   aryutt[254]= - 1303./36.*aryutt[12] - 113./6.*aryutt[21] + 1./4.*
   aryutt[254] - 53*aryutt[41];
   aryutt[254]=aryutt[16]*aryutt[254];
   aryutt[255]=aryutt[3]*aryutt[15]*aryutt[189];
   aryutt[138]=1./4.*aryutt[245] + 1./4.*aryutt[255] + 1./2.*
   aryutt[138] + 1./6.*aryutt[254] + aryutt[253] + aryutt[252] - 253./
   72.*aryutt[20] - 149./8.*aryutt[18] - 3803./144.*aryutt[19] + 
   aryutt[251] + 1./3.*aryutt[248] + aryutt[250];
   aryutt[245]=1./9. + aryutt[222];
   aryutt[245]=aryutt[21]*aryutt[245];
   aryutt[106]=aryutt[106] + 59./12. + aryutt[21];
   aryutt[106]=aryutt[12]*aryutt[106];
   aryutt[106]=1./12.*aryutt[106] + aryutt[245] + 13./48.*aryutt[13] + 
   59./144.*aryutt[45] - 1./3.*aryutt[87] - 5./48.*aryutt[75] - 19./48.
   *aryutt[90] + 277./576.*aryutt[141] - 1./6.*aryutt[76] + 3./4. + 1./
   9.*aryutt[58];
   aryutt[106]=MMZ*aryutt[106];
   aryutt[245]= - 109./24.*aryutt[20] - 25./16.*aryutt[18] - 5./16.*
   aryutt[19] + 37./24.*aryutt[33] + 5./6.*aryutt[34] - aryutt[96] + 5./
   8.*aryutt[93];
   aryutt[248]= - 5./6.*aryutt[12];
   aryutt[254]=aryutt[248] + 47./9. + aryutt[21];
   aryutt[254]=aryutt[16]*aryutt[254];
   aryutt[255]= - 5./8.*aryutt[12] - 77./12. + aryutt[21];
   aryutt[255]=aryutt[15]*aryutt[255];
   aryutt[256]=203./8.*aryutt[12] + 119./4. - 11*aryutt[21];
   aryutt[256]=aryutt[17]*aryutt[256];
   aryutt[106]=aryutt[106] + 1./36.*aryutt[256] + 1./12.*aryutt[255] + 
   1./3.*aryutt[245] + 1./16.*aryutt[254];
   aryutt[106]=MMZ*aryutt[106];
   aryutt[245]=pow(MMZ,2);
   aryutt[254]= - aryutt[245]*aryutt[141];
   aryutt[255]=pow(aryutt[17],2);
   aryutt[254]=1./3.*aryutt[255] + 63./32.*aryutt[254];
   aryutt[254]=MMZ*aryutt[254];
   aryutt[256]=aryutt[6]*aryutt[196]*aryutt[141];
   aryutt[254]=aryutt[254] + 5./8.*aryutt[256];
   aryutt[254]=aryutt[6]*aryutt[254];
   aryutt[256]=1./3.*aryutt[16];
   aryutt[257]=aryutt[256] + 3./2.*aryutt[15];
   aryutt[257]=aryutt[15]*aryutt[257];
   aryutt[258]= - 29./9.*aryutt[16] + 11*aryutt[15];
   aryutt[258]=5*aryutt[258] - 11./6.*aryutt[17];
   aryutt[258]=aryutt[17]*aryutt[258];
   aryutt[257]=5*aryutt[257] + aryutt[258];
   aryutt[106]=aryutt[254] + 1./32.*aryutt[257] + aryutt[106];
   aryutt[106]=aryutt[6]*aryutt[106];
   aryutt[254]= - 1./2.*aryutt[74];
   aryutt[257]= - 1./2.*aryutt[72];
   aryutt[258]=1./3.*aryutt[79];
   aryutt[169]=aryutt[258] + aryutt[169] + aryutt[257] + 13./9. + 
   aryutt[254];
   aryutt[241]=aryutt[244] + 11./2. + aryutt[241];
   aryutt[244]= - 1./2.*aryutt[8];
   aryutt[241]=1./3.*aryutt[241] + aryutt[244];
   aryutt[241]=aryutt[12]*aryutt[241];
   aryutt[259]=1./8.*aryutt[51];
   aryutt[169]=1./6.*aryutt[241] + 1./24.*aryutt[13] - 1./8.*aryutt[84]
    + aryutt[259] + 1./4.*aryutt[169] + 1./3.*aryutt[45];
   aryutt[169]=1./4.*MMH*aryutt[169];
   aryutt[106]=aryutt[106] + aryutt[169] + 1./2.*aryutt[138] + 
   aryutt[119];
   aryutt[106]=aryutt[6]*aryutt[106];
   aryutt[119]=13./2.*aryutt[39];
   aryutt[138]= - 3./2.*aryutt[38];
   aryutt[241]=3./16.*aryutt[21];
   aryutt[260]=11./8.*aryutt[12];
   aryutt[261]=aryutt[260] + aryutt[241] + aryutt[128] + aryutt[138] + 
   179./8.*aryutt[41] + aryutt[122] + aryutt[119] + 79./3. + 
   aryutt[139];
   aryutt[261]=aryutt[15]*aryutt[261];
   aryutt[262]=19*aryutt[39];
   aryutt[263]=759./2. + aryutt[262];
   aryutt[264]=3*aryutt[37];
   aryutt[263]=1./2.*aryutt[263] + aryutt[264];
   aryutt[265]=53*aryutt[41];
   aryutt[263]=19./8.*aryutt[12] + 1./2.*aryutt[263] + aryutt[265];
   aryutt[263]=aryutt[16]*aryutt[263];
   aryutt[266]= - 11*aryutt[50] - 3./2.*aryutt[48];
   aryutt[267]=5./3.*aryutt[42];
   aryutt[268]=5./3. + aryutt[264];
   aryutt[268]=1./8.*aryutt[14]*aryutt[268];
   aryutt[269]=13./24.*aryutt[183];
   aryutt[270]= - 1./2.*aryutt[15];
   aryutt[271]= - aryutt[16] + aryutt[270];
   aryutt[272]=9./8.*aryutt[3]*aryutt[15]*aryutt[271];
   aryutt[261]=aryutt[272] + 1./2.*aryutt[261] + 1./2.*aryutt[263] + 
   aryutt[269] + aryutt[268] + 21./2.*aryutt[18] + 33./2.*aryutt[19] + 
   3./4.*aryutt[266] + aryutt[267];
   aryutt[261]=aryutt[3]*aryutt[261];
   aryutt[263]=1./3.*aryutt[36];
   aryutt[266]=aryutt[263] + aryutt[148];
   aryutt[266]=aryutt[12]*aryutt[266];
   aryutt[273]=aryutt[72] - 119./12. + aryutt[74];
   aryutt[274]= - 1./3.*aryutt[79];
   aryutt[275]= - 1./3.*aryutt[36];
   aryutt[149]=aryutt[266] + aryutt[275] + 5./24.*aryutt[13] - 5./2.*
   aryutt[84] + aryutt[149] + aryutt[274] + 1./8.*aryutt[273] + 5./3.*
   aryutt[88];
   aryutt[149]=aryutt[99]*aryutt[149];
   aryutt[266]= - 1./2.*aryutt[64] - aryutt[60];
   aryutt[149]=1./3.*aryutt[266] + 1./2.*aryutt[149];
   aryutt[149]=5./4.*MMH*aryutt[149];
   aryutt[227]= - aryutt[69] + aryutt[227];
   aryutt[266]=1./4.*aryutt[67];
   aryutt[273]=11./4.*aryutt[64];
   aryutt[276]=11./2.*aryutt[60];
   aryutt[277]=1./2.*aryutt[63];
   aryutt[227]=aryutt[276] + aryutt[273] + aryutt[266] + aryutt[277] + 
   35./8.*aryutt[65] + 2*aryutt[227] + 1./16.*aryutt[66];
   aryutt[278]= - 1 - aryutt[37];
   aryutt[279]=1./2.*aryutt[21];
   aryutt[280]=5 + aryutt[264];
   aryutt[280]=aryutt[8]*aryutt[280];
   aryutt[278]=aryutt[280] + 3*aryutt[278] + aryutt[279];
   aryutt[278]=aryutt[3]*aryutt[278];
   aryutt[280]=1./4.*aryutt[278];
   aryutt[203]=3./4.*aryutt[203];
   aryutt[281]=3./4.*aryutt[112];
   aryutt[227]=aryutt[280] + aryutt[281] + 1./3.*aryutt[227] + 
   aryutt[203];
   aryutt[227]=MMt*aryutt[227];
   aryutt[282]=aryutt[185] + 29*aryutt[101];
   aryutt[282]=11./3.*aryutt[114] + 1./3.*aryutt[282] + aryutt[113];
   aryutt[282]=1./4.*aryutt[282];
   aryutt[283]=4133./24. + 297*aryutt[41];
   aryutt[284]= - 5*aryutt[36];
   aryutt[283]= - 225./32.*aryutt[21] + 1./4.*aryutt[283] + aryutt[284]
   ;
   aryutt[283]=1./4.*aryutt[283] + 2./3.*aryutt[12];
   aryutt[283]=aryutt[99]*aryutt[283];
   aryutt[283]=aryutt[282] + aryutt[283];
   aryutt[283]=aryutt[15]*aryutt[283];
   aryutt[285]= - 193 - 297./2.*aryutt[41];
   aryutt[229]=aryutt[229] + 1./2.*aryutt[285] + aryutt[228];
   aryutt[229]=aryutt[99]*aryutt[229];
   aryutt[285]= - 3./4.*aryutt[101];
   aryutt[229]=aryutt[285] + aryutt[229];
   aryutt[286]=1 + aryutt[140];
   aryutt[286]=aryutt[3]*aryutt[286];
   aryutt[229]=1./2.*aryutt[229] + aryutt[286];
   aryutt[229]=aryutt[17]*aryutt[229];
   aryutt[287]=aryutt[116] + 7./3.*aryutt[18] + aryutt[111] + 11./6.*
   aryutt[42];
   aryutt[287]=1./4.*aryutt[101]*aryutt[287];
   aryutt[288]=4913./24. - 113*aryutt[39];
   aryutt[265]=1./4.*aryutt[288] + aryutt[265];
   aryutt[265]=1./4.*aryutt[265] + 13./3.*aryutt[21];
   aryutt[265]=aryutt[21]*aryutt[265];
   aryutt[288]= - 337./32.*aryutt[80] + 1./4.*aryutt[43] - 113./256.*
   aryutt[76] + 43./32.*aryutt[83] - 31./12.*aryutt[58] - 13183./1024.
    + 17./3.*aryutt[81];
   aryutt[289]= - 157./192.*aryutt[51];
   aryutt[290]=7./24.*aryutt[234];
   aryutt[291]=aryutt[19]*aryutt[98];
   aryutt[292]=1./16.*aryutt[291];
   aryutt[293]= - aryutt[18]*aryutt[98];
   aryutt[294]=13./72.*aryutt[38];
   aryutt[295]= - 53./36.*aryutt[36];
   aryutt[296]=11./24.*aryutt[163];
   aryutt[297]= - 1 + aryutt[118];
   aryutt[297]=aryutt[8]*aryutt[297];
   aryutt[298]=1./2.*aryutt[297];
   aryutt[299]=397./64.*aryutt[21] + 11./2.*aryutt[104] + 13./3.*
   aryutt[36] + 7./6.*aryutt[38] + 29./4.*aryutt[41] + 35./3. - 31./2.*
   aryutt[39];
   aryutt[299]=79./48.*aryutt[12] + 1./3.*aryutt[299] + 9./16.*
   aryutt[8];
   aryutt[299]=aryutt[12]*aryutt[299];
   aryutt[300]=aryutt[16]*aryutt[98];
   aryutt[301]=3./8.*aryutt[300];
   aryutt[302]= - 27./2.*aryutt[50] - 7*aryutt[94];
   aryutt[302]= - 25*aryutt[92] - 39*aryutt[93] + 11*aryutt[302] + 225./
   4.*aryutt[96];
   aryutt[302]=1887./8.*aryutt[16] + 25./2.*aryutt[14] + 5279./48.*
   aryutt[20] + 10033./96.*aryutt[18] + 451./2.*aryutt[19] - 69./8.*
   aryutt[33] - 225./16.*aryutt[34] + 5./4.*aryutt[32] + 1./2.*
   aryutt[302] + 5*aryutt[48];
   aryutt[302]=aryutt[99]*aryutt[302];
   aryutt[303]= - 1./6.*aryutt[24];
   aryutt[304]= - 9./8.*aryutt[35];
   aryutt[305]= - 1./4.*aryutt[40];
   aryutt[103]=aryutt[103] + aryutt[227] + aryutt[106] + aryutt[149] + 
   aryutt[110] + aryutt[229] + aryutt[261] + aryutt[283] + 1./4.*
   aryutt[302] + aryutt[301] + 1./4.*aryutt[299] + aryutt[298] + 1./6.*
   aryutt[265] + aryutt[296] + aryutt[295] + aryutt[294] + aryutt[287]
    - 197./48.*aryutt[41] + aryutt[37] + 1./32.*aryutt[293] + 131./48.*
   aryutt[39] + 9055./2304.*aryutt[13] + aryutt[305] + aryutt[292] + 
   aryutt[290] + aryutt[303] + 7./12.*aryutt[25] - 17./32.*aryutt[23]
    + aryutt[84] + aryutt[289] + 7./9.*aryutt[45] - 11./24.*aryutt[79]
    - 145./288.*aryutt[88] + aryutt[304] + 301./384.*aryutt[87] - 19./
   192.*aryutt[72] - 19./192.*aryutt[74] + 11./48.*aryutt[75] + 11./12.
   *aryutt[77] - 79./128.*aryutt[90] + 1./3.*aryutt[288] + 1./8.*
   aryutt[141];
   aryutt[103]=aryutt[210]*aryutt[103];
   aryutt[106]=17041./6. + 289*aryutt[39];
   aryutt[110]= - 17*aryutt[41];
   aryutt[106]=1./3.*aryutt[106] + aryutt[110];
   aryutt[227]= - 17*aryutt[36];
   aryutt[106]= - 239./27.*aryutt[12] + 17./6.*aryutt[163] + 1./12.*
   aryutt[106] + aryutt[227];
   aryutt[106]=aryutt[12]*aryutt[106];
   aryutt[229]= - aryutt[3]*aryutt[12];
   aryutt[171]=17./4.*aryutt[229] + 17./2.*aryutt[176] - 257./9.*
   aryutt[66] + 17*aryutt[171];
   aryutt[171]=MMZ*aryutt[171];
   aryutt[117]=aryutt[117] + aryutt[156];
   aryutt[229]= - 1./3. + aryutt[131];
   aryutt[229]=aryutt[15]*aryutt[229];
   aryutt[117]=1./3.*aryutt[117] + aryutt[229];
   aryutt[117]=aryutt[3]*aryutt[117];
   aryutt[229]= - 637./9.*aryutt[75] - 7027./12. + 17*aryutt[80];
   aryutt[229]= - 17./2.*aryutt[84] + 721./9.*aryutt[45] + 17./2.*
   aryutt[79] - 221./48.*aryutt[88] + 1./8.*aryutt[229] + aryutt[87];
   aryutt[229]=1./4.*aryutt[229] - 4./9.*aryutt[13];
   aryutt[170]=1./3.*aryutt[170] + aryutt[18];
   aryutt[170]=aryutt[101]*aryutt[170];
   aryutt[106]=1./12.*aryutt[171] + 17./48.*aryutt[178] + 17./16.*
   aryutt[117] + 17./24.*aryutt[177] + 1./8.*aryutt[106] + 17./16.*
   aryutt[104] + 1./3.*aryutt[229] + 17./8.*aryutt[170];
   aryutt[106]=MMZ*aryutt[106];
   aryutt[110]=aryutt[110] + 2507./2. + 289./3.*aryutt[39];
   aryutt[110]=17./3.*aryutt[162] + 221./72.*aryutt[12] + 17./3.*
   aryutt[163] + 1./12.*aryutt[110] + aryutt[227];
   aryutt[110]=aryutt[15]*aryutt[110];
   aryutt[117]=aryutt[50] - aryutt[94];
   aryutt[117]=17./8.*aryutt[183] + 17./24.*aryutt[14] + 6065./144.*
   aryutt[20] + 977./12.*aryutt[18] + aryutt[19] + 221./96.*aryutt[42]
    + aryutt[154] - 2297./288.*aryutt[33] + 17./4.*aryutt[48] - 119./24.
   *aryutt[92] - 967./48.*aryutt[93] + 17./16.*aryutt[117] + aryutt[96]
   ;
   aryutt[162]= - 839./3. - 289./8.*aryutt[39];
   aryutt[162]=1./3.*aryutt[162] + 17./8.*aryutt[41];
   aryutt[170]=aryutt[3]*aryutt[15];
   aryutt[162]=17./8.*aryutt[170] + 1475./432.*aryutt[12] + 1./3.*
   aryutt[162] + 17./2.*aryutt[36];
   aryutt[162]=aryutt[17]*aryutt[162];
   aryutt[171]= - 1./2.*aryutt[39];
   aryutt[177]=1253./144.*aryutt[12] - 5./3. + aryutt[171];
   aryutt[177]=aryutt[16]*aryutt[177];
   aryutt[229]=pow(aryutt[15],2);
   aryutt[261]= - aryutt[3]*aryutt[229];
   aryutt[110]=1./2.*aryutt[162] + 17./16.*aryutt[261] + 1./4.*
   aryutt[110] + 1./3.*aryutt[117] + 1./2.*aryutt[177];
   aryutt[117]= - 17./3.*aryutt[79] + 187./6.*aryutt[88] + 17./2.*
   aryutt[72] - 41./9. + 17./2.*aryutt[74];
   aryutt[162]= - 17*aryutt[38];
   aryutt[177]= - 7./4. + aryutt[162];
   aryutt[177]=17./4.*aryutt[8] + 1./3.*aryutt[177] + 17*aryutt[36];
   aryutt[177]=aryutt[12]*aryutt[177];
   aryutt[265]= - 1./3.*aryutt[45];
   aryutt[117]=1./6.*aryutt[177] - 17./48.*aryutt[13] + 17./16.*
   aryutt[84] - 17./16.*aryutt[51] + 1./8.*aryutt[117] + aryutt[265];
   aryutt[117]=MMH*aryutt[117];
   aryutt[177]=3211./3. + 139./2.*aryutt[12];
   aryutt[177]=aryutt[15]*aryutt[177];
   aryutt[283]=1607 - 5557./2.*aryutt[12];
   aryutt[283]=aryutt[17]*aryutt[283];
   aryutt[177]=1./6.*aryutt[283] + 1./2.*aryutt[177] + 1927./3.*
   aryutt[20] + 535./2.*aryutt[18] - 107*aryutt[93] - 803./3.*
   aryutt[33];
   aryutt[283]= - 107*aryutt[13] - 1285./3.*aryutt[45] - 3211./3. + 107
   *aryutt[75];
   aryutt[288]= - 257./9. - 5./4.*aryutt[12];
   aryutt[288]=aryutt[12]*aryutt[288];
   aryutt[283]=1./3.*aryutt[283] + 5*aryutt[288];
   aryutt[283]=MMZ*aryutt[283];
   aryutt[177]=1./3.*aryutt[177] + 1./2.*aryutt[283];
   aryutt[177]=MMZ*aryutt[177];
   aryutt[283]= - 13265./3.*aryutt[15] + 1381./2.*aryutt[17];
   aryutt[283]=aryutt[17]*aryutt[283];
   aryutt[283]= - 931./2.*aryutt[229] + aryutt[283];
   aryutt[177]=1./12.*aryutt[283] + aryutt[177];
   aryutt[177]=aryutt[6]*aryutt[177];
   aryutt[106]=1./48.*aryutt[177] + 1./4.*aryutt[117] + 1./2.*
   aryutt[110] + aryutt[106];
   aryutt[106]=aryutt[6]*aryutt[106];
   aryutt[110]= - 17*aryutt[39];
   aryutt[117]=7./2. + aryutt[110];
   aryutt[117]=aryutt[41] + 1./2.*aryutt[117] + aryutt[264];
   aryutt[117]= - 7./48.*aryutt[12] + 1./4.*aryutt[117] + aryutt[201];
   aryutt[117]=aryutt[15]*aryutt[117];
   aryutt[177]= - aryutt[50] - 9./2.*aryutt[48];
   aryutt[283]=95./9. + aryutt[264];
   aryutt[283]=aryutt[14]*aryutt[283];
   aryutt[288]= - 65./6. + aryutt[39];
   aryutt[288]=aryutt[16]*aryutt[288];
   aryutt[117]=9./16.*aryutt[261] + aryutt[117] + 1./4.*aryutt[288] + 
   41./72.*aryutt[156] + 1./8.*aryutt[283] + 5./2.*aryutt[18] + 
   aryutt[187] + 1./4.*aryutt[177] + 5./9.*aryutt[42];
   aryutt[117]=aryutt[3]*aryutt[117];
   aryutt[177]= - 1./32.*aryutt[21];
   aryutt[148]=25./24.*aryutt[12] + aryutt[148] + aryutt[177] + 1./12.*
   aryutt[41] + aryutt[36];
   aryutt[148]=aryutt[12]*aryutt[148];
   aryutt[187]=993./4. + 1./3.*aryutt[76];
   aryutt[187]=1./4.*aryutt[87] + 1./6.*aryutt[72] + 1./6.*aryutt[74]
    - 25./9.*aryutt[75] + 1./8.*aryutt[187] - 1./9.*aryutt[80];
   aryutt[148]=1./3.*aryutt[148] + 1./96.*aryutt[21] + aryutt[275] - 1./
   36.*aryutt[41] + 289./288.*aryutt[13] - 1./24.*aryutt[51] + 
   aryutt[274] + 1./4.*aryutt[187] + aryutt[88];
   aryutt[148]=aryutt[99]*aryutt[148];
   aryutt[187]= - 11 - 17./2.*aryutt[39];
   aryutt[261]=aryutt[3]*aryutt[150];
   aryutt[172]=3./2.*aryutt[261] - 7./36.*aryutt[12] + aryutt[172] + 7./
   12.*aryutt[41] + 1./3.*aryutt[187] + aryutt[37];
   aryutt[172]=aryutt[3]*aryutt[172];
   aryutt[143]=7./18.*aryutt[13] - 7./9.*aryutt[84] + aryutt[144] - 7./
   18.*aryutt[79] + 7./18.*aryutt[88] - 37./9. + aryutt[143];
   aryutt[143]=aryutt[101]*aryutt[143];
   aryutt[144]=1./2. - aryutt[37];
   aryutt[144]=aryutt[153]*aryutt[144];
   aryutt[187]=MMZ*aryutt[144];
   aryutt[261]= - 1 - aryutt[87];
   aryutt[274]=aryutt[98]*aryutt[261];
   aryutt[274]=7./32.*aryutt[274] - 17./4.*aryutt[60] + 11*aryutt[66]
    - 17./8.*aryutt[64];
   aryutt[283]= - aryutt[8]*aryutt[101];
   aryutt[143]=3./4.*aryutt[187] + 1./4.*aryutt[172] + 25./16.*
   aryutt[148] + 7./72.*aryutt[176] + 3./8.*aryutt[283] + 1./9.*
   aryutt[274] + 1./4.*aryutt[143];
   aryutt[143]=MMZ*aryutt[143];
   aryutt[148]=aryutt[275] + aryutt[200];
   aryutt[148]=aryutt[12]*aryutt[148];
   aryutt[172]= - aryutt[72] + 119./12. - aryutt[74];
   aryutt[148]=aryutt[148] + aryutt[263] - 5./24.*aryutt[13] + 5./2.*
   aryutt[84] + aryutt[259] + aryutt[258] + 1./8.*aryutt[172] - 5./3.*
   aryutt[88];
   aryutt[148]=aryutt[99]*aryutt[148];
   aryutt[165]=aryutt[165] + aryutt[60];
   aryutt[148]=1./3.*aryutt[165] + 1./2.*aryutt[148];
   aryutt[148]=MMH*aryutt[148];
   aryutt[165]= - 5209./24. + aryutt[41];
   aryutt[165]=aryutt[177] + 1./12.*aryutt[165] + aryutt[36];
   aryutt[165]=aryutt[99]*aryutt[165];
   aryutt[172]=aryutt[185] + 47*aryutt[101];
   aryutt[113]=25./3.*aryutt[165] + 7./9.*aryutt[176] + 1./9.*
   aryutt[172] + aryutt[113];
   aryutt[113]=aryutt[15]*aryutt[113];
   aryutt[165]= - aryutt[68] - 437./54.*aryutt[66];
   aryutt[172]= - 1./6.*aryutt[67];
   aryutt[165]= - 7./3.*aryutt[60] - 7./6.*aryutt[64] + aryutt[172] + 1.
   /4.*aryutt[165] - 17./3.*aryutt[63];
   aryutt[112]=1./2.*aryutt[278] + 3./2.*aryutt[112] + 1./3.*
   aryutt[165] + aryutt[204];
   aryutt[112]=MMt*aryutt[112];
   aryutt[165]=aryutt[39] + aryutt[12];
   aryutt[165]=aryutt[15]*aryutt[165];
   aryutt[176]=aryutt[12]*aryutt[39];
   aryutt[185]=MMZ*aryutt[176];
   aryutt[187]= - aryutt[17]*aryutt[39];
   aryutt[200]= - aryutt[16]*aryutt[12];
   aryutt[165]=aryutt[185] + aryutt[187] + aryutt[200] + aryutt[165];
   aryutt[165]=aryutt[6]*aryutt[165];
   aryutt[176]=7./27.*aryutt[176] + aryutt[261] + 1./27.*aryutt[39];
   aryutt[185]= - 1./3. + aryutt[171];
   aryutt[185]=aryutt[15]*aryutt[185];
   aryutt[185]=aryutt[256] + aryutt[185];
   aryutt[185]=aryutt[3]*aryutt[185];
   aryutt[204]= - MMZ*aryutt[3]*aryutt[39];
   aryutt[165]=17./108.*aryutt[165] + 1./3.*aryutt[204] + 1./4.*
   aryutt[176] + aryutt[185];
   aryutt[176]=pow(aryutt[2],2);
   aryutt[165]=aryutt[176]*aryutt[165];
   aryutt[185]= - 1./2.*aryutt[50] + aryutt[94];
   aryutt[204]=1./4.*aryutt[96];
   aryutt[123]=aryutt[123] + 25./3.*aryutt[93] + 1./3.*aryutt[185] + 
   aryutt[204];
   aryutt[123]= - 11./24.*aryutt[16] + aryutt[124] - 1675./144.*
   aryutt[20] - 2869./288.*aryutt[18] - 1./6.*aryutt[19] + 193./72.*
   aryutt[33] - 1./16.*aryutt[34] + aryutt[125] + 1./2.*aryutt[123] - 
   aryutt[48];
   aryutt[123]=aryutt[99]*aryutt[123];
   aryutt[124]= - aryutt[90] - 41./18.*aryutt[80] + 331./144.*
   aryutt[76] - 56965./1728. + aryutt[83];
   aryutt[185]= - 9./2.*aryutt[35];
   aryutt[124]=7./18.*aryutt[79] - 595./216.*aryutt[88] + aryutt[185]
    - 967./288.*aryutt[87] + 143./144.*aryutt[72] + 143./144.*
   aryutt[74] - 191./324.*aryutt[75] + 1./4.*aryutt[124] - 7./9.*
   aryutt[77];
   aryutt[124]= - 31./576.*aryutt[51] + 1./4.*aryutt[124] - 139./81.*
   aryutt[45];
   aryutt[111]=aryutt[116] + 61./9.*aryutt[18] + aryutt[111] - 7./18.*
   aryutt[42];
   aryutt[111]=aryutt[101]*aryutt[111];
   aryutt[258]= - 3329./2. + 119*aryutt[39];
   aryutt[258]=1./3.*aryutt[258] - 7./2.*aryutt[41];
   aryutt[258]= - 227./96.*aryutt[21] + 7./3.*aryutt[163] + 1./6.*
   aryutt[258] - 7*aryutt[36];
   aryutt[258]= - 235./648.*aryutt[12] + 1./3.*aryutt[258] - 31./8.*
   aryutt[8];
   aryutt[258]=aryutt[12]*aryutt[258];
   aryutt[261]=67./3. - aryutt[41];
   aryutt[261]=1./32.*aryutt[21] + 1./12.*aryutt[261] - aryutt[36];
   aryutt[261]=aryutt[99]*aryutt[261];
   aryutt[261]=aryutt[285] + 25./3.*aryutt[261];
   aryutt[261]=1./2.*aryutt[261] + aryutt[286];
   aryutt[261]=aryutt[17]*aryutt[261];
   aryutt[274]=25./16. - aryutt[39];
   aryutt[274]=aryutt[21]*aryutt[274];
   aryutt[278]=1./6.*aryutt[38];
   aryutt[106]=1./8.*aryutt[165] + 1./4.*aryutt[112] + 1./3.*
   aryutt[106] + 25./24.*aryutt[148] + aryutt[143] + 1./2.*aryutt[261]
    + 1./2.*aryutt[117] + 1./8.*aryutt[113] + 25./24.*aryutt[123] + 1./
   16.*aryutt[300] + 1./16.*aryutt[258] + 1./4.*aryutt[297] + 1./96.*
   aryutt[274] + 7./144.*aryutt[104] - 79./144.*aryutt[36] + 
   aryutt[278] + 1./8.*aryutt[111] + 5./1728.*aryutt[41] + aryutt[199]
    + 35./864.*aryutt[39] - 2983./41472.*aryutt[13] - 1./8.*aryutt[40]
    + 1./96.*aryutt[291] + 7./144.*aryutt[234] - 1./36.*aryutt[24] - 1./
   32.*aryutt[23] + 1./2.*aryutt[124] - 1./3.*aryutt[84];
   aryutt[106]=aryutt[176]*aryutt[106];
   aryutt[111]=14./27.*aryutt[57] + aryutt[54];
   aryutt[112]=1./8.*aryutt[155] - 13./16.*aryutt[87] - 1 - 11./48.*
   aryutt[90];
   aryutt[112]=aryutt[98]*aryutt[112];
   aryutt[111]=1./2.*aryutt[112] + aryutt[226] + aryutt[225] + 
   aryutt[172] - 7./2.*aryutt[66] + 2*aryutt[111] + 1./4.*aryutt[68];
   aryutt[112]=aryutt[16]*aryutt[146];
   aryutt[112]=aryutt[112] + aryutt[150];
   aryutt[112]=aryutt[3]*aryutt[112];
   aryutt[113]= - aryutt[38] + 31./6.*aryutt[41] + aryutt[37] + 19./4.*
   aryutt[39] - 59./18. + aryutt[40];
   aryutt[112]=3./4.*aryutt[112] + 11./24.*aryutt[12] + aryutt[222] + 1.
   /2.*aryutt[113] + aryutt[201];
   aryutt[112]=aryutt[3]*aryutt[112];
   aryutt[113]= - 224603./36. + 7*aryutt[76];
   aryutt[113]=275./3.*aryutt[75] + 5./4.*aryutt[90] + 1./8.*
   aryutt[113] - 37./3.*aryutt[80];
   aryutt[113]=7./4.*aryutt[87] + aryutt[257] + 1./3.*aryutt[113] + 
   aryutt[254];
   aryutt[117]=7./96.*aryutt[21];
   aryutt[113]=aryutt[117] + aryutt[36] - 37./36.*aryutt[41] - 9091./
   864.*aryutt[13] + aryutt[259] + aryutt[79] + 1./4.*aryutt[113] - 3*
   aryutt[88];
   aryutt[123]=10./3. + 37./32.*aryutt[41];
   aryutt[123]= - 275./576.*aryutt[12] + 1./64.*aryutt[8] - 7./768.*
   aryutt[21] + 1./9.*aryutt[123] - 1./8.*aryutt[36];
   aryutt[123]=aryutt[12]*aryutt[123];
   aryutt[113]=1./8.*aryutt[113] + aryutt[123];
   aryutt[113]=aryutt[99]*aryutt[113];
   aryutt[123]=1 + aryutt[87];
   aryutt[123]=aryutt[151]*aryutt[123];
   aryutt[123]=1./48.*aryutt[123] + 3*aryutt[144];
   aryutt[123]=MMZ*aryutt[123];
   aryutt[111]=1./2.*aryutt[123] + aryutt[112] + 1./48.*aryutt[160] + 5
   *aryutt[113] + 1./32.*aryutt[159] + aryutt[236] + 3./4.*aryutt[283]
    + 1./3.*aryutt[111] + aryutt[224];
   aryutt[111]=MMZ*aryutt[111];
   aryutt[112]=46471./216. + aryutt[76];
   aryutt[112]=143./18.*aryutt[80] + 1./2.*aryutt[112] + 1./3.*
   aryutt[43];
   aryutt[113]=1./3.*aryutt[90];
   aryutt[123]=1./4.*aryutt[84];
   aryutt[112]=607./432.*aryutt[13] + 1./8.*aryutt[24] - 68./27.*
   aryutt[25] - 1./12.*aryutt[23] + aryutt[123] - 457./54.*aryutt[45]
    + aryutt[193] + 13./96.*aryutt[88] + 1./12.*aryutt[87] - 377./432.*
   aryutt[75] + 1./4.*aryutt[112] + aryutt[113];
   aryutt[124]=2*aryutt[57] - aryutt[54];
   aryutt[124]=34*aryutt[124] + 263./2.*aryutt[66];
   aryutt[124]=1./9.*aryutt[124] - 1./4.*aryutt[67];
   aryutt[143]=aryutt[3]*aryutt[12];
   aryutt[114]=1./8.*aryutt[143] + 1./4.*aryutt[114] + 1./3.*
   aryutt[124] + 1./2.*aryutt[174];
   aryutt[114]=MMZ*aryutt[114];
   aryutt[124]=aryutt[205] + aryutt[183];
   aryutt[143]= - 1 - 17./6.*aryutt[12];
   aryutt[143]=aryutt[16]*aryutt[143];
   aryutt[144]=1./3. + aryutt[194];
   aryutt[144]=aryutt[15]*aryutt[144];
   aryutt[124]=aryutt[144] + 1./3.*aryutt[124] + aryutt[143];
   aryutt[124]=aryutt[3]*aryutt[124];
   aryutt[143]=17*aryutt[38];
   aryutt[144]=aryutt[284] + aryutt[143] - 821./12.*aryutt[41] - 3659./
   6. - 73*aryutt[39];
   aryutt[144]=1./3.*aryutt[144] + aryutt[104];
   aryutt[144]=1./3.*aryutt[144] - 3./8.*aryutt[21];
   aryutt[144]=1./8.*aryutt[144] + 59./81.*aryutt[12];
   aryutt[144]=aryutt[12]*aryutt[144];
   aryutt[146]= - aryutt[21]*aryutt[39];
   aryutt[112]=1./3.*aryutt[114] + 5./8.*aryutt[178] + 1./8.*
   aryutt[124] + aryutt[249] + aryutt[144] + 1./24.*aryutt[146] + 1./8.
   *aryutt[163] + 1./3.*aryutt[112] + aryutt[243];
   aryutt[112]=MMZ*aryutt[112];
   aryutt[114]= - 5./4.*aryutt[36] + 17./4.*aryutt[38] - 821./48.*
   aryutt[41] - 917./3. - 73./4.*aryutt[39];
   aryutt[105]=1./2.*aryutt[235] - 4013./432.*aryutt[12] - 35./96.*
   aryutt[21] + 1./3.*aryutt[114] + aryutt[105];
   aryutt[105]=aryutt[15]*aryutt[105];
   aryutt[114]=277./12.*aryutt[12] + 19./8.*aryutt[21] + aryutt[228] + 
   aryutt[162] + 821./12.*aryutt[41] + 25201./36. + 73*aryutt[39];
   aryutt[124]=17./3.*aryutt[16] + 15*aryutt[15];
   aryutt[124]=aryutt[3]*aryutt[124];
   aryutt[114]=1./9.*aryutt[114] + 1./2.*aryutt[124];
   aryutt[114]=aryutt[17]*aryutt[114];
   aryutt[124]=295*aryutt[50] - 337*aryutt[94];
   aryutt[124]=1993./27.*aryutt[93] + 1./27.*aryutt[124] + 3*aryutt[96]
   ;
   aryutt[124]=1./2.*aryutt[124] + 7./9.*aryutt[92];
   aryutt[124]=1./2.*aryutt[124] - 1./3.*aryutt[48];
   aryutt[144]=43./6. - 37*aryutt[39];
   aryutt[144]= - 151./18.*aryutt[12] + 1./2.*aryutt[144] + aryutt[21];
   aryutt[144]=aryutt[16]*aryutt[144];
   aryutt[148]= - 13./3.*aryutt[16] + aryutt[188];
   aryutt[148]=aryutt[3]*aryutt[15]*aryutt[148];
   aryutt[105]=1./4.*aryutt[114] + 1./4.*aryutt[148] + 1./3.*
   aryutt[105] + 1./12.*aryutt[144] + aryutt[253] + aryutt[252] - 12533.
   /648.*aryutt[20] - 5899./216.*aryutt[18] + 23./48.*aryutt[19] + 
   aryutt[251] - 13./36.*aryutt[5] + 547./144.*aryutt[33] + 1./2.*
   aryutt[124] - 1./9.*aryutt[34];
   aryutt[114]=173*aryutt[13] + 1043./3.*aryutt[45] + 4157./3. - 173*
   aryutt[75];
   aryutt[124]=1043./27. + 5./4.*aryutt[12];
   aryutt[124]=aryutt[12]*aryutt[124];
   aryutt[114]=1./9.*aryutt[114] + aryutt[124];
   aryutt[114]=MMZ*aryutt[114];
   aryutt[124]=173*aryutt[93] + 781./3.*aryutt[33];
   aryutt[124]= - 875./9.*aryutt[20] - 865./18.*aryutt[18] + 1./9.*
   aryutt[124] + 11./2.*aryutt[19];
   aryutt[144]= - 11./3. + 9./2.*aryutt[12];
   aryutt[144]=aryutt[16]*aryutt[144];
   aryutt[148]= - 4157./3. - 301./2.*aryutt[12];
   aryutt[148]=aryutt[15]*aryutt[148];
   aryutt[150]= - 529 + 4259./2.*aryutt[12];
   aryutt[150]=aryutt[17]*aryutt[150];
   aryutt[114]=1./6.*aryutt[114] + 1./162.*aryutt[150] + 1./54.*
   aryutt[148] + 1./3.*aryutt[124] + 1./2.*aryutt[144];
   aryutt[114]=MMZ*aryutt[114];
   aryutt[124]=5*aryutt[16];
   aryutt[144]=aryutt[124] + 1429./18.*aryutt[15];
   aryutt[144]=aryutt[15]*aryutt[144];
   aryutt[148]=925./54.*aryutt[17] - 9*aryutt[16] + 11495./81.*
   aryutt[15];
   aryutt[148]=aryutt[17]*aryutt[148];
   aryutt[144]=1./3.*aryutt[144] + aryutt[148];
   aryutt[114]=1./4.*aryutt[144] + aryutt[114];
   aryutt[144]=aryutt[6]*MMZ*aryutt[255];
   aryutt[114]=1./8.*aryutt[114] + 17./27.*aryutt[144];
   aryutt[114]=aryutt[6]*aryutt[114];
   aryutt[105]=aryutt[114] + aryutt[169] + 1./2.*aryutt[105] + 
   aryutt[112];
   aryutt[105]=aryutt[6]*aryutt[105];
   aryutt[112]=aryutt[260] + aryutt[241] + aryutt[128] + aryutt[138] + 
   35./8.*aryutt[41] + aryutt[122] + aryutt[119] + 43./3. + aryutt[139]
   ;
   aryutt[112]=aryutt[15]*aryutt[112];
   aryutt[114]=197./6. + aryutt[262];
   aryutt[114]= - 7./12.*aryutt[12] + 1./2.*aryutt[114] + aryutt[264];
   aryutt[114]=aryutt[16]*aryutt[114];
   aryutt[119]= - aryutt[50] + aryutt[219];
   aryutt[112]=aryutt[272] + 1./2.*aryutt[112] + 1./4.*aryutt[114] + 
   aryutt[269] + aryutt[268] + 9./2.*aryutt[18] + aryutt[126] + 9./4.*
   aryutt[119] + aryutt[267];
   aryutt[112]=aryutt[3]*aryutt[112];
   aryutt[114]=28*aryutt[54] - 997./16.*aryutt[66];
   aryutt[114]=aryutt[276] + aryutt[273] + aryutt[266] + aryutt[277] + 
   1./27.*aryutt[114] + 19./8.*aryutt[65];
   aryutt[114]=aryutt[280] + aryutt[281] + 1./3.*aryutt[114] + 
   aryutt[203];
   aryutt[114]=MMt*aryutt[114];
   aryutt[119]=271./8.*aryutt[76] - 253733./864. + 7*aryutt[83];
   aryutt[119]= - 473./72.*aryutt[80] + 1./8.*aryutt[119] + aryutt[43];
   aryutt[119]= - 11./6.*aryutt[79] - 145./72.*aryutt[88] + aryutt[185]
    + 439./288.*aryutt[87] - 19./48.*aryutt[72] - 19./48.*aryutt[74] + 
   809./324.*aryutt[75] + 11./3.*aryutt[77] - 109./96.*aryutt[90] + 1./
   3.*aryutt[119] + 1./2.*aryutt[141];
   aryutt[126]=5153./8. + 37./3.*aryutt[41];
   aryutt[126]= - 7./96.*aryutt[21] + 1./12.*aryutt[126] - aryutt[36];
   aryutt[126]=1./4.*aryutt[126] + 10./27.*aryutt[12];
   aryutt[126]=aryutt[99]*aryutt[126];
   aryutt[126]=aryutt[282] + 5*aryutt[126];
   aryutt[126]=aryutt[15]*aryutt[126];
   aryutt[128]=13./2.*aryutt[36] + 7./4.*aryutt[38] - 71./24.*
   aryutt[41] - 173./9. - 29./4.*aryutt[39];
   aryutt[104]=269./128.*aryutt[21] + 1./3.*aryutt[128] + 11./4.*
   aryutt[104];
   aryutt[104]=85./2592.*aryutt[12] + 1./3.*aryutt[104] + 9./32.*
   aryutt[8];
   aryutt[104]=aryutt[12]*aryutt[104];
   aryutt[128]= - 37./2.*aryutt[50] + 47*aryutt[94];
   aryutt[128]= - 275./3.*aryutt[93] + 1./3.*aryutt[128] + 7./4.*
   aryutt[96];
   aryutt[128]=1./3.*aryutt[128] - 5*aryutt[92];
   aryutt[138]=1./4.*aryutt[32];
   aryutt[128]= - 437./72.*aryutt[16] + 5./2.*aryutt[14] + 4985./144.*
   aryutt[20] + 29797./864.*aryutt[18] - 19./6.*aryutt[19] - 1499./216.
   *aryutt[33] - 7./48.*aryutt[34] + aryutt[138] + 1./2.*aryutt[128] + 
   aryutt[48];
   aryutt[128]=aryutt[99]*aryutt[128];
   aryutt[144]= - 643./3. - 37./2.*aryutt[41];
   aryutt[117]=aryutt[117] + 1./18.*aryutt[144] + aryutt[36];
   aryutt[117]=aryutt[99]*aryutt[117];
   aryutt[117]=aryutt[285] + 5*aryutt[117];
   aryutt[117]=1./2.*aryutt[117] + aryutt[286];
   aryutt[117]=aryutt[17]*aryutt[117];
   aryutt[144]=49./24. + aryutt[110];
   aryutt[144]=1./8.*aryutt[144] - aryutt[21];
   aryutt[144]=aryutt[21]*aryutt[144];
   aryutt[148]= - 1./18.*aryutt[24];
   aryutt[104]=aryutt[106] + aryutt[114] + aryutt[105] + aryutt[149] + 
   aryutt[111] + aryutt[117] + aryutt[112] + aryutt[126] + 5./4.*
   aryutt[128] + aryutt[301] + 1./2.*aryutt[104] + aryutt[298] + 1./12.
   *aryutt[144] + aryutt[296] + aryutt[295] + aryutt[294] + aryutt[287]
    + 1./144.*aryutt[41] + aryutt[37] + 1./96.*aryutt[136] + 41./144.*
   aryutt[39] + 2461./6912.*aryutt[13] + aryutt[305] + aryutt[292] + 
   aryutt[290] + aryutt[148] - 103./324.*aryutt[25] - 5./32.*aryutt[23]
    + aryutt[84] + aryutt[289] + 1./4.*aryutt[119] - 1./9.*aryutt[45];
   aryutt[104]=aryutt[176]*aryutt[104];
   aryutt[105]=3*aryutt[97];
   aryutt[106]=3./4.*aryutt[71];
   aryutt[111]= - 3*aryutt[86];
   aryutt[112]=1./2.*aryutt[77];
   aryutt[114]=1./2.*aryutt[85];
   aryutt[117]=27./4.*aryutt[35];
   aryutt[119]=1./4.*aryutt[79];
   aryutt[126]= - 21./4.*aryutt[51];
   aryutt[128]= - 1./4.*aryutt[13];
   aryutt[136]= - 1./4.*aryutt[84];
   aryutt[144]= - 1./6.*aryutt[38];
   aryutt[149]= - 1./12.*aryutt[36] + aryutt[144] + aryutt[140] + 
   aryutt[128] + aryutt[139] + aryutt[135] + aryutt[136] + aryutt[126]
    + aryutt[119] + aryutt[117] + aryutt[114] + aryutt[112] + 
   aryutt[134] + aryutt[137] + aryutt[78] + aryutt[111] + aryutt[106]
    - 293./32. + aryutt[105];
   aryutt[150]= - 11./3. - 9*aryutt[35];
   aryutt[150]=aryutt[129] + 1./2.*aryutt[150] - aryutt[40];
   aryutt[155]=1./2.*aryutt[150];
   aryutt[159]=1./4.*aryutt[8];
   aryutt[160]=aryutt[159] + aryutt[155] + aryutt[275];
   aryutt[160]=aryutt[8]*aryutt[160];
   aryutt[162]=1./2.*aryutt[18];
   aryutt[163]= - aryutt[14] - aryutt[20] + aryutt[162] + aryutt[92] + 
   aryutt[219];
   aryutt[163]=aryutt[99]*aryutt[163];
   aryutt[165]=aryutt[38] + aryutt[107];
   aryutt[165]=aryutt[12]*aryutt[165];
   aryutt[169]= - 1 + aryutt[107];
   aryutt[169]=aryutt[15]*aryutt[99]*aryutt[169];
   aryutt[174]= - 1./2. - aryutt[36];
   aryutt[178]=aryutt[17]*aryutt[99]*aryutt[174];
   aryutt[185]=1 + aryutt[84];
   aryutt[193]=MMH*aryutt[99]*aryutt[185];
   aryutt[149]=1./4.*aryutt[193] + 1./2.*aryutt[178] + 1./2.*
   aryutt[169] + 1./2.*aryutt[163] + 1./12.*aryutt[165] + 1./2.*
   aryutt[149] + aryutt[160];
   aryutt[149]=MMH*aryutt[149];
   aryutt[160]= - 1./4.*aryutt[47] - aryutt[91];
   aryutt[165]= - 1./2.*aryutt[32];
   aryutt[160]=aryutt[165] - 3./8.*aryutt[48] + 1./4.*aryutt[92] - 1./4.
   *aryutt[93] - 1./2.*aryutt[94] - aryutt[49] + 3./2.*aryutt[160] + 
   aryutt[95];
   aryutt[203]=25./8.*aryutt[19];
   aryutt[205]=aryutt[116] + 11./8.*aryutt[18] + aryutt[203] + 
   aryutt[160] + 27./8.*aryutt[42];
   aryutt[219]= - 1./8.*aryutt[21];
   aryutt[222]=1./6.*aryutt[36];
   aryutt[224]= - 1./24.*aryutt[12];
   aryutt[225]=aryutt[224] + 7./16.*aryutt[8] + aryutt[219] + 
   aryutt[222] + 1./16. + aryutt[38];
   aryutt[225]=aryutt[16]*aryutt[225];
   aryutt[226]=aryutt[199] + aryutt[40] + 55./3. + 9./2.*aryutt[35];
   aryutt[226]=1./8.*aryutt[226] + aryutt[263];
   aryutt[226]=aryutt[14]*aryutt[226];
   aryutt[228]= - aryutt[8]*aryutt[14];
   aryutt[234]=1./4.*aryutt[228];
   aryutt[235]=aryutt[36] + 5./6. + aryutt[38];
   aryutt[235]= - 1./6.*aryutt[12] + 1./2.*aryutt[235] - 5./3.*
   aryutt[8];
   aryutt[235]=aryutt[15]*aryutt[235];
   aryutt[182]=1./2.*aryutt[235] + aryutt[225] + 1./8.*aryutt[183] + 
   aryutt[234] + aryutt[182] + 1./2.*aryutt[205] + aryutt[226];
   aryutt[205]=1./16.*aryutt[37] + 1./8.*aryutt[40] + 1./3. + 9./16.*
   aryutt[35];
   aryutt[225]= - 7./16.*aryutt[8];
   aryutt[226]=aryutt[225] + aryutt[205] + aryutt[222];
   aryutt[226]=aryutt[17]*aryutt[226];
   aryutt[149]=1./4.*aryutt[149] + 1./2.*aryutt[182] + aryutt[226];
   aryutt[149]=MMH*aryutt[149];
   aryutt[182]= - aryutt[16]*aryutt[38];
   aryutt[180]=aryutt[180] - aryutt[36];
   aryutt[226]=aryutt[15]*aryutt[180];
   aryutt[235]=1./2.*aryutt[38];
   aryutt[236]=aryutt[235] + aryutt[36];
   aryutt[236]=aryutt[17]*aryutt[236];
   aryutt[226]=1./3.*aryutt[236] + aryutt[182] + 1./3.*aryutt[226];
   aryutt[226]=MMH*aryutt[226];
   aryutt[236]= - aryutt[14] + aryutt[256];
   aryutt[241]=1./3.*aryutt[15];
   aryutt[236]=1./2.*aryutt[236] + aryutt[241];
   aryutt[236]=aryutt[15]*aryutt[236];
   aryutt[243]= - 1./3.*aryutt[16];
   aryutt[249]=aryutt[14] + aryutt[243];
   aryutt[251]= - 1./3.*aryutt[15];
   aryutt[249]=1./2.*aryutt[249] + aryutt[251];
   aryutt[249]=aryutt[17]*aryutt[249];
   aryutt[252]= - aryutt[14] + aryutt[16];
   aryutt[253]=aryutt[16]*aryutt[252];
   aryutt[226]=aryutt[226] + aryutt[249] + aryutt[253] + aryutt[236];
   aryutt[226]=aryutt[6]*MMH*aryutt[226];
   aryutt[236]=pow(aryutt[14],2);
   aryutt[249]= - 3./2.*aryutt[236];
   aryutt[254]=11*aryutt[14] - 137./2.*aryutt[16];
   aryutt[254]=aryutt[16]*aryutt[254];
   aryutt[257]=107*aryutt[14] + 151*aryutt[16];
   aryutt[257]=1./3.*aryutt[257] - 11./2.*aryutt[15];
   aryutt[257]=aryutt[15]*aryutt[257];
   aryutt[254]=aryutt[257] + aryutt[249] + aryutt[254];
   aryutt[257]=11./3.*aryutt[15] - aryutt[14] - 17./2.*aryutt[16];
   aryutt[257]=1./4.*aryutt[257] + aryutt[17];
   aryutt[257]=aryutt[17]*aryutt[257];
   aryutt[254]=1./16.*aryutt[254] + aryutt[257];
   aryutt[149]=1./8.*aryutt[226] + 1./2.*aryutt[254] + aryutt[149];
   aryutt[226]=101./64. - aryutt[71];
   aryutt[226]=3./64.*aryutt[141] - 5./16.*aryutt[73] + 1./2.*
   aryutt[226] + aryutt[86];
   aryutt[254]=15 + 11./4.*aryutt[21];
   aryutt[257]=3*aryutt[8];
   aryutt[254]=1./2.*aryutt[254] + aryutt[257];
   aryutt[254]=aryutt[8]*aryutt[254];
   aryutt[258]=aryutt[181] + 25*aryutt[42] - 9*aryutt[14];
   aryutt[259]=aryutt[8]*aryutt[14];
   aryutt[258]=1./8.*aryutt[258] + aryutt[259];
   aryutt[258]=aryutt[3]*aryutt[258];
   aryutt[260]=3./2.*aryutt[21];
   aryutt[261]=1 + aryutt[260];
   aryutt[261]=aryutt[21]*aryutt[261];
   aryutt[239]=aryutt[239] + 3*aryutt[59] - aryutt[62];
   aryutt[239]=MMH*aryutt[239];
   aryutt[262]=MMt*aryutt[62];
   aryutt[266]=1./2.*aryutt[262] + 1./4.*aryutt[239] + 1./2.*
   aryutt[258] + 1./4.*aryutt[254] + 1./64.*aryutt[261] + 13./32.*
   aryutt[51] + 1./16.*aryutt[85] + 1./2.*aryutt[226] + aryutt[44];
   aryutt[266]=MMt*aryutt[266];
   aryutt[267]= - 5*aryutt[91] - aryutt[95];
   aryutt[268]=1./4.*aryutt[93];
   aryutt[192]=7./4.*aryutt[42] + aryutt[192] - 3*aryutt[32] - 3./2.*
   aryutt[92] + aryutt[268] + aryutt[204] + 3./4.*aryutt[267] - 
   aryutt[94];
   aryutt[204]=25./2.*aryutt[20];
   aryutt[267]=35./2.*aryutt[14];
   aryutt[269]=aryutt[238] + aryutt[267] + aryutt[204] - 7*aryutt[18]
    + aryutt[192] - 21./2.*aryutt[19];
   aryutt[272]=11 - 5*aryutt[21];
   aryutt[273]= - 7*aryutt[8];
   aryutt[272]=9./4.*aryutt[272] + aryutt[273];
   aryutt[272]=aryutt[16]*aryutt[272];
   aryutt[269]=1./2.*aryutt[272] + 1./2.*aryutt[269] + aryutt[259];
   aryutt[272]=47./2. - 31./3.*aryutt[21];
   aryutt[272]=1./16.*aryutt[272] + 5./3.*aryutt[8];
   aryutt[272]=aryutt[15]*aryutt[272];
   aryutt[269]=1./4.*aryutt[269] + aryutt[272];
   aryutt[272]= - 31./2. - 9*aryutt[97];
   aryutt[272]=aryutt[72] + aryutt[74] - aryutt[75] - 5./2.*aryutt[77]
    - 3./2.*aryutt[82] + aryutt[137] + aryutt[147] + 5./4.*aryutt[73]
    + 1./2.*aryutt[272] - 5*aryutt[86];
   aryutt[136]= - 1./8.*aryutt[37] + aryutt[305] + 3./8.*aryutt[23] + 
   aryutt[136] + 31./16.*aryutt[51] + aryutt[304] - 1./4.*aryutt[85] + 
   1./4.*aryutt[272] - aryutt[44];
   aryutt[272]= - 1./4. + aryutt[36];
   aryutt[272]=aryutt[21]*aryutt[272];
   aryutt[272]=aryutt[136] + 1./6.*aryutt[272];
   aryutt[274]=9*aryutt[35];
   aryutt[276]=11./3. + aryutt[274];
   aryutt[276]=aryutt[199] + 1./2.*aryutt[276] + aryutt[40];
   aryutt[280]=1./8.*aryutt[276];
   aryutt[281]=aryutt[244] + aryutt[280] + aryutt[263];
   aryutt[281]=aryutt[8]*aryutt[281];
   aryutt[282]=aryutt[8] + aryutt[131];
   aryutt[282]=aryutt[12]*aryutt[282];
   aryutt[284]=1./8.*aryutt[282];
   aryutt[285]=aryutt[64] - 9*aryutt[59] - aryutt[62];
   aryutt[285]=1./2.*aryutt[285] + aryutt[60];
   aryutt[285]=MMH*aryutt[285];
   aryutt[286]=1./8.*aryutt[285];
   aryutt[272]=aryutt[286] + aryutt[284] + 1./2.*aryutt[272] + 
   aryutt[281];
   aryutt[272]=MMH*aryutt[272];
   aryutt[257]=aryutt[257] + 7 + aryutt[260];
   aryutt[260]= - aryutt[3]*aryutt[14];
   aryutt[281]=1./8.*aryutt[257] + aryutt[260];
   aryutt[281]=aryutt[17]*aryutt[281];
   aryutt[269]=aryutt[266] + 1./2.*aryutt[272] + 1./2.*aryutt[269] + 
   aryutt[281];
   aryutt[269]=MMt*aryutt[269];
   aryutt[272]=1./3.*aryutt[38];
   aryutt[224]=aryutt[224] - 17./24.*aryutt[8] - 1./24.*aryutt[21] + 
   aryutt[275] + 1./8. + aryutt[272];
   aryutt[224]=aryutt[16]*aryutt[224];
   aryutt[287]=1./3.*aryutt[216];
   aryutt[288]= - aryutt[38] + aryutt[36];
   aryutt[289]=aryutt[288] + aryutt[287];
   aryutt[290]=aryutt[8]*aryutt[214];
   aryutt[291]=1./3.*aryutt[290];
   aryutt[215]=1./12.*aryutt[215] + 1./4.*aryutt[289] + aryutt[291];
   aryutt[215]=MMH*aryutt[215];
   aryutt[289]=aryutt[14]*aryutt[288];
   aryutt[292]= - 1./2. + aryutt[38];
   aryutt[293]=1./6.*aryutt[12] + 17./6.*aryutt[8] + 1./6.*aryutt[21]
    + aryutt[292] - aryutt[36];
   aryutt[293]=aryutt[15]*aryutt[293];
   aryutt[294]=aryutt[17]*aryutt[288];
   aryutt[215]=1./2.*aryutt[215] + 1./3.*aryutt[294] + 1./4.*
   aryutt[293] + 1./3.*aryutt[289] + aryutt[224];
   aryutt[215]=MMH*aryutt[215];
   aryutt[224]=17*aryutt[14] + 5./2.*aryutt[16];
   aryutt[224]=aryutt[16]*aryutt[224];
   aryutt[289]= - 17*aryutt[14];
   aryutt[293]=aryutt[289] - 11./2.*aryutt[16];
   aryutt[293]=1./3.*aryutt[293] + aryutt[15];
   aryutt[293]=aryutt[15]*aryutt[293];
   aryutt[295]=aryutt[16] - aryutt[15];
   aryutt[296]=aryutt[17]*aryutt[295];
   aryutt[224]=17./3.*aryutt[296] + 1./3.*aryutt[224] + aryutt[293];
   aryutt[293]=aryutt[16]*aryutt[288];
   aryutt[297]=aryutt[15]*aryutt[288];
   aryutt[298]=aryutt[17]*aryutt[214];
   aryutt[293]=1./2.*aryutt[298] + aryutt[293] + 1./2.*aryutt[297];
   aryutt[293]=MMH*aryutt[293];
   aryutt[297]= - aryutt[16] - aryutt[15];
   aryutt[297]=aryutt[15]*aryutt[297];
   aryutt[297]=aryutt[186] + 1./2.*aryutt[297];
   aryutt[299]= - aryutt[16] + aryutt[15];
   aryutt[300]=aryutt[17]*aryutt[299];
   aryutt[301]=1./2.*aryutt[300];
   aryutt[293]=aryutt[293] + aryutt[297] + aryutt[301];
   aryutt[293]=aryutt[6]*MMH*aryutt[293];
   aryutt[215]=1./12.*aryutt[293] + 1./4.*aryutt[224] + aryutt[215];
   aryutt[224]=aryutt[202] + aryutt[8];
   aryutt[293]=aryutt[16]*aryutt[224];
   aryutt[302]=aryutt[130] - aryutt[8];
   aryutt[304]=aryutt[15]*aryutt[302];
   aryutt[293]=aryutt[293] + aryutt[304];
   aryutt[304]=aryutt[21]*aryutt[288];
   aryutt[305]=aryutt[8]*aryutt[288];
   aryutt[306]=1./4.*aryutt[304] + aryutt[305];
   aryutt[306]=MMH*aryutt[306];
   aryutt[293]=17./4.*aryutt[293] + aryutt[306];
   aryutt[293]=MMt*aryutt[293];
   aryutt[215]=1./2.*aryutt[215] + 1./3.*aryutt[293];
   aryutt[215]=aryutt[210]*aryutt[215];
   aryutt[149]=1./2.*aryutt[215] + 1./2.*aryutt[149] + aryutt[269];
   aryutt[149]=aryutt[210]*aryutt[149];
   aryutt[215]= - 1./3.*aryutt[38];
   aryutt[150]=1./2.*aryutt[8] + aryutt[275] + aryutt[150] + 
   aryutt[215];
   aryutt[150]=aryutt[8]*aryutt[150];
   aryutt[269]= - aryutt[21]*aryutt[38];
   aryutt[150]=aryutt[150] + 1./12.*aryutt[269] - 17./54.*aryutt[36] + 
   5./54.*aryutt[38] + aryutt[140] - 11./12.*aryutt[13] + aryutt[139]
    + aryutt[135] + 13./12.*aryutt[84] + aryutt[126] + 11./12.*
   aryutt[79] + aryutt[117] + aryutt[114] + aryutt[112] + aryutt[134]
    + aryutt[137] + aryutt[78] + aryutt[111] + aryutt[106] - 301./32.
    + aryutt[105];
   aryutt[293]= - 1./2.*aryutt[18];
   aryutt[167]=aryutt[14] + aryutt[20] + aryutt[293] - aryutt[92] + 
   aryutt[167];
   aryutt[167]=aryutt[99]*aryutt[167];
   aryutt[306]=aryutt[38] + 13./8.*aryutt[36];
   aryutt[306]=aryutt[12]*aryutt[306];
   aryutt[307]= - 1./2.*aryutt[36];
   aryutt[308]=1 + aryutt[307];
   aryutt[308]=aryutt[15]*aryutt[99]*aryutt[308];
   aryutt[309]=1./2. + aryutt[36];
   aryutt[309]=aryutt[17]*aryutt[99]*aryutt[309];
   aryutt[310]= - 1 - aryutt[84];
   aryutt[310]=MMH*aryutt[99]*aryutt[310];
   aryutt[150]=5./24.*aryutt[310] + 5./12.*aryutt[309] + 5./12.*
   aryutt[308] + 5./12.*aryutt[167] + 1./4.*aryutt[150] + 1./27.*
   aryutt[306];
   aryutt[150]=MMH*aryutt[150];
   aryutt[167]=307./9. + aryutt[274];
   aryutt[167]=aryutt[199] + 1./2.*aryutt[167] + aryutt[40];
   aryutt[167]=aryutt[263] + 1./4.*aryutt[167] + aryutt[272];
   aryutt[167]=aryutt[14]*aryutt[167];
   aryutt[167]=47./72.*aryutt[183] + 1./2.*aryutt[228] + 5./24.*
   aryutt[181] + aryutt[167] + aryutt[116] - 31./24.*aryutt[18] + 
   aryutt[203] + aryutt[160] + 97./24.*aryutt[42];
   aryutt[274]=5./8.*aryutt[8] - 5./8.*aryutt[21] + aryutt[107] - 1./72.
    + 5*aryutt[38];
   aryutt[306]= - 1./9.*aryutt[12];
   aryutt[274]=1./2.*aryutt[274] + aryutt[306];
   aryutt[274]=aryutt[16]*aryutt[274];
   aryutt[308]=331./54. + aryutt[201];
   aryutt[308]= - 125./216.*aryutt[12] + 1./4.*aryutt[308] - aryutt[8];
   aryutt[308]=aryutt[15]*aryutt[308];
   aryutt[167]=1./2.*aryutt[308] + 1./2.*aryutt[167] + 1./3.*
   aryutt[274];
   aryutt[274]=aryutt[225] + 1./12.*aryutt[36] + aryutt[205] + 1./12.*
   aryutt[38];
   aryutt[274]=aryutt[17]*aryutt[274];
   aryutt[150]=1./2.*aryutt[150] + 1./2.*aryutt[167] + aryutt[274];
   aryutt[150]=MMH*aryutt[150];
   aryutt[167]=1./8.*aryutt[85];
   aryutt[226]=aryutt[262] + 1./2.*aryutt[239] + aryutt[258] + 1./2.*
   aryutt[254] + 1./32.*aryutt[261] + 13./16.*aryutt[51] + aryutt[167]
    + aryutt[226] + 2*aryutt[44];
   aryutt[226]=MMt*aryutt[226];
   aryutt[239]= - aryutt[8] + aryutt[263] + 1./4.*aryutt[276] + 
   aryutt[272];
   aryutt[239]=aryutt[8]*aryutt[239];
   aryutt[254]=aryutt[292] + aryutt[36];
   aryutt[254]=aryutt[21]*aryutt[254];
   aryutt[239]=1./4.*aryutt[285] + 1./4.*aryutt[282] + aryutt[239] + 
   aryutt[136] + 1./12.*aryutt[254];
   aryutt[239]=MMH*aryutt[239];
   aryutt[254]=aryutt[238] + aryutt[267] + aryutt[204] + 69./4.*
   aryutt[18] + aryutt[192] - 307./36.*aryutt[19];
   aryutt[258]=205 - 143./2.*aryutt[21];
   aryutt[258]=1./3.*aryutt[258] - 5*aryutt[8];
   aryutt[258]=aryutt[16]*aryutt[258];
   aryutt[254]=1./6.*aryutt[258] + 1./2.*aryutt[254] + aryutt[259];
   aryutt[257]=1./4.*aryutt[257] + 2*aryutt[260];
   aryutt[257]=aryutt[17]*aryutt[257];
   aryutt[258]= - 9 + 37./6.*aryutt[21];
   aryutt[258]=1./16.*aryutt[258] + aryutt[8];
   aryutt[258]=aryutt[15]*aryutt[258];
   aryutt[226]=aryutt[226] + 1./2.*aryutt[239] + aryutt[257] + 1./4.*
   aryutt[254] + aryutt[258];
   aryutt[226]=MMt*aryutt[226];
   aryutt[230]=aryutt[38] + aryutt[230];
   aryutt[239]=aryutt[15]*aryutt[230];
   aryutt[254]= - aryutt[38] + 5./8.*aryutt[36];
   aryutt[254]=aryutt[17]*aryutt[254];
   aryutt[239]=1./9.*aryutt[254] + 5./4.*aryutt[182] + 1./9.*
   aryutt[239];
   aryutt[239]=MMH*aryutt[239];
   aryutt[243]=5./24.*aryutt[15] + 1./8.*aryutt[14] + aryutt[243];
   aryutt[243]=aryutt[15]*aryutt[243];
   aryutt[254]=283./24.*aryutt[15] - 1./8.*aryutt[14] + aryutt[256];
   aryutt[254]=aryutt[17]*aryutt[254];
   aryutt[239]=aryutt[239] + 1./3.*aryutt[254] + 5./4.*aryutt[253] + 1./
   3.*aryutt[243];
   aryutt[239]=aryutt[6]*MMH*aryutt[239];
   aryutt[243]=53*aryutt[14] - 779./6.*aryutt[16];
   aryutt[243]=aryutt[16]*aryutt[243];
   aryutt[243]=aryutt[249] + 1./3.*aryutt[243];
   aryutt[254]=37./2.*aryutt[14] + 61*aryutt[16];
   aryutt[254]=1./3.*aryutt[254] - 15./4.*aryutt[15];
   aryutt[254]=aryutt[15]*aryutt[254];
   aryutt[243]=1./2.*aryutt[243] + aryutt[254];
   aryutt[254]= - 29./3.*aryutt[15] - aryutt[14] - 35./6.*aryutt[16];
   aryutt[254]=1./4.*aryutt[254] + aryutt[17];
   aryutt[254]=aryutt[17]*aryutt[254];
   aryutt[243]=1./8.*aryutt[243] + aryutt[254];
   aryutt[150]=aryutt[226] + 1./12.*aryutt[239] + 1./2.*aryutt[243] + 
   aryutt[150];
   aryutt[149]=aryutt[150] + aryutt[149];
   aryutt[149]=aryutt[210]*aryutt[149];
   aryutt[226]=1./6.*aryutt[269];
   aryutt[105]=aryutt[226] + 5./108.*aryutt[36] + 19./54.*aryutt[38] + 
   aryutt[140] + 7./36.*aryutt[13] + aryutt[139] + aryutt[135] - 41./36.
   *aryutt[84] + aryutt[126] - 7./36.*aryutt[79] + aryutt[117] + 
   aryutt[114] + aryutt[112] + aryutt[134] + aryutt[137] + aryutt[78]
    + aryutt[111] + aryutt[106] - 863./96. + aryutt[105];
   aryutt[106]=aryutt[159] + aryutt[155] + aryutt[215];
   aryutt[106]=aryutt[8]*aryutt[106];
   aryutt[111]=aryutt[272] + aryutt[307];
   aryutt[111]=aryutt[12]*aryutt[111];
   aryutt[105]=25./36.*aryutt[193] + 25./18.*aryutt[178] + 25./18.*
   aryutt[169] + 25./18.*aryutt[163] + 7./36.*aryutt[111] + 1./2.*
   aryutt[105] + aryutt[106];
   aryutt[105]=MMH*aryutt[105];
   aryutt[106]=aryutt[116] + 227./72.*aryutt[18] + aryutt[203] + 
   aryutt[160] + 211./72.*aryutt[42];
   aryutt[111]=aryutt[212] + aryutt[157] + 277./27. + 9./4.*aryutt[35];
   aryutt[111]=1./4.*aryutt[111] + aryutt[272];
   aryutt[111]=aryutt[14]*aryutt[111];
   aryutt[106]=7./108.*aryutt[156] + aryutt[234] + 1./12.*aryutt[181]
    + 1./2.*aryutt[106] + aryutt[111];
   aryutt[111]=aryutt[225] + aryutt[205] + aryutt[278];
   aryutt[111]=aryutt[17]*aryutt[111];
   aryutt[114]= - 7./144.*aryutt[12];
   aryutt[116]=aryutt[114] + 85./32.*aryutt[8] + aryutt[219] - 29./288.
    + aryutt[38];
   aryutt[116]=aryutt[16]*aryutt[116];
   aryutt[117]= - 17./27. - aryutt[38];
   aryutt[117]=7./36.*aryutt[12] + 13./9.*aryutt[8] + 1./2.*aryutt[117]
    + aryutt[36];
   aryutt[117]=aryutt[15]*aryutt[117];
   aryutt[105]=1./4.*aryutt[105] + aryutt[111] + 1./4.*aryutt[117] + 1./
   2.*aryutt[106] + 1./3.*aryutt[116];
   aryutt[105]=MMH*aryutt[105];
   aryutt[106]=aryutt[272] - aryutt[36];
   aryutt[111]=aryutt[15]*aryutt[106];
   aryutt[116]=aryutt[215] + aryutt[36];
   aryutt[116]=aryutt[17]*aryutt[116];
   aryutt[111]=17./12.*aryutt[116] + aryutt[182] + 17./12.*aryutt[111];
   aryutt[111]=MMH*aryutt[111];
   aryutt[116]= - 1./2.*aryutt[16];
   aryutt[117]= - aryutt[14] + aryutt[116];
   aryutt[126]=1./3.*aryutt[117] + aryutt[188];
   aryutt[126]=aryutt[15]*aryutt[126];
   aryutt[134]=aryutt[14] + 1./2.*aryutt[16];
   aryutt[135]=17./3.*aryutt[134] + aryutt[195];
   aryutt[135]=aryutt[17]*aryutt[135];
   aryutt[111]=aryutt[111] + 1./6.*aryutt[135] + aryutt[253] + 17./6.*
   aryutt[126];
   aryutt[111]=aryutt[6]*MMH*aryutt[111];
   aryutt[126]= - 7*aryutt[14] - 53./6.*aryutt[16];
   aryutt[126]=aryutt[16]*aryutt[126];
   aryutt[135]= - 11*aryutt[14];
   aryutt[137]= - 19./2.*aryutt[15] + aryutt[135] + 31*aryutt[16];
   aryutt[137]=aryutt[15]*aryutt[137];
   aryutt[126]=aryutt[137] + aryutt[249] + 17./3.*aryutt[126];
   aryutt[137]=1./9.*aryutt[15] - aryutt[14] - 115./6.*aryutt[16];
   aryutt[137]=1./4.*aryutt[137] + aryutt[17];
   aryutt[137]=aryutt[17]*aryutt[137];
   aryutt[126]=1./16.*aryutt[126] + aryutt[137];
   aryutt[105]=1./12.*aryutt[111] + 1./2.*aryutt[126] + aryutt[105];
   aryutt[111]=aryutt[238] + aryutt[267] + aryutt[204] - 41./18.*
   aryutt[18] + aryutt[192] - 59./9.*aryutt[19];
   aryutt[126]=107 - 65*aryutt[21];
   aryutt[126]=7./12.*aryutt[126] - 85*aryutt[8];
   aryutt[126]=aryutt[16]*aryutt[126];
   aryutt[111]=1./6.*aryutt[126] + 1./2.*aryutt[111] + aryutt[259];
   aryutt[126]= - 1./16. + 11./3.*aryutt[21];
   aryutt[126]=1./2.*aryutt[126] - 13./3.*aryutt[8];
   aryutt[126]=aryutt[15]*aryutt[126];
   aryutt[111]=1./4.*aryutt[111] + 1./3.*aryutt[126];
   aryutt[126]= - 1./4. + aryutt[38];
   aryutt[126]=aryutt[21]*aryutt[126];
   aryutt[136]=aryutt[136] + 1./6.*aryutt[126];
   aryutt[137]=aryutt[244] + aryutt[280] + aryutt[272];
   aryutt[137]=aryutt[8]*aryutt[137];
   aryutt[136]=aryutt[286] + aryutt[284] + 1./2.*aryutt[136] + 
   aryutt[137];
   aryutt[136]=MMH*aryutt[136];
   aryutt[111]=aryutt[266] + 1./2.*aryutt[136] + 1./2.*aryutt[111] + 
   aryutt[281];
   aryutt[111]=MMt*aryutt[111];
   aryutt[136]=aryutt[14] - aryutt[16];
   aryutt[137]=aryutt[16]*aryutt[136];
   aryutt[139]=aryutt[15]*aryutt[252];
   aryutt[137]=1./4.*aryutt[296] + 1./4.*aryutt[137] + aryutt[139];
   aryutt[139]=aryutt[14] + 7*aryutt[183];
   aryutt[140]= - 7./9.*aryutt[12];
   aryutt[155]=aryutt[140] - 1./9. - aryutt[8];
   aryutt[155]=aryutt[16]*aryutt[155];
   aryutt[139]=1./9.*aryutt[139] + aryutt[155];
   aryutt[155]=1./6.*aryutt[8];
   aryutt[160]= - aryutt[38] + aryutt[155];
   aryutt[160]=aryutt[15]*aryutt[160];
   aryutt[163]=aryutt[12]*aryutt[38];
   aryutt[169]=aryutt[38] + 7*aryutt[163];
   aryutt[169]=MMH*aryutt[169];
   aryutt[139]=1./54.*aryutt[169] + 1./6.*aryutt[139] + aryutt[160];
   aryutt[139]=MMH*aryutt[139];
   aryutt[160]=aryutt[15]*aryutt[38];
   aryutt[169]= - aryutt[17]*aryutt[38];
   aryutt[160]=aryutt[160] + aryutt[169];
   aryutt[160]=MMH*aryutt[160];
   aryutt[178]=aryutt[15]*aryutt[136];
   aryutt[182]=aryutt[17]*aryutt[252];
   aryutt[160]=aryutt[160] + aryutt[178] + aryutt[182];
   aryutt[160]=aryutt[6]*MMH*aryutt[160];
   aryutt[178]=1./3.*aryutt[21];
   aryutt[192]= - 1./2. + aryutt[178];
   aryutt[193]= - 1./3.*aryutt[8];
   aryutt[192]=1./2.*aryutt[192] + aryutt[193];
   aryutt[192]=aryutt[15]*aryutt[192];
   aryutt[195]= - 1./3.*aryutt[21];
   aryutt[203]=1./2. + aryutt[195];
   aryutt[204]=1./3.*aryutt[8];
   aryutt[203]=1./2.*aryutt[203] + aryutt[204];
   aryutt[203]=aryutt[16]*aryutt[203];
   aryutt[205]= - aryutt[19] + aryutt[18];
   aryutt[192]=aryutt[192] + 1./4.*aryutt[205] + aryutt[203];
   aryutt[192]=MMt*aryutt[192];
   aryutt[137]=1./2.*aryutt[192] + 17./216.*aryutt[160] + 1./3.*
   aryutt[137] + 1./4.*aryutt[139];
   aryutt[137]=aryutt[176]*aryutt[137];
   aryutt[105]=1./4.*aryutt[137] + 1./2.*aryutt[105] + aryutt[111];
   aryutt[105]=aryutt[176]*aryutt[105];
   aryutt[105]=aryutt[150] + aryutt[105];
   aryutt[105]=aryutt[176]*aryutt[105];
   aryutt[111]=aryutt[8]*aryutt[38];
   aryutt[137]= - 1 - aryutt[85];
   aryutt[139]=aryutt[137] + 1./6.*aryutt[111];
   aryutt[139]=MMH*aryutt[139];
   aryutt[150]= - 1./2.*aryutt[19];
   aryutt[127]=aryutt[150] - aryutt[95] + aryutt[127];
   aryutt[144]=1 + aryutt[144];
   aryutt[160]=aryutt[14]*aryutt[144];
   aryutt[192]= - 1./12.*aryutt[8] + 1 - 1./12.*aryutt[38];
   aryutt[192]=aryutt[16]*aryutt[192];
   aryutt[139]=1./2.*aryutt[139] + 1./6.*aryutt[169] + aryutt[192] + 1./
   12.*aryutt[259] + aryutt[127] + aryutt[160];
   aryutt[139]=MMH*aryutt[139];
   aryutt[160]=7*aryutt[14];
   aryutt[192]= - 5*aryutt[16];
   aryutt[205]=aryutt[160] + aryutt[192];
   aryutt[205]=aryutt[16]*aryutt[205];
   aryutt[205]=aryutt[182] - aryutt[236] + 1./2.*aryutt[205];
   aryutt[139]=1./6.*aryutt[205] + aryutt[139];
   aryutt[139]=MMH*aryutt[139];
   aryutt[111]=aryutt[137] + 1./3.*aryutt[111];
   aryutt[111]=MMH*aryutt[111];
   aryutt[205]=1 + aryutt[215];
   aryutt[205]=aryutt[14]*aryutt[205];
   aryutt[144]=aryutt[144] - 1./6.*aryutt[8];
   aryutt[144]=aryutt[16]*aryutt[144];
   aryutt[215]=1./6.*aryutt[259];
   aryutt[111]=1./2.*aryutt[111] + 1./3.*aryutt[169] + aryutt[144] + 
   aryutt[215] + aryutt[127] + aryutt[205];
   aryutt[111]=MMH*aryutt[111];
   aryutt[116]=aryutt[14] + aryutt[116];
   aryutt[116]=aryutt[16]*aryutt[116];
   aryutt[116]=1./2.*aryutt[182] - 1./2.*aryutt[236] + aryutt[116];
   aryutt[111]=1./3.*aryutt[116] + 1./2.*aryutt[111];
   aryutt[111]=MMH*aryutt[111];
   aryutt[116]=1./3.*aryutt[181];
   aryutt[144]= - aryutt[14] + aryutt[42] - aryutt[19];
   aryutt[169]=1./2.*aryutt[144] + aryutt[116];
   aryutt[182]= - aryutt[23] + aryutt[85] + 1./4. + aryutt[82];
   aryutt[205]=aryutt[182] + 1./3.*aryutt[269];
   aryutt[225]= - 1./3.*aryutt[8]*aryutt[38];
   aryutt[205]=1./4.*aryutt[205] + aryutt[225];
   aryutt[205]=MMH*aryutt[205];
   aryutt[234]=1./3.*aryutt[228];
   aryutt[169]=aryutt[205] + aryutt[203] + 1./2.*aryutt[169] + 
   aryutt[234];
   aryutt[169]=MMt*MMH*aryutt[169];
   aryutt[111]=aryutt[111] + aryutt[169];
   aryutt[111]=aryutt[176]*aryutt[111];
   aryutt[169]=aryutt[144] + 5./6.*aryutt[181];
   aryutt[203]=1 - 5./6.*aryutt[21];
   aryutt[203]=1./2.*aryutt[203] + aryutt[204];
   aryutt[203]=aryutt[16]*aryutt[203];
   aryutt[205]=aryutt[182] + aryutt[226];
   aryutt[205]=1./2.*aryutt[205] + aryutt[225];
   aryutt[205]=MMH*aryutt[205];
   aryutt[169]=aryutt[205] + aryutt[203] + 1./2.*aryutt[169] + 
   aryutt[234];
   aryutt[169]=MMt*MMH*aryutt[169];
   aryutt[111]=aryutt[111] + aryutt[139] + aryutt[169];
   aryutt[111]=aryutt[176]*aryutt[111];
   aryutt[139]=aryutt[8]*aryutt[36];
   aryutt[139]=3*aryutt[137] + 1./3.*aryutt[139];
   aryutt[139]=MMH*aryutt[139];
   aryutt[169]=3 + aryutt[275];
   aryutt[169]=aryutt[14]*aryutt[169];
   aryutt[203]= - 1./6.*aryutt[36];
   aryutt[205]=3 + aryutt[203];
   aryutt[205]=aryutt[16]*aryutt[205];
   aryutt[225]= - aryutt[17]*aryutt[36];
   aryutt[226]= - aryutt[15]*aryutt[8];
   aryutt[234]=1./6.*aryutt[226];
   aryutt[139]=1./2.*aryutt[139] + 1./3.*aryutt[225] + aryutt[234] + 
   aryutt[205] + aryutt[215] + 3*aryutt[127] + aryutt[169];
   aryutt[139]=MMH*aryutt[139];
   aryutt[169]=1./3.*aryutt[14] - 3./8.*aryutt[16];
   aryutt[169]=aryutt[16]*aryutt[169];
   aryutt[205]= - aryutt[14] + aryutt[15];
   aryutt[205]=aryutt[17]*aryutt[205];
   aryutt[134]=aryutt[15]*aryutt[134];
   aryutt[139]=1./4.*aryutt[139] + 1./12.*aryutt[205] + 1./12.*
   aryutt[134] - 1./12.*aryutt[236] + aryutt[169];
   aryutt[139]=MMH*aryutt[139];
   aryutt[169]= - aryutt[21]*aryutt[36];
   aryutt[169]=3*aryutt[182] + 1./3.*aryutt[169];
   aryutt[205]= - aryutt[8]*aryutt[36];
   aryutt[169]=1./4.*aryutt[169] + 1./3.*aryutt[205];
   aryutt[169]=MMH*aryutt[169];
   aryutt[205]=1 - aryutt[21];
   aryutt[205]=aryutt[16]*aryutt[205];
   aryutt[215]=aryutt[15]*aryutt[224];
   aryutt[116]=1./2.*aryutt[169] + 1./6.*aryutt[215] + 3./8.*
   aryutt[205] + 1./6.*aryutt[228] + 3./8.*aryutt[144] + aryutt[116];
   aryutt[116]=MMt*MMH*aryutt[116];
   aryutt[116]=aryutt[139] + aryutt[116];
   aryutt[111]=aryutt[116] + 1./2.*aryutt[111];
   aryutt[111]=aryutt[176]*aryutt[111];
   aryutt[139]=aryutt[155] + aryutt[203] + 1 + aryutt[278];
   aryutt[139]=aryutt[16]*aryutt[139];
   aryutt[155]=aryutt[275] + 1 + aryutt[272];
   aryutt[155]=aryutt[14]*aryutt[155];
   aryutt[137]=aryutt[137] + 1./3.*aryutt[305];
   aryutt[137]=MMH*aryutt[137];
   aryutt[127]=1./2.*aryutt[137] + 1./3.*aryutt[298] + aryutt[234] + 
   aryutt[139] + aryutt[127] + aryutt[155];
   aryutt[127]=MMH*aryutt[127];
   aryutt[137]=1./4.*aryutt[14];
   aryutt[139]=aryutt[137] - aryutt[16];
   aryutt[139]=aryutt[16]*aryutt[139];
   aryutt[139]=aryutt[301] + aryutt[139] + 1./2.*aryutt[134];
   aryutt[127]=1./3.*aryutt[139] + 1./2.*aryutt[127];
   aryutt[127]=MMH*aryutt[127];
   aryutt[139]=aryutt[14]*aryutt[214];
   aryutt[155]=aryutt[214] + aryutt[8];
   aryutt[155]=aryutt[16]*aryutt[155];
   aryutt[169]=MMH*aryutt[305];
   aryutt[139]=1./2.*aryutt[169] + aryutt[298] + 1./2.*aryutt[226] + 
   aryutt[139] + 1./2.*aryutt[155];
   aryutt[139]=MMH*aryutt[139];
   aryutt[117]=aryutt[16]*aryutt[117];
   aryutt[117]=aryutt[139] + aryutt[300] + aryutt[117] + aryutt[134];
   aryutt[117]=MMH*aryutt[117];
   aryutt[134]=1./4.*aryutt[216] + aryutt[290];
   aryutt[134]=MMH*aryutt[134];
   aryutt[139]=aryutt[16]*aryutt[302];
   aryutt[134]=aryutt[134] + aryutt[139] + aryutt[215];
   aryutt[134]=MMt*MMH*aryutt[134];
   aryutt[117]=1./2.*aryutt[117] + aryutt[134];
   aryutt[117]=aryutt[210]*aryutt[117];
   aryutt[134]=aryutt[144] + aryutt[181];
   aryutt[139]=aryutt[193] + 1./4. + aryutt[195];
   aryutt[139]=aryutt[16]*aryutt[139];
   aryutt[144]=aryutt[182] + aryutt[287];
   aryutt[144]=1./4.*aryutt[144] + aryutt[291];
   aryutt[144]=MMH*aryutt[144];
   aryutt[134]=aryutt[144] + 1./3.*aryutt[215] + 1./4.*aryutt[134] + 
   aryutt[139];
   aryutt[134]=MMt*MMH*aryutt[134];
   aryutt[117]=1./3.*aryutt[117] + aryutt[127] + aryutt[134];
   aryutt[117]=aryutt[210]*aryutt[117];
   aryutt[116]=aryutt[116] + 1./2.*aryutt[117];
   aryutt[116]=aryutt[210]*aryutt[116];
   aryutt[111]=aryutt[111] + aryutt[116];
   aryutt[111]=aryutt[9]*aryutt[111];
   aryutt[105]=1./4.*aryutt[111] + aryutt[105] + aryutt[149];
   aryutt[105]=aryutt[9]*aryutt[105];
   aryutt[111]=77./36.*aryutt[76] + 691./144. + 15*aryutt[83];
   aryutt[111]= - 1./4.*aryutt[90] - 1./6.*aryutt[141] + 1./4.*
   aryutt[111] - 1./3.*aryutt[80];
   aryutt[116]=1./2.*aryutt[75];
   aryutt[111]=239./144.*aryutt[87] + 7./36.*aryutt[72] + 7./36.*
   aryutt[74] + aryutt[116] + 1./2.*aryutt[111] + aryutt[77];
   aryutt[111]= - 23./72.*aryutt[88] + 1./2.*aryutt[111] + 17./9.*
   aryutt[44];
   aryutt[117]= - 17./3.*aryutt[39];
   aryutt[127]=7./4. + aryutt[117];
   aryutt[127]=1./2.*aryutt[127] + aryutt[178];
   aryutt[127]=aryutt[21]*aryutt[127];
   aryutt[134]= - 9 - 7./36.*aryutt[21];
   aryutt[139]=1./4.*aryutt[12];
   aryutt[134]=aryutt[139] + 1./8.*aryutt[134] - 7./9.*aryutt[8];
   aryutt[134]=aryutt[12]*aryutt[134];
   aryutt[144]=7*aryutt[156] + 61*aryutt[42] + aryutt[160];
   aryutt[149]= - 1 + aryutt[202];
   aryutt[149]=1./2.*aryutt[149] + aryutt[8];
   aryutt[155]=aryutt[15]*aryutt[149];
   aryutt[160]=3*aryutt[155];
   aryutt[144]=1./18.*aryutt[144] + aryutt[160];
   aryutt[144]=aryutt[3]*aryutt[144];
   aryutt[169]= - 1./2.*aryutt[68] - aryutt[65];
   aryutt[169]=1./2.*aryutt[169] - 7./3.*aryutt[63];
   aryutt[178]= - 1./3.*aryutt[67];
   aryutt[169]=1./2.*aryutt[169] + aryutt[178];
   aryutt[169]=MMt*aryutt[169];
   aryutt[182]=77./2. + aryutt[110];
   aryutt[182]=1./3.*aryutt[182] - aryutt[37];
   aryutt[182]=aryutt[8]*aryutt[182];
   aryutt[193]=13./6.*aryutt[63] - aryutt[64];
   aryutt[193]=1./2.*aryutt[193] - aryutt[60];
   aryutt[193]=MMH*aryutt[193];
   aryutt[195]=aryutt[101] + aryutt[283];
   aryutt[195]=aryutt[15]*aryutt[195];
   aryutt[111]=1./6.*aryutt[169] + 1./4.*aryutt[193] + 1./4.*
   aryutt[144] + 3./8.*aryutt[195] + 1./4.*aryutt[134] + 1./8.*
   aryutt[182] + 1./16.*aryutt[127] + 1./8.*aryutt[37] - 55./384.*
   aryutt[13] + aryutt[148] - 1./18.*aryutt[25] - 15./32.*aryutt[23] - 
   7./72.*aryutt[84] - 43./144.*aryutt[51] + 1./2.*aryutt[111] - 4./9.*
   aryutt[45];
   aryutt[111]=MMt*aryutt[111];
   aryutt[127]=MMH*aryutt[163];
   aryutt[127]=aryutt[127] + aryutt[183] + aryutt[200];
   aryutt[127]=MMH*aryutt[127];
   aryutt[134]=aryutt[15]*aryutt[299];
   aryutt[127]=aryutt[127] + aryutt[134] + aryutt[296];
   aryutt[127]=aryutt[6]*aryutt[127];
   aryutt[134]=1./3.*aryutt[146] - aryutt[23] + aryutt[87] + 1./4. + 
   aryutt[83];
   aryutt[144]= - aryutt[8]*aryutt[39];
   aryutt[134]=1./4.*aryutt[134] + 1./3.*aryutt[144];
   aryutt[134]=MMt*aryutt[134];
   aryutt[144]=aryutt[150] + 1./2.*aryutt[50] - aryutt[96];
   aryutt[146]= - 1 + aryutt[171];
   aryutt[146]=aryutt[14]*aryutt[146];
   aryutt[140]=aryutt[140] + 89./9. - aryutt[39];
   aryutt[140]=aryutt[16]*aryutt[140];
   aryutt[148]=55 + 7*aryutt[12];
   aryutt[148]=aryutt[15]*aryutt[148];
   aryutt[150]=aryutt[3]*aryutt[15]*aryutt[295];
   aryutt[163]=aryutt[8]*aryutt[39];
   aryutt[163]= - aryutt[38] + 1./4.*aryutt[163];
   aryutt[163]=MMH*aryutt[163];
   aryutt[127]=aryutt[134] + 17./108.*aryutt[127] + 1./3.*aryutt[163]
    + 1./6.*aryutt[187] + 1./2.*aryutt[150] + 1./108.*aryutt[148] + 1./
   12.*aryutt[140] + 1./2.*aryutt[144] + 1./3.*aryutt[146];
   aryutt[127]=aryutt[176]*aryutt[127];
   aryutt[134]= - 7./2.*aryutt[92] + aryutt[48];
   aryutt[138]= - 1./3.*aryutt[18] + 1./3.*aryutt[134] + aryutt[138];
   aryutt[138]=7./12.*aryutt[14] + 1./2.*aryutt[138] + 1./3.*aryutt[20]
   ;
   aryutt[138]=aryutt[99]*aryutt[138];
   aryutt[140]=aryutt[12]*aryutt[36];
   aryutt[144]=1./6.*aryutt[140] + aryutt[203] + 1./6.*aryutt[13] - 7./
   3.*aryutt[84] - 1./6.*aryutt[79] - 31./16. + 1./3.*aryutt[88];
   aryutt[144]=aryutt[99]*aryutt[144];
   aryutt[144]= - 1./3.*aryutt[63] + 25./16.*aryutt[144];
   aryutt[144]=MMH*aryutt[144];
   aryutt[146]= - 3*aryutt[77];
   aryutt[148]= - 97./36.*aryutt[72] - 97./36.*aryutt[74] - 2213./432.
    + aryutt[146];
   aryutt[148]=25./36.*aryutt[79] + 91./72.*aryutt[88] + 1./4.*
   aryutt[148] - 17./9.*aryutt[44];
   aryutt[150]=17*aryutt[39];
   aryutt[163]= - 163./6. + aryutt[150];
   aryutt[163]=1./12.*aryutt[163] + aryutt[37];
   aryutt[163]=aryutt[8]*aryutt[163];
   aryutt[169]=151./4. + aryutt[206];
   aryutt[169]=1./3.*aryutt[169] + 7./2.*aryutt[36];
   aryutt[169]=1./3.*aryutt[169] - 23./16.*aryutt[8];
   aryutt[169]=aryutt[12]*aryutt[169];
   aryutt[171]=17./4. - aryutt[36];
   aryutt[171]=aryutt[15]*aryutt[99]*aryutt[171];
   aryutt[182]=7./8. + aryutt[36];
   aryutt[182]=aryutt[17]*aryutt[99]*aryutt[182];
   aryutt[138]=1./6.*aryutt[144] + 25./36.*aryutt[182] + 25./72.*
   aryutt[171] + 25./12.*aryutt[138] + 1./24.*aryutt[169] + 1./8.*
   aryutt[163] + 29./216.*aryutt[36] - 5./108.*aryutt[38] - 3./16.*
   aryutt[37] - 1./384.*aryutt[13] + 209./1728.*aryutt[84] + 205./1152.
   *aryutt[51] + 1./8.*aryutt[148] + 1./9.*aryutt[45];
   aryutt[138]=MMH*aryutt[138];
   aryutt[144]=47./6. + aryutt[110];
   aryutt[144]=1./3.*aryutt[144] - aryutt[37];
   aryutt[144]=aryutt[14]*aryutt[144];
   aryutt[110]=661./9. + aryutt[110];
   aryutt[110]=403./36.*aryutt[12] + 1./2.*aryutt[110] - aryutt[21];
   aryutt[110]=aryutt[16]*aryutt[110];
   aryutt[148]=3./4.*aryutt[48];
   aryutt[163]=aryutt[99]*aryutt[186];
   aryutt[110]=25./18.*aryutt[163] + 1./12.*aryutt[110] + 11./72.*
   aryutt[156] + 1./4.*aryutt[144] + 383./288.*aryutt[20] - 35507./5184.
   *aryutt[18] - 103./48.*aryutt[19] - 103./432.*aryutt[42] - 1./36.*
   aryutt[5] + 107./144.*aryutt[33] - 155./288.*aryutt[34] - 23./18.*
   aryutt[32] + aryutt[148] + 11./36.*aryutt[92] - 7./12.*aryutt[93] - 
   413./144.*aryutt[96] + 15./8.*aryutt[50] + 1./3.*aryutt[94];
   aryutt[144]= - 27*aryutt[35];
   aryutt[169]=29653./648. + aryutt[144];
   aryutt[187]= - 3*aryutt[40];
   aryutt[169]=aryutt[118] + 1./2.*aryutt[169] + aryutt[187];
   aryutt[169]=1./4.*aryutt[169] + aryutt[38];
   aryutt[193]= - 3*aryutt[101];
   aryutt[203]=aryutt[193] + 625./9.*aryutt[99];
   aryutt[203]=aryutt[15]*aryutt[203];
   aryutt[205]= - aryutt[16]*aryutt[98];
   aryutt[169]=1./8.*aryutt[203] + 1./48.*aryutt[205] + 13./1296.*
   aryutt[12] - 197./144.*aryutt[8] - 103./1152.*aryutt[21] + 1./2.*
   aryutt[169] - aryutt[36];
   aryutt[169]=aryutt[15]*aryutt[169];
   aryutt[117]=19./2. + aryutt[117];
   aryutt[117]=79./36.*aryutt[12] + 17./18.*aryutt[8] + 97./144.*
   aryutt[21] + 1./2.*aryutt[117] - aryutt[37];
   aryutt[203]=aryutt[193] + 4825./54.*aryutt[99];
   aryutt[203]=aryutt[15]*aryutt[203];
   aryutt[206]=aryutt[99]*aryutt[16];
   aryutt[117]=1./2.*aryutt[203] + 1./2.*aryutt[117] + 25./9.*
   aryutt[206];
   aryutt[203]=9./8.*aryutt[15];
   aryutt[215]= - 17./9.*aryutt[14] + aryutt[203];
   aryutt[215]=aryutt[3]*aryutt[215];
   aryutt[216]= - aryutt[17]*aryutt[99];
   aryutt[117]=1075./216.*aryutt[216] + 1./2.*aryutt[117] + aryutt[215]
   ;
   aryutt[117]=aryutt[17]*aryutt[117];
   aryutt[215]= - aryutt[14] - 25./4.*aryutt[16];
   aryutt[215]=3*aryutt[215] - 7./2.*aryutt[15];
   aryutt[215]=aryutt[3]*aryutt[15]*aryutt[215];
   aryutt[110]=aryutt[117] + 1./4.*aryutt[215] + 1./2.*aryutt[110] + 
   aryutt[169];
   aryutt[117]=13./3.*aryutt[92] - aryutt[48];
   aryutt[117]=1./3.*aryutt[117] + aryutt[165];
   aryutt[169]=aryutt[211] + 85./24.*aryutt[8] + 17./3.*aryutt[36] - 37.
   /4. - 17./9.*aryutt[38];
   aryutt[169]=aryutt[15]*aryutt[169];
   aryutt[117]=1./2.*aryutt[169] + 17./72.*aryutt[200] + 17./36.*
   aryutt[156] - 187./72.*aryutt[14] - 5./18.*aryutt[20] + 323./144.*
   aryutt[18] + 17./8.*aryutt[117] - 1./3.*aryutt[33];
   aryutt[143]= - 125./8. + aryutt[143];
   aryutt[143]= - 85./8.*aryutt[8] + 1./3.*aryutt[143] + aryutt[227];
   aryutt[143]=1./4.*aryutt[143] - aryutt[12];
   aryutt[143]=aryutt[17]*aryutt[143];
   aryutt[106]=aryutt[12]*aryutt[106];
   aryutt[106]=11./6.*aryutt[185] + aryutt[106];
   aryutt[106]=MMH*aryutt[106];
   aryutt[106]=17./48.*aryutt[106] + 1./2.*aryutt[117] + 1./3.*
   aryutt[143];
   aryutt[106]=MMH*aryutt[106];
   aryutt[117]= - 5423./3.*aryutt[15] + aryutt[289] + 1085*aryutt[16];
   aryutt[117]=aryutt[15]*aryutt[117];
   aryutt[117]= - 11*aryutt[186] + 1./12.*aryutt[117];
   aryutt[143]= - 85*aryutt[14] - 1481./2.*aryutt[16];
   aryutt[143]=355./2.*aryutt[17] + 1./4.*aryutt[143] - 1531./3.*
   aryutt[15];
   aryutt[143]=aryutt[17]*aryutt[143];
   aryutt[117]=1./2.*aryutt[117] + 1./3.*aryutt[143];
   aryutt[106]=1./12.*aryutt[117] + aryutt[106];
   aryutt[106]=aryutt[6]*aryutt[106];
   aryutt[106]=1./8.*aryutt[127] + aryutt[111] + 1./6.*aryutt[106] + 1./
   2.*aryutt[110] + aryutt[138];
   aryutt[106]=aryutt[176]*aryutt[106];
   aryutt[110]=11./24.*aryutt[58] + 8185./768. - aryutt[81];
   aryutt[111]=11./12.*aryutt[84];
   aryutt[115]=aryutt[115] - aryutt[19];
   aryutt[115]=aryutt[100]*aryutt[115];
   aryutt[110]=1./16.*aryutt[166] + aryutt[212] + 1./8.*aryutt[115] + 
   73./192.*aryutt[13] + 1./4.*aryutt[40] + aryutt[303] + 25./12.*
   aryutt[25] + 15./32.*aryutt[23] + aryutt[111] - 403./288.*aryutt[51]
    + aryutt[45] + 19./24.*aryutt[88] + 3./16.*aryutt[85] + 10./9.*
   aryutt[44] + 13./48.*aryutt[87] - 11./24.*aryutt[72] - 11./24.*
   aryutt[74] + aryutt[191] + aryutt[112] - 33./64.*aryutt[90] - 13./36.
   *aryutt[141] - 15./64.*aryutt[89] + 1./6.*aryutt[80] - 11./64.*
   aryutt[76] + 3./4.*aryutt[78] - 9./32.*aryutt[73] + 1./3.*
   aryutt[110] - 15./32.*aryutt[83];
   aryutt[112]=295./9. - aryutt[40];
   aryutt[117]= - 1./12.*aryutt[41];
   aryutt[127]=13./3.*aryutt[39];
   aryutt[138]=7./16.*aryutt[21];
   aryutt[112]=aryutt[138] - aryutt[36] - aryutt[38] + aryutt[117] + 
   aryutt[129] + 1./2.*aryutt[112] + aryutt[127];
   aryutt[112]=aryutt[8]*aryutt[112];
   aryutt[143]=13*aryutt[39];
   aryutt[169]= - 1./4.*aryutt[41] + 329./12. + aryutt[143];
   aryutt[185]=235./72.*aryutt[21];
   aryutt[169]=aryutt[185] + aryutt[246] - aryutt[36] + 1./3.*
   aryutt[169] - aryutt[38];
   aryutt[169]=aryutt[21]*aryutt[169];
   aryutt[191]=31./3.*aryutt[42] - 5*aryutt[14];
   aryutt[191]=5./3.*aryutt[191] + aryutt[181];
   aryutt[191]=1./2.*aryutt[191] + 11./3.*aryutt[183];
   aryutt[149]=3*aryutt[16]*aryutt[149];
   aryutt[160]=aryutt[160] + 1./2.*aryutt[191] + aryutt[149];
   aryutt[160]=1./2.*aryutt[3]*aryutt[160];
   aryutt[191]=125./3. + 33./4.*aryutt[21];
   aryutt[191]=aryutt[139] + 1./8.*aryutt[191] + 11./3.*aryutt[8];
   aryutt[191]=1./2.*aryutt[12]*aryutt[191];
   aryutt[161]=7*aryutt[100] + aryutt[161];
   aryutt[200]= - aryutt[8]*aryutt[100];
   aryutt[161]=1./2.*aryutt[161] + 3*aryutt[200];
   aryutt[161]=aryutt[16]*aryutt[161];
   aryutt[200]=1./4.*aryutt[161];
   aryutt[215]=1./4.*aryutt[70];
   aryutt[224]=1./3.*aryutt[52];
   aryutt[225]=1./2.*aryutt[65] - 1./8.*aryutt[68] - aryutt[69] + 
   aryutt[224] - 16./3.*aryutt[53] + aryutt[215];
   aryutt[172]=aryutt[172] + 11./6.*aryutt[63] + 1./3.*aryutt[225] + 1./
   4.*aryutt[61];
   aryutt[172]=MMt*aryutt[172];
   aryutt[225]= - 4./3.*aryutt[56] + aryutt[53];
   aryutt[225]= - 1./2.*aryutt[60] - 1./4.*aryutt[64] - 9./8.*
   aryutt[63] + 4./3.*aryutt[225] - 3./8.*aryutt[61];
   aryutt[225]=MMH*aryutt[225];
   aryutt[195]=3./4.*aryutt[195];
   aryutt[112]=aryutt[172] + aryutt[225] + aryutt[160] + aryutt[195] + 
   aryutt[200] + aryutt[191] + 1./2.*aryutt[112] + aryutt[110] + 1./8.*
   aryutt[169];
   aryutt[112]=MMt*aryutt[112];
   aryutt[169]= - 541./54. + aryutt[144];
   aryutt[226]= - 5*aryutt[38];
   aryutt[227]= - 13./2.*aryutt[8];
   aryutt[169]= - 925./108.*aryutt[12] + aryutt[227] + 45./8.*
   aryutt[21] + aryutt[120] - aryutt[36] + aryutt[226] + aryutt[117] + 
   aryutt[118] + aryutt[127] + 1./2.*aryutt[169] + aryutt[187];
   aryutt[234]=51./16.*aryutt[98] - aryutt[100];
   aryutt[234]=aryutt[16]*aryutt[234];
   aryutt[169]=1./4.*aryutt[169] + aryutt[234];
   aryutt[169]=aryutt[16]*aryutt[169];
   aryutt[238]= - 3*aryutt[73];
   aryutt[239]= - 1289./12. + aryutt[238];
   aryutt[243]=9./4.*aryutt[89];
   aryutt[146]= - 19./12.*aryutt[72] - 19./12.*aryutt[74] + aryutt[146]
    + aryutt[197] + aryutt[243] + 1./4.*aryutt[239] + aryutt[147];
   aryutt[146]= - 5./12.*aryutt[79] - 23./24.*aryutt[88] - 5./24.*
   aryutt[85] + 1./4.*aryutt[146] - 7./3.*aryutt[44];
   aryutt[146]=1./24.*aryutt[126] + 11./27.*aryutt[36] + 13./108.*
   aryutt[38] - 3./8.*aryutt[37] + 13./64.*aryutt[13] - 3./8.*
   aryutt[40] - 4./9.*aryutt[25] + 1./32.*aryutt[23] - 277./288.*
   aryutt[84] + 139./144.*aryutt[51] + 1./4.*aryutt[146] + aryutt[265];
   aryutt[239]= - 13./6.*aryutt[39];
   aryutt[244]= - 1./2.*aryutt[21];
   aryutt[249]=aryutt[244] + aryutt[107] + aryutt[235] + 1./24.*
   aryutt[41] + aryutt[37] + aryutt[239] - 365./36. + aryutt[40];
   aryutt[249]=aryutt[8]*aryutt[249];
   aryutt[252]=7./2.*aryutt[92] - aryutt[48];
   aryutt[253]=1./3.*aryutt[18];
   aryutt[252]=aryutt[253] + 1./3.*aryutt[252] + aryutt[125];
   aryutt[252]= - 7./12.*aryutt[14] + 1./2.*aryutt[252] - 1./3.*
   aryutt[20];
   aryutt[252]=5./2.*aryutt[99]*aryutt[252];
   aryutt[254]=2*aryutt[56] - aryutt[53];
   aryutt[254]=2./3.*aryutt[254] + aryutt[277];
   aryutt[256]= - aryutt[12]*aryutt[36];
   aryutt[222]=1./6.*aryutt[256] + aryutt[222] - 1./6.*aryutt[13] + 7./
   3.*aryutt[84] + 1./6.*aryutt[79] + 31./16. - 1./3.*aryutt[88];
   aryutt[222]=aryutt[99]*aryutt[222];
   aryutt[222]=1./3.*aryutt[254] + 5./16.*aryutt[222];
   aryutt[222]=MMH*aryutt[222];
   aryutt[254]= - 83./2. - 7./3.*aryutt[38];
   aryutt[254]=1./2.*aryutt[254] - 13./3.*aryutt[36];
   aryutt[254]=1./9.*aryutt[254] - 21./16.*aryutt[8];
   aryutt[254]=1./4.*aryutt[12]*aryutt[254];
   aryutt[256]= - 17./4. + aryutt[36];
   aryutt[256]=5./12.*aryutt[15]*aryutt[99]*aryutt[256];
   aryutt[257]= - 7./8. - aryutt[36];
   aryutt[257]=5./6.*aryutt[17]*aryutt[99]*aryutt[257];
   aryutt[249]=aryutt[222] + aryutt[257] + aryutt[256] + aryutt[252] + 
   aryutt[254] + aryutt[146] + 1./4.*aryutt[249];
   aryutt[249]=MMH*aryutt[249];
   aryutt[258]=349./144.*aryutt[21];
   aryutt[259]=109./18.*aryutt[8];
   aryutt[260]=13./12.*aryutt[12];
   aryutt[261]=23./4.*aryutt[98] - aryutt[100];
   aryutt[261]=3*aryutt[16]*aryutt[261];
   aryutt[262]=aryutt[261] + aryutt[260] + aryutt[259] + aryutt[258] - 
   aryutt[36] - aryutt[38] + aryutt[117] - aryutt[37] + aryutt[127] + 
   157./12. - aryutt[40];
   aryutt[265]=aryutt[193] - 9095./54.*aryutt[99];
   aryutt[265]=aryutt[15]*aryutt[265];
   aryutt[262]=1./2.*aryutt[265] + 1./2.*aryutt[262] + 235./9.*
   aryutt[206];
   aryutt[265]=aryutt[203] - 10./9.*aryutt[14] + 9./8.*aryutt[16];
   aryutt[265]=aryutt[3]*aryutt[265];
   aryutt[262]=755./216.*aryutt[216] + 1./2.*aryutt[262] + aryutt[265];
   aryutt[262]=aryutt[17]*aryutt[262];
   aryutt[266]= - 5935./648. + aryutt[144];
   aryutt[266]=aryutt[118] + 1./2.*aryutt[266] + aryutt[187];
   aryutt[267]=aryutt[193] - 3325./27.*aryutt[99];
   aryutt[267]=aryutt[15]*aryutt[267];
   aryutt[274]=1./12.*aryutt[8];
   aryutt[275]=1./4.*aryutt[205];
   aryutt[266]=1./2.*aryutt[267] + aryutt[275] - 739./324.*aryutt[12]
    + aryutt[274] + 169./288.*aryutt[21] + aryutt[217] + 1./2.*
   aryutt[266] + aryutt[38];
   aryutt[266]=aryutt[15]*aryutt[266];
   aryutt[267]=11./16. + aryutt[44];
   aryutt[230]=aryutt[12]*aryutt[230];
   aryutt[167]=1./18.*aryutt[230] + aryutt[204] + 1./4.*aryutt[269] - 
   11./48.*aryutt[84] + 1./3.*aryutt[267] + aryutt[167];
   aryutt[167]=MMH*aryutt[167];
   aryutt[204]=1./3.*aryutt[95] - 1./4.*aryutt[49];
   aryutt[230]=19./24.*aryutt[19];
   aryutt[267]=1./72.*aryutt[183] + 1./6.*aryutt[152] - 1./12.*
   aryutt[14] - 25./9.*aryutt[20] - 19./72.*aryutt[18] + aryutt[230] + 
   1./3.*aryutt[33] + 25./72.*aryutt[32] + 1./12.*aryutt[48] + 
   aryutt[204] - 13./36.*aryutt[92];
   aryutt[269]=1./12.*aryutt[21];
   aryutt[276]=5./16.*aryutt[8];
   aryutt[277]=aryutt[276] + aryutt[269] - 5./16. + aryutt[272];
   aryutt[278]=aryutt[277] - 1./54.*aryutt[12];
   aryutt[278]=aryutt[16]*aryutt[278];
   aryutt[237]= - 97./4.*aryutt[8] + aryutt[240] - 67./4. + aryutt[237]
   ;
   aryutt[237]=1./12.*aryutt[237] + aryutt[12];
   aryutt[237]=aryutt[17]*aryutt[237];
   aryutt[240]=5./108.*aryutt[12] - 5./12.*aryutt[8] + 5./27.*
   aryutt[36] + 1./2. - 17./27.*aryutt[38];
   aryutt[240]=aryutt[15]*aryutt[240];
   aryutt[167]=1./3.*aryutt[167] + 1./3.*aryutt[237] + 1./4.*
   aryutt[240] + 1./2.*aryutt[267] + aryutt[278];
   aryutt[167]=MMH*aryutt[167];
   aryutt[237]= - aryutt[14] + 121./9.*aryutt[16];
   aryutt[237]=aryutt[16]*aryutt[237];
   aryutt[240]=8893./36.*aryutt[15] + aryutt[137] - 209./3.*aryutt[16];
   aryutt[240]=aryutt[15]*aryutt[240];
   aryutt[237]=1./4.*aryutt[237] + 1./3.*aryutt[240];
   aryutt[240]=5*aryutt[14];
   aryutt[267]=aryutt[240] - 5971./6.*aryutt[16];
   aryutt[267]=101./6.*aryutt[17] + 1./4.*aryutt[267] + 3335./9.*
   aryutt[15];
   aryutt[267]=aryutt[17]*aryutt[267];
   aryutt[237]=1./2.*aryutt[237] + 1./3.*aryutt[267];
   aryutt[237]=1./6.*aryutt[237] + aryutt[167];
   aryutt[267]=1./3.*aryutt[6]*MMH*aryutt[255];
   aryutt[237]=1./2.*aryutt[237] + aryutt[267];
   aryutt[237]=aryutt[6]*aryutt[237];
   aryutt[278]=745./24. - aryutt[40];
   aryutt[117]= - aryutt[36] - aryutt[38] + aryutt[117] + aryutt[129]
    + 1./2.*aryutt[278] + aryutt[127];
   aryutt[117]=aryutt[14]*aryutt[117];
   aryutt[278]= - 3*aryutt[14];
   aryutt[280]=aryutt[278] - 5./4.*aryutt[16];
   aryutt[280]=aryutt[16]*aryutt[280];
   aryutt[236]= - 4./9.*aryutt[236];
   aryutt[281]=9./2.*aryutt[15] + aryutt[278] + 17./4.*aryutt[16];
   aryutt[281]=1./4.*aryutt[15]*aryutt[281];
   aryutt[280]=aryutt[281] + aryutt[236] + 1./4.*aryutt[280];
   aryutt[280]=aryutt[3]*aryutt[280];
   aryutt[282]= - 1./24.*aryutt[95] + aryutt[49];
   aryutt[148]=aryutt[148] - 13./12.*aryutt[92] + aryutt[268] + 139./48.
   *aryutt[96] + 7./6.*aryutt[94] + aryutt[282] - 11./16.*aryutt[50];
   aryutt[148]= - 763./288.*aryutt[19] - 83./144.*aryutt[42] - 1./24.*
   aryutt[5] - 31./96.*aryutt[33] - 341./576.*aryutt[34] + 1./2.*
   aryutt[148] - 8./9.*aryutt[32];
   aryutt[268]=2137./576.*aryutt[20];
   aryutt[283]=3./64.*aryutt[152];
   aryutt[228]=2./9.*aryutt[228];
   aryutt[284]=17./48.*aryutt[156];
   aryutt[106]=aryutt[106] + aryutt[112] + aryutt[237] + aryutt[249] + 
   aryutt[262] + aryutt[280] + 1./4.*aryutt[266] + 235./36.*aryutt[163]
    + 1./2.*aryutt[169] + aryutt[284] + aryutt[228] + aryutt[283] + 1./
   4.*aryutt[117] + aryutt[268] + aryutt[148] - 1649./1152.*aryutt[18];
   aryutt[106]=aryutt[176]*aryutt[106];
   aryutt[112]= - 9./2.*aryutt[73] + 101./96. - 13*aryutt[83];
   aryutt[112]= - 35./16.*aryutt[87] - 1./4.*aryutt[72] - 1./4.*
   aryutt[74] + aryutt[116] + aryutt[77] - 39./16.*aryutt[90] + 1./4.*
   aryutt[141] - 15./16.*aryutt[89] - 1./2.*aryutt[80] + 5./32.*
   aryutt[76] + 1./4.*aryutt[112] + aryutt[132];
   aryutt[102]= - 55*aryutt[41] - 617./12. + aryutt[102];
   aryutt[102]=aryutt[219] + aryutt[246] - aryutt[36] + 1./4.*
   aryutt[102] - aryutt[38];
   aryutt[102]=aryutt[21]*aryutt[102];
   aryutt[116]=55./2.*aryutt[39];
   aryutt[117]= - 55./2.*aryutt[41];
   aryutt[132]=aryutt[117] + aryutt[129] + aryutt[116] - 51./4. - 
   aryutt[40];
   aryutt[132]=aryutt[138] - aryutt[36] + 1./2.*aryutt[132] - 
   aryutt[38];
   aryutt[132]=aryutt[8]*aryutt[132];
   aryutt[169]=1./2.*aryutt[183] + 1./2.*aryutt[181] + 3*aryutt[42] - 
   aryutt[14];
   aryutt[149]=3./2.*aryutt[155] + 1./2.*aryutt[169] + aryutt[149];
   aryutt[149]=aryutt[3]*aryutt[149];
   aryutt[155]=7 + 17./4.*aryutt[21];
   aryutt[139]=aryutt[139] + 1./8.*aryutt[155] + aryutt[8];
   aryutt[139]=aryutt[12]*aryutt[139];
   aryutt[155]=1./2.*aryutt[25];
   aryutt[169]= - aryutt[64] - 3*aryutt[61] - 1./2.*aryutt[63];
   aryutt[169]=1./2.*aryutt[169] - aryutt[60];
   aryutt[169]=MMH*aryutt[169];
   aryutt[181]=1./2.*aryutt[68];
   aryutt[183]=aryutt[181] - aryutt[65];
   aryutt[183]=aryutt[63] + 1./2.*aryutt[183] + aryutt[61];
   aryutt[183]=MMt*aryutt[183];
   aryutt[102]=1./2.*aryutt[183] + 1./2.*aryutt[169] + aryutt[149] + 
   aryutt[195] + 1./2.*aryutt[161] + 1./2.*aryutt[139] + aryutt[132] + 
   1./4.*aryutt[102] + 1./8.*aryutt[166] + aryutt[212] + 1./4.*
   aryutt[115] + 3./64.*aryutt[13] + aryutt[157] + aryutt[155] + 13./8.
   *aryutt[23] + aryutt[123] - 21./16.*aryutt[51] + 1./8.*aryutt[88] + 
   3./8.*aryutt[85] + 1./2.*aryutt[112] + aryutt[44];
   aryutt[102]=MMt*aryutt[102];
   aryutt[112]= - 85./24. + aryutt[238];
   aryutt[112]= - 9./8.*aryutt[72] - 9./8.*aryutt[74] - 3./2.*
   aryutt[77] + aryutt[197] + aryutt[243] + 1./4.*aryutt[112] + 
   aryutt[147];
   aryutt[115]=1./4.*aryutt[23];
   aryutt[112]=aryutt[118] + 5./16.*aryutt[13] + aryutt[187] + 
   aryutt[115] - 23./24.*aryutt[84] + 51./16.*aryutt[51] + aryutt[119]
    + 3./8.*aryutt[88] - 5./12.*aryutt[85] + 1./2.*aryutt[112] - 
   aryutt[44];
   aryutt[119]=aryutt[244] + aryutt[107] + aryutt[235] + 55./8.*
   aryutt[41] + aryutt[199] - 55./8.*aryutt[39] + 55./16. + aryutt[40];
   aryutt[119]=aryutt[8]*aryutt[119];
   aryutt[123]=1./2.*aryutt[140] + aryutt[307] + aryutt[173] - 7*
   aryutt[84] - 1./2.*aryutt[79] - 93./16. + aryutt[88];
   aryutt[123]=MMH*aryutt[99]*aryutt[123];
   aryutt[132]= - 13./8.*aryutt[8];
   aryutt[139]=1./3.*aryutt[174] + aryutt[132];
   aryutt[139]=aryutt[12]*aryutt[139];
   aryutt[134]= - aryutt[18] + aryutt[134] + 3./4.*aryutt[32];
   aryutt[134]=7./4.*aryutt[14] + 1./2.*aryutt[134] + aryutt[20];
   aryutt[134]=aryutt[99]*aryutt[134];
   aryutt[112]=1./16.*aryutt[123] + 1./2.*aryutt[182] + 1./4.*
   aryutt[171] + 1./2.*aryutt[134] + 1./8.*aryutt[139] + 1./2.*
   aryutt[119] + 1./12.*aryutt[126] + 5./12.*aryutt[36] + 1./4.*
   aryutt[112] + aryutt[272];
   aryutt[112]=MMH*aryutt[112];
   aryutt[119]= - 199./3. + aryutt[144];
   aryutt[123]= - 55./4.*aryutt[41];
   aryutt[119]= - 71./8.*aryutt[12] + aryutt[227] - 133./24.*aryutt[21]
    + aryutt[120] - aryutt[36] + aryutt[226] + aryutt[123] + 
   aryutt[118] + aryutt[220] + 1./2.*aryutt[119] + aryutt[187];
   aryutt[119]=1./4.*aryutt[119] + aryutt[234];
   aryutt[119]=aryutt[16]*aryutt[119];
   aryutt[126]= - aryutt[101] + 9./2.*aryutt[99];
   aryutt[126]=aryutt[15]*aryutt[126];
   aryutt[123]=3./2.*aryutt[126] + 33*aryutt[206] + aryutt[261] + 7./8.
   *aryutt[12] + aryutt[159] + 17./32.*aryutt[21] - aryutt[36] - 
   aryutt[38] + aryutt[123] + aryutt[129] + aryutt[220] - 85./8. - 
   aryutt[40];
   aryutt[126]=aryutt[203] - aryutt[14] + 9./4.*aryutt[16];
   aryutt[126]=aryutt[3]*aryutt[126];
   aryutt[134]=9./8.*aryutt[216];
   aryutt[123]=aryutt[134] + 1./2.*aryutt[123] + aryutt[126];
   aryutt[123]=aryutt[17]*aryutt[123];
   aryutt[126]=aryutt[8]*aryutt[207];
   aryutt[139]=aryutt[12]*aryutt[288];
   aryutt[126]=1./6.*aryutt[139] + 1./2.*aryutt[126] + aryutt[214] + 1./
   6.*aryutt[304];
   aryutt[126]=MMH*aryutt[126];
   aryutt[140]=aryutt[223] + aryutt[214] + aryutt[202];
   aryutt[140]=aryutt[16]*aryutt[140];
   aryutt[147]=aryutt[184] + aryutt[214] + aryutt[244];
   aryutt[147]=aryutt[15]*aryutt[147];
   aryutt[139]=aryutt[304] + 1./2.*aryutt[139];
   aryutt[139]=MMH*aryutt[139];
   aryutt[139]=1./4.*aryutt[139] + 1./2.*aryutt[294] + aryutt[140] + 1./
   2.*aryutt[147];
   aryutt[139]=MMH*aryutt[139];
   aryutt[140]=aryutt[16] + aryutt[15];
   aryutt[147]=aryutt[15]*aryutt[140];
   aryutt[147]=1./2.*aryutt[296] - aryutt[186] + 1./2.*aryutt[147];
   aryutt[139]=17./4.*aryutt[147] + aryutt[139];
   aryutt[139]=aryutt[6]*aryutt[139];
   aryutt[147]=aryutt[179] - 7 + aryutt[164];
   aryutt[147]=aryutt[36] + 1./4.*aryutt[147] - aryutt[38];
   aryutt[147]=aryutt[211] + 3*aryutt[147] + 17./12.*aryutt[21];
   aryutt[147]=aryutt[16]*aryutt[147];
   aryutt[149]=aryutt[3]*aryutt[297];
   aryutt[157]=aryutt[14]*aryutt[198];
   aryutt[161]=aryutt[209] - 17./24.*aryutt[21] + aryutt[36] + 21./8.
    - aryutt[38];
   aryutt[161]=aryutt[15]*aryutt[161];
   aryutt[126]=1./3.*aryutt[139] + aryutt[126] + aryutt[218] + 17./2.*
   aryutt[149] + aryutt[161] + aryutt[157] + 1./2.*aryutt[147];
   aryutt[139]=aryutt[8]*aryutt[198];
   aryutt[139]=1./4.*aryutt[208] + aryutt[139];
   aryutt[139]=MMt*aryutt[139];
   aryutt[126]=1./2.*aryutt[126] + aryutt[139];
   aryutt[126]=aryutt[210]*aryutt[126];
   aryutt[139]=aryutt[221] + 1./4.*aryutt[152] - 5./8.*aryutt[14] - 1./
   6.*aryutt[20] + 19./48.*aryutt[18] + aryutt[230] - 3./16.*aryutt[32]
    - 1./8.*aryutt[48] + aryutt[204] + 13./24.*aryutt[92];
   aryutt[147]=aryutt[277] + 1./48.*aryutt[12];
   aryutt[147]=aryutt[16]*aryutt[147];
   aryutt[149]=1./12.*aryutt[12] + aryutt[276] + aryutt[269] - 7./8. + 
   aryutt[263];
   aryutt[149]=aryutt[15]*aryutt[149];
   aryutt[111]=aryutt[111] + 5./4. + 1./3.*aryutt[85];
   aryutt[157]= - aryutt[38] + aryutt[307];
   aryutt[157]=aryutt[21]*aryutt[157];
   aryutt[161]=aryutt[12]*aryutt[180];
   aryutt[111]=1./6.*aryutt[161] + 1./2.*aryutt[111] + 1./3.*
   aryutt[157];
   aryutt[111]=MMH*aryutt[111];
   aryutt[157]= - 13./16. - aryutt[36];
   aryutt[157]=1./3.*aryutt[157] - 5./16.*aryutt[8];
   aryutt[157]=aryutt[17]*aryutt[157];
   aryutt[111]=1./4.*aryutt[111] + 1./2.*aryutt[157] + 1./2.*
   aryutt[149] + 1./2.*aryutt[139] + aryutt[147];
   aryutt[111]=MMH*aryutt[111];
   aryutt[139]= - aryutt[14] - 197*aryutt[16];
   aryutt[139]=aryutt[16]*aryutt[139];
   aryutt[147]= - 85*aryutt[15] - aryutt[14] - 13*aryutt[16];
   aryutt[147]=aryutt[15]*aryutt[147];
   aryutt[139]=aryutt[139] + 1./2.*aryutt[147];
   aryutt[147]= - 5./3.*aryutt[14] - 423./2.*aryutt[16];
   aryutt[147]=5./2.*aryutt[17] + 1./4.*aryutt[147] + aryutt[15];
   aryutt[147]=aryutt[17]*aryutt[147];
   aryutt[139]=1./12.*aryutt[139] + aryutt[147];
   aryutt[111]=1./4.*aryutt[139] + aryutt[111];
   aryutt[111]=aryutt[6]*aryutt[111];
   aryutt[139]= - 3113./24. + aryutt[144];
   aryutt[139]=aryutt[118] + 1./2.*aryutt[139] + aryutt[187];
   aryutt[147]= - aryutt[101] + 3*aryutt[99];
   aryutt[147]=aryutt[15]*aryutt[147];
   aryutt[132]=3./4.*aryutt[147] + 3./8.*aryutt[205] + 23./24.*
   aryutt[12] + aryutt[132] - 205./192.*aryutt[21] + 1./4.*aryutt[139]
    - aryutt[36];
   aryutt[132]=aryutt[15]*aryutt[132];
   aryutt[116]=aryutt[117] + aryutt[129] + aryutt[116] - 305./24. - 
   aryutt[40];
   aryutt[116]= - aryutt[36] + 1./2.*aryutt[116] - aryutt[38];
   aryutt[116]=aryutt[14]*aryutt[116];
   aryutt[117]= - aryutt[14] - 23./4.*aryutt[16];
   aryutt[117]=aryutt[16]*aryutt[117];
   aryutt[139]=25./2.*aryutt[15] + aryutt[278] + 13./4.*aryutt[16];
   aryutt[139]=aryutt[15]*aryutt[139];
   aryutt[117]=3*aryutt[117] + 1./2.*aryutt[139];
   aryutt[117]=aryutt[3]*aryutt[117];
   aryutt[102]=aryutt[126] + aryutt[102] + aryutt[111] + aryutt[112] + 
   aryutt[123] + 1./2.*aryutt[117] + 1./2.*aryutt[132] + 33./4.*
   aryutt[163] + aryutt[119] + 3./16.*aryutt[156] + 3./32.*aryutt[152]
    + 1./2.*aryutt[116] + 27./64.*aryutt[20] - 1901./384.*aryutt[18] - 
   169./32.*aryutt[19] - 149./96.*aryutt[42] - 1./8.*aryutt[5] + 3./32.
   *aryutt[33] - 23./64.*aryutt[34] + aryutt[125] + 3./8.*aryutt[48] - 
   1./8.*aryutt[92] - 1./8.*aryutt[93] + 247./32.*aryutt[96] + 
   aryutt[94] + aryutt[282] - 15./8.*aryutt[50];
   aryutt[102]=aryutt[210]*aryutt[102];
   aryutt[111]=367./9. - aryutt[40];
   aryutt[112]=47./12.*aryutt[41];
   aryutt[111]=aryutt[138] - aryutt[36] - aryutt[38] + aryutt[112] + 
   aryutt[129] + 1./2.*aryutt[111] + aryutt[127];
   aryutt[111]=aryutt[8]*aryutt[111];
   aryutt[116]=47./4.*aryutt[41] + 473./12. + aryutt[143];
   aryutt[116]=aryutt[185] + aryutt[246] - aryutt[36] + 1./3.*
   aryutt[116] - aryutt[38];
   aryutt[116]=aryutt[21]*aryutt[116];
   aryutt[110]=aryutt[172] + aryutt[225] + aryutt[160] + aryutt[195] + 
   aryutt[200] + aryutt[191] + 1./2.*aryutt[111] + aryutt[110] + 1./8.*
   aryutt[116];
   aryutt[110]=MMt*aryutt[110];
   aryutt[111]=3635./54. + aryutt[144];
   aryutt[111]= - 5149./108.*aryutt[12] + aryutt[227] + 77./8.*
   aryutt[21] + aryutt[120] - aryutt[36] + aryutt[226] + aryutt[112] + 
   aryutt[118] + aryutt[127] + 1./2.*aryutt[111] + aryutt[187];
   aryutt[111]=1./4.*aryutt[111] + aryutt[234];
   aryutt[111]=aryutt[16]*aryutt[111];
   aryutt[107]=aryutt[244] + aryutt[107] + aryutt[235] - 47./24.*
   aryutt[41] + aryutt[37] + aryutt[239] - 437./36. + aryutt[40];
   aryutt[107]=aryutt[8]*aryutt[107];
   aryutt[107]=aryutt[222] + aryutt[257] + aryutt[256] + aryutt[252] + 
   aryutt[254] + aryutt[146] + 1./4.*aryutt[107];
   aryutt[107]=MMH*aryutt[107];
   aryutt[116]=aryutt[261] + aryutt[260] + aryutt[259] + aryutt[258] - 
   aryutt[36] - aryutt[38] + aryutt[112] - aryutt[37] + aryutt[127] + 
   205./12. - aryutt[40];
   aryutt[117]= - aryutt[99]*aryutt[16];
   aryutt[119]=aryutt[193] - 271./6.*aryutt[99];
   aryutt[119]=aryutt[15]*aryutt[119];
   aryutt[116]=1./2.*aryutt[119] + 1./2.*aryutt[116] + 77*aryutt[117];
   aryutt[116]=aryutt[134] + 1./2.*aryutt[116] + aryutt[265];
   aryutt[116]=aryutt[17]*aryutt[116];
   aryutt[117]= - 6101./216. + aryutt[144];
   aryutt[117]=aryutt[118] + 1./2.*aryutt[117] + aryutt[187];
   aryutt[118]=aryutt[193] - 85./3.*aryutt[99];
   aryutt[118]=aryutt[15]*aryutt[118];
   aryutt[117]=1./2.*aryutt[118] + aryutt[275] + 95./108.*aryutt[12] + 
   aryutt[274] + 211./96.*aryutt[21] + aryutt[217] + 1./2.*aryutt[117]
    + aryutt[38];
   aryutt[117]=aryutt[15]*aryutt[117];
   aryutt[118]= - aryutt[14] + 563./3.*aryutt[16];
   aryutt[118]=aryutt[16]*aryutt[118];
   aryutt[119]=575./12.*aryutt[15] + aryutt[137] - 737./3.*aryutt[16];
   aryutt[119]=aryutt[15]*aryutt[119];
   aryutt[118]=1./4.*aryutt[118] + 1./3.*aryutt[119];
   aryutt[119]=aryutt[240] - 211./6.*aryutt[16];
   aryutt[119]=1./4.*aryutt[119] + 269./3.*aryutt[15];
   aryutt[120]= - 1./2.*aryutt[17];
   aryutt[119]=1./9.*aryutt[119] + aryutt[120];
   aryutt[119]=aryutt[17]*aryutt[119];
   aryutt[118]=1./6.*aryutt[118] + aryutt[119];
   aryutt[118]=1./2.*aryutt[118] + aryutt[167];
   aryutt[118]=1./2.*aryutt[118] + aryutt[267];
   aryutt[118]=aryutt[6]*aryutt[118];
   aryutt[119]=937./24. - aryutt[40];
   aryutt[112]= - aryutt[36] - aryutt[38] + aryutt[112] + aryutt[129]
    + 1./2.*aryutt[119] + aryutt[127];
   aryutt[112]=aryutt[14]*aryutt[112];
   aryutt[119]=aryutt[278] + 91./4.*aryutt[16];
   aryutt[119]=aryutt[16]*aryutt[119];
   aryutt[119]=aryutt[281] + aryutt[236] + 1./4.*aryutt[119];
   aryutt[119]=aryutt[3]*aryutt[119];
   aryutt[123]= - aryutt[99]*aryutt[186];
   aryutt[102]=1./2.*aryutt[102] + aryutt[110] + aryutt[118] + 
   aryutt[107] + aryutt[116] + aryutt[119] + 1./4.*aryutt[117] + 77./4.
   *aryutt[123] + 1./2.*aryutt[111] + aryutt[284] + aryutt[228] + 
   aryutt[283] + 1./4.*aryutt[112] + aryutt[268] + aryutt[148] - 161./
   1152.*aryutt[18];
   aryutt[102]=aryutt[210]*aryutt[102];
   aryutt[107]=aryutt[15]*aryutt[99];
   aryutt[110]=aryutt[17]*aryutt[99];
   aryutt[110]=8./9.*aryutt[110] + aryutt[206] + 8./3.*aryutt[107];
   aryutt[110]=aryutt[17]*aryutt[110];
   aryutt[111]= - 1 + aryutt[12];
   aryutt[112]=1./3.*aryutt[111] + 4*aryutt[107];
   aryutt[112]=aryutt[15]*aryutt[112];
   aryutt[116]=aryutt[16]*aryutt[111];
   aryutt[116]=2./9.*aryutt[18] + 5*aryutt[116];
   aryutt[110]=4*aryutt[110] + 8./9.*aryutt[112] + 1./3.*aryutt[116] + 
   2*aryutt[163];
   aryutt[112]=aryutt[124] - 16./3.*aryutt[15];
   aryutt[112]=aryutt[15]*aryutt[112];
   aryutt[116]= - 4./3.*aryutt[17] + aryutt[16] - 11./3.*aryutt[15];
   aryutt[116]=aryutt[17]*aryutt[116];
   aryutt[112]=4*aryutt[116] + aryutt[186] + aryutt[112];
   aryutt[112]=aryutt[6]*aryutt[112];
   aryutt[110]=2*aryutt[110] + 1./3.*aryutt[112];
   aryutt[102]=aryutt[105] + aryutt[102] + 8./3.*aryutt[110] + 
   aryutt[106];
   aryutt[102]=aryutt[9]*aryutt[102];
   aryutt[105]= - 1./2.*aryutt[70];
   aryutt[106]= - 2./3.*aryutt[52];
   aryutt[110]= - 2*aryutt[69];
   aryutt[112]=2./3.*aryutt[67] + 2*aryutt[65] + aryutt[181] + 
   aryutt[110] + aryutt[105] + aryutt[106];
   aryutt[116]=pow(CW,2);
   aryutt[112]=aryutt[116]*aryutt[112];
   aryutt[105]=64./27.*aryutt[54] + aryutt[110] + aryutt[106] - 128./27.
   *aryutt[57] + aryutt[105];
   aryutt[106]=aryutt[130] - 4./9.*aryutt[12];
   aryutt[106]=aryutt[3]*aryutt[106];
   aryutt[105]=aryutt[106] + aryutt[175] + aryutt[168] + 1./3.*
   aryutt[112] + 17./18.*aryutt[67] + 5./2.*aryutt[65] - 128./81.*
   aryutt[66] + 1./3.*aryutt[105] + 3./4.*aryutt[68];
   aryutt[105]=MMZ*aryutt[105];
   aryutt[106]=20*aryutt[116];
   aryutt[110]=8*aryutt[116];
   aryutt[112]=49./3. + aryutt[110];
   aryutt[112]=aryutt[41]*aryutt[112];
   aryutt[117]= - 2*aryutt[21];
   aryutt[118]= - 32./9.*aryutt[12] + aryutt[117] + 2*aryutt[112] + 313.
   /9. + aryutt[106];
   aryutt[118]=aryutt[12]*aryutt[118];
   aryutt[119]=4./3.*aryutt[116] + aryutt[90] + 2./3.*aryutt[80] - 38./
   9. - aryutt[81];
   aryutt[119]=aryutt[116]*aryutt[119];
   aryutt[123]=53./12. + aryutt[116];
   aryutt[123]=aryutt[41]*aryutt[123];
   aryutt[123]=aryutt[123] + aryutt[39] + 22./3. + aryutt[116];
   aryutt[123]=aryutt[21]*aryutt[123];
   aryutt[124]= - 701./27. - aryutt[81];
   aryutt[124]=2*aryutt[124] + 25./6.*aryutt[58];
   aryutt[125]= - 2./9.*aryutt[24];
   aryutt[126]= - 1 - 2./3.*aryutt[116];
   aryutt[126]=aryutt[13]*aryutt[126];
   aryutt[127]= - aryutt[3]*aryutt[15];
   aryutt[132]= - aryutt[17]*aryutt[3];
   aryutt[105]=aryutt[105] + 8./9.*aryutt[132] + 4./9.*aryutt[127] + 2./
   9.*aryutt[118] + aryutt[123] + 4./3.*aryutt[126] + 4./3.*aryutt[119]
    + aryutt[125] + 146./81.*aryutt[25] + 5./4.*aryutt[23] + 64./81.*
   aryutt[45] + 101./36.*aryutt[87] + 64./81.*aryutt[75] + 3*aryutt[90]
    - 29./24.*aryutt[141] - 28./27.*aryutt[80] + 2./9.*aryutt[76] + 1./
   3.*aryutt[124] - 5./4.*aryutt[83];
   aryutt[105]=MMZ*aryutt[105];
   aryutt[118]= - 203./6. - aryutt[58];
   aryutt[118]=139./384.*aryutt[141] + 1./3.*aryutt[118] + aryutt[76];
   aryutt[119]=29./4.*aryutt[141] + 2 - aryutt[58];
   aryutt[119]=1./3.*aryutt[119] + aryutt[90];
   aryutt[119]=aryutt[116]*aryutt[119];
   aryutt[123]= - 1 - aryutt[116];
   aryutt[123]=1./9.*aryutt[123] + aryutt[177];
   aryutt[123]=aryutt[21]*aryutt[123];
   aryutt[118]= - 16./27.*aryutt[12] + aryutt[123] - 25./27.*aryutt[13]
    + 1./3.*aryutt[119] - 16./27.*aryutt[45] + 2./3.*aryutt[87] + 16./
   27.*aryutt[75] + 1./3.*aryutt[118] + 1./2.*aryutt[90];
   aryutt[118]=MMZ*aryutt[118];
   aryutt[119]=104./27. - aryutt[116];
   aryutt[119]=aryutt[20]*aryutt[119];
   aryutt[119]=2*aryutt[119] + 40./9.*aryutt[18] - 16./9.*aryutt[33] - 
   2./3.*aryutt[34] + aryutt[96] - 16./9.*aryutt[93];
   aryutt[123]=85 + 16*aryutt[12];
   aryutt[123]=aryutt[15]*aryutt[123];
   aryutt[124]= - 16./3.*aryutt[12];
   aryutt[126]=aryutt[124] - 26./3. + aryutt[21];
   aryutt[126]=aryutt[17]*aryutt[126];
   aryutt[119]=2./3.*aryutt[126] + 1./9.*aryutt[123] + 2*aryutt[119] - 
   5./3.*aryutt[16];
   aryutt[118]=1./3.*aryutt[119] + aryutt[118];
   aryutt[118]=MMZ*aryutt[118];
   aryutt[119]=aryutt[116]*aryutt[141];
   aryutt[119]=aryutt[141] + 1./2.*aryutt[119];
   aryutt[119]=aryutt[245]*aryutt[119];
   aryutt[119]= - 32./27.*aryutt[255] + 21./16.*aryutt[119];
   aryutt[119]=MMZ*aryutt[119];
   aryutt[123]= - aryutt[116]*aryutt[141];
   aryutt[123]= - aryutt[141] + 1./2.*aryutt[123];
   aryutt[123]=aryutt[116]*aryutt[123];
   aryutt[123]= - 3./2.*aryutt[141] + aryutt[123];
   aryutt[123]=aryutt[6]*aryutt[196]*aryutt[123];
   aryutt[119]=aryutt[119] + 5./16.*aryutt[123];
   aryutt[119]=aryutt[6]*aryutt[119];
   aryutt[123]= - 8./9.*aryutt[17] + aryutt[16] - 80./27.*aryutt[15];
   aryutt[123]=aryutt[17]*aryutt[123];
   aryutt[123]= - 8./9.*aryutt[229] + aryutt[123];
   aryutt[118]=aryutt[119] + 4./3.*aryutt[123] + aryutt[118];
   aryutt[118]=aryutt[6]*aryutt[118];
   aryutt[119]= - 23./3.*aryutt[50] - 2*aryutt[94];
   aryutt[123]=aryutt[116]*aryutt[50];
   aryutt[126]= - 4*aryutt[116];
   aryutt[127]=29./3. + aryutt[126];
   aryutt[127]=aryutt[19]*aryutt[127];
   aryutt[119]=2./3.*aryutt[127] + 4./3.*aryutt[123] - 128./27.*
   aryutt[33] - 2*aryutt[34] - 88./9.*aryutt[93] + 1./3.*aryutt[119] + 
   4*aryutt[96];
   aryutt[112]=160./9.*aryutt[12] + 5./4.*aryutt[21] + 4*aryutt[112] + 
   2009./12. + 40*aryutt[116];
   aryutt[112]=aryutt[15]*aryutt[112];
   aryutt[127]= - 16*aryutt[116];
   aryutt[134]= - 49./3. - 8*aryutt[116];
   aryutt[134]=aryutt[41]*aryutt[134];
   aryutt[134]=2*aryutt[134] - 265./3. + aryutt[127];
   aryutt[134]=32./27.*aryutt[12] + 1./3.*aryutt[134] + 2*aryutt[21];
   aryutt[134]=aryutt[17]*aryutt[134];
   aryutt[137]=69./4. - 8./9.*aryutt[116];
   aryutt[137]=aryutt[18]*aryutt[137];
   aryutt[138]=217./9. + aryutt[127];
   aryutt[138]=aryutt[20]*aryutt[138];
   aryutt[139]= - 37./3. + aryutt[127];
   aryutt[139]=40./9.*aryutt[12] + aryutt[21] + 1./9.*aryutt[139] + 
   aryutt[41];
   aryutt[139]=aryutt[16]*aryutt[139];
   aryutt[105]=aryutt[118] + aryutt[105] + 2./3.*aryutt[134] + 1./9.*
   aryutt[112] + aryutt[139] + 4./9.*aryutt[138] + 2./3.*aryutt[119] + 
   aryutt[137];
   aryutt[105]=aryutt[6]*aryutt[105];
   aryutt[112]=4*aryutt[116] - aryutt[90] + 23./3. + aryutt[80];
   aryutt[112]=aryutt[116]*aryutt[112];
   aryutt[118]= - 64./3.*aryutt[75] - 10*aryutt[90] + 104./3.*
   aryutt[80] + 4211./18. - 4*aryutt[76];
   aryutt[119]=11./9. + aryutt[126];
   aryutt[119]=aryutt[13]*aryutt[119];
   aryutt[134]=13./3. + 2*aryutt[116];
   aryutt[134]=aryutt[41]*aryutt[134];
   aryutt[137]= - 2*aryutt[116];
   aryutt[138]= - 13./3. + aryutt[137];
   aryutt[138]=aryutt[41]*aryutt[138];
   aryutt[139]= - 16./9. + aryutt[138];
   aryutt[139]=8./3.*aryutt[12] + 2*aryutt[139] + aryutt[21];
   aryutt[139]=aryutt[12]*aryutt[139];
   aryutt[112]=4./3.*aryutt[139] - 4./3.*aryutt[21] + 8./3.*aryutt[134]
    + 4./3.*aryutt[119] + 16./3.*aryutt[112] + 1./3.*aryutt[118] - 8*
   aryutt[87];
   aryutt[112]=aryutt[99]*aryutt[112];
   aryutt[113]=1./8. + aryutt[113];
   aryutt[113]=aryutt[151]*aryutt[113];
   aryutt[109]=aryutt[153]*aryutt[109];
   aryutt[109]=1./4.*aryutt[113] + 3*aryutt[109];
   aryutt[109]=MMZ*aryutt[109];
   aryutt[113]= - 3*aryutt[39] + 37./9. + aryutt[137];
   aryutt[118]= - 53./3. + aryutt[126];
   aryutt[118]=aryutt[41]*aryutt[118];
   aryutt[113]= - 8./9.*aryutt[12] + 2*aryutt[113] + aryutt[118];
   aryutt[113]=aryutt[3]*aryutt[113];
   aryutt[109]=aryutt[109] + aryutt[113] + aryutt[112] + aryutt[178] + 
   1./3.*aryutt[65] - 1./4.*aryutt[68] - 1./3.*aryutt[69] + aryutt[224]
    - 256./81.*aryutt[57] + aryutt[215];
   aryutt[109]=MMZ*aryutt[109];
   aryutt[112]= - 8./3. - aryutt[116];
   aryutt[112]=aryutt[19]*aryutt[112];
   aryutt[112]=4*aryutt[112] + 2*aryutt[123] + 16./9.*aryutt[33] + 
   aryutt[34] + 16./3.*aryutt[93] - 2*aryutt[96] + 13./3.*aryutt[50] + 
   aryutt[94];
   aryutt[113]= - 29 + aryutt[126];
   aryutt[113]=aryutt[18]*aryutt[113];
   aryutt[118]= - 29./3. + aryutt[126];
   aryutt[118]=aryutt[16]*aryutt[118];
   aryutt[112]=2*aryutt[118] - 260./9.*aryutt[20] + 2*aryutt[112] + 
   aryutt[113];
   aryutt[112]=aryutt[99]*aryutt[112];
   aryutt[110]=aryutt[117] + 4*aryutt[134] + 341./9. + aryutt[110];
   aryutt[110]=aryutt[17]*aryutt[99]*aryutt[110];
   aryutt[113]=4./3.*aryutt[80] + 1 - 4./3.*aryutt[81];
   aryutt[113]=aryutt[116]*aryutt[113];
   aryutt[117]= - 181./9. + aryutt[127];
   aryutt[117]=aryutt[13]*aryutt[117];
   aryutt[113]=2./3.*aryutt[117] + 8*aryutt[113] + 8./27.*aryutt[25] + 
   16./27.*aryutt[45] - 32./27.*aryutt[75] + 110./9.*aryutt[80] + 5 - 
   98./9.*aryutt[81];
   aryutt[117]= - 65./3. + aryutt[127];
   aryutt[117]=aryutt[41]*aryutt[117];
   aryutt[113]=2*aryutt[113] + 1./3.*aryutt[117];
   aryutt[116]=5./3. + aryutt[116];
   aryutt[116]=aryutt[41]*aryutt[116];
   aryutt[106]=16./9.*aryutt[12] + 4*aryutt[116] + 241./9. + 
   aryutt[106];
   aryutt[106]=aryutt[12]*aryutt[106];
   aryutt[116]= - 1 - aryutt[41];
   aryutt[116]=aryutt[21]*aryutt[116];
   aryutt[117]= - 16./9.*aryutt[12] + aryutt[21] - 259./9. + 2*
   aryutt[138];
   aryutt[117]=aryutt[15]*aryutt[99]*aryutt[117];
   aryutt[118]= - 5 - 3*aryutt[41];
   aryutt[118]=aryutt[3]*aryutt[16]*aryutt[118];
   aryutt[119]=aryutt[66] + aryutt[55] - 2*aryutt[54];
   aryutt[119]=MMt*aryutt[119];
   aryutt[102]=aryutt[102] + aryutt[103] + aryutt[104] + 128./81.*
   aryutt[119] + aryutt[105] + aryutt[109] + 8./3.*aryutt[110] + 2*
   aryutt[118] + 8./3.*aryutt[117] + 8./3.*aryutt[112] + 4./9.*
   aryutt[106] + 1./3.*aryutt[113] + 1./2.*aryutt[116];
   aryutt[102]=aryutt[46]*aryutt[102];
   aryutt[103]=3./2.*aryutt[20] + aryutt[162] + aryutt[19] + 
   aryutt[250] - aryutt[34] - 1./4.*aryutt[33];
   aryutt[104]=7./2.*aryutt[10];
   aryutt[105]= - 3*aryutt[30];
   aryutt[106]=aryutt[104] + aryutt[105] + 3 + 11./2.*aryutt[29];
   aryutt[106]=aryutt[15]*aryutt[106];
   aryutt[109]=aryutt[10] + 3./2. + 2*aryutt[29];
   aryutt[109]=aryutt[16]*aryutt[109];
   aryutt[103]=1./4.*aryutt[106] + 3*aryutt[103] + aryutt[109];
   aryutt[103]=aryutt[3]*aryutt[103];
   aryutt[106]= - 1 + aryutt[29];
   aryutt[109]=3./4.*aryutt[106] + aryutt[10];
   aryutt[109]=aryutt[12]*aryutt[109];
   aryutt[110]= - 3./4.*aryutt[29];
   aryutt[109]=aryutt[109] - aryutt[10] + aryutt[110] + aryutt[13] - 
   aryutt[24] + 3./2. - aryutt[25];
   aryutt[109]=aryutt[99]*aryutt[109];
   aryutt[112]=25./4.*aryutt[29];
   aryutt[113]=17./4.*aryutt[10];
   aryutt[116]=aryutt[113] - 11 + aryutt[112];
   aryutt[116]=aryutt[3]*aryutt[116];
   aryutt[109]=3./16.*aryutt[109] + aryutt[116];
   aryutt[109]=MMZ*aryutt[109];
   aryutt[116]=1 - aryutt[29];
   aryutt[117]=aryutt[232] + 3./4.*aryutt[116] - aryutt[10];
   aryutt[117]=aryutt[99]*aryutt[117];
   aryutt[118]=15./2.*aryutt[12] + 31 + 15*aryutt[21];
   aryutt[118]=aryutt[3]*aryutt[118];
   aryutt[117]=3./2.*aryutt[117] + aryutt[118];
   aryutt[117]=aryutt[17]*aryutt[117];
   aryutt[118]=1./8.*aryutt[13] - 1./8.*aryutt[24] + 9./8.*aryutt[25]
    + 1 - 3./8.*aryutt[45];
   aryutt[118]= - 9./8.*aryutt[10] + 5./2.*aryutt[30] + 3*aryutt[118]
    - 65./32.*aryutt[29];
   aryutt[119]= - 1./8.*aryutt[30];
   aryutt[123]=1./12.*aryutt[10] + aryutt[119] - 2./3. - 3./8.*
   aryutt[29];
   aryutt[123]=aryutt[21]*aryutt[123];
   aryutt[126]= - 239 - 13*aryutt[29];
   aryutt[126]=11./12.*aryutt[10] + 1./48.*aryutt[126] - aryutt[30];
   aryutt[126]=aryutt[12]*aryutt[126];
   aryutt[127]=1./4.*aryutt[33];
   aryutt[134]= - 5./2.*aryutt[20] + 3./4.*aryutt[18] + aryutt[127] - 
   aryutt[5];
   aryutt[134]=aryutt[99]*aryutt[134];
   aryutt[137]=aryutt[10] - 1 + 3./4.*aryutt[29];
   aryutt[137]=aryutt[15]*aryutt[99]*aryutt[137];
   aryutt[103]=1./2.*aryutt[109] + 1./4.*aryutt[117] + aryutt[103] + 3./
   16.*aryutt[137] + 3./16.*aryutt[134] + 1./8.*aryutt[126] + 1./4.*
   aryutt[118] + aryutt[123];
   aryutt[103]=aryutt[27]*aryutt[103];
   aryutt[109]=aryutt[173] - 1./2.*aryutt[24] + aryutt[155] + 3./8. - 
   aryutt[22];
   aryutt[117]= - 13./2.*aryutt[10] + 3*aryutt[109] + 5*aryutt[11];
   aryutt[117]=aryutt[31]*aryutt[117];
   aryutt[109]= - 13./6.*aryutt[10] + aryutt[109] + 5./3.*aryutt[11];
   aryutt[109]=aryutt[1]*aryutt[109];
   aryutt[109]=aryutt[117] + aryutt[109];
   aryutt[117]=aryutt[10] + 1./2. - aryutt[11];
   aryutt[118]=aryutt[31]*aryutt[117];
   aryutt[117]=aryutt[1]*aryutt[117];
   aryutt[117]=aryutt[118] + 1./3.*aryutt[117];
   aryutt[117]=aryutt[21]*aryutt[117];
   aryutt[118]=1 - aryutt[11];
   aryutt[123]=aryutt[118] + 5./2.*aryutt[10];
   aryutt[126]=aryutt[31]*aryutt[123];
   aryutt[123]=aryutt[1]*aryutt[123];
   aryutt[123]=aryutt[126] + 1./3.*aryutt[123];
   aryutt[123]=aryutt[12]*aryutt[123];
   aryutt[126]= - aryutt[5] + aryutt[18];
   aryutt[134]=1./2.*aryutt[126] - aryutt[20];
   aryutt[137]=aryutt[31]*aryutt[134];
   aryutt[134]=aryutt[1]*aryutt[134];
   aryutt[138]=3*aryutt[137] + aryutt[134];
   aryutt[138]=aryutt[99]*aryutt[138];
   aryutt[139]=aryutt[31]*aryutt[10];
   aryutt[141]=aryutt[1]*aryutt[10];
   aryutt[143]=3*aryutt[139] + aryutt[141];
   aryutt[144]=aryutt[15]*aryutt[99]*aryutt[143];
   aryutt[109]=1./2.*aryutt[144] + aryutt[138] + 1./2.*aryutt[123] + 1./
   2.*aryutt[109] + aryutt[117];
   aryutt[117]=aryutt[145] + aryutt[133] + aryutt[26] - aryutt[23];
   aryutt[123]=aryutt[31]*aryutt[117];
   aryutt[117]=aryutt[1]*aryutt[117];
   aryutt[133]= - 1 - 3*aryutt[10];
   aryutt[133]=aryutt[31]*aryutt[133];
   aryutt[138]= - 1./3. - aryutt[10];
   aryutt[144]=aryutt[1]*aryutt[138];
   aryutt[133]=aryutt[133] + aryutt[144];
   aryutt[133]=aryutt[21]*aryutt[133];
   aryutt[117]=aryutt[133] + 3*aryutt[123] + aryutt[117];
   aryutt[123]= - 1 + aryutt[11];
   aryutt[123]=1./4.*aryutt[123] - aryutt[10];
   aryutt[133]=aryutt[31]*aryutt[123];
   aryutt[123]=aryutt[1]*aryutt[123];
   aryutt[123]=aryutt[133] + 1./3.*aryutt[123];
   aryutt[133]=aryutt[12]*aryutt[123];
   aryutt[117]=1./2.*aryutt[117] + aryutt[133];
   aryutt[117]=MMZ*aryutt[117];
   aryutt[133]=1./2.*aryutt[20] + aryutt[293] - aryutt[19] + aryutt[7]
    + aryutt[233];
   aryutt[145]=aryutt[31]*aryutt[133];
   aryutt[133]=aryutt[1]*aryutt[133];
   aryutt[133]=3*aryutt[145] + aryutt[133];
   aryutt[145]= - 1 - aryutt[11];
   aryutt[145]=1./2.*aryutt[145] - aryutt[10];
   aryutt[146]=aryutt[31]*aryutt[145];
   aryutt[145]=aryutt[1]*aryutt[145];
   aryutt[145]=aryutt[146] + 1./3.*aryutt[145];
   aryutt[145]=aryutt[16]*aryutt[145];
   aryutt[123]=aryutt[15]*aryutt[123];
   aryutt[146]= - 7./2. - aryutt[11];
   aryutt[146]=1./4.*aryutt[146] + aryutt[10];
   aryutt[147]=aryutt[31]*aryutt[146];
   aryutt[146]=aryutt[1]*aryutt[146];
   aryutt[146]=aryutt[147] + 1./3.*aryutt[146];
   aryutt[146]=aryutt[17]*aryutt[146];
   aryutt[117]=aryutt[117] + aryutt[146] + aryutt[123] + 1./2.*
   aryutt[133] + aryutt[145];
   aryutt[123]= - 1 + 1./6.*aryutt[29];
   aryutt[133]=5./6.*aryutt[10];
   aryutt[145]= - 1./2.*aryutt[30];
   aryutt[123]=aryutt[248] - 5./3.*aryutt[21] + aryutt[133] + 7*
   aryutt[123] + aryutt[145];
   aryutt[123]=aryutt[17]*aryutt[123];
   aryutt[146]= - 1./3.*aryutt[29];
   aryutt[147]= - 1./6.*aryutt[10];
   aryutt[148]= - 1./4.*aryutt[30];
   aryutt[149]=aryutt[147] + aryutt[148] - 1./4. + aryutt[146];
   aryutt[149]=aryutt[16]*aryutt[149];
   aryutt[146]=1./2. + aryutt[146];
   aryutt[151]= - 5./3.*aryutt[10];
   aryutt[146]=aryutt[151] + 7*aryutt[146] + aryutt[30];
   aryutt[146]=aryutt[12]*aryutt[146];
   aryutt[128]=aryutt[128] + 1./4.*aryutt[24] - aryutt[23] + 1./2.*
   aryutt[45] + 1./4. + aryutt[43];
   aryutt[152]= - 7./6.*aryutt[10] - 1 - 11./6.*aryutt[29];
   aryutt[152]=aryutt[21]*aryutt[152];
   aryutt[128]=1./2.*aryutt[146] + 3*aryutt[128] + aryutt[152];
   aryutt[128]=MMZ*aryutt[128];
   aryutt[127]= - 1./4.*aryutt[20] - aryutt[19] + 1./4.*aryutt[5] + 
   aryutt[34] + aryutt[127];
   aryutt[146]=aryutt[151] + aryutt[30] - 1 - 7./3.*aryutt[29];
   aryutt[146]=aryutt[15]*aryutt[146];
   aryutt[123]=1./4.*aryutt[128] + 1./4.*aryutt[123] + 1./8.*
   aryutt[146] + 3./4.*aryutt[127] + aryutt[149];
   aryutt[123]=aryutt[27]*aryutt[123];
   aryutt[127]=1 - aryutt[12];
   aryutt[127]=aryutt[17]*aryutt[127];
   aryutt[128]= - aryutt[12] - 1 - aryutt[45];
   aryutt[128]=MMZ*aryutt[128];
   aryutt[146]=1./2.*aryutt[128] + aryutt[127] + aryutt[188] - 1./2.*
   aryutt[33] + aryutt[20];
   aryutt[146]=MMZ*aryutt[146];
   aryutt[149]= - aryutt[15] + 3./4.*aryutt[17];
   aryutt[149]=aryutt[17]*aryutt[149];
   aryutt[149]=aryutt[149] + 1./2.*aryutt[146];
   aryutt[149]=aryutt[6]*aryutt[27]*aryutt[149];
   aryutt[117]=3./8.*aryutt[149] + 1./2.*aryutt[117] + aryutt[123];
   aryutt[117]=aryutt[6]*aryutt[117];
   aryutt[123]=1./2.*aryutt[29] - aryutt[30];
   aryutt[149]=aryutt[123] + 1./2.*aryutt[10];
   aryutt[152]=aryutt[21]*aryutt[149];
   aryutt[155]=aryutt[12]*aryutt[149];
   aryutt[152]=aryutt[152] + 1./2.*aryutt[155];
   aryutt[152]=MMZ*aryutt[152];
   aryutt[155]=aryutt[16]*aryutt[149];
   aryutt[156]=aryutt[15]*aryutt[149];
   aryutt[157]= - 1./2.*aryutt[10];
   aryutt[160]=aryutt[157] - 1./2.*aryutt[29] + aryutt[30];
   aryutt[161]=aryutt[17]*aryutt[160];
   aryutt[152]=aryutt[152] + 1./2.*aryutt[161] + aryutt[155] + 1./2.*
   aryutt[156];
   aryutt[152]=aryutt[27]*aryutt[152];
   aryutt[155]= - aryutt[11] + aryutt[10];
   aryutt[156]=aryutt[31]*aryutt[155];
   aryutt[155]=aryutt[1]*aryutt[155];
   aryutt[161]=aryutt[156] + 1./3.*aryutt[155];
   aryutt[163]=aryutt[21]*aryutt[161];
   aryutt[164]=aryutt[12]*aryutt[161];
   aryutt[163]=aryutt[163] + 1./2.*aryutt[164];
   aryutt[163]=MMZ*aryutt[163];
   aryutt[164]=aryutt[16]*aryutt[161];
   aryutt[166]=aryutt[15]*aryutt[161];
   aryutt[167]=aryutt[11] - aryutt[10];
   aryutt[168]=aryutt[31]*aryutt[167];
   aryutt[167]=aryutt[1]*aryutt[167];
   aryutt[169]=aryutt[168] + 1./3.*aryutt[167];
   aryutt[171]=aryutt[17]*aryutt[169];
   aryutt[152]=aryutt[152] + aryutt[163] + 1./2.*aryutt[171] + 
   aryutt[164] + 1./2.*aryutt[166];
   aryutt[152]=aryutt[6]*aryutt[152];
   aryutt[163]=aryutt[12]*aryutt[169];
   aryutt[164]=aryutt[21]*aryutt[169];
   aryutt[163]=aryutt[163] + aryutt[164] + 3*aryutt[156] + aryutt[155];
   aryutt[166]=3*aryutt[168] + aryutt[167];
   aryutt[172]=aryutt[16]*aryutt[166];
   aryutt[173]=aryutt[15]*aryutt[166];
   aryutt[172]=aryutt[172] + 1./2.*aryutt[173];
   aryutt[172]=aryutt[3]*aryutt[172];
   aryutt[173]= - 1./4.*aryutt[10];
   aryutt[110]=aryutt[173] + aryutt[110] + aryutt[30];
   aryutt[174]=aryutt[21]*aryutt[110];
   aryutt[175]=3./2.*aryutt[30];
   aryutt[177]=aryutt[157] - aryutt[29] + aryutt[175];
   aryutt[177]=aryutt[12]*aryutt[177];
   aryutt[174]=1./2.*aryutt[177] + 3./2.*aryutt[149] + aryutt[174];
   aryutt[177]=aryutt[16]*aryutt[160];
   aryutt[178]=aryutt[15]*aryutt[160];
   aryutt[177]=aryutt[177] + 1./2.*aryutt[178];
   aryutt[177]=aryutt[3]*aryutt[177];
   aryutt[178]=MMZ*aryutt[3]*aryutt[160];
   aryutt[174]=3*aryutt[178] + 1./2.*aryutt[174] + 3*aryutt[177];
   aryutt[174]=aryutt[27]*aryutt[174];
   aryutt[166]=MMZ*aryutt[3]*aryutt[166];
   aryutt[177]=aryutt[29] - aryutt[30];
   aryutt[178]=MMt*aryutt[27]*aryutt[3]*aryutt[177];
   aryutt[152]=3./2.*aryutt[178] + 1./2.*aryutt[152] + aryutt[174] + 
   aryutt[166] + 1./4.*aryutt[163] + aryutt[172];
   aryutt[152]=aryutt[210]*aryutt[152];
   aryutt[163]=3*aryutt[10];
   aryutt[166]=aryutt[118] + aryutt[163];
   aryutt[172]=aryutt[31]*aryutt[166];
   aryutt[166]=aryutt[1]*aryutt[166];
   aryutt[166]=3*aryutt[172] + aryutt[166];
   aryutt[166]=aryutt[15]*aryutt[166];
   aryutt[154]=aryutt[162] + aryutt[19] - aryutt[7] + aryutt[154];
   aryutt[162]=aryutt[31]*aryutt[154];
   aryutt[154]=aryutt[1]*aryutt[154];
   aryutt[172]=1./2. + aryutt[10];
   aryutt[174]=aryutt[31]*aryutt[172];
   aryutt[172]=aryutt[1]*aryutt[172];
   aryutt[178]=3*aryutt[174] + aryutt[172];
   aryutt[178]=aryutt[16]*aryutt[178];
   aryutt[154]=1./4.*aryutt[166] + aryutt[178] + 3*aryutt[162] + 
   aryutt[154];
   aryutt[154]=aryutt[3]*aryutt[154];
   aryutt[143]=aryutt[12]*aryutt[143];
   aryutt[162]= - aryutt[10] + aryutt[13] - aryutt[25] - aryutt[24];
   aryutt[166]=aryutt[31]*aryutt[162];
   aryutt[162]=aryutt[1]*aryutt[162];
   aryutt[143]=aryutt[143] + 3*aryutt[166] + aryutt[162];
   aryutt[143]=aryutt[99]*aryutt[143];
   aryutt[178]= - 11 + 21./2.*aryutt[10];
   aryutt[178]=aryutt[31]*aryutt[178];
   aryutt[104]= - 11./3. + aryutt[104];
   aryutt[179]=aryutt[1]*aryutt[104];
   aryutt[178]=aryutt[178] + aryutt[179];
   aryutt[178]=aryutt[3]*aryutt[178];
   aryutt[143]=1./8.*aryutt[143] + aryutt[178];
   aryutt[143]=MMZ*aryutt[143];
   aryutt[178]= - 9*aryutt[28];
   aryutt[105]=aryutt[105] + 23./2.*aryutt[29] + 35./2. + aryutt[178];
   aryutt[105]=3./2.*aryutt[12] + 1./2.*aryutt[105] + aryutt[231];
   aryutt[105]=aryutt[3]*aryutt[105];
   aryutt[180]=aryutt[199] - 1./2. + aryutt[40];
   aryutt[180]=aryutt[17]*aryutt[153]*aryutt[180];
   aryutt[105]=1./2.*aryutt[105] + 9*aryutt[180];
   aryutt[105]=aryutt[27]*aryutt[105];
   aryutt[180]=MMt*aryutt[27]*aryutt[153]*aryutt[28];
   aryutt[105]=aryutt[105] + 9*aryutt[180];
   aryutt[105]=MMt*aryutt[105];
   aryutt[181]=1 - aryutt[10];
   aryutt[182]=aryutt[31]*aryutt[181];
   aryutt[181]=aryutt[1]*aryutt[181];
   aryutt[183]=3*aryutt[182] + aryutt[181];
   aryutt[183]=aryutt[17]*aryutt[99]*aryutt[183];
   aryutt[103]=1./2.*aryutt[152] + aryutt[105] + aryutt[117] + 
   aryutt[103] + 1./2.*aryutt[143] + 1./4.*aryutt[183] + 1./4.*
   aryutt[109] + aryutt[154];
   aryutt[103]=aryutt[210]*aryutt[103];
   aryutt[105]=aryutt[126] - aryutt[20];
   aryutt[109]=aryutt[31]*aryutt[105];
   aryutt[105]=aryutt[1]*aryutt[105];
   aryutt[117]=aryutt[109] + 1./3.*aryutt[105];
   aryutt[143]= - 1./3. + aryutt[10];
   aryutt[152]=aryutt[31]*aryutt[143];
   aryutt[143]=aryutt[1]*aryutt[143];
   aryutt[143]=5./3.*aryutt[152] + aryutt[143];
   aryutt[152]=aryutt[16]*aryutt[143];
   aryutt[154]=17./4.*aryutt[11];
   aryutt[183]=13./3. + aryutt[154];
   aryutt[113]=1./3.*aryutt[183] + aryutt[113];
   aryutt[113]=aryutt[31]*aryutt[113];
   aryutt[183]=67./12.*aryutt[10] + 1 + 17./12.*aryutt[11];
   aryutt[183]=aryutt[1]*aryutt[183];
   aryutt[113]=aryutt[113] + 1./3.*aryutt[183];
   aryutt[183]=aryutt[15]*aryutt[113];
   aryutt[184]= - 17./4.*aryutt[11];
   aryutt[185]=41./3. + aryutt[184];
   aryutt[185]=1./3.*aryutt[185] - 17./4.*aryutt[10];
   aryutt[185]=aryutt[31]*aryutt[185];
   aryutt[186]= - 67./12.*aryutt[10] + 5 - 17./12.*aryutt[11];
   aryutt[186]=aryutt[1]*aryutt[186];
   aryutt[185]=aryutt[185] + 1./3.*aryutt[186];
   aryutt[185]=aryutt[17]*aryutt[185];
   aryutt[117]=1./2.*aryutt[185] + 1./2.*aryutt[183] + 2*aryutt[117] + 
   aryutt[152];
   aryutt[152]= - aryutt[26] + aryutt[23];
   aryutt[125]=2./9.*aryutt[13] + 1./4.*aryutt[152] + aryutt[125];
   aryutt[125]=aryutt[1]*aryutt[125];
   aryutt[183]=7./27. + aryutt[11];
   aryutt[183]=1./2.*aryutt[183] + 19./9.*aryutt[10];
   aryutt[183]=aryutt[31]*aryutt[183];
   aryutt[185]= - 1./3. + aryutt[11];
   aryutt[185]=1./6.*aryutt[185] + aryutt[10];
   aryutt[185]=aryutt[1]*aryutt[185];
   aryutt[183]=aryutt[183] + aryutt[185];
   aryutt[183]=aryutt[21]*aryutt[183];
   aryutt[113]=aryutt[12]*aryutt[113];
   aryutt[152]=2./3.*aryutt[13] + 3./4.*aryutt[152] - 2./3.*aryutt[24];
   aryutt[152]=aryutt[31]*aryutt[152];
   aryutt[113]=1./6.*aryutt[113] + 1./2.*aryutt[183] + aryutt[152] + 
   aryutt[125];
   aryutt[113]=MMZ*aryutt[113];
   aryutt[125]=1./12. + aryutt[29];
   aryutt[152]=1./4.*aryutt[30];
   aryutt[125]=5./18.*aryutt[10] + 7./9.*aryutt[125] + aryutt[152];
   aryutt[125]=aryutt[21]*aryutt[125];
   aryutt[183]=19./18.*aryutt[13] - 19./18.*aryutt[24] + 3*aryutt[23]
    - aryutt[45] + 11./18. - 3*aryutt[43];
   aryutt[185]=17./3.*aryutt[30];
   aryutt[186]=aryutt[185] - 89./9. + 25./2.*aryutt[29];
   aryutt[186]=1./3.*aryutt[186] + 3./2.*aryutt[10];
   aryutt[186]=aryutt[12]*aryutt[186];
   aryutt[125]=1./8.*aryutt[186] + 1./4.*aryutt[183] + aryutt[125];
   aryutt[125]=MMZ*aryutt[125];
   aryutt[183]=11*aryutt[33];
   aryutt[186]=aryutt[183] - 19*aryutt[5];
   aryutt[188]= - 5*aryutt[18];
   aryutt[186]= - 41./2.*aryutt[20] + 1./2.*aryutt[186] + aryutt[188];
   aryutt[191]=4*aryutt[29];
   aryutt[193]=aryutt[10] - 5./3. + aryutt[191];
   aryutt[195]=aryutt[16]*aryutt[193];
   aryutt[186]=1./4.*aryutt[186] + aryutt[195];
   aryutt[195]= - 25./2.*aryutt[29];
   aryutt[196]= - 17./3.*aryutt[30];
   aryutt[197]=aryutt[196] + 11./9. + aryutt[195];
   aryutt[197]=1./3.*aryutt[197] - 3./2.*aryutt[10];
   aryutt[197]=19./12.*aryutt[12] + 1./8.*aryutt[197] + 14./9.*
   aryutt[21];
   aryutt[197]=aryutt[17]*aryutt[197];
   aryutt[198]=17./12.*aryutt[30] - 17./9. + 25./8.*aryutt[29];
   aryutt[198]=1./3.*aryutt[198] + 3./8.*aryutt[10];
   aryutt[198]=aryutt[15]*aryutt[198];
   aryutt[125]=aryutt[125] + aryutt[197] + 1./9.*aryutt[186] + 1./2.*
   aryutt[198];
   aryutt[125]=aryutt[27]*aryutt[125];
   aryutt[111]=aryutt[17]*aryutt[111];
   aryutt[186]=aryutt[12] + 1 + aryutt[45];
   aryutt[186]=MMZ*aryutt[186];
   aryutt[197]=1./2.*aryutt[33];
   aryutt[111]=1./2.*aryutt[186] + aryutt[111] + aryutt[270] + 
   aryutt[197] - aryutt[20];
   aryutt[111]=MMZ*aryutt[111];
   aryutt[186]=aryutt[241] - 1./4.*aryutt[17];
   aryutt[186]=aryutt[17]*aryutt[186];
   aryutt[111]=aryutt[186] + 1./6.*aryutt[111];
   aryutt[111]=aryutt[6]*aryutt[27]*aryutt[111];
   aryutt[113]=29./12.*aryutt[111] + aryutt[125] + 1./3.*aryutt[117] + 
   aryutt[113];
   aryutt[113]=aryutt[6]*aryutt[113];
   aryutt[117]= - 11*aryutt[45];
   aryutt[125]=323./12.*aryutt[29] + 49./3.*aryutt[13] - 49./3.*
   aryutt[24] + 25./3.*aryutt[25] - 31./2. + aryutt[117];
   aryutt[186]=1./3.*aryutt[30];
   aryutt[125]=1./2.*aryutt[125] + aryutt[186];
   aryutt[198]= - 1./3.*aryutt[10];
   aryutt[125]=1./2.*aryutt[125] + aryutt[198];
   aryutt[199]=aryutt[147] + 89./18. + aryutt[191];
   aryutt[199]=aryutt[21]*aryutt[199];
   aryutt[200]=2809./9. + 187*aryutt[29];
   aryutt[200]=1./4.*aryutt[200] - aryutt[30];
   aryutt[200]=1./4.*aryutt[200] + aryutt[10];
   aryutt[200]=aryutt[12]*aryutt[200];
   aryutt[125]=1./4.*aryutt[200] + 1./4.*aryutt[125] + aryutt[199];
   aryutt[199]=5./3. - 29./8.*aryutt[29];
   aryutt[200]= - 3./4.*aryutt[30];
   aryutt[199]=7./24.*aryutt[10] + 1./3.*aryutt[199] + aryutt[200];
   aryutt[199]=aryutt[15]*aryutt[199];
   aryutt[202]=aryutt[33] - aryutt[5];
   aryutt[202]=1./2.*aryutt[202] - aryutt[20];
   aryutt[203]= - 4*aryutt[29];
   aryutt[204]= - aryutt[10] + 5./3. + aryutt[203];
   aryutt[205]=aryutt[16]*aryutt[204];
   aryutt[205]=aryutt[199] + aryutt[202] + 2./3.*aryutt[205];
   aryutt[205]=aryutt[3]*aryutt[205];
   aryutt[206]=5 - 13./4.*aryutt[29];
   aryutt[207]= - 7*aryutt[10];
   aryutt[206]=3*aryutt[206] + aryutt[207];
   aryutt[206]=aryutt[15]*aryutt[99]*aryutt[206];
   aryutt[208]= - 7 + 39*aryutt[29];
   aryutt[209]=7*aryutt[10];
   aryutt[208]= - 39./4.*aryutt[12] + 1./4.*aryutt[208] + aryutt[209];
   aryutt[208]=aryutt[99]*aryutt[208];
   aryutt[211]= - 5./2.*aryutt[12];
   aryutt[212]= - 53./3. + aryutt[211];
   aryutt[212]=aryutt[3]*aryutt[212];
   aryutt[208]=1./2.*aryutt[208] + aryutt[212];
   aryutt[208]=aryutt[17]*aryutt[208];
   aryutt[207]=39./4.*aryutt[116] + aryutt[207];
   aryutt[207]=aryutt[12]*aryutt[207];
   aryutt[207]=aryutt[207] + aryutt[209] + 39./4.*aryutt[29] - 7*
   aryutt[13] + 7*aryutt[24] - 51./2. + 7*aryutt[25];
   aryutt[207]=aryutt[99]*aryutt[207];
   aryutt[195]=43./3. + aryutt[195];
   aryutt[195]= - 2./3.*aryutt[10] + 1./3.*aryutt[195] - 2*aryutt[30];
   aryutt[195]=aryutt[3]*aryutt[195];
   aryutt[195]=1./16.*aryutt[207] + aryutt[195];
   aryutt[195]=MMZ*aryutt[195];
   aryutt[207]=7./2.*aryutt[20] - 1./4.*aryutt[18] - 3./4.*aryutt[33]
    + aryutt[5];
   aryutt[207]=aryutt[99]*aryutt[207];
   aryutt[125]=aryutt[195] + 1./2.*aryutt[208] + aryutt[205] + 1./8.*
   aryutt[206] + 1./3.*aryutt[125] + 7./8.*aryutt[207];
   aryutt[125]=aryutt[27]*aryutt[125];
   aryutt[195]=1./3. - aryutt[10];
   aryutt[205]=aryutt[31]*aryutt[195];
   aryutt[195]=aryutt[1]*aryutt[195];
   aryutt[195]=5./3.*aryutt[205] + aryutt[195];
   aryutt[205]=aryutt[21]*aryutt[195];
   aryutt[206]=7./4.*aryutt[11];
   aryutt[207]=41./3. + aryutt[206];
   aryutt[207]=1./3.*aryutt[207] + aryutt[173];
   aryutt[207]=aryutt[31]*aryutt[207];
   aryutt[208]= - 19./12.*aryutt[10] + 5 + 7./12.*aryutt[11];
   aryutt[208]=aryutt[1]*aryutt[208];
   aryutt[207]=aryutt[207] + 1./3.*aryutt[208];
   aryutt[207]=aryutt[12]*aryutt[207];
   aryutt[208]=1./24.*aryutt[11];
   aryutt[209]=35./24.*aryutt[10] + aryutt[208] + aryutt[13] - 
   aryutt[24] + 7./6. + aryutt[25];
   aryutt[209]=aryutt[31]*aryutt[209];
   aryutt[212]=83./24.*aryutt[10] + aryutt[208] + aryutt[13] - 
   aryutt[24] + 1./2. + aryutt[25];
   aryutt[212]=aryutt[1]*aryutt[212];
   aryutt[207]=1./2.*aryutt[207] + 1./2.*aryutt[205] + aryutt[209] + 1./
   3.*aryutt[212];
   aryutt[209]=aryutt[5] - aryutt[18];
   aryutt[212]=2*aryutt[20];
   aryutt[214]=aryutt[209] + aryutt[212];
   aryutt[215]=aryutt[31]*aryutt[214];
   aryutt[214]=aryutt[1]*aryutt[214];
   aryutt[217]=aryutt[215] + 1./3.*aryutt[214];
   aryutt[217]=aryutt[99]*aryutt[217];
   aryutt[195]=aryutt[16]*aryutt[195];
   aryutt[218]= - 11./12.*aryutt[10] + 5./9. - 3./4.*aryutt[11];
   aryutt[218]=aryutt[31]*aryutt[218];
   aryutt[219]= - 1./4.*aryutt[11];
   aryutt[220]=1./3. + aryutt[219];
   aryutt[221]=aryutt[220] - 3./4.*aryutt[10];
   aryutt[221]=aryutt[1]*aryutt[221];
   aryutt[218]=aryutt[218] + aryutt[221];
   aryutt[218]=aryutt[15]*aryutt[218];
   aryutt[195]=2*aryutt[195] + aryutt[218];
   aryutt[195]=aryutt[3]*aryutt[195];
   aryutt[221]= - aryutt[31]*aryutt[10];
   aryutt[222]= - aryutt[1]*aryutt[10];
   aryutt[223]=aryutt[221] + 1./3.*aryutt[222];
   aryutt[224]=aryutt[12]*aryutt[223];
   aryutt[225]=aryutt[10] - aryutt[13] + aryutt[25] + aryutt[24];
   aryutt[226]=aryutt[31]*aryutt[225];
   aryutt[225]=aryutt[1]*aryutt[225];
   aryutt[224]=aryutt[224] + aryutt[226] + 1./3.*aryutt[225];
   aryutt[224]=aryutt[99]*aryutt[224];
   aryutt[227]= - 2*aryutt[11];
   aryutt[228]=17./3. + aryutt[227];
   aryutt[229]= - 5./2.*aryutt[10];
   aryutt[228]=1./3.*aryutt[228] + aryutt[229];
   aryutt[228]=aryutt[1]*aryutt[228];
   aryutt[227]= - 29./6.*aryutt[10] + 43./9. + aryutt[227];
   aryutt[227]=aryutt[31]*aryutt[227];
   aryutt[227]=aryutt[227] + aryutt[228];
   aryutt[227]=aryutt[3]*aryutt[227];
   aryutt[224]=aryutt[224] + aryutt[227];
   aryutt[224]=MMZ*aryutt[224];
   aryutt[227]= - 3*aryutt[28];
   aryutt[228]=1./2.*aryutt[30];
   aryutt[230]= - aryutt[12] + aryutt[228] - 77./3.*aryutt[29] - 37./6.
    + aryutt[227];
   aryutt[230]=aryutt[3]*aryutt[230];
   aryutt[231]= - 1 + aryutt[264];
   aryutt[231]=3*aryutt[17]*aryutt[153]*aryutt[231];
   aryutt[230]=1./2.*aryutt[230] + aryutt[231];
   aryutt[230]=aryutt[27]*aryutt[230];
   aryutt[232]=6*aryutt[180];
   aryutt[230]=aryutt[230] + aryutt[232];
   aryutt[230]=MMt*aryutt[230];
   aryutt[223]=aryutt[15]*aryutt[99]*aryutt[223];
   aryutt[233]= - 1 + aryutt[10];
   aryutt[234]=aryutt[31]*aryutt[233];
   aryutt[233]=aryutt[1]*aryutt[233];
   aryutt[235]=aryutt[234] + 1./3.*aryutt[233];
   aryutt[235]=aryutt[17]*aryutt[99]*aryutt[235];
   aryutt[103]=aryutt[103] + aryutt[230] + aryutt[113] + aryutt[125] + 
   aryutt[224] + 4*aryutt[235] + aryutt[195] + 2*aryutt[223] + 1./3.*
   aryutt[207] + 2*aryutt[217];
   aryutt[103]=aryutt[210]*aryutt[103];
   aryutt[113]= - 13*aryutt[39];
   aryutt[125]=11./3.*aryutt[10] - 20./9. + 3*aryutt[11];
   aryutt[125]=aryutt[31]*aryutt[125];
   aryutt[163]=aryutt[163] - 4./3. + aryutt[11];
   aryutt[163]=aryutt[1]*aryutt[163];
   aryutt[195]= - 55./4.*aryutt[12];
   aryutt[140]=9./2.*aryutt[3]*aryutt[140];
   aryutt[207]= - 15./8.*aryutt[21];
   aryutt[217]=aryutt[140] + aryutt[195] + aryutt[207] + aryutt[163] + 
   aryutt[125] + aryutt[201] + aryutt[213] + 1./4.*aryutt[41] + 
   aryutt[142] + aryutt[113] - 299./12. + aryutt[187];
   aryutt[217]=aryutt[17]*aryutt[3]*aryutt[217];
   aryutt[223]= - 17*aryutt[29];
   aryutt[224]=121./3. + aryutt[223];
   aryutt[230]= - 5./9.*aryutt[10];
   aryutt[235]=3*aryutt[30];
   aryutt[224]=aryutt[230] + 1./9.*aryutt[224] + aryutt[235];
   aryutt[224]=aryutt[21]*aryutt[224];
   aryutt[183]= - 25*aryutt[20] - 11*aryutt[18] + aryutt[242] + 3*
   aryutt[34] + aryutt[183];
   aryutt[236]=aryutt[175] - 139./6.*aryutt[29] - 59./3. + aryutt[178];
   aryutt[236]=aryutt[15]*aryutt[236];
   aryutt[237]= - 1 - aryutt[28];
   aryutt[237]=3*aryutt[237] - aryutt[30];
   aryutt[237]=aryutt[16]*aryutt[237];
   aryutt[237]=1./2.*aryutt[236] + aryutt[183] + 3./2.*aryutt[237];
   aryutt[237]=aryutt[3]*aryutt[237];
   aryutt[238]= - 20./3. + 17./2.*aryutt[29];
   aryutt[133]=aryutt[133] + 1./3.*aryutt[238] + aryutt[235];
   aryutt[133]=aryutt[27]*aryutt[17]*aryutt[3]*aryutt[133];
   aryutt[238]=13759./2. + 511*aryutt[29];
   aryutt[239]= - 5./2.*aryutt[30];
   aryutt[238]=1./27.*aryutt[238] + aryutt[239];
   aryutt[238]=aryutt[12]*aryutt[238];
   aryutt[240]=10./3. - 17./4.*aryutt[29];
   aryutt[240]= - 5./36.*aryutt[10] + 1./9.*aryutt[240] + aryutt[145];
   aryutt[240]=aryutt[8]*aryutt[240];
   aryutt[241]= - 3./8.*aryutt[43];
   aryutt[242]= - 49./27. + aryutt[241];
   aryutt[243]=3./16.*aryutt[23];
   aryutt[245]= - 19./144.*aryutt[30];
   aryutt[217]=aryutt[133] + aryutt[217] + 1./2.*aryutt[237] + 1./24.*
   aryutt[238] + aryutt[240] + 1./16.*aryutt[224] + aryutt[245] - 73./
   216.*aryutt[29] + 400./81.*aryutt[25] + aryutt[243] + 1./2.*
   aryutt[242] + 400./81.*aryutt[45];
   aryutt[217]=aryutt[27]*aryutt[217];
   aryutt[224]= - 11./9.*aryutt[10] + 47./27. - aryutt[11];
   aryutt[224]=aryutt[31]*aryutt[224];
   aryutt[237]=7./3. - aryutt[11];
   aryutt[237]=1./3.*aryutt[237] - aryutt[10];
   aryutt[237]=aryutt[1]*aryutt[237];
   aryutt[224]=aryutt[224] + aryutt[237];
   aryutt[237]=aryutt[21]*aryutt[224];
   aryutt[238]= - aryutt[31]*aryutt[22];
   aryutt[242]= - aryutt[1]*aryutt[22];
   aryutt[238]=aryutt[238] + 1./3.*aryutt[242];
   aryutt[237]=aryutt[238] + aryutt[237];
   aryutt[242]= - 1./2.*aryutt[11];
   aryutt[246]=2./3. + aryutt[242];
   aryutt[246]=1./3.*aryutt[246] + aryutt[157];
   aryutt[246]=aryutt[1]*aryutt[246];
   aryutt[242]= - 11./18.*aryutt[10] + 10./27. + aryutt[242];
   aryutt[242]=aryutt[31]*aryutt[242];
   aryutt[242]=aryutt[242] + aryutt[246];
   aryutt[246]=aryutt[8]*aryutt[242];
   aryutt[237]=1./8.*aryutt[237] + aryutt[246];
   aryutt[246]=aryutt[16]*aryutt[28];
   aryutt[248]=aryutt[15]*aryutt[28];
   aryutt[249]=aryutt[246] + aryutt[248];
   aryutt[249]=aryutt[3]*aryutt[249];
   aryutt[250]= - 11*aryutt[12] - 3./2.*aryutt[21] + 133./3. + 
   aryutt[239];
   aryutt[249]=1./2.*aryutt[250] + 9*aryutt[249];
   aryutt[249]=MMt*aryutt[27]*aryutt[3]*aryutt[249];
   aryutt[217]=aryutt[249] + aryutt[237] + aryutt[217];
   aryutt[217]=MMt*aryutt[217];
   aryutt[250]=107./4. - 259*aryutt[29];
   aryutt[252]= - 19./4.*aryutt[30];
   aryutt[250]=aryutt[229] + 1./9.*aryutt[250] + aryutt[252];
   aryutt[250]=887./324.*aryutt[12] + 1./9.*aryutt[250] + aryutt[279];
   aryutt[247]=aryutt[247] - 41./3.*aryutt[15];
   aryutt[247]=aryutt[3]*aryutt[247];
   aryutt[254]=2173./54.*aryutt[99] - 9*aryutt[3];
   aryutt[254]=aryutt[17]*aryutt[254];
   aryutt[256]= - aryutt[15]*aryutt[99];
   aryutt[247]=5./2.*aryutt[254] + 1./2.*aryutt[247] + 1./2.*
   aryutt[250] + 1685./27.*aryutt[256];
   aryutt[247]=aryutt[17]*aryutt[247];
   aryutt[250]=5./3. - 17./8.*aryutt[29];
   aryutt[250]=1./9.*aryutt[250];
   aryutt[254]= - 5./72.*aryutt[10];
   aryutt[257]=aryutt[254] + aryutt[250] + aryutt[152];
   aryutt[257]=aryutt[16]*aryutt[257];
   aryutt[250]=aryutt[254] + aryutt[250] + aryutt[148];
   aryutt[250]=aryutt[14]*aryutt[250];
   aryutt[258]=1663 + 4775./6.*aryutt[29];
   aryutt[259]= - 17*aryutt[30];
   aryutt[258]=1./3.*aryutt[258] + aryutt[259];
   aryutt[258]=aryutt[15]*aryutt[258];
   aryutt[253]= - aryutt[33] + aryutt[253];
   aryutt[253]=11./2.*aryutt[253];
   aryutt[260]=aryutt[253] + 7291./81.*aryutt[20];
   aryutt[247]=1./2.*aryutt[247] + 1./144.*aryutt[258] + 1./2.*
   aryutt[257] + 1./16.*aryutt[260] + aryutt[250];
   aryutt[247]=aryutt[27]*aryutt[247];
   aryutt[150]=3./2. + aryutt[150];
   aryutt[150]=9./2.*aryutt[170] + 35./12.*aryutt[12] + 1./2.*
   aryutt[150] + aryutt[142];
   aryutt[150]=aryutt[17]*aryutt[3]*aryutt[150];
   aryutt[257]=103 + 35./2.*aryutt[29];
   aryutt[257]=1./9.*aryutt[257] - aryutt[30];
   aryutt[258]= - 7*aryutt[29];
   aryutt[260]= - 13 + aryutt[258];
   aryutt[260]=7./144.*aryutt[260] - aryutt[30];
   aryutt[260]=aryutt[12]*aryutt[260];
   aryutt[257]=1./24.*aryutt[257] + aryutt[260];
   aryutt[260]= - 11./3. + aryutt[178];
   aryutt[260]=aryutt[175] + 1./2.*aryutt[260] + 7./3.*aryutt[29];
   aryutt[260]=aryutt[15]*aryutt[260];
   aryutt[261]= - aryutt[33] + aryutt[18];
   aryutt[261]=1./2.*aryutt[261] + aryutt[20];
   aryutt[260]=7./3.*aryutt[261] + 1./2.*aryutt[260];
   aryutt[260]=aryutt[3]*aryutt[260];
   aryutt[150]=aryutt[150] + 1./3.*aryutt[257] + aryutt[260];
   aryutt[150]=aryutt[27]*aryutt[150];
   aryutt[257]=7./3.*aryutt[12] + 11./3. + aryutt[30];
   aryutt[260]=aryutt[3]*aryutt[248];
   aryutt[257]=1./2.*aryutt[257] + 9*aryutt[260];
   aryutt[257]=MMt*aryutt[27]*aryutt[3]*aryutt[257];
   aryutt[150]=aryutt[150] + aryutt[257];
   aryutt[150]=MMt*aryutt[150];
   aryutt[257]=7*aryutt[29];
   aryutt[260]=323./27. + aryutt[257];
   aryutt[185]= - 193./18.*aryutt[12] + 1./2.*aryutt[260] + aryutt[185]
   ;
   aryutt[170]=67./2.*aryutt[170] + 1./4.*aryutt[185] + 725./9.*
   aryutt[107];
   aryutt[185]= - 4775./54.*aryutt[99] - 17*aryutt[3];
   aryutt[185]=aryutt[17]*aryutt[185];
   aryutt[170]=1./3.*aryutt[170] + 1./2.*aryutt[185];
   aryutt[170]=aryutt[17]*aryutt[170];
   aryutt[185]= - 221 - 301./2.*aryutt[29];
   aryutt[185]=1./6.*aryutt[185] + aryutt[259];
   aryutt[185]=aryutt[15]*aryutt[185];
   aryutt[260]=aryutt[33] + 77./9.*aryutt[18];
   aryutt[260]=1./2.*aryutt[260] - aryutt[20];
   aryutt[185]=7./2.*aryutt[260] + 1./3.*aryutt[185];
   aryutt[170]=1./12.*aryutt[185] + aryutt[170];
   aryutt[170]=aryutt[27]*aryutt[170];
   aryutt[185]= - 29./2.*aryutt[15] + 19*aryutt[17];
   aryutt[185]=aryutt[6]*aryutt[27]*aryutt[17]*aryutt[185];
   aryutt[170]=aryutt[170] + 17./27.*aryutt[185];
   aryutt[185]=1 + aryutt[39];
   aryutt[185]=aryutt[17]*aryutt[3]*aryutt[185];
   aryutt[260]= - aryutt[12]*aryutt[30];
   aryutt[185]=17./36.*aryutt[260] + aryutt[185];
   aryutt[185]=aryutt[27]*aryutt[185];
   aryutt[261]=aryutt[3]*aryutt[30];
   aryutt[262]=MMt*aryutt[27]*aryutt[261];
   aryutt[185]=aryutt[185] + aryutt[262];
   aryutt[185]=MMt*aryutt[185];
   aryutt[262]= - aryutt[27]*aryutt[17]*aryutt[12];
   aryutt[185]=17./36.*aryutt[262] + aryutt[185];
   aryutt[185]=aryutt[176]*aryutt[185];
   aryutt[150]=1./2.*aryutt[185] + 1./2.*aryutt[170] + aryutt[150];
   aryutt[150]=aryutt[176]*aryutt[150];
   aryutt[170]=1./3.*aryutt[220] + aryutt[173];
   aryutt[170]=aryutt[1]*aryutt[170];
   aryutt[173]= - 11./36.*aryutt[10] + 5./27. + aryutt[219];
   aryutt[173]=aryutt[31]*aryutt[173];
   aryutt[170]=aryutt[173] + aryutt[170];
   aryutt[173]=aryutt[14]*aryutt[170];
   aryutt[185]=aryutt[16]*aryutt[224];
   aryutt[170]=aryutt[17]*aryutt[170];
   aryutt[170]=aryutt[170] + aryutt[173] + 1./8.*aryutt[185];
   aryutt[173]=1./4.*aryutt[11];
   aryutt[185]=11./36.*aryutt[10] - 5./27. + aryutt[173];
   aryutt[185]=aryutt[31]*aryutt[185];
   aryutt[173]= - 1./3. + aryutt[173];
   aryutt[173]=1./3.*aryutt[173] + 1./4.*aryutt[10];
   aryutt[173]=aryutt[1]*aryutt[173];
   aryutt[173]=aryutt[185] + aryutt[173];
   aryutt[173]=aryutt[8]*aryutt[173];
   aryutt[185]= - 5./3. + 17./8.*aryutt[29];
   aryutt[185]=5./72.*aryutt[10] + 1./9.*aryutt[185] + aryutt[152];
   aryutt[185]=aryutt[27]*aryutt[8]*aryutt[185];
   aryutt[173]=aryutt[173] + aryutt[185];
   aryutt[173]=1./2.*MMH*aryutt[173];
   aryutt[185]=aryutt[16] + 3049./81.*aryutt[15];
   aryutt[185]=1./2.*aryutt[185] - 2165./81.*aryutt[17];
   aryutt[185]=aryutt[6]*aryutt[27]*aryutt[17]*aryutt[185];
   aryutt[150]=1./2.*aryutt[150] + aryutt[217] + 1./2.*aryutt[185] + 
   aryutt[173] + aryutt[170] + aryutt[247];
   aryutt[150]=aryutt[176]*aryutt[150];
   aryutt[113]=aryutt[140] + aryutt[195] + aryutt[207] + aryutt[163] + 
   aryutt[125] + aryutt[201] + aryutt[213] - 47./4.*aryutt[41] + 
   aryutt[142] + aryutt[113] - 443./12. + aryutt[187];
   aryutt[113]=aryutt[17]*aryutt[3]*aryutt[113];
   aryutt[125]= - 91./3. + aryutt[178];
   aryutt[125]=aryutt[200] + 1./4.*aryutt[125] - 16./3.*aryutt[29];
   aryutt[125]=aryutt[16]*aryutt[125];
   aryutt[125]=1./4.*aryutt[236] + 1./2.*aryutt[183] + aryutt[125];
   aryutt[125]=aryutt[3]*aryutt[125];
   aryutt[140]=aryutt[230] + aryutt[235] - 71./27. - 9*aryutt[29];
   aryutt[140]=aryutt[21]*aryutt[140];
   aryutt[142]=7./3. + aryutt[241];
   aryutt[163]=aryutt[239] + 391./6. - 19*aryutt[29];
   aryutt[163]=aryutt[12]*aryutt[163];
   aryutt[113]=aryutt[133] + aryutt[113] + aryutt[125] + 1./24.*
   aryutt[163] + aryutt[240] + 1./16.*aryutt[140] + aryutt[245] + 125./
   72.*aryutt[29] + 16./9.*aryutt[25] + aryutt[243] + 1./2.*aryutt[142]
    + 16./9.*aryutt[45];
   aryutt[113]=aryutt[27]*aryutt[113];
   aryutt[113]=aryutt[249] + aryutt[237] + aryutt[113];
   aryutt[113]=MMt*aryutt[113];
   aryutt[110]=aryutt[16]*aryutt[110];
   aryutt[125]= - aryutt[29] + aryutt[235];
   aryutt[125]=aryutt[194] + aryutt[21] + 1./2.*aryutt[125] - 
   aryutt[10];
   aryutt[125]=aryutt[17]*aryutt[125];
   aryutt[133]=aryutt[14]*aryutt[160];
   aryutt[140]= - aryutt[29] + aryutt[30];
   aryutt[142]=aryutt[15]*aryutt[140];
   aryutt[110]=1./2.*aryutt[125] + 1./4.*aryutt[142] + aryutt[133] + 
   aryutt[110];
   aryutt[110]=aryutt[27]*aryutt[110];
   aryutt[125]=aryutt[14]*aryutt[169];
   aryutt[133]=aryutt[16]*aryutt[169];
   aryutt[142]=aryutt[8]*aryutt[161];
   aryutt[163]=aryutt[27]*aryutt[8]*aryutt[149];
   aryutt[163]=aryutt[142] + aryutt[163];
   aryutt[163]=MMH*aryutt[163];
   aryutt[110]=1./2.*aryutt[163] + aryutt[110] + aryutt[171] + 
   aryutt[125] + 1./2.*aryutt[133];
   aryutt[125]=aryutt[235] - aryutt[10];
   aryutt[125]=aryutt[21]*aryutt[125];
   aryutt[125]=3*aryutt[140] + aryutt[125];
   aryutt[133]=aryutt[16]*aryutt[177];
   aryutt[163]=aryutt[15]*aryutt[177];
   aryutt[133]=aryutt[133] + 1./2.*aryutt[163];
   aryutt[133]=aryutt[3]*aryutt[133];
   aryutt[160]=aryutt[8]*aryutt[160];
   aryutt[163]=aryutt[12]*aryutt[29];
   aryutt[125]=3./2.*aryutt[133] + 1./8.*aryutt[163] + 1./8.*
   aryutt[125] + aryutt[160];
   aryutt[108]=aryutt[108] - 1 + 33*aryutt[39];
   aryutt[108]=aryutt[156] + aryutt[36] + 1./4.*aryutt[108] - 
   aryutt[38];
   aryutt[108]=3*aryutt[108] + aryutt[155];
   aryutt[108]=aryutt[17]*aryutt[3]*aryutt[108];
   aryutt[133]=aryutt[27]*aryutt[17]*aryutt[3]*aryutt[149];
   aryutt[108]=3*aryutt[133] + 1./2.*aryutt[125] + aryutt[108];
   aryutt[108]=aryutt[27]*aryutt[108];
   aryutt[125]=aryutt[8]*aryutt[169];
   aryutt[133]=1./4.*aryutt[164] + aryutt[125];
   aryutt[149]= - aryutt[3]*aryutt[30];
   aryutt[155]=MMt*aryutt[27]*aryutt[149];
   aryutt[108]=3./4.*aryutt[155] + 1./2.*aryutt[133] + aryutt[108];
   aryutt[108]=MMt*aryutt[108];
   aryutt[108]=1./4.*aryutt[110] + aryutt[108];
   aryutt[108]=aryutt[210]*aryutt[108];
   aryutt[110]= - 55 + 49*aryutt[29];
   aryutt[133]= - 7./2.*aryutt[99] - 3*aryutt[3];
   aryutt[133]=aryutt[17]*aryutt[133];
   aryutt[155]=1./6.*aryutt[10];
   aryutt[110]=3./2.*aryutt[133] + 25*aryutt[190] + 3*aryutt[107] + 35./
   24.*aryutt[12] + 11./6.*aryutt[21] + aryutt[155] + 1./24.*
   aryutt[110] - aryutt[30];
   aryutt[110]=aryutt[17]*aryutt[110];
   aryutt[133]= - aryutt[33] + aryutt[188];
   aryutt[133]=1./2.*aryutt[133] + aryutt[20];
   aryutt[155]=aryutt[155] + 5./6.*aryutt[29] - aryutt[30];
   aryutt[156]=aryutt[14]*aryutt[155];
   aryutt[133]=3./8.*aryutt[133] + aryutt[156];
   aryutt[156]=1./8.*aryutt[30];
   aryutt[160]=1./48.*aryutt[10] + aryutt[156] - 2./3. - 9./16.*
   aryutt[29];
   aryutt[160]=aryutt[16]*aryutt[160];
   aryutt[163]= - 29 - 49./2.*aryutt[29];
   aryutt[163]=aryutt[15]*aryutt[163];
   aryutt[110]=1./4.*aryutt[110] + 1./96.*aryutt[163] + 1./4.*
   aryutt[133] + aryutt[160];
   aryutt[110]=aryutt[27]*aryutt[110];
   aryutt[129]=aryutt[168] + aryutt[36] + aryutt[38] + 55./4.*
   aryutt[41] + aryutt[129] - 55./4.*aryutt[39] + 57./8. - aryutt[40];
   aryutt[129]=9./2.*aryutt[190] - 15./8.*aryutt[12] + aryutt[207] + 3*
   aryutt[129] + aryutt[167];
   aryutt[129]=aryutt[17]*aryutt[3]*aryutt[129];
   aryutt[133]=1./3.*aryutt[10];
   aryutt[160]=aryutt[133] + aryutt[235] + 25./3. + aryutt[257];
   aryutt[160]=aryutt[21]*aryutt[160];
   aryutt[163]=aryutt[23] - 5 - aryutt[43];
   aryutt[163]=3*aryutt[163] - 31./2.*aryutt[29];
   aryutt[160]=1./2.*aryutt[160] + 1./2.*aryutt[163] - aryutt[30];
   aryutt[155]=aryutt[8]*aryutt[155];
   aryutt[163]=1 + aryutt[29];
   aryutt[164]=29./6.*aryutt[163] + aryutt[30];
   aryutt[164]=aryutt[12]*aryutt[164];
   aryutt[155]=1./8.*aryutt[164] + 1./4.*aryutt[160] + aryutt[155];
   aryutt[160]=7 + aryutt[178];
   aryutt[164]=aryutt[200] + 1./4.*aryutt[160] + aryutt[191];
   aryutt[164]=aryutt[16]*aryutt[164];
   aryutt[167]=aryutt[293] - aryutt[19] + aryutt[34] + aryutt[197];
   aryutt[167]=1./2.*aryutt[167] - aryutt[20];
   aryutt[160]=aryutt[160] + 13*aryutt[29];
   aryutt[160]=aryutt[15]*aryutt[160];
   aryutt[160]=1./8.*aryutt[160] + 3*aryutt[167] + aryutt[164];
   aryutt[160]=aryutt[3]*aryutt[160];
   aryutt[164]=aryutt[157] - 5./2.*aryutt[29] + aryutt[235];
   aryutt[164]=aryutt[27]*aryutt[17]*aryutt[3]*aryutt[164];
   aryutt[129]=aryutt[164] + aryutt[129] + 1./2.*aryutt[155] + 
   aryutt[160];
   aryutt[129]=aryutt[27]*aryutt[129];
   aryutt[118]=aryutt[118] + aryutt[10];
   aryutt[155]=aryutt[31]*aryutt[118];
   aryutt[118]=aryutt[1]*aryutt[118];
   aryutt[118]=aryutt[155] + 1./3.*aryutt[118];
   aryutt[155]=aryutt[21]*aryutt[118];
   aryutt[155]=aryutt[238] + aryutt[155];
   aryutt[142]=1./4.*aryutt[155] + aryutt[142];
   aryutt[131]=aryutt[131] + aryutt[244] + 7./2. - aryutt[30];
   aryutt[155]=aryutt[246] + 1./2.*aryutt[248];
   aryutt[155]=aryutt[3]*aryutt[155];
   aryutt[131]=1./2.*aryutt[131] + 3*aryutt[155];
   aryutt[131]=MMt*aryutt[27]*aryutt[3]*aryutt[131];
   aryutt[129]=3*aryutt[131] + 1./2.*aryutt[142] + aryutt[129];
   aryutt[129]=MMt*aryutt[129];
   aryutt[131]=aryutt[14]*aryutt[161];
   aryutt[118]=aryutt[16]*aryutt[118];
   aryutt[142]=aryutt[17]*aryutt[161];
   aryutt[118]=aryutt[142] + aryutt[131] + 1./2.*aryutt[118];
   aryutt[131]=aryutt[147] - 5./6.*aryutt[29] + aryutt[30];
   aryutt[131]=aryutt[27]*aryutt[8]*aryutt[131];
   aryutt[125]=aryutt[125] + aryutt[131];
   aryutt[125]=MMH*aryutt[125];
   aryutt[131]=5*aryutt[271] + 7*aryutt[17];
   aryutt[131]=aryutt[6]*aryutt[27]*aryutt[17]*aryutt[131];
   aryutt[108]=aryutt[108] + aryutt[129] + 1./12.*aryutt[131] + 1./8.*
   aryutt[125] + 1./4.*aryutt[118] + aryutt[110];
   aryutt[108]=aryutt[210]*aryutt[108];
   aryutt[110]= - 23./2.*aryutt[21] + aryutt[229] + aryutt[252] + 761./
   12. - 43*aryutt[29];
   aryutt[110]=1./3.*aryutt[110] - 43./4.*aryutt[12];
   aryutt[118]= - 37*aryutt[16] - 41*aryutt[15];
   aryutt[118]=aryutt[3]*aryutt[118];
   aryutt[125]=11./2.*aryutt[99] - 5*aryutt[3];
   aryutt[125]=aryutt[17]*aryutt[125];
   aryutt[110]=9./2.*aryutt[125] + 1./6.*aryutt[118] + 1./6.*
   aryutt[110] + 15*aryutt[256];
   aryutt[110]=aryutt[17]*aryutt[110];
   aryutt[118]=53./9. + 37./8.*aryutt[29];
   aryutt[118]=aryutt[254] + 1./3.*aryutt[118] + aryutt[152];
   aryutt[118]=aryutt[16]*aryutt[118];
   aryutt[125]=aryutt[253] + 355./9.*aryutt[20];
   aryutt[129]=aryutt[196] + 71 + 101./2.*aryutt[29];
   aryutt[129]=aryutt[15]*aryutt[129];
   aryutt[110]=1./2.*aryutt[110] + 1./48.*aryutt[129] + 1./2.*
   aryutt[118] + 1./16.*aryutt[125] + aryutt[250];
   aryutt[110]=aryutt[27]*aryutt[110];
   aryutt[118]=41*aryutt[16] + 97*aryutt[15];
   aryutt[118]=1./2.*aryutt[118] - 77*aryutt[17];
   aryutt[118]=aryutt[6]*aryutt[27]*aryutt[17]*aryutt[118];
   aryutt[108]=aryutt[108] + aryutt[113] + 1./18.*aryutt[118] + 
   aryutt[173] + aryutt[170] + aryutt[110];
   aryutt[108]=aryutt[210]*aryutt[108];
   aryutt[110]=27*aryutt[35];
   aryutt[113]=3*aryutt[28];
   aryutt[118]=aryutt[113] - 365./12. + aryutt[110];
   aryutt[118]=13./36.*aryutt[30] - 7./36.*aryutt[29] + aryutt[122] + 1.
   /2.*aryutt[118] + aryutt[121];
   aryutt[125]= - 2*aryutt[15] + 12*aryutt[14] - 5./2.*aryutt[16];
   aryutt[125]=aryutt[3]*aryutt[125];
   aryutt[118]=9*aryutt[132] + aryutt[125] + aryutt[306] + 68./9.*
   aryutt[8] - 97./144.*aryutt[21] - aryutt[36] + 1./2.*aryutt[118] - 
   aryutt[38];
   aryutt[118]=aryutt[17]*aryutt[118];
   aryutt[125]=aryutt[19] + 9*aryutt[42] + 3*aryutt[32] - aryutt[34];
   aryutt[129]=3./8.*aryutt[125];
   aryutt[131]= - 3./4.*aryutt[20];
   aryutt[142]=aryutt[156] - 7./72.*aryutt[29] - 20./9. + 3./8.*
   aryutt[28];
   aryutt[142]=aryutt[14]*aryutt[142];
   aryutt[147]=1 - 7./2.*aryutt[29];
   aryutt[147]=1./9.*aryutt[147];
   aryutt[152]=aryutt[147] + 11./2.*aryutt[30];
   aryutt[152]=1./8.*aryutt[16]*aryutt[152];
   aryutt[155]= - aryutt[15]*aryutt[30];
   aryutt[160]=1./18.*aryutt[155];
   aryutt[161]=aryutt[27]*aryutt[255]*aryutt[3];
   aryutt[164]=26./3.*aryutt[161];
   aryutt[167]=aryutt[164] + aryutt[118] + aryutt[160] + aryutt[152] + 
   aryutt[142] + aryutt[131] + aryutt[129] + 400./81.*aryutt[18];
   aryutt[167]=aryutt[27]*aryutt[167];
   aryutt[168]= - 4*aryutt[20];
   aryutt[169]= - 3./2.*aryutt[30];
   aryutt[170]= - 1 + aryutt[169];
   aryutt[170]=aryutt[16]*aryutt[170];
   aryutt[171]= - aryutt[19] - 9./2.*aryutt[42] + 3./2.*aryutt[32] + 
   aryutt[34];
   aryutt[173]= - 1 - 1./2.*aryutt[28];
   aryutt[173]=aryutt[14]*aryutt[173];
   aryutt[170]=1./2.*aryutt[170] + aryutt[173] + aryutt[171] + 
   aryutt[168];
   aryutt[170]=aryutt[3]*aryutt[170];
   aryutt[178]=aryutt[169] + 7./6.*aryutt[29] + 26./3. + 9*aryutt[28];
   aryutt[178]=aryutt[3]*aryutt[178];
   aryutt[183]= - aryutt[17]*aryutt[153];
   aryutt[178]=aryutt[178] + 9*aryutt[183];
   aryutt[178]=aryutt[27]*aryutt[17]*aryutt[178];
   aryutt[185]= - 23./16. - aryutt[43];
   aryutt[115]= - 7./4.*aryutt[28] + 3./8.*aryutt[25] + aryutt[115] + 1.
   /4.*aryutt[185] - 3*aryutt[44];
   aryutt[115]=3*aryutt[115];
   aryutt[147]=aryutt[147] + aryutt[30];
   aryutt[147]=aryutt[21]*aryutt[147];
   aryutt[185]=15*aryutt[28];
   aryutt[187]=aryutt[30] - 7./9.*aryutt[29] - 403./9. + aryutt[185];
   aryutt[187]=aryutt[8]*aryutt[187];
   aryutt[188]= - 5./4.*aryutt[21];
   aryutt[190]=aryutt[188] - 3./2. - aryutt[28];
   aryutt[194]=aryutt[190] - 6*aryutt[8];
   aryutt[194]=aryutt[17]*aryutt[3]*aryutt[194];
   aryutt[147]=aryutt[178] + 3*aryutt[194] + 3*aryutt[170] + 1./9.*
   aryutt[260] + 1./4.*aryutt[187] + 1./8.*aryutt[147] + aryutt[115] - 
   5./72.*aryutt[30];
   aryutt[147]=aryutt[27]*aryutt[147];
   aryutt[170]= - 3 - aryutt[28];
   aryutt[170]=aryutt[8]*aryutt[170];
   aryutt[170]=aryutt[170] + aryutt[244] + 5./2. + aryutt[28];
   aryutt[170]=aryutt[3]*aryutt[170];
   aryutt[153]= - aryutt[27]*aryutt[17]*aryutt[153]*aryutt[28];
   aryutt[178]=aryutt[170] + 12*aryutt[153];
   aryutt[178]=MMt*aryutt[27]*aryutt[178];
   aryutt[147]=aryutt[147] + 3*aryutt[178];
   aryutt[147]=MMt*aryutt[147];
   aryutt[178]= - 3./2.*aryutt[28];
   aryutt[119]=aryutt[119] + 7./72.*aryutt[29] + 101./9. + aryutt[178];
   aryutt[119]=aryutt[8]*aryutt[119];
   aryutt[187]=1./4.*aryutt[28] + 1./2. + aryutt[44];
   aryutt[187]=9*aryutt[187];
   aryutt[119]=aryutt[187] + aryutt[119];
   aryutt[119]=1./2.*MMH*aryutt[27]*aryutt[119];
   aryutt[167]=aryutt[147] + aryutt[167] + aryutt[119];
   aryutt[167]=MMt*aryutt[167];
   aryutt[194]=aryutt[113] - 1091./36. + aryutt[110];
   aryutt[194]=35./36.*aryutt[30] + aryutt[122] + 1./2.*aryutt[194] + 
   aryutt[121];
   aryutt[195]=27./8.*aryutt[15] + 6*aryutt[14] + 73./8.*aryutt[16];
   aryutt[195]=aryutt[3]*aryutt[195];
   aryutt[132]=9./2.*aryutt[132];
   aryutt[114]=aryutt[132] + aryutt[195] + aryutt[114] + 4*aryutt[8] + 
   aryutt[130] + 1./4.*aryutt[194] - aryutt[38];
   aryutt[114]=aryutt[17]*aryutt[114];
   aryutt[130]=3./16.*aryutt[28];
   aryutt[194]=aryutt[156] - 1 + aryutt[130];
   aryutt[194]=aryutt[14]*aryutt[194];
   aryutt[195]=1 + 5*aryutt[30];
   aryutt[195]=aryutt[16]*aryutt[195];
   aryutt[125]=1./2.*aryutt[125] - aryutt[20];
   aryutt[125]=3./8.*aryutt[125];
   aryutt[161]=3*aryutt[161];
   aryutt[114]=aryutt[161] + aryutt[114] + 17./144.*aryutt[155] + 1./16.
   *aryutt[195] + aryutt[125] + aryutt[194];
   aryutt[114]=aryutt[27]*aryutt[114];
   aryutt[194]=1./2. + aryutt[30];
   aryutt[194]=aryutt[21]*aryutt[194];
   aryutt[195]= - 43 + aryutt[185];
   aryutt[195]=1./2.*aryutt[195] + aryutt[30];
   aryutt[195]=aryutt[8]*aryutt[195];
   aryutt[194]=7./72.*aryutt[260] + 1./2.*aryutt[195] + 1./4.*
   aryutt[194] + aryutt[115] - 19./72.*aryutt[30];
   aryutt[195]= - 2*aryutt[20];
   aryutt[171]=1./2.*aryutt[173] + 1./2.*aryutt[171] + aryutt[195];
   aryutt[173]= - 1 - aryutt[30];
   aryutt[173]=aryutt[16]*aryutt[173];
   aryutt[196]=aryutt[15]*aryutt[30];
   aryutt[173]=1./8.*aryutt[196] + aryutt[171] + 1./4.*aryutt[173];
   aryutt[173]=aryutt[3]*aryutt[173];
   aryutt[197]=aryutt[145] + 1 + 3./2.*aryutt[28];
   aryutt[197]=aryutt[3]*aryutt[197];
   aryutt[197]=aryutt[197] + 3./2.*aryutt[183];
   aryutt[197]=aryutt[27]*aryutt[17]*aryutt[197];
   aryutt[190]=1./2.*aryutt[190] - 3*aryutt[8];
   aryutt[190]=3*aryutt[17]*aryutt[3]*aryutt[190];
   aryutt[173]=3*aryutt[197] + aryutt[190] + 1./2.*aryutt[194] + 3*
   aryutt[173];
   aryutt[173]=aryutt[27]*aryutt[173];
   aryutt[153]=1./2.*aryutt[170] + 6*aryutt[153];
   aryutt[153]=3*MMt*aryutt[27]*aryutt[153];
   aryutt[170]=aryutt[173] + aryutt[153];
   aryutt[170]=MMt*aryutt[170];
   aryutt[173]=aryutt[148] + 11 + aryutt[178];
   aryutt[173]=aryutt[8]*aryutt[173];
   aryutt[173]=aryutt[187] + aryutt[173];
   aryutt[173]=MMH*aryutt[27]*aryutt[173];
   aryutt[114]=aryutt[170] + aryutt[114] + 1./4.*aryutt[173];
   aryutt[114]=MMt*aryutt[114];
   aryutt[170]=17*aryutt[30];
   aryutt[173]= - 7*aryutt[12] - 1 + aryutt[170];
   aryutt[194]= - aryutt[16] + 5./2.*aryutt[15];
   aryutt[194]=aryutt[3]*aryutt[194];
   aryutt[173]=1./36.*aryutt[173] + aryutt[194];
   aryutt[173]=aryutt[17]*aryutt[173];
   aryutt[173]=17./36.*aryutt[155] + aryutt[173];
   aryutt[173]=aryutt[27]*aryutt[173];
   aryutt[194]= - aryutt[30] + 7*aryutt[260];
   aryutt[197]=aryutt[3]*aryutt[196];
   aryutt[194]=1./18.*aryutt[194] + 3*aryutt[197];
   aryutt[194]=MMt*aryutt[27]*aryutt[194];
   aryutt[173]=aryutt[173] + 1./2.*aryutt[194];
   aryutt[173]=MMt*aryutt[173];
   aryutt[194]= - aryutt[15] + aryutt[17];
   aryutt[194]=aryutt[27]*aryutt[17]*aryutt[194];
   aryutt[173]=17./36.*aryutt[194] + aryutt[173];
   aryutt[173]=aryutt[176]*aryutt[173];
   aryutt[194]=9 + aryutt[273];
   aryutt[194]=aryutt[17]*aryutt[194];
   aryutt[165]=1./2.*aryutt[14] + aryutt[165] + aryutt[20];
   aryutt[165]=9*aryutt[165];
   aryutt[194]=aryutt[165] + aryutt[194];
   aryutt[194]=aryutt[27]*aryutt[194];
   aryutt[197]= - aryutt[8] - 1 - aryutt[44];
   aryutt[197]=9./2.*MMH*aryutt[27]*aryutt[197];
   aryutt[194]=aryutt[194] + aryutt[197];
   aryutt[194]=MMH*aryutt[194];
   aryutt[135]=143./18.*aryutt[17] - 17./18.*aryutt[15] + aryutt[135]
    + aryutt[16];
   aryutt[135]=aryutt[27]*aryutt[17]*aryutt[135];
   aryutt[135]=aryutt[135] + 1./2.*aryutt[194];
   aryutt[114]=1./4.*aryutt[173] + 1./8.*aryutt[135] + aryutt[114];
   aryutt[114]=aryutt[176]*aryutt[114];
   aryutt[135]= - 107*aryutt[14] + 19./2.*aryutt[16];
   aryutt[173]=1./2.*aryutt[135] - 809./9.*aryutt[15];
   aryutt[194]=19./2.*aryutt[17];
   aryutt[173]=1./3.*aryutt[173] + aryutt[194];
   aryutt[173]=aryutt[27]*aryutt[17]*aryutt[173];
   aryutt[200]=9 - 55./9.*aryutt[8];
   aryutt[200]=aryutt[17]*aryutt[200];
   aryutt[200]=aryutt[165] + aryutt[200];
   aryutt[200]=aryutt[27]*aryutt[200];
   aryutt[200]=aryutt[200] + aryutt[197];
   aryutt[200]=1./4.*MMH*aryutt[200];
   aryutt[173]=1./3.*aryutt[173] + aryutt[200];
   aryutt[114]=aryutt[114] + 1./2.*aryutt[173] + aryutt[167];
   aryutt[114]=aryutt[176]*aryutt[114];
   aryutt[118]=aryutt[164] + aryutt[118] + aryutt[160] + aryutt[152] + 
   aryutt[142] + aryutt[131] + aryutt[129] + 16./9.*aryutt[18];
   aryutt[118]=aryutt[27]*aryutt[118];
   aryutt[118]=aryutt[147] + aryutt[118] + aryutt[119];
   aryutt[118]=MMt*aryutt[118];
   aryutt[110]=aryutt[113] - 281./12. + aryutt[110];
   aryutt[110]=aryutt[148] + 19./6.*aryutt[29] + aryutt[122] + 1./2.*
   aryutt[110] + aryutt[121];
   aryutt[113]=2*aryutt[14] + 1./8.*aryutt[16];
   aryutt[113]=3*aryutt[113] - 43./8.*aryutt[15];
   aryutt[113]=aryutt[3]*aryutt[113];
   aryutt[119]= - 1./16.*aryutt[12];
   aryutt[110]=aryutt[132] + aryutt[113] + aryutt[119] + 16./3.*
   aryutt[8] + 1./48.*aryutt[21] + 1./4.*aryutt[110] - aryutt[36];
   aryutt[110]=aryutt[17]*aryutt[110];
   aryutt[113]=19./6.*aryutt[163] + aryutt[235];
   aryutt[113]=aryutt[16]*aryutt[113];
   aryutt[121]=19./24.*aryutt[29] - 1./3. + aryutt[130];
   aryutt[121]=aryutt[14]*aryutt[121];
   aryutt[122]= - aryutt[27]*aryutt[255]*aryutt[3];
   aryutt[110]=5*aryutt[122] + aryutt[110] + 1./16.*aryutt[196] + 1./8.
   *aryutt[113] + aryutt[125] + aryutt[121];
   aryutt[110]=aryutt[27]*aryutt[110];
   aryutt[113]=aryutt[21]*aryutt[163];
   aryutt[121]= - 97./3. + aryutt[185];
   aryutt[121]=1./2.*aryutt[121] + 19./3.*aryutt[29];
   aryutt[121]=aryutt[8]*aryutt[121];
   aryutt[113]=1./8.*aryutt[260] + 1./2.*aryutt[121] + 19./24.*
   aryutt[113] + aryutt[115] + aryutt[156];
   aryutt[115]= - 1./2. - aryutt[30];
   aryutt[115]=aryutt[16]*aryutt[115];
   aryutt[115]=1./8.*aryutt[155] + aryutt[171] + 1./2.*aryutt[115];
   aryutt[115]=aryutt[3]*aryutt[115];
   aryutt[121]= - 19./2.*aryutt[29];
   aryutt[125]=aryutt[121] - 5 + 9./2.*aryutt[28];
   aryutt[125]=aryutt[3]*aryutt[125];
   aryutt[125]=aryutt[125] + 9./2.*aryutt[183];
   aryutt[125]=aryutt[27]*aryutt[17]*aryutt[125];
   aryutt[113]=aryutt[125] + aryutt[190] + 1./2.*aryutt[113] + 3*
   aryutt[115];
   aryutt[113]=aryutt[27]*aryutt[113];
   aryutt[113]=aryutt[113] + aryutt[153];
   aryutt[113]=MMt*aryutt[113];
   aryutt[115]= - 19./12.*aryutt[29] + 29./3. + aryutt[178];
   aryutt[115]=aryutt[8]*aryutt[115];
   aryutt[115]=aryutt[187] + aryutt[115];
   aryutt[115]=MMH*aryutt[27]*aryutt[115];
   aryutt[110]=aryutt[113] + aryutt[110] + 1./4.*aryutt[115];
   aryutt[110]=MMt*aryutt[110];
   aryutt[113]=aryutt[169] + 3./2. + aryutt[29];
   aryutt[115]=aryutt[192] + 31./8.*aryutt[15];
   aryutt[115]=aryutt[3]*aryutt[115];
   aryutt[113]=aryutt[115] + aryutt[119] + aryutt[158] - aryutt[36] + 1.
   /8.*aryutt[113] + aryutt[38];
   aryutt[113]=aryutt[17]*aryutt[113];
   aryutt[115]=aryutt[14]*aryutt[177];
   aryutt[119]=aryutt[29] + aryutt[30];
   aryutt[119]=aryutt[16]*aryutt[119];
   aryutt[115]=1./2.*aryutt[196] + aryutt[115] + 1./2.*aryutt[119];
   aryutt[113]=1./8.*aryutt[115] + aryutt[113];
   aryutt[113]=aryutt[27]*aryutt[113];
   aryutt[115]=aryutt[21]*aryutt[123];
   aryutt[115]=aryutt[175] + aryutt[115];
   aryutt[119]= - aryutt[16]*aryutt[30];
   aryutt[123]=aryutt[119] + 1./2.*aryutt[155];
   aryutt[123]=aryutt[3]*aryutt[123];
   aryutt[125]=aryutt[8]*aryutt[177];
   aryutt[115]=3*aryutt[123] + 1./4.*aryutt[260] + 1./2.*aryutt[115] + 
   aryutt[125];
   aryutt[123]=aryutt[27]*aryutt[17]*aryutt[3]*aryutt[140];
   aryutt[115]=1./2.*aryutt[115] + 3*aryutt[123];
   aryutt[115]=MMt*aryutt[27]*aryutt[115];
   aryutt[123]=MMH*aryutt[27]*aryutt[8]*aryutt[140];
   aryutt[113]=1./2.*aryutt[115] + aryutt[113] + 1./16.*aryutt[123];
   aryutt[113]=MMt*aryutt[113];
   aryutt[115]=aryutt[189] + aryutt[120];
   aryutt[115]=aryutt[27]*aryutt[17]*aryutt[115];
   aryutt[113]=1./8.*aryutt[115] + aryutt[113];
   aryutt[113]=aryutt[210]*aryutt[113];
   aryutt[115]=9 - 37./3.*aryutt[8];
   aryutt[115]=aryutt[17]*aryutt[115];
   aryutt[115]=aryutt[165] + aryutt[115];
   aryutt[115]=aryutt[27]*aryutt[115];
   aryutt[115]=aryutt[115] + aryutt[197];
   aryutt[115]=MMH*aryutt[115];
   aryutt[123]= - 17./2.*aryutt[14] + 7*aryutt[16];
   aryutt[123]=71./12.*aryutt[17] + 1./3.*aryutt[123] + 1./4.*
   aryutt[15];
   aryutt[123]=aryutt[27]*aryutt[17]*aryutt[123];
   aryutt[115]=aryutt[123] + 1./4.*aryutt[115];
   aryutt[110]=aryutt[113] + 1./4.*aryutt[115] + aryutt[110];
   aryutt[110]=aryutt[210]*aryutt[110];
   aryutt[113]=aryutt[194] + 1./6.*aryutt[135] - 11*aryutt[15];
   aryutt[113]=aryutt[27]*aryutt[17]*aryutt[113];
   aryutt[113]=1./3.*aryutt[113] + aryutt[200];
   aryutt[110]=aryutt[110] + 1./2.*aryutt[113] + aryutt[118];
   aryutt[110]=aryutt[210]*aryutt[110];
   aryutt[113]=aryutt[14]*aryutt[30];
   aryutt[115]=1./2.*aryutt[16]*aryutt[30];
   aryutt[118]= - aryutt[20] - aryutt[34] + aryutt[19];
   aryutt[123]=aryutt[115] + 3*aryutt[118] + aryutt[113];
   aryutt[125]=aryutt[8] + aryutt[188] + 3 + aryutt[228];
   aryutt[125]=aryutt[17]*aryutt[125];
   aryutt[123]=1./2.*aryutt[123] + aryutt[125];
   aryutt[122]=3*aryutt[122];
   aryutt[123]=1./2.*aryutt[123] + aryutt[122];
   aryutt[123]=aryutt[27]*aryutt[123];
   aryutt[125]=aryutt[21]*aryutt[30];
   aryutt[129]= - aryutt[43] + aryutt[23];
   aryutt[129]=1./2.*aryutt[129] - aryutt[25];
   aryutt[130]=3*aryutt[129];
   aryutt[131]=aryutt[8]*aryutt[30];
   aryutt[132]=aryutt[131] + aryutt[130] + 1./4.*aryutt[125];
   aryutt[135]=3*aryutt[27]*aryutt[17]*aryutt[149];
   aryutt[132]=1./2.*aryutt[132] + aryutt[135];
   aryutt[132]=MMt*aryutt[27]*aryutt[132];
   aryutt[140]= - aryutt[8]*aryutt[30];
   aryutt[142]=1./8.*MMH*aryutt[27]*aryutt[140];
   aryutt[123]=aryutt[132] + aryutt[123] + aryutt[142];
   aryutt[123]=MMt*aryutt[123];
   aryutt[132]=aryutt[14] - 3./2.*aryutt[16];
   aryutt[132]=3*aryutt[132] + aryutt[17];
   aryutt[132]=aryutt[27]*aryutt[17]*aryutt[132];
   aryutt[147]=aryutt[38] - 1./4.*aryutt[8];
   aryutt[147]=MMH*aryutt[27]*aryutt[17]*aryutt[147];
   aryutt[132]=1./2.*aryutt[132] + aryutt[147];
   aryutt[123]=1./2.*aryutt[132] + aryutt[123];
   aryutt[123]=MMt*aryutt[123];
   aryutt[132]=3./2.*aryutt[118];
   aryutt[113]=aryutt[115] + aryutt[132] + aryutt[113];
   aryutt[115]= - aryutt[21] + 3 + aryutt[30];
   aryutt[115]=1./2.*aryutt[115] + aryutt[8];
   aryutt[115]=aryutt[17]*aryutt[115];
   aryutt[113]=1./2.*aryutt[113] + aryutt[115];
   aryutt[113]=1./2.*aryutt[113] + aryutt[122];
   aryutt[113]=aryutt[27]*aryutt[113];
   aryutt[115]=aryutt[130] + 1./2.*aryutt[125];
   aryutt[115]=1./2.*aryutt[115] + aryutt[131];
   aryutt[115]=1./2.*aryutt[115] + aryutt[135];
   aryutt[115]=MMt*aryutt[27]*aryutt[115];
   aryutt[113]=aryutt[115] + aryutt[113] + aryutt[142];
   aryutt[113]=MMt*aryutt[113];
   aryutt[115]=3*aryutt[136] + aryutt[17];
   aryutt[115]=aryutt[27]*aryutt[17]*aryutt[115];
   aryutt[115]=1./2.*aryutt[115] + aryutt[147];
   aryutt[113]=1./2.*aryutt[115] + aryutt[113];
   aryutt[113]=aryutt[176]*MMt*aryutt[113];
   aryutt[113]=aryutt[123] + aryutt[113];
   aryutt[113]=aryutt[176]*aryutt[113];
   aryutt[115]=1 + aryutt[244];
   aryutt[115]=aryutt[17]*aryutt[115];
   aryutt[115]=1./2.*aryutt[118] + aryutt[115];
   aryutt[115]=aryutt[27]*aryutt[115];
   aryutt[118]=MMt*aryutt[27]*aryutt[129];
   aryutt[115]=aryutt[115] + aryutt[118];
   aryutt[115]=MMt*aryutt[115];
   aryutt[118]= - aryutt[15] + aryutt[14] - 9./4.*aryutt[16];
   aryutt[118]=aryutt[27]*aryutt[17]*aryutt[118];
   aryutt[122]=MMH*aryutt[27]*aryutt[17]*aryutt[36];
   aryutt[115]=9./2.*aryutt[115] + aryutt[118] + aryutt[122];
   aryutt[115]=1./2.*MMt*aryutt[115];
   aryutt[113]=aryutt[115] + aryutt[113];
   aryutt[113]=aryutt[176]*aryutt[113];
   aryutt[118]= - aryutt[14] + 3./2.*aryutt[16];
   aryutt[118]=aryutt[120] + 1./2.*aryutt[118] - aryutt[15];
   aryutt[118]=aryutt[27]*aryutt[17]*aryutt[118];
   aryutt[122]=aryutt[288] + aryutt[159];
   aryutt[122]=MMH*aryutt[27]*aryutt[17]*aryutt[122];
   aryutt[118]=aryutt[118] + aryutt[122];
   aryutt[123]= - aryutt[14]*aryutt[30];
   aryutt[119]=1./2.*aryutt[119];
   aryutt[125]=aryutt[123] + aryutt[119];
   aryutt[129]= - aryutt[30] + aryutt[244];
   aryutt[129]=1./2.*aryutt[129] - aryutt[8];
   aryutt[129]=aryutt[17]*aryutt[129];
   aryutt[125]=1./2.*aryutt[125] + aryutt[129];
   aryutt[125]=1./2.*aryutt[125] + aryutt[161];
   aryutt[125]=aryutt[27]*aryutt[125];
   aryutt[129]= - aryutt[21]*aryutt[30];
   aryutt[135]=1./4.*aryutt[129] + aryutt[140];
   aryutt[136]=3*aryutt[27]*aryutt[17]*aryutt[261];
   aryutt[135]=1./2.*aryutt[135] + aryutt[136];
   aryutt[135]=MMt*aryutt[27]*aryutt[135];
   aryutt[131]=1./8.*MMH*aryutt[27]*aryutt[131];
   aryutt[125]=aryutt[135] + aryutt[125] + aryutt[131];
   aryutt[125]=MMt*aryutt[125];
   aryutt[118]=1./2.*aryutt[118] + aryutt[125];
   aryutt[118]=aryutt[210]*MMt*aryutt[118];
   aryutt[119]=aryutt[119] + aryutt[132] + aryutt[123];
   aryutt[123]=3 - aryutt[30];
   aryutt[123]= - aryutt[8] + 1./2.*aryutt[123] - aryutt[21];
   aryutt[123]=aryutt[17]*aryutt[123];
   aryutt[119]=1./2.*aryutt[119] + aryutt[123];
   aryutt[119]=1./2.*aryutt[119] + aryutt[161];
   aryutt[119]=aryutt[27]*aryutt[119];
   aryutt[123]=aryutt[130] + 1./2.*aryutt[129];
   aryutt[123]=1./2.*aryutt[123] + aryutt[140];
   aryutt[123]=1./2.*aryutt[123] + aryutt[136];
   aryutt[123]=MMt*aryutt[27]*aryutt[123];
   aryutt[119]=aryutt[123] + aryutt[119] + aryutt[131];
   aryutt[119]=MMt*aryutt[119];
   aryutt[120]=aryutt[120] - 1./2.*aryutt[14] - aryutt[15];
   aryutt[120]=aryutt[27]*aryutt[17]*aryutt[120];
   aryutt[120]=aryutt[120] + aryutt[122];
   aryutt[119]=1./2.*aryutt[120] + aryutt[119];
   aryutt[119]=MMt*aryutt[119];
   aryutt[118]=aryutt[119] + aryutt[118];
   aryutt[118]=aryutt[210]*aryutt[118];
   aryutt[115]=aryutt[115] + aryutt[118];
   aryutt[115]=aryutt[210]*aryutt[115];
   aryutt[113]=aryutt[113] + aryutt[115];
   aryutt[113]=aryutt[9]*aryutt[113];
   aryutt[115]=aryutt[27]*aryutt[17]*aryutt[15];
   aryutt[118]= - MMt*aryutt[27]*aryutt[18];
   aryutt[115]=aryutt[115] + aryutt[118];
   aryutt[110]=1./2.*aryutt[113] + aryutt[110] + 1024./81.*aryutt[115]
    + aryutt[114];
   aryutt[110]=aryutt[9]*aryutt[110];
   aryutt[113]= - aryutt[29] - aryutt[12];
   aryutt[107]=16*aryutt[216] + 1./3.*aryutt[113] + 10*aryutt[107];
   aryutt[107]=aryutt[17]*aryutt[107];
   aryutt[113]= - aryutt[18] - 8./3.*aryutt[20];
   aryutt[114]= - 1 - 1./9.*aryutt[29];
   aryutt[114]=aryutt[15]*aryutt[114];
   aryutt[107]=2./3.*aryutt[107] + 2./3.*aryutt[113] + aryutt[114];
   aryutt[107]=aryutt[27]*aryutt[107];
   aryutt[113]= - 5*aryutt[15] + 11*aryutt[17];
   aryutt[113]=aryutt[6]*aryutt[27]*aryutt[17]*aryutt[113];
   aryutt[114]= - 17 - aryutt[29];
   aryutt[114]=aryutt[12]*aryutt[114];
   aryutt[114]=aryutt[114] + aryutt[29] - 8*aryutt[25] + 1 - 8*
   aryutt[45];
   aryutt[114]=MMt*aryutt[27]*aryutt[114];
   aryutt[107]=2./9.*aryutt[114] + aryutt[107] + 2./9.*aryutt[113];
   aryutt[107]=aryutt[110] + aryutt[108] + 64./9.*aryutt[107] + 
   aryutt[150];
   aryutt[107]=aryutt[9]*aryutt[107];
   aryutt[108]=2059./3. - 1253*aryutt[29];
   aryutt[108]= - 35./6.*aryutt[10] + 1./24.*aryutt[108] + aryutt[259];
   aryutt[108]=aryutt[12]*aryutt[108];
   aryutt[110]=491./9. + 301./8.*aryutt[45];
   aryutt[113]= - 5*aryutt[5];
   aryutt[114]=41./4.*aryutt[33] + aryutt[113];
   aryutt[114]= - 61./6.*aryutt[20] + 1./3.*aryutt[114] - 7./4.*
   aryutt[18];
   aryutt[114]=aryutt[99]*aryutt[114];
   aryutt[112]=5./3.*aryutt[10] - 29./3. + aryutt[112];
   aryutt[112]=aryutt[15]*aryutt[99]*aryutt[112];
   aryutt[108]=25./12.*aryutt[112] + 25./12.*aryutt[114] + 1./12.*
   aryutt[108] + 25./216.*aryutt[10] - 185./864.*aryutt[29] - 35./72.*
   aryutt[13] + 35./72.*aryutt[24] + 1./9.*aryutt[110] + 5./8.*
   aryutt[25];
   aryutt[110]= - 17*aryutt[33] + aryutt[113];
   aryutt[112]=5*aryutt[10] + 11 + 17*aryutt[29];
   aryutt[112]=aryutt[15]*aryutt[112];
   aryutt[110]=1./2.*aryutt[112] + 17*aryutt[20] + 1./2.*aryutt[110] + 
   11*aryutt[18];
   aryutt[110]=aryutt[3]*aryutt[110];
   aryutt[112]= - 7 - 25*aryutt[29];
   aryutt[112]=25./4.*aryutt[12] + 1./4.*aryutt[112] + aryutt[151];
   aryutt[112]=aryutt[99]*aryutt[112];
   aryutt[113]=101./3. + 85./2.*aryutt[12];
   aryutt[113]=aryutt[3]*aryutt[113];
   aryutt[112]=25./6.*aryutt[112] + aryutt[113];
   aryutt[112]=aryutt[17]*aryutt[112];
   aryutt[113]=5./4.*aryutt[106] + aryutt[133];
   aryutt[113]=aryutt[12]*aryutt[113];
   aryutt[113]=5./3.*aryutt[113] + aryutt[230] - 25./12.*aryutt[29] + 5.
   /9.*aryutt[13] - 5./9.*aryutt[24] + 11./2. - 5./9.*aryutt[25];
   aryutt[113]=aryutt[99]*aryutt[113];
   aryutt[114]=35./4.*aryutt[10] - 121./3. + 119./4.*aryutt[29];
   aryutt[114]=aryutt[3]*aryutt[114];
   aryutt[113]=25./16.*aryutt[113] + 1./3.*aryutt[114];
   aryutt[113]=MMZ*aryutt[113];
   aryutt[108]=aryutt[113] + 1./2.*aryutt[112] + 1./2.*aryutt[108] + 
   aryutt[110];
   aryutt[108]=aryutt[27]*aryutt[108];
   aryutt[110]=5./27.*aryutt[10] - 7./9.*aryutt[13] + 7./9.*aryutt[24]
    - 23./324. + aryutt[25];
   aryutt[110]=aryutt[31]*aryutt[110];
   aryutt[112]=5./9.*aryutt[10] - 7./3.*aryutt[13] + 7./3.*aryutt[24]
    - 23./108. + 3*aryutt[25];
   aryutt[112]=aryutt[1]*aryutt[112];
   aryutt[110]=11./3.*aryutt[110] + aryutt[112];
   aryutt[112]= - 1./3. + aryutt[157];
   aryutt[113]=aryutt[31]*aryutt[112];
   aryutt[112]=aryutt[1]*aryutt[112];
   aryutt[112]=11./9.*aryutt[113] + aryutt[112];
   aryutt[112]=aryutt[12]*aryutt[112];
   aryutt[110]=1./2.*aryutt[110] + 7./3.*aryutt[112];
   aryutt[112]=11./9.*aryutt[137] + aryutt[134];
   aryutt[112]=aryutt[99]*aryutt[112];
   aryutt[113]=11./9.*aryutt[139] + aryutt[141];
   aryutt[114]=aryutt[15]*aryutt[99]*aryutt[113];
   aryutt[110]=25./6.*aryutt[114] + 1./2.*aryutt[110] + 25./3.*
   aryutt[112];
   aryutt[112]=11./3.*aryutt[174] + 3*aryutt[172];
   aryutt[112]=aryutt[15]*aryutt[112];
   aryutt[114]=aryutt[31]*aryutt[126];
   aryutt[115]=aryutt[1]*aryutt[126];
   aryutt[112]=aryutt[112] + 11./3.*aryutt[114] + 3*aryutt[115];
   aryutt[112]=aryutt[3]*aryutt[112];
   aryutt[113]=aryutt[12]*aryutt[113];
   aryutt[113]=aryutt[113] + 11./9.*aryutt[166] + aryutt[162];
   aryutt[113]=aryutt[99]*aryutt[113];
   aryutt[104]=aryutt[31]*aryutt[104];
   aryutt[104]=11./9.*aryutt[104] + aryutt[179];
   aryutt[104]=aryutt[3]*aryutt[104];
   aryutt[104]=25./24.*aryutt[113] + aryutt[104];
   aryutt[104]=MMZ*aryutt[104];
   aryutt[113]=11./9.*aryutt[182] + aryutt[181];
   aryutt[113]=aryutt[17]*aryutt[99]*aryutt[113];
   aryutt[104]=1./3.*aryutt[108] + aryutt[104] + 25./6.*aryutt[113] + 1.
   /2.*aryutt[110] + aryutt[112];
   aryutt[108]=aryutt[31]*aryutt[138];
   aryutt[110]=11./9.*aryutt[108] + aryutt[144];
   aryutt[112]=aryutt[12]*aryutt[110];
   aryutt[113]=aryutt[24] - aryutt[13];
   aryutt[114]=aryutt[31]*aryutt[113];
   aryutt[113]=aryutt[1]*aryutt[113];
   aryutt[112]=aryutt[112] + 11./9.*aryutt[114] + aryutt[113];
   aryutt[112]=MMZ*aryutt[112];
   aryutt[110]=aryutt[15]*aryutt[110];
   aryutt[115]= - 7./6. + aryutt[10];
   aryutt[118]=aryutt[31]*aryutt[115];
   aryutt[115]=aryutt[1]*aryutt[115];
   aryutt[115]=11./9.*aryutt[118] + aryutt[115];
   aryutt[115]=aryutt[17]*aryutt[115];
   aryutt[118]=aryutt[209] + aryutt[20];
   aryutt[119]=aryutt[31]*aryutt[118];
   aryutt[118]=aryutt[1]*aryutt[118];
   aryutt[110]=aryutt[112] + aryutt[115] + aryutt[110] + 11./9.*
   aryutt[119] + aryutt[118];
   aryutt[112]=aryutt[251] + 1./4.*aryutt[17];
   aryutt[112]=aryutt[17]*aryutt[112];
   aryutt[112]=aryutt[112] + 1./6.*aryutt[146];
   aryutt[112]=aryutt[6]*aryutt[27]*aryutt[112];
   aryutt[115]=aryutt[33] + 5*aryutt[5];
   aryutt[118]=aryutt[229] + 13./3. - 17./2.*aryutt[29];
   aryutt[118]=aryutt[15]*aryutt[118];
   aryutt[115]=1./6.*aryutt[118] + 1./4.*aryutt[20] + 1./12.*
   aryutt[115] + aryutt[18];
   aryutt[118]= - 13./8.*aryutt[12] + 5./16.*aryutt[10] - 2./3. + 17./
   16.*aryutt[29];
   aryutt[118]=aryutt[17]*aryutt[118];
   aryutt[115]=1./4.*aryutt[115] + 1./3.*aryutt[118];
   aryutt[118]= - 5*aryutt[10] + 83./3. + aryutt[223];
   aryutt[118]=aryutt[12]*aryutt[118];
   aryutt[118]=1./18.*aryutt[118] - 5./18.*aryutt[13] + 5./18.*
   aryutt[24] + 1./18. + aryutt[45];
   aryutt[118]=MMZ*aryutt[118];
   aryutt[115]=1./3.*aryutt[115] + 1./8.*aryutt[118];
   aryutt[115]=aryutt[27]*aryutt[115];
   aryutt[110]=17./24.*aryutt[112] + 1./8.*aryutt[110] + aryutt[115];
   aryutt[110]=aryutt[6]*aryutt[110];
   aryutt[112]=17./3.*aryutt[12] + aryutt[30] + 49./18.*aryutt[29] - 71.
   /18. + aryutt[227];
   aryutt[112]=aryutt[3]*aryutt[112];
   aryutt[112]=1./2.*aryutt[112] + aryutt[231];
   aryutt[112]=aryutt[27]*aryutt[112];
   aryutt[112]=1./2.*aryutt[112] + 3*aryutt[180];
   aryutt[112]=MMt*aryutt[112];
   aryutt[104]=aryutt[112] + 1./2.*aryutt[104] + 17./3.*aryutt[110];
   aryutt[104]=aryutt[176]*aryutt[104];
   aryutt[110]=11./9.*aryutt[215] + aryutt[214];
   aryutt[110]=aryutt[99]*aryutt[110];
   aryutt[112]=517./36.*aryutt[10] + 625./27. + aryutt[206];
   aryutt[112]=aryutt[31]*aryutt[112];
   aryutt[115]=173./3. + aryutt[206];
   aryutt[115]=1./3.*aryutt[115] + 47./4.*aryutt[10];
   aryutt[115]=aryutt[1]*aryutt[115];
   aryutt[112]=aryutt[112] + aryutt[115];
   aryutt[112]=aryutt[12]*aryutt[112];
   aryutt[115]=11./9.*aryutt[221] + aryutt[222];
   aryutt[118]=aryutt[15]*aryutt[99]*aryutt[115];
   aryutt[119]=11*aryutt[13] - 11*aryutt[24] + 131./6. + 11*aryutt[25];
   aryutt[119]= - 143./24.*aryutt[10] + 5./9.*aryutt[119] + 1./8.*
   aryutt[11];
   aryutt[119]=aryutt[31]*aryutt[119];
   aryutt[120]=aryutt[208] + 5*aryutt[13] - 5*aryutt[24] + 179./18. + 5
   *aryutt[25];
   aryutt[120]=1./3.*aryutt[120] - 13./8.*aryutt[10];
   aryutt[120]=aryutt[1]*aryutt[120];
   aryutt[110]=10*aryutt[118] + 10*aryutt[110] + 1./6.*aryutt[112] + 1./
   3.*aryutt[119] + aryutt[120];
   aryutt[112]=67./12.*aryutt[29] + 377./27.*aryutt[13] - 377./27.*
   aryutt[24] + 161./27.*aryutt[25] - 11951./162. + aryutt[117];
   aryutt[112]=1./2.*aryutt[112] + aryutt[186];
   aryutt[117]=54977./3. + 16313*aryutt[29];
   aryutt[117]=1./108.*aryutt[117] - aryutt[30];
   aryutt[117]=1./4.*aryutt[117] + 43./27.*aryutt[10];
   aryutt[117]=aryutt[12]*aryutt[117];
   aryutt[118]= - 523./4.*aryutt[33] + 25*aryutt[5];
   aryutt[118]=623./18.*aryutt[20] + 1./9.*aryutt[118] + 47./4.*
   aryutt[18];
   aryutt[118]=aryutt[99]*aryutt[118];
   aryutt[119]= - 25./3.*aryutt[10] + 337./3. - 275./4.*aryutt[29];
   aryutt[119]=aryutt[15]*aryutt[99]*aryutt[119];
   aryutt[112]=5./6.*aryutt[119] + 5./2.*aryutt[118] + aryutt[117] + 1./
   2.*aryutt[112] + aryutt[151];
   aryutt[116]=11./4.*aryutt[116] + aryutt[198];
   aryutt[116]=aryutt[12]*aryutt[116];
   aryutt[117]=25./9.*aryutt[10];
   aryutt[116]=25./3.*aryutt[116] + aryutt[117] + 275./12.*aryutt[29]
    - 25./9.*aryutt[13] + 25./9.*aryutt[24] - 133./2. + 25./9.*
   aryutt[25];
   aryutt[116]=aryutt[99]*aryutt[116];
   aryutt[118]=10./3. + aryutt[121];
   aryutt[118]=4./9.*aryutt[10] + 1./9.*aryutt[118] + aryutt[145];
   aryutt[118]=aryutt[3]*aryutt[118];
   aryutt[116]=5./48.*aryutt[116] + aryutt[118];
   aryutt[116]=MMZ*aryutt[116];
   aryutt[118]=aryutt[202] + aryutt[199];
   aryutt[118]=aryutt[3]*aryutt[118];
   aryutt[119]=47 + 275./3.*aryutt[29];
   aryutt[117]= - 275./12.*aryutt[12] + 1./4.*aryutt[119] + aryutt[117]
   ;
   aryutt[117]=aryutt[99]*aryutt[117];
   aryutt[119]= - 47./9. + aryutt[211];
   aryutt[119]=aryutt[3]*aryutt[119];
   aryutt[117]=5./6.*aryutt[117] + aryutt[119];
   aryutt[117]=aryutt[17]*aryutt[117];
   aryutt[112]=aryutt[116] + 1./2.*aryutt[117] + 1./12.*aryutt[112] + 
   aryutt[118];
   aryutt[112]=aryutt[27]*aryutt[112];
   aryutt[116]= - aryutt[24] + aryutt[13];
   aryutt[117]=aryutt[31]*aryutt[116];
   aryutt[116]=aryutt[1]*aryutt[116];
   aryutt[116]=11./9.*aryutt[117] + aryutt[116];
   aryutt[117]=1067./36.*aryutt[10] + 245./27. + aryutt[154];
   aryutt[117]=aryutt[31]*aryutt[117];
   aryutt[118]=73./3. + aryutt[154];
   aryutt[118]=1./3.*aryutt[118] + 97./4.*aryutt[10];
   aryutt[118]=aryutt[1]*aryutt[118];
   aryutt[117]=aryutt[117] + aryutt[118];
   aryutt[118]=aryutt[12]*aryutt[117];
   aryutt[116]=10*aryutt[116] + 1./2.*aryutt[118];
   aryutt[116]=MMZ*aryutt[116];
   aryutt[105]=11./9.*aryutt[109] + aryutt[105];
   aryutt[109]=aryutt[15]*aryutt[117];
   aryutt[117]= - 1067./36.*aryutt[10] + 745./27. + aryutt[184];
   aryutt[117]=aryutt[31]*aryutt[117];
   aryutt[118]=197./3. + aryutt[184];
   aryutt[118]=1./3.*aryutt[118] - 97./4.*aryutt[10];
   aryutt[118]=aryutt[1]*aryutt[118];
   aryutt[117]=aryutt[117] + aryutt[118];
   aryutt[117]=aryutt[17]*aryutt[117];
   aryutt[105]=aryutt[116] + 1./2.*aryutt[117] + 10*aryutt[105] + 1./2.
   *aryutt[109];
   aryutt[109]= - 601./3. - 391./2.*aryutt[29];
   aryutt[109]= - 179./18.*aryutt[10] + 5./9.*aryutt[109] + aryutt[259]
   ;
   aryutt[109]=1./2.*aryutt[109] + 185*aryutt[12];
   aryutt[109]=aryutt[17]*aryutt[109];
   aryutt[116]=611*aryutt[33] - 107*aryutt[5];
   aryutt[116]= - 443./6.*aryutt[20] + 1./18.*aryutt[116] - 37*
   aryutt[18];
   aryutt[117]= - 793./3. + 1955./8.*aryutt[29];
   aryutt[117]=179./72.*aryutt[10] + 1./9.*aryutt[117] + 17./4.*
   aryutt[30];
   aryutt[117]=aryutt[15]*aryutt[117];
   aryutt[109]=1./2.*aryutt[109] + 1./2.*aryutt[116] + aryutt[117];
   aryutt[116]= - 365./3. + 391./2.*aryutt[29];
   aryutt[116]=179./18.*aryutt[10] + 5./9.*aryutt[116] + aryutt[170];
   aryutt[116]=aryutt[12]*aryutt[116];
   aryutt[116]=1./18.*aryutt[116] + 107./162.*aryutt[13] - 107./162.*
   aryutt[24] + 611./162. - aryutt[45];
   aryutt[116]=MMZ*aryutt[116];
   aryutt[109]=1./9.*aryutt[109] + 1./2.*aryutt[116];
   aryutt[109]=aryutt[27]*aryutt[109];
   aryutt[105]=773./108.*aryutt[111] + 1./9.*aryutt[105] + 1./2.*
   aryutt[109];
   aryutt[105]=aryutt[6]*aryutt[105];
   aryutt[109]=aryutt[12]*aryutt[115];
   aryutt[109]=aryutt[109] + 11./9.*aryutt[226] + aryutt[225];
   aryutt[109]=aryutt[99]*aryutt[109];
   aryutt[111]=aryutt[3]*aryutt[242];
   aryutt[109]=5./3.*aryutt[109] + aryutt[111];
   aryutt[109]=MMZ*aryutt[109];
   aryutt[111]= - aryutt[12] + aryutt[228] - 119./9.*aryutt[29] + 113./
   18. + aryutt[227];
   aryutt[111]=aryutt[3]*aryutt[111];
   aryutt[111]=1./2.*aryutt[111] + aryutt[231];
   aryutt[111]=aryutt[27]*aryutt[111];
   aryutt[111]=aryutt[111] + aryutt[232];
   aryutt[111]=MMt*aryutt[111];
   aryutt[115]=aryutt[3]*aryutt[218];
   aryutt[116]=11./9.*aryutt[234] + aryutt[233];
   aryutt[116]=aryutt[17]*aryutt[99]*aryutt[116];
   aryutt[104]=aryutt[104] + aryutt[111] + aryutt[105] + aryutt[112] + 
   aryutt[109] + 20./3.*aryutt[116] + 1./3.*aryutt[110] + aryutt[115];
   aryutt[104]=aryutt[176]*aryutt[104];
   aryutt[105]=61./9. - 4*aryutt[45];
   aryutt[109]= - 2./3.*aryutt[29];
   aryutt[105]=aryutt[133] + aryutt[109] + 1./3.*aryutt[13] - 1./3.*
   aryutt[24] + 2*aryutt[105] - 23./3.*aryutt[25];
   aryutt[110]= - 5./3. + aryutt[29];
   aryutt[110]=2*aryutt[110] + aryutt[133];
   aryutt[110]=aryutt[15]*aryutt[99]*aryutt[110];
   aryutt[111]= - 1 - aryutt[29];
   aryutt[111]=aryutt[21]*aryutt[111];
   aryutt[112]= - 97./3. + aryutt[258];
   aryutt[112]=2*aryutt[112] - aryutt[10];
   aryutt[112]=aryutt[12]*aryutt[112];
   aryutt[115]=4*aryutt[33] - aryutt[5];
   aryutt[115]= - 10./3.*aryutt[20] + 1./3.*aryutt[115] - aryutt[18];
   aryutt[115]=aryutt[99]*aryutt[115];
   aryutt[116]=2*aryutt[12] + aryutt[198] - 1 - 2*aryutt[29];
   aryutt[116]=aryutt[99]*aryutt[116];
   aryutt[116]=2*aryutt[116] + aryutt[3];
   aryutt[116]=aryutt[17]*aryutt[116];
   aryutt[105]=4*aryutt[116] + 4*aryutt[110] + 4*aryutt[115] + 2./9.*
   aryutt[112] + 2./3.*aryutt[105] + aryutt[111];
   aryutt[106]=2*aryutt[106] + aryutt[133];
   aryutt[106]=aryutt[12]*aryutt[106];
   aryutt[106]=1./3.*aryutt[106] - 1./9.*aryutt[10] + aryutt[109] + 1./
   9.*aryutt[13] - 1./9.*aryutt[24] + 2 - 1./9.*aryutt[25];
   aryutt[106]=aryutt[99]*aryutt[106];
   aryutt[109]=aryutt[3]*aryutt[193];
   aryutt[106]=4*aryutt[106] + 1./3.*aryutt[109];
   aryutt[106]=MMZ*aryutt[106];
   aryutt[105]=2./3.*aryutt[105] + aryutt[106];
   aryutt[105]=aryutt[27]*aryutt[105];
   aryutt[106]= - 4./3. - aryutt[10];
   aryutt[109]=aryutt[31]*aryutt[106];
   aryutt[106]=aryutt[1]*aryutt[106];
   aryutt[106]=5./3.*aryutt[109] + aryutt[106];
   aryutt[106]=aryutt[12]*aryutt[106];
   aryutt[109]=aryutt[10] + aryutt[13] - aryutt[24] + 16./3. + 
   aryutt[25];
   aryutt[110]=aryutt[31]*aryutt[109];
   aryutt[109]=aryutt[1]*aryutt[109];
   aryutt[106]=aryutt[106] + 5./3.*aryutt[110] + aryutt[109];
   aryutt[109]=aryutt[126] + aryutt[195];
   aryutt[110]=aryutt[31]*aryutt[109];
   aryutt[109]=aryutt[1]*aryutt[109];
   aryutt[109]=5./3.*aryutt[110] + aryutt[109];
   aryutt[109]=aryutt[99]*aryutt[109];
   aryutt[110]=5./3.*aryutt[139] + aryutt[141];
   aryutt[111]=aryutt[15]*aryutt[99]*aryutt[110];
   aryutt[112]=5./3.*aryutt[182] + aryutt[181];
   aryutt[112]=aryutt[17]*aryutt[99]*aryutt[112];
   aryutt[106]=4*aryutt[112] + 2*aryutt[111] + 1./3.*aryutt[106] + 2*
   aryutt[109];
   aryutt[109]=aryutt[12]*aryutt[110];
   aryutt[109]=aryutt[109] + 5./3.*aryutt[166] + aryutt[162];
   aryutt[109]=aryutt[99]*aryutt[109];
   aryutt[110]=aryutt[3]*aryutt[143];
   aryutt[109]=4./3.*aryutt[109] + aryutt[110];
   aryutt[109]=MMZ*aryutt[109];
   aryutt[105]=aryutt[105] + 4./3.*aryutt[106] + aryutt[109];
   aryutt[106]=5./3.*aryutt[108] + aryutt[144];
   aryutt[108]=aryutt[15]*aryutt[106];
   aryutt[109]= - 11./3. + aryutt[10];
   aryutt[110]=aryutt[31]*aryutt[109];
   aryutt[109]=aryutt[1]*aryutt[109];
   aryutt[109]=5./3.*aryutt[110] + aryutt[109];
   aryutt[109]=aryutt[17]*aryutt[109];
   aryutt[110]=aryutt[209] + 4*aryutt[20];
   aryutt[111]=aryutt[31]*aryutt[110];
   aryutt[110]=aryutt[1]*aryutt[110];
   aryutt[108]=2*aryutt[109] + 2*aryutt[108] + 5./3.*aryutt[111] + 
   aryutt[110];
   aryutt[109]=5./3.*aryutt[114] + aryutt[113];
   aryutt[106]=aryutt[12]*aryutt[106];
   aryutt[106]=32./9.*aryutt[106] + 16./9.*aryutt[109] + aryutt[205];
   aryutt[106]=MMZ*aryutt[106];
   aryutt[109]= - 4*aryutt[33] + aryutt[5];
   aryutt[110]= - aryutt[10] + 7./3. + aryutt[203];
   aryutt[110]=aryutt[15]*aryutt[110];
   aryutt[109]=2./3.*aryutt[110] + aryutt[168] + 1./3.*aryutt[109] + 
   aryutt[18];
   aryutt[110]=aryutt[10] + 17./3. + aryutt[191];
   aryutt[110]=aryutt[124] + 4./9.*aryutt[110] - aryutt[21];
   aryutt[110]=aryutt[17]*aryutt[110];
   aryutt[109]=2./3.*aryutt[109] + aryutt[110];
   aryutt[110]= - aryutt[13] - 4 + aryutt[24];
   aryutt[111]=aryutt[21]*aryutt[204];
   aryutt[112]= - aryutt[10] + 1./3. + aryutt[203];
   aryutt[112]=aryutt[12]*aryutt[112];
   aryutt[110]=32./9.*aryutt[112] + 16./9.*aryutt[110] + aryutt[111];
   aryutt[110]=MMZ*aryutt[110];
   aryutt[109]=8*aryutt[109] + aryutt[110];
   aryutt[109]=aryutt[27]*aryutt[109];
   aryutt[110]=aryutt[128] + 2*aryutt[127] + aryutt[15] - aryutt[33] + 
   aryutt[212];
   aryutt[110]=MMZ*aryutt[110];
   aryutt[111]=aryutt[251] + aryutt[17];
   aryutt[111]=aryutt[17]*aryutt[111];
   aryutt[110]=4*aryutt[111] + 1./3.*aryutt[110];
   aryutt[110]=aryutt[6]*aryutt[27]*aryutt[110];
   aryutt[106]=32./9.*aryutt[110] + 1./3.*aryutt[109] + 16./9.*
   aryutt[108] + aryutt[106];
   aryutt[106]=aryutt[6]*aryutt[106];
   aryutt[108]=MMt*aryutt[27]*aryutt[3]*aryutt[163];
   aryutt[105]=32./3.*aryutt[108] + 4*aryutt[105] + aryutt[106];

      yuttret = aryutt[102] + aryutt[103] + aryutt[104] + 1./3.*
      aryutt[105] + aryutt[107];
      return yuttret;
}
