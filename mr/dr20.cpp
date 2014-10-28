#include <dr.hpp>
std::complex<long double>
dr::dr20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> ardr[208], drret;

    ardr[1]=double(nL + nH);
    ardr[2]=pow(CW,-1);
    ardr[3]=pow(MMH,-1);
    ardr[4]=pow(MMZ,-1);
    ardr[5]=pow(SW,-1);
    ardr[6]=Tsil::I2(0,0,MMZ,mu2);
    ardr[7]=Tsil::I2(0,0,MMW,mu2);
    ardr[8]=Tsil::A(MMH,mu2);
    ardr[9]=std::real(Tsil::B(0,0,MMZ,mu2));
    ardr[10]=std::real(Tsil::B(0,0,MMW,mu2));
    ardr[11]=Tsil::A(MMZ,mu2);
    ardr[12]=Tsil::A(MMW,mu2);
    ardr[13]=Tsil::Aeps(MMZ,mu2);
    ardr[14]=Tsil::Aeps(MMW,mu2);
    ardr[15]=double(nH);
    ardr[16]=Tsil::A(MMt,mu2);
    ardr[17]=Tsil::B(MMt,MMt,MMH,mu2);
    ardr[18]=Tsil::B(MMt,MMt,MMZ,mu2);
    ardr[19]=Tsil::B(0,MMt,MMW,mu2);
    ardr[20]=double(nL);
    ardr[21]=Tsil::I2(MMH,MMt,MMt,mu2);
    ardr[22]=Tsil::I2(MMZ,MMt,MMt,mu2);
    ardr[23]=pow(MMt,-1);
    ardr[24]=Tsil::I2(MMW,MMt,MMt,mu2);
    ardr[25]=Tsil::I2(0,MMH,MMt,mu2);
    ardr[26]=Tsil::I2(0,MMZ,MMt,mu2);
    ardr[27]=Tsil::I2(0,MMW,MMt,mu2);
    ardr[28]=Tsil::I2(0,0,MMt,mu2);
    ardr[29]=Tsil::B(MMH,MMH,MMH,mu2);
    ardr[30]=Tsil::B(MMH,MMt,MMt,mu2);
    ardr[31]=Tsil::B(MMZ,MMH,MMZ,mu2);
    ardr[32]=Tsil::B(MMZ,MMZ,MMH,mu2);
    ardr[33]=Tsil::B(MMZ,MMt,MMt,mu2);
    ardr[34]=Tsil::B(MMW,MMH,MMW,mu2);
    ardr[35]=Tsil::B(MMW,MMZ,MMW,mu2);
    ardr[36]=Tsil::B(MMW,MMW,MMH,mu2);
    ardr[37]=Tsil::B(MMW,MMW,MMZ,mu2);
    ardr[38]=std::real(Tsil::B(0,MMW,MMt,mu2));
    ardr[39]=Tsil::Aeps(MMH,mu2);
    ardr[40]=Tsil::Aeps(MMt,mu2);
    ardr[41]=double(boson);
    ardr[42]=Tsil::I2(MMH,MMH,MMH,mu2);
    ardr[43]=Tsil::I2(MMH,MMZ,MMZ,mu2);
    ardr[44]=Tsil::I2(MMH,MMW,MMW,mu2);
    ardr[45]=Tsil::I2(MMW,MMH,MMH,mu2);
    ardr[46]=Tsil::I2(MMW,MMZ,MMH,mu2);
    ardr[47]=Tsil::I2(MMW,MMZ,MMZ,mu2);
    ardr[48]=Tsil::I2(MMW,MMW,MMZ,mu2);
    ardr[49]=Tsil::I2(MMW,MMW,MMW,mu2);
    ardr[50]=Tsil::I2(0,MMZ,MMH,mu2);
    ardr[51]=Tsil::I2(0,MMW,MMH,mu2);
    ardr[52]=Tsil::I2(0,MMW,MMZ,mu2);
    ardr[53]=Tsil::I2(0,0,MMH,mu2);
    ardr[54]=1/(4*MMt - MMZ);
    ardr[55]=1/( - 4*MMW + MMH);
    ardr[56]=1/( - MMW + MMH);
    ardr[57]=1/( - MMW + 4*MMZ);
    ardr[58]=1/( - 4*MMZ + MMH);
   ardr[59]=55./2.*ardr[37];
   ardr[60]= - 55./2.*ardr[35];
   ardr[61]=2*ardr[31];
   ardr[62]=ardr[10] - ardr[9];
   ardr[63]=ardr[20]*ardr[62];
   ardr[64]=2*ardr[34];
   ardr[65]=2*ardr[63] + ardr[64] + ardr[61] + ardr[60] - ardr[32] + 
   ardr[59] - 2*ardr[36] - ardr[33] + 57./4. - ardr[38];
   ardr[62]=ardr[1]*ardr[62];
   ardr[65]=3*ardr[65] + 2*ardr[62];
   ardr[65]=ardr[16]*ardr[65];
   ardr[66]= - 2*ardr[40];
   ardr[67]=1./2.*ardr[22];
   ardr[68]= - 1./2.*ardr[13];
   ardr[69]= - ardr[14] + ardr[68] + ardr[66] + ardr[67] + ardr[27];
   ardr[70]= - 9*ardr[17];
   ardr[71]=23./2. + ardr[70];
   ardr[72]=8*ardr[18];
   ardr[73]= - 3./2.*ardr[19];
   ardr[74]=ardr[73] + 1./2.*ardr[71] + ardr[72];
   ardr[74]=ardr[12]*ardr[74];
   ardr[71]=ardr[71] + 13*ardr[18];
   ardr[71]=ardr[11]*ardr[71];
   ardr[65]=ardr[65] + ardr[74] + 3*ardr[69] + 1./4.*ardr[71];
   ardr[65]=ardr[3]*ardr[65];
   ardr[69]= - 9./8.*ardr[17];
   ardr[71]= - 3./8.*ardr[19];
   ardr[74]=ardr[71] + 2*ardr[18] + 2 + ardr[69];
   ardr[74]=ardr[12]*ardr[74];
   ardr[75]= - 1./4.*ardr[25];
   ardr[76]= - ardr[14] + 1./2.*ardr[39] - 5./4.*ardr[40] - 1./4.*
   ardr[21] + ardr[75] + ardr[24];
   ardr[77]=3./2.*ardr[76];
   ardr[78]= - 11./4. - 2*ardr[18];
   ardr[78]=ardr[8]*ardr[78];
   ardr[79]= - 15./8.*ardr[16];
   ardr[74]=ardr[79] + ardr[74] + ardr[77] + ardr[78];
   ardr[74]=ardr[56]*ardr[74];
   ardr[78]= - 3*ardr[17];
   ardr[80]=1./4.*ardr[19];
   ardr[81]=ardr[78] + ardr[80];
   ardr[81]=ardr[8]*ardr[81];
   ardr[82]= - 1./4.*ardr[19];
   ardr[83]=3*ardr[17];
   ardr[84]=ardr[83] + ardr[82];
   ardr[84]=ardr[12]*ardr[84];
   ardr[85]=3*ardr[21];
   ardr[75]=ardr[84] + ardr[81] + 11./4.*ardr[14] - 11./4.*ardr[39] + 1.
   /4.*ardr[27] + ardr[85] + ardr[75] - 3*ardr[24];
   ardr[75]=ardr[56]*ardr[75];
   ardr[75]=ardr[75] + ardr[82] - 1./4. + ardr[83];
   ardr[75]=ardr[56]*ardr[75];
   ardr[81]=ardr[11]*ardr[17];
   ardr[84]=ardr[12]*ardr[17];
   ardr[86]=ardr[81] + 2*ardr[84];
   ardr[86]=ardr[3]*ardr[86];
   ardr[87]= - ardr[33] + 7 - ardr[38];
   ardr[86]=3*ardr[86] + 1./2.*ardr[87] - ardr[19];
   ardr[86]=ardr[3]*ardr[86];
   ardr[86]=1./2.*ardr[75] + ardr[86];
   ardr[86]=MMt*ardr[86];
   ardr[87]=9*ardr[38];
   ardr[88]= - 165*ardr[37] + 9*ardr[33] - 325./4. + ardr[87];
   ardr[88]=165./2.*ardr[35] + 1./2.*ardr[88] + ardr[70];
   ardr[89]= - 3./4.*ardr[31];
   ardr[90]= - ardr[10] + ardr[9];
   ardr[91]=ardr[20]*ardr[90];
   ardr[90]=ardr[1]*ardr[90];
   ardr[92]= - 3./4.*ardr[34];
   ardr[93]=3./4.*ardr[19];
   ardr[65]=3*ardr[86] + ardr[65] + ardr[74] + 1./4.*ardr[90] + 3./4.*
   ardr[91] + ardr[93] + ardr[92] + ardr[89] + 1./8.*ardr[88] - 3*
   ardr[18];
   ardr[65]=MMt*ardr[65];
   ardr[74]=55./4.*ardr[35];
   ardr[86]=1./4.*ardr[38];
   ardr[88]=ardr[91] - ardr[34] - ardr[31] + ardr[74] - 55./4.*ardr[37]
    + 1./4.*ardr[33] - 9 + ardr[86];
   ardr[94]= - 5./2.*ardr[54] - ardr[23];
   ardr[94]=ardr[11]*ardr[94];
   ardr[95]= - ardr[12]*ardr[23];
   ardr[96]=5*ardr[54] + ardr[23];
   ardr[96]=ardr[16]*ardr[96];
   ardr[88]=3./4.*ardr[96] + 3*ardr[95] + 3./2.*ardr[94] + 3*ardr[88]
    + ardr[90];
   ardr[88]=ardr[16]*ardr[88];
   ardr[94]=5*ardr[18] + ardr[9];
   ardr[96]= - 3*ardr[19];
   ardr[94]=1./2.*ardr[94] + ardr[96];
   ardr[97]= - 5*ardr[18] - ardr[9];
   ardr[98]=6*ardr[19];
   ardr[99]=ardr[97] + ardr[98];
   ardr[99]=ardr[3]*ardr[16]*ardr[99];
   ardr[99]=1./4.*ardr[94] + ardr[99];
   ardr[99]=MMt*ardr[99];
   ardr[100]=ardr[16]*ardr[94];
   ardr[99]=1./2.*ardr[100] + ardr[99];
   ardr[99]=ardr[15]*ardr[99];
   ardr[94]=ardr[8]*ardr[94];
   ardr[100]= - 7./2.*ardr[9] + 133./96. - 11*ardr[18];
   ardr[100]=1./3.*ardr[100] + ardr[93];
   ardr[100]=ardr[11]*ardr[100];
   ardr[101]=1./2.*ardr[11];
   ardr[102]=ardr[101] + ardr[12];
   ardr[103]= - 3*ardr[16];
   ardr[102]=31*ardr[102] + ardr[103];
   ardr[102]=ardr[3]*ardr[16]*ardr[102];
   ardr[104]= - 1./3.*ardr[24];
   ardr[105]=3./8.*ardr[21];
   ardr[106]=5./12.*ardr[7];
   ardr[107]=ardr[27] + ardr[22] - ardr[24] - ardr[26];
   ardr[107]=EPAIR2*ardr[107];
   ardr[108]=1./6.*ardr[14];
   ardr[109]=19*ardr[18];
   ardr[110]=7*ardr[9] - 13./3. + ardr[109];
   ardr[110]=ardr[12]*ardr[110];
   ardr[97]=1./6.*ardr[97] + ardr[19];
   ardr[97]=MMH*ardr[97];
   ardr[111]= - 31*ardr[8] + 25*ardr[12];
   ardr[111]=ardr[56]*ardr[16]*ardr[111];
   ardr[112]= - 1./4.*ardr[27];
   ardr[65]=ardr[99] + ardr[65] + 1./2.*ardr[102] + 1./8.*ardr[111] + 1.
   /8.*ardr[97] + 1./2.*ardr[88] + 1./24.*ardr[110] + ardr[100] + 1./4.
   *ardr[94] + ardr[108] + 3./4.*ardr[107] - 151./96.*ardr[13] - 25./8.
   *ardr[40] + ardr[112] - 7./24.*ardr[6] + 173./96.*ardr[22] + 
   ardr[106] + 1./16.*ardr[26] + ardr[105] - 3./8.*ardr[25] + ardr[104]
   ;
   ardr[65]=ardr[15]*ardr[65];
   ardr[88]=143*ardr[35] - 111./2. - 209*ardr[37];
   ardr[94]=ardr[8]*ardr[58];
   ardr[97]= - ardr[11]*ardr[58];
   ardr[99]=5*ardr[31];
   ardr[100]= - 11*ardr[34];
   ardr[88]=3*ardr[97] + 3*ardr[94] + ardr[100] + 3./4.*ardr[88] + 
   ardr[99];
   ardr[88]=ardr[11]*ardr[88];
   ardr[102]= - 11*ardr[35];
   ardr[107]=ardr[102] + 15./2. + 77*ardr[37];
   ardr[107]= - 31*ardr[34] + 3*ardr[107] + 55*ardr[31];
   ardr[110]= - ardr[8]*ardr[55];
   ardr[111]=ardr[12]*ardr[55];
   ardr[107]=3*ardr[111] + 1./4.*ardr[107] + 3*ardr[110];
   ardr[107]=ardr[12]*ardr[107];
   ardr[113]=5*ardr[34];
   ardr[114]=ardr[37] - ardr[35];
   ardr[115]=ardr[113] + 33*ardr[114] - 5*ardr[31];
   ardr[115]=ardr[8]*ardr[115];
   ardr[116]=ardr[12] - ardr[8] - ardr[11];
   ardr[116]=ardr[12]*ardr[116];
   ardr[117]=ardr[11]*ardr[8];
   ardr[116]=ardr[117] + ardr[116];
   ardr[116]=ardr[56]*ardr[116];
   ardr[118]= - 1./2.*ardr[51];
   ardr[119]=3./8.*ardr[44] + 1./2.*ardr[50] + 99./8.*ardr[48] + 
   ardr[118] - 11./16.*ardr[47] - 187./16.*ardr[49] + 3*ardr[52];
   ardr[120]= - ardr[37] + ardr[35];
   ardr[121]=ardr[34] + 33./16.*ardr[120] - ardr[31];
   ardr[121]=MMH*ardr[121];
   ardr[122]= - 1./2.*ardr[7];
   ardr[123]=1./2.*ardr[6];
   ardr[124]= - 9./8.*ardr[43];
   ardr[88]=ardr[116] + ardr[121] + 1./2.*ardr[107] + 1./2.*ardr[88] + 
   3./8.*ardr[115] + 95./4.*ardr[14] - 167./4.*ardr[13] + ardr[124] + 
   ardr[123] + 3*ardr[119] + ardr[122];
   ardr[107]=2*ardr[12];
   ardr[115]= - ardr[11] + ardr[107];
   ardr[115]=ardr[12]*ardr[115];
   ardr[116]=pow(ardr[11],2);
   ardr[115]= - ardr[116] + ardr[115];
   ardr[115]=ardr[3]*ardr[115];
   ardr[88]=1./2.*ardr[88] + ardr[115];
   ardr[88]=ardr[41]*ardr[88];
   ardr[115]=11*ardr[37];
   ardr[119]=ardr[102] - 1./4. + ardr[115];
   ardr[119]=ardr[63] + ardr[34] + 3./4.*ardr[119] - ardr[31];
   ardr[121]=3./8.*ardr[54] - ardr[23];
   ardr[121]=ardr[11]*ardr[121];
   ardr[125]=ardr[12]*ardr[23];
   ardr[126]= - ardr[16]*ardr[54];
   ardr[119]=9./8.*ardr[126] + 3*ardr[125] + 3*ardr[121] + 3*ardr[119]
    + ardr[62];
   ardr[119]=ardr[16]*ardr[119];
   ardr[121]= - 33*ardr[37];
   ardr[125]=33*ardr[35] - 1 + ardr[121];
   ardr[127]= - 2*ardr[34];
   ardr[61]=2*ardr[91] + ardr[127] + 1./2.*ardr[125] + ardr[61];
   ardr[61]=3*ardr[61] + 2*ardr[90];
   ardr[61]=ardr[16]*ardr[61];
   ardr[125]=ardr[18] - ardr[19];
   ardr[128]=ardr[12]*ardr[125];
   ardr[129]=ardr[11]*ardr[125];
   ardr[130]=1./2.*ardr[129] + ardr[128];
   ardr[61]=3./2.*ardr[130] + ardr[61];
   ardr[61]=ardr[3]*ardr[61];
   ardr[115]= - ardr[18] + ardr[102] + 3./2. + ardr[115];
   ardr[115]=ardr[63] + ardr[93] + ardr[34] + 3./4.*ardr[115] - 
   ardr[31];
   ardr[130]= - ardr[18] + ardr[19];
   ardr[131]=ardr[8]*ardr[130];
   ardr[128]=ardr[131] + ardr[128];
   ardr[128]=ardr[56]*ardr[128];
   ardr[115]=3./2.*ardr[128] + 3*ardr[115] + ardr[62];
   ardr[128]= - MMt*ardr[3]*ardr[19];
   ardr[61]=3./2.*ardr[128] + 1./4.*ardr[115] + ardr[61];
   ardr[61]=MMt*ardr[61];
   ardr[115]=1./2.*ardr[9];
   ardr[128]=ardr[96] + ardr[115] - 3 - 7./2.*ardr[18];
   ardr[128]=ardr[12]*ardr[128];
   ardr[131]=ardr[18] + ardr[9];
   ardr[132]=ardr[131] - 2*ardr[19];
   ardr[132]=ardr[3]*ardr[16]*ardr[132];
   ardr[133]= - ardr[18] - ardr[9];
   ardr[133]=1./2.*ardr[133] + ardr[19];
   ardr[132]=1./4.*ardr[133] + ardr[132];
   ardr[132]=MMt*ardr[132];
   ardr[134]=ardr[16]*ardr[133];
   ardr[132]=1./2.*ardr[134] + ardr[132];
   ardr[132]=ardr[15]*ardr[132];
   ardr[134]=183./32.*ardr[13] - 9./16.*ardr[40] + 21./4.*ardr[27] - 5./
   2.*ardr[6] + 17./32.*ardr[22] + 7./4.*ardr[7] - ardr[24] - 15./4.*
   ardr[26];
   ardr[135]= - 3*ardr[14];
   ardr[136]=ardr[8]*ardr[133];
   ardr[137]=11*ardr[9] + 39./8. + ardr[109];
   ardr[137]=1./2.*ardr[137] - 9*ardr[19];
   ardr[137]=ardr[11]*ardr[137];
   ardr[131]=1./2.*ardr[131] - ardr[19];
   ardr[138]=MMH*ardr[131];
   ardr[61]=3*ardr[132] + ardr[61] + 1./8.*ardr[138] + 1./2.*ardr[119]
    + 1./4.*ardr[128] + 1./4.*ardr[137] + 3./4.*ardr[136] + 1./2.*
   ardr[134] + ardr[135];
   ardr[61]=ardr[15]*ardr[61];
   ardr[119]= - 1 - ardr[10];
   ardr[128]=5*ardr[9];
   ardr[119]=3*ardr[119] + ardr[128];
   ardr[132]=ardr[20]*ardr[119];
   ardr[119]=ardr[1]*ardr[119];
   ardr[119]=3*ardr[132] + ardr[119];
   ardr[119]=ardr[11]*ardr[119];
   ardr[132]=3*ardr[63] + ardr[62];
   ardr[134]=ardr[8]*ardr[132];
   ardr[136]= - ardr[9] + 3 - ardr[10];
   ardr[137]=ardr[20]*ardr[136];
   ardr[136]=ardr[1]*ardr[136];
   ardr[136]=3*ardr[137] + ardr[136];
   ardr[136]=ardr[12]*ardr[136];
   ardr[137]=ardr[91] + 1./3.*ardr[90];
   ardr[137]=MMH*ardr[137];
   ardr[138]= - ardr[14] + ardr[13] + ardr[7] - ardr[6];
   ardr[139]=ardr[20]*ardr[138];
   ardr[138]=ardr[1]*ardr[138];
   ardr[119]=1./4.*ardr[137] + 1./2.*ardr[136] + 1./2.*ardr[119] + 1./2.
   *ardr[134] + 3*ardr[139] + ardr[138];
   ardr[134]= - 1./2.*ardr[24];
   ardr[136]=ardr[27] + ardr[123] + ardr[67] + ardr[122] + ardr[134] - 
   ardr[26];
   ardr[137]=ardr[11]*ardr[133];
   ardr[138]=ardr[12]*ardr[131];
   ardr[136]=ardr[138] + 1./2.*ardr[136] + ardr[137];
   ardr[136]=ardr[15]*ardr[136];
   ardr[138]=ardr[49] + ardr[47];
   ardr[138]=1./2.*ardr[138] - ardr[48];
   ardr[139]=ardr[34] + 33./4.*ardr[114] - ardr[31];
   ardr[140]=ardr[11]*ardr[139];
   ardr[141]= - ardr[34] + 33./4.*ardr[120] + ardr[31];
   ardr[141]=ardr[12]*ardr[141];
   ardr[138]=ardr[141] + 99./8.*ardr[138] + ardr[140];
   ardr[138]=ardr[41]*ardr[138];
   ardr[140]=ardr[11]*ardr[132];
   ardr[90]=3*ardr[91] + ardr[90];
   ardr[91]=ardr[12]*ardr[90];
   ardr[91]=3*ardr[136] + 3*ardr[138] + ardr[140] + ardr[91];
   ardr[136]=pow(ardr[5],2);
   ardr[91]=ardr[136]*ardr[91];
   ardr[61]=1./2.*ardr[91] + ardr[61] + 1./2.*ardr[119] + ardr[88];
   ardr[61]=ardr[136]*ardr[61];
   ardr[88]=541 - 135*ardr[29];
   ardr[88]= - 15./2.*ardr[32] + 71*ardr[37] + 1./2.*ardr[88] - 15*
   ardr[36];
   ardr[91]= - 9./4.*ardr[34];
   ardr[119]=5./2.*ardr[110];
   ardr[138]=5./2.*ardr[111];
   ardr[88]=ardr[138] + ardr[119] + ardr[91] + 37./8.*ardr[31] + 1./8.*
   ardr[88] - 13*ardr[35];
   ardr[88]=ardr[12]*ardr[88];
   ardr[141]= - 3./2.*ardr[32];
   ardr[142]= - 13./4.*ardr[34];
   ardr[143]= - 3*ardr[36];
   ardr[144]=ardr[142] - 13./4.*ardr[31] + 165./4.*ardr[35] + ardr[141]
    - 165./4.*ardr[37] - 139./4. + ardr[143];
   ardr[144]=ardr[8]*ardr[144];
   ardr[145]= - 27*ardr[29];
   ardr[146]=49 + ardr[145];
   ardr[59]=ardr[60] - 13./2.*ardr[32] + ardr[59] + 1./2.*ardr[146] - 
   13*ardr[36];
   ardr[59]=19./8.*ardr[34] + 1./8.*ardr[59] + ardr[31];
   ardr[59]=MMH*ardr[59];
   ardr[60]= - 3*ardr[49] + 11*ardr[45] - 9./2.*ardr[42];
   ardr[146]=1./4.*ardr[60];
   ardr[147]=1./2.*ardr[53];
   ardr[148]=3*ardr[48] + 3*ardr[51] - 151./8.*ardr[47] + ardr[147] + 
   ardr[146] - 3*ardr[52];
   ardr[149]=985./2. + ardr[145];
   ardr[149]=1./2.*ardr[149] + ardr[143];
   ardr[150]=3./4.*ardr[34];
   ardr[149]=ardr[150] - 1./2.*ardr[31] - 259./16.*ardr[35] - 3./16.*
   ardr[32] + 1./8.*ardr[149] + 43*ardr[37];
   ardr[149]=ardr[11]*ardr[149];
   ardr[151]= - ardr[8] + 1./4.*ardr[11];
   ardr[151]=ardr[11]*ardr[151];
   ardr[152]=pow(ardr[8],2);
   ardr[151]= - 7./4.*ardr[152] + 9*ardr[151];
   ardr[153]= - 5./8.*ardr[12] - ardr[8] + 15./4.*ardr[11];
   ardr[153]=ardr[12]*ardr[153];
   ardr[151]=1./4.*ardr[151] + ardr[153];
   ardr[151]=ardr[56]*ardr[151];
   ardr[153]=23./8.*ardr[8] + 2*ardr[11];
   ardr[153]=ardr[11]*ardr[153];
   ardr[154]=23*ardr[8];
   ardr[155]=ardr[154] - 49./2.*ardr[11];
   ardr[155]=1./4.*ardr[155] - 18*ardr[12];
   ardr[155]=ardr[12]*ardr[155];
   ardr[153]=ardr[153] + ardr[155];
   ardr[153]=ardr[3]*ardr[153];
   ardr[59]=ardr[153] + ardr[151] + 1./2.*ardr[59] + 1./2.*ardr[88] + 
   ardr[149] + 1./4.*ardr[144] + 25./4.*ardr[14] + 151./8.*ardr[13] - 
   13./4.*ardr[39] - 7./4.*ardr[43] - 1./8.*ardr[6] + 21./8.*ardr[46]
    + 3./16.*ardr[7] - 23./16.*ardr[44] + 1./2.*ardr[148] + ardr[50];
   ardr[59]=ardr[41]*ardr[59];
   ardr[88]=ardr[8]*ardr[90];
   ardr[144]= - 7*ardr[14] - 29*ardr[13] + 7*ardr[7] + 29*ardr[6];
   ardr[144]=ardr[20]*ardr[144];
   ardr[148]= - 5*ardr[13];
   ardr[149]=ardr[14] + ardr[148] - ardr[7] + 5*ardr[6];
   ardr[149]=ardr[1]*ardr[149];
   ardr[144]=ardr[88] + 1./3.*ardr[144] + ardr[149];
   ardr[149]= - 29./6.*ardr[9] + 5./9. + 3./4.*ardr[10];
   ardr[149]=ardr[20]*ardr[149];
   ardr[151]=1./4.*ardr[10];
   ardr[153]= - 5./2.*ardr[9] + 1./3. + ardr[151];
   ardr[153]=ardr[1]*ardr[153];
   ardr[149]=ardr[149] + ardr[153];
   ardr[149]=ardr[11]*ardr[149];
   ardr[153]= - 5./3. + 13./4.*ardr[9];
   ardr[153]=ardr[20]*ardr[153];
   ardr[155]= - 1./3. + 5./4.*ardr[9];
   ardr[155]=ardr[1]*ardr[155];
   ardr[153]=1./3.*ardr[153] + ardr[155];
   ardr[153]=ardr[12]*ardr[153];
   ardr[62]=ardr[63] + 1./3.*ardr[62];
   ardr[62]=MMH*ardr[62];
   ardr[59]=ardr[61] + ardr[65] + ardr[59] + 1./8.*ardr[62] + ardr[153]
    + 1./4.*ardr[144] + ardr[149];
   ardr[59]=ardr[136]*ardr[59];
   ardr[61]= - 3./2.*ardr[26];
   ardr[62]=ardr[61] + 3*ardr[28] + ardr[134];
   ardr[62]=3./2.*ardr[27] + 3*ardr[62] + 11./2.*ardr[22];
   ardr[63]= - 7*ardr[40];
   ardr[62]=ardr[68] + 1./2.*ardr[62] + ardr[63];
   ardr[62]=EPAIR2*ardr[62];
   ardr[65]=17*ardr[18];
   ardr[134]=ardr[128] - 53./9. + ardr[65];
   ardr[134]=ardr[12]*ardr[134];
   ardr[144]= - 5./3.*ardr[28] - 3./4.*ardr[25];
   ardr[149]= - 7./4.*ardr[26];
   ardr[153]= - 5./8.*ardr[9] + 5./3. - 17./8.*ardr[18];
   ardr[153]=1./3.*ardr[153] - 3./4.*ardr[19];
   ardr[155]=ardr[8]*ardr[153];
   ardr[156]=3./2.*ardr[19];
   ardr[157]=ardr[156] - 1./12.*ardr[9] + 23./12.*ardr[18] + 181./288.
    + EPAIR2;
   ardr[157]=ardr[11]*ardr[157];
   ardr[62]=1./8.*ardr[134] + 1./2.*ardr[157] + ardr[155] + ardr[108]
    + ardr[62] - 11./64.*ardr[13] + 9./32.*ardr[40] + ardr[112] + 
   ardr[123] + 91./64.*ardr[22] + ardr[106] + ardr[149] + ardr[105] + 1.
   /2.*ardr[144] + ardr[104];
   ardr[104]= - 3*ardr[38];
   ardr[105]= - 22*ardr[33];
   ardr[106]= - 6*ardr[36];
   ardr[108]= - 6*ardr[32];
   ardr[134]= - 26*ardr[35];
   ardr[144]=6*ardr[31];
   ardr[155]=6*ardr[34];
   ardr[157]=11./3.*ardr[9] - 20./9. + 3*ardr[10];
   ardr[157]=2*ardr[20]*ardr[157];
   ardr[158]=3*ardr[9];
   ardr[159]=ardr[158] - 4./3. + ardr[10];
   ardr[159]=2*ardr[1]*ardr[159];
   ardr[160]=ardr[159] + ardr[157] + ardr[155] + ardr[144] + ardr[134]
    + ardr[108] - 47./2.*ardr[37] + ardr[106] + ardr[105] - 157./2. + 
   ardr[104];
   ardr[160]=ardr[16]*ardr[160];
   ardr[161]=ardr[156] - 139./6.*ardr[18] - 91./6. + ardr[70];
   ardr[161]=ardr[11]*ardr[161];
   ardr[161]=1./2.*ardr[161] + ardr[135] - 11*ardr[13] - 25*ardr[40] + 
   11*ardr[22] + 3*ardr[27];
   ardr[162]= - 155./6. + ardr[70];
   ardr[162]=ardr[73] + 1./2.*ardr[162] - 32./3.*ardr[18];
   ardr[162]=ardr[12]*ardr[162];
   ardr[160]=ardr[160] + ardr[161] + ardr[162];
   ardr[160]=ardr[3]*ardr[160];
   ardr[162]= - 1067./6. + ardr[87];
   ardr[163]=33*ardr[33];
   ardr[164]=47./2.*ardr[37];
   ardr[162]=ardr[70] + ardr[164] + 1./2.*ardr[162] + ardr[163];
   ardr[69]=ardr[71] - 8./3.*ardr[18] - 8./3. + ardr[69];
   ardr[69]=ardr[12]*ardr[69];
   ardr[72]=23./4. + ardr[72];
   ardr[72]=ardr[8]*ardr[72];
   ardr[69]=ardr[79] + ardr[69] + ardr[77] + 1./3.*ardr[72];
   ardr[69]=ardr[56]*ardr[69];
   ardr[72]=ardr[81] + ardr[84];
   ardr[72]=ardr[3]*ardr[72];
   ardr[72]=18*ardr[72] - 5./2.*ardr[19] - 11*ardr[33] + 133./3. - 3./2.
   *ardr[38];
   ardr[72]=ardr[3]*ardr[72];
   ardr[72]=3./2.*ardr[75] + ardr[72];
   ardr[72]=MMt*ardr[72];
   ardr[75]=7*EPAIR2;
   ardr[77]=13./4.*ardr[35];
   ardr[79]=11./16.*ardr[19];
   ardr[84]= - 11./12.*ardr[9] + 5./9. - 3./4.*ardr[10];
   ardr[84]=ardr[20]*ardr[84];
   ardr[165]= - 3./4.*ardr[9];
   ardr[166]=ardr[165] + 1./3. - 1./4.*ardr[10];
   ardr[166]=ardr[1]*ardr[166];
   ardr[69]=ardr[72] + ardr[160] + ardr[69] + ardr[166] + ardr[84] + 
   ardr[79] + ardr[92] + ardr[89] + 71./16.*ardr[18] + ardr[77] + 1./8.
   *ardr[162] + ardr[75];
   ardr[69]=MMt*ardr[69];
   ardr[160]=3*ardr[38];
   ardr[162]=1595./12. + ardr[160];
   ardr[167]=11*ardr[33];
   ardr[162]=ardr[164] + 1./2.*ardr[162] + ardr[167];
   ardr[164]=1./16.*ardr[54] - 1./3.*ardr[23];
   ardr[168]=EPAIR2*ardr[23];
   ardr[164]=5*ardr[164] + 2*ardr[168];
   ardr[164]=ardr[16]*ardr[164];
   ardr[169]= - 10*EPAIR2;
   ardr[170]=13./2.*ardr[35];
   ardr[171]= - 3./2.*ardr[34];
   ardr[172]= - 3./2.*ardr[10];
   ardr[173]= - 11./6.*ardr[9] + 10./9. + ardr[172];
   ardr[173]=ardr[20]*ardr[173];
   ardr[174]= - 1./2.*ardr[10];
   ardr[175]= - 3./2.*ardr[9] + 2./3. + ardr[174];
   ardr[175]=ardr[1]*ardr[175];
   ardr[168]=ardr[168] - 5./8.*ardr[54] + ardr[23];
   ardr[168]=1./2.*ardr[11]*ardr[168];
   ardr[95]=3./2.*ardr[95];
   ardr[176]= - 3./2.*ardr[31];
   ardr[162]=ardr[164] + ardr[95] + ardr[168] + ardr[175] + ardr[173]
    + ardr[171] + ardr[176] + ardr[170] + 1./4.*ardr[162] + ardr[169];
   ardr[162]=ardr[16]*ardr[162];
   ardr[177]= - 5./4.*ardr[9] + 10./3. - 17./4.*ardr[18];
   ardr[177]=1./3.*ardr[177] + ardr[73];
   ardr[177]=ardr[16]*ardr[177];
   ardr[178]=ardr[128] - 40./3. + ardr[65];
   ardr[98]=1./3.*ardr[178] + ardr[98];
   ardr[98]=ardr[3]*ardr[16]*ardr[98];
   ardr[98]=ardr[153] + ardr[98];
   ardr[98]=MMt*ardr[98];
   ardr[98]=ardr[177] + ardr[98];
   ardr[153]=ardr[15]*ardr[98];
   ardr[177]=5./8.*ardr[9] - 5./3. + 17./8.*ardr[18];
   ardr[177]=1./9.*ardr[177] + ardr[80];
   ardr[177]=1./2.*MMH*ardr[177];
   ardr[178]= - 11*ardr[11] - 19*ardr[12];
   ardr[179]= - 15*ardr[16];
   ardr[178]=1./6.*ardr[178] + ardr[179];
   ardr[178]=ardr[3]*ardr[16]*ardr[178];
   ardr[180]=19*ardr[8] - 37*ardr[12];
   ardr[180]=ardr[56]*ardr[16]*ardr[180];
   ardr[69]=ardr[153] + ardr[69] + ardr[178] + 1./24.*ardr[180] + 
   ardr[177] + ardr[62] + ardr[162];
   ardr[69]=ardr[15]*ardr[69];
   ardr[153]= - 527./3. + ardr[145];
   ardr[162]= - 25./2.*ardr[35];
   ardr[178]= - 9*ardr[31];
   ardr[153]=ardr[113] + ardr[178] + ardr[162] + ardr[141] - 209./4.*
   ardr[37] + 1./2.*ardr[153] + ardr[143];
   ardr[180]= - ardr[8]*ardr[58];
   ardr[181]=ardr[11]*ardr[58];
   ardr[153]=ardr[181] + 1./4.*ardr[153] + ardr[180];
   ardr[153]=ardr[11]*ardr[153];
   ardr[182]= - 27./16.*ardr[29];
   ardr[183]= - 13./8.*ardr[36];
   ardr[184]= - 23./16.*ardr[32];
   ardr[185]= - 13./12.*ardr[35];
   ardr[186]=11./6.*ardr[31];
   ardr[187]=17./8.*ardr[34];
   ardr[188]=ardr[187] + ardr[186] + ardr[185] + ardr[184] - 47./48.*
   ardr[37] + ardr[183] - 7./3. + ardr[182];
   ardr[188]=MMH*ardr[188];
   ardr[189]= - 3*ardr[32];
   ardr[190]=13*ardr[35];
   ardr[191]= - 3*ardr[31];
   ardr[192]=ardr[142] + ardr[191] + ardr[190] + ardr[189] + 47./4.*
   ardr[37] + 115./6. + ardr[143];
   ardr[192]=ardr[8]*ardr[192];
   ardr[60]=ardr[60] - 55*ardr[52];
   ardr[60]=ardr[192] + 7*ardr[14] + 265./4.*ardr[13] - 65./4.*ardr[39]
    - 45./4.*ardr[43] - 3*ardr[6] + 11*ardr[46] + 63./4.*ardr[7] - 23./
   4.*ardr[44] + 11./2.*ardr[50] + 15./4.*ardr[48] + 27./2.*ardr[51] - 
   67./4.*ardr[47] + 1./2.*ardr[60] - ardr[53];
   ardr[192]= - 3*ardr[152];
   ardr[193]= - 41*ardr[8] + 17./2.*ardr[11];
   ardr[193]=ardr[11]*ardr[193];
   ardr[194]= - ardr[12] - 49*ardr[8] + 87./2.*ardr[11];
   ardr[194]=ardr[12]*ardr[194];
   ardr[193]=ardr[194] + ardr[192] + ardr[193];
   ardr[193]=ardr[56]*ardr[193];
   ardr[110]=5./4.*ardr[111] + 5./4.*ardr[110] - 9./8.*ardr[34] + 9./4.
   *ardr[31] - 59./16.*ardr[35] - 15./32.*ardr[32] - 45./16.*ardr[37]
    - 15./16.*ardr[36] + 47./3. - 135./32.*ardr[29];
   ardr[110]=ardr[12]*ardr[110];
   ardr[194]=ardr[154] - 67./2.*ardr[11];
   ardr[194]=1./4.*ardr[194];
   ardr[107]=ardr[194] + ardr[107];
   ardr[107]=ardr[12]*ardr[107];
   ardr[195]=5*ardr[8] - 3./2.*ardr[11];
   ardr[195]=3./4.*ardr[11]*ardr[195];
   ardr[107]=ardr[195] + ardr[107];
   ardr[107]=ardr[3]*ardr[107];
   ardr[60]=ardr[107] + 1./8.*ardr[193] + 1./2.*ardr[188] + ardr[110]
    + 1./4.*ardr[60] + ardr[153];
   ardr[60]=ardr[41]*ardr[60];
   ardr[107]=ardr[7] - ardr[14];
   ardr[107]=ardr[20]*ardr[107];
   ardr[110]= - ardr[7] + ardr[14];
   ardr[110]=ardr[1]*ardr[110];
   ardr[107]=7./3.*ardr[107] + ardr[110];
   ardr[110]=ardr[84] + ardr[166];
   ardr[110]=ardr[8]*ardr[110];
   ardr[153]=1./2.*ardr[10];
   ardr[188]= - 3*EPAIR2;
   ardr[193]=3./2.*ardr[9] + ardr[153] + 5./6. + ardr[188];
   ardr[193]=ardr[1]*ardr[193];
   ardr[196]=3*EPAIR2;
   ardr[197]=11./6.*ardr[9] + 3./2.*ardr[10] + 13./18. + ardr[196];
   ardr[197]=ardr[20]*ardr[197];
   ardr[193]=ardr[197] + ardr[193];
   ardr[193]=ardr[11]*ardr[193];
   ardr[197]= - 10./9. + 11./4.*ardr[9];
   ardr[197]=ardr[20]*ardr[197];
   ardr[198]= - 2./3. + 9./4.*ardr[9];
   ardr[198]=ardr[1]*ardr[198];
   ardr[197]=ardr[197] + ardr[198];
   ardr[197]=ardr[12]*ardr[197];
   ardr[198]=11./36.*ardr[9] - 5./27. + ardr[151];
   ardr[198]=ardr[20]*ardr[198];
   ardr[151]= - 1./3. + ardr[151];
   ardr[151]=1./3.*ardr[151] + 1./4.*ardr[9];
   ardr[151]=ardr[1]*ardr[151];
   ardr[151]=ardr[198] + ardr[151];
   ardr[151]=MMH*ardr[151];
   ardr[107]=1./2.*ardr[151] + ardr[197] + 1./2.*ardr[193] + 1./4.*
   ardr[107] + ardr[110];
   ardr[59]=ardr[59] + ardr[69] + ardr[107] + ardr[60];
   ardr[59]=ardr[136]*ardr[59];
   ardr[60]=ardr[159] + ardr[157] + ardr[155] + ardr[144] + ardr[134]
    + ardr[108] + 1./2.*ardr[37] + ardr[106] + ardr[105] - 109./2. + 
   ardr[104];
   ardr[60]=ardr[16]*ardr[60];
   ardr[69]= - 1./2. - ardr[17];
   ardr[69]=3*ardr[69] - ardr[19];
   ardr[69]=ardr[12]*ardr[69];
   ardr[60]=ardr[60] + ardr[161] + 3./2.*ardr[69];
   ardr[60]=ardr[3]*ardr[60];
   ardr[69]= - 4321./18. + ardr[87];
   ardr[87]= - 1./2.*ardr[37];
   ardr[69]=ardr[70] + ardr[87] + 1./2.*ardr[69] + ardr[163];
   ardr[104]=ardr[78] - ardr[19];
   ardr[105]=1./4.*ardr[12]*ardr[104];
   ardr[76]= - 5./4.*ardr[16] + ardr[105] + ardr[76] - 1./2.*ardr[8];
   ardr[76]=ardr[56]*ardr[76];
   ardr[60]=ardr[72] + ardr[60] + 3./2.*ardr[76] + ardr[166] + ardr[84]
    + ardr[79] + ardr[92] + ardr[89] + 511./144.*ardr[18] + ardr[77] + 
   1./8.*ardr[69] + ardr[75];
   ardr[60]=MMt*ardr[60];
   ardr[69]=2801./36. + ardr[160];
   ardr[69]=ardr[87] + 1./2.*ardr[69] + ardr[167];
   ardr[69]=ardr[164] + ardr[95] + ardr[168] + ardr[175] + ardr[173] + 
   ardr[171] + ardr[176] + ardr[170] + 1./4.*ardr[69] + ardr[169];
   ardr[69]=ardr[16]*ardr[69];
   ardr[70]=5./6. + ardr[70];
   ardr[72]=7./3.*ardr[18];
   ardr[70]=ardr[156] + 1./2.*ardr[70] + ardr[72];
   ardr[70]=ardr[11]*ardr[70];
   ardr[75]=7*ardr[33];
   ardr[76]=17./4. + ardr[75];
   ardr[77]=17./2.*ardr[35];
   ardr[76]=ardr[77] + 1./3.*ardr[76] + ardr[189];
   ardr[76]=ardr[16]*ardr[76];
   ardr[79]= - 1./2.*ardr[22];
   ardr[84]=1./2.*ardr[13] + ardr[79] + ardr[40];
   ardr[70]=ardr[76] + 7./3.*ardr[84] + 1./2.*ardr[70];
   ardr[70]=ardr[3]*ardr[70];
   ardr[76]= - 7*ardr[33];
   ardr[84]= - 17*ardr[35];
   ardr[89]= - ardr[19] - 49./9.*ardr[18] + ardr[84] + 515./36. + 
   ardr[76];
   ardr[75]=11 + ardr[75];
   ardr[75]=1./3.*ardr[75] + ardr[19];
   ardr[81]=ardr[3]*ardr[81];
   ardr[75]=1./2.*ardr[75] + 9*ardr[81];
   ardr[75]=MMt*ardr[3]*ardr[75];
   ardr[70]=ardr[75] + 1./16.*ardr[89] + ardr[70];
   ardr[70]=MMt*ardr[70];
   ardr[75]= - 215./12. + ardr[76];
   ardr[75]=1./3.*ardr[75] + ardr[84];
   ardr[76]=25./4.*ardr[54] - 17./3.*ardr[23];
   ardr[76]=ardr[11]*ardr[76];
   ardr[81]= - 25./2.*ardr[54] + 31./3.*ardr[23];
   ardr[81]=ardr[16]*ardr[81];
   ardr[75]=1./2.*ardr[81] + 1./2.*ardr[75] + ardr[76];
   ardr[75]=ardr[16]*ardr[75];
   ardr[76]= - 103*ardr[26] + 439./4.*ardr[22];
   ardr[81]= - 5*ardr[9];
   ardr[84]=ardr[81] - 167./8. - 17*ardr[18];
   ardr[84]=ardr[11]*ardr[84];
   ardr[76]=ardr[84] - 115./8.*ardr[13] - 233./4.*ardr[40] + 1./2.*
   ardr[76] + 11*ardr[6];
   ardr[75]=1./6.*ardr[76] + ardr[75];
   ardr[76]=77./2.*ardr[11] - 17*ardr[16];
   ardr[76]=ardr[3]*ardr[16]*ardr[76];
   ardr[75]=1./2.*ardr[75] + 1./3.*ardr[76];
   ardr[76]=1 + ardr[35];
   ardr[76]=ardr[3]*ardr[16]*ardr[76];
   ardr[84]=MMt*ardr[3]*ardr[19];
   ardr[76]=ardr[84] - 1./8.*ardr[35] + ardr[76];
   ardr[76]=MMt*ardr[76];
   ardr[84]= - ardr[16]*ardr[35];
   ardr[76]=1./4.*ardr[84] + ardr[76];
   ardr[84]=pow(ardr[2],2);
   ardr[76]=ardr[84]*ardr[76];
   ardr[70]=1./2.*ardr[76] + 1./2.*ardr[75] + ardr[70];
   ardr[70]=ardr[84]*ardr[70];
   ardr[75]=3*ardr[12];
   ardr[76]= - 5*ardr[8];
   ardr[89]=ardr[76] + ardr[75];
   ardr[89]=ardr[56]*ardr[16]*ardr[89];
   ardr[95]= - 11./3.*ardr[11] + 15*ardr[12];
   ardr[95]=1./2.*ardr[95] + ardr[179];
   ardr[95]=ardr[3]*ardr[16]*ardr[95];
   ardr[60]=ardr[70] + ardr[60] + ardr[95] + 3./8.*ardr[89] + ardr[177]
    + ardr[62] + ardr[69];
   ardr[60]=ardr[84]*ardr[60];
   ardr[62]=ardr[15]*ardr[84]*ardr[98];
   ardr[60]=ardr[60] + ardr[62];
   ardr[60]=ardr[15]*ardr[60];
   ardr[62]=9*ardr[30];
   ardr[69]=3./4.*ardr[38];
   ardr[70]=27./2.*ardr[29];
   ardr[89]=3*ardr[36];
   ardr[95]=3./2.*ardr[32];
   ardr[98]=ardr[95] + ardr[89] + ardr[70] + ardr[69] + 197./8. + 
   ardr[62];
   ardr[98]=1./2.*ardr[98] + ardr[127];
   ardr[98]=ardr[16]*ardr[98];
   ardr[106]=15*ardr[30];
   ardr[108]=9./4.*ardr[38];
   ardr[110]= - 15*ardr[17];
   ardr[134]=11./4.*ardr[19] + ardr[110] + ardr[108] - 721./16. + 
   ardr[106];
   ardr[151]=1./2. + ardr[17];
   ardr[155]=ardr[8]*ardr[151];
   ardr[156]= - 1./2.*ardr[27] + 1./2.*ardr[25] - ardr[21];
   ardr[156]=1./4.*ardr[14] + 1./4.*ardr[39] + 1./2.*ardr[156] + 
   ardr[40];
   ardr[155]=ardr[156] + ardr[155];
   ardr[157]=3*ardr[16];
   ardr[159]= - ardr[12]*ardr[19];
   ardr[155]=ardr[157] + 3*ardr[155] + 1./4.*ardr[159];
   ardr[155]=ardr[56]*ardr[155];
   ardr[134]=1./2.*ardr[134] + 3*ardr[155];
   ardr[155]= - 2 + ardr[78];
   ardr[155]=ardr[8]*ardr[155];
   ardr[63]=ardr[155] - ardr[14] - 3*ardr[39] + ardr[63] + ardr[85] + 
   ardr[27];
   ardr[85]= - 1 - ardr[19];
   ardr[85]=ardr[12]*ardr[85];
   ardr[155]=ardr[11]*ardr[19];
   ardr[161]= - 4*ardr[30];
   ardr[163]= - ardr[38] - 9./2. + ardr[161];
   ardr[163]=ardr[16]*ardr[163];
   ardr[85]=ardr[163] + 1./2.*ardr[85] + ardr[63] + 1./4.*ardr[155];
   ardr[85]=ardr[3]*ardr[85];
   ardr[164]= - 1./2.*ardr[38] + 9./2. - 2*ardr[30];
   ardr[164]=ardr[3]*ardr[164];
   ardr[164]= - 3./2.*ardr[56] + ardr[164];
   ardr[164]=3*MMt*ardr[164];
   ardr[85]=ardr[164] + 1./2.*ardr[134] + 3*ardr[85];
   ardr[85]=MMt*ardr[85];
   ardr[134]= - 9./4.*ardr[30];
   ardr[166]=15./4.*ardr[17];
   ardr[167]=ardr[82] + ardr[34] + ardr[166] - 5 + ardr[134];
   ardr[167]=MMH*ardr[167];
   ardr[168]= - 7./24.*ardr[18];
   ardr[169]=ardr[71] - 3 + ardr[168];
   ardr[169]=ardr[11]*ardr[169];
   ardr[170]= - ardr[28] + ardr[27];
   ardr[171]=EPAIR2*ardr[170];
   ardr[177]= - 3./4.*EPAIR2;
   ardr[179]=11./3. + ardr[177];
   ardr[179]=ardr[14]*ardr[179];
   ardr[193]=ardr[80] + 1 + 3./4.*ardr[17];
   ardr[193]=ardr[8]*ardr[193];
   ardr[197]= - 3./4.*ardr[16] + 21./8.*ardr[8] - 2*ardr[12];
   ardr[197]=ardr[56]*ardr[16]*ardr[197];
   ardr[198]= - ardr[8] + 3./2.*ardr[11];
   ardr[199]=3*ardr[198] + 13*ardr[12];
   ardr[199]=ardr[3]*ardr[16]*ardr[199];
   ardr[200]= - 2*ardr[28];
   ardr[85]=ardr[85] + 5./2.*ardr[199] + 3*ardr[197] + 1./4.*ardr[167]
    + ardr[98] - 157./24.*ardr[12] + ardr[169] + 3./2.*ardr[193] + 
   ardr[179] + 3./4.*ardr[171] + 1./12.*ardr[13] + 45./16.*ardr[39] + 
   16*ardr[40] + ardr[112] - 139./48.*ardr[22] + 45./16.*ardr[26] - 9./
   8.*ardr[21] - 41./12.*ardr[24] + ardr[200] - 27./16.*ardr[25];
   ardr[85]=MMt*ardr[85];
   ardr[98]=ardr[3]*ardr[155];
   ardr[82]=ardr[82] + 3*ardr[98];
   ardr[82]=MMt*ardr[82];
   ardr[98]= - ardr[11] + ardr[12];
   ardr[155]=ardr[98] - ardr[16];
   ardr[167]=5./2.*ardr[11] - ardr[12];
   ardr[167]=ardr[3]*ardr[16]*ardr[167];
   ardr[82]=1./2.*ardr[82] + 1./8.*ardr[155] + ardr[167];
   ardr[82]=MMt*ardr[82];
   ardr[98]=ardr[16]*ardr[98];
   ardr[82]=1./4.*ardr[98] + ardr[82];
   ardr[82]=ardr[84]*ardr[82];
   ardr[155]= - 1./16.*ardr[8] - ardr[11];
   ardr[167]= - 119./3. + ardr[196];
   ardr[167]=ardr[12]*ardr[167];
   ardr[155]= - 77./16.*ardr[16] + 3*ardr[155] + 1./4.*ardr[167];
   ardr[155]=ardr[16]*ardr[155];
   ardr[167]= - 3./4.*ardr[30];
   ardr[169]=ardr[167] + ardr[34];
   ardr[169]=ardr[16]*ardr[169];
   ardr[171]=ardr[40] + ardr[25] - ardr[21];
   ardr[179]=9./8.*ardr[171];
   ardr[169]=ardr[179] + ardr[169];
   ardr[169]=MMH*ardr[169];
   ardr[82]=1./2.*ardr[82] + ardr[85] + ardr[155] + 1./2.*ardr[169];
   ardr[82]=ardr[84]*ardr[82];
   ardr[85]= - 2*ardr[31];
   ardr[62]=ardr[127] + ardr[85] + ardr[95] + ardr[89] + ardr[70] + 
   ardr[69] - 173./24. + ardr[62];
   ardr[62]=ardr[16]*ardr[62];
   ardr[69]=1./2.*ardr[34];
   ardr[70]=1./2.*ardr[31];
   ardr[127]= - 1./8.*ardr[19];
   ardr[155]=ardr[127] + ardr[69] + ardr[70] + 7./72.*ardr[18] + 
   ardr[166] - 43./9. + ardr[134];
   ardr[155]=MMH*ardr[155];
   ardr[151]=3*ardr[151];
   ardr[169]=1./8.*ardr[19];
   ardr[193]=ardr[151] + ardr[169];
   ardr[193]=ardr[8]*ardr[193];
   ardr[156]=3*ardr[156];
   ardr[193]=ardr[157] + 3./8.*ardr[159] + ardr[156] + ardr[193];
   ardr[193]=ardr[56]*ardr[193];
   ardr[197]=ardr[110] + ardr[108] + 773./48. + ardr[106];
   ardr[199]= - 1 + ardr[73];
   ardr[199]=ardr[12]*ardr[199];
   ardr[201]= - 2*ardr[38] - 9 - 8*ardr[30];
   ardr[201]=ardr[16]*ardr[201];
   ardr[199]=ardr[201] + 2*ardr[63] + ardr[199];
   ardr[199]=ardr[3]*ardr[199];
   ardr[161]= - ardr[38] + 9 + ardr[161];
   ardr[161]=ardr[3]*ardr[161];
   ardr[161]= - 3*ardr[56] + ardr[161];
   ardr[161]=MMt*ardr[161];
   ardr[161]=3*ardr[161] + 3*ardr[199] + 3*ardr[193] + 1./2.*ardr[197]
    + 2*ardr[19];
   ardr[161]=MMt*ardr[161];
   ardr[193]=ardr[28] - 9./2.*ardr[25];
   ardr[197]= - 5*ardr[28];
   ardr[199]=ardr[197] + 3*ardr[26];
   ardr[199]= - 3./2.*ardr[13] + 1./2.*ardr[199] + ardr[27];
   ardr[199]=EPAIR2*ardr[199];
   ardr[201]= - 9./2.*ardr[21];
   ardr[202]=45./4.*ardr[39];
   ardr[193]=3*ardr[199] - 149./8.*ardr[13] + ardr[202] + 17./6.*
   ardr[40] - ardr[27] + 105./8.*ardr[22] + 11./2.*ardr[26] + ardr[201]
    + 3./2.*ardr[193] - 29./3.*ardr[24];
   ardr[199]=9./4.*ardr[17];
   ardr[168]=3./8.*ardr[19] + ardr[168] + 7./3. + ardr[199];
   ardr[168]=ardr[8]*ardr[168];
   ardr[203]=413./2. + 73*ardr[18];
   ardr[203]=1./3.*ardr[203] + ardr[96];
   ardr[203]=ardr[11]*ardr[203];
   ardr[204]=43*ardr[8] - 33*ardr[12];
   ardr[204]=1./4.*ardr[204] + ardr[103];
   ardr[204]=ardr[56]*ardr[16]*ardr[204];
   ardr[205]= - 15*ardr[8];
   ardr[101]=28*ardr[12] + ardr[205] + ardr[101];
   ardr[101]=ardr[3]*ardr[16]*ardr[101];
   ardr[206]=16./3. - 3./2.*EPAIR2;
   ardr[206]=ardr[14]*ardr[206];
   ardr[207]= - 7./3. + 1./2.*ardr[18];
   ardr[207]=7*ardr[207] + 9./2.*ardr[19];
   ardr[207]=ardr[12]*ardr[207];
   ardr[62]=ardr[161] + ardr[101] + 3./2.*ardr[204] + 1./2.*ardr[155]
    + ardr[62] + 1./4.*ardr[207] + 1./8.*ardr[203] + ardr[168] + 1./2.*
   ardr[193] + ardr[206];
   ardr[62]=MMt*ardr[62];
   ardr[69]=ardr[69] + ardr[70] + 1./9. + ardr[167];
   ardr[69]=ardr[16]*ardr[69];
   ardr[69]=ardr[179] + ardr[69];
   ardr[69]=MMH*ardr[69];
   ardr[101]=9*EPAIR2;
   ardr[155]=605./12. + ardr[101];
   ardr[155]=ardr[11]*ardr[155];
   ardr[155]= - 25./6.*ardr[8] + ardr[155];
   ardr[161]= - 167./12. + ardr[196];
   ardr[161]=ardr[12]*ardr[161];
   ardr[155]=227./24.*ardr[16] + 1./2.*ardr[155] + ardr[161];
   ardr[155]=ardr[16]*ardr[155];
   ardr[62]=ardr[62] + 1./2.*ardr[155] + ardr[69];
   ardr[69]=ardr[62] + ardr[82];
   ardr[69]=ardr[84]*ardr[69];
   ardr[82]= - 1./2. + ardr[19];
   ardr[155]= - ardr[19] + 1./2. + ardr[83];
   ardr[155]=ardr[3]*ardr[16]*ardr[155];
   ardr[161]=pow(ardr[3],2);
   ardr[167]= - MMt*ardr[161]*ardr[16]*ardr[17];
   ardr[82]=12*ardr[167] + 1./8.*ardr[82] + ardr[155];
   ardr[82]=MMt*ardr[82];
   ardr[155]=pow(ardr[16],2);
   ardr[168]= - ardr[3]*ardr[155];
   ardr[179]=ardr[16]*ardr[19];
   ardr[82]=ardr[82] + 1./4.*ardr[179] + ardr[168];
   ardr[82]=MMt*ardr[82];
   ardr[193]=1./4.*ardr[155];
   ardr[82]=ardr[193] + ardr[82];
   ardr[82]=ardr[84]*ardr[82];
   ardr[203]= - 25 - 7*ardr[18];
   ardr[204]=3*ardr[19];
   ardr[203]=1./3.*ardr[203] + ardr[204];
   ardr[72]=ardr[96] + ardr[72] + 25./3. + 18*ardr[17];
   ardr[72]=ardr[3]*ardr[16]*ardr[72];
   ardr[72]=72*ardr[167] + 1./8.*ardr[203] + ardr[72];
   ardr[72]=MMt*ardr[72];
   ardr[203]=ardr[93] - 2 - 7./12.*ardr[18];
   ardr[203]=ardr[16]*ardr[203];
   ardr[72]=ardr[72] + ardr[203] + 2./3.*ardr[168];
   ardr[72]=MMt*ardr[72];
   ardr[72]=1./6.*ardr[155] + ardr[72];
   ardr[82]=ardr[72] + 3*ardr[82];
   ardr[82]=ardr[15]*ardr[84]*ardr[82];
   ardr[69]=ardr[69] + ardr[82];
   ardr[69]=ardr[15]*ardr[69];
   ardr[72]=ardr[15]*ardr[72];
   ardr[62]=ardr[62] + ardr[72];
   ardr[62]=ardr[15]*ardr[62];
   ardr[72]=ardr[151] + ardr[80];
   ardr[72]=ardr[8]*ardr[72];
   ardr[72]=ardr[157] + 1./2.*ardr[159] + ardr[156] + ardr[72];
   ardr[72]=ardr[56]*ardr[72];
   ardr[82]=21./4.*ardr[19] + ardr[110] + ardr[108] + 287./16. + 
   ardr[106];
   ardr[72]=1./2.*ardr[82] + 3*ardr[72];
   ardr[82]= - 1./2. - ardr[19];
   ardr[106]=ardr[12]*ardr[82];
   ardr[108]= - ardr[11]*ardr[19];
   ardr[110]=1./4.*ardr[108];
   ardr[63]=ardr[163] + ardr[106] + ardr[63] + ardr[110];
   ardr[63]=ardr[3]*ardr[63];
   ardr[63]=ardr[164] + 1./2.*ardr[72] + 3*ardr[63];
   ardr[63]=MMt*ardr[63];
   ardr[72]=ardr[28] - 9./4.*ardr[25];
   ardr[106]=ardr[13] - ardr[26] + ardr[27];
   ardr[106]=EPAIR2*ardr[106];
   ardr[72]=3*ardr[106] - 155./12.*ardr[13] + ardr[202] + ardr[40] - 
   ardr[27] + 67./6.*ardr[22] + 7./4.*ardr[26] + ardr[201] + 3*ardr[72]
    - 17./3.*ardr[24];
   ardr[106]=19./4.*ardr[18];
   ardr[151]=ardr[106] + 7 + ardr[199];
   ardr[151]=ardr[8]*ardr[151];
   ardr[156]= - 77 - 43*ardr[18];
   ardr[157]=9*ardr[19];
   ardr[156]=1./3.*ardr[156] + ardr[157];
   ardr[156]=ardr[12]*ardr[156];
   ardr[163]=1./2.*ardr[32];
   ardr[86]=ardr[163] + ardr[36] + 9./2.*ardr[29] + ardr[86] - 17./8.
    + 3*ardr[30];
   ardr[86]=3./2.*ardr[86] + ardr[85];
   ardr[86]=ardr[16]*ardr[86];
   ardr[134]=ardr[31] - 19./12.*ardr[18] + ardr[166] - 19./3. + 
   ardr[134];
   ardr[134]=MMH*ardr[134];
   ardr[164]=5./3. + ardr[177];
   ardr[164]=ardr[14]*ardr[164];
   ardr[166]= - 239./16. - 22*ardr[18];
   ardr[166]=ardr[11]*ardr[166];
   ardr[103]=ardr[103] + 11*ardr[8] - 17./2.*ardr[12];
   ardr[103]=ardr[56]*ardr[16]*ardr[103];
   ardr[177]=39*ardr[12] + ardr[205] - 43./2.*ardr[11];
   ardr[177]=ardr[3]*ardr[16]*ardr[177];
   ardr[63]=ardr[63] + 1./2.*ardr[177] + 3./4.*ardr[103] + 1./4.*
   ardr[134] + ardr[86] + 1./8.*ardr[156] + 1./3.*ardr[166] + 1./2.*
   ardr[151] + 1./4.*ardr[72] + ardr[164];
   ardr[63]=MMt*ardr[63];
   ardr[72]=29./2. + ardr[109];
   ardr[86]= - 19*ardr[18];
   ardr[103]=ardr[86] - 29./2. + 9*ardr[17];
   ardr[103]=ardr[3]*ardr[16]*ardr[103];
   ardr[72]=36*ardr[167] + 1./8.*ardr[72] + ardr[103];
   ardr[72]=MMt*ardr[72];
   ardr[103]=6 + ardr[106];
   ardr[103]=ardr[16]*ardr[103];
   ardr[72]=ardr[72] + ardr[103] + 19*ardr[168];
   ardr[72]=MMt*ardr[72];
   ardr[72]=19./4.*ardr[155] + ardr[72];
   ardr[72]=ardr[15]*ardr[72];
   ardr[103]= - 59./12. + ardr[188];
   ardr[103]=ardr[11]*ardr[103];
   ardr[106]= - 137./6. + ardr[196];
   ardr[106]=ardr[12]*ardr[106];
   ardr[103]=23./2.*ardr[16] + ardr[106] + 29./4.*ardr[8] + ardr[103];
   ardr[103]=ardr[16]*ardr[103];
   ardr[70]=ardr[70] - 1./3. - 3./8.*ardr[30];
   ardr[70]=ardr[16]*ardr[70];
   ardr[70]=9./16.*ardr[171] + ardr[70];
   ardr[70]=MMH*ardr[70];
   ardr[63]=ardr[72] + ardr[63] + 1./4.*ardr[103] + ardr[70];
   ardr[63]=ardr[15]*ardr[63];
   ardr[70]=ardr[8]*ardr[19];
   ardr[72]=ardr[70] + ardr[159];
   ardr[72]=ardr[56]*ardr[72];
   ardr[72]=1./2.*ardr[72] + 1 + ardr[93];
   ardr[103]=1./2.*ardr[108];
   ardr[106]=ardr[103] + ardr[159];
   ardr[106]=ardr[3]*ardr[106];
   ardr[72]=1./2.*ardr[72] + ardr[106];
   ardr[72]=MMt*ardr[72];
   ardr[64]=ardr[64] - 9./16. + ardr[85];
   ardr[64]=ardr[16]*ardr[64];
   ardr[85]=ardr[8]*ardr[125];
   ardr[106]=ardr[204] + 3 + 29*ardr[18];
   ardr[106]=ardr[11]*ardr[106];
   ardr[109]=ardr[157] - 15./2. - 41*ardr[18];
   ardr[109]=ardr[12]*ardr[109];
   ardr[80]=ardr[80] - ardr[34] - 1./4.*ardr[18] + ardr[31];
   ardr[80]=MMH*ardr[80];
   ardr[134]=ardr[8] - ardr[12];
   ardr[134]=ardr[56]*ardr[16]*ardr[134];
   ardr[151]=13./2.*ardr[11] - 11*ardr[12];
   ardr[151]=ardr[3]*ardr[16]*ardr[151];
   ardr[64]=3./2.*ardr[72] + 1./2.*ardr[151] + 3./8.*ardr[134] + 1./4.*
   ardr[80] + ardr[64] + 1./8.*ardr[109] + 1./8.*ardr[106] + 3./8.*
   ardr[85] - 2*ardr[14] + 23./16.*ardr[13] - 3./4.*ardr[40] - 29./16.*
   ardr[22] + 2*ardr[24] + 3./8.*ardr[26];
   ardr[64]=MMt*ardr[64];
   ardr[72]=ardr[3]*ardr[16]*ardr[130];
   ardr[72]=1./8.*ardr[125] + ardr[72];
   ardr[72]=MMt*ardr[72];
   ardr[80]=ardr[16]*ardr[125];
   ardr[72]=1./4.*ardr[80] + ardr[72];
   ardr[72]=ardr[15]*MMt*ardr[72];
   ardr[80]=3./16.*ardr[16] + 25./16.*ardr[11] - ardr[12];
   ardr[80]=ardr[16]*ardr[80];
   ardr[85]=ardr[31] - ardr[34];
   ardr[106]=MMH*ardr[16]*ardr[85];
   ardr[64]=3*ardr[72] + ardr[64] + ardr[80] + 1./2.*ardr[106];
   ardr[64]=ardr[15]*ardr[64];
   ardr[72]=3./2.*ardr[34];
   ardr[80]=ardr[72] + 1./3. + ardr[176];
   ardr[80]=ardr[11]*ardr[80];
   ardr[106]=ardr[51] - ardr[50];
   ardr[109]= - 3*ardr[46];
   ardr[134]=ardr[8]*ardr[85];
   ardr[151]=17./4.*ardr[34] - 1./3. - 17./4.*ardr[31];
   ardr[151]=ardr[12]*ardr[151];
   ardr[156]=MMH*ardr[85];
   ardr[80]=5./24.*ardr[156] + ardr[151] + ardr[80] + 5./4.*ardr[134]
    + 5*ardr[14] + ardr[148] + 11./2.*ardr[43] + ardr[109] + 3*
   ardr[106] - 5./2.*ardr[44];
   ardr[80]=MMH*ardr[80];
   ardr[106]=ardr[8] - ardr[11];
   ardr[134]=ardr[11]*ardr[106];
   ardr[148]=49./4.*ardr[12] - ardr[8] + 27./4.*ardr[11];
   ardr[148]=ardr[12]*ardr[148];
   ardr[80]=ardr[80] + ardr[134] + ardr[148];
   ardr[80]=ardr[41]*ardr[80];
   ardr[134]=ardr[11]*ardr[85];
   ardr[148]= - ardr[31] + ardr[34];
   ardr[151]=ardr[12]*ardr[148];
   ardr[157]=1./2.*ardr[43] + 1./2.*ardr[44] - ardr[46];
   ardr[134]=ardr[151] + 3./2.*ardr[157] + ardr[134];
   ardr[134]=MMH*ardr[134];
   ardr[151]=1./2.*ardr[12];
   ardr[159]= - ardr[11] + ardr[151];
   ardr[159]=ardr[12]*ardr[159];
   ardr[159]=1./2.*ardr[116] + ardr[159];
   ardr[134]=103./8.*ardr[159] + ardr[134];
   ardr[134]=ardr[41]*ardr[134];
   ardr[164]= - ardr[27] - ardr[22] + ardr[24] + ardr[26];
   ardr[130]=ardr[12]*ardr[130];
   ardr[129]=ardr[130] + 1./2.*ardr[164] + ardr[129];
   ardr[129]=ardr[15]*MMt*ardr[129];
   ardr[129]=ardr[134] + 3./2.*ardr[129];
   ardr[129]=ardr[136]*ardr[129];
   ardr[64]=1./2.*ardr[129] + 1./4.*ardr[80] + ardr[64];
   ardr[64]=ardr[136]*ardr[64];
   ardr[80]=45*ardr[29];
   ardr[129]= - 557./6. + ardr[80];
   ardr[130]=5*ardr[36];
   ardr[134]=5./2.*ardr[32];
   ardr[129]=ardr[134] + 1./2.*ardr[129] + ardr[130];
   ardr[164]= - 7./4.*ardr[34] + 1./4.*ardr[129] - 1./3.*ardr[31];
   ardr[164]=MMH*ardr[164];
   ardr[166]=1./4.*ardr[32] + 1./2.*ardr[36] + 5 + 9./4.*ardr[29];
   ardr[167]= - 1./4.*ardr[34];
   ardr[166]=ardr[167] + 3./2.*ardr[166] + ardr[31];
   ardr[166]=ardr[8]*ardr[166];
   ardr[171]= - 9./4.*ardr[42];
   ardr[177]=7./8.*ardr[44] - 17./2.*ardr[50] - 11./2.*ardr[51] + 
   ardr[171] + 5*ardr[53];
   ardr[188]=5*ardr[39];
   ardr[199]= - 13./8.*ardr[13];
   ardr[201]=ardr[34] + 1./3. + ardr[31];
   ardr[201]=ardr[11]*ardr[201];
   ardr[202]=5./4. - ardr[31];
   ardr[202]=3*ardr[202] + 11./2.*ardr[34];
   ardr[202]=ardr[12]*ardr[202];
   ardr[164]=1./4.*ardr[164] + 1./2.*ardr[202] + 1./4.*ardr[201] + 1./2.
   *ardr[166] + 23./8.*ardr[14] + ardr[199] + ardr[188] + 55./16.*
   ardr[43] + 1./2.*ardr[177] - ardr[46];
   ardr[164]=MMH*ardr[164];
   ardr[166]=15*ardr[8] - 199./4.*ardr[11];
   ardr[166]=ardr[11]*ardr[166];
   ardr[177]= - 59./8.*ardr[12] - 9*ardr[8] + 91./2.*ardr[11];
   ardr[177]=ardr[12]*ardr[177];
   ardr[166]=1./2.*ardr[166] + ardr[177];
   ardr[164]=1./2.*ardr[166] + ardr[164];
   ardr[164]=ardr[41]*ardr[164];
   ardr[63]=ardr[64] + 1./2.*ardr[164] + ardr[63];
   ardr[63]=ardr[136]*ardr[63];
   ardr[64]=121./2. + 27*ardr[29];
   ardr[64]=ardr[95] + 1./2.*ardr[64] + ardr[89];
   ardr[64]=ardr[150] + 1./2.*ardr[64] + ardr[31];
   ardr[64]=1./4.*ardr[8]*ardr[64];
   ardr[129]=ardr[129] - 11./6.*ardr[31];
   ardr[150]= - 3*ardr[34];
   ardr[129]=1./2.*ardr[129] + ardr[150];
   ardr[129]=1./8.*MMH*ardr[129];
   ardr[164]=29./8.*ardr[43] - 7./2.*ardr[46] + 17./8.*ardr[44] - 1./2.
   *ardr[50] + ardr[171] - 5*ardr[51];
   ardr[166]=17./8.*ardr[14];
   ardr[99]=47./12. + ardr[99];
   ardr[99]=1./2.*ardr[99] - ardr[34];
   ardr[99]=1./4.*ardr[11]*ardr[99];
   ardr[171]=27./4.*ardr[34] + 77./12. + ardr[191];
   ardr[171]=1./4.*ardr[12]*ardr[171];
   ardr[164]=ardr[129] + ardr[171] + ardr[99] + ardr[64] + ardr[166] + 
   ardr[199] + 1./2.*ardr[164] + ardr[188];
   ardr[164]=MMH*ardr[164];
   ardr[177]=1./4.*ardr[152];
   ardr[199]= - 25*ardr[11];
   ardr[201]=11./2.*ardr[8] + ardr[199];
   ardr[201]=ardr[11]*ardr[201];
   ardr[202]= - 25*ardr[8] - 71./2.*ardr[11];
   ardr[202]=ardr[12]*ardr[202];
   ardr[201]=1./2.*ardr[202] + ardr[177] + ardr[201];
   ardr[164]=1./4.*ardr[201] + ardr[164];
   ardr[164]=ardr[41]*ardr[164];
   ardr[62]=ardr[63] + ardr[164] + ardr[62];
   ardr[62]=ardr[136]*ardr[62];
   ardr[63]= - 9./8.*ardr[42];
   ardr[164]= - 9./8.*ardr[13];
   ardr[201]= - 1./4.*ardr[50];
   ardr[64]=ardr[129] + ardr[171] + ardr[99] + ardr[64] + ardr[166] + 
   ardr[164] + ardr[188] + 29./16.*ardr[43] - 9./4.*ardr[46] + 17./16.*
   ardr[44] + ardr[201] + ardr[63] - 2*ardr[51];
   ardr[64]=MMH*ardr[64];
   ardr[80]= - 461./6. + ardr[80];
   ardr[99]= - 5*ardr[34];
   ardr[80]=ardr[99] - 7./3.*ardr[31] + ardr[134] + 1./2.*ardr[80] + 
   ardr[130];
   ardr[80]=MMH*ardr[80];
   ardr[129]=5 + 3*ardr[29];
   ardr[129]=ardr[163] + 3./2.*ardr[129] + ardr[36];
   ardr[129]=3./4.*ardr[129] + ardr[34];
   ardr[129]=ardr[8]*ardr[129];
   ardr[118]=ardr[118] - 3./4.*ardr[42] - ardr[53];
   ardr[118]=ardr[129] + 7./4.*ardr[14] - 5./4.*ardr[13] + 9*ardr[39]
    + 19./8.*ardr[43] - 7*ardr[46] + 27./8.*ardr[44] + 3*ardr[118] + 7./
   2.*ardr[50];
   ardr[92]=ardr[92] + 31./48. + ardr[31];
   ardr[92]=ardr[11]*ardr[92];
   ardr[129]=7./3. + 5./8.*ardr[34];
   ardr[129]=ardr[12]*ardr[129];
   ardr[80]=1./16.*ardr[80] + ardr[129] + 1./2.*ardr[118] + ardr[92];
   ardr[80]=MMH*ardr[80];
   ardr[92]= - 3*ardr[8];
   ardr[118]=ardr[92] + 61./8.*ardr[11];
   ardr[118]=ardr[11]*ardr[118];
   ardr[129]=443./4.*ardr[12] - 27*ardr[8] + 43./2.*ardr[11];
   ardr[129]=ardr[12]*ardr[129];
   ardr[118]=1./2.*ardr[129] + ardr[177] + ardr[118];
   ardr[129]= - ardr[8] + ardr[199];
   ardr[75]=1./8.*ardr[129] + ardr[75];
   ardr[75]=ardr[12]*ardr[75];
   ardr[129]=1./12. - ardr[34];
   ardr[129]=ardr[11]*ardr[129];
   ardr[130]=MMH*ardr[34];
   ardr[134]=1./24.*ardr[130];
   ardr[129]=ardr[134] - 1./12.*ardr[12] + 1./2.*ardr[129] + 1./24.*
   ardr[8] + ardr[14] + ardr[50] - ardr[46];
   ardr[129]=MMH*ardr[129];
   ardr[166]=ardr[92] + ardr[11];
   ardr[166]=ardr[11]*ardr[166];
   ardr[75]=1./2.*ardr[129] + 1./8.*ardr[166] + ardr[75];
   ardr[75]=ardr[84]*ardr[75];
   ardr[75]=ardr[75] + 1./2.*ardr[118] + ardr[80];
   ardr[75]=ardr[84]*ardr[75];
   ardr[80]=9./2.*ardr[8] - 23*ardr[11];
   ardr[80]=ardr[11]*ardr[80];
   ardr[118]=ardr[76] - 11./2.*ardr[11];
   ardr[118]=ardr[12]*ardr[118];
   ardr[80]=5./2.*ardr[118] + ardr[177] + ardr[80];
   ardr[64]=1./2.*ardr[75] + 1./4.*ardr[80] + ardr[64];
   ardr[64]=ardr[84]*ardr[64];
   ardr[64]= - 1./4.*ardr[116] + ardr[64];
   ardr[64]=ardr[41]*ardr[64];
   ardr[75]=1./4.*ardr[46];
   ardr[80]=1./6. - ardr[34];
   ardr[80]=ardr[8]*ardr[80];
   ardr[118]=ardr[11]*ardr[34];
   ardr[118]=1./48.*ardr[130] - 1./48.*ardr[12] + 1./8.*ardr[118] + 1./
   8.*ardr[80] + 3./4.*ardr[14] + ardr[75] - ardr[44] + ardr[51] + 
   ardr[201];
   ardr[118]=MMH*ardr[118];
   ardr[117]= - ardr[152] + ardr[117];
   ardr[129]= - 3*ardr[11];
   ardr[76]=ardr[76] + ardr[129];
   ardr[76]=1./8.*ardr[76] + ardr[12];
   ardr[76]=ardr[12]*ardr[76];
   ardr[76]=ardr[118] + 1./8.*ardr[117] + ardr[76];
   ardr[76]=ardr[84]*MMH*ardr[76];
   ardr[117]= - 1./6. + ardr[150];
   ardr[117]=ardr[12]*ardr[117];
   ardr[118]=ardr[11]*ardr[148];
   ardr[80]=ardr[134] + 1./4.*ardr[117] + 1./4.*ardr[118] + 1./4.*
   ardr[80] + 3./2.*ardr[14] - ardr[43] + 3*ardr[46] - 7./2.*ardr[44]
    + 5./2.*ardr[51] - ardr[50];
   ardr[80]=MMH*ardr[80];
   ardr[117]=1./4.*ardr[8] + ardr[11];
   ardr[117]=ardr[11]*ardr[117];
   ardr[130]= - 3./4.*ardr[8] - ardr[11];
   ardr[134]=3*ardr[130] + 17./4.*ardr[12];
   ardr[134]=ardr[12]*ardr[134];
   ardr[80]=ardr[80] + ardr[134] - 1./4.*ardr[152] + ardr[117];
   ardr[80]=MMH*ardr[80];
   ardr[76]=1./2.*ardr[80] + ardr[76];
   ardr[76]=ardr[84]*ardr[76];
   ardr[80]=3*ardr[31];
   ardr[100]=ardr[80] + ardr[100];
   ardr[100]=ardr[12]*ardr[100];
   ardr[117]=1./6. - ardr[31];
   ardr[117]=ardr[8]*ardr[117];
   ardr[134]= - 1./6. - ardr[31];
   ardr[134]=1./2.*ardr[134] + ardr[34];
   ardr[134]=ardr[11]*ardr[134];
   ardr[166]=MMH*ardr[31];
   ardr[100]=1./96.*ardr[166] + 1./16.*ardr[100] + 1./8.*ardr[134] + 1./
   16.*ardr[117] + 7./8.*ardr[13] + ardr[124] + 3./2.*ardr[46] - 5./4.*
   ardr[44] + ardr[51] - 1./8.*ardr[50];
   ardr[100]=MMH*ardr[100];
   ardr[117]=ardr[92] + 5*ardr[11];
   ardr[124]=ardr[11]*ardr[117];
   ardr[124]= - ardr[152] + 3*ardr[124];
   ardr[134]=15./4.*ardr[12] - 11./8.*ardr[8] + ardr[129];
   ardr[134]=ardr[12]*ardr[134];
   ardr[124]=1./8.*ardr[124] + ardr[134];
   ardr[100]=1./2.*ardr[124] + ardr[100];
   ardr[100]=MMH*ardr[100];
   ardr[76]=ardr[100] + 1./2.*ardr[76];
   ardr[76]=ardr[41]*ardr[84]*ardr[76];
   ardr[124]=1./6. + ardr[80];
   ardr[99]=ardr[124] + ardr[99];
   ardr[99]=ardr[12]*ardr[99];
   ardr[134]= - ardr[14] + ardr[13] + ardr[46] - ardr[43];
   ardr[152]=ardr[8]*ardr[148];
   ardr[166]=ardr[34] - 1./6. + ardr[31];
   ardr[166]=ardr[11]*ardr[166];
   ardr[99]=1./12.*ardr[156] + 1./2.*ardr[99] + 1./2.*ardr[166] + 3*
   ardr[134] + 1./2.*ardr[152];
   ardr[99]=MMH*ardr[99];
   ardr[134]=ardr[11]*ardr[198];
   ardr[166]=ardr[8] - 7./2.*ardr[11];
   ardr[166]=1./2.*ardr[166] + ardr[12];
   ardr[166]=ardr[12]*ardr[166];
   ardr[99]=1./2.*ardr[99] + 1./2.*ardr[134] + ardr[166];
   ardr[99]=ardr[41]*MMH*ardr[99];
   ardr[134]= - ardr[8]*ardr[19];
   ardr[166]=ardr[14] - ardr[13] + ardr[26] - ardr[27];
   ardr[171]=ardr[12]*ardr[19];
   ardr[177]=5./4.*ardr[171] + ardr[110] + ardr[166] + 1./4.*ardr[134];
   ardr[188]=MMH*ardr[19];
   ardr[177]=3*ardr[177] + 1./8.*ardr[188];
   ardr[177]=MMt*ardr[177];
   ardr[202]=ardr[11] - ardr[12];
   ardr[203]= - ardr[34] + 1./8. + ardr[31];
   ardr[203]=ardr[16]*ardr[203];
   ardr[204]=MMH*ardr[148];
   ardr[203]=1./8.*ardr[204] + 1./8.*ardr[202] + ardr[203];
   ardr[203]=MMH*ardr[203];
   ardr[117]=ardr[117] + 7*ardr[12];
   ardr[117]=ardr[16]*ardr[117];
   ardr[117]=ardr[177] + 1./4.*ardr[117] + ardr[203];
   ardr[117]=MMt*ardr[117];
   ardr[177]=ardr[16]*ardr[202];
   ardr[148]=MMH*ardr[16]*ardr[148];
   ardr[148]=ardr[177] + ardr[148];
   ardr[148]=MMH*ardr[148];
   ardr[117]=1./4.*ardr[148] + ardr[117];
   ardr[177]=ardr[3]*ardr[179];
   ardr[127]=ardr[127] + ardr[177];
   ardr[127]=MMt*ardr[127];
   ardr[82]=ardr[16]*ardr[82];
   ardr[177]=ardr[3]*ardr[155];
   ardr[82]=ardr[127] + 1./4.*ardr[82] + ardr[177];
   ardr[82]=MMt*ardr[82];
   ardr[82]= - 1./4.*ardr[155] + ardr[82];
   ardr[82]=3*ardr[15]*MMt*ardr[82];
   ardr[117]=1./2.*ardr[117] + ardr[82];
   ardr[117]=ardr[15]*ardr[117];
   ardr[85]=ardr[12]*ardr[85];
   ardr[85]=1./2.*ardr[85] + 1./2.*ardr[118] - 1./2.*ardr[43] - 1./2.*
   ardr[44] + ardr[46];
   ardr[85]=MMH*ardr[85];
   ardr[85]=3./2.*ardr[159] + ardr[85];
   ardr[85]=ardr[41]*MMH*ardr[85];
   ardr[118]=ardr[171] + 1./2.*ardr[166] + ardr[108];
   ardr[118]=MMt*ardr[118];
   ardr[98]=1./2.*ardr[98] + ardr[118];
   ardr[98]=ardr[15]*MMt*ardr[98];
   ardr[85]=ardr[85] + 3*ardr[98];
   ardr[85]=ardr[136]*ardr[85];
   ardr[85]=1./4.*ardr[85] + 1./4.*ardr[99] + ardr[117];
   ardr[85]=ardr[136]*ardr[85];
   ardr[98]= - 3*ardr[44] + ardr[51] + ardr[50];
   ardr[99]= - 1./6. + ardr[34];
   ardr[99]=ardr[11]*ardr[99];
   ardr[117]=1./8.*ardr[124] - ardr[34];
   ardr[117]=ardr[12]*ardr[117];
   ardr[98]=1./96.*ardr[156] + 1./2.*ardr[117] + 1./16.*ardr[99] + 1./
   16.*ardr[152] - 3./8.*ardr[14] + 5./8.*ardr[13] - 7./8.*ardr[43] + 1.
   /8.*ardr[98] + ardr[46];
   ardr[98]=MMH*ardr[98];
   ardr[99]=ardr[92] + 11./2.*ardr[11];
   ardr[99]=ardr[11]*ardr[99];
   ardr[117]=13./16.*ardr[12] - 1./8.*ardr[8] - ardr[11];
   ardr[117]=ardr[12]*ardr[117];
   ardr[98]=ardr[98] + 1./8.*ardr[99] + ardr[117];
   ardr[98]=ardr[41]*MMH*ardr[98];
   ardr[99]=ardr[130] + 31./4.*ardr[12];
   ardr[99]=ardr[16]*ardr[99];
   ardr[99]=ardr[99] + ardr[203];
   ardr[117]=ardr[103] + 1./2.*ardr[134] - 1./2.*ardr[14] + ardr[68] + 
   1./2.*ardr[27] - ardr[28] + 1./2.*ardr[26];
   ardr[117]=1./4.*ardr[117] + ardr[171];
   ardr[117]=3*ardr[117] + 1./16.*ardr[188];
   ardr[117]=MMt*ardr[117];
   ardr[99]=1./2.*ardr[99] + ardr[117];
   ardr[99]=MMt*ardr[99];
   ardr[82]=ardr[82] + 1./8.*ardr[148] + ardr[99];
   ardr[82]=ardr[15]*ardr[82];
   ardr[82]=ardr[85] + ardr[98] + ardr[82];
   ardr[82]=ardr[136]*ardr[82];
   ardr[85]=ardr[41]*ardr[100];
   ardr[98]= - ardr[8] + ardr[11];
   ardr[99]=ardr[16]*ardr[31];
   ardr[100]= - MMH*ardr[31];
   ardr[99]=1./8.*ardr[100] + 1./8.*ardr[98] + ardr[99];
   ardr[99]=MMH*ardr[99];
   ardr[100]=ardr[106] + 75./4.*ardr[12];
   ardr[100]=ardr[16]*ardr[100];
   ardr[99]=ardr[100] + ardr[99];
   ardr[100]=11./8.*ardr[171] + ardr[110] - 7./4.*ardr[14] - 1./4.*
   ardr[13] + 7./4.*ardr[27] + ardr[200] + 1./4.*ardr[26];
   ardr[100]=MMt*ardr[100];
   ardr[99]=1./2.*ardr[99] + 3*ardr[100];
   ardr[99]=MMt*ardr[99];
   ardr[98]=ardr[16]*ardr[98];
   ardr[100]= - MMH*ardr[16]*ardr[31];
   ardr[98]=ardr[98] + ardr[100];
   ardr[98]=MMH*ardr[98];
   ardr[98]=1./8.*ardr[98] + ardr[99];
   ardr[99]=ardr[15]*ardr[98];
   ardr[82]=ardr[82] + ardr[85] + ardr[99];
   ardr[82]=ardr[136]*ardr[82];
   ardr[61]=3./2.*ardr[171] + ardr[103] + 1./2.*ardr[70] - 13./2.*
   ardr[14] + 3./2.*ardr[13] + 13./2.*ardr[27] + ardr[197] + ardr[61];
   ardr[85]= - MMH*ardr[19];
   ardr[61]=3*ardr[61] + 1./4.*ardr[85];
   ardr[61]=MMt*ardr[61];
   ardr[99]= - ardr[8] + ardr[12];
   ardr[100]= - 1./8. + ardr[34];
   ardr[100]=ardr[16]*ardr[100];
   ardr[103]= - MMH*ardr[34];
   ardr[100]=1./8.*ardr[103] + 1./8.*ardr[99] + ardr[100];
   ardr[100]=MMH*ardr[100];
   ardr[103]=11*ardr[12] + 7./4.*ardr[8] + ardr[129];
   ardr[103]=ardr[16]*ardr[103];
   ardr[61]=1./2.*ardr[61] + ardr[103] + ardr[100];
   ardr[61]=MMt*ardr[61];
   ardr[99]=ardr[16]*ardr[99];
   ardr[103]= - MMH*ardr[16]*ardr[34];
   ardr[99]=ardr[99] + ardr[103];
   ardr[99]=MMH*ardr[99];
   ardr[61]=1./4.*ardr[99] + ardr[61];
   ardr[103]=7*ardr[8] + ardr[129];
   ardr[103]=1./4.*ardr[103] + 5*ardr[12];
   ardr[103]=ardr[16]*ardr[103];
   ardr[100]=ardr[103] + ardr[100];
   ardr[70]=1./8.*ardr[108] + 1./8.*ardr[70] + ardr[170] - ardr[14];
   ardr[70]=3*ardr[70] + 1./16.*ardr[85];
   ardr[70]=MMt*ardr[70];
   ardr[70]=1./2.*ardr[100] + ardr[70];
   ardr[70]=MMt*ardr[70];
   ardr[70]=1./8.*ardr[99] + ardr[70];
   ardr[70]=ardr[84]*ardr[70];
   ardr[61]=1./2.*ardr[61] + ardr[70];
   ardr[61]=ardr[84]*ardr[61];
   ardr[61]=ardr[98] + ardr[61];
   ardr[61]=ardr[84]*ardr[61];
   ardr[70]= - ardr[3]*ardr[16]*ardr[19];
   ardr[70]=ardr[169] + ardr[70];
   ardr[70]=MMt*ardr[70];
   ardr[85]=1./2. + ardr[19];
   ardr[85]=ardr[16]*ardr[85];
   ardr[70]=ardr[70] + 1./4.*ardr[85] + ardr[168];
   ardr[70]=MMt*ardr[70];
   ardr[70]=ardr[193] + ardr[70];
   ardr[70]=MMt*ardr[70];
   ardr[85]=ardr[84]*ardr[70];
   ardr[70]=ardr[70] + ardr[85];
   ardr[70]=ardr[15]*ardr[70]*pow(ardr[2],4);
   ardr[61]=ardr[61] + 3*ardr[70];
   ardr[61]=ardr[15]*ardr[61];
   ardr[70]=ardr[11] - 1./2.*ardr[12];
   ardr[70]=ardr[12]*ardr[70];
   ardr[85]=MMH*ardr[157];
   ardr[70]=ardr[85] - 1./2.*ardr[116] + ardr[70];
   ardr[70]=ardr[70]*pow(MMH,2);
   ardr[85]=ardr[84]*ardr[70];
   ardr[85]=ardr[70] + 1./4.*ardr[85];
   ardr[85]=ardr[41]*ardr[84]*ardr[85];
   ardr[70]=ardr[41]*ardr[70];
   ardr[98]=ardr[136]*ardr[70];
   ardr[98]=ardr[70] + 1./2.*ardr[98];
   ardr[98]=ardr[136]*ardr[98];
   ardr[98]=3./2.*ardr[70] + ardr[98];
   ardr[98]=ardr[136]*ardr[98];
   ardr[70]=ardr[70] + 1./2.*ardr[98];
   ardr[70]=ardr[136]*ardr[70];
   ardr[70]=ardr[85] + ardr[70];
   ardr[70]=ardr[4]*ardr[70];
   ardr[61]=1./4.*ardr[70] + ardr[82] + ardr[76] + ardr[61];
   ardr[61]=ardr[4]*ardr[61];
   ardr[61]=ardr[61] + ardr[62] + ardr[64] + ardr[69];
   ardr[61]=ardr[4]*ardr[61];
   ardr[62]= - 47./3. + ardr[145];
   ardr[62]=63./16.*ardr[181] + 63./16.*ardr[180] + ardr[113] + 
   ardr[178] + ardr[162] + ardr[141] - 17./4.*ardr[37] + 1./2.*ardr[62]
    + ardr[143];
   ardr[62]=ardr[11]*ardr[62];
   ardr[64]=ardr[187] + ardr[186] + ardr[185] + ardr[184] + 1./48.*
   ardr[37] + ardr[183] - 1./3. + ardr[182];
   ardr[64]=MMH*ardr[64];
   ardr[69]=ardr[142] + ardr[191] + ardr[190] + ardr[189] - 1./4.*
   ardr[37] + 43./6. + ardr[143];
   ardr[69]=ardr[8]*ardr[69];
   ardr[70]= - 135./16.*ardr[29];
   ardr[76]=9./2.*ardr[31];
   ardr[82]=ardr[138] + ardr[119] + ardr[91] + ardr[76] - 59./8.*
   ardr[35] - 15./16.*ardr[32] + 3./8.*ardr[37] - 15./8.*ardr[36] + 1./
   3. + ardr[70];
   ardr[82]=ardr[12]*ardr[82];
   ardr[85]= - 281*ardr[8] + 69*ardr[11];
   ardr[85]=ardr[11]*ardr[85];
   ardr[91]=ardr[198] - ardr[12];
   ardr[91]=ardr[12]*ardr[91];
   ardr[85]=25*ardr[91] + ardr[192] + 1./8.*ardr[85];
   ardr[85]=ardr[56]*ardr[85];
   ardr[62]=1./4.*ardr[85] + ardr[64] + ardr[82] + 1./2.*ardr[62] + 1./
   2.*ardr[69] - 69./2.*ardr[14] + 133./8.*ardr[13] - 69./8.*ardr[39]
    - 45./8.*ardr[43] - 3./2.*ardr[6] + 7*ardr[46] + 35./2.*ardr[7] - 
   23./8.*ardr[44] + 11./4.*ardr[50] + 159./8.*ardr[48] + 23./4.*
   ardr[51] - 13./4.*ardr[47] - 1./2.*ardr[53] + ardr[146] - 27*
   ardr[52];
   ardr[64]=ardr[12]*ardr[11];
   ardr[64]= - ardr[116] + ardr[64];
   ardr[64]=ardr[3]*ardr[64];
   ardr[69]= - 1./4.*ardr[35];
   ardr[82]= - 1./3. + ardr[69];
   ardr[82]=ardr[8]*ardr[82];
   ardr[85]=1./6. + ardr[35];
   ardr[85]=ardr[11]*ardr[85];
   ardr[91]=1./8.*ardr[35] - ardr[34];
   ardr[91]=MMH*ardr[91];
   ardr[64]=1./4.*ardr[64] + 1./6.*ardr[91] + 7./48.*ardr[12] + 1./8.*
   ardr[85] + 1./2.*ardr[82] + 13./4.*ardr[14] + ardr[123] + ardr[75]
    + ardr[201] + 5./2.*ardr[52] - 3*ardr[48];
   ardr[64]=ardr[84]*ardr[64];
   ardr[82]= - 9*ardr[29];
   ardr[85]=5./2. + ardr[82];
   ardr[85]=1./2.*ardr[85] - ardr[36];
   ardr[77]=ardr[77] + ardr[141] + 3*ardr[85] + ardr[87];
   ardr[77]=1./2.*ardr[77] - 7*ardr[31];
   ardr[85]=15./8.*ardr[181];
   ardr[77]=ardr[85] + 15./8.*ardr[180] + 1./2.*ardr[77] + ardr[34];
   ardr[77]=ardr[11]*ardr[77];
   ardr[87]= - 17./2.*ardr[35];
   ardr[91]=ardr[87] - 11./6. + ardr[189];
   ardr[91]=ardr[8]*ardr[91];
   ardr[98]= - 23*ardr[8] + 9./2.*ardr[11];
   ardr[98]=ardr[11]*ardr[98];
   ardr[99]=3./2.*ardr[8] + 19*ardr[11];
   ardr[99]=1./2.*ardr[99] - ardr[12];
   ardr[99]=ardr[12]*ardr[99];
   ardr[98]=1./2.*ardr[98] + ardr[99];
   ardr[98]=ardr[56]*ardr[98];
   ardr[99]=ardr[154] + ardr[199];
   ardr[99]=ardr[11]*ardr[99];
   ardr[100]= - ardr[12]*ardr[11];
   ardr[99]=1./2.*ardr[99] + 33*ardr[100];
   ardr[99]=ardr[3]*ardr[99];
   ardr[100]= - 3*ardr[35];
   ardr[103]=125./6. + ardr[100];
   ardr[103]=ardr[12]*ardr[103];
   ardr[106]=17./6.*ardr[35] + 31./6. - 5*ardr[32];
   ardr[106]= - 1./2.*ardr[34] + 1./4.*ardr[106] + 7./3.*ardr[31];
   ardr[106]=MMH*ardr[106];
   ardr[64]=1./2.*ardr[64] + 1./4.*ardr[99] + 1./4.*ardr[98] + 1./4.*
   ardr[106] + 1./16.*ardr[103] + 1./2.*ardr[77] + 1./8.*ardr[91] - 11./
   4.*ardr[14] - 31./16.*ardr[13] - 17./16.*ardr[39] + 3./16.*ardr[43]
    - 9./8.*ardr[6] + 9./4.*ardr[46] - 11./8.*ardr[50] - 5./16.*
   ardr[48] + 1./8.*ardr[52] + ardr[47];
   ardr[64]=ardr[84]*ardr[64];
   ardr[77]=ardr[194] - 10*ardr[12];
   ardr[77]=ardr[12]*ardr[77];
   ardr[77]=ardr[195] + ardr[77];
   ardr[77]=ardr[3]*ardr[77];
   ardr[62]=ardr[64] + 1./2.*ardr[62] + ardr[77];
   ardr[62]=ardr[84]*ardr[62];
   ardr[64]=pow(CW,2);
   ardr[77]=9./2. - ardr[64];
   ardr[77]=ardr[52]*ardr[77];
   ardr[91]= - 1 + 1./4.*ardr[64];
   ardr[91]=ardr[47]*ardr[91];
   ardr[98]= - 5 + ardr[64];
   ardr[98]=ardr[7]*ardr[98];
   ardr[62]=ardr[62] + ardr[68] + ardr[75] + 1./4.*ardr[98] - 1./4.*
   ardr[51] + 1./2.*ardr[77] + ardr[91];
   ardr[62]=ardr[41]*ardr[62];
   ardr[68]= - 1 - ardr[9];
   ardr[75]=ardr[20]*ardr[68];
   ardr[68]=ardr[1]*ardr[68];
   ardr[68]=11./3.*ardr[75] + 3*ardr[68];
   ardr[68]=ardr[11]*ardr[68];
   ardr[75]=ardr[6] - ardr[13];
   ardr[77]=ardr[20]*ardr[75];
   ardr[75]=ardr[1]*ardr[75];
   ardr[68]=ardr[68] + 11./3.*ardr[77] + 3*ardr[75];
   ardr[68]=ardr[84]*ardr[68];
   ardr[68]=ardr[107] + 1./4.*ardr[68];
   ardr[68]=ardr[84]*ardr[68];
   ardr[59]=ardr[61] + ardr[59] + ardr[60] + ardr[68] + ardr[62];
   ardr[59]=ardr[4]*ardr[59];
   ardr[60]=9*ardr[31];
   ardr[61]=3*ardr[32];
   ardr[62]=ardr[150] + ardr[60] + ardr[190] + ardr[61] + 179./4.*
   ardr[37] + 125./2. + ardr[89];
   ardr[62]=ardr[11]*ardr[62];
   ardr[68]=53*ardr[37];
   ardr[75]=19./4.*ardr[35];
   ardr[77]=ardr[75] + ardr[95] + 2201./24. + ardr[68];
   ardr[77]=ardr[12]*ardr[77];
   ardr[91]= - 1 + 2*ardr[36];
   ardr[91]=2*ardr[91] - ardr[32];
   ardr[98]=pow(MMZ,2);
   ardr[91]=ardr[98]*ardr[91];
   ardr[99]=ardr[141] - 1 + ardr[89];
   ardr[99]=ardr[12]*MMZ*ardr[99];
   ardr[103]=1 + ardr[141];
   ardr[106]=MMZ*ardr[103];
   ardr[107]=ardr[11]*ardr[106];
   ardr[91]=ardr[99] + ardr[91] + ardr[107];
   ardr[91]=ardr[3]*ardr[91];
   ardr[99]=3*ardr[39];
   ardr[108]=17*ardr[13] + ardr[99] - 11*ardr[48] - 3*ardr[43];
   ardr[108]=1./2.*ardr[108] + 11*ardr[14];
   ardr[110]= - 4*ardr[34] + ardr[144] + 151./4.*ardr[35] + ardr[32] + 
   205./2.*ardr[37] - 835./18. - ardr[36];
   ardr[110]=MMZ*ardr[110];
   ardr[113]=9*ardr[32];
   ardr[116]=11 + ardr[113];
   ardr[116]=ardr[8]*ardr[116];
   ardr[117]=1./2.*ardr[116];
   ardr[62]=3*ardr[91] + ardr[77] + 1./2.*ardr[62] + ardr[117] + 3*
   ardr[108] + ardr[110];
   ardr[62]=ardr[3]*ardr[62];
   ardr[77]= - 1./4.*ardr[31];
   ardr[68]=ardr[77] + ardr[75] + ardr[95] + 441./4. + ardr[68];
   ardr[68]=ardr[12]*ardr[68];
   ardr[91]=3./8.*ardr[36] + 1 + 3./16.*ardr[29];
   ardr[108]= - 11./16.*ardr[32];
   ardr[110]=55./8.*ardr[35];
   ardr[118]= - 3./8.*ardr[34];
   ardr[91]=ardr[118] + ardr[110] + 3*ardr[91] + ardr[108];
   ardr[91]=MMZ*ardr[91];
   ardr[91]=ardr[91] + ardr[92];
   ardr[91]=ardr[12]*ardr[91];
   ardr[92]=11./2.*ardr[47] - ardr[45] - 23./2.*ardr[49];
   ardr[92]=21./2.*ardr[14] - 9./2.*ardr[39] - 11./16.*ardr[43] + 55./8.
   *ardr[46] + 25./16.*ardr[44] - 55./8.*ardr[48] + 1./8.*ardr[92] - 3*
   ardr[51];
   ardr[92]=1./2.*ardr[92] - 3*MMZ;
   ardr[92]=MMZ*ardr[92];
   ardr[119]= - 3./8.*ardr[36] + 1 - 3./16.*ardr[29];
   ardr[119]=3./8.*ardr[34] - 55./8.*ardr[35] + 3*ardr[119] + 11./16.*
   ardr[32];
   ardr[119]=ardr[8]*MMZ*ardr[119];
   ardr[91]=1./2.*ardr[91] + ardr[92] + 1./2.*ardr[119];
   ardr[91]=ardr[56]*ardr[91];
   ardr[92]=ardr[118] + ardr[110] + ardr[108] + 9./8.*ardr[36] - 17 + 9.
   /16.*ardr[29];
   ardr[92]=MMZ*ardr[92];
   ardr[108]=ardr[7] + 3./2.*ardr[50] - 91./8.*ardr[48] - 15*ardr[51]
    - 1./2.*ardr[47] - 3./2.*ardr[52] - ardr[53];
   ardr[110]=1./16.*ardr[31] + 47./16.*ardr[35] - 9./4.*ardr[32] - 12
    - 53./4.*ardr[37];
   ardr[110]=ardr[8]*ardr[110];
   ardr[68]=3*ardr[91] + 1./4.*ardr[68] - 15./16.*ardr[11] + ardr[110]
    + 3./2.*ardr[92] + 115./8.*ardr[14] - 15./16.*ardr[13] + 19./16.*
   ardr[39] + 17./16.*ardr[43] + 1./2.*ardr[108] + 5*ardr[46];
   ardr[68]=ardr[56]*ardr[68];
   ardr[91]=1853./48. + ardr[82];
   ardr[92]= - ardr[43]*ardr[58];
   ardr[108]=ardr[39]*ardr[58];
   ardr[110]=ardr[13]*ardr[58];
   ardr[118]=2*ardr[32];
   ardr[119]= - 23./8.*ardr[31];
   ardr[124]=11./24.*ardr[34];
   ardr[127]=ardr[55] - 2*ardr[58];
   ardr[127]=MMZ*ardr[127];
   ardr[62]=ardr[62] + ardr[68] + 2*ardr[181] + ardr[94] + 3*ardr[127]
    + ardr[124] + ardr[119] - 39./8.*ardr[35] + ardr[118] + 2*ardr[110]
    - 405./32.*ardr[37] + ardr[108] - 1./2.*ardr[36] + 1./4.*ardr[91]
    + ardr[92];
   ardr[62]=ardr[41]*ardr[62];
   ardr[68]= - ardr[22] - ardr[6];
   ardr[91]=2*ardr[14];
   ardr[68]=ardr[91] + ardr[13] + 3*ardr[40] + 1./2.*ardr[68] - 2*
   ardr[27];
   ardr[127]=7./2.*ardr[9];
   ardr[96]=ardr[96] + ardr[127] + 3 + 11./2.*ardr[18];
   ardr[96]=ardr[11]*ardr[96];
   ardr[129]=4*ardr[18];
   ardr[130]=2*ardr[9];
   ardr[134]=ardr[130] + 3 + ardr[129];
   ardr[134]=ardr[12]*ardr[134];
   ardr[138]=10 + ardr[160];
   ardr[141]=3*ardr[33];
   ardr[138]=2*ardr[138] + ardr[141];
   ardr[138]=ardr[16]*ardr[138];
   ardr[142]=17./4.*ardr[9] - 11 + 25./4.*ardr[18];
   ardr[142]=MMZ*ardr[142];
   ardr[68]=ardr[138] + ardr[134] + 1./2.*ardr[96] + 3*ardr[68] + 
   ardr[142];
   ardr[68]=ardr[3]*ardr[68];
   ardr[96]=ardr[14] - ardr[39] + ardr[25] - ardr[27];
   ardr[96]=MMZ*ardr[96];
   ardr[134]= - MMZ*ardr[19];
   ardr[138]=ardr[8]*ardr[134];
   ardr[142]=MMZ*ardr[19];
   ardr[145]=ardr[12]*ardr[142];
   ardr[96]=ardr[145] + ardr[96] + ardr[138];
   ardr[96]=ardr[56]*ardr[96];
   ardr[138]= - 3*ardr[27] + ardr[21] + ardr[25] - ardr[24];
   ardr[138]= - ardr[39] + 1./2.*ardr[138] + ardr[40];
   ardr[145]= - 1 + ardr[19];
   ardr[145]=MMZ*ardr[145];
   ardr[138]=1./2.*ardr[145] + 1./2.*ardr[138] + ardr[14];
   ardr[93]=ardr[93] - ardr[18] - 1./2.*ardr[9];
   ardr[93]=ardr[8]*ardr[93];
   ardr[115]=ardr[115] + 3./4. + ardr[18];
   ardr[115]=ardr[12]*ardr[115];
   ardr[93]=3./4.*ardr[96] + 3./4.*ardr[16] + ardr[115] + 3./2.*
   ardr[138] + ardr[93];
   ardr[93]=ardr[56]*ardr[93];
   ardr[96]=ardr[27] - 3*ardr[21] - ardr[25] + 3*ardr[24];
   ardr[83]=ardr[83] + ardr[19];
   ardr[83]=ardr[8]*ardr[83];
   ardr[83]=ardr[105] + 1./4.*ardr[83] - ardr[14] + 1./4.*ardr[96] + 
   ardr[39];
   ardr[83]=ardr[56]*ardr[83];
   ardr[83]=1./4.*ardr[104] + ardr[83];
   ardr[83]=ardr[56]*ardr[83];
   ardr[96]= - 1 + ardr[36];
   ardr[96]=2*ardr[96] + ardr[32];
   ardr[96]=ardr[3]*ardr[16]*ardr[96];
   ardr[96]=9*ardr[96] + ardr[73] + 23./4.*ardr[18] - 9./2.*ardr[17] + 
   3./2.*ardr[33] + 11 + ardr[160];
   ardr[96]=ardr[3]*ardr[96];
   ardr[104]=MMt*ardr[161]*ardr[17];
   ardr[83]=18*ardr[104] + 3./2.*ardr[83] + ardr[96];
   ardr[83]=MMt*ardr[83];
   ardr[96]= - ardr[54]*ardr[22];
   ardr[105]=ardr[40]*ardr[54];
   ardr[115]=15./2.*ardr[105] + 15./4.*ardr[96] - 9./2.*ardr[33] + 149./
   8. - 9*ardr[38];
   ardr[138]= - 1./2.*ardr[6];
   ardr[145]= - ardr[27] + ardr[138] + ardr[79] + 1./2.*ardr[7] + 1./2.
   *ardr[24] + ardr[26];
   ardr[145]=ardr[23]*ardr[145];
   ardr[145]=3./2. + ardr[145];
   ardr[145]=EPAIR2*ardr[145];
   ardr[79]=ardr[138] + ardr[26] + ardr[79];
   ardr[138]=pow(ardr[23],2);
   ardr[79]=ardr[138]*ardr[79];
   ardr[79]= - 15./16.*ardr[54] + ardr[79];
   ardr[79]=MMZ*ardr[79];
   ardr[146]= - ardr[38] - 1./2.*ardr[33];
   ardr[146]=ardr[23]*ardr[146];
   ardr[148]=ardr[16]*ardr[138];
   ardr[146]=1./2.*ardr[148] + 5./4.*ardr[54] + ardr[146];
   ardr[146]=ardr[16]*ardr[146];
   ardr[112]=15./8.*ardr[40] + ardr[112] + 5./8.*ardr[6] + ardr[22] + 5.
   /4.*ardr[7] - ardr[24] - 13./8.*ardr[26];
   ardr[112]=ardr[23]*ardr[112];
   ardr[148]=ardr[13]*ardr[54];
   ardr[152]=ardr[11]*ardr[54];
   ardr[68]=ardr[83] + ardr[68] + ardr[93] + 3./4.*ardr[146] + 15./32.*
   ardr[152] + 3./4.*ardr[79] + ardr[71] + ardr[165] - 3./2.*ardr[18]
    + 3./2.*ardr[145] + 15./32.*ardr[148] + 1./8.*ardr[115] + ardr[112]
   ;
   ardr[68]=ardr[15]*ardr[68];
   ardr[71]=ardr[76] + 33./2.*ardr[35] + ardr[95] - 825./8.*ardr[37] + 
   125./6. + ardr[89];
   ardr[71]=MMZ*ardr[71];
   ardr[76]= - 1./2.*ardr[32];
   ardr[79]=ardr[76] + 1 - ardr[36];
   ardr[83]=ardr[98]*ardr[79];
   ardr[79]=MMZ*ardr[79];
   ardr[93]=ardr[11]*ardr[79];
   ardr[79]=ardr[12]*ardr[79];
   ardr[79]=ardr[79] + ardr[83] + 1./2.*ardr[93];
   ardr[79]=ardr[3]*ardr[79];
   ardr[83]=ardr[95] + 11./2. + ardr[89];
   ardr[83]=ardr[8]*ardr[83];
   ardr[74]= - ardr[34] + ardr[31] + ardr[74] + ardr[163] - 121./4.*
   ardr[37] - 97./8. + ardr[36];
   ardr[74]=ardr[11]*ardr[74];
   ardr[93]= - 2*ardr[7];
   ardr[95]= - 9./4.*ardr[43];
   ardr[112]=11./2.*ardr[35] + ardr[163] - 22*ardr[37] - 175./8. + 
   ardr[36];
   ardr[112]=ardr[12]*ardr[112];
   ardr[71]=9*ardr[79] + 3*ardr[112] + 3./2.*ardr[74] + 3./2.*ardr[83]
    + ardr[71] - 81./2.*ardr[14] - 81./4.*ardr[13] + 15./4.*ardr[39] + 
   ardr[95] - ardr[6] + ardr[93] - 9./2.*ardr[44] + ardr[50] + 2*
   ardr[51] + 99./4.*ardr[48];
   ardr[71]=ardr[3]*ardr[71];
   ardr[74]=ardr[89] - 59 - 27./2.*ardr[29];
   ardr[74]=ardr[167] + ardr[77] + 33./4.*ardr[35] + 3./8.*ardr[32] + 1.
   /4.*ardr[74] + ardr[121];
   ardr[74]=ardr[12]*ardr[74];
   ardr[63]= - 47./2.*ardr[14] - 21./4.*ardr[13] + 43./4.*ardr[39] - 25.
   /16.*ardr[43] - 33./8.*ardr[46] + ardr[122] - 67./16.*ardr[44] + 33./
   4.*ardr[48] + 17./8.*ardr[47] + ardr[147] + 45./8.*ardr[49] + 
   ardr[45] + ardr[63];
   ardr[77]= - 9./2.*ardr[43] - 33*ardr[46] - 25./2.*ardr[44] + 33*
   ardr[48] + 9./2.*ardr[47] + ardr[45] + 23./2.*ardr[49];
   ardr[77]=ardr[135] + 1./16.*ardr[77] + ardr[99];
   ardr[77]=MMZ*ardr[77];
   ardr[76]=ardr[76] - 1./2.*ardr[29] - ardr[36];
   ardr[76]=ardr[34] + 3*ardr[76] + ardr[102];
   ardr[76]=MMZ*ardr[76];
   ardr[76]=ardr[151] + 3./4.*ardr[76] - ardr[8];
   ardr[76]=ardr[12]*ardr[76];
   ardr[79]=ardr[163] + 1./2.*ardr[29] + ardr[36];
   ardr[79]= - ardr[34] + 3*ardr[79] + 11*ardr[35];
   ardr[79]=MMZ*ardr[79];
   ardr[79]=3./2.*ardr[79] + ardr[8];
   ardr[79]=ardr[8]*ardr[79];
   ardr[76]=1./4.*ardr[76] + ardr[77] + 1./8.*ardr[79];
   ardr[76]=ardr[56]*ardr[76];
   ardr[77]=1./16.*ardr[34] - 11./16.*ardr[35] - 3./32.*ardr[32] - 3./
   16.*ardr[36] + 1 - 3./32.*ardr[29];
   ardr[77]=MMZ*ardr[77];
   ardr[79]=11./16.*ardr[34] + 1./8.*ardr[31] - 165./16.*ardr[35] - 33./
   32.*ardr[32] + 33./2.*ardr[37] - 33./16.*ardr[36] + 1 + 27./32.*
   ardr[29];
   ardr[79]=ardr[8]*ardr[79];
   ardr[63]=3*ardr[76] + 1./2.*ardr[74] - 21./8.*ardr[11] + ardr[79] + 
   1./2.*ardr[63] + 9*ardr[77];
   ardr[63]=ardr[56]*ardr[63];
   ardr[74]=ardr[39]*ardr[55];
   ardr[76]= - ardr[44]*ardr[55];
   ardr[70]= - 99./16.*ardr[35] + 27./16.*ardr[32] + 99./2.*ardr[37] + 
   ardr[74] + 27./8.*ardr[36] + ardr[70] - 13511./192. + ardr[76];
   ardr[74]=MMZ*ardr[55];
   ardr[77]=ardr[8]*ardr[55];
   ardr[79]=ardr[14]*ardr[55];
   ardr[63]=ardr[71] + ardr[63] + ardr[111] + 1./2.*ardr[77] + 3./2.*
   ardr[74] - 25./16.*ardr[34] + ardr[176] + 1./2.*ardr[70] + ardr[79];
   ardr[63]=ardr[41]*ardr[63];
   ardr[70]=ardr[54]*ardr[22];
   ardr[71]=3./2. + ardr[70];
   ardr[74]= - ardr[40]*ardr[54];
   ardr[71]=1./2.*ardr[71] + ardr[74];
   ardr[77]=ardr[8]*ardr[131];
   ardr[83]=ardr[12]*ardr[133];
   ardr[77]=ardr[77] + ardr[83];
   ardr[77]=ardr[56]*ardr[77];
   ardr[99]=ardr[6] - ardr[22] + ardr[24] - ardr[7];
   ardr[99]=ardr[23]*ardr[99];
   ardr[102]= - ardr[13]*ardr[54];
   ardr[112]=MMZ*ardr[54];
   ardr[115]= - ardr[11]*ardr[54];
   ardr[71]=ardr[77] + 3./8.*ardr[126] + 3./16.*ardr[115] + 9./32.*
   ardr[112] + ardr[73] + 3./4.*ardr[9] + 3./4.*ardr[18] + 3./16.*
   ardr[102] + 3./8.*ardr[71] + ardr[99];
   ardr[77]=MMZ*ardr[133];
   ardr[77]=ardr[83] + ardr[77] + 1./2.*ardr[137];
   ardr[77]=ardr[3]*ardr[77];
   ardr[83]=MMt*ardr[3]*ardr[125];
   ardr[71]=1./2.*ardr[83] + 1./4.*ardr[71] + ardr[77];
   ardr[71]=ardr[15]*ardr[71];
   ardr[77]=ardr[43]*ardr[58];
   ardr[83]=ardr[55] - ardr[58];
   ardr[99]=ardr[39]*ardr[83];
   ardr[112]= - ardr[13]*ardr[58];
   ardr[76]=1./2.*ardr[79] + 99./32.*ardr[35] + 1./2.*ardr[112] - 99./
   32.*ardr[37] + 1./4.*ardr[99] + 1./4.*ardr[77] + 3 + 1./4.*ardr[76];
   ardr[77]= - ardr[43] + ardr[44] - ardr[49] + ardr[47];
   ardr[79]= - ardr[34] + 11*ardr[120] + ardr[31];
   ardr[79]=ardr[8]*ardr[79];
   ardr[99]=ardr[34] + 11*ardr[114] - ardr[31];
   ardr[99]=ardr[12]*ardr[99];
   ardr[77]=3*ardr[99] + 11*ardr[77] + 3*ardr[79];
   ardr[77]=ardr[56]*ardr[77];
   ardr[72]=ardr[72] + ardr[176] - 99./8.*ardr[35] - 2 + 99./8.*
   ardr[37];
   ardr[72]=ardr[11]*ardr[72];
   ardr[79]=MMZ*ardr[139];
   ardr[99]=3*ardr[34] + ardr[191] - 99./4.*ardr[35] + 2 + 99./4.*
   ardr[37];
   ardr[99]=ardr[12]*ardr[99];
   ardr[72]=ardr[99] + 3*ardr[79] + ardr[72];
   ardr[72]=ardr[3]*ardr[72];
   ardr[79]= - ardr[55] + ardr[58];
   ardr[79]=MMZ*ardr[79];
   ardr[83]=ardr[8]*ardr[83];
   ardr[72]=ardr[72] + 3./16.*ardr[77] + 3./2.*ardr[111] + 3./2.*
   ardr[97] + 3./4.*ardr[83] + 9./2.*ardr[79] - 17./8.*ardr[34] + 3*
   ardr[76] + 17./8.*ardr[31];
   ardr[72]=ardr[41]*ardr[72];
   ardr[76]=MMZ*ardr[132];
   ardr[77]=ardr[12]*ardr[132];
   ardr[76]=ardr[77] + ardr[76] + 1./2.*ardr[140];
   ardr[76]=ardr[3]*ardr[76];
   ardr[77]=ardr[88] + ardr[77];
   ardr[77]=ardr[56]*ardr[77];
   ardr[77]=3./2.*ardr[90] + ardr[77];
   ardr[71]=3*ardr[71] + ardr[72] + 1./4.*ardr[77] + ardr[76];
   ardr[71]=ardr[136]*ardr[71];
   ardr[72]=ardr[14] - ardr[39] + ardr[53] - ardr[7];
   ardr[76]=ardr[20]*ardr[72];
   ardr[72]=ardr[1]*ardr[72];
   ardr[72]=3*ardr[76] + ardr[72];
   ardr[72]=MMZ*ardr[72];
   ardr[76]= - ardr[20]*ardr[10];
   ardr[77]= - ardr[1]*ardr[10];
   ardr[76]=3*ardr[76] + ardr[77];
   ardr[76]=MMZ*ardr[76];
   ardr[77]=ardr[8]*ardr[76];
   ardr[79]=ardr[20]*ardr[10];
   ardr[83]=ardr[1]*ardr[10];
   ardr[79]=3*ardr[79] + ardr[83];
   ardr[79]=MMZ*ardr[79];
   ardr[83]=ardr[12]*ardr[79];
   ardr[72]=ardr[83] + ardr[72] + ardr[77];
   ardr[72]=ardr[56]*ardr[72];
   ardr[77]= - 1./2.*ardr[39];
   ardr[83]=ardr[14] + ardr[77] + ardr[147] - ardr[7];
   ardr[88]=ardr[20]*ardr[83];
   ardr[83]=ardr[1]*ardr[83];
   ardr[90]= - 1 + ardr[10];
   ardr[97]=ardr[20]*ardr[90];
   ardr[90]=ardr[1]*ardr[90];
   ardr[90]=3*ardr[97] + ardr[90];
   ardr[90]=MMZ*ardr[90];
   ardr[97]=ardr[153] - ardr[9];
   ardr[99]=ardr[20]*ardr[97];
   ardr[97]=ardr[1]*ardr[97];
   ardr[97]=3*ardr[99] + ardr[97];
   ardr[97]=ardr[8]*ardr[97];
   ardr[99]=1./2. + ardr[9];
   ardr[111]=ardr[20]*ardr[99];
   ardr[99]=ardr[1]*ardr[99];
   ardr[112]=3*ardr[111] + ardr[99];
   ardr[112]=ardr[12]*ardr[112];
   ardr[72]=1./2.*ardr[72] + ardr[112] + ardr[97] + 1./2.*ardr[90] + 3*
   ardr[88] + ardr[83];
   ardr[72]=ardr[56]*ardr[72];
   ardr[83]= - 9*ardr[9] + ardr[172] + 193./8. + ardr[101];
   ardr[83]=ardr[20]*ardr[83];
   ardr[88]= - 3*ardr[9] + ardr[174] + 193./24. + ardr[196];
   ardr[88]=ardr[1]*ardr[88];
   ardr[83]=ardr[83] + ardr[88];
   ardr[72]=1./2.*ardr[83] + ardr[72];
   ardr[83]=1 - ardr[10];
   ardr[88]=ardr[83] + ardr[158];
   ardr[90]=ardr[20]*ardr[88];
   ardr[88]=ardr[1]*ardr[88];
   ardr[88]=3*ardr[90] + ardr[88];
   ardr[88]=ardr[11]*ardr[88];
   ardr[90]=ardr[91] + ardr[13] + ardr[93] - ardr[6];
   ardr[91]=ardr[20]*ardr[90];
   ardr[90]=ardr[1]*ardr[90];
   ardr[93]= - 11 + 21./2.*ardr[9];
   ardr[93]=ardr[20]*ardr[93];
   ardr[97]= - 11./3. + ardr[127];
   ardr[101]=ardr[1]*ardr[97];
   ardr[93]=ardr[93] + ardr[101];
   ardr[93]=MMZ*ardr[93];
   ardr[112]=1 + ardr[130];
   ardr[114]=ardr[20]*ardr[112];
   ardr[112]=ardr[1]*ardr[112];
   ardr[112]=3*ardr[114] + ardr[112];
   ardr[112]=ardr[12]*ardr[112];
   ardr[88]=ardr[112] + 1./2.*ardr[88] + ardr[93] + 3*ardr[91] + 
   ardr[90];
   ardr[88]=ardr[3]*ardr[88];
   ardr[63]=ardr[71] + ardr[68] + ardr[63] + 1./2.*ardr[72] + ardr[88];
   ardr[63]=ardr[136]*ardr[63];
   ardr[68]=ardr[22] - ardr[6];
   ardr[66]=ardr[68] + ardr[66];
   ardr[71]=7./4.*ardr[9] + 10./3. - 29./4.*ardr[18];
   ardr[71]=1./3.*ardr[71] + ardr[73];
   ardr[71]=ardr[11]*ardr[71];
   ardr[72]= - 4*ardr[9] + 86./3. - 25*ardr[18];
   ardr[72]=1./3.*ardr[72] - 4*ardr[19];
   ardr[72]=MMZ*ardr[72];
   ardr[73]= - ardr[9] + 5./3. - 4*ardr[18];
   ardr[73]=ardr[12]*ardr[73];
   ardr[88]= - 22./3. - ardr[33];
   ardr[88]=ardr[16]*ardr[88];
   ardr[72]=2*ardr[88] + 4./3.*ardr[73] + ardr[71] + ardr[66] + 
   ardr[72];
   ardr[72]=ardr[3]*ardr[72];
   ardr[88]= - ardr[14] + ardr[39] - ardr[25] + ardr[27];
   ardr[88]=MMZ*ardr[88];
   ardr[90]=ardr[8]*ardr[142];
   ardr[91]=ardr[12]*ardr[134];
   ardr[88]=ardr[91] + ardr[88] + ardr[90];
   ardr[88]=ardr[56]*ardr[88];
   ardr[90]=1 - ardr[19];
   ardr[90]=MMZ*ardr[90];
   ardr[91]=ardr[9] - 5./3. + ardr[129];
   ardr[93]=ardr[8]*ardr[91];
   ardr[73]=3./4.*ardr[88] + 1./3.*ardr[73] + 3./4.*ardr[90] + 1./3.*
   ardr[93];
   ardr[73]=ardr[56]*ardr[73];
   ardr[88]=5./8.*ardr[96];
   ardr[90]=5./4.*ardr[105];
   ardr[93]=ardr[90] + ardr[88] - 1223./144. + ardr[141];
   ardr[96]=1./2.*ardr[19];
   ardr[105]= - 2 + ardr[61];
   ardr[105]=ardr[3]*ardr[16]*ardr[105];
   ardr[112]=6*ardr[105];
   ardr[114]=ardr[112] + ardr[96] - 77./3.*ardr[18] + ardr[78] - 14./3.
    - ardr[33];
   ardr[114]=ardr[3]*ardr[114];
   ardr[120]=12*ardr[104];
   ardr[114]=ardr[114] + ardr[120];
   ardr[114]=MMt*ardr[114];
   ardr[121]= - 1./4.*ardr[40] + 3./4.*ardr[6] + ardr[149] + ardr[22];
   ardr[121]=ardr[23]*ardr[121];
   ardr[122]=5./64.*ardr[148];
   ardr[68]=1./2.*ardr[68] - ardr[40];
   ardr[68]=ardr[23]*ardr[68];
   ardr[68]=11./4. + ardr[68];
   ardr[68]=1./2.*EPAIR2*ardr[68];
   ardr[125]= - 15./128.*MMZ*ardr[54];
   ardr[126]=5./64.*ardr[152];
   ardr[129]=ardr[23]*ardr[33];
   ardr[129]=5./8.*ardr[54] + ardr[129];
   ardr[129]=1./4.*ardr[16]*ardr[129];
   ardr[72]=ardr[114] + ardr[72] + ardr[73] + ardr[129] + ardr[126] + 
   ardr[125] + ardr[169] + 1./16.*ardr[9] + 25./16.*ardr[18] + ardr[68]
    + ardr[122] + 1./8.*ardr[93] + ardr[121];
   ardr[72]=ardr[15]*ardr[72];
   ardr[73]= - ardr[14] + ardr[39] - ardr[53] + ardr[7];
   ardr[93]=ardr[20]*ardr[73];
   ardr[73]=ardr[1]*ardr[73];
   ardr[73]=3*ardr[93] + ardr[73];
   ardr[73]=MMZ*ardr[73];
   ardr[79]=ardr[8]*ardr[79];
   ardr[76]=ardr[12]*ardr[76];
   ardr[73]=ardr[76] + ardr[73] + ardr[79];
   ardr[73]=ardr[56]*ardr[73];
   ardr[76]=ardr[20]*ardr[83];
   ardr[79]=ardr[1]*ardr[83];
   ardr[76]=3*ardr[76] + ardr[79];
   ardr[76]=MMZ*ardr[76];
   ardr[79]= - 1./3. + ardr[9];
   ardr[83]=ardr[20]*ardr[79];
   ardr[79]=ardr[1]*ardr[79];
   ardr[79]=5./3.*ardr[83] + ardr[79];
   ardr[83]=ardr[8]*ardr[79];
   ardr[93]=1./3. - ardr[9];
   ardr[114]=ardr[20]*ardr[93];
   ardr[93]=ardr[1]*ardr[93];
   ardr[93]=5./3.*ardr[114] + ardr[93];
   ardr[93]=ardr[12]*ardr[93];
   ardr[73]=1./4.*ardr[73] + ardr[93] + 1./4.*ardr[76] + ardr[83];
   ardr[73]=ardr[56]*ardr[73];
   ardr[76]= - 2*ardr[10];
   ardr[83]=17./3. + ardr[76];
   ardr[81]=2./3.*ardr[83] + ardr[81];
   ardr[81]=ardr[1]*ardr[81];
   ardr[76]=43./9. + ardr[76];
   ardr[76]=2*ardr[76] - 29./3.*ardr[9];
   ardr[76]=ardr[20]*ardr[76];
   ardr[76]=ardr[76] + ardr[81];
   ardr[76]=MMZ*ardr[76];
   ardr[81]=ardr[173] + ardr[175];
   ardr[81]=ardr[11]*ardr[81];
   ardr[76]=4*ardr[93] + ardr[76] + ardr[81];
   ardr[76]=ardr[3]*ardr[76];
   ardr[83]=7./2.*EPAIR2;
   ardr[93]=13./2.*ardr[9] + ardr[153] - 43./9. + ardr[83];
   ardr[93]=ardr[20]*ardr[93];
   ardr[114]= - 7./2.*EPAIR2;
   ardr[130]=1./6.*ardr[10];
   ardr[131]=29./6.*ardr[9] + ardr[130] - 11./3. + ardr[114];
   ardr[131]=ardr[1]*ardr[131];
   ardr[93]=ardr[93] + ardr[131];
   ardr[62]=ardr[63] + ardr[72] + ardr[62] + ardr[76] + 1./4.*ardr[93]
    + ardr[73];
   ardr[62]=ardr[136]*ardr[62];
   ardr[63]=ardr[128] + 11 + ardr[65];
   ardr[63]=ardr[11]*ardr[63];
   ardr[65]= - 17*ardr[22] - 5*ardr[6];
   ardr[72]=35./4.*ardr[9] - 121./3. + 119./4.*ardr[18];
   ardr[72]=MMZ*ardr[72];
   ardr[73]=64./3. + 17*ardr[33];
   ardr[73]=ardr[16]*ardr[73];
   ardr[63]=ardr[73] + 1./2.*ardr[63] + 1./3.*ardr[72] + 11*ardr[13] + 
   1./2.*ardr[65] + 17*ardr[40];
   ardr[63]=ardr[3]*ardr[63];
   ardr[65]= - 11./3. + 17./2.*ardr[33];
   ardr[65]=3*ardr[105] + ardr[96] + 49./36.*ardr[18] + 1./3.*ardr[65]
    - 3./2.*ardr[17];
   ardr[65]=ardr[3]*ardr[65];
   ardr[65]=ardr[65] + 6*ardr[104];
   ardr[65]=MMt*ardr[65];
   ardr[70]=25./2.*ardr[74] + 25./4.*ardr[70] + 3989./216. - 17*
   ardr[33];
   ardr[67]=ardr[123] - ardr[26] + ardr[67];
   ardr[67]=ardr[138]*ardr[67];
   ardr[67]=75./32.*ardr[54] + 1./3.*ardr[67];
   ardr[67]=MMZ*ardr[67];
   ardr[72]= - ardr[23]*ardr[33];
   ardr[73]= - ardr[16]*ardr[138];
   ardr[72]=1./3.*ardr[73] - 25./4.*ardr[54] + 17./3.*ardr[72];
   ardr[72]=ardr[16]*ardr[72];
   ardr[73]=13./8.*ardr[26] - 2*ardr[22];
   ardr[73]=19./24.*ardr[40] + 1./3.*ardr[73] + 1./8.*ardr[6];
   ardr[73]=ardr[23]*ardr[73];
   ardr[63]=ardr[65] + 1./3.*ardr[63] + 1./8.*ardr[72] + 25./64.*
   ardr[115] + 1./4.*ardr[67] - 35./144.*ardr[9] - 119./144.*ardr[18]
    + 25./64.*ardr[102] + 1./16.*ardr[70] + ardr[73];
   ardr[63]=ardr[84]*ardr[63];
   ardr[65]=8*ardr[9] + 20./3. + ardr[86];
   ardr[65]=1./9.*ardr[65] - ardr[19];
   ardr[65]=MMZ*ardr[65];
   ardr[67]= - 10./9. - ardr[33];
   ardr[67]=ardr[16]*ardr[67];
   ardr[65]=2*ardr[67] + ardr[71] + ardr[66] + ardr[65];
   ardr[65]=ardr[3]*ardr[65];
   ardr[66]=ardr[90] + ardr[88] - 3029./432. + ardr[141];
   ardr[67]=ardr[112] + ardr[96] - 119./9.*ardr[18] + ardr[78] + 70./9.
    - ardr[33];
   ardr[67]=ardr[3]*ardr[67];
   ardr[67]=ardr[67] + ardr[120];
   ardr[67]=MMt*ardr[67];
   ardr[63]=ardr[63] + ardr[67] + ardr[65] + ardr[129] + ardr[126] + 
   ardr[125] + ardr[169] - 7./144.*ardr[9] + 161./144.*ardr[18] + 
   ardr[68] + ardr[122] + 1./8.*ardr[66] + ardr[121];
   ardr[63]=ardr[84]*ardr[63];
   ardr[65]=MMZ*ardr[91];
   ardr[65]=ardr[65] + 8*ardr[16];
   ardr[65]=ardr[3]*ardr[65];
   ardr[66]=1 + ardr[18];
   ardr[66]=MMt*ardr[3]*ardr[66];
   ardr[65]=ardr[65] + 8*ardr[66];
   ardr[63]=8./9.*ardr[65] + ardr[63];
   ardr[63]=ardr[15]*ardr[63];
   ardr[65]=ardr[87] + ardr[61] + 85./12. + ardr[37];
   ardr[65]=1./4.*ardr[65] + ardr[80];
   ardr[65]=ardr[11]*ardr[65];
   ardr[66]=ardr[98]*ardr[103];
   ardr[66]=ardr[66] + 3./2.*ardr[107];
   ardr[66]=ardr[3]*ardr[66];
   ardr[67]= - 151./3. + 7./2.*ardr[37];
   ardr[67]=7*ardr[31] - 17./6.*ardr[35] + 1./6.*ardr[67] + ardr[32];
   ardr[67]=MMZ*ardr[67];
   ardr[68]= - 61./3. + ardr[35];
   ardr[68]=ardr[12]*ardr[68];
   ardr[65]=ardr[66] + 1./4.*ardr[68] + ardr[65] + 1./4.*ardr[116] + 1./
   2.*ardr[67] + 1./2.*ardr[14] + 19./4.*ardr[13] + 5./4.*ardr[39] + 
   ardr[95] - ardr[6] - 1./4.*ardr[48] + ardr[50];
   ardr[65]=ardr[3]*ardr[65];
   ardr[66]=15./2.*ardr[92] - 1381./144. + ardr[82];
   ardr[66]= - 49./6.*ardr[31] + 5./24.*ardr[35] + 13./4.*ardr[32] + 15.
   /2.*ardr[110] - 7./24.*ardr[37] + 15./4.*ardr[108] + 1./2.*ardr[66]
    - ardr[36];
   ardr[67]= - 7./2. + ardr[189];
   ardr[67]=ardr[8]*ardr[67];
   ardr[68]=3 + 1./4.*ardr[35];
   ardr[68]=ardr[12]*ardr[68];
   ardr[67]=1./4.*ardr[68] - 9./8.*ardr[11] + 3./8.*ardr[67] + 9./8.*
   MMZ - 19./16.*ardr[14] + ardr[164] + 5./8.*ardr[39] + 9./16.*
   ardr[43] - 3./16.*ardr[46] - ardr[50] + ardr[52] + 3./16.*ardr[48];
   ardr[67]=ardr[56]*ardr[67];
   ardr[68]= - MMZ*ardr[35];
   ardr[70]= - 1./2.*ardr[35];
   ardr[71]= - 1./3. + ardr[70];
   ardr[71]=ardr[11]*ardr[71];
   ardr[68]=1./3.*ardr[12] + 1./3.*ardr[68] + ardr[71];
   ardr[68]=ardr[3]*ardr[68];
   ardr[68]=1./24.*ardr[35] + ardr[68];
   ardr[68]=ardr[84]*ardr[68];
   ardr[71]= - MMZ*ardr[58];
   ardr[65]=1./4.*ardr[68] + ardr[65] + ardr[67] + ardr[85] + 15./16.*
   ardr[94] + 45./8.*ardr[71] + 1./4.*ardr[66] + 1./3.*ardr[34];
   ardr[65]=ardr[84]*ardr[65];
   ardr[60]=ardr[150] + ardr[60] + ardr[190] + ardr[61] + 35./4.*
   ardr[37] + 77./2. + ardr[89];
   ardr[60]=ardr[11]*ardr[60];
   ardr[66]=2 + ardr[189];
   ardr[66]=ardr[98]*ardr[66];
   ardr[67]=ardr[12]*ardr[106];
   ardr[66]=3*ardr[67] + ardr[66] + 3*ardr[107];
   ardr[66]=ardr[3]*ardr[66];
   ardr[67]= - ardr[34] + ardr[144] + ardr[75] + ardr[32] + 31./6.*
   ardr[37] - 7 + ardr[36];
   ardr[67]=MMZ*ardr[67];
   ardr[68]=19./2.*ardr[35];
   ardr[72]=ardr[68] + 19./4. + ardr[61];
   ardr[72]=ardr[12]*ardr[72];
   ardr[73]=3*ardr[13] + ardr[39] - ardr[48] - ardr[43];
   ardr[73]=1./2.*ardr[73] + ardr[14];
   ardr[60]=ardr[66] + 1./2.*ardr[72] + 1./2.*ardr[60] + ardr[117] + 9*
   ardr[73] + ardr[67];
   ardr[60]=ardr[3]*ardr[60];
   ardr[66]=ardr[113] + ardr[70];
   ardr[66]=ardr[8]*MMZ*ardr[66];
   ardr[67]= - 9*ardr[32];
   ardr[70]=1./2.*ardr[35];
   ardr[72]=ardr[67] + ardr[70];
   ardr[72]=ardr[12]*MMZ*ardr[72];
   ardr[73]= - 17./2.*ardr[14] + 17./2.*ardr[39] - 9*ardr[43] + 1./2.*
   ardr[46] + 9*ardr[47] - 1./2.*ardr[48];
   ardr[73]=MMZ*ardr[73];
   ardr[66]=ardr[72] + ardr[73] + ardr[66];
   ardr[66]=ardr[56]*ardr[66];
   ardr[70]=ardr[70] + 71./4. + ardr[67];
   ardr[70]=MMZ*ardr[70];
   ardr[67]=ardr[69] - 151./16. + ardr[67];
   ardr[67]=ardr[8]*ardr[67];
   ardr[61]=ardr[68] + 23./2. + ardr[61];
   ardr[61]=ardr[12]*ardr[61];
   ardr[61]=1./4.*ardr[66] + 1./4.*ardr[61] - 53./16.*ardr[11] + 1./2.*
   ardr[67] + 1./4.*ardr[70] - 9./4.*ardr[14] - 53./16.*ardr[13] + 13./
   32.*ardr[39] + 83./32.*ardr[43] + ardr[109] - ardr[47] + 25./8.*
   ardr[48];
   ardr[61]=ardr[56]*ardr[61];
   ardr[66]=63./16.*ardr[92] - 1907./48. + ardr[82];
   ardr[66]=63./16.*ardr[110] - 119./48.*ardr[37] + 63./32.*ardr[108]
    + 1./2.*ardr[66] - ardr[36];
   ardr[60]=ardr[65] + ardr[60] + 1./2.*ardr[61] + 63./32.*ardr[181] + 
   63./64.*ardr[94] + 189./32.*ardr[71] + ardr[124] + ardr[119] - 23./8.
   *ardr[35] + 1./2.*ardr[66] + ardr[118];
   ardr[60]=ardr[84]*ardr[60];
   ardr[61]= - ardr[43]*ardr[57];
   ardr[65]=ardr[39]*ardr[57];
   ardr[66]=ardr[13]*ardr[57];
   ardr[67]= - MMZ*ardr[57];
   ardr[61]=3./32.*ardr[67] + ardr[100] + 1./32.*ardr[66] + 1./64.*
   ardr[65] + 21 + 1./64.*ardr[61];
   ardr[61]=MMZ*ardr[61];
   ardr[65]=MMZ*ardr[57];
   ardr[66]= - ardr[8]*ardr[57];
   ardr[67]=ardr[11]*ardr[57];
   ardr[66]=1./2.*ardr[67] + ardr[65] + 1./2.*ardr[66];
   ardr[66]=ardr[11]*ardr[66];
   ardr[67]=1 + ardr[64];
   ardr[67]=MMZ*ardr[67];
   ardr[67]=3*ardr[67] - 5./2.*ardr[14] + ardr[77] - ardr[46] + 3./2.*
   ardr[51] + ardr[48];
   ardr[67]=MMZ*ardr[67];
   ardr[68]= - 3./2. + ardr[35];
   ardr[68]=ardr[8]*MMZ*ardr[68];
   ardr[69]= - 3./2. - ardr[35];
   ardr[69]=ardr[12]*MMZ*ardr[69];
   ardr[67]=ardr[69] + ardr[67] + ardr[68];
   ardr[67]=ardr[56]*ardr[67];
   ardr[68]=1 + ardr[37];
   ardr[65]=3*ardr[68] + 1./64.*ardr[65];
   ardr[65]=ardr[8]*ardr[65];
   ardr[68]= - 1 - ardr[37];
   ardr[68]=ardr[12]*ardr[68];
   ardr[61]=3*ardr[67] + 3*ardr[68] + 1./32.*ardr[66] + ardr[61] + 
   ardr[65];
   ardr[61]=ardr[56]*ardr[61];
   ardr[65]=ardr[58]*ardr[57];
   ardr[66]=ardr[43]*ardr[65];
   ardr[67]= - ardr[58]*ardr[57];
   ardr[68]=ardr[39]*ardr[67];
   ardr[66]=ardr[66] + ardr[68];
   ardr[68]=ardr[13]*ardr[67];
   ardr[69]=MMZ*ardr[65];
   ardr[66]=3*ardr[69] + 1./2.*ardr[66] + ardr[68];
   ardr[66]=MMZ*ardr[66];
   ardr[68]=MMZ*ardr[67];
   ardr[65]=ardr[8]*ardr[65];
   ardr[67]=ardr[11]*ardr[67];
   ardr[65]=1./2.*ardr[67] + ardr[68] + 1./2.*ardr[65];
   ardr[65]=ardr[11]*ardr[65];
   ardr[67]=11./3. - 2*ardr[64];
   ardr[64]= - 53./3. - 4*ardr[64];
   ardr[64]=ardr[37]*ardr[64];
   ardr[64]= - 6*ardr[35] + 2*ardr[67] + ardr[64];
   ardr[64]=MMZ*ardr[64];
   ardr[67]=2 + ardr[143];
   ardr[67]=ardr[3]*ardr[98]*ardr[67];
   ardr[69]= - 5 - 3*ardr[37];
   ardr[69]=ardr[12]*ardr[69];
   ardr[64]=ardr[67] + ardr[64] + 2*ardr[69];
   ardr[64]=ardr[3]*ardr[64];
   ardr[67]=ardr[8]*ardr[68];
   ardr[60]=ardr[60] + 2*ardr[64] + ardr[61] + 1./32.*ardr[65] + 1./64.
   *ardr[67] + 1./32.*ardr[66] + 81./4. + ardr[37];
   ardr[60]=ardr[41]*ardr[60];
   ardr[61]=ardr[20]*ardr[97];
   ardr[61]=11./9.*ardr[61] + ardr[101];
   ardr[61]=MMZ*ardr[61];
   ardr[64]=11./3.*ardr[111] + 3*ardr[99];
   ardr[64]=ardr[11]*ardr[64];
   ardr[65]= - ardr[6] + ardr[13];
   ardr[66]=ardr[20]*ardr[65];
   ardr[65]=ardr[1]*ardr[65];
   ardr[61]=ardr[64] + ardr[61] + 11./3.*ardr[66] + 3*ardr[65];
   ardr[61]=ardr[3]*ardr[61];
   ardr[64]=97./12. - 7*ardr[9];
   ardr[65]=ardr[20]*ardr[64];
   ardr[64]=ardr[1]*ardr[64];
   ardr[64]=11./9.*ardr[65] + ardr[64];
   ardr[61]=1./8.*ardr[64] + ardr[61];
   ardr[61]=ardr[84]*ardr[61];
   ardr[64]=ardr[127] + ardr[130] - 29./9. + ardr[114];
   ardr[64]=ardr[1]*ardr[64];
   ardr[65]=77./18.*ardr[9] + ardr[153] - 109./27. + ardr[83];
   ardr[65]=ardr[20]*ardr[65];
   ardr[64]=ardr[65] + ardr[64];
   ardr[65]= - 11./9.*ardr[9] + 20./27. - ardr[10];
   ardr[65]=ardr[20]*ardr[65];
   ardr[66]=4./3. - ardr[10];
   ardr[66]=1./3.*ardr[66] - ardr[9];
   ardr[66]=ardr[1]*ardr[66];
   ardr[65]=ardr[65] + ardr[66];
   ardr[65]=MMZ*ardr[65];
   ardr[65]=ardr[65] + ardr[81];
   ardr[65]=ardr[3]*ardr[65];
   ardr[61]=ardr[61] + 1./4.*ardr[64] + ardr[65];
   ardr[61]=ardr[84]*ardr[61];
   ardr[64]=ardr[3]*MMZ*ardr[79];

      drret = ardr[59] + ardr[60] + ardr[61] + ardr[62] + ardr[63] + 8./
      3.*ardr[64];
      return drret;
}
