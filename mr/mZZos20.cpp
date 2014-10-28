#include <ZZ.hpp>
std::complex<long double>
ZZ<MS>::m20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armZZos[266], mZZosret;

    armZZos[1]=double(nL + nH);
    armZZos[2]=pow(mmZ,-1);
    armZZos[3]=pow(s,-1);
    armZZos[4]=pow(c,-1);
    armZZos[5]=pow(mmH,-1);
    armZZos[6]=std::real(Tsil::B(0,0,mmZ,mu2));
    armZZos[7]=Tsil::I2(0,0,mmW,mu2);
    armZZos[8]=Tsil::I2(0,0,mmZ,mu2);
    armZZos[9]=Tsil::B(mmW,mmW,mmZ,mu2);
    armZZos[10]=Tsil::B(mmZ,mmH,mmZ,mu2);
    armZZos[11]=Tsil::B(0,mmH,mmZ,mu2);
    armZZos[12]=Tsil::Beps(mmW,mmW,mmZ,mu2);
    armZZos[13]=Tsil::Beps(mmZ,mmH,mmZ,mu2);
    armZZos[14]=Tsil::A(mmW,mu2);
    armZZos[15]=Tsil::A(mmZ,mu2);
    armZZos[16]=Tsil::A(mmH,mu2);
    armZZos[17]=Tsil::Aeps(mmW,mu2);
    armZZos[18]=Tsil::Aeps(mmZ,mu2);
    armZZos[19]=Tsil::Aeps(mmH,mu2);
    armZZos[20]=std::real(Tsil::B(0,mmW,mmZ,mu2));
    armZZos[21]=protW0W00->M(0);
    armZZos[22]=prot0000Z->M(0);
    armZZos[23]=prot0000W->M(0);
    armZZos[24]=prot00000->M(0);
    armZZos[25]=protHZ00->Uxzuv(0);
    armZZos[26]=protW0W00->Uzxyv(0);
    armZZos[27]=protHZ00->Txuv(0);
    armZZos[28]=prot0000Z->Tvxu(0);
    armZZos[29]=protW0W00->Txuv(0);
    armZZos[30]=double(nH);
    armZZos[31]=pow(mmt,-1);
    armZZos[32]=Tsil::B(mmt,mmt,mmZ,mu2);
    armZZos[33]=Tsil::A(mmt,mu2);
    armZZos[34]=double(nL);
    armZZos[35]=Tsil::I2(mmZ,mmt,mmt,mu2);
    armZZos[36]=Tsil::I2(mmH,mmt,mmt,mu2);
    armZZos[37]=Tsil::I2(0,mmW,mmt,mu2);
    armZZos[38]=Tsil::Beps(mmt,mmt,mmZ,mu2);
    armZZos[39]=Tsil::Aeps(mmt,mu2);
    armZZos[40]=protWtWt0->M(0);
    armZZos[41]=protW0W0t->M(0);
    armZZos[42]=prottZtHt->M(0);
    armZZos[43]=protttttH->M(0);
    armZZos[44]=protttttZ->M(0);
    armZZos[45]=prottttt0->M(0);
    armZZos[46]=prott0t0W->M(0);
    armZZos[47]=prottttt0->Vzxyv(0);
    armZZos[48]=prottZtHt->Uuyxv(0);
    armZZos[49]=protWtWt0->Uzxyv(0);
    armZZos[50]=protttttH->Uzxyv(0);
    armZZos[51]=protttttZ->Uzxyv(0);
    armZZos[52]=protWtWt0->Uuyxv(0);
    armZZos[53]=prottZtHt->Tuxv(0);
    armZZos[54]=protWtWt0->Tzyv(0);
    armZZos[55]=protWtWt0->Tyzv(0);
    armZZos[56]=prottZtHt->Txuv(0);
    armZZos[57]=prottZtHt->Suxv(0);
    armZZos[58]=prottZtHt->Svyz(0);
    armZZos[59]=protWtWt0->Svyz(0);
    armZZos[60]=prottttt0->Suxv(0);
    armZZos[61]=prottZtHt->Uyuzv(0);
    armZZos[62]=double(boson);
    armZZos[63]=Tsil::I2(mmW,mmW,mmZ,mu2);
    armZZos[64]=Tsil::I2(mmH,mmW,mmW,mu2);
    armZZos[65]=Tsil::I2(mmH,mmZ,mmZ,mu2);
    armZZos[66]=Tsil::I2(mmH,mmH,mmH,mu2);
    armZZos[67]=protZHHZZ->M(0);
    armZZos[68]=protZZHHH->M(0);
    armZZos[69]=protZWHWW->M(0);
    armZZos[70]=protWWWWH->M(0);
    armZZos[71]=protWWWWZ->M(0);
    armZZos[72]=protWWWW0->M(0);
    armZZos[73]=protZHHZZ->Uzxyv(0);
    armZZos[74]=protZWHWW->Uzxyv(0);
    armZZos[75]=protZZHHH->Uyuzv(0);
    armZZos[76]=protZHHZZ->Uuyxv(0);
    armZZos[77]=protZWHWW->Uuyxv(0);
    armZZos[78]=protZWHWW->Tzyv(0);
    armZZos[79]=protZWHWW->Tyzv(0);
    armZZos[80]=protZZHHH->Suxv(0);
    armZZos[81]=protZWHWW->Svyz(0);
    armZZos[82]=protWZWHW->Svyz(0);
    armZZos[83]=protWWWW0->Suxv(0);
    armZZos[84]=protWWWW0->Vzxyv(0);
    armZZos[85]=protZWHWW->Uxzuv(0);
    armZZos[86]=protZWHWW->Uyuzv(0);
    armZZos[87]=1/(mmt - mmW);
    armZZos[88]=1/(4*mmt - mmZ);
    armZZos[89]=1/(mmH - 4*mmZ);
    armZZos[90]=1/(mmH - 4*mmW);
   armZZos[91]=99./2.*armZZos[9];
   armZZos[92]=59./3. + armZZos[91];
   armZZos[93]=1./3.*armZZos[1];
   armZZos[94]= - 3*armZZos[10];
   armZZos[95]= - armZZos[1] - 3*armZZos[34];
   armZZos[96]=armZZos[6]*armZZos[95];
   armZZos[92]=armZZos[96] + armZZos[34] + armZZos[93] + 1./2.*
   armZZos[92] + armZZos[94];
   armZZos[92]=armZZos[32]*armZZos[88]*armZZos[92];
   armZZos[97]=2./3.*armZZos[10];
   armZZos[98]=1 - armZZos[10];
   armZZos[99]=armZZos[6]*armZZos[98];
   armZZos[100]=armZZos[32]*armZZos[98];
   armZZos[101]= - 1./2.*armZZos[48];
   armZZos[102]=1./2.*armZZos[100] + 1./2.*armZZos[99] + armZZos[97] + 
   armZZos[13] - 3*armZZos[53] + armZZos[101] + 3./2.*armZZos[11] - 1./
   2.*armZZos[25] - 25./24. + armZZos[27];
   armZZos[102]=armZZos[89]*armZZos[102];
   armZZos[103]=3*armZZos[13];
   armZZos[104]= - 99./4.*armZZos[9] + armZZos[103] - 9./2.*armZZos[56]
    - 3*armZZos[48] + 141./16.*armZZos[38] + 243./32.*armZZos[55] + 3./
   4.*armZZos[50] + 81./8.*armZZos[54] - 9./2.*armZZos[51] - 7./3. - 81.
   /16.*armZZos[52];
   armZZos[105]=armZZos[1] + 3*armZZos[34];
   armZZos[106]=armZZos[6]*armZZos[105];
   armZZos[107]=3*armZZos[10];
   armZZos[104]=1./2.*armZZos[106] - 1./2.*armZZos[34] - 1./6.*
   armZZos[1] + 1./2.*armZZos[104] + armZZos[107];
   armZZos[104]=armZZos[88]*armZZos[104];
   armZZos[108]=pow(armZZos[87],2);
   armZZos[109]= - armZZos[6]*armZZos[108];
   armZZos[109]=armZZos[108] + armZZos[109];
   armZZos[109]=armZZos[33]*armZZos[109];
   armZZos[110]=armZZos[54] - armZZos[55];
   armZZos[110]=armZZos[108]*armZZos[110];
   armZZos[111]=mmZ*armZZos[110];
   armZZos[112]= - 1./2.*armZZos[55] - 33./16. - armZZos[54];
   armZZos[112]=armZZos[87]*armZZos[112];
   armZZos[113]= - armZZos[39]*armZZos[108];
   armZZos[114]=armZZos[17]*armZZos[108];
   armZZos[115]= - 3./2.*armZZos[38] + 1./2.*armZZos[51] - 9./8. + 
   armZZos[52];
   armZZos[115]=armZZos[31]*armZZos[115];
   armZZos[116]= - armZZos[14]*armZZos[108];
   armZZos[117]=armZZos[14]*armZZos[108];
   armZZos[118]=armZZos[87] + armZZos[117];
   armZZos[118]=armZZos[6]*armZZos[118];
   armZZos[119]=3*armZZos[6];
   armZZos[120]=3*armZZos[32] + 7 + armZZos[119];
   armZZos[120]=armZZos[5]*armZZos[120];
   armZZos[92]=3./4.*armZZos[111] + armZZos[102] + 1./2.*armZZos[120]
    + 3./4.*armZZos[109] + 1./16.*armZZos[92] + 1./8.*armZZos[104] + 3./
   4.*armZZos[118] + 3./4.*armZZos[116] + 1./2.*armZZos[115] + 3./4.*
   armZZos[114] + 3./4.*armZZos[113] + armZZos[112] + 4*armZZos[46] - 
   armZZos[22] + 6*armZZos[40] - armZZos[44] + 6*armZZos[41];
   armZZos[92]=mmZ*armZZos[92];
   armZZos[102]= - 2*armZZos[17];
   armZZos[104]=1./2.*armZZos[8];
   armZZos[111]=1./2.*armZZos[35];
   armZZos[112]=armZZos[102] - 3*armZZos[39] + armZZos[111] + 2*
   armZZos[37] + armZZos[104];
   armZZos[115]=3*armZZos[14];
   armZZos[120]= - 1./2.*armZZos[16];
   armZZos[121]=armZZos[115] + armZZos[120];
   armZZos[122]=armZZos[121] + 3./2.*armZZos[15];
   armZZos[123]=armZZos[6]*armZZos[122];
   armZZos[122]=armZZos[32]*armZZos[122];
   armZZos[124]= - 1./2.*armZZos[19];
   armZZos[125]= - 3*armZZos[18];
   armZZos[126]= - 9*armZZos[14];
   armZZos[127]=1./2.*armZZos[16];
   armZZos[128]=3*armZZos[33];
   armZZos[112]=armZZos[128] + 1./2.*armZZos[122] + 1./2.*armZZos[123]
    - 9./2.*armZZos[15] + armZZos[127] + armZZos[126] + armZZos[125] + 
   3*armZZos[112] + armZZos[124];
   armZZos[112]=armZZos[5]*armZZos[112];
   armZZos[122]=armZZos[120] + armZZos[15];
   armZZos[123]=armZZos[6]*armZZos[122];
   armZZos[129]=2 - armZZos[10];
   armZZos[129]=armZZos[33]*armZZos[129];
   armZZos[130]=armZZos[32]*armZZos[122];
   armZZos[131]=1./2.*armZZos[19];
   armZZos[132]=2*armZZos[129] + 1./2.*armZZos[130] + 1./2.*
   armZZos[123] - 2./3.*armZZos[15] + 4./3.*armZZos[16] - armZZos[18]
    + armZZos[131] + armZZos[111] + armZZos[104] - armZZos[57];
   armZZos[132]=armZZos[89]*armZZos[132];
   armZZos[133]=739./864. - 3*armZZos[49];
   armZZos[134]=5*armZZos[28];
   armZZos[135]=23./24.*armZZos[11];
   armZZos[136]=5./8.*armZZos[48];
   armZZos[137]=pow(Pi,2);
   armZZos[133]=armZZos[136] - 5./4.*armZZos[137] + 33*armZZos[12] + 
   1037./128.*armZZos[38] - 141./256.*armZZos[55] + armZZos[135] + 3./
   32.*armZZos[50] - 207./64.*armZZos[54] - 49./16.*armZZos[51] + 
   armZZos[25] - 657./128.*armZZos[52] + 11*armZZos[133] + armZZos[134]
   ;
   armZZos[138]= - 491./9. + armZZos[91];
   armZZos[139]=armZZos[14]*armZZos[31];
   armZZos[140]=13./8.*armZZos[34] - 5 + 13./24.*armZZos[1];
   armZZos[140]=armZZos[6]*armZZos[140];
   armZZos[141]=9./32.*armZZos[14] - armZZos[15];
   armZZos[141]=armZZos[88]*armZZos[141];
   armZZos[142]=armZZos[15]*armZZos[31];
   armZZos[138]=5./4.*armZZos[32] + armZZos[141] + armZZos[140] + 1./2.
   *armZZos[142] - 37./24.*armZZos[34] - 37./72.*armZZos[1] - 3./8.*
   armZZos[10] + 1./16.*armZZos[138] + armZZos[139];
   armZZos[138]=armZZos[32]*armZZos[138];
   armZZos[140]=3*armZZos[61];
   armZZos[141]= - 1 + armZZos[10];
   armZZos[143]=armZZos[32]*armZZos[141];
   armZZos[144]= - 2*armZZos[10];
   armZZos[145]=1./2.*armZZos[48];
   armZZos[146]=3*armZZos[53];
   armZZos[147]=1./2.*armZZos[143] + armZZos[144] - 7./2.*armZZos[13]
    - armZZos[56] + armZZos[146] + armZZos[145] + 7./4. + armZZos[140];
   armZZos[147]=armZZos[89]*armZZos[147];
   armZZos[148]=19./2. - armZZos[32];
   armZZos[149]= - armZZos[5]*armZZos[33];
   armZZos[150]=21*armZZos[149];
   armZZos[148]=1./2.*armZZos[148] + armZZos[150];
   armZZos[148]=armZZos[5]*armZZos[148];
   armZZos[151]=pow(armZZos[5],2);
   armZZos[152]= - mmt*armZZos[151];
   armZZos[153]= - 2*armZZos[42];
   armZZos[147]=18*armZZos[152] + armZZos[147] + 3*armZZos[148] + 
   armZZos[153] - 4*armZZos[46] - 7*armZZos[40] + 2*armZZos[44] - 7./2.
   *armZZos[41];
   armZZos[147]=mmt*armZZos[147];
   armZZos[148]=99*armZZos[9];
   armZZos[154]= - 149./3. + armZZos[148];
   armZZos[155]=1./6.*armZZos[1];
   armZZos[154]=1./2.*armZZos[96] + 1./2.*armZZos[34] + armZZos[155] + 
   1./8.*armZZos[154] + armZZos[94];
   armZZos[154]=armZZos[88]*armZZos[154];
   armZZos[156]=armZZos[6]*armZZos[87];
   armZZos[157]= - armZZos[31] + 27./8.*armZZos[88];
   armZZos[157]=armZZos[32]*armZZos[157];
   armZZos[154]=1./4.*armZZos[157] + 1./2.*armZZos[154] + armZZos[156]
    + 5./4.*armZZos[87] + 3*armZZos[31];
   armZZos[154]=armZZos[33]*armZZos[154];
   armZZos[157]=1./2.*armZZos[17];
   armZZos[158]=1./4.*armZZos[35];
   armZZos[159]= - armZZos[59] + 1./2.*armZZos[37];
   armZZos[160]=armZZos[157] + armZZos[158] + armZZos[159] - 1./2.*
   armZZos[58];
   armZZos[160]=armZZos[31]*armZZos[160];
   armZZos[161]=armZZos[59] - 1./2.*armZZos[37];
   armZZos[161]=9./8.*armZZos[161] + armZZos[58];
   armZZos[162]= - 1./2.*armZZos[57];
   armZZos[163]=1./4.*armZZos[36];
   armZZos[161]=armZZos[163] + 3*armZZos[161] + armZZos[162];
   armZZos[164]= - 5./4.*armZZos[18];
   armZZos[161]=7./8.*armZZos[15] + 3./8.*armZZos[16] - 45./16.*
   armZZos[14] + armZZos[164] - 59./32.*armZZos[39] + 3./4.*
   armZZos[161] - armZZos[35];
   armZZos[161]=armZZos[88]*armZZos[161];
   armZZos[165]= - armZZos[14]*armZZos[87];
   armZZos[166]= - 37./4.*armZZos[34] - 37./12.*armZZos[1] - 43./9. + 
   armZZos[165];
   armZZos[167]=armZZos[34] + 5./8. + armZZos[93];
   armZZos[167]=armZZos[6]*armZZos[167];
   armZZos[166]=1./4.*armZZos[166] + armZZos[167];
   armZZos[166]=armZZos[6]*armZZos[166];
   armZZos[167]=armZZos[38] - 3 - armZZos[50];
   armZZos[167]= - armZZos[13] + 5./2.*armZZos[56] - 3./2.*armZZos[53]
    + 3./4.*armZZos[167] + armZZos[48];
   armZZos[167]=1./2.*armZZos[167] - armZZos[10];
   armZZos[167]=armZZos[88]*armZZos[167];
   armZZos[168]=1./2. + armZZos[10];
   armZZos[168]=armZZos[32]*armZZos[88]*armZZos[168];
   armZZos[167]=armZZos[167] + 1./2.*armZZos[168];
   armZZos[167]=mmH*armZZos[167];
   armZZos[169]= - armZZos[87]*armZZos[59];
   armZZos[170]=armZZos[39]*armZZos[87];
   armZZos[171]=armZZos[17]*armZZos[87];
   armZZos[172]=armZZos[18]*armZZos[31];
   armZZos[173]=5*armZZos[87] + armZZos[31];
   armZZos[173]=armZZos[14]*armZZos[173];
   armZZos[174]=119./144. - armZZos[137];
   armZZos[175]=armZZos[1]*armZZos[174];
   armZZos[174]=armZZos[34]*armZZos[174];
   armZZos[92]=armZZos[147] + 1./8.*armZZos[167] + armZZos[92] + 
   armZZos[132] + armZZos[112] + armZZos[154] + 1./2.*armZZos[138] + 1./
   2.*armZZos[161] + armZZos[166] + 1./4.*armZZos[142] + armZZos[174]
    + 1./3.*armZZos[175] - 47./24.*armZZos[10] + 1./2.*armZZos[173] + 
   1837./64.*armZZos[9] + 1./4.*armZZos[172] - 13./16.*armZZos[13] + 
   armZZos[160] + 1./2.*armZZos[171] + 7./4.*armZZos[170] + 3./2.*
   armZZos[169] - 9./32.*armZZos[56] + 1./2.*armZZos[133] - armZZos[53]
   ;
   armZZos[112]=pow(armZZos[3],2);
   armZZos[92]=armZZos[112]*armZZos[92];
   armZZos[132]= - 297*armZZos[9];
   armZZos[133]= - 1685./9. + armZZos[132];
   armZZos[138]=5*armZZos[10];
   armZZos[133]=1./4.*armZZos[133] + armZZos[138];
   armZZos[147]= - 1./9.*armZZos[1];
   armZZos[154]= - 1./3.*armZZos[34];
   armZZos[160]=armZZos[93] + armZZos[34];
   armZZos[161]=armZZos[6]*armZZos[160];
   armZZos[133]=armZZos[161] + armZZos[154] + 1./8.*armZZos[133] + 
   armZZos[147];
   armZZos[133]=armZZos[88]*armZZos[133];
   armZZos[166]=armZZos[32]*armZZos[88];
   armZZos[133]=armZZos[133] + 1./3.*armZZos[166];
   armZZos[133]=armZZos[32]*armZZos[133];
   armZZos[167]=2303./9. + 225*armZZos[52];
   armZZos[167]= - 1423./24.*armZZos[38] - 675./16.*armZZos[55] - 5./2.
   *armZZos[50] - 531./8.*armZZos[54] + 1./8.*armZZos[167] + 39*
   armZZos[51];
   armZZos[167]=297./4.*armZZos[9] - 5*armZZos[13] + 15./2.*armZZos[56]
    + 1./2.*armZZos[167] + 5*armZZos[48];
   armZZos[167]=1./2.*armZZos[167] - 5*armZZos[10];
   armZZos[173]=1./9.*armZZos[1];
   armZZos[174]= - 1./3.*armZZos[1];
   armZZos[175]=armZZos[174] - armZZos[34];
   armZZos[176]=armZZos[6]*armZZos[175];
   armZZos[167]=armZZos[176] + 1./3.*armZZos[34] + 1./4.*armZZos[167]
    + armZZos[173];
   armZZos[167]=armZZos[88]*armZZos[167];
   armZZos[177]=armZZos[6]*armZZos[108];
   armZZos[177]= - armZZos[108] + armZZos[177];
   armZZos[177]=armZZos[33]*armZZos[177];
   armZZos[178]=1./3. - armZZos[6];
   armZZos[179]= - 3*armZZos[32];
   armZZos[180]=armZZos[178] + armZZos[179];
   armZZos[180]=armZZos[5]*armZZos[180];
   armZZos[181]=1./3.*armZZos[143] + 1./3.*armZZos[99] + 2*armZZos[53]
    + 1./3.*armZZos[48] + armZZos[11] - 1./3.*armZZos[25] + 9./4. + 2./
   3.*armZZos[27];
   armZZos[181]=armZZos[89]*armZZos[181];
   armZZos[182]= - armZZos[54] + armZZos[55];
   armZZos[182]=mmZ*armZZos[108]*armZZos[182];
   armZZos[183]= - 4*armZZos[45] - armZZos[24];
   armZZos[183]=1./3.*armZZos[183] + 10*armZZos[44];
   armZZos[184]=17./3.*armZZos[55] + 23 + 25./3.*armZZos[54];
   armZZos[184]=armZZos[87]*armZZos[184];
   armZZos[185]=armZZos[39]*armZZos[108];
   armZZos[108]= - armZZos[17]*armZZos[108];
   armZZos[186]=13./6.*armZZos[38] + 1./2.*armZZos[54] - 5./6.*
   armZZos[51] + 17./8. - 4./3.*armZZos[52];
   armZZos[186]=armZZos[31]*armZZos[186];
   armZZos[187]= - armZZos[87] + armZZos[116];
   armZZos[187]=armZZos[6]*armZZos[187];
   armZZos[108]=7./4.*armZZos[182] + armZZos[181] + armZZos[180] + 
   armZZos[177] + armZZos[133] + armZZos[167] + armZZos[187] + 
   armZZos[117] + armZZos[186] + armZZos[108] + armZZos[185] + 1./4.*
   armZZos[184] - 8*armZZos[46] + 2./3.*armZZos[22] - 16*armZZos[40] + 
   1./3.*armZZos[183] - 12*armZZos[41];
   armZZos[108]=mmZ*armZZos[108];
   armZZos[117]= - 7279./9. - 1619*armZZos[9];
   armZZos[133]= - armZZos[14]*armZZos[31];
   armZZos[117]=1./16.*armZZos[117] + 5*armZZos[133];
   armZZos[133]=5./4.*armZZos[10];
   armZZos[117]=1./3.*armZZos[117] + armZZos[133];
   armZZos[154]=armZZos[154] + 7./2. + armZZos[147];
   armZZos[154]=armZZos[6]*armZZos[154];
   armZZos[167]= - armZZos[15]*armZZos[31];
   armZZos[177]=5./6.*armZZos[167];
   armZZos[180]= - 43./32.*armZZos[14] + 11./3.*armZZos[15];
   armZZos[180]=armZZos[88]*armZZos[180];
   armZZos[117]= - 7./4.*armZZos[32] + armZZos[180] + armZZos[154] + 
   armZZos[177] + 13./9.*armZZos[34] + 1./2.*armZZos[117] + 13./27.*
   armZZos[1];
   armZZos[117]=armZZos[32]*armZZos[117];
   armZZos[132]=1795./9. + armZZos[132];
   armZZos[154]= - 4./9.*armZZos[1];
   armZZos[132]=4*armZZos[161] - 4./3.*armZZos[34] + armZZos[154] + 1./
   8.*armZZos[132] + armZZos[138];
   armZZos[132]=armZZos[88]*armZZos[132];
   armZZos[161]= - armZZos[6]*armZZos[87];
   armZZos[180]=5*armZZos[31];
   armZZos[182]=armZZos[180] - 143./8.*armZZos[88];
   armZZos[182]=armZZos[32]*armZZos[182];
   armZZos[132]=1./6.*armZZos[182] + armZZos[132] + 1./3.*armZZos[161]
    - 8./3.*armZZos[87] - 9*armZZos[31];
   armZZos[132]=armZZos[33]*armZZos[132];
   armZZos[161]= - armZZos[35] + armZZos[8] + 2*armZZos[57];
   armZZos[182]= - 2 + armZZos[10];
   armZZos[182]=armZZos[33]*armZZos[182];
   armZZos[183]=armZZos[127] - armZZos[15];
   armZZos[184]=armZZos[32]*armZZos[183];
   armZZos[161]=4./3.*armZZos[182] + 1./3.*armZZos[184] + 1./3.*
   armZZos[123] + 1./3.*armZZos[161] - armZZos[16];
   armZZos[161]=armZZos[89]*armZZos[161];
   armZZos[182]= - armZZos[38] + 3 + armZZos[50];
   armZZos[185]= - 1./3.*armZZos[48];
   armZZos[186]=1./3.*armZZos[13];
   armZZos[187]=1./2.*armZZos[53];
   armZZos[182]=armZZos[186] - 5./6.*armZZos[56] + armZZos[187] + 1./4.
   *armZZos[182] + armZZos[185];
   armZZos[188]=1./3.*armZZos[10];
   armZZos[182]=1./2.*armZZos[182] + armZZos[188];
   armZZos[182]=armZZos[88]*armZZos[182];
   armZZos[189]= - 1./2. - armZZos[10];
   armZZos[189]=armZZos[32]*armZZos[88]*armZZos[189];
   armZZos[182]=armZZos[182] + 1./6.*armZZos[189];
   armZZos[182]=5./4.*mmH*armZZos[182];
   armZZos[189]=4./3.*armZZos[47] + armZZos[45];
   armZZos[190]= - 7*armZZos[44];
   armZZos[189]=11*armZZos[41] + 2*armZZos[189] + armZZos[190];
   armZZos[191]=2*armZZos[40];
   armZZos[189]=2./3.*armZZos[42] + 2*armZZos[46] + 1./3.*armZZos[189]
    + armZZos[191];
   armZZos[149]=42*armZZos[149];
   armZZos[192]=armZZos[149] - 29./2. - 9*armZZos[32];
   armZZos[192]=armZZos[5]*armZZos[192];
   armZZos[143]=11./3.*armZZos[143] + 4./3.*armZZos[10] - 29./3.*
   armZZos[13] + 2./3.*armZZos[56] + 22*armZZos[53] + 11./3.*
   armZZos[48] + 65./6. + 6*armZZos[61];
   armZZos[143]=armZZos[89]*armZZos[143];
   armZZos[193]=12*armZZos[152];
   armZZos[189]=armZZos[193] + armZZos[143] + 2*armZZos[189] + 
   armZZos[192];
   armZZos[189]=mmt*armZZos[189];
   armZZos[192]=2*armZZos[39];
   armZZos[194]=armZZos[192] + armZZos[8] - armZZos[35];
   armZZos[195]= - 1./3.*armZZos[16];
   armZZos[196]=armZZos[15] - armZZos[14] + armZZos[195];
   armZZos[196]=armZZos[6]*armZZos[196];
   armZZos[197]= - 5*armZZos[14];
   armZZos[198]= - armZZos[15] + armZZos[197] + 1./3.*armZZos[16];
   armZZos[198]=armZZos[32]*armZZos[198];
   armZZos[196]= - 6*armZZos[33] + 1./2.*armZZos[198] + 1./2.*
   armZZos[196] + armZZos[194] + 4*armZZos[14];
   armZZos[196]=armZZos[5]*armZZos[196];
   armZZos[198]=armZZos[14]*armZZos[87];
   armZZos[198]=armZZos[198] - 139./9. - 47*armZZos[9];
   armZZos[199]= - 2./3.*armZZos[1];
   armZZos[200]= - 2*armZZos[34] - 5./4. + armZZos[199];
   armZZos[200]=armZZos[6]*armZZos[200];
   armZZos[198]=armZZos[200] + 23./3.*armZZos[34] + 1./4.*armZZos[198]
    + 23./9.*armZZos[1];
   armZZos[198]=armZZos[6]*armZZos[198];
   armZZos[200]=1./3.*armZZos[25];
   armZZos[201]= - 5./32.*armZZos[50];
   armZZos[202]=23./72.*armZZos[11];
   armZZos[203]= - 28*armZZos[12];
   armZZos[204]=5./12.*armZZos[137];
   armZZos[205]=7./24.*armZZos[48];
   armZZos[206]=2./3.*armZZos[53];
   armZZos[207]=15./16.*armZZos[56];
   armZZos[208]=armZZos[87]*armZZos[59];
   armZZos[209]= - armZZos[39]*armZZos[87];
   armZZos[210]= - armZZos[17]*armZZos[87];
   armZZos[211]=2*armZZos[59] - armZZos[37];
   armZZos[212]= - 5./2.*armZZos[17] - armZZos[39] - 5./2.*armZZos[35]
    + 4*armZZos[211] + 5*armZZos[58];
   armZZos[212]=armZZos[31]*armZZos[212];
   armZZos[213]= - armZZos[18]*armZZos[31];
   armZZos[214]= - 11*armZZos[31];
   armZZos[215]= - 37./2.*armZZos[87] + armZZos[214];
   armZZos[215]=armZZos[14]*armZZos[215];
   armZZos[216]=2*armZZos[137];
   armZZos[217]= - 11./3. + armZZos[216];
   armZZos[218]=armZZos[1]*armZZos[217];
   armZZos[217]=armZZos[34]*armZZos[217];
   armZZos[219]=75./8.*armZZos[159] - 13*armZZos[58];
   armZZos[219]= - 5./4.*armZZos[36] + 3*armZZos[219] + 5./2.*
   armZZos[57];
   armZZos[219]= - 59./24.*armZZos[15] - 5./8.*armZZos[16] + 67./8.*
   armZZos[14] + 73./12.*armZZos[18] + 307./32.*armZZos[39] + 1./4.*
   armZZos[219] + 11./3.*armZZos[35];
   armZZos[219]=armZZos[88]*armZZos[219];
   armZZos[220]= - 5./8.*armZZos[13];
   armZZos[221]= - 5./4.*armZZos[10];
   armZZos[92]=armZZos[92] + armZZos[189] + armZZos[182] + armZZos[108]
    + armZZos[161] + armZZos[196] + armZZos[132] + armZZos[117] + 
   armZZos[219] + 1./3.*armZZos[198] + armZZos[177] + 1./3.*
   armZZos[217] + 1./9.*armZZos[218] + armZZos[221] + 1./6.*
   armZZos[215] - 1463./32.*armZZos[9] + 5./6.*armZZos[213] + 
   armZZos[220] + 1./3.*armZZos[212] + 11./12.*armZZos[210] + 4./3.*
   armZZos[209] + 2*armZZos[208] + armZZos[207] + armZZos[206] + 
   armZZos[205] + armZZos[204] + armZZos[203] - 2383./384.*armZZos[38]
    + 535./768.*armZZos[55] + armZZos[202] + armZZos[201] + 5287./384.*
   armZZos[54] + 157./48.*armZZos[51] + armZZos[200] + 1315./384.*
   armZZos[52] - 5./3.*armZZos[28] + 10601./3456. + 28*armZZos[49];
   armZZos[92]=armZZos[112]*armZZos[92];
   armZZos[108]= - 311./9. + 1./2.*armZZos[52];
   armZZos[108]=1./2.*armZZos[50] + 1./4.*armZZos[108] - 25./3.*
   armZZos[51];
   armZZos[108]=185./72.*armZZos[38] + 1./3.*armZZos[108] - 1./16.*
   armZZos[55];
   armZZos[117]= - 1./2.*armZZos[56];
   armZZos[108]=1./36.*armZZos[9] + armZZos[186] + armZZos[117] + 1./2.
   *armZZos[108] + armZZos[185];
   armZZos[132]=11./9.*armZZos[34];
   armZZos[177]=armZZos[1] + armZZos[132];
   armZZos[185]=armZZos[6]*armZZos[177];
   armZZos[108]=1./6.*armZZos[185] - 11./162.*armZZos[34] - 1./18.*
   armZZos[1] + 1./2.*armZZos[108] + armZZos[188];
   armZZos[108]=armZZos[88]*armZZos[108];
   armZZos[189]=19./3. - 1./4.*armZZos[9];
   armZZos[196]= - armZZos[1] - 11./9.*armZZos[34];
   armZZos[198]=armZZos[6]*armZZos[196];
   armZZos[189]=armZZos[198] + 11./27.*armZZos[34] + armZZos[93] + 1./3.
   *armZZos[189] - armZZos[10];
   armZZos[189]=armZZos[32]*armZZos[88]*armZZos[189];
   armZZos[208]=5*armZZos[27];
   armZZos[209]= - 5./2.*armZZos[25] - 761./24. + armZZos[208];
   armZZos[99]=17./6.*armZZos[100] + 5./6.*armZZos[99] + 22./9.*
   armZZos[10] + 11./3.*armZZos[13] - 17*armZZos[53] - 17./6.*
   armZZos[48] + 1./3.*armZZos[209] + 5./2.*armZZos[11];
   armZZos[99]=armZZos[89]*armZZos[99];
   armZZos[209]= - 257*armZZos[44] - 17*armZZos[22];
   armZZos[212]= - 1./3.*armZZos[38] - 1./4. + 1./3.*armZZos[51];
   armZZos[212]=armZZos[31]*armZZos[212];
   armZZos[209]=1./9.*armZZos[209] + 107./4.*armZZos[212];
   armZZos[212]=5*armZZos[6];
   armZZos[215]=17*armZZos[32] + 77./3. + armZZos[212];
   armZZos[215]=armZZos[5]*armZZos[215];
   armZZos[99]=1./3.*armZZos[99] + 1./18.*armZZos[215] + 25./48.*
   armZZos[189] + 1./9.*armZZos[209] + 25./8.*armZZos[108];
   armZZos[99]=mmZ*armZZos[99];
   armZZos[108]= - 11*armZZos[18];
   armZZos[123]=34*armZZos[129] + 17./2.*armZZos[130] + 5./2.*
   armZZos[123] - 22./3.*armZZos[15] + 71./3.*armZZos[16] + 
   armZZos[108] + 11./2.*armZZos[19] + 17./2.*armZZos[35] + 5./2.*
   armZZos[8] - 17*armZZos[57];
   armZZos[123]=armZZos[89]*armZZos[123];
   armZZos[129]=7*armZZos[32];
   armZZos[189]=395./2. + armZZos[129];
   armZZos[150]=1./18.*armZZos[189] + armZZos[150];
   armZZos[150]=armZZos[5]*armZZos[150];
   armZZos[100]=7./18.*armZZos[100] - 34./9.*armZZos[10] - 47./18.*
   armZZos[13] - 17./9.*armZZos[56] - 7./3.*armZZos[53] - 7./18.*
   armZZos[48] - 25./36. + armZZos[140];
   armZZos[100]=armZZos[89]*armZZos[100];
   armZZos[140]=22*armZZos[44] - 17*armZZos[42];
   armZZos[100]=6*armZZos[152] + armZZos[100] + 2./9.*armZZos[140] + 
   armZZos[150];
   armZZos[100]=mmt*armZZos[100];
   armZZos[140]=armZZos[163] + armZZos[162] + 1./8.*armZZos[159] + 25./
   3.*armZZos[58];
   armZZos[150]=1./8.*armZZos[16];
   armZZos[140]=37./72.*armZZos[15] + armZZos[150] + 1./48.*armZZos[14]
    - 47./36.*armZZos[18] - 941./288.*armZZos[39] + 1./4.*armZZos[140]
    - 7./9.*armZZos[35];
   armZZos[140]=armZZos[88]*armZZos[140];
   armZZos[152]=421 + 197./4.*armZZos[9];
   armZZos[152]= - 6655./81.*armZZos[34] - 605./9.*armZZos[1] + 1./9.*
   armZZos[152] - 25*armZZos[10];
   armZZos[152]=197./12.*armZZos[185] + 1./4.*armZZos[152] + 107./9.*
   armZZos[142];
   armZZos[162]=1./32.*armZZos[14] - 7./3.*armZZos[15];
   armZZos[162]=armZZos[88]*armZZos[162];
   armZZos[152]=1285./108.*armZZos[32] + 1./2.*armZZos[152] + 25./3.*
   armZZos[162];
   armZZos[152]=armZZos[32]*armZZos[152];
   armZZos[162]= - 1./3.*armZZos[9];
   armZZos[163]= - 49 + armZZos[162];
   armZZos[155]=1./2.*armZZos[198] + 11./54.*armZZos[34] + armZZos[155]
    + 1./8.*armZZos[163] - armZZos[10];
   armZZos[155]=armZZos[88]*armZZos[155];
   armZZos[163]= - 107*armZZos[31];
   armZZos[189]=armZZos[163] + 4825./8.*armZZos[88];
   armZZos[189]=armZZos[32]*armZZos[189];
   armZZos[155]=1./18.*armZZos[189] + 107./3.*armZZos[31] + 25*
   armZZos[155];
   armZZos[155]=armZZos[33]*armZZos[155];
   armZZos[189]=5*armZZos[8] + 17*armZZos[35];
   armZZos[108]=11./6.*armZZos[16] + armZZos[108] - 11./6.*armZZos[19]
    + 1./2.*armZZos[189] - 17*armZZos[39];
   armZZos[189]=armZZos[195] + armZZos[15];
   armZZos[195]=armZZos[6]*armZZos[189];
   armZZos[189]=armZZos[32]*armZZos[189];
   armZZos[108]=17./9.*armZZos[33] + 17./12.*armZZos[189] + 5./12.*
   armZZos[195] + 1./3.*armZZos[108] - 11./2.*armZZos[15];
   armZZos[108]=armZZos[5]*armZZos[108];
   armZZos[189]= - 86381./192. + 85*armZZos[28];
   armZZos[189]=115./72.*armZZos[11] + 25./32.*armZZos[50] - 10721./432.
   *armZZos[51] + 5./3.*armZZos[25] + 1./27.*armZZos[189] + 25./128.*
   armZZos[52];
   armZZos[189]=61./72.*armZZos[48] - 85./324.*armZZos[137] + 82393./
   10368.*armZZos[38] + 1./3.*armZZos[189] - 25./256.*armZZos[55];
   armZZos[195]=5./2.*armZZos[9];
   armZZos[209]=1 + armZZos[195];
   armZZos[209]= - 1045./72.*armZZos[34] + 1./9.*armZZos[209] - 95./8.*
   armZZos[1];
   armZZos[132]=armZZos[132] + 17./72. + armZZos[1];
   armZZos[132]=armZZos[6]*armZZos[132];
   armZZos[132]=1./2.*armZZos[209] + 5./3.*armZZos[132];
   armZZos[132]=armZZos[6]*armZZos[132];
   armZZos[209]=1./3.*armZZos[38];
   armZZos[215]=armZZos[209] - 1 - 1./3.*armZZos[50];
   armZZos[215]= - 1./9.*armZZos[13] + 5./18.*armZZos[56] - 1./6.*
   armZZos[53] + 1./4.*armZZos[215] + 1./9.*armZZos[48];
   armZZos[215]=1./2.*armZZos[215] - 1./9.*armZZos[10];
   armZZos[215]=armZZos[88]*armZZos[215];
   armZZos[168]=armZZos[215] + 1./18.*armZZos[168];
   armZZos[168]=mmH*armZZos[168];
   armZZos[215]=armZZos[39] - armZZos[58] + armZZos[111];
   armZZos[215]=armZZos[31]*armZZos[215];
   armZZos[217]=1183./144. - 5*armZZos[137];
   armZZos[218]=armZZos[1]*armZZos[217];
   armZZos[217]=armZZos[34]*armZZos[217];
   armZZos[99]=armZZos[100] + 25./8.*armZZos[168] + armZZos[99] + 1./9.
   *armZZos[123] + armZZos[108] + 1./6.*armZZos[155] + 1./6.*
   armZZos[152] + 25./6.*armZZos[140] + 1./3.*armZZos[132] + 107./108.*
   armZZos[142] + 11./81.*armZZos[217] + 1./9.*armZZos[218] - 391./216.
   *armZZos[10] - 655./5184.*armZZos[9] + 107./108.*armZZos[172] - 101./
   144.*armZZos[13] + 107./54.*armZZos[215] - 25./32.*armZZos[56] + 1./
   2.*armZZos[189] - 17./9.*armZZos[53];
   armZZos[100]=pow(armZZos[4],2);
   armZZos[99]=armZZos[100]*armZZos[99];
   armZZos[108]=7321./9. + 7*armZZos[52];
   armZZos[108]= - 5./8.*armZZos[54] + 1./8.*armZZos[108] + 275./3.*
   armZZos[51];
   armZZos[123]= - 1./2.*armZZos[50];
   armZZos[108]= - 5915./216.*armZZos[38] - 7./16.*armZZos[55] + 1./3.*
   armZZos[108] + armZZos[123];
   armZZos[108]=37./36.*armZZos[9] - armZZos[13] + 3./2.*armZZos[56] + 
   1./2.*armZZos[108] + armZZos[48];
   armZZos[108]=1./2.*armZZos[108] - armZZos[10];
   armZZos[108]=1./3.*armZZos[198] + 11./81.*armZZos[34] + 1./4.*
   armZZos[108] + armZZos[173];
   armZZos[108]=armZZos[88]*armZZos[108];
   armZZos[132]= - 37*armZZos[9];
   armZZos[140]= - 1483./3. + armZZos[132];
   armZZos[140]=1./36.*armZZos[140] + armZZos[10];
   armZZos[140]=1./3.*armZZos[185] - 11./81.*armZZos[34] + 1./8.*
   armZZos[140] + armZZos[147];
   armZZos[140]=armZZos[88]*armZZos[140];
   armZZos[140]=armZZos[140] + 5./27.*armZZos[166];
   armZZos[140]=armZZos[32]*armZZos[140];
   armZZos[147]= - 5*armZZos[24];
   armZZos[152]=22*armZZos[22] + 526*armZZos[44] - 68*armZZos[45] + 
   armZZos[147];
   armZZos[155]=armZZos[209] + 1./4. - 1./3.*armZZos[51];
   armZZos[155]=armZZos[31]*armZZos[155];
   armZZos[152]=173./2.*armZZos[155] + 1./9.*armZZos[152] - armZZos[46]
   ;
   armZZos[155]=armZZos[6] - armZZos[32];
   armZZos[155]=armZZos[5]*armZZos[155];
   armZZos[108]=armZZos[181] + 1./3.*armZZos[155] + 5*armZZos[140] + 1./
   9.*armZZos[152] + 5*armZZos[108];
   armZZos[108]=mmZ*armZZos[108];
   armZZos[140]= - 33941./9. + 223*armZZos[9];
   armZZos[139]=1./48.*armZZos[140] + 11*armZZos[139];
   armZZos[133]=1./3.*armZZos[139] + armZZos[133];
   armZZos[139]=173./54.*armZZos[167];
   armZZos[140]= - 55./27.*armZZos[34] + 1./2. - 5./3.*armZZos[1];
   armZZos[140]=armZZos[6]*armZZos[140];
   armZZos[152]= - 103./32.*armZZos[14] + 71./3.*armZZos[15];
   armZZos[152]=armZZos[88]*armZZos[152];
   armZZos[133]= - 1015./324.*armZZos[32] + 5./9.*armZZos[152] + 1./3.*
   armZZos[140] + armZZos[139] + 715./243.*armZZos[34] + 1./2.*
   armZZos[133] + 65./27.*armZZos[1];
   armZZos[133]=armZZos[32]*armZZos[133];
   armZZos[140]=68./27.*armZZos[47] + armZZos[45];
   armZZos[152]= - 1./3.*armZZos[46];
   armZZos[155]=2*armZZos[42];
   armZZos[140]=armZZos[155] + armZZos[152] + armZZos[41] + 2*
   armZZos[140] + armZZos[190];
   armZZos[168]= - 5./2. - armZZos[32];
   armZZos[149]=11./3.*armZZos[168] + armZZos[149];
   armZZos[149]=armZZos[5]*armZZos[149];
   armZZos[140]=armZZos[193] + armZZos[143] + 2./3.*armZZos[140] + 
   armZZos[149];
   armZZos[140]=mmt*armZZos[140];
   armZZos[132]=1631 + armZZos[132];
   armZZos[132]=4./3.*armZZos[185] - 44./81.*armZZos[34] + armZZos[154]
    + 1./72.*armZZos[132] + armZZos[10];
   armZZos[132]=armZZos[88]*armZZos[132];
   armZZos[143]=armZZos[156] - armZZos[87] - 173*armZZos[31];
   armZZos[149]=173./9.*armZZos[31] - 655./8.*armZZos[88];
   armZZos[149]=armZZos[32]*armZZos[149];
   armZZos[132]=1./6.*armZZos[149] + 1./9.*armZZos[143] + 5*
   armZZos[132];
   armZZos[132]=armZZos[33]*armZZos[132];
   armZZos[143]=5*armZZos[14] - armZZos[16];
   armZZos[143]=1./3.*armZZos[143] + armZZos[15];
   armZZos[143]=armZZos[6]*armZZos[143];
   armZZos[149]=17*armZZos[14];
   armZZos[156]=armZZos[149] + armZZos[16];
   armZZos[156]=1./3.*armZZos[156] - armZZos[15];
   armZZos[156]=armZZos[32]*armZZos[156];
   armZZos[143]= - 2./3.*armZZos[33] + 1./2.*armZZos[156] + 1./2.*
   armZZos[143] + armZZos[194] - 44./9.*armZZos[14];
   armZZos[143]=armZZos[5]*armZZos[143];
   armZZos[156]= - 41./27. + armZZos[9];
   armZZos[156]=5*armZZos[156] + 1./3.*armZZos[165];
   armZZos[165]= - 2*armZZos[1];
   armZZos[168]= - 22./9.*armZZos[34] - 55./36. + armZZos[165];
   armZZos[168]=armZZos[6]*armZZos[168];
   armZZos[156]=1./3.*armZZos[168] + 187./27.*armZZos[34] + 1./4.*
   armZZos[156] + 17./3.*armZZos[1];
   armZZos[156]=armZZos[6]*armZZos[156];
   armZZos[168]=21661./384. - 11*armZZos[28];
   armZZos[168]=565./384.*armZZos[54] + 9359./432.*armZZos[51] + 
   armZZos[25] + 5./27.*armZZos[168] + 163./128.*armZZos[52];
   armZZos[173]= - 1./2.*armZZos[35];
   armZZos[181]= - armZZos[39] + armZZos[58] + armZZos[173];
   armZZos[181]=173./9.*armZZos[181] + 11./2.*armZZos[17];
   armZZos[181]=armZZos[31]*armZZos[181];
   armZZos[189]=1./6.*armZZos[87] + armZZos[214];
   armZZos[189]=armZZos[14]*armZZos[189];
   armZZos[190]= - 109./9. + armZZos[216];
   armZZos[193]=armZZos[1]*armZZos[190];
   armZZos[190]=armZZos[34]*armZZos[190];
   armZZos[194]=7./8.*armZZos[159] - 275./3.*armZZos[58];
   armZZos[194]= - 1./4.*armZZos[36] + 1./3.*armZZos[194] + 1./2.*
   armZZos[57];
   armZZos[194]= - 311./216.*armZZos[15] - 1./8.*armZZos[16] + 31./72.*
   armZZos[14] + 541./108.*armZZos[18] + 9773./864.*armZZos[39] + 1./4.
   *armZZos[194] + 71./27.*armZZos[35];
   armZZos[194]=armZZos[88]*armZZos[194];
   armZZos[99]=armZZos[99] + armZZos[140] + armZZos[182] + armZZos[108]
    + armZZos[161] + armZZos[143] + armZZos[132] + armZZos[133] + 5*
   armZZos[194] + 1./3.*armZZos[156] + armZZos[139] + 11./81.*
   armZZos[190] + 1./9.*armZZos[193] + armZZos[221] + 1./6.*
   armZZos[189] - 1525./864.*armZZos[9] + 173./54.*armZZos[213] + 
   armZZos[220] + 1./3.*armZZos[181] + 1./36.*armZZos[210] + 1./9.*
   armZZos[170] + armZZos[207] + armZZos[206] + armZZos[205] + 55./324.
   *armZZos[137] - 68053./10368.*armZZos[38] - 59./2304.*armZZos[55] + 
   armZZos[202] + 1./3.*armZZos[168] + armZZos[201];
   armZZos[99]=armZZos[100]*armZZos[99];
   armZZos[108]=2*armZZos[17];
   armZZos[132]=armZZos[108] + 8*armZZos[39] - 2*armZZos[37] - 3*
   armZZos[36];
   armZZos[133]=6*armZZos[15];
   armZZos[139]=armZZos[32]*armZZos[16];
   armZZos[128]=armZZos[128] + 11./3.*armZZos[139] + armZZos[133] - 11./
   3.*armZZos[16] + 3*armZZos[132] + 38./3.*armZZos[19];
   armZZos[128]=armZZos[5]*armZZos[128];
   armZZos[132]= - 15*armZZos[49];
   armZZos[140]= - 11./2.*armZZos[52];
   armZZos[143]=armZZos[140] + 7963./324. + armZZos[132];
   armZZos[156]= - 2*armZZos[15];
   armZZos[161]=2*armZZos[33] + armZZos[156] + armZZos[19] - 
   armZZos[36] + armZZos[192];
   armZZos[161]=3*armZZos[89]*armZZos[161];
   armZZos[168]= - 2*armZZos[46] + 11*armZZos[43] - armZZos[41] + 
   armZZos[191];
   armZZos[181]= - 9*armZZos[5];
   armZZos[168]= - 6*armZZos[89] + 1./3.*armZZos[168] + armZZos[181];
   armZZos[168]=mmt*armZZos[168];
   armZZos[182]=4*armZZos[61];
   armZZos[189]=2*armZZos[51];
   armZZos[190]= - 11./3.*armZZos[50];
   armZZos[193]=73./24.*armZZos[38];
   armZZos[194]=15./4.*armZZos[12];
   armZZos[198]=11./3.*armZZos[53];
   armZZos[201]=145./36.*armZZos[56];
   armZZos[205]= - 4*armZZos[13];
   armZZos[206]=armZZos[6]*armZZos[9];
   armZZos[207]=1./3.*armZZos[206];
   armZZos[209]=4./3.*armZZos[32] - armZZos[6] + 23./3. + armZZos[9];
   armZZos[209]=2./3.*armZZos[32]*armZZos[209];
   armZZos[210]= - 9./4.*armZZos[43] - 4*armZZos[42];
   armZZos[210]=mmH*armZZos[210];
   armZZos[143]=armZZos[168] + armZZos[210] + armZZos[161] + 
   armZZos[128] + armZZos[209] + armZZos[207] + armZZos[144] + 
   armZZos[205] + armZZos[201] + armZZos[198] + armZZos[194] + 
   armZZos[193] - 439./144.*armZZos[55] + armZZos[190] - 67./9.*
   armZZos[54] + armZZos[189] + 1./4.*armZZos[143] + armZZos[182];
   armZZos[143]=mmt*armZZos[143];
   armZZos[213]= - 11*armZZos[9];
   armZZos[214]= - 11735./162. + armZZos[213];
   armZZos[215]=2*armZZos[10];
   armZZos[216]= - 1./3.*armZZos[6];
   armZZos[217]= - 29./36.*armZZos[32];
   armZZos[214]=armZZos[217] + armZZos[216] + 1./4.*armZZos[214] + 
   armZZos[215];
   armZZos[214]=armZZos[33]*armZZos[214];
   armZZos[218]=23./6.*armZZos[13];
   armZZos[219]= - 3*armZZos[61];
   armZZos[222]= - 11./3.*armZZos[53];
   armZZos[223]=armZZos[218] - 7./3.*armZZos[56] + armZZos[222] - 5./6.
   *armZZos[48] + 19./12.*armZZos[38] - 19./12.*armZZos[50] + 151./72.
    + armZZos[219];
   armZZos[224]= - 11 - 59./2.*armZZos[10];
   armZZos[224]=1./6.*armZZos[224] - armZZos[32];
   armZZos[224]=armZZos[32]*armZZos[224];
   armZZos[225]=mmH*armZZos[43];
   armZZos[223]=1./3.*armZZos[225] + 1./3.*armZZos[224] + 1./2.*
   armZZos[223] - 5./9.*armZZos[10];
   armZZos[223]=mmH*armZZos[223];
   armZZos[224]= - 4*armZZos[14];
   armZZos[225]= - 3*armZZos[16];
   armZZos[226]=armZZos[15] + armZZos[224] + armZZos[225];
   armZZos[226]=4*armZZos[226] - 3*armZZos[33];
   armZZos[226]=armZZos[5]*armZZos[33]*armZZos[226];
   armZZos[227]=151./16.*armZZos[37] + 4*armZZos[60] + 395./8.*
   armZZos[59];
   armZZos[228]= - 61*armZZos[16];
   armZZos[229]= - 475./27.*armZZos[15] - 1283./18.*armZZos[14] + 
   armZZos[228];
   armZZos[229]=armZZos[32]*armZZos[229];
   armZZos[230]=67./36.*armZZos[57];
   armZZos[231]= - 7./3.*armZZos[17];
   armZZos[232]= - 107./18.*armZZos[19];
   armZZos[233]=11./9.*armZZos[16];
   armZZos[234]=armZZos[6]*armZZos[14];
   armZZos[235]=1./3.*armZZos[234];
   armZZos[236]=armZZos[156] - armZZos[33];
   armZZos[236]=3*armZZos[89]*armZZos[33]*armZZos[236];
   armZZos[237]=9./4.*armZZos[36];
   armZZos[143]=armZZos[143] + armZZos[223] + armZZos[236] + 
   armZZos[226] + armZZos[214] + 1./12.*armZZos[229] + armZZos[235] - 
   581./324.*armZZos[15] + armZZos[233] - 919./108.*armZZos[14] + 349./
   324.*armZZos[18] + armZZos[232] + armZZos[231] - 2027./162.*
   armZZos[39] - 673./324.*armZZos[35] + armZZos[237] + armZZos[230] + 
   1./9.*armZZos[227] + armZZos[58];
   armZZos[143]=mmt*armZZos[143];
   armZZos[214]= - armZZos[32]*armZZos[16];
   armZZos[227]= - 3./2.*armZZos[36];
   armZZos[229]=armZZos[17] + 4*armZZos[39] - armZZos[37] + 
   armZZos[227];
   armZZos[229]=3*armZZos[229];
   armZZos[238]=3*armZZos[15];
   armZZos[239]=3./2.*armZZos[33];
   armZZos[214]=armZZos[239] + 7./18.*armZZos[214] + armZZos[238] + 7./
   18.*armZZos[16] + armZZos[229] + 37./9.*armZZos[19];
   armZZos[214]=armZZos[5]*armZZos[214];
   armZZos[240]= - 7./3.*armZZos[43] - armZZos[41] - armZZos[40];
   armZZos[241]= - 2./3.*armZZos[46];
   armZZos[240]=1./2.*armZZos[240] + armZZos[241];
   armZZos[242]= - 3*armZZos[89];
   armZZos[240]=armZZos[242] + 1./3.*armZZos[240] - 9./2.*armZZos[5];
   armZZos[240]=mmt*armZZos[240];
   armZZos[243]= - 1./2.*armZZos[36];
   armZZos[244]=armZZos[33] - armZZos[15] + armZZos[131] + armZZos[243]
    + armZZos[39];
   armZZos[244]=3*armZZos[89]*armZZos[244];
   armZZos[245]=1433./324. + 3*armZZos[49];
   armZZos[245]=5*armZZos[245] + 77./36.*armZZos[52];
   armZZos[246]= - 1./2.*armZZos[9];
   armZZos[247]= - 4./3.*armZZos[32] - 2./3.*armZZos[6] + 13./3. + 
   armZZos[246];
   armZZos[247]=armZZos[32]*armZZos[247];
   armZZos[248]=13./24.*armZZos[43] + armZZos[153];
   armZZos[248]=mmH*armZZos[248];
   armZZos[249]=2*armZZos[61];
   armZZos[250]= - 1./2.*armZZos[54];
   armZZos[251]= - 2*armZZos[13];
   armZZos[206]=armZZos[240] + armZZos[248] + armZZos[244] + 
   armZZos[214] + 1./3.*armZZos[247] + 1./6.*armZZos[206] - armZZos[10]
    + armZZos[251] + 595./216.*armZZos[56] - 7./18.*armZZos[53] - 15./4.
   *armZZos[12] - 277./144.*armZZos[38] + 967./288.*armZZos[55] + 7./18.
   *armZZos[50] + armZZos[250] + armZZos[51] + 1./4.*armZZos[245] + 
   armZZos[249];
   armZZos[206]=mmt*armZZos[206];
   armZZos[214]=29./18.*armZZos[13] - 31./9.*armZZos[56] + 67./9.*
   armZZos[53] + 25./18.*armZZos[48] + 97./36.*armZZos[38] - 97./36.*
   armZZos[50] + 1933./216. + armZZos[219];
   armZZos[214]=1./2.*armZZos[214] + 85./27.*armZZos[10];
   armZZos[240]=7 + 103./2.*armZZos[10];
   armZZos[240]=1./12.*armZZos[240] + armZZos[32];
   armZZos[240]=armZZos[32]*armZZos[240];
   armZZos[245]= - mmH*armZZos[43];
   armZZos[214]=1./9.*armZZos[245] + 1./2.*armZZos[214] + 1./9.*
   armZZos[240];
   armZZos[214]=mmH*armZZos[214];
   armZZos[240]=47*armZZos[59] + 275./6.*armZZos[37];
   armZZos[240]=661./36.*armZZos[57] + 1./8.*armZZos[240] - 7*
   armZZos[58];
   armZZos[240]=869./324.*armZZos[18] - 251./54.*armZZos[19] - 61./36.*
   armZZos[17] - 1957./648.*armZZos[39] - 113./324.*armZZos[35] + 1./3.
   *armZZos[240] + armZZos[237];
   armZZos[245]= - 2*armZZos[14];
   armZZos[225]=armZZos[15] + armZZos[245] + armZZos[225];
   armZZos[247]= - 3./2.*armZZos[33];
   armZZos[225]=2*armZZos[225] + armZZos[247];
   armZZos[225]=armZZos[5]*armZZos[33]*armZZos[225];
   armZZos[248]=17*armZZos[16];
   armZZos[252]= - 155./9.*armZZos[15] + 53./6.*armZZos[14] + 
   armZZos[248];
   armZZos[252]=armZZos[32]*armZZos[252];
   armZZos[253]= - 6551./162. + 15*armZZos[9];
   armZZos[253]=43./72.*armZZos[32] - 1./18.*armZZos[6] + 1./4.*
   armZZos[253] + armZZos[10];
   armZZos[253]=armZZos[33]*armZZos[253];
   armZZos[254]= - 1./2.*armZZos[33];
   armZZos[255]= - armZZos[15] + armZZos[254];
   armZZos[255]=3*armZZos[89]*armZZos[33]*armZZos[255];
   armZZos[206]=armZZos[206] + armZZos[214] + armZZos[255] + 
   armZZos[225] + armZZos[253] + 1./72.*armZZos[252] + 1./18.*
   armZZos[234] - 781./648.*armZZos[15] - 71./27.*armZZos[16] + 1./2.*
   armZZos[240] - 43./27.*armZZos[14];
   armZZos[206]=mmt*armZZos[206];
   armZZos[214]= - 971./8. + 5*armZZos[11];
   armZZos[225]= - armZZos[6]*armZZos[10];
   armZZos[240]= - armZZos[88]*armZZos[39];
   armZZos[252]=1./2. - armZZos[10];
   armZZos[252]=armZZos[33]*armZZos[88]*armZZos[252];
   armZZos[253]= - armZZos[32]*armZZos[10];
   armZZos[214]=25*armZZos[252] + 211./24.*armZZos[253] + 25./2.*
   armZZos[240] + 5./3.*armZZos[225] - 215./36.*armZZos[10] + 25./8.*
   armZZos[13] - 25./8.*armZZos[56] - 523./12.*armZZos[53] + 1./3.*
   armZZos[214] - 25./8.*armZZos[48];
   armZZos[214]=mmH*armZZos[214];
   armZZos[256]= - 25./12.*armZZos[18];
   armZZos[257]=215./54.*armZZos[15] + 5741./108.*armZZos[16] + 
   armZZos[256] + 227./18.*armZZos[39] + 25./12.*armZZos[35] - 455./9.*
   armZZos[57] + 43./2.*armZZos[36];
   armZZos[258]= - armZZos[16] + armZZos[15];
   armZZos[259]=armZZos[6]*armZZos[258];
   armZZos[260]=armZZos[32]*armZZos[258];
   armZZos[257]=211./72.*armZZos[260] + 1./2.*armZZos[257] + 5./9.*
   armZZos[259];
   armZZos[261]=1879./32. + 29*armZZos[10];
   armZZos[122]=armZZos[88]*armZZos[122];
   armZZos[262]= - armZZos[33]*armZZos[88];
   armZZos[261]=25./16.*armZZos[262] + 1./3.*armZZos[261] + 25./8.*
   armZZos[122];
   armZZos[261]=armZZos[33]*armZZos[261];
   armZZos[214]=1./24.*armZZos[214] + 1./8.*armZZos[257] + 1./3.*
   armZZos[261];
   armZZos[214]=mmH*armZZos[214];
   armZZos[257]=551./2.*armZZos[14] + 391*armZZos[16];
   armZZos[257]= - 265./3.*armZZos[33] + 1./2.*armZZos[257] - 431./3.*
   armZZos[15];
   armZZos[257]=armZZos[33]*armZZos[257];
   armZZos[214]=1./36.*armZZos[257] + armZZos[214];
   armZZos[257]= - 1 + armZZos[9];
   armZZos[261]=armZZos[33]*armZZos[257];
   armZZos[263]= - armZZos[12] + armZZos[55] + 1./4. + armZZos[49];
   armZZos[263]=mmt*armZZos[263];
   armZZos[261]=armZZos[263] + armZZos[39] + armZZos[261];
   armZZos[261]=armZZos[100]*mmt*armZZos[261];
   armZZos[206]=1./4.*armZZos[261] + 1./3.*armZZos[214] + armZZos[206];
   armZZos[206]=armZZos[100]*armZZos[206];
   armZZos[214]= - 17./8. + armZZos[11];
   armZZos[261]=armZZos[88]*armZZos[39];
   armZZos[263]= - 1./2. + armZZos[10];
   armZZos[263]=armZZos[33]*armZZos[88]*armZZos[263];
   armZZos[225]=1./3.*armZZos[225];
   armZZos[264]=armZZos[32]*armZZos[10];
   armZZos[136]=5*armZZos[263] + 23./24.*armZZos[264] + 5./2.*
   armZZos[261] + armZZos[225] + armZZos[221] + armZZos[220] + 5./8.*
   armZZos[56] - 1./12.*armZZos[53] + 1./3.*armZZos[214] + armZZos[136]
   ;
   armZZos[136]=mmH*armZZos[136];
   armZZos[214]=1./2.*armZZos[36];
   armZZos[220]=5./6.*armZZos[15] - 1./12.*armZZos[16] + 5./12.*
   armZZos[18] - 151./18.*armZZos[39] - 5./12.*armZZos[35] - 5./9.*
   armZZos[57] + armZZos[214];
   armZZos[221]=armZZos[16] - armZZos[15];
   armZZos[261]=armZZos[32]*armZZos[221];
   armZZos[220]=23./72.*armZZos[261] + 1./2.*armZZos[220] + 1./9.*
   armZZos[259];
   armZZos[263]=19./16. + armZZos[144];
   armZZos[183]=armZZos[88]*armZZos[183];
   armZZos[265]=armZZos[33]*armZZos[88];
   armZZos[183]=5./8.*armZZos[265] + 7./3.*armZZos[263] + 5./4.*
   armZZos[183];
   armZZos[183]=armZZos[33]*armZZos[183];
   armZZos[136]=1./12.*armZZos[136] + 1./4.*armZZos[220] + 1./3.*
   armZZos[183];
   armZZos[136]=mmH*armZZos[136];
   armZZos[183]= - 23*armZZos[16];
   armZZos[220]= - 167./6.*armZZos[14] + armZZos[183];
   armZZos[220]= - 395./9.*armZZos[33] + 1./2.*armZZos[220] - 331./9.*
   armZZos[15];
   armZZos[220]=armZZos[33]*armZZos[220];
   armZZos[143]=armZZos[206] + armZZos[143] + 1./18.*armZZos[220] + 
   armZZos[136];
   armZZos[143]=armZZos[100]*armZZos[143];
   armZZos[132]=armZZos[140] + 289./12. + armZZos[132];
   armZZos[128]=armZZos[168] + armZZos[210] + armZZos[161] + 
   armZZos[128] + armZZos[209] + armZZos[207] + armZZos[144] + 
   armZZos[205] + armZZos[201] + armZZos[198] + armZZos[194] + 
   armZZos[193] - 301./48.*armZZos[55] + armZZos[190] - 32./3.*
   armZZos[54] + armZZos[189] + 1./4.*armZZos[132] + armZZos[182];
   armZZos[128]=mmt*armZZos[128];
   armZZos[132]= - 1295./18. + armZZos[213];
   armZZos[132]=armZZos[217] + armZZos[216] + 1./4.*armZZos[132] + 
   armZZos[215];
   armZZos[132]=armZZos[33]*armZZos[132];
   armZZos[140]=4./3.*armZZos[60];
   armZZos[161]=205./16.*armZZos[37] + armZZos[140] - 23./8.*
   armZZos[59];
   armZZos[168]=61./3.*armZZos[15] - 2305./6.*armZZos[14] + 
   armZZos[228];
   armZZos[168]=armZZos[32]*armZZos[168];
   armZZos[128]=armZZos[128] + armZZos[223] + armZZos[236] + 
   armZZos[226] + armZZos[132] + 1./12.*armZZos[168] + armZZos[235] + 
   163./36.*armZZos[15] + armZZos[233] - 1013./36.*armZZos[14] + 
   armZZos[256] + armZZos[232] + armZZos[231] - 113./6.*armZZos[39] + 
   13./12.*armZZos[35] + armZZos[237] + armZZos[230] + 1./3.*
   armZZos[161] + armZZos[58];
   armZZos[128]=mmt*armZZos[128];
   armZZos[132]=armZZos[239] + 1./2.*armZZos[139] + armZZos[238] + 
   armZZos[120] + armZZos[229] + 5*armZZos[19];
   armZZos[132]=armZZos[5]*armZZos[132];
   armZZos[139]=armZZos[181] + armZZos[43] + armZZos[41] - armZZos[40];
   armZZos[139]=1./2.*armZZos[139] + armZZos[242];
   armZZos[139]=mmt*armZZos[139];
   armZZos[161]=5./8.*armZZos[52] + 223./12. - 13*armZZos[49];
   armZZos[168]= - armZZos[6]*armZZos[9];
   armZZos[181]=1 + armZZos[246];
   armZZos[181]=armZZos[32]*armZZos[181];
   armZZos[153]= - 1./8.*armZZos[43] + armZZos[153];
   armZZos[153]=mmH*armZZos[153];
   armZZos[123]=armZZos[139] + armZZos[153] + armZZos[244] + 
   armZZos[132] + armZZos[181] + 1./2.*armZZos[168] - armZZos[10] + 
   armZZos[251] + 59./24.*armZZos[56] + armZZos[187] + 13./2.*
   armZZos[12] - 13./16.*armZZos[38] - 365./32.*armZZos[55] + 
   armZZos[123] + armZZos[250] + armZZos[51] + 1./2.*armZZos[161] + 
   armZZos[249];
   armZZos[123]=mmt*armZZos[123];
   armZZos[132]=129*armZZos[59] - 157./2.*armZZos[37];
   armZZos[139]= - 1./4.*armZZos[35];
   armZZos[153]=5./4.*armZZos[18];
   armZZos[132]= - 25./2.*armZZos[14] + armZZos[153] - 31./6.*
   armZZos[19] - 19./4.*armZZos[17] - 841./24.*armZZos[39] + 
   armZZos[139] + armZZos[237] + 53./12.*armZZos[57] + 1./8.*
   armZZos[132] - armZZos[58];
   armZZos[161]=5./2.*armZZos[13];
   armZZos[145]=armZZos[161] - 3*armZZos[56] + armZZos[146] + 
   armZZos[145] + 9./4.*armZZos[38] - 9./4.*armZZos[50] + 149./24. + 
   armZZos[219];
   armZZos[168]=5./3.*armZZos[10];
   armZZos[181]= - 1./2.*armZZos[10];
   armZZos[182]= - 1 + armZZos[181];
   armZZos[182]=armZZos[32]*armZZos[182];
   armZZos[145]=1./6.*armZZos[182] + 1./2.*armZZos[145] + armZZos[168];
   armZZos[145]=mmH*armZZos[145];
   armZZos[182]=2*armZZos[14];
   armZZos[189]=armZZos[182] - armZZos[16];
   armZZos[189]=3*armZZos[189] + armZZos[15];
   armZZos[189]=2*armZZos[189] + armZZos[247];
   armZZos[189]=armZZos[5]*armZZos[33]*armZZos[189];
   armZZos[190]= - 7*armZZos[16];
   armZZos[193]= - 1./3.*armZZos[15] - 65./2.*armZZos[14] + 
   armZZos[190];
   armZZos[193]=armZZos[32]*armZZos[193];
   armZZos[194]=3./8.*armZZos[32] - 1./2.*armZZos[6] - 15./2.*
   armZZos[9] + armZZos[10];
   armZZos[194]=armZZos[33]*armZZos[194];
   armZZos[123]=armZZos[123] + 1./2.*armZZos[145] + armZZos[255] + 
   armZZos[189] + armZZos[194] + 1./8.*armZZos[193] + 1./2.*
   armZZos[234] - 7./24.*armZZos[15] + 1./2.*armZZos[132] - 4./3.*
   armZZos[16];
   armZZos[123]=mmt*armZZos[123];
   armZZos[132]= - 67./8. + armZZos[11];
   armZZos[132]=armZZos[252] + 11./24.*armZZos[253] + 1./2.*
   armZZos[240] + armZZos[225] - 31./36.*armZZos[10] + 1./8.*
   armZZos[13] - 1./8.*armZZos[56] - 35./12.*armZZos[53] + 1./3.*
   armZZos[132] - 1./8.*armZZos[48];
   armZZos[132]=mmH*armZZos[132];
   armZZos[145]= - 1./4.*armZZos[18];
   armZZos[189]=31./18.*armZZos[15] + 373./36.*armZZos[16] + 
   armZZos[145] - 5./6.*armZZos[39] + armZZos[158] - 31./3.*armZZos[57]
    + 9./2.*armZZos[36];
   armZZos[189]=11./24.*armZZos[260] + 1./2.*armZZos[189] + 1./3.*
   armZZos[259];
   armZZos[193]=143./32. + armZZos[10];
   armZZos[122]=1./16.*armZZos[262] + 1./3.*armZZos[193] + 1./8.*
   armZZos[122];
   armZZos[122]=armZZos[33]*armZZos[122];
   armZZos[122]=1./8.*armZZos[132] + 1./8.*armZZos[189] + armZZos[122];
   armZZos[122]=mmH*armZZos[122];
   armZZos[132]= - 79./2.*armZZos[14] + 23./3.*armZZos[16];
   armZZos[132]= - 7./3.*armZZos[33] + 1./2.*armZZos[132] - 5./3.*
   armZZos[15];
   armZZos[132]=armZZos[33]*armZZos[132];
   armZZos[122]=armZZos[123] + 1./4.*armZZos[132] + armZZos[122];
   armZZos[122]=armZZos[112]*armZZos[122];
   armZZos[123]= - 1933./2.*armZZos[14] + armZZos[183];
   armZZos[123]=13*armZZos[33] + 1./2.*armZZos[123] + 77*armZZos[15];
   armZZos[123]=armZZos[33]*armZZos[123];
   armZZos[122]=armZZos[122] + armZZos[128] + 1./18.*armZZos[123] + 
   armZZos[136];
   armZZos[122]=armZZos[112]*armZZos[122];
   armZZos[123]= - 7./9. - armZZos[61];
   armZZos[128]=1./2.*armZZos[10];
   armZZos[132]=mmH*armZZos[42];
   armZZos[123]=armZZos[132] + armZZos[128] + armZZos[161] - 91./36.*
   armZZos[56] + 7./9.*armZZos[53] + armZZos[50] + 5./2.*armZZos[123]
    - armZZos[51];
   armZZos[123]=mmH*armZZos[123];
   armZZos[136]= - 469./2.*armZZos[59] + 115*armZZos[37];
   armZZos[136]= - 59./9.*armZZos[36] + 65./18.*armZZos[57] + 1./9.*
   armZZos[136] + 1./2.*armZZos[58];
   armZZos[183]= - 109 - 239./3.*armZZos[55];
   armZZos[183]=1./4.*armZZos[183] + 23./3.*armZZos[56];
   armZZos[183]=mmt*armZZos[183];
   armZZos[189]=5*armZZos[39];
   armZZos[193]=9./4.*armZZos[19];
   armZZos[194]= - 11./4.*armZZos[15];
   armZZos[123]=1./3.*armZZos[183] + armZZos[123] + 79./6.*armZZos[33]
    + armZZos[194] - 47./36.*armZZos[16] + 775./36.*armZZos[14] + 
   armZZos[145] + armZZos[193] - 33./4.*armZZos[17] + 1./2.*
   armZZos[136] + armZZos[189];
   armZZos[123]=mmt*armZZos[123];
   armZZos[136]= - 29./36.*armZZos[13] - 25./36.*armZZos[56] + 
   armZZos[222] - 7./36.*armZZos[48] - 107./36. + armZZos[61];
   armZZos[136]=35./432.*armZZos[253] + 1./4.*armZZos[136] - 17./27.*
   armZZos[10];
   armZZos[136]=mmH*armZZos[136];
   armZZos[183]=7./36.*armZZos[35] + 43./18.*armZZos[36] - armZZos[58]
    - 41./9.*armZZos[57];
   armZZos[198]=11 - armZZos[10];
   armZZos[198]=armZZos[33]*armZZos[198];
   armZZos[136]=armZZos[136] + 1./4.*armZZos[198] + 35./432.*
   armZZos[260] + 17./27.*armZZos[15] + 445./432.*armZZos[16] + 29./144.
   *armZZos[18] - 3./8.*armZZos[19] + 1./4.*armZZos[183] + 2./9.*
   armZZos[39];
   armZZos[136]=mmH*armZZos[136];
   armZZos[183]= - 7*armZZos[15];
   armZZos[198]=armZZos[183] + 33*armZZos[14] - armZZos[16];
   armZZos[201]= - 5*armZZos[33];
   armZZos[198]=1./2.*armZZos[198] + armZZos[201];
   armZZos[198]=armZZos[33]*armZZos[198];
   armZZos[123]=1./2.*armZZos[123] + 1./4.*armZZos[198] + armZZos[136];
   armZZos[123]=mmt*armZZos[123];
   armZZos[136]=3./2.*armZZos[14];
   armZZos[198]= - 1 - armZZos[55];
   armZZos[198]=mmt*armZZos[198];
   armZZos[198]=1./2.*armZZos[198] + armZZos[33] + armZZos[136] + 
   armZZos[159] - 1./2.*armZZos[17];
   armZZos[198]=mmt*armZZos[198];
   armZZos[205]=armZZos[33]*armZZos[14];
   armZZos[198]=1./2.*armZZos[205] + armZZos[198];
   armZZos[198]=armZZos[100]*mmt*armZZos[198];
   armZZos[205]= - 143./2.*armZZos[16] + 109*armZZos[15];
   armZZos[205]=1./3.*armZZos[205] - 25./2.*armZZos[33];
   armZZos[205]=armZZos[33]*armZZos[205];
   armZZos[206]=25./2. - 109./3.*armZZos[10];
   armZZos[206]=armZZos[33]*armZZos[206];
   armZZos[206]= - 25./2.*armZZos[39] + armZZos[206];
   armZZos[206]=mmH*armZZos[206];
   armZZos[205]=armZZos[205] + armZZos[206];
   armZZos[205]=mmH*armZZos[205];
   armZZos[123]=1./2.*armZZos[198] + 1./72.*armZZos[205] + armZZos[123]
   ;
   armZZos[123]=armZZos[100]*armZZos[123];
   armZZos[198]= - 1./3. - armZZos[61];
   armZZos[198]=armZZos[132] + armZZos[128] + armZZos[161] + 23./12.*
   armZZos[56] + armZZos[222] + armZZos[50] + 5./2.*armZZos[198] - 
   armZZos[51];
   armZZos[198]=mmH*armZZos[198];
   armZZos[205]= - 49*armZZos[59] + 23*armZZos[37];
   armZZos[205]= - 85./3.*armZZos[57] + 1./3.*armZZos[205] + 
   armZZos[58];
   armZZos[205]=1./2.*armZZos[205] + 7./3.*armZZos[36];
   armZZos[206]= - 27 - 13./3.*armZZos[55];
   armZZos[206]=1./2.*armZZos[206] - 19./3.*armZZos[56];
   armZZos[206]=mmt*armZZos[206];
   armZZos[198]=armZZos[206] + armZZos[198] + 22*armZZos[33] + 
   armZZos[194] + 91./12.*armZZos[16] + 40./3.*armZZos[14] + 
   armZZos[145] + armZZos[193] - 9*armZZos[17] + 1./2.*armZZos[205] + 
   armZZos[189];
   armZZos[198]=mmt*armZZos[198];
   armZZos[205]= - 23./12.*armZZos[13] + 5./12.*armZZos[56] + 
   armZZos[146] + 11./12.*armZZos[48] + 31./12. + armZZos[61];
   armZZos[206]=2./9.*armZZos[10];
   armZZos[205]=55./72.*armZZos[264] + 1./2.*armZZos[205] + 
   armZZos[206];
   armZZos[205]=mmH*armZZos[205];
   armZZos[207]= - 11./12.*armZZos[35] + 1./6.*armZZos[36] - 
   armZZos[58] + 13./3.*armZZos[57];
   armZZos[209]= - 9 - armZZos[10];
   armZZos[209]=armZZos[33]*armZZos[209];
   armZZos[205]=armZZos[205] + 1./2.*armZZos[209] + 55./72.*
   armZZos[261] - 2./9.*armZZos[15] - 185./72.*armZZos[16] + 23./24.*
   armZZos[18] - 3./4.*armZZos[19] + 1./2.*armZZos[207] - 2./3.*
   armZZos[39];
   armZZos[205]=mmH*armZZos[205];
   armZZos[207]= - 5./2.*armZZos[33] - 7./4.*armZZos[15] + 9*
   armZZos[14] - 1./4.*armZZos[16];
   armZZos[207]=armZZos[33]*armZZos[207];
   armZZos[198]=armZZos[198] + armZZos[207] + armZZos[205];
   armZZos[198]=mmt*armZZos[198];
   armZZos[205]=19./2.*armZZos[16];
   armZZos[207]= - 17*armZZos[15];
   armZZos[209]=armZZos[205] + armZZos[207];
   armZZos[209]=1./3.*armZZos[209] + 5./2.*armZZos[33];
   armZZos[209]=armZZos[33]*armZZos[209];
   armZZos[210]= - 5./2. + 17./3.*armZZos[10];
   armZZos[210]=armZZos[33]*armZZos[210];
   armZZos[210]=5./2.*armZZos[39] + armZZos[210];
   armZZos[210]=mmH*armZZos[210];
   armZZos[209]=armZZos[209] + armZZos[210];
   armZZos[209]=mmH*armZZos[209];
   armZZos[198]=1./12.*armZZos[209] + armZZos[198];
   armZZos[123]=armZZos[198] + armZZos[123];
   armZZos[123]=armZZos[100]*armZZos[123];
   armZZos[209]= - 3 - 5*armZZos[61];
   armZZos[128]=armZZos[132] + armZZos[128] + armZZos[161] - 3./4.*
   armZZos[56] - armZZos[53] + armZZos[50] + 1./2.*armZZos[209] - 
   armZZos[51];
   armZZos[128]=mmH*armZZos[128];
   armZZos[161]=9./4.*armZZos[16];
   armZZos[209]=11 + 35*armZZos[55];
   armZZos[209]=1./4.*armZZos[209] - armZZos[56];
   armZZos[209]=mmt*armZZos[209];
   armZZos[128]=armZZos[209] + armZZos[128] - 21./2.*armZZos[33] + 
   armZZos[194] + armZZos[161] - 31./4.*armZZos[14] + armZZos[145] + 
   armZZos[193] - 39./4.*armZZos[17] + armZZos[189] + armZZos[227] - 7./
   4.*armZZos[57] + 1./4.*armZZos[58] + 71./4.*armZZos[59] - 9*
   armZZos[37];
   armZZos[128]=mmt*armZZos[128];
   armZZos[189]=3./2.*armZZos[36];
   armZZos[139]=5./12.*armZZos[16] + armZZos[153] - 3./2.*armZZos[19]
    + armZZos[139] + armZZos[189] - armZZos[58] - armZZos[57];
   armZZos[153]= - 5./4.*armZZos[13] - 1./4.*armZZos[56] - armZZos[53]
    + 1./4.*armZZos[48] - 3./4. + armZZos[61];
   armZZos[193]= - 1./3.*armZZos[10];
   armZZos[153]=5./48.*armZZos[264] + 1./4.*armZZos[153] + armZZos[193]
   ;
   armZZos[153]=mmH*armZZos[153];
   armZZos[194]=3 - armZZos[10];
   armZZos[194]=armZZos[33]*armZZos[194];
   armZZos[139]=armZZos[153] + 1./4.*armZZos[194] + 5./48.*armZZos[261]
    + 1./4.*armZZos[139] + 1./3.*armZZos[15];
   armZZos[139]=mmH*armZZos[139];
   armZZos[153]=armZZos[183] + 39*armZZos[14] - armZZos[16];
   armZZos[153]=1./2.*armZZos[153] + armZZos[201];
   armZZos[153]=armZZos[33]*armZZos[153];
   armZZos[128]=1./2.*armZZos[128] + 1./4.*armZZos[153] + armZZos[139];
   armZZos[128]=mmt*armZZos[128];
   armZZos[139]= - 7./2.*armZZos[16];
   armZZos[153]=armZZos[139] + 5*armZZos[15];
   armZZos[153]=1./3.*armZZos[153] + armZZos[254];
   armZZos[153]=armZZos[33]*armZZos[153];
   armZZos[183]= - 5./3.*armZZos[10];
   armZZos[194]=1./2. + armZZos[183];
   armZZos[194]=armZZos[33]*armZZos[194];
   armZZos[194]= - 1./2.*armZZos[39] + armZZos[194];
   armZZos[194]=mmH*armZZos[194];
   armZZos[153]=armZZos[153] + armZZos[194];
   armZZos[153]=mmH*armZZos[153];
   armZZos[128]=1./8.*armZZos[153] + armZZos[128];
   armZZos[128]=armZZos[112]*armZZos[128];
   armZZos[128]=armZZos[198] + armZZos[128];
   armZZos[128]=armZZos[112]*armZZos[128];
   armZZos[123]=armZZos[123] + armZZos[128];
   armZZos[123]=armZZos[2]*armZZos[123];
   armZZos[128]= - armZZos[18] + armZZos[35] - 2*armZZos[39];
   armZZos[153]=8./27.*armZZos[15];
   armZZos[194]=armZZos[14] + 4./27.*armZZos[15];
   armZZos[194]=armZZos[32]*armZZos[194];
   armZZos[128]=8./27.*mmt - 8./27.*armZZos[33] + armZZos[194] + 
   armZZos[153] + 4./27.*armZZos[128] + armZZos[14];
   armZZos[128]=mmt*armZZos[128];
   armZZos[153]=4./27.*armZZos[33] + armZZos[14] + armZZos[153];
   armZZos[153]=armZZos[33]*armZZos[153];
   armZZos[128]=armZZos[153] + armZZos[128];
   armZZos[122]=armZZos[123] + armZZos[122] + 128./3.*armZZos[128] + 
   armZZos[143];
   armZZos[122]=armZZos[2]*armZZos[122];
   armZZos[123]= - 869 + 31*armZZos[9];
   armZZos[128]= - 1./2.*armZZos[87] + 11*armZZos[31];
   armZZos[128]=armZZos[14]*armZZos[128];
   armZZos[123]=1./12.*armZZos[123] + armZZos[128];
   armZZos[128]=13*armZZos[10];
   armZZos[143]=44./9.*armZZos[34] - 1./2. + 4*armZZos[1];
   armZZos[143]=armZZos[6]*armZZos[143];
   armZZos[153]= - 103./8.*armZZos[14] - 257./3.*armZZos[15];
   armZZos[153]=armZZos[88]*armZZos[153];
   armZZos[194]=173*armZZos[31] - 4925./6.*armZZos[88];
   armZZos[194]=armZZos[33]*armZZos[194];
   armZZos[198]=100./9.*armZZos[1];
   armZZos[201]=1100./81.*armZZos[34];
   armZZos[123]=1./6.*armZZos[194] - 4297./432.*armZZos[32] + 5./3.*
   armZZos[153] + 5./3.*armZZos[143] + 173./6.*armZZos[142] + 
   armZZos[201] + armZZos[198] + 1./2.*armZZos[123] + armZZos[128];
   armZZos[123]=armZZos[33]*armZZos[123];
   armZZos[143]=271./8.*armZZos[52] - 767./108. + 7*armZZos[49];
   armZZos[153]=22*armZZos[61];
   armZZos[194]= - 19./8.*armZZos[50];
   armZZos[209]= - 11./2.*armZZos[48];
   armZZos[210]=23*armZZos[53];
   armZZos[215]=61./6.*armZZos[56];
   armZZos[143]=armZZos[215] + armZZos[210] + armZZos[209] - 7./4.*
   armZZos[12] - 5771./288.*armZZos[38] + 181./192.*armZZos[55] + 
   armZZos[194] + 313./96.*armZZos[54] + 809./54.*armZZos[51] + 1./4.*
   armZZos[143] + armZZos[153];
   armZZos[216]=7*armZZos[14];
   armZZos[217]=13*armZZos[16];
   armZZos[220]=armZZos[216] + armZZos[217];
   armZZos[222]= - 11*armZZos[15];
   armZZos[220]=1./3.*armZZos[220] + armZZos[222];
   armZZos[220]=armZZos[32]*armZZos[220];
   armZZos[223]= - 6*armZZos[10];
   armZZos[225]=armZZos[223] + 128./3. - 9*armZZos[9];
   armZZos[225]=armZZos[33]*armZZos[225];
   armZZos[226]=11*armZZos[18] + 13./6.*armZZos[19] + 3*armZZos[17] + 
   25*armZZos[39] - 3*armZZos[37] - 11*armZZos[35];
   armZZos[227]= - 13./6.*armZZos[16];
   armZZos[228]= - armZZos[14] - armZZos[15];
   armZZos[229]=36*armZZos[5]*armZZos[33]*armZZos[228];
   armZZos[220]=armZZos[229] + armZZos[225] + 1./2.*armZZos[220] + 
   armZZos[133] + armZZos[227] + armZZos[226] + 62./3.*armZZos[14];
   armZZos[220]=armZZos[5]*armZZos[220];
   armZZos[225]= - 1271./3. - 173*armZZos[9];
   armZZos[230]=11*armZZos[10];
   armZZos[225]=1./12.*armZZos[225] + armZZos[230];
   armZZos[231]=110./9.*armZZos[34] - 1 + 10*armZZos[1];
   armZZos[231]=armZZos[6]*armZZos[231];
   armZZos[198]= - 61./27.*armZZos[32] + 2./3.*armZZos[231] + 
   armZZos[201] + 1./2.*armZZos[225] + armZZos[198];
   armZZos[198]=armZZos[32]*armZZos[198];
   armZZos[201]=2*armZZos[47] + armZZos[45];
   armZZos[201]=56*armZZos[201] - 997./4.*armZZos[44];
   armZZos[225]=44*armZZos[42];
   armZZos[201]=armZZos[225] + armZZos[46] + armZZos[43] + 1./27.*
   armZZos[201] + 19./2.*armZZos[40];
   armZZos[231]=armZZos[103] + 11./3.*armZZos[56] + 20./3. + 
   armZZos[219];
   armZZos[231]=2*armZZos[89]*armZZos[231];
   armZZos[228]=armZZos[5]*armZZos[228];
   armZZos[228]=18*armZZos[228] - 35 + 6*armZZos[10];
   armZZos[228]=armZZos[5]*armZZos[228];
   armZZos[201]=armZZos[231] + 1./3.*armZZos[201] + armZZos[228];
   armZZos[201]=mmt*armZZos[201];
   armZZos[232]= - 3*armZZos[58];
   armZZos[233]=armZZos[232] + 11./3.*armZZos[57];
   armZZos[233]= - 8./3.*armZZos[33] + 11./3.*armZZos[184] + 5./3.*
   armZZos[15] - 17./2.*armZZos[16] + 29./3.*armZZos[18] - 29./6.*
   armZZos[19] + 6*armZZos[39] - 11./3.*armZZos[35] + 2*armZZos[233] + 
   3*armZZos[36];
   armZZos[233]=armZZos[89]*armZZos[233];
   armZZos[234]=220./27.*armZZos[34] + 7./4. + 20./3.*armZZos[1];
   armZZos[234]=armZZos[6]*armZZos[234];
   armZZos[235]= - 11./2.*armZZos[13];
   armZZos[236]= - 10./3.*mmH*armZZos[42];
   armZZos[143]=armZZos[201] + armZZos[236] + armZZos[233] + 
   armZZos[220] + 1./3.*armZZos[198] + 1./3.*armZZos[234] + 1100./243.*
   armZZos[34] + 100./27.*armZZos[1] - 293./72.*armZZos[9] + 1./3.*
   armZZos[143] + armZZos[235];
   armZZos[143]=mmt*armZZos[143];
   armZZos[198]= - 11./12.*armZZos[16];
   armZZos[201]=3*armZZos[39];
   armZZos[130]=68./9.*armZZos[33] + 7./18.*armZZos[130] + 95./18.*
   armZZos[15] + armZZos[198] + 47./18.*armZZos[18] - 47./36.*
   armZZos[19] + armZZos[201] + 7./18.*armZZos[35] + armZZos[189] + 
   armZZos[232] - 7./9.*armZZos[57];
   armZZos[130]=armZZos[89]*armZZos[130];
   armZZos[220]= - 5./3. - armZZos[9];
   armZZos[234]= - 3*armZZos[1];
   armZZos[237]=armZZos[234] - 11./3.*armZZos[34];
   armZZos[237]=armZZos[6]*armZZos[237];
   armZZos[220]=2*armZZos[237] + 88./9.*armZZos[34] + 8*armZZos[1] + 1./
   2.*armZZos[220] + armZZos[94];
   armZZos[220]=armZZos[33]*armZZos[220];
   armZZos[111]=armZZos[111] - armZZos[39];
   armZZos[111]=41./12.*armZZos[16] - 7./2.*armZZos[18] + 7*
   armZZos[111] - 41./12.*armZZos[19];
   armZZos[237]= - 41./3.*armZZos[16] + 7*armZZos[15];
   armZZos[237]=armZZos[32]*armZZos[237];
   armZZos[238]= - armZZos[5]*armZZos[33]*armZZos[15];
   armZZos[111]=18*armZZos[238] + armZZos[220] + 1./12.*armZZos[237] + 
   1./3.*armZZos[111] + 13*armZZos[15];
   armZZos[111]=armZZos[5]*armZZos[111];
   armZZos[220]= - armZZos[5]*armZZos[15];
   armZZos[220]=9*armZZos[220] - 14./3. + armZZos[107];
   armZZos[220]=armZZos[5]*armZZos[220];
   armZZos[237]=armZZos[103] - 7./9.*armZZos[56] + 20./9. + 
   armZZos[219];
   armZZos[237]=armZZos[89]*armZZos[237];
   armZZos[152]=armZZos[152] - 17./3.*armZZos[43] - 437./108.*
   armZZos[44] - armZZos[41];
   armZZos[152]=1./2.*armZZos[152] - 14./3.*armZZos[42];
   armZZos[152]=armZZos[237] + 1./3.*armZZos[152] + armZZos[220];
   armZZos[152]=mmt*armZZos[152];
   armZZos[220]= - 137./4.*armZZos[9];
   armZZos[237]= - 751./3. + armZZos[220];
   armZZos[237]= - 539./27.*armZZos[34] - 49./3.*armZZos[1] + 1./3.*
   armZZos[237] - 7*armZZos[10];
   armZZos[238]=77./36.*armZZos[34] - 1 + 7./4.*armZZos[1];
   armZZos[238]=armZZos[6]*armZZos[238];
   armZZos[237]= - 157./18.*armZZos[32] + 1./4.*armZZos[237] + 
   armZZos[238];
   armZZos[237]=armZZos[32]*armZZos[237];
   armZZos[238]=331./144.*armZZos[52] - 47./648. + armZZos[49];
   armZZos[239]= - 1./4.*armZZos[12];
   armZZos[240]=671./27.*armZZos[34] + 61./3.*armZZos[1] + 5./2. + 
   armZZos[9];
   armZZos[240]=armZZos[6]*armZZos[240];
   armZZos[111]=armZZos[152] + 25./9.*armZZos[132] + armZZos[130] + 
   armZZos[111] + 1./9.*armZZos[237] + 1./12.*armZZos[240] - 2915./972.
   *armZZos[34] - 265./108.*armZZos[1] - 47./432.*armZZos[9] + 7./12.*
   armZZos[13] - 227./108.*armZZos[56] - 31./18.*armZZos[53] + 7./36.*
   armZZos[48] + armZZos[239] - 5071./5184.*armZZos[38] - 257./1152.*
   armZZos[55] + 143./144.*armZZos[50] - 191./324.*armZZos[51] + 1./4.*
   armZZos[238] - 7./9.*armZZos[61];
   armZZos[111]=mmt*armZZos[111];
   armZZos[130]=1127 + 61./3.*armZZos[9];
   armZZos[152]=671./9.*armZZos[34] - 1 + 61*armZZos[1];
   armZZos[152]=armZZos[6]*armZZos[152];
   armZZos[237]=1./8.*armZZos[14] + 19./3.*armZZos[15];
   armZZos[237]=armZZos[88]*armZZos[237];
   armZZos[163]=armZZos[163] + 3875./6.*armZZos[88];
   armZZos[163]=armZZos[33]*armZZos[163];
   armZZos[130]=1./2.*armZZos[163] + 263./144.*armZZos[32] + 25*
   armZZos[237] + 1./2.*armZZos[152] + 107./2.*armZZos[167] - 2915./54.
   *armZZos[34] - 265./6.*armZZos[1] + 1./8.*armZZos[130] - 41*
   armZZos[10];
   armZZos[130]=armZZos[33]*armZZos[130];
   armZZos[152]= - 1741./72.*armZZos[39] - 1./3.*armZZos[35] + 25./4.*
   armZZos[36] - 259./18.*armZZos[57] + 215./27.*armZZos[8] + 25./8.*
   armZZos[159] - 59./9.*armZZos[58];
   armZZos[159]=185./216.*armZZos[15] + 4./9.*armZZos[14] + 5./8.*
   armZZos[16];
   armZZos[159]=armZZos[6]*armZZos[159];
   armZZos[163]= - 37./18.*armZZos[15] + 1483./144.*armZZos[14] + 
   armZZos[248];
   armZZos[163]=armZZos[32]*armZZos[163];
   armZZos[130]=1./6.*armZZos[130] + 1./4.*armZZos[163] + armZZos[159]
    - 3157./432.*armZZos[15] - 67./18.*armZZos[16] - 2095./864.*
   armZZos[14] + 311./216.*armZZos[18] + 73./24.*armZZos[19] + 1./8.*
   armZZos[152] + 5./3.*armZZos[17];
   armZZos[152]=armZZos[57] + armZZos[243];
   armZZos[159]= - 1./6.*armZZos[15] + armZZos[198] - 1./12.*
   armZZos[18] + 13./6.*armZZos[39] + armZZos[152] + 1./12.*armZZos[35]
   ;
   armZZos[159]=armZZos[88]*armZZos[159];
   armZZos[163]=7./2. + 109*armZZos[10];
   armZZos[198]=armZZos[88]*armZZos[258];
   armZZos[163]=1./3.*armZZos[163] + 25./4.*armZZos[198];
   armZZos[163]=armZZos[32]*armZZos[163];
   armZZos[208]= - 5*armZZos[25] + 1811./72. + armZZos[208];
   armZZos[237]=armZZos[6]*armZZos[141];
   armZZos[159]=1./18.*armZZos[163] + 25./6.*armZZos[159] + 5./27.*
   armZZos[237] + 347./81.*armZZos[10] + 19./18.*armZZos[13] + 125./36.
   *armZZos[56] + 523./108.*armZZos[53] + armZZos[101] - 37./8.*
   armZZos[38] - 5./54.*armZZos[11] + 1./9.*armZZos[208] + 37./8.*
   armZZos[50];
   armZZos[163]=1./6.*armZZos[13];
   armZZos[208]=7 - armZZos[48];
   armZZos[188]=armZZos[188] + armZZos[163] - 1./6.*armZZos[56] + 1./6.
   *armZZos[208] + armZZos[53];
   armZZos[188]=armZZos[88]*armZZos[188];
   armZZos[238]= - armZZos[32]*armZZos[88]*armZZos[10];
   armZZos[188]=armZZos[188] + 1./6.*armZZos[238];
   armZZos[188]=mmH*armZZos[188];
   armZZos[240]= - 23./16. + armZZos[10];
   armZZos[240]=armZZos[33]*armZZos[88]*armZZos[240];
   armZZos[159]=25./48.*armZZos[188] + 1./4.*armZZos[159] + 25./9.*
   armZZos[240];
   armZZos[159]=mmH*armZZos[159];
   armZZos[188]=armZZos[15] - 1./3.*armZZos[33];
   armZZos[188]=armZZos[5]*armZZos[33]*armZZos[188];
   armZZos[242]= - armZZos[16] + 2*armZZos[15];
   armZZos[242]=armZZos[89]*armZZos[33]*armZZos[242];
   armZZos[111]=armZZos[111] + 1./2.*armZZos[159] + 17./9.*armZZos[242]
    + 1./3.*armZZos[130] + 17./2.*armZZos[188];
   armZZos[111]=armZZos[100]*armZZos[111];
   armZZos[130]= - 11./2. - 17*armZZos[10];
   armZZos[159]=armZZos[88]*armZZos[221];
   armZZos[130]=1./3.*armZZos[130] + 5./4.*armZZos[159];
   armZZos[130]=armZZos[32]*armZZos[130];
   armZZos[159]= - armZZos[25] + 115./24. + armZZos[27];
   armZZos[188]=1./6.*armZZos[15] + 11./12.*armZZos[16] + 1./12.*
   armZZos[18] - 13./6.*armZZos[39] - 1./12.*armZZos[35] - armZZos[57]
    + armZZos[214];
   armZZos[188]=armZZos[88]*armZZos[188];
   armZZos[130]=1./6.*armZZos[130] + 5./2.*armZZos[188] + 1./9.*
   armZZos[237] + armZZos[168] + 5./6.*armZZos[13] - 25./12.*
   armZZos[56] + 1./36.*armZZos[53] + armZZos[101] + 3./8.*armZZos[38]
    - 1./18.*armZZos[11] + 1./3.*armZZos[159] - 3./8.*armZZos[50];
   armZZos[159]= - 7 + armZZos[48];
   armZZos[159]=armZZos[193] - 1./6.*armZZos[13] + 1./6.*armZZos[56] + 
   1./6.*armZZos[159] - armZZos[53];
   armZZos[159]=armZZos[88]*armZZos[159];
   armZZos[168]=armZZos[32]*armZZos[88]*armZZos[10];
   armZZos[159]=armZZos[159] + 1./6.*armZZos[168];
   armZZos[159]=mmH*armZZos[159];
   armZZos[168]=23./16. - armZZos[10];
   armZZos[168]=armZZos[33]*armZZos[88]*armZZos[168];
   armZZos[130]=5./16.*armZZos[159] + 1./4.*armZZos[130] + 5./3.*
   armZZos[168];
   armZZos[130]=mmH*armZZos[130];
   armZZos[159]= - 3*armZZos[15];
   armZZos[168]= - 15*armZZos[33];
   armZZos[188]=armZZos[168] + 35./3.*armZZos[14] + armZZos[159];
   armZZos[188]=armZZos[5]*armZZos[33]*armZZos[188];
   armZZos[193]=47./8.*armZZos[57] - 83./36.*armZZos[8] - 139./12.*
   armZZos[58] + 425./64.*armZZos[37] + 68./3.*armZZos[60] + 151./32.*
   armZZos[59];
   armZZos[214]=1./4.*armZZos[16];
   armZZos[221]= - 101./324.*armZZos[15] + 46./27.*armZZos[14] + 
   armZZos[214];
   armZZos[221]=armZZos[6]*armZZos[221];
   armZZos[243]= - 19./18.*armZZos[15];
   armZZos[244]=armZZos[243] + 1151./432.*armZZos[14] - armZZos[16];
   armZZos[244]=armZZos[32]*armZZos[244];
   armZZos[247]= - 5./16.*armZZos[36];
   armZZos[248]= - 7./12.*armZZos[35];
   armZZos[249]= - 5./12.*armZZos[16];
   armZZos[156]=armZZos[16] + armZZos[156];
   armZZos[156]=2./3.*armZZos[89]*armZZos[33]*armZZos[156];
   armZZos[250]= - 1./4.*armZZos[19];
   armZZos[111]=armZZos[111] + armZZos[143] + armZZos[130] + 
   armZZos[156] + armZZos[188] + 1./3.*armZZos[123] + 1./2.*
   armZZos[244] + armZZos[221] - 1477./648.*armZZos[15] + armZZos[249]
    - 3169./648.*armZZos[14] + 469./324.*armZZos[18] + armZZos[250] - 
   armZZos[17] + 2099./2592.*armZZos[39] + armZZos[248] + 1./9.*
   armZZos[193] + armZZos[247];
   armZZos[111]=armZZos[100]*armZZos[111];
   armZZos[123]= - 113./8.*armZZos[52] + 23./36. + 43*armZZos[49];
   armZZos[123]=armZZos[215] + armZZos[210] + armZZos[209] - 43./4.*
   armZZos[12] - 89./96.*armZZos[38] + 5671./64.*armZZos[55] + 
   armZZos[194] + 19./32.*armZZos[54] + 11./2.*armZZos[51] + 1./4.*
   armZZos[123] + armZZos[153];
   armZZos[143]=armZZos[222] - 19*armZZos[14] + 13./3.*armZZos[16];
   armZZos[143]=armZZos[32]*armZZos[143];
   armZZos[153]=armZZos[223] + 56./3. - 33*armZZos[9];
   armZZos[153]=armZZos[33]*armZZos[153];
   armZZos[133]=armZZos[229] + armZZos[153] + 1./2.*armZZos[143] + 
   armZZos[133] + armZZos[227] + armZZos[226] + 10*armZZos[14];
   armZZos[133]=armZZos[5]*armZZos[133];
   armZZos[143]= - 2*armZZos[47] - armZZos[45];
   armZZos[153]=armZZos[225] + armZZos[46] + armZZos[43] + 35./2.*
   armZZos[40] + 8./3.*armZZos[143] + 1./4.*armZZos[44];
   armZZos[153]=armZZos[231] + 1./3.*armZZos[153] + armZZos[228];
   armZZos[153]=mmt*armZZos[153];
   armZZos[188]=4./3.*armZZos[1];
   armZZos[193]=4*armZZos[34];
   armZZos[194]=armZZos[193] + armZZos[188] + 7./4. + 2*armZZos[9];
   armZZos[194]=armZZos[6]*armZZos[194];
   armZZos[209]= - 1607./9. - 799*armZZos[9];
   armZZos[209]=1./4.*armZZos[209] + armZZos[230];
   armZZos[210]=2./3.*armZZos[1];
   armZZos[215]=2*armZZos[34] - 1 + armZZos[210];
   armZZos[215]=armZZos[6]*armZZos[215];
   armZZos[209]=11./3.*armZZos[32] + 2*armZZos[215] + 20./3.*
   armZZos[34] + 1./2.*armZZos[209] + 20./9.*armZZos[1];
   armZZos[209]=armZZos[32]*armZZos[209];
   armZZos[215]=20./27.*armZZos[1];
   armZZos[221]=20./9.*armZZos[34];
   armZZos[123]=armZZos[153] + armZZos[236] + armZZos[233] + 
   armZZos[133] + 1./3.*armZZos[209] + 1./3.*armZZos[194] + 
   armZZos[221] + armZZos[215] - 895./24.*armZZos[9] + 1./3.*
   armZZos[123] + armZZos[235];
   armZZos[123]=mmt*armZZos[123];
   armZZos[133]= - 23./3. + armZZos[148];
   armZZos[153]=8./3.*armZZos[1];
   armZZos[96]=2*armZZos[96] + 8*armZZos[34] + armZZos[153] + 1./2.*
   armZZos[133] + armZZos[94];
   armZZos[96]=armZZos[33]*armZZos[96];
   armZZos[133]=armZZos[17] + armZZos[192] - armZZos[37] + armZZos[173]
   ;
   armZZos[192]= - 3*armZZos[14];
   armZZos[194]= - 3./2.*armZZos[15] + armZZos[192] + armZZos[120];
   armZZos[194]=armZZos[32]*armZZos[194];
   armZZos[209]=armZZos[245] - armZZos[15];
   armZZos[222]=armZZos[5]*armZZos[33]*armZZos[209];
   armZZos[96]=18*armZZos[222] + armZZos[96] + 1./2.*armZZos[194] + 9*
   armZZos[15] + armZZos[214] + 18*armZZos[14] + 3./2.*armZZos[18] + 3*
   armZZos[133] + armZZos[250];
   armZZos[96]=armZZos[5]*armZZos[96];
   armZZos[133]=4*armZZos[33] + 1./2.*armZZos[184] + 7./2.*armZZos[15]
    - 9./4.*armZZos[16] + 7./2.*armZZos[18] - 7./4.*armZZos[19] + 
   armZZos[201] + armZZos[173] + armZZos[189] + armZZos[232] + 
   armZZos[57];
   armZZos[133]=armZZos[89]*armZZos[133];
   armZZos[103]=armZZos[103] + armZZos[56] + 4 + armZZos[219];
   armZZos[103]=armZZos[89]*armZZos[103];
   armZZos[173]=3*armZZos[41];
   armZZos[184]=3*armZZos[46] - armZZos[43] + 3*armZZos[40] - 5./4.*
   armZZos[44] + armZZos[173];
   armZZos[189]=armZZos[5]*armZZos[209];
   armZZos[189]=3*armZZos[189] - 3 + armZZos[10];
   armZZos[189]=armZZos[5]*armZZos[189];
   armZZos[103]=armZZos[103] + 3*armZZos[189] + 1./2.*armZZos[184] + 
   armZZos[155];
   armZZos[103]=mmt*armZZos[103];
   armZZos[155]=53./9. - 35./2.*armZZos[9];
   armZZos[155]=7./3.*armZZos[34] + 7./9.*armZZos[1] + 1./2.*
   armZZos[155] + armZZos[10];
   armZZos[184]= - 1./4.*armZZos[34] + 1 - 1./12.*armZZos[1];
   armZZos[184]=armZZos[6]*armZZos[184];
   armZZos[155]= - 1./2.*armZZos[32] + 1./4.*armZZos[155] + 
   armZZos[184];
   armZZos[155]=armZZos[32]*armZZos[155];
   armZZos[184]=5*armZZos[34] + 5./3.*armZZos[1] - 1./2. - 7*armZZos[9]
   ;
   armZZos[184]=armZZos[6]*armZZos[184];
   armZZos[96]=armZZos[103] + armZZos[132] + armZZos[133] + armZZos[96]
    + armZZos[155] + 1./4.*armZZos[184] - 17./12.*armZZos[34] - 17./36.
   *armZZos[1] - 57./16.*armZZos[9] - 3./4.*armZZos[13] - 7./12.*
   armZZos[56] + armZZos[187] - 1./4.*armZZos[48] + 12*armZZos[12] - 3./
   64.*armZZos[38] - 3533./128.*armZZos[55] + 7./16.*armZZos[50] + 25./
   32.*armZZos[54] + 1./4.*armZZos[51] + armZZos[61] - 41./64.*
   armZZos[52] + 149./288. - 12*armZZos[49];
   armZZos[96]=mmt*armZZos[96];
   armZZos[103]=6331./9. - 453*armZZos[9];
   armZZos[132]= - 3./2.*armZZos[87] + armZZos[31];
   armZZos[132]=armZZos[14]*armZZos[132];
   armZZos[133]=armZZos[34] - 1 + armZZos[93];
   armZZos[133]=armZZos[6]*armZZos[133];
   armZZos[155]=9./8.*armZZos[14] + armZZos[15];
   armZZos[155]=armZZos[88]*armZZos[155];
   armZZos[184]= - 3*armZZos[31] + 25./2.*armZZos[88];
   armZZos[184]=armZZos[33]*armZZos[184];
   armZZos[103]=1./2.*armZZos[184] - 9./16.*armZZos[32] + armZZos[155]
    + 5./2.*armZZos[133] + 3./2.*armZZos[167] - 17./6.*armZZos[34] - 17.
   /18.*armZZos[1] - armZZos[10] + 1./8.*armZZos[103] + armZZos[132];
   armZZos[103]=armZZos[33]*armZZos[103];
   armZZos[132]= - 1./2.*armZZos[15] - 11./4.*armZZos[16] + 
   armZZos[145] + 13./2.*armZZos[39] + 3*armZZos[152] + armZZos[158];
   armZZos[132]=armZZos[88]*armZZos[132];
   armZZos[133]= - 1./2. + armZZos[138];
   armZZos[133]=1./3.*armZZos[133] + 1./4.*armZZos[198];
   armZZos[133]=armZZos[32]*armZZos[133];
   armZZos[101]=1./2.*armZZos[133] + 1./2.*armZZos[132] + 1./3.*
   armZZos[237] + 43./9.*armZZos[10] + 3./2.*armZZos[13] + 5./4.*
   armZZos[56] + 35./12.*armZZos[53] + armZZos[101] - 21./8.*
   armZZos[38] - 1./6.*armZZos[11] + 21./8.*armZZos[50] - armZZos[25]
    + 235./72. + armZZos[27];
   armZZos[117]=armZZos[10] + 1./2.*armZZos[13] + armZZos[117] + 1./2.*
   armZZos[208] + armZZos[146];
   armZZos[117]=armZZos[88]*armZZos[117];
   armZZos[117]=armZZos[117] + 1./2.*armZZos[238];
   armZZos[117]=mmH*armZZos[117];
   armZZos[101]=1./16.*armZZos[117] + 1./4.*armZZos[101] + armZZos[240]
   ;
   armZZos[101]=mmH*armZZos[101];
   armZZos[117]=armZZos[14] + 1./2.*armZZos[15];
   armZZos[132]=3*armZZos[117] + armZZos[254];
   armZZos[132]=armZZos[5]*armZZos[33]*armZZos[132];
   armZZos[133]= - 943*armZZos[59] + 399./2.*armZZos[37];
   armZZos[133]= - 1011./8.*armZZos[39] - armZZos[35] + 3./4.*
   armZZos[36] - 11./6.*armZZos[57] + 11*armZZos[8] + 1./8.*
   armZZos[133] - armZZos[58];
   armZZos[138]=armZZos[14] + armZZos[150];
   armZZos[138]=3*armZZos[138] + 31./24.*armZZos[15];
   armZZos[138]=armZZos[6]*armZZos[138];
   armZZos[146]= - 5./6.*armZZos[15] + 121./16.*armZZos[14] + 3*
   armZZos[16];
   armZZos[146]=armZZos[32]*armZZos[146];
   armZZos[96]=armZZos[96] + 1./2.*armZZos[101] + armZZos[242] + 3*
   armZZos[132] + 1./2.*armZZos[103] + 1./4.*armZZos[146] + 
   armZZos[138] - 685./144.*armZZos[15] - 13./12.*armZZos[16] + 1321./
   96.*armZZos[14] + 11./8.*armZZos[18] + 5./8.*armZZos[19] + 1./8.*
   armZZos[133] + 10*armZZos[17];
   armZZos[96]=armZZos[112]*armZZos[96];
   armZZos[101]= - 16559./9. - 499*armZZos[9];
   armZZos[103]=1./2.*armZZos[87] - 5./3.*armZZos[31];
   armZZos[103]=armZZos[14]*armZZos[103];
   armZZos[101]=1./12.*armZZos[101] + armZZos[103];
   armZZos[103]=armZZos[193] + 1./2. + armZZos[188];
   armZZos[103]=armZZos[6]*armZZos[103];
   armZZos[132]= - 43./8.*armZZos[14] - 29./3.*armZZos[15];
   armZZos[132]=armZZos[88]*armZZos[132];
   armZZos[133]=armZZos[180] - 149./6.*armZZos[88];
   armZZos[133]=armZZos[33]*armZZos[133];
   armZZos[101]=1./2.*armZZos[133] - 193./144.*armZZos[32] + 
   armZZos[132] + 1./3.*armZZos[103] + 5./2.*armZZos[142] + 
   armZZos[221] + armZZos[215] + 1./2.*armZZos[101] + 13./3.*
   armZZos[10];
   armZZos[101]=armZZos[33]*armZZos[101];
   armZZos[103]=armZZos[168] + armZZos[14] + armZZos[159];
   armZZos[103]=armZZos[5]*armZZos[33]*armZZos[103];
   armZZos[132]= - 29./36.*armZZos[15] - 14./9.*armZZos[14] + 
   armZZos[214];
   armZZos[132]=armZZos[6]*armZZos[132];
   armZZos[133]=armZZos[243] - 2987./144.*armZZos[14] - armZZos[16];
   armZZos[133]=armZZos[32]*armZZos[133];
   armZZos[96]=armZZos[96] + armZZos[123] + armZZos[130] + armZZos[156]
    + armZZos[103] + armZZos[101] + 1./2.*armZZos[133] + armZZos[132]
    + 11./8.*armZZos[15] + armZZos[249] - 2449./72.*armZZos[14] - 17./
   12.*armZZos[18] + armZZos[250] - 26./3.*armZZos[17] + 2107./288.*
   armZZos[39] + armZZos[248] + armZZos[247] + 47./72.*armZZos[57] - 3./
   4.*armZZos[8] + 13./12.*armZZos[58] - 1949./192.*armZZos[37] + 
   armZZos[140] + 991./32.*armZZos[59];
   armZZos[96]=armZZos[112]*armZZos[96];
   armZZos[101]=armZZos[14] + 32./9.*armZZos[15];
   armZZos[101]=armZZos[88]*armZZos[101];
   armZZos[103]= - armZZos[31] + 8./3.*armZZos[88];
   armZZos[103]=armZZos[33]*armZZos[103];
   armZZos[123]=pow(c,2);
   armZZos[130]=107./9. + 8*armZZos[123];
   armZZos[130]=armZZos[9]*armZZos[130];
   armZZos[132]= - armZZos[1] - 5./3.*armZZos[34];
   armZZos[133]=armZZos[6]*armZZos[132];
   armZZos[101]=32./3.*armZZos[103] + 16*armZZos[101] + 16./3.*
   armZZos[133] + 32./3.*armZZos[167] - 80./27.*armZZos[34] - 16./9.*
   armZZos[1] + 4*armZZos[130] + 239./9. + 32*armZZos[123];
   armZZos[101]=armZZos[33]*armZZos[101];
   armZZos[103]= - 32./3.*armZZos[60] - 29*armZZos[59];
   armZZos[103]= - 172./9.*armZZos[18] - 128./9.*armZZos[39] + 20./9.*
   armZZos[8] + 64./3.*armZZos[58] + 2*armZZos[103] + 29*armZZos[37];
   armZZos[138]=6*armZZos[14];
   armZZos[140]=armZZos[14] + 5./27.*armZZos[15];
   armZZos[140]=armZZos[6]*armZZos[140];
   armZZos[146]=armZZos[32]*armZZos[14];
   armZZos[101]=2./3.*armZZos[101] + 8*armZZos[146] + 4./3.*
   armZZos[140] + 140./81.*armZZos[15] + 1./9.*armZZos[103] + 
   armZZos[138];
   armZZos[103]=4*armZZos[123];
   armZZos[140]=137./27. + armZZos[103];
   armZZos[146]= - 20./27.*armZZos[34];
   armZZos[150]=4./3.*armZZos[133];
   armZZos[140]=8./9.*armZZos[32] + armZZos[150] + armZZos[146] + 
   armZZos[154] + 2*armZZos[140] + armZZos[130];
   armZZos[140]=armZZos[32]*armZZos[140];
   armZZos[152]=256./9.*armZZos[38] - 58*armZZos[55] + 29*armZZos[54]
    + 131 - 512./9.*armZZos[51];
   armZZos[143]=armZZos[143] + armZZos[44];
   armZZos[143]=mmt*armZZos[143];
   armZZos[130]=512./27.*armZZos[143] + 16*armZZos[140] + 64./3.*
   armZZos[133] - 320./27.*armZZos[34] - 64./9.*armZZos[1] + 16*
   armZZos[130] + 1./3.*armZZos[152] + 128*armZZos[123];
   armZZos[130]=mmt*armZZos[130];
   armZZos[96]=armZZos[122] + armZZos[96] + armZZos[111] + 2*
   armZZos[101] + 1./3.*armZZos[130];
   armZZos[96]=armZZos[2]*armZZos[96];
   armZZos[101]=2*armZZos[123];
   armZZos[111]=13./3. + armZZos[101];
   armZZos[111]=2*armZZos[9]*armZZos[111];
   armZZos[122]=4./9.*armZZos[1];
   armZZos[130]=20./27.*armZZos[34];
   armZZos[133]=armZZos[150] + armZZos[130] + armZZos[122] + 
   armZZos[111] + 127./9. + armZZos[103];
   armZZos[133]=armZZos[88]*armZZos[133];
   armZZos[140]= - armZZos[32]*armZZos[88];
   armZZos[133]=armZZos[133] + 16./9.*armZZos[140];
   armZZos[133]=armZZos[32]*armZZos[133];
   armZZos[143]=5*armZZos[54] - 64./3.*armZZos[51] - 209./9. - 4*
   armZZos[52];
   armZZos[152]=armZZos[123]*armZZos[54];
   armZZos[155]= - 2*armZZos[123];
   armZZos[156]= - 13./3. + armZZos[155];
   armZZos[156]=armZZos[9]*armZZos[156];
   armZZos[158]=armZZos[1] + 5./3.*armZZos[34];
   armZZos[158]=armZZos[6]*armZZos[158];
   armZZos[143]=8./9.*armZZos[158] - 40./81.*armZZos[34] - 8./27.*
   armZZos[1] + 4./3.*armZZos[156] + 8./3.*armZZos[152] + 196./27.*
   armZZos[38] + 1./3.*armZZos[143] + 2*armZZos[55];
   armZZos[143]=armZZos[88]*armZZos[143];
   armZZos[152]=armZZos[123]*armZZos[110];
   armZZos[110]=armZZos[110] + 1./3.*armZZos[152];
   armZZos[110]=mmZ*armZZos[110];
   armZZos[152]= - 64*armZZos[44] + 16*armZZos[45] + armZZos[24];
   armZZos[152]=17./9.*armZZos[46] - 16./81.*armZZos[22] + 5*
   armZZos[40] + 4./81.*armZZos[152] + armZZos[173];
   armZZos[156]=2./3.*armZZos[46] + armZZos[41] + armZZos[191];
   armZZos[159]= - armZZos[55] - 4 - armZZos[54];
   armZZos[159]=armZZos[87]*armZZos[159];
   armZZos[156]=4*armZZos[156] + armZZos[159];
   armZZos[156]=armZZos[123]*armZZos[156];
   armZZos[159]= - 1 - armZZos[54];
   armZZos[159]=armZZos[123]*armZZos[159];
   armZZos[159]=4*armZZos[159] - 100./9.*armZZos[38] + 2*armZZos[54] + 
   64./9.*armZZos[51] - 19./3. + 4*armZZos[52];
   armZZos[159]=armZZos[31]*armZZos[159];
   armZZos[167]= - 8./3.*armZZos[55] - 43./4. - 10./3.*armZZos[54];
   armZZos[167]=armZZos[87]*armZZos[167];
   armZZos[168]=4*armZZos[32] - 20./3. + armZZos[6];
   armZZos[168]=armZZos[5]*armZZos[168];
   armZZos[109]=armZZos[110] + 8./9.*armZZos[168] + 1./3.*armZZos[109]
    + 4./3.*armZZos[133] + 2*armZZos[143] + 1./3.*armZZos[118] + 1./3.*
   armZZos[116] + 2./3.*armZZos[159] + 1./3.*armZZos[156] + 1./3.*
   armZZos[114] + 1./3.*armZZos[113] + 2*armZZos[152] + 1./3.*
   armZZos[167];
   armZZos[109]=mmZ*armZZos[109];
   armZZos[101]=77./27. + armZZos[101];
   armZZos[101]=armZZos[9]*armZZos[101];
   armZZos[110]= - 10./3.*armZZos[34] - 1 + armZZos[165];
   armZZos[110]=armZZos[6]*armZZos[110];
   armZZos[113]=armZZos[14] - 32./9.*armZZos[15];
   armZZos[113]=armZZos[88]*armZZos[113];
   armZZos[101]=32./81.*armZZos[32] + 4./3.*armZZos[113] + 2./9.*
   armZZos[110] + 32./27.*armZZos[142] + armZZos[146] + armZZos[154] + 
   2*armZZos[101] + 583./81. + armZZos[103];
   armZZos[101]=armZZos[32]*armZZos[101];
   armZZos[110]= - 7 + 11*armZZos[32];
   armZZos[110]=armZZos[33]*armZZos[110];
   armZZos[113]=pow(armZZos[33],2);
   armZZos[114]=armZZos[5]*armZZos[113];
   armZZos[116]=36*armZZos[114];
   armZZos[110]=armZZos[110] + armZZos[116];
   armZZos[110]=2*armZZos[5]*armZZos[110];
   armZZos[118]=5113 + 2111*armZZos[32];
   armZZos[118]=armZZos[32]*armZZos[118];
   armZZos[118]=1501 + 1./2.*armZZos[118];
   armZZos[133]=mmt*armZZos[151]*armZZos[33];
   armZZos[143]=72*armZZos[133];
   armZZos[118]=armZZos[143] + 1./162.*armZZos[118] + armZZos[110];
   armZZos[118]=mmt*armZZos[118];
   armZZos[146]=1501 + 5113./4.*armZZos[32];
   armZZos[146]=armZZos[33]*armZZos[146];
   armZZos[152]=4*armZZos[114];
   armZZos[118]=armZZos[118] + 1./81.*armZZos[146] + armZZos[152];
   armZZos[118]=mmt*armZZos[118];
   armZZos[146]= - 7*armZZos[32];
   armZZos[154]= - 61 + armZZos[146];
   armZZos[154]=armZZos[33]*armZZos[154];
   armZZos[116]=1./3.*armZZos[154] + armZZos[116];
   armZZos[116]=armZZos[5]*armZZos[116];
   armZZos[146]= - 41 + armZZos[146];
   armZZos[146]=armZZos[32]*armZZos[146];
   armZZos[146]= - 17 + 1./2.*armZZos[146];
   armZZos[133]=36*armZZos[133];
   armZZos[116]=armZZos[133] + 7./324.*armZZos[146] + armZZos[116];
   armZZos[116]=mmt*armZZos[116];
   armZZos[146]= - 17 - 41./4.*armZZos[32];
   armZZos[146]=armZZos[33]*armZZos[146];
   armZZos[154]= - armZZos[5]*armZZos[113];
   armZZos[146]=7./54.*armZZos[146] + 34*armZZos[154];
   armZZos[116]=1./3.*armZZos[146] + armZZos[116];
   armZZos[116]=mmt*armZZos[116];
   armZZos[116]= - 119./324.*armZZos[113] + armZZos[116];
   armZZos[116]=armZZos[100]*armZZos[116];
   armZZos[116]=armZZos[116] + 1501./162.*armZZos[113] + armZZos[118];
   armZZos[116]=armZZos[100]*armZZos[116];
   armZZos[118]=113 + armZZos[129];
   armZZos[118]=armZZos[32]*armZZos[118];
   armZZos[118]=53 + 1./2.*armZZos[118];
   armZZos[110]=armZZos[143] + 1./18.*armZZos[118] + armZZos[110];
   armZZos[110]=mmt*armZZos[110];
   armZZos[118]=53 + 113./4.*armZZos[32];
   armZZos[118]=armZZos[33]*armZZos[118];
   armZZos[110]=armZZos[110] + 1./9.*armZZos[118] + armZZos[152];
   armZZos[110]=mmt*armZZos[110];
   armZZos[118]= - 5 + armZZos[32];
   armZZos[118]=armZZos[33]*armZZos[118];
   armZZos[114]=armZZos[118] + 12*armZZos[114];
   armZZos[114]=armZZos[5]*armZZos[114];
   armZZos[118]=1 - armZZos[32];
   armZZos[118]=armZZos[32]*armZZos[118];
   armZZos[118]=1 + 1./2.*armZZos[118];
   armZZos[114]=armZZos[133] + 1./4.*armZZos[118] + 3*armZZos[114];
   armZZos[114]=mmt*armZZos[114];
   armZZos[118]=1 + 1./4.*armZZos[32];
   armZZos[118]=armZZos[33]*armZZos[118];
   armZZos[114]=armZZos[114] + 1./2.*armZZos[118] + 6*armZZos[154];
   armZZos[114]=mmt*armZZos[114];
   armZZos[114]=1./4.*armZZos[113] + armZZos[114];
   armZZos[114]=armZZos[112]*armZZos[114];
   armZZos[110]=armZZos[114] + 53./18.*armZZos[113] + armZZos[110];
   armZZos[110]=armZZos[112]*armZZos[110];
   armZZos[114]= - 1 - armZZos[32];
   armZZos[114]=armZZos[33]*armZZos[114];
   armZZos[118]= - 2 - armZZos[32];
   armZZos[118]=armZZos[32]*armZZos[118];
   armZZos[118]= - 1 + armZZos[118];
   armZZos[118]=mmt*armZZos[118];
   armZZos[114]=2*armZZos[114] + armZZos[118];
   armZZos[114]=mmt*armZZos[114];
   armZZos[113]= - armZZos[113] + armZZos[114];
   armZZos[110]=armZZos[110] + 512./27.*armZZos[113] + armZZos[116];
   armZZos[110]=armZZos[2]*armZZos[110];
   armZZos[113]= - 583./3. + 61./2.*armZZos[6];
   armZZos[113]=5*armZZos[113] + 1073./8.*armZZos[32];
   armZZos[113]=1./3.*armZZos[113] + 425*armZZos[262];
   armZZos[113]=armZZos[33]*armZZos[113];
   armZZos[114]= - 1567./6. + 61*armZZos[6];
   armZZos[116]= - 49./4.*armZZos[32] + 2057./12. + 35*armZZos[6];
   armZZos[116]=armZZos[32]*armZZos[116];
   armZZos[114]=5*armZZos[114] + armZZos[116];
   armZZos[116]= - 5*armZZos[6];
   armZZos[118]= - 17*armZZos[32] + 88./3. + armZZos[116];
   armZZos[118]=armZZos[5]*armZZos[33]*armZZos[118];
   armZZos[114]=1./216.*armZZos[114] + armZZos[118];
   armZZos[114]=mmt*armZZos[114];
   armZZos[113]=1./36.*armZZos[113] + armZZos[114];
   armZZos[113]=armZZos[100]*armZZos[113];
   armZZos[114]=100./3. + 13./4.*armZZos[6];
   armZZos[114]=11*armZZos[114] + 10793./16.*armZZos[32];
   armZZos[114]=1./3.*armZZos[114] + 755./2.*armZZos[265];
   armZZos[114]=armZZos[33]*armZZos[114];
   armZZos[118]=6911./4.*armZZos[32] + 30089./12. - 19*armZZos[6];
   armZZos[118]=armZZos[32]*armZZos[118];
   armZZos[118]=armZZos[118] + 2005./6. + 143*armZZos[6];
   armZZos[129]= - armZZos[6] + armZZos[32];
   armZZos[129]=2*armZZos[5]*armZZos[33]*armZZos[129];
   armZZos[118]=1./324.*armZZos[118] + armZZos[129];
   armZZos[118]=mmt*armZZos[118];
   armZZos[113]=1./3.*armZZos[113] + 1./27.*armZZos[114] + armZZos[118]
   ;
   armZZos[113]=armZZos[100]*armZZos[113];
   armZZos[114]= - 1./4.*armZZos[32] + 89./12. - armZZos[6];
   armZZos[114]=armZZos[32]*armZZos[114];
   armZZos[114]=armZZos[114] - 59./6. + armZZos[212];
   armZZos[118]= - 3*armZZos[6];
   armZZos[133]=armZZos[179] + 8 + armZZos[118];
   armZZos[133]=armZZos[5]*armZZos[33]*armZZos[133];
   armZZos[114]=1./8.*armZZos[114] + armZZos[133];
   armZZos[114]=mmt*armZZos[114];
   armZZos[133]=3*armZZos[262] + 17./8.*armZZos[32] - 17./3. + 5./2.*
   armZZos[6];
   armZZos[133]=armZZos[33]*armZZos[133];
   armZZos[114]=1./4.*armZZos[133] + armZZos[114];
   armZZos[114]=armZZos[112]*armZZos[114];
   armZZos[133]=289./16.*armZZos[32] + 20 + 23./4.*armZZos[6];
   armZZos[133]=1./9.*armZZos[133] + 9./2.*armZZos[265];
   armZZos[133]=armZZos[33]*armZZos[133];
   armZZos[143]=79./2. + 23*armZZos[6];
   armZZos[146]=37./4. + 1./3.*armZZos[6];
   armZZos[146]=5*armZZos[146] + 199./12.*armZZos[32];
   armZZos[146]=armZZos[32]*armZZos[146];
   armZZos[143]=1./3.*armZZos[143] + armZZos[146];
   armZZos[129]=1./12.*armZZos[143] + armZZos[129];
   armZZos[129]=mmt*armZZos[129];
   armZZos[114]=armZZos[114] + armZZos[133] + armZZos[129];
   armZZos[114]=armZZos[112]*armZZos[114];
   armZZos[129]= - 17./3. - armZZos[6];
   armZZos[133]= - 2*armZZos[32];
   armZZos[129]=1./3.*armZZos[129] + armZZos[133];
   armZZos[129]=armZZos[32]*armZZos[129];
   armZZos[129]=1./3.*armZZos[178] + armZZos[129];
   armZZos[129]=mmt*armZZos[129];
   armZZos[143]=8*armZZos[262] - 8*armZZos[32] - 5./3. - armZZos[6];
   armZZos[143]=armZZos[33]*armZZos[143];
   armZZos[129]=1./3.*armZZos[143] + armZZos[129];
   armZZos[110]=armZZos[110] + armZZos[114] + 64./9.*armZZos[129] + 
   armZZos[113];
   armZZos[110]=armZZos[2]*armZZos[110];
   armZZos[113]= - 1 - armZZos[6];
   armZZos[113]=1./3.*armZZos[113] + armZZos[133];
   armZZos[113]=armZZos[32]*armZZos[113];
   armZZos[114]=5./3. - armZZos[6];
   armZZos[114]=armZZos[88]*armZZos[114];
   armZZos[114]=armZZos[114] + 8*armZZos[140];
   armZZos[114]=armZZos[33]*armZZos[114];
   armZZos[129]=17./3. - armZZos[6];
   armZZos[129]=armZZos[88]*armZZos[129];
   armZZos[129]=1./3.*armZZos[129] + 2*armZZos[140];
   armZZos[129]=armZZos[32]*armZZos[129];
   armZZos[133]=1./3. + armZZos[6];
   armZZos[133]=armZZos[88]*armZZos[133];
   armZZos[129]=1./3.*armZZos[133] + armZZos[129];
   armZZos[129]=mmZ*armZZos[129];
   armZZos[133]=14./3. - armZZos[6];
   armZZos[113]=armZZos[129] + 4./3.*armZZos[114] + 2./9.*armZZos[133]
    + armZZos[113];
   armZZos[114]=7*armZZos[6];
   armZZos[129]=283./24. + armZZos[114];
   armZZos[129]=armZZos[6]*armZZos[129];
   armZZos[133]=6263./4.*armZZos[32] + 22225./12. - 89*armZZos[6];
   armZZos[133]=armZZos[32]*armZZos[133];
   armZZos[129]=1./8.*armZZos[133] + armZZos[129] - 31649./144. - 7*
   armZZos[137];
   armZZos[133]= - 44./3. + 25./4.*armZZos[6];
   armZZos[133]=armZZos[88]*armZZos[133];
   armZZos[133]=armZZos[133] + 1127./16.*armZZos[166];
   armZZos[133]=armZZos[33]*armZZos[133];
   armZZos[129]=1./6.*armZZos[129] + 5*armZZos[133];
   armZZos[133]= - 131./24. + armZZos[6];
   armZZos[133]=armZZos[6]*armZZos[133];
   armZZos[143]= - 8803./12. + 197*armZZos[6];
   armZZos[143]=5*armZZos[143] + 3623./4.*armZZos[32];
   armZZos[143]=armZZos[32]*armZZos[143];
   armZZos[133]=1./8.*armZZos[143] + 25*armZZos[133] + 32551./144. - 25
   *armZZos[137];
   armZZos[143]=11./3. - 5./2.*armZZos[6];
   armZZos[143]=armZZos[88]*armZZos[143];
   armZZos[143]=armZZos[143] + 109./8.*armZZos[140];
   armZZos[143]=armZZos[33]*armZZos[143];
   armZZos[133]=1./3.*armZZos[133] + 25*armZZos[143];
   armZZos[116]=211./12. + armZZos[116];
   armZZos[116]=armZZos[88]*armZZos[116];
   armZZos[116]=1./3.*armZZos[116] + 25./4.*armZZos[140];
   armZZos[116]=armZZos[32]*armZZos[116];
   armZZos[143]=7./6. + armZZos[212];
   armZZos[143]=armZZos[88]*armZZos[143];
   armZZos[116]=1./3.*armZZos[143] + armZZos[116];
   armZZos[116]=mmZ*armZZos[116];
   armZZos[116]=1./3.*armZZos[133] + 25./8.*armZZos[116];
   armZZos[116]=armZZos[100]*armZZos[116];
   armZZos[133]= - 2273./12. + 25*armZZos[6];
   armZZos[133]=armZZos[88]*armZZos[133];
   armZZos[133]=1./3.*armZZos[133] + 275./4.*armZZos[166];
   armZZos[133]=armZZos[32]*armZZos[133];
   armZZos[143]= - 101./6. - 25*armZZos[6];
   armZZos[143]=armZZos[88]*armZZos[143];
   armZZos[133]=1./3.*armZZos[143] + armZZos[133];
   armZZos[133]=mmZ*armZZos[133];
   armZZos[116]=1./4.*armZZos[116] + 1./3.*armZZos[129] + 5./16.*
   armZZos[133];
   armZZos[116]=armZZos[100]*armZZos[116];
   armZZos[113]=16*armZZos[113] + armZZos[116];
   armZZos[116]=1 - 3./2.*armZZos[6];
   armZZos[116]=armZZos[88]*armZZos[116];
   armZZos[116]=armZZos[116] + 15./8.*armZZos[140];
   armZZos[116]=armZZos[33]*armZZos[116];
   armZZos[118]=11./4. + armZZos[118];
   armZZos[118]=armZZos[88]*armZZos[118];
   armZZos[118]=armZZos[118] + 9./4.*armZZos[140];
   armZZos[118]=armZZos[32]*armZZos[118];
   armZZos[119]= - 1./2. + armZZos[119];
   armZZos[119]=armZZos[88]*armZZos[119];
   armZZos[118]=armZZos[119] + armZZos[118];
   armZZos[118]=mmZ*armZZos[118];
   armZZos[119]= - 71./24. + armZZos[6];
   armZZos[119]=armZZos[6]*armZZos[119];
   armZZos[129]=23./4.*armZZos[32] - 287./12. + 13*armZZos[6];
   armZZos[129]=armZZos[32]*armZZos[129];
   armZZos[116]=1./8.*armZZos[118] + armZZos[116] + 1./8.*armZZos[129]
    + armZZos[119] + 247./144. - armZZos[137];
   armZZos[116]=armZZos[112]*armZZos[116];
   armZZos[118]=65./8. - armZZos[6];
   armZZos[118]=armZZos[6]*armZZos[118];
   armZZos[119]=127./4.*armZZos[32] + 643./4. - 17*armZZos[6];
   armZZos[119]=armZZos[32]*armZZos[119];
   armZZos[118]=1./8.*armZZos[119] + armZZos[118] - 209./16. + 
   armZZos[137];
   armZZos[119]= - 4./3. + 7./4.*armZZos[6];
   armZZos[119]=armZZos[88]*armZZos[119];
   armZZos[119]=armZZos[119] + 57./16.*armZZos[166];
   armZZos[119]=armZZos[33]*armZZos[119];
   armZZos[114]= - 127./12. + armZZos[114];
   armZZos[114]=armZZos[88]*armZZos[114];
   armZZos[114]=armZZos[114] + 39./4.*armZZos[166];
   armZZos[114]=armZZos[32]*armZZos[114];
   armZZos[129]=5./6. - 7*armZZos[6];
   armZZos[129]=armZZos[88]*armZZos[129];
   armZZos[114]=armZZos[129] + armZZos[114];
   armZZos[114]=mmZ*armZZos[114];
   armZZos[114]=1./4.*armZZos[116] + 1./16.*armZZos[114] + 1./18.*
   armZZos[118] + armZZos[119];
   armZZos[114]=armZZos[112]*armZZos[114];
   armZZos[110]=armZZos[110] + 1./9.*armZZos[113] + armZZos[114];
   armZZos[110]=armZZos[30]*armZZos[110];
   armZZos[111]=armZZos[150] + armZZos[130] + armZZos[122] + 
   armZZos[111] - 293./9. + armZZos[103];
   armZZos[111]=armZZos[88]*armZZos[111];
   armZZos[113]= - armZZos[31] + 2*armZZos[88];
   armZZos[113]=armZZos[32]*armZZos[113];
   armZZos[111]=128./27.*armZZos[113] + 16./3.*armZZos[111] + 
   armZZos[87] + 328./9.*armZZos[31];
   armZZos[111]=armZZos[33]*armZZos[111];
   armZZos[113]=14./9. + armZZos[123];
   armZZos[113]=armZZos[9]*armZZos[113];
   armZZos[113]= - 10./9.*armZZos[34] + armZZos[199] + armZZos[113] + 
   26./27. + armZZos[123];
   armZZos[113]=2*armZZos[113] + 5./27.*armZZos[6];
   armZZos[113]=armZZos[6]*armZZos[113];
   armZZos[114]= - 11 - 29./3.*armZZos[54];
   armZZos[114]=armZZos[123]*armZZos[114];
   armZZos[116]=1 + armZZos[32];
   armZZos[116]=armZZos[5]*armZZos[116];
   armZZos[116]=32./9.*armZZos[116] + armZZos[241] + 2./3.*armZZos[40]
    - 256./81.*armZZos[47] - armZZos[41];
   armZZos[116]=mmt*armZZos[116];
   armZZos[118]=23./9.*armZZos[39] + 16./9.*armZZos[35] - 32./9.*
   armZZos[58] - 2*armZZos[59] + armZZos[37];
   armZZos[118]=armZZos[31]*armZZos[118];
   armZZos[119]= - 4*armZZos[123];
   armZZos[122]=617./9. + armZZos[119];
   armZZos[122]=armZZos[9]*armZZos[122];
   armZZos[129]=armZZos[87] + 16./3.*armZZos[31];
   armZZos[129]=armZZos[14]*armZZos[129];
   armZZos[130]= - 64./9.*armZZos[18] - 14*armZZos[39] - 32./9.*
   armZZos[35] + armZZos[211] + 32./3.*armZZos[58];
   armZZos[130]=16./27.*armZZos[15] + 1./3.*armZZos[130] - armZZos[14];
   armZZos[130]=armZZos[88]*armZZos[130];
   armZZos[133]=armZZos[5]*armZZos[33];
   armZZos[92]=armZZos[110] + armZZos[96] + armZZos[92] + armZZos[99]
    + 2*armZZos[116] + armZZos[109] + 64./9.*armZZos[133] + 
   armZZos[111] + 4*armZZos[101] + 16*armZZos[130] + 4./3.*armZZos[113]
    + 128./27.*armZZos[142] + 560./243.*armZZos[34] + 112./81.*
   armZZos[1] + armZZos[129] + 2./9.*armZZos[122] + 128./27.*
   armZZos[172] + 8./3.*armZZos[118] + 2./3.*armZZos[114] + 1./3.*
   armZZos[171] + 1./3.*armZZos[170] + 2./3.*armZZos[169] - 20./81.*
   armZZos[137] + 10*armZZos[12] + 520./81.*armZZos[38] - 5./3.*
   armZZos[55] - 107./9.*armZZos[54] - 640./81.*armZZos[51] - 8./9.*
   armZZos[52] + 80./81.*armZZos[28] - 2243./162. - 10*armZZos[49];
   armZZos[92]=armZZos[30]*armZZos[92];
   armZZos[96]=armZZos[224] + armZZos[190];
   armZZos[96]=1./3.*armZZos[10]*armZZos[96];
   armZZos[99]=10*armZZos[80];
   armZZos[101]= - 4*armZZos[82];
   armZZos[109]=43./12.*armZZos[65];
   armZZos[110]=6*armZZos[17];
   armZZos[111]= - 157 + armZZos[213];
   armZZos[111]=armZZos[16]*armZZos[111];
   armZZos[113]=43./3.*armZZos[10];
   armZZos[114]=armZZos[113] + 163./3. + 11*armZZos[9];
   armZZos[114]=armZZos[15]*armZZos[114];
   armZZos[116]=9*armZZos[66];
   armZZos[118]= - 11*armZZos[81];
   armZZos[111]=1./2.*armZZos[114] + armZZos[96] + 1./4.*armZZos[111]
    + 100./3.*armZZos[14] - 56./3.*armZZos[18] - 293./6.*armZZos[19] + 
   armZZos[110] + armZZos[109] + armZZos[101] + armZZos[64] + 
   armZZos[99] + 11./2.*armZZos[63] + armZZos[116] + armZZos[118];
   armZZos[111]=armZZos[89]*armZZos[111];
   armZZos[114]= - 5 + armZZos[213];
   armZZos[114]=1./2.*armZZos[114] + armZZos[144];
   armZZos[114]=armZZos[10]*armZZos[114];
   armZZos[122]= - 11*armZZos[74];
   armZZos[129]= - 55*armZZos[79] + 541./4. + armZZos[122];
   armZZos[130]=18*armZZos[75];
   armZZos[133]= - 4*armZZos[73];
   armZZos[140]=2*armZZos[85];
   armZZos[142]=11./2.*armZZos[9];
   armZZos[143]=3*armZZos[76];
   armZZos[114]=armZZos[114] + armZZos[142] - 27./2.*armZZos[13] - 33*
   armZZos[78] + armZZos[140] + armZZos[133] + armZZos[130] + 1./2.*
   armZZos[129] + armZZos[143];
   armZZos[114]=armZZos[89]*armZZos[114];
   armZZos[129]= - armZZos[14] + armZZos[15];
   armZZos[129]=armZZos[5]*armZZos[129];
   armZZos[91]=21./2.*armZZos[129] + armZZos[181] + 1015./9. + 
   armZZos[91];
   armZZos[91]=armZZos[5]*armZZos[91];
   armZZos[129]= - 2*armZZos[84] - 1./2.*armZZos[72];
   armZZos[129]=11*armZZos[129] + 255./8.*armZZos[71];
   armZZos[146]= - 2*armZZos[67];
   armZZos[150]=44*armZZos[12] - 88*armZZos[78] - 2*armZZos[85] + 44*
   armZZos[79] - 41./4. - 44*armZZos[86];
   armZZos[150]=armZZos[90]*armZZos[150];
   armZZos[152]=armZZos[13]*armZZos[90];
   armZZos[154]=armZZos[10]*armZZos[90];
   armZZos[156]= - mmZ*armZZos[151];
   armZZos[159]=armZZos[9]*armZZos[90];
   armZZos[91]=15*armZZos[156] + armZZos[114] + armZZos[91] + 2*
   armZZos[154] + 44*armZZos[159] + 2*armZZos[152] + armZZos[150] + 
   armZZos[146] - 22*armZZos[70] + 3*armZZos[129] - 55*armZZos[69];
   armZZos[91]=mmZ*armZZos[91];
   armZZos[114]=11*armZZos[63] + 3./2.*armZZos[65];
   armZZos[114]= - 7*armZZos[18] - 25./4.*armZZos[19] + 1./2.*
   armZZos[114] - 11*armZZos[17];
   armZZos[129]= - 4./3. + 15./4.*armZZos[9];
   armZZos[129]=armZZos[14]*armZZos[129];
   armZZos[150]=19 - 33./2.*armZZos[9];
   armZZos[150]=armZZos[16]*armZZos[150];
   armZZos[152]=armZZos[136] - 2*armZZos[16];
   armZZos[152]=armZZos[10]*armZZos[152];
   armZZos[154]=armZZos[107] + 41./3. + 33./2.*armZZos[9];
   armZZos[154]=armZZos[15]*armZZos[154];
   armZZos[156]=9*armZZos[5]*armZZos[15]*armZZos[117];
   armZZos[114]=armZZos[156] + 1./2.*armZZos[154] + armZZos[152] + 
   armZZos[150] + 3*armZZos[114] + 11*armZZos[129];
   armZZos[114]=armZZos[5]*armZZos[114];
   armZZos[129]= - 3*armZZos[68];
   armZZos[150]=4./3.*armZZos[67];
   armZZos[154]=armZZos[150] + 67./12.*armZZos[70] + armZZos[129] + 89./
   6.*armZZos[69];
   armZZos[154]=mmH*armZZos[154];
   armZZos[166]=152789./324. - 173*armZZos[77];
   armZZos[167]= - 8./3.*armZZos[76];
   armZZos[168]=9./2.*armZZos[75];
   armZZos[169]= - 79./12.*armZZos[19] + 227./6.*armZZos[17] - 43./6.*
   armZZos[64] + 2*armZZos[82];
   armZZos[169]=armZZos[90]*armZZos[169];
   armZZos[170]=97./2. + 21*armZZos[9];
   armZZos[170]=armZZos[9]*armZZos[170];
   armZZos[171]= - 56./3.*armZZos[90] + 55./2.*armZZos[159];
   armZZos[171]=armZZos[14]*armZZos[171];
   armZZos[172]= - armZZos[9]*armZZos[90];
   armZZos[173]=251./3.*armZZos[90] + 55*armZZos[172];
   armZZos[173]=armZZos[16]*armZZos[173];
   armZZos[178]= - armZZos[14]*armZZos[90];
   armZZos[179]=1./2.*armZZos[178];
   armZZos[180]= - armZZos[10] + armZZos[179] - 1 - 11./6.*armZZos[9];
   armZZos[180]=armZZos[10]*armZZos[180];
   armZZos[184]= - armZZos[18]*armZZos[90];
   armZZos[187]=armZZos[15]*armZZos[90];
   armZZos[91]=armZZos[154] + armZZos[91] + armZZos[111] + armZZos[114]
    + armZZos[187] + armZZos[180] + 1./4.*armZZos[173] + armZZos[171]
    + 7*armZZos[170] + 2*armZZos[184] + armZZos[218] + armZZos[169] + 
   1891./24.*armZZos[12] - 77./2.*armZZos[78] - 28./3.*armZZos[85] + 
   armZZos[168] + armZZos[167] + 51./2.*armZZos[79] + 22./3.*
   armZZos[86] + 5./8.*armZZos[166] + 11./3.*armZZos[74];
   armZZos[91]=armZZos[62]*armZZos[91];
   armZZos[111]=14065./36. + 33*armZZos[74];
   armZZos[111]=armZZos[143] + 1./2.*armZZos[111] + 33*armZZos[79];
   armZZos[114]=9*armZZos[75];
   armZZos[154]= - 2*armZZos[73];
   armZZos[166]=3*armZZos[9];
   armZZos[169]=25./9. + armZZos[166];
   armZZos[169]=11./4.*armZZos[169] - armZZos[10];
   armZZos[169]=armZZos[10]*armZZos[169];
   armZZos[111]=armZZos[169] - 33./4.*armZZos[9] - 73./4.*armZZos[13]
    + 99./2.*armZZos[78] + 3./2.*armZZos[85] + armZZos[154] + 1./2.*
   armZZos[111] + armZZos[114];
   armZZos[111]=armZZos[89]*armZZos[111];
   armZZos[94]=armZZos[94] - 1847./6. - 99*armZZos[9];
   armZZos[117]=armZZos[5]*armZZos[117];
   armZZos[94]=1./2.*armZZos[94] + 63*armZZos[117];
   armZZos[94]=armZZos[5]*armZZos[94];
   armZZos[117]=11./4.*armZZos[70] - 153./16.*armZZos[71] + 11*
   armZZos[69];
   armZZos[169]=armZZos[85] - 11*armZZos[79] - 5./2. + 11*armZZos[86];
   armZZos[169]= - 11./2.*armZZos[12] + 1./2.*armZZos[169] + 11*
   armZZos[78];
   armZZos[169]=armZZos[90]*armZZos[169];
   armZZos[170]= - armZZos[13]*armZZos[90];
   armZZos[171]= - armZZos[10]*armZZos[90];
   armZZos[151]=mmZ*armZZos[151];
   armZZos[94]=45./2.*armZZos[151] + armZZos[111] + 1./2.*armZZos[94]
    + 3./2.*armZZos[171] + 33./2.*armZZos[172] + 3./2.*armZZos[170] + 3
   *armZZos[169] + 3*armZZos[117] - armZZos[67];
   armZZos[94]=mmZ*armZZos[94];
   armZZos[111]= - 3*armZZos[9];
   armZZos[117]=17./9. + armZZos[111];
   armZZos[117]=11*armZZos[117] + armZZos[113];
   armZZos[117]=armZZos[15]*armZZos[117];
   armZZos[169]=3*armZZos[66];
   armZZos[173]= - 11./2.*armZZos[63] + armZZos[169] + 11*armZZos[81];
   armZZos[180]=5*armZZos[80];
   armZZos[188]= - 2*armZZos[82];
   armZZos[189]=43./24.*armZZos[65];
   armZZos[190]= - 3073./9. + 33*armZZos[9];
   armZZos[190]=armZZos[16]*armZZos[190];
   armZZos[138]=armZZos[138] - 7./6.*armZZos[16];
   armZZos[138]=armZZos[10]*armZZos[138];
   armZZos[117]=1./4.*armZZos[117] + armZZos[138] + 1./8.*armZZos[190]
    - 69./2.*armZZos[14] + 5./3.*armZZos[18] - 181./6.*armZZos[19] + 5./
   2.*armZZos[17] + armZZos[189] + armZZos[188] + 3./4.*armZZos[64] + 3.
   /2.*armZZos[173] + armZZos[180];
   armZZos[117]=armZZos[89]*armZZos[117];
   armZZos[138]=1./2.*armZZos[65];
   armZZos[173]=armZZos[138] - 11*armZZos[63] + armZZos[64];
   armZZos[190]= - 11./2.*armZZos[9];
   armZZos[191]= - 1 + armZZos[190];
   armZZos[191]=armZZos[14]*armZZos[191];
   armZZos[173]=1./2.*armZZos[191] + 5./2.*armZZos[18] + armZZos[19] + 
   1./4.*armZZos[173] + 5*armZZos[17];
   armZZos[191]=armZZos[14] + 1./4.*armZZos[15];
   armZZos[191]=armZZos[15]*armZZos[191];
   armZZos[193]=pow(armZZos[14],2);
   armZZos[191]=armZZos[193] + armZZos[191];
   armZZos[191]=armZZos[5]*armZZos[191];
   armZZos[148]= - 85 + armZZos[148];
   armZZos[148]=armZZos[16]*armZZos[148];
   armZZos[194]=armZZos[136] - armZZos[16];
   armZZos[194]=armZZos[10]*armZZos[194];
   armZZos[198]=armZZos[107] - 7 - 99./2.*armZZos[9];
   armZZos[198]=armZZos[15]*armZZos[198];
   armZZos[148]=9*armZZos[191] + 1./4.*armZZos[198] + armZZos[194] + 9*
   armZZos[173] + 1./8.*armZZos[148];
   armZZos[148]=armZZos[5]*armZZos[148];
   armZZos[173]=armZZos[14]*armZZos[90];
   armZZos[191]=armZZos[181] + 3./2.*armZZos[173] + 41./9. + 
   armZZos[142];
   armZZos[191]=armZZos[10]*armZZos[191];
   armZZos[194]=48019./162. + 435*armZZos[77];
   armZZos[122]=1./2.*armZZos[194] + armZZos[122];
   armZZos[122]=1./2.*armZZos[122] - 11*armZZos[86];
   armZZos[194]=21./4.*armZZos[19] - 37./2.*armZZos[17] + 3*armZZos[64]
    + armZZos[188];
   armZZos[194]=armZZos[90]*armZZos[194];
   armZZos[198]=armZZos[18]*armZZos[90];
   armZZos[201]= - 327*armZZos[9];
   armZZos[208]= - 1697./3. + armZZos[201];
   armZZos[208]=armZZos[9]*armZZos[208];
   armZZos[209]=5*armZZos[90] + 11./2.*armZZos[172];
   armZZos[209]=armZZos[14]*armZZos[209];
   armZZos[211]= - 15*armZZos[90] + 11*armZZos[159];
   armZZos[211]=armZZos[16]*armZZos[211];
   armZZos[212]= - armZZos[15]*armZZos[90];
   armZZos[215]=2./3.*armZZos[67];
   armZZos[218]=armZZos[215] - 13./4.*armZZos[70] - 3./2.*armZZos[68]
    - 13*armZZos[69];
   armZZos[218]=mmH*armZZos[218];
   armZZos[219]=9./4.*armZZos[75];
   armZZos[94]=armZZos[218] + armZZos[94] + armZZos[117] + armZZos[148]
    + armZZos[212] + armZZos[191] + 3./4.*armZZos[211] + 3*armZZos[209]
    + 1./16.*armZZos[208] + 2*armZZos[198] - 7./6.*armZZos[13] + 
   armZZos[194] - 391./8.*armZZos[12] + 33*armZZos[78] + 3*armZZos[85]
    + armZZos[219] + 1./2.*armZZos[122] - 4./3.*armZZos[76];
   armZZos[94]=armZZos[62]*armZZos[94];
   armZZos[117]=37*armZZos[29] + 11053./432. - 33*armZZos[26];
   armZZos[122]= - 7./3.*armZZos[10];
   armZZos[148]=4./9. - armZZos[137];
   armZZos[191]=armZZos[1]*armZZos[148];
   armZZos[148]=armZZos[34]*armZZos[148];
   armZZos[117]=armZZos[148] + 2./3.*armZZos[191] + armZZos[122] + 121./
   4.*armZZos[9] - armZZos[13] + 5./4.*armZZos[137] + 33./2.*
   armZZos[12] + armZZos[135] + armZZos[25] + 1./2.*armZZos[117] + 
   armZZos[134];
   armZZos[117]=armZZos[34]*armZZos[117];
   armZZos[134]= - 1./4.*armZZos[15];
   armZZos[194]= - 1./2.*armZZos[82];
   armZZos[198]=armZZos[134] + armZZos[136] + 1./2.*armZZos[18] + 
   armZZos[194] + armZZos[17];
   armZZos[198]=armZZos[89]*armZZos[198];
   armZZos[208]=1./2.*armZZos[82] - armZZos[17];
   armZZos[208]=armZZos[90]*armZZos[208];
   armZZos[209]=armZZos[90] - armZZos[89];
   armZZos[209]=mmZ*armZZos[209];
   armZZos[184]=3./16.*armZZos[209] + armZZos[198] + 1./4.*armZZos[187]
    + 3./2.*armZZos[178] + armZZos[208] + 1./2.*armZZos[184];
   armZZos[184]=armZZos[112]*armZZos[62]*armZZos[184];
   armZZos[187]=37./3.*armZZos[29] + 11053./1296. - 11*armZZos[26];
   armZZos[187]=1./9.*armZZos[191] - 7./9.*armZZos[10] + 121./12.*
   armZZos[9] - 1./3.*armZZos[13] + armZZos[204] + 11./2.*armZZos[12]
    + armZZos[202] + armZZos[200] + 1./2.*armZZos[187] + 5./3.*
   armZZos[28];
   armZZos[187]=armZZos[1]*armZZos[187];
   armZZos[198]=2*armZZos[7];
   armZZos[102]=armZZos[102] + armZZos[198] + armZZos[8];
   armZZos[126]=armZZos[127] + armZZos[126] + armZZos[125] + 3*
   armZZos[102] + armZZos[124];
   armZZos[126]=armZZos[34]*armZZos[126];
   armZZos[200]= - 1./6.*armZZos[19];
   armZZos[202]=1./6.*armZZos[16];
   armZZos[102]=armZZos[202] + armZZos[192] - armZZos[18] + 
   armZZos[102] + armZZos[200];
   armZZos[102]=armZZos[1]*armZZos[102];
   armZZos[192]=armZZos[34]*armZZos[121];
   armZZos[204]= - 1./6.*armZZos[16];
   armZZos[208]=armZZos[14] + armZZos[204];
   armZZos[208]=armZZos[1]*armZZos[208];
   armZZos[105]=armZZos[15]*armZZos[105];
   armZZos[105]=1./2.*armZZos[105] + armZZos[208] + armZZos[192];
   armZZos[105]=armZZos[6]*armZZos[105];
   armZZos[95]=armZZos[15]*armZZos[95];
   armZZos[95]=armZZos[105] + 3./2.*armZZos[95] + armZZos[102] + 
   armZZos[126];
   armZZos[95]=armZZos[5]*armZZos[95];
   armZZos[102]=armZZos[204] - armZZos[18] + armZZos[8] + armZZos[131];
   armZZos[105]=armZZos[1]*armZZos[102];
   armZZos[102]=armZZos[34]*armZZos[102];
   armZZos[126]= - armZZos[1]*armZZos[16];
   armZZos[131]= - armZZos[34]*armZZos[16];
   armZZos[192]=1./3.*armZZos[126] + armZZos[131];
   armZZos[204]=armZZos[15]*armZZos[160];
   armZZos[208]=1./2.*armZZos[192] + armZZos[204];
   armZZos[208]=armZZos[6]*armZZos[208];
   armZZos[175]=armZZos[15]*armZZos[175];
   armZZos[208]=armZZos[208] + 2./3.*armZZos[175] + 1./3.*armZZos[105]
    + armZZos[102];
   armZZos[208]=armZZos[89]*armZZos[208];
   armZZos[209]= - armZZos[25] + 7./3. + 2*armZZos[27];
   armZZos[186]=armZZos[206] + armZZos[186] + 1./3.*armZZos[209] + 
   armZZos[11];
   armZZos[206]=armZZos[1]*armZZos[186];
   armZZos[97]=armZZos[97] + armZZos[13] + armZZos[209] + 3*armZZos[11]
   ;
   armZZos[209]=armZZos[34]*armZZos[97];
   armZZos[211]=armZZos[1]*armZZos[98];
   armZZos[98]=armZZos[34]*armZZos[98];
   armZZos[212]=1./3.*armZZos[211] + armZZos[98];
   armZZos[212]=armZZos[6]*armZZos[212];
   armZZos[206]=armZZos[212] + armZZos[206] + armZZos[209];
   armZZos[206]=armZZos[89]*armZZos[206];
   armZZos[209]=armZZos[21] + 1./3.*armZZos[23];
   armZZos[209]=2*armZZos[209] - 1./3.*armZZos[22];
   armZZos[209]=armZZos[1]*armZZos[209];
   armZZos[212]=3*armZZos[21] + armZZos[23];
   armZZos[212]=2*armZZos[212] - armZZos[22];
   armZZos[212]=armZZos[34]*armZZos[212];
   armZZos[209]=armZZos[209] + armZZos[212];
   armZZos[106]=7./2.*armZZos[160] + armZZos[106];
   armZZos[106]=armZZos[5]*armZZos[106];
   armZZos[106]=armZZos[206] + 2*armZZos[209] + armZZos[106];
   armZZos[106]=mmZ*armZZos[106];
   armZZos[206]= - 7 - armZZos[1];
   armZZos[206]=armZZos[1]*armZZos[206];
   armZZos[165]= - 7 + armZZos[165];
   armZZos[165]=1./3.*armZZos[165] - armZZos[34];
   armZZos[165]=armZZos[34]*armZZos[165];
   armZZos[165]=1./9.*armZZos[206] + armZZos[165];
   armZZos[93]= - 5./4. + armZZos[93];
   armZZos[93]=armZZos[1]*armZZos[93];
   armZZos[209]=armZZos[34] - 5./4. + armZZos[210];
   armZZos[209]=armZZos[34]*armZZos[209];
   armZZos[93]=1./3.*armZZos[93] + armZZos[209];
   armZZos[93]=armZZos[6]*armZZos[93];
   armZZos[93]=5./3.*armZZos[165] + armZZos[93];
   armZZos[93]=armZZos[6]*armZZos[93];
   armZZos[93]=3*armZZos[184] + armZZos[94] + armZZos[106] + 
   armZZos[208] + armZZos[95] + armZZos[93] + armZZos[187] + 
   armZZos[117];
   armZZos[93]=armZZos[112]*armZZos[93];
   armZZos[94]=28*armZZos[26];
   armZZos[95]= - 37./4.*armZZos[20];
   armZZos[106]= - 7./9. + armZZos[137];
   armZZos[117]=armZZos[1]*armZZos[106];
   armZZos[106]=armZZos[34]*armZZos[106];
   armZZos[165]=armZZos[106] + 2./3.*armZZos[117] - 55*armZZos[9] - 
   armZZos[137] + armZZos[203] - 10*armZZos[28] - 81./2.*armZZos[29] + 
   armZZos[95] - 7159./216. + armZZos[94];
   armZZos[165]=armZZos[34]*armZZos[165];
   armZZos[184]= - 28*armZZos[21];
   armZZos[187]= - 8*armZZos[23];
   armZZos[203]=4*armZZos[22];
   armZZos[208]=armZZos[203] + armZZos[187] + armZZos[184] - 
   armZZos[24];
   armZZos[208]=armZZos[1]*armZZos[208];
   armZZos[184]=armZZos[203] + armZZos[187] + armZZos[184] - 5./9.*
   armZZos[24];
   armZZos[184]=armZZos[34]*armZZos[184];
   armZZos[160]=1./3.*armZZos[160] + 4*armZZos[176];
   armZZos[160]=armZZos[5]*armZZos[160];
   armZZos[160]=armZZos[160] + 1./3.*armZZos[208] + armZZos[184];
   armZZos[160]=mmZ*armZZos[160];
   armZZos[94]=armZZos[95] - 7531./216. + armZZos[94];
   armZZos[94]=1./9.*armZZos[117] - 55./3.*armZZos[9] - 1./3.*
   armZZos[137] - 28./3.*armZZos[12] - 10./3.*armZZos[28] + 1./3.*
   armZZos[94] - 27./2.*armZZos[29];
   armZZos[94]=armZZos[1]*armZZos[94];
   armZZos[95]=1 + armZZos[174];
   armZZos[95]=armZZos[1]*armZZos[95];
   armZZos[174]= - armZZos[34] + 1 + armZZos[199];
   armZZos[174]=armZZos[34]*armZZos[174];
   armZZos[95]=1./3.*armZZos[95] + armZZos[174];
   armZZos[95]=armZZos[6]*armZZos[95];
   armZZos[174]= - 23./2.*armZZos[9];
   armZZos[176]=8./9.*armZZos[1] - 40./9. + armZZos[174];
   armZZos[176]=armZZos[1]*armZZos[176];
   armZZos[184]=8./3.*armZZos[34] + 16./9.*armZZos[1] - 34./9. + 
   armZZos[174];
   armZZos[184]=armZZos[34]*armZZos[184];
   armZZos[95]=armZZos[95] + 1./3.*armZZos[176] + armZZos[184];
   armZZos[95]=armZZos[6]*armZZos[95];
   armZZos[176]=armZZos[1]*armZZos[14];
   armZZos[184]=armZZos[34]*armZZos[14];
   armZZos[187]=1./3.*armZZos[176] + armZZos[184];
   armZZos[199]= - armZZos[1]*armZZos[14];
   armZZos[203]= - armZZos[34]*armZZos[14];
   armZZos[208]=armZZos[199] + 3*armZZos[203];
   armZZos[208]=armZZos[6]*armZZos[208];
   armZZos[187]=4*armZZos[187] + armZZos[208];
   armZZos[187]=armZZos[5]*armZZos[187];
   armZZos[91]=armZZos[93] + armZZos[91] + armZZos[160] + armZZos[187]
    + armZZos[95] + armZZos[94] + armZZos[165];
   armZZos[91]=armZZos[112]*armZZos[91];
   armZZos[93]=armZZos[113] + 115./3. + armZZos[166];
   armZZos[93]=armZZos[15]*armZZos[93];
   armZZos[94]=1./2.*armZZos[63] + armZZos[169] - armZZos[81];
   armZZos[95]= - 181 + armZZos[111];
   armZZos[95]=armZZos[16]*armZZos[95];
   armZZos[93]=1./2.*armZZos[93] + armZZos[96] + 1./4.*armZZos[95] + 52.
   /3.*armZZos[14] - 44./3.*armZZos[18] - 305./6.*armZZos[19] + 
   armZZos[110] + armZZos[109] + armZZos[101] + armZZos[64] + 3*
   armZZos[94] + armZZos[99];
   armZZos[93]=armZZos[89]*armZZos[93];
   armZZos[94]=11 + armZZos[111];
   armZZos[94]=1./2.*armZZos[94] + armZZos[144];
   armZZos[94]=armZZos[10]*armZZos[94];
   armZZos[95]= - 3*armZZos[74];
   armZZos[96]= - 1./3.*armZZos[79] + 2483./12. + armZZos[95];
   armZZos[99]=3./2.*armZZos[9];
   armZZos[94]=armZZos[94] + armZZos[99] - 35./2.*armZZos[13] - 9*
   armZZos[78] + armZZos[140] + armZZos[133] + armZZos[130] + 1./2.*
   armZZos[96] + armZZos[143];
   armZZos[94]=armZZos[89]*armZZos[94];
   armZZos[96]=armZZos[14] + armZZos[15];
   armZZos[96]=armZZos[5]*armZZos[96];
   armZZos[96]=21*armZZos[96] - armZZos[10] + 35./3. + armZZos[166];
   armZZos[96]=armZZos[5]*armZZos[96];
   armZZos[101]= - armZZos[72] + 23./4.*armZZos[71];
   armZZos[101]=1./2.*armZZos[101] - armZZos[69];
   armZZos[94]=5*armZZos[151] + armZZos[94] + 1./2.*armZZos[96] + 1./3.
   *armZZos[101] + armZZos[146];
   armZZos[94]=mmZ*armZZos[94];
   armZZos[96]=armZZos[63] + armZZos[138];
   armZZos[96]=1./2.*armZZos[96] - armZZos[17];
   armZZos[101]= - 3 + 1./4.*armZZos[9];
   armZZos[101]=armZZos[14]*armZZos[101];
   armZZos[109]=2 + armZZos[246];
   armZZos[109]=armZZos[16]*armZZos[109];
   armZZos[110]=armZZos[107] - 67./3. + 9./2.*armZZos[9];
   armZZos[110]=armZZos[15]*armZZos[110];
   armZZos[96]=armZZos[156] + 1./2.*armZZos[110] + armZZos[152] + 5./3.
   *armZZos[109] + armZZos[101] - 9*armZZos[18] + 9*armZZos[96] - 37./
   12.*armZZos[19];
   armZZos[96]=armZZos[5]*armZZos[96];
   armZZos[101]= - 1./12.*armZZos[70];
   armZZos[109]=armZZos[150] + armZZos[101] + armZZos[129] - 1./2.*
   armZZos[69];
   armZZos[109]=mmH*armZZos[109];
   armZZos[110]=70297./324. - 33*armZZos[77];
   armZZos[113]=armZZos[17] + armZZos[124];
   armZZos[113]=armZZos[90]*armZZos[113];
   armZZos[130]= - 1./2. - 4./3.*armZZos[9];
   armZZos[130]=armZZos[9]*armZZos[130];
   armZZos[133]= - 5*armZZos[90] + 1./2.*armZZos[159];
   armZZos[133]=armZZos[14]*armZZos[133];
   armZZos[144]=armZZos[90] + armZZos[172];
   armZZos[144]=armZZos[16]*armZZos[144];
   armZZos[146]=3*armZZos[173];
   armZZos[150]=armZZos[146] - 3 + 1./3.*armZZos[9];
   armZZos[150]=1./2.*armZZos[150] - armZZos[10];
   armZZos[150]=armZZos[10]*armZZos[150];
   armZZos[93]=armZZos[109] + armZZos[94] + armZZos[93] + armZZos[96]
    + armZZos[150] + 1./12.*armZZos[144] + 1./3.*armZZos[133] + 1./3.*
   armZZos[130] - 8./3.*armZZos[13] + 1./6.*armZZos[113] + 33./8.*
   armZZos[12] - 19./6.*armZZos[78] - 5./6.*armZZos[85] + armZZos[168]
    + armZZos[167] + 5./18.*armZZos[79] + 1./8.*armZZos[110] + 5./3.*
   armZZos[74];
   armZZos[93]=armZZos[62]*armZZos[93];
   armZZos[94]=armZZos[182] + armZZos[139];
   armZZos[94]=armZZos[10]*armZZos[94];
   armZZos[96]=1./6.*armZZos[63] + armZZos[116] - 1./3.*armZZos[81];
   armZZos[109]= - 1769./3. - armZZos[9];
   armZZos[109]=armZZos[16]*armZZos[109];
   armZZos[110]=43*armZZos[10] + 401./3. + armZZos[9];
   armZZos[110]=armZZos[15]*armZZos[110];
   armZZos[94]=1./12.*armZZos[110] + 1./3.*armZZos[94] + 1./24.*
   armZZos[109] + 11./3.*armZZos[14] - 20./3.*armZZos[18] - 105./4.*
   armZZos[19] + armZZos[108] + armZZos[189] + armZZos[188] + 
   armZZos[64] + 1./2.*armZZos[96] + armZZos[180];
   armZZos[94]=armZZos[89]*armZZos[94];
   armZZos[96]=7757./12. - armZZos[74];
   armZZos[96]=1./6.*armZZos[96] + armZZos[143];
   armZZos[109]=79./3. - armZZos[9];
   armZZos[109]=1./12.*armZZos[109] - armZZos[10];
   armZZos[109]=armZZos[10]*armZZos[109];
   armZZos[110]=1./12.*armZZos[9];
   armZZos[96]=armZZos[109] + armZZos[110] - 125./12.*armZZos[13] - 1./
   2.*armZZos[78] + armZZos[140] + armZZos[154] + 1./2.*armZZos[96] + 
   armZZos[114];
   armZZos[96]=armZZos[89]*armZZos[96];
   armZZos[109]=149./6. + armZZos[9];
   armZZos[113]=armZZos[5]*armZZos[15];
   armZZos[109]=21*armZZos[113] + 1./3.*armZZos[109] - armZZos[10];
   armZZos[109]=armZZos[5]*armZZos[109];
   armZZos[96]=5./2.*armZZos[151] + armZZos[96] + 1./4.*armZZos[109] - 
   3./16.*armZZos[71] - armZZos[67];
   armZZos[96]=mmZ*armZZos[96];
   armZZos[109]=armZZos[63] + 9./2.*armZZos[65];
   armZZos[113]=31 - armZZos[9];
   armZZos[113]=armZZos[16]*armZZos[113];
   armZZos[114]= - 5*armZZos[18];
   armZZos[109]=1./12.*armZZos[113] + 1./3.*armZZos[14] + armZZos[114]
    - 7./3.*armZZos[19] + 1./2.*armZZos[109] - armZZos[17];
   armZZos[107]=armZZos[107] - 31./3. + 1./2.*armZZos[9];
   armZZos[107]=armZZos[15]*armZZos[107];
   armZZos[113]=armZZos[5]*pow(armZZos[15],2);
   armZZos[130]= - armZZos[10]*armZZos[16];
   armZZos[107]=9./4.*armZZos[113] + 1./4.*armZZos[107] + 1./2.*
   armZZos[109] + armZZos[130];
   armZZos[107]=armZZos[5]*armZZos[107];
   armZZos[109]= - 1./3.*armZZos[69];
   armZZos[113]=armZZos[129] + armZZos[109];
   armZZos[113]=1./2.*armZZos[113] + armZZos[215];
   armZZos[113]=mmH*armZZos[113];
   armZZos[129]=1./2.*armZZos[74] + 41801./432. - 7*armZZos[77];
   armZZos[129]=1./2.*armZZos[129] - 4*armZZos[76];
   armZZos[133]=17./3. + 19*armZZos[9];
   armZZos[133]=armZZos[9]*armZZos[133];
   armZZos[139]=7./9. + armZZos[181];
   armZZos[139]=armZZos[10]*armZZos[139];
   armZZos[143]=1./2.*armZZos[85];
   armZZos[94]=armZZos[113] + armZZos[96] + armZZos[94] + armZZos[107]
    + armZZos[139] + 1./144.*armZZos[133] - 3./2.*armZZos[13] + 7./6.*
   armZZos[12] - 1./6.*armZZos[78] + armZZos[143] + 1./3.*armZZos[129]
    + armZZos[219];
   armZZos[94]=armZZos[62]*armZZos[94];
   armZZos[96]=armZZos[191] + armZZos[122] - 5./36.*armZZos[9] - 
   armZZos[13] - 15./4.*armZZos[137] + armZZos[135] + armZZos[25] - 
   1861./864. + 15*armZZos[28];
   armZZos[96]=armZZos[1]*armZZos[96];
   armZZos[107]= - 9911./96. + 685*armZZos[28];
   armZZos[107]=121./9.*armZZos[148] + 22*armZZos[191] - 77./3.*
   armZZos[10] - 55./36.*armZZos[9] - 11*armZZos[13] - 685./36.*
   armZZos[137] + 253./24.*armZZos[11] + 1./9.*armZZos[107] + 11*
   armZZos[25];
   armZZos[107]=armZZos[34]*armZZos[107];
   armZZos[113]=armZZos[127] + armZZos[125] + 3*armZZos[8] + 
   armZZos[124];
   armZZos[113]=armZZos[1]*armZZos[113];
   armZZos[122]=armZZos[202] - armZZos[18] + armZZos[8] + armZZos[200];
   armZZos[122]=armZZos[34]*armZZos[122];
   armZZos[124]=armZZos[126] + 11./9.*armZZos[131];
   armZZos[125]=3*armZZos[1];
   armZZos[126]=armZZos[125] + 11./3.*armZZos[34];
   armZZos[126]=armZZos[15]*armZZos[126];
   armZZos[126]=armZZos[124] + armZZos[126];
   armZZos[126]=armZZos[6]*armZZos[126];
   armZZos[127]= - 9*armZZos[1] - 11*armZZos[34];
   armZZos[127]=armZZos[15]*armZZos[127];
   armZZos[113]=1./2.*armZZos[126] + 1./2.*armZZos[127] + armZZos[113]
    + 11./3.*armZZos[122];
   armZZos[113]=armZZos[5]*armZZos[113];
   armZZos[122]=armZZos[15]*armZZos[177];
   armZZos[126]=1./2.*armZZos[124] + armZZos[122];
   armZZos[126]=armZZos[6]*armZZos[126];
   armZZos[127]=armZZos[15]*armZZos[196];
   armZZos[102]=armZZos[126] + 2./3.*armZZos[127] + armZZos[105] + 11./
   9.*armZZos[102];
   armZZos[102]=armZZos[89]*armZZos[102];
   armZZos[105]=armZZos[34]*armZZos[186];
   armZZos[97]=armZZos[1]*armZZos[97];
   armZZos[98]=armZZos[211] + 11./9.*armZZos[98];
   armZZos[98]=armZZos[6]*armZZos[98];
   armZZos[97]=armZZos[98] + armZZos[97] + 11./3.*armZZos[105];
   armZZos[97]=armZZos[89]*armZZos[97];
   armZZos[98]= - armZZos[1]*armZZos[22];
   armZZos[105]= - armZZos[34]*armZZos[22];
   armZZos[98]=3*armZZos[98] + 137./81.*armZZos[105];
   armZZos[105]=7./6.*armZZos[177] + armZZos[185];
   armZZos[105]=armZZos[5]*armZZos[105];
   armZZos[97]=armZZos[97] + 2*armZZos[98] + armZZos[105];
   armZZos[97]=mmZ*armZZos[97];
   armZZos[98]=13./3. + armZZos[9];
   armZZos[98]=1./2.*armZZos[98] - 5*armZZos[1];
   armZZos[98]=armZZos[1]*armZZos[98];
   armZZos[105]=7 + 11./3.*armZZos[9];
   armZZos[105]= - 605./27.*armZZos[34] + 1./2.*armZZos[105] - 110./3.*
   armZZos[1];
   armZZos[105]=armZZos[34]*armZZos[105];
   armZZos[98]=armZZos[98] + 1./3.*armZZos[105];
   armZZos[105]=15./4. + armZZos[1];
   armZZos[105]=armZZos[1]*armZZos[105];
   armZZos[126]=121./9.*armZZos[34] + 685./36. + 22*armZZos[1];
   armZZos[126]=armZZos[34]*armZZos[126];
   armZZos[105]=armZZos[105] + 1./9.*armZZos[126];
   armZZos[105]=armZZos[6]*armZZos[105];
   armZZos[98]=1./3.*armZZos[98] + armZZos[105];
   armZZos[98]=armZZos[6]*armZZos[98];
   armZZos[105]= - 1./3.*armZZos[77];
   armZZos[126]=1./3.*armZZos[12] + 1./4. + armZZos[105];
   armZZos[126]=armZZos[100]*armZZos[62]*armZZos[126];
   armZZos[94]=1./8.*armZZos[126] + armZZos[94] + armZZos[97] + 
   armZZos[102] + armZZos[113] + armZZos[98] + armZZos[96] + 1./9.*
   armZZos[107];
   armZZos[94]=armZZos[100]*armZZos[94];
   armZZos[96]=167./108. - armZZos[20];
   armZZos[96]=121./27.*armZZos[106] + 22./3.*armZZos[117] - 65./9.*
   armZZos[9] + 329./27.*armZZos[137] - 1370./27.*armZZos[28] + 1./6.*
   armZZos[96] + armZZos[29];
   armZZos[96]=armZZos[34]*armZZos[96];
   armZZos[97]= - 37./12. - armZZos[20];
   armZZos[97]=armZZos[117] - 17./9.*armZZos[9] + 7*armZZos[137] - 30*
   armZZos[28] + 11./18.*armZZos[97] + armZZos[29];
   armZZos[97]=armZZos[1]*armZZos[97];
   armZZos[98]=armZZos[153] - 53./9. + armZZos[99];
   armZZos[98]=armZZos[1]*armZZos[98];
   armZZos[102]= - 121./9.*armZZos[34] - 329./9. - 22*armZZos[1];
   armZZos[102]=armZZos[34]*armZZos[102];
   armZZos[102]=armZZos[206] + 1./9.*armZZos[102];
   armZZos[102]=armZZos[6]*armZZos[102];
   armZZos[106]=968./81.*armZZos[34] + 176./9.*armZZos[1] - 377./27. + 
   armZZos[142];
   armZZos[106]=armZZos[34]*armZZos[106];
   armZZos[98]=armZZos[102] + armZZos[98] + 1./3.*armZZos[106];
   armZZos[98]=armZZos[6]*armZZos[98];
   armZZos[102]=armZZos[199] + 11./9.*armZZos[203];
   armZZos[106]=3*armZZos[176] + 11./3.*armZZos[184];
   armZZos[106]=armZZos[6]*armZZos[106];
   armZZos[102]=4*armZZos[102] + armZZos[106];
   armZZos[102]=armZZos[5]*armZZos[102];
   armZZos[106]=armZZos[147] - armZZos[23];
   armZZos[106]=1./3.*armZZos[106] + 12*armZZos[22];
   armZZos[106]=armZZos[1]*armZZos[106];
   armZZos[107]=548./9.*armZZos[22] - 73./9.*armZZos[24] - armZZos[23];
   armZZos[107]=armZZos[34]*armZZos[107];
   armZZos[106]=armZZos[106] + 1./9.*armZZos[107];
   armZZos[106]=mmZ*armZZos[106];
   armZZos[93]=armZZos[94] + armZZos[93] + armZZos[106] + armZZos[102]
    + armZZos[98] + armZZos[97] + 1./3.*armZZos[96];
   armZZos[93]=armZZos[100]*armZZos[93];
   armZZos[94]=3*armZZos[68];
   armZZos[96]=armZZos[94] + armZZos[69];
   armZZos[97]=armZZos[96] + 7./12.*armZZos[70];
   armZZos[98]=1./3.*armZZos[67];
   armZZos[97]=1./2.*armZZos[97] + armZZos[98];
   armZZos[97]=mmH*armZZos[97];
   armZZos[102]= - 1./4. - armZZos[9];
   armZZos[102]=armZZos[9]*armZZos[102];
   armZZos[106]= - 47./9. + armZZos[195];
   armZZos[106]=1./2.*armZZos[106] + armZZos[10];
   armZZos[106]=armZZos[10]*armZZos[106];
   armZZos[107]=5./3.*armZZos[76];
   armZZos[113]=1./6.*armZZos[85];
   armZZos[117]= - 1./3.*armZZos[73];
   armZZos[97]=armZZos[97] + 1./2.*armZZos[106] + 1./6.*armZZos[102] - 
   31./8.*armZZos[13] - 13./12.*armZZos[12] + 5./12.*armZZos[78] + 
   armZZos[113] + armZZos[117] + armZZos[219] + armZZos[107] + 1./24.*
   armZZos[79] + 17./12.*armZZos[86] + 1./8.*armZZos[74] + 143./32. + 
   armZZos[105];
   armZZos[97]=mmH*armZZos[97];
   armZZos[102]=armZZos[115] - armZZos[16];
   armZZos[102]=armZZos[16]*armZZos[102];
   armZZos[105]= - 5./4.*armZZos[15];
   armZZos[106]=armZZos[105] + 13./2.*armZZos[14] + armZZos[16];
   armZZos[106]=armZZos[15]*armZZos[106];
   armZZos[126]=armZZos[106] + 33./4.*armZZos[193] + armZZos[102];
   armZZos[126]=armZZos[5]*armZZos[126];
   armZZos[127]=5*armZZos[16];
   armZZos[129]= - 4./3.*armZZos[14] + armZZos[127];
   armZZos[129]=2*armZZos[129] + armZZos[134];
   armZZos[129]=armZZos[15]*armZZos[129];
   armZZos[131]=armZZos[245] + 23./2.*armZZos[16];
   armZZos[131]=armZZos[16]*armZZos[131];
   armZZos[129]=armZZos[129] - 2*armZZos[193] + 1./3.*armZZos[131];
   armZZos[129]=armZZos[89]*armZZos[129];
   armZZos[131]= - 4955./3. + 613./2.*armZZos[9];
   armZZos[131]=1./3.*armZZos[131] + 11*armZZos[173];
   armZZos[131]=armZZos[14]*armZZos[131];
   armZZos[133]=1./8.*armZZos[9];
   armZZos[134]= - 2./3. + armZZos[133];
   armZZos[134]=43*armZZos[134] + armZZos[179];
   armZZos[134]=armZZos[16]*armZZos[134];
   armZZos[135]= - 5*armZZos[9];
   armZZos[139]=7./4. + armZZos[135];
   armZZos[139]=7./3.*armZZos[139] + armZZos[146];
   armZZos[139]=1./2.*armZZos[139] + armZZos[10];
   armZZos[139]=armZZos[15]*armZZos[139];
   armZZos[142]=9./4.*armZZos[66];
   armZZos[144]=11./12.*armZZos[80];
   armZZos[147]=1./6.*armZZos[65];
   armZZos[148]=3./2.*armZZos[16];
   armZZos[150]=23./3.*armZZos[14] + armZZos[148];
   armZZos[150]=armZZos[10]*armZZos[150];
   armZZos[97]=armZZos[97] + armZZos[129] + armZZos[126] + armZZos[139]
    + armZZos[150] + 1./3.*armZZos[134] + 1./12.*armZZos[131] - 33./4.*
   armZZos[18] - 511./72.*armZZos[19] - 67./12.*armZZos[17] + 
   armZZos[147] + 151./12.*armZZos[82] - 2./3.*armZZos[64] + 
   armZZos[144] - 67./12.*armZZos[63] - 41./72.*armZZos[81] - 10./9.*
   armZZos[83] + armZZos[142];
   armZZos[97]=armZZos[62]*armZZos[97];
   armZZos[101]=armZZos[101] + armZZos[94] + armZZos[109];
   armZZos[101]=1./2.*armZZos[101] + armZZos[98];
   armZZos[101]=mmH*armZZos[101];
   armZZos[95]=3089./108. + armZZos[95];
   armZZos[95]=1./3.*armZZos[79] + 1./2.*armZZos[95] + armZZos[86];
   armZZos[109]= - 361./3. + armZZos[174];
   armZZos[109]=1./18.*armZZos[109] + armZZos[10];
   armZZos[109]=armZZos[10]*armZZos[109];
   armZZos[95]=armZZos[101] + 1./2.*armZZos[109] + 5./72.*armZZos[9] - 
   35./8.*armZZos[13] + armZZos[239] + 1./36.*armZZos[78] + 7./6.*
   armZZos[85] + armZZos[117] + armZZos[219] + 1./4.*armZZos[95] + 
   armZZos[107];
   armZZos[95]=mmH*armZZos[95];
   armZZos[101]=armZZos[105] + armZZos[136] + armZZos[16];
   armZZos[101]=armZZos[15]*armZZos[101];
   armZZos[109]= - 1./2.*armZZos[193];
   armZZos[126]=pow(armZZos[16],2);
   armZZos[101]=armZZos[101] + armZZos[109] - armZZos[126];
   armZZos[101]=armZZos[5]*armZZos[101];
   armZZos[131]= - 1./8.*armZZos[15];
   armZZos[134]=armZZos[131] - 8./3.*armZZos[14] + armZZos[127];
   armZZos[134]=armZZos[15]*armZZos[134];
   armZZos[136]=armZZos[14] + 23./4.*armZZos[16];
   armZZos[136]=armZZos[16]*armZZos[136];
   armZZos[134]=armZZos[134] - armZZos[193] + 1./3.*armZZos[136];
   armZZos[134]=armZZos[89]*armZZos[134];
   armZZos[136]=11./3.*armZZos[80];
   armZZos[139]=armZZos[136] - 15./2.*armZZos[63] + armZZos[116] - 1./
   18.*armZZos[81];
   armZZos[139]=1./2.*armZZos[139] + armZZos[64];
   armZZos[139]=armZZos[147] + 1./2.*armZZos[139] + 7./3.*armZZos[82];
   armZZos[152]= - 1285./3. + 83*armZZos[9];
   armZZos[152]=armZZos[14]*armZZos[152];
   armZZos[133]= - 9 + armZZos[133];
   armZZos[133]=armZZos[16]*armZZos[133];
   armZZos[153]=3./4.*armZZos[16];
   armZZos[154]=5./3.*armZZos[14] + armZZos[153];
   armZZos[154]=armZZos[10]*armZZos[154];
   armZZos[156]=61./3. + armZZos[220];
   armZZos[156]=1./18.*armZZos[156] + armZZos[10];
   armZZos[156]=armZZos[15]*armZZos[156];
   armZZos[95]=1./2.*armZZos[95] + armZZos[134] + 1./2.*armZZos[101] + 
   1./2.*armZZos[156] + armZZos[154] + 1./2.*armZZos[133] + 1./72.*
   armZZos[152] - 41./48.*armZZos[18] - 221./48.*armZZos[19] + 1./2.*
   armZZos[139] - 2./3.*armZZos[17];
   armZZos[95]=armZZos[62]*armZZos[95];
   armZZos[101]=4./3.*armZZos[14] + armZZos[153];
   armZZos[101]=armZZos[1]*armZZos[101];
   armZZos[133]=20./9.*armZZos[14] + 11./4.*armZZos[16];
   armZZos[133]=armZZos[34]*armZZos[133];
   armZZos[134]=91*armZZos[1] + 1403./27.*armZZos[34];
   armZZos[134]=armZZos[15]*armZZos[134];
   armZZos[101]=1./12.*armZZos[134] + armZZos[101] + 1./3.*armZZos[133]
   ;
   armZZos[101]=armZZos[6]*armZZos[101];
   armZZos[133]=armZZos[1]*armZZos[141];
   armZZos[134]=armZZos[34]*armZZos[141];
   armZZos[139]=armZZos[133] + 11./9.*armZZos[134];
   armZZos[139]=armZZos[6]*armZZos[139];
   armZZos[141]=13./18.*armZZos[10] + 1./4.*armZZos[13] - 1./24.*
   armZZos[11] - 1./4.*armZZos[25] + 4./9. + 1./4.*armZZos[27];
   armZZos[152]=armZZos[1]*armZZos[141];
   armZZos[141]=armZZos[34]*armZZos[141];
   armZZos[139]=1./12.*armZZos[139] + armZZos[152] + 11./9.*
   armZZos[141];
   armZZos[139]=mmH*armZZos[139];
   armZZos[154]= - 1./2.*armZZos[18] - armZZos[17] - 1./2.*armZZos[63]
    + armZZos[82];
   armZZos[156]= - 1 + 1./6.*armZZos[9];
   armZZos[156]=armZZos[14]*armZZos[156];
   armZZos[160]= - 1 - armZZos[9];
   armZZos[160]=armZZos[15]*armZZos[160];
   armZZos[154]=1./6.*armZZos[160] + 1./3.*armZZos[154] + armZZos[156];
   armZZos[154]=armZZos[100]*armZZos[62]*armZZos[154];
   armZZos[156]=1./4.*armZZos[19];
   armZZos[160]= - 23./24.*armZZos[16];
   armZZos[165]=armZZos[160] - 23./18.*armZZos[14] + 29./4.*armZZos[18]
    + armZZos[156] + 31./4.*armZZos[8] + armZZos[17];
   armZZos[165]=armZZos[1]*armZZos[165];
   armZZos[166]= - 253./72.*armZZos[16] - 109./54.*armZZos[14] + 1271./
   108.*armZZos[18] + 11./12.*armZZos[19] + 1469./108.*armZZos[8] + 
   armZZos[17];
   armZZos[166]=armZZos[34]*armZZos[166];
   armZZos[167]= - 359*armZZos[1] - 623./3.*armZZos[34];
   armZZos[167]=armZZos[15]*armZZos[167];
   armZZos[95]=1./4.*armZZos[154] + armZZos[95] + armZZos[139] + 
   armZZos[101] + 1./18.*armZZos[167] + armZZos[165] + 1./3.*
   armZZos[166];
   armZZos[95]=armZZos[100]*armZZos[95];
   armZZos[101]= - 637./18.*armZZos[14] - 685./9.*armZZos[18] + 
   armZZos[108] + armZZos[198] - 685./9.*armZZos[8];
   armZZos[101]=armZZos[34]*armZZos[101];
   armZZos[139]=armZZos[176] + 29./27.*armZZos[184];
   armZZos[154]=armZZos[234] - 137./81.*armZZos[34];
   armZZos[154]=armZZos[15]*armZZos[154];
   armZZos[139]=4*armZZos[139] + 5*armZZos[154];
   armZZos[139]=armZZos[6]*armZZos[139];
   armZZos[154]= - 29./6.*armZZos[14] - 15*armZZos[18] + 2./3.*
   armZZos[17] + 2./3.*armZZos[7] - 15*armZZos[8];
   armZZos[154]=armZZos[1]*armZZos[154];
   armZZos[125]=armZZos[125] + 137./81.*armZZos[34];
   armZZos[125]=armZZos[15]*armZZos[125];
   armZZos[95]=armZZos[95] + armZZos[97] + armZZos[139] + 13*
   armZZos[125] + armZZos[154] + 1./9.*armZZos[101];
   armZZos[95]=armZZos[100]*armZZos[95];
   armZZos[96]=armZZos[96] - 17./12.*armZZos[70];
   armZZos[96]=1./2.*armZZos[96] + armZZos[98];
   armZZos[96]=mmH*armZZos[96];
   armZZos[97]=223./96. - armZZos[77];
   armZZos[101]=7./12. + armZZos[9];
   armZZos[101]=armZZos[9]*armZZos[101];
   armZZos[125]=73./3. + 71./2.*armZZos[9];
   armZZos[125]=1./6.*armZZos[125] + armZZos[10];
   armZZos[125]=armZZos[10]*armZZos[125];
   armZZos[96]=armZZos[96] + 1./2.*armZZos[125] + 1./2.*armZZos[101] - 
   39./8.*armZZos[13] + 1./4.*armZZos[12] + 197./12.*armZZos[78] + 
   armZZos[113] + armZZos[117] + armZZos[219] + armZZos[107] - 89./8.*
   armZZos[79] + 1./12.*armZZos[86] + 1./3.*armZZos[97] + 9./8.*
   armZZos[74];
   armZZos[96]=mmH*armZZos[96];
   armZZos[97]=armZZos[106] + 81./4.*armZZos[193] + armZZos[102];
   armZZos[97]=armZZos[5]*armZZos[97];
   armZZos[101]=28531./12. + 677*armZZos[9];
   armZZos[101]=1./3.*armZZos[101] + armZZos[146];
   armZZos[101]=1./2.*armZZos[101] + armZZos[10];
   armZZos[101]=armZZos[15]*armZZos[101];
   armZZos[102]=59*armZZos[173] - 78599./9. + 3487./2.*armZZos[9];
   armZZos[102]=armZZos[14]*armZZos[102];
   armZZos[106]=25./6.*armZZos[178] - 107./3. + 49./8.*armZZos[9];
   armZZos[106]=armZZos[16]*armZZos[106];
   armZZos[96]=armZZos[96] + armZZos[129] + armZZos[97] + armZZos[101]
    + armZZos[150] + armZZos[106] + 1./12.*armZZos[102] - 3467./12.*
   armZZos[18] - 101./24.*armZZos[19] - 7099./12.*armZZos[17] + 
   armZZos[147] + 2087./12.*armZZos[82] - 40./3.*armZZos[64] + 
   armZZos[144] + 455./4.*armZZos[63] + 605./24.*armZZos[81] - 2*
   armZZos[83] + armZZos[142];
   armZZos[96]=armZZos[62]*armZZos[96];
   armZZos[94]=11./12.*armZZos[70] + armZZos[94] + 11./3.*armZZos[69];
   armZZos[94]=1./2.*armZZos[94] + armZZos[98];
   armZZos[94]=mmH*armZZos[94];
   armZZos[97]=5./2.*armZZos[76] - 29./8.*armZZos[86] - 29./16.*
   armZZos[74] + 205./576. - 2*armZZos[77];
   armZZos[98]= - 7./8. + armZZos[162];
   armZZos[98]=armZZos[9]*armZZos[98];
   armZZos[101]= - 1247./9. - 161./2.*armZZos[9];
   armZZos[101]=1./6.*armZZos[101] + armZZos[10];
   armZZos[101]=armZZos[10]*armZZos[101];
   armZZos[102]= - 1./6.*armZZos[73];
   armZZos[94]=1./2.*armZZos[94] + 1./4.*armZZos[101] + 1./2.*
   armZZos[98] - 137./48.*armZZos[13] + 15./8.*armZZos[12] - 95./8.*
   armZZos[78] + 5./3.*armZZos[85] + armZZos[102] + 1./3.*armZZos[97]
    + 9./8.*armZZos[75];
   armZZos[94]=mmH*armZZos[94];
   armZZos[97]=armZZos[16]*armZZos[121];
   armZZos[98]=armZZos[105] - 49./2.*armZZos[14] + armZZos[16];
   armZZos[98]=armZZos[15]*armZZos[98];
   armZZos[97]=1./2.*armZZos[98] - 23./2.*armZZos[193] + armZZos[97];
   armZZos[97]=armZZos[5]*armZZos[97];
   armZZos[98]=armZZos[131] - 8*armZZos[14] + armZZos[127];
   armZZos[98]=armZZos[15]*armZZos[98];
   armZZos[101]=armZZos[115] + 23./12.*armZZos[16];
   armZZos[101]=armZZos[16]*armZZos[101];
   armZZos[98]=armZZos[98] - armZZos[193] + armZZos[101];
   armZZos[98]=armZZos[89]*armZZos[98];
   armZZos[101]=armZZos[136] - 489./2.*armZZos[63] + armZZos[116] - 197.
   /2.*armZZos[81];
   armZZos[101]=1./4.*armZZos[101] + 11*armZZos[64];
   armZZos[105]=13*armZZos[178] + 45421./18. + armZZos[201];
   armZZos[105]=armZZos[14]*armZZos[105];
   armZZos[106]=517./9. - 147./4.*armZZos[9];
   armZZos[106]=1./4.*armZZos[106] + armZZos[146];
   armZZos[106]=armZZos[16]*armZZos[106];
   armZZos[107]= - 36203./27. - 467./2.*armZZos[9];
   armZZos[107]=armZZos[10] + 1./4.*armZZos[107] + armZZos[178];
   armZZos[107]=armZZos[15]*armZZos[107];
   armZZos[113]=29*armZZos[14] + armZZos[148];
   armZZos[113]=armZZos[10]*armZZos[113];
   armZZos[115]=1./12.*armZZos[65];
   armZZos[94]=armZZos[94] + armZZos[98] + armZZos[97] + 1./2.*
   armZZos[107] + 1./2.*armZZos[113] + armZZos[106] + 1./4.*
   armZZos[105] + 3023./16.*armZZos[18] - 151./16.*armZZos[19] + 1571./
   4.*armZZos[17] + armZZos[115] + 1./2.*armZZos[101] - 159*armZZos[82]
   ;
   armZZos[94]=armZZos[62]*armZZos[94];
   armZZos[97]=8./3.*armZZos[14] + armZZos[214];
   armZZos[97]=armZZos[1]*armZZos[97];
   armZZos[98]=8*armZZos[14] + armZZos[153];
   armZZos[98]=armZZos[34]*armZZos[98];
   armZZos[97]=31./12.*armZZos[204] + armZZos[97] + armZZos[98];
   armZZos[97]=armZZos[6]*armZZos[97];
   armZZos[98]=1./3.*armZZos[133] + armZZos[134];
   armZZos[98]=armZZos[6]*armZZos[98];
   armZZos[98]=1./12.*armZZos[98] + 1./3.*armZZos[152] + armZZos[141];
   armZZos[98]=mmH*armZZos[98];
   armZZos[101]= - armZZos[15]*armZZos[14];
   armZZos[105]=armZZos[109] + armZZos[101];
   armZZos[105]=armZZos[89]*armZZos[105];
   armZZos[106]=armZZos[193]*armZZos[90];
   armZZos[107]=armZZos[15]*armZZos[173];
   armZZos[105]=armZZos[105] + 1./2.*armZZos[106] + armZZos[107];
   armZZos[105]=armZZos[112]*armZZos[62]*armZZos[105];
   armZZos[104]=armZZos[7] + armZZos[104];
   armZZos[104]=armZZos[156] + 11./2.*armZZos[104] + 13*armZZos[17];
   armZZos[106]= - 23./72.*armZZos[16] - 163./36.*armZZos[14] + 1./3.*
   armZZos[104] + 3./4.*armZZos[18];
   armZZos[106]=armZZos[1]*armZZos[106];
   armZZos[104]=armZZos[160] - 163./12.*armZZos[14] + armZZos[104] + 9./
   4.*armZZos[18];
   armZZos[104]=armZZos[34]*armZZos[104];
   armZZos[94]=3./2.*armZZos[105] + armZZos[94] + armZZos[98] + 
   armZZos[97] + 125./18.*armZZos[175] + armZZos[106] + armZZos[104];
   armZZos[94]=armZZos[112]*armZZos[94];
   armZZos[97]=1./3.*armZZos[199] + armZZos[203];
   armZZos[97]=28./3.*armZZos[97] + 5*armZZos[175];
   armZZos[97]=armZZos[6]*armZZos[97];
   armZZos[98]=239./18.*armZZos[14] + armZZos[114] - 8*armZZos[17] - 14
   *armZZos[7] - 5*armZZos[8];
   armZZos[104]=armZZos[1]*armZZos[98];
   armZZos[98]=armZZos[34]*armZZos[98];
   armZZos[94]=armZZos[94] + armZZos[96] + armZZos[97] + 13*
   armZZos[204] + 1./3.*armZZos[104] + armZZos[98];
   armZZos[94]=armZZos[112]*armZZos[94];
   armZZos[96]=619./9. + 59*armZZos[9];
   armZZos[97]= - 11./3.*armZZos[10];
   armZZos[96]=1./4.*armZZos[96] + armZZos[97];
   armZZos[96]=armZZos[10]*armZZos[96];
   armZZos[98]=1./4.*armZZos[76] + 13./2.*armZZos[86] + 13./4.*
   armZZos[74] + 9877./384. + armZZos[77];
   armZZos[98]=1./3.*armZZos[98] - 7./8.*armZZos[75];
   armZZos[104]=pow(armZZos[9],2);
   armZZos[105]=1./24.*armZZos[104];
   armZZos[106]= - 1./4.*armZZos[70];
   armZZos[107]= - armZZos[67] - armZZos[69] + armZZos[106];
   armZZos[107]=mmH*armZZos[107];
   armZZos[96]=1./6.*armZZos[107] + 1./12.*armZZos[96] + armZZos[105]
    + 17./48.*armZZos[13] - 5./4.*armZZos[12] + 85./12.*armZZos[78] - 1.
   /6.*armZZos[85] + 1./2.*armZZos[98] + armZZos[117];
   armZZos[96]=mmH*armZZos[96];
   armZZos[98]= - 5*armZZos[63] + armZZos[169] + 169./3.*armZZos[81];
   armZZos[98]= - 25./3.*armZZos[64] + 1./2.*armZZos[98] + armZZos[80];
   armZZos[107]=1./3.*armZZos[82];
   armZZos[109]= - 227./6. + armZZos[135];
   armZZos[109]=armZZos[14]*armZZos[109];
   armZZos[113]= - 4505./27. + 41*armZZos[9];
   armZZos[113]=armZZos[16]*armZZos[113];
   armZZos[114]= - 21*armZZos[14] - 5./2.*armZZos[16];
   armZZos[114]=armZZos[10]*armZZos[114];
   armZZos[116]= - 167./27. - 21./2.*armZZos[9];
   armZZos[116]=1./2.*armZZos[116] + 11./9.*armZZos[10];
   armZZos[116]=armZZos[15]*armZZos[116];
   armZZos[96]=armZZos[96] + 1./4.*armZZos[116] + 1./4.*armZZos[114] + 
   1./16.*armZZos[113] + 1./4.*armZZos[109] + armZZos[145] + 17./48.*
   armZZos[19] - 28./3.*armZZos[17] + 7./48.*armZZos[65] + 1./4.*
   armZZos[98] + armZZos[107];
   armZZos[96]=mmH*armZZos[96];
   armZZos[98]=9*armZZos[16];
   armZZos[109]= - 5165./3.*armZZos[14] + armZZos[98];
   armZZos[113]=11./9.*armZZos[15];
   armZZos[109]=1./4.*armZZos[109] + armZZos[113];
   armZZos[109]=armZZos[15]*armZZos[109];
   armZZos[114]= - armZZos[14] + 13./12.*armZZos[16];
   armZZos[114]=armZZos[16]*armZZos[114];
   armZZos[109]=armZZos[109] - 1109./4.*armZZos[193] + armZZos[114];
   armZZos[96]=1./2.*armZZos[109] + armZZos[96];
   armZZos[96]=armZZos[62]*armZZos[96];
   armZZos[109]= - armZZos[1]*armZZos[10];
   armZZos[114]= - armZZos[34]*armZZos[10];
   armZZos[116]=1./3.*armZZos[109] + armZZos[114];
   armZZos[116]=armZZos[6]*armZZos[116];
   armZZos[121]=armZZos[11] + armZZos[183];
   armZZos[125]=armZZos[1]*armZZos[121];
   armZZos[121]=armZZos[34]*armZZos[121];
   armZZos[116]=armZZos[116] + 1./3.*armZZos[125] + armZZos[121];
   armZZos[116]=mmH*armZZos[116];
   armZZos[129]=armZZos[192] + 5./2.*armZZos[204];
   armZZos[131]=armZZos[192] + armZZos[204];
   armZZos[131]=armZZos[6]*armZZos[131];
   armZZos[116]=1./2.*armZZos[116] + 1./3.*armZZos[129] + 1./2.*
   armZZos[131];
   armZZos[116]=mmH*armZZos[116];
   armZZos[96]=1./6.*armZZos[116] + armZZos[96];
   armZZos[96]=armZZos[112]*armZZos[96];
   armZZos[116]= - 311./192. + armZZos[77];
   armZZos[129]=5./4.*armZZos[86];
   armZZos[131]=1./6.*armZZos[76];
   armZZos[133]= - 7./4.*armZZos[75];
   armZZos[116]=armZZos[133] + armZZos[131] + 43./12.*armZZos[79] + 
   armZZos[129] + 1./3.*armZZos[116] - 5./4.*armZZos[74];
   armZZos[134]= - 41./3. - 13*armZZos[9];
   armZZos[97]=1./2.*armZZos[134] + armZZos[97];
   armZZos[97]=armZZos[10]*armZZos[97];
   armZZos[134]= - 2./3.*armZZos[73];
   armZZos[135]= - 7./24.*armZZos[85];
   armZZos[136]= - 19./24.*armZZos[12];
   armZZos[104]= - 1./12.*armZZos[104];
   armZZos[139]= - armZZos[69] + 1./2.*armZZos[70];
   armZZos[139]=1./2.*armZZos[139] - armZZos[67];
   armZZos[139]=1./3.*mmH*armZZos[139];
   armZZos[97]=armZZos[139] + 1./6.*armZZos[97] + armZZos[104] + 19./8.
   *armZZos[13] + armZZos[136] - 89./24.*armZZos[78] + armZZos[135] + 1.
   /2.*armZZos[116] + armZZos[134];
   armZZos[97]=mmH*armZZos[97];
   armZZos[116]=11./6.*armZZos[63] + armZZos[169] - 89./3.*armZZos[81];
   armZZos[141]=7./12.*armZZos[65];
   armZZos[142]= - 19./4.*armZZos[19];
   armZZos[144]=1319./3. + armZZos[213];
   armZZos[144]=armZZos[14]*armZZos[144];
   armZZos[145]=215./3. + armZZos[190];
   armZZos[145]=armZZos[16]*armZZos[145];
   armZZos[146]= - 91./9.*armZZos[14] - 5*armZZos[16];
   armZZos[146]=1./2.*armZZos[10]*armZZos[146];
   armZZos[116]=armZZos[146] + 1./6.*armZZos[145] + 1./12.*armZZos[144]
    - 49./12.*armZZos[18] + armZZos[142] - 79./12.*armZZos[17] + 
   armZZos[141] + armZZos[82] + 61./12.*armZZos[64] + 1./2.*
   armZZos[116] + armZZos[80];
   armZZos[144]=11./6.*armZZos[10] + 7 + 11./4.*armZZos[9];
   armZZos[144]=armZZos[15]*armZZos[144];
   armZZos[97]=armZZos[97] + 1./2.*armZZos[116] + 1./3.*armZZos[144];
   armZZos[97]=mmH*armZZos[97];
   armZZos[116]=armZZos[149] + armZZos[217];
   armZZos[116]=armZZos[16]*armZZos[116];
   armZZos[144]=9791./6.*armZZos[193] + armZZos[116];
   armZZos[145]=13453./18.*armZZos[14] + armZZos[98];
   armZZos[145]=1./4.*armZZos[145] + armZZos[113];
   armZZos[145]=armZZos[15]*armZZos[145];
   armZZos[97]=armZZos[97] + 1./12.*armZZos[144] + armZZos[145];
   armZZos[97]=armZZos[62]*armZZos[97];
   armZZos[96]=armZZos[97] + armZZos[96];
   armZZos[96]=armZZos[112]*armZZos[96];
   armZZos[97]=2185./192. + armZZos[77];
   armZZos[97]=armZZos[133] + armZZos[131] + 13./4.*armZZos[79] + 
   armZZos[129] + 1./3.*armZZos[97] - 1./4.*armZZos[74];
   armZZos[129]=7./9. - armZZos[9];
   armZZos[129]=1./2.*armZZos[129] - 11./9.*armZZos[10];
   armZZos[129]=armZZos[10]*armZZos[129];
   armZZos[97]=armZZos[139] + 1./2.*armZZos[129] + armZZos[104] + 15./8.
   *armZZos[13] + armZZos[136] - 3./8.*armZZos[78] + armZZos[135] + 1./
   2.*armZZos[97] + armZZos[134];
   armZZos[97]=mmH*armZZos[97];
   armZZos[104]= - 1./6.*armZZos[63] + armZZos[169] + armZZos[118];
   armZZos[118]=575./3. + armZZos[213];
   armZZos[118]=armZZos[14]*armZZos[118];
   armZZos[99]=41./9. + armZZos[99];
   armZZos[99]=armZZos[16]*armZZos[99];
   armZZos[99]=armZZos[146] + 1./2.*armZZos[99] + 1./12.*armZZos[118]
    - 37./12.*armZZos[18] + armZZos[142] - 103./12.*armZZos[17] + 
   armZZos[141] + armZZos[82] + 29./12.*armZZos[64] + 1./2.*
   armZZos[104] + armZZos[80];
   armZZos[104]=11./18.*armZZos[10] + 1 + armZZos[110];
   armZZos[104]=armZZos[15]*armZZos[104];
   armZZos[97]=armZZos[97] + 1./2.*armZZos[99] + armZZos[104];
   armZZos[97]=mmH*armZZos[97];
   armZZos[99]=1151./6.*armZZos[193] + armZZos[116];
   armZZos[98]=1357./18.*armZZos[14] + armZZos[98];
   armZZos[98]=1./4.*armZZos[98] + armZZos[113];
   armZZos[98]=armZZos[15]*armZZos[98];
   armZZos[97]=armZZos[97] + 1./12.*armZZos[99] + armZZos[98];
   armZZos[97]=armZZos[62]*armZZos[97];
   armZZos[98]= - 7./2.*armZZos[75] + 1./3.*armZZos[76] + armZZos[79]
    + 2621./288. - 5*armZZos[86];
   armZZos[99]=65./3. - armZZos[9];
   armZZos[99]=1./4.*armZZos[99] - 11*armZZos[10];
   armZZos[99]=armZZos[10]*armZZos[99];
   armZZos[104]=armZZos[106] - armZZos[67];
   armZZos[104]=mmH*armZZos[104];
   armZZos[98]=1./6.*armZZos[104] + 1./36.*armZZos[99] + armZZos[105]
    + 41./48.*armZZos[13] + 5./8.*armZZos[12] - 37./72.*armZZos[78] - 1.
   /8.*armZZos[85] + 1./8.*armZZos[98] + armZZos[117];
   armZZos[98]=mmH*armZZos[98];
   armZZos[99]=armZZos[169] - 13./9.*armZZos[81];
   armZZos[99]=7./6.*armZZos[64] + armZZos[80] + 1./2.*armZZos[99] - 1./
   3.*armZZos[63];
   armZZos[99]=armZZos[164] - 115./24.*armZZos[19] - 25./36.*
   armZZos[17] + 7./24.*armZZos[65] + 1./2.*armZZos[99] + armZZos[107];
   armZZos[104]=83./3. + armZZos[190];
   armZZos[104]=1./2.*armZZos[104] + armZZos[230];
   armZZos[104]=armZZos[15]*armZZos[104];
   armZZos[105]= - 1./9. + 5./8.*armZZos[9];
   armZZos[105]=armZZos[14]*armZZos[105];
   armZZos[106]= - 187./3. - 79*armZZos[9];
   armZZos[106]=armZZos[16]*armZZos[106];
   armZZos[107]= - 11./9.*armZZos[14] - 5./4.*armZZos[16];
   armZZos[107]=armZZos[10]*armZZos[107];
   armZZos[98]=armZZos[98] + 1./36.*armZZos[104] + 1./2.*armZZos[107]
    + 1./144.*armZZos[106] + 1./2.*armZZos[99] + armZZos[105];
   armZZos[98]=mmH*armZZos[98];
   armZZos[99]=17./3.*armZZos[14];
   armZZos[104]=armZZos[99] + 13./2.*armZZos[16];
   armZZos[104]=armZZos[16]*armZZos[104];
   armZZos[104]=49./3.*armZZos[193] + 1./2.*armZZos[104];
   armZZos[105]=armZZos[113] + 49./9.*armZZos[14] + armZZos[161];
   armZZos[105]=armZZos[15]*armZZos[105];
   armZZos[104]=1./3.*armZZos[104] + armZZos[105];
   armZZos[98]=1./2.*armZZos[104] + armZZos[98];
   armZZos[98]=armZZos[62]*armZZos[98];
   armZZos[104]=armZZos[109] + 11./9.*armZZos[114];
   armZZos[104]=armZZos[6]*armZZos[104];
   armZZos[104]=armZZos[104] + armZZos[125] + 11./9.*armZZos[121];
   armZZos[104]=mmH*armZZos[104];
   armZZos[105]=armZZos[124] + 5./2.*armZZos[122];
   armZZos[106]=armZZos[124] + armZZos[122];
   armZZos[106]=armZZos[6]*armZZos[106];
   armZZos[104]=1./2.*armZZos[104] + 1./3.*armZZos[105] + 1./2.*
   armZZos[106];
   armZZos[104]=mmH*armZZos[104];
   armZZos[105]=armZZos[14]*armZZos[257];
   armZZos[106]=1 - armZZos[9];
   armZZos[106]=armZZos[16]*armZZos[106];
   armZZos[107]=armZZos[12] - armZZos[78] - 1./4. - armZZos[86];
   armZZos[107]=mmH*armZZos[107];
   armZZos[105]=armZZos[107] + armZZos[106] + armZZos[105] + 
   armZZos[17] - armZZos[19];
   armZZos[105]=mmH*armZZos[105];
   armZZos[106]=armZZos[15]*armZZos[14];
   armZZos[105]=1./3.*armZZos[105] + armZZos[193] + armZZos[106];
   armZZos[105]=armZZos[100]*armZZos[62]*armZZos[105];
   armZZos[98]=1./8.*armZZos[105] + 1./6.*armZZos[104] + armZZos[98];
   armZZos[98]=armZZos[100]*armZZos[98];
   armZZos[97]=armZZos[97] + armZZos[98];
   armZZos[97]=armZZos[100]*armZZos[97];
   armZZos[98]=pow(armZZos[10],2);
   armZZos[98]= - 1./9.*armZZos[98] + armZZos[163] + 1 + armZZos[102];
   armZZos[98]=mmH*armZZos[98];
   armZZos[102]= - 5./2.*armZZos[19] + armZZos[80] + armZZos[138];
   armZZos[102]=1./2.*armZZos[102] - armZZos[18];
   armZZos[104]=1 + 7./9.*armZZos[10];
   armZZos[104]=armZZos[15]*armZZos[104];
   armZZos[98]=1./2.*armZZos[98] + 1./4.*armZZos[104] + 7./36.*
   armZZos[130] + 1./3.*armZZos[102] + armZZos[120];
   armZZos[98]=mmH*armZZos[98];
   armZZos[102]=19*armZZos[16] + 25./4.*armZZos[15];
   armZZos[102]=armZZos[15]*armZZos[102];
   armZZos[102]=armZZos[126] + 1./4.*armZZos[102];
   armZZos[98]=1./9.*armZZos[102] + armZZos[98];
   armZZos[98]=armZZos[62]*armZZos[98]*pow(mmH,2);
   armZZos[102]=armZZos[100]*armZZos[98];
   armZZos[102]=armZZos[98] + 1./2.*armZZos[102];
   armZZos[102]=armZZos[100]*armZZos[102];
   armZZos[104]=armZZos[112]*armZZos[98];
   armZZos[98]=armZZos[98] + 1./2.*armZZos[104];
   armZZos[98]=armZZos[112]*armZZos[98];
   armZZos[98]=armZZos[102] + armZZos[98];
   armZZos[98]=armZZos[2]*armZZos[98];
   armZZos[102]=7./3.*armZZos[76];
   armZZos[104]=3*armZZos[75];
   armZZos[105]=armZZos[104] + 5./2. + armZZos[102];
   armZZos[106]= - 17./6.*armZZos[13];
   armZZos[107]= - 1 + armZZos[128];
   armZZos[107]=armZZos[10]*armZZos[107];
   armZZos[109]=1./18.*armZZos[107];
   armZZos[110]= - armZZos[68] - 1./3.*armZZos[67];
   armZZos[110]=mmH*armZZos[110];
   armZZos[113]=1./2.*armZZos[110];
   armZZos[105]=armZZos[113] + armZZos[109] + armZZos[106] + 1./2.*
   armZZos[78] + armZZos[143] + 1./4.*armZZos[105] + armZZos[73];
   armZZos[105]=mmH*armZZos[105];
   armZZos[114]= - 55./18.*armZZos[16] - 49./6.*armZZos[14] + 7./12.*
   armZZos[18] - 23./6.*armZZos[19] + 37./6.*armZZos[17] + armZZos[115]
    + armZZos[194] - 1./4.*armZZos[80] + 5./4.*armZZos[66] + 
   armZZos[81];
   armZZos[115]=25./4.*armZZos[16];
   armZZos[116]=11*armZZos[14] + armZZos[115];
   armZZos[116]=armZZos[10]*armZZos[116];
   armZZos[117]= - 47./12. - armZZos[10];
   armZZos[117]=armZZos[15]*armZZos[117];
   armZZos[118]=1./12.*armZZos[117];
   armZZos[105]=1./2.*armZZos[105] + armZZos[118] + 1./2.*armZZos[114]
    + 1./9.*armZZos[116];
   armZZos[105]=mmH*armZZos[105];
   armZZos[114]= - 13*armZZos[14];
   armZZos[116]=armZZos[207] + armZZos[114] + armZZos[205];
   armZZos[116]=armZZos[15]*armZZos[116];
   armZZos[120]= - 67./3.*armZZos[14] + armZZos[127];
   armZZos[120]=armZZos[16]*armZZos[120];
   armZZos[120]=43*armZZos[193] + armZZos[120];
   armZZos[116]=1./2.*armZZos[120] + 1./3.*armZZos[116];
   armZZos[105]=1./6.*armZZos[116] + armZZos[105];
   armZZos[105]=1./2.*armZZos[62]*mmH*armZZos[105];
   armZZos[116]= - 19./2. + armZZos[76];
   armZZos[116]=7./3.*armZZos[116] + armZZos[104];
   armZZos[106]=armZZos[113] + armZZos[109] + armZZos[106] - 17./3.*
   armZZos[78] + armZZos[143] + 1./4.*armZZos[116] + armZZos[73];
   armZZos[106]=mmH*armZZos[106];
   armZZos[109]=armZZos[216] + 25./36.*armZZos[16];
   armZZos[109]=armZZos[10]*armZZos[109];
   armZZos[99]=1./2.*armZZos[106] + armZZos[118] + armZZos[109] + 167./
   36.*armZZos[16] + armZZos[99] + 7./24.*armZZos[18] - 23./12.*
   armZZos[19] + 17./3.*armZZos[17] + 1./24.*armZZos[65] - 1./4.*
   armZZos[82] + 37./12.*armZZos[64] - 1./8.*armZZos[80] + 5./8.*
   armZZos[66] - 17./3.*armZZos[81];
   armZZos[99]=mmH*armZZos[99];
   armZZos[106]=armZZos[14] + 5./16.*armZZos[16];
   armZZos[106]=armZZos[16]*armZZos[106];
   armZZos[106]=37./8.*armZZos[193] + armZZos[106];
   armZZos[109]= - 17./9.*armZZos[15] + armZZos[114] + 19./18.*
   armZZos[16];
   armZZos[109]=armZZos[15]*armZZos[109];
   armZZos[99]=1./4.*armZZos[99] + 1./3.*armZZos[106] + 1./8.*
   armZZos[109];
   armZZos[99]=armZZos[112]*armZZos[62]*mmH*armZZos[99];
   armZZos[99]=armZZos[105] + armZZos[99];
   armZZos[99]=armZZos[112]*armZZos[99];
   armZZos[106]=armZZos[207] + armZZos[197] + armZZos[205];
   armZZos[106]=armZZos[15]*armZZos[106];
   armZZos[109]= - 31./3.*armZZos[14] + 5./2.*armZZos[16];
   armZZos[109]=armZZos[16]*armZZos[109];
   armZZos[106]=1./18.*armZZos[106] + 3*armZZos[193] + 1./6.*
   armZZos[109];
   armZZos[102]=armZZos[104] + 33./2. + armZZos[102];
   armZZos[102]=armZZos[143] + 1./4.*armZZos[102] + armZZos[73];
   armZZos[102]=1./8.*armZZos[110] + 1./72.*armZZos[107] - 17./24.*
   armZZos[13] + 1./4.*armZZos[102] + armZZos[78];
   armZZos[102]=mmH*armZZos[102];
   armZZos[104]=armZZos[216] + armZZos[115];
   armZZos[104]=armZZos[10]*armZZos[104];
   armZZos[102]=1./2.*armZZos[102] + 1./48.*armZZos[117] + 1./36.*
   armZZos[104] - 181./144.*armZZos[16] - 21./8.*armZZos[14] + 7./96.*
   armZZos[18] - 23./48.*armZZos[19] + 5./8.*armZZos[17] + 1./96.*
   armZZos[65] - 1./16.*armZZos[82] - 7./16.*armZZos[64] - 1./32.*
   armZZos[80] + 5./32.*armZZos[66] + armZZos[81];
   armZZos[102]=mmH*armZZos[102];
   armZZos[102]=1./4.*armZZos[106] + armZZos[102];
   armZZos[102]=armZZos[62]*mmH*armZZos[102];
   armZZos[104]= - armZZos[16]*armZZos[14];
   armZZos[104]=armZZos[193] + armZZos[104];
   armZZos[106]=1 + armZZos[78];
   armZZos[106]=mmH*armZZos[106];
   armZZos[106]=1./2.*armZZos[106] - armZZos[16] - 5./2.*armZZos[14] + 
   armZZos[157] + armZZos[81] - 1./2.*armZZos[64];
   armZZos[106]=mmH*armZZos[106];
   armZZos[104]=1./2.*armZZos[104] + armZZos[106];
   armZZos[104]=armZZos[100]*armZZos[62]*mmH*armZZos[104];
   armZZos[102]=armZZos[102] + 1./12.*armZZos[104];
   armZZos[102]=armZZos[100]*armZZos[102];
   armZZos[102]=armZZos[105] + armZZos[102];
   armZZos[100]=armZZos[100]*armZZos[102];
   armZZos[98]=1./2.*armZZos[98] + armZZos[100] + armZZos[99];
   armZZos[98]=armZZos[2]*armZZos[98];
   armZZos[99]= - 5*armZZos[193] + 4*armZZos[101];
   armZZos[99]=armZZos[62]*armZZos[99];
   armZZos[96]=armZZos[98] + armZZos[96] + 8*armZZos[99] + armZZos[97];
   armZZos[96]=armZZos[2]*armZZos[96];
   armZZos[97]=4*armZZos[18] + armZZos[7] + 4*armZZos[8];
   armZZos[97]=5*armZZos[97] - armZZos[14];
   armZZos[97]=armZZos[1]*armZZos[97];
   armZZos[98]= - 13./3.*armZZos[14] + 340./27.*armZZos[18] + 11*
   armZZos[7] + 340./27.*armZZos[8];
   armZZos[98]=armZZos[34]*armZZos[98];
   armZZos[99]= - armZZos[1] - 17./27.*armZZos[34];
   armZZos[99]=armZZos[15]*armZZos[99];
   armZZos[97]=52*armZZos[99] + armZZos[97] + armZZos[98];
   armZZos[98]=armZZos[1] + 17./27.*armZZos[34];
   armZZos[98]=armZZos[15]*armZZos[98];
   armZZos[98]=5./3.*armZZos[98] + armZZos[176] + 5./3.*armZZos[184];
   armZZos[98]=armZZos[6]*armZZos[98];
   armZZos[97]=1./3.*armZZos[97] + 4*armZZos[98];
   armZZos[98]= - armZZos[63] + armZZos[108];
   armZZos[98]=armZZos[123]*armZZos[98];
   armZZos[99]= - 23./3.*armZZos[63] + 1./9.*armZZos[83] - armZZos[81];
   armZZos[100]=20./3. + armZZos[123];
   armZZos[100]=armZZos[18]*armZZos[100];
   armZZos[98]=4*armZZos[100] + 4*armZZos[98] + 494./9.*armZZos[17] - 
   34./3.*armZZos[82] + 2*armZZos[99] + armZZos[64];
   armZZos[99]=403./9. + armZZos[119];
   armZZos[100]= - 24*armZZos[123];
   armZZos[101]= - 425./9. + armZZos[100];
   armZZos[101]=armZZos[9]*armZZos[101];
   armZZos[99]=2*armZZos[99] + armZZos[101];
   armZZos[99]=armZZos[14]*armZZos[99];
   armZZos[101]= - 23./3. + armZZos[155];
   armZZos[101]=armZZos[9]*armZZos[101];
   armZZos[101]=armZZos[101] - 21 + armZZos[119];
   armZZos[101]=armZZos[15]*armZZos[101];
   armZZos[98]=4*armZZos[101] + 4*armZZos[16] + 2*armZZos[98] + 
   armZZos[99];
   armZZos[99]=5 + 14*armZZos[79];
   armZZos[99]=1./3.*armZZos[99] - 4*armZZos[78];
   armZZos[99]=mmH*armZZos[99];
   armZZos[98]=2*armZZos[98] + armZZos[99];
   armZZos[98]=armZZos[62]*armZZos[98];
   armZZos[94]=armZZos[96] + armZZos[94] + armZZos[95] + 2*armZZos[97]
    + armZZos[98];
   armZZos[94]=armZZos[2]*armZZos[94];
   armZZos[95]= - armZZos[12] + 2*armZZos[78] - armZZos[79] + 2 + 
   armZZos[86];
   armZZos[95]=armZZos[123]*armZZos[95];
   armZZos[96]= - 82*armZZos[79] + 29 + 82*armZZos[86];
   armZZos[95]=8*armZZos[95] - 82./3.*armZZos[12] + 164./3.*armZZos[78]
    + 1./3.*armZZos[96] + armZZos[140];
   armZZos[95]=armZZos[90]*armZZos[95];
   armZZos[96]= - armZZos[71] + 2*armZZos[84] + armZZos[72];
   armZZos[96]=armZZos[123]*armZZos[96];
   armZZos[97]= - 46*armZZos[71] + 58*armZZos[84] + 23*armZZos[72];
   armZZos[96]=4*armZZos[96] + 1./3.*armZZos[97] + armZZos[70];
   armZZos[96]=armZZos[123]*armZZos[96];
   armZZos[97]= - 41./3. + armZZos[119];
   armZZos[98]=armZZos[9]*armZZos[97];
   armZZos[98]=armZZos[98] - 35./3. + armZZos[119];
   armZZos[98]=armZZos[5]*armZZos[98];
   armZZos[99]=100./3.*armZZos[84] + 7*armZZos[72];
   armZZos[97]=armZZos[9]*armZZos[90]*armZZos[97];
   armZZos[101]=1 + armZZos[79];
   armZZos[101]=armZZos[89]*armZZos[101];
   armZZos[95]=10*armZZos[151] + 8*armZZos[101] + 2*armZZos[98] + 2*
   armZZos[171] + 2*armZZos[97] + 2*armZZos[170] + armZZos[95] + 4*
   armZZos[96] + 41./3.*armZZos[70] + 16*armZZos[69] + 2*armZZos[99] - 
   65*armZZos[71];
   armZZos[95]=mmZ*armZZos[95];
   armZZos[96]= - 767./9. + armZZos[100];
   armZZos[96]=armZZos[123]*armZZos[96];
   armZZos[97]= - 361./9. - 12*armZZos[123];
   armZZos[97]=armZZos[123]*armZZos[97];
   armZZos[97]= - 592./9. + armZZos[97];
   armZZos[97]=armZZos[9]*armZZos[97];
   armZZos[96]=2*armZZos[97] - 2731./9. + 2*armZZos[96];
   armZZos[96]=armZZos[9]*armZZos[96];
   armZZos[97]= - 14*armZZos[123] - 2*armZZos[12] - 2*armZZos[79] - 517.
   /9. + 4*armZZos[77];
   armZZos[97]=armZZos[123]*armZZos[97];
   armZZos[98]= - 5 + armZZos[111];
   armZZos[98]=armZZos[14]*armZZos[98];
   armZZos[99]=armZZos[16]*armZZos[257];
   armZZos[98]=armZZos[99] + armZZos[19] + armZZos[98];
   armZZos[98]=armZZos[5]*armZZos[98];
   armZZos[99]= - 3055./3. + 244*armZZos[77];
   armZZos[100]=armZZos[64] - 4*armZZos[17];
   armZZos[100]=armZZos[90]*armZZos[100];
   armZZos[101]=armZZos[14]*armZZos[172];
   armZZos[102]= - 2*armZZos[90] + armZZos[159];
   armZZos[102]=armZZos[16]*armZZos[102];
   armZZos[104]= - mmH*armZZos[70];
   armZZos[95]=10./3.*armZZos[104] + armZZos[95] + 4*armZZos[98] + 4*
   armZZos[102] + 8*armZZos[101] + armZZos[96] + 4*armZZos[100] + 4*
   armZZos[97] - 58*armZZos[12] + 8*armZZos[78] - 190./9.*armZZos[79]
    + 1./3.*armZZos[99] - 2*armZZos[86];
   armZZos[95]=armZZos[62]*armZZos[95];
   armZZos[96]=11*armZZos[29] + 53./3. + 11*armZZos[20];
   armZZos[96]=armZZos[123]*armZZos[96];
   armZZos[97]= - 5*armZZos[26];
   armZZos[98]=5*armZZos[12];
   armZZos[99]=217./9. + armZZos[103];
   armZZos[99]=armZZos[9]*armZZos[99];
   armZZos[96]=200./243.*armZZos[34] + 80./81.*armZZos[1] + 5./9.*
   armZZos[99] + 1./3.*armZZos[96] - 134./81.*armZZos[137] + 
   armZZos[98] + 680./81.*armZZos[28] + 101./9.*armZZos[29] + 55./9.*
   armZZos[20] + 808./81. + armZZos[97];
   armZZos[96]=armZZos[34]*armZZos[96];
   armZZos[99]=2*armZZos[21];
   armZZos[100]=armZZos[99] + armZZos[24];
   armZZos[101]=armZZos[123]*armZZos[21];
   armZZos[100]=2*armZZos[101] - 16*armZZos[22] + 4*armZZos[100] + 
   armZZos[23];
   armZZos[100]=armZZos[1]*armZZos[100];
   armZZos[101]=armZZos[21] + 2./9.*armZZos[23];
   armZZos[101]=armZZos[123]*armZZos[101];
   armZZos[99]=armZZos[99] + 17./81.*armZZos[24];
   armZZos[99]=2*armZZos[101] - 272./81.*armZZos[22] + 4*armZZos[99] + 
   17./9.*armZZos[23];
   armZZos[99]=armZZos[34]*armZZos[99];
   armZZos[101]=4./3.*armZZos[132] + armZZos[158];
   armZZos[101]=armZZos[5]*armZZos[101];
   armZZos[99]=4./3.*armZZos[101] + 1./3.*armZZos[100] + armZZos[99];
   armZZos[99]=mmZ*armZZos[99];
   armZZos[97]= - 10*armZZos[137] + armZZos[98] + 40*armZZos[28] + 13*
   armZZos[29] + 17./3.*armZZos[20] + 136./9. + armZZos[97];
   armZZos[98]=5./3.*armZZos[29] + 3 + 5./3.*armZZos[20];
   armZZos[98]=armZZos[123]*armZZos[98];
   armZZos[100]=157./9. + armZZos[103];
   armZZos[100]=armZZos[9]*armZZos[100];
   armZZos[97]=8./27.*armZZos[1] + 1./3.*armZZos[100] + 1./3.*
   armZZos[97] + armZZos[98];
   armZZos[97]=armZZos[1]*armZZos[97];
   armZZos[98]=8./9. + armZZos[123];
   armZZos[98]=armZZos[9]*armZZos[98];
   armZZos[98]= - 2./9.*armZZos[1] + armZZos[98] + 11./9. + 
   armZZos[123];
   armZZos[98]=armZZos[1]*armZZos[98];
   armZZos[100]=5*armZZos[123];
   armZZos[101]=52./9. + armZZos[100];
   armZZos[101]=armZZos[9]*armZZos[101];
   armZZos[100]= - 50./27.*armZZos[34] - 20./9.*armZZos[1] + 
   armZZos[101] + 127./27. + armZZos[100];
   armZZos[100]=armZZos[34]*armZZos[100];
   armZZos[98]=armZZos[98] + 1./3.*armZZos[100];
   armZZos[100]=5*armZZos[1] + 67./27.*armZZos[34];
   armZZos[100]=armZZos[6]*armZZos[100];
   armZZos[98]=2*armZZos[98] + 1./3.*armZZos[100];
   armZZos[98]=armZZos[6]*armZZos[98];
   armZZos[96]=armZZos[99] + 2*armZZos[98] + armZZos[97] + armZZos[96];

      mZZosret = armZZos[91] + armZZos[92] + armZZos[93] + armZZos[94]
       + armZZos[95] + 2*armZZos[96];
      return mZZosret;
}
