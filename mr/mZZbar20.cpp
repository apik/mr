#include <ZZ.hpp>
std::complex<long double>
ZZ<OS>::m20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armZZbar[391], mZZbarret;

    armZZbar[1]=double(nL + nH);
    armZZbar[2]=pow(CW,-1);
    armZZbar[3]=pow(MMH,-1);
    armZZbar[4]=pow(MMZ,-1);
    armZZbar[5]=pow(SW,-1);
    armZZbar[6]=std::real(Tsil::B(0,0,MMZ,mu2));
    armZZbar[7]=std::real(Tsil::B(0,0,MMW,mu2));
    armZZbar[8]=Tsil::I2(0,0,MMZ,mu2);
    armZZbar[9]=Tsil::I2(0,0,MMW,mu2);
    armZZbar[10]=Tsil::B(MMZ,MMH,MMZ,mu2);
    armZZbar[11]=Tsil::B(MMW,MMH,MMW,mu2);
    armZZbar[12]=Tsil::B(MMW,MMZ,MMW,mu2);
    armZZbar[13]=Tsil::B(MMW,MMW,MMZ,mu2);
    armZZbar[14]=Tsil::B(0,MMH,MMZ,mu2);
    armZZbar[15]=Tsil::Beps(MMZ,MMH,MMZ,mu2);
    armZZbar[16]=Tsil::Beps(MMW,MMW,MMZ,mu2);
    armZZbar[17]=Tsil::A(MMH,mu2);
    armZZbar[18]=Tsil::A(MMZ,mu2);
    armZZbar[19]=Tsil::A(MMW,mu2);
    armZZbar[20]=Tsil::Aeps(MMH,mu2);
    armZZbar[21]=Tsil::Aeps(MMZ,mu2);
    armZZbar[22]=Tsil::Aeps(MMW,mu2);
    armZZbar[23]=std::real(Tsil::B(0,MMW,MMZ,mu2));
    armZZbar[24]=protW0W00->M(0);
    armZZbar[25]=prot0000Z->M(0);
    armZZbar[26]=prot0000W->M(0);
    armZZbar[27]=prot00000->M(0);
    armZZbar[28]=protHZ00->Uxzuv(0);
    armZZbar[29]=protW0W00->Uzxyv(0);
    armZZbar[30]=protHZ00->Txuv(0);
    armZZbar[31]=prot0000Z->Tvxu(0);
    armZZbar[32]=protW0W00->Txuv(0);
    armZZbar[33]=double(nH);
    armZZbar[34]=pow(MMt,-1);
    armZZbar[35]=Tsil::B(MMt,MMt,MMZ,mu2);
    armZZbar[36]=Tsil::B(0,MMt,MMW,mu2);
    armZZbar[37]=Tsil::A(MMt,mu2);
    armZZbar[38]=Tsil::B(MMt,MMt,MMH,mu2);
    armZZbar[39]=double(nL);
    armZZbar[40]=Tsil::I2(MMH,MMt,MMt,mu2);
    armZZbar[41]=Tsil::I2(MMZ,MMt,MMt,mu2);
    armZZbar[42]=Tsil::I2(0,MMW,MMt,mu2);
    armZZbar[43]=Tsil::B(MMH,MMH,MMH,mu2);
    armZZbar[44]=Tsil::B(MMH,MMt,MMt,mu2);
    armZZbar[45]=Tsil::B(MMZ,MMZ,MMH,mu2);
    armZZbar[46]=Tsil::B(MMZ,MMt,MMt,mu2);
    armZZbar[47]=Tsil::B(MMW,MMW,MMH,mu2);
    armZZbar[48]=std::real(Tsil::B(0,MMW,MMt,mu2));
    armZZbar[49]=Tsil::Beps(MMt,MMt,MMZ,mu2);
    armZZbar[50]=Tsil::Aeps(MMt,mu2);
    armZZbar[51]=protWtWt0->M(0);
    armZZbar[52]=protW0W0t->M(0);
    armZZbar[53]=prottZtHt->M(0);
    armZZbar[54]=protttttH->M(0);
    armZZbar[55]=protttttZ->M(0);
    armZZbar[56]=prottttt0->M(0);
    armZZbar[57]=prott0t0W->M(0);
    armZZbar[58]=prottttt0->Vzxyv(0);
    armZZbar[59]=prottZtHt->Uuyxv(0);
    armZZbar[60]=protWtWt0->Uzxyv(0);
    armZZbar[61]=protttttH->Uzxyv(0);
    armZZbar[62]=protttttZ->Uzxyv(0);
    armZZbar[63]=protWtWt0->Uuyxv(0);
    armZZbar[64]=prottZtHt->Tuxv(0);
    armZZbar[65]=protWtWt0->Tzyv(0);
    armZZbar[66]=protWtWt0->Tyzv(0);
    armZZbar[67]=prottZtHt->Txuv(0);
    armZZbar[68]=prottZtHt->Suxv(0);
    armZZbar[69]=prottZtHt->Svyz(0);
    armZZbar[70]=protWtWt0->Svyz(0);
    armZZbar[71]=prottttt0->Suxv(0);
    armZZbar[72]=prottZtHt->Uyuzv(0);
    armZZbar[73]=double(boson);
    armZZbar[74]=Tsil::I2(MMH,MMH,MMH,mu2);
    armZZbar[75]=Tsil::I2(MMH,MMZ,MMZ,mu2);
    armZZbar[76]=Tsil::I2(MMH,MMW,MMW,mu2);
    armZZbar[77]=Tsil::I2(MMW,MMW,MMZ,mu2);
    armZZbar[78]=protZHHZZ->M(0);
    armZZbar[79]=protZZHHH->M(0);
    armZZbar[80]=protZWHWW->M(0);
    armZZbar[81]=protWWWWH->M(0);
    armZZbar[82]=protWWWWZ->M(0);
    armZZbar[83]=protWWWW0->M(0);
    armZZbar[84]=protZHHZZ->Uzxyv(0);
    armZZbar[85]=protZWHWW->Uzxyv(0);
    armZZbar[86]=protZZHHH->Uyuzv(0);
    armZZbar[87]=protZHHZZ->Uuyxv(0);
    armZZbar[88]=protZWHWW->Uuyxv(0);
    armZZbar[89]=protZWHWW->Tzyv(0);
    armZZbar[90]=protZWHWW->Tyzv(0);
    armZZbar[91]=protZZHHH->Suxv(0);
    armZZbar[92]=protZWHWW->Svyz(0);
    armZZbar[93]=protWZWHW->Svyz(0);
    armZZbar[94]=protWWWW0->Suxv(0);
    armZZbar[95]=protWWWW0->Vzxyv(0);
    armZZbar[96]=protZWHWW->Uxzuv(0);
    armZZbar[97]=protZWHWW->Uyuzv(0);
    armZZbar[98]=1/(MMt - MMW);
    armZZbar[99]=1/(4*MMt - MMZ);
    armZZbar[100]=1/( - 4*MMW + MMH);
    armZZbar[101]=1/( - 4*MMZ + MMH);
   armZZbar[102]= - 161./27. + 5*armZZbar[12];
   armZZbar[103]= - 2*armZZbar[13];
   armZZbar[102]=1./8.*armZZbar[102] + armZZbar[103];
   armZZbar[104]=1./9.*armZZbar[1];
   armZZbar[105]=1./3.*armZZbar[39];
   armZZbar[106]=armZZbar[105] - 1./8. + armZZbar[104];
   armZZbar[106]=armZZbar[6]*armZZbar[106];
   armZZbar[107]=1./4.*armZZbar[19]*armZZbar[98];
   armZZbar[108]=pow(armZZbar[98],2);
   armZZbar[109]= - armZZbar[19]*armZZbar[108];
   armZZbar[110]= - armZZbar[98] + armZZbar[109];
   armZZbar[111]=3./2.*armZZbar[110] + armZZbar[101];
   armZZbar[111]=MMZ*armZZbar[111];
   armZZbar[112]= - 25./9. - armZZbar[7];
   armZZbar[112]=7./3.*armZZbar[35] + 1./2.*armZZbar[112] - 
   armZZbar[36];
   armZZbar[113]=armZZbar[1]*armZZbar[112];
   armZZbar[112]=armZZbar[39]*armZZbar[112];
   armZZbar[114]= - MMZ*armZZbar[101];
   armZZbar[115]=7./3. + armZZbar[114];
   armZZbar[115]=armZZbar[10]*armZZbar[115];
   armZZbar[116]= - armZZbar[101]*armZZbar[17];
   armZZbar[117]=armZZbar[18]*armZZbar[101];
   armZZbar[118]= - 1./2.*armZZbar[11];
   armZZbar[102]=5*armZZbar[106] + 1./2.*armZZbar[117] + 1./2.*
   armZZbar[115] + armZZbar[112] + 1./3.*armZZbar[113] + 5./2.*
   armZZbar[35] + 1./2.*armZZbar[111] + 1./4.*armZZbar[116] + 
   armZZbar[107] + 11*armZZbar[102] + armZZbar[118];
   armZZbar[102]=armZZbar[6]*armZZbar[102];
   armZZbar[106]=2*armZZbar[22];
   armZZbar[111]= - 1./2.*armZZbar[8];
   armZZbar[112]= - 1./2.*armZZbar[41];
   armZZbar[113]=armZZbar[106] + 3*armZZbar[50] + armZZbar[112] - 2*
   armZZbar[42] + armZZbar[111];
   armZZbar[115]=11*armZZbar[19];
   armZZbar[119]=1./2.*armZZbar[17];
   armZZbar[120]=31./2.*MMZ + armZZbar[115] + armZZbar[119];
   armZZbar[120]=armZZbar[35]*armZZbar[120];
   armZZbar[121]= - 3./2.*armZZbar[36];
   armZZbar[122]=7./2.*armZZbar[35];
   armZZbar[123]=armZZbar[122] + 1 + armZZbar[121];
   armZZbar[123]=armZZbar[18]*armZZbar[123];
   armZZbar[124]=7*armZZbar[19];
   armZZbar[125]=5*armZZbar[18] + 23./2.*MMZ + armZZbar[124] + 
   armZZbar[119];
   armZZbar[125]=armZZbar[6]*armZZbar[125];
   armZZbar[126]=1./2.*armZZbar[20];
   armZZbar[127]=3*armZZbar[21];
   armZZbar[128]= - 12*MMZ;
   armZZbar[129]=3*armZZbar[46];
   armZZbar[130]=6*armZZbar[48] + 37./2. + armZZbar[129];
   armZZbar[130]=armZZbar[37]*armZZbar[130];
   armZZbar[131]= - 1./2.*armZZbar[17];
   armZZbar[132]=2*armZZbar[19];
   armZZbar[113]=armZZbar[130] + 1./2.*armZZbar[125] + armZZbar[123] + 
   1./2.*armZZbar[120] + armZZbar[128] + armZZbar[131] + armZZbar[132]
    + armZZbar[127] + 3*armZZbar[113] + armZZbar[126];
   armZZbar[113]=armZZbar[3]*armZZbar[113];
   armZZbar[120]= - armZZbar[19] - MMZ;
   armZZbar[123]= - 1./2.*armZZbar[18];
   armZZbar[125]=armZZbar[120] + armZZbar[123];
   armZZbar[125]=armZZbar[6]*armZZbar[125];
   armZZbar[130]=MMZ*armZZbar[36];
   armZZbar[133]=armZZbar[35]*armZZbar[120];
   armZZbar[134]= - 1./2.*armZZbar[35];
   armZZbar[135]=armZZbar[36] + armZZbar[134];
   armZZbar[136]=armZZbar[18]*armZZbar[135];
   armZZbar[137]=armZZbar[19]*armZZbar[36];
   armZZbar[125]=1./2.*armZZbar[125] + 1./2.*armZZbar[136] + 1./2.*
   armZZbar[133] + armZZbar[137] + armZZbar[130];
   armZZbar[125]=armZZbar[3]*armZZbar[125];
   armZZbar[130]=33*armZZbar[13] + 59./9. - 33./2.*armZZbar[12];
   armZZbar[130]=1./2.*armZZbar[130] + armZZbar[11];
   armZZbar[133]=1 + armZZbar[7];
   armZZbar[133]= - armZZbar[35] + 1./2.*armZZbar[133] + armZZbar[36];
   armZZbar[136]=armZZbar[1]*armZZbar[133];
   armZZbar[133]=armZZbar[39]*armZZbar[133];
   armZZbar[138]= - 1./3.*armZZbar[1];
   armZZbar[139]=armZZbar[138] - armZZbar[39];
   armZZbar[140]=armZZbar[6]*armZZbar[139];
   armZZbar[133]=armZZbar[140] - armZZbar[10] + armZZbar[133] + 1./2.*
   armZZbar[130] + 1./3.*armZZbar[136];
   armZZbar[133]=armZZbar[6]*armZZbar[133];
   armZZbar[136]= - armZZbar[7] - armZZbar[36];
   armZZbar[141]=1./3. + armZZbar[7];
   armZZbar[142]=armZZbar[35]*armZZbar[141];
   armZZbar[136]=1./3.*armZZbar[136] + 1./2.*armZZbar[142];
   armZZbar[142]=armZZbar[1]*armZZbar[136];
   armZZbar[136]=armZZbar[39]*armZZbar[136];
   armZZbar[143]= - 59./9. - 33./2.*armZZbar[13];
   armZZbar[143]=armZZbar[36]*armZZbar[143];
   armZZbar[144]=armZZbar[35]*armZZbar[130];
   armZZbar[145]= - armZZbar[36] + armZZbar[35];
   armZZbar[146]=MMt*armZZbar[3]*armZZbar[145];
   armZZbar[147]=armZZbar[12] - armZZbar[13];
   armZZbar[148]=1./3. + armZZbar[36];
   armZZbar[149]=armZZbar[148] - armZZbar[35];
   armZZbar[150]=armZZbar[10]*armZZbar[149];
   armZZbar[151]= - 1./3.*armZZbar[11];
   armZZbar[125]=3./2.*armZZbar[146] + 3*armZZbar[125] + armZZbar[133]
    + armZZbar[150] + armZZbar[136] + 1./3.*armZZbar[142] + 1./2.*
   armZZbar[144] + 1./2.*armZZbar[143] + 11./4.*armZZbar[147] + 
   armZZbar[151];
   armZZbar[133]=pow(armZZbar[5],2);
   armZZbar[125]=armZZbar[133]*armZZbar[125];
   armZZbar[136]= - 7 + 27./4.*armZZbar[63];
   armZZbar[142]= - 1./2.*armZZbar[61];
   armZZbar[136]= - 47./8.*armZZbar[49] - 81./16.*armZZbar[66] + 
   armZZbar[142] - 27./4.*armZZbar[65] + 1./2.*armZZbar[136] + 3*
   armZZbar[62];
   armZZbar[143]= - 1./4.*armZZbar[44];
   armZZbar[136]= - armZZbar[15] + armZZbar[143] + 27./16.*armZZbar[48]
    + 3./4.*armZZbar[46] + 3./2.*armZZbar[67] + 1./2.*armZZbar[136] + 
   armZZbar[59];
   armZZbar[136]=armZZbar[99]*armZZbar[136];
   armZZbar[144]= - armZZbar[65] + armZZbar[66];
   armZZbar[144]=armZZbar[108]*armZZbar[144];
   armZZbar[146]=MMZ*armZZbar[144];
   armZZbar[152]=1./2.*armZZbar[66] + 33./16. + armZZbar[65];
   armZZbar[152]=armZZbar[98]*armZZbar[152];
   armZZbar[153]=armZZbar[50]*armZZbar[108];
   armZZbar[154]= - armZZbar[22]*armZZbar[108];
   armZZbar[155]=3./2.*armZZbar[49] - 1./2.*armZZbar[62] + 9./8. - 
   armZZbar[63];
   armZZbar[155]=armZZbar[34]*armZZbar[155];
   armZZbar[156]=armZZbar[19]*armZZbar[108];
   armZZbar[157]=1./2.*armZZbar[59];
   armZZbar[158]= - armZZbar[15] + 3*armZZbar[64] + armZZbar[157] - 3./
   2.*armZZbar[14] + 1./2.*armZZbar[28] + 3./8. - armZZbar[30];
   armZZbar[158]=armZZbar[101]*armZZbar[158];
   armZZbar[136]=3./4.*armZZbar[146] + armZZbar[158] + 3./4.*
   armZZbar[156] + 3./16.*armZZbar[136] + 1./2.*armZZbar[155] + 3./4.*
   armZZbar[154] + 3./4.*armZZbar[153] + armZZbar[152] - 4*armZZbar[57]
    + armZZbar[25] - 6*armZZbar[51] + armZZbar[55] - 6*armZZbar[52];
   armZZbar[136]=MMZ*armZZbar[136];
   armZZbar[146]= - 3*armZZbar[72];
   armZZbar[152]= - 1./2.*armZZbar[59];
   armZZbar[155]= - 3*armZZbar[64];
   armZZbar[158]= - 3*armZZbar[38];
   armZZbar[159]=armZZbar[158] + 7./2.*armZZbar[15] + armZZbar[67] + 
   armZZbar[155] + armZZbar[152] + 1./4. + armZZbar[146];
   armZZbar[159]=armZZbar[101]*armZZbar[159];
   armZZbar[160]=47./2. + armZZbar[129];
   armZZbar[161]= - 1 + armZZbar[47];
   armZZbar[161]=2*armZZbar[161] + armZZbar[45];
   armZZbar[161]=armZZbar[3]*armZZbar[37]*armZZbar[161];
   armZZbar[162]=17./4.*armZZbar[35];
   armZZbar[160]=9*armZZbar[161] + armZZbar[162] + armZZbar[121] - 9./2.
   *armZZbar[38] + 1./2.*armZZbar[160] + 3*armZZbar[48];
   armZZbar[160]=armZZbar[3]*armZZbar[160];
   armZZbar[161]=armZZbar[101]*armZZbar[38];
   armZZbar[163]=3*armZZbar[161];
   armZZbar[164]=armZZbar[35]*armZZbar[101];
   armZZbar[165]=armZZbar[163] + 1./2.*armZZbar[164];
   armZZbar[165]=armZZbar[10]*armZZbar[165];
   armZZbar[166]=pow(armZZbar[3],2);
   armZZbar[167]=MMt*armZZbar[166]*armZZbar[38];
   armZZbar[168]= - armZZbar[35]*armZZbar[101];
   armZZbar[169]=1./2.*armZZbar[168];
   armZZbar[170]=2*armZZbar[53];
   armZZbar[159]=18*armZZbar[167] + armZZbar[160] + armZZbar[165] + 
   armZZbar[169] + armZZbar[159] + armZZbar[170] + 4*armZZbar[57] + 7*
   armZZbar[51] - 2*armZZbar[55] + 7./2.*armZZbar[52];
   armZZbar[159]=MMt*armZZbar[159];
   armZZbar[160]=3*armZZbar[44];
   armZZbar[165]=armZZbar[160] - 81./4.*armZZbar[48] + 151./2. - 9*
   armZZbar[46];
   armZZbar[165]=armZZbar[99]*armZZbar[165];
   armZZbar[171]= - MMZ*armZZbar[108];
   armZZbar[172]=armZZbar[34] - 9./4.*armZZbar[99];
   armZZbar[172]=armZZbar[35]*armZZbar[172];
   armZZbar[173]=MMZ*armZZbar[108];
   armZZbar[174]= - armZZbar[98] + 3./4.*armZZbar[173];
   armZZbar[174]=armZZbar[6]*armZZbar[174];
   armZZbar[175]=armZZbar[10]*armZZbar[99];
   armZZbar[165]=armZZbar[174] + 3./4.*armZZbar[175] + 1./4.*
   armZZbar[172] + 3./4.*armZZbar[171] - 2*armZZbar[101] + 1./16.*
   armZZbar[165] - 5./4.*armZZbar[98] - 3*armZZbar[34];
   armZZbar[165]=armZZbar[37]*armZZbar[165];
   armZZbar[172]= - 5*armZZbar[31];
   armZZbar[174]= - 23./24.*armZZbar[14];
   armZZbar[176]= - 5./8.*armZZbar[59];
   armZZbar[177]=pow(Pi,2);
   armZZbar[178]=armZZbar[176] + 5./4.*armZZbar[177] - 33*armZZbar[16]
    - 1037./128.*armZZbar[49] + 141./256.*armZZbar[66] + armZZbar[174]
    - 3./32.*armZZbar[61] + 207./64.*armZZbar[65] + 49./16.*
   armZZbar[62] - armZZbar[28] + 657./128.*armZZbar[63] + armZZbar[172]
    + 1295./216. + 33*armZZbar[60];
   armZZbar[179]=armZZbar[160] + 111./4.*armZZbar[48] - 16243./27. + 15
   *armZZbar[46];
   armZZbar[180]=55*armZZbar[12];
   armZZbar[179]=1./8.*armZZbar[179] + armZZbar[180];
   armZZbar[181]= - 1 - armZZbar[46];
   armZZbar[182]= - 9./4.*armZZbar[48];
   armZZbar[183]=armZZbar[181] + armZZbar[182];
   armZZbar[183]=3*armZZbar[183] + armZZbar[44];
   armZZbar[183]=armZZbar[99]*armZZbar[183];
   armZZbar[183]=3./32.*armZZbar[183] + armZZbar[101];
   armZZbar[183]=MMZ*armZZbar[183];
   armZZbar[184]= - armZZbar[19]*armZZbar[34];
   armZZbar[185]=1./2.*armZZbar[116];
   armZZbar[179]= - 5./4.*armZZbar[35] + armZZbar[183] + armZZbar[185]
    + armZZbar[184] - armZZbar[11] + 1./4.*armZZbar[179] - 55*
   armZZbar[13];
   armZZbar[179]=armZZbar[35]*armZZbar[179];
   armZZbar[183]= - 11./9. - armZZbar[7];
   armZZbar[186]=armZZbar[35]*armZZbar[183];
   armZZbar[187]=1./3.*armZZbar[7];
   armZZbar[188]=1./3.*armZZbar[36];
   armZZbar[186]=1./2.*armZZbar[186] + armZZbar[188] + armZZbar[187] + 
   2./9. + armZZbar[177];
   armZZbar[189]=armZZbar[1]*armZZbar[186];
   armZZbar[186]=armZZbar[39]*armZZbar[186];
   armZZbar[190]= - 1./4.*armZZbar[41];
   armZZbar[191]= - 1./2.*armZZbar[22];
   armZZbar[192]=armZZbar[70] - 1./2.*armZZbar[42];
   armZZbar[193]=armZZbar[191] + armZZbar[190] + armZZbar[192] + 1./2.*
   armZZbar[69];
   armZZbar[193]=armZZbar[34]*armZZbar[193];
   armZZbar[194]= - armZZbar[70] + 1./2.*armZZbar[42];
   armZZbar[194]=9./8.*armZZbar[194] - armZZbar[69];
   armZZbar[195]=1./2.*armZZbar[68];
   armZZbar[196]= - 1./4.*armZZbar[40];
   armZZbar[194]=armZZbar[196] + 3*armZZbar[194] + armZZbar[195];
   armZZbar[197]=5./4.*armZZbar[21];
   armZZbar[194]=armZZbar[197] + 59./32.*armZZbar[50] + 3./4.*
   armZZbar[194] + armZZbar[41];
   armZZbar[194]=armZZbar[99]*armZZbar[194];
   armZZbar[198]= - 1./2.*armZZbar[20];
   armZZbar[199]= - armZZbar[17] + armZZbar[21] + armZZbar[198] + 
   armZZbar[112] + armZZbar[111] + armZZbar[68];
   armZZbar[199]=armZZbar[101]*armZZbar[199];
   armZZbar[200]= - MMZ*armZZbar[99];
   armZZbar[201]=11./3. + armZZbar[114];
   armZZbar[201]=armZZbar[35]*armZZbar[201];
   armZZbar[202]=23./48. - armZZbar[36];
   armZZbar[201]=1./2.*armZZbar[201] + armZZbar[202] + 3./16.*
   armZZbar[200];
   armZZbar[201]=armZZbar[10]*armZZbar[201];
   armZZbar[203]= - 3*armZZbar[49] + 7 + 3*armZZbar[61];
   armZZbar[203]=armZZbar[15] + 3./4.*armZZbar[44] - 5./2.*armZZbar[67]
    + 3./2.*armZZbar[64] + 1./4.*armZZbar[203] - armZZbar[59];
   armZZbar[203]=armZZbar[99]*armZZbar[203];
   armZZbar[204]= - armZZbar[35]*armZZbar[99]*armZZbar[44];
   armZZbar[203]=armZZbar[175] + armZZbar[203] + 3./4.*armZZbar[204];
   armZZbar[203]=MMH*armZZbar[203];
   armZZbar[205]= - armZZbar[34] - 7./8.*armZZbar[99];
   armZZbar[206]= - armZZbar[34] + 9./8.*armZZbar[99];
   armZZbar[206]=1./2.*armZZbar[206] + armZZbar[101];
   armZZbar[206]=armZZbar[35]*armZZbar[206];
   armZZbar[205]=1./2.*armZZbar[205] + armZZbar[206];
   armZZbar[205]=armZZbar[18]*armZZbar[205];
   armZZbar[206]=armZZbar[98]*armZZbar[70];
   armZZbar[207]= - armZZbar[50]*armZZbar[98];
   armZZbar[208]= - armZZbar[22]*armZZbar[98];
   armZZbar[209]= - armZZbar[21]*armZZbar[34];
   armZZbar[210]=11./2.*armZZbar[13];
   armZZbar[211]=113./9. + armZZbar[210];
   armZZbar[211]=armZZbar[36]*armZZbar[211];
   armZZbar[212]=81./32.*armZZbar[99] - 5*armZZbar[98] - armZZbar[34];
   armZZbar[212]=armZZbar[19]*armZZbar[212];
   armZZbar[213]= - armZZbar[17]*armZZbar[99];
   armZZbar[214]= - 3./64.*armZZbar[44];
   armZZbar[215]= - 55./12.*armZZbar[12];
   armZZbar[216]=1./3.*armZZbar[11];
   armZZbar[102]=armZZbar[125] + armZZbar[159] + 1./16.*armZZbar[203]
    + armZZbar[113] + armZZbar[165] + armZZbar[102] + 1./2.*
   armZZbar[205] + armZZbar[201] + armZZbar[186] + 1./3.*armZZbar[189]
    + 1./2.*armZZbar[179] + armZZbar[136] + armZZbar[199] + 3./16.*
   armZZbar[213] + 1./2.*armZZbar[212] + 1./2.*armZZbar[211] + 1./2.*
   armZZbar[194] + armZZbar[216] - 77./12.*armZZbar[13] + armZZbar[215]
    + 1./4.*armZZbar[209] + 13./16.*armZZbar[15] + armZZbar[193] + 
   armZZbar[214] + 1./2.*armZZbar[208] + 7./4.*armZZbar[207] + 3./2.*
   armZZbar[206] - 303./256.*armZZbar[48] - 39./64.*armZZbar[46] + 9./
   32.*armZZbar[67] + 1./2.*armZZbar[178] + armZZbar[64];
   armZZbar[102]=armZZbar[133]*armZZbar[102];
   armZZbar[113]= - 2./3.*armZZbar[30];
   armZZbar[125]=1./3.*armZZbar[28];
   armZZbar[136]= - 1./3.*armZZbar[59];
   armZZbar[159]= - 2*armZZbar[64] + armZZbar[136] - armZZbar[14] + 
   armZZbar[125] - 9./4. + armZZbar[113];
   armZZbar[159]=armZZbar[101]*armZZbar[159];
   armZZbar[165]=armZZbar[65] - armZZbar[66];
   armZZbar[165]=MMZ*armZZbar[108]*armZZbar[165];
   armZZbar[178]=4*armZZbar[56] + armZZbar[27];
   armZZbar[178]=1./3.*armZZbar[178] - 10*armZZbar[55];
   armZZbar[179]= - 17./3.*armZZbar[66] - 23 - 25./3.*armZZbar[65];
   armZZbar[179]=armZZbar[98]*armZZbar[179];
   armZZbar[186]= - armZZbar[50]*armZZbar[108];
   armZZbar[108]=armZZbar[22]*armZZbar[108];
   armZZbar[189]= - 13./6.*armZZbar[49] - 1./2.*armZZbar[65] + 5./6.*
   armZZbar[62] - 17./8. + 4./3.*armZZbar[63];
   armZZbar[189]=armZZbar[34]*armZZbar[189];
   armZZbar[193]=209 - 225*armZZbar[63];
   armZZbar[193]=1423./24.*armZZbar[49] + 675./16.*armZZbar[66] + 5./2.
   *armZZbar[61] + 531./8.*armZZbar[65] + 1./8.*armZZbar[193] - 39*
   armZZbar[62];
   armZZbar[193]=5*armZZbar[15] + 5./4.*armZZbar[44] - 225./16.*
   armZZbar[48] - 39./4.*armZZbar[46] - 15./2.*armZZbar[67] + 1./2.*
   armZZbar[193] - 5*armZZbar[59];
   armZZbar[193]=armZZbar[99]*armZZbar[193];
   armZZbar[108]=7./4.*armZZbar[165] + armZZbar[159] + armZZbar[109] + 
   1./8.*armZZbar[193] + armZZbar[189] + armZZbar[108] + armZZbar[186]
    + 1./4.*armZZbar[179] + 8*armZZbar[57] - 2./3.*armZZbar[25] + 16*
   armZZbar[51] + 1./3.*armZZbar[178] + 12*armZZbar[52];
   armZZbar[108]=MMZ*armZZbar[108];
   armZZbar[109]= - armZZbar[19]*armZZbar[98];
   armZZbar[109]=armZZbar[185] + 1./4.*armZZbar[109] - 5./6.*
   armZZbar[11] + 289./4.*armZZbar[13] + 1967./54. + armZZbar[12];
   armZZbar[165]=31./3. - armZZbar[7];
   armZZbar[178]= - 8./9.*armZZbar[35];
   armZZbar[165]=armZZbar[178] + 5./54.*armZZbar[165] - armZZbar[36];
   armZZbar[165]=armZZbar[1]*armZZbar[165];
   armZZbar[179]= - 5*armZZbar[7];
   armZZbar[186]=41 + armZZbar[179];
   armZZbar[189]= - 11*armZZbar[36];
   armZZbar[186]=1./2.*armZZbar[186] + armZZbar[189];
   armZZbar[186]=1./3.*armZZbar[186] - 8*armZZbar[35];
   armZZbar[186]=armZZbar[39]*armZZbar[186];
   armZZbar[193]=1./3.*armZZbar[101] + armZZbar[98] + armZZbar[156];
   armZZbar[193]=MMZ*armZZbar[193];
   armZZbar[194]=2 + armZZbar[114];
   armZZbar[194]=1./3.*armZZbar[10]*armZZbar[194];
   armZZbar[199]=1./3.*armZZbar[117];
   armZZbar[201]= - 4./3.*armZZbar[1];
   armZZbar[203]= - 4*armZZbar[39] + 5./4. + armZZbar[201];
   armZZbar[203]=armZZbar[6]*armZZbar[203];
   armZZbar[205]= - 7./2.*armZZbar[35];
   armZZbar[109]=1./3.*armZZbar[203] + armZZbar[199] + armZZbar[194] + 
   1./3.*armZZbar[186] + armZZbar[165] + armZZbar[205] + 1./3.*
   armZZbar[109] + armZZbar[193];
   armZZbar[109]=armZZbar[6]*armZZbar[109];
   armZZbar[165]= - 11*armZZbar[19] + armZZbar[17];
   armZZbar[186]=13./4.*armZZbar[18];
   armZZbar[165]=armZZbar[186] + 1./2.*armZZbar[165] - 7*MMZ;
   armZZbar[165]=armZZbar[6]*armZZbar[165];
   armZZbar[193]= - 2*armZZbar[50];
   armZZbar[203]=armZZbar[193] - armZZbar[8] + armZZbar[41];
   armZZbar[211]=29./9.*armZZbar[19];
   armZZbar[212]= - 2*armZZbar[36];
   armZZbar[217]=49./9. + armZZbar[212];
   armZZbar[217]=MMZ*armZZbar[217];
   armZZbar[218]= - 47*armZZbar[19] - armZZbar[17];
   armZZbar[218]=1./2.*armZZbar[218] - 34*MMZ;
   armZZbar[218]=armZZbar[35]*armZZbar[218];
   armZZbar[219]= - 35./12.*armZZbar[35] + 10./9. + armZZbar[121];
   armZZbar[219]=armZZbar[18]*armZZbar[219];
   armZZbar[220]= - 2*armZZbar[46];
   armZZbar[221]= - 71./3. + armZZbar[220];
   armZZbar[221]=armZZbar[37]*armZZbar[221];
   armZZbar[165]=armZZbar[221] + 1./3.*armZZbar[165] + armZZbar[219] + 
   1./3.*armZZbar[218] + 2*armZZbar[217] + armZZbar[203] + 
   armZZbar[211];
   armZZbar[165]=armZZbar[3]*armZZbar[165];
   armZZbar[217]=39*armZZbar[46];
   armZZbar[218]=225./4.*armZZbar[48];
   armZZbar[221]= - 5*armZZbar[44];
   armZZbar[222]=armZZbar[221] + armZZbar[218] - 1921./6. + 
   armZZbar[217];
   armZZbar[222]=armZZbar[99]*armZZbar[222];
   armZZbar[223]=4./3.*armZZbar[101];
   armZZbar[224]= - armZZbar[34] + 17./4.*armZZbar[99];
   armZZbar[224]=armZZbar[35]*armZZbar[224];
   armZZbar[225]= - armZZbar[10]*armZZbar[99];
   armZZbar[226]=5./2.*armZZbar[225];
   armZZbar[227]=1./3.*armZZbar[98] + armZZbar[171];
   armZZbar[227]=armZZbar[6]*armZZbar[227];
   armZZbar[222]=armZZbar[227] + armZZbar[226] + 5./6.*armZZbar[224] + 
   armZZbar[173] + armZZbar[223] + 1./8.*armZZbar[222] + 8./3.*
   armZZbar[98] + 9*armZZbar[34];
   armZZbar[222]=armZZbar[37]*armZZbar[222];
   armZZbar[224]=3*armZZbar[45];
   armZZbar[227]= - 2 + armZZbar[224];
   armZZbar[227]=armZZbar[3]*armZZbar[37]*armZZbar[227];
   armZZbar[228]=6*armZZbar[227];
   armZZbar[229]=1./2.*armZZbar[36];
   armZZbar[230]=armZZbar[228] - 104./3.*armZZbar[35] + armZZbar[229]
    + armZZbar[158] - 73./6. - armZZbar[46];
   armZZbar[230]=armZZbar[3]*armZZbar[230];
   armZZbar[231]= - 4./3.*armZZbar[58] - armZZbar[56];
   armZZbar[232]=7*armZZbar[55];
   armZZbar[231]= - 11*armZZbar[52] + 2*armZZbar[231] + armZZbar[232];
   armZZbar[233]= - 2*armZZbar[51];
   armZZbar[231]= - 2./3.*armZZbar[53] - 2*armZZbar[57] + 1./3.*
   armZZbar[231] + armZZbar[233];
   armZZbar[161]=6*armZZbar[161] + 11./3.*armZZbar[164];
   armZZbar[161]=armZZbar[10]*armZZbar[161];
   armZZbar[234]= - 6*armZZbar[38];
   armZZbar[235]=armZZbar[234] + 29./3.*armZZbar[15] - 2./3.*
   armZZbar[67] - 22*armZZbar[64] - 11./3.*armZZbar[59] - 73./6. - 6*
   armZZbar[72];
   armZZbar[235]=armZZbar[101]*armZZbar[235];
   armZZbar[236]=12*armZZbar[167];
   armZZbar[237]=11./3.*armZZbar[168];
   armZZbar[230]=armZZbar[236] + armZZbar[230] + armZZbar[161] + 
   armZZbar[237] + 2*armZZbar[231] + armZZbar[235];
   armZZbar[230]=MMt*armZZbar[230];
   armZZbar[231]=armZZbar[221] + armZZbar[218] + 58175./81. + 31*
   armZZbar[46];
   armZZbar[238]=armZZbar[19]*armZZbar[34];
   armZZbar[239]=armZZbar[101]*armZZbar[17];
   armZZbar[240]=1./3.*armZZbar[239];
   armZZbar[241]=35./3.*armZZbar[12];
   armZZbar[242]= - 17./9.*armZZbar[11];
   armZZbar[231]=armZZbar[240] + 5./3.*armZZbar[238] + armZZbar[242] + 
   487./6.*armZZbar[13] + 1./16.*armZZbar[231] + armZZbar[241];
   armZZbar[217]=armZZbar[221] + armZZbar[218] + 71 + armZZbar[217];
   armZZbar[217]=armZZbar[99]*armZZbar[217];
   armZZbar[218]= - 1./3.*armZZbar[101];
   armZZbar[217]=1./32.*armZZbar[217] + armZZbar[218];
   armZZbar[217]=MMZ*armZZbar[217];
   armZZbar[238]=7./4. + 1./3.*armZZbar[200];
   armZZbar[238]=armZZbar[35]*armZZbar[238];
   armZZbar[217]=armZZbar[238] + 1./2.*armZZbar[231] + armZZbar[217];
   armZZbar[217]=armZZbar[35]*armZZbar[217];
   armZZbar[231]=armZZbar[49] - 7./3. - armZZbar[61];
   armZZbar[238]=1./3.*armZZbar[59];
   armZZbar[243]= - 1./3.*armZZbar[15];
   armZZbar[244]= - 1./2.*armZZbar[64];
   armZZbar[231]=armZZbar[243] + armZZbar[143] + 5./6.*armZZbar[67] + 
   armZZbar[244] + 1./4.*armZZbar[231] + armZZbar[238];
   armZZbar[231]=armZZbar[99]*armZZbar[231];
   armZZbar[245]=armZZbar[35]*armZZbar[99]*armZZbar[44];
   armZZbar[231]=1./3.*armZZbar[225] + armZZbar[231] + 1./4.*
   armZZbar[245];
   armZZbar[231]=5./8.*MMH*armZZbar[231];
   armZZbar[245]= - 29./27. - armZZbar[177];
   armZZbar[246]=11./9.*armZZbar[7];
   armZZbar[247]= - 17*armZZbar[7];
   armZZbar[248]=65./3. + armZZbar[247];
   armZZbar[248]=armZZbar[35]*armZZbar[248];
   armZZbar[245]=1./6.*armZZbar[248] + 11./9.*armZZbar[36] + 2*
   armZZbar[245] + armZZbar[246];
   armZZbar[245]=armZZbar[39]*armZZbar[245];
   armZZbar[249]=MMZ*armZZbar[99];
   armZZbar[250]=MMZ*armZZbar[101];
   armZZbar[251]= - 2 + armZZbar[250];
   armZZbar[251]=armZZbar[35]*armZZbar[251];
   armZZbar[252]=23./24. - armZZbar[36];
   armZZbar[251]=1./3.*armZZbar[251] + armZZbar[252] + 5./8.*
   armZZbar[249];
   armZZbar[251]=armZZbar[10]*armZZbar[251];
   armZZbar[253]= - 37./27. - armZZbar[177];
   armZZbar[253]=2*armZZbar[253] + armZZbar[246];
   armZZbar[248]=1./18.*armZZbar[248] + 1./3.*armZZbar[253] + 
   armZZbar[36];
   armZZbar[248]=armZZbar[1]*armZZbar[248];
   armZZbar[253]=5./3.*armZZbar[34];
   armZZbar[254]=armZZbar[253] - 39./8.*armZZbar[99];
   armZZbar[254]=1./2.*armZZbar[254] + armZZbar[218];
   armZZbar[254]=armZZbar[35]*armZZbar[254];
   armZZbar[255]=5*armZZbar[34] + 59./8.*armZZbar[99];
   armZZbar[254]=1./6.*armZZbar[255] + armZZbar[254];
   armZZbar[254]=armZZbar[18]*armZZbar[254];
   armZZbar[255]= - 1./3.*armZZbar[28];
   armZZbar[256]=5./32.*armZZbar[61];
   armZZbar[257]= - 23./72.*armZZbar[14];
   armZZbar[258]=28*armZZbar[16];
   armZZbar[259]= - 5./12.*armZZbar[177];
   armZZbar[260]= - 7./24.*armZZbar[59];
   armZZbar[261]= - 2./3.*armZZbar[64];
   armZZbar[262]= - 15./16.*armZZbar[67];
   armZZbar[263]= - armZZbar[98]*armZZbar[70];
   armZZbar[264]=armZZbar[50]*armZZbar[98];
   armZZbar[265]=armZZbar[22]*armZZbar[98];
   armZZbar[266]=5./32.*armZZbar[44];
   armZZbar[267]= - 2*armZZbar[70] + armZZbar[42];
   armZZbar[268]=5./2.*armZZbar[22] + armZZbar[50] + 5./2.*armZZbar[41]
    + 4*armZZbar[267] - 5*armZZbar[69];
   armZZbar[268]=armZZbar[34]*armZZbar[268];
   armZZbar[269]=armZZbar[21]*armZZbar[34];
   armZZbar[270]= - 37./18.*armZZbar[12];
   armZZbar[271]=11./27.*armZZbar[11];
   armZZbar[272]=75./8.*armZZbar[192] + 13*armZZbar[69];
   armZZbar[272]=5./4.*armZZbar[40] + 3*armZZbar[272] - 5./2.*
   armZZbar[68];
   armZZbar[272]= - 73./12.*armZZbar[21] - 307./32.*armZZbar[50] + 1./4.
   *armZZbar[272] - 11./3.*armZZbar[41];
   armZZbar[272]=armZZbar[99]*armZZbar[272];
   armZZbar[273]=45./2.*armZZbar[13];
   armZZbar[274]=61./9. + armZZbar[273];
   armZZbar[274]=armZZbar[36]*armZZbar[274];
   armZZbar[275]=11*armZZbar[34];
   armZZbar[276]=37./2.*armZZbar[98] + armZZbar[275];
   armZZbar[276]=1./3.*armZZbar[276] - 225./16.*armZZbar[99];
   armZZbar[276]=armZZbar[19]*armZZbar[276];
   armZZbar[277]=armZZbar[17]*armZZbar[99];
   armZZbar[278]=5./8.*armZZbar[277];
   armZZbar[279]=armZZbar[41] - armZZbar[8] - 2*armZZbar[68];
   armZZbar[279]=1./3.*armZZbar[279] + armZZbar[17];
   armZZbar[279]=armZZbar[101]*armZZbar[279];
   armZZbar[280]=5./8.*armZZbar[15];
   armZZbar[102]=armZZbar[102] + armZZbar[230] + armZZbar[231] + 
   armZZbar[165] + armZZbar[222] + armZZbar[109] + armZZbar[254] + 
   armZZbar[251] + 1./3.*armZZbar[245] + 1./3.*armZZbar[248] + 
   armZZbar[217] + armZZbar[108] + armZZbar[279] + armZZbar[278] + 1./2.
   *armZZbar[276] + 1./2.*armZZbar[274] + armZZbar[272] + armZZbar[271]
    + 415./36.*armZZbar[13] + armZZbar[270] + 5./6.*armZZbar[269] + 
   armZZbar[280] + 1./3.*armZZbar[268] + armZZbar[266] + 11./12.*
   armZZbar[265] + 4./3.*armZZbar[264] + 2*armZZbar[263] + 287./128.*
   armZZbar[48] + 41./32.*armZZbar[46] + armZZbar[262] + armZZbar[261]
    + armZZbar[260] + armZZbar[259] + armZZbar[258] + 2383./384.*
   armZZbar[49] - 535./768.*armZZbar[66] + armZZbar[257] + 
   armZZbar[256] - 5287./384.*armZZbar[65] - 157./48.*armZZbar[62] + 
   armZZbar[255] - 1315./384.*armZZbar[63] + 5./3.*armZZbar[31] - 
   620713./31104. - 28*armZZbar[60];
   armZZbar[102]=armZZbar[133]*armZZbar[102];
   armZZbar[108]= - 7./3. - 85*armZZbar[12];
   armZZbar[108]=1./2.*armZZbar[108] + 5*armZZbar[13];
   armZZbar[108]=1./3.*armZZbar[108] + 5*armZZbar[116];
   armZZbar[108]=1./2.*armZZbar[108] + 5*armZZbar[250];
   armZZbar[109]= - 1 + 17./9.*armZZbar[35];
   armZZbar[165]=armZZbar[1]*armZZbar[109];
   armZZbar[109]=armZZbar[39]*armZZbar[109];
   armZZbar[217]=1 + 1./3.*armZZbar[114];
   armZZbar[222]=armZZbar[10]*armZZbar[217];
   armZZbar[230]=11./9.*armZZbar[39];
   armZZbar[245]=armZZbar[230] - 17./72. + armZZbar[1];
   armZZbar[245]=armZZbar[6]*armZZbar[245];
   armZZbar[108]=5./9.*armZZbar[245] + 5./18.*armZZbar[117] + 5./6.*
   armZZbar[222] + 11./9.*armZZbar[109] + 1./18.*armZZbar[108] + 
   armZZbar[165];
   armZZbar[108]=armZZbar[6]*armZZbar[108];
   armZZbar[109]= - 1./2. + 17*armZZbar[46];
   armZZbar[109]=armZZbar[122] + armZZbar[36] + 1./3.*armZZbar[109] + 
   armZZbar[158];
   armZZbar[109]=1./2.*armZZbar[109] + 3*armZZbar[227];
   armZZbar[109]=armZZbar[3]*armZZbar[109];
   armZZbar[146]=armZZbar[158] + 47./18.*armZZbar[15] + 17./9.*
   armZZbar[67] + 7./3.*armZZbar[64] + 7./18.*armZZbar[59] + 161./36.
    + armZZbar[146];
   armZZbar[146]=armZZbar[101]*armZZbar[146];
   armZZbar[163]=armZZbar[163] + 7./18.*armZZbar[168];
   armZZbar[163]=armZZbar[10]*armZZbar[163];
   armZZbar[165]= - 22*armZZbar[55] + 17*armZZbar[53];
   armZZbar[164]=7./18.*armZZbar[164];
   armZZbar[109]=6*armZZbar[167] + armZZbar[109] + armZZbar[163] + 
   armZZbar[164] + 2./9.*armZZbar[165] + armZZbar[146];
   armZZbar[109]=MMt*armZZbar[109];
   armZZbar[146]= - 5*armZZbar[8] - 17*armZZbar[41];
   armZZbar[163]=11*armZZbar[21];
   armZZbar[146]= - 44./3.*MMZ - 11./6.*armZZbar[17] + armZZbar[163] + 
   11./6.*armZZbar[20] + 1./2.*armZZbar[146] + 17*armZZbar[50];
   armZZbar[165]=1./9.*armZZbar[17] + MMZ;
   armZZbar[167]=armZZbar[35]*armZZbar[165];
   armZZbar[168]=11./9. + armZZbar[162];
   armZZbar[168]=armZZbar[18]*armZZbar[168];
   armZZbar[222]=armZZbar[165] + armZZbar[18];
   armZZbar[222]=armZZbar[6]*armZZbar[222];
   armZZbar[227]=15./2. + 17./3.*armZZbar[46];
   armZZbar[227]=armZZbar[37]*armZZbar[227];
   armZZbar[146]=armZZbar[227] + 5./4.*armZZbar[222] + armZZbar[168] + 
   1./3.*armZZbar[146] + 17./4.*armZZbar[167];
   armZZbar[146]=armZZbar[3]*armZZbar[146];
   armZZbar[167]=103./3. - 1./2.*armZZbar[63];
   armZZbar[142]=armZZbar[142] + 1./4.*armZZbar[167] + 25./3.*
   armZZbar[62];
   armZZbar[142]= - 185./72.*armZZbar[49] + 1./3.*armZZbar[142] + 1./16.
   *armZZbar[66];
   armZZbar[167]=1./2.*armZZbar[67];
   armZZbar[142]=armZZbar[243] - 1./12.*armZZbar[44] - 1./48.*
   armZZbar[48] + 25./36.*armZZbar[46] + armZZbar[167] + 1./2.*
   armZZbar[142] + armZZbar[238];
   armZZbar[142]=armZZbar[99]*armZZbar[142];
   armZZbar[168]=257*armZZbar[55] + 17*armZZbar[25];
   armZZbar[222]=1./3.*armZZbar[49] + 1./4. - 1./3.*armZZbar[62];
   armZZbar[222]=armZZbar[34]*armZZbar[222];
   armZZbar[168]=1./9.*armZZbar[168] + 107./4.*armZZbar[222];
   armZZbar[222]= - 1./2.*armZZbar[14] + 1./6.*armZZbar[28] + 13./8. - 
   1./3.*armZZbar[30];
   armZZbar[222]= - 11./3.*armZZbar[15] + 17*armZZbar[64] + 5*
   armZZbar[222] + 17./6.*armZZbar[59];
   armZZbar[222]=armZZbar[101]*armZZbar[222];
   armZZbar[142]=1./3.*armZZbar[222] + 1./9.*armZZbar[168] + 25./16.*
   armZZbar[142];
   armZZbar[142]=MMZ*armZZbar[142];
   armZZbar[168]= - armZZbar[49] + 7./3. + armZZbar[61];
   armZZbar[222]=1./4.*armZZbar[44];
   armZZbar[227]=1./3.*armZZbar[15];
   armZZbar[136]=armZZbar[227] + armZZbar[222] - 5./6.*armZZbar[67] + 1.
   /2.*armZZbar[64] + 1./4.*armZZbar[168] + armZZbar[136];
   armZZbar[136]=armZZbar[99]*armZZbar[136];
   armZZbar[168]=1./3.*armZZbar[175];
   armZZbar[136]=armZZbar[168] + armZZbar[136] + 1./4.*armZZbar[204];
   armZZbar[136]=MMH*armZZbar[136];
   armZZbar[195]=armZZbar[196] + armZZbar[195] + 1./8.*armZZbar[192] - 
   25./3.*armZZbar[69];
   armZZbar[195]=47./36.*armZZbar[21] + 941./288.*armZZbar[50] + 1./4.*
   armZZbar[195] + 7./9.*armZZbar[41];
   armZZbar[195]=armZZbar[99]*armZZbar[195];
   armZZbar[196]= - 4841./9. - 217*armZZbar[46];
   armZZbar[204]=25*armZZbar[44];
   armZZbar[238]=25./4.*armZZbar[48];
   armZZbar[196]=armZZbar[204] + 1./3.*armZZbar[196] + armZZbar[238];
   armZZbar[245]= - 289./9.*armZZbar[12];
   armZZbar[196]=1./8.*armZZbar[196] + armZZbar[245];
   armZZbar[196]=17./3.*armZZbar[116] + 1./2.*armZZbar[196] + 17./9.*
   armZZbar[13];
   armZZbar[248]=1./4.*armZZbar[48];
   armZZbar[181]=armZZbar[44] + 25./3.*armZZbar[181] + armZZbar[248];
   armZZbar[181]=armZZbar[99]*armZZbar[181];
   armZZbar[181]=25./32.*armZZbar[181] + 17./3.*armZZbar[101];
   armZZbar[181]=MMZ*armZZbar[181];
   armZZbar[181]= - 1285./108.*armZZbar[35] + 1./2.*armZZbar[196] + 
   armZZbar[181];
   armZZbar[181]=armZZbar[35]*armZZbar[181];
   armZZbar[196]=841./6. - 25*armZZbar[46];
   armZZbar[196]=armZZbar[44] + 1./3.*armZZbar[196] + armZZbar[248];
   armZZbar[196]=armZZbar[99]*armZZbar[196];
   armZZbar[248]= - 107./3.*armZZbar[34];
   armZZbar[196]=armZZbar[248] + 25./8.*armZZbar[196];
   armZZbar[254]=107./3.*armZZbar[34] - 625./4.*armZZbar[99];
   armZZbar[254]=armZZbar[35]*armZZbar[254];
   armZZbar[196]=25./4.*armZZbar[175] + 1./12.*armZZbar[254] + 1./2.*
   armZZbar[196] - 34./3.*armZZbar[101];
   armZZbar[196]=armZZbar[37]*armZZbar[196];
   armZZbar[254]=87691./192. - 85*armZZbar[31];
   armZZbar[254]= - 115./72.*armZZbar[14] - 25./32.*armZZbar[61] + 
   10721./432.*armZZbar[62] - 5./3.*armZZbar[28] + 1./27.*armZZbar[254]
    - 25./128.*armZZbar[63];
   armZZbar[254]= - 61./72.*armZZbar[59] + 85./324.*armZZbar[177] - 
   82393./10368.*armZZbar[49] + 1./3.*armZZbar[254] + 25./256.*
   armZZbar[66];
   armZZbar[263]= - 49./9. + 25*armZZbar[200];
   armZZbar[264]=armZZbar[35]*armZZbar[217];
   armZZbar[263]=1./8.*armZZbar[263] + 17*armZZbar[264];
   armZZbar[263]=armZZbar[10]*armZZbar[263];
   armZZbar[264]= - 107*armZZbar[34] - 925./8.*armZZbar[99];
   armZZbar[248]=armZZbar[248] + 625./8.*armZZbar[99];
   armZZbar[248]=1./2.*armZZbar[248] + 17*armZZbar[101];
   armZZbar[248]=armZZbar[35]*armZZbar[248];
   armZZbar[248]=1./6.*armZZbar[264] + armZZbar[248];
   armZZbar[248]=armZZbar[18]*armZZbar[248];
   armZZbar[264]= - armZZbar[6]*armZZbar[12];
   armZZbar[268]= - armZZbar[35]*armZZbar[12];
   armZZbar[264]=5./2.*armZZbar[264] + 11./3.*armZZbar[12] + 17./2.*
   armZZbar[268];
   armZZbar[272]=pow(armZZbar[2],2);
   armZZbar[264]=armZZbar[272]*armZZbar[264];
   armZZbar[274]= - armZZbar[50] + armZZbar[69] + armZZbar[112];
   armZZbar[274]=armZZbar[34]*armZZbar[274];
   armZZbar[276]= - armZZbar[19]*armZZbar[99];
   armZZbar[163]= - 20*armZZbar[17] + armZZbar[163] - 11./2.*
   armZZbar[20] - 17./2.*armZZbar[41] - 5./2.*armZZbar[8] + 17*
   armZZbar[68];
   armZZbar[163]=armZZbar[101]*armZZbar[163];
   armZZbar[281]= - 17./3.*armZZbar[35] + 22./9. + 5*armZZbar[177];
   armZZbar[282]=armZZbar[1]*armZZbar[281];
   armZZbar[281]=armZZbar[39]*armZZbar[281];
   armZZbar[108]=1./108.*armZZbar[264] + armZZbar[109] + 25./48.*
   armZZbar[136] + armZZbar[146] + 1./3.*armZZbar[196] + armZZbar[108]
    + 1./18.*armZZbar[248] + 1./6.*armZZbar[263] + 11./81.*
   armZZbar[281] + 1./9.*armZZbar[282] + 1./6.*armZZbar[181] + 
   armZZbar[142] + 1./9.*armZZbar[163] + 25./48.*armZZbar[213] + 25./
   192.*armZZbar[276] + 25./6.*armZZbar[195] - 11./162.*armZZbar[13] + 
   187./324.*armZZbar[12] + 107./108.*armZZbar[209] + 101./144.*
   armZZbar[15] + 107./54.*armZZbar[274] - 25./192.*armZZbar[44] - 25./
   768.*armZZbar[48] - 2749./1728.*armZZbar[46] + 25./32.*armZZbar[67]
    + 1./2.*armZZbar[254] + 17./9.*armZZbar[64];
   armZZbar[108]=armZZbar[272]*armZZbar[108];
   armZZbar[109]=5*armZZbar[19];
   armZZbar[136]=armZZbar[109] + armZZbar[17];
   armZZbar[142]=armZZbar[186] + 1./2.*armZZbar[136] + 11./3.*MMZ;
   armZZbar[142]=armZZbar[6]*armZZbar[142];
   armZZbar[146]= - 11./9.*armZZbar[19];
   armZZbar[163]=20./27. - armZZbar[36];
   armZZbar[163]=MMZ*armZZbar[163];
   armZZbar[181]=17*armZZbar[19];
   armZZbar[195]=armZZbar[181] - armZZbar[17];
   armZZbar[195]=1./2.*armZZbar[195] - 22./3.*MMZ;
   armZZbar[195]=armZZbar[35]*armZZbar[195];
   armZZbar[196]= - 53./9. + armZZbar[220];
   armZZbar[196]=armZZbar[37]*armZZbar[196];
   armZZbar[142]=armZZbar[196] + 1./3.*armZZbar[142] + armZZbar[219] + 
   1./3.*armZZbar[195] + armZZbar[163] + armZZbar[203] + armZZbar[146];
   armZZbar[142]=armZZbar[3]*armZZbar[142];
   armZZbar[163]=11*armZZbar[12];
   armZZbar[107]=armZZbar[107] - 5./2.*armZZbar[11] + 89./12.*
   armZZbar[13] + 283./18. + armZZbar[163];
   armZZbar[107]=armZZbar[134] + armZZbar[250] + 1./3.*armZZbar[107] + 
   armZZbar[185];
   armZZbar[195]=433./9. + armZZbar[179];
   armZZbar[195]= - 440./9.*armZZbar[35] + 1./2.*armZZbar[195] + 
   armZZbar[189];
   armZZbar[195]=armZZbar[39]*armZZbar[195];
   armZZbar[196]=347./3. + armZZbar[179];
   armZZbar[196]= - 40./9.*armZZbar[35] + 1./54.*armZZbar[196] - 
   armZZbar[36];
   armZZbar[196]=armZZbar[1]*armZZbar[196];
   armZZbar[203]= - 44./9.*armZZbar[39] + 55./36. - 4*armZZbar[1];
   armZZbar[203]=armZZbar[6]*armZZbar[203];
   armZZbar[107]=1./9.*armZZbar[203] + armZZbar[199] + armZZbar[194] + 
   1./9.*armZZbar[195] + 1./3.*armZZbar[107] + armZZbar[196];
   armZZbar[107]=armZZbar[6]*armZZbar[107];
   armZZbar[194]= - 3431./2. + 275*armZZbar[46];
   armZZbar[195]=7./4.*armZZbar[48];
   armZZbar[194]=1./3.*armZZbar[194] + armZZbar[195];
   armZZbar[194]=1./3.*armZZbar[194] - armZZbar[44];
   armZZbar[194]=armZZbar[99]*armZZbar[194];
   armZZbar[196]=173*armZZbar[34];
   armZZbar[199]=armZZbar[98] + armZZbar[196];
   armZZbar[203]= - 173*armZZbar[34];
   armZZbar[219]=armZZbar[203] + 3325./4.*armZZbar[99];
   armZZbar[219]=armZZbar[35]*armZZbar[219];
   armZZbar[248]= - armZZbar[6]*armZZbar[98];
   armZZbar[194]=1./9.*armZZbar[248] + armZZbar[226] + 1./54.*
   armZZbar[219] + armZZbar[223] + 1./9.*armZZbar[199] + 5./8.*
   armZZbar[194];
   armZZbar[194]=armZZbar[37]*armZZbar[194];
   armZZbar[199]= - 68./27.*armZZbar[58] - armZZbar[56];
   armZZbar[219]=1./3.*armZZbar[57];
   armZZbar[223]= - 2*armZZbar[53];
   armZZbar[199]=armZZbar[223] + armZZbar[219] - armZZbar[52] + 2*
   armZZbar[199] + armZZbar[232];
   armZZbar[226]=armZZbar[228] - 152./9.*armZZbar[35] + armZZbar[229]
    + armZZbar[158] + 101./18. - armZZbar[46];
   armZZbar[226]=armZZbar[3]*armZZbar[226];
   armZZbar[161]=armZZbar[236] + armZZbar[226] + armZZbar[161] + 
   armZZbar[237] + 2./3.*armZZbar[199] + armZZbar[235];
   armZZbar[161]=MMt*armZZbar[161];
   armZZbar[199]=41407./9. + 1303*armZZbar[46];
   armZZbar[226]=35./4.*armZZbar[48];
   armZZbar[199]=1./3.*armZZbar[199] + armZZbar[226];
   armZZbar[199]=1./3.*armZZbar[199] + armZZbar[221];
   armZZbar[228]=73./9.*armZZbar[12];
   armZZbar[184]=armZZbar[240] + 11./3.*armZZbar[184] + armZZbar[242]
    + 263./54.*armZZbar[13] + 1./16.*armZZbar[199] + armZZbar[228];
   armZZbar[199]=29 + 55./3.*armZZbar[46];
   armZZbar[195]=5*armZZbar[199] + armZZbar[195];
   armZZbar[195]=1./3.*armZZbar[195] - armZZbar[44];
   armZZbar[195]=armZZbar[99]*armZZbar[195];
   armZZbar[195]=5./32.*armZZbar[195] + armZZbar[218];
   armZZbar[195]=MMZ*armZZbar[195];
   armZZbar[199]=203./12. + 5*armZZbar[200];
   armZZbar[199]=armZZbar[35]*armZZbar[199];
   armZZbar[184]=5./27.*armZZbar[199] + 1./2.*armZZbar[184] + 
   armZZbar[195];
   armZZbar[184]=armZZbar[35]*armZZbar[184];
   armZZbar[195]=5*armZZbar[27];
   armZZbar[199]= - 22*armZZbar[25] - 526*armZZbar[55] + 68*
   armZZbar[56] + armZZbar[195];
   armZZbar[200]= - 1./3.*armZZbar[49] - 1./4. + 1./3.*armZZbar[62];
   armZZbar[200]=armZZbar[34]*armZZbar[200];
   armZZbar[199]=173./2.*armZZbar[200] + 1./9.*armZZbar[199] + 
   armZZbar[57];
   armZZbar[200]= - 2203./3. - 7*armZZbar[63];
   armZZbar[200]=5./8.*armZZbar[65] + 1./8.*armZZbar[200] - 275./3.*
   armZZbar[62];
   armZZbar[218]=1./2.*armZZbar[61];
   armZZbar[200]=5915./216.*armZZbar[49] + 7./16.*armZZbar[66] + 1./3.*
   armZZbar[200] + armZZbar[218];
   armZZbar[200]=armZZbar[15] + armZZbar[222] - 7./48.*armZZbar[48] - 
   275./36.*armZZbar[46] - 3./2.*armZZbar[67] + 1./2.*armZZbar[200] - 
   armZZbar[59];
   armZZbar[200]=armZZbar[99]*armZZbar[200];
   armZZbar[159]=armZZbar[159] + 1./9.*armZZbar[199] + 5./8.*
   armZZbar[200];
   armZZbar[159]=MMZ*armZZbar[159];
   armZZbar[199]= - 85./27. - armZZbar[177];
   armZZbar[199]=2*armZZbar[199] + armZZbar[246];
   armZZbar[200]=257./3. + armZZbar[247];
   armZZbar[200]=armZZbar[35]*armZZbar[200];
   armZZbar[199]=1./18.*armZZbar[200] + 1./3.*armZZbar[199] + 
   armZZbar[36];
   armZZbar[199]=armZZbar[1]*armZZbar[199];
   armZZbar[200]= - 31./9. - armZZbar[177];
   armZZbar[200]=armZZbar[36] + 2./3.*armZZbar[200] + armZZbar[7];
   armZZbar[222]=1033./27. + armZZbar[247];
   armZZbar[222]=armZZbar[35]*armZZbar[222];
   armZZbar[200]=11./3.*armZZbar[200] + 1./2.*armZZbar[222];
   armZZbar[200]=armZZbar[39]*armZZbar[200];
   armZZbar[196]=armZZbar[196] + 1555./8.*armZZbar[99];
   armZZbar[222]=173./3.*armZZbar[34] - 1375./8.*armZZbar[99];
   armZZbar[222]=1./6.*armZZbar[222] - armZZbar[101];
   armZZbar[222]=armZZbar[35]*armZZbar[222];
   armZZbar[196]=1./18.*armZZbar[196] + armZZbar[222];
   armZZbar[196]=armZZbar[18]*armZZbar[196];
   armZZbar[222]= - 35939./128. + 55*armZZbar[31];
   armZZbar[222]= - 565./384.*armZZbar[65] - 9359./432.*armZZbar[62] - 
   armZZbar[28] + 1./27.*armZZbar[222] - 163./128.*armZZbar[63];
   armZZbar[232]=1./2.*armZZbar[41];
   armZZbar[235]=armZZbar[50] - armZZbar[69] + armZZbar[232];
   armZZbar[235]=173./9.*armZZbar[235] - 11./2.*armZZbar[22];
   armZZbar[235]=armZZbar[34]*armZZbar[235];
   armZZbar[236]=7./8.*armZZbar[192] + 275./3.*armZZbar[69];
   armZZbar[236]=1./4.*armZZbar[40] + 1./3.*armZZbar[236] - 1./2.*
   armZZbar[68];
   armZZbar[236]= - 541./108.*armZZbar[21] - 9773./864.*armZZbar[50] + 
   1./4.*armZZbar[236] - 71./27.*armZZbar[41];
   armZZbar[236]=armZZbar[99]*armZZbar[236];
   armZZbar[240]= - 1./2.*armZZbar[13];
   armZZbar[246]=1./3. + armZZbar[240];
   armZZbar[248]=armZZbar[36]*armZZbar[246];
   armZZbar[254]= - 35./16.*armZZbar[99] - 1./6.*armZZbar[98] + 
   armZZbar[275];
   armZZbar[254]=armZZbar[19]*armZZbar[254];
   armZZbar[107]=armZZbar[108] + armZZbar[161] + armZZbar[231] + 
   armZZbar[142] + armZZbar[194] + armZZbar[107] + 1./3.*armZZbar[196]
    + armZZbar[251] + 1./9.*armZZbar[200] + 1./3.*armZZbar[199] + 
   armZZbar[184] + armZZbar[159] + armZZbar[279] + armZZbar[278] + 1./6.
   *armZZbar[254] + 1./6.*armZZbar[248] + 5*armZZbar[236] + 
   armZZbar[271] - 289./324.*armZZbar[13] - 95./54.*armZZbar[12] + 173./
   54.*armZZbar[269] + armZZbar[280] + 1./3.*armZZbar[235] + 
   armZZbar[266] + 1./36.*armZZbar[265] + 1./9.*armZZbar[207] - 35./384.
   *armZZbar[48] + 2131./864.*armZZbar[46] + armZZbar[262] + 
   armZZbar[261] + armZZbar[260] - 55./324.*armZZbar[177] + 68053./
   10368.*armZZbar[49] + 59./2304.*armZZbar[66] + armZZbar[257] + 1./3.
   *armZZbar[222] + armZZbar[256];
   armZZbar[107]=armZZbar[272]*armZZbar[107];
   armZZbar[108]=pow(CW,2);
   armZZbar[142]=armZZbar[108]*armZZbar[144];
   armZZbar[142]=armZZbar[144] + 1./3.*armZZbar[142];
   armZZbar[142]=MMZ*armZZbar[142];
   armZZbar[144]=64*armZZbar[55] - 16*armZZbar[56] - armZZbar[27];
   armZZbar[159]= - 3*armZZbar[52];
   armZZbar[144]= - 17./9.*armZZbar[57] + 16./81.*armZZbar[25] - 5*
   armZZbar[51] + 4./81.*armZZbar[144] + armZZbar[159];
   armZZbar[161]= - 2./3.*armZZbar[57] - armZZbar[52] + armZZbar[233];
   armZZbar[184]=armZZbar[66] + 4 + armZZbar[65];
   armZZbar[184]=armZZbar[98]*armZZbar[184];
   armZZbar[161]=4*armZZbar[161] + armZZbar[184];
   armZZbar[161]=armZZbar[108]*armZZbar[161];
   armZZbar[184]=1 + armZZbar[65];
   armZZbar[184]=armZZbar[108]*armZZbar[184];
   armZZbar[184]=4*armZZbar[184] + 100./9.*armZZbar[49] - 2*
   armZZbar[65] - 64./9.*armZZbar[62] + 19./3. - 4*armZZbar[63];
   armZZbar[184]=armZZbar[34]*armZZbar[184];
   armZZbar[194]= - 1 - armZZbar[65];
   armZZbar[194]=armZZbar[108]*armZZbar[194];
   armZZbar[196]= - 5*armZZbar[65] + 64./3.*armZZbar[62] + 49./3. + 4*
   armZZbar[63];
   armZZbar[194]=8./3.*armZZbar[194] + 4./3.*armZZbar[48] + 32./9.*
   armZZbar[46] - 196./27.*armZZbar[49] + 1./3.*armZZbar[196] - 2*
   armZZbar[66];
   armZZbar[194]=armZZbar[99]*armZZbar[194];
   armZZbar[196]=8./3.*armZZbar[66] + 43./4. + 10./3.*armZZbar[65];
   armZZbar[196]=armZZbar[98]*armZZbar[196];
   armZZbar[142]=armZZbar[142] + 1./3.*armZZbar[156] + 2*armZZbar[194]
    + 2./3.*armZZbar[184] + 1./3.*armZZbar[154] + 1./3.*armZZbar[161]
    + 1./3.*armZZbar[153] + 2*armZZbar[144] + 1./3.*armZZbar[196];
   armZZbar[142]=MMZ*armZZbar[142];
   armZZbar[144]= - 11./9. - armZZbar[108];
   armZZbar[153]= - 8*armZZbar[108];
   armZZbar[154]= - 79./3. + armZZbar[153];
   armZZbar[154]=armZZbar[13]*armZZbar[154];
   armZZbar[144]=8*armZZbar[144] + armZZbar[154];
   armZZbar[110]=MMZ*armZZbar[110];
   armZZbar[154]=2*armZZbar[35];
   armZZbar[156]= - 1 + armZZbar[154];
   armZZbar[161]=armZZbar[1]*armZZbar[156];
   armZZbar[156]=armZZbar[39]*armZZbar[156];
   armZZbar[184]=40./3.*armZZbar[39] - 5./3. + 8*armZZbar[1];
   armZZbar[184]=armZZbar[6]*armZZbar[184];
   armZZbar[194]=8./3.*armZZbar[35];
   armZZbar[110]=4./9.*armZZbar[184] + 320./27.*armZZbar[156] + 64./9.*
   armZZbar[161] + armZZbar[194] + 4./3.*armZZbar[144] + armZZbar[110];
   armZZbar[110]=armZZbar[6]*armZZbar[110];
   armZZbar[144]= - 59./9. + armZZbar[220];
   armZZbar[156]= - 113./3. - 16*armZZbar[108];
   armZZbar[156]=armZZbar[13]*armZZbar[156];
   armZZbar[161]= - 2 - armZZbar[46];
   armZZbar[161]=8./3.*armZZbar[161] - armZZbar[48];
   armZZbar[161]=MMZ*armZZbar[99]*armZZbar[161];
   armZZbar[184]= - 2./3. + armZZbar[249];
   armZZbar[184]=armZZbar[35]*armZZbar[184];
   armZZbar[144]=8./9.*armZZbar[184] + armZZbar[161] + 1./3.*
   armZZbar[156] - 16./3.*armZZbar[108] + 4./3.*armZZbar[144] - 
   armZZbar[48];
   armZZbar[144]=armZZbar[35]*armZZbar[144];
   armZZbar[148]=armZZbar[148] + armZZbar[134];
   armZZbar[148]=armZZbar[35]*armZZbar[148];
   armZZbar[156]= - 1./2.*armZZbar[6];
   armZZbar[149]=armZZbar[149] + armZZbar[156];
   armZZbar[149]=armZZbar[6]*armZZbar[149];
   armZZbar[161]= - 1./3.*armZZbar[36];
   armZZbar[148]=1./2.*armZZbar[149] + armZZbar[161] + 1./2.*
   armZZbar[148];
   armZZbar[148]=armZZbar[133]*armZZbar[148];
   armZZbar[149]= - 11./9. - armZZbar[36];
   armZZbar[184]=2./3.*armZZbar[35];
   armZZbar[196]=1./2.*armZZbar[149] + armZZbar[184];
   armZZbar[196]=armZZbar[35]*armZZbar[196];
   armZZbar[199]= - 7./9. - armZZbar[36];
   armZZbar[199]=1./3.*armZZbar[6] + 1./2.*armZZbar[199] + armZZbar[35]
   ;
   armZZbar[199]=armZZbar[6]*armZZbar[199];
   armZZbar[148]=armZZbar[148] + armZZbar[199] + armZZbar[196] + 
   armZZbar[188] + 1./9. + 1./4.*armZZbar[177];
   armZZbar[148]=armZZbar[133]*armZZbar[148];
   armZZbar[196]= - 5*armZZbar[36];
   armZZbar[199]= - 17*armZZbar[35];
   armZZbar[200]= - 7./2.*armZZbar[6] + armZZbar[199] + 29./3. + 
   armZZbar[196];
   armZZbar[200]=armZZbar[6]*armZZbar[200];
   armZZbar[220]= - 17*armZZbar[36];
   armZZbar[222]= - 31./2.*armZZbar[35] + 65./3. + armZZbar[220];
   armZZbar[222]=armZZbar[35]*armZZbar[222];
   armZZbar[200]=1./2.*armZZbar[200] + 1./2.*armZZbar[222] + 11./3.*
   armZZbar[36] - 29./9. - 1./2.*armZZbar[177];
   armZZbar[148]=1./9.*armZZbar[200] + armZZbar[148];
   armZZbar[148]=armZZbar[133]*armZZbar[148];
   armZZbar[162]= - 11./3. + armZZbar[162];
   armZZbar[162]=armZZbar[35]*armZZbar[162];
   armZZbar[200]= - 11./3. + 17./2.*armZZbar[35];
   armZZbar[222]=armZZbar[200] + 5./4.*armZZbar[6];
   armZZbar[222]=armZZbar[6]*armZZbar[222];
   armZZbar[162]=5*armZZbar[222] + 17*armZZbar[162] + 121./9. + 25./4.*
   armZZbar[177];
   armZZbar[162]=armZZbar[272]*armZZbar[162];
   armZZbar[222]= - 341./9. + 7./2.*armZZbar[177];
   armZZbar[222]=1./3.*armZZbar[222] + 11*armZZbar[36];
   armZZbar[231]= - 791./18.*armZZbar[35] + 1033./27. + armZZbar[220];
   armZZbar[231]=armZZbar[35]*armZZbar[231];
   armZZbar[235]=1./18.*armZZbar[6] - 89./9.*armZZbar[35] + 133./27. + 
   armZZbar[196];
   armZZbar[235]=armZZbar[6]*armZZbar[235];
   armZZbar[162]=1./9.*armZZbar[162] + 1./2.*armZZbar[235] + 1./3.*
   armZZbar[222] + 1./2.*armZZbar[231];
   armZZbar[162]=armZZbar[272]*armZZbar[162];
   armZZbar[222]= - 5./3. + armZZbar[154];
   armZZbar[222]=armZZbar[35]*armZZbar[222];
   armZZbar[231]=4*armZZbar[35];
   armZZbar[235]= - 5./3. + armZZbar[231];
   armZZbar[236]=2*armZZbar[235] + armZZbar[6];
   armZZbar[236]=armZZbar[6]*armZZbar[236];
   armZZbar[222]=armZZbar[236] + 25./9. + 8*armZZbar[222];
   armZZbar[162]=16./9.*armZZbar[222] + armZZbar[162];
   armZZbar[148]=1./9.*armZZbar[162] + armZZbar[148];
   armZZbar[148]=armZZbar[33]*armZZbar[148];
   armZZbar[162]= - 4*armZZbar[46];
   armZZbar[222]=85./3. + armZZbar[162];
   armZZbar[222]=2./3.*armZZbar[222] - armZZbar[48];
   armZZbar[222]=armZZbar[99]*armZZbar[222];
   armZZbar[236]=armZZbar[34] - 4*armZZbar[99];
   armZZbar[236]=armZZbar[35]*armZZbar[236];
   armZZbar[173]=armZZbar[6]*armZZbar[173];
   armZZbar[171]=1./3.*armZZbar[173] + 128./27.*armZZbar[236] + 1./3.*
   armZZbar[171] + 32./3.*armZZbar[222] - armZZbar[98] - 328./9.*
   armZZbar[34];
   armZZbar[171]=armZZbar[37]*armZZbar[171];
   armZZbar[173]=167./3. + 29*armZZbar[65];
   armZZbar[173]=armZZbar[108]*armZZbar[173];
   armZZbar[222]= - armZZbar[34] - armZZbar[99];
   armZZbar[236]= - 1./3.*armZZbar[34] + armZZbar[99];
   armZZbar[236]=armZZbar[35]*armZZbar[236];
   armZZbar[222]=1./3.*armZZbar[222] + armZZbar[236];
   armZZbar[222]=armZZbar[18]*armZZbar[222];
   armZZbar[236]=armZZbar[35]*MMZ;
   armZZbar[248]=armZZbar[6]*MMZ;
   armZZbar[236]=8*armZZbar[37] + armZZbar[248] - 5./3.*MMZ + 4*
   armZZbar[236];
   armZZbar[236]=armZZbar[3]*armZZbar[236];
   armZZbar[248]=1 + armZZbar[35];
   armZZbar[249]=armZZbar[3]*armZZbar[248];
   armZZbar[251]=2./3.*armZZbar[57];
   armZZbar[249]=64./9.*armZZbar[249] + armZZbar[251] - 2./3.*
   armZZbar[51] + 256./81.*armZZbar[58] + armZZbar[52];
   armZZbar[249]=MMt*armZZbar[249];
   armZZbar[254]= - 23./9.*armZZbar[50] - 16./9.*armZZbar[41] + 32./9.*
   armZZbar[69] + 2*armZZbar[70] - armZZbar[42];
   armZZbar[254]=armZZbar[34]*armZZbar[254];
   armZZbar[256]=29./3. + 10*armZZbar[108];
   armZZbar[260]=armZZbar[13]*armZZbar[256];
   armZZbar[261]=64./9.*armZZbar[21] + 14*armZZbar[50] + 32./9.*
   armZZbar[41] + armZZbar[267] - 32./3.*armZZbar[69];
   armZZbar[261]=armZZbar[99]*armZZbar[261];
   armZZbar[262]= - 1 - armZZbar[13];
   armZZbar[263]=armZZbar[36]*armZZbar[262];
   armZZbar[264]=32./3.*armZZbar[99] - armZZbar[98] - 16./3.*
   armZZbar[34];
   armZZbar[264]=armZZbar[19]*armZZbar[264];
   armZZbar[265]= - 4*armZZbar[35];
   armZZbar[266]=5./3. + armZZbar[265];
   armZZbar[267]=armZZbar[1]*armZZbar[266];
   armZZbar[266]=armZZbar[39]*armZZbar[266];
   armZZbar[102]=armZZbar[148] + armZZbar[102] + armZZbar[107] + 2*
   armZZbar[249] + 16./9.*armZZbar[236] + armZZbar[171] + 1./3.*
   armZZbar[110] + 128./9.*armZZbar[222] + 160./243.*armZZbar[266] + 32.
   /81.*armZZbar[267] + 8./3.*armZZbar[144] + armZZbar[142] + 
   armZZbar[264] + 6*armZZbar[263] + 16./3.*armZZbar[261] + 16./27.*
   armZZbar[260] + 128./27.*armZZbar[209] + 8./3.*armZZbar[254] + 1./3.
   *armZZbar[208] + 2./9.*armZZbar[173] + 1./3.*armZZbar[207] + 2./3.*
   armZZbar[206] - 8./3.*armZZbar[48] - 64./27.*armZZbar[46] + 20./81.*
   armZZbar[177] - 10*armZZbar[16] - 520./81.*armZZbar[49] + 5./3.*
   armZZbar[66] + 107./9.*armZZbar[65] + 640./81.*armZZbar[62] + 8./9.*
   armZZbar[63] - 80./81.*armZZbar[31] + 2687./162. + 10*armZZbar[60];
   armZZbar[102]=armZZbar[33]*armZZbar[102];
   armZZbar[107]= - 2*armZZbar[22];
   armZZbar[110]=armZZbar[107] - 8*armZZbar[50] + 2*armZZbar[42] + 3*
   armZZbar[40];
   armZZbar[121]= - 1 + armZZbar[121];
   armZZbar[121]=armZZbar[19]*armZZbar[121];
   armZZbar[142]= - 7./3. + armZZbar[158];
   armZZbar[142]=armZZbar[17]*armZZbar[142];
   armZZbar[144]= - 1 - armZZbar[38];
   armZZbar[148]=armZZbar[18]*armZZbar[144];
   armZZbar[171]= - 8*armZZbar[44] - 3 - 2*armZZbar[48];
   armZZbar[171]=armZZbar[37]*armZZbar[171];
   armZZbar[173]= - armZZbar[35]*armZZbar[17];
   armZZbar[110]=3*armZZbar[171] + 6*armZZbar[148] + 11./3.*
   armZZbar[173] + armZZbar[142] + 3*armZZbar[121] + 3*armZZbar[110] - 
   38./3.*armZZbar[20];
   armZZbar[110]=armZZbar[3]*armZZbar[110];
   armZZbar[121]=armZZbar[1]*armZZbar[36];
   armZZbar[142]=armZZbar[39]*armZZbar[36];
   armZZbar[148]=1./3.*armZZbar[142];
   armZZbar[171]=armZZbar[148] + armZZbar[121] - armZZbar[13] + 
   armZZbar[154];
   armZZbar[171]=1./3.*armZZbar[6]*armZZbar[171];
   armZZbar[206]=15*armZZbar[60];
   armZZbar[207]=11./2.*armZZbar[63];
   armZZbar[208]=armZZbar[207] - 7963./324. + armZZbar[206];
   armZZbar[209]= - armZZbar[17]*armZZbar[38];
   armZZbar[222]=armZZbar[209] - armZZbar[20] + armZZbar[40] + 
   armZZbar[193];
   armZZbar[222]=3*armZZbar[101]*armZZbar[222];
   armZZbar[236]= - 4*armZZbar[44];
   armZZbar[249]=3 - armZZbar[48];
   armZZbar[254]=armZZbar[249] + armZZbar[236];
   armZZbar[254]=armZZbar[3]*armZZbar[254];
   armZZbar[233]=2*armZZbar[57] - 11*armZZbar[54] + armZZbar[52] + 
   armZZbar[233];
   armZZbar[233]=3*armZZbar[254] + 1./3.*armZZbar[233] + 6*
   armZZbar[101];
   armZZbar[233]=MMt*armZZbar[233];
   armZZbar[254]= - 4*armZZbar[72];
   armZZbar[260]= - 2*armZZbar[62];
   armZZbar[261]=11./3.*armZZbar[61];
   armZZbar[263]= - 73./24.*armZZbar[49];
   armZZbar[264]= - 15./4.*armZZbar[16];
   armZZbar[266]= - 11./3.*armZZbar[64];
   armZZbar[267]= - 145./36.*armZZbar[67];
   armZZbar[269]=1./8.*armZZbar[48];
   armZZbar[274]=7./2.*armZZbar[44];
   armZZbar[275]=4*armZZbar[15];
   armZZbar[276]= - 55./18. + 4*armZZbar[13];
   armZZbar[276]=armZZbar[36]*armZZbar[276];
   armZZbar[178]=armZZbar[178] - 2./3.*armZZbar[13] - 11./2.*
   armZZbar[44] - 46./9. - 11./8.*armZZbar[48];
   armZZbar[178]=armZZbar[35]*armZZbar[178];
   armZZbar[278]= - armZZbar[1]*armZZbar[36];
   armZZbar[279]=1./9.*armZZbar[278];
   armZZbar[281]= - armZZbar[39]*armZZbar[36];
   armZZbar[282]=1./27.*armZZbar[281];
   armZZbar[283]=3*armZZbar[38];
   armZZbar[284]=1 + armZZbar[283];
   armZZbar[285]=2*armZZbar[10]*armZZbar[284];
   armZZbar[286]=1 + armZZbar[38];
   armZZbar[286]=armZZbar[18]*armZZbar[101]*armZZbar[286];
   armZZbar[287]=6*armZZbar[286];
   armZZbar[288]= - armZZbar[37]*armZZbar[101];
   armZZbar[289]=6*armZZbar[288];
   armZZbar[290]=9./4.*armZZbar[54] + 4*armZZbar[53];
   armZZbar[290]=MMH*armZZbar[290];
   armZZbar[208]=armZZbar[233] + armZZbar[290] + armZZbar[110] + 
   armZZbar[289] + armZZbar[171] + armZZbar[287] + armZZbar[285] + 
   armZZbar[282] + armZZbar[279] + armZZbar[178] + armZZbar[222] + 
   armZZbar[276] + armZZbar[234] + armZZbar[275] + armZZbar[274] + 
   armZZbar[269] + armZZbar[267] + armZZbar[266] + armZZbar[264] + 
   armZZbar[263] + 439./144.*armZZbar[66] + armZZbar[261] + 67./9.*
   armZZbar[65] + armZZbar[260] + 1./4.*armZZbar[208] + armZZbar[254];
   armZZbar[208]=MMt*armZZbar[208];
   armZZbar[291]= - 5*armZZbar[48];
   armZZbar[292]=5381./81. + armZZbar[291];
   armZZbar[293]=7*armZZbar[44];
   armZZbar[294]=27*armZZbar[43];
   armZZbar[292]=armZZbar[293] + 1./4.*armZZbar[292] + armZZbar[294];
   armZZbar[295]=27./4.*armZZbar[13];
   armZZbar[296]= - 2*armZZbar[11];
   armZZbar[297]= - 35./18.*armZZbar[35];
   armZZbar[298]= - 1./9.*armZZbar[1];
   armZZbar[299]= - 1./27.*armZZbar[39];
   armZZbar[300]=4*armZZbar[10];
   armZZbar[301]=6*armZZbar[117];
   armZZbar[302]=armZZbar[105] + 1 + armZZbar[1];
   armZZbar[302]=1./3.*armZZbar[6]*armZZbar[302];
   armZZbar[303]=armZZbar[37]*armZZbar[101];
   armZZbar[304]=3*armZZbar[303];
   armZZbar[305]=3*armZZbar[47];
   armZZbar[306]=3./2.*armZZbar[45];
   armZZbar[292]=armZZbar[304] + armZZbar[302] + armZZbar[301] + 
   armZZbar[300] + armZZbar[299] + armZZbar[298] + armZZbar[297] + 
   armZZbar[296] + armZZbar[295] + armZZbar[306] + 1./2.*armZZbar[292]
    + armZZbar[305];
   armZZbar[292]=armZZbar[37]*armZZbar[292];
   armZZbar[307]= - 23./6.*armZZbar[15];
   armZZbar[308]=3*armZZbar[72];
   armZZbar[309]=11./3.*armZZbar[64];
   armZZbar[310]=armZZbar[283] + armZZbar[307] + armZZbar[143] + 7./3.*
   armZZbar[67] + armZZbar[309] + 5./6.*armZZbar[59] - 19./12.*
   armZZbar[49] + 19./12.*armZZbar[61] - 353./216. + armZZbar[308];
   armZZbar[311]=73./27. + 11./2.*armZZbar[44];
   armZZbar[311]=1./2.*armZZbar[311] + 7./27.*armZZbar[11];
   armZZbar[312]=1./3.*armZZbar[35];
   armZZbar[311]=1./2.*armZZbar[311] + armZZbar[312];
   armZZbar[311]=armZZbar[35]*armZZbar[311];
   armZZbar[313]=59./18.*armZZbar[35] + armZZbar[161] + 7./9. + 
   armZZbar[158];
   armZZbar[313]=armZZbar[10]*armZZbar[313];
   armZZbar[314]=17./27.*armZZbar[11];
   armZZbar[315]= - 1./12.*armZZbar[36];
   armZZbar[316]= - MMH*armZZbar[54];
   armZZbar[310]=1./3.*armZZbar[316] + 1./2.*armZZbar[313] + 
   armZZbar[311] + armZZbar[315] + 1./2.*armZZbar[310] + armZZbar[314];
   armZZbar[310]=MMH*armZZbar[310];
   armZZbar[311]=11./3.*armZZbar[17];
   armZZbar[181]=armZZbar[181] + armZZbar[311];
   armZZbar[181]=armZZbar[35]*armZZbar[181];
   armZZbar[313]=armZZbar[229] - 11./3. + armZZbar[283];
   armZZbar[313]=1./2.*armZZbar[17]*armZZbar[313];
   armZZbar[316]= - 151./16.*armZZbar[42] - 4*armZZbar[71] - 395./8.*
   armZZbar[70];
   armZZbar[317]=13*armZZbar[36];
   armZZbar[318]=1229./6. + armZZbar[317];
   armZZbar[318]=armZZbar[19]*armZZbar[318];
   armZZbar[319]=161./27.*armZZbar[35] + 673./27. + armZZbar[229];
   armZZbar[319]=armZZbar[18]*armZZbar[319];
   armZZbar[320]= - 67./36.*armZZbar[68];
   armZZbar[321]=7./3.*armZZbar[22];
   armZZbar[322]=107./18.*armZZbar[20];
   armZZbar[323]= - armZZbar[6]*armZZbar[19];
   armZZbar[324]=1./3.*armZZbar[323];
   armZZbar[325]= - 3*armZZbar[37] - 4*armZZbar[18] - 9./2.*
   armZZbar[19] + armZZbar[17];
   armZZbar[325]=3*armZZbar[3]*armZZbar[37]*armZZbar[325];
   armZZbar[326]= - 9./4.*armZZbar[40];
   armZZbar[181]=armZZbar[208] + armZZbar[310] + armZZbar[325] + 
   armZZbar[292] + armZZbar[324] + 1./6.*armZZbar[319] + 1./3.*
   armZZbar[181] + armZZbar[313] + 1./12.*armZZbar[318] - 349./324.*
   armZZbar[21] + armZZbar[322] + armZZbar[321] + 2027./162.*
   armZZbar[50] + 673./324.*armZZbar[41] + armZZbar[326] + 
   armZZbar[320] + 1./9.*armZZbar[316] - armZZbar[69];
   armZZbar[181]=MMt*armZZbar[181];
   armZZbar[208]=1./4.*armZZbar[36];
   armZZbar[292]=armZZbar[144] + armZZbar[208];
   armZZbar[292]=armZZbar[18]*armZZbar[292];
   armZZbar[316]= - 1 - armZZbar[36];
   armZZbar[318]=armZZbar[19]*armZZbar[316];
   armZZbar[319]= - 61./9. + armZZbar[158];
   armZZbar[319]=armZZbar[17]*armZZbar[319];
   armZZbar[327]=armZZbar[35]*armZZbar[17];
   armZZbar[328]=3./2.*armZZbar[40];
   armZZbar[329]= - armZZbar[22] - 4*armZZbar[50] + armZZbar[42] + 
   armZZbar[328];
   armZZbar[329]=3*armZZbar[329];
   armZZbar[236]=armZZbar[236] - 3./2. - armZZbar[48];
   armZZbar[236]=3*armZZbar[37]*armZZbar[236];
   armZZbar[292]=armZZbar[236] + 3*armZZbar[292] + 7./18.*armZZbar[327]
    + 1./2.*armZZbar[319] + 3./2.*armZZbar[318] + armZZbar[329] - 37./9.
   *armZZbar[20];
   armZZbar[292]=armZZbar[3]*armZZbar[292];
   armZZbar[318]=armZZbar[240] + armZZbar[184];
   armZZbar[318]=11./18.*armZZbar[142] + 1./3.*armZZbar[318] + 1./2.*
   armZZbar[121];
   armZZbar[318]=armZZbar[6]*armZZbar[318];
   armZZbar[319]=7./3.*armZZbar[54] + armZZbar[52] + armZZbar[51];
   armZZbar[251]=1./2.*armZZbar[319] + armZZbar[251];
   armZZbar[319]=3*armZZbar[101];
   armZZbar[249]=1./2.*armZZbar[249] - 2*armZZbar[44];
   armZZbar[249]=3*armZZbar[3]*armZZbar[249];
   armZZbar[251]=armZZbar[249] + 1./3.*armZZbar[251] + armZZbar[319];
   armZZbar[251]=MMt*armZZbar[251];
   armZZbar[330]=1./2.*armZZbar[40];
   armZZbar[209]=1./2.*armZZbar[209] + armZZbar[198] + armZZbar[330] - 
   armZZbar[50];
   armZZbar[209]=3*armZZbar[101]*armZZbar[209];
   armZZbar[331]=armZZbar[284] + armZZbar[229];
   armZZbar[331]=armZZbar[10]*armZZbar[331];
   armZZbar[332]= - 1433./324. - 3*armZZbar[60];
   armZZbar[332]=5*armZZbar[332] - 77./36.*armZZbar[63];
   armZZbar[333]=1./2.*armZZbar[13];
   armZZbar[334]= - 7./9. + armZZbar[333];
   armZZbar[334]=1./4.*armZZbar[36]*armZZbar[334];
   armZZbar[335]=4./3.*armZZbar[35] + armZZbar[333] + 7./4.*
   armZZbar[44] - 13./3. + 7./16.*armZZbar[48];
   armZZbar[335]=armZZbar[35]*armZZbar[335];
   armZZbar[336]= - 13./24.*armZZbar[54] + armZZbar[170];
   armZZbar[336]=MMH*armZZbar[336];
   armZZbar[337]= - 2*armZZbar[72];
   armZZbar[338]=1./2.*armZZbar[65];
   armZZbar[339]=2*armZZbar[15];
   armZZbar[286]=3*armZZbar[286];
   armZZbar[340]=3*armZZbar[288];
   armZZbar[251]=armZZbar[251] + armZZbar[336] + armZZbar[292] + 
   armZZbar[340] + armZZbar[318] + armZZbar[286] + armZZbar[331] + 11./
   54.*armZZbar[281] + 1./6.*armZZbar[278] + 1./3.*armZZbar[335] + 
   armZZbar[209] + armZZbar[334] + armZZbar[158] + armZZbar[339] + 61./
   12.*armZZbar[44] + 43./48.*armZZbar[48] - 595./216.*armZZbar[67] + 7.
   /18.*armZZbar[64] + 15./4.*armZZbar[16] + 277./144.*armZZbar[49] - 
   967./288.*armZZbar[66] - 7./18.*armZZbar[61] + armZZbar[338] - 
   armZZbar[62] + 1./4.*armZZbar[332] + armZZbar[337];
   armZZbar[251]=MMt*armZZbar[251];
   armZZbar[292]= - 47*armZZbar[70] - 275./6.*armZZbar[42];
   armZZbar[292]= - 661./36.*armZZbar[68] + 1./8.*armZZbar[292] + 7*
   armZZbar[69];
   armZZbar[318]= - 983./36. + armZZbar[36];
   armZZbar[318]=armZZbar[19]*armZZbar[318];
   armZZbar[331]=armZZbar[36] + 515./27. + armZZbar[283];
   armZZbar[331]=armZZbar[17]*armZZbar[331];
   armZZbar[332]=11*armZZbar[17];
   armZZbar[335]= - 205./3.*armZZbar[19] + armZZbar[332];
   armZZbar[335]=armZZbar[35]*armZZbar[335];
   armZZbar[336]=71./54.*armZZbar[35] + 173./27. + armZZbar[36];
   armZZbar[336]=armZZbar[18]*armZZbar[336];
   armZZbar[292]=1./9.*armZZbar[323] + 1./6.*armZZbar[336] + 1./36.*
   armZZbar[335] + 1./2.*armZZbar[331] + 1./6.*armZZbar[318] - 869./324.
   *armZZbar[21] + 251./54.*armZZbar[20] + 61./36.*armZZbar[22] + 1957./
   648.*armZZbar[50] + 113./324.*armZZbar[41] + 1./3.*armZZbar[292] + 
   armZZbar[326];
   armZZbar[318]=armZZbar[283] - 29./18.*armZZbar[15] - 43./12.*
   armZZbar[44] + 31./9.*armZZbar[67] - 67./9.*armZZbar[64] - 25./18.*
   armZZbar[59] - 97./36.*armZZbar[49] + 97./36.*armZZbar[61] - 1897./
   216. + armZZbar[308];
   armZZbar[143]= - 1./9. + armZZbar[143];
   armZZbar[331]=1./9.*armZZbar[11];
   armZZbar[143]=1./2.*armZZbar[143] + armZZbar[331];
   armZZbar[335]= - 1./3.*armZZbar[35];
   armZZbar[143]=7./2.*armZZbar[143] + armZZbar[335];
   armZZbar[143]=armZZbar[35]*armZZbar[143];
   armZZbar[336]= - 229./27. + armZZbar[158];
   armZZbar[336]= - 131./108.*armZZbar[35] + 1./2.*armZZbar[336] + 
   armZZbar[161];
   armZZbar[336]=armZZbar[10]*armZZbar[336];
   armZZbar[341]=MMH*armZZbar[54];
   armZZbar[143]=1./9.*armZZbar[341] + 1./2.*armZZbar[336] + 1./3.*
   armZZbar[143] + armZZbar[315] + 1./4.*armZZbar[318] + armZZbar[314];
   armZZbar[143]=MMH*armZZbar[143];
   armZZbar[315]=2659./27. + armZZbar[238];
   armZZbar[315]=61./3.*armZZbar[44] + 1./3.*armZZbar[315] + 
   armZZbar[294];
   armZZbar[315]= - 29./4.*armZZbar[13] + armZZbar[306] + 1./2.*
   armZZbar[315] + armZZbar[305];
   armZZbar[318]=armZZbar[230] + 1./9. + armZZbar[1];
   armZZbar[318]=armZZbar[6]*armZZbar[318];
   armZZbar[117]=3*armZZbar[117];
   armZZbar[336]=3./2.*armZZbar[303];
   armZZbar[315]=armZZbar[336] + 1./2.*armZZbar[318] + armZZbar[117] + 
   7./2.*armZZbar[10] - 11./54.*armZZbar[39] - 1./6.*armZZbar[1] - 11./
   36.*armZZbar[35] + 1./2.*armZZbar[315] + armZZbar[296];
   armZZbar[315]=armZZbar[37]*armZZbar[315];
   armZZbar[318]=6*armZZbar[19] + armZZbar[119];
   armZZbar[318]= - 9./2.*armZZbar[37] + 3*armZZbar[318] + 11./4.*
   armZZbar[18];
   armZZbar[318]=armZZbar[3]*armZZbar[37]*armZZbar[318];
   armZZbar[143]=armZZbar[251] + armZZbar[143] + armZZbar[318] + 1./2.*
   armZZbar[292] + armZZbar[315];
   armZZbar[143]=MMt*armZZbar[143];
   armZZbar[251]=971./8. - 5*armZZbar[14];
   armZZbar[251]= - 25./8.*armZZbar[15] + 25./8.*armZZbar[67] + 523./12.
   *armZZbar[64] + 1./3.*armZZbar[251] + 25./8.*armZZbar[59];
   armZZbar[292]=11./9.*armZZbar[11];
   armZZbar[315]= - 49./48. + 17*armZZbar[35];
   armZZbar[315]=armZZbar[10]*armZZbar[315];
   armZZbar[318]=armZZbar[118] + armZZbar[10];
   armZZbar[341]=armZZbar[6]*armZZbar[318];
   armZZbar[342]=armZZbar[99]*armZZbar[50];
   armZZbar[343]= - armZZbar[99] + armZZbar[175];
   armZZbar[343]=armZZbar[37]*armZZbar[343];
   armZZbar[344]= - armZZbar[35]*armZZbar[11];
   armZZbar[251]=25./4.*armZZbar[343] + 5./3.*armZZbar[341] + 1./3.*
   armZZbar[315] + 17./6.*armZZbar[344] + 25./4.*armZZbar[342] + 1./2.*
   armZZbar[251] + armZZbar[292];
   armZZbar[251]=MMH*armZZbar[251];
   armZZbar[315]=25./12.*armZZbar[21];
   armZZbar[341]=armZZbar[315] - 227./18.*armZZbar[50] - 25./12.*
   armZZbar[41] + 455./9.*armZZbar[68] - 43./2.*armZZbar[40];
   armZZbar[345]=49./48. + armZZbar[199];
   armZZbar[345]=armZZbar[18]*armZZbar[345];
   armZZbar[346]=armZZbar[19] + armZZbar[17];
   armZZbar[347]=1./2.*armZZbar[346] - armZZbar[18];
   armZZbar[347]=armZZbar[6]*armZZbar[347];
   armZZbar[348]=armZZbar[35]*armZZbar[346];
   armZZbar[341]=5./9.*armZZbar[347] + 1./9.*armZZbar[345] + 17./18.*
   armZZbar[348] - 1511./108.*armZZbar[17] + 1./4.*armZZbar[341] - 11./
   27.*armZZbar[19];
   armZZbar[345]= - 1843./18. - 25*armZZbar[44];
   armZZbar[347]=17./9.*armZZbar[11];
   armZZbar[349]= - armZZbar[18]*armZZbar[99];
   armZZbar[350]=armZZbar[37]*armZZbar[99];
   armZZbar[345]=25./48.*armZZbar[350] + 25./48.*armZZbar[349] - 175./
   36.*armZZbar[10] + 1./16.*armZZbar[345] + armZZbar[347];
   armZZbar[345]=armZZbar[37]*armZZbar[345];
   armZZbar[251]=1./12.*armZZbar[251] + 1./4.*armZZbar[341] + 
   armZZbar[345];
   armZZbar[251]=MMH*armZZbar[251];
   armZZbar[341]=1205./24.*armZZbar[37] + 427./24.*armZZbar[18] - 149*
   armZZbar[19] - 91./8.*armZZbar[17];
   armZZbar[341]=armZZbar[37]*armZZbar[341];
   armZZbar[251]=1./9.*armZZbar[341] + armZZbar[251];
   armZZbar[341]=armZZbar[121] + 11./9.*armZZbar[142];
   armZZbar[341]=armZZbar[6]*armZZbar[341];
   armZZbar[345]= - 1./3. + armZZbar[333];
   armZZbar[345]=armZZbar[36]*armZZbar[345];
   armZZbar[351]=1./3.*armZZbar[345] + armZZbar[16] - armZZbar[66] - 1./
   4. - armZZbar[60];
   armZZbar[352]=armZZbar[18]*armZZbar[36];
   armZZbar[353]=armZZbar[3]*armZZbar[352];
   armZZbar[278]=1./3.*armZZbar[278];
   armZZbar[354]=armZZbar[10]*armZZbar[36];
   armZZbar[341]=3./2.*armZZbar[353] + armZZbar[341] + armZZbar[354] + 
   11./27.*armZZbar[281] + 1./2.*armZZbar[351] + armZZbar[278];
   armZZbar[341]=MMt*armZZbar[341];
   armZZbar[351]= - 11./18.*armZZbar[39] - 1./2.*armZZbar[1] + 2./3. - 
   5./8.*armZZbar[13];
   armZZbar[230]=armZZbar[1] + armZZbar[230];
   armZZbar[353]=armZZbar[6]*armZZbar[230];
   armZZbar[355]=1./2.*armZZbar[10];
   armZZbar[351]=1./2.*armZZbar[353] + 1./3.*armZZbar[351] + 
   armZZbar[355];
   armZZbar[351]=armZZbar[37]*armZZbar[351];
   armZZbar[353]=armZZbar[35]*armZZbar[19];
   armZZbar[205]= - 17 + armZZbar[205];
   armZZbar[205]=armZZbar[18]*armZZbar[205];
   armZZbar[205]=1./27.*armZZbar[205] + 7./54.*armZZbar[353] - 
   armZZbar[50] + 17./27.*armZZbar[19];
   armZZbar[356]= - armZZbar[19] + 5./2.*armZZbar[18];
   armZZbar[356]=armZZbar[3]*armZZbar[37]*armZZbar[356];
   armZZbar[205]=1./2.*armZZbar[341] + 1./2.*armZZbar[356] + 1./4.*
   armZZbar[205] + armZZbar[351];
   armZZbar[205]=MMt*armZZbar[205];
   armZZbar[341]=11./3.*armZZbar[11];
   armZZbar[351]= - armZZbar[6]*armZZbar[11];
   armZZbar[351]=5./2.*armZZbar[351] + armZZbar[341] + 17./2.*
   armZZbar[344];
   armZZbar[351]=MMH*armZZbar[351];
   armZZbar[356]=armZZbar[19] - armZZbar[17];
   armZZbar[357]=armZZbar[6]*armZZbar[356];
   armZZbar[358]=armZZbar[35]*armZZbar[356];
   armZZbar[359]= - armZZbar[19] + armZZbar[17];
   armZZbar[351]=armZZbar[351] + 5./2.*armZZbar[357] + 11./3.*
   armZZbar[359] + 17./2.*armZZbar[358];
   armZZbar[351]=MMH*armZZbar[351];
   armZZbar[357]=armZZbar[19] - armZZbar[18];
   armZZbar[360]=armZZbar[37]*armZZbar[357];
   armZZbar[351]=17*armZZbar[360] + armZZbar[351];
   armZZbar[205]=1./108.*armZZbar[351] + armZZbar[205];
   armZZbar[205]=armZZbar[272]*armZZbar[205];
   armZZbar[143]=armZZbar[205] + 1./3.*armZZbar[251] + armZZbar[143];
   armZZbar[143]=armZZbar[272]*armZZbar[143];
   armZZbar[205]=5*armZZbar[44];
   armZZbar[251]= - 299./54. + armZZbar[205];
   armZZbar[351]=armZZbar[18]*armZZbar[99];
   armZZbar[361]= - armZZbar[37]*armZZbar[99];
   armZZbar[362]=25./18.*armZZbar[10];
   armZZbar[251]=5./24.*armZZbar[361] + 5./24.*armZZbar[351] + 
   armZZbar[362] + 1./8.*armZZbar[251] + armZZbar[314];
   armZZbar[251]=armZZbar[37]*armZZbar[251];
   armZZbar[361]=17./8. - armZZbar[14];
   armZZbar[176]=armZZbar[280] - 5./8.*armZZbar[67] + 1./12.*
   armZZbar[64] + 1./3.*armZZbar[361] + armZZbar[176];
   armZZbar[280]= - armZZbar[99]*armZZbar[50];
   armZZbar[176]=5./4.*armZZbar[280] + 1./2.*armZZbar[176] + 1./27.*
   armZZbar[11];
   armZZbar[252]=armZZbar[252] + armZZbar[335];
   armZZbar[252]=armZZbar[10]*armZZbar[252];
   armZZbar[280]=armZZbar[216] + armZZbar[355];
   armZZbar[280]=armZZbar[6]*armZZbar[280];
   armZZbar[355]=armZZbar[99] + armZZbar[225];
   armZZbar[355]=armZZbar[37]*armZZbar[355];
   armZZbar[176]=5./8.*armZZbar[355] + 1./6.*armZZbar[280] + 1./4.*
   armZZbar[252] + 1./2.*armZZbar[176] + 1./9.*armZZbar[344];
   armZZbar[176]=MMH*armZZbar[176];
   armZZbar[252]= - 1./2.*armZZbar[40];
   armZZbar[280]= - 5./12.*armZZbar[21] + 151./18.*armZZbar[50] + 5./12.
   *armZZbar[41] + 5./9.*armZZbar[68] + armZZbar[252];
   armZZbar[355]= - 5./54. - armZZbar[36];
   armZZbar[355]=armZZbar[17]*armZZbar[355];
   armZZbar[280]=1./6.*armZZbar[355] + 1./4.*armZZbar[280] - 1./81.*
   armZZbar[19];
   armZZbar[355]= - 23./24. + armZZbar[36];
   armZZbar[361]=armZZbar[355] + armZZbar[312];
   armZZbar[361]=armZZbar[18]*armZZbar[361];
   armZZbar[363]= - armZZbar[19] + 5./2.*armZZbar[17];
   armZZbar[123]=1./3.*armZZbar[363] + armZZbar[123];
   armZZbar[123]=armZZbar[6]*armZZbar[123];
   armZZbar[363]=armZZbar[19] - 7./4.*armZZbar[17];
   armZZbar[363]=armZZbar[35]*armZZbar[363];
   armZZbar[123]=1./3.*armZZbar[176] + armZZbar[251] + 1./18.*
   armZZbar[123] + 1./12.*armZZbar[361] + 1./2.*armZZbar[280] + 1./27.*
   armZZbar[363];
   armZZbar[123]=MMH*armZZbar[123];
   armZZbar[176]=385./27.*armZZbar[37] + 1427./27.*armZZbar[18] + 103*
   armZZbar[19] - 53./3.*armZZbar[17];
   armZZbar[176]=armZZbar[37]*armZZbar[176];
   armZZbar[143]=armZZbar[143] + armZZbar[181] + 1./12.*armZZbar[176]
    + armZZbar[123];
   armZZbar[143]=armZZbar[272]*armZZbar[143];
   armZZbar[176]=armZZbar[207] - 289./12. + armZZbar[206];
   armZZbar[110]=armZZbar[233] + armZZbar[290] + armZZbar[110] + 
   armZZbar[289] + armZZbar[171] + armZZbar[287] + armZZbar[285] + 
   armZZbar[282] + armZZbar[279] + armZZbar[178] + armZZbar[222] + 
   armZZbar[276] + armZZbar[234] + armZZbar[275] + armZZbar[274] + 
   armZZbar[269] + armZZbar[267] + armZZbar[266] + armZZbar[264] + 
   armZZbar[263] + 301./48.*armZZbar[66] + armZZbar[261] + 32./3.*
   armZZbar[65] + armZZbar[260] + 1./4.*armZZbar[176] + armZZbar[254];
   armZZbar[110]=MMt*armZZbar[110];
   armZZbar[171]=589./9. + armZZbar[291];
   armZZbar[171]=armZZbar[293] + 1./4.*armZZbar[171] + armZZbar[294];
   armZZbar[171]=armZZbar[304] + armZZbar[302] + armZZbar[301] + 
   armZZbar[300] + armZZbar[299] + armZZbar[298] + armZZbar[297] + 
   armZZbar[296] + armZZbar[295] + armZZbar[306] + 1./2.*armZZbar[171]
    + armZZbar[305];
   armZZbar[171]=armZZbar[37]*armZZbar[171];
   armZZbar[176]= - 4./3.*armZZbar[71];
   armZZbar[178]= - 205./16.*armZZbar[42] + armZZbar[176] + 23./8.*
   armZZbar[70];
   armZZbar[181]=4855./18. + armZZbar[317];
   armZZbar[181]=armZZbar[19]*armZZbar[181];
   armZZbar[206]=43./3.*armZZbar[19] + armZZbar[17];
   armZZbar[206]=armZZbar[35]*armZZbar[206];
   armZZbar[207]= - 13*armZZbar[35] - 13 + armZZbar[229];
   armZZbar[207]=armZZbar[18]*armZZbar[207];
   armZZbar[110]=armZZbar[110] + armZZbar[310] + armZZbar[325] + 
   armZZbar[171] + armZZbar[324] + 1./6.*armZZbar[207] + 11./9.*
   armZZbar[206] + armZZbar[313] + 1./12.*armZZbar[181] + armZZbar[315]
    + armZZbar[322] + armZZbar[321] + 113./6.*armZZbar[50] - 13./12.*
   armZZbar[41] + armZZbar[326] + armZZbar[320] + 1./3.*armZZbar[178]
    - armZZbar[69];
   armZZbar[110]=MMt*armZZbar[110];
   armZZbar[171]= - 1./4.*armZZbar[36];
   armZZbar[178]=armZZbar[144] + armZZbar[171];
   armZZbar[178]=armZZbar[18]*armZZbar[178];
   armZZbar[181]= - 1./2. - armZZbar[36];
   armZZbar[181]=armZZbar[19]*armZZbar[181];
   armZZbar[206]= - 5 + armZZbar[158];
   armZZbar[206]=armZZbar[17]*armZZbar[206];
   armZZbar[173]=armZZbar[236] + 3*armZZbar[178] + 1./2.*armZZbar[173]
    + 1./2.*armZZbar[206] + 3*armZZbar[181] + armZZbar[329] - 5*
   armZZbar[20];
   armZZbar[173]=armZZbar[3]*armZZbar[173];
   armZZbar[178]= - armZZbar[54] - armZZbar[52] + armZZbar[51];
   armZZbar[178]=armZZbar[249] + 1./2.*armZZbar[178] + armZZbar[319];
   armZZbar[178]=MMt*armZZbar[178];
   armZZbar[181]= - 1./2.*armZZbar[36];
   armZZbar[206]=armZZbar[284] + armZZbar[181];
   armZZbar[206]=armZZbar[10]*armZZbar[206];
   armZZbar[207]=armZZbar[281] + armZZbar[13] + armZZbar[278];
   armZZbar[207]=armZZbar[6]*armZZbar[207];
   armZZbar[222]= - 5./8.*armZZbar[63] - 223./12. + 13*armZZbar[60];
   armZZbar[233]=1./9. + armZZbar[210];
   armZZbar[233]=armZZbar[36]*armZZbar[233];
   armZZbar[234]=armZZbar[333] - 3./4.*armZZbar[44] - 1 - 3./16.*
   armZZbar[48];
   armZZbar[234]=armZZbar[35]*armZZbar[234];
   armZZbar[170]=1./8.*armZZbar[54] + armZZbar[170];
   armZZbar[170]=MMH*armZZbar[170];
   armZZbar[142]=armZZbar[178] + armZZbar[170] + armZZbar[173] + 
   armZZbar[340] + 1./2.*armZZbar[207] + armZZbar[286] + armZZbar[206]
    + 1./6.*armZZbar[142] + 1./18.*armZZbar[121] + armZZbar[234] + 
   armZZbar[209] + 5./4.*armZZbar[233] + armZZbar[158] + armZZbar[339]
    + 15./4.*armZZbar[44] + 9./16.*armZZbar[48] - 59./24.*armZZbar[67]
    + armZZbar[244] - 13./2.*armZZbar[16] + 13./16.*armZZbar[49] + 365./
   32.*armZZbar[66] + armZZbar[218] + armZZbar[338] - armZZbar[62] + 1./
   2.*armZZbar[222] + armZZbar[337];
   armZZbar[142]=MMt*armZZbar[142];
   armZZbar[170]= - 161./9. + 3./2.*armZZbar[48];
   armZZbar[170]=15*armZZbar[44] + 1./2.*armZZbar[170] + armZZbar[294];
   armZZbar[173]=115./4.*armZZbar[13];
   armZZbar[104]=armZZbar[10] + armZZbar[105] + armZZbar[104] - 3./2.*
   armZZbar[35] + armZZbar[173] + armZZbar[306] + 1./2.*armZZbar[170]
    + armZZbar[305];
   armZZbar[105]= - armZZbar[39] + 1 + armZZbar[138];
   armZZbar[105]=armZZbar[6]*armZZbar[105];
   armZZbar[104]=armZZbar[336] + 1./2.*armZZbar[105] + 1./2.*
   armZZbar[104] + armZZbar[117];
   armZZbar[104]=armZZbar[37]*armZZbar[104];
   armZZbar[105]= - 5./2.*armZZbar[15];
   armZZbar[117]=armZZbar[283] + armZZbar[105] - 9./4.*armZZbar[44] + 3
   *armZZbar[67] + armZZbar[155] + armZZbar[152] - 9./4.*armZZbar[49]
    + 9./4.*armZZbar[61] - 563./72. + armZZbar[308];
   armZZbar[138]= - 4./9. + 3./16.*armZZbar[44];
   armZZbar[138]=armZZbar[35]*armZZbar[138];
   armZZbar[152]= - 73./18.*armZZbar[35] - 59./9. + armZZbar[158];
   armZZbar[152]=armZZbar[10]*armZZbar[152];
   armZZbar[117]=1./4.*armZZbar[152] + 1./4.*armZZbar[117] + 
   armZZbar[138];
   armZZbar[117]=MMH*armZZbar[117];
   armZZbar[138]= - 129*armZZbar[70] + 157./2.*armZZbar[42];
   armZZbar[152]=1./4.*armZZbar[41];
   armZZbar[170]= - 5./4.*armZZbar[21];
   armZZbar[138]=armZZbar[170] + 31./6.*armZZbar[20] + 19./4.*
   armZZbar[22] + 841./24.*armZZbar[50] + armZZbar[152] + armZZbar[326]
    - 53./12.*armZZbar[68] + 1./8.*armZZbar[138] + armZZbar[69];
   armZZbar[178]= - 275./48. + 4*armZZbar[36];
   armZZbar[178]=armZZbar[19]*armZZbar[178];
   armZZbar[206]=47./3. + armZZbar[283];
   armZZbar[206]=armZZbar[17]*armZZbar[206];
   armZZbar[207]= - 115*armZZbar[19] + 47*armZZbar[17];
   armZZbar[207]=armZZbar[35]*armZZbar[207];
   armZZbar[209]=109 - 43./2.*armZZbar[35];
   armZZbar[209]=armZZbar[18]*armZZbar[209];
   armZZbar[115]=armZZbar[115] + armZZbar[17];
   armZZbar[115]= - 9*armZZbar[37] + 3*armZZbar[115] - 59./2.*
   armZZbar[18];
   armZZbar[115]=armZZbar[3]*armZZbar[37]*armZZbar[115];
   armZZbar[104]=armZZbar[142] + armZZbar[117] + 1./2.*armZZbar[115] + 
   armZZbar[104] + 1./2.*armZZbar[323] + 1./36.*armZZbar[209] + 1./24.*
   armZZbar[207] + 1./4.*armZZbar[206] + 1./2.*armZZbar[138] + 
   armZZbar[178];
   armZZbar[104]=MMt*armZZbar[104];
   armZZbar[115]=67./8. - armZZbar[14];
   armZZbar[115]= - 1./8.*armZZbar[15] + 1./8.*armZZbar[67] + 35./12.*
   armZZbar[64] + 1./3.*armZZbar[115] + 1./8.*armZZbar[59];
   armZZbar[117]=armZZbar[35]*armZZbar[11];
   armZZbar[138]= - 1./9.*armZZbar[11];
   armZZbar[115]=1./6.*armZZbar[117] + 1./4.*armZZbar[342] + 1./2.*
   armZZbar[115] + armZZbar[138];
   armZZbar[142]=1./4.*armZZbar[202] + armZZbar[312];
   armZZbar[142]=armZZbar[10]*armZZbar[142];
   armZZbar[178]=1./4.*armZZbar[11] + 1./3.*armZZbar[10];
   armZZbar[178]=armZZbar[6]*armZZbar[178];
   armZZbar[115]=1./16.*armZZbar[343] + 1./6.*armZZbar[178] + 1./4.*
   armZZbar[115] + 1./3.*armZZbar[142];
   armZZbar[115]=MMH*armZZbar[115];
   armZZbar[142]=1./4.*armZZbar[21];
   armZZbar[178]=armZZbar[142] + 5./6.*armZZbar[50] + armZZbar[190] + 
   31./3.*armZZbar[68] - 9./2.*armZZbar[40];
   armZZbar[202]= - 107./12. - armZZbar[36];
   armZZbar[202]=armZZbar[17]*armZZbar[202];
   armZZbar[206]= - armZZbar[19] + armZZbar[311];
   armZZbar[207]=armZZbar[35]*armZZbar[206];
   armZZbar[178]=1./6.*armZZbar[207] + 1./3.*armZZbar[202] + 1./4.*
   armZZbar[178] + 1./9.*armZZbar[19];
   armZZbar[202]= - 545./18. - 3*armZZbar[44];
   armZZbar[202]=1./4.*armZZbar[350] + 1./4.*armZZbar[349] + 1./4.*
   armZZbar[202] - 41./9.*armZZbar[10];
   armZZbar[202]=armZZbar[37]*armZZbar[202];
   armZZbar[207]= - 23./48. + armZZbar[36];
   armZZbar[209]=1./4.*armZZbar[207] + armZZbar[335];
   armZZbar[209]=armZZbar[18]*armZZbar[209];
   armZZbar[218]= - armZZbar[19] + 7./3.*armZZbar[17];
   armZZbar[222]= - 1./3.*armZZbar[18];
   armZZbar[218]=1./4.*armZZbar[218] + armZZbar[222];
   armZZbar[218]=armZZbar[6]*armZZbar[218];
   armZZbar[115]=armZZbar[115] + 1./4.*armZZbar[202] + 1./6.*
   armZZbar[218] + 1./4.*armZZbar[178] + 1./3.*armZZbar[209];
   armZZbar[115]=MMH*armZZbar[115];
   armZZbar[178]=armZZbar[278] + armZZbar[281];
   armZZbar[178]=armZZbar[6]*armZZbar[178];
   armZZbar[202]= - armZZbar[19]*armZZbar[36];
   armZZbar[209]= - armZZbar[18]*armZZbar[36];
   armZZbar[218]=armZZbar[202] + 1./2.*armZZbar[209];
   armZZbar[218]=armZZbar[3]*armZZbar[218];
   armZZbar[233]=59./9. + 33./2.*armZZbar[13];
   armZZbar[234]=armZZbar[36]*armZZbar[233];
   armZZbar[236]= - armZZbar[10]*armZZbar[36];
   armZZbar[121]=3*armZZbar[218] + armZZbar[178] + armZZbar[236] + 
   armZZbar[148] + 1./2.*armZZbar[234] + 1./9.*armZZbar[121];
   armZZbar[121]=MMt*armZZbar[121];
   armZZbar[148]=5*armZZbar[36];
   armZZbar[178]=17./3. + armZZbar[148];
   armZZbar[178]=armZZbar[19]*armZZbar[178];
   armZZbar[218]= - 47./6.*armZZbar[19] + armZZbar[17];
   armZZbar[218]=armZZbar[35]*armZZbar[218];
   armZZbar[249]=19./2.*armZZbar[35] - 17 - armZZbar[36];
   armZZbar[249]=armZZbar[18]*armZZbar[249];
   armZZbar[251]= - armZZbar[17]*armZZbar[36];
   armZZbar[178]=1./3.*armZZbar[249] + armZZbar[218] + armZZbar[178] + 
   armZZbar[251];
   armZZbar[218]= - 1./2. + armZZbar[11];
   armZZbar[249]=armZZbar[35]*armZZbar[218];
   armZZbar[254]= - armZZbar[35] + 1 + armZZbar[229];
   armZZbar[260]=armZZbar[10]*armZZbar[254];
   armZZbar[249]=armZZbar[260] + 1./2.*armZZbar[249] - armZZbar[11] + 
   armZZbar[208];
   armZZbar[249]=MMH*armZZbar[249];
   armZZbar[233]=1./2.*armZZbar[140] - 5./2.*armZZbar[10] + 1./6.*
   armZZbar[39] + 1./18.*armZZbar[1] + 1./4.*armZZbar[233] + 2*
   armZZbar[11];
   armZZbar[233]=armZZbar[37]*armZZbar[233];
   armZZbar[260]= - 10*armZZbar[19] + 31./4.*armZZbar[18];
   armZZbar[260]=armZZbar[3]*armZZbar[37]*armZZbar[260];
   armZZbar[121]=1./2.*armZZbar[121] + 1./3.*armZZbar[249] + 
   armZZbar[260] + 1./4.*armZZbar[178] + armZZbar[233];
   armZZbar[121]=MMt*armZZbar[121];
   armZZbar[178]= - armZZbar[19] - armZZbar[17];
   armZZbar[233]=armZZbar[35]*armZZbar[178];
   armZZbar[249]=armZZbar[35] - 1./3. - armZZbar[36];
   armZZbar[260]=armZZbar[18]*armZZbar[249];
   armZZbar[261]=1./2.*armZZbar[178] + armZZbar[18];
   armZZbar[261]=armZZbar[6]*armZZbar[261];
   armZZbar[263]=1./3.*armZZbar[19];
   armZZbar[264]=armZZbar[17]*armZZbar[36];
   armZZbar[233]=armZZbar[261] + armZZbar[260] + 1./2.*armZZbar[233] + 
   armZZbar[263] + armZZbar[264];
   armZZbar[260]=1./2.*armZZbar[11];
   armZZbar[261]=armZZbar[260] - armZZbar[10];
   armZZbar[261]=armZZbar[6]*armZZbar[261];
   armZZbar[117]=armZZbar[261] + armZZbar[150] + armZZbar[151] + 1./2.*
   armZZbar[117];
   armZZbar[117]=MMH*armZZbar[117];
   armZZbar[150]= - armZZbar[11] + armZZbar[10];
   armZZbar[261]=armZZbar[37]*armZZbar[150];
   armZZbar[117]=1./4.*armZZbar[117] + 1./4.*armZZbar[233] + 
   armZZbar[261];
   armZZbar[117]=MMH*armZZbar[117];
   armZZbar[117]=17./4.*armZZbar[360] + armZZbar[117];
   armZZbar[117]=1./3.*armZZbar[117] + armZZbar[121];
   armZZbar[117]=armZZbar[133]*armZZbar[117];
   armZZbar[121]=7./2.*armZZbar[17];
   armZZbar[233]=23./6.*armZZbar[37] + 227./18.*armZZbar[18] - 89./3.*
   armZZbar[19] + armZZbar[121];
   armZZbar[233]=armZZbar[37]*armZZbar[233];
   armZZbar[104]=armZZbar[117] + armZZbar[104] + 1./4.*armZZbar[233] + 
   armZZbar[115];
   armZZbar[104]=armZZbar[133]*armZZbar[104];
   armZZbar[115]=2207./3.*armZZbar[19] - 53*armZZbar[17];
   armZZbar[115]= - 71./3.*armZZbar[37] + 1./3.*armZZbar[115] - 23*
   armZZbar[18];
   armZZbar[115]=armZZbar[37]*armZZbar[115];
   armZZbar[104]=armZZbar[104] + armZZbar[110] + 1./12.*armZZbar[115]
    + armZZbar[123];
   armZZbar[104]=armZZbar[133]*armZZbar[104];
   armZZbar[110]=9*armZZbar[38];
   armZZbar[115]= - 16*armZZbar[35];
   armZZbar[117]=armZZbar[115] - 16 + armZZbar[110];
   armZZbar[117]=armZZbar[3]*armZZbar[37]*armZZbar[117];
   armZZbar[123]=8 + armZZbar[229];
   armZZbar[171]= - 4./3.*armZZbar[35] + 4./3. + armZZbar[171];
   armZZbar[171]=armZZbar[35]*armZZbar[171];
   armZZbar[166]= - MMt*armZZbar[166]*armZZbar[37]*armZZbar[38];
   armZZbar[233]=36*armZZbar[166];
   armZZbar[261]= - armZZbar[6]*armZZbar[36];
   armZZbar[266]=1./4.*armZZbar[261];
   armZZbar[117]=armZZbar[233] + armZZbar[117] + armZZbar[266] + 1./3.*
   armZZbar[123] + armZZbar[171];
   armZZbar[117]=MMt*armZZbar[117];
   armZZbar[123]=armZZbar[156] + 11 + 13./6.*armZZbar[35];
   armZZbar[123]=armZZbar[37]*armZZbar[123];
   armZZbar[171]=pow(armZZbar[37],2);
   armZZbar[267]= - armZZbar[3]*armZZbar[171];
   armZZbar[269]=16*armZZbar[267];
   armZZbar[117]=armZZbar[117] + 1./2.*armZZbar[123] + armZZbar[269];
   armZZbar[117]=MMt*armZZbar[117];
   armZZbar[123]=1 + armZZbar[134];
   armZZbar[123]=armZZbar[35]*armZZbar[123];
   armZZbar[134]=armZZbar[36] - armZZbar[35];
   armZZbar[134]=armZZbar[3]*armZZbar[37]*armZZbar[134];
   armZZbar[123]=3*armZZbar[134] + armZZbar[266] + armZZbar[161] + 1./2.
   *armZZbar[123];
   armZZbar[123]=MMt*armZZbar[123];
   armZZbar[134]=1./2.*armZZbar[35];
   armZZbar[266]=armZZbar[156] + armZZbar[134] + 1./3. - armZZbar[36];
   armZZbar[266]=armZZbar[37]*armZZbar[266];
   armZZbar[123]=1./2.*armZZbar[266] + armZZbar[123];
   armZZbar[123]=armZZbar[133]*MMt*armZZbar[123];
   armZZbar[117]=armZZbar[123] + 8./3.*armZZbar[171] + armZZbar[117];
   armZZbar[117]=armZZbar[133]*armZZbar[117];
   armZZbar[123]=5./2.*armZZbar[36];
   armZZbar[266]= - 7./6.*armZZbar[35] - 113./3. + armZZbar[123];
   armZZbar[266]=armZZbar[35]*armZZbar[266];
   armZZbar[261]=1./3.*armZZbar[261];
   armZZbar[266]=armZZbar[261] + 1./2.*armZZbar[266] - 16 + 49./18.*
   armZZbar[36];
   armZZbar[274]=32./3. + armZZbar[110];
   armZZbar[275]= - 3*armZZbar[36];
   armZZbar[274]=73./3.*armZZbar[35] + 2*armZZbar[274] + armZZbar[275];
   armZZbar[274]=armZZbar[3]*armZZbar[37]*armZZbar[274];
   armZZbar[166]=72*armZZbar[166];
   armZZbar[266]=armZZbar[166] + 1./3.*armZZbar[266] + armZZbar[274];
   armZZbar[266]=MMt*armZZbar[266];
   armZZbar[276]= - 17./3. + armZZbar[229];
   armZZbar[276]= - armZZbar[6] + 17*armZZbar[276] - 109./2.*
   armZZbar[35];
   armZZbar[276]=armZZbar[37]*armZZbar[276];
   armZZbar[278]=armZZbar[3]*armZZbar[171];
   armZZbar[279]=64*armZZbar[278];
   armZZbar[276]=1./3.*armZZbar[276] + armZZbar[279];
   armZZbar[266]=1./3.*armZZbar[276] + armZZbar[266];
   armZZbar[266]=MMt*armZZbar[266];
   armZZbar[117]=armZZbar[117] - 16./3.*armZZbar[171] + armZZbar[266];
   armZZbar[117]=armZZbar[133]*armZZbar[117];
   armZZbar[148]=34./3. + armZZbar[148];
   armZZbar[266]=49./108.*armZZbar[35] + 175./54. + 2*armZZbar[36];
   armZZbar[266]=armZZbar[35]*armZZbar[266];
   armZZbar[276]=armZZbar[6]*armZZbar[36];
   armZZbar[148]=5./12.*armZZbar[276] + 4./9.*armZZbar[148] + 
   armZZbar[266];
   armZZbar[110]= - 7./3.*armZZbar[35] + armZZbar[275] - 16./3. + 
   armZZbar[110];
   armZZbar[110]=armZZbar[3]*armZZbar[37]*armZZbar[110];
   armZZbar[110]=armZZbar[233] + 1./3.*armZZbar[148] + armZZbar[110];
   armZZbar[110]=MMt*armZZbar[110];
   armZZbar[148]=5./2.*armZZbar[6];
   armZZbar[233]=17*armZZbar[36];
   armZZbar[266]=armZZbar[148] + 503./18.*armZZbar[35] + 511./9. + 
   armZZbar[233];
   armZZbar[266]=armZZbar[37]*armZZbar[266];
   armZZbar[266]=1./6.*armZZbar[266] + armZZbar[269];
   armZZbar[110]=1./3.*armZZbar[266] + armZZbar[110];
   armZZbar[110]=MMt*armZZbar[110];
   armZZbar[148]=armZZbar[200] + armZZbar[148];
   armZZbar[148]=armZZbar[37]*armZZbar[148];
   armZZbar[200]=armZZbar[35]*armZZbar[36];
   armZZbar[266]=5./2.*armZZbar[276] - 11./3.*armZZbar[36] + 17./2.*
   armZZbar[200];
   armZZbar[266]=MMt*armZZbar[266];
   armZZbar[266]=armZZbar[148] + armZZbar[266];
   armZZbar[266]=armZZbar[272]*MMt*armZZbar[266];
   armZZbar[110]=1./18.*armZZbar[266] + 136./81.*armZZbar[171] + 
   armZZbar[110];
   armZZbar[110]=armZZbar[272]*armZZbar[110];
   armZZbar[266]= - 944./3. + 49./2.*armZZbar[36];
   armZZbar[269]= - 613./27. + armZZbar[229];
   armZZbar[269]=5*armZZbar[269] - 2111./54.*armZZbar[35];
   armZZbar[269]=armZZbar[35]*armZZbar[269];
   armZZbar[261]=armZZbar[261] + 1./9.*armZZbar[266] + 1./2.*
   armZZbar[269];
   armZZbar[166]=armZZbar[166] + 1./3.*armZZbar[261] + armZZbar[274];
   armZZbar[166]=MMt*armZZbar[166];
   armZZbar[261]=17./2.*armZZbar[36];
   armZZbar[266]= - armZZbar[6] - 3029./18.*armZZbar[35] - 1891./9. + 
   armZZbar[261];
   armZZbar[266]=armZZbar[37]*armZZbar[266];
   armZZbar[266]=1./3.*armZZbar[266] + armZZbar[279];
   armZZbar[166]=1./3.*armZZbar[266] + armZZbar[166];
   armZZbar[166]=MMt*armZZbar[166];
   armZZbar[110]=armZZbar[110] - 944./81.*armZZbar[171] + armZZbar[166]
   ;
   armZZbar[110]=armZZbar[272]*armZZbar[110];
   armZZbar[166]=armZZbar[37]*armZZbar[248];
   armZZbar[266]=2 + armZZbar[35];
   armZZbar[269]=armZZbar[35]*armZZbar[266];
   armZZbar[269]=1 + armZZbar[269];
   armZZbar[269]=MMt*armZZbar[269];
   armZZbar[166]=2*armZZbar[166] + armZZbar[269];
   armZZbar[166]=MMt*armZZbar[166];
   armZZbar[166]=armZZbar[171] + armZZbar[166];
   armZZbar[110]=armZZbar[117] + 1024./81.*armZZbar[166] + 
   armZZbar[110];
   armZZbar[110]=armZZbar[33]*armZZbar[110];
   armZZbar[117]=armZZbar[21] - armZZbar[41] + 2*armZZbar[50];
   armZZbar[166]= - armZZbar[35]*armZZbar[19];
   armZZbar[269]= - 2 - armZZbar[35];
   armZZbar[269]=armZZbar[18]*armZZbar[269];
   armZZbar[117]= - 4./9.*MMt + 4./9.*armZZbar[37] + 2./9.*
   armZZbar[269] + armZZbar[166] + 2./9.*armZZbar[117] - armZZbar[19];
   armZZbar[117]=MMt*armZZbar[117];
   armZZbar[269]= - 2./9.*armZZbar[37] - armZZbar[19] - 4./9.*
   armZZbar[18];
   armZZbar[269]=armZZbar[37]*armZZbar[269];
   armZZbar[117]=armZZbar[269] + armZZbar[117];
   armZZbar[104]=armZZbar[110] + armZZbar[104] + 256./9.*armZZbar[117]
    + armZZbar[143];
   armZZbar[104]=armZZbar[33]*armZZbar[104];
   armZZbar[110]= - 1./2. + armZZbar[158];
   armZZbar[117]=armZZbar[110] + armZZbar[161];
   armZZbar[117]=armZZbar[10]*armZZbar[117];
   armZZbar[143]=7./9. + armZZbar[72];
   armZZbar[269]= - 1./6.*armZZbar[36];
   armZZbar[274]= - MMH*armZZbar[53];
   armZZbar[117]=armZZbar[274] + armZZbar[117] + armZZbar[269] + 
   armZZbar[105] + 91./36.*armZZbar[67] - 7./9.*armZZbar[64] - 
   armZZbar[61] + 5./2.*armZZbar[143] + armZZbar[62];
   armZZbar[117]=MMH*armZZbar[117];
   armZZbar[143]=armZZbar[229] + 47./36. + armZZbar[158];
   armZZbar[143]=armZZbar[17]*armZZbar[143];
   armZZbar[276]=11./4. + armZZbar[283];
   armZZbar[279]=1./6.*armZZbar[36];
   armZZbar[280]=armZZbar[276] + armZZbar[279];
   armZZbar[280]=armZZbar[18]*armZZbar[280];
   armZZbar[281]=469./2.*armZZbar[70] - 115*armZZbar[42];
   armZZbar[281]=59./9.*armZZbar[40] - 65./18.*armZZbar[68] + 1./9.*
   armZZbar[281] - 1./2.*armZZbar[69];
   armZZbar[282]= - 775./6. + armZZbar[317];
   armZZbar[282]=armZZbar[19]*armZZbar[282];
   armZZbar[284]=109 + 239./3.*armZZbar[66];
   armZZbar[284]=1./4.*armZZbar[284] - 23./3.*armZZbar[67];
   armZZbar[284]=MMt*armZZbar[284];
   armZZbar[285]= - 5*armZZbar[50];
   armZZbar[286]= - 9./4.*armZZbar[20];
   armZZbar[117]=1./3.*armZZbar[284] + armZZbar[117] - 79./6.*
   armZZbar[37] + armZZbar[280] + armZZbar[143] + 1./6.*armZZbar[282]
    + armZZbar[142] + armZZbar[286] + 33./4.*armZZbar[22] + 1./2.*
   armZZbar[281] + armZZbar[285];
   armZZbar[117]=MMt*armZZbar[117];
   armZZbar[143]=armZZbar[188] + 229./54. + armZZbar[283];
   armZZbar[143]=1./2.*armZZbar[143] + 7./27.*armZZbar[35];
   armZZbar[143]=armZZbar[10]*armZZbar[143];
   armZZbar[143]=armZZbar[143] + 7./54.*armZZbar[344] - 17./27.*
   armZZbar[11] + 29./36.*armZZbar[15] + 25./36.*armZZbar[67] + 
   armZZbar[309] + 7./36.*armZZbar[59] + 107./36. - armZZbar[72];
   armZZbar[143]=MMH*armZZbar[143];
   armZZbar[161]=armZZbar[161] - 229./54. + armZZbar[158];
   armZZbar[161]=1./2.*armZZbar[161] - 7./27.*armZZbar[35];
   armZZbar[161]=armZZbar[18]*armZZbar[161];
   armZZbar[280]= - 7./36.*armZZbar[41] - 43./18.*armZZbar[40] + 
   armZZbar[69] + 41./9.*armZZbar[68];
   armZZbar[281]=armZZbar[279] - 139./27. + 3./2.*armZZbar[38];
   armZZbar[281]=armZZbar[17]*armZZbar[281];
   armZZbar[282]= - 11./6.*armZZbar[10] - 17./3. + armZZbar[11];
   armZZbar[282]=armZZbar[37]*armZZbar[282];
   armZZbar[143]=1./4.*armZZbar[143] + 1./2.*armZZbar[282] + 1./4.*
   armZZbar[161] + 7./216.*armZZbar[348] + 1./4.*armZZbar[281] + 17./
   108.*armZZbar[19] - 29./144.*armZZbar[21] + 3./8.*armZZbar[20] + 1./
   4.*armZZbar[280] - 2./9.*armZZbar[50];
   armZZbar[143]=MMH*armZZbar[143];
   armZZbar[161]=47./3.*armZZbar[18] - 85./3.*armZZbar[19] - 
   armZZbar[17];
   armZZbar[280]=5*armZZbar[37];
   armZZbar[161]=1./2.*armZZbar[161] + armZZbar[280];
   armZZbar[161]=armZZbar[37]*armZZbar[161];
   armZZbar[117]=1./2.*armZZbar[117] + 1./4.*armZZbar[161] + 
   armZZbar[143];
   armZZbar[117]=MMt*armZZbar[117];
   armZZbar[143]=armZZbar[181] + armZZbar[236];
   armZZbar[143]=MMH*armZZbar[143];
   armZZbar[161]= - 3 + armZZbar[188];
   armZZbar[161]=armZZbar[19]*armZZbar[161];
   armZZbar[281]=1 + armZZbar[66];
   armZZbar[281]=MMt*armZZbar[281];
   armZZbar[143]=1./2.*armZZbar[281] + 1./3.*armZZbar[143] - 
   armZZbar[37] + 1./6.*armZZbar[352] + 1./2.*armZZbar[264] + 1./2.*
   armZZbar[161] + armZZbar[192] + 1./2.*armZZbar[22];
   armZZbar[143]=MMt*armZZbar[143];
   armZZbar[161]=17*armZZbar[356] + 7./2.*armZZbar[358];
   armZZbar[281]= - 1./3.*armZZbar[10];
   armZZbar[282]=armZZbar[281] - 1./6. + armZZbar[11];
   armZZbar[282]=armZZbar[37]*armZZbar[282];
   armZZbar[284]= - 17*armZZbar[11] + 7./2.*armZZbar[344];
   armZZbar[284]=MMH*armZZbar[284];
   armZZbar[161]=1./54.*armZZbar[284] + 1./54.*armZZbar[161] + 
   armZZbar[282];
   armZZbar[161]=MMH*armZZbar[161];
   armZZbar[282]=1./12.*armZZbar[18] - 2./3.*armZZbar[19] + 3./4.*
   armZZbar[17];
   armZZbar[282]=armZZbar[37]*armZZbar[282];
   armZZbar[143]=1./2.*armZZbar[143] + armZZbar[282] + 1./2.*
   armZZbar[161];
   armZZbar[143]=MMt*armZZbar[143];
   armZZbar[161]=armZZbar[37]*armZZbar[356];
   armZZbar[282]= - MMH*armZZbar[37]*armZZbar[11];
   armZZbar[161]=armZZbar[161] + armZZbar[282];
   armZZbar[161]=MMH*armZZbar[161];
   armZZbar[143]=17./108.*armZZbar[161] + armZZbar[143];
   armZZbar[143]=armZZbar[272]*armZZbar[143];
   armZZbar[161]=2*armZZbar[17];
   armZZbar[282]= - 175./16.*armZZbar[18] + 17./4.*armZZbar[19] + 
   armZZbar[161];
   armZZbar[282]=1./3.*armZZbar[282] + 25./16.*armZZbar[37];
   armZZbar[282]=armZZbar[37]*armZZbar[282];
   armZZbar[284]=175./12.*armZZbar[10];
   armZZbar[287]=armZZbar[284] - 25./4. - 17./3.*armZZbar[11];
   armZZbar[287]=armZZbar[37]*armZZbar[287];
   armZZbar[287]=25./4.*armZZbar[50] + armZZbar[287];
   armZZbar[287]=MMH*armZZbar[287];
   armZZbar[282]=armZZbar[282] + 1./4.*armZZbar[287];
   armZZbar[282]=MMH*armZZbar[282];
   armZZbar[117]=armZZbar[143] + 1./9.*armZZbar[282] + armZZbar[117];
   armZZbar[117]=armZZbar[272]*armZZbar[117];
   armZZbar[143]=1./3. + armZZbar[72];
   armZZbar[282]=armZZbar[10]*armZZbar[110];
   armZZbar[143]=armZZbar[274] + armZZbar[282] + armZZbar[105] - 23./12.
   *armZZbar[67] + armZZbar[309] - armZZbar[61] + 5./2.*armZZbar[143]
    + armZZbar[62];
   armZZbar[143]=MMH*armZZbar[143];
   armZZbar[282]=49*armZZbar[70] - 23*armZZbar[42];
   armZZbar[282]=85./3.*armZZbar[68] + 1./3.*armZZbar[282] - 
   armZZbar[69];
   armZZbar[282]=1./2.*armZZbar[282] - 7./3.*armZZbar[40];
   armZZbar[287]= - 40 + 19*armZZbar[36];
   armZZbar[287]=armZZbar[19]*armZZbar[287];
   armZZbar[289]= - 91./12. + armZZbar[158];
   armZZbar[289]=armZZbar[17]*armZZbar[289];
   armZZbar[290]=armZZbar[18]*armZZbar[276];
   armZZbar[291]=27 + 13./3.*armZZbar[66];
   armZZbar[291]=1./2.*armZZbar[291] + 19./3.*armZZbar[67];
   armZZbar[291]=MMt*armZZbar[291];
   armZZbar[143]=armZZbar[291] + armZZbar[143] - 22*armZZbar[37] + 
   armZZbar[290] + armZZbar[289] + 1./3.*armZZbar[287] + armZZbar[142]
    + armZZbar[286] + 9*armZZbar[22] + 1./2.*armZZbar[282] + 
   armZZbar[285];
   armZZbar[143]=MMt*armZZbar[143];
   armZZbar[282]= - 11./9.*armZZbar[35];
   armZZbar[287]=armZZbar[282] + armZZbar[279] - 7./18. + armZZbar[283]
   ;
   armZZbar[287]=armZZbar[10]*armZZbar[287];
   armZZbar[289]=23./12.*armZZbar[15] - 5./12.*armZZbar[67] + 
   armZZbar[155] - 11./12.*armZZbar[59] - 31./12. - armZZbar[72];
   armZZbar[287]=1./4.*armZZbar[287] + 2./27.*armZZbar[344] + 1./2.*
   armZZbar[289] - 2./27.*armZZbar[11];
   armZZbar[287]=MMH*armZZbar[287];
   armZZbar[289]=armZZbar[279] + 283./27. + armZZbar[283];
   armZZbar[289]=armZZbar[17]*armZZbar[289];
   armZZbar[290]=11./9.*armZZbar[35] + armZZbar[269] + 7./18. + 
   armZZbar[158];
   armZZbar[290]=armZZbar[18]*armZZbar[290];
   armZZbar[291]=11./12.*armZZbar[41] - 1./6.*armZZbar[40] + 
   armZZbar[69] - 13./3.*armZZbar[68];
   armZZbar[293]=armZZbar[132] - 41./4.*armZZbar[17];
   armZZbar[293]=armZZbar[35]*armZZbar[293];
   armZZbar[294]=9./2. - armZZbar[10];
   armZZbar[294]=armZZbar[37]*armZZbar[294];
   armZZbar[287]=armZZbar[287] + armZZbar[294] + 1./4.*armZZbar[290] + 
   1./27.*armZZbar[293] + 1./4.*armZZbar[289] + 2./27.*armZZbar[19] - 
   23./24.*armZZbar[21] + 3./4.*armZZbar[20] + 1./2.*armZZbar[291] + 2./
   3.*armZZbar[50];
   armZZbar[287]=MMH*armZZbar[287];
   armZZbar[186]=5./2.*armZZbar[37] + armZZbar[186] - 8./3.*
   armZZbar[19] - 5./4.*armZZbar[17];
   armZZbar[186]=armZZbar[37]*armZZbar[186];
   armZZbar[143]=armZZbar[143] + armZZbar[186] + armZZbar[287];
   armZZbar[143]=MMt*armZZbar[143];
   armZZbar[186]=armZZbar[132] - 23./4.*armZZbar[17];
   armZZbar[186]=1./3.*armZZbar[186] + 25./8.*armZZbar[18];
   armZZbar[186]=1./3.*armZZbar[186] - 5./8.*armZZbar[37];
   armZZbar[186]=armZZbar[37]*armZZbar[186];
   armZZbar[287]= - 25./24.*armZZbar[10] + 5./8. - 2./9.*armZZbar[11];
   armZZbar[287]=armZZbar[37]*armZZbar[287];
   armZZbar[287]= - 5./8.*armZZbar[50] + armZZbar[287];
   armZZbar[287]=MMH*armZZbar[287];
   armZZbar[186]=armZZbar[186] + armZZbar[287];
   armZZbar[186]=MMH*armZZbar[186];
   armZZbar[143]=1./3.*armZZbar[186] + armZZbar[143];
   armZZbar[117]=armZZbar[143] + armZZbar[117];
   armZZbar[117]=armZZbar[272]*armZZbar[117];
   armZZbar[110]=armZZbar[110] + armZZbar[188];
   armZZbar[110]=armZZbar[10]*armZZbar[110];
   armZZbar[186]=3 + 5*armZZbar[72];
   armZZbar[105]=armZZbar[274] + armZZbar[110] + armZZbar[279] + 
   armZZbar[105] + 3./4.*armZZbar[67] + armZZbar[64] - armZZbar[61] + 1.
   /2.*armZZbar[186] + armZZbar[62];
   armZZbar[105]=MMH*armZZbar[105];
   armZZbar[110]= - 3./4. - armZZbar[38];
   armZZbar[110]=3*armZZbar[110] + armZZbar[181];
   armZZbar[110]=armZZbar[17]*armZZbar[110];
   armZZbar[186]=armZZbar[276] + armZZbar[269];
   armZZbar[186]=armZZbar[18]*armZZbar[186];
   armZZbar[269]=31./2. + 21*armZZbar[36];
   armZZbar[269]=armZZbar[19]*armZZbar[269];
   armZZbar[276]= - 11 - 35*armZZbar[66];
   armZZbar[276]=1./4.*armZZbar[276] + armZZbar[67];
   armZZbar[276]=MMt*armZZbar[276];
   armZZbar[105]=armZZbar[276] + armZZbar[105] + 21./2.*armZZbar[37] + 
   armZZbar[186] + armZZbar[110] + 1./2.*armZZbar[269] + armZZbar[142]
    + armZZbar[286] + 39./4.*armZZbar[22] + armZZbar[285] + 
   armZZbar[328] + 7./4.*armZZbar[68] - 1./4.*armZZbar[69] - 71./4.*
   armZZbar[70] + 9*armZZbar[42];
   armZZbar[105]=MMt*armZZbar[105];
   armZZbar[110]= - 3./2.*armZZbar[40];
   armZZbar[186]= - 1./3.*armZZbar[19];
   armZZbar[152]=armZZbar[186] + armZZbar[170] + 3./2.*armZZbar[20] + 
   armZZbar[152] + armZZbar[110] + armZZbar[69] + armZZbar[68];
   armZZbar[170]=1./9. + 3./4.*armZZbar[38];
   armZZbar[170]=armZZbar[17]*armZZbar[170];
   armZZbar[269]=armZZbar[19] + 13./3.*armZZbar[17];
   armZZbar[269]=armZZbar[35]*armZZbar[269];
   armZZbar[152]=1./12.*armZZbar[269] + 1./2.*armZZbar[152] + 
   armZZbar[170];
   armZZbar[170]= - 1./12.*armZZbar[10];
   armZZbar[269]=armZZbar[170] - 2./3. + armZZbar[118];
   armZZbar[269]=armZZbar[37]*armZZbar[269];
   armZZbar[276]=1./6.*armZZbar[344] + armZZbar[216] + 5./4.*
   armZZbar[15] + 1./4.*armZZbar[67] + armZZbar[64] - 1./4.*
   armZZbar[59] + 3./4. - armZZbar[72];
   armZZbar[285]=59./18. + armZZbar[283];
   armZZbar[285]=1./8.*armZZbar[285] + 2./9.*armZZbar[35];
   armZZbar[285]=armZZbar[10]*armZZbar[285];
   armZZbar[276]=1./4.*armZZbar[276] + armZZbar[285];
   armZZbar[276]=MMH*armZZbar[276];
   armZZbar[285]= - 59./18. + armZZbar[158];
   armZZbar[285]=1./8.*armZZbar[285] - 2./9.*armZZbar[35];
   armZZbar[285]=armZZbar[18]*armZZbar[285];
   armZZbar[152]=armZZbar[276] + armZZbar[269] + 1./2.*armZZbar[152] + 
   armZZbar[285];
   armZZbar[152]=MMH*armZZbar[152];
   armZZbar[269]=31./3.*armZZbar[18] + armZZbar[124] - 9*armZZbar[17];
   armZZbar[269]=1./2.*armZZbar[269] + armZZbar[280];
   armZZbar[269]=armZZbar[37]*armZZbar[269];
   armZZbar[105]=1./2.*armZZbar[105] + 1./4.*armZZbar[269] + 
   armZZbar[152];
   armZZbar[105]=MMt*armZZbar[105];
   armZZbar[152]=armZZbar[18]*armZZbar[254];
   armZZbar[152]=armZZbar[152] + 1./2.*armZZbar[348] - armZZbar[19] + 1.
   /2.*armZZbar[251];
   armZZbar[254]=1./6. - armZZbar[11];
   armZZbar[269]=2./3.*armZZbar[10];
   armZZbar[254]=1./2.*armZZbar[254] + armZZbar[269];
   armZZbar[254]=armZZbar[37]*armZZbar[254];
   armZZbar[276]=armZZbar[35] - 1 + armZZbar[181];
   armZZbar[280]=armZZbar[10]*armZZbar[276];
   armZZbar[285]=armZZbar[280] + armZZbar[11] + 1./2.*armZZbar[344];
   armZZbar[285]=MMH*armZZbar[285];
   armZZbar[152]=1./12.*armZZbar[285] + 1./12.*armZZbar[152] + 
   armZZbar[254];
   armZZbar[152]=MMH*armZZbar[152];
   armZZbar[137]=1./3.*armZZbar[209] + 5*armZZbar[137] + armZZbar[251];
   armZZbar[254]=armZZbar[229] + armZZbar[354];
   armZZbar[254]=MMH*armZZbar[254];
   armZZbar[137]=1./2.*armZZbar[137] + 1./3.*armZZbar[254];
   armZZbar[137]=MMt*armZZbar[137];
   armZZbar[254]=armZZbar[124] - armZZbar[17];
   armZZbar[285]=armZZbar[254] - 7./3.*armZZbar[18];
   armZZbar[285]=armZZbar[37]*armZZbar[285];
   armZZbar[137]=1./2.*armZZbar[137] + 1./4.*armZZbar[285] + 
   armZZbar[152];
   armZZbar[137]=MMt*armZZbar[137];
   armZZbar[152]= - armZZbar[19] + armZZbar[18];
   armZZbar[152]=armZZbar[37]*armZZbar[152];
   armZZbar[285]=armZZbar[11] - armZZbar[10];
   armZZbar[286]=MMH*armZZbar[37]*armZZbar[285];
   armZZbar[152]=armZZbar[152] + armZZbar[286];
   armZZbar[152]=MMH*armZZbar[152];
   armZZbar[137]=1./12.*armZZbar[152] + armZZbar[137];
   armZZbar[137]=armZZbar[133]*armZZbar[137];
   armZZbar[152]=armZZbar[206] - 41./12.*armZZbar[18];
   armZZbar[152]=1./3.*armZZbar[152] + 1./4.*armZZbar[37];
   armZZbar[152]=armZZbar[37]*armZZbar[152];
   armZZbar[206]=41./36.*armZZbar[10] - 1./4. + armZZbar[216];
   armZZbar[206]=armZZbar[37]*armZZbar[206];
   armZZbar[206]=1./4.*armZZbar[50] + armZZbar[206];
   armZZbar[206]=MMH*armZZbar[206];
   armZZbar[152]=armZZbar[152] + armZZbar[206];
   armZZbar[152]=MMH*armZZbar[152];
   armZZbar[105]=armZZbar[137] + 1./4.*armZZbar[152] + armZZbar[105];
   armZZbar[105]=armZZbar[133]*armZZbar[105];
   armZZbar[105]=armZZbar[143] + armZZbar[105];
   armZZbar[105]=armZZbar[133]*armZZbar[105];
   armZZbar[137]=1 + armZZbar[36];
   armZZbar[122]=17*armZZbar[137] + armZZbar[122];
   armZZbar[122]=armZZbar[37]*armZZbar[122];
   armZZbar[143]=armZZbar[233] + 7./2.*armZZbar[200];
   armZZbar[152]= - armZZbar[3]*armZZbar[37]*armZZbar[36];
   armZZbar[143]=1./18.*armZZbar[143] + 3*armZZbar[152];
   armZZbar[143]=MMt*armZZbar[143];
   armZZbar[122]=armZZbar[143] + 1./18.*armZZbar[122] + 3*armZZbar[267]
   ;
   armZZbar[122]=MMt*armZZbar[122];
   armZZbar[122]=17./18.*armZZbar[171] + armZZbar[122];
   armZZbar[122]=MMt*armZZbar[122];
   armZZbar[143]=armZZbar[272]*armZZbar[122];
   armZZbar[122]=armZZbar[122] + armZZbar[143];
   armZZbar[122]=armZZbar[272]*armZZbar[122];
   armZZbar[137]=armZZbar[137] + armZZbar[35];
   armZZbar[137]=armZZbar[37]*armZZbar[137];
   armZZbar[143]=armZZbar[36] + armZZbar[200];
   armZZbar[143]=MMt*armZZbar[143];
   armZZbar[137]=armZZbar[137] + armZZbar[143];
   armZZbar[137]=MMt*armZZbar[137];
   armZZbar[137]=armZZbar[171] + armZZbar[137];
   armZZbar[137]=4./9.*MMt*armZZbar[137];
   armZZbar[122]=armZZbar[137] + armZZbar[122];
   armZZbar[122]=armZZbar[272]*armZZbar[122];
   armZZbar[143]=armZZbar[316] + armZZbar[134];
   armZZbar[143]=armZZbar[37]*armZZbar[143];
   armZZbar[152]= - armZZbar[36] + 1./2.*armZZbar[200];
   armZZbar[200]=armZZbar[3]*armZZbar[37]*armZZbar[36];
   armZZbar[152]=1./2.*armZZbar[152] + 3*armZZbar[200];
   armZZbar[152]=MMt*armZZbar[152];
   armZZbar[143]=armZZbar[152] + 1./2.*armZZbar[143] + 3*armZZbar[278];
   armZZbar[143]=MMt*armZZbar[143];
   armZZbar[143]= - 1./2.*armZZbar[171] + armZZbar[143];
   armZZbar[143]=MMt*armZZbar[143];
   armZZbar[152]=armZZbar[133]*armZZbar[143];
   armZZbar[143]=armZZbar[143] + armZZbar[152];
   armZZbar[143]=armZZbar[133]*armZZbar[143];
   armZZbar[137]=armZZbar[137] + armZZbar[143];
   armZZbar[137]=armZZbar[133]*armZZbar[137];
   armZZbar[122]=armZZbar[122] + armZZbar[137];
   armZZbar[122]=armZZbar[33]*armZZbar[122];
   armZZbar[105]=armZZbar[122] + armZZbar[117] + armZZbar[105];
   armZZbar[105]=armZZbar[33]*armZZbar[105];
   armZZbar[117]= - 5./9.*armZZbar[10];
   armZZbar[122]=armZZbar[117] - 1 + armZZbar[331];
   armZZbar[122]=armZZbar[18]*armZZbar[122];
   armZZbar[137]=1 - 1./18.*armZZbar[11];
   armZZbar[137]=armZZbar[17]*armZZbar[137];
   armZZbar[143]=1./4.*armZZbar[19] + armZZbar[17];
   armZZbar[143]=armZZbar[10]*armZZbar[143];
   armZZbar[152]= - 1./2.*armZZbar[75];
   armZZbar[171]=5./2.*armZZbar[20] - armZZbar[91] + armZZbar[152];
   armZZbar[171]=1./2.*armZZbar[171] + armZZbar[21];
   armZZbar[171]=1./3.*armZZbar[171];
   armZZbar[122]=1./4.*armZZbar[122] + 1./9.*armZZbar[143] + 
   armZZbar[171] + 1./2.*armZZbar[137];
   armZZbar[122]=armZZbar[73]*armZZbar[122];
   armZZbar[137]=1./6.*armZZbar[84];
   armZZbar[143]= - 1./6.*armZZbar[15];
   armZZbar[200]=armZZbar[143] - 1 + armZZbar[137];
   armZZbar[150]=armZZbar[10]*armZZbar[150];
   armZZbar[150]=armZZbar[200] + 1./18.*armZZbar[150];
   armZZbar[150]=MMH*armZZbar[73]*armZZbar[150];
   armZZbar[122]=armZZbar[122] + 1./2.*armZZbar[150];
   armZZbar[122]=MMH*armZZbar[122];
   armZZbar[150]=1./6.*armZZbar[19];
   armZZbar[206]=armZZbar[150] - armZZbar[17];
   armZZbar[206]=armZZbar[17]*armZZbar[206];
   armZZbar[233]= - 29./64.*armZZbar[18] - 1./16.*armZZbar[19] - 
   armZZbar[17];
   armZZbar[233]=armZZbar[18]*armZZbar[233];
   armZZbar[206]=1./8.*armZZbar[206] + 1./3.*armZZbar[233];
   armZZbar[206]=armZZbar[73]*armZZbar[206];
   armZZbar[122]=1./3.*armZZbar[206] + 1./4.*armZZbar[122];
   armZZbar[206]=pow(MMH,2);
   armZZbar[122]=armZZbar[206]*armZZbar[122];
   armZZbar[233]= - armZZbar[17]*armZZbar[11];
   armZZbar[267]=armZZbar[10]*armZZbar[356];
   armZZbar[278]=armZZbar[18]*armZZbar[11];
   armZZbar[233]=armZZbar[278] + armZZbar[233] + armZZbar[267];
   armZZbar[233]=armZZbar[73]*armZZbar[233];
   armZZbar[267]= - MMH*armZZbar[73]*armZZbar[10]*armZZbar[11];
   armZZbar[233]=armZZbar[233] + armZZbar[267];
   armZZbar[233]=MMH*armZZbar[233];
   armZZbar[267]=armZZbar[17]*armZZbar[356];
   armZZbar[278]=armZZbar[18]*armZZbar[359];
   armZZbar[267]=armZZbar[267] + armZZbar[278];
   armZZbar[267]=armZZbar[73]*armZZbar[267];
   armZZbar[233]=armZZbar[267] + armZZbar[233];
   armZZbar[233]=armZZbar[272]*armZZbar[206]*armZZbar[233];
   armZZbar[122]=armZZbar[122] + 1./144.*armZZbar[233];
   armZZbar[122]=armZZbar[272]*armZZbar[122];
   armZZbar[233]=armZZbar[10]*armZZbar[17];
   armZZbar[267]= - 1./4. - 1./9.*armZZbar[10];
   armZZbar[267]=armZZbar[18]*armZZbar[267];
   armZZbar[233]=armZZbar[267] + 1./9.*armZZbar[233] + armZZbar[171] + 
   armZZbar[119];
   armZZbar[233]=armZZbar[73]*armZZbar[233];
   armZZbar[267]=pow(armZZbar[10],2);
   armZZbar[267]=armZZbar[200] + 1./36.*armZZbar[267];
   armZZbar[267]=MMH*armZZbar[73]*armZZbar[267];
   armZZbar[233]=armZZbar[233] + 1./2.*armZZbar[267];
   armZZbar[233]=MMH*armZZbar[233];
   armZZbar[267]= - 2*armZZbar[17];
   armZZbar[278]=armZZbar[267] - 31./32.*armZZbar[18];
   armZZbar[278]=armZZbar[18]*armZZbar[278];
   armZZbar[286]=pow(armZZbar[17],2);
   armZZbar[278]= - 11./16.*armZZbar[286] + armZZbar[278];
   armZZbar[278]=armZZbar[73]*armZZbar[278];
   armZZbar[233]=1./9.*armZZbar[278] + 1./2.*armZZbar[233];
   armZZbar[233]=armZZbar[206]*armZZbar[233];
   armZZbar[122]=armZZbar[233] + armZZbar[122];
   armZZbar[122]=armZZbar[272]*armZZbar[122];
   armZZbar[278]=armZZbar[281] - 1 + armZZbar[138];
   armZZbar[278]=armZZbar[18]*armZZbar[278];
   armZZbar[281]=1 + 1./18.*armZZbar[11];
   armZZbar[281]=armZZbar[17]*armZZbar[281];
   armZZbar[287]= - 1./4.*armZZbar[19];
   armZZbar[289]=armZZbar[287] + armZZbar[17];
   armZZbar[289]=armZZbar[10]*armZZbar[289];
   armZZbar[171]=1./4.*armZZbar[278] + 1./9.*armZZbar[289] + 
   armZZbar[171] + 1./2.*armZZbar[281];
   armZZbar[171]=armZZbar[73]*armZZbar[171];
   armZZbar[278]=armZZbar[10]*armZZbar[11];
   armZZbar[200]=armZZbar[200] + 1./18.*armZZbar[278];
   armZZbar[200]=MMH*armZZbar[73]*armZZbar[200];
   armZZbar[171]=armZZbar[171] + 1./2.*armZZbar[200];
   armZZbar[171]=MMH*armZZbar[171];
   armZZbar[200]= - armZZbar[19] - 5*armZZbar[17];
   armZZbar[200]=armZZbar[17]*armZZbar[200];
   armZZbar[281]=1./16.*armZZbar[19] - armZZbar[17];
   armZZbar[281]=1./3.*armZZbar[281] - 11./64.*armZZbar[18];
   armZZbar[281]=armZZbar[18]*armZZbar[281];
   armZZbar[200]=1./48.*armZZbar[200] + armZZbar[281];
   armZZbar[200]=armZZbar[73]*armZZbar[200];
   armZZbar[171]=1./3.*armZZbar[200] + 1./4.*armZZbar[171];
   armZZbar[171]=armZZbar[206]*armZZbar[171];
   armZZbar[200]= - armZZbar[17]*armZZbar[19];
   armZZbar[281]=armZZbar[346] - armZZbar[18];
   armZZbar[281]=armZZbar[18]*armZZbar[281];
   armZZbar[200]=armZZbar[200] + armZZbar[281];
   armZZbar[200]=armZZbar[73]*armZZbar[200];
   armZZbar[281]=armZZbar[17]*armZZbar[11];
   armZZbar[289]=armZZbar[10]*armZZbar[178];
   armZZbar[289]=armZZbar[281] + armZZbar[289];
   armZZbar[290]=armZZbar[18]*armZZbar[318];
   armZZbar[289]=1./2.*armZZbar[289] + armZZbar[290];
   armZZbar[289]=armZZbar[73]*armZZbar[289];
   armZZbar[285]=MMH*armZZbar[73]*armZZbar[10]*armZZbar[285];
   armZZbar[285]=armZZbar[289] + 1./2.*armZZbar[285];
   armZZbar[285]=MMH*armZZbar[285];
   armZZbar[200]=1./2.*armZZbar[200] + armZZbar[285];
   armZZbar[200]=armZZbar[133]*armZZbar[206]*armZZbar[200];
   armZZbar[171]=armZZbar[171] + 1./72.*armZZbar[200];
   armZZbar[171]=armZZbar[133]*armZZbar[171];
   armZZbar[171]=armZZbar[233] + armZZbar[171];
   armZZbar[171]=armZZbar[133]*armZZbar[171];
   armZZbar[200]=MMH*armZZbar[354];
   armZZbar[200]=armZZbar[200] + armZZbar[264] + armZZbar[209];
   armZZbar[200]=MMt*MMH*armZZbar[200];
   armZZbar[206]=armZZbar[17] - armZZbar[18];
   armZZbar[206]=armZZbar[37]*armZZbar[206];
   armZZbar[209]=MMH*armZZbar[37]*armZZbar[10];
   armZZbar[206]=armZZbar[206] + armZZbar[209];
   armZZbar[206]=MMH*armZZbar[206];
   armZZbar[200]=armZZbar[206] + armZZbar[200];
   armZZbar[200]=MMt*armZZbar[200];
   armZZbar[206]=armZZbar[272]*armZZbar[200];
   armZZbar[200]=armZZbar[200] + armZZbar[206];
   armZZbar[200]=armZZbar[200]*pow(armZZbar[2],4);
   armZZbar[206]=MMH*armZZbar[236];
   armZZbar[206]=armZZbar[206] + armZZbar[251] + armZZbar[352];
   armZZbar[206]=MMt*MMH*armZZbar[206];
   armZZbar[209]= - armZZbar[17] + armZZbar[18];
   armZZbar[209]=armZZbar[37]*armZZbar[209];
   armZZbar[233]= - MMH*armZZbar[37]*armZZbar[10];
   armZZbar[209]=armZZbar[209] + armZZbar[233];
   armZZbar[209]=MMH*armZZbar[209];
   armZZbar[206]=armZZbar[209] + armZZbar[206];
   armZZbar[206]=MMt*armZZbar[206];
   armZZbar[209]=armZZbar[133]*armZZbar[206];
   armZZbar[206]=armZZbar[206] + armZZbar[209];
   armZZbar[206]=armZZbar[206]*pow(armZZbar[5],4);
   armZZbar[200]=armZZbar[200] + armZZbar[206];
   armZZbar[200]=armZZbar[33]*armZZbar[200];
   armZZbar[122]=1./24.*armZZbar[200] + armZZbar[122] + armZZbar[171];
   armZZbar[122]=armZZbar[4]*armZZbar[122];
   armZZbar[171]=1./3.*armZZbar[17];
   armZZbar[200]=armZZbar[287] + armZZbar[171];
   armZZbar[200]=armZZbar[10]*armZZbar[200];
   armZZbar[206]=7 + armZZbar[118];
   armZZbar[206]=armZZbar[19]*armZZbar[206];
   armZZbar[209]=7./3. - armZZbar[11];
   armZZbar[209]=armZZbar[17]*armZZbar[209];
   armZZbar[233]= - armZZbar[11] - 1./2.*armZZbar[10];
   armZZbar[233]=armZZbar[18]*armZZbar[233];
   armZZbar[191]=1./6.*armZZbar[233] + armZZbar[200] + 1./2.*
   armZZbar[209] + 1./3.*armZZbar[206] + armZZbar[191] - armZZbar[92]
    + 1./2.*armZZbar[76];
   armZZbar[191]=armZZbar[73]*armZZbar[191];
   armZZbar[200]=armZZbar[216] - 1 - armZZbar[89];
   armZZbar[200]=1./2.*armZZbar[200] + 1./3.*armZZbar[278];
   armZZbar[200]=MMH*armZZbar[73]*armZZbar[200];
   armZZbar[191]=armZZbar[191] + armZZbar[200];
   armZZbar[191]=MMH*armZZbar[191];
   armZZbar[200]=1./3.*armZZbar[18];
   armZZbar[206]=armZZbar[200] + armZZbar[263] - armZZbar[17];
   armZZbar[206]=armZZbar[18]*armZZbar[206];
   armZZbar[209]=pow(armZZbar[19],2);
   armZZbar[233]= - 1./3.*armZZbar[209];
   armZZbar[236]=11./6.*armZZbar[19] - armZZbar[17];
   armZZbar[236]=armZZbar[17]*armZZbar[236];
   armZZbar[206]=1./4.*armZZbar[206] + armZZbar[233] + 1./2.*
   armZZbar[236];
   armZZbar[206]=armZZbar[73]*armZZbar[206];
   armZZbar[191]=armZZbar[206] + armZZbar[191];
   armZZbar[191]=armZZbar[272]*MMH*armZZbar[191];
   armZZbar[206]=9./4.*armZZbar[43];
   armZZbar[236]=1./2.*armZZbar[47];
   armZZbar[251]=1./4.*armZZbar[45];
   armZZbar[278]=armZZbar[138] + armZZbar[251] + armZZbar[236] + 31./3.
    + armZZbar[206];
   armZZbar[278]=armZZbar[17]*armZZbar[278];
   armZZbar[285]= - 9*armZZbar[43];
   armZZbar[287]=35./18. + armZZbar[285];
   armZZbar[289]= - 1./2.*armZZbar[45];
   armZZbar[287]=armZZbar[289] + 1./2.*armZZbar[287] - armZZbar[47];
   armZZbar[290]=1./6.*armZZbar[10];
   armZZbar[287]=armZZbar[290] + 1./2.*armZZbar[287] + armZZbar[151];
   armZZbar[287]=armZZbar[18]*armZZbar[287];
   armZZbar[291]=47 - 13./4.*armZZbar[11];
   armZZbar[291]=armZZbar[19]*armZZbar[291];
   armZZbar[293]= - 47./2.*armZZbar[19] + armZZbar[17];
   armZZbar[293]=armZZbar[10]*armZZbar[293];
   armZZbar[278]=1./8.*armZZbar[287] + 1./72.*armZZbar[293] + 1./8.*
   armZZbar[278] + 1./18.*armZZbar[291] - 7./96.*armZZbar[21] + 23./48.
   *armZZbar[20] - 5./8.*armZZbar[22] - 1./96.*armZZbar[75] + 1./16.*
   armZZbar[93] + 7./16.*armZZbar[76] + 1./32.*armZZbar[91] - 5./32.*
   armZZbar[74] - armZZbar[92];
   armZZbar[278]=armZZbar[73]*armZZbar[278];
   armZZbar[287]=9./8.*armZZbar[43];
   armZZbar[291]=1./4.*armZZbar[47];
   armZZbar[293]=1./8.*armZZbar[45];
   armZZbar[294]=armZZbar[293] + armZZbar[291] + 1./9. + armZZbar[287];
   armZZbar[294]= - 11./72.*armZZbar[10] + 1./2.*armZZbar[294] + 
   armZZbar[331];
   armZZbar[294]=armZZbar[10]*armZZbar[294];
   armZZbar[295]= - 7./3.*armZZbar[87];
   armZZbar[297]= - 3*armZZbar[86];
   armZZbar[298]=armZZbar[297] - 33./2. + armZZbar[295];
   armZZbar[299]= - 1./2.*armZZbar[96];
   armZZbar[298]=armZZbar[299] + 1./4.*armZZbar[298] - armZZbar[84];
   armZZbar[294]=armZZbar[294] + 1./36.*armZZbar[11] + 17./24.*
   armZZbar[15] + 1./4.*armZZbar[298] - armZZbar[89];
   armZZbar[294]=armZZbar[73]*armZZbar[294];
   armZZbar[298]=armZZbar[79] + 1./3.*armZZbar[78];
   armZZbar[298]=MMH*armZZbar[73]*armZZbar[298];
   armZZbar[294]=armZZbar[294] + 1./8.*armZZbar[298];
   armZZbar[294]=MMH*armZZbar[294];
   armZZbar[278]=armZZbar[278] + 1./2.*armZZbar[294];
   armZZbar[278]=MMH*armZZbar[278];
   armZZbar[294]=armZZbar[17]*armZZbar[19];
   armZZbar[300]= - 41./9.*armZZbar[209] + 3./2.*armZZbar[294];
   armZZbar[301]=3./8.*armZZbar[19] - 5./9.*armZZbar[17];
   armZZbar[301]=1./2.*armZZbar[301] + 2./9.*armZZbar[18];
   armZZbar[301]=armZZbar[18]*armZZbar[301];
   armZZbar[300]=1./8.*armZZbar[300] + armZZbar[301];
   armZZbar[300]=armZZbar[73]*armZZbar[300];
   armZZbar[278]=armZZbar[300] + armZZbar[278];
   armZZbar[278]=MMH*armZZbar[278];
   armZZbar[191]=armZZbar[278] + 1./12.*armZZbar[191];
   armZZbar[191]=armZZbar[272]*armZZbar[191];
   armZZbar[278]=9./2.*armZZbar[43];
   armZZbar[300]=1./2.*armZZbar[45];
   armZZbar[301]=armZZbar[300] + armZZbar[47] + 5./9. + armZZbar[278];
   armZZbar[302]= - 7./18.*armZZbar[10];
   armZZbar[301]=armZZbar[302] + 1./2.*armZZbar[301] + armZZbar[331];
   armZZbar[301]=armZZbar[10]*armZZbar[301];
   armZZbar[295]=armZZbar[297] - 5./2. + armZZbar[295];
   armZZbar[304]=17./6.*armZZbar[15];
   armZZbar[295]=armZZbar[301] + armZZbar[304] - 1./2.*armZZbar[89] + 
   armZZbar[299] + 1./4.*armZZbar[295] - armZZbar[84];
   armZZbar[295]=armZZbar[73]*armZZbar[295];
   armZZbar[295]=armZZbar[295] + 1./2.*armZZbar[298];
   armZZbar[295]=MMH*armZZbar[295];
   armZZbar[301]=1./2.*armZZbar[93];
   armZZbar[309]= - 1./12.*armZZbar[75];
   armZZbar[310]= - 7./12.*armZZbar[21] + 23./6.*armZZbar[20] - 37./6.*
   armZZbar[22] + armZZbar[309] + armZZbar[301] + 1./4.*armZZbar[91] - 
   5./4.*armZZbar[74] - armZZbar[92];
   armZZbar[311]=armZZbar[300] + armZZbar[47] + 59./9. + armZZbar[278];
   armZZbar[311]=1./2.*armZZbar[311] + armZZbar[331];
   armZZbar[311]=armZZbar[17]*armZZbar[311];
   armZZbar[313]=49./4. - 19./3.*armZZbar[11];
   armZZbar[313]=armZZbar[19]*armZZbar[313];
   armZZbar[310]=1./2.*armZZbar[311] + 1./2.*armZZbar[310] + 1./3.*
   armZZbar[313];
   armZZbar[311]=31./18. + armZZbar[285];
   armZZbar[311]=armZZbar[289] + 1./2.*armZZbar[311] - armZZbar[47];
   armZZbar[311]=7./9.*armZZbar[10] + 1./2.*armZZbar[311] + 
   armZZbar[138];
   armZZbar[311]=armZZbar[18]*armZZbar[311];
   armZZbar[313]= - armZZbar[19] - 1./8.*armZZbar[17];
   armZZbar[313]=armZZbar[10]*armZZbar[313];
   armZZbar[310]=1./4.*armZZbar[311] + 1./2.*armZZbar[310] + 1./9.*
   armZZbar[313];
   armZZbar[310]=armZZbar[73]*armZZbar[310];
   armZZbar[295]=armZZbar[310] + 1./4.*armZZbar[295];
   armZZbar[295]=MMH*armZZbar[295];
   armZZbar[310]=3*armZZbar[19];
   armZZbar[171]=armZZbar[310] + armZZbar[171];
   armZZbar[171]=armZZbar[17]*armZZbar[171];
   armZZbar[171]= - 53./9.*armZZbar[209] + armZZbar[171];
   armZZbar[311]= - 1./2.*armZZbar[19] - armZZbar[17];
   armZZbar[311]=5*armZZbar[311] + 17./4.*armZZbar[18];
   armZZbar[311]=armZZbar[18]*armZZbar[311];
   armZZbar[171]=1./4.*armZZbar[171] + 1./9.*armZZbar[311];
   armZZbar[171]=armZZbar[73]*armZZbar[171];
   armZZbar[171]=1./2.*armZZbar[171] + armZZbar[295];
   armZZbar[171]=MMH*armZZbar[171];
   armZZbar[191]=armZZbar[171] + armZZbar[191];
   armZZbar[191]=armZZbar[272]*armZZbar[191];
   armZZbar[295]=1./6.*armZZbar[11] + armZZbar[293] + armZZbar[291] - 
   41./9. + armZZbar[287];
   armZZbar[295]=armZZbar[17]*armZZbar[295];
   armZZbar[311]=1./2. - 3*armZZbar[43];
   armZZbar[311]=armZZbar[289] + 3./2.*armZZbar[311] - armZZbar[47];
   armZZbar[311]=armZZbar[362] + 1./2.*armZZbar[311] + armZZbar[331];
   armZZbar[311]=armZZbar[18]*armZZbar[311];
   armZZbar[313]= - 7*armZZbar[11];
   armZZbar[315]= - 101./9. + armZZbar[313];
   armZZbar[315]=armZZbar[19]*armZZbar[315];
   armZZbar[316]= - 137./4.*armZZbar[19] - armZZbar[17];
   armZZbar[316]=armZZbar[10]*armZZbar[316];
   armZZbar[295]=1./2.*armZZbar[311] + 1./9.*armZZbar[316] + 
   armZZbar[295] + 1./2.*armZZbar[315] - 7./24.*armZZbar[21] + 23./12.*
   armZZbar[20] - 17./3.*armZZbar[22] - 1./24.*armZZbar[75] + 1./4.*
   armZZbar[93] - 37./12.*armZZbar[76] + 1./8.*armZZbar[91] - 5./8.*
   armZZbar[74] + 17./3.*armZZbar[92];
   armZZbar[295]=armZZbar[73]*armZZbar[295];
   armZZbar[311]=19./2. - armZZbar[87];
   armZZbar[311]=7./3.*armZZbar[311] + armZZbar[297];
   armZZbar[304]=armZZbar[138] + armZZbar[304] + 17./3.*armZZbar[89] + 
   armZZbar[299] + 1./4.*armZZbar[311] - armZZbar[84];
   armZZbar[311]=armZZbar[251] + armZZbar[236] + 1./3. + armZZbar[206];
   armZZbar[170]=armZZbar[170] + 1./2.*armZZbar[311] + armZZbar[138];
   armZZbar[170]=armZZbar[10]*armZZbar[170];
   armZZbar[170]=1./2.*armZZbar[304] + armZZbar[170];
   armZZbar[170]=armZZbar[73]*armZZbar[170];
   armZZbar[170]=armZZbar[170] + 1./4.*armZZbar[298];
   armZZbar[170]=MMH*armZZbar[170];
   armZZbar[170]=armZZbar[295] + armZZbar[170];
   armZZbar[170]=MMH*armZZbar[170];
   armZZbar[295]= - 23./2.*armZZbar[19] + armZZbar[17];
   armZZbar[295]=armZZbar[17]*armZZbar[295];
   armZZbar[298]=121./2.*armZZbar[19] + armZZbar[18];
   armZZbar[298]=armZZbar[18]*armZZbar[298];
   armZZbar[295]=1./24.*armZZbar[298] - 2*armZZbar[209] + 1./8.*
   armZZbar[295];
   armZZbar[295]=armZZbar[73]*armZZbar[295];
   armZZbar[170]=1./3.*armZZbar[295] + 1./4.*armZZbar[170];
   armZZbar[170]=MMH*armZZbar[170];
   armZZbar[295]= - 11*armZZbar[17];
   armZZbar[298]=19*armZZbar[18] - 49*armZZbar[19] + armZZbar[295];
   armZZbar[298]=armZZbar[18]*armZZbar[298];
   armZZbar[298]=1./6.*armZZbar[298] + 5*armZZbar[209] + 11./6.*
   armZZbar[294];
   armZZbar[298]=armZZbar[73]*armZZbar[298];
   armZZbar[304]=1./3. - 5*armZZbar[11];
   armZZbar[304]=armZZbar[19]*armZZbar[304];
   armZZbar[311]= - 1./3.*armZZbar[17];
   armZZbar[315]=17./2.*armZZbar[19] + armZZbar[311];
   armZZbar[315]=armZZbar[10]*armZZbar[315];
   armZZbar[281]=armZZbar[315] + armZZbar[304] + 1./3.*armZZbar[281];
   armZZbar[304]= - 1./3. + armZZbar[11];
   armZZbar[304]=1./3.*armZZbar[304] - 3./2.*armZZbar[10];
   armZZbar[304]=armZZbar[18]*armZZbar[304];
   armZZbar[281]=1./3.*armZZbar[281] + armZZbar[304];
   armZZbar[281]=armZZbar[73]*armZZbar[281];
   armZZbar[304]=armZZbar[10] + 1./4. - armZZbar[11];
   armZZbar[304]=armZZbar[10]*armZZbar[304];
   armZZbar[304]= - 1./4.*armZZbar[11] + armZZbar[304];
   armZZbar[304]=MMH*armZZbar[73]*armZZbar[304];
   armZZbar[281]=1./4.*armZZbar[281] + 1./9.*armZZbar[304];
   armZZbar[281]=MMH*armZZbar[281];
   armZZbar[281]=1./12.*armZZbar[298] + armZZbar[281];
   armZZbar[281]=armZZbar[133]*MMH*armZZbar[281];
   armZZbar[170]=armZZbar[170] + 1./2.*armZZbar[281];
   armZZbar[170]=armZZbar[133]*armZZbar[170];
   armZZbar[170]=armZZbar[171] + armZZbar[170];
   armZZbar[170]=armZZbar[133]*armZZbar[170];
   armZZbar[105]=armZZbar[122] + armZZbar[105] + armZZbar[191] + 
   armZZbar[170];
   armZZbar[105]=armZZbar[4]*armZZbar[105];
   armZZbar[122]= - 739./64. - armZZbar[88];
   armZZbar[170]= - 5./4.*armZZbar[97];
   armZZbar[171]= - 1./6.*armZZbar[87];
   armZZbar[191]=7./4.*armZZbar[86];
   armZZbar[122]=armZZbar[191] + armZZbar[171] - 13./4.*armZZbar[90] + 
   armZZbar[170] + 1./3.*armZZbar[122] + 1./4.*armZZbar[85];
   armZZbar[281]=53./27. + armZZbar[285];
   armZZbar[298]= - 3./4.*armZZbar[45];
   armZZbar[304]=13./18.*armZZbar[12];
   armZZbar[315]= - 7./18.*armZZbar[11];
   armZZbar[281]=armZZbar[315] + 19./36.*armZZbar[13] + armZZbar[304]
    + armZZbar[298] + 1./4.*armZZbar[281] - armZZbar[47];
   armZZbar[316]=4./9.*armZZbar[10];
   armZZbar[281]=1./2.*armZZbar[281] + armZZbar[316];
   armZZbar[281]=armZZbar[10]*armZZbar[281];
   armZZbar[317]=2./3.*armZZbar[84];
   armZZbar[318]=7./24.*armZZbar[96];
   armZZbar[319]=19./24.*armZZbar[16];
   armZZbar[320]=pow(armZZbar[13],2);
   armZZbar[321]=1./12.*armZZbar[320];
   armZZbar[322]=49./36. + armZZbar[103];
   armZZbar[322]=1./3.*armZZbar[11]*armZZbar[322];
   armZZbar[122]=armZZbar[281] + armZZbar[322] + armZZbar[321] + 
   armZZbar[293] - 15./8.*armZZbar[15] + armZZbar[291] + armZZbar[287]
    + armZZbar[319] + 3./8.*armZZbar[89] + armZZbar[318] + 1./2.*
   armZZbar[122] + armZZbar[317];
   armZZbar[122]=armZZbar[73]*armZZbar[122];
   armZZbar[281]=armZZbar[1]*armZZbar[11];
   armZZbar[323]=armZZbar[39]*armZZbar[11];
   armZZbar[324]=armZZbar[281] + 1./3.*armZZbar[323];
   armZZbar[325]=1./3. - armZZbar[7];
   armZZbar[326]=armZZbar[1]*armZZbar[325];
   armZZbar[325]=armZZbar[39]*armZZbar[325];
   armZZbar[325]=1./3.*armZZbar[326] + armZZbar[325];
   armZZbar[325]=armZZbar[10]*armZZbar[325];
   armZZbar[326]= - armZZbar[1]*armZZbar[11];
   armZZbar[328]= - armZZbar[39]*armZZbar[11];
   armZZbar[329]=armZZbar[326] + 1./3.*armZZbar[328];
   armZZbar[329]=armZZbar[6]*armZZbar[329];
   armZZbar[324]=1./3.*armZZbar[329] + 1./9.*armZZbar[324] + 1./2.*
   armZZbar[325];
   armZZbar[324]=1./6.*armZZbar[324];
   armZZbar[329]=armZZbar[80] - 1./2.*armZZbar[81];
   armZZbar[329]=1./2.*armZZbar[329] + armZZbar[78];
   armZZbar[329]=1./3.*MMH*armZZbar[73]*armZZbar[329];
   armZZbar[122]=armZZbar[329] + armZZbar[324] + armZZbar[122];
   armZZbar[122]=MMH*armZZbar[122];
   armZZbar[333]=9*armZZbar[43];
   armZZbar[335]= - armZZbar[45] + 151./27. + armZZbar[333];
   armZZbar[336]=13./9.*armZZbar[12];
   armZZbar[335]=armZZbar[216] - 37./9.*armZZbar[13] + 1./2.*
   armZZbar[335] + armZZbar[336];
   armZZbar[335]=armZZbar[17]*armZZbar[335];
   armZZbar[337]= - 13./9.*armZZbar[12];
   armZZbar[338]=5./9.*armZZbar[11];
   armZZbar[339]=armZZbar[302] + armZZbar[338] - 7./18.*armZZbar[13] + 
   armZZbar[337] + armZZbar[45] - 329./54. + armZZbar[47];
   armZZbar[339]=armZZbar[18]*armZZbar[339];
   armZZbar[340]=11*armZZbar[92];
   armZZbar[342]= - 3*armZZbar[74];
   armZZbar[343]=1./6.*armZZbar[77] + armZZbar[342] + armZZbar[340];
   armZZbar[344]=9*armZZbar[13];
   armZZbar[348]= - 1847./27. + armZZbar[344];
   armZZbar[350]=53./9.*armZZbar[11];
   armZZbar[348]=1./2.*armZZbar[348] + armZZbar[350];
   armZZbar[348]=armZZbar[19]*armZZbar[348];
   armZZbar[352]= - 7./12.*armZZbar[75];
   armZZbar[354]=19./4.*armZZbar[20];
   armZZbar[358]=25./6.*armZZbar[19] + armZZbar[17];
   armZZbar[358]=1./3.*armZZbar[10]*armZZbar[358];
   armZZbar[335]=1./2.*armZZbar[339] + armZZbar[358] + 1./2.*
   armZZbar[335] + 1./2.*armZZbar[348] + 37./12.*armZZbar[21] + 
   armZZbar[354] + 103./12.*armZZbar[22] + armZZbar[352] - armZZbar[93]
    - 29./12.*armZZbar[76] + 1./2.*armZZbar[343] - armZZbar[91];
   armZZbar[335]=armZZbar[73]*armZZbar[335];
   armZZbar[339]=1 - armZZbar[7];
   armZZbar[339]=armZZbar[17]*armZZbar[339];
   armZZbar[339]=armZZbar[186] + 1./2.*armZZbar[339];
   armZZbar[339]=armZZbar[1]*armZZbar[339];
   armZZbar[343]= - 1./3. + armZZbar[7];
   armZZbar[348]=armZZbar[1]*armZZbar[343];
   armZZbar[343]=armZZbar[39]*armZZbar[343];
   armZZbar[343]=1./3.*armZZbar[348] + armZZbar[343];
   armZZbar[348]=armZZbar[18]*armZZbar[343];
   armZZbar[360]=armZZbar[1]*armZZbar[356];
   armZZbar[356]=armZZbar[39]*armZZbar[356];
   armZZbar[361]=armZZbar[360] + 1./3.*armZZbar[356];
   armZZbar[361]=armZZbar[6]*armZZbar[361];
   armZZbar[362]=11./27. - armZZbar[7];
   armZZbar[362]=armZZbar[17]*armZZbar[362];
   armZZbar[362]= - 1./27.*armZZbar[19] + 1./2.*armZZbar[362];
   armZZbar[362]=armZZbar[39]*armZZbar[362];
   armZZbar[339]=1./3.*armZZbar[361] + 1./2.*armZZbar[348] + 1./3.*
   armZZbar[339] + armZZbar[362];
   armZZbar[339]=1./3.*armZZbar[339];
   armZZbar[335]=armZZbar[339] + armZZbar[335];
   armZZbar[122]=1./2.*armZZbar[335] + armZZbar[122];
   armZZbar[122]=MMH*armZZbar[122];
   armZZbar[335]= - 3./2.*armZZbar[45];
   armZZbar[348]=1./9.*armZZbar[13] - 17./18.*armZZbar[12] + 
   armZZbar[335] - armZZbar[47] - 29./27. - 9./2.*armZZbar[43];
   armZZbar[361]= - 7./9.*armZZbar[11];
   armZZbar[348]=13./9.*armZZbar[10] + 1./2.*armZZbar[348] + 
   armZZbar[361];
   armZZbar[348]=armZZbar[10]*armZZbar[348];
   armZZbar[362]=7./2.*armZZbar[86] - 1./3.*armZZbar[87] - armZZbar[90]
    - 2653./288. + 5*armZZbar[97];
   armZZbar[363]= - 5./9. + armZZbar[240];
   armZZbar[364]=armZZbar[11]*armZZbar[363];
   armZZbar[365]=1./3.*armZZbar[84];
   armZZbar[366]=9./16.*armZZbar[43];
   armZZbar[367]=1./8.*armZZbar[47];
   armZZbar[368]=1./16.*armZZbar[45];
   armZZbar[320]= - 1./24.*armZZbar[320];
   armZZbar[348]=1./4.*armZZbar[348] + 1./24.*armZZbar[364] + 
   armZZbar[320] + armZZbar[368] - 41./48.*armZZbar[15] + armZZbar[367]
    + armZZbar[366] - 5./8.*armZZbar[16] + 37./72.*armZZbar[89] + 1./8.
   *armZZbar[96] + 1./8.*armZZbar[362] + armZZbar[365];
   armZZbar[348]=armZZbar[73]*armZZbar[348];
   armZZbar[362]= - armZZbar[14] + armZZbar[216];
   armZZbar[364]=armZZbar[1]*armZZbar[362];
   armZZbar[362]=armZZbar[39]*armZZbar[362];
   armZZbar[369]=armZZbar[10]*armZZbar[230];
   armZZbar[370]=1./3.*armZZbar[369];
   armZZbar[362]=armZZbar[370] + armZZbar[364] + 11./9.*armZZbar[362];
   armZZbar[364]=armZZbar[326] + 11./9.*armZZbar[328];
   armZZbar[369]=1./2.*armZZbar[364] + armZZbar[369];
   armZZbar[369]=armZZbar[6]*armZZbar[369];
   armZZbar[362]=1./2.*armZZbar[362] + armZZbar[369];
   armZZbar[369]=1./4.*armZZbar[81];
   armZZbar[371]=armZZbar[369] + armZZbar[78];
   armZZbar[371]=MMH*armZZbar[73]*armZZbar[371];
   armZZbar[348]=1./6.*armZZbar[371] + 1./6.*armZZbar[362] + 
   armZZbar[348];
   armZZbar[348]=MMH*armZZbar[348];
   armZZbar[362]=167./9. - 29./2.*armZZbar[13];
   armZZbar[341]=1./2.*armZZbar[362] + armZZbar[341];
   armZZbar[341]=armZZbar[19]*armZZbar[341];
   armZZbar[362]=77./36.*armZZbar[13] - 17./36.*armZZbar[12] - 1./4.*
   armZZbar[45] + armZZbar[236] + 97./27. + armZZbar[206];
   armZZbar[362]=1./2.*armZZbar[362] + armZZbar[216];
   armZZbar[362]=armZZbar[17]*armZZbar[362];
   armZZbar[371]=armZZbar[342] + 13./9.*armZZbar[92];
   armZZbar[371]= - 7./6.*armZZbar[76] - armZZbar[91] + 1./2.*
   armZZbar[371] + 1./3.*armZZbar[77];
   armZZbar[372]=armZZbar[10]*armZZbar[19];
   armZZbar[373]=25./9.*armZZbar[10] - 5./9.*armZZbar[11] + 5./9.*
   armZZbar[13] + 17./18.*armZZbar[12] - 83./54. + armZZbar[45];
   armZZbar[373]=armZZbar[18]*armZZbar[373];
   armZZbar[374]= - 1./3.*armZZbar[93];
   armZZbar[197]=1./4.*armZZbar[373] + 23./9.*armZZbar[372] + 
   armZZbar[362] + 1./6.*armZZbar[341] + armZZbar[197] + 115./24.*
   armZZbar[20] + 25./36.*armZZbar[22] - 7./24.*armZZbar[75] + 1./2.*
   armZZbar[371] + armZZbar[374];
   armZZbar[197]=armZZbar[73]*armZZbar[197];
   armZZbar[341]= - armZZbar[1] - 11./9.*armZZbar[39];
   armZZbar[362]=armZZbar[18]*armZZbar[341];
   armZZbar[371]=armZZbar[1]*armZZbar[178];
   armZZbar[178]=armZZbar[39]*armZZbar[178];
   armZZbar[373]=armZZbar[362] + armZZbar[371] + 11./9.*armZZbar[178];
   armZZbar[375]=armZZbar[1]*armZZbar[346];
   armZZbar[346]=armZZbar[39]*armZZbar[346];
   armZZbar[346]=armZZbar[375] + 11./9.*armZZbar[346];
   armZZbar[346]=1./2.*armZZbar[346] + armZZbar[362];
   armZZbar[346]=armZZbar[6]*armZZbar[346];
   armZZbar[346]=1./6.*armZZbar[373] + armZZbar[346];
   armZZbar[197]=1./3.*armZZbar[346] + armZZbar[197];
   armZZbar[197]=1./2.*armZZbar[197] + armZZbar[348];
   armZZbar[197]=MMH*armZZbar[197];
   armZZbar[346]=armZZbar[281] + 11./9.*armZZbar[323];
   armZZbar[348]=armZZbar[6]*armZZbar[364];
   armZZbar[364]=armZZbar[11]*armZZbar[246];
   armZZbar[364]=1./3.*armZZbar[364] - armZZbar[16] + armZZbar[89] + 1./
   4. + armZZbar[97];
   armZZbar[373]= - 1./12.*armZZbar[12] - armZZbar[11];
   armZZbar[373]=armZZbar[10]*armZZbar[373];
   armZZbar[364]=1./2.*armZZbar[364] + armZZbar[373];
   armZZbar[364]=armZZbar[73]*armZZbar[364];
   armZZbar[346]=armZZbar[364] + 1./3.*armZZbar[346] + armZZbar[348];
   armZZbar[346]=MMH*armZZbar[346];
   armZZbar[348]=armZZbar[1]*armZZbar[359];
   armZZbar[359]=armZZbar[39]*armZZbar[359];
   armZZbar[348]=armZZbar[348] + 11./9.*armZZbar[359];
   armZZbar[356]=armZZbar[360] + 11./9.*armZZbar[356];
   armZZbar[356]=armZZbar[6]*armZZbar[356];
   armZZbar[348]=1./3.*armZZbar[348] + armZZbar[356];
   armZZbar[356]=armZZbar[263] + armZZbar[131];
   armZZbar[356]=armZZbar[10]*armZZbar[356];
   armZZbar[246]=armZZbar[19]*armZZbar[246];
   armZZbar[246]=5./3.*armZZbar[246] - armZZbar[22] + armZZbar[20];
   armZZbar[359]=5./16.*armZZbar[13] - 1./3. - 1./16.*armZZbar[12];
   armZZbar[359]=armZZbar[17]*armZZbar[359];
   armZZbar[246]=1./2.*armZZbar[356] + 1./8.*armZZbar[246] + 1./3.*
   armZZbar[359];
   armZZbar[356]=1 + 1./2.*armZZbar[12];
   armZZbar[356]=1./9.*armZZbar[356] - armZZbar[11];
   armZZbar[356]=1./2.*armZZbar[356] + 1./9.*armZZbar[10];
   armZZbar[356]=armZZbar[18]*armZZbar[356];
   armZZbar[246]=1./3.*armZZbar[246] + 1./4.*armZZbar[356];
   armZZbar[246]=armZZbar[73]*armZZbar[246];
   armZZbar[246]=1./12.*armZZbar[346] + 1./12.*armZZbar[348] + 
   armZZbar[246];
   armZZbar[246]=MMH*armZZbar[246];
   armZZbar[346]= - armZZbar[17] - 1./12.*armZZbar[18];
   armZZbar[346]=armZZbar[18]*armZZbar[346];
   armZZbar[233]=1./2.*armZZbar[346] + armZZbar[233] + 1./8.*
   armZZbar[294];
   armZZbar[233]=armZZbar[73]*armZZbar[233];
   armZZbar[233]=1./3.*armZZbar[233] + armZZbar[246];
   armZZbar[233]=armZZbar[272]*armZZbar[233];
   armZZbar[246]= - 91./18.*armZZbar[19] - armZZbar[17];
   armZZbar[246]=armZZbar[17]*armZZbar[246];
   armZZbar[346]= - 97./6.*armZZbar[18] - 235./6.*armZZbar[19] - 7*
   armZZbar[17];
   armZZbar[346]=armZZbar[18]*armZZbar[346];
   armZZbar[246]=1./3.*armZZbar[346] - 233./18.*armZZbar[209] + 
   armZZbar[246];
   armZZbar[246]=armZZbar[73]*armZZbar[246];
   armZZbar[197]=armZZbar[233] + 1./4.*armZZbar[246] + armZZbar[197];
   armZZbar[197]=armZZbar[272]*armZZbar[197];
   armZZbar[233]= - 1453./3.*armZZbar[19] + 7*armZZbar[17];
   armZZbar[246]= - 17*armZZbar[18];
   armZZbar[233]=1./3.*armZZbar[233] + armZZbar[246];
   armZZbar[233]=armZZbar[18]*armZZbar[233];
   armZZbar[346]=29./6.*armZZbar[19] - armZZbar[17];
   armZZbar[346]=armZZbar[17]*armZZbar[346];
   armZZbar[233]=1./4.*armZZbar[233] - 1471./36.*armZZbar[209] + 
   armZZbar[346];
   armZZbar[233]=armZZbar[73]*armZZbar[233];
   armZZbar[122]=armZZbar[197] + 1./2.*armZZbar[233] + armZZbar[122];
   armZZbar[122]=armZZbar[272]*armZZbar[122];
   armZZbar[197]= - 1./3.*armZZbar[88];
   armZZbar[170]=armZZbar[191] + armZZbar[171] - 43./12.*armZZbar[90]
    + armZZbar[170] + 5./4.*armZZbar[85] + 31./64. + armZZbar[197];
   armZZbar[171]=305./27. + armZZbar[285];
   armZZbar[171]=armZZbar[315] + 67./36.*armZZbar[13] + armZZbar[304]
    + armZZbar[298] + 1./4.*armZZbar[171] - armZZbar[47];
   armZZbar[171]=1./2.*armZZbar[171] + armZZbar[316];
   armZZbar[171]=armZZbar[10]*armZZbar[171];
   armZZbar[170]=armZZbar[171] + armZZbar[322] + armZZbar[321] + 
   armZZbar[293] - 19./8.*armZZbar[15] + armZZbar[291] + armZZbar[287]
    + armZZbar[319] + 89./24.*armZZbar[89] + armZZbar[318] + 1./2.*
   armZZbar[170] + armZZbar[317];
   armZZbar[170]=armZZbar[73]*armZZbar[170];
   armZZbar[170]=armZZbar[329] + armZZbar[324] + armZZbar[170];
   armZZbar[170]=MMH*armZZbar[170];
   armZZbar[171]= - armZZbar[45] - 929./27. + armZZbar[333];
   armZZbar[171]=armZZbar[216] - 13./9.*armZZbar[13] + 1./2.*
   armZZbar[171] + armZZbar[336];
   armZZbar[171]=armZZbar[17]*armZZbar[171];
   armZZbar[191]=armZZbar[302] + armZZbar[338] - 55./18.*armZZbar[13]
    + armZZbar[337] + armZZbar[45] - 581./54. + armZZbar[47];
   armZZbar[191]=armZZbar[18]*armZZbar[191];
   armZZbar[233]= - 4079./27. + armZZbar[344];
   armZZbar[233]=1./2.*armZZbar[233] + armZZbar[350];
   armZZbar[233]=armZZbar[19]*armZZbar[233];
   armZZbar[293]= - 11./6.*armZZbar[77] + armZZbar[342] + 89./3.*
   armZZbar[92];
   armZZbar[171]=1./2.*armZZbar[191] + armZZbar[358] + 1./2.*
   armZZbar[171] + 1./2.*armZZbar[233] + 49./12.*armZZbar[21] + 
   armZZbar[354] + 79./12.*armZZbar[22] + armZZbar[352] - armZZbar[93]
    - 61./12.*armZZbar[76] + 1./2.*armZZbar[293] - armZZbar[91];
   armZZbar[171]=armZZbar[73]*armZZbar[171];
   armZZbar[171]=armZZbar[339] + armZZbar[171];
   armZZbar[170]=1./2.*armZZbar[171] + armZZbar[170];
   armZZbar[170]=MMH*armZZbar[170];
   armZZbar[171]= - 9./8.*armZZbar[43];
   armZZbar[191]= - 3./4.*armZZbar[47];
   armZZbar[233]= - 3./8.*armZZbar[45];
   armZZbar[290]=armZZbar[290] - 23./6.*armZZbar[13] + 55./24.*
   armZZbar[12] + armZZbar[233] + armZZbar[191] - 101./27. + 
   armZZbar[171];
   armZZbar[290]=armZZbar[10]*armZZbar[290];
   armZZbar[197]=7./8.*armZZbar[86] - 1./12.*armZZbar[87] - 13./6.*
   armZZbar[97] - 13./12.*armZZbar[85] - 1101./128. + armZZbar[197];
   armZZbar[293]= - 11./2.*armZZbar[13];
   armZZbar[302]= - 1./9. + armZZbar[293];
   armZZbar[302]=armZZbar[11]*armZZbar[302];
   armZZbar[197]=1./2.*armZZbar[290] + 5./24.*armZZbar[302] + 
   armZZbar[320] + armZZbar[368] - 17./48.*armZZbar[15] + armZZbar[367]
    + armZZbar[366] + 5./4.*armZZbar[16] - 85./12.*armZZbar[89] + 1./6.
   *armZZbar[96] + 1./2.*armZZbar[197] + armZZbar[365];
   armZZbar[197]=armZZbar[73]*armZZbar[197];
   armZZbar[290]= - armZZbar[14] + armZZbar[151];
   armZZbar[302]=armZZbar[1]*armZZbar[290];
   armZZbar[290]=armZZbar[39]*armZZbar[290];
   armZZbar[290]=1./3.*armZZbar[302] + armZZbar[290];
   armZZbar[302]= - 1./2.*armZZbar[7];
   armZZbar[304]=1./3. + armZZbar[302];
   armZZbar[315]=armZZbar[1]*armZZbar[304];
   armZZbar[304]=armZZbar[39]*armZZbar[304];
   armZZbar[304]=1./3.*armZZbar[315] + armZZbar[304];
   armZZbar[304]=armZZbar[10]*armZZbar[304];
   armZZbar[281]=1./3.*armZZbar[281] + armZZbar[323];
   armZZbar[281]=1./2.*armZZbar[281];
   armZZbar[316]=1./3.*armZZbar[1] + armZZbar[39];
   armZZbar[317]=armZZbar[10]*armZZbar[316];
   armZZbar[318]=armZZbar[281] + armZZbar[317];
   armZZbar[318]=armZZbar[6]*armZZbar[318];
   armZZbar[290]=armZZbar[318] + 1./2.*armZZbar[290] + armZZbar[304];
   armZZbar[304]=armZZbar[78] + armZZbar[80] + armZZbar[369];
   armZZbar[304]=MMH*armZZbar[73]*armZZbar[304];
   armZZbar[197]=1./6.*armZZbar[304] + 1./6.*armZZbar[290] + 
   armZZbar[197];
   armZZbar[197]=MMH*armZZbar[197];
   armZZbar[290]= - 211./6.*armZZbar[13] + 55./6.*armZZbar[12] + 
   armZZbar[289] - armZZbar[47] + 2255./27. + armZZbar[278];
   armZZbar[290]=1./2.*armZZbar[290] + armZZbar[151];
   armZZbar[290]=armZZbar[17]*armZZbar[290];
   armZZbar[304]=armZZbar[215] + armZZbar[300] + 397./108. + 
   armZZbar[47];
   armZZbar[304]= - 4./9.*armZZbar[10] + 5./24.*armZZbar[11] + 1./4.*
   armZZbar[304] + 2*armZZbar[13];
   armZZbar[304]=armZZbar[18]*armZZbar[304];
   armZZbar[318]=5*armZZbar[77] + armZZbar[342] - 169./3.*armZZbar[92];
   armZZbar[318]=25./3.*armZZbar[76] + 1./2.*armZZbar[318] - 
   armZZbar[91];
   armZZbar[173]=1123./9. + armZZbar[173];
   armZZbar[173]=1./2.*armZZbar[173] + 19*armZZbar[11];
   armZZbar[173]=armZZbar[19]*armZZbar[173];
   armZZbar[319]=49./2.*armZZbar[19] + armZZbar[17];
   armZZbar[319]=armZZbar[10]*armZZbar[319];
   armZZbar[173]=armZZbar[304] + 1./6.*armZZbar[319] + 1./4.*
   armZZbar[290] + 1./6.*armZZbar[173] + armZZbar[142] - 17./48.*
   armZZbar[20] + 28./3.*armZZbar[22] - 7./48.*armZZbar[75] + 1./4.*
   armZZbar[318] + armZZbar[374];
   armZZbar[173]=armZZbar[73]*armZZbar[173];
   armZZbar[290]= - 1./3. + armZZbar[302];
   armZZbar[290]=armZZbar[17]*armZZbar[290];
   armZZbar[150]=armZZbar[150] + armZZbar[290];
   armZZbar[290]=armZZbar[1]*armZZbar[150];
   armZZbar[150]=armZZbar[39]*armZZbar[150];
   armZZbar[304]=1./2.*armZZbar[7];
   armZZbar[318]= - 1./3. + armZZbar[304];
   armZZbar[319]=armZZbar[1]*armZZbar[318];
   armZZbar[318]=armZZbar[39]*armZZbar[318];
   armZZbar[318]=1./3.*armZZbar[319] + armZZbar[318];
   armZZbar[318]=armZZbar[18]*armZZbar[318];
   armZZbar[150]=armZZbar[318] + 1./3.*armZZbar[290] + armZZbar[150];
   armZZbar[186]=armZZbar[186] + armZZbar[17];
   armZZbar[290]=armZZbar[1]*armZZbar[186];
   armZZbar[186]=armZZbar[39]*armZZbar[186];
   armZZbar[186]=1./3.*armZZbar[290] + armZZbar[186];
   armZZbar[290]=armZZbar[18]*armZZbar[139];
   armZZbar[186]=1./2.*armZZbar[186] + 1./3.*armZZbar[290];
   armZZbar[186]=armZZbar[6]*armZZbar[186];
   armZZbar[150]=1./3.*armZZbar[150] + armZZbar[186];
   armZZbar[150]=armZZbar[197] + 1./2.*armZZbar[150] + armZZbar[173];
   armZZbar[150]=MMH*armZZbar[150];
   armZZbar[173]=armZZbar[17]*armZZbar[7];
   armZZbar[186]=armZZbar[263] + armZZbar[173];
   armZZbar[197]=armZZbar[1]*armZZbar[186];
   armZZbar[186]=armZZbar[39]*armZZbar[186];
   armZZbar[263]= - 1./3. - armZZbar[7];
   armZZbar[318]=armZZbar[1]*armZZbar[263];
   armZZbar[319]=armZZbar[39]*armZZbar[263];
   armZZbar[318]=1./3.*armZZbar[318] + armZZbar[319];
   armZZbar[319]=armZZbar[18]*armZZbar[318];
   armZZbar[186]=armZZbar[319] + 1./3.*armZZbar[197] + armZZbar[186];
   armZZbar[178]=1./3.*armZZbar[371] + armZZbar[178];
   armZZbar[197]=armZZbar[18]*armZZbar[316];
   armZZbar[178]=1./2.*armZZbar[178] + armZZbar[197];
   armZZbar[178]=armZZbar[6]*armZZbar[178];
   armZZbar[178]=1./2.*armZZbar[186] + armZZbar[178];
   armZZbar[186]= - armZZbar[12] + armZZbar[13];
   armZZbar[316]=11./4.*armZZbar[186] + armZZbar[151];
   armZZbar[316]=armZZbar[17]*armZZbar[316];
   armZZbar[319]=13./12.*armZZbar[11] + 1./27. + 11./16.*armZZbar[13];
   armZZbar[319]=armZZbar[19]*armZZbar[319];
   armZZbar[320]= - 59./3.*armZZbar[19] + armZZbar[17];
   armZZbar[320]=armZZbar[10]*armZZbar[320];
   armZZbar[321]=41./72.*armZZbar[10] - 1./72.*armZZbar[11] - 11./8.*
   armZZbar[13] - 1./27. + 11./16.*armZZbar[12];
   armZZbar[321]=armZZbar[18]*armZZbar[321];
   armZZbar[316]=armZZbar[321] + 1./12.*armZZbar[320] + armZZbar[319]
    + 1./4.*armZZbar[316];
   armZZbar[316]=armZZbar[73]*armZZbar[316];
   armZZbar[319]=1./3.*armZZbar[326] + armZZbar[328];
   armZZbar[320]=armZZbar[1]*armZZbar[141];
   armZZbar[141]=armZZbar[39]*armZZbar[141];
   armZZbar[321]=1./3.*armZZbar[320] + armZZbar[141];
   armZZbar[321]=armZZbar[10]*armZZbar[321];
   armZZbar[319]=1./3.*armZZbar[319] + armZZbar[321];
   armZZbar[322]=armZZbar[10]*armZZbar[139];
   armZZbar[281]=armZZbar[281] + armZZbar[322];
   armZZbar[281]=armZZbar[6]*armZZbar[281];
   armZZbar[281]=1./2.*armZZbar[319] + armZZbar[281];
   armZZbar[319]= - 11./2.*armZZbar[12];
   armZZbar[323]=11*armZZbar[13];
   armZZbar[324]=armZZbar[323] + 47./27. + armZZbar[319];
   armZZbar[117]=armZZbar[117] + 1./4.*armZZbar[324] + armZZbar[338];
   armZZbar[117]=armZZbar[10]*armZZbar[117];
   armZZbar[324]= - 47./27. + armZZbar[293];
   armZZbar[324]=armZZbar[11]*armZZbar[324];
   armZZbar[117]=1./4.*armZZbar[324] + armZZbar[117];
   armZZbar[117]=armZZbar[73]*armZZbar[117];
   armZZbar[117]=1./3.*armZZbar[281] + armZZbar[117];
   armZZbar[117]=MMH*armZZbar[117];
   armZZbar[117]=1./2.*armZZbar[117] + 1./6.*armZZbar[178] + 
   armZZbar[316];
   armZZbar[117]=MMH*armZZbar[117];
   armZZbar[178]= - 91*armZZbar[209] + 17*armZZbar[294];
   armZZbar[281]= - 17*armZZbar[17];
   armZZbar[294]=281./3.*armZZbar[19] + armZZbar[281];
   armZZbar[222]=1./8.*armZZbar[294] + armZZbar[222];
   armZZbar[222]=armZZbar[18]*armZZbar[222];
   armZZbar[178]=1./8.*armZZbar[178] + armZZbar[222];
   armZZbar[178]=armZZbar[73]*armZZbar[178];
   armZZbar[117]=1./3.*armZZbar[178] + armZZbar[117];
   armZZbar[117]=armZZbar[133]*armZZbar[117];
   armZZbar[178]= - 7./2.*armZZbar[19] - armZZbar[17];
   armZZbar[178]=armZZbar[17]*armZZbar[178];
   armZZbar[178]=1679./3.*armZZbar[209] + armZZbar[178];
   armZZbar[121]=2633./3.*armZZbar[19] + armZZbar[121];
   armZZbar[121]=1./4.*armZZbar[121] - 7./9.*armZZbar[18];
   armZZbar[121]=armZZbar[18]*armZZbar[121];
   armZZbar[121]=1./4.*armZZbar[178] + armZZbar[121];
   armZZbar[121]=armZZbar[73]*armZZbar[121];
   armZZbar[117]=armZZbar[117] + armZZbar[121] + armZZbar[150];
   armZZbar[117]=armZZbar[133]*armZZbar[117];
   armZZbar[121]= - 2059./3.*armZZbar[19] + armZZbar[17];
   armZZbar[121]=7./3.*armZZbar[121] + armZZbar[246];
   armZZbar[121]=armZZbar[18]*armZZbar[121];
   armZZbar[121]=1./4.*armZZbar[121] - 10687./36.*armZZbar[209] + 
   armZZbar[346];
   armZZbar[121]=armZZbar[73]*armZZbar[121];
   armZZbar[117]=armZZbar[117] + 1./2.*armZZbar[121] + armZZbar[170];
   armZZbar[117]=armZZbar[133]*armZZbar[117];
   armZZbar[121]=armZZbar[18]*armZZbar[19];
   armZZbar[121]=armZZbar[209] + armZZbar[121];
   armZZbar[121]=armZZbar[73]*armZZbar[121];
   armZZbar[104]=armZZbar[105] + armZZbar[104] + armZZbar[117] + 32*
   armZZbar[121] + armZZbar[122];
   armZZbar[104]=armZZbar[4]*armZZbar[104];
   armZZbar[105]= - 22*armZZbar[46];
   armZZbar[117]= - 6*armZZbar[45];
   armZZbar[121]= - 26*armZZbar[12];
   armZZbar[122]= - 4./3. + armZZbar[7];
   armZZbar[122]=2*armZZbar[1]*armZZbar[122];
   armZZbar[150]=3*armZZbar[7];
   armZZbar[170]= - 20./9. + armZZbar[150];
   armZZbar[170]=2*armZZbar[39]*armZZbar[170];
   armZZbar[178]= - 12*armZZbar[10];
   armZZbar[222]=3*armZZbar[1] + 11./3.*armZZbar[39];
   armZZbar[222]=2*armZZbar[6]*armZZbar[222];
   armZZbar[246]= - 3*armZZbar[48];
   armZZbar[294]= - 6*armZZbar[47];
   armZZbar[316]=6*armZZbar[11];
   armZZbar[324]=armZZbar[222] + armZZbar[178] + armZZbar[170] + 
   armZZbar[122] + armZZbar[316] - 17./2.*armZZbar[13] + armZZbar[121]
    + armZZbar[117] + armZZbar[294] + armZZbar[246] - 95./2. + 
   armZZbar[105];
   armZZbar[324]=armZZbar[37]*armZZbar[324];
   armZZbar[326]= - 9*armZZbar[38];
   armZZbar[328]=armZZbar[275] + 7./3. + armZZbar[326];
   armZZbar[328]=armZZbar[19]*armZZbar[328];
   armZZbar[329]= - 13*armZZbar[17];
   armZZbar[124]=armZZbar[124] + armZZbar[329];
   armZZbar[124]=armZZbar[35]*armZZbar[124];
   armZZbar[338]= - 11*armZZbar[21] - 13./6.*armZZbar[20] - 3*
   armZZbar[22] - 25*armZZbar[50] + 3*armZZbar[42] + 11*armZZbar[41];
   armZZbar[339]=13./6.*armZZbar[17];
   armZZbar[344]= - 205./6.*armZZbar[35] + 3./2.*armZZbar[36] - 65./3.
    + armZZbar[326];
   armZZbar[344]=1./2.*armZZbar[18]*armZZbar[344];
   armZZbar[124]=armZZbar[324] + armZZbar[344] + 1./6.*armZZbar[124] + 
   armZZbar[339] + armZZbar[338] + 1./2.*armZZbar[328];
   armZZbar[124]=armZZbar[3]*armZZbar[124];
   armZZbar[324]=157./108. - armZZbar[60];
   armZZbar[324]=7*armZZbar[324] - 271./8.*armZZbar[63];
   armZZbar[328]= - 22*armZZbar[72];
   armZZbar[346]=19./8.*armZZbar[61];
   armZZbar[348]=11./2.*armZZbar[59];
   armZZbar[350]= - 23*armZZbar[64];
   armZZbar[352]= - 61./6.*armZZbar[67];
   armZZbar[324]=101./32.*armZZbar[48] + 1699./72.*armZZbar[46] + 
   armZZbar[352] + armZZbar[350] + armZZbar[348] + 7./4.*armZZbar[16]
    + 5771./288.*armZZbar[49] - 181./192.*armZZbar[66] + armZZbar[346]
    - 313./96.*armZZbar[65] - 809./54.*armZZbar[62] + 1./4.*
   armZZbar[324] + armZZbar[328];
   armZZbar[354]= - 3*armZZbar[15];
   armZZbar[356]=armZZbar[283] + armZZbar[354] - 11./3.*armZZbar[67] - 
   20./3. + armZZbar[308];
   armZZbar[356]=2*armZZbar[101]*armZZbar[356];
   armZZbar[358]= - 2*armZZbar[58] - armZZbar[56];
   armZZbar[358]=56*armZZbar[358] + 997./4.*armZZbar[55];
   armZZbar[359]= - 44*armZZbar[53];
   armZZbar[358]=armZZbar[359] - armZZbar[57] - armZZbar[54] + 1./27.*
   armZZbar[358] - 19./2.*armZZbar[51];
   armZZbar[360]=armZZbar[19]*armZZbar[38];
   armZZbar[364]=armZZbar[18]*armZZbar[38];
   armZZbar[366]=armZZbar[360] + armZZbar[364];
   armZZbar[366]=armZZbar[3]*armZZbar[366];
   armZZbar[144]=armZZbar[10]*armZZbar[144];
   armZZbar[366]=18*armZZbar[366] + 6*armZZbar[144] - 5./2.*
   armZZbar[36] + 6*armZZbar[38] - 3./2.*armZZbar[48] + 151./3. - 11*
   armZZbar[46];
   armZZbar[366]=armZZbar[3]*armZZbar[366];
   armZZbar[367]= - armZZbar[10]*armZZbar[101]*armZZbar[38];
   armZZbar[368]=6*armZZbar[367];
   armZZbar[358]=armZZbar[366] + armZZbar[368] + 1./3.*armZZbar[358] + 
   armZZbar[356];
   armZZbar[358]=MMt*armZZbar[358];
   armZZbar[369]=4213./3. + 511*armZZbar[46];
   armZZbar[369]=1./9.*armZZbar[369] - 1./4.*armZZbar[48];
   armZZbar[369]=1./3.*armZZbar[369] + armZZbar[221];
   armZZbar[371]=11./3.*armZZbar[239];
   armZZbar[369]=armZZbar[371] + armZZbar[361] + 301./54.*armZZbar[13]
    + 1./4.*armZZbar[369] + 29./9.*armZZbar[12];
   armZZbar[369]=1./2.*armZZbar[369] + 61./81.*armZZbar[35];
   armZZbar[369]=armZZbar[35]*armZZbar[369];
   armZZbar[373]=284./3. + armZZbar[247];
   armZZbar[374]= - 7*armZZbar[7];
   armZZbar[375]=487./3. + armZZbar[374];
   armZZbar[375]=armZZbar[35]*armZZbar[375];
   armZZbar[375]=1./18.*armZZbar[375] + 1./9.*armZZbar[373] + 
   armZZbar[181];
   armZZbar[375]=armZZbar[1]*armZZbar[375];
   armZZbar[376]=1132./27. + armZZbar[247];
   armZZbar[377]=1823./27. + armZZbar[374];
   armZZbar[377]=armZZbar[35]*armZZbar[377];
   armZZbar[378]= - 11./6.*armZZbar[36];
   armZZbar[377]=1./2.*armZZbar[377] + armZZbar[376] + armZZbar[378];
   armZZbar[377]=armZZbar[39]*armZZbar[377];
   armZZbar[379]= - 7./4. + armZZbar[184];
   armZZbar[380]= - 80./9.*armZZbar[35] - 89./9. + armZZbar[229];
   armZZbar[381]=armZZbar[1]*armZZbar[380];
   armZZbar[380]=armZZbar[39]*armZZbar[380];
   armZZbar[379]=11./9.*armZZbar[380] + 1./3.*armZZbar[379] + 
   armZZbar[381];
   armZZbar[379]=armZZbar[6]*armZZbar[379];
   armZZbar[380]=3*armZZbar[69];
   armZZbar[381]=armZZbar[380] - 11./3.*armZZbar[68];
   armZZbar[382]=55./6. + armZZbar[283];
   armZZbar[382]=armZZbar[17]*armZZbar[382];
   armZZbar[381]=armZZbar[382] - 29./3.*armZZbar[21] + 29./6.*
   armZZbar[20] - 6*armZZbar[50] + 11./3.*armZZbar[41] + 2*
   armZZbar[381] - 3*armZZbar[40];
   armZZbar[381]=armZZbar[101]*armZZbar[381];
   armZZbar[382]= - 11./3.*armZZbar[35] + armZZbar[229] + 1./6. + 
   armZZbar[283];
   armZZbar[382]=armZZbar[10]*armZZbar[382];
   armZZbar[383]= - 1 - 2*armZZbar[38];
   armZZbar[383]=armZZbar[101]*armZZbar[383];
   armZZbar[237]=3*armZZbar[383] + armZZbar[237];
   armZZbar[237]=armZZbar[18]*armZZbar[237];
   armZZbar[383]=5./8.*armZZbar[44];
   armZZbar[384]=11./2.*armZZbar[15];
   armZZbar[303]=8./3.*armZZbar[303];
   armZZbar[385]=10./3.*MMH*armZZbar[53];
   armZZbar[124]=armZZbar[358] + armZZbar[385] + armZZbar[124] + 
   armZZbar[303] + armZZbar[379] + armZZbar[237] + armZZbar[382] + 1./9.
   *armZZbar[377] + 1./3.*armZZbar[375] + armZZbar[369] + armZZbar[381]
    + armZZbar[334] + armZZbar[242] + 967./108.*armZZbar[13] + 
   armZZbar[228] + armZZbar[158] + armZZbar[384] + 1./3.*armZZbar[324]
    + armZZbar[383];
   armZZbar[124]=MMt*armZZbar[124];
   armZZbar[324]=29 + 7*armZZbar[46];
   armZZbar[334]=armZZbar[3]*armZZbar[364];
   armZZbar[324]=9*armZZbar[334] + 3*armZZbar[144] + armZZbar[229] + 1./
   6.*armZZbar[324] + armZZbar[283];
   armZZbar[324]=armZZbar[3]*armZZbar[324];
   armZZbar[334]=armZZbar[283] + armZZbar[354] + 7./9.*armZZbar[67] - 
   20./9. + armZZbar[308];
   armZZbar[334]=armZZbar[101]*armZZbar[334];
   armZZbar[219]=armZZbar[219] + 17./3.*armZZbar[54] + 437./108.*
   armZZbar[55] + armZZbar[52];
   armZZbar[219]=1./2.*armZZbar[219] + 14./3.*armZZbar[53];
   armZZbar[358]=3*armZZbar[367];
   armZZbar[219]=armZZbar[324] + armZZbar[358] + 1./3.*armZZbar[219] + 
   armZZbar[334];
   armZZbar[219]=MMt*armZZbar[219];
   armZZbar[324]= - 35./18. + armZZbar[283];
   armZZbar[324]=armZZbar[17]*armZZbar[324];
   armZZbar[334]= - 3*armZZbar[50];
   armZZbar[324]=1./2.*armZZbar[324] - 47./18.*armZZbar[21] + 47./36.*
   armZZbar[20] + armZZbar[334] - 7./18.*armZZbar[41] + armZZbar[110]
    + armZZbar[380] + 7./9.*armZZbar[68];
   armZZbar[324]=armZZbar[101]*armZZbar[324];
   armZZbar[367]= - 5./2. - armZZbar[13];
   armZZbar[312]=1./4.*armZZbar[367] + armZZbar[312];
   armZZbar[367]=7./9.*armZZbar[35];
   armZZbar[369]=armZZbar[367] + 25./9. + armZZbar[229];
   armZZbar[375]=armZZbar[1]*armZZbar[369];
   armZZbar[377]=armZZbar[39]*armZZbar[369];
   armZZbar[312]=11./9.*armZZbar[377] + 1./3.*armZZbar[312] + 
   armZZbar[375];
   armZZbar[312]=armZZbar[6]*armZZbar[312];
   armZZbar[375]= - 3*armZZbar[1] - 11./3.*armZZbar[39];
   armZZbar[377]=armZZbar[6]*armZZbar[375];
   armZZbar[379]= - 3*armZZbar[45];
   armZZbar[377]=2*armZZbar[377] - 9*armZZbar[10] + 22./9.*armZZbar[39]
    + 2*armZZbar[1] + armZZbar[240] + 17./2.*armZZbar[12] + 
   armZZbar[379] + 4 + 7./3.*armZZbar[46];
   armZZbar[377]=armZZbar[37]*armZZbar[377];
   armZZbar[112]=armZZbar[112] + armZZbar[50];
   armZZbar[112]=41./12.*armZZbar[327] - 41./12.*armZZbar[17] + 7./2.*
   armZZbar[21] + 7*armZZbar[112] + 41./12.*armZZbar[20];
   armZZbar[327]=3*armZZbar[36];
   armZZbar[386]=7*armZZbar[35] + armZZbar[327] + 23./3. + 
   armZZbar[326];
   armZZbar[386]=armZZbar[18]*armZZbar[386];
   armZZbar[112]=armZZbar[377] + 1./3.*armZZbar[112] + 1./4.*
   armZZbar[386];
   armZZbar[112]=armZZbar[3]*armZZbar[112];
   armZZbar[377]=1013 - 49*armZZbar[46];
   armZZbar[204]=armZZbar[204] + 1./9.*armZZbar[377] + armZZbar[238];
   armZZbar[204]=1./2.*armZZbar[204] - 119./9.*armZZbar[12];
   armZZbar[204]=7./3.*armZZbar[116] + 1./2.*armZZbar[204] + 43./9.*
   armZZbar[13];
   armZZbar[204]=1./2.*armZZbar[204] + 157./27.*armZZbar[35];
   armZZbar[204]=armZZbar[35]*armZZbar[204];
   armZZbar[238]=armZZbar[36] + 31./6. + armZZbar[283];
   armZZbar[238]=1./2.*armZZbar[238] + armZZbar[367];
   armZZbar[238]=armZZbar[10]*armZZbar[238];
   armZZbar[367]= - 1./2. - armZZbar[38];
   armZZbar[367]=3*armZZbar[101]*armZZbar[367];
   armZZbar[164]=armZZbar[367] + armZZbar[164];
   armZZbar[164]=armZZbar[18]*armZZbar[164];
   armZZbar[377]= - 331./144.*armZZbar[63] + 41./72. - armZZbar[60];
   armZZbar[386]=1./4.*armZZbar[16];
   armZZbar[387]= - 25./9. + armZZbar[181];
   armZZbar[388]=armZZbar[387] - 7./9.*armZZbar[35];
   armZZbar[389]=armZZbar[1]*armZZbar[388];
   armZZbar[388]=armZZbar[39]*armZZbar[388];
   armZZbar[390]= - 3./2.*armZZbar[38];
   armZZbar[112]=armZZbar[219] + 25./9.*armZZbar[274] + armZZbar[112]
    + 68./9.*armZZbar[288] + armZZbar[312] + armZZbar[164] + 
   armZZbar[238] + 11./27.*armZZbar[388] + 1./3.*armZZbar[389] + 1./6.*
   armZZbar[204] + armZZbar[324] + 1./12.*armZZbar[345] + 13./27.*
   armZZbar[13] - 289./108.*armZZbar[12] + armZZbar[390] - 7./12.*
   armZZbar[15] - 25./48.*armZZbar[44] - 25./192.*armZZbar[48] - 301./
   432.*armZZbar[46] + 227./108.*armZZbar[67] + 31./18.*armZZbar[64] - 
   7./36.*armZZbar[59] + armZZbar[386] + 5071./5184.*armZZbar[49] + 257.
   /1152.*armZZbar[66] - 143./144.*armZZbar[61] + 191./324.*
   armZZbar[62] + 1./4.*armZZbar[377] + 7./9.*armZZbar[72];
   armZZbar[112]=MMt*armZZbar[112];
   armZZbar[164]=5*armZZbar[28] - 677./72. - 5*armZZbar[30];
   armZZbar[164]= - 19./18.*armZZbar[15] + 25./24.*armZZbar[44] - 125./
   36.*armZZbar[67] - 523./108.*armZZbar[64] + armZZbar[157] + 37./8.*
   armZZbar[49] + 5./54.*armZZbar[14] + 1./9.*armZZbar[164] - 37./8.*
   armZZbar[61];
   armZZbar[204]= - 17./9. - 25./16.*armZZbar[44];
   armZZbar[204]=1./2.*armZZbar[204] + armZZbar[347];
   armZZbar[204]=armZZbar[35]*armZZbar[204];
   armZZbar[219]=3./2. - 1./3.*armZZbar[44];
   armZZbar[219]=armZZbar[99]*armZZbar[219];
   armZZbar[219]=1./4.*armZZbar[219] + 1./9.*armZZbar[225];
   armZZbar[219]=armZZbar[37]*armZZbar[219];
   armZZbar[238]= - 7 + armZZbar[59];
   armZZbar[143]=armZZbar[143] + 1./6.*armZZbar[67] + 1./6.*
   armZZbar[238] - armZZbar[64];
   armZZbar[143]=armZZbar[99]*armZZbar[143];
   armZZbar[143]=armZZbar[143] + 1./6.*armZZbar[225];
   armZZbar[143]=MMH*armZZbar[143];
   armZZbar[312]= - armZZbar[68] + armZZbar[330];
   armZZbar[324]=1./12.*armZZbar[21] - 13./6.*armZZbar[50] + 
   armZZbar[312] - 1./12.*armZZbar[41];
   armZZbar[324]=armZZbar[99]*armZZbar[324];
   armZZbar[330]=49./6. - 187*armZZbar[35];
   armZZbar[330]=armZZbar[10]*armZZbar[330];
   armZZbar[345]=armZZbar[218] - 11./4.*armZZbar[10];
   armZZbar[345]=armZZbar[6]*armZZbar[345];
   armZZbar[347]= - 11./81.*armZZbar[11];
   armZZbar[143]=25./96.*armZZbar[143] + 25./4.*armZZbar[219] + 5./54.*
   armZZbar[345] + 25./576.*armZZbar[351] + 1./216.*armZZbar[330] + 1./
   6.*armZZbar[204] + 25./48.*armZZbar[277] + 25./48.*armZZbar[324] + 1.
   /8.*armZZbar[164] + armZZbar[347];
   armZZbar[143]=MMH*armZZbar[143];
   armZZbar[164]= - 3229./6. - 35*armZZbar[46];
   armZZbar[164]=armZZbar[205] + 1./9.*armZZbar[164] + 5./4.*
   armZZbar[48];
   armZZbar[164]=563./108.*armZZbar[35] + 25./9.*armZZbar[13] + 5./4.*
   armZZbar[164] + armZZbar[245];
   armZZbar[204]=107*armZZbar[34] - 2825./6.*armZZbar[99];
   armZZbar[204]=armZZbar[18]*armZZbar[204];
   armZZbar[205]=275./9.*armZZbar[39] + 1./4. + 25*armZZbar[1];
   armZZbar[205]=armZZbar[6]*armZZbar[205];
   armZZbar[219]=107./4.*armZZbar[34] - 350./3.*armZZbar[99];
   armZZbar[219]=armZZbar[37]*armZZbar[219];
   armZZbar[164]=1./3.*armZZbar[219] + 1./3.*armZZbar[205] + 1./12.*
   armZZbar[204] + armZZbar[284] - 275./81.*armZZbar[39] + 1./4.*
   armZZbar[164] - 25./9.*armZZbar[1];
   armZZbar[164]=armZZbar[37]*armZZbar[164];
   armZZbar[192]=1741./72.*armZZbar[50] + 1./3.*armZZbar[41] - 25./4.*
   armZZbar[40] + 259./18.*armZZbar[68] - 215./27.*armZZbar[8] + 25./8.
   *armZZbar[192] + 59./9.*armZZbar[69];
   armZZbar[204]= - 635./6.*armZZbar[19] + armZZbar[281];
   armZZbar[204]=armZZbar[35]*armZZbar[204];
   armZZbar[205]=571 - 233*armZZbar[35];
   armZZbar[205]=armZZbar[18]*armZZbar[205];
   armZZbar[219]=5*armZZbar[17];
   armZZbar[245]= - 167./3.*armZZbar[19] + armZZbar[219];
   armZZbar[245]=1./8.*armZZbar[245] - 25./9.*armZZbar[18];
   armZZbar[245]=armZZbar[6]*armZZbar[245];
   armZZbar[164]=armZZbar[164] + 1./3.*armZZbar[245] + 1./288.*
   armZZbar[205] + 1./12.*armZZbar[204] + 2./3.*armZZbar[17] + 7661./
   1728.*armZZbar[19] - 311./216.*armZZbar[21] - 73./24.*armZZbar[20]
    + 1./8.*armZZbar[192] - 5./3.*armZZbar[22];
   armZZbar[192]=11./3. - 17./2.*armZZbar[35];
   armZZbar[192]=armZZbar[18]*armZZbar[192];
   armZZbar[204]=armZZbar[6]*armZZbar[357];
   armZZbar[205]= - armZZbar[37]*armZZbar[12];
   armZZbar[192]=17*armZZbar[205] + 5./2.*armZZbar[204] + armZZbar[192]
    - 11./3.*armZZbar[19] + 17./2.*armZZbar[353];
   armZZbar[204]= - 17*armZZbar[12];
   armZZbar[205]=armZZbar[204] + 7./2.*armZZbar[268];
   armZZbar[245]=1 + armZZbar[12];
   armZZbar[245]=armZZbar[3]*armZZbar[37]*armZZbar[245];
   armZZbar[268]=MMt*armZZbar[3]*armZZbar[36];
   armZZbar[205]=armZZbar[268] + 1./54.*armZZbar[205] + armZZbar[245];
   armZZbar[205]=MMt*armZZbar[205];
   armZZbar[192]=1./54.*armZZbar[192] + armZZbar[205];
   armZZbar[192]=armZZbar[272]*armZZbar[192];
   armZZbar[205]=7*armZZbar[18] - 17./6.*armZZbar[37];
   armZZbar[205]=armZZbar[3]*armZZbar[37]*armZZbar[205];
   armZZbar[112]=1./2.*armZZbar[192] + armZZbar[112] + armZZbar[143] + 
   1./3.*armZZbar[164] + armZZbar[205];
   armZZbar[112]=armZZbar[272]*armZZbar[112];
   armZZbar[143]=armZZbar[28] - 1085./216. - armZZbar[30];
   armZZbar[143]= - 5./6.*armZZbar[15] - 5./8.*armZZbar[44] + 25./12.*
   armZZbar[67] - 1./36.*armZZbar[64] + armZZbar[157] - 3./8.*
   armZZbar[49] + 1./18.*armZZbar[14] + 1./3.*armZZbar[143] + 3./8.*
   armZZbar[61];
   armZZbar[164]=23./27. + armZZbar[383];
   armZZbar[164]=1./2.*armZZbar[164] + armZZbar[314];
   armZZbar[164]=armZZbar[35]*armZZbar[164];
   armZZbar[192]= - 9./2. + armZZbar[44];
   armZZbar[192]=armZZbar[99]*armZZbar[192];
   armZZbar[168]=1./4.*armZZbar[192] + armZZbar[168];
   armZZbar[168]=armZZbar[37]*armZZbar[168];
   armZZbar[192]=7 - armZZbar[59];
   armZZbar[192]=1./6.*armZZbar[15] - 1./6.*armZZbar[67] + 1./6.*
   armZZbar[192] + armZZbar[64];
   armZZbar[192]=armZZbar[99]*armZZbar[192];
   armZZbar[175]=armZZbar[192] + 1./6.*armZZbar[175];
   armZZbar[175]=MMH*armZZbar[175];
   armZZbar[192]= - 1./12.*armZZbar[21] + 13./6.*armZZbar[50] + 1./12.*
   armZZbar[41] + armZZbar[68] + armZZbar[252];
   armZZbar[192]=armZZbar[99]*armZZbar[192];
   armZZbar[205]=armZZbar[355] + 7./12.*armZZbar[35];
   armZZbar[205]=armZZbar[10]*armZZbar[205];
   armZZbar[245]= - 1./2. + 5*armZZbar[11];
   armZZbar[245]=1./3.*armZZbar[245] - 7./2.*armZZbar[10];
   armZZbar[245]=armZZbar[6]*armZZbar[245];
   armZZbar[143]=5./16.*armZZbar[175] + 5./2.*armZZbar[168] + 1./18.*
   armZZbar[245] + 5./96.*armZZbar[349] + 1./3.*armZZbar[205] + 1./2.*
   armZZbar[164] + 5./8.*armZZbar[213] + armZZbar[279] + 5./8.*
   armZZbar[192] + 1./4.*armZZbar[143] + armZZbar[347];
   armZZbar[143]=MMH*armZZbar[143];
   armZZbar[164]=28877./6. + 1105*armZZbar[46];
   armZZbar[164]=1./9.*armZZbar[164] + armZZbar[226];
   armZZbar[164]=1./3.*armZZbar[164] + armZZbar[221];
   armZZbar[168]=1./2.*armZZbar[98] - 11*armZZbar[34];
   armZZbar[168]=armZZbar[19]*armZZbar[168];
   armZZbar[175]=armZZbar[1]*armZZbar[373];
   armZZbar[192]=armZZbar[39]*armZZbar[376];
   armZZbar[203]=armZZbar[203] + 6695./6.*armZZbar[99];
   armZZbar[203]=armZZbar[18]*armZZbar[203];
   armZZbar[205]= - 979./9.*armZZbar[39] + 5./2. - 89*armZZbar[1];
   armZZbar[205]=armZZbar[6]*armZZbar[205];
   armZZbar[213]= - 173./2.*armZZbar[34] + 1420./3.*armZZbar[99];
   armZZbar[213]=armZZbar[37]*armZZbar[213];
   armZZbar[226]= - 25./6.*armZZbar[10];
   armZZbar[164]=1./9.*armZZbar[213] + 1./9.*armZZbar[205] + 1./18.*
   armZZbar[203] + armZZbar[226] + 1./9.*armZZbar[192] + 1./27.*
   armZZbar[175] + 851./648.*armZZbar[35] + 1./6.*armZZbar[168] + 
   armZZbar[242] + 481./108.*armZZbar[13] + 1./8.*armZZbar[164] + 
   armZZbar[228];
   armZZbar[164]=armZZbar[37]*armZZbar[164];
   armZZbar[168]= - 11*armZZbar[18];
   armZZbar[175]=13*armZZbar[19];
   armZZbar[192]=armZZbar[175] + armZZbar[168];
   armZZbar[203]= - 15*armZZbar[37];
   armZZbar[192]=2./3.*armZZbar[192] + armZZbar[203];
   armZZbar[192]=armZZbar[3]*armZZbar[37]*armZZbar[192];
   armZZbar[205]= - 47./8.*armZZbar[68] + 83./36.*armZZbar[8] + 139./12.
   *armZZbar[69] - 425./64.*armZZbar[42] - 68./3.*armZZbar[71] - 151./
   32.*armZZbar[70];
   armZZbar[213]= - 169./48. - armZZbar[36];
   armZZbar[213]=armZZbar[19]*armZZbar[213];
   armZZbar[228]= - 11./4.*armZZbar[17];
   armZZbar[245]=37*armZZbar[19] + armZZbar[228];
   armZZbar[245]=armZZbar[35]*armZZbar[245];
   armZZbar[252]=2131./72.*armZZbar[35] + 7049./216. - armZZbar[36];
   armZZbar[252]=armZZbar[18]*armZZbar[252];
   armZZbar[268]=armZZbar[254] + 281./18.*armZZbar[18];
   armZZbar[268]=armZZbar[6]*armZZbar[268];
   armZZbar[281]=5./16.*armZZbar[40];
   armZZbar[284]=7./12.*armZZbar[41];
   armZZbar[314]=85./54. - armZZbar[36];
   armZZbar[314]=1./2.*armZZbar[17]*armZZbar[314];
   armZZbar[324]=1./4.*armZZbar[20];
   armZZbar[112]=armZZbar[112] + armZZbar[124] + armZZbar[143] + 
   armZZbar[192] + armZZbar[164] + 1./18.*armZZbar[268] + 1./6.*
   armZZbar[252] + 1./9.*armZZbar[245] + armZZbar[314] + 1./6.*
   armZZbar[213] - 469./324.*armZZbar[21] + armZZbar[324] + 
   armZZbar[22] - 2099./2592.*armZZbar[50] + armZZbar[284] + 1./9.*
   armZZbar[205] + armZZbar[281];
   armZZbar[112]=armZZbar[272]*armZZbar[112];
   armZZbar[105]=armZZbar[222] + armZZbar[178] + armZZbar[170] + 
   armZZbar[122] + armZZbar[316] - 113./2.*armZZbar[13] + armZZbar[121]
    + armZZbar[117] + armZZbar[294] + armZZbar[246] - 191./2. + 
   armZZbar[105];
   armZZbar[105]=armZZbar[37]*armZZbar[105];
   armZZbar[117]=armZZbar[275] - 121./3. + armZZbar[326];
   armZZbar[117]=armZZbar[19]*armZZbar[117];
   armZZbar[121]= - 121*armZZbar[19] + armZZbar[329];
   armZZbar[121]=armZZbar[35]*armZZbar[121];
   armZZbar[105]=armZZbar[105] + armZZbar[344] + 1./6.*armZZbar[121] + 
   armZZbar[339] + armZZbar[338] + 1./2.*armZZbar[117];
   armZZbar[105]=armZZbar[3]*armZZbar[105];
   armZZbar[117]=113./8.*armZZbar[63] + 30167./108. - 43*armZZbar[60];
   armZZbar[117]=armZZbar[352] + armZZbar[350] + armZZbar[348] + 43./4.
   *armZZbar[16] + 89./96.*armZZbar[49] - 5671./64.*armZZbar[66] + 
   armZZbar[346] - 19./32.*armZZbar[65] - 11./2.*armZZbar[62] + 1./4.*
   armZZbar[117] + armZZbar[328];
   armZZbar[121]=2*armZZbar[58] + armZZbar[56];
   armZZbar[122]=armZZbar[359] - armZZbar[57] - armZZbar[54] - 35./2.*
   armZZbar[51] + 8./3.*armZZbar[121] - 1./4.*armZZbar[55];
   armZZbar[122]=armZZbar[366] + armZZbar[368] + 1./3.*armZZbar[122] + 
   armZZbar[356];
   armZZbar[122]=MMt*armZZbar[122];
   armZZbar[124]=armZZbar[221] - 43./4.*armZZbar[48] + 10453./81. - 19*
   armZZbar[46];
   armZZbar[124]=armZZbar[371] + armZZbar[361] + 581./6.*armZZbar[13]
    + 1./4.*armZZbar[124] + 31./3.*armZZbar[12];
   armZZbar[124]=1./2.*armZZbar[124] + armZZbar[282];
   armZZbar[124]=armZZbar[35]*armZZbar[124];
   armZZbar[164]=92./3. + armZZbar[247];
   armZZbar[170]=103./3. + armZZbar[374];
   armZZbar[170]=armZZbar[35]*armZZbar[170];
   armZZbar[178]=1./18.*armZZbar[170] + 1./9.*armZZbar[164] + 
   armZZbar[181];
   armZZbar[178]=armZZbar[1]*armZZbar[178];
   armZZbar[192]=76./3. + armZZbar[247];
   armZZbar[170]=1./2.*armZZbar[170] + armZZbar[192] + armZZbar[378];
   armZZbar[170]=armZZbar[39]*armZZbar[170];
   armZZbar[103]=armZZbar[154] - 7./4. + armZZbar[103];
   armZZbar[154]= - 59 + 11./2.*armZZbar[36];
   armZZbar[154]=1./3.*armZZbar[154] + armZZbar[115];
   armZZbar[154]=armZZbar[39]*armZZbar[154];
   armZZbar[205]= - 16./9.*armZZbar[35] - 25./9. + armZZbar[229];
   armZZbar[205]=armZZbar[1]*armZZbar[205];
   armZZbar[103]=1./3.*armZZbar[154] + 1./3.*armZZbar[103] + 
   armZZbar[205];
   armZZbar[103]=armZZbar[6]*armZZbar[103];
   armZZbar[154]=armZZbar[36]*armZZbar[363];
   armZZbar[103]=armZZbar[122] + armZZbar[385] + armZZbar[105] + 
   armZZbar[303] + armZZbar[103] + armZZbar[237] + armZZbar[382] + 1./9.
   *armZZbar[170] + 1./3.*armZZbar[178] + armZZbar[124] + armZZbar[381]
    + 23./4.*armZZbar[154] + armZZbar[242] + 863./12.*armZZbar[13] + 
   armZZbar[241] + armZZbar[158] + armZZbar[384] + armZZbar[383] - 9./
   32.*armZZbar[48] + 1./3.*armZZbar[117] + 25./8.*armZZbar[46];
   armZZbar[103]=MMt*armZZbar[103];
   armZZbar[105]= - 3*armZZbar[46];
   armZZbar[117]=1 + armZZbar[150];
   armZZbar[117]=armZZbar[39]*armZZbar[117];
   armZZbar[122]= - 3*armZZbar[10];
   armZZbar[124]= - armZZbar[1] - 3*armZZbar[39];
   armZZbar[154]=armZZbar[6]*armZZbar[124];
   armZZbar[117]=4*armZZbar[154] + armZZbar[122] + 2*armZZbar[117] + 2*
   armZZbar[320] + armZZbar[316] + 132*armZZbar[13] - 165./2.*
   armZZbar[12] + armZZbar[379] + armZZbar[294] + armZZbar[246] + 379./
   6. + armZZbar[105];
   armZZbar[117]=armZZbar[37]*armZZbar[117];
   armZZbar[154]= - armZZbar[22] + armZZbar[193] + armZZbar[42] + 
   armZZbar[232];
   armZZbar[158]=13 + armZZbar[326];
   armZZbar[170]=armZZbar[158] + armZZbar[275];
   armZZbar[170]=armZZbar[19]*armZZbar[170];
   armZZbar[175]=armZZbar[175] + armZZbar[119];
   armZZbar[175]=armZZbar[35]*armZZbar[175];
   armZZbar[178]=5*armZZbar[35];
   armZZbar[158]=1./2.*armZZbar[158] + armZZbar[178];
   armZZbar[158]=armZZbar[18]*armZZbar[158];
   armZZbar[117]=armZZbar[117] + 1./2.*armZZbar[158] + 1./2.*
   armZZbar[175] - 1./4.*armZZbar[17] + 1./2.*armZZbar[170] - 3./2.*
   armZZbar[21] + 3*armZZbar[154] + armZZbar[324];
   armZZbar[117]=armZZbar[3]*armZZbar[117];
   armZZbar[154]=5./2. + armZZbar[283];
   armZZbar[154]=armZZbar[17]*armZZbar[154];
   armZZbar[110]=1./2.*armZZbar[154] - 7./2.*armZZbar[21] + 7./4.*
   armZZbar[20] + armZZbar[334] + armZZbar[232] + armZZbar[110] + 
   armZZbar[380] - armZZbar[68];
   armZZbar[110]=armZZbar[101]*armZZbar[110];
   armZZbar[154]=armZZbar[283] + armZZbar[354] - armZZbar[67] - 4 + 
   armZZbar[308];
   armZZbar[154]=armZZbar[101]*armZZbar[154];
   armZZbar[158]=2*armZZbar[360] + armZZbar[364];
   armZZbar[158]=armZZbar[3]*armZZbar[158];
   armZZbar[170]= - armZZbar[48] + 9 - armZZbar[46];
   armZZbar[144]=3*armZZbar[158] + armZZbar[144] - armZZbar[36] + 1./2.
   *armZZbar[170] + armZZbar[38];
   armZZbar[144]=armZZbar[3]*armZZbar[144];
   armZZbar[158]= - 3*armZZbar[57] + armZZbar[54] - 3*armZZbar[51] + 5./
   4.*armZZbar[55] + armZZbar[159];
   armZZbar[144]=3*armZZbar[144] + armZZbar[358] + armZZbar[154] + 1./2.
   *armZZbar[158] + armZZbar[223];
   armZZbar[144]=MMt*armZZbar[144];
   armZZbar[105]=armZZbar[160] + armZZbar[182] - 1799./27. + 
   armZZbar[105];
   armZZbar[105]=1./2.*armZZbar[105] - 55*armZZbar[12];
   armZZbar[105]=armZZbar[35] + 1./2.*armZZbar[239] + armZZbar[11] + 1./
   4.*armZZbar[105] - 19*armZZbar[13];
   armZZbar[105]=armZZbar[35]*armZZbar[105];
   armZZbar[154]=1./2. + 7*armZZbar[13];
   armZZbar[158]=14 + armZZbar[178];
   armZZbar[159]=armZZbar[1]*armZZbar[158];
   armZZbar[158]=armZZbar[39]*armZZbar[158];
   armZZbar[154]=1./3.*armZZbar[158] + 1./9.*armZZbar[159] + 1./4.*
   armZZbar[154] - armZZbar[35];
   armZZbar[154]=armZZbar[6]*armZZbar[154];
   armZZbar[158]= - 13./9. + armZZbar[7];
   armZZbar[158]=armZZbar[35]*armZZbar[158];
   armZZbar[158]=armZZbar[183] + 1./2.*armZZbar[158];
   armZZbar[159]=armZZbar[1]*armZZbar[158];
   armZZbar[158]=armZZbar[39]*armZZbar[158];
   armZZbar[170]=35./6. + armZZbar[283];
   armZZbar[170]=1./2.*armZZbar[170] + armZZbar[194];
   armZZbar[170]=armZZbar[10]*armZZbar[170];
   armZZbar[169]=armZZbar[367] + armZZbar[169];
   armZZbar[169]=armZZbar[18]*armZZbar[169];
   armZZbar[175]=55./4.*armZZbar[12];
   armZZbar[178]= - 3 + armZZbar[210];
   armZZbar[178]=armZZbar[36]*armZZbar[178];
   armZZbar[105]=armZZbar[144] + armZZbar[274] + armZZbar[117] + 4*
   armZZbar[288] + armZZbar[154] + armZZbar[169] + armZZbar[170] + 
   armZZbar[158] + 1./3.*armZZbar[159] + 1./2.*armZZbar[105] + 
   armZZbar[110] + 1./2.*armZZbar[178] - armZZbar[11] - 203./4.*
   armZZbar[13] + armZZbar[175] + armZZbar[390] + 3./4.*armZZbar[15] - 
   3./16.*armZZbar[44] + 33./64.*armZZbar[48] + 9./16.*armZZbar[46] + 7.
   /12.*armZZbar[67] + armZZbar[244] + 1./4.*armZZbar[59] - 12*
   armZZbar[16] + 3./64.*armZZbar[49] + 3533./128.*armZZbar[66] - 7./16.
   *armZZbar[61] - 25./32.*armZZbar[65] - 1./4.*armZZbar[62] - 
   armZZbar[72] + 41./64.*armZZbar[63] - 17197./864. + 12*armZZbar[60];
   armZZbar[105]=MMt*armZZbar[105];
   armZZbar[110]=33*armZZbar[12];
   armZZbar[117]= - 33*armZZbar[13];
   armZZbar[144]=armZZbar[117] - 1 + armZZbar[110];
   armZZbar[144]=1./2.*armZZbar[144] + armZZbar[296];
   armZZbar[154]= - armZZbar[1]*armZZbar[7];
   armZZbar[158]= - armZZbar[39]*armZZbar[7];
   armZZbar[159]=armZZbar[1] + 3*armZZbar[39];
   armZZbar[169]=armZZbar[6]*armZZbar[159];
   armZZbar[144]=2*armZZbar[169] + 6*armZZbar[10] + 6*armZZbar[158] + 3
   *armZZbar[144] + 2*armZZbar[154];
   armZZbar[144]=armZZbar[37]*armZZbar[144];
   armZZbar[145]=armZZbar[18]*armZZbar[145];
   armZZbar[145]=1./2.*armZZbar[145] + armZZbar[202] + armZZbar[353];
   armZZbar[144]=3./2.*armZZbar[145] + armZZbar[144];
   armZZbar[144]=armZZbar[3]*armZZbar[144];
   armZZbar[145]=33./2.*armZZbar[12];
   armZZbar[169]=armZZbar[117] - 59./9. + armZZbar[145];
   armZZbar[169]=1./2.*armZZbar[169] - armZZbar[11];
   armZZbar[169]=armZZbar[35]*armZZbar[169];
   armZZbar[170]=armZZbar[35]*armZZbar[263];
   armZZbar[170]=1./2.*armZZbar[170] + armZZbar[7] + armZZbar[279];
   armZZbar[178]=armZZbar[1]*armZZbar[170];
   armZZbar[170]=armZZbar[39]*armZZbar[170];
   armZZbar[182]=armZZbar[1]*armZZbar[276];
   armZZbar[193]=armZZbar[39]*armZZbar[276];
   armZZbar[182]=1./3.*armZZbar[182] + armZZbar[193];
   armZZbar[182]=armZZbar[6]*armZZbar[182];
   armZZbar[186]=33./4.*armZZbar[186] + armZZbar[11];
   armZZbar[193]= - MMt*armZZbar[3]*armZZbar[36];
   armZZbar[144]=3./2.*armZZbar[193] + armZZbar[144] + armZZbar[182] + 
   armZZbar[280] + armZZbar[170] + 1./3.*armZZbar[178] + 1./2.*
   armZZbar[169] + armZZbar[186] + 1./4.*armZZbar[234];
   armZZbar[144]=MMt*armZZbar[144];
   armZZbar[169]= - 17./18. + armZZbar[196];
   armZZbar[169]=armZZbar[19]*armZZbar[169];
   armZZbar[170]=47./6.*armZZbar[19] - armZZbar[17];
   armZZbar[178]=armZZbar[35]*armZZbar[170];
   armZZbar[182]= - 19./4.*armZZbar[35] + 17./6. + armZZbar[36];
   armZZbar[182]=armZZbar[18]*armZZbar[182];
   armZZbar[193]=armZZbar[170] - 19./6.*armZZbar[18];
   armZZbar[193]=armZZbar[6]*armZZbar[193];
   armZZbar[169]=1./2.*armZZbar[193] + 1./3.*armZZbar[182] + 1./2.*
   armZZbar[178] + armZZbar[169] + armZZbar[264];
   armZZbar[178]=1./2. - armZZbar[11];
   armZZbar[182]=armZZbar[35]*armZZbar[178];
   armZZbar[193]=armZZbar[10]*armZZbar[249];
   armZZbar[196]=1./2.*armZZbar[178] + armZZbar[10];
   armZZbar[196]=armZZbar[6]*armZZbar[196];
   armZZbar[181]=armZZbar[196] + armZZbar[193] + 1./2.*armZZbar[182] + 
   armZZbar[216] + armZZbar[181];
   armZZbar[181]=MMH*armZZbar[181];
   armZZbar[182]=armZZbar[1]*armZZbar[7];
   armZZbar[193]=armZZbar[39]*armZZbar[7];
   armZZbar[140]=armZZbar[140] - armZZbar[10] + armZZbar[193] + 
   armZZbar[186] + 1./3.*armZZbar[182];
   armZZbar[140]=armZZbar[37]*armZZbar[140];
   armZZbar[140]=armZZbar[144] + 1./3.*armZZbar[181] + 1./2.*
   armZZbar[169] + armZZbar[140];
   armZZbar[140]=armZZbar[133]*armZZbar[140];
   armZZbar[144]=armZZbar[160] + 15./4.*armZZbar[48] - 53359./54. + 
   armZZbar[129];
   armZZbar[144]=1./4.*armZZbar[144] + armZZbar[180];
   armZZbar[160]=3./2.*armZZbar[98] - armZZbar[34];
   armZZbar[160]=armZZbar[19]*armZZbar[160];
   armZZbar[169]=armZZbar[1]*armZZbar[183];
   armZZbar[180]=armZZbar[39]*armZZbar[183];
   armZZbar[181]=3*armZZbar[34] - 11./2.*armZZbar[99];
   armZZbar[181]=armZZbar[18]*armZZbar[181];
   armZZbar[183]=14./3.*armZZbar[39] + 5./4. + 14./9.*armZZbar[1];
   armZZbar[183]=armZZbar[6]*armZZbar[183];
   armZZbar[196]=3./4.*armZZbar[34] - 2*armZZbar[99];
   armZZbar[196]=armZZbar[37]*armZZbar[196];
   armZZbar[144]=armZZbar[196] + armZZbar[183] + 1./4.*armZZbar[181] + 
   41./12.*armZZbar[10] + armZZbar[180] + 1./3.*armZZbar[169] + 3./16.*
   armZZbar[35] + 1./2.*armZZbar[160] - armZZbar[11] + 1./4.*
   armZZbar[144] - 26*armZZbar[13];
   armZZbar[144]=armZZbar[37]*armZZbar[144];
   armZZbar[142]=armZZbar[142] - 13./2.*armZZbar[50] + 3*armZZbar[312]
    + armZZbar[190];
   armZZbar[142]=armZZbar[99]*armZZbar[142];
   armZZbar[160]=1./6.*armZZbar[14];
   armZZbar[142]=1./2.*armZZbar[142] - 3./2.*armZZbar[15] + 3./8.*
   armZZbar[44] - 5./4.*armZZbar[67] - 35./12.*armZZbar[64] + 
   armZZbar[157] + 21./8.*armZZbar[49] + armZZbar[160] - 21./8.*
   armZZbar[61] + armZZbar[28] - 109./72. - armZZbar[30];
   armZZbar[142]=3./8.*armZZbar[277] + 1./4.*armZZbar[142] + 
   armZZbar[188];
   armZZbar[155]= - 1./2.*armZZbar[15] + armZZbar[167] + 1./2.*
   armZZbar[238] + armZZbar[155];
   armZZbar[155]=armZZbar[99]*armZZbar[155];
   armZZbar[155]=armZZbar[155] + 1./2.*armZZbar[225];
   armZZbar[155]=MMH*armZZbar[155];
   armZZbar[157]=9./2. - armZZbar[44];
   armZZbar[157]=armZZbar[99]*armZZbar[157];
   armZZbar[157]=3./4.*armZZbar[157] + armZZbar[225];
   armZZbar[157]=armZZbar[37]*armZZbar[157];
   armZZbar[167]= - 2./9. + armZZbar[214];
   armZZbar[167]=armZZbar[35]*armZZbar[167];
   armZZbar[169]=armZZbar[207] - 41./24.*armZZbar[35];
   armZZbar[169]=armZZbar[10]*armZZbar[169];
   armZZbar[180]= - 1 - 25./8.*armZZbar[10];
   armZZbar[180]=armZZbar[6]*armZZbar[180];
   armZZbar[142]=1./32.*armZZbar[155] + 1./4.*armZZbar[157] + 1./9.*
   armZZbar[180] + 1./64.*armZZbar[351] + 1./3.*armZZbar[169] + 1./2.*
   armZZbar[142] + armZZbar[167];
   armZZbar[142]=MMH*armZZbar[142];
   armZZbar[155]= - 109*armZZbar[19] + armZZbar[219];
   armZZbar[155]=1./2.*armZZbar[155] + armZZbar[200];
   armZZbar[155]=armZZbar[6]*armZZbar[155];
   armZZbar[157]=943*armZZbar[70] - 399./2.*armZZbar[42];
   armZZbar[157]=1011./8.*armZZbar[50] + armZZbar[41] - 3./4.*
   armZZbar[40] + 11./6.*armZZbar[68] - 11*armZZbar[8] + 1./8.*
   armZZbar[157] + armZZbar[69];
   armZZbar[167]= - 607./32. + armZZbar[189];
   armZZbar[167]=armZZbar[19]*armZZbar[167];
   armZZbar[169]=1./2. - armZZbar[36];
   armZZbar[169]=armZZbar[17]*armZZbar[169];
   armZZbar[180]= - 101./4.*armZZbar[19] + armZZbar[17];
   armZZbar[180]=armZZbar[35]*armZZbar[180];
   armZZbar[181]=517./48.*armZZbar[35] + 307./16. - armZZbar[36];
   armZZbar[181]=armZZbar[18]*armZZbar[181];
   armZZbar[183]=armZZbar[132] + armZZbar[18];
   armZZbar[183]=7*armZZbar[183] - 3./2.*armZZbar[37];
   armZZbar[183]=armZZbar[3]*armZZbar[37]*armZZbar[183];
   armZZbar[105]=armZZbar[140] + armZZbar[105] + armZZbar[142] + 
   armZZbar[183] + armZZbar[144] + 1./12.*armZZbar[155] + 1./6.*
   armZZbar[181] + 1./6.*armZZbar[180] + 1./2.*armZZbar[169] + 1./2.*
   armZZbar[167] - 11./8.*armZZbar[21] - 5./8.*armZZbar[20] + 1./8.*
   armZZbar[157] - 10*armZZbar[22];
   armZZbar[105]=armZZbar[133]*armZZbar[105];
   armZZbar[129]=armZZbar[221] - 31./4.*armZZbar[48] + 140957./162. + 
   armZZbar[129];
   armZZbar[140]= - 1./2.*armZZbar[98] + armZZbar[253];
   armZZbar[140]=armZZbar[19]*armZZbar[140];
   armZZbar[142]=armZZbar[1]*armZZbar[164];
   armZZbar[144]=armZZbar[39]*armZZbar[192];
   armZZbar[155]= - armZZbar[34] + 35./6.*armZZbar[99];
   armZZbar[155]=armZZbar[18]*armZZbar[155];
   armZZbar[157]= - 59./3.*armZZbar[39] - 1./2. - 25./3.*armZZbar[1];
   armZZbar[157]=armZZbar[6]*armZZbar[157];
   armZZbar[164]= - 5./2.*armZZbar[34] + 44./3.*armZZbar[99];
   armZZbar[164]=armZZbar[37]*armZZbar[164];
   armZZbar[129]=armZZbar[164] + 1./3.*armZZbar[157] + 5./2.*
   armZZbar[155] + armZZbar[226] + 1./9.*armZZbar[144] + 1./27.*
   armZZbar[142] - 133./72.*armZZbar[35] + 1./2.*armZZbar[140] + 
   armZZbar[242] + 665./12.*armZZbar[13] + 1./8.*armZZbar[129] + 
   armZZbar[241];
   armZZbar[129]=armZZbar[37]*armZZbar[129];
   armZZbar[140]= - 19*armZZbar[19] + armZZbar[168];
   armZZbar[140]=2./3.*armZZbar[140] + armZZbar[203];
   armZZbar[140]=armZZbar[3]*armZZbar[37]*armZZbar[140];
   armZZbar[142]=69551./432. + 35*armZZbar[36];
   armZZbar[142]=armZZbar[19]*armZZbar[142];
   armZZbar[144]=367./3.*armZZbar[19] + armZZbar[228];
   armZZbar[144]=armZZbar[35]*armZZbar[144];
   armZZbar[155]= - 253./72. - armZZbar[36];
   armZZbar[155]=1./3.*armZZbar[155] + 41./8.*armZZbar[35];
   armZZbar[155]=armZZbar[18]*armZZbar[155];
   armZZbar[157]=49./2.*armZZbar[18] + 269./3.*armZZbar[19] - 
   armZZbar[17];
   armZZbar[157]=armZZbar[6]*armZZbar[157];
   armZZbar[103]=armZZbar[105] + armZZbar[103] + armZZbar[143] + 
   armZZbar[140] + armZZbar[129] + 1./18.*armZZbar[157] + 1./2.*
   armZZbar[155] + 1./9.*armZZbar[144] + armZZbar[314] + 1./6.*
   armZZbar[142] + 17./12.*armZZbar[21] + armZZbar[324] + 26./3.*
   armZZbar[22] - 2107./288.*armZZbar[50] + armZZbar[284] + 
   armZZbar[281] - 47./72.*armZZbar[68] + 3./4.*armZZbar[8] - 13./12.*
   armZZbar[69] + 1949./192.*armZZbar[42] + armZZbar[176] - 991./32.*
   armZZbar[70];
   armZZbar[103]=armZZbar[133]*armZZbar[103];
   armZZbar[105]= - 43./3. + armZZbar[162];
   armZZbar[129]= - 4*armZZbar[108];
   armZZbar[140]= - 13./3. - 2*armZZbar[108];
   armZZbar[140]=armZZbar[13]*armZZbar[140];
   armZZbar[105]= - 2./3.*armZZbar[35] + 2*armZZbar[140] + 1./3.*
   armZZbar[105] + armZZbar[129];
   armZZbar[105]=armZZbar[35]*armZZbar[105];
   armZZbar[142]=armZZbar[1]*armZZbar[248];
   armZZbar[143]=armZZbar[39]*armZZbar[248];
   armZZbar[142]=armZZbar[142] + 5./3.*armZZbar[143];
   armZZbar[142]=armZZbar[6]*armZZbar[142];
   armZZbar[143]=151 + 512./3.*armZZbar[62];
   armZZbar[144]= - 1 - armZZbar[35];
   armZZbar[155]=armZZbar[1]*armZZbar[144];
   armZZbar[144]=armZZbar[39]*armZZbar[144];
   armZZbar[121]=armZZbar[121] - armZZbar[55];
   armZZbar[121]=MMt*armZZbar[121];
   armZZbar[105]=512./9.*armZZbar[121] + 256./3.*armZZbar[142] + 1280./
   27.*armZZbar[144] + 256./9.*armZZbar[155] + 64*armZZbar[105] + 128*
   armZZbar[140] - 256*armZZbar[108] - 256./3.*armZZbar[46] - 256./9.*
   armZZbar[49] + 58*armZZbar[66] + 1./3.*armZZbar[143] - 29*
   armZZbar[65];
   armZZbar[105]=MMt*armZZbar[105];
   armZZbar[121]= - 167 - 64*armZZbar[46];
   armZZbar[142]=armZZbar[34] - 20./3.*armZZbar[99];
   armZZbar[142]=armZZbar[18]*armZZbar[142];
   armZZbar[143]=armZZbar[1] + 5./3.*armZZbar[39];
   armZZbar[143]=armZZbar[6]*armZZbar[143];
   armZZbar[144]=armZZbar[34] - 16./3.*armZZbar[99];
   armZZbar[144]=armZZbar[37]*armZZbar[144];
   armZZbar[121]=32*armZZbar[144] + 64./3.*armZZbar[143] + 32*
   armZZbar[142] - 320./27.*armZZbar[39] - 64./9.*armZZbar[1] - 64./3.*
   armZZbar[35] + 32*armZZbar[140] + 1./3.*armZZbar[121] - 64*
   armZZbar[108];
   armZZbar[121]=armZZbar[37]*armZZbar[121];
   armZZbar[140]=32./3.*armZZbar[71] + 29*armZZbar[70];
   armZZbar[115]= - 83./3. + armZZbar[115];
   armZZbar[115]=armZZbar[18]*armZZbar[115];
   armZZbar[142]= - 4*armZZbar[19] - 5./9.*armZZbar[18];
   armZZbar[142]=armZZbar[6]*armZZbar[142];
   armZZbar[115]=2*armZZbar[121] + 4*armZZbar[142] + 4./3.*
   armZZbar[115] + 64*armZZbar[166] - 94./3.*armZZbar[19] + 172./9.*
   armZZbar[21] + 128./9.*armZZbar[50] - 20./9.*armZZbar[8] - 64./3.*
   armZZbar[69] + 2*armZZbar[140] - 29*armZZbar[42];
   armZZbar[105]=2*armZZbar[115] + armZZbar[105];
   armZZbar[115]= - 113./2.*armZZbar[36];
   armZZbar[121]=1132./9. + armZZbar[115];
   armZZbar[140]= - 1501./54.*armZZbar[35] - 1691./81. + armZZbar[208];
   armZZbar[140]=armZZbar[35]*armZZbar[140];
   armZZbar[123]= - 259./9.*armZZbar[35] - 223./9. + armZZbar[123];
   armZZbar[123]=armZZbar[6]*armZZbar[123];
   armZZbar[121]=1./6.*armZZbar[123] + 1./9.*armZZbar[121] + 
   armZZbar[140];
   armZZbar[123]= - 20./9. + armZZbar[327];
   armZZbar[123]= - 1./3.*armZZbar[6] + 2*armZZbar[123] + 23./3.*
   armZZbar[35];
   armZZbar[123]=armZZbar[3]*armZZbar[37]*armZZbar[123];
   armZZbar[121]=1./3.*armZZbar[121] + armZZbar[123];
   armZZbar[121]=MMt*armZZbar[121];
   armZZbar[140]=armZZbar[6]*armZZbar[369];
   armZZbar[142]=119./9.*armZZbar[35] + 1121./27. + armZZbar[261];
   armZZbar[142]=armZZbar[35]*armZZbar[142];
   armZZbar[140]=5./2.*armZZbar[140] + 11./3.*armZZbar[387] + 1./2.*
   armZZbar[142];
   armZZbar[142]= - 5*armZZbar[6] + 22./3. + armZZbar[199];
   armZZbar[142]=armZZbar[3]*armZZbar[37]*armZZbar[142];
   armZZbar[140]=1./3.*armZZbar[140] + armZZbar[142];
   armZZbar[140]=MMt*armZZbar[140];
   armZZbar[140]=25./27.*armZZbar[148] + armZZbar[140];
   armZZbar[140]=armZZbar[272]*armZZbar[140];
   armZZbar[142]= - 223./18.*armZZbar[6] - 1735./18.*armZZbar[35] + 
   1132./27. + armZZbar[220];
   armZZbar[142]=armZZbar[37]*armZZbar[142];
   armZZbar[121]=1./3.*armZZbar[140] + 1./9.*armZZbar[142] + 
   armZZbar[121];
   armZZbar[121]=armZZbar[272]*armZZbar[121];
   armZZbar[140]=35./9. + armZZbar[36];
   armZZbar[140]=1./2.*armZZbar[140] + armZZbar[184];
   armZZbar[140]=armZZbar[35]*armZZbar[140];
   armZZbar[142]= - 2*armZZbar[6] + armZZbar[265] + 1 + armZZbar[327];
   armZZbar[142]=armZZbar[3]*armZZbar[37]*armZZbar[142];
   armZZbar[143]=armZZbar[6]*armZZbar[266];
   armZZbar[140]=2*armZZbar[142] + armZZbar[143] + armZZbar[149] + 
   armZZbar[140];
   armZZbar[140]=MMt*armZZbar[140];
   armZZbar[134]=armZZbar[134] - 2./3. - 3./4.*armZZbar[36];
   armZZbar[134]=armZZbar[35]*armZZbar[134];
   armZZbar[142]=armZZbar[6]*armZZbar[276];
   armZZbar[143]=armZZbar[6] + armZZbar[212] + armZZbar[35];
   armZZbar[143]=armZZbar[3]*armZZbar[37]*armZZbar[143];
   armZZbar[134]=3*armZZbar[143] + 1./2.*armZZbar[142] + 7./6.*
   armZZbar[36] + armZZbar[134];
   armZZbar[134]=MMt*armZZbar[134];
   armZZbar[135]=armZZbar[135] + armZZbar[156];
   armZZbar[135]=armZZbar[37]*armZZbar[135];
   armZZbar[134]=armZZbar[135] + armZZbar[134];
   armZZbar[134]=armZZbar[133]*armZZbar[134];
   armZZbar[135]=2*armZZbar[6] + armZZbar[149] + armZZbar[194];
   armZZbar[135]=armZZbar[37]*armZZbar[135];
   armZZbar[134]=armZZbar[134] + armZZbar[135] + armZZbar[140];
   armZZbar[134]=armZZbar[133]*armZZbar[134];
   armZZbar[115]=76 + armZZbar[115];
   armZZbar[135]= - 53./6.*armZZbar[35] - 67./9. + armZZbar[208];
   armZZbar[135]=armZZbar[35]*armZZbar[135];
   armZZbar[140]= - 43./3.*armZZbar[35] - 13 + 5./6.*armZZbar[36];
   armZZbar[140]=armZZbar[6]*armZZbar[140];
   armZZbar[115]=1./2.*armZZbar[140] + 1./9.*armZZbar[115] + 
   armZZbar[135];
   armZZbar[115]=1./3.*armZZbar[115] + armZZbar[123];
   armZZbar[115]=MMt*armZZbar[115];
   armZZbar[123]= - 79./2.*armZZbar[35] + 76./3. + armZZbar[220];
   armZZbar[123]=1./3.*armZZbar[123] - 13./2.*armZZbar[6];
   armZZbar[123]=armZZbar[37]*armZZbar[123];
   armZZbar[115]=armZZbar[134] + 1./3.*armZZbar[123] + armZZbar[115];
   armZZbar[115]=armZZbar[133]*armZZbar[115];
   armZZbar[123]=7./3. + armZZbar[231];
   armZZbar[123]=armZZbar[35]*armZZbar[123];
   armZZbar[134]=armZZbar[6]*armZZbar[248];
   armZZbar[123]=armZZbar[134] - 5./3. + armZZbar[123];
   armZZbar[123]=MMt*armZZbar[123];
   armZZbar[134]=armZZbar[235] + armZZbar[6];
   armZZbar[134]=armZZbar[37]*armZZbar[134];
   armZZbar[123]=armZZbar[134] + armZZbar[123];
   armZZbar[115]=armZZbar[115] + 256./81.*armZZbar[123] + armZZbar[121]
   ;
   armZZbar[115]=armZZbar[33]*armZZbar[115];
   armZZbar[103]=armZZbar[115] + armZZbar[103] + 1./9.*armZZbar[105] + 
   armZZbar[112];
   armZZbar[103]=armZZbar[33]*armZZbar[103];
   armZZbar[105]=35./27. + armZZbar[333];
   armZZbar[112]=2./3.*armZZbar[11];
   armZZbar[115]= - 7./6.*armZZbar[10];
   armZZbar[105]=armZZbar[115] + armZZbar[112] - 103./72.*armZZbar[13]
    + armZZbar[337] + armZZbar[306] + 1./4.*armZZbar[105] + 
   armZZbar[47];
   armZZbar[105]=armZZbar[10]*armZZbar[105];
   armZZbar[121]= - 1543./96. + armZZbar[88];
   armZZbar[123]= - 23./12. + armZZbar[13];
   armZZbar[123]=armZZbar[13]*armZZbar[123];
   armZZbar[134]=49./9. + armZZbar[13];
   armZZbar[134]=armZZbar[11]*armZZbar[134];
   armZZbar[135]= - 5./3.*armZZbar[87];
   armZZbar[140]= - 9./4.*armZZbar[86];
   armZZbar[142]= - 1./6.*armZZbar[96];
   armZZbar[143]= - 9./4.*armZZbar[43];
   armZZbar[144]= - 13./18.*armZZbar[12];
   armZZbar[105]=armZZbar[105] + 1./12.*armZZbar[134] + 1./6.*
   armZZbar[123] + armZZbar[144] + armZZbar[298] + 31./8.*armZZbar[15]
    - armZZbar[47] + armZZbar[143] + 13./12.*armZZbar[16] - 5./12.*
   armZZbar[89] + armZZbar[142] + armZZbar[365] + armZZbar[140] + 
   armZZbar[135] - 1./24.*armZZbar[90] - 17./12.*armZZbar[97] + 1./3.*
   armZZbar[121] - 1./8.*armZZbar[85];
   armZZbar[105]=armZZbar[73]*armZZbar[105];
   armZZbar[121]=armZZbar[10]*armZZbar[343];
   armZZbar[123]=1./2. + armZZbar[11];
   armZZbar[134]=armZZbar[1]*armZZbar[123];
   armZZbar[123]=armZZbar[39]*armZZbar[123];
   armZZbar[123]=armZZbar[134] + 11./9.*armZZbar[123];
   armZZbar[123]=armZZbar[6]*armZZbar[123];
   armZZbar[134]= - armZZbar[11] - 2./3. + armZZbar[304];
   armZZbar[134]=armZZbar[1]*armZZbar[134];
   armZZbar[148]= - 11./27.*armZZbar[11] - 10./27. + armZZbar[304];
   armZZbar[148]=armZZbar[39]*armZZbar[148];
   armZZbar[121]=armZZbar[123] + armZZbar[121] + 1./3.*armZZbar[134] + 
   armZZbar[148];
   armZZbar[121]=1./3.*armZZbar[121];
   armZZbar[123]= - 3*armZZbar[79];
   armZZbar[134]=armZZbar[123] - armZZbar[80];
   armZZbar[148]=armZZbar[134] - 7./12.*armZZbar[81];
   armZZbar[149]= - 1./3.*armZZbar[78];
   armZZbar[148]=1./2.*armZZbar[148] + armZZbar[149];
   armZZbar[148]=MMH*armZZbar[73]*armZZbar[148];
   armZZbar[105]=armZZbar[148] + armZZbar[121] + armZZbar[105];
   armZZbar[105]=MMH*armZZbar[105];
   armZZbar[148]=3*armZZbar[85];
   armZZbar[155]= - 3137./108. + armZZbar[148];
   armZZbar[155]= - 1./3.*armZZbar[90] + 1./2.*armZZbar[155] - 
   armZZbar[97];
   armZZbar[156]= - 1./2.*armZZbar[47];
   armZZbar[157]= - 11./3. + armZZbar[13];
   armZZbar[157]=armZZbar[11]*armZZbar[157];
   armZZbar[162]=17./36.*armZZbar[12];
   armZZbar[155]=1./18.*armZZbar[157] - 1./9.*armZZbar[13] + 
   armZZbar[162] + armZZbar[298] + 35./8.*armZZbar[15] + armZZbar[156]
    + armZZbar[143] + armZZbar[386] - 1./36.*armZZbar[89] - 7./6.*
   armZZbar[96] + armZZbar[365] + armZZbar[140] + 1./4.*armZZbar[155]
    + armZZbar[135];
   armZZbar[157]=3./4.*armZZbar[45];
   armZZbar[164]= - 11./12.*armZZbar[10] + armZZbar[216] + 13./144.*
   armZZbar[13] + armZZbar[162] + armZZbar[157] + armZZbar[291] + 52./
   27. + armZZbar[287];
   armZZbar[164]=armZZbar[10]*armZZbar[164];
   armZZbar[155]=1./2.*armZZbar[155] + armZZbar[164];
   armZZbar[155]=armZZbar[73]*armZZbar[155];
   armZZbar[164]=armZZbar[1]*armZZbar[218];
   armZZbar[166]=armZZbar[39]*armZZbar[218];
   armZZbar[167]=armZZbar[10]*armZZbar[341];
   armZZbar[164]=11./4.*armZZbar[167] + armZZbar[164] + 11./9.*
   armZZbar[166];
   armZZbar[164]=armZZbar[6]*armZZbar[164];
   armZZbar[166]=1./3.*armZZbar[80];
   armZZbar[168]=1./12.*armZZbar[81];
   armZZbar[169]=armZZbar[168] + armZZbar[123] + armZZbar[166];
   armZZbar[169]=1./2.*armZZbar[169] + armZZbar[149];
   armZZbar[169]=MMH*armZZbar[73]*armZZbar[169];
   armZZbar[160]= - armZZbar[15] + armZZbar[160] + armZZbar[28] - 7./9.
    - armZZbar[30];
   armZZbar[176]=1./4.*armZZbar[160] + armZZbar[138];
   armZZbar[180]=armZZbar[1]*armZZbar[176];
   armZZbar[176]=armZZbar[39]*armZZbar[176];
   armZZbar[155]=1./2.*armZZbar[169] + armZZbar[155] + 1./3.*
   armZZbar[164] + 1./9.*armZZbar[167] + armZZbar[180] + 11./9.*
   armZZbar[176];
   armZZbar[155]=MMH*armZZbar[155];
   armZZbar[164]= - 17./12.*armZZbar[12];
   armZZbar[167]=1./24.*armZZbar[13] + armZZbar[164] + armZZbar[300] + 
   armZZbar[236] + 77./9. + armZZbar[206];
   armZZbar[167]=armZZbar[17]*armZZbar[167];
   armZZbar[169]=armZZbar[132] - 53./12.*armZZbar[17];
   armZZbar[169]=armZZbar[101]*armZZbar[169];
   armZZbar[176]= - armZZbar[18]*armZZbar[101];
   armZZbar[180]=17./24.*armZZbar[176];
   armZZbar[181]=armZZbar[180] - 7./3.*armZZbar[10] + armZZbar[169] + 
   armZZbar[260] + 67./72.*armZZbar[13] - 17./72.*armZZbar[12] - 13./16.
   *armZZbar[45] - 7./8.*armZZbar[47] - 29./27. - 63./16.*armZZbar[43];
   armZZbar[181]=armZZbar[18]*armZZbar[181];
   armZZbar[183]= - 9*armZZbar[74];
   armZZbar[184]= - 11./3.*armZZbar[91];
   armZZbar[188]=armZZbar[184] + 15./2.*armZZbar[77] + armZZbar[183] + 
   1./18.*armZZbar[92];
   armZZbar[188]=1./2.*armZZbar[188] - armZZbar[76];
   armZZbar[189]= - 1./6.*armZZbar[75];
   armZZbar[188]=armZZbar[189] + 1./2.*armZZbar[188] - 7./3.*
   armZZbar[93];
   armZZbar[190]=armZZbar[18] - 37*armZZbar[19] + armZZbar[219];
   armZZbar[190]=armZZbar[18]*armZZbar[190];
   armZZbar[190]=armZZbar[286] + 1./4.*armZZbar[190];
   armZZbar[190]=armZZbar[3]*armZZbar[190];
   armZZbar[192]=armZZbar[209] - 2*armZZbar[286];
   armZZbar[192]=armZZbar[101]*armZZbar[192];
   armZZbar[194]= - 191./2.*armZZbar[13] + 1225./3. - 29*armZZbar[12];
   armZZbar[194]=armZZbar[19]*armZZbar[194];
   armZZbar[196]= - 43./4.*armZZbar[19] + armZZbar[17];
   armZZbar[196]=armZZbar[10]*armZZbar[196];
   armZZbar[167]=1./2.*armZZbar[190] + armZZbar[181] + 1./3.*
   armZZbar[196] + armZZbar[192] + 1./2.*armZZbar[167] + 1./72.*
   armZZbar[194] + 41./48.*armZZbar[21] + 221./48.*armZZbar[20] + 1./2.
   *armZZbar[188] + 2./3.*armZZbar[22];
   armZZbar[167]=armZZbar[73]*armZZbar[167];
   armZZbar[181]= - armZZbar[1]*armZZbar[19];
   armZZbar[188]= - armZZbar[39]*armZZbar[19];
   armZZbar[190]=armZZbar[18]*armZZbar[230];
   armZZbar[194]=armZZbar[190] + armZZbar[181] + 11./9.*armZZbar[188];
   armZZbar[196]=armZZbar[1]*armZZbar[19];
   armZZbar[199]=armZZbar[39]*armZZbar[19];
   armZZbar[196]=armZZbar[362] + armZZbar[196] + 11./9.*armZZbar[199];
   armZZbar[196]=armZZbar[6]*armZZbar[196];
   armZZbar[194]=1./3.*armZZbar[194] + armZZbar[196];
   armZZbar[196]= - 5./2.*armZZbar[13] + 71./3. - armZZbar[12];
   armZZbar[196]=armZZbar[19]*armZZbar[196];
   armZZbar[199]= - 1./2.*armZZbar[12];
   armZZbar[200]= - 1 + armZZbar[199];
   armZZbar[200]=armZZbar[17]*armZZbar[200];
   armZZbar[202]=5./4.*armZZbar[13] + 5./3. + armZZbar[199];
   armZZbar[202]=1./3.*armZZbar[202] - armZZbar[10];
   armZZbar[202]=armZZbar[18]*armZZbar[202];
   armZZbar[196]=armZZbar[202] + armZZbar[372] + armZZbar[200] + 1./6.*
   armZZbar[196] + 1./2.*armZZbar[21] + armZZbar[22] + 1./2.*
   armZZbar[77] - armZZbar[93];
   armZZbar[200]=armZZbar[3]*armZZbar[18]*armZZbar[357];
   armZZbar[196]=1./3.*armZZbar[196] + 1./2.*armZZbar[200];
   armZZbar[196]=armZZbar[73]*armZZbar[196];
   armZZbar[200]=armZZbar[10]*armZZbar[12];
   armZZbar[200]=1./3.*armZZbar[200] + 1./6.*armZZbar[12] - 
   armZZbar[11];
   armZZbar[200]=MMH*armZZbar[73]*armZZbar[200];
   armZZbar[194]=1./3.*armZZbar[200] + 1./3.*armZZbar[194] + 
   armZZbar[196];
   armZZbar[194]=armZZbar[272]*armZZbar[194];
   armZZbar[196]= - 43./3.*armZZbar[19] + armZZbar[17];
   armZZbar[196]=armZZbar[1]*armZZbar[196];
   armZZbar[200]= - 377./3.*armZZbar[19] + armZZbar[332];
   armZZbar[200]=armZZbar[39]*armZZbar[200];
   armZZbar[196]=armZZbar[196] + 1./9.*armZZbar[200];
   armZZbar[200]= - 23*armZZbar[1] - 359./27.*armZZbar[39];
   armZZbar[200]=armZZbar[18]*armZZbar[200];
   armZZbar[196]=1./4.*armZZbar[196] + 1./3.*armZZbar[200];
   armZZbar[196]=armZZbar[6]*armZZbar[196];
   armZZbar[200]= - 1./4.*armZZbar[20];
   armZZbar[202]=1./8.*armZZbar[17] + 67./36.*armZZbar[19] - 29./4.*
   armZZbar[21] + armZZbar[200] - 31./4.*armZZbar[8] - armZZbar[22];
   armZZbar[202]=armZZbar[1]*armZZbar[202];
   armZZbar[203]=11./24.*armZZbar[17] + 449./108.*armZZbar[19] - 1271./
   108.*armZZbar[21] - 11./12.*armZZbar[20] - 1469./108.*armZZbar[8] - 
   armZZbar[22];
   armZZbar[203]=armZZbar[39]*armZZbar[203];
   armZZbar[205]=713*armZZbar[1] + 3683./9.*armZZbar[39];
   armZZbar[205]=armZZbar[18]*armZZbar[205];
   armZZbar[155]=1./4.*armZZbar[194] + armZZbar[155] + armZZbar[167] + 
   armZZbar[196] + 1./36.*armZZbar[205] + armZZbar[202] + 1./3.*
   armZZbar[203];
   armZZbar[155]=armZZbar[272]*armZZbar[155];
   armZZbar[167]= - 63*armZZbar[43];
   armZZbar[194]= - 391./27. + armZZbar[167];
   armZZbar[196]= - 7*armZZbar[47];
   armZZbar[202]= - 13./2.*armZZbar[45];
   armZZbar[194]=armZZbar[202] + 1./2.*armZZbar[194] + armZZbar[196];
   armZZbar[194]=107./9.*armZZbar[13] + 1./2.*armZZbar[194] + 
   armZZbar[336];
   armZZbar[203]=4*armZZbar[19] - 53./6.*armZZbar[17];
   armZZbar[203]=armZZbar[101]*armZZbar[203];
   armZZbar[205]= - 5./3.*armZZbar[10];
   armZZbar[176]=17./12.*armZZbar[176];
   armZZbar[207]= - armZZbar[19]*armZZbar[100];
   armZZbar[208]=3./2.*armZZbar[207];
   armZZbar[194]=armZZbar[176] + armZZbar[205] + armZZbar[203] + 
   armZZbar[208] + 1./2.*armZZbar[194] + armZZbar[216];
   armZZbar[194]=armZZbar[18]*armZZbar[194];
   armZZbar[212]= - 27./2.*armZZbar[43];
   armZZbar[213]= - 3*armZZbar[47];
   armZZbar[214]=armZZbar[335] + armZZbar[213] + 557./3. + 
   armZZbar[212];
   armZZbar[218]= - 3*armZZbar[12];
   armZZbar[214]=11./12.*armZZbar[207] - 7./6.*armZZbar[11] - 503./72.*
   armZZbar[13] + 1./4.*armZZbar[214] + armZZbar[218];
   armZZbar[214]=armZZbar[19]*armZZbar[214];
   armZZbar[220]=431./9. + armZZbar[333];
   armZZbar[221]=13./3.*armZZbar[12];
   armZZbar[220]=armZZbar[221] + 1./2.*armZZbar[220] + armZZbar[45];
   armZZbar[222]=armZZbar[19]*armZZbar[100];
   armZZbar[220]=1./6.*armZZbar[222] + armZZbar[118] + 1./2.*
   armZZbar[220] - 1./3.*armZZbar[13];
   armZZbar[220]=armZZbar[17]*armZZbar[220];
   armZZbar[223]= - 3./4.*armZZbar[19] + armZZbar[17];
   armZZbar[223]=armZZbar[17]*armZZbar[223];
   armZZbar[136]=5*armZZbar[136] + 17*armZZbar[18];
   armZZbar[136]=1./4.*armZZbar[18]*armZZbar[136];
   armZZbar[225]=armZZbar[136] - 1./2.*armZZbar[209] + armZZbar[223];
   armZZbar[225]=armZZbar[3]*armZZbar[225];
   armZZbar[226]=2*armZZbar[192];
   armZZbar[228]= - 9./4.*armZZbar[74];
   armZZbar[229]= - 11./12.*armZZbar[91];
   armZZbar[230]= - 7*armZZbar[19] + armZZbar[119];
   armZZbar[230]=1./3.*armZZbar[10]*armZZbar[230];
   armZZbar[194]=armZZbar[225] + armZZbar[194] + armZZbar[230] + 
   armZZbar[226] + armZZbar[220] + armZZbar[214] + 33./4.*armZZbar[21]
    + 511./72.*armZZbar[20] + 67./12.*armZZbar[22] + armZZbar[189] - 
   151./12.*armZZbar[93] + 2./3.*armZZbar[76] + armZZbar[229] + 67./12.
   *armZZbar[77] + 41./72.*armZZbar[92] + 10./9.*armZZbar[94] + 
   armZZbar[228];
   armZZbar[194]=armZZbar[73]*armZZbar[194];
   armZZbar[214]= - 2*armZZbar[9];
   armZZbar[220]=685./9.*armZZbar[21] + armZZbar[107] + armZZbar[214]
    + 685./9.*armZZbar[8];
   armZZbar[225]= - 1./3.*armZZbar[7];
   armZZbar[231]= - 1 + armZZbar[225];
   armZZbar[231]=armZZbar[19]*armZZbar[231];
   armZZbar[232]=10./27. + armZZbar[302];
   armZZbar[232]=armZZbar[17]*armZZbar[232];
   armZZbar[220]=armZZbar[232] + 1./9.*armZZbar[220] + 1./2.*
   armZZbar[231];
   armZZbar[220]=armZZbar[39]*armZZbar[220];
   armZZbar[231]=55./9.*armZZbar[19] - armZZbar[17];
   armZZbar[231]=armZZbar[1]*armZZbar[231];
   armZZbar[234]=armZZbar[39]*armZZbar[254];
   armZZbar[235]=293*armZZbar[1] + 1631./9.*armZZbar[39];
   armZZbar[235]=armZZbar[18]*armZZbar[235];
   armZZbar[231]=1./9.*armZZbar[235] + armZZbar[231] + 11./9.*
   armZZbar[234];
   armZZbar[231]=armZZbar[6]*armZZbar[231];
   armZZbar[225]=7 + armZZbar[225];
   armZZbar[225]=armZZbar[19]*armZZbar[225];
   armZZbar[234]=2./3. + armZZbar[302];
   armZZbar[235]=armZZbar[17]*armZZbar[234];
   armZZbar[225]=1./3.*armZZbar[235] + 1./6.*armZZbar[225] + 15*
   armZZbar[21] - 2./3.*armZZbar[22] - 2./3.*armZZbar[9] + 15*
   armZZbar[8];
   armZZbar[225]=armZZbar[1]*armZZbar[225];
   armZZbar[237]= - 1064./3. + armZZbar[302];
   armZZbar[237]=armZZbar[1]*armZZbar[237];
   armZZbar[238]= - 1820./27. + armZZbar[302];
   armZZbar[238]=armZZbar[39]*armZZbar[238];
   armZZbar[237]=1./3.*armZZbar[237] + armZZbar[238];
   armZZbar[237]=armZZbar[18]*armZZbar[237];
   armZZbar[105]=armZZbar[155] + armZZbar[105] + armZZbar[194] + 1./2.*
   armZZbar[231] + 1./3.*armZZbar[237] + armZZbar[225] + armZZbar[220];
   armZZbar[105]=armZZbar[272]*armZZbar[105];
   armZZbar[155]= - 469./27. + armZZbar[333];
   armZZbar[112]=armZZbar[115] + armZZbar[112] - 367./72.*armZZbar[13]
    + armZZbar[337] + armZZbar[306] + 1./4.*armZZbar[155] + 
   armZZbar[47];
   armZZbar[112]=armZZbar[10]*armZZbar[112];
   armZZbar[115]= - 959./96. + armZZbar[88];
   armZZbar[155]= - 167./36. - armZZbar[13];
   armZZbar[155]=armZZbar[13]*armZZbar[155];
   armZZbar[194]= - 167./9. - 23*armZZbar[13];
   armZZbar[194]=armZZbar[11]*armZZbar[194];
   armZZbar[112]=armZZbar[112] + 1./12.*armZZbar[194] + 1./2.*
   armZZbar[155] + armZZbar[144] + armZZbar[298] + 39./8.*armZZbar[15]
    - armZZbar[47] + armZZbar[143] - 1./4.*armZZbar[16] - 197./12.*
   armZZbar[89] + armZZbar[142] + armZZbar[365] + armZZbar[140] + 
   armZZbar[135] + 89./8.*armZZbar[90] - 1./12.*armZZbar[97] + 1./3.*
   armZZbar[115] - 9./8.*armZZbar[85];
   armZZbar[112]=armZZbar[73]*armZZbar[112];
   armZZbar[115]=armZZbar[134] + 17./12.*armZZbar[81];
   armZZbar[115]=1./2.*armZZbar[115] + armZZbar[149];
   armZZbar[115]=MMH*armZZbar[73]*armZZbar[115];
   armZZbar[112]=armZZbar[115] + armZZbar[121] + armZZbar[112];
   armZZbar[112]=MMH*armZZbar[112];
   armZZbar[115]= - 85759./27. + armZZbar[167];
   armZZbar[115]=armZZbar[202] + 1./2.*armZZbar[115] + armZZbar[196];
   armZZbar[115]= - 2335./9.*armZZbar[13] + 1./2.*armZZbar[115] + 
   armZZbar[336];
   armZZbar[115]=armZZbar[176] + armZZbar[205] + armZZbar[203] + 
   armZZbar[208] + 1./2.*armZZbar[115] + armZZbar[216];
   armZZbar[115]=armZZbar[18]*armZZbar[115];
   armZZbar[121]=armZZbar[335] + armZZbar[213] + 26863./9. + 
   armZZbar[212];
   armZZbar[121]=59./12.*armZZbar[207] + 29./6.*armZZbar[11] - 1423./8.
   *armZZbar[13] + 1./4.*armZZbar[121] - 247./3.*armZZbar[12];
   armZZbar[121]=armZZbar[19]*armZZbar[121];
   armZZbar[134]=529./3. + armZZbar[333];
   armZZbar[134]=armZZbar[221] + 1./2.*armZZbar[134] + armZZbar[45];
   armZZbar[134]=25./6.*armZZbar[222] + armZZbar[118] + 1./2.*
   armZZbar[134] + 4./3.*armZZbar[13];
   armZZbar[134]=armZZbar[17]*armZZbar[134];
   armZZbar[135]=armZZbar[136] + 47./2.*armZZbar[209] + armZZbar[223];
   armZZbar[135]=armZZbar[3]*armZZbar[135];
   armZZbar[115]=armZZbar[135] + armZZbar[115] + armZZbar[230] + 
   armZZbar[226] + armZZbar[134] + armZZbar[121] + 3467./12.*
   armZZbar[21] + 101./24.*armZZbar[20] + 7099./12.*armZZbar[22] + 
   armZZbar[189] - 2087./12.*armZZbar[93] + 40./3.*armZZbar[76] + 
   armZZbar[229] - 455./4.*armZZbar[77] - 605./24.*armZZbar[92] + 2*
   armZZbar[94] + armZZbar[228];
   armZZbar[115]=armZZbar[73]*armZZbar[115];
   armZZbar[118]=3*armZZbar[207] + armZZbar[118] - 95./16.*armZZbar[13]
    + 55./8.*armZZbar[12] + armZZbar[251] + armZZbar[236] - 164./9. + 
   armZZbar[287];
   armZZbar[118]=armZZbar[17]*armZZbar[118];
   armZZbar[121]=137 - 7./2.*armZZbar[43];
   armZZbar[121]=107*armZZbar[13] + 55./3.*armZZbar[12] + armZZbar[202]
    + 9*armZZbar[121] - 13*armZZbar[47];
   armZZbar[121]=armZZbar[222] + 1./4.*armZZbar[121] + armZZbar[151];
   armZZbar[121]=armZZbar[180] + armZZbar[269] + 1./2.*armZZbar[121] + 
   armZZbar[169];
   armZZbar[121]=armZZbar[18]*armZZbar[121];
   armZZbar[134]= - 14143./3. - 27*armZZbar[43];
   armZZbar[134]=6091./12.*armZZbar[13] + 253./2.*armZZbar[12] + 
   armZZbar[335] + 1./2.*armZZbar[134] + armZZbar[213];
   armZZbar[134]=13./2.*armZZbar[222] + 1./2.*armZZbar[134] - 13*
   armZZbar[11];
   armZZbar[134]=armZZbar[19]*armZZbar[134];
   armZZbar[135]= - 3./2.*armZZbar[19];
   armZZbar[136]=armZZbar[135] + armZZbar[17];
   armZZbar[136]=armZZbar[17]*armZZbar[136];
   armZZbar[142]=33*armZZbar[18] - 9*armZZbar[19] + armZZbar[219];
   armZZbar[142]=armZZbar[18]*armZZbar[142];
   armZZbar[136]=1./4.*armZZbar[142] - 33*armZZbar[209] + armZZbar[136]
   ;
   armZZbar[136]=armZZbar[3]*armZZbar[136];
   armZZbar[142]=armZZbar[184] + 489./2.*armZZbar[77] + armZZbar[183]
    + 197./2.*armZZbar[92];
   armZZbar[142]=1./4.*armZZbar[142] - 11*armZZbar[76];
   armZZbar[143]= - 53./2.*armZZbar[19] + armZZbar[311];
   armZZbar[143]=armZZbar[10]*armZZbar[143];
   armZZbar[118]=1./2.*armZZbar[136] + armZZbar[121] + 1./2.*
   armZZbar[143] + armZZbar[192] + armZZbar[118] + 1./2.*armZZbar[134]
    - 3023./16.*armZZbar[21] + 151./16.*armZZbar[20] - 1571./4.*
   armZZbar[22] + armZZbar[309] + 1./2.*armZZbar[142] + 159*
   armZZbar[93];
   armZZbar[118]=armZZbar[73]*armZZbar[118];
   armZZbar[121]=491./24.*armZZbar[13] - 55./6.*armZZbar[12] + 
   armZZbar[306] + armZZbar[305] + 503./27. + armZZbar[206];
   armZZbar[121]= - 1./4.*armZZbar[10] + 1./2.*armZZbar[121] + 
   armZZbar[216];
   armZZbar[121]=armZZbar[10]*armZZbar[121];
   armZZbar[134]= - 5./2.*armZZbar[87] + 29./8.*armZZbar[97] + 29./16.*
   armZZbar[85] + 1667./576. + 2*armZZbar[88];
   armZZbar[136]=115./4. + armZZbar[13];
   armZZbar[136]=armZZbar[13]*armZZbar[136];
   armZZbar[142]= - 1 + armZZbar[210];
   armZZbar[142]=armZZbar[11]*armZZbar[142];
   armZZbar[121]=armZZbar[121] + 1./3.*armZZbar[142] + 1./6.*
   armZZbar[136] - 55./24.*armZZbar[12] + armZZbar[233] + 137./48.*
   armZZbar[15] + armZZbar[191] + armZZbar[171] - 15./8.*armZZbar[16]
    + 95./8.*armZZbar[89] - 5./3.*armZZbar[96] + armZZbar[137] + 1./3.*
   armZZbar[134] - 9./8.*armZZbar[86];
   armZZbar[121]=armZZbar[73]*armZZbar[121];
   armZZbar[134]=1./2.*armZZbar[160] + armZZbar[187];
   armZZbar[136]=armZZbar[1]*armZZbar[134];
   armZZbar[134]=armZZbar[39]*armZZbar[134];
   armZZbar[134]=1./3.*armZZbar[136] + armZZbar[134];
   armZZbar[136]= - 2./3. + armZZbar[7];
   armZZbar[137]=armZZbar[1]*armZZbar[136];
   armZZbar[136]=armZZbar[39]*armZZbar[136];
   armZZbar[136]=1./3.*armZZbar[137] + armZZbar[136];
   armZZbar[136]=armZZbar[10]*armZZbar[136];
   armZZbar[137]=armZZbar[139] + 11./4.*armZZbar[322];
   armZZbar[137]=armZZbar[6]*armZZbar[137];
   armZZbar[123]= - 11./12.*armZZbar[81] + armZZbar[123] - 11./3.*
   armZZbar[80];
   armZZbar[123]=1./2.*armZZbar[123] + armZZbar[149];
   armZZbar[123]=MMH*armZZbar[73]*armZZbar[123];
   armZZbar[121]=1./2.*armZZbar[123] + armZZbar[121] + 1./3.*
   armZZbar[137] + 1./2.*armZZbar[134] + 1./3.*armZZbar[136];
   armZZbar[121]=MMH*armZZbar[121];
   armZZbar[123]=armZZbar[1]*armZZbar[170];
   armZZbar[134]=armZZbar[39]*armZZbar[170];
   armZZbar[123]=19./6.*armZZbar[290] + 1./3.*armZZbar[123] + 
   armZZbar[134];
   armZZbar[123]=armZZbar[6]*armZZbar[123];
   armZZbar[134]= - 17./18. + armZZbar[179];
   armZZbar[134]=armZZbar[19]*armZZbar[134];
   armZZbar[134]=armZZbar[134] + armZZbar[173];
   armZZbar[136]=armZZbar[1]*armZZbar[134];
   armZZbar[134]=armZZbar[39]*armZZbar[134];
   armZZbar[137]=17./6. + armZZbar[7];
   armZZbar[139]=armZZbar[1]*armZZbar[137];
   armZZbar[137]=armZZbar[39]*armZZbar[137];
   armZZbar[137]=1./3.*armZZbar[139] + armZZbar[137];
   armZZbar[137]=armZZbar[18]*armZZbar[137];
   armZZbar[123]=armZZbar[123] + 1./3.*armZZbar[137] + 1./3.*
   armZZbar[136] + armZZbar[134];
   armZZbar[134]= - 517./2.*armZZbar[13] - 1057./27. + 165*armZZbar[12]
   ;
   armZZbar[134]=armZZbar[208] + 1./4.*armZZbar[134] + armZZbar[313];
   armZZbar[134]=armZZbar[19]*armZZbar[134];
   armZZbar[136]=armZZbar[17]*armZZbar[186];
   armZZbar[137]=armZZbar[101]*armZZbar[209];
   armZZbar[139]=59./6.*armZZbar[19] - armZZbar[17];
   armZZbar[139]=armZZbar[10]*armZZbar[139];
   armZZbar[134]=armZZbar[139] + 3./2.*armZZbar[137] + armZZbar[134] + 
   armZZbar[136];
   armZZbar[136]=209./2.*armZZbar[13] + 1057./27. - 11*armZZbar[12];
   armZZbar[137]=armZZbar[101]*armZZbar[19];
   armZZbar[136]= - 13./12.*armZZbar[10] + 3./2.*armZZbar[137] + 
   armZZbar[208] + 1./8.*armZZbar[136] + armZZbar[151];
   armZZbar[136]=armZZbar[18]*armZZbar[136];
   armZZbar[137]= - armZZbar[19] - armZZbar[18];
   armZZbar[137]=armZZbar[18]*armZZbar[137];
   armZZbar[137]=armZZbar[209] + 1./2.*armZZbar[137];
   armZZbar[137]=armZZbar[3]*armZZbar[137];
   armZZbar[134]=17./4.*armZZbar[137] + 1./2.*armZZbar[134] + 
   armZZbar[136];
   armZZbar[134]=armZZbar[73]*armZZbar[134];
   armZZbar[136]=armZZbar[1]*armZZbar[178];
   armZZbar[137]=armZZbar[39]*armZZbar[178];
   armZZbar[136]=2*armZZbar[317] + 1./3.*armZZbar[136] + armZZbar[137];
   armZZbar[136]=armZZbar[6]*armZZbar[136];
   armZZbar[137]=armZZbar[302] + armZZbar[216];
   armZZbar[139]=armZZbar[1]*armZZbar[137];
   armZZbar[137]=armZZbar[39]*armZZbar[137];
   armZZbar[142]=armZZbar[10]*armZZbar[318];
   armZZbar[136]=armZZbar[136] + armZZbar[142] + 1./3.*armZZbar[139] + 
   armZZbar[137];
   armZZbar[137]=127./27. + armZZbar[323];
   armZZbar[137]=armZZbar[11]*armZZbar[137];
   armZZbar[137]=11./2.*armZZbar[147] + armZZbar[137];
   armZZbar[139]= - 127./27. + armZZbar[163];
   armZZbar[139]=1./2.*armZZbar[139] - 11*armZZbar[13];
   armZZbar[139]=armZZbar[269] + 1./2.*armZZbar[139] - 2./3.*
   armZZbar[11];
   armZZbar[139]=armZZbar[10]*armZZbar[139];
   armZZbar[137]=1./4.*armZZbar[137] + armZZbar[139];
   armZZbar[137]=armZZbar[73]*armZZbar[137];
   armZZbar[136]=1./3.*armZZbar[136] + armZZbar[137];
   armZZbar[136]=MMH*armZZbar[136];
   armZZbar[123]=armZZbar[136] + 1./2.*armZZbar[123] + armZZbar[134];
   armZZbar[123]=armZZbar[133]*armZZbar[123];
   armZZbar[111]= - armZZbar[9] + armZZbar[111];
   armZZbar[111]=armZZbar[200] + 11./2.*armZZbar[111] - 13*armZZbar[22]
   ;
   armZZbar[134]=13 - 11./3.*armZZbar[7];
   armZZbar[134]=armZZbar[19]*armZZbar[134];
   armZZbar[136]=1./4. - armZZbar[7];
   armZZbar[136]=armZZbar[17]*armZZbar[136];
   armZZbar[134]=1./6.*armZZbar[136] + 1./2.*armZZbar[134] + 1./3.*
   armZZbar[111] - 3./4.*armZZbar[21];
   armZZbar[134]=armZZbar[1]*armZZbar[134];
   armZZbar[137]= - 11*armZZbar[7];
   armZZbar[139]=39 + armZZbar[137];
   armZZbar[139]=armZZbar[19]*armZZbar[139];
   armZZbar[111]=1./2.*armZZbar[136] + 1./2.*armZZbar[139] + 
   armZZbar[111] - 9./4.*armZZbar[21];
   armZZbar[111]=armZZbar[39]*armZZbar[111];
   armZZbar[136]=71./2. - armZZbar[7];
   armZZbar[139]=armZZbar[1]*armZZbar[136];
   armZZbar[136]=armZZbar[39]*armZZbar[136];
   armZZbar[136]=1./3.*armZZbar[139] + armZZbar[136];
   armZZbar[136]=armZZbar[18]*armZZbar[136];
   armZZbar[139]= - 17*armZZbar[19] + armZZbar[17];
   armZZbar[142]=armZZbar[1]*armZZbar[139];
   armZZbar[139]=armZZbar[39]*armZZbar[139];
   armZZbar[139]=armZZbar[142] + 3*armZZbar[139];
   armZZbar[139]=1./2.*armZZbar[139] + 1./3.*armZZbar[197];
   armZZbar[139]=armZZbar[6]*armZZbar[139];
   armZZbar[111]=armZZbar[123] + armZZbar[121] + armZZbar[118] + 1./2.*
   armZZbar[139] + 1./6.*armZZbar[136] + armZZbar[134] + armZZbar[111];
   armZZbar[111]=armZZbar[133]*armZZbar[111];
   armZZbar[118]=35*armZZbar[7];
   armZZbar[121]= - 449./3. + armZZbar[118];
   armZZbar[121]=armZZbar[19]*armZZbar[121];
   armZZbar[123]=5*armZZbar[21];
   armZZbar[134]=armZZbar[123] + 8*armZZbar[22] + 14*armZZbar[9] + 5*
   armZZbar[8];
   armZZbar[121]=armZZbar[235] + armZZbar[134] + 1./6.*armZZbar[121];
   armZZbar[121]=armZZbar[1]*armZZbar[121];
   armZZbar[118]= - 3937./27. + armZZbar[118];
   armZZbar[118]=armZZbar[19]*armZZbar[118];
   armZZbar[118]=armZZbar[232] + armZZbar[134] + 1./6.*armZZbar[118];
   armZZbar[118]=armZZbar[39]*armZZbar[118];
   armZZbar[134]=37./3.*armZZbar[19] - armZZbar[17];
   armZZbar[134]=armZZbar[1]*armZZbar[134];
   armZZbar[136]=895./3.*armZZbar[19] + armZZbar[295];
   armZZbar[136]=armZZbar[39]*armZZbar[136];
   armZZbar[139]=53*armZZbar[1] + 119*armZZbar[39];
   armZZbar[139]=armZZbar[18]*armZZbar[139];
   armZZbar[134]=1./9.*armZZbar[139] + armZZbar[134] + 1./9.*
   armZZbar[136];
   armZZbar[134]=armZZbar[6]*armZZbar[134];
   armZZbar[136]= - 128./3. + armZZbar[302];
   armZZbar[136]=armZZbar[1]*armZZbar[136];
   armZZbar[139]= - 364./9. + armZZbar[302];
   armZZbar[139]=armZZbar[39]*armZZbar[139];
   armZZbar[136]=1./3.*armZZbar[136] + armZZbar[139];
   armZZbar[136]=armZZbar[18]*armZZbar[136];
   armZZbar[111]=armZZbar[111] + armZZbar[112] + armZZbar[115] + 1./2.*
   armZZbar[134] + 1./3.*armZZbar[136] + 1./3.*armZZbar[121] + 
   armZZbar[118];
   armZZbar[111]=armZZbar[133]*armZZbar[111];
   armZZbar[112]=4*armZZbar[108];
   armZZbar[115]=55./3. + armZZbar[112];
   armZZbar[115]=armZZbar[13]*armZZbar[115];
   armZZbar[115]=armZZbar[115] + 45 + 8*armZZbar[108];
   armZZbar[115]=armZZbar[18]*armZZbar[115];
   armZZbar[118]=23./3.*armZZbar[77] - 1./9.*armZZbar[94] + 
   armZZbar[92];
   armZZbar[121]=armZZbar[108]*armZZbar[77];
   armZZbar[134]= - 247./9. + armZZbar[129];
   armZZbar[134]=armZZbar[22]*armZZbar[134];
   armZZbar[136]= - 20./3. - armZZbar[108];
   armZZbar[136]=armZZbar[21]*armZZbar[136];
   armZZbar[139]=113./9. + armZZbar[112];
   armZZbar[139]=armZZbar[13]*armZZbar[139];
   armZZbar[139]=armZZbar[139] - 233./9. + 3*armZZbar[12];
   armZZbar[139]=armZZbar[19]*armZZbar[139];
   armZZbar[115]=armZZbar[115] + armZZbar[267] + 2*armZZbar[139] + 4*
   armZZbar[136] + 2*armZZbar[134] + 4*armZZbar[121] + 34./3.*
   armZZbar[93] + 2*armZZbar[118] - armZZbar[76];
   armZZbar[115]=armZZbar[73]*armZZbar[115];
   armZZbar[118]=armZZbar[181] + 5./3.*armZZbar[188];
   armZZbar[121]= - armZZbar[1] - 17./27.*armZZbar[39];
   armZZbar[121]=armZZbar[18]*armZZbar[121];
   armZZbar[118]=4*armZZbar[118] + 5*armZZbar[121];
   armZZbar[118]=armZZbar[6]*armZZbar[118];
   armZZbar[121]= - 4*armZZbar[21] - armZZbar[9] - 4*armZZbar[8];
   armZZbar[121]=5*armZZbar[121] + 31./3.*armZZbar[19];
   armZZbar[121]=armZZbar[1]*armZZbar[121];
   armZZbar[134]=179./9.*armZZbar[19] - 340./27.*armZZbar[21] - 11*
   armZZbar[9] - 340./27.*armZZbar[8];
   armZZbar[134]=armZZbar[39]*armZZbar[134];
   armZZbar[136]=armZZbar[1] + 17./27.*armZZbar[39];
   armZZbar[136]=armZZbar[18]*armZZbar[136];
   armZZbar[118]=4*armZZbar[118] + 52*armZZbar[136] + armZZbar[121] + 
   armZZbar[134];
   armZZbar[115]=1./3.*armZZbar[118] + 2*armZZbar[115];
   armZZbar[118]= - 5 - 14*armZZbar[90];
   armZZbar[118]=1./3.*armZZbar[118] + 4*armZZbar[89];
   armZZbar[118]=MMH*armZZbar[73]*armZZbar[118];
   armZZbar[103]=armZZbar[104] + armZZbar[103] + armZZbar[111] + 
   armZZbar[105] + 2*armZZbar[115] + armZZbar[118];
   armZZbar[103]=armZZbar[4]*armZZbar[103];
   armZZbar[104]=11*armZZbar[85];
   armZZbar[105]=55*armZZbar[90] - 501./4. + armZZbar[104];
   armZZbar[111]= - 3*armZZbar[87];
   armZZbar[115]= - 18*armZZbar[86];
   armZZbar[118]=4*armZZbar[84];
   armZZbar[121]= - 2*armZZbar[96];
   armZZbar[134]= - 18*armZZbar[43];
   armZZbar[136]= - 2*armZZbar[47];
   armZZbar[105]=armZZbar[210] + armZZbar[335] + 27./2.*armZZbar[15] + 
   armZZbar[136] + armZZbar[134] + 33*armZZbar[89] + armZZbar[121] + 
   armZZbar[118] + armZZbar[115] + 1./2.*armZZbar[105] + armZZbar[111];
   armZZbar[105]=armZZbar[101]*armZZbar[105];
   armZZbar[139]=2*armZZbar[95] + 1./2.*armZZbar[83];
   armZZbar[139]=11*armZZbar[139] - 255./8.*armZZbar[82];
   armZZbar[142]=2*armZZbar[78];
   armZZbar[143]= - 44*armZZbar[16] + 88*armZZbar[89] + 2*armZZbar[96]
    - 44*armZZbar[90] + 41./4. + 44*armZZbar[97];
   armZZbar[143]=armZZbar[100]*armZZbar[143];
   armZZbar[144]= - armZZbar[15]*armZZbar[100];
   armZZbar[147]= - armZZbar[13]*armZZbar[100];
   armZZbar[105]=armZZbar[105] + 44*armZZbar[147] + 2*armZZbar[144] + 
   armZZbar[143] + armZZbar[142] + 22*armZZbar[81] + 3*armZZbar[139] + 
   55*armZZbar[80];
   armZZbar[105]=MMZ*armZZbar[105];
   armZZbar[139]= - 11./4.*armZZbar[13] + armZZbar[157] + armZZbar[47]
    + 445./12. + armZZbar[333];
   armZZbar[139]=armZZbar[17]*armZZbar[139];
   armZZbar[143]= - 10*armZZbar[91];
   armZZbar[144]=4*armZZbar[93];
   armZZbar[149]= - 43./12.*armZZbar[75];
   armZZbar[155]= - 6*armZZbar[22];
   armZZbar[139]=armZZbar[139] - 32*armZZbar[19] + 56./3.*armZZbar[21]
    + 293./6.*armZZbar[20] + armZZbar[155] + armZZbar[149] + 
   armZZbar[144] - armZZbar[76] + armZZbar[143] - 11./2.*armZZbar[77]
    + armZZbar[183] + armZZbar[340];
   armZZbar[139]=armZZbar[101]*armZZbar[139];
   armZZbar[160]=1./2. + 6*armZZbar[43];
   armZZbar[167]=2*armZZbar[47];
   armZZbar[160]=armZZbar[306] + 3*armZZbar[160] + armZZbar[167];
   armZZbar[169]=armZZbar[160] + armZZbar[293];
   armZZbar[169]=armZZbar[101]*armZZbar[169];
   armZZbar[169]= - 2*armZZbar[100] + armZZbar[169];
   armZZbar[169]=MMZ*armZZbar[169];
   armZZbar[170]=179./9. + armZZbar[333];
   armZZbar[171]= - 2*armZZbar[45];
   armZZbar[173]=1 + armZZbar[114];
   armZZbar[173]=2*armZZbar[10]*armZZbar[173];
   armZZbar[169]=armZZbar[173] + armZZbar[169] + 1./2.*armZZbar[222] - 
   armZZbar[11] + armZZbar[323] + armZZbar[221] + 1./2.*armZZbar[170]
    + armZZbar[171];
   armZZbar[169]=armZZbar[10]*armZZbar[169];
   armZZbar[170]=403./6. + armZZbar[305];
   armZZbar[176]=13./2.*armZZbar[12];
   armZZbar[178]=9*armZZbar[10];
   armZZbar[179]= - 3./2.*armZZbar[11];
   armZZbar[170]=armZZbar[178] + armZZbar[179] + 245./8.*armZZbar[13]
    + armZZbar[176] + 1./2.*armZZbar[170] + armZZbar[224];
   armZZbar[170]=armZZbar[18]*armZZbar[170];
   armZZbar[180]=19 + armZZbar[224];
   armZZbar[180]=MMZ*armZZbar[180];
   armZZbar[181]=9./2.*armZZbar[19];
   armZZbar[161]=1./2.*armZZbar[180] + armZZbar[181] + armZZbar[161];
   armZZbar[161]=armZZbar[10]*armZZbar[161];
   armZZbar[180]=19./4.*armZZbar[12];
   armZZbar[184]=377./4.*armZZbar[13] + armZZbar[180] + 392./3. + 
   armZZbar[306];
   armZZbar[184]=armZZbar[19]*armZZbar[184];
   armZZbar[188]=armZZbar[323] - 11 + armZZbar[300];
   armZZbar[188]=armZZbar[17]*armZZbar[188];
   armZZbar[189]=armZZbar[335] - 1 + armZZbar[305];
   armZZbar[189]=armZZbar[19]*armZZbar[189];
   armZZbar[192]= - 1 + armZZbar[167];
   armZZbar[192]=2*armZZbar[192] - armZZbar[45];
   armZZbar[192]=MMZ*armZZbar[192];
   armZZbar[189]=armZZbar[189] + armZZbar[192];
   armZZbar[189]=MMZ*armZZbar[189];
   armZZbar[192]=1 + armZZbar[335];
   armZZbar[194]=armZZbar[18]*MMZ*armZZbar[192];
   armZZbar[189]=armZZbar[189] + armZZbar[194];
   armZZbar[189]=armZZbar[3]*armZZbar[189];
   armZZbar[196]= - 11*armZZbar[77] - 3./2.*armZZbar[75];
   armZZbar[196]=7*armZZbar[21] + 25./4.*armZZbar[20] + 1./2.*
   armZZbar[196] + 11*armZZbar[22];
   armZZbar[197]= - 4*armZZbar[11] + 152*armZZbar[13] + 151./4.*
   armZZbar[12] + armZZbar[289] - 223./9. - armZZbar[47];
   armZZbar[197]=MMZ*armZZbar[197];
   armZZbar[170]=3*armZZbar[189] + armZZbar[170] + armZZbar[161] + 
   armZZbar[197] + 3./2.*armZZbar[188] + 3*armZZbar[196] + 
   armZZbar[184];
   armZZbar[170]=armZZbar[3]*armZZbar[170];
   armZZbar[184]=armZZbar[335] + armZZbar[136] - 143./6. + 
   armZZbar[134];
   armZZbar[188]=armZZbar[184] + armZZbar[210];
   armZZbar[188]=armZZbar[101]*armZZbar[188];
   armZZbar[189]= - armZZbar[10]*armZZbar[101];
   armZZbar[196]=3./2.*armZZbar[189];
   armZZbar[188]=armZZbar[196] - armZZbar[100] + armZZbar[188];
   armZZbar[188]=armZZbar[18]*armZZbar[188];
   armZZbar[197]= - 185773./324. + 173*armZZbar[88];
   armZZbar[200]=8./3.*armZZbar[87];
   armZZbar[202]= - 9./2.*armZZbar[86];
   armZZbar[203]= - 27./4.*armZZbar[43];
   armZZbar[205]=79./12.*armZZbar[20] - 227./6.*armZZbar[22] + 43./6.*
   armZZbar[76] - 2*armZZbar[93];
   armZZbar[205]=armZZbar[100]*armZZbar[205];
   armZZbar[206]=7./4.*armZZbar[45];
   armZZbar[207]= - 5619./4.*armZZbar[13] - 115177./54. - 363*
   armZZbar[12];
   armZZbar[207]=armZZbar[13]*armZZbar[207];
   armZZbar[209]=67./9. + armZZbar[273];
   armZZbar[209]=armZZbar[11]*armZZbar[209];
   armZZbar[212]=56./3.*armZZbar[100] + 55./2.*armZZbar[147];
   armZZbar[212]=armZZbar[19]*armZZbar[212];
   armZZbar[219]=armZZbar[13]*armZZbar[100];
   armZZbar[220]= - 251./3.*armZZbar[100] + 55*armZZbar[219];
   armZZbar[220]=armZZbar[17]*armZZbar[220];
   armZZbar[223]=armZZbar[21]*armZZbar[100];
   armZZbar[105]=armZZbar[170] + armZZbar[188] + armZZbar[169] + 
   armZZbar[105] + armZZbar[139] + 1./4.*armZZbar[220] + armZZbar[212]
    + 1./2.*armZZbar[209] + 1./4.*armZZbar[207] + 149./108.*
   armZZbar[12] + armZZbar[206] + 2*armZZbar[223] + armZZbar[307] + 
   armZZbar[156] + armZZbar[205] + armZZbar[203] - 1891./24.*
   armZZbar[16] + 77./2.*armZZbar[89] + 28./3.*armZZbar[96] + 
   armZZbar[202] + armZZbar[200] - 51./2.*armZZbar[90] - 22./3.*
   armZZbar[97] + 5./8.*armZZbar[197] - 11./3.*armZZbar[85];
   armZZbar[105]=armZZbar[73]*armZZbar[105];
   armZZbar[139]= - 363./4.*armZZbar[13] + armZZbar[145] + 
   armZZbar[306] - 491./6. + armZZbar[305];
   armZZbar[139]=armZZbar[19]*armZZbar[139];
   armZZbar[145]= - 99./2.*armZZbar[13];
   armZZbar[169]=armZZbar[145] + armZZbar[306] + 115./2. + 
   armZZbar[305];
   armZZbar[169]=armZZbar[17]*armZZbar[169];
   armZZbar[110]= - 1023./4.*armZZbar[13] + armZZbar[110] + 
   armZZbar[306] + 133./4. + armZZbar[305];
   armZZbar[110]=MMZ*armZZbar[110];
   armZZbar[170]=15./4.*armZZbar[10] + armZZbar[179] - 231./4.*
   armZZbar[13] + 165./8.*armZZbar[12] + armZZbar[306] - 475./24. + 
   armZZbar[305];
   armZZbar[170]=armZZbar[18]*armZZbar[170];
   armZZbar[188]=armZZbar[300] + 13./2. + armZZbar[47];
   armZZbar[188]=MMZ*armZZbar[188];
   armZZbar[181]=3./2.*armZZbar[188] + armZZbar[181] + armZZbar[17];
   armZZbar[181]=armZZbar[10]*armZZbar[181];
   armZZbar[188]=armZZbar[289] + 1 - armZZbar[47];
   armZZbar[197]=armZZbar[19]*armZZbar[188];
   armZZbar[188]=MMZ*armZZbar[188];
   armZZbar[197]=armZZbar[197] + armZZbar[188];
   armZZbar[197]=MMZ*armZZbar[197];
   armZZbar[188]=armZZbar[18]*armZZbar[188];
   armZZbar[188]=armZZbar[197] + 1./2.*armZZbar[188];
   armZZbar[188]=armZZbar[3]*armZZbar[188];
   armZZbar[197]=armZZbar[152] + 11*armZZbar[77] - armZZbar[76];
   armZZbar[197]= - 5./2.*armZZbar[21] - armZZbar[20] + 1./4.*
   armZZbar[197] - 5*armZZbar[22];
   armZZbar[110]=9*armZZbar[188] + armZZbar[170] + armZZbar[181] + 1./2.
   *armZZbar[110] + 1./4.*armZZbar[169] + 9*armZZbar[197] + 
   armZZbar[139];
   armZZbar[110]=armZZbar[3]*armZZbar[110];
   armZZbar[139]= - 555./4. - 11*armZZbar[85];
   armZZbar[139]= - armZZbar[87] + 1./2.*armZZbar[139] - 11*
   armZZbar[90];
   armZZbar[139]=1./2.*armZZbar[139] + armZZbar[297];
   armZZbar[169]=2*armZZbar[84];
   armZZbar[170]= - 3./2.*armZZbar[47];
   armZZbar[181]= - 33./4.*armZZbar[13];
   armZZbar[139]=armZZbar[181] + armZZbar[298] + 73./4.*armZZbar[15] + 
   armZZbar[170] + armZZbar[285] - 99./2.*armZZbar[89] - 3./2.*
   armZZbar[96] + 3*armZZbar[139] + armZZbar[169];
   armZZbar[139]=armZZbar[101]*armZZbar[139];
   armZZbar[188]= - 11./4.*armZZbar[81] + 153./16.*armZZbar[82] - 11*
   armZZbar[80];
   armZZbar[197]= - armZZbar[96] + 11*armZZbar[90] + 5./2. - 11*
   armZZbar[97];
   armZZbar[197]=11./2.*armZZbar[16] + 1./2.*armZZbar[197] - 11*
   armZZbar[89];
   armZZbar[197]=armZZbar[100]*armZZbar[197];
   armZZbar[205]=armZZbar[15]*armZZbar[100];
   armZZbar[139]=armZZbar[139] + 33./2.*armZZbar[219] + 3./2.*
   armZZbar[205] + 3*armZZbar[197] + 3*armZZbar[188] + armZZbar[78];
   armZZbar[139]=MMZ*armZZbar[139];
   armZZbar[188]=1./4. + 3*armZZbar[43];
   armZZbar[197]=11./4.*armZZbar[13] + armZZbar[251] + armZZbar[188] + 
   armZZbar[236];
   armZZbar[197]=armZZbar[101]*armZZbar[197];
   armZZbar[197]=1./2.*armZZbar[100] + armZZbar[197];
   armZZbar[197]=MMZ*armZZbar[197];
   armZZbar[207]= - 625./9. + armZZbar[333];
   armZZbar[209]=1./2. + armZZbar[114];
   armZZbar[209]=armZZbar[10]*armZZbar[209];
   armZZbar[175]=armZZbar[209] + 3*armZZbar[197] + armZZbar[208] - 
   armZZbar[11] - 121./4.*armZZbar[13] + armZZbar[175] - armZZbar[45]
    + 1./4.*armZZbar[207] + armZZbar[136];
   armZZbar[175]=armZZbar[10]*armZZbar[175];
   armZZbar[197]=33./4.*armZZbar[13] + armZZbar[157] + 3./2.*
   armZZbar[47] + 1105./12. + armZZbar[333];
   armZZbar[197]=armZZbar[17]*armZZbar[197];
   armZZbar[207]=11./2.*armZZbar[77] + armZZbar[342] - 11*armZZbar[92];
   armZZbar[209]= - 5*armZZbar[91];
   armZZbar[212]=2*armZZbar[93];
   armZZbar[220]= - 43./24.*armZZbar[75];
   armZZbar[197]=1./2.*armZZbar[197] + 57./2.*armZZbar[19] - 5./3.*
   armZZbar[21] + 181./6.*armZZbar[20] - 5./2.*armZZbar[22] + 
   armZZbar[220] + armZZbar[212] - 3./4.*armZZbar[76] + 3./2.*
   armZZbar[207] + armZZbar[209];
   armZZbar[197]=armZZbar[101]*armZZbar[197];
   armZZbar[170]=armZZbar[181] + armZZbar[298] + armZZbar[170] - 149./
   12. + armZZbar[285];
   armZZbar[170]=armZZbar[101]*armZZbar[170];
   armZZbar[181]=3./4.*armZZbar[189];
   armZZbar[170]=armZZbar[181] + armZZbar[100] + armZZbar[170];
   armZZbar[170]=armZZbar[18]*armZZbar[170];
   armZZbar[189]= - 15821./162. - 435*armZZbar[88];
   armZZbar[104]=1./2.*armZZbar[189] + armZZbar[104];
   armZZbar[104]=1./2.*armZZbar[104] + 11*armZZbar[97];
   armZZbar[189]= - 21./4.*armZZbar[20] + 37./2.*armZZbar[22] - 3*
   armZZbar[76] + armZZbar[212];
   armZZbar[189]=armZZbar[100]*armZZbar[189];
   armZZbar[207]= - armZZbar[21]*armZZbar[100];
   armZZbar[225]=1797*armZZbar[13] + 16301./9. - 1089./2.*armZZbar[12];
   armZZbar[225]=armZZbar[13]*armZZbar[225];
   armZZbar[210]=95./9. + armZZbar[210];
   armZZbar[210]=armZZbar[11]*armZZbar[210];
   armZZbar[226]= - 5*armZZbar[100] + 11./2.*armZZbar[219];
   armZZbar[226]=armZZbar[19]*armZZbar[226];
   armZZbar[228]=15*armZZbar[100] + 11*armZZbar[147];
   armZZbar[228]=armZZbar[17]*armZZbar[228];
   armZZbar[104]=armZZbar[110] + armZZbar[170] + armZZbar[175] + 
   armZZbar[139] + armZZbar[197] + 3./4.*armZZbar[228] + 3*
   armZZbar[226] + 1./2.*armZZbar[210] + 1./8.*armZZbar[225] - 5027./72.
   *armZZbar[12] + 5./8.*armZZbar[45] + 2*armZZbar[207] + 7./6.*
   armZZbar[15] + 5./4.*armZZbar[47] + armZZbar[189] - 45./8.*
   armZZbar[43] + 391./8.*armZZbar[16] - 33*armZZbar[89] - 3*
   armZZbar[96] + armZZbar[140] + 1./2.*armZZbar[104] + 4./3.*
   armZZbar[87];
   armZZbar[104]=armZZbar[73]*armZZbar[104];
   armZZbar[110]=armZZbar[119] + armZZbar[21] - armZZbar[8] + 
   armZZbar[198];
   armZZbar[110]=armZZbar[101]*armZZbar[110];
   armZZbar[139]= - 3*armZZbar[24] - armZZbar[26];
   armZZbar[139]=2*armZZbar[139] + armZZbar[25];
   armZZbar[170]= - armZZbar[15] - 3*armZZbar[14] + armZZbar[28] - 3 - 
   2*armZZbar[30];
   armZZbar[170]=armZZbar[101]*armZZbar[170];
   armZZbar[139]=2*armZZbar[139] + armZZbar[170];
   armZZbar[139]=MMZ*armZZbar[139];
   armZZbar[175]=1./9. + armZZbar[177];
   armZZbar[187]=armZZbar[175] + armZZbar[187];
   armZZbar[189]=armZZbar[1]*armZZbar[187];
   armZZbar[187]=armZZbar[39]*armZZbar[187];
   armZZbar[197]= - 37*armZZbar[32] - 3989./432. + 33*armZZbar[29];
   armZZbar[198]= - 7./3. + armZZbar[7];
   armZZbar[198]=armZZbar[13]*armZZbar[198];
   armZZbar[139]=armZZbar[187] + 2./3.*armZZbar[189] + armZZbar[139] + 
   armZZbar[110] + armZZbar[216] + 11./4.*armZZbar[198] + 113./18.*
   armZZbar[7] + armZZbar[215] + armZZbar[15] - 5./4.*armZZbar[177] - 
   33./2.*armZZbar[16] + armZZbar[174] - armZZbar[28] + 1./2.*
   armZZbar[197] + armZZbar[172];
   armZZbar[139]=armZZbar[39]*armZZbar[139];
   armZZbar[113]=armZZbar[243] - armZZbar[14] + armZZbar[125] - 1 + 
   armZZbar[113];
   armZZbar[113]=armZZbar[101]*armZZbar[113];
   armZZbar[125]= - armZZbar[24] - 1./3.*armZZbar[26];
   armZZbar[125]=2*armZZbar[125] + 1./3.*armZZbar[25];
   armZZbar[125]=2*armZZbar[125] + armZZbar[113];
   armZZbar[125]=MMZ*armZZbar[125];
   armZZbar[172]= - 37./3.*armZZbar[32] - 3989./1296. + 11*armZZbar[29]
   ;
   armZZbar[125]=1./9.*armZZbar[189] + armZZbar[125] + 1./3.*
   armZZbar[110] + armZZbar[331] + 11./12.*armZZbar[198] + 113./54.*
   armZZbar[7] - 55./36.*armZZbar[12] + armZZbar[227] + armZZbar[259]
    - 11./2.*armZZbar[16] + armZZbar[257] + armZZbar[255] + 1./2.*
   armZZbar[172] - 5./3.*armZZbar[31];
   armZZbar[125]=armZZbar[1]*armZZbar[125];
   armZZbar[172]=99./2.*armZZbar[13] - 17 - 99./2.*armZZbar[12];
   armZZbar[187]=3*armZZbar[11];
   armZZbar[122]=armZZbar[122] + 1./2.*armZZbar[172] + armZZbar[187];
   armZZbar[122]=armZZbar[18]*armZZbar[122];
   armZZbar[172]=99*armZZbar[13] + 17 - 99*armZZbar[12];
   armZZbar[172]=1./4.*armZZbar[172] + armZZbar[187];
   armZZbar[172]=armZZbar[19]*armZZbar[172];
   armZZbar[186]=MMZ*armZZbar[186];
   armZZbar[187]=armZZbar[10]*armZZbar[120];
   armZZbar[122]=1./2.*armZZbar[122] + 3*armZZbar[187] + armZZbar[172]
    + 3*armZZbar[186];
   armZZbar[122]=armZZbar[3]*armZZbar[122];
   armZZbar[172]= - 1./2.*armZZbar[93] + armZZbar[22];
   armZZbar[172]=armZZbar[100]*armZZbar[172];
   armZZbar[172]=armZZbar[172] + 1./2.*armZZbar[223];
   armZZbar[135]=armZZbar[135] - 1./2.*armZZbar[21] + armZZbar[301] - 
   armZZbar[22];
   armZZbar[135]=armZZbar[101]*armZZbar[135];
   armZZbar[145]=armZZbar[145] - 59./3. + 99./2.*armZZbar[12];
   armZZbar[145]=armZZbar[13]*armZZbar[145];
   armZZbar[186]= - 3./2.*armZZbar[13];
   armZZbar[187]= - 7./9. + armZZbar[186];
   armZZbar[187]=armZZbar[11]*armZZbar[187];
   armZZbar[189]= - armZZbar[100] + armZZbar[101];
   armZZbar[197]=MMZ*armZZbar[189];
   armZZbar[198]=3*armZZbar[13];
   armZZbar[207]=armZZbar[198] + 7./9. - 3./2.*armZZbar[12];
   armZZbar[207]= - armZZbar[10] + 11./2.*armZZbar[207] + armZZbar[11];
   armZZbar[207]=armZZbar[10]*armZZbar[207];
   armZZbar[189]=armZZbar[18]*armZZbar[189];
   armZZbar[122]=armZZbar[122] + 3./4.*armZZbar[189] + armZZbar[207] + 
   9./16.*armZZbar[197] + 3*armZZbar[135] + 9./2.*armZZbar[222] + 11./2.
   *armZZbar[187] + 11./8.*armZZbar[145] + 3*armZZbar[172] + 649./24.*
   armZZbar[12];
   armZZbar[122]=armZZbar[73]*armZZbar[122];
   armZZbar[135]=armZZbar[1]*armZZbar[120];
   armZZbar[120]=armZZbar[39]*armZZbar[120];
   armZZbar[124]=armZZbar[18]*armZZbar[124];
   armZZbar[120]=1./2.*armZZbar[124] + armZZbar[135] + 3*armZZbar[120];
   armZZbar[120]=armZZbar[6]*armZZbar[120];
   armZZbar[124]=armZZbar[19]*armZZbar[7];
   armZZbar[135]=MMZ*armZZbar[7];
   armZZbar[124]=armZZbar[124] + armZZbar[135];
   armZZbar[135]=armZZbar[1]*armZZbar[124];
   armZZbar[124]=armZZbar[39]*armZZbar[124];
   armZZbar[145]=armZZbar[182] + 3*armZZbar[193];
   armZZbar[145]=armZZbar[18]*armZZbar[145];
   armZZbar[120]=armZZbar[120] + 1./2.*armZZbar[145] + armZZbar[135] + 
   3*armZZbar[124];
   armZZbar[120]=armZZbar[3]*armZZbar[120];
   armZZbar[124]=11./2.*armZZbar[12] - 59./9.*armZZbar[7];
   armZZbar[135]= - 3*armZZbar[7];
   armZZbar[145]= - 1 + armZZbar[135];
   armZZbar[145]=armZZbar[13]*armZZbar[145];
   armZZbar[145]=armZZbar[124] + 11./2.*armZZbar[145];
   armZZbar[145]=1./3.*armZZbar[158] + 2./9.*armZZbar[154] + 1./2.*
   armZZbar[145] + armZZbar[151];
   armZZbar[145]=armZZbar[39]*armZZbar[145];
   armZZbar[158]=armZZbar[323] + 59./27. + armZZbar[319];
   armZZbar[158]=1./9.*armZZbar[320] + 1./2.*armZZbar[158] + 
   armZZbar[216];
   armZZbar[158]=armZZbar[1]*armZZbar[158];
   armZZbar[130]=armZZbar[141] + armZZbar[130] + 2./3.*armZZbar[320];
   armZZbar[130]=armZZbar[39]*armZZbar[130];
   armZZbar[141]= - 2./3.*armZZbar[1];
   armZZbar[172]=armZZbar[141] - armZZbar[39];
   armZZbar[172]=armZZbar[39]*armZZbar[172];
   armZZbar[182]=pow(armZZbar[1],2);
   armZZbar[172]= - 1./9.*armZZbar[182] + armZZbar[172];
   armZZbar[172]=armZZbar[6]*armZZbar[172];
   armZZbar[130]=armZZbar[172] + 2*armZZbar[322] + armZZbar[158] + 
   armZZbar[130];
   armZZbar[130]=armZZbar[6]*armZZbar[130];
   armZZbar[158]=armZZbar[13]*armZZbar[263];
   armZZbar[124]=1./3.*armZZbar[124] + 11./2.*armZZbar[158];
   armZZbar[124]=1./27.*armZZbar[154] + 1./2.*armZZbar[124] + 
   armZZbar[138];
   armZZbar[124]=armZZbar[1]*armZZbar[124];
   armZZbar[120]=armZZbar[122] + armZZbar[120] + armZZbar[130] + 
   armZZbar[321] + armZZbar[124] + armZZbar[145];
   armZZbar[120]=armZZbar[133]*armZZbar[120];
   armZZbar[122]= - 41./3. + armZZbar[163];
   armZZbar[117]=5./6.*armZZbar[122] + armZZbar[117];
   armZZbar[124]= - 1 - armZZbar[7];
   armZZbar[130]=armZZbar[1]*armZZbar[124];
   armZZbar[117]=1./9.*armZZbar[130] + 1./3.*armZZbar[250] + 1./6.*
   armZZbar[116] + 1./2.*armZZbar[117] + armZZbar[151];
   armZZbar[117]=armZZbar[1]*armZZbar[117];
   armZZbar[122]=5./2.*armZZbar[122] - 99*armZZbar[13];
   armZZbar[124]=armZZbar[39]*armZZbar[124];
   armZZbar[122]=armZZbar[124] + 2./3.*armZZbar[130] + armZZbar[250] + 
   armZZbar[185] + 1./2.*armZZbar[122] - armZZbar[11];
   armZZbar[122]=armZZbar[39]*armZZbar[122];
   armZZbar[124]=armZZbar[1]*armZZbar[217];
   armZZbar[130]=3 + armZZbar[114];
   armZZbar[138]=armZZbar[39]*armZZbar[130];
   armZZbar[124]=armZZbar[124] + armZZbar[138];
   armZZbar[124]=armZZbar[10]*armZZbar[124];
   armZZbar[138]=armZZbar[1]*armZZbar[101];
   armZZbar[145]=armZZbar[39]*armZZbar[101];
   armZZbar[151]=1./3.*armZZbar[138] + armZZbar[145];
   armZZbar[151]=armZZbar[18]*armZZbar[151];
   armZZbar[154]=5./4. + 2./3.*armZZbar[1];
   armZZbar[154]=armZZbar[1]*armZZbar[154];
   armZZbar[158]=4./3.*armZZbar[1];
   armZZbar[163]=2*armZZbar[39] + 5./4. + armZZbar[158];
   armZZbar[163]=armZZbar[39]*armZZbar[163];
   armZZbar[154]=1./3.*armZZbar[154] + armZZbar[163];
   armZZbar[154]=armZZbar[6]*armZZbar[154];
   armZZbar[117]=armZZbar[154] + armZZbar[151] + armZZbar[124] + 
   armZZbar[117] + armZZbar[122];
   armZZbar[117]=armZZbar[6]*armZZbar[117];
   armZZbar[106]=armZZbar[106] + armZZbar[214] - armZZbar[8];
   armZZbar[122]=armZZbar[128] + armZZbar[131] + armZZbar[132] + 
   armZZbar[127] + 3*armZZbar[106] + armZZbar[126];
   armZZbar[122]=armZZbar[39]*armZZbar[122];
   armZZbar[124]=1./6.*armZZbar[20];
   armZZbar[128]= - 1./6.*armZZbar[17];
   armZZbar[132]= - 4*MMZ;
   armZZbar[106]=armZZbar[132] + armZZbar[128] + 2./3.*armZZbar[19] + 
   armZZbar[21] + armZZbar[106] + armZZbar[124];
   armZZbar[106]=armZZbar[1]*armZZbar[106];
   armZZbar[151]=9./2.*MMZ + armZZbar[310] + 1./6.*armZZbar[17];
   armZZbar[151]=armZZbar[1]*armZZbar[151];
   armZZbar[119]=27./2.*MMZ + 9*armZZbar[19] + armZZbar[119];
   armZZbar[119]=armZZbar[39]*armZZbar[119];
   armZZbar[154]=armZZbar[18]*armZZbar[159];
   armZZbar[119]=2*armZZbar[154] + armZZbar[151] + armZZbar[119];
   armZZbar[119]=armZZbar[6]*armZZbar[119];
   armZZbar[151]= - 3./2.*armZZbar[7];
   armZZbar[154]=1 + armZZbar[151];
   armZZbar[154]=armZZbar[39]*armZZbar[154];
   armZZbar[154]=armZZbar[315] + armZZbar[154];
   armZZbar[154]=armZZbar[18]*armZZbar[154];
   armZZbar[106]=armZZbar[119] + armZZbar[154] + armZZbar[106] + 
   armZZbar[122];
   armZZbar[106]=armZZbar[3]*armZZbar[106];
   armZZbar[119]=2./3. - armZZbar[7];
   armZZbar[122]=armZZbar[1]*armZZbar[119];
   armZZbar[119]=armZZbar[39]*armZZbar[119];
   armZZbar[119]=1./3.*armZZbar[122] + armZZbar[119];
   armZZbar[119]=armZZbar[10]*armZZbar[119];
   armZZbar[122]= - 2./3.*armZZbar[78];
   armZZbar[154]=armZZbar[122] + 13./4.*armZZbar[81] + 3./2.*
   armZZbar[79] + 13*armZZbar[80];
   armZZbar[154]=MMH*armZZbar[73]*armZZbar[154];
   armZZbar[104]=armZZbar[120] + armZZbar[154] + armZZbar[104] + 
   armZZbar[106] + armZZbar[117] + armZZbar[119] + armZZbar[125] + 
   armZZbar[139];
   armZZbar[104]=armZZbar[133]*armZZbar[104];
   armZZbar[106]=28*armZZbar[24];
   armZZbar[117]=8*armZZbar[26];
   armZZbar[119]= - 4*armZZbar[25];
   armZZbar[120]=armZZbar[119] + armZZbar[117] + armZZbar[106] + 5./9.*
   armZZbar[27];
   armZZbar[120]=MMZ*armZZbar[120];
   armZZbar[125]= - 28*armZZbar[29];
   armZZbar[139]=37./4.*armZZbar[23];
   armZZbar[154]=83./9. + 9*armZZbar[7];
   armZZbar[154]=armZZbar[13]*armZZbar[154];
   armZZbar[159]=19./27.*armZZbar[7] - 37./81. - armZZbar[177];
   armZZbar[159]=armZZbar[1]*armZZbar[159];
   armZZbar[163]=11./27.*armZZbar[7] - 29./81. - armZZbar[177];
   armZZbar[163]=armZZbar[39]*armZZbar[163];
   armZZbar[120]=armZZbar[163] + 2./3.*armZZbar[159] + armZZbar[120] + 
   armZZbar[271] + 5./4.*armZZbar[154] + 61./18.*armZZbar[7] + 
   armZZbar[270] + armZZbar[177] + armZZbar[258] + 10*armZZbar[31] + 81.
   /2.*armZZbar[32] + armZZbar[139] + 24539./1944. + armZZbar[125];
   armZZbar[120]=armZZbar[39]*armZZbar[120];
   armZZbar[106]=armZZbar[119] + armZZbar[117] + armZZbar[106] + 
   armZZbar[27];
   armZZbar[106]=MMZ*armZZbar[106];
   armZZbar[117]=armZZbar[139] + 893./72. + armZZbar[125];
   armZZbar[119]=83./27. + armZZbar[150];
   armZZbar[119]=armZZbar[13]*armZZbar[119];
   armZZbar[125]=armZZbar[7] - 5./9. - armZZbar[177];
   armZZbar[125]=armZZbar[1]*armZZbar[125];
   armZZbar[106]=1./9.*armZZbar[125] + 1./3.*armZZbar[106] + 
   armZZbar[216] + 5./4.*armZZbar[119] + 61./54.*armZZbar[7] - 125./54.
   *armZZbar[12] + 1./3.*armZZbar[177] + 28./3.*armZZbar[16] + 10./3.*
   armZZbar[31] + 1./3.*armZZbar[117] + 27./2.*armZZbar[32];
   armZZbar[106]=armZZbar[1]*armZZbar[106];
   armZZbar[117]= - 5*armZZbar[19] - 19./3.*MMZ;
   armZZbar[117]=armZZbar[1]*armZZbar[117];
   armZZbar[119]= - 29*armZZbar[19] - 41*MMZ;
   armZZbar[119]=armZZbar[39]*armZZbar[119];
   armZZbar[125]=1./2.*armZZbar[18]*armZZbar[375];
   armZZbar[117]=armZZbar[125] + armZZbar[117] + 1./3.*armZZbar[119];
   armZZbar[117]=armZZbar[6]*armZZbar[117];
   armZZbar[119]= - 2*armZZbar[7];
   armZZbar[139]=19./3. + armZZbar[119];
   armZZbar[139]=MMZ*armZZbar[139];
   armZZbar[109]=armZZbar[109] + 2*armZZbar[139];
   armZZbar[109]=armZZbar[1]*armZZbar[109];
   armZZbar[119]=49./9. + armZZbar[119];
   armZZbar[119]=MMZ*armZZbar[119];
   armZZbar[119]=armZZbar[211] + 2*armZZbar[119];
   armZZbar[119]=armZZbar[39]*armZZbar[119];
   armZZbar[139]=armZZbar[1]*armZZbar[234];
   armZZbar[150]=10./9. + armZZbar[151];
   armZZbar[150]=armZZbar[39]*armZZbar[150];
   armZZbar[139]=armZZbar[139] + armZZbar[150];
   armZZbar[139]=armZZbar[18]*armZZbar[139];
   armZZbar[109]=armZZbar[117] + armZZbar[139] + 1./3.*armZZbar[109] + 
   armZZbar[119];
   armZZbar[109]=armZZbar[3]*armZZbar[109];
   armZZbar[117]=194*armZZbar[13];
   armZZbar[119]= - 19*armZZbar[7];
   armZZbar[150]=55./3. + armZZbar[119];
   armZZbar[150]=armZZbar[1]*armZZbar[150];
   armZZbar[137]=47./3. + armZZbar[137];
   armZZbar[137]=armZZbar[39]*armZZbar[137];
   armZZbar[137]=1./3.*armZZbar[137] + 2./9.*armZZbar[150] - 11./3.*
   armZZbar[11] + armZZbar[117] + 2615./27. + 37./2.*armZZbar[12];
   armZZbar[137]=armZZbar[39]*armZZbar[137];
   armZZbar[117]=armZZbar[117] + 349./3. + 125./2.*armZZbar[12];
   armZZbar[150]=7./9. - armZZbar[7];
   armZZbar[150]=armZZbar[1]*armZZbar[150];
   armZZbar[117]=1./3.*armZZbar[150] + 1./9.*armZZbar[117] - 
   armZZbar[11];
   armZZbar[117]=armZZbar[1]*armZZbar[117];
   armZZbar[150]= - 1 + armZZbar[141];
   armZZbar[150]=armZZbar[1]*armZZbar[150];
   armZZbar[151]= - 2*armZZbar[39] - 1 + armZZbar[201];
   armZZbar[151]=armZZbar[39]*armZZbar[151];
   armZZbar[150]=1./3.*armZZbar[150] + armZZbar[151];
   armZZbar[150]=armZZbar[6]*armZZbar[150];
   armZZbar[117]=armZZbar[150] + armZZbar[117] + 1./3.*armZZbar[137];
   armZZbar[117]=armZZbar[6]*armZZbar[117];
   armZZbar[137]=3*armZZbar[79];
   armZZbar[150]= - 4./3.*armZZbar[78];
   armZZbar[151]=armZZbar[150] - 67./12.*armZZbar[81] + armZZbar[137]
    - 89./6.*armZZbar[80];
   armZZbar[151]=MMH*armZZbar[73]*armZZbar[151];
   armZZbar[104]=armZZbar[104] + armZZbar[151] + armZZbar[105] + 
   armZZbar[109] + armZZbar[117] + armZZbar[325] + armZZbar[106] + 
   armZZbar[120];
   armZZbar[104]=armZZbar[133]*armZZbar[104];
   armZZbar[105]=1./3.*armZZbar[90] - 2555./12. + armZZbar[148];
   armZZbar[106]=3./2.*armZZbar[13];
   armZZbar[105]=armZZbar[106] + armZZbar[335] + 35./2.*armZZbar[15] + 
   armZZbar[136] + armZZbar[134] + 9*armZZbar[89] + armZZbar[121] + 
   armZZbar[118] + armZZbar[115] + 1./2.*armZZbar[105] + armZZbar[111];
   armZZbar[105]=armZZbar[101]*armZZbar[105];
   armZZbar[109]=armZZbar[83] - 23./4.*armZZbar[82];
   armZZbar[109]=1./2.*armZZbar[109] + armZZbar[80];
   armZZbar[105]=armZZbar[105] + 1./3.*armZZbar[109] + armZZbar[142];
   armZZbar[105]=MMZ*armZZbar[105];
   armZZbar[109]=armZZbar[160] + armZZbar[186];
   armZZbar[109]=MMZ*armZZbar[101]*armZZbar[109];
   armZZbar[109]=armZZbar[173] + armZZbar[109] + armZZbar[208] - 
   armZZbar[11] + armZZbar[198] + armZZbar[221] + armZZbar[171] + 40./9.
    + armZZbar[278];
   armZZbar[109]=armZZbar[10]*armZZbar[109];
   armZZbar[115]=187./6. + armZZbar[305];
   armZZbar[115]=armZZbar[178] + armZZbar[179] + 53./8.*armZZbar[13] + 
   armZZbar[176] + 1./2.*armZZbar[115] + armZZbar[224];
   armZZbar[115]=armZZbar[18]*armZZbar[115];
   armZZbar[117]=1./4.*armZZbar[13];
   armZZbar[118]=armZZbar[117] + armZZbar[180] + 8 + armZZbar[306];
   armZZbar[118]=armZZbar[19]*armZZbar[118];
   armZZbar[120]= - armZZbar[11] + 20./3.*armZZbar[13] + armZZbar[180]
    + armZZbar[289] - 25./6. + armZZbar[47];
   armZZbar[120]=MMZ*armZZbar[120];
   armZZbar[133]=armZZbar[19]*armZZbar[192];
   armZZbar[134]=2 + armZZbar[379];
   armZZbar[134]=MMZ*armZZbar[134];
   armZZbar[133]=3*armZZbar[133] + armZZbar[134];
   armZZbar[133]=MMZ*armZZbar[133];
   armZZbar[133]=armZZbar[133] + 3*armZZbar[194];
   armZZbar[133]=armZZbar[3]*armZZbar[133];
   armZZbar[134]= - armZZbar[77] + armZZbar[152];
   armZZbar[134]=1./2.*armZZbar[134] + armZZbar[22];
   armZZbar[142]=5./3.*armZZbar[13] - 5./3. + armZZbar[306];
   armZZbar[142]=armZZbar[17]*armZZbar[142];
   armZZbar[115]=armZZbar[133] + armZZbar[115] + armZZbar[161] + 
   armZZbar[120] + 1./2.*armZZbar[142] + armZZbar[118] + 9*armZZbar[21]
    + 9*armZZbar[134] + 37./12.*armZZbar[20];
   armZZbar[115]=armZZbar[3]*armZZbar[115];
   armZZbar[118]= - 3./4.*armZZbar[13] + armZZbar[157] + armZZbar[47]
    + 565./12. + armZZbar[333];
   armZZbar[118]=armZZbar[17]*armZZbar[118];
   armZZbar[120]= - 1./2.*armZZbar[77] + armZZbar[342] + armZZbar[92];
   armZZbar[118]=armZZbar[118] - 16*armZZbar[19] + 44./3.*armZZbar[21]
    + 305./6.*armZZbar[20] + armZZbar[155] + armZZbar[149] + 
   armZZbar[144] - armZZbar[76] + 3*armZZbar[120] + armZZbar[143];
   armZZbar[118]=armZZbar[101]*armZZbar[118];
   armZZbar[120]=armZZbar[184] + armZZbar[106];
   armZZbar[120]=armZZbar[101]*armZZbar[120];
   armZZbar[120]=armZZbar[120] + armZZbar[196];
   armZZbar[120]=armZZbar[18]*armZZbar[120];
   armZZbar[133]= - 75473./324. + 33*armZZbar[88];
   armZZbar[134]= - armZZbar[22] + armZZbar[126];
   armZZbar[134]=armZZbar[100]*armZZbar[134];
   armZZbar[142]=29./2. - 19*armZZbar[12];
   armZZbar[142]=1./9.*armZZbar[142] + 13./4.*armZZbar[13];
   armZZbar[142]=armZZbar[13]*armZZbar[142];
   armZZbar[143]=7./3. + armZZbar[240];
   armZZbar[143]=armZZbar[11]*armZZbar[143];
   armZZbar[144]=5*armZZbar[100] + 1./2.*armZZbar[147];
   armZZbar[144]=armZZbar[19]*armZZbar[144];
   armZZbar[148]= - armZZbar[100] + armZZbar[219];
   armZZbar[148]=armZZbar[17]*armZZbar[148];
   armZZbar[105]=armZZbar[115] + armZZbar[120] + armZZbar[109] + 
   armZZbar[105] + armZZbar[118] + 1./12.*armZZbar[148] + 1./3.*
   armZZbar[144] + 1./6.*armZZbar[143] + 1./4.*armZZbar[142] + 77./108.
   *armZZbar[12] + armZZbar[206] + 8./3.*armZZbar[15] + armZZbar[156]
    + 1./6.*armZZbar[134] + armZZbar[203] - 33./8.*armZZbar[16] + 19./6.
   *armZZbar[89] + 5./6.*armZZbar[96] + armZZbar[202] + armZZbar[200]
    - 5./18.*armZZbar[90] + 1./8.*armZZbar[133] - 5./3.*armZZbar[85];
   armZZbar[105]=armZZbar[73]*armZZbar[105];
   armZZbar[109]= - 877./4. + 1./3.*armZZbar[85];
   armZZbar[109]=1./2.*armZZbar[109] + armZZbar[111];
   armZZbar[111]=1./12.*armZZbar[13];
   armZZbar[109]=armZZbar[111] + armZZbar[298] + 125./12.*armZZbar[15]
    + armZZbar[136] + armZZbar[285] + 1./2.*armZZbar[89] + 
   armZZbar[121] + armZZbar[169] + 1./2.*armZZbar[109] - 9*armZZbar[86]
   ;
   armZZbar[109]=armZZbar[101]*armZZbar[109];
   armZZbar[109]=armZZbar[109] + 3./16.*armZZbar[82] + armZZbar[78];
   armZZbar[109]=MMZ*armZZbar[109];
   armZZbar[115]= - armZZbar[77] - 9./2.*armZZbar[75];
   armZZbar[115]=armZZbar[123] + 7./3.*armZZbar[20] + 1./2.*
   armZZbar[115] + armZZbar[22];
   armZZbar[118]=pow(MMZ,2);
   armZZbar[120]=armZZbar[118]*armZZbar[192];
   armZZbar[120]=armZZbar[120] + 3./2.*armZZbar[194];
   armZZbar[120]=armZZbar[3]*armZZbar[120];
   armZZbar[123]= - 8./3. + 1./4.*armZZbar[12];
   armZZbar[123]=armZZbar[19]*armZZbar[123];
   armZZbar[133]=1./3.*armZZbar[13] - 1./3. + armZZbar[224];
   armZZbar[133]=armZZbar[17]*armZZbar[133];
   armZZbar[106]=armZZbar[106] - 17./3.*armZZbar[12] - 21./2. - 
   armZZbar[45];
   armZZbar[106]=MMZ*armZZbar[106];
   armZZbar[134]=7 + armZZbar[45];
   armZZbar[134]=MMZ*armZZbar[134];
   armZZbar[134]=armZZbar[17] + 3./4.*armZZbar[134];
   armZZbar[134]=armZZbar[10]*armZZbar[134];
   armZZbar[142]=21./2.*armZZbar[10] + 3./4.*armZZbar[13] - 17./4.*
   armZZbar[12] + 13./12. + armZZbar[224];
   armZZbar[142]=armZZbar[18]*armZZbar[142];
   armZZbar[106]=armZZbar[120] + 1./2.*armZZbar[142] + armZZbar[134] + 
   1./4.*armZZbar[106] + 1./8.*armZZbar[133] + 1./2.*armZZbar[115] + 
   armZZbar[123];
   armZZbar[106]=armZZbar[3]*armZZbar[106];
   armZZbar[115]= - 1./6.*armZZbar[77] + armZZbar[183] + 1./3.*
   armZZbar[92];
   armZZbar[120]=605./12. + armZZbar[333];
   armZZbar[120]= - 1./24.*armZZbar[13] + 3./8.*armZZbar[45] + 1./2.*
   armZZbar[120] + armZZbar[47];
   armZZbar[120]=armZZbar[17]*armZZbar[120];
   armZZbar[107]=armZZbar[120] - 13./3.*armZZbar[19] + 20./3.*
   armZZbar[21] + 105./4.*armZZbar[20] + armZZbar[107] + armZZbar[220]
    + armZZbar[212] - armZZbar[76] + 1./2.*armZZbar[115] + 
   armZZbar[209];
   armZZbar[107]=armZZbar[101]*armZZbar[107];
   armZZbar[115]= - 1./12.*armZZbar[13] + armZZbar[157] + 3*
   armZZbar[188] + armZZbar[167];
   armZZbar[115]=MMZ*armZZbar[101]*armZZbar[115];
   armZZbar[120]= - 61./9. + armZZbar[333];
   armZZbar[120]=1./2.*armZZbar[120] + armZZbar[47];
   armZZbar[114]=3./2. + armZZbar[114];
   armZZbar[114]=armZZbar[10]*armZZbar[114];
   armZZbar[114]=armZZbar[114] + armZZbar[115] + armZZbar[117] + 
   armZZbar[164] + 1./2.*armZZbar[120] - armZZbar[45];
   armZZbar[114]=armZZbar[10]*armZZbar[114];
   armZZbar[111]=armZZbar[111] + armZZbar[298] + armZZbar[136] - 155./
   12. + armZZbar[285];
   armZZbar[111]=armZZbar[101]*armZZbar[111];
   armZZbar[111]=armZZbar[111] + armZZbar[181];
   armZZbar[111]=armZZbar[18]*armZZbar[111];
   armZZbar[115]= - 6283./432. + armZZbar[88];
   armZZbar[115]=7*armZZbar[115] - 1./2.*armZZbar[85];
   armZZbar[115]=1./2.*armZZbar[115] + 4*armZZbar[87];
   armZZbar[117]= - 17./16.*armZZbar[13] - 1./3. - 19./16.*armZZbar[12]
   ;
   armZZbar[117]=armZZbar[13]*armZZbar[117];
   armZZbar[106]=armZZbar[106] + armZZbar[111] + armZZbar[114] + 
   armZZbar[109] + armZZbar[107] + armZZbar[216] + 1./9.*armZZbar[117]
    + 23./216.*armZZbar[12] + 7./8.*armZZbar[45] + 3./2.*armZZbar[15]
    + armZZbar[191] - 27./8.*armZZbar[43] - 7./6.*armZZbar[16] + 1./6.*
   armZZbar[89] + armZZbar[299] + 1./3.*armZZbar[115] + armZZbar[140];
   armZZbar[106]=armZZbar[73]*armZZbar[106];
   armZZbar[107]=6*armZZbar[25] + armZZbar[170];
   armZZbar[107]=MMZ*armZZbar[107];
   armZZbar[109]=armZZbar[1]*armZZbar[175];
   armZZbar[107]=armZZbar[109] + armZZbar[107] + armZZbar[110] - 1./18.
   *armZZbar[13] + armZZbar[162] + armZZbar[15] + 15./4.*armZZbar[177]
    + armZZbar[174] - armZZbar[28] + 1901./864. - 15*armZZbar[31];
   armZZbar[107]=armZZbar[1]*armZZbar[107];
   armZZbar[111]= - 25./3. + armZZbar[204];
   armZZbar[111]=1./2.*armZZbar[111] + armZZbar[13];
   armZZbar[111]=1./3.*armZZbar[111] + armZZbar[116];
   armZZbar[111]=armZZbar[141] + 1./2.*armZZbar[111] + armZZbar[250];
   armZZbar[111]=armZZbar[1]*armZZbar[111];
   armZZbar[114]= - 115./3. - 187*armZZbar[12];
   armZZbar[114]=1./2.*armZZbar[114] + armZZbar[323];
   armZZbar[114]=1./3.*armZZbar[114] + 11*armZZbar[116];
   armZZbar[114]= - 242./27.*armZZbar[39] - 44./3.*armZZbar[1] + 1./2.*
   armZZbar[114] + 11*armZZbar[250];
   armZZbar[114]=armZZbar[39]*armZZbar[114];
   armZZbar[115]=armZZbar[1]*armZZbar[130];
   armZZbar[116]=armZZbar[39]*armZZbar[217];
   armZZbar[115]=armZZbar[115] + 11./3.*armZZbar[116];
   armZZbar[115]=armZZbar[10]*armZZbar[115];
   armZZbar[116]=armZZbar[138] + 11./9.*armZZbar[145];
   armZZbar[116]=armZZbar[18]*armZZbar[116];
   armZZbar[117]= - 15./4. + armZZbar[1];
   armZZbar[117]=armZZbar[1]*armZZbar[117];
   armZZbar[120]=121./9.*armZZbar[39] - 685./36. + 22*armZZbar[1];
   armZZbar[120]=armZZbar[39]*armZZbar[120];
   armZZbar[117]=armZZbar[117] + 1./9.*armZZbar[120];
   armZZbar[117]=armZZbar[6]*armZZbar[117];
   armZZbar[111]=armZZbar[117] + armZZbar[116] + armZZbar[115] + 
   armZZbar[111] + 1./9.*armZZbar[114];
   armZZbar[111]=armZZbar[6]*armZZbar[111];
   armZZbar[114]=10351./96. - 685*armZZbar[31];
   armZZbar[110]=11*armZZbar[110] - 11./18.*armZZbar[13] + 187./36.*
   armZZbar[12] + 11*armZZbar[15] + 685./36.*armZZbar[177] - 253./24.*
   armZZbar[14] + 1./9.*armZZbar[114] - 11*armZZbar[28];
   armZZbar[113]=274./27.*armZZbar[25] + 11*armZZbar[113];
   armZZbar[113]=MMZ*armZZbar[113];
   armZZbar[114]=armZZbar[39]*armZZbar[175];
   armZZbar[109]=121./27.*armZZbar[114] + 22./3.*armZZbar[109] + 1./3.*
   armZZbar[110] + armZZbar[113];
   armZZbar[109]=armZZbar[39]*armZZbar[109];
   armZZbar[110]=armZZbar[132] + armZZbar[131] + armZZbar[127] - 3*
   armZZbar[8] + armZZbar[126];
   armZZbar[110]=armZZbar[1]*armZZbar[110];
   armZZbar[113]=armZZbar[17] + 9*MMZ;
   armZZbar[113]=armZZbar[1]*armZZbar[113];
   armZZbar[114]=armZZbar[39]*armZZbar[165];
   armZZbar[115]=9*armZZbar[1] + 11*armZZbar[39];
   armZZbar[115]=armZZbar[18]*armZZbar[115];
   armZZbar[113]=armZZbar[115] + armZZbar[113] + 11*armZZbar[114];
   armZZbar[113]=armZZbar[6]*armZZbar[113];
   armZZbar[114]= - 4./3.*MMZ + armZZbar[128] + armZZbar[21] - 
   armZZbar[8] + armZZbar[124];
   armZZbar[114]=armZZbar[39]*armZZbar[114];
   armZZbar[110]=1./2.*armZZbar[113] + armZZbar[190] + armZZbar[110] + 
   11./3.*armZZbar[114];
   armZZbar[110]=armZZbar[3]*armZZbar[110];
   armZZbar[113]=armZZbar[1]*armZZbar[12];
   armZZbar[114]=armZZbar[39]*armZZbar[12];
   armZZbar[113]=armZZbar[113] + 11./9.*armZZbar[114];
   armZZbar[114]= - armZZbar[1]*armZZbar[12];
   armZZbar[115]= - armZZbar[39]*armZZbar[12];
   armZZbar[114]=armZZbar[114] + 11./9.*armZZbar[115];
   armZZbar[114]=armZZbar[6]*armZZbar[114];
   armZZbar[113]=1./3.*armZZbar[113] + armZZbar[114];
   armZZbar[114]= - MMZ*armZZbar[12];
   armZZbar[114]=armZZbar[19] + armZZbar[114];
   armZZbar[115]= - 1./3. + armZZbar[199];
   armZZbar[115]=armZZbar[18]*armZZbar[115];
   armZZbar[114]=1./3.*armZZbar[114] + armZZbar[115];
   armZZbar[114]=armZZbar[3]*armZZbar[114];
   armZZbar[115]= - armZZbar[13]*armZZbar[12];
   armZZbar[115]=1./18.*armZZbar[115] + 1./27.*armZZbar[12] - 1./3.*
   armZZbar[16] - 1./4. + 1./3.*armZZbar[88];
   armZZbar[116]= - armZZbar[10]*armZZbar[12];
   armZZbar[114]=armZZbar[114] + 1./2.*armZZbar[115] + 1./3.*
   armZZbar[116];
   armZZbar[114]=armZZbar[73]*armZZbar[114];
   armZZbar[113]=1./3.*armZZbar[113] + armZZbar[114];
   armZZbar[113]=armZZbar[272]*armZZbar[113];
   armZZbar[114]=armZZbar[137] + armZZbar[166];
   armZZbar[114]=1./2.*armZZbar[114] + armZZbar[122];
   armZZbar[114]=MMH*armZZbar[73]*armZZbar[114];
   armZZbar[106]=1./4.*armZZbar[113] + armZZbar[114] + armZZbar[106] + 
   armZZbar[110] + armZZbar[111] + armZZbar[370] + armZZbar[107] + 1./3.
   *armZZbar[109];
   armZZbar[106]=armZZbar[272]*armZZbar[106];
   armZZbar[107]= - 713./36. + armZZbar[23];
   armZZbar[109]= - 289./27. - armZZbar[7];
   armZZbar[109]=armZZbar[13]*armZZbar[109];
   armZZbar[110]= - 548./9.*armZZbar[25] + 73./9.*armZZbar[27] + 
   armZZbar[26];
   armZZbar[110]=MMZ*armZZbar[110];
   armZZbar[111]= - 11*armZZbar[177];
   armZZbar[113]=19./9.*armZZbar[7] - 85./27. + armZZbar[111];
   armZZbar[113]=armZZbar[1]*armZZbar[113];
   armZZbar[111]= - 31./9. + armZZbar[111];
   armZZbar[111]=1./3.*armZZbar[111] + armZZbar[7];
   armZZbar[111]=armZZbar[39]*armZZbar[111];
   armZZbar[107]=11./9.*armZZbar[111] + 2./3.*armZZbar[113] + 1./3.*
   armZZbar[110] + armZZbar[292] + 1./4.*armZZbar[109] + 1./6.*
   armZZbar[7] - 95./18.*armZZbar[12] - 329./27.*armZZbar[177] + 1370./
   27.*armZZbar[31] + 1./6.*armZZbar[107] - armZZbar[32];
   armZZbar[107]=armZZbar[39]*armZZbar[107];
   armZZbar[109]=689./36. + 11*armZZbar[23];
   armZZbar[110]= - 89./3. - armZZbar[7];
   armZZbar[110]=armZZbar[13]*armZZbar[110];
   armZZbar[111]=armZZbar[195] + armZZbar[26];
   armZZbar[111]=1./3.*armZZbar[111] - 12*armZZbar[25];
   armZZbar[111]=MMZ*armZZbar[111];
   armZZbar[113]=1./9.*armZZbar[7] - 7./27. - armZZbar[177];
   armZZbar[113]=armZZbar[1]*armZZbar[113];
   armZZbar[109]=armZZbar[113] + armZZbar[111] + armZZbar[216] + 1./36.
   *armZZbar[110] + 1./54.*armZZbar[7] - 77./54.*armZZbar[12] - 7*
   armZZbar[177] + 30*armZZbar[31] + 1./18.*armZZbar[109] - 
   armZZbar[32];
   armZZbar[109]=armZZbar[1]*armZZbar[109];
   armZZbar[110]=armZZbar[310] - MMZ;
   armZZbar[110]=armZZbar[1]*armZZbar[110];
   armZZbar[111]=armZZbar[19] - 1./3.*MMZ;
   armZZbar[111]=armZZbar[39]*armZZbar[111];
   armZZbar[110]=armZZbar[125] + armZZbar[110] + 11./3.*armZZbar[111];
   armZZbar[110]=armZZbar[6]*armZZbar[110];
   armZZbar[111]=20./27. - armZZbar[7];
   armZZbar[111]=MMZ*armZZbar[111];
   armZZbar[111]=armZZbar[146] + armZZbar[111];
   armZZbar[111]=armZZbar[39]*armZZbar[111];
   armZZbar[113]=4./3. - armZZbar[7];
   armZZbar[113]=MMZ*armZZbar[113];
   armZZbar[113]= - armZZbar[19] + 1./3.*armZZbar[113];
   armZZbar[113]=armZZbar[1]*armZZbar[113];
   armZZbar[110]=armZZbar[110] + armZZbar[139] + armZZbar[113] + 
   armZZbar[111];
   armZZbar[110]=armZZbar[3]*armZZbar[110];
   armZZbar[111]=151./3. + armZZbar[119];
   armZZbar[111]=armZZbar[1]*armZZbar[111];
   armZZbar[113]=53./27. - armZZbar[7];
   armZZbar[113]=armZZbar[39]*armZZbar[113];
   armZZbar[111]=11*armZZbar[113] + 2./3.*armZZbar[111] - 11*
   armZZbar[11] + 88./3.*armZZbar[13] + 784./9. + 95./2.*armZZbar[12];
   armZZbar[111]=armZZbar[39]*armZZbar[111];
   armZZbar[113]=7 - 2*armZZbar[1];
   armZZbar[113]=armZZbar[1]*armZZbar[113];
   armZZbar[114]= - 242./9.*armZZbar[39] + 329./9. - 44*armZZbar[1];
   armZZbar[114]=armZZbar[39]*armZZbar[114];
   armZZbar[113]=armZZbar[113] + 1./9.*armZZbar[114];
   armZZbar[113]=armZZbar[6]*armZZbar[113];
   armZZbar[114]=13./3. - armZZbar[7];
   armZZbar[114]=armZZbar[1]*armZZbar[114];
   armZZbar[114]=1./3.*armZZbar[114] - armZZbar[11] + 8./3.*
   armZZbar[13] + 10 + 77./18.*armZZbar[12];
   armZZbar[114]=armZZbar[1]*armZZbar[114];
   armZZbar[111]=armZZbar[113] + armZZbar[114] + 1./9.*armZZbar[111];
   armZZbar[111]=armZZbar[6]*armZZbar[111];
   armZZbar[113]=armZZbar[150] + armZZbar[168] + armZZbar[137] + 1./2.*
   armZZbar[80];
   armZZbar[113]=MMH*armZZbar[73]*armZZbar[113];
   armZZbar[105]=armZZbar[106] + armZZbar[113] + armZZbar[105] + 
   armZZbar[110] + armZZbar[111] + armZZbar[325] + armZZbar[109] + 1./3.
   *armZZbar[107];
   armZZbar[105]=armZZbar[272]*armZZbar[105];
   armZZbar[106]=armZZbar[16] - 2*armZZbar[89] + armZZbar[90] - 2 - 
   armZZbar[97];
   armZZbar[106]=armZZbar[108]*armZZbar[106];
   armZZbar[107]=82*armZZbar[90] - 29 - 82*armZZbar[97];
   armZZbar[106]=8*armZZbar[106] + 82./3.*armZZbar[16] - 164./3.*
   armZZbar[89] + 1./3.*armZZbar[107] + armZZbar[121];
   armZZbar[106]=armZZbar[100]*armZZbar[106];
   armZZbar[107]=armZZbar[82] - 2*armZZbar[95] - armZZbar[83];
   armZZbar[107]=armZZbar[108]*armZZbar[107];
   armZZbar[109]=46*armZZbar[82] - 58*armZZbar[95] - 23*armZZbar[83];
   armZZbar[107]=4*armZZbar[107] + 1./3.*armZZbar[109] - armZZbar[81];
   armZZbar[107]=armZZbar[108]*armZZbar[107];
   armZZbar[109]= - 100./3.*armZZbar[95] - 7*armZZbar[83];
   armZZbar[110]=41./3. + armZZbar[112];
   armZZbar[110]=armZZbar[13]*armZZbar[100]*armZZbar[110];
   armZZbar[111]= - 1 - armZZbar[90];
   armZZbar[111]=armZZbar[101]*armZZbar[111];
   armZZbar[106]=8*armZZbar[111] + 2*armZZbar[110] + 2*armZZbar[205] + 
   armZZbar[106] + 4*armZZbar[107] - 41./3.*armZZbar[81] - 16*
   armZZbar[80] + 2*armZZbar[109] + 65*armZZbar[82];
   armZZbar[106]=MMZ*armZZbar[106];
   armZZbar[107]=149./9. + armZZbar[112];
   armZZbar[107]=armZZbar[108]*armZZbar[107];
   armZZbar[109]=12*armZZbar[108];
   armZZbar[110]=119./3. + armZZbar[109];
   armZZbar[110]=armZZbar[12]*armZZbar[110];
   armZZbar[107]=armZZbar[110] + 1495./9. + 4*armZZbar[107];
   armZZbar[110]=44./9. + armZZbar[108];
   armZZbar[110]=armZZbar[108]*armZZbar[110];
   armZZbar[110]=1813./9. + 16*armZZbar[110];
   armZZbar[110]=armZZbar[13]*armZZbar[110];
   armZZbar[107]=2*armZZbar[107] + armZZbar[110];
   armZZbar[107]=armZZbar[13]*armZZbar[107];
   armZZbar[110]= - 47./3. + armZZbar[129];
   armZZbar[110]=armZZbar[13]*armZZbar[110];
   armZZbar[110]=armZZbar[110] + armZZbar[218] + 1./3. + armZZbar[129];
   armZZbar[110]=MMZ*armZZbar[110];
   armZZbar[111]= - 5 - 3*armZZbar[13];
   armZZbar[111]=armZZbar[19]*armZZbar[111];
   armZZbar[112]=1 - armZZbar[13];
   armZZbar[112]=armZZbar[17]*armZZbar[112];
   armZZbar[110]=armZZbar[110] + armZZbar[112] - armZZbar[20] + 2*
   armZZbar[111];
   armZZbar[111]=2 + armZZbar[213];
   armZZbar[111]=armZZbar[3]*armZZbar[118]*armZZbar[111];
   armZZbar[110]=2*armZZbar[110] + armZZbar[111];
   armZZbar[110]=armZZbar[3]*armZZbar[110];
   armZZbar[111]=109./9. - armZZbar[88];
   armZZbar[111]=6*armZZbar[108] + armZZbar[16] + 2*armZZbar[111] + 
   armZZbar[90];
   armZZbar[111]=armZZbar[108]*armZZbar[111];
   armZZbar[112]=1117 - 244*armZZbar[88];
   armZZbar[113]= - armZZbar[76] + 4*armZZbar[22];
   armZZbar[113]=armZZbar[100]*armZZbar[113];
   armZZbar[109]=71./3. + armZZbar[109];
   armZZbar[109]=armZZbar[12]*armZZbar[109];
   armZZbar[114]=armZZbar[11]*armZZbar[262];
   armZZbar[115]=armZZbar[19]*armZZbar[219];
   armZZbar[116]=2*armZZbar[100] + armZZbar[147];
   armZZbar[116]=armZZbar[17]*armZZbar[116];
   armZZbar[117]=armZZbar[10]*MMZ*armZZbar[100];
   armZZbar[106]=2*armZZbar[110] + 2*armZZbar[117] + armZZbar[106] + 4*
   armZZbar[116] + 8*armZZbar[115] + 6*armZZbar[114] + armZZbar[107] + 
   2*armZZbar[109] + 4*armZZbar[113] + 8*armZZbar[111] + 58*
   armZZbar[16] - 8*armZZbar[89] + 190./9.*armZZbar[90] + 1./3.*
   armZZbar[112] + 2*armZZbar[97];
   armZZbar[106]=armZZbar[73]*armZZbar[106];
   armZZbar[107]= - 11*armZZbar[32] - 19./9. - 11*armZZbar[23];
   armZZbar[107]=armZZbar[108]*armZZbar[107];
   armZZbar[109]=8./27.*armZZbar[256] + armZZbar[135];
   armZZbar[109]=armZZbar[13]*armZZbar[109];
   armZZbar[110]= - armZZbar[24] - 2./9.*armZZbar[26];
   armZZbar[110]=armZZbar[108]*armZZbar[110];
   armZZbar[111]= - 2*armZZbar[24];
   armZZbar[112]=armZZbar[111] - 17./81.*armZZbar[27];
   armZZbar[110]=2*armZZbar[110] + 272./81.*armZZbar[25] + 4*
   armZZbar[112] - 17./9.*armZZbar[26];
   armZZbar[110]=MMZ*armZZbar[110];
   armZZbar[112]=5*armZZbar[29];
   armZZbar[113]= - 5*armZZbar[16];
   armZZbar[107]=200./729.*armZZbar[39] + 80./243.*armZZbar[1] + 
   armZZbar[110] + armZZbar[109] + armZZbar[135] + 1./3.*armZZbar[107]
    + 134./81.*armZZbar[177] + armZZbar[113] - 680./81.*armZZbar[31] - 
   101./9.*armZZbar[32] - 55./9.*armZZbar[23] - 511./81. + 
   armZZbar[112];
   armZZbar[107]=armZZbar[39]*armZZbar[107];
   armZZbar[109]= - 5*armZZbar[32] + 1./3. - 5*armZZbar[23];
   armZZbar[109]=armZZbar[108]*armZZbar[109];
   armZZbar[109]=armZZbar[109] + 10*armZZbar[177] + armZZbar[113] - 40*
   armZZbar[31] - 13*armZZbar[32] - 17./3.*armZZbar[23] - 127./9. + 
   armZZbar[112];
   armZZbar[110]=armZZbar[111] - armZZbar[27];
   armZZbar[111]= - armZZbar[108]*armZZbar[24];
   armZZbar[110]=2*armZZbar[111] + 16*armZZbar[25] + 4*armZZbar[110] - 
   armZZbar[26];
   armZZbar[110]=MMZ*armZZbar[110];
   armZZbar[111]=7./3. + 2*armZZbar[108];
   armZZbar[111]=8./9.*armZZbar[111] - armZZbar[7];
   armZZbar[111]=armZZbar[13]*armZZbar[111];
   armZZbar[109]=8./81.*armZZbar[1] + 1./3.*armZZbar[110] + 
   armZZbar[111] + 1./3.*armZZbar[109] - armZZbar[7];
   armZZbar[109]=armZZbar[1]*armZZbar[109];
   armZZbar[110]= - 2./3. - armZZbar[108];
   armZZbar[111]= - 49./3. + armZZbar[153];
   armZZbar[111]=armZZbar[13]*armZZbar[111];
   armZZbar[110]= - 8./9.*armZZbar[1] + 8*armZZbar[110] + armZZbar[111]
   ;
   armZZbar[110]=armZZbar[1]*armZZbar[110];
   armZZbar[108]= - 34./9. - 5*armZZbar[108];
   armZZbar[111]= - 61./3. + armZZbar[153];
   armZZbar[111]=armZZbar[13]*armZZbar[111];
   armZZbar[108]= - 200./27.*armZZbar[39] - 80./9.*armZZbar[1] + 8*
   armZZbar[108] + 5*armZZbar[111];
   armZZbar[108]=armZZbar[39]*armZZbar[108];
   armZZbar[111]= - 5 + armZZbar[158];
   armZZbar[111]=armZZbar[1]*armZZbar[111];
   armZZbar[112]=100./3.*armZZbar[39] - 67./3. + 40*armZZbar[1];
   armZZbar[112]=armZZbar[39]*armZZbar[112];
   armZZbar[111]=armZZbar[111] + 1./9.*armZZbar[112];
   armZZbar[111]=armZZbar[6]*armZZbar[111];
   armZZbar[108]=armZZbar[111] + armZZbar[110] + 1./3.*armZZbar[108];
   armZZbar[108]=armZZbar[6]*armZZbar[108];
   armZZbar[110]= - armZZbar[1]*MMZ;
   armZZbar[111]= - armZZbar[39]*MMZ;
   armZZbar[110]=armZZbar[110] + 5./3.*armZZbar[111];
   armZZbar[111]=armZZbar[1]*MMZ;
   armZZbar[112]=armZZbar[39]*MMZ;
   armZZbar[111]=armZZbar[111] + 5./3.*armZZbar[112];
   armZZbar[111]=armZZbar[6]*armZZbar[111];
   armZZbar[110]=1./3.*armZZbar[110] + armZZbar[111];
   armZZbar[110]=armZZbar[3]*armZZbar[110];
   armZZbar[107]=8./3.*armZZbar[110] + 2./3.*armZZbar[108] + 
   armZZbar[109] + armZZbar[107];
   armZZbar[108]=MMH*armZZbar[73]*armZZbar[81];

      mZZbarret = armZZbar[102] + armZZbar[103] + armZZbar[104] + 
      armZZbar[105] + armZZbar[106] + 2*armZZbar[107] + 10./3.*
      armZZbar[108];
      return mZZbarret;
}
