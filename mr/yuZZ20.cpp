#include <ZZ.hpp>
std::complex<long double>
ZZ<OS>::my20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuZZ[391], yuZZret;

    aryuZZ[1]=double(nL + nH);
    aryuZZ[2]=pow(CW,-1);
    aryuZZ[3]=pow(MMH,-1);
    aryuZZ[4]=pow(MMZ,-1);
    aryuZZ[5]=pow(SW,-1);
    aryuZZ[6]=std::real(Tsil::B(0,0,MMZ,mu2));
    aryuZZ[7]=std::real(Tsil::B(0,0,MMW,mu2));
    aryuZZ[8]=Tsil::I2(0,0,MMZ,mu2);
    aryuZZ[9]=Tsil::I2(0,0,MMW,mu2);
    aryuZZ[10]=Tsil::B(MMZ,MMH,MMZ,mu2);
    aryuZZ[11]=Tsil::B(MMW,MMH,MMW,mu2);
    aryuZZ[12]=Tsil::B(MMW,MMZ,MMW,mu2);
    aryuZZ[13]=Tsil::B(MMW,MMW,MMZ,mu2);
    aryuZZ[14]=Tsil::B(0,MMH,MMZ,mu2);
    aryuZZ[15]=Tsil::Beps(MMZ,MMH,MMZ,mu2);
    aryuZZ[16]=Tsil::Beps(MMW,MMW,MMZ,mu2);
    aryuZZ[17]=Tsil::A(MMH,mu2);
    aryuZZ[18]=Tsil::A(MMZ,mu2);
    aryuZZ[19]=Tsil::A(MMW,mu2);
    aryuZZ[20]=Tsil::Aeps(MMH,mu2);
    aryuZZ[21]=Tsil::Aeps(MMZ,mu2);
    aryuZZ[22]=Tsil::Aeps(MMW,mu2);
    aryuZZ[23]=std::real(Tsil::B(0,MMW,MMZ,mu2));
    aryuZZ[24]=protW0W00->M(0);
    aryuZZ[25]=prot0000Z->M(0);
    aryuZZ[26]=prot0000W->M(0);
    aryuZZ[27]=prot00000->M(0);
    aryuZZ[28]=protHZ00->Uxzuv(0);
    aryuZZ[29]=protW0W00->Uzxyv(0);
    aryuZZ[30]=protHZ00->Txuv(0);
    aryuZZ[31]=prot0000Z->Tvxu(0);
    aryuZZ[32]=protW0W00->Txuv(0);
    aryuZZ[33]=double(nH);
    aryuZZ[34]=pow(MMt,-1);
    aryuZZ[35]=Tsil::B(MMt,MMt,MMZ,mu2);
    aryuZZ[36]=Tsil::B(0,MMt,MMW,mu2);
    aryuZZ[37]=Tsil::A(MMt,mu2);
    aryuZZ[38]=Tsil::B(MMt,MMt,MMH,mu2);
    aryuZZ[39]=double(nL);
    aryuZZ[40]=Tsil::I2(MMH,MMt,MMt,mu2);
    aryuZZ[41]=Tsil::I2(MMZ,MMt,MMt,mu2);
    aryuZZ[42]=Tsil::I2(0,MMW,MMt,mu2);
    aryuZZ[43]=Tsil::B(MMH,MMH,MMH,mu2);
    aryuZZ[44]=Tsil::B(MMH,MMt,MMt,mu2);
    aryuZZ[45]=Tsil::B(MMZ,MMZ,MMH,mu2);
    aryuZZ[46]=Tsil::B(MMZ,MMt,MMt,mu2);
    aryuZZ[47]=Tsil::B(MMW,MMW,MMH,mu2);
    aryuZZ[48]=std::real(Tsil::B(0,MMW,MMt,mu2));
    aryuZZ[49]=Tsil::Beps(MMt,MMt,MMZ,mu2);
    aryuZZ[50]=Tsil::Aeps(MMt,mu2);
    aryuZZ[51]=protWtWt0->M(0);
    aryuZZ[52]=protW0W0t->M(0);
    aryuZZ[53]=prottZtHt->M(0);
    aryuZZ[54]=protttttH->M(0);
    aryuZZ[55]=protttttZ->M(0);
    aryuZZ[56]=prottttt0->M(0);
    aryuZZ[57]=prott0t0W->M(0);
    aryuZZ[58]=prottttt0->Vzxyv(0);
    aryuZZ[59]=prottZtHt->Uuyxv(0);
    aryuZZ[60]=protWtWt0->Uzxyv(0);
    aryuZZ[61]=protttttH->Uzxyv(0);
    aryuZZ[62]=protttttZ->Uzxyv(0);
    aryuZZ[63]=protWtWt0->Uuyxv(0);
    aryuZZ[64]=prottZtHt->Tuxv(0);
    aryuZZ[65]=protWtWt0->Tzyv(0);
    aryuZZ[66]=protWtWt0->Tyzv(0);
    aryuZZ[67]=prottZtHt->Txuv(0);
    aryuZZ[68]=prottZtHt->Suxv(0);
    aryuZZ[69]=prottZtHt->Svyz(0);
    aryuZZ[70]=protWtWt0->Svyz(0);
    aryuZZ[71]=prottttt0->Suxv(0);
    aryuZZ[72]=prottZtHt->Uyuzv(0);
    aryuZZ[73]=double(boson);
    aryuZZ[74]=Tsil::I2(MMH,MMH,MMH,mu2);
    aryuZZ[75]=Tsil::I2(MMH,MMZ,MMZ,mu2);
    aryuZZ[76]=Tsil::I2(MMH,MMW,MMW,mu2);
    aryuZZ[77]=Tsil::I2(MMW,MMW,MMZ,mu2);
    aryuZZ[78]=protZHHZZ->M(0);
    aryuZZ[79]=protZZHHH->M(0);
    aryuZZ[80]=protZWHWW->M(0);
    aryuZZ[81]=protWWWWH->M(0);
    aryuZZ[82]=protWWWWZ->M(0);
    aryuZZ[83]=protWWWW0->M(0);
    aryuZZ[84]=protZHHZZ->Uzxyv(0);
    aryuZZ[85]=protZWHWW->Uzxyv(0);
    aryuZZ[86]=protZZHHH->Uyuzv(0);
    aryuZZ[87]=protZHHZZ->Uuyxv(0);
    aryuZZ[88]=protZWHWW->Uuyxv(0);
    aryuZZ[89]=protZWHWW->Tzyv(0);
    aryuZZ[90]=protZWHWW->Tyzv(0);
    aryuZZ[91]=protZZHHH->Suxv(0);
    aryuZZ[92]=protZWHWW->Svyz(0);
    aryuZZ[93]=protWZWHW->Svyz(0);
    aryuZZ[94]=protWWWW0->Suxv(0);
    aryuZZ[95]=protWWWW0->Vzxyv(0);
    aryuZZ[96]=protZWHWW->Uxzuv(0);
    aryuZZ[97]=protZWHWW->Uyuzv(0);
    aryuZZ[98]=1/(MMt - MMW);
    aryuZZ[99]=1/(4*MMt - MMZ);
    aryuZZ[100]=1/( - 4*MMW + MMH);
    aryuZZ[101]=1/( - 4*MMZ + MMH);
   aryuZZ[102]= - 161./27. + 5*aryuZZ[12];
   aryuZZ[103]= - 2*aryuZZ[13];
   aryuZZ[102]=1./8.*aryuZZ[102] + aryuZZ[103];
   aryuZZ[104]=1./9.*aryuZZ[1];
   aryuZZ[105]=1./3.*aryuZZ[39];
   aryuZZ[106]=aryuZZ[105] - 1./8. + aryuZZ[104];
   aryuZZ[106]=aryuZZ[6]*aryuZZ[106];
   aryuZZ[107]=1./4.*aryuZZ[19]*aryuZZ[98];
   aryuZZ[108]=pow(aryuZZ[98],2);
   aryuZZ[109]= - aryuZZ[19]*aryuZZ[108];
   aryuZZ[110]= - aryuZZ[98] + aryuZZ[109];
   aryuZZ[111]=3./2.*aryuZZ[110] + aryuZZ[101];
   aryuZZ[111]=MMZ*aryuZZ[111];
   aryuZZ[112]= - 25./9. - aryuZZ[7];
   aryuZZ[112]=7./3.*aryuZZ[35] + 1./2.*aryuZZ[112] - aryuZZ[36];
   aryuZZ[113]=aryuZZ[1]*aryuZZ[112];
   aryuZZ[112]=aryuZZ[39]*aryuZZ[112];
   aryuZZ[114]= - MMZ*aryuZZ[101];
   aryuZZ[115]=7./3. + aryuZZ[114];
   aryuZZ[115]=aryuZZ[10]*aryuZZ[115];
   aryuZZ[116]= - aryuZZ[101]*aryuZZ[17];
   aryuZZ[117]=aryuZZ[18]*aryuZZ[101];
   aryuZZ[118]= - 1./2.*aryuZZ[11];
   aryuZZ[102]=5*aryuZZ[106] + 1./2.*aryuZZ[117] + 1./2.*aryuZZ[115] + 
   aryuZZ[112] + 1./3.*aryuZZ[113] + 5./2.*aryuZZ[35] + 1./2.*
   aryuZZ[111] + 1./4.*aryuZZ[116] + aryuZZ[107] + 11*aryuZZ[102] + 
   aryuZZ[118];
   aryuZZ[102]=aryuZZ[6]*aryuZZ[102];
   aryuZZ[106]=2*aryuZZ[22];
   aryuZZ[111]= - 1./2.*aryuZZ[8];
   aryuZZ[112]= - 1./2.*aryuZZ[41];
   aryuZZ[113]=aryuZZ[106] + 3*aryuZZ[50] + aryuZZ[112] - 2*aryuZZ[42]
    + aryuZZ[111];
   aryuZZ[115]=11*aryuZZ[19];
   aryuZZ[119]=1./2.*aryuZZ[17];
   aryuZZ[120]=31./2.*MMZ + aryuZZ[115] + aryuZZ[119];
   aryuZZ[120]=aryuZZ[35]*aryuZZ[120];
   aryuZZ[121]= - 3./2.*aryuZZ[36];
   aryuZZ[122]=7./2.*aryuZZ[35];
   aryuZZ[123]=aryuZZ[122] + 1 + aryuZZ[121];
   aryuZZ[123]=aryuZZ[18]*aryuZZ[123];
   aryuZZ[124]=7*aryuZZ[19];
   aryuZZ[125]=5*aryuZZ[18] + 23./2.*MMZ + aryuZZ[124] + aryuZZ[119];
   aryuZZ[125]=aryuZZ[6]*aryuZZ[125];
   aryuZZ[126]=1./2.*aryuZZ[20];
   aryuZZ[127]=3*aryuZZ[21];
   aryuZZ[128]= - 12*MMZ;
   aryuZZ[129]=3*aryuZZ[46];
   aryuZZ[130]=6*aryuZZ[48] + 37./2. + aryuZZ[129];
   aryuZZ[130]=aryuZZ[37]*aryuZZ[130];
   aryuZZ[131]= - 1./2.*aryuZZ[17];
   aryuZZ[132]=2*aryuZZ[19];
   aryuZZ[113]=aryuZZ[130] + 1./2.*aryuZZ[125] + aryuZZ[123] + 1./2.*
   aryuZZ[120] + aryuZZ[128] + aryuZZ[131] + aryuZZ[132] + aryuZZ[127]
    + 3*aryuZZ[113] + aryuZZ[126];
   aryuZZ[113]=aryuZZ[3]*aryuZZ[113];
   aryuZZ[120]= - aryuZZ[19] - MMZ;
   aryuZZ[123]= - 1./2.*aryuZZ[18];
   aryuZZ[125]=aryuZZ[120] + aryuZZ[123];
   aryuZZ[125]=aryuZZ[6]*aryuZZ[125];
   aryuZZ[130]=MMZ*aryuZZ[36];
   aryuZZ[133]=aryuZZ[35]*aryuZZ[120];
   aryuZZ[134]= - 1./2.*aryuZZ[35];
   aryuZZ[135]=aryuZZ[36] + aryuZZ[134];
   aryuZZ[136]=aryuZZ[18]*aryuZZ[135];
   aryuZZ[137]=aryuZZ[19]*aryuZZ[36];
   aryuZZ[125]=1./2.*aryuZZ[125] + 1./2.*aryuZZ[136] + 1./2.*
   aryuZZ[133] + aryuZZ[137] + aryuZZ[130];
   aryuZZ[125]=aryuZZ[3]*aryuZZ[125];
   aryuZZ[130]=33*aryuZZ[13] + 59./9. - 33./2.*aryuZZ[12];
   aryuZZ[130]=1./2.*aryuZZ[130] + aryuZZ[11];
   aryuZZ[133]=1 + aryuZZ[7];
   aryuZZ[133]= - aryuZZ[35] + 1./2.*aryuZZ[133] + aryuZZ[36];
   aryuZZ[136]=aryuZZ[1]*aryuZZ[133];
   aryuZZ[133]=aryuZZ[39]*aryuZZ[133];
   aryuZZ[138]= - 1./3.*aryuZZ[1];
   aryuZZ[139]=aryuZZ[138] - aryuZZ[39];
   aryuZZ[140]=aryuZZ[6]*aryuZZ[139];
   aryuZZ[133]=aryuZZ[140] - aryuZZ[10] + aryuZZ[133] + 1./2.*
   aryuZZ[130] + 1./3.*aryuZZ[136];
   aryuZZ[133]=aryuZZ[6]*aryuZZ[133];
   aryuZZ[136]= - aryuZZ[7] - aryuZZ[36];
   aryuZZ[141]=1./3. + aryuZZ[7];
   aryuZZ[142]=aryuZZ[35]*aryuZZ[141];
   aryuZZ[136]=1./3.*aryuZZ[136] + 1./2.*aryuZZ[142];
   aryuZZ[142]=aryuZZ[1]*aryuZZ[136];
   aryuZZ[136]=aryuZZ[39]*aryuZZ[136];
   aryuZZ[143]= - 59./9. - 33./2.*aryuZZ[13];
   aryuZZ[143]=aryuZZ[36]*aryuZZ[143];
   aryuZZ[144]=aryuZZ[35]*aryuZZ[130];
   aryuZZ[145]= - aryuZZ[36] + aryuZZ[35];
   aryuZZ[146]=MMt*aryuZZ[3]*aryuZZ[145];
   aryuZZ[147]=aryuZZ[12] - aryuZZ[13];
   aryuZZ[148]=1./3. + aryuZZ[36];
   aryuZZ[149]=aryuZZ[148] - aryuZZ[35];
   aryuZZ[150]=aryuZZ[10]*aryuZZ[149];
   aryuZZ[151]= - 1./3.*aryuZZ[11];
   aryuZZ[125]=3./2.*aryuZZ[146] + 3*aryuZZ[125] + aryuZZ[133] + 
   aryuZZ[150] + aryuZZ[136] + 1./3.*aryuZZ[142] + 1./2.*aryuZZ[144] + 
   1./2.*aryuZZ[143] + 11./4.*aryuZZ[147] + aryuZZ[151];
   aryuZZ[133]=pow(aryuZZ[5],2);
   aryuZZ[125]=aryuZZ[133]*aryuZZ[125];
   aryuZZ[136]= - 7 + 27./4.*aryuZZ[63];
   aryuZZ[142]= - 1./2.*aryuZZ[61];
   aryuZZ[136]= - 47./8.*aryuZZ[49] - 81./16.*aryuZZ[66] + aryuZZ[142]
    - 27./4.*aryuZZ[65] + 1./2.*aryuZZ[136] + 3*aryuZZ[62];
   aryuZZ[143]= - 1./4.*aryuZZ[44];
   aryuZZ[136]= - aryuZZ[15] + aryuZZ[143] + 27./16.*aryuZZ[48] + 3./4.
   *aryuZZ[46] + 3./2.*aryuZZ[67] + 1./2.*aryuZZ[136] + aryuZZ[59];
   aryuZZ[136]=aryuZZ[99]*aryuZZ[136];
   aryuZZ[144]= - aryuZZ[65] + aryuZZ[66];
   aryuZZ[144]=aryuZZ[108]*aryuZZ[144];
   aryuZZ[146]=MMZ*aryuZZ[144];
   aryuZZ[152]=1./2.*aryuZZ[66] + 33./16. + aryuZZ[65];
   aryuZZ[152]=aryuZZ[98]*aryuZZ[152];
   aryuZZ[153]=aryuZZ[50]*aryuZZ[108];
   aryuZZ[154]= - aryuZZ[22]*aryuZZ[108];
   aryuZZ[155]=3./2.*aryuZZ[49] - 1./2.*aryuZZ[62] + 9./8. - aryuZZ[63]
   ;
   aryuZZ[155]=aryuZZ[34]*aryuZZ[155];
   aryuZZ[156]=aryuZZ[19]*aryuZZ[108];
   aryuZZ[157]=1./2.*aryuZZ[59];
   aryuZZ[158]= - aryuZZ[15] + 3*aryuZZ[64] + aryuZZ[157] - 3./2.*
   aryuZZ[14] + 1./2.*aryuZZ[28] + 3./8. - aryuZZ[30];
   aryuZZ[158]=aryuZZ[101]*aryuZZ[158];
   aryuZZ[136]=3./4.*aryuZZ[146] + aryuZZ[158] + 3./4.*aryuZZ[156] + 3./
   16.*aryuZZ[136] + 1./2.*aryuZZ[155] + 3./4.*aryuZZ[154] + 3./4.*
   aryuZZ[153] + aryuZZ[152] - 4*aryuZZ[57] + aryuZZ[25] - 6*aryuZZ[51]
    + aryuZZ[55] - 6*aryuZZ[52];
   aryuZZ[136]=MMZ*aryuZZ[136];
   aryuZZ[146]= - 3*aryuZZ[72];
   aryuZZ[152]= - 1./2.*aryuZZ[59];
   aryuZZ[155]= - 3*aryuZZ[64];
   aryuZZ[158]= - 3*aryuZZ[38];
   aryuZZ[159]=aryuZZ[158] + 7./2.*aryuZZ[15] + aryuZZ[67] + 
   aryuZZ[155] + aryuZZ[152] + 1./4. + aryuZZ[146];
   aryuZZ[159]=aryuZZ[101]*aryuZZ[159];
   aryuZZ[160]=47./2. + aryuZZ[129];
   aryuZZ[161]= - 1 + aryuZZ[47];
   aryuZZ[161]=2*aryuZZ[161] + aryuZZ[45];
   aryuZZ[161]=aryuZZ[3]*aryuZZ[37]*aryuZZ[161];
   aryuZZ[162]=17./4.*aryuZZ[35];
   aryuZZ[160]=9*aryuZZ[161] + aryuZZ[162] + aryuZZ[121] - 9./2.*
   aryuZZ[38] + 1./2.*aryuZZ[160] + 3*aryuZZ[48];
   aryuZZ[160]=aryuZZ[3]*aryuZZ[160];
   aryuZZ[161]=aryuZZ[101]*aryuZZ[38];
   aryuZZ[163]=3*aryuZZ[161];
   aryuZZ[164]=aryuZZ[35]*aryuZZ[101];
   aryuZZ[165]=aryuZZ[163] + 1./2.*aryuZZ[164];
   aryuZZ[165]=aryuZZ[10]*aryuZZ[165];
   aryuZZ[166]=pow(aryuZZ[3],2);
   aryuZZ[167]=MMt*aryuZZ[166]*aryuZZ[38];
   aryuZZ[168]= - aryuZZ[35]*aryuZZ[101];
   aryuZZ[169]=1./2.*aryuZZ[168];
   aryuZZ[170]=2*aryuZZ[53];
   aryuZZ[159]=18*aryuZZ[167] + aryuZZ[160] + aryuZZ[165] + aryuZZ[169]
    + aryuZZ[159] + aryuZZ[170] + 4*aryuZZ[57] + 7*aryuZZ[51] - 2*
   aryuZZ[55] + 7./2.*aryuZZ[52];
   aryuZZ[159]=MMt*aryuZZ[159];
   aryuZZ[160]=3*aryuZZ[44];
   aryuZZ[165]=aryuZZ[160] - 81./4.*aryuZZ[48] + 151./2. - 9*aryuZZ[46]
   ;
   aryuZZ[165]=aryuZZ[99]*aryuZZ[165];
   aryuZZ[171]= - MMZ*aryuZZ[108];
   aryuZZ[172]=aryuZZ[34] - 9./4.*aryuZZ[99];
   aryuZZ[172]=aryuZZ[35]*aryuZZ[172];
   aryuZZ[173]=MMZ*aryuZZ[108];
   aryuZZ[174]= - aryuZZ[98] + 3./4.*aryuZZ[173];
   aryuZZ[174]=aryuZZ[6]*aryuZZ[174];
   aryuZZ[175]=aryuZZ[10]*aryuZZ[99];
   aryuZZ[165]=aryuZZ[174] + 3./4.*aryuZZ[175] + 1./4.*aryuZZ[172] + 3./
   4.*aryuZZ[171] - 2*aryuZZ[101] + 1./16.*aryuZZ[165] - 5./4.*
   aryuZZ[98] - 3*aryuZZ[34];
   aryuZZ[165]=aryuZZ[37]*aryuZZ[165];
   aryuZZ[172]= - 5*aryuZZ[31];
   aryuZZ[174]= - 23./24.*aryuZZ[14];
   aryuZZ[176]= - 5./8.*aryuZZ[59];
   aryuZZ[177]=pow(Pi,2);
   aryuZZ[178]=aryuZZ[176] + 5./4.*aryuZZ[177] - 33*aryuZZ[16] - 1037./
   128.*aryuZZ[49] + 141./256.*aryuZZ[66] + aryuZZ[174] - 3./32.*
   aryuZZ[61] + 207./64.*aryuZZ[65] + 49./16.*aryuZZ[62] - aryuZZ[28]
    + 657./128.*aryuZZ[63] + aryuZZ[172] + 1295./216. + 33*aryuZZ[60];
   aryuZZ[179]=aryuZZ[160] + 111./4.*aryuZZ[48] - 16243./27. + 15*
   aryuZZ[46];
   aryuZZ[180]=55*aryuZZ[12];
   aryuZZ[179]=1./8.*aryuZZ[179] + aryuZZ[180];
   aryuZZ[181]= - 1 - aryuZZ[46];
   aryuZZ[182]= - 9./4.*aryuZZ[48];
   aryuZZ[183]=aryuZZ[181] + aryuZZ[182];
   aryuZZ[183]=3*aryuZZ[183] + aryuZZ[44];
   aryuZZ[183]=aryuZZ[99]*aryuZZ[183];
   aryuZZ[183]=3./32.*aryuZZ[183] + aryuZZ[101];
   aryuZZ[183]=MMZ*aryuZZ[183];
   aryuZZ[184]= - aryuZZ[19]*aryuZZ[34];
   aryuZZ[185]=1./2.*aryuZZ[116];
   aryuZZ[179]= - 5./4.*aryuZZ[35] + aryuZZ[183] + aryuZZ[185] + 
   aryuZZ[184] - aryuZZ[11] + 1./4.*aryuZZ[179] - 55*aryuZZ[13];
   aryuZZ[179]=aryuZZ[35]*aryuZZ[179];
   aryuZZ[183]= - 11./9. - aryuZZ[7];
   aryuZZ[186]=aryuZZ[35]*aryuZZ[183];
   aryuZZ[187]=1./3.*aryuZZ[7];
   aryuZZ[188]=1./3.*aryuZZ[36];
   aryuZZ[186]=1./2.*aryuZZ[186] + aryuZZ[188] + aryuZZ[187] + 2./9. + 
   aryuZZ[177];
   aryuZZ[189]=aryuZZ[1]*aryuZZ[186];
   aryuZZ[186]=aryuZZ[39]*aryuZZ[186];
   aryuZZ[190]= - 1./4.*aryuZZ[41];
   aryuZZ[191]= - 1./2.*aryuZZ[22];
   aryuZZ[192]=aryuZZ[70] - 1./2.*aryuZZ[42];
   aryuZZ[193]=aryuZZ[191] + aryuZZ[190] + aryuZZ[192] + 1./2.*
   aryuZZ[69];
   aryuZZ[193]=aryuZZ[34]*aryuZZ[193];
   aryuZZ[194]= - aryuZZ[70] + 1./2.*aryuZZ[42];
   aryuZZ[194]=9./8.*aryuZZ[194] - aryuZZ[69];
   aryuZZ[195]=1./2.*aryuZZ[68];
   aryuZZ[196]= - 1./4.*aryuZZ[40];
   aryuZZ[194]=aryuZZ[196] + 3*aryuZZ[194] + aryuZZ[195];
   aryuZZ[197]=5./4.*aryuZZ[21];
   aryuZZ[194]=aryuZZ[197] + 59./32.*aryuZZ[50] + 3./4.*aryuZZ[194] + 
   aryuZZ[41];
   aryuZZ[194]=aryuZZ[99]*aryuZZ[194];
   aryuZZ[198]= - 1./2.*aryuZZ[20];
   aryuZZ[199]= - aryuZZ[17] + aryuZZ[21] + aryuZZ[198] + aryuZZ[112]
    + aryuZZ[111] + aryuZZ[68];
   aryuZZ[199]=aryuZZ[101]*aryuZZ[199];
   aryuZZ[200]= - MMZ*aryuZZ[99];
   aryuZZ[201]=11./3. + aryuZZ[114];
   aryuZZ[201]=aryuZZ[35]*aryuZZ[201];
   aryuZZ[202]=23./48. - aryuZZ[36];
   aryuZZ[201]=1./2.*aryuZZ[201] + aryuZZ[202] + 3./16.*aryuZZ[200];
   aryuZZ[201]=aryuZZ[10]*aryuZZ[201];
   aryuZZ[203]= - 3*aryuZZ[49] + 7 + 3*aryuZZ[61];
   aryuZZ[203]=aryuZZ[15] + 3./4.*aryuZZ[44] - 5./2.*aryuZZ[67] + 3./2.
   *aryuZZ[64] + 1./4.*aryuZZ[203] - aryuZZ[59];
   aryuZZ[203]=aryuZZ[99]*aryuZZ[203];
   aryuZZ[204]= - aryuZZ[35]*aryuZZ[99]*aryuZZ[44];
   aryuZZ[203]=aryuZZ[175] + aryuZZ[203] + 3./4.*aryuZZ[204];
   aryuZZ[203]=MMH*aryuZZ[203];
   aryuZZ[205]= - aryuZZ[34] - 7./8.*aryuZZ[99];
   aryuZZ[206]= - aryuZZ[34] + 9./8.*aryuZZ[99];
   aryuZZ[206]=1./2.*aryuZZ[206] + aryuZZ[101];
   aryuZZ[206]=aryuZZ[35]*aryuZZ[206];
   aryuZZ[205]=1./2.*aryuZZ[205] + aryuZZ[206];
   aryuZZ[205]=aryuZZ[18]*aryuZZ[205];
   aryuZZ[206]=aryuZZ[98]*aryuZZ[70];
   aryuZZ[207]= - aryuZZ[50]*aryuZZ[98];
   aryuZZ[208]= - aryuZZ[22]*aryuZZ[98];
   aryuZZ[209]= - aryuZZ[21]*aryuZZ[34];
   aryuZZ[210]=11./2.*aryuZZ[13];
   aryuZZ[211]=113./9. + aryuZZ[210];
   aryuZZ[211]=aryuZZ[36]*aryuZZ[211];
   aryuZZ[212]=81./32.*aryuZZ[99] - 5*aryuZZ[98] - aryuZZ[34];
   aryuZZ[212]=aryuZZ[19]*aryuZZ[212];
   aryuZZ[213]= - aryuZZ[17]*aryuZZ[99];
   aryuZZ[214]= - 3./64.*aryuZZ[44];
   aryuZZ[215]= - 55./12.*aryuZZ[12];
   aryuZZ[216]=1./3.*aryuZZ[11];
   aryuZZ[102]=aryuZZ[125] + aryuZZ[159] + 1./16.*aryuZZ[203] + 
   aryuZZ[113] + aryuZZ[165] + aryuZZ[102] + 1./2.*aryuZZ[205] + 
   aryuZZ[201] + aryuZZ[186] + 1./3.*aryuZZ[189] + 1./2.*aryuZZ[179] + 
   aryuZZ[136] + aryuZZ[199] + 3./16.*aryuZZ[213] + 1./2.*aryuZZ[212]
    + 1./2.*aryuZZ[211] + 1./2.*aryuZZ[194] + aryuZZ[216] - 77./12.*
   aryuZZ[13] + aryuZZ[215] + 1./4.*aryuZZ[209] + 13./16.*aryuZZ[15] + 
   aryuZZ[193] + aryuZZ[214] + 1./2.*aryuZZ[208] + 7./4.*aryuZZ[207] + 
   3./2.*aryuZZ[206] - 303./256.*aryuZZ[48] - 39./64.*aryuZZ[46] + 9./
   32.*aryuZZ[67] + 1./2.*aryuZZ[178] + aryuZZ[64];
   aryuZZ[102]=aryuZZ[133]*aryuZZ[102];
   aryuZZ[113]= - 2./3.*aryuZZ[30];
   aryuZZ[125]=1./3.*aryuZZ[28];
   aryuZZ[136]= - 1./3.*aryuZZ[59];
   aryuZZ[159]= - 2*aryuZZ[64] + aryuZZ[136] - aryuZZ[14] + aryuZZ[125]
    - 9./4. + aryuZZ[113];
   aryuZZ[159]=aryuZZ[101]*aryuZZ[159];
   aryuZZ[165]=aryuZZ[65] - aryuZZ[66];
   aryuZZ[165]=MMZ*aryuZZ[108]*aryuZZ[165];
   aryuZZ[178]=4*aryuZZ[56] + aryuZZ[27];
   aryuZZ[178]=1./3.*aryuZZ[178] - 10*aryuZZ[55];
   aryuZZ[179]= - 17./3.*aryuZZ[66] - 23 - 25./3.*aryuZZ[65];
   aryuZZ[179]=aryuZZ[98]*aryuZZ[179];
   aryuZZ[186]= - aryuZZ[50]*aryuZZ[108];
   aryuZZ[108]=aryuZZ[22]*aryuZZ[108];
   aryuZZ[189]= - 13./6.*aryuZZ[49] - 1./2.*aryuZZ[65] + 5./6.*
   aryuZZ[62] - 17./8. + 4./3.*aryuZZ[63];
   aryuZZ[189]=aryuZZ[34]*aryuZZ[189];
   aryuZZ[193]=209 - 225*aryuZZ[63];
   aryuZZ[193]=1423./24.*aryuZZ[49] + 675./16.*aryuZZ[66] + 5./2.*
   aryuZZ[61] + 531./8.*aryuZZ[65] + 1./8.*aryuZZ[193] - 39*aryuZZ[62];
   aryuZZ[193]=5*aryuZZ[15] + 5./4.*aryuZZ[44] - 225./16.*aryuZZ[48] - 
   39./4.*aryuZZ[46] - 15./2.*aryuZZ[67] + 1./2.*aryuZZ[193] - 5*
   aryuZZ[59];
   aryuZZ[193]=aryuZZ[99]*aryuZZ[193];
   aryuZZ[108]=7./4.*aryuZZ[165] + aryuZZ[159] + aryuZZ[109] + 1./8.*
   aryuZZ[193] + aryuZZ[189] + aryuZZ[108] + aryuZZ[186] + 1./4.*
   aryuZZ[179] + 8*aryuZZ[57] - 2./3.*aryuZZ[25] + 16*aryuZZ[51] + 1./3.
   *aryuZZ[178] + 12*aryuZZ[52];
   aryuZZ[108]=MMZ*aryuZZ[108];
   aryuZZ[109]= - aryuZZ[19]*aryuZZ[98];
   aryuZZ[109]=aryuZZ[185] + 1./4.*aryuZZ[109] - 5./6.*aryuZZ[11] + 289.
   /4.*aryuZZ[13] + 1967./54. + aryuZZ[12];
   aryuZZ[165]=31./3. - aryuZZ[7];
   aryuZZ[178]= - 8./9.*aryuZZ[35];
   aryuZZ[165]=aryuZZ[178] + 5./54.*aryuZZ[165] - aryuZZ[36];
   aryuZZ[165]=aryuZZ[1]*aryuZZ[165];
   aryuZZ[179]= - 5*aryuZZ[7];
   aryuZZ[186]=41 + aryuZZ[179];
   aryuZZ[189]= - 11*aryuZZ[36];
   aryuZZ[186]=1./2.*aryuZZ[186] + aryuZZ[189];
   aryuZZ[186]=1./3.*aryuZZ[186] - 8*aryuZZ[35];
   aryuZZ[186]=aryuZZ[39]*aryuZZ[186];
   aryuZZ[193]=1./3.*aryuZZ[101] + aryuZZ[98] + aryuZZ[156];
   aryuZZ[193]=MMZ*aryuZZ[193];
   aryuZZ[194]=2 + aryuZZ[114];
   aryuZZ[194]=1./3.*aryuZZ[10]*aryuZZ[194];
   aryuZZ[199]=1./3.*aryuZZ[117];
   aryuZZ[201]= - 4./3.*aryuZZ[1];
   aryuZZ[203]= - 4*aryuZZ[39] + 5./4. + aryuZZ[201];
   aryuZZ[203]=aryuZZ[6]*aryuZZ[203];
   aryuZZ[205]= - 7./2.*aryuZZ[35];
   aryuZZ[109]=1./3.*aryuZZ[203] + aryuZZ[199] + aryuZZ[194] + 1./3.*
   aryuZZ[186] + aryuZZ[165] + aryuZZ[205] + 1./3.*aryuZZ[109] + 
   aryuZZ[193];
   aryuZZ[109]=aryuZZ[6]*aryuZZ[109];
   aryuZZ[165]= - 11*aryuZZ[19] + aryuZZ[17];
   aryuZZ[186]=13./4.*aryuZZ[18];
   aryuZZ[165]=aryuZZ[186] + 1./2.*aryuZZ[165] - 7*MMZ;
   aryuZZ[165]=aryuZZ[6]*aryuZZ[165];
   aryuZZ[193]= - 2*aryuZZ[50];
   aryuZZ[203]=aryuZZ[193] - aryuZZ[8] + aryuZZ[41];
   aryuZZ[211]=29./9.*aryuZZ[19];
   aryuZZ[212]= - 2*aryuZZ[36];
   aryuZZ[217]=49./9. + aryuZZ[212];
   aryuZZ[217]=MMZ*aryuZZ[217];
   aryuZZ[218]= - 47*aryuZZ[19] - aryuZZ[17];
   aryuZZ[218]=1./2.*aryuZZ[218] - 34*MMZ;
   aryuZZ[218]=aryuZZ[35]*aryuZZ[218];
   aryuZZ[219]= - 35./12.*aryuZZ[35] + 10./9. + aryuZZ[121];
   aryuZZ[219]=aryuZZ[18]*aryuZZ[219];
   aryuZZ[220]= - 2*aryuZZ[46];
   aryuZZ[221]= - 71./3. + aryuZZ[220];
   aryuZZ[221]=aryuZZ[37]*aryuZZ[221];
   aryuZZ[165]=aryuZZ[221] + 1./3.*aryuZZ[165] + aryuZZ[219] + 1./3.*
   aryuZZ[218] + 2*aryuZZ[217] + aryuZZ[203] + aryuZZ[211];
   aryuZZ[165]=aryuZZ[3]*aryuZZ[165];
   aryuZZ[217]=39*aryuZZ[46];
   aryuZZ[218]=225./4.*aryuZZ[48];
   aryuZZ[221]= - 5*aryuZZ[44];
   aryuZZ[222]=aryuZZ[221] + aryuZZ[218] - 1921./6. + aryuZZ[217];
   aryuZZ[222]=aryuZZ[99]*aryuZZ[222];
   aryuZZ[223]=4./3.*aryuZZ[101];
   aryuZZ[224]= - aryuZZ[34] + 17./4.*aryuZZ[99];
   aryuZZ[224]=aryuZZ[35]*aryuZZ[224];
   aryuZZ[225]= - aryuZZ[10]*aryuZZ[99];
   aryuZZ[226]=5./2.*aryuZZ[225];
   aryuZZ[227]=1./3.*aryuZZ[98] + aryuZZ[171];
   aryuZZ[227]=aryuZZ[6]*aryuZZ[227];
   aryuZZ[222]=aryuZZ[227] + aryuZZ[226] + 5./6.*aryuZZ[224] + 
   aryuZZ[173] + aryuZZ[223] + 1./8.*aryuZZ[222] + 8./3.*aryuZZ[98] + 9
   *aryuZZ[34];
   aryuZZ[222]=aryuZZ[37]*aryuZZ[222];
   aryuZZ[224]=3*aryuZZ[45];
   aryuZZ[227]= - 2 + aryuZZ[224];
   aryuZZ[227]=aryuZZ[3]*aryuZZ[37]*aryuZZ[227];
   aryuZZ[228]=6*aryuZZ[227];
   aryuZZ[229]=1./2.*aryuZZ[36];
   aryuZZ[230]=aryuZZ[228] - 104./3.*aryuZZ[35] + aryuZZ[229] + 
   aryuZZ[158] - 73./6. - aryuZZ[46];
   aryuZZ[230]=aryuZZ[3]*aryuZZ[230];
   aryuZZ[231]= - 4./3.*aryuZZ[58] - aryuZZ[56];
   aryuZZ[232]=7*aryuZZ[55];
   aryuZZ[231]= - 11*aryuZZ[52] + 2*aryuZZ[231] + aryuZZ[232];
   aryuZZ[233]= - 2*aryuZZ[51];
   aryuZZ[231]= - 2./3.*aryuZZ[53] - 2*aryuZZ[57] + 1./3.*aryuZZ[231]
    + aryuZZ[233];
   aryuZZ[161]=6*aryuZZ[161] + 11./3.*aryuZZ[164];
   aryuZZ[161]=aryuZZ[10]*aryuZZ[161];
   aryuZZ[234]= - 6*aryuZZ[38];
   aryuZZ[235]=aryuZZ[234] + 29./3.*aryuZZ[15] - 2./3.*aryuZZ[67] - 22*
   aryuZZ[64] - 11./3.*aryuZZ[59] - 73./6. - 6*aryuZZ[72];
   aryuZZ[235]=aryuZZ[101]*aryuZZ[235];
   aryuZZ[236]=12*aryuZZ[167];
   aryuZZ[237]=11./3.*aryuZZ[168];
   aryuZZ[230]=aryuZZ[236] + aryuZZ[230] + aryuZZ[161] + aryuZZ[237] + 
   2*aryuZZ[231] + aryuZZ[235];
   aryuZZ[230]=MMt*aryuZZ[230];
   aryuZZ[231]=aryuZZ[221] + aryuZZ[218] + 58175./81. + 31*aryuZZ[46];
   aryuZZ[238]=aryuZZ[19]*aryuZZ[34];
   aryuZZ[239]=aryuZZ[101]*aryuZZ[17];
   aryuZZ[240]=1./3.*aryuZZ[239];
   aryuZZ[241]=35./3.*aryuZZ[12];
   aryuZZ[242]= - 17./9.*aryuZZ[11];
   aryuZZ[231]=aryuZZ[240] + 5./3.*aryuZZ[238] + aryuZZ[242] + 487./6.*
   aryuZZ[13] + 1./16.*aryuZZ[231] + aryuZZ[241];
   aryuZZ[217]=aryuZZ[221] + aryuZZ[218] + 71 + aryuZZ[217];
   aryuZZ[217]=aryuZZ[99]*aryuZZ[217];
   aryuZZ[218]= - 1./3.*aryuZZ[101];
   aryuZZ[217]=1./32.*aryuZZ[217] + aryuZZ[218];
   aryuZZ[217]=MMZ*aryuZZ[217];
   aryuZZ[238]=7./4. + 1./3.*aryuZZ[200];
   aryuZZ[238]=aryuZZ[35]*aryuZZ[238];
   aryuZZ[217]=aryuZZ[238] + 1./2.*aryuZZ[231] + aryuZZ[217];
   aryuZZ[217]=aryuZZ[35]*aryuZZ[217];
   aryuZZ[231]=aryuZZ[49] - 7./3. - aryuZZ[61];
   aryuZZ[238]=1./3.*aryuZZ[59];
   aryuZZ[243]= - 1./3.*aryuZZ[15];
   aryuZZ[244]= - 1./2.*aryuZZ[64];
   aryuZZ[231]=aryuZZ[243] + aryuZZ[143] + 5./6.*aryuZZ[67] + 
   aryuZZ[244] + 1./4.*aryuZZ[231] + aryuZZ[238];
   aryuZZ[231]=aryuZZ[99]*aryuZZ[231];
   aryuZZ[245]=aryuZZ[35]*aryuZZ[99]*aryuZZ[44];
   aryuZZ[231]=1./3.*aryuZZ[225] + aryuZZ[231] + 1./4.*aryuZZ[245];
   aryuZZ[231]=5./8.*MMH*aryuZZ[231];
   aryuZZ[245]= - 29./27. - aryuZZ[177];
   aryuZZ[246]=11./9.*aryuZZ[7];
   aryuZZ[247]= - 17*aryuZZ[7];
   aryuZZ[248]=65./3. + aryuZZ[247];
   aryuZZ[248]=aryuZZ[35]*aryuZZ[248];
   aryuZZ[245]=1./6.*aryuZZ[248] + 11./9.*aryuZZ[36] + 2*aryuZZ[245] + 
   aryuZZ[246];
   aryuZZ[245]=aryuZZ[39]*aryuZZ[245];
   aryuZZ[249]=MMZ*aryuZZ[99];
   aryuZZ[250]=MMZ*aryuZZ[101];
   aryuZZ[251]= - 2 + aryuZZ[250];
   aryuZZ[251]=aryuZZ[35]*aryuZZ[251];
   aryuZZ[252]=23./24. - aryuZZ[36];
   aryuZZ[251]=1./3.*aryuZZ[251] + aryuZZ[252] + 5./8.*aryuZZ[249];
   aryuZZ[251]=aryuZZ[10]*aryuZZ[251];
   aryuZZ[253]= - 37./27. - aryuZZ[177];
   aryuZZ[253]=2*aryuZZ[253] + aryuZZ[246];
   aryuZZ[248]=1./18.*aryuZZ[248] + 1./3.*aryuZZ[253] + aryuZZ[36];
   aryuZZ[248]=aryuZZ[1]*aryuZZ[248];
   aryuZZ[253]=5./3.*aryuZZ[34];
   aryuZZ[254]=aryuZZ[253] - 39./8.*aryuZZ[99];
   aryuZZ[254]=1./2.*aryuZZ[254] + aryuZZ[218];
   aryuZZ[254]=aryuZZ[35]*aryuZZ[254];
   aryuZZ[255]=5*aryuZZ[34] + 59./8.*aryuZZ[99];
   aryuZZ[254]=1./6.*aryuZZ[255] + aryuZZ[254];
   aryuZZ[254]=aryuZZ[18]*aryuZZ[254];
   aryuZZ[255]= - 1./3.*aryuZZ[28];
   aryuZZ[256]=5./32.*aryuZZ[61];
   aryuZZ[257]= - 23./72.*aryuZZ[14];
   aryuZZ[258]=28*aryuZZ[16];
   aryuZZ[259]= - 5./12.*aryuZZ[177];
   aryuZZ[260]= - 7./24.*aryuZZ[59];
   aryuZZ[261]= - 2./3.*aryuZZ[64];
   aryuZZ[262]= - 15./16.*aryuZZ[67];
   aryuZZ[263]= - aryuZZ[98]*aryuZZ[70];
   aryuZZ[264]=aryuZZ[50]*aryuZZ[98];
   aryuZZ[265]=aryuZZ[22]*aryuZZ[98];
   aryuZZ[266]=5./32.*aryuZZ[44];
   aryuZZ[267]= - 2*aryuZZ[70] + aryuZZ[42];
   aryuZZ[268]=5./2.*aryuZZ[22] + aryuZZ[50] + 5./2.*aryuZZ[41] + 4*
   aryuZZ[267] - 5*aryuZZ[69];
   aryuZZ[268]=aryuZZ[34]*aryuZZ[268];
   aryuZZ[269]=aryuZZ[21]*aryuZZ[34];
   aryuZZ[270]= - 37./18.*aryuZZ[12];
   aryuZZ[271]=11./27.*aryuZZ[11];
   aryuZZ[272]=75./8.*aryuZZ[192] + 13*aryuZZ[69];
   aryuZZ[272]=5./4.*aryuZZ[40] + 3*aryuZZ[272] - 5./2.*aryuZZ[68];
   aryuZZ[272]= - 73./12.*aryuZZ[21] - 307./32.*aryuZZ[50] + 1./4.*
   aryuZZ[272] - 11./3.*aryuZZ[41];
   aryuZZ[272]=aryuZZ[99]*aryuZZ[272];
   aryuZZ[273]=45./2.*aryuZZ[13];
   aryuZZ[274]=61./9. + aryuZZ[273];
   aryuZZ[274]=aryuZZ[36]*aryuZZ[274];
   aryuZZ[275]=11*aryuZZ[34];
   aryuZZ[276]=37./2.*aryuZZ[98] + aryuZZ[275];
   aryuZZ[276]=1./3.*aryuZZ[276] - 225./16.*aryuZZ[99];
   aryuZZ[276]=aryuZZ[19]*aryuZZ[276];
   aryuZZ[277]=aryuZZ[17]*aryuZZ[99];
   aryuZZ[278]=5./8.*aryuZZ[277];
   aryuZZ[279]=aryuZZ[41] - aryuZZ[8] - 2*aryuZZ[68];
   aryuZZ[279]=1./3.*aryuZZ[279] + aryuZZ[17];
   aryuZZ[279]=aryuZZ[101]*aryuZZ[279];
   aryuZZ[280]=5./8.*aryuZZ[15];
   aryuZZ[102]=aryuZZ[102] + aryuZZ[230] + aryuZZ[231] + aryuZZ[165] + 
   aryuZZ[222] + aryuZZ[109] + aryuZZ[254] + aryuZZ[251] + 1./3.*
   aryuZZ[245] + 1./3.*aryuZZ[248] + aryuZZ[217] + aryuZZ[108] + 
   aryuZZ[279] + aryuZZ[278] + 1./2.*aryuZZ[276] + 1./2.*aryuZZ[274] + 
   aryuZZ[272] + aryuZZ[271] + 415./36.*aryuZZ[13] + aryuZZ[270] + 5./6.
   *aryuZZ[269] + aryuZZ[280] + 1./3.*aryuZZ[268] + aryuZZ[266] + 11./
   12.*aryuZZ[265] + 4./3.*aryuZZ[264] + 2*aryuZZ[263] + 287./128.*
   aryuZZ[48] + 41./32.*aryuZZ[46] + aryuZZ[262] + aryuZZ[261] + 
   aryuZZ[260] + aryuZZ[259] + aryuZZ[258] + 2383./384.*aryuZZ[49] - 
   535./768.*aryuZZ[66] + aryuZZ[257] + aryuZZ[256] - 5287./384.*
   aryuZZ[65] - 157./48.*aryuZZ[62] + aryuZZ[255] - 1315./384.*
   aryuZZ[63] + 5./3.*aryuZZ[31] - 620713./31104. - 28*aryuZZ[60];
   aryuZZ[102]=aryuZZ[133]*aryuZZ[102];
   aryuZZ[108]= - 7./3. - 85*aryuZZ[12];
   aryuZZ[108]=1./2.*aryuZZ[108] + 5*aryuZZ[13];
   aryuZZ[108]=1./3.*aryuZZ[108] + 5*aryuZZ[116];
   aryuZZ[108]=1./2.*aryuZZ[108] + 5*aryuZZ[250];
   aryuZZ[109]= - 1 + 17./9.*aryuZZ[35];
   aryuZZ[165]=aryuZZ[1]*aryuZZ[109];
   aryuZZ[109]=aryuZZ[39]*aryuZZ[109];
   aryuZZ[217]=1 + 1./3.*aryuZZ[114];
   aryuZZ[222]=aryuZZ[10]*aryuZZ[217];
   aryuZZ[230]=11./9.*aryuZZ[39];
   aryuZZ[245]=aryuZZ[230] - 17./72. + aryuZZ[1];
   aryuZZ[245]=aryuZZ[6]*aryuZZ[245];
   aryuZZ[108]=5./9.*aryuZZ[245] + 5./18.*aryuZZ[117] + 5./6.*
   aryuZZ[222] + 11./9.*aryuZZ[109] + 1./18.*aryuZZ[108] + aryuZZ[165];
   aryuZZ[108]=aryuZZ[6]*aryuZZ[108];
   aryuZZ[109]= - 1./2. + 17*aryuZZ[46];
   aryuZZ[109]=aryuZZ[122] + aryuZZ[36] + 1./3.*aryuZZ[109] + 
   aryuZZ[158];
   aryuZZ[109]=1./2.*aryuZZ[109] + 3*aryuZZ[227];
   aryuZZ[109]=aryuZZ[3]*aryuZZ[109];
   aryuZZ[146]=aryuZZ[158] + 47./18.*aryuZZ[15] + 17./9.*aryuZZ[67] + 7.
   /3.*aryuZZ[64] + 7./18.*aryuZZ[59] + 161./36. + aryuZZ[146];
   aryuZZ[146]=aryuZZ[101]*aryuZZ[146];
   aryuZZ[163]=aryuZZ[163] + 7./18.*aryuZZ[168];
   aryuZZ[163]=aryuZZ[10]*aryuZZ[163];
   aryuZZ[165]= - 22*aryuZZ[55] + 17*aryuZZ[53];
   aryuZZ[164]=7./18.*aryuZZ[164];
   aryuZZ[109]=6*aryuZZ[167] + aryuZZ[109] + aryuZZ[163] + aryuZZ[164]
    + 2./9.*aryuZZ[165] + aryuZZ[146];
   aryuZZ[109]=MMt*aryuZZ[109];
   aryuZZ[146]= - 5*aryuZZ[8] - 17*aryuZZ[41];
   aryuZZ[163]=11*aryuZZ[21];
   aryuZZ[146]= - 44./3.*MMZ - 11./6.*aryuZZ[17] + aryuZZ[163] + 11./6.
   *aryuZZ[20] + 1./2.*aryuZZ[146] + 17*aryuZZ[50];
   aryuZZ[165]=1./9.*aryuZZ[17] + MMZ;
   aryuZZ[167]=aryuZZ[35]*aryuZZ[165];
   aryuZZ[168]=11./9. + aryuZZ[162];
   aryuZZ[168]=aryuZZ[18]*aryuZZ[168];
   aryuZZ[222]=aryuZZ[165] + aryuZZ[18];
   aryuZZ[222]=aryuZZ[6]*aryuZZ[222];
   aryuZZ[227]=15./2. + 17./3.*aryuZZ[46];
   aryuZZ[227]=aryuZZ[37]*aryuZZ[227];
   aryuZZ[146]=aryuZZ[227] + 5./4.*aryuZZ[222] + aryuZZ[168] + 1./3.*
   aryuZZ[146] + 17./4.*aryuZZ[167];
   aryuZZ[146]=aryuZZ[3]*aryuZZ[146];
   aryuZZ[167]=103./3. - 1./2.*aryuZZ[63];
   aryuZZ[142]=aryuZZ[142] + 1./4.*aryuZZ[167] + 25./3.*aryuZZ[62];
   aryuZZ[142]= - 185./72.*aryuZZ[49] + 1./3.*aryuZZ[142] + 1./16.*
   aryuZZ[66];
   aryuZZ[167]=1./2.*aryuZZ[67];
   aryuZZ[142]=aryuZZ[243] - 1./12.*aryuZZ[44] - 1./48.*aryuZZ[48] + 25.
   /36.*aryuZZ[46] + aryuZZ[167] + 1./2.*aryuZZ[142] + aryuZZ[238];
   aryuZZ[142]=aryuZZ[99]*aryuZZ[142];
   aryuZZ[168]=257*aryuZZ[55] + 17*aryuZZ[25];
   aryuZZ[222]=1./3.*aryuZZ[49] + 1./4. - 1./3.*aryuZZ[62];
   aryuZZ[222]=aryuZZ[34]*aryuZZ[222];
   aryuZZ[168]=1./9.*aryuZZ[168] + 107./4.*aryuZZ[222];
   aryuZZ[222]= - 1./2.*aryuZZ[14] + 1./6.*aryuZZ[28] + 13./8. - 1./3.*
   aryuZZ[30];
   aryuZZ[222]= - 11./3.*aryuZZ[15] + 17*aryuZZ[64] + 5*aryuZZ[222] + 
   17./6.*aryuZZ[59];
   aryuZZ[222]=aryuZZ[101]*aryuZZ[222];
   aryuZZ[142]=1./3.*aryuZZ[222] + 1./9.*aryuZZ[168] + 25./16.*
   aryuZZ[142];
   aryuZZ[142]=MMZ*aryuZZ[142];
   aryuZZ[168]= - aryuZZ[49] + 7./3. + aryuZZ[61];
   aryuZZ[222]=1./4.*aryuZZ[44];
   aryuZZ[227]=1./3.*aryuZZ[15];
   aryuZZ[136]=aryuZZ[227] + aryuZZ[222] - 5./6.*aryuZZ[67] + 1./2.*
   aryuZZ[64] + 1./4.*aryuZZ[168] + aryuZZ[136];
   aryuZZ[136]=aryuZZ[99]*aryuZZ[136];
   aryuZZ[168]=1./3.*aryuZZ[175];
   aryuZZ[136]=aryuZZ[168] + aryuZZ[136] + 1./4.*aryuZZ[204];
   aryuZZ[136]=MMH*aryuZZ[136];
   aryuZZ[195]=aryuZZ[196] + aryuZZ[195] + 1./8.*aryuZZ[192] - 25./3.*
   aryuZZ[69];
   aryuZZ[195]=47./36.*aryuZZ[21] + 941./288.*aryuZZ[50] + 1./4.*
   aryuZZ[195] + 7./9.*aryuZZ[41];
   aryuZZ[195]=aryuZZ[99]*aryuZZ[195];
   aryuZZ[196]= - 4841./9. - 217*aryuZZ[46];
   aryuZZ[204]=25*aryuZZ[44];
   aryuZZ[238]=25./4.*aryuZZ[48];
   aryuZZ[196]=aryuZZ[204] + 1./3.*aryuZZ[196] + aryuZZ[238];
   aryuZZ[245]= - 289./9.*aryuZZ[12];
   aryuZZ[196]=1./8.*aryuZZ[196] + aryuZZ[245];
   aryuZZ[196]=17./3.*aryuZZ[116] + 1./2.*aryuZZ[196] + 17./9.*
   aryuZZ[13];
   aryuZZ[248]=1./4.*aryuZZ[48];
   aryuZZ[181]=aryuZZ[44] + 25./3.*aryuZZ[181] + aryuZZ[248];
   aryuZZ[181]=aryuZZ[99]*aryuZZ[181];
   aryuZZ[181]=25./32.*aryuZZ[181] + 17./3.*aryuZZ[101];
   aryuZZ[181]=MMZ*aryuZZ[181];
   aryuZZ[181]= - 1285./108.*aryuZZ[35] + 1./2.*aryuZZ[196] + 
   aryuZZ[181];
   aryuZZ[181]=aryuZZ[35]*aryuZZ[181];
   aryuZZ[196]=841./6. - 25*aryuZZ[46];
   aryuZZ[196]=aryuZZ[44] + 1./3.*aryuZZ[196] + aryuZZ[248];
   aryuZZ[196]=aryuZZ[99]*aryuZZ[196];
   aryuZZ[248]= - 107./3.*aryuZZ[34];
   aryuZZ[196]=aryuZZ[248] + 25./8.*aryuZZ[196];
   aryuZZ[254]=107./3.*aryuZZ[34] - 625./4.*aryuZZ[99];
   aryuZZ[254]=aryuZZ[35]*aryuZZ[254];
   aryuZZ[196]=25./4.*aryuZZ[175] + 1./12.*aryuZZ[254] + 1./2.*
   aryuZZ[196] - 34./3.*aryuZZ[101];
   aryuZZ[196]=aryuZZ[37]*aryuZZ[196];
   aryuZZ[254]=87691./192. - 85*aryuZZ[31];
   aryuZZ[254]= - 115./72.*aryuZZ[14] - 25./32.*aryuZZ[61] + 10721./432.
   *aryuZZ[62] - 5./3.*aryuZZ[28] + 1./27.*aryuZZ[254] - 25./128.*
   aryuZZ[63];
   aryuZZ[254]= - 61./72.*aryuZZ[59] + 85./324.*aryuZZ[177] - 82393./
   10368.*aryuZZ[49] + 1./3.*aryuZZ[254] + 25./256.*aryuZZ[66];
   aryuZZ[263]= - 49./9. + 25*aryuZZ[200];
   aryuZZ[264]=aryuZZ[35]*aryuZZ[217];
   aryuZZ[263]=1./8.*aryuZZ[263] + 17*aryuZZ[264];
   aryuZZ[263]=aryuZZ[10]*aryuZZ[263];
   aryuZZ[264]= - 107*aryuZZ[34] - 925./8.*aryuZZ[99];
   aryuZZ[248]=aryuZZ[248] + 625./8.*aryuZZ[99];
   aryuZZ[248]=1./2.*aryuZZ[248] + 17*aryuZZ[101];
   aryuZZ[248]=aryuZZ[35]*aryuZZ[248];
   aryuZZ[248]=1./6.*aryuZZ[264] + aryuZZ[248];
   aryuZZ[248]=aryuZZ[18]*aryuZZ[248];
   aryuZZ[264]= - aryuZZ[6]*aryuZZ[12];
   aryuZZ[268]= - aryuZZ[35]*aryuZZ[12];
   aryuZZ[264]=5./2.*aryuZZ[264] + 11./3.*aryuZZ[12] + 17./2.*
   aryuZZ[268];
   aryuZZ[272]=pow(aryuZZ[2],2);
   aryuZZ[264]=aryuZZ[272]*aryuZZ[264];
   aryuZZ[274]= - aryuZZ[50] + aryuZZ[69] + aryuZZ[112];
   aryuZZ[274]=aryuZZ[34]*aryuZZ[274];
   aryuZZ[276]= - aryuZZ[19]*aryuZZ[99];
   aryuZZ[163]= - 20*aryuZZ[17] + aryuZZ[163] - 11./2.*aryuZZ[20] - 17./
   2.*aryuZZ[41] - 5./2.*aryuZZ[8] + 17*aryuZZ[68];
   aryuZZ[163]=aryuZZ[101]*aryuZZ[163];
   aryuZZ[281]= - 17./3.*aryuZZ[35] + 22./9. + 5*aryuZZ[177];
   aryuZZ[282]=aryuZZ[1]*aryuZZ[281];
   aryuZZ[281]=aryuZZ[39]*aryuZZ[281];
   aryuZZ[108]=1./108.*aryuZZ[264] + aryuZZ[109] + 25./48.*aryuZZ[136]
    + aryuZZ[146] + 1./3.*aryuZZ[196] + aryuZZ[108] + 1./18.*
   aryuZZ[248] + 1./6.*aryuZZ[263] + 11./81.*aryuZZ[281] + 1./9.*
   aryuZZ[282] + 1./6.*aryuZZ[181] + aryuZZ[142] + 1./9.*aryuZZ[163] + 
   25./48.*aryuZZ[213] + 25./192.*aryuZZ[276] + 25./6.*aryuZZ[195] - 11.
   /162.*aryuZZ[13] + 187./324.*aryuZZ[12] + 107./108.*aryuZZ[209] + 
   101./144.*aryuZZ[15] + 107./54.*aryuZZ[274] - 25./192.*aryuZZ[44] - 
   25./768.*aryuZZ[48] - 2749./1728.*aryuZZ[46] + 25./32.*aryuZZ[67] + 
   1./2.*aryuZZ[254] + 17./9.*aryuZZ[64];
   aryuZZ[108]=aryuZZ[272]*aryuZZ[108];
   aryuZZ[109]=5*aryuZZ[19];
   aryuZZ[136]=aryuZZ[109] + aryuZZ[17];
   aryuZZ[142]=aryuZZ[186] + 1./2.*aryuZZ[136] + 11./3.*MMZ;
   aryuZZ[142]=aryuZZ[6]*aryuZZ[142];
   aryuZZ[146]= - 11./9.*aryuZZ[19];
   aryuZZ[163]=20./27. - aryuZZ[36];
   aryuZZ[163]=MMZ*aryuZZ[163];
   aryuZZ[181]=17*aryuZZ[19];
   aryuZZ[195]=aryuZZ[181] - aryuZZ[17];
   aryuZZ[195]=1./2.*aryuZZ[195] - 22./3.*MMZ;
   aryuZZ[195]=aryuZZ[35]*aryuZZ[195];
   aryuZZ[196]= - 53./9. + aryuZZ[220];
   aryuZZ[196]=aryuZZ[37]*aryuZZ[196];
   aryuZZ[142]=aryuZZ[196] + 1./3.*aryuZZ[142] + aryuZZ[219] + 1./3.*
   aryuZZ[195] + aryuZZ[163] + aryuZZ[203] + aryuZZ[146];
   aryuZZ[142]=aryuZZ[3]*aryuZZ[142];
   aryuZZ[163]=11*aryuZZ[12];
   aryuZZ[107]=aryuZZ[107] - 5./2.*aryuZZ[11] + 89./12.*aryuZZ[13] + 
   283./18. + aryuZZ[163];
   aryuZZ[107]=aryuZZ[134] + aryuZZ[250] + 1./3.*aryuZZ[107] + 
   aryuZZ[185];
   aryuZZ[195]=433./9. + aryuZZ[179];
   aryuZZ[195]= - 440./9.*aryuZZ[35] + 1./2.*aryuZZ[195] + aryuZZ[189];
   aryuZZ[195]=aryuZZ[39]*aryuZZ[195];
   aryuZZ[196]=347./3. + aryuZZ[179];
   aryuZZ[196]= - 40./9.*aryuZZ[35] + 1./54.*aryuZZ[196] - aryuZZ[36];
   aryuZZ[196]=aryuZZ[1]*aryuZZ[196];
   aryuZZ[203]= - 44./9.*aryuZZ[39] + 55./36. - 4*aryuZZ[1];
   aryuZZ[203]=aryuZZ[6]*aryuZZ[203];
   aryuZZ[107]=1./9.*aryuZZ[203] + aryuZZ[199] + aryuZZ[194] + 1./9.*
   aryuZZ[195] + 1./3.*aryuZZ[107] + aryuZZ[196];
   aryuZZ[107]=aryuZZ[6]*aryuZZ[107];
   aryuZZ[194]= - 3431./2. + 275*aryuZZ[46];
   aryuZZ[195]=7./4.*aryuZZ[48];
   aryuZZ[194]=1./3.*aryuZZ[194] + aryuZZ[195];
   aryuZZ[194]=1./3.*aryuZZ[194] - aryuZZ[44];
   aryuZZ[194]=aryuZZ[99]*aryuZZ[194];
   aryuZZ[196]=173*aryuZZ[34];
   aryuZZ[199]=aryuZZ[98] + aryuZZ[196];
   aryuZZ[203]= - 173*aryuZZ[34];
   aryuZZ[219]=aryuZZ[203] + 3325./4.*aryuZZ[99];
   aryuZZ[219]=aryuZZ[35]*aryuZZ[219];
   aryuZZ[248]= - aryuZZ[6]*aryuZZ[98];
   aryuZZ[194]=1./9.*aryuZZ[248] + aryuZZ[226] + 1./54.*aryuZZ[219] + 
   aryuZZ[223] + 1./9.*aryuZZ[199] + 5./8.*aryuZZ[194];
   aryuZZ[194]=aryuZZ[37]*aryuZZ[194];
   aryuZZ[199]= - 68./27.*aryuZZ[58] - aryuZZ[56];
   aryuZZ[219]=1./3.*aryuZZ[57];
   aryuZZ[223]= - 2*aryuZZ[53];
   aryuZZ[199]=aryuZZ[223] + aryuZZ[219] - aryuZZ[52] + 2*aryuZZ[199]
    + aryuZZ[232];
   aryuZZ[226]=aryuZZ[228] - 152./9.*aryuZZ[35] + aryuZZ[229] + 
   aryuZZ[158] + 101./18. - aryuZZ[46];
   aryuZZ[226]=aryuZZ[3]*aryuZZ[226];
   aryuZZ[161]=aryuZZ[236] + aryuZZ[226] + aryuZZ[161] + aryuZZ[237] + 
   2./3.*aryuZZ[199] + aryuZZ[235];
   aryuZZ[161]=MMt*aryuZZ[161];
   aryuZZ[199]=41407./9. + 1303*aryuZZ[46];
   aryuZZ[226]=35./4.*aryuZZ[48];
   aryuZZ[199]=1./3.*aryuZZ[199] + aryuZZ[226];
   aryuZZ[199]=1./3.*aryuZZ[199] + aryuZZ[221];
   aryuZZ[228]=73./9.*aryuZZ[12];
   aryuZZ[184]=aryuZZ[240] + 11./3.*aryuZZ[184] + aryuZZ[242] + 263./54.
   *aryuZZ[13] + 1./16.*aryuZZ[199] + aryuZZ[228];
   aryuZZ[199]=29 + 55./3.*aryuZZ[46];
   aryuZZ[195]=5*aryuZZ[199] + aryuZZ[195];
   aryuZZ[195]=1./3.*aryuZZ[195] - aryuZZ[44];
   aryuZZ[195]=aryuZZ[99]*aryuZZ[195];
   aryuZZ[195]=5./32.*aryuZZ[195] + aryuZZ[218];
   aryuZZ[195]=MMZ*aryuZZ[195];
   aryuZZ[199]=203./12. + 5*aryuZZ[200];
   aryuZZ[199]=aryuZZ[35]*aryuZZ[199];
   aryuZZ[184]=5./27.*aryuZZ[199] + 1./2.*aryuZZ[184] + aryuZZ[195];
   aryuZZ[184]=aryuZZ[35]*aryuZZ[184];
   aryuZZ[195]=5*aryuZZ[27];
   aryuZZ[199]= - 22*aryuZZ[25] - 526*aryuZZ[55] + 68*aryuZZ[56] + 
   aryuZZ[195];
   aryuZZ[200]= - 1./3.*aryuZZ[49] - 1./4. + 1./3.*aryuZZ[62];
   aryuZZ[200]=aryuZZ[34]*aryuZZ[200];
   aryuZZ[199]=173./2.*aryuZZ[200] + 1./9.*aryuZZ[199] + aryuZZ[57];
   aryuZZ[200]= - 2203./3. - 7*aryuZZ[63];
   aryuZZ[200]=5./8.*aryuZZ[65] + 1./8.*aryuZZ[200] - 275./3.*
   aryuZZ[62];
   aryuZZ[218]=1./2.*aryuZZ[61];
   aryuZZ[200]=5915./216.*aryuZZ[49] + 7./16.*aryuZZ[66] + 1./3.*
   aryuZZ[200] + aryuZZ[218];
   aryuZZ[200]=aryuZZ[15] + aryuZZ[222] - 7./48.*aryuZZ[48] - 275./36.*
   aryuZZ[46] - 3./2.*aryuZZ[67] + 1./2.*aryuZZ[200] - aryuZZ[59];
   aryuZZ[200]=aryuZZ[99]*aryuZZ[200];
   aryuZZ[159]=aryuZZ[159] + 1./9.*aryuZZ[199] + 5./8.*aryuZZ[200];
   aryuZZ[159]=MMZ*aryuZZ[159];
   aryuZZ[199]= - 85./27. - aryuZZ[177];
   aryuZZ[199]=2*aryuZZ[199] + aryuZZ[246];
   aryuZZ[200]=257./3. + aryuZZ[247];
   aryuZZ[200]=aryuZZ[35]*aryuZZ[200];
   aryuZZ[199]=1./18.*aryuZZ[200] + 1./3.*aryuZZ[199] + aryuZZ[36];
   aryuZZ[199]=aryuZZ[1]*aryuZZ[199];
   aryuZZ[200]= - 31./9. - aryuZZ[177];
   aryuZZ[200]=aryuZZ[36] + 2./3.*aryuZZ[200] + aryuZZ[7];
   aryuZZ[222]=1033./27. + aryuZZ[247];
   aryuZZ[222]=aryuZZ[35]*aryuZZ[222];
   aryuZZ[200]=11./3.*aryuZZ[200] + 1./2.*aryuZZ[222];
   aryuZZ[200]=aryuZZ[39]*aryuZZ[200];
   aryuZZ[196]=aryuZZ[196] + 1555./8.*aryuZZ[99];
   aryuZZ[222]=173./3.*aryuZZ[34] - 1375./8.*aryuZZ[99];
   aryuZZ[222]=1./6.*aryuZZ[222] - aryuZZ[101];
   aryuZZ[222]=aryuZZ[35]*aryuZZ[222];
   aryuZZ[196]=1./18.*aryuZZ[196] + aryuZZ[222];
   aryuZZ[196]=aryuZZ[18]*aryuZZ[196];
   aryuZZ[222]= - 35939./128. + 55*aryuZZ[31];
   aryuZZ[222]= - 565./384.*aryuZZ[65] - 9359./432.*aryuZZ[62] - 
   aryuZZ[28] + 1./27.*aryuZZ[222] - 163./128.*aryuZZ[63];
   aryuZZ[232]=1./2.*aryuZZ[41];
   aryuZZ[235]=aryuZZ[50] - aryuZZ[69] + aryuZZ[232];
   aryuZZ[235]=173./9.*aryuZZ[235] - 11./2.*aryuZZ[22];
   aryuZZ[235]=aryuZZ[34]*aryuZZ[235];
   aryuZZ[236]=7./8.*aryuZZ[192] + 275./3.*aryuZZ[69];
   aryuZZ[236]=1./4.*aryuZZ[40] + 1./3.*aryuZZ[236] - 1./2.*aryuZZ[68];
   aryuZZ[236]= - 541./108.*aryuZZ[21] - 9773./864.*aryuZZ[50] + 1./4.*
   aryuZZ[236] - 71./27.*aryuZZ[41];
   aryuZZ[236]=aryuZZ[99]*aryuZZ[236];
   aryuZZ[240]= - 1./2.*aryuZZ[13];
   aryuZZ[246]=1./3. + aryuZZ[240];
   aryuZZ[248]=aryuZZ[36]*aryuZZ[246];
   aryuZZ[254]= - 35./16.*aryuZZ[99] - 1./6.*aryuZZ[98] + aryuZZ[275];
   aryuZZ[254]=aryuZZ[19]*aryuZZ[254];
   aryuZZ[107]=aryuZZ[108] + aryuZZ[161] + aryuZZ[231] + aryuZZ[142] + 
   aryuZZ[194] + aryuZZ[107] + 1./3.*aryuZZ[196] + aryuZZ[251] + 1./9.*
   aryuZZ[200] + 1./3.*aryuZZ[199] + aryuZZ[184] + aryuZZ[159] + 
   aryuZZ[279] + aryuZZ[278] + 1./6.*aryuZZ[254] + 1./6.*aryuZZ[248] + 
   5*aryuZZ[236] + aryuZZ[271] - 289./324.*aryuZZ[13] - 95./54.*
   aryuZZ[12] + 173./54.*aryuZZ[269] + aryuZZ[280] + 1./3.*aryuZZ[235]
    + aryuZZ[266] + 1./36.*aryuZZ[265] + 1./9.*aryuZZ[207] - 35./384.*
   aryuZZ[48] + 2131./864.*aryuZZ[46] + aryuZZ[262] + aryuZZ[261] + 
   aryuZZ[260] - 55./324.*aryuZZ[177] + 68053./10368.*aryuZZ[49] + 59./
   2304.*aryuZZ[66] + aryuZZ[257] + 1./3.*aryuZZ[222] + aryuZZ[256];
   aryuZZ[107]=aryuZZ[272]*aryuZZ[107];
   aryuZZ[108]=pow(CW,2);
   aryuZZ[142]=aryuZZ[108]*aryuZZ[144];
   aryuZZ[142]=aryuZZ[144] + 1./3.*aryuZZ[142];
   aryuZZ[142]=MMZ*aryuZZ[142];
   aryuZZ[144]=64*aryuZZ[55] - 16*aryuZZ[56] - aryuZZ[27];
   aryuZZ[159]= - 3*aryuZZ[52];
   aryuZZ[144]= - 17./9.*aryuZZ[57] + 16./81.*aryuZZ[25] - 5*aryuZZ[51]
    + 4./81.*aryuZZ[144] + aryuZZ[159];
   aryuZZ[161]= - 2./3.*aryuZZ[57] - aryuZZ[52] + aryuZZ[233];
   aryuZZ[184]=aryuZZ[66] + 4 + aryuZZ[65];
   aryuZZ[184]=aryuZZ[98]*aryuZZ[184];
   aryuZZ[161]=4*aryuZZ[161] + aryuZZ[184];
   aryuZZ[161]=aryuZZ[108]*aryuZZ[161];
   aryuZZ[184]=1 + aryuZZ[65];
   aryuZZ[184]=aryuZZ[108]*aryuZZ[184];
   aryuZZ[184]=4*aryuZZ[184] + 100./9.*aryuZZ[49] - 2*aryuZZ[65] - 64./
   9.*aryuZZ[62] + 19./3. - 4*aryuZZ[63];
   aryuZZ[184]=aryuZZ[34]*aryuZZ[184];
   aryuZZ[194]= - 1 - aryuZZ[65];
   aryuZZ[194]=aryuZZ[108]*aryuZZ[194];
   aryuZZ[196]= - 5*aryuZZ[65] + 64./3.*aryuZZ[62] + 49./3. + 4*
   aryuZZ[63];
   aryuZZ[194]=8./3.*aryuZZ[194] + 4./3.*aryuZZ[48] + 32./9.*aryuZZ[46]
    - 196./27.*aryuZZ[49] + 1./3.*aryuZZ[196] - 2*aryuZZ[66];
   aryuZZ[194]=aryuZZ[99]*aryuZZ[194];
   aryuZZ[196]=8./3.*aryuZZ[66] + 43./4. + 10./3.*aryuZZ[65];
   aryuZZ[196]=aryuZZ[98]*aryuZZ[196];
   aryuZZ[142]=aryuZZ[142] + 1./3.*aryuZZ[156] + 2*aryuZZ[194] + 2./3.*
   aryuZZ[184] + 1./3.*aryuZZ[154] + 1./3.*aryuZZ[161] + 1./3.*
   aryuZZ[153] + 2*aryuZZ[144] + 1./3.*aryuZZ[196];
   aryuZZ[142]=MMZ*aryuZZ[142];
   aryuZZ[144]= - 11./9. - aryuZZ[108];
   aryuZZ[153]= - 8*aryuZZ[108];
   aryuZZ[154]= - 79./3. + aryuZZ[153];
   aryuZZ[154]=aryuZZ[13]*aryuZZ[154];
   aryuZZ[144]=8*aryuZZ[144] + aryuZZ[154];
   aryuZZ[110]=MMZ*aryuZZ[110];
   aryuZZ[154]=2*aryuZZ[35];
   aryuZZ[156]= - 1 + aryuZZ[154];
   aryuZZ[161]=aryuZZ[1]*aryuZZ[156];
   aryuZZ[156]=aryuZZ[39]*aryuZZ[156];
   aryuZZ[184]=40./3.*aryuZZ[39] - 5./3. + 8*aryuZZ[1];
   aryuZZ[184]=aryuZZ[6]*aryuZZ[184];
   aryuZZ[194]=8./3.*aryuZZ[35];
   aryuZZ[110]=4./9.*aryuZZ[184] + 320./27.*aryuZZ[156] + 64./9.*
   aryuZZ[161] + aryuZZ[194] + 4./3.*aryuZZ[144] + aryuZZ[110];
   aryuZZ[110]=aryuZZ[6]*aryuZZ[110];
   aryuZZ[144]= - 59./9. + aryuZZ[220];
   aryuZZ[156]= - 113./3. - 16*aryuZZ[108];
   aryuZZ[156]=aryuZZ[13]*aryuZZ[156];
   aryuZZ[161]= - 2 - aryuZZ[46];
   aryuZZ[161]=8./3.*aryuZZ[161] - aryuZZ[48];
   aryuZZ[161]=MMZ*aryuZZ[99]*aryuZZ[161];
   aryuZZ[184]= - 2./3. + aryuZZ[249];
   aryuZZ[184]=aryuZZ[35]*aryuZZ[184];
   aryuZZ[144]=8./9.*aryuZZ[184] + aryuZZ[161] + 1./3.*aryuZZ[156] - 16.
   /3.*aryuZZ[108] + 4./3.*aryuZZ[144] - aryuZZ[48];
   aryuZZ[144]=aryuZZ[35]*aryuZZ[144];
   aryuZZ[148]=aryuZZ[148] + aryuZZ[134];
   aryuZZ[148]=aryuZZ[35]*aryuZZ[148];
   aryuZZ[156]= - 1./2.*aryuZZ[6];
   aryuZZ[149]=aryuZZ[149] + aryuZZ[156];
   aryuZZ[149]=aryuZZ[6]*aryuZZ[149];
   aryuZZ[161]= - 1./3.*aryuZZ[36];
   aryuZZ[148]=1./2.*aryuZZ[149] + aryuZZ[161] + 1./2.*aryuZZ[148];
   aryuZZ[148]=aryuZZ[133]*aryuZZ[148];
   aryuZZ[149]= - 11./9. - aryuZZ[36];
   aryuZZ[184]=2./3.*aryuZZ[35];
   aryuZZ[196]=1./2.*aryuZZ[149] + aryuZZ[184];
   aryuZZ[196]=aryuZZ[35]*aryuZZ[196];
   aryuZZ[199]= - 7./9. - aryuZZ[36];
   aryuZZ[199]=1./3.*aryuZZ[6] + 1./2.*aryuZZ[199] + aryuZZ[35];
   aryuZZ[199]=aryuZZ[6]*aryuZZ[199];
   aryuZZ[148]=aryuZZ[148] + aryuZZ[199] + aryuZZ[196] + aryuZZ[188] + 
   1./9. + 1./4.*aryuZZ[177];
   aryuZZ[148]=aryuZZ[133]*aryuZZ[148];
   aryuZZ[196]= - 5*aryuZZ[36];
   aryuZZ[199]= - 17*aryuZZ[35];
   aryuZZ[200]= - 7./2.*aryuZZ[6] + aryuZZ[199] + 29./3. + aryuZZ[196];
   aryuZZ[200]=aryuZZ[6]*aryuZZ[200];
   aryuZZ[220]= - 17*aryuZZ[36];
   aryuZZ[222]= - 31./2.*aryuZZ[35] + 65./3. + aryuZZ[220];
   aryuZZ[222]=aryuZZ[35]*aryuZZ[222];
   aryuZZ[200]=1./2.*aryuZZ[200] + 1./2.*aryuZZ[222] + 11./3.*
   aryuZZ[36] - 29./9. - 1./2.*aryuZZ[177];
   aryuZZ[148]=1./9.*aryuZZ[200] + aryuZZ[148];
   aryuZZ[148]=aryuZZ[133]*aryuZZ[148];
   aryuZZ[162]= - 11./3. + aryuZZ[162];
   aryuZZ[162]=aryuZZ[35]*aryuZZ[162];
   aryuZZ[200]= - 11./3. + 17./2.*aryuZZ[35];
   aryuZZ[222]=aryuZZ[200] + 5./4.*aryuZZ[6];
   aryuZZ[222]=aryuZZ[6]*aryuZZ[222];
   aryuZZ[162]=5*aryuZZ[222] + 17*aryuZZ[162] + 121./9. + 25./4.*
   aryuZZ[177];
   aryuZZ[162]=aryuZZ[272]*aryuZZ[162];
   aryuZZ[222]= - 341./9. + 7./2.*aryuZZ[177];
   aryuZZ[222]=1./3.*aryuZZ[222] + 11*aryuZZ[36];
   aryuZZ[231]= - 791./18.*aryuZZ[35] + 1033./27. + aryuZZ[220];
   aryuZZ[231]=aryuZZ[35]*aryuZZ[231];
   aryuZZ[235]=1./18.*aryuZZ[6] - 89./9.*aryuZZ[35] + 133./27. + 
   aryuZZ[196];
   aryuZZ[235]=aryuZZ[6]*aryuZZ[235];
   aryuZZ[162]=1./9.*aryuZZ[162] + 1./2.*aryuZZ[235] + 1./3.*
   aryuZZ[222] + 1./2.*aryuZZ[231];
   aryuZZ[162]=aryuZZ[272]*aryuZZ[162];
   aryuZZ[222]= - 5./3. + aryuZZ[154];
   aryuZZ[222]=aryuZZ[35]*aryuZZ[222];
   aryuZZ[231]=4*aryuZZ[35];
   aryuZZ[235]= - 5./3. + aryuZZ[231];
   aryuZZ[236]=2*aryuZZ[235] + aryuZZ[6];
   aryuZZ[236]=aryuZZ[6]*aryuZZ[236];
   aryuZZ[222]=aryuZZ[236] + 25./9. + 8*aryuZZ[222];
   aryuZZ[162]=16./9.*aryuZZ[222] + aryuZZ[162];
   aryuZZ[148]=1./9.*aryuZZ[162] + aryuZZ[148];
   aryuZZ[148]=aryuZZ[33]*aryuZZ[148];
   aryuZZ[162]= - 4*aryuZZ[46];
   aryuZZ[222]=85./3. + aryuZZ[162];
   aryuZZ[222]=2./3.*aryuZZ[222] - aryuZZ[48];
   aryuZZ[222]=aryuZZ[99]*aryuZZ[222];
   aryuZZ[236]=aryuZZ[34] - 4*aryuZZ[99];
   aryuZZ[236]=aryuZZ[35]*aryuZZ[236];
   aryuZZ[173]=aryuZZ[6]*aryuZZ[173];
   aryuZZ[171]=1./3.*aryuZZ[173] + 128./27.*aryuZZ[236] + 1./3.*
   aryuZZ[171] + 32./3.*aryuZZ[222] - aryuZZ[98] - 328./9.*aryuZZ[34];
   aryuZZ[171]=aryuZZ[37]*aryuZZ[171];
   aryuZZ[173]=167./3. + 29*aryuZZ[65];
   aryuZZ[173]=aryuZZ[108]*aryuZZ[173];
   aryuZZ[222]= - aryuZZ[34] - aryuZZ[99];
   aryuZZ[236]= - 1./3.*aryuZZ[34] + aryuZZ[99];
   aryuZZ[236]=aryuZZ[35]*aryuZZ[236];
   aryuZZ[222]=1./3.*aryuZZ[222] + aryuZZ[236];
   aryuZZ[222]=aryuZZ[18]*aryuZZ[222];
   aryuZZ[236]=aryuZZ[35]*MMZ;
   aryuZZ[248]=aryuZZ[6]*MMZ;
   aryuZZ[236]=8*aryuZZ[37] + aryuZZ[248] - 5./3.*MMZ + 4*aryuZZ[236];
   aryuZZ[236]=aryuZZ[3]*aryuZZ[236];
   aryuZZ[248]=1 + aryuZZ[35];
   aryuZZ[249]=aryuZZ[3]*aryuZZ[248];
   aryuZZ[251]=2./3.*aryuZZ[57];
   aryuZZ[249]=64./9.*aryuZZ[249] + aryuZZ[251] - 2./3.*aryuZZ[51] + 
   256./81.*aryuZZ[58] + aryuZZ[52];
   aryuZZ[249]=MMt*aryuZZ[249];
   aryuZZ[254]= - 23./9.*aryuZZ[50] - 16./9.*aryuZZ[41] + 32./9.*
   aryuZZ[69] + 2*aryuZZ[70] - aryuZZ[42];
   aryuZZ[254]=aryuZZ[34]*aryuZZ[254];
   aryuZZ[256]=29./3. + 10*aryuZZ[108];
   aryuZZ[260]=aryuZZ[13]*aryuZZ[256];
   aryuZZ[261]=64./9.*aryuZZ[21] + 14*aryuZZ[50] + 32./9.*aryuZZ[41] + 
   aryuZZ[267] - 32./3.*aryuZZ[69];
   aryuZZ[261]=aryuZZ[99]*aryuZZ[261];
   aryuZZ[262]= - 1 - aryuZZ[13];
   aryuZZ[263]=aryuZZ[36]*aryuZZ[262];
   aryuZZ[264]=32./3.*aryuZZ[99] - aryuZZ[98] - 16./3.*aryuZZ[34];
   aryuZZ[264]=aryuZZ[19]*aryuZZ[264];
   aryuZZ[265]= - 4*aryuZZ[35];
   aryuZZ[266]=5./3. + aryuZZ[265];
   aryuZZ[267]=aryuZZ[1]*aryuZZ[266];
   aryuZZ[266]=aryuZZ[39]*aryuZZ[266];
   aryuZZ[102]=aryuZZ[148] + aryuZZ[102] + aryuZZ[107] + 2*aryuZZ[249]
    + 16./9.*aryuZZ[236] + aryuZZ[171] + 1./3.*aryuZZ[110] + 128./9.*
   aryuZZ[222] + 160./243.*aryuZZ[266] + 32./81.*aryuZZ[267] + 8./3.*
   aryuZZ[144] + aryuZZ[142] + aryuZZ[264] + 6*aryuZZ[263] + 16./3.*
   aryuZZ[261] + 16./27.*aryuZZ[260] + 128./27.*aryuZZ[209] + 8./3.*
   aryuZZ[254] + 1./3.*aryuZZ[208] + 2./9.*aryuZZ[173] + 1./3.*
   aryuZZ[207] + 2./3.*aryuZZ[206] - 8./3.*aryuZZ[48] - 64./27.*
   aryuZZ[46] + 20./81.*aryuZZ[177] - 10*aryuZZ[16] - 520./81.*
   aryuZZ[49] + 5./3.*aryuZZ[66] + 107./9.*aryuZZ[65] + 640./81.*
   aryuZZ[62] + 8./9.*aryuZZ[63] - 80./81.*aryuZZ[31] + 2687./162. + 10
   *aryuZZ[60];
   aryuZZ[102]=aryuZZ[33]*aryuZZ[102];
   aryuZZ[107]= - 2*aryuZZ[22];
   aryuZZ[110]=aryuZZ[107] - 8*aryuZZ[50] + 2*aryuZZ[42] + 3*aryuZZ[40]
   ;
   aryuZZ[121]= - 1 + aryuZZ[121];
   aryuZZ[121]=aryuZZ[19]*aryuZZ[121];
   aryuZZ[142]= - 7./3. + aryuZZ[158];
   aryuZZ[142]=aryuZZ[17]*aryuZZ[142];
   aryuZZ[144]= - 1 - aryuZZ[38];
   aryuZZ[148]=aryuZZ[18]*aryuZZ[144];
   aryuZZ[171]= - 8*aryuZZ[44] - 3 - 2*aryuZZ[48];
   aryuZZ[171]=aryuZZ[37]*aryuZZ[171];
   aryuZZ[173]= - aryuZZ[35]*aryuZZ[17];
   aryuZZ[110]=3*aryuZZ[171] + 6*aryuZZ[148] + 11./3.*aryuZZ[173] + 
   aryuZZ[142] + 3*aryuZZ[121] + 3*aryuZZ[110] - 38./3.*aryuZZ[20];
   aryuZZ[110]=aryuZZ[3]*aryuZZ[110];
   aryuZZ[121]=aryuZZ[1]*aryuZZ[36];
   aryuZZ[142]=aryuZZ[39]*aryuZZ[36];
   aryuZZ[148]=1./3.*aryuZZ[142];
   aryuZZ[171]=aryuZZ[148] + aryuZZ[121] - aryuZZ[13] + aryuZZ[154];
   aryuZZ[171]=1./3.*aryuZZ[6]*aryuZZ[171];
   aryuZZ[206]=15*aryuZZ[60];
   aryuZZ[207]=11./2.*aryuZZ[63];
   aryuZZ[208]=aryuZZ[207] - 7963./324. + aryuZZ[206];
   aryuZZ[209]= - aryuZZ[17]*aryuZZ[38];
   aryuZZ[222]=aryuZZ[209] - aryuZZ[20] + aryuZZ[40] + aryuZZ[193];
   aryuZZ[222]=3*aryuZZ[101]*aryuZZ[222];
   aryuZZ[236]= - 4*aryuZZ[44];
   aryuZZ[249]=3 - aryuZZ[48];
   aryuZZ[254]=aryuZZ[249] + aryuZZ[236];
   aryuZZ[254]=aryuZZ[3]*aryuZZ[254];
   aryuZZ[233]=2*aryuZZ[57] - 11*aryuZZ[54] + aryuZZ[52] + aryuZZ[233];
   aryuZZ[233]=3*aryuZZ[254] + 1./3.*aryuZZ[233] + 6*aryuZZ[101];
   aryuZZ[233]=MMt*aryuZZ[233];
   aryuZZ[254]= - 4*aryuZZ[72];
   aryuZZ[260]= - 2*aryuZZ[62];
   aryuZZ[261]=11./3.*aryuZZ[61];
   aryuZZ[263]= - 73./24.*aryuZZ[49];
   aryuZZ[264]= - 15./4.*aryuZZ[16];
   aryuZZ[266]= - 11./3.*aryuZZ[64];
   aryuZZ[267]= - 145./36.*aryuZZ[67];
   aryuZZ[269]=1./8.*aryuZZ[48];
   aryuZZ[274]=7./2.*aryuZZ[44];
   aryuZZ[275]=4*aryuZZ[15];
   aryuZZ[276]= - 55./18. + 4*aryuZZ[13];
   aryuZZ[276]=aryuZZ[36]*aryuZZ[276];
   aryuZZ[178]=aryuZZ[178] - 2./3.*aryuZZ[13] - 11./2.*aryuZZ[44] - 46./
   9. - 11./8.*aryuZZ[48];
   aryuZZ[178]=aryuZZ[35]*aryuZZ[178];
   aryuZZ[278]= - aryuZZ[1]*aryuZZ[36];
   aryuZZ[279]=1./9.*aryuZZ[278];
   aryuZZ[281]= - aryuZZ[39]*aryuZZ[36];
   aryuZZ[282]=1./27.*aryuZZ[281];
   aryuZZ[283]=3*aryuZZ[38];
   aryuZZ[284]=1 + aryuZZ[283];
   aryuZZ[285]=2*aryuZZ[10]*aryuZZ[284];
   aryuZZ[286]=1 + aryuZZ[38];
   aryuZZ[286]=aryuZZ[18]*aryuZZ[101]*aryuZZ[286];
   aryuZZ[287]=6*aryuZZ[286];
   aryuZZ[288]= - aryuZZ[37]*aryuZZ[101];
   aryuZZ[289]=6*aryuZZ[288];
   aryuZZ[290]=9./4.*aryuZZ[54] + 4*aryuZZ[53];
   aryuZZ[290]=MMH*aryuZZ[290];
   aryuZZ[208]=aryuZZ[233] + aryuZZ[290] + aryuZZ[110] + aryuZZ[289] + 
   aryuZZ[171] + aryuZZ[287] + aryuZZ[285] + aryuZZ[282] + aryuZZ[279]
    + aryuZZ[178] + aryuZZ[222] + aryuZZ[276] + aryuZZ[234] + 
   aryuZZ[275] + aryuZZ[274] + aryuZZ[269] + aryuZZ[267] + aryuZZ[266]
    + aryuZZ[264] + aryuZZ[263] + 439./144.*aryuZZ[66] + aryuZZ[261] + 
   67./9.*aryuZZ[65] + aryuZZ[260] + 1./4.*aryuZZ[208] + aryuZZ[254];
   aryuZZ[208]=MMt*aryuZZ[208];
   aryuZZ[291]= - 5*aryuZZ[48];
   aryuZZ[292]=5381./81. + aryuZZ[291];
   aryuZZ[293]=7*aryuZZ[44];
   aryuZZ[294]=27*aryuZZ[43];
   aryuZZ[292]=aryuZZ[293] + 1./4.*aryuZZ[292] + aryuZZ[294];
   aryuZZ[295]=27./4.*aryuZZ[13];
   aryuZZ[296]= - 2*aryuZZ[11];
   aryuZZ[297]= - 35./18.*aryuZZ[35];
   aryuZZ[298]= - 1./9.*aryuZZ[1];
   aryuZZ[299]= - 1./27.*aryuZZ[39];
   aryuZZ[300]=4*aryuZZ[10];
   aryuZZ[301]=6*aryuZZ[117];
   aryuZZ[302]=aryuZZ[105] + 1 + aryuZZ[1];
   aryuZZ[302]=1./3.*aryuZZ[6]*aryuZZ[302];
   aryuZZ[303]=aryuZZ[37]*aryuZZ[101];
   aryuZZ[304]=3*aryuZZ[303];
   aryuZZ[305]=3*aryuZZ[47];
   aryuZZ[306]=3./2.*aryuZZ[45];
   aryuZZ[292]=aryuZZ[304] + aryuZZ[302] + aryuZZ[301] + aryuZZ[300] + 
   aryuZZ[299] + aryuZZ[298] + aryuZZ[297] + aryuZZ[296] + aryuZZ[295]
    + aryuZZ[306] + 1./2.*aryuZZ[292] + aryuZZ[305];
   aryuZZ[292]=aryuZZ[37]*aryuZZ[292];
   aryuZZ[307]= - 23./6.*aryuZZ[15];
   aryuZZ[308]=3*aryuZZ[72];
   aryuZZ[309]=11./3.*aryuZZ[64];
   aryuZZ[310]=aryuZZ[283] + aryuZZ[307] + aryuZZ[143] + 7./3.*
   aryuZZ[67] + aryuZZ[309] + 5./6.*aryuZZ[59] - 19./12.*aryuZZ[49] + 
   19./12.*aryuZZ[61] - 353./216. + aryuZZ[308];
   aryuZZ[311]=73./27. + 11./2.*aryuZZ[44];
   aryuZZ[311]=1./2.*aryuZZ[311] + 7./27.*aryuZZ[11];
   aryuZZ[312]=1./3.*aryuZZ[35];
   aryuZZ[311]=1./2.*aryuZZ[311] + aryuZZ[312];
   aryuZZ[311]=aryuZZ[35]*aryuZZ[311];
   aryuZZ[313]=59./18.*aryuZZ[35] + aryuZZ[161] + 7./9. + aryuZZ[158];
   aryuZZ[313]=aryuZZ[10]*aryuZZ[313];
   aryuZZ[314]=17./27.*aryuZZ[11];
   aryuZZ[315]= - 1./12.*aryuZZ[36];
   aryuZZ[316]= - MMH*aryuZZ[54];
   aryuZZ[310]=1./3.*aryuZZ[316] + 1./2.*aryuZZ[313] + aryuZZ[311] + 
   aryuZZ[315] + 1./2.*aryuZZ[310] + aryuZZ[314];
   aryuZZ[310]=MMH*aryuZZ[310];
   aryuZZ[311]=11./3.*aryuZZ[17];
   aryuZZ[181]=aryuZZ[181] + aryuZZ[311];
   aryuZZ[181]=aryuZZ[35]*aryuZZ[181];
   aryuZZ[313]=aryuZZ[229] - 11./3. + aryuZZ[283];
   aryuZZ[313]=1./2.*aryuZZ[17]*aryuZZ[313];
   aryuZZ[316]= - 151./16.*aryuZZ[42] - 4*aryuZZ[71] - 395./8.*
   aryuZZ[70];
   aryuZZ[317]=13*aryuZZ[36];
   aryuZZ[318]=1229./6. + aryuZZ[317];
   aryuZZ[318]=aryuZZ[19]*aryuZZ[318];
   aryuZZ[319]=161./27.*aryuZZ[35] + 673./27. + aryuZZ[229];
   aryuZZ[319]=aryuZZ[18]*aryuZZ[319];
   aryuZZ[320]= - 67./36.*aryuZZ[68];
   aryuZZ[321]=7./3.*aryuZZ[22];
   aryuZZ[322]=107./18.*aryuZZ[20];
   aryuZZ[323]= - aryuZZ[6]*aryuZZ[19];
   aryuZZ[324]=1./3.*aryuZZ[323];
   aryuZZ[325]= - 3*aryuZZ[37] - 4*aryuZZ[18] - 9./2.*aryuZZ[19] + 
   aryuZZ[17];
   aryuZZ[325]=3*aryuZZ[3]*aryuZZ[37]*aryuZZ[325];
   aryuZZ[326]= - 9./4.*aryuZZ[40];
   aryuZZ[181]=aryuZZ[208] + aryuZZ[310] + aryuZZ[325] + aryuZZ[292] + 
   aryuZZ[324] + 1./6.*aryuZZ[319] + 1./3.*aryuZZ[181] + aryuZZ[313] + 
   1./12.*aryuZZ[318] - 349./324.*aryuZZ[21] + aryuZZ[322] + 
   aryuZZ[321] + 2027./162.*aryuZZ[50] + 673./324.*aryuZZ[41] + 
   aryuZZ[326] + aryuZZ[320] + 1./9.*aryuZZ[316] - aryuZZ[69];
   aryuZZ[181]=MMt*aryuZZ[181];
   aryuZZ[208]=1./4.*aryuZZ[36];
   aryuZZ[292]=aryuZZ[144] + aryuZZ[208];
   aryuZZ[292]=aryuZZ[18]*aryuZZ[292];
   aryuZZ[316]= - 1 - aryuZZ[36];
   aryuZZ[318]=aryuZZ[19]*aryuZZ[316];
   aryuZZ[319]= - 61./9. + aryuZZ[158];
   aryuZZ[319]=aryuZZ[17]*aryuZZ[319];
   aryuZZ[327]=aryuZZ[35]*aryuZZ[17];
   aryuZZ[328]=3./2.*aryuZZ[40];
   aryuZZ[329]= - aryuZZ[22] - 4*aryuZZ[50] + aryuZZ[42] + aryuZZ[328];
   aryuZZ[329]=3*aryuZZ[329];
   aryuZZ[236]=aryuZZ[236] - 3./2. - aryuZZ[48];
   aryuZZ[236]=3*aryuZZ[37]*aryuZZ[236];
   aryuZZ[292]=aryuZZ[236] + 3*aryuZZ[292] + 7./18.*aryuZZ[327] + 1./2.
   *aryuZZ[319] + 3./2.*aryuZZ[318] + aryuZZ[329] - 37./9.*aryuZZ[20];
   aryuZZ[292]=aryuZZ[3]*aryuZZ[292];
   aryuZZ[318]=aryuZZ[240] + aryuZZ[184];
   aryuZZ[318]=11./18.*aryuZZ[142] + 1./3.*aryuZZ[318] + 1./2.*
   aryuZZ[121];
   aryuZZ[318]=aryuZZ[6]*aryuZZ[318];
   aryuZZ[319]=7./3.*aryuZZ[54] + aryuZZ[52] + aryuZZ[51];
   aryuZZ[251]=1./2.*aryuZZ[319] + aryuZZ[251];
   aryuZZ[319]=3*aryuZZ[101];
   aryuZZ[249]=1./2.*aryuZZ[249] - 2*aryuZZ[44];
   aryuZZ[249]=3*aryuZZ[3]*aryuZZ[249];
   aryuZZ[251]=aryuZZ[249] + 1./3.*aryuZZ[251] + aryuZZ[319];
   aryuZZ[251]=MMt*aryuZZ[251];
   aryuZZ[330]=1./2.*aryuZZ[40];
   aryuZZ[209]=1./2.*aryuZZ[209] + aryuZZ[198] + aryuZZ[330] - 
   aryuZZ[50];
   aryuZZ[209]=3*aryuZZ[101]*aryuZZ[209];
   aryuZZ[331]=aryuZZ[284] + aryuZZ[229];
   aryuZZ[331]=aryuZZ[10]*aryuZZ[331];
   aryuZZ[332]= - 1433./324. - 3*aryuZZ[60];
   aryuZZ[332]=5*aryuZZ[332] - 77./36.*aryuZZ[63];
   aryuZZ[333]=1./2.*aryuZZ[13];
   aryuZZ[334]= - 7./9. + aryuZZ[333];
   aryuZZ[334]=1./4.*aryuZZ[36]*aryuZZ[334];
   aryuZZ[335]=4./3.*aryuZZ[35] + aryuZZ[333] + 7./4.*aryuZZ[44] - 13./
   3. + 7./16.*aryuZZ[48];
   aryuZZ[335]=aryuZZ[35]*aryuZZ[335];
   aryuZZ[336]= - 13./24.*aryuZZ[54] + aryuZZ[170];
   aryuZZ[336]=MMH*aryuZZ[336];
   aryuZZ[337]= - 2*aryuZZ[72];
   aryuZZ[338]=1./2.*aryuZZ[65];
   aryuZZ[339]=2*aryuZZ[15];
   aryuZZ[286]=3*aryuZZ[286];
   aryuZZ[340]=3*aryuZZ[288];
   aryuZZ[251]=aryuZZ[251] + aryuZZ[336] + aryuZZ[292] + aryuZZ[340] + 
   aryuZZ[318] + aryuZZ[286] + aryuZZ[331] + 11./54.*aryuZZ[281] + 1./6.
   *aryuZZ[278] + 1./3.*aryuZZ[335] + aryuZZ[209] + aryuZZ[334] + 
   aryuZZ[158] + aryuZZ[339] + 61./12.*aryuZZ[44] + 43./48.*aryuZZ[48]
    - 595./216.*aryuZZ[67] + 7./18.*aryuZZ[64] + 15./4.*aryuZZ[16] + 
   277./144.*aryuZZ[49] - 967./288.*aryuZZ[66] - 7./18.*aryuZZ[61] + 
   aryuZZ[338] - aryuZZ[62] + 1./4.*aryuZZ[332] + aryuZZ[337];
   aryuZZ[251]=MMt*aryuZZ[251];
   aryuZZ[292]= - 47*aryuZZ[70] - 275./6.*aryuZZ[42];
   aryuZZ[292]= - 661./36.*aryuZZ[68] + 1./8.*aryuZZ[292] + 7*
   aryuZZ[69];
   aryuZZ[318]= - 983./36. + aryuZZ[36];
   aryuZZ[318]=aryuZZ[19]*aryuZZ[318];
   aryuZZ[331]=aryuZZ[36] + 515./27. + aryuZZ[283];
   aryuZZ[331]=aryuZZ[17]*aryuZZ[331];
   aryuZZ[332]=11*aryuZZ[17];
   aryuZZ[335]= - 205./3.*aryuZZ[19] + aryuZZ[332];
   aryuZZ[335]=aryuZZ[35]*aryuZZ[335];
   aryuZZ[336]=71./54.*aryuZZ[35] + 173./27. + aryuZZ[36];
   aryuZZ[336]=aryuZZ[18]*aryuZZ[336];
   aryuZZ[292]=1./9.*aryuZZ[323] + 1./6.*aryuZZ[336] + 1./36.*
   aryuZZ[335] + 1./2.*aryuZZ[331] + 1./6.*aryuZZ[318] - 869./324.*
   aryuZZ[21] + 251./54.*aryuZZ[20] + 61./36.*aryuZZ[22] + 1957./648.*
   aryuZZ[50] + 113./324.*aryuZZ[41] + 1./3.*aryuZZ[292] + aryuZZ[326];
   aryuZZ[318]=aryuZZ[283] - 29./18.*aryuZZ[15] - 43./12.*aryuZZ[44] + 
   31./9.*aryuZZ[67] - 67./9.*aryuZZ[64] - 25./18.*aryuZZ[59] - 97./36.
   *aryuZZ[49] + 97./36.*aryuZZ[61] - 1897./216. + aryuZZ[308];
   aryuZZ[143]= - 1./9. + aryuZZ[143];
   aryuZZ[331]=1./9.*aryuZZ[11];
   aryuZZ[143]=1./2.*aryuZZ[143] + aryuZZ[331];
   aryuZZ[335]= - 1./3.*aryuZZ[35];
   aryuZZ[143]=7./2.*aryuZZ[143] + aryuZZ[335];
   aryuZZ[143]=aryuZZ[35]*aryuZZ[143];
   aryuZZ[336]= - 229./27. + aryuZZ[158];
   aryuZZ[336]= - 131./108.*aryuZZ[35] + 1./2.*aryuZZ[336] + 
   aryuZZ[161];
   aryuZZ[336]=aryuZZ[10]*aryuZZ[336];
   aryuZZ[341]=MMH*aryuZZ[54];
   aryuZZ[143]=1./9.*aryuZZ[341] + 1./2.*aryuZZ[336] + 1./3.*
   aryuZZ[143] + aryuZZ[315] + 1./4.*aryuZZ[318] + aryuZZ[314];
   aryuZZ[143]=MMH*aryuZZ[143];
   aryuZZ[315]=2659./27. + aryuZZ[238];
   aryuZZ[315]=61./3.*aryuZZ[44] + 1./3.*aryuZZ[315] + aryuZZ[294];
   aryuZZ[315]= - 29./4.*aryuZZ[13] + aryuZZ[306] + 1./2.*aryuZZ[315]
    + aryuZZ[305];
   aryuZZ[318]=aryuZZ[230] + 1./9. + aryuZZ[1];
   aryuZZ[318]=aryuZZ[6]*aryuZZ[318];
   aryuZZ[117]=3*aryuZZ[117];
   aryuZZ[336]=3./2.*aryuZZ[303];
   aryuZZ[315]=aryuZZ[336] + 1./2.*aryuZZ[318] + aryuZZ[117] + 7./2.*
   aryuZZ[10] - 11./54.*aryuZZ[39] - 1./6.*aryuZZ[1] - 11./36.*
   aryuZZ[35] + 1./2.*aryuZZ[315] + aryuZZ[296];
   aryuZZ[315]=aryuZZ[37]*aryuZZ[315];
   aryuZZ[318]=6*aryuZZ[19] + aryuZZ[119];
   aryuZZ[318]= - 9./2.*aryuZZ[37] + 3*aryuZZ[318] + 11./4.*aryuZZ[18];
   aryuZZ[318]=aryuZZ[3]*aryuZZ[37]*aryuZZ[318];
   aryuZZ[143]=aryuZZ[251] + aryuZZ[143] + aryuZZ[318] + 1./2.*
   aryuZZ[292] + aryuZZ[315];
   aryuZZ[143]=MMt*aryuZZ[143];
   aryuZZ[251]=971./8. - 5*aryuZZ[14];
   aryuZZ[251]= - 25./8.*aryuZZ[15] + 25./8.*aryuZZ[67] + 523./12.*
   aryuZZ[64] + 1./3.*aryuZZ[251] + 25./8.*aryuZZ[59];
   aryuZZ[292]=11./9.*aryuZZ[11];
   aryuZZ[315]= - 49./48. + 17*aryuZZ[35];
   aryuZZ[315]=aryuZZ[10]*aryuZZ[315];
   aryuZZ[318]=aryuZZ[118] + aryuZZ[10];
   aryuZZ[341]=aryuZZ[6]*aryuZZ[318];
   aryuZZ[342]=aryuZZ[99]*aryuZZ[50];
   aryuZZ[343]= - aryuZZ[99] + aryuZZ[175];
   aryuZZ[343]=aryuZZ[37]*aryuZZ[343];
   aryuZZ[344]= - aryuZZ[35]*aryuZZ[11];
   aryuZZ[251]=25./4.*aryuZZ[343] + 5./3.*aryuZZ[341] + 1./3.*
   aryuZZ[315] + 17./6.*aryuZZ[344] + 25./4.*aryuZZ[342] + 1./2.*
   aryuZZ[251] + aryuZZ[292];
   aryuZZ[251]=MMH*aryuZZ[251];
   aryuZZ[315]=25./12.*aryuZZ[21];
   aryuZZ[341]=aryuZZ[315] - 227./18.*aryuZZ[50] - 25./12.*aryuZZ[41]
    + 455./9.*aryuZZ[68] - 43./2.*aryuZZ[40];
   aryuZZ[345]=49./48. + aryuZZ[199];
   aryuZZ[345]=aryuZZ[18]*aryuZZ[345];
   aryuZZ[346]=aryuZZ[19] + aryuZZ[17];
   aryuZZ[347]=1./2.*aryuZZ[346] - aryuZZ[18];
   aryuZZ[347]=aryuZZ[6]*aryuZZ[347];
   aryuZZ[348]=aryuZZ[35]*aryuZZ[346];
   aryuZZ[341]=5./9.*aryuZZ[347] + 1./9.*aryuZZ[345] + 17./18.*
   aryuZZ[348] - 1511./108.*aryuZZ[17] + 1./4.*aryuZZ[341] - 11./27.*
   aryuZZ[19];
   aryuZZ[345]= - 1843./18. - 25*aryuZZ[44];
   aryuZZ[347]=17./9.*aryuZZ[11];
   aryuZZ[349]= - aryuZZ[18]*aryuZZ[99];
   aryuZZ[350]=aryuZZ[37]*aryuZZ[99];
   aryuZZ[345]=25./48.*aryuZZ[350] + 25./48.*aryuZZ[349] - 175./36.*
   aryuZZ[10] + 1./16.*aryuZZ[345] + aryuZZ[347];
   aryuZZ[345]=aryuZZ[37]*aryuZZ[345];
   aryuZZ[251]=1./12.*aryuZZ[251] + 1./4.*aryuZZ[341] + aryuZZ[345];
   aryuZZ[251]=MMH*aryuZZ[251];
   aryuZZ[341]=1205./24.*aryuZZ[37] + 427./24.*aryuZZ[18] - 149*
   aryuZZ[19] - 91./8.*aryuZZ[17];
   aryuZZ[341]=aryuZZ[37]*aryuZZ[341];
   aryuZZ[251]=1./9.*aryuZZ[341] + aryuZZ[251];
   aryuZZ[341]=aryuZZ[121] + 11./9.*aryuZZ[142];
   aryuZZ[341]=aryuZZ[6]*aryuZZ[341];
   aryuZZ[345]= - 1./3. + aryuZZ[333];
   aryuZZ[345]=aryuZZ[36]*aryuZZ[345];
   aryuZZ[351]=1./3.*aryuZZ[345] + aryuZZ[16] - aryuZZ[66] - 1./4. - 
   aryuZZ[60];
   aryuZZ[352]=aryuZZ[18]*aryuZZ[36];
   aryuZZ[353]=aryuZZ[3]*aryuZZ[352];
   aryuZZ[278]=1./3.*aryuZZ[278];
   aryuZZ[354]=aryuZZ[10]*aryuZZ[36];
   aryuZZ[341]=3./2.*aryuZZ[353] + aryuZZ[341] + aryuZZ[354] + 11./27.*
   aryuZZ[281] + 1./2.*aryuZZ[351] + aryuZZ[278];
   aryuZZ[341]=MMt*aryuZZ[341];
   aryuZZ[351]= - 11./18.*aryuZZ[39] - 1./2.*aryuZZ[1] + 2./3. - 5./8.*
   aryuZZ[13];
   aryuZZ[230]=aryuZZ[1] + aryuZZ[230];
   aryuZZ[353]=aryuZZ[6]*aryuZZ[230];
   aryuZZ[355]=1./2.*aryuZZ[10];
   aryuZZ[351]=1./2.*aryuZZ[353] + 1./3.*aryuZZ[351] + aryuZZ[355];
   aryuZZ[351]=aryuZZ[37]*aryuZZ[351];
   aryuZZ[353]=aryuZZ[35]*aryuZZ[19];
   aryuZZ[205]= - 17 + aryuZZ[205];
   aryuZZ[205]=aryuZZ[18]*aryuZZ[205];
   aryuZZ[205]=1./27.*aryuZZ[205] + 7./54.*aryuZZ[353] - aryuZZ[50] + 
   17./27.*aryuZZ[19];
   aryuZZ[356]= - aryuZZ[19] + 5./2.*aryuZZ[18];
   aryuZZ[356]=aryuZZ[3]*aryuZZ[37]*aryuZZ[356];
   aryuZZ[205]=1./2.*aryuZZ[341] + 1./2.*aryuZZ[356] + 1./4.*
   aryuZZ[205] + aryuZZ[351];
   aryuZZ[205]=MMt*aryuZZ[205];
   aryuZZ[341]=11./3.*aryuZZ[11];
   aryuZZ[351]= - aryuZZ[6]*aryuZZ[11];
   aryuZZ[351]=5./2.*aryuZZ[351] + aryuZZ[341] + 17./2.*aryuZZ[344];
   aryuZZ[351]=MMH*aryuZZ[351];
   aryuZZ[356]=aryuZZ[19] - aryuZZ[17];
   aryuZZ[357]=aryuZZ[6]*aryuZZ[356];
   aryuZZ[358]=aryuZZ[35]*aryuZZ[356];
   aryuZZ[359]= - aryuZZ[19] + aryuZZ[17];
   aryuZZ[351]=aryuZZ[351] + 5./2.*aryuZZ[357] + 11./3.*aryuZZ[359] + 
   17./2.*aryuZZ[358];
   aryuZZ[351]=MMH*aryuZZ[351];
   aryuZZ[357]=aryuZZ[19] - aryuZZ[18];
   aryuZZ[360]=aryuZZ[37]*aryuZZ[357];
   aryuZZ[351]=17*aryuZZ[360] + aryuZZ[351];
   aryuZZ[205]=1./108.*aryuZZ[351] + aryuZZ[205];
   aryuZZ[205]=aryuZZ[272]*aryuZZ[205];
   aryuZZ[143]=aryuZZ[205] + 1./3.*aryuZZ[251] + aryuZZ[143];
   aryuZZ[143]=aryuZZ[272]*aryuZZ[143];
   aryuZZ[205]=5*aryuZZ[44];
   aryuZZ[251]= - 299./54. + aryuZZ[205];
   aryuZZ[351]=aryuZZ[18]*aryuZZ[99];
   aryuZZ[361]= - aryuZZ[37]*aryuZZ[99];
   aryuZZ[362]=25./18.*aryuZZ[10];
   aryuZZ[251]=5./24.*aryuZZ[361] + 5./24.*aryuZZ[351] + aryuZZ[362] + 
   1./8.*aryuZZ[251] + aryuZZ[314];
   aryuZZ[251]=aryuZZ[37]*aryuZZ[251];
   aryuZZ[361]=17./8. - aryuZZ[14];
   aryuZZ[176]=aryuZZ[280] - 5./8.*aryuZZ[67] + 1./12.*aryuZZ[64] + 1./
   3.*aryuZZ[361] + aryuZZ[176];
   aryuZZ[280]= - aryuZZ[99]*aryuZZ[50];
   aryuZZ[176]=5./4.*aryuZZ[280] + 1./2.*aryuZZ[176] + 1./27.*
   aryuZZ[11];
   aryuZZ[252]=aryuZZ[252] + aryuZZ[335];
   aryuZZ[252]=aryuZZ[10]*aryuZZ[252];
   aryuZZ[280]=aryuZZ[216] + aryuZZ[355];
   aryuZZ[280]=aryuZZ[6]*aryuZZ[280];
   aryuZZ[355]=aryuZZ[99] + aryuZZ[225];
   aryuZZ[355]=aryuZZ[37]*aryuZZ[355];
   aryuZZ[176]=5./8.*aryuZZ[355] + 1./6.*aryuZZ[280] + 1./4.*
   aryuZZ[252] + 1./2.*aryuZZ[176] + 1./9.*aryuZZ[344];
   aryuZZ[176]=MMH*aryuZZ[176];
   aryuZZ[252]= - 1./2.*aryuZZ[40];
   aryuZZ[280]= - 5./12.*aryuZZ[21] + 151./18.*aryuZZ[50] + 5./12.*
   aryuZZ[41] + 5./9.*aryuZZ[68] + aryuZZ[252];
   aryuZZ[355]= - 5./54. - aryuZZ[36];
   aryuZZ[355]=aryuZZ[17]*aryuZZ[355];
   aryuZZ[280]=1./6.*aryuZZ[355] + 1./4.*aryuZZ[280] - 1./81.*
   aryuZZ[19];
   aryuZZ[355]= - 23./24. + aryuZZ[36];
   aryuZZ[361]=aryuZZ[355] + aryuZZ[312];
   aryuZZ[361]=aryuZZ[18]*aryuZZ[361];
   aryuZZ[363]= - aryuZZ[19] + 5./2.*aryuZZ[17];
   aryuZZ[123]=1./3.*aryuZZ[363] + aryuZZ[123];
   aryuZZ[123]=aryuZZ[6]*aryuZZ[123];
   aryuZZ[363]=aryuZZ[19] - 7./4.*aryuZZ[17];
   aryuZZ[363]=aryuZZ[35]*aryuZZ[363];
   aryuZZ[123]=1./3.*aryuZZ[176] + aryuZZ[251] + 1./18.*aryuZZ[123] + 1.
   /12.*aryuZZ[361] + 1./2.*aryuZZ[280] + 1./27.*aryuZZ[363];
   aryuZZ[123]=MMH*aryuZZ[123];
   aryuZZ[176]=385./27.*aryuZZ[37] + 1427./27.*aryuZZ[18] + 103*
   aryuZZ[19] - 53./3.*aryuZZ[17];
   aryuZZ[176]=aryuZZ[37]*aryuZZ[176];
   aryuZZ[143]=aryuZZ[143] + aryuZZ[181] + 1./12.*aryuZZ[176] + 
   aryuZZ[123];
   aryuZZ[143]=aryuZZ[272]*aryuZZ[143];
   aryuZZ[176]=aryuZZ[207] - 289./12. + aryuZZ[206];
   aryuZZ[110]=aryuZZ[233] + aryuZZ[290] + aryuZZ[110] + aryuZZ[289] + 
   aryuZZ[171] + aryuZZ[287] + aryuZZ[285] + aryuZZ[282] + aryuZZ[279]
    + aryuZZ[178] + aryuZZ[222] + aryuZZ[276] + aryuZZ[234] + 
   aryuZZ[275] + aryuZZ[274] + aryuZZ[269] + aryuZZ[267] + aryuZZ[266]
    + aryuZZ[264] + aryuZZ[263] + 301./48.*aryuZZ[66] + aryuZZ[261] + 
   32./3.*aryuZZ[65] + aryuZZ[260] + 1./4.*aryuZZ[176] + aryuZZ[254];
   aryuZZ[110]=MMt*aryuZZ[110];
   aryuZZ[171]=589./9. + aryuZZ[291];
   aryuZZ[171]=aryuZZ[293] + 1./4.*aryuZZ[171] + aryuZZ[294];
   aryuZZ[171]=aryuZZ[304] + aryuZZ[302] + aryuZZ[301] + aryuZZ[300] + 
   aryuZZ[299] + aryuZZ[298] + aryuZZ[297] + aryuZZ[296] + aryuZZ[295]
    + aryuZZ[306] + 1./2.*aryuZZ[171] + aryuZZ[305];
   aryuZZ[171]=aryuZZ[37]*aryuZZ[171];
   aryuZZ[176]= - 4./3.*aryuZZ[71];
   aryuZZ[178]= - 205./16.*aryuZZ[42] + aryuZZ[176] + 23./8.*aryuZZ[70]
   ;
   aryuZZ[181]=4855./18. + aryuZZ[317];
   aryuZZ[181]=aryuZZ[19]*aryuZZ[181];
   aryuZZ[206]=43./3.*aryuZZ[19] + aryuZZ[17];
   aryuZZ[206]=aryuZZ[35]*aryuZZ[206];
   aryuZZ[207]= - 13*aryuZZ[35] - 13 + aryuZZ[229];
   aryuZZ[207]=aryuZZ[18]*aryuZZ[207];
   aryuZZ[110]=aryuZZ[110] + aryuZZ[310] + aryuZZ[325] + aryuZZ[171] + 
   aryuZZ[324] + 1./6.*aryuZZ[207] + 11./9.*aryuZZ[206] + aryuZZ[313]
    + 1./12.*aryuZZ[181] + aryuZZ[315] + aryuZZ[322] + aryuZZ[321] + 
   113./6.*aryuZZ[50] - 13./12.*aryuZZ[41] + aryuZZ[326] + aryuZZ[320]
    + 1./3.*aryuZZ[178] - aryuZZ[69];
   aryuZZ[110]=MMt*aryuZZ[110];
   aryuZZ[171]= - 1./4.*aryuZZ[36];
   aryuZZ[178]=aryuZZ[144] + aryuZZ[171];
   aryuZZ[178]=aryuZZ[18]*aryuZZ[178];
   aryuZZ[181]= - 1./2. - aryuZZ[36];
   aryuZZ[181]=aryuZZ[19]*aryuZZ[181];
   aryuZZ[206]= - 5 + aryuZZ[158];
   aryuZZ[206]=aryuZZ[17]*aryuZZ[206];
   aryuZZ[173]=aryuZZ[236] + 3*aryuZZ[178] + 1./2.*aryuZZ[173] + 1./2.*
   aryuZZ[206] + 3*aryuZZ[181] + aryuZZ[329] - 5*aryuZZ[20];
   aryuZZ[173]=aryuZZ[3]*aryuZZ[173];
   aryuZZ[178]= - aryuZZ[54] - aryuZZ[52] + aryuZZ[51];
   aryuZZ[178]=aryuZZ[249] + 1./2.*aryuZZ[178] + aryuZZ[319];
   aryuZZ[178]=MMt*aryuZZ[178];
   aryuZZ[181]= - 1./2.*aryuZZ[36];
   aryuZZ[206]=aryuZZ[284] + aryuZZ[181];
   aryuZZ[206]=aryuZZ[10]*aryuZZ[206];
   aryuZZ[207]=aryuZZ[281] + aryuZZ[13] + aryuZZ[278];
   aryuZZ[207]=aryuZZ[6]*aryuZZ[207];
   aryuZZ[222]= - 5./8.*aryuZZ[63] - 223./12. + 13*aryuZZ[60];
   aryuZZ[233]=1./9. + aryuZZ[210];
   aryuZZ[233]=aryuZZ[36]*aryuZZ[233];
   aryuZZ[234]=aryuZZ[333] - 3./4.*aryuZZ[44] - 1 - 3./16.*aryuZZ[48];
   aryuZZ[234]=aryuZZ[35]*aryuZZ[234];
   aryuZZ[170]=1./8.*aryuZZ[54] + aryuZZ[170];
   aryuZZ[170]=MMH*aryuZZ[170];
   aryuZZ[142]=aryuZZ[178] + aryuZZ[170] + aryuZZ[173] + aryuZZ[340] + 
   1./2.*aryuZZ[207] + aryuZZ[286] + aryuZZ[206] + 1./6.*aryuZZ[142] + 
   1./18.*aryuZZ[121] + aryuZZ[234] + aryuZZ[209] + 5./4.*aryuZZ[233]
    + aryuZZ[158] + aryuZZ[339] + 15./4.*aryuZZ[44] + 9./16.*aryuZZ[48]
    - 59./24.*aryuZZ[67] + aryuZZ[244] - 13./2.*aryuZZ[16] + 13./16.*
   aryuZZ[49] + 365./32.*aryuZZ[66] + aryuZZ[218] + aryuZZ[338] - 
   aryuZZ[62] + 1./2.*aryuZZ[222] + aryuZZ[337];
   aryuZZ[142]=MMt*aryuZZ[142];
   aryuZZ[170]= - 161./9. + 3./2.*aryuZZ[48];
   aryuZZ[170]=15*aryuZZ[44] + 1./2.*aryuZZ[170] + aryuZZ[294];
   aryuZZ[173]=115./4.*aryuZZ[13];
   aryuZZ[104]=aryuZZ[10] + aryuZZ[105] + aryuZZ[104] - 3./2.*
   aryuZZ[35] + aryuZZ[173] + aryuZZ[306] + 1./2.*aryuZZ[170] + 
   aryuZZ[305];
   aryuZZ[105]= - aryuZZ[39] + 1 + aryuZZ[138];
   aryuZZ[105]=aryuZZ[6]*aryuZZ[105];
   aryuZZ[104]=aryuZZ[336] + 1./2.*aryuZZ[105] + 1./2.*aryuZZ[104] + 
   aryuZZ[117];
   aryuZZ[104]=aryuZZ[37]*aryuZZ[104];
   aryuZZ[105]= - 5./2.*aryuZZ[15];
   aryuZZ[117]=aryuZZ[283] + aryuZZ[105] - 9./4.*aryuZZ[44] + 3*
   aryuZZ[67] + aryuZZ[155] + aryuZZ[152] - 9./4.*aryuZZ[49] + 9./4.*
   aryuZZ[61] - 563./72. + aryuZZ[308];
   aryuZZ[138]= - 4./9. + 3./16.*aryuZZ[44];
   aryuZZ[138]=aryuZZ[35]*aryuZZ[138];
   aryuZZ[152]= - 73./18.*aryuZZ[35] - 59./9. + aryuZZ[158];
   aryuZZ[152]=aryuZZ[10]*aryuZZ[152];
   aryuZZ[117]=1./4.*aryuZZ[152] + 1./4.*aryuZZ[117] + aryuZZ[138];
   aryuZZ[117]=MMH*aryuZZ[117];
   aryuZZ[138]= - 129*aryuZZ[70] + 157./2.*aryuZZ[42];
   aryuZZ[152]=1./4.*aryuZZ[41];
   aryuZZ[170]= - 5./4.*aryuZZ[21];
   aryuZZ[138]=aryuZZ[170] + 31./6.*aryuZZ[20] + 19./4.*aryuZZ[22] + 
   841./24.*aryuZZ[50] + aryuZZ[152] + aryuZZ[326] - 53./12.*aryuZZ[68]
    + 1./8.*aryuZZ[138] + aryuZZ[69];
   aryuZZ[178]= - 275./48. + 4*aryuZZ[36];
   aryuZZ[178]=aryuZZ[19]*aryuZZ[178];
   aryuZZ[206]=47./3. + aryuZZ[283];
   aryuZZ[206]=aryuZZ[17]*aryuZZ[206];
   aryuZZ[207]= - 115*aryuZZ[19] + 47*aryuZZ[17];
   aryuZZ[207]=aryuZZ[35]*aryuZZ[207];
   aryuZZ[209]=109 - 43./2.*aryuZZ[35];
   aryuZZ[209]=aryuZZ[18]*aryuZZ[209];
   aryuZZ[115]=aryuZZ[115] + aryuZZ[17];
   aryuZZ[115]= - 9*aryuZZ[37] + 3*aryuZZ[115] - 59./2.*aryuZZ[18];
   aryuZZ[115]=aryuZZ[3]*aryuZZ[37]*aryuZZ[115];
   aryuZZ[104]=aryuZZ[142] + aryuZZ[117] + 1./2.*aryuZZ[115] + 
   aryuZZ[104] + 1./2.*aryuZZ[323] + 1./36.*aryuZZ[209] + 1./24.*
   aryuZZ[207] + 1./4.*aryuZZ[206] + 1./2.*aryuZZ[138] + aryuZZ[178];
   aryuZZ[104]=MMt*aryuZZ[104];
   aryuZZ[115]=67./8. - aryuZZ[14];
   aryuZZ[115]= - 1./8.*aryuZZ[15] + 1./8.*aryuZZ[67] + 35./12.*
   aryuZZ[64] + 1./3.*aryuZZ[115] + 1./8.*aryuZZ[59];
   aryuZZ[117]=aryuZZ[35]*aryuZZ[11];
   aryuZZ[138]= - 1./9.*aryuZZ[11];
   aryuZZ[115]=1./6.*aryuZZ[117] + 1./4.*aryuZZ[342] + 1./2.*
   aryuZZ[115] + aryuZZ[138];
   aryuZZ[142]=1./4.*aryuZZ[202] + aryuZZ[312];
   aryuZZ[142]=aryuZZ[10]*aryuZZ[142];
   aryuZZ[178]=1./4.*aryuZZ[11] + 1./3.*aryuZZ[10];
   aryuZZ[178]=aryuZZ[6]*aryuZZ[178];
   aryuZZ[115]=1./16.*aryuZZ[343] + 1./6.*aryuZZ[178] + 1./4.*
   aryuZZ[115] + 1./3.*aryuZZ[142];
   aryuZZ[115]=MMH*aryuZZ[115];
   aryuZZ[142]=1./4.*aryuZZ[21];
   aryuZZ[178]=aryuZZ[142] + 5./6.*aryuZZ[50] + aryuZZ[190] + 31./3.*
   aryuZZ[68] - 9./2.*aryuZZ[40];
   aryuZZ[202]= - 107./12. - aryuZZ[36];
   aryuZZ[202]=aryuZZ[17]*aryuZZ[202];
   aryuZZ[206]= - aryuZZ[19] + aryuZZ[311];
   aryuZZ[207]=aryuZZ[35]*aryuZZ[206];
   aryuZZ[178]=1./6.*aryuZZ[207] + 1./3.*aryuZZ[202] + 1./4.*
   aryuZZ[178] + 1./9.*aryuZZ[19];
   aryuZZ[202]= - 545./18. - 3*aryuZZ[44];
   aryuZZ[202]=1./4.*aryuZZ[350] + 1./4.*aryuZZ[349] + 1./4.*
   aryuZZ[202] - 41./9.*aryuZZ[10];
   aryuZZ[202]=aryuZZ[37]*aryuZZ[202];
   aryuZZ[207]= - 23./48. + aryuZZ[36];
   aryuZZ[209]=1./4.*aryuZZ[207] + aryuZZ[335];
   aryuZZ[209]=aryuZZ[18]*aryuZZ[209];
   aryuZZ[218]= - aryuZZ[19] + 7./3.*aryuZZ[17];
   aryuZZ[222]= - 1./3.*aryuZZ[18];
   aryuZZ[218]=1./4.*aryuZZ[218] + aryuZZ[222];
   aryuZZ[218]=aryuZZ[6]*aryuZZ[218];
   aryuZZ[115]=aryuZZ[115] + 1./4.*aryuZZ[202] + 1./6.*aryuZZ[218] + 1./
   4.*aryuZZ[178] + 1./3.*aryuZZ[209];
   aryuZZ[115]=MMH*aryuZZ[115];
   aryuZZ[178]=aryuZZ[278] + aryuZZ[281];
   aryuZZ[178]=aryuZZ[6]*aryuZZ[178];
   aryuZZ[202]= - aryuZZ[19]*aryuZZ[36];
   aryuZZ[209]= - aryuZZ[18]*aryuZZ[36];
   aryuZZ[218]=aryuZZ[202] + 1./2.*aryuZZ[209];
   aryuZZ[218]=aryuZZ[3]*aryuZZ[218];
   aryuZZ[233]=59./9. + 33./2.*aryuZZ[13];
   aryuZZ[234]=aryuZZ[36]*aryuZZ[233];
   aryuZZ[236]= - aryuZZ[10]*aryuZZ[36];
   aryuZZ[121]=3*aryuZZ[218] + aryuZZ[178] + aryuZZ[236] + aryuZZ[148]
    + 1./2.*aryuZZ[234] + 1./9.*aryuZZ[121];
   aryuZZ[121]=MMt*aryuZZ[121];
   aryuZZ[148]=5*aryuZZ[36];
   aryuZZ[178]=17./3. + aryuZZ[148];
   aryuZZ[178]=aryuZZ[19]*aryuZZ[178];
   aryuZZ[218]= - 47./6.*aryuZZ[19] + aryuZZ[17];
   aryuZZ[218]=aryuZZ[35]*aryuZZ[218];
   aryuZZ[249]=19./2.*aryuZZ[35] - 17 - aryuZZ[36];
   aryuZZ[249]=aryuZZ[18]*aryuZZ[249];
   aryuZZ[251]= - aryuZZ[17]*aryuZZ[36];
   aryuZZ[178]=1./3.*aryuZZ[249] + aryuZZ[218] + aryuZZ[178] + 
   aryuZZ[251];
   aryuZZ[218]= - 1./2. + aryuZZ[11];
   aryuZZ[249]=aryuZZ[35]*aryuZZ[218];
   aryuZZ[254]= - aryuZZ[35] + 1 + aryuZZ[229];
   aryuZZ[260]=aryuZZ[10]*aryuZZ[254];
   aryuZZ[249]=aryuZZ[260] + 1./2.*aryuZZ[249] - aryuZZ[11] + 
   aryuZZ[208];
   aryuZZ[249]=MMH*aryuZZ[249];
   aryuZZ[233]=1./2.*aryuZZ[140] - 5./2.*aryuZZ[10] + 1./6.*aryuZZ[39]
    + 1./18.*aryuZZ[1] + 1./4.*aryuZZ[233] + 2*aryuZZ[11];
   aryuZZ[233]=aryuZZ[37]*aryuZZ[233];
   aryuZZ[260]= - 10*aryuZZ[19] + 31./4.*aryuZZ[18];
   aryuZZ[260]=aryuZZ[3]*aryuZZ[37]*aryuZZ[260];
   aryuZZ[121]=1./2.*aryuZZ[121] + 1./3.*aryuZZ[249] + aryuZZ[260] + 1./
   4.*aryuZZ[178] + aryuZZ[233];
   aryuZZ[121]=MMt*aryuZZ[121];
   aryuZZ[178]= - aryuZZ[19] - aryuZZ[17];
   aryuZZ[233]=aryuZZ[35]*aryuZZ[178];
   aryuZZ[249]=aryuZZ[35] - 1./3. - aryuZZ[36];
   aryuZZ[260]=aryuZZ[18]*aryuZZ[249];
   aryuZZ[261]=1./2.*aryuZZ[178] + aryuZZ[18];
   aryuZZ[261]=aryuZZ[6]*aryuZZ[261];
   aryuZZ[263]=1./3.*aryuZZ[19];
   aryuZZ[264]=aryuZZ[17]*aryuZZ[36];
   aryuZZ[233]=aryuZZ[261] + aryuZZ[260] + 1./2.*aryuZZ[233] + 
   aryuZZ[263] + aryuZZ[264];
   aryuZZ[260]=1./2.*aryuZZ[11];
   aryuZZ[261]=aryuZZ[260] - aryuZZ[10];
   aryuZZ[261]=aryuZZ[6]*aryuZZ[261];
   aryuZZ[117]=aryuZZ[261] + aryuZZ[150] + aryuZZ[151] + 1./2.*
   aryuZZ[117];
   aryuZZ[117]=MMH*aryuZZ[117];
   aryuZZ[150]= - aryuZZ[11] + aryuZZ[10];
   aryuZZ[261]=aryuZZ[37]*aryuZZ[150];
   aryuZZ[117]=1./4.*aryuZZ[117] + 1./4.*aryuZZ[233] + aryuZZ[261];
   aryuZZ[117]=MMH*aryuZZ[117];
   aryuZZ[117]=17./4.*aryuZZ[360] + aryuZZ[117];
   aryuZZ[117]=1./3.*aryuZZ[117] + aryuZZ[121];
   aryuZZ[117]=aryuZZ[133]*aryuZZ[117];
   aryuZZ[121]=7./2.*aryuZZ[17];
   aryuZZ[233]=23./6.*aryuZZ[37] + 227./18.*aryuZZ[18] - 89./3.*
   aryuZZ[19] + aryuZZ[121];
   aryuZZ[233]=aryuZZ[37]*aryuZZ[233];
   aryuZZ[104]=aryuZZ[117] + aryuZZ[104] + 1./4.*aryuZZ[233] + 
   aryuZZ[115];
   aryuZZ[104]=aryuZZ[133]*aryuZZ[104];
   aryuZZ[115]=2207./3.*aryuZZ[19] - 53*aryuZZ[17];
   aryuZZ[115]= - 71./3.*aryuZZ[37] + 1./3.*aryuZZ[115] - 23*aryuZZ[18]
   ;
   aryuZZ[115]=aryuZZ[37]*aryuZZ[115];
   aryuZZ[104]=aryuZZ[104] + aryuZZ[110] + 1./12.*aryuZZ[115] + 
   aryuZZ[123];
   aryuZZ[104]=aryuZZ[133]*aryuZZ[104];
   aryuZZ[110]=9*aryuZZ[38];
   aryuZZ[115]= - 16*aryuZZ[35];
   aryuZZ[117]=aryuZZ[115] - 16 + aryuZZ[110];
   aryuZZ[117]=aryuZZ[3]*aryuZZ[37]*aryuZZ[117];
   aryuZZ[123]=8 + aryuZZ[229];
   aryuZZ[171]= - 4./3.*aryuZZ[35] + 4./3. + aryuZZ[171];
   aryuZZ[171]=aryuZZ[35]*aryuZZ[171];
   aryuZZ[166]= - MMt*aryuZZ[166]*aryuZZ[37]*aryuZZ[38];
   aryuZZ[233]=36*aryuZZ[166];
   aryuZZ[261]= - aryuZZ[6]*aryuZZ[36];
   aryuZZ[266]=1./4.*aryuZZ[261];
   aryuZZ[117]=aryuZZ[233] + aryuZZ[117] + aryuZZ[266] + 1./3.*
   aryuZZ[123] + aryuZZ[171];
   aryuZZ[117]=MMt*aryuZZ[117];
   aryuZZ[123]=aryuZZ[156] + 11 + 13./6.*aryuZZ[35];
   aryuZZ[123]=aryuZZ[37]*aryuZZ[123];
   aryuZZ[171]=pow(aryuZZ[37],2);
   aryuZZ[267]= - aryuZZ[3]*aryuZZ[171];
   aryuZZ[269]=16*aryuZZ[267];
   aryuZZ[117]=aryuZZ[117] + 1./2.*aryuZZ[123] + aryuZZ[269];
   aryuZZ[117]=MMt*aryuZZ[117];
   aryuZZ[123]=1 + aryuZZ[134];
   aryuZZ[123]=aryuZZ[35]*aryuZZ[123];
   aryuZZ[134]=aryuZZ[36] - aryuZZ[35];
   aryuZZ[134]=aryuZZ[3]*aryuZZ[37]*aryuZZ[134];
   aryuZZ[123]=3*aryuZZ[134] + aryuZZ[266] + aryuZZ[161] + 1./2.*
   aryuZZ[123];
   aryuZZ[123]=MMt*aryuZZ[123];
   aryuZZ[134]=1./2.*aryuZZ[35];
   aryuZZ[266]=aryuZZ[156] + aryuZZ[134] + 1./3. - aryuZZ[36];
   aryuZZ[266]=aryuZZ[37]*aryuZZ[266];
   aryuZZ[123]=1./2.*aryuZZ[266] + aryuZZ[123];
   aryuZZ[123]=aryuZZ[133]*MMt*aryuZZ[123];
   aryuZZ[117]=aryuZZ[123] + 8./3.*aryuZZ[171] + aryuZZ[117];
   aryuZZ[117]=aryuZZ[133]*aryuZZ[117];
   aryuZZ[123]=5./2.*aryuZZ[36];
   aryuZZ[266]= - 7./6.*aryuZZ[35] - 113./3. + aryuZZ[123];
   aryuZZ[266]=aryuZZ[35]*aryuZZ[266];
   aryuZZ[261]=1./3.*aryuZZ[261];
   aryuZZ[266]=aryuZZ[261] + 1./2.*aryuZZ[266] - 16 + 49./18.*
   aryuZZ[36];
   aryuZZ[274]=32./3. + aryuZZ[110];
   aryuZZ[275]= - 3*aryuZZ[36];
   aryuZZ[274]=73./3.*aryuZZ[35] + 2*aryuZZ[274] + aryuZZ[275];
   aryuZZ[274]=aryuZZ[3]*aryuZZ[37]*aryuZZ[274];
   aryuZZ[166]=72*aryuZZ[166];
   aryuZZ[266]=aryuZZ[166] + 1./3.*aryuZZ[266] + aryuZZ[274];
   aryuZZ[266]=MMt*aryuZZ[266];
   aryuZZ[276]= - 17./3. + aryuZZ[229];
   aryuZZ[276]= - aryuZZ[6] + 17*aryuZZ[276] - 109./2.*aryuZZ[35];
   aryuZZ[276]=aryuZZ[37]*aryuZZ[276];
   aryuZZ[278]=aryuZZ[3]*aryuZZ[171];
   aryuZZ[279]=64*aryuZZ[278];
   aryuZZ[276]=1./3.*aryuZZ[276] + aryuZZ[279];
   aryuZZ[266]=1./3.*aryuZZ[276] + aryuZZ[266];
   aryuZZ[266]=MMt*aryuZZ[266];
   aryuZZ[117]=aryuZZ[117] - 16./3.*aryuZZ[171] + aryuZZ[266];
   aryuZZ[117]=aryuZZ[133]*aryuZZ[117];
   aryuZZ[148]=34./3. + aryuZZ[148];
   aryuZZ[266]=49./108.*aryuZZ[35] + 175./54. + 2*aryuZZ[36];
   aryuZZ[266]=aryuZZ[35]*aryuZZ[266];
   aryuZZ[276]=aryuZZ[6]*aryuZZ[36];
   aryuZZ[148]=5./12.*aryuZZ[276] + 4./9.*aryuZZ[148] + aryuZZ[266];
   aryuZZ[110]= - 7./3.*aryuZZ[35] + aryuZZ[275] - 16./3. + aryuZZ[110]
   ;
   aryuZZ[110]=aryuZZ[3]*aryuZZ[37]*aryuZZ[110];
   aryuZZ[110]=aryuZZ[233] + 1./3.*aryuZZ[148] + aryuZZ[110];
   aryuZZ[110]=MMt*aryuZZ[110];
   aryuZZ[148]=5./2.*aryuZZ[6];
   aryuZZ[233]=17*aryuZZ[36];
   aryuZZ[266]=aryuZZ[148] + 503./18.*aryuZZ[35] + 511./9. + 
   aryuZZ[233];
   aryuZZ[266]=aryuZZ[37]*aryuZZ[266];
   aryuZZ[266]=1./6.*aryuZZ[266] + aryuZZ[269];
   aryuZZ[110]=1./3.*aryuZZ[266] + aryuZZ[110];
   aryuZZ[110]=MMt*aryuZZ[110];
   aryuZZ[148]=aryuZZ[200] + aryuZZ[148];
   aryuZZ[148]=aryuZZ[37]*aryuZZ[148];
   aryuZZ[200]=aryuZZ[35]*aryuZZ[36];
   aryuZZ[266]=5./2.*aryuZZ[276] - 11./3.*aryuZZ[36] + 17./2.*
   aryuZZ[200];
   aryuZZ[266]=MMt*aryuZZ[266];
   aryuZZ[266]=aryuZZ[148] + aryuZZ[266];
   aryuZZ[266]=aryuZZ[272]*MMt*aryuZZ[266];
   aryuZZ[110]=1./18.*aryuZZ[266] + 136./81.*aryuZZ[171] + aryuZZ[110];
   aryuZZ[110]=aryuZZ[272]*aryuZZ[110];
   aryuZZ[266]= - 944./3. + 49./2.*aryuZZ[36];
   aryuZZ[269]= - 613./27. + aryuZZ[229];
   aryuZZ[269]=5*aryuZZ[269] - 2111./54.*aryuZZ[35];
   aryuZZ[269]=aryuZZ[35]*aryuZZ[269];
   aryuZZ[261]=aryuZZ[261] + 1./9.*aryuZZ[266] + 1./2.*aryuZZ[269];
   aryuZZ[166]=aryuZZ[166] + 1./3.*aryuZZ[261] + aryuZZ[274];
   aryuZZ[166]=MMt*aryuZZ[166];
   aryuZZ[261]=17./2.*aryuZZ[36];
   aryuZZ[266]= - aryuZZ[6] - 3029./18.*aryuZZ[35] - 1891./9. + 
   aryuZZ[261];
   aryuZZ[266]=aryuZZ[37]*aryuZZ[266];
   aryuZZ[266]=1./3.*aryuZZ[266] + aryuZZ[279];
   aryuZZ[166]=1./3.*aryuZZ[266] + aryuZZ[166];
   aryuZZ[166]=MMt*aryuZZ[166];
   aryuZZ[110]=aryuZZ[110] - 944./81.*aryuZZ[171] + aryuZZ[166];
   aryuZZ[110]=aryuZZ[272]*aryuZZ[110];
   aryuZZ[166]=aryuZZ[37]*aryuZZ[248];
   aryuZZ[266]=2 + aryuZZ[35];
   aryuZZ[269]=aryuZZ[35]*aryuZZ[266];
   aryuZZ[269]=1 + aryuZZ[269];
   aryuZZ[269]=MMt*aryuZZ[269];
   aryuZZ[166]=2*aryuZZ[166] + aryuZZ[269];
   aryuZZ[166]=MMt*aryuZZ[166];
   aryuZZ[166]=aryuZZ[171] + aryuZZ[166];
   aryuZZ[110]=aryuZZ[117] + 1024./81.*aryuZZ[166] + aryuZZ[110];
   aryuZZ[110]=aryuZZ[33]*aryuZZ[110];
   aryuZZ[117]=aryuZZ[21] - aryuZZ[41] + 2*aryuZZ[50];
   aryuZZ[166]= - aryuZZ[35]*aryuZZ[19];
   aryuZZ[269]= - 2 - aryuZZ[35];
   aryuZZ[269]=aryuZZ[18]*aryuZZ[269];
   aryuZZ[117]= - 4./9.*MMt + 4./9.*aryuZZ[37] + 2./9.*aryuZZ[269] + 
   aryuZZ[166] + 2./9.*aryuZZ[117] - aryuZZ[19];
   aryuZZ[117]=MMt*aryuZZ[117];
   aryuZZ[269]= - 2./9.*aryuZZ[37] - aryuZZ[19] - 4./9.*aryuZZ[18];
   aryuZZ[269]=aryuZZ[37]*aryuZZ[269];
   aryuZZ[117]=aryuZZ[269] + aryuZZ[117];
   aryuZZ[104]=aryuZZ[110] + aryuZZ[104] + 256./9.*aryuZZ[117] + 
   aryuZZ[143];
   aryuZZ[104]=aryuZZ[33]*aryuZZ[104];
   aryuZZ[110]= - 1./2. + aryuZZ[158];
   aryuZZ[117]=aryuZZ[110] + aryuZZ[161];
   aryuZZ[117]=aryuZZ[10]*aryuZZ[117];
   aryuZZ[143]=7./9. + aryuZZ[72];
   aryuZZ[269]= - 1./6.*aryuZZ[36];
   aryuZZ[274]= - MMH*aryuZZ[53];
   aryuZZ[117]=aryuZZ[274] + aryuZZ[117] + aryuZZ[269] + aryuZZ[105] + 
   91./36.*aryuZZ[67] - 7./9.*aryuZZ[64] - aryuZZ[61] + 5./2.*
   aryuZZ[143] + aryuZZ[62];
   aryuZZ[117]=MMH*aryuZZ[117];
   aryuZZ[143]=aryuZZ[229] + 47./36. + aryuZZ[158];
   aryuZZ[143]=aryuZZ[17]*aryuZZ[143];
   aryuZZ[276]=11./4. + aryuZZ[283];
   aryuZZ[279]=1./6.*aryuZZ[36];
   aryuZZ[280]=aryuZZ[276] + aryuZZ[279];
   aryuZZ[280]=aryuZZ[18]*aryuZZ[280];
   aryuZZ[281]=469./2.*aryuZZ[70] - 115*aryuZZ[42];
   aryuZZ[281]=59./9.*aryuZZ[40] - 65./18.*aryuZZ[68] + 1./9.*
   aryuZZ[281] - 1./2.*aryuZZ[69];
   aryuZZ[282]= - 775./6. + aryuZZ[317];
   aryuZZ[282]=aryuZZ[19]*aryuZZ[282];
   aryuZZ[284]=109 + 239./3.*aryuZZ[66];
   aryuZZ[284]=1./4.*aryuZZ[284] - 23./3.*aryuZZ[67];
   aryuZZ[284]=MMt*aryuZZ[284];
   aryuZZ[285]= - 5*aryuZZ[50];
   aryuZZ[286]= - 9./4.*aryuZZ[20];
   aryuZZ[117]=1./3.*aryuZZ[284] + aryuZZ[117] - 79./6.*aryuZZ[37] + 
   aryuZZ[280] + aryuZZ[143] + 1./6.*aryuZZ[282] + aryuZZ[142] + 
   aryuZZ[286] + 33./4.*aryuZZ[22] + 1./2.*aryuZZ[281] + aryuZZ[285];
   aryuZZ[117]=MMt*aryuZZ[117];
   aryuZZ[143]=aryuZZ[188] + 229./54. + aryuZZ[283];
   aryuZZ[143]=1./2.*aryuZZ[143] + 7./27.*aryuZZ[35];
   aryuZZ[143]=aryuZZ[10]*aryuZZ[143];
   aryuZZ[143]=aryuZZ[143] + 7./54.*aryuZZ[344] - 17./27.*aryuZZ[11] + 
   29./36.*aryuZZ[15] + 25./36.*aryuZZ[67] + aryuZZ[309] + 7./36.*
   aryuZZ[59] + 107./36. - aryuZZ[72];
   aryuZZ[143]=MMH*aryuZZ[143];
   aryuZZ[161]=aryuZZ[161] - 229./54. + aryuZZ[158];
   aryuZZ[161]=1./2.*aryuZZ[161] - 7./27.*aryuZZ[35];
   aryuZZ[161]=aryuZZ[18]*aryuZZ[161];
   aryuZZ[280]= - 7./36.*aryuZZ[41] - 43./18.*aryuZZ[40] + aryuZZ[69]
    + 41./9.*aryuZZ[68];
   aryuZZ[281]=aryuZZ[279] - 139./27. + 3./2.*aryuZZ[38];
   aryuZZ[281]=aryuZZ[17]*aryuZZ[281];
   aryuZZ[282]= - 11./6.*aryuZZ[10] - 17./3. + aryuZZ[11];
   aryuZZ[282]=aryuZZ[37]*aryuZZ[282];
   aryuZZ[143]=1./4.*aryuZZ[143] + 1./2.*aryuZZ[282] + 1./4.*
   aryuZZ[161] + 7./216.*aryuZZ[348] + 1./4.*aryuZZ[281] + 17./108.*
   aryuZZ[19] - 29./144.*aryuZZ[21] + 3./8.*aryuZZ[20] + 1./4.*
   aryuZZ[280] - 2./9.*aryuZZ[50];
   aryuZZ[143]=MMH*aryuZZ[143];
   aryuZZ[161]=47./3.*aryuZZ[18] - 85./3.*aryuZZ[19] - aryuZZ[17];
   aryuZZ[280]=5*aryuZZ[37];
   aryuZZ[161]=1./2.*aryuZZ[161] + aryuZZ[280];
   aryuZZ[161]=aryuZZ[37]*aryuZZ[161];
   aryuZZ[117]=1./2.*aryuZZ[117] + 1./4.*aryuZZ[161] + aryuZZ[143];
   aryuZZ[117]=MMt*aryuZZ[117];
   aryuZZ[143]=aryuZZ[181] + aryuZZ[236];
   aryuZZ[143]=MMH*aryuZZ[143];
   aryuZZ[161]= - 3 + aryuZZ[188];
   aryuZZ[161]=aryuZZ[19]*aryuZZ[161];
   aryuZZ[281]=1 + aryuZZ[66];
   aryuZZ[281]=MMt*aryuZZ[281];
   aryuZZ[143]=1./2.*aryuZZ[281] + 1./3.*aryuZZ[143] - aryuZZ[37] + 1./
   6.*aryuZZ[352] + 1./2.*aryuZZ[264] + 1./2.*aryuZZ[161] + aryuZZ[192]
    + 1./2.*aryuZZ[22];
   aryuZZ[143]=MMt*aryuZZ[143];
   aryuZZ[161]=17*aryuZZ[356] + 7./2.*aryuZZ[358];
   aryuZZ[281]= - 1./3.*aryuZZ[10];
   aryuZZ[282]=aryuZZ[281] - 1./6. + aryuZZ[11];
   aryuZZ[282]=aryuZZ[37]*aryuZZ[282];
   aryuZZ[284]= - 17*aryuZZ[11] + 7./2.*aryuZZ[344];
   aryuZZ[284]=MMH*aryuZZ[284];
   aryuZZ[161]=1./54.*aryuZZ[284] + 1./54.*aryuZZ[161] + aryuZZ[282];
   aryuZZ[161]=MMH*aryuZZ[161];
   aryuZZ[282]=1./12.*aryuZZ[18] - 2./3.*aryuZZ[19] + 3./4.*aryuZZ[17];
   aryuZZ[282]=aryuZZ[37]*aryuZZ[282];
   aryuZZ[143]=1./2.*aryuZZ[143] + aryuZZ[282] + 1./2.*aryuZZ[161];
   aryuZZ[143]=MMt*aryuZZ[143];
   aryuZZ[161]=aryuZZ[37]*aryuZZ[356];
   aryuZZ[282]= - MMH*aryuZZ[37]*aryuZZ[11];
   aryuZZ[161]=aryuZZ[161] + aryuZZ[282];
   aryuZZ[161]=MMH*aryuZZ[161];
   aryuZZ[143]=17./108.*aryuZZ[161] + aryuZZ[143];
   aryuZZ[143]=aryuZZ[272]*aryuZZ[143];
   aryuZZ[161]=2*aryuZZ[17];
   aryuZZ[282]= - 175./16.*aryuZZ[18] + 17./4.*aryuZZ[19] + aryuZZ[161]
   ;
   aryuZZ[282]=1./3.*aryuZZ[282] + 25./16.*aryuZZ[37];
   aryuZZ[282]=aryuZZ[37]*aryuZZ[282];
   aryuZZ[284]=175./12.*aryuZZ[10];
   aryuZZ[287]=aryuZZ[284] - 25./4. - 17./3.*aryuZZ[11];
   aryuZZ[287]=aryuZZ[37]*aryuZZ[287];
   aryuZZ[287]=25./4.*aryuZZ[50] + aryuZZ[287];
   aryuZZ[287]=MMH*aryuZZ[287];
   aryuZZ[282]=aryuZZ[282] + 1./4.*aryuZZ[287];
   aryuZZ[282]=MMH*aryuZZ[282];
   aryuZZ[117]=aryuZZ[143] + 1./9.*aryuZZ[282] + aryuZZ[117];
   aryuZZ[117]=aryuZZ[272]*aryuZZ[117];
   aryuZZ[143]=1./3. + aryuZZ[72];
   aryuZZ[282]=aryuZZ[10]*aryuZZ[110];
   aryuZZ[143]=aryuZZ[274] + aryuZZ[282] + aryuZZ[105] - 23./12.*
   aryuZZ[67] + aryuZZ[309] - aryuZZ[61] + 5./2.*aryuZZ[143] + 
   aryuZZ[62];
   aryuZZ[143]=MMH*aryuZZ[143];
   aryuZZ[282]=49*aryuZZ[70] - 23*aryuZZ[42];
   aryuZZ[282]=85./3.*aryuZZ[68] + 1./3.*aryuZZ[282] - aryuZZ[69];
   aryuZZ[282]=1./2.*aryuZZ[282] - 7./3.*aryuZZ[40];
   aryuZZ[287]= - 40 + 19*aryuZZ[36];
   aryuZZ[287]=aryuZZ[19]*aryuZZ[287];
   aryuZZ[289]= - 91./12. + aryuZZ[158];
   aryuZZ[289]=aryuZZ[17]*aryuZZ[289];
   aryuZZ[290]=aryuZZ[18]*aryuZZ[276];
   aryuZZ[291]=27 + 13./3.*aryuZZ[66];
   aryuZZ[291]=1./2.*aryuZZ[291] + 19./3.*aryuZZ[67];
   aryuZZ[291]=MMt*aryuZZ[291];
   aryuZZ[143]=aryuZZ[291] + aryuZZ[143] - 22*aryuZZ[37] + aryuZZ[290]
    + aryuZZ[289] + 1./3.*aryuZZ[287] + aryuZZ[142] + aryuZZ[286] + 9*
   aryuZZ[22] + 1./2.*aryuZZ[282] + aryuZZ[285];
   aryuZZ[143]=MMt*aryuZZ[143];
   aryuZZ[282]= - 11./9.*aryuZZ[35];
   aryuZZ[287]=aryuZZ[282] + aryuZZ[279] - 7./18. + aryuZZ[283];
   aryuZZ[287]=aryuZZ[10]*aryuZZ[287];
   aryuZZ[289]=23./12.*aryuZZ[15] - 5./12.*aryuZZ[67] + aryuZZ[155] - 
   11./12.*aryuZZ[59] - 31./12. - aryuZZ[72];
   aryuZZ[287]=1./4.*aryuZZ[287] + 2./27.*aryuZZ[344] + 1./2.*
   aryuZZ[289] - 2./27.*aryuZZ[11];
   aryuZZ[287]=MMH*aryuZZ[287];
   aryuZZ[289]=aryuZZ[279] + 283./27. + aryuZZ[283];
   aryuZZ[289]=aryuZZ[17]*aryuZZ[289];
   aryuZZ[290]=11./9.*aryuZZ[35] + aryuZZ[269] + 7./18. + aryuZZ[158];
   aryuZZ[290]=aryuZZ[18]*aryuZZ[290];
   aryuZZ[291]=11./12.*aryuZZ[41] - 1./6.*aryuZZ[40] + aryuZZ[69] - 13./
   3.*aryuZZ[68];
   aryuZZ[293]=aryuZZ[132] - 41./4.*aryuZZ[17];
   aryuZZ[293]=aryuZZ[35]*aryuZZ[293];
   aryuZZ[294]=9./2. - aryuZZ[10];
   aryuZZ[294]=aryuZZ[37]*aryuZZ[294];
   aryuZZ[287]=aryuZZ[287] + aryuZZ[294] + 1./4.*aryuZZ[290] + 1./27.*
   aryuZZ[293] + 1./4.*aryuZZ[289] + 2./27.*aryuZZ[19] - 23./24.*
   aryuZZ[21] + 3./4.*aryuZZ[20] + 1./2.*aryuZZ[291] + 2./3.*aryuZZ[50]
   ;
   aryuZZ[287]=MMH*aryuZZ[287];
   aryuZZ[186]=5./2.*aryuZZ[37] + aryuZZ[186] - 8./3.*aryuZZ[19] - 5./4.
   *aryuZZ[17];
   aryuZZ[186]=aryuZZ[37]*aryuZZ[186];
   aryuZZ[143]=aryuZZ[143] + aryuZZ[186] + aryuZZ[287];
   aryuZZ[143]=MMt*aryuZZ[143];
   aryuZZ[186]=aryuZZ[132] - 23./4.*aryuZZ[17];
   aryuZZ[186]=1./3.*aryuZZ[186] + 25./8.*aryuZZ[18];
   aryuZZ[186]=1./3.*aryuZZ[186] - 5./8.*aryuZZ[37];
   aryuZZ[186]=aryuZZ[37]*aryuZZ[186];
   aryuZZ[287]= - 25./24.*aryuZZ[10] + 5./8. - 2./9.*aryuZZ[11];
   aryuZZ[287]=aryuZZ[37]*aryuZZ[287];
   aryuZZ[287]= - 5./8.*aryuZZ[50] + aryuZZ[287];
   aryuZZ[287]=MMH*aryuZZ[287];
   aryuZZ[186]=aryuZZ[186] + aryuZZ[287];
   aryuZZ[186]=MMH*aryuZZ[186];
   aryuZZ[143]=1./3.*aryuZZ[186] + aryuZZ[143];
   aryuZZ[117]=aryuZZ[143] + aryuZZ[117];
   aryuZZ[117]=aryuZZ[272]*aryuZZ[117];
   aryuZZ[110]=aryuZZ[110] + aryuZZ[188];
   aryuZZ[110]=aryuZZ[10]*aryuZZ[110];
   aryuZZ[186]=3 + 5*aryuZZ[72];
   aryuZZ[105]=aryuZZ[274] + aryuZZ[110] + aryuZZ[279] + aryuZZ[105] + 
   3./4.*aryuZZ[67] + aryuZZ[64] - aryuZZ[61] + 1./2.*aryuZZ[186] + 
   aryuZZ[62];
   aryuZZ[105]=MMH*aryuZZ[105];
   aryuZZ[110]= - 3./4. - aryuZZ[38];
   aryuZZ[110]=3*aryuZZ[110] + aryuZZ[181];
   aryuZZ[110]=aryuZZ[17]*aryuZZ[110];
   aryuZZ[186]=aryuZZ[276] + aryuZZ[269];
   aryuZZ[186]=aryuZZ[18]*aryuZZ[186];
   aryuZZ[269]=31./2. + 21*aryuZZ[36];
   aryuZZ[269]=aryuZZ[19]*aryuZZ[269];
   aryuZZ[276]= - 11 - 35*aryuZZ[66];
   aryuZZ[276]=1./4.*aryuZZ[276] + aryuZZ[67];
   aryuZZ[276]=MMt*aryuZZ[276];
   aryuZZ[105]=aryuZZ[276] + aryuZZ[105] + 21./2.*aryuZZ[37] + 
   aryuZZ[186] + aryuZZ[110] + 1./2.*aryuZZ[269] + aryuZZ[142] + 
   aryuZZ[286] + 39./4.*aryuZZ[22] + aryuZZ[285] + aryuZZ[328] + 7./4.*
   aryuZZ[68] - 1./4.*aryuZZ[69] - 71./4.*aryuZZ[70] + 9*aryuZZ[42];
   aryuZZ[105]=MMt*aryuZZ[105];
   aryuZZ[110]= - 3./2.*aryuZZ[40];
   aryuZZ[186]= - 1./3.*aryuZZ[19];
   aryuZZ[152]=aryuZZ[186] + aryuZZ[170] + 3./2.*aryuZZ[20] + 
   aryuZZ[152] + aryuZZ[110] + aryuZZ[69] + aryuZZ[68];
   aryuZZ[170]=1./9. + 3./4.*aryuZZ[38];
   aryuZZ[170]=aryuZZ[17]*aryuZZ[170];
   aryuZZ[269]=aryuZZ[19] + 13./3.*aryuZZ[17];
   aryuZZ[269]=aryuZZ[35]*aryuZZ[269];
   aryuZZ[152]=1./12.*aryuZZ[269] + 1./2.*aryuZZ[152] + aryuZZ[170];
   aryuZZ[170]= - 1./12.*aryuZZ[10];
   aryuZZ[269]=aryuZZ[170] - 2./3. + aryuZZ[118];
   aryuZZ[269]=aryuZZ[37]*aryuZZ[269];
   aryuZZ[276]=1./6.*aryuZZ[344] + aryuZZ[216] + 5./4.*aryuZZ[15] + 1./
   4.*aryuZZ[67] + aryuZZ[64] - 1./4.*aryuZZ[59] + 3./4. - aryuZZ[72];
   aryuZZ[285]=59./18. + aryuZZ[283];
   aryuZZ[285]=1./8.*aryuZZ[285] + 2./9.*aryuZZ[35];
   aryuZZ[285]=aryuZZ[10]*aryuZZ[285];
   aryuZZ[276]=1./4.*aryuZZ[276] + aryuZZ[285];
   aryuZZ[276]=MMH*aryuZZ[276];
   aryuZZ[285]= - 59./18. + aryuZZ[158];
   aryuZZ[285]=1./8.*aryuZZ[285] - 2./9.*aryuZZ[35];
   aryuZZ[285]=aryuZZ[18]*aryuZZ[285];
   aryuZZ[152]=aryuZZ[276] + aryuZZ[269] + 1./2.*aryuZZ[152] + 
   aryuZZ[285];
   aryuZZ[152]=MMH*aryuZZ[152];
   aryuZZ[269]=31./3.*aryuZZ[18] + aryuZZ[124] - 9*aryuZZ[17];
   aryuZZ[269]=1./2.*aryuZZ[269] + aryuZZ[280];
   aryuZZ[269]=aryuZZ[37]*aryuZZ[269];
   aryuZZ[105]=1./2.*aryuZZ[105] + 1./4.*aryuZZ[269] + aryuZZ[152];
   aryuZZ[105]=MMt*aryuZZ[105];
   aryuZZ[152]=aryuZZ[18]*aryuZZ[254];
   aryuZZ[152]=aryuZZ[152] + 1./2.*aryuZZ[348] - aryuZZ[19] + 1./2.*
   aryuZZ[251];
   aryuZZ[254]=1./6. - aryuZZ[11];
   aryuZZ[269]=2./3.*aryuZZ[10];
   aryuZZ[254]=1./2.*aryuZZ[254] + aryuZZ[269];
   aryuZZ[254]=aryuZZ[37]*aryuZZ[254];
   aryuZZ[276]=aryuZZ[35] - 1 + aryuZZ[181];
   aryuZZ[280]=aryuZZ[10]*aryuZZ[276];
   aryuZZ[285]=aryuZZ[280] + aryuZZ[11] + 1./2.*aryuZZ[344];
   aryuZZ[285]=MMH*aryuZZ[285];
   aryuZZ[152]=1./12.*aryuZZ[285] + 1./12.*aryuZZ[152] + aryuZZ[254];
   aryuZZ[152]=MMH*aryuZZ[152];
   aryuZZ[137]=1./3.*aryuZZ[209] + 5*aryuZZ[137] + aryuZZ[251];
   aryuZZ[254]=aryuZZ[229] + aryuZZ[354];
   aryuZZ[254]=MMH*aryuZZ[254];
   aryuZZ[137]=1./2.*aryuZZ[137] + 1./3.*aryuZZ[254];
   aryuZZ[137]=MMt*aryuZZ[137];
   aryuZZ[254]=aryuZZ[124] - aryuZZ[17];
   aryuZZ[285]=aryuZZ[254] - 7./3.*aryuZZ[18];
   aryuZZ[285]=aryuZZ[37]*aryuZZ[285];
   aryuZZ[137]=1./2.*aryuZZ[137] + 1./4.*aryuZZ[285] + aryuZZ[152];
   aryuZZ[137]=MMt*aryuZZ[137];
   aryuZZ[152]= - aryuZZ[19] + aryuZZ[18];
   aryuZZ[152]=aryuZZ[37]*aryuZZ[152];
   aryuZZ[285]=aryuZZ[11] - aryuZZ[10];
   aryuZZ[286]=MMH*aryuZZ[37]*aryuZZ[285];
   aryuZZ[152]=aryuZZ[152] + aryuZZ[286];
   aryuZZ[152]=MMH*aryuZZ[152];
   aryuZZ[137]=1./12.*aryuZZ[152] + aryuZZ[137];
   aryuZZ[137]=aryuZZ[133]*aryuZZ[137];
   aryuZZ[152]=aryuZZ[206] - 41./12.*aryuZZ[18];
   aryuZZ[152]=1./3.*aryuZZ[152] + 1./4.*aryuZZ[37];
   aryuZZ[152]=aryuZZ[37]*aryuZZ[152];
   aryuZZ[206]=41./36.*aryuZZ[10] - 1./4. + aryuZZ[216];
   aryuZZ[206]=aryuZZ[37]*aryuZZ[206];
   aryuZZ[206]=1./4.*aryuZZ[50] + aryuZZ[206];
   aryuZZ[206]=MMH*aryuZZ[206];
   aryuZZ[152]=aryuZZ[152] + aryuZZ[206];
   aryuZZ[152]=MMH*aryuZZ[152];
   aryuZZ[105]=aryuZZ[137] + 1./4.*aryuZZ[152] + aryuZZ[105];
   aryuZZ[105]=aryuZZ[133]*aryuZZ[105];
   aryuZZ[105]=aryuZZ[143] + aryuZZ[105];
   aryuZZ[105]=aryuZZ[133]*aryuZZ[105];
   aryuZZ[137]=1 + aryuZZ[36];
   aryuZZ[122]=17*aryuZZ[137] + aryuZZ[122];
   aryuZZ[122]=aryuZZ[37]*aryuZZ[122];
   aryuZZ[143]=aryuZZ[233] + 7./2.*aryuZZ[200];
   aryuZZ[152]= - aryuZZ[3]*aryuZZ[37]*aryuZZ[36];
   aryuZZ[143]=1./18.*aryuZZ[143] + 3*aryuZZ[152];
   aryuZZ[143]=MMt*aryuZZ[143];
   aryuZZ[122]=aryuZZ[143] + 1./18.*aryuZZ[122] + 3*aryuZZ[267];
   aryuZZ[122]=MMt*aryuZZ[122];
   aryuZZ[122]=17./18.*aryuZZ[171] + aryuZZ[122];
   aryuZZ[122]=MMt*aryuZZ[122];
   aryuZZ[143]=aryuZZ[272]*aryuZZ[122];
   aryuZZ[122]=aryuZZ[122] + aryuZZ[143];
   aryuZZ[122]=aryuZZ[272]*aryuZZ[122];
   aryuZZ[137]=aryuZZ[137] + aryuZZ[35];
   aryuZZ[137]=aryuZZ[37]*aryuZZ[137];
   aryuZZ[143]=aryuZZ[36] + aryuZZ[200];
   aryuZZ[143]=MMt*aryuZZ[143];
   aryuZZ[137]=aryuZZ[137] + aryuZZ[143];
   aryuZZ[137]=MMt*aryuZZ[137];
   aryuZZ[137]=aryuZZ[171] + aryuZZ[137];
   aryuZZ[137]=4./9.*MMt*aryuZZ[137];
   aryuZZ[122]=aryuZZ[137] + aryuZZ[122];
   aryuZZ[122]=aryuZZ[272]*aryuZZ[122];
   aryuZZ[143]=aryuZZ[316] + aryuZZ[134];
   aryuZZ[143]=aryuZZ[37]*aryuZZ[143];
   aryuZZ[152]= - aryuZZ[36] + 1./2.*aryuZZ[200];
   aryuZZ[200]=aryuZZ[3]*aryuZZ[37]*aryuZZ[36];
   aryuZZ[152]=1./2.*aryuZZ[152] + 3*aryuZZ[200];
   aryuZZ[152]=MMt*aryuZZ[152];
   aryuZZ[143]=aryuZZ[152] + 1./2.*aryuZZ[143] + 3*aryuZZ[278];
   aryuZZ[143]=MMt*aryuZZ[143];
   aryuZZ[143]= - 1./2.*aryuZZ[171] + aryuZZ[143];
   aryuZZ[143]=MMt*aryuZZ[143];
   aryuZZ[152]=aryuZZ[133]*aryuZZ[143];
   aryuZZ[143]=aryuZZ[143] + aryuZZ[152];
   aryuZZ[143]=aryuZZ[133]*aryuZZ[143];
   aryuZZ[137]=aryuZZ[137] + aryuZZ[143];
   aryuZZ[137]=aryuZZ[133]*aryuZZ[137];
   aryuZZ[122]=aryuZZ[122] + aryuZZ[137];
   aryuZZ[122]=aryuZZ[33]*aryuZZ[122];
   aryuZZ[105]=aryuZZ[122] + aryuZZ[117] + aryuZZ[105];
   aryuZZ[105]=aryuZZ[33]*aryuZZ[105];
   aryuZZ[117]= - 5./9.*aryuZZ[10];
   aryuZZ[122]=aryuZZ[117] - 1 + aryuZZ[331];
   aryuZZ[122]=aryuZZ[18]*aryuZZ[122];
   aryuZZ[137]=1 - 1./18.*aryuZZ[11];
   aryuZZ[137]=aryuZZ[17]*aryuZZ[137];
   aryuZZ[143]=1./4.*aryuZZ[19] + aryuZZ[17];
   aryuZZ[143]=aryuZZ[10]*aryuZZ[143];
   aryuZZ[152]= - 1./2.*aryuZZ[75];
   aryuZZ[171]=5./2.*aryuZZ[20] - aryuZZ[91] + aryuZZ[152];
   aryuZZ[171]=1./2.*aryuZZ[171] + aryuZZ[21];
   aryuZZ[171]=1./3.*aryuZZ[171];
   aryuZZ[122]=1./4.*aryuZZ[122] + 1./9.*aryuZZ[143] + aryuZZ[171] + 1./
   2.*aryuZZ[137];
   aryuZZ[122]=aryuZZ[73]*aryuZZ[122];
   aryuZZ[137]=1./6.*aryuZZ[84];
   aryuZZ[143]= - 1./6.*aryuZZ[15];
   aryuZZ[200]=aryuZZ[143] - 1 + aryuZZ[137];
   aryuZZ[150]=aryuZZ[10]*aryuZZ[150];
   aryuZZ[150]=aryuZZ[200] + 1./18.*aryuZZ[150];
   aryuZZ[150]=MMH*aryuZZ[73]*aryuZZ[150];
   aryuZZ[122]=aryuZZ[122] + 1./2.*aryuZZ[150];
   aryuZZ[122]=MMH*aryuZZ[122];
   aryuZZ[150]=1./6.*aryuZZ[19];
   aryuZZ[206]=aryuZZ[150] - aryuZZ[17];
   aryuZZ[206]=aryuZZ[17]*aryuZZ[206];
   aryuZZ[233]= - 29./64.*aryuZZ[18] - 1./16.*aryuZZ[19] - aryuZZ[17];
   aryuZZ[233]=aryuZZ[18]*aryuZZ[233];
   aryuZZ[206]=1./8.*aryuZZ[206] + 1./3.*aryuZZ[233];
   aryuZZ[206]=aryuZZ[73]*aryuZZ[206];
   aryuZZ[122]=1./3.*aryuZZ[206] + 1./4.*aryuZZ[122];
   aryuZZ[206]=pow(MMH,2);
   aryuZZ[122]=aryuZZ[206]*aryuZZ[122];
   aryuZZ[233]= - aryuZZ[17]*aryuZZ[11];
   aryuZZ[267]=aryuZZ[10]*aryuZZ[356];
   aryuZZ[278]=aryuZZ[18]*aryuZZ[11];
   aryuZZ[233]=aryuZZ[278] + aryuZZ[233] + aryuZZ[267];
   aryuZZ[233]=aryuZZ[73]*aryuZZ[233];
   aryuZZ[267]= - MMH*aryuZZ[73]*aryuZZ[10]*aryuZZ[11];
   aryuZZ[233]=aryuZZ[233] + aryuZZ[267];
   aryuZZ[233]=MMH*aryuZZ[233];
   aryuZZ[267]=aryuZZ[17]*aryuZZ[356];
   aryuZZ[278]=aryuZZ[18]*aryuZZ[359];
   aryuZZ[267]=aryuZZ[267] + aryuZZ[278];
   aryuZZ[267]=aryuZZ[73]*aryuZZ[267];
   aryuZZ[233]=aryuZZ[267] + aryuZZ[233];
   aryuZZ[233]=aryuZZ[272]*aryuZZ[206]*aryuZZ[233];
   aryuZZ[122]=aryuZZ[122] + 1./144.*aryuZZ[233];
   aryuZZ[122]=aryuZZ[272]*aryuZZ[122];
   aryuZZ[233]=aryuZZ[10]*aryuZZ[17];
   aryuZZ[267]= - 1./4. - 1./9.*aryuZZ[10];
   aryuZZ[267]=aryuZZ[18]*aryuZZ[267];
   aryuZZ[233]=aryuZZ[267] + 1./9.*aryuZZ[233] + aryuZZ[171] + 
   aryuZZ[119];
   aryuZZ[233]=aryuZZ[73]*aryuZZ[233];
   aryuZZ[267]=pow(aryuZZ[10],2);
   aryuZZ[267]=aryuZZ[200] + 1./36.*aryuZZ[267];
   aryuZZ[267]=MMH*aryuZZ[73]*aryuZZ[267];
   aryuZZ[233]=aryuZZ[233] + 1./2.*aryuZZ[267];
   aryuZZ[233]=MMH*aryuZZ[233];
   aryuZZ[267]= - 2*aryuZZ[17];
   aryuZZ[278]=aryuZZ[267] - 31./32.*aryuZZ[18];
   aryuZZ[278]=aryuZZ[18]*aryuZZ[278];
   aryuZZ[286]=pow(aryuZZ[17],2);
   aryuZZ[278]= - 11./16.*aryuZZ[286] + aryuZZ[278];
   aryuZZ[278]=aryuZZ[73]*aryuZZ[278];
   aryuZZ[233]=1./9.*aryuZZ[278] + 1./2.*aryuZZ[233];
   aryuZZ[233]=aryuZZ[206]*aryuZZ[233];
   aryuZZ[122]=aryuZZ[233] + aryuZZ[122];
   aryuZZ[122]=aryuZZ[272]*aryuZZ[122];
   aryuZZ[278]=aryuZZ[281] - 1 + aryuZZ[138];
   aryuZZ[278]=aryuZZ[18]*aryuZZ[278];
   aryuZZ[281]=1 + 1./18.*aryuZZ[11];
   aryuZZ[281]=aryuZZ[17]*aryuZZ[281];
   aryuZZ[287]= - 1./4.*aryuZZ[19];
   aryuZZ[289]=aryuZZ[287] + aryuZZ[17];
   aryuZZ[289]=aryuZZ[10]*aryuZZ[289];
   aryuZZ[171]=1./4.*aryuZZ[278] + 1./9.*aryuZZ[289] + aryuZZ[171] + 1./
   2.*aryuZZ[281];
   aryuZZ[171]=aryuZZ[73]*aryuZZ[171];
   aryuZZ[278]=aryuZZ[10]*aryuZZ[11];
   aryuZZ[200]=aryuZZ[200] + 1./18.*aryuZZ[278];
   aryuZZ[200]=MMH*aryuZZ[73]*aryuZZ[200];
   aryuZZ[171]=aryuZZ[171] + 1./2.*aryuZZ[200];
   aryuZZ[171]=MMH*aryuZZ[171];
   aryuZZ[200]= - aryuZZ[19] - 5*aryuZZ[17];
   aryuZZ[200]=aryuZZ[17]*aryuZZ[200];
   aryuZZ[281]=1./16.*aryuZZ[19] - aryuZZ[17];
   aryuZZ[281]=1./3.*aryuZZ[281] - 11./64.*aryuZZ[18];
   aryuZZ[281]=aryuZZ[18]*aryuZZ[281];
   aryuZZ[200]=1./48.*aryuZZ[200] + aryuZZ[281];
   aryuZZ[200]=aryuZZ[73]*aryuZZ[200];
   aryuZZ[171]=1./3.*aryuZZ[200] + 1./4.*aryuZZ[171];
   aryuZZ[171]=aryuZZ[206]*aryuZZ[171];
   aryuZZ[200]= - aryuZZ[17]*aryuZZ[19];
   aryuZZ[281]=aryuZZ[346] - aryuZZ[18];
   aryuZZ[281]=aryuZZ[18]*aryuZZ[281];
   aryuZZ[200]=aryuZZ[200] + aryuZZ[281];
   aryuZZ[200]=aryuZZ[73]*aryuZZ[200];
   aryuZZ[281]=aryuZZ[17]*aryuZZ[11];
   aryuZZ[289]=aryuZZ[10]*aryuZZ[178];
   aryuZZ[289]=aryuZZ[281] + aryuZZ[289];
   aryuZZ[290]=aryuZZ[18]*aryuZZ[318];
   aryuZZ[289]=1./2.*aryuZZ[289] + aryuZZ[290];
   aryuZZ[289]=aryuZZ[73]*aryuZZ[289];
   aryuZZ[285]=MMH*aryuZZ[73]*aryuZZ[10]*aryuZZ[285];
   aryuZZ[285]=aryuZZ[289] + 1./2.*aryuZZ[285];
   aryuZZ[285]=MMH*aryuZZ[285];
   aryuZZ[200]=1./2.*aryuZZ[200] + aryuZZ[285];
   aryuZZ[200]=aryuZZ[133]*aryuZZ[206]*aryuZZ[200];
   aryuZZ[171]=aryuZZ[171] + 1./72.*aryuZZ[200];
   aryuZZ[171]=aryuZZ[133]*aryuZZ[171];
   aryuZZ[171]=aryuZZ[233] + aryuZZ[171];
   aryuZZ[171]=aryuZZ[133]*aryuZZ[171];
   aryuZZ[200]=MMH*aryuZZ[354];
   aryuZZ[200]=aryuZZ[200] + aryuZZ[264] + aryuZZ[209];
   aryuZZ[200]=MMt*MMH*aryuZZ[200];
   aryuZZ[206]=aryuZZ[17] - aryuZZ[18];
   aryuZZ[206]=aryuZZ[37]*aryuZZ[206];
   aryuZZ[209]=MMH*aryuZZ[37]*aryuZZ[10];
   aryuZZ[206]=aryuZZ[206] + aryuZZ[209];
   aryuZZ[206]=MMH*aryuZZ[206];
   aryuZZ[200]=aryuZZ[206] + aryuZZ[200];
   aryuZZ[200]=MMt*aryuZZ[200];
   aryuZZ[206]=aryuZZ[272]*aryuZZ[200];
   aryuZZ[200]=aryuZZ[200] + aryuZZ[206];
   aryuZZ[200]=aryuZZ[200]*pow(aryuZZ[2],4);
   aryuZZ[206]=MMH*aryuZZ[236];
   aryuZZ[206]=aryuZZ[206] + aryuZZ[251] + aryuZZ[352];
   aryuZZ[206]=MMt*MMH*aryuZZ[206];
   aryuZZ[209]= - aryuZZ[17] + aryuZZ[18];
   aryuZZ[209]=aryuZZ[37]*aryuZZ[209];
   aryuZZ[233]= - MMH*aryuZZ[37]*aryuZZ[10];
   aryuZZ[209]=aryuZZ[209] + aryuZZ[233];
   aryuZZ[209]=MMH*aryuZZ[209];
   aryuZZ[206]=aryuZZ[209] + aryuZZ[206];
   aryuZZ[206]=MMt*aryuZZ[206];
   aryuZZ[209]=aryuZZ[133]*aryuZZ[206];
   aryuZZ[206]=aryuZZ[206] + aryuZZ[209];
   aryuZZ[206]=aryuZZ[206]*pow(aryuZZ[5],4);
   aryuZZ[200]=aryuZZ[200] + aryuZZ[206];
   aryuZZ[200]=aryuZZ[33]*aryuZZ[200];
   aryuZZ[122]=1./24.*aryuZZ[200] + aryuZZ[122] + aryuZZ[171];
   aryuZZ[122]=aryuZZ[4]*aryuZZ[122];
   aryuZZ[171]=1./3.*aryuZZ[17];
   aryuZZ[200]=aryuZZ[287] + aryuZZ[171];
   aryuZZ[200]=aryuZZ[10]*aryuZZ[200];
   aryuZZ[206]=7 + aryuZZ[118];
   aryuZZ[206]=aryuZZ[19]*aryuZZ[206];
   aryuZZ[209]=7./3. - aryuZZ[11];
   aryuZZ[209]=aryuZZ[17]*aryuZZ[209];
   aryuZZ[233]= - aryuZZ[11] - 1./2.*aryuZZ[10];
   aryuZZ[233]=aryuZZ[18]*aryuZZ[233];
   aryuZZ[191]=1./6.*aryuZZ[233] + aryuZZ[200] + 1./2.*aryuZZ[209] + 1./
   3.*aryuZZ[206] + aryuZZ[191] - aryuZZ[92] + 1./2.*aryuZZ[76];
   aryuZZ[191]=aryuZZ[73]*aryuZZ[191];
   aryuZZ[200]=aryuZZ[216] - 1 - aryuZZ[89];
   aryuZZ[200]=1./2.*aryuZZ[200] + 1./3.*aryuZZ[278];
   aryuZZ[200]=MMH*aryuZZ[73]*aryuZZ[200];
   aryuZZ[191]=aryuZZ[191] + aryuZZ[200];
   aryuZZ[191]=MMH*aryuZZ[191];
   aryuZZ[200]=1./3.*aryuZZ[18];
   aryuZZ[206]=aryuZZ[200] + aryuZZ[263] - aryuZZ[17];
   aryuZZ[206]=aryuZZ[18]*aryuZZ[206];
   aryuZZ[209]=pow(aryuZZ[19],2);
   aryuZZ[233]= - 1./3.*aryuZZ[209];
   aryuZZ[236]=11./6.*aryuZZ[19] - aryuZZ[17];
   aryuZZ[236]=aryuZZ[17]*aryuZZ[236];
   aryuZZ[206]=1./4.*aryuZZ[206] + aryuZZ[233] + 1./2.*aryuZZ[236];
   aryuZZ[206]=aryuZZ[73]*aryuZZ[206];
   aryuZZ[191]=aryuZZ[206] + aryuZZ[191];
   aryuZZ[191]=aryuZZ[272]*MMH*aryuZZ[191];
   aryuZZ[206]=9./4.*aryuZZ[43];
   aryuZZ[236]=1./2.*aryuZZ[47];
   aryuZZ[251]=1./4.*aryuZZ[45];
   aryuZZ[278]=aryuZZ[138] + aryuZZ[251] + aryuZZ[236] + 31./3. + 
   aryuZZ[206];
   aryuZZ[278]=aryuZZ[17]*aryuZZ[278];
   aryuZZ[285]= - 9*aryuZZ[43];
   aryuZZ[287]=35./18. + aryuZZ[285];
   aryuZZ[289]= - 1./2.*aryuZZ[45];
   aryuZZ[287]=aryuZZ[289] + 1./2.*aryuZZ[287] - aryuZZ[47];
   aryuZZ[290]=1./6.*aryuZZ[10];
   aryuZZ[287]=aryuZZ[290] + 1./2.*aryuZZ[287] + aryuZZ[151];
   aryuZZ[287]=aryuZZ[18]*aryuZZ[287];
   aryuZZ[291]=47 - 13./4.*aryuZZ[11];
   aryuZZ[291]=aryuZZ[19]*aryuZZ[291];
   aryuZZ[293]= - 47./2.*aryuZZ[19] + aryuZZ[17];
   aryuZZ[293]=aryuZZ[10]*aryuZZ[293];
   aryuZZ[278]=1./8.*aryuZZ[287] + 1./72.*aryuZZ[293] + 1./8.*
   aryuZZ[278] + 1./18.*aryuZZ[291] - 7./96.*aryuZZ[21] + 23./48.*
   aryuZZ[20] - 5./8.*aryuZZ[22] - 1./96.*aryuZZ[75] + 1./16.*
   aryuZZ[93] + 7./16.*aryuZZ[76] + 1./32.*aryuZZ[91] - 5./32.*
   aryuZZ[74] - aryuZZ[92];
   aryuZZ[278]=aryuZZ[73]*aryuZZ[278];
   aryuZZ[287]=9./8.*aryuZZ[43];
   aryuZZ[291]=1./4.*aryuZZ[47];
   aryuZZ[293]=1./8.*aryuZZ[45];
   aryuZZ[294]=aryuZZ[293] + aryuZZ[291] + 1./9. + aryuZZ[287];
   aryuZZ[294]= - 11./72.*aryuZZ[10] + 1./2.*aryuZZ[294] + aryuZZ[331];
   aryuZZ[294]=aryuZZ[10]*aryuZZ[294];
   aryuZZ[295]= - 7./3.*aryuZZ[87];
   aryuZZ[297]= - 3*aryuZZ[86];
   aryuZZ[298]=aryuZZ[297] - 33./2. + aryuZZ[295];
   aryuZZ[299]= - 1./2.*aryuZZ[96];
   aryuZZ[298]=aryuZZ[299] + 1./4.*aryuZZ[298] - aryuZZ[84];
   aryuZZ[294]=aryuZZ[294] + 1./36.*aryuZZ[11] + 17./24.*aryuZZ[15] + 1.
   /4.*aryuZZ[298] - aryuZZ[89];
   aryuZZ[294]=aryuZZ[73]*aryuZZ[294];
   aryuZZ[298]=aryuZZ[79] + 1./3.*aryuZZ[78];
   aryuZZ[298]=MMH*aryuZZ[73]*aryuZZ[298];
   aryuZZ[294]=aryuZZ[294] + 1./8.*aryuZZ[298];
   aryuZZ[294]=MMH*aryuZZ[294];
   aryuZZ[278]=aryuZZ[278] + 1./2.*aryuZZ[294];
   aryuZZ[278]=MMH*aryuZZ[278];
   aryuZZ[294]=aryuZZ[17]*aryuZZ[19];
   aryuZZ[300]= - 41./9.*aryuZZ[209] + 3./2.*aryuZZ[294];
   aryuZZ[301]=3./8.*aryuZZ[19] - 5./9.*aryuZZ[17];
   aryuZZ[301]=1./2.*aryuZZ[301] + 2./9.*aryuZZ[18];
   aryuZZ[301]=aryuZZ[18]*aryuZZ[301];
   aryuZZ[300]=1./8.*aryuZZ[300] + aryuZZ[301];
   aryuZZ[300]=aryuZZ[73]*aryuZZ[300];
   aryuZZ[278]=aryuZZ[300] + aryuZZ[278];
   aryuZZ[278]=MMH*aryuZZ[278];
   aryuZZ[191]=aryuZZ[278] + 1./12.*aryuZZ[191];
   aryuZZ[191]=aryuZZ[272]*aryuZZ[191];
   aryuZZ[278]=9./2.*aryuZZ[43];
   aryuZZ[300]=1./2.*aryuZZ[45];
   aryuZZ[301]=aryuZZ[300] + aryuZZ[47] + 5./9. + aryuZZ[278];
   aryuZZ[302]= - 7./18.*aryuZZ[10];
   aryuZZ[301]=aryuZZ[302] + 1./2.*aryuZZ[301] + aryuZZ[331];
   aryuZZ[301]=aryuZZ[10]*aryuZZ[301];
   aryuZZ[295]=aryuZZ[297] - 5./2. + aryuZZ[295];
   aryuZZ[304]=17./6.*aryuZZ[15];
   aryuZZ[295]=aryuZZ[301] + aryuZZ[304] - 1./2.*aryuZZ[89] + 
   aryuZZ[299] + 1./4.*aryuZZ[295] - aryuZZ[84];
   aryuZZ[295]=aryuZZ[73]*aryuZZ[295];
   aryuZZ[295]=aryuZZ[295] + 1./2.*aryuZZ[298];
   aryuZZ[295]=MMH*aryuZZ[295];
   aryuZZ[301]=1./2.*aryuZZ[93];
   aryuZZ[309]= - 1./12.*aryuZZ[75];
   aryuZZ[310]= - 7./12.*aryuZZ[21] + 23./6.*aryuZZ[20] - 37./6.*
   aryuZZ[22] + aryuZZ[309] + aryuZZ[301] + 1./4.*aryuZZ[91] - 5./4.*
   aryuZZ[74] - aryuZZ[92];
   aryuZZ[311]=aryuZZ[300] + aryuZZ[47] + 59./9. + aryuZZ[278];
   aryuZZ[311]=1./2.*aryuZZ[311] + aryuZZ[331];
   aryuZZ[311]=aryuZZ[17]*aryuZZ[311];
   aryuZZ[313]=49./4. - 19./3.*aryuZZ[11];
   aryuZZ[313]=aryuZZ[19]*aryuZZ[313];
   aryuZZ[310]=1./2.*aryuZZ[311] + 1./2.*aryuZZ[310] + 1./3.*
   aryuZZ[313];
   aryuZZ[311]=31./18. + aryuZZ[285];
   aryuZZ[311]=aryuZZ[289] + 1./2.*aryuZZ[311] - aryuZZ[47];
   aryuZZ[311]=7./9.*aryuZZ[10] + 1./2.*aryuZZ[311] + aryuZZ[138];
   aryuZZ[311]=aryuZZ[18]*aryuZZ[311];
   aryuZZ[313]= - aryuZZ[19] - 1./8.*aryuZZ[17];
   aryuZZ[313]=aryuZZ[10]*aryuZZ[313];
   aryuZZ[310]=1./4.*aryuZZ[311] + 1./2.*aryuZZ[310] + 1./9.*
   aryuZZ[313];
   aryuZZ[310]=aryuZZ[73]*aryuZZ[310];
   aryuZZ[295]=aryuZZ[310] + 1./4.*aryuZZ[295];
   aryuZZ[295]=MMH*aryuZZ[295];
   aryuZZ[310]=3*aryuZZ[19];
   aryuZZ[171]=aryuZZ[310] + aryuZZ[171];
   aryuZZ[171]=aryuZZ[17]*aryuZZ[171];
   aryuZZ[171]= - 53./9.*aryuZZ[209] + aryuZZ[171];
   aryuZZ[311]= - 1./2.*aryuZZ[19] - aryuZZ[17];
   aryuZZ[311]=5*aryuZZ[311] + 17./4.*aryuZZ[18];
   aryuZZ[311]=aryuZZ[18]*aryuZZ[311];
   aryuZZ[171]=1./4.*aryuZZ[171] + 1./9.*aryuZZ[311];
   aryuZZ[171]=aryuZZ[73]*aryuZZ[171];
   aryuZZ[171]=1./2.*aryuZZ[171] + aryuZZ[295];
   aryuZZ[171]=MMH*aryuZZ[171];
   aryuZZ[191]=aryuZZ[171] + aryuZZ[191];
   aryuZZ[191]=aryuZZ[272]*aryuZZ[191];
   aryuZZ[295]=1./6.*aryuZZ[11] + aryuZZ[293] + aryuZZ[291] - 41./9. + 
   aryuZZ[287];
   aryuZZ[295]=aryuZZ[17]*aryuZZ[295];
   aryuZZ[311]=1./2. - 3*aryuZZ[43];
   aryuZZ[311]=aryuZZ[289] + 3./2.*aryuZZ[311] - aryuZZ[47];
   aryuZZ[311]=aryuZZ[362] + 1./2.*aryuZZ[311] + aryuZZ[331];
   aryuZZ[311]=aryuZZ[18]*aryuZZ[311];
   aryuZZ[313]= - 7*aryuZZ[11];
   aryuZZ[315]= - 101./9. + aryuZZ[313];
   aryuZZ[315]=aryuZZ[19]*aryuZZ[315];
   aryuZZ[316]= - 137./4.*aryuZZ[19] - aryuZZ[17];
   aryuZZ[316]=aryuZZ[10]*aryuZZ[316];
   aryuZZ[295]=1./2.*aryuZZ[311] + 1./9.*aryuZZ[316] + aryuZZ[295] + 1./
   2.*aryuZZ[315] - 7./24.*aryuZZ[21] + 23./12.*aryuZZ[20] - 17./3.*
   aryuZZ[22] - 1./24.*aryuZZ[75] + 1./4.*aryuZZ[93] - 37./12.*
   aryuZZ[76] + 1./8.*aryuZZ[91] - 5./8.*aryuZZ[74] + 17./3.*aryuZZ[92]
   ;
   aryuZZ[295]=aryuZZ[73]*aryuZZ[295];
   aryuZZ[311]=19./2. - aryuZZ[87];
   aryuZZ[311]=7./3.*aryuZZ[311] + aryuZZ[297];
   aryuZZ[304]=aryuZZ[138] + aryuZZ[304] + 17./3.*aryuZZ[89] + 
   aryuZZ[299] + 1./4.*aryuZZ[311] - aryuZZ[84];
   aryuZZ[311]=aryuZZ[251] + aryuZZ[236] + 1./3. + aryuZZ[206];
   aryuZZ[170]=aryuZZ[170] + 1./2.*aryuZZ[311] + aryuZZ[138];
   aryuZZ[170]=aryuZZ[10]*aryuZZ[170];
   aryuZZ[170]=1./2.*aryuZZ[304] + aryuZZ[170];
   aryuZZ[170]=aryuZZ[73]*aryuZZ[170];
   aryuZZ[170]=aryuZZ[170] + 1./4.*aryuZZ[298];
   aryuZZ[170]=MMH*aryuZZ[170];
   aryuZZ[170]=aryuZZ[295] + aryuZZ[170];
   aryuZZ[170]=MMH*aryuZZ[170];
   aryuZZ[295]= - 23./2.*aryuZZ[19] + aryuZZ[17];
   aryuZZ[295]=aryuZZ[17]*aryuZZ[295];
   aryuZZ[298]=121./2.*aryuZZ[19] + aryuZZ[18];
   aryuZZ[298]=aryuZZ[18]*aryuZZ[298];
   aryuZZ[295]=1./24.*aryuZZ[298] - 2*aryuZZ[209] + 1./8.*aryuZZ[295];
   aryuZZ[295]=aryuZZ[73]*aryuZZ[295];
   aryuZZ[170]=1./3.*aryuZZ[295] + 1./4.*aryuZZ[170];
   aryuZZ[170]=MMH*aryuZZ[170];
   aryuZZ[295]= - 11*aryuZZ[17];
   aryuZZ[298]=19*aryuZZ[18] - 49*aryuZZ[19] + aryuZZ[295];
   aryuZZ[298]=aryuZZ[18]*aryuZZ[298];
   aryuZZ[298]=1./6.*aryuZZ[298] + 5*aryuZZ[209] + 11./6.*aryuZZ[294];
   aryuZZ[298]=aryuZZ[73]*aryuZZ[298];
   aryuZZ[304]=1./3. - 5*aryuZZ[11];
   aryuZZ[304]=aryuZZ[19]*aryuZZ[304];
   aryuZZ[311]= - 1./3.*aryuZZ[17];
   aryuZZ[315]=17./2.*aryuZZ[19] + aryuZZ[311];
   aryuZZ[315]=aryuZZ[10]*aryuZZ[315];
   aryuZZ[281]=aryuZZ[315] + aryuZZ[304] + 1./3.*aryuZZ[281];
   aryuZZ[304]= - 1./3. + aryuZZ[11];
   aryuZZ[304]=1./3.*aryuZZ[304] - 3./2.*aryuZZ[10];
   aryuZZ[304]=aryuZZ[18]*aryuZZ[304];
   aryuZZ[281]=1./3.*aryuZZ[281] + aryuZZ[304];
   aryuZZ[281]=aryuZZ[73]*aryuZZ[281];
   aryuZZ[304]=aryuZZ[10] + 1./4. - aryuZZ[11];
   aryuZZ[304]=aryuZZ[10]*aryuZZ[304];
   aryuZZ[304]= - 1./4.*aryuZZ[11] + aryuZZ[304];
   aryuZZ[304]=MMH*aryuZZ[73]*aryuZZ[304];
   aryuZZ[281]=1./4.*aryuZZ[281] + 1./9.*aryuZZ[304];
   aryuZZ[281]=MMH*aryuZZ[281];
   aryuZZ[281]=1./12.*aryuZZ[298] + aryuZZ[281];
   aryuZZ[281]=aryuZZ[133]*MMH*aryuZZ[281];
   aryuZZ[170]=aryuZZ[170] + 1./2.*aryuZZ[281];
   aryuZZ[170]=aryuZZ[133]*aryuZZ[170];
   aryuZZ[170]=aryuZZ[171] + aryuZZ[170];
   aryuZZ[170]=aryuZZ[133]*aryuZZ[170];
   aryuZZ[105]=aryuZZ[122] + aryuZZ[105] + aryuZZ[191] + aryuZZ[170];
   aryuZZ[105]=aryuZZ[4]*aryuZZ[105];
   aryuZZ[122]= - 739./64. - aryuZZ[88];
   aryuZZ[170]= - 5./4.*aryuZZ[97];
   aryuZZ[171]= - 1./6.*aryuZZ[87];
   aryuZZ[191]=7./4.*aryuZZ[86];
   aryuZZ[122]=aryuZZ[191] + aryuZZ[171] - 13./4.*aryuZZ[90] + 
   aryuZZ[170] + 1./3.*aryuZZ[122] + 1./4.*aryuZZ[85];
   aryuZZ[281]=53./27. + aryuZZ[285];
   aryuZZ[298]= - 3./4.*aryuZZ[45];
   aryuZZ[304]=13./18.*aryuZZ[12];
   aryuZZ[315]= - 7./18.*aryuZZ[11];
   aryuZZ[281]=aryuZZ[315] + 19./36.*aryuZZ[13] + aryuZZ[304] + 
   aryuZZ[298] + 1./4.*aryuZZ[281] - aryuZZ[47];
   aryuZZ[316]=4./9.*aryuZZ[10];
   aryuZZ[281]=1./2.*aryuZZ[281] + aryuZZ[316];
   aryuZZ[281]=aryuZZ[10]*aryuZZ[281];
   aryuZZ[317]=2./3.*aryuZZ[84];
   aryuZZ[318]=7./24.*aryuZZ[96];
   aryuZZ[319]=19./24.*aryuZZ[16];
   aryuZZ[320]=pow(aryuZZ[13],2);
   aryuZZ[321]=1./12.*aryuZZ[320];
   aryuZZ[322]=49./36. + aryuZZ[103];
   aryuZZ[322]=1./3.*aryuZZ[11]*aryuZZ[322];
   aryuZZ[122]=aryuZZ[281] + aryuZZ[322] + aryuZZ[321] + aryuZZ[293] - 
   15./8.*aryuZZ[15] + aryuZZ[291] + aryuZZ[287] + aryuZZ[319] + 3./8.*
   aryuZZ[89] + aryuZZ[318] + 1./2.*aryuZZ[122] + aryuZZ[317];
   aryuZZ[122]=aryuZZ[73]*aryuZZ[122];
   aryuZZ[281]=aryuZZ[1]*aryuZZ[11];
   aryuZZ[323]=aryuZZ[39]*aryuZZ[11];
   aryuZZ[324]=aryuZZ[281] + 1./3.*aryuZZ[323];
   aryuZZ[325]=1./3. - aryuZZ[7];
   aryuZZ[326]=aryuZZ[1]*aryuZZ[325];
   aryuZZ[325]=aryuZZ[39]*aryuZZ[325];
   aryuZZ[325]=1./3.*aryuZZ[326] + aryuZZ[325];
   aryuZZ[325]=aryuZZ[10]*aryuZZ[325];
   aryuZZ[326]= - aryuZZ[1]*aryuZZ[11];
   aryuZZ[328]= - aryuZZ[39]*aryuZZ[11];
   aryuZZ[329]=aryuZZ[326] + 1./3.*aryuZZ[328];
   aryuZZ[329]=aryuZZ[6]*aryuZZ[329];
   aryuZZ[324]=1./3.*aryuZZ[329] + 1./9.*aryuZZ[324] + 1./2.*
   aryuZZ[325];
   aryuZZ[324]=1./6.*aryuZZ[324];
   aryuZZ[329]=aryuZZ[80] - 1./2.*aryuZZ[81];
   aryuZZ[329]=1./2.*aryuZZ[329] + aryuZZ[78];
   aryuZZ[329]=1./3.*MMH*aryuZZ[73]*aryuZZ[329];
   aryuZZ[122]=aryuZZ[329] + aryuZZ[324] + aryuZZ[122];
   aryuZZ[122]=MMH*aryuZZ[122];
   aryuZZ[333]=9*aryuZZ[43];
   aryuZZ[335]= - aryuZZ[45] + 151./27. + aryuZZ[333];
   aryuZZ[336]=13./9.*aryuZZ[12];
   aryuZZ[335]=aryuZZ[216] - 37./9.*aryuZZ[13] + 1./2.*aryuZZ[335] + 
   aryuZZ[336];
   aryuZZ[335]=aryuZZ[17]*aryuZZ[335];
   aryuZZ[337]= - 13./9.*aryuZZ[12];
   aryuZZ[338]=5./9.*aryuZZ[11];
   aryuZZ[339]=aryuZZ[302] + aryuZZ[338] - 7./18.*aryuZZ[13] + 
   aryuZZ[337] + aryuZZ[45] - 329./54. + aryuZZ[47];
   aryuZZ[339]=aryuZZ[18]*aryuZZ[339];
   aryuZZ[340]=11*aryuZZ[92];
   aryuZZ[342]= - 3*aryuZZ[74];
   aryuZZ[343]=1./6.*aryuZZ[77] + aryuZZ[342] + aryuZZ[340];
   aryuZZ[344]=9*aryuZZ[13];
   aryuZZ[348]= - 1847./27. + aryuZZ[344];
   aryuZZ[350]=53./9.*aryuZZ[11];
   aryuZZ[348]=1./2.*aryuZZ[348] + aryuZZ[350];
   aryuZZ[348]=aryuZZ[19]*aryuZZ[348];
   aryuZZ[352]= - 7./12.*aryuZZ[75];
   aryuZZ[354]=19./4.*aryuZZ[20];
   aryuZZ[358]=25./6.*aryuZZ[19] + aryuZZ[17];
   aryuZZ[358]=1./3.*aryuZZ[10]*aryuZZ[358];
   aryuZZ[335]=1./2.*aryuZZ[339] + aryuZZ[358] + 1./2.*aryuZZ[335] + 1./
   2.*aryuZZ[348] + 37./12.*aryuZZ[21] + aryuZZ[354] + 103./12.*
   aryuZZ[22] + aryuZZ[352] - aryuZZ[93] - 29./12.*aryuZZ[76] + 1./2.*
   aryuZZ[343] - aryuZZ[91];
   aryuZZ[335]=aryuZZ[73]*aryuZZ[335];
   aryuZZ[339]=1 - aryuZZ[7];
   aryuZZ[339]=aryuZZ[17]*aryuZZ[339];
   aryuZZ[339]=aryuZZ[186] + 1./2.*aryuZZ[339];
   aryuZZ[339]=aryuZZ[1]*aryuZZ[339];
   aryuZZ[343]= - 1./3. + aryuZZ[7];
   aryuZZ[348]=aryuZZ[1]*aryuZZ[343];
   aryuZZ[343]=aryuZZ[39]*aryuZZ[343];
   aryuZZ[343]=1./3.*aryuZZ[348] + aryuZZ[343];
   aryuZZ[348]=aryuZZ[18]*aryuZZ[343];
   aryuZZ[360]=aryuZZ[1]*aryuZZ[356];
   aryuZZ[356]=aryuZZ[39]*aryuZZ[356];
   aryuZZ[361]=aryuZZ[360] + 1./3.*aryuZZ[356];
   aryuZZ[361]=aryuZZ[6]*aryuZZ[361];
   aryuZZ[362]=11./27. - aryuZZ[7];
   aryuZZ[362]=aryuZZ[17]*aryuZZ[362];
   aryuZZ[362]= - 1./27.*aryuZZ[19] + 1./2.*aryuZZ[362];
   aryuZZ[362]=aryuZZ[39]*aryuZZ[362];
   aryuZZ[339]=1./3.*aryuZZ[361] + 1./2.*aryuZZ[348] + 1./3.*
   aryuZZ[339] + aryuZZ[362];
   aryuZZ[339]=1./3.*aryuZZ[339];
   aryuZZ[335]=aryuZZ[339] + aryuZZ[335];
   aryuZZ[122]=1./2.*aryuZZ[335] + aryuZZ[122];
   aryuZZ[122]=MMH*aryuZZ[122];
   aryuZZ[335]= - 3./2.*aryuZZ[45];
   aryuZZ[348]=1./9.*aryuZZ[13] - 17./18.*aryuZZ[12] + aryuZZ[335] - 
   aryuZZ[47] - 29./27. - 9./2.*aryuZZ[43];
   aryuZZ[361]= - 7./9.*aryuZZ[11];
   aryuZZ[348]=13./9.*aryuZZ[10] + 1./2.*aryuZZ[348] + aryuZZ[361];
   aryuZZ[348]=aryuZZ[10]*aryuZZ[348];
   aryuZZ[362]=7./2.*aryuZZ[86] - 1./3.*aryuZZ[87] - aryuZZ[90] - 2653./
   288. + 5*aryuZZ[97];
   aryuZZ[363]= - 5./9. + aryuZZ[240];
   aryuZZ[364]=aryuZZ[11]*aryuZZ[363];
   aryuZZ[365]=1./3.*aryuZZ[84];
   aryuZZ[366]=9./16.*aryuZZ[43];
   aryuZZ[367]=1./8.*aryuZZ[47];
   aryuZZ[368]=1./16.*aryuZZ[45];
   aryuZZ[320]= - 1./24.*aryuZZ[320];
   aryuZZ[348]=1./4.*aryuZZ[348] + 1./24.*aryuZZ[364] + aryuZZ[320] + 
   aryuZZ[368] - 41./48.*aryuZZ[15] + aryuZZ[367] + aryuZZ[366] - 5./8.
   *aryuZZ[16] + 37./72.*aryuZZ[89] + 1./8.*aryuZZ[96] + 1./8.*
   aryuZZ[362] + aryuZZ[365];
   aryuZZ[348]=aryuZZ[73]*aryuZZ[348];
   aryuZZ[362]= - aryuZZ[14] + aryuZZ[216];
   aryuZZ[364]=aryuZZ[1]*aryuZZ[362];
   aryuZZ[362]=aryuZZ[39]*aryuZZ[362];
   aryuZZ[369]=aryuZZ[10]*aryuZZ[230];
   aryuZZ[370]=1./3.*aryuZZ[369];
   aryuZZ[362]=aryuZZ[370] + aryuZZ[364] + 11./9.*aryuZZ[362];
   aryuZZ[364]=aryuZZ[326] + 11./9.*aryuZZ[328];
   aryuZZ[369]=1./2.*aryuZZ[364] + aryuZZ[369];
   aryuZZ[369]=aryuZZ[6]*aryuZZ[369];
   aryuZZ[362]=1./2.*aryuZZ[362] + aryuZZ[369];
   aryuZZ[369]=1./4.*aryuZZ[81];
   aryuZZ[371]=aryuZZ[369] + aryuZZ[78];
   aryuZZ[371]=MMH*aryuZZ[73]*aryuZZ[371];
   aryuZZ[348]=1./6.*aryuZZ[371] + 1./6.*aryuZZ[362] + aryuZZ[348];
   aryuZZ[348]=MMH*aryuZZ[348];
   aryuZZ[362]=167./9. - 29./2.*aryuZZ[13];
   aryuZZ[341]=1./2.*aryuZZ[362] + aryuZZ[341];
   aryuZZ[341]=aryuZZ[19]*aryuZZ[341];
   aryuZZ[362]=77./36.*aryuZZ[13] - 17./36.*aryuZZ[12] - 1./4.*
   aryuZZ[45] + aryuZZ[236] + 97./27. + aryuZZ[206];
   aryuZZ[362]=1./2.*aryuZZ[362] + aryuZZ[216];
   aryuZZ[362]=aryuZZ[17]*aryuZZ[362];
   aryuZZ[371]=aryuZZ[342] + 13./9.*aryuZZ[92];
   aryuZZ[371]= - 7./6.*aryuZZ[76] - aryuZZ[91] + 1./2.*aryuZZ[371] + 1.
   /3.*aryuZZ[77];
   aryuZZ[372]=aryuZZ[10]*aryuZZ[19];
   aryuZZ[373]=25./9.*aryuZZ[10] - 5./9.*aryuZZ[11] + 5./9.*aryuZZ[13]
    + 17./18.*aryuZZ[12] - 83./54. + aryuZZ[45];
   aryuZZ[373]=aryuZZ[18]*aryuZZ[373];
   aryuZZ[374]= - 1./3.*aryuZZ[93];
   aryuZZ[197]=1./4.*aryuZZ[373] + 23./9.*aryuZZ[372] + aryuZZ[362] + 1.
   /6.*aryuZZ[341] + aryuZZ[197] + 115./24.*aryuZZ[20] + 25./36.*
   aryuZZ[22] - 7./24.*aryuZZ[75] + 1./2.*aryuZZ[371] + aryuZZ[374];
   aryuZZ[197]=aryuZZ[73]*aryuZZ[197];
   aryuZZ[341]= - aryuZZ[1] - 11./9.*aryuZZ[39];
   aryuZZ[362]=aryuZZ[18]*aryuZZ[341];
   aryuZZ[371]=aryuZZ[1]*aryuZZ[178];
   aryuZZ[178]=aryuZZ[39]*aryuZZ[178];
   aryuZZ[373]=aryuZZ[362] + aryuZZ[371] + 11./9.*aryuZZ[178];
   aryuZZ[375]=aryuZZ[1]*aryuZZ[346];
   aryuZZ[346]=aryuZZ[39]*aryuZZ[346];
   aryuZZ[346]=aryuZZ[375] + 11./9.*aryuZZ[346];
   aryuZZ[346]=1./2.*aryuZZ[346] + aryuZZ[362];
   aryuZZ[346]=aryuZZ[6]*aryuZZ[346];
   aryuZZ[346]=1./6.*aryuZZ[373] + aryuZZ[346];
   aryuZZ[197]=1./3.*aryuZZ[346] + aryuZZ[197];
   aryuZZ[197]=1./2.*aryuZZ[197] + aryuZZ[348];
   aryuZZ[197]=MMH*aryuZZ[197];
   aryuZZ[346]=aryuZZ[281] + 11./9.*aryuZZ[323];
   aryuZZ[348]=aryuZZ[6]*aryuZZ[364];
   aryuZZ[364]=aryuZZ[11]*aryuZZ[246];
   aryuZZ[364]=1./3.*aryuZZ[364] - aryuZZ[16] + aryuZZ[89] + 1./4. + 
   aryuZZ[97];
   aryuZZ[373]= - 1./12.*aryuZZ[12] - aryuZZ[11];
   aryuZZ[373]=aryuZZ[10]*aryuZZ[373];
   aryuZZ[364]=1./2.*aryuZZ[364] + aryuZZ[373];
   aryuZZ[364]=aryuZZ[73]*aryuZZ[364];
   aryuZZ[346]=aryuZZ[364] + 1./3.*aryuZZ[346] + aryuZZ[348];
   aryuZZ[346]=MMH*aryuZZ[346];
   aryuZZ[348]=aryuZZ[1]*aryuZZ[359];
   aryuZZ[359]=aryuZZ[39]*aryuZZ[359];
   aryuZZ[348]=aryuZZ[348] + 11./9.*aryuZZ[359];
   aryuZZ[356]=aryuZZ[360] + 11./9.*aryuZZ[356];
   aryuZZ[356]=aryuZZ[6]*aryuZZ[356];
   aryuZZ[348]=1./3.*aryuZZ[348] + aryuZZ[356];
   aryuZZ[356]=aryuZZ[263] + aryuZZ[131];
   aryuZZ[356]=aryuZZ[10]*aryuZZ[356];
   aryuZZ[246]=aryuZZ[19]*aryuZZ[246];
   aryuZZ[246]=5./3.*aryuZZ[246] - aryuZZ[22] + aryuZZ[20];
   aryuZZ[359]=5./16.*aryuZZ[13] - 1./3. - 1./16.*aryuZZ[12];
   aryuZZ[359]=aryuZZ[17]*aryuZZ[359];
   aryuZZ[246]=1./2.*aryuZZ[356] + 1./8.*aryuZZ[246] + 1./3.*
   aryuZZ[359];
   aryuZZ[356]=1 + 1./2.*aryuZZ[12];
   aryuZZ[356]=1./9.*aryuZZ[356] - aryuZZ[11];
   aryuZZ[356]=1./2.*aryuZZ[356] + 1./9.*aryuZZ[10];
   aryuZZ[356]=aryuZZ[18]*aryuZZ[356];
   aryuZZ[246]=1./3.*aryuZZ[246] + 1./4.*aryuZZ[356];
   aryuZZ[246]=aryuZZ[73]*aryuZZ[246];
   aryuZZ[246]=1./12.*aryuZZ[346] + 1./12.*aryuZZ[348] + aryuZZ[246];
   aryuZZ[246]=MMH*aryuZZ[246];
   aryuZZ[346]= - aryuZZ[17] - 1./12.*aryuZZ[18];
   aryuZZ[346]=aryuZZ[18]*aryuZZ[346];
   aryuZZ[233]=1./2.*aryuZZ[346] + aryuZZ[233] + 1./8.*aryuZZ[294];
   aryuZZ[233]=aryuZZ[73]*aryuZZ[233];
   aryuZZ[233]=1./3.*aryuZZ[233] + aryuZZ[246];
   aryuZZ[233]=aryuZZ[272]*aryuZZ[233];
   aryuZZ[246]= - 91./18.*aryuZZ[19] - aryuZZ[17];
   aryuZZ[246]=aryuZZ[17]*aryuZZ[246];
   aryuZZ[346]= - 97./6.*aryuZZ[18] - 235./6.*aryuZZ[19] - 7*aryuZZ[17]
   ;
   aryuZZ[346]=aryuZZ[18]*aryuZZ[346];
   aryuZZ[246]=1./3.*aryuZZ[346] - 233./18.*aryuZZ[209] + aryuZZ[246];
   aryuZZ[246]=aryuZZ[73]*aryuZZ[246];
   aryuZZ[197]=aryuZZ[233] + 1./4.*aryuZZ[246] + aryuZZ[197];
   aryuZZ[197]=aryuZZ[272]*aryuZZ[197];
   aryuZZ[233]= - 1453./3.*aryuZZ[19] + 7*aryuZZ[17];
   aryuZZ[246]= - 17*aryuZZ[18];
   aryuZZ[233]=1./3.*aryuZZ[233] + aryuZZ[246];
   aryuZZ[233]=aryuZZ[18]*aryuZZ[233];
   aryuZZ[346]=29./6.*aryuZZ[19] - aryuZZ[17];
   aryuZZ[346]=aryuZZ[17]*aryuZZ[346];
   aryuZZ[233]=1./4.*aryuZZ[233] - 1471./36.*aryuZZ[209] + aryuZZ[346];
   aryuZZ[233]=aryuZZ[73]*aryuZZ[233];
   aryuZZ[122]=aryuZZ[197] + 1./2.*aryuZZ[233] + aryuZZ[122];
   aryuZZ[122]=aryuZZ[272]*aryuZZ[122];
   aryuZZ[197]= - 1./3.*aryuZZ[88];
   aryuZZ[170]=aryuZZ[191] + aryuZZ[171] - 43./12.*aryuZZ[90] + 
   aryuZZ[170] + 5./4.*aryuZZ[85] + 31./64. + aryuZZ[197];
   aryuZZ[171]=305./27. + aryuZZ[285];
   aryuZZ[171]=aryuZZ[315] + 67./36.*aryuZZ[13] + aryuZZ[304] + 
   aryuZZ[298] + 1./4.*aryuZZ[171] - aryuZZ[47];
   aryuZZ[171]=1./2.*aryuZZ[171] + aryuZZ[316];
   aryuZZ[171]=aryuZZ[10]*aryuZZ[171];
   aryuZZ[170]=aryuZZ[171] + aryuZZ[322] + aryuZZ[321] + aryuZZ[293] - 
   19./8.*aryuZZ[15] + aryuZZ[291] + aryuZZ[287] + aryuZZ[319] + 89./24.
   *aryuZZ[89] + aryuZZ[318] + 1./2.*aryuZZ[170] + aryuZZ[317];
   aryuZZ[170]=aryuZZ[73]*aryuZZ[170];
   aryuZZ[170]=aryuZZ[329] + aryuZZ[324] + aryuZZ[170];
   aryuZZ[170]=MMH*aryuZZ[170];
   aryuZZ[171]= - aryuZZ[45] - 929./27. + aryuZZ[333];
   aryuZZ[171]=aryuZZ[216] - 13./9.*aryuZZ[13] + 1./2.*aryuZZ[171] + 
   aryuZZ[336];
   aryuZZ[171]=aryuZZ[17]*aryuZZ[171];
   aryuZZ[191]=aryuZZ[302] + aryuZZ[338] - 55./18.*aryuZZ[13] + 
   aryuZZ[337] + aryuZZ[45] - 581./54. + aryuZZ[47];
   aryuZZ[191]=aryuZZ[18]*aryuZZ[191];
   aryuZZ[233]= - 4079./27. + aryuZZ[344];
   aryuZZ[233]=1./2.*aryuZZ[233] + aryuZZ[350];
   aryuZZ[233]=aryuZZ[19]*aryuZZ[233];
   aryuZZ[293]= - 11./6.*aryuZZ[77] + aryuZZ[342] + 89./3.*aryuZZ[92];
   aryuZZ[171]=1./2.*aryuZZ[191] + aryuZZ[358] + 1./2.*aryuZZ[171] + 1./
   2.*aryuZZ[233] + 49./12.*aryuZZ[21] + aryuZZ[354] + 79./12.*
   aryuZZ[22] + aryuZZ[352] - aryuZZ[93] - 61./12.*aryuZZ[76] + 1./2.*
   aryuZZ[293] - aryuZZ[91];
   aryuZZ[171]=aryuZZ[73]*aryuZZ[171];
   aryuZZ[171]=aryuZZ[339] + aryuZZ[171];
   aryuZZ[170]=1./2.*aryuZZ[171] + aryuZZ[170];
   aryuZZ[170]=MMH*aryuZZ[170];
   aryuZZ[171]= - 9./8.*aryuZZ[43];
   aryuZZ[191]= - 3./4.*aryuZZ[47];
   aryuZZ[233]= - 3./8.*aryuZZ[45];
   aryuZZ[290]=aryuZZ[290] - 23./6.*aryuZZ[13] + 55./24.*aryuZZ[12] + 
   aryuZZ[233] + aryuZZ[191] - 101./27. + aryuZZ[171];
   aryuZZ[290]=aryuZZ[10]*aryuZZ[290];
   aryuZZ[197]=7./8.*aryuZZ[86] - 1./12.*aryuZZ[87] - 13./6.*aryuZZ[97]
    - 13./12.*aryuZZ[85] - 1101./128. + aryuZZ[197];
   aryuZZ[293]= - 11./2.*aryuZZ[13];
   aryuZZ[302]= - 1./9. + aryuZZ[293];
   aryuZZ[302]=aryuZZ[11]*aryuZZ[302];
   aryuZZ[197]=1./2.*aryuZZ[290] + 5./24.*aryuZZ[302] + aryuZZ[320] + 
   aryuZZ[368] - 17./48.*aryuZZ[15] + aryuZZ[367] + aryuZZ[366] + 5./4.
   *aryuZZ[16] - 85./12.*aryuZZ[89] + 1./6.*aryuZZ[96] + 1./2.*
   aryuZZ[197] + aryuZZ[365];
   aryuZZ[197]=aryuZZ[73]*aryuZZ[197];
   aryuZZ[290]= - aryuZZ[14] + aryuZZ[151];
   aryuZZ[302]=aryuZZ[1]*aryuZZ[290];
   aryuZZ[290]=aryuZZ[39]*aryuZZ[290];
   aryuZZ[290]=1./3.*aryuZZ[302] + aryuZZ[290];
   aryuZZ[302]= - 1./2.*aryuZZ[7];
   aryuZZ[304]=1./3. + aryuZZ[302];
   aryuZZ[315]=aryuZZ[1]*aryuZZ[304];
   aryuZZ[304]=aryuZZ[39]*aryuZZ[304];
   aryuZZ[304]=1./3.*aryuZZ[315] + aryuZZ[304];
   aryuZZ[304]=aryuZZ[10]*aryuZZ[304];
   aryuZZ[281]=1./3.*aryuZZ[281] + aryuZZ[323];
   aryuZZ[281]=1./2.*aryuZZ[281];
   aryuZZ[316]=1./3.*aryuZZ[1] + aryuZZ[39];
   aryuZZ[317]=aryuZZ[10]*aryuZZ[316];
   aryuZZ[318]=aryuZZ[281] + aryuZZ[317];
   aryuZZ[318]=aryuZZ[6]*aryuZZ[318];
   aryuZZ[290]=aryuZZ[318] + 1./2.*aryuZZ[290] + aryuZZ[304];
   aryuZZ[304]=aryuZZ[78] + aryuZZ[80] + aryuZZ[369];
   aryuZZ[304]=MMH*aryuZZ[73]*aryuZZ[304];
   aryuZZ[197]=1./6.*aryuZZ[304] + 1./6.*aryuZZ[290] + aryuZZ[197];
   aryuZZ[197]=MMH*aryuZZ[197];
   aryuZZ[290]= - 211./6.*aryuZZ[13] + 55./6.*aryuZZ[12] + aryuZZ[289]
    - aryuZZ[47] + 2255./27. + aryuZZ[278];
   aryuZZ[290]=1./2.*aryuZZ[290] + aryuZZ[151];
   aryuZZ[290]=aryuZZ[17]*aryuZZ[290];
   aryuZZ[304]=aryuZZ[215] + aryuZZ[300] + 397./108. + aryuZZ[47];
   aryuZZ[304]= - 4./9.*aryuZZ[10] + 5./24.*aryuZZ[11] + 1./4.*
   aryuZZ[304] + 2*aryuZZ[13];
   aryuZZ[304]=aryuZZ[18]*aryuZZ[304];
   aryuZZ[318]=5*aryuZZ[77] + aryuZZ[342] - 169./3.*aryuZZ[92];
   aryuZZ[318]=25./3.*aryuZZ[76] + 1./2.*aryuZZ[318] - aryuZZ[91];
   aryuZZ[173]=1123./9. + aryuZZ[173];
   aryuZZ[173]=1./2.*aryuZZ[173] + 19*aryuZZ[11];
   aryuZZ[173]=aryuZZ[19]*aryuZZ[173];
   aryuZZ[319]=49./2.*aryuZZ[19] + aryuZZ[17];
   aryuZZ[319]=aryuZZ[10]*aryuZZ[319];
   aryuZZ[173]=aryuZZ[304] + 1./6.*aryuZZ[319] + 1./4.*aryuZZ[290] + 1./
   6.*aryuZZ[173] + aryuZZ[142] - 17./48.*aryuZZ[20] + 28./3.*
   aryuZZ[22] - 7./48.*aryuZZ[75] + 1./4.*aryuZZ[318] + aryuZZ[374];
   aryuZZ[173]=aryuZZ[73]*aryuZZ[173];
   aryuZZ[290]= - 1./3. + aryuZZ[302];
   aryuZZ[290]=aryuZZ[17]*aryuZZ[290];
   aryuZZ[150]=aryuZZ[150] + aryuZZ[290];
   aryuZZ[290]=aryuZZ[1]*aryuZZ[150];
   aryuZZ[150]=aryuZZ[39]*aryuZZ[150];
   aryuZZ[304]=1./2.*aryuZZ[7];
   aryuZZ[318]= - 1./3. + aryuZZ[304];
   aryuZZ[319]=aryuZZ[1]*aryuZZ[318];
   aryuZZ[318]=aryuZZ[39]*aryuZZ[318];
   aryuZZ[318]=1./3.*aryuZZ[319] + aryuZZ[318];
   aryuZZ[318]=aryuZZ[18]*aryuZZ[318];
   aryuZZ[150]=aryuZZ[318] + 1./3.*aryuZZ[290] + aryuZZ[150];
   aryuZZ[186]=aryuZZ[186] + aryuZZ[17];
   aryuZZ[290]=aryuZZ[1]*aryuZZ[186];
   aryuZZ[186]=aryuZZ[39]*aryuZZ[186];
   aryuZZ[186]=1./3.*aryuZZ[290] + aryuZZ[186];
   aryuZZ[290]=aryuZZ[18]*aryuZZ[139];
   aryuZZ[186]=1./2.*aryuZZ[186] + 1./3.*aryuZZ[290];
   aryuZZ[186]=aryuZZ[6]*aryuZZ[186];
   aryuZZ[150]=1./3.*aryuZZ[150] + aryuZZ[186];
   aryuZZ[150]=aryuZZ[197] + 1./2.*aryuZZ[150] + aryuZZ[173];
   aryuZZ[150]=MMH*aryuZZ[150];
   aryuZZ[173]=aryuZZ[17]*aryuZZ[7];
   aryuZZ[186]=aryuZZ[263] + aryuZZ[173];
   aryuZZ[197]=aryuZZ[1]*aryuZZ[186];
   aryuZZ[186]=aryuZZ[39]*aryuZZ[186];
   aryuZZ[263]= - 1./3. - aryuZZ[7];
   aryuZZ[318]=aryuZZ[1]*aryuZZ[263];
   aryuZZ[319]=aryuZZ[39]*aryuZZ[263];
   aryuZZ[318]=1./3.*aryuZZ[318] + aryuZZ[319];
   aryuZZ[319]=aryuZZ[18]*aryuZZ[318];
   aryuZZ[186]=aryuZZ[319] + 1./3.*aryuZZ[197] + aryuZZ[186];
   aryuZZ[178]=1./3.*aryuZZ[371] + aryuZZ[178];
   aryuZZ[197]=aryuZZ[18]*aryuZZ[316];
   aryuZZ[178]=1./2.*aryuZZ[178] + aryuZZ[197];
   aryuZZ[178]=aryuZZ[6]*aryuZZ[178];
   aryuZZ[178]=1./2.*aryuZZ[186] + aryuZZ[178];
   aryuZZ[186]= - aryuZZ[12] + aryuZZ[13];
   aryuZZ[316]=11./4.*aryuZZ[186] + aryuZZ[151];
   aryuZZ[316]=aryuZZ[17]*aryuZZ[316];
   aryuZZ[319]=13./12.*aryuZZ[11] + 1./27. + 11./16.*aryuZZ[13];
   aryuZZ[319]=aryuZZ[19]*aryuZZ[319];
   aryuZZ[320]= - 59./3.*aryuZZ[19] + aryuZZ[17];
   aryuZZ[320]=aryuZZ[10]*aryuZZ[320];
   aryuZZ[321]=41./72.*aryuZZ[10] - 1./72.*aryuZZ[11] - 11./8.*
   aryuZZ[13] - 1./27. + 11./16.*aryuZZ[12];
   aryuZZ[321]=aryuZZ[18]*aryuZZ[321];
   aryuZZ[316]=aryuZZ[321] + 1./12.*aryuZZ[320] + aryuZZ[319] + 1./4.*
   aryuZZ[316];
   aryuZZ[316]=aryuZZ[73]*aryuZZ[316];
   aryuZZ[319]=1./3.*aryuZZ[326] + aryuZZ[328];
   aryuZZ[320]=aryuZZ[1]*aryuZZ[141];
   aryuZZ[141]=aryuZZ[39]*aryuZZ[141];
   aryuZZ[321]=1./3.*aryuZZ[320] + aryuZZ[141];
   aryuZZ[321]=aryuZZ[10]*aryuZZ[321];
   aryuZZ[319]=1./3.*aryuZZ[319] + aryuZZ[321];
   aryuZZ[322]=aryuZZ[10]*aryuZZ[139];
   aryuZZ[281]=aryuZZ[281] + aryuZZ[322];
   aryuZZ[281]=aryuZZ[6]*aryuZZ[281];
   aryuZZ[281]=1./2.*aryuZZ[319] + aryuZZ[281];
   aryuZZ[319]= - 11./2.*aryuZZ[12];
   aryuZZ[323]=11*aryuZZ[13];
   aryuZZ[324]=aryuZZ[323] + 47./27. + aryuZZ[319];
   aryuZZ[117]=aryuZZ[117] + 1./4.*aryuZZ[324] + aryuZZ[338];
   aryuZZ[117]=aryuZZ[10]*aryuZZ[117];
   aryuZZ[324]= - 47./27. + aryuZZ[293];
   aryuZZ[324]=aryuZZ[11]*aryuZZ[324];
   aryuZZ[117]=1./4.*aryuZZ[324] + aryuZZ[117];
   aryuZZ[117]=aryuZZ[73]*aryuZZ[117];
   aryuZZ[117]=1./3.*aryuZZ[281] + aryuZZ[117];
   aryuZZ[117]=MMH*aryuZZ[117];
   aryuZZ[117]=1./2.*aryuZZ[117] + 1./6.*aryuZZ[178] + aryuZZ[316];
   aryuZZ[117]=MMH*aryuZZ[117];
   aryuZZ[178]= - 91*aryuZZ[209] + 17*aryuZZ[294];
   aryuZZ[281]= - 17*aryuZZ[17];
   aryuZZ[294]=281./3.*aryuZZ[19] + aryuZZ[281];
   aryuZZ[222]=1./8.*aryuZZ[294] + aryuZZ[222];
   aryuZZ[222]=aryuZZ[18]*aryuZZ[222];
   aryuZZ[178]=1./8.*aryuZZ[178] + aryuZZ[222];
   aryuZZ[178]=aryuZZ[73]*aryuZZ[178];
   aryuZZ[117]=1./3.*aryuZZ[178] + aryuZZ[117];
   aryuZZ[117]=aryuZZ[133]*aryuZZ[117];
   aryuZZ[178]= - 7./2.*aryuZZ[19] - aryuZZ[17];
   aryuZZ[178]=aryuZZ[17]*aryuZZ[178];
   aryuZZ[178]=1679./3.*aryuZZ[209] + aryuZZ[178];
   aryuZZ[121]=2633./3.*aryuZZ[19] + aryuZZ[121];
   aryuZZ[121]=1./4.*aryuZZ[121] - 7./9.*aryuZZ[18];
   aryuZZ[121]=aryuZZ[18]*aryuZZ[121];
   aryuZZ[121]=1./4.*aryuZZ[178] + aryuZZ[121];
   aryuZZ[121]=aryuZZ[73]*aryuZZ[121];
   aryuZZ[117]=aryuZZ[117] + aryuZZ[121] + aryuZZ[150];
   aryuZZ[117]=aryuZZ[133]*aryuZZ[117];
   aryuZZ[121]= - 2059./3.*aryuZZ[19] + aryuZZ[17];
   aryuZZ[121]=7./3.*aryuZZ[121] + aryuZZ[246];
   aryuZZ[121]=aryuZZ[18]*aryuZZ[121];
   aryuZZ[121]=1./4.*aryuZZ[121] - 10687./36.*aryuZZ[209] + aryuZZ[346]
   ;
   aryuZZ[121]=aryuZZ[73]*aryuZZ[121];
   aryuZZ[117]=aryuZZ[117] + 1./2.*aryuZZ[121] + aryuZZ[170];
   aryuZZ[117]=aryuZZ[133]*aryuZZ[117];
   aryuZZ[121]=aryuZZ[18]*aryuZZ[19];
   aryuZZ[121]=aryuZZ[209] + aryuZZ[121];
   aryuZZ[121]=aryuZZ[73]*aryuZZ[121];
   aryuZZ[104]=aryuZZ[105] + aryuZZ[104] + aryuZZ[117] + 32*aryuZZ[121]
    + aryuZZ[122];
   aryuZZ[104]=aryuZZ[4]*aryuZZ[104];
   aryuZZ[105]= - 22*aryuZZ[46];
   aryuZZ[117]= - 6*aryuZZ[45];
   aryuZZ[121]= - 26*aryuZZ[12];
   aryuZZ[122]= - 4./3. + aryuZZ[7];
   aryuZZ[122]=2*aryuZZ[1]*aryuZZ[122];
   aryuZZ[150]=3*aryuZZ[7];
   aryuZZ[170]= - 20./9. + aryuZZ[150];
   aryuZZ[170]=2*aryuZZ[39]*aryuZZ[170];
   aryuZZ[178]= - 12*aryuZZ[10];
   aryuZZ[222]=3*aryuZZ[1] + 11./3.*aryuZZ[39];
   aryuZZ[222]=2*aryuZZ[6]*aryuZZ[222];
   aryuZZ[246]= - 3*aryuZZ[48];
   aryuZZ[294]= - 6*aryuZZ[47];
   aryuZZ[316]=6*aryuZZ[11];
   aryuZZ[324]=aryuZZ[222] + aryuZZ[178] + aryuZZ[170] + aryuZZ[122] + 
   aryuZZ[316] - 17./2.*aryuZZ[13] + aryuZZ[121] + aryuZZ[117] + 
   aryuZZ[294] + aryuZZ[246] - 95./2. + aryuZZ[105];
   aryuZZ[324]=aryuZZ[37]*aryuZZ[324];
   aryuZZ[326]= - 9*aryuZZ[38];
   aryuZZ[328]=aryuZZ[275] + 7./3. + aryuZZ[326];
   aryuZZ[328]=aryuZZ[19]*aryuZZ[328];
   aryuZZ[329]= - 13*aryuZZ[17];
   aryuZZ[124]=aryuZZ[124] + aryuZZ[329];
   aryuZZ[124]=aryuZZ[35]*aryuZZ[124];
   aryuZZ[338]= - 11*aryuZZ[21] - 13./6.*aryuZZ[20] - 3*aryuZZ[22] - 25
   *aryuZZ[50] + 3*aryuZZ[42] + 11*aryuZZ[41];
   aryuZZ[339]=13./6.*aryuZZ[17];
   aryuZZ[344]= - 205./6.*aryuZZ[35] + 3./2.*aryuZZ[36] - 65./3. + 
   aryuZZ[326];
   aryuZZ[344]=1./2.*aryuZZ[18]*aryuZZ[344];
   aryuZZ[124]=aryuZZ[324] + aryuZZ[344] + 1./6.*aryuZZ[124] + 
   aryuZZ[339] + aryuZZ[338] + 1./2.*aryuZZ[328];
   aryuZZ[124]=aryuZZ[3]*aryuZZ[124];
   aryuZZ[324]=157./108. - aryuZZ[60];
   aryuZZ[324]=7*aryuZZ[324] - 271./8.*aryuZZ[63];
   aryuZZ[328]= - 22*aryuZZ[72];
   aryuZZ[346]=19./8.*aryuZZ[61];
   aryuZZ[348]=11./2.*aryuZZ[59];
   aryuZZ[350]= - 23*aryuZZ[64];
   aryuZZ[352]= - 61./6.*aryuZZ[67];
   aryuZZ[324]=101./32.*aryuZZ[48] + 1699./72.*aryuZZ[46] + aryuZZ[352]
    + aryuZZ[350] + aryuZZ[348] + 7./4.*aryuZZ[16] + 5771./288.*
   aryuZZ[49] - 181./192.*aryuZZ[66] + aryuZZ[346] - 313./96.*
   aryuZZ[65] - 809./54.*aryuZZ[62] + 1./4.*aryuZZ[324] + aryuZZ[328];
   aryuZZ[354]= - 3*aryuZZ[15];
   aryuZZ[356]=aryuZZ[283] + aryuZZ[354] - 11./3.*aryuZZ[67] - 20./3.
    + aryuZZ[308];
   aryuZZ[356]=2*aryuZZ[101]*aryuZZ[356];
   aryuZZ[358]= - 2*aryuZZ[58] - aryuZZ[56];
   aryuZZ[358]=56*aryuZZ[358] + 997./4.*aryuZZ[55];
   aryuZZ[359]= - 44*aryuZZ[53];
   aryuZZ[358]=aryuZZ[359] - aryuZZ[57] - aryuZZ[54] + 1./27.*
   aryuZZ[358] - 19./2.*aryuZZ[51];
   aryuZZ[360]=aryuZZ[19]*aryuZZ[38];
   aryuZZ[364]=aryuZZ[18]*aryuZZ[38];
   aryuZZ[366]=aryuZZ[360] + aryuZZ[364];
   aryuZZ[366]=aryuZZ[3]*aryuZZ[366];
   aryuZZ[144]=aryuZZ[10]*aryuZZ[144];
   aryuZZ[366]=18*aryuZZ[366] + 6*aryuZZ[144] - 5./2.*aryuZZ[36] + 6*
   aryuZZ[38] - 3./2.*aryuZZ[48] + 151./3. - 11*aryuZZ[46];
   aryuZZ[366]=aryuZZ[3]*aryuZZ[366];
   aryuZZ[367]= - aryuZZ[10]*aryuZZ[101]*aryuZZ[38];
   aryuZZ[368]=6*aryuZZ[367];
   aryuZZ[358]=aryuZZ[366] + aryuZZ[368] + 1./3.*aryuZZ[358] + 
   aryuZZ[356];
   aryuZZ[358]=MMt*aryuZZ[358];
   aryuZZ[369]=4213./3. + 511*aryuZZ[46];
   aryuZZ[369]=1./9.*aryuZZ[369] - 1./4.*aryuZZ[48];
   aryuZZ[369]=1./3.*aryuZZ[369] + aryuZZ[221];
   aryuZZ[371]=11./3.*aryuZZ[239];
   aryuZZ[369]=aryuZZ[371] + aryuZZ[361] + 301./54.*aryuZZ[13] + 1./4.*
   aryuZZ[369] + 29./9.*aryuZZ[12];
   aryuZZ[369]=1./2.*aryuZZ[369] + 61./81.*aryuZZ[35];
   aryuZZ[369]=aryuZZ[35]*aryuZZ[369];
   aryuZZ[373]=284./3. + aryuZZ[247];
   aryuZZ[374]= - 7*aryuZZ[7];
   aryuZZ[375]=487./3. + aryuZZ[374];
   aryuZZ[375]=aryuZZ[35]*aryuZZ[375];
   aryuZZ[375]=1./18.*aryuZZ[375] + 1./9.*aryuZZ[373] + aryuZZ[181];
   aryuZZ[375]=aryuZZ[1]*aryuZZ[375];
   aryuZZ[376]=1132./27. + aryuZZ[247];
   aryuZZ[377]=1823./27. + aryuZZ[374];
   aryuZZ[377]=aryuZZ[35]*aryuZZ[377];
   aryuZZ[378]= - 11./6.*aryuZZ[36];
   aryuZZ[377]=1./2.*aryuZZ[377] + aryuZZ[376] + aryuZZ[378];
   aryuZZ[377]=aryuZZ[39]*aryuZZ[377];
   aryuZZ[379]= - 7./4. + aryuZZ[184];
   aryuZZ[380]= - 80./9.*aryuZZ[35] - 89./9. + aryuZZ[229];
   aryuZZ[381]=aryuZZ[1]*aryuZZ[380];
   aryuZZ[380]=aryuZZ[39]*aryuZZ[380];
   aryuZZ[379]=11./9.*aryuZZ[380] + 1./3.*aryuZZ[379] + aryuZZ[381];
   aryuZZ[379]=aryuZZ[6]*aryuZZ[379];
   aryuZZ[380]=3*aryuZZ[69];
   aryuZZ[381]=aryuZZ[380] - 11./3.*aryuZZ[68];
   aryuZZ[382]=55./6. + aryuZZ[283];
   aryuZZ[382]=aryuZZ[17]*aryuZZ[382];
   aryuZZ[381]=aryuZZ[382] - 29./3.*aryuZZ[21] + 29./6.*aryuZZ[20] - 6*
   aryuZZ[50] + 11./3.*aryuZZ[41] + 2*aryuZZ[381] - 3*aryuZZ[40];
   aryuZZ[381]=aryuZZ[101]*aryuZZ[381];
   aryuZZ[382]= - 11./3.*aryuZZ[35] + aryuZZ[229] + 1./6. + aryuZZ[283]
   ;
   aryuZZ[382]=aryuZZ[10]*aryuZZ[382];
   aryuZZ[383]= - 1 - 2*aryuZZ[38];
   aryuZZ[383]=aryuZZ[101]*aryuZZ[383];
   aryuZZ[237]=3*aryuZZ[383] + aryuZZ[237];
   aryuZZ[237]=aryuZZ[18]*aryuZZ[237];
   aryuZZ[383]=5./8.*aryuZZ[44];
   aryuZZ[384]=11./2.*aryuZZ[15];
   aryuZZ[303]=8./3.*aryuZZ[303];
   aryuZZ[385]=10./3.*MMH*aryuZZ[53];
   aryuZZ[124]=aryuZZ[358] + aryuZZ[385] + aryuZZ[124] + aryuZZ[303] + 
   aryuZZ[379] + aryuZZ[237] + aryuZZ[382] + 1./9.*aryuZZ[377] + 1./3.*
   aryuZZ[375] + aryuZZ[369] + aryuZZ[381] + aryuZZ[334] + aryuZZ[242]
    + 967./108.*aryuZZ[13] + aryuZZ[228] + aryuZZ[158] + aryuZZ[384] + 
   1./3.*aryuZZ[324] + aryuZZ[383];
   aryuZZ[124]=MMt*aryuZZ[124];
   aryuZZ[324]=29 + 7*aryuZZ[46];
   aryuZZ[334]=aryuZZ[3]*aryuZZ[364];
   aryuZZ[324]=9*aryuZZ[334] + 3*aryuZZ[144] + aryuZZ[229] + 1./6.*
   aryuZZ[324] + aryuZZ[283];
   aryuZZ[324]=aryuZZ[3]*aryuZZ[324];
   aryuZZ[334]=aryuZZ[283] + aryuZZ[354] + 7./9.*aryuZZ[67] - 20./9. + 
   aryuZZ[308];
   aryuZZ[334]=aryuZZ[101]*aryuZZ[334];
   aryuZZ[219]=aryuZZ[219] + 17./3.*aryuZZ[54] + 437./108.*aryuZZ[55]
    + aryuZZ[52];
   aryuZZ[219]=1./2.*aryuZZ[219] + 14./3.*aryuZZ[53];
   aryuZZ[358]=3*aryuZZ[367];
   aryuZZ[219]=aryuZZ[324] + aryuZZ[358] + 1./3.*aryuZZ[219] + 
   aryuZZ[334];
   aryuZZ[219]=MMt*aryuZZ[219];
   aryuZZ[324]= - 35./18. + aryuZZ[283];
   aryuZZ[324]=aryuZZ[17]*aryuZZ[324];
   aryuZZ[334]= - 3*aryuZZ[50];
   aryuZZ[324]=1./2.*aryuZZ[324] - 47./18.*aryuZZ[21] + 47./36.*
   aryuZZ[20] + aryuZZ[334] - 7./18.*aryuZZ[41] + aryuZZ[110] + 
   aryuZZ[380] + 7./9.*aryuZZ[68];
   aryuZZ[324]=aryuZZ[101]*aryuZZ[324];
   aryuZZ[367]= - 5./2. - aryuZZ[13];
   aryuZZ[312]=1./4.*aryuZZ[367] + aryuZZ[312];
   aryuZZ[367]=7./9.*aryuZZ[35];
   aryuZZ[369]=aryuZZ[367] + 25./9. + aryuZZ[229];
   aryuZZ[375]=aryuZZ[1]*aryuZZ[369];
   aryuZZ[377]=aryuZZ[39]*aryuZZ[369];
   aryuZZ[312]=11./9.*aryuZZ[377] + 1./3.*aryuZZ[312] + aryuZZ[375];
   aryuZZ[312]=aryuZZ[6]*aryuZZ[312];
   aryuZZ[375]= - 3*aryuZZ[1] - 11./3.*aryuZZ[39];
   aryuZZ[377]=aryuZZ[6]*aryuZZ[375];
   aryuZZ[379]= - 3*aryuZZ[45];
   aryuZZ[377]=2*aryuZZ[377] - 9*aryuZZ[10] + 22./9.*aryuZZ[39] + 2*
   aryuZZ[1] + aryuZZ[240] + 17./2.*aryuZZ[12] + aryuZZ[379] + 4 + 7./3.
   *aryuZZ[46];
   aryuZZ[377]=aryuZZ[37]*aryuZZ[377];
   aryuZZ[112]=aryuZZ[112] + aryuZZ[50];
   aryuZZ[112]=41./12.*aryuZZ[327] - 41./12.*aryuZZ[17] + 7./2.*
   aryuZZ[21] + 7*aryuZZ[112] + 41./12.*aryuZZ[20];
   aryuZZ[327]=3*aryuZZ[36];
   aryuZZ[386]=7*aryuZZ[35] + aryuZZ[327] + 23./3. + aryuZZ[326];
   aryuZZ[386]=aryuZZ[18]*aryuZZ[386];
   aryuZZ[112]=aryuZZ[377] + 1./3.*aryuZZ[112] + 1./4.*aryuZZ[386];
   aryuZZ[112]=aryuZZ[3]*aryuZZ[112];
   aryuZZ[377]=1013 - 49*aryuZZ[46];
   aryuZZ[204]=aryuZZ[204] + 1./9.*aryuZZ[377] + aryuZZ[238];
   aryuZZ[204]=1./2.*aryuZZ[204] - 119./9.*aryuZZ[12];
   aryuZZ[204]=7./3.*aryuZZ[116] + 1./2.*aryuZZ[204] + 43./9.*
   aryuZZ[13];
   aryuZZ[204]=1./2.*aryuZZ[204] + 157./27.*aryuZZ[35];
   aryuZZ[204]=aryuZZ[35]*aryuZZ[204];
   aryuZZ[238]=aryuZZ[36] + 31./6. + aryuZZ[283];
   aryuZZ[238]=1./2.*aryuZZ[238] + aryuZZ[367];
   aryuZZ[238]=aryuZZ[10]*aryuZZ[238];
   aryuZZ[367]= - 1./2. - aryuZZ[38];
   aryuZZ[367]=3*aryuZZ[101]*aryuZZ[367];
   aryuZZ[164]=aryuZZ[367] + aryuZZ[164];
   aryuZZ[164]=aryuZZ[18]*aryuZZ[164];
   aryuZZ[377]= - 331./144.*aryuZZ[63] + 41./72. - aryuZZ[60];
   aryuZZ[386]=1./4.*aryuZZ[16];
   aryuZZ[387]= - 25./9. + aryuZZ[181];
   aryuZZ[388]=aryuZZ[387] - 7./9.*aryuZZ[35];
   aryuZZ[389]=aryuZZ[1]*aryuZZ[388];
   aryuZZ[388]=aryuZZ[39]*aryuZZ[388];
   aryuZZ[390]= - 3./2.*aryuZZ[38];
   aryuZZ[112]=aryuZZ[219] + 25./9.*aryuZZ[274] + aryuZZ[112] + 68./9.*
   aryuZZ[288] + aryuZZ[312] + aryuZZ[164] + aryuZZ[238] + 11./27.*
   aryuZZ[388] + 1./3.*aryuZZ[389] + 1./6.*aryuZZ[204] + aryuZZ[324] + 
   1./12.*aryuZZ[345] + 13./27.*aryuZZ[13] - 289./108.*aryuZZ[12] + 
   aryuZZ[390] - 7./12.*aryuZZ[15] - 25./48.*aryuZZ[44] - 25./192.*
   aryuZZ[48] - 301./432.*aryuZZ[46] + 227./108.*aryuZZ[67] + 31./18.*
   aryuZZ[64] - 7./36.*aryuZZ[59] + aryuZZ[386] + 5071./5184.*
   aryuZZ[49] + 257./1152.*aryuZZ[66] - 143./144.*aryuZZ[61] + 191./324.
   *aryuZZ[62] + 1./4.*aryuZZ[377] + 7./9.*aryuZZ[72];
   aryuZZ[112]=MMt*aryuZZ[112];
   aryuZZ[164]=5*aryuZZ[28] - 677./72. - 5*aryuZZ[30];
   aryuZZ[164]= - 19./18.*aryuZZ[15] + 25./24.*aryuZZ[44] - 125./36.*
   aryuZZ[67] - 523./108.*aryuZZ[64] + aryuZZ[157] + 37./8.*aryuZZ[49]
    + 5./54.*aryuZZ[14] + 1./9.*aryuZZ[164] - 37./8.*aryuZZ[61];
   aryuZZ[204]= - 17./9. - 25./16.*aryuZZ[44];
   aryuZZ[204]=1./2.*aryuZZ[204] + aryuZZ[347];
   aryuZZ[204]=aryuZZ[35]*aryuZZ[204];
   aryuZZ[219]=3./2. - 1./3.*aryuZZ[44];
   aryuZZ[219]=aryuZZ[99]*aryuZZ[219];
   aryuZZ[219]=1./4.*aryuZZ[219] + 1./9.*aryuZZ[225];
   aryuZZ[219]=aryuZZ[37]*aryuZZ[219];
   aryuZZ[238]= - 7 + aryuZZ[59];
   aryuZZ[143]=aryuZZ[143] + 1./6.*aryuZZ[67] + 1./6.*aryuZZ[238] - 
   aryuZZ[64];
   aryuZZ[143]=aryuZZ[99]*aryuZZ[143];
   aryuZZ[143]=aryuZZ[143] + 1./6.*aryuZZ[225];
   aryuZZ[143]=MMH*aryuZZ[143];
   aryuZZ[312]= - aryuZZ[68] + aryuZZ[330];
   aryuZZ[324]=1./12.*aryuZZ[21] - 13./6.*aryuZZ[50] + aryuZZ[312] - 1./
   12.*aryuZZ[41];
   aryuZZ[324]=aryuZZ[99]*aryuZZ[324];
   aryuZZ[330]=49./6. - 187*aryuZZ[35];
   aryuZZ[330]=aryuZZ[10]*aryuZZ[330];
   aryuZZ[345]=aryuZZ[218] - 11./4.*aryuZZ[10];
   aryuZZ[345]=aryuZZ[6]*aryuZZ[345];
   aryuZZ[347]= - 11./81.*aryuZZ[11];
   aryuZZ[143]=25./96.*aryuZZ[143] + 25./4.*aryuZZ[219] + 5./54.*
   aryuZZ[345] + 25./576.*aryuZZ[351] + 1./216.*aryuZZ[330] + 1./6.*
   aryuZZ[204] + 25./48.*aryuZZ[277] + 25./48.*aryuZZ[324] + 1./8.*
   aryuZZ[164] + aryuZZ[347];
   aryuZZ[143]=MMH*aryuZZ[143];
   aryuZZ[164]= - 3229./6. - 35*aryuZZ[46];
   aryuZZ[164]=aryuZZ[205] + 1./9.*aryuZZ[164] + 5./4.*aryuZZ[48];
   aryuZZ[164]=563./108.*aryuZZ[35] + 25./9.*aryuZZ[13] + 5./4.*
   aryuZZ[164] + aryuZZ[245];
   aryuZZ[204]=107*aryuZZ[34] - 2825./6.*aryuZZ[99];
   aryuZZ[204]=aryuZZ[18]*aryuZZ[204];
   aryuZZ[205]=275./9.*aryuZZ[39] + 1./4. + 25*aryuZZ[1];
   aryuZZ[205]=aryuZZ[6]*aryuZZ[205];
   aryuZZ[219]=107./4.*aryuZZ[34] - 350./3.*aryuZZ[99];
   aryuZZ[219]=aryuZZ[37]*aryuZZ[219];
   aryuZZ[164]=1./3.*aryuZZ[219] + 1./3.*aryuZZ[205] + 1./12.*
   aryuZZ[204] + aryuZZ[284] - 275./81.*aryuZZ[39] + 1./4.*aryuZZ[164]
    - 25./9.*aryuZZ[1];
   aryuZZ[164]=aryuZZ[37]*aryuZZ[164];
   aryuZZ[192]=1741./72.*aryuZZ[50] + 1./3.*aryuZZ[41] - 25./4.*
   aryuZZ[40] + 259./18.*aryuZZ[68] - 215./27.*aryuZZ[8] + 25./8.*
   aryuZZ[192] + 59./9.*aryuZZ[69];
   aryuZZ[204]= - 635./6.*aryuZZ[19] + aryuZZ[281];
   aryuZZ[204]=aryuZZ[35]*aryuZZ[204];
   aryuZZ[205]=571 - 233*aryuZZ[35];
   aryuZZ[205]=aryuZZ[18]*aryuZZ[205];
   aryuZZ[219]=5*aryuZZ[17];
   aryuZZ[245]= - 167./3.*aryuZZ[19] + aryuZZ[219];
   aryuZZ[245]=1./8.*aryuZZ[245] - 25./9.*aryuZZ[18];
   aryuZZ[245]=aryuZZ[6]*aryuZZ[245];
   aryuZZ[164]=aryuZZ[164] + 1./3.*aryuZZ[245] + 1./288.*aryuZZ[205] + 
   1./12.*aryuZZ[204] + 2./3.*aryuZZ[17] + 7661./1728.*aryuZZ[19] - 311.
   /216.*aryuZZ[21] - 73./24.*aryuZZ[20] + 1./8.*aryuZZ[192] - 5./3.*
   aryuZZ[22];
   aryuZZ[192]=11./3. - 17./2.*aryuZZ[35];
   aryuZZ[192]=aryuZZ[18]*aryuZZ[192];
   aryuZZ[204]=aryuZZ[6]*aryuZZ[357];
   aryuZZ[205]= - aryuZZ[37]*aryuZZ[12];
   aryuZZ[192]=17*aryuZZ[205] + 5./2.*aryuZZ[204] + aryuZZ[192] - 11./3.
   *aryuZZ[19] + 17./2.*aryuZZ[353];
   aryuZZ[204]= - 17*aryuZZ[12];
   aryuZZ[205]=aryuZZ[204] + 7./2.*aryuZZ[268];
   aryuZZ[245]=1 + aryuZZ[12];
   aryuZZ[245]=aryuZZ[3]*aryuZZ[37]*aryuZZ[245];
   aryuZZ[268]=MMt*aryuZZ[3]*aryuZZ[36];
   aryuZZ[205]=aryuZZ[268] + 1./54.*aryuZZ[205] + aryuZZ[245];
   aryuZZ[205]=MMt*aryuZZ[205];
   aryuZZ[192]=1./54.*aryuZZ[192] + aryuZZ[205];
   aryuZZ[192]=aryuZZ[272]*aryuZZ[192];
   aryuZZ[205]=7*aryuZZ[18] - 17./6.*aryuZZ[37];
   aryuZZ[205]=aryuZZ[3]*aryuZZ[37]*aryuZZ[205];
   aryuZZ[112]=1./2.*aryuZZ[192] + aryuZZ[112] + aryuZZ[143] + 1./3.*
   aryuZZ[164] + aryuZZ[205];
   aryuZZ[112]=aryuZZ[272]*aryuZZ[112];
   aryuZZ[143]=aryuZZ[28] - 1085./216. - aryuZZ[30];
   aryuZZ[143]= - 5./6.*aryuZZ[15] - 5./8.*aryuZZ[44] + 25./12.*
   aryuZZ[67] - 1./36.*aryuZZ[64] + aryuZZ[157] - 3./8.*aryuZZ[49] + 1./
   18.*aryuZZ[14] + 1./3.*aryuZZ[143] + 3./8.*aryuZZ[61];
   aryuZZ[164]=23./27. + aryuZZ[383];
   aryuZZ[164]=1./2.*aryuZZ[164] + aryuZZ[314];
   aryuZZ[164]=aryuZZ[35]*aryuZZ[164];
   aryuZZ[192]= - 9./2. + aryuZZ[44];
   aryuZZ[192]=aryuZZ[99]*aryuZZ[192];
   aryuZZ[168]=1./4.*aryuZZ[192] + aryuZZ[168];
   aryuZZ[168]=aryuZZ[37]*aryuZZ[168];
   aryuZZ[192]=7 - aryuZZ[59];
   aryuZZ[192]=1./6.*aryuZZ[15] - 1./6.*aryuZZ[67] + 1./6.*aryuZZ[192]
    + aryuZZ[64];
   aryuZZ[192]=aryuZZ[99]*aryuZZ[192];
   aryuZZ[175]=aryuZZ[192] + 1./6.*aryuZZ[175];
   aryuZZ[175]=MMH*aryuZZ[175];
   aryuZZ[192]= - 1./12.*aryuZZ[21] + 13./6.*aryuZZ[50] + 1./12.*
   aryuZZ[41] + aryuZZ[68] + aryuZZ[252];
   aryuZZ[192]=aryuZZ[99]*aryuZZ[192];
   aryuZZ[205]=aryuZZ[355] + 7./12.*aryuZZ[35];
   aryuZZ[205]=aryuZZ[10]*aryuZZ[205];
   aryuZZ[245]= - 1./2. + 5*aryuZZ[11];
   aryuZZ[245]=1./3.*aryuZZ[245] - 7./2.*aryuZZ[10];
   aryuZZ[245]=aryuZZ[6]*aryuZZ[245];
   aryuZZ[143]=5./16.*aryuZZ[175] + 5./2.*aryuZZ[168] + 1./18.*
   aryuZZ[245] + 5./96.*aryuZZ[349] + 1./3.*aryuZZ[205] + 1./2.*
   aryuZZ[164] + 5./8.*aryuZZ[213] + aryuZZ[279] + 5./8.*aryuZZ[192] + 
   1./4.*aryuZZ[143] + aryuZZ[347];
   aryuZZ[143]=MMH*aryuZZ[143];
   aryuZZ[164]=28877./6. + 1105*aryuZZ[46];
   aryuZZ[164]=1./9.*aryuZZ[164] + aryuZZ[226];
   aryuZZ[164]=1./3.*aryuZZ[164] + aryuZZ[221];
   aryuZZ[168]=1./2.*aryuZZ[98] - 11*aryuZZ[34];
   aryuZZ[168]=aryuZZ[19]*aryuZZ[168];
   aryuZZ[175]=aryuZZ[1]*aryuZZ[373];
   aryuZZ[192]=aryuZZ[39]*aryuZZ[376];
   aryuZZ[203]=aryuZZ[203] + 6695./6.*aryuZZ[99];
   aryuZZ[203]=aryuZZ[18]*aryuZZ[203];
   aryuZZ[205]= - 979./9.*aryuZZ[39] + 5./2. - 89*aryuZZ[1];
   aryuZZ[205]=aryuZZ[6]*aryuZZ[205];
   aryuZZ[213]= - 173./2.*aryuZZ[34] + 1420./3.*aryuZZ[99];
   aryuZZ[213]=aryuZZ[37]*aryuZZ[213];
   aryuZZ[226]= - 25./6.*aryuZZ[10];
   aryuZZ[164]=1./9.*aryuZZ[213] + 1./9.*aryuZZ[205] + 1./18.*
   aryuZZ[203] + aryuZZ[226] + 1./9.*aryuZZ[192] + 1./27.*aryuZZ[175]
    + 851./648.*aryuZZ[35] + 1./6.*aryuZZ[168] + aryuZZ[242] + 481./108.
   *aryuZZ[13] + 1./8.*aryuZZ[164] + aryuZZ[228];
   aryuZZ[164]=aryuZZ[37]*aryuZZ[164];
   aryuZZ[168]= - 11*aryuZZ[18];
   aryuZZ[175]=13*aryuZZ[19];
   aryuZZ[192]=aryuZZ[175] + aryuZZ[168];
   aryuZZ[203]= - 15*aryuZZ[37];
   aryuZZ[192]=2./3.*aryuZZ[192] + aryuZZ[203];
   aryuZZ[192]=aryuZZ[3]*aryuZZ[37]*aryuZZ[192];
   aryuZZ[205]= - 47./8.*aryuZZ[68] + 83./36.*aryuZZ[8] + 139./12.*
   aryuZZ[69] - 425./64.*aryuZZ[42] - 68./3.*aryuZZ[71] - 151./32.*
   aryuZZ[70];
   aryuZZ[213]= - 169./48. - aryuZZ[36];
   aryuZZ[213]=aryuZZ[19]*aryuZZ[213];
   aryuZZ[228]= - 11./4.*aryuZZ[17];
   aryuZZ[245]=37*aryuZZ[19] + aryuZZ[228];
   aryuZZ[245]=aryuZZ[35]*aryuZZ[245];
   aryuZZ[252]=2131./72.*aryuZZ[35] + 7049./216. - aryuZZ[36];
   aryuZZ[252]=aryuZZ[18]*aryuZZ[252];
   aryuZZ[268]=aryuZZ[254] + 281./18.*aryuZZ[18];
   aryuZZ[268]=aryuZZ[6]*aryuZZ[268];
   aryuZZ[281]=5./16.*aryuZZ[40];
   aryuZZ[284]=7./12.*aryuZZ[41];
   aryuZZ[314]=85./54. - aryuZZ[36];
   aryuZZ[314]=1./2.*aryuZZ[17]*aryuZZ[314];
   aryuZZ[324]=1./4.*aryuZZ[20];
   aryuZZ[112]=aryuZZ[112] + aryuZZ[124] + aryuZZ[143] + aryuZZ[192] + 
   aryuZZ[164] + 1./18.*aryuZZ[268] + 1./6.*aryuZZ[252] + 1./9.*
   aryuZZ[245] + aryuZZ[314] + 1./6.*aryuZZ[213] - 469./324.*aryuZZ[21]
    + aryuZZ[324] + aryuZZ[22] - 2099./2592.*aryuZZ[50] + aryuZZ[284]
    + 1./9.*aryuZZ[205] + aryuZZ[281];
   aryuZZ[112]=aryuZZ[272]*aryuZZ[112];
   aryuZZ[105]=aryuZZ[222] + aryuZZ[178] + aryuZZ[170] + aryuZZ[122] + 
   aryuZZ[316] - 113./2.*aryuZZ[13] + aryuZZ[121] + aryuZZ[117] + 
   aryuZZ[294] + aryuZZ[246] - 191./2. + aryuZZ[105];
   aryuZZ[105]=aryuZZ[37]*aryuZZ[105];
   aryuZZ[117]=aryuZZ[275] - 121./3. + aryuZZ[326];
   aryuZZ[117]=aryuZZ[19]*aryuZZ[117];
   aryuZZ[121]= - 121*aryuZZ[19] + aryuZZ[329];
   aryuZZ[121]=aryuZZ[35]*aryuZZ[121];
   aryuZZ[105]=aryuZZ[105] + aryuZZ[344] + 1./6.*aryuZZ[121] + 
   aryuZZ[339] + aryuZZ[338] + 1./2.*aryuZZ[117];
   aryuZZ[105]=aryuZZ[3]*aryuZZ[105];
   aryuZZ[117]=113./8.*aryuZZ[63] + 30167./108. - 43*aryuZZ[60];
   aryuZZ[117]=aryuZZ[352] + aryuZZ[350] + aryuZZ[348] + 43./4.*
   aryuZZ[16] + 89./96.*aryuZZ[49] - 5671./64.*aryuZZ[66] + aryuZZ[346]
    - 19./32.*aryuZZ[65] - 11./2.*aryuZZ[62] + 1./4.*aryuZZ[117] + 
   aryuZZ[328];
   aryuZZ[121]=2*aryuZZ[58] + aryuZZ[56];
   aryuZZ[122]=aryuZZ[359] - aryuZZ[57] - aryuZZ[54] - 35./2.*
   aryuZZ[51] + 8./3.*aryuZZ[121] - 1./4.*aryuZZ[55];
   aryuZZ[122]=aryuZZ[366] + aryuZZ[368] + 1./3.*aryuZZ[122] + 
   aryuZZ[356];
   aryuZZ[122]=MMt*aryuZZ[122];
   aryuZZ[124]=aryuZZ[221] - 43./4.*aryuZZ[48] + 10453./81. - 19*
   aryuZZ[46];
   aryuZZ[124]=aryuZZ[371] + aryuZZ[361] + 581./6.*aryuZZ[13] + 1./4.*
   aryuZZ[124] + 31./3.*aryuZZ[12];
   aryuZZ[124]=1./2.*aryuZZ[124] + aryuZZ[282];
   aryuZZ[124]=aryuZZ[35]*aryuZZ[124];
   aryuZZ[164]=92./3. + aryuZZ[247];
   aryuZZ[170]=103./3. + aryuZZ[374];
   aryuZZ[170]=aryuZZ[35]*aryuZZ[170];
   aryuZZ[178]=1./18.*aryuZZ[170] + 1./9.*aryuZZ[164] + aryuZZ[181];
   aryuZZ[178]=aryuZZ[1]*aryuZZ[178];
   aryuZZ[192]=76./3. + aryuZZ[247];
   aryuZZ[170]=1./2.*aryuZZ[170] + aryuZZ[192] + aryuZZ[378];
   aryuZZ[170]=aryuZZ[39]*aryuZZ[170];
   aryuZZ[103]=aryuZZ[154] - 7./4. + aryuZZ[103];
   aryuZZ[154]= - 59 + 11./2.*aryuZZ[36];
   aryuZZ[154]=1./3.*aryuZZ[154] + aryuZZ[115];
   aryuZZ[154]=aryuZZ[39]*aryuZZ[154];
   aryuZZ[205]= - 16./9.*aryuZZ[35] - 25./9. + aryuZZ[229];
   aryuZZ[205]=aryuZZ[1]*aryuZZ[205];
   aryuZZ[103]=1./3.*aryuZZ[154] + 1./3.*aryuZZ[103] + aryuZZ[205];
   aryuZZ[103]=aryuZZ[6]*aryuZZ[103];
   aryuZZ[154]=aryuZZ[36]*aryuZZ[363];
   aryuZZ[103]=aryuZZ[122] + aryuZZ[385] + aryuZZ[105] + aryuZZ[303] + 
   aryuZZ[103] + aryuZZ[237] + aryuZZ[382] + 1./9.*aryuZZ[170] + 1./3.*
   aryuZZ[178] + aryuZZ[124] + aryuZZ[381] + 23./4.*aryuZZ[154] + 
   aryuZZ[242] + 863./12.*aryuZZ[13] + aryuZZ[241] + aryuZZ[158] + 
   aryuZZ[384] + aryuZZ[383] - 9./32.*aryuZZ[48] + 1./3.*aryuZZ[117] + 
   25./8.*aryuZZ[46];
   aryuZZ[103]=MMt*aryuZZ[103];
   aryuZZ[105]= - 3*aryuZZ[46];
   aryuZZ[117]=1 + aryuZZ[150];
   aryuZZ[117]=aryuZZ[39]*aryuZZ[117];
   aryuZZ[122]= - 3*aryuZZ[10];
   aryuZZ[124]= - aryuZZ[1] - 3*aryuZZ[39];
   aryuZZ[154]=aryuZZ[6]*aryuZZ[124];
   aryuZZ[117]=4*aryuZZ[154] + aryuZZ[122] + 2*aryuZZ[117] + 2*
   aryuZZ[320] + aryuZZ[316] + 132*aryuZZ[13] - 165./2.*aryuZZ[12] + 
   aryuZZ[379] + aryuZZ[294] + aryuZZ[246] + 379./6. + aryuZZ[105];
   aryuZZ[117]=aryuZZ[37]*aryuZZ[117];
   aryuZZ[154]= - aryuZZ[22] + aryuZZ[193] + aryuZZ[42] + aryuZZ[232];
   aryuZZ[158]=13 + aryuZZ[326];
   aryuZZ[170]=aryuZZ[158] + aryuZZ[275];
   aryuZZ[170]=aryuZZ[19]*aryuZZ[170];
   aryuZZ[175]=aryuZZ[175] + aryuZZ[119];
   aryuZZ[175]=aryuZZ[35]*aryuZZ[175];
   aryuZZ[178]=5*aryuZZ[35];
   aryuZZ[158]=1./2.*aryuZZ[158] + aryuZZ[178];
   aryuZZ[158]=aryuZZ[18]*aryuZZ[158];
   aryuZZ[117]=aryuZZ[117] + 1./2.*aryuZZ[158] + 1./2.*aryuZZ[175] - 1./
   4.*aryuZZ[17] + 1./2.*aryuZZ[170] - 3./2.*aryuZZ[21] + 3*aryuZZ[154]
    + aryuZZ[324];
   aryuZZ[117]=aryuZZ[3]*aryuZZ[117];
   aryuZZ[154]=5./2. + aryuZZ[283];
   aryuZZ[154]=aryuZZ[17]*aryuZZ[154];
   aryuZZ[110]=1./2.*aryuZZ[154] - 7./2.*aryuZZ[21] + 7./4.*aryuZZ[20]
    + aryuZZ[334] + aryuZZ[232] + aryuZZ[110] + aryuZZ[380] - 
   aryuZZ[68];
   aryuZZ[110]=aryuZZ[101]*aryuZZ[110];
   aryuZZ[154]=aryuZZ[283] + aryuZZ[354] - aryuZZ[67] - 4 + aryuZZ[308]
   ;
   aryuZZ[154]=aryuZZ[101]*aryuZZ[154];
   aryuZZ[158]=2*aryuZZ[360] + aryuZZ[364];
   aryuZZ[158]=aryuZZ[3]*aryuZZ[158];
   aryuZZ[170]= - aryuZZ[48] + 9 - aryuZZ[46];
   aryuZZ[144]=3*aryuZZ[158] + aryuZZ[144] - aryuZZ[36] + 1./2.*
   aryuZZ[170] + aryuZZ[38];
   aryuZZ[144]=aryuZZ[3]*aryuZZ[144];
   aryuZZ[158]= - 3*aryuZZ[57] + aryuZZ[54] - 3*aryuZZ[51] + 5./4.*
   aryuZZ[55] + aryuZZ[159];
   aryuZZ[144]=3*aryuZZ[144] + aryuZZ[358] + aryuZZ[154] + 1./2.*
   aryuZZ[158] + aryuZZ[223];
   aryuZZ[144]=MMt*aryuZZ[144];
   aryuZZ[105]=aryuZZ[160] + aryuZZ[182] - 1799./27. + aryuZZ[105];
   aryuZZ[105]=1./2.*aryuZZ[105] - 55*aryuZZ[12];
   aryuZZ[105]=aryuZZ[35] + 1./2.*aryuZZ[239] + aryuZZ[11] + 1./4.*
   aryuZZ[105] - 19*aryuZZ[13];
   aryuZZ[105]=aryuZZ[35]*aryuZZ[105];
   aryuZZ[154]=1./2. + 7*aryuZZ[13];
   aryuZZ[158]=14 + aryuZZ[178];
   aryuZZ[159]=aryuZZ[1]*aryuZZ[158];
   aryuZZ[158]=aryuZZ[39]*aryuZZ[158];
   aryuZZ[154]=1./3.*aryuZZ[158] + 1./9.*aryuZZ[159] + 1./4.*
   aryuZZ[154] - aryuZZ[35];
   aryuZZ[154]=aryuZZ[6]*aryuZZ[154];
   aryuZZ[158]= - 13./9. + aryuZZ[7];
   aryuZZ[158]=aryuZZ[35]*aryuZZ[158];
   aryuZZ[158]=aryuZZ[183] + 1./2.*aryuZZ[158];
   aryuZZ[159]=aryuZZ[1]*aryuZZ[158];
   aryuZZ[158]=aryuZZ[39]*aryuZZ[158];
   aryuZZ[170]=35./6. + aryuZZ[283];
   aryuZZ[170]=1./2.*aryuZZ[170] + aryuZZ[194];
   aryuZZ[170]=aryuZZ[10]*aryuZZ[170];
   aryuZZ[169]=aryuZZ[367] + aryuZZ[169];
   aryuZZ[169]=aryuZZ[18]*aryuZZ[169];
   aryuZZ[175]=55./4.*aryuZZ[12];
   aryuZZ[178]= - 3 + aryuZZ[210];
   aryuZZ[178]=aryuZZ[36]*aryuZZ[178];
   aryuZZ[105]=aryuZZ[144] + aryuZZ[274] + aryuZZ[117] + 4*aryuZZ[288]
    + aryuZZ[154] + aryuZZ[169] + aryuZZ[170] + aryuZZ[158] + 1./3.*
   aryuZZ[159] + 1./2.*aryuZZ[105] + aryuZZ[110] + 1./2.*aryuZZ[178] - 
   aryuZZ[11] - 203./4.*aryuZZ[13] + aryuZZ[175] + aryuZZ[390] + 3./4.*
   aryuZZ[15] - 3./16.*aryuZZ[44] + 33./64.*aryuZZ[48] + 9./16.*
   aryuZZ[46] + 7./12.*aryuZZ[67] + aryuZZ[244] + 1./4.*aryuZZ[59] - 12
   *aryuZZ[16] + 3./64.*aryuZZ[49] + 3533./128.*aryuZZ[66] - 7./16.*
   aryuZZ[61] - 25./32.*aryuZZ[65] - 1./4.*aryuZZ[62] - aryuZZ[72] + 41.
   /64.*aryuZZ[63] - 17197./864. + 12*aryuZZ[60];
   aryuZZ[105]=MMt*aryuZZ[105];
   aryuZZ[110]=33*aryuZZ[12];
   aryuZZ[117]= - 33*aryuZZ[13];
   aryuZZ[144]=aryuZZ[117] - 1 + aryuZZ[110];
   aryuZZ[144]=1./2.*aryuZZ[144] + aryuZZ[296];
   aryuZZ[154]= - aryuZZ[1]*aryuZZ[7];
   aryuZZ[158]= - aryuZZ[39]*aryuZZ[7];
   aryuZZ[159]=aryuZZ[1] + 3*aryuZZ[39];
   aryuZZ[169]=aryuZZ[6]*aryuZZ[159];
   aryuZZ[144]=2*aryuZZ[169] + 6*aryuZZ[10] + 6*aryuZZ[158] + 3*
   aryuZZ[144] + 2*aryuZZ[154];
   aryuZZ[144]=aryuZZ[37]*aryuZZ[144];
   aryuZZ[145]=aryuZZ[18]*aryuZZ[145];
   aryuZZ[145]=1./2.*aryuZZ[145] + aryuZZ[202] + aryuZZ[353];
   aryuZZ[144]=3./2.*aryuZZ[145] + aryuZZ[144];
   aryuZZ[144]=aryuZZ[3]*aryuZZ[144];
   aryuZZ[145]=33./2.*aryuZZ[12];
   aryuZZ[169]=aryuZZ[117] - 59./9. + aryuZZ[145];
   aryuZZ[169]=1./2.*aryuZZ[169] - aryuZZ[11];
   aryuZZ[169]=aryuZZ[35]*aryuZZ[169];
   aryuZZ[170]=aryuZZ[35]*aryuZZ[263];
   aryuZZ[170]=1./2.*aryuZZ[170] + aryuZZ[7] + aryuZZ[279];
   aryuZZ[178]=aryuZZ[1]*aryuZZ[170];
   aryuZZ[170]=aryuZZ[39]*aryuZZ[170];
   aryuZZ[182]=aryuZZ[1]*aryuZZ[276];
   aryuZZ[193]=aryuZZ[39]*aryuZZ[276];
   aryuZZ[182]=1./3.*aryuZZ[182] + aryuZZ[193];
   aryuZZ[182]=aryuZZ[6]*aryuZZ[182];
   aryuZZ[186]=33./4.*aryuZZ[186] + aryuZZ[11];
   aryuZZ[193]= - MMt*aryuZZ[3]*aryuZZ[36];
   aryuZZ[144]=3./2.*aryuZZ[193] + aryuZZ[144] + aryuZZ[182] + 
   aryuZZ[280] + aryuZZ[170] + 1./3.*aryuZZ[178] + 1./2.*aryuZZ[169] + 
   aryuZZ[186] + 1./4.*aryuZZ[234];
   aryuZZ[144]=MMt*aryuZZ[144];
   aryuZZ[169]= - 17./18. + aryuZZ[196];
   aryuZZ[169]=aryuZZ[19]*aryuZZ[169];
   aryuZZ[170]=47./6.*aryuZZ[19] - aryuZZ[17];
   aryuZZ[178]=aryuZZ[35]*aryuZZ[170];
   aryuZZ[182]= - 19./4.*aryuZZ[35] + 17./6. + aryuZZ[36];
   aryuZZ[182]=aryuZZ[18]*aryuZZ[182];
   aryuZZ[193]=aryuZZ[170] - 19./6.*aryuZZ[18];
   aryuZZ[193]=aryuZZ[6]*aryuZZ[193];
   aryuZZ[169]=1./2.*aryuZZ[193] + 1./3.*aryuZZ[182] + 1./2.*
   aryuZZ[178] + aryuZZ[169] + aryuZZ[264];
   aryuZZ[178]=1./2. - aryuZZ[11];
   aryuZZ[182]=aryuZZ[35]*aryuZZ[178];
   aryuZZ[193]=aryuZZ[10]*aryuZZ[249];
   aryuZZ[196]=1./2.*aryuZZ[178] + aryuZZ[10];
   aryuZZ[196]=aryuZZ[6]*aryuZZ[196];
   aryuZZ[181]=aryuZZ[196] + aryuZZ[193] + 1./2.*aryuZZ[182] + 
   aryuZZ[216] + aryuZZ[181];
   aryuZZ[181]=MMH*aryuZZ[181];
   aryuZZ[182]=aryuZZ[1]*aryuZZ[7];
   aryuZZ[193]=aryuZZ[39]*aryuZZ[7];
   aryuZZ[140]=aryuZZ[140] - aryuZZ[10] + aryuZZ[193] + aryuZZ[186] + 1.
   /3.*aryuZZ[182];
   aryuZZ[140]=aryuZZ[37]*aryuZZ[140];
   aryuZZ[140]=aryuZZ[144] + 1./3.*aryuZZ[181] + 1./2.*aryuZZ[169] + 
   aryuZZ[140];
   aryuZZ[140]=aryuZZ[133]*aryuZZ[140];
   aryuZZ[144]=aryuZZ[160] + 15./4.*aryuZZ[48] - 53359./54. + 
   aryuZZ[129];
   aryuZZ[144]=1./4.*aryuZZ[144] + aryuZZ[180];
   aryuZZ[160]=3./2.*aryuZZ[98] - aryuZZ[34];
   aryuZZ[160]=aryuZZ[19]*aryuZZ[160];
   aryuZZ[169]=aryuZZ[1]*aryuZZ[183];
   aryuZZ[180]=aryuZZ[39]*aryuZZ[183];
   aryuZZ[181]=3*aryuZZ[34] - 11./2.*aryuZZ[99];
   aryuZZ[181]=aryuZZ[18]*aryuZZ[181];
   aryuZZ[183]=14./3.*aryuZZ[39] + 5./4. + 14./9.*aryuZZ[1];
   aryuZZ[183]=aryuZZ[6]*aryuZZ[183];
   aryuZZ[196]=3./4.*aryuZZ[34] - 2*aryuZZ[99];
   aryuZZ[196]=aryuZZ[37]*aryuZZ[196];
   aryuZZ[144]=aryuZZ[196] + aryuZZ[183] + 1./4.*aryuZZ[181] + 41./12.*
   aryuZZ[10] + aryuZZ[180] + 1./3.*aryuZZ[169] + 3./16.*aryuZZ[35] + 1.
   /2.*aryuZZ[160] - aryuZZ[11] + 1./4.*aryuZZ[144] - 26*aryuZZ[13];
   aryuZZ[144]=aryuZZ[37]*aryuZZ[144];
   aryuZZ[142]=aryuZZ[142] - 13./2.*aryuZZ[50] + 3*aryuZZ[312] + 
   aryuZZ[190];
   aryuZZ[142]=aryuZZ[99]*aryuZZ[142];
   aryuZZ[160]=1./6.*aryuZZ[14];
   aryuZZ[142]=1./2.*aryuZZ[142] - 3./2.*aryuZZ[15] + 3./8.*aryuZZ[44]
    - 5./4.*aryuZZ[67] - 35./12.*aryuZZ[64] + aryuZZ[157] + 21./8.*
   aryuZZ[49] + aryuZZ[160] - 21./8.*aryuZZ[61] + aryuZZ[28] - 109./72.
    - aryuZZ[30];
   aryuZZ[142]=3./8.*aryuZZ[277] + 1./4.*aryuZZ[142] + aryuZZ[188];
   aryuZZ[155]= - 1./2.*aryuZZ[15] + aryuZZ[167] + 1./2.*aryuZZ[238] + 
   aryuZZ[155];
   aryuZZ[155]=aryuZZ[99]*aryuZZ[155];
   aryuZZ[155]=aryuZZ[155] + 1./2.*aryuZZ[225];
   aryuZZ[155]=MMH*aryuZZ[155];
   aryuZZ[157]=9./2. - aryuZZ[44];
   aryuZZ[157]=aryuZZ[99]*aryuZZ[157];
   aryuZZ[157]=3./4.*aryuZZ[157] + aryuZZ[225];
   aryuZZ[157]=aryuZZ[37]*aryuZZ[157];
   aryuZZ[167]= - 2./9. + aryuZZ[214];
   aryuZZ[167]=aryuZZ[35]*aryuZZ[167];
   aryuZZ[169]=aryuZZ[207] - 41./24.*aryuZZ[35];
   aryuZZ[169]=aryuZZ[10]*aryuZZ[169];
   aryuZZ[180]= - 1 - 25./8.*aryuZZ[10];
   aryuZZ[180]=aryuZZ[6]*aryuZZ[180];
   aryuZZ[142]=1./32.*aryuZZ[155] + 1./4.*aryuZZ[157] + 1./9.*
   aryuZZ[180] + 1./64.*aryuZZ[351] + 1./3.*aryuZZ[169] + 1./2.*
   aryuZZ[142] + aryuZZ[167];
   aryuZZ[142]=MMH*aryuZZ[142];
   aryuZZ[155]= - 109*aryuZZ[19] + aryuZZ[219];
   aryuZZ[155]=1./2.*aryuZZ[155] + aryuZZ[200];
   aryuZZ[155]=aryuZZ[6]*aryuZZ[155];
   aryuZZ[157]=943*aryuZZ[70] - 399./2.*aryuZZ[42];
   aryuZZ[157]=1011./8.*aryuZZ[50] + aryuZZ[41] - 3./4.*aryuZZ[40] + 11.
   /6.*aryuZZ[68] - 11*aryuZZ[8] + 1./8.*aryuZZ[157] + aryuZZ[69];
   aryuZZ[167]= - 607./32. + aryuZZ[189];
   aryuZZ[167]=aryuZZ[19]*aryuZZ[167];
   aryuZZ[169]=1./2. - aryuZZ[36];
   aryuZZ[169]=aryuZZ[17]*aryuZZ[169];
   aryuZZ[180]= - 101./4.*aryuZZ[19] + aryuZZ[17];
   aryuZZ[180]=aryuZZ[35]*aryuZZ[180];
   aryuZZ[181]=517./48.*aryuZZ[35] + 307./16. - aryuZZ[36];
   aryuZZ[181]=aryuZZ[18]*aryuZZ[181];
   aryuZZ[183]=aryuZZ[132] + aryuZZ[18];
   aryuZZ[183]=7*aryuZZ[183] - 3./2.*aryuZZ[37];
   aryuZZ[183]=aryuZZ[3]*aryuZZ[37]*aryuZZ[183];
   aryuZZ[105]=aryuZZ[140] + aryuZZ[105] + aryuZZ[142] + aryuZZ[183] + 
   aryuZZ[144] + 1./12.*aryuZZ[155] + 1./6.*aryuZZ[181] + 1./6.*
   aryuZZ[180] + 1./2.*aryuZZ[169] + 1./2.*aryuZZ[167] - 11./8.*
   aryuZZ[21] - 5./8.*aryuZZ[20] + 1./8.*aryuZZ[157] - 10*aryuZZ[22];
   aryuZZ[105]=aryuZZ[133]*aryuZZ[105];
   aryuZZ[129]=aryuZZ[221] - 31./4.*aryuZZ[48] + 140957./162. + 
   aryuZZ[129];
   aryuZZ[140]= - 1./2.*aryuZZ[98] + aryuZZ[253];
   aryuZZ[140]=aryuZZ[19]*aryuZZ[140];
   aryuZZ[142]=aryuZZ[1]*aryuZZ[164];
   aryuZZ[144]=aryuZZ[39]*aryuZZ[192];
   aryuZZ[155]= - aryuZZ[34] + 35./6.*aryuZZ[99];
   aryuZZ[155]=aryuZZ[18]*aryuZZ[155];
   aryuZZ[157]= - 59./3.*aryuZZ[39] - 1./2. - 25./3.*aryuZZ[1];
   aryuZZ[157]=aryuZZ[6]*aryuZZ[157];
   aryuZZ[164]= - 5./2.*aryuZZ[34] + 44./3.*aryuZZ[99];
   aryuZZ[164]=aryuZZ[37]*aryuZZ[164];
   aryuZZ[129]=aryuZZ[164] + 1./3.*aryuZZ[157] + 5./2.*aryuZZ[155] + 
   aryuZZ[226] + 1./9.*aryuZZ[144] + 1./27.*aryuZZ[142] - 133./72.*
   aryuZZ[35] + 1./2.*aryuZZ[140] + aryuZZ[242] + 665./12.*aryuZZ[13]
    + 1./8.*aryuZZ[129] + aryuZZ[241];
   aryuZZ[129]=aryuZZ[37]*aryuZZ[129];
   aryuZZ[140]= - 19*aryuZZ[19] + aryuZZ[168];
   aryuZZ[140]=2./3.*aryuZZ[140] + aryuZZ[203];
   aryuZZ[140]=aryuZZ[3]*aryuZZ[37]*aryuZZ[140];
   aryuZZ[142]=69551./432. + 35*aryuZZ[36];
   aryuZZ[142]=aryuZZ[19]*aryuZZ[142];
   aryuZZ[144]=367./3.*aryuZZ[19] + aryuZZ[228];
   aryuZZ[144]=aryuZZ[35]*aryuZZ[144];
   aryuZZ[155]= - 253./72. - aryuZZ[36];
   aryuZZ[155]=1./3.*aryuZZ[155] + 41./8.*aryuZZ[35];
   aryuZZ[155]=aryuZZ[18]*aryuZZ[155];
   aryuZZ[157]=49./2.*aryuZZ[18] + 269./3.*aryuZZ[19] - aryuZZ[17];
   aryuZZ[157]=aryuZZ[6]*aryuZZ[157];
   aryuZZ[103]=aryuZZ[105] + aryuZZ[103] + aryuZZ[143] + aryuZZ[140] + 
   aryuZZ[129] + 1./18.*aryuZZ[157] + 1./2.*aryuZZ[155] + 1./9.*
   aryuZZ[144] + aryuZZ[314] + 1./6.*aryuZZ[142] + 17./12.*aryuZZ[21]
    + aryuZZ[324] + 26./3.*aryuZZ[22] - 2107./288.*aryuZZ[50] + 
   aryuZZ[284] + aryuZZ[281] - 47./72.*aryuZZ[68] + 3./4.*aryuZZ[8] - 
   13./12.*aryuZZ[69] + 1949./192.*aryuZZ[42] + aryuZZ[176] - 991./32.*
   aryuZZ[70];
   aryuZZ[103]=aryuZZ[133]*aryuZZ[103];
   aryuZZ[105]= - 43./3. + aryuZZ[162];
   aryuZZ[129]= - 4*aryuZZ[108];
   aryuZZ[140]= - 13./3. - 2*aryuZZ[108];
   aryuZZ[140]=aryuZZ[13]*aryuZZ[140];
   aryuZZ[105]= - 2./3.*aryuZZ[35] + 2*aryuZZ[140] + 1./3.*aryuZZ[105]
    + aryuZZ[129];
   aryuZZ[105]=aryuZZ[35]*aryuZZ[105];
   aryuZZ[142]=aryuZZ[1]*aryuZZ[248];
   aryuZZ[143]=aryuZZ[39]*aryuZZ[248];
   aryuZZ[142]=aryuZZ[142] + 5./3.*aryuZZ[143];
   aryuZZ[142]=aryuZZ[6]*aryuZZ[142];
   aryuZZ[143]=151 + 512./3.*aryuZZ[62];
   aryuZZ[144]= - 1 - aryuZZ[35];
   aryuZZ[155]=aryuZZ[1]*aryuZZ[144];
   aryuZZ[144]=aryuZZ[39]*aryuZZ[144];
   aryuZZ[121]=aryuZZ[121] - aryuZZ[55];
   aryuZZ[121]=MMt*aryuZZ[121];
   aryuZZ[105]=512./9.*aryuZZ[121] + 256./3.*aryuZZ[142] + 1280./27.*
   aryuZZ[144] + 256./9.*aryuZZ[155] + 64*aryuZZ[105] + 128*aryuZZ[140]
    - 256*aryuZZ[108] - 256./3.*aryuZZ[46] - 256./9.*aryuZZ[49] + 58*
   aryuZZ[66] + 1./3.*aryuZZ[143] - 29*aryuZZ[65];
   aryuZZ[105]=MMt*aryuZZ[105];
   aryuZZ[121]= - 167 - 64*aryuZZ[46];
   aryuZZ[142]=aryuZZ[34] - 20./3.*aryuZZ[99];
   aryuZZ[142]=aryuZZ[18]*aryuZZ[142];
   aryuZZ[143]=aryuZZ[1] + 5./3.*aryuZZ[39];
   aryuZZ[143]=aryuZZ[6]*aryuZZ[143];
   aryuZZ[144]=aryuZZ[34] - 16./3.*aryuZZ[99];
   aryuZZ[144]=aryuZZ[37]*aryuZZ[144];
   aryuZZ[121]=32*aryuZZ[144] + 64./3.*aryuZZ[143] + 32*aryuZZ[142] - 
   320./27.*aryuZZ[39] - 64./9.*aryuZZ[1] - 64./3.*aryuZZ[35] + 32*
   aryuZZ[140] + 1./3.*aryuZZ[121] - 64*aryuZZ[108];
   aryuZZ[121]=aryuZZ[37]*aryuZZ[121];
   aryuZZ[140]=32./3.*aryuZZ[71] + 29*aryuZZ[70];
   aryuZZ[115]= - 83./3. + aryuZZ[115];
   aryuZZ[115]=aryuZZ[18]*aryuZZ[115];
   aryuZZ[142]= - 4*aryuZZ[19] - 5./9.*aryuZZ[18];
   aryuZZ[142]=aryuZZ[6]*aryuZZ[142];
   aryuZZ[115]=2*aryuZZ[121] + 4*aryuZZ[142] + 4./3.*aryuZZ[115] + 64*
   aryuZZ[166] - 94./3.*aryuZZ[19] + 172./9.*aryuZZ[21] + 128./9.*
   aryuZZ[50] - 20./9.*aryuZZ[8] - 64./3.*aryuZZ[69] + 2*aryuZZ[140] - 
   29*aryuZZ[42];
   aryuZZ[105]=2*aryuZZ[115] + aryuZZ[105];
   aryuZZ[115]= - 113./2.*aryuZZ[36];
   aryuZZ[121]=1132./9. + aryuZZ[115];
   aryuZZ[140]= - 1501./54.*aryuZZ[35] - 1691./81. + aryuZZ[208];
   aryuZZ[140]=aryuZZ[35]*aryuZZ[140];
   aryuZZ[123]= - 259./9.*aryuZZ[35] - 223./9. + aryuZZ[123];
   aryuZZ[123]=aryuZZ[6]*aryuZZ[123];
   aryuZZ[121]=1./6.*aryuZZ[123] + 1./9.*aryuZZ[121] + aryuZZ[140];
   aryuZZ[123]= - 20./9. + aryuZZ[327];
   aryuZZ[123]= - 1./3.*aryuZZ[6] + 2*aryuZZ[123] + 23./3.*aryuZZ[35];
   aryuZZ[123]=aryuZZ[3]*aryuZZ[37]*aryuZZ[123];
   aryuZZ[121]=1./3.*aryuZZ[121] + aryuZZ[123];
   aryuZZ[121]=MMt*aryuZZ[121];
   aryuZZ[140]=aryuZZ[6]*aryuZZ[369];
   aryuZZ[142]=119./9.*aryuZZ[35] + 1121./27. + aryuZZ[261];
   aryuZZ[142]=aryuZZ[35]*aryuZZ[142];
   aryuZZ[140]=5./2.*aryuZZ[140] + 11./3.*aryuZZ[387] + 1./2.*
   aryuZZ[142];
   aryuZZ[142]= - 5*aryuZZ[6] + 22./3. + aryuZZ[199];
   aryuZZ[142]=aryuZZ[3]*aryuZZ[37]*aryuZZ[142];
   aryuZZ[140]=1./3.*aryuZZ[140] + aryuZZ[142];
   aryuZZ[140]=MMt*aryuZZ[140];
   aryuZZ[140]=25./27.*aryuZZ[148] + aryuZZ[140];
   aryuZZ[140]=aryuZZ[272]*aryuZZ[140];
   aryuZZ[142]= - 223./18.*aryuZZ[6] - 1735./18.*aryuZZ[35] + 1132./27.
    + aryuZZ[220];
   aryuZZ[142]=aryuZZ[37]*aryuZZ[142];
   aryuZZ[121]=1./3.*aryuZZ[140] + 1./9.*aryuZZ[142] + aryuZZ[121];
   aryuZZ[121]=aryuZZ[272]*aryuZZ[121];
   aryuZZ[140]=35./9. + aryuZZ[36];
   aryuZZ[140]=1./2.*aryuZZ[140] + aryuZZ[184];
   aryuZZ[140]=aryuZZ[35]*aryuZZ[140];
   aryuZZ[142]= - 2*aryuZZ[6] + aryuZZ[265] + 1 + aryuZZ[327];
   aryuZZ[142]=aryuZZ[3]*aryuZZ[37]*aryuZZ[142];
   aryuZZ[143]=aryuZZ[6]*aryuZZ[266];
   aryuZZ[140]=2*aryuZZ[142] + aryuZZ[143] + aryuZZ[149] + aryuZZ[140];
   aryuZZ[140]=MMt*aryuZZ[140];
   aryuZZ[134]=aryuZZ[134] - 2./3. - 3./4.*aryuZZ[36];
   aryuZZ[134]=aryuZZ[35]*aryuZZ[134];
   aryuZZ[142]=aryuZZ[6]*aryuZZ[276];
   aryuZZ[143]=aryuZZ[6] + aryuZZ[212] + aryuZZ[35];
   aryuZZ[143]=aryuZZ[3]*aryuZZ[37]*aryuZZ[143];
   aryuZZ[134]=3*aryuZZ[143] + 1./2.*aryuZZ[142] + 7./6.*aryuZZ[36] + 
   aryuZZ[134];
   aryuZZ[134]=MMt*aryuZZ[134];
   aryuZZ[135]=aryuZZ[135] + aryuZZ[156];
   aryuZZ[135]=aryuZZ[37]*aryuZZ[135];
   aryuZZ[134]=aryuZZ[135] + aryuZZ[134];
   aryuZZ[134]=aryuZZ[133]*aryuZZ[134];
   aryuZZ[135]=2*aryuZZ[6] + aryuZZ[149] + aryuZZ[194];
   aryuZZ[135]=aryuZZ[37]*aryuZZ[135];
   aryuZZ[134]=aryuZZ[134] + aryuZZ[135] + aryuZZ[140];
   aryuZZ[134]=aryuZZ[133]*aryuZZ[134];
   aryuZZ[115]=76 + aryuZZ[115];
   aryuZZ[135]= - 53./6.*aryuZZ[35] - 67./9. + aryuZZ[208];
   aryuZZ[135]=aryuZZ[35]*aryuZZ[135];
   aryuZZ[140]= - 43./3.*aryuZZ[35] - 13 + 5./6.*aryuZZ[36];
   aryuZZ[140]=aryuZZ[6]*aryuZZ[140];
   aryuZZ[115]=1./2.*aryuZZ[140] + 1./9.*aryuZZ[115] + aryuZZ[135];
   aryuZZ[115]=1./3.*aryuZZ[115] + aryuZZ[123];
   aryuZZ[115]=MMt*aryuZZ[115];
   aryuZZ[123]= - 79./2.*aryuZZ[35] + 76./3. + aryuZZ[220];
   aryuZZ[123]=1./3.*aryuZZ[123] - 13./2.*aryuZZ[6];
   aryuZZ[123]=aryuZZ[37]*aryuZZ[123];
   aryuZZ[115]=aryuZZ[134] + 1./3.*aryuZZ[123] + aryuZZ[115];
   aryuZZ[115]=aryuZZ[133]*aryuZZ[115];
   aryuZZ[123]=7./3. + aryuZZ[231];
   aryuZZ[123]=aryuZZ[35]*aryuZZ[123];
   aryuZZ[134]=aryuZZ[6]*aryuZZ[248];
   aryuZZ[123]=aryuZZ[134] - 5./3. + aryuZZ[123];
   aryuZZ[123]=MMt*aryuZZ[123];
   aryuZZ[134]=aryuZZ[235] + aryuZZ[6];
   aryuZZ[134]=aryuZZ[37]*aryuZZ[134];
   aryuZZ[123]=aryuZZ[134] + aryuZZ[123];
   aryuZZ[115]=aryuZZ[115] + 256./81.*aryuZZ[123] + aryuZZ[121];
   aryuZZ[115]=aryuZZ[33]*aryuZZ[115];
   aryuZZ[103]=aryuZZ[115] + aryuZZ[103] + 1./9.*aryuZZ[105] + 
   aryuZZ[112];
   aryuZZ[103]=aryuZZ[33]*aryuZZ[103];
   aryuZZ[105]=35./27. + aryuZZ[333];
   aryuZZ[112]=2./3.*aryuZZ[11];
   aryuZZ[115]= - 7./6.*aryuZZ[10];
   aryuZZ[105]=aryuZZ[115] + aryuZZ[112] - 103./72.*aryuZZ[13] + 
   aryuZZ[337] + aryuZZ[306] + 1./4.*aryuZZ[105] + aryuZZ[47];
   aryuZZ[105]=aryuZZ[10]*aryuZZ[105];
   aryuZZ[121]= - 1543./96. + aryuZZ[88];
   aryuZZ[123]= - 23./12. + aryuZZ[13];
   aryuZZ[123]=aryuZZ[13]*aryuZZ[123];
   aryuZZ[134]=49./9. + aryuZZ[13];
   aryuZZ[134]=aryuZZ[11]*aryuZZ[134];
   aryuZZ[135]= - 5./3.*aryuZZ[87];
   aryuZZ[140]= - 9./4.*aryuZZ[86];
   aryuZZ[142]= - 1./6.*aryuZZ[96];
   aryuZZ[143]= - 9./4.*aryuZZ[43];
   aryuZZ[144]= - 13./18.*aryuZZ[12];
   aryuZZ[105]=aryuZZ[105] + 1./12.*aryuZZ[134] + 1./6.*aryuZZ[123] + 
   aryuZZ[144] + aryuZZ[298] + 31./8.*aryuZZ[15] - aryuZZ[47] + 
   aryuZZ[143] + 13./12.*aryuZZ[16] - 5./12.*aryuZZ[89] + aryuZZ[142]
    + aryuZZ[365] + aryuZZ[140] + aryuZZ[135] - 1./24.*aryuZZ[90] - 17./
   12.*aryuZZ[97] + 1./3.*aryuZZ[121] - 1./8.*aryuZZ[85];
   aryuZZ[105]=aryuZZ[73]*aryuZZ[105];
   aryuZZ[121]=aryuZZ[10]*aryuZZ[343];
   aryuZZ[123]=1./2. + aryuZZ[11];
   aryuZZ[134]=aryuZZ[1]*aryuZZ[123];
   aryuZZ[123]=aryuZZ[39]*aryuZZ[123];
   aryuZZ[123]=aryuZZ[134] + 11./9.*aryuZZ[123];
   aryuZZ[123]=aryuZZ[6]*aryuZZ[123];
   aryuZZ[134]= - aryuZZ[11] - 2./3. + aryuZZ[304];
   aryuZZ[134]=aryuZZ[1]*aryuZZ[134];
   aryuZZ[148]= - 11./27.*aryuZZ[11] - 10./27. + aryuZZ[304];
   aryuZZ[148]=aryuZZ[39]*aryuZZ[148];
   aryuZZ[121]=aryuZZ[123] + aryuZZ[121] + 1./3.*aryuZZ[134] + 
   aryuZZ[148];
   aryuZZ[121]=1./3.*aryuZZ[121];
   aryuZZ[123]= - 3*aryuZZ[79];
   aryuZZ[134]=aryuZZ[123] - aryuZZ[80];
   aryuZZ[148]=aryuZZ[134] - 7./12.*aryuZZ[81];
   aryuZZ[149]= - 1./3.*aryuZZ[78];
   aryuZZ[148]=1./2.*aryuZZ[148] + aryuZZ[149];
   aryuZZ[148]=MMH*aryuZZ[73]*aryuZZ[148];
   aryuZZ[105]=aryuZZ[148] + aryuZZ[121] + aryuZZ[105];
   aryuZZ[105]=MMH*aryuZZ[105];
   aryuZZ[148]=3*aryuZZ[85];
   aryuZZ[155]= - 3137./108. + aryuZZ[148];
   aryuZZ[155]= - 1./3.*aryuZZ[90] + 1./2.*aryuZZ[155] - aryuZZ[97];
   aryuZZ[156]= - 1./2.*aryuZZ[47];
   aryuZZ[157]= - 11./3. + aryuZZ[13];
   aryuZZ[157]=aryuZZ[11]*aryuZZ[157];
   aryuZZ[162]=17./36.*aryuZZ[12];
   aryuZZ[155]=1./18.*aryuZZ[157] - 1./9.*aryuZZ[13] + aryuZZ[162] + 
   aryuZZ[298] + 35./8.*aryuZZ[15] + aryuZZ[156] + aryuZZ[143] + 
   aryuZZ[386] - 1./36.*aryuZZ[89] - 7./6.*aryuZZ[96] + aryuZZ[365] + 
   aryuZZ[140] + 1./4.*aryuZZ[155] + aryuZZ[135];
   aryuZZ[157]=3./4.*aryuZZ[45];
   aryuZZ[164]= - 11./12.*aryuZZ[10] + aryuZZ[216] + 13./144.*
   aryuZZ[13] + aryuZZ[162] + aryuZZ[157] + aryuZZ[291] + 52./27. + 
   aryuZZ[287];
   aryuZZ[164]=aryuZZ[10]*aryuZZ[164];
   aryuZZ[155]=1./2.*aryuZZ[155] + aryuZZ[164];
   aryuZZ[155]=aryuZZ[73]*aryuZZ[155];
   aryuZZ[164]=aryuZZ[1]*aryuZZ[218];
   aryuZZ[166]=aryuZZ[39]*aryuZZ[218];
   aryuZZ[167]=aryuZZ[10]*aryuZZ[341];
   aryuZZ[164]=11./4.*aryuZZ[167] + aryuZZ[164] + 11./9.*aryuZZ[166];
   aryuZZ[164]=aryuZZ[6]*aryuZZ[164];
   aryuZZ[166]=1./3.*aryuZZ[80];
   aryuZZ[168]=1./12.*aryuZZ[81];
   aryuZZ[169]=aryuZZ[168] + aryuZZ[123] + aryuZZ[166];
   aryuZZ[169]=1./2.*aryuZZ[169] + aryuZZ[149];
   aryuZZ[169]=MMH*aryuZZ[73]*aryuZZ[169];
   aryuZZ[160]= - aryuZZ[15] + aryuZZ[160] + aryuZZ[28] - 7./9. - 
   aryuZZ[30];
   aryuZZ[176]=1./4.*aryuZZ[160] + aryuZZ[138];
   aryuZZ[180]=aryuZZ[1]*aryuZZ[176];
   aryuZZ[176]=aryuZZ[39]*aryuZZ[176];
   aryuZZ[155]=1./2.*aryuZZ[169] + aryuZZ[155] + 1./3.*aryuZZ[164] + 1./
   9.*aryuZZ[167] + aryuZZ[180] + 11./9.*aryuZZ[176];
   aryuZZ[155]=MMH*aryuZZ[155];
   aryuZZ[164]= - 17./12.*aryuZZ[12];
   aryuZZ[167]=1./24.*aryuZZ[13] + aryuZZ[164] + aryuZZ[300] + 
   aryuZZ[236] + 77./9. + aryuZZ[206];
   aryuZZ[167]=aryuZZ[17]*aryuZZ[167];
   aryuZZ[169]=aryuZZ[132] - 53./12.*aryuZZ[17];
   aryuZZ[169]=aryuZZ[101]*aryuZZ[169];
   aryuZZ[176]= - aryuZZ[18]*aryuZZ[101];
   aryuZZ[180]=17./24.*aryuZZ[176];
   aryuZZ[181]=aryuZZ[180] - 7./3.*aryuZZ[10] + aryuZZ[169] + 
   aryuZZ[260] + 67./72.*aryuZZ[13] - 17./72.*aryuZZ[12] - 13./16.*
   aryuZZ[45] - 7./8.*aryuZZ[47] - 29./27. - 63./16.*aryuZZ[43];
   aryuZZ[181]=aryuZZ[18]*aryuZZ[181];
   aryuZZ[183]= - 9*aryuZZ[74];
   aryuZZ[184]= - 11./3.*aryuZZ[91];
   aryuZZ[188]=aryuZZ[184] + 15./2.*aryuZZ[77] + aryuZZ[183] + 1./18.*
   aryuZZ[92];
   aryuZZ[188]=1./2.*aryuZZ[188] - aryuZZ[76];
   aryuZZ[189]= - 1./6.*aryuZZ[75];
   aryuZZ[188]=aryuZZ[189] + 1./2.*aryuZZ[188] - 7./3.*aryuZZ[93];
   aryuZZ[190]=aryuZZ[18] - 37*aryuZZ[19] + aryuZZ[219];
   aryuZZ[190]=aryuZZ[18]*aryuZZ[190];
   aryuZZ[190]=aryuZZ[286] + 1./4.*aryuZZ[190];
   aryuZZ[190]=aryuZZ[3]*aryuZZ[190];
   aryuZZ[192]=aryuZZ[209] - 2*aryuZZ[286];
   aryuZZ[192]=aryuZZ[101]*aryuZZ[192];
   aryuZZ[194]= - 191./2.*aryuZZ[13] + 1225./3. - 29*aryuZZ[12];
   aryuZZ[194]=aryuZZ[19]*aryuZZ[194];
   aryuZZ[196]= - 43./4.*aryuZZ[19] + aryuZZ[17];
   aryuZZ[196]=aryuZZ[10]*aryuZZ[196];
   aryuZZ[167]=1./2.*aryuZZ[190] + aryuZZ[181] + 1./3.*aryuZZ[196] + 
   aryuZZ[192] + 1./2.*aryuZZ[167] + 1./72.*aryuZZ[194] + 41./48.*
   aryuZZ[21] + 221./48.*aryuZZ[20] + 1./2.*aryuZZ[188] + 2./3.*
   aryuZZ[22];
   aryuZZ[167]=aryuZZ[73]*aryuZZ[167];
   aryuZZ[181]= - aryuZZ[1]*aryuZZ[19];
   aryuZZ[188]= - aryuZZ[39]*aryuZZ[19];
   aryuZZ[190]=aryuZZ[18]*aryuZZ[230];
   aryuZZ[194]=aryuZZ[190] + aryuZZ[181] + 11./9.*aryuZZ[188];
   aryuZZ[196]=aryuZZ[1]*aryuZZ[19];
   aryuZZ[199]=aryuZZ[39]*aryuZZ[19];
   aryuZZ[196]=aryuZZ[362] + aryuZZ[196] + 11./9.*aryuZZ[199];
   aryuZZ[196]=aryuZZ[6]*aryuZZ[196];
   aryuZZ[194]=1./3.*aryuZZ[194] + aryuZZ[196];
   aryuZZ[196]= - 5./2.*aryuZZ[13] + 71./3. - aryuZZ[12];
   aryuZZ[196]=aryuZZ[19]*aryuZZ[196];
   aryuZZ[199]= - 1./2.*aryuZZ[12];
   aryuZZ[200]= - 1 + aryuZZ[199];
   aryuZZ[200]=aryuZZ[17]*aryuZZ[200];
   aryuZZ[202]=5./4.*aryuZZ[13] + 5./3. + aryuZZ[199];
   aryuZZ[202]=1./3.*aryuZZ[202] - aryuZZ[10];
   aryuZZ[202]=aryuZZ[18]*aryuZZ[202];
   aryuZZ[196]=aryuZZ[202] + aryuZZ[372] + aryuZZ[200] + 1./6.*
   aryuZZ[196] + 1./2.*aryuZZ[21] + aryuZZ[22] + 1./2.*aryuZZ[77] - 
   aryuZZ[93];
   aryuZZ[200]=aryuZZ[3]*aryuZZ[18]*aryuZZ[357];
   aryuZZ[196]=1./3.*aryuZZ[196] + 1./2.*aryuZZ[200];
   aryuZZ[196]=aryuZZ[73]*aryuZZ[196];
   aryuZZ[200]=aryuZZ[10]*aryuZZ[12];
   aryuZZ[200]=1./3.*aryuZZ[200] + 1./6.*aryuZZ[12] - aryuZZ[11];
   aryuZZ[200]=MMH*aryuZZ[73]*aryuZZ[200];
   aryuZZ[194]=1./3.*aryuZZ[200] + 1./3.*aryuZZ[194] + aryuZZ[196];
   aryuZZ[194]=aryuZZ[272]*aryuZZ[194];
   aryuZZ[196]= - 43./3.*aryuZZ[19] + aryuZZ[17];
   aryuZZ[196]=aryuZZ[1]*aryuZZ[196];
   aryuZZ[200]= - 377./3.*aryuZZ[19] + aryuZZ[332];
   aryuZZ[200]=aryuZZ[39]*aryuZZ[200];
   aryuZZ[196]=aryuZZ[196] + 1./9.*aryuZZ[200];
   aryuZZ[200]= - 23*aryuZZ[1] - 359./27.*aryuZZ[39];
   aryuZZ[200]=aryuZZ[18]*aryuZZ[200];
   aryuZZ[196]=1./4.*aryuZZ[196] + 1./3.*aryuZZ[200];
   aryuZZ[196]=aryuZZ[6]*aryuZZ[196];
   aryuZZ[200]= - 1./4.*aryuZZ[20];
   aryuZZ[202]=1./8.*aryuZZ[17] + 67./36.*aryuZZ[19] - 29./4.*
   aryuZZ[21] + aryuZZ[200] - 31./4.*aryuZZ[8] - aryuZZ[22];
   aryuZZ[202]=aryuZZ[1]*aryuZZ[202];
   aryuZZ[203]=11./24.*aryuZZ[17] + 449./108.*aryuZZ[19] - 1271./108.*
   aryuZZ[21] - 11./12.*aryuZZ[20] - 1469./108.*aryuZZ[8] - aryuZZ[22];
   aryuZZ[203]=aryuZZ[39]*aryuZZ[203];
   aryuZZ[205]=713*aryuZZ[1] + 3683./9.*aryuZZ[39];
   aryuZZ[205]=aryuZZ[18]*aryuZZ[205];
   aryuZZ[155]=1./4.*aryuZZ[194] + aryuZZ[155] + aryuZZ[167] + 
   aryuZZ[196] + 1./36.*aryuZZ[205] + aryuZZ[202] + 1./3.*aryuZZ[203];
   aryuZZ[155]=aryuZZ[272]*aryuZZ[155];
   aryuZZ[167]= - 63*aryuZZ[43];
   aryuZZ[194]= - 391./27. + aryuZZ[167];
   aryuZZ[196]= - 7*aryuZZ[47];
   aryuZZ[202]= - 13./2.*aryuZZ[45];
   aryuZZ[194]=aryuZZ[202] + 1./2.*aryuZZ[194] + aryuZZ[196];
   aryuZZ[194]=107./9.*aryuZZ[13] + 1./2.*aryuZZ[194] + aryuZZ[336];
   aryuZZ[203]=4*aryuZZ[19] - 53./6.*aryuZZ[17];
   aryuZZ[203]=aryuZZ[101]*aryuZZ[203];
   aryuZZ[205]= - 5./3.*aryuZZ[10];
   aryuZZ[176]=17./12.*aryuZZ[176];
   aryuZZ[207]= - aryuZZ[19]*aryuZZ[100];
   aryuZZ[208]=3./2.*aryuZZ[207];
   aryuZZ[194]=aryuZZ[176] + aryuZZ[205] + aryuZZ[203] + aryuZZ[208] + 
   1./2.*aryuZZ[194] + aryuZZ[216];
   aryuZZ[194]=aryuZZ[18]*aryuZZ[194];
   aryuZZ[212]= - 27./2.*aryuZZ[43];
   aryuZZ[213]= - 3*aryuZZ[47];
   aryuZZ[214]=aryuZZ[335] + aryuZZ[213] + 557./3. + aryuZZ[212];
   aryuZZ[218]= - 3*aryuZZ[12];
   aryuZZ[214]=11./12.*aryuZZ[207] - 7./6.*aryuZZ[11] - 503./72.*
   aryuZZ[13] + 1./4.*aryuZZ[214] + aryuZZ[218];
   aryuZZ[214]=aryuZZ[19]*aryuZZ[214];
   aryuZZ[220]=431./9. + aryuZZ[333];
   aryuZZ[221]=13./3.*aryuZZ[12];
   aryuZZ[220]=aryuZZ[221] + 1./2.*aryuZZ[220] + aryuZZ[45];
   aryuZZ[222]=aryuZZ[19]*aryuZZ[100];
   aryuZZ[220]=1./6.*aryuZZ[222] + aryuZZ[118] + 1./2.*aryuZZ[220] - 1./
   3.*aryuZZ[13];
   aryuZZ[220]=aryuZZ[17]*aryuZZ[220];
   aryuZZ[223]= - 3./4.*aryuZZ[19] + aryuZZ[17];
   aryuZZ[223]=aryuZZ[17]*aryuZZ[223];
   aryuZZ[136]=5*aryuZZ[136] + 17*aryuZZ[18];
   aryuZZ[136]=1./4.*aryuZZ[18]*aryuZZ[136];
   aryuZZ[225]=aryuZZ[136] - 1./2.*aryuZZ[209] + aryuZZ[223];
   aryuZZ[225]=aryuZZ[3]*aryuZZ[225];
   aryuZZ[226]=2*aryuZZ[192];
   aryuZZ[228]= - 9./4.*aryuZZ[74];
   aryuZZ[229]= - 11./12.*aryuZZ[91];
   aryuZZ[230]= - 7*aryuZZ[19] + aryuZZ[119];
   aryuZZ[230]=1./3.*aryuZZ[10]*aryuZZ[230];
   aryuZZ[194]=aryuZZ[225] + aryuZZ[194] + aryuZZ[230] + aryuZZ[226] + 
   aryuZZ[220] + aryuZZ[214] + 33./4.*aryuZZ[21] + 511./72.*aryuZZ[20]
    + 67./12.*aryuZZ[22] + aryuZZ[189] - 151./12.*aryuZZ[93] + 2./3.*
   aryuZZ[76] + aryuZZ[229] + 67./12.*aryuZZ[77] + 41./72.*aryuZZ[92]
    + 10./9.*aryuZZ[94] + aryuZZ[228];
   aryuZZ[194]=aryuZZ[73]*aryuZZ[194];
   aryuZZ[214]= - 2*aryuZZ[9];
   aryuZZ[220]=685./9.*aryuZZ[21] + aryuZZ[107] + aryuZZ[214] + 685./9.
   *aryuZZ[8];
   aryuZZ[225]= - 1./3.*aryuZZ[7];
   aryuZZ[231]= - 1 + aryuZZ[225];
   aryuZZ[231]=aryuZZ[19]*aryuZZ[231];
   aryuZZ[232]=10./27. + aryuZZ[302];
   aryuZZ[232]=aryuZZ[17]*aryuZZ[232];
   aryuZZ[220]=aryuZZ[232] + 1./9.*aryuZZ[220] + 1./2.*aryuZZ[231];
   aryuZZ[220]=aryuZZ[39]*aryuZZ[220];
   aryuZZ[231]=55./9.*aryuZZ[19] - aryuZZ[17];
   aryuZZ[231]=aryuZZ[1]*aryuZZ[231];
   aryuZZ[234]=aryuZZ[39]*aryuZZ[254];
   aryuZZ[235]=293*aryuZZ[1] + 1631./9.*aryuZZ[39];
   aryuZZ[235]=aryuZZ[18]*aryuZZ[235];
   aryuZZ[231]=1./9.*aryuZZ[235] + aryuZZ[231] + 11./9.*aryuZZ[234];
   aryuZZ[231]=aryuZZ[6]*aryuZZ[231];
   aryuZZ[225]=7 + aryuZZ[225];
   aryuZZ[225]=aryuZZ[19]*aryuZZ[225];
   aryuZZ[234]=2./3. + aryuZZ[302];
   aryuZZ[235]=aryuZZ[17]*aryuZZ[234];
   aryuZZ[225]=1./3.*aryuZZ[235] + 1./6.*aryuZZ[225] + 15*aryuZZ[21] - 
   2./3.*aryuZZ[22] - 2./3.*aryuZZ[9] + 15*aryuZZ[8];
   aryuZZ[225]=aryuZZ[1]*aryuZZ[225];
   aryuZZ[237]= - 1064./3. + aryuZZ[302];
   aryuZZ[237]=aryuZZ[1]*aryuZZ[237];
   aryuZZ[238]= - 1820./27. + aryuZZ[302];
   aryuZZ[238]=aryuZZ[39]*aryuZZ[238];
   aryuZZ[237]=1./3.*aryuZZ[237] + aryuZZ[238];
   aryuZZ[237]=aryuZZ[18]*aryuZZ[237];
   aryuZZ[105]=aryuZZ[155] + aryuZZ[105] + aryuZZ[194] + 1./2.*
   aryuZZ[231] + 1./3.*aryuZZ[237] + aryuZZ[225] + aryuZZ[220];
   aryuZZ[105]=aryuZZ[272]*aryuZZ[105];
   aryuZZ[155]= - 469./27. + aryuZZ[333];
   aryuZZ[112]=aryuZZ[115] + aryuZZ[112] - 367./72.*aryuZZ[13] + 
   aryuZZ[337] + aryuZZ[306] + 1./4.*aryuZZ[155] + aryuZZ[47];
   aryuZZ[112]=aryuZZ[10]*aryuZZ[112];
   aryuZZ[115]= - 959./96. + aryuZZ[88];
   aryuZZ[155]= - 167./36. - aryuZZ[13];
   aryuZZ[155]=aryuZZ[13]*aryuZZ[155];
   aryuZZ[194]= - 167./9. - 23*aryuZZ[13];
   aryuZZ[194]=aryuZZ[11]*aryuZZ[194];
   aryuZZ[112]=aryuZZ[112] + 1./12.*aryuZZ[194] + 1./2.*aryuZZ[155] + 
   aryuZZ[144] + aryuZZ[298] + 39./8.*aryuZZ[15] - aryuZZ[47] + 
   aryuZZ[143] - 1./4.*aryuZZ[16] - 197./12.*aryuZZ[89] + aryuZZ[142]
    + aryuZZ[365] + aryuZZ[140] + aryuZZ[135] + 89./8.*aryuZZ[90] - 1./
   12.*aryuZZ[97] + 1./3.*aryuZZ[115] - 9./8.*aryuZZ[85];
   aryuZZ[112]=aryuZZ[73]*aryuZZ[112];
   aryuZZ[115]=aryuZZ[134] + 17./12.*aryuZZ[81];
   aryuZZ[115]=1./2.*aryuZZ[115] + aryuZZ[149];
   aryuZZ[115]=MMH*aryuZZ[73]*aryuZZ[115];
   aryuZZ[112]=aryuZZ[115] + aryuZZ[121] + aryuZZ[112];
   aryuZZ[112]=MMH*aryuZZ[112];
   aryuZZ[115]= - 85759./27. + aryuZZ[167];
   aryuZZ[115]=aryuZZ[202] + 1./2.*aryuZZ[115] + aryuZZ[196];
   aryuZZ[115]= - 2335./9.*aryuZZ[13] + 1./2.*aryuZZ[115] + aryuZZ[336]
   ;
   aryuZZ[115]=aryuZZ[176] + aryuZZ[205] + aryuZZ[203] + aryuZZ[208] + 
   1./2.*aryuZZ[115] + aryuZZ[216];
   aryuZZ[115]=aryuZZ[18]*aryuZZ[115];
   aryuZZ[121]=aryuZZ[335] + aryuZZ[213] + 26863./9. + aryuZZ[212];
   aryuZZ[121]=59./12.*aryuZZ[207] + 29./6.*aryuZZ[11] - 1423./8.*
   aryuZZ[13] + 1./4.*aryuZZ[121] - 247./3.*aryuZZ[12];
   aryuZZ[121]=aryuZZ[19]*aryuZZ[121];
   aryuZZ[134]=529./3. + aryuZZ[333];
   aryuZZ[134]=aryuZZ[221] + 1./2.*aryuZZ[134] + aryuZZ[45];
   aryuZZ[134]=25./6.*aryuZZ[222] + aryuZZ[118] + 1./2.*aryuZZ[134] + 4.
   /3.*aryuZZ[13];
   aryuZZ[134]=aryuZZ[17]*aryuZZ[134];
   aryuZZ[135]=aryuZZ[136] + 47./2.*aryuZZ[209] + aryuZZ[223];
   aryuZZ[135]=aryuZZ[3]*aryuZZ[135];
   aryuZZ[115]=aryuZZ[135] + aryuZZ[115] + aryuZZ[230] + aryuZZ[226] + 
   aryuZZ[134] + aryuZZ[121] + 3467./12.*aryuZZ[21] + 101./24.*
   aryuZZ[20] + 7099./12.*aryuZZ[22] + aryuZZ[189] - 2087./12.*
   aryuZZ[93] + 40./3.*aryuZZ[76] + aryuZZ[229] - 455./4.*aryuZZ[77] - 
   605./24.*aryuZZ[92] + 2*aryuZZ[94] + aryuZZ[228];
   aryuZZ[115]=aryuZZ[73]*aryuZZ[115];
   aryuZZ[118]=3*aryuZZ[207] + aryuZZ[118] - 95./16.*aryuZZ[13] + 55./8.
   *aryuZZ[12] + aryuZZ[251] + aryuZZ[236] - 164./9. + aryuZZ[287];
   aryuZZ[118]=aryuZZ[17]*aryuZZ[118];
   aryuZZ[121]=137 - 7./2.*aryuZZ[43];
   aryuZZ[121]=107*aryuZZ[13] + 55./3.*aryuZZ[12] + aryuZZ[202] + 9*
   aryuZZ[121] - 13*aryuZZ[47];
   aryuZZ[121]=aryuZZ[222] + 1./4.*aryuZZ[121] + aryuZZ[151];
   aryuZZ[121]=aryuZZ[180] + aryuZZ[269] + 1./2.*aryuZZ[121] + 
   aryuZZ[169];
   aryuZZ[121]=aryuZZ[18]*aryuZZ[121];
   aryuZZ[134]= - 14143./3. - 27*aryuZZ[43];
   aryuZZ[134]=6091./12.*aryuZZ[13] + 253./2.*aryuZZ[12] + aryuZZ[335]
    + 1./2.*aryuZZ[134] + aryuZZ[213];
   aryuZZ[134]=13./2.*aryuZZ[222] + 1./2.*aryuZZ[134] - 13*aryuZZ[11];
   aryuZZ[134]=aryuZZ[19]*aryuZZ[134];
   aryuZZ[135]= - 3./2.*aryuZZ[19];
   aryuZZ[136]=aryuZZ[135] + aryuZZ[17];
   aryuZZ[136]=aryuZZ[17]*aryuZZ[136];
   aryuZZ[142]=33*aryuZZ[18] - 9*aryuZZ[19] + aryuZZ[219];
   aryuZZ[142]=aryuZZ[18]*aryuZZ[142];
   aryuZZ[136]=1./4.*aryuZZ[142] - 33*aryuZZ[209] + aryuZZ[136];
   aryuZZ[136]=aryuZZ[3]*aryuZZ[136];
   aryuZZ[142]=aryuZZ[184] + 489./2.*aryuZZ[77] + aryuZZ[183] + 197./2.
   *aryuZZ[92];
   aryuZZ[142]=1./4.*aryuZZ[142] - 11*aryuZZ[76];
   aryuZZ[143]= - 53./2.*aryuZZ[19] + aryuZZ[311];
   aryuZZ[143]=aryuZZ[10]*aryuZZ[143];
   aryuZZ[118]=1./2.*aryuZZ[136] + aryuZZ[121] + 1./2.*aryuZZ[143] + 
   aryuZZ[192] + aryuZZ[118] + 1./2.*aryuZZ[134] - 3023./16.*aryuZZ[21]
    + 151./16.*aryuZZ[20] - 1571./4.*aryuZZ[22] + aryuZZ[309] + 1./2.*
   aryuZZ[142] + 159*aryuZZ[93];
   aryuZZ[118]=aryuZZ[73]*aryuZZ[118];
   aryuZZ[121]=491./24.*aryuZZ[13] - 55./6.*aryuZZ[12] + aryuZZ[306] + 
   aryuZZ[305] + 503./27. + aryuZZ[206];
   aryuZZ[121]= - 1./4.*aryuZZ[10] + 1./2.*aryuZZ[121] + aryuZZ[216];
   aryuZZ[121]=aryuZZ[10]*aryuZZ[121];
   aryuZZ[134]= - 5./2.*aryuZZ[87] + 29./8.*aryuZZ[97] + 29./16.*
   aryuZZ[85] + 1667./576. + 2*aryuZZ[88];
   aryuZZ[136]=115./4. + aryuZZ[13];
   aryuZZ[136]=aryuZZ[13]*aryuZZ[136];
   aryuZZ[142]= - 1 + aryuZZ[210];
   aryuZZ[142]=aryuZZ[11]*aryuZZ[142];
   aryuZZ[121]=aryuZZ[121] + 1./3.*aryuZZ[142] + 1./6.*aryuZZ[136] - 55.
   /24.*aryuZZ[12] + aryuZZ[233] + 137./48.*aryuZZ[15] + aryuZZ[191] + 
   aryuZZ[171] - 15./8.*aryuZZ[16] + 95./8.*aryuZZ[89] - 5./3.*
   aryuZZ[96] + aryuZZ[137] + 1./3.*aryuZZ[134] - 9./8.*aryuZZ[86];
   aryuZZ[121]=aryuZZ[73]*aryuZZ[121];
   aryuZZ[134]=1./2.*aryuZZ[160] + aryuZZ[187];
   aryuZZ[136]=aryuZZ[1]*aryuZZ[134];
   aryuZZ[134]=aryuZZ[39]*aryuZZ[134];
   aryuZZ[134]=1./3.*aryuZZ[136] + aryuZZ[134];
   aryuZZ[136]= - 2./3. + aryuZZ[7];
   aryuZZ[137]=aryuZZ[1]*aryuZZ[136];
   aryuZZ[136]=aryuZZ[39]*aryuZZ[136];
   aryuZZ[136]=1./3.*aryuZZ[137] + aryuZZ[136];
   aryuZZ[136]=aryuZZ[10]*aryuZZ[136];
   aryuZZ[137]=aryuZZ[139] + 11./4.*aryuZZ[322];
   aryuZZ[137]=aryuZZ[6]*aryuZZ[137];
   aryuZZ[123]= - 11./12.*aryuZZ[81] + aryuZZ[123] - 11./3.*aryuZZ[80];
   aryuZZ[123]=1./2.*aryuZZ[123] + aryuZZ[149];
   aryuZZ[123]=MMH*aryuZZ[73]*aryuZZ[123];
   aryuZZ[121]=1./2.*aryuZZ[123] + aryuZZ[121] + 1./3.*aryuZZ[137] + 1./
   2.*aryuZZ[134] + 1./3.*aryuZZ[136];
   aryuZZ[121]=MMH*aryuZZ[121];
   aryuZZ[123]=aryuZZ[1]*aryuZZ[170];
   aryuZZ[134]=aryuZZ[39]*aryuZZ[170];
   aryuZZ[123]=19./6.*aryuZZ[290] + 1./3.*aryuZZ[123] + aryuZZ[134];
   aryuZZ[123]=aryuZZ[6]*aryuZZ[123];
   aryuZZ[134]= - 17./18. + aryuZZ[179];
   aryuZZ[134]=aryuZZ[19]*aryuZZ[134];
   aryuZZ[134]=aryuZZ[134] + aryuZZ[173];
   aryuZZ[136]=aryuZZ[1]*aryuZZ[134];
   aryuZZ[134]=aryuZZ[39]*aryuZZ[134];
   aryuZZ[137]=17./6. + aryuZZ[7];
   aryuZZ[139]=aryuZZ[1]*aryuZZ[137];
   aryuZZ[137]=aryuZZ[39]*aryuZZ[137];
   aryuZZ[137]=1./3.*aryuZZ[139] + aryuZZ[137];
   aryuZZ[137]=aryuZZ[18]*aryuZZ[137];
   aryuZZ[123]=aryuZZ[123] + 1./3.*aryuZZ[137] + 1./3.*aryuZZ[136] + 
   aryuZZ[134];
   aryuZZ[134]= - 517./2.*aryuZZ[13] - 1057./27. + 165*aryuZZ[12];
   aryuZZ[134]=aryuZZ[208] + 1./4.*aryuZZ[134] + aryuZZ[313];
   aryuZZ[134]=aryuZZ[19]*aryuZZ[134];
   aryuZZ[136]=aryuZZ[17]*aryuZZ[186];
   aryuZZ[137]=aryuZZ[101]*aryuZZ[209];
   aryuZZ[139]=59./6.*aryuZZ[19] - aryuZZ[17];
   aryuZZ[139]=aryuZZ[10]*aryuZZ[139];
   aryuZZ[134]=aryuZZ[139] + 3./2.*aryuZZ[137] + aryuZZ[134] + 
   aryuZZ[136];
   aryuZZ[136]=209./2.*aryuZZ[13] + 1057./27. - 11*aryuZZ[12];
   aryuZZ[137]=aryuZZ[101]*aryuZZ[19];
   aryuZZ[136]= - 13./12.*aryuZZ[10] + 3./2.*aryuZZ[137] + aryuZZ[208]
    + 1./8.*aryuZZ[136] + aryuZZ[151];
   aryuZZ[136]=aryuZZ[18]*aryuZZ[136];
   aryuZZ[137]= - aryuZZ[19] - aryuZZ[18];
   aryuZZ[137]=aryuZZ[18]*aryuZZ[137];
   aryuZZ[137]=aryuZZ[209] + 1./2.*aryuZZ[137];
   aryuZZ[137]=aryuZZ[3]*aryuZZ[137];
   aryuZZ[134]=17./4.*aryuZZ[137] + 1./2.*aryuZZ[134] + aryuZZ[136];
   aryuZZ[134]=aryuZZ[73]*aryuZZ[134];
   aryuZZ[136]=aryuZZ[1]*aryuZZ[178];
   aryuZZ[137]=aryuZZ[39]*aryuZZ[178];
   aryuZZ[136]=2*aryuZZ[317] + 1./3.*aryuZZ[136] + aryuZZ[137];
   aryuZZ[136]=aryuZZ[6]*aryuZZ[136];
   aryuZZ[137]=aryuZZ[302] + aryuZZ[216];
   aryuZZ[139]=aryuZZ[1]*aryuZZ[137];
   aryuZZ[137]=aryuZZ[39]*aryuZZ[137];
   aryuZZ[142]=aryuZZ[10]*aryuZZ[318];
   aryuZZ[136]=aryuZZ[136] + aryuZZ[142] + 1./3.*aryuZZ[139] + 
   aryuZZ[137];
   aryuZZ[137]=127./27. + aryuZZ[323];
   aryuZZ[137]=aryuZZ[11]*aryuZZ[137];
   aryuZZ[137]=11./2.*aryuZZ[147] + aryuZZ[137];
   aryuZZ[139]= - 127./27. + aryuZZ[163];
   aryuZZ[139]=1./2.*aryuZZ[139] - 11*aryuZZ[13];
   aryuZZ[139]=aryuZZ[269] + 1./2.*aryuZZ[139] - 2./3.*aryuZZ[11];
   aryuZZ[139]=aryuZZ[10]*aryuZZ[139];
   aryuZZ[137]=1./4.*aryuZZ[137] + aryuZZ[139];
   aryuZZ[137]=aryuZZ[73]*aryuZZ[137];
   aryuZZ[136]=1./3.*aryuZZ[136] + aryuZZ[137];
   aryuZZ[136]=MMH*aryuZZ[136];
   aryuZZ[123]=aryuZZ[136] + 1./2.*aryuZZ[123] + aryuZZ[134];
   aryuZZ[123]=aryuZZ[133]*aryuZZ[123];
   aryuZZ[111]= - aryuZZ[9] + aryuZZ[111];
   aryuZZ[111]=aryuZZ[200] + 11./2.*aryuZZ[111] - 13*aryuZZ[22];
   aryuZZ[134]=13 - 11./3.*aryuZZ[7];
   aryuZZ[134]=aryuZZ[19]*aryuZZ[134];
   aryuZZ[136]=1./4. - aryuZZ[7];
   aryuZZ[136]=aryuZZ[17]*aryuZZ[136];
   aryuZZ[134]=1./6.*aryuZZ[136] + 1./2.*aryuZZ[134] + 1./3.*
   aryuZZ[111] - 3./4.*aryuZZ[21];
   aryuZZ[134]=aryuZZ[1]*aryuZZ[134];
   aryuZZ[137]= - 11*aryuZZ[7];
   aryuZZ[139]=39 + aryuZZ[137];
   aryuZZ[139]=aryuZZ[19]*aryuZZ[139];
   aryuZZ[111]=1./2.*aryuZZ[136] + 1./2.*aryuZZ[139] + aryuZZ[111] - 9./
   4.*aryuZZ[21];
   aryuZZ[111]=aryuZZ[39]*aryuZZ[111];
   aryuZZ[136]=71./2. - aryuZZ[7];
   aryuZZ[139]=aryuZZ[1]*aryuZZ[136];
   aryuZZ[136]=aryuZZ[39]*aryuZZ[136];
   aryuZZ[136]=1./3.*aryuZZ[139] + aryuZZ[136];
   aryuZZ[136]=aryuZZ[18]*aryuZZ[136];
   aryuZZ[139]= - 17*aryuZZ[19] + aryuZZ[17];
   aryuZZ[142]=aryuZZ[1]*aryuZZ[139];
   aryuZZ[139]=aryuZZ[39]*aryuZZ[139];
   aryuZZ[139]=aryuZZ[142] + 3*aryuZZ[139];
   aryuZZ[139]=1./2.*aryuZZ[139] + 1./3.*aryuZZ[197];
   aryuZZ[139]=aryuZZ[6]*aryuZZ[139];
   aryuZZ[111]=aryuZZ[123] + aryuZZ[121] + aryuZZ[118] + 1./2.*
   aryuZZ[139] + 1./6.*aryuZZ[136] + aryuZZ[134] + aryuZZ[111];
   aryuZZ[111]=aryuZZ[133]*aryuZZ[111];
   aryuZZ[118]=35*aryuZZ[7];
   aryuZZ[121]= - 449./3. + aryuZZ[118];
   aryuZZ[121]=aryuZZ[19]*aryuZZ[121];
   aryuZZ[123]=5*aryuZZ[21];
   aryuZZ[134]=aryuZZ[123] + 8*aryuZZ[22] + 14*aryuZZ[9] + 5*aryuZZ[8];
   aryuZZ[121]=aryuZZ[235] + aryuZZ[134] + 1./6.*aryuZZ[121];
   aryuZZ[121]=aryuZZ[1]*aryuZZ[121];
   aryuZZ[118]= - 3937./27. + aryuZZ[118];
   aryuZZ[118]=aryuZZ[19]*aryuZZ[118];
   aryuZZ[118]=aryuZZ[232] + aryuZZ[134] + 1./6.*aryuZZ[118];
   aryuZZ[118]=aryuZZ[39]*aryuZZ[118];
   aryuZZ[134]=37./3.*aryuZZ[19] - aryuZZ[17];
   aryuZZ[134]=aryuZZ[1]*aryuZZ[134];
   aryuZZ[136]=895./3.*aryuZZ[19] + aryuZZ[295];
   aryuZZ[136]=aryuZZ[39]*aryuZZ[136];
   aryuZZ[139]=53*aryuZZ[1] + 119*aryuZZ[39];
   aryuZZ[139]=aryuZZ[18]*aryuZZ[139];
   aryuZZ[134]=1./9.*aryuZZ[139] + aryuZZ[134] + 1./9.*aryuZZ[136];
   aryuZZ[134]=aryuZZ[6]*aryuZZ[134];
   aryuZZ[136]= - 128./3. + aryuZZ[302];
   aryuZZ[136]=aryuZZ[1]*aryuZZ[136];
   aryuZZ[139]= - 364./9. + aryuZZ[302];
   aryuZZ[139]=aryuZZ[39]*aryuZZ[139];
   aryuZZ[136]=1./3.*aryuZZ[136] + aryuZZ[139];
   aryuZZ[136]=aryuZZ[18]*aryuZZ[136];
   aryuZZ[111]=aryuZZ[111] + aryuZZ[112] + aryuZZ[115] + 1./2.*
   aryuZZ[134] + 1./3.*aryuZZ[136] + 1./3.*aryuZZ[121] + aryuZZ[118];
   aryuZZ[111]=aryuZZ[133]*aryuZZ[111];
   aryuZZ[112]=4*aryuZZ[108];
   aryuZZ[115]=55./3. + aryuZZ[112];
   aryuZZ[115]=aryuZZ[13]*aryuZZ[115];
   aryuZZ[115]=aryuZZ[115] + 45 + 8*aryuZZ[108];
   aryuZZ[115]=aryuZZ[18]*aryuZZ[115];
   aryuZZ[118]=23./3.*aryuZZ[77] - 1./9.*aryuZZ[94] + aryuZZ[92];
   aryuZZ[121]=aryuZZ[108]*aryuZZ[77];
   aryuZZ[134]= - 247./9. + aryuZZ[129];
   aryuZZ[134]=aryuZZ[22]*aryuZZ[134];
   aryuZZ[136]= - 20./3. - aryuZZ[108];
   aryuZZ[136]=aryuZZ[21]*aryuZZ[136];
   aryuZZ[139]=113./9. + aryuZZ[112];
   aryuZZ[139]=aryuZZ[13]*aryuZZ[139];
   aryuZZ[139]=aryuZZ[139] - 233./9. + 3*aryuZZ[12];
   aryuZZ[139]=aryuZZ[19]*aryuZZ[139];
   aryuZZ[115]=aryuZZ[115] + aryuZZ[267] + 2*aryuZZ[139] + 4*
   aryuZZ[136] + 2*aryuZZ[134] + 4*aryuZZ[121] + 34./3.*aryuZZ[93] + 2*
   aryuZZ[118] - aryuZZ[76];
   aryuZZ[115]=aryuZZ[73]*aryuZZ[115];
   aryuZZ[118]=aryuZZ[181] + 5./3.*aryuZZ[188];
   aryuZZ[121]= - aryuZZ[1] - 17./27.*aryuZZ[39];
   aryuZZ[121]=aryuZZ[18]*aryuZZ[121];
   aryuZZ[118]=4*aryuZZ[118] + 5*aryuZZ[121];
   aryuZZ[118]=aryuZZ[6]*aryuZZ[118];
   aryuZZ[121]= - 4*aryuZZ[21] - aryuZZ[9] - 4*aryuZZ[8];
   aryuZZ[121]=5*aryuZZ[121] + 31./3.*aryuZZ[19];
   aryuZZ[121]=aryuZZ[1]*aryuZZ[121];
   aryuZZ[134]=179./9.*aryuZZ[19] - 340./27.*aryuZZ[21] - 11*aryuZZ[9]
    - 340./27.*aryuZZ[8];
   aryuZZ[134]=aryuZZ[39]*aryuZZ[134];
   aryuZZ[136]=aryuZZ[1] + 17./27.*aryuZZ[39];
   aryuZZ[136]=aryuZZ[18]*aryuZZ[136];
   aryuZZ[118]=4*aryuZZ[118] + 52*aryuZZ[136] + aryuZZ[121] + 
   aryuZZ[134];
   aryuZZ[115]=1./3.*aryuZZ[118] + 2*aryuZZ[115];
   aryuZZ[118]= - 5 - 14*aryuZZ[90];
   aryuZZ[118]=1./3.*aryuZZ[118] + 4*aryuZZ[89];
   aryuZZ[118]=MMH*aryuZZ[73]*aryuZZ[118];
   aryuZZ[103]=aryuZZ[104] + aryuZZ[103] + aryuZZ[111] + aryuZZ[105] + 
   2*aryuZZ[115] + aryuZZ[118];
   aryuZZ[103]=aryuZZ[4]*aryuZZ[103];
   aryuZZ[104]=11*aryuZZ[85];
   aryuZZ[105]=55*aryuZZ[90] - 501./4. + aryuZZ[104];
   aryuZZ[111]= - 3*aryuZZ[87];
   aryuZZ[115]= - 18*aryuZZ[86];
   aryuZZ[118]=4*aryuZZ[84];
   aryuZZ[121]= - 2*aryuZZ[96];
   aryuZZ[134]= - 18*aryuZZ[43];
   aryuZZ[136]= - 2*aryuZZ[47];
   aryuZZ[105]=aryuZZ[210] + aryuZZ[335] + 27./2.*aryuZZ[15] + 
   aryuZZ[136] + aryuZZ[134] + 33*aryuZZ[89] + aryuZZ[121] + 
   aryuZZ[118] + aryuZZ[115] + 1./2.*aryuZZ[105] + aryuZZ[111];
   aryuZZ[105]=aryuZZ[101]*aryuZZ[105];
   aryuZZ[139]=2*aryuZZ[95] + 1./2.*aryuZZ[83];
   aryuZZ[139]=11*aryuZZ[139] - 255./8.*aryuZZ[82];
   aryuZZ[142]=2*aryuZZ[78];
   aryuZZ[143]= - 44*aryuZZ[16] + 88*aryuZZ[89] + 2*aryuZZ[96] - 44*
   aryuZZ[90] + 41./4. + 44*aryuZZ[97];
   aryuZZ[143]=aryuZZ[100]*aryuZZ[143];
   aryuZZ[144]= - aryuZZ[15]*aryuZZ[100];
   aryuZZ[147]= - aryuZZ[13]*aryuZZ[100];
   aryuZZ[105]=aryuZZ[105] + 44*aryuZZ[147] + 2*aryuZZ[144] + 
   aryuZZ[143] + aryuZZ[142] + 22*aryuZZ[81] + 3*aryuZZ[139] + 55*
   aryuZZ[80];
   aryuZZ[105]=MMZ*aryuZZ[105];
   aryuZZ[139]= - 11./4.*aryuZZ[13] + aryuZZ[157] + aryuZZ[47] + 445./
   12. + aryuZZ[333];
   aryuZZ[139]=aryuZZ[17]*aryuZZ[139];
   aryuZZ[143]= - 10*aryuZZ[91];
   aryuZZ[144]=4*aryuZZ[93];
   aryuZZ[149]= - 43./12.*aryuZZ[75];
   aryuZZ[155]= - 6*aryuZZ[22];
   aryuZZ[139]=aryuZZ[139] - 32*aryuZZ[19] + 56./3.*aryuZZ[21] + 293./6.
   *aryuZZ[20] + aryuZZ[155] + aryuZZ[149] + aryuZZ[144] - aryuZZ[76]
    + aryuZZ[143] - 11./2.*aryuZZ[77] + aryuZZ[183] + aryuZZ[340];
   aryuZZ[139]=aryuZZ[101]*aryuZZ[139];
   aryuZZ[160]=1./2. + 6*aryuZZ[43];
   aryuZZ[167]=2*aryuZZ[47];
   aryuZZ[160]=aryuZZ[306] + 3*aryuZZ[160] + aryuZZ[167];
   aryuZZ[169]=aryuZZ[160] + aryuZZ[293];
   aryuZZ[169]=aryuZZ[101]*aryuZZ[169];
   aryuZZ[169]= - 2*aryuZZ[100] + aryuZZ[169];
   aryuZZ[169]=MMZ*aryuZZ[169];
   aryuZZ[170]=179./9. + aryuZZ[333];
   aryuZZ[171]= - 2*aryuZZ[45];
   aryuZZ[173]=1 + aryuZZ[114];
   aryuZZ[173]=2*aryuZZ[10]*aryuZZ[173];
   aryuZZ[169]=aryuZZ[173] + aryuZZ[169] + 1./2.*aryuZZ[222] - 
   aryuZZ[11] + aryuZZ[323] + aryuZZ[221] + 1./2.*aryuZZ[170] + 
   aryuZZ[171];
   aryuZZ[169]=aryuZZ[10]*aryuZZ[169];
   aryuZZ[170]=403./6. + aryuZZ[305];
   aryuZZ[176]=13./2.*aryuZZ[12];
   aryuZZ[178]=9*aryuZZ[10];
   aryuZZ[179]= - 3./2.*aryuZZ[11];
   aryuZZ[170]=aryuZZ[178] + aryuZZ[179] + 245./8.*aryuZZ[13] + 
   aryuZZ[176] + 1./2.*aryuZZ[170] + aryuZZ[224];
   aryuZZ[170]=aryuZZ[18]*aryuZZ[170];
   aryuZZ[180]=19 + aryuZZ[224];
   aryuZZ[180]=MMZ*aryuZZ[180];
   aryuZZ[181]=9./2.*aryuZZ[19];
   aryuZZ[161]=1./2.*aryuZZ[180] + aryuZZ[181] + aryuZZ[161];
   aryuZZ[161]=aryuZZ[10]*aryuZZ[161];
   aryuZZ[180]=19./4.*aryuZZ[12];
   aryuZZ[184]=377./4.*aryuZZ[13] + aryuZZ[180] + 392./3. + aryuZZ[306]
   ;
   aryuZZ[184]=aryuZZ[19]*aryuZZ[184];
   aryuZZ[188]=aryuZZ[323] - 11 + aryuZZ[300];
   aryuZZ[188]=aryuZZ[17]*aryuZZ[188];
   aryuZZ[189]=aryuZZ[335] - 1 + aryuZZ[305];
   aryuZZ[189]=aryuZZ[19]*aryuZZ[189];
   aryuZZ[192]= - 1 + aryuZZ[167];
   aryuZZ[192]=2*aryuZZ[192] - aryuZZ[45];
   aryuZZ[192]=MMZ*aryuZZ[192];
   aryuZZ[189]=aryuZZ[189] + aryuZZ[192];
   aryuZZ[189]=MMZ*aryuZZ[189];
   aryuZZ[192]=1 + aryuZZ[335];
   aryuZZ[194]=aryuZZ[18]*MMZ*aryuZZ[192];
   aryuZZ[189]=aryuZZ[189] + aryuZZ[194];
   aryuZZ[189]=aryuZZ[3]*aryuZZ[189];
   aryuZZ[196]= - 11*aryuZZ[77] - 3./2.*aryuZZ[75];
   aryuZZ[196]=7*aryuZZ[21] + 25./4.*aryuZZ[20] + 1./2.*aryuZZ[196] + 
   11*aryuZZ[22];
   aryuZZ[197]= - 4*aryuZZ[11] + 152*aryuZZ[13] + 151./4.*aryuZZ[12] + 
   aryuZZ[289] - 223./9. - aryuZZ[47];
   aryuZZ[197]=MMZ*aryuZZ[197];
   aryuZZ[170]=3*aryuZZ[189] + aryuZZ[170] + aryuZZ[161] + aryuZZ[197]
    + 3./2.*aryuZZ[188] + 3*aryuZZ[196] + aryuZZ[184];
   aryuZZ[170]=aryuZZ[3]*aryuZZ[170];
   aryuZZ[184]=aryuZZ[335] + aryuZZ[136] - 143./6. + aryuZZ[134];
   aryuZZ[188]=aryuZZ[184] + aryuZZ[210];
   aryuZZ[188]=aryuZZ[101]*aryuZZ[188];
   aryuZZ[189]= - aryuZZ[10]*aryuZZ[101];
   aryuZZ[196]=3./2.*aryuZZ[189];
   aryuZZ[188]=aryuZZ[196] - aryuZZ[100] + aryuZZ[188];
   aryuZZ[188]=aryuZZ[18]*aryuZZ[188];
   aryuZZ[197]= - 185773./324. + 173*aryuZZ[88];
   aryuZZ[200]=8./3.*aryuZZ[87];
   aryuZZ[202]= - 9./2.*aryuZZ[86];
   aryuZZ[203]= - 27./4.*aryuZZ[43];
   aryuZZ[205]=79./12.*aryuZZ[20] - 227./6.*aryuZZ[22] + 43./6.*
   aryuZZ[76] - 2*aryuZZ[93];
   aryuZZ[205]=aryuZZ[100]*aryuZZ[205];
   aryuZZ[206]=7./4.*aryuZZ[45];
   aryuZZ[207]= - 5619./4.*aryuZZ[13] - 115177./54. - 363*aryuZZ[12];
   aryuZZ[207]=aryuZZ[13]*aryuZZ[207];
   aryuZZ[209]=67./9. + aryuZZ[273];
   aryuZZ[209]=aryuZZ[11]*aryuZZ[209];
   aryuZZ[212]=56./3.*aryuZZ[100] + 55./2.*aryuZZ[147];
   aryuZZ[212]=aryuZZ[19]*aryuZZ[212];
   aryuZZ[219]=aryuZZ[13]*aryuZZ[100];
   aryuZZ[220]= - 251./3.*aryuZZ[100] + 55*aryuZZ[219];
   aryuZZ[220]=aryuZZ[17]*aryuZZ[220];
   aryuZZ[223]=aryuZZ[21]*aryuZZ[100];
   aryuZZ[105]=aryuZZ[170] + aryuZZ[188] + aryuZZ[169] + aryuZZ[105] + 
   aryuZZ[139] + 1./4.*aryuZZ[220] + aryuZZ[212] + 1./2.*aryuZZ[209] + 
   1./4.*aryuZZ[207] + 149./108.*aryuZZ[12] + aryuZZ[206] + 2*
   aryuZZ[223] + aryuZZ[307] + aryuZZ[156] + aryuZZ[205] + aryuZZ[203]
    - 1891./24.*aryuZZ[16] + 77./2.*aryuZZ[89] + 28./3.*aryuZZ[96] + 
   aryuZZ[202] + aryuZZ[200] - 51./2.*aryuZZ[90] - 22./3.*aryuZZ[97] + 
   5./8.*aryuZZ[197] - 11./3.*aryuZZ[85];
   aryuZZ[105]=aryuZZ[73]*aryuZZ[105];
   aryuZZ[139]= - 363./4.*aryuZZ[13] + aryuZZ[145] + aryuZZ[306] - 491./
   6. + aryuZZ[305];
   aryuZZ[139]=aryuZZ[19]*aryuZZ[139];
   aryuZZ[145]= - 99./2.*aryuZZ[13];
   aryuZZ[169]=aryuZZ[145] + aryuZZ[306] + 115./2. + aryuZZ[305];
   aryuZZ[169]=aryuZZ[17]*aryuZZ[169];
   aryuZZ[110]= - 1023./4.*aryuZZ[13] + aryuZZ[110] + aryuZZ[306] + 133.
   /4. + aryuZZ[305];
   aryuZZ[110]=MMZ*aryuZZ[110];
   aryuZZ[170]=15./4.*aryuZZ[10] + aryuZZ[179] - 231./4.*aryuZZ[13] + 
   165./8.*aryuZZ[12] + aryuZZ[306] - 475./24. + aryuZZ[305];
   aryuZZ[170]=aryuZZ[18]*aryuZZ[170];
   aryuZZ[188]=aryuZZ[300] + 13./2. + aryuZZ[47];
   aryuZZ[188]=MMZ*aryuZZ[188];
   aryuZZ[181]=3./2.*aryuZZ[188] + aryuZZ[181] + aryuZZ[17];
   aryuZZ[181]=aryuZZ[10]*aryuZZ[181];
   aryuZZ[188]=aryuZZ[289] + 1 - aryuZZ[47];
   aryuZZ[197]=aryuZZ[19]*aryuZZ[188];
   aryuZZ[188]=MMZ*aryuZZ[188];
   aryuZZ[197]=aryuZZ[197] + aryuZZ[188];
   aryuZZ[197]=MMZ*aryuZZ[197];
   aryuZZ[188]=aryuZZ[18]*aryuZZ[188];
   aryuZZ[188]=aryuZZ[197] + 1./2.*aryuZZ[188];
   aryuZZ[188]=aryuZZ[3]*aryuZZ[188];
   aryuZZ[197]=aryuZZ[152] + 11*aryuZZ[77] - aryuZZ[76];
   aryuZZ[197]= - 5./2.*aryuZZ[21] - aryuZZ[20] + 1./4.*aryuZZ[197] - 5
   *aryuZZ[22];
   aryuZZ[110]=9*aryuZZ[188] + aryuZZ[170] + aryuZZ[181] + 1./2.*
   aryuZZ[110] + 1./4.*aryuZZ[169] + 9*aryuZZ[197] + aryuZZ[139];
   aryuZZ[110]=aryuZZ[3]*aryuZZ[110];
   aryuZZ[139]= - 555./4. - 11*aryuZZ[85];
   aryuZZ[139]= - aryuZZ[87] + 1./2.*aryuZZ[139] - 11*aryuZZ[90];
   aryuZZ[139]=1./2.*aryuZZ[139] + aryuZZ[297];
   aryuZZ[169]=2*aryuZZ[84];
   aryuZZ[170]= - 3./2.*aryuZZ[47];
   aryuZZ[181]= - 33./4.*aryuZZ[13];
   aryuZZ[139]=aryuZZ[181] + aryuZZ[298] + 73./4.*aryuZZ[15] + 
   aryuZZ[170] + aryuZZ[285] - 99./2.*aryuZZ[89] - 3./2.*aryuZZ[96] + 3
   *aryuZZ[139] + aryuZZ[169];
   aryuZZ[139]=aryuZZ[101]*aryuZZ[139];
   aryuZZ[188]= - 11./4.*aryuZZ[81] + 153./16.*aryuZZ[82] - 11*
   aryuZZ[80];
   aryuZZ[197]= - aryuZZ[96] + 11*aryuZZ[90] + 5./2. - 11*aryuZZ[97];
   aryuZZ[197]=11./2.*aryuZZ[16] + 1./2.*aryuZZ[197] - 11*aryuZZ[89];
   aryuZZ[197]=aryuZZ[100]*aryuZZ[197];
   aryuZZ[205]=aryuZZ[15]*aryuZZ[100];
   aryuZZ[139]=aryuZZ[139] + 33./2.*aryuZZ[219] + 3./2.*aryuZZ[205] + 3
   *aryuZZ[197] + 3*aryuZZ[188] + aryuZZ[78];
   aryuZZ[139]=MMZ*aryuZZ[139];
   aryuZZ[188]=1./4. + 3*aryuZZ[43];
   aryuZZ[197]=11./4.*aryuZZ[13] + aryuZZ[251] + aryuZZ[188] + 
   aryuZZ[236];
   aryuZZ[197]=aryuZZ[101]*aryuZZ[197];
   aryuZZ[197]=1./2.*aryuZZ[100] + aryuZZ[197];
   aryuZZ[197]=MMZ*aryuZZ[197];
   aryuZZ[207]= - 625./9. + aryuZZ[333];
   aryuZZ[209]=1./2. + aryuZZ[114];
   aryuZZ[209]=aryuZZ[10]*aryuZZ[209];
   aryuZZ[175]=aryuZZ[209] + 3*aryuZZ[197] + aryuZZ[208] - aryuZZ[11]
    - 121./4.*aryuZZ[13] + aryuZZ[175] - aryuZZ[45] + 1./4.*aryuZZ[207]
    + aryuZZ[136];
   aryuZZ[175]=aryuZZ[10]*aryuZZ[175];
   aryuZZ[197]=33./4.*aryuZZ[13] + aryuZZ[157] + 3./2.*aryuZZ[47] + 
   1105./12. + aryuZZ[333];
   aryuZZ[197]=aryuZZ[17]*aryuZZ[197];
   aryuZZ[207]=11./2.*aryuZZ[77] + aryuZZ[342] - 11*aryuZZ[92];
   aryuZZ[209]= - 5*aryuZZ[91];
   aryuZZ[212]=2*aryuZZ[93];
   aryuZZ[220]= - 43./24.*aryuZZ[75];
   aryuZZ[197]=1./2.*aryuZZ[197] + 57./2.*aryuZZ[19] - 5./3.*aryuZZ[21]
    + 181./6.*aryuZZ[20] - 5./2.*aryuZZ[22] + aryuZZ[220] + aryuZZ[212]
    - 3./4.*aryuZZ[76] + 3./2.*aryuZZ[207] + aryuZZ[209];
   aryuZZ[197]=aryuZZ[101]*aryuZZ[197];
   aryuZZ[170]=aryuZZ[181] + aryuZZ[298] + aryuZZ[170] - 149./12. + 
   aryuZZ[285];
   aryuZZ[170]=aryuZZ[101]*aryuZZ[170];
   aryuZZ[181]=3./4.*aryuZZ[189];
   aryuZZ[170]=aryuZZ[181] + aryuZZ[100] + aryuZZ[170];
   aryuZZ[170]=aryuZZ[18]*aryuZZ[170];
   aryuZZ[189]= - 15821./162. - 435*aryuZZ[88];
   aryuZZ[104]=1./2.*aryuZZ[189] + aryuZZ[104];
   aryuZZ[104]=1./2.*aryuZZ[104] + 11*aryuZZ[97];
   aryuZZ[189]= - 21./4.*aryuZZ[20] + 37./2.*aryuZZ[22] - 3*aryuZZ[76]
    + aryuZZ[212];
   aryuZZ[189]=aryuZZ[100]*aryuZZ[189];
   aryuZZ[207]= - aryuZZ[21]*aryuZZ[100];
   aryuZZ[225]=1797*aryuZZ[13] + 16301./9. - 1089./2.*aryuZZ[12];
   aryuZZ[225]=aryuZZ[13]*aryuZZ[225];
   aryuZZ[210]=95./9. + aryuZZ[210];
   aryuZZ[210]=aryuZZ[11]*aryuZZ[210];
   aryuZZ[226]= - 5*aryuZZ[100] + 11./2.*aryuZZ[219];
   aryuZZ[226]=aryuZZ[19]*aryuZZ[226];
   aryuZZ[228]=15*aryuZZ[100] + 11*aryuZZ[147];
   aryuZZ[228]=aryuZZ[17]*aryuZZ[228];
   aryuZZ[104]=aryuZZ[110] + aryuZZ[170] + aryuZZ[175] + aryuZZ[139] + 
   aryuZZ[197] + 3./4.*aryuZZ[228] + 3*aryuZZ[226] + 1./2.*aryuZZ[210]
    + 1./8.*aryuZZ[225] - 5027./72.*aryuZZ[12] + 5./8.*aryuZZ[45] + 2*
   aryuZZ[207] + 7./6.*aryuZZ[15] + 5./4.*aryuZZ[47] + aryuZZ[189] - 45.
   /8.*aryuZZ[43] + 391./8.*aryuZZ[16] - 33*aryuZZ[89] - 3*aryuZZ[96]
    + aryuZZ[140] + 1./2.*aryuZZ[104] + 4./3.*aryuZZ[87];
   aryuZZ[104]=aryuZZ[73]*aryuZZ[104];
   aryuZZ[110]=aryuZZ[119] + aryuZZ[21] - aryuZZ[8] + aryuZZ[198];
   aryuZZ[110]=aryuZZ[101]*aryuZZ[110];
   aryuZZ[139]= - 3*aryuZZ[24] - aryuZZ[26];
   aryuZZ[139]=2*aryuZZ[139] + aryuZZ[25];
   aryuZZ[170]= - aryuZZ[15] - 3*aryuZZ[14] + aryuZZ[28] - 3 - 2*
   aryuZZ[30];
   aryuZZ[170]=aryuZZ[101]*aryuZZ[170];
   aryuZZ[139]=2*aryuZZ[139] + aryuZZ[170];
   aryuZZ[139]=MMZ*aryuZZ[139];
   aryuZZ[175]=1./9. + aryuZZ[177];
   aryuZZ[187]=aryuZZ[175] + aryuZZ[187];
   aryuZZ[189]=aryuZZ[1]*aryuZZ[187];
   aryuZZ[187]=aryuZZ[39]*aryuZZ[187];
   aryuZZ[197]= - 37*aryuZZ[32] - 3989./432. + 33*aryuZZ[29];
   aryuZZ[198]= - 7./3. + aryuZZ[7];
   aryuZZ[198]=aryuZZ[13]*aryuZZ[198];
   aryuZZ[139]=aryuZZ[187] + 2./3.*aryuZZ[189] + aryuZZ[139] + 
   aryuZZ[110] + aryuZZ[216] + 11./4.*aryuZZ[198] + 113./18.*aryuZZ[7]
    + aryuZZ[215] + aryuZZ[15] - 5./4.*aryuZZ[177] - 33./2.*aryuZZ[16]
    + aryuZZ[174] - aryuZZ[28] + 1./2.*aryuZZ[197] + aryuZZ[172];
   aryuZZ[139]=aryuZZ[39]*aryuZZ[139];
   aryuZZ[113]=aryuZZ[243] - aryuZZ[14] + aryuZZ[125] - 1 + aryuZZ[113]
   ;
   aryuZZ[113]=aryuZZ[101]*aryuZZ[113];
   aryuZZ[125]= - aryuZZ[24] - 1./3.*aryuZZ[26];
   aryuZZ[125]=2*aryuZZ[125] + 1./3.*aryuZZ[25];
   aryuZZ[125]=2*aryuZZ[125] + aryuZZ[113];
   aryuZZ[125]=MMZ*aryuZZ[125];
   aryuZZ[172]= - 37./3.*aryuZZ[32] - 3989./1296. + 11*aryuZZ[29];
   aryuZZ[125]=1./9.*aryuZZ[189] + aryuZZ[125] + 1./3.*aryuZZ[110] + 
   aryuZZ[331] + 11./12.*aryuZZ[198] + 113./54.*aryuZZ[7] - 55./36.*
   aryuZZ[12] + aryuZZ[227] + aryuZZ[259] - 11./2.*aryuZZ[16] + 
   aryuZZ[257] + aryuZZ[255] + 1./2.*aryuZZ[172] - 5./3.*aryuZZ[31];
   aryuZZ[125]=aryuZZ[1]*aryuZZ[125];
   aryuZZ[172]=99./2.*aryuZZ[13] - 17 - 99./2.*aryuZZ[12];
   aryuZZ[187]=3*aryuZZ[11];
   aryuZZ[122]=aryuZZ[122] + 1./2.*aryuZZ[172] + aryuZZ[187];
   aryuZZ[122]=aryuZZ[18]*aryuZZ[122];
   aryuZZ[172]=99*aryuZZ[13] + 17 - 99*aryuZZ[12];
   aryuZZ[172]=1./4.*aryuZZ[172] + aryuZZ[187];
   aryuZZ[172]=aryuZZ[19]*aryuZZ[172];
   aryuZZ[186]=MMZ*aryuZZ[186];
   aryuZZ[187]=aryuZZ[10]*aryuZZ[120];
   aryuZZ[122]=1./2.*aryuZZ[122] + 3*aryuZZ[187] + aryuZZ[172] + 3*
   aryuZZ[186];
   aryuZZ[122]=aryuZZ[3]*aryuZZ[122];
   aryuZZ[172]= - 1./2.*aryuZZ[93] + aryuZZ[22];
   aryuZZ[172]=aryuZZ[100]*aryuZZ[172];
   aryuZZ[172]=aryuZZ[172] + 1./2.*aryuZZ[223];
   aryuZZ[135]=aryuZZ[135] - 1./2.*aryuZZ[21] + aryuZZ[301] - 
   aryuZZ[22];
   aryuZZ[135]=aryuZZ[101]*aryuZZ[135];
   aryuZZ[145]=aryuZZ[145] - 59./3. + 99./2.*aryuZZ[12];
   aryuZZ[145]=aryuZZ[13]*aryuZZ[145];
   aryuZZ[186]= - 3./2.*aryuZZ[13];
   aryuZZ[187]= - 7./9. + aryuZZ[186];
   aryuZZ[187]=aryuZZ[11]*aryuZZ[187];
   aryuZZ[189]= - aryuZZ[100] + aryuZZ[101];
   aryuZZ[197]=MMZ*aryuZZ[189];
   aryuZZ[198]=3*aryuZZ[13];
   aryuZZ[207]=aryuZZ[198] + 7./9. - 3./2.*aryuZZ[12];
   aryuZZ[207]= - aryuZZ[10] + 11./2.*aryuZZ[207] + aryuZZ[11];
   aryuZZ[207]=aryuZZ[10]*aryuZZ[207];
   aryuZZ[189]=aryuZZ[18]*aryuZZ[189];
   aryuZZ[122]=aryuZZ[122] + 3./4.*aryuZZ[189] + aryuZZ[207] + 9./16.*
   aryuZZ[197] + 3*aryuZZ[135] + 9./2.*aryuZZ[222] + 11./2.*aryuZZ[187]
    + 11./8.*aryuZZ[145] + 3*aryuZZ[172] + 649./24.*aryuZZ[12];
   aryuZZ[122]=aryuZZ[73]*aryuZZ[122];
   aryuZZ[135]=aryuZZ[1]*aryuZZ[120];
   aryuZZ[120]=aryuZZ[39]*aryuZZ[120];
   aryuZZ[124]=aryuZZ[18]*aryuZZ[124];
   aryuZZ[120]=1./2.*aryuZZ[124] + aryuZZ[135] + 3*aryuZZ[120];
   aryuZZ[120]=aryuZZ[6]*aryuZZ[120];
   aryuZZ[124]=aryuZZ[19]*aryuZZ[7];
   aryuZZ[135]=MMZ*aryuZZ[7];
   aryuZZ[124]=aryuZZ[124] + aryuZZ[135];
   aryuZZ[135]=aryuZZ[1]*aryuZZ[124];
   aryuZZ[124]=aryuZZ[39]*aryuZZ[124];
   aryuZZ[145]=aryuZZ[182] + 3*aryuZZ[193];
   aryuZZ[145]=aryuZZ[18]*aryuZZ[145];
   aryuZZ[120]=aryuZZ[120] + 1./2.*aryuZZ[145] + aryuZZ[135] + 3*
   aryuZZ[124];
   aryuZZ[120]=aryuZZ[3]*aryuZZ[120];
   aryuZZ[124]=11./2.*aryuZZ[12] - 59./9.*aryuZZ[7];
   aryuZZ[135]= - 3*aryuZZ[7];
   aryuZZ[145]= - 1 + aryuZZ[135];
   aryuZZ[145]=aryuZZ[13]*aryuZZ[145];
   aryuZZ[145]=aryuZZ[124] + 11./2.*aryuZZ[145];
   aryuZZ[145]=1./3.*aryuZZ[158] + 2./9.*aryuZZ[154] + 1./2.*
   aryuZZ[145] + aryuZZ[151];
   aryuZZ[145]=aryuZZ[39]*aryuZZ[145];
   aryuZZ[158]=aryuZZ[323] + 59./27. + aryuZZ[319];
   aryuZZ[158]=1./9.*aryuZZ[320] + 1./2.*aryuZZ[158] + aryuZZ[216];
   aryuZZ[158]=aryuZZ[1]*aryuZZ[158];
   aryuZZ[130]=aryuZZ[141] + aryuZZ[130] + 2./3.*aryuZZ[320];
   aryuZZ[130]=aryuZZ[39]*aryuZZ[130];
   aryuZZ[141]= - 2./3.*aryuZZ[1];
   aryuZZ[172]=aryuZZ[141] - aryuZZ[39];
   aryuZZ[172]=aryuZZ[39]*aryuZZ[172];
   aryuZZ[182]=pow(aryuZZ[1],2);
   aryuZZ[172]= - 1./9.*aryuZZ[182] + aryuZZ[172];
   aryuZZ[172]=aryuZZ[6]*aryuZZ[172];
   aryuZZ[130]=aryuZZ[172] + 2*aryuZZ[322] + aryuZZ[158] + aryuZZ[130];
   aryuZZ[130]=aryuZZ[6]*aryuZZ[130];
   aryuZZ[158]=aryuZZ[13]*aryuZZ[263];
   aryuZZ[124]=1./3.*aryuZZ[124] + 11./2.*aryuZZ[158];
   aryuZZ[124]=1./27.*aryuZZ[154] + 1./2.*aryuZZ[124] + aryuZZ[138];
   aryuZZ[124]=aryuZZ[1]*aryuZZ[124];
   aryuZZ[120]=aryuZZ[122] + aryuZZ[120] + aryuZZ[130] + aryuZZ[321] + 
   aryuZZ[124] + aryuZZ[145];
   aryuZZ[120]=aryuZZ[133]*aryuZZ[120];
   aryuZZ[122]= - 41./3. + aryuZZ[163];
   aryuZZ[117]=5./6.*aryuZZ[122] + aryuZZ[117];
   aryuZZ[124]= - 1 - aryuZZ[7];
   aryuZZ[130]=aryuZZ[1]*aryuZZ[124];
   aryuZZ[117]=1./9.*aryuZZ[130] + 1./3.*aryuZZ[250] + 1./6.*
   aryuZZ[116] + 1./2.*aryuZZ[117] + aryuZZ[151];
   aryuZZ[117]=aryuZZ[1]*aryuZZ[117];
   aryuZZ[122]=5./2.*aryuZZ[122] - 99*aryuZZ[13];
   aryuZZ[124]=aryuZZ[39]*aryuZZ[124];
   aryuZZ[122]=aryuZZ[124] + 2./3.*aryuZZ[130] + aryuZZ[250] + 
   aryuZZ[185] + 1./2.*aryuZZ[122] - aryuZZ[11];
   aryuZZ[122]=aryuZZ[39]*aryuZZ[122];
   aryuZZ[124]=aryuZZ[1]*aryuZZ[217];
   aryuZZ[130]=3 + aryuZZ[114];
   aryuZZ[138]=aryuZZ[39]*aryuZZ[130];
   aryuZZ[124]=aryuZZ[124] + aryuZZ[138];
   aryuZZ[124]=aryuZZ[10]*aryuZZ[124];
   aryuZZ[138]=aryuZZ[1]*aryuZZ[101];
   aryuZZ[145]=aryuZZ[39]*aryuZZ[101];
   aryuZZ[151]=1./3.*aryuZZ[138] + aryuZZ[145];
   aryuZZ[151]=aryuZZ[18]*aryuZZ[151];
   aryuZZ[154]=5./4. + 2./3.*aryuZZ[1];
   aryuZZ[154]=aryuZZ[1]*aryuZZ[154];
   aryuZZ[158]=4./3.*aryuZZ[1];
   aryuZZ[163]=2*aryuZZ[39] + 5./4. + aryuZZ[158];
   aryuZZ[163]=aryuZZ[39]*aryuZZ[163];
   aryuZZ[154]=1./3.*aryuZZ[154] + aryuZZ[163];
   aryuZZ[154]=aryuZZ[6]*aryuZZ[154];
   aryuZZ[117]=aryuZZ[154] + aryuZZ[151] + aryuZZ[124] + aryuZZ[117] + 
   aryuZZ[122];
   aryuZZ[117]=aryuZZ[6]*aryuZZ[117];
   aryuZZ[106]=aryuZZ[106] + aryuZZ[214] - aryuZZ[8];
   aryuZZ[122]=aryuZZ[128] + aryuZZ[131] + aryuZZ[132] + aryuZZ[127] + 
   3*aryuZZ[106] + aryuZZ[126];
   aryuZZ[122]=aryuZZ[39]*aryuZZ[122];
   aryuZZ[124]=1./6.*aryuZZ[20];
   aryuZZ[128]= - 1./6.*aryuZZ[17];
   aryuZZ[132]= - 4*MMZ;
   aryuZZ[106]=aryuZZ[132] + aryuZZ[128] + 2./3.*aryuZZ[19] + 
   aryuZZ[21] + aryuZZ[106] + aryuZZ[124];
   aryuZZ[106]=aryuZZ[1]*aryuZZ[106];
   aryuZZ[151]=9./2.*MMZ + aryuZZ[310] + 1./6.*aryuZZ[17];
   aryuZZ[151]=aryuZZ[1]*aryuZZ[151];
   aryuZZ[119]=27./2.*MMZ + 9*aryuZZ[19] + aryuZZ[119];
   aryuZZ[119]=aryuZZ[39]*aryuZZ[119];
   aryuZZ[154]=aryuZZ[18]*aryuZZ[159];
   aryuZZ[119]=2*aryuZZ[154] + aryuZZ[151] + aryuZZ[119];
   aryuZZ[119]=aryuZZ[6]*aryuZZ[119];
   aryuZZ[151]= - 3./2.*aryuZZ[7];
   aryuZZ[154]=1 + aryuZZ[151];
   aryuZZ[154]=aryuZZ[39]*aryuZZ[154];
   aryuZZ[154]=aryuZZ[315] + aryuZZ[154];
   aryuZZ[154]=aryuZZ[18]*aryuZZ[154];
   aryuZZ[106]=aryuZZ[119] + aryuZZ[154] + aryuZZ[106] + aryuZZ[122];
   aryuZZ[106]=aryuZZ[3]*aryuZZ[106];
   aryuZZ[119]=2./3. - aryuZZ[7];
   aryuZZ[122]=aryuZZ[1]*aryuZZ[119];
   aryuZZ[119]=aryuZZ[39]*aryuZZ[119];
   aryuZZ[119]=1./3.*aryuZZ[122] + aryuZZ[119];
   aryuZZ[119]=aryuZZ[10]*aryuZZ[119];
   aryuZZ[122]= - 2./3.*aryuZZ[78];
   aryuZZ[154]=aryuZZ[122] + 13./4.*aryuZZ[81] + 3./2.*aryuZZ[79] + 13*
   aryuZZ[80];
   aryuZZ[154]=MMH*aryuZZ[73]*aryuZZ[154];
   aryuZZ[104]=aryuZZ[120] + aryuZZ[154] + aryuZZ[104] + aryuZZ[106] + 
   aryuZZ[117] + aryuZZ[119] + aryuZZ[125] + aryuZZ[139];
   aryuZZ[104]=aryuZZ[133]*aryuZZ[104];
   aryuZZ[106]=28*aryuZZ[24];
   aryuZZ[117]=8*aryuZZ[26];
   aryuZZ[119]= - 4*aryuZZ[25];
   aryuZZ[120]=aryuZZ[119] + aryuZZ[117] + aryuZZ[106] + 5./9.*
   aryuZZ[27];
   aryuZZ[120]=MMZ*aryuZZ[120];
   aryuZZ[125]= - 28*aryuZZ[29];
   aryuZZ[139]=37./4.*aryuZZ[23];
   aryuZZ[154]=83./9. + 9*aryuZZ[7];
   aryuZZ[154]=aryuZZ[13]*aryuZZ[154];
   aryuZZ[159]=19./27.*aryuZZ[7] - 37./81. - aryuZZ[177];
   aryuZZ[159]=aryuZZ[1]*aryuZZ[159];
   aryuZZ[163]=11./27.*aryuZZ[7] - 29./81. - aryuZZ[177];
   aryuZZ[163]=aryuZZ[39]*aryuZZ[163];
   aryuZZ[120]=aryuZZ[163] + 2./3.*aryuZZ[159] + aryuZZ[120] + 
   aryuZZ[271] + 5./4.*aryuZZ[154] + 61./18.*aryuZZ[7] + aryuZZ[270] + 
   aryuZZ[177] + aryuZZ[258] + 10*aryuZZ[31] + 81./2.*aryuZZ[32] + 
   aryuZZ[139] + 24539./1944. + aryuZZ[125];
   aryuZZ[120]=aryuZZ[39]*aryuZZ[120];
   aryuZZ[106]=aryuZZ[119] + aryuZZ[117] + aryuZZ[106] + aryuZZ[27];
   aryuZZ[106]=MMZ*aryuZZ[106];
   aryuZZ[117]=aryuZZ[139] + 893./72. + aryuZZ[125];
   aryuZZ[119]=83./27. + aryuZZ[150];
   aryuZZ[119]=aryuZZ[13]*aryuZZ[119];
   aryuZZ[125]=aryuZZ[7] - 5./9. - aryuZZ[177];
   aryuZZ[125]=aryuZZ[1]*aryuZZ[125];
   aryuZZ[106]=1./9.*aryuZZ[125] + 1./3.*aryuZZ[106] + aryuZZ[216] + 5./
   4.*aryuZZ[119] + 61./54.*aryuZZ[7] - 125./54.*aryuZZ[12] + 1./3.*
   aryuZZ[177] + 28./3.*aryuZZ[16] + 10./3.*aryuZZ[31] + 1./3.*
   aryuZZ[117] + 27./2.*aryuZZ[32];
   aryuZZ[106]=aryuZZ[1]*aryuZZ[106];
   aryuZZ[117]= - 5*aryuZZ[19] - 19./3.*MMZ;
   aryuZZ[117]=aryuZZ[1]*aryuZZ[117];
   aryuZZ[119]= - 29*aryuZZ[19] - 41*MMZ;
   aryuZZ[119]=aryuZZ[39]*aryuZZ[119];
   aryuZZ[125]=1./2.*aryuZZ[18]*aryuZZ[375];
   aryuZZ[117]=aryuZZ[125] + aryuZZ[117] + 1./3.*aryuZZ[119];
   aryuZZ[117]=aryuZZ[6]*aryuZZ[117];
   aryuZZ[119]= - 2*aryuZZ[7];
   aryuZZ[139]=19./3. + aryuZZ[119];
   aryuZZ[139]=MMZ*aryuZZ[139];
   aryuZZ[109]=aryuZZ[109] + 2*aryuZZ[139];
   aryuZZ[109]=aryuZZ[1]*aryuZZ[109];
   aryuZZ[119]=49./9. + aryuZZ[119];
   aryuZZ[119]=MMZ*aryuZZ[119];
   aryuZZ[119]=aryuZZ[211] + 2*aryuZZ[119];
   aryuZZ[119]=aryuZZ[39]*aryuZZ[119];
   aryuZZ[139]=aryuZZ[1]*aryuZZ[234];
   aryuZZ[150]=10./9. + aryuZZ[151];
   aryuZZ[150]=aryuZZ[39]*aryuZZ[150];
   aryuZZ[139]=aryuZZ[139] + aryuZZ[150];
   aryuZZ[139]=aryuZZ[18]*aryuZZ[139];
   aryuZZ[109]=aryuZZ[117] + aryuZZ[139] + 1./3.*aryuZZ[109] + 
   aryuZZ[119];
   aryuZZ[109]=aryuZZ[3]*aryuZZ[109];
   aryuZZ[117]=194*aryuZZ[13];
   aryuZZ[119]= - 19*aryuZZ[7];
   aryuZZ[150]=55./3. + aryuZZ[119];
   aryuZZ[150]=aryuZZ[1]*aryuZZ[150];
   aryuZZ[137]=47./3. + aryuZZ[137];
   aryuZZ[137]=aryuZZ[39]*aryuZZ[137];
   aryuZZ[137]=1./3.*aryuZZ[137] + 2./9.*aryuZZ[150] - 11./3.*
   aryuZZ[11] + aryuZZ[117] + 2615./27. + 37./2.*aryuZZ[12];
   aryuZZ[137]=aryuZZ[39]*aryuZZ[137];
   aryuZZ[117]=aryuZZ[117] + 349./3. + 125./2.*aryuZZ[12];
   aryuZZ[150]=7./9. - aryuZZ[7];
   aryuZZ[150]=aryuZZ[1]*aryuZZ[150];
   aryuZZ[117]=1./3.*aryuZZ[150] + 1./9.*aryuZZ[117] - aryuZZ[11];
   aryuZZ[117]=aryuZZ[1]*aryuZZ[117];
   aryuZZ[150]= - 1 + aryuZZ[141];
   aryuZZ[150]=aryuZZ[1]*aryuZZ[150];
   aryuZZ[151]= - 2*aryuZZ[39] - 1 + aryuZZ[201];
   aryuZZ[151]=aryuZZ[39]*aryuZZ[151];
   aryuZZ[150]=1./3.*aryuZZ[150] + aryuZZ[151];
   aryuZZ[150]=aryuZZ[6]*aryuZZ[150];
   aryuZZ[117]=aryuZZ[150] + aryuZZ[117] + 1./3.*aryuZZ[137];
   aryuZZ[117]=aryuZZ[6]*aryuZZ[117];
   aryuZZ[137]=3*aryuZZ[79];
   aryuZZ[150]= - 4./3.*aryuZZ[78];
   aryuZZ[151]=aryuZZ[150] - 67./12.*aryuZZ[81] + aryuZZ[137] - 89./6.*
   aryuZZ[80];
   aryuZZ[151]=MMH*aryuZZ[73]*aryuZZ[151];
   aryuZZ[104]=aryuZZ[104] + aryuZZ[151] + aryuZZ[105] + aryuZZ[109] + 
   aryuZZ[117] + aryuZZ[325] + aryuZZ[106] + aryuZZ[120];
   aryuZZ[104]=aryuZZ[133]*aryuZZ[104];
   aryuZZ[105]=1./3.*aryuZZ[90] - 2555./12. + aryuZZ[148];
   aryuZZ[106]=3./2.*aryuZZ[13];
   aryuZZ[105]=aryuZZ[106] + aryuZZ[335] + 35./2.*aryuZZ[15] + 
   aryuZZ[136] + aryuZZ[134] + 9*aryuZZ[89] + aryuZZ[121] + aryuZZ[118]
    + aryuZZ[115] + 1./2.*aryuZZ[105] + aryuZZ[111];
   aryuZZ[105]=aryuZZ[101]*aryuZZ[105];
   aryuZZ[109]=aryuZZ[83] - 23./4.*aryuZZ[82];
   aryuZZ[109]=1./2.*aryuZZ[109] + aryuZZ[80];
   aryuZZ[105]=aryuZZ[105] + 1./3.*aryuZZ[109] + aryuZZ[142];
   aryuZZ[105]=MMZ*aryuZZ[105];
   aryuZZ[109]=aryuZZ[160] + aryuZZ[186];
   aryuZZ[109]=MMZ*aryuZZ[101]*aryuZZ[109];
   aryuZZ[109]=aryuZZ[173] + aryuZZ[109] + aryuZZ[208] - aryuZZ[11] + 
   aryuZZ[198] + aryuZZ[221] + aryuZZ[171] + 40./9. + aryuZZ[278];
   aryuZZ[109]=aryuZZ[10]*aryuZZ[109];
   aryuZZ[115]=187./6. + aryuZZ[305];
   aryuZZ[115]=aryuZZ[178] + aryuZZ[179] + 53./8.*aryuZZ[13] + 
   aryuZZ[176] + 1./2.*aryuZZ[115] + aryuZZ[224];
   aryuZZ[115]=aryuZZ[18]*aryuZZ[115];
   aryuZZ[117]=1./4.*aryuZZ[13];
   aryuZZ[118]=aryuZZ[117] + aryuZZ[180] + 8 + aryuZZ[306];
   aryuZZ[118]=aryuZZ[19]*aryuZZ[118];
   aryuZZ[120]= - aryuZZ[11] + 20./3.*aryuZZ[13] + aryuZZ[180] + 
   aryuZZ[289] - 25./6. + aryuZZ[47];
   aryuZZ[120]=MMZ*aryuZZ[120];
   aryuZZ[133]=aryuZZ[19]*aryuZZ[192];
   aryuZZ[134]=2 + aryuZZ[379];
   aryuZZ[134]=MMZ*aryuZZ[134];
   aryuZZ[133]=3*aryuZZ[133] + aryuZZ[134];
   aryuZZ[133]=MMZ*aryuZZ[133];
   aryuZZ[133]=aryuZZ[133] + 3*aryuZZ[194];
   aryuZZ[133]=aryuZZ[3]*aryuZZ[133];
   aryuZZ[134]= - aryuZZ[77] + aryuZZ[152];
   aryuZZ[134]=1./2.*aryuZZ[134] + aryuZZ[22];
   aryuZZ[142]=5./3.*aryuZZ[13] - 5./3. + aryuZZ[306];
   aryuZZ[142]=aryuZZ[17]*aryuZZ[142];
   aryuZZ[115]=aryuZZ[133] + aryuZZ[115] + aryuZZ[161] + aryuZZ[120] + 
   1./2.*aryuZZ[142] + aryuZZ[118] + 9*aryuZZ[21] + 9*aryuZZ[134] + 37./
   12.*aryuZZ[20];
   aryuZZ[115]=aryuZZ[3]*aryuZZ[115];
   aryuZZ[118]= - 3./4.*aryuZZ[13] + aryuZZ[157] + aryuZZ[47] + 565./12.
    + aryuZZ[333];
   aryuZZ[118]=aryuZZ[17]*aryuZZ[118];
   aryuZZ[120]= - 1./2.*aryuZZ[77] + aryuZZ[342] + aryuZZ[92];
   aryuZZ[118]=aryuZZ[118] - 16*aryuZZ[19] + 44./3.*aryuZZ[21] + 305./6.
   *aryuZZ[20] + aryuZZ[155] + aryuZZ[149] + aryuZZ[144] - aryuZZ[76]
    + 3*aryuZZ[120] + aryuZZ[143];
   aryuZZ[118]=aryuZZ[101]*aryuZZ[118];
   aryuZZ[120]=aryuZZ[184] + aryuZZ[106];
   aryuZZ[120]=aryuZZ[101]*aryuZZ[120];
   aryuZZ[120]=aryuZZ[120] + aryuZZ[196];
   aryuZZ[120]=aryuZZ[18]*aryuZZ[120];
   aryuZZ[133]= - 75473./324. + 33*aryuZZ[88];
   aryuZZ[134]= - aryuZZ[22] + aryuZZ[126];
   aryuZZ[134]=aryuZZ[100]*aryuZZ[134];
   aryuZZ[142]=29./2. - 19*aryuZZ[12];
   aryuZZ[142]=1./9.*aryuZZ[142] + 13./4.*aryuZZ[13];
   aryuZZ[142]=aryuZZ[13]*aryuZZ[142];
   aryuZZ[143]=7./3. + aryuZZ[240];
   aryuZZ[143]=aryuZZ[11]*aryuZZ[143];
   aryuZZ[144]=5*aryuZZ[100] + 1./2.*aryuZZ[147];
   aryuZZ[144]=aryuZZ[19]*aryuZZ[144];
   aryuZZ[148]= - aryuZZ[100] + aryuZZ[219];
   aryuZZ[148]=aryuZZ[17]*aryuZZ[148];
   aryuZZ[105]=aryuZZ[115] + aryuZZ[120] + aryuZZ[109] + aryuZZ[105] + 
   aryuZZ[118] + 1./12.*aryuZZ[148] + 1./3.*aryuZZ[144] + 1./6.*
   aryuZZ[143] + 1./4.*aryuZZ[142] + 77./108.*aryuZZ[12] + aryuZZ[206]
    + 8./3.*aryuZZ[15] + aryuZZ[156] + 1./6.*aryuZZ[134] + aryuZZ[203]
    - 33./8.*aryuZZ[16] + 19./6.*aryuZZ[89] + 5./6.*aryuZZ[96] + 
   aryuZZ[202] + aryuZZ[200] - 5./18.*aryuZZ[90] + 1./8.*aryuZZ[133] - 
   5./3.*aryuZZ[85];
   aryuZZ[105]=aryuZZ[73]*aryuZZ[105];
   aryuZZ[109]= - 877./4. + 1./3.*aryuZZ[85];
   aryuZZ[109]=1./2.*aryuZZ[109] + aryuZZ[111];
   aryuZZ[111]=1./12.*aryuZZ[13];
   aryuZZ[109]=aryuZZ[111] + aryuZZ[298] + 125./12.*aryuZZ[15] + 
   aryuZZ[136] + aryuZZ[285] + 1./2.*aryuZZ[89] + aryuZZ[121] + 
   aryuZZ[169] + 1./2.*aryuZZ[109] - 9*aryuZZ[86];
   aryuZZ[109]=aryuZZ[101]*aryuZZ[109];
   aryuZZ[109]=aryuZZ[109] + 3./16.*aryuZZ[82] + aryuZZ[78];
   aryuZZ[109]=MMZ*aryuZZ[109];
   aryuZZ[115]= - aryuZZ[77] - 9./2.*aryuZZ[75];
   aryuZZ[115]=aryuZZ[123] + 7./3.*aryuZZ[20] + 1./2.*aryuZZ[115] + 
   aryuZZ[22];
   aryuZZ[118]=pow(MMZ,2);
   aryuZZ[120]=aryuZZ[118]*aryuZZ[192];
   aryuZZ[120]=aryuZZ[120] + 3./2.*aryuZZ[194];
   aryuZZ[120]=aryuZZ[3]*aryuZZ[120];
   aryuZZ[123]= - 8./3. + 1./4.*aryuZZ[12];
   aryuZZ[123]=aryuZZ[19]*aryuZZ[123];
   aryuZZ[133]=1./3.*aryuZZ[13] - 1./3. + aryuZZ[224];
   aryuZZ[133]=aryuZZ[17]*aryuZZ[133];
   aryuZZ[106]=aryuZZ[106] - 17./3.*aryuZZ[12] - 21./2. - aryuZZ[45];
   aryuZZ[106]=MMZ*aryuZZ[106];
   aryuZZ[134]=7 + aryuZZ[45];
   aryuZZ[134]=MMZ*aryuZZ[134];
   aryuZZ[134]=aryuZZ[17] + 3./4.*aryuZZ[134];
   aryuZZ[134]=aryuZZ[10]*aryuZZ[134];
   aryuZZ[142]=21./2.*aryuZZ[10] + 3./4.*aryuZZ[13] - 17./4.*aryuZZ[12]
    + 13./12. + aryuZZ[224];
   aryuZZ[142]=aryuZZ[18]*aryuZZ[142];
   aryuZZ[106]=aryuZZ[120] + 1./2.*aryuZZ[142] + aryuZZ[134] + 1./4.*
   aryuZZ[106] + 1./8.*aryuZZ[133] + 1./2.*aryuZZ[115] + aryuZZ[123];
   aryuZZ[106]=aryuZZ[3]*aryuZZ[106];
   aryuZZ[115]= - 1./6.*aryuZZ[77] + aryuZZ[183] + 1./3.*aryuZZ[92];
   aryuZZ[120]=605./12. + aryuZZ[333];
   aryuZZ[120]= - 1./24.*aryuZZ[13] + 3./8.*aryuZZ[45] + 1./2.*
   aryuZZ[120] + aryuZZ[47];
   aryuZZ[120]=aryuZZ[17]*aryuZZ[120];
   aryuZZ[107]=aryuZZ[120] - 13./3.*aryuZZ[19] + 20./3.*aryuZZ[21] + 
   105./4.*aryuZZ[20] + aryuZZ[107] + aryuZZ[220] + aryuZZ[212] - 
   aryuZZ[76] + 1./2.*aryuZZ[115] + aryuZZ[209];
   aryuZZ[107]=aryuZZ[101]*aryuZZ[107];
   aryuZZ[115]= - 1./12.*aryuZZ[13] + aryuZZ[157] + 3*aryuZZ[188] + 
   aryuZZ[167];
   aryuZZ[115]=MMZ*aryuZZ[101]*aryuZZ[115];
   aryuZZ[120]= - 61./9. + aryuZZ[333];
   aryuZZ[120]=1./2.*aryuZZ[120] + aryuZZ[47];
   aryuZZ[114]=3./2. + aryuZZ[114];
   aryuZZ[114]=aryuZZ[10]*aryuZZ[114];
   aryuZZ[114]=aryuZZ[114] + aryuZZ[115] + aryuZZ[117] + aryuZZ[164] + 
   1./2.*aryuZZ[120] - aryuZZ[45];
   aryuZZ[114]=aryuZZ[10]*aryuZZ[114];
   aryuZZ[111]=aryuZZ[111] + aryuZZ[298] + aryuZZ[136] - 155./12. + 
   aryuZZ[285];
   aryuZZ[111]=aryuZZ[101]*aryuZZ[111];
   aryuZZ[111]=aryuZZ[111] + aryuZZ[181];
   aryuZZ[111]=aryuZZ[18]*aryuZZ[111];
   aryuZZ[115]= - 6283./432. + aryuZZ[88];
   aryuZZ[115]=7*aryuZZ[115] - 1./2.*aryuZZ[85];
   aryuZZ[115]=1./2.*aryuZZ[115] + 4*aryuZZ[87];
   aryuZZ[117]= - 17./16.*aryuZZ[13] - 1./3. - 19./16.*aryuZZ[12];
   aryuZZ[117]=aryuZZ[13]*aryuZZ[117];
   aryuZZ[106]=aryuZZ[106] + aryuZZ[111] + aryuZZ[114] + aryuZZ[109] + 
   aryuZZ[107] + aryuZZ[216] + 1./9.*aryuZZ[117] + 23./216.*aryuZZ[12]
    + 7./8.*aryuZZ[45] + 3./2.*aryuZZ[15] + aryuZZ[191] - 27./8.*
   aryuZZ[43] - 7./6.*aryuZZ[16] + 1./6.*aryuZZ[89] + aryuZZ[299] + 1./
   3.*aryuZZ[115] + aryuZZ[140];
   aryuZZ[106]=aryuZZ[73]*aryuZZ[106];
   aryuZZ[107]=6*aryuZZ[25] + aryuZZ[170];
   aryuZZ[107]=MMZ*aryuZZ[107];
   aryuZZ[109]=aryuZZ[1]*aryuZZ[175];
   aryuZZ[107]=aryuZZ[109] + aryuZZ[107] + aryuZZ[110] - 1./18.*
   aryuZZ[13] + aryuZZ[162] + aryuZZ[15] + 15./4.*aryuZZ[177] + 
   aryuZZ[174] - aryuZZ[28] + 1901./864. - 15*aryuZZ[31];
   aryuZZ[107]=aryuZZ[1]*aryuZZ[107];
   aryuZZ[111]= - 25./3. + aryuZZ[204];
   aryuZZ[111]=1./2.*aryuZZ[111] + aryuZZ[13];
   aryuZZ[111]=1./3.*aryuZZ[111] + aryuZZ[116];
   aryuZZ[111]=aryuZZ[141] + 1./2.*aryuZZ[111] + aryuZZ[250];
   aryuZZ[111]=aryuZZ[1]*aryuZZ[111];
   aryuZZ[114]= - 115./3. - 187*aryuZZ[12];
   aryuZZ[114]=1./2.*aryuZZ[114] + aryuZZ[323];
   aryuZZ[114]=1./3.*aryuZZ[114] + 11*aryuZZ[116];
   aryuZZ[114]= - 242./27.*aryuZZ[39] - 44./3.*aryuZZ[1] + 1./2.*
   aryuZZ[114] + 11*aryuZZ[250];
   aryuZZ[114]=aryuZZ[39]*aryuZZ[114];
   aryuZZ[115]=aryuZZ[1]*aryuZZ[130];
   aryuZZ[116]=aryuZZ[39]*aryuZZ[217];
   aryuZZ[115]=aryuZZ[115] + 11./3.*aryuZZ[116];
   aryuZZ[115]=aryuZZ[10]*aryuZZ[115];
   aryuZZ[116]=aryuZZ[138] + 11./9.*aryuZZ[145];
   aryuZZ[116]=aryuZZ[18]*aryuZZ[116];
   aryuZZ[117]= - 15./4. + aryuZZ[1];
   aryuZZ[117]=aryuZZ[1]*aryuZZ[117];
   aryuZZ[120]=121./9.*aryuZZ[39] - 685./36. + 22*aryuZZ[1];
   aryuZZ[120]=aryuZZ[39]*aryuZZ[120];
   aryuZZ[117]=aryuZZ[117] + 1./9.*aryuZZ[120];
   aryuZZ[117]=aryuZZ[6]*aryuZZ[117];
   aryuZZ[111]=aryuZZ[117] + aryuZZ[116] + aryuZZ[115] + aryuZZ[111] + 
   1./9.*aryuZZ[114];
   aryuZZ[111]=aryuZZ[6]*aryuZZ[111];
   aryuZZ[114]=10351./96. - 685*aryuZZ[31];
   aryuZZ[110]=11*aryuZZ[110] - 11./18.*aryuZZ[13] + 187./36.*
   aryuZZ[12] + 11*aryuZZ[15] + 685./36.*aryuZZ[177] - 253./24.*
   aryuZZ[14] + 1./9.*aryuZZ[114] - 11*aryuZZ[28];
   aryuZZ[113]=274./27.*aryuZZ[25] + 11*aryuZZ[113];
   aryuZZ[113]=MMZ*aryuZZ[113];
   aryuZZ[114]=aryuZZ[39]*aryuZZ[175];
   aryuZZ[109]=121./27.*aryuZZ[114] + 22./3.*aryuZZ[109] + 1./3.*
   aryuZZ[110] + aryuZZ[113];
   aryuZZ[109]=aryuZZ[39]*aryuZZ[109];
   aryuZZ[110]=aryuZZ[132] + aryuZZ[131] + aryuZZ[127] - 3*aryuZZ[8] + 
   aryuZZ[126];
   aryuZZ[110]=aryuZZ[1]*aryuZZ[110];
   aryuZZ[113]=aryuZZ[17] + 9*MMZ;
   aryuZZ[113]=aryuZZ[1]*aryuZZ[113];
   aryuZZ[114]=aryuZZ[39]*aryuZZ[165];
   aryuZZ[115]=9*aryuZZ[1] + 11*aryuZZ[39];
   aryuZZ[115]=aryuZZ[18]*aryuZZ[115];
   aryuZZ[113]=aryuZZ[115] + aryuZZ[113] + 11*aryuZZ[114];
   aryuZZ[113]=aryuZZ[6]*aryuZZ[113];
   aryuZZ[114]= - 4./3.*MMZ + aryuZZ[128] + aryuZZ[21] - aryuZZ[8] + 
   aryuZZ[124];
   aryuZZ[114]=aryuZZ[39]*aryuZZ[114];
   aryuZZ[110]=1./2.*aryuZZ[113] + aryuZZ[190] + aryuZZ[110] + 11./3.*
   aryuZZ[114];
   aryuZZ[110]=aryuZZ[3]*aryuZZ[110];
   aryuZZ[113]=aryuZZ[1]*aryuZZ[12];
   aryuZZ[114]=aryuZZ[39]*aryuZZ[12];
   aryuZZ[113]=aryuZZ[113] + 11./9.*aryuZZ[114];
   aryuZZ[114]= - aryuZZ[1]*aryuZZ[12];
   aryuZZ[115]= - aryuZZ[39]*aryuZZ[12];
   aryuZZ[114]=aryuZZ[114] + 11./9.*aryuZZ[115];
   aryuZZ[114]=aryuZZ[6]*aryuZZ[114];
   aryuZZ[113]=1./3.*aryuZZ[113] + aryuZZ[114];
   aryuZZ[114]= - MMZ*aryuZZ[12];
   aryuZZ[114]=aryuZZ[19] + aryuZZ[114];
   aryuZZ[115]= - 1./3. + aryuZZ[199];
   aryuZZ[115]=aryuZZ[18]*aryuZZ[115];
   aryuZZ[114]=1./3.*aryuZZ[114] + aryuZZ[115];
   aryuZZ[114]=aryuZZ[3]*aryuZZ[114];
   aryuZZ[115]= - aryuZZ[13]*aryuZZ[12];
   aryuZZ[115]=1./18.*aryuZZ[115] + 1./27.*aryuZZ[12] - 1./3.*
   aryuZZ[16] - 1./4. + 1./3.*aryuZZ[88];
   aryuZZ[116]= - aryuZZ[10]*aryuZZ[12];
   aryuZZ[114]=aryuZZ[114] + 1./2.*aryuZZ[115] + 1./3.*aryuZZ[116];
   aryuZZ[114]=aryuZZ[73]*aryuZZ[114];
   aryuZZ[113]=1./3.*aryuZZ[113] + aryuZZ[114];
   aryuZZ[113]=aryuZZ[272]*aryuZZ[113];
   aryuZZ[114]=aryuZZ[137] + aryuZZ[166];
   aryuZZ[114]=1./2.*aryuZZ[114] + aryuZZ[122];
   aryuZZ[114]=MMH*aryuZZ[73]*aryuZZ[114];
   aryuZZ[106]=1./4.*aryuZZ[113] + aryuZZ[114] + aryuZZ[106] + 
   aryuZZ[110] + aryuZZ[111] + aryuZZ[370] + aryuZZ[107] + 1./3.*
   aryuZZ[109];
   aryuZZ[106]=aryuZZ[272]*aryuZZ[106];
   aryuZZ[107]= - 713./36. + aryuZZ[23];
   aryuZZ[109]= - 289./27. - aryuZZ[7];
   aryuZZ[109]=aryuZZ[13]*aryuZZ[109];
   aryuZZ[110]= - 548./9.*aryuZZ[25] + 73./9.*aryuZZ[27] + aryuZZ[26];
   aryuZZ[110]=MMZ*aryuZZ[110];
   aryuZZ[111]= - 11*aryuZZ[177];
   aryuZZ[113]=19./9.*aryuZZ[7] - 85./27. + aryuZZ[111];
   aryuZZ[113]=aryuZZ[1]*aryuZZ[113];
   aryuZZ[111]= - 31./9. + aryuZZ[111];
   aryuZZ[111]=1./3.*aryuZZ[111] + aryuZZ[7];
   aryuZZ[111]=aryuZZ[39]*aryuZZ[111];
   aryuZZ[107]=11./9.*aryuZZ[111] + 2./3.*aryuZZ[113] + 1./3.*
   aryuZZ[110] + aryuZZ[292] + 1./4.*aryuZZ[109] + 1./6.*aryuZZ[7] - 95.
   /18.*aryuZZ[12] - 329./27.*aryuZZ[177] + 1370./27.*aryuZZ[31] + 1./6.
   *aryuZZ[107] - aryuZZ[32];
   aryuZZ[107]=aryuZZ[39]*aryuZZ[107];
   aryuZZ[109]=689./36. + 11*aryuZZ[23];
   aryuZZ[110]= - 89./3. - aryuZZ[7];
   aryuZZ[110]=aryuZZ[13]*aryuZZ[110];
   aryuZZ[111]=aryuZZ[195] + aryuZZ[26];
   aryuZZ[111]=1./3.*aryuZZ[111] - 12*aryuZZ[25];
   aryuZZ[111]=MMZ*aryuZZ[111];
   aryuZZ[113]=1./9.*aryuZZ[7] - 7./27. - aryuZZ[177];
   aryuZZ[113]=aryuZZ[1]*aryuZZ[113];
   aryuZZ[109]=aryuZZ[113] + aryuZZ[111] + aryuZZ[216] + 1./36.*
   aryuZZ[110] + 1./54.*aryuZZ[7] - 77./54.*aryuZZ[12] - 7*aryuZZ[177]
    + 30*aryuZZ[31] + 1./18.*aryuZZ[109] - aryuZZ[32];
   aryuZZ[109]=aryuZZ[1]*aryuZZ[109];
   aryuZZ[110]=aryuZZ[310] - MMZ;
   aryuZZ[110]=aryuZZ[1]*aryuZZ[110];
   aryuZZ[111]=aryuZZ[19] - 1./3.*MMZ;
   aryuZZ[111]=aryuZZ[39]*aryuZZ[111];
   aryuZZ[110]=aryuZZ[125] + aryuZZ[110] + 11./3.*aryuZZ[111];
   aryuZZ[110]=aryuZZ[6]*aryuZZ[110];
   aryuZZ[111]=20./27. - aryuZZ[7];
   aryuZZ[111]=MMZ*aryuZZ[111];
   aryuZZ[111]=aryuZZ[146] + aryuZZ[111];
   aryuZZ[111]=aryuZZ[39]*aryuZZ[111];
   aryuZZ[113]=4./3. - aryuZZ[7];
   aryuZZ[113]=MMZ*aryuZZ[113];
   aryuZZ[113]= - aryuZZ[19] + 1./3.*aryuZZ[113];
   aryuZZ[113]=aryuZZ[1]*aryuZZ[113];
   aryuZZ[110]=aryuZZ[110] + aryuZZ[139] + aryuZZ[113] + aryuZZ[111];
   aryuZZ[110]=aryuZZ[3]*aryuZZ[110];
   aryuZZ[111]=151./3. + aryuZZ[119];
   aryuZZ[111]=aryuZZ[1]*aryuZZ[111];
   aryuZZ[113]=53./27. - aryuZZ[7];
   aryuZZ[113]=aryuZZ[39]*aryuZZ[113];
   aryuZZ[111]=11*aryuZZ[113] + 2./3.*aryuZZ[111] - 11*aryuZZ[11] + 88./
   3.*aryuZZ[13] + 784./9. + 95./2.*aryuZZ[12];
   aryuZZ[111]=aryuZZ[39]*aryuZZ[111];
   aryuZZ[113]=7 - 2*aryuZZ[1];
   aryuZZ[113]=aryuZZ[1]*aryuZZ[113];
   aryuZZ[114]= - 242./9.*aryuZZ[39] + 329./9. - 44*aryuZZ[1];
   aryuZZ[114]=aryuZZ[39]*aryuZZ[114];
   aryuZZ[113]=aryuZZ[113] + 1./9.*aryuZZ[114];
   aryuZZ[113]=aryuZZ[6]*aryuZZ[113];
   aryuZZ[114]=13./3. - aryuZZ[7];
   aryuZZ[114]=aryuZZ[1]*aryuZZ[114];
   aryuZZ[114]=1./3.*aryuZZ[114] - aryuZZ[11] + 8./3.*aryuZZ[13] + 10
    + 77./18.*aryuZZ[12];
   aryuZZ[114]=aryuZZ[1]*aryuZZ[114];
   aryuZZ[111]=aryuZZ[113] + aryuZZ[114] + 1./9.*aryuZZ[111];
   aryuZZ[111]=aryuZZ[6]*aryuZZ[111];
   aryuZZ[113]=aryuZZ[150] + aryuZZ[168] + aryuZZ[137] + 1./2.*
   aryuZZ[80];
   aryuZZ[113]=MMH*aryuZZ[73]*aryuZZ[113];
   aryuZZ[105]=aryuZZ[106] + aryuZZ[113] + aryuZZ[105] + aryuZZ[110] + 
   aryuZZ[111] + aryuZZ[325] + aryuZZ[109] + 1./3.*aryuZZ[107];
   aryuZZ[105]=aryuZZ[272]*aryuZZ[105];
   aryuZZ[106]=aryuZZ[16] - 2*aryuZZ[89] + aryuZZ[90] - 2 - aryuZZ[97];
   aryuZZ[106]=aryuZZ[108]*aryuZZ[106];
   aryuZZ[107]=82*aryuZZ[90] - 29 - 82*aryuZZ[97];
   aryuZZ[106]=8*aryuZZ[106] + 82./3.*aryuZZ[16] - 164./3.*aryuZZ[89]
    + 1./3.*aryuZZ[107] + aryuZZ[121];
   aryuZZ[106]=aryuZZ[100]*aryuZZ[106];
   aryuZZ[107]=aryuZZ[82] - 2*aryuZZ[95] - aryuZZ[83];
   aryuZZ[107]=aryuZZ[108]*aryuZZ[107];
   aryuZZ[109]=46*aryuZZ[82] - 58*aryuZZ[95] - 23*aryuZZ[83];
   aryuZZ[107]=4*aryuZZ[107] + 1./3.*aryuZZ[109] - aryuZZ[81];
   aryuZZ[107]=aryuZZ[108]*aryuZZ[107];
   aryuZZ[109]= - 100./3.*aryuZZ[95] - 7*aryuZZ[83];
   aryuZZ[110]=41./3. + aryuZZ[112];
   aryuZZ[110]=aryuZZ[13]*aryuZZ[100]*aryuZZ[110];
   aryuZZ[111]= - 1 - aryuZZ[90];
   aryuZZ[111]=aryuZZ[101]*aryuZZ[111];
   aryuZZ[106]=8*aryuZZ[111] + 2*aryuZZ[110] + 2*aryuZZ[205] + 
   aryuZZ[106] + 4*aryuZZ[107] - 41./3.*aryuZZ[81] - 16*aryuZZ[80] + 2*
   aryuZZ[109] + 65*aryuZZ[82];
   aryuZZ[106]=MMZ*aryuZZ[106];
   aryuZZ[107]=149./9. + aryuZZ[112];
   aryuZZ[107]=aryuZZ[108]*aryuZZ[107];
   aryuZZ[109]=12*aryuZZ[108];
   aryuZZ[110]=119./3. + aryuZZ[109];
   aryuZZ[110]=aryuZZ[12]*aryuZZ[110];
   aryuZZ[107]=aryuZZ[110] + 1495./9. + 4*aryuZZ[107];
   aryuZZ[110]=44./9. + aryuZZ[108];
   aryuZZ[110]=aryuZZ[108]*aryuZZ[110];
   aryuZZ[110]=1813./9. + 16*aryuZZ[110];
   aryuZZ[110]=aryuZZ[13]*aryuZZ[110];
   aryuZZ[107]=2*aryuZZ[107] + aryuZZ[110];
   aryuZZ[107]=aryuZZ[13]*aryuZZ[107];
   aryuZZ[110]= - 47./3. + aryuZZ[129];
   aryuZZ[110]=aryuZZ[13]*aryuZZ[110];
   aryuZZ[110]=aryuZZ[110] + aryuZZ[218] + 1./3. + aryuZZ[129];
   aryuZZ[110]=MMZ*aryuZZ[110];
   aryuZZ[111]= - 5 - 3*aryuZZ[13];
   aryuZZ[111]=aryuZZ[19]*aryuZZ[111];
   aryuZZ[112]=1 - aryuZZ[13];
   aryuZZ[112]=aryuZZ[17]*aryuZZ[112];
   aryuZZ[110]=aryuZZ[110] + aryuZZ[112] - aryuZZ[20] + 2*aryuZZ[111];
   aryuZZ[111]=2 + aryuZZ[213];
   aryuZZ[111]=aryuZZ[3]*aryuZZ[118]*aryuZZ[111];
   aryuZZ[110]=2*aryuZZ[110] + aryuZZ[111];
   aryuZZ[110]=aryuZZ[3]*aryuZZ[110];
   aryuZZ[111]=109./9. - aryuZZ[88];
   aryuZZ[111]=6*aryuZZ[108] + aryuZZ[16] + 2*aryuZZ[111] + aryuZZ[90];
   aryuZZ[111]=aryuZZ[108]*aryuZZ[111];
   aryuZZ[112]=1117 - 244*aryuZZ[88];
   aryuZZ[113]= - aryuZZ[76] + 4*aryuZZ[22];
   aryuZZ[113]=aryuZZ[100]*aryuZZ[113];
   aryuZZ[109]=71./3. + aryuZZ[109];
   aryuZZ[109]=aryuZZ[12]*aryuZZ[109];
   aryuZZ[114]=aryuZZ[11]*aryuZZ[262];
   aryuZZ[115]=aryuZZ[19]*aryuZZ[219];
   aryuZZ[116]=2*aryuZZ[100] + aryuZZ[147];
   aryuZZ[116]=aryuZZ[17]*aryuZZ[116];
   aryuZZ[117]=aryuZZ[10]*MMZ*aryuZZ[100];
   aryuZZ[106]=2*aryuZZ[110] + 2*aryuZZ[117] + aryuZZ[106] + 4*
   aryuZZ[116] + 8*aryuZZ[115] + 6*aryuZZ[114] + aryuZZ[107] + 2*
   aryuZZ[109] + 4*aryuZZ[113] + 8*aryuZZ[111] + 58*aryuZZ[16] - 8*
   aryuZZ[89] + 190./9.*aryuZZ[90] + 1./3.*aryuZZ[112] + 2*aryuZZ[97];
   aryuZZ[106]=aryuZZ[73]*aryuZZ[106];
   aryuZZ[107]= - 11*aryuZZ[32] - 19./9. - 11*aryuZZ[23];
   aryuZZ[107]=aryuZZ[108]*aryuZZ[107];
   aryuZZ[109]=8./27.*aryuZZ[256] + aryuZZ[135];
   aryuZZ[109]=aryuZZ[13]*aryuZZ[109];
   aryuZZ[110]= - aryuZZ[24] - 2./9.*aryuZZ[26];
   aryuZZ[110]=aryuZZ[108]*aryuZZ[110];
   aryuZZ[111]= - 2*aryuZZ[24];
   aryuZZ[112]=aryuZZ[111] - 17./81.*aryuZZ[27];
   aryuZZ[110]=2*aryuZZ[110] + 272./81.*aryuZZ[25] + 4*aryuZZ[112] - 17.
   /9.*aryuZZ[26];
   aryuZZ[110]=MMZ*aryuZZ[110];
   aryuZZ[112]=5*aryuZZ[29];
   aryuZZ[113]= - 5*aryuZZ[16];
   aryuZZ[107]=200./729.*aryuZZ[39] + 80./243.*aryuZZ[1] + aryuZZ[110]
    + aryuZZ[109] + aryuZZ[135] + 1./3.*aryuZZ[107] + 134./81.*
   aryuZZ[177] + aryuZZ[113] - 680./81.*aryuZZ[31] - 101./9.*aryuZZ[32]
    - 55./9.*aryuZZ[23] - 511./81. + aryuZZ[112];
   aryuZZ[107]=aryuZZ[39]*aryuZZ[107];
   aryuZZ[109]= - 5*aryuZZ[32] + 1./3. - 5*aryuZZ[23];
   aryuZZ[109]=aryuZZ[108]*aryuZZ[109];
   aryuZZ[109]=aryuZZ[109] + 10*aryuZZ[177] + aryuZZ[113] - 40*
   aryuZZ[31] - 13*aryuZZ[32] - 17./3.*aryuZZ[23] - 127./9. + 
   aryuZZ[112];
   aryuZZ[110]=aryuZZ[111] - aryuZZ[27];
   aryuZZ[111]= - aryuZZ[108]*aryuZZ[24];
   aryuZZ[110]=2*aryuZZ[111] + 16*aryuZZ[25] + 4*aryuZZ[110] - 
   aryuZZ[26];
   aryuZZ[110]=MMZ*aryuZZ[110];
   aryuZZ[111]=7./3. + 2*aryuZZ[108];
   aryuZZ[111]=8./9.*aryuZZ[111] - aryuZZ[7];
   aryuZZ[111]=aryuZZ[13]*aryuZZ[111];
   aryuZZ[109]=8./81.*aryuZZ[1] + 1./3.*aryuZZ[110] + aryuZZ[111] + 1./
   3.*aryuZZ[109] - aryuZZ[7];
   aryuZZ[109]=aryuZZ[1]*aryuZZ[109];
   aryuZZ[110]= - 2./3. - aryuZZ[108];
   aryuZZ[111]= - 49./3. + aryuZZ[153];
   aryuZZ[111]=aryuZZ[13]*aryuZZ[111];
   aryuZZ[110]= - 8./9.*aryuZZ[1] + 8*aryuZZ[110] + aryuZZ[111];
   aryuZZ[110]=aryuZZ[1]*aryuZZ[110];
   aryuZZ[108]= - 34./9. - 5*aryuZZ[108];
   aryuZZ[111]= - 61./3. + aryuZZ[153];
   aryuZZ[111]=aryuZZ[13]*aryuZZ[111];
   aryuZZ[108]= - 200./27.*aryuZZ[39] - 80./9.*aryuZZ[1] + 8*
   aryuZZ[108] + 5*aryuZZ[111];
   aryuZZ[108]=aryuZZ[39]*aryuZZ[108];
   aryuZZ[111]= - 5 + aryuZZ[158];
   aryuZZ[111]=aryuZZ[1]*aryuZZ[111];
   aryuZZ[112]=100./3.*aryuZZ[39] - 67./3. + 40*aryuZZ[1];
   aryuZZ[112]=aryuZZ[39]*aryuZZ[112];
   aryuZZ[111]=aryuZZ[111] + 1./9.*aryuZZ[112];
   aryuZZ[111]=aryuZZ[6]*aryuZZ[111];
   aryuZZ[108]=aryuZZ[111] + aryuZZ[110] + 1./3.*aryuZZ[108];
   aryuZZ[108]=aryuZZ[6]*aryuZZ[108];
   aryuZZ[110]= - aryuZZ[1]*MMZ;
   aryuZZ[111]= - aryuZZ[39]*MMZ;
   aryuZZ[110]=aryuZZ[110] + 5./3.*aryuZZ[111];
   aryuZZ[111]=aryuZZ[1]*MMZ;
   aryuZZ[112]=aryuZZ[39]*MMZ;
   aryuZZ[111]=aryuZZ[111] + 5./3.*aryuZZ[112];
   aryuZZ[111]=aryuZZ[6]*aryuZZ[111];
   aryuZZ[110]=1./3.*aryuZZ[110] + aryuZZ[111];
   aryuZZ[110]=aryuZZ[3]*aryuZZ[110];
   aryuZZ[107]=8./3.*aryuZZ[110] + 2./3.*aryuZZ[108] + aryuZZ[109] + 
   aryuZZ[107];
   aryuZZ[108]=MMH*aryuZZ[73]*aryuZZ[81];

      yuZZret = aryuZZ[102] + aryuZZ[103] + aryuZZ[104] + aryuZZ[105]
       + aryuZZ[106] + 2*aryuZZ[107] + 10./3.*aryuZZ[108];
      return yuZZret;
}
