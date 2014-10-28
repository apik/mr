#include <bb.hpp>
std::complex<long double>
bb::m20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armbbbar[299], mbbbarret;

    armbbbar[1]=double(nL + nH);
    armbbbar[2]=pow(CW,-1);
    armbbbar[3]=pow(MMH,-1);
    armbbbar[4]=pow(MMZ,-1);
    armbbbar[5]=pow(SW,-1);
    armbbbar[6]=Tsil::I2(0,0,MMZ,mu2);
    armbbbar[7]=Tsil::I2(0,0,MMW,mu2);
    armbbbar[8]=Tsil::A(MMH,mu2);
    armbbbar[9]=std::real(Tsil::B(0,0,MMZ,mu2));
    armbbbar[10]=std::real(Tsil::B(0,0,MMW,mu2));
    armbbbar[11]=Tsil::A(MMZ,mu2);
    armbbbar[12]=Tsil::A(MMW,mu2);
    armbbbar[13]=Tsil::A(MMt,mu2);
    armbbbar[14]=Tsil::A(MMb,mu2);
    armbbbar[15]=pow(MMb,-1);
    armbbbar[16]=Tsil::Aeps(MMZ,mu2);
    armbbbar[17]=Tsil::Aeps(MMW,mu2);
    armbbbar[18]=Tsil::Aeps(MMb,mu2);
    armbbbar[19]=prot0bb0b->Tvxu(0);
    armbbbar[20]=double(nL);
    armbbbar[21]=double(boson);
    armbbbar[22]=det(MMW,MMH,MMt);
    armbbbar[23]=Tsil::I2(MMW,MMH,MMt,mu2);
    armbbbar[24]=Tsil::Aeps(MMH,mu2);
    armbbbar[25]=Tsil::Aeps(MMt,mu2);
    armbbbar[26]=det(MMW,MMZ,MMt);
    armbbbar[27]=Tsil::I2(MMW,MMZ,MMt,mu2);
    armbbbar[28]=Tsil::I2(MMH,MMH,MMH,mu2);
    armbbbar[29]=Tsil::I2(MMH,MMZ,MMZ,mu2);
    armbbbar[30]=Tsil::I2(MMH,MMW,MMW,mu2);
    armbbbar[31]=Tsil::I2(MMH,MMt,MMt,mu2);
    armbbbar[32]=Tsil::I2(MMZ,MMt,MMt,mu2);
    armbbbar[33]=Tsil::I2(MMW,MMW,MMZ,mu2);
    armbbbar[34]=Tsil::I2(0,MMZ,MMH,mu2);
    armbbbar[35]=Tsil::I2(0,MMW,MMt,mu2);
    armbbbar[36]=Tsil::I2(0,0,MMH,mu2);
    armbbbar[37]=Tsil::B(MMH,MMH,MMH,mu2);
    armbbbar[38]=Tsil::B(MMH,MMt,MMt,mu2);
    armbbbar[39]=Tsil::B(MMZ,MMH,MMZ,mu2);
    armbbbar[40]=Tsil::B(MMZ,MMZ,MMH,mu2);
    armbbbar[41]=Tsil::B(MMZ,MMt,MMt,mu2);
    armbbbar[42]=Tsil::B(MMW,MMH,MMW,mu2);
    armbbbar[43]=Tsil::B(MMW,MMZ,MMW,mu2);
    armbbbar[44]=Tsil::B(MMW,MMW,MMH,mu2);
    armbbbar[45]=Tsil::B(MMW,MMW,MMZ,mu2);
    armbbbar[46]=Tsil::B(MMt,MMt,MMH,mu2);
    armbbbar[47]=Tsil::B(MMt,MMt,MMZ,mu2);
    armbbbar[48]=Tsil::B(0,MMt,MMW,mu2);
    armbbbar[49]=std::real(Tsil::B(MMt,MMt,MMH,mu2));
    armbbbar[50]=std::real(Tsil::B(0,MMW,MMt,mu2));
    armbbbar[51]=pow(MMt,-1);
    armbbbar[52]=1/(MMt - MMW);
    armbbbar[53]=Tsil::I2(0,0,MMt,mu2);
    armbbbar[54]=1/(4*MMt - MMZ);
    armbbbar[55]=1/( - 4*MMW + MMZ);
    armbbbar[56]=1/( - 4*MMW + MMH);
    armbbbar[57]=1/( - 4*MMZ + MMH);
   armbbbar[58]= - 9*armbbbar[37];
   armbbbar[59]= - armbbbar[57]*armbbbar[29];
   armbbbar[60]=1./2.*armbbbar[59] - 331./32. + armbbbar[58];
   armbbbar[60]=3*armbbbar[60] - 13./2.*armbbbar[50];
   armbbbar[61]=armbbbar[25] + armbbbar[17] + armbbbar[24] - 
   armbbbar[23];
   armbbbar[62]=armbbbar[22]*armbbbar[61];
   armbbbar[63]=armbbbar[8]*armbbbar[22];
   armbbbar[64]=armbbbar[62] + armbbbar[63];
   armbbbar[65]= - MMH*armbbbar[22];
   armbbbar[66]=5./4.*armbbbar[64] + armbbbar[65];
   armbbbar[66]=MMH*armbbbar[66];
   armbbbar[67]= - armbbbar[54]*armbbbar[32];
   armbbbar[68]=armbbbar[24]*armbbbar[57];
   armbbbar[69]=7./8.*armbbbar[54];
   armbbbar[70]=3*armbbbar[57] + armbbbar[69];
   armbbbar[71]=armbbbar[16]*armbbbar[70];
   armbbbar[72]=armbbbar[25]*armbbbar[54];
   armbbbar[73]= - armbbbar[30] + armbbbar[24];
   armbbbar[73]=1./2.*armbbbar[73] + armbbbar[17];
   armbbbar[73]=armbbbar[56]*armbbbar[73];
   armbbbar[74]=1./2.*armbbbar[57] + armbbbar[56];
   armbbbar[74]=armbbbar[8]*armbbbar[74];
   armbbbar[70]=armbbbar[11]*armbbbar[70];
   armbbbar[75]=MMH*armbbbar[22];
   armbbbar[76]=armbbbar[56] + 5./4.*armbbbar[75];
   armbbbar[76]=armbbbar[12]*armbbbar[76];
   armbbbar[77]= - 9./4.*armbbbar[48];
   armbbbar[78]= - 3*armbbbar[38];
   armbbbar[79]=3*armbbbar[44];
   armbbbar[80]=3./2.*armbbbar[40];
   armbbbar[81]= - 3*armbbbar[39];
   armbbbar[82]=33./2.*armbbbar[45];
   armbbbar[83]= - 1./2.*armbbbar[9];
   armbbbar[60]=3./2.*armbbbar[76] + 3./2.*armbbbar[66] + 1./4.*
   armbbbar[70] + 3./4.*armbbbar[74] + 3./2.*armbbbar[73] + 
   armbbbar[83] - 5./4.*armbbbar[42] + armbbbar[77] + 7./16.*
   armbbbar[72] + armbbbar[81] + 9./4.*armbbbar[47] - 165./16.*
   armbbbar[43] + 1./4.*armbbbar[71] + 3./8.*armbbbar[68] + 
   armbbbar[80] + armbbbar[82] + armbbbar[78] + 7./32.*armbbbar[67] - 7.
   /16.*armbbbar[41] + 1./4.*armbbbar[60] + armbbbar[79];
   armbbbar[66]= - armbbbar[25] - armbbbar[17] - armbbbar[24] + 
   armbbbar[23];
   armbbbar[70]=armbbbar[56]*armbbbar[66];
   armbbbar[71]=21./4. + armbbbar[70];
   armbbbar[71]=armbbbar[56]*armbbbar[71];
   armbbbar[73]=pow(armbbbar[56],2);
   armbbbar[74]= - armbbbar[8]*armbbbar[73];
   armbbbar[76]= - armbbbar[12]*armbbbar[73];
   armbbbar[84]=3./4.*armbbbar[3];
   armbbbar[85]= - armbbbar[13]*armbbbar[73];
   armbbbar[85]=armbbbar[85] + armbbbar[84] + armbbbar[76] + 
   armbbbar[71] + armbbbar[74];
   armbbbar[86]=MMZ*armbbbar[73];
   armbbbar[85]=1./2.*armbbbar[85] + 3*armbbbar[86];
   armbbbar[85]=MMZ*armbbbar[85];
   armbbbar[86]=armbbbar[8]*armbbbar[73];
   armbbbar[87]= - 1./4.*armbbbar[56] + armbbbar[86];
   armbbbar[88]= - 3./2.*armbbbar[3];
   armbbbar[76]=armbbbar[88] + 1./2.*armbbbar[87] + armbbbar[76];
   armbbbar[76]=armbbbar[13]*armbbbar[76];
   armbbbar[87]= - 3*armbbbar[41];
   armbbbar[89]=armbbbar[87] + 149./2. - 9*armbbbar[50];
   armbbbar[90]=33./8.*armbbbar[43];
   armbbbar[91]= - 1./4.*armbbbar[25] - 9./4.*armbbbar[17] + 1./4.*
   armbbbar[23] + armbbbar[30] - 5./4.*armbbbar[24];
   armbbbar[91]=armbbbar[56]*armbbbar[91];
   armbbbar[92]= - 9./4.*armbbbar[56] + armbbbar[86];
   armbbbar[92]=armbbbar[12]*armbbbar[92];
   armbbbar[93]=1./4.*armbbbar[11] + armbbbar[12];
   armbbbar[93]=armbbbar[3]*armbbbar[93];
   armbbbar[94]= - armbbbar[8]*armbbbar[56];
   armbbbar[76]=3*armbbbar[85] + 3./2.*armbbbar[76] + 9./4.*
   armbbbar[93] + 3./4.*armbbbar[92] + 15./16.*armbbbar[94] + 3./4.*
   armbbbar[91] + 7./8.*armbbbar[9] + armbbbar[77] + 3./4.*armbbbar[39]
    + 21./8.*armbbbar[47] + armbbbar[90] - 363./16.*armbbbar[45] + 1./8.
   *armbbbar[89] + armbbbar[78];
   armbbbar[76]=MMZ*armbbbar[76];
   armbbbar[77]= - 3*armbbbar[50];
   armbbbar[85]= - 113./4. + armbbbar[77];
   armbbbar[85]=5*armbbbar[85] + armbbbar[87];
   armbbbar[87]=armbbbar[8]*armbbbar[56];
   armbbbar[89]= - armbbbar[12]*armbbbar[56];
   armbbbar[91]=1./2.*armbbbar[9];
   armbbbar[92]= - 3*armbbbar[48];
   armbbbar[85]=3./4.*armbbbar[89] + 3./4.*armbbbar[87] + armbbbar[91]
    + armbbbar[92] + 17./4.*armbbbar[47] + armbbbar[90] - 33./2.*
   armbbbar[45] + 1./16.*armbbbar[85] + armbbbar[78];
   armbbbar[85]=armbbbar[12]*armbbbar[85];
   armbbbar[90]=3*armbbbar[50];
   armbbbar[93]=armbbbar[41] - 1./4. + armbbbar[90];
   armbbbar[95]=3*armbbbar[38];
   armbbbar[96]= - 55./4.*armbbbar[43];
   armbbbar[97]=3*armbbbar[48];
   armbbbar[93]=armbbbar[42] + armbbbar[97] - armbbbar[39] - 7./2.*
   armbbbar[47] + armbbbar[96] + 121./4.*armbbbar[45] + 1./2.*
   armbbbar[93] + armbbbar[95];
   armbbbar[98]= - 3*armbbbar[12];
   armbbbar[99]= - 1./2.*armbbbar[11];
   armbbbar[100]=armbbbar[99] + armbbbar[98];
   armbbbar[100]=armbbbar[3]*armbbbar[100];
   armbbbar[101]=3*armbbbar[89];
   armbbbar[102]=armbbbar[13]*armbbbar[3];
   armbbbar[93]=9*armbbbar[102] + 9./2.*armbbbar[100] + armbbbar[101]
    + 3./2.*armbbbar[87] + 3*armbbbar[93] - 7./2.*armbbbar[9];
   armbbbar[93]=armbbbar[13]*armbbbar[93];
   armbbbar[100]= - armbbbar[41] - 49./2. + armbbbar[77];
   armbbbar[103]= - 121./4.*armbbbar[45];
   armbbbar[100]=armbbbar[42] + armbbbar[92] + armbbbar[39] + 7./2.*
   armbbbar[47] - 11./4.*armbbbar[43] + armbbbar[103] + 1./2.*
   armbbbar[100] - 5*armbbbar[38];
   armbbbar[104]=7./2.*armbbbar[9];
   armbbbar[100]=armbbbar[101] + 9./2.*armbbbar[87] + 3*armbbbar[100]
    + armbbbar[104];
   armbbbar[100]=armbbbar[12]*armbbbar[100];
   armbbbar[101]=11*armbbbar[33];
   armbbbar[105]=1./2.*armbbbar[30];
   armbbbar[106]= - 3./8.*armbbbar[31] + armbbbar[101] + armbbbar[105];
   armbbbar[107]=5 + armbbbar[95];
   armbbbar[107]=1./4.*armbbbar[107] - armbbbar[42];
   armbbbar[107]=MMH*armbbbar[107];
   armbbbar[108]=1./2.*armbbbar[11];
   armbbbar[109]=armbbbar[108] + armbbbar[12];
   armbbbar[110]=armbbbar[3]*armbbbar[12]*armbbbar[109];
   armbbbar[111]= - 3./4.*armbbbar[8];
   armbbbar[93]=1./2.*armbbbar[93] + 9./4.*armbbbar[110] + 1./2.*
   armbbbar[100] + armbbbar[107] - 33./4.*armbbbar[11] + armbbbar[111]
    + 535./16.*armbbbar[25] - 823./16.*armbbbar[17] - 33./4.*
   armbbbar[16] + 3./8.*armbbbar[23] - 3./4.*armbbbar[24] - 287./16.*
   armbbbar[27] + 3*armbbbar[106] - 109./16.*armbbbar[32];
   armbbbar[100]= - 3./2.*armbbbar[17];
   armbbbar[106]= - 1./2.*armbbbar[25];
   armbbbar[107]=armbbbar[106] + armbbbar[100] + 1./2.*armbbbar[23] + 
   armbbbar[105] - armbbbar[24];
   armbbbar[107]=armbbbar[56]*armbbbar[107];
   armbbbar[112]=armbbbar[3]*armbbbar[12];
   armbbbar[88]= - armbbbar[56] + armbbbar[88];
   armbbbar[88]=armbbbar[13]*armbbbar[88];
   armbbbar[113]= - 33./4.*armbbbar[43];
   armbbbar[88]=1./2.*armbbbar[88] + 3./4.*armbbbar[112] + 3./2.*
   armbbbar[89] + armbbbar[94] + armbbbar[107] + armbbbar[42] + 
   armbbbar[113] + 309./32. - armbbbar[38];
   armbbbar[89]=MMZ*armbbbar[56];
   armbbbar[88]=1./2.*armbbbar[88] + 3*armbbbar[89];
   armbbbar[88]=MMZ*armbbbar[88];
   armbbbar[88]=1./2.*armbbbar[93] + 3*armbbbar[88];
   armbbbar[88]=MMZ*armbbbar[88];
   armbbbar[89]=pow(MMH,2);
   armbbbar[93]= - armbbbar[89]*armbbbar[42];
   armbbbar[114]= - armbbbar[13]*MMH;
   armbbbar[115]=1./2.*armbbbar[114];
   armbbbar[116]=armbbbar[12]*MMH;
   armbbbar[117]=armbbbar[115] + armbbbar[93] + armbbbar[116];
   armbbbar[117]=armbbbar[13]*armbbbar[117];
   armbbbar[118]= - armbbbar[30] + armbbbar[31];
   armbbbar[118]= - armbbbar[25] + 1./2.*armbbbar[118] + armbbbar[17];
   armbbbar[118]=armbbbar[89]*armbbbar[118];
   armbbbar[119]=armbbbar[89]*armbbbar[42];
   armbbbar[120]= - armbbbar[12]*MMH;
   armbbbar[119]=armbbbar[119] + 1./2.*armbbbar[120];
   armbbbar[119]=armbbbar[12]*armbbbar[119];
   armbbbar[117]=armbbbar[117] + armbbbar[118] + armbbbar[119];
   armbbbar[118]= - 1./2.*armbbbar[30];
   armbbbar[119]= - 1./2.*armbbbar[17];
   armbbbar[121]=1./2.*armbbbar[25];
   armbbbar[122]=armbbbar[121] + armbbbar[119] + 3./2.*armbbbar[23] + 
   armbbbar[118] - armbbbar[31];
   armbbbar[122]=MMH*armbbbar[122];
   armbbbar[123]=5*armbbbar[33] - armbbbar[32];
   armbbbar[100]=3./2.*armbbbar[25] + armbbbar[100] + 1./4.*
   armbbbar[123] - armbbbar[27];
   armbbbar[113]=armbbbar[42] + armbbbar[113] - 3./16. - armbbbar[38];
   armbbbar[113]=armbbbar[12]*armbbbar[113];
   armbbbar[123]= - armbbbar[42] + 33./4.*armbbbar[43] + 3./16. + 
   armbbbar[38];
   armbbbar[123]=armbbbar[13]*armbbbar[123];
   armbbbar[100]=armbbbar[123] + 11./2.*armbbbar[100] + armbbbar[113];
   armbbbar[100]=MMZ*armbbbar[100];
   armbbbar[113]=3./4.*armbbbar[38];
   armbbbar[123]=armbbbar[113] - armbbbar[42];
   armbbbar[123]=MMH*armbbbar[123];
   armbbbar[124]=armbbbar[123] + 5./16.*armbbbar[12];
   armbbbar[124]=armbbbar[12]*armbbbar[124];
   armbbbar[125]= - 3./4.*armbbbar[38];
   armbbbar[126]=armbbbar[125] + armbbbar[42];
   armbbbar[126]=MMH*armbbbar[126];
   armbbbar[127]=7./8.*armbbbar[13] + armbbbar[126] - 19./16.*
   armbbbar[12];
   armbbbar[127]=armbbbar[13]*armbbbar[127];
   armbbbar[100]=3*armbbbar[100] + armbbbar[127] + 1./2.*armbbbar[122]
    + armbbbar[124];
   armbbbar[100]=MMZ*armbbbar[100];
   armbbbar[122]=armbbbar[12] - 1./2.*armbbbar[13];
   armbbbar[122]=armbbbar[13]*armbbbar[122];
   armbbbar[124]=pow(armbbbar[12],2);
   armbbbar[122]= - 1./2.*armbbbar[124] + armbbbar[122];
   armbbbar[127]=pow(MMZ,2);
   armbbbar[122]=armbbbar[52]*armbbbar[127]*armbbbar[122];
   armbbbar[100]=9./16.*armbbbar[122] + 1./4.*armbbbar[117] + 
   armbbbar[100];
   armbbbar[100]=armbbbar[52]*armbbbar[100];
   armbbbar[117]= - 7*armbbbar[30];
   armbbbar[128]=armbbbar[117] - 15./2.*armbbbar[31];
   armbbbar[129]=1 + armbbbar[42];
   armbbbar[129]=MMH*armbbbar[129];
   armbbbar[128]=armbbbar[129] - armbbbar[8] - 3./4.*armbbbar[25] - 5./
   4.*armbbbar[17] + 33./4.*armbbbar[23] + 1./2.*armbbbar[128] - 
   armbbbar[24];
   armbbbar[128]=MMH*armbbbar[128];
   armbbbar[129]=3*armbbbar[8];
   armbbbar[130]=armbbbar[129] + 239./4.*armbbbar[11];
   armbbbar[131]= - 1 + armbbbar[95];
   armbbbar[131]=3./4.*armbbbar[131] - armbbbar[39];
   armbbbar[131]=1./2.*armbbbar[131] - armbbbar[42];
   armbbbar[131]=MMH*armbbbar[131];
   armbbbar[130]= - 353./32.*armbbbar[12] + 1./4.*armbbbar[130] + 
   armbbbar[131];
   armbbbar[130]=armbbbar[12]*armbbbar[130];
   armbbbar[131]= - 1 + armbbbar[78];
   armbbbar[131]=1./4.*armbbbar[131] + armbbbar[39];
   armbbbar[131]=MMH*armbbbar[131];
   armbbbar[131]= - 51./8.*armbbbar[13] + 79./8.*armbbbar[12] - 107./8.
   *armbbbar[11] + armbbbar[131];
   armbbbar[131]=armbbbar[13]*armbbbar[131];
   armbbbar[128]=1./2.*armbbbar[131] + 1./4.*armbbbar[128] + 
   armbbbar[130];
   armbbbar[88]=1./2.*armbbbar[100] + 1./2.*armbbbar[128] + 
   armbbbar[88];
   armbbbar[88]=armbbbar[52]*armbbbar[88];
   armbbbar[100]= - armbbbar[41] - 61./4. - armbbbar[50];
   armbbbar[100]=1./4.*armbbbar[100] + armbbbar[38];
   armbbbar[128]=3./4.*armbbbar[94];
   armbbbar[130]=armbbbar[12]*armbbbar[56];
   armbbbar[100]=3./4.*armbbbar[130] + armbbbar[128] + armbbbar[97] + 3
   *armbbbar[100] - 13./2.*armbbbar[47];
   armbbbar[131]= - armbbbar[3]*armbbbar[12];
   armbbbar[100]=9./2.*armbbbar[102] + 1./2.*armbbbar[100] + 9*
   armbbbar[131];
   armbbbar[100]=armbbbar[13]*armbbbar[100];
   armbbbar[101]=armbbbar[101] - armbbbar[35];
   armbbbar[101]= - 95./16.*armbbbar[16] + 33./16.*armbbbar[23] - 15./
   16.*armbbbar[24] - 59./16.*armbbbar[27] - 11./4.*armbbbar[32] - 17./
   8.*armbbbar[31] + 9./8.*armbbbar[101] + armbbbar[30];
   armbbbar[132]=23./4. + armbbbar[95];
   armbbbar[132]= - armbbbar[42] + 1./2.*armbbbar[132] - armbbbar[39];
   armbbbar[132]=MMH*armbbbar[132];
   armbbbar[76]=1./2.*armbbbar[88] + 1./2.*armbbbar[76] + 1./4.*
   armbbbar[100] + 9./16.*armbbbar[110] + 1./2.*armbbbar[85] + 1./8.*
   armbbbar[132] - 5./32.*armbbbar[11] - 9./64.*armbbbar[8] + 25./8.*
   armbbbar[25] + 1./4.*armbbbar[101] - 6*armbbbar[17];
   armbbbar[76]=armbbbar[52]*armbbbar[76];
   armbbbar[85]=armbbbar[63] + armbbbar[65];
   armbbbar[85]=3./4.*MMH*armbbbar[85];
   armbbbar[88]=165*armbbbar[43] - 165*armbbbar[45] + 3*armbbbar[41] - 
   319./4. + armbbbar[90];
   armbbbar[88]=1./2.*armbbbar[88] + 5*armbbbar[47];
   armbbbar[100]= - 3./2.*armbbbar[48];
   armbbbar[101]= - 3*armbbbar[42];
   armbbbar[110]= - armbbbar[11]*armbbbar[54];
   armbbbar[132]=3./4.*armbbbar[12]*armbbbar[75];
   armbbbar[88]=armbbbar[132] + armbbbar[85] + 7./4.*armbbbar[110] + 
   armbbbar[91] + armbbbar[101] + armbbbar[100] + 1./2.*armbbbar[88] + 
   armbbbar[81];
   armbbbar[133]= - 3*armbbbar[3];
   armbbbar[69]=armbbbar[69] + armbbbar[133];
   armbbbar[69]=armbbbar[13]*armbbbar[69];
   armbbbar[134]=79./4.*armbbbar[12] - armbbbar[8] + 97./8.*
   armbbbar[11];
   armbbbar[134]=armbbbar[3]*armbbbar[134];
   armbbbar[69]=armbbbar[69] + 1./2.*armbbbar[88] + armbbbar[134];
   armbbbar[69]=armbbbar[13]*armbbbar[69];
   armbbbar[88]=75./4.*armbbbar[24];
   armbbbar[134]= - 5./4.*armbbbar[23];
   armbbbar[135]=9*armbbbar[46];
   armbbbar[136]=19*armbbbar[47] + 167./8. + armbbbar[135];
   armbbbar[136]=armbbbar[8]*armbbbar[136];
   armbbbar[136]=armbbbar[136] + 91./2.*armbbbar[25] + 21./2.*
   armbbbar[17] + 33./4.*armbbbar[16] + armbbbar[134] + armbbbar[88] - 
   33./4.*armbbbar[27] - armbbbar[35] - 35./2.*armbbbar[31];
   armbbbar[137]=1./3.*armbbbar[11];
   armbbbar[138]= - 7./4.*armbbbar[38] - 127./12. + armbbbar[135];
   armbbbar[139]=armbbbar[138] + 1./3.*armbbbar[39];
   armbbbar[139]=MMH*armbbbar[139];
   armbbbar[140]=1./4. - armbbbar[48];
   armbbbar[140]=armbbbar[12]*armbbbar[140];
   armbbbar[141]=armbbbar[3]*armbbbar[12]*armbbbar[8];
   armbbbar[136]=1./4.*armbbbar[141] + 3./4.*armbbbar[140] + 1./8.*
   armbbbar[139] + 1./8.*armbbbar[136] + armbbbar[137];
   armbbbar[139]=27*armbbbar[37];
   armbbbar[140]=3./4.*armbbbar[50];
   armbbbar[142]=armbbbar[140] + 667./48. + armbbbar[139];
   armbbbar[143]=15./2.*armbbbar[38];
   armbbbar[142]=19./4.*armbbbar[47] + armbbbar[80] + armbbbar[143] + 1.
   /2.*armbbbar[142] + armbbbar[79];
   armbbbar[144]=5*armbbbar[8];
   armbbbar[145]=armbbbar[144] - 43*armbbbar[11];
   armbbbar[145]=1./2.*armbbbar[145] + 3*armbbbar[12];
   armbbbar[145]=armbbbar[3]*armbbbar[145];
   armbbbar[146]= - armbbbar[13]*armbbbar[3];
   armbbbar[142]=29./4.*armbbbar[146] + 1./4.*armbbbar[145] + 1./4.*
   armbbbar[142] - armbbbar[39];
   armbbbar[142]=armbbbar[13]*armbbbar[142];
   armbbbar[145]= - 9*armbbbar[49];
   armbbbar[147]=7./16.*armbbbar[50];
   armbbbar[148]=19./4.*armbbbar[38];
   armbbbar[149]=3./8.*armbbbar[48] + 19./24.*armbbbar[47] + 
   armbbbar[148] + armbbbar[147] + 1477./384. + armbbbar[145];
   armbbbar[150]=13*armbbbar[35] + 27*armbbbar[31];
   armbbbar[151]= - 1./2.*armbbbar[23];
   armbbbar[150]=armbbbar[151] + 1./2.*armbbbar[150] - 13*armbbbar[24];
   armbbbar[152]= - 3*armbbbar[17];
   armbbbar[153]= - 5./2. - 3*armbbbar[49];
   armbbbar[153]=armbbbar[8]*armbbbar[153];
   armbbbar[150]=3./2.*armbbbar[153] - 33./2.*armbbbar[25] + 1./2.*
   armbbbar[150] + armbbbar[152];
   armbbbar[153]= - 1./2. - armbbbar[48];
   armbbbar[153]=armbbbar[12]*armbbbar[153];
   armbbbar[154]= - armbbbar[11]*armbbbar[48];
   armbbbar[153]=3*armbbbar[153] + armbbbar[150] + 3./4.*armbbbar[154];
   armbbbar[153]=armbbbar[3]*armbbbar[153];
   armbbbar[149]=1./2.*armbbbar[149] + armbbbar[153];
   armbbbar[153]=armbbbar[77] - 155./8. + armbbbar[135];
   armbbbar[153]= - 19./2.*armbbbar[47] + 1./2.*armbbbar[153] - 6*
   armbbbar[38];
   armbbbar[153]=armbbbar[3]*armbbbar[153];
   armbbbar[155]=pow(armbbbar[3],2);
   armbbbar[156]= - armbbbar[13]*armbbbar[155];
   armbbbar[153]=armbbbar[153] + 9./2.*armbbbar[156];
   armbbbar[153]=armbbbar[13]*armbbbar[153];
   armbbbar[157]=3 - 1./2.*armbbbar[50];
   armbbbar[158]=1./2.*armbbbar[157] - armbbbar[38];
   armbbbar[158]=armbbbar[3]*armbbbar[158];
   armbbbar[159]= - armbbbar[13]*armbbbar[155]*armbbbar[49];
   armbbbar[158]=armbbbar[158] + 6*armbbbar[159];
   armbbbar[158]=3*MMt*armbbbar[158];
   armbbbar[149]=armbbbar[158] + 1./2.*armbbbar[149] + armbbbar[153];
   armbbbar[149]=MMt*armbbbar[149];
   armbbbar[136]=armbbbar[149] + 1./2.*armbbbar[136] + armbbbar[142];
   armbbbar[136]=MMt*armbbbar[136];
   armbbbar[142]=31 + armbbbar[139];
   armbbbar[142]=armbbbar[80] + 1./2.*armbbbar[142] + armbbbar[79];
   armbbbar[142]=1./4.*armbbbar[142] + armbbbar[39];
   armbbbar[142]=armbbbar[8]*armbbbar[142];
   armbbbar[117]= - 13./2.*armbbbar[31] + armbbbar[117] - 5./2.*
   armbbbar[29] - 15./2.*armbbbar[28] + armbbbar[34];
   armbbbar[149]=armbbbar[42] + 3./2. + armbbbar[39];
   armbbbar[149]=armbbbar[11]*armbbbar[149];
   armbbbar[117]=1./2.*armbbbar[149] + armbbbar[142] + 7./8.*
   armbbbar[25] + 9./8.*armbbbar[17] + armbbbar[16] + 19./8.*
   armbbbar[23] + 1./4.*armbbbar[117] + 7*armbbbar[24];
   armbbbar[142]= - 1./32.*armbbbar[42] + 3./32.*armbbbar[40] + 3./16.*
   armbbbar[44] - 1 + 27./32.*armbbbar[37];
   armbbbar[142]=MMH*armbbbar[142];
   armbbbar[117]=1./4.*armbbbar[117] + armbbbar[142];
   armbbbar[117]=MMH*armbbbar[117];
   armbbbar[142]=pow(armbbbar[8],2);
   armbbbar[149]=3./8.*armbbbar[142];
   armbbbar[153]=13*armbbbar[8] + armbbbar[99];
   armbbbar[153]=armbbbar[11]*armbbbar[153];
   armbbbar[153]=armbbbar[149] + armbbbar[153];
   armbbbar[160]=7*armbbbar[8];
   armbbbar[161]=armbbbar[160] + 41*armbbbar[11];
   armbbbar[162]=13./16. + armbbbar[42];
   armbbbar[162]=MMH*armbbbar[162];
   armbbbar[161]= - 1./8.*armbbbar[12] + 1./16.*armbbbar[161] + 
   armbbbar[162];
   armbbbar[161]=armbbbar[12]*armbbbar[161];
   armbbbar[162]=23*armbbbar[8] + 7./2.*armbbbar[11];
   armbbbar[163]=1 + armbbbar[78];
   armbbbar[164]=1./4.*armbbbar[163];
   armbbbar[165]=armbbbar[164] + armbbbar[39];
   armbbbar[165]=MMH*armbbbar[165];
   armbbbar[162]=1./4.*armbbbar[162] + armbbbar[165];
   armbbbar[165]= - 5*armbbbar[12];
   armbbbar[162]=215./32.*armbbbar[13] + 1./2.*armbbbar[162] + 
   armbbbar[165];
   armbbbar[162]=armbbbar[13]*armbbbar[162];
   armbbbar[117]=1./2.*armbbbar[162] + 1./2.*armbbbar[161] + 1./8.*
   armbbbar[153] + armbbbar[117];
   armbbbar[153]= - 3./2.*armbbbar[23];
   armbbbar[119]=armbbbar[121] + armbbbar[119] + armbbbar[153] + 
   armbbbar[30] + 1./2.*armbbbar[31];
   armbbbar[119]=armbbbar[89]*armbbbar[119];
   armbbbar[161]= - MMH*armbbbar[8];
   armbbbar[114]=armbbbar[114] + armbbbar[161] + 3*armbbbar[116];
   armbbbar[114]=armbbbar[13]*armbbbar[114];
   armbbbar[162]=MMH*armbbbar[42];
   armbbbar[166]=3./2.*armbbbar[8] + armbbbar[162];
   armbbbar[166]=MMH*armbbbar[166];
   armbbbar[166]=1./2.*armbbbar[166] + armbbbar[120];
   armbbbar[166]=armbbbar[12]*armbbbar[166];
   armbbbar[166]=1./4.*armbbbar[114] + 1./2.*armbbbar[119] + 
   armbbbar[166];
   armbbbar[166]=armbbbar[52]*armbbbar[166];
   armbbbar[167]= - 3*armbbbar[8] - 11./2.*armbbbar[11];
   armbbbar[168]=armbbbar[39] - armbbbar[42];
   armbbbar[169]=MMH*armbbbar[168];
   armbbbar[170]= - 3./4.*armbbbar[13];
   armbbbar[167]=armbbbar[170] + armbbbar[12] + 1./4.*armbbbar[167] + 
   armbbbar[169];
   armbbbar[167]=armbbbar[13]*armbbbar[167];
   armbbbar[171]=3*armbbbar[11] - 1./6.*MMH;
   armbbbar[171]=armbbbar[12]*armbbbar[171];
   armbbbar[172]= - armbbbar[39] + armbbbar[42];
   armbbbar[173]=MMH*armbbbar[172];
   armbbbar[174]=armbbbar[11] + armbbbar[173];
   armbbbar[174]=MMH*armbbbar[174];
   armbbbar[171]=1./6.*armbbbar[174] + armbbbar[171];
   armbbbar[167]=1./8.*armbbbar[171] + armbbbar[167];
   armbbbar[171]=armbbbar[16] + armbbbar[35] - armbbbar[27];
   armbbbar[175]= - armbbbar[8]*armbbbar[48];
   armbbbar[171]=1./2.*armbbbar[171] + armbbbar[175];
   armbbbar[176]= - 1./2. + armbbbar[92];
   armbbbar[176]=1./8.*armbbbar[176] + 3*armbbbar[102];
   armbbbar[176]=armbbbar[13]*armbbbar[176];
   armbbbar[177]=armbbbar[3]*armbbbar[48];
   armbbbar[178]=armbbbar[13]*armbbbar[177];
   armbbbar[178]= - 1./16.*armbbbar[48] + 3*armbbbar[178];
   armbbbar[178]=MMt*armbbbar[178];
   armbbbar[171]=armbbbar[178] + 3./8.*armbbbar[171] + armbbbar[176];
   armbbbar[171]=MMt*armbbbar[171];
   armbbbar[167]=1./2.*armbbbar[167] + armbbbar[171];
   armbbbar[167]=MMt*armbbbar[167];
   armbbbar[171]=armbbbar[11]*armbbbar[8];
   armbbbar[179]=MMH*armbbbar[8]*armbbbar[172];
   armbbbar[179]=armbbbar[171] + armbbbar[179];
   armbbbar[179]=MMH*armbbbar[179];
   armbbbar[180]=armbbbar[174] + armbbbar[120];
   armbbbar[181]=armbbbar[13]*armbbbar[180];
   armbbbar[161]=armbbbar[12]*armbbbar[161];
   armbbbar[161]=armbbbar[181] + armbbbar[179] + armbbbar[161];
   armbbbar[161]=1./16.*armbbbar[161];
   armbbbar[167]=armbbbar[161] + armbbbar[167];
   armbbbar[167]=armbbbar[4]*armbbbar[167];
   armbbbar[117]=1./2.*armbbbar[167] + 1./8.*armbbbar[166] + 1./2.*
   armbbbar[117] + armbbbar[136];
   armbbbar[117]=armbbbar[4]*armbbbar[117];
   armbbbar[136]= - 3*armbbbar[44];
   armbbbar[166]= - 3./2.*armbbbar[40];
   armbbbar[167]=armbbbar[91] + armbbbar[101] + armbbbar[92] + 
   armbbbar[81] + 5./2.*armbbbar[47] + 165./4.*armbbbar[43] + 
   armbbbar[166] - 165./4.*armbbbar[45] - 397./16. + armbbbar[136];
   armbbbar[167]=armbbbar[8]*armbbbar[167];
   armbbbar[179]=1./4. - 3*armbbbar[37];
   armbbbar[181]= - 1./2.*armbbbar[40];
   armbbbar[179]=armbbbar[181] + 3./2.*armbbbar[179] - armbbbar[44];
   armbbbar[182]= - armbbbar[8]*armbbbar[57];
   armbbbar[183]=armbbbar[11]*armbbbar[57];
   armbbbar[179]=3./8.*armbbbar[183] + 3./8.*armbbbar[182] + 3./4.*
   armbbbar[179] - armbbbar[39];
   armbbbar[179]=armbbbar[11]*armbbbar[179];
   armbbbar[66]=armbbbar[22]*armbbbar[66];
   armbbbar[184]= - armbbbar[8]*armbbbar[22];
   armbbbar[185]=armbbbar[66] + armbbbar[184];
   armbbbar[185]=1./2.*armbbbar[185] + armbbbar[75];
   armbbbar[185]=MMH*armbbbar[185];
   armbbbar[186]=1./2.*armbbbar[39] + armbbbar[181] + 5./8. - 
   armbbbar[44];
   armbbbar[186]=3./4.*armbbbar[185] + 3./2.*armbbbar[186] + 
   armbbbar[42];
   armbbbar[186]=MMH*armbbbar[186];
   armbbbar[187]= - 3*armbbbar[35] + 13*armbbbar[29] + armbbbar[34] - 3
   *armbbbar[6];
   armbbbar[187]=1./2.*armbbbar[187] + 11*armbbbar[30];
   armbbbar[187]=15./32.*armbbbar[32] + 1./4.*armbbbar[187] - 
   armbbbar[31];
   armbbbar[167]=1./2.*armbbbar[186] + 1./2.*armbbbar[179] + 1./4.*
   armbbbar[167] + 15./16.*armbbbar[25] - 75./32.*armbbbar[17] - 47./64.
   *armbbbar[16] + 25./32.*armbbbar[23] - 81./32.*armbbbar[24] + 1./2.*
   armbbbar[187] - armbbbar[27];
   armbbbar[179]= - 27*armbbbar[37];
   armbbbar[186]= - 3./2.*armbbbar[50];
   armbbbar[187]=armbbbar[186] + 259./8. + armbbbar[179];
   armbbbar[187]=armbbbar[166] + armbbbar[78] + 1./2.*armbbbar[187] + 
   armbbbar[136];
   armbbbar[188]=armbbbar[184] + 1./2.*armbbbar[65];
   armbbbar[188]=3./16.*MMH*armbbbar[188];
   armbbbar[189]= - 15./16.*armbbbar[48];
   armbbbar[190]= - 1./2.*armbbbar[42];
   armbbbar[191]=3./16.*armbbbar[94];
   armbbbar[192]=3./16.*armbbbar[130];
   armbbbar[187]=armbbbar[192] + armbbbar[188] + armbbbar[191] + 
   armbbbar[190] + armbbbar[189] + 1./8.*armbbbar[187] + armbbbar[47];
   armbbbar[187]=armbbbar[12]*armbbbar[187];
   armbbbar[193]=111./8. - armbbbar[50];
   armbbbar[96]=armbbbar[96] + armbbbar[181] + 55./4.*armbbbar[45] - 1./
   2.*armbbbar[41] + 1./2.*armbbbar[193] - armbbbar[44];
   armbbbar[193]=3*armbbbar[39];
   armbbbar[109]=armbbbar[3]*armbbbar[109];
   armbbbar[194]=3*armbbbar[42];
   armbbbar[96]=9./2.*armbbbar[109] + armbbbar[83] + armbbbar[194] + 
   armbbbar[97] + armbbbar[193] + 3*armbbbar[96] - 5./2.*armbbbar[47];
   armbbbar[96]=armbbbar[3]*armbbbar[96];
   armbbbar[109]= - armbbbar[11]*armbbbar[26];
   armbbbar[195]=3./8.*armbbbar[75];
   armbbbar[196]= - armbbbar[26] - 9./8.*armbbbar[22];
   armbbbar[196]=armbbbar[12]*armbbbar[196];
   armbbbar[196]=armbbbar[196] + armbbbar[195] + 3./8.*armbbbar[184] + 
   armbbbar[109];
   armbbbar[96]=1./4.*armbbbar[196] + armbbbar[96];
   armbbbar[96]=armbbbar[13]*armbbbar[96];
   armbbbar[64]=1./2.*armbbbar[64] + armbbbar[65];
   armbbbar[64]=3*MMH*armbbbar[64];
   armbbbar[196]=7*armbbbar[41];
   armbbbar[197]=55*armbbbar[43] - 55*armbbbar[45] + armbbbar[196] - 
   179 + 7*armbbbar[50];
   armbbbar[197]=1./2.*armbbbar[197] - 43./3.*armbbbar[47];
   armbbbar[198]=1./6.*armbbbar[9];
   armbbbar[197]=armbbbar[64] + armbbbar[198] - armbbbar[42] - 7*
   armbbbar[48] + 1./2.*armbbbar[197] - armbbbar[39];
   armbbbar[199]=armbbbar[11]*armbbbar[26];
   armbbbar[200]=3./16.*armbbbar[75] + 3./4.*armbbbar[63] + 
   armbbbar[199];
   armbbbar[200]=armbbbar[12]*armbbbar[200];
   armbbbar[197]=1./8.*armbbbar[197] + armbbbar[200];
   armbbbar[200]= - 9*armbbbar[46];
   armbbbar[201]=59./8. + armbbbar[200];
   armbbbar[202]=4*armbbbar[47];
   armbbbar[203]= - 3./4.*armbbbar[48];
   armbbbar[204]=armbbbar[203] + 1./4.*armbbbar[201] + armbbbar[202];
   armbbbar[204]=armbbbar[12]*armbbbar[204];
   armbbbar[205]= - 3*armbbbar[16] + armbbbar[151] + 1./2.*armbbbar[24]
    + 13./2.*armbbbar[35] + 3*armbbbar[32];
   armbbbar[205]=1./2.*armbbbar[205] + armbbbar[152];
   armbbbar[201]=armbbbar[201] + 13*armbbbar[47];
   armbbbar[201]=armbbbar[11]*armbbbar[201];
   armbbbar[201]=armbbbar[204] + 1./8.*armbbbar[201] + 1./2.*
   armbbbar[205] - 3*armbbbar[25];
   armbbbar[201]=armbbbar[3]*armbbbar[201];
   armbbbar[204]=armbbbar[11]*armbbbar[49];
   armbbbar[205]=armbbbar[12]*armbbbar[49];
   armbbbar[206]=1./2.*armbbbar[204] + armbbbar[205];
   armbbbar[206]=armbbbar[3]*armbbbar[206];
   armbbbar[207]= - armbbbar[41] + 7 - armbbbar[50];
   armbbbar[207]=1./2.*armbbbar[207] - armbbbar[48];
   armbbbar[206]=1./2.*armbbbar[207] + 3*armbbbar[206];
   armbbbar[206]=MMt*armbbbar[3]*armbbbar[206];
   armbbbar[96]=3*armbbbar[206] + armbbbar[96] + 1./2.*armbbbar[197] + 
   armbbbar[201];
   armbbbar[96]=MMt*armbbbar[96];
   armbbbar[105]=armbbbar[106] + 1./2.*armbbbar[17] + armbbbar[153] + 
   armbbbar[105] + armbbbar[31];
   armbbbar[89]=1./2.*armbbbar[89]*armbbbar[105];
   armbbbar[153]= - armbbbar[8] + armbbbar[11];
   armbbbar[197]= - armbbbar[39] - armbbbar[42];
   armbbbar[197]=MMH*armbbbar[197];
   armbbbar[197]=armbbbar[153] + armbbbar[197];
   armbbbar[197]=MMH*armbbbar[197];
   armbbbar[197]=armbbbar[115] + 1./4.*armbbbar[197] + armbbbar[116];
   armbbbar[197]=armbbbar[13]*armbbbar[197];
   armbbbar[201]=armbbbar[39] + armbbbar[194];
   armbbbar[201]=MMH*armbbbar[201];
   armbbbar[129]=armbbbar[201] + armbbbar[129] - armbbbar[11];
   armbbbar[129]=MMH*armbbbar[129];
   armbbbar[129]=1./4.*armbbbar[129] + armbbbar[120];
   armbbbar[129]=armbbbar[12]*armbbbar[129];
   armbbbar[129]=armbbbar[197] + armbbbar[89] + armbbbar[129];
   armbbbar[129]=armbbbar[52]*armbbbar[129];
   armbbbar[197]= - 5*armbbbar[31];
   armbbbar[201]=9*armbbbar[23] - 1./2.*armbbbar[24] - 7./2.*
   armbbbar[30] + armbbbar[197];
   armbbbar[201]=armbbbar[121] + 1./2.*armbbbar[201] - armbbbar[17];
   armbbbar[206]=1./4.*armbbbar[8];
   armbbbar[207]=1 + armbbbar[39];
   armbbbar[207]=1./2.*armbbbar[207] + armbbbar[42];
   armbbbar[207]=MMH*armbbbar[207];
   armbbbar[207]=1./2.*armbbbar[207] - 1./4.*armbbbar[11] + 
   armbbbar[201] + armbbbar[206];
   armbbbar[207]=MMH*armbbbar[207];
   armbbbar[208]=13./4.*armbbbar[8] + 59*armbbbar[11];
   armbbbar[209]= - 1 + armbbbar[38];
   armbbbar[209]=3./4.*armbbbar[209] - armbbbar[42];
   armbbbar[209]=MMH*armbbbar[209];
   armbbbar[208]= - 65./4.*armbbbar[12] + 1./2.*armbbbar[208] + 
   armbbbar[209];
   armbbbar[208]=armbbbar[12]*armbbbar[208];
   armbbbar[210]=1./2.*armbbbar[8];
   armbbbar[211]=7*armbbbar[13] - 105./4.*armbbbar[12] + armbbbar[210]
    - 17*armbbbar[11];
   armbbbar[211]=armbbbar[13]*armbbbar[211];
   armbbbar[129]=armbbbar[129] + 1./2.*armbbbar[211] + armbbbar[207] + 
   armbbbar[208];
   armbbbar[129]=armbbbar[52]*armbbbar[129];
   armbbbar[207]= - armbbbar[8] + 41./2.*armbbbar[11];
   armbbbar[207]=armbbbar[11]*armbbbar[207];
   armbbbar[208]=5./2.*armbbbar[11];
   armbbbar[211]= - 81./2.*armbbbar[12] - armbbbar[8] + armbbbar[208];
   armbbbar[211]=armbbbar[12]*armbbbar[211];
   armbbbar[207]=1./2.*armbbbar[207] + armbbbar[211];
   armbbbar[207]=armbbbar[3]*armbbbar[207];
   armbbbar[69]=armbbbar[117] + 1./8.*armbbbar[129] + armbbbar[96] + 1./
   4.*armbbbar[69] + 1./8.*armbbbar[207] + 1./2.*armbbbar[167] + 
   armbbbar[187];
   armbbbar[69]=armbbbar[4]*armbbbar[69];
   armbbbar[96]=143./8. + armbbbar[200];
   armbbbar[117]=3./2.*armbbbar[41];
   armbbbar[96]=armbbbar[100] + 23./4.*armbbbar[47] + armbbbar[117] + 1.
   /2.*armbbbar[96] + armbbbar[90];
   armbbbar[96]=armbbbar[3]*armbbbar[96];
   armbbbar[129]= - armbbbar[27] + armbbbar[16];
   armbbbar[167]=armbbbar[25] + armbbbar[129] + armbbbar[17];
   armbbbar[167]=armbbbar[26]*armbbbar[167];
   armbbbar[187]=3./8.*armbbbar[63];
   armbbbar[207]=armbbbar[26] + 3./8.*armbbbar[22];
   armbbbar[211]=armbbbar[12]*armbbbar[207];
   armbbbar[96]=armbbbar[96] + armbbbar[211] + 21./8.*armbbbar[65] + 
   armbbbar[199] + armbbbar[187] + armbbbar[167] + 3./8.*armbbbar[62];
   armbbbar[211]= - armbbbar[22]*armbbbar[56];
   armbbbar[212]=armbbbar[12]*armbbbar[211];
   armbbbar[213]=armbbbar[22]*armbbbar[56];
   armbbbar[214]=armbbbar[8]*armbbbar[213];
   armbbbar[207]=3./16.*armbbbar[212] + armbbbar[207] + 3./16.*
   armbbbar[214];
   armbbbar[215]=1./2.*armbbbar[40];
   armbbbar[216]=armbbbar[215] - 1./2. + armbbbar[44];
   armbbbar[216]=armbbbar[155]*armbbbar[216];
   armbbbar[207]=1./2.*armbbbar[207] + 9*armbbbar[216];
   armbbbar[207]=armbbbar[13]*armbbbar[207];
   armbbbar[216]=MMt*armbbbar[155]*armbbbar[49];
   armbbbar[96]=9*armbbbar[216] + 1./2.*armbbbar[96] + armbbbar[207];
   armbbbar[96]=MMt*armbbbar[96];
   armbbbar[207]= - 3*armbbbar[47];
   armbbbar[217]=armbbbar[207] - 99./2.*armbbbar[43] - 17 + 99./2.*
   armbbbar[45];
   armbbbar[218]= - 3./2.*armbbbar[9];
   armbbbar[217]=armbbbar[218] + armbbbar[194] + armbbbar[97] + 1./2.*
   armbbbar[217] + armbbbar[81];
   armbbbar[217]=armbbbar[11]*armbbbar[217];
   armbbbar[219]= - 99*armbbbar[43] + 17 + 99*armbbbar[45];
   armbbbar[219]=1./2.*armbbbar[219];
   armbbbar[207]=armbbbar[219] + armbbbar[207];
   armbbbar[207]=armbbbar[218] + armbbbar[194] + armbbbar[97] + 1./2.*
   armbbbar[207] + armbbbar[81];
   armbbbar[207]=armbbbar[12]*armbbbar[207];
   armbbbar[207]=1./2.*armbbbar[217] + armbbbar[207];
   armbbbar[207]=armbbbar[3]*armbbbar[207];
   armbbbar[217]=3*armbbbar[47];
   armbbbar[219]=armbbbar[219] + armbbbar[217];
   armbbbar[218]=armbbbar[218] + armbbbar[194] - 9./2.*armbbbar[48] + 1.
   /2.*armbbbar[219] + armbbbar[81];
   armbbbar[218]=armbbbar[12]*armbbbar[218];
   armbbbar[219]=armbbbar[45] - armbbbar[43];
   armbbbar[220]=armbbbar[83] + armbbbar[42] + 33./4.*armbbbar[219] - 
   armbbbar[39];
   armbbbar[221]=armbbbar[12]*armbbbar[220];
   armbbbar[222]= - armbbbar[45] + armbbbar[43];
   armbbbar[223]=armbbbar[91] - armbbbar[42] + 33./4.*armbbbar[222] + 
   armbbbar[39];
   armbbbar[223]=armbbbar[13]*armbbbar[223];
   armbbbar[221]=armbbbar[221] + armbbbar[223];
   armbbbar[221]=MMZ*armbbbar[221];
   armbbbar[223]= - 17./4.*armbbbar[11];
   armbbbar[224]=armbbbar[223] + armbbbar[169];
   armbbbar[225]=17./4.*armbbbar[12];
   armbbbar[226]=armbbbar[224] + armbbbar[225];
   armbbbar[226]=armbbbar[12]*armbbbar[226];
   armbbbar[227]=3./2.*armbbbar[13] - 23./4.*armbbbar[12] + 17./4.*
   armbbbar[11] + armbbbar[173];
   armbbbar[227]=armbbbar[13]*armbbbar[227];
   armbbbar[221]=3*armbbbar[221] + armbbbar[226] + armbbbar[227];
   armbbbar[221]=armbbbar[52]*armbbbar[221];
   armbbbar[220]=MMZ*armbbbar[220];
   armbbbar[226]=armbbbar[97] - 1 - armbbbar[47];
   armbbbar[226]=armbbbar[13]*armbbbar[226];
   armbbbar[218]=armbbbar[221] + 3*armbbbar[220] + 3./2.*armbbbar[226]
    + armbbbar[224] + armbbbar[218];
   armbbbar[218]=armbbbar[52]*armbbbar[218];
   armbbbar[220]=11./2.*armbbbar[222] + armbbbar[47];
   armbbbar[221]= - 9*armbbbar[48];
   armbbbar[220]=3./2.*armbbbar[9] - 11*armbbbar[42] + armbbbar[221] + 
   9./2.*armbbbar[220] + 11*armbbbar[39];
   armbbbar[222]=33./2.*armbbbar[219] - armbbbar[47];
   armbbbar[222]=armbbbar[83] + armbbbar[42] + armbbbar[48] + 1./2.*
   armbbbar[222] - armbbbar[39];
   armbbbar[226]=MMZ*armbbbar[3]*armbbbar[222];
   armbbbar[227]=armbbbar[47] - armbbbar[48];
   armbbbar[228]=MMt*armbbbar[3]*armbbbar[227];
   armbbbar[207]=1./4.*armbbbar[218] + 3*armbbbar[226] + 3./2.*
   armbbbar[228] + 1./8.*armbbbar[220] + armbbbar[207];
   armbbbar[168]=armbbbar[8]*armbbbar[168];
   armbbbar[218]=armbbbar[42] - 1./4. - armbbbar[39];
   armbbbar[218]=armbbbar[11]*armbbbar[218];
   armbbbar[168]=1./8.*armbbbar[169] + armbbbar[168] + 1./2.*
   armbbbar[218];
   armbbbar[168]=MMH*armbbbar[168];
   armbbbar[218]=17./2.*armbbbar[8] + armbbbar[11];
   armbbbar[220]=armbbbar[42] + 1./8. - armbbbar[39];
   armbbbar[220]=MMH*armbbbar[220];
   armbbbar[218]= - armbbbar[12] + 1./2.*armbbbar[218] + armbbbar[220];
   armbbbar[218]=armbbbar[12]*armbbbar[218];
   armbbbar[220]=armbbbar[224] + 11./4.*armbbbar[12];
   armbbbar[220]=armbbbar[13]*armbbbar[220];
   armbbbar[224]= - 17./2.*armbbbar[8] + armbbbar[11];
   armbbbar[224]=armbbbar[11]*armbbbar[224];
   armbbbar[168]=armbbbar[220] + armbbbar[218] + 1./2.*armbbbar[224] + 
   armbbbar[168];
   armbbbar[218]=armbbbar[8]*armbbbar[227];
   armbbbar[220]=17./12. + armbbbar[92];
   armbbbar[220]=armbbbar[12]*armbbbar[220];
   armbbbar[218]=armbbbar[220] + 1./3.*armbbbar[169] + 3*armbbbar[218]
    - 17./12.*armbbbar[11];
   armbbbar[220]=1./2. + armbbbar[47];
   armbbbar[165]=31./8.*armbbbar[11] + armbbbar[165];
   armbbbar[165]=armbbbar[3]*armbbbar[165];
   armbbbar[165]=armbbbar[165] + armbbbar[42] - 3./16.*armbbbar[48] + 3.
   /16.*armbbbar[220] - armbbbar[39];
   armbbbar[165]=armbbbar[13]*armbbbar[165];
   armbbbar[220]= - armbbbar[12]*armbbbar[48];
   armbbbar[220]=1./2.*armbbbar[154] + armbbbar[220];
   armbbbar[220]=armbbbar[3]*armbbbar[220];
   armbbbar[224]=1./2.*armbbbar[47] + armbbbar[48];
   armbbbar[220]=1./4.*armbbbar[224] + 3*armbbbar[220];
   armbbbar[224]= - armbbbar[47] + armbbbar[48];
   armbbbar[224]=armbbbar[13]*armbbbar[3]*armbbbar[224];
   armbbbar[220]=1./2.*armbbbar[220] + 3*armbbbar[224];
   armbbbar[220]=MMt*armbbbar[220];
   armbbbar[165]=1./2.*armbbbar[220] + 1./16.*armbbbar[218] + 
   armbbbar[165];
   armbbbar[165]=MMt*armbbbar[165];
   armbbbar[170]=armbbbar[170] + armbbbar[12] + armbbbar[169] + 
   armbbbar[111] - armbbbar[11];
   armbbbar[170]=armbbbar[13]*armbbbar[170];
   armbbbar[170]=1./48.*armbbbar[180] + armbbbar[170];
   armbbbar[175]=armbbbar[178] + 3./8.*armbbbar[175] + armbbbar[176];
   armbbbar[175]=MMt*armbbbar[175];
   armbbbar[170]=1./2.*armbbbar[170] + armbbbar[175];
   armbbbar[170]=MMt*armbbbar[170];
   armbbbar[161]=armbbbar[161] + armbbbar[170];
   armbbbar[161]=armbbbar[4]*armbbbar[161];
   armbbbar[170]=armbbbar[12]*armbbbar[180];
   armbbbar[175]=armbbbar[52]*armbbbar[170];
   armbbbar[161]=1./2.*armbbbar[161] + 1./32.*armbbbar[175] + 1./8.*
   armbbbar[168] + armbbbar[165];
   armbbbar[161]=armbbbar[4]*armbbbar[161];
   armbbbar[165]=armbbbar[8]*armbbbar[222];
   armbbbar[168]= - armbbbar[42] + 25./16. + armbbbar[39];
   armbbbar[168]=armbbbar[11]*armbbbar[168];
   armbbbar[165]=3./4.*armbbbar[173] + 3./2.*armbbbar[165] + 
   armbbbar[168];
   armbbbar[82]= - armbbbar[47] - 33./2.*armbbbar[43] - 1 + 
   armbbbar[82];
   armbbbar[168]=3./2.*armbbbar[48];
   armbbbar[82]=armbbbar[83] + armbbbar[42] + armbbbar[168] + 1./2.*
   armbbbar[82] - armbbbar[39];
   armbbbar[82]=armbbbar[13]*armbbbar[82];
   armbbbar[173]=pow(armbbbar[11],2);
   armbbbar[175]= - 1./2.*armbbbar[173];
   armbbbar[176]=armbbbar[99] + armbbbar[12];
   armbbbar[176]=armbbbar[12]*armbbbar[176];
   armbbbar[176]=armbbbar[175] + armbbbar[176];
   armbbbar[176]=armbbbar[3]*armbbbar[176];
   armbbbar[178]= - 25./4. + armbbbar[217];
   armbbbar[178]= - armbbbar[42] - 9./8.*armbbbar[48] + 1./8.*
   armbbbar[178] + armbbbar[39];
   armbbbar[178]=armbbbar[12]*armbbbar[178];
   armbbbar[82]=3./4.*armbbbar[82] + 17./4.*armbbbar[176] + 1./2.*
   armbbbar[165] + armbbbar[178];
   armbbbar[165]= - armbbbar[42] - 1./4. + armbbbar[39];
   armbbbar[165]=MMH*armbbbar[165];
   armbbbar[165]=armbbbar[225] + armbbbar[223] + armbbbar[165];
   armbbbar[165]=armbbbar[12]*armbbbar[165];
   armbbbar[169]= - armbbbar[11] + armbbbar[169];
   armbbbar[169]=MMH*armbbbar[169];
   armbbbar[169]=armbbbar[169] + armbbbar[116];
   armbbbar[169]=armbbbar[13]*armbbbar[169];
   armbbbar[169]=armbbbar[170] + armbbbar[169];
   armbbbar[169]=armbbbar[52]*armbbbar[169];
   armbbbar[170]= - armbbbar[12] + 1./2.*armbbbar[13];
   armbbbar[170]=armbbbar[13]*armbbbar[170];
   armbbbar[165]=1./4.*armbbbar[169] + 3*armbbbar[170] + 1./4.*
   armbbbar[174] + armbbbar[165];
   armbbbar[165]=armbbbar[52]*armbbbar[165];
   armbbbar[169]= - 1./4.*armbbbar[9];
   armbbbar[174]=armbbbar[169] + 1./2.*armbbbar[42] - 1./4.*
   armbbbar[48] - 1./2.*armbbbar[39] + 33./8.*armbbbar[219] - 
   armbbbar[47];
   armbbbar[176]=armbbbar[11]*armbbbar[227];
   armbbbar[178]=armbbbar[12]*armbbbar[227];
   armbbbar[176]=1./2.*armbbbar[176] + armbbbar[178];
   armbbbar[176]=armbbbar[3]*armbbbar[176];
   armbbbar[174]=1./2.*armbbbar[174] + 3*armbbbar[176];
   armbbbar[176]=33*armbbbar[43];
   armbbbar[178]=armbbbar[176] - 1 - 33*armbbbar[45];
   armbbbar[178]=1./2.*armbbbar[178] + armbbbar[47];
   armbbbar[91]=armbbbar[91] - armbbbar[42] - armbbbar[48] + 1./2.*
   armbbbar[178] + armbbbar[39];
   armbbbar[91]=armbbbar[13]*armbbbar[3]*armbbbar[91];
   armbbbar[178]= - armbbbar[3]*armbbbar[48];
   armbbbar[180]=MMt*armbbbar[178];
   armbbbar[91]=3./4.*armbbbar[180] + 1./4.*armbbbar[174] + 3*
   armbbbar[91];
   armbbbar[91]=MMt*armbbbar[91];
   armbbbar[82]=armbbbar[161] + 1./8.*armbbbar[165] + 1./2.*
   armbbbar[82] + armbbbar[91];
   armbbbar[82]=armbbbar[4]*armbbbar[82];
   armbbbar[82]=1./2.*armbbbar[207] + armbbbar[82];
   armbbbar[91]=pow(armbbbar[5],2);
   armbbbar[82]=armbbbar[91]*armbbbar[82];
   armbbbar[161]=armbbbar[22]*armbbbar[73];
   armbbbar[165]=armbbbar[8]*armbbbar[161];
   armbbbar[174]=3*armbbbar[165];
   armbbbar[180]=7./4.*armbbbar[211] + armbbbar[174];
   armbbbar[207]=armbbbar[12]*armbbbar[180];
   armbbbar[217]=3 + 7./4.*armbbbar[70];
   armbbbar[217]=armbbbar[22]*armbbbar[217];
   armbbbar[218]=armbbbar[8]*armbbbar[211];
   armbbbar[207]=armbbbar[207] + armbbbar[217] + 7./4.*armbbbar[218];
   armbbbar[219]=armbbbar[181] + 3./4. - armbbbar[44];
   armbbbar[219]=armbbbar[155]*armbbbar[219];
   armbbbar[220]= - armbbbar[22]*armbbbar[73];
   armbbbar[222]=armbbbar[12]*armbbbar[220];
   armbbbar[180]=1./2.*armbbbar[180] + 3*armbbbar[222];
   armbbbar[180]=armbbbar[13]*armbbbar[180];
   armbbbar[180]=1./2.*armbbbar[180] + 1./4.*armbbbar[207] + 
   armbbbar[219];
   armbbbar[61]=armbbbar[56]*armbbbar[61];
   armbbbar[207]= - 9./4. + armbbbar[61];
   armbbbar[207]=armbbbar[22]*armbbbar[56]*armbbbar[207];
   armbbbar[219]=armbbbar[12]*armbbbar[161];
   armbbbar[223]=armbbbar[13]*armbbbar[161];
   armbbbar[224]=armbbbar[223] + armbbbar[219] + armbbbar[207] + 
   armbbbar[165];
   armbbbar[224]=MMt*armbbbar[224];
   armbbbar[226]=13./4. + armbbbar[70];
   armbbbar[226]=armbbbar[22]*armbbbar[56]*armbbbar[226];
   armbbbar[227]=armbbbar[8]*armbbbar[220];
   armbbbar[228]=armbbbar[13]*armbbbar[220];
   armbbbar[226]=armbbbar[228] + armbbbar[222] + armbbbar[226] + 
   armbbbar[227];
   armbbbar[229]=MMt*armbbbar[220];
   armbbbar[230]=MMZ*armbbbar[161];
   armbbbar[226]=9*armbbbar[230] + 3./2.*armbbbar[226] + armbbbar[229];
   armbbbar[226]=MMZ*armbbbar[226];
   armbbbar[180]=3*armbbbar[226] + 3*armbbbar[180] + 1./2.*
   armbbbar[224];
   armbbbar[180]=MMZ*armbbbar[180];
   armbbbar[224]= - 3*armbbbar[57] - 7./16.*armbbbar[54];
   armbbbar[226]= - 1./2.*armbbbar[22] + armbbbar[214];
   armbbbar[226]=armbbbar[12]*armbbbar[226];
   armbbbar[65]=3./2.*armbbbar[226] + armbbbar[65] + 3./4.*
   armbbbar[184] + 3./4.*armbbbar[66] + 1./2.*armbbbar[224] - 3*
   armbbbar[56];
   armbbbar[66]=1./2. - armbbbar[44];
   armbbbar[224]=armbbbar[66] + armbbbar[181];
   armbbbar[226]=armbbbar[11]*armbbbar[224];
   armbbbar[224]=armbbbar[12]*armbbbar[224];
   armbbbar[224]=1./2.*armbbbar[226] + armbbbar[224];
   armbbbar[224]=armbbbar[3]*armbbbar[224];
   armbbbar[224]=9*armbbbar[224] + 17./4.*armbbbar[9] + 9./2.*
   armbbbar[39] + 25./4.*armbbbar[47] + 33./2.*armbbbar[43] + 
   armbbbar[80] - 825./8.*armbbbar[45] + 661./48. + armbbbar[79];
   armbbbar[224]=armbbbar[3]*armbbbar[224];
   armbbbar[212]=armbbbar[212] - armbbbar[22] + armbbbar[218];
   armbbbar[212]=armbbbar[13]*armbbbar[212];
   armbbbar[65]=9./16.*armbbbar[212] + 3./4.*armbbbar[65] + 
   armbbbar[224];
   armbbbar[212]=1./4.*armbbbar[213] + armbbbar[227];
   armbbbar[224]=armbbbar[12]*armbbbar[212];
   armbbbar[212]=1./2.*armbbbar[212] + armbbbar[219];
   armbbbar[212]=armbbbar[13]*armbbbar[212];
   armbbbar[226]= - 3 + 1./2.*armbbbar[61];
   armbbbar[226]=armbbbar[22]*armbbbar[226];
   armbbbar[212]=1./4.*armbbbar[212] + 1./8.*armbbbar[224] + 1./32.*
   armbbbar[214] - armbbbar[26] + 1./16.*armbbbar[226];
   armbbbar[212]=MMt*armbbbar[212];
   armbbbar[65]=3./2.*armbbbar[180] + 1./2.*armbbbar[65] + 3*
   armbbbar[212];
   armbbbar[65]=MMZ*armbbbar[65];
   armbbbar[103]=55./4.*armbbbar[43] + armbbbar[215] + armbbbar[103] - 
   155./16. + armbbbar[44];
   armbbbar[103]=armbbbar[104] + armbbbar[101] + armbbbar[92] + 
   armbbbar[193] + 3*armbbbar[103] + 11./2.*armbbbar[47];
   armbbbar[103]=armbbbar[11]*armbbbar[103];
   armbbbar[180]=armbbbar[34] - 7*armbbbar[6];
   armbbbar[212]= - 27./4.*armbbbar[29];
   armbbbar[215]=99*armbbbar[33] + armbbbar[180] + armbbbar[212];
   armbbbar[215]= - 147./2.*armbbbar[17] - 147./4.*armbbbar[16] + 
   armbbbar[23] + 69./8.*armbbbar[24] - 3*armbbbar[32] - 27./4.*
   armbbbar[30] + 1./2.*armbbbar[215] - 13*armbbbar[35];
   armbbbar[224]=armbbbar[80] + 17./4. + armbbbar[79];
   armbbbar[224]=armbbbar[8]*armbbbar[224];
   armbbbar[103]=1./2.*armbbbar[103] + 3./4.*armbbbar[224] + 1./2.*
   armbbbar[215] + 9*armbbbar[25];
   armbbbar[215]= - armbbbar[11] - armbbbar[12];
   armbbbar[215]=armbbbar[12]*armbbbar[215];
   armbbbar[224]= - 1./4.*armbbbar[173] + armbbbar[215];
   armbbbar[224]=armbbbar[3]*armbbbar[224];
   armbbbar[226]= - 311./16. + armbbbar[44];
   armbbbar[226]=11./4.*armbbbar[43] + 1./4.*armbbbar[40] + 1./2.*
   armbbbar[226] - 11*armbbbar[45];
   armbbbar[226]=armbbbar[9] + 3*armbbbar[226] + 2*armbbbar[47];
   armbbbar[226]=armbbbar[12]*armbbbar[226];
   armbbbar[103]=9./8.*armbbbar[224] + 1./2.*armbbbar[103] + 
   armbbbar[226];
   armbbbar[103]=armbbbar[3]*armbbbar[103];
   armbbbar[224]=armbbbar[12]*armbbbar[22];
   armbbbar[224]=9*armbbbar[224] + 15./2.*armbbbar[75] + 7./4.*
   armbbbar[54] + 9*armbbbar[184];
   armbbbar[90]=armbbbar[117] + 115./16. + armbbbar[90];
   armbbbar[90]=armbbbar[3]*armbbbar[90];
   armbbbar[90]=1./16.*armbbbar[224] + armbbbar[90];
   armbbbar[90]=armbbbar[13]*armbbbar[90];
   armbbbar[60]=armbbbar[82] + armbbbar[69] + armbbbar[76] + 
   armbbbar[65] + armbbbar[96] + armbbbar[90] + 1./4.*armbbbar[60] + 
   armbbbar[103];
   armbbbar[60]=armbbbar[91]*armbbbar[60];
   armbbbar[65]=11./2.*armbbbar[41];
   armbbbar[69]=47./4.*armbbbar[45];
   armbbbar[76]=13*armbbbar[43];
   armbbbar[82]= - 17./6.*armbbbar[47];
   armbbbar[90]= - 5./6.*armbbbar[9];
   armbbbar[96]=armbbbar[11]*armbbbar[54];
   armbbbar[103]=armbbbar[15]*armbbbar[14];
   armbbbar[117]=armbbbar[132] + armbbbar[85] + 1./6.*armbbbar[96] + 
   armbbbar[90] + armbbbar[101] + armbbbar[100] + armbbbar[81] + 
   armbbbar[82] + armbbbar[76] + armbbbar[69] + armbbbar[65] + 
   armbbbar[140] + 3895./72. + armbbbar[103];
   armbbbar[224]= - armbbbar[8] - 125./12.*armbbbar[11];
   armbbbar[226]=armbbbar[224] - 211./12.*armbbbar[12];
   armbbbar[226]=armbbbar[3]*armbbbar[226];
   armbbbar[117]=1./2.*armbbbar[117] + armbbbar[226];
   armbbbar[226]= - 1./24.*armbbbar[54];
   armbbbar[229]=armbbbar[226] - 15*armbbbar[3];
   armbbbar[229]=armbbbar[13]*armbbbar[229];
   armbbbar[117]=1./2.*armbbbar[117] + armbbbar[229];
   armbbbar[117]=armbbbar[13]*armbbbar[117];
   armbbbar[229]= - armbbbar[15]*armbbbar[14];
   armbbbar[230]= - 3*armbbbar[40];
   armbbbar[231]= - 13*armbbbar[43];
   armbbbar[232]=17./6.*armbbbar[47];
   armbbbar[233]=5./6.*armbbbar[9];
   armbbbar[234]=armbbbar[11] + armbbbar[12];
   armbbbar[234]=9./2.*armbbbar[3]*armbbbar[234];
   armbbbar[235]= - 11*armbbbar[41];
   armbbbar[236]=armbbbar[234] + armbbbar[233] + armbbbar[194] + 
   armbbbar[97] + armbbbar[193] + armbbbar[232] + armbbbar[231] + 
   armbbbar[230] - 47./4.*armbbbar[45] + armbbbar[235] + armbbbar[136]
    + armbbbar[186] - 2971./72. + armbbbar[229];
   armbbbar[236]=armbbbar[3]*armbbbar[236];
   armbbbar[184]=3*armbbbar[184];
   armbbbar[237]=3*armbbbar[75];
   armbbbar[238]= - 9*armbbbar[22];
   armbbbar[239]=35./3.*armbbbar[26] + armbbbar[238];
   armbbbar[239]=armbbbar[12]*armbbbar[239];
   armbbbar[239]=armbbbar[239] + armbbbar[237] + armbbbar[184] + 29./3.
   *armbbbar[199];
   armbbbar[236]=1./32.*armbbbar[239] + armbbbar[236];
   armbbbar[236]=armbbbar[13]*armbbbar[236];
   armbbbar[239]=6409./72. + armbbbar[103];
   armbbbar[240]=7./4.*armbbbar[50];
   armbbbar[241]=77./6.*armbbbar[41];
   armbbbar[242]=13./3.*armbbbar[43];
   armbbbar[243]= - 157./18.*armbbbar[48];
   armbbbar[244]= - 5./18.*armbbbar[9];
   armbbbar[239]=armbbbar[64] + armbbbar[244] - armbbbar[42] + 
   armbbbar[243] - armbbbar[39] + 139./9.*armbbbar[47] + armbbbar[242]
    + 47./12.*armbbbar[45] + armbbbar[241] + 1./3.*armbbbar[239] + 
   armbbbar[240];
   armbbbar[245]=armbbbar[168] - 139./6.*armbbbar[47] - 463./24. + 
   armbbbar[200];
   armbbbar[245]=armbbbar[11]*armbbbar[245];
   armbbbar[246]= - 1./4.*armbbbar[23];
   armbbbar[152]=1./2.*armbbbar[245] - 25*armbbbar[25] + armbbbar[152]
    - 11*armbbbar[16] + armbbbar[246] + 1./4.*armbbbar[24] + 13./4.*
   armbbbar[35] + 11*armbbbar[32];
   armbbbar[245]= - 719./24. + armbbbar[200];
   armbbbar[247]= - 16./3.*armbbbar[47];
   armbbbar[203]=armbbbar[203] + 1./4.*armbbbar[245] + armbbbar[247];
   armbbbar[203]=armbbbar[12]*armbbbar[203];
   armbbbar[203]=1./2.*armbbbar[152] + armbbbar[203];
   armbbbar[203]=armbbbar[3]*armbbbar[203];
   armbbbar[75]=3./32.*armbbbar[75] + armbbbar[187] + 2./3.*
   armbbbar[109];
   armbbbar[75]=armbbbar[12]*armbbbar[75];
   armbbbar[187]= - 5./2.*armbbbar[48] + armbbbar[235] + 133./3. + 
   armbbbar[186];
   armbbbar[205]=armbbbar[204] + armbbbar[205];
   armbbbar[205]=armbbbar[3]*armbbbar[205];
   armbbbar[187]=1./2.*armbbbar[187] + 9*armbbbar[205];
   armbbbar[187]=MMt*armbbbar[3]*armbbbar[187];
   armbbbar[75]=armbbbar[187] + armbbbar[236] + armbbbar[203] + 1./16.*
   armbbbar[239] + armbbbar[75];
   armbbbar[75]=MMt*armbbbar[75];
   armbbbar[69]=armbbbar[90] + armbbbar[101] + armbbbar[92] + 
   armbbbar[81] + armbbbar[82] + armbbbar[76] + armbbbar[230] + 
   armbbbar[69] + armbbbar[136] + 2245./72. + armbbbar[103];
   armbbbar[69]=armbbbar[8]*armbbbar[69];
   armbbbar[203]=armbbbar[168] - 7./6.*armbbbar[47] + 53./24. + 
   armbbbar[135];
   armbbbar[203]=armbbbar[8]*armbbbar[203];
   armbbbar[205]= - 73./6.*armbbbar[35] - 35*armbbbar[31];
   armbbbar[236]=49./6.*armbbbar[17];
   armbbbar[205]=armbbbar[203] + 335./6.*armbbbar[25] + armbbbar[236]
    + 43./6.*armbbbar[16] + armbbbar[134] + armbbbar[88] - 5./6.*
   armbbbar[27] + 1./2.*armbbbar[205] - 19./3.*armbbbar[32];
   armbbbar[239]= - armbbbar[48] - 53./8. + armbbbar[47];
   armbbbar[239]=armbbbar[11]*armbbbar[239];
   armbbbar[245]=1./6.*armbbbar[42] + armbbbar[138] + 1./6.*
   armbbbar[39];
   armbbbar[245]=1./2.*MMH*armbbbar[245];
   armbbbar[248]=13./3. + armbbbar[221];
   armbbbar[248]=1./4.*armbbbar[12]*armbbbar[248];
   armbbbar[205]=armbbbar[141] + armbbbar[248] + armbbbar[245] + 1./2.*
   armbbbar[205] + 1./3.*armbbbar[239];
   armbbbar[239]=armbbbar[140] + 893./144. + armbbbar[139];
   armbbbar[249]= - 7./24.*armbbbar[47];
   armbbbar[239]=armbbbar[249] + armbbbar[80] + armbbbar[143] + 1./2.*
   armbbbar[239] + armbbbar[79];
   armbbbar[250]=3./16.*armbbbar[48];
   armbbbar[251]= - 7./2.*armbbbar[12] + 5./2.*armbbbar[8] - 
   armbbbar[11];
   armbbbar[251]=1./2.*armbbbar[3]*armbbbar[251];
   armbbbar[252]=25./6.*armbbbar[102];
   armbbbar[239]=armbbbar[252] + armbbbar[251] - armbbbar[42] + 
   armbbbar[250] + 1./2.*armbbbar[239] - armbbbar[39];
   armbbbar[239]=armbbbar[13]*armbbbar[239];
   armbbbar[253]= - 7./144.*armbbbar[47];
   armbbbar[254]=5./144.*armbbbar[48];
   armbbbar[255]=armbbbar[254] + armbbbar[253] + armbbbar[148] + 
   armbbbar[147] + 559./1152. + armbbbar[145];
   armbbbar[256]= - 1 + armbbbar[100];
   armbbbar[256]=armbbbar[12]*armbbbar[256];
   armbbbar[256]=armbbbar[150] + 3./2.*armbbbar[256];
   armbbbar[256]=armbbbar[3]*armbbbar[256];
   armbbbar[77]=armbbbar[100] + 7./6.*armbbbar[47] - 12*armbbbar[38] + 
   armbbbar[77] - 17./24. + armbbbar[135];
   armbbbar[77]=armbbbar[3]*armbbbar[77];
   armbbbar[77]=armbbbar[77] + 9*armbbbar[156];
   armbbbar[77]=armbbbar[13]*armbbbar[77];
   armbbbar[135]= - 2*armbbbar[38];
   armbbbar[157]=armbbbar[157] + armbbbar[135];
   armbbbar[157]=armbbbar[3]*armbbbar[157];
   armbbbar[157]=armbbbar[157] + 12*armbbbar[159];
   armbbbar[157]=3*MMt*armbbbar[157];
   armbbbar[159]=armbbbar[157] + armbbbar[77] + 1./2.*armbbbar[255] + 
   armbbbar[256];
   armbbbar[159]=MMt*armbbbar[159];
   armbbbar[159]=armbbbar[159] + 1./4.*armbbbar[205] + armbbbar[239];
   armbbbar[159]=MMt*armbbbar[159];
   armbbbar[205]=71./2. + armbbbar[139];
   armbbbar[205]=armbbbar[80] + 1./2.*armbbbar[205] + armbbbar[79];
   armbbbar[205]=armbbbar[42] + 1./2.*armbbbar[205] + armbbbar[39];
   armbbbar[205]=armbbbar[8]*armbbbar[205];
   armbbbar[239]=27./8.*armbbbar[37];
   armbbbar[255]=3./4.*armbbbar[44];
   armbbbar[257]=3./8.*armbbbar[40];
   armbbbar[258]=1./108.*armbbbar[42] - 103./432.*armbbbar[39] + 
   armbbbar[257] + armbbbar[255] - 13./3. + armbbbar[239];
   armbbbar[258]=MMH*armbbbar[258];
   armbbbar[259]=499./12. + 35*armbbbar[39];
   armbbbar[259]=1./8.*armbbbar[259] - armbbbar[42];
   armbbbar[259]=armbbbar[11]*armbbbar[259];
   armbbbar[205]=1./2.*armbbbar[258] + 1./9.*armbbbar[259] + 1./4.*
   armbbbar[205] + 7./16.*armbbbar[25] + 9./16.*armbbbar[17] + 1./4.*
   armbbbar[16] + 19./16.*armbbbar[23] + 11./3.*armbbbar[24] - 13./16.*
   armbbbar[31] - 7./8.*armbbbar[30] - 11./48.*armbbbar[29] + 5./24.*
   armbbbar[34] - 1./3.*armbbbar[36] - 15./16.*armbbbar[28];
   armbbbar[205]=1./2.*MMH*armbbbar[205];
   armbbbar[114]=1./2.*armbbbar[114];
   armbbbar[258]=armbbbar[144] + 3*armbbbar[162];
   armbbbar[258]=MMH*armbbbar[258];
   armbbbar[258]=armbbbar[258] + 7*armbbbar[120];
   armbbbar[258]=armbbbar[12]*armbbbar[258];
   armbbbar[258]=armbbbar[114] + armbbbar[119] + 1./4.*armbbbar[258];
   armbbbar[259]=armbbbar[52]*armbbbar[258];
   armbbbar[260]=armbbbar[16] + armbbbar[34] - armbbbar[29];
   armbbbar[261]= - armbbbar[8]*armbbbar[39];
   armbbbar[172]=1./18.*armbbbar[11]*armbbbar[172];
   armbbbar[261]=armbbbar[172] + 1./3.*armbbbar[260] + 1./8.*
   armbbbar[261];
   armbbbar[261]=MMH*armbbbar[261];
   armbbbar[262]= - 1./4.*armbbbar[142];
   armbbbar[263]= - 5./4.*armbbbar[8] + 7./3.*armbbbar[11];
   armbbbar[263]=armbbbar[11]*armbbbar[263];
   armbbbar[263]=armbbbar[262] + 1./3.*armbbbar[263];
   armbbbar[261]=1./2.*armbbbar[263] + armbbbar[261];
   armbbbar[261]=MMH*armbbbar[261];
   armbbbar[263]= - MMH*armbbbar[39];
   armbbbar[263]=armbbbar[153] + armbbbar[263];
   armbbbar[263]=MMH*armbbbar[263];
   armbbbar[264]= - armbbbar[12]*armbbbar[11];
   armbbbar[265]=1./4.*armbbbar[263] + 5*armbbbar[264];
   armbbbar[266]= - armbbbar[16] - armbbbar[35] + armbbbar[27];
   armbbbar[267]=5./4.*armbbbar[266] + armbbbar[154];
   armbbbar[267]=MMt*armbbbar[267];
   armbbbar[268]=MMH*armbbbar[39];
   armbbbar[268]=armbbbar[268] + armbbbar[8] - 11./12.*armbbbar[11];
   armbbbar[268]=armbbbar[13]*armbbbar[268];
   armbbbar[265]=1./3.*armbbbar[267] + 1./12.*armbbbar[265] + 
   armbbbar[268];
   armbbbar[265]=MMt*armbbbar[265];
   armbbbar[263]=armbbbar[13]*armbbbar[263];
   armbbbar[267]= - armbbbar[12]*MMH*armbbbar[11];
   armbbbar[261]=armbbbar[265] + 1./8.*armbbbar[263] + armbbbar[261] + 
   1./18.*armbbbar[267];
   armbbbar[263]=armbbbar[4]*armbbbar[261];
   armbbbar[99]= - armbbbar[8] + armbbbar[99];
   armbbbar[163]=armbbbar[42] + 1./2.*armbbbar[163] + armbbbar[39];
   armbbbar[163]=MMH*armbbbar[163];
   armbbbar[265]= - 59./6.*armbbbar[12];
   armbbbar[99]=15./8.*armbbbar[13] + armbbbar[265] + 43./6.*
   armbbbar[99] + armbbbar[163];
   armbbbar[99]=armbbbar[13]*armbbbar[99];
   armbbbar[267]=19*armbbbar[8];
   armbbbar[268]=armbbbar[267] + 1403./9.*armbbbar[11];
   armbbbar[269]=161./54. + armbbbar[194];
   armbbbar[269]=MMH*armbbbar[269];
   armbbbar[268]=151./6.*armbbbar[12] + 1./4.*armbbbar[268] + 
   armbbbar[269];
   armbbbar[268]=armbbbar[12]*armbbbar[268];
   armbbbar[270]=3./64.*armbbbar[142];
   armbbbar[271]=armbbbar[8] - 7./12.*armbbbar[11];
   armbbbar[271]=armbbbar[11]*armbbbar[271];
   armbbbar[99]=1./4.*armbbbar[263] + 1./8.*armbbbar[259] + 
   armbbbar[159] + 1./8.*armbbbar[99] + 1./8.*armbbbar[268] + 
   armbbbar[205] + armbbbar[270] + 1./3.*armbbbar[271];
   armbbbar[99]=armbbbar[4]*armbbbar[99];
   armbbbar[159]=2551./108. + armbbbar[179];
   armbbbar[159]= - 11*armbbbar[43] + armbbbar[166] + 11*armbbbar[45]
    + 1./2.*armbbbar[159] + armbbbar[136];
   armbbbar[159]= - 11./3.*armbbbar[39] + 1./2.*armbbbar[159] - 1./3.*
   armbbbar[47];
   armbbbar[259]=1./3.*armbbbar[48];
   armbbbar[263]=5./6.*armbbbar[42];
   armbbbar[268]= - 1./6.*armbbbar[9];
   armbbbar[271]=11./8.*armbbbar[182];
   armbbbar[272]=11./8.*armbbbar[183];
   armbbbar[159]=armbbbar[272] + armbbbar[271] + armbbbar[268] + 
   armbbbar[263] + 1./2.*armbbbar[159] + armbbbar[259];
   armbbbar[159]=armbbbar[11]*armbbbar[159];
   armbbbar[273]=armbbbar[186] - 57761./216. + armbbbar[179];
   armbbbar[273]=armbbbar[166] + armbbbar[78] + 1./2.*armbbbar[273] + 
   armbbbar[136];
   armbbbar[273]=armbbbar[192] + armbbbar[188] + armbbbar[191] + 
   armbbbar[190] + armbbbar[189] + 1./8.*armbbbar[273] - 4./3.*
   armbbbar[47];
   armbbbar[273]=armbbbar[12]*armbbbar[273];
   armbbbar[93]=armbbbar[115] + 1./4.*armbbbar[93] + armbbbar[116];
   armbbbar[93]=armbbbar[13]*armbbbar[93];
   armbbbar[115]=armbbbar[8] + 3./2.*armbbbar[162];
   armbbbar[115]=MMH*armbbbar[115];
   armbbbar[115]=1./2.*armbbbar[115] + armbbbar[120];
   armbbbar[115]=armbbbar[12]*armbbbar[115];
   armbbbar[89]=armbbbar[93] + armbbbar[89] + armbbbar[115];
   armbbbar[93]=armbbbar[52]*armbbbar[89];
   armbbbar[115]=1./2. + armbbbar[42];
   armbbbar[115]=MMH*armbbbar[115];
   armbbbar[115]=armbbbar[201] + 1./2.*armbbbar[115];
   armbbbar[115]=MMH*armbbbar[115];
   armbbbar[201]=25./24.*armbbbar[12] + 1./8.*armbbbar[209] + 13./64.*
   armbbbar[8] + armbbbar[11];
   armbbbar[201]=armbbbar[12]*armbbbar[201];
   armbbbar[274]= - 145./2.*armbbbar[12] + armbbbar[8] - 7*armbbbar[11]
   ;
   armbbbar[274]=1./2.*armbbbar[274] + 17*armbbbar[13];
   armbbbar[274]=armbbbar[13]*armbbbar[274];
   armbbbar[93]=1./8.*armbbbar[93] + 1./16.*armbbbar[274] + 1./8.*
   armbbbar[115] + armbbbar[201];
   armbbbar[93]=armbbbar[52]*armbbbar[93];
   armbbbar[201]=2*armbbbar[18];
   armbbbar[274]=17./16.*armbbbar[34];
   armbbbar[275]=23./16.*armbbbar[29];
   armbbbar[276]=armbbbar[275] - 1./16.*armbbbar[6] + armbbbar[201] + 
   armbbbar[274];
   armbbbar[185]=3./2.*armbbbar[185] + 59./54.*armbbbar[42] + 211./54.*
   armbbbar[39] + armbbbar[230] + 37./24. + armbbbar[136];
   armbbbar[185]=1./8.*MMH*armbbbar[185];
   armbbbar[277]= - 7./3.*armbbbar[8] + 3./2.*armbbbar[11];
   armbbbar[277]=armbbbar[11]*armbbbar[277];
   armbbbar[153]=armbbbar[153] + 79./2.*armbbbar[12];
   armbbbar[153]=armbbbar[12]*armbbbar[153];
   armbbbar[153]=armbbbar[277] + armbbbar[153];
   armbbbar[153]=armbbbar[3]*armbbbar[153];
   armbbbar[278]=11./16.*armbbbar[30];
   armbbbar[279]= - 1./4.*armbbbar[31];
   armbbbar[280]=pow(armbbbar[14],2);
   armbbbar[281]= - armbbbar[15]*armbbbar[280];
   armbbbar[282]= - 319./192.*armbbbar[24];
   armbbbar[283]=25./64.*armbbbar[23];
   armbbbar[69]=armbbbar[99] + armbbbar[93] + armbbbar[75] + 1./2.*
   armbbbar[117] + 1./8.*armbbbar[153] + armbbbar[273] + armbbbar[185]
    + 1./2.*armbbbar[159] + 1./8.*armbbbar[69] - 435./32.*armbbbar[17]
    - 11./2.*armbbbar[16] + armbbbar[283] + armbbbar[282] - 179./192.*
   armbbbar[27] - 23./64.*armbbbar[32] + 1./3.*armbbbar[281] + 
   armbbbar[279] + armbbbar[278] + 169./96.*armbbbar[35] + 1./3.*
   armbbbar[276] + 11./2.*armbbbar[33];
   armbbbar[69]=armbbbar[4]*armbbbar[69];
   armbbbar[75]=armbbbar[121] + 3./2.*armbbbar[17] + armbbbar[151] + 
   armbbbar[118] + armbbbar[24];
   armbbbar[75]=armbbbar[56]*armbbbar[75];
   armbbbar[93]= - 21*armbbbar[41];
   armbbbar[99]= - 2225./12. + armbbbar[93];
   armbbbar[117]=armbbbar[56] + armbbbar[3];
   armbbbar[117]=armbbbar[13]*armbbbar[117];
   armbbbar[118]= - MMZ*armbbbar[56];
   armbbbar[75]=27*armbbbar[118] + 3./2.*armbbbar[117] + 3./2.*
   armbbbar[131] + 9./2.*armbbbar[130] + 3*armbbbar[87] + 3*
   armbbbar[75] + armbbbar[101] + armbbbar[176] + 1./8.*armbbbar[99] + 
   armbbbar[95];
   armbbbar[75]=MMZ*armbbbar[75];
   armbbbar[99]=armbbbar[221] + 69./2.*armbbbar[47] + 47./4.*
   armbbbar[43] - 119*armbbbar[45] - 9*armbbbar[38] + 11./4.*
   armbbbar[41] - 9./2.*armbbbar[50] - 9659./144. + armbbbar[229];
   armbbbar[117]=3./2.*armbbbar[130];
   armbbbar[99]=9./2.*armbbbar[146] + 27./4.*armbbbar[112] + 
   armbbbar[117] + armbbbar[128] + 1./2.*armbbbar[99] + 5./3.*
   armbbbar[9];
   armbbbar[99]=armbbbar[13]*armbbbar[99];
   armbbbar[118]=9./2.*armbbbar[50];
   armbbbar[130]=9*armbbbar[48] - 69./2.*armbbbar[47] + 283./4.*
   armbbbar[43] + 119*armbbbar[45] + 15*armbbbar[38] - 53./4.*
   armbbbar[41] + armbbbar[118] + 13469./144. + armbbbar[103];
   armbbbar[94]=armbbbar[117] + 9./4.*armbbbar[94] - 5./3.*armbbbar[9]
    + 1./2.*armbbbar[130] + armbbbar[101];
   armbbbar[94]=armbbbar[12]*armbbbar[94];
   armbbbar[117]= - 5 + armbbbar[78];
   armbbbar[117]=1./4.*armbbbar[117] + armbbbar[42];
   armbbbar[117]=MMH*armbbbar[117];
   armbbbar[130]= - armbbbar[3]*armbbbar[124];
   armbbbar[94]=armbbbar[99] + 9./4.*armbbbar[130] + armbbbar[94] + 
   armbbbar[117] + 37./4.*armbbbar[11] + 3./4.*armbbbar[8] - 1373./24.*
   armbbbar[25] + 559./8.*armbbbar[17] + 37./4.*armbbbar[16] - 3./8.*
   armbbbar[23] + 3./4.*armbbbar[24] + 259./12.*armbbbar[27] + 361./24.
   *armbbbar[32] + 9./8.*armbbbar[31] - 3./2.*armbbbar[30] - 367./8.*
   armbbbar[33] + 11./3.*armbbbar[35];
   armbbbar[75]=1./2.*armbbbar[94] + armbbbar[75];
   armbbbar[75]=MMZ*armbbbar[75];
   armbbbar[93]=17./2. + armbbbar[93];
   armbbbar[93]=armbbbar[101] + armbbbar[176] + 1./8.*armbbbar[93] + 
   armbbbar[95];
   armbbbar[93]=armbbbar[12]*armbbbar[93];
   armbbbar[94]= - 1./2. + 21*armbbbar[41];
   armbbbar[78]=armbbbar[194] - 33*armbbbar[43] + 1./8.*armbbbar[94] + 
   armbbbar[78];
   armbbbar[78]=armbbbar[13]*armbbbar[78];
   armbbbar[94]= - 239./4.*armbbbar[25] + 247./4.*armbbbar[17] + 135./4.
   *armbbbar[27] + 27./2.*armbbbar[32] - 189./4.*armbbbar[33] - 
   armbbbar[35];
   armbbbar[78]= - MMZ + armbbbar[78] + 1./2.*armbbbar[94] + 
   armbbbar[93];
   armbbbar[78]=MMZ*armbbbar[78];
   armbbbar[93]=MMH*armbbbar[105];
   armbbbar[94]=armbbbar[126] + 419./48.*armbbbar[12];
   armbbbar[94]=armbbbar[12]*armbbbar[94];
   armbbbar[95]=49./6.*armbbbar[13] + armbbbar[123] - 859./48.*
   armbbbar[12];
   armbbbar[95]=armbbbar[13]*armbbbar[95];
   armbbbar[93]=armbbbar[95] + 1./2.*armbbbar[93] + armbbbar[94];
   armbbbar[78]=1./2.*armbbbar[93] + armbbbar[78];
   armbbbar[78]=MMZ*armbbbar[78];
   armbbbar[93]=1./2.*armbbbar[124] + armbbbar[170];
   armbbbar[93]=armbbbar[52]*armbbbar[127]*armbbbar[93];
   armbbbar[78]=armbbbar[78] + 9./16.*armbbbar[93];
   armbbbar[78]=armbbbar[52]*armbbbar[78];
   armbbbar[93]=17*armbbbar[12];
   armbbbar[94]= - 57./8.*armbbbar[11] + armbbbar[93];
   armbbbar[94]=armbbbar[12]*armbbbar[94];
   armbbbar[95]=235./3.*armbbbar[13] + 49*armbbbar[11] - 173*
   armbbbar[12];
   armbbbar[95]=armbbbar[13]*armbbbar[95];
   armbbbar[94]=armbbbar[94] + 1./8.*armbbbar[95];
   armbbbar[75]=armbbbar[78] + 1./2.*armbbbar[94] + armbbbar[75];
   armbbbar[75]=armbbbar[52]*armbbbar[75];
   armbbbar[78]= - 21./4. + armbbbar[61];
   armbbbar[78]=armbbbar[56]*armbbbar[78];
   armbbbar[94]=armbbbar[12]*armbbbar[73];
   armbbbar[95]=armbbbar[13]*armbbbar[73];
   armbbbar[99]= - MMZ*armbbbar[73];
   armbbbar[78]=9*armbbbar[99] + armbbbar[95] - 1./2.*armbbbar[3] + 
   armbbbar[94] + armbbbar[78] + armbbbar[86];
   armbbbar[78]=MMZ*armbbbar[78];
   armbbbar[86]=armbbbar[118] + 8423./144. + armbbbar[103];
   armbbbar[95]=1./4.*armbbbar[56] + armbbbar[74];
   armbbbar[94]=3./2.*armbbbar[3] + 1./2.*armbbbar[95] + armbbbar[94];
   armbbbar[94]=armbbbar[13]*armbbbar[94];
   armbbbar[95]=1./4.*armbbbar[25] + 9./4.*armbbbar[17] + armbbbar[246]
    - armbbbar[30] + 5./4.*armbbbar[24];
   armbbbar[95]=armbbbar[56]*armbbbar[95];
   armbbbar[99]=9./4.*armbbbar[56] + armbbbar[74];
   armbbbar[99]=armbbbar[12]*armbbbar[99];
   armbbbar[78]=3./2.*armbbbar[78] + 3./4.*armbbbar[94] + 3./4.*
   armbbbar[131] + 3./8.*armbbbar[99] + 15./32.*armbbbar[87] + 3./8.*
   armbbbar[95] - 5./12.*armbbbar[9] - 3./8.*armbbbar[42] + 9./8.*
   armbbbar[48] - 69./16.*armbbbar[47] + 59./16.*armbbbar[43] + 119./8.
   *armbbbar[45] + 3./2.*armbbbar[38] + 1./8.*armbbbar[86] - 
   armbbbar[41];
   armbbbar[78]=MMZ*armbbbar[78];
   armbbbar[86]=53*armbbbar[45];
   armbbbar[87]=19./4.*armbbbar[43];
   armbbbar[94]=armbbbar[87] + armbbbar[86] - 19./4.*armbbbar[41] + 
   2897./144. + armbbbar[103];
   armbbbar[94]=armbbbar[268] + 1./8.*armbbbar[94] + armbbbar[247];
   armbbbar[94]=armbbbar[12]*armbbbar[94];
   armbbbar[95]= - 199./6. + armbbbar[41];
   armbbbar[99]=7./3.*armbbbar[47];
   armbbbar[95]=1./16.*armbbbar[95] + armbbbar[99];
   armbbbar[95]=armbbbar[13]*armbbbar[95];
   armbbbar[105]= - 1./8.*armbbbar[33] + armbbbar[35];
   armbbbar[105]= - 35./24.*armbbbar[11] - 93./8.*armbbbar[25] - 63./8.
   *armbbbar[17] + 5./4.*armbbbar[16] - 3./8.*armbbbar[27] + 11*
   armbbbar[105] + 1./2.*armbbbar[32];
   armbbbar[117]=armbbbar[12]*armbbbar[11];
   armbbbar[118]=armbbbar[3]*armbbbar[117];
   armbbbar[75]=1./2.*armbbbar[75] + armbbbar[78] + armbbbar[95] + 9./
   32.*armbbbar[118] + 1./4.*armbbbar[105] + armbbbar[94];
   armbbbar[75]=armbbbar[52]*armbbbar[75];
   armbbbar[78]= - 29./6.*armbbbar[47];
   armbbbar[94]=9*armbbbar[39];
   armbbbar[95]=7./6.*armbbbar[9];
   armbbbar[105]=3*armbbbar[40];
   armbbbar[123]=armbbbar[95] + armbbbar[101] + armbbbar[92] + 
   armbbbar[94] + armbbbar[78] + armbbbar[76] + armbbbar[105] + 179./4.
   *armbbbar[45] + armbbbar[79] + 4213./72. + armbbbar[103];
   armbbbar[123]=armbbbar[11]*armbbbar[123];
   armbbbar[86]=armbbbar[87] + armbbbar[80] + armbbbar[86] + 14171./144.
    + armbbbar[103];
   armbbbar[126]= - 2./3.*armbbbar[9];
   armbbbar[86]=armbbbar[126] + 1./2.*armbbbar[86] - 8./3.*armbbbar[47]
   ;
   armbbbar[86]=armbbbar[12]*armbbbar[86];
   armbbbar[130]=1./3.*armbbbar[180] + armbbbar[212];
   armbbbar[131]=armbbbar[130] - 33*armbbbar[33];
   armbbbar[151]=77./24.*armbbbar[24];
   armbbbar[131]=33*armbbbar[17] + 93./4.*armbbbar[16] + armbbbar[151]
    + 1./2.*armbbbar[131] + armbbbar[32];
   armbbbar[153]=armbbbar[175] + armbbbar[264];
   armbbbar[153]=9./8.*armbbbar[3]*armbbbar[153];
   armbbbar[159]=17./2. + 9*armbbbar[40];
   armbbbar[159]=1./8.*armbbbar[8]*armbbbar[159];
   armbbbar[86]=armbbbar[153] + armbbbar[86] + 1./4.*armbbbar[123] + 
   armbbbar[159] + 1./2.*armbbbar[131] - armbbbar[25];
   armbbbar[86]=armbbbar[3]*armbbbar[86];
   armbbbar[123]=9./4. + armbbbar[70];
   armbbbar[123]=armbbbar[22]*armbbbar[56]*armbbbar[123];
   armbbbar[123]=armbbbar[228] + armbbbar[222] + armbbbar[123] + 
   armbbbar[227];
   armbbbar[123]=MMt*armbbbar[123];
   armbbbar[131]= - 13./4. + armbbbar[61];
   armbbbar[131]=armbbbar[22]*armbbbar[56]*armbbbar[131];
   armbbbar[131]=armbbbar[223] + armbbbar[219] + armbbbar[131] + 
   armbbbar[165];
   armbbbar[170]=MMt*armbbbar[161];
   armbbbar[131]=3./2.*armbbbar[131] + armbbbar[170];
   armbbbar[170]=MMZ*armbbbar[220];
   armbbbar[131]=1./2.*armbbbar[131] + 6*armbbbar[170];
   armbbbar[131]=MMZ*armbbbar[131];
   armbbbar[170]=7./4.*armbbbar[213] + 3*armbbbar[227];
   armbbbar[175]=armbbbar[12]*armbbbar[170];
   armbbbar[61]= - 3 + 7./4.*armbbbar[61];
   armbbbar[61]=armbbbar[22]*armbbbar[61];
   armbbbar[61]=armbbbar[175] + armbbbar[61] + 7./4.*armbbbar[214];
   armbbbar[175]=armbbbar[181] - 3./4. + 2*armbbbar[44];
   armbbbar[175]=armbbbar[155]*armbbbar[175];
   armbbbar[170]=1./2.*armbbbar[170] + 3*armbbbar[219];
   armbbbar[170]=armbbbar[13]*armbbbar[170];
   armbbbar[61]=9*armbbbar[131] + 1./2.*armbbbar[123] + 3./2.*
   armbbbar[170] + 3./4.*armbbbar[61] + armbbbar[175];
   armbbbar[61]=MMZ*armbbbar[61];
   armbbbar[123]= - 33*armbbbar[57];
   armbbbar[131]=1./2.*armbbbar[22] + armbbbar[218];
   armbbbar[131]=armbbbar[12]*armbbbar[131];
   armbbbar[62]=9./2.*armbbbar[131] + armbbbar[237] + 9./4.*
   armbbbar[63] + 9./4.*armbbbar[62] + 9*armbbbar[56] + armbbbar[123]
    + 1./16.*armbbbar[54];
   armbbbar[131]=armbbbar[166] - 1./2. + armbbbar[79];
   armbbbar[131]=armbbbar[12]*armbbbar[131];
   armbbbar[170]=1 + armbbbar[230];
   armbbbar[175]=armbbbar[11]*armbbbar[170];
   armbbbar[131]=1./2.*armbbbar[175] + armbbbar[131];
   armbbbar[131]=armbbbar[3]*armbbbar[131];
   armbbbar[176]= - 25./3.*armbbbar[47] + 151./4.*armbbbar[43] + 
   armbbbar[40] + 205./2.*armbbbar[45] - armbbbar[44] - 907./24. + 
   armbbbar[103];
   armbbbar[126]=3./2.*armbbbar[131] + armbbbar[126] - 2*armbbbar[42]
    - 2*armbbbar[48] + 1./2.*armbbbar[176] + armbbbar[193];
   armbbbar[126]=armbbbar[3]*armbbbar[126];
   armbbbar[131]=1./4.*armbbbar[211] + armbbbar[165];
   armbbbar[176]=armbbbar[12]*armbbbar[131];
   armbbbar[131]=1./2.*armbbbar[131] + armbbbar[222];
   armbbbar[131]=armbbbar[13]*armbbbar[131];
   armbbbar[70]=3 + 1./2.*armbbbar[70];
   armbbbar[70]=armbbbar[22]*armbbbar[70];
   armbbbar[70]=3./4.*armbbbar[131] + 3./8.*armbbbar[176] + 3./32.*
   armbbbar[218] + 7*armbbbar[26] + 3./16.*armbbbar[70];
   armbbbar[70]=MMt*armbbbar[70];
   armbbbar[131]=armbbbar[12]*armbbbar[213];
   armbbbar[131]=armbbbar[131] + armbbbar[22] + armbbbar[214];
   armbbbar[131]=armbbbar[13]*armbbbar[131];
   armbbbar[61]=3*armbbbar[61] + armbbbar[70] + 9./32.*armbbbar[131] + 
   1./8.*armbbbar[62] + armbbbar[126];
   armbbbar[61]=MMZ*armbbbar[61];
   armbbbar[62]=11./2.*armbbbar[59];
   armbbbar[70]=1./2.*armbbbar[229] + armbbbar[62] + 10783./288. + 
   armbbbar[58];
   armbbbar[126]=7./24.*armbbbar[41];
   armbbbar[131]=armbbbar[54]*armbbbar[32];
   armbbbar[70]=497./48.*armbbbar[45] + 1./48.*armbbbar[131] + 
   armbbbar[126] + 1./2.*armbbbar[70] - armbbbar[44];
   armbbbar[176]=armbbbar[27] - armbbbar[16];
   armbbbar[181]= - armbbbar[25] + armbbbar[176] - armbbbar[17];
   armbbbar[181]=armbbbar[26]*armbbbar[181];
   armbbbar[213]= - armbbbar[12]*armbbbar[26];
   armbbbar[214]=armbbbar[213] + armbbbar[181] + armbbbar[109];
   armbbbar[219]= - 3*armbbbar[46];
   armbbbar[221]=1./2.*armbbbar[48];
   armbbbar[222]=armbbbar[221] - 77./3.*armbbbar[47] - armbbbar[41] - 
   145./24. + armbbbar[219];
   armbbbar[222]=armbbbar[3]*armbbbar[222];
   armbbbar[214]=11./6.*armbbbar[214] + armbbbar[222];
   armbbbar[222]= - 1 + armbbbar[105];
   armbbbar[222]=armbbbar[155]*armbbbar[222];
   armbbbar[223]=3*armbbbar[222];
   armbbbar[227]= - 11./12.*armbbbar[26] + armbbbar[223];
   armbbbar[227]=armbbbar[13]*armbbbar[227];
   armbbbar[228]=6*armbbbar[216];
   armbbbar[214]=armbbbar[228] + 1./2.*armbbbar[214] + armbbbar[227];
   armbbbar[214]=MMt*armbbbar[214];
   armbbbar[227]=11*armbbbar[57];
   armbbbar[226]=armbbbar[227] + armbbbar[226];
   armbbbar[246]=armbbbar[16]*armbbbar[226];
   armbbbar[226]=armbbbar[11]*armbbbar[226];
   armbbbar[247]=3*armbbbar[213] - 1./12.*armbbbar[54] + 3*
   armbbbar[199];
   armbbbar[273]= - 203./24. - armbbbar[41];
   armbbbar[273]=armbbbar[3]*armbbbar[273];
   armbbbar[247]=1./8.*armbbbar[247] + armbbbar[273];
   armbbbar[247]=armbbbar[13]*armbbbar[247];
   armbbbar[273]=11./16.*armbbbar[68];
   armbbbar[276]= - 247./144.*armbbbar[39];
   armbbbar[284]= - armbbbar[25]*armbbbar[54];
   armbbbar[285]=31./144.*armbbbar[48];
   armbbbar[286]=55./144.*armbbbar[42];
   armbbbar[287]=armbbbar[8]*armbbbar[57];
   armbbbar[288]=11./16.*armbbbar[287];
   armbbbar[60]=armbbbar[60] + armbbbar[69] + armbbbar[75] + 
   armbbbar[61] + armbbbar[214] + armbbbar[247] + armbbbar[86] + 1./8.*
   armbbbar[226] + armbbbar[288] - 79./288.*armbbbar[9] + armbbbar[286]
    + armbbbar[285] + 1./96.*armbbbar[284] + armbbbar[276] - 679./288.*
   armbbbar[47] - 25./24.*armbbbar[43] + 1./8.*armbbbar[246] + 
   armbbbar[273] + 1./4.*armbbbar[70] + armbbbar[40];
   armbbbar[60]=armbbbar[91]*armbbbar[60];
   armbbbar[61]=4./3.*armbbbar[41];
   armbbbar[69]= - 3./4. + armbbbar[61];
   armbbbar[70]=pow(CW,2);
   armbbbar[69]=armbbbar[70]*armbbbar[69];
   armbbbar[75]=337./128. + armbbbar[41];
   armbbbar[86]=armbbbar[32] + armbbbar[27];
   armbbbar[86]=1./2.*armbbbar[86] - armbbbar[16];
   armbbbar[86]=armbbbar[55]*armbbbar[86];
   armbbbar[214]= - armbbbar[17]*armbbbar[55];
   armbbbar[226]= - armbbbar[25]*armbbbar[55];
   armbbbar[246]= - armbbbar[11]*armbbbar[55];
   armbbbar[247]= - armbbbar[55] - 9*armbbbar[56];
   armbbbar[247]=armbbbar[12]*armbbbar[247];
   armbbbar[289]= - armbbbar[3] - armbbbar[55] - armbbbar[56];
   armbbbar[289]=armbbbar[13]*armbbbar[289];
   armbbbar[290]=1 + 1./2.*armbbbar[70];
   armbbbar[291]=armbbbar[56]*armbbbar[290];
   armbbbar[291]=1./8.*armbbbar[55] + 3*armbbbar[291];
   armbbbar[291]=MMZ*armbbbar[291];
   armbbbar[292]= - 3*armbbbar[70];
   armbbbar[293]= - 41./4. + armbbbar[292];
   armbbbar[293]=armbbbar[43]*armbbbar[293];
   armbbbar[294]=3./4.*armbbbar[42];
   armbbbar[69]=3*armbbbar[291] + 3./8.*armbbbar[289] + 3./8.*
   armbbbar[112] + 1./8.*armbbbar[247] + 1./4.*armbbbar[246] + 
   armbbbar[128] + 3./4.*armbbbar[107] + armbbbar[294] + 3./8.*
   armbbbar[226] + 1./8.*armbbbar[214] + 1./4.*armbbbar[86] + 
   armbbbar[293] + armbbbar[69] + 5./3.*armbbbar[75] + armbbbar[125];
   armbbbar[69]=MMZ*armbbbar[69];
   armbbbar[75]=3*armbbbar[70];
   armbbbar[86]=53./4. + armbbbar[75];
   armbbbar[86]=armbbbar[45]*armbbbar[86];
   armbbbar[107]= - 3 - 4./3.*armbbbar[70];
   armbbbar[107]=armbbbar[47]*armbbbar[107];
   armbbbar[128]=armbbbar[11]*armbbbar[55];
   armbbbar[247]= - armbbbar[12]*armbbbar[55];
   armbbbar[289]= - armbbbar[13]*armbbbar[55];
   armbbbar[86]=1./2.*armbbbar[289] + 1./2.*armbbbar[247] + 3./4.*
   armbbbar[128] + armbbbar[268] + armbbbar[107] + 1./2.*armbbbar[86]
    + 1./6.*armbbbar[70] - 2./3.*armbbbar[41] + 25./3. + 1./8.*
   armbbbar[103];
   armbbbar[86]=armbbbar[13]*armbbbar[86];
   armbbbar[107]=1./8.*armbbbar[229];
   armbbbar[268]= - 53./4. + armbbbar[292];
   armbbbar[268]=1./2.*armbbbar[45]*armbbbar[268];
   armbbbar[291]=3 + 4./3.*armbbbar[70];
   armbbbar[291]=armbbbar[47]*armbbbar[291];
   armbbbar[292]=1./4.*armbbbar[128] + armbbbar[198] + armbbbar[291] - 
   3*armbbbar[43] + armbbbar[268] - 1./6.*armbbbar[70] + 2*armbbbar[41]
    - 22./3. + armbbbar[107];
   armbbbar[292]=armbbbar[12]*armbbbar[292];
   armbbbar[295]= - 55./96.*armbbbar[25] + 55./96.*armbbbar[17] + 119./
   96.*armbbbar[27] - 29./32.*armbbbar[33] - 1./3.*armbbbar[32];
   armbbbar[296]=17./3.*armbbbar[41];
   armbbbar[297]=armbbbar[296] + 3*armbbbar[43];
   armbbbar[297]=armbbbar[12]*armbbbar[297];
   armbbbar[298]=armbbbar[296] - armbbbar[43];
   armbbbar[298]=armbbbar[13]*armbbbar[298];
   armbbbar[295]=1./32.*armbbbar[298] + 1./3.*armbbbar[295] + 1./32.*
   armbbbar[297];
   armbbbar[297]=pow(armbbbar[2],2);
   armbbbar[295]=armbbbar[297]*armbbbar[295];
   armbbbar[298]=71./2.*armbbbar[33] - 11*armbbbar[35];
   armbbbar[298]= - 5./2.*armbbbar[11] + 349./12.*armbbbar[25] - 277./
   12.*armbbbar[17] - 5./2.*armbbbar[16] - 83./12.*armbbbar[27] + 1./2.
   *armbbbar[298] - 25./3.*armbbbar[32];
   armbbbar[69]=armbbbar[69] + armbbbar[295] + armbbbar[86] + 1./6.*
   armbbbar[298] + armbbbar[292];
   armbbbar[69]=MMZ*armbbbar[69];
   armbbbar[86]= - 1./4. + armbbbar[61];
   armbbbar[86]=armbbbar[70]*armbbbar[86];
   armbbbar[292]= - 5./64. + 1./3.*armbbbar[41];
   armbbbar[86]=armbbbar[294] + armbbbar[293] + armbbbar[86] + 5*
   armbbbar[292] + armbbbar[125];
   armbbbar[86]=armbbbar[12]*armbbbar[86];
   armbbbar[125]= - 1./4. - 4./3.*armbbbar[41];
   armbbbar[125]=armbbbar[70]*armbbbar[125];
   armbbbar[75]=41./4. + armbbbar[75];
   armbbbar[75]=armbbbar[43]*armbbbar[75];
   armbbbar[75]= - 3./4.*armbbbar[42] + armbbbar[75] + armbbbar[125] + 
   armbbbar[113] - 7./64. - 5./3.*armbbbar[41];
   armbbbar[75]=armbbbar[13]*armbbbar[75];
   armbbbar[113]=1./4.*armbbbar[35];
   armbbbar[125]= - 1./3.*armbbbar[27] - 2./3.*armbbbar[32] + 
   armbbbar[33] + armbbbar[113];
   armbbbar[125]=armbbbar[70]*armbbbar[125];
   armbbbar[292]=1 + armbbbar[70];
   armbbbar[292]=armbbbar[70]*armbbbar[292];
   armbbbar[292]=1 + armbbbar[292];
   armbbbar[292]=MMZ*armbbbar[292];
   armbbbar[293]= - 53 - 23./2.*armbbbar[70];
   armbbbar[293]=armbbbar[17]*armbbbar[293];
   armbbbar[294]=25 + 17./4.*armbbbar[70];
   armbbbar[294]=armbbbar[25]*armbbbar[294];
   armbbbar[75]=1./2.*armbbbar[292] + armbbbar[75] + armbbbar[86] + 1./
   3.*armbbbar[294] + 1./6.*armbbbar[293] + armbbbar[125] - 49./12.*
   armbbbar[27] - 9./4.*armbbbar[32] + 19./3.*armbbbar[33] + 
   armbbbar[113];
   armbbbar[75]=MMZ*armbbbar[75];
   armbbbar[86]= - 2*armbbbar[13];
   armbbbar[113]=armbbbar[225] + armbbbar[86];
   armbbbar[113]=armbbbar[13]*armbbbar[113];
   armbbbar[75]=armbbbar[75] - 2*armbbbar[124] + armbbbar[113];
   armbbbar[75]=MMZ*armbbbar[75];
   armbbbar[75]=armbbbar[75] + 9./64.*armbbbar[122];
   armbbbar[75]=armbbbar[52]*armbbbar[75];
   armbbbar[98]=armbbbar[137] + armbbbar[98];
   armbbbar[98]=armbbbar[12]*armbbbar[98];
   armbbbar[113]=23./12.*armbbbar[11] + 5*armbbbar[12];
   armbbbar[113]=armbbbar[12]*armbbbar[113];
   armbbbar[122]= - 19*armbbbar[13] + 97*armbbbar[11] - 161*
   armbbbar[12];
   armbbbar[122]=armbbbar[13]*armbbbar[122];
   armbbbar[113]=armbbbar[113] + 1./12.*armbbbar[122];
   armbbbar[113]=armbbbar[297]*armbbbar[113];
   armbbbar[93]= - armbbbar[11] + armbbbar[93];
   armbbbar[93]=1./2.*armbbbar[93] - 4*armbbbar[13];
   armbbbar[93]=armbbbar[13]*armbbbar[93];
   armbbbar[69]=armbbbar[75] + armbbbar[69] + 1./24.*armbbbar[113] + 1./
   2.*armbbbar[98] + 1./3.*armbbbar[93];
   armbbbar[69]=armbbbar[52]*armbbbar[69];
   armbbbar[75]= - 3*armbbbar[73];
   armbbbar[93]=pow(armbbbar[55],2);
   armbbbar[98]= - armbbbar[93] + armbbbar[75];
   armbbbar[98]=armbbbar[12]*armbbbar[98];
   armbbbar[75]= - 5*armbbbar[93] + armbbbar[75];
   armbbbar[75]=armbbbar[13]*armbbbar[75];
   armbbbar[113]= - 3./2.*armbbbar[16] + armbbbar[32] + 1./2.*
   armbbbar[27];
   armbbbar[113]=armbbbar[55]*armbbbar[113];
   armbbbar[113]=13./8. + armbbbar[113];
   armbbbar[113]=armbbbar[55]*armbbbar[113];
   armbbbar[122]= - armbbbar[17]*armbbbar[93];
   armbbbar[124]= - armbbbar[25]*armbbbar[93];
   armbbbar[125]= - armbbbar[11]*armbbbar[93];
   armbbbar[71]=1./2.*armbbbar[75] + armbbbar[84] + 1./2.*armbbbar[98]
    + 3./2.*armbbbar[125] + 3./2.*armbbbar[74] + 3./2.*armbbbar[71] + 5.
   /2.*armbbbar[124] + armbbbar[113] + 1./2.*armbbbar[122];
   armbbbar[74]=armbbbar[73]*armbbbar[290];
   armbbbar[74]=1./8.*armbbbar[93] + armbbbar[74];
   armbbbar[74]=MMZ*armbbbar[74];
   armbbbar[71]=1./2.*armbbbar[71] + 9*armbbbar[74];
   armbbbar[71]=MMZ*armbbbar[71];
   armbbbar[74]=armbbbar[11]*armbbbar[93];
   armbbbar[75]= - 7./4.*armbbbar[55] + 5*armbbbar[74];
   armbbbar[84]= - armbbbar[12]*armbbbar[93];
   armbbbar[98]= - armbbbar[13]*armbbbar[93];
   armbbbar[75]=2*armbbbar[98] + 1./2.*armbbbar[75] + armbbbar[84];
   armbbbar[75]=armbbbar[13]*armbbbar[75];
   armbbbar[84]=3*armbbbar[112] + 55./18. + armbbbar[43];
   armbbbar[84]=armbbbar[297]*armbbbar[84];
   armbbbar[98]= - 13./8.*armbbbar[16] + armbbbar[32] + 5./8.*
   armbbbar[27];
   armbbbar[98]=armbbbar[55]*armbbbar[98];
   armbbbar[74]= - 5./12.*armbbbar[55] + armbbbar[74];
   armbbbar[74]=armbbbar[12]*armbbbar[74];
   armbbbar[61]=armbbbar[71] + 1./16.*armbbbar[84] + armbbbar[75] + 1./
   2.*armbbbar[74] + 13./24.*armbbbar[246] + armbbbar[198] + 7./8.*
   armbbbar[226] + 5./24.*armbbbar[214] + 1./3.*armbbbar[98] + 
   armbbbar[291] - 3./2.*armbbbar[43] + armbbbar[268] - 29./12.*
   armbbbar[70] + armbbbar[61] - 98./9. + armbbbar[107];
   armbbbar[61]=MMZ*armbbbar[61];
   armbbbar[71]= - 13./8. - 4*armbbbar[47];
   armbbbar[71]=4./3.*armbbbar[289] + 5./6.*armbbbar[247] + 1./3.*
   armbbbar[71] + 3./2.*armbbbar[128];
   armbbbar[71]=armbbbar[13]*armbbbar[71];
   armbbbar[74]=59./6.*armbbbar[25] - 59./6.*armbbbar[17] + 47./6.*
   armbbbar[27] + armbbbar[33] - 53./6.*armbbbar[32];
   armbbbar[75]=armbbbar[12]*armbbbar[43];
   armbbbar[84]= - armbbbar[13]*armbbbar[41];
   armbbbar[74]=17./6.*armbbbar[84] + 1./3.*armbbbar[74] + 1./2.*
   armbbbar[75];
   armbbbar[74]=armbbbar[297]*armbbbar[74];
   armbbbar[75]= - 403./12. + armbbbar[196];
   armbbbar[75]=1./3.*armbbbar[75] + 19*armbbbar[43];
   armbbbar[75]=armbbbar[12]*armbbbar[75];
   armbbbar[84]=1./9. + armbbbar[41];
   armbbbar[84]=armbbbar[13]*armbbbar[84];
   armbbbar[74]=1./2.*armbbbar[74] + 1./2.*armbbbar[84] + 9./4.*
   armbbbar[118] + 1./4.*armbbbar[75] - 23./9.*armbbbar[11] - 49./18.*
   armbbbar[25] - 1./3.*armbbbar[17] - 55./36.*armbbbar[16] - 1./6.*
   armbbbar[27] + 1./4.*armbbbar[33] + 13./9.*armbbbar[32];
   armbbbar[74]=armbbbar[297]*armbbbar[74];
   armbbbar[75]= - 3*armbbbar[45];
   armbbbar[84]=59./12. + armbbbar[75];
   armbbbar[84]=2./3.*armbbbar[128] + 1./2.*armbbbar[84] + 8./3.*
   armbbbar[47];
   armbbbar[84]=armbbbar[12]*armbbbar[84];
   armbbbar[61]=armbbbar[69] + armbbbar[61] + 1./8.*armbbbar[74] + 
   armbbbar[71] + armbbbar[84] + 1./8.*armbbbar[11] - 13./24.*
   armbbbar[25] + 67./24.*armbbbar[17] + 1./24.*armbbbar[16] + 7./24.*
   armbbbar[27] + 2./3.*armbbbar[32] - armbbbar[33] - 13./12.*
   armbbbar[35];
   armbbbar[61]=armbbbar[52]*armbbbar[61];
   armbbbar[69]=armbbbar[95] + armbbbar[101] + armbbbar[92] + 
   armbbbar[94] + armbbbar[78] + armbbbar[76] + armbbbar[105] + 35./4.*
   armbbbar[45] + armbbbar[79] + 2629./72. + armbbbar[103];
   armbbbar[69]=armbbbar[11]*armbbbar[69];
   armbbbar[71]=armbbbar[130] - 9*armbbbar[33];
   armbbbar[71]=9*armbbbar[17] + 45./4.*armbbbar[16] + armbbbar[151] + 
   1./2.*armbbbar[71] + armbbbar[32];
   armbbbar[74]=19./2.*armbbbar[43] + 353./24. + armbbbar[105];
   armbbbar[74]=armbbbar[12]*armbbbar[74];
   armbbbar[69]=armbbbar[153] + 1./4.*armbbbar[74] + 1./4.*armbbbar[69]
    + armbbbar[159] + 1./2.*armbbbar[71] - armbbbar[25];
   armbbbar[69]=armbbbar[3]*armbbbar[69];
   armbbbar[71]= - 17./2.*armbbbar[43];
   armbbbar[74]=armbbbar[71] + armbbbar[105] + 281./24. + armbbbar[45];
   armbbbar[74]=1./2.*armbbbar[74] + 17./3.*armbbbar[47];
   armbbbar[74]=armbbbar[233] + 1./2.*armbbbar[74] + armbbbar[193];
   armbbbar[74]=armbbbar[11]*armbbbar[74];
   armbbbar[78]= - armbbbar[33] + 5./9.*armbbbar[180] + armbbbar[212];
   armbbbar[78]=armbbbar[17] + 175./12.*armbbbar[16] + 223./72.*
   armbbbar[24] + 1./2.*armbbbar[78] - 17./3.*armbbbar[32];
   armbbbar[84]= - 17./3. + 1./2.*armbbbar[43];
   armbbbar[84]=armbbbar[12]*armbbbar[84];
   armbbbar[94]= - armbbbar[3]*armbbbar[173];
   armbbbar[74]=9./16.*armbbbar[94] + 1./2.*armbbbar[84] + armbbbar[74]
    + armbbbar[159] + 1./2.*armbbbar[78] + 17./3.*armbbbar[25];
   armbbbar[74]=armbbbar[3]*armbbbar[74];
   armbbbar[59]=25./6.*armbbbar[59] - 33791./7776. + armbbbar[58];
   armbbbar[59]= - 55./216.*armbbbar[45] + 125./144.*armbbbar[131] - 
   119./72.*armbbbar[41] + 1./2.*armbbbar[59] - armbbbar[44];
   armbbbar[78]=armbbbar[57] - 5./24.*armbbbar[54];
   armbbbar[84]=armbbbar[16]*armbbbar[78];
   armbbbar[78]=armbbbar[11]*armbbbar[78];
   armbbbar[94]=431./24. + 17*armbbbar[41];
   armbbbar[94]=armbbbar[3]*armbbbar[94];
   armbbbar[94]= - 125./96.*armbbbar[54] + armbbbar[94];
   armbbbar[94]=armbbbar[13]*armbbbar[94];
   armbbbar[95]=1./3.*armbbbar[42];
   armbbbar[59]=1./3.*armbbbar[94] + armbbbar[74] + 25./24.*
   armbbbar[78] + 25./48.*armbbbar[287] - 275./1296.*armbbbar[9] + 
   armbbbar[95] + 125./288.*armbbbar[284] - 139./72.*armbbbar[39] - 935.
   /1296.*armbbbar[47] + 527./864.*armbbbar[43] + 25./24.*armbbbar[84]
    + 25./48.*armbbbar[68] + 1./4.*armbbbar[59] + armbbbar[40];
   armbbbar[68]=armbbbar[48] + 49./18.*armbbbar[47] + armbbbar[296] - 
   275./72. + armbbbar[219];
   armbbbar[68]=armbbbar[3]*armbbbar[68];
   armbbbar[74]=armbbbar[13]*armbbbar[222];
   armbbbar[68]=1./2.*armbbbar[68] + 3*armbbbar[74];
   armbbbar[68]=1./2.*armbbbar[68] + 3*armbbbar[216];
   armbbbar[68]=MMt*armbbbar[68];
   armbbbar[74]= - 1./3. - 1./2.*armbbbar[43];
   armbbbar[74]=armbbbar[11]*armbbbar[74];
   armbbbar[74]=armbbbar[74] + 1./3.*armbbbar[12];
   armbbbar[74]=armbbbar[3]*armbbbar[74];
   armbbbar[74]=31./216.*armbbbar[43] + armbbbar[74];
   armbbbar[74]=armbbbar[297]*armbbbar[74];
   armbbbar[59]=1./8.*armbbbar[74] + 1./2.*armbbbar[59] + armbbbar[68];
   armbbbar[59]=armbbbar[297]*armbbbar[59];
   armbbbar[58]=31./54.*armbbbar[229] + armbbbar[62] - 187739./7776. + 
   armbbbar[58];
   armbbbar[58]= - 983./432.*armbbbar[45] + 125./144.*armbbbar[67] + 
   armbbbar[126] + 1./2.*armbbbar[58] - armbbbar[44];
   armbbbar[62]=armbbbar[167] + armbbbar[199];
   armbbbar[68]=armbbbar[12]*armbbbar[26];
   armbbbar[74]=armbbbar[62] + armbbbar[68];
   armbbbar[78]=armbbbar[221] - 119./9.*armbbbar[47] - armbbbar[41] + 
   461./72. + armbbbar[219];
   armbbbar[78]=armbbbar[3]*armbbbar[78];
   armbbbar[74]=19./48.*armbbbar[74] + armbbbar[78];
   armbbbar[78]=19./96.*armbbbar[26] + armbbbar[223];
   armbbbar[78]=armbbbar[13]*armbbbar[78];
   armbbbar[74]=armbbbar[228] + 1./2.*armbbbar[74] + armbbbar[78];
   armbbbar[74]=MMt*armbbbar[74];
   armbbbar[78]=125./3.*armbbbar[54] + 19*armbbbar[199];
   armbbbar[78]=1./3.*armbbbar[78] + 27*armbbbar[213];
   armbbbar[84]= - 161./72. - armbbbar[41];
   armbbbar[84]=armbbbar[3]*armbbbar[84];
   armbbbar[78]=1./32.*armbbbar[78] + armbbbar[84];
   armbbbar[78]=armbbbar[13]*armbbbar[78];
   armbbbar[84]=armbbbar[227] + 125./72.*armbbbar[54];
   armbbbar[94]=armbbbar[16]*armbbbar[84];
   armbbbar[84]=armbbbar[11]*armbbbar[84];
   armbbbar[98]=armbbbar[12]*armbbbar[199];
   armbbbar[58]=armbbbar[59] + armbbbar[74] + armbbbar[78] + 
   armbbbar[69] + 31./48.*armbbbar[98] + 1./8.*armbbbar[84] + 
   armbbbar[288] - 295./2592.*armbbbar[9] + armbbbar[286] + 
   armbbbar[285] + 125./288.*armbbbar[72] + armbbbar[276] + 449./2592.*
   armbbbar[47] - 67./72.*armbbbar[43] + 1./8.*armbbbar[94] + 
   armbbbar[273] + 1./4.*armbbbar[58] + armbbbar[40];
   armbbbar[58]=armbbbar[297]*armbbbar[58];
   armbbbar[59]= - 1./4.*armbbbar[45];
   armbbbar[65]=armbbbar[132] + armbbbar[85] + 125./18.*armbbbar[110]
    + armbbbar[90] + armbbbar[101] + armbbbar[100] + armbbbar[81] + 
   armbbbar[82] + armbbbar[76] + armbbbar[59] + armbbbar[65] + 
   armbbbar[140] + 21407./648. + armbbbar[103];
   armbbbar[69]=armbbbar[224] + 15./4.*armbbbar[12];
   armbbbar[69]=armbbbar[3]*armbbbar[69];
   armbbbar[65]=1./2.*armbbbar[65] + armbbbar[69];
   armbbbar[69]=25./72.*armbbbar[54] + armbbbar[133];
   armbbbar[69]=armbbbar[13]*armbbbar[69];
   armbbbar[65]=1./2.*armbbbar[65] + 5*armbbbar[69];
   armbbbar[65]=armbbbar[13]*armbbbar[65];
   armbbbar[69]=armbbbar[234] + armbbbar[233] + armbbbar[194] + 
   armbbbar[97] + armbbbar[193] + armbbbar[232] + armbbbar[231] + 
   armbbbar[230] + 1./4.*armbbbar[45] + armbbbar[235] + armbbbar[136]
    + armbbbar[186] - 2107./72. + armbbbar[229];
   armbbbar[69]=armbbbar[3]*armbbbar[69];
   armbbbar[72]=43./3.*armbbbar[26] + armbbbar[238];
   armbbbar[72]=armbbbar[12]*armbbbar[72];
   armbbbar[72]=armbbbar[72] + armbbbar[237] + armbbbar[184] + 19./3.*
   armbbbar[109];
   armbbbar[69]=1./32.*armbbbar[72] + armbbbar[69];
   armbbbar[69]=armbbbar[13]*armbbbar[69];
   armbbbar[72]=15115./216. + armbbbar[103];
   armbbbar[64]=armbbbar[64] + armbbbar[244] - armbbbar[42] + 
   armbbbar[243] - armbbbar[39] + 1763./81.*armbbbar[47] + 
   armbbbar[242] - 1./12.*armbbbar[45] + armbbbar[241] + 1./3.*
   armbbbar[72] + armbbbar[240];
   armbbbar[63]=armbbbar[195] + 3./2.*armbbbar[63] + armbbbar[109];
   armbbbar[63]=armbbbar[12]*armbbbar[63];
   armbbbar[63]=1./4.*armbbbar[64] + armbbbar[63];
   armbbbar[64]= - armbbbar[48] - 23./8. + armbbbar[219];
   armbbbar[64]=armbbbar[12]*armbbbar[64];
   armbbbar[64]=armbbbar[152] + 3./2.*armbbbar[64];
   armbbbar[64]=armbbbar[3]*armbbbar[64];
   armbbbar[63]=1./2.*armbbbar[63] + armbbbar[64];
   armbbbar[63]=armbbbar[187] + 1./2.*armbbbar[63] + armbbbar[69];
   armbbbar[63]=MMt*armbbbar[63];
   armbbbar[59]=armbbbar[90] + armbbbar[101] + armbbbar[92] + 
   armbbbar[81] + armbbbar[82] + armbbbar[76] + armbbbar[230] + 
   armbbbar[59] + armbbbar[136] + 1381./72. + armbbbar[103];
   armbbbar[59]=armbbbar[8]*armbbbar[59];
   armbbbar[64]=10303./216. + armbbbar[179];
   armbbbar[64]=armbbbar[192] + armbbbar[188] + armbbbar[191] + 
   armbbbar[190] + armbbbar[189] - 3./16.*armbbbar[40] - 3./8.*
   armbbbar[38] - 3./8.*armbbbar[44] - 3./32.*armbbbar[50] + 1./16.*
   armbbbar[64] + 7./3.*armbbbar[229];
   armbbbar[64]=armbbbar[12]*armbbbar[64];
   armbbbar[69]= - 79./24. + armbbbar[200];
   armbbbar[69]=armbbbar[168] + 1./2.*armbbbar[69] + armbbbar[99];
   armbbbar[69]=armbbbar[11]*armbbbar[69];
   armbbbar[72]= - armbbbar[32] + armbbbar[16];
   armbbbar[72]=1./2.*armbbbar[72] + armbbbar[25];
   armbbbar[69]=7./3.*armbbbar[72] + 1./2.*armbbbar[69];
   armbbbar[69]=armbbbar[3]*armbbbar[69];
   armbbbar[72]=47./8. + armbbbar[196];
   armbbbar[74]=armbbbar[3]*armbbbar[11];
   armbbbar[72]=9./2.*armbbbar[74] + 17./2.*armbbbar[43] + 1./3.*
   armbbbar[72] + armbbbar[230];
   armbbbar[72]=armbbbar[13]*armbbbar[3]*armbbbar[72];
   armbbbar[76]= - 7*armbbbar[41];
   armbbbar[78]= - 97./18. + armbbbar[76];
   armbbbar[81]= - 17*armbbbar[43];
   armbbbar[78]=7./3.*armbbbar[78] + armbbbar[81];
   armbbbar[78]= - 31./3.*armbbbar[48] + 1./2.*armbbbar[78] - 385./27.*
   armbbbar[47];
   armbbbar[82]=11 + armbbbar[196];
   armbbbar[82]=1./3.*armbbbar[82] + armbbbar[48];
   armbbbar[84]=armbbbar[3]*armbbbar[204];
   armbbbar[82]=1./2.*armbbbar[82] + 9*armbbbar[84];
   armbbbar[82]=MMt*armbbbar[3]*armbbbar[82];
   armbbbar[69]=armbbbar[82] + armbbbar[72] + 1./48.*armbbbar[78] + 
   armbbbar[69];
   armbbbar[69]=MMt*armbbbar[69];
   armbbbar[72]= - armbbbar[173] + armbbbar[117];
   armbbbar[72]=armbbbar[3]*armbbbar[72];
   armbbbar[78]= - 1./3. - 1./4.*armbbbar[43];
   armbbbar[78]=armbbbar[8]*armbbbar[78];
   armbbbar[82]=31./24. + armbbbar[43];
   armbbbar[82]=armbbbar[11]*armbbbar[82];
   armbbbar[84]= - armbbbar[13]*armbbbar[43];
   armbbbar[85]= - MMH*armbbbar[42];
   armbbbar[72]=1./4.*armbbbar[84] + 1./2.*armbbbar[72] + 41./216.*
   armbbbar[12] + 1./3.*armbbbar[85] + armbbbar[78] + 1./9.*
   armbbbar[82];
   armbbbar[78]=1 + armbbbar[43];
   armbbbar[78]=armbbbar[13]*armbbbar[3]*armbbbar[78];
   armbbbar[82]=MMt*armbbbar[177];
   armbbbar[78]=armbbbar[82] - 1./48.*armbbbar[43] + armbbbar[78];
   armbbbar[78]=MMt*armbbbar[78];
   armbbbar[72]=1./2.*armbbbar[72] + armbbbar[78];
   armbbbar[72]=armbbbar[297]*armbbbar[72];
   armbbbar[76]= - 4051./108. + armbbbar[76];
   armbbbar[74]=227./3.*armbbbar[74] + 125./9.*armbbbar[96] + 1./3.*
   armbbbar[76] + armbbbar[81];
   armbbbar[76]= - 125./24.*armbbbar[54] - 17*armbbbar[3];
   armbbbar[76]=armbbbar[13]*armbbbar[76];
   armbbbar[74]=1./8.*armbbbar[74] + 1./3.*armbbbar[76];
   armbbbar[74]=armbbbar[13]*armbbbar[74];
   armbbbar[76]=401./324. + armbbbar[179];
   armbbbar[76]=1./2.*armbbbar[76] + armbbbar[136];
   armbbbar[76]=17./9.*armbbbar[43] - 3./4.*armbbbar[40] + 1./2.*
   armbbbar[76] - 1./9.*armbbbar[45];
   armbbbar[76]=1./2.*armbbbar[76] - 17./27.*armbbbar[47];
   armbbbar[76]=25./96.*armbbbar[183] + 25./96.*armbbbar[182] - 5./108.
   *armbbbar[9] + 1./4.*armbbbar[42] + 1./4.*armbbbar[76] - 2./3.*
   armbbbar[39];
   armbbbar[76]=armbbbar[11]*armbbbar[76];
   armbbbar[71]=armbbbar[71] - 17./8. + armbbbar[230];
   armbbbar[71]=armbbbar[8]*armbbbar[71];
   armbbbar[78]=85*armbbbar[29] + 37*armbbbar[34] - 107./9.*armbbbar[6]
   ;
   armbbbar[78]=841./192.*armbbbar[32] + 1./16.*armbbbar[78] + 
   armbbbar[33];
   armbbbar[78]= - 389./144.*armbbbar[25] - 43./96.*armbbbar[17] - 9745.
   /1728.*armbbbar[16] - 61./24.*armbbbar[24] + 1./3.*armbbbar[78] - 7./
   32.*armbbbar[27];
   armbbbar[71]=1./3.*armbbbar[78] + 1./8.*armbbbar[71];
   armbbbar[78]= - 17./9.*armbbbar[8] - 27./2.*armbbbar[11];
   armbbbar[78]=armbbbar[11]*armbbbar[78];
   armbbbar[78]=armbbbar[78] + 39*armbbbar[264];
   armbbbar[78]=armbbbar[3]*armbbbar[78];
   armbbbar[81]= - 49./54.*armbbbar[42] + 59./27.*armbbbar[39] + 5./9.
    + armbbbar[166];
   armbbbar[81]=MMH*armbbbar[81];
   armbbbar[69]=1./4.*armbbbar[72] + 1./2.*armbbbar[69] + 1./4.*
   armbbbar[74] + 1./16.*armbbbar[78] + 349./576.*armbbbar[12] + 1./8.*
   armbbbar[81] + 1./2.*armbbbar[71] + armbbbar[76];
   armbbbar[69]=armbbbar[297]*armbbbar[69];
   armbbbar[71]= - 1235./36. + armbbbar[179];
   armbbbar[71]=armbbbar[272] + armbbbar[271] - 1./54.*armbbbar[9] + 
   armbbbar[263] + armbbbar[259] - 11./6.*armbbbar[39] + 23./54.*
   armbbbar[47] - 17./12.*armbbbar[43] - 3./8.*armbbbar[40] - 17./36.*
   armbbbar[45] - 3./4.*armbbbar[44] + 1./8.*armbbbar[71] + 1./9.*
   armbbbar[229];
   armbbbar[71]=armbbbar[11]*armbbbar[71];
   armbbbar[72]=137./32.*armbbbar[35] + 7./3.*armbbbar[33] + 
   armbbbar[275] - 139./432.*armbbbar[6] + 2./9.*armbbbar[18] + 
   armbbbar[274];
   armbbbar[74]= - 17./2.*armbbbar[12] - armbbbar[8] + 5*armbbbar[11];
   armbbbar[74]=armbbbar[12]*armbbbar[74];
   armbbbar[74]=armbbbar[277] + armbbbar[74];
   armbbbar[74]=armbbbar[3]*armbbbar[74];
   armbbbar[59]=armbbbar[69] + armbbbar[63] + 1./2.*armbbbar[65] + 1./8.
   *armbbbar[74] + armbbbar[64] + armbbbar[185] + 1./2.*armbbbar[71] + 
   1./8.*armbbbar[59] - 17./18.*armbbbar[25] - 1019./288.*armbbbar[17]
    - 193./162.*armbbbar[16] + armbbbar[283] + armbbbar[282] - 697./576.
   *armbbbar[27] + 241./576.*armbbbar[32] + 1./27.*armbbbar[281] + 
   armbbbar[279] + 1./3.*armbbbar[72] + armbbbar[278];
   armbbbar[59]=armbbbar[297]*armbbbar[59];
   armbbbar[63]= - 5./18.*armbbbar[35] - armbbbar[31];
   armbbbar[63]=7./2.*armbbbar[63] - 5./9.*armbbbar[32];
   armbbbar[63]=armbbbar[203] + 877./18.*armbbbar[25] + armbbbar[236]
    + 29./6.*armbbbar[16] + armbbbar[134] + armbbbar[88] + 5*
   armbbbar[63] - 37./18.*armbbbar[27];
   armbbbar[64]=419./8. + 73*armbbbar[47];
   armbbbar[64]=1./9.*armbbbar[64] - armbbbar[48];
   armbbbar[64]=armbbbar[11]*armbbbar[64];
   armbbbar[63]=armbbbar[141] + armbbbar[248] + armbbbar[245] + 1./2.*
   armbbbar[63] + 1./3.*armbbbar[64];
   armbbbar[64]=armbbbar[140] + 127./48. + armbbbar[139];
   armbbbar[64]=armbbbar[249] + armbbbar[80] + armbbbar[143] + 1./2.*
   armbbbar[64] + armbbbar[79];
   armbbbar[64]=armbbbar[252] + armbbbar[251] - armbbbar[42] + 
   armbbbar[250] + 1./2.*armbbbar[64] - armbbbar[39];
   armbbbar[64]=armbbbar[13]*armbbbar[64];
   armbbbar[65]=armbbbar[254] + armbbbar[253] + armbbbar[148] + 
   armbbbar[147] + 869./384. + armbbbar[145];
   armbbbar[65]=armbbbar[157] + armbbbar[77] + 1./2.*armbbbar[65] + 
   armbbbar[256];
   armbbbar[65]=MMt*armbbbar[65];
   armbbbar[63]=armbbbar[65] + 1./4.*armbbbar[63] + armbbbar[64];
   armbbbar[63]=MMt*armbbbar[63];
   armbbbar[64]= - 11./9.*armbbbar[35] + armbbbar[197];
   armbbbar[65]=armbbbar[48] + 13./8. + 3*armbbbar[46];
   armbbbar[65]=armbbbar[8]*armbbbar[65];
   armbbbar[64]=3./2.*armbbbar[65] + 2441./108.*armbbbar[25] + 149./36.
   *armbbbar[17] + 401./216.*armbbbar[16] - 5./8.*armbbbar[23] + 75./8.
   *armbbbar[24] - 11./8.*armbbbar[27] + 7./4.*armbbbar[64] - 13./27.*
   armbbbar[32];
   armbbbar[65]=armbbbar[138] + armbbbar[95];
   armbbbar[65]=MMH*armbbbar[65];
   armbbbar[69]= - 11./8. - 7./3.*armbbbar[47];
   armbbbar[69]=1./3.*armbbbar[69] - armbbbar[48];
   armbbbar[69]=armbbbar[11]*armbbbar[69];
   armbbbar[71]= - 7./6. + armbbbar[92];
   armbbbar[71]=armbbbar[12]*armbbbar[71];
   armbbbar[64]=1./2.*armbbbar[141] + 1./4.*armbbbar[71] + 1./4.*
   armbbbar[65] + 1./2.*armbbbar[64] + 1./3.*armbbbar[69];
   armbbbar[65]=armbbbar[140] + 731./432. + armbbbar[139];
   armbbbar[65]=3./4.*armbbbar[48] + armbbbar[80] + armbbbar[143] + 1./
   2.*armbbbar[65] + armbbbar[79];
   armbbbar[69]=armbbbar[144] + 31*armbbbar[11];
   armbbbar[69]=1./4.*armbbbar[69] + 19*armbbbar[12];
   armbbbar[69]=armbbbar[3]*armbbbar[69];
   armbbbar[65]=3./4.*armbbbar[102] + 1./2.*armbbbar[69] + 1./4.*
   armbbbar[65] - armbbbar[42];
   armbbbar[65]=armbbbar[13]*armbbbar[65];
   armbbbar[69]= - 11./36.*armbbbar[48] + armbbbar[148] + armbbbar[147]
    + 11245./3456. + armbbbar[145];
   armbbbar[71]= - 1 - armbbbar[48];
   armbbbar[71]=armbbbar[12]*armbbbar[71];
   armbbbar[72]=armbbbar[11]*armbbbar[48];
   armbbbar[71]=3./2.*armbbbar[71] + armbbbar[150] + 3./4.*armbbbar[72]
   ;
   armbbbar[71]=armbbbar[3]*armbbbar[71];
   armbbbar[69]=1./2.*armbbbar[69] + armbbbar[71];
   armbbbar[71]= - 3./8. + armbbbar[46];
   armbbbar[71]=3*armbbbar[71] - armbbbar[50];
   armbbbar[74]= - 1./2.*armbbbar[48];
   armbbbar[71]=armbbbar[74] + 1./2.*armbbbar[71] + armbbbar[135];
   armbbbar[71]=armbbbar[3]*armbbbar[71];
   armbbbar[71]=armbbbar[71] + 3./2.*armbbbar[156];
   armbbbar[71]=armbbbar[13]*armbbbar[71];
   armbbbar[69]=armbbbar[158] + 1./2.*armbbbar[69] + 3*armbbbar[71];
   armbbbar[69]=MMt*armbbbar[69];
   armbbbar[64]=armbbbar[69] + 1./4.*armbbbar[64] + armbbbar[65];
   armbbbar[64]=MMt*armbbbar[64];
   armbbbar[65]=armbbbar[42] + armbbbar[257] + armbbbar[255] + 13./3.
    + armbbbar[239];
   armbbbar[65]=armbbbar[8]*armbbbar[65];
   armbbbar[69]= - 17./2.*armbbbar[42] + 241./24. + 13*armbbbar[39];
   armbbbar[69]=armbbbar[11]*armbbbar[69];
   armbbbar[71]=31./432.*armbbbar[42] - 55./432.*armbbbar[39] + 3./16.*
   armbbbar[40] + 3./8.*armbbbar[44] - 19./9. + 27./16.*armbbbar[37];
   armbbbar[71]=MMH*armbbbar[71];
   armbbbar[65]=1./2.*armbbbar[71] + 1./36.*armbbbar[69] + 1./4.*
   armbbbar[65] + 7./32.*armbbbar[25] + 9./32.*armbbbar[17] + 1./6.*
   armbbbar[16] + 19./32.*armbbbar[23] + 65./36.*armbbbar[24] - 13./32.
   *armbbbar[31] - 7./16.*armbbbar[30] - 37./288.*armbbbar[29] + 13./
   144.*armbbbar[34] - 1./9.*armbbbar[36] - 15./32.*armbbbar[28];
   armbbbar[65]=MMH*armbbbar[65];
   armbbbar[69]= - 47*armbbbar[8] + 161./18.*armbbbar[11];
   armbbbar[69]=armbbbar[11]*armbbbar[69];
   armbbbar[69]=armbbbar[149] + 1./9.*armbbbar[69];
   armbbbar[71]= - 65*armbbbar[8] + 437./9.*armbbbar[11];
   armbbbar[76]=293./216. + armbbbar[42];
   armbbbar[76]=MMH*armbbbar[76];
   armbbbar[71]=419./36.*armbbbar[12] + 1./8.*armbbbar[71] + 
   armbbbar[76];
   armbbbar[71]=armbbbar[12]*armbbbar[71];
   armbbbar[76]= - 9*armbbbar[8] - 139./6.*armbbbar[11];
   armbbbar[77]=armbbbar[164] + armbbbar[42];
   armbbbar[77]=MMH*armbbbar[77];
   armbbbar[76]=2141./432.*armbbbar[13] - 257./18.*armbbbar[12] + 1./4.
   *armbbbar[76] + armbbbar[77];
   armbbbar[76]=armbbbar[13]*armbbbar[76];
   armbbbar[65]=1./4.*armbbbar[76] + 1./4.*armbbbar[71] + 1./8.*
   armbbbar[69] + armbbbar[65];
   armbbbar[69]= - armbbbar[11]*armbbbar[42];
   armbbbar[69]=31./108.*armbbbar[162] + 31./108.*armbbbar[8] + 
   armbbbar[69];
   armbbbar[69]=MMH*armbbbar[69];
   armbbbar[71]=armbbbar[111] + 1./9.*armbbbar[11];
   armbbbar[71]=armbbbar[11]*armbbbar[71];
   armbbbar[76]= - 31./108.*MMH + armbbbar[210] + 7./9.*armbbbar[11];
   armbbbar[76]=armbbbar[12]*armbbbar[76];
   armbbbar[77]= - armbbbar[11] + armbbbar[12];
   armbbbar[78]=armbbbar[13]*armbbbar[77];
   armbbbar[69]=1./4.*armbbbar[78] + 1./2.*armbbbar[76] + armbbbar[71]
    + 1./2.*armbbbar[69];
   armbbbar[71]=armbbbar[208] - armbbbar[12];
   armbbbar[71]=armbbbar[3]*armbbbar[71];
   armbbbar[71]= - 31./72. + armbbbar[71];
   armbbbar[71]=armbbbar[13]*armbbbar[71];
   armbbbar[72]=armbbbar[3]*armbbbar[72];
   armbbbar[72]= - 31./36.*armbbbar[48] + 3*armbbbar[72];
   armbbbar[72]=MMt*armbbbar[72];
   armbbbar[71]=1./2.*armbbbar[72] + 1./48.*armbbbar[77] + armbbbar[71]
   ;
   armbbbar[71]=MMt*armbbbar[71];
   armbbbar[69]=1./2.*armbbbar[69] + armbbbar[71];
   armbbbar[69]=armbbbar[297]*armbbbar[69];
   armbbbar[64]=1./4.*armbbbar[69] + 1./2.*armbbbar[65] + armbbbar[64];
   armbbbar[64]=armbbbar[297]*armbbbar[64];
   armbbbar[65]= - 43*armbbbar[8] + 245./18.*armbbbar[11];
   armbbbar[65]=391./72.*armbbbar[13] + armbbbar[265] + 1./6.*
   armbbbar[65] + armbbbar[163];
   armbbbar[65]=armbbbar[13]*armbbbar[65];
   armbbbar[69]=armbbbar[267] + 27*armbbbar[11];
   armbbbar[69]= - 19./2.*armbbbar[12] + 1./4.*armbbbar[69] + 
   armbbbar[269];
   armbbbar[69]=armbbbar[12]*armbbbar[69];
   armbbbar[71]=armbbbar[8] - 59./108.*armbbbar[11];
   armbbbar[71]=armbbbar[11]*armbbbar[71];
   armbbbar[63]=armbbbar[64] + armbbbar[63] + 1./8.*armbbbar[65] + 1./8.
   *armbbbar[69] + armbbbar[205] + armbbbar[270] + 1./3.*armbbbar[71];
   armbbbar[63]=armbbbar[297]*armbbbar[63];
   armbbbar[64]=armbbbar[8] + 1./2.*armbbbar[162];
   armbbbar[64]=MMH*armbbbar[64];
   armbbbar[64]=armbbbar[64] + 3./2.*armbbbar[120];
   armbbbar[64]=armbbbar[12]*armbbbar[64];
   armbbbar[64]=armbbbar[114] + armbbbar[119] + armbbbar[64];
   armbbbar[64]=armbbbar[297]*armbbbar[64];
   armbbbar[64]=armbbbar[258] + 1./2.*armbbbar[64];
   armbbbar[64]=armbbbar[52]*armbbbar[297]*armbbbar[64];
   armbbbar[65]=armbbbar[160] - 11./18.*armbbbar[11];
   armbbbar[69]=3./4.*armbbbar[13];
   armbbbar[65]=armbbbar[69] - armbbbar[12] + 1./4.*armbbbar[65] + 
   armbbbar[162];
   armbbbar[65]=armbbbar[13]*armbbbar[65];
   armbbbar[71]= - 13./3.*armbbbar[11] + 1./2.*MMH;
   armbbbar[71]=armbbbar[12]*armbbbar[71];
   armbbbar[72]= - armbbbar[8] + armbbbar[85];
   armbbbar[72]=MMH*armbbbar[72];
   armbbbar[71]=1./2.*armbbbar[72] + armbbbar[71];
   armbbbar[65]=1./24.*armbbbar[71] + armbbbar[65];
   armbbbar[71]=armbbbar[8]*armbbbar[48];
   armbbbar[76]=13./18.*armbbbar[266] + 3*armbbbar[71];
   armbbbar[77]=1./3.*armbbbar[154];
   armbbbar[76]=1./4.*armbbbar[76] + armbbbar[77];
   armbbbar[78]=1./2. + armbbbar[97];
   armbbbar[78]=1./8.*armbbbar[78] + 3*armbbbar[146];
   armbbbar[78]=armbbbar[13]*armbbbar[78];
   armbbbar[79]=armbbbar[13]*armbbbar[178];
   armbbbar[79]=1./16.*armbbbar[48] + 3*armbbbar[79];
   armbbbar[79]=MMt*armbbbar[79];
   armbbbar[76]=armbbbar[79] + 1./2.*armbbbar[76] + armbbbar[78];
   armbbbar[76]=MMt*armbbbar[76];
   armbbbar[65]=1./2.*armbbbar[65] + armbbbar[76];
   armbbbar[65]=MMt*armbbbar[65];
   armbbbar[72]=armbbbar[72] + armbbbar[116];
   armbbbar[69]=armbbbar[69] - armbbbar[12] + armbbbar[162] + 7./4.*
   armbbbar[8] - 1./3.*armbbbar[11];
   armbbbar[69]=armbbbar[13]*armbbbar[69];
   armbbbar[69]=1./48.*armbbbar[72] + armbbbar[69];
   armbbbar[71]=3./4.*armbbbar[71] + armbbbar[77];
   armbbbar[71]=armbbbar[79] + 1./2.*armbbbar[71] + armbbbar[78];
   armbbbar[71]=MMt*armbbbar[71];
   armbbbar[69]=1./2.*armbbbar[69] + armbbbar[71];
   armbbbar[69]=MMt*armbbbar[69];
   armbbbar[71]= - armbbbar[8]*armbbbar[42];
   armbbbar[76]=armbbbar[11]*armbbbar[42];
   armbbbar[76]=1./4.*armbbbar[71] + 1./9.*armbbbar[76];
   armbbbar[76]=MMH*armbbbar[76];
   armbbbar[76]=armbbbar[76] + armbbbar[262] + 1./9.*armbbbar[171];
   armbbbar[76]=MMH*armbbbar[76];
   armbbbar[72]=armbbbar[13]*armbbbar[72];
   armbbbar[77]=armbbbar[206] - 1./9.*armbbbar[11];
   armbbbar[77]=armbbbar[12]*MMH*armbbbar[77];
   armbbbar[76]=1./4.*armbbbar[72] + armbbbar[76] + armbbbar[77];
   armbbbar[69]=1./4.*armbbbar[76] + armbbbar[69];
   armbbbar[69]=armbbbar[297]*armbbbar[69];
   armbbbar[71]=armbbbar[172] + 1./9.*armbbbar[260] + 1./8.*
   armbbbar[71];
   armbbbar[71]=MMH*armbbbar[71];
   armbbbar[76]= - 1./3.*armbbbar[8] + armbbbar[108];
   armbbbar[76]=armbbbar[11]*armbbbar[76];
   armbbbar[71]=armbbbar[71] - 1./8.*armbbbar[142] + 1./3.*armbbbar[76]
   ;
   armbbbar[71]=MMH*armbbbar[71];
   armbbbar[71]=1./8.*armbbbar[72] + armbbbar[71] + 1./2.*armbbbar[77];
   armbbbar[65]=armbbbar[69] + 1./2.*armbbbar[71] + armbbbar[65];
   armbbbar[65]=armbbbar[297]*armbbbar[65];
   armbbbar[65]=1./2.*armbbbar[261] + armbbbar[65];
   armbbbar[65]=armbbbar[4]*armbbbar[297]*armbbbar[65];
   armbbbar[69]=armbbbar[11] + armbbbar[13];
   armbbbar[69]=armbbbar[13]*armbbbar[69];
   armbbbar[71]=2*MMt + armbbbar[86] + armbbbar[11] - 2*armbbbar[25] + 
   armbbbar[32] - armbbbar[16];
   armbbbar[71]=MMt*armbbbar[71];
   armbbbar[69]=80./9.*armbbbar[71] + 80./9.*armbbbar[69] + 1./27.*
   armbbbar[173] + 10*armbbbar[215];
   armbbbar[63]=1./2.*armbbbar[65] + 1./8.*armbbbar[64] + 1./3.*
   armbbbar[69] + armbbbar[63];
   armbbbar[63]=armbbbar[4]*armbbbar[63];
   armbbbar[64]=17./3.*armbbbar[11] - 5./4.*armbbbar[12];
   armbbbar[64]=armbbbar[12]*armbbbar[64];
   armbbbar[65]= - 55*armbbbar[11] - 47./2.*armbbbar[12];
   armbbbar[65]=1./4.*armbbbar[65] + 13*armbbbar[13];
   armbbbar[65]=armbbbar[13]*armbbbar[65];
   armbbbar[64]=1./2.*armbbbar[64] + 1./3.*armbbbar[65];
   armbbbar[64]=armbbbar[297]*armbbbar[64];
   armbbbar[65]= - 35./6.*armbbbar[12] + 1./2.*armbbbar[209] + 13./16.*
   armbbbar[8] + 53./9.*armbbbar[11];
   armbbbar[65]=armbbbar[12]*armbbbar[65];
   armbbbar[69]= - 371./6.*armbbbar[12] + armbbbar[8] - 7./9.*
   armbbbar[11];
   armbbbar[69]=1./2.*armbbbar[69] + 67./3.*armbbbar[13];
   armbbbar[69]=armbbbar[13]*armbbbar[69];
   armbbbar[64]=1./3.*armbbbar[64] + 1./4.*armbbbar[69] + 1./2.*
   armbbbar[115] + armbbbar[65];
   armbbbar[64]=armbbbar[297]*armbbbar[64];
   armbbbar[65]=armbbbar[52]*armbbbar[297]*armbbbar[89];
   armbbbar[64]=armbbbar[64] + 1./2.*armbbbar[65];
   armbbbar[64]=armbbbar[52]*armbbbar[64];
   armbbbar[65]= - 310./9. + armbbbar[229];
   armbbbar[65]= - 2*armbbbar[45] + 1./3.*armbbbar[65] - 10*
   armbbbar[70];
   armbbbar[65]=armbbbar[11]*armbbbar[65];
   armbbbar[69]= - armbbbar[13]*armbbbar[54];
   armbbbar[69]=2*armbbbar[69] - 23./3. + 2*armbbbar[96];
   armbbbar[69]=armbbbar[13]*armbbbar[69];
   armbbbar[71]= - 32*armbbbar[18] + 5./3.*armbbbar[6];
   armbbbar[72]=armbbbar[15]*armbbbar[280];
   armbbbar[76]= - armbbbar[70]*armbbbar[33];
   armbbbar[77]=1489./108. + 10*armbbbar[70];
   armbbbar[77]=armbbbar[16]*armbbbar[77];
   armbbbar[78]=3 + 5./3.*armbbbar[70];
   armbbbar[78]=armbbbar[17]*armbbbar[78];
   armbbbar[79]=13 + 5*armbbbar[70];
   armbbbar[79]=armbbbar[12]*armbbbar[79];
   armbbbar[80]= - 19 + armbbbar[202];
   armbbbar[80]=MMt*armbbbar[80];
   armbbbar[59]=armbbbar[63] + 1./4.*armbbbar[64] + armbbbar[59] + 8./
   27.*armbbbar[80] + 8./9.*armbbbar[69] + 4./3.*armbbbar[79] + 1./3.*
   armbbbar[65] - 56./27.*armbbbar[25] + 4*armbbbar[78] + 1./3.*
   armbbbar[77] + 10./3.*armbbbar[76] + 11./36.*armbbbar[27] + 28./27.*
   armbbbar[32] + 16./27.*armbbbar[72] - 11./36.*armbbbar[35] + 1./27.*
   armbbbar[71] - 6*armbbbar[33];
   armbbbar[59]=armbbbar[4]*armbbbar[59];
   armbbbar[63]= - 1 - 1./2.*armbbbar[70];
   armbbbar[64]=armbbbar[24]*armbbbar[63];
   armbbbar[65]=armbbbar[23]*armbbbar[290];
   armbbbar[69]=armbbbar[17]*armbbbar[63];
   armbbbar[71]=armbbbar[25]*armbbbar[63];
   armbbbar[64]=armbbbar[71] + armbbbar[69] + armbbbar[64] + 
   armbbbar[65];
   armbbbar[64]=armbbbar[56]*armbbbar[64];
   armbbbar[64]=13./4.*armbbbar[290] + armbbbar[64];
   armbbbar[64]=armbbbar[22]*armbbbar[56]*armbbbar[64];
   armbbbar[65]=armbbbar[55]*armbbbar[176];
   armbbbar[69]=3./4. + armbbbar[65];
   armbbbar[69]=armbbbar[55]*armbbbar[69];
   armbbbar[69]=armbbbar[124] + armbbbar[69] + armbbbar[122];
   armbbbar[69]=armbbbar[26]*armbbbar[69];
   armbbbar[63]=armbbbar[22]*armbbbar[73]*armbbbar[63];
   armbbbar[71]=armbbbar[8]*armbbbar[63];
   armbbbar[76]= - armbbbar[26]*armbbbar[93];
   armbbbar[63]=1./8.*armbbbar[76] + 3*armbbbar[63];
   armbbbar[77]=armbbbar[12]*armbbbar[63];
   armbbbar[78]=armbbbar[13]*armbbbar[63];
   armbbbar[63]=MMt*armbbbar[63];
   armbbbar[79]=armbbbar[70]*armbbbar[290];
   armbbbar[79]=3./2. + armbbbar[79];
   armbbbar[73]=armbbbar[22]*armbbbar[73]*armbbbar[79];
   armbbbar[79]=armbbbar[26]*armbbbar[93];
   armbbbar[73]=1./32.*armbbbar[79] + 3*armbbbar[73];
   armbbbar[73]=MMZ*armbbbar[73];
   armbbbar[80]= - armbbbar[93]*armbbbar[54];
   armbbbar[81]=armbbbar[11]*armbbbar[76];
   armbbbar[63]=9*armbbbar[73] + armbbbar[63] + 3./2.*armbbbar[78] + 3./
   2.*armbbbar[77] + 3./16.*armbbbar[81] + 9./2.*armbbbar[71] + 9./2.*
   armbbbar[64] + armbbbar[80] + 3./16.*armbbbar[69];
   armbbbar[63]=MMZ*armbbbar[63];
   armbbbar[64]=armbbbar[76] + 3*armbbbar[220];
   armbbbar[64]=armbbbar[12]*armbbbar[64];
   armbbbar[69]=armbbbar[93]*armbbbar[54];
   armbbbar[71]= - armbbbar[26]*armbbbar[55];
   armbbbar[73]=9./32.*armbbbar[71];
   armbbbar[76]=armbbbar[11]*armbbbar[79];
   armbbbar[64]=9./4.*armbbbar[64] + 9./8.*armbbbar[76] + 27./8.*
   armbbbar[165] + 63./32.*armbbbar[211] + 4*armbbbar[69] + 
   armbbbar[73];
   armbbbar[64]=armbbbar[13]*armbbbar[64];
   armbbbar[77]=armbbbar[55]*armbbbar[129];
   armbbbar[78]= - 5./4. + armbbbar[77];
   armbbbar[78]=armbbbar[55]*armbbbar[78];
   armbbbar[82]=armbbbar[17]*armbbbar[93];
   armbbbar[84]=armbbbar[25]*armbbbar[93];
   armbbbar[78]=armbbbar[84] + armbbbar[78] + armbbbar[82];
   armbbbar[78]=armbbbar[26]*armbbbar[78];
   armbbbar[82]=armbbbar[79] + 3*armbbbar[161];
   armbbbar[84]=armbbbar[12]*armbbbar[82];
   armbbbar[82]=armbbbar[13]*armbbbar[82];
   armbbbar[78]=armbbbar[82] + armbbbar[84] + armbbbar[76] + 
   armbbbar[174] + armbbbar[78] + 3*armbbbar[207];
   armbbbar[78]=MMt*armbbbar[78];
   armbbbar[82]=23./6. + armbbbar[70];
   armbbbar[82]=armbbbar[70]*armbbbar[82];
   armbbbar[65]=9./8.*armbbbar[226] + 9./8.*armbbbar[214] + 9./8.*
   armbbbar[65] - 83./48. + armbbbar[82];
   armbbbar[65]=armbbbar[26]*armbbbar[65];
   armbbbar[82]=armbbbar[71] + 7*armbbbar[211];
   armbbbar[76]=armbbbar[76] + 1./4.*armbbbar[82] + armbbbar[174];
   armbbbar[76]=armbbbar[12]*armbbbar[76];
   armbbbar[82]=1./2. - armbbbar[40];
   armbbbar[82]=armbbbar[155]*armbbbar[82];
   armbbbar[84]=armbbbar[297]*armbbbar[82];
   armbbbar[82]=3./2.*armbbbar[84] + 19./24.*armbbbar[26] + 3*
   armbbbar[82];
   armbbbar[82]=armbbbar[297]*armbbbar[82];
   armbbbar[84]=armbbbar[16]*armbbbar[54];
   armbbbar[67]=armbbbar[67] + armbbbar[84];
   armbbbar[67]=armbbbar[55]*armbbbar[67];
   armbbbar[84]= - armbbbar[54] + armbbbar[67];
   armbbbar[84]=armbbbar[55]*armbbbar[84];
   armbbbar[85]=armbbbar[25]*armbbbar[69];
   armbbbar[84]=armbbbar[84] + 2*armbbbar[85];
   armbbbar[73]=2*armbbbar[69] + armbbbar[73];
   armbbbar[73]=armbbbar[11]*armbbbar[73];
   armbbbar[66]=armbbbar[155]*armbbbar[66];
   armbbbar[63]=3*armbbbar[63] + 1./2.*armbbbar[82] + 1./4.*
   armbbbar[78] + armbbbar[64] + 3*armbbbar[66] + 9./8.*armbbbar[76] + 
   armbbbar[73] + 63./32.*armbbbar[218] + 9./8.*armbbbar[217] + 2*
   armbbbar[84] + 1./4.*armbbbar[65];
   armbbbar[63]=MMZ*armbbbar[63];
   armbbbar[64]= - 881./72. + armbbbar[103];
   armbbbar[64]= - 19./9.*armbbbar[47] + armbbbar[87] + armbbbar[40] + 
   31./6.*armbbbar[45] + 1./3.*armbbbar[64] + armbbbar[44];
   armbbbar[65]=armbbbar[12]*armbbbar[170];
   armbbbar[65]=armbbbar[175] + armbbbar[65];
   armbbbar[65]=armbbbar[3]*armbbbar[65];
   armbbbar[66]=4./9.*armbbbar[9];
   armbbbar[64]=3./4.*armbbbar[65] + armbbbar[66] + armbbbar[190] + 
   armbbbar[74] + 1./2.*armbbbar[64] + armbbbar[193];
   armbbbar[64]=armbbbar[3]*armbbbar[64];
   armbbbar[65]=19./12.*armbbbar[213] + 19./12.*armbbbar[109] + 19./12.
   *armbbbar[181] + armbbbar[123] - 125./48.*armbbbar[54];
   armbbbar[73]= - 3301./18. + 7*armbbbar[45];
   armbbbar[74]=armbbbar[3]*armbbbar[175];
   armbbbar[73]=3./2.*armbbbar[74] + 35./18.*armbbbar[9] + 7*
   armbbbar[39] + 119./18.*armbbbar[47] - 17./6.*armbbbar[43] + 1./12.*
   armbbbar[73] + armbbbar[40];
   armbbbar[73]=armbbbar[3]*armbbbar[73];
   armbbbar[74]= - armbbbar[57] + 5./48.*armbbbar[54];
   armbbbar[76]= - armbbbar[297]*armbbbar[3]*armbbbar[43];
   armbbbar[73]=1./6.*armbbbar[76] + 25./4.*armbbbar[74] + armbbbar[73]
   ;
   armbbbar[73]=armbbbar[297]*armbbbar[73];
   armbbbar[74]= - armbbbar[13]*armbbbar[26];
   armbbbar[76]= - MMt*armbbbar[26];
   armbbbar[64]=1./4.*armbbbar[73] + 19./48.*armbbbar[76] + 19./96.*
   armbbbar[74] + 1./8.*armbbbar[65] + armbbbar[64];
   armbbbar[64]=armbbbar[297]*armbbbar[64];
   armbbbar[65]=armbbbar[55]*armbbbar[54];
   armbbbar[73]= - 1./6. - armbbbar[70];
   armbbbar[74]=armbbbar[26]*armbbbar[73];
   armbbbar[76]=1./8.*armbbbar[74];
   armbbbar[78]=armbbbar[11]*armbbbar[80];
   armbbbar[80]=armbbbar[12]*armbbbar[71];
   armbbbar[69]=armbbbar[13]*armbbbar[69];
   armbbbar[69]=8*armbbbar[69] + 9./8.*armbbbar[80] + 8*armbbbar[78] + 
   8./3.*armbbbar[65] + armbbbar[76];
   armbbbar[69]=armbbbar[13]*armbbbar[69];
   armbbbar[78]=armbbbar[70]*armbbbar[27];
   armbbbar[80]=armbbbar[16]*armbbbar[73];
   armbbbar[82]=armbbbar[17]*armbbbar[73];
   armbbbar[73]=armbbbar[25]*armbbbar[73];
   armbbbar[73]=armbbbar[73] + armbbbar[82] + armbbbar[80] + 1./6.*
   armbbbar[27] + armbbbar[78];
   armbbbar[73]=armbbbar[26]*armbbbar[73];
   armbbbar[78]=365./18. + armbbbar[229];
   armbbbar[80]= - 4*armbbbar[70];
   armbbbar[82]= - 53./3. + armbbbar[80];
   armbbbar[82]=armbbbar[45]*armbbbar[82];
   armbbbar[66]=armbbbar[66] + 16./9.*armbbbar[47] - 6*armbbbar[43] + 
   armbbbar[82] + 1./3.*armbbbar[78] + armbbbar[80];
   armbbbar[66]=armbbbar[3]*armbbbar[66];
   armbbbar[78]=armbbbar[17]*armbbbar[55];
   armbbbar[80]=armbbbar[25]*armbbbar[55];
   armbbbar[77]=5./6.*armbbbar[80] + 5./6.*armbbbar[78] + 5./6.*
   armbbbar[77] - 85./6. - armbbbar[70];
   armbbbar[77]=armbbbar[26]*armbbbar[77];
   armbbbar[78]=armbbbar[26]*armbbbar[55];
   armbbbar[80]=armbbbar[11]*armbbbar[78];
   armbbbar[77]=armbbbar[77] + 5./6.*armbbbar[80];
   armbbbar[81]=5./12.*armbbbar[78] + armbbbar[81];
   armbbbar[82]=armbbbar[12]*armbbbar[81];
   armbbbar[77]=1./2.*armbbbar[77] + armbbbar[82];
   armbbbar[79]=armbbbar[12]*armbbbar[79];
   armbbbar[79]=1./2.*armbbbar[81] + armbbbar[79];
   armbbbar[79]=armbbbar[13]*armbbbar[79];
   armbbbar[77]=1./2.*armbbbar[77] + armbbbar[79];
   armbbbar[77]=MMt*armbbbar[77];
   armbbbar[79]=armbbbar[25]*armbbbar[65];
   armbbbar[67]=4*armbbbar[79] + armbbbar[54] + 2*armbbbar[67];
   armbbbar[76]=4./3.*armbbbar[65] + armbbbar[76];
   armbbbar[76]=armbbbar[11]*armbbbar[76];
   armbbbar[74]=armbbbar[74] + 9*armbbbar[80];
   armbbbar[74]=armbbbar[12]*armbbbar[74];
   armbbbar[63]=armbbbar[63] + armbbbar[64] + armbbbar[77] + 
   armbbbar[69] + armbbbar[66] + 1./8.*armbbbar[74] + armbbbar[76] + 2./
   3.*armbbbar[67] + 1./8.*armbbbar[73];
   armbbbar[63]=MMZ*armbbbar[63];
   armbbbar[64]=armbbbar[11]*armbbbar[71];
   armbbbar[66]=armbbbar[12]*armbbbar[78];
   armbbbar[66]=5./3.*armbbbar[66] + 1./4.*armbbbar[26] + 1./3.*
   armbbbar[64];
   armbbbar[66]=armbbbar[13]*armbbbar[66];
   armbbbar[64]=1./8.*armbbbar[26] + 2./3.*armbbbar[64];
   armbbbar[64]=armbbbar[12]*armbbbar[64];
   armbbbar[67]=1 + armbbbar[47];
   armbbbar[67]=armbbbar[3]*armbbbar[67];
   armbbbar[62]=1./2.*armbbbar[66] + 32./9.*armbbbar[67] + 1./8.*
   armbbbar[62] + armbbbar[64];
   armbbbar[62]=MMt*armbbbar[62];
   armbbbar[64]=10*armbbbar[72] - 13*armbbbar[18] - 101./3.*
   armbbbar[14];
   armbbbar[64]=armbbbar[15]*armbbbar[64];
   armbbbar[64]=armbbbar[64] + 1849./18. - 13*armbbbar[19];
   armbbbar[64]=1./3.*armbbbar[64] + 4*armbbbar[131];
   armbbbar[66]= - 1 - 2*armbbbar[70];
   armbbbar[66]=armbbbar[70]*armbbbar[66];
   armbbbar[64]=1./3.*armbbbar[64] + 10*armbbbar[66];
   armbbbar[66]=armbbbar[15]*armbbbar[14]*armbbbar[51];
   armbbbar[66]= - 13./3.*armbbbar[51] + 2*armbbbar[66];
   armbbbar[66]=1./3.*armbbbar[66] - armbbbar[54];
   armbbbar[67]= - armbbbar[55]*armbbbar[54];
   armbbbar[67]=16*armbbbar[67] + 5./8.*armbbbar[26];
   armbbbar[67]=armbbbar[11]*armbbbar[67];
   armbbbar[66]=8./3.*armbbbar[66] + armbbbar[67];
   armbbbar[67]=pow(armbbbar[51],2);
   armbbbar[65]= - 1./9.*armbbbar[67] + 2*armbbbar[65];
   armbbbar[65]=armbbbar[13]*armbbbar[65];
   armbbbar[65]=8./3.*armbbbar[65] + 32./9.*armbbbar[3] + 1./3.*
   armbbbar[66] + 1./8.*armbbbar[68];
   armbbbar[65]=armbbbar[13]*armbbbar[65];
   armbbbar[66]= - 5 + armbbbar[75];
   armbbbar[66]=armbbbar[12]*armbbbar[66];
   armbbbar[66]=armbbbar[137] + 2*armbbbar[66];
   armbbbar[66]=armbbbar[3]*armbbbar[66];
   armbbbar[67]= - 1 - 1./3.*armbbbar[70];
   armbbbar[67]=armbbbar[45]*armbbbar[67];
   armbbbar[68]= - armbbbar[16]*armbbbar[54];
   armbbbar[69]=armbbbar[12]*armbbbar[109];
   armbbbar[58]=armbbbar[60] + armbbbar[59] + armbbbar[61] + 
   armbbbar[63] + armbbbar[58] + armbbbar[62] + armbbbar[65] + 
   armbbbar[66] + 1./3.*armbbbar[69] + 4./9.*armbbbar[110] + 4./27.*
   armbbbar[9] + 8./9.*armbbbar[284] + 52./27.*armbbbar[47] + 4./9.*
   armbbbar[68] + 1./3.*armbbbar[64] + 4*armbbbar[67];
   armbbbar[58]=armbbbar[21]*armbbbar[58];
   armbbbar[59]=1./2. + armbbbar[9];
   armbbbar[60]=armbbbar[20]*armbbbar[59];
   armbbbar[59]=armbbbar[1]*armbbbar[59];
   armbbbar[61]=3*armbbbar[60] + armbbbar[59];
   armbbbar[61]=armbbbar[12]*armbbbar[61];
   armbbbar[62]= - 1 + armbbbar[9];
   armbbbar[63]=armbbbar[20]*armbbbar[62];
   armbbbar[62]=armbbbar[1]*armbbbar[62];
   armbbbar[62]=3*armbbbar[63] + armbbbar[62];
   armbbbar[62]=MMZ*armbbbar[62];
   armbbbar[63]= - armbbbar[7] + armbbbar[17];
   armbbbar[64]=armbbbar[20]*armbbbar[63];
   armbbbar[63]=armbbbar[1]*armbbbar[63];
   armbbbar[65]=3*armbbbar[20] + armbbbar[1];
   armbbbar[65]=armbbbar[13]*armbbbar[65];
   armbbbar[62]=3./2.*armbbbar[62] + 1./2.*armbbbar[65] + armbbbar[61]
    + 3*armbbbar[64] + armbbbar[63];
   armbbbar[63]=3*armbbbar[9];
   armbbbar[64]=1 + armbbbar[10];
   armbbbar[65]=armbbbar[64] + armbbbar[63];
   armbbbar[66]=armbbbar[20]*armbbbar[65];
   armbbbar[65]=armbbbar[1]*armbbbar[65];
   armbbbar[65]=3*armbbbar[66] + armbbbar[65];
   armbbbar[65]=armbbbar[12]*armbbbar[65];
   armbbbar[66]= - 3*armbbbar[9];
   armbbbar[64]=armbbbar[64] + armbbbar[66];
   armbbbar[67]=armbbbar[20]*armbbbar[64];
   armbbbar[64]=armbbbar[1]*armbbbar[64];
   armbbbar[64]=3*armbbbar[67] + armbbbar[64];
   armbbbar[64]=armbbbar[13]*armbbbar[64];
   armbbbar[67]= - 1 + armbbbar[10];
   armbbbar[68]=armbbbar[20]*armbbbar[67];
   armbbbar[67]=armbbbar[1]*armbbbar[67];
   armbbbar[67]=3*armbbbar[68] + armbbbar[67];
   armbbbar[67]=1./2.*MMZ*armbbbar[67];
   armbbbar[68]=armbbbar[106] + armbbbar[17] + 1./2.*armbbbar[53] - 
   armbbbar[7];
   armbbbar[69]=armbbbar[20]*armbbbar[68];
   armbbbar[68]=armbbbar[1]*armbbbar[68];
   armbbbar[64]=armbbbar[67] + 1./4.*armbbbar[64] + 1./4.*armbbbar[65]
    + 3*armbbbar[69] + armbbbar[68];
   armbbbar[64]=MMZ*armbbbar[64];
   armbbbar[65]=armbbbar[20]*armbbbar[10];
   armbbbar[68]=armbbbar[1]*armbbbar[10];
   armbbbar[65]=3*armbbbar[65] + armbbbar[68];
   armbbbar[68]=armbbbar[12]*armbbbar[65];
   armbbbar[69]= - armbbbar[20]*armbbbar[10];
   armbbbar[70]= - armbbbar[1]*armbbbar[10];
   armbbbar[69]=3*armbbbar[69] + armbbbar[70];
   armbbbar[70]=armbbbar[13]*armbbbar[69];
   armbbbar[71]= - armbbbar[25] + armbbbar[17] + armbbbar[53] - 
   armbbbar[7];
   armbbbar[72]=armbbbar[20]*armbbbar[71];
   armbbbar[71]=armbbbar[1]*armbbbar[71];
   armbbbar[68]=armbbbar[70] + armbbbar[68] + 3*armbbbar[72] + 
   armbbbar[71];
   armbbbar[68]=1./2.*armbbbar[52]*armbbbar[127]*armbbbar[68];
   armbbbar[64]=armbbbar[64] + armbbbar[68];
   armbbbar[64]=armbbbar[52]*armbbbar[64];
   armbbbar[62]=1./2.*armbbbar[62] + armbbbar[64];
   armbbbar[62]=armbbbar[52]*armbbbar[62];
   armbbbar[64]=armbbbar[10] - armbbbar[9];
   armbbbar[70]=armbbbar[20]*armbbbar[64];
   armbbbar[64]=armbbbar[1]*armbbbar[64];
   armbbbar[71]=3*armbbbar[70] + armbbbar[64];
   armbbbar[72]=armbbbar[12]*armbbbar[71];
   armbbbar[73]= - armbbbar[10] + armbbbar[9];
   armbbbar[74]=armbbbar[20]*armbbbar[73];
   armbbbar[73]=armbbbar[1]*armbbbar[73];
   armbbbar[75]=3*armbbbar[74] + armbbbar[73];
   armbbbar[76]=armbbbar[13]*armbbbar[75];
   armbbbar[77]=armbbbar[72] + armbbbar[76];
   armbbbar[77]=armbbbar[52]*MMZ*armbbbar[77];
   armbbbar[78]=MMZ*armbbbar[71];
   armbbbar[77]=armbbbar[77] + armbbbar[72] + armbbbar[78];
   armbbbar[77]=armbbbar[52]*armbbbar[77];
   armbbbar[78]=armbbbar[11]*armbbbar[71];
   armbbbar[72]=1./2.*armbbbar[78] + armbbbar[72];
   armbbbar[72]=armbbbar[3]*armbbbar[72];
   armbbbar[78]=armbbbar[3]*armbbbar[71];
   armbbbar[79]=MMZ*armbbbar[78];
   armbbbar[72]=1./4.*armbbbar[77] + armbbbar[79] + 1./8.*armbbbar[75]
    + armbbbar[72];
   armbbbar[77]=armbbbar[8]*armbbbar[71];
   armbbbar[71]=armbbbar[13]*armbbbar[71];
   armbbbar[71]=armbbbar[77] + armbbbar[71];
   armbbbar[64]=armbbbar[70] + 1./3.*armbbbar[64];
   armbbbar[70]=armbbbar[13]*armbbbar[3]*armbbbar[75];
   armbbbar[64]=1./16.*armbbbar[64] + armbbbar[70];
   armbbbar[64]=MMt*armbbbar[64];
   armbbbar[64]=1./8.*armbbbar[71] + armbbbar[64];
   armbbbar[64]=armbbbar[4]*armbbbar[64];
   armbbbar[64]=1./2.*armbbbar[72] + armbbbar[64];
   armbbbar[64]=armbbbar[91]*armbbbar[64];
   armbbbar[70]=1 - armbbbar[10];
   armbbbar[71]=armbbbar[70] + armbbbar[63];
   armbbbar[72]=armbbbar[20]*armbbbar[71];
   armbbbar[71]=armbbbar[1]*armbbbar[71];
   armbbbar[71]=3*armbbbar[72] + armbbbar[71];
   armbbbar[71]=armbbbar[11]*armbbbar[71];
   armbbbar[72]=armbbbar[17] + 1./2.*armbbbar[16] - armbbbar[7] - 1./2.
   *armbbbar[6];
   armbbbar[77]=armbbbar[20]*armbbbar[72];
   armbbbar[72]=armbbbar[1]*armbbbar[72];
   armbbbar[61]=armbbbar[61] + 1./4.*armbbbar[71] + 3*armbbbar[77] + 
   armbbbar[72];
   armbbbar[61]=armbbbar[3]*armbbbar[61];
   armbbbar[71]=armbbbar[8]*armbbbar[75];
   armbbbar[71]=armbbbar[71] + armbbbar[76];
   armbbbar[72]=armbbbar[74] + 1./3.*armbbbar[73];
   armbbbar[73]=armbbbar[13]*armbbbar[78];
   armbbbar[72]=1./16.*armbbbar[72] + armbbbar[73];
   armbbbar[72]=MMt*armbbbar[72];
   armbbbar[71]=1./8.*armbbbar[71] + armbbbar[72];
   armbbbar[71]=armbbbar[4]*armbbbar[71];
   armbbbar[72]=1./4. + armbbbar[10];
   armbbbar[72]=1./2.*armbbbar[72] - armbbbar[9];
   armbbbar[73]=armbbbar[20]*armbbbar[72];
   armbbbar[72]=armbbbar[1]*armbbbar[72];
   armbbbar[72]=3*armbbbar[73] + armbbbar[72];
   armbbbar[73]= - 11 + 21./2.*armbbbar[9];
   armbbbar[73]=armbbbar[20]*armbbbar[73];
   armbbbar[74]= - 11./3. + armbbbar[104];
   armbbbar[75]=armbbbar[1]*armbbbar[74];
   armbbbar[73]=armbbbar[73] + armbbbar[75];
   armbbbar[73]=MMZ*armbbbar[3]*armbbbar[73];
   armbbbar[61]=armbbbar[64] + armbbbar[71] + 1./2.*armbbbar[62] + 1./2.
   *armbbbar[73] + 1./8.*armbbbar[72] + armbbbar[61];
   armbbbar[61]=armbbbar[91]*armbbbar[61];
   armbbbar[62]= - 3*armbbbar[10];
   armbbbar[64]= - 19./3.*armbbbar[9];
   armbbbar[71]=armbbbar[64] - 7./18. + armbbbar[62];
   armbbbar[71]=armbbbar[20]*armbbbar[71];
   armbbbar[72]=armbbbar[66] + 1./6. - armbbbar[10];
   armbbbar[72]=armbbbar[1]*armbbbar[72];
   armbbbar[71]=armbbbar[71] + armbbbar[72];
   armbbbar[71]=armbbbar[12]*armbbbar[71];
   armbbbar[72]= - 47./6. + 19*armbbbar[9];
   armbbbar[72]=armbbbar[20]*armbbbar[72];
   armbbbar[73]= - 7./6. + armbbbar[63];
   armbbbar[73]=armbbbar[1]*armbbbar[73];
   armbbbar[72]=1./3.*armbbbar[72] + armbbbar[73];
   armbbbar[72]=armbbbar[13]*armbbbar[72];
   armbbbar[73]=armbbbar[20]*armbbbar[70];
   armbbbar[70]=armbbbar[1]*armbbbar[70];
   armbbbar[70]=3*armbbbar[73] + armbbbar[70];
   armbbbar[70]=MMZ*armbbbar[70];
   armbbbar[73]=armbbbar[121] - armbbbar[17] - 1./2.*armbbbar[53] + 
   armbbbar[7];
   armbbbar[76]=armbbbar[20]*armbbbar[73];
   armbbbar[73]=armbbbar[1]*armbbbar[73];
   armbbbar[70]=armbbbar[70] + 1./2.*armbbbar[72] + 1./2.*armbbbar[71]
    + 3*armbbbar[76] + armbbbar[73];
   armbbbar[70]=MMZ*armbbbar[70];
   armbbbar[69]=armbbbar[12]*armbbbar[69];
   armbbbar[65]=armbbbar[13]*armbbbar[65];
   armbbbar[71]=armbbbar[25] - armbbbar[17] - armbbbar[53] + 
   armbbbar[7];
   armbbbar[72]=armbbbar[20]*armbbbar[71];
   armbbbar[71]=armbbbar[1]*armbbbar[71];
   armbbbar[65]=armbbbar[65] + armbbbar[69] + 3*armbbbar[72] + 
   armbbbar[71];
   armbbbar[65]=armbbbar[52]*armbbbar[127]*armbbbar[65];
   armbbbar[65]=armbbbar[70] + armbbbar[65];
   armbbbar[65]=armbbbar[52]*armbbbar[65];
   armbbbar[62]=101./9. + armbbbar[62];
   armbbbar[62]=1./2.*armbbbar[62] + armbbbar[64];
   armbbbar[62]=armbbbar[20]*armbbbar[62];
   armbbbar[64]=13./3. - armbbbar[10];
   armbbbar[64]=1./2.*armbbbar[64] + armbbbar[66];
   armbbbar[64]=armbbbar[1]*armbbbar[64];
   armbbbar[62]=armbbbar[62] + armbbbar[64];
   armbbbar[62]=MMZ*armbbbar[62];
   armbbbar[64]=1./3. - armbbbar[9];
   armbbbar[66]=armbbbar[20]*armbbbar[64];
   armbbbar[64]=armbbbar[1]*armbbbar[64];
   armbbbar[64]=5./3.*armbbbar[66] + armbbbar[64];
   armbbbar[66]=armbbbar[12]*armbbbar[64];
   armbbbar[62]=armbbbar[65] + armbbbar[66] + 1./2.*armbbbar[62];
   armbbbar[62]=armbbbar[52]*armbbbar[62];
   armbbbar[65]= - 11./12.*armbbbar[9] + 5./9. - 3./4.*armbbbar[10];
   armbbbar[65]=armbbbar[20]*armbbbar[65];
   armbbbar[69]= - 1./4.*armbbbar[10];
   armbbbar[70]=1./3. + armbbbar[69];
   armbbbar[71]=armbbbar[70] - 3./4.*armbbbar[9];
   armbbbar[71]=armbbbar[1]*armbbbar[71];
   armbbbar[65]=armbbbar[65] + armbbbar[71];
   armbbbar[71]=armbbbar[8]*armbbbar[65];
   armbbbar[72]= - 1./3.*armbbbar[9] - 1 + 1./3.*armbbbar[10];
   armbbbar[73]=armbbbar[20]*armbbbar[72];
   armbbbar[72]=armbbbar[1]*armbbbar[72];
   armbbbar[72]=armbbbar[73] + 1./3.*armbbbar[72];
   armbbbar[72]=armbbbar[11]*armbbbar[72];
   armbbbar[73]=armbbbar[13]*armbbbar[65];
   armbbbar[72]=armbbbar[73] + armbbbar[71] + armbbbar[72];
   armbbbar[70]=1./3.*armbbbar[70] + armbbbar[169];
   armbbbar[70]=armbbbar[1]*armbbbar[70];
   armbbbar[69]= - 11./36.*armbbbar[9] + 5./27. + armbbbar[69];
   armbbbar[69]=armbbbar[20]*armbbbar[69];
   armbbbar[69]=armbbbar[69] + armbbbar[70];
   armbbbar[70]=11./3.*armbbbar[9] - 20./9. + 3*armbbbar[10];
   armbbbar[70]=armbbbar[20]*armbbbar[70];
   armbbbar[63]=armbbbar[63] - 4./3. + armbbbar[10];
   armbbbar[63]=armbbbar[1]*armbbbar[63];
   armbbbar[63]=armbbbar[70] + armbbbar[63];
   armbbbar[63]=armbbbar[13]*armbbbar[3]*armbbbar[63];
   armbbbar[63]=1./4.*armbbbar[69] + armbbbar[63];
   armbbbar[63]=MMt*armbbbar[63];
   armbbbar[69]=1./2.*armbbbar[72] + armbbbar[63];
   armbbbar[69]=armbbbar[4]*armbbbar[69];
   armbbbar[65]=armbbbar[11]*armbbbar[65];
   armbbbar[66]=armbbbar[65] + 2*armbbbar[66];
   armbbbar[66]=armbbbar[3]*armbbbar[66];
   armbbbar[70]=31./4.*armbbbar[10];
   armbbbar[72]= - 43./4.*armbbbar[9] + 13 + armbbbar[70];
   armbbbar[72]=armbbbar[20]*armbbbar[72];
   armbbbar[76]=5./36.*armbbbar[9] + 1 + 31./36.*armbbbar[10];
   armbbbar[76]=armbbbar[1]*armbbbar[76];
   armbbbar[72]=1./3.*armbbbar[72] + armbbbar[76];
   armbbbar[76]= - 2*armbbbar[10];
   armbbbar[77]= - 29./6.*armbbbar[9] + 43./9. + armbbbar[76];
   armbbbar[77]=armbbbar[20]*armbbbar[77];
   armbbbar[76]=17./3. + armbbbar[76];
   armbbbar[76]=1./3.*armbbbar[76] - 5./2.*armbbbar[9];
   armbbbar[76]=armbbbar[1]*armbbbar[76];
   armbbbar[76]=armbbbar[77] + armbbbar[76];
   armbbbar[76]=MMZ*armbbbar[3]*armbbbar[76];
   armbbbar[61]=armbbbar[61] + armbbbar[69] + 1./2.*armbbbar[62] + 
   armbbbar[76] + 1./12.*armbbbar[72] + armbbbar[66];
   armbbbar[61]=armbbbar[91]*armbbbar[61];
   armbbbar[62]= - 1./3. + armbbbar[9];
   armbbbar[66]=armbbbar[20]*armbbbar[62];
   armbbbar[62]=armbbbar[1]*armbbbar[62];
   armbbbar[62]=5./3.*armbbbar[66] + armbbbar[62];
   armbbbar[66]=armbbbar[12]*armbbbar[62];
   armbbbar[64]=armbbbar[13]*armbbbar[64];
   armbbbar[64]=armbbbar[67] + armbbbar[66] + armbbbar[64];
   armbbbar[64]=MMZ*armbbbar[64];
   armbbbar[64]=armbbbar[64] + armbbbar[68];
   armbbbar[64]=armbbbar[52]*armbbbar[64];
   armbbbar[66]=MMZ*armbbbar[62];
   armbbbar[64]=armbbbar[66] + armbbbar[64];
   armbbbar[64]=armbbbar[52]*armbbbar[64];
   armbbbar[66]=11./9.*armbbbar[9] - 53./27. + armbbbar[10];
   armbbbar[66]=armbbbar[20]*armbbbar[66];
   armbbbar[67]= - 13./3. + armbbbar[10];
   armbbbar[67]=1./3.*armbbbar[67] + armbbbar[9];
   armbbbar[67]=armbbbar[1]*armbbbar[67];
   armbbbar[66]=armbbbar[66] + armbbbar[67];
   armbbbar[66]=armbbbar[11]*armbbbar[66];
   armbbbar[66]=armbbbar[73] + armbbbar[71] + 1./3.*armbbbar[66];
   armbbbar[67]= - 1 - armbbbar[9];
   armbbbar[68]=armbbbar[20]*armbbbar[67];
   armbbbar[67]=armbbbar[1]*armbbbar[67];
   armbbbar[67]=11./9.*armbbbar[68] + armbbbar[67];
   armbbbar[67]=armbbbar[11]*armbbbar[67];
   armbbbar[68]=armbbbar[6] - armbbbar[16];
   armbbbar[69]=armbbbar[20]*armbbbar[68];
   armbbbar[68]=armbbbar[1]*armbbbar[68];
   armbbbar[67]=armbbbar[67] + 11./9.*armbbbar[69] + armbbbar[68];
   armbbbar[67]=armbbbar[297]*armbbbar[67];
   armbbbar[63]=1./6.*armbbbar[67] + 1./2.*armbbbar[66] + armbbbar[63];
   armbbbar[63]=armbbbar[297]*armbbbar[63];
   armbbbar[66]=5./3.*armbbbar[20] + armbbbar[1];
   armbbbar[66]=armbbbar[11]*armbbbar[66];
   armbbbar[66]=2./3.*armbbbar[66] + 5./3.*armbbbar[69] + armbbbar[68];
   armbbbar[63]=4./9.*armbbbar[66] + armbbbar[63];
   armbbbar[63]=armbbbar[4]*armbbbar[63];
   armbbbar[59]=11./3.*armbbbar[60] + 3*armbbbar[59];
   armbbbar[59]=armbbbar[11]*armbbbar[59];
   armbbbar[60]= - armbbbar[6] + armbbbar[16];
   armbbbar[66]=armbbbar[20]*armbbbar[60];
   armbbbar[60]=armbbbar[1]*armbbbar[60];
   armbbbar[59]=armbbbar[59] + 11./3.*armbbbar[66] + 3*armbbbar[60];
   armbbbar[59]=armbbbar[3]*armbbbar[59];
   armbbbar[60]=631./12. - 55*armbbbar[9];
   armbbbar[66]=armbbbar[20]*armbbbar[60];
   armbbbar[60]=armbbbar[1]*armbbbar[60];
   armbbbar[60]=11./9.*armbbbar[66] + armbbbar[60];
   armbbbar[59]=1./72.*armbbbar[60] + armbbbar[59];
   armbbbar[59]=armbbbar[297]*armbbbar[59];
   armbbbar[60]=77./36.*armbbbar[9] + 43./27. + armbbbar[70];
   armbbbar[60]=armbbbar[20]*armbbbar[60];
   armbbbar[66]=23./3. + armbbbar[70];
   armbbbar[66]=1./3.*armbbbar[66] + 7./4.*armbbbar[9];
   armbbbar[66]=armbbbar[1]*armbbbar[66];
   armbbbar[60]=armbbbar[60] + armbbbar[66];
   armbbbar[65]=armbbbar[3]*armbbbar[65];
   armbbbar[59]=1./2.*armbbbar[59] + 1./36.*armbbbar[60] + armbbbar[65]
   ;
   armbbbar[59]=armbbbar[297]*armbbbar[59];
   armbbbar[60]= - 1./2.*armbbbar[10];
   armbbbar[65]=2./3. + armbbbar[60];
   armbbbar[65]=1./3.*armbbbar[65] + armbbbar[83];
   armbbbar[65]=armbbbar[1]*armbbbar[65];
   armbbbar[60]= - 11./18.*armbbbar[9] + 10./27. + armbbbar[60];
   armbbbar[60]=armbbbar[20]*armbbbar[60];
   armbbbar[60]=armbbbar[60] + armbbbar[65];
   armbbbar[60]=armbbbar[3]*armbbbar[60];
   armbbbar[65]=armbbbar[20]*armbbbar[74];
   armbbbar[65]=11./9.*armbbbar[65] + armbbbar[75];
   armbbbar[65]=armbbbar[297]*armbbbar[3]*armbbbar[65];
   armbbbar[60]=armbbbar[60] + 1./2.*armbbbar[65];
   armbbbar[60]=armbbbar[297]*armbbbar[60];
   armbbbar[62]=armbbbar[3]*armbbbar[62];
   armbbbar[60]=4./3.*armbbbar[62] + armbbbar[60];
   armbbbar[60]=MMZ*armbbbar[60];
   armbbbar[62]=armbbbar[201] - 13./3.*armbbbar[14];
   armbbbar[62]=armbbbar[15]*armbbbar[62];
   armbbbar[62]=4*armbbbar[9] + 2*armbbbar[62] + 19./3. + 4*
   armbbbar[19];
   armbbbar[65]=armbbbar[20]*armbbbar[62];
   armbbbar[62]=armbbbar[1]*armbbbar[62];
   armbbbar[62]=5./3.*armbbbar[65] + armbbbar[62];

      mbbbarret = armbbbar[58] + armbbbar[59] + armbbbar[60] + 
      armbbbar[61] + 1./9.*armbbbar[62] + armbbbar[63] + 1./2.*
      armbbbar[64];
      return mbbbarret;
}
