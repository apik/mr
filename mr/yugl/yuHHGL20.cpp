#include <HH.hpp>
namespace mr
{
  long double HH<OS>::ygl20(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<long double> aryuHHGL[59], yuHHGLret;

    aryuHHGL[1]=double(boson);
    aryuHHGL[2]=pow(SW,-1);
    aryuHHGL[3]=pow(MMH,-1);
    aryuHHGL[4]=pow(MMW,-1);
    aryuHHGL[5]=Tsil::I2(MMH,MMH,MMH,mu2);
    aryuHHGL[6]=Tsil::I2(MMH,MMt,MMt,mu2);
    aryuHHGL[7]=Tsil::I2(0,MMH,MMt,mu2);
    aryuHHGL[8]=Tsil::I2(0,0,MMH,mu2);
    aryuHHGL[9]=Tsil::I2(0,0,MMt,mu2);
    aryuHHGL[10]=Tsil::B(MMH,MMH,MMH,mu2);
    aryuHHGL[11]=Tsil::B(MMt,MMt,MMH,mu2);
    aryuHHGL[12]=Tsil::A(MMH,mu2);
    aryuHHGL[13]=Tsil::A(MMt,mu2);
    aryuHHGL[14]=std::real(Tsil::B(0,0,MMH,mu2));
    aryuHHGL[15]=Tsil::B(MMH,MMt,MMt,mu2);
    aryuHHGL[16]=std::real(Tsil::B(0,0,MMt,mu2));
    aryuHHGL[17]=Tsil::B(0,MMt,MMH,mu2);
    aryuHHGL[18]=Tsil::B(0,0,MMH,mu2);
    aryuHHGL[19]=Tsil::Beps(MMH,MMH,MMH,mu2);
    aryuHHGL[20]=Tsil::Beps(MMt,MMt,MMH,mu2);
    aryuHHGL[21]=pow(MMt,-1);
    aryuHHGL[22]=Tsil::Aeps(MMH,mu2);
    aryuHHGL[23]=Tsil::Aeps(MMt,mu2);
    aryuHHGL[24]=prottttt0->Suxv(0);
    aryuHHGL[25]=protHHHHH->M(0);
    aryuHHGL[26]=protHtHtt->M(0);
    aryuHHGL[27]=protttttH->M(0);
    aryuHHGL[28]=protHHHHH->Uzxyv(0);
    aryuHHGL[29]=protHtHtt->Uzxyv(0);
    aryuHHGL[30]=protHtHtt->Uuyxv(0);
    aryuHHGL[31]=protHtHtt->Svyz(0);
    aryuHHGL[32]=prot0H0H0->M(0);
    aryuHHGL[33]=prot0t0tt->M(0);
    aryuHHGL[34]=prot0t0t0->M(0);
    aryuHHGL[35]=prot0000H->M(0);
    aryuHHGL[36]=prot0H0H0->Uuyxv(0);
    aryuHHGL[37]=prot0t0t0->Uuyxv(0);
    aryuHHGL[38]=prot0H0H0->Tyzv(0);
    aryuHHGL[39]=prot0t0t0->Tyzv(0);
    aryuHHGL[40]=1/(4*MMt - MMH);
    aryuHHGL[41]=aryuHHGL[15]*aryuHHGL[3];
    aryuHHGL[42]=aryuHHGL[16]*aryuHHGL[3];
    aryuHHGL[43]=aryuHHGL[41] + 1./4.*aryuHHGL[42];
    aryuHHGL[44]=aryuHHGL[11]*aryuHHGL[3];
    aryuHHGL[45]=pow(aryuHHGL[3],2);
    aryuHHGL[46]=aryuHHGL[45]*aryuHHGL[12];
    aryuHHGL[47]=2*aryuHHGL[45];
    aryuHHGL[47]= - aryuHHGL[13]*aryuHHGL[47];
    aryuHHGL[47]=aryuHHGL[47] + 13./4.*aryuHHGL[44] + 119./16.*
      aryuHHGL[3] - 4*aryuHHGL[46] - aryuHHGL[43];
    aryuHHGL[47]=aryuHHGL[13]*aryuHHGL[47];
    aryuHHGL[48]=aryuHHGL[29] - aryuHHGL[19];
    aryuHHGL[49]= - 9./2.*aryuHHGL[23] + aryuHHGL[6] - 5./4.*
      aryuHHGL[22];
    aryuHHGL[49]=aryuHHGL[3]*aryuHHGL[49];
    aryuHHGL[50]= - 5*aryuHHGL[3] - 21./2.*aryuHHGL[46];
    aryuHHGL[50]=aryuHHGL[12]*aryuHHGL[50];
    aryuHHGL[49]=5./4.*aryuHHGL[20] - 1 + 1./4.*aryuHHGL[50] + 
      aryuHHGL[49] + 9./2.*aryuHHGL[48];
    aryuHHGL[50]=1./2.*aryuHHGL[14];
    aryuHHGL[51]=1./8.*aryuHHGL[16];
    aryuHHGL[52]=aryuHHGL[12]*aryuHHGL[3];
    aryuHHGL[53]= - 1./16.*aryuHHGL[11] + aryuHHGL[51] - aryuHHGL[50] + 
      aryuHHGL[52] + 7./4.*aryuHHGL[15];
    aryuHHGL[53]=aryuHHGL[11]*aryuHHGL[53];
    aryuHHGL[54]=1./4.*MMH;
    aryuHHGL[55]=1./2.*aryuHHGL[27] + 9*aryuHHGL[26] - aryuHHGL[33];
    aryuHHGL[55]=aryuHHGL[55]*aryuHHGL[54];
    aryuHHGL[56]=9./2.*aryuHHGL[10];
    aryuHHGL[57]= - 1./2. - aryuHHGL[11];
    aryuHHGL[57]=aryuHHGL[57]*aryuHHGL[56];
    aryuHHGL[58]= - 1./8.*aryuHHGL[24] - 17./8.*aryuHHGL[31] - 3./8.*
      aryuHHGL[7];
    aryuHHGL[58]=aryuHHGL[3]*aryuHHGL[58];
    aryuHHGL[47]=aryuHHGL[57] + aryuHHGL[55] + 9./16.*aryuHHGL[17] + 
      aryuHHGL[47] + aryuHHGL[53] - 3./16.*aryuHHGL[16] + 1./2.*
      aryuHHGL[49] - aryuHHGL[15] + aryuHHGL[58];
    aryuHHGL[49]=pow(aryuHHGL[2],4);
    aryuHHGL[47]=aryuHHGL[49]*aryuHHGL[47];
    aryuHHGL[42]=aryuHHGL[42] - aryuHHGL[3];
    aryuHHGL[41]= - 7./8.*aryuHHGL[44] - 5*aryuHHGL[41] - 2*aryuHHGL[46]
      - 5./4.*aryuHHGL[42];
    aryuHHGL[41]=aryuHHGL[11]*aryuHHGL[41];
    aryuHHGL[42]=aryuHHGL[6] - aryuHHGL[22];
    aryuHHGL[44]=aryuHHGL[23] - aryuHHGL[12];
    aryuHHGL[42]=9./2.*aryuHHGL[13] + 4*aryuHHGL[44] - 2*aryuHHGL[42];
    aryuHHGL[42]=aryuHHGL[45]*aryuHHGL[42];
    aryuHHGL[44]=5./16.*aryuHHGL[17] - 3./8.*aryuHHGL[39] - 13./4.*
      aryuHHGL[20] + 19./16.;
    aryuHHGL[44]=aryuHHGL[3]*aryuHHGL[44];
    aryuHHGL[41]=aryuHHGL[41] + aryuHHGL[27] - 1./2.*aryuHHGL[34] - 6*
      aryuHHGL[26] + aryuHHGL[42] + aryuHHGL[44] + aryuHHGL[43];
    aryuHHGL[41]=aryuHHGL[41]*aryuHHGL[49];
    aryuHHGL[42]=1./2.*aryuHHGL[39];
    aryuHHGL[43]= - 1./2.*aryuHHGL[17] - aryuHHGL[42] - 9./2.;
    aryuHHGL[43]=aryuHHGL[45]*aryuHHGL[43];
    aryuHHGL[44]=aryuHHGL[27]*aryuHHGL[3];
    aryuHHGL[43]= - 2*aryuHHGL[44] + aryuHHGL[43];
    aryuHHGL[44]=aryuHHGL[49]*MMt;
    aryuHHGL[43]=aryuHHGL[43]*aryuHHGL[44];
    aryuHHGL[41]=aryuHHGL[41] + aryuHHGL[43];
    aryuHHGL[41]=MMt*aryuHHGL[41];
    aryuHHGL[41]=aryuHHGL[41] + aryuHHGL[47];
    aryuHHGL[41]=MMt*aryuHHGL[41];
    aryuHHGL[43]=aryuHHGL[13]*aryuHHGL[3];
    aryuHHGL[43]= - 11./4.*aryuHHGL[43] - 3./8.*aryuHHGL[11] - 
      aryuHHGL[51] - aryuHHGL[50] - 1./4.*aryuHHGL[15] + 57./16. - 2*
      aryuHHGL[52];
    aryuHHGL[43]=aryuHHGL[13]*aryuHHGL[43];
    aryuHHGL[46]=aryuHHGL[15] + aryuHHGL[18];
    aryuHHGL[46]=1./4.*aryuHHGL[20] + 17./16. - aryuHHGL[48] + 1./8.*
      aryuHHGL[46];
    aryuHHGL[47]=1./2.*aryuHHGL[11];
    aryuHHGL[48]= - 1./8.*aryuHHGL[11] + 3./2.*aryuHHGL[14] - 1 - 1./2.*
      aryuHHGL[15];
    aryuHHGL[48]=aryuHHGL[48]*aryuHHGL[47];
    aryuHHGL[50]=1./4.*aryuHHGL[14];
    aryuHHGL[42]= - 3./4.*aryuHHGL[17] + aryuHHGL[42] + aryuHHGL[48] + 3
      *aryuHHGL[46] - aryuHHGL[50];
    aryuHHGL[46]=1./2.*MMH;
    aryuHHGL[42]=aryuHHGL[42]*aryuHHGL[46];
    aryuHHGL[48]= - 1 + 15./4.*aryuHHGL[52];
    aryuHHGL[48]=aryuHHGL[12]*aryuHHGL[48];
    aryuHHGL[48]= - 3./2.*aryuHHGL[6] + 1./2.*aryuHHGL[48] + 3*
      aryuHHGL[22];
    aryuHHGL[51]=1./4.*aryuHHGL[11];
    aryuHHGL[52]= - aryuHHGL[12]*aryuHHGL[51];
    aryuHHGL[48]= - aryuHHGL[24] + aryuHHGL[52] + 1./2.*aryuHHGL[48] + 5
      *aryuHHGL[23];
    aryuHHGL[52]=3./8. + aryuHHGL[11];
    aryuHHGL[52]=MMH*aryuHHGL[52];
    aryuHHGL[52]= - 3./2.*aryuHHGL[13] + aryuHHGL[52];
    aryuHHGL[52]=aryuHHGL[10]*aryuHHGL[52];
    aryuHHGL[42]=3./2.*aryuHHGL[52] + aryuHHGL[42] + 9./16.*aryuHHGL[7]
      + aryuHHGL[43] + 1./2.*aryuHHGL[48];
    aryuHHGL[42]=aryuHHGL[49]*aryuHHGL[42];
    aryuHHGL[43]=aryuHHGL[49]*aryuHHGL[3];
    aryuHHGL[48]=pow(MMt,2);
    aryuHHGL[52]=aryuHHGL[48]*aryuHHGL[43];
    aryuHHGL[53]=aryuHHGL[49]*MMH;
    aryuHHGL[55]=3./8.*aryuHHGL[53];
    aryuHHGL[52]= - aryuHHGL[55] + 2*aryuHHGL[52];
    aryuHHGL[52]=aryuHHGL[30]*aryuHHGL[52];
    aryuHHGL[41]=aryuHHGL[52] + aryuHHGL[41] + aryuHHGL[42];
    aryuHHGL[41]=MMt*aryuHHGL[41];
    aryuHHGL[42]=aryuHHGL[47] - 1;
    aryuHHGL[42]=aryuHHGL[42]*aryuHHGL[51];
    aryuHHGL[51]= - 11 + aryuHHGL[14];
    aryuHHGL[50]=aryuHHGL[51]*aryuHHGL[50];
    aryuHHGL[51]=aryuHHGL[42] + 1./8.;
    aryuHHGL[51]=aryuHHGL[40]*aryuHHGL[51];
    aryuHHGL[52]=9./2.*aryuHHGL[25] + aryuHHGL[32];
    aryuHHGL[51]=3*aryuHHGL[52] + aryuHHGL[51];
    aryuHHGL[51]=MMH*aryuHHGL[51];
    aryuHHGL[52]= - 3./2.*aryuHHGL[28] + 1./2.*aryuHHGL[18] - 31./8. + 3
      *aryuHHGL[19];
    aryuHHGL[57]= - 1 + aryuHHGL[11];
    aryuHHGL[57]=aryuHHGL[13]*aryuHHGL[40]*aryuHHGL[57];
    aryuHHGL[42]=aryuHHGL[51] + aryuHHGL[57] + aryuHHGL[42] + 3*
      aryuHHGL[52] + aryuHHGL[50];
    aryuHHGL[42]=aryuHHGL[42]*aryuHHGL[54];
    aryuHHGL[50]=aryuHHGL[6] - aryuHHGL[12];
    aryuHHGL[51]=aryuHHGL[8] + aryuHHGL[22];
    aryuHHGL[50]=1./2.*aryuHHGL[51] + 3*aryuHHGL[50];
    aryuHHGL[51]=3./4.*aryuHHGL[21] + aryuHHGL[40];
    aryuHHGL[51]=aryuHHGL[13]*aryuHHGL[51];
    aryuHHGL[47]=aryuHHGL[51] - aryuHHGL[47] - 7./2. + aryuHHGL[15];
    aryuHHGL[47]=aryuHHGL[13]*aryuHHGL[47];
    aryuHHGL[51]=aryuHHGL[14]*aryuHHGL[12];
    aryuHHGL[42]=aryuHHGL[42] - 3./4.*aryuHHGL[7] + 1./2.*aryuHHGL[47]
      + 7./8.*aryuHHGL[51] + 1./4.*aryuHHGL[50] + aryuHHGL[23];
    aryuHHGL[42]=MMH*aryuHHGL[42];
    aryuHHGL[47]=pow(aryuHHGL[12],2);
    aryuHHGL[50]=5*aryuHHGL[12] + 9*aryuHHGL[13];
    aryuHHGL[50]=aryuHHGL[13]*aryuHHGL[50];
    aryuHHGL[47]=33./4.*aryuHHGL[47] + aryuHHGL[50];
    aryuHHGL[50]=pow(MMH,2);
    aryuHHGL[51]=aryuHHGL[36]*aryuHHGL[50];
    aryuHHGL[42]= - 9./8.*aryuHHGL[51] + 1./4.*aryuHHGL[47] + 
      aryuHHGL[42];
    aryuHHGL[42]=aryuHHGL[42]*aryuHHGL[49];
    aryuHHGL[47]=3./4. + aryuHHGL[14];
    aryuHHGL[46]=aryuHHGL[47]*aryuHHGL[46];
    aryuHHGL[46]=aryuHHGL[46] + 3./4.*aryuHHGL[12] - aryuHHGL[13];
    aryuHHGL[46]=aryuHHGL[46]*aryuHHGL[53];
    aryuHHGL[47]=aryuHHGL[50]*aryuHHGL[49];
    aryuHHGL[51]=aryuHHGL[10]*aryuHHGL[47];
    aryuHHGL[46]=aryuHHGL[46] + 7./8.*aryuHHGL[51];
    aryuHHGL[46]=aryuHHGL[46]*aryuHHGL[56];
    aryuHHGL[51]=aryuHHGL[5]*aryuHHGL[55];
    aryuHHGL[42]=aryuHHGL[51] + aryuHHGL[42] + aryuHHGL[46];
    aryuHHGL[45]= - aryuHHGL[45]*aryuHHGL[44];
    aryuHHGL[43]=7./4.*aryuHHGL[43] + aryuHHGL[45];
    aryuHHGL[43]=aryuHHGL[9]*aryuHHGL[43];
    aryuHHGL[44]=aryuHHGL[44]*aryuHHGL[3];
    aryuHHGL[45]= - 1./2.*aryuHHGL[49] + aryuHHGL[44];
    aryuHHGL[45]=aryuHHGL[37]*aryuHHGL[45];
    aryuHHGL[43]=5./4.*aryuHHGL[45] + 1./2.*aryuHHGL[43];
    aryuHHGL[43]=aryuHHGL[48]*aryuHHGL[43];
    aryuHHGL[45]=81./32.*S2 + 11./32.*aryuHHGL[38];
    aryuHHGL[45]=aryuHHGL[50]*aryuHHGL[45];
    aryuHHGL[46]=aryuHHGL[35]*pow(MMH,3);
    aryuHHGL[45]=1./32.*aryuHHGL[46] + aryuHHGL[45];
    aryuHHGL[45]=aryuHHGL[49]*aryuHHGL[45];
    aryuHHGL[41]=aryuHHGL[45] + 1./4.*aryuHHGL[42] + aryuHHGL[43] + 
      aryuHHGL[41];
    aryuHHGL[42]= - 35./2.*aryuHHGL[49] - 29*aryuHHGL[44];
    aryuHHGL[42]=MMt*aryuHHGL[42];
    aryuHHGL[42]=9*aryuHHGL[53] + aryuHHGL[42];
    aryuHHGL[42]=MMt*aryuHHGL[42];
    aryuHHGL[42]=37./6.*aryuHHGL[47] + aryuHHGL[42];
    aryuHHGL[42]=aryuHHGL[42]*pow(Pi,2);
    aryuHHGL[41]=3*aryuHHGL[41] + 1./32.*aryuHHGL[42];

    yuHHGLret = aryuHHGL[41]*pow(aryuHHGL[4],2)*aryuHHGL[1];
    return yuHHGLret.real();
  }
} // namespace mr
