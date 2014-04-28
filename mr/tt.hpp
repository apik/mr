#ifndef __TOP_HPP__
#define __TOP_HPP__

#include "tsil.hpp"
#include "mr.hpp"

class tt
{

  long double MMt, MMH, MMW, MMZ, mmu;
  // alpha*alphaS
  Tsil* protWt000;
  Tsil* prot0ttHt;
  Tsil* prot0ttZt;
  Tsil* prot0tt0t;
  TsilSTU* prottH0H;
  TsilSTU* prottZ0Z;
  // alpha^2
  Tsil* protHHttH;
  Tsil* protHZttZ;
  Tsil* protHWt0W;
  Tsil* protHttHt;
  Tsil* protHttZt;
  Tsil* protZZttH;
  Tsil* protZWt0W;
  Tsil* protZttZt;
  Tsil* protZ0tW0;
  Tsil* protWW00Z;
  Tsil* protW00tW;
  Tsil* prot00WW0;
  TsilSTU* prot0W00;
  TsilST* prot000;

  TsilST* protos[20];
public:
  tt()
  {
  }

  tt(long double MMt_,long double MMH_,long double MMW_,long double MMZ_,long double mmu_);// :
    // MMt(MMt_), MMH(MMH_), MMW(MMW_), MMZ(MMZ_), mmu(mmu_);

// #include "tt_x.hpp"
  
  std::complex<long double> Mfin1(TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParameters (&result, x, y, z, u, v, mmu); 
    TSIL_Evaluate( &result, MMt); 
    return TSIL_GetFunction(&result, "M");
  }

  std::complex<long double> Vfin1(TSIL_REAL x, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParametersSTU (&result, x, z, u, v, mmu); 
    TSIL_Evaluate( &result, MMt); 
    return TSIL_GetBoldFunction(&result, "Vxzuv",0);
  }

  std::complex<long double> Ufin1(TSIL_REAL x, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParametersSTU (&result, x, z, u, v, mmu); 
    TSIL_Evaluate( &result, MMt); 
    return TSIL_GetBoldFunction(&result, "Uxzuv",0);
  }

  std::complex<long double> Tfin1(TSIL_REAL x, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParametersST (&result, x, u, v, mmu); 
    TSIL_Evaluate( &result, MMt); 
    return TSIL_GetBoldFunction(&result, "Txuv",0);
  }

  std::complex<long double> Sfin1(TSIL_REAL x, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParametersST (&result, x, u, v, mmu); 
    TSIL_Evaluate( &result, MMt); 
    return TSIL_GetBoldFunction(&result, "Sxuv",0);
  }

  
  void test()
  {
    std::complex<long double> diffMfin[23];
    std::complex<long double> diffVfin[2];
    std::complex<long double> diffUfin[17];
    std::complex<long double> diffTfin[10];
    std::complex<long double> diffSfin[8];
    std::complex<long double> diffIfin[9];

    diffMfin[0] =
      - protHHttH->M(0)
      + Mfin1(MMH,MMH,MMt,MMt,MMH)
      ;

    diffMfin[1] =
      - protHZttZ->M(0)
      + Mfin1(MMH,MMZ,MMt,MMt,MMZ)
      ;

    diffMfin[2] =
      - protHWt0W->M(0)
      + Mfin1(MMH,MMW,MMt,0,MMW)
      ;

    diffMfin[3] =
      - protHttHt->M(0)
      + Mfin1(MMH,MMt,MMt,MMH,MMt)
      ;

    diffMfin[4] =
      - protHttZt->M(0)
      + Mfin1(MMH,MMt,MMt,MMZ,MMt)
      ;

    diffMfin[5] =
      - protHZttZ->M(0)
      + Mfin1(MMZ,MMH,MMt,MMt,MMZ)
      ;

    diffMfin[6] =
      - protZZttH->M(0)
      + Mfin1(MMZ,MMZ,MMt,MMt,MMH)
      ;

    diffMfin[7] =
      - protZWt0W->M(0)
      + Mfin1(MMZ,MMW,MMt,0,MMW)
      ;

    diffMfin[8] =
      - protHttZt->M(0)
      + Mfin1(MMZ,MMt,MMt,MMH,MMt)
      ;

    diffMfin[9] =
      - protZttZt->M(0)
      + Mfin1(MMZ,MMt,MMt,MMZ,MMt)
      ;

    diffMfin[10] =
      - protZ0tW0->M(0)
      + Mfin1(MMZ,0,MMt,MMW,0)
      ;

    diffMfin[11] =
      - protHWt0W->M(0)
      + Mfin1(MMW,MMH,0,MMt,MMW)
      ;

    diffMfin[12] =
      - protZWt0W->M(0)
      + Mfin1(MMW,MMZ,0,MMt,MMW)
      ;

    diffMfin[13] =
      - protWW00Z->M(0)
      + Mfin1(MMW,MMW,0,0,MMZ)
      ;

    diffMfin[14] =
      - protZ0tW0->M(0)
      + Mfin1(MMW,MMt,0,MMZ,0)
      ;

    diffMfin[15] =
      - protWt000->M(0)
      + Mfin1(MMW,MMt,0,0,0)
      ;

    diffMfin[16] =
      - protW00tW->M(0)
      + Mfin1(MMW,0,0,MMt,MMW)
      ;

    diffMfin[17] =
      - protW00tW->M(0)
      + Mfin1(0,MMW,MMt,0,MMW)
      ;

    diffMfin[18] =
      - prot0ttHt->M(0)
      + Mfin1(0,MMt,MMt,MMH,MMt)
      ;

    diffMfin[19] =
      - prot0ttZt->M(0)
      + Mfin1(0,MMt,MMt,MMZ,MMt)
      ;

    diffMfin[20] =
      - prot0tt0t->M(0)
      + Mfin1(0,MMt,MMt,0,MMt)
      ;

    diffMfin[21] =
      - prot00WW0->M(0)
      + Mfin1(0,0,MMW,MMW,0)
      ;

    diffMfin[22] =
      - protWt000->M(0)
      + Mfin1(0,0,MMt,MMW,0)
      ;
    
    for(int i = 0; i < 23; i++)
      std::cout << "Test diffMfin[" << i << "]= " << diffMfin[i] << std::endl;

    diffVfin[0] =
      - prottH0H->Vxzuv(0)
      + Vfin1(MMt,MMH,0,MMH)
      ;

    diffVfin[1] =
      - prottZ0Z->Vxzuv(0)
      + Vfin1(MMt,MMZ,0,MMZ)
      ;

    for(int i = 0; i < 2; i++)
      std::cout << "Test diffVfin[" << i << "]= " << diffVfin[i] << std::endl;



    diffUfin[0] =
      - protHHttH->Uxzuv(0)
      + Ufin1(MMH,MMt,MMH,MMt)
      ;

    diffUfin[1] =
      - protHZttZ->Uxzuv(0)
      + Ufin1(MMH,MMt,MMZ,MMt)
      ;

    diffUfin[2] =
      - protHWt0W->Uxzuv(0)
      + Ufin1(MMH,MMt,0,MMW)
      ;

    diffUfin[3] =
      - protHttZt->Uuyxv(0)
      + Ufin1(MMZ,MMt,MMH,MMt)
      ;

    diffUfin[4] =
      - protHZttZ->Uyuzv(0)
      + Ufin1(MMZ,MMt,MMZ,MMt)
      ;

    diffUfin[5] =
      - protZ0tW0->Uxzuv(0)
      + Ufin1(MMZ,MMt,0,MMW)
      ;

    diffUfin[6] =
      - protHHttH->Uuyxv(0)
      + Ufin1(MMt,MMH,MMH,MMH)
      ;

    diffUfin[7] =
      - protHZttZ->Uzxyv(0)
      + Ufin1(MMt,MMH,MMZ,MMZ)
      ;

    diffUfin[8] =
      - protHWt0W->Uzxyv(0)
      + Ufin1(MMt,MMH,MMW,MMW)
      ;

    diffUfin[9] =
      - protZZttH->Uzxyv(0)
      + Ufin1(MMt,MMZ,MMZ,MMH)
      ;

    diffUfin[10] =
      - protZWt0W->Uzxyv(0)
      + Ufin1(MMt,MMZ,MMW,MMW)
      ;

    diffUfin[11] =
      - protZ0tW0->Uzxyv(0)
      + Ufin1(MMt,MMZ,0,0)
      ;

    diffUfin[12] =
      - protW00tW->Uuyxv(0)
      + Ufin1(MMt,0,MMW,MMW)
      ;

    diffUfin[13] =
      - protHWt0W->Uuyxv(0)
      + Ufin1(0,MMW,MMW,MMH)
      ;

    diffUfin[14] =
      - protZWt0W->Uuyxv(0)
      + Ufin1(0,MMW,MMW,MMZ)
      ;

    diffUfin[15] =
      - protWt000->Uzxyv(0)
      + Ufin1(0,MMW,0,MMt)
      ;

    // diffUfin[16] =
    //   - prot00WW0->Uuyxv(0)
    //   + Ufin1(0,MMW,0,0)
    //   ;

    diffUfin[16] =
      - prot0W00->Uxzuv(0)
      + Ufin1(0,MMW,0,0)
      ;

    for(int i = 0; i < 17; i++)
      std::cout << "Test diffUfin[" << i << "]= " << diffUfin[i] << std::endl;

    diffTfin[0] =
      - protHZttZ->Txuv(0)
      + Tfin1(MMH,MMZ,MMt)
      ;

    diffTfin[1] =
      - protHWt0W->Txuv(0)
      + Tfin1(MMH,MMW,0)
      ;

    diffTfin[2] =
      - protHttHt->Txuv(0)
      + Tfin1(MMH,MMt,MMH)
      ;

    diffTfin[3] =
      - prot0ttHt->Tuxv(0)
      + Tfin1(MMH,MMt,0)
      ;

    diffTfin[4] =
      - protZWt0W->Txuv(0)
      + Tfin1(MMZ,MMW,0)
      ;

    diffTfin[5] =
      - protHttZt->Tuxv(0)
      + Tfin1(MMZ,MMt,MMH)
      ;

    diffTfin[6] =
      - prot0ttZt->Tuxv(0)
      + Tfin1(MMZ,MMt,0)
      ;

    diffTfin[7] =
      - protHWt0W->Tvxu(0)
      + Tfin1(MMW,MMH,0)
      ;

    diffTfin[8] =
      - protZWt0W->Tvxu(0)
      + Tfin1(MMW,MMZ,0)
      ;

    diffTfin[9] =
      - protWt000->Txuv(0)
      + Tfin1(MMW,0,0)
      ;

    for(int i = 0; i < 10; i++)
      std::cout << "Test diffTfin[" << i << "]= " << diffTfin[i] << std::endl;


    diffSfin[0] =
      - protHHttH->Suxv(0)
      + Sfin1(MMH,MMH,MMt)
      ;

    diffSfin[1] =
      - protHZttZ->Suxv(0)
      + Sfin1(MMZ,MMH,MMt)
      ;

    diffSfin[2] =
      - protHZttZ->Svyz(0)
      + Sfin1(MMZ,MMZ,MMt)
      ;

    diffSfin[3] =
      - protHWt0W->Svyz(0)
      + Sfin1(MMW,MMW,MMt)
      ;

    diffSfin[4] =
      - protHWt0W->Suxv(0)
      + Sfin1(0,MMW,MMH)
      ;

    diffSfin[5] =
      - protZ0tW0->Suxv(0)
      + Sfin1(0,MMW,MMZ)
      ;

    diffSfin[6] =
      - protWt000->Svyz(0)
      + Sfin1(0,0,MMt)
      ;

    diffSfin[7] =
      - prot000->Suxv(0)
      + Sfin1(0,0,0)
      ;

    for(int i = 0; i < 8; i++)
      std::cout << "Test diffSfin[" << i << "]= " << diffSfin[i] << std::endl;


    // diffIfin[0] =
    //   - Tsil::I2(MMH,MMH,MMH,mu2)
    //   + Ifin1(MMH,MMH,MMH)
    //   ;

    // diffIfin[1] =
    //   - Tsil::I2(MMH,MMt,MMt,mu2)
    //   + Ifin1(MMH,MMt,MMt)
    //   ;

    // diffIfin[2] =
    //   - Tsil::I2(MMZ,MMZ,MMH,mu2)
    //   + Ifin1(MMZ,MMZ,MMH)
    //   ;

    // diffIfin[3] =
    //   - Tsil::I2(MMZ,MMt,MMt,mu2)
    //   + Ifin1(MMZ,MMt,MMt)
    //   ;

    // diffIfin[4] =
    //   - Tsil::I2(MMW,MMW,MMH,mu2)
    //   + Ifin1(MMW,MMW,MMH)
    //   ;

    // diffIfin[5] =
    //   - Tsil::I2(MMW,MMW,MMZ,mu2)
    //   + Ifin1(MMW,MMW,MMZ)
    //   ;

    // diffIfin[6] =
    //   - Tsil::I2(0,MMW,MMt,mu2)
    //   + Ifin1(0,MMW,MMt)
    //   ;

    // diffIfin[7] =
    //   - Tsil::I2(0,0,MMZ,mu2)
    //   + Ifin1(0,0,MMZ)
    //   ;

    // diffIfin[8] =
    //   - Tsil::I2(0,0,MMW,mu2)
    //   + Ifin1(0,0,MMW)
    //   ;  
  }
  std::complex<long double> mtopas()
  {    return mtopas(MMt,MMH, MMW, MMZ, mmu);  }
  std::complex<long double> mtopas(long double MMt,long double MMH,long double MMW,long double MMZ,long double mmu);

  
  std::complex<long double> mtopal()
  {    return mtopal(MMt,MMH, MMW, MMZ, mmu);  }  
  std::complex<long double> mtopal(long double MMt,long double MMH,long double MMW,long double MMZ,long double mmu);

  
  std::complex<long double> mtopalas()
  {    return mtopalas(MMt,MMH, MMW, MMZ, mmu);  }
  std::complex<long double> mtopalas(long double MMt,long double MMH,long double MMW,long double MMZ,long double mmu);
  

  std::complex<long double> mtopalal()
  {    return mtopalal(MMt,MMH, MMW, MMZ, mmu);  }  
  std::complex<long double> mtopalal(long double MMt,long double MMH,long double MMW,long double MMZ,long double mmu);


  // std::complex<long double> ttx()
  // {    
  //   std::complex<long double> res = 0;
  //   res += tt_x0(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x1(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x2(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x3(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x4(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x5(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x6(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x7(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x8(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x9(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x10(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x11(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x12(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x13(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x14(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x15(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x16(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x17(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x18(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x19(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x20(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x21(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x22(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x23(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x24(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x25(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x26(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x27(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x28(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x29(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x30(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x31(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x32(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x33(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x34(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x35(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x36(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x37(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x38(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x39(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x40(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x41(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x42(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x43(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x44(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x45(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x46(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x47(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x48(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x49(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x50(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x51(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x52(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x53(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x54(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x55(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x56(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x57(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x58(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x59(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x60(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x61(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x62(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x63(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x64(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x65(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x66(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x67(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x68(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x69(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x70(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x71(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x72(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x73(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x74(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x75(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x76(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x77(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x78(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x79(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x80(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x81(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x82(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x83(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x84(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x85(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x86(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x87(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x88(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x89(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x90(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x91(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x92(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x93(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x94(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x95(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x96(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x97(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x98(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x99(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x100(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x101(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x102(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x103(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x104(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x105(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x106(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x107(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x108(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x109(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x110(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x111(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x112(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x113(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x114(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x115(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x116(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x117(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x118(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x119(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x120(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x121(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x122(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x123(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x124(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x125(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x126(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x127(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x128(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x129(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x130(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x131(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x132(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x133(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x134(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x135(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x136(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x137(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x138(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x139(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x140(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x141(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x142(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x143(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x144(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x145(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x146(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x147(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x148(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x149(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x150(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x151(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x152(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x153(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x154(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x155(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x156(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x157(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x158(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x159(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x160(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x161(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x162(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x163(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x164(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x165(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x166(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x167(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x168(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x169(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x170(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x171(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x172(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x173(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x174(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x175(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x176(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x177(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x178(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x179(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x180(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x181(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x182(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x183(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x184(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x185(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x186(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x187(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x188(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x189(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x190(MMt, MMH, MMW, MMZ, mmu);
  //   res += tt_x191(MMt, MMH, MMW, MMZ, mmu);
  //   return res;
  // }  


  std::complex<long double> dalpha(long double MMt,long double MMH)
  {
    return 0;
  }
  
};


#endif  //  __TOP_HPP__
