/* Miscellaneous useful constants.  Users who #include this file
   should be careful to insure that there are no name conflicts with
   their own variables or functions. */

#ifndef TSIL_GLOBAL_H
#define TSIL_GLOBAL_H

/* Status codes (status) */
enum {UNEVALUATED, ANALYTIC, REAXIS, CONTOUR};

/* Function codes (whichFns) */
enum {ALL, STUM=0, STU, ST};

/* Size codes */
enum {LONG_DOUBLE, DOUBLE};

/* Number of each type of function (could just be hardwired) */
#define NUM_U_FUNCS 4
#define NUM_V_FUNCS 4
#define NUM_T_FUNCS 6
#define NUM_S_FUNCS 2
#define NUM_B_FUNCS 2

/* Enums for indexing */
enum {Bxz, Byu};
enum { xz,  yu};

enum {Svyz, Suxv};
/* enum {Tvyz, Tyzv, Tzyv, Tuxv, Txuv, Tvxu}; */
/* enum { vyz,  yzv,  zyv,  uxv,  xuv,  vxu}; */

/* Original: */
enum {Tvyz, Tuxv, Tyzv, Txuv, Tzyv, Tvxu};
enum { vyz,  uxv,  yzv,  xuv,  zyv,  vxu};


enum {Uzxyv, Uuyxv, Uxzuv, Uyuzv};
enum {Vzxyv, Vuyxv, Vxzuv, Vyuzv};
enum { zxyv,  uyxv,  xzuv,  yuzv};

#endif /* tsil_global.h */
