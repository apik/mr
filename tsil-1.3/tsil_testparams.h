/* Parameters used in test program: */

#ifndef TSIL_TESTPARAMS_H
#define TSIL_TESTPARAMS_H

/* NOTE: Generic cases should usually exceed these relative precisions
   with ease; these pass/warn/fail parameters are aimed to be
   "friendly" to the more difficult cases. Even so, you should expect
   to see some WARNs, and possibly FAILs, when running the test suite,
   depending on your platform.  */

#if defined(TSIL_SIZE_DOUBLE)
#define TSIL_PASS    3.e-8
#define TSIL_WARN    3.e-4
#define TSIL_PASS_V  3.e-5 /* The V functions have weaker requirements. */
#define TSIL_WARN_V  1.e-3

#else /* then it must be LONG */
#define TSIL_PASS    1.e-9
#define TSIL_WARN    1.e-6
#define TSIL_PASS_V  1.e-6  /* The V functions have weaker requirements. */
#define TSIL_WARN_V  1.e-4

#endif

#endif /* tsil_testparams.h */
