/* General header for "internal" use */

#include "tsil.h"        /* Contains types and user API prototypes */
#include "tsil_global.h" /* Misc. global objects */
#include "tsil_funcs.h"  /* Contains remaining functions */
#include "tsil_params.h" /* Parameters, e.g. for integration */
#include <fcntl.h>
#include <unistd.h>

#ifndef PI
#define PI 4.0L*TSIL_ATAN(1.0L)
#endif

#define Zeta2 1.644934066848226436472415166646025189219L
#define Zeta3 1.202056903159594285399738161511449990765L
#define TSIL_Infinity ((1.0L+ 1.0L*I)/0.0L)

enum {FALSE, TRUE};
enum {NO, YES};
