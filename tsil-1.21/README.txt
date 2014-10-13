                     ********************************
		            Welcome to TSIL
                     ********************************

Copyright (C) 2005 S.P. Martin and D.G. Robertson

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.  See the file LICENSE.txt for further
details.

Contents of this file:

I.   Overview
II.  Building TSIL
III. Using TSIL
IV.  The TSIL API
V.   Numerical Integration in TSIL
VI.  Using TSIL with Fortran


************************************************************************
I. Overview
************************************************************************

TSIL is a library of utilities for the numerical calculation of
dimensionally regularized two-loop self-energy integrals.  A
convenient basis for these functions is given by the integrals
obtained at the end of O.V. Tarasov's recurrence relation algorithm.
TSIL computes the values of all of these basis functions, for
arbitrary input masses and external momentum.  When analytical
expressions in terms of polylogarithms are available, they are
used. Otherwise, the evaluation proceeds by a Runge-Kutta integration
of the coupled first-order differential equations for the basis
integrals, using the external momentum invariant as the independent
variable. The code is written in C, and may be linked from C/C++ or
Fortran.

Authors:
S.P. Martin [spmartin AT niu.edu] 
Department of Physics
Northern Illinois University
DeKalb, IL 60115
USA

D.G. Robertson [drobertson AT otterbein.edu]
Department of Physics
Otterbein University
Westerville, OH 43081
USA

To cite this program, please reference the paper (referred to in this
document as MR05):

"TSIL: a program for the calculation of two-loop self-energy
integrals", by S.P. Martin and D.G. Robertson, [hep-ph/0501132].

Also, please reference the paper containing some of the results on
which it is based:

"Evaluation of two-loop self-energy basis integrals using differential
equations," by S.P. Martin
Phys. Rev. D 68, 075002 (2003) [hep-ph/0307101].

TSIL is available from:

       http://www.niu.edu/spmartin/TSIL
       http://faculty.otterbein.edu/drobertson/TSIL

Version number: 1.21


************************************************************************
II. Building TSIL
************************************************************************

TSIL can be compiled on any system that supports the GNU Compiler
Collection (gcc), the Intel C compiler (icc), or a similar C compiler
with support for complex mathematics.

To compile TSIL, edit the Makefile and choose:

1. The size of basic data types.  This is controlled by the compiler
   flag

	 -DTSIL_SIZE_<size>

   where <size> may be LONG or DOUBLE:

	 LONG	    Basic floating point type is long double
	 DOUBLE	    Basic floating point type is double

   Simply uncomment the line for the type you wish.  (If neither flag
   is given, LONG is selected.)  LONG is strongly recommended on
   systems where it is available.  (However, the gcc currently available for 
   Mac OS X will not work with LONG.) There is a speed penalty due to 
   the use of long double intrinsic functions, but it is minor and execution 
   times are in any case low (less than one second) on modern hardware.

   In your own programs that use TSIL, you may declare variables as

         TSIL_REAL
         TSIL_COMPLEX

   which will automatically correspond to the type selected when the
   library was compiled.  Note also that macros of the form

         TSIL_CLOG

   will automatically select the appropriate intrinsic function (in
   this case, either clog or clogl).  For a full list of intrinsics
   with this behavior, see the file tsil.h.

2. Compiler and optimization flags.

   Several sets are pre-defined in the Makefile; simply uncomment the
   appropriate one for your system if present.  TSIL is currently
   known to compile with gcc (under Linux or Mac OS X) and icc.  Other
   C compilers should work provided that complex mathematics is
   supported, but in this case you will need to explicitly set the
   compiler name and optimization flags.

   If you succeed in building TSIL on a new platform, the authors
   would be grateful to know it, and to learn of any special measures
   that were needed to compile it.

3. Install directories, if desired.

   You can set TSIL_LIBDIR and TSIL_INCDIR to point to directories
   where you would like the library and the TSIL header file,
   respectively, to be placed after compilation. If TSIL_LIBDIR and
   TSIL_INCDIR are not set, then the library and the header file will be 
   left after compilation in the directory where the sources reside, and 
   they can be moved by hand to an appropriate place. Standard directories
   that are automatically searched by compilers and linkers typically
   include /usr/lib and /usr/include, but you will need root access to
   write to these directories.  If you specify other directories not
   on the standard search path, note that when compiling your own code
   it will be necessary to specify these directories using the options
   -I<dir> and -L<dir>.  See the compiler/linker man pages for
   complete details. 

Once these choices have been made, simply type

       make

to build the library, along with the associated test and sample
programs.  After this command is complete, and after TSIL has been
tested with satisfactory results (see below), you can type

       make install

to install the library and header files in the specified locations.
Congratulations, TSIL is ready for action!

The end product intended for the user consists of the files:

       libtsil.a     The static TSIL archive (will be placed in
		     TSIL_LIBDIR upon "make install")

       tsil.h        The TSIL header file, must be included in any 
		     user code that uses TSIL (will be placed in
		     TSIL_INCDIR upon "make install")

       tsil	     Executable for basic computation (see below 
		     for details)

       tsil.tst	     Test program (see below for details)

In addition, the file

       tsil_fort.inc

may be useful if you are planning to call TSIL from Fortran.  See
section VII of this document for additional information on using TSIL
with Fortran.

It is strongly recommended that you run the test program after
compiling TSIL, to insure that correct results are being obtained.
The program

       tsil.tst

compares the output of TSIL_Evaluate to predefined results in 320 data
files, which are located in the directory TestData.  These data files
include cases representing all known analytic results as well as cases
requiring integration that have thresholds and pseudo-thresholds at
s=0 and at the final s.

tsil.tst takes a list of filenames as command line arguments, and
outputs, for each test case (i.e. for each file) either PASS, WARN or
FAIL.  If a WARN or FAIL results, the responsible functions are
printed, with both the expected and obtained values.  Assuming you
have just compiled TSIL, the easiest way to run the entire suite of
tests is via the command

       ./tsil.tst TestData/* > foo.txt

with the output redirected to the file foo.txt (recommended here due
to the large number of tests).  The end of this file will contain a
summary with the total number of PASS, FAIL and WARN results.

The pass/fail/warn criteria are controlled by macros TSIL_PASS,
TSIL_WARN, TSIL_PASS_V and TSIL_WARN_V, defined in tsil_testparams.h.
The first sets the maximum relative error allowed for the test to
pass; the second sets a lower error threshold below which the test is
deemed to fail.  A relative error between these two values results in
a warning.  (TSIL_PASS_V and TSIL_WARN_V define slightly less
stringent requirements for the V functions.)  As an example, for long
double data the default values are

       TSIL_PASS   = 1.e-9
       TSIL_WARN   = 1.e-6
       TSIL_PASS_V = 1.e-6
       TSIL_WARN_V = 1.e-4

Generic cases should usually exceed these relative precisions with
ease; these pass/warn/fail parameters are aimed to be "friendly" to
the more difficult cases. Even so, you may see some WARNs, and
possibly even FAILs, when running the test suite, depending on your
platform.

Users should not need to do this to insure correct functionality, but
the test program can be configured to evaluate subset cases (STU or ST
functions) by uncommenting the appropriate flag in the Makefile. In
this case the test program evaluates the selected subset case for each
set of input parameters, comparing the results to values taken from
the full data files.

The make command also produces the executable

       tsil

which implements the most basic TSIL calculation: it takes as
command-line arguments x, y, z, u, v, s, and Q^2 and prints the values
of all integral functions together with timing and other information.
As an example,

       ./tsil 1 2 3 4 5 10 1

evaluates all functions for x=1.0, y=2.0, z=3.0, u=4.0, v=5.0, s=10.0,
and Q^2=1.0 and prints the results to stdout.


************************************************************************
III. Using TSIL
************************************************************************

To use TSIL functions in your code, you must:

1. Include the header file tsil.h in any source file that makes use of
   TSIL data structures or functions, e.g. by adding the line

	#include "tsil.h"

   This is appropriate if the file tsil.h is located in the directory
   where the code is being compiled; if it has been placed in a
   standard location such as /usr/include, then

	#include <tsil.h>

   would work.  If it is a nonstandard directory <inc_dir>, the
   compiler option

        -I<inc_dir>

   will generically be necessary.  See the compiler man pages on your
   system for further details.

   Some potentially useful constants are also defined in the file

	tsil_global.h 

   which may be included by users wishing to access them.  Note that
   the constant names defined herein are not prefixed by TSIL_ and are
   therefore potentially subject to namespace collisions with objects
   in the user code.  Users wishing to make use of these constants
   should be careful to insure that their own variable and function
   names do not match any of those defined in tsil_global.h.

2. Link to the library at the end of the compilation process.  This is
   accomplished via the (linker) flag

	-ltsil

   If libtsil.a is not in a standard location (including the case
   where it is in the current directory), you will generally need to
   add the flag

        -L<lib_dir>

   where <lib_dir> is the directory in which libtsil.a may be found.
   If this is the current directory, then

	 -L.

   may be used.  Again, consult the compiler man pages for complete
   details on making user libraries available to the linker.

Complete details regarding the TSIL functions are given in section IV
of this document.  In the rest of the section we will discuss some
general points and exhibit some simple examples.

The basic data object in TSIL is a C struct with type name

       TSIL_DATA

that contains the parameter values x, y, z, u, v and Q^2 as well as
the 15 basis functions of types B, S, T, U, M.  Each basis function is
itself a struct containing its value, arguments, and various
unchanging coefficients used in computing its derivative.  Also
contained in the basic data struct are values of the integrals Tbar,
V, and bold versions of S, T, U, and V.  Definitions of all datatypes
are contained in the header file tsil.h.

In any program that calls TSIL functions requiring Runge-Kutta
evaluation, at least one of these high-level data objects must be
declared, e.g.

	TSIL_DATA foo;

More than one such object, and arrays of such objects, are allowed.
Users can of course access the items in the struct directly, though it
is recommended that the provided user interface routines be used.
These allow one to extract values of individual functions (or all of
them), set the values of the external parameters, and so on.  See
section VI for additional details.

Note that, as discussed above, the types TSIL_REAL and TSIL_COMPLEX
are available to insure that intrinsic floating-point types match
those used in compiling the library.

In the simplest application of TSIL, the parameters x, y, z, u, v and
Q^2 will be set using TSIL_SetParameters, the integrals for real s
evaluated using TSIL_Evaluate, and the results either printed using
TSIL_PrintData, or perhaps extracted by the calling program with the
command TSIL_GetFunction.  Generic code for this would look like:

	TSIL_DATA    foo;
	TSIL_REAL    x, y, z, u, v, s, qq;
	TSIL_COMPLEX integral1, integral2;
	...
	TSIL_SetParameters (&foo, x, y, z, u, v, qq);
	TSIL_Evaluate (&foo, s);
	integral1 = TSIL_GetFunction (&foo, <string1>);
	integral2 = TSIL_GetFunction (&foo, <string2>);
	...

where the strings <string1>, <string2> can each be one of

         "M", "Uzxyv", "Uuyxv", "Uxzuv", "Uyuzv", "Tvyz", "Tuxv",
	 "Tyzv", "Txuv", "Tzyv", "Tvxu", "Svyz", or "Suxv"

according to which of the basis integrals is desired. The strings can
also be one of

          "Vzxyv", "Vuyxv", "Vxzuv", "Vyuzv", "Bxz", or "Byu"

to access the functions V and the one-loop B functions.

The "bold" variants of the S, T, U and V functions can be accessed in
a similar way, e.g.

          integral3 = TSIL_GetBoldFunction (&foo, <string3>, n);

would return the coefficient of 1/\epsilon^n (for n=0,1, or 2) in the
bold-faced function corresponding to an appropriate <string3> from the
list above.

TSIL_Evaluate first decides whether the case at hand is known
analytically; if so, the basis functions are computed directly.  If
not, numerical integration is performed.  TSIL_Evaluate returns 1
(TRUE) for successful execution or 0 (FALSE) for error execution.  (A
warning message is printed if the external parameters correspond to
the unnatural threshold case discussed in [MR05].)  The data object
further contains a status parameter, accessible via the function
TSIL_GetStatus, which indicates how the master integral evaluation was
performed: either analytic, numerical integration along real axis, or
numerical integration along the displaced contour.

All integrals that are analytically known in terms of polylogarithms
can also be evaluated directly, without TSIL_SetParameters or
TSIL_Evaluate or TSIL_GetFunction. For example,

       TSIL_Manalytic (x,y,z,u,v,s,&res);

will return the int value 1 and set the variable res equal to
M(x,y,z,u,v) for the appropriate s, if it is analytically available,
and otherwise will return 0. Here x,y,z,u,v,qq are of type TSIL_REAL,
and s and res are of type TSIL_COMPLEX.  The functions

       TSIL_Sanalytic
       TSIL_Tanalytic
       TSIL_Tbaranalytic
       TSIL_Uanalytic
       TSIL_Vanalytic

have analogous behavior, except that they carry an additional argument
qq of type TSIL_REAL for the renormalization scale squared Q^2. For
example,

       TSIL_Uanalytic (x,y,z,u,s,qq,&res)}

will return the int value 1 and set the variable res equal to
U(x,y,z,u) for the appropriate s and Q^2, if it is analytically
available, and otherwise will return 0.

The other analytic functions (i.e. those that are known for arbitrary
parameter values and values of s) assign without pointers, for example

       res = TSIL_Bp (x,y,s,qq);

will set result equal to B(x',y) computed analytically for the
appropriate s and Q^2.

In addition to the evaluation for generic parameters described above,
TSIL provides functions for direct analytical evaluation of the vacuum
integrals A(x) and I(x,y,z), the one-loop integral B(x,y), as well as
various derivatives of these.

The standard output function is TSIL_PrintData, which prints all
function values on stdout.  An alternate format, designed so that
captured output can serve as valid input files for Mathematica, is
given by TSIL_PrintDataM.  Additional utilities allow the user to
extract individual basis functions or sets of functions to arrays.
Note that warning and error messages appear on stderr so they may be
redirected by the shell and examined separately.

A complete example program that uses TSIL to compute the two-loop pole
mass of a scalar particle with both cubic and quartic self-
interactions is given in

     scalarpole.c,

included with the distribution.  (Detailed background on this
calculation may be found in ref. [MR05].)  It takes parameter values
m^2, g, lambda and Q^2 as command-line inputs, in that order.  It then
computes the required basis functions at s=m^2, assembles the one- and
two-loop pole mass squared values as outlined in [MR05], and prints
the results to stdout.  Most of the basic functionality available in
TSIL is exhibited in scalarpole.c.

To compile this program using gcc, assuming tsil.h and libtsil.a are
present in the current directory, use e.g.

     gcc -o spole -DTSIL_SIZE_<size> scalarpole.c -L. -ltsil -lm

where <size> is either LONG or DOUBLE, and matches the size chosen
when libtsil.a was compiled.  If you used the default size when
compiling libtsil.a, then you may omit this flag.  This command
produces the executable spole, which can then be run as, e.g.

     ./spole 1 2 3 1

Note the use in scalarpole.c of TSIL_A, TSIL_I2p, TSIL_B, TSIL_Bp, and
TSIL_dBds to evaluate the functions A(x), I(x',x,x), B(x,x), B(x',x),
and the partial derivative of B(x,x) with respect to s, respectively.
In the evaluation of pi2, we arbitrarily chose to use

     TSIL_GetFunction(&result, "Bxz") 

where TSIL_B could also have been used.

For convenience there is also a struct of type TSIL_RESULT, which
contains only the parameter values and the results for the B, S, T,
Tbar, U, V, and M functions.  A function TSIL_CopyResult takes an
evlauated TSIL_DATA struct and copies the results into a specified
TSIL_RESULT. There is also a function TSIL_PermuteResult that can
permute the arguments of a TSIL_RESULT to produce another
TSIL_RESULT. See section IV below for complete details.


************************************************************************
IV. The TSIL Application Programmer Interface
************************************************************************

In this section we give the signatures of all TSIL functions, and
describe their operation.

Basic Evaluation and Extraction Functions
=========================================

Before evaluation, one of the following three functions should be called
to set the relevant squared mass and renormalization scale parameters.

1a. int TSIL_SetParameters (TSIL_DATA *foo,
                            TSIL_REAL x,
			    TSIL_REAL y,
			    TSIL_REAL z,
			    TSIL_REAL u,
			    TSIL_REAL v,
                            TSIL_REAL qq)
  
   Sets parameter values x, y, z, u, v, qq in the data object *foo.
   Return value is currently ignored.

1b. int TSIL_SetParametersSTU (TSIL_DATA *foo,
                               TSIL_REAL x,
			       TSIL_REAL z,
			       TSIL_REAL u,
			       TSIL_REAL v,
			       TSIL_REAL qq)
  
   Sets parameter values x, z, u, v, qq in the data object *foo, and
   selects evaluation of the STU subset case. Return value is
   currently ignored.

1c. int TSIL_SetParametersST (TSIL_DATA *foo,
                              TSIL_REAL x,
			      TSIL_REAL u,
			      TSIL_REAL v,
                              TSIL_REAL qq)
  
   Sets parameter values x, u, v, qq in the data object *foo and
   selects evaluation of the ST subset case. Return value is currently
   ignored.

   Subsequent calls of TSIL_SetParameters, TSIL_SetParametersSTU, or
   TSIL_SetParametersST to change one or more of x, y, z, u, v, qq are
   allowed, and to select the default (1a) or subset (1b) or (1c)
   evaluation modes. Each such call also resets the default values of
   parameters related to numerical integration (see section VII
   below).

---------------------------------------------------------------------

2. int TSIL_Evaluate (TSIL_DATA *foo, TSIL_REAL s)

   Evaluates all functions in *foo, including bold variants, at the
   specified value of s.  If *foo was initialized with
   TSIL_SetParametersSTU() or TSIL_SetParametersST(), then only the
   relevant subset of the basis functions is computed.  Return value
   is 1 (TRUE) for successful execution, 0 (FALSE) for error
   execution.

---------------------------------------------------------------------

3. int TSIL_GetStatus (TSIL_DATA *foo)

   Returns the evaluation status of *foo; either

	   0  Unevaluated
	   1  Evaluated by analytical formula
	   2  Evaluated by numerical integration along real s axis
	   3  Evaluated by numerical integration along displaced contour

   Note that a set of enum'ed constants (UNEVALUATED, ANALYTIC,
   REAXIS, CONTOUR) representing these are available in
   tsil_global.h.

---------------------------------------------------------------------

4. void TSIL_GetData (TSIL_DATA *foo,
	              const char *str,
		      TSIL_COMPLEX *res)

   Extracts a collection of functions to an array res[] according to
   the string str provided.  This may be one of

       "M", "U", "T", "S", "V", "B", or "TBAR"

   The user is responsible for insuring that the array res[] is of the
   correct size to hold all of the specified values.  The macros
   NUM_U_FUNCS, NUM_V_FUNCS, NUM_T_FUNCS, NUM_S_FUNCS, NUM_B_FUNCS are
   available and may be used in dimensioning such arrays.

   Example:
   ========
   TSIL_DATA foo;
   TSIL_COMPLEX uvals[NUM_U_FUNCS];
   ...
   TSIL_GetData (&foo, "U", uvals);

   Note that in subset evaluation (STU or ST), calling this function
   will generate an error message.  Instead, TSIL_GetFunction should
   be used to extract results.

---------------------------------------------------------------------

5. void TSIL_GetBoldData (TSIL_DATA *foo,
                          const char *str,
	 		  TSIL_COMPLEX res[][3])

   Similar to TSIL_ExtractData, but extracts a collection of "bold"
   functions to a two-dimensional array res[][3] according to the
   string provided.  As above, this may be one of

       "M", "U", "T", "S", "V", "B", or "TBAR"

   The user is responsible for insuring that the array res[][3] is of
   the correct size to hold all of the specified values.  

   Example:
   ========
   TSIL_DATA foo;
   TSIL_COMPLEX bolduvals[NUM_U_FUNCS][3];
   ...
   TSIL_GetBoldData (&foo, "U", bolduvals);

   Note that in subset evaluation (STU or ST), calling this function
   will generate an error message.  Instead, TSIL_GetBoldFunction
   should be used to extract results.

---------------------------------------------------------------------

6. TSIL_COMPLEX TSIL_GetFunction (TSIL_DATA *foo, const char *str)

   Returns a single function from *foo according to the string
   provided.  str may be any of

         "M",
	 "Uzxyv", "Uuyxv", "Uxzuv", "Uyuzv",
	 "Vzxyv", "Vuyxv", "Vxzuv", "Vyuzv",
	 "Tvyz", "Tuxv", "Tyzv", "Txuv", "Tzyv", "Tvxu",
	 "Svyz", "Suxv",
	 "Bxz", or "Byu"

   Example:
   ========
   TSIL_DATA foo;
   TSIL_COMPLEX tyzv;
   ...
   tyzv = TSIL_GetFunction (&foo, "Tyzv");

   Note that where there is only a single function of a given type
   available, only the first character in the specification string is
   relevant.  For example, for STU evaluation the single U-type function
   Uxzuv can be extracted using "Uxzuv" or "U".

---------------------------------------------------------------------

7. TSIL_COMPLEX TSIL_GetBoldFunction (TSIL_DATA *foo, 
                                      const char *str,
				      int n)			      

   Like TSIL_GetFunction but returns the coefficient of 1/epsilon^n in
   the "bold" variant function specified by the string str.  As above,
   str may be any of

         "M",
	 "Uzxyv", "Uuyxv", "Uxzuv", "Uyuzv",
	 "Vzxyv", "Vuyxv", "Vxzuv", "Vyuzv",
	 "Tvyz", "Tuxv", "Tyzv", "Txuv", "Tzyv", "Tvxu",
	 "Svyz", "Suxv",
	 "Bxz", or "Byu"

   Example:
   ========
   TSIL_DATA foo;
   TSIL_COMPLEX tyzv2;
   ...
   /* Returns power of 1/epsilon^2 in bold Tyzv: */
   tyzv2 = TSIL_GetBoldFunction (&foo, "Tyzv", 2);

   Note that where there is only a single function of a given type
   available, only the first character in the specification string is
   relevant.  E.g., for STU evaluation the single U-type function
   Uxzuv can be extracted using "Uxzuv" or "U".

---------------------------------------------------------------------

8. void TSIL_CopyResult (TSIL_DATA *foo, TSIL_RESULT *bar)

   Copies the parameters and integral values from an evaluated
   TSIL_DATA struct into a smaller TSIL_RESULT struct, for
   convenience.

---------------------------------------------------------------------

9. void TSIL_PermuteResult (TSIL_RESULT *in, int perm, TSIL_RESULT *out)

   Permutes the arguments of a TSIL_RESULT *in to produce a new
   TSIL_RESULT *out. The integral functions are reorganized
   appropriately.

   Allowed values of perm are:

   0 (or NOSWAP)  - No permutation, just copies the result
   1 (or XYandZU) - Permute x <-> y and z <-> u
   2 (or XZandYU) - Permute x <-> z and y <-> u
   3 (or XUandYZ) - Permute x <-> u and y <-> z

---------------------------------------------------------------------

I/O and Related Functions
=========================

1. void TSIL_PrintStatus (TSIL_DATA *foo)

   Prints to stdout the evaluation status of *foo, whether
   unevaluated, evaluated by analytical formula, evaluated by
   numerical integration along real s axis, or evaluated by numerical
   integration along displaced contour.

---------------------------------------------------------------------

2. void TSIL_PrintData (TSIL_DATA *foo)

   Prints to stdout the values of all integral functions associated
   with *foo, including Tbar, V, and "bold" variants. In STU or
   ST evaluation only the appropriate subset of functions is
   displayed.

---------------------------------------------------------------------

3. void TSIL_WriteData (FILE *fp, TSIL_DATA *foo)

   Like TSIL_PrintData, but writes to the file pointer fp rather than
   stdout.  The user is responsible for properly initializing the file
   pointer.

---------------------------------------------------------------------

4. void TSIL_PrintDataM (TSIL_DATA *foo)

   Same as TSIL_PrintData, but uses a format that is valid Mathematica
   input.

---------------------------------------------------------------------

5. void TSIL_WriteDataM (FILE *fp, TSIL_DATA *foo)

   Same as TSIL_WriteData, but uses a format that is valid Mathematica
   input.

---------------------------------------------------------------------

6. void TSIL_cprintf (TSIL_COMPLEX)

   Prints a complex value to stdout.  (No newline is appended.)

---------------------------------------------------------------------

7. void TSIL_cprintfM (TSIL_COMPLEX)

   Prints a complex value to stdout in a form digestible by
   Mathematica.  (No newline is appended.)

---------------------------------------------------------------------

8. void TSIL_Error (char *func, char *msg, int errcode)

   Prints an error message to stderr and exits with status errcode.
   The message contains the function in which the error was generated
   (func) and an error-specific message (msg) of the user's choice.

---------------------------------------------------------------------

9. void TSIL_Warn (char *func, char *msg)

   Prints a warning message to stderr; execution continues normally.
   The message contains the function in which the error was generated
   (func) and an error-specific message (msg) of the user's choice.

---------------------------------------------------------------------

Utilities
=========

1. void TSIL_ResetStepSizeParams (TSIL_DATA *foo,
			          TSIL_REAL precisionGoal,
                                  int nStepsStart,
				  int nStepsMaxCon,
				  int nStepsMaxVar,
				  int nStepsMin)

   Allows modification of integration parameters.  In *foo, sets
   precisionGoal, nStepsStart, nStepsMaxCon, nStepsMaxVar, nStepsMin
   to the specified values.  A subsequent call to TSIL_SetParameters
   will cause these to reset to their default values.

---------------------------------------------------------------------

2. int TSIL_IsInfinite (TSIL_COMPLEX z)

   Returns 1 (TRUE) or 0 (FALSE) according to whether z is infinite.

---------------------------------------------------------------------

3. int TSIL_DataSize (void)

   Returns a code indicating the size of intrinsic datatypes.
   Possible values are:

	    0   long double
	    1   double

   Note that there are enum'ed constants (LONG_DOUBLE, DOUBLE)
   representing these available in tsil_global.h.

---------------------------------------------------------------------

4. void TSIL_PrintInfo (void)

   Prints intrinsic datatype used by TSIL to stdout.

---------------------------------------------------------------------

5. void TSIL_PrintVersion (void)

   Prints the current version of TSIL to stdout.

---------------------------------------------------------------------

Analytic Cases
==============

1. TSIL_COMPLEX TSIL_Dilog (TSIL_COMPLEX z)

   Returns the dilogarithm of z.

---------------------------------------------------------------------

2. TSIL_COMPLEX TSIL_Trilog (TSIL_COMPLEX z)

   Returns the trilogarithm of z.

---------------------------------------------------------------------

3. TSIL_REAL TSIL_A (TSIL_REAL x, TSIL_REAL qq)

   Returns the one-loop vacuum function A(x) with renormalization
   scale squared equal to qq.

---------------------------------------------------------------------

4. TSIL_REAL TSIL_Ap (TSIL_REAL x, TSIL_REAL qq)

   Returns the derivative of the one-loop vacuum function A(x) with
   respect to x, with renormalization scale squared equal to qq.

---------------------------------------------------------------------

5. TSIL_COMPLEX TSIL_Aeps (TSIL_REAL x, TSIL_REAL qq)

   Returns A_epsilon(x) with renormalization scale squared equal to
   qq.

---------------------------------------------------------------------

6. TSIL_COMPLEX TSIL_B (TSIL_REAL x,
		        TSIL_REAL y,
			TSIL_COMPLEX s,
			TSIL_REAL qq)

   Returns the one-loop self-energy function B(x,y) with external
   momentum invariant s and renormalization scale squared equal to qq.

---------------------------------------------------------------------

7. TSIL_COMPLEX TSIL_Bp (TSIL_REAL x,
			 TSIL_REAL y,
			 TSIL_COMPLEX s,
			 TSIL_REAL qq)

   Returns the function B(x',y) for external momentum invariant s and
   renormalization scale squared equal to qq.

---------------------------------------------------------------------

8. TSIL_COMPLEX TSIL_dBds (TSIL_REAL x,
			   TSIL_REAL y,
			   TSIL_COMPLEX s, 
			   TSIL_REAL qq)

   Returns the derivative of B(x,y) with respect to s, for external
   momentum invariant s and renormalization scale squared equal to qq.

---------------------------------------------------------------------

9. TSIL_COMPLEX TSIL_Beps (TSIL_REAL x,
			   TSIL_REAL y,
			   TSIL_COMPLEX s, 
                           TSIL_REAL qq)

   Returns the function B_epsilon(x,y) for external momentum invariant
   s and renormalization scale squared equal to qq.

---------------------------------------------------------------------

10. TSIL_COMPLEX TSIL_I2 (TSIL_REAL x,
     			  TSIL_REAL y,
			  TSIL_REAL z,
			  TSIL_REAL qq)

   Returns the two-loop vacuum integral I(x,y,z), with renormalization
   scale squared equal to qq.

---------------------------------------------------------------------

11. TSIL_COMPLEX TSIL_I2p (TSIL_REAL x,
			   TSIL_REAL y,
			   TSIL_REAL z,
			   TSIL_REAL qq)

   Returns the first derivative I(x',y,z), with renormalization scale
   squared equal to qq.

---------------------------------------------------------------------

12. TSIL_COMPLEX TSIL_I2p2 (TSIL_REAL x,
			    TSIL_REAL y,
			    TSIL_REAL z,
			    TSIL_REAL qq)

   Returns the second derivative I(x'',y,z), with renormalization
   scale squared equal to qq.

---------------------------------------------------------------------

13. TSIL_COMPLEX TSIL_I2pp (TSIL_REAL x,
			    TSIL_REAL y,
			    TSIL_REAL z,
			    TSIL_REAL qq)

   Returns the mixed second derivative I(x',y',z), with
   renormalization scale squared equal to qq.

---------------------------------------------------------------------

14. TSIL_COMPLEX TSIL_I2p3 (TSIL_REAL x,
			    TSIL_REAL y,
			    TSIL_REAL z,
			    TSIL_REAL qq)

   Returns the third derivative I(x''',y,z), with renormalization
   scale squared equal to qq.

---------------------------------------------------------------------

15. int TSIL_Sanalytic (TSIL_REAL x,
		        TSIL_REAL y,
			TSIL_REAL z,
			TSIL_COMPLEX s,
			TSIL_REAL qq,
			TSIL_COMPLEX *res)

    If an analytic result for S(x,y,z) for external momentum invariant
    s is available, return 1 (TRUE) and evaluate it for
    renormalization scale squared qq.  The result is placed in *res.

    If no analytical result for S is available, return 0 (FALSE).  In
    this case *res is unchanged.

---------------------------------------------------------------------

16. int TSIL_Tanalytic (TSIL_REAL x,
		        TSIL_REAL y,
			TSIL_REAL z,
			TSIL_COMPLEX s,
			TSIL_REAL qq,
			TSIL_COMPLEX *res)

    If an analytic result for T(x,y,z) for external momentum invariant
    s is available, return 1 (TRUE) and evaluate it for
    renormalization scale squared qq.  The result is placed in *res.

    If no analytical result for T is available, return 0 (FALSE).  In
    this case *res is unchanged.

---------------------------------------------------------------------

17. int TSIL_Tbaranalytic (TSIL_REAL x,
			   TSIL_REAL y,
			   TSIL_REAL z,
			   TSIL_COMPLEX s,
			   TSIL_REAL qq,
			   TSIL_COMPLEX *res)

    If an analytic result for Tbar(x,y,z) for external momentum
    invariant s is available, return 1 (TRUE) and evaluate it for
    renormalization scale squared qq.  The result is placed in *res.

    If no analytical result for Tbar is available, return 0 (FALSE).
    In this case *res is unchanged.

---------------------------------------------------------------------

18. int TSIL_Uanalytic (TSIL_REAL x,
    			TSIL_REAL y,
			TSIL_REAL z,
			TSIL_REAL u,
			TSIL_COMPLEX s,
			TSIL_REAL qq,
			TSIL_COMPLEX *res)

    If an analytic result for U(x,y,z,u) for external momentum
    invariant s is available, return 1 (TRUE) and evaluate it for
    renormalization scale squared qq.  The result is placed in *res.

    If no analytical result for U is available, return 0 (FALSE).  In
    this case *res is unchanged.

---------------------------------------------------------------------

19. int TSIL_Vanalytic (TSIL_REAL x,
    			TSIL_REAL y,
			TSIL_REAL z,
			TSIL_REAL u,
			TSIL_COMPLEX s,
			TSIL_REAL qq,
			TSIL_COMPLEX *res)

    If an analytic result for V(x,y,z,u) for external momentum
    invariant s is available, return 1 (TRUE) and evaluate it for
    renormalization scale squared qq.  The result is placed in *res.

    If no analytical result for V is available, return 0 (FALSE).  In
    this case *res is unchanged.

---------------------------------------------------------------------

20. int TSIL_Manalytic (TSIL_REAL x,
    			TSIL_REAL y,
			TSIL_REAL z,
			TSIL_REAL u,
			TSIL_REAL v,
			TSIL_COMPLEX s,
			TSIL_REAL qq,
			TSIL_COMPLEX *res)

    If an analytic result for M(x,y,z,u,v) for external momentum
    invariant s is available, return 1 (TRUE) and evaluate it for
    renormalization scale squared qq.  The result is placed in *res.
    If no analytical result for M is available, return 0 (FALSE).  In
    this case *res is unchanged.

---------------------------------------------------------------------


Fortran Interface
=================

subroutine tsilfevaluate (REAL x,
			  REAL y,
			  REAL z,
			  REAL u,
			  REAL v,
			  REAL qq,
			  REAL s)

Wrapper for TSIL_Evaluate, to be called from Fortran.  See "Using TSIL
with Fortran" below for complete details.

---------------------------------------------------------------------


************************************************************************
V. Numerical Integration in TSIL
************************************************************************

As discussed in [MR05], we generically use a 6-stage, 5th-order
embedded Runge-Kutta scheme with coefficients given by Cash and Karp.
This gives an error estimate for each dependent variable at each step
and thus allows estimation of the step size necessary to achieve a
given precision.

In addition, whenever the final value of s is at (or near) a threshold
or pseudo-threshold, we switch for the final leg of the numerical
integration to a 5-stage, 4th-order scheme.  This is slower and does
not provide an error estimate, but it avoids evaluation of basis
function derivatives at the endpoint, where they may be singular.

The "nearness" criterion is controlled by the constant THRESH_CUTOFF,
which is defined in tsil_params.h and has default value 0.025.  (This
is in dimensionless units; prior to numerical integration all
dimensionful quantities are rescaled by the largest of x,y,z,u,v, and
|s|.)  That is, if the final value of s is within 0.025 of the nearest
threshold or pseudo-threshold, the special 5-stage, 4th order
integration is used for the last leg of integration, otherwise the
normal 6-stage scheme is used throughout.  Users can adjust this value
if they find their results to be excessively sensitive to thresholds
or pseudo-thresholds outside this limit.

For the generic rectangular integration contour described in [MR05],
the displacement along the imaginary axis of the complex s plane is
set by IM_DISPL, also defined in tsil_params.h.  Its default value is
0.2, again in dimensionless units.  Of course all results should be
independent of this value, so users can change this as a probe of
numerical stability if unusual results are observed.

Finally let us describe the parameters associated with adaptive
stepsize control.  These are also all defined in tsil_params.h, and,
along with the intrinsic data size chosen, they have the most
significant effect on execution speed and accuracy.

These are realized as members of the TSIL_DATA struct, with names:

* precisionGoal: This is \delta_P in eq. (3.14) of [MR05].  (We use a
                 safety factor S=0.9.)  If the maximum estimated error
                 for any dependent variable exceeds \delta_P
                 multiplied by the increment of that variable for that
                 step, and also exceeds the relative precision of the
                 computer arithmetic times the absolute value of that
                 variable, then the step is retried with a smaller
                 step size, unless the step size would become smaller
                 than specified below. Also, after a successful step,
                 the size for the next step is chosen according to
                 eq. (3.14), unless it would exceed the
                 amount specified below.  (Defaults: 10^{-12} for long
                 double data, 5 x 10^{-11} for double data.)

* nStepsStart:   For each leg of the contour, the initial step size is
                 chosen so that there would be this many steps if the
                 step size did not change. (Default: 500)

* nStepsMin:     The maximum allowed step size on a leg of the contour
                 with dimensionless (rescaled) independent variable 
		 length L is given by L/nStepsMin. (Default: 100)

* nStepsMaxCon,
  nStepsMaxVar:  The minimum allowed step size on a leg of the contour
		 with dimensionless independent variable length L is 
		 given by L/(nStepsMaxCon + L*nStepsMaxVar).
		 (Defaults: 10000, 10000)

The step size is not allowed to increase by more than a factor of 1.5
or decrease by more than a factor of 2 after each step or attempted
step.  Note that by setting precisionGoal to 0, one can arrange that
the total number of steps on each leg tends to nStepsMaxCon +
L*nStepsMaxVar. If instead one sets precisionGoal to a very large
number, the number of steps will tend to nStepsMin.

The default values have been found to give good results for a wide
variety of different choices of input parameters, for the integration
variables used in the program.  However, to deal with exceptional
situations, they can be reset at run time with the function
TSIL_ResetStepSizeParams, after calling TSIL_SetParameters and before
calling TSIL_Evaluate.

Note that in general, long double data (typically with 63 or more bits
of relative precision) gives results with relative accuracies better
than 10^{-10} for generic cases, but sometimes somewhat worse in cases
with large mass hierarchies, and in some particularly difficult cases
significantly worse.  (The function V(x,y,z,u) for very small but
non-zero y can be particularly sensitive to roundoff errors, since the
individual terms in its evaluation are proportional to 1/y and yet it
is only logarithmically divergent as y approaches 0.)  The user should
consider modifying the default parameters of the program if
significant sensitivity to parameters is expected (or observed), or if
speed is an overriding concern.


************************************************************************
VI. Using TSIL with Fortran
************************************************************************

It is possible to use TSIL with Fortran, and some basic utilities for
this are included with the package.  Essential functionality is
provided by a wrapper function tsilfevaluate (defined in the file
fevaluate.c), which is called as a subroutine from a Fortran program.
This subroutine implements the most general TSIL calculation: it takes
as arguments x,y,z,u,v,Q^2, and s and returns the values of all basis
functions, including Tbar, V, and "bold" functions.

These results are returned to the calling program in a COMMON block,
which corresponds to a special C struct used by tsilfevaluate.  This
COMMON block contains a number of pre-defined arrays that hold the
various function values.  Definitions of the COMMON block and
subsidiary arrays are given in the header file

       tsil_fort.inc

which users should INCLUDE in their Fortran program.

IMPORTANT NOTE: It is crucial that the type sizes match between the
calling Fortran program and TSIL, because the COMMON block must have
the same memory "footprint" as the C struct defined in the wrapper
code.  In addition, the arguments to tsilfevaluate must be
consistently defined in both fevaluate.c and the calling Fortran
program.  Type mismatches will often result in non-fatal errors.
Users should be aware that corresponding types may not actually exist,
depending on platform and compiler.  (For example, on a Pentium (IA32)
CPU, gcc supports 8-byte doubles and 12-byte long doubles, but g77
only supports REAL*8.)  It is strongly recommended that exact type
sizes for your system be determined before building TSIL for use with
Fortran.

The Fortran include file

       tsil_fort.inc

assumes that the basic floating-point type in the Fortran program is
REAL*8, assumed to be equivalent to C type double.  (This is true on
all systems of which the authors are aware, but anomalies may exist.)
The associated complex type is of course COMPLEX*16.  These types must
match those defined in the file

       tsil_fortran.h

which contains definitions for the C wrapper program (fevaluate.c).
Specifically, the macro

       TSIL_REAL_F

should be set to the C type corresponding to the Fortran floating-
point type used (default: double).  Furthermore,

       TSIL_COMPLEX_F

should be set to the C type corresponding to the Fortran complex type
(default: double complex).  Finally,

       TSIL_REAL_SIZE

should be set to the number of bytes in the basic floating-point type
(default: 8 for REAL*8/double).  This value is used when tsilfevaluate
is called to test whether there may be a type mismatch.  (If a
possible mismatch is detected, a warning message is printed but
execution continues.)  However, this test is rather simple and should
not be relied upon to detect all possible errors.

Note that in principle these types need not match the basic types used
in the main TSIL routines.  When the wrapper calls TSIL_Evaluate, the
arguments are cast to the appropriate type used in the main TSIL
library (i.e., to TSIL_REAL), and upon return the results are cast
back to the type TSIL_COMPLEX_F, appropriate for return to the calling
Fortran program.  Thus it can be that the Fortran program uses a
"smaller" type, which is cast by the wrapper to a more precise type
for evaluation, and then cast back to the earlier type for return.

For example, in the case mentioned parenthetically above (IA32 CPU
with gcc/g77) the 8-byte type is the only option for the struct/COMMON
block that enables communication between Fortran and the wrapper
program.  Thus the Fortran code would use REAL*8, and TSIL_REAL_F and
TSIL_COMPLEX_F would be set to double and double complex,
respectively.  However, the main TSIL library could be built with
either double or long double.  In the latter case, the input
parameters would be cast to long double for evaluation, giving
additional accuracy in the results.


Sample Code
===========

A Fortran program fragment that uses these utilities is shown below. 

       PROGRAM ftest
c      Includes array and COMMON definitions:
       INCLUDE 'tsil_fort.inc'

c      (Code setting values of x,y,z,u,v,qq,s not shown)
       ...
c      Evaluate basis integrals:
       CALL tsilfevaluate(x,y,z,u,v,qq,s)

c      Print a representative value:
       PRINT *, U(xzuv)
       ...

Note in the last line the use of the INTEGER xzuv.  This is one of a
set of integer variables, also defined in tsil_fort.inc, that allow
items in the COMMON block arrays to be referred to by name.

A complete functional sample fortran program is available in ftest.f. It 
can be compiled by

     g77 ftest.f -L. -ltsil

Running ./a.out then prints out the values of the integral functions
for the case x, y, z, u, v, s, qq = 1, 2, 3, 4, 5, 10, 1


Compiling with g77 under Mac OS X
=================================

Note that when using g77 under Mac OS X, the flags

     -lmx -lcc_dynamic

are required in addition to the flags that link the tsil archive.
Thus to compile the test program ftest.f one would typically use

     g77 ftest.f -lmx -lcc_dynamic -L. -ltsil

assuming libtsil.a is in the current directory (.).  See the g77
documentation for additional information.


Other General Fortran Issues
============================

The provided wrapper code can serve as a model for users wishing to
write their own interface routines with additional functionality.  In
addition to the type size issue discussed above, the following general
points should be noted:

* In Fortran arguments are passed by reference, i.e. as pointers in C.
  Thus any C function that is called from Fortran should have pointer
  arguments.

* Fortran compilers typically append an underscore to all external
  names, to identify them as Fortran constructs.  (This is important
  because of the call-by-reference convention.)  Thus the Fortran
  subroutine tsilfevaluate is known as tsilfevaluate_ in C, and the
  COMMON block results corresponds to the struct results_ in C.
  Fortran names that already contain an underscore typically have
  two underscores appended.  Unfortunately, these conventions are
  not uniform among Fortran compilers, so testing is required in any
  particular case.  The utilities provided with TSIL have been tested
  with the GNU compilers gcc/g77 and the Intel compilers icc/ifort.

* Fortran and C store the elements of multi-dimensional arrays in a
  different order (column-major and row-major ordering, respectively).
  Thus a two-dimensional array generated in C must be reorganized for
  proper indexing in Fortran.  This reorganization amounts to
  "transposing" the arrays, interpreted as matrices.  Note also that
  indexing conventions differ: Fortran array indices start at 1 while
  C indices start at 0.

************************************************************************
End of README.txt
************************************************************************
