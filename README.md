mr [![Build Status](https://travis-ci.org/apik/mr.svg)](https://travis-ci.org/apik/mr)
==

**Mathching &amp; Running** is a set of tools for NNLO Standard Model
vacuum stability analysis. It includes three-loop
renormalization-group equations for Standard Model couplings with
four-loop additions and two-loop thershold corrections.


## Installation

for git version

    $ git clone https://github.com/apik/mr.git
    $ cd mr
    $ autoreconf -i
    $ ./configure
    $ make && make install

or for prepared package

    $ curl -OL https://github.com/apik/mr/releases/download/v1.3.2/mr-1.3.2.tar.gz
    $ tar -zxf mr-1.3.2.tar.gz
    $ cd mr-1.01
    $ ./configure
    $ make && make install

To enable *Mathematica* interface configure **mr** with full path to
    directory containing both of **mcc** and **mprep** tools. For example 

    $ ./configure --with-mcc-path=/Applications/Mathematica.app/SystemFiles/Links/MathLink/DeveloperKit/MacOSX-x86-64/CompilerAdditions

After successfull compilation `libmr` library will be prepared in directory `mr/` and after successfull installation all header files will be placed to `$PREFIX/include` and library to `$PREFIX/lib`. Prefix may be set at coonfigure step with `./configure --prefix=<PATH>`.

## Usage

Just include in your code

    #include <mr.hpp>
    using namespace mr;

And compile with flags

    $ pkg-config --cflags --libs mr

Sample applications are in `examples` directory


### Tests

Due to differen architectures for package installation there is a set
of tests. Running

     $ ./testsuite

in `tests` directory prepare results of numerous tests.

* * * * *
### Package structure
```
boost               - OdeInt library (part of boost)
examples            - examples of library usage in C++ applications
math                - Mathematica interface
 |_ mrimpl.cpp      - mathlink functions implementation
 |_ mr.m            - test file for Mathematica interface
 |_ mr.tm           - mathlink definitions
mr                  - mr library source files
 |_ aGF/            - EW corrections to relation between Gf and aEM
 |_ dr/             - EW corrections to delta-r
 |_ gl/             - gaugeless limit for MS mass in terms of OS 
 |_ mMSmOS          - MS masses in tems of OS masses
 |_ mOSmMS/         - OS masses in tems of MS masses
 |_ yu/             - EW corrections to couplings
 |_ yugl/           - EW corrections to couplings in gaugless limit
 |_ HH.cpp          - corrections to X_H and Y_H source
 |_ HH.hpp          - corrections to X_H and Y_H header
 |_ WW.cpp          - corrections to X_W and Y_W source
 |_ WW.hpp          - corrections to X_W and Y_W source
 |_ ZZ.cpp          - corrections to X_Z and Y_Z source
 |_ ZZ.hpp          - corrections to X_Z and Y_Z source
 |_ alphaGF.cpp     - running aEM from Gf and pole masses source 
 |_ alphaGF.hpp     - running aEM from Gf and pole masses header
 |_ alphas.cpp      - QCD coupling evolution and matching source 
 |_ alphas.hpp      - QCD coupling evolution and matching header
 |_ base.hpp        - base class for all EW corrections classes
 |_ bb.cpp          - corrections to X_b and Y_b source
 |_ bb.hpp          - corrections to X_b and Y_b source
 |_ constants.hpp   - numerical constants and PDG input
 |_ dr.cpp          - Corrections to delta-r source
 |_ dr.hpp          - Corrections to delta-r header
 |_ logger.hpp      - simple logging facility
 |_ mOS2mMSQCD.cpp  - QCD corrections to relation between OS and MS masses
 |_ mr.hpp          - main header to include in user progrm
 |_ operators.hpp   - arithmetic operations with complex types
 |_ p2ms.cpp        - extraction running couplings from OS input source
 |_ p2ms.hpp        - extraction running couplings from OS input header
 |_ smRGE.cpp       - SM beta functions and evolution source
 |_ smRGE.hpp       - SM beta functions and evolution header
 |_ sminput.hpp     - simple type for manipulation with input
 |_ tdecl.hpp       - declaration of types used in library
 |_ timer.hpp       - simple timer implementation
 |_ tsil.cpp        - TSIL integral library interface source
 |_ tsil.hpp        - TSIL integral library interface header
 |_ tt.cpp          - Corrections to X_t and Y_t source
 |_ tt.hpp          - Corrections to X_t and Y_t header
tests               - tests for package consistency check
```

* * * * *
## References

### SM running

*  *[Anomalous dimensions of gauge fields and gauge coupling
beta-functions in the Standard Model at three loops.](http://inspirehep.net/record/1193366)
By A.V. Bednyakov, A.F. Pikelner, V.N. Velizhanin.
JHEP 1301 (2013) 017.*
*  *[Yukawa coupling beta-functions in the Standard Model at three
loops.](http://inspirehep.net/record/1208862)
By A.V. Bednyakov, A.F. Pikelner, V.N. Velizhanin.
Phys.Lett. B722 (2013) 336-340.*
*  *[Higgs self-coupling beta-function in the Standard Model at three
loops.](http://inspirehep.net/record/1224266)
By A.V. Bednyakov, A.F. Pikelner, V.N. Velizhanin.
Nucl.Phys. B875 (2013) 552-565.*

### Threshold corrections
*  *[Two-loop electroweak threshold corrections in the Standard Model.](http://inspirehep.net/record/1351233)
By Bernd A. Kniehl, Andrey F. Pikelner, Oleg L. Veretin.
[arXiv:1503.02138 [hep-ph]].*

### Additional tools 

Two-loop massive self-energy integrals evaluated numerically using TSIL
library:

* *[TSIL: A Program for the calculation of two-loop self-energy
integrals.](http://inspirehep.net/record/675010)
By Stephen P. Martin, David G. Robertson.
Comput.Phys.Commun. 174 (2006) 133-151.*

Effective potential up to three loops in SM in gaugless limit:

*  *[Three-loop Standard Model effective potential at leading order in
strong and top Yukawa couplings.](http://inspirehep.net/record/1262358)
By Stephen P. Martin.
Phys.Rev. D89 (2014) 1, 013003.*

Four-loop QCD corrections to strong and Yukawa couplings running:

*  *[The Four loop beta function in quantum chromodynamics.](http://inspirehep.net/record/439866)
By T. van Ritbergen, J.A.M. Vermaseren, S.A. Larin.
Phys.Lett. B400 (1997) 379-384.*
*   *[The four loop quark mass anomalous dimension and the invariant
quark mass.](http://inspirehep.net/record/441078)
By J.A.M. Vermaseren, S.A. Larin, T. van Ritbergen.
Phys.Lett. B405 (1997) 327-333.*
