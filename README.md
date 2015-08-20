mr
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
    $ make

or for prepared package

    $ curl -OL
    $ tar -xf 
    $ cd mr-*.*
    $ ./configure
    $ make

To enable *Mathematica* interface coigure **mr** with full path to
    directory containing both of **mcc** and **mprep** tools. For example `/Applications/Mathematica.app/SystemFiles/Links/MathLink/DeveloperKit/MacOSX-x86-64/CompilerAdditions` 

    $ ./configure
    --with-mcc-path=/Applications/Mathematica.app/SystemFiles/Links/MathLink/DeveloperKit/MacOSX-x86-64/CompilerAdditions

## Usage

Just include in your code

    #include <mr.hpp>
    using namespace mr;

* * * * *

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
