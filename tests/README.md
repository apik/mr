##Test set

Test prepared with the help of
[catch](https://github.com/philsquared/Catch) library. 

To list available tests use 
```
./testsuite -l
```
To run all test silently use without any arguments
```
./testsuite
```
To see not only errors but also succesfully passed tests use
```
./testsuite -s
```
Each test is supplied with its own tag, which can be used to run
selected test. For example to run tests with tag `[MS]` use
```
./testsuite "[MS]"
```
For other options use help with
```
./testsuite -h
```

##Available tests

* `tst-MWeqMZ.cpp`, tags:`[W]` `[Z]`

    Test for equality of corrections to MW and MZ in limit of equal masses MW=MZ, Mb=Mt

* `tst-RunDec.cpp`, tags:`[as]`

    Comaprison of running and decoupling aof strong coupling constant
    with `RunDec` package

* `tst-mPlanck.cpp`, tags: `[1307.3536v2]` `[1307.3536v4]` `[g1]` `[g2]` `[gs]` `[lam]` `[yt]`

   >   [Investigating the near-criticality of the Higgs boson.](https://inspirehep.net/record/1242456) 
   >   By Dario Buttazzo, Giuseppe Degrassi, Pier Paolo Giardino, Gian F. Giudice, Filippo Sala, Alberto Salvio, Alessandro Strumia.
   >   JHEP 1312 (2013) 089.
    
    Test of RGE running for all SM couplings up to a Plack
    scale. Different versions of arXiv versions of paper have
    different input at top mass scale. 

* `tst-BKKS_1205.2893.cpp`, tags: `[1205.1893]` `[Higgs]` `[top]`
   
   >   [Higgs Boson Mass and New Physics.](https://inspirehep.net/record/1114503)
   >   By Fedor Bezrukov, Mikhail Yu. Kalmykov, Bernd A. Kniehl, Mikhail Shaposhnikov.
   >   JHEP 1210 (2012) 140.
    
    Comparison mixed EWxQCD corrections for Higgs and top masses and Yukawa top and Higgs self-coupling


* `tst-JKV_0105304.cpp`, tags: `[0105304]` `[W]` `[Z]`
    
   >   [MS versus pole masses of gauge bosons: Electroweak bosonic two loop corrections.](https://inspirehep.net/record/557405)
   >   By F. Jegerlehner, M.Yu. Kalmykov, O. Veretin.
   >   Nucl.Phys. B641 (2002) 285-326.
   
    Comparison of pure bosonic corrections to masses of gauge bosons
    MW and MZ

* `tst-JKV_0212319.cpp`, tags: `[0212319]` `[W]` `[Z]`

   >   [MS-bar versus pole masses of gauge bosons. 2. Two loop electroweak fermion corrections.](https://inspirehep.net/record/605355)
   >   By F. Jegerlehner, Mikhail Yu. Kalmykov, O. Veretin.
   >   Nucl.Phys. B658 (2003) 49-112.

    Comparison of bosonic and fermionic corrections to masses of gauge bosons
    MW and MZ

* `tst-mOSmMS_QCD.cpp`, tags: `[MS]` `[OS]` `[QCD]`

   >   [Quark Mass Relations to Four-Loop Order in Perturbative QCD.](https://inspirehep.net/record/1342942)
   >   By Peter Marquard, Alexander V. Smirnov, Vladimir A. Smirnov, Matthias Steinhauser.
   >   Phys.Rev.Lett. 114 (2015) 14, 142002.

    Check for correctnes of logarithm part controlled by RGE of four-loop relation
    between OS and MS masses in QCD
