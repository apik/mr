---
layout: default
---
	
[mr:Matching and Running](http://apik.github.io/mr) is a C++ package for NNLO Standard Model stability analysis. It includes full two-loop electroweak threshold corrections connecting input in terms of pole masses with running couplings and three-loop renormalization group equations withh additional higher order QCD corrections for evolution of running couplings up to a needed scale.

<iframe src="http://ghbtns.com/github-btn.html?user=apik&amp;repo=mr&amp;type=watch&amp;count=true&amp;size=large"
  allowtransparency="true" frameborder="0" scrolling="0" width="170" height="30"></iframe><br/>


## Sample plots


<figure>
  <img src="plots/gauge123.png" alt="Running SM couplings up to a Planck scale">
  <figcaption>Running SM couplings up to a Planck scale. <a href="plots/gauge123.svg">[SVG]</a>,<a href="plots/gauge123.svg">[EPS]</a></figcaption>	
</figure>	      	


## Installation

From github repository using git:

~~~
$ git clone https://github.com/apik/mr.git
$ cd mr
~~~

Downloading tarball and extracting:

~~~
$ git clone https://github.com/apik/mr.git
$ cd mr
~~~

{% highlight ruby %}
def print_hi(name)
puts "Hi, #{name}"
end
print_hi('Tom')
#=> prints 'Hi, Tom' to STDOUT.
{% endhighlight %}


{% highlight c++ %}
int aa(int b,char c)
{
return 0;
}
// comment
{% endhighlight %}


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
Nucl.Phys. B896 (2015) 19-51.*

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

## Contacts


### License

[GPLv3 License](https://www.gnu.org/licenses/gpl.html)

<!-- <div class="github-fork-ribbon-wrapper right fixed" style="width: 150px;height: 150px;position: fixed;overflow: hidden;top: 0;z-index: 9999;pointer-events: none;right: 0;"><div class="github-fork-ribbon" style="position: absolute;padding: 2px 0;background-color: #333;background-image: linear-gradient(to bottom, rgba(0, 0, 0, 0), rgba(0, 0, 0, 0.15));-webkit-box-shadow: 0 2px 3px 0 rgba(0, 0, 0, 0.5);-moz-box-shadow: 0 2px 3px 0 rgba(0, 0, 0, 0.5);box-shadow: 0 2px 3px 0 rgba(0, 0, 0, 0.5);z-index: 9999;pointer-events: auto;top: 42px;right: -43px;-webkit-transform: rotate(45deg);-moz-transform: rotate(45deg);-ms-transform: rotate(45deg);-o-transform: rotate(45deg);transform: rotate(45deg);"><a href="https://github.com/chibicode/solo" style="font: 700 13px &quot;Helvetica Neue&quot;, Helvetica, Arial, sans-serif;color: #fff;text-decoration: none;text-shadow: 0 -1px rgba(0, 0, 0, 0.5);text-align: center;width: 200px;line-height: 20px;display: inline-block;padding: 2px 0;border-width: 1px 0;border-style: dotted;border-color: rgba(255, 255, 255, 0.7);">Fork me on GitHub</a></div></div> -->
