#include <algorithm>
#include <vector>

typedef std::vector<long double> state_type;

// template < int pocoa1, int pocoa2, int pocoas, int pocoat, int pocoab, int pocoatau, int pocolam > 
class BetaSMFull
{
  // Maximal numer of loops available for beta-functions
  static const size_t maxBetaOrder = 3;

  static const size_t maxPoco = maxBetaOrder + 2;

  long double be1[maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2];
  long double be2[maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2];
  long double be3[maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2];
  long double be4[maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2];
  long double be5[maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2];
  long double be6[maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2];
  long double be7[maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2][maxBetaOrder+2];
  
  int pocoa1, pocoa2, pocoas, pocoat, pocoab, pocoatau, pocolam;
  size_t maxPower;
public:

  BetaSMFull(int pocoa1_, int pocoa2_, int pocoas_, int pocoat_, int pocoab_, int pocoatau_, int pocolam_);

  void operator() (const state_type &, state_type &, const double);
  long double betaQCD(long double);
};

// template<size_t pocoGauge=3, size_t pocoYukawa=3, size_t pocoLam=3> 
// class BetaSM : public BetaSMFull<pocoGauge, pocoGauge, pocoGauge, pocoYukawa, pocoYukawa, pocoYukawa, pocoLam>
// {
//   // public:
//   //   template <class T2> CDerived(const T2 &x) : CBase<T>(x) {;}
//   //   template <class T2> CDerived (const CBase<T2> &x) : CBase<T>(x) {;}
//   //   ~CDerived() {;}
// };


