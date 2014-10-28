#ifndef __BASE_HPP__
#define __BASE_HPP__

class PoleMass
{
    // 
  // Pure QCD part m_ij=mY_ij by definition
  // 
  // virtual std::complex<long double> m01();

  // virtual std::complex<long double> m02(size_t nL = 5);

  // virtual std::complex<long double> m03(size_t nL = 5);

  // 
  // Mass corrections
  // 
public:
  virtual std::complex<long double> m10(size_t nL = 2, size_t nH = 1, size_t boson = 1) = 0;
  
  virtual std::complex<long double> m11(size_t nL = 2, size_t nH = 1, size_t boson = 1) = 0;

  virtual std::complex<long double> m20(size_t nL = 2, size_t nH = 1, size_t boson = 1) = 0;


  // Gaugeless limit
  // virtual std::complex<long double> mgl01(size_t nL = 2, size_t nH = 1);

  // virtual std::complex<long double> mgl10(size_t nL = 2, size_t nH = 1);
  
  // virtual std::complex<long double> mgl11(size_t nL = 2, size_t nH = 1);

  // virtual std::complex<long double> mgl20(size_t nL = 2, size_t nH = 1);
  

  // 
  // mass definition using Yukawa couplings 
  // mY=y/sqrt(2*sqrt(2)*GF)
  // 
  // virtual std::complex<long double> my01() = 0;

  virtual std::complex<long double> my10(size_t nL = 2, size_t nH = 1, size_t boson = 1) = 0;

  virtual std::complex<long double> my11(size_t nL = 2, size_t nH = 1, size_t boson = 1) = 0;

  virtual std::complex<long double> my20(size_t nL = 2, size_t nH = 1, size_t boson = 1) = 0;


  // Gaugeless limit
  // virtual std::complex<long double> mygl01();

  // virtual std::complex<long double> mygl10(size_t nL = 2, size_t nH = 1);

  // virtual std::complex<long double> mygl11(size_t nL = 2, size_t nH = 1);

  // virtual std::complex<long double> mygl20(size_t nL = 2, size_t nH = 1);


};

// class RunningMass
// {
//   virtual void  m01()
// };

#endif  // __BASE_HPP__
