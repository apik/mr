#include <tt.hpp>
std::complex<long double> tt::mgl11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mttgl[14];

    mttgl[1]=pow(SW,-1);
    mttgl[2]=pow(MMH,-1);
    mttgl[3]=pow(MMW,-1);
    mttgl[4]=Tsil::B(MMH,MMt,MMt,mu2);
    mttgl[5]=Tsil::A(MMt,mu2);
    mttgl[6]=Tsil::A(MMH,mu2);
    mttgl[7]=pow(MMt,-1);
    mttgl[8]=std::real(Tsil::B(0,0,MMt,mu2));
   mttgl[9]=3*mttgl[5];
   mttgl[10]= - MMt - mttgl[9];
   mttgl[9]=mttgl[10]*mttgl[9];
   mttgl[10]=pow(MMt,2);
   mttgl[9]=mttgl[10] + mttgl[9];
   mttgl[9]=mttgl[2]*mttgl[9];
   mttgl[10]=1./3.*MMt;
   mttgl[11]=2*mttgl[5];
   mttgl[12]=mttgl[10] + mttgl[11];
   mttgl[12]=mttgl[8]*mttgl[12];
   mttgl[13]= - MMH + 2*MMt;
   mttgl[13]=mttgl[4]*mttgl[13];
   mttgl[11]=mttgl[11] + 3*mttgl[6];
   mttgl[11]=mttgl[7]*mttgl[11];
   mttgl[11]=1./3. + 6*mttgl[4] + mttgl[11];
   mttgl[11]=mttgl[5]*mttgl[11];
   mttgl[9]=mttgl[12] + 8*mttgl[9] + mttgl[11] + mttgl[6] + mttgl[10]
    + mttgl[13];

      return mttgl[9]*mttgl[3]*pow(mttgl[1],2);
}
