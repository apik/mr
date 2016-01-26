// Example of coefficients in expressions 
// for conversion OS mass to MS and 
// MS couplings in terms of OS masses

#include "mr.hpp"

using namespace mr;

int main (int argc, char *argv[])
{
  try
    {
      // loglevel = logINFO;
      
      // Input: Pole masses and Fermi constant in OS scheme
      OSinput oi(pdg2014::Mb, pdg2014::MW, pdg2014::MZ, pdg2014::MH, pdg2014::Mt);
      
      bb<OS> Xbb(oi, oi.MMt());
      WW<OS> XWW(oi, oi.MMt());
      ZZ<OS> XZZ(oi, oi.MMt());
      HH<OS> XHH(oi, oi.MMt());
      tt<OS> Xtt(oi, oi.MMt());





      const size_t fw = 12;
      std::cout << std::setprecision(2) << std::fixed;

      // OS masses to MS constants
    
      std::cout << "MS parameters in terms of OS" << std::endl;
      std::cout << " [ Y_ij OS masses to MS constants, eq.36 from 1503.02138 ]" << std::endl;
      
      std::cout << std::endl << " W-boson " << std::endl;
      
      std::cout << "     as     " << "     al     " << "    al*as   " 
                << "    al^2    " << "    as^2    " << "    as^3    " << std::endl;
      std::cout << std::setw(fw) << "-" 
                << std::setw(fw) << XWW.y10() 
                << std::setw(fw) << XWW.y11() 
                << std::setw(fw) << XWW.y20() 
                << std::setw(fw) << "-"
                << std::setw(fw) << "-"
                << std::endl;


      std::cout << std::endl << " Z-boson " << std::endl;

      std::cout << "     as     " << "     al     " << "    al*as   " 
                << "    al^2    " << "    as^2    " << "    as^3    " << std::endl;
      std::cout << std::setw(fw) << "-" 
                << std::setw(fw) << XZZ.y10() 
                << std::setw(fw) << XZZ.y11() 
                << std::setw(fw) << XZZ.y20() 
                << std::setw(fw) << "-"
                << std::setw(fw) << "-"
                << std::endl;
      

      std::cout << std::endl << " H-boson " << std::endl;

      std::cout << "     as     " << "     al     " << "    al*as   " 
                << "    al^2    " << "    as^2    " << "    as^3    " << std::endl;
      std::cout << std::setw(fw) << "-" 
                << std::setw(fw) << XHH.y10() 
                << std::setw(fw) << XHH.y11() 
                << std::setw(fw) << XHH.y20() 
                << std::setw(fw) << "-"
                << std::setw(fw) << "-"
                << std::endl;


      std::cout << std::endl << " b-quark " << std::endl;

      std::cout << "     as     " << "     al     " << "    al*as   " 
                << "    al^2    " << "    as^2    " << "    as^3    " << std::endl;
      std::cout << std::setw(fw) << Xbb.y01() 
                << std::setw(fw) << Xbb.y10() 
                << std::setw(fw) << Xbb.y11() 
                << std::setw(fw) << Xbb.y20() 
                << std::setw(fw) << Xbb.y02() 
                << std::setw(fw) << Xbb.y03() << std::endl;


      std::cout << std::endl << " t-quark " << std::endl;

      std::cout << "     as     " << "     al     " << "    al*as   " 
                << "    al^2    " << "    as^2    " << "    as^3    " << std::endl;
      std::cout << std::setw(fw) << Xtt.y01() 
                << std::setw(fw) << Xtt.y10() 
                << std::setw(fw) << Xtt.y11() 
                << std::setw(fw) << Xtt.y20() 
                << std::setw(fw) << Xtt.y02() 
                << std::setw(fw) << Xtt.y03() << std::endl;

      
      // OS masses to MS masses
      
      std::cout << std::endl << std::endl << " [ X_ij OS masses to MS for bosons, eq.39 from 1503.02138 ]" << std::endl;

      std::cout << std::endl << " W-boson " << std::endl;
      
      std::cout << "     as     " << "     al     " << "    al*as   " 
                << "    al^2    " << "    as^2    " << "    as^3    " << std::endl;
      std::cout << std::setw(fw) << "-" 
                << std::setw(fw) << XWW.x10() 
                << std::setw(fw) << XWW.x11() 
                << std::setw(fw) << XWW.x20() 
                << std::setw(fw) << "-"
                << std::setw(fw) << "-"
                << std::endl;


      std::cout << std::endl << " Z-boson " << std::endl;

      std::cout << "     as     " << "     al     " << "    al*as   " 
                << "    al^2    " << "    as^2    " << "    as^3    " << std::endl;
      std::cout << std::setw(fw) << "-" 
                << std::setw(fw) << XZZ.x10() 
                << std::setw(fw) << XZZ.x11() 
                << std::setw(fw) << XZZ.x20() 
                << std::setw(fw) << "-"
                << std::setw(fw) << "-"
                << std::endl;
      

      std::cout << std::endl << " H-boson " << std::endl;

      std::cout << "     as     " << "     al     " << "    al*as   " 
                << "    al^2    " << "    as^2    " << "    as^3    " << std::endl;
      std::cout << std::setw(fw) << "-" 
                << std::setw(fw) << XHH.x10() 
                << std::setw(fw) << XHH.x11() 
                << std::setw(fw) << XHH.x20() 
                << std::setw(fw) << "-"
                << std::setw(fw) << "-"
                << std::endl;

      std::cout << " \n[ X_ij OS masses to MS for fermions, eq.40 from 1503.02138 ]" << std::endl;

      std::cout << " b-quark " << std::endl;

      std::cout << "     as     " << "     al     " << "    al*as   " 
                << "    al^2    " << "    as^2    " << "    as^3    " << std::endl;
      std::cout << std::setw(fw) << Xbb.x01() 
                << std::setw(fw) << Xbb.x10() 
                << std::setw(fw) << Xbb.x11() 
                << std::setw(fw) << Xbb.x20() 
                << std::setw(fw) << Xbb.x02() 
                << std::setw(fw) << Xbb.x03() << std::endl;


      std::cout << std::endl << " t-quark " << std::endl;

      std::cout << "     as     " << "     al     " << "    al*as   " 
                << "    al^2    " << "    as^2    " << "    as^3    " << std::endl;
      std::cout << std::setw(fw) << Xtt.x01() 
                << std::setw(fw) << Xtt.x10() 
                << std::setw(fw) << Xtt.x11() 
                << std::setw(fw) << Xtt.x20() 
                << std::setw(fw) << Xtt.x02() 
                << std::setw(fw) << Xtt.x03() << std::endl;








      // MS mass to OS, inverse relations

      AlphaS as(oi);

      // From the same OS input ass in first part obtain MS masses 
      P2MS<AlphaGF> pMSmt(oi,pdg2014::Gf, as(oi.Mt()), oi.Mt(), order::all);
    
    
      MSinput mi(pMSmt.getMSpar());
      
      bb<MS> xbb(mi, oi.MMt());
      WW<MS> xWW(mi, oi.MMt());
      ZZ<MS> xZZ(mi, oi.MMt());
      HH<MS> xHH(mi, oi.MMt());
      tt<MS> xtt(mi, oi.MMt());


      std::cout << std::setprecision(2) << std::fixed;
      
      std::cout << "\n\n\nInverse relations for OS mass in terms of running MS masses" << std::endl;
      std::cout << " [ x_ij MS masses to OS for bosons ]" << std::endl;
      
      std::cout << std::endl << " W-boson " << std::endl;
      
      std::cout << "     as     " << "     al     " << "    al*as   " 
                << "    al^2    " << "    as^2    " << "    as^3    " << std::endl;
      std::cout << std::setw(fw) << "-" 
                << std::setw(fw) << xWW.x10() 
                << std::setw(fw) << xWW.x11() 
                << std::setw(fw) << xWW.x20() 
                << std::setw(fw) << "-"
                << std::setw(fw) << "-"
                << std::endl;


      std::cout << std::endl << " Z-boson " << std::endl;

      std::cout << "     as     " << "     al     " << "    al*as   " 
                << "    al^2    " << "    as^2    " << "    as^3    " << std::endl;
      std::cout << std::setw(fw) << "-" 
                << std::setw(fw) << xZZ.x10() 
                << std::setw(fw) << xZZ.x11() 
                << std::setw(fw) << xZZ.x20() 
                << std::setw(fw) << "-"
                << std::setw(fw) << "-"
                << std::endl;
      

      std::cout << std::endl << " H-boson " << std::endl;

      std::cout << "     as     " << "     al     " << "    al*as   " 
                << "    al^2    " << "    as^2    " << "    as^3    " << std::endl;
      std::cout << std::setw(fw) << "-" 
                << std::setw(fw) << xHH.x10() 
                << std::setw(fw) << xHH.x11() 
                << std::setw(fw) << xHH.x20() 
                << std::setw(fw) << "-"
                << std::setw(fw) << "-"
                << std::endl;

      std::cout << " \n[ x_ij MS masses to OS for fermions ]" << std::endl;

      std::cout << " b-quark " << std::endl;

      std::cout << "     as     " << "     al     " << "    al*as   " 
                << "    al^2    " << "    as^2    " << "    as^3    " << std::endl;
      std::cout << std::setw(fw) << xbb.x01() 
                << std::setw(fw) << xbb.x10() 
                << std::setw(fw) << xbb.x11() 
                << std::setw(fw) << xbb.x20() 
                << std::setw(fw) << xbb.x02() 
                << std::setw(fw) << xbb.x03() << std::endl;


      std::cout << std::endl << " t-quark " << std::endl;

      std::cout << "     as     " << "     al     " << "    al*as   " 
                << "    al^2    " << "    as^2    " << "    as^3    " << std::endl;
      std::cout << std::setw(fw) << xtt.x01() 
                << std::setw(fw) << xtt.x10() 
                << std::setw(fw) << xtt.x11() 
                << std::setw(fw) << xtt.x20() 
                << std::setw(fw) << xtt.x02() 
                << std::setw(fw) << xtt.x03() << std::endl;

      

    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

