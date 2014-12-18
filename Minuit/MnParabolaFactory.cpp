#include "Minuit/MnParabolaFactory.h"
#include "Minuit/MnParabola.h"
#include "Minuit/MnParabolaPoint.h"

// #include <iostream>

MnParabola MnParabolaFactory::operator()(const MnParabolaPoint& p1, 
					 const MnParabolaPoint& p2, 
					 const MnParabolaPoint& p3) const {
  double x1 = p1.x();
  double x2 = p2.x();
  double x3 = p3.x();
  double dx12 = x1-x2;
  double dx13 = x1-x3;
  double dx23 = x2-x3;

//   std::cout<<"MnParabolaFactory x1, x2, x3: "<<x1<<" "<<x2<<" "<<x3<<std::endl;

  double xm = (x1+x2+x3)/3.;
  x1 -= xm;
  x2 -= xm;
  x3 -= xm;
  
  double y1 = p1.y();
  double y2 = p2.y();
  double y3 = p3.y();
//   std::cout<<"MnParabolaFactory y1, y2, y3: "<<y1<<" "<<y2<<" "<<y3<<std::endl;

  double a = y1/(dx12*dx13) - y2/(dx12*dx23) + y3/(dx13*dx23);
  double b = -y1*(x2+x3)/(dx12*dx13) + y2*(x1+x3)/(dx12*dx23) - y3*(x1+x2)/(dx13*dx23);
  double c = y1 - a*x1*x1 - b*x1;

  c += xm*(xm*a - b);
  b -= 2.*xm*a;

//   std::cout<<"a,b,c= "<<a<<" "<<b<<" "<<c<<std::endl;
  return MnParabola(a, b, c);
}

MnParabola MnParabolaFactory::operator()(const MnParabolaPoint& p1, 
					 double dxdy1, 
					 const MnParabolaPoint& p2) const {
  
  double x1 = p1.x();
  double xx1 = x1*x1;
  double x2 = p2.x();
  double xx2 = x2*x2;
  double y1 = p1.y();
  double y12 = p1.y() - p2.y();

  double det = xx1-xx2 - 2.*x1*(x1-x2);
  double a = -(        y12 +   (x2-x1)*dxdy1)/det;
  double b = -( -2.*x1*y12 + (xx1-xx2)*dxdy1)/det;
  double c = y1 - a*xx1 - b*x1;

  return MnParabola(a, b, c);
}

