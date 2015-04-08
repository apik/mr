//
// MR - 2-loop matching and 3-loop Running, including full 2-loop EW corrections
// Copyright (C) 2014 Andrey Pikelner <pikelner@theor.jinr.ru>
//
// This file is part of MR.
//
// MR is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MR is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MR.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef __TOOLS_HPP__
#define __TOOLS_HPP__

#include "sminput.hpp"

#include <Eigen/Dense>

#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

// Generic functor
template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{
  typedef _Scalar Scalar;
  enum {
    InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
  };
  typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

  int m_inputs, m_values;

  Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
  Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }

};

// struct my_functor : Functor<double>
// {
//   my_functor(void): Functor<double>(2,2) {}
//   int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
//   {
//     // Implement y = 10*(x0+3)^2 + (x1-5)^2
//     fvec(0) = 10.0*pow(x(0)+3.0,4) +  pow(x(1)-5.0,4);
//     fvec(1) = 0;

//     return 0;
//   }
// };


long double critMH(const OSinput&, long double);

long double critMt(const OSinput&, long double);


long double critMH2(const OSinput&, long double);

long double critMH_scaleNotFixed(const OSinput&, long double);

long double critMt_scaleNotFixed(const OSinput&, long double);

#endif  // __TOOLS_HPP__
