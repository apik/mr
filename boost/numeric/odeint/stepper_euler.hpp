/* Boost odeint/stepper_euler.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 
 This file includes the explicit euler solver for
 ordinary differential equations.

 It solves any ODE dx/dt = f(x,t) via
 x(t+dt) = x(t) + dt*f(x,t)

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_EULER_HPP
#define BOOST_NUMERIC_ODEINT_STEPPER_EULER_HPP

#include <boost/numeric/odeint/detail/iterator_algebra.hpp>
#include <boost/numeric/odeint/container_traits.hpp>

namespace boost {
namespace numeric {
namespace odeint {

    template<
        class Container ,
        class Time = double ,
        class Traits = container_traits< Container >
        >
    class stepper_euler
    {
        //
        // provide basic typedefs
        //
    public:

        typedef unsigned short order_type;
        typedef Time time_type;
        typedef Traits traits_type;
        typedef typename traits_type::container_type container_type;
        typedef container_type state_type;
        typedef typename traits_type::value_type value_type;


        //
        // private members
        //
    private:

        state_type m_dxdt;



        //
        // public interface
        //
    public:

        // return the order of the stepper
        order_type order_step() const { return 1; }


        // standard constructor, m_dxdt is not adjusted
        stepper_euler( void )
        {
        }



        // contructor, which adjusts m_dxdt
        stepper_euler( const state_type &x )
        {
            adjust_size( x );
        }



        // adjust the size of m_dxdt
        void adjust_size( const state_type &x )
        {
            traits_type::adjust_size( x , m_dxdt );
        }



        // performs one step with the knowledge of dxdt(t)
        template< class DynamicalSystem >
        void do_step( DynamicalSystem &system ,
		      state_type &x ,
		      const state_type &dxdt ,
		      time_type t ,
		      time_type dt )
        {
            // x = x + dt*dxdt
            detail::it_algebra::increment( traits_type::begin(x) ,
                                           traits_type::end(x) ,
                                           traits_type::begin(dxdt) , 
                                           dt );
        }



        // performs one step
        template< class DynamicalSystem >
        void do_step( DynamicalSystem &system ,
                      state_type &x ,
                      time_type t ,
                      time_type dt )
        {
        	system( x , m_dxdt , t );
            do_step( system , x , m_dxdt , t , dt );
        }

    };



} // namespace odeint
} // namespace numeric
} // namespace boost


#endif // BOOST_NUMERIC_ODEINT_STEPPER_EULER_HPP
