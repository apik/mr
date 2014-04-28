/* Boost odeint/stepper_rk4.hpp header file
 
 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 
 This file includes the explicit runge kutta solver for
 ordinary differential equations.

 It solves any ODE dx/dt = f(x,t).

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_STEPPER_RK4_CLASSICAL_HPP
#define BOOST_NUMERIC_ODEINT_STEPPER_RK4_CLASSICAL_HPP

#include <boost/numeric/odeint/container_traits.hpp>
#include <boost/numeric/odeint/detail/iterator_algebra.hpp>

namespace boost {
namespace numeric {
namespace odeint {

    template<
        class Container ,
        class Time = double ,
        class Traits = container_traits< Container >
        >
    class stepper_rk4_classical
    {

        // provide basic typedefs
    public:

        typedef unsigned short order_type;
        typedef Time time_type;
        typedef Traits traits_type;
        typedef typename traits_type::container_type container_type;
        typedef container_type state_type;
        typedef typename traits_type::value_type value_type;
//        typedef typename traits_type::iterator iterator;
//        typedef typename traits_type::const_iterator const_iterator;






        // private members
    private:

        state_type m_dxdt;
        state_type m_dxt;
        state_type m_dxm;
        state_type m_xt;

        



        // public interface
    public:

        order_type order_step() const { return 4; }

	// standard constructor, internal containers are not initialized
	stepper_rk4_classical( void )
	{
	}

	// constructor, which adjusts the internal containers
	stepper_rk4_classical( const state_type &x )
	{
	    adjust_size( x );
	}

	void adjust_size( const state_type &x )
	{
	    traits_type::adjust_size( x , m_dxdt );
	    traits_type::adjust_size( x , m_dxt );
	    traits_type::adjust_size( x , m_dxm );
	    traits_type::adjust_size( x , m_xt );
	}

        template< class DynamicalSystem >
        void do_step( DynamicalSystem &system ,
                      state_type &x ,
                      const state_type &dxdt ,
                      time_type t ,
                      time_type dt )
        {
            using namespace detail::it_algebra;

            const time_type val2 = time_type( 2.0 );


            time_type  dh = time_type( 0.5 ) * dt;
            time_type th = t + dh;

            //m_xt = x + dh*dxdt
            scale_sum( traits_type::begin(m_xt) , traits_type::end(m_xt) ,
                       static_cast< time_type >(1.0) ,
                       traits_type::begin(x) , 
                       dh, 
                       traits_type::begin(dxdt) );


            //m_xt = x + dh*m_dxdt
            system( m_xt , m_dxt , th );
            scale_sum( traits_type::begin(m_xt) , traits_type::end(m_xt) ,
                        static_cast< time_type >(1.0),
                        traits_type::begin(x) , 
                        dh ,
                        traits_type::begin(m_dxt) );


            //m_xt = x + dt*m_dxm ; m_dxm += m_dxt
            system( m_xt , m_dxm , th );
            assign_sum_increment( traits_type::begin(m_xt) , traits_type::end(m_xt) ,
                                  traits_type::begin(x) , traits_type::begin(m_dxm) ,
                                  traits_type::begin(m_dxt) , dt );


            //x = dt/6 * ( m_dxdt + m_dxt + val2*m_dxm )
            system( m_xt , m_dxt , t + dt );
            increment_sum_sum( traits_type::begin(x) , traits_type::end(x) , 
                               traits_type::begin(dxdt) ,
                               traits_type::begin(m_dxt) ,
                               traits_type::begin(m_dxm) ,
                               dt /  time_type( 6.0 ) , val2 );
        }



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


#endif // BOOST_NUMERIC_ODEINT_STEPPER_RK4_CLASSICAL_HPP
