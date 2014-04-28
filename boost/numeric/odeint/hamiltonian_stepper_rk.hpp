/*
 boost header: numeric/odeint/hamiltonian_stepper_rk.hpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#ifndef BOOST_NUMERIC_ODEINT_HAMILTONIAN_STEPPER_RK_HPP_INCLUDED
#define BOOST_NUMERIC_ODEINT_HAMILTONIAN_STEPPER_RK_HPP_INCLUDED

#include <stdexcept>
#include <utility>

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
    class hamiltonian_stepper_rk
    {
        // provide basic typedefs
    public:

        typedef unsigned short order_type;
        typedef Time time_type;
        typedef Traits traits_type;
        typedef typename traits_type::container_type container_type;
        typedef std::pair< container_type , container_type > state_type;
        typedef typename traits_type::value_type value_type;
        typedef typename traits_type::iterator iterator;
        typedef typename traits_type::const_iterator const_iterator;

    private:

        container_type m_dqdt , m_dpdt;

/*
	rk_a[0]=0.40518861839525227722;
	rk_a[1]=-0.28714404081652408900;
	rk_a[2]=0.5-(rk_a[0]+rk_a[1]);
	rk_a[3]=rk_a[2];
	rk_a[4]=rk_a[1];
	rk_a[5]=rk_a[0];

	rk_b[0]=-3.0/73.0;
	rk_b[1]=17.0/59.0;
	rk_b[2]=1.0-2.0*(rk_b[0]+rk_b[1]);
	rk_b[3]=rk_b[1];
	rk_b[4]=rk_b[0];
	rk_b[5]=0.0;
*/

    public:

        template< class SystemFunction >
        void do_step( SystemFunction func , 
                      state_type &state ,
                      time_type t ,
                      time_type dt )
        {
            const size_t order = 6;

            const time_type rk_a[order] = {
                static_cast<time_type>( 0.40518861839525227722 ) ,
                static_cast<time_type>( -0.28714404081652408900 ) ,
                static_cast<time_type>( 0.3819554224212718118 ) ,
                static_cast<time_type>( 0.3819554224212718118 ) ,
                static_cast<time_type>( -0.28714404081652408900 ) ,
                static_cast<time_type>( 0.40518861839525227722 )
            };
            const time_type rk_b[order] = {
                static_cast<time_type>( -3.0/73.0 ) ,
                static_cast<time_type>( 17.0/59.0 ) ,
                static_cast<time_type>( 0.50592059438123984212 ) ,
                static_cast<time_type>( 17.0/59.0 ) ,
                static_cast<time_type>( -3.0/73.0 ) ,
                static_cast<time_type>( 0.0 )
            };
            
            container_type &q = state.first;
            container_type &p = state.second;

            if( !traits_type::same_size( q , p ) )
            {
                std::string msg( "hamiltonian_stepper_rk::do_step(): " );
                msg += "q and p have different sizes";
                throw std::invalid_argument( msg );
            }

            traits_type::adjust_size( p , m_dqdt );
            traits_type::adjust_size( p , m_dpdt );

            for( size_t l=0 ; l<order ; ++l )
            {
                func.first( p , m_dqdt );
                detail::it_algebra::increment( traits_type::begin(q) ,
                                               traits_type::end(q) ,
                                               traits_type::begin(m_dqdt) ,
                                               rk_a[l]*dt );
                func.second( q , m_dpdt );
                detail::it_algebra::increment( traits_type::begin(p) ,
                                               traits_type::end(p) ,
                                               traits_type::begin(m_dpdt) ,
                                               rk_b[l]*dt );
            }
        }
        
    };



    template<
        class Container ,
        class Time = double ,
        class Traits = container_traits< Container >
        >
    class hamiltonian_stepper_rk_qfunc
    {
        // provide basic typedefs
    public:

        typedef unsigned short order_type;
        typedef Time time_type;
        typedef Traits traits_type;
        typedef typename traits_type::container_type container_type;
        typedef std::pair< container_type , container_type > state_type;
        typedef typename traits_type::value_type value_type;
        typedef typename traits_type::iterator iterator;
        typedef typename traits_type::const_iterator const_iterator;

    private:

	container_type m_dpdt;

/*
	rk_a[0]=0.40518861839525227722;
	rk_a[1]=-0.28714404081652408900;
	rk_a[2]=0.5-(rk_a[0]+rk_a[1]);
	rk_a[3]=rk_a[2];
	rk_a[4]=rk_a[1];
	rk_a[5]=rk_a[0];

	rk_b[0]=-3.0/73.0;
	rk_b[1]=17.0/59.0;
	rk_b[2]=1.0-2.0*(rk_b[0]+rk_b[1]);
	rk_b[3]=rk_b[1];
	rk_b[4]=rk_b[0];
	rk_b[5]=0.0;
*/

    public:

        template< class CoordinateFunction >
        void do_step( CoordinateFunction qfunc ,
                      state_type &state ,
                      time_type t ,
                      time_type dt )
        {
            const size_t order = 6;

            const time_type rk_a[order] = {
                static_cast<time_type>( 0.40518861839525227722 ) ,
                static_cast<time_type>( -0.28714404081652408900 ) ,
                static_cast<time_type>( 0.3819554224212718118 ) ,
                static_cast<time_type>( 0.3819554224212718118 ) ,
                static_cast<time_type>( -0.28714404081652408900 ) ,
                static_cast<time_type>( 0.40518861839525227722 )
            };
            const time_type rk_b[order] = {
                static_cast<time_type>( -3.0/73.0 ) ,
                static_cast<time_type>( 17.0/59.0 ) ,
                static_cast<time_type>( 0.50592059438123984212 ) ,
                static_cast<time_type>( 17.0/59.0 ) ,
                static_cast<time_type>( -3.0/73.0 ) ,
                static_cast<time_type>( 0.0 )
            };
            
            container_type &q = state.first;
            container_type &p = state.second;

            if( !traits_type::same_size( q , p ) )
            {
                std::string msg( "hamiltonian_stepper_rk::do_step(): " );
                msg += "q and p have different sizes";
                throw std::invalid_argument( msg );
            }

            traits_type::adjust_size( p , m_dpdt );

            for( size_t l=0 ; l<order ; ++l )
            {
                detail::it_algebra::increment( traits_type::begin(q) ,
                                               traits_type::end(q) ,
                                               traits_type::begin(p) ,
                                               rk_a[l]*dt );
                qfunc( q , m_dpdt );
                detail::it_algebra::increment( traits_type::begin(p) ,
                                               traits_type::end(p) ,
                                               traits_type::begin(m_dpdt) ,
                                               rk_b[l]*dt );
            }
        }
        
    };


} // namespace odeint
} // namespace numeric
} // namespace boost





#endif //BOOST_NUMERIC_ODEINT_HAMILTONIAN_STEPPER_RK_HPP_INCLUDED
