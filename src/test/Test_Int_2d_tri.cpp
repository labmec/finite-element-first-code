//
//  TestTopology.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/6/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//


//#include <gmm/gmm.h>


#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "tintegracao2d_tri.h"


double func2(double x, double y);

BOOST_AUTO_TEST_SUITE(int_tri_2d_test)


BOOST_AUTO_TEST_CASE(known_result)
{
	
    TIntegracao2d_tri intrule(4);
    int npoints = intrule.NPoints();
    double weight, result = 0;
    std::pair <double,double> x;
    for (int i = 0; i < npoints; i++) 
    {
        intrule.Point(i, x, weight);
        result += func2(x.first, x.second)*weight;
    }
    double exact = 1.;
    
	
    BOOST_CHECK_SMALL(result-exact,1E-3);
}


BOOST_AUTO_TEST_SUITE_END()


double func2(double x, double y)
{
  return (2*x + 3*y + 5*x*x*x*x);
}