//
//  TestTopology.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/6/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//


//#include <gmm/gmm.h>
#include "tintegracao2d_quad.h"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

double func(double x, double y);

BOOST_AUTO_TEST_SUITE(int_quad_2d_test)

BOOST_AUTO_TEST_CASE(test_npoints)
{
	
	TIntegracao2d_quad intrule(1);
	int npoints = intrule.NPoints();
	
	BOOST_CHECK_EQUAL(npoints,1);
	
	TIntegracao2d_quad intrule2(2);
	npoints = intrule2.NPoints();
	
	BOOST_CHECK_EQUAL(npoints,4);
	
	TIntegracao2d_quad intrule3(4);
	npoints = intrule3.NPoints();
	
	BOOST_CHECK_EQUAL(npoints,9);
}

BOOST_AUTO_TEST_CASE(known_result)
{
	
		TIntegracao2d_quad intrule(4);
		int npoints = intrule.NPoints();
		double weight, result = 0;
		std::pair <double,double> x;
		for (int i = 0; i < npoints; i++) 
		{
			intrule.Point(i, x, weight);
			result += func(x.first, x.second)*weight;
		}
		double exact = 2.197;
		
	
		BOOST_CHECK_SMALL(result-exact,1E-3);
}

// soma dos pesos deve ser 4
BOOST_AUTO_TEST_CASE(weightsum)
{
	TIntegracao2d_quad intrule(1);
	int npoints = intrule.NPoints();
	double weight, sum = 0;
	std::pair <double,double> x;
	for (int i = 0 ; i < npoints ; i++) 
	{
		intrule.Point(i, x, weight);
		sum += weight;
	}
	BOOST_CHECK_SMALL(sum-4, 1E-8);
	
	TIntegracao2d_quad intrule2(2);
	npoints = intrule2.NPoints();
	sum = 0;
	for (int i = 0 ; i < npoints ; i++) 
	{
		intrule2.Point(i, x, weight);
		sum += weight;
	}
	BOOST_CHECK_SMALL(sum-4, 1E-8);
	
	TIntegracao2d_quad intrule3(4);
	npoints = intrule3.NPoints();
	sum = 0;
	for (int i = 0 ; i < npoints ; i++) 
	{
		intrule3.Point(i, x, weight);
		sum += weight;
	}
	BOOST_CHECK_SMALL(sum-4, 1E-8);
	
}

BOOST_AUTO_TEST_CASE(pontoforadoescopo)
{
	TIntegracao2d_quad unid(3);
	std::pair<double,double> ponto(0,0);
	double peso;
	int np = unid.NPoints();
	BOOST_CHECK_THROW(unid.Point(np,ponto,peso), femsc_exception);
	BOOST_CHECK_THROW(unid.Point(-1,ponto,peso), femsc_exception);
}


BOOST_AUTO_TEST_SUITE_END()

double func(double x, double y)
{
  return (1+x)/(2+y);
}