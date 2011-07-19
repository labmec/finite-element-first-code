/*
 *  Test_Elemento_2d.cpp
 *  Neo_Fem_sc
 *
 *  Created by Nathan Shauer on 7/18/11.
 *  Copyright 2011 Unicamp. All rights reserved.
 *
 */


#include "telemento2d_quad.h"
#include "pzfmatrix.h"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_SUITE(Elementoquad_2d_test)

BOOST_AUTO_TEST_CASE(shape_functions)
{
	std::vector <int> nodes(4,0);
	telemento2d_quad obj(1,2,nodes);
	int Nshape=(2+1)*(2+1);
	int dim=2;
	std::vector<double> phi(Nshape), knownphi(Nshape,0), point(2,0);
	TPZFMatrix dphi(dim,Nshape);
	
	//Teste 1 - lado esquerdo em baixo
	point[0] = -1;
	point[1] = -1;
	obj.Shape(point, phi, dphi);
	knownphi[0] = 1;
	for (int is = 0 ; is < Nshape ; is++) 
	{
		BOOST_CHECK_EQUAL(phi[is],knownphi[is]);
	}
	
	//Teste 2 - lado direito em baixo
	point[0] = 1;
	point[1] = -1;
	obj.Shape(point, phi, dphi);
	knownphi[0] = 0;
	knownphi[2] = 1;
	for (int is = 0 ; is < Nshape ; is++) 
	{
		BOOST_CHECK_EQUAL(phi[is],knownphi[is]);
	}
	
	//Teste 3 - lado esquerdo em cima
	point[0] = -1;
	point[1] = 1;
	obj.Shape(point, phi, dphi);
	knownphi[2] = 0;
	knownphi[6] = 1;
	for (int is = 0 ; is < Nshape ; is++) 
	{
		BOOST_CHECK_EQUAL(phi[is],knownphi[is]);
	}
	
	//Teste 4 - lado direito em cima
	point[0] = 1;
	point[1] = 1;
	obj.Shape(point, phi, dphi);
	knownphi[6] = 0;
	knownphi[8] = 1;
	for (int is = 0 ; is < Nshape ; is++) 
	{
		BOOST_CHECK_EQUAL(phi[is],knownphi[is]);
	}
	
}

/*
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
	
	TIntegracao2d_quad intrule(6);
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
*/

BOOST_AUTO_TEST_SUITE_END()
