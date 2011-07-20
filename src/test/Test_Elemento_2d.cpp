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
#include "tmalha.h"
#include "tmaterial1d.h"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_SUITE(Elementoquad_2d_test)

//Inacabado
BOOST_AUTO_TEST_CASE(jacobian)
{
	TMalha malha;
	int morder = 2;
	std::vector <int> nodes(9,0);
	for (int i = 0 ; i < 9 ; i++) 
	{
		nodes[i] = i;
	}
	
	std::vector<double> co(2,0);
	for (double i = 0; i <= 1 ; i += 0.5) 
	{
		for (double j = 0 ; j <= 1; j = j + 0.5)
		{
			co[0] = j;
			co[1] = i;
			TNo no(co);
			malha.insertNode(no);
		}
	}
	
	int id = 1, k = 1, c = 0, b = 1, f = 1;
	TMaterial1d *mat = new TMaterial1d(id,k,c,b,f);
	malha.insertMaterial(mat);
	telemento2d_quad *elem = new telemento2d_quad(1,morder,nodes);
	malha.insertElement(elem);
	TPZFMatrix stiff, rhs;
	std::vector <double> ponto(2,0);
	double detjac;
	elem->JacobTest(malha,stiff,rhs,ponto,detjac);
	BOOST_CHECK_EQUAL(detjac,0.25);
}

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

BOOST_AUTO_TEST_SUITE_END()
