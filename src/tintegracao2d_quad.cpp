//
//  TIntegracao2d_quad.cpp
//  Fem_sc
//
//  Created by Tiago Costa on 09/05/11.
//  Copyright 2011 Unicamp. All rights reserved.
//

#include "TIntegracao2d_quad.h"

#include "nr3.h"
//#include <math.h>
#include <iostream>
//#include "gauss_wgts.h"

//**************Integracao pela Quadratura de Gauss-Legendre n-point***************


// IMPLEMENTAR


static void gauleg(const double x1, const double x2, std::vector <double> &x, std::vector <double> &w)
{
	const double EPS=1.0e-16;
	double z1,z,xm,xl,pp,p3,p2,p1;
	int n=x.size();
	int m=(n+1)/2;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for (int i=0;i<m;i++) {
		z=cos(3.14159265358979*(i+0.75)/(n+0.5)); // cos do math.h eh impreciso, USEI o do Numerical Recipes
		do {
			p1=1.0;
			p2=0.0;
			for (int j=0;j<n;j++) {
				p3=p2;
				p2=p1;
				p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
			}
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (abs(z-z1) > EPS);
		x[i]=xm-xl*z;
		x[n-1-i]=xm+xl*z;
		w[i]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n-1-i]=w[i];
	}
}



TIntegracao2d_quad::TIntegracao2d_quad(int order)
{
    fOrder = order;
		int i, j, np1d, np, xmin = -1., xmax = 1.;
		std::vector <double> tempvec, temppesos;
	
    if (order%2==1)
    {
        np1d = (order+1)/2;
				np = np1d*np1d;
    }
    else 
    {
        np1d = (order/2+1);
				np = np1d*np1d;
    }
		tempvec.resize(np1d);
		temppesos.resize(np1d);
    fPontos.resize(np);
    fPesos.resize(np);

    gauleg(xmin,xmax,tempvec,temppesos);
		for (i = 0 ; i < np ; i = i + np1d)
		{
			for (j = 0 ; j < np1d ; j++)
			{
				fPontos[i+j].first = tempvec[j];
				fPontos[i+j].second = tempvec[i/np1d];
				fPesos[i+j] = temppesos[j]*temppesos[i/np1d];
			}
			
		}
}

int TIntegracao2d_quad::NPoints()
{
    return fPontos.size();
}


void TIntegracao2d_quad::Point(int ip, std::pair<double,double> &x, double &weight)
{
    if (ip < 0 || ip >= fPontos.size()) 
    {
        femsc_exception toto;
        throw toto;
    }
		x = fPontos[ip];
    weight = fPesos[ip];
}
/*
**************Integracao pela Quadratura de Gauss-Hermite****************************************
**************Método adequado para integrar funçoes do tipo: ∫ ∞. −∞ e−x2 f(x) dx ***************
**************Não funcionou bem para Integrar a função de teste! Erro muito grande!**************

TIntegrategauher::TIntegrategauher(int order)
{
    fOrder = order;
    int np = order;
    fPontos.resize(np);
    fPesos.resize(np);
    gauher(fPontos,fPesos);
}
int TIntegrategauher::NPoints()
{
    return fPontos.size();
    
}

double TIntegrategauher::Point(int ip, std::vector<double> &loc, double &weight)
{
    if(loc.size() != 1)
    {
        femsc_exception toto;
        throw toto;
    }
    if (ip<0) 
    {
        femsc_exception toto;
        throw toto;
    }
    loc[0]=fPontos[ip];
    weight=fPesos[ip];
    return 0;
}

//**************Integracao pela Quadratura de Gauss-Lobatto****************************************
//**************São incluídos para a integração os pontos extremos e o ponto medio do intervalo ***************

TIntegratelobatto::TIntegratelobatto(int order)
{
    fOrder = order;
    int np = order;
    fPontos.resize(np);
    fPesos.resize(np);
    lobatto(fa, fb,1,-1,1,fPontos,fPesos);
}
int TIntegratelobatto::NPoints()
{
    return fPontos.size();
    
}

double TIntegratelobatto::Point(int ip, std::vector<double> &loc, double &weight)
{
    if(loc.size() != 1)
    {
        femsc_exception toto;
        throw toto;
    }
    if (ip<0) 
    {
        femsc_exception toto;
        throw toto;
    }
    loc[0]=fPontos[ip];
    weight=fPesos[ip];
    return 0;
}*/

