//
//  tintegracao.cpp
//  Fem_sc
//
//  Created by Tiago Costa on 09/05/11.
//  Copyright 2011 Unicamp. All rights reserved.
//

#include "tintegracao.h"

#include "nr3.h"
#include "gauss_wgts.h"
#include "gauss_wgts2.h"

//**************Integracao pela Quadratura de Gauss-Legendre n-point***************

TIntegracao::TIntegracao(int order, double xmin, double xmax)
{
    VecDoub x1;
    fOrder = order;
    int np;
    if (order%2==1)
    {
        np = (order+1)/2;   
    }
    else 
    {
        np = (order/2+1);
    }
    fPontos.resize(np);
    fPesos.resize(np);
    gauleg(xmin,xmax,fPontos,fPesos);
}

int TIntegracao::NPoints()
{
    return fPontos.size();
    
}

void TIntegracao::Point(int ip, double &x, double &weight)
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

