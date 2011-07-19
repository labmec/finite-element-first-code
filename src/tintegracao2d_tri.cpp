//
//  TIntegracao2d_tri.cpp
//  Fem_sc
//
//  Created by Tiago Costa on 09/05/11.
//  Copyright 2011 Unicamp. All rights reserved.
//

#include "TIntegracao2d_tri.h"
#include <iostream>


//**************Integracao pela Quadratura de Gauss-Legendre n-point***************


TIntegracao2d_tri::TIntegracao2d_tri(int order)
{
    fOrder = order;
    int np, xmin = -1., xmax = 1.;
    switch(order)// está correto: order*order?
    {
        case 1:
            np = 1;
            fPontos.resize(np);
            fPesos.resize(np);
            fPontos[0].first  = 0.333333333333333;
            fPontos[0].second = 0.333333333333333;
            fPesos[0] = 1.;
            break;
        case 2:
            np = 3;
            fPontos.resize(np);
            fPesos.resize(np);
            fPontos[0].first  = 0.166666666666667;
            fPontos[0].second = 0.166666666666667;
            fPesos[0]         = 0.333333333333333;
            fPontos[1].first  = 0.666666666666667;
            fPontos[1].second = 0.166666666666667;
            fPesos[1]         = 0.333333333333333;
            fPontos[2].first  = 0.166666666666667;
            fPontos[2].second = 0.666666666666667;
            fPesos[2]         = 0.333333333333333;
            break;
        case 3:
            np = 4;
            fPontos.resize(np);
            fPesos.resize(np);
            fPontos[0].first  = 0.2;
            fPontos[0].second = 0.2;
            fPesos[0]         = 0.520833333333333;
            fPontos[1].first  = 0.6;
            fPontos[1].second = 0.2;
            fPesos[1]         = 0.520833333333333;
            fPontos[2].first  = 0.2;
            fPontos[2].second = 0.6;
            fPesos[2]         = 0.520833333333333;
            fPontos[3].first  = 0.333333333333333;
            fPontos[3].second = 0.333333333333333;
            fPesos[3]         =-0.56250;
            break;
        case 4:
            np = 6;
            fPontos.resize(np);
            fPesos.resize(np);
            fPontos[0].first  = 0.091576213509771;
            fPontos[0].second = 0.091576213509771;
            fPesos[0]         = 0.109951743655322;
            fPontos[1].first  = 0.816847572980459;
            fPontos[1].second = 0.091576213509771;
            fPesos[1]         = 0.109951743655322;
            fPontos[2].first  = 0.091576213509771;
            fPontos[2].second = 0.816847572980459;
            fPesos[2]         = 0.109951743655322;
            fPontos[3].first  = 0.445948490915965;
            fPontos[3].second = 0.108103018168070;
            fPesos[3]         = 0.223381589678011;
            fPontos[4].first  = 0.445948490915965;
            fPontos[4].second = 0.445948490915965;
            fPesos[4]         = 0.223381589678011;
            fPontos[5].first  = 0.108103018168070;
            fPontos[5].second = 0.445948490915965;
            fPesos[5]         = 0.223381589678011;
            break;
        case 5:
            np = 7;
            fPontos.resize(np);
            fPesos.resize(np);
            fPontos[0].first  = 0.101286507323456;
            fPontos[0].second = 0.101286507323456;
            fPesos[0]         = 0.125939180544827;
            fPontos[1].first  = 0.797426985353087;
            fPontos[1].second = 0.101286507323456;
            fPesos[1]         = 0.125939180544827;
            fPontos[2].first  = 0.101286507323456;
            fPontos[2].second = 0.797426985353087;
            fPesos[2]         = 0.125939180544827;
            fPontos[3].first  = 0.470142064105115;
            fPontos[3].second = 0.059715871789770;
            fPesos[3]         = 0.132394152788506;
            fPontos[4].first  = 0.470142064105115;
            fPontos[4].second = 0.470142064105115;
            fPesos[4]         = 0.132394152788506;
            fPontos[5].first  = 0.059715871789770;
            fPontos[5].second = 0.470142064105115;
            fPesos[5]         = 0.132394152788506;
            fPontos[6].first  = 0.333333333333333;
            fPontos[6].second = 0.333333333333333;
            fPesos[6]         = 0.225030000300000;
            break;
        case 6:
            np = 11;
            fPontos.resize(np);
            fPesos.resize(np);
            fPontos[0].first   = 0.063089014491502;
            fPontos[0].second  = 0.063089014491502;
            fPesos[0]          = 0.050844960370207;
            fPontos[1].first   = 0.873821971016996;
            fPontos[1].second  = 0.063089014491502;
            fPesos[1]          = 0.050844960370207;
            fPontos[2].first   = 0.063089014491502;
            fPontos[2].second  = 0.873821971016996;
            fPesos[2]          = 0.050844960370207;
            fPontos[3].first   = 0.249286745170910;
            fPontos[3].second  = 0.249286745170910;
            fPesos[3]          = 0.116786275726379;
            fPontos[4].first   = 0.501426509658179;
            fPontos[4].second  = 0.249286745170910;
            fPesos[4]          = 0.116786275726379;
            fPontos[5].first   = 0.249286745170910;
            fPontos[5].second  = 0.501426509658179;
            fPesos[5]          = 0.116786275726379;
            fPontos[6].first   = 0.310352451033785;
            fPontos[6].second  = 0.053145049844816;
            fPesos[6]          = 0.082851075618374;
            fPontos[7].first   = 0.636502499121399;
            fPontos[7].second  = 0.053145049844816;
            fPesos[7]          = 0.082851075618374;
            fPontos[8].first   = 0.636502499121399;
            fPontos[8].second  = 0.310352451033785;
            fPesos[8]          = 0.082851075618374;
            fPontos[9].first   = 0.310352451033785;
            fPontos[9].second  = 0.636502499121399;
            fPesos[9]          = 0.082851075618374;
            fPontos[10].first  = 0.053145049844816;
            fPontos[10].second = 0.636502499121399;
            fPesos[10]         = 0.082851075618374;
            fPontos[11].first  = 0.053145049844816;
            fPontos[11].second = 0.310352451033785;
            fPesos[11]         = 0.082851075618374;
            break;
        }
   }

int TIntegracao2d_tri::NPoints()
{
    return fPontos.size();    
}

void TIntegracao2d_tri::Point(int ip, std::pair<double,double> &x, double &weight)
{
    if (ip < 0 || ip >= fPontos.size()) 
    {
        femsc_exception toto;
        throw toto;
    }
    x = fPontos[ip];
    weight = fPesos[ip]/2.; //Pois a integracao numerica de triangulos pede a multiplicacao do somatorio final por 1/2.
}

