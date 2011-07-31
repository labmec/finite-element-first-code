/***************************************************************************
 *   Copyright (C) 2005 by Philippe R. B. Devloo                           *
 *   phil@corona                                                           *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef TELEMENTO1D_H
#define TELEMENTO1D_H

#include <telemento.h>


/**
@author Philippe R. B. Devloo
*/
class TElemento1d : public TElemento
{
public:
    TElemento1d(int matid, int order, std::vector< int >& nodes);
    
    TElemento1d();

    ~TElemento1d();

    virtual MElementType getElType();
    virtual void CalcStiff(TMalha &malha, TPZFMatrix& stiff, TPZFMatrix& rhs);
    virtual void Jacobian(std::vector< double >& point, TPZFMatrix& jacobian, TPZFMatrix& jacinv, double &detjac, TMalha& malha);
	//virtual void RealCoord(std::vector <double> &point);

    
    /**
     * Calcula os valores das funcoes de forma e suas derivadas
     * @point ponto onde calcular as funcoes de forma
     * @phi valores das funcoes de forma
     * @dphi valores das derivadas das funcoes de forma
     */
    virtual void Shape(std::vector<double> &point, std::vector<double> &phi, TPZFMatrix &dphi);
    
    /**
     * Calcula o erro do elemento
     * @param exact funcao que calcula o valor exato
     * @param energy [out] erro na norma de energia
     * @param l2 [out] erro na norma l2
     */
virtual void Error(TPZFMatrix &solution, TMalha &malha, void (f)(std::vector<double> &,double &, std::vector<double> &), double &energy, double &l2);

	virtual void RealCoord(std::vector <double> &point, TMalha &malha);
	


};

#endif
