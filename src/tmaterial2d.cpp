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
#include "tmaterial2d.h"
#include "pzfmatrix.h"

TMaterial2d::TMaterial2d(int id, double K, double C, double B, double F) : TMaterial(id), fK(K), fC(C), fB(B), fF(F)
{
}



TMaterial2d::~TMaterial2d()
{
}


void TMaterial2d::Print(std::ostream& out) const
{
    TMaterial::Print(out);
    out << "Coeficient values K " << fK << " C " << fC << " B " << fB << " F " << fF << std::endl;
}

  /**
   * Calcula o valor da contribuição da equação variacional no ponto dado
   * na matriz de rigidez do elemento e no vetor de carga
   * @param pt [in]: ponto de integração de Gauss
   * @param weight [in]: peso de integração
   * @param phiVal [in] : valor da função teste no ponto dado
   * @param dphi [in] : valor das derivadas da função de forma no ponto de integração
   * @param elementK [inout]: matriz de rigidez do elemento
   * @param elementF [inout]: vetor de carga do elemento
   */
void TMaterial2d::Contribute (std::vector<double> &point, double  weight,
                           std::vector<double> & philVal,
                           TPZFMatrix & dphi,TPZFMatrix & elementK,
                           TPZFMatrix & elementF) const
{
	
  int i, j, nshape;
	nshape = philVal.size();
	
	if (fForcing) 
	{
		std::vector <double> force;
		fForcing(point, force);
		for(i=0; i<nshape; i++)
		{
			elementF(i,0) += weight*philVal[i]*force[0];
			for(j=0; j<nshape; j++)
			{
				elementK(i,j) += dphi(0,i)*dphi(0,j)*fK*weight+
				philVal[i]*dphi(0,j)*weight*fC+
				philVal[i]*philVal[j]*fB*weight;
			}
		}
		
	}
	else
	{
		for(i=0; i<nshape; i++)
		{
			elementF(i,0) += weight*philVal[i]*fF;
			for(j=0; j<nshape; j++)
			{
				elementK(i,j) += dphi(0,i)*dphi(0,j)*fK*weight + dphi(1,i)*dphi(1,j)*fK*weight;//+
				//philVal[i]*philVal[j]*fB*weight;
			}
		}
		
	}
}

   /**
    * Calcula a contribuicao para o erro da solucao
    * @param weight [in] peso do ponto de integracao
    * @param sol [in] valor da solucao
    * @param deriv [in] valor da derivada da solucao
    * @param function [in] ponteiro para funcao que calcula o valor exato
    * @param energy [in/out] contribuicao para norma da energia
    * @param l2 [in/out] contribuicao para norma em L2
    *
    */
void TMaterial2d::ContributeErrorSquare(std::vector<double> &x, double weight, double sol, std::vector<double> &deriv,
   	void (function)(std::vector<double>& x, double &val, std::vector<double>&der), double &energy, double &l2)
{
	double solexact;
	int dimension = deriv.size();
	std::vector <double> derivsolexact(dimension);
	function(x,solexact,derivsolexact);
	int id;
	for(id=0; id<dimension; id++)
	{
		std::cout << deriv[id] << "\t"<< derivsolexact[id] << std::endl;
	  energy += weight*(deriv[id]-derivsolexact[id])*(deriv[id]-derivsolexact[id]);
		//energy += weight*(deriv[id]);
	}
	//energy += weight*(sol-solexact)*(sol-solexact);
	l2 += weight*(sol-solexact)*(sol-solexact);
	
} 
