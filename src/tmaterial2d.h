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
#ifndef TMATERIAL2D_H
#define TMATERIAL2D_H

#include <tmaterial.h>

/**
This class implements the variational formulation of a 2D partial differential equation

@author Philippe R. B. Devloo
*/
class TMaterial2d : public TMaterial
{
public:
    
    TMaterial2d(int id, double K, double C, double B, double F);

    ~TMaterial2d();

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
  virtual void Contribute (std::vector<double> &point, double  weight,
                           std::vector<double> & philVal,
                           TPZFMatrix & dphi,TPZFMatrix & elementK,
                           TPZFMatrix & elementF) const; 
   /**
    * Calcula a contribuicao para o erro da solucao
    * @param x [in] localizacao do ponto
    * @param weight [in] peso do ponto de integracao
    * @param sol [in] valor da solucao
    * @param deriv [in] valor da derivada da solucao
    * @param function [in] ponteiro para funcao que calcula o valor exato
    * @param energy [in/out] contribuicao para norma da energia
    * @param l2 [in/out] contribuicao para norma em L2
    *
    */
   virtual void ContributeErrorSquare(std::vector<double> &x, double weight, double sol, std::vector<double> &deriv,
   	void (function)(std::vector<double>& x, double &val, std::vector<double>&der), double &energy, double &l2)  ;                      

    virtual void Print(std::ostream& out) const;

protected:

  /**
   * Definition of the differential equation coeficients according to the book of Becker, Carey and Oden
   */
  double fK, fC, fB, fF;
  
};

#endif
