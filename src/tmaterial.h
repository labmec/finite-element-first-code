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
#ifndef TMATERIAL_H
#define TMATERIAL_H

#include <vector>
#include <iostream>
#include "pzfmatrix.h"

class TPZFMatrix;

/**
Implementa a interface t�pica de um material


@author Philippe R. B. Devloo
*/
class TMaterial{
public:
    /**
     * Construtor vazio
     */
    TMaterial(int id);
    
    /**
     * Construtor de c�pia
     */
    TMaterial (const TMaterial & copy);

    /**
     * Destrutor
     */
    virtual ~TMaterial();
    
public: //fun��es de c�lculo

  /**
   * Calcula o valor da contribui��o da equa��o variacional no ponto dado
   * na matriz de rigidez do elemento e no vetor de carga
   * @param pt [in]: ponto de integra��o de Gauss
   * @param weight [in]: peso de integra��o
   * @param phiVal [in] : valor da fun��o teste no ponto dado
   * @param dphi [in] : valor das derivadas da fun��o de forma no ponto de integra��o
   * @param elementK [inout]: matriz de rigidez do elemento
   * @param elementF [inout]: vetor de carga do elemento
   */
  virtual void Contribute (std::vector<double> &point, double  weight,
                           std::vector<double> & philVal,
                           TPZFMatrix & dphi,TPZFMatrix & elementK,
                           TPZFMatrix & elementF) const = 0; 

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
   	void (function)(std::vector<double>& x, double &val, std::vector<double>&der), double &energy, double &l2) = 0;
	
	
public://Diversos
  /**
   * Imprime os dados da classe
   */
  virtual void Print(std::ostream &out = std::cout) const;
	
	virtual void SetForcingFunction(void (*func)(std::vector<double> &point, std::vector<double> &force))
	{
		fForcing = func;
	}
  
  /**
   * retorna o identificador do material
   */
   int Id() 
   {
    return fId;
   }
  
protected:
  /**
   * Identificador do material
   */
   int fId;
   
	void (*fForcing)(std::vector<double> &point, std::vector<double> &force);
};

#endif
