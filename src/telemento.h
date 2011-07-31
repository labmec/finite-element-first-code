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
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *t Wt WIIt WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GN    
    /**
    
    TNos * getNode(i);
U General Public License for more details.          t WI                *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef TELEMENTO_H
#define TELEMENTO_H

class TMalha;

class TPZMatrix;
class TPZFMatrix;
#include <vector>
#include <iostream>
#include "pzfmatrix.h"
/**
Classe virtual que implementa a interface b�sica de um elemento

@author Philippe R. B. Devloo
*/
/**
 * Defini��o dos tipos de elementos
 */
 
enum MElementType { EPoint, ELinear, ETriangle, EQuadrilateral };

class TElemento{
public:
    /**
     * Construtor vazio - inicializar com valores que indiquem falta de 
     * inicializa��o
     */
    TElemento();
    
    /**
     *  Elemento criado com
     * @matid : indice do material
     * @order : ordem de interpolacao
     * @nodes : indices dos nos
     * nota que existe correspondencia entre o numero de nos e a ordem do elemento
     */
     TElemento(int matid, int order, std::vector<int> &nodes);
    
    /**
     * Destrutor padr�o
     */
    ~TElemento();
    
    /**
     * Metodo para acesso aos nos do elemento
     */
     const std::vector<int> &getNodeVec() const
     {
       return fNodes;
     }
    
public: //Acesso a leitura dos dados da classe
    /**
     * Retorna o tipo de elemento
     */
    virtual MElementType getElType()=0;

public: //M�todos de c�lculo
    
    static int main();
    
public:
	


    /**
     * C�lcula o valor da fun��o de forma em um dado ponto
     * @param pt [in] ponto onde se quer calcular o valor das fun��es de forma
     * @param phiValue [out] valor de cada fun��o de forma do elemento no ponto
     */
    static void Shape1d (int order, std::vector<double> & pt,std::vector<double> & phi, TPZFMatrix& dphi);
    
    /**
     * Calcula a matriz de rigidez local
     */
    virtual void CalcStiff (TMalha &malha, TPZFMatrix & stiff, TPZFMatrix &rhs) = 0;
    
    /**
     * Calculo do jacobiano
     * @point : ponto em qual calculamos o jacobiano do mapeamento
     * @jacobian : matriz jacobiana (a dimensao depende da dimensao do elemento
     * @jacinv : inverso da matriz jacobiana
     * @malha : objeto malha necessaria para relacionar os indices dos nos com os nos reais
     * 
     */
     virtual void Jacobian(std::vector<double> &point, TPZFMatrix &jacobian, TPZFMatrix &jacinv, double &detjac, TMalha &malha) = 0;
     
     /**
      * return the string corresponding to the type
      */
    static std::string TypeName(MElementType type);
    
    /**
     * Print the element
     */
    void Print(std::ostream &out = std::cout);
    
    /**
     * Calcula os valores das funcoes de forma e suas derivadas
     * @point ponto onde calcular as funcoes de forma
     * @phi valores das funcoes de forma
     * @dphi valores das derivadas das funcoes de forma
     */
    virtual void Shape(std::vector<double> &point, std::vector<double> &phi, TPZFMatrix &dphi) = 0;
    
    /**
     * Calcula o erro do elemento
     * @param exact funcao que calcula o valor exato
     * @param energy [out] erro na norma de energia
     * @param l2 [out] erro na norma l2
     */
	virtual void Error(TPZFMatrix &solution, TMalha &malha, void (f)(std::vector<double> &,double &, std::vector<double> &), double &energy, double &l2) = 0;
    
	
	/** Calcula as coord reais do ponto de integracao
	 
	 **/
	virtual void RealCoord(std::vector <double> &point, TMalha &malha) = 0;
	
	
protected:
    
    /**
     * Vetor de ponteiros para os n�s. Os n�s s�o armazenados
     * na malha
     */
    std::vector<int> fNodes;
    
    /**
     * Refer�ncia para o �ndice do material na malha
     */
    int fMaterialId;
    
    /**
     * Ordem de interpola��o p
     */
    int fPorder;
		
	  /**
	  * Funcoes base
	  */
		std::vector <double> fphi;
		
	  /**
	  * Derivada das funcoes base
	  */
		TPZFMatrix fdphi;
    
};

#endif
