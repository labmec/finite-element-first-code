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
#include "telemento.h"
#include "tnos.h"
#include "tmalha.h"
#include "pzfmatrix.h"

TElemento::TElemento() : fNodes(), fMaterialId(-1), fPorder(-1)
{
  TPZFMatrix a(4,4,0.);
  
}

/**
  *  Elemento criado com
  * @matid : indice do material
  * @order : ordem de interpolacao
  * @nodes : indices dos nos i dont know
  * nota que existe correspondencia entre o numero de nos e a ordem do elemento
  */
TElemento::TElemento(int matid, int order, std::vector<int> &nodes) :
  fNodes(nodes), fMaterialId(matid), fPorder(order)
{
}


TElemento::~TElemento()
{
}

    /**
     * Cálcula o valor da função de forma em um dado ponto
     * @param pt [in] ponto onde se quer calcular o valor das funções de forma
     * @param phiValue [out] valor de cada função de forma do elemento no ponto
     */
void TElemento::Shape1d (int order, std::vector<double> & pt,std::vector<double> & phi, TPZFMatrix& dphi)
{
  std::cout << __PRETTY_FUNCTION__ << "\nOrdem : " << order << " Ponto de entrada ";
  int i,j,k;
  for(i=0; i< pt.size(); i++) std::cout << pt[i] << ' ';
  std::cout << std::endl;
  
  phi.resize(order+1);
  dphi.Redim(1,order+1);
  for(i=0; i<= order; i++)
  {
    phi[i] = 1.;
    dphi(0,i) = 0.;
  }
  for(i=0; i<= order; i++)
  {
    double pointi = -1.+(i*2./order);
    for(j=0; j<=order; j++)
    {
      double pointj = -1.+(j*2./order);
      if(j!=i) 
      {
        phi[i] *= (pt[0]-pointj)/(pointi-pointj);
        double prod = 1/(pointi-pointj);
        for(k=0; k <= order; k++)
        {
          double pointk = -1.+(k*2./order);
          if(k!=j && k!=i) prod *= (pt[0]-pointk)/(pointi-pointk);
        }
        dphi(0,i) += prod;
      }
    }
  }
}

int TElemento::main()
{
  std::vector<double> pt(1), phi(3);
  TPZFMatrix dphi(1,3);
  pt[0] = 0.5;
  Shape1d(3,pt,phi,dphi);
}



/*!
    \fn TElemento::Print(std::ostream &out = cout)
 */
void TElemento::Print(std::ostream &out)
{
  out << "ElType " << TypeName(getElType()) << " matid " << fMaterialId << " order " << fPorder << " nodes ";
  int i;
  for(i=0; i<fNodes.size(); i++) out << fNodes[i] << ' ';
  out << std::endl;
}


/*!
    \fn TElemento::TypeName(MElementType type)
 */
std::string TElemento::TypeName(MElementType type)
{
  std::string result;
    switch (type)
    {
      case EPoint:
      result = "Point";
      break;
      case ELinear:
      result = "Linear";
      break;
      case EQuadrilateral:
      result = "Quadrilateral";
      break;
      case  ETriangle:
      result = "Triangle";
      break;
      default:
      result = "Unknown";
    }
    return result;
}


