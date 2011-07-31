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
#include "telemento0d.h"
#include "tmaterial.h"
#include "tmalha.h"
#include "pzfmatrix.h"

TElemento0d::TElemento0d()
{
}


TElemento0d::TElemento0d(int matid, int order, std::vector< int >& nodes): TElemento(matid, order, nodes)
{
}


TElemento0d::~TElemento0d()
{
}


std::string TElemento0d::TypeName(MElementType type)
{
    return TElemento::TypeName(type);
}

void TElemento0d::CalcStiff(TMalha &malha, TPZFMatrix& stiff, TPZFMatrix& rhs)
{
  stiff.Redim(1,1);
  rhs.Redim(1,1);
  std::vector<double> phi(1,1.);
  TPZFMatrix dphi(0,1);
  std::vector<double> point(0);
  Shape(point,phi,dphi);
  double weight = 1.;
  TMaterial *mat = malha.getMaterial(this->fMaterialId);
  mat->Contribute(point,weight,phi,dphi,stiff,rhs);
}

void TElemento0d::Jacobian(std::vector< double >& point, TPZFMatrix& jacobian, TPZFMatrix& jacinv, double &detjac, TMalha& malha)
{
  jacobian.Resize(0,0);
  jacinv.Resize(0,0);
  detjac = 1;
}



/*!
    \fn TElemento0d::getElType()
 */
MElementType TElemento0d::getElType()
{
    return EPoint;
}

    /**
     * Calcula os valores das funcoes de forma e suas derivadas
     * @point ponto onde calcular as funcoes de forma
     * @phi valores das funcoes de forma
     * @dphi valores das derivadas das funcoes de forma
     */
void TElemento0d::Shape(std::vector<double> &point, std::vector<double> &phi, TPZFMatrix &dphi) 
{
  phi.resize(1);
  phi[0] = 1.;
  dphi.Resize(0,1);
}

void TElemento0d::RealCoord(std::vector <double> &point, TMalha &malha)
{
	
}
