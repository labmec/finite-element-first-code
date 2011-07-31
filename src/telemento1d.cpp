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
#include "telemento1d.h"

#include "tnos.h"
#include "tmalha.h"
#include "tmaterial.h"
#include "pzfmatrix.h"
#include "pzquad.h"
#include "tintegracao.h"

using namespace std; 

TElemento1d::TElemento1d(int matid, int order, std::vector< int >& nodes): TElemento(matid, order, nodes)
{
}


TElemento1d::TElemento1d(): TElemento()
{
}


TElemento1d::~TElemento1d()
{
}


MElementType TElemento1d::getElType()
{
    return ELinear;
}

void TElemento1d::CalcStiff(TMalha &malha, TPZFMatrix& stiff, TPZFMatrix& rhs)
{
	int Nshape = fPorder+1, dim = 1;
  std::vector<double> phi(Nshape),pointstl(1,0.);
  TPZFMatrix dphi(1,Nshape), dphix(1,Nshape);
  TPZFMatrix jac(1,1),jacinv(1,1);
  // TPZInt1d intrule(fPorder+fPorder);
    
  //*************** Chama a Regra de Integracao Numerica ***************//
    
  TIntegracao intrule(fPorder+fPorder); // A ordem do polinomio que pode ser integrado eh (2x num de ptos - 1)
  stiff.Redim(Nshape,Nshape);
  rhs.Redim(Nshape,1);
  TMaterial *mat = malha.getMaterial(fMaterialId);
  
  //int npoints = intrule.NPoints();
    int npoints = intrule.NPoints(); //Retorna o tamanho de fPontos
    
  //  cout << "numpontospz="<<npoints<<"numpointsmeu="<<npoints2<<endl;  !!!!! Teste para checar se o num de pontos da nova regra de integracao estah igual ao anterior !!!!
  
    TPZVec<double> point(1);
    
//  double weight;
    double x,weight;
    int ip;
    for(ip = 0; ip<npoints; ip++)
		{
		 // intrule.Point(ip,point,weight);
				
			intrule.Point(ip, x, weight); // Retorna os pontos e os pesos obtidos pela Regra de Integracao
				
			pointstl[0] = x;
			
			this->fphi.resize(Nshape);
			this->fdphi.Redim(dim, Nshape);
			// compute the value of the shape functions
			Shape(pointstl,phi,dphi);
			for (int ishape = 0; ishape < Nshape; ishape++) 
			{
				fphi[ishape] = phi[ishape];
				for (int idim = 0 ; idim < dim; idim++) 
				{
					fdphi(idim,ishape) = dphi(idim,ishape);				
				}
			}
			
			double detjac;
			// compute the jacobian
			Jacobian(pointstl,jac,jacinv,detjac,malha);
			// compute the derivative of the shapefunctions with respect to x
			int i;
			for(i=0; i<phi.size(); i++)
			{
				dphix(0,i) = jacinv(0,0)*dphi(0,i);
			}
			dphix.Print("Derivada real da funcao de forma");
				
			weight *= fabs(detjac);  // !!! Ajusta o peso do ponto de integracao ao determinante do jacobiano !!!
			cout << "Peso " << weight << endl;
			
			RealCoord(pointstl,malha); // Acha as coord reais do ponto de integracao. Para usar na forcingfunction depois
			
			// accumulate the contribution to the stiffness matrix
			mat->Contribute(pointstl,weight,phi,dphix,stiff,rhs);
		}
}
//**********************************************************************//

    /**
     * Calculo do jacobiano
     * @point : ponto em qual calculamos o jacobiano do mapeamento
     * @jacobian : matriz jacobiana (a dimensao depende da dimensao do elemento
     * @jacinv : inverso da matriz jacobiana
     * @malha : objeto malha necessaria para relacionar os indices dos nos com os nos reais
     * 
     */
void TElemento1d::Jacobian(std::vector<double> &point, TPZFMatrix &jacobian, TPZFMatrix &jacinv, double &detjac, TMalha &malha)
{
  if(!fNodes.size()) return;
  //TNo &no1 = malha.getNode(fNodes[0]);
  //TNo &no2 = malha.getNode(fNodes[fPorder]);
  
	jacobian.Resize(1,1);
  jacinv.Resize(1,1);
	jacobian.Zero();
	jacinv.Zero();
	double dxdn = 0, dydn = 0;
	for (int i = 0; i < fPorder+1 ; i++) 
	{
		TNo &no = malha.getNode(fNodes[i]);
		dxdn += no.Co(0)*fdphi(0,i);
		dydn += no.Co(1)*fdphi(0,i);
	}
	
	jacobian(0,0) = sqrt(dxdn*dxdn + dydn*dydn);
	jacinv = 1/jacobian(0,0);
	detjac = jacobian(0,0);
	
	/* MODO ANTIGO
  double delx = sqrt(
    (no1.Co(0)-no2.Co(0))*(no1.Co(0)-no2.Co(0))+
    (no1.Co(1)-no2.Co(1))*(no1.Co(1)-no2.Co(1))
    );
  jacobian.Resize(1,1);
  jacinv.Resize(1,1);
  jacobian(0,0) = delx/2.;
  jacinv(0,0) = 2./delx;
  detjac = jacobian(0,0);*/
	 
}


    /**
     * Calcula os valores das funcoes de forma e suas derivadas
     * @point ponto onde calcular as funcoes de forma
     * @phi valores das funcoes de forma
     * @dphi valores das derivadas das funcoes de forma
     */
void TElemento1d::Shape(std::vector<double> &point, std::vector<double> &phi, TPZFMatrix &dphi)
{
  TElemento::Shape1d(fPorder, point, phi, dphi);
}

    /**
     * Calcula o erro do elemento
     * @param exact funcao que calcula o valor exato
     * @param energy [out] erro na norma de energia
     * @param l2 [out] erro na norma l2
     */
void TElemento1d::Error(TPZFMatrix &solution, TMalha &malha, void (exact)(std::vector<double> &,double &, std::vector<double> &), double &energy, double &l2)
{
  std::vector<double> phi(this->fPorder+1),pointstl(1,0.);
  TPZFMatrix dphi(1,fPorder+1), dphix(1,fPorder+1);
  TPZFMatrix jac(1,1),jacinv(1,1);
  TIntegracao intrule(fPorder+fPorder);
  
  TMaterial *mat = malha.getMaterial(fMaterialId);
  TNo &no1 = malha.getNode(fNodes[0]);
  TNo &no2 = malha.getNode(fNodes[fPorder]);
  
  int npoints = intrule.NPoints();
  double point;
  std::vector<double> x(1,0.);
  double weight;
  energy = 0.;
  l2 = 0.;
  int numinterval = 50;
  int interval;
  for(interval = 0; interval<numinterval; interval++)
  {
			int ip;
			for(ip = 0; ip<npoints; ip++)
			{
				intrule.Point(ip,point,weight);
				pointstl[0] = -1.+interval*2./numinterval+(point+1.)/(numinterval);
				// compute the jacobian
				double detjac;
				Jacobian(pointstl,jac,jacinv,detjac,malha);
				
				x[0] = no1.Co(0) + (1.+pointstl[0])*(no2.Co(0)-no1.Co(0))/2.;
				// compute the value of the shape functions
				Shape(pointstl,phi,dphi);
				// compute the derivative of the shapefunctions with respect to x
				int i;
				for(i=0; i<phi.size(); i++)
				{
					dphix(0,i) = jacinv(0,0)*dphi(0,i);
				}
				dphix.Print("Derivada real da funcao de forma");
				// adjust the weight of the integration point with the determinant of the jacobian
				weight *= fabs(detjac)/numinterval;
				cout << "Peso " << weight << endl;
				double sol = 0.;
				std::vector<double> derivsol(1,0.);
				for(i=0; i<phi.size(); i++)
				{
					sol += phi[i]*solution(fNodes[i],0);
					derivsol[0] += dphix(0,i)*solution(fNodes[i],0);
				}
				// accumulate the contribution to the stiffness matrix
				mat->ContributeErrorSquare(x,weight,sol,derivsol,exact,energy,l2);
			}
   }
   energy = sqrt(energy);
   l2 = sqrt(l2);

}

void TElemento1d::RealCoord(std::vector <double> &point, TMalha &malha)
{
	double x, y;
	point.resize(2,0);
	for (int i = 0 ; i < fphi.size(); i++) 
	{
		TNo &no = malha.getNode(fNodes[i]);
		x += no.Co(0)*fphi[i];
		y += no.Co(1)*fphi[i];
	}
	point[0] = x;
	point[1] = y;
}
