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
#include "telemento2d_quad.h"

#include "tnos.h"
#include "tmalha.h"
#include "tmaterial.h"
#include "pzfmatrix.h"
#include "pzquad.h"
#include "tintegracao2d_quad.h"

using namespace std; 

telemento2d_quad::telemento2d_quad(int matid, int order, std::vector< int >& nodes): TElemento(matid, order, nodes)
{
}


telemento2d_quad::telemento2d_quad(): TElemento()
{
}


telemento2d_quad::~telemento2d_quad()
{
}


MElementType telemento2d_quad::getElType()
{
    return EQuadrilateral;
}



void telemento2d_quad::CalcStiff(TMalha &malha, TPZFMatrix& stiff, TPZFMatrix& rhs)
{
	int Nshape=(this->fPorder+1)*(this->fPorder+1);
	int dim=2;
  std::vector<double> phi(Nshape),pointstl(dim,0.);
  TPZFMatrix dphi(dim,Nshape), gradphi(dim,Nshape);
  TPZFMatrix jac(dim,dim),jacinv(dim,dim);
  // TPZInt1d intrule(fPorder+fPorder);
    
  //*************** Chama a Regra de Integracao Numerica ***************//
    
  TIntegracao2d_quad intrule(fPorder+fPorder); // A ordem do polinomio que pode ser integrado eh (2x num de ptos - 1)
  stiff.Redim(Nshape,Nshape);
  rhs.Redim(Nshape,1);
  TMaterial *mat = malha.getMaterial(fMaterialId);
  
  //int npoints = intrule.NPoints();
    int npoints = intrule.NPoints(); //Retorna o tamanho de fPontos
    
  //  cout << "numpontospz="<<npoints<<"numpointsmeu="<<npoints2<<endl;  !!!!! Teste para checar se o num de pontos da nova regra de integracao estah igual ao anterior !!!!
  
	std::pair<double,double> x;
    
//  double weight;
    double weight;
    int ip;
	for(ip = 0; ip<npoints; ip++)
  {
		//intrule.Point(ip,point,weight);
      
		// AGORA RECEBE UM PAIR NO X
    intrule.Point(ip, x, weight); // Retorna os pontos e os pesos obtidos pela Regra de Integracao
    
		pointstl[0] = x.first;
		pointstl[1] = x.second;

		// compute the value of the shape functions
    Shape(pointstl,phi,dphi);
		
		// Acredito que de pra fazer melhor do que definir dois novos atributos
		this->fphi.resize(Nshape);
		this->fdphi.Redim(dim, Nshape);
		
		// prenchendo fdphi e fphi
		for (int ishape = 0; ishape < Nshape; ishape++) 
		{
			fphi[ishape] = phi[ishape];
			for (int idim = 0 ; idim < dim; idim++) 
			{
				fdphi(idim,ishape) = dphi(idim,ishape);				
			}
		}
		
    // compute the jacobian
    double detjac;
    Jacobian(pointstl,jac,jacinv,detjac,malha);

		
    // compute the derivative of the shapefunctions with respect to x
    int i,j,idim,jdim;
    for(i=0; i<phi.size(); i++)
		{
			for (idim=0; idim<dim; idim++) 
			{
				gradphi(idim,i)=0;
				for (jdim=0; jdim<dim; jdim++) 
				{
					gradphi(idim,i) += jacinv(jdim,idim)*dphi(jdim,i); // Aqui estava jacinv(idim,jdim), troquei, pois acho que o certo eh o contrario
				}
			}
		}
    gradphi.Print("Derivada real da funcao de forma");
      
    weight *= fabs(detjac);  // !!! Ajusta o peso do ponto de integracao ao determinante do jacobiano !!!

    cout << "Peso " << weight << endl;
		
		RealCoord(pointstl,malha); // Acha as coord reais do ponto de integracao. Para usar na forcingfunction depois
		
    // accumulate the contribution to the stiffness matrix
    mat->Contribute(pointstl,weight,phi,gradphi,stiff,rhs);
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


void telemento2d_quad::Jacobian(std::vector<double> &point, TPZFMatrix &jacobian, TPZFMatrix &jacinv, double &detjac, TMalha &malha)
{
  if(!fNodes.size()) return;
  jacobian.Resize(2,2);	
	jacinv.Resize(2,2);
	jacobian.Zero();
	jacinv.Zero();
	int in;
	for (in = 0 ; in < fNodes.size() ; in++) 
	{
		TNo &no = malha.getNode(fNodes[in]);
		//std::cout << "no x = " << no.Co(0) << "\t" << "no y = " << no.Co(1) << std::endl;
		jacobian(0,0) += no.Co(0)*fdphi(0,in);
		jacobian(0,1) += no.Co(0)*fdphi(1,in);
		jacobian(1,0) += no.Co(1)*fdphi(0,in);
		jacobian(1,1) += no.Co(1)*fdphi(1,in);
	}
	//jacobian.Print("Jac:");
	detjac = jacobian(0,0)*jacobian(1,1)-jacobian(0,1)*jacobian(1,0);
	double invdetjac = 1/detjac;
	jacinv(0,0) = invdetjac*jacobian(1,1);
	jacinv(0,1) = -invdetjac*jacobian(0,1);
	jacinv(1,0) = -invdetjac*jacobian(1,0);
	jacinv(1,1) = invdetjac*jacobian(0,0);
	//jacinv.Print("JacInv:");
	
  /*TNo &no1 = malha.getNode(fNodes[0]);
  TNo &no2 = malha.getNode(fNodes[fPorder]);
  double delx = sqrt(
    (no1.Co(0)-no2.Co(0))*(no1.Co(0)-no2.Co(0))+
    (no1.Co(1)-no2.Co(1))*(no1.Co(1)-no2.Co(1))
    );

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


void telemento2d_quad::Shape(std::vector<double> &point, std::vector<double> &phi, TPZFMatrix &dphi)
{
	std::vector<double> phi1d_a(fPorder+1), phi1d_b(fPorder+1), pt(1,0);
	TPZFMatrix dphi1d_a(0,fPorder+1), dphi1d_b(0,fPorder+1);
	pt[0] = point[0];
  TElemento::Shape1d(fPorder, pt, phi1d_a, dphi1d_a);
	pt[0] = point[1];
	TElemento::Shape1d(fPorder, pt, phi1d_b, dphi1d_b);
	
	/*int is, js; MUITO LEGAL, PENA QUE NAO FICA ANTI HORARIO
	for (is = 0 ; is < (fPorder+1) ; is++)
	{
		for (js = 0 ; js < (fPorder+1) ; js++ ) 
		{
			phi[is*(fPorder+1)+js] = phi1d_a[js]*phi1d_b[is];
			dphi(0,is*(fPorder+1)+js) = dphi1d_a[js]*phi1d_b[is];
			dphi(1,is*(fPorder+1)+js) = phi1d_a[js]*dphi1d_b[is];
		}
	}*/
	
	
	
	int k = 3, ns = fPorder;
	
	// Fiz um monte de blablabla mas no final soh funciona pra p=2
	phi[0] = phi1d_a[0]*phi1d_b[0];	
	phi[1] = phi1d_a[ns]*phi1d_b[0];		
	phi[2] = phi1d_a[ns]*phi1d_b[ns];		
	phi[3] = phi1d_a[0]*phi1d_b[ns];		
	dphi(0,0) = dphi1d_a[0]*phi1d_b[0];
	dphi(0,1) = dphi1d_a[ns]*phi1d_b[0];
	dphi(0,2) = dphi1d_a[ns]*phi1d_b[ns];
	dphi(0,3) = dphi1d_a[0]*phi1d_b[ns];
	dphi(1,0) = phi1d_a[0]*dphi1d_b[0];
	dphi(1,1) = phi1d_a[ns]*dphi1d_b[0];
	dphi(1,2) = phi1d_a[ns]*dphi1d_b[ns];
	dphi(1,3) = phi1d_a[0]*dphi1d_b[ns];
	
	int im;
	for (im = 1; im < ns; im++) 
	{
		k++;
		phi[k] = phi1d_a[im]*phi1d_b[0];
		dphi(0,k) = dphi1d_a[im]*phi1d_b[0];
		dphi(1,k) = phi1d_a[im]*dphi1d_b[0];

	}
	for (im = 1; im < ns; im++) 
	{
		k++;
		phi[k] = phi1d_a[ns]*phi1d_b[im];
		dphi(0,k) = dphi1d_a[ns]*phi1d_b[im];
		dphi(1,k) = phi1d_a[ns]*dphi1d_b[im];
	}
	for (im = ns-1; im > 0; im--) 
	{
		k++;
		phi[k] = phi1d_a[im]*phi1d_b[ns];
		dphi(0,k) = dphi1d_a[im]*phi1d_b[ns];
		dphi(1,k) = phi1d_a[im]*dphi1d_b[ns];
	}
	for (im = ns-1; im > 0; im--) 
	{
		k++;
		phi[k] = phi1d_a[0]*phi1d_b[im];
		dphi(0,k) = dphi1d_a[0]*phi1d_b[im];
		dphi(1,k) = phi1d_a[0]*dphi1d_b[im];
	}
	if (k == 7)
	{
		k++;
		phi[k] = phi1d_a[im]*phi1d_b[im];
		dphi(0,k) = dphi1d_a[0]*phi1d_b[im];
		dphi(1,k) = phi1d_a[0]*dphi1d_b[im];
	}

}


    /**
     * Calcula o erro do elemento
     * @param exact funcao que calcula o valor exato
     * @param energy [out] erro na norma de energia
     * @param l2 [out] erro na norma l2
     */
void telemento2d_quad::Error(TPZFMatrix &solution, TMalha &malha, void (exact)(std::vector<double> &,double &, std::vector<double> &), double &energy, double &l2)
{
	int Nshape=(this->fPorder+1)*(this->fPorder+1);
	int dim=2;
  std::vector<double> phi(Nshape),pointstl(dim,0.);
  TPZFMatrix dphi(dim,Nshape), gradphi(dim,Nshape);
  TPZFMatrix jac(dim,dim),jacinv(dim,dim);
  // TPZInt1d intrule(fPorder+fPorder);
	
  //*************** Chama a Regra de Integracao Numerica ***************//
	
  TIntegracao2d_quad intrule(6); // A ordem do polinomio que pode ser integrado eh (2x num de ptos - 1)
  TMaterial *mat = malha.getMaterial(fMaterialId);
  
  //int npoints = intrule.NPoints();
	int npoints = intrule.NPoints(); //Retorna o tamanho de fPontos
	
  //  cout << "numpontospz="<<npoints<<"numpointsmeu="<<npoints2<<endl;  !!!!! Teste para checar se o num de pontos da nova regra de integracao estah igual ao anterior !!!!
  
	std::pair<double,double> x;
	
	//  double weight;
	double weight;
	int ip;
	for(ip = 0; ip<npoints; ip++)
  {
		//intrule.Point(ip,point,weight);
		
		// AGORA RECEBE UM PAIR NO X
    intrule.Point(ip, x, weight); // Retorna os pontos e os pesos obtidos pela Regra de Integracao
    
		pointstl[0] = x.first;
		pointstl[1] = x.second;
		
		// compute the value of the shape functions
    Shape(pointstl,phi,dphi);
		
		// Acredito que de pra fazer melhor do que definir dois novos atributos
		this->fphi.resize(Nshape);
		this->fdphi.Redim(dim, Nshape);
		
		// prenchendo fdphi e fphi
		for (int ishape = 0; ishape < Nshape; ishape++) 
		{
			fphi[ishape] = phi[ishape];
			for (int idim = 0 ; idim < dim; idim++) 
			{
				fdphi(idim,ishape) = dphi(idim,ishape);				
			}
		}
		
    // compute the jacobian
    double detjac;
    Jacobian(pointstl,jac,jacinv,detjac,malha);
		
		
    // compute the derivative of the shapefunctions with respect to x
    int i,j,idim,jdim;
    for(i=0; i<phi.size(); i++)
		{
			for (idim=0; idim<dim; idim++) 
			{
				gradphi(idim,i)=0;
				for (jdim=0; jdim<dim; jdim++) 
				{
					gradphi(idim,i) += jacinv(jdim,idim)*dphi(jdim,i); // Aqui estava jacinv(idim,jdim), troquei, pois acho que o certo eh o contrario
				}
			}
		}

		
    weight *= fabs(detjac);  // !!! Ajusta o peso do ponto de integracao ao determinante do jacobiano !!!
		
		
		RealCoord(pointstl,malha); // Acha as coord reais do ponto de integracao. Para usar na forcingfunction depois
		
    // accumulate the contribution to the stiffness matrix
		//TMaterial2d::ContributeErrorSquare(std::vector<double> &x, double weight, double sol, std::vector<double> &deriv,
   	//void (function)(std::vector<double>& x, double &val, std::vector<double>&der), double &energy, double &l2
		double sol;
		std::vector <double> deriv;
		this->ComputeSol(solution, sol, deriv);
		mat->ContributeErrorSquare(pointstl, weight, sol, deriv, exact, energy, l2);
    //mat->Contribute(pointstl,weight,phi,gradphi,stiff,rhs);
  }
	
}

// Only made this method to test the jacobian in a Unit Test (better ways?)
void telemento2d_quad::JacobTest(TMalha &malha, TPZFMatrix& stiff, TPZFMatrix& rhs, std::vector <double> &ponto, double &detjacob)
{
	int Nshape=(this->fPorder+1)*(this->fPorder+1);
	int dim=2;
  std::vector<double> phi(Nshape),pointstl(dim,0.);
  TPZFMatrix dphi(dim,Nshape), gradphi(dim,Nshape);
  TPZFMatrix jac(dim,dim),jacinv(dim,dim);

	
	pointstl[0] = ponto[0];
	pointstl[1] = ponto[1];
	
	Shape(pointstl,phi,dphi);
	
	// Acredito que de pra fazer melhor do que definir dois novos atributos
	this->fphi.resize(Nshape);
	this->fdphi.Redim(dim, Nshape);
	for (int ishape = 0; ishape < Nshape; ishape++) 
	{
		fphi[ishape] = phi[ishape];
		for (int idim = 0 ; idim < dim; idim++) 
		{
			fdphi(idim,ishape) = dphi(idim,ishape);				
		}
	}
	
	// compute the jacobian
	double detjac;
	Jacobian(pointstl,jac,jacinv,detjac,malha);
	detjacob = detjac;
	
}

void telemento2d_quad::ComputeSol(TPZFMatrix& sol, double &val, std::vector <double> &deriv)
{
	int n = fNodes.size();
	val = 0;
	deriv.resize(2,0);
	for (int i = 0; i < n; i++) 
	{
		val += sol(fNodes[i])*fphi[i];
		deriv[0] = sol(fNodes[i])*fdphi(0,i);
		deriv[1] = sol(fNodes[i])*fdphi(1,i);		
	}
}

void telemento2d_quad::RealCoord(std::vector <double> &point, TMalha &malha)
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
