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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <iostream>
#include <cstdlib>
#include <math.h>
#include "tnos.h"
#include "telemento.h"
#include "telemento0d.h"
#include "telemento1d.h"
#include "tmalha.h"
#include "tanalysis.h"
#include "tmaterial1d.h"
#include "tmaterialbc.h"
#include "pzquad.h"
#include "pzvec.h"
#include "telemento2d_quad.h"
#include "tmaterial2d.h"
#include "telemento2d_tri.h"
//#include "log4cxx.h"


//using namespace std;

double Integrate(int ord, double xmin, double xmax,double (*g)(double x));

double func(double x, double y);

void ReadMesh(TMalha &malha, std::string &filename);

void ReadMeshGid(TMalha &malha, std::string &Filename);

void ReadMeshGidTri(TMalha &malha, std::string &Filename);

void ForcingFunction(std::vector<double> &point, std::vector<double> &force);

void Exact(std::vector<double> &x, double &val, std::vector<double> &deriv);

void ExactQuad(std::vector <double> &x, double &val, std::vector <double> &deriv);


int main(int argc, char *argv[])
{
/*  std::cout << "Hello, world!" << std::endl;
  
  std::cout << "Testando classe nó" << std::endl;
  TNo::main();
  std::cout << "Fim do teste da classe nó" << std::endl;

  std::cout << "Testando classe elemento" << std::endl;
  TElemento::main();
  std::cout << "Fim do teste da classe elemento" << std::endl;
  
  std::cout << "Testando uma malha uni dimensional" << std::endl;
  TMalha::main();
  std::cout << "Fim do teste da class malha" << std::endl;
  
  double result;
  int ord;
  for (ord = 1; ord < 10; ord++)
  {
    result = Integrate(ord,0.,3.,func);
    std::cout << "O resultado da integracao " << result << std::endl;
  }*/

	
			
  if(argc < 2) 
  {
      std::cout << "Numero de argumentos menor que dois";
        return -1;
  }
  TMalha NumeroUm;
  std::string filename(argv[1]);
  //ReadMeshGidTri(NumeroUm,filename);
  ReadMeshGid(NumeroUm,filename);
	//ReadMesh(NumeroUm,filename);
  TAnalysis analysis(&NumeroUm);
  analysis.Run();
	double energy, l2;
	
	analysis.Error(ExactQuad, energy, l2);
	std::cout << energy << "\t" << l2 ;

    
  return EXIT_SUCCESS;
}

double Integrate(int ord, double xmin, double xmax,double (*g)(double x))
{
    double result;
  /* -- 
    
  // escala devido ao tamanho do intervalo
  double mult = (xmax-xmin)/2.;
  // objecto tipo regra de integracao
  TPZInt1d rule(ord);
  // numero de pontos de integracao
  int np = rule.NPoints();
  std::cout << "Numero de pontos " << np << std::endl;
  // estrutura de dados para receber a posicao do ponto e o peso
  TPZVec<double> point(1);
  double w, result = 0.;
  int p;
  // loop sobre os pontos de integracao
  for(p=0; p<np; p++)
  {
    // pegue a posicao do ponto e o peso
    rule.Point(p,point,w);
    // calcula a coordenada x correspondente
    double x = xmin + (point[0]+1.)*(xmax-xmin)/2.;
    // calcula o valor da funcao, integra
    result += g(x)*w*mult;
    
  }*/
  return result;
  
}

double func(double x, double y)
{
  return (1+x)/(2+y);
}



void ReadMesh(TMalha &malha, std::string &filename)
{
   int nmat,nbc,npoint; 
   std::ifstream input(filename.c_str());
   input >> nmat >> nbc >> npoint;
   int im;
   for(im=0; im<nmat; im++) // lendo materiais
   {
		 int id;
		 double k,c,b,f;
		 input >> id >> k >> c >> b >> f;
		 if(f == -1000)
		 {
			 TMaterial1d *mat1d = new TMaterial1d(id,k,c,b,f);
			 mat1d->SetForcingFunction(ForcingFunction);
			 malha.insertMaterial(mat1d);
		 }
		 else 
		 {
			 TMaterial1d *mat1d = new TMaterial1d(id,k,c,b,f);
			 malha.insertMaterial(mat1d);
		 }
   }
   for(im=0; im<nbc; im++) // lendo cond de contorno
   {
   	int id, type;
		 double contrstiff, contrrhs;
		 input >> id >> type >> contrstiff >> contrrhs;
		 TMaterialBC *matbc = new TMaterialBC(id,type,contrstiff,contrrhs);
		 malha.insertMaterial(matbc);
   }
   for(im=0; im < npoint; im++) // lendo os pontos com cond de contorno no meio da malha
   {
   	int mat;
	int node, order = 0;
   	input >> mat >> node;
	std::vector<int> nodes(1,node);
	TElemento0d *elem = new TElemento0d(mat,order,nodes);
	malha.insertElement(elem);
   }
   double x0;
   int matleft, matright;
   input >> x0 >> matleft >> matright; // lendo primeiro noh, matleft matright
   std::vector<double> co(1,x0);
   TNo no(co);
   malha.insertNode(no);
   int count = 0;
   while(input)
   {
   	double x1;
	int nelem, mat, morder;
	input >> x1 >> nelem >> mat >> morder; // lendo proximo nohs
	if(!input) break;
	int numnewnodes = nelem*morder;
	int firstnode = malha.getNodeVec().size();
	malha.getNodeVec().resize(firstnode+numnewnodes);
	int in;
	for(in=1; in<=numnewnodes; in++)
	{
		co[0] = x0+in*(x1-x0)/numnewnodes;
		TNo no(co);
		malha.getNode(firstnode+in-1) = no;
	}
	int el;
	count = firstnode-1;
	for(el=0; el<nelem; el++)
	{
		std::vector<int> nodes(morder+1);
		int in;
		for(in=0; in<=morder; in++)
		{
			nodes[in] = count++;
		}
		TElemento1d *elem = new TElemento1d(mat,morder,nodes);
		malha.insertElement(elem);
		count--;
	}
   }
   int node, order = 0;
   std::vector<int> nodes(1,0);
   TElemento0d *elem = new TElemento0d(matleft,order,nodes); // condicao da esquerda
   malha.insertElement(elem);
   nodes[0] = count;
   elem = new TElemento0d(matright,order,nodes); // condicao da direita
   malha.insertElement(elem);
}

void ReadMeshGid(TMalha &malha, std::string &FileName) 
{
	int numnodes=-1;
	int numelements=-1;
	int numelements1d=-1;
	int id = 1;
	double k = 1, c = 0, b = 0, f = 1;
	TMaterial2d *mat2d = new TMaterial2d(id,k,c,b,f);
	malha.insertMaterial(mat2d);
	
	//string FileName;
	//FileName = "8andares02.txt";
	
	{
		bool countnodes = false;
		bool countelements = false;
		
		std::ifstream read (FileName.c_str());
		
		while(read)
		{
			char buf[1024];
			read.getline(buf, 1024);
			std::string str(buf);
			if(str == "Coordinates") countnodes = true;
			if(str == "end coordinates") countnodes = false;
			if(countnodes) numnodes++;
			
			if(str == "Elements") countelements = true;
			if(str == "end elements")
			{
				countelements = false;
				while (read) 
				{
					char buf[1024];
					read.getline(buf, 1024);
					std::string str(buf);
					if(str == "Coordinates")
					{
						numnodes--;
						countnodes = true;					
					}
					if(str == "end coordinates") countnodes = false;
					if(countnodes) numnodes++;
					
					if(str == "Elements") countelements = true;
					if(str == "end elements") countelements = false;
					if(countelements) numelements1d++;	
				}
				
			}
			if(countelements) numelements++;
			
			
		}
	}


	int nodeId = 0, elementId = 0, matElId = 0;
	std::ifstream read;
	read.open(FileName.c_str());
	
	double nodecoordX , nodecoordY , nodecoordZ ;
	
	char buf[1024];
	read.getline(buf, 1024);
	read.getline(buf, 1024);
	std::string str(buf);
	int in;
	std::vector <double> co(3,0);
	for(in=0; in<numnodes; in++)
	{ 
		read >> nodeId;
		nodeId--;
		read >> nodecoordX;
		read >> nodecoordY;
		read >> nodecoordZ;
		co[0] = nodecoordX;
		co[1] = nodecoordY;
		co[2] = nodecoordZ;
		TNo no(co);
		malha.insertNode(no);		
	}
	


	{
		read.close();
		read.open(FileName.c_str());
		
		int l , m = numnodes+5;
		for(l=0; l<m; l++)
		{
			read.getline(buf, 1024);
		}
		// Soh funciona para porder 1
		int mat = 1, morder = 1;
		std::vector <int> nodes(4,0);
		int el;
		for(el=0; el<numelements; el++)
		{
			read >> elementId;
			read >> nodes[0]; //node 1
			read >> nodes[1]; //node 2
			read >> nodes[2]; //node 3
			read >> nodes[3]; //node 4
			
			// O GID comeca com 1 na contagem dos nodes, e nao zero como no programa, assim o node 1 na verdade Ž o node 0
			nodes[0]--;
			nodes[1]--;
			nodes[2]--;
			nodes[3]--; 
			
			telemento2d_quad *elem = new telemento2d_quad(mat,morder,nodes);
			malha.insertElement(elem);
		}
		
		//Dirichilet em baixo
		double id = 2, type = 0, contrstiff = 0, contrrhs = 0;
		TMaterialBC *matbc1 = new TMaterialBC(id,type,contrstiff,contrrhs);
		malha.insertMaterial(matbc1);
		
		// Neumann
		id = 3;
		type = 1;
		contrstiff = 0;
		contrrhs = 0;
		TMaterialBC *matbc2 = new TMaterialBC(id,type,contrstiff,contrrhs);
		malha.insertMaterial(matbc2);
		
		// Dirichilet em cima 
		id = 4;
		type = 0;
		contrstiff = 0;
		contrrhs = -2;
		TMaterialBC *matbc3 = new TMaterialBC(id,type,contrstiff,contrrhs);
		malha.insertMaterial(matbc3);


		for(l=0; l<7; l++)
		{
			read.getline(buf, 1024);
		}
		
		nodes.resize(2,0);
		mat = 2;
		for(el=0; el<numelements1d; el++)
		{
			read >> elementId;
			read >> nodes[0]; //node 1
			read >> nodes[1]; //node 2
			
			// O GID comeca com 1 na contagem dos nodes, e nao zero como no programa, assim o node 1 na verdade Ž o node 0
			nodes[0]--;
			nodes[1]--;
			
			TNo &no1 = malha.getNode(nodes[0]);
			TNo	&no2 = malha.getNode(nodes[1]);
			
			if (no1.Co(1) == 0 && no2.Co(1) == 0) 
			{
				mat = 2;
				TElemento1d *elem = 
				new TElemento1d(mat,morder,nodes);
				malha.insertElement(elem);
			}
			else if (no1.Co(1) == 2 && no2.Co(1) == 2)
			{
				mat = 4;
				TElemento1d *elem = new TElemento1d(mat,morder,nodes);
				malha.insertElement(elem);
			}

		}
		
	}
	
}


void ForcingFunction(std::vector<double> &point, std::vector<double> &force)
{
	force.resize(1);
	force[0] = point[0];	
}

void Exact(std::vector<double> &x, double &val, std::vector<double> &deriv)
{
	val = -1 + cos(x[0]) - sin(x[0])/tan(1) + sin(x[0])/sin(1);
	deriv[0] = -cos(x[0])/tan(1) + cos(x[0])/sin(1) - sin(x[0]);
}

void ExactQuad(std::vector <double> &x, double &val, std::vector <double> &deriv)
{
	deriv.resize(2,0);
	val = -x[1]*x[1]/2;
	deriv[0] = 0;
	deriv[1] = -x[1];
	
}

void ReadMeshGidTri(TMalha &malha, std::string &FileName)
{
	int numnodes=-1;
	int numelements=-1;
	int numelements1d=-1;
	int id = 1;
	double k = 1, c = 0, b = 0, f = -1;
	TMaterial2d *mat2d = new TMaterial2d(id,k,c,b,f);
	malha.insertMaterial(mat2d);
	
	//string FileName;
	//FileName = "8andares02.txt";
	
	{
		bool countnodes = false;
		bool countelements = false;
		
		std::ifstream read (FileName.c_str());
		
		while(read)
		{
			char buf[1024];
			read.getline(buf, 1024);
			std::string str(buf);
			if(str == "Coordinates") countnodes = true;
			if(str == "end coordinates") countnodes = false;
			if(countnodes) numnodes++;
			
			if(str == "Elements") countelements = true;
			if(str == "end elements")
			{
				countelements = false;
				while (read) 
				{
					char buf[1024];
					read.getline(buf, 1024);
					std::string str(buf);
					if(str == "Coordinates")
					{
						numnodes--;
						countnodes = true;					
					}
					if(str == "end coordinates") countnodes = false;
					if(countnodes) numnodes++;
					
					if(str == "Elements") countelements = true;
					if(str == "end elements") countelements = false;
					if(countelements) numelements1d++;	
				}
				
			}
			if(countelements) numelements++;
			
			
		}
	}
	
	
	int nodeId = 0, elementId = 0, matElId = 0;
	std::ifstream read;
	read.open(FileName.c_str());
	
	double nodecoordX , nodecoordY , nodecoordZ ;
	
	char buf[1024];
	read.getline(buf, 1024);
	read.getline(buf, 1024);
	std::string str(buf);
	int in;
	std::vector <double> co(3,0);
	for(in=0; in<numnodes; in++)
	{ 
		read >> nodeId;
		nodeId--;
		read >> nodecoordX;
		read >> nodecoordY;
		read >> nodecoordZ;
		co[0] = nodecoordX;
		co[1] = nodecoordY;
		co[2] = nodecoordZ;
		TNo no(co);
		malha.insertNode(no);		
	}
	
	
	
	{
		read.close();
		read.open(FileName.c_str());
		
		int l , m = numnodes+5;
		for(l=0; l<m; l++)
		{
			read.getline(buf, 1024);
		}
		// Soh funciona para porder 1
		int mat = 1, morder = 1;
		std::vector <int> nodes(3,0);
		int el;
		for(el=0; el<numelements; el++)
		{
			read >> elementId;
			read >> nodes[0]; //node 1
			read >> nodes[1]; //node 2
			read >> nodes[2]; //node 3
			
			// O GID comeca com 1 na contagem dos nodes, e nao zero como no programa, assim o node 1 na verdade Ž o node 0
			nodes[0]--;
			nodes[1]--;
			nodes[2]--;
			
			telemento2d_tri *elem = new telemento2d_tri(mat,morder,nodes);
			malha.insertElement(elem);
		}
		
		//Dirichilet em baixo
		double id = 2, type = 0, contrstiff = 0, contrrhs = 0;
		TMaterialBC *matbc1 = new TMaterialBC(id,type,contrstiff,contrrhs);
		malha.insertMaterial(matbc1);
		
		// Neumann
		id = 3;
		type = 1;
		contrstiff = 0;
		contrrhs = 0;
		TMaterialBC *matbc2 = new TMaterialBC(id,type,contrstiff,contrrhs);
		malha.insertMaterial(matbc2);
		
		// Dirichilet em cima 
		id = 4;
		type = 0;
		contrstiff = 0;
		contrrhs = 2;
		TMaterialBC *matbc3 = new TMaterialBC(id,type,contrstiff,contrrhs);
		malha.insertMaterial(matbc3);
		
		
		for(l=0; l<7; l++)
		{
			read.getline(buf, 1024);
		}
		
		nodes.resize(2,0);
		mat = 2;
		for(el=0; el<numelements1d; el++)
		{
			read >> elementId;
			read >> nodes[0]; //node 1
			read >> nodes[1]; //node 2
			
			// O GID comeca com 1 na contagem dos nodes, e nao zero como no programa, assim o node 1 na verdade Ž o node 0
			nodes[0]--;
			nodes[1]--;
			
			TNo &no1 = malha.getNode(nodes[0]);
			TNo	&no2 = malha.getNode(nodes[1]);
			
			if (no1.Co(1) == 0 && no2.Co(1) == 0) 
			{
				mat = 2;
				TElemento1d *elem = 
				new TElemento1d(mat,morder,nodes);
				malha.insertElement(elem);
			}
			else if (no1.Co(1) == 2 && no2.Co(1) == 2)
			{
				mat = 4;
				TElemento1d *elem = new TElemento1d(mat,morder,nodes);
				malha.insertElement(elem);
			}
			
		}
		
	}
	
}
