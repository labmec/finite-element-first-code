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
#include "tmalha.h"
#include "telemento.h"
#include "tmaterial.h"
#include "telemento1d.h"
#include "telemento0d.h"
#include "tmaterialbc.h"
#include "tmaterial1d.h"
#include "pzfmatrix.h"

TMalha::TMalha() : fNodes(), fElements(), fMaterials()
{
}


TMalha::~TMalha()
{
  int i;
  for (i=0;i<fElements.size();i++)
  {
    if (fElements[i]) delete fElements[i];
  }
  fElements.resize(0);
  
  std::map<int, TMaterial *>::iterator it;
  for (it=fMaterials.begin(); it != fMaterials.end(); it++)
  {
    if ((*it).second) delete (*it).second;
  }
  fMaterials.clear();
  
  fNodes.resize(0);

}

    
int TMalha::insertNode(TNo &node)
{
  int result = fNodes.size();
  fNodes.push_back(node);
  return result;
}
        
int TMalha::insertMaterial(TMaterial *mat)
{
  if (!mat)
  {
    std::cout << "ERROR - int TMalha::insertMaterial - tentando inserir material nulo"
              << std::endl;
    return -1;
  }
  int id = mat->Id();
  if(fMaterials.find(id) != fMaterials.end())
  {
    std::cout << "ERROR - int TMalha::insertMaterial - inserindo material duplicado"
              << std::endl;
    delete fMaterials[id];
  }
  fMaterials[id]= mat;
  
  return id;
}

int TMalha::insertElement (TElemento *el)
{
  if (!el)
  {
    std::cout << "ERROR - int TMalha::insertElement - tentando inserir elemento nulo"
              << std::endl;
    return -1;
  }
  int result = fElements.size();
  fElements.push_back(el);
  return result;
}
      
std::vector <TNo> & TMalha::getNodeVec()
{
  return fNodes;
}
    
std::map <int, TMaterial *> & TMalha::getMaterialMap()
{
  return fMaterials;
}
    
std::vector <TElemento *> & TMalha::getElementVec()
{
  return fElements;
}
    

TNo &TMalha::getNode(int id)
{
  static TNo gZero;
  if(id < 0 || id >= fNodes.size())
  {
    return gZero;
  }
  return fNodes[id];
}

    
TMaterial * TMalha::getMaterial (int i)
{
    return fMaterials[i];
}
    

TElemento * TMalha::getElement (int i)
{
    return fElements[i];
}
    

void TMalha::Print(std::ostream &out)
{
  int i;
  out << "Impressão da malha" << std::endl;
  out << "Vetor de Nós - número de nós = " << fNodes.size() << std::endl;
  for (i=0;i<fNodes.size();i++) 
  {
    out << "Ind " << i << ' ';
    fNodes[i].Print(out);
  }
  out << "Vetor de Elementos - número de elementos = " << fElements.size() << std::endl;
  for (i=0; i<this->fElements.size(); i++) 
  {
    out << "Ind " << i << ' ';
    fElements[i]->Print(out);
  }
  out << "Vetor de Materiais - número de materiais = " << fMaterials.size() << std::endl;
  std::map<int, TMaterial *>::iterator it;
  for (it=fMaterials.begin(); it != this->fMaterials.end(); it++) 
  {
    (*it).second->Print(out);
  }
}

void TMalha::main()
{
  TMalha malha;
  
  int i;
  
  // criacao dos nos
  int nnodes = 10;
  std::vector<double> co(2,0.);
  for(i=0; i<nnodes; i++)
  {
    co[0] = i*1./(nnodes-1);
    TNo no(co);
    malha.insertNode(no);
  }
  
  // criacao dos elementos uni dimensionais
  std::vector<int> nodeindex(2,0);
  int matindex = 0;
  int order = 1;
  for(i=0; i<nnodes-1; i++)
  {
    nodeindex[0] = i;
    nodeindex[1] = i+1;
    TElemento1d *el = new TElemento1d(matindex,order,nodeindex);
    malha.insertElement(el);
  }
  // criacao do primeiro elemento pontual
  nodeindex.resize(1);
  nodeindex[0] = 0;
  matindex = 1;
  TElemento0d *el = new TElemento0d(matindex,0,nodeindex);
  malha.insertElement(el);
  // criacao do segundo elemento pontual
  matindex = 2;
  nodeindex[0] = nnodes-1;
  el = new TElemento0d(matindex,0,nodeindex);
  malha.insertElement(el);
  
  double K(1.),C(0.),B(1.),F(2.);
  TMaterial1d *mat1d = new TMaterial1d(0,K,C,B,F);
  malha.insertMaterial(mat1d);
  
  int bctype(0);
  matindex = 1;
  double contrstiff(0.),contrrhs(0.);
  TMaterialBC *matbc0 = new TMaterialBC(matindex,bctype,contrstiff,contrrhs);
  malha.insertMaterial(matbc0);
  bctype = 1;
  matindex = 2;
  contrrhs = 1.;
  TMaterialBC *matbc1 = new TMaterialBC(matindex,bctype,contrstiff,contrrhs);
  malha.insertMaterial(matbc1);

  // impressao da malha
  malha.Print();
  
  int elem;
  TPZFMatrix stiff,rhs;
  for(elem=0; elem< malha.fElements.size(); elem++)
  {
    malha.fElements[elem]->CalcStiff(malha,stiff,rhs);
    stiff.Print("Matriz de rigidez");
    rhs.Print("Vetor de carga");
  }

}
