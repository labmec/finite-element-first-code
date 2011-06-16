/***************************************************************************
 *   Copyright (C) 2005 by Philippe R. B. Devloo                           *
 *   phil@fec.unicamp.br                                                   *
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
#include "tanalysis.h"
#include "pzfmatrix.h"

TAnalysis::TAnalysis(TMalha *malha) : fMalha(malha), fSolution()
{
}


TAnalysis::~TAnalysis()
{
}

    // Builds the global stiffness matrix and right hand side
void TAnalysis::Run()
{
   if(!fMalha) return;
   int neq = fMalha->getNodeVec().size();
   TPZFMatrix globstiff(neq,neq,0.);
   TPZFMatrix globrhs(neq,1,0.);
   Assemble(globstiff,globrhs);
   fSolution = globrhs;
   globstiff.Solve_LU(&fSolution);
   
   fSolution.Print("The solution ");
   
}

    // assembly method
void TAnalysis::Assemble(TPZMatrix &stiff, TPZFMatrix &rhs)
{
   int nelem = fMalha->getElementVec().size();
   int iel;
   TPZFMatrix locstiff,locrhs;
   for(iel=0; iel<nelem; iel++)
   {
     fMalha->getElement(iel)->CalcStiff(*fMalha,locstiff,locrhs);
     const std::vector<int> &nodes = fMalha->getElement(iel)->getNodeVec();
     int nnodes = nodes.size();
     int in,jn;
     for(in=0; in<nnodes; in++)
     {
     	rhs(nodes[in],0) += locrhs(in,0);
	for(jn=0; jn<nnodes; jn++)
	{
		stiff(nodes[in],nodes[jn]) += locstiff(in,jn);
	}
     }
   }
}


    // computes the error in energy and l2 norm
void TAnalysis::Error(void (exact) (std::vector<double> &x, double &val, std::vector<double> &deriv), double &energy, double &l2)
{
   int nelem = fMalha->getElementVec().size();
   int iel;
   energy = 0.;
   l2 = 0.;
   if(!fSolution.Rows()) return;

   for(iel=0; iel<nelem; iel++)
   {
     double energyel = 0.;
     double l2el = 0.;
     fMalha->getElement(iel)->Error(fSolution,*fMalha,exact,energyel,l2el);
     energy += energyel*energyel;
     l2 += l2el*l2el;
   }
   energy = sqrt(energy);
   l2 = sqrt(l2);
}
