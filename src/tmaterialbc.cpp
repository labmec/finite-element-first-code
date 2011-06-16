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
#include "tmaterialbc.h"
#include "pzfmatrix.h"

static double BigNumber = 1.e12;

TMaterialBC::    TMaterialBC(int id, int bctype, double contrstiff, double contrrhs) : TMaterial(id), fBCType(bctype),
  fContrStiff(contrstiff), fContrRhs(contrrhs)
{
}


TMaterialBC::~TMaterialBC()
{
}


void TMaterialBC::Print(std::ostream& out) const
{
    TMaterial::Print(out);
    out << "BC type " << fBCType;
    out << " BC values stiff " << fContrStiff << " rhs " << fContrRhs << std::endl;
}

void TMaterialBC::Contribute (double  weight,
                           std::vector<double> & philVal,
                           TPZFMatrix & dphi,TPZFMatrix & elementK,
                           TPZFMatrix & elementF) const 
  {
    int i, nshape, j;
    nshape = philVal.size();
    switch(fBCType) 
    {
      case 0: // DIRICHLET
      {
        for(i=0; i<nshape; i++)
        {
          elementF(i,0) += philVal[i]*weight*BigNumber*fContrRhs;
          for(j=0; j<nshape; j++)
          {
            elementK(i,j) += philVal[i]*philVal[j]*weight*BigNumber;
          }
        }
      }
      break;
      case 1: // Neumann
      {
        for(i=0; i<nshape; i++)
        {
          elementF(i,0) += weight*philVal[i]*fContrRhs;
        }
      }
      break;
      case 2: // Mixed BC
      {
        for(i=0; i<nshape; i++)
        {
          elementF(i,0) += philVal[i]*weight*fContrRhs;
          for(j=0; j<nshape; j++)
          {
            elementK(i,j) += philVal[i]*philVal[j]*weight*fContrStiff;
          }
        }
      }
      break;
      default:
        std::cout << __PRETTY_FUNCTION__ << " wrong boundary condition type " << fBCType << std::endl;
        break;
    }
  }

