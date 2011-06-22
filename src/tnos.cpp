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
#include "tnos.h"

TNo::TNo()
{
  fCoordenadas[0] = 0.;
  fCoordenadas[1] = 0.;
}


TNo::~TNo()
{
}

TNo::TNo (std::vector<double> &coordenadas)
{
  fCoordenadas[0] = 0.;
  fCoordenadas[1] = 0.;
  int i = 0;
  for(i=0; i<coordenadas.size() && i<2;i++)
  {
    fCoordenadas[i]=coordenadas[i];
  }
  
/*  
  std::vector<double>::iterator it;
  for(it=coordenadas.begin(); it!=coordenadas.end() && i<2; it++,i++)
  {
    fCoordenadas[i] = (*it);
  }
*/  
}
    
TNo::TNo (const TNo &copy)
{ 
  int i;
  for(i=0; i<2; i++)
  {
    fCoordenadas[i] = copy.fCoordenadas[i];
  }
}


void TNo::setData(std::vector<double> &coordenadas)
{
  fCoordenadas[0] = 0.;
  fCoordenadas[1] = 0.;
  int i = 0;
  std::vector<double>::iterator it;
  for(it=coordenadas.begin(); it!=coordenadas.end() && i<2; it++,i++)
  {
    fCoordenadas[i] = (*it);
  }  
}

double &TNo::Co(int i)
{
  static double gZero = 0;
  if (i< 0 || i>1) 
  {
    std::cout << "ERRO - double TNo::getCoordenada(int i)\n"
              << "\tIndice de coordenada " << i << " fora das dimensões do nó!"
              << "\tRetornando valor nulo"
              << std::endl;
    gZero = 0.;
    return gZero;
  }
  return fCoordenadas[i];
}
       
    

void TNo::Print(std::ostream &out)
{
  out << "Nó: coordenadas ( ";
  int i;
  out << fCoordenadas[0];
  for (i=1;i<2;i++)
  {
     out << " , " << fCoordenadas[i];
  }
  
  out  << " );" << std::endl;
}

int TNo::main()
{
  std::cout << "Testando a classe no!" << std::endl;
  
  std::cout << "\tcriando vetor de 2 double (0.,0.)" << std::endl;
  std::vector <double> coord (2,0.);
  
  std::cout << "\tcriando o noh 0 com identificador 100 coordenadas do vetor acima..." <<std::endl;
  TNo no0 (coord);
  
  std::cout << "\to resultado eh: " << std::endl;
  no0.Print();
  
  std::cout << std::endl;
  
  std::cout << "Criando um vetor de nos. No loop sao definidos os dados dos nohs" << std::endl;
  std::cout << "Entre com numero de linhas e de colunas :";
  std::cout.flush();
  int linha,coluna;
  std::cin >> linha >> coluna;
  
  std::vector <TNo> nos (linha*coluna);
  int i,j;
  for (i=0;i<linha;i++)
  {
    coord[0] = (double) i;
    for (j=0;j<coluna;j++)
    {
      coord[1] = (double) j; 
      nos[i*coluna + j].setData(coord);
    }
  }
  
  std::cout << "os nos criados foram:" << std::endl;
  for (i=0;i<nos.size();i++)
  {
    nos[i].Print();
  }
}
