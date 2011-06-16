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
#ifndef TMALHA_H
#define TMALHA_H

#include "tnos.h"
//#include "telemento.h"
//#include "tmaterial.h"

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
/**
Implementa uma malha de elementos finitos

@author Philippe R. B. Devloo
*/

class TMaterial;
class TElemento;

class TMalha{
public:
    TMalha();

    ~TMalha();
    
public: //inser��o de n�s, materiais e elementos
    /**
     * Insere um novo n� na malha. O Retorno � o �ndice do elemento no vetor de n�s;
     */
    int insertNode(TNo &node);
    
    /**
     * Insere um material na malha. O Retorno � o �ndice do material no vetor de materiais
     */
    int insertMaterial(TMaterial *mat);
    
    /**
     * Insere um elemento na malha. O Retorno � o �ndice do elemento no vetor de elementos
     */
    int insertElement (TElemento *el);
  

public: //Acesso a leitura dos dados da malha
    /**
     * Acesso ao vetor de nos
     */
    std::vector <TNo> & getNodeVec() ;
    
    /**
     * Acesso ao vetor de Materiais
     */
    std::map <int ,TMaterial *> & getMaterialMap();
    
    /**
     * Acesso ao vetor de elementos
     */
    std::vector <TElemento *> & getElementVec();
    
    /**
     * Retorna um ponteiro para o n� com identificador id
     */
    TNo &getNode(int id);
    
    /**
     * Retorna um ponteiro para o material com identificador id
     */
    TMaterial * getMaterial (int id);
    
    /**
     * Retorna um ponteiro para o elemento com identificador id;
     */
    TElemento * getElement (int id);
    
public://diversos
    /**
     * Imprime a malha
     */
    void Print(std::ostream &out=std::cout);
    
    /**
     * Metodo para ilustrar a criacao de uma malha
     */
     static void main();
  
protected:
    /**
     * Vetor de n�s
     */
    std::vector<TNo> fNodes;
    
    /**
     * Vetor de elementos
     */
    std::vector<TElemento *> fElements;
    
    /**
     * Vetor de materiais
     */
    std::map<int, TMaterial *>  fMaterials;

private:
    
    

};

#endif
