//
//  TIntegracao2d_tri.h
//  Fem_sc
//
//  Created by Tiago Costa on 09/05/11.
//  Copyright 2011 Unicamp. All rights reserved.
//
#ifndef TIntegracao2d_triH
#define TIntegracao2d_triH

#include <vector>

class femsc_exception{};

class TIntegracao2d_tri
{
    /// ordem da regra de integracao
    int fOrder;
    
    /// coordenadas dos pontos de integracao
    std::vector<std::pair<double,double> > fPontos;
    
    /// pesos dos pontos de integracao
    std::vector<double> fPesos;
  
public:
    
    TIntegracao2d_tri(int order);
    
    int NPoints();
    
    void Point(int ip, std::vector<double> &x, double &weight);
};
/*
class TIntegrategauher
{
    /// ordem da regra de integracao
    int fOrder;
    
    /// coordenadas dos pontos de integracao
    std::vector<double> fPontos;
    
    /// pesos dos pontos de integracao
    std::vector<double> fPesos;
    
public:
    
    TIntegrategauher(int order);
    
    int NPoints();
    
    double Point(int ip, std::vector<double> &loc,double &weight);
};

class TIntegratelobatto
{
    /// ordem da regra de integracao
    int fOrder;
    
    /// coordenadas dos pontos de integracao
    std::vector<double> fPontos;
    
    /// pesos dos pontos de integracao
    std::vector<double> fPesos;
    
    /// coordenadas dos pontos de integracao
    std::vector<double> fa;
    
    /// pesos dos pontos de integracao
    std::vector<double> fb;
    
public:
    
    TIntegratelobatto(int order);
    
    int NPoints();
    
    double Point(int ip, std::vector<double> &loc,double &weight);
};*/

#endif