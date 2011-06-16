//
//  tintegracao.h
//  Fem_sc
//
//  Created by Tiago Costa on 09/05/11.
//  Copyright 2011 Unicamp. All rights reserved.
//
#ifndef TINTEGRACAOH
#define TINTEGRACAOH

#include <vector>

class femsc_exception{};

class TIntegracao
{
    /// ordem da regra de integracao
    int fOrder;
    
    /// coordenadas dos pontos de integracao
    std::vector<double> fPontos;
    
    /// pesos dos pontos de integracao
    std::vector<double> fPesos;
  
public:
    
    TIntegracao(int order, double xmin = -1., double xmax = 1.);
    
    int NPoints();
    
    void Point(int ip, double &x, double &weight);
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