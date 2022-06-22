//
//  CQILTE.h
//  LTE_PHY
//
//  Created by Vaishnavi  on 25/12/19.
//  Copyright © 2019 none. All rights reserved.
//

#ifndef CQILTE_h
#define CQILTE_h
#include "Modulation.hpp"
struct CQI{
    unsigned int id;
    Modulationscheme modulation;
    double coderate;
    double efficiency;
};


CQI CQI_Table[]=
{   {0,BPSK(),0,0},
    {1,QPSK(),78.0/1024.0,.1523},
    {2,QPSK(),120.0/1024.0,.2344},
    {3,QPSK(),193.0/1024.0,.3770},
    {4,QPSK(),308.0/1024.0,.6016},
    {5,QPSK(),449.0/1024.0,.8770},
    {6,QPSK(),602.0/1024.0,1.1758},
    {7,QAM16(),378.0/1024.0,1.4766},
    {8,QAM16(),490.0/1024.0,1.9141},
    {9,QAM16(),616.0/1024.0,2.4063},
    {10,QAM64(),466.0/1024.0,2.7305},
    {11,QAM64(),567.0/1024.0,3.3223},
    {12,QAM64(),666.0/1024.0,3.9023},
    {13,QAM64(),772.0/1024.0,4.5234},
    {14,QAM64(),873.0/1024.0,5.1152},
    {15,QAM64(),948.0/1024.0,5.5547}
};
#endif /* CQILTE_h */
