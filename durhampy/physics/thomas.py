#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 15:48:53 2020

@author: mialy94
"""

import math 
import numpy as np

def miniMax(chaine):
    test = True
    try:
        #print(len(np.shape(chaine)))
        test = (len(np.shape(chaine))==1)
        #print(test)
        assert test

    except:
        raise AssertionError('Thomas, t as insere autre chose qu une liste')
    test = True
    try : 
            if(np.shape(chaine)[0]<1):
                test = False
            assert test
    except:
            raise AssertionError('Thomas, cette liste est vide')        

    out=[chaine[0],chaine[0]] #[min,max]
    for el in chaine:
        if(el<out[0]): out[0]=el
        if(el>out[1]): out[1]=el
    
    return out

        
        