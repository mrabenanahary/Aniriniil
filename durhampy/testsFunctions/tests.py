#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 09:14:38 2019

@author: mialy94
=================================================================
= Auteur : Mialy RABENANAHARY
= 
= tests.py
= 
= Module contenant toutes les fonctions réalisant dans le main.py
= des tests, notamment sur les données et modèles manipulés
= dans ce script Python
=
=================================================================

"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
from scipy.optimize import curve_fit
import csv
import os
from astropy.io import ascii
from matplotlib.ticker import AutoMinorLocator
import matplotlib.text as mtext
#from durhampy.data.globals import *

#variables globales
#global input_dir

"================FONCTIONS======================"

def verifieModele(actual_dir,                                                   # nom du répertoire parent direct à celui de travail actuel 
                  inputDir='cj-n4-b1',                                          # choix du modèle cx-ny-bz traité
                  vit=10,                                                       # choix de la vitesse de choc traité dans la nomenclature v"a"-t"b"
                  ag=100,                                                       # choix de l'âge du choc traité dans la nomenclature v"a"-t"b"
                  nom_de_fonction = 'extraire_temperature',                     # nom de la fonction mère d'où est appelée verifieModele()
                  return_info_mhd = False,                                      # booléen switchant le retour ou non par la fonction verifieModele()
                                                                                # de la liste des paramètres du modèle actuel tels qu'extraits de
                                                                                # chaque fichier d'input 'info_mhd.out' sorti par le code de Paris-Durham
                  simpilified_log10nH_and_Bbeta = True,                         # bool à activer si la notation des dossiers cx-ny-bz ont ny et bz écrits sous forme d'entiers sans décimales (y et z étant forcément entiers dans le modèle)
                  simpilified_vs = True,                                        # bool comme simpilified_log10nH_and_Bbeta mais pour la vitesse de choc v_s
                  simpilified_timeJ = True,                                     # bool comme simpilified_log10nH_and_Bbeta mais pour l'âge du choc t
                  param_tester = {'shock parameters':                           # dictionnaires des paramètres à inclure dans la série des 3 tests effectués par cette fonction et impliquant l'extraction de ces paramètres à partir de 'info_mhd.out'
                    ['shocktype',                                               # par défaut, param_tester = {'shock parameters': ['shocktype, 'nH(cm-3)', 'B(microGauss)', 'Vs(km.s-1)'], 'integration parameters': ['timeJ']}
                     'nH(cm-3)',                                                # s'il vous faut rajouter des variables à tester, il faudra les renseigner dans l'input de cette fonction via les dictionnaires valeurs_param_tester et valeurs_param_tester
                     'B(microGauss)',                                           
                     'Vs(km.s-1)'],                                             
                                 'integration parameters':
                    ['timeJ'],
                    'environmental parameters': ['Zeta(s-1)']
                    },
                  valeurs_param_tester = {'input files':{},
                                          'shock parameters': {},
                                          'environmental parameters': {},
                                          'grains parameters': {},
                                          'excitation & cooling': {},
                                          'integration parameters': {},
                                          'DVODE parameters': {},
                                          'outputs parameters': {},
                                          'H2 molecule': {},
                                          'elemental abundances (gas + mantles + PAH)': {}
                                          },                                  # même si ce dictionnaire est vide, verifieModele() effectuera toujours en test 2 et 3 les tests sur cx-ny-bz et v'a'-t'b'   
                  print_tests=False):                                         # si True, affiche dans la console des indications détaillées et longues des résultats des tests, sinon si False, se contente d'afficher 'Test n°{1,2,3} réussi/échoué' à chaque modèle
    """
    
        \n Fonction qui effectue une série de 3 tests dans le répertoire de travail et modèle actuels (i.e. type de choc,
        densité pré-choc, champ magnétique, vitesse de choc, âge de choc, etc...) :
            
            - 1er Test : verifieModele() teste que dans le modèle de travail actuel, les paramètres
            de type de choc cx, de densité pré-choc ny, de champ magnétique bz, de vitesse de choc v'a', d' âge de choc
            t'b' qui sont envoyés en input de la fonction verifieModele(**kwargs) (via les arguments inputDir, vit, et ag)
            correspondent bien au dossier courant où travail la fonction mère (celle qui fait appel à verifieModele()), dont
            les deux derniers niveaux d'arborescence doivent forcément être écrites sous le format ../../cx-ny-bz/v'a'-t'b'
            ce test donne lieu à une levée d'exception sur un booléen testé et en cas d'échec du test, retourne un message 
            d'erreur détaillé
            
            - 2ème Test : verifieModele() teste que dans le modèle de travail actuel, les paramètres
            (tels que le type de choc cx, le densité pré-choc ny, le champ magnétique bz, le vitesse de choc v'a', l' âge de choc
            t'b') qui sont extraits du fichier 'info_mhd.out' (sorti à la fin de chaque run du code de choc de Paris-Durham )
            et qui sont choisis par l'utilisateur (via le dictionnaire en input param_tester)
            correspondent bien aux paramètres passés en input de la fonction verifieModele() (via les dictionnaires en input param_tester
            et valeurs_param_tester) via la fonction mère qui l'appelle (par défaut, si ces dictionnaires ne sont pas re-précisés en input,
            seuls les paramètres cx, ny, bz, v'a' et t'b' sont testés)
            
            - 2ème Test : verifieModele() teste que dans le modèle de travail actuel, les paramètres
            (tels que le type de choc cx, le densité pré-choc ny, le champ magnétique bz, le vitesse de choc v'a', l' âge de choc
            t'b') qui sont extraits du fichier 'info_mhd.out' (sorti à la fin de chaque run du code de choc de Paris-Durham )
            et qui sont choisis par l'utilisateur (via le dictionnaire en input param_tester)
            correspondent bien aux paramètres passés en input de la fonction verifieModele() (via les dictionnaires en input param_tester
            et valeurs_param_tester) via la fonction mère qui l'appelle (par défaut, si ces dictionnaires ne sont pas re-précisés en input,
            seuls les paramètres cx, ny, bz, v'a' et t'b' sont testés). 
            Ce test donne lieu à une levée d'exception sur un booléen testé et en cas d'échec du test, retourne un message 
            d'erreur détaillé
            
            - 3ème test: Même test que le test 2, si ce n'est que les paramètres extraits de 'info_mhd.out' sont ici
            cette fois comparées aux deux derniers niveaux d'arborescence du répertoire de travail actuel : donc, ici
            seuls les paramètres cx, ny, bz, v'a' et t'b' dans 'info_mhd.out'  sont testés. 
            Ce test donne lieu à une levée d'exception sur un booléen testé et en cas d'échec du test, retourne un message 
            d'erreur détaillé
            
            
        Grâce à cette série de 3 tests, verifieModele() permets de vérifier sans ambiguité (en cas d'échec et de succès des tests)
        qu'au moment de son appel, la fonction mère travaille bien dans le bon dossier et sur le bon modèle de choc 
        (avec les paramètres consistants et concordants entre le répertoire
        de travail courant, les paramètres input lors du run du code de choc de Paris-Durham, ainsi que les paramètres Python de modèle avec
        lesquels travaille la fonction mère).
            
        (Dernière vérification de cette fonction le 02/12/2019)
            
    Utilisation:
            
        * si  return_info_mhd == False: 
            verifieModele(**kwargs)
        * si  return_info_mhd == True: 
            output = verifieModele(**kwargs)
            
    Paramètres d'entrée:
    ===================
        Paramètres obligatoires:
        ------------------------
        
        actual_dir : nom du répertoire parent direct (dernier niveau d'arborescence) à celui de travail courant 
        
        Paramètres facultatifs, **kwargs:
        ---------------------------------
            
        inputDir: nom du dossier de format cx-ny-bz (normalement au niveau du deuxième dernier niveau d'arborescence de travail courant)
        renseigné depuis la fonction mère dans verifieModele(). Doit être un string de format cx-ny-bz et vaut par défaut 'cj-n4-b1'
        
        vit: vitesse de choc traitée par la fonction à partir 
            d'un dossier obligatoirement sous le format de nom v"a"-t"b". vit 
            doit être un nombre réel et vaut par défaut 10
        
        ag: âge du choc traité par la fonction à partir 
            des dossiers obligatoirement sous le format de nom v"a"-t"b". vit_dir 
            ne peut prendre qu'une valeur  réelle et vaut par défaut 100.
        
        nom_de_fonction : string indiquant le nom de la fonction mère d'où est appelée verifieModele(). 
            Par défaut, vaut 'extraire temperature'
        
        return_info_mhd : booléen switchant le retour ou non par la fonction verifieModele() de la liste des 
            paramètres du modèle actuel tels qu'extraits de chaque fichier d'input 'info_mhd.out' sorti par le code de Paris-Durham.
            Par défaut, vaut False.
        
        simpilified_log10nH_and_Bbeta : booléen à activer si la notation des dossiers cx-ny-bz ont ny et bz écrits 
            sous forme d'entiers sans décimales (y et z étant forcément entiers dans le modèle). Par défaut, vaut True.
            
        simpilified_vs : booléen comme simpilified_log10nH_and_Bbeta mais pour la vitesse de choc v_s. Par défaut True.
        
        simpilified_timeJ : booléen comme simpilified_log10nH_and_Bbeta mais pour l'âge du choc t. Par défaut True.
        
        param_tester : dictionnaire de format { nom_de_categorie : [nom_de_parametre] } } des noms de paramètres à inclure dans la série 
            des 3 tests effectués par cette fonction et impliquant l'extraction de ces paramètres à partir de 'info_mhd.out'.
            Par défaut, param_tester = {'shock parameters': ['shocktype, 'nH(cm-3)', 'B(microGauss)', 'Vs(km.s-1)'], 'integration parameters': ['timeJ']}
            S'il vous faut rajouter des variables à tester, il faudra les renseigner dans l'input de cette fonction 
            via les dictionnaires valeurs_param_tester et valeurs_param_tester.

        valeurs_param_tester : dictionnaire de format { nom_de_categorie : {nom_du_parametre : valeur_du_parametre } } } des noms et valeurs des paramètres à inclure dans la série 
            des 3 tests effectués par cette fonction et impliquant l'extraction de ces paramètres à partir de 'info_mhd.out'.
            Par défaut, valeurs_param_tester = {'input files':{},
                                          'shock parameters': {},
                                          'environmental parameters': {},
                                          'grains parameters': {},
                                          'excitation & cooling': {},
                                          'integration parameters': {},
                                          'DVODE parameters': {},
                                          'outputs parameters': {},
                                          'H2 molecule': {},
                                          'elemental abundances (gas + mantles + PAH)': {}
                                          }, 
            mais la fonction verifieModele() inclut toujours systématiquement dans les tests 2 et 3 les paramètres 
             {'shock parameters': ['shocktype, 'nH(cm-3)', 'B(microGauss)', 'Vs(km.s-1)'], 'integration parameters': ['timeJ']}.
            S'il vous faut rajouter des variables à tester, il faudra les renseigner dans l'input de cette fonction 
            via les dictionnaires valeurs_param_tester et valeurs_param_tester.
        
        print_tests : booléen. Si True, la fonction verifieModele() affiche dans la console des indications détaillées et  longues des résultats
            des tests. Si False, verifieModele() se contente d'afficher 'Test n°{1,2,3} réussi/échoué' pour le modèle, dossier et info_mhd.out courants donnés
            
            
    Paramètres de sortie:
    ===================
        
        parameters_Dict_file: SEULEMENT SI return_info_mhd == True,
            verifieModele() retourne en sortie un dictionnaire listant tous les paramètres extraits du fichier 'info_mhd.out'
            du modèle et dossier courants, classés comme dans ce fichier (nom de la catégorie de paramètres, puis nom du paramètre, puis sa valeur).
            Il s'agit donc des paramètres renseignés dans le fichier input 'input_mhd.in' utilisé par le code de choc de Paris Durham
            pour lancer le run de ce modèle.
            
            
    e.g. d'utilisation : 
    ===================
            
            - input_dir = 'cj-n4-b1'
            - modele = 'v10-t100', vitesse=10, age = 100
            - nomDeFonction = 'extraire_raie'
            - print_tests_results = True
            - info_mhd_H2 = verifieModele(actual_dir=modele,inputDir=input_dir,vit=vitesse,ag=age, nom_de_fonction=nomDeFonction, return_info_mhd=True, print_tests=print_tests_results, 
                    param_tester =      {'shock parameters':                           
                                        ['shocktype',                                  
                                         'nH(cm-3)',                                   
                                         'B(microGauss)',                                           
                                         'Vs(km.s-1)'],                                             
                                                     'integration parameters':
                                        ['timeJ'],
                                        'environmental parameters': ['Zeta(s-1)']
                                        },
                    valeurs_param_tester = {'input files':{},
                                          'shock parameters': {},
                                          'environmental parameters': {'Zeta(s-1)':5.00E-17},
                                          'grains parameters': {},
                                          'excitation & cooling': {},
                                          'integration parameters': {},
                                          'DVODE parameters': {},
                                          'outputs parameters': {},
                                          'H2 molecule': {},
                                          'elemental abundances (gas + mantles + PAH)': {}
                                          }          )
            - verifieModele(actual_dir=modele,inputDir=input_dir,vit=vitesse,ag=age, nom_de_fonction=nomDeFonction, print_tests=False, 
                    param_tester =      {'shock parameters':                           
                                        ['shocktype',                                  
                                         'nH(cm-3)',                                   
                                         'B(microGauss)',                                           
                                         'Vs(km.s-1)'],                                             
                                                     'integration parameters':
                                        ['timeJ'],
                                        'environmental parameters': ['Zeta(s-1)']
                                        },
                    valeurs_param_tester = {'input files':{},
                                          'shock parameters': {},
                                          'environmental parameters': {'Zeta(s-1)':5.00E-17},
                                          'grains parameters': {},
                                          'excitation & cooling': {},
                                          'integration parameters': {},
                                          'DVODE parameters': {},
                                          'outputs parameters': {},
                                          'H2 molecule': {},
                                          'elemental abundances (gas + mantles + PAH)': {}
                                          }          )
    """     
    
    """Test n°1: Tester que les paramètres {cx-ny-bz; v_s; t} input dans la fonction de test (i.e.
    renseignés au niveau de la fonction-mère) correspond bien au dossier actuel où travaille
    cette fonction mère ==> 1 booléen propre + message d'erreur"""
    
    try:
        #booléens de test:
        va_tb_ok = True #tout va bien
        cx_ny_bz_ok = True #tout va bien !
        cx_ny_bz_va_tb_ok = True #Tout va bien! 
        
        #extraction des chemins
        chemin      = os.getcwd()
        chemin_list = chemin.split('/')
        
        #1ères variables qui seront testées et qui nous sont retournées par os.getcwd()
        va_tb       = chemin_list[-1]                                           # nom du premier sous-dossier du niveau inférieur au répertoire actuel : si tout est bien fait, doit être de format v'a'-t'b'
        cx_ny_bz    = chemin_list[-2]                                           # nom du deuxième sous-dossier du niveau inférieur au répertoire actuel : si tout est bien fait, doit être de format cx-ny-bz
        
        #2ème série de variables qui seront testées et qui nous sont récupérées en input de verifieModele()
        va_tb_tester = 'v{0:}-t{1:}'.format(vit,ag)
        
        error_msg = []                                                          # initialisation de la liste des potentiels messages d'erreurs à afficher à l'issu de chacun des 3 tests
        if(not print_tests) : print(">>>>Fonction mère : ",nom_de_fonction,"() : modele => ", inputDir, '/', va_tb_tester)    #si print_tests == False : affiche dans quel modèles sont réalisés chaque test d'après les inputs de verifieModele()
        if(print_tests) : print("\n>>>>>>>>>Test n°1: Tester que les paramètre {cx-ny-bz; v_s; t} input dans la fonction de test (i.e. renseignées au niveua de la fonction-mère) correspond bien au dossier actuel où travaille cette fonction mère ==> 1 booléen propre + message d'erreur<<<<<<<<<<<<<<\n""")
        
        #test sur v'a'-t'b'
        try:
            assert va_tb_tester == va_tb
        except AssertionError:
            va_tb_ok = False 
            error_msg.append("1]Attention, dans la fonction {3:}(), l'utilisateur entre un modele \n nommé {0:}, alors que le dossier de travail à cet endroit est \n ../{1:} ; le chemin complet est : {2:}\n".format(va_tb_tester,va_tb,chemin,nom_de_fonction))
        
        #test sur cx-ny-bz
        try:
            assert cx_ny_bz == inputDir
        except AssertionError:
            cx_ny_bz_ok = False 
            error_msg.append("2]Attention, dans la fonction {3:}(), l'utilisateur entre un modele \n nommé {0:}, alors que le dossier de travail au niveau supérieur est \n ../../{1:} ; le chemin complet est : {2:}\n".format(inputDir,cx_ny_bz,chemin,nom_de_fonction))
        
        #test global
        cx_ny_bz_va_tb_ok = (va_tb_ok and cx_ny_bz_ok)
        assert cx_ny_bz_va_tb_ok == True
    
    except AssertionError:
        if(print_tests):
            for messages in error_msg:
                print(messages)
            raise AssertionError(">>FAIRE ATTENTION AU MESSAGE D'ERREUR CI-DESSUS |^|\n")
        else:
            raise AssertionError('Test 1 échoué\n')
    else:
        if(print_tests): print("3]Dans la fonction {3:}(), dans le modele \n nommé {0:}, aucun problème de correspondance entre les paramètres input dans la fonction mère \n et le chemin de travail actuel \n ../../{1:} ; le chemin complet est : {2:}\n".format(inputDir+'/'+va_tb_tester,cx_ny_bz+'/'+va_tb,chemin,nom_de_fonction))
        else: print('Test 1 réussi\n')




    """Test n°2: Tester que les paramètres {cx-ny-bz; v_s; t} du code de choc retrouvés dans le finchier input 'info_mhd.out' 
    correspond bien aux paramètres input dans la fonction de test (~test 1) ==> 1 booléen propre + message d'erreur"""
    
    # extraction des paramètres issus du fichier 'info_mhd.out'
    parameters_Dict_file={'date':'',
                     'time':'',
                     'input files':{},
                     'shock parameters':{},
                     'environmental parameters':{},                        
                     'grains parameters':{},  
                     'integration parameters':{},  
                     'excitation & cooling':{},
                     'DVODE parameters':{},  
                     'outputs parameters':{}, 
                     'H2 molecule':{}, 
                     'elemental abundances (gas + mantles + PAH)':{}
        }
    with open('info_mhd.out','r') as input_mhd_file:
        #extraction de : la date
        line_loop = input_mhd_file.readline()
        dateTemp = ((line_loop.split('date :')[1]).split('time :')[0]).split()
        
        
        date = ''
        for elements in dateTemp:
            date+=elements
        parameters_Dict_file['date']=date
        
        #extraction de : l'heure
        heureTemp = ((line_loop.split('date :')[1]).split('time :')[1]).split()
        heure = ''
        for elements in heureTemp:
            heure+=elements
        parameters_Dict_file['time']=heure
        
        #extraction des paramètres de : input files
        line_loop = input_mhd_file.readline()                                       #ligne 2
        line_loop = input_mhd_file.readline()                                       #ligne 3
        
        for i in range(0,5):
            line_loop = input_mhd_file.readline()
            key = (line_loop.split(':'))[0].replace(' ','')
            value_of_key = (line_loop.split(':'))[1].replace('\n','')
            (parameters_Dict_file['input files'])[key]=value_of_key
        
        
        #extraction des paramètres de : shock parameters
        line_loop = input_mhd_file.readline()                                       #ligne 9
        line_loop = input_mhd_file.readline()                                       #ligne 10
        
        for i in range(0,8):
            line_loop = input_mhd_file.readline()
            key = (line_loop.split(':'))[0].replace(' ','')
            value_of_key = (line_loop.split(':'))[1].replace('\n','')
            if(i==0) : value_of_key = value_of_key.replace(' ','')
            else: value_of_key = np.double(value_of_key)
            (parameters_Dict_file['shock parameters'])[key]=value_of_key
        
        #extraction des paramètres de : environmental parameters
        line_loop = input_mhd_file.readline()                                       #ligne 19
        line_loop = input_mhd_file.readline()                                       #ligne 20
        
        for i in range(0,11):
            line_loop = input_mhd_file.readline()
            key = (line_loop.split(':'))[0].replace(' ','')
            value_of_key = (line_loop.split(':'))[1].replace('\n','')
            value_of_key = np.double(value_of_key)
            (parameters_Dict_file['environmental parameters'])[key]=value_of_key
            
        #extraction des paramètres de : grains parameters
        line_loop = input_mhd_file.readline()                                       #ligne 32
        line_loop = input_mhd_file.readline()                                       #ligne 33
        
        for i in range(0,7):
            line_loop = input_mhd_file.readline()
            key = (line_loop.split(':'))[0].replace(' ','')
            value_of_key = (line_loop.split(':'))[1].replace('\n','')
            value_of_key = np.double(value_of_key)
            (parameters_Dict_file['grains parameters'])[key]=value_of_key
            
        #extraction des paramètres de : excitation & cooling
        line_loop = input_mhd_file.readline()                                       #ligne 41
        line_loop = input_mhd_file.readline()                                       #ligne 42
        
        for i in range(0,9):
            line_loop = input_mhd_file.readline()
            key = (line_loop.split(':'))[0].replace(' ','')
            value_of_key = (line_loop.split(':'))[1].replace('\n','')
            if(i==4): value_of_key = value_of_key
            else: value_of_key = np.double(value_of_key)
            (parameters_Dict_file['excitation & cooling'])[key]=value_of_key
            
        #extraction des paramètres de : integration parameters
        line_loop = input_mhd_file.readline()                                       #ligne 52
        line_loop = input_mhd_file.readline()                                       #ligne 53
        
        for i in range(0,4):
            line_loop = input_mhd_file.readline()
            key = (line_loop.split(':'))[0].replace(' ','')
            value_of_key = (line_loop.split(':'))[1].replace('\n','')
            value_of_key = np.double(value_of_key)
            (parameters_Dict_file['integration parameters'])[key]=value_of_key
            
        #extraction des paramètres de : DVODE parameters 
        line_loop = input_mhd_file.readline()                                       #ligne 58
        line_loop = input_mhd_file.readline()                                       #ligne 59
        
        for i in range(0,9):
            line_loop = input_mhd_file.readline()
            key = (line_loop.split(':'))[0].replace(' ','')
            value_of_key = (line_loop.split(':'))[1].replace('\n','')
            value_of_key = np.double(value_of_key)
            (parameters_Dict_file['DVODE parameters'])[key]=value_of_key
            
        #extraction des paramètres de : outputs parameters
        line_loop = input_mhd_file.readline()                                       #ligne 69
        line_loop = input_mhd_file.readline()                                       #ligne 70
        
        for i in range(0,9):
            line_loop = input_mhd_file.readline()
            key = (line_loop.split(':'))[0].replace(' ','')
            value_of_key = (line_loop.split(':'))[1].replace('\n','')
            if(i in [5,6,7,8]): value_of_key = value_of_key
            else: value_of_key = np.double(value_of_key)
            (parameters_Dict_file['outputs parameters'])[key]=value_of_key
            
        #extraction des paramètres de : H2 molecule 
        line_loop = input_mhd_file.readline()                                       #ligne 80
        line_loop = input_mhd_file.readline()                                       #ligne 81
        
        for i in range(0,8):
            line_loop = input_mhd_file.readline()
            key = (line_loop.split(':'))[0].replace(' ','')
            value_of_key = (line_loop.split(':'))[1].replace('\n','')
            if(i in [4]): value_of_key = value_of_key
            else: 
                if(i in [1,7]): value_of_key = np.double(value_of_key.split()[0]) #print(value_of_key.split()[0])
                else : value_of_key = np.double(value_of_key)
            (parameters_Dict_file['H2 molecule'])[key]=value_of_key
            
        #extraction des paramètres de : elemental abundances (gas + mantles + PAH)
        line_loop = input_mhd_file.readline()                                       #ligne 90
        line_loop = input_mhd_file.readline()                                       #ligne 91
        
        for i in range(0,12):
            line_loop = input_mhd_file.readline()
            key = (line_loop.split(':'))[0].replace(' ','')
            value_of_key = (line_loop.split(':'))[1].replace('\n','')
            value_of_key = np.double(value_of_key.split()[0]) #print(value_of_key.split()[0])
            (parameters_Dict_file['elemental abundances (gas + mantles + PAH)'])[key]=value_of_key
            
    # extraction des paramètres issus de la fonction mère passés en input de cette fonction-ci de test
    
    #type de choc
    type_de_choc = ''
    type_de_choc_Temp = inputDir.split('-')[0][0]
    if(type_de_choc_Temp=='c'): type_de_choc_Temp = type_de_choc_Temp.replace('c','C')
    if(type_de_choc_Temp=='j'): type_de_choc_Temp = type_de_choc_Temp.replace('j','J')
    if(np.abs(parameters_Dict_file['integration parameters']['timeJ']-9.99E+99)/9.99E+99>0.10): #tolérance de 10%
        type_de_choc=type_de_choc_Temp#+'j'
    else:
        type_de_choc=type_de_choc_Temp
    valeurs_param_tester['shock parameters']['shocktype'] = type_de_choc
    
    #log10-densité pré-choc (nHinit en cm-3)
    log10_nH = np.double(inputDir.split('-')[1][1:])
    ##on convertit le format entier des de log10_nH en entier ou en décimaux
    reste_log10_nH = log10_nH - np.round(log10_nH)
    if(simpilified_log10nH_and_Bbeta):
        if (np.abs(reste_log10_nH)/np.round(log10_nH)<1e-3): #tolérance de 0.1%
            log10_nH = np.int(np.round(log10_nH))    
    valeurs_param_tester['shock parameters']['nH(cm-3)'] = 10**log10_nH
    
    #paramètre b de champs magnétique:
    Bbeta = np.double(inputDir.split('-')[2][1:]) 
    reste_Bbeta = Bbeta - np.round(Bbeta)
    if(simpilified_log10nH_and_Bbeta):
        if (np.abs(reste_Bbeta)/np.round(Bbeta)<1e-3): #tolérance de 0.1%
            Bbeta = np.int(np.round(Bbeta)) 
    
        
    Bfield = Bbeta * np.sqrt(10**log10_nH)                                      #Bfield = Bbeta * sqrt(nH) (micro Gauss)
    valeurs_param_tester['shock parameters']['B(microGauss)'] = Bfield
    #print(">>>>Bbeta = ", Bbeta , ">>>>>Bbeta infomhd = ", parameters_Dict_file['shock parameters'] ['B(microGauss)']/np.sqrt(parameters_Dict_file['shock parameters'] ['nH(cm-3)']))
    #print(">>>>Bfield = ", valeurs_param_tester['shock parameters']['B(microGauss)']  , ">>>>>Bbeta infomhd = ", parameters_Dict_file['shock parameters'] ['B(microGauss)'] )
    
    #vitesse du choc:
    vitesse_du_choc = vit
    valeurs_param_tester['shock parameters']['Vs(km.s-1)'] = vitesse_du_choc

    #âge du choc:
    age_du_choc = ag
    valeurs_param_tester[ 'integration parameters']['timeJ'] = age_du_choc     
    
    """Test n°2 à proprement parler ====================="""
    serie_de_test2 = {}
    for category in list(param_tester.keys()):
        serie_de_test2[category] = {}
        for key in param_tester[category]:
            if(key=='B(microGauss)'):
                serie_de_test2[category][key] = (np.abs(valeurs_param_tester[category][key]-parameters_Dict_file[category][key])/parameters_Dict_file[category][key]<0.03)
                #étrangement, même si on entre Bbeta = 1 dans input_mhd.in dans le code de Paris-Durham, la valeur retournée de Bfield dans info_mhd.out peut 
                #être arrondie à l'entier prêt, ce que je ne comprends pas : peut-être le code de Paris-Durham effectue une approximation à un niveau ou
                #bien seulement au moment de sortir les valeurs de info_mhd.out
                #c'est par exemple le cas pour les modèles cj-n5-b1
            else:
                serie_de_test2[category][key] = (valeurs_param_tester[category][key]==parameters_Dict_file[category][key])
            
    
    verite = True #si tout va bien, verite reste True
    for category in list(param_tester.keys()):
        for key in param_tester[category]:
            verite = verite and serie_de_test2[category][key]
            
    if(print_tests) : print(""">>>>>>>>>>Test 2: Tester que les paramètres {cx-ny-bz; v_s; t} du code de choc retrouvés dans le finchier input 'info_mhd.out'    correspond bien aux paramètres input dans la fonction de test (~test 1) ==> 1 booléen propre + message d'erreur\n""")
    try:
        assert verite #teste si tous les paramètres du fichier 'info_mhd.out' renseignés et ceux input depuis la fonction mère sont correspondants
    
    except AssertionError:
        errormsg_part = ''
        for category in list(param_tester.keys()):
            errormsg_part = errormsg_part + '\n)' + category +' : \n\n'
            for key in param_tester[category]:                
                errormsg_part = errormsg_part + '>> ' + key + ' :\n' 
                errormsg_part = errormsg_part + np.str(valeurs_param_tester[category][key]) + ' (inputs fonction mere)'+'\n'
                errormsg_part = errormsg_part + np.str(parameters_Dict_file[category][key]) + ' (inputs dans \'info_mhd.out\')' + '\n'
        if(print_tests): raise AssertionError("4]Attention, dans la fonction {0:}(), l'utilisateur entre un modele \n nommé ../../{1:}, dans lequel les inputs comparés sont \n {2:} ".format(nom_de_fonction,inputDir+'/'+va_tb_tester,errormsg_part))
        else: raise AssertionError('Test 2 échoué\n')
    else:
        if(print_tests) : print("5]Dans la fonction {3:}(), dans le modele \n nommé {0:}, aucun problème de correspondance entre les paramètres input testés dans la fonction mère \n et le fichier 'info_mhd.out' ; le chemin complet est : {2:}\n".format(inputDir+'/'+va_tb_tester,cx_ny_bz+'/'+va_tb,chemin,nom_de_fonction))
        else: print('Test 2 réussi\n')




    
    """Test 3: Meme chose que le test 2 mais où les paramètres tirés de 'info_mhd.out' sont comparées au dossier actuel d'output où travaille la fonction mère
    (donc seuls cx-ny-bz et v'a'-t'b' peuvent être testés)"""
        
    #conversion de l'input depuis la fonction mère pour le test:
    chemin3      = os.getcwd()
    chemin_list3 = chemin3.split('/')
    
    #1ère série de variables testées retournées par os.getcwd()
    va_tb3       = chemin_list3[-1] # nom du premier sous-dossier du niveau inférieur au répertoire actuel : si tout est bien fait, doit être de format v'a'-t'b'
    cx_ny_bz3    = chemin_list3[-2] # nom du deuxième sous-dossier du niveau inférieur au répertoire actuel : si tout est bien fait, doit être de format cx-ny-bz
        
    #type de choc
    verite3 = True #si tout va bien, ce bool reste True
    
    type_de_choc3 = ''
    type_de_choc_Temp3 = parameters_Dict_file['shock parameters']['shocktype']
    if(type_de_choc_Temp3=='C'): type_de_choc_Temp3 = type_de_choc_Temp3.replace('C','c')
    if(type_de_choc_Temp3=='J'): type_de_choc_Temp3 = type_de_choc_Temp3.replace('J','j')
    if(np.abs(parameters_Dict_file['integration parameters']['timeJ']-9.99E+99)/9.99E+99>0.10): #tolérance de 10%
        type_de_choc3=type_de_choc_Temp3+'j'
    else:
        type_de_choc3=type_de_choc_Temp3
    
    
    #log10-densité pré-choc (nHinit en cm-3)
    log10_nH3 = np.log10(parameters_Dict_file['shock parameters']['nH(cm-3)'])
    ##on convertit le format entier des de log10_nH en entier ou en décimaux
    reste_log10_nH3 = log10_nH3 - np.round(log10_nH3)
    if(simpilified_log10nH_and_Bbeta):
        if (np.abs(reste_log10_nH3)/np.round(log10_nH3)<1e-2): #tolérance de 1%
            log10_nH3 = np.int(np.round(log10_nH3))    
    
    
    #paramètre b de champs magnétique:
    Bfield3 = parameters_Dict_file['shock parameters']['B(microGauss)']
    Bbeta3 = Bfield3 / np.sqrt(parameters_Dict_file['shock parameters']['nH(cm-3)']) #Bbeta_(no units) = Bfield_(micro Gauss) /  sqrt(nH_{cm-3}) 
    reste_Bbeta3 = Bbeta3 - np.round(Bbeta3)
    if(simpilified_log10nH_and_Bbeta):
        if (np.abs(reste_Bbeta3)/np.round(Bbeta3)<1e-2): #tolérance de 1%
            Bbeta3 = np.int(np.round(Bbeta3)) 
    
    cx_ny_bz_file3 = type_de_choc3 + '-n' + np.str(log10_nH3) + '-b' + np.str(Bbeta3)
    verite3 = verite3 and (cx_ny_bz_file3==cx_ny_bz3)    

    #vitesse du choc:
    vitesse_du_choc3 = parameters_Dict_file['shock parameters']['Vs(km.s-1)']
    reste_vitesse_du_choc3 = vitesse_du_choc3 - np.round(vitesse_du_choc3)
    if(simpilified_vs):
        if (np.abs(reste_vitesse_du_choc3)/np.round(vitesse_du_choc3)<1e-2): #tolérance de 1%
            vitesse_du_choc3 = np.int(np.round(vitesse_du_choc3)) 

    #âge du choc:
    age_du_choc3 = parameters_Dict_file['integration parameters']['timeJ']
    reste_age_du_choc3 = age_du_choc3 - np.round(age_du_choc3)
    #print(np.str(age_du_choc3), va_tb3)
    if(simpilified_timeJ):
        if (np.abs(reste_age_du_choc3)/np.round(age_du_choc3)<1e-2): #tolérance de 1%
            age_du_choc3 = np.int(np.round(age_du_choc3)) 
    va_tb_file3  = 'v' + np.str(vitesse_du_choc3) + '-t' + np.str(age_du_choc3)
    
    verite3 = verite3 and (va_tb_file3==va_tb3)    
    
    #message d'erreur a afficher si erreur detectée
    errormsg_part3 = ''
    errormsg_part3 = errormsg_part3 + 'Chemin dossier actuel : ../../' + cx_ny_bz3 +'/' + va_tb3 + '\n\n'
    
    errormsg_part3 = errormsg_part3 + '\n)shock parameters : \n\n'
    
    errormsg_part3 = errormsg_part3 + '>> shock type : ' 
    errormsg_part3 = errormsg_part3 + np.str(type_de_choc3) + ' (inputs dans \'info_mhd.out\')' + '\n'
    errormsg_part3 = errormsg_part3 + '>> log10_nH(cm-3) : ' 
    errormsg_part3 = errormsg_part3 + np.str(log10_nH3) + ' (inputs dans \'info_mhd.out\')' + '\n'
    errormsg_part3 = errormsg_part3 + '>> Bbeta : ' 
    errormsg_part3 = errormsg_part3 + np.str(Bbeta3) + ' (inputs dans \'info_mhd.out\')' + '\n'
    errormsg_part3 = errormsg_part3 + '>> Vs (km.s-1) : ' 
    errormsg_part3 = errormsg_part3 + np.str(vitesse_du_choc3) + ' (inputs dans \'info_mhd.out\')' + '\n'    
    
    errormsg_part3 = errormsg_part3 + '\n)integration parameters : \n\n'
    
    errormsg_part3 = errormsg_part3 + '>> timeJ : ' 
    errormsg_part3 = errormsg_part3 + np.str(age_du_choc3) + ' (inputs dans \'info_mhd.out\')' + '\n'
    
    
    if(print_tests) :print(""">>>>>>>>>>Test 3: Meme chose que le test 2 mais où les paramètres tirés de 'info_mhd.out' sont comparées au dossier actuel d'output où travaille la fonction mère (donc seuls cx-ny-bz et v'a'-t'b' peuvent être testés)\n""")
    try:
        assert verite3 #teste si tous les paramètres du fichier 'info_mhd.out' renseignés et ceux dans le chemin du dossier actuel de travail sont correspondants
    
    except AssertionError:
        if(print_tests): raise AssertionError("6]Attention, dans la fonction {0:}(), lon travaille sur le modele \n nommé ../../{1:}, dans lequel l'input file comparé au dossier actuel donne \n {2:} ".format(nom_de_fonction,inputDir+'/'+va_tb_tester,errormsg_part3))
        else: raise AssertionError('Test 3 échoué\n')
    else:
        if(print_tests): print("7]Dans la fonction {2:}(), dans le modele \n nommé {0:}, aucun problème de \n correspondance entre les paramètres déduits du chemin de du dossier de travail actuel  \n et issus du fichier 'info_mhd.out' ; le chemin complet est : {1:}\n".format(cx_ny_bz3+'/'+va_tb3,chemin3,nom_de_fonction))      
        else: print('Test 3 réussi\n')
        
    if(return_info_mhd==True):
        return parameters_Dict_file
    
"""END FUNCTION verifieModele()"""