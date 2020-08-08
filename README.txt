# Options barières dans le LIBOR Market Model

Nous présentons ici notre projet intitulé “Options barières dans le LIBOR Market Model”.


# Rapport

Un rapport détaillé au format pdf est disponible dans le dossier rapport


# Installation du programme

Le programme permettant de calculer le prix des options peut être installé à l’aide des commandes:

cd src
make all clean

L’executable sera automatiquement copié dans le dossier mère (où se trouve le fichier README.txt).
On pourra revenir à ce dossier à l’aide de la commande:

cd ..


En cas de probleme a la lecture du Makefile, les instructions suivantes permettent de compiler le projet (une fois dans le doosier src)

g++  -I ./libs/boost_1_71_0 -I ./libs/Eigen -o Monte_Carlo.o -c Monte_Carlo.cpp -std=gnu++14 -O2
g++  -I ./libs/boost_1_71_0 -I ./libs/Eigen -o main.o -c main.cpp -std=gnu++14 -O2
g++  -o barrier_option_libor Monte_Carlo.o main.o -lm 


On peut ensuite copier l'executable barrier_option_libor avec la commande:
cp barrier_option_libor ../

# Utilisation du programme

Le nom de l’exécutable est barrier_option_libor
Il se lance à l’aide de la commande:

./barrier_option_libor

Les menus à l’intérieur du logiciel sont explicites.
Avant chaque simulation du prix, il est demandé de saisir les paramètres du modèle. 
On peut ensuite réaliser plusieurs simulations avec ces paramètres. 
Entre chaque simulation du prix, il est possible de changer le pas de discrétisation
ainsi que le nombre de simulations Monte-Carlo. 
Cependant, il n’est pas possible de changer les autres paramètres sans revenir au menu initial.


# Listing des fichiers

- rapport/Oulefki_Szymanski.pdf : rapport du projet
- README.txt 

- src/libs : dossier contenant les deux librairies utilisées dans ce projet: Eigen et Boost

- src/Makefile 

Le code source du projet est composé des fichiers suivants:
- src/GeneralCapletPricing.hpp
- src/GeneralSwaptionPricing.hpp
- src/Header.hpp
- src/main.cpp
- src/Monte_Carlo.cpp
- src/Monte_Carlo.hpp
- src/Timer.hpp
- src/VanillaBarrier.hpp
