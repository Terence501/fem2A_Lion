# Éléments Finis en C++

#### Objectif

Écrire [un code aux éléments finis en C++](course/code.md) pour résoudre des 
[problèmes de Poisson](course/poisson.md) 2D avec des éléments triangulaires
linéaires. 

### Explications organisation fichiers

main.cpp : Programme principal permettant de lancer Tests et Simulations
test.h : Fichier où sont définis les tests
simu.h : Fichier où sont définis les simulations
fem.h et fem.cpp : Fichiers où sont définis et déclarés les différents fonctions utiles dans la résolution par éléments finis

Les autres fichiers ont été utilisés tels quels et n'ont pas été modifiés lors du projet.

### Execution sous linux

1) Dans le fichier main.cpp entre les lignes 59 et 63 choisir en les décommentant les simulations que l'on souhaite effectuer. Il est également possible de choisir le maillage grossier (square.mesh) ou fin (square_fine.mesh).
Ou choisir les tests en changeant le False en True dans les lignes 31 à 38.

2) Dans la console, se placer dans le repertoire fem2A_Lion et lancer la commande make;
3) Dans la console, lancer la commande build/fem2a -h pour afficher l'aide sur les options. Cela permettera de choisir de lancer les simulations choisies ou les tests choisis.

Les différents tests sont expliqués dans le rapport joins.

Pour les simulations, si on veut afficher les résultats, il suffit de lancer la commande medit sur le fichier .mesh et appuyer sur la touche m pour afficher les résultats contenus dans le fichier .bb devant se trouver dans le même répertoire.


