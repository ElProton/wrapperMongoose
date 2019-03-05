#### Installation

L'installation classique de Mongoose est toujours diponible.

make         # Builds Mongoose (uses CMake) and runs the demo
sudo make install    # installs the mongoose library in /usr/local/lib

L'executable sera toujours placé produit ici : "build/bin/mongoose"
Celui çi est un appel à edge_cut basique directement en C++. (Voir l'exemple présent dans "Doc/Mongoose_UserGuide.pdf")


#### Usage

Pour lancer un test depuis l'executable C++:

./build/bin/mongoose filename.mtx

Les matrices de tests du projet sont dans le dossier Matrix


## Version C
## Installation

Il y a un Makefile dans mongooseApplication.

Il utilise la libmongoose, qui est dans build/lib ou bien dans /usr/local/lib (si on a fait sudo make install). Comme elle n'est pas dans un emplacement standard, pour lancer l'exécutable, il faudra dire :

export LD_LIBRARY_PATH = [path où se trouve libmongoose.so]

Tous les fichiers sont désormais centralisés dans mongooseApplication/

## Usage

L'executable prend la matrice sur son entrée stdin.

cat Matrix/file.mtx | ./connectorTest

Celui ci dispose aussi d'une option pouvant être mise à A ou B.(par défaut c'est A qui est choisie)

A : appelle direct à edge_cut depuis la matrice en entrée. (passage par la structure spasm et ensuite GraphC)
B : execute la recherche de partition modulaire sur la matrice d'entrée. Le résultat de cette partition est envoyé à edge_cut. (passage par la structure spasm et ensuite GraphC)


#### TODO

- corriger le passage par l'option B.

- Modifier modules.c afin d'utiliser les permutations factorisantes dans l'algo de modular_partition. (sur ce sujet voir aussi : https://github.com/antonovvk/decmod/blob/master/dm.c implémentation C de la décomposition modulaire par permutation factorisante)
