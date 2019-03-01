#### Installation

L'installation calssique de Mongoose est toujours diponible.

make         # Builds Mongoose (uses CMake) and runs the demo


La commande suivante ne génère plus l'executable classique du projet Mongoose.

sudo make install 

L'executable sera toujours placé produit ici : "build/bin/mongoose"
Celui çi est un appel à edge_cut basique directement en C++. (Voir l'exemple présent dans "Doc/Mongoose_UserGuide.pdf")


#### Usage

Pour lancer un test depuis l'executable C++:

./build/bin/mongoose filename.mtx

Les matrices de tests du projet sont dans le dossier Matrix



## Version C
## Installation

Pour l'instant aucun Makefile n'est présent pour automatiser la compilation.

L'ensemble des .o nécéssaire sont dans le dossier "objects". Recompiler avec g++ les fichiers .cpp de Mongoose si il une modification a été apporté.
(Il est nécéssaire d'ajouter l'option -c à la compilation pour ne pas faire intervenir le linkage à cette étape)
g++ -c Source/classe.cpp -o2 objects/classe.o


Le fichier "mongooseApplication/connectorTest.c" se compile lui à l'aide du compilateur gcc. Le .o résultant est aussi stocké dans le dossier "objects".

Une fois l'ensemble des .o mis à jours. Lancer la commande suivante pour créer l'executable :

g++ objects/*.o -o connectorTest


## Usage

L'executable prend la matrice sur son entrée stdin.

cat Matrix/file.mtx | ./connectorTest

Celui ci dispose aussi d'une option pouvant être mise à A ou B.(par défaut c'est A qui est choisie)

A : appelle direct à edge_cut depuis la matrice en entrée. (passage par la structure spasm et ensuite GraphC)
B : execute la recherche de partition modulaire sur la matrice d'entrée. Le résultat de cette partition est envoyé à edge_cut. (passage par la structure spasm et ensuite GraphC)


#### TODO

- Changer l'option de compilation des fichiers du dossier "objects" et de l'executable final.

- Ajouter des timers dans mongoose.cpp, conectorTest.c et Mongoose_EdgeCut_Connector.cpp

- Modifier modules.c afin d'utiliser les permutations factorisantes dans l'algo de modular_partition. (sur ce sujet voir aussi : https://github.com/antonovvk/decmod/blob/master/dm.c implémentation C de la décomposition modulaire par permutation factorisante)

- Créer un header à modules.c afin de ne pas importer directement le .c dans connectorTest.c


