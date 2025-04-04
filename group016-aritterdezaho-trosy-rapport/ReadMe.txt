# myFem

## Description
L'exécutable `myFem` nécessite deux fichiers en entrée : un fichier de maillage (`meshFile`) et un fichier décrivant le problème (`problemFile`).  
Deux paramètres supplémentaires permettent de spécifier le type de solveur et la méthode de renumérotation.

## Compilation

Pour compiler le projet, suivez les étapes suivantes :

1. Naviguez dans le répertoire `build` :
   ```bash
   cd build
   ```

2. Exécutez `cmake` pour générer les fichiers de construction :
   ```bash
   cmake ..
   ```

3. Compilez le projet avec `make` :
   ```bash
   make
   ```

L'exécutable `myFem` sera généré dans le répertoire `build`.

## Utilisation

```bash
./myFem [meshFile] [problemFile] [solver] [renumType]
```

### Paramètres

- **`meshFile`** : chemin vers le fichier de maillage (ex : `../data/mesh.txt`)
- **`problemFile`** : chemin vers le fichier décrivant le problème (ex : `../data/problem.txt`)
- **`solver`** : type de solveur à utiliser
- **`renumType`** : méthode de renumérotation à appliquer

### Valeurs possibles

#### `solver`
- `B` : solveur bande
- `F` : solveur plein
- `G` : solveur par gradient conjugué

#### `renumType`
- `X` : renumérotation selon l'axe X
- `Y` : renumérotation selon l'axe Y
- `0` : pas de renumérotation

## Exemple

```bash
./myFem ../data/mesh_precis.txt ../data/problem_inconel600.txt B X 
```

## Structure du projet

├── CMakeLists.txt          # Fichier de configuration CMake pour le projet  
├── ReadMe.txt              # Documentation du projet  
├── data                    # Répertoire contenant les données d'entrée  
│   ├── mesh.txt            # Fichier de maillage avec tous les éléments de même taille  
│   ├── mesh_precis.txt     # Fichier de maillage adapté à la géométrie du problème  
│   ├── problem_inconel600.txt # Fichier décrivant le problème pour l'Inconel 600  
│   ├── problem_molybdene.txt  # Fichier décrivant le problème pour le Molybdène  
│   └── problem_inconel_lourd.txt  # Fichier décrivant le problème pour l'Inconel 600 avec une masse 10X plus élevée  
├── src                     
│   ├── fem.c               # Implémentation des fonctions principales FEM  
│   ├── fem.h               # Déclarations des fonctions et structures FEM  
│   ├── main.c              # Point d'entrée principal du programme  
│   └── homework.c          # Fichier pour les fonctions spécifiques au projet  
└── build                   # Répertoire pour les fichiers de construction (généré par CMake)

