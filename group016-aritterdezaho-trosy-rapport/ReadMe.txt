# myFem

## Description

l'executable `myFem`nécessite deux fichiers en entrée : un fichier de maillage (`meshFile`) et un fichier décrivant le problème (`problemFile`).  
Deux paramètres supplémentaires permettent de spécifier le type de solveur et la méthode de renumérotation.

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
./myFem ../data/mesh.txt ../data/problem.txt B X
```

