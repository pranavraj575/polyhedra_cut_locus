[//]: <> (hello)
# Geodesic Complexity Visualizations
Visualize relevant features of a polytope consisting of 2-d faces
(tested with [Python 3.8.10](https://www.python.org/downloads/release/python-3810/))

[//]: <> (python cut_locus.py -s tetra --legend --font 14 --width 7.5 --height 5.5)
[//]: <> (python unfolding.py -s cube --legend --font 11 --width 7.5 --height 5.5)
[//]: <> (python unfolding.py -s octa --voronoi --legend --font 11)
* **Interactive Cut Locus Visualization:**

  ![](https://github.com/pranavraj575/mitochondria/blob/main/images/tetrahedron/demo_cut_locus.gif)
* **Unfolding Visualization:**
  
  ![](https://github.com/pranavraj575/mitochondria/blob/main/images/cube/demo_unfold.gif)
* **Voronoi Star Unfolding Visualization:**

  ![](https://github.com/pranavraj575/mitochondria/blob/main/images/octahedron/voronoi_star_demo.gif)


|**TABLE OF CONTENTS**||
| ------------- | ------------- |
| [Installation Instructions](#installation-assuming-python-is-installed) | instructions for how to set up code |
| [Installation Test](#installation-test) | test example (tetrahedron cut locus) | 
| [Cut Locus](#cut-locus) | Instructions for running interactive cut locus visualization |
| [Unfolding](#unfolding) | Instructions for running interactive unfolding visualization |
| [Implemented Shapes](#implemented-shapes) | List of implemented polyhedra and examples of their cut loci |

## Installation (assuming [Python](https://www.python.org/downloads/release/python-3810/) is installed)
* **Option 1**: [Download zip](https://github.com/pranavraj575/mitochondria/archive/refs/heads/main.zip), then extract all
  
  run the following in terminal/command prompt to install Python package:

  (replace `<name of folder>` with path of folder you extracted it to)
  ```bash
  cd <name of folder>
  pip3 install -e .
  ```
* **Option 2**: Clone Repository, then install the python package (assumes [Git is installed](https://github.com/git-guides/install-git))
  ```bash
  git clone https://github.com/pranavraj575/mitochondria
  pip3 install -e mitochondria
  ```

## Installation Test

Run `cut_locus.py` from terminal/command prompt in the mitochondria folder 
  
  (move to correct folder with ```cd <name of folder>```)

**Example**: 
```bash
python3 cut_locus.py --shape tetrahedron --center-pt --legend 
```
To see the possible arguments, run the following:
```bash
python3 cut_locus.py -h
```
**Note**: try replacing `python3` with `python` if you get error "Python was not found"

## Cut Locus:

Run `cut_locus.py` from terminal/command prompt in the mitochondria folder 
  
  (move to correct folder with ```cd <name of folder>```)

**Example**: 
```bash
python3 cut_locus.py --shape tetrahedron --center_pt --legend 
```

## Unfolding:

Run `unfolding.py` from terminal/command prompt in the mitochondria folder

**Example**: 
```bash
python3 unfolding.py --shape cube
```

### Voronoi Star Unfolding:

Same as unfolding, with the additional argument `--voronoi-star`

**Example**: 
```bash
python3 unfolding.py --shape octahedron --voronoi-star
```

## Implemented shapes:
  * **Tetrahedron**:
    ```bash
    python3 cut_locus.py --shape tetrahedron
    python3 unfolding.py --shape tetrahedron --voronoi-star
    ```
    [\\]: <> (python cut_locus.py -s tetra --source 3 --point-x -.5 --point-y -0.86602540378 --save images/tetrahedron/tetra_locus_three.png --height 6 --width 9 --legend --no-show --font 14)
    [\\]: <> (python unfolding.py -s tetra --source 0 --point-x .5 --point-y 0.86602540378 --label-unwrapping --save images/tetrahedron/initial_voronoi_star.png --height 5.5 --width 6.5 --font 20 --voronoi-star --ignore-points --no-show)
    <img src="https://github.com/pranavraj575/mitochondria/blob/main/images/tetrahedron/tetra_locus_three.png" width=420 />
    <img src="https://github.com/pranavraj575/mitochondria/blob/main/images/tetrahedron/initial_voronoi_star.png" width=420 />

  * **Cube**:
    ```bash
    python3 cut_locus.py --shape cube
    ```
    [\\]: <> (python cut_locus.py --shape cube --point-x .8 --point-y .2 --source 2 --legend --width 7 --height 5.5 --save images/cube/demo_cut_locus.png --no-show)
    <img src="https://github.com/pranavraj575/mitochondria/blob/main/images/cube/demo_cut_locus.png" width=690 />
    
  * **Octahedron**:
    ```bash
    python3 cut_locus.py --shape octahedron
    ```
    [\\]: <> (python cut_locus.py -s octa --source 6  --point-x 0.45 --point-y .5  --height 6 --width 12 --font 14 --save images\octahedron\unlabeled_locus_tau1.png --legend --no-show)
    <img src="https://github.com/pranavraj575/mitochondria/blob/main/images/octahedron/unlabeled_locus_tau1.png" width=690 />
    
  * **Dodecahedron**:
    ```bash
    python3 cut_locus.py --shape dodecahedron --click
    python3 unfolding.py --shape dodecahedron --click --voronoi-star
    ```
    [\\]: <> (python cut_locus.py -s dodeca --source 11 --point-x 0 --point-y 0 --height 10 --width 12 --font 17 --save images/dodecahedron/dodeca_locus_zero.png --no-show)
    [\\]: <> (python unfolding.py -s dodeca --source 0 --point-x 0 --point-y 0 --height 10 --width 11 --font 17 --voronoi-star --save images/dodecahedron/dodeca_star_zero.png --ignore-points --no-show)
    <img src="https://github.com/pranavraj575/mitochondria/blob/main/images/dodecahedron/dodeca_locus_zero.png" width=420 />
    <img src="https://github.com/pranavraj575/mitochondria/blob/main/images/dodecahedron/dodeca_star_zero.png" width=420 />
    
  * **Icosahedron**:
    ```bash
    python3 cut_locus.py --shape icosahedron --click
    python3 unfolding.py --shape icosahedron --click --voronoi-star
    ```
    [\\]: <> (python cut_locus.py -s icosa --source 17 --point-x 0 --point-y 0 --height 10 --width 15 --font 17 --save images/icosahedron/icosa_locus_zero.png --no-show)
    [\\]: <> (python unfolding.py -s icosa --source 0 --point-x 0 --point-y 0 --height 10 --width 12 --font 17 --voronoi-star  --save images/icosahedron/icosa_star_zero.png --no-show)
    <img src="https://github.com/pranavraj575/mitochondria/blob/main/images/icosahedron/icosa_locus_zero.png" width=420 />
    <img src="https://github.com/pranavraj575/mitochondria/blob/main/images/icosahedron/icosa_star_zero.png" width=420 />
    
  * **Prism**:
    
    replace `<n>` with the n-gon you want (`3<=n`)
    ```bash
    python3 cut_locus.py --shape prism --n <n>
    ```
    [\\]: <> (python cut_locus.py --shape prism --n 5 --legend --source 5 --point-x -.2 --point-y .69 --height 6 --width 9 --save images/display_images/pentagonal_prism.png --no-show)
    <img src="https://github.com/pranavraj575/mitochondria/blob/main/images/display_images/pentagonal_prism.png" width=690 />
    
  * **Antiprism**:
    
    replace `<n>` with the n-gon you want (`2<=n`)
    ```bash
    python3 cut_locus.py --shape antiprism --n <n>
    ```
    [\\]: <> (python cut_locus.py --shape antiprism --n 5 --legend --source 0 --point-x -1.3 --point-y 1.5 --height 8 --width 10 --save images/display_images/pentagonal_antiprism.png --no-show)
    <img src="https://github.com/pranavraj575/mitochondria/blob/main/images/display_images/pentagonal_antiprism.png" width=690 />
    
  * **Pyramid**:
    
    replace `<n>` with the n-gon you want (`3<=n<=5`)
    ```bash
    python3 cut_locus.py --shape pyramid --n <n>
    ```
    [\\]: <> (python cut_locus.py --shape pyramid --n 5 --legend --source 0 --point-x -.69 --point-y .420 --height 5 --width 12 --save images/display_images/pentagonal_pyramid.png --no-show)
    <img src="https://github.com/pranavraj575/mitochondria/blob/main/images/display_images/pentagonal_pyramid.png" width=690 />
    
  * **Elongated Pyramid**:
    
    replace `<n>` with the n-gon you want (`2<=n<=5`)
    ```bash
    python3 cut_locus.py --shape longpyramid --n <n>
    ```
    [\\]: <> (python cut_locus.py --shape longpyramid --n 5 --legend --source 0 --point-x .4 --point-y -.1 --height 5 --width 8 --save images/display_images/elongated_pentagonal_pyramid.png --no-show)
    <img src="https://github.com/pranavraj575/mitochondria/blob/main/images/display_images/elongated_pentagonal_pyramid.png" width=690 />
    
  * **Bipyramid**:
    
    replace `<n>` with the n-gon you want (`2<=n<=5`)
    ```bash
    python3 cut_locus.py --shape bipyramid --n <n>
    ```
    [\\]: <> (python cut_locus.py --shape bipyramid --n 5 --legend --source 0 --point-x .5 --point-y -.8 --height 5 --width 12 --save images/display_images/pentagonal_bipyramid.png --no-show)
    <img src="https://github.com/pranavraj575/mitochondria/blob/main/images/display_images/pentagonal_bipyramid.png" width=690 />
    
  * **Elongated Bipyramid**:
    
    replace `<n>` with the n-gon you want (`3<=n<=5`)
    ```bash
    python3 cut_locus.py --shape longbipyramid --n <n>
    ```
    [\\]: <> (python cut_locus.py --shape longbipyramid --n 4 --legend --source 0 --point-x .5 --point-y -.8 --height 5 --width 7 --save images/display_images/elongated_octahedron.png --no-show)
    <img src="https://github.com/pranavraj575/mitochondria/blob/main/images/display_images/elongated_octahedron.png" width=690 />
    
  * **Torus**:
    ```bash
    python3 cut_locus.py --shape torus
    ```
    [\\]: <> (python cut_locus.py --shape torus --point-x .69 --point-y .420 --source 00 --height 5.5 --width 5.5 --save images/2torus/sample.png --no-show)
    <img src="https://github.com/pranavraj575/mitochondria/blob/main/images/2torus/sample.png" width=690 />
  
  * **Mirror** (two n-gons glued to each other):

    replace `<n>` with the n-gon you want (`3<=n`)
    ```bash
    python3 cut_locus.py --shape mirror --n <n>
    ```
    [\\]: <> (python cut_locus.py --shape mirror --legend -n 6 --point-x .5 --point-y .5 --source 0 --height 4 --width 8 --save images/display_images/hexagonal_mirror.png --no-show)
    <img src="https://github.com/pranavraj575/mitochondria/blob/main/images/display_images/hexagonal_mirror.png" width=690 />
