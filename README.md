# Geodesic Complexity Visualizations
Visualize relevant features of a polytope consisting of 2-d faces
(tested with [Python 3.8.10](https://www.python.org/downloads/release/python-3810/))

* **Interactive Cut Locus Visualization:**
  
  ![](https://github.com/pranavraj575/mitochondria/blob/main/images/tetrahedron/demo_cut_locus.gif)
  
* **Unfolding Visualization:**
  
  ![](https://github.com/pranavraj575/mitochondria/blob/main/images/cube/demo_unfold.gif)


|**TABLE OF CONTENTS**||
| ------------- | ------------- |
| [Installation Instructions](#installation-assuming-python-is-installed) | instructions for how to set up code |
| [Installation Test](#installation-test) | test example (tetrahedron cut locus) | 
| [Cut Locus Visualization](#cut-locus-visualization) | Instructions for running interactive cut locus visualization |
| [Unfolding Visualization](#unfolding-visualization) | Instructions for running interactive unfolding visualization |
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

Run `voronoi_interactive.py` from terminal/command prompt in the mitochondria folder 
  
  (move to correct folder with ```cd <name of folder>```)

**Example**: 
```bash
python3 voronoi_interactive.py --shape tetrahedron --center-pt --legend 
```
To see the possible arguments, run the following:
```bash
python3 voronoi_interactive.py -h
```
**Note**: try replacing `python3` with `python` if you get error "Python was not found"

## Cut Locus Visualization:

Run `voronoi_interactive.py` from terminal/command prompt in the mitochondria folder 
  
  (move to correct folder with ```cd <name of folder>```)

**Example**: 
```bash
python3 voronoi_interactive.py --shape tetrahedron --center_pt --legend 
```

## Unfolding Visualization:

Run `unfolding.py` from terminal/command prompt in the mitochondria folder 
  
  (move to correct folder with ```cd <name of folder>```)

**Example**: 
```bash
python3 unfolding.py --shape cube
```

## Implemented shapes:
  * **Tetrahedron**:
    ```bash
    python3 voronoi_interactive.py --shape tetrahedron
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/tetrahedron/p_(-0.4%2C%20-0.6)_face_3.png)
    
  * **Cube**:
    ```bash
    python3 voronoi_interactive.py --shape cube
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/cube/p_(0.8%2C%200.2)_face_1.png)
    
  * **Octahedron**:
    ```bash
    python3 voronoi_interactive.py --shape octahedron
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/octahedron/p_(0.0649519052838329%2C%200.16250000000000003)_face_1.png)
    
  * **Dodecahedron**:
    ```bash
    python3 voronoi_interactive.py --shape dodecahedron --click
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/dodecahedron/p_(0.1902113032590307%2C%200.1368033988749895)_face_0.png)
    
  * **Icosahedron**:
    ```bash
    python3 voronoi_interactive.py --shape icosahedron --click
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/icosahedron/p_(0.10825317547305484%2C%200.1375)_face_1.png)
    
  * **Prism**:
    
    replace `<n>` with the n-gon you want (`3<=n`)
    ```bash
    python3 voronoi_interactive.py --shape prism --n <n>
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/display_images/pentagonal_prism.png)
    
  * **Antiprism**:
    
    replace `<n>` with the n-gon you want (`2<=n`)
    ```bash
    python3 voronoi_interactive.py --shape antiprism --n <n>
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/antiprisms/5/p_(-0.21650635094610957%2C%20-1.034680636892047)_face_0.png)
    
  * **Pyramid**:
    
    replace `<n>` with the n-gon you want (`3<=n<=5`)
    ```bash
    python3 voronoi_interactive.py --shape pyramid --n <n>
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/display_images/pentagonal_pyramid.png)
    
  * **Elongated Pyramid**:
    
    replace `<n>` with the n-gon you want (`2<=n<=5`)
    ```bash
    python3 voronoi_interactive.py --shape longpyramid --n <n>
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/display_images/elongated_pentagonal_pyramid.png)
    
  * **Bipyramid**:
    
    replace `<n>` with the n-gon you want (`2<=n<=5`)
    ```bash
    python3 voronoi_interactive.py --shape bipyramid --n <n>
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/display_images/pentagonal_bipyramid.png)
    
  * **Elongated Bipyramid**:
    
    replace `<n>` with the n-gon you want (`3<=n<=5`)
    ```bash
    python3 voronoi_interactive.py --shape longbipyramid --n <n>
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/display_images/elongated_octahedron.png)
    
  * **Torus**:
    ```bash
    python3 voronoi_interactive.py --shape torus
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/2torus/sample.png)
  
  * **Mirror** (two n-gons glued to each other):
    

    replace `<n>` with the n-gon you want (`3<=n`)
    ```bash
    python3 voronoi_interactive.py --shape mirror --n <n>
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/display_images/hexagonal_mirror.png)
