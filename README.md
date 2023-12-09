# Cut Locus Visualization
Visualize the cut locus of points on a polytope consisting of 2-d faces

tested with [Python 3.8.10](https://www.python.org/downloads/release/python-3810/)

![](https://github.com/pranavraj575/mitochondria/blob/main/images/tetrahedron/demo.gif)
## Installation (assuming [Python](https://www.python.org/downloads/release/python-3810/) is installed)
* Option 1: [Download zip](https://github.com/pranavraj575/mitochondria/archive/refs/heads/main.zip), then extract all
  run the following in terminal/command prompt to install Python package: (replace `<name of folder>` with path of folder you extracted it to)
  ```bash
  cd <name of folder>
  pip3 install -e .
  ```
* Option 2: Clone Repository, then install the python package
  ```bash
  git clone https://github.com/pranavraj575/mitochondria
  pip3 install -e mitochondria
  ```
  
## Run Cut Locus visualization

Run `voronoi_interactive.py` in the mitochondria folder

Example:
```bash
python3 voronoi_interactive.py --shape tetrahedron --center_pt --legend 
```
To see the possible arguments, run the following:
```bash
python3 voronoi_interactive.py -h
```

### Implemented shapes:
  * Tetrahedron:
    ```bash
    python3 voronoi_interactive.py --shape tetrahedron
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/tetrahedron/p_(-0.4%2C%20-0.6)_face_3.png)
    
  * Cube:
    ```bash
    python3 voronoi_interactive.py --shape cube
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/cube/p_(0.8%2C%200.2)_face_1.png)
    
  * Octahedron:
    ```bash
    python3 voronoi_interactive.py --shape octahedron
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/octahedron/p_(0.0649519052838329%2C%200.16250000000000003)_face_1.png)
    
  * Dodecahedron:
    ```bash
    python3 voronoi_interactive.py --shape dodecahedron --click
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/dodecahedron/p_(0.1902113032590307%2C%200.1368033988749895)_face_0.png)
    
  * Icosahedron:
    ```bash
    python3 voronoi_interactive.py --shape icosahedron --click
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/icosahedron/p_(0.10825317547305484%2C%200.1375)_face_1.png)
    
  * Prism:
    (replace `<n>` with the n-gon you want
    ```bash
    python3 voronoi_interactive.py --shape prism --n <n>
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/display_images/pentagonal_prism.png)
    
  * Antiprism:
    (replace `<n>` with the n-gon you want
    ```bash
    python3 voronoi_interactive.py --shape antiprism --n <n>
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/antiprisms/5/p_(-0.21650635094610957%2C%20-1.034680636892047)_face_0.png)
    
  * Pyramid:
    (replace `<n>` with the n-gon you want (`3<=n<=6`)
    ```bash
    python3 voronoi_interactive.py --shape pyramid --n <n>
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/display_images/pentagonal_pyramid.png)
    
  * Elongated Pyramid:
    (replace `<n>` with the n-gon you want (`3<=n<=6`)
    ```bash
    python3 voronoi_interactive.py --shape longpyramid --n <n>
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/display_images/elongated_pentagonal_pyramid.png)
    
  * Bipyramid:
    (replace `<n>` with the n-gon you want (`3<=n<=6`)
    ```bash
    python3 voronoi_interactive.py --shape bipyramid --n <n>
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/display_images/pentagonal_bipyramid.png)
    
  * Elongated Bipyramid:
    (replace `<n>` with the n-gon you want (`3<=n<=6`)
    ```bash
    python3 voronoi_interactive.py --shape longbipyramid --n <n>
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/display_images/elongated_octahedron.png)
    
  * Torus:
    ```bash
    python3 voronoi_interactive.py --shape torus
    ```
    ![](https://github.com/pranavraj575/mitochondria/blob/main/images/2torus/sample.png)
