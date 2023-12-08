# mitochondria
Visualize the cut locus of points on a polytope consisting of 2-d faces

tested with [Python 3.8.10](https://www.python.org/downloads/release/python-3810/)

## Installation
Clone Repository, then install the python package
```bash
git clone https://github.com/pranavraj575/mitochondria
pip3 install -e mitochondria
```
## Run Cut Locus visualization

Run `voronoi_interactive.py` in the mitochondria folder

Example:
```bash
python3 voronoi_interactive.py --shape antiprism --n 4 --center_pt --legend 
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
  * Cube:
    ```bash
    python3 voronoi_interactive.py --shape cube
    ```
  * Octahedron:
    ```bash
    python3 voronoi_interactive.py --shape octahedron --click
    ```
  * Dodecahedron:
    ```bash
    python3 voronoi_interactive.py --shape dodecahedron --click
    ```
  * Icosahedron:
    ```bash
    python3 voronoi_interactive.py --shape icosahedron --click
    ```
  * Antiprism:
    (replace `<n>` with the n-gon you want
    ```bash
    python3 voronoi_interactive.py --shape antiprism --n <n> --click
    ```
  * Torus:
    ```bash
    python3 voronoi_interactive.py --shape torus
    ```
