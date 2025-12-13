[//]: <> (hello)
# Polyhedra Cut Locus Calculation
Visualize cut loci of the 2-d surfaces of polyhedra

Tested with [Python 3.8.10](https://www.python.org/downloads/release/python-3810/) on Ubuntu (versions 20.04, 22.04, and 24.04) and on Windows 11

[//]: <> (python cut_locus.py -s tetrahedron --legend --font 14 --display-dims 7.5 5.5)
[//]: <> (python unfolding.py -s cube --legend --font 11 --display-dims 6.5 5.5)
[//]: <> (python unfolding.py -s octa --voronoi --legend --font 11)
* **Interactive Cut Locus Visualization:**

  ![](https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/tetrahedron/demo_cut_locus.gif)
* **Path Unfolding Visualization:**
  
  ![](https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/cube/demo_unfold.gif)
* **Voronoi Star Unfolding Visualization:**

  ![](https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/octahedron/voronoi_star_demo.gif)

<table>
    <tr>
        <td colspan="2" align="center"><strong>TABLE OF CONTENTS</strong></td>
    </tr>
    <tr>
        <td><a href=#installation-assuming-python-is-installed>Installation Instructions</a></td>
        <td>instructions for how to set up code</td>
    </tr>
    <tr>
        <td><a href=#installation-test>Installation Test</a></td>
        <td>test example (tetrahedron cut locus)</td>
    </tr>
    <tr>
        <td><a href=#cut-locus>Cut Locus</a></td>
        <td>instructions for running interactive cut locus visualization</td>
    </tr>
    <tr>
        <td><a href=#path-unfolding>Path Unfolding</a></td>
        <td>instructions for running interactive path unfolding visualization</td>
    </tr>
    <tr>
        <td><a href=#implemented-shapes>Implemented Shapes</a></td>
        <td>list of implemented polyhedra and examples of their cut loci</td>
    </tr>
    <tr>
        <td><a href=https://github.com/pranavraj575/polyhedra_cut_locus/tree/main/src#to-implement-new-shapes>Implement New Shapes</a></td>
        <td>instructions to implement new shapes</td>
    </tr>
</table>

## Installation (assuming [Python](https://www.python.org/downloads/release/python-3810/) is installed)
* **Option 1**: [Download zip](https://github.com/pranavraj575/polyhedra_cut_locus/archive/refs/heads/main.zip), then extract all
  
  run the following in terminal/command prompt to install Python package:

  (replace `<name of folder>` with path of folder you extracted it to)
  ```bash
  cd <name of folder>
  pip3 install -e .
  ```
* **Option 2**: Clone Repository, then install the python package (assumes [Git is installed](https://github.com/git-guides/install-git))
  ```bash
  git clone https://github.com/pranavraj575/polyhedra_cut_locus
  pip3 install -e polyhedra_cut_locus
  ```

## Installation Test

Run `cut_locus.py` from terminal/command prompt in the polyhedra_cut_locus folder 
  
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

Run `cut_locus.py` from terminal/command prompt in the polyhedra_cut_locus folder 
  
  (move to correct folder with ```cd <name of folder>```)

**Example**: 
```bash
python3 cut_locus.py --shape tetrahedron --center-pt --legend 
```

## Path Unfolding:

Run `unfolding.py` from terminal/command prompt in the polyhedra_cut_locus folder

**Example**: 
```bash
python3 unfolding.py --shape cube
```

### Voronoi Star Unfolding:

Same as path unfolding, with the additional argument `--voronoi-star`

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
    [\\]: <> (python cut_locus.py -s tetrahedron --source 3 --point -.5 -0.86602540378 --save images/tetrahedron/tetra_locus_three.png --display-dims 9 6 --legend --no-show --font 14)
    [\\]: <> (python unfolding.py -s tetrahedron --source 0 --point .5 0.86602540378 --label-unwrapping --save images/tetrahedron/initial_voronoi_star.png --display-dims 6.5 5.5 --font 20 --voronoi-star --ignore-points --no-show)
    <img src="https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/tetrahedron/tetra_locus_three.png" width=420 />
    <img src="https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/tetrahedron/initial_voronoi_star.png" width=420 />

  * **Cube**:
    ```bash
    python3 cut_locus.py --shape cube
    ```
    [\\]: <> (python cut_locus.py --shape cube --point .8 .2 --source 2 --legend --display-dims 7 5.5 --save images/cube/demo_cut_locus.png --no-show)
    <img src="https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/cube/demo_cut_locus.png" width=690 />
    
  * **Octahedron**:
    ```bash
    python3 cut_locus.py --shape octahedron
    ```
    [\\]: <> (python cut_locus.py -s octahedron --source 6  --point 0.45 .5  --display-dims 12 6 --font 14 --save images\octahedron\unlabeled_locus_tau1.png --legend --no-show)
    <img src="https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/octahedron/unlabeled_locus_tau1.png" width=690 />
    
  * **Dodecahedron**:
    ```bash
    python3 cut_locus.py --shape dodecahedron --click
    python3 unfolding.py --shape dodecahedron --click --voronoi-star
    ```
    [\\]: <> (python cut_locus.py -s dodecahedron --source 11 --point 0 0 --display-dims 12 10 --font 17 --save images/dodecahedron/dodeca_locus_zero.png --no-show)
    [\\]: <> (python unfolding.py -s dodecahedron --source 0 --point 0 0 --display-dims 11 10 --font 17 --voronoi-star --save images/dodecahedron/dodeca_star_zero.png --ignore-points --no-show)
    <img src="https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/dodecahedron/dodeca_locus_zero.png" width=420 />
    <img src="https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/dodecahedron/dodeca_star_zero.png" width=420 />
    
  * **Icosahedron**:
    ```bash
    python3 cut_locus.py --shape icosahedron --click
    python3 unfolding.py --shape icosahedron --click --voronoi-star
    ```
    [\\]: <> (python cut_locus.py -s icosahedron --source 17 --point 0 0 --display-dims 15 10 --font 17 --save images/icosahedron/icosa_locus_zero.png --no-show)
    [\\]: <> (python unfolding.py -s icosahedron --source 0 --point 0 0 --display-dims 12 10 --font 17 --voronoi-star  --save images/icosahedron/icosa_star_zero.png --no-show)
    <img src="https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/icosahedron/icosa_locus_zero.png" width=420 />
    <img src="https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/icosahedron/icosa_star_zero.png" width=420 />
    
  * **Truncated Tetrahedron**:
    ```bash
    python3 cut_locus.py --shape icosahedron --click
    python3 unfolding.py --shape icosahedron --click --voronoi-star
    ```
    [\\]: <> (python cut_locus.py --shape trunc-tetrahedron --legend --source 0 --point 0 0 --display-dims 6 8 --save images/display_images/truncated_tetrahedron.png --no-show)
    <img src="https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/icosahedron/icosa_locus_zero.png" width=420 />
    <img src="https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/icosahedron/icosa_star_zero.png" width=420 />
    
    
  * **Prism**:
    
    replace `<n>` with the n-gon you want (`3<=n`)
    ```bash
    python3 cut_locus.py --shape prism --n <n>
    ```
    [\\]: <> (python cut_locus.py --shape prism --n 5 --legend --source 5 --point -.2 .69 --display-dims 9 6 --save images/display_images/pentagonal_prism.png --no-show)
    <img src="https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/display_images/pentagonal_prism.png" width=690 />
    
  * **Antiprism**:
    
    replace `<n>` with the n-gon you want (`2<=n`)
    ```bash
    python3 cut_locus.py --shape antiprism --n <n>
    ```
    [\\]: <> (python cut_locus.py --shape antiprism --n 5 --legend --source 0 --point -1.3 1.5 --display-dims 10 8 --save images/display_images/pentagonal_antiprism.png --no-show)
    <img src="https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/display_images/pentagonal_antiprism.png" width=690 />
    
  * **Pyramid**:
    
    replace `<n>` with the n-gon you want (`3<=n<=5`)
    ```bash
    python3 cut_locus.py --shape pyramid --n <n>
    ```
    [\\]: <> (python cut_locus.py --shape pyramid --n 5 --legend --source 0 --point -.69 .420 --display-dims 12 5 --save images/display_images/pentagonal_pyramid.png --no-show)
    <img src="https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/display_images/pentagonal_pyramid.png" width=690 />
    
  * **Elongated Pyramid**:
    
    replace `<n>` with the n-gon you want (`2<=n<=5`)
    ```bash
    python3 cut_locus.py --shape longpyramid --n <n>
    ```
    [\\]: <> (python cut_locus.py --shape longpyramid --n 5 --legend --source 0 --point .4 -.1 --display-dims 8 5 --save images/display_images/elongated_pentagonal_pyramid.png --no-show)
    <img src="https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/display_images/elongated_pentagonal_pyramid.png" width=690 />
    
  * **Bipyramid**:
    
    replace `<n>` with the n-gon you want (`2<=n<=5`)
    ```bash
    python3 cut_locus.py --shape bipyramid --n <n>
    ```
    [\\]: <> (python cut_locus.py --shape bipyramid --n 5 --legend --source 0 --point .5 -.8 --display-dims 12 5 --save images/display_images/pentagonal_bipyramid.png --no-show)
    <img src="https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/display_images/pentagonal_bipyramid.png" width=690 />
    
  * **Elongated Bipyramid**:
    
    replace `<n>` with the n-gon you want (`3<=n<=5`)
    ```bash
    python3 cut_locus.py --shape longbipyramid --n <n>
    ```
    [\\]: <> (python cut_locus.py --shape longbipyramid --n 4 --legend --source 0 --point .5 -.8 --display-dims 7 5 --save images/display_images/elongated_octahedron.png --no-show)
    <img src="https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/display_images/elongated_octahedron.png" width=690 />
    
  * **Torus**:
    ```bash
    python3 cut_locus.py --shape torus
    ```
    [\\]: <> (python cut_locus.py --shape torus --point .69 .420 --source 00 --display-dims 5.5 5.5 --save images/2torus/sample.png --no-show)
    <img src="https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/2torus/sample.png" width=690 />
  
  * **Mirror** (two n-gons glued to each other):

    replace `<n>` with the n-gon you want (`3<=n`)
    ```bash
    python3 cut_locus.py --shape mirror --n <n>
    ```
    [\\]: <> (python cut_locus.py --shape mirror --legend -n 6 --point .5 .5 --source 0 --display-dims 8 4 --save images/display_images/hexagonal_mirror.png --no-show)
    <img src="https://github.com/pranavraj575/polyhedra_cut_locus/blob/main/images/display_images/hexagonal_mirror.png" width=690 />
