# writeMesh (write finite element mesh to inp, bdf, and msh files)



**writeMesh** is a MATLAB project for writing MATLAB finite element mesh to inp, bdf, and msh files. writeMesh is a sub-project of [Im2mesh package](https://github.com/mjx888/im2mesh). 

writeMesh supports mesh with single or multiple phases. Note that phase is also known as part, domain, or physical surface. In finite element modeling of composite materials, each phase in the mesh represents a distinct material component.

**Mesh file formats:**

- `inp` file (Abaqus)
- `bdf` file (Nastran bulk data, compatible with COMSOL) 
- `msh` file (Gmsh mesh file format)
- For other formats, you can import the generated `msh` file into Gmsh and then export.

**Notes:**

- **writeMesh currently only woks for 2D mesh.** I will extend it to 3D mesh later this year.
- Tested in MATLAB R2017b
- No dependency on MATLAB Toolbox

## How to start

After downloading writeMesh, I suggest you start with writeMesh manual in the folder. 4 examples are provided. 

- Examples are live script `mlx` files (`demo01.mlx` ~ `demo04.mlx`). If you find some text in the `mlx` file is missing, please read the `html` file instead.
- Examples are also available as `html` files in the folder "demo_html".

**Examples:**

- [demo01](https://mjx888.github.io/im2mesh_demo_html/demo01.html) - Triangular mesh with single phase
- [demo02](https://mjx888.github.io/im2mesh_demo_html/demo02.html) - Triangular mesh with multiple phases
- [demo03](https://mjx888.github.io/im2mesh_demo_html/demo03.html) - Quadrilateral mesh
- [demo04](https://mjx888.github.io/im2mesh_demo_html/demo04.html) - Triangular mesh with isolated domain.

## Cite as

The development of writeMesh is nontrivial. If you use writeMesh for your project, please cite Im2mesh package as follows.

Ma, J., & Li, Y. (2025). Im2mesh: A MATLAB/Octave package for generating finite element mesh based on 2D multi-phase image (2.1.5). Zenodo. https://doi.org/10.5281/zenodo.14847059

## Acknowledgments

Great thanks Dr. Yang Lu for providing valuable suggestions.
