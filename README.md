# writeMesh (write finite element mesh to inp, bdf, and msh files)



**writeMesh** is a MATLAB project for writing MATLAB finite element mesh to inp, bdf, and msh files. writeMesh is a part of [Im2mesh package](https://github.com/mjx888/im2mesh). Tested in MATLAB R2017b.

writeMesh supports mesh with single or multiple phases. Note that phase is also known as part, domain, or physical surface. In finite element modeling of composite materials, each phase in the mesh represents a distinct material component.

writeMesh currently only woks for 2D mesh. I will extend it to 3D mesh later this year.



## Mesh file formats

- `inp` file (Abaqus)
- `bdf` file (Nastran bulk data, compatible with COMSOL) 
- `msh` file (Gmsh mesh file format)
- For other formats, you can import the generated `msh` file into Gmsh and then export.



## How to start

4 examples are provided. Examples are live script `mlx` files (`demo01.mlx` ~ `demo04.mlx`). If you find some text in the `mlx` file is missing, read the `html` file instead. Examples are also available as `html` files in the folder "demo_html".

**Examples:**

- [demo01](https://mjx888.github.io/writeMesh/demo01.html) - Triangular mesh with single phase
- [demo02](https://mjx888.github.io/writeMesh/demo02.html) - Triangular mesh with multiple phases
- [demo03](https://mjx888.github.io/writeMesh/demo03.html) - Quadrilateral mesh
- [demo04](https://mjx888.github.io/writeMesh/demo04.html) - Mesh with quadratic elements



## Function and input arguments

**Functions:**

```matlab
printInp2d( vert, ele, tnum, ele_type, precision, file_name );
printBdf2d( vert, ele, tnum, ele_type, precision, file_name );
printMsh( vert, ele, tnum, conn, precision, file_name );
```

In the above functions, the 1st and 2nd input arguments are required. The other input arguments are optional. The typical usage of function `printInp2d` is shown below. 

```matlab
printInp2d( vert, ele );
printInp2d( vert, ele, [], [], [], file_name );
printInp2d( vert, ele, [], [], precision, file_name );
printInp2d( vert, ele, tnum, [], precision, file_name );
printInp2d( vert, ele, tnum, ele_type, precision, file_name );
```

The usage of other functions is similar.



**Input arguments:**

- `vert`: Mesh nodes. It’s a Nn-by-2 matrix, where Nn is the number of nodes in the mesh. Each row of `vert` contains the x, y coordinates for that mesh node.
- `ele`: Mesh elements. For linear triangular mesh, it s a Ne-by-3 matrix. For linear quadrilateral mesh, it s a Ne-by-4 matrix. For quadratic element, the number of column in `ele` doubles. Ne is the number of elements in the mesh. Each row in `ele` contains the indices of the nodes for that mesh element. 
- `tnum`: Label of phase, which corresponds to physical surface tag in Gmsh. `tnum` is a Ne-by-1 array, where Ne is the number of elements. If your mesh has multiple part or multiple phase, you can use `tnum` for labelling. `tnum(j,1) = k;` means the j-th element belongs to the k-th phase. If you mesh has single phase, you can set `tnum` as an empty array.
- `conn`: C-by-2 array of constraining edges, where each row defines an edge. You can set `conn` as an empty array.
- `ele_type`: element type. When omitted, `ele_typ` will be determined based on the size of variable `ele`.
- `precision`: number of digits to the right of the decimal point when writing node coordinates. When omitted, `precision=8;`
- `file_name`: file name of exported file, such as `'aaa.inp'`, `'D:\aaa.inp'`. When omitted, `file_name='test.suffix';`



**Notes:**

- Function `printInp2d` works for triangular and quadrilateral mesh; works for linear and quadratic element.
- Function `printBdf2d` does not work for quadratic element.
- Function `printMsh` does not work for quadrilateral mesh. 
- Function `printMsh` can be modified to support quadratic triangular element. I'm just not sure whether this is necessary. If you prefer `msh` file with quadratic triangular element, feel free to let me know.
- I have tested the `msh` file generated by function `printMsh` in Gmsh 4.13.1. I didn't see any errors. The `msh` file generated by function `printMsh` is MSH file format version 4.1. Please refer to the Gmsh manual about details.



## Cite as

The development of writeMesh is nontrivial. If you use writeMesh for your project, please cite Im2mesh package as follows.

Ma, J., & Li, Y. (2025). Im2mesh: A MATLAB/Octave package for generating finite element mesh based on 2D multi-phase image (2.1.5). Zenodo. https://doi.org/10.5281/zenodo.14847059



## Acknowledgments

Great thanks Dr. Yang Lu for providing valuable suggestions.
