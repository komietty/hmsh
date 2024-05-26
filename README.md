# hmsh -- Tiny halfedge structured mesh library
Hmsh is header only tiny library for generating halfedge structured mesh from basic triangle mesh. 
This structure contains many basic properties needed for geometry processing, such as gradient, hodgestar, homology generator, and so on.

The only dependancy is Eigen, which is imported automatically by cmake fetchcontent.

Documents for the basic properties will be added soon. 
For the time being, please check the code in `include/hmsh/hmsh.h` and `include/hmsh/hgen.h` 
If you are familier with basic idea of halfedge structure, you will find the concept of this library.