# ImWG2D
This code implements the Projected Kirchhoff migration for a homogeneous 2D waveguide with 
Dirichlet/Dirichlet boundary conditions. The method was introduced in [1] for the case where the array spans the
whole depth of the waveguide. 

Implemented in the code are 4 different scatterers with zero volume:

- a point scatterer,
- a semicircle 'facing' the array',
- a vertical screen with unit reflectivity,
- a horizontal screen with unit reflectivity.

Any scatterer shape or multiple scatterers should be supported by the code, assuming its response matrix may be written by using the 
Born approximation. Closed surfaces that occupy a positive volume inside the waveguide are not supported.

 Call as :

     
     [IKMT,z1,x1] = IKMt(frqs,J)
     
frqs: the frequency (or frequencies) we use in the computation
       (enter an array for multiple frequencies).

J   : selective imaging using subspace projection on the Jth singular 
       vector.
       Enter 0 or leave blank for no projection,
       Enter a vector if you want to project on more than one singular 
       vector.
 
 Examples of use:
 
       [IKMT,z1,x1] = IKMt(73); Uses a single frequency and no projection.
       [IKMt,z1,x2] = IKMt(31.875:3.75:118.125); Uses multiple frequencies
                      and no projection.
       [IKMt,z1,x1] = IKMt(73,1); Uses a single frequency and projection
                       on the first singular vector.
       [IKMt,z1,x1] = IKMt(73,[1,3,5]); Uses a single frequency and 
                      projection on the first, third and fifth singular 
                      vectors (the result of each projection is summed.


 [1] C. Tsogka, D. A. Mitsoudis, and S. Papadimitropoulos. 
 Selective imaging of extended reflectors in two-dimensional waveguides. 
 SIAM J. Imaging Sci., 6(4):2714{2739, 2013.
