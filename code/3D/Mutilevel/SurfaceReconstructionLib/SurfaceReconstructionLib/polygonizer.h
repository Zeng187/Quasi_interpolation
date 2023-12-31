/**********************************************************************

polygonizer.h

This is Jules Bloomenthal's implicit surface polygonizer from GRAPHICS 
GEMS IV. Bloomenthal's polygonizer is still used and the present code
is simply the original code morphed into C++.

J. Andreas B�rentzen 2003.

**********************************************************************/





#ifndef POLYGONIZER_H
#define POLYGONIZER_H

#include <vector>

enum ToTetraHedralize
	{
		TET = 0,  // use tetrahedral decomposition 
		NOTET = 1  // no tetrahedral decomposition  */
	};

/** The implicit function class represents the implicit function we wish 
		to polygonize. Derive a class from this one and implement your 
		implicit primitive in the eval function. Eval takes x,y,z coordinates and
		returns a value. We assume that the surface is the zero level set 
		and that the negative values are outside. This an arbitrary choice
		which does not make the code less general. */
class ImplicitFunction
{
 public:
  virtual double eval(double,double,double) = 0;
};

struct GPOINT { double x, y, z;	};

typedef GPOINT VERTEX;
typedef GPOINT NORMAL;

/** TRIANGLE struct contains the indices of the vertices comprising 
		the triangle */
struct TRIANGLE
{
  int v0,v1,v2;
};

/** Polygonizer is the class used to perform polygonization.*/
class Polygonizer
{
  std::vector<NORMAL> gnormals;  
  
  ImplicitFunction* func;
  double size;
  int bounds;

 public:
	 std::vector<VERTEX> gvertices;
	 std::vector<TRIANGLE> gtriangles;
	 double xlength, ylength, zlength;
	
	/** Constructor of Polygonizer. The first argument is the ImplicitFunction
			that we wish to polygonize. The second argument is the size of the 
			polygonizing cell. The final arg. is the limit to how far away we will
			look for components of the implicit surface. */
  Polygonizer(ImplicitFunction* _func, double _size, int _bounds,double _xlength,double _ylength,double _zlength):
  func(_func), size(_size), bounds(_bounds),xlength(_xlength),ylength(_ylength),zlength(_zlength) {}

	/** March erases the triangles gathered so far and builds a new 
			polygonization. The first argument indicates whether the primitive
			cell is a tetrahedron (true) or a cube (fale). The final x,y,z 
			arguments indicate a point near the surface. */
  void march(bool tetra, double x, double y, double z);

	/** Return number of triangles generated after the polygonization.
			Call this function only when march has been called. */
  int no_triangles() const
  {
    return gtriangles.size();
  }

	/** Return number of vertices generated after the polygonization.
			Call this function only when march has been called. */
  int no_vertices() const
  {
    return gvertices.size();
  }
	
	/** Return number of normals generated after the polygonization.
			Of course the result of calling this function is the same as
			no_vertices.
			Call this function only when march has been called. */
  int no_normals() const
  {
    return gnormals.size();
  }
	
	/// Return triangle with index i. 
	TRIANGLE& get_triangle(int i) 
  {
    return gtriangles[i];
  }

	/// Return vertex with index i. 
	VERTEX& get_vertex(int i) 
  {
    return gvertices[i];
  }
	
	
	/// Return normal with index i. 
	NORMAL& get_normal(int i) 
  {
    return gnormals[i];
  }


};



#endif
