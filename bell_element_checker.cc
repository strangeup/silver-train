// A c++ program to test a single element
#include <fenv.h> 
//Generic routines
#include "generic.h"
// My New code to check
#include "C1_basis/C1_curved_elements.h"

using namespace oomph;
using namespace MathematicalConstants;


//#############################################################################//
// Checker class
//#############################################################################//
namespace bell_element_tests {

/// Pretty print a dense matrix, optional argument for omission of dofs
void print_dense_matrix(const DenseMatrix<double>& matrix, const double tol=1e-15)
 {
 const unsigned nrow = matrix.nrow(); 
 const unsigned ncol = matrix.ncol(); 
  //Output any non zeros
  oomph_info<<"\n";
  for(unsigned i=0;i<nrow;++i)
   {
    for(unsigned j=0;j<ncol;++j)
     {
      oomph_info<<(std::abs(matrix(i,j))>tol? matrix(i,j): 0.0 )<<(j!=ncol-1 ? " ":"\n");
     }
   }
 }

/// Class to test the BellElement basis class
class BellElementTestBasis : 
 protected MyShape::BellElementBasis
{
public:
 /// \short Shorthand for a vector of vectors containining the vertices
 typedef Vector<Vector<double> > VertexList;
 
 /// New type of function pointer
 typedef void (*ExactSolnFctPt)(const Vector<double>& x, Vector<double>& p);
 
 /// Constructor
 BellElementTestBasis()
  {
   // Fill in vertices with local vertices 
   Vertices_local = VertexList (3,Vector<double> (2));
   fill_in_local_dofs_position(Vertices_local);
   Vertices = Vertices_local;
  }
 
 /// Destructor
 ~BellElementTestBasis()
  {/*Do nothing everything initialised on upgrade */}
 
 // Private copy and assign - so this should cause a compilation error
private:
 /// Broken copy constructor
 BellElementTestBasis(BellElementTestBasis& dummy) 
  { 
   BrokenCopy::broken_copy("BellElementTestBasis");
  } 
 /// Broken assignment operator
 void operator=(const BellElementTestBasis&) 
  {
   BrokenCopy::broken_assign("BellElementTestBasis");
  }

/// Private data
VertexList Vertices;
VertexList Vertices_local;

public:
 /// Get vertices
 void get_vertices(VertexList& vertices) const
   { vertices = Vertices; }

 const VertexList get_vertices() const
   { return Vertices; }

 /// Get vertices
 void set_vertices(const VertexList& vertices)
   { Vertices = vertices; }

 /// Fill in area coordinates
 void fill_in_area_coordinates(const Vector<double>& s, Shape& omega) const
   {
    omega(0) = s[0];
    omega(1) = s[1];
    omega(2) = 1.0-s[0]-s[1];
   }

 /// Interpolated position
 void interpolated_position(const Vector<double>& s, Vector<double>& x) const
   {
    // Initialize
    unsigned nnode = 3, dim =2 ;
    Shape phi(nnode);
    x = Vector<double>(2,0.0);
    // Now get the area coodinates
    fill_in_area_coordinates(s,phi);
    // Loop over nodes
    for(unsigned i=0;i<nnode;++i)   
     {
      // Loop over components
      for(unsigned alpha=0;alpha<dim;++alpha)
       {
        // Interpolated the position
        x[alpha] += Vertices[i][alpha] * phi(i);
       }
     }
   }

 // Return local coordinates of vertices
 void fill_in_local_dofs_position(VertexList& local_vertices) const
  {
   // Number of nodes
   const unsigned nnode = 3, dim =2;
   // Override via assignment
   local_vertices = VertexList(nnode,Vector<double>(dim,0.0));
   // Modify two of the values
   local_vertices[0][0] =1.0; local_vertices[1][1] =1.0;
  }
 

/// Pretty print a dense matrix, optional argument for omission of dofs
void print_matrix_identity_diff(const DenseMatrix<double>& matrix, const double tol) const
 {
 const unsigned nrow = matrix.nrow(); 
 const unsigned ncol = matrix.ncol(); 
  //Output any non zeros
  oomph_info<<"\n";
  for(unsigned i=0;i<nrow;++i)
   {
    for(unsigned j=0;j<ncol;++j)
     {
      // Identity difference
      const double diff =matrix(i,j)- (i==j ? 1.0 :0.0);
      // Nonzero difference
      if(std::abs(diff) > tol)
       {
        oomph_info<<"Entry ("<<i<<","<<j<<") differs. Difference:"<< diff<<std::endl;
       }
     }
   }
 }

private:
// Return the global dofs of a function on an element
void check_basis_global_dofs(const VertexList& vertices) const
 {
  // Three nodes, 6 hermite dofs
  const  unsigned nnode = 3, ntype =6, dim =2;

  // Initialise the matrix to store the dofs of the basis
  DenseMatrix<double> gdofs(nnode*ntype,nnode*ntype,0.0);

  // Initialise the shape
  Shape psi(nnode,ntype);
  DShape dpsi(nnode,ntype,dim);
  DShape d2psi(nnode,ntype,dim*dim-1);

  // Get the dofs
  for (unsigned i=0; i<nnode; ++i)
   {
    // Fill in basis
    this->d2_basis_eulerian(Vertices_local[i],vertices,psi,dpsi,d2psi);
    // Loop over basis
    for(unsigned ii=0; ii <nnode; ++ii)
     {
      // Loop over basis type  
      for(unsigned jj=0; jj <ntype; ++jj)
       {
        // Flatpack
        const unsigned ishape = ntype*ii + jj;
        // Zeroth derivative of basis
        gdofs(ishape,ntype*i + 0) = psi(ii,jj);
        // First derivatives of basis
        for(unsigned alpha = 0;alpha<dim;++alpha)
         {gdofs(ishape,ntype*i + 1 + alpha) = dpsi(ii,jj,alpha);}
        // Second derivatives of basis
        for(unsigned l = 0;l<dim*dim-1;++l)
         {gdofs(ishape,ntype*i + 3 + l) = d2psi(ii,jj,l);}
       } 
     } 
   }

  // Now output
  print_dense_matrix(gdofs,0); 
  print_matrix_identity_diff(gdofs,0); 
 }

/// Get degrees of freedom
void get_global_dofs(const ExactSolnFctPt& exact_fpt, Vector<double>& gdofs)const 
  {
  // Three nodes, 6 hermite dofs
  const  unsigned nnode = 3, ntype =6;
  // Fill in
  // Loop over basis
  for(unsigned i=0; i <nnode; ++i)
   {
    // Fill in exact_f
    Vector<double> exact_f (6,0.0);
    (*exact_fpt)(Vertices[i],exact_f);
    // Zeroth derivative of basis
    for(unsigned l=0 ; l <ntype ; ++l)
     {gdofs[ntype*i + l] = exact_f[l];}
   } 
  }

public:
// Return the global dofs of a function on an element
double check_analytic_function(const ExactSolnFctPt& exact_fpt) const
 {
  // Three nodes, 6 hermite dofs
  const  unsigned nnode = 3, ntype =6, dim =2;
  // Create integral
  oomph::TGauss<2,5> integral;
  // Now get the traces at several values of s
  const unsigned n_ipoints=integral.nweight();
  // Get the global dofs
  Vector<double> global_dofs (nnode*ntype,0.0);
  get_global_dofs(exact_fpt,global_dofs);
  
  double integrated_squared_error(0.0),integrated_area(0.0);
  // Loop integral points
  for(unsigned ipt=0 ; ipt< n_ipoints; ++ipt)
   {
    // The corresponding local coordinates
    Vector<double> s(2,0.0);
    s[0] = integral.knot(ipt,0);
    s[1] = integral.knot(ipt,1);
    // Initialise approximate function
    double approx_f = 0.0; 
    // The integral weights
    const double weight=integral.weight(ipt);
    // Fill in basis
    // Initialise the shape
    Shape psi(nnode,ntype);
    DShape dpsi(nnode,ntype,dim);
    DShape d2psi(nnode,ntype,dim*dim-1);
    const double J =this->d2_basis_eulerian(s,Vertices,psi,dpsi,d2psi);

    // Now Sum over shape functions
    for(unsigned inod = 0; inod<nnode; ++inod)
     {
      for(unsigned itype =0 ; itype<ntype;++itype)
       {
        // Interpolated f
        approx_f+=global_dofs[ntype*inod+itype]*psi(inod,itype);
       }
     }

    // Get exact
    Vector<double> exact_f(6,0.0),x(2,0.0);
    interpolated_position(s,x);
    (*exact_fpt)(x,exact_f);

    // Calculate L2 norm
    integrated_squared_error+=(exact_f[0]-approx_f)*(exact_f[0]-approx_f)*weight*J;
    integrated_area+=weight*J;
   }
  // Return the result
  return sqrt(integrated_squared_error/integrated_area);
 }

/// Check the basis degrees of freedom on reference element with delta property
void check_local_basis_dofs() const
 {
  // Use the function on these
  check_basis_global_dofs(Vertices_local);
  }

/// Check the basis degrees of freedom on physical element with delta property
void check_basis_global_dofs() const
 {
  // Use the function on these
  check_basis_global_dofs(Vertices);
  }
  
};

}

namespace TestSoln {
void get_p4(const Vector<double>& x, Vector<double>& p)
{
 // Copy of the position Vector
 DenseMatrix<double> a(5,5,0.0);
 a(0,0)= 0.2; 
 a(0,1)= 1.0;  a(1,0)=-0.4; 
 a(1,1)= 0.2;  a(2,0)=-0.9; a(0,2)=-0.1; 
 a(3,0)= 0.4;  a(0,3)= 0.6; a(1,2)= 0.3; a(2,1)=-0.025;
 a(4,0)=-0.7;  a(0,4)= 0.9; a(3,1)= 0.6; a(1,3)=-0.015; a(2,2)=0.3;
 // Zero the vector
 // Now loop over coefficients
 for(unsigned i=0; i<5;++i)
  {
   for(unsigned j=0; j<5-i;++j)
    {
     // Function
     p[0]+=a(i,j)*pow(x[0],i)*pow(x[1],j);
     // Values
     if(i>0)
       p[1]+=a(i,j)*i*pow(x[0],i-1)*pow(x[1],j);
     if(j>0)
       p[2]+=a(i,j)*j*pow(x[0],i)*pow(x[1],j-1);
     // Second derivatives
     if(i>1)
       p[3]+=a(i,j)*i*(i-1)*pow(x[0],i-2)*pow(x[1],j);
     if(i>0 && j>0)
       p[4]+=a(i,j)*j*i*pow(x[0],i-1)*pow(x[1],j-1);
     if(j>1)
       p[5]+=a(i,j)*j*(j-1)*pow(x[0],i)*pow(x[1],j-2);
    }
  }
}

void get_p5(const Vector<double>& x, Vector<double>& p)
{
 // Copy of the position Vector
 DenseMatrix<double> a(6,6,0.0);
 a(0,0)= 0.2; 
 a(0,1)= 1.0;  a(1,0)=-0.4; 
 a(1,1)= 0.2;  a(2,0)=-0.9; a(0,2)=-0.1; 
 a(3,0)= 0.4;  a(0,3)= 0.6; a(1,2)= 0.3; a(2,1)=-0.025;
 a(4,0)=-0.7;  a(0,4)= 0.9; a(3,1)= 0.6; a(1,3)=-0.015; a(2,2)= 0.3;
 a(5,0)=-0.2;  a(0,5)= 0.3; a(4,1)=-0.6; a(1,4)= 0.150; a(3,2)=-0.2; a(2,3)=0.1;

 // Zero the vector
 // Now loop over coefficients
 for(unsigned i=0; i<6;++i)
  {
   for(unsigned j=0; j<6-i;++j)
    {
     // Function
     p[0]+=a(i,j)*pow(x[0],i)*pow(x[1],j);
     // Values
     if(i>0)
       p[1]+=a(i,j)*i*pow(x[0],i-1)*pow(x[1],j);
     if(j>0)
       p[2]+=a(i,j)*j*pow(x[0],i)*pow(x[1],j-1);
     // Second derivatives
     if(i>1)
       p[3]+=a(i,j)*i*(i-1)*pow(x[0],i-2)*pow(x[1],j);
     if(i>0 && j>0)
       p[4]+=a(i,j)*j*i*pow(x[0],i-1)*pow(x[1],j-1);
     if(j>1)
       p[5]+=a(i,j)*j*(j-1)*pow(x[0],i)*pow(x[1],j-2);
    }
  }
}

}

//#############################################################################//
// Main function - performs all the tests.
//#############################################################################//
int main(int argc, char **argv)
{
 // Use convenient namespace
 using namespace bell_element_tests; 
 feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

 // Set up command line arguments
 oomph::CommandLineArgs::setup(argc,argv);
 BellElementTestBasis mybasis;

 // Now define some vertices  
 Vector<Vector<double> > vertices(3,Vector<double>(2));
 vertices[0][0] = 1.2;
 vertices[0][1] = 0.3;
 vertices[1][0] = 0.2;
 vertices[1][1] = 0.8;
 vertices[2][0] = 0.6;
 vertices[2][1] =-0.3;
 mybasis.set_vertices(vertices);

 // Check on an element
 oomph_info <<"Checking on reference element."<<std::endl;
 oomph_info <<"Dofs of shape functions gives:"<<std::endl;
 // Check on the reference basis
 mybasis.check_local_basis_dofs();

 // Check the dofs on the vertices 
 oomph_info <<"Checking on element with vertices:"<<mybasis.get_vertices()<<std::endl;
 oomph_info <<"Dofs of shape functions gives:"<<std::endl;
 mybasis.check_basis_global_dofs();

 // Output some points on the reference triangle
 oomph_info <<"Some points on the reference triangle:"<<std::endl;
 const unsigned n_points_per_side = 6, dim = 2; 
 oomph_info<<"Graphics[{";
 for(unsigned i=0;i<n_points_per_side;++i)
  {
   for(unsigned j=0;j<n_points_per_side-i;++j)
    {
     Vector<double> s(dim,0.0);
     s[0] = i * 1.0/(n_points_per_side-1);
     s[1] = j * 1.0/(n_points_per_side-1);
     oomph_info << "Point[{"<<s[0]<<","<<s[1]<<"}]"<<((i==6-1 && j== 0) ? "}" :",");
    }
  }
 oomph_info<<"]\n"<<std::endl;

 // Output the same points on the physical triangle
 oomph_info <<"The same points on the physical triangle:"<<std::endl;
 oomph_info<<"Graphics[{";
 for(unsigned i=0;i<n_points_per_side;++i)
  {
   for(unsigned j=0;j<n_points_per_side-i;++j)
    {
     Vector<double> s(dim,0.0),x(dim,0.0);
     s[0] = i * 1.0/(n_points_per_side-1);
     s[1] = j * 1.0/(n_points_per_side-1);
     mybasis.interpolated_position(s,x);
     oomph_info << "Point[{"<<x[0]<<","<<x[1]<<"}]"<<((i==6-1 && j== 0) ? "}" :",");
    }
  }
 oomph_info<<"]\n"<<std::endl;


 // Check the differences with a known function
 oomph_info <<"Checking interpolation element with vertices:"<<mybasis.get_vertices()<<std::endl;
 oomph_info <<"L2 norm of p4 polynomial, which should be exact gives:";
 oomph_info <<mybasis.check_analytic_function(&TestSoln::get_p4)<<std::endl;
 oomph_info <<std::endl;

 // Check the differences with a known function
 oomph_info <<"Checking convergence of interpolation on element"<<std::endl;
 oomph_info <<"L2 norm of p5 polynomial is"<<std::endl;

 // Set up parameters of heirarchy
 const unsigned n_length_steps = 7, nnode =3;
 double length_scale = sqrt(3);
 vertices[0][0] = sqrt(3)/2;
 vertices[0][1] = 1./2;
 vertices[1][0] = -sqrt(3)/2;
 vertices[1][1] =  1./2;
 vertices[2][0] = 0;
 vertices[0][1] =-1 ;

 // Open trace
 std::ofstream trace_file;
 trace_file.open("RESLT/trace_bell.dat");
 // Loop over the the 'size' of the element 
 for(unsigned istep = 0 ; istep < n_length_steps;  ++istep)
  {
   // New length_scale
   Vector<Vector<double> > test_vertices(nnode,Vector<double>(dim));
   length_scale /= sqrt(2.0);
   // Fill in new vertices
   for(unsigned i=0;i<dim*nnode;++i)
     {test_vertices[i/2][i%2] = vertices[i/2][i%2]*length_scale;}
   // Update basis
   mybasis.set_vertices(test_vertices);
   // Now Test
   oomph_info <<length_scale<< " "<<sqrt(3)*pow(length_scale,2)/4.<<" "<<mybasis.check_analytic_function(&TestSoln::get_p5)<<std::endl;
   trace_file<<length_scale<< " "<<sqrt(3)*pow(length_scale,2)/4.<<" "<<mybasis.check_analytic_function(&TestSoln::get_p5)<<std::endl;
  } 
 trace_file.close();
}
