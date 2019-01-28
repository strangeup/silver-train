// A c++ program to test a single element

#include <fenv.h> 
//Generic routines
#include "generic.h"

// My New code to check
#include "single_element_tester.h"

using namespace oomph;

using namespace MathematicalConstants;

//#############################################################################//
// Parametric Boundaries
//#############################################################################//
// Parametric function for boundary part 1
void parametric_edge_quad(const double& s, Vector<double>& x)
{ x[0] =-s;  x[1] =-s*s+0.75;};
// Derivative of parametric function
void d_parametric_edge_quad(const double& s, Vector<double>& dx)
{ dx[0] =-1;  dx[1] =-2*s;};

// Parametric function to describe the single (straight) curved side
void chi_0(const double& s, Vector<double>& chi)
{
 // x equal to s
 chi[0]=s;

 // y equal to Sqrt[3]/2
 chi[1]=sqrt(3)/2.;
}

void d_chi_0(const double& s, Vector<double>& d_chi)
{
 // x equal to s
 d_chi[0]=1;

 // y equal to Sqrt[3]/2
 d_chi[1]=0;
}

// Parametric function to describe the single (quadratic) curved side
void chi_2(const double& s, Vector<double>& chi)
{
 // x equal to s
 chi[0]=s;

 // y equal to simplest quadratic that I could think of
 chi[1]=-s*s+1/4.+sqrt(3)/2.;
}

// Parametric function to describe the single (quadratic) curved side
void d_chi_2(const double& s, Vector<double>& d_chi)
{
 // derivative of y component 
 d_chi[0]=1;

 // derivative of x component
 d_chi[1]=-2*s;
}

/// Parametric function to describe the single (cubic) curved side
void chi_circle(const double& s, Vector<double>& chi)
{
 // x equal to s
 chi[0]=std::cos(s);

 // y equal to simplest cubic I could think of 
 chi[1]=std::sin(s);
}

/// Parametric function to describe the single (cubic) curved side
void d_chi_circle(const double& s, Vector<double>& d_chi)
{
 // x equal to s
 d_chi[0]=-std::sin(s);

 // y equal to simplest cubic I could think of 
 d_chi[1]= std::cos(s);
}
/// Parametric function to describe the single (cubic) curved side
void chi_3(const double& s, Vector<double>& chi)
{
 // x equal to s
 chi[0]=s;

 // y equal to simplest cubic I could think of 
 chi[1]=s*(s*s-1/4.)+sqrt(3)/2.;
}

/// Parametric function to describe the single (cubic) curved side
void d_chi_3(const double& s, Vector<double>& d_chi)
{
 // x equal to s
 d_chi[0]=1;

 // y equal to simplest cubic I could think of 
 d_chi[1]=3*s*s-1/4.;
}

/// Exact solution as a Vector
void get_p0(const Vector<double>& x, Vector<double>& p)
{
 // Copy of the position Vector
 const double a0=1.0;
 Vector <double> xp(x); 
 p[0]=a0;
}

/// Exact solution as a Vector
void get_p1(const Vector<double>& x, Vector<double>& p)
{
 // Copy of the position Vector
 const double a00=0.2, a01=1.0, a10=- 0.4;
 Vector <double> xp(x); 
 p[0]=a00 + a10*x[0] + a01*x[1];
 p[1]=a10;
 p[2]=a01;
}

//Exact solution for constant pressure, circular domain and resting boundary conditions
void get_p2_clamp(const Vector<double>& x, Vector<double>& w)
{
 w[0]= +x[0]*x[0]-3./4.+x[1];
 w[1]= 2*x[0];
 w[2]= 1.0;
 w[3]= 2;
 w[4]= 0;
 w[5]= 0;
}

void get_p2(const Vector<double>& x, Vector<double>& p)
{
 // Copy of the position Vector
 DenseMatrix<double> a(3,3,0.0);
 a(0,0)= 0.2; a(0,1)= 1.0; a(1,0)=- 0.4;
 a(1,1)= 0.2, a(2,0)=-0.9; a(0,2)=-0.1;

 // Zero the vector
 // Now loop over coefficients
 for(unsigned i=0; i<3;++i)
  {
   for(unsigned j=0; j<3-i;++j)
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
     if(j>1)
       p[4]+=a(i,j)*j*(j-1)*pow(x[0],i)*pow(x[1],j-2);
     if(i>0 && j>0)
       p[5]+=a(i,j)*j*i*pow(x[0],i-1)*pow(x[1],j-1);
    }
  }
}

void get_p3(const Vector<double>& x, Vector<double>& p)
{
 // Copy of the position Vector
 DenseMatrix<double> a(4,4,0.0);
 a(0,0)= 0.2; a(0,1)= 1.0; a(1,0)=-0.4; a(1,1)=0.2; 
 a(2,0)=-0.9; a(0,2)=-0.1; a(3,0)= 0.4;
 a(0,3)= 0.6; a(1,2)= 0.3; a(2,1)=-0.025;
 // Zero the vector
 // Now loop over coefficients
 for(unsigned i=0; i<4;++i)
  {
   for(unsigned j=0; j<4-i;++j)
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
     if(j>1)
       p[4]+=a(i,j)*j*(j-1)*pow(x[0],i)*pow(x[1],j-2);
     if(i>0 && j>0)
       p[5]+=a(i,j)*j*i*pow(x[0],i-1)*pow(x[1],j-1);
    }
  }
}

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
     if(j>1)
       p[4]+=a(i,j)*j*(j-1)*pow(x[0],i)*pow(x[1],j-2);
     if(i>0 && j>0)
       p[5]+=a(i,j)*j*i*pow(x[0],i-1)*pow(x[1],j-1);
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
 a(4,0)=-0.7;  a(0,4)= 0.9; a(3,1)= 0.6; a(1,3)=-0.015; a(2,2)=0.3;
 a(5,0)=-0.7;  a(0,5)= 0.9; a(4,1)= 0.6; a(1,4)=-0.015; a(3,2)=0.3; a(2,3)=0.3;

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
     if(j>1)
       p[4]+=a(i,j)*j*(j-1)*pow(x[0],i)*pow(x[1],j-2);
     if(i>0 && j>0)
       p[5]+=a(i,j)*j*i*pow(x[0],i-1)*pow(x[1],j-1);
    }
  }
}

/// Exact solution as a Vector
void get_radialp4(const Vector<double>& x, Vector<double>& p)
{
 p[0] = std::pow(1.0-x[0]*x[0]-x[1]*x[1],2) ;
 p[1] = -4*x[0]*(1.0-x[0]*x[0]-x[1]*x[1]);
 p[2] = -4*x[1]*(1.0-x[0]*x[0]-x[1]*x[1]);
 p[3] = -4*(1.0-3*x[0]*x[0]-x[1]*x[1]);
 p[4] = -4*(1.0-x[0]*x[0]-3*x[1]*x[1]);
 p[5] = 8*x[1]*x[0];
}


// Parametric function for boundary part 0
// We want to make sure s=2Pi is as far away as possible from the boundary
void parametric_edge_0(const double& s, Vector<double>& x)
{x[0] = std::cos(s-Pi/2);   x[1] =std::sin(s-Pi/2); };

// Derivative of parametric function
void d_parametric_edge_0(const double& s, Vector<double>& dx)
{dx[0] =-std::sin(s-Pi/2.);  dx[1] =std::cos(s-Pi/2);};

// Parametric function for boundary part 1
void parametric_edge_1(const double& s, Vector<double>& x)
{x[0] =-s;  x[1] = 0.5;};
// Derivative of parametric function
void d_parametric_edge_1(const double& s, Vector<double>& dx)
{dx[0] =-1;  dx[1] = 0.0;};
// Get s from x
double get_s_1(const Vector<double>& x){return -x[0];};

// Typedef
typedef void (*Parametric_fct_pt)(const double& s, Vector<double>& v);

//#############################################################################//
// Main function - performs all the tests.
//#############################################################################//
int main(int argc, char **argv)
{
 // Use convenient namespace
 using namespace c1_curved_checks; 
 // Set up command line arguments
 oomph::CommandLineArgs::setup(argc,argv);
 
 // Set up the parametric function
 CurvilineCircleTop parametric_curve;
 double s_ubar(-0.5), s_obar(0.5);

 // Vertices for test triangle
 Vector<Vector<double> > vertices(3,Vector<double>(2,0.0));
 (*chi_pt)(s_ubar, vertices[0]);
 (*chi_pt)(s_obar, vertices[1]);
 vertices[2][0] = 0.0; 
 vertices[2][1] = 0.2;

 // Output the vertices
 oomph_info<<"Vertices"<<vertices<<"\n";

 // Set up the problem
 BernadouElementTestBasis bernadou_test_basis;
 bernadou_test_basis.upgrade_element(vertices,s_ubar,s_obar,parametric_curve);

 // Output the element to Mathematica
 bernadou_test_basis.output_to_mathematica_graphics();

 // Do some rudimentary checks and the Jacobian and Hessian
 bernadou_test_basis.check_jacobian_and_hessian();
 
 // Do some checks on the constants so we know they are functioning correctly
 bernadou_test_basis.check_constant_consistency();

 // Check we get the identity matrix on this bit
 oomph_info<<"\nDofs of shape functions (and dshape), should gve the identity"
         <<"matrix. Any difference will be recorded.";
 bernadou_test_basis.check_basic_shape(1e-12);

 // Now check the submatrix 
 bernadou_test_basis.check_b2l_submatrix_3(1e-20);

 // Check 1d hermite shape functions
 bernadou_test_basis.check_1d_hermite_shape(1e-15);

 // Check traces of the shape functions
 oomph_info<<"Check a radial function (p4)\n";
 bernadou_test_basis.check_traces(get_radialp4, 1e-16,false,false);

 oomph_info<<"Check a degree 5 bivariate polynomial (p5)\n";
 bernadou_test_basis.check_traces(get_p5, 1e-16,false,true);

 oomph_info<<"Check a degree 4 bivariate polynomial (p4)\n";
 bernadou_test_basis.check_traces(get_p4, 1e-16,false,false);

 oomph_info<<"Check a degree 3 bivariate polynomial (p3)\n";
 bernadou_test_basis.check_traces(get_p3, 1e-16,false,false);

 oomph_info<<"Check a degree 2 bivariate polynomial (p2)\n";
 bernadou_test_basis.check_traces(get_p2, 1e-16,false,false);

 oomph_info<<"Check a degree 1 bivariate polynomial (p1)\n";
 bernadou_test_basis.check_traces(get_p1, 1e-16,false,false);

 oomph_info<<"Check a degree 0 bivariate polynomial (p0)\n";
 bernadou_test_basis.check_traces(get_p0, 1e-16,false,false);

 // Check g3 trace
 oomph_info<<"\nCheck g3 trace:\n";
 // This has it's own check
 bernadou_test_basis.check_g3_trace(0.0);

 oomph_info<<"\nCheck f3 trace:\n";
 bernadou_test_basis.check_f3_trace(0.0);

 // Check we can construct the Submatrices
 oomph_info<<"\nCheck delta property:\n";
 bernadou_test_basis.check_shape(1e-16);
 oomph_info<<"\nCheck interpolation of p0:\n";
 bernadou_test_basis.check_function(get_p0,1e-16);
 oomph_info<<"\nCheck interpolation of p1:\n";
 bernadou_test_basis.check_function(get_p1,1e-16);
 oomph_info<<"\nCheck interpolation of p2:\n";
 bernadou_test_basis.check_function(get_p2,1e-16);
 oomph_info<<"\nCheck interpolation of p3:\n";
 bernadou_test_basis.check_function(get_p3,1e-16);
 oomph_info<<"\nCheck interpolation of p4:\n";
 bernadou_test_basis.check_function(get_p4,1e-16);
 oomph_info<<"\nCheck interpolation of radial p4:\n";
 bernadou_test_basis.check_function(get_radialp4,1e-16);

 // Check get shape
 bernadou_test_basis.check_get_shape(1e-11);

 // Check we can interpolate p0 - p5
 oomph_info<<"\nCheck interpolation of p0:\n";
 bernadou_test_basis.check_basis_function(get_p0,1e-16);

 oomph_info<<"\nCheck interpolation of p1:\n";
 bernadou_test_basis.check_basis_function(get_p1,1e-16);

 oomph_info<<"\nCheck interpolation of p2:\n";
 bernadou_test_basis.check_basis_function(get_p2,1e-16);

 oomph_info<<"\nCheck interpolation of p3:\n";
 bernadou_test_basis.check_basis_function(get_p3,1e-16);

 oomph_info<<"\nCheck interpolation of p4:\n";
 bernadou_test_basis.check_basis_function(get_p4,1e-16);

 oomph_info<<"\nCheck interpolation of radial p4:\n";
 bernadou_test_basis.check_basis_function(get_radialp4,1e-16);
}

 

