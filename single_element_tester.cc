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

// Helper
void print_at_three_points(std::string fcthead, const  Parametric_fct_pt& chi_p
  , const double& s_ubar, const double& s_obar)
{
 Vector<double> v_l(2,0.0), v_m(2,0.0),v_u(2,0.0);
 chi_p(s_ubar,v_l);
 chi_p((s_obar+s_ubar)/2,v_m);
 chi_p(s_obar,v_u);
 // Now call the functions at 0 and 1
 std::cout << fcthead << v_l<< " ; "<< v_u 
           << " ; "  << v_m <<"\n";
}
//#############################################################################//
// Main function - performs all the tests.
//#############################################################################//
int main(int argc, char **argv){
 using namespace c1_curved_checks; 
 //set up
 oomph::CommandLineArgs::setup(argc,argv);
 
 // Set up the parametric function
 void (*chi_pt)(const double& s, Vector<double>& v) = &parametric_edge_quad;
 void (*d_chi_pt)(const double& s, Vector<double>& v) = &d_parametric_edge_quad;
 double s_ubar(-0.5), s_obar(0.5);

 // Vertices for test triangle
 Vector<Vector<double> > vertices(3,Vector<double>(2,0.0));
 (*chi_pt)(s_ubar, vertices[0]);
 (*chi_pt)(s_obar, vertices[1]);
 vertices[2][0] = 0.0; 
 vertices[2][1] = 0.2;

 // ORDER OF BOUNDARY
 const unsigned boundary_order =3;

 // Do all the tests - even though we expect some to fail.
 bool do_all_tests=true;

 std::cout<<"Vertices"<<vertices<<"\n";

 // Now call the functions at 0 and 1
 print_at_three_points("X0 ",&chi_0,s_ubar,s_obar);
 print_at_three_points("X0'",&d_chi_0,s_ubar,s_obar);
 print_at_three_points("X2 ",&chi_2,s_ubar,s_obar);
 print_at_three_points("X2'",&d_chi_2,s_ubar,s_obar);
 print_at_three_points("X3 ",&chi_3,s_ubar,s_obar);
 print_at_three_points("X3'",&d_chi_3,s_ubar,s_obar);
// std::cout << "X0 "  << chi_0(s_ubar)<< " ; "<< chi_0(s_obar) 
//           << " ; "  << chi_0((s_obar+s_ubar)/2)<<"\n";
// std::cout << "X0'"  << d_chi_0(s_ubar)<< " ; "<< d_chi_0(s_obar)
//           << " ; "  << d_chi_0((s_obar+s_ubar)/2)<<"\n";
//
// std::cout << "X2 "  << chi_2(s_ubar)<< " ; "<< chi_2(s_obar)
//           << " ; "  << chi_2((s_obar+s_ubar)/2)<<"\n";
// std::cout << "X2'"  << d_chi_2(s_ubar)<< " ; "<< d_chi_2(s_obar)
//           << " ; "  << d_chi_2((s_obar+s_ubar)/2)<<"\n";
//
// std::cout << "X3 "  << chi_3(s_ubar)<< " ; "<< chi_3(s_obar) 
//           << " ; "  << chi_3((s_obar+s_ubar)/2)<<"\n";
// std::cout << "X3'"  << d_chi_3(s_ubar)<< " ; "<< d_chi_3(s_obar)
//           << " ; "  << d_chi_3((s_obar+s_ubar)/2)<<"\n";

 // Set up the problem
 CurvilineCircleTop parametric_curve;
 CurvedElementChecker test; //(vertices,s_ubar,s_obar);
 test.upgrade_element(vertices,s_ubar,s_obar,parametric_curve);
// test.get_chi_fct_pt() = chi_pt;
// test.get_d_chi_fct_pt() = d_chi_pt;

 // Output the element to Mathematica
 test.output_to_mathematica_graphics();

 // Do some rudimentary checks and the Jacobian and Hessian
 test.check_jacobian_and_hessian();
 
 // Do some checks on the constants so we know they are functioning correctly
 test.check_constant_consistency();

 // Check we get the identity matrix on this bit
 std::cout<<"\nDofs of shape functions (and dshape), should gve the identity"
         <<"matrix. Any difference will be recorded.";
 test.check_basic_shape(1e-12);

 // Now check the submatrix 
 test.check_b2l_submatrix_3(1e-20);


 // Check 1d hermite shape functions
 test.check_1d_hermite_shape(1e-15);

 // Check traces of the shape functions
 std::cout<<"Check a radial function (p4)\n";
 test.check_traces(get_radialp4, 1e-16,false,false);

 std::cout<<"Check a degree 5 bivariate polynomial (p5)\n";
 test.check_traces(get_p5, 1e-16,false,true);

 std::cout<<"Check a degree 4 bivariate polynomial (p4)\n";
 test.check_traces(get_p4, 1e-16,false,false);

 std::cout<<"Check a degree 3 bivariate polynomial (p3)\n";
 test.check_traces(get_p3, 1e-16,false,false);

 std::cout<<"Check a degree 2 bivariate polynomial (p2)\n";
 test.check_traces(get_p2, 1e-16,false,false);

 std::cout<<"Check a degree 1 bivariate polynomial (p1)\n";
 test.check_traces(get_p1, 1e-16,false,false);

 std::cout<<"Check a degree 0 bivariate polynomial (p0)\n";
 test.check_traces(get_p0, 1e-16,false,false);

 // Check g3 trace
 std::cout<<"\nCheck g3 trace:\n";
 // This has it's own check
 test.check_g3_trace(0.0);

 std::cout<<"\nCheck f3 trace:\n";
 test.check_f3_trace(0.0);

 // Check we can construct the Submatrices
 std::cout<<"\nCheck delta property:\n";
 test.check_shape(1e-16);
 std::cout<<"\nCheck interpolation of p0:\n";
 test.check_function(get_p0,1e-16);
 std::cout<<"\nCheck interpolation of p1:\n";
 test.check_function(get_p1,1e-16);
 std::cout<<"\nCheck interpolation of p2:\n";
 test.check_function(get_p2,1e-16);
 std::cout<<"\nCheck interpolation of p3:\n";
 test.check_function(get_p3,1e-16);
 std::cout<<"\nCheck interpolation of p4:\n";
 test.check_function(get_p4,1e-16);
 std::cout<<"\nCheck interpolation of radial p4:\n";
 test.check_function(get_radialp4,1e-16);

 // Check get shape
 test.check_get_shape(1e-11);

 // Check we can interpolate p0 - p5
 std::cout<<"\nCheck interpolation of p0:\n";
 test.check_basis_function(get_p0,1e-16);

 std::cout<<"\nCheck interpolation of p1:\n";
 test.check_basis_function(get_p1,1e-16);

 std::cout<<"\nCheck interpolation of p2:\n";
 test.check_basis_function(get_p2,1e-16);

 std::cout<<"\nCheck interpolation of p3:\n";
 test.check_basis_function(get_p3,1e-16);

 std::cout<<"\nCheck interpolation of p4:\n";
 test.check_basis_function(get_p4,1e-16);

 std::cout<<"\nCheck interpolation of radial p4:\n";
 test.check_basis_function(get_radialp4,1e-16);
}

 

