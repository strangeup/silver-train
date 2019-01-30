// A c++ program to test a single element

#include <fenv.h> 
//Generic routines
#include "generic.h"

// My New code to check
#include "single_element_tester.h"
#include "test_boundaries.h"

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
// Derivative of parametric function
void d2_parametric_edge_quad(const double& s, Vector<double>& dx)
{ dx[0] =0;  dx[1] =-2;};

// Parametric function for boundary part 1
void parametric_edge_quint(const double& s, Vector<double>& x)
{ x[0] =-s;/*-s*s*s;*/  x[1] =-s*s*s*(s*s - 1)+2;};
// Derivative of parametric function
void d_parametric_edge_quint(const double& s, Vector<double>& dx)
//{ dx[0] =-1;/*-3*s*s;*/  dx[1] =-s*s*s*s*(7*s*s - 5);};
{ dx[0] =-1;/*-3*s*s;*/  dx[1] =-s*s*(5*s*s - 3);};
// Derivative of parametric function
void d2_parametric_edge_quint(const double& s, Vector<double>& dx)
//{ dx[0] =0;/*-6*s;*/  dx[1] =-s*s*s*(35*s*s - 20) ;};
{ dx[0] =0;/*-6*s;*/  dx[1] =-s*s*(20*s) + 6*s;};

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

/// Exact solution as a Vector
void get_clamped_exact(const Vector<double>& x, Vector<double>& w)
{
 const double p_mag=64;
 w[0]= p_mag*pow(x[0]*x[0]+x[1]*x[1]-1,2)/64;
 w[1]= x[0]*p_mag*(x[0]*x[0]+x[1]*x[1]-1)/16;
 w[2]= x[1]*p_mag*(x[0]*x[0]+x[1]*x[1]-1)/16;
 w[3]= p_mag*(3*x[0]*x[0]+x[1]*x[1]-1)/16;
 w[5]= p_mag*x[0]*x[1]/8;
 w[4]= p_mag*(x[0]*x[0]+3*x[1]*x[1]-1)/16;
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

void get_p8(const Vector<double>& x_alpha, Vector<double>& p)
{
 const double x=x_alpha[0], y=x_alpha[1];
 // Copy of the position Vector
 p[0] =  pow(x,4)*pow(y,4);
 p[1] =  4*pow(x,3)*pow(y,4);
 p[2] =  4*pow(x,4)*pow(y,3);
 p[3] =  12*pow(x,2)*pow(y,4);
 p[5] =  16*pow(x,3)*pow(y,3);
 p[4] =  12*pow(x,4)*pow(y,2);
}

void get_p8_on_radius(const Vector<double>& x_alpha, Vector<double>& p)
{
 const double x=x_alpha[0], y=x_alpha[1], r=sqrt(x*x+y*y), t=atan2(x,y);
 // Copy of the position Vector
 p[0] =  pow(r*cos(t),4)*pow(r*sin(t),4);
 p[1] =  8*pow(r,7)*pow(cos(t),4)*pow(sin(t),4);
 p[2] = -4*pow(r,8)*pow(cos(t),3)*pow(sin(t),5)
        +4*pow(r,8)*pow(cos(t),5)*pow(sin(t),3);
 p[3] =  56*pow(r,6)*pow(cos(t),4)*pow(sin(t),4);
 p[4] = -32*pow(r,7)*pow(cos(t),3)*pow(sin(t),5)
        +32*pow(r,7)*pow(cos(t),5)*pow(sin(t),3);
 p[5] = +12*pow(r,8)*pow(cos(t),2)*pow(sin(t),6)
        -16*pow(r,8)*pow(cos(t),4)*pow(sin(t),4)
        -16*pow(r,8)*pow(cos(t),4)*pow(sin(t),4)
        +12*pow(r,8)*pow(cos(t),6)*pow(sin(t),2);
}

void get_solution_and_radial_derivative
 (const Vector<double>& x_alpha, Vector<double>& w)
{
 // Copy of the position Vector
 // Solution r^4 Cos(3\theta) /35
 const double x=x_alpha[0], y=x_alpha[1];

 // Check for r = 0
 if(fabs(x)<1e-14 && fabs(y)<1e-14)
  {
   // All zero at centre
   w=Vector<double>(6,0.0);
   return;
  }
 // Otherwise
  w[0] = sqrt(x*x+y*y)*(pow(x,3)-3*x*y*y);

  w[1] = x*(pow(x, 3) - 3*x*pow(y, 2))*std::pow(pow(x, 2) + pow(y, 2), -0.5) +
(3*pow(x, 2) - 3*pow(y, 2))*std::pow(pow(x, 2) + pow(y, 2), 0.5);

  w[2] = y*(pow(x, 3) - 3*x*pow(y, 2))*std::pow(pow(x, 2) + pow(y, 2), -0.5) - 
6*x*y*pow(pow(x, 2) + std::pow(y, 2), 0.5); 

  w[3] =2*x*(3*pow(x, 2) - 3*pow(y, 2))*std::pow(pow(x, 2) + pow(y, 2), -0.5) + 
(pow(x, 3) - 3*x*pow(y, 2))*(-(pow(x, 2)*std::pow(pow(x, 2) + pow(y, 2),-1.5)) 
+ std::pow(pow(x, 2) + pow(y, 2), -0.5)) + 6*x*std::pow(pow(x, 2) + pow(y, 2), 0.5);

  w[5] =-(x*y*(pow(x, 3) - 3*x*pow(y, 2))*std::pow(pow(x, 2) + pow(y, 2), -1.5)) - 
6*y*pow(x, 2)*std::pow(pow(x, 2) + pow(y, 2), -0.5) + y*(3*pow(x, 2) - 
3*pow(y, 2))*std::pow(pow(x, 2) + pow(y, 2), -0.5) - 6*y*std::pow(pow(x, 2) + 
pow(y, 2), 0.5);

  w[4] =-12*x*pow(y, 2)*std::pow(pow(x, 2) + pow(y, 2), -0.5) + (pow(x, 3) - 
3*x*pow(y, 2))*(-(pow(y, 2)*std::pow(pow(x, 2) + pow(y, 2), -1.5)) + 
std::pow(pow(x, 2) + pow(y, 2), -0.5)) - 6*x*std::pow(pow(x, 2) + pow(y, 2), 
0.5); 

 // Get the radial derivative 
// w[6] = 4*x*(x*x+y*y);

}

void get_sine_cosine(const Vector<double>& x, Vector<double>& w)
{
 w[0]= sin(x[0]*Pi)*cos(x[1]*Pi);
 w[1]= cos(x[0]*Pi)*cos(x[1]*Pi)*Pi;
 w[2]=-sin(x[0]*Pi)*sin(x[1]*Pi)*Pi;
 w[3]=-sin(x[0]*Pi)*cos(x[1]*Pi)*Pi*Pi;
 w[4]=-cos(x[0]*Pi)*sin(x[1]*Pi)*Pi*Pi;
 w[5]=-sin(x[0]*Pi)*cos(x[1]*Pi)*Pi*Pi;
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

// Get the tip of an isoceles triangle
void get_isosceles_tip(const double& s, const Vector<double>& x1, 
const Vector<double>& x2, Vector<double>& x3)
 {
  // Definitions
  const double x1x2_mag_2 =pow(x2[0]-x1[0],2)+pow(x2[0]-x1[0],2);
  const double alt3 = sqrt(s*s-(x1x2_mag_2)/4.);
  // Fill in x3
  x3[0] = (x1[0] + x2[0]) / 2. - alt3 * (x2[1]-x1[1]) / sqrt(x1x2_mag_2) ;
  x3[1] = (x1[1] + x2[1]) / 2. + alt3 * (x2[0]-x1[0]) / sqrt(x1x2_mag_2) ;
 }

// Get the tip of an isoceles triangle
void get_right_handed_equilateral_vertex(const Vector<double>& x1, 
const Vector<double>& x2, Vector<double>& x3)
 {
  // x component
  x3[0] = 0.5*(x2[0]+x1[0])-sqrt(3)/2.*(x2[1]-x1[1]);
  x3[1] = 0.5*(x2[1]+x1[1])-sqrt(3)/2.*(x1[0]-x2[0]);
 }

// Get the tip of an isoceles triangle
double distance_between_two_points(const Vector<double>& x1, 
const Vector<double>& x2)
 {
  // x1. x2
 return std::sqrt(std::pow(x1[0]-x2[0],2)+std::pow(x1[1]-x2[1],2));
 }

//#############################################################################//
// Main function - performs all the tests.
//#############################################################################//
int main(int argc, char **argv){
 using namespace c1_curved_checks; 
 //set up
 oomph::CommandLineArgs::setup(argc,argv);
 oomph::CommandLineArgs::specify_command_line_flag("--higher_order");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 

 // Higher order flag - this code is a bit dodgy
 const bool do_higher_order = CommandLineArgs::command_line_flag_has_been_set("--higher_order");

 // Set up the parametric function
 CurvilineCircleRight parametric_curve;
 // Constants for the hierarchy of triangles 
 const double centre = 0.8 , s_start = 0.4;//Pi/8.;

 // Open a mathematica module
 std::ofstream outfile;
 char filename[100];
 sprintf(filename,"RESLT/tests.m");
 outfile.open(filename);
 // Now output the p5 basis (hard coded) to mathematica script
 outfile << "(* p5 basis polynomials *)\n";
 outfile << "p5[y_]:= { \
 {-(-1 + y)^3 (1 + 3 y + 6 y^2)},\
 {y^3 (10 - 15 y + 6 y^2)},\
 {-(-1 + y)^3 y (1 + 3 y)},\
 {-(-1 + y) y^3 (-4 + 3 y)},\
 {-(1/2) (-1 + y)^3 y^2},\
 {1/2 (-1 + y)^2 y^3} \
}\n";

 // Find the parametric dofs for the mathematica script
 Vector<Vector<double> > parametric_dofs(6,Vector<double> (2));
 parametric_curve.position(Vector<double>(1,0),parametric_dofs[0]);
 parametric_curve.position(Vector<double>(1,1),parametric_dofs[1]);
 parametric_curve.dposition(Vector<double>(1,0),parametric_dofs[2]);
 parametric_curve.dposition(Vector<double>(1,1),parametric_dofs[3]);
 parametric_curve.d2position(Vector<double>(1,0),parametric_dofs[4]);
 parametric_curve.d2position(Vector<double>(1,1),parametric_dofs[5]);

 // Initialize graphics in mathematica plotting script
 outfile << "<<JavaGraphics`\n";
 // Now output the exact chi and dchi
 outfile << "(* The Edge Parametric function *)\n";
 outfile << "\nchi1 = {";
 for(unsigned i=0; i<6; ++i)
  {outfile << parametric_dofs[i][0] << (i==5?"}":","); } 

 outfile << ".p5[s] \nchi2 = {";
 for(unsigned i=0; i<6; ++i)
  {outfile << parametric_dofs[i][1] << (i==5?"}":","); } 
 outfile << ".p5[s] \n";
 
 // Establish an order of  colours in the script so that it's appealing 
 Vector<std::string> color_order(7);
 color_order[0] ="Red";
 color_order[1] ="Orange";
 color_order[2] ="Yellow";
 color_order[3] ="Green";
 color_order[4] ="Blue";
 color_order[5] ="Purple";
 color_order[6] ="Magenta";
 outfile << "Show[{";
 outfile << "Graphics[{";

 // Initialise number of elements
 const unsigned nel =7;
 // Open a new file
 std::ofstream trace_file;
 trace_file.open("RESLT/trace.dat");
 // Now loop over the elements
 for(unsigned iel=0 ; iel< nel ; ++iel)
  {
   // Definitions 
   double delta_s = s_start/std::pow(3./2.,iel+1);
   double s_ubar(centre - delta_s/2.), s_obar(centre + delta_s/2.);

   // Now get the vertices and then plot the elements
   Vector<Vector<double> > vertices(3,Vector<double>(2,0.0));
   parametric_curve.position(Vector<double>(1,s_ubar),vertices[0]);
   parametric_curve.position(Vector<double>(1,s_obar),vertices[1]);

   // Get the length of the sides and the final vertice to construct
   // the equilateral triangle
   double side_length = distance_between_two_points(vertices[0],vertices[1]);
   get_right_handed_equilateral_vertex(vertices[0],vertices[1],vertices[2]);
//   get_isosceles_tip(delta_s, vertices[0],vertices[1],vertices[2]);
   // Now output
   outfile <<"Opacity[0.4],"<< color_order[iel % 7] <<",";
   output_polygon_to_mathematica(outfile,3,vertices);
   outfile <<"\n" <<(iel==nel-1?"}]":","); 

   // Set up the problem
   //BernadouElementTestBasis<5> test;
   //   if(iel==0)
   //    test_basis_3.check_d_matrix(&get_p8);
 
    // Do some checks on the constants so we know they are functioning correctly
    //    test_basis_3.check_constant_consistency();
    if(do_higher_order)
    {
     // Construct the element
     BernadouElementTestBasis<5> test_basis_5;
     test_basis_5.upgrade_element(vertices,s_ubar,s_obar,parametric_curve);
     oomph_info<<(iel==0?"The answer is:\n[":"")
       //<< delta_s <<","<<
       << side_length <<","<<
          test_basis_5.check_function_norm(&get_solution_and_radial_derivative) 
       //mapping_on_trace_norm()
        <<(iel==nel-1? "]\n":";");
    // Output to trace file
    trace_file << side_length <<" "<<
        test_basis_5.check_function_norm(&get_solution_and_radial_derivative) 
        <<(iel==nel-1? "":"\n");
    }
    else
    {
     // Construct the element
     BernadouElementTestBasis<3> test_basis_3;
     test_basis_3.upgrade_element(vertices,s_ubar,s_obar,parametric_curve);
     oomph_info<<(iel==0?"The answer is:\n[":"")<< side_length <<","
         <<test_basis_3.check_function_norm(&get_solution_and_radial_derivative) 
         <<(iel==nel-1? "]\n":";");
    // Output to trace file
    trace_file << side_length <<" "
         <<test_basis_3.check_function_norm(&get_solution_and_radial_derivative) 
         <<(iel==nel-1? "":"\n");
    }
  }

 outfile << ",\nParametricPlot[Join[chi1,chi2],{s,"<<
   -s_start + centre<<","<<s_start + centre<<"}]\n";
 // End show
 outfile << "}]";
 // Close file
 outfile.close();
 trace_file.close();
}

 

