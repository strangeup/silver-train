#ifndef C1_CURVED_CHECKS_HEADER
#define C1_CURVED_CHECKS_HEADER

#include <fenv.h>
//Generic routines
#include "generic.h"

// My New code to check
#include "C1_curved_elements.h"

namespace oomph {
namespace c1_curved_checks {
// Output polygon as mathematica code
void output_polygon_to_mathematica(std::ostream& os, const unsigned& n_vertices,
const Vector<Vector <double> >& vertices)
{
 os <<"Polygon[{\n";
 // Output all of the vertices
 for (unsigned i=0;i<n_vertices;++i)
  {os << "{"<<vertices[i][0]<<","<<vertices[i][1]<<"}"
      <<(i+1==n_vertices ?"\n":",");}
 os<<"}]";
}

// Output polygon as mathematica code
void output_arrow_to_mathematica(std::ostream& os, const Vector<double> start, 
 const Vector<double> end)
{
 os <<"Arrow[{\n";
 // Output all of the vertices
 os<<"{"<<start[0]<<","<<start[1]<<"},{"<<end[0]<<","<<end[1]<<"}\n}]";
}

/// 6th Order Accurate First Finite Difference
// (centred on the point fpt[2] with with spacing h)
double dfdx_fd4(const Vector<double>& fpt, const double& h)
 {
  // Central differences with stencil over 5 pts
  return (-(fpt[4]-fpt[0])/12. + 2.*(fpt[3]-fpt[1])/3.) / h;
 }

/// 4th Order Accurate Second Finite Difference
// (centred on the point fpt[3] with with spacing h)
double d2fdx2_fd6(const Vector<double>& fpt, const double& h)
 {
  // Central differences with stencil over 7 pts
  return (+(fpt[0]+fpt[6])/90.-3*(fpt[1]+fpt[5])/20. + 3.*(fpt[2]+fpt[4])/2.
          -49./18.*fpt[3]) /(h*h);
 }

/// 6th Order Accurate First Finite Difference
// (centred on the point fpt[3] with with spacing h)
double dfdx_fd6(const Vector<double>& fpt, const double& h)
 {
 // Central differences with stencil over 7 pts
 return (+(fpt[6]-fpt[0])/60. - 3.*(fpt[5]-fpt[1])/20. + 3.*(fpt[4]-fpt[2])/4.)
         / h;
 }

/// 4th Order Cross derivative, coefficients exact for fifth order polynomials
double d2fdxdy_fd4_coeff(const int m, const int n)
 {
  double coeff;
  // Fill in the coeffcients symmetrically
  if (abs(m)==1 && abs(n) == 1)
   {coeff = 4./9.;}
  else if (abs(m)==1 && abs(n) == 2)
   {coeff =-1./18.;}
  else if (abs(m)==2 && abs(n) == 1)
   {coeff =-1./18.;}
  else if (abs(m)==2 && abs(n) == 2)
   {coeff = 1./144.;}
  else
   {return 0;}
 // Invert the top left and bottom right
 if(m*n<0)
  {coeff*=-1;}
 return coeff;
 }

/// 4th Order Cross derivative, exact for fifth order polynomials
double d2fdxdy_fd4(const DenseMatrix<double>& fpt, const double& h)
 {
 double value=0.0;
  for(unsigned i=0;i<5;++i)
   {
    for(unsigned j=0;j<5;++j)
     {
      // Sum up contributions
      value+=d2fdxdy_fd4_coeff(int(i-2),int(j-2))*fpt(i,j)/(h*h);
     }
   }
  return value;
 }

//#############################################################################//
// New types and shorthands
//#############################################################################//
// New type of function pointer
typedef void (*ExactSolnFctPt)(const Vector<double>& x, Vector<double>& p);

// 1D analytic functions
typedef void (*analytic_1d_function)(const double& x, Vector<double>& p);

//#############################################################################//
// Generic function definitions
//#############################################################################//
// Two by two inverse
DenseMatrix<double> invert_two_by_two(const DenseMatrix<double>& jac)
 {
  //Initialise inverse and det
  DenseMatrix<double> inverse(2,2,0.0);
  double det=jac(0,0)*jac(1,1)-jac(1,0)*jac(0,1);
  // Now inverse
  inverse(0,0)= jac(1,1)/det;
  inverse(1,1)= jac(0,0)/det;
  inverse(0,1)=-jac(0,1)/det;
  inverse(1,0)=-jac(1,0)/det;
  return inverse;
 }



//#############################################################################//
// Checker class
//#############################################################################//
class CurvedElementChecker : public MyC1CurvedElements::TestElement<3>

{
public:
// Set edge to be two - so that no rotation etc. has to occur
  CurvedElementChecker(const VertexList& verts, const double& su, const double&so)
    : MyC1CurvedElements::TestElement<3>(verts,su,so) { set_edge(two);}
//#############################################################################//
// the element specific function definitions
//#############################################################################//
// Return the global dofs of a function on an element
Vector<double> get_global_dofs(const ExactSolnFctPt& w)
 {
  //Initialise
  Vector<Vector<double> > ai(3,Vector<double>(2,0.0)),
                          ei(3,Vector<double>(2,0.25)),
                          w_at_xi(6,Vector<double> (6,0.0));
  ai=get_vertices();
  Vector<double> gdofs(21,0.0);


  // Fill in ei
  ei[0][0]=0.5;
  ei[1][1]=0.5;

  // Get the dofs
  for (unsigned i=0; i<3; ++i)
   {
    // Fill in Vectors on nodes
    (*w)(ai[i],w_at_xi[i]);

    // Fill in Vectors at internal points
    Vector<double> x(2); f_k(ei[i],x);
    (*w)(x,w_at_xi[3+i]);
   }

  // Now rearrange them
  for (unsigned i=0; i< 3; ++i)
   {
   // First three dofs are on nodes
   gdofs[i]=w_at_xi[i][0];

   // Now fill in the next six first derivatives
   gdofs[3+2*i]=w_at_xi[i][1];
   gdofs[4+2*i]=w_at_xi[i][2];

   // Finally fill in the second derivatives
    gdofs[ 9+3*i]=w_at_xi[i][3];
    gdofs[10+3*i]=w_at_xi[i][5];
    gdofs[11+3*i]=w_at_xi[i][4];
   // Internal dofs
   gdofs[18+i]=w_at_xi[3+i][0];
   }
  return gdofs;
 }

// Return the global dofs of a function on an element
Vector<double> get_basis_global_dofs(const ExactSolnFctPt& w)
 {
  //Initialise
  Vector<Vector<double> > ai(3,Vector<double>(2,0.0)),
                          ei(3,Vector<double>(2,0.25)),
                          w_at_xi(6,Vector<double> (6,0.0));
  ai=get_vertices();
  Vector<double> gdofs(21,0.0);

  // Fill in ei
  ei[0][0]=0.5;
  ei[1][1]=0.5;

  // Get the dofs
  for (unsigned i=0; i<3; ++i)
   {
    // Fill in Vectors on nodes
    (*w)(ai[i],w_at_xi[i]);

    // Fill in Vectors at internal points
    Vector<double> x(2); f_k(ei[i],x);
    (*w)(x,w_at_xi[3+i]);
   }

  // Now rearrange them
  for (unsigned i=0; i< 3; ++i)
   {
   // First three dofs are on nodes
   gdofs[6*i+0]=w_at_xi[i][0];
   gdofs[6*i+1]=w_at_xi[i][1];
   gdofs[6*i+2]=w_at_xi[i][2];
   gdofs[6*i+3]=w_at_xi[i][3];
   gdofs[6*i+5]=w_at_xi[i][4];
   gdofs[6*i+4]=w_at_xi[i][5];
   // The get dofs are the wrong way around

   // Internal dofs
   gdofs[18+i]=w_at_xi[3+i][0];
   }
  return gdofs;
 }
// These next two won't be pretty.

// Return the global dofs of a function on an element
void check_basic_shape(const double& tol)
 {
  //Initialise
  Vector<double> bdofs(36,0.0);
  Vector<Vector<double> > ai(3,Vector<double>(2,0.0)),
    ei(3,Vector<double>(2,0.25)), ni(3,Vector<double>(2,0.0)),
    bi(3,Vector<double>(2,0.5)),di(6,Vector<double>(2,0.0));

  // fill in ai (local)
  ai[0][0]=1.0;  ai[1][1]=1.0;

  // Fill in ei
  ei[0][0]=0.5;  ei[1][1]=0.5;

  // Fill in bi
  bi[0][0]=0.0;  bi[1][1]=0.0;

  // Fill in the normal vectors
  ni[0][0]=-1.0;  ni[0][1]= 0.0;
  ni[1][0]= 0.0;  ni[1][1]=-1.0;
  ni[2][0]= std::sqrt(2.)/2.;
  ni[2][1]= std::sqrt(2.)/2.;

  // Fill in di
  di[0][0]=0.0;   di[0][1]=0.75;
  di[1][0]=0.0;   di[1][1]=0.25;
  di[2][0]=0.25;  di[2][1]=0.0;
  di[3][0]=0.75;  di[3][1]=0.0;
  di[4][0]=0.75;  di[4][1]=0.25;
  di[5][0]=0.25;  di[5][1]=0.75;

  // Get the shape functions
  Shape p7 (36), m7 (36);
  DShape dp7 (36,2), dm7 (36,2), d2p7 (36,3), d2m7 (36,3);
  DenseMatrix<double> a_matrix (36,36,0.0);
  monomial_to_basic_matrix(a_matrix);

  // Shape function dofs
  DenseMatrix<double> shape_dofs(36,36,0.0);

  // Do the ai
  for(unsigned i=0; i<3;++i)
   {
    // Construct shape functions
    full_basis_monomials(ai[i],p7);

    for (unsigned j=0; j<36; ++j)
     {
      m7(j)=0.0;
      for (unsigned k=0; k<36; ++k)
       {
        // Fill in basis at ai
        m7(j)+=a_matrix(j,k)*p7(k);
       }
      // Fill in matrix entries
      shape_dofs(i,j)=m7(j);
     }
   }

  // Do the ai
  for(unsigned i=0; i<3;++i)
   {
    // Construct shape functions
    dfull_basis_monomials(ai[i],dp7);
 //   for(unsigned k=0;k<36;++k)
 //       std::cout<<dp7(k,0)<<" "<<dp7(k,1)<<"\n";

    for (unsigned j=0; j<36; ++j)
     {
      dm7(j,0)=0.0;
      dm7(j,1)=0.0;
      for (unsigned k=0; k<36; ++k)
       {
        // Fill in basis at ai
        dm7(j,0)+=a_matrix(j,k)*dp7(k,0);
        dm7(j,1)+=a_matrix(j,k)*dp7(k,1);
       }
     }

    // Fill in matrix entries
    for (unsigned j=0; j<36; ++j)
      shape_dofs(3+2*i,j)=dm7(j,0);
    for (unsigned j=0; j<36; ++j)
      shape_dofs(4+2*i,j)=dm7(j,1);
   }

  // Do the ai
  for(unsigned i=0; i<3;++i)
   {
    // Construct shape functions
    d2full_basis_monomials(ai[i],d2p7);
  //  for(unsigned k=0;k<36;++k)
  //      std::cout<<d2p7(k,0)<<" "<<d2p7(k,1)<<" "<<d2p7(k,2)<<"\n";

    for (unsigned j=0; j<36; ++j)
     {
      d2m7(j,0)=0.0;
      d2m7(j,1)=0.0;
      d2m7(j,2)=0.0;
      for (unsigned k=0; k<36; ++k)
       {
        // Fill in basis at ai
        d2m7(j,0)+=a_matrix(j,k)*d2p7(k,0);
        d2m7(j,1)+=a_matrix(j,k)*d2p7(k,1);
        d2m7(j,2)+=a_matrix(j,k)*d2p7(k,2);
       }
      // Fill in matrix entries
      shape_dofs(9+3*i,j)=d2m7(j,0);
      shape_dofs(10+3*i,j)=d2m7(j,1);
      shape_dofs(11+3*i,j)=d2m7(j,2);
     }
   }

  // Do the bi
  for(unsigned i=0; i<3;++i)
   {
    // Construct shape functions
    dfull_basis_monomials(bi[i],dp7);

    for (unsigned j=0; j<36; ++j)
     {
      dm7(j,0)=0.0;
      dm7(j,1)=0.0;
      for (unsigned k=0; k<36; ++k)
       {
        // Fill in basis at ai
        dm7(j,0)+=a_matrix(j,k)*dp7(k,0);
        dm7(j,1)+=a_matrix(j,k)*dp7(k,1);
       }
      // Fill in matrix entries
      shape_dofs(18+i,j)=dm7(j,0)*ni[i][0]+dm7(j,1)*ni[i][1];
     }
   }

  // Do the dis
  for(unsigned i=0; i<6;++i)
   {
    // Construct shape functions
    full_basis_monomials(di[i],p7);

    for (unsigned j=0; j<36; ++j)
     {
      m7(j)=0.0;
      for (unsigned k=0; k<36; ++k)
       {
        // Fill in basis at ai
        m7(j)+=a_matrix(j,k)*p7(k);
       }
      // Fill in matrix entries
      shape_dofs(21+i,j)=m7(j);
     }
   }

  // Do the bi
  for(unsigned i=0; i<6;++i)
   {
    // Construct shape functions
    dfull_basis_monomials(di[i],dp7);

    for (unsigned j=0; j<36; ++j)
     {
      dm7(j,0)=0.0;
      dm7(j,1)=0.0;
      for (unsigned k=0; k<36; ++k)
       {
        // Fill in basis at ai
        dm7(j,0)+=a_matrix(j,k)*dp7(k,0);
        dm7(j,1)+=a_matrix(j,k)*dp7(k,1);
       }
      // Fill in matrix entries
      shape_dofs(27+i,j)=dm7(j,0)*ni[i/2][0]+dm7(j,1)*ni[i/2][1];
     }
   }

  // Do the eis
  for(unsigned i=0; i<3;++i)
   {
    // Construct shape functions
    full_basis_monomials(ei[i],p7);

    for (unsigned j=0; j<36; ++j)
     {
      m7(j)=0.0;
      for (unsigned k=0; k<36; ++k)
       {
        // Fill in basis at ai
        m7(j)+=a_matrix(j,k)*p7(k);
       }
    // Fill in matrix entries
    shape_dofs(33+i,j)=m7(j);
    }
   }

  //Output any non zeros
  std::cout<<"\n";
  for(unsigned i=0;i<21;++i)
    for(unsigned j=0;j<21;++j)
      std::cout<<(fabs(shape_dofs(i,j))>1e-11? shape_dofs(i,j): 0.0 )<<(j!=20 ? " ":"\n");

  //Output any non zeros
  for(unsigned i=0;i<36;++i)
    for(unsigned j=0;j<36;++j)
      if(std::abs(shape_dofs(i,j)-((i==j)?1.0:0.0))>tol)
        std::cout/*<<std::scientific*/
                 <<"Nonzero difference at ("<<i<<","<<j<<"): "
                 <<shape_dofs(i,j)-((i==j)?1.0:0.0)<<"\n";

 }


// degree 5 analytic function
void p5_1d(const double& x, Vector<double>& p)
 {
  // Clear vector
  p=Vector<double>(3,0.0);
  // Initialise constants
  Vector<double> a(6,0.0);
  a[0]=0.1,a[1]=-0.1,a[2]=0.375,a[3]=0.6,a[4]=1.0,a[5]=-0.2;
  for(unsigned i=0; i<6;++i)
   {
    // Polynomial
    p[0] += a[i]*pow(x,i);
    // Derivatives
    if(i>0)
     p[1] += i*a[i]*pow(x,i-1);
    if(i>1)
     p[2] += i*(i-1)*a[i]*pow(x,i-2);
   }
 }

// degree 3 analytic function
void p3_1d(const double& x, Vector<double>& p)
 {
  // Clear vector
  p=Vector<double>(3,0.0);
  // Initialise constants
  Vector<double> a(6,0.0);
  a[0]=0.3,a[1]=-0.5,a[2]=0.675,a[3]=1.0;
  for(unsigned i=0; i<4;++i)
   {
    // Polynomial
    p[0] += a[i]*pow(x,i);
    // Derivatives
    if(i>0)
     p[1] += i*a[i]*pow(x,i-1);
    if(i>1)
     p[2] += i*(i-1)*a[i]*pow(x,i-2);
   }
 }

// Return the global dofs of a function on an element
void check_1d_hermite_shape(const double& tol)
 {
  // Get the shape functions
  Shape p6(6);
  Vector<double> vertices(2,0.0);
  vertices[0]=0.0;  vertices[1]=1.0;

  // Now loop over nodes
  DenseMatrix<double> shape_matrix_5(6,6,0.0);
  DenseMatrix<double> d_shape_matrix_5(6,6,0.0);
  DenseMatrix<double> shape_matrix_3(4,4,0.0);
  for(unsigned i=0;i<2;++i)
   {
    // Initialise
    Shape psi5(6),psi3(4);
    DShape dpsi5(6,1),dpsi3(4,1);

    //Fill in
    hermite_shape_1d_5(vertices[i],psi5);
    hermite_shape_1d_3(vertices[i],psi3);

    d_hermite_shape_1d_5(vertices[i],dpsi5);

   for(unsigned j=0;j<6;++j)
      {
       // Fill in shape matrix row 1 and 4
       shape_matrix_5(i,j) = psi5[j];
       d_shape_matrix_5(i+2,j) = dpsi5(j,0);

       // Initialise the stencil
       double h=0.01, s=vertices[i];
       Vector<double> fd_points_5(7,0.0);
       Vector<double> d_fd_points_5(7,0.0);

       // Loop finite difference points
       for(unsigned l=0; l<7;++l)
        {
         int p=l-3;
         Shape chi(6);
         DShape dchi(6,1);

         // Fill in fd values
         hermite_shape_1d_5(s+p*h,chi);
         fd_points_5[l]=chi[j];
         d_hermite_shape_1d_5(s+p*h,dchi);
         d_fd_points_5[l]=dchi(j,0);
        }

      // Fill in matrix entries
      shape_matrix_5(i+2,j)=dfdx_fd6(fd_points_5,h);
      shape_matrix_5(i+4,j)=d2fdx2_fd6(fd_points_5,h);
      d_shape_matrix_5(i+4,j)=dfdx_fd6(d_fd_points_5,h);
      }

    for(unsigned j=0;j<4;++j)
      {
       shape_matrix_3(i,j) = psi3[j];

       // Initialise the stencil
       double h=0.01, s=vertices[i];
       Vector<double> fd_points_3(7,0.0);

       // Loop finite difference points
       for(unsigned l=0; l<7;++l)
        {
         int p=l-3;
         Shape chi(4);
         // Fill in fd values
         hermite_shape_1d_3(s+p*h,chi);
         fd_points_3[l]=chi[j];
        }

      // Fill in matrix entries
      shape_matrix_3(i+2,j)=dfdx_fd6(fd_points_3,h);
      }
   }

  // Output
  std::cout<<"\nshape_matrix_5 should be identity matrix to machine precision:\n";
  for(unsigned i=0;i<6;i++)
   for(unsigned j=0;j<6;j++)
     if(std::abs(shape_matrix_5(i,j)-(i==j? 1.0 :0.0))>tol)
      std::cout<<"shape matrix - I  nonzero at ("<<i<<","<<j<<"): "
               <<shape_matrix_5(i,j)-(i==j? 1.0 :0.0)<<"\n";

  std::cout<<"\nd_shape_matrix_5 should be I3 in bottom right to machine"
           <<" precision:\n";
  for(unsigned i=0;i<6;i++)
    for(unsigned j=0;j<6;j++)
     if(std::abs(d_shape_matrix_5(i,j)-(i==j && j>1 ? 1.0 :0.0))>tol)
      std::cout<<"d_shape matrix - I  nonzero at ("<<i<<","<<j<<"): "
               <<d_shape_matrix_5(i,j)-(i==j && j>1 ? 1.0 :0.0)<<"\n";

  std::cout<<"\nshape_matrix_3 should be identity matrix to machine precision:\n";
  for(unsigned i=0;i<4;i++)
   for(unsigned j=0;j<4;j++)
     if(std::abs(shape_matrix_3(i,j)-(i==j? 1.0 :0.0))>tol)
      std::cout<<"shape matrix - I  nonzero at ("<<i<<","<<j<<"): "
               <<shape_matrix_3(i,j)-(i==j? 1.0 :0.0)<<"\n";

  // Now check we can represent a generic degree 5 polynomial
  Vector<double> p5_at_zero(3,0.0), p5_at_one(3,1.0), p3_at_zero(3,0.0),
    p3_at_one(3,1.0);

  // Fill in the vectors
  p5_1d(0,p5_at_zero);
  p5_1d(1,p5_at_one);
  p3_1d(0,p3_at_zero);
  p3_1d(1,p3_at_one);

  // Fill in the dofs
  Vector<double> p5_dofs(6,0.0);
  p5_dofs[0]=p5_at_zero[0];
  p5_dofs[1]=p5_at_one[0];
  p5_dofs[2]=p5_at_zero[1];
  p5_dofs[3]=p5_at_one[1];
  p5_dofs[4]=p5_at_zero[2];
  p5_dofs[5]=p5_at_one[2];

  Vector<double> p3_dofs(4,0.0);
  p3_dofs[0]=p3_at_zero[0];
  p3_dofs[1]=p3_at_one[0];
  p3_dofs[2]=p3_at_zero[1];
  p3_dofs[3]=p3_at_one[1];

  unsigned n_points=51;
  for(unsigned i=0; i<n_points; ++i)
   {
    // The 'nodes'
    double s=i*1.0/(n_points-1);

    //Initialise
    double p5_apx=0.0,p3_apx=0.0;
    Shape psi_5(6),psi_3(4);
    hermite_shape_1d_5(s,psi_5);
    hermite_shape_1d_3(s,psi_3);

    // Get exact
    Vector<double> p5_ex(3);
    Vector<double> p3_ex(3);
    p5_1d(s,p5_ex);
    p3_1d(s,p3_ex);

    // Get aprox
    for(unsigned j=0; j< 6; ++j)
     {
      p5_apx+=p5_dofs[j]*psi_5[j];
      if(j<4)
       p3_apx+=p3_dofs[j]*psi_3[j];
     }
    if(p5_apx-p5_ex[0]>tol)
     std::cout<<"\nNonzero diff between analytic and aprox p5:\n"
             <<"s \t Aprox \t True \t Diff\n"
             <<s<<"\t"<<p5_apx<<"\t"<<p5_ex[0]<<"\t"<<p5_apx-p5_ex[0]<<"\n";

    if(p3_apx-p3_ex[0]>tol)
     std::cout<<"\nNonzero diff between analytic and aprox p3:\n"
             <<"s \t Aprox \t True \t Diff\n"
             <<s<<"\t"<<p3_apx<<"\t"<<p3_ex[0]<<"\t"<<p3_apx-p3_ex[0]<<"\n";
   }
 }

// Check the submatrix B2 using the Jacobian and Hessian
void check_b2l_submatrix_3(const double& tol)
{
// Get the B3 matrix
DenseMatrix<double> B3 (21,9,0.0);
basic_to_local_submatrix_3(B3);

// Output
//std::cout<<"B3:\n";
//for(unsigned i=0;i<21;++i)
// for(unsigned j=0;j<9;++j)
//   std::cout<<B3(i,j)<<(j==8?"\n":" ");

// Get the d matrix
DenseMatrix<double> d (21,21,0.0);
local_to_global_matrix(d);

// Get the vertices
Vector<Vector<double> > vertices(3,Vector<double>(2,0.0));
vertices=get_vertices();


Vector<Vector<double> > local_vertices(3,Vector<double>(2,0.0));
local_vertices[0][0]=1.0; local_vertices[1][1]=1.0;

// Now compare the parts at the different nodes
for(unsigned n=0;n<3;++n)
 {
  // Get Jacobian
  DenseMatrix<double> jacobian(2,2,0.0);
  get_basic_jacobian(local_vertices[n],jacobian);
  // We need these for comparison
  DenseMatrix<double> Datnode(5,5,0.0);
  DenseMatrix<double> Batnode(5,3,0.0);

// We have checked Jacobian and Hessian
//  std::cout<<"Jacobian: \n";
//  std::cout<<jacobian(0,0)<<" "<<jacobian(0,1)<<"\n"
//           <<jacobian(1,0)<<" "<<jacobian(1,1)<<"\n";

  // Get Hessian
  RankThreeTensor<double> hessian(2,2,2);
  get_basic_hessian(local_vertices[n],hessian);

//  std::cout<<"Hessian:\n";
//  std::cout<<hessian(0,0,0)<<" "<<hessian(0,0,1)<<" "<<hessian(0,1,1)<<"\n"
//           <<hessian(1,0,0)<<" "<<hessian(1,0,1)<<" "<<hessian(1,1,1)<<"\n";

  // Construct M (5 x 3)
  DenseMatrix<double> M(5,3,0.0);
  for(unsigned i=0; i<2; ++i)
    for(unsigned j=0; j<2; ++j)
      for(unsigned k=j; k<2; ++k)
        M(i,j+k)+=hessian(i,j,k);

  // Construct M (5 x 3)
  for(unsigned i=0; i<2; ++i)
    for(unsigned j=0; j<2; ++j)
      for(unsigned k=0; k<2; ++k)
        for(unsigned l=j; l<2; ++l)
          M(2+i+k,j+l)+=jacobian(i,j)*jacobian(k,l);

  // Fill in Batnode
  // b31
  for(unsigned i=0; i<2; ++i)
    for(unsigned j=0; j<3; ++j)
      Batnode(i,j)=B3(3 + 2*n+i ,3*n+j);

  // b32 pt 1
  for(unsigned i=0; i<2; ++i)
    for(unsigned j=0; j<3; ++j)
      Batnode(2+i,j)=B3(9 + 2*n+i ,3*n+j);

  // b32 pt 2
  for(unsigned j=0; j<3; ++j)
      Batnode(4,j)=B3(15 + n ,3*n+j);

  // Fill in Datnode
  for(unsigned i=0; i<2; ++i)
    for(unsigned j=0; j<2; ++j)
      Datnode(i,j)=d(3 + 2*n+i ,3 + 2*n+j);

  for(unsigned i=0; i<3; ++i)
    for(unsigned j=0; j<2; ++j)
      Datnode(2+i,2+j)=d(9 + 3*n+i ,9 + 2*n+j);

  for(unsigned j=0; j<3; ++j)
      Datnode(2+j,4)=d(9+3*n+j,15 + n );

//  //Output
//  std::cout<<"Datnode:\n";
//  for(unsigned i=0; i<5; ++i)
//    for(unsigned j=0; j<5; ++j)
//      std::cout<<Datnode(i,j)<<(j==4?"\n":" ");;
//
//  //Output
//  std::cout<<"Batnode:\n";
//  for(unsigned i=0; i<5; ++i)
//    for(unsigned j=0; j<3; ++j)
//      std::cout<<Batnode(i,j)<<(j==2?"\n":" ");;
//
//  // Output M - we have checked this matrix
//  std::cout<<"M"<<n<<": \n";
//  for(unsigned i=0;i<5;++i)
//   for(unsigned j=0;j<3;++j)
//     std::cout<<M(i,j)<<(j==2?"\n":" ");


  // Now perform the est!
  for(unsigned i=0;i<5;++i)
   {
    for(unsigned j=0;j<3;++j)
     {
      // Residual for the check
      double residual = M(i,j);
      for(unsigned k=0;k<5;++k)
       {
        residual-=Datnode(i,k)*Batnode(k,j);
       }
      //Output the resulting matrix
      if(std::abs(residual)>tol)
       std::cout<<"Non zero result for test at node "<< n
                <<" Dik Bkj element ("<<i<<","<<j<<")"
                <<", with value : "<<residual <<"\n";
     }
   }
 }
}


// Check the Traces against analytic function
void check_traces(const ExactSolnFctPt&
get_analyticfunction, const double& tol, bool ignore_side_3, bool ignore_g )
{
 // Get the dofs
 Vector<double> dofs (21,0.0);
 dofs=get_global_dofs(get_analyticfunction);

 // Now get D matrix
 DenseMatrix<double> d (21,21,0.0);
 local_to_global_matrix(d);

 // Get the local dofs
 Vector<double> local_dofs(21,0.0);
 for(unsigned i=0; i<21;++i)
  {
   for(unsigned j=0; j<21;++j)
    {
     local_dofs[j]+=dofs[i]*d(i,j);
    }
  }

 // Now get the traces at several values of s
 unsigned n_points=21;
 Vector<double> values_f1(n_points,0.0), values_f2(n_points,0.0),
                values_f3(n_points,0.0), exact_values_f3(n_points,0.0),
                exact_values_f1(n_points,0.0), exact_values_f2(n_points,0.0);

 Vector<double> values_g1(n_points,0.0), values_g2(n_points,0.0),
                values_g3(n_points,0.0), exact_values_g3(n_points,0.0),
                exact_values_g1(n_points,0.0), exact_values_g2(n_points,0.0);

 for(unsigned i=0; i<n_points;++i)
  {
   // Initialise position vector
   Vector<double> s(2,0.0);
   double t=1.0*i*1./(n_points-1);

   // Test on side 1
   s[0]=t;
   s[1]=0.0;

   // Sum up at this point
   Vector<double> f1(21,0.0),g1(21,0.0);
   f1=f_1(t);
   g1=g_1(t);

   // Loop
   for (unsigned j=0; j<21; ++j)
    {
     values_f1[i]+=local_dofs[j]*f1[j];
     values_g1[i]+=local_dofs[j]*g1[j];
    }

   // Get the exact result
   Vector<double> pex(6,0.0),x(2);
   f_k(s,x);  
   (*get_analyticfunction)(x,pex);
   exact_values_f1[i]=pex[0];

  // Report any differences
   if(std::abs(values_f1[i] - exact_values_f1[i])> tol)
     std::cout<<"Nonzero difference for trace on side 2 (f_1):\n"
              <<values_f1[i] - exact_values_f1[i]<<" at "<< s<<"\n";

   // Do the g bit
   if(!ignore_g)
    {
     Vector<double> pex(6,0.0),x(2);
     f_k(s,x);  
     (*get_analyticfunction)(x,pex);
     for(unsigned j=0; j<2;++j)
      exact_values_g1[i]+=pex[1+j]*altitude_vector_2(j);
     // Check the difference
     if(std::abs(values_g1[i] - exact_values_g1[i])> tol)
       std::cout<<"Nonzero difference for trace on side 2 (g_1):\n"
                <<values_g1[i]<< " "<<exact_values_g1[i]<<" "
                <<values_g1[i] - exact_values_g1[i]<<" at "<< s<<"\n";
    }

   // Test on side 2
   s[0]=0.0;
   s[1]=t;

  // Sum up at this point
   Vector<double> f2(21,0.0), g2(21,0.0);
   f2=f_2(t);
   g2=g_2(t);
   // Loop
   for (unsigned j=0; j<21; ++j)
    {
     values_f2[i]+=local_dofs[j]*f2[j];
     values_g2[i]+=local_dofs[j]*g2[j];
    }

   // Get the exact result
   pex=Vector<double>(6,0.0);
   f_k(s,x);  
   (*get_analyticfunction)(x,pex);
   exact_values_f2[i]=pex[0];

   // Report differences
   if(std::abs(values_f2[i] - exact_values_f2[i])> tol)
     std::cout<<"Nonzero difference for trace on side 1 (f_2):\n"
              <<values_f2[i] - exact_values_f2[i]<<" at "<< s<<"\n";

   if(!ignore_g)
    {
     Vector<double> pex(6,0.0);
     f_k(s,x);  
     (*get_analyticfunction)(x,pex);
     for(unsigned j=0; j<2;++j)
      exact_values_g2[i]+=pex[1+j]*altitude_vector_1(j);
     // Check the difference
     if(fabs(values_g2[i] - exact_values_g2[i])> tol)
       std::cout<<"Nonzero difference for trace on side 1 (g_2):\n"
                <<values_g2[i]<< " "<<exact_values_g2[i]<<" "
                <<values_g2[i] - exact_values_g2[i]<<" at "<< s<<"\n";
    }

   // Test on side 3
   if(!ignore_side_3)
   {
   s[0]=t;
   s[1]=1-t;

   Vector<double> f3(21,0.0);
   f3=f_3(t);
   // Loop
   for (unsigned j=0; j<21; ++j)
    for (unsigned k=0; k<21; ++k)
    values_f3[i]+=dofs[j]*d(j,k)*f3[k];

   // Get the exact result
   pex=Vector<double>(6,0.0);
   f_k(s,x);
   (*get_analyticfunction)(x,pex);
   exact_values_f3[i]=pex[0];

   if(std::abs(values_f3[i] - exact_values_f3[i])> tol)
     std::cout<<"Nonzero difference for trace on side 1 (f_3):\n"
                <<values_f3[i]<< " "<<exact_values_f3[i]<<" "
                <<values_f3[i] - exact_values_f3[i]<<" at "<< s<<"\n";
   }
  }
 // Does it agree with our analytic function
//  std::cout<<exact_values_g2<<"\n";
//  std::cout<<values_g2<<"\n";

}

// Check the Traces against what we would expect
void check_g3_trace(const double& tol)
{
 // Now we check G3 using the other submatrices
 // Get the B3 matrix
 DenseMatrix<double> B2 (21,6,0.0);
 basic_to_local_submatrix_2(B2);

 DenseMatrix<double> B3 (21,9,0.0);
 basic_to_local_submatrix_3(B3);

 const unsigned n_points=11;
 // Loop over some points
 for(unsigned i=0; i<n_points;++i)
  {
   // Initialise position vector
   double t=1.0*i*1./(n_points-1);
   // Now initialise the B3 matrix
   Vector<double> g3(21,0.0);
   g3=g_3(t);

   // Get 1d shape
   Shape psi(4);
   hermite_shape_1d_3(t,psi);

   // Now initialise g3
   Vector<double> g3_exact(21,0.0);
   for(unsigned j=0;j<21;++j)
    {
     // Normal derivative at node 0
     g3_exact[j]+=(-0.5*B2(j,0)-0.5*B2(j,1))*psi[1];
     // Normal derivative at node 1
     g3_exact[j]+=(-0.5*B2(j,2)-0.5*B2(j,3))*psi[0];

     // D2 w (a0)(n,t) derivative at node 0
     g3_exact[j]+=(-0.5*B3(j,0)+0.5*B3(j,2))*psi[3];
     // D2 w (a1)(n,t) derivative at node 1
     g3_exact[j]+=(-0.5*B3(j,3)+0.5*B3(j,5))*psi[2];
    }

  // Now compare exact with aprox.
  for(unsigned j=0;j<21;++j)
   {
    if(std::abs(g3_exact[j]-g3[j])>0)
      std::cout<<"Non zero difference for g3_exact at dof: "<<i
               << " diff: "<<g3_exact[j]-g3[j]<<"\n";
   }
  }
}

// Return the global dofs of a function on an element
void check_shape(const double& tol)
 {
  //Initialise
  Vector<Vector<double> > ei(3,Vector<double>(2,0.25)),
                          ai(3,Vector<double>(2,0.0));
  Vector<double> bdofs(36,0.0);

  // fill in ai (local)
  ai[0][0]=1.0; ai[1][1]=1.0;

  // Fill in ei
  ei[0][0]=0.5 ; ei[1][1]=0.5;

  // Get the shape functions
  Shape p7 (36), m7 (36);
  DShape dp7 (36,2), dm7 (36,2);
  DShape d2p7 (36,3), d2m7 (36,3);

  // Get the matrices
  DenseMatrix<double> a_matrix (36,36,0.0), b_matrix (21,36,0.0),
   d_matrix (21,21,0.0),  conversion_matrix (21,36,0.0),
   gl2basic_matrix (21,36,0.0);

  monomial_to_basic_matrix(a_matrix);
  basic_to_local_matrix(b_matrix);
  local_to_global_matrix(d_matrix);

  // Fill in conversion matrix
  for(unsigned i=0;i<21;++i)
   for(unsigned j=0;j<36;++j)
     for(unsigned k=0;k<21;++k)
      for(unsigned l=0;l<36;++l)
        conversion_matrix(i,j)+=d_matrix(i,k)*b_matrix(k,l)*a_matrix(l,j);

//  // Fill in conversion matrix
//  for(unsigned i=0;i<21;++i)
//    for(unsigned k=0;k<21;++k)
//     for(unsigned l=0;l<36;++l)
//       gl2basic_matrix(i,l)+=d_matrix(i,k)*b_matrix(k,l);
//
//   for(unsigned i=0;i<21;++i)
//     for(unsigned l=0;l<36;++l)
//     std::cout<<gl2basic_matrix(i,l)<<(l==35 ? "\n":" ");

  // Shape function dofs
  DenseMatrix<double> shape_dofs(21,21,0.0);

  // Do the ai
  for(unsigned i=0; i<3;++i)
   {
    // Construct shape functions
    full_basis_monomials(ai[i],p7);

    for (unsigned j=0; j<21; ++j)
     {
      m7(j)=0.0;
      for (unsigned k=0; k<36; ++k)
       {
        // Fill in basis at ai
        m7(j)+=conversion_matrix(j,k)*p7(k);
       }
      // Fill in matrix entries
      shape_dofs(i,j)=m7(j);
     }
   }

  // Do the ai
  for(unsigned i=0; i<3;++i)
   {
    // Construct shape functions
    dfull_basis_monomials(ai[i],dp7);

    for (unsigned j=0; j<21; ++j)
     {
      // Initialise
      Vector<double> dm7j_ds(2,0.0);
      dm7(j,0)=0.0;
      dm7(j,1)=0.0;
      for (unsigned k=0; k<36; ++k)
       {
        // Fill in basis at ai
        dm7j_ds[0]+=conversion_matrix(j,k)*dp7(k,0);
        dm7j_ds[1]+=conversion_matrix(j,k)*dp7(k,1);
       }

      // Now transform
      DenseMatrix<double> jacobian(2,2,0.0), inv_jacobian(2,2,0.0);
      get_basic_jacobian(ai[i],jacobian);
      // Now invert
      inv_jacobian=invert_two_by_two(jacobian);

      // Transform
      for(unsigned l=0; l<2;++l)
       {
        for(unsigned k=0; k<2;++k)
        {
         dm7(j,k)+=inv_jacobian(l,k)*dm7j_ds[l];
        }
       }
    }
    // Fill in matrix entries
    for (unsigned j=0; j<21; ++j)
      shape_dofs(3+2*i,j)=dm7(j,0);
    for (unsigned j=0; j<21; ++j)
      shape_dofs(4+2*i,j)=dm7(j,1);
   }

  // Do the ai
  for(unsigned i=0; i<3;++i)
   {
    // Construct shape functions
    d2full_basis_monomials(ai[i],d2p7);
    dfull_basis_monomials(ai[i],dp7);
 //   std::cout<<ai[i]<<"\n";
   // for(unsigned k=0;k<36;++k)
   //     std::cout<<d2p7(k,0)<<" "<<d2p7(k,1)<<" "<<d2p7(k,2)<<"\n";
    for (unsigned j=0; j<21; ++j)
     {
      //Initialise
      DenseMatrix<double> d2m7j_ds2(2,2,0.0);
      Vector<double> dm7j_ds(2,0.0);
      d2m7(j,0)=0.0;
      d2m7(j,1)=0.0;
      d2m7(j,2)=0.0;
      for (unsigned k=0; k<36; ++k)
       {
        // Fill in basis at ai
        dm7j_ds[0]+=conversion_matrix(j,k)*dp7(k,0);
        dm7j_ds[1]+=conversion_matrix(j,k)*dp7(k,1);
        d2m7j_ds2(0,0)+=conversion_matrix(j,k)*d2p7(k,0);
        d2m7j_ds2(1,0)+=conversion_matrix(j,k)*d2p7(k,1);
        d2m7j_ds2(0,1)+=conversion_matrix(j,k)*d2p7(k,1);
        d2m7j_ds2(1,1)+=conversion_matrix(j,k)*d2p7(k,2);
       }

      // Now transform
      DenseMatrix<double> jacobian(2,2,0.0), inv_jacobian(2,2,0.0);
      RankThreeTensor<double> hessian(2,2,2,0.0), d_inv_jac_ds(2,2,2,0.0);
      get_basic_jacobian(ai[i],jacobian);
      get_basic_hessian(ai[i],hessian);
//      for(unsigned alpha=0;alpha<2;++alpha)
//       for(unsigned beta =0;beta<2 ;++beta )
//         std::cout<<"("<<hessian(alpha,beta,0)<<","<<hessian(alpha,beta,1)<<")"
//                  <<(beta==1 ? "\n":" ");
      // Now invert
      inv_jacobian=invert_two_by_two(jacobian);

      //Fill in values
       for(unsigned alpha=0;alpha<2;++alpha)
        for(unsigned beta =0;beta<2 ;++beta )
         for(unsigned gamma=0;gamma<2;++gamma)
          for(unsigned delta=0;delta<2;++delta)
            for(unsigned zeta=0 ;zeta<2; ++zeta)
             d_inv_jac_ds(alpha,delta,zeta)-=inv_jacobian(alpha,beta)
                  *hessian(beta,gamma,zeta)*inv_jacobian(gamma,delta);

//      for(unsigned alpha=0;alpha<2;++alpha)
//       for(unsigned beta =0;beta<2 ;++beta )
//         std::cout<<"("<<d_inv_jac_ds(alpha,beta,0)<<","<<d_inv_jac_ds(alpha,beta,1)<<")"
//                  <<(beta==1 ? "\n":" ");

//       std::cout<<"{{"<<d2m7j_ds2(0,0)<<","<<d2m7j_ds2(0,1)<<"},"<<
//                   "{"<<d2m7j_ds2(1,0)<<","<<d2m7j_ds2(1,1)<<"}}\n";
       // Now find the global derivatives at local coordinates
       for(unsigned alpha=0;alpha<2;++alpha)
        for(unsigned beta =alpha;beta<2 ;++beta )
         for(unsigned gamma=0;gamma<2;++gamma)
          for(unsigned delta=0;delta<2;++delta)
            d2m7(j,alpha + beta)+=d2m7j_ds2(gamma,delta)*inv_jacobian(gamma,alpha)
                *inv_jacobian(delta,beta)
             + dm7j_ds[gamma]*d_inv_jac_ds(gamma,beta,delta)
               *inv_jacobian(delta,alpha);
      // Fill in matrix entries
      shape_dofs(9+3*i,j) =d2m7(j,0);
      shape_dofs(10+3*i,j)=d2m7(j,1);
      shape_dofs(11+3*i,j)=d2m7(j,2);
     }
   }

//  std::cout<<"The m7: \n[";
//  for(unsigned i=0;i<36;++i)
//      std::cout<<p7[i]<<(i==36?"]":",");
//  std::cout<<"\n";

  // Do the eis
  for(unsigned i=0; i<3;++i)
   {
    // Construct shape functions
    full_basis_monomials(ei[i],p7);

    for (unsigned j=0; j<21; ++j)
     {
      m7(j)=0.0;
      for (unsigned k=0; k<36; ++k)
       {
        // Fill in basis at ai
        m7(j)+=conversion_matrix(j,k)*p7(k);
       }
    // Fill in matrix entries
    shape_dofs(18+i,j)=m7(j);
    }
   }

  //Output any non zeros
  for(unsigned i=0;i<21;++i)
    for(unsigned j=0;j<21;++j)
      std::cout<<(shape_dofs(i,j)>tol? shape_dofs(i,j): 0.0 )<<(j!=20 ? " ":"\n");

  for(unsigned i=0;i<21;++i)
    for(unsigned j=0;j<21;++j)
      if(std::abs(shape_dofs(i,j)-((i==j)?1.0:0.0))>tol)
        std::cout/*<<std::scientific*/
                 <<"Nonzero difference at ("<<i<<","<<j<<"): "
                 <<shape_dofs(i,j)-((i==j)?1.0:0.0)<<"\n";

 }

// Return the global dofs of a function on an element
void check_get_shape(const double& tol)
 {
  // Initialise shape
  Shape psi(3,6);
  DShape dpsi(3,6,2);
  DShape d2psi(3,6,3);

  // Initialise bubble shape
  Shape bpsi(3,1);
  DShape dbpsi(3,1,2);
  DShape d2bpsi(3,1,3);

  // Initialise local vertices
  Vector<Vector<double> > lvertices(3,Vector<double>(2,0.0));
  lvertices[0][0]=1.0;
  lvertices[1][1]=1.0;

  Vector<Vector<double> > linternpts(3,Vector<double>(2,0.25));
  linternpts[0][0]=0.5;
  linternpts[1][1]=0.5;

  DenseMatrix<double> shape_matrix(21,21,0.0);

  // Now check the shape
  // Loop over the shape functions
  for(unsigned i=0;i<21;++i)
   {
    // The nodal basis functions
    if(i<18)
     {
      // Now we have the indices
      unsigned n = i / 6;
      unsigned l = i % 6;
      // basis at nth node
      // Loop over the dofs
      for(unsigned ii=0;ii<3;++ii)
       {
        // At the nodes the Hermite dofs
        unsigned nn=ii;
        d2_shape_dx2(lvertices[nn],psi,bpsi,dpsi,dbpsi,d2psi,d2bpsi);
        d_shape_dx(lvertices[nn],psi,bpsi,dpsi,dbpsi);
        shape(lvertices[nn],psi,bpsi);

        shape_matrix(i,6*nn)  =  psi(n,l);
        shape_matrix(i,6*nn+1)= dpsi(n,l,0);
        shape_matrix(i,6*nn+2)= dpsi(n,l,1);
        shape_matrix(i,6*nn+3)= d2psi(n,l,0);
        shape_matrix(i,6*nn+4)= d2psi(n,l,1);
        shape_matrix(i,6*nn+5)= d2psi(n,l,2);

        // At the interior points the bubble dofs
        unsigned kk = ii;
        d2_shape_dx2(linternpts[kk],psi,bpsi,dpsi,dbpsi,d2psi,d2bpsi);
        d_shape_dx(linternpts[kk],psi,bpsi,dpsi,dbpsi);
        shape(linternpts[kk],psi,bpsi);
        shape_matrix(i,18+kk)  =  psi(n,l);
       }
     }
    // The bubble basis functions
    else
     {
      // The points k
      unsigned k = i -18;
      // basis at kth internal point
      for(unsigned ii=0;ii<3;++ii)
       {
        // At the nodes the Hermite dofs
        unsigned nn=ii;
        d2_shape_dx2(lvertices[nn],psi,bpsi,dpsi,dbpsi,d2psi,d2bpsi);
        d_shape_dx(lvertices[nn],psi,bpsi,dpsi,dbpsi);
        shape(lvertices[nn],psi,bpsi);
        shape_matrix(i,6*nn)  = bpsi(k,0);
        shape_matrix(i,6*nn+1)= dbpsi(k,0,0);
        shape_matrix(i,6*nn+2)= dbpsi(k,0,1);
        shape_matrix(i,6*nn+3)= d2bpsi(k,0,0);
        shape_matrix(i,6*nn+4)= d2bpsi(k,0,1);
        shape_matrix(i,6*nn+5)= d2bpsi(k,0,2);

        // At the interior points the bubble dofs
        unsigned kk = ii;
        d2_shape_dx2(linternpts[kk],psi,bpsi,dpsi,dbpsi,d2psi,d2bpsi);
        d_shape_dx(linternpts[kk],psi,bpsi,dpsi,dbpsi);
        shape(linternpts[kk],psi,bpsi);
        shape_matrix(i,18+kk)  =  bpsi(k,0);
       }
     }
   }
//  for(unsigned i=0;i<3;++i)
//   {
//    // Fill in the Bell dofs
//    d2_shape_dx2(lvertices[i],psi,bpsi,dpsi,dbpsi,d2psi,d2psi);
//    d_shape_dx(lvertices[i],psi,bpsi,dpsi,dbpsi);
//   // shape(lvertices[i],psi,bpsi);
//    for( unsigned k=0;k<18; ++k)
//     {
//      unsigned kk = k / 6;
//      unsigned ll = k % 6;
//      shape_matrix(6*i,k) = psi(kk, ll);
//      shape_matrix(6*i+1,k) = dpsi(kk, ll,0);
//      shape_matrix(6*i+2,k) = dpsi(kk, ll,1);
//      shape_matrix(6*i+3,k) = d2psi(kk, ll,0);
//      shape_matrix(6*i+4,k) = d2psi(kk, ll,1);
//      shape_matrix(6*i+5,k) = d2psi(kk, ll,2);
//     }
//    for( unsigned k=0;k<3; ++k)
//     {
//      shape_matrix(6*i,18+k) = bpsi(k,0);
//      shape_matrix(6*i+1,18+k) = dbpsi(k,0,0);
//      shape_matrix(6*i+2,18+k) = dbpsi(k,0,1);
//      shape_matrix(6*i+3,18+k) = d2bpsi(k,0,0);
//      shape_matrix(6*i+4,18+k) = d2bpsi(k,0,1);
//      shape_matrix(6*i+5,18+k) = d2bpsi(k,0,2);
//     }
//
//    // Fill in the other dofs
//    d2_shape_dx2(linternpts[i],psi,bpsi,dpsi,dbpsi,d2psi,d2psi);
//    d_shape_dx(linternpts[i],psi,bpsi,dpsi,dbpsi);
//    shape(linternpts[i],psi,bpsi);
//    for( unsigned k=0;k<18; ++k)
//     {
//      unsigned kk = k / 6;
//      unsigned ll = k % 6;
//      shape_matrix(i+18,k) = psi(kk, ll);
//     }
//    for( unsigned k=0;k<3; ++k)
//     {
//      shape_matrix(i+18,18+k) = bpsi(k, 0);
//     }
//   }

  //Output
  for(unsigned i=0;i<21;++i)
   {
    for( unsigned k=0;k<21; ++k)
     {
       std::cout<<(fabs(shape_matrix(i,k))>tol ? shape_matrix(i,k):0.0)
                <<(k==20?"\n":" ");
     }
   }
 }

// Return the global dofs of a function on an element
void check_function(const ExactSolnFctPt&
get_analyticfunction, const double& tol)
 {
  // Get the matrices
  DenseMatrix<double> a_matrix (36,36,0.0);
  DenseMatrix<double> b_matrix (21,36,0.0);
  DenseMatrix<double> d_matrix (21,21,0.0);
  DenseMatrix<double> conversion_matrix (21,36,0.0);
  DenseMatrix<double> gl2basic_matrix (21,36,0.0);
  monomial_to_basic_matrix(a_matrix);
  basic_to_local_matrix(b_matrix);
  local_to_global_matrix(d_matrix);

  //Get the dofs
  Vector<double> dofs(21,0.0);
  dofs=get_global_dofs(get_analyticfunction);

  // Fill in conversion matrix
  for(unsigned i=0;i<21;++i)
   for(unsigned j=0;j<36;++j)
     for(unsigned k=0;k<21;++k)
      for(unsigned l=0;l<36;++l)
        conversion_matrix(i,j)+=d_matrix(i,k)*b_matrix(k,l)*a_matrix(l,j);

 //Loop over grid points
 const unsigned n_points=5;
 for(unsigned k=0;k<n_points;++k)
  {
   double s0(1.0/(n_points-1)*k);
   for(unsigned l=0;l<n_points-k;++l)
    {
     double s1(1.0/(n_points-1)*l);
     // The position on the basic element
     Vector<double> s(2,0.0); s[0]=s0; s[1]=s1;

    // Now Sum over shape functions
    double aprox_w(0.0);
    Shape p7 (36);
    full_basis_monomials(s,p7);

    for(unsigned i=0;i<21;++i)
     {
     // Get the shape functions
      // Loop over the conversion matrix
      for (unsigned k=0; k<36; ++k)
       {
        // Sum over
        aprox_w+=dofs[i]*conversion_matrix(i,k)*p7(k);
       }
     }

     // Compare
     Vector<double> w_exact(6,0.0),x(2);
     f_k(s,x);
     (*get_analyticfunction)(x,w_exact);
      if( fabs(aprox_w - w_exact[0])>tol )
        std::cout<<"Non zero difference at: "<<s<<" "<<aprox_w-w_exact[0]<<"\n";
    }
   }

}

// Return the global dofs of a function on an element
void check_basis_function(const ExactSolnFctPt&
get_analyticfunction, const double& tol)
 {
  // Get the matrices
  //Get the dofs
  Vector<double> dofs(21,0.0);
  dofs=get_basis_global_dofs(get_analyticfunction);
 //Loop over grid points
 const unsigned n_points=5;
 for(unsigned k=0;k<n_points;++k)
  {
   double s0(1.0/(n_points-1)*k);
   for(unsigned l=0;l<n_points-k;++l)
    {
     double s1(1.0/(n_points-1)*l);
     // The position on the basic element
     Vector<double> s(2,0.0); s[0]=s0; s[1]=s1;
     Shape psi(3,6),psib(3,1);
     DShape dpsi(3,6,2),dbpsi(3,1,2);
     DShape d2psi(3,6,3),d2bpsi(3,1,3);
     d2_shape_dx2(s,psi,psib,dpsi,dbpsi, d2psi,d2bpsi);
     d_shape_dx(s,psi,psib,dpsi,dbpsi);
     shape(s,psi,psib);
      // Now Sum over shape functions
      Vector<double> aprox_w(6,0.0);

       // Get the shape functions
       // Loop over the conversion matrix
       // Sum over
      for(unsigned k=0;k<18;++k)
        {
        aprox_w[0]+=dofs[k]*psi(k / 6,k % 6);
        aprox_w[1]+=dofs[k]*dpsi(k / 6,k % 6,0);
        aprox_w[2]+=dofs[k]*dpsi(k / 6,k % 6,1);
        aprox_w[3]+=dofs[k]*d2psi(k / 6,k % 6,0);
        aprox_w[4]+=dofs[k]*d2psi(k / 6,k % 6,2);
        aprox_w[5]+=dofs[k]*d2psi(k / 6,k % 6,1);
        }
      for(unsigned k=0;k<3;++k)
        {
        aprox_w[0]+=dofs[18+k]*psib(k ,0);
        aprox_w[1]+=dofs[18+k]*dbpsi(k ,0,0);
        aprox_w[2]+=dofs[18+k]*dbpsi(k ,0,1);
        aprox_w[3]+=dofs[18+k]*d2bpsi(k ,0,0);
        aprox_w[4]+=dofs[18+k]*d2bpsi(k ,0,2);
        aprox_w[5]+=dofs[18+k]*d2bpsi(k ,0,1);
        }
       // Compare
       Vector<double> w_exact(6,0.0),x(2); f_k(s,x);
       (*get_analyticfunction)(x,w_exact);
       for(unsigned i=0;i<6;++i)
        if( fabs(aprox_w[i]- w_exact[i])>tol )
          std::cout<<"Non zero difference at: "<<i<<","<<s<<" "<<aprox_w[i]-w_exact[i]<<"\n";
    }
   }
}

void output_to_mathematica_graphics()
 {
 // Get copies of the private data
 VertexList vertices = get_vertices(); 
 double sobar(get_s_obar()), subar(get_s_ubar());

 // Open the file
 std::ofstream some_file;
 char filename[100];

 sprintf(filename,"RESLT/tests.m");
 some_file.open(filename);

 // Output a mathematica plotting script
 some_file << "<<JavaGraphics`\n";

 // TESTS FOR PSI AND CHI
 some_file << "Print[\"Triangle and normal using \\[Chi](s) (red)\"]\n";
 some_file << "Show[Graphics[{Red,\n";
 some_file <<"Polygon[{\n";
 // The Curved Domain
 for (unsigned i=0;i<21;++i)
 {
  double s = i*(sobar-subar)/20. + subar;
  // Get chi
  Vector<double> mychi(2); chi(s,mychi);
  some_file<<"{"<<std::fixed
           <<(mychi[0]) <<","<<(mychi[1]) <<"},";
 }
 some_file<<"{"<<vertices[2][0]<<","<<vertices[2][1]<<"}}],"<<"\n";

 // The Normals
 for (unsigned i=0;i<21;++i)
 {
  double s = i*(sobar-subar)/20. + subar;
  // Get chi
  Vector<double> mychi(2); chi(s,mychi);
  Vector<double> myd_chi(2); d_chi(s,myd_chi);
  some_file<<"Arrow[{{"<<std::fixed
           <<(mychi[0])<<","<<(mychi[1])<<"},"
           <<"{"<<(mychi[0] -0.2*myd_chi[1])
           <<","<<(mychi[1] +0.2*myd_chi[0]) << "}}],\n";
 }
 some_file<<"{Point[{0,0}]}}]]"<<"\n";


 some_file << "Print[\"Triangle and normal using \\[Psi](s) (blue)\"]\n";
 some_file << "Show[Graphics[{Blue,\n";
 some_file <<"Polygon[{\n";

 // The Curved Domain
 for (unsigned i=0;i<21;++i)
 {
  double s1 =1 -i/20.;
  // Get psi
  Vector<double> mypsi_h(2); psi_h(s1,mypsi_h);
  some_file<<"{"<<std::fixed
           <<(mypsi_h[0]) <<","<<(mypsi_h[1]) <<"},";
 }
 some_file<<"{"<<vertices[2][0]<<","<<vertices[2][1]<<"}}],"<<"\n";

 // The Normals
 for (unsigned i=0;i<21;++i)
 {
  double s =1 -i/20.;
  Vector<double> mypsi(2); psi(s,mypsi);
  Vector<double> myd_psi(2); d_psi(s,myd_psi);
  some_file<<"Arrow[{{"<<std::fixed
           <<(mypsi[0])<<","<<(mypsi[1])<<"},"
           <<"{"<<(mypsi[0] +0.2*myd_psi[1])
           <<","<<( mypsi[1] -0.2*myd_psi[0]) << "}}],\n";
 }
 some_file<<"{Point[{0,0}]}}]]"<<"\n";
  
 // Show Vectors
 some_file << "Print[\"Triangle and normal using \\[Chi](s) (black)\"]\n";
 some_file << "Show[Graphics[{Black,\n";
 // Side 0
 some_file << "Arrow[{\n{"<<vertices[0][0]<<","<<vertices[0][1]<<"},\n";
 some_file << "{"<<vertices[0][0]+A1(0)<<","<<vertices[0][1]+A1(1)
           <<"}}],\n";
 // Side 1
 some_file << "Arrow[{\n{"<<vertices[1][0]<<","<<vertices[1][1]<<"},\n";
 some_file << "{"<<vertices[1][0]+B2(0)<<","<<vertices[1][1]+B2(1)
           <<"}}],\n";

 // Tangents
 some_file << "Arrow[{\n{"<<vertices[0][0]<<","<<vertices[0][1]<<"},\n";
 some_file << "{"<<vertices[0][0]+A2(0)<<","<<vertices[0][1]+A2(1)
           <<"}}],\n";

 some_file << "Arrow[{\n{"<<vertices[1][0]<<","<<vertices[1][1]<<"},\n";
 some_file << "{"<<vertices[1][0]+B1(0)<<","<<vertices[1][1]+B1(1)
           <<"}}],\n";
 // Side 2
 some_file << "Arrow[{\n{"<<vertices[0][0]<<","<<vertices[0][1]<<"},\n";
 some_file << "{"<<vertices[0][0]+A1(0)-B2(0)<<","
           <<vertices[0][1]+A1(1)-B2(1)<<"}}]\n";
 some_file<<"}]]"<<"\n";

 some_file << "Print[\"Triangle and normal using \\[Psi](s) (blue)\"]\n";
 some_file << "Show[Graphics[{Blue,\n";
 some_file <<"Polygon[{\n";

 // The Curved Domain
 for (unsigned i=0;i<21;++i)
 {
  // Get psi
  double s1 =1 -i/20.;
  Vector<double> mypsi_h(2); psi_h(s1,mypsi_h);
  some_file<<"{"<<std::fixed
           <<(mypsi_h[0]) <<","<<(mypsi_h[1]) <<"},";
 }
 some_file<<"{"<<vertices[2][0]<<","<<vertices[2][1]<<"}}],"<<"\n";

 // The Normals
 for (unsigned i=0;i<21;++i)
 {
  double s =1 -i/20.;
  Vector<double> mypsi(2); psi(s,mypsi);
  Vector<double> myd_psi(2); d_psi(s,myd_psi);
  some_file<<"Arrow[{{"<<std::fixed
           <<(mypsi[0])<<","<<(mypsi[1])<<"},"
           <<"{"<<(mypsi[0] +0.2*myd_psi[1])
           <<","<<(mypsi[1] -0.2*myd_psi[0]) << "}}],\n";
 }
 some_file<<"{Point[{0,0}]}}]]"<<"\n";

// // Ouput A1 A2 B1 B2
// std::cout<<"A1:\n";
// std::cout<<A1()<<"\n";
// std::cout<<"A2:\n";
// std::cout<<A2()<<"\n";
// std::cout<<"B1:\n";
// std::cout<<B1()<<"\n";
// std::cout<<"B2:\n";
// std::cout<<B2()<<"\n";

 // TESTS FOR FK
 some_file << "Print[\"A number of points on the reference triangle (Blue)\"]\n";
 some_file << "Show[Graphics[{Blue,\n";
 some_file <<"Line[{\n";

 // The Curved Domain
 for (double s0=1.0;s0>=0.0-1e-11;s0-=0.1)
 {
  for (double s1=1.0-s0;s1>=0.0-1e-11;s1-=0.1)
   {
    some_file<<"{"<<std::fixed
             <<s0<<","<<s1 <<"},\n";
   }
 }
 some_file<<"{0,0}}],";
 some_file<<"Line[{{1,0},{0,1},{0,0},{1,0}}]}]]"<<"\n";

 // Should be affine for vertices
 Vector<double> s(2,0.0),x(2); s[0]=1; s[1]=0; f_k(s,x);
 std::cout<<"These should be the vertices:"<<x<<"\n";
 s[0]=0; s[1]=1; f_k(s,x);
 std::cout<<"             "<<x<<"\n";
 s[0]=0; s[1]=0; f_k(s,x);
 std::cout<<"             "<<x<<"\n";

 some_file << "Print[\"The same points on the curved triangle (Red)\"]\n";
 some_file << "Show[Graphics[{Red,\n";
 some_file <<"Point[{\n";

 for (double i=0;i<11;++i)
 {
  // Initialise shape
  Vector<double> s(2,0.0);
  s[0]= 0.1*i;
  for (double j=0;i+j<11;++j)
   {
    s[1]= 0.1*j; f_k(s,x);
    // Output to file
    some_file<<"{"<<std::fixed
             <<x[0]<<","<<x[1] <<"},\n";
   }
 }
 some_file<<"{0,0}}],";
 some_file<<"Line[{{"
          <<vertices[0][0]<<","<<vertices[0][1]<<"},{"
          <<vertices[1][0]<<","<<vertices[1][1]<<"},{"
          <<vertices[2][0]<<","<<vertices[2][1]<<"},{"
          <<vertices[0][0]<<","<<vertices[0][1]<<"}"
          <<"}]}]]"<<"\n";
 some_file.close();
}

void check_jacobian_and_hessian()
 {
 // Now some tests on the Jacobian
 Vector<double> fk1_s0(5,0.0),fk2_s0(5,0.0),fk1_s1(5,0.0),fk2_s1(5,0.0);
 Vector<double> locus(2,0.0);
 locus[0]=.19; locus[1]=0.7910;

 // Grid spacing
 double h=0.1;

 for (unsigned i=0; i<3;++i) //2 pts in each direction
 {
  // The negative points in s0 direction
  Vector<double> s=locus,x(2);
  s[0]-=i*h; f_k(s,x);
  fk1_s0[2-i]=x[0]; 
  fk2_s0[2-i]=x[1];
  // Now the positive points in s0 direction
  s= locus;
  s[0]+=i*h; f_k(s,x);
  fk1_s0[2+i]=x[0];
  fk2_s0[2+i]=x[1];
  // The negative points in s1 direction
  s=locus;
  s[1]-=i*h; f_k(s,x);
  fk1_s1[2-i]=x[0];
  fk2_s1[2-i]=x[1];
  // Now the positive points in s1 direction
  s= locus;
  s[1]+=i*h; f_k(s,x);
  fk1_s1[2+i]=x[0]; 
  fk2_s1[2+i]=x[1];
 }

 // Check jacobian and Hessian are correct
 DenseMatrix<double> jacobian(2,2,0.0);
 get_basic_jacobian(locus,jacobian);
  
 std::cout<<"The Jacobian entries compared to FD jacobian:\n"
          <<jacobian(0,0)-dfdx_fd4(fk1_s0,h)<<" "
          <<jacobian(0,1)-dfdx_fd4(fk1_s1,h)<<"\n"
          <<jacobian(1,0)-dfdx_fd4(fk2_s0,h)<<" "
          <<jacobian(1,1)-dfdx_fd4(fk2_s1,h)<<"\n";

 // Now some tests on the Hessian
 Vector<double> J11_s0(5,0.0), J11_s1(5,0.0), J12_s0(5,0.0), J12_s1(5,0.0),
                J21_s0(5,0.0), J21_s1(5,0.0), J22_s0(5,0.0), J22_s1(5,0.0);

 for (unsigned i=0; i<3;++i) //2 pts in each direction
 {
  // The negative points in s0 direction
  Vector<double> s=locus;
  s[0]-=i*h;
  get_basic_jacobian(s,jacobian);
  J11_s0[2-i]=jacobian(0,0); 
  J21_s0[2-i]=jacobian(1,0);
  J12_s0[2-i]=jacobian(0,1); 
  J22_s0[2-i]=jacobian(1,1);
  // Now the positive points in s0 direction
  s= locus;
  s[0]+=i*h;
  get_basic_jacobian(s,jacobian);
  J11_s0[2+i]=jacobian(0,0); 
  J21_s0[2+i]=jacobian(1,0);
  J12_s0[2+i]=jacobian(0,1); 
  J22_s0[2+i]=jacobian(1,1);
  // The negative points in s1 direction
  s=locus;
  s[1]-=i*h;
  get_basic_jacobian(s,jacobian);
  J11_s1[2-i]=jacobian(0,0); 
  J21_s1[2-i]=jacobian(1,0);
  J12_s1[2-i]=jacobian(0,1); 
  J22_s1[2-i]=jacobian(1,1);
  // Now the positive points in s1 direction
  s= locus;
  s[1]+=i*h;
  get_basic_jacobian(s,jacobian);
  J11_s1[2+i]=jacobian(0,0); 
  J21_s1[2+i]=jacobian(1,0);
  J12_s1[2+i]=jacobian(0,1); 
  J22_s1[2+i]=jacobian(1,1);
 }
 
 // Output Hessian entries
 RankThreeTensor<double> hess(2,2,2,0);
 get_basic_hessian(locus,hess);
 std::cout<< "Hessian Entries minus finite differenced:\n";
 std::cout<< hess(0,0,0)-dfdx_fd4(J11_s0,h)<< "," 
          << hess(0,1,0)-dfdx_fd4(J12_s0,h)<<"," 
          << hess(1,0,0)-dfdx_fd4(J21_s0,h)<< "," 
          << hess(1,1,0)-dfdx_fd4(J22_s0,h)<<"\n"
          << hess(0,0,1)-dfdx_fd4(J11_s1,h)<< "," 
          << hess(0,1,1)-dfdx_fd4(J12_s1,h)<<"," 
          << hess(1,0,1)-dfdx_fd4(J21_s1,h)<< "," 
          << hess(1,1,1)-dfdx_fd4(J22_s1,h)<<"\n";
 }

void check_constant_consistency()
 {
 // Check the consistency of the constants
 std::cout<<"\nCheck constant consistency - using defining relations (these "
          <<"should be zero):\n";

 // Check a tilde alpha
 std::cout<<"(";
 for(unsigned i=0; i<2;++i)
 std::cout<< B2(i) - A1(i) -a_tilde_1()*A1(i)
             -a_tilde_2()*A2(i) <<(i==1? ")":","); 
 std::cout<<"\n";

 // check a tilde tilde alpha
 std::cout<<"(";
 for(unsigned i=0; i<2;++i)
 std::cout<< -B1(i) -a_tildetilde_1()*A1(i)
             -a_tildetilde_2()*A2(i) <<(i==1? ")":","); 
 std::cout<<"\n";

 // Check b tilde alpha
 std::cout<<"(";
 for(unsigned i=0; i<2;++i)
 std::cout<< -B2(i) + A1(i) -b_tilde_1()*B1(i)
             -b_tilde_2()*B2(i) <<(i==1? ")":","); 
 std::cout<<"\n";

 // check b tilde tilde alpha
 std::cout<<"(";
 for(unsigned i=0; i<2;++i)
 std::cout<< A2(i) -b_tildetilde_1()*B1(i)
             -b_tildetilde_2()*B2(i) <<(i==1? ")":","); 
 std::cout<<"\n";

 // check c tilde alpha
 std::cout<<"(";
 for(unsigned i=0; i<2;++i)
 std::cout<< A2(i) +c_tilde_1()*B2(i)
             +c_tilde_2()*A1(i) <<(i==1? ")":","); 
 std::cout<<"\n";

 // check b tilde tilde alpha
 std::cout<<"(";
 for(unsigned i=0; i<2;++i)
 std::cout<< B1(i) +c_tildetilde_1()*B2(i)
             +c_tildetilde_2()*A1(i) <<(i==1? ")":","); 
 std::cout<<"\n";

 // Check the eccentricity parameters
 //First we check that the normals are indeed normal
 std::cout<<"\nThe altitude vectors:\n";
 std::cout<<"["<< altitude_vector_1(0)<<","<<altitude_vector_1(1)
          <<"]\n";
 std::cout<<"["<< altitude_vector_2(0)<<","<<altitude_vector_2(1)
          <<"]\n";
 std::cout<<"\nThe dot products with the tangents should be zero:\n";
 std::cout<< B2(0)*altitude_vector_1(0)
           + B2(1)*altitude_vector_1(1)<<"\n";
 std::cout<< A1(0)*altitude_vector_2(0)
           + A1(1)*altitude_vector_2(1)<<"\n";
 
 std::cout<<"\nThe altitudes:\n"
          << altitude_1()<<" "<<altitude_2()<<"\n"; 
 std::cout<<"\nThe altitudes:\n";
 std::cout<<sqrt(pow(altitude_vector_1(0),2)
                   +pow(altitude_vector_1(1),2))<<" "
          <<sqrt(pow(altitude_vector_2(0),2)
                   +pow(altitude_vector_2(1),2))<<"\n";
 std::cout<<"Check the length is equal to the altitude vector length:\n";
 std::cout<< pow(altitude_vector_1(0),2)+pow(altitude_vector_1(1),2)
           - pow(altitude_1(),2)<<"\n"; 
 std::cout<< pow(altitude_vector_2(0),2)+pow(altitude_vector_2(1),2)
           - pow(altitude_2(),2)<<"\n"; 
 std::cout<<"Now check the eccentricity parameters (2 tests each): \n";

 // The eccentricity parameter check
 {
 double B2B2(0.0), v2B2(0.0), A1A1(0.0), v1A1(0.0);
 for(unsigned i=0; i<2; ++i)
  {
   // B2 dot B2
   B2B2 += B2(i)*B2(i);
   A1A1 += A1(i)*A1(i);
   // v2 dot B2
   v2B2 +=  ( B2(i)-A1(i)-altitude_vector_1(i))*B2(i);
   v1A1 +=  (-B2(i)+A1(i)-altitude_vector_2(i))*A1(i);
  }
  // Now calculate ratio 1 and 2
  double r1(1-sqrt(A1A1/B2B2-pow(altitude_1(),2)/B2B2));
  double r2(1-sqrt(B2B2/A1A1-pow(altitude_2(),2)/A1A1));
  std::cout<<"\nTest eta_1 (should give zeros):\n";
  std::cout<< eta_1()-2*v2B2/B2B2+1<<"\n";
  std::cout<< eta_1()-2*r1+1<<"\n";

  std::cout<<"\nTest eta_2 (should give zeros):\n";
  std::cout<< eta_2()+2*v1A1/A1A1-1<<"\n";
  std::cout<< eta_2()+2*r2-1<<"\n";
 }
 }

// Check the Traces against what we would expect
void check_f3_trace(const double& tol)
{
 // Now we check G3 using the other submatrices
 // Get the B3 matrix
 DenseMatrix<double> B2 (21,6,0.0);
 basic_to_local_submatrix_2(B2);

 DenseMatrix<double> B3 (21,9,0.0);
 basic_to_local_submatrix_3(B3);

 const unsigned n_points=11;
 // Loop over some points
 for(unsigned i=0; i<n_points;++i)
  {
   // Initialise position vector
   double t=1.0*i*1./(n_points-1);
   // Now initialise the B3 matrix
   Vector<double> f3(21,0.0);
   f3=f_3(t);

   // Get 1d shape
   Shape psi(6);
   hermite_shape_1d_5(t,psi);

   // Now initialise f3
   Vector<double> f3_exact(21,0.0);
   for(unsigned j=0;j<21;++j)
    {
     // Function at node 0
     if(j==0)
      {f3_exact[j]+= psi[1];}
     // Function at node 1
     if(j==1)
      {f3_exact[j]+= psi[0];}

     // D2 w (a0)(n,t) derivative at node 0
     f3_exact[j]-=(-B2(j,0)+B2(j,1))*psi[3];
     // D2 w (a1)(n,t) derivative at node 1
     f3_exact[j]+=(B2(j,2)-B2(j,3))*psi[2];

     // D2 w (a0)(n,t) derivative at node 0
     f3_exact[j]+=(B3(j,0)-2*B3(j,1)+B3(j,2))*psi[5];
     // D2 w (a1)(n,t) derivative at node 1
     f3_exact[j]+=(B3(j,3)-2*B3(j,4)+B3(j,5))*psi[4];
    }

  // Now compare exact with aprox.
  for(unsigned j=0;j<21;++j)
   {
    if(std::abs(f3_exact[j]-f3[j])>0)
      std::cout<<"Non zero difference for f3_exact at dof: "<<j
               << " diff: "<<f3_exact[j]-f3[j]<<"\n";
   }
  }
}
};

} // end c1_curved_checks namespace
} //end namespace expansion (oomph)
#endif // C1_CURVED_CHECKS_HEADER
