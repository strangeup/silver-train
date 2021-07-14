#include "generic.h"

double integrate_unit_triangle(Integral* my_integral_pt){
    // Create a TElement<2,2> to provide the local-to-global mapping
    unsigned const NNODE_1D=2, DIM=2, n_position_type=1;
    TElement<DIM,NNODE_1D> t_element;
    // Vertices for the unit triangle
    auto const vertices = std::vector<std::vector<double>> ({{1,0},{0,1},{0,0}});
    // Test data for equilateral triangle, to verify that the local-gloabl mapping is correct
    // auto vertices = std::vector<std::vector<double>> ({{-0.5,0},{0.5,0},{0,0.866025404}});
    if (vertices.size() != t_element.nnode()) {
      oomph_info<<"Number of nodes of element does not match length of vertices vector"<<std::endl;
      throw std::logic_error("Number of nodes of element does not match length of vertices vector");
    } 
    
    // Assign the vertices
    for (unsigned i_vtx = 0; i_vtx < t_element.nvertex_node(); ++i_vtx)
    {
      auto vnode = t_element.construct_node(i_vtx);
      auto vertex = vertices[i_vtx];
      for(unsigned x_itr = 0; x_itr<DIM; ++x_itr ) {
        vnode->x(x_itr) = vertex[x_itr];
      }
    }
    
    // Integrate it
    auto integral_value = 0.0;
    // auto my_integral_pt = &tg4;
    for(unsigned i_ipt=0; i_ipt<my_integral_pt->nweight(); ++i_ipt)
    {
      auto s = Vector<double>(DIM);
      for(unsigned x_itr = 0; x_itr<DIM; ++x_itr ) {
       s[x_itr] = my_integral_pt->knot(i_ipt,x_itr);
      }
      
      #ifdef VERBOSE
      for(unsigned x_itr = 0; x_itr<DIM; ++x_itr ) {
        std::cout<<s[x_itr]<<" ";
      }
      std::cout<<std::endl;
      #endif
    
      // Fill in Jacobian
      Shape psi(t_element.nnode(), n_position_type);
      DShape dpsi_ds(t_element.nnode(), n_position_type, DIM);
      t_element.dshape_local(s, psi, dpsi_ds);
      DenseMatrix<double> jacobian(DIM, DIM, 0.0), ijacobian(DIM, DIM, 0.0);
    
      for (unsigned l=0; l < t_element.nnode();  ++l)
      {
      for (unsigned x_itr =0; x_itr < DIM ; ++x_itr)
      {
        for (unsigned x2_itr =0; x2_itr < DIM ; ++x2_itr)
        {
          for (unsigned i_type=0; i_type < n_position_type; ++i_type) {
            jacobian(x_itr,x2_itr) += t_element.raw_nodal_position_gen(l, i_type ,x2_itr)*dpsi_ds(l, i_type ,x_itr);
          }
        }
      }
      }
    
      //Get the integral weight
      double const w = my_integral_pt->weight(i_ipt);
      double const J = t_element.invert_jacobian_mapping(jacobian, ijacobian);
      
      //Premultiply the weights and the Jacobian
      double W = w*J;
      integral_value += W;
    }
    
    oomph_info << "Integral value is: " << integral_value << std::endl;
    return integral_value;
}

// Main
int main()
{
using namespace oomph;
// Use Polymorphism to create a list of integral pt
Vector<Integral*> integral_pt_list;
TGauss<2,2> tg2;
integral_pt_list.push_back(&tg2);
TGauss<2,3> tg3;
integral_pt_list.push_back(&tg3);
TGauss<2,4> tg4;
integral_pt_list.push_back(&tg4);
TGauss<2,5> tg5;
integral_pt_list.push_back(&tg5);
TGauss<2,9> tg9;
integral_pt_list.push_back(&tg9);
TGauss<2,13> tg13;
integral_pt_list.push_back(&tg13);
TGauss<2,16> tg16;
integral_pt_list.push_back(&tg16);

// Open Trace file
std::ofstream Trace_file;
Trace_file.open("RESLT/trace_integral.dat");

// Now loop over list
for(unsigned i_ipt=0;i_ipt<integral_pt_list.size();++i_ipt)
 {
  // Initialize
  double triangle_area = 0.0;
  const unsigned nweight = integral_pt_list[i_ipt]->nweight();
  // Loop weights
  for(unsigned iweight=0;iweight<nweight;++iweight)
   {
    // Get weight and add
    triangle_area+=integral_pt_list[i_ipt]->weight(iweight);
   }
  oomph_info << "Area of unit triangle by sum: " << triangle_area <<std::endl;
  // Now do it properly
  auto triangle_area2 = integrate_unit_triangle(integral_pt_list[i_ipt]);
  oomph_info << "Area of unit triangle by integral: " << triangle_area2 <<std::endl;
  Trace_file <<std::setprecision(15)<< triangle_area2 <<std::endl;
 }

// Close Tracefile
Trace_file.close();
}
