#include "generic.h"

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
  oomph_info << "Area of unit triangle: " << triangle_area <<std::endl;
  Trace_file <<std::setprecision(15)<< triangle_area <<std::endl;
 }
// Close Tracefile
Trace_file.close();
}
