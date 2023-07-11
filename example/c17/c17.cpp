// This program demonstrates the basic usage of OpenTimer C++ API
// to get the critical path information of c17 benchmark.
//
// Design : c17.v
// SPEF   : c17.spef
// SDC    : c17.sdc
// Library: c17_Early.lib (MIN), c17_Late.lib (MAX)

#include <ot/timer/timer.hpp>

int main(int argc, char *argv[]) {

  // create a timer object
  ot::Timer timer;
  
  // Read design
  timer.read_celllib("c17_Early.lib", ot::MIN)
       .read_celllib("c17_Late.lib", ot::MAX)
       .read_verilog("c17.v")
       .read_spef("c17.spef")
       .read_sdc("c17.sdc");

  // get the top-5 worst critical paths
  auto paths = timer.report_timing(1);

  // for(size_t i=0; i<paths.size(); ++i) {
  //   std::cout << "----- Critical Path " << i << " -----\n";
  //   std::cout << paths[i] << '\n';
  // }

  // timer.dump_at(std::cout);

  // dump the timing graph to dot format for debugging
  // timer.dump_graph(std::cout);

  return 0;
}



















