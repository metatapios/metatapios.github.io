#include "progress_schedule.hpp"
#include <cassert>
#include <iostream>
int main(){
  constexpr double poll=600.0,heartbeat=1800.0;
  assert(!should_write_progress(10,10,poll,heartbeat));
  assert(!should_write_progress(10,10,2*poll,heartbeat));
  assert(should_write_progress(10,10,3*poll,heartbeat));
  assert(should_write_progress(11,10,poll,heartbeat));
  assert(should_write_progress(10,-1,0,heartbeat));
  assert(should_write_progress(10,10,5,heartbeat,true));
  std::cout<<"10-minute poll / 30-minute heartbeat policy: PASS\n";
}
