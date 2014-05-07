#ifndef __TIMER_HPP__
#define __TIMER_HPP__

#include <iostream>
#include <ctime>
class Timer
{
  std::clock_t    start;
public:
  Timer():start(std::clock())
  {
  }
   
  void elapsed()
  {
    std::cout << "Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
  }
};

#endif  // __TIMER_HPP__
