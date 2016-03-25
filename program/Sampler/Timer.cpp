#include "Timer.h"

Timer::Timer(class System* system)
{
  my_system = system;
}

int Timer::elapsedTimeSeconds()
{
  int count = std::chrono::duration_cast<seconds_type>(my_end - my_start).count();
  return count;
}

int Timer::elapsedTimeMilli()
{
  int count = std::chrono::duration_cast<milli_type>(my_end - my_start).count();
  return count;
}

int Timer::elapsedTimeMicro()
{
  int count = std::chrono::duration_cast<micro_type>(my_end - my_start).count();
  return count;
}

