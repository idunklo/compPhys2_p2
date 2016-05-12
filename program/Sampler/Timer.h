#pragma once
#include "../System/System.h"

class Timer
{
  public:

    Timer (class System* system);

    void    startTimer	      ()  {my_start = clock::now();}
    void    stopTimer	      ()  {my_end = clock::now();}
    int	    elapsedTimeSeconds();
    int     elapsedTimeMilli  ();
    int     elapsedTimeMicro  ();

  protected:

    class System* my_system = nullptr;

    typedef std::chrono::high_resolution_clock	clock;
    typedef std::chrono::seconds        seconds_type;
    typedef std::chrono::milliseconds   milli_type;
    typedef std::chrono::microseconds   micro_type;

    std::chrono::duration<int> seconds;
    std::chrono::duration<int> milli;
    std::chrono::duration<int> micro;

    clock::time_point	      my_start;
    clock::time_point	      my_end;
};



