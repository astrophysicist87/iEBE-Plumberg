#ifndef CPSTOPWATCH_H
#define CPSTOPWATCH_H

#include <ctime>

class CPStopwatch
{
	private:
		time_t start, end;
		time_t current_elapsed_time;
	public:
		CPStopwatch() {start=clock(); end=0; current_elapsed_time=0;}
		void Start() {start=clock();}
		void Stop() {end=clock(); current_elapsed_time+=(end - start);}
		void Reset() {start=clock(); end=0; current_elapsed_time=0;}
		double printTime() {return ((double)current_elapsed_time) / CLOCKS_PER_SEC;}
};

#endif

/*-----------------------------------------------------------------------
  Usage:
  Declare a class as:
    Stopwatch sw;
  Then the time a piece of code takes can be recorded as:
    sw.tic();
    **** code ****
    sw.toc();
  And the result can be outputted using:
    cout << sw.takeTime() << endl;
-----------------------------------------------------------------------*/
