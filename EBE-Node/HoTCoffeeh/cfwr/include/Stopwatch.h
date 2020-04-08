#ifndef STOPWATCH_H
#define STOPWATCH_H

#include <ctime>

class Stopwatch
{
	private:
		time_t start, end;
		time_t current_elapsed_time;
	public:
		Stopwatch() {start=clock(); end=0; current_elapsed_time=0;}
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
    sw.Start();
    **** code ****
    sw.Stop();
    **** more code ****
    sw.Start();
    **** still more code ****
    sw.Stop();
  And the result can be output using:
    cout << sw.printTime() << endl;
  To reset:
    sw.Reset();
-----------------------------------------------------------------------*/
