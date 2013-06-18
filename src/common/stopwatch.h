/*
 * stopwatch.h
 *
 *  Created on: 17-apr-2009
 *      Author: s030858
 */

#ifndef STOPWATCH_H_
#define STOPWATCH_H_

#include <time.h>

class StopWatch
{
private:
	clock_t _starttime;
	clock_t _stoptime;

public:
	StopWatch();
	~StopWatch();

	void start();
	void stop();
	double getElapsedTime();
};

inline double StopWatch::getElapsedTime()
{
	return (_stoptime - _starttime) / (double) CLOCKS_PER_SEC;
}

#endif /* STOPWATCH_H_ */
