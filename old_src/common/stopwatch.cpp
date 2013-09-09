/*
 * stopwatch.cpp
 *
 *  Created on: 17-apr-2009
 *      Author: M. El-Kebir
 */

#include "stopwatch.h"

StopWatch::StopWatch()
	: _starttime(clock())
	, _stoptime(clock())
{
}

StopWatch::~StopWatch()
{
}

void StopWatch::start()
{
	_starttime = clock();
}

void StopWatch::stop()
{
	_stoptime = clock();
}

