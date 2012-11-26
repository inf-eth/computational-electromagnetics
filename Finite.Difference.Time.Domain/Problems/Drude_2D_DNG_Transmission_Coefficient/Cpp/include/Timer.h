// Reference: http://stackoverflow.com/questions/2215744/c-programming-how-can-i-get-execution-time-of-a-program-in-milliseconds
// gprof, which is part of the GNU toolkit, is an option. 
// Most POSIX systems will have it installed, and it's available under Cygwin for Windows. 
// Tracking the time yourself using gettimeofday() works fine, but it's the performance 
// equivalent of using print statements for debugging. It's good if you just want a quick 
// and dirty solution, but it's not quite as elegant as using proper tools.
// To use gprof, you must specify the -pg option when compiling with gcc as in:

// gcc -o prg source.c -pg

//Then you can run gprof on the generated program as follows:

// gprof prg > gprof.out

// By default, gprof will generate the overall runtime of your program, as well as the 
// amount of time spent in each function, the number of times each function was called,
// the average time spent in each function call, and similar information.

// There are a large number of options you can set with gprof. If you're interested,
// there is more information in the man pages or through Google.

#ifndef TIMER_H_
#define TIMER_H_

#define NULL 0

#if defined __linux__ || defined __CYGWIN__
#define __int64 long long
#include <sys/time.h>
#else
#include <Windows.h>
#endif
__int64 GetTimeMs64();
__int64 GetTimeus64();

#endif // #ifndef TIMER_H_
