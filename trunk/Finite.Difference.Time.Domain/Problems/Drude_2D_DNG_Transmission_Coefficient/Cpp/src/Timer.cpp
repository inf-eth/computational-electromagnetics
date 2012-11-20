#include "Timer.h"
__int64 GetTimeMs64()
{
#if defined __linux__ || defined __CYGWIN__
	/* Linux */
	struct timeval tv;

	gettimeofday(&tv, NULL);

	__int64 ret = tv.tv_usec;
	/* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
	ret /= 1000;

	/* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
	ret += (tv.tv_sec * 1000);

	return ret;
#else
	/* Windows */
	FILETIME ft;
	LARGE_INTEGER li;

	/* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
	* to a LARGE_INTEGER structure. */
	GetSystemTimeAsFileTime(&ft);
	li.LowPart = ft.dwLowDateTime;
	li.HighPart = ft.dwHighDateTime;

	__int64 ret = li.QuadPart;
	ret -= 116444736000000000L; /* Convert from file time to UNIX epoch time. */
	ret /= 10000; /* From 100 nano seconds (10^-7) to 1 millisecond (10^-3) intervals */

	return ret;
#endif
}

/* Returns the amount of microseconds elapsed since the UNIX epoch. Works on both
* windows and linux. */

__int64 GetTimeus64()
{
#if defined __linux__ || defined __CYGWIN__
	/* Linux */
	struct timeval tv;

	gettimeofday(&tv, NULL);

	__int64 ret = tv.tv_usec;

	/* Adds the seconds (10^0) after converting them to microseconds (10^-6) */
	ret += (tv.tv_sec * 1000000);

	return ret;
#else
	/* Windows */
	FILETIME ft;
	LARGE_INTEGER li;

	/* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
	* to a LARGE_INTEGER structure. */
	GetSystemTimeAsFileTime(&ft);
	li.LowPart = ft.dwLowDateTime;
	li.HighPart = ft.dwHighDateTime;

	__int64 ret = li.QuadPart;
	ret -= 116444736000000000L; /* Convert from file time to UNIX epoch time. */
	ret /= 10; /* From 100 nano seconds (10^-7) to 1 microsecond (10^-6) intervals */

	return ret;
#endif
}
