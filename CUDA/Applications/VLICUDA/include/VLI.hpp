#ifndef __VLI_H
#define __VLI_H

#include <vector>
#include <iostream>
using std::vector;
using std::ostream;
using std::istream;

#include <cutil.h>

class CVLI
{
	bool Sign;			// 0=positive, 1=negative
	vector<short> Number;
	void Format ();		// Properly format the Number.
	bool IsZero () const;		// Is this zero?
	CVLI Abs () const;
	CVLI Power10 (int n);

	public:
	CVLI ();
	CVLI (bool, vector<short> &);		// VLI from sign and number.
	// VCLI from unsigned int.
	CVLI (const unsigned int);

	// Check for primality.
	bool CheckPrime (bool=false);
	CVLI Sqrt (CVLI&, int = 5);	// Default number of iterations is 5. More iterations = better accuracy.

	// Comparison operators.
	bool operator == (const CVLI&);
	bool operator != (const CVLI&);
	bool operator < (const CVLI&);
	bool operator > (const CVLI&);
	bool operator <= (const CVLI&);
	bool operator >= (const CVLI&);

	// Pre-post increment.
	CVLI operator ++ ();
	CVLI operator ++ (int);

	// Addition and Subtraction.
	CVLI operator + (const CVLI&);
	CVLI operator - (const CVLI&);

	// Multiplication and division.
	CVLI operator * (const CVLI&);
	CVLI operator / (const CVLI&);

	// modulus
	CVLI operator % (const CVLI&);

	// Input and output.
	friend ostream& operator << (ostream&, const CVLI&);
	friend istream& operator >> (istream&, CVLI&);
};
#endif // __VLI_H
