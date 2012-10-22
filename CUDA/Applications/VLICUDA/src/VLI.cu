#include <VLI.hpp>
#include <string>
#include <sstream>
using std::string;
using std::stringstream;

CVLI::CVLI (): Sign(0)
{
}

CVLI::CVLI (bool pSign, vector<short> &pNumber): Sign(pSign), Number(pNumber)
{
	Format ();
}

CVLI::CVLI (const unsigned int x)
{
	stringstream VLIstream;
	VLIstream << x;
	VLIstream >> (*this);
}

void CVLI::Format ()
{
	// Remove any preceding zeros.
	for (int i=Number.size()-1; i > 0; i--)
	{
		if (Number[i] == 0)
			Number.pop_back();
		else
			break;
	}
	// Zero should have a positive sign.
	if (Number.size() == 1 && Number[0] == 0)
		Sign = 0;
}

bool CVLI::IsZero () const
{
	if (Sign == 0 && Number.size () == 1 && Number[0] == 0)
		return true;
	else
		return false;
}

CVLI CVLI::Abs () const
{
	CVLI rAbs = (*this);
	rAbs.Sign = 0;
	return rAbs;
}

CVLI CVLI::Power10 (int n)
{
	CVLI temp;
	temp.Number.insert (temp.Number.begin(), (short)1);
	for (int i=0; i<n; i++)
		temp.Number.insert (temp.Number.begin(), (short)0);
	return temp;
}
// Using twice of rough estimation as break-point
// http://en.wikipedia.org/wiki/Methods_of_computing_square_roots

bool CVLI::CheckPrime (bool CheckBaseCase)
{
	// Base case conditions upto 11.
	if (CheckBaseCase == true)
	{
		if (Number.size() == 1)
		{
			if (Number[0] == 1 || Number[0] == 2 || Number[0] == 3 || Number[0] == 5 || Number[0] == 7)
				return true;
		}
		else if (Number.size() == 2 && Number[0] == 1 && Number[1] == 1)
			return true;
	}
	// Division by 2 and 3.
	CVLI two, three;
	two.Number.push_back (2);
	three.Number.push_back (3);
	if ( ((*this) % two).IsZero() == true || ((*this) % three).IsZero() == true )
		return false;

	CVLI UpperBound;
	CVLI temp;

	CVLI i;
	CVLI one;
	CVLI six;
	CVLI Check;
	bool odd = false;

	i.Number.push_back (1);
	one.Number.push_back (1);
	six.Number.push_back (6);
	Check.Number.push_back (5);

	UpperBound = ++(Sqrt ((*this)));
	while (Check < UpperBound)
	{
		if (((*this) % Check).IsZero() == true)
			return false;

		odd == false ? Check = (six * (i++))+one : Check = (six * i)-one;
		odd = !odd;
	}
	return true;
}

CVLI CVLI::Sqrt (CVLI& VLI, int Iterations)
{
	int n;
	int D = VLI.Number.size();
	CVLI x;
	CVLI two;
	CVLI six;
	two.Number.push_back (2);
	six.Number.push_back (6);
	
	// Determine guess.
	if (D % 2 == 1)
	{
		n = (D-1)/2;
		x = two * Power10 (n);
	}
	else
	{
		n = (D-2)/2;
		x = six * Power10 (n);
	}
	for (int i=0; i<Iterations; i++)
	{
		x = (x+(VLI/x))/two;
	}
	return x;
}
// Comparison operators.
bool CVLI::operator== (const CVLI& pVLI)
{
	if (Sign == pVLI.Sign && Number == pVLI.Number)
		return true;
	else
		return false;
}

bool CVLI::operator!= (const CVLI& pVLI)
{
	if ((*this) == pVLI)
		return false;
	else 
		return true;
}
bool CVLI::operator< (const CVLI& pVLI)
{
	// Equality check.
	if ((*this) == pVLI)
		return false;

	// Zero check.
	if ((*this).IsZero() || pVLI.IsZero())
	{
		if ((*this).IsZero())
			return !pVLI.Sign;
		else
			return Sign;
	}

	// Sign check.
	if (Sign != pVLI.Sign)
		return Sign;

	// Both positive.
	if (Sign == 0)
	{
		// Check length.
		if (Number.size() != pVLI.Number.size())
		{
			return (Number.size() < pVLI.Number.size());
		}
		else
		{
			for (int i=0; i < (int)Number.size(); i++)
			{
				if (Number[Number.size()-1-i] == pVLI.Number[pVLI.Number.size()-1-i])
					continue;
				else
					return (Number[Number.size()-1-i] < pVLI.Number[pVLI.Number.size()-1-i]);
			}
		}
	}
	// Else if both negative.
	else
	{
		// Check length.
		if (Number.size() != pVLI.Number.size())
		{
			return !(Number.size() < pVLI.Number.size());
		}
		else
		{
			for (int i=0; i < (int)Number.size(); i++)
			{
				if (Number[Number.size()-1-i] == pVLI.Number[pVLI.Number.size()-1-i])
					continue;
				else
					return !(Number[Number.size()-1-i] < pVLI.Number[pVLI.Number.size()-1-i]);
			}
		}
	}
	return false;
}

bool CVLI::operator> (const CVLI& pVLI)
{
	// Equality check.
	if ((*this) == pVLI)
		return false;

	// Zero check.
	if ((*this).IsZero() || pVLI.IsZero())
	{
		if ((*this).IsZero())
			return pVLI.Sign;
		else
			return !Sign;
	}

	// Sign check.
	if (Sign != pVLI.Sign)
		return !Sign;

	// Both positive.
	if (Sign == 0)
	{
		// Check length.
		if (Number.size() != pVLI.Number.size())
		{
			return (Number.size() > pVLI.Number.size());
		}
		else
		{
			for (int i=0; i < (int)Number.size(); i++)
			{
				if (Number[Number.size()-1-i] == pVLI.Number[pVLI.Number.size()-1-i])
					continue;
				else
					return (Number[Number.size()-1-i] > pVLI.Number[pVLI.Number.size()-1-i]);
			}
		}
	}
	// Else if both negative.
	else
	{
		// Check length.
		if (Number.size() != pVLI.Number.size())
		{
			return !(Number.size() > pVLI.Number.size());
		}
		else
		{
			for (int i=0; i < (int)Number.size(); i++)
			{
				if (Number[Number.size()-1-i] == pVLI.Number[pVLI.Number.size()-1-i])
					continue;
				else
					return !(Number[Number.size()-1-i] > pVLI.Number[pVLI.Number.size()-1-i]);
			}
		}
	}
	return false;
}

bool CVLI::operator<= (const CVLI& pVLI)
{
	if ((*this) < pVLI || (*this) == pVLI)
		return true;
	else
		return false;
}

bool CVLI::operator>= (const CVLI& pVLI)
{
	if ((*this) > pVLI || (*this) == pVLI)
		return true;
	else
		return false;
}

// Pre-post increment.
CVLI CVLI::operator++ ()
{
	Number[0]++;
	// Carry Pass.
	short Carry = 0;
	for (int i=0; i < (int)Number.size(); i++)
	{
		Number[i] += Carry;
		if (Number[i] > 9)
		{
			Carry = 1;
			Number[i] -= 10;
		}
		else
			Carry = 0;
	}
	if (Carry == 1)
		Number.push_back (1);
	return (*this);
}

// Pre-post increment.
CVLI CVLI::operator++ (int dummy)
{
	CVLI temp = (*this);
	++(*this);
	return temp;
}

// Addition and subtraction
CVLI CVLI::operator+ (const CVLI& pVLI)
{
	CVLI Sum;
	// Same signs, simple addition.
	if (Sign == pVLI.Sign)
	{
		Sum.Sign = Sign;
		Number.size()<pVLI.Number.size() ? Sum.Number = pVLI.Number : Sum.Number = Number;
		// First Pass.
		for (int i=0; i < (int)(Number.size()<pVLI.Number.size() ? Number.size() : pVLI.Number.size()); i++)
			Sum.Number[i] = Number.size()<pVLI.Number.size() ? Number[i] + Sum.Number[i]: pVLI.Number[i] + Sum.Number[i];
		// Carry Pass.
		short Carry = 0;
		for (int i=0; i < (int)(Number.size()<pVLI.Number.size() ? pVLI.Number.size() : Number.size()); i++)
		{
			Sum.Number[i] += Carry;
			if (Sum.Number[i] > 9)
			{
				Carry = 1;
				Sum.Number[i] -= 10;
			}
			else
				Carry = 0;
		}
		if (Carry == 1)
			Sum.Number.push_back (1);
		return Sum;
	}
	// Subtraction in case of different signs.
	else
	{
		Sum.Sign = (*this).Abs() > pVLI.Abs() ? (*this).Sign : pVLI.Sign;
		(*this).Abs() < pVLI.Abs() ? Sum.Number = pVLI.Number : Sum.Number = Number;
		// First Pass.
		for (int i=0; i < (int)((*this).Abs() < pVLI.Abs() ? Number.size() : pVLI.Number.size()); i++)
			Sum.Number[i] = (*this).Abs() < pVLI.Abs() ? Sum.Number[i] - Number[i]: Sum.Number[i] - pVLI.Number[i];
		// Carry Pass.
		short Carry = 0;
		for (int i=0; i < (int)((*this).Abs() < pVLI.Abs() ? pVLI.Number.size() : Number.size()); i++)
		{
			Sum.Number[i] += Carry;
			if (Sum.Number[i] < 0)
			{
				Carry = -1;
				Sum.Number[i] += 10;
			}
			else
				Carry = 0;
		}
		Sum.Format();
		return Sum;
	}
}

CVLI CVLI::operator- (const CVLI& pVLI)
{
	CVLI temp;
	temp.Sign = !pVLI.Sign;
	temp.Number = pVLI.Number;
	return ((*this)+temp);
}

// Multiplication
CVLI CVLI::operator* (const CVLI& pVLI)
{
	CVLI Product, temp;
	int z;
	Product.Number.push_back (0);

	//  this
	//x pVLI
	//-------
	//Product
	//-------

	for (int i=0; i < (int)pVLI.Number.size(); i++)
	{
		temp.Number.clear();
		for (z=0; z<i; z++)
			temp.Number.push_back (0);

		// Raw multiplication
		for (int j=0; j < (int)Number.size(); j++)
		{
			temp.Number.push_back (Number[j] * pVLI.Number[i]);
		}

		// Carry manipulation.
		short Carry = 0;
		for (int k=0; k < (int)temp.Number.size(); k++)
		{
			temp.Number[k] += Carry;
			if (temp.Number[k] > 9)
			{
				Carry = temp.Number[k] / 10;
				temp.Number[k] %= 10;
			}
			else
				Carry = 0;
		}
		if (Carry != 0)
			temp.Number.push_back (Carry);
		temp.Format();
		Product = Product + temp;
	}
	Product.Sign = Sign ^ pVLI.Sign;
	return Product;
}

// Division
CVLI CVLI::operator/ (const CVLI& pVLI)
{
	CVLI Quotient;
	CVLI Dividend;
	CVLI Remainder;
	CVLI Product;
	CVLI Divisor = pVLI;
	Divisor.Sign = 0;

	bool ZeroFlag = false;
	bool InsertedFlag = false;
	for (int i=0; i<(int)Number.size(); i++)
	{
		if (ZeroFlag == true && Number[Number.size()-1-i] == (short)0)
		{
			Quotient.Number.insert (Quotient.Number.begin(), (short)0);
			continue;
		}
		else
		{
			Dividend.Number.insert (Dividend.Number.begin(), 1, Number[Number.size()-1-i]);
			ZeroFlag = false;
		}
		if (Dividend < Divisor)
		{
			if (InsertedFlag == true)
				Quotient.Number.insert (Quotient.Number.begin(), (short)0);
			continue;
		}
		else
		{
			int j = 1;
			Product = Divisor;
			for (; j<10; j++)
			{
				if (Product == Dividend)
					break;
				if ((Product+Divisor) > Dividend)
					break;
				Product = Product + Divisor;
			}
			Remainder = Dividend - Product;
			Quotient.Number.insert (Quotient.Number.begin(), 1, (short)j);
			InsertedFlag = true;
			if (Remainder.IsZero() == true)
			{
				ZeroFlag = true;
				Dividend.Number.clear();
			}
			else
				Dividend = Remainder;
		}
	}
	Quotient.Sign = Sign ^ pVLI.Sign;
	return Quotient;
}

// modulus
CVLI CVLI::operator% (const CVLI& pVLI)
{
	CVLI Quotient;
	CVLI Dividend;
	CVLI Remainder;
	CVLI Product;
	CVLI Divisor = pVLI;
	Divisor.Sign = 0;

	bool ZeroFlag = false;
	bool InsertedFlag = false;
	for (int i=0; i<(int)Number.size(); i++)
	{
		if (ZeroFlag == true && Number[Number.size()-1-i] == (short)0)
		{
			Quotient.Number.insert (Quotient.Number.begin(), (short)0);
			continue;
		}
		else
		{
			Dividend.Number.insert (Dividend.Number.begin(), 1, Number[Number.size()-1-i]);
			ZeroFlag = false;
		}
		if (Dividend < Divisor)
		{
			if (InsertedFlag == true)
				Quotient.Number.insert (Quotient.Number.begin(), (short)0);
			continue;
		}
		else
		{
			int j = 1;
			Product = Divisor;
			for (; j<10; j++)
			{
				if (Product == Dividend)
					break;
				if ((Product+Divisor) > Dividend)
					break;
				Product = Product + Divisor;
			}
			Remainder = Dividend - Product;
			Quotient.Number.insert (Quotient.Number.begin(), 1, (short)j);
			InsertedFlag = true;
			if (Remainder.IsZero() == true)
			{
				ZeroFlag = true;
				Dividend.Number.clear();
			}
			else
				Dividend = Remainder;
		}
	}
	if (Dividend.Number.empty() == false)
		Remainder = Dividend;
	return Remainder;
}

// << and >> operators.
ostream& operator<<(ostream& out, const CVLI& pVLI)
{
	out << (pVLI.Sign?"-":"");
	for(int i= pVLI.Number.size()-1; i >= 0; i--)
		out << pVLI.Number[i];
	return out;
}
istream& operator >> (istream& in, CVLI& pVLI)
{
	pVLI.Number.clear();
	string input;
	in >> input;

	// Sign check.
	if (input.size() > 0 && input[0]=='-')
	{
		pVLI.Sign = 1;
		input.replace (0, 1, "");
	}
	else
		pVLI.Sign = 0;

	// Remove any preceding zeros.
	for (int i=pVLI.Sign; i < (int)input.size(); i++)
	{
		if (input[0] == '0')
			input.replace (0, 1, "");
		else
			break;
	}

	for (int i=input.size()-1; i >= 0; i--)
	{
		if ((unsigned short)(input.c_str()[i]-48) <= 9U)
			pVLI.Number.push_back ((unsigned short)(input.c_str()[i]-48));
		else
			continue;
	}
	// If zero string.
	if (input.size() == 0)
	{
		pVLI.Number.push_back ((unsigned short)(0));
		pVLI.Sign = 0;
	}
	return in;
}
