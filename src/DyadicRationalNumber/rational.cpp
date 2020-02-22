#include "DyadicRationalNumber/rational.h"

std::ostream& operator<<(std::ostream& os, DyadicRationalNumber obj) {
	os << obj.getGMP();
	return os;
}
