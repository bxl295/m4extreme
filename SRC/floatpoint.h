//
//
//

#ifndef _M4EXTREME_UTILITY_FLOATPOINT_
#define _M4EXTREME_UTILITY_FLOATPOINT_

#include <iostream>

namespace m4extreme {

    class floatpoint {
    public:

        floatpoint() : _x(0.0) {
        }

        floatpoint(double x) : _x(x) {
        }

        ~floatpoint() {
        }

        floatpoint & operator=(double y) {
            _x = y;
            return *this;
        }

        bool operator<(double y) const {
            if (_isnan()) {
	        std::cerr << "floating point is not-a-number. throw out exeception!" << std::endl;
                throw (1);
            } else {
                return _x < y;
            }
        }

	bool operator>(double y) const {
            if (_isnan()) {
	        std::cerr << "floating point is not-a-number. throw out exeception!" << std::endl;
                throw (1);
            } else {
                return _x > y;
            }
        }

	std::ostream & print ( std::ostream & os ) const {
	  os << _x;
	  return os;
	}

	double value() {
	  return _x;
	}

    private:

        bool _isnan() const {
            return _x != _x;
        }


    private:
        double _x;
    };

}

std::ostream & operator << (std::ostream & os, const m4extreme::floatpoint & err) {
  return err.print(os);
}

#endif //_M4EXTREME_UTILITY_FLOATPOINT_
