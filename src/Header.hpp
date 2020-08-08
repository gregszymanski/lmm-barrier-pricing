#ifndef Header_h
#define Header_h

// STL

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <array>
#include <random>
#include <algorithm>
#include <string>
#include <functional>
#include <chrono>
#include <utility>


using namespace std;
using namespace chrono;



// Boost

#include <boost/math/distributions/normal.hpp>

using boost::math::normal;




// Eigen

#include <Eigen/Dense>

using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::LLT;
using Eigen::DiagonalMatrix;
using Eigen::Dynamic;
using Eigen::ArrayXXd;



typedef
Matrix<function<double(double)>, Dynamic, Dynamic>
Matrix_dtod;





// Macros affichage

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
#define clear_terminal "\n\n\n\n\n\n\n\n\n\n\n"
#else
#define clear_terminal "\033[H\033[2J"
#endif

// Quantités

#define max_double numeric_limits<double>::max()
#define min_double numeric_limits<double>::min()

#define max_int numeric_limits<int>::max()
#define min_int numeric_limits<int>::min()





#endif /* Header_h */
