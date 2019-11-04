
//#include <boost/multiprecision/cpp_dec_float.hpp>
//#include <boost/multiprecision/cpp_complex.hpp>

#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/multiprecision/float128.hpp>
#include <iostream>
#include <quadmath.h>
#include <random>

int main()
{
    //boost::multiprecision::cpp_complex<100> input(10, 10);
    boost::multiprecision::cpp_complex_quad input(10, 10);
    boost::multiprecision::float128 output;
    output = (boost::multiprecision::float128)abs(input);
    double output_2;
    //output_2 = double(real(input));
    std::cout << std::setprecision(std::numeric_limits<boost::multiprecision::float128>::max_digits10) <<"output 1 = " << output << std::endl;
    std::cout << "max value = " << std::numeric_limits<boost::multiprecision::cpp_complex_quad>::max() << std::endl;

    int seed = 10;
    std::mt19937_64 generator(seed);
    std::uniform_real_distribution<double> uniform(0, 1);
    std::normal_distribution<double> gaussian(0.0, 4.0);

    for(int i = 0; i < 9; ++i)
    {
        std::cout << uniform(generator) << " ";
    }
    std::cout << std::endl;

    for(int i = 0; i < 10; ++i)
    {
        std::cout << gaussian(generator) << " ";
    }
    std::cout << std::endl;
    return 0;
}
