#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <complex>

template <class T> class matrix
{
    private:
    int size_r, size_c;
    double d;
    bool square;
    T* storage;
    T* LU;

    void calc_LU();

    public:
    void set_element(int r, int c, T val);
    T access_element(int r, int c); {return storage[r + size_r*c];}

    //constructors and destructors
    matrix(int size_r, int size_c);
    ~matrix();
}

#endif //MATRIX_H_INCLUDED