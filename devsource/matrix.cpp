#include "matrix.hpp"

template <class T>matrix::matrix(int r, int c)
{   
    size_r = r;
    size_c = c;
    T storage = new T[size_r*size_c]; //original matrix
    T LU = new LU[size_r*size_c]; //triangular matrix
    if (size_r == size_c)
    {
        square = true;
    }
    else
    {
        square = false;
    }
    d = 1.;
    
}

template <class T>matrix::~matrix()
{
    delete[] storage;
    delete[] LU;
}

void template<class T>matrix::set_element(int r, int c, T val);
{
    storage[r + r_size*c] = val;
}

void template<class T>matrix::calc_LU();
{
    //LU decomposition basically copied from numerical recipies 3rd edition
    const double TINY = 1.0e-40;
    int r, rmax, c, k;
    double BIG, temp;
    vector<double> vv(size_r);
    vector<int> index(size_r);

    for(r = 0; r < pow(size_r, 2); ++r)
    {
        LU[r] = storage[r];
    }

    for (r = 0; i < size_r; ++r)
    {
        BIG = 0.;
        for(c = 0; c <size_r; ++c)
        {
            if((temp = std::abs(LU[r + size_r*c])) > BIG)
            {
                BIG = temp;
            }
        }
        if (BIG == 0.) 
        {
            printf("Singular matrix \n");
            exit(7);
        }
        vv[r] = 1./BIG;
    }

    for(k = 0; k < size_r; ++k)
    {
        BIG = 0.;
        for (r = k; r < size_r; ++r)
        {
            temp = vv[r]*std::abs(LU[r + size_r*k]);
            if (temp > BIG)
            {
                BIG = temp;
                rmax = r;
            }
        }
        if (k != rmax)
        {
            for (c = 0; c < size_r; ++c)
            {
                temp = LU[rmax + c*size_r];
                LU[rmax + c*size_r] = LU[k + c*size_r];
                LU[k + c*size_r] = temp;
            }
            d *= -1.;
            vv[rmax] = vv[k];
        }
        index[k] = rmax;
        if(LU[k + k*size_r])
        {
            LU[k + k*size_r] = TINY;
        }

        for (r = k + 1; r < size_r; ++r)
        {
            //continue from here page 53
        }
    }
}
