#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <complex>

template <class T> class matrix
{
    private:
    int size_r, size_c;
    double d;
    std::vector<int> index;
    bool square, LUcheck, det_check;
    T* storage;
    T* LU;

    void calc_LU();

    public:
    void set_element(int r, int c, T val);
    T get_element(int r, int c) const {return storage[r + size_r*c];}
    void solve(T x[], T b[]);
    matrix<T> operator * (matrix const &obj);

    //constructors and destructors
    matrix(int size_r, int size_c);
    matrix(const matrix& obj);
    ~matrix();
};

template <class T>matrix<T>::matrix(int r, int c)
{   
    size_r = r;
    size_c = c;
    storage = new T[size_r*size_c]; //original matrix
    LU = new T[size_r*size_c]; //triangular matrix
    if (size_r == size_c)
    {
        square = true;
        index.resize(size_r, 0);
    }
    else
    {
        square = false;
    }
    LUcheck = false;
    det_check = false;
    d = 1.;
    
}

template <class T>matrix<T>::~matrix()
{
    delete[] storage;
    delete[] LU;
}

template <class T>matrix<T>::matrix(const matrix& obj) : size_r(obj.size_r),
size_c(obj.size_c),
d(obj.d),
index(index),
square(obj.square),
LUcheck(obj.LUcheck),
det_check(obj.det_check),
storage(new T[obj.size_c*obj.size_r]),
LU(new T[obj.size_c*obj.size_r])
{
    for(int i = 0; i < size_c*size_r; ++i)
    {
        storage[i] = obj.storage[i];
        LU[i] = obj.LU[i];
    }
}

template <class T>matrix<T>matrix<T>::operator * (matrix const &obj)
{
    T element;
    matrix return_val(size_r, obj.size_c);
    
    if (size_c != obj.size_r)
    {
        printf("invalid matrix multiplication \n");
        exit(8);
    }
    for (int r = 0; r < size_r; ++r)
    {
        for (int c = 0; c < obj.size_c; ++c)
        {
            element = 0;
            for (int s = 0; s < size_c; ++s)
            {
                element += get_element(r, s)*obj.get_element(s, c);
            }
            return_val.set_element(r, c, element);
        }
    }
    return return_val;
}

template<class T> void matrix<T>::set_element(int r, int c, T val)
{
    storage[r + size_r*c] = val;
    LUcheck = false;
    det_check = false;
}

template<class T> void matrix<T>::calc_LU()
{
    //LU decomposition basically copied from numerical recipies 3rd edition, as a result I can't really comment on it
    if (!square)
    {
        printf("Cannot LU decompose non square matrix \n");
        exit(9);
    }
    const double TINY = 1.0e-40;
    int r, rmax, c, k;
    double BIG, temp;
    std::vector<double> vv(size_r);
    std::vector<int> index(size_r);

    for(r = 0; r < pow(size_r, 2); ++r)
    {
        LU[r] = storage[r];
    }

    for (r = 0; r < size_r; ++r)
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
            temp = LU[r + size_r*k] /= LU[k + size_r*k];
            for (c = k + 1; c < size_r; ++c)
            {
                LU[r + c*size_r] -= temp*LU[k + size_r*c];
            }
        }
    }
    LUcheck = true;
}

template<class T> void matrix<T>::solve(T x[], T b[])
{
    //solves a system of equations of the form matrix * x = b for known b and matrix.
    //again pretty much copied from page 53 of numerical recipies 3rd ed.
    if (LUcheck == false)
    {
        calc_LU();
    }
    int i, ii = 0, ip, k;
    T sum;
    for(i = 0; i < size_r; ++i)
    {
        x[i] = b[i];
    }
    for (i = 0; i < size_r; ++i)
    {
        ip = index[i];
        sum = x[ip];
        x[ip] = x[i];
        if (ii != 0)
        {
            for(k = ii - 1; k < i; ++k)
            {
                sum -= LU[i + size_r*k]*x[k];
            }
        }
        else if(sum != 0.0)
        {
            ii = i + 1;
        }
        x[i] = sum;
    }
    for (i = size_r - 1; i >= 0; --i)
    {
        sum = x[i];
        for (k = i + 1; k < size_r; ++k)
        {
            sum -= LU[i + size_r*k]*x[k];
        }
        x[i] = sum/LU[i + i*size_r];
    }
}

#endif //MATRIX_H_INCLUDED