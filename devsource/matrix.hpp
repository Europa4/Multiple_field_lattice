#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <complex>
#include <vector>

template <class T> class matrix
{
    private:
    int size_r, size_c;
    double d;
    std::vector<int> index;
    bool square, LUcheck, det_check;
    T* storage;
    T* LU;
    T det;

    void calc_LU();
    void calc_det();

    public:
    void set_element(int r, int c, T val);
    T get_element(int r, int c) const {return storage[r + size_r*c];}
    void solve(T x[], T b[]);
    matrix solve(T b[]);
    void resize(int new_r, int new_c);
    T get_det();
    matrix conj();
    matrix transpose();
    matrix forward_vector_multiplication(T vec[]);
    matrix backward_vector_multiplication(T vec[]);

    //operators
    matrix<T> operator * (matrix const &obj);
    matrix<T> operator = (matrix const &obj);
    matrix<T> operator + (matrix const &obj);

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
        for(int i = 0; i < size_r*size_r; ++i)
        {
            LU[i] = 0;
        }
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
index(obj.index),
square(obj.square),
LUcheck(obj.LUcheck),
det_check(obj.det_check),
storage(new T[obj.size_c*obj.size_r]),
LU(new T[obj.size_c*obj.size_r]),
det(obj.det)
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
    matrix<T> return_val(size_r, obj.size_c);
    
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

template <class T> matrix<T> matrix<T>::operator = (matrix const &obj)
{
    size_r = obj.size_r;
    size_c = obj.size_c;
    index = obj.index;
    square = obj.square;
    LUcheck = obj.LUcheck;
    det_check = obj.det_check;
    det = obj.det;
    LU = new T(size_r*size_c);
    storage = new T(size_r*size_c);
    for (int i = 0; i < size_r*size_c; ++i)
    {
        LU[i] = obj.LU[i];
        storage[i] = obj.storage[i];
    }
    return *this;
}

template <class T> matrix<T> matrix<T>::operator + (matrix const &obj)
{
    if (!(size_r == obj.size_r && size_c == obj.size_c))
    {
        printf("Error in matrix class, invalid matrix addition \n");
        exit(10);
    }
    matrix<T> return_val(size_r, size_c);
    for (int r = 0; r < size_r; ++r)
    {
        for (int c = 0; c < size_c; ++c)
        {
            return_val.set_element(r, c, get_element(r, c) + obj.get_element(r, c));
        }
    }
    return return_val;
}

template<class T> matrix<T> matrix<T>::conj()
{
    matrix<T> ret(size_c, size_r);
    for (int r = 0; r < size_r; ++r)
    {
        for (int c = 0; c < size_c; ++c)
        {
            ret.set_element(c, r, std::conj(get_element(r, c)));
        }
    }
    return ret;
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
    T BIG, temp;
    std::vector<double> vv(size_r, 0);
    std::vector<int> index(size_r, 0);
    
    d = 1.;

    for(r = 0; r < pow(size_r, 2); ++r)
    {
        LU[r] = storage[r];
    }

    for (r = 0; r < size_r; ++r)
    {
        BIG = 0.;
        for(c = 0; c < size_r; ++c)
        {
            if(std::abs(temp = LU[r + size_r*c]) > std::abs(BIG))
            {
                BIG = temp;
            }
        }
        if (BIG == 0.) 
        {
            printf("Singular matrix \n");
            exit(7);
        }
        vv[r] = 1./std::abs(BIG);
    }

    for(k = 0; k < size_r; ++k)
    {
        BIG = 0.;
        for (r = k; r < size_r; ++r)
        {
            temp = vv[r]*std::abs(LU[r + size_r*k]);
            if (std::abs(temp) > std::abs(BIG))
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
        if(std::abs(LU[k + k*size_r]) == 0)
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

template<class T> matrix<T> matrix<T>::solve(T b[])
{
    //This is just a wrapper for the other function, but allows it to be returned as a matrix. Can be useful, as in our main project.
    T x[size_c];
    matrix<T> return_val(size_r, 1);
    solve(x, b);
    for (int i = 0; i < size_r; ++i)
    {
        return_val.set_element(i, 0, x[i]);
    }
    return return_val;
}

template<class T> void matrix<T>::calc_det()
{
    if(LUcheck == false)
    {
        calc_LU();
    }
    det = d;
    for (int i = 0; i < size_r; ++i)
    {
        det *= LU[i + i*size_r];
    }
    det_check = true;
}

template<class T> T matrix<T>::get_det()
{
    if(!det_check)
    {
        calc_det();
    }
    return det;
}

template<class T> matrix<T> matrix<T>::forward_vector_multiplication(T vec[])
{
    matrix<T> return_val(1, size_r);
    T element;
    for (int i = 0; i < size_c; ++i)
    {
        element = 0;
        for (int k = 0; k < size_r; ++k)
        {
            element += get_element(i, k)*vec[k];
        }
        return_val.set_element(0, i, element);
    }
    return return_val;
}

template<class T> matrix<T> matrix<T>::backward_vector_multiplication(T vec[])
{
    matrix<T> return_val(size_r, 1);
    T element;
    for (int i = 0; i < size_c; ++i)
    {
        element = 0;
        for (int k = 0; k < size_r; ++k)
        {
            element += vec[k]*get_element(k, i);
        }
        return_val.set_element(0, i, element);
    }
    return return_val;
}

template<class T> matrix<T> matrix<T>::transpose()
{
    matrix<T> return_val(size_c, size_r);
    for(int r = 0; r < size_r; ++r)
    {
        for (int c = 0; c < size_c; ++c)
        {
            return_val.set_element(c, r, get_element(r, c));
        }
    }
    return return_val;
}

template<class T> void matrix<T>::resize(int new_r, int new_c)
{
    
    if (new_r >= size_r && new_c >= size_c)
    {
        T* temp = new T[size_r*size_c];
        for (int i = 0; i < size_r*size_c; ++i)
        {
            temp[i] = storage[i];
        }
    }

    delete[] storage;
    storage = new T[new_r*new_c];

    if (new_r >= size_r && new_c >= size_c)
    {
        for(int r = 0; r < new_r; ++r)
        {
            for (int c = 0; c < new_c; ++c)
            {
                storage[r + new_r*c] = temp[r + size_r*c];
            }
        }
        delete[] temp;
    }

    size_r = new_r;
    size_c = new_c;

    if(size_r == size_c)
    {
        square = true;
    }

    det_check = false;
    LU_check = false;
}

#endif //MATRIX_H_INCLUDED