#include "../header/finite_field.hpp"
#include "../header/polynomial_ring.hpp"
#include <algorithm>
#include <iostream>

template<typename T>
class idx_pair
{
    private:
    unsigned int idx;
    T value;
    public:
    idx_pair() {}
    idx_pair(unsigned int idx_, T value_)
    : idx(idx_), value(value_)
    {}
    unsigned int get_idx() {return idx;}
    T get_value() {return value;}
    bool operator<(const idx_pair & comp)
    {
        return (value < comp.value);
    }
};

bool search_idx(polynomial_ring<62, 2> searcher, unsigned int* result, 
idx_pair<polynomial_ring<62,2>>* space, unsigned int size)
{
    unsigned int lo = 0; unsigned int hi = size-1; unsigned int mid;
    while(hi - lo > 1)
    {
        mid = (hi+lo)/2;
        if(space[mid].get_value() < searcher)
        lo = mid;
        else
        hi = mid;
    }
    if(space[hi].get_value() == searcher)
    {
        *result = space[hi].get_idx();
        return true;
    }
    if(space[lo].get_value() == searcher)
    {
        *result = space[lo].get_idx();
        return true;
    }
    return false;
}

polynomial_ring<62,2>cache_space(idx_pair<polynomial_ring<62, 2>>* space, 
polynomial_ring<62, 2> P, unsigned int m)
{

    finite_field<2> c[2] = {finite_field<2>(0), finite_field<2>(1)};
    polynomial_ring<62, 2> x(1, c);
    polynomial_ring<62, 2> pow_x(0, c+1);
    polynomial_ring<62, 2> inv_xm;

    for(int i=0; i<m; i++)
    {
        space[i] = idx_pair<polynomial_ring<62, 2>>(i, pow_x);
        pow_x = x*pow_x; pow_x = pow_x % P;
    }
    inv_xm = mod_pow(pow_x, 4294967294ul, P);
    std::sort(space, space+m);
    return inv_xm;
}

unsigned int baby_giant(polynomial_ring<62, 2> a, polynomial_ring<62, 2> & P,
polynomial_ring<62, 2> & inv_xm, idx_pair<polynomial_ring<62, 2>>* space, unsigned int m)
{
    unsigned int* idx = new unsigned int;
    unsigned int L = 0;
    unsigned int max_giant = (4294967295ul/m)+1;
    for(int i=0; i<=max_giant; i++)
    {
        if(search_idx(a, idx, space, m))
        {L = i*m+(*idx); break;}
        a = (a*inv_xm) % P;
    }
    return L;
}

int poly_str_to_poly(char* poly_str, finite_field<2>* buffer)
{
    int i=0, tmp = 0, deg = 0;
    while(poly_str[i] != '\0')
    {
        if(poly_str[i] != 'x' && poly_str[i] != '^')
        {
            if(poly_str[i+1] == '+' || poly_str[i+1] == '\0')
            {
                tmp = poly_str[i] - '0';
                buffer[tmp] = 1;
                if(deg < tmp)
                deg = tmp;
                if(poly_str[i+1] == '\0')
                break;
                i += 2;
            }
            else
            {
                tmp = 10*(poly_str[i]-'0')+(poly_str[i+1]-'0');
                buffer[tmp] = 1;
                if(deg < tmp)
                deg = tmp;
                if(poly_str[i+2] == '\0')
                break;
                i += 3;
            }
        }
        i++;
    }
    return deg;
}

const int m = 1<<16;
polynomial_ring<62, 2> a[200];
idx_pair<polynomial_ring<62, 2>> space[m];
finite_field<2> buffer[33];
char poly_str[161];
unsigned int result[200];


int main(void)
{
    int T = 0, deg;
    polynomial_ring<62, 2> P, inv_xm;
    std::ios_base::sync_with_stdio(false);
    std::cin.tie();
    for(int i=0; i<=32; i++)
    buffer[i] = 0;
    buffer[32] = 1; buffer[26] = 1;
    buffer[15] = 1; buffer[7] = 1;
    buffer[0] = 1;
    P = polynomial_ring<62, 2>(32, buffer);
    inv_xm = cache_space(space, P, m);
    while(true)
    {
        //std::cin >> poly_str;
        //if(poly_str[0] == '0')
        //break;
        for(int i=0; i<=32; i++)
        buffer[i] = 0;
        deg = 31;
        buffer[31] = 1; buffer[25] = 1; buffer[14] = 1; buffer[6] = 1;
        //deg = poly_str_to_poly(poly_str, buffer);
        a[T] = polynomial_ring<62, 2>(deg, buffer);
        result[T] = baby_giant(a[T], P, inv_xm, space, m);
        T++;
        if(T>=200)
        break;
    }
    for(int i=0; i<T; i++)
    std::cout << result[i] << '\n';
    return 0;
}