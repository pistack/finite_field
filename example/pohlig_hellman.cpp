#include "../header/finite_field.hpp"
#include "../header/polynomial_ring.hpp"
#include <algorithm>
#include <iostream>

const int length = 62;

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

bool search_idx(polynomial_ring<length, 2> searcher, unsigned int* result, 
idx_pair<polynomial_ring<length, 2>>* space, unsigned int size)
{
    unsigned int lo = 0, hi = size-1, mid;
    while(hi - lo > 1)
    {
        mid = (hi+lo)>>1;
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

template<unsigned int sub_grp_size, unsigned int grp_size>
polynomial_ring<length, 2>cache_space(idx_pair<polynomial_ring<length, 2>>* space, 
polynomial_ring<length, 2> P, unsigned int m)
{

    finite_field<2> c[2] = {finite_field<2>(0), finite_field<2>(1)};
    polynomial_ring<length, 2> x(1, c);
    polynomial_ring<length, 2> pow_x(0, c+1);
    polynomial_ring<length, 2> inv_xm;
    
    x = mod_pow(x, grp_size/sub_grp_size, P);
    for(int i=0; i<m; i++)
    {
        space[i] = idx_pair<polynomial_ring<length, 2>>(i, pow_x);
        pow_x = x*pow_x; pow_x = pow_x % P;
    }
    inv_xm = mod_pow(pow_x, sub_grp_size-1, P);
    std::sort(space, space+m);
    return inv_xm;
}

template<unsigned int sub_grp_size, unsigned int grp_size>
unsigned int baby_giant(polynomial_ring<length, 2> a, polynomial_ring<length, 2> & P,
polynomial_ring<length, 2> & inv_xm, idx_pair<polynomial_ring<length, 2>>* space, unsigned int m)
{
    unsigned int* idx = new unsigned int;
    unsigned int L = 0;
    unsigned int max_giant = (sub_grp_size/m)+1;
    a = mod_pow(a, grp_size/sub_grp_size, P);
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

const unsigned int grp_size = 4294967295;
const unsigned int sub_grp_1 = 3, sub_grp_2 = 5, sub_grp_3 = 17, sub_grp_4 = 257, sub_grp_5 = 65537;
const unsigned int m_1 = 2, m_2 = 3, m_3 = 5, m_4 = 17, m_5 = 259;
const unsigned int f_1 = 1*(grp_size/sub_grp_1); 
const unsigned int f_2 = 4*(grp_size/sub_grp_2);
const unsigned int f_3 = 2*(grp_size/sub_grp_3);
const unsigned int f_4 = 64*(grp_size/sub_grp_4);
const unsigned int f_5 = 32768*(grp_size/sub_grp_5);
polynomial_ring<length, 2> a[200];
idx_pair<polynomial_ring<length, 2>> space_1[m_1], space_2[m_2], space_3[m_3], space_4[m_4], space_5[m_5];
finite_field<2> buffer[33];
char poly_str[161];
unsigned int result[200];


int main(void)
{
    int T = 0, deg;
    unsigned int l_1, l_2, l_3, l_4, l_5;
    unsigned long long int tmp;
    polynomial_ring<length, 2> P;
    polynomial_ring<length, 2> inv_xm_1, inv_xm_2, inv_xm_3, inv_xm_4, inv_xm_5;
    std::ios_base::sync_with_stdio(false);
    std::cin.tie();
    for(int i=0; i<=32; i++)
    buffer[i] = 0;
    buffer[32] = 1; buffer[26] = 1;
    buffer[15] = 1; buffer[7] = 1;
    buffer[0] = 1;
    P = polynomial_ring<length, 2>(32, buffer);
    inv_xm_1 = cache_space<sub_grp_1, grp_size>(space_1, P, m_1);
    inv_xm_2 = cache_space<sub_grp_2, grp_size>(space_2, P, m_2);
    inv_xm_3 = cache_space<sub_grp_3, grp_size>(space_3, P, m_3);
    inv_xm_4 = cache_space<sub_grp_4, grp_size>(space_4, P, m_4);
    inv_xm_5 = cache_space<sub_grp_5, grp_size>(space_5, P, m_5);
    while(true)
    {
        std::cin >> poly_str;
        if(poly_str[0] == '0')
        break;
        for(int i=0; i<=32; i++)
        buffer[i] = 0;
        deg = poly_str_to_poly(poly_str, buffer);
        a[T] = polynomial_ring<length, 2>(deg, buffer);
        T++;
    }
    for(int i=0; i<T; i++)
    {
        l_1 = baby_giant<sub_grp_1, grp_size>(a[i], P, inv_xm_1, space_1, m_1);
        l_2 = baby_giant<sub_grp_2, grp_size>(a[i], P, inv_xm_2, space_2, m_2);
        l_3 = baby_giant<sub_grp_3, grp_size>(a[i], P, inv_xm_3, space_3, m_3);
        l_4 = baby_giant<sub_grp_4, grp_size>(a[i], P, inv_xm_4, space_4, m_4);
        l_5 = baby_giant<sub_grp_5, grp_size>(a[i], P, inv_xm_5, space_5, m_5);
        result[i] = (1ull*l_1*f_1+1ull*l_2*f_2+1ull*l_3*f_3+1ull*l_4*f_4+1ull*l_5*f_5) % grp_size;
    }
    for(int i=0; i<T; i++)
    std::cout << result[i] << '\n';
    return 0;
}