#include "../header/finite_field.hpp"
#include "../header/polynomial_ring.hpp"
#include "../header/matrix.hpp"
#include <algorithm>
#include <iostream>

typedef unsigned long long int ull;

const int length = 126;

template<typename T>
class idx_pair
{
    private:
    ull idx;
    T value;
    public:
    idx_pair() {}
    idx_pair(ull idx_, T value_)
    : idx(idx_), value(value_)
    {}
    ull get_idx() {return idx;}
    T get_value() {return value;}
    bool operator<(const idx_pair & comp)
    {
        return (value < comp.value);
    }
};

bool search_idx(polynomial_ring<length, 2> searcher, ull* result, 
idx_pair<polynomial_ring<length, 2>>* space, ull size)
{
    ull lo = 0, hi = size-1, mid;
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

template<ull sub_grp_size, ull grp_size>
polynomial_ring<length, 2>cache_space(idx_pair<polynomial_ring<length, 2>>* space, 
polynomial_ring<length, 2> P, ull m)
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

template<ull sub_grp_size, ull grp_size>
ull baby_giant(polynomial_ring<length, 2> a, polynomial_ring<length, 2> & P,
polynomial_ring<length, 2> & inv_xm, idx_pair<polynomial_ring<length, 2>>* space, ull m)
{
    ull* idx = new ull;
    ull L = 0;
    ull max_giant = (sub_grp_size/m)+1;
    a = mod_pow(a, grp_size/sub_grp_size, P);
    for(int i=0; i<=max_giant; i++)
    {
        if(search_idx(a, idx, space, m))
        {L = i*m+(*idx); break;}
        a = (a*inv_xm) % P;
    }
    return L;
}

int int_to_poly(ull poly_int, finite_field<2>* poly)
{
    int tmp, deg=0;
    ull mask = 1ull<<63ull;
    for(int i=63; i>=0; i--)
    {
        if(poly_int>=mask)
        {
            poly_int ^= mask;
            poly[i] = 1; tmp = i;
            if(tmp>deg)
            deg = tmp;
        }
        else
        poly[i] = 0;
        mask >>= 1;
    }
    return deg;
}

ull xorshift64(ull x)
{
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    return x;
}

void xor_shift64_gen(ull x0, finite_field<2>* mat)
{
    for(int i=0; i<64*64; i += 64)
    {
        int_to_poly(x0, mat+i);
        x0 = xorshift64(x0);
    }
    return;
}

const ull grp_size = 18446744073709551615ull;
const ull sub_grp_1 = 3, sub_grp_2 = 5, sub_grp_3 = 17;
const ull sub_grp_4 = 257, sub_grp_5 = 641, sub_grp_6 = 65537;
const ull sub_grp_7 = 6700417;
const ull m_1 = 2, m_2 = 3, m_3 = 5;
const ull m_4 = 17, m_5 = 26, m_6 = 259;
const ull m_7 = 2589;
const ull f_1 = 2*(grp_size/sub_grp_1); 
const ull f_2 = 2*(grp_size/sub_grp_2);
const ull f_3 = 1*(grp_size/sub_grp_3);
const ull f_4 = 32*(grp_size/sub_grp_4);
const ull f_5 = 590*(grp_size/sub_grp_5);
const ull f_6 = 16384*(grp_size/sub_grp_6);
const ull f_7 = 3883315*(grp_size/sub_grp_7);
idx_pair<polynomial_ring<length, 2>> space_1[m_1], space_2[m_2], space_3[m_3];
idx_pair<polynomial_ring<length, 2>> space_4[m_4], space_5[m_5], space_6[m_6];
idx_pair<polynomial_ring<length, 2>> space_7[m_7];
finite_field<2> buffer[65], poly_mat[64*64], y[64];


int main(void)
{
    int T = 0, deg;
    ull l_1, l_2, l_3, l_4, l_5, l_6, l_7;
    ull x0, xf, result;
    polynomial_ring<length, 2> a;
    polynomial_ring<length, 2> P;
    polynomial_ring<length, 2> inv_xm_1, inv_xm_2, inv_xm_3;
    polynomial_ring<length, 2> inv_xm_4, inv_xm_5, inv_xm_6;
    polynomial_ring<length, 2> inv_xm_7;
    gauss_elimination<finite_field<2>> ge;

    std::ios_base::sync_with_stdio(false);
    std::cin.tie();

    for(int i=0; i<=64; i++)
    buffer[i] = 0;
    buffer[64] = 1; buffer[56] = 1; buffer[53] = 1; buffer[52] = 1; 
    buffer[53] = 1; buffer[52] = 1; buffer[51] = 1; buffer[50] = 1; 
    buffer[49] = 1; buffer[47] = 1; buffer[46] = 1; buffer[44] = 1; 
    buffer[42] = 1; buffer[39] = 1; buffer[37] = 1; buffer[33] = 1; 
    buffer[32] = 1; buffer[30] = 1; buffer[28] = 1; buffer[27] = 1; 
    buffer[23] = 1; buffer[20] = 1; buffer[16] = 1; buffer[13] = 1; 
    buffer[12] = 1; buffer[9] = 1; buffer[0] = 1;
    P = polynomial_ring<length, 2>(64, buffer);
    inv_xm_1 = cache_space<sub_grp_1, grp_size>(space_1, P, m_1);
    inv_xm_2 = cache_space<sub_grp_2, grp_size>(space_2, P, m_2);
    inv_xm_3 = cache_space<sub_grp_3, grp_size>(space_3, P, m_3);
    inv_xm_4 = cache_space<sub_grp_4, grp_size>(space_4, P, m_4);
    inv_xm_5 = cache_space<sub_grp_5, grp_size>(space_5, P, m_5);
    inv_xm_6 = cache_space<sub_grp_6, grp_size>(space_6, P, m_6);
    inv_xm_7 = cache_space<sub_grp_7, grp_size>(space_7, P, m_7);

    std::cin >> x0 >> xf;
    xor_shift64_gen(x0, poly_mat);
    int_to_poly(xf, y);
    ge.transpose(poly_mat, 64);
    ge.solve(poly_mat, y, 64);
    deg = 63;
    while(y[deg] == 0)
    deg--;
    if(deg < 0) deg = 0;
    a = polynomial_ring<length, 2>(deg, y);
    l_1 = baby_giant<sub_grp_1, grp_size>(a, P, inv_xm_1, space_1, m_1);
    l_2 = baby_giant<sub_grp_2, grp_size>(a, P, inv_xm_2, space_2, m_2);
    l_3 = baby_giant<sub_grp_3, grp_size>(a, P, inv_xm_3, space_3, m_3);
    l_4 = baby_giant<sub_grp_4, grp_size>(a, P, inv_xm_4, space_4, m_4);
    l_5 = baby_giant<sub_grp_5, grp_size>(a, P, inv_xm_5, space_5, m_5);
    l_6 = baby_giant<sub_grp_6, grp_size>(a, P, inv_xm_6, space_6, m_6);
    l_7 = baby_giant<sub_grp_7, grp_size>(a, P, inv_xm_7, space_7, m_7);
    result = ((__int128_t)l_1*f_1+(__int128_t)l_2*f_2+(__int128_t)l_3*f_3+
    (__int128_t)l_4*f_4+(__int128_t)l_5*f_5+(__int128_t)l_6*f_6+
    (__int128_t)l_7*f_7) % grp_size;
    std::cout << result+1 << '\n';
    return 0;
}