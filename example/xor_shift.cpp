#include <iostream>
#include "../header/finite_field.hpp"
#include "../header/berlekamp_massey.hpp"
#include "../header/kitamasa.hpp"

typedef long long int ll;
u_int64_t xorshift64(uint64_t x)
{
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    return x;
}

void gen(finite_field<2>* guess, int max_guess)
{
    u_int64_t tmp = 1;
    guess[0] = finite_field<2>(tmp % 2);
    for(int i=1; i<max_guess; i++)
    {
        tmp = xorshift64(tmp);
        guess[i] = finite_field<2>(tmp % 2);
    }
    return;
}

constexpr int max_guess = 300;
finite_field<2> guess[max_guess];

int main(void)
{
    int* deg = new int;
    berlekamp_massey<finite_field<2>> bm;
    kitamasa<finite_field<2>> kt;
    finite_field<2>* char_poly;
    gen(guess, max_guess);
    bm.allocate(900);
    kt.allocate(1200);
    char_poly = bm(guess, deg, max_guess);
    std::cout << *deg << '\n';
    for(int i=0; i<=(*deg); i++)
    std::cout << char_poly[i] << ' ';
    std::cout << '\n';
    for(int i=0; i<=(*deg); i++)
    {
        if(char_poly[i].get_val())
        std::cout << (*deg)-i << ' ';
    }
    std::cout << '\n';
    for(int i=0; i<max_guess; i++)
    std::cout << i << ' ' << guess[i] << ' ' << kt(i, guess, char_poly, *deg) << '\n';
    return 0;
}