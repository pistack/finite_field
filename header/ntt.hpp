#ifndef NTT_H
#define NTT_H

#include "finite_field.hpp"

template<int k, int p>
class ntt
{
    private:
    finite_field<p> w = 0, w_inv = 0;


    public:
    constexpr ntt();
    void operator()(finite_field<p>* ary, 
    int k_idx, bool inverse) const;
};

template<int k>
class arb_ntt
{
    private:
    static constexpr ntt<k, 998244353> ntt_p1 = ntt<k, 998244353>();
    static constexpr ntt<k, 897581057> ntt_p2 = ntt<k, 897581057>();
    static constexpr ntt<k, 880803841> ntt_p3 = ntt<k, 880803841>();
    
    int pow(long long int base, int idx, int p);
    int wisdom(int size, int p);
    void fft(int* a, finite_field<998244353>* a1, finite_field<897581057>* a2, finite_field<880803841>* a3,
    int* b, finite_field<998244353>* b1, finite_field<897581057>* b2, finite_field<880803841>* b3,
    int* b_inv, int size, int k_idx, int p);

    public:
    arb_ntt() {};
    void operator()(int* ary, int size, int p, bool inverse);
};

#include "ntt.tpp"

#endif