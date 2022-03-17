#ifndef POLY_H
#include "finite_field.hpp"
#include "ntt.hpp"
#define POLY_H

template<int k, int p>
class polynomial_ring
{
    static ntt<k, 998244353> ntt_1;
    static ntt<k, 249561089> ntt_2;
    static ntt<k, 104857601> ntt_3; 
    int deg;
    finite_field<p> coeff[1<<k];

    public:
    polynomial_ring() {};
    polynomial_ring(int c) {deg = 0; coeff[0] = c % p;}
    polynomial_ring(int d, finite_field<p>* c);
    polynomial_ring(const polynomial_ring<k, p> & poly);
    bool operator<(const polynomial_ring<k, p> & poly);
    bool operator>(const polynomial_ring<k, p> & poly);
    bool operator==(int c)
    {return (deg == 0) && (coeff[0] == c); }
    bool operator!=(int c)
    {return (deg > 0) || (coeff[0] != c);}
    bool operator==(finite_field<p> c)
    {return (deg == 0) && (coeff[0] == c); }
    bool operator!=(finite_field<p> c)
    {return (deg > 0) || (coeff[0] != c);}
    bool operator==(const polynomial_ring<k, p> & poly);
    bool operator!=(const polynomial_ring<k, p> & poly);
    polynomial_ring<k, p> & operator=(const polynomial_ring<k, p> & assign);
    const polynomial_ring<k, p> operator+(const polynomial_ring<k, p> & add) const;
    const polynomial_ring<k, p> operator-(const polynomial_ring<k, p> & sub) const; 
    const polynomial_ring<k, p> operator*(const polynomial_ring<k, p> & mult) const;
    const polynomial_ring<k, p> operator/(const polynomial_ring<k, p> & div) const;
    const polynomial_ring<k, p> operator%(const polynomial_ring<k, p> & mod) const;
    int get_deg() const {return deg;}
    finite_field<p>* to_ary(int* d) const;
    friend std::ostream & operator<<(std::ostream & os, 
    const polynomial_ring<k, p> & poly)
    {
        for(int i=0; i<=poly.deg; i++)
        os << poly.coeff[i] << ' ';
        return os;
    }
    friend std::istream & operator>>(std::istream & is,
    polynomial_ring<k, p> & poly)
    {
        is >> poly.deg;
        for(int i=0; i<=poly.deg; i++)
        is >> poly.coeff[i];
        return is;
    }
};

template<int max_deg, int p>
polynomial_ring<max_deg, p> mod_pow(polynomial_ring<max_deg, p> base, 
unsigned long long int idx, polynomial_ring<max_deg, p> & mod);

template<int max_deg, int p>
bool is_prime(polynomial_ring<max_deg, p> poly);

#include "polynomial_ring.tpp"

#endif