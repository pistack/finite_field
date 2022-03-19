#ifndef POLY_H
#define POLY_H
#include "finite_field.hpp"
#include "ntt.hpp"

template<int max_deg, int p>
class polynomial_ring
{
    int deg = -1;
    finite_field<p>* coeff;

    public:
    polynomial_ring() {}
    polynomial_ring(const int & c);
    polynomial_ring(const finite_field<p> & c);
    polynomial_ring(int d, finite_field<p>* c);
    polynomial_ring(const polynomial_ring<max_deg, p> & poly);
    ~polynomial_ring();
    bool operator<(const polynomial_ring<max_deg, p> & poly) const;
    bool operator>(const polynomial_ring<max_deg, p> & poly) const;
    bool operator==(const int & c) const 
    {return (deg == 0) && (coeff[0] == c); }
    bool operator!=(const int & c) const 
    {return (deg > 0) || (coeff[0] != c);}
    bool operator==(const finite_field<p> & c) const
    {return (deg == 0) && (coeff[0] == c); }
    bool operator!=(const finite_field<p> & c) const
    {return (deg > 0) || (coeff[0] != c);}
    bool operator==(const polynomial_ring<max_deg, p> & poly) const;
    bool operator!=(const polynomial_ring<max_deg, p> & poly) const;
    polynomial_ring<max_deg, p> & operator=(const polynomial_ring<max_deg, p> & assign);
    const polynomial_ring<max_deg, p> operator+(const polynomial_ring<max_deg, p> & add) const;
    const polynomial_ring<max_deg, p> operator-(const polynomial_ring<max_deg, p> & sub) const; 
    const polynomial_ring<max_deg, p> operator*(const polynomial_ring<max_deg, p> & mult) const;
    const polynomial_ring<max_deg, p> operator/(const polynomial_ring<max_deg, p> & div) const;
    const polynomial_ring<max_deg, p> operator%(const polynomial_ring<max_deg, p> & mod) const;
    bool is_unit() const {return (deg == 0) && (coeff[0] != 0);}
    int get_deg() const {return deg;}
    finite_field<p>* to_ary(int* d) const;
    friend std::ostream & operator<<(std::ostream & os, 
    const polynomial_ring<max_deg, p> & poly)
    {
        for(int i=0; i<=poly.deg; i++)
        os << poly.coeff[i] << ' ';
        return os;
    }
    friend std::istream & operator>>(std::istream & is,
    polynomial_ring<max_deg, p> & poly)
    {
        if(poly.deg == -1)
        poly.coeff = new finite_field<p>[max_deg+1];
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