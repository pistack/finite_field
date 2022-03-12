#ifndef POLY_H
#include "finite_field.hpp"
#define POLY_H

template<int max_deg, int p>
class polynomial_ring
{
    int deg;
    finite_field<p> coeff[max_deg+1];

    public:
    polynomial_ring() {};
    polynomial_ring(int d_, finite_field<p>* c_);
    polynomial_ring(const polynomial_ring<max_deg, p> & poly);
    bool operator==(const polynomial_ring<max_deg, p> & poly);
    polynomial_ring<max_deg, p> & operator=(const polynomial_ring<max_deg, p> & assign);
    const polynomial_ring<max_deg, p> operator+(const polynomial_ring<max_deg, p> & add) const;
    const polynomial_ring<max_deg, p> operator-(const polynomial_ring<max_deg, p> & mult) const; 
    const polynomial_ring<max_deg, p> operator*(const polynomial_ring<max_deg, p> & mult) const;
    const polynomial_ring<max_deg, p> operator/(const polynomial_ring<max_deg, p> & div) const;
    const polynomial_ring<max_deg, p> operator%(const polynomial_ring<max_deg, p> & mod) const;
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
        is >> poly.deg;
        for(int i=0; i<=poly.deg; i++)
        is >> poly.coeff[i];
        return is;
    }
};

template<int max_deg>
class polynomial_ring<max_deg, 2>
{
    int deg;
    finite_field<2> coeff[max_deg+1];

    public:
    polynomial_ring() {};
    polynomial_ring(int d_, finite_field<2>* c_);
    polynomial_ring(const polynomial_ring<max_deg, 2> & poly);
    bool operator==(const polynomial_ring<max_deg, 2> & poly);
    bool operator<(const polynomial_ring<max_deg, 2> & poly);
    polynomial_ring<max_deg, 2> & operator=(const polynomial_ring<max_deg, 2> & assign);
    const polynomial_ring<max_deg, 2> operator+(const polynomial_ring<max_deg, 2> & add) const;
    const polynomial_ring<max_deg, 2> operator-(const polynomial_ring<max_deg, 2> & mult) const; 
    const polynomial_ring<max_deg, 2> operator*(const polynomial_ring<max_deg, 2> & mult) const;
    const polynomial_ring<max_deg, 2> operator/(const polynomial_ring<max_deg, 2> & div) const;
    const polynomial_ring<max_deg, 2> operator%(const polynomial_ring<max_deg, 2> & mod) const;
    friend std::ostream & operator<<(std::ostream & os, 
    const polynomial_ring<max_deg, 2> & poly)
    {
        for(int i=0; i<=poly.deg; i++)
        os << poly.coeff[i] << ' ';
        return os;
    }
    friend std::istream & operator>>(std::istream & is,
    polynomial_ring<max_deg, 2> & poly)
    {
        is >> poly.deg;
        for(int i=0; i<=poly.deg; i++)
        is >> poly.coeff[i];
        return is;
    }
};

#include "polynomial_ring.tpp"

#endif