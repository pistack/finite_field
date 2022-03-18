#ifndef FINITE_FIELD_H
#define FINITE_FIELD_H
#include <iostream>

template<int p>
class finite_field
{
    private:
    int value;
    public:
    finite_field() {}
    constexpr finite_field(const int & val)
    : value(val) {}
    constexpr finite_field(const finite_field<p> & copy)
    : value(copy.value) {}
    constexpr int get_val() const {return value;}
    constexpr bool operator<(const finite_field<p> & comp)
    {return (this->value < comp.value);}
    constexpr bool operator>(const finite_field<p> & comp)
    {return (this->value > comp.value);}
    constexpr bool operator==(const int & comp)
    {return (value == comp);}
    constexpr bool operator==(const finite_field<p> & comp)
    {return (value == comp.value);}
    constexpr bool operator!=(const int & comp)
    {return (value != comp);}
    constexpr bool operator!=(const finite_field<p> & comp)
    {return (value != comp.value);}
    constexpr finite_field<p> & operator+=(const finite_field<p> & add)
    {this -> value += add.value; this -> value %= p; return *this;}
    constexpr finite_field<p> & operator-=(const finite_field<p> & sub)
    {this -> value += p-sub.value; this -> value %= p; return *this;}
    constexpr finite_field<p> & operator*=(const finite_field<p> & mult)
    {this -> value = (1ll*(this->value)*mult.value) % p; return *this;}
    constexpr finite_field<p> & operator*=(const int & mult)
    {
        this -> value = (1ll*(this->value)*mult) % p; 
        this -> value = (this -> value < 0 ? this -> value = p + this -> value : this -> value);
        return *this;
    }
    constexpr finite_field<p> & operator/=(const finite_field<p> & div);
    constexpr finite_field<p> & operator=(const int & assign)
    {this -> value = assign; return *this;}
    constexpr finite_field<p> & operator=(const finite_field<p> & assign)
    {if(this != &assign){value = assign.value;} return *this;}
    constexpr finite_field<p> operator+(const finite_field<p> & add) const
    {return finite_field<p>(*this) += add;}
    constexpr finite_field<p> operator-(const finite_field<p> & sub) const
    {return finite_field<p>(*this) -= sub;}
    constexpr finite_field<p> operator*(const int & mult) const
    {return finite_field<p>(*this) *= mult;}
    constexpr finite_field<p> operator*(const finite_field<p> & mult) const
    {return finite_field<p>(*this) *= mult;}
    constexpr finite_field<p> operator/(const finite_field<p> & div) const
    {return finite_field<p>(*this) /= div;}
    friend std::ostream & operator<<(std::ostream & os, 
    const finite_field<p> & f_p)
    {os << f_p.value; return os;}
    friend std::istream & operator>>(std::istream & is,
    finite_field<p> & f_p)
    {is >> f_p.value; return is;}
};

template<>
class finite_field<2>
{
    private:
    bool value;
    public:
    finite_field() {}
    constexpr finite_field(const int & val)
    : value(val) {}
    constexpr finite_field(const finite_field<2> & copy)
    : value(copy.value) {}
    constexpr bool get_val() const {return value;}
    constexpr bool operator<(const finite_field<2> & comp)
    {return comp.value && (this->value ^ comp.value);}
    constexpr bool operator>(const finite_field<2> & comp)
    {return !comp.value && (this->value ^ comp.value);}
    constexpr bool operator==(const int & comp)
    {return (value == comp);}
    constexpr bool operator==(const finite_field<2> & comp)
    {return !(this -> value ^ comp.value);}
    constexpr bool operator!=(const int & comp)
    {return (value != comp);}
    constexpr bool operator!=(const finite_field<2> & comp)
    {return this->value ^ comp.value;}
    constexpr finite_field<2> & operator+=(const finite_field<2> & add)
    {this -> value ^= add.value; return *this;}
    constexpr finite_field<2> & operator-=(const finite_field<2> & sub)
    {this -> value ^= sub.value; return *this;}
    constexpr finite_field<2> & operator*=(const int & mult)
    {this -> value = this -> value && (mult != 0); return *this;}
    constexpr finite_field<2> & operator*=(const finite_field<2> & mult)
    {this -> value = this -> value && mult.value; return *this;}
    constexpr finite_field<2> & operator/=(const finite_field<2> & div)
    {this -> value = this -> value && div.value; return *this;}
    constexpr finite_field<2> & operator=(const int & assign)
    {this -> value = assign; return *this;}
    constexpr finite_field<2> & operator=(const finite_field<2> & assign)
    {if(this != &assign){value = assign.value;} return *this;}
    constexpr finite_field<2> operator+(const finite_field<2> & add) const
    {return finite_field<2>(*this) += add;}
    constexpr finite_field<2> operator-(const finite_field<2> & sub) const
    {return finite_field<2>(*this) -= sub;}
    constexpr finite_field<2> operator*(const int & mult) const
    {return finite_field<2>(*this) *= mult;}
    constexpr finite_field<2> operator*(const finite_field<2> & mult) const
    {return finite_field<2>(*this) *= mult;}
    constexpr finite_field<2> operator/(const finite_field<2> & div) const
    {return finite_field<2>(*this) /= div;}
    friend std::ostream & operator<<(std::ostream & os, 
    const finite_field<2> & f_p)
    {os << f_p.value; return os;}
    friend std::istream & operator>>(std::istream & is,
    finite_field<2> & f_p)
    {is >> f_p.value; return is;}
};

template<typename mon>
constexpr mon pow(mon base, unsigned long long int idx);

template<typename ED>
ED gcd(ED a, ED b);

template<typename ED, typename... EDs>
ED gcd(ED a, ED b, EDs... cs);

template<int p>
constexpr finite_field<p> quad_non_res();

template<int p>
class tonelli_shanks
{
    private:
    static constexpr finite_field<p> non_residual = quad_non_res<p>();
    public:
    tonelli_shanks() {}
    constexpr finite_field<p> operator()(finite_field<p> x);
};

template<>
class tonelli_shanks<2>
{
    public:
    constexpr finite_field<2> operator()(const finite_field<2> & x)
    {return finite_field<2>(x);}
};

#include "finite_field.tpp"

#endif