#ifndef FINITE_FIELD_H
#define FINITE_FIELD_H
#include <random>
#include <iostream>

template<int p>
class finite_field
{
    private:
    int value;
    public:
    finite_field() {}
    finite_field(int val)
    : value(val) {}
    finite_field(const finite_field<p> & copy)
    : value(copy.value) {}
    int get_val() const {return value;}
    bool operator<(const finite_field<p> & comp)
    {
        int a = (this->value < 0 ? p+this->value : this->value);
        int b = (comp.value < 0 ? p+comp.value : comp.value);
        return a<b;
    }
    bool operator==(const int & comp)
    {return (value == comp);}
    bool operator==(const finite_field<p> & comp)
    {return ((this -> value - comp.value) % p == 0);}
    bool operator!=(const int & comp)
    {return (value != comp);}
    bool operator!=(const finite_field<p> & comp)
    {return ((this -> value - comp.value) % p != 0);}
    finite_field<p> & operator+=(const finite_field<p> & add)
    {this -> value += add.value; this -> value %= p; return *this;}
    finite_field<p> & operator-=(const finite_field<p> & sub)
    {this -> value -= sub.value; this -> value %= p; return *this;}
    finite_field<p> & operator*=(const finite_field<p> & mult)
    {this -> value = (1ll*(this->value)*mult.value) % p; return *this;}
    finite_field<p> & operator/=(const finite_field<p> & div);
    finite_field<p> & operator=(const int & assign)
    {this -> value = assign; return *this;}
    finite_field<p> & operator=(const finite_field<p> & assign)
    {if(this != &assign){value = assign.value;} return *this;}
    const finite_field<p> operator+(const finite_field<p> & add) const
    {return finite_field<p>(*this) += add;}
    const finite_field<p> operator-(const finite_field<p> & sub) const
    {return finite_field<p>(*this) -= sub;}
    const finite_field<p> operator*(const finite_field<p> & mult) const
    {return finite_field<p>(*this) *= mult;}
    const finite_field<p> operator/(const finite_field<p> & div) const
    {return finite_field<p>(*this) /= div;}
    friend std::ostream & operator<<(std::ostream & os, 
    const finite_field<p> & f_p)
    {os << (f_p.value<0 ? p+f_p.value : f_p.value); return os;}
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
    finite_field(int val)
    : value(val) {}
    finite_field(const finite_field<2> & copy)
    : value(copy.value) {}
    bool get_val() const {return value;}
    bool operator<(const finite_field<2> & comp)
    {return comp.value && (this->value ^ comp.value);}
    bool operator==(const int & comp)
    {return (value == comp);}
    bool operator==(const finite_field<2> & comp)
    {return !(this -> value ^ comp.value);}
    bool operator!=(const int & comp)
    {return (value != comp);}
    bool operator!=(const finite_field<2> & comp)
    {return this->value ^ comp.value;}
    finite_field<2> & operator+=(const finite_field<2> & add)
    {this -> value ^= add.value; return *this;}
    finite_field<2> & operator-=(const finite_field<2> & sub)
    {this -> value ^= sub.value; return *this;}
    finite_field<2> & operator*=(const finite_field<2> & mult)
    {this -> value = this -> value && mult.value; return *this;}
    finite_field<2> & operator/=(const finite_field<2> & div)
    {this -> value = this -> value && div.value; return *this;}
    finite_field<2> & operator=(const int & assign)
    {this -> value = assign; return *this;}
    finite_field<2> & operator=(const finite_field<2> & assign)
    {if(this != &assign){value = assign.value;} return *this;}
    const finite_field<2> operator+(const finite_field<2> & add) const
    {return finite_field<2>(*this) += add;}
    const finite_field<2> operator-(const finite_field<2> & sub) const
    {return finite_field<2>(*this) -= sub;}
    const finite_field<2> operator*(const finite_field<2> & mult) const
    {return finite_field<2>(*this) *= mult;}
    const finite_field<2> operator/(const finite_field<2> & div) const
    {return finite_field<2>(*this) /= div;}
    friend std::ostream & operator<<(std::ostream & os, 
    const finite_field<2> & f_p)
    {os << f_p.value; return os;}
    friend std::istream & operator>>(std::istream & is,
    finite_field<2> & f_p)
    {is >> f_p.value; return is;}
};

template<typename mon>
mon pow(mon base, unsigned long long int idx);

template<int p>
class tonelli_shanks
{
    private:
    std::random_device rd;
    std::mt19937 rng{rd()};
    std::uniform_int_distribution<int> dist = \
    std::uniform_int_distribution<int>(2, p-1);
    finite_field<p> non_residual;
    void gen_non_res();
    public:
    tonelli_shanks() {gen_non_res();}
    finite_field<p> operator()(finite_field<p> x);
};

template<>
class tonelli_shanks<2>
{
    public:
    finite_field<2> operator()(finite_field<2> x)
    {return x;}
};

#include "finite_field.tpp"

#endif