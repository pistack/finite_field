#ifndef KITAMASA_H
#define KITAMASA_H
template<typename ring>
class kitamasa
{
    private:
    int mult(ring* x, ring* y, ring* r, ring* mod, 
    int deg_x, int deg_y, int deg_mod);
    public:
    kitamasa() {};
    ring operator()(unsigned long long int n, 
    ring* init_term, ring* char_poly, int deg);
};
#include "kitamasa.tpp"
#endif