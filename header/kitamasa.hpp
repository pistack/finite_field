#ifndef KITAMASA_H
#define KITAMASA_H
template<typename ring>
class kitamasa
{
    private:
    int allocated_space = 0;
    ring *space, *a, *b, *r;
    int mult(ring* x, ring* y, ring* mod, 
    int deg_x, int deg_y, int deg_mod);
    public:
    kitamasa() {};
    void allocate(int size);
    void allocate(ring* out_space, int size);
    ring operator()(unsigned long long int n, 
    ring* init_term, ring* char_poly, int deg);
    ~kitamasa();
};
#include "kitamasa.tpp"
#endif