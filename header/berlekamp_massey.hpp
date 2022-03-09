#ifndef BERLEKAMP_MASSEY_H
#define BERLEKAMP_MASSEY_H
template<typename field>
class berlekamp_massey
{
    private:
    field *space, *c, *t1, *t2;
    int allocated_size = 0;
    public:
    berlekamp_massey() {};
    void allocate(int size);
    void allocate(field* out_space, int size);
    field* operator()(field* s, int* deg, int size);
    ~berlekamp_massey ();
};

#include "berlekamp_massey.tpp"
#endif