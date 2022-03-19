#ifndef BERLEKAMP_MASSEY_H
#define BERLEKAMP_MASSEY_H
template<typename field>
field* berlekamp_massey(field s, int* deg, int size);

#include "berlekamp_massey.tpp"
#endif