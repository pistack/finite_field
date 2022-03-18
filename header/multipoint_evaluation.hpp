#ifndef MULT_PTS_EV_H
#define MULT_PTS_EV_H

#include "ntt.hpp"

template<int k, int p>
class multipoint_evaluation
{
    private:
    static constexpr ntt<k+2, p> ntt_p = ntt<k+2, p>();
    int allocated_idx = -1;
    finite_field<p> *space;
    finite_field<p> *poly_tree;
    finite_field<p> *tmp, *tmp2, *tmp3;
    finite_field<p> *ext_f, *ext_g, *collect;
    int construct(finite_field<p>* x, int num_pts);
    void inverse(finite_field<p>* p1, finite_field<p>* p1_rev,
    finite_field<p>* p1_inv, finite_field<p>* p1_inv_pre,
    int deg_p1, int deg_inv);
    void remainder(finite_field<p>* f, finite_field<p>* g, 
    finite_field<p>* r, finite_field<p>* supp1, 
    finite_field<p>* supp2, int deg_f, int deg_g, int k_idx);
    public:
    multipoint_evaluation() {};
    void allocate(int k_idx);
    void allocate(finite_field<p>* out_space, int k_idx);
    finite_field<p>* operator()(finite_field<p>* f, int deg, 
    finite_field<p>* x, int num_pts);
    ~multipoint_evaluation();
};

#include "multipoint_evaluation.tpp"

#endif