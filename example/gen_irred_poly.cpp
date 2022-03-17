#include <random>
#include <iostream>
#include "../header/finite_field.hpp"
#include "../header/polynomial_ring.hpp"

const int max_deg = 2048;
finite_field<2> coeff[max_deg+1];
int main(void)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 1);
    int deg; polynomial_ring<max_deg, 2> poly;
    std::cin >> deg;
    do
    {
        std::cout << "tic" << '\n';
        for(int i=0; i<deg; i++)
        coeff[i] = dist(gen);
        coeff[deg] = 1;
        poly = polynomial_ring<max_deg, 2>(deg, coeff);
        std::cout << "toc" << '\n';
    }
    while(!is_prime(poly));
    std::cout << poly << '\n';
}

