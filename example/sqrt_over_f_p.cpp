#include <iostream>
#include "../header/finite_field.hpp"


constexpr tonelli_shanks<17> sqrt_17 = tonelli_shanks<17>();
constexpr finite_field<17> test = sqrt_17(finite_field<17>(8));

int main(void)
{
    std::cout << test << '\n';
    return 0;
}