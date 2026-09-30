#pragma once
#include <iostream>
#include <vector>

namespace lib {

class Task {
public:
    virtual void load(std::istream& in) = 0;
    virtual double f(const std::vector<int>& x) = 0;
    virtual std::vector<int> neighbour(const std::vector<int>& x) = 0;
    virtual std::vector<int> random() = 0;
};

class Solver {
public:
    virtual std::vector<int> solve() = 0;
};

}
