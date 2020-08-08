#include "Monte_Carlo.hpp"

normal Monte_Carlo_Results::nd = normal(0,1);

Monte_Carlo_Results::Monte_Carlo_Results(int n, double sx, double sx2) : n(n), sum_x(sx), sum_x2(sx2) {};

void Monte_Carlo_Results::add(double value) {
    sum_x  += value;
    sum_x2 += value * value;
    ++n;
}

double Monte_Carlo_Results::mean() const {
    return sum_x / n;
}

double Monte_Carlo_Results::var() const {
    return sum_x2 / n - pow(mean(), 2);
}

Monte_Carlo_Results& Monte_Carlo_Results::operator+=(Monte_Carlo_Results const&m) {
    n += m.n;
    sum_x += m.sum_x;
    sum_x2 += m.sum_x2;
    return *this;
}

Monte_Carlo_Results operator+(Monte_Carlo_Results const&m1, Monte_Carlo_Results const&m2) {
    return Monte_Carlo_Results(m1.n + m2.n,
                   m1.sum_x + m2.sum_x,
                   m1.sum_x2 + m2.sum_x2);
}

int Monte_Carlo_Results::nbSim() const {
    return n;
}

double Monte_Carlo_Results::confidence(double alpha) const {
    return 2 * quantile(nd, (1+alpha)/2) * sqrt(var()/n);
}

ostream & operator<<(ostream &out, Monte_Carlo_Results const&m) {
    out << "Simulations = " << m.n << endl;
    out << "Mean        = " << m.mean() << endl;
    out << "Variance    = " << m.var() << endl;
    out << "Confidence  = " << m.confidence() << endl;
    return out;
}



