#ifndef Monte_Carlo_hpp
#define Monte_Carlo_hpp

#include "Header.hpp"


class Monte_Carlo_Results {
public:
    Monte_Carlo_Results(int n=0, double sx = 0, double sx2 = 0);
    void add(double value);
    double mean() const;
    double var() const;
    double confidence(double alpha = 0.95) const;
    int nbSim() const;
    
    Monte_Carlo_Results& operator+=(Monte_Carlo_Results const&);
    friend Monte_Carlo_Results operator+(Monte_Carlo_Results const&, Monte_Carlo_Results const&);
    
    friend ostream & operator<<(ostream &out, Monte_Carlo_Results const&);
    
    
protected:
    double sum_x;
    double sum_x2;
    int n;
    static normal nd;
};

//
template <class RandomVariable>
Monte_Carlo_Results monte_carlo(RandomVariable &&rv, int n = 1) {
    double sum_x = 0;
    double sum_x2 = 0;
    double u = 0;
    
    for(int i = 0; i < n; ++i) {
        u = rv();
        sum_x += u;
        sum_x2 += u*u;
    }
    
    return {n, sum_x, sum_x2};
}


template <class RandomVariable>
Monte_Carlo_Results monte_carlo(RandomVariable &&rv, double prec, double alpha=0.95, int batch_size = 1000, int max_iter = 100000) {
    auto res = monte_carlo(rv, batch_size);
    for(int i=1; i<max_iter && res.confidence(alpha) > prec; ++i)
        res += monte_carlo(rv, batch_size);
    return res;
}



#endif /* Monte_Carlo_hpp */
