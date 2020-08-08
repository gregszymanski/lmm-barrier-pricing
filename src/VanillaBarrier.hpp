#ifndef VanillaBarrier_hpp
#define VanillaBarrier_hpp


#include "Header.hpp"


/*
 * Abstract class VanillaBarrier 
 * 
 * Used to simulate one path of a Vanilla Barrier (CALL/PUT) (UP/DOWN) (IN/OUT)
 * We can use this abstract class to get several implementations
 * and compare the resultss
*/




template <class Generator>
class VanillaBarrier {
public:
    // All the parameters of the simulation
    struct params_type {
        double s0;
        double barrier;
        double strike;
        
        double T;
        
        double sigma;
        double r;
        
        double n; // Discretisation du B&S
    };

    // Enums type to describe the kind of vanilla Barrier
    enum Vanilla {CALL = false, PUT = true};
    enum Barrier {UP = true, DOWN = false};
    enum TypeBar {IN = true, OUT = false};
    
    // Constructor
    VanillaBarrier(Vanilla v, Barrier b, TypeBar tb, params_type p) {
        vanilla = v;
        barrier = b;
        typebarrier = tb;
        params = p;
        step = params.T / params.n;
        nd = normal_distribution<>(0,1);
    }
    
    // Simulation of a Vanilla Barrier payoff
    virtual double operator()(Generator &gen) {
        current = params.s0;
        flag = check_barrier(current);
        
        for(int i = 0; i<params.n; ++i) {
            update(gen);
        }
        
        return payoff(current);
    };
    
    // Return the value given by the closed formula
    double trueValue() const {
        
        return 0;
    }
    
protected:
    normal_distribution<> nd;
    
    virtual void update(Generator &gen) = 0;
    
    bool check_barrier(double value) {
        return (barrier == UP and value >= params.barrier)
            or (barrier == DOWN and value <= params.barrier);
    };
    
    // We return the payoff (and we check if the barrier is triggered)
    double payoff(double value) const {
        double diff = value - params.strike;
        
        if(flag == typebarrier) {
            if(vanilla == CALL and diff > 0)
                return diff;
            if(vanilla == PUT and diff < 0)
                return -diff;
        }
        
        return 0;
    };

    
    Vanilla vanilla;
    Barrier barrier;
    TypeBar typebarrier;
    
    params_type params;
    
    double step;
    double current;
    bool flag;
};



template <class Generator>
class BrownianBridgeVanillaBarrier : public VanillaBarrier<Generator> {
private:
    using typename VanillaBarrier<Generator>::Vanilla;
    using typename VanillaBarrier<Generator>::Barrier;
    using typename VanillaBarrier<Generator>::TypeBar;
    using typename VanillaBarrier<Generator>::params_type;
    
    using VanillaBarrier<Generator>::current;
    using VanillaBarrier<Generator>::flag;
    using VanillaBarrier<Generator>::step;
    using VanillaBarrier<Generator>::barrier;
    using VanillaBarrier<Generator>::params;
    using VanillaBarrier<Generator>::nd;
    
    using VanillaBarrier<Generator>::check_barrier;
    using VanillaBarrier<Generator>::payoff;
    
public:
    BrownianBridgeVanillaBarrier(Vanilla v,
                           Barrier b,
                           TypeBar tb,
                           params_type p,
                           bool antit = false):
    VanillaBarrier<Generator>(v, b, tb, p) {
        nu = params.r/params.sigma - params.sigma/2;
        unif = uniform_real_distribution<>(0,1);
        sgn = barrier ? 1 : -1;
        anti = antit;
    };
    
    // Simulation of a Vanilla Barrier payoff
    virtual double operator()(Generator &gen) {
        
        double tmp_normal = nd(gen);
        double tmp_unif = unif(gen);
        
        if(anti)
            return (internal_payoff(tmp_normal,tmp_unif) +
                    internal_payoff(-tmp_normal,1-tmp_unif)) / 2;
        else
            return internal_payoff(tmp_normal,tmp_unif);
    };
    
protected:
    void update(Generator &gen)  {
        // We do not need any update
    };
    
    double internal_payoff(double normal, double uniform) {
        double y = sqrt(params.T) * normal + nu * params.T;
        double z = (y + sgn * sqrt(y*y - 2 * params.T * log(uniform)))/2;
        
        y = params.s0 * exp(params.sigma * y);
        z = params.s0 * exp(params.sigma * z);
        
        flag = check_barrier(z);
        
        return payoff(y);
    };
    
    uniform_real_distribution<> unif;
    double nu;
    int sgn;
    
    bool anti;
};


template <class Generator>
class explicitVanillaBarrier : public VanillaBarrier<Generator> {
private:
    using typename VanillaBarrier<Generator>::Vanilla;
    using typename VanillaBarrier<Generator>::Barrier;
    using typename VanillaBarrier<Generator>::TypeBar;
    using typename VanillaBarrier<Generator>::params_type;
    
    using VanillaBarrier<Generator>::current;
    using VanillaBarrier<Generator>::flag;
    using VanillaBarrier<Generator>::step;
    using VanillaBarrier<Generator>::params;
    using VanillaBarrier<Generator>::nd;
    
    using VanillaBarrier<Generator>::check_barrier;
    
public:
    explicitVanillaBarrier(Vanilla v,
                           Barrier b,
                           TypeBar tb,
                           params_type p):
    VanillaBarrier<Generator>(v, b, tb, p) {};
    
    
protected:
    void update(Generator &gen)  {
        // We update the value of the process using usual Log normal formula
        double s = params.sigma;
        double tmp = step * (params.r - s*s/2)
            + sqrt(step) * s * nd(gen);
        current = current * exp(tmp);
        
        // We check if we touch the barrier
        flag = flag || check_barrier(current);
    };
    
    
};


template <class Generator>
class EulerVanillaBarrier : public VanillaBarrier<Generator> {
private:
    using typename VanillaBarrier<Generator>::Vanilla;
    using typename VanillaBarrier<Generator>::Barrier;
    using typename VanillaBarrier<Generator>::TypeBar;
    using typename VanillaBarrier<Generator>::params_type;
    
    using VanillaBarrier<Generator>::current;
    using VanillaBarrier<Generator>::flag;
    using VanillaBarrier<Generator>::step;
    using VanillaBarrier<Generator>::params;
    using VanillaBarrier<Generator>::nd;
    
    using VanillaBarrier<Generator>::check_barrier;
    
public:
    EulerVanillaBarrier(Vanilla v,
                        Barrier b,
                        TypeBar tb,
                        params_type p):
    VanillaBarrier<Generator>(v, b, tb, p) {};
    
    
protected:
    void update(Generator &gen)  {
        // We update the value of the process using the Euler scheme
        double drift_inc = params.r * step;
        double vol_inc = params.sigma * sqrt(step) * nd(gen);
        
        current = current * (1 + drift_inc + vol_inc);
        
        // We check if we touch the barrier
        flag = flag || check_barrier(current);
    };
    
    
};


template <class Generator>
class RWVanillaBarrier : public VanillaBarrier<Generator> {
private:
    using typename VanillaBarrier<Generator>::Vanilla;
    using typename VanillaBarrier<Generator>::Barrier;
    using typename VanillaBarrier<Generator>::TypeBar;
    using typename VanillaBarrier<Generator>::params_type;
    
    using VanillaBarrier<Generator>::current;
    using VanillaBarrier<Generator>::flag;
    using VanillaBarrier<Generator>::step;
    using VanillaBarrier<Generator>::params;
    using VanillaBarrier<Generator>::nd;
    
    using VanillaBarrier<Generator>::check_barrier;
    
public:
    RWVanillaBarrier(Vanilla v,
                        Barrier b,
                        TypeBar tb,
                        params_type p):
    VanillaBarrier<Generator>(v, b, tb, p) {
        unif = uniform_int_distribution<>(0,1);
    };
    
    
protected:
    void update(Generator &gen)  {
        // We update the value of the process using the Euler scheme
        double drift_inc = params.r * step;
        double vol_inc = params.sigma * sqrt(step) * (unif(gen)*2-1);
        
        current = current * (1 + drift_inc + vol_inc);
        
        // We check if we touch the barrier
        flag = flag || check_barrier(current);
    };
    
    uniform_int_distribution<> unif;
    
    
};




#endif /* VanillaBarrier_hpp */
