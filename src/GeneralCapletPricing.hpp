#ifndef GeneralCapletPricing_hpp
#define GeneralCapletPricing_hpp

#include "Header.hpp"

template <class Generator>
class GeneralCapletPricing {
public:
    // All the parameters involved in the simulation
    struct params_type {
        function<double(double)> volatitly;
        function<double(double, double)> F;
        
        double L0;
        double barrier;
        
        double strike;
        
        int n;
        
        double T;
        
        double optimalF(double t, double lnl) {
            normal nd{};

            double inst_vol = volatitly(t);
            
            t = T-t;
            double sqrtt = sqrt(t);
            double v = inst_vol * sqrtt;
            double v2 = v*v;
            double logstrike = log(strike);
            double logbarrier = log(barrier);
            double rate = exp(lnl);
            
            double delta = 0;
            double sum = 0;
            
            delta = (lnl - logstrike + v2 / 2)/v;
            sum += cdf(nd, delta);
            
            delta = (lnl - logbarrier + v2 / 2)/v;
            sum -= cdf(nd, delta);
            
            delta = (2 * logbarrier - logstrike - lnl - v2 / 2)/v;
            sum += strike/barrier * cdf(nd, delta);
            
            delta = (logbarrier - lnl - v2 / 2)/v;
            sum -= strike/barrier * cdf(nd, delta);
            
            delta = (lnl - logbarrier + v2 / 2)/v;
            sum -= 2 * pdf(nd, delta) / v;
            sum += 2 * strike/barrier * pdf(nd, delta) / v;
            
            return - inst_vol * rate *  sum;

        }
    };
    
    
    
    // Result of a simulation
    struct result {
        double tau; // Stopping time
        double L; // Value of the Libor when the simulation stopped
        double Z; // Used to reduce variance
    };
    
    
    GeneralCapletPricing(params_type p) {
        params = p;
        
        step = params.T/params.n;
        sqrt_step = sqrt(step);
        
        log_barrier = log(params.barrier);

        unif = uniform_int_distribution<>(0,1);
        unif_real = uniform_real_distribution<>(0,1);
        rv_normal = normal_distribution<>(0,1);
    };
    
    // Simulate one instance of the caplet using LLM with antithetic method
    double antithetic_rw_simulation(Generator &gen) {
        bool flag_anti = true;
        bool flag_last = true;
        
        int i = 0;
        
        double current = log(params.L0);
        double currentZ = 0;
        
        double currentAntithetic = current;
        double currentAntitheticZ = 0;
        
        // The or is to ensure both finish
        while((flag_anti || flag_last) && i<params.n) {
            int y = unif(gen) * 2 - 1;
            double time = i*step;
            double vol = params.volatitly(time);
            double lambdash = - vol * vol * step / 2 + vol * sqrt_step;
            // Represent lambda*sqrt(h) because lambda never appears alone
            
            
            
            // Algorithm 2:
            // If we are close to the barrier,
            // we check whether we escape or not
            // Then we update the current values
            
            // We do this for the usual simulation
            // and then for the antithetic one
            
            if(flag_last) {
                if(current + lambdash >= log_barrier) {
                    double p = lambdash / (log_barrier - current + lambdash);
                    double z = unif_real(gen);
                    
                    if(z<p) {
                        flag_last = false;
                        last = {.tau = time, .L = params.barrier, .Z = currentZ};
                    }
                    
                    current = current - lambdash;
                }
                
                currentZ = currentZ + params.F(time, current)* sqrt_step * y;
                current = current - step * vol * vol / 2 + sqrt_step * vol * y;
            }
            
            
            if(flag_anti) {
                if(currentAntithetic + lambdash >= log_barrier) {
                    double p = lambdash
                    / (log_barrier - currentAntithetic + lambdash);
                    double z = unif_real(gen);
                    
                    if(z<p) {
                        flag_anti = false;
                        lastAntithetic = {.tau = time,
                            .L = params.barrier,
                            .Z = currentAntitheticZ};
                    }
                    
                    currentAntithetic = currentAntithetic - lambdash;
                }
                
                currentAntitheticZ = currentAntitheticZ
                - params.F(time, currentAntithetic) * sqrt_step * y;
                currentAntithetic = currentAntithetic
                - step * vol * vol / 2 - sqrt_step * vol * y;
                
            }
            
            ++i;
        }
        
        
        // We do the same simulations only on the remaining part !
        
        if(flag_last) {
            last = {.tau = params.T, .L = exp(current), .Z = currentZ};
        }
        if(flag_anti) {
            lastAntithetic = {.tau = params.T, .L = exp(currentAntithetic), .Z = currentAntitheticZ};
        }
        
        return payoff();
    }
    
    // Simulate one instance of the caplet using LLM without antithetic method
    double slow_rw_simulation(Generator &gen) {
        bool flag_last = true;
        
        int i = 0;
        
        double current = log(params.L0);
        double currentZ = 0;
        
        
        while(flag_last && i<params.n) {
            int y = unif(gen) * 2 - 1;
            double time = i*step;
            double vol = params.volatitly(time);
            double lambdash = - vol * vol * step / 2 + vol * sqrt_step;
            // Represent lambda*sqrt(h) because lambda never appears alone
            
            // Algorithm 2:
            // If we are close to the barrier, we check whether we escape or not
            
            if(current + lambdash >= log_barrier) {
                flag_last = false;
                last = {.tau = time, .L = params.barrier, .Z = currentZ};
                current = current - lambdash;
            }
            
            // We update the current values
            
            currentZ = currentZ + params.F(time, current)* sqrt_step * y;
            current = current - step * vol * vol / 2 + sqrt_step * vol * y;
            
            ++i;
        }
        
        
        
        
        // There is no remaining part since antithetic can't finish
        // before the classic part
        
        if(flag_last) {
            last = {.tau = params.T, .L = exp(current), .Z = currentZ};
        }
        
        return payoff_without_anti();
    }
    
    // Simulate one instance of the caplet using LLM with antithetic method
    double antithetic_slow_rw_simulation(Generator &gen) {
        bool flag_anti = true;
        bool flag_last = true;
        
        int i = 0;
        
        double current = log(params.L0);
        double currentZ = 0;
        
        double currentAntithetic = current;
        double currentAntitheticZ = 0;
        
        // The or is to ensure both finish
        while((flag_anti || flag_last) && i<params.n) {
            int y = unif(gen) * 2 - 1;
            double time = i*step;
            double vol = params.volatitly(time);
            double lambdash = - vol * vol * step / 2 + vol * sqrt_step;
            // Represent lambda*sqrt(h) because lambda never appears alone
            
            
            
            // Algorithm 2:
            // If we are close to the barrier,
            // we check whether we escape or not
            // Then we update the current values
            
            // We do this for the usual simulation
            // and then for the antithetic one
            
            if(flag_last) {
                if(current + lambdash >= log_barrier) {
                    flag_last = false;
                    last = {.tau = time, .L = params.barrier, .Z = currentZ};
                    current = current - lambdash;
                }
                
                currentZ = currentZ + params.F(time, current)* sqrt_step * y;
                current = current - step * vol * vol / 2 + sqrt_step * vol * y;
            }
            
            
            if(flag_anti) {
                if(currentAntithetic + lambdash >= log_barrier) {
                    flag_anti = false;
                    lastAntithetic = {.tau = time,
                        .L = params.barrier,
                        .Z = currentAntitheticZ};
                    currentAntithetic = currentAntithetic - lambdash;
                }
                
                currentAntitheticZ = currentAntitheticZ
                - params.F(time, currentAntithetic) * sqrt_step * y;
                currentAntithetic = currentAntithetic
                - step * vol * vol / 2 - sqrt_step * vol * y;
                
            }
            
            ++i;
        }
        
        
        // We do the same simulations only on the remaining part !
        
        if(flag_last) {
            last = {.tau = params.T, .L = exp(current), .Z = currentZ};
        }
        if(flag_anti) {
            lastAntithetic = {.tau = params.T, .L = exp(currentAntithetic), .Z = currentAntitheticZ};
        }
        
        return payoff();
    }
    
    // Simulate one instance of the caplet using LLM without antithetic method
    double rw_simulation(Generator &gen) {
        bool flag_last = true;
        
        int i = 0;
        
        double current = log(params.L0);
        double currentZ = 0;
        
        
        while(flag_last && i<params.n) {
            int y = unif(gen) * 2 - 1;
            double time = i*step;
            double vol = params.volatitly(time);
            double lambdash = - vol * vol * step / 2 + vol * sqrt_step;
            // Represent lambda*sqrt(h) because lambda never appears alone
            
            // Algorithm 2:
            // If we are close to the barrier, we check whether we escape or not
            
            if(current + lambdash >= log_barrier) {
                double p = lambdash / (log_barrier - current + lambdash);
                double z = unif_real(gen);
                
                if(z<p) {
                    flag_last = false;
                    last = {.tau = time, .L = params.barrier, .Z = currentZ};
                }
                
                current = current - lambdash;
            }
            
            // We update the current values
            
            currentZ = currentZ + params.F(time, current)* sqrt_step * y;
            current = current - step * vol * vol / 2 + sqrt_step * vol * y;
            
            ++i;
        }
        
        
        
        
        // There is no remaining part since antithetic can't finish
        // before the classic part
        
        if(flag_last) {
            last = {.tau = params.T, .L = exp(current), .Z = currentZ};
        }
        
        return payoff_without_anti();
    }
    
    // Simulate one instance of the caplet using LLM without antithetic method usign Brownian motions
    double normal_simulation(Generator &gen) {
        bool flag_last = true;
        
        int i = 0;
        
        double current = log(params.L0);
        double currentZ = 0;
        
        
        while(flag_last && i<params.n) {
            double y = rv_normal(gen);
            double time = i*step;
            double vol = params.volatitly(time);
            
            // We update the current values
            
            currentZ = currentZ + params.F(time, current)* sqrt_step * y;
            current = current - step * vol * vol / 2 + sqrt_step * vol * y;
            
            if(current >= log_barrier) {
                flag_last = false;
                last = {.tau = time, .L = params.barrier, .Z = currentZ};
            }
            
            ++i;
        }
        
        // There is no remaining part since antithetic can't finish
        // before the classic part
        
        if(flag_last) {
            last = {.tau = params.T, .L = exp(current), .Z = currentZ};
        }
        
        return payoff_without_anti();
    }
    
    // Simulate one instance of the caplet using LLM with antithetic method  usign Brownian motions
    double antithetic_normal_simulation(Generator &gen) {
        bool flag_anti = true;
        bool flag_last = true;
        
        int i = 0;
        
        double current = log(params.L0);
        double currentZ = 0;
        
        double currentAntithetic = current;
        double currentAntitheticZ = 0;
        
        // The or is to ensure both finish
        while((flag_anti || flag_last) && i<params.n) {
            double y = rv_normal(gen);
            double time = i*step;
            double vol = params.volatitly(time);
            // Represent lambda*sqrt(h) because lambda never appears alone
            
            
            
            // Algorithm 2:
            // If we are close to the barrier,
            // we check whether we escape or not
            // Then we update the current values
            
            // We do this for the usual simulation
            // and then for the antithetic one
            
            if(flag_last) {
                currentZ = currentZ + params.F(time, current)* sqrt_step * y;
                current = current - step * vol * vol / 2 + sqrt_step * vol * y;
                
                if(current >= log_barrier) {
                    flag_last = false;
                    last = {.tau = time, .L = params.barrier, .Z = currentZ};
                }
            }
            
            
            if(flag_anti) {
                
                currentAntitheticZ = currentAntitheticZ
                - params.F(time, currentAntithetic) * sqrt_step * y;
                currentAntithetic = currentAntithetic
                - step * vol * vol / 2 - sqrt_step * vol * y;
                
                if(currentAntithetic >= log_barrier) {
                    flag_anti = false;
                    lastAntithetic = {.tau = time, .L = params.barrier, .Z = currentAntitheticZ};
                }
                
            }
            
            ++i;
        }
        
        
        // We do the same simulations only on the remaining part !
        
        if(flag_last) {
            last = {.tau = params.T, .L = exp(current), .Z = currentZ};
        }
        if(flag_anti) {
            lastAntithetic = {.tau = params.T, .L = exp(currentAntithetic), .Z = currentAntitheticZ};
        }
        
        return payoff();
    }
    
    // Simulate one instance of the caplet using LLM without antithetic method usign Brownian motions
    double bridge_normal_simulation(Generator &gen) {
        bool flag_last = true;
        
        int i = 0;
        
        double current = log(params.L0);
        double currentZ = 0;
        
        
        while(flag_last && i<params.n) {
            double y = rv_normal(gen);
            double u = unif_real(gen);
            double time = i*step;
            double vol = params.volatitly(time);
            
            // We update the current values
            
            currentZ = currentZ + params.F(time, current)* sqrt_step * y;
            double new_current = current - step * vol * vol / 2 + sqrt_step * vol * y;
            
            // Maximum computed using the Gaussian bridge method
            double maximum = bridge(current, new_current, vol, u);
            
            if(maximum >= log_barrier) {
                flag_last = false;
                last = {.tau = time, .L = params.barrier, .Z = currentZ};
            }
            
            current = new_current;
            
            ++i;
        }
        
        // There is no remaining part since antithetic can't finish
        // before the classic part
        
        if(flag_last) {
            last = {.tau = params.T, .L = exp(current), .Z = currentZ};
        }
        
        return payoff_without_anti();
    }

    // Simulate one instance of the caplet using LLM without antithetic method usign Brownian motions
    double antithetic_bridge_normal_simulation(Generator &gen) {
        bool flag_anti = true;
        bool flag_last = true;
        
        int i = 0;
        
        double current = log(params.L0);
        double currentZ = 0;
        
        double currentAntithetic = current;
        double currentAntitheticZ = 0;
        
        // The or is to ensure both finish
        while((flag_anti || flag_last) && i<params.n) {
            double y = rv_normal(gen);
            double u = unif_real(gen);
            double time = i*step;
            double vol = params.volatitly(time);
            // Represent lambda*sqrt(h) because lambda never appears alone
            
            
            
            // Algorithm 2:
            // If we are close to the barrier,
            // we check whether we escape or not
            // Then we update the current values
            
            // We do this for the usual simulation
            // and then for the antithetic one
            
            if(flag_last) {
                currentZ = currentZ + params.F(time, current)* sqrt_step * y;
                double new_current = current - step * vol * vol / 2 + sqrt_step * vol * y;
                
                // Maximum computed using the Gaussian bridge method
                double maximum = bridge(current, new_current, vol, u);
                if(maximum >= log_barrier) {
                    flag_last = false;
                    last = {.tau = time, .L = params.barrier, .Z = currentZ};
                }
                
                current = new_current;
            }
            
            
            if(flag_anti) {
                
                currentAntitheticZ = currentAntitheticZ
                - params.F(time, currentAntithetic) * sqrt_step * y;
                double new_current = currentAntithetic
                - step * vol * vol / 2 - sqrt_step * vol * y;
                
                double maximum = bridge(currentAntithetic, new_current, vol, 1-u);
                if(maximum >= log_barrier) {
                    flag_anti = false;
                    lastAntithetic = {.tau = time, .L = params.barrier, .Z = currentAntitheticZ};
                }
                
                currentAntithetic = new_current;
            }
            
            ++i;
        }
        
        
        // We do the same simulations only on the remaining part !
        
        if(flag_last) {
            last = {.tau = params.T, .L = exp(current), .Z = currentZ};
        }
        if(flag_anti) {
            lastAntithetic = {.tau = params.T, .L = exp(currentAntithetic), .Z = currentAntitheticZ};
        }
        
        return payoff();
    }

    
    
    
    
    
    
    
    result getLast() const {
        return last;
    };
    
    result getAntithetic() const {
        return lastAntithetic;
    };
    
    
    double payoff() const
    {
        double last_swap = last.L - params.strike;
        if(last.tau < params.T or last_swap <= 0)
            last_swap = 0;
        
        double lastAntithetic_swap = lastAntithetic.L - params.strike;
        if(lastAntithetic.tau < params.T or lastAntithetic_swap <= 0)
            lastAntithetic_swap = 0;
        
        return (last_swap + lastAntithetic_swap + last.Z + lastAntithetic.Z) / 2;
    };
    
    double payoff_without_anti() const
    {
        double last_swap = last.L - params.strike;
        if(last.tau < params.T or last_swap <= 0)
            last_swap = 0;
        
        
        return (last_swap + last.Z);
    };
    
    

protected:
    uniform_int_distribution<> unif;
    uniform_real_distribution<> unif_real;
    normal_distribution<> rv_normal;
    
    params_type params;
    
    result last;
    result lastAntithetic;
    
    double step;
    double sqrt_step;
    
    double log_barrier;
    
    // Simulation of the max of the log process
    // using Brownian Bridge method
    double bridge(double x, double y, double vol, double u) {
        double maximum = (x - y) * (x - y);
        maximum -= 2 * step * vol * vol * log(u);
        maximum = sqrt(maximum);
        
        return (maximum + x + y) / 2;
    }
};



#endif /* GeneralCapletPricing_hpp */
