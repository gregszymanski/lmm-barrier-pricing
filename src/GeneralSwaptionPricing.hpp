#ifndef GeneralSwaptionPricing_h
#define GeneralSwaptionPricing_h

#include "Header.hpp"

// These tree funcitions extends classical mathematics functions
// to apply to VectorXd element wise
inline VectorXd exp(VectorXd tmp) {
    const long n = tmp.size();
    for(int i=0; i<n; ++i) {
        tmp(i) = exp(tmp(i));
    }
    return tmp;
};

inline VectorXd square(VectorXd tmp) {
    const long n = tmp.size();
    for(int i=0; i<n; ++i) {
        tmp(i) = tmp(i)*tmp(i);
    }
    return tmp;
};

inline VectorXd log(VectorXd tmp) {
    const long n = tmp.size();
    for(int i=0; i<n; ++i) {
        tmp(i) = log(tmp(i));
    }
    return tmp;
};

// Inner product of two VectorXd
inline double dot(VectorXd tmp1, VectorXd tmp2) {
    const long n = tmp1.size();
    double r = 0;
    for(int i=0; i<n; ++i) {
        r += tmp1(i) * tmp2(i);
    }
    return r;
};







template <class Generator>
class GeneralSwaptionPricing {
public:
    
    MatrixXd somme_y;
    int TOTAL;
    
    // All the parameters involved in the simulation
    struct params_type {
        int dimension;
        
        function<VectorXd(double)> volatitly;
        function<VectorXd(double, VectorXd)> F;
        
        VectorXd initial;
        double barrier;
        
        double strike;
        
        MatrixXd correlations;
        
        int n;
        
        double T0;
        double delta;
        
        
        // Can be used to improve performance of the simulation (to reduce variance)
        // The parameter useVege indicates whether we should differentiate
        // the Rebonato's approximation of the Swap variance with respect to
        // the LIBOR's rates.
        // By default, it is false because empirically we had better results
        // when we consider this variance as fixed and then we introduce this
        // approximation (after the differentiation)
        VectorXd optimalF_constantvol(double t, VectorXd log_libor, bool useVega = false) {
            normal nd{};

            double swap_rate;
            double log_swap_rate;
            double factor;
            VectorXd gradient_swap_rate = MatrixXd::Zero(dimension,1);
            VectorXd gradient_factor = MatrixXd::Zero(dimension,1);
            VectorXd gradient_rebonatoVolatility = MatrixXd::Zero(dimension,1);
            VectorXd omegas = MatrixXd::Zero(dimension, 1);
            
            VectorXd sums = MatrixXd::Zero(dimension,1);
            VectorXd products = MatrixXd::Zero(dimension,1);
            
            VectorXd libor = exp(log_libor);
            VectorXd constant_vol = volatitly(t);
            
            double remaining = T0 - t;
            
            double logstrike = log(strike);
            double logbarrier = log(barrier);

            double rebonatoVolatility = 0;
            double rebonatoVolatility2 = 0;
            double sum = 0;
            double product = 1;
            double delta_bs = 0;
            double delta_bs_op = 0;
            double rebonatoApprox = 0;
            double diffRebonatoApprox = 0;
            double vegaRebonatoApprox = 0;
            

            
            for(int j = 0; j < dimension; ++j) {
                product /= 1 + delta * libor(j);
                // Now, product contains
                // \Pi _ {l = 0}^{j} \frac{1}{1 + delta * libor(j)}
                
                sum += product;
                // Now, sum contains
                // \sum_{m = 0}^{j} \Pi _ {l = 0}^{m} \frac{1}{1 + delta * libor(j)}
                
                sums(j) = sum;
                products(j) = product;
            }
            
            //cout << "products " << endl << products << endl;
            //cout << "sums     " << endl << sums << endl;
            
            
            VectorXd rest_sums = MatrixXd::Zero(dimension,1);
            for(int j = 0; j < dimension; ++j) {
                rest_sums(j) = sum - sums(j);
            }

            double sum2 = sum * sum;

            
            for(int j = 0; j < dimension; ++j) {
                double tmp_sum = 0;
                if(j > 0)
                    tmp_sum = sums(j-1);
                
                gradient_swap_rate(j) = (tmp_sum * product + (sum - tmp_sum))
                / ( (1 + delta * libor(j)) * sum2);
                
                gradient_factor(j) = - delta * (sum - tmp_sum) / (1 + delta * libor(j));
                
                omegas(j) = products(j) / (delta * sum);
            }
            
            swap_rate = (1 - product) / (delta * sum);
            factor = sum;
            log_swap_rate = log(swap_rate);
            
            //cout << "swap_rate " << swap_rate << endl;
            

            // Calcul de la volatilité de rebenato
            for(int i = 0; i<dimension; ++i) {
                for(int j=0; j<dimension; ++j) {
                    rebonatoVolatility2 += omegas(i) * omegas(j) * libor(i) * libor(j) * constant_vol(i) * constant_vol(j) * correlations(i,j);
                }
            }
            rebonatoVolatility2 *= remaining/(swap_rate*swap_rate);
            rebonatoVolatility = sqrt(rebonatoVolatility2);
            
            
            
            // Calcul de l'approximation de rebenato
            delta_bs = (log_swap_rate - logstrike + rebonatoVolatility2 / 2)/rebonatoVolatility;
            rebonatoApprox += swap_rate * cdf(nd, delta_bs);
            
            delta_bs = (log_swap_rate - logbarrier + rebonatoVolatility2 / 2)/rebonatoVolatility;
            rebonatoApprox -= swap_rate * cdf(nd, delta_bs);
            
            delta_bs = (log_swap_rate - logstrike - rebonatoVolatility2 / 2)/rebonatoVolatility;
            rebonatoApprox -= strike * cdf(nd, delta_bs);
            
            delta_bs = (log_swap_rate - logbarrier - rebonatoVolatility2 / 2)/rebonatoVolatility;
            rebonatoApprox += strike * cdf(nd, delta_bs);
            
            
            delta_bs = (2* logbarrier - logstrike - log_swap_rate + rebonatoVolatility2 / 2)/rebonatoVolatility;
            rebonatoApprox -= barrier * cdf(nd, delta_bs);
            
            delta_bs = (logbarrier - log_swap_rate + rebonatoVolatility2 / 2)/rebonatoVolatility;
            rebonatoApprox += barrier * cdf(nd, delta_bs);
            
            delta_bs = (2* logbarrier - logstrike - log_swap_rate - rebonatoVolatility2 / 2)/rebonatoVolatility;
            rebonatoApprox += strike * swap_rate / barrier * cdf(nd, delta_bs);
            
            delta_bs = (logbarrier - log_swap_rate - rebonatoVolatility2 / 2)/rebonatoVolatility;
            rebonatoApprox -= strike * swap_rate / barrier * cdf(nd, delta_bs);
            
            
            // Vega de l'approximation de rebenato
            
            delta_bs    = (log_swap_rate - logstrike + rebonatoVolatility2 / 2)/rebonatoVolatility;
            delta_bs_op = (log_swap_rate - logstrike - rebonatoVolatility2 / 2)/rebonatoVolatility;
            vegaRebonatoApprox += swap_rate * pdf(nd, delta_bs) * delta_bs_op;
            
            delta_bs    = (log_swap_rate - logbarrier + rebonatoVolatility2 / 2)/rebonatoVolatility;
            delta_bs_op = (log_swap_rate - logbarrier - rebonatoVolatility2 / 2)/rebonatoVolatility;
            vegaRebonatoApprox -= swap_rate * pdf(nd, delta_bs) * delta_bs_op;
            
            delta_bs    = (log_swap_rate - logstrike - rebonatoVolatility2 / 2)/rebonatoVolatility;
            delta_bs_op = (log_swap_rate - logstrike + rebonatoVolatility2 / 2)/rebonatoVolatility;
            vegaRebonatoApprox -= strike * pdf(nd, delta_bs) * delta_bs_op;
            
            delta_bs    = (log_swap_rate - logbarrier - rebonatoVolatility2 / 2)/rebonatoVolatility;
            delta_bs_op = (log_swap_rate - logbarrier + rebonatoVolatility2 / 2)/rebonatoVolatility;
            vegaRebonatoApprox += strike * pdf(nd, delta_bs) * delta_bs_op;
            
            
            delta_bs    = (2* logbarrier - logstrike - log_swap_rate + rebonatoVolatility2 / 2)/rebonatoVolatility;
            delta_bs_op = (2* logbarrier - logstrike - log_swap_rate - rebonatoVolatility2 / 2)/rebonatoVolatility;
            vegaRebonatoApprox -= barrier * pdf(nd, delta_bs) * delta_bs_op;
            
            delta_bs    = (logbarrier - log_swap_rate + rebonatoVolatility2 / 2)/rebonatoVolatility;
            delta_bs_op = (logbarrier - log_swap_rate - rebonatoVolatility2 / 2)/rebonatoVolatility;
            vegaRebonatoApprox += barrier * pdf(nd, delta_bs) * delta_bs_op;
            
            delta_bs    = (2* logbarrier - logstrike - log_swap_rate - rebonatoVolatility2 / 2)/rebonatoVolatility;
            delta_bs_op = (2* logbarrier - logstrike - log_swap_rate + rebonatoVolatility2 / 2)/rebonatoVolatility;
            vegaRebonatoApprox += strike * swap_rate / barrier * pdf(nd, delta_bs) * delta_bs_op;
            
            delta_bs    = (logbarrier - log_swap_rate - rebonatoVolatility2 / 2)/rebonatoVolatility;
            delta_bs_op = (logbarrier - log_swap_rate + rebonatoVolatility2 / 2)/rebonatoVolatility;
            vegaRebonatoApprox -= strike * swap_rate / barrier * pdf(nd, delta_bs) * delta_bs_op;
            
            vegaRebonatoApprox = - vegaRebonatoApprox / rebonatoVolatility;
            

            // Calcul de la differentielle de l'approximation de rebonato
            delta_bs = (log_swap_rate - logstrike + rebonatoVolatility2 / 2)/rebonatoVolatility;
            diffRebonatoApprox += cdf(nd, delta_bs);
            
            
            delta_bs = (log_swap_rate - logbarrier + rebonatoVolatility2 / 2)/rebonatoVolatility;
            diffRebonatoApprox -= cdf(nd, delta_bs);
            
            
            delta_bs = (2 * logbarrier - logstrike - log_swap_rate - rebonatoVolatility2 / 2)/rebonatoVolatility;
            diffRebonatoApprox += strike/barrier * cdf(nd, delta_bs);
            
            
            delta_bs = (logbarrier - log_swap_rate - rebonatoVolatility2 / 2)/rebonatoVolatility;
            diffRebonatoApprox -= strike/barrier * cdf(nd, delta_bs);
            
            
            delta_bs = (log_swap_rate - logbarrier + rebonatoVolatility2 / 2)/rebonatoVolatility;
            diffRebonatoApprox -= 2 * pdf(nd, delta_bs) / rebonatoVolatility;
            diffRebonatoApprox += 2 * strike/barrier * pdf(nd, delta_bs) / rebonatoVolatility;

            
            // Calcul du gradient de la volatilité de Rebonato
            
            for(int k=0; k<dimension; ++k) {
                for(int i = 0; i < dimension; ++i) {
                    gradient_rebonatoVolatility(k) += 2 * omegas(i) * omegas(k) * correlations(i,k) * constant_vol(i) * constant_vol(k) * T0 * libor(i);
                }
            }

            gradient_rebonatoVolatility /= 2 * swap_rate * swap_rate * rebonatoVolatility;
            

            
            // Calcul du gradient
            
            VectorXd gradient_V_swaption = MatrixXd::Zero(dimension,1);
            /*for(int i = 0; i<dimension; ++i) {
                gradient_V_swaption(i) = gradient_factor(i) * rebonatoApprox + factor * gradient_swap_rate(i) * diffRebonatoApprox;
            }*/
            if(useVega)
                gradient_V_swaption = rebonatoApprox * gradient_factor + factor * (diffRebonatoApprox * gradient_swap_rate + vegaRebonatoApprox * gradient_rebonatoVolatility);
            else
                gradient_V_swaption = rebonatoApprox * gradient_factor + factor * (diffRebonatoApprox * gradient_swap_rate);
            

            // Calcul du résultat
            VectorXd result = MatrixXd::Zero(dimension,1);
            for(int i = 0; i<dimension; ++i) {
                result(i) = - constant_vol(i) * libor(i) * gradient_V_swaption(i);
            }
            
            
            return result;
        }
    };

    
    
    
    // Result of a simulation
    struct result {
        double tau; // Stopping time
        
        VectorXd L; // Value of the Libor when the simulation stopped
        double Z; // Used to reduce variance
    };
    
    /*
     double currentZ = 0;
     double currentZAntithetic = 0;

    
    currentZ = currentZ +  sqrt_step * dot(params.F(time, current), y);
    currentZAntithetic = currentZAntithetic - sqrt_step * dot(params.F(time, currentAntithetic), y);
*/
    
    
    
    // Constructor
    GeneralSwaptionPricing(params_type p) {
        params = p;
        
        step = params.T0/params.n;
        sqrt_step = sqrt(step);
        
        log_barrier = log(params.barrier);
        
        unif = uniform_int_distribution<>(0,1);
        unif_real = uniform_real_distribution<>(0,1);
        rv_normal = normal_distribution<>(0,1);
        
        
        
        LLT<MatrixXd> decomp(params.correlations);
        pseudoroot = decomp.matrixL(); // Used to simulate the correlated noise

        somme_y = MatrixXd::Zero(params.dimension, params.dimension);
        TOTAL = 0;
    };
    
    // Simulation with Gaussian noise and no antithetic
    double normal_simulation(Generator &gen) {
        bool flag_last = true;
        int i = 0;
        
        VectorXd current = log(params.initial);
        double currentZ = 0;
        
        while(flag_last && i<params.n) {
            VectorXd y = pseudoroot * vec_normal(gen);
            //somme_y += y * y.transpose();
            //TOTAL++;
            
            double time = i*step;
            VectorXd vol = params.volatitly(time);
            VectorXd vol2 = square(vol);
            
            VectorXd drift = VectorXd::Zero(params.dimension);
            for(int j = 0; j < params.dimension; ++j) {
                double coef = exp(current(j));
                coef = coef / (1 + params.delta * coef) * vol(j);
                for(int i = j; i < params.dimension; ++i) {
                    drift(i) += coef * params.correlations(i,j);
                }
            }
            
            for(int i = 0; i < params.dimension; ++i) {
                drift(i) *= step * vol(i) * params.delta;
            }
            
            VectorXd cpt = y;
            for(int i=0; i<params.dimension; ++i) {
                cpt(i) *= vol(i);
            }
            
            
            VectorXd update = current + drift - step/2 * vol2 + sqrt_step * cpt ;
            
            
            if(internalSwapRate(exp(update)) > params.barrier) {
                last = {.tau = time, .L = exp(current), .Z = currentZ};
                return payoff();
            }
            
            
            currentZ = currentZ +  sqrt_step * dot(params.F(time, current), y);
            current = update;
            
            i = i + 1;
        }
        
        if(flag_last) {
            last = {.tau = params.T0, .L = exp(current), .Z = currentZ};
        }
        
        return payoff();
    }

    // Simulation with Gaussian noise and antithetic
    double normal_antithetic_simulation(Generator &gen) {
        bool flag_last = true;
        bool flag_lastAntithetic = true;
        int i = 0;
        
        VectorXd current = log(params.initial);
        VectorXd currentAntithetic = current;
        
        double currentZ = 0;
        double currentZAntithetic = 0;
        
        while((flag_last or flag_lastAntithetic) && i<params.n) {
            VectorXd y = pseudoroot * vec_normal(gen);
            
            double time = i*step;
            VectorXd vol = params.volatitly(time);
            VectorXd vol2 = square(vol);
            
            VectorXd drift = VectorXd::Zero(params.dimension);
            VectorXd driftAntithetic = VectorXd::Zero(params.dimension);
            
            for(int j = 0; j < params.dimension; ++j) {
                double coef = exp(current(j));
                double coefAntithetic = exp(currentAntithetic(j));
                
                coef = coef / (1 + params.delta * coef) * vol(j);
                coefAntithetic = coefAntithetic / (1 + params.delta * coefAntithetic) * vol(j);

                for(int i = j; i < params.dimension; ++i) {
                    drift(i) += coef * params.correlations(i,j);
                    driftAntithetic(i) += coefAntithetic * params.correlations(i,j);
                }
            }
            
            for(int i = 0; i < params.dimension; ++i) {
                drift(i) *= step * vol(i) * params.delta;
                driftAntithetic(i) *= step * vol(i) * params.delta;
            }
            
            VectorXd cpt = y;
            for(int i=0; i<params.dimension; ++i) {
                cpt(i) *= vol(i);
            }
            
            
            VectorXd update = current + drift - step/2 * vol2 + sqrt_step * cpt ;
            VectorXd updateAntithetic = currentAntithetic + driftAntithetic - step/2 * vol2 - sqrt_step * cpt ;
            
            currentZ = currentZ +  sqrt_step * dot(params.F(time, current), y);
            currentZAntithetic = currentZAntithetic - sqrt_step * dot(params.F(time, currentAntithetic), y);

            if(flag_last) {
                if(internalSwapRate(exp(update)) > params.barrier) {
                    flag_last = false;
                    last = {.tau = time, .L = exp(current), .Z = currentZ};
                }
                current = update;
            }
            if(flag_lastAntithetic) {
                if(internalSwapRate(exp(updateAntithetic)) > params.barrier) {
                    flag_lastAntithetic = false;
                    lastAntithetic = {.tau = time, .L = exp(currentAntithetic), .Z = currentZAntithetic};
                }
                currentAntithetic = updateAntithetic;
            }
            
            
            i = i + 1;
        }
        
        if(flag_last) {
            last = {.tau = params.T0, .L = exp(current), .Z = currentZ};
        }

        if(flag_lastAntithetic) {
            lastAntithetic = {.tau = params.T0, .L = exp(currentAntithetic), .Z = currentZAntithetic};
        }

        return payoffAntithetic();
    }

    // Simulation with RW(order=1/2) noise and no antithetic
    double slow_rw_simulation(Generator &gen) {
        bool flag_last = true;
        int i = 0;
        
        VectorXd current = log(params.initial);
        double currentZ = 0;

        while(flag_last && i<params.n) {
            VectorXd y = pseudoroot * vec_rw(gen);
            somme_y += y * y.transpose();
            TOTAL++;
            
            double time = i*step;
            
            VectorXd vol = params.volatitly(time);
            VectorXd vol2 = square(vol);
            
            
            double volmax = vol.maxCoeff();
            double libormax = current.maxCoeff();
            double futuremax = libormax + volmax*volmax * step * params.dimension + volmax * sqrt_step *sqrt(params.dimension);
            
            if(futuremax > log_barrier) {
                VectorXd diff_max = vol;
                for(int i = 0; i<params.dimension; ++i)
                    diff_max(i) = diff_max(i) * (sqrt_step *sqrt(i+1) + volmax * (i+1) * step);
                VectorXd updatemax = current + diff_max -vol2/2*step;
                if(internalSwapRate(exp(updatemax)) > params.barrier)
                {
                    last = {.tau = time, .L = exp(current), .Z = currentZ};
                    return payoff();
                }
            }
            
            
            VectorXd drift = VectorXd::Zero(params.dimension);
            for(int j = 0; j < params.dimension; ++j) {
                double coef = exp(current(j));
                coef = coef / (1 + params.delta * coef) * vol(j);
                for(int i = j; i < params.dimension; ++i) {
                    drift(i) += coef * params.correlations(i,j);
                }
            }
            
            for(int i = 0; i < params.dimension; ++i) {
                drift(i) *= step * vol(i) * params.delta;
            }
            
            VectorXd cpt = y;
            for(int i=0; i<params.dimension; ++i) {
                cpt(i) *= vol(i);
            }
            
            
            VectorXd update = current + drift - step/2 * vol2 + sqrt_step * cpt ;
            currentZ = currentZ +  sqrt_step * dot(params.F(time, current), y);

            
            
            current = update;
            
            i = i + 1;
        }
        
        if(flag_last) {
            last = {.tau = params.T0, .L = exp(current), .Z = currentZ};
        }
        
        return payoff();
    }
    
    // Simulation with RW(order=1/2) noise and antithetic
    double slow_rw_antithetic_simulation(Generator &gen) {
        bool flag_last = true;
        bool flag_anti = true;
        int i = 0;
        
        VectorXd current = log(params.initial);
        VectorXd currentAntithetic = current;

        double currentZ = 0;
        double currentZAntithetic = 0;

        while((flag_anti || flag_last) && i<params.n) {
            VectorXd y = pseudoroot * vec_rw(gen);
            somme_y += y * y.transpose();
            TOTAL++;
            
            double time = i*step;
            
            VectorXd vol = params.volatitly(time);
            VectorXd vol2 = square(vol);
            
            
            double volmax = vol.maxCoeff();
            
            VectorXd cpt = y;
            for(int i=0; i<params.dimension; ++i) {
                cpt(i) *= vol(i);
            }
            
            if(flag_last) {
                double libormax = current.maxCoeff();
                double futuremax = libormax + volmax*volmax * step * params.dimension + volmax * sqrt_step *sqrt(params.dimension);
                if(futuremax > log_barrier) {
                    VectorXd diff_max = vol;
                    for(int i = 0; i<params.dimension; ++i)
                        diff_max(i) = diff_max(i) * (sqrt_step *sqrt(i+1) + volmax * (i+1) * step);
                    VectorXd updatemax = current + diff_max -vol2/2*step;
                    if(internalSwapRate(exp(updatemax)) > params.barrier) {
                        flag_last = false;
                        last = {.tau = time, .L = exp(current), .Z = currentZ};
                    }
                }
                
                VectorXd drift = VectorXd::Zero(params.dimension);
                for(int j = 0; j < params.dimension; ++j) {
                    double coef = exp(current(j));
                    coef = coef / (1 + params.delta * coef) * vol(j);
                    for(int i = j; i < params.dimension; ++i) {
                        drift(i) += coef * params.correlations(i,j);
                    }
                }
                
                for(int i = 0; i < params.dimension; ++i) {
                    drift(i) *= step * vol(i) * params.delta;
                }
                
                currentZ = currentZ +  sqrt_step * dot(params.F(time, current), y);

                current = current + drift - step/2 * vol2 + sqrt_step * cpt ;
            }

            if(flag_anti) {
                double libormaxAntithetic = currentAntithetic.maxCoeff();
                double futuremaxAntithetic = libormaxAntithetic + volmax*volmax * step * params.dimension + volmax * sqrt_step *sqrt(params.dimension);
                if(futuremaxAntithetic > log_barrier) {
                    VectorXd diff_maxAntithetic = vol;
                    for(int i = 0; i<params.dimension; ++i)
                        diff_maxAntithetic(i) = diff_maxAntithetic(i) * (sqrt_step *sqrt(i+1) + volmax * (i+1) * step);
                    VectorXd updatemaxAntithetic = currentAntithetic + diff_maxAntithetic -vol2/2*step;
                    if(internalSwapRate(exp(updatemaxAntithetic)) > params.barrier) {
                        flag_anti = false;
                        lastAntithetic = {.tau = time, .L = exp(currentAntithetic), .Z = currentZAntithetic};
                    }
                }
                
                VectorXd driftAntithetic = VectorXd::Zero(params.dimension);
                for(int j = 0; j < params.dimension; ++j) {
                    double coefAntithetic = exp(currentAntithetic(j));
                    coefAntithetic = coefAntithetic / (1 + params.delta * coefAntithetic) * vol(j);
                    for(int i = j; i < params.dimension; ++i) {
                        driftAntithetic(i) += coefAntithetic * params.correlations(i,j);
                    }
                }
                
                for(int i = 0; i < params.dimension; ++i) {
                    driftAntithetic(i) *= step * vol(i) * params.delta;
                }
                
                currentZAntithetic = currentZAntithetic - sqrt_step * dot(params.F(time, currentAntithetic), y);
                currentAntithetic = currentAntithetic + driftAntithetic - step/2 * vol2 - sqrt_step * cpt ;
            }

            i = i + 1;
        }
        
        if(flag_last) {
            last = {.tau = params.T0, .L = exp(current), .Z = currentZ};
        }

        if(flag_anti) {
            lastAntithetic = {.tau = params.T0, .L = exp(currentAntithetic), .Z = currentZAntithetic};
        }

        return payoffAntithetic();
    }
    
    // Simulation with RW(order=1) noise and no antithetic
    double rw_simulation(Generator &gen) {
        bool flag_last = true;
        int i = 0;
        
        VectorXd current = log(params.initial);
        double currentZ = 0;

        while(flag_last && i<params.n) {
            VectorXd y = pseudoroot * vec_rw(gen);
            
            double time = i*step;
            
            VectorXd vol = params.volatitly(time);
            VectorXd vol2 = square(vol);
            
            
            double volmax = vol.maxCoeff();
            double libormax = current.maxCoeff();
            double futuremax = libormax + volmax*volmax * step * params.dimension + volmax * sqrt_step *sqrt(params.dimension);
            
            if(futuremax > log_barrier) {
                VectorXd diff_max = vol;
                for(int i = 0; i<params.dimension; ++i)
                    diff_max(i) = diff_max(i) * (sqrt_step *sqrt(i+1) + volmax * (i+1) * step);
                VectorXd updatemax = current + diff_max - vol2/2*step;
                if(internalSwapRate(exp(updatemax)) > params.barrier) {
                    VectorXd projected = projection(current);
                    
                    double lambdash = sqrt(params.dimension) * (volmax*volmax * step * params.dimension + volmax * sqrt_step *sqrt(params.dimension));
                    double p = lambdash / ((current - projected).norm() + lambdash) ;
                    double z = unif_real(gen);

                    
                    if(z<p) {
                        flag_last = false;
                        last = {.tau = time, .L = exp(current), .Z = currentZ};
                    }
                    
                    current = current + lambdash * (current - projected) / (current - projected).norm();
                    

                }
            }
            
            
            VectorXd drift = VectorXd::Zero(params.dimension);
            for(int j = 0; j < params.dimension; ++j) {
                double coef = exp(current(j));
                coef = coef / (1 + params.delta * coef) * vol(j);
                for(int i = j; i < params.dimension; ++i) {
                    drift(i) += coef * params.correlations(i,j);
                }
            }
            
            for(int i = 0; i < params.dimension; ++i) {
                drift(i) *= step * vol(i) * params.delta;
            }
            
            VectorXd cpt = y;
            for(int i=0; i<params.dimension; ++i) {
                cpt(i) *= vol(i);
            }
            
            
            VectorXd update = current + drift - step/2 * vol2 + sqrt_step * cpt ;
            currentZ = currentZ + sqrt_step * dot(params.F(time, current), y);
            
            current = update;
            
            i = i + 1;
        }
        
        if(flag_last) {
            last = {.tau = params.T0, .L = exp(current), .Z = currentZ};
        }
        
        return payoff();
    }

    // Simulation with RW(order=1) noise and antithetic
    double rw_antithetic_simulation(Generator &gen) {
        bool flag_last = true;
        bool flag_anti = true;
        int i = 0;
        
        VectorXd current = log(params.initial);
        VectorXd currentAntithetic = current;
        
        double currentZ = 0;
        double currentZAntithetic = 0;

        while((flag_anti || flag_last) && i<params.n) {
            VectorXd y = pseudoroot * vec_rw(gen);
            
            double time = i*step;
            
            VectorXd vol = params.volatitly(time);
            VectorXd vol2 = square(vol);
            double volmax = vol.maxCoeff();
            double lambdash = sqrt(params.dimension) * (volmax*volmax * step * params.dimension + volmax * sqrt_step *sqrt(params.dimension));

            VectorXd cpt = y;
            for(int i=0; i<params.dimension; ++i) {
                cpt(i) *= vol(i);
            }
            
            if(flag_last) {
                double libormax = current.maxCoeff();
                double futuremax = libormax + volmax*volmax * step * params.dimension + volmax * sqrt_step *sqrt(params.dimension);
                if(futuremax > log_barrier) {
                    VectorXd diff_max = vol;
                    for(int i = 0; i<params.dimension; ++i)
                        diff_max(i) = diff_max(i) * (sqrt_step *sqrt(i+1) + volmax * (i+1) * step);
                    VectorXd updatemax = current + diff_max -vol2/2*step;
                    if(internalSwapRate(exp(updatemax)) > params.barrier) {
                        VectorXd projected = projection(current);
                        double p = lambdash / ((current - projected).norm() + lambdash) ;
                        double z = unif_real(gen);
                        
                        
                        if(z<p) {
                            flag_last = false;
                            last = {.tau = time, .L = exp(current), .Z = currentZ};
                        }
                        
                        current = current + lambdash * (current - projected) / (current - projected).norm();

                        }
                }
                
                VectorXd drift = VectorXd::Zero(params.dimension);
                for(int j = 0; j < params.dimension; ++j) {
                    double coef = exp(current(j));
                    coef = coef / (1 + params.delta * coef) * vol(j);
                    for(int i = j; i < params.dimension; ++i) {
                        drift(i) += coef * params.correlations(i,j);
                    }
                }
                
                for(int i = 0; i < params.dimension; ++i) {
                    drift(i) *= step * vol(i) * params.delta;
                }
                
                currentZ = currentZ +  sqrt_step * dot(params.F(time, current), y);

                current = current + drift - step/2 * vol2 + sqrt_step * cpt ;
            }
            
            if(flag_anti) {
                double libormaxAntithetic = currentAntithetic.maxCoeff();
                double futuremaxAntithetic = libormaxAntithetic + volmax*volmax * step * params.dimension + volmax * sqrt_step *sqrt(params.dimension);
                if(futuremaxAntithetic > log_barrier) {
                    VectorXd diff_maxAntithetic = vol;
                    for(int i = 0; i<params.dimension; ++i)
                        diff_maxAntithetic(i) = diff_maxAntithetic(i) * (sqrt_step *sqrt(i+1) + volmax * (i+1) * step);
                    VectorXd updatemaxAntithetic = currentAntithetic + diff_maxAntithetic -vol2/2*step;
                    if(internalSwapRate(exp(updatemaxAntithetic)) > params.barrier) {
                        VectorXd projected = projection(currentAntithetic);
                        double p = lambdash / ((currentAntithetic - projected).norm() + lambdash) ;
                        double z = unif_real(gen);
                        
                        
                        if(z<p) {
                            flag_anti = false;
                            lastAntithetic = {.tau = time, .L = exp(currentAntithetic), .Z = currentZAntithetic};;
                        }
                        
                        currentAntithetic = currentAntithetic + lambdash * (currentAntithetic - projected) / (currentAntithetic - projected).norm();
                        
                    }

                }
                
                VectorXd driftAntithetic = VectorXd::Zero(params.dimension);
                for(int j = 0; j < params.dimension; ++j) {
                    double coefAntithetic = exp(currentAntithetic(j));
                    coefAntithetic = coefAntithetic / (1 + params.delta * coefAntithetic) * vol(j);
                    for(int i = j; i < params.dimension; ++i) {
                        driftAntithetic(i) += coefAntithetic * params.correlations(i,j);
                    }
                }
                
                for(int i = 0; i < params.dimension; ++i) {
                    driftAntithetic(i) *= step * vol(i) * params.delta;
                }
                
                currentZAntithetic = currentZAntithetic - sqrt_step * dot(params.F(time, currentAntithetic), y);
                currentAntithetic = currentAntithetic + driftAntithetic - step/2 * vol2 - sqrt_step * cpt ;
            }
            
            i = i + 1;
        }
        
        if(flag_last) {
            last = {.tau = params.T0, .L = exp(current), .Z = currentZ};
        }
        
        if(flag_anti) {
            lastAntithetic = {.tau = params.T0, .L = exp(currentAntithetic), .Z = currentZAntithetic};
        }
        
        return payoffAntithetic();
    }
    

    
    
    double payoff() {
        if(last.tau < params.T0)
            return last.Z;
        
        double sum = 0;
        double product = 1;
        for(int j = 0; j < params.dimension; ++j) {
            product /= 1 + params.delta * last.L(j);
            sum += product;
        }
        
        double swap_rate = (1 - product) / (params.delta * sum);
        
        double tmp = swap_rate - params.strike;
        if(tmp <= 0)
            return last.Z;
        
        return tmp * sum + last.Z;
    }

    double payoffAntithetic() {
        double tmp = payoff();
        last = lastAntithetic;
        return (tmp + payoff())/2;
    }

protected:
    uniform_int_distribution<> unif;
    uniform_real_distribution<> unif_real;
    normal_distribution<> rv_normal;

    
    VectorXd vec_normal(Generator &gen) {
        VectorXd r = VectorXd::Zero(params.dimension);
        for(int i = 0; i<params.dimension; ++i) {
            r(i) = rv_normal(gen);
        }
        return r;
    }

    VectorXd vec_rw(Generator &gen) {
        VectorXd r = VectorXd::Zero(params.dimension);
        for(int i = 0; i<params.dimension; ++i) {
            r(i) = unif(gen) * 2 - 1;
        }
        return r;
    }

    
    
    
    // Result of the computation of the first term in the least square minimisation
    struct evaluation {
        double exp_value;
        double value;
        VectorXd gradient;
    };
    
    // We do not use guess(0) which is then supposed to be equal to psi(guess) to be on the border! (voir pdf)
    evaluation psi(VectorXd guess) {
        
        VectorXd exp_guess = exp(guess);
        VectorXd products = MatrixXd::Constant(params.dimension, 1, 1.0);
        VectorXd sums = MatrixXd::Zero(params.dimension, 1);
        
        evaluation res;
        
        products(params.dimension - 1) = 1 + params.delta * exp_guess(params.dimension - 1);
        for(int i = params.dimension - 2; i >= 0; --i) {
            products(i) = products(i+1) * (1 + params.delta * exp_guess(i));
            sums(i) = sums(i+1) + products(i+1);
        }
        
        
        res.gradient = VectorXd(params.dimension);
        for(int i(1); i<params.dimension; ++i) {
            res.gradient(i) = exp_guess(i);
            res.gradient(i) = sums(i);
            res.gradient(i) = products(1);
            res.gradient(i) = - exp_guess(i) * ( params.delta * params.barrier * (1 + sums(i)) + 1 )/ ((1 + params.delta * exp_guess(i)) * products(1));
        }
        
        /*
        cout << "Barriere: " << params.barrier << endl;
        cout << "Sommes: " << sums << endl;
        cout << "Produits: " << products << endl;
        cout << "exp(psi(x)): " << (params.barrier * params.delta * (1+sums(0))+1)/(products(1)* params.delta) - 1/params.delta << endl;
        cout << "gradient: " << res.gradient << endl;
        */
        /*if(products(1) == 0) {
            cout << "PRODUCTS 1 == 0!!!!!!!" << endl;
            cout << guess << endl;
            
        }*/
        
        res.exp_value = (params.barrier * params.delta * (1+sums(0))+1)/(products(1)* params.delta) - 1/params.delta;
        res.value = log(res.exp_value);
        
        //cout << res.exp_value << endl;
        
        return res;
    }

    
    // Projection based on the gradient descent algorithm
    // Stable while "VectorXd current"'s elements are quite close
    VectorXd projection(VectorXd current) {
        VectorXd guess = current;
        VectorXd previous = current;
        
        double precision = 0.01;
        double gamma = 0.01;
        int i = 0;
        int max_iter = 1000;
        evaluation e;
        //double spread = precision / 10;
        
        VectorXd tmp_g;
        
        double lastprevious = psi(previous).exp_value;
        
        do {
            
            ++i;

            e = psi(guess);

            
            guess(0) = e.value;
            
            
            tmp_g = 2* (e.gradient / e.exp_value * (e.value - current(0)) + (guess - current));
            
            
            tmp_g(0) = 0;
            

            if(e.exp_value <= 0) {
                if(lastprevious <= 0) {
                    return current;
                }

                e = psi(previous);
                previous(0) = e.value;
                return previous;
            }
            else {
                previous = guess;
                lastprevious = psi(previous).exp_value;
                guess = guess - gamma * tmp_g;
            }
            

            
            /*{
                cout << " --- ITERATION " << i << " ---" << endl;
                cout << "Gradient: " << endl << e.gradient << endl;
                cout << "     exp: " << endl << exp(e.gradient) << endl;
                cout << "value: " << e.value << endl;
                cout << "  exp: " << exp(e.value) << endl;
                cout << "Guess: " << endl << guess << endl;
                cout << "  exp: " << endl << exp(guess) << endl;
                cout << "current: " << endl << current << endl;
                cout << "    exp: " << endl << exp(current) << endl;
                cout << "tmp_g: " << endl << tmp_g << endl;
                cout << "  exp: " << endl << exp(tmp_g) << endl;
                cout << "Diffs: " << e.value - current(0) << endl;
                cout << endl << endl << endl << endl;

            }*/
            


        } while(i < max_iter and tmp_g.norm() > precision);
        
        e = psi(guess);
        guess(0) = e.value;

        
        //cout << "Nb_iterations: " << i << endl;
        //cout << "Dernier gradient: " << tmp_g.norm() << endl;

        return guess;
    }
    
    
    
    
    
    
    Eigen::MatrixXd pseudoroot;
    
    params_type params;
    
    result last;
    result lastAntithetic;
    
    double step;
    double sqrt_step;
    
    double log_barrier;
    
    
    // Compute the swap rate associated to a given rate curve
    double internalSwapRate(VectorXd libor) {
        double sum = 0;
        double product = 1;
        for(int j = 0; j < params.dimension; ++j) {
            product /= 1 + params.delta * libor(j);
            sum += product;
        }
        return (1 - product) / (params.delta * sum);
    }
    
    
    
    
};

#endif /* GeneralSwaptionPricing_h */
