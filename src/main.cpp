#include "Header.hpp"

#include "Monte_Carlo.hpp"
#include "Timer.hpp"

#include "VanillaBarrier.hpp"
#include "GeneralCapletPricing.hpp"
#include "GeneralSwaptionPricing.hpp"


// Function to print option in a user menu
void print_option(int number, string name) {
    cout << "#" << number << ": " << name << endl;
}



// Description of the three following functions:
// Get an integer or a double and check if it is in good bounds
// The parameter incl indicates whether bounds are inclusive of exclusive

int get_answer(int min, int max) {
    int result = -1;
    cout << "Saisissez votre choix: " << flush;
    
    std::string s = "";
    
    while(result == -1) {
        try
        {
            cin >> s;
            result = std::stoi(s);
            if(result < min or result > max) {
                std::cout << "Erreur: l'entrée n'est pas un des choix proposés" << endl;
                std::cout << "Veuillez entrer votre choix: " << flush;
                result = -1;
            }
        }
        catch (std::invalid_argument const &e)
        {
            std::cout << "Erreur: l'entree est invalide" << endl;
            std::cout << "Veuillez entrer votre choix: " << flush;
            result = -1;
        }
        catch (std::out_of_range const &e)
        {
            std::cout << "Erreur: l'entrée n'est pas un des choix proposés" << endl;
            std::cout << "Veuillez entrer votre choix: " << flush;
            result = -1;
        }
    }
    
    return result;
}


double get_double(double min, double max, std::string add_str="", bool incl = true) {
    double result = 0;
    bool finished = false;
    
    std::string s = "";
    
    while(not finished) {
        try
        {
            cout << "Saisissez votre choix: " << add_str << flush;
            cin >> s;
            result = std::stod(s);
            if(result < min or result > max or ((not incl) and (result == min or result == max))) {
                std::cout << "Erreur: l'entree depasse les bornes autorisees" << endl;
            }
            else
                finished = true;
        }
        catch (std::invalid_argument const &e)
        {
            std::cout << "Erreur: l'entree est invalide" << endl;
            std::cout << "Veuillez entrer votre choix: " << add_str << endl;
        }
        catch (std::out_of_range const &e)
        {
            std::cout << "Erreur: l'entree depasse les bornes autorisees" << endl;
        }
    }
    
    return result;
}


double get_int(int min, int max, std::string add_str="", bool incl = true) {
    int result = 0;
    bool finished = false;
    
    std::string s = "";
    
    while(not finished) {
        try
        {
            cout << "Saisissez votre choix: " << add_str << flush;
            cin >> s;
            result = std::stoi(s);
            if(result < min or result > max or ((not incl) and (result == min or result == max))) {
                std::cout << "Erreur: l'entree depasse les bornes autorisees" << endl;
            }
            else
                finished = true;
        }
        catch (std::invalid_argument const &e)
        {
            std::cout << "Erreur: l'entree est invalide" << endl;
            std::cout << "Veuillez entrer votre choix: " << add_str << endl;
        }
        catch (std::out_of_range const &e)
        {
            std::cout << "Erreur: l'entree depasse les bornes autorisees" << endl;
        }
    }
    
    return result;
}






// Function to simulate the Black and Scholes barrier vanilla options
void user_simultation_BS(mt19937_64 &gen) {
    int answer = 0;
    
    VanillaBarrier<mt19937_64>::params_type p {
        .s0 = 100,
        .barrier = 95,
        .strike = 100,
        
        .T = 1,
        
        .sigma = 0.2,
        .r = 0.00,
        
        .n = 50 // Discretisation du B&S
    };

    
    cout << clear_terminal << flush;
    cout << "Commencons pas saisir les parametres du modele" << endl;
    
    cout << endl;
    cout << "Choix de la valeur initiale" << endl;
    p.s0 = get_double(0, max_double, "S_0 = ", false);
    
    cout << endl;
    cout << "Choix de la barriere" << endl;
    p.barrier = get_double(0, max_double, "U = ", false);

    cout << endl;
    cout << "Choix du strike" << endl;
    p.strike = get_double(0, max_double, "K = ", false);

    cout << endl;
    cout << "Choix de la maturite" << endl;
    p.T = get_double(0, max_double, "T = ", false);

    cout << endl;
    cout << "Choix de la volatilite" << endl;
    p.sigma = get_double(0, max_double, "sigma = ", false);

    cout << endl;
    cout << "Choix du taux d'interet" << endl;
    p.r = get_double(-1, max_double, "r = ", false);

    cout << endl;
    cout << "Type d'option" << endl;
    print_option(0, "Call");
    print_option(1, "Put");
    auto paramCP = VanillaBarrier<mt19937_64>::CALL;
    answer = get_answer(0,1);
    if(answer == 1)
        paramCP = VanillaBarrier<mt19937_64>::PUT;
    
    cout << endl;
    cout << "Type d'option barriere" << endl;
    print_option(0, "Down");
    print_option(1, "Up");
    auto paramUD = VanillaBarrier<mt19937_64>::DOWN;
    answer = get_answer(0,1);
    if(answer == 1)
        paramUD = VanillaBarrier<mt19937_64>::UP;
    
    cout << endl;
    cout << "Type d'option activante/désactivante" << endl;
    print_option(0, "In");
    print_option(1, "Out");
    auto paramIO = VanillaBarrier<mt19937_64>::IN;
    answer = get_answer(0,1);
    if(answer == 1)
        paramIO = VanillaBarrier<mt19937_64>::OUT;
    
    cout << endl;
    cout << "Nombre de pas de discretisation" << endl;
    p.n = get_int(1, max_int, "n = ", true);

    cout << endl;
    cout << "Nombre de simulation Monte Carlo (au moins 2)" << endl;
    int simulations = get_int(2, max_int, "M = ", true);


    bool finished = false;
    Timer t;

    
    do {
        cout << endl << endl;

        cout << "Menu de simulation" << endl;
        
        print_option(0, "Exit");
        print_option(1, "Simulation exacte basée sur un pont Brownien");
        print_option(2, "Simulation en discrétisant la formule explicite (sans pont Brownien)");
        print_option(3, "Simulation en discrétisant l'EDS (sans pont Brownien)");
        print_option(4, "Simulation en discrétisant l'EDS avec des variable uniforme sur {-1, 1}");
        print_option(5, "Changer le nombre de pas de discretisation");
        print_option(6, "Changer le nombre de simulations MC");
        
        int answer = get_answer(0,6);
        
        
        switch (answer) {
            case 0:
            {
                finished = true;
                break;
            }
                
            case 1:
            {
                BrownianBridgeVanillaBarrier<mt19937_64> vb_bb
                (paramCP, paramUD, paramIO, p);
                
                auto variable = bind(vb_bb, ref(gen));
                
                t.reset();
                auto res = monte_carlo(variable, simulations);
                t.stop();
                
                cout << t << endl;
                cout << res << endl;
                
                break;
            }
                

            case 2:
            {
                explicitVanillaBarrier<mt19937_64> vb_explicit
                (paramCP, paramUD, paramIO, p);
                
                auto variable = bind(vb_explicit, ref(gen));
                
                t.reset();
                auto res = monte_carlo(variable, simulations);
                t.stop();
                
                cout << t << endl;
                cout << res << endl;
                
                break;
            }

            case 3:
            {
                EulerVanillaBarrier<mt19937_64> vb_euler
                (paramCP, paramUD, paramIO, p);
                
                auto variable = bind(vb_euler, ref(gen));
                
                t.reset();
                auto res = monte_carlo(variable, simulations);
                t.stop();
                
                cout << t << endl;
                cout << res << endl;
                
                break;
            }

            case 4:
            {
                RWVanillaBarrier<mt19937_64> vb_rw
                (paramCP, paramUD, paramIO, p);
                
                auto variable = bind(vb_rw, ref(gen));
                
                t.reset();
                auto res = monte_carlo(variable, simulations);
                t.stop();
                
                cout << t << endl;
                cout << res << endl;
                
                break;
            }
                
            case 5: {
                cout << endl;
                cout << "Nombre de pas de discretisation" << endl;
                p.n = get_int(1, max_int, "n = ", true);
                
                break;
            }
                
            case 6: {
                cout << endl;
                cout << "Nombre de simulation Monte Carlo (au moins 2)" << endl;
                simulations = get_int(2, max_int, "M = ", true);
                
                break;
            }

            default:
                finished = true;
                break;
        }
    }
    
    while(not finished);
    
}




// Function to simulate the caplets in the libor market model
void user_simultation_caplet(mt19937_64 &gen) {
    double L0 = 0.15;
    double barrier = 0.20;
    double strike = 0.05;
    
    
    double T = 10;
    
    double volatility = 0.25;
    
    int n = 1;

    
    
    
    cout << clear_terminal << flush;
    cout << "Commencons pas saisir les parametres du modele" << endl;
    
    cout << endl;
    cout << "Choix de la valeur initiale" << endl;
    L0 = get_double(0, max_double, "L_0 = ", false);
    
    cout << endl;
    cout << "Choix de la barriere" << endl;
    barrier = get_double(0, max_double, "H = ", false);
    
    cout << endl;
    cout << "Choix du strike" << endl;
    strike = get_double(0, max_double, "K = ", false);
    
    cout << endl;
    cout << "Choix de la maturite" << endl;
    T = get_double(0, max_double, "T = ", false);
    
    cout << endl;
    cout << "Choix de la volatilite" << endl;
    volatility = get_double(0, max_double, "sigma = ", false);
    
    
    cout << endl;
    cout << "Nombre de pas de discretisation" << endl;
    n = get_int(1, max_int, "n = ", true);
    
    cout << endl;
    cout << "Nombre de simulation Monte Carlo (au moins 2)" << endl;
    int simulations = get_int(2, max_int, "M = ", true);
    
    
    
    GeneralCapletPricing<mt19937_64>::params_type params {
        .volatitly = [volatility](double tmp) {return volatility;},
        .F = [](double t, double lnl) { return 0.0;},
        
        .L0 = L0,
        .barrier = barrier,
        .strike = strike,
        
        .n = n,
        
        .T = T
    };

    normal nd{};
    
    auto var_red = [&nd, barrier, strike, T, volatility](double t, double lnl) {
        t = T-t;
        double sqrtt = sqrt(t);
        double v = volatility * sqrtt;
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
        
        return - volatility * rate *  sum;
    };
    
    
    bool finished = false;
    Timer t;
    
    GeneralCapletPricing<mt19937_64> llm1d(params);
    
    
    do {
        cout << endl << endl;
        cout << "Menu de simulation" << endl;
        
        print_option(0, "Exit");
        
        print_option(1, "Simulation classique");
        print_option(2, "Simulation classique antithetique");
        print_option(3, "Algorithme des marches aléatoires d'ordre 1/2");
        print_option(4, "Algorithme des marches aléatoires antithetique d'ordre 1/2");
        print_option(5, "Algorithme des marches aléatoires d'ordre 1");
        print_option(6, "Algorithme des marches aléatoires antithetique d'ordre 1");
        print_option(7, "Algorithme basé sur un pont Brownien (avec discrétisation)");
        print_option(8, "Algorithme basé sur un pont Brownien antithetique (avec discrétisation)");
        
        print_option(9, "Reduction de variance (activation/desactivation");
        print_option(10, "Changer le nombre de pas de discretisation");
        print_option(11, "Changer le nombre de simulations MC");
        
        int answer = get_answer(0,11);
        
        cout << endl;

        
        switch (answer) {
            case 0:
                finished = true;
                break;
                
            case 1:
            {
                auto variable = [&llm1d, &gen]() {
                    return llm1d.normal_simulation(gen);
                };
                
                t.reset();
                auto res = monte_carlo(variable, simulations);
                t.stop();
                
                cout << t << endl;
                cout << res << endl;
                
                break;
            }

                
            case 2:
            {
                auto variable = [&llm1d, &gen]() {
                    return llm1d.antithetic_normal_simulation(gen);
                };
                
                
                t.reset();
                auto res = monte_carlo(variable, simulations);
                t.stop();
                
                cout << t << endl;
                cout << res << endl;
                
                break;
            }

                
            case 3:
            {
                auto variable = [&llm1d, &gen]() {
                    return llm1d.slow_rw_simulation(gen);
                };
                
                t.reset();
                auto res = monte_carlo(variable, simulations);
                t.stop();
                
                cout << t << endl;
                cout << res << endl;
                
                break;
            }

                
            case 4:
            {
                auto variable = [&llm1d, &gen]() {
                    return llm1d.antithetic_slow_rw_simulation(gen);
                };
                
                t.reset();
                auto res = monte_carlo(variable, simulations);
                t.stop();
                
                cout << t << endl;
                cout << res << endl;
                
                break;
            }

                
            case 5:
            {
                auto variable = [&llm1d, &gen]() {
                    return llm1d.rw_simulation(gen);
                };
                
                t.reset();
                auto res = monte_carlo(variable, simulations);
                t.stop();
                
                cout << t << endl;
                cout << res << endl;
                
                break;
            }

                
            case 6:
            {
                auto variable = [&llm1d, &gen]() {
                    return llm1d.antithetic_rw_simulation(gen);
                };
                
                t.reset();
                auto res = monte_carlo(variable, simulations);
                t.stop();
                
                cout << t << endl;
                cout << res << endl;
                
                break;
            }

                
                
            case 7:
            {
                auto variable = [&llm1d, &gen]() {
                    return llm1d.bridge_normal_simulation(gen);
                };
                
                t.reset();
                auto res = monte_carlo(variable, simulations);
                t.stop();
                
                cout << t << endl;
                cout << res << endl;
                
                break;
            }
                

                
            case 8:
            {
                auto variable = [&llm1d, &gen]() {
                    return llm1d.antithetic_bridge_normal_simulation(gen);
                };
                
                t.reset();
                auto res = monte_carlo(variable, simulations);
                t.stop();
                
                cout << t << endl;
                cout << res << endl;
                
                break;
            }
                

                
                
                
                
            case 9: {
                cout << "Voulez vous activer la réduction de variance (basée sur la fonction F)" << endl;
                
                print_option(1, "Oui");
                print_option(2, "Non");
                
                int answer = get_answer(1,2);
                
                switch (answer) {
                    case 1:
                        params.F = var_red;
                        break;
                        
                    case 2:
                        params.F = [](double t, double lnl) { return 0.0;};
                        break;
                        
                    default:
                        finished = true;
                        break;
                }
                
                llm1d = GeneralCapletPricing<mt19937_64> (params);
                
                break;
            }
                
            case 10: {
                cout << "Nombre de pas de discretisation" << endl;
                n = get_int(1, max_int, "n = ", true);
                params.n = n;
                
                llm1d = GeneralCapletPricing<mt19937_64> (params);

                break;
            }
                
            case 11: {
                cout << "Nombre de simulation Monte Carlo (au moins 2)" << endl;
                simulations = get_int(2, max_int, "M = ", true);
                
                break;
            }
                
            default:
                finished = true;
                break;
        }
    }
    while(not finished);
    
}




// Function to simulate the swaption in the libor market model
void user_simultation_swaption(mt19937_64 &gen) {
    cout << clear_terminal << flush;
    
    cout << "Commencons pas saisir les parametres du modele" << endl;
    
    
    cout << endl;
    cout << "Dimension (nombre d'echeances considerees)" << endl;
    int dimension = get_int(0, max_int, "dim = ", false);
    
    cout << endl;
    cout << "Temps entre deux echanges du swap" << endl;
    double delta = get_double(0, max_double, "delta = ", false);
    
    VectorXd initial = MatrixXd::Constant(dimension,1,0.05);
    
    cout << endl;
    cout << "Choix des " << dimension << " valeurs initiales" << endl;
    for(int i = 0; i<dimension; ++i) {
        initial(i) = get_double(0, max_double, "L_0(" + to_string(i+1) + ") = ", false);
    }
    
    cout << endl;
    cout << "Choix de la barriere" << endl;
    double barrier = get_double(0, max_double, "H = ", false);

    cout << endl;
    cout << "Choix du strike" << endl;
    double strike = get_double(0, max_double, "K = ", false);
    
    cout << endl;
    cout << "Choix de la maturite (début du swap)" << endl;
    double T0 = get_double(0, max_double, "T0 = ", false);
    
    cout << endl;
    cout << "Choix de la volatilite de chaque composante (supposée constante)" << endl;
    VectorXd volatility = MatrixXd::Constant(dimension,1,0.05);
    for(int i = 0; i<dimension; ++i) {
        volatility(i) = get_double(0, max_double, "sigma(" + to_string(i+1) + ") = ", false);
    }
    
    cout << endl;
    cout << "Choix du parametre beta de la matrice de variance covariance" << endl;
    double beta = get_double(0, max_double, "beta = ", false);

    
    cout << endl;
    cout << "Nombre de pas de discretisation" << endl;
    int n = get_int(1, max_int, "n = ", true);
    
    cout << endl;
    cout << "Nombre de simulation Monte Carlo (au moins 2)" << endl;
    int simulations = get_int(2, max_int, "M = ", true);


    // Compute correlation matrix according user's beta
    MatrixXd corr = MatrixXd::Zero(dimension, dimension);
    for(int i = 0; i<dimension; ++i) {
        for(int j = 0; j<i; ++j) {
            corr(i,j) = corr(j,i) = exp(beta*(j-i)*delta);
        }
        corr(i,i) = 1;
    }

    
    GeneralSwaptionPricing<mt19937_64>::params_type params {
        .dimension = dimension,
        
        .volatitly = [dimension, &volatility] (double t) { return volatility; },
        .F = [&dimension](double t, VectorXd c) { return MatrixXd::Zero(dimension,1);},
        
        .initial = initial,
        .barrier = barrier,
        
        .strike = strike,
        .correlations = corr,
        
        
        .n = n,
        .T0 = T0,
        .delta = delta
        
    };
    
    bool finished = false;
    Timer t;
    
    GeneralSwaptionPricing<mt19937_64> gsp(params);
    
    do {
        cout << endl << endl;
        cout << "Menu de simulation" << endl;
        
        print_option(0, "Exit");
        
        print_option(1, "Simulation classique");
        print_option(2, "Simulation classique antithetique");
        print_option(3, "Algorithme des marches aléatoires d'ordre 1/2");
        print_option(4, "Algorithme des marches aléatoires antithetique d'ordre 1/2");
        print_option(5, "Algorithme des marches aléatoires d'ordre 1");
        print_option(6, "Algorithme des marches aléatoires antithetique d'ordre 1");
        
        print_option(7, "Reduction de variance (activation/desactivation");
        print_option(8, "Changer le nombre de pas de discretisation");
        print_option(9, "Changer le nombre de simulations MC");
        
        int answer = get_answer(0,9);
        
        cout << endl;
        
        
        switch (answer) {
            case 0:
                finished = true;
                break;
                
            case 1:
            {
                auto variable = [&gsp, &gen]() {
                    return gsp.normal_simulation(gen);
                };
                
                t.reset();
                auto res = monte_carlo(variable, simulations);
                t.stop();
                
                cout << t << endl;
                cout << res << endl;
                
                break;
            }
                
                
            case 2:
            {
                auto variable = [&gsp, &gen]() {
                    return gsp.normal_antithetic_simulation(gen);
                };
                
                
                t.reset();
                auto res = monte_carlo(variable, simulations);
                t.stop();
                
                cout << t << endl;
                cout << res << endl;
                
                break;
            }
                
                
            case 3:
            {
                auto variable = [&gsp, &gen]() {
                    return gsp.slow_rw_simulation(gen);
                };
                
                t.reset();
                auto res = monte_carlo(variable, simulations);
                t.stop();
                
                cout << t << endl;
                cout << res << endl;
                
                break;
            }
                
                
            case 4:
            {
                auto variable = [&gsp, &gen]() {
                    return gsp.slow_rw_antithetic_simulation(gen);
                };
                
                t.reset();
                auto res = monte_carlo(variable, simulations);
                t.stop();
                
                cout << t << endl;
                cout << res << endl;
                
                break;
            }
                
                
            case 5:
            {
                auto variable = [&gsp, &gen]() {
                    return gsp.rw_simulation(gen);
                };
                
                t.reset();
                auto res = monte_carlo(variable, simulations);
                t.stop();
                
                cout << t << endl;
                cout << res << endl;
                
                break;
            }
                
                
            case 6:
            {
                auto variable = [&gsp, &gen]() {
                    return gsp.rw_antithetic_simulation(gen);
                };
                
                t.reset();
                auto res = monte_carlo(variable, simulations);
                t.stop();
                
                cout << t << endl;
                cout << res << endl;
                
                break;
            }
                
                
                
                
                
            case 7: {
                cout << "Voulez vous activer la réduction de variance (basée sur la fonction F)" << endl;
                
                print_option(1, "Oui");
                print_option(2, "Non");
                
                int answer = get_answer(1,2);
                
                switch (answer) {
                    case 1:
                        params.F = [&params](double t, VectorXd c) { return params.optimalF_constantvol(t, c); };
                        break;
                        
                    case 2:
                        params.F = [&dimension](double t, VectorXd c) { return MatrixXd::Zero(dimension,1);};
                        break;
                        
                    default:
                        finished = true;
                        break;
                }
                
                gsp = GeneralSwaptionPricing<mt19937_64>(params);
                
                break;
            }
                
            case 8: {
                cout << "Nombre de pas de discretisation" << endl;
                n = get_int(1, max_int, "n = ", true);
                params.n = n;
                
                gsp = GeneralSwaptionPricing<mt19937_64>(params);
                
                break;
            }
                
            case 9: {
                cout << "Nombre de simulation Monte Carlo (au moins 2)" << endl;
                simulations = get_int(2, max_int, "M = ", true);
                
                break;
            }
                
            default:
                finished = true;
                break;
        }
    }
    while(not finished);
    
}









int main(int argc, const char * argv[]) {
    // Seeds of the randomness in the program
    random_device rd;
    auto seed = rd();
    mt19937_64 gen(seed);

    

    bool finished = false;
    
    do {
        cout << clear_terminal << flush;
        
        cout << "Menu principal" << endl;
        
        print_option(0, "Exit");
        print_option(1, "Simulation d'une option barriere Black and Scholes");
        print_option(2, "Simulation d'un caplet dans le modele LMM");
        print_option(3, "Simulation d'un swaption dans le modele LMM");
        
        int answer = get_answer(0,3);
        
        switch (answer) {
            case 0:
                finished = true;
                break;

            case 1:
                user_simultation_BS(gen);
                break;

            case 2:
                user_simultation_caplet(gen);
                break;
                
            case 3:
                user_simultation_swaption(gen);
                break;

            default:
                finished = true;
                break;
        }
    }
    while(not finished);
    
    return 0;
}




