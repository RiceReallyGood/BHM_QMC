#include "dwa.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

int main()
{
    double t = 0.031;
    double beta = 1.0 / 0.01;
    double U = 1;
    double mu = 0.4 * U;
    unsigned short dim = 3;
    unsigned int Nl = 20;
    unsigned short max_state = 40;
    double time_of_flight = 0;
    double harmonic_shape = 0.0 * U;
    directed_worm_algorithm::count_type thermalization_sweeps = 5000000;
    directed_worm_algorithm::count_type measurement_sweeps = 5000000;
    unsigned short measurement_time_steps = beta / 0.1;
    /* directed_worm_algorithm(double beta_, double t_, double U_, double mu_, unsigned short dim_, unsigned int Nl_,
                            unsigned short max_state_, count_type _thermalization_sweeps, count_type _measurement_sweeps,
                            unsigned short measurement_time_steps_, double time_of_flight_, double harmonic_shape_); */
    directed_worm_algorithm mysystem(beta, t, U, mu, dim, Nl, max_state,
                                     thermalization_sweeps, measurement_sweeps, measurement_time_steps, time_of_flight, harmonic_shape);
    mysystem.simulation();
    mysystem.print_result();
    mysystem.print_greenfunction({0, 0, 0});
    mysystem.print_tof_greenfunction();

    //can also print to a file
    for (int k = 0; k < 10; k++)
    {
        ostringstream dir, filename;
        dir << "./t" << t << "/k" << k << "/";
        system(("mkdir -p " + dir.str()).c_str());
        filename << "greent" << t << "k" << k << ".txt";
        ofstream outgreen(dir.str() + filename.str());
        mysystem.print_greenfunction({k, 0, 0}, outgreen);
        /* ofstream outtofgreen("tof_green0.034.txt");
        mysystem.print_tof_greenfunction(outtofgreen); */
    }
}
