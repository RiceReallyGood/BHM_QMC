#ifndef DWA_HPP
#define DWA_HPP

//2019.06.26: change the definition of green's function from -<a_{i}(\tau)a_{j}^{\dagger}> to <a_{i}^{\dagger}a_{j}>

#include <vector>
#include <cmath>
#include <random>
#include "worldlines.hpp"
#include <iomanip>
#include <complex>

class directed_worm_algorithm
{
public:
    typedef std::uint64_t count_type;
    typedef worldlines::line line;
    typedef worldlines::location_type location_type;
    directed_worm_algorithm(double beta_, double t_, double U_, double mu_, unsigned short dim_, unsigned int Nl_,
                            unsigned short max_state_, count_type _thermalization_sweeps, count_type _measurement_sweeps,
                            unsigned short measurement_time_steps_, double time_of_flight_, double harmonic_shape_);

    void simulation();
    void print_result(std::ostream &os = std::cout) const;
    void print_greenfunction(const std::vector<int> &momentum, std::ostream &os = std::cout) const;
    void print_tof_greenfunction(std::ostream &os = std::cout) const;

private:
    void sweep();
    bool wormhead_propagates_till_collision_with_wormtail();
    void insert_jump_or_bounce(double diagonal_energy_relative_);
    void delete_relink_jump_or_bounce();

    void do_diagonal_measurement();
    void do_greenfunction_measurement(double newtime_);

    double diagonal_energy_relative(unsigned int site_, unsigned short state_, bool forward_, bool creation_) const
    {
        //std::cout << "calculate the diagonal_energy_relative of site " << site_ << " state " << state_ << " " << (creation_ ? "creation " : "annihilation ") << (forward_ ? "forward" : "backward") << std::endl;
        double _diagonal_energy = _onsite_energy[site_][state_];
        double _diagonal_energy_before = _onsite_energy[site_][creation_ ? state_ - 1 : state_ + 1];
        return forward_ ? _diagonal_energy_before - std::min(_diagonal_energy, _diagonal_energy_before) + _energy_offset
                        : _diagonal_energy - std::min(_diagonal_energy, _diagonal_energy_before) + _energy_offset;
    }

    std::vector<int> site2coordinate(unsigned int site_) const
    {
        std::vector<int> coordinate_;
        coordinate_.reserve(_dim);
        while (site_)
        {
            coordinate_.push_back(site_ % _Nl);
            site_ /= _Nl;
        }
        coordinate_.resize(_dim);
        return coordinate_;
    }

    unsigned int coordinate2site(std::vector<int> coordinate_) const
    {
        unsigned int site_ = 0;
        while (!coordinate_.empty())
        {
            site_ = site_ * _Nl + coordinate_.back();
            coordinate_.pop_back();
        }
        return site_;
    }

    unsigned int displacement(unsigned int site1_, unsigned int site2_) const
    {
        std::vector<int> displacement_coordinate;
        displacement_coordinate.reserve(_dim);
        for (int i = 0; i < _dim; i++)
        {
            displacement_coordinate.push_back((_Nl + (site1_ % _Nl) - (site2_ % _Nl)) % _Nl);
            site1_ /= _Nl, site2_ /= _Nl;
        }
        return coordinate2site(displacement_coordinate);
    }

    double square_displacement(unsigned int site1_, unsigned int site2_) const
    {
        double sd = 0;
        for (int i = 0; i < _dim; i++)
        {
            double coordinate1 = site1_ % _Nl;
            double coordinate2 = site2_ % _Nl;
            sd += coordinate1 * coordinate1 - coordinate2 * coordinate2 - _Nl * (coordinate1 - coordinate2);
            site1_ /= _Nl, site2_ /= _Nl;
        }
        return sd;
    }

    int periodic(int coordinate_) { return (coordinate_ + _Nl) % _Nl; };
    bool is_thermalized() const { return (_sweep_counter > _thermalization_sweeps); }

    //parameters
    double _beta; //one over temperature
    double _hopping_strength;
    double _onsite_coupling;
    double _energy_offset = 0.1;

    unsigned short _dim;
    unsigned int _Nl;
    unsigned int _num_sites;
    unsigned short _max_state;

    //bool _is_homogeneous;

    std::vector<double> _chemical_potential;            //-agr0: site
    std::vector<std::vector<double>> _onsite_energy;    //-arg0: site  arg1: state
    std::vector<std::vector<unsigned int>> _neighbours; //-arg0: site  arg1: which neighbour
    std::vector<double> _measurement_times;
    std::vector<std::vector<double>> _greenfunction;
    std::vector<std::complex<double>> _fourier_coeff;
    //regarding TOF
    std::vector<std::complex<double>> _tof_greenfunction;
    double _time_of_flight;
    double _time_of_flight_phase;

    worldlines wl;
    wormpair worm;

    count_type _thermalization_sweeps;
    count_type _measurement_sweeps;
    count_type _sweep_counter;

    std::default_random_engine _dre;
    std::uniform_real_distribution<double> _urd;

    double _average_particle_num;
    double _winding_number2_in_x_direction; //winding number square in x direction
    unsigned short _measurement_time_steps;
};

// ==================================================
// directed_worm_algorithm member functions
// ==================================================

directed_worm_algorithm::directed_worm_algorithm(double beta_, double t_, double U_, double mu_, unsigned short dim_, unsigned int Nl_,
                                                 unsigned short max_state_, count_type thermalization_sweeps_, count_type measurement_sweeps_,
                                                 unsigned short measurement_time_steps_, double time_of_flight_, double harmonic_shape_)
    : _beta(beta_),
      _hopping_strength(t_ * _beta),
      _onsite_coupling(U_ * _beta),
      _dim(dim_),
      _Nl(Nl_),
      _max_state(max_state_),
      _thermalization_sweeps(thermalization_sweeps_),
      _measurement_sweeps(measurement_sweeps_),
      _sweep_counter(0),
      _num_sites(std::pow(_Nl, _dim)),
      _chemical_potential(_num_sites, mu_ * _beta),
      _measurement_time_steps(measurement_time_steps_),
      _greenfunction(_num_sites, std::vector<double>(measurement_time_steps_ + 1, 0)),
      _time_of_flight(time_of_flight_),
      _time_of_flight_phase(time_of_flight_),
      //_time_of_flight_phase(87. * 1.67 * 5.32 * 5.32 / (2 * 1.0545718 * time_of_flight_ * 10000000)),
      _tof_greenfunction(_num_sites, 0),
      _onsite_energy(_num_sites, std::vector<double>(_max_state + 1, 0)),
      _neighbours(_num_sites),
      wl(_num_sites),
      _dre(time(NULL)), //time(NULL)
      _average_particle_num(0),
      _winding_number2_in_x_direction(0)
{
    //initiallize _measurement_times
    double time_step = 1. / measurement_time_steps_;
    _measurement_times.reserve(measurement_time_steps_);
    for (int i = 0; i <= measurement_time_steps_; ++i)
        _measurement_times.push_back(i * time_step);

    //initiallize _chemical_potential
    for (unsigned int site = 0; site < _num_sites; ++site)
    {
        std::vector<int> coo = site2coordinate(site);
        double distance2 = 0;
        for (int i = 0; i < _dim; ++i)
            distance2 += (coo[i] - _Nl / 2.) * (coo[i] - _Nl / 2.);
        _chemical_potential[site] -= _beta * harmonic_shape_ * distance2;
    }

    //initiallize _onsite_energy
    for (unsigned int site = 0; site < _num_sites; ++site)
    {
        for (unsigned short state = 1; state <= _max_state; ++state)
            _onsite_energy[site][state] = 0.5 * _onsite_coupling * state * (state - 1) - _chemical_potential[site] * state;
    }

    //initiallize _neighbours
    for (unsigned int site = 0; site < _num_sites; ++site)
    {
        std::vector<int> this_site_neighbours_coordinate = site2coordinate(site);
        for (unsigned short dir = 0; dir < _dim; ++dir)
        {
            this_site_neighbours_coordinate[dir] = periodic(this_site_neighbours_coordinate[dir] - 1);
            _neighbours[site].push_back(coordinate2site(this_site_neighbours_coordinate));
            this_site_neighbours_coordinate[dir] = periodic(this_site_neighbours_coordinate[dir] + 2);
            _neighbours[site].push_back(coordinate2site(this_site_neighbours_coordinate));
            this_site_neighbours_coordinate[dir] = periodic(this_site_neighbours_coordinate[dir] - 1);
        }
    }

    //initiallize _fourier_coeff
    std::complex<double> I(0, 1);
    double pi = 3.1415926535897932385;
    for (int i = 0; i < Nl_; ++i)
        _fourier_coeff.push_back(exp(I * 2. * pi * static_cast<double>(i) / static_cast<double>(_Nl)));
}

void directed_worm_algorithm::simulation()
{
    while (_sweep_counter < _thermalization_sweeps + _measurement_sweeps)
    {
        ++_sweep_counter;
        sweep();
        std::cout << _sweep_counter << " " << wl.total_interactions() << std::endl;
        //if (!wl.is_valid(_max_state))
        //    return;
        //wl.print();
    }
}

void directed_worm_algorithm::sweep()
{
    unsigned int _site = _num_sites * _urd(_dre);
    double _time = _urd(_dre);
    bool _forward = _urd(_dre) < 0.5;
    bool _creation = _urd(_dre) < 0.5;

    location_type _location = wl.location(_site, _time);
    unsigned short _state = wl.state_before(_location);

    /* if ((wormpair::increasing(_forward, _creation) && _state == _max_state) ||
        (!wormpair::increasing(_forward, _creation) && _state == 0) ||
        (!wl.location_is_kink_unoccupied(_location, _time)))
    {
        --_sweep_counter;
        return;
    } */

    if (!wl.location_is_kink_unoccupied(_location, _time))
    {
        --_sweep_counter;
        return;
    }

    if (_state == 0 && !wormpair::increasing(_forward, _creation))
    {
        if (is_thermalized())
            do_diagonal_measurement();
        return;
    }

    if (_state == _max_state && wormpair::increasing(_forward, _creation))
    {
        if (is_thermalized())
        {
            _greenfunction[0][0] += 2 * (_state + 1);
            do_diagonal_measurement();
        }
        return;
    }
    //std::cout << "sweep counter = " << _sweep_counter << " insert a wormpair at site " << _site << " time " << _time
    //          << " wormhead type is " << (_creation ? "creation " : "annihilation ") << (_forward ? "forward" : "backward") << std::endl;

    worm = wormpair(_location, kink(_site, _time, _state), _forward, _creation);

    if (is_thermalized())
    {
        if (worm.increasing())
            _greenfunction[0].back() += worm.wormpair_state();
        else
        {
            _greenfunction[0][0] += worm.wormpair_state();
            _tof_greenfunction[0] += worm.wormpair_state();
        }
    }

    while (wormhead_propagates_till_collision_with_wormtail())
        ;
    //std::cout << "reach a diagonal configuration" << std::endl;

    /* if (is_thermalized() && !worm.increasing())
        _greenfunction[0][0] += worm.wormpair_state(); */
    if (is_thermalized())
    {
        if (!worm.increasing())
            _greenfunction[0].back() += worm.wormpair_state();
        else
        {
            _greenfunction[0][0] += worm.wormpair_state();
            _tof_greenfunction[0] += worm.wormpair_state();
        }
    }

    worm.wormhead_annihilates_wormtail();
    if (is_thermalized())
        do_diagonal_measurement();
    return;
}

bool directed_worm_algorithm::wormhead_propagates_till_collision_with_wormtail()
{
    double _time2next = worm.time2next();
    double _diagonal_energy_relative = diagonal_energy_relative(worm.wormhead_site(), worm.wormhead_state(), worm.forward(), worm.creation());
    double _deltatime = -std::log(1. - _urd(_dre)) / _diagonal_energy_relative;
    bool _halted = (_deltatime >= _time2next);
    if (_halted)
        _deltatime = _time2next;

    double _newtime = worm.forward() ? worm.wormhead_time() + _deltatime : worm.wormhead_time() - _deltatime;
    bool _winding_over_time = (_newtime < 0. || _newtime >= 1.);

    if (_halted)
        _newtime = worm.forward() ? worm.next_time() - std::numeric_limits<double>::epsilon() : worm.next_time() + std::numeric_limits<double>::epsilon();
    else
        _newtime = wormpair::mod_one(_newtime);

    if (is_thermalized())
        do_greenfunction_measurement(_newtime);

    worm.wormhead_moves_to_new_time(_newtime, _winding_over_time);

    if (!_halted)
    {
        insert_jump_or_bounce(_diagonal_energy_relative);
        return true;
    }
    else //halted
    {
        if (worm.wormhead_reaches_wormtail())
        {
            //std::cout << "wormhead reaches wormtail" << std::endl;
            return false;
        }

        else if (worm.wormhead_is_same_type_as_next())
        {
            //std::cout << "wormhead is same type as next" << std::endl;
            worm.wormhead_crosses_vertex();
            return true;
        }

        else
        {
            //std::cout << "wormhead is different type with next" << std::endl;
            delete_relink_jump_or_bounce();
            return true;
        }
    }
}

void directed_worm_algorithm::do_greenfunction_measurement(double newtime_)
{
    static const std::complex<double> I(0, 1);
    unsigned int green_site; // = worm.creation() ? displacement(worm.wormtail_site(), worm.wormhead_site()) : displacement(worm.wormhead_site(), worm.wormtail_site());
    double start_green_time, end_green_time;
    bool green_time_forward, green_time_winding_over_time;
    if (worm.creation())
    {
        green_site = displacement(worm.wormhead_site(), worm.wormtail_site());
        start_green_time = wormpair::mod_one(worm.wormhead_time() - worm.wormtail_time());
        end_green_time = wormpair::mod_one(newtime_ - worm.wormtail_time());
        green_time_forward = worm.forward();
    }
    else
    {
        green_site = displacement(worm.wormtail_site(), worm.wormhead_site());
        start_green_time = wormpair::mod_one(worm.wormtail_time() - worm.wormhead_time());
        end_green_time = wormpair::mod_one(worm.wormtail_time() - newtime_);
        green_time_forward = !worm.forward();
    }
    green_time_winding_over_time = (green_time_forward && (end_green_time < start_green_time)) ||
                                   (!green_time_forward && (end_green_time > start_green_time));

    if (green_time_winding_over_time)
    {
        double tof_green_phase;
        if (worm.creation())
            tof_green_phase = _time_of_flight_phase * square_displacement(worm.wormhead_site(), worm.wormtail_site());
        else
            tof_green_phase = _time_of_flight_phase * square_displacement(worm.wormtail_site(), worm.wormhead_site());

        _tof_greenfunction[green_site] += std::exp(-I * tof_green_phase) * static_cast<double>(worm.wormpair_state());
    }

    //std::cout << "start time = " << start_green_time << " end green time = " << end_green_time << std::endl;
    int measure_time_steps = _measurement_times.size();
    if (green_time_forward)
    {
        int start_index = std::upper_bound(_measurement_times.begin(), _measurement_times.end(), start_green_time) - _measurement_times.begin();
        int end_index = std::upper_bound(_measurement_times.begin(), _measurement_times.end(), end_green_time) - _measurement_times.begin();
        if (green_time_winding_over_time)
        {
            for (int index_ = start_index; index_ < measure_time_steps; ++index_)
                _greenfunction[green_site][index_] += worm.wormpair_state();
            for (int index_ = 0; index_ < end_index; ++index_)
                _greenfunction[green_site][index_] += worm.wormpair_state();
        }
        else
        {
            for (int index_ = start_index; index_ < end_index; ++index_)
                _greenfunction[green_site][index_] += worm.wormpair_state();
        }
        //std::cout << "start index = " << start_index << " end index = " << end_index << " direction " << (green_time_forward ? "forward" : "backward") << std::endl;
    }
    else
    {
        int start_index = std::lower_bound(_measurement_times.begin(), _measurement_times.end(), start_green_time) - _measurement_times.begin();
        int end_index = std::lower_bound(_measurement_times.begin(), _measurement_times.end(), end_green_time) - _measurement_times.begin();
        if (green_time_winding_over_time)
        {
            for (int index_ = start_index - 1; index_ >= 0; --index_)
                _greenfunction[green_site][index_] += worm.wormpair_state();
            for (int index_ = measure_time_steps - 1; index_ >= end_index; --index_)
                _greenfunction[green_site][index_] += worm.wormpair_state();
        }
        else
        {
            for (int index_ = start_index - 1; index_ >= end_index; --index_)
                _greenfunction[green_site][index_] += worm.wormpair_state();
        }
        //std::cout << "start index = " << start_index << " end index = " << end_index << " direction " << (green_time_forward ? "forward" : "backward") << std::endl;
    }
}

void directed_worm_algorithm::insert_jump_or_bounce(double diagonal_energy_relative_)
{
    unsigned int _neighbour_site = _neighbours[worm.wormhead_site()][_urd(_dre) * 2 * _dim];
    location_type _neighborlocation = wl.location(_neighbour_site, worm.wormhead_time());
    unsigned short _targetstate = (_neighborlocation.second - 1)->state();

    if ((worm.increasing() && _targetstate == _max_state) ||
        (!worm.increasing() && _targetstate == 0) ||
        (!wl.location_is_kink_unoccupied(_neighborlocation, worm.wormhead_time())))
    {
        worm.wormhead_turns_around();
        //std::cout << "wormhead turns around" << std::endl;
        return;
    }
    else
    {
        double _new_weight = _hopping_strength * (worm.increasing() ? _targetstate + 1 : _targetstate);
        if (_new_weight >= diagonal_energy_relative_ || _urd(_dre) < (_new_weight / diagonal_energy_relative_))
        {
            worm.wormhead_inserts_vertex_and_jumps_to_new_site(_neighborlocation);
            return;
        }
        else
        {
            worm.wormhead_turns_around();
            //std::cout << "wormhead turns around" << std::endl;
            return;
        }
    }
}

void directed_worm_algorithm::delete_relink_jump_or_bounce()
{
    unsigned int _source_site = worm.next_partnersite();
    location_type _sourcelocation = wl.location(_source_site, worm.next_time());
    unsigned int _target_site = _neighbours[_source_site][_urd(_dre) * 2 * _dim];

    double _weight = _hopping_strength * (worm.creation() ? worm.wormhead_state() : worm.next_state());

    if (_target_site == worm.wormhead_site()) //delete or bounce
    {
        double _new_weight = diagonal_energy_relative(_source_site, _sourcelocation.second->state(), !worm.forward(), worm.creation());
        if (_new_weight >= _weight || _urd(_dre) < (_new_weight / _weight))
        {
            worm.wormhead_deletes_vertex_and_jumps_to_new_site(_sourcelocation);
            return;
        }
        else
        {
            worm.wormhead_turns_around();
            //std::cout << "wormhead turns around" << std::endl;
            return;
        }
    }
    else //relink
    {
        location_type _targetlocation = wl.location(_target_site, worm.next_time());
        unsigned short _target_state = (_targetlocation.second - 1)->state();
        if ((!worm.increasing() && _target_state == _max_state) ||
            (worm.increasing() && _target_state == 0) ||
            !wl.location_is_kink_unoccupied(_targetlocation, worm.next_time()))
        {
            worm.wormhead_turns_around();
            //std::cout << "wormhead turns around" << std::endl;
            return;
        }

        double _new_weight = _hopping_strength * (!worm.increasing() ? _target_state + 1 : _target_site);
        if (_new_weight >= _weight || _urd(_dre) < (_new_weight / _weight))
        {
            worm.wormhead_relinks_vertex_and_jumps_to_new_site(_sourcelocation, _targetlocation);
            return;
        }
        else
        {
            worm.wormhead_turns_around();
            //std::cout << "wormhead turns around" << std::endl;
            return;
        }
    }
}

void directed_worm_algorithm::do_diagonal_measurement()
{
    //measure average particle num
    double average_n = 0;
    for (unsigned int site = 0; site < _num_sites; ++site)
        average_n += wl.site_state(site);
    average_n /= _num_sites;

    _average_particle_num += average_n;

    //measure winding number square in x direction
    double _net_number_of_directed_hops = 0;
    for (unsigned int other_dir = 0; other_dir < std::pow(_Nl, _dim - 1); ++other_dir)
        _net_number_of_directed_hops += wl.net_number_of_directed_hops(other_dir * _Nl, other_dir * _Nl + 1);
    _winding_number2_in_x_direction += _net_number_of_directed_hops * _net_number_of_directed_hops;
}

void directed_worm_algorithm::print_result(std::ostream &os) const
{
    std::cout << "average particle number = " << _average_particle_num / _measurement_sweeps << std::endl;
    std::cout << "superfluid density = " << _winding_number2_in_x_direction * _Nl * _Nl / (_beta * _num_sites * _average_particle_num) << std::endl;
}

void directed_worm_algorithm::print_greenfunction(const std::vector<int> &momentum, std::ostream &os) const
{
    std::ios::fmtflags oldFlags = os.flags();
    double time_step = _beta / _measurement_time_steps;
    for (int time_slice = 0; time_slice <= _measurement_time_steps; ++time_slice)
    {
        std::complex<double> gktau = 0;
        for (unsigned int site = 0; site < _num_sites; ++site)
        {
            std::vector<int> coordinate = site2coordinate(site);
            std::complex<double> coeff(1, 0);
            for (int dir = 0; dir < _dim; ++dir)
                coeff *= _fourier_coeff[(coordinate[dir] * momentum[dir]) % _Nl];
            gktau += _greenfunction[site][time_slice] * coeff;
        }
        os << setiosflags(std::ios::fixed | std::ios::right)
           << std::setprecision(2) << std::setw(5) << time_step * time_slice << " "
           << std::setprecision(6) << std::setw(10) << gktau.real() / _measurement_sweeps << " "
           << std::setw(10) << gktau.imag() / _measurement_sweeps << std::endl;
    }
    os.flags(oldFlags);
}

/* void directed_worm_algorithm::print_greenfunction(const std::vector<int> &momentum, std::ostream &os) const
{
    std::ios::fmtflags oldFlags = os.flags();
    std::vector<int> coordinate(_dim, 0);
    while (true)
    {
        unsigned int site = coordinate2site(coordinate);
        for (int i = 0; i < _dim; ++i)
            os << std::setw(2) << coordinate[i] << " ";
        os << setiosflags(std::ios::fixed | std::ios::right) << std::setprecision(6)
           << std::setw(10) << _greenfunction[site][0] / _measurement_sweeps << std::endl;

        //next coordinate
        int dir = _dim - 1;
        while (dir >= 0 && coordinate[dir] == _Nl - 1)
            coordinate[dir] = 0, --dir;
        if (dir == -1)
            break;
        ++coordinate[dir];
    }
    os.flags(oldFlags);
} */

void directed_worm_algorithm::print_tof_greenfunction(std::ostream &os) const
{
    std::ios::fmtflags oldFlags = os.flags();
    std::vector<int> coordinate(_dim, 0);
    while (true)
    {
        unsigned int site = coordinate2site(coordinate);
        for (int i = 0; i < _dim; ++i)
            os << std::setw(2) << coordinate[i] << " ";
        os << setiosflags(std::ios::fixed | std::ios::right) << std::setprecision(6)
           << std::setw(10) << _tof_greenfunction[site].real() / _measurement_sweeps << " "
           << std::setw(10) << _tof_greenfunction[site].imag() / _measurement_sweeps << std::endl;

        //next coordinate
        int dir = _dim - 1;
        while (dir >= 0 && coordinate[dir] == _Nl - 1)
            coordinate[dir] = 0, --dir;
        if (dir == -1)
            break;
        ++coordinate[dir];
    }
    os.flags(oldFlags);
}

#endif
