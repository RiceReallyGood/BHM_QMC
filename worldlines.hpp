#ifndef WORLDLINES_HPP
#define WORLDLINES_HPP

#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>

// ==================================================
// _node struct
// ==================================================
struct _Node
{
    _Node() = default;
    _Node(unsigned int partnersite_, double time_ = 0, unsigned short state_ = 0) : _partnersite(partnersite_), _time(time_), _state(state_) {}
    _Node(const _Node &node_) : _partnersite(node_._partnersite), _time(node_._time), _state(node_._state){};
    unsigned int _partnersite;
    double _time;
    unsigned short _state;
};

// ==================================================
// kink class
// ==================================================
class kink
{
public:
    kink() : _node(new _Node) {}
    kink(unsigned int partnersite_, double time_ = 0, unsigned short state_ = 0) : _node(new _Node(partnersite_, time_, state_)) {}
    kink(const kink &k) : _node(new _Node(*k._node)) {}
    kink(kink &&k) noexcept : _node(k._node) { k._node = nullptr; }
    kink &operator=(const kink &rhs)
    {
        *_node = *rhs._node;
        return *this;
    }
    kink &operator=(kink &&rhs)
    {
        if (this != &rhs)
        {
            delete _node;
            _node = rhs._node;
            rhs._node = nullptr;
        }
        return *this;
    }
    ~kink() { delete _node; }

    unsigned int partnersite() const { return _node->_partnersite; }
    double time() const { return _node->_time; }
    unsigned short state() const { return _node->_state; }

    void set_partnersite(unsigned int partnersite_) { _node->_partnersite = partnersite_; }
    void set_time(double time_) { _node->_time = time_; }
    void set_state(unsigned short state_) { _node->_state = state_; }

    void state_increment() { ++(_node->_state); }
    void state_decrement() { --(_node->_state); }

    bool operator>(const double &time_) const { return _node->_time > time_; }
    bool operator<(const double &time_) const { return _node->_time < time_; }
    bool operator>=(const double &time_) const { return _node->_time >= time_; }
    bool operator<=(const double &time_) const { return _node->_time <= time_; }

    kink &operator+=(const double &delta_time_)
    {
        _node->_time += delta_time_;
        return *this;
    }

    kink &operator-=(const double &delta_time_)
    {
        _node->_time -= delta_time_;
        return *this;
    }

    friend void swap_state(kink &obj1_, kink &obj2_) { std::swap(obj1_._node->_state, obj2_._node->_state); }
    friend std::ostream &operator<<(std::ostream &out, kink const &obj)
    {
        out << "\t" << obj.partnersite() << "\t" << obj.time() << "\t" << obj.state();
        return out;
    }

private:
    _Node *_node;
};

// ==================================================
// worldlines class
// ==================================================

class worldlines
{
public:
    typedef std::vector<kink> line;
    typedef std::vector<line> lines;
    typedef std::pair<lines::iterator, line::iterator> location_type;

    worldlines() {}
    worldlines(unsigned int num_sites_)
    {
        _worldlines.resize(num_sites_);
        for (unsigned int site = 0; site < num_sites_; ++site)
            _worldlines[site].push_back(kink(site));
    }

    unsigned int num_sites() const { return _worldlines.size(); }
    unsigned int num_kinks(unsigned int site_) const { return _worldlines[site_].size() - 1; }
    unsigned int total_interactions() const;

    unsigned int site_state(unsigned int site_) const { return _worldlines[site_][0].state(); }

    location_type location(unsigned int site_, double time_)
    {
        return std::make_pair(_worldlines.begin() + site_, std::lower_bound(_worldlines[site_].begin() + 1, _worldlines[site_].end(), time_));
    }

    unsigned short state_before(const location_type &location_) const { return (location_.second - 1)->state(); }
    unsigned short state(const location_type &location_) const { return (location_.second == location_.first->end() ? location_.first->begin()->state() : location_.second->state()); }

    bool location_is_kink_unoccupied(location_type const &location_, double time_) const { return (location_.second == location_.first->end() || location_.second->time() != time_); }
    int net_number_of_directed_hops(unsigned int site_, unsigned int partnersite_) const;

    bool is_valid(unsigned short max_state);
    void print() const;

private:
    lines _worldlines;
};

unsigned int worldlines::total_interactions() const
{
    unsigned int ti = 0;
    for (int site = 0; site < num_sites(); ++site)
        ti += num_kinks(site);
    return ti / 2;
}

int worldlines::net_number_of_directed_hops(unsigned int site_, unsigned int partnersite_) const
{
    int net_number = 0;
    line::const_iterator it = _worldlines[site_].begin();
    ++it;
    for (; it != _worldlines[site_].end(); ++it)
        if (it->partnersite() == partnersite_)
            (it->state() > (it - 1)->state()) ? ++net_number : --net_number;
    return net_number;
}

void worldlines::print() const
{
    for (unsigned int site = 0; site < num_sites(); site++)
    {
        std::cout << "site " << site << ": ";
        for (kink _kink : _worldlines[site])
            std::cout << _kink.time() << "->";
        std::cout << std::endl;
    }
}

bool worldlines::is_valid(unsigned short max_state)
{
    bool valid = true;

    //testing vertex state
    for (unsigned int site = 0; site < _worldlines.size(); ++site)
    {
        if (_worldlines[site][0].state() > max_state)
            valid = false;
        for (unsigned int i = 1; i < _worldlines[site].size(); ++i)
        {
            if (_worldlines[site][i].state() > max_state)
                valid = false;
            short this_state_increment = _worldlines[site][i].state() - _worldlines[site][i - 1].state();
            if ((this_state_increment != 1) && (this_state_increment != -1))
                valid = false;
        }
        if (!valid)
        {
            std::cout << "\nError: testing vertex state fails...\n";
            std::cout << "site " << site << " : ";
            for (unsigned int i = 0; i < _worldlines[site].size(); ++i)
                std::cout << _worldlines[site][i].state() << "  ";
            std::cout << "\n";
            return false;
        }
    }

    // testing vertex time
    for (unsigned int site = 0; site < _worldlines.size(); ++site)
    {
        if (_worldlines[site][0].time() != 0)
            valid = false;
        for (unsigned int i = 1; i < _worldlines[site].size(); ++i)
        {
            if (_worldlines[site][i].time() < 0 || _worldlines[site][i].time() > 1)
                valid = false;
            if (_worldlines[site][i].time() <= _worldlines[site][i - 1].time())
                valid = false;
        }
        if (!valid)
        {
            std::cout << "\nError: testing vertex time fails...\n";
            std::cout << "site " << site << " : ";
            for (unsigned int i = 0; i < _worldlines[site].size(); ++i)
                std::cout << _worldlines[site][i].time() << "  ";
            std::cout << "\n";
            return false;
        }
    }

    // testing vertex siteindicator
    for (unsigned int site = 0; site < _worldlines.size(); ++site)
    {
        if (_worldlines[site][0].partnersite() != site)
            valid = false;
        for (unsigned int i = 1; i < _worldlines[site].size(); ++i)
            if (_worldlines[site][i].partnersite() >= _worldlines.size())
                valid = false;

        if (!valid)
        {
            std::cout << "\nError: testing vertex siteindicator fails...\n";
            std::cout << "site " << site << " : ";
            for (unsigned int i = 0; i < _worldlines[site].size(); ++i)
                std::cout << _worldlines[site][i].partnersite() << "  ";
            std::cout << "\n";
            return false;
        }
    }

    // testing vertex pairing
    for (unsigned int site = 0; site < _worldlines.size(); ++site)
    {
        if (_worldlines[site].size() > 1)
        {
            for (unsigned int i = 1; i < _worldlines[site].size(); ++i)
            {
                kink linkedto(*(location(_worldlines[site][i].partnersite(), _worldlines[site][i].time()).second));

                if (linkedto.time() != _worldlines[site][i].time())
                    valid = false;

                if (linkedto.partnersite() != site)
                    valid = false;
            }
        }

        if (!valid)
        {
            std::cout << "\nError: testing vertex paring fails...\n";
            std::cout << "site " << site << "\n";
            for (unsigned int i = 1; i < _worldlines[site].size(); ++i)
            {
                std::cout << "kink : " << _worldlines[site][i] << " , linkedto : " << *(location(_worldlines[site][i].partnersite(), _worldlines[site][i].time()).second) << "\n";
            }
            return false;
        }
    }

    return valid;
}

// ==================================================
// wormpair class
// ==================================================

class wormpair
{
public:
    typedef worldlines::lines::iterator linesiterator;
    typedef worldlines::line::iterator lineiterator;
    typedef worldlines::location_type location_type;

    wormpair() {}
    wormpair(location_type location_, const kink &kink_, bool forward_, bool creation_);

    unsigned int wormhead_site() const { return _wormhead.partnersite(); }
    unsigned int wormtail_site() const { return _wormtail.partnersite(); }

    double wormhead_time() const { return _wormhead.time(); }
    double wormtail_time() const { return _wormtail.time(); }

    unsigned short wormhead_state() const { return _wormhead.state(); }
    unsigned short wormtail_state() const { return _wormtail.state(); }
    unsigned short wormpair_state() const { return _wormpair_state; }
    unsigned short state_before_head() const { return _creation ? _wormhead.state() - 1 : _wormhead.state() + 1; }

    void set_next_iterator()
    {
        _next = _forward ? (_location.second == _location.first->end() ? _location.first->begin() + 1 : _location.second)
                         : (_location.second == _location.first->begin() + 1 ? _location.first->end() - 1 : _location.second - 1);
    }

    bool forward() const { return _forward; }
    bool creation() const { return _creation; }
    bool increasing() const { return (_forward && !_creation) || (!_forward && _creation); }
    static bool increasing(bool forward_, bool creation_) { return (forward_ && !creation_) || (!forward_ && creation_); }

    unsigned int next_partnersite() const { return _next->partnersite(); }
    double next_time() const { return _next->time(); }
    unsigned short next_state() const { return _next->state(); }

    static double mod_one(double value_) { return value_ - std::floor(value_); }
    double time2wormtail() const { return _forward ? wormpair::mod_one(_wormtail.time() - _wormhead.time()) : wormpair::mod_one(_wormhead.time() - _wormtail.time()); }
    double time2next() const { return _forward ? wormpair::mod_one(_next->time() - _wormhead.time()) : wormpair::mod_one(_wormhead.time() - _next->time()); }

    bool wormhead_reaches_wormtail() const { return _next->partnersite() == wormhead_site(); }
    bool wormhead_is_same_type_as_next() const
    {
        return _forward ? ((_creation && _next->state() == wormhead_state() + 1) || (!_creation && _next->state() + 1 == wormhead_state()))
                        : ((_creation && _next->state() == (_next - 1)->state() + 1) || (!_creation && _next->state() + 1 == (_next - 1)->state()));
    }

    inline void wormhead_turns_around()
    {
        _forward = !_forward;
        set_next_iterator();
    };

    inline void wormhead_moves_to_new_time(double time_, bool winding_over_time_ = false);
    inline void wormhead_inserts_vertex_and_jumps_to_new_site(const location_type &targetlocation_);
    inline void wormhead_deletes_vertex_and_jumps_to_new_site(const location_type &sourcelocation_);
    inline void wormhead_relinks_vertex_and_jumps_to_new_site(const location_type &sourcelocation_, const location_type &targetlocation_);
    inline void wormhead_crosses_vertex();
    inline void wormhead_annihilates_wormtail() { _location.first->erase(_next); }

private:
    unsigned short _wormpair_state;

    kink _wormtail;
    kink _wormhead;
    bool _forward;
    bool _creation;

    location_type _location;
    lineiterator _next;
};

// ==================================================
// wormpair member functions
// ==================================================

wormpair::wormpair(location_type location_, const kink &kink_, bool forward_, bool creation_)
    : _location(location_),
      _wormhead(kink_),
      _wormtail(kink_),
      _forward(forward_),
      _creation(creation_)
{
    unsigned short state_between = increasing() ? kink_.state() + 1 : kink_.state() - 1;
    _forward ? _wormtail.set_state(state_between) : _wormhead.set_state(state_between);

    _wormhead += forward_ ? std::numeric_limits<double>::epsilon() : -std::numeric_limits<double>::epsilon();

    _wormpair_state = (_creation ? wormhead_state() : wormtail_state());

    _location.second = location_.first->insert(location_.second, _wormtail);
    if (_forward)
        ++_location.second;

    set_next_iterator();
}

inline void wormpair::wormhead_moves_to_new_time(double time_, bool winding_over_time_)
{
    //std::cout << "wormhead moves " << (_forward ? "forward " : "backward ") << "from " << wormhead_time() << " to " << time_ << std::endl;
    if (winding_over_time_)
    {
        increasing() ? _location.first->begin()->state_increment() : _location.first->begin()->state_decrement();
        _location.second = _forward ? _location.first->begin() + 1 : _location.first->end();
        //std::cout << "winding site = " << _location.first->begin()->partnersite() << " state = " << _location.first->begin()->state() << std::endl;
    }
    _wormhead.set_time(time_);
}

inline void wormpair::wormhead_inserts_vertex_and_jumps_to_new_site(const location_type &targetlocation_)
{
    //std::cout << "insert a vertex" << (forward() ? "forward" : "backward") << " from " << wormhead_site() << " to " << targetlocation_.first->begin()->partnersite() << std::endl;
    unsigned int target_site = targetlocation_.first->begin()->partnersite();
    unsigned short target_state = (targetlocation_.second - 1)->state();
    unsigned short newstate = increasing() ? target_state + 1 : target_state - 1;
    _location.first->insert(_location.second, kink(target_site, wormhead_time(), wormhead_state()));
    _location = targetlocation_;
    _location.second = targetlocation_.first->insert(targetlocation_.second, kink(wormhead_site(), wormhead_time(), target_state));
    if (_forward)
    {
        _location.second->set_state(newstate);
        _wormhead = kink(target_site, wormhead_time() + std::numeric_limits<double>::epsilon(), target_state);
        ++_location.second;
    }
    else
    {
        _wormhead = kink(target_site, wormhead_time() - std::numeric_limits<double>::epsilon(), newstate);
    }
    set_next_iterator();
}

inline void wormpair::wormhead_deletes_vertex_and_jumps_to_new_site(const location_type &sourcelocation_)
{
    //std::cout << "delete a vertex from site " << wormhead_site() << " to site " << sourcelocation_.first->begin()->partnersite() << "and jump to site " << sourcelocation_.first->begin()->partnersite() << std::endl;
    _location.first->erase(_next);

    _wormhead = *sourcelocation_.second;
    _wormhead.set_partnersite(sourcelocation_.first->begin()->partnersite());

    _location = sourcelocation_;
    _location.second = _location.first->erase(_location.second);

    set_next_iterator();
}

inline void wormpair::wormhead_relinks_vertex_and_jumps_to_new_site(const location_type &sourcelocation_, const location_type &targetlocation_)
{
    //std::cout << "wormhead relinks vertex from " << wormhead_site() << "and " << sourcelocation_.first->begin()->partnersite()
    //          << " to " << sourcelocation_.first->begin()->partnersite() << " and " << targetlocation_.first->begin()->partnersite() << std::endl;
    _location.first->erase(_next);

    sourcelocation_.second->set_partnersite(targetlocation_.first->begin()->partnersite());

    _forward = !_forward;
    unsigned short targer_state = (targetlocation_.second - 1)->state();
    unsigned short newstate = increasing() ? targer_state + 1 : targer_state - 1;

    _location = targetlocation_;
    if (_forward)
    {
        _location.second = targetlocation_.first->insert(targetlocation_.second, kink(sourcelocation_.first->begin()->partnersite(), sourcelocation_.second->time(), newstate));
        _wormhead = kink(_location.first->begin()->partnersite(), sourcelocation_.second->time() + std::numeric_limits<double>::epsilon(), targer_state);
        ++_location.second;
    }
    else
    {
        _location.second = targetlocation_.first->insert(targetlocation_.second, kink(sourcelocation_.first->begin()->partnersite(), sourcelocation_.second->time(), targer_state));
        _wormhead = kink(_location.first->begin()->partnersite(), sourcelocation_.second->time() - std::numeric_limits<double>::epsilon(), newstate);
    }

    set_next_iterator();
}

inline void wormpair::wormhead_crosses_vertex()
{
    //std::cout << "wormhead crosses a vertex" << std::endl;
    swap_state(_wormhead, *_next);
    if (_forward)
    {
        _wormhead.set_time(next_time() + std::numeric_limits<double>::epsilon());
        ++_location.second;
    }
    else
    {
        _wormhead.set_time(next_time() - std::numeric_limits<double>::epsilon());
        --_location.second;
    }

    set_next_iterator();
}

#endif