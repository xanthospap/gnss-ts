#ifndef __TIMESERIES__
#define __TIMESERIES__

#include <vector>
#include <type_traits>
#include "datetime.hpp"
#include "flags.hpp"

namespace ts {

struct cmp_record {
    double    val_;
    double    sigma_;
    flag_type flag_;
};

class ts_cmp {

using std::vector<cmp_record> = recvec;

public:
    ts_cmp() noexcept : {}

    ts_cmp(const ts_cmp& rhs)
        noexcept( std::is_nothrow_copy_constructible<recvec>::value )
        : records_{rhs.records_} {}

    ts_cmp& operator=(const ts_cmp& rhs)
        noexcept( std::is_nothrow_copy_assignable<recvec>::value )
    {
        if ( this != &rhs ) {
            epochs_  = rhs.epochs_;
        }
        return *this;
    }
    
    ts_cmp(ts_cmp&& rhs)
        noexcept( std::is_nothrow_move_constructible<recvec>::value )
        : records_{std::move(rhs.records_) } {}

    ts_cmp& operator=(ts_cmp&& rhs)
        noexcept( std::is_nothrow_move_assignable<recvec>::value )
    {
        records_ = std::move(rhs.records_);
        return *this;
    }
    
    inline
    dtvec*& dt_ptr() 
    noexcept
    {
        return epochs_;
    }

    inline
    cmp_record& at(std::size_t idx)
    {
        return records_[idx];
    }

    inline
    void push_back(const cmp_record& rec)
    {
        records_.push_back( rec );
    }

    inline
    void push_back(double v, double s = 1.0e0, flag_type f = 0)
    {
        records_.emplace_back( v, s, f );
    }

    /*
    std::size_t
    average( double& average, double& sigma_2 )
    {
        double K, tmp;
        double average_n  { records_[0].val_   };
        double sigma_n    { records_[0].sigma_ };
        double n          { 1 };
        std::size_t n_max { records_.size() };
        for (int i=1; i<n_max; ++i) {
            K   = 1.0 / (n + 1.0);
            tmp = K * ( records_[i] - average_n );
            average_n += tmp;
            sigma_n   += (tmp * tmp);
            sigma_n   *= ( 1.0 - K );
            n = i;
        }
        average = average_n;
        sigma_2 = sigma_n;
    }
    */
    
private:
    recvec  records_;
};

} // end namespace

#endif
