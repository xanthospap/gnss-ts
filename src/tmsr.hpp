#ifndef __TIMESERIES__
#define __TIMESERIES__

#include <vector>
#include <type_traits>
#include "tscmp.hpp"

namespace ts {

class time_series {

using std::vector<ngpt::datetime> dtvec;
using ts::ts_cmp;

public:
    time_series() {};
    time_series(const time_series& ts)         = default;
    time_series(time_series&& ts)              = default;
    time_series& operator=(const time_series&) = default;
    time_series& operator=(time_series&&)      = default;

    inline std::size_t size() const noexcept
    {
        return epochs_.size();
    }

    std::size_t read_from_cts(const char*);

private:
    dtvec  epochs_;
    ts_cmp xcmp_;
    ts_cmp ycmp_;
    ts_cmp zcmp_;
};

} // end namespace

#endif
