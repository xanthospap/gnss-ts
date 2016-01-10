#ifndef __DATETIME_DUMMY__
#define __DATETIME_DUMMY__

namespace ngpt {

constexpr long     _JAN61980_  { 44244     };
constexpr long     _JAN11901_  { 15385     };
constexpr double   _SECPERDAY_ { 86400.0e0 };

const long LeapMonths[13] =   { 0,  31,  60,  91, 121, 152, 182,
                              213, 244, 274, 305, 335, 366 };

const long NormalMonths[13] = { 0,  31,  59,  90, 120, 151, 181,
                              212, 243, 273, 304, 334, 365 };

class datetime {

public:
    datetime() noexcept
        : mjd_(0), fday_(.0)
    {}

    datetime(int year, int month, int day, 
             int hour = 0, int minute = 0, double sec = 0)
    noexcept;

    datetime(const datetime&) noexcept = default;
    datetime(datetime&&)      noexcept = default;

    ~datetime() noexcept = default;

    bool operator==(const datetime& d) 
    const noexcept
    {
        return ( mjd_ == d.mjd_ && fday_ == d.fday_ );
    }

private:
    long   mjd_;
    double fday_;
};

} // end of namespace

#endif
