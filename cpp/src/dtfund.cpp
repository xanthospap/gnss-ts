#include "dtfund.hpp"

/// Definition for static month array (short names).
constexpr const char* ngpt::month::short_names[];

/// Definition for static month array (long names).
constexpr const char* ngpt::month::long_names[];

/// Number of days past at the end of non-leap and leap years.
constexpr static long month_day[2][13] = {
    {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365},
    {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366}
};

/// Transform a year, month, day of month to modified_julian_day.
ngpt::modified_julian_day
ngpt::cal2mjd(ngpt::year y, ngpt::month m, ngpt::day_of_month d)
{
    long mjd {cal2mjd( y.as_underlying_type(), m.as_underlying_type(), 
        d.as_underlying_type()) };
    return ngpt::modified_julian_day{mjd};
}

/// Cast a modified_julian_day to year, day_of_year
std::tuple<ngpt::year, ngpt::day_of_year>
ngpt::modified_julian_day::to_ydoy() const noexcept
{
    long days_fr_jan1_1901 { m_mjd - ngpt::jan11901 };
    long num_four_yrs      { days_fr_jan1_1901/1461L };
    long years_so_far      { 1901L + 4*num_four_yrs };
    long days_left         { days_fr_jan1_1901 - 1461*num_four_yrs };
    long delta_yrs         { days_left/365 - days_left/1460 };
    day_of_year yday {static_cast<day_of_year::underlying_type>
                                (days_left - 365*delta_yrs + 1)};
    ngpt::year y {static_cast<year::underlying_type>(years_so_far + delta_yrs)};

    return std::make_tuple(y, yday);
}

/// Cast a modified_julian_day to year, month, day_of_month
std::tuple<ngpt::year, ngpt::month, ngpt::day_of_month>
ngpt::modified_julian_day::to_ymd() const noexcept
{
    ngpt::day_of_year doy;
    ngpt::year y;
    std::tie(y, doy) = this->to_ydoy();
    long yday  { static_cast<long>(doy.as_underlying_type()) };
    long leap  { ((y.as_underlying_type()%4L) == 0) };
    long guess { static_cast<long>(yday*0.032) };
    long more  { ((yday-month_day[leap][guess+1]) > 0) };

    ngpt::month mon {static_cast<ngpt::month::underlying_type>(guess + more + 1)};
    ngpt::day_of_month dom {static_cast<ngpt::day_of_month::underlying_type>
        (yday-month_day[leap][guess+more]) };

    return std::make_tuple(y, mon, dom);
}

/// Convert a pair of Year, Day of year to a Modified Julian Day
/// Reference: http://www.ngs.noaa.gov/gps-toolbox/bwr-c.txt :: ydhms_to_mjd
ngpt::modified_julian_day
ngpt::ydoy2ymd(ngpt::year yr, ngpt::day_of_year doy) noexcept
{
    long iyr { static_cast<long>(yr.as_underlying_type()) };
    long idy { static_cast<long>(doy.as_underlying_type()) };
    return modified_julian_day {((iyr-1901)/4)*1461 + ((iyr-1901)%4)*365 +
        idy - 1 + ngpt::jan11901};
}

/// Calendar date (i.e. year-month-day) to Modified Julian Date.
///
/// \return    The Modified Julian Date (as \c long).
///
/// \throw     A runtime_error if the month and/or day is invalid.
///
/// Reference: iauCal2jd
///
long
ngpt::cal2mjd(int iy, int im, int id)
{
    // Month lengths in days
    static constexpr int mtab[] = 
        {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    // Validate month
    if ( im < 1 || im > 12 ) {
        throw std::out_of_range("ngpt::cal2mjd -> Invalid Month.");
    }

    // If February in a leap year, 1, otherwise 0
    int ly ( (im == 2) && /*ngpt::*/is_leap(iy) );

    // Validate day, taking into account leap years
    if ( (id < 1) || (id > (mtab[im-1] + ly))) {
        throw std::out_of_range("ngpt::cal2mjd -> Invalid Day of Month.");
    }

    // Compute mjd
    int  my    { (im-14) / 12 };
    long iypmy { static_cast<long>(iy + my) };

        
    return  (1461L * (iypmy + 4800L)) / 4L
            + (367L * static_cast<long>(im - 2 - 12 * my)) / 12L
            - (3L * ((iypmy + 4900L) / 100L)) / 4L
            + static_cast<long>(id) - 2432076L;
}
