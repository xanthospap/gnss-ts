#ifndef __DTCALENDAR_NGPT__HPP__
#define __DTCALENDAR_NGPT__HPP__

#include "dtfund.hpp"
#ifdef DEBUG
    #include <iostream>
    #include <iomanip>  // std::setprecision
    #include <sstream>  // std::ostringstream
#endif

namespace ngpt {

/*
 * A datetime class. Holds (integral) days as MJD and fraction of day as any
 * of the is_of_sec_type class (i.e. seconds/milli/micro).
 */
template<class S,
        typename = std::enable_if_t<S::is_of_sec_type>
        >
    class datetime {
public:

    /// Only allow S parameter to be of sec type (seconds/milli/nano).
    static_assert( S::is_of_sec_type, "" );
    
    /// Default (zero) constructor.
    explicit constexpr datetime() noexcept : m_mjd(0), m_sec(0) {};
    
    /// Constructor from year, month, day of month and sec type.
    explicit constexpr datetime(year y, month m, day_of_month d, S s)
    : m_mjd{cal2mjd(y, m, d)}, m_sec{s}
    {}

    /// Constructor from year, month, day of month and any sec type (T),
    /// convertible to S.
    template<class T,
            typename = std::enable_if_t<T::is_of_sec_type>,
            typename = std::enable_if_t<
                std::is_same<S, decltype(static_cast<S>(T{}))>::value,
                bool
                >
            >
        explicit datetime(year y, month m, day_of_month d, T t)
        : m_mjd{cal2mjd(y, m, d)}, m_sec{S(t)}
    {}

    /// Constructor from year, month, day of month and fractional seconds.
    explicit
    datetime(year y, month m, day_of_month d, hours hr, minutes mn,
        double fsecs)
    : m_mjd{cal2mjd(y, m, d) }, m_sec{hr, mn, fsecs}
    {}

    /// Constructor from modified julian day, hours, minutes and 
    /// micro- or milli- or seconds.
    explicit
    datetime(modified_julian_day mjd, hours hr=hours(), minutes mn=minutes(),
        S sec=S())
    : m_mjd{mjd}, m_sec{hr, mn, sec}
    {}
    
    /// Constructor from year, month, day of month, hours, minutes and
    /// second type S.
    explicit
    datetime(year y, month m, day_of_month d, hours hr=hours(),
        minutes mn=minutes(), S sec=S())
    : m_mjd{cal2mjd(y, m, d)}, m_sec{hr, mn, sec}
    {}
    
    /// Constructor from year, month, day of month, hours, minutes and
    /// micro- or milli- or seconds (if T can be cast to S).
    template<class T,
            typename = std::enable_if_t<T::is_of_sec_type>
            >
    explicit
        datetime(year y, month m, day_of_month d, hours hr, minutes mn, T sec)
        : m_mjd{cal2mjd(y, m, d)}, m_sec{hr, mn, S(sec)}
    {}

    /// Get the MJDay.
    constexpr modified_julian_day  mjd() const noexcept { return m_mjd; }
    constexpr modified_julian_day& mjd() noexcept { return m_mjd; }

    #ifdef DEBUG
    /// Get the sec type:
    constexpr typename S::underlying_type secs() const noexcept
    { return m_sec.as_underlying_type(); }
    #endif

    /// Add any second type T, convertible to S.
    template<class T>
        constexpr void add_seconds(T t) noexcept
    { 
        m_sec += (S)t;
        if ( m_sec.more_than_day() ) {
            this->normalize();
        }
        return;
    }
    
    /// Subtract any second type T, convertible to S.
    template<class T>
        constexpr void remove_seconds(T t) noexcept
    { 
        m_sec -= (S)t;
        if ( m_sec < (S)0 ) {
            while ( m_sec < (S)0 ) {
                --m_mjd;
                m_sec += (S)S::max_in_day;
            }
        }
        return;
    }

    /// Return the difference of two datetimes as second type S.
    constexpr S delta_sec(const datetime& d) const noexcept
    {
        return (m_sec - d.m_sec) + mjd_sec_diff<S>(m_mjd, d.m_mjd);
    }

    /// Overload equality operator.
    constexpr bool operator==(const datetime& d) const noexcept
    { return m_mjd == d.m_mjd && m_sec == d.m_sec; }

    /// Overload ">" operator.
    constexpr bool operator>(const datetime& d) const noexcept
    { return m_mjd > d.m_mjd || (m_mjd == d.m_mjd && m_sec > d.m_sec); }
    
    /// Overload ">=" operator.
    constexpr bool operator>=(const datetime& d) const noexcept
    { return m_mjd > d.m_mjd || (m_mjd == d.m_mjd && m_sec >= d.m_sec); }
    
    /// Overload "<" operator.
    constexpr bool operator<(const datetime& d) const noexcept
    { return m_mjd < d.m_mjd || (m_mjd == d.m_mjd && m_sec < d.m_sec); }
    
    /// Overload "<=" operator.
    constexpr bool operator<=(const datetime& d) const noexcept
    { return m_mjd < d.m_mjd || (m_mjd == d.m_mjd && m_sec <= d.m_sec); }

    /// Reset the seconds/milli/nano after removing whole days.
    constexpr void normalize() noexcept
    {
        // TODO what happens if m_sec < 0 ??
        modified_julian_day add {static_cast<long>(m_sec.remove_days())};
        m_mjd += add;
        return;
    }

    /// Cast to double (i.e. fractional) Modified Julian Date.
    constexpr double as_mjd() const noexcept
    {
        return static_cast<double>(m_mjd.as_underlying_type())
            + m_sec.fractional_days();
    }

    /// Cast to year, month, day of month
    constexpr std::tuple<year, month, day_of_month>
    as_ymd() const noexcept
    {
        // TODO should i normalize the (this) date before calling to_ymd()
        year y;
        month m;
        day_of_month d;
        std::tie(y, m, d) = m_mjd.to_ymd();
        return std::make_tuple(y, m, d);
    }

    /// Cast to year, day_of_year
    constexpr std::tuple<year, day_of_year>
    as_ydoy() const noexcept
    {
        // TODO should i normalize the (this) date before calling to_ydoy()
        year y;
        day_of_year d;
        std::tie(y, d) = m_mjd.to_ydoy();
        return std::make_tuple(y, d);
    }

    /// Convert the *seconds to hours, minutes, seconds and S
    constexpr std::tuple<hours, minutes, seconds, long>
    as_hmsf() const noexcept { return m_sec.to_hmsf(); }

#ifdef DEBUG
    std::string stringify() const
    {
        auto ymd { this->as_ymd() };
        auto hms { this->as_hmsf() };
        S st { std::get<3>(hms) };
        double fsec { st.as_underlying_type() * S::template sec_ifactor<double>() };
        std::ostringstream out;
        out << std::fixed << std::setprecision(9) << fsec;

        return std::string {
                     std::to_string( std::get<0>(ymd).as_underlying_type() )
             + "/" + std::to_string( std::get<1>(ymd).as_underlying_type() )
             + "/" + std::to_string( std::get<2>(ymd).as_underlying_type() )
             + " " + std::to_string( std::get<0>(hms).as_underlying_type() )
             + ":" + std::to_string( std::get<1>(hms).as_underlying_type() )
             + ":" + std::to_string( std::get<2>(hms).as_underlying_type() )
             + "+" + out.str()
        };
    }
#endif

private:
    modified_julian_day m_mjd;  ///< Modified Julian Day
    S                   m_sec;  ///< Fraction of day in milli/nano/seconds

}; // end class datetime

} // end namespace

#endif // define DATETIME
