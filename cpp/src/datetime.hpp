#ifndef __DATETIME_NGPT__HPP__
#define __DATETIME_NGPT__HPP__

#include <ostream>
#include <iomanip>
#include <cstdio>
#include <limits>
#include <cassert>
#include <cmath>
#include <tuple>
#include <string>
#ifdef DEBUG
    #include <iostream>
#endif

namespace ngpt {

/// Check if long is big enough to hold two days in microseconds.
static_assert( 86400L  * 1000000L * 2 < std::numeric_limits<long>::max(),
    "FUCK! Long is not big enough to hold two days in microseconds" );

/// Jan 1st 1980
constexpr long jan61980 { 44244L };

constexpr long jan11901 { 15385L };

/// Seconds per day.
constexpr double sec_per_day { 86400.0e0 };

/// Days per Julian year.
constexpr double days_in_julian_year { 365.25e0 };

/// Days per Julian century.
constexpr double days_in_julian_cent { 36525.0e0 };

/// Reference epoch (J2000.0), Julian Date.
constexpr double j2000_jd { 2451545.0e0 };

/// Reference epoch (J2000.0), Modified Julian Date.
constexpr double j2000_mjd { 51544.5e0 };

/// Julian Date of Modified Julian Date zero.
constexpr double mjd0_jd { 2400000.5e0 };

/// TT minus TAI (s)
constexpr double tt_minus_tai { 32.184e0 };

/// Time-Scales
/// enum class time_scale : char { tai, tt, utc, ut1 };

/// Calendar date to MJD.
long cal2mjd(int, int, int);

/// MJD to calendar date.
void mjd2cal(long, int&, int&, int&) noexcept;

/// Convert hours, minutes, seconds into fractional days.
double hms2fd(int, int, double) noexcept;

/// Decompose fractional days to hours, minutes, seconds and fractional seconds
/// with a given precision.
void fd2hms(double, int, int ihmsf[4]);

/// Return true if given year is leap.
inline bool is_leap(int iy) noexcept
{ 
    return !(iy%4) && (iy%100 || !(iy%400));
}

/// Forward declerations
class year;
class month;
class day_of_month;
class day_of_year;
class day;
class modified_julian_day;
class julian_day;
class seconds;
class milliseconds;
class microseconds;

/// A wrapper class for years.
class year {

public:
    /// Years are represented as integers.
    typedef int underlying_type;

    /// Constructor.
    explicit constexpr year (underlying_type i=0) noexcept : m_year(i) {};

    /// Get the underlying int.
    constexpr underlying_type as_underlying_type() const noexcept
    { return m_year; }

private:
    /// The year as underlying type.
    underlying_type m_year;

}; // class year

/// A wrapper class for months.
class month {

public:
    /// Months are represented as int.
    typedef int underlying_type;

    /// Constructor.
    explicit constexpr month (underlying_type i=1) noexcept : m_month(i) {};

    /// Get the underlying int.
    constexpr underlying_type as_underlying_type() const noexcept
    { return m_month; }

    const char* short_name() const noexcept { return short_names[m_month-1]; }
    const char* long_name()  const noexcept { return long_names[m_month-1];  }

private:
    /// Decleration of short month names. Note that we do need a definition
    /// in the .cpp file.
    constexpr static const char* short_names[] = {
        "Jan", "Feb", "Mar", "Apr", "May", "Jun",
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
    };
    
    /// Decleration of long month names. Note that we do need a definition
    /// in the .cpp file.
    constexpr static const char* long_names[] = {
        "January", "February", "March", "April", "May", "June",
        "July", "August", "September", "October", "November", "December"
    };

    /// The month as underlying_type.
    underlying_type m_month;

}; // class month

/*
/// A wrapper class for days (in general!).
class day {

public:
    /// Days are represented as ints.
    typedef int underlying_type;
    
    /// Constructor.
    explicit constexpr day(underlying_type i=0) noexcept : m_day(i) {};
    
    /// Get the underlying int.
    constexpr underlying_type as_underlying_type() const noexcept
    { return m_day; }

private:
    /// The day as underlying_type.
    underlying_type m_day;
};
*/

/// A wrapper class for day of month.
class day_of_month {

public:
    /// Days are represented as ints.
    typedef int underlying_type;
    
    /// Constructor.
    explicit constexpr day_of_month(underlying_type i=0) noexcept : m_dom(i) {};
    
    /// Get the underlying int.
    constexpr underlying_type as_underlying_type() const noexcept
    { return m_dom; }

private:
    /// The day of month as underlying_type.
    underlying_type m_dom;   

}; // class day_of_month

/// A wrapper class for Modified Julian Day (only the integer part is 
/// considered).
class modified_julian_day {
public:
    /// MJDs are represented as long ints.
    typedef long underlying_type;
    
    /// Constructor.
    explicit constexpr modified_julian_day(underlying_type i=0) noexcept
    : m_mjd(i) 
    {};
    
    /// Get the underlying long int.
    constexpr underlying_type as_underlying_type() const noexcept
    { return m_mjd; }
    
    /// Define addition (between MJDs).
    constexpr void operator+=(const modified_julian_day& d) noexcept
    { m_mjd += d.m_mjd; }
    
    /// Operator - (subtraction).
    constexpr modified_julian_day operator-(const modified_julian_day& mjd)
    const noexcept
    { return modified_julian_day{m_mjd-mjd.m_mjd}; }

    /// Operator == (equality).
    constexpr bool operator==(const modified_julian_day& d) const noexcept
    { return m_mjd == d.m_mjd; }

    /// Operator > (greater than).
    constexpr bool operator>(const modified_julian_day& d) const noexcept
    { return m_mjd > d.m_mjd; }

    /// Operator >= (greater or equal to).
    constexpr bool operator>=(const modified_julian_day& d) const noexcept
    { return m_mjd >= d.m_mjd; }

    /// Operator < (less than).
    constexpr bool operator<(const modified_julian_day& d) const noexcept
    { return m_mjd < d.m_mjd; }

    /// Operator <= (less or equal to).
    constexpr bool operator<=(const modified_julian_day& d) const noexcept
    { return m_mjd <= d.m_mjd; }

    /// Operator -- (minus 1 mjd).
    constexpr modified_julian_day& operator--() noexcept
    {
        --m_mjd;
        return *this;
    }

    /// Cast to year, day_of_year
    year to_ydoy(day_of_year&) const noexcept;
    
    /// Cast to year, month, day_of_month
    year to_ymd(month&, day_of_month&) const noexcept;
    
private:
    /// The modified julian day as underlying type.
    underlying_type m_mjd;

}; // class modified_julian_day

/// A wrapper class for day of year
class day_of_year {

public:
    /// Day of year represented as int.
    typedef int underlying_type;
    
    /// Constructor.
    explicit constexpr day_of_year(underlying_type i=0) noexcept : m_doy(i) {};
    
    /// Cast to underlying type
    constexpr underlying_type as_underlying_type() const noexcept
    { return m_doy; }

private:
    /// The day of year as underlying type.
    underlying_type m_doy;   
};

/// A wrapper class for hours.
class hours {
public:
    /// Hours are represented by ints.
    typedef int underlying_type;
    
    /// Constructor
    explicit constexpr hours(underlying_type i=0) noexcept : m_hours(i) {};
    
    /// Pass the underlying type
    constexpr underlying_type as_underlying_type() const noexcept
    { return m_hours; }

private:
    /// The hours as underlying type.
    underlying_type m_hours;

}; // class hours

/// A wrapper class for minutes.
class minutes {
public:
    /// Minutes are represented by ints
    typedef int underlying_type;
    
    /// Constructor
    explicit constexpr minutes(underlying_type i=0) noexcept : m_min(i) {};

    /// Pass the underlying type
    constexpr underlying_type as_underlying_type() const noexcept
    { return m_min; }

private:
    /// The minutes as underlying type.
    underlying_type m_min;

}; // class minutes

/// A wrapper class for seconds.
class seconds {
public:
    /// Seconds are represented as long ints.
    typedef long underlying_type;
    
    /// Seconds is a subdivision of seconds.
    static constexpr bool is_of_sec_type { true };
    
    /// Max seconds in day.
    static constexpr underlying_type max_in_day { 86400L };
/*
    template<typename T>
    static constexpr T sec_factor() noexcept
    { return static_cast<T>(1); }
    
    template<typename T>
    static constexpr T sec_ifactor() noexcept
    { return static_cast<T>(1); }
*/  
    /// Constructor.
    explicit constexpr seconds(underlying_type i=0L) noexcept : m_sec(i) {};

    /// Constructor from hours, minutes, seconds.
    explicit constexpr seconds(hours h, minutes m, seconds c) noexcept
    : m_sec {  c.as_underlying_type()
             + m.as_underlying_type()*60L
             + h.as_underlying_type()*3600L}
    {}

    /// Constructor from hours, minutes, fractional seconds.
    explicit constexpr seconds(hours h, minutes m, double fs) noexcept
    : m_sec{  static_cast<long>(fs)
            + m.as_underlying_type()*60L
            + h.as_underlying_type()*3600L}
    {}
    
    /// Addition operator between seconds.
    constexpr void operator+=(const seconds& sc) noexcept { m_sec+=sc.m_sec; }

    /// Subtraction operator between seconds.
    constexpr void operator-=(const seconds& sc) noexcept { m_sec-=sc.m_sec; }
    
    /// Overload - operator (subtraction).
    constexpr seconds operator-(const seconds& n) const noexcept
    { return seconds{m_sec - n.m_sec}; }

    /// Overload + operator (addition).
    constexpr seconds operator+(const seconds& sec) const noexcept
    { return seconds{m_sec+sec.m_sec}; }

    /// Overload operator '/'
    constexpr seconds operator/(const seconds& n) const noexcept
    { return seconds{m_sec/n.m_sec}; }
  
    /// Equality operator.  
    constexpr bool operator==(const seconds& d) const noexcept
    { return m_sec == d.m_sec; }

    /// Greater than operator.
    constexpr bool operator>(const seconds& d) const noexcept
    { return m_sec > d.m_sec; }

    /// Greater or equal to operator.
    constexpr bool operator>=(const seconds& d) const noexcept
    { return m_sec >= d.m_sec; }
    
    /// Less than operator.
    constexpr bool operator<(const seconds& d) const noexcept
    { return m_sec < d.m_sec; }

    /// Less or equal to operator.
    constexpr bool operator<=(const seconds& d) const noexcept
    { return m_sec <= d.m_sec; }
    
    /// Do the secods sum up to more than one day?
    constexpr bool more_than_day() const noexcept { return m_sec>max_in_day; }
    
    /// Get the underlying type numeric.
    constexpr underlying_type as_underlying_type() const noexcept
    { return m_sec; }
    
    /// If the seconds sum up to more (or equal to) one day, remove the integer
    /// days (and return them); reset the seconds to seconds of the new day.
    constexpr int remove_days() noexcept
    {
        int d { static_cast<int>(m_sec/max_in_day) };
        m_sec%=max_in_day;
        return d;
    }
    
    /// If the seconds sum up to more (or equal to) one day, return the 
    /// (integral) number of days.
    constexpr int to_days() const noexcept
    {
        return static_cast<int>(m_sec/max_in_day);
    }
    
    /// Interpret the seconds as fractional days.
    constexpr double fractional_days() const noexcept
    {
        return static_cast<double>(m_sec)/static_cast<double>(max_in_day);
    }
    
    /// Cast to double (i.e. fractional seconds).
    constexpr double to_fractional_seconds() const noexcept
    { return static_cast<double>(m_sec); }

    /// Translate to hours, minutes, seconds and fractional seconds
    constexpr std::tuple<hours, minutes, seconds, double>
    to_hmsf() const noexcept
    {
        return std::make_tuple(hours  {static_cast<int>(m_sec/3600L)},
                               minutes{static_cast<int>((m_sec%3600L)/60L)},
                               seconds{(m_sec%3600L)%60L},
                               0.0e0);
    }

private:
    /// Cast to any arithmetic type.
    template<typename T,
             typename=std::enable_if_t<std::is_arithmetic<T>::value>
             >
        constexpr T cast_to() const noexcept
    { return static_cast<T>(m_sec); }

    /// The seconds as underlying type.
    underlying_type m_sec;

}; // class seconds


/// A wrapper class for milliseconds (i.e. 10**-3 sec).
class milliseconds {

public:
    /// MilliSeconds are represented as long ints.
    typedef long underlying_type;
    
    /// MilliSeconds are a subdivision of seconds.
    static constexpr bool is_of_sec_type { true };
    
    /// Max milliseconds in one day.
    static constexpr long max_in_day { 86400L * 1000L };

/*
    template<typename T>
    static constexpr T sec_factor() noexcept
    { return static_cast<T>(1000); }
    template<typename T>
    static constexpr T sec_ifactor() noexcept
    { return ((static_cast<T>(1))/1000); }
*/
    
    /// Constructor.
    explicit constexpr milliseconds(underlying_type i=0L) noexcept
    : m_msec(i)
    {};
    
    explicit constexpr milliseconds(hours h, minutes m, milliseconds c) noexcept
    : m_msec { c.as_underlying_type()
        + m.as_underlying_type()*60L  *1000L
        + h.as_underlying_type()*3600L*1000L}
    {}
    
    /// Constructor from hours, minutes, fractional seconds
    explicit constexpr milliseconds(hours h, minutes m, double fs) noexcept
    : m_msec{ static_cast<long>(fs*1000.0e0)
        + (m.as_underlying_type()*60L
        + h.as_underlying_type()*3600L) * 1000L}
    {}
    
    /// Milliseconds can be cast to seconds (with a loss of precission).
    constexpr explicit operator seconds() const
    { return seconds{m_msec/1000L}; }
    
    /// Addition operator.
    constexpr milliseconds operator+(const milliseconds& sec) const noexcept
    { return milliseconds{m_msec+sec.m_msec}; }
    
    /// Addition operator.
    constexpr void operator+=(const milliseconds& ms) noexcept
    { m_msec+=ms.m_msec; }

    /// Subtraction operator.
    constexpr void operator-=(const milliseconds& ms) noexcept
    { m_msec-=ms.m_msec; }
    
    /// Subtraction operator.
    constexpr milliseconds operator-(const milliseconds& n) const noexcept
    { return milliseconds{m_msec - n.m_msec}; }
    
    /// Division operator between milliseconds.
    constexpr milliseconds operator/(const milliseconds& sc) noexcept
    { return milliseconds{m_msec/sc.m_msec}; }
    
    /// Equality operator.
    constexpr bool operator==(const milliseconds& d) const noexcept
    { return m_msec == d.m_msec; }

    /// Greater than operator.
    constexpr bool operator>(const milliseconds& d) const noexcept
    { return m_msec > d.m_msec; }

    /// Greater or equal to operator.
    constexpr bool operator>=(const milliseconds& d) const noexcept
    { return m_msec >= d.m_msec; }

    /// Less than operator.
    constexpr bool operator<(const milliseconds& d) const noexcept
    { return m_msec < d.m_msec; }

    /// Less or equal to operator.
    constexpr bool operator<=(const milliseconds& d) const noexcept
    { return m_msec <= d.m_msec; }
    
    /// Do the milliseconds sum up to more than one day ?
    constexpr bool more_than_day() const noexcept { return m_msec>max_in_day; }
    
    /// Get the milliseconds cast to the underlying type.
    constexpr underlying_type as_underlying_type() const noexcept
    { return m_msec; }
    
    /// If the milliseconds sum up to more (or equal to) one day, remove the 
    /// integral days (and return them); reset the milliseconds to milliseconds
    /// of the new day.
    constexpr int remove_days() noexcept
    {
        int day { static_cast<int>(m_msec/max_in_day) };
        m_msec %= max_in_day;
        return day;
    }
    
    /// Return the milliseconds as integral day(s).
    constexpr int to_days() const noexcept
    {
        return int{static_cast<int>(m_msec/max_in_day)};
    }
    
    /// Cast to fractional days.
    constexpr double fractional_days() const noexcept
    {
        return static_cast<double>(m_msec)/static_cast<double>(max_in_day);
    }
    
    /// Cast to fractional seconds
    constexpr double to_fractional_seconds() const noexcept
    { return static_cast<double>(m_msec)*1.0e-3; }
    
    /// Resolve to (integer) seconds and fractional seconds.
    constexpr seconds resolve_sec(double& fraction) const noexcept
    {
        seconds sec { m_msec/1000L };
        fraction = static_cast<double>(m_msec%1000L)*1e-3;
        return sec;
    }
    
    /// Cast to any arithmetic type.
    template<typename T,
             typename=std::enable_if_t<std::is_arithmetic<T>::value>
             >
    constexpr T cast_to() const noexcept
    { return static_cast<T>( m_msec ); }
    
    /// Translate to hours, minutes, seconds and fraction of seconds
    constexpr std::tuple<hours, minutes, seconds, double>
    to_hmsf() const noexcept
    {
        long hr { m_msec/3600000L                    };  // hours
        long mn { (m_msec%3600000L)/60000L           };  // minutes
        long sc { ((m_msec%3600000L)%60000L)/1000L   };  // seconds
        long ms { m_msec-(hr*3600L+mn*60L+ sc)*1000L };  // milliseconds
        return std::make_tuple( hours  { static_cast<hours::underlying_type>(hr) },
                                minutes{ static_cast<minutes::underlying_type>(mn) },
                                seconds{ sc },
                                static_cast<double>(ms) );
    }

private:
    /// Milliseconds as underlying type.
    underlying_type m_msec;

}; /// class milliseconds

/// A wrapper class for microseconds (i.e 10**-6 sec.).
class microseconds {

public:
    /// Nanoseconds are represented as long integers.
    typedef long underlying_type;
    
    /// Nanoseconds is a subdivision of seconds.
    static constexpr bool is_of_sec_type { true };
    
    /// Max microseconds in day.
    static constexpr long max_in_day { 86400L * 1000000L };

    template<typename T>
    static constexpr T sec_factor() noexcept
    { return static_cast<T>(1000000); }
    template<typename T>
    static constexpr T sec_ifactor() noexcept
    { return ((static_cast<T>(1))/1000000); }
    
    /// Constructor.
    explicit constexpr microseconds(underlying_type i=0L) noexcept : s(i) {};
    
    explicit constexpr microseconds(hours h, minutes m, microseconds c) noexcept
        : s { c.as_underlying_type()
            + m.as_underlying_type()*60L   *1000L * 1000L
            + h.as_underlying_type()*3600L *1000L * 1000L}
    {}
    
    /// Constructor from hours, minutes, fractional seconds
    /*explicit constexpr microseconds(hours h, minutes m, double fs)
    noexcept
        : s(  static_cast<long>(fs*1e6)
            + (m.as_underlying_type()*60L
            + h.as_underlying_type()*3600L) * 1000000L
            )
    {}*/
    
    /// Nanoseconds can be cast to milliseconds will a loss of accuracy.
    constexpr explicit operator milliseconds() const
    { return milliseconds(s/1000L); }
    
    /// Nanoseconds can be cast to seconds will a loss of accuracy.
    constexpr explicit operator seconds() const { return seconds(s/1000000L); }
    
    /// Overload operatpr "+=" between microseconds.
    constexpr void operator+=(const microseconds& ns) noexcept { s+=ns.s; }
    constexpr void operator-=(const microseconds& ns) noexcept { s-=ns.s; }

    constexpr microseconds operator+(const microseconds& sec) const noexcept
    { return microseconds(s+sec.s); }

    /// Overload '-' operator
    constexpr microseconds operator-(const microseconds& n) const noexcept
    { return microseconds(s - n.s); }
    
    /// Overload operatpr "/" between microseconds.
    constexpr microseconds operator/(const microseconds& sc) noexcept
    { return microseconds(s/sc.s); }
    
    constexpr bool operator==(const microseconds& d) const noexcept
    { return s == d.s; }

    constexpr bool operator>(const microseconds& d) const noexcept
    { return s > d.s; }

    constexpr bool operator>=(const microseconds& d) const noexcept
    { return s >= d.s; }

    constexpr bool operator<(const microseconds& d) const noexcept
    { return s < d.s; }

    constexpr bool operator<=(const microseconds& d) const noexcept
    { return s <= d.s; }
    
    /// Do the microseconds sum up to more than one day?
    constexpr bool more_than_day() const noexcept { return s>max_in_day; }
    
    /// Cast to underlying type.
    constexpr underlying_type as_underlying_type() const noexcept { return s; }
    
    /// Access (set/get) the underlying type.
    constexpr underlying_type& assign() noexcept { return s; }
    
    /// If the microseconds sum up to more (or equal to) one day, remove the
    /// integral days (and return them); reset the microseconds to microseconds
    /// of the new day.
    constexpr day remove_days() noexcept {
        day d ( static_cast<int>(s/max_in_day) );
        s %= max_in_day;
        return d;
    }
    
    /// Cast to days.
    constexpr day to_days() const noexcept {
        return day(static_cast<day::underlying_type>(s/max_in_day));
    }
    
    /// Cast to fractional days.
    constexpr double fractional_days() const noexcept {
        return static_cast<double>(s)/static_cast<double>(max_in_day);
    }
    
    /// Cast to fractional seconds
    constexpr double to_fractional_seconds() const noexcept
    { return static_cast<double>(s)*1.0e-6; }
    
    
    /// Translate to hours, minutes, seconds and fraction of seconds.
    constexpr std::tuple<hours, minutes, seconds, double>
    to_hmsf() const noexcept
    {
        long hr { m_msec/3600000000L                       };  // hours
        long mn { (m_msec%3600000000L)/60000000L           };  // minutes
        long sc { ((m_msec%3600000000L)%60000000L)/1000000L};  // seconds
        long ns { m_msec-(hr*3600L+mn*60L+sc)*1000000L     };  // nanoseconds
        return std::make_tuple( hours  { static_cast<hours::underlying_type>(hr) },
                                minutes{ static_cast<minutes::underlying_type>(mn) },
                                seconds{ sc },
                                static_cast<double>(ns) );
    }

private:
    /// Cast to any arithmetic type.
    template<typename T,
             typename=std::enable_if_t<std::is_arithmetic<T>::value>
             >
    constexpr T cast_to() const noexcept
    { return static_cast<T>(m_msec); }

    /// Microseconds as long ints.
    underlying_type m_msec;

}; // class microseconds

/// Calendar date (i.e. year, momth, day) to MJDay.
modified_julian_day cal2mjd(year, month, day_of_month);

/// Valid output formats
enum class datetime_output_format : char
{
    ymd, ymdhms, ydhms, gps, ydoy, jd, mjd
};

/*
 * A datetime class. Holds (integral) days as MJD and fraction of day as any
 * of the is_of_sec_type class (i.e. seconds/milli/nano).
 */
template<class S>
class datev2 {
public:

    /// Only allow S parameter to be of sec type (seconds/milli/nano).
    static_assert( S::is_of_sec_type, "" );
    
    /// Default (zero) constructor.
    explicit constexpr datev2() noexcept : mjd_(0), sect_(0) {};
    
    /// Constructor from year, month, day of month and any sec type.
    explicit constexpr datev2(year y, month m, day_of_month d, S s)
        : mjd_ (cal2mjd(y, m, d)),
          sect_(s)
        {}

    template<class T>
    explicit datev2(year y, month m, day_of_month d, T t)
        : mjd_ (cal2mjd(y, m, d)),
          sect_(S(t))
    {}

    explicit
    datev2(year y, month m, day_of_month d, hours hr, minutes mn, double fsecs)
        : mjd_ ( cal2mjd(y, m, d) ),
          sect_( hr, mn, fsecs )
    {}

    explicit
    datev2(modified_julian_day mjd, hours hr=hours(), minutes mn=minutes(), S sec=S())
        : mjd_ {mjd}, 
          sect_{hr, mn, sec}
    {}
    
    explicit
    datev2(year y, month m, day_of_month d,
           hours hr=hours(), minutes mn=minutes(), S sec=S())
        : mjd_ {cal2mjd(y, m, d)}, 
          sect_{hr, mn, sec}
    {}

    constexpr modified_julian_day mjd() const noexcept { return mjd_; }

    template<class T>
    constexpr void add_seconds(T t) noexcept
    { 
        sect_ += (S)t;
        if ( sect_.more_than_day() ) this->normalize();
        return;
    }
    
    template<class T>
    constexpr void remove_seconds(T t) noexcept
    { 
        sect_ -= (S)t;
        if ( sect_ < (S)0 ) {
            while ( sect_ < (S)0 ) {
                --mjd_;
                sect_ += (S)S::max_in_day;
            }
        }
        return;
    }

    constexpr S delta_sec(const datev2& d) const noexcept
    { return (sect_ - d.sect_) + to_sec_type<S>(mjd_ - d.mjd_); }

    /// Overload equality operator.
    constexpr bool operator==(const datev2& d) const noexcept
    { return mjd_ == d.mjd_ && sect_ == d.sect_; }

    /// Overload ">" operator.
    constexpr bool operator>(const datev2& d) const noexcept
    { return mjd_ > d.mjd_ || (mjd_ == d.mjd_ && sect_ > d.sect_); }
    
    /// Overload ">=" operator.
    constexpr bool operator>=(const datev2& d) const noexcept
    { return mjd_ > d.mjd_ || (mjd_ == d.mjd_ && sect_ >= d.sect_); }
    
    /// Overload "<" operator.
    constexpr bool operator<(const datev2& d) const noexcept
    { return mjd_ < d.mjd_ || (mjd_ == d.mjd_ && sect_ < d.sect_); }
    
    /// Overload "<=" operator.
    constexpr bool operator<=(const datev2& d) const noexcept
    { return mjd_ < d.mjd_ || (mjd_ == d.mjd_ && sect_ <= d.sect_); }

    /// Reset the seconds/milli/nano after removing whole days.
    constexpr void normalize() noexcept
    {
        mjd_ += sect_.remove_days();
        return;
    }

    /// Cast to double Modified Julian Date.
    constexpr double as_mjd() const noexcept
    {
        return static_cast<double>(mjd_.as_underlying_type())
                                + sect_.fractional_days();
    }

    /// Cast to gps_datetime.
    constexpr gps_datetime as_gps_datetime() const noexcept
    {
        long week   = (mjd_.as_underlying_type() - jan61980)/7L;
        double secs = sect_.as_fractional_seconds();
        secs       += static_cast<double>(
                    ((mjd_.as_underlying_type() - jan61980) - week*7L ) * 86400L);
        return gps_datetime( week, secs );
    }

    /// Cast to year, month, day of month
    /*constexpr*/ std::tuple<year, month, day_of_month>
    as_ymd() const noexcept
    {
        year y;
        month m;
        day_of_month d;
        y = mjd_.to_ymd(m, d);
        return std::make_tuple( y, m, d );
    }

    /// Cast to year, day_of_year
    /*constexpr*/ std::tuple<year, day_of_year> as_ydoy() const noexcept
    {
        year y;
        day_of_year d;
        y = mjd_.to_ydoy( d );
        return std::make_tuple( y, d );
    }

    /// Convert the *seconds to hours, minutes, seconds and fractional seconds
    /*constexpr*/ std::tuple<hours, minutes, seconds, long>
    as_hms() const noexcept
    { return sect_.to_hms(); }

    std::string stringify() const
    {
        auto ymd { this->as_ymd() };
        auto hms { this->as_hms() };
        /*double fsec = S::template sec_ifactor<double>() * std::get<3>(hms);*/
        return std::string {
                     std::to_string( std::get<0>(ymd).as_underlying_type() )
             + "/" + std::to_string( std::get<1>(ymd).as_underlying_type() )
             + "/" + std::to_string( std::get<2>(ymd).as_underlying_type() )
             + " " + std::to_string( std::get<0>(hms).as_underlying_type() )
             + ":" + std::to_string( std::get<1>(hms).as_underlying_type() )
             + ":" + std::to_string( std::get<2>(hms).as_underlying_type() )
             + "." + std::to_string( S::template sec_ifactor<int>() * std::get<3>(hms) )
            };
    }
    
    /// Overload operator "<<" TODO
    /*
    friend std::ostream& operator<<(std::ostream& o, const datev2& d)
    {
        o << d.stringify();
        return o;
    }
    */

private:
    modified_julian_day mjd_;  ///< Modified Julian Day
    S                   sect_; ///< Fraction of day in milli/nano/seconds
};


/*
/// Cast a modified_julian_day to year, day_of_year
constexpr ngpt::year
ngpt::modified_julian_day::to_ydoy(day_of_year& d)
const noexcept
{
    long days_fr_jan1_1901 = m - ngpt::jan11901;
    long num_four_yrs      = days_fr_jan1_1901/1461L;
    long years_so_far      = 1901L + 4*num_four_yrs;
    long days_left         = days_fr_jan1_1901 - 1461*num_four_yrs;
    long delta_yrs         = days_left/365 - days_left/1460;
    day_of_year yday (static_cast<day_of_year::underlying_type>
                                (days_left - 365*delta_yrs + 1));
    d = yday;
    return year (years_so_far + delta_yrs);
}
*/
/*
constexpr ngpt::year
ngpt::modified_julian_day::to_ymd(ngpt::month& mm, ngpt::day_of_month& dd)
const noexcept
{
    ngpt::day_of_year doy;
    auto y     = this->to_ydoy( doy );
    long yday  = static_cast<long>( doy.as_underlying_type() );
    long leap  = ( y.as_underlying_type()%4L == 0 );
    long guess = yday*0.032;
    long more  = (( yday - month_day[leap][guess+1] ) > 0);
    mm.assign() = static_cast<ngpt::month::underlying_type>
                                    (guess + more + 1);
    dd.assign() = static_cast<ngpt::day_of_month::underlying_type>
                                    (yday - month_day[leap][guess+more]);
    return y;
}
*/

} // end namespace

#endif // define DATETIME
