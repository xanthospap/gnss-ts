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
{ return !(iy%4) && (iy%100 || !(iy%400)); }

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

/// Format options
namespace datetime_format_options {
    enum class year_digits  : char { two_digit, four_digit };
    enum class month_format : char { two_digit, short_name, long_name };
}

/// Time systems
enum class time_system : char {
    gps,
    glonass,
    galileo,
    tai,
    utc,
    ut,
    qzss
};

/// Forward declare cout operator for year ...
/*
template<datetime_format_options::year_digits F = datetime_format_options::year_digits::four_digit>
std::ostream& operator<<(std::ostream&, const year&);
*/

/// ... and month
/*
template<datetime_format_options::month_format F = datetime_format_options::month_format::two_digit>
std::ostream& operator<<(std::ostream&, const month&);
*/

/// A wrapper class for years.
class year {
    int y;
public:
    /// Years are represented as integers.
    typedef int underlying_type;

    /// Constructor.
    explicit constexpr year (underlying_type i=0) noexcept : y(i) {};

    /// Get the underlying int.
    constexpr underlying_type as_underlying_type() const noexcept
    { return y; }

    /// overload operator "<<"
    /*
    template<datetime_format_options::year_digits F>
    friend std::ostream& operator<<(std::ostream&, const year&);
    */
};

/// Declare the '<<' operator(s) (definition in .cpp)
/*
template<datetime_format_options::year_digits F>
std::ostream& operator<<(std::ostream&, const year&);
*/

/// Specialization of 2-digit year
/*
template<>
std::ostream& operator<<<datetime_format_options::year_digits::two_digit>
(std::ostream&, const year&);
*/

/// A wrapper class for months.
class month {
    int m;
    
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
public:
    /// Months are represented as int.
    typedef int underlying_type;

    /// Constructor.
    explicit constexpr month (underlying_type i=1) noexcept : m(i) {};

    /// Get the underlying int.
    constexpr underlying_type as_underlying_type() const noexcept
    { return m; }

    /// Access (get/set) the underlying type
    constexpr underlying_type& assign() noexcept { return m; }

    const char* short_name() const noexcept { return short_names[m-1]; }
    const char* long_name()  const noexcept { return long_names[m-1];  }
    
    /// overload operator "<<"
    /*
    template<datetime_format_options::month_format F>
    friend std::ostream& operator<<(std::ostream&, const month&);
    */
};

/*
template<datetime_format_options::month_format F>
std::ostream& operator<<(std::ostream& o, const month& mm);

template<>
std::ostream& operator<<<datetime_format_options::month_format::short_name>
(std::ostream& o, const month& mm);

template<>
std::ostream& operator<<<datetime_format_options::month_format::long_name>
(std::ostream& o, const month& mm);
*/

/// A wrapper class for days (in general!).
class day {
    int d;
public:
    /// Days are represented as int.
    typedef int underlying_type;
    
    /// Constructor.
    explicit constexpr day(underlying_type i=0) noexcept : d(i) {};
    
    /// Get the underlying int.
    constexpr underlying_type as_underlying_type() const noexcept
    { return d; }
};

/// A wrapper class for day of month.
class day_of_month {
    int d;
public:
    /// Days are represented as int.
    typedef int underlying_type;
    
    /// Constructor.
    explicit constexpr day_of_month(underlying_type i=0) noexcept : d(i) {};
    
    /// Get the underlying int.
    constexpr underlying_type as_underlying_type() const noexcept
    { return d; }
    
    /// Access (get/set) the underlying type
    constexpr underlying_type& assign() noexcept { return d; }
    
    /// Overload << operator
    friend std::ostream& operator<<(std::ostream& o, const day_of_month& dm)
    { o << std::setw(2) << dm.d; return o; }
};

/// A wrapper class for Modified Julian Days.
class modified_julian_day {
    long m;
public:
    /// MJDs are represented as long ints.
    typedef long underlying_type;
    
    /// Constructor.
    explicit constexpr modified_julian_day(underlying_type i=0) noexcept
        : m(i) 
    {};
    
    /// Get the underlying long int.
    constexpr underlying_type as_underlying_type() const noexcept
    { return m; }
    
    /// Access (get/set) the underlying long.
    constexpr underlying_type& assign() noexcept { return m; }
    
    /// Define addition (between MJDs).
    constexpr void operator+=(const modified_julian_day& d) noexcept
    { m += d.m; }
    
    /// Define addition (between an MJDs and a day).
    constexpr void operator+=(const day& d) noexcept
    { m += static_cast<underlying_type>(d.as_underlying_type()); }

    /// Overload - operator
    constexpr modified_julian_day operator-(const modified_julian_day& mjd) const noexcept
    { return modified_julian_day(m-mjd.m); }

    constexpr bool operator==(const modified_julian_day& d) const noexcept
    { return m == d.m; }

    constexpr bool operator>(const modified_julian_day& d) const noexcept
    { return m > d.m; }

    constexpr bool operator>=(const modified_julian_day& d) const noexcept
    { return m >= d.m; }

    constexpr bool operator<(const modified_julian_day& d) const noexcept
    { return m < d.m; }

    constexpr bool operator<=(const modified_julian_day& d) const noexcept
    { return m <= d.m; }

    constexpr modified_julian_day& operator--() noexcept
    { --m; return *this; }

    /// Cast to year, day_of_year
    /*constexpr*/ year to_ydoy(day_of_year&) const noexcept;
    
    /// Cast to year, month, day_of_month
    /*constexpr*/ year to_ymd(month&, day_of_month&) const noexcept;
    
    /// overload operator <<
    friend std::ostream&
    operator<<(std::ostream& o, const modified_julian_day& mjd)
    { o << std::setw(5) << mjd.m; return o; }
};
    
/// Trasnform to sec_type
template<typename S>
S to_sec_type(const modified_julian_day& mjd) noexcept
{ return S(mjd.as_underlying_type() * S::max_in_day); }

/// A wrapper class for Julian Days.
/// TODO JD have a fraction part!
class julian_day {
    long j;
public:
    /// MJDs are represented as long ints.
    typedef long underlying_type;
    
    explicit constexpr julian_day(long i=0) noexcept : j(i) {};
    
    constexpr long as_long() const noexcept { return j; }
    /*
    constexpr long& assign() noexcept { return j; }
    constexpr void operator+=(const julian_day& d) noexcept
    { j += d.j; }
    constexpr void operator+=(const day& d) noexcept
    { j += d.as_int(); }
    */
};

/// A wrapper class for day of year
class day_of_year {
    int d;
public:
    /// DOY represented by ints.
    typedef int underlying_type;
    
    /// Constructor.
    explicit constexpr day_of_year(underlying_type i=0) noexcept : d(i) {};
    
    /// Cast to underlying type
    constexpr underlying_type as_underlying_type() const noexcept
    { return d; }
    
    /// overload << operator
    friend std::ostream& operator<<(std::ostream& o, const day_of_year&  doy)
    { o << std::setw(3) << doy.d; return o; }
};

class hours {
    int h;
public:
    /// Hours are represented by ints.
    typedef int underlying_type;
    
    /// Constructor
    explicit constexpr hours(underlying_type i=0) noexcept : h(i) {};
    
    /// Pass the underlying type
    constexpr underlying_type as_underlying_type() const noexcept
    { return h; }

    friend std::ostream& operator<<(std::ostream& o, const hours& hr)
    {
        o << std::setw(2) << hr.h;
        return o;
    }
};

class minutes {
    int m;
public:
    /// Minutes are represented by ints
    typedef int underlying_type;
    
    /// Constructor
    explicit constexpr minutes(underlying_type i=0) noexcept : m(i) {};

    /// Pass the underlying type
    constexpr underlying_type as_underlying_type() const noexcept
    { return m; }

    friend std::ostream& operator<<(std::ostream& o, const minutes& mn)
    {
        o << std::setw(2) << mn.m;
        return o;
    }
};

/// A wrapper class for seconds.
class seconds {
    long s;
public:
    /// Seconds are represented as long ints.
    typedef long underlying_type;
    
    /// Seconds is a subdivision of seconds.
    static constexpr bool is_of_sec_type { true };
    
    /// Max seconds in day.
    static constexpr underlying_type max_in_day { 86400L };
    
    template<typename T>
    static constexpr T sec_factor() noexcept
    { return static_cast<T>(1); }
    template<typename T>
    static constexpr T sec_ifactor() noexcept
    { return static_cast<T>(1); }
    
    /// Constructor
    explicit constexpr seconds(underlying_type i=0L) noexcept : s(i) {};

    /// Constructor from hours, minutes, seconds
    explicit constexpr seconds(hours h, minutes m, seconds c) noexcept
        : s {  c.as_underlying_type()
             + m.as_underlying_type()*60L
             + h.as_underlying_type()*3600L}
    {}

    /// Constructor from hours, minutes, fractional seconds
    explicit constexpr seconds(hours h, minutes m, double fs) noexcept
        : s(  static_cast<long>(fs)
            + m.as_underlying_type()*60L
            + h.as_underlying_type()*3600L
            )
    {}
    
    /// Addition operator between seconds.
    constexpr void operator+=(const seconds& sc) noexcept { s+=sc.s; }
    constexpr void operator-=(const seconds& sc) noexcept { s-=sc.s; }
    
    /// Overload '-' operator
    constexpr seconds operator-(const seconds& n) const noexcept
    { return seconds(s - n.s); }

    constexpr seconds operator+(const seconds& sec) const noexcept
    { return seconds(s+sec.s); }

    /// Overload operator '/'
    constexpr seconds operator/(const seconds& n) const noexcept
    { return seconds(s/n.s); }
    
    constexpr bool operator==(const seconds& d) const noexcept
    { return s == d.s; }

    constexpr bool operator>(const seconds& d) const noexcept
    { return s > d.s; }

    constexpr bool operator>=(const seconds& d) const noexcept
    { return s >= d.s; }

    constexpr bool operator<(const seconds& d) const noexcept
    { return s < d.s; }

    constexpr bool operator<=(const seconds& d) const noexcept
    { return s <= d.s; }
    
    /// Do the secods sum up to more than one day?
    constexpr bool more_than_day() const noexcept { return s>max_in_day; }
    
    /// Get the underlying type numeric.
    constexpr underlying_type as_underlying_type() const noexcept { return s; }
    
    /// Access (get/set) the underlying type (long).
    constexpr underlying_type& assign() noexcept { return s; }
    
    /// If the seconds sum up to more (or equal to) one day, remove the integer
    /// days (and return them); reset the seconds to seconds of the new day.
    constexpr day remove_days() noexcept {
        day d ( static_cast<day::underlying_type>(s/max_in_day) );
        s %= max_in_day;
        return d;
    }
    
    /// Return the integral number of days.
    constexpr day to_days() const noexcept {
        return day(static_cast<day::underlying_type>(s/max_in_day));
    }
    
    /// Interpret the seconds as fractional days.
    constexpr double fractional_days() const noexcept {
        return static_cast<double>(s)/static_cast<double>(max_in_day);
    }
    
    /// Cast to double (i.e. fractional seconds)
    constexpr double to_fractional_seconds() const noexcept
    { return static_cast<double>(s); }

    /// Resolve to (integer) seconds and fractional seconds.
    constexpr seconds resolve_sec(double& fraction) const noexcept
    {
        fraction = 0.0e0;
        return *this;
    }
    
    /// Translate to hours, minutes, seconds and fractional seconds
    constexpr std::tuple<hours, minutes, seconds, long> to_hms() const noexcept
    {
        return std::make_tuple(hours  {static_cast<int>(s / 3600L)},
                               minutes{static_cast<int>((s % 3600L) / 60L)},
                               seconds{(s % 3600) % 60},
                               0L );
    }

    /// Cast to any arithmetic type.
    template<typename T,
             typename=std::enable_if_t<std::is_arithmetic<T>::value>
             >
    constexpr T cast_to() const noexcept
    { return static_cast<T>( s ); }

    friend std::ostream& operator<<(std::ostream& o, const seconds& sec)
    {
        o << sec.s << "sec";
        return o;
    }
};

/// A wrapper class to represent a datetime in GPSTime, i.e. gps week and
/// seconds of (gps) week.
class gps_datetime {
    long   week_;
    double sec_of_week_;
public:
    explicit constexpr gps_datetime(long w, double s) noexcept
        : week_(w), sec_of_week_(s) {};
};

/// A wrapper class for milliseconds.
class milliseconds {
    long s;
public:
    /// MilliSeconds are represented as long ints.
    typedef long underlying_type;
    
    /// MilliSeconds are a subdivision of seconds.
    static constexpr bool is_of_sec_type { true };
    
    /// Max milliseconds in one day.
    static constexpr long max_in_day { 86400L * 1000L };
    
    template<typename T>
    static constexpr T sec_factor() noexcept
    { return static_cast<T>(1000); }
    template<typename T>
    static constexpr T sec_ifactor() noexcept
    { return ((static_cast<T>(1))/1000); }
    
    /// Cinstructor.
    explicit constexpr milliseconds(underlying_type i=0L) noexcept : s(i) {};
    
    explicit constexpr milliseconds(hours h, minutes m, milliseconds c) noexcept
        : s { c.as_underlying_type()
            + m.as_underlying_type()*60L  *1000L
            + h.as_underlying_type()*3600L*1000L}
    {}
    
    /// Constructor from hours, minutes, fractional seconds
    explicit constexpr milliseconds(hours h, minutes m, double fs)
    noexcept
        : s(  static_cast<long>(fs*1000.0e0)
            + (m.as_underlying_type()*60L
            + h.as_underlying_type()*3600L) * 1000L
            )
    {}
    
    /// Milliseconds can be cast to seconds (with a loss of precission).
    constexpr explicit operator seconds() const { return seconds(s/1000L); }
    
    constexpr milliseconds operator+(const milliseconds& sec) const noexcept
    { return milliseconds(s+sec.s); }
    
    /// Overload operator "+=" between milliseconds.
    constexpr void operator+=(const milliseconds& ms) noexcept { s+=ms.s; }
    constexpr void operator-=(const milliseconds& ms) noexcept { s-=ms.s; }
    
    /// Overload '-' operator
    constexpr milliseconds operator-(const milliseconds& n) const noexcept
    { return milliseconds(s - n.s); }
    
    /// Overload operator "/" between milliseconds.
    constexpr milliseconds operator/(const milliseconds& sc) noexcept
    { return milliseconds(s/sc.s); }
    
    constexpr bool operator==(const milliseconds& d) const noexcept
    { return s == d.s; }

    constexpr bool operator>(const milliseconds& d) const noexcept
    { return s > d.s; }

    constexpr bool operator>=(const milliseconds& d) const noexcept
    { return s >= d.s; }

    constexpr bool operator<(const milliseconds& d) const noexcept
    { return s < d.s; }

    constexpr bool operator<=(const milliseconds& d) const noexcept
    { return s <= d.s; }
    
    /// Do the milliseconds sum up to more than one day ?
    constexpr bool more_than_day() const noexcept { return s>max_in_day; }
    
    /// Get the milliseconds cast to the underlying type.
    constexpr underlying_type as_underlying_type() const noexcept { return s; }
    
    /// Access (get/set) the underlying type (long int).
    constexpr underlying_type& assign() noexcept { return s; }
    
    /// If the milliseconds sum up to more (or equal to) one day, remove the 
    /// integral days (and return them); reset the milliseconds to milliseconds
    /// of the new day.
    constexpr day remove_days() noexcept {
        day d ( static_cast<day::underlying_type>(s/max_in_day) );
        s %= max_in_day;
        return d;
    }
    
    /// Return the milliseconds as whole day(s) .
    constexpr day to_days() const noexcept {
        return day(static_cast<day::underlying_type>(s/max_in_day));
    }
    
    /// Cast to fractional days.
    constexpr double fractional_days() const noexcept {
        return static_cast<double>(s)/static_cast<double>(max_in_day);
    }
    
    /// Cast to fractional seconds
    constexpr double to_fractional_seconds() const noexcept
    { return static_cast<double>(s)*1.0e-3; }
    
    /// Resolve to (integer) seconds and fractional seconds.
    constexpr seconds resolve_sec(double& fraction) const noexcept
    {
        seconds sec ( s / 1000L );
        fraction = static_cast<double>( s % 1000L ) * 1e-3;
        return sec;
    }
    
    /// Cast to any arithmetic type.
    template<typename T,
             typename=std::enable_if_t<std::is_arithmetic<T>::value>
             >
    constexpr T cast_to() const noexcept
    { return static_cast<T>( s ); }
    
    /// Translate to hours minutes and milliseconds
    constexpr std::tuple<hours, minutes, seconds, long> to_hms() const noexcept
    {
        long hr { s / 3600000L                           };  // hours
        long mn { (s % 3600000L) / 60000L                };  // minutes
        long sc { ((s % 3600000L) % 60000L) / 1000L      };  // seconds
        long ms { s - ( hr*3600L + mn*60L + sc ) * 1000L };  // milliseconds
        return std::make_tuple( hours  { (int)hr },
                                minutes{ (int)mn },
                                seconds{ sc },
                                ms
                              );
    }

    friend std::ostream& operator<<(std::ostream& o, const milliseconds& ms)
    {
        o << ms.s << "millisec";
        return o;
    }
    
};

/// A wrapper class for microseconds.
class microseconds {
    long s;
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
    
    /// Cast to any arithmetic type.
    template<typename T,
             typename=std::enable_if_t<std::is_arithmetic<T>::value>
             >
    constexpr T cast_to() const noexcept
    { return static_cast<T>( s ); }
    
    /// Translate to hours minutes and nanoseconds
    constexpr std::tuple<hours, minutes, seconds, long> to_hms() const noexcept
    {
        long hr { s / 3600000000L                           };  // hours
        long mn { (s % 3600000000L) / 60000000L             };  // minutes
        long sc { ((s % 3600000000L) % 60000000L) / 1000000L};  // seconds
        long ns { s - ( hr*3600L + mn*60L + sc ) * 1000000L };  // nanoseconds
        return std::make_tuple( hours  { (int)hr },
                                minutes{ (int)mn },
                                seconds{ sc },
                                ns
                               );
    }
};

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
