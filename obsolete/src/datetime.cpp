#include "datetime.hpp"

using ngpt::datetime;

datetime::datetime(int year, int month, int day, 
    int hour, int minute, double sec) 
noexcept
{
  long doy;

  if ( year % 4 ) {
      doy = NormalMonths[month - 1] + day;
  } else {
      doy = LeapMonths[month - 1] + day;
  }

  mjd_ = ((year - 1901)/4)*1461 + 
         ((year - 1901)%4)*365 + 
         doy - 1 + ngpt::_JAN11901_;

  fday_ = ((sec/60.0 + minute)/60.0 + hour)/24.0;
}
