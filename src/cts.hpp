#ifndef __CTS_TS__
#define __CTS_TS__

#include "datetime.hpp"


enum class bp_flag : uint_8 {
  outlier  = 0x01,
  skip     = 0x02,
  erthq    = 0x04,
  velchg   = 0x08,
  ofst     = 0x10
};

using T = std::underlying_type_t<bp_flag>;

inline bp_flag operator|(bp_flag lhs, bp_flag rhs)

{
    return (bp_flag)(static_cast<T>(lhs) | static_cast<T>(rhs));
}

inline bp_flag& operator|=(bp_flag& lhs, bp_flag rhs)
{
    lhs = (bp_flag)(static_cast<T>(lhs) | static_cast<T>(rhs));
    return lhs;
}


class flag {
public:
  flag() noexcept : bits_{0x00} {};
  
  uint_8 
  is_outlier() 
  noexcept const
  {
      return bits_ & static_cast<uint_8>( bp_flag::outlier );
  }

  uint_8&
  set_outlier()
  noexcept const
  {
      return bits_ |= 1 << static_cast<uint_8>( bp_flag::outlier );
  }

  uint_8&
  unset_outlier()
  noexcept const
  {
      return bits_ |= 1 << static_cast<uint_8>( bp_flag::outlier );
  }


private:
  uint_8 bits_;
};



/*
class ts {

using ngpt::datetime;

public:
    ts( std::vector<datetime>* e = nullptr, std::size_t hint = 0 )
        : epochs_{ e }
    {
      if ( hint ) {
        vals_.reserve( hint );
        sigmas_.reserve( hint );
      }
    }

    ~ts() noexcept {
      epochs_ = nullptr;
    }

    friend void swap(ts& lhs, ts& rhs) noexcept
    {
        std::swap(lhs.epochs_, rhs.epochs_);
        std::swap(lhs.vals_,   rhs.vals_);
        std::swap(lhs.sigmas_, rhs.sigmas_);
    }

    ts(const ts& rhs)
      noexcept(std::is_nothrow_copy_constructible<std::vector<double>>::value)
      : epochs_{ rhs.epochs_},
        vals_  { rhs.vals_  },
        sigmas_{ rhs.sigmas_}
    {}
    
    ts(ts&& rhs)
      noexcept(std::is_nothrow_move_constructible<std::vector<double>>::value)
      : epochs_{ std::move(rhs.epochs_) },
        vals_  { std::move(rhs.vals_)   },
        sigmas_{ std::move(rhs.sigmas_) }
    {
        rhs.epochs_ = nullptr;
    }

    ts& operator=(ts rhs) noexcept
    {
        swap(*this, rhs);
        return this;
    }

    std::vector<datetime>* epoch_ptr() noexcept
    {
        return epochs_;
    }

    void push_back(double val, double sigma)
    {
        vals_.push_back( val );
        sigmas_.push_back( sigma );
    }

    double& value(std::size_t i)
    {
        vals_[i];
    }
    
    double& sigma(std::size_t i)
    {
        sigmas_[i];
    }
    
private:
    std::vector<datetime>* epochs_;
    std::vector<double>    vals_;
    std::vector<double>    sigmas_;
};
*/

#endif
