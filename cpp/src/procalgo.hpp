#ifndef __NGPT_CRD_GENALGO__
#define __NGPT_CRD_GENALGO__

// standard headers
#include <stdexcept>
#ifdef DEBUG
#include <iostream>
#endif
// gtms headers
#include "timeseries.hpp"
#include "earthquake_cat.hpp"
#include "period.hpp"

namespace ngpt
{

template<class T>
  ngpt::ts_model<T>
identify_harmonics(ngpt::timeseries<T, ngpt::pt_marker>& ts,
    ngpt::ts_model<T>& model, double Ut=1e-3, double min_amplitude=0e0)
{
  // Time-span of the time-series in days and years.
  auto   tdif = ts.last_epoch().delta_date(ts.first_epoch());
  double ddif = tdif.days().as_underlying_type() / 365.25e0;
  double div  = 1e0 / ddif;
  std::size_t N = ts.data_pts() - ts.skipped_pts();
  double ofac{4},
         hifac{div/div},
         *px,
         *py,
         prob,
         *mempool;
  int    nout (0.5*ofac*hifac*N+1),
         jmax;
  try {
    mempool = new double[2*nout];
  } catch (std::bad_alloc&) {
    std::cerr<<"\n[ERROR] Cannot allocate memmory for identify_harmonics";
    return model;
  }
  px = mempool;
  py = mempool + nout;
  auto   nmodel {model},
         amodel {model};
  double nstddev,
         astddev,
         factor;
  auto nts{ts};
  // first fit.
  nts = nts.qr_ls_solve(amodel, astddev, 1e-3, false, true);
  std::size_t dummy_it = 0;
  while ( dummy_it < 10 ) {
    ngpt::lomb_scargle_period(nts, ofac, hifac, px, py, nout, nout, jmax, prob);
    std::cout<<"\n\tDominant frequency in time-series: "<<px[jmax]<<" (at: "<<jmax<<")"
      <<"; this is a period of "<<1e0/px[jmax]<<" days -- posibility: "<<prob<<" --";
    nmodel.add_period( 1e0/px[jmax] );
    auto res_ts = nts.qr_ls_solve(nmodel, nstddev, 1e-3, false, false);
    factor = astddev / nstddev;
    double ampl = nmodel.harmonic_with_period(1e0/px[jmax])->amplitude();
    if ((factor-1e0) > Ut && ampl > min_amplitude) {
      nts = res_ts;
      amodel.add_period( 1e0/px[jmax] );
      astddev = nstddev;
      nmodel = amodel;
      std::cout<<"\n\tFrequency added to the model; factor was: "<<(factor-1e0);
    } else {
      break;
    }
    ++dummy_it;
  }

  delete[] mempool;

  if ( dummy_it >= 50 ) {
    std::cerr<<"\n[ERROR] Harmonic Analysis is corrupt! More than 50 frequencies identified!";
  }

  return amodel;
}

/// @brief Filter out earthquakes from a model.
///
/// Given a time-series and a model, try to filter out (from the model) any
/// earthquakes that produce statisticaly insignificant offsets. We will start
/// with a no-earthquake model, and iteratively add all earthquakes included
/// in the a-priori (input) model only keeping the ones for which the offset
/// induced is statisticaly significant. By 'statisticaly significant' we mean
/// that the a-priori std. deviation and the a-posteriori std. deviation (i.e.
/// before anf after including the earthquake in the model) satisfy the
/// relation:
/// (a_priori_stddev/a_posteriori_stddev) - 1 > U_t
///
/// @parameter[in] ts    The input time-series (this will remain const).
/// @parameter[in] model The input, a-priori model; the final estimated 
///                      (returned) model will include all earthquakes from 
///                      this model that produce statisticaly significant
///                      offset values.
/// @parameter[in] Ut    Critical value for the test of significance (see
///                      description). Default value is .001
/// @return              A new model; besides the earthquakes, all other member
///                      variables will be indentical to the a-priori model. Yet,
///                      this model will only include a subset of the
///                      earthquakes; the ones with statisticaly significant
///                      offset values.
template<class T>
  ngpt::ts_model<T>
filter_earthquakes(ngpt::timeseries<T, ngpt::pt_marker>& ts,
    ngpt::ts_model<T>& model, double Ut=1e-3)
{
#ifdef DEBUG
  std::cout<<"\n[DEBUG] Filtering earthquake events based on their statistical significance";
#endif
  // get the list of earthquakes from the input model.
  std::vector<ngpt::md_earthquake<T>> erthqk_vec {model.earthquakes()};

  // no earthquakes; quick return
  if (!erthqk_vec.size()) return model;

  double stddev_a, // previous std. dev
         stddev_n; // next (this) std. dev

  // Initial estimate to get the approximate earthquake offsets (no outlier
  // marking)
  ts.qr_ls_solve(model, stddev_n, 1e-3, false, false);
  erthqk_vec = model.earthquakes();

  // Get the a-priori earthquake vector and sort it according to each 
  // earthquake's resulting offset (as estimated with the a-priori model).
  std::sort(erthqk_vec.begin(), erthqk_vec.end(),
      [](const ngpt::md_earthquake<T>& a,
        const ngpt::md_earthquake<T>& b)
      {return a.abs_offset() > b.abs_offset();}
      );
#ifdef DEBUG
  std::cout<<"\n[DEBUG] Initial list of earthquakes, based on estimated offset:";
  for (auto i : erthqk_vec) {
    std::cout<<"\n\tEpoch: "<<ngpt::strftime_ymd_hms(i.start())<<" offset: "<<i.offset();
  }
#endif

  // Initialize a new model (identical to the a-priori) with no earthquakes.
  ngpt::ts_model<ngpt::milliseconds> amodel{model},
                                     nmodel{model};
  nmodel.clear_earthquakes();
  amodel.clear_earthquakes();

  // initial guess is a no-earthquake model.
  ts.qr_ls_solve(amodel, stddev_a, 1e-3, false, false);
  // iterators to the (sorted) earthquakes vector.
  typename std::vector<ngpt::md_earthquake<T>>::iterator
      it = erthqk_vec.begin(),
      it_end = erthqk_vec.end();
  double factor;

  // Add all earthquakes from list to the model (iteratively) and check the
  // post-fit residuals; if needed add the earthquake to the (final) model.
#ifdef DEBUG
  std::cout<<"\n[DEBUG] Start testing for significant offsets:";
#endif
  for (; it != it_end; ++it) {
    // add (the next) earthquake to the model and check the residuals.
#ifdef DEBUG
    std::cout<<"\n\tTesting earthquake: "<<ngpt::strftime_ymd_hms(it->start())<<", Jump Value:"<<it->offset();
#endif
    nmodel.add_earthquake(*it);
    ts.qr_ls_solve(nmodel, stddev_n, 1e-3, false, false);
#ifdef DEBUG
    std::cout<<"\n\t\tResulting std. dev: "<<stddev_n<<" (alpha model std.dev: "<<stddev_a<<")";
#endif
    factor = stddev_a / stddev_n;
    auto sz = nmodel.earthquakes().size();
    const auto last_eq = nmodel.earthquakes()[sz-1];
#ifdef DEBUG
    std::cout<<"\n\t\tOffset estimated at: "<<last_eq.abs_offset();
#endif
    if ((factor-1e0) > Ut && last_eq.abs_offset()>1e-3) {
      amodel = nmodel;
      stddev_a = stddev_n;
      std::cout<<"\n\t[DEBUG] Earthquake at "<<ngpt::strftime_ymd_hms(it->start())
        <<" added to the model; factor is "<<factor
        <<", value is: "<<it->a1();
    } else {
      // nmodel.erase_earthquake_at( it->start() );
      std::cout<<"\n\t[DEBUG] Earthquake at "<<ngpt::strftime_ymd_hms(it->start())
        <<" removed from the model; factor is "
        <<factor<<", value is: "<<it->a1();
      nmodel = amodel;
      stddev_n = stddev_a;
    }
#ifdef DEBUG
    std::cout<<"\n\tEq's in model: "<<amodel.earthquakes().size();
#endif
  }

  return amodel;
}

template<class T>
  ngpt::ts_model<T>
  try_earthquakes(ngpt::timeseries<T, ngpt::pt_marker>& ts, 
    ngpt::ts_model<T>& model, double min_mag=0e0/*, const event_list<T>* evns=nullptr*/, double Ut=1e-3)
{
    
  // get the list of vectors from the input model.
  std::vector<ngpt::md_earthquake<T>> erthqk_vec {model.earthquakes()};

  // no earthquakes; quick return
  if ( !erthqk_vec.size() ) return model;
    
  double stddev_a, // previous std. dev
         stddev_n, // next (this) std. dev
         stddev_n1, stddev_n2;
    
  // Initial estimate to get the approximate earthquake offsets (no outlier
  // marking)
  ts.qr_ls_solve(model, stddev_a, 1e-3, false, false);

  // Get the a-priori earthquake vector and sort it according to each 
  // earthquake's resulting offset (as estimated with the a-priori model).
  try {
  std::sort(erthqk_vec.begin(), erthqk_vec.end(),
      [](const ngpt::md_earthquake<T>& a,
        const ngpt::md_earthquake<T>& b)
      {return a.abs_offset() > b.abs_offset();}
      );
  } catch (std::exception& e) {
    std::cerr<<"\n ---- Caught exception id std::sort ----";
    return model;
  }
  
  /*
  try {
  if (min_mag && evns) {
    erthqk_vec.erase(std::remove_if(erthqk_vec.begin(), erthqk_vec.end(),
          [min_mag, &evns](const ngpt::md_earthquake<T>& e){return e.magnitude(*evns) < min_mag;}),
        erthqk_vec.end());
  }
  } catch (std::exception& e) {
    std::cerr<<"\n ---- Caught exception id vec.erase ----";
    return model;
  }
  */
  if (min_mag) {
    std::cerr<<"\n[WARNING] Skipping minimum magnitude for now .... (try_earthquakes)";
  }

  // Initialize a new model (identical to the a-priori) with earthquakes.
  ngpt::ts_model<ngpt::milliseconds> amodel{model},
                                     nmodel{model},
                                     logmdl{model},
                                     expmdl{model};
    
  // iterators to the (sorted) earthquakes vector.
  typename std::vector<ngpt::md_earthquake<T>>::iterator
      it = erthqk_vec.begin(),
      it_end = erthqk_vec.end();
  double factor;

  // Add all earthquakes from list to the model (iteratively) and check the
  // post-fit residuals; if needed add the earthquake to the (final) model.
  double pre_stddev;
  std::cout<<"\n[DEBUG] Start testing for significant offsets:";
  for (; it != it_end; ++it) {
    // add (the next) earthquake to the model and check the residuals.
    double ofst = it->a1();
    logmdl = amodel;
    expmdl = amodel;
    logmdl.change_earthquake_at(it->start(), psd_model::log, ofst, 0.15e0);
    expmdl.change_earthquake_at(it->start(), psd_model::exp, ofst, 0.15e0);
    pre_stddev = std::numeric_limits<double>::max();
    for (int iter=0; iter<10; iter++) {
      ts.qr_ls_solve(logmdl, stddev_n1, 1e-3, false, false);
      if (std::abs(stddev_n1 - pre_stddev)>1e-6) {
        pre_stddev = stddev_n1;
      } else {
        break;
      }
    }
    pre_stddev = std::numeric_limits<double>::max();
    for (int iter=0; iter<10; iter++) {
      ts.qr_ls_solve(expmdl, stddev_n2, 1e-3, false, false);
      if (std::abs(stddev_n2 - pre_stddev)>1e-6) {
        pre_stddev = stddev_n2;
      } else {
        break;
      }
    }
    if (stddev_n1 < stddev_n2) {
      stddev_n = stddev_n1;
      nmodel = logmdl;
    } else {
      stddev_n = stddev_n2;
      nmodel = expmdl;
    }
    factor = stddev_a / stddev_n;
    if ((factor-1e0) > Ut && ) {
      amodel = nmodel;
      stddev_a = stddev_n;
      std::cout<<"\n\t[DEBUG] Earthquake at "<<ngpt::strftime_ymd_hms(it->start())
        <<" added to the model; factor is "<<factor
        <<", value is: "<<it->a1();
    } else {
      // nmodel.change_earthquake_at(it->start(), psd_model::pwl, ofst);
      std::cout<<"\n\t[DEBUG] Earthquake at "<<ngpt::strftime_ymd_hms(it->start())
        <<" removed from the model; factor is "
        <<factor<<", value is: "<<it->a1();
    }
  }

  return amodel;
}

}// namespace ngpt

#endif
