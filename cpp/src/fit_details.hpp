#ifndef __NGPT_FIT_DETAILS_HPP__
#define __NGPT_FIT_DETAILS_HPP__

// Eigen headers
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/QR"

// ggdatetime headers
// #include "ggdatetime/dtcalendar.hpp"

// gtms headers
#include "timeseries.hpp"
#include "model.hpp"

namespace ngpt
{

template<class T, class F>
    void
    construct_ls_matrices(
            const ngpt::timeseries<T, F>& ts,
            const ngpt::ts_model<T>&      model,
            Eigen::MatrixXd&        A,
            Eigen::VectorXd&        b,
            datetime<T>*            central_epoch = nullptr)
{
    // is the timeseries valid ?
    assert( ts.epoch_ptr() != nullptr && ts.epochs() == ts.data_pts() );

    // Number of parameters to be estimated:
    std::size_t parameters = model.parameters();

    // Number of observations (ommiting the ones to be skipped)
    std::size_t observations = ts.data_pts() - ts.data_skipped();
        
    // Can we estimate ?
    if ( observations < parameters ) {
        throw std::runtime_error("[ERROR] Too few observations to perform LS"
                "\n\tNumber of parameters:   "+std::to_string(parameters)+"\n\tNumber of observations: "+std::to_string(observations));
    }

    // Set the central epoch. If given, just assign it, else
    // compute the mean epoch; all deltatimes are computed as differences
    // from this (mean) epoch.
    double mean_epoch { central_epoch ? 
                       *central_epoch : 
                        ts.central_epoch().as_mjd() };

    // TODO would it be better to fill this column-wise??
    // Go ahead and form the A and b matrices.
    // We're gonna need some variables ...
    double dt, weight;
    std::size_t idx{0},       // index (i.e. row of A and b matrices)
                counter{0};   // index of the current data point (m_data)
    datetime<T> current_epoch;// the current epoch
    timeseries_const_iterator<T, F> it     = ts.cbegin(),
                                    it_end = ts.cend();
    data_point<F> entry;

    for (it; it != it_end; ++it) {
        entry         = it.data();
        current_epoch = it.epoch();
        // only include data that are not marked as 'skip'
        if ( !entry.skip() ) {
            dt = current_epoch.as_mjd() - mean_epoch;
            weight = sigma0 / entry.sigma();
            model.assign_row(A, b, current_epoch, weight, dt, entry.value(), idx);
            ++idx;
        }
        ++counter;
    }
#ifdef DEBUG
        assert(counter == ts.data_pts());
#endif

    return;
}

} // end namespace ngpt

#endif
