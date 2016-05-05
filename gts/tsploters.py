import datetime
import timeseries as ts
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def ts_plot(ts, y_erbar=False):
    ##  this is nice for large time-series
    mdates_series = mdates.date2num( ts.epoch_array )
    min_mdt = mdates.date2num( ts.epoch_array[0]  - datetime.timedelta( days=1.0 ) )
    max_mdt = mdates.date2num( ts.epoch_array[-1] + datetime.timedelta( days=1.0 ) )
    
    # a temporary ts with all outliers
    outliers_ts = ts.collect_outliers()
    mdates_out = mdates.date2num( outliers_ts.epoch_array )

    fig, axs = plt.subplots(3, sharex=True )
    axs[0].set_xlim( min_mdt, max_mdt )

    if not y_erbar:
        axs[0].plot_date( x=mdates_series, y=ts.x_array, fmt="-" )
    else:
        axs[0].errorbar(x=mdates_series, y=ts.x_array, xerr=None, fmt='-', 
                        yerr=ts.sx_array, ecolor='0.1')
    axs[0].set_title( 'Station %s'%ts.station )
    axs[0].set_ylabel( 'DNorth (m)' )
    axs[0].plot_date( x=mdates_out, y=outliers_ts.x_array, fmt="ro", markersize=4 )
    axs[0].set_ylim([np.amin(outliers_ts.x_array), np.amax(outliers_ts.x_array)])

    if not y_erbar:
       axs[1].plot_date( x=mdates_series, y=ts.y_array, fmt="-" )
    else:
        axs[1].errorbar(x=mdates_series, y=ts.y_array, xerr=None, fmt='-', 
                        yerr=ts.sy_array, ecolor='0.1')
    axs[1].set_ylabel( 'DEast (m)' )
    axs[1].plot_date( x=mdates_out, y=outliers_ts.y_array, fmt="ro", markersize=4)
    axs[1].set_ylim([np.amin(outliers_ts.y_array), np.amax(outliers_ts.y_array)])

    if not y_erbar:
        axs[2].plot_date( x=mdates_series, y=ts.z_array, fmt="-" )
    else:
        axs[2].errorbar(x=mdates_series, y=ts.z_array, xerr=None, fmt='-', 
                        yerr=ts.sz_array, ecolor='0.1')
    axs[2].set_ylabel( 'DUp (m)' )
    axs[2].plot_date( x=mdates_out, y=outliers_ts.z_array, fmt="ro", markersize=4 )
    axs[2].set_ylim([np.amin(outliers_ts.z_array), np.amax(outliers_ts.z_array)])

    axs[0].xaxis.set_major_locator( mdates.YearLocator() )
    axs[0].xaxis.set_major_formatter( mdates.DateFormatter('%Y') )
    axs[0].xaxis.set_minor_locator( mdates.MonthLocator(bymonth=[1,3,6,9,12]) )
    #axs[0].xaxis.set_minor_formatter( mdates.DateFormatter('%b') )

    fig.autofmt_xdate()
    plt.savefig( '%s.eps'%ts.station )
    plt.show()

    return
