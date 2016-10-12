#include <stdexcept>
#include <cmath>
#include <random>
#include <iostream>
constexpr double TWOPID       = 3.14159265*2;

void
avevar(double data[], unsigned long n, double& ave, double& var)
{
    unsigned long j;
    double s,ep;

    for (ave=0.0,j=0;j<n;j++) ave += data[j];
    ave /= n;
    var=ep=0.0;
    for (j=0;j<n;j++) {
        s=data[j]-ave;
        ep += s;
        var += s*s;
    }
    var=(var-ep*ep/n)/(n-1);
}

///
/// Given n data points with abscissas x[1..n] (which need not be equally
/// spaced) and ordinates y[1..n], and given a desired oversampling factor of
/// ofac (a typical value being 4 or larger), this routine fills array px[1..np]
/// with an increasing sequence of frequencies (not angular frequencies) up to
/// hifac times the “average” Nyquist frequency, and fills array py[1..np] with
/// the values of the Lomb normalized periodogram at those frequencies. The
/// arrays x and y are not altered. np, the dimension of px and py, must be
/// large enough to contain the output, or an error results. The routine also
/// returns jmax such that py[jmax] is the maximum element in py, and prob, an
/// estimate of the significance of that maximum against the hypothesis of
/// random noise. A small value of prob indicates that a significant periodic
/// signal is present.
///
/// \param[in]  x     size n
/// \param[in]  y     size n
/// \param[in]  n     size of x and y
/// \param[in]  ofac  oversampling factor
/// \param[in]  hifac px includes frequencies up to hifac times the “average”
///                   Nyquist frequency
/// \param[out] px    increasing sequence of frequencies (not angular), size np
/// \param[out] py    values of the Lomb normalized periodogram at px
///                   frequencies, size np
/// \param[in]  np    size of px and py arrays
/// \param[out] nout  the number of different frequencies returned by the
///                   function : ofac*hifac*N/2
/// \param[out] jmax  py[jmax] is the maximum element in py
/// \param[out] prob  an estimate of the significance of the maximum (i.e.
///                   py[jmax]) against the hypothesis of random noise
///
void
period(double x[], double y[], int n, double ofac, double hifac, double px[], 
        double py[], int np, int& nout, int& jmax, double& prob)
{
    double ave,c,cc,cwtau,effm,expy,pnow,pymax,s,ss,sumc,sumcy,sums,sumsh,
           sumsy,swtau,var,wtau,xave,xdif,xmax,xmin,yy;
    double arg,wtemp,*wi,*wpi,*wpr,*wr;
    
    nout = 0.5*ofac*hifac*n;
    if (nout > np) {
        throw std::out_of_range{"output arrays too short in period"};
    }

    // Get mean and variance of the input data.
    avevar(y, n, ave, var);
    if (var == 0.0) {
        throw std::runtime_error{"zero variance in period"};
    }

    // Allocate memory
    try {
        wi  = new double[n];
        wpi = new double[n];
        wpr = new double[n];
        wr  = new double [n];
    } catch (std::bad_alloc&) {
        throw 1;
    }

    // xmax = xmin = x[0];
    xmin = x[0];
    xmax = x[n-1];

    xdif = xmax - xmin;
    xave = 0.5*(xmax+xmin);
    pymax= 0.0;
    pnow = 1.0/(xdif*ofac); // starting frequency
    // Initialize values for the trigonometric recurrences at each data point.
    for (int j = 0; j < n; j++) { 
        arg    = TWOPID*((x[j]-xave)*pnow);
        wpr[j] = -2.0*sin(0.5*arg)*sin(0.5*arg);
        wpi[j] = std::sin(arg);
        wr[j]  = std::cos(arg);
        wi[j]  = wpi[j];
    }

    // Main loop over the frequencies to be evaluated
    for (int i = 0; i < nout; i++) {
        px[i] = pnow;
        sumsh = sumc = 0.0;
        // First, loop over the data to get τ and related quantities.
        for (int j = 0; j < n; j++) {
            c = wr[j];
            s = wi[j];
            sumsh += s*c;
            sumc  += (c-s)*(c+s);
        }
        wtau  = 0.5*std::atan2(2.0*sumsh, sumc);
        swtau = std::sin(wtau);
        cwtau = std::cos(wtau);
        sums  = sumc = sumsy = sumcy = 0.0;
        // Then, loop over the data again to get the periodogram value.
        for (int j = 0; j < n; j++) {
            s  = wi[j];
            c  = wr[j];
            ss = s*cwtau-c*swtau;
            cc = c*cwtau+s*swtau;
            sums += ss*ss;
            sumc += cc*cc;
            yy = y[j]-ave;
            sumsy += yy*ss;
            sumcy += yy*cc;
            // Update the trigonometric recurrences.
            wr[j] = ((wtemp=wr[j])*wpr[j]-wi[j]*wpi[j])+wr[j];
            wi[j] = (wi[j]*wpr[j]+wtemp*wpi[j])+wi[j];
        }
        py[i] = 0.5*(sumcy*sumcy/sumc+sumsy*sumsy/sums)/var;
        if (py[i] >= pymax)
            pymax=py[(jmax=i)];
        pnow += 1.0/(ofac*xdif); // The next frequency.
    }

    // Evaluate statistical significance of the maximum.
    expy = std::exp(-pymax);
    effm = 2.0*nout/ofac;
    prob = effm*expy;
    if (prob > 0.01)
        prob=1.0-std::pow(1.0-expy,effm);

    // Free memory
    delete[] wr;
    delete[] wpr;
    delete[] wpi;
    delete[] wi;

    // The end
    return;
}


constexpr double MJD_START    = 53005.5e0;
constexpr double DAYS_IN_YEAR = 365.25e0;
constexpr double D2PI         = 3.14159265*2;
constexpr int    N            = 365*7;

int main()
{
    double time[N], vals[N];
    double mean_epoch = MJD_START + N/2;

    std::vector<double> frequencies;
    frequencies.push_back(1/(DAYS_IN_YEAR/2));
    frequencies.push_back(1/(DAYS_IN_YEAR/4));
    
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::normal_distribution<double> d(0e0, .5e0);
    for (int i = 0; i < N; i++) {
        time[i] = MJD_START + (double)i;
        vals[i] = d(gen)
            + .3*std::sin(D2PI*frequencies[0]*(mean_epoch-time[i]))
            + .3*std::cos(D2PI*frequencies[0]*(mean_epoch-time[i]))
            + .8*std::sin(D2PI*frequencies[1]*(mean_epoch-time[i]))
            + .8*std::cos(D2PI*frequencies[1]*(mean_epoch-time[i]));
    }

    //for (int i = 0; i < N; i++) {
    //    std::cout<<"\n"<<time[i]<<" "<<vals[i];
    //}
    
    double ofac{4}, hifac{.999999999};
    int nout = 0.5*ofac*hifac*N + 1;
    double px[nout], py[nout], prob;
    int jmax;

    period(time, vals, N, ofac, hifac, px, py, nout, nout, jmax, prob);

    for (int i = 0; i < nout; ++i) {
        std::cout << "\n" << px[i] << " " << py[i];
    }
    std::cerr<<"\nFrequencies present: "<<frequencies[0]<<" and "<<frequencies[1];
    std::cerr<<"\nAs angular: "<<D2PI*frequencies[0]<<" and "<<D2PI*frequencies[1];
    std::cerr<<"\nMax frequency = "<<px[jmax]<<" with probability: "<<prob;
    std::cerr<<"\nAs yearly fraction: 1/"<<DAYS_IN_YEAR*px[jmax]<<"\n";

    return 0;
}
