#include <ot/liberty/delay.hpp>
#include <ot/static/logger.hpp>
#include <ot/headerdef.hpp>
#include <cmath>

namespace ot
{

    // Constructor
    Statisical_delay::Statisical_delay(float nom, float ms, float std, float skew)
    {
        assert(DATA_MODE == 0);
        _nominal = nom;
        _mean_shift = ms;
        _std = std;
        _skew = skew;
    }

    Statisical_delay::Statisical_delay(const std::vector<float> &Dist, float St)
    {
        assert(DATA_MODE == 1);
        dist = Dist;
        st = St;
    }

    Statisical_delay::Statisical_delay(const Statisical_delay &x)
    {
        if (DATA_MODE == 0)
        {
            _nominal = x._nominal;
            _mean_shift = x._mean_shift;
            _std = x._std;
            _skew = x._skew;
        }
        else if (DATA_MODE == 1)
        {
            dist = x.dist;
            st = x.st;
        }
    }

    Statisical_delay::Statisical_delay(float St)
    {
        if (DATA_MODE == 0)
        {
            _nominal = St;
            _mean_shift = 0;
            _std = 0;
            _skew = 0;
        }
        else if (DATA_MODE == 1)
        {
            std::vector<float> temp(1, 1);
            dist = temp;
            st = St;
        }
    }

    Statisical_delay::Statisical_delay()
    {
        if (DATA_MODE == 0)
        {
            _nominal = 0;
            _mean_shift = 0;
            _std = 0;
            _skew = 0;
        }
        else if (DATA_MODE == 1)
        {
            std::vector<float> temp(1, 1);
            dist = temp;
            st = 0;
        }
    }

    // for debugging
    void Statisical_delay::display()
    {
        if (DATA_MODE == 0)
        {
            OT_LOGD("Statisical delay\n"
                    "DATA_MODE: ",
                    DATA_MODE, "\n",
                    "nom:  ", nominal(), "\n",
                    "ms:   ", meanshift(), "\n",
                    "std:  ", stdev(), "\n",
                    "skew: ", skew(), "\n\n");
        }
        else if (DATA_MODE == 1)
        {
            OT_LOGD("Statisical delay\n"
                    "DATA_MODE: ",
                    DATA_MODE, "\n",
                    "st:   ", st, "\n",
                    "ms:   ", meanshift(), "\n",
                    "std:  ", stdev(), "\n",
                    "skew: ", skew(), "\n\n");
        }
    }

    // Operator
    Statisical_delay Statisical_delay::operator*(float s) const
    {
        if (DATA_MODE == 0)
        {
            return Statisical_delay(_nominal * s, _mean_shift * s, _std * s, _skew * s);
        }
        else if (DATA_MODE == 1)
        {
            Statisical_delay temp(*this);
            temp.time_unit = temp.time_unit * s;
            return temp;
        }
    }

    // Operator
    Statisical_delay Statisical_delay::operator+(const Statisical_delay &x) const
    {
        if (DATA_MODE == 0)
        {
            float nominal = _nominal + x.nominal();
            float meanshift = _mean_shift + x.meanshift();
            float stdev = sqrt(_std * _std + x.stdev() * x.stdev());
            float skew = _skew + x.skew();

            return Statisical_delay(nominal, meanshift, stdev, skew);
        }
        else if (DATA_MODE == 1)
        {
            Statisical_delay temp1 = *this;
            Statisical_delay temp2 = x;
            float res_st = temp1.st + temp2.st;
            // std::cout << "convolution st: " << res_st << std::endl;

            std::vector<float> res;
            if (CONVOLVE_MODE == 0)
            {
                res = convolve(temp1.dist, temp2.dist);
            }
            else if (CONVOLVE_MODE == 1)
            {
                res = fft_convolve(temp1.dist, temp2.dist);
            }

            // std::cout << "convolution result: ";
            // for (int i = 0; i < res.size(); i++) {
            // 	std::cout << res[i] << " ";
            // }
            // std::cout << std::endl;

            Statisical_delay temp(res, res_st);
            return temp;
        }
    }

    // Operator
    Statisical_delay Statisical_delay::operator-(const Statisical_delay &x) const
    {
        if (DATA_MODE == 0)
        {
            float nominal = _nominal - x.nominal();
            float meanshift = _mean_shift - x.meanshift();
            float stdev = sqrt(_std * _std + x.stdev() * x.stdev());
            float skew = _skew - x.skew();

            return Statisical_delay(nominal, meanshift, stdev, skew);
        }
        else if (DATA_MODE == 1)
        {
            Statisical_delay temp = x;
            temp.negate();
            Statisical_delay res = *this + temp;
            return res;
        }
    }

    // Operator
    Statisical_delay Statisical_delay::operator+(float f) const
    {
        if (DATA_MODE == 0)
        {
            return Statisical_delay(_nominal + f, _mean_shift, _std, _skew);
        }
        else if (DATA_MODE == 1)
        {
            Statisical_delay temp = *this;
            temp.st += f;
            return temp;
        }
    }

    // Operator
    Statisical_delay Statisical_delay::operator-(float f) const
    {
        if (DATA_MODE == 0)
        {
            return Statisical_delay(_nominal - f, _mean_shift, _std, _skew);
        }
        else if (DATA_MODE == 1)
        {
            Statisical_delay temp = *this;
            temp.st -= f;
            return temp;
        }
    }

    std::vector<float> convolve(std::vector<float> &a, std::vector<float> &b)
    {
        int n = a.size(), m = b.size();
        std::vector<float> res(n + m - 1);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                res[i + j] += a[i] * b[j];
            }
        }
        return res;
    }

    // inplace FFT
    // sign 1 => forward FFT
    // sign -1 => inverse FFT
    void FFT(std::vector<Complex> &x, int sign)
    {
        int n = x.size();
        if (n <= 1)
            return;
        std::vector<Complex> even(n / 2), odd(n / 2);
        for (int i = 0; i < n / 2; i++)
        {
            even[i] = x[2 * i];
            odd[i] = x[2 * i + 1];
        }
        FFT(even, sign);
        FFT(odd, sign);
        Complex w(1, 0);
        Complex wn(cos(2 * M_PI / n), sign * sin(2 * M_PI / n));
        for (int i = 0; i < n / 2; i++)
        {
            x[i] = even[i] + w * odd[i];
            x[i + n / 2] = even[i] - w * odd[i];
            w *= wn;
        }
    }

    std::vector<float> fft_convolve(std::vector<float> &a, std::vector<float> &b)
    {
        int nn = a.size() + b.size() - 1;
        int n = 1;
        while (n < nn)
            n <<= 1; // find FFT min length
        std::vector<Complex> x(n), y(n);
        for (int i = 0; i < a.size(); i++)
            x[i] = Complex(a[i], 0);
        for (int i = 0; i < b.size(); i++)
            y[i] = Complex(b[i], 0);
        FFT(x, 1);
        FFT(y, 1);
        for (int i = 0; i < n; i++)
            x[i] *= y[i]; // multiplication
        FFT(x, -1);       // inverse FFT
        std::vector<float> res(nn);
        for (int i = 0; i < nn; i++)
            res[i] = x[i].real() / n; // take the real part and normalize
        for (int i = nn; i < n; i++)
            if (x[i].real() != 0)
                std::cout << "error! fft convolution has non-zero value outside res_dist (not in n+m-1 values)";
        return res;
    }

    Statisical_delay &Statisical_delay::operator=(const float &St)
    {
        if (DATA_MODE == 0)
        {
            _nominal = St;
            _mean_shift = 0;
            _std = 0;
            _skew = 0;
            return *this;
        }
        else if (DATA_MODE == 1)
        {
            std::vector<float> temp(1, 1);
            dist = temp;
            st = St;
            return *this;
        }
    }

    Statisical_delay &Statisical_delay::operator=(const Statisical_delay &x)
    {
        if (DATA_MODE == 0)
        {
            _nominal = x._nominal;
            _mean_shift = x._mean_shift;
            _std = x._std;
            _skew = x._skew;
            return *this;
        }
        else if (DATA_MODE == 1)
        {
            dist = x.dist;
            st = x.st;
            return *this;
        }
    }

    void Statisical_delay::dist_norm()
    {
        assert(DATA_MODE == 1);
        float sum = 0;
        for (auto &n : dist)
            sum += n;
        for (auto &n : dist)
            n /= sum;
    }

    void Statisical_delay::negate()
    {
        assert(DATA_MODE == 1);
        // set start time to negative end time
        st = -get_et();
        // reverse the distribution
        std::reverse(dist.begin(), dist.end());
    }

    void Statisical_delay::shrink()
    {
        assert(DATA_MODE == 1);
        // get new start time and end time, skipping zeros (or low values)
        float threshold = SHRINK_THRESHOLD;
        int new_st_index = 0;
        while (dist[new_st_index] <= threshold)
        {
            new_st_index++;
            st += TIME_STEP;
        }
        dist.erase(dist.begin(), dist.begin() + new_st_index);
        while (dist.back() <= threshold)
        {
            dist.pop_back();
        }
        dist_norm();
    }

    void Statisical_delay::print_info() const
    {
        assert(DATA_MODE == 1);
        std::cout << "st: " << st << std::endl;
        std::cout << "dist: ";
        for (int i = 0; i < dist.size(); i++)
        {
            std::cout << dist[i] << " ";
        }
        std::cout << std::endl;
    }

    const float Statisical_delay::nominal() const
    {
        assert(DATA_MODE == 0);
        return _nominal;
    }

    const float Statisical_delay::meanshift() const
    {
        assert(DATA_MODE == 0);
        return _mean_shift;
    }

    const float Statisical_delay::mean() const
    {
        if (DATA_MODE == 0)
        {
            return _nominal + _mean_shift;
        }
        else if (DATA_MODE == 1)
        {
            // sum = E[X] w/o st
            float sum = 0;
            float p_sum = 0;
            for (int i = 0; i < dist.size(); i++)
            {
                sum += dist[i] * TIME_STEP * i;
                p_sum += dist[i];
            }
            sum /= p_sum;
            return sum + st;
        }
    }

    const float Statisical_delay::stdev() const
    {
        if (DATA_MODE == 0)
        {
            return _std;
        }
        else if (DATA_MODE == 1)
        {
            return sqrt(variance());
        }
    }

    const float Statisical_delay::variance() const
    {
        if (DATA_MODE == 0)
        {
            return _std * _std;
        }
        else if (DATA_MODE == 1)
        {
            // sum = E[X^2] with st
            float sum = 0;
            float p_sum = 0;
            for (int i = 0; i < dist.size(); i++)
            {
                sum += dist[i] * pow(TIME_STEP * i + st, 2);
                p_sum += dist[i];
            }
            sum /= p_sum;
            return sum - pow(mean(), 2);
        }
    }

    // Return skewness of delay distribution
    const float Statisical_delay::skew() const
    {
        assert(DATA_MODE == 0);
        return _skew;
    }

    const float Statisical_delay::max() const
    {
        if (DATA_MODE == 0)
        {
            return _nominal + _mean_shift + 3 * _std;
        }
        else if (DATA_MODE == 1)
        {
            return accum(0.9987);
        }
    }

    const float Statisical_delay::min() const
    {
        if (DATA_MODE == 0)
        {
            return _nominal + _mean_shift - 3 * _std;
        }
        else if (DATA_MODE == 1)
        {
            return accum(0.0013);
        }
    }

    const float Statisical_delay::accum(float p = 0.5) const
    {
        assert(DATA_MODE == 1);
        if (p <= 0 || p >= 1)
        {
            std::cout << "error, in Statisical_delay::accum(p), p < 0 || p > 1" << std::endl;
            return 0;
        }
        Statisical_delay temp = *this;
        temp.dist_norm();
        float p_sum = 0;
        for (int i = 0; i < temp.dist.size(); i++)
        {
            p_sum += temp.dist[i];
            if (p_sum >= p)
            {
                return TIME_STEP * i + temp.st;
            }
        }
        std::cout << "error, in Statisical_delay::accum(p), overflow" << std::endl;
        return 0;
    }

    Statisical_delay Statisical_delay_mm(Statisical_delay &x, Statisical_delay &y, bool mode)
    {
        assert(DATA_MODE == 1);
        float st = std::min(x.st, y.st);
        float et = std::max(x.get_et(), y.get_et());
        int length = std::round((et - st) / TIME_STEP) + 1;
        // std::cout << "min/max info (st, et, length): " << st << ", " << et << ", " << length << std::endl;

        int x_offset = std::round((x.st - st) / TIME_STEP);
        int y_offset = std::round((y.st - st) / TIME_STEP);
        std::vector<float> res_dist(length, 0);
        std::vector<int> index_order;
        float px_sum = 0; // x's probability sum
        float py_sum = 0;
        if (mode == Max)
        {
            for (int i = 0; i < length; i++)
                index_order.push_back(i);
        }
        else if (mode == Min)
        {
            for (int i = length - 1; i >= 0; i--)
                index_order.push_back(i);
        }
        for (const auto &i : index_order)
        {
            float xp = 0; // x's probability at this time
            float yp = 0;
            if (i >= x_offset && i - x_offset < x.dist.size())
                xp = x.dist[i - x_offset];
            if (i >= y_offset && i - y_offset < x.dist.size())
                yp = y.dist[i - y_offset];
            res_dist[i] = py_sum * xp + px_sum * yp + xp * yp;
            px_sum = px_sum + xp;
            py_sum = py_sum + yp;
        }

        Statisical_delay res(res_dist, st);
        res.shrink();
        return res;
    }

    std::ostream &operator<<(std::ostream &os, const Statisical_delay &obj)
    {
        // write obj to stream
        os << obj.mean();
        return os;
    }

} // end of namespace ot. -----------------------------------------------------------------------
