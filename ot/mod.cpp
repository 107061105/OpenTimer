#include <ot/mod.hpp>


std::vector<double> convolve(std::vector<double>& a, std::vector<double>& b) {
    int n = a.size(), m = b.size();
    std::vector<double> res(n+m-1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            res[i+j] += a[i] * b[j];
        }
    }
    return res;
}

// inplace FFT
// sign 1 => forward FFT
// sign -1 => inverse FFT
void FFT(std::vector<Complex>& x, int sign) {
    int n = x.size();
    if (n <= 1) return;
    std::vector<Complex> even(n/2), odd(n/2);
    for (int i = 0; i < n/2; i++) {
        even[i] = x[2*i];
        odd[i] = x[2*i+1];
    }
    FFT(even, sign);
    FFT(odd, sign);
    Complex w(1, 0);
    Complex wn(cos(2*M_PI/n), sign*sin(2*M_PI/n));
    for (int i = 0; i < n/2; i++) {
        x[i] = even[i] + w*odd[i];
        x[i+n/2] = even[i] - w*odd[i];
        w *= wn;
    }
}

std::vector<double> fft_convolve(std::vector<double>& a, std::vector<double>& b) {
    int nn = a.size() + b.size() - 1;
    int n = 1;
    while (n < nn) n <<= 1; // find FFT min length
    std::vector<Complex> x(n), y(n);
    for (int i = 0; i < a.size(); i++) x[i] = Complex(a[i], 0);
    for (int i = 0; i < b.size(); i++) y[i] = Complex(b[i], 0);
    FFT(x, 1);
    FFT(y, 1);
    for (int i = 0; i < n; i++) x[i] *= y[i]; // multiplication
    FFT(x, -1); // inverse FFT
    std::vector<double> res(nn);
    for (int i = 0; i < nn; i++) res[i] = x[i].real() / n; // take the real part and normalize
    for (int i = nn; i < n; i++) if (x[i].real() != 0) std::cout << "error! fft convolution has non-zero value outside res_dist (not in n+m-1 values)";
    return res;
}


float_mod::float_mod(const std::vector<double>& Dist, double St) {
	dist = Dist;
	st = St;
}

float_mod::float_mod(double St) {
    std::vector<double> temp(1, 1);
	dist = temp;
	st = St;
}

float_mod::float_mod() {
    std::vector<double> temp(1, 1);
	dist = temp;
	st = 0.0;
}

float_mod& float_mod::operator=(const double& St) {
    std::vector<double> temp(1, 1);
	dist = temp;
	st = St;
    return *this;
}

float_mod& float_mod::operator=(const float_mod& x) {
	dist = x.dist;
	st = x.st;
    return *this;
}

float_mod float_mod::operator+(const float_mod& x) const {
    float_mod temp1 = *this;
    float_mod temp2 = x;
	double res_st = temp1.st + temp2.st;
	// std::cout << "convolution st: " << res_st << std::endl;

	std::vector<double> res;
	if (CONVOLVE_MODE == 0) {
		res = convolve(temp1.dist, temp2.dist);
	} else if (CONVOLVE_MODE == 1) {
		res = fft_convolve(temp1.dist, temp2.dist);
	}
	
	// std::cout << "convolution result: ";
	// for (int i = 0; i < res.size(); i++) {
	// 	std::cout << res[i] << " ";
	// }
	// std::cout << std::endl;

	float_mod temp(res, res_st);
	return temp;
}

float_mod float_mod::operator-(const float_mod& x) const {
    float_mod temp = x;
	temp.negate();
    float_mod res = *this + temp;
    return res;
}

float_mod float_mod::operator+(float x) const {
    float_mod temp = *this;
	temp.st += x;
    return temp;
}

float_mod float_mod::operator-(float x) const {
    float_mod temp = *this;
	temp.st -= x;
    return temp;
}

void float_mod::dist_norm() {
    double sum = 0;
    for (auto& n : dist) sum += n;
	for (auto& n : dist) n /= sum;
}

void float_mod::negate() {
	// set start time to negative end time
	st = -get_et();
    // reverse the distribution
	std::reverse(dist.begin(), dist.end());
}

void float_mod::shrink() {
	// get new start time and end time, skipping zeros (or low values)
    double threshold = SHRINK_THRESHOLD;
    int new_st_index = 0;
    while (dist[new_st_index] <= threshold) {
        new_st_index++;
        st += TIME_STEP;
    }
    dist.erase(dist.begin(), dist.begin() + new_st_index);
    while (dist.back() <= threshold) {
        dist.pop_back();
    }
    dist_norm();
}

void float_mod::print_info() const {
	std::cout << "st: " << st << std::endl;
	std::cout << "dist: ";
	for (int i = 0; i < dist.size(); i++) {
		std::cout << dist[i] << " ";
	}
	std::cout << std::endl;
}

const float float_mod::mean() const {
    // sum = E[X] w/o st
    double sum = 0;
    double p_sum = 0;
    for (int i = 0; i < dist.size(); i++) {
        sum += dist[i] * TIME_STEP * i;
        p_sum += dist[i];
    }
    sum /= p_sum;
    return sum + st;
}

const float float_mod::variance() const {
    // sum = E[X^2] with st
    double sum = 0;
    double p_sum = 0;
    for (int i = 0; i < dist.size(); i++) {
        sum += dist[i] * pow(TIME_STEP * i + st, 2);
        p_sum += dist[i];
    }
    sum /= p_sum;
    return sum - pow(mean(), 2);
}

const float float_mod::accum(float p = 0.5) const {
    if (p <= 0 || p >= 1) {
        std::cout << "error, in float_mod::accum(p), p < 0 || p > 1" << std::endl;
        return 0;
    }
    float_mod temp = *this;
    temp.dist_norm();
    double p_sum = 0;
    for (int i = 0; i < temp.dist.size(); i++) {
        p_sum += temp.dist[i];
        if (p_sum >= p) {
            return TIME_STEP * i + temp.st;
        }
    }
    std::cout << "error, in float_mod::accum(p), overflow" << std::endl;
    return 0;
}

float_mod float_mod::sacle_time(float s) {
    return time_unit * s;
}

float_mod float_mod_mm(float_mod& x, float_mod&y, bool mode) {
    double st = std::min(x.st, y.st);
    double et = std::max(x.get_et(), y.get_et());
    int length = std::round((et - st) / TIME_STEP) + 1;
    // std::cout << "min/max info (st, et, length): " << st << ", " << et << ", " << length << std::endl;

    int x_offset = std::round((x.st-st)/TIME_STEP);
    int y_offset = std::round((y.st-st)/TIME_STEP);
    std::vector<double> res_dist(length, 0);
    std::vector<int> index_order;
    double px_sum = 0; // x's probability sum
    double py_sum = 0;
    if (mode == Max) {
        for (int i = 0; i < length; i++) index_order.push_back(i);
    } else if (mode == Min) {
        for (int i = length - 1; i >= 0; i--) index_order.push_back(i);
    }
    for (const auto& i : index_order) {
        double xp = 0; // x's probability at this time
        double yp = 0;
        if (i >= x_offset && i - x_offset < x.dist.size()) xp = x.dist[i - x_offset];
        if (i >= y_offset && i - y_offset < x.dist.size()) yp = y.dist[i - y_offset];
        res_dist[i] = py_sum * xp + px_sum * yp + xp * yp;
        px_sum = px_sum + xp;
        py_sum = py_sum + yp;
    }

    float_mod res(res_dist, st);
    res.shrink();
    return res;
}

std::ostream& operator<<(std::ostream& os, const float_mod& obj) {
    // write obj to stream
    os << obj.mean();
    return os;
}
