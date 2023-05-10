#ifndef OT_LIBERTY_DELAY_HPP_
#define OT_LIBERTY_DELAY_HPP_

#include <ot/headerdef.hpp>
#include <vector>
#include <complex>
#include <algorithm>
#include <math.h>
#include <iostream>

#define DATA_MODE 0
#define CONVOLVE_MODE 0
#define TIME_STEP 0.001
#define SHRINK_THRESHOLD 1e-6;

typedef std::complex<float> Complex;

enum min_max_mode
{
	Min,
	Max
};

namespace ot
{

	class Statisical_delay;

	std::vector<float> convolve(std::vector<float> &, std::vector<float> &);

	// inplace FFT
	// sign 1 => forward FFT
	// sign -1 => inverse FFT
	void FFT(std::vector<Complex> &, int);

	std::vector<float> fft_convolve(std::vector<float> &, std::vector<float> &);

	Statisical_delay Statisical_delay_mm(Statisical_delay &x, Statisical_delay &y, bool mode);

	std::ostream &operator<<(std::ostream &os, const Statisical_delay &obj);

	// Class Statisical_delay
	class Statisical_delay
	{

	public:
		Statisical_delay(float, float, float, float);
		Statisical_delay(const std::vector<float> &Dist, float St = 0.0);
		Statisical_delay(const Statisical_delay &x);
		// Statisical_delay(Statisical_delay &&) = default;
		Statisical_delay(float);
		Statisical_delay();

		const float nominal() const;
		const float meanshift() const;
		const float mean() const;
		const float stdev() const;
		const float variance() const;
		const float skew() const;
		const float max() const;
		const float min() const;
		const float accum(float) const;

		// for debugging
		void display();

		Statisical_delay operator*(float s) const;
		Statisical_delay operator+(float v) const;
		Statisical_delay operator-(float v) const;
		Statisical_delay operator+(const Statisical_delay &) const;
		Statisical_delay operator-(const Statisical_delay &) const;

		Statisical_delay &operator=(const float &);
		Statisical_delay &operator=(const Statisical_delay &);
		Statisical_delay &operator=(Statisical_delay &&) = default;

		void dist_norm();
		void shrink();
		void negate();
		void print_info() const;

		// const float upper() const {return accum(0.9987);}
		// const float lower() const {return accum(0.0013);}

	private:
		float _nominal;
		float _mean_shift;
		float _std;
		float _skew;

		// lien
		std::vector<float> dist;
		float st;
		float get_et()
		{
			assert(DATA_MODE == 1);
			return st + TIME_STEP * (dist.size() - 1);
		}

		float time_unit = 1; // 1: ns

		friend Statisical_delay Statisical_delay_mm(Statisical_delay &, Statisical_delay &, bool);
	};

}; // end of namespace ot. -----------------------------------------------------------------------

#endif
