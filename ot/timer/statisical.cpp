#include <ot/timer/pin.hpp>
#include <ot/timer/arc.hpp>
#include <ot/timer/net.hpp>
#include <ot/timer/test.hpp>

namespace ot {

float findQuantile(const std::vector<float>& samples) {
    // Check if the vector is not empty
    if (samples.empty()) {
        std::cerr << "Error: Empty vector." << std::endl;
        return std::numeric_limits<float>::quiet_NaN();
    }

    // Create a copy of the vector and sort it
    std::vector<float> sortedSamples = samples;
    std::sort(sortedSamples.begin(), sortedSamples.end());

    // Calculate the index corresponding to the 0.135% quantile
    size_t index = static_cast<size_t>((1.0 - 0.00135) * sortedSamples.size());

    // Return the value at the calculated index
    return sortedSamples[index];
}

}