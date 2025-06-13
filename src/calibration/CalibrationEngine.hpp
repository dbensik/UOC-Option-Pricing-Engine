// CalibrationEngine.hpp
#ifndef CALIBRATION_ENGINE_HPP
#define CALIBRATION_ENGINE_HPP

#include <vector>
#include <iostream>
#include "UpOutCallOption.hpp"

namespace Calibration {

    void gridSearch(
        std::vector<double>& initialTheta,
        std::vector<double>& optimalTheta,
        const std::vector<double>& strikes,
        const std::vector<double>& barriers,
        const std::vector<std::vector<double>>& marketPrices,
        std::vector<std::vector<double>>& modelPrices);

    void getPriceGrid(
        const std::vector<double>& theta,
        const std::vector<double>& strikes,
        const std::vector<double>& barriers,
        std::vector<std::vector<double>>& modelPrices);

    double ssre(
        const std::vector<std::vector<double>>& modelPrices,
        const std::vector<std::vector<double>>& marketPrices);

    void ssre(
        const std::vector<std::vector<double>>& modelPrices,
        const std::vector<std::vector<double>>& marketPrices,
        std::vector<double>& diffVector);

    double objective(std::vector<double>& simplexTheta);

    int nelderMeadExample();
    int nelderMead1(const std::vector<double>& theta);
    int nelderMead2(const std::vector<double>& theta);

} // namespace Calibration

#endif
