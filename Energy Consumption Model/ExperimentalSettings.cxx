#include "ExperimentalSettings.h"
#include <tuple>
#include <limits>
#include <cmath>
#include <algorithm> // For min/max functions
#include "nlohmann/json.hpp" // to read and work with JSON files 
#include <fstream> // for file I/O
#include <vector>
#include <string>



ExperimentalSettings::ExperimentalSettings() {
    // Initialize with default values or leave empty
    initializeParameters();
}

void ExperimentalSettings::initializeParameters() {
    // Set initial values as in the constructor
    lambda1 = deg2rad(32.5876); 
    tau1 = deg2rad(-96.935);
    lambda2 = deg2rad(32.802); 
    tau2 = deg2rad(-96.6058);
    // setBounds(V_cruise, lambda1, lambda2, tau1, tau2);

    waypointSets = {{lambda1, tau1}, {lambda2, tau2}};
    altitudeSet = {500}; //{500,1000,2000,3000};
    VcruiseSet = {55,57,58};//,50,51,52,53,54,55,56}; //{50.0,55,60,65};
}


void ExperimentalSettings::setMissionParameters(double lambda1, double tau1, double lambda2, double tau2) {
    this->lambda1 = lambda1;
    this->tau1 = tau1;
    this->lambda2 = lambda2;
    this->tau2 = tau2;
}

void ExperimentalSettings::setAircraftParameters(double m, double g, double R_Earth, double rho, double A_rotor, double sigma, double C_d_mean, double F_P, double kappa, double Omega) {
    this->m = m;
    this->g = g;
    this->R_Earth = R_Earth;
    this->rho = rho;
    this->A_rotor = A_rotor;
    this->sigma = sigma;
    this->C_d_mean = C_d_mean;
    this->F_P = F_P;
    this->kappa = kappa;
    this->Omega = Omega;
}

void ExperimentalSettings::setBounds(double V_cruise, double lambda1, double lambda2, double tau1, double tau2) {
    // Example logic for setting bounds, adjust as necessary
    velocityBounds = {V_cruise, V_cruise};
    psiBounds = {calculatePsi(V_cruise, lambda1, lambda2, tau1, tau2) - M_PI/6, 
                 calculatePsi(V_cruise, lambda1, lambda2, tau1, tau2) + M_PI/6};
    lambdaBounds = {std::min(lambda1, lambda2) -0.2, std::max(lambda1, lambda2) +0.2};
    tauBounds = {std::min(tau1, tau2) - 0.2, std::max(tau1, tau2) + 0.2};
    timeBounds = {5,
                  calculateTravelTime(lambda1, tau1, lambda2, tau2, V_cruise) + 1000};
}

std::pair<double, double> ExperimentalSettings::getVelocityBounds() const {
    return velocityBounds;
}

std::pair<double, double> ExperimentalSettings::getPsiBounds() const {
    return psiBounds;
}

std::pair<double, double> ExperimentalSettings::getLambdaBounds() const {
    return lambdaBounds;
}

std::pair<double, double> ExperimentalSettings::getTauBounds() const {
    return tauBounds;
}
std::pair<double, double> ExperimentalSettings::getTimeBounds() const {
    return timeBounds;
}

std::tuple<double, double, double, double> ExperimentalSettings::getMissionParameters() const {
    return std::make_tuple(lambda1, tau1, lambda2, tau2);
}

std::tuple<double, double, double, double, double, double, double, double, double, double> ExperimentalSettings::getAircraftParameters() const {
    return std::make_tuple(m, g, R_Earth, rho, A_rotor, sigma, C_d_mean, F_P, kappa, Omega);
}

void ExperimentalSettings::setWaypointSets(const std::vector<std::pair<double, double>>& waypoints) {
    waypointSets = waypoints;
}

void ExperimentalSettings::setAltitudeSet(const std::vector<double>& altitudes) {
    altitudeSet = altitudes;
}

void ExperimentalSettings::setVcruiseSet(const std::vector<double>& speeds) {
    VcruiseSet = speeds;
}

const std::vector<std::pair<double, double>>& ExperimentalSettings::getWaypointSets() const {
    return waypointSets;
}
const std::vector<double>& ExperimentalSettings::getAltitudeSet() const {
    return altitudeSet;
}
const std::vector<double>& ExperimentalSettings::getVcruiseSet() const {
    return VcruiseSet;
}

double ExperimentalSettings::calculateTravelTime(double lat1, double lon1, double lat2, double lon2, double velocity) {
    const double R = 6371000; // Earth's radius in meters
    // Convert latitude and longitude from radians to degrees
    double lat1_deg = lat1 * 180 / M_PI;
    double lon1_deg = lon1 * 180 / M_PI;
    double lat2_deg = lat2 * 180 / M_PI;
    double lon2_deg = lon2 * 180 / M_PI;

    // Calculate the great circle distance
    double dlat = deg2rad(lat2_deg - lat1_deg);
    double dlon = deg2rad(lon2_deg - lon1_deg);

    double a = sin(dlat / 2) * sin(dlat / 2) + cos(lat1) * cos(lat2) * sin(dlon / 2) * sin(dlon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    double distance = R * c;

    double travelTime = distance / velocity;
    return travelTime;
}


double ExperimentalSettings::calculatePsi(double V, double lambda1, double lambda2, double tau1, double tau2) {
    double numerator = sin(tau2 - tau1) * cos(lambda2);
    double denominator = sin(lambda2) * cos(lambda1) - sin(lambda1) * cos(lambda2) * cos(tau2 - tau1);

    if (denominator != 0.0) {
        double psi = atan2(numerator, denominator);
        if (psi < 0.0) {
            psi += 2.0 * M_PI;
        }
        return psi;
    } else {
        // Handle the case where the denominator is zero (undefined)
        return std::numeric_limits<double>::quiet_NaN(); 
    }
}

double ExperimentalSettings::deg2rad(double degrees) {
    return degrees * M_PI / 180;
}

double ExperimentalSettings::rad2deg(double radians) {
    return radians * 180.0 / M_PI;
}


std::pair<double, double> ExperimentalSettings::calculateWindComponents(double lambda, double tau) {
    // Constants for wind model
    const double a_speed = 2.60080096e+02, b_speed = -1.23048100e+02, c_speed = 5.82544218e+00, d_speed = 4.76339623e+04;
    const double a_dir = 7.01187622e+03, b_dir = -4.01908890e+03, c_dir = 1.46096591e+02, d_dir = 1.27238444e+06;
    const double radian2degree = 180.0 / M_PI;

    // Calculate wind speed
    double windSpeed = a_speed * tau * radian2degree + b_speed * lambda * radian2degree + c_speed * tau * radian2degree * lambda * radian2degree + d_speed;

    // Calculate wind direction in radians
    double windDirRadians = (a_dir * tau * radian2degree + b_dir * lambda * radian2degree + c_dir * tau * radian2degree * lambda * radian2degree + d_dir) * M_PI / 180.0;

    // Calculate wind components
    double W_n_temp = windSpeed * sin(windDirRadians);
    double W_e_temp = windSpeed * cos(windDirRadians);

    double minWind = -12.5;
    double maxWind = -11;

    double W_n = ((W_n_temp < minWind) ? minWind : ((W_n_temp > maxWind) ? maxWind : W_n_temp));
    double W_e = ((W_e_temp < minWind) ? minWind : ((W_e_temp > maxWind) ? maxWind : W_e_temp));

    return std::make_pair(W_n, W_e);
}

double ExperimentalSettings::calculateCruiseSpeed(double lambda, double tau) {
    // Define the waypoints and speed limits
    const double lambda1 = deg2rad(32.5876);
    const double tau1 = deg2rad(-96.935);
    const double lambda2 = deg2rad(32.802);
    const double tau2 = deg2rad(-96.6058);
    const double maxSpeed = 65.0;
    const double minSpeed = 55.0;

    // Calculate relative position between waypoints
    double totalDistance = sqrt(pow(lambda2 - lambda1, 2) + pow(tau2 - tau1, 2));
    double currentDistance = sqrt(pow(lambda - lambda1, 2) + pow(tau - tau1, 2));
    double progress = currentDistance / totalDistance;

    // Calculate speed based on trapezoidal pattern
    double speed = minSpeed + (maxSpeed - minSpeed) * progress;
    speed = std::min(std::max(speed, minSpeed), maxSpeed); // Ensure speed is within bounds

    return speed;
}
