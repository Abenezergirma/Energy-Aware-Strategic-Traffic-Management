#ifndef EXPERIMENTAL_SETTINGS_H
#define EXPERIMENTAL_SETTINGS_H

#include <vector>
#include <utility>
#include <cmath>

class ExperimentalSettings {
public:
    // Constructor
    ExperimentalSettings();

    // Setters for parameters
    void setMissionParameters(double lambda1, double tau1, double lambda2, double tau2);
    void setAircraftParameters(double m, double g, double R_Earth, double rho, double A_rotor, double sigma, double C_d_mean, double F_P, double kappa, double Omega);
    void setWaypointSets(const std::vector<std::pair<double, double>>& waypoints);
    void setAltitudeSet(const std::vector<double>& altitudes);
    void setVcruiseSet(const std::vector<double>& speeds);
    void initializeParameters(); 
    // Method to set bounds
    void setBounds(double V_cruise, double lambda1, double lambda2, double tau1, double tau2);

    // Getters for bounds
    std::pair<double, double> getVelocityBounds() const;
    std::pair<double, double> getPsiBounds() const;
    std::pair<double, double> getLambdaBounds() const;
    std::pair<double, double> getTauBounds() const;
    std::pair<double, double> getTimeBounds() const;
    // Getter methods for mission parameters
    std::tuple<double, double, double, double> getMissionParameters() const;
    // Getter methods for aircraft parameters
    std::tuple<double, double, double, double, double, double, double, double, double, double> getAircraftParameters() const;

    const std::vector<std::pair<double, double>>& getWaypointSets() const;
    const std::vector<double>& getAltitudeSet() const;
    const std::vector<double>& getVcruiseSet() const;
    std::pair<double, double> calculateWindComponents(double lambda, double tau);

    // Calculation functions
    double calculateTravelTime(double lat1, double lon1, double lat2, double lon2, double velocity);
    double calculatePsi(double V, double lambda1, double lambda2, double tau1, double tau2);
    double calculateCruiseSpeed(double lambda, double tau);

    static double deg2rad(double degrees);
    static double rad2deg(double radians);

    // Getters for parameters (if needed)

private:

    std::vector<std::pair<double, double>> waypointSets;
    std::vector<double> altitudeSet;
    std::vector<double> VcruiseSet;

    // Mission related parameters
    double lambda1, tau1, lambda2, tau2;
    // Aircraft parameters
    double m, g, R_Earth, rho, A_rotor, sigma, C_d_mean, F_P, kappa, Omega;

    // Bounds
    std::pair<double, double> velocityBounds;
    std::pair<double, double> psiBounds;
    std::pair<double, double> lambdaBounds;
    std::pair<double, double> tauBounds;
    std::pair<double, double> timeBounds;
};

#endif // EXPERIMENTAL_SETTINGS_H
