
#include "psopt.h"
#include <cmath>
#include "ExperimentalSettings.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <map>
#include <iomanip> // for std::setw and std::left



// Function declarations
double deg2rad(double degrees);

// Define constants

// General paramters 
const double m = 6480/2.204;        // Aircraft mass
const double g = 9.81;       // Gravity
const double R_Earth =  6378100;  // Earth's radius
const double W_n = 0.0;      // Wind in the North direction // consider the time here when there is wind 
const double W_e = 0.0;      // Wind in the East direction
const double pi = 3.1415;

//Aircraft related paramters 
const double rho = 1.293e-3; // change this as a function of altitude 
const double A_rotor = 3.14 * 4 * 4; // R is 4m, so A_rotor = pi*R^2
const double sigma = 0.055; 
const double C_d_mean = 0.0089; // mean drag coefficient;
const double F_P = 0.97; // profile loss factor
const double kappa = 1.75; // proportionality constant
const double Omega = 30.12; // rotational velocity of the rotor blades

// Mission related paramters 
double lambda1; 
double tau1;
double lambda2;
double tau2;
double h;

// initial state related parameters 
double V_cruise;
double lambda_lower;
double lambda_upper;
double tau_lower;
double tau_upper;

// Calculate travel time
double travelTime_lower;
double travelTime_upper;
double psi_lower;
double psi_upper;


double deg2rad(double degrees) {
    return degrees * M_PI / 180;
}

// Instatiate the confuguration class globally 
ExperimentalSettings settings;


adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
   return 0;
}


adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                       adouble& time, adouble* xad, int iphase, Workspace* workspace) {

    adouble T_rotor = controls[0]/4;  
    adouble V = states[0];
    adouble theta = controls[1];
    adouble lambda = states[2];
    adouble tau = states[3];
    adouble alpha = theta;

    // Compute induced velocity in hover
    adouble vh = sqrt(T_rotor / (2 * rho * A_rotor)); 

    //computing  induced velocity using numerical method
    adouble vi = vh; 

    int i; 
    for (int i = 0; i < 50; i++){
        // eqn(9)
         vi = vh/sqrt((V * cos(alpha))*(V * cos(alpha)) + (V * sin(alpha) + vi)*(V * sin(alpha)) + vi); 
    }  

    // Compute the integrand of the performance index J
    adouble J_integrand = (4 * kappa * (T_rotor * vi) +
                        T_rotor * V * sin(alpha) + 
                        rho * A_rotor * pow(Omega*4, 3) * sigma * C_d_mean * F_P / 8.0)/1000000; // cost function will be energy in Mega Joules


    return J_integrand;
}



void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
   adouble x0 = initial_states[ 0 ];
   adouble x1 = initial_states[ 1 ];
   adouble x2 = initial_states[ 2 ];
   adouble x3 = initial_states[ 3 ];
   adouble y0 = final_states[ 0 ];
   // adouble y1 = final_states[ 1 ];
   adouble y2 = final_states[ 2 ];
   adouble y3 = final_states[ 3 ];


   e[ 0 ] = x0;
   e[ 1 ] = x1;
   e[ 2 ] = x2;
   e[ 3 ] = x3;
   e[ 4 ] = y0;
   // e[ 5 ] = y1;
   e[ 5 ] = y2;
   e[ 6 ] = y3;

}

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)

{    // Extract states and controls
    adouble V = states[0]; // V airspeed 
    adouble psi = states[1];
    adouble lambda = states[2];
    adouble tau = states[3];
    
    // adouble t = time;  
    adouble T = controls[0];
    adouble theta = controls[1];
    adouble phi = controls[2];

    // auto V = settings.calculateCruiseSpeed(lambda.value(), tau.value());
    // std::cout << "V_in: " << V_in << std::endl;



    // Calculate drag force (D)
    adouble D = 1.1984 * (rho * V * V) / 2.0;

    bool includeWind = false;// set this to true to include wind 

    if (includeWind){
        auto windComponents = settings.calculateWindComponents(lambda.value(), tau.value());
        // std::cout << "W_n: " << windComponents.first << ", W_e: " << windComponents.second << std::endl;

        double W_n = windComponents.first;
        double W_e = windComponents.second;
    }
    else {
        double W_n = 0;
        double W_e = 0;
    }
    // std::cout << "W_n = " << W_n <<  std::endl;
    // std::cout << "W_e = " << W_e << std::endl;

    // Define system dynamics
    derivatives[0] = (T*cos(phi)*sin(theta) - D)/m;
    derivatives[1] = (T*sin(phi))/(m*V);
    derivatives[2] = (V*cos(psi) + W_n)/(R_Earth + h);
    derivatives[3] = (V*sin(psi) + W_e)/((R_Earth + h)*cos(lambda));

    // Define path constraints
    path[0] = T*cos(phi)*cos(theta) - m*g;
    
    path[1] = (V*sin(psi) + W_e) * (sin(lambda2)*cos(lambda) - sin(lambda)*cos(lambda2)*cos(tau2 - tau))
              - (V*cos(psi) + W_n) * (sin(tau2 - tau)*cos(lambda2));
}



void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
  // No linkages as this is a single phase problem
}


int main() {

    const std::string filename = "flight_plans.json";
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file " + filename);
    }

    nlohmann::json json;
    file >> json;

    std::map<std::string, std::vector<std::pair<double, double>>> allWaypointSets;

    for (const auto& [flightPlanName, waypoints] : json.items()) {
        std::vector<std::pair<double, double>> waypointSets;
        for (const auto& wp : waypoints) {
            double latitude = deg2rad(wp.at(0));  // Accessing by index if waypoints are arrays
            double longitude = deg2rad(wp.at(1));
            waypointSets.emplace_back(latitude, longitude);
        }
        allWaypointSets[flightPlanName] = waypointSets;
    }

    struct SpeedCost {
        double speed;
        double cost;
        SpeedCost(double s, double c) : speed(s), cost(c) {}
    };


    // Declare optimal cost saving variables 
    struct SegmentResult {
        std::pair<std::pair<double, double>, std::pair<double, double>> waypoints; 
        std::vector<SpeedCost> speedCosts; // This will hold multiple speed-cost pairs
    };

    struct FlightPlanResult {
        std::string name;
        std::vector<SegmentResult> segments;
        std::map<double, double> totalCostPerSpeed; // Map to hold the total cost for each speed
    };

    std::vector<FlightPlanResult> allResults;


    // Iterating over all elements in allWaypointSets
    for (const auto& [flightPlanName, waypointSets] : allWaypointSets) {
        std::cout << "Flight Plan: " << flightPlanName << "\n";
        std::cout << "Waypoints:\n";
        
        FlightPlanResult currentFlightPlanResult;
        currentFlightPlanResult.name = flightPlanName;        
        // Loop through each waypoint set and optimize
        for (size_t i = 1; i < waypointSets.size(); ++i) {
        SegmentResult currentSegment;
        currentSegment.waypoints = std::make_pair(waypointSets[i - 1], waypointSets[i]);

        int fileIndex = 0; // Initialize a counter for file names

        lambda1 = waypointSets[i-1].first;
        tau1 = waypointSets[i-1].second;
        // Here, waypointSets[i] represents the next waypoint
        lambda2 = waypointSets[i].first;
        tau2 = waypointSets[i].second;

        std::cout << "Latitude 1: " << lambda1 << ", Longitude 1: " << lambda2 << "\n";
        std::cout << "Latitude 2: " << tau1 << ", Longitude 2: " << tau2 << "\n";

        for (const auto& altitude : settings.getAltitudeSet()) {
            h = altitude;

        for (double speed : settings.getVcruiseSet()) {
            double V_cruise = speed;

        // Set bounds
        settings.setBounds(V_cruise, lambda1, lambda2, tau1, tau2);
        // Retrieve bounds and use them in the problem setup
        auto velocityBounds = settings.getVelocityBounds();
        auto psiBounds = settings.getPsiBounds();
        auto lambdaBounds = settings.getLambdaBounds();
        auto tauBounds = settings.getTauBounds();
        auto timeBounds = settings.getTimeBounds();

        const double lambda_lower = lambdaBounds.first;
        const double lambda_upper = lambdaBounds.second;
        const double tau_lower = tauBounds.first;
        const double tau_upper = tauBounds.second;
        const double psi_lower = psiBounds.first;
        const double psi_upper = psiBounds.second;
        const double travelTime_lower = timeBounds.first;
        const double travelTime_upper = timeBounds.second;
        const double V_lower = velocityBounds.first;
        const double V_upper = velocityBounds.second;


        Alg  algorithm;
        Sol  solution;
        Prob problem;

        problem.name = "Aircraft Dynamics";


        // Generate unique file names based on lambda1, tau1, lambda2, and tau2
        std::ostringstream SummaryFileName;
        SummaryFileName << "eVTOLEnergy_" << fileIndex << "nowind" << "alt_"<< h << "speed_" << V_cruise<<".txt";

        problem.outfilename = SummaryFileName.str().c_str();
        problem.nphases = 1;
        problem.nlinkages = 0;

        // Define bounds and dimensions for the problem
        // This includes bounds on states, controls, and path constraints
        psopt_level1_setup(problem);
        problem.phases(1).nstates   		= 4;
        problem.phases(1).ncontrols 		= 3;
        problem.phases(1).nevents   		= 7;
        problem.phases(1).npath     		= 2;
        problem.phases(1).nodes                     << 120;

        psopt_level2_setup(problem, algorithm);

        problem.phases(1).bounds.lower.states <<    V_lower, psi_lower, lambda_lower, tau_lower; // bound should be based on the path 

        problem.phases(1).bounds.upper.states <<      V_upper, psi_upper, lambda_upper, tau_upper;

        problem.phases(1).bounds.lower.controls <<  20e3, -pi/8, -pi/8;

        problem.phases(1).bounds.upper.controls <<  45e3, pi/8, pi/8;

        problem.phases(1).bounds.lower.events << V_lower, psi_lower, lambda1, tau1, V_lower,  lambda2, tau2; // assign psi mannualy 
        
        problem.phases(1).bounds.upper.events << V_upper, psi_upper, lambda1, tau1, V_upper,  lambda2, tau2; 

        problem.phases(1).bounds.upper.path         << 0.0,0.0;

        problem.phases(1).bounds.lower.path         << 0.0,0.0;

        problem.phases(1).bounds.lower.StartTime    = 0.0;

        problem.phases(1).bounds.upper.StartTime    = 0.0;

        problem.phases(1).bounds.lower.EndTime      = travelTime_lower; 

        problem.phases(1).bounds.upper.EndTime      = travelTime_upper; 

        problem.integrand_cost 	= &integrand_cost;
        problem.endpoint_cost 	   = &endpoint_cost;
        problem.dae             	= &dae;
        problem.events 		      = &events;
        problem.linkages		      = &linkages;

        int nnodes    			               = problem.phases(1).nodes(0);

        MatrixXd x_guess    		            =  zeros(4,nnodes);

        x_guess.row(0) 			               =  V_lower * ones(1, nnodes);
        x_guess.row(1) 			               =  ((psi_lower+psi_upper)/2) * ones(1, nnodes);
        x_guess.row(2) 			               =  lambda1 * ones(1, nnodes);
        x_guess.row(3) 			               =  tau1 * ones(1, nnodes);

        MatrixXd control_guess                 = zeros(3,nnodes);
        control_guess.row(0)                  = 30e3 * ones(1, nnodes);
        control_guess.row(1)                  = (pi/8) * ones(1, nnodes);
        control_guess.row(2)                  = (pi/8) * ones(1, nnodes);//zeros(1, nnodes);

        problem.phases(1).guess.controls      = control_guess;
        problem.phases(1).guess.states        = x_guess;
        problem.phases(1).guess.time          = linspace(travelTime_lower,travelTime_upper,nnodes);

        algorithm.nlp_iter_max                = 1000;
        algorithm.nlp_tolerance               = 1.e-6;
        algorithm.nlp_method                  = "IPOPT";
        algorithm.scaling                     = "automatic";
        algorithm.derivatives                 = "automatic";
        algorithm.jac_sparsity_ratio          = 0.20;
        algorithm.collocation_method          = "Legendre";
        algorithm.diff_matrix                 = "central-differences";
        algorithm.mesh_refinement             = "automatic";
        algorithm.mr_max_increment_factor     = 0.3;
        algorithm.mr_max_iterations           = 3;
        algorithm.defect_scaling              = "jacobian-based";

        psopt(solution, problem, algorithm);

    ////////////////////////////////////////////////////////////////////////////
    ///////////  Extract relevant variables from solution structure   //////////
    ////////////////////////////////////////////////////////////////////////////


        MatrixXd x = solution.get_states_in_phase(1);
        MatrixXd u = solution.get_controls_in_phase(1);
        MatrixXd t = solution.get_time_in_phase(1);


    ////////////////////////////////////////////////////////////////////////////
    ///////////  Save solution data to files if desired ////////////////////////
    ////////////////////////////////////////////////////////////////////////////

       // Save(x,"x.dat");
       // Save(u,"u.dat");
       // Save(t,"t.dat");

        // Generate unique file names based on lambda1, tau1, lambda2, and tau2
        std::ostringstream xFileName, uFileName, tFileName;
        xFileName << "x_" << "no_wind" << "alt_"<< h << "speed_" << V_cruise<<".dat";
        uFileName << "u_" << "no_wind" << "alt_"<< h << "speed_" << V_cruise<<".dat";
        tFileName << "t_" << "no_wind" << "alt_"<< h << "speed_" << V_cruise<<".dat";

        Save(x, xFileName.str().c_str());
        Save(u, uFileName.str().c_str());
        Save(t, tFileName.str().c_str());
    ////////////////////////////////////////////////////////////////////////////
    ///////////  Plot some results if desired (requires gnuplot) ///////////////
    ////////////////////////////////////////////////////////////////////////////
if (false){
        plot(t,x.row(0),problem.name+": state V", "time (s)", "V","V");

        plot(t,x.row(1),problem.name+": state psi", "time (s)", "psi","psi");

        plot(t,x.row(2),problem.name+": state lambda", "time (s)", "lambda","lambda");

        plot(t,x.row(3),problem.name+": state tau", "time (s)", "tau","tau");

        plot(t,u.row(0),problem.name+": control T","time (s)", "control T", "T");

        plot(t,u.row(1),problem.name+": control theta","time (s)", "control theta", "theta");

        plot(t,u.row(2),problem.name+": control phi","time (s)", "control phi", "phi");

        plot(t,x.row(0),problem.name+": state V", "time (s)", "state","V",
                                            "pdf", "eVTOL_state1.pdf");

        plot(t,x.row(1),problem.name+": state psi", "time (s)", "state","psi",
                                            "pdf", "eVTOL_state2.pdf");

        plot(t,x.row(2),problem.name+": state lamabda", "time (s)", "state","lamabda",
                                            "pdf", "eVTOL_state3.pdf");

        plot(t,x.row(3),problem.name+": state tau", "time (s)", "state","tau",
                                            "pdf", "eVTOL_state4.pdf");

        plot(t,u.row(0),problem.name+": control T","time (s)", "control", "T",
                                            "pdf", "eVTOL_control1.pdf");

        plot(t,u.row(1),problem.name+": control theta","time (s)", "control", "theta",
                                            "pdf", "eVTOL_control2.pdf");

        plot(t,u.row(2),problem.name+": control phi","time (s)", "control", "phi",
                                            "pdf", "eVTOL_control3.pdf");
}
        // Increment the file index
        fileIndex++;
            
        double optimal_cost = solution.get_cost();/* calculate or retrieve your optimal_cost here */
        currentSegment.speedCosts.emplace_back(speed, optimal_cost);
        currentFlightPlanResult.totalCostPerSpeed[speed] += optimal_cost;

        // flightPlanOptimalCost += optimal_cost;
        // optimal_costs.push_back(optimal_cost);     // Store the cost

    }
        }
                currentFlightPlanResult.segments.push_back(currentSegment);
    }
    allResults.push_back(currentFlightPlanResult);

      }
    for (const auto& flightPlan : allResults) {
        std::cout << "Flight Plan Name: " << flightPlan.name << "\n";
        std::cout << std::left << std::setw(50) << "Segment (initial & final waypoints)"
                << std::setw(15) << "Speed"
                << "Optimal Cost\n";
        std::cout << std::string(80, '-') << "\n"; // Separator line

        for (const auto& segment : flightPlan.segments) {
            bool firstSpeedCost = true;
            for (const auto& speedCost : segment.speedCosts) {
                if (firstSpeedCost) {
                    std::cout << std::left << std::setw(50)
                            << "(" + std::to_string(segment.waypoints.first.first) + ", " 
                            + std::to_string(segment.waypoints.first.second) + ") -> ("
                            + std::to_string(segment.waypoints.second.first) + ", "
                            + std::to_string(segment.waypoints.second.second) + ")";
                    firstSpeedCost = false;
                } else {
                    std::cout << std::left << std::setw(50) << ""; // Empty for subsequent lines
                }
                std::cout << std::setw(15) << speedCost.speed << speedCost.cost << "\n";
            }
            std::cout << "\n";
        }
        // Print summary for the flight plan
        std::cout << "\nSummary for " << flightPlan.name << ":\n";
        std::cout << "First waypoint: (" << flightPlan.segments.front().waypoints.first.first << ", " 
                << flightPlan.segments.front().waypoints.first.second << ")\n";
        std::cout << "Last waypoint: (" << flightPlan.segments.back().waypoints.second.first << ", " 
                << flightPlan.segments.back().waypoints.second.second << ")\n";

        std::cout << std::left << std::setw(15) << "Speed"
                << "Total Optimal Cost\n";
        for (const auto& [speed, totalCost] : flightPlan.totalCostPerSpeed) {
            std::cout << std::setw(15) << speed
                    << totalCost << "\n";
        }
        std::cout << "\n";
        }
        std::ofstream outFile("output.csv");
        for (const auto& flightPlan : allResults) {
        // Write flight plan header
        outFile << "Flight Plan Name: " << flightPlan.name << "\n";

        // Write table headers
        outFile << "Segment (initial & final waypoints),Speed,Optimal Cost\n";

        for (const auto& segment : flightPlan.segments) {
            for (const auto& [speed, cost] : segment.speedCosts) {
                outFile << "\"" << segment.waypoints.first.first << ", " << segment.waypoints.first.second
                        << " -> " << segment.waypoints.second.first << ", " << segment.waypoints.second.second << "\","
                        << speed << "," << cost << "\n";
            }
        }

        // Write summary for the flight plan
        outFile << "Summary\n";
        outFile << "First waypoint," << flightPlan.segments.front().waypoints.first.first << ", " 
                << flightPlan.segments.front().waypoints.first.second << "\n";
        outFile << "Last waypoint," << flightPlan.segments.back().waypoints.second.first << ", " 
                << flightPlan.segments.back().waypoints.second.second << "\n";

        outFile << "Speed,Total Optimal Cost\n";
        for (const auto& [speed, totalCost] : flightPlan.totalCostPerSpeed) {
            outFile << speed << "," << totalCost << "\n";
        }
        outFile << "\n";
    }

    // Close the file
    outFile.close();
}
