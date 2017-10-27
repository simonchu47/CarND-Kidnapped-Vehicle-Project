/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

#define ALMOST_ZERO 0.0001

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

        num_particles = 10;
        default_random_engine gen;
        // Creates a normal (Gaussian) distribution for x, y and theta
	normal_distribution<double> dist_x(x, std[0]);
        normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
        
        for (int i = 0; i < num_particles; ++i) {
                Particle p;
                p.id = i;
                p.x = dist_x(gen);
                p.y = dist_y(gen);
                p.theta = dist_theta(gen);
                p.weight = 1.0;
                particles.push_back(p);
        }

        last_max_possible_x = x;
        last_max_possible_y = y;
        
        is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
        default_random_engine gen;

        //double x_cur, y_cur, theta_cur;
        double delta_x, delta_y;
        double delta_theta = yaw_rate*delta_t;
        double updated_theta;

        
        for (vector<Particle>::iterator it = particles.begin();
             it != particles.end(); ++it) {
            updated_theta = (*it).theta + delta_theta;
            if (fabs(yaw_rate) > ALMOST_ZERO) {
                double v_by_yaw_rate = velocity/yaw_rate;
                delta_x = v_by_yaw_rate*(sin(updated_theta) - sin((*it).theta));
                delta_y = v_by_yaw_rate*(-cos(updated_theta) + cos((*it).theta));
            } else {
                delta_x = velocity*cos((*it).theta)*delta_t;
                delta_y = velocity*sin((*it).theta)*delta_t;
            }
            normal_distribution<double> dist_x((*it).x + delta_x, std_pos[0]);
            normal_distribution<double> dist_y((*it).y + delta_y, std_pos[1]);
	    normal_distribution<double> dist_theta(updated_theta, std_pos[2]);
            (*it).x = dist_x(gen);
            (*it).y = dist_y(gen);
            (*it).theta = dist_theta(gen);

        }

        /**
        double x_cur, y_cur, theta_cur;
        if (fabs(yaw_rate) > ALMOST_ZERO) {
            double v_by_yaw_rate = velocity/yaw_rate;
            for (int i = 0; i < particles.size(); ++i) {
                x_cur = particles.at(i).x;
                y_cur = particles.at(i).y;
                theta_cur = particles.at(i).theta;
                updated_theta = theta_cur + delta_theta;
                delta_x = v_by_yaw_rate*(sin(updated_theta) - sin(theta_cur));
                delta_y = v_by_yaw_rate*(-cos(updated_theta) + cos(theta_cur));
                normal_distribution<double> dist_x(x_cur + delta_x, std_pos[0]);
                normal_distribution<double> dist_y(y_cur + delta_y, std_pos[1]);
	        normal_distribution<double> dist_theta(updated_theta, std_pos[2]);
                particles.at(i).x = dist_x(gen);
                particles.at(i).y = dist_y(gen);
                particles.at(i).theta = dist_theta(gen);
            }
        } else {
            //TODO: when yaw rate is almost zero...
            for (int i = 0; i < particles.size(); ++i) { 
                x_cur = particles.at(i).x;
                y_cur = particles.at(i).y;
                theta_cur = particles.at(i).theta;
                updated_theta = theta_cur + delta_theta;
                delta_x = velocity*cos(theta_cur)*delta_t;
                delta_y = velocity*sin(theta_cur)*delta_t;
                normal_distribution<double> dist_x(x_cur + delta_x, std_pos[0]);
                normal_distribution<double> dist_y(y_cur + delta_y, std_pos[1]);
	        normal_distribution<double> dist_theta(updated_theta, std_pos[2]);
                particles.at(i).x = dist_x(gen);
                particles.at(i).y = dist_y(gen);
                particles.at(i).theta = dist_theta(gen);
            }
        }**/
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
        //double obs_x, obs_y;
        //double landmark_x, landmark_y;
        double min_dist, temp_dist;
        for (vector<LandmarkObs>::iterator it_o = observations.begin();
             it_o != observations.end(); ++it_o) {
            min_dist = dist((*it_o).x, (*it_o).y, predicted.at(0).x, predicted.at(0).y);
            (*it_o).id = predicted.at(0).id;
            for (vector<LandmarkObs>::iterator it_l = predicted.begin() + 1;
                 it_l != predicted.end(); ++it_l) {
                temp_dist = dist((*it_o).x, (*it_o).y, (*it_l).x, (*it_l).y);
                if (temp_dist < min_dist) {
                    min_dist = temp_dist;
                    (*it_o).id = (*it_l).id;
                }

            }
        }

        /**
        double obs_x, obs_y;
        double landmark_x, landmark_y;
        for (int i = 0; i < observations.size(); ++i) {
            obs_x = observations.at(i).x;
            obs_y = observations.at(i).y;
            landmark_x = predicted.at(0).x;
            landmark_y = predicted.at(0).y;
            min_dist = dist(obs_x, obs_y, landmark_x, landmark_y);
            observations.at(i).id = predicted.at(0).id;
            for(int j = 1; j < predicted.size(); ++j) {
                landmark_x = predicted.at(j).x;
                landmark_y = predicted.at(j).y;
                temp_dist = dist(obs_x, obs_y, landmark_x, landmark_y);
                if (temp_dist < min_dist) {
                    min_dist = temp_dist;
                    observations.at(i).id = predicted.at(j).id;
                }
            }
        }**/
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
        
        // Prepare the vector of map landmarks
        //std::vector<LandmarkObs> landmarks;
        if (landmarks.empty()) {
            LandmarkObs single_landmark;
            int index = 0;
            for (int i = 0; i < map_landmarks.landmark_list.size(); ++i) {
                //single_landmark.id = map_landmarks.landmark_list.at(i).id_i;
                //if (dist(last_max_possible_x, last_max_possible_y, map_landmarks.landmark_list.at(i).x_f, map_landmarks.landmark_list.at(i).y_f) < sensor_range) {
                    single_landmark.id = index;
                    single_landmark.x = map_landmarks.landmark_list.at(i).x_f;
                    single_landmark.y = map_landmarks.landmark_list.at(i).y_f;
                    landmarks.push_back(single_landmark);
                    index++;
                //}
            }
        }
        

        if (!landmarks.empty()) {
            // Transform the measurement based on each particle
            //std::vector<double> weight_raw;
            double weight_temp; 
            //double sum_weight_raw = 0.0;
            double max_weight = -1.0;
            weights.clear();
            for (int i = 0; i < particles.size(); ++i) {
                std::vector<LandmarkObs> transformed_obs = transformObservation(observations, particles.at(i).x, particles.at(i).y, particles.at(i).theta);
                dataAssociation(landmarks, transformed_obs);
                //cout << transformed_obs.at(0).x << " ," << transformed_obs.at(0).y << endl;
                weight_temp = calculateWeight(transformed_obs, std_landmark, landmarks);
                //weight_raw.push_back(calculateWeight(transformed_obs, std_landmark, landmarks));
                //weight_raw.push_back(weight_temp);
                //sum_weight_raw += weight_temp;
                if (weight_temp > max_weight) {
                    max_weight = weight_temp;
                    last_max_possible_x = particles.at(i).x;
                    last_max_possible_y = particles.at(i).y;
                }
                particles.at(i).weight = weight_temp;
                weights.push_back(weight_temp);
            }
        }

        /**
        double sum_weight_raw = 0.0;
        for (int i = 0; i < weight_raw.size(); ++i) {
            sum_weight_raw += weight_raw.at(i);
        }**/

        /**
        weights.clear();
        for (int i = 0; i < particles.size(); ++i) { 
            particles.at(i).weight = weight_raw.at(i)/sum_weight_raw;
            weights.push_back(particles.at(i).weight);
        }**/
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
        
        std::initializer_list<double> weight_init;
        std::random_device rd;
        std::mt19937 gen(rd());
        //std::discrete_distribution<> d(weight_init);
        std::discrete_distribution<> d(weights.begin(), weights.end());
        //std::discrete_distribution<> d({0.4, 0.1, 0.4, 0.1});

        std::vector<Particle> resample;
        for (int i = 0; i < num_particles; ++i) {
           resample.push_back(particles.at(d(gen)));
        }
        
        particles = resample;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

std::vector<LandmarkObs> ParticleFilter::transformObservation(const std::vector<LandmarkObs> &observations, double x, double y, double theta) {
    vector<LandmarkObs> transformedObs;
    LandmarkObs transformed, obs;    
    for (int i = 0; i < observations.size(); ++i) {
        obs = observations.at(i);
        transformed.x = obs.x*cos(theta) - obs.y*sin(theta) + x;
        transformed.y = obs.x*sin(theta) + obs.y*cos(theta) + y;
        transformedObs.push_back(transformed);
    }
    return transformedObs;
}

double ParticleFilter::calculateWeight(const std::vector<LandmarkObs> &observations, double std_landmark[], const std::vector<LandmarkObs> &map_landmarks) {
   double probability = 1.0;
   double std_x_std_y = std_landmark[0]*std_landmark[1];
   double std_x_2 = std_landmark[0]*std_landmark[0];
   double std_y_2 = std_landmark[1]*std_landmark[1];
   double dif_x, dif_y, param;
   LandmarkObs obs;
   for (int i = 0; i < observations.size(); ++i) {
       obs = observations.at(i);
       dif_x = obs.x - map_landmarks.at(obs.id).x;
       dif_y = obs.y - map_landmarks.at(obs.id).y;
       param = dif_x*dif_x/2/std_x_2 + dif_y*dif_y/2/std_y_2;
       probability *= exp(-param)/2/M_PI/std_x_std_y;
   }
   return probability; 
}
