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

        num_particles = 20;
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
            updated_theta = it->theta + delta_theta;
            if (fabs(yaw_rate) > ALMOST_ZERO) {
                double v_by_yaw_rate = velocity/yaw_rate;
                delta_x = v_by_yaw_rate*(sin(updated_theta) - sin(it->theta));
                delta_y = v_by_yaw_rate*(-cos(updated_theta) + cos(it->theta));
            } else {
                delta_x = velocity*cos(it->theta)*delta_t;
                delta_y = velocity*sin(it->theta)*delta_t;
            }
            normal_distribution<double> dist_x(it->x + delta_x, std_pos[0]);
            normal_distribution<double> dist_y(it->y + delta_y, std_pos[1]);
	    normal_distribution<double> dist_theta(updated_theta, std_pos[2]);
            it->x = dist_x(gen);
            it->y = dist_y(gen);
            it->theta = dist_theta(gen);

        }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
        double min_dist, temp_dist;
        for (vector<LandmarkObs>::iterator it_o = observations.begin();
             it_o != observations.end(); ++it_o) {
            min_dist = dist(it_o->x, it_o->y, predicted.at(0).x, predicted.at(0).y);
            it_o->id = predicted.at(0).id;
            for (vector<LandmarkObs>::iterator it_l = predicted.begin() + 1;
                 it_l != predicted.end(); ++it_l) {
                temp_dist = dist(it_o->x, it_o->y, it_l->x, it_l->y);
                if (temp_dist < min_dist) {
                    min_dist = temp_dist;
                    it_o->id = it_l->id;
                }
            }
        }
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

        // Transform the measurement based on each particle
        double weight_temp; 
        weights.clear();
        for (vector<Particle>::iterator it_p = particles.begin();
             it_p != particles.end(); ++it_p) {
            // Prepare the vector of map landmarks
            std::vector<LandmarkObs> landmarks;
            for (vector<Map::single_landmark_s>::const_iterator it_l = map_landmarks.landmark_list.begin();
                 it_l != map_landmarks.landmark_list.end(); ++it_l) {
                if (dist(it_p->x, it_p->y, it_l->x_f, it_l->y_f) < sensor_range) {
                    LandmarkObs single_landmark;
                    single_landmark.id = it_l->id_i;
                    single_landmark.x = it_l->x_f;
                    single_landmark.y = it_l->y_f;
                    landmarks.push_back(single_landmark);
                }
            }
            
            std::vector<LandmarkObs> transformed_obs = transformObservation(observations, it_p->x, it_p->y, it_p->theta);
            dataAssociation(landmarks, transformed_obs);
            std::vector<int> associations;
            std::vector<double> sense_x;
            std::vector<double> sense_y;
            for (vector<LandmarkObs>::iterator it_o = transformed_obs.begin();
                 it_o != transformed_obs.end(); ++it_o) {
                associations.push_back(it_o->id);
                sense_x.push_back(it_o->x);
                sense_y.push_back(it_o->y);
            }
            SetAssociations(*it_p, associations, sense_x, sense_y);
            weight_temp = calculateWeight(transformed_obs, std_landmark, map_landmarks);
            it_p->weight = weight_temp;
            weights.push_back(weight_temp);
        }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
        
        std::initializer_list<double> weight_init;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::discrete_distribution<> d(weights.begin(), weights.end());
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
    for (vector<LandmarkObs>::const_iterator it = observations.begin();
         it != observations.end(); ++it) {
        LandmarkObs transformed;
        transformed.x = it->x*cos(theta) - it->y*sin(theta) + x;
        transformed.y = it->x*sin(theta) + it->y*cos(theta) + y;
        transformedObs.push_back(transformed);
    }
    return transformedObs;
}

double ParticleFilter::calculateWeight(const std::vector<LandmarkObs> &observations, double std_landmark[], const Map map_landmarks) {
   double probability = 1.0;
   double std_x_std_y = std_landmark[0]*std_landmark[1];
   double std_x_2 = std_landmark[0]*std_landmark[0];
   double std_y_2 = std_landmark[1]*std_landmark[1];
   double dif_x, dif_y, param; 
   for (vector<LandmarkObs>::const_iterator it = observations.begin();
        it != observations.end(); ++it) {
       dif_x = it->x - map_landmarks.landmark_list.at(it->id - 1).x_f;
       dif_y = it->y - map_landmarks.landmark_list.at(it->id - 1).y_f;
       param = dif_x*dif_x/2/std_x_2 + dif_y*dif_y/2/std_y_2;
       probability *= exp(-param)/2/M_PI/std_x_std_y;
   }
   return probability; 
}
