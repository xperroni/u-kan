/*
 * Copyright (c) Helio Perroni Filho <xperroni@gmail.com>
 *
 * This file is part of U-KAN.
 *
 * U-KAN is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * U-KAN is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with U-KAN. If not, see <http://www.gnu.org/licenses/>.
 */

#include "model_ctrv.h"
#include "sensors.h"
#include "ukf.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace ukan;

using namespace std;

// Laser measurement position variance across the X axis.
static const double s2_px = 0.0225;

// Laser measurement position variance across the Y axis.
static const double s2_py = 0.02;

// Radar direction variance.
static const double s2_r = 0.1;

// Radar position variance.
static const double s2_d = 0.1;

// Radar speed variance.
static const double s2_v = 0.1;

// Linear acceleration variance.
static const double s2_a = 0.04;

// Radial acceleration variance.
static const double s2_u = 0.04;

// Initial state variance.
static const double s2_P0 = 0.01;

/**
 * @brief Make sure user has provided input and output files.
 */
void check_arguments(int argc, char* argv[]) {
  if (argc == 2) {
    cerr << "Please include an output file." << endl;
  }
  else if (argc > 3) {
    cerr << "Too many arguments." << endl;
  }

  if (argc != 3) {
    cerr << "Usage instructions: " << argv[0] << " path/to/input.txt output.txt" << endl;
    exit(EXIT_FAILURE);
  }
}

/**
 * @brief Try to open given file, quit program if unsuccessful.
 */
template<class T> void open(T &file, const char *path) {
  file.open(path);
  if (!file.is_open()) {
    cerr << "Cannot open file \"" << path << '"' << endl;
    exit(EXIT_FAILURE);
  }
}

/**
 * @brief Run program.
 */
int main(int argc, char* argv[]) {
  check_arguments(argc, argv);

  ifstream data;
  open(data, argv[1]);

  ofstream out;
  open(out, argv[2]);

  vector<Measurement> measurements;
  vector<State> ground_truth;

  // Initalize model.
  Model model = new ctrv::Model(s2_a, s2_u, s2_P0);

  // Initalize sensors array.
  Sensors sensors(model, s2_px, s2_py, s2_d, s2_r, s2_v);

  // Read input measurements and ground truth.
  for (;;) {
    Measurement z = sensors(data);
    State g = model->newState();
    data >> g;
    if (data.eof())
      break;

    measurements.push_back(z);
    ground_truth.push_back(g);
  }

  // Create a Kalman filter instance.
  UKF filter(model);

  // State estimates.
  vector<State> estimates;

  // Produce state estimates for loaded measurements and
  // record them along ground truth values.
  for (size_t k = 0, n = measurements.size(); k < n; ++k) {
    Measurement z = measurements[k];
    State x = filter(z);

    cerr << "x:" << endl << filter.x << endl;
    cerr << "P:" << endl << filter.P << endl;

    out << x << '\t' << z << '\t' << ground_truth[k] << endl;

    estimates.push_back(x);
  }

  // Compute overall estimate accuracy using Root Mean Squared Error (RMSE).
  cerr << "Accuracy - RMSE:" << endl << fixed << setprecision(2) << RMSE(estimates, ground_truth) << endl;

  out.close();
  data.close();

  return 0;
}
