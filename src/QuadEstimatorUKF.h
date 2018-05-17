#pragma once

#include "BaseQuadEstimator.h"
#include "matrix/math.hpp"
#include "Math/Quaternion.h"

using matrix::Vector;
using matrix::Matrix;
using matrix::SquareMatrix;

#include "Eigen/Dense"
#include "Eigen/SVD"
using Eigen::MatrixXf;
using Eigen::VectorXf;

class QuadEstimatorUKF : public BaseQuadEstimator
{
public:
  QuadEstimatorUKF(string config, string name);
  virtual ~QuadEstimatorUKF();

  virtual void Init();

  virtual void Predict(float dt, V3F accel, V3F gyro);

  void ComputeSigmaPointsAndWeights();

  // helper functions for Predict
  VectorXf PredictState(VectorXf curState, float dt, V3F accel, V3F gyro);
  MatrixXf GetRbgPrime(float roll, float pitch, float yaw);

  virtual void UpdateFromIMU(V3F accel, V3F gyro);
  virtual void UpdateFromGPS(V3F pos, V3F vel);
  virtual void UpdateFromBaro(float z) {};
	virtual void UpdateFromMag(float magYaw);

  static const int QUAD_UKF_NUM_STATES = 7;
  static const int QUAD_UKF_NUM_SIGMA_POINTS = 2 * QUAD_UKF_NUM_STATES + 1;

  // process covariance
	MatrixXf Q;

	// GPS measurement covariance
	MatrixXf R_GPS;

	// Magnetometer measurement covariance
	MatrixXf R_Mag;

  // attitude filter state
  float pitchEst, rollEst;
  float accelPitch, accelRoll; // raw pitch/roll angles as calculated from last accelerometer.. purely for graphing.
	V3F accelG;
	V3F lastGyro;

  // generic EKF update
  // z: measurement
  // R: observation error model covariance 
  // zFromX: measurement prediction based on current state
  void Update(VectorXf& z, MatrixXf& R, MatrixXf& zFromSigmaPoints);

  // UKF state and covariance
	VectorXf ukfState;
	MatrixXf ukfCov;

  // UKF sigma points
  MatrixXf ukfSigmaPoints;
  // UKF mean weights
  VectorXf ukfMeanWeights;
  // UKF covariance weights
  VectorXf ukfCovWeights;

  // UKF parameters
  // the following are free parameters
  float alpha;
  float beta;
  float kappa;
  // the following are derived parameters
  float gamma;
  float lambda;

  // params
  float attitudeTau;
  float dtIMU;

  // Access functions for graphing variables
  virtual bool GetData(const string& name, float& ret) const;
  virtual vector<string> GetFields() const;
  string _name;

	// error vs ground truth (trueError = estimated-actual)
	virtual void UpdateTrueError(V3F truePos, V3F trueVel, SLR::Quaternion<float> trueAtt);
	VectorXf trueError;
	float pitchErr, rollErr, maxEuler;

	float posErrorMag, velErrorMag;

	virtual V3F EstimatedPosition() 
	{
		return V3F(ukfState(0), ukfState(1), ukfState(2));
	}

	virtual V3F EstimatedVelocity()
	{
		return V3F(ukfState(3), ukfState(4), ukfState(5));
	}

	virtual Quaternion<float> EstimatedAttitude()
	{
		return Quaternion<float>::FromEuler123_RPY(rollEst, pitchEst, ukfState(6));
	}

	virtual V3F EstimatedOmega()
	{
		return lastGyro;
	}

	float CovConditionNumber() const;
};
