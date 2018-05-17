#include "Common.h"
#include "QuadEstimatorUKF.h"
#include "Utility/SimpleConfig.h"
#include "Utility/StringUtils.h"
#include "Math/Quaternion.h"

#include <iostream>
using namespace std;

using namespace SLR;

const int QuadEstimatorUKF::QUAD_UKF_NUM_STATES;
const int QuadEstimatorUKF::QUAD_UKF_NUM_SIGMA_POINTS;

QuadEstimatorUKF::QuadEstimatorUKF(string config, string name)
  : BaseQuadEstimator(config),
  Q(QUAD_UKF_NUM_STATES, QUAD_UKF_NUM_STATES),
  R_GPS(6, 6),
  R_Mag(1, 1),
  ukfState(QUAD_UKF_NUM_STATES),
  ukfCov(QUAD_UKF_NUM_STATES, QUAD_UKF_NUM_STATES),
  trueError(QUAD_UKF_NUM_STATES),
  ukfSigmaPoints(QUAD_UKF_NUM_STATES, QUAD_UKF_NUM_SIGMA_POINTS),
  ukfMeanWeights(QUAD_UKF_NUM_SIGMA_POINTS),
  ukfCovWeights(QUAD_UKF_NUM_SIGMA_POINTS)
{
  _name = name;
  Init();
}

QuadEstimatorUKF::~QuadEstimatorUKF()
{

}

void QuadEstimatorUKF::Init()
{
  ParamsHandle paramSys = SimpleConfig::GetInstance();

  paramSys->GetFloatVector(_config + ".InitState", ukfState);

  VectorXf initStdDevs(QUAD_UKF_NUM_STATES);
  paramSys->GetFloatVector(_config + ".InitStdDevs", initStdDevs);
  ukfCov.setIdentity();
  for (int i = 0; i < QUAD_UKF_NUM_STATES; i++)
  {
    ukfCov(i, i) = initStdDevs(i) * initStdDevs(i);
  }

  // complementary filter params
  attitudeTau = paramSys->Get(_config + ".AttitudeTau", .1f);
  dtIMU = paramSys->Get(_config + ".dtIMU", .002f);

  pitchEst = 0;
  rollEst = 0;
  
  // GPS measurement model covariance
  R_GPS.setZero();
  R_GPS(0, 0) = R_GPS(1, 1) = powf(paramSys->Get(_config + ".GPSPosXYStd", 0), 2);
  R_GPS(2, 2) = powf(paramSys->Get(_config + ".GPSPosZStd", 0), 2);
  R_GPS(3, 3) = R_GPS(4, 4) = powf(paramSys->Get(_config + ".GPSVelXYStd", 0), 2);
  R_GPS(5, 5) = powf(paramSys->Get(_config + ".GPSVelZStd", 0), 2);

  // magnetometer measurement model covariance
  R_Mag.setZero();
  R_Mag(0, 0) = powf(paramSys->Get(_config + ".MagYawStd", 0), 2);

  // load the transition model covariance
  Q.setZero();
  Q(0, 0) = Q(1, 1) = powf(paramSys->Get(_config + ".QPosXYStd", 0), 2);
  Q(2, 2) = powf(paramSys->Get(_config + ".QPosZStd", 0), 2);
  Q(3, 3) = Q(4, 4) = powf(paramSys->Get(_config + ".QVelXYStd", 0), 2);
  Q(5, 5) = powf(paramSys->Get(_config + ".QVelZStd", 0), 2);
  Q(6, 6) = powf(paramSys->Get(_config + ".QYawStd", 0), 2);
  Q *= dtIMU;

  // initialize sigma points matrix
  // for N states, we will generate 2 * N + 1 sigma points
  ukfSigmaPoints.setZero();
  ukfMeanWeights.setZero();
  ukfCovWeights.setZero();

  // initialize UKF parameters
  // In UKF there're 3 free parameters: alpha, beta and kappa
  // It's too hard for me to figure out the meaning of different set of values.
  // I'll used the set used in UKF exercise jupyter notebook, that is:
  // alpha = 1
  // beta = 2
  // kappa = 3 - N
  alpha = 1.0f;
  beta = 2.0f;
  kappa = 3 - QUAD_UKF_NUM_STATES;
  lambda = alpha * alpha * (QUAD_UKF_NUM_STATES + kappa) - QUAD_UKF_NUM_STATES;
  gamma = sqrtf(QUAD_UKF_NUM_STATES + lambda);

  rollErr = pitchErr = maxEuler = 0;
  posErrorMag = velErrorMag = 0;
}

void QuadEstimatorUKF::UpdateFromIMU(V3F accel, V3F gyro)
{
  // Improve a complementary filter-type attitude filter
  // 
  // Currently a small-angle approximation integration method is implemented
  // The integrated (predicted) value is then updated in a complementary filter style with attitude information from accelerometers
  // 
  // Implement a better integration method that uses the current attitude estimate (rollEst, pitchEst and ukfState(6))
  // to integrate the body rates into new Euler angles.
  //
  // HINTS:
  //  - there are several ways to go about this, including:
  //    1) create a rotation matrix based on your current Euler angles, integrate that, convert back to Euler angles
  //    OR 
  //    2) use the Quaternion<float> class, which has a handy FromEuler123_RPY function for creating a quaternion from Euler Roll/PitchYaw
  //       (Quaternion<float> also has a IntegrateBodyRate function, though this uses quaternions, not Euler angles)

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  // SMALL ANGLE GYRO INTEGRATION:
  // (replace the code below)
  // make sure you comment it out when you add your own code -- otherwise e.g. you might integrate yaw twice
  
  Quaternion<float> q = Quaternion<float>::FromEuler123_RPY(rollEst, pitchEst, ukfState(6));
  q.IntegrateBodyRate(gyro, dtIMU);
  
  V3D eulerRPY = q.ToEulerRPY();
  float predictedRoll = static_cast<float>(eulerRPY[0]);
  float predictedPitch = static_cast<float>(eulerRPY[1]);
  ukfState(6) = static_cast<float>(eulerRPY[2]);
  
  /*
  float e[9] = {
    1, sinf(rollEst) * tanf(pitchEst), cosf(rollEst) * tan(pitchEst),
    0, cosf(rollEst), - sinf(rollEst),
    0, sinf(rollEst) / cosf(pitchEst), cosf(rollEst) / cosf(pitchEst)
  };
  V3F eulerDerivatives = Mat3x3F(e) * gyro;
  V3F newEuler = V3F(rollEst, pitchEst, ukfState(6)) + eulerDerivatives * dtIMU;
  float predictedRoll = newEuler[0];
  float predictedPitch = newEuler[1];
  ukfState(6) = newEuler[2];
  */
  // normalize yaw to -pi .. pi
  if (ukfState(6) > F_PI) ukfState(6) -= 2.f*F_PI;
  if (ukfState(6) < -F_PI) ukfState(6) += 2.f*F_PI;

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  // CALCULATE UPDATE
  accelRoll = atan2f(accel.y, accel.z);
  accelPitch = atan2f(-accel.x, 9.81f);

  // FUSE INTEGRATION AND UPDATE
  rollEst = attitudeTau / (attitudeTau + dtIMU) * (predictedRoll)+dtIMU / (attitudeTau + dtIMU) * accelRoll;
  pitchEst = attitudeTau / (attitudeTau + dtIMU) * (predictedPitch)+dtIMU / (attitudeTau + dtIMU) * accelPitch;

  lastGyro = gyro;
}

void QuadEstimatorUKF::UpdateTrueError(V3F truePos, V3F trueVel, Quaternion<float> trueAtt)
{
  VectorXf trueState(QUAD_UKF_NUM_STATES);
  trueState(0) = truePos.x;
  trueState(1) = truePos.y;
  trueState(2) = truePos.z;
  trueState(3) = trueVel.x;
  trueState(4) = trueVel.y;
  trueState(5) = trueVel.z;
  trueState(6) = trueAtt.Yaw();

  trueError = ukfState - trueState;
  if (trueError(6) > F_PI) trueError(6) -= 2.f*F_PI;
  if (trueError(6) < -F_PI) trueError(6) += 2.f*F_PI;

  pitchErr = pitchEst - trueAtt.Pitch();
  rollErr = rollEst - trueAtt.Roll();
  maxEuler = MAX(fabs(pitchErr), MAX(fabs(rollErr), fabs(trueError(6))));

  posErrorMag = truePos.dist(V3F(ukfState(0), ukfState(1), ukfState(2)));
  velErrorMag = trueVel.dist(V3F(ukfState(3), ukfState(4), ukfState(5)));
}

void QuadEstimatorUKF::ComputeSigmaPointsAndWeights() {
  // first, compute the square root of covariance matrix
  // ... that said, the most common implemention is perform Cholesky Decomposition, and take the lower triangular
  // matrix as the square root.
  // 
  // i.e. P = L * L^T, and take L as the square root
  MatrixXf covSqrt = ukfCov.llt().matrixL();

  // sigma points are determined by the following rules:
  // 1. mu_t as the first sigma point
  // 2. mu_t + gamma * S_i as the 1 .. N sigma points, S_i means the ith column of the square root of covariance matrix
  // 3. mu_t - gamma * S_{i-N} as the N + 1 .. 2N sigma points.
  //
  // all parameters (alpha, beta, kappa, lambda, gamma) are initialized in Init()
  ukfSigmaPoints.col(0) = ukfState;
  ukfSigmaPoints.block<QUAD_UKF_NUM_STATES, QUAD_UKF_NUM_STATES>(0, 1) = (gamma * covSqrt).colwise() + ukfState;
  ukfSigmaPoints.block<QUAD_UKF_NUM_STATES, QUAD_UKF_NUM_STATES>(0, QUAD_UKF_NUM_STATES + 1) = (-gamma * covSqrt).colwise() + ukfState;

  // mean weights:
  // lambda / (N+lambda) for i = 0
  // 1 / (2 * (N+lambda)) for i > 0
  float k = QUAD_UKF_NUM_STATES + lambda;
  ukfMeanWeights.setConstant(1 / (2 * k));
  ukfMeanWeights(0) = lambda / k;

  // covariance weights:
  // lamda / (N+lambda) + ( 1 - alpha ^ 2 + beta ) for i = 0
  // 1 / (2 * (N+lambda)) for i > 0
  ukfCovWeights.setConstant(1 / (2 * k));
  ukfCovWeights(0) = lambda / k + (1 - powf(alpha, 2) + beta);
}

VectorXf QuadEstimatorUKF::PredictState(VectorXf curState, float dt, V3F accel, V3F gyro)
{
  assert(curState.size() == QUAD_UKF_NUM_STATES);
  VectorXf predictedState = curState;
  // Predict the current state forward by time dt using current accelerations and body rates as input
  // INPUTS: 
  //   curState: starting state
  //   dt: time step to predict forward by [s]
  //   accel: acceleration of the vehicle, in body frame, *not including gravity* [m/s2]
  //   gyro: body rates of the vehicle, in body frame [rad/s]
  //   
  // OUTPUT:
  //   return the predicted state as a vector

  // HINTS 
  // - dt is the time duration for which you should predict. It will be very short (on the order of 1ms)
  //   so simplistic integration methods are fine here
  // - we've created an Attitude Quaternion for you from the current state. Use 
  //   attitude.Rotate_BtoI(<V3F>) to rotate a vector from body frame to inertial frame
  // - the yaw integral is already done in the IMU update. Be sure not to integrate it again here

  Quaternion<float> attitude = Quaternion<float>::FromEuler123_RPY(rollEst, pitchEst, curState(6));

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  
  VectorXf term2(QUAD_UKF_NUM_STATES);
  V3F v = attitude.Rotate_BtoI(accel);
  term2 << 0, 0, 0, v[0], v[1], v[2], 0;
  term2 *= dt;
  VectorXf term1(predictedState);
  term1[0] += term1[3] * dt;
  term1[1] += term1[4] * dt;
  term1[2] += term1[5] * dt;
  term1[5] -= static_cast<float>(CONST_GRAVITY * dt);
  predictedState = term1 + term2;

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return predictedState;
}

MatrixXf QuadEstimatorUKF::GetRbgPrime(float roll, float pitch, float yaw)
{
  // first, figure out the Rbg_prime
  MatrixXf RbgPrime(3, 3);
  RbgPrime.setZero();

  // Return the partial derivative of the Rbg rotation matrix with respect to yaw. We call this RbgPrime.
  // INPUTS: 
  //   roll, pitch, yaw: Euler angles at which to calculate RbgPrime
  //   
  // OUTPUT:
  //   return the 3x3 matrix representing the partial derivative at the given point

  // HINTS
  // - this is just a matter of putting the right sin() and cos() functions in the right place.
  //   make sure you write clear code and triple-check your math
  // - You can also do some numerical partial derivatives in a unit test scheme to check 
  //   that your calculations are reasonable

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  float cosPhi = cosf(roll);
  float sinPhi = sinf(roll);
  float cosTheta = cosf(pitch);
  float sinTheta = sinf(pitch);
  float cosPsi = cosf(yaw);
  float sinPsi = sinf(yaw);

  RbgPrime << -cosTheta * sinPsi, -sinPhi * sinTheta * sinPsi - cosTheta * cosPsi, -cosPhi * sinTheta * sinPsi + sinPhi * cosPsi,
               cosTheta * cosPsi, sinPhi * sinTheta * cosPsi - cosPhi * sinPsi, cosPhi * sinTheta * cosPsi + sinPhi * sinPsi,
               0, 0, 0;

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return RbgPrime;
}

void QuadEstimatorUKF::Predict(float dt, V3F accel, V3F gyro)
{
  // Compute sigma points with current mean and covariance
  ComputeSigmaPointsAndWeights();
  // Predict next values of sigma points
  for (int i = 0; i < QUAD_UKF_NUM_SIGMA_POINTS; ++i) {
    ukfSigmaPoints.col(i) = PredictState(ukfSigmaPoints.col(i), dt, accel, gyro);
  }
  ukfState = ukfSigmaPoints * ukfMeanWeights;
  //cout << "sigma points:"<<endl<< ukfSigmaPoints << endl;
  //cout << "weights:"<<endl<<ukfMeanWeights << endl;
  cout << "new mean"<<endl<<ukfState.transpose() << endl;
  MatrixXf diff_sp_mu = ukfSigmaPoints.colwise() - ukfState;
  ukfCov = (ukfCovWeights.replicate(1, QUAD_UKF_NUM_STATES).transpose().array() * diff_sp_mu.array()).matrix() * diff_sp_mu.transpose() + Q;
  //cout << "cov" << ukfCov << endl;
  //cout << endl;
}

void QuadEstimatorUKF::UpdateFromGPS(V3F pos, V3F vel)
{
  VectorXf z(6);
  MatrixXf zFromX(6, QUAD_UKF_NUM_SIGMA_POINTS);
  z(0) = pos.x;
  z(1) = pos.y;
  z(2) = pos.z;
  z(3) = vel.x;
  z(4) = vel.y;
  z(5) = vel.z;

  MatrixXf hPrime(6, QUAD_UKF_NUM_STATES);
  hPrime.setZero();

  // GPS UPDATE
  // Hints: 
  //  - The GPS measurement covariance is available in member variable R_GPS
  //  - this is a very simple update
  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  int diagSize = QUAD_UKF_NUM_STATES - 1;
  hPrime.topLeftCorner(diagSize, diagSize) = MatrixXf::Identity(diagSize, diagSize);
  zFromX = hPrime * ukfSigmaPoints;

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  Update(z, R_GPS, zFromX);
}

void QuadEstimatorUKF::UpdateFromMag(float magYaw)
{
  VectorXf z(1);
  MatrixXf zFromX(1, QUAD_UKF_NUM_SIGMA_POINTS);
  z(0) = magYaw;

  MatrixXf hPrime(1, QUAD_UKF_NUM_STATES);
  hPrime.setZero();

  // MAGNETOMETER UPDATE
  // Hints: 
  //  - Your current estimated yaw can be found in the state vector: ukfState(6)
  //  - Make sure to normalize the difference between your measured and estimated yaw
  //    (you don't want to update your yaw the long way around the circle)
  //  - The magnetomer measurement covariance is available in member variable R_Mag
  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  hPrime(0, 6) = 1;
  zFromX = hPrime * ukfSigmaPoints;
  for (int i = 0; i < QUAD_UKF_NUM_SIGMA_POINTS; ++i) {
    if (magYaw - zFromX(0, i) > F_PI) {
      zFromX(0, i) += 2 * F_PI;
    }
    if (magYaw - zFromX(0, i) < -F_PI) {
      zFromX(0, i) -= 2 * F_PI;
    }
  }

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  Update(z, R_Mag, zFromX);
}

// Execute an EKF update step
// z: measurement
// R: observation error model covariance
// zFromSigmaPoints: measurement prediction based on sigma points
void QuadEstimatorUKF::Update(VectorXf& z, MatrixXf& R, MatrixXf& zFromSigmaPoints)
{
  // 1. compute new mean mu_t_bar before measurement
  // 2. compute new covariance sigma_t_bar before measure
  // these two steps are done in Predict.
  // 3. compute z mean
  VectorXf mu_z = zFromSigmaPoints * ukfMeanWeights;
  // 4. compute convariance for measurement
  MatrixXf diff_zsp_mu = zFromSigmaPoints.colwise() - mu_z;
  MatrixXf sigma_z = (ukfCovWeights.replicate(1, zFromSigmaPoints.rows()).transpose().array() * diff_zsp_mu.array()).matrix() * diff_zsp_mu.transpose() + R;
  // 5. compute convariance between prediction and measure
  MatrixXf diff_sp_mu = ukfSigmaPoints.colwise() - ukfState;
  MatrixXf sigma_xz = (ukfCovWeights.replicate(1, QUAD_UKF_NUM_STATES).transpose().array() * diff_sp_mu.array()).matrix() * diff_zsp_mu.transpose();
  // 6. compute kalman gain
  MatrixXf K = sigma_xz * sigma_z.inverse();
  // 7. compute new mean after measurement
  ukfState = ukfState + K * (z - mu_z);
  // 8. compute new covariance after measurement
  ukfCov = ukfCov - K * sigma_z * K.transpose();
}

// Calculate the condition number of the EKF ovariance matrix (useful for numerical diagnostics)
// The condition number provides a measure of how similar the magnitudes of the error metric beliefs 
// about the different states are. If the magnitudes are very far apart, numerical issues will start to come up.
float QuadEstimatorUKF::CovConditionNumber() const
{
  MatrixXf m(7, 7);
  for (int i = 0; i < 7; i++)
  {
    for (int j = 0; j < 7; j++)
    {
      m(i, j) = ukfCov(i, j);
    }
  }

  Eigen::JacobiSVD<MatrixXf> svd(m);
  float cond = svd.singularValues()(0)
    / svd.singularValues()(svd.singularValues().size() - 1);
  return cond;
}

// Access functions for graphing variables
bool QuadEstimatorUKF::GetData(const string& name, float& ret) const
{
  if (name.find_first_of(".") == string::npos) return false;
  string leftPart = LeftOf(name, '.');
  string rightPart = RightOf(name, '.');

  if (ToUpper(leftPart) == ToUpper(_name))
  {
#define GETTER_HELPER(A,B) if (SLR::ToUpper(rightPart) == SLR::ToUpper(A)){ ret=(B); return true; }
    GETTER_HELPER("Est.roll", rollEst);
    GETTER_HELPER("Est.pitch", pitchEst);

    GETTER_HELPER("Est.x", ukfState(0));
    GETTER_HELPER("Est.y", ukfState(1));
    GETTER_HELPER("Est.z", ukfState(2));
    GETTER_HELPER("Est.vx", ukfState(3));
    GETTER_HELPER("Est.vy", ukfState(4));
    GETTER_HELPER("Est.vz", ukfState(5));
    GETTER_HELPER("Est.yaw", ukfState(6));

    GETTER_HELPER("Est.S.x", sqrtf(ukfCov(0, 0)));
    GETTER_HELPER("Est.S.y", sqrtf(ukfCov(1, 1)));
    GETTER_HELPER("Est.S.z", sqrtf(ukfCov(2, 2)));
    GETTER_HELPER("Est.S.vx", sqrtf(ukfCov(3, 3)));
    GETTER_HELPER("Est.S.vy", sqrtf(ukfCov(4, 4)));
    GETTER_HELPER("Est.S.vz", sqrtf(ukfCov(5, 5)));
    GETTER_HELPER("Est.S.yaw", sqrtf(ukfCov(6, 6)));

    // diagnostic variables
    GETTER_HELPER("Est.D.AccelPitch", accelPitch);
    GETTER_HELPER("Est.D.AccelRoll", accelRoll);

    GETTER_HELPER("Est.D.ax_g", accelG[0]);
    GETTER_HELPER("Est.D.ay_g", accelG[1]);
    GETTER_HELPER("Est.D.az_g", accelG[2]);

    GETTER_HELPER("Est.E.x", trueError(0));
    GETTER_HELPER("Est.E.y", trueError(1));
    GETTER_HELPER("Est.E.z", trueError(2));
    GETTER_HELPER("Est.E.vx", trueError(3));
    GETTER_HELPER("Est.E.vy", trueError(4));
    GETTER_HELPER("Est.E.vz", trueError(5));
    GETTER_HELPER("Est.E.yaw", trueError(6));
    GETTER_HELPER("Est.E.pitch", pitchErr);
    GETTER_HELPER("Est.E.roll", rollErr);
    GETTER_HELPER("Est.E.MaxEuler", maxEuler);

    GETTER_HELPER("Est.E.pos", posErrorMag);
    GETTER_HELPER("Est.E.vel", velErrorMag);

    GETTER_HELPER("Est.D.covCond", CovConditionNumber());
#undef GETTER_HELPER
  }
  return false;
};

vector<string> QuadEstimatorUKF::GetFields() const
{
  vector<string> ret = BaseQuadEstimator::GetFields();
  ret.push_back(_name + ".Est.roll");
  ret.push_back(_name + ".Est.pitch");

  ret.push_back(_name + ".Est.x");
  ret.push_back(_name + ".Est.y");
  ret.push_back(_name + ".Est.z");
  ret.push_back(_name + ".Est.vx");
  ret.push_back(_name + ".Est.vy");
  ret.push_back(_name + ".Est.vz");
  ret.push_back(_name + ".Est.yaw");

  ret.push_back(_name + ".Est.S.x");
  ret.push_back(_name + ".Est.S.y");
  ret.push_back(_name + ".Est.S.z");
  ret.push_back(_name + ".Est.S.vx");
  ret.push_back(_name + ".Est.S.vy");
  ret.push_back(_name + ".Est.S.vz");
  ret.push_back(_name + ".Est.S.yaw");

  ret.push_back(_name + ".Est.E.x");
  ret.push_back(_name + ".Est.E.y");
  ret.push_back(_name + ".Est.E.z");
  ret.push_back(_name + ".Est.E.vx");
  ret.push_back(_name + ".Est.E.vy");
  ret.push_back(_name + ".Est.E.vz");
  ret.push_back(_name + ".Est.E.yaw");
  ret.push_back(_name + ".Est.E.pitch");
  ret.push_back(_name + ".Est.E.roll");

  ret.push_back(_name + ".Est.E.pos");
  ret.push_back(_name + ".Est.E.vel");

  ret.push_back(_name + ".Est.E.maxEuler");

  ret.push_back(_name + ".Est.D.covCond");

  // diagnostic variables
  ret.push_back(_name + ".Est.D.AccelPitch");
  ret.push_back(_name + ".Est.D.AccelRoll");
  ret.push_back(_name + ".Est.D.ax_g");
  ret.push_back(_name + ".Est.D.ay_g");
  ret.push_back(_name + ".Est.D.az_g");
  return ret;
};
