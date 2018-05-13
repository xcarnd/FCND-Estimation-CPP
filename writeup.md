# Project: Building an Estimator

## Writeup

### 1. Provide a Writeup / README that includes all the rubric points and how you addressed each one. You can submit your writeup as markdown or pdf.

You are reading it! Below I'll describe how I addressed each rubric point and where my code each point is handled.

## Implement Estimator

### 1. Determine the standard deviation of the measurement noise of both GPS X data and Accelerometer X data.

By running the corresponding scene in the simulator, I got two log files. They contained the measurement for GPS x location and acceleometer x velocity. By reading it with Python's pandas package, I can easily figure out the standard deviation of the measurements: `MeasuredStdDev_GPSPosXY` should be around 0.7, and `MeasuredStdDev_AccelXY` should be around 0.49.

### 2. Implement a better rate gyro attitude integration scheme in the UpdateFromIMU() function.

Codes for a better rate gyro attitude integration can be found in `QuadEstimatorEKF.cpp`, L96 ~ L102. 

The way I integrate rate gyro attitude is as following: first, calculate quaternion for the orientation of the drone from the current euler angle estimate, then use `Quaternion::IntegrateBodyRate` to integrate body rate into the quaternion, and lastly convert quaternion back to euler angle. This is the new euler angle estimates.

### 3. Implement all of the elements of the prediction step for the estimator.

See `QuadEstimatorEKF::Predict`, `QuadEstimatorEKF::PredictState` and `QuadEstimatorEKF::GetRbgPrime` for the implementation, in L156 ~ L289. The main functionalities of each function are as follow:

1. `QuadEstimatorEKF::PredictState` is used to perform EKF state prediction step. By correctly implement the transition model in *Estimation for Quadrotors*, it can perform the prediction.

2. `QuadEstimatorEKF::GetRbgPrime` is used for calculating the first order derivative `RbgPrime` of `Rbg`, which is the rotation matrix converting from body frame to global frame. `RbgPrime` is then used in the calculation of Jacobian matrix, the key for performing covariace prediction in EKF.

3. `QuadEstimatorEKF::Predict` is first calling `QuadEstimatorEKF::PredictState` to make prediction of the current state, then perform another prediction of the corresponding covariance.

### 4. Implement the magnetometer update.

Since EKF update routine is already implemented and given to us, the implementation for magnetometer update is simple. All I have to do the realize the measurement model for magnetometer, then call `QuadEstimatorEKF::Update` to do the state update and covariance update. Codes can be found in L318 ~ L341 in `QuadEstimatorEKF::UpdateFromMag`. 

One thing to be kept in mind: the difference of measurement and estimated yaw must be normalized.

### 5. Implement the GPS update.

Similar to the magnetometer update, I only have to implement the measurement model for GPS. Codes can be found in L291 ~ L316.
