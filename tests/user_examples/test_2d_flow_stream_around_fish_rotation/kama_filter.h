#ifndef KAMA_FILTER_H
#define KAMA_FILTER_H

#include "base_data_type.h"

namespace SPH
{
    class KalmanFilter {
    public:
        KalmanFilter(Real dt, Real covariance) : dt_(dt), covariance_(covariance)
        {
            //State matrix
            P_ = Mat2d::Identity();
            H_ = Mat2d::Identity();
            R_ = Mat2d::Identity();
            x_ = Vec2d::Zero();
        }

        void filter(const Real z)
        {
            //Transition matrix corresponding to the system state
            A_ = Mat2d{ {1.0, dt_},{0.0, 1.0} };
            //Dynamics of the deterministic disturbance and projection onto the system state
            B_ = Vec2d(0.5 * dt_ * dt_, dt_);
            //Process noise
            Q_ = Mat2d{ {covariance_, 0.0},{0.0, covariance_} };

            // predict state
            Vec2d x_hat_minus = A_ * x_;// +B_ * aw;
            Mat2d P_minus = A_ * P_ * A_.transpose() + Q_;

            // calculate kalmann coefficient
            Mat2d K = P_minus * H_.transpose() * (H_ * P_minus * H_.transpose() + R_).inverse();

            // update the state and covariance
            x_ = x_hat_minus + K * Vec2d(z - (H_ * x_hat_minus)(0), 0);
            P_ = (Mat2d::Identity() - K * H_) * P_minus * (Mat2d::Identity() - K * H_).transpose() + K * R_ * K.transpose();
        }

        Vec2d get_state() const
        {
            return x_;
        }

    private:
        Real dt_, aw, covariance_;
        Vec2d x_;
        Mat2d P_, Q_, R_, A_, H_;
        Vec2d B_;
    };
}
#endif // KAMA_FILTER_H