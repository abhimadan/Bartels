#ifndef SIM_linear_tri3dstvk_dq_H
#define SIM_linear_tri3dstvk_dq_H

#include <Eigen/Dense>
#include <EigenTypes.h>

#include <dpsi_stvk_dF.h>

namespace sim {

template<typename HessianType, typename DefoType, typename DerivedV, 
         typename NormalType, typename DNormalType,
         typename Scalar, typename DerivedParams>
void linear_tri3d_stvk_dq(Eigen::DenseBase<HessianType> &out, const Eigen::MatrixBase<DerivedV> &q, const Eigen::Ref<const Eigen::RowVectorXi> element,  
                                    const Eigen::MatrixBase<DefoType> &dXinv,
                                    const Eigen::MatrixBase<NormalType> &N, 
                                    const Eigen::MatrixBase<NormalType> &n,
                                    const Eigen::MatrixBase<DNormalType> &dndq,
                                    const Eigen::MatrixBase<DerivedParams> &params, const Scalar &volume);
}

#ifndef SIM_STATIC_LIBRARY
#   include <../src/linear_tri3d_stvk_dq.cpp>
#endif

#endif
