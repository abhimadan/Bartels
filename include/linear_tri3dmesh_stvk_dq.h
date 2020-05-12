#ifndef SIM_linear_tri3dmeshstvk_dq_H
#define SIM_linear_tri3dmeshstvk_dq_H

#include <Eigen/Dense>
#include <EigenTypes.h>

#include <assemble.h>
#include <linear_tri3d_stvk_dq.h>

namespace sim {

template<typename DerivedRet,
        typename DerivedV, typename DerivedQ, typename DefoType, 
         typename NormalType,
         typename DNormalType,
         typename DerivedVol, 
         typename DerivedParam, class ElementMatrixCallback = decltype(default_callback)>
void linear_tri3dmesh_stvk_dq(Eigen::VectorXx<DerivedRet> &H, Eigen::DenseBase<DerivedV> &V,  Eigen::Ref<const Eigen::MatrixXi> E,
                                        const Eigen::MatrixBase<DerivedQ> &q,
                                        const Eigen::MatrixBase<DefoType> &dphidX, 
                                        const Eigen::MatrixBase<NormalType> &N,
                                        const Eigen::MatrixBase<NormalType> &n, 
                                        const Eigen::MatrixBase<DNormalType> &dndq,
                                        const Eigen::MatrixBase<DerivedVol>  &volume, 
                                        const Eigen::MatrixBase<DerivedParam> &params,
                                        const ElementMatrixCallback func = default_callback);



}

#ifndef SIM_STATIC_LIBRARY
#   include <../src/linear_tri3dmesh_stvk_dq.cpp>
#endif

#endif
