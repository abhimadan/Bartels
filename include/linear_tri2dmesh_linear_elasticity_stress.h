#ifndef linear_tri2dmesh_linear_elasticity_stress_H
#define linear_tri2dmesh_linear_elasticity_stress_H

#include <Eigen/Dense>
#include <EigenTypes.h>

#include <d2psi_linear_elasticity_de2.h>

namespace sim {

template<typename DerivedRet, typename DerivedV, typename DerivedQ,
         typename DefoType, typename DerivedParam>
void linear_tri2dmesh_linear_elasticity_stress(
    Eigen::PlainObjectBase<DerivedRet> &stress, const Eigen::MatrixBase<DerivedV> &V,
    Eigen::Ref<const Eigen::MatrixXi> E, const Eigen::MatrixBase<DerivedQ> &q,
    const Eigen::MatrixBase<DefoType> &dphidX,
    const Eigen::MatrixBase<DerivedParam> &params);

}

#ifndef SIM_STATIC_LIBRARY
# include <../src/linear_tri2dmesh_linear_elasticity_stress.cpp>
#endif

#endif
