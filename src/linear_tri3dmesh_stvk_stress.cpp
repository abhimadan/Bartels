#ifdef SIM_STATIC_LIBRARY
# include <../include/linear_tri3dmesh_stvk_stress.cpp>
#endif

template<typename DerivedRet, typename DerivedV, typename DerivedQ,
         typename DefoType, typename DerivedParam>
void sim::linear_tri3dmesh_stvk_stress(
    Eigen::PlainObjectBase<DerivedRet> &stress, const Eigen::MatrixBase<DerivedV> &V,
    Eigen::Ref<const Eigen::MatrixXi> E, const Eigen::MatrixBase<DerivedQ> &q,
    const Eigen::MatrixBase<DefoType> &dphidX,
    const Eigen::MatrixBase<DerivedParam> &params) {
  stress.resize(E.rows(), 6);

  for(unsigned int ii=0; ii<E.rows(); ++ii) {
      Eigen::Vector9x<typename DerivedQ::Scalar> ue;
      ue << q.segment(3*E(ii,0),3), q.segment(3*E(ii,1),3), q.segment(3*E(ii,2),3);

      Eigen::Matrix<typename DefoType::Scalar, 9,9> B = sim::flatten_multiply_right<Eigen::Matrix<typename DefoType::Scalar, 3,3> >(sim::unflatten<3,3>(dphidX.row(ii)));

      Eigen::Vector9x<typename DerivedRet::Scalar> dF;
      sim::dpsi_stvk_dF(dF, (sim::unflatten<3,3>(ue)*sim::unflatten<3,3>(dphidX.row(ii))).eval(), params.row(ii));

      Eigen::Matrix3x<typename DerivedRet::Scalar> sigma = sim::unflatten<3,3>(dF)*((sim::unflatten<3,3>(ue)*sim::unflatten<3,3>(dphidX.row(ii))).eval()).transpose();
      stress.row(ii) << sigma(0,0), sigma(1,1), sigma(2,2), sigma(1,2), sigma(0,2), sigma(0,1);
  }
}
