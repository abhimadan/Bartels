#ifdef SIM_STATIC_LIBRARY
# include <../include/linear_tri2dmesh_linear_elasticity_stress.h>
#endif

template<typename DerivedRet, typename DerivedV, typename DerivedQ,
         typename DefoType, typename DerivedParam>
void sim::linear_tri2dmesh_linear_elasticity_stress(
    Eigen::PlainObjectBase<DerivedRet> &stress, const Eigen::MatrixBase<DerivedV> &V,
    Eigen::Ref<const Eigen::MatrixXi> E, const Eigen::MatrixBase<DerivedQ> &q,
    const Eigen::MatrixBase<DefoType> &dphidX,
    const Eigen::MatrixBase<DerivedParam> &params) {
  stress.resize(E.rows(), 6);

  Eigen::Matrix3x<double> F = Eigen::Matrix3x<double>::Zero();
  Eigen::Matrix6x<double> dF2; 
  Eigen::Matrix<double, 9,6> P;
  Eigen::Matrix<double, 6,9> P2; //collect up off diagonal shears to reduce from 6 to 9
  Eigen::Vector6x<double> ue;

  P<<1, 0, 0, 0, 0, 0,
     0, 1, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0,
     0, 0, 1, 0, 0, 0,
     0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 1,
     0, 0, 0, 0, 0, 0;

  P2<<1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 1, 0, 1, 0,
      0, 0, 1, 0, 0, 0, 1, 0, 0,
      0, 1, 0, 1, 0, 0, 0, 0, 0;

  for(unsigned int ii=0; ii<E.rows(); ++ii) {

      F.block(0,0,3,2) = sim::unflatten<3,2>(dphidX.row(ii)); 

      Eigen::Matrix<double, 6,6> B = P2*sim::flatten_multiply_right<Eigen::Matrix3d>(F)*P;

      ue << q.segment(2*E(ii,0),2), q.segment(2*E(ii,1),2), q.segment(2*E(ii,2),2);

      //grab per element positions
      sim::d2psi_linear_elasticity_de2(dF2, params.row(ii));

      stress.row(ii) = -(dF2*B*ue).transpose();
  }
}
