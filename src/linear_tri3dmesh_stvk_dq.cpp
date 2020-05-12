#ifdef SIM_STATIC_LIBRARY
# include<../include/linear_tri3dmesh_stvk_dq.h>
#endif

template<typename DerivedRet, typename DerivedV, typename DerivedQ, typename DefoType, 
         typename NormalType, typename DNormalType, typename DerivedVol, 
         typename DerivedParam, typename ElementMatrixCallback>
void sim::linear_tri3dmesh_stvk_dq(Eigen::VectorXx<DerivedRet> &out, Eigen::DenseBase<DerivedV> &V,  Eigen::Ref<const Eigen::MatrixXi> E,
                                        const Eigen::MatrixBase<DerivedQ> &q, 
                                        const Eigen::MatrixBase<DefoType> &dphidX, 
                                        const Eigen::MatrixBase<NormalType> &N,
                                        const Eigen::MatrixBase<NormalType> &n,
                                        const Eigen::MatrixBase<DNormalType> &dndq,
                                        const Eigen::MatrixBase<DerivedVol>  &volume, 
                                        const Eigen::MatrixBase<DerivedParam> &params,
                                        const ElementMatrixCallback func) {

    auto assemble_func = [&q, &func](auto &H,  auto &e, 
                            const auto &dphidX, const auto &N, const auto &n, const auto &dndq,
                            const auto &volume, const auto &params) 
                           { 
                             linear_tri3d_stvk_dq(H, q, e, sim::unflatten<3,3>(dphidX), N, n, dndq, params, volume(0));
                             func(H); //callback stuff
                           };
    

    Eigen::Vector9x<DerivedRet> Htmp;
    sim::assemble(out, 3*V.rows(), E, E, assemble_func, Htmp, dXinv, N, n, dndq, volume, params);
}
