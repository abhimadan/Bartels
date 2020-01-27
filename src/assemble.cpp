#ifdef SIM_STATIC_LIBRARY
# include<../include/assemble.h>
#endif

//utility objects
struct relay_params_struct {
  template <typename Func, typename ...Params>
  inline relay_params_struct(Func &func, Params && ... params) {
    func(params ...);
  }
};

template<typename Func, typename Ret, typename ...Params>
inline void relay_params(unsigned int eid, Func &func, Ret && tmp, Params && ... params) {
  std::cout<<"Relay: "<<tmp<<"\n";
  relay_params_struct{ func, tmp, (params.row(eid))...};
}

//assembler code
template<typename Func, typename ...Params, typename DerivedRet,  typename DerivedTmp>
void sim::assemble(
                Eigen::SparseMatrix<DerivedRet> &assembled, 
                unsigned int rows, unsigned int cols, 
                Eigen::Ref<Eigen::MatrixXi> E_from,  
                Eigen::Ref<Eigen::MatrixXi> E_to,  
                Func func, Eigen::MatrixBase<DerivedTmp> &tmp,
                Params && ... params) {

    using Scalar = typename DerivedTmp::Scalar;

    //init triplet list
    std::vector<Eigen::Triplet<Scalar>> triplet_list;
    unsigned int block_size_to, block_size_from;

    assembled.resize(rows, cols); 

    std::cout<<"Rows: "<<assembled.rows()<<" Cols: "<<assembled.cols()<<"\n";

    std::cout<<"A "<<tmp<<" \n";
    //iterate through element hypergraph
    for(unsigned int ie=0; ie < E_from.rows(); ++ie) {
        
        std::cout<<"B "<<tmp<<" \n";
        //evaluate func for each hyper edge
        //func(tmp, E_from.row(ie), params...); //run function on element hyper edge 
        relay_params(ie, func,tmp, E_from, params...);
        
        std::cout<<"B1 \n";
        block_size_to = tmp.rows()/E_to.cols();
        block_size_from = tmp.cols()/E_from.cols(); 
        std::cout<<"B2 \n"<<block_size_to<<"\n"<<block_size_from<<"\n";

        for(unsigned int iblock_to=0; iblock_to <E_to.cols(); ++iblock_to) {
          for(unsigned int iblock_from=0; iblock_from<E_from.cols(); ++iblock_from) {
            
            for(unsigned int iblk_r=0; iblk_r < block_size_to; ++iblk_r) {
              for(unsigned int iblk_c=0; iblk_c < block_size_from; ++iblk_c) {
                 triplet_list.push_back(Eigen::Triplet<Scalar>(block_size_to*E_to(ie, iblock_to) + iblk_r, 
                                                       block_size_from*E_from(ie, iblock_from) + iblk_c,
                                                       tmp(iblk_r, iblk_c)));

                 std::cout<<"("<<block_size_to*E_to(ie, iblock_to) + iblk_r<<","<<block_size_from*E_from(ie, iblock_from) + iblk_c<<") \n";
              }
            } 
          }
        }
          std::cout<<"C \n";
    }

    assembled.setFromTriplets(triplet_list.begin(), triplet_list.end());

    std::cout<<"D \n";

 }

template<typename Func, typename ...Params, typename DerivedRet,  typename DerivedTmpA, typename DerivedTmpB>
void sim::assemble(
                Eigen::SparseMatrixBase<DerivedRet> &assembled, 
                Eigen::Ref<Eigen::MatrixXi> E,  
                Func func, Params && ... params, 
                Flatten_Multiply<DerivedTmpA, DerivedTmpB> &tmp) {

}