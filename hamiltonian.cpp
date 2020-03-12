namespace hamiltonian {

    double molecular::nnrep( common& com) {
      /*
        Get the nuclear-nuclear repulsion
      */
      int i, j, natm = com.natm() ;
      double cx = d0, cy = d0, cz = d0 ;
      double r2 = d0, n_rep = d0 ;
      Eigen::MatrixXd c (natm, 3) ;
      Eigen::VectorXd a (natm) ;
      c = com.getC() ;
      a = com.getA() ;
    
      for (i = 0; i < natm; i++){
        for (j = i+1; j < natm; j++){
          cx = c(j, 0) - c(i, 0) ;
          cy = c(j, 1) - c(i, 1) ;
          cz = c(j, 2) - c(i, 2) ;
          r2 = std::pow(cx, d2) + std::pow(cy, d2) + std::pow(cz, d2);
          n_rep += a(i)*a(j)/std::sqrt(r2);
          }
        }
     
      com.nrep(n_rep) ;
    
      std::cout << " Nuclear repulsion is " << n_rep << std::endl ;
    
      a.resize(0) ;
      c.resize(0, 0) ;
    
      return ;
    
      } ;

}
