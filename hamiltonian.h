namespace hamiltonian {

  class molecular {
    private:
      std::vector<double> atomic_number;
      std::vector<double> coordinates;
  
    public:
      molecular(atomic_number, coordinates);
      double nnrep(void);
  
  }
}
