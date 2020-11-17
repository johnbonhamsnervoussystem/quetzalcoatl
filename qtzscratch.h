#include <string>
#include <unordered_map>
#include <Eigen/Core>

class QtzScratch {
  private:
    std::unordered_map<std::string, Eigen::MatrixXd*> _doublemat;
  public:
    Eigen::MatrixXd* dblmat(std::string key);
    void new_dblmat(std::string key, int dim);
    bool is_member(std::string key);
    void print(std::string key);
};
