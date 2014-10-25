/**
 * @file mtl_test.cpp
 * Test script for interfacing with MTL4 and it's linear solvers.
 */

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>


class IdentityMatrix {
public:
  
  template <typename VectorIn, typename VectorOut, typename Assign>
  void mult( const VectorIn& v, VectorOut& w, Assign) const {
    auto v_iter = v.begin();
    auto w_iter = w.begin();

    while (v_iter != v.end() ) {
      Assign::apply(*w_iter, *v_iter);
      v_iter++;
      w_iter++;
    }
  }
  
  template<typename Vector>
  mtl::vec::mat_cvec_multiplier<IdentityMatrix, Vector>
  operator*(const Vector& v) const {
    return mtl::vec::mat_cvec_multiplier
                           <IdentityMatrix, Vector>(*this, v);
  }


  IdentityMatrix(std::size_t num_columns, std::size_t num_rows) : matrix_params_({num_columns, num_rows}) {};

  std::size_t num_rows() const {
    return matrix_params_.num_rows;
  }

  std::size_t num_cols() const {
    return matrix_params_.num_cols;
  }

  private:

  struct matrix_params {

    std::size_t num_rows;
    std::size_t num_cols;
  };

  matrix_params matrix_params_;

};


  namespace mtl {
    namespace ashape {


      template<>
      struct ashape_aux<IdentityMatrix> {
        typedef nonscal type;
      };
    }

    template <>
    struct Collection<IdentityMatrix> {
      typedef double value_type;
      typedef unsigned size_type;
    };
  }

  inline std::size_t size(const IdentityMatrix& A) {
    return A.num_rows()*A.num_cols();
  }

  inline std::size_t num_rows(const IdentityMatrix& A) {
    return A.num_rows();
  }

  inline std::size_t num_cols(const IdentityMatrix& A) {
    return A.num_cols();
  }

int main()
{

  std::size_t N = 1000;
  IdentityMatrix I = IdentityMatrix(N,N);
  mtl::dense_vector<double> x(N, 1.0), b(N);
  b = I * x; x = 0;

  itl::noisy_iteration<double> iter(b, 500, 1.e-6);

  itl::cg(I, x, b, iter);

  return 0;
}
