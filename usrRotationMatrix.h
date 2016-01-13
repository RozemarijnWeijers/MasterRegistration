#include <armadillo>
#include <math.h>

#define PI 3.14159265

using namespace arma;


class RotationMatrix
{

  public:

  RotationMatrix();
  ~RotationMatrix();
  void Set3Angels( double[3] );
  void ShowMatrix();

  mat matrix;
  double matrixDouble[9];

};
