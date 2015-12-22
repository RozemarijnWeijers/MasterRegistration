#ifndef TRANSFORMMATRIX_H
#define TRANSFORMMATRIX_H

#include "igtlImageMessage.h"
#include <armadillo>

using namespace arma;

class TransformMatrix
{

  public:

  TransformMatrix();
  ~TransformMatrix();

  void SetTransformFromDouble( double[16] );
  void SetTransformFromIGT( igtl::ImageMessage::Pointer );
  void SetOriginInTransform( double[3] );
  void SetDimensionsForIGTMatrix( int[3] );
  void SetSpacingForIGTMatrix( float[3] );

  igtl::Matrix4x4 IGTMatrix;

  //private:

  mat matrix;
  int dimensions[3];
  float spacing[3];

  void SetIGTTransformFromMat( mat );

};

#endif //TRANSFORMMATRIX_H
