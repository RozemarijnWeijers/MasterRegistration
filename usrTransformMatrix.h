#ifndef TRANSFORMMATRIX_H
#define TRANSFORMMATRIX_H

#include "igtlImageMessage.h"
#include <armadillo>
#include <math.h>

#define PI 3.14159265

using namespace arma;

class TransformMatrix
{

  public:

  TransformMatrix();
  ~TransformMatrix();

  void SetTransformFromDouble( double[16] );
  void SetTransformFromIGT( igtl::ImageMessage::Pointer );
  void SetOriginInTransform( float[3] );
  void SetCentreOriginInTransform( float[3] );
  void SetDirectionInTransform( double[9] );
  void SetDimensionsForIGTMatrix( int[3] );
  void SetSpacingForIGTMatrix( float[3] );
  void SetDirectionFrom3Angles( double[3] );
  void SetIGTTransformFromMat();
  void ShowMatrix();
  void MultiplyWith( mat multiplyMatrix , bool NeedsToSetIGTMatrix = true );

  igtl::Matrix4x4 IGTMatrix;

  mat matrix;
  int dimensions[3];
  float spacing[3];

  private:

  bool spacingSetCheck;
  bool dimensionsSetCheck;

};

#endif //TRANSFORMMATRIX_H
