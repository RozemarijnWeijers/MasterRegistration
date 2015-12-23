#include "usrTransformMatrix.h"

TransformMatrix::TransformMatrix()
{

  //initialize
  this->matrix.eye(4,4);
  this->spacingSetCheck = false;
  this->dimensionsSetCheck = false;

}

TransformMatrix::~TransformMatrix()
{
}

void TransformMatrix::SetTransformFromDouble( double transformMatrix[16] )
{

  this->matrix << transformMatrix[0] << transformMatrix[1] << transformMatrix[2] << transformMatrix[3] << endr
               << transformMatrix[4] << transformMatrix[5] << transformMatrix[6] << transformMatrix[7] << endr
               << transformMatrix[8] << transformMatrix[9] << transformMatrix[10] << transformMatrix[11] << endr
               << transformMatrix[12] << transformMatrix[13] << transformMatrix[14] << transformMatrix[15] << endr;

  this->SetIGTTransformFromMat( this->matrix );

  return;

}

void TransformMatrix::SetTransformFromIGT( igtl::ImageMessage::Pointer imgMsg )
{

  double originITK[3];
  imgMsg->GetMatrix( this->IGTMatrix );
  imgMsg->GetSpacing( this->spacing );
  imgMsg->GetDimensions( this->dimensions );
  for ( int i=0; i<3; i++)
  {
    originITK[i] = this->IGTMatrix[i][3]-((this->dimensions[i]-1)*spacing[i]/2);
  }

  float transformD[16];
  transformD[0] = this->IGTMatrix[0][0]; transformD[1] = this->IGTMatrix[0][1]; transformD[2] = this->IGTMatrix[0][2]; transformD[3] = originITK[0];
  transformD[4] = this->IGTMatrix[1][0]; transformD[5] = this->IGTMatrix[1][1]; transformD[6] = this->IGTMatrix[1][2]; transformD[7] = originITK[1];
  transformD[8] = this->IGTMatrix[2][0] ; transformD[9] = this->IGTMatrix[2][1];transformD[10] = this->IGTMatrix[2][2]; transformD[11] = originITK[2];
  transformD[12] = this->IGTMatrix[3][0]; transformD[13] = this->IGTMatrix[3][1];transformD[14] = this->IGTMatrix[3][2]; transformD[15] = this->IGTMatrix[3][3];

  this->SetTransformFromDouble( (double*)transformD );

  return;

}

void TransformMatrix::SetIGTTransformFromMat( mat matrix )
{

  //dimensions and spacing!!
  if ( this->spacingSetCheck && this->dimensionsSetCheck == true )
  {
      double originIGT[3];
      for ( int i=0; i<3; i++)
      {
        originIGT[i] = matrix(i,3)+((this->dimensions[i]-1)*this->spacing[i]/2);
      }

      this->IGTMatrix[0][0] = matrix(0,0); this->IGTMatrix[0][1] = matrix(0,1); this->IGTMatrix[0][2] = matrix(0,2); this->IGTMatrix[0][3] = originIGT[0];
      this->IGTMatrix[1][0] = matrix(1,0); this->IGTMatrix[1][1] = matrix(1,1); this->IGTMatrix[1][2] = matrix(1,2); this->IGTMatrix[1][3] = originIGT[1];
      this->IGTMatrix[2][0] = matrix(2,0); this->IGTMatrix[2][1] = matrix(2,1); this->IGTMatrix[2][2] = matrix(2,2); this->IGTMatrix[2][3] = originIGT[2];
      this->IGTMatrix[3][0] = matrix(3,0); this->IGTMatrix[3][1] = matrix(3,1); this->IGTMatrix[3][2] = matrix(3,2); this->IGTMatrix[3][3] = matrix(3,3);
  }
  else
  {
    std::cerr << "Spacing ans dimensions not set" << std::endl;
  }
  return;

}

void TransformMatrix::SetOriginInTransform( float origin[3] )
{

  //dimensions!!
  this->matrix(0,3) = origin[0];
  this->matrix(1,3) = origin[1];
  this->matrix(2,3) = origin[2];

  if ( this->spacingSetCheck && this->dimensionsSetCheck == true )
  {
    SetIGTTransformFromMat( this->matrix );
  }
  else
  {
    std::cerr << "Spacing ans dimensions not set" << std::endl;
  }

  return;

}

void TransformMatrix::SetDirectionInTransform( double dir[9] )
{
  this->matrix(0,0) = dir[0];
  if ( this->spacingSetCheck && this->dimensionsSetCheck == true )
  {
    SetIGTTransformFromMat( this->matrix );
  }
  else
  {
    std::cerr << "Spacing ans dimensions not set" << std::endl;
  }

  return;
}

void TransformMatrix::SetDimensionsForIGTMatrix( int dim[3] )
{

  this->dimensions[0] = dim[0];
  this->dimensions[1] = dim[1];
  this->dimensions[2] = dim[2];
  this->dimensionsSetCheck = true;

  return;

}

void TransformMatrix::SetSpacingForIGTMatrix( float spac[3] )
{

  this->spacing[0] = spac[0];
  this->spacing[1] = spac[1];
  this->spacing[2] = spac[2];
  this->spacingSetCheck = true;

  return;

}
