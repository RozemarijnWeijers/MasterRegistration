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

  this->SetIGTTransformFromMat();

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

mat MultiplyMatrixAWithB( mat matrixA, mat matrixB )
{

  mat tempMatrix;
  tempMatrix = matrixB * matrixA;

  return tempMatrix;

}

void TransformMatrix::SetIGTTransformFromMat()
{

  if ( this->spacingSetCheck && this->dimensionsSetCheck == true )
  {
      mat tempMat1;
      mat tempMat2;
      mat tempMat3;
      tempMat1 = eye(4,4);
      tempMat2 = this->matrix;
      double originIGT[3];
      for ( int i=0; i<3; i++)
      {
        originIGT[i] = (this->dimensions[i]-1)*this->spacing[i]/2;//matrix(i,3)+((this->dimensions[i]-1)*this->spacing[i]/2);
      }
      tempMat1(0,3) = originIGT[0];
      tempMat1(1,3) = originIGT[1];
      tempMat1(2,3) = originIGT[2];
      tempMat3 = MultiplyMatrixAWithB( tempMat1, tempMat2 );

      this->IGTMatrix[0][0] = this->matrix(0,0); this->IGTMatrix[0][1] = this->matrix(0,1); this->IGTMatrix[0][2] = this->matrix(0,2); this->IGTMatrix[0][3] = tempMat3(0,3);
      this->IGTMatrix[1][0] = this->matrix(1,0); this->IGTMatrix[1][1] = this->matrix(1,1); this->IGTMatrix[1][2] = this->matrix(1,2); this->IGTMatrix[1][3] = tempMat3(1,3);
      this->IGTMatrix[2][0] = this->matrix(2,0); this->IGTMatrix[2][1] = this->matrix(2,1); this->IGTMatrix[2][2] = this->matrix(2,2); this->IGTMatrix[2][3] = tempMat3(2,3);
      this->IGTMatrix[3][0] = this->matrix(3,0); this->IGTMatrix[3][1] = this->matrix(3,1); this->IGTMatrix[3][2] = this->matrix(3,2); this->IGTMatrix[3][3] = this->matrix(3,3);
  }
  else
  {
    std::cerr << "Spacing ans dimensions not set" << std::endl;
  }

  return;

}

void TransformMatrix::SetOriginInTransform( float origin[3], bool RAS )
{

  this->matrix(0,3) = origin[0];
  this->matrix(1,3) = origin[1];
  this->matrix(2,3) = origin[2];
  //dimensions!!
  if (RAS == false)
  {
      this->matrix(0,3) = -origin[0];//LPS
      this->matrix(1,3) = -origin[1];//LPS
      this->matrix(2,3) = origin[2];
  }

  if ( this->spacingSetCheck && this->dimensionsSetCheck == true )
  {
    SetIGTTransformFromMat();
  }
  else
  {
    std::cerr << "Spacing ans dimensions not set" << std::endl;
  }

  return;

}

void TransformMatrix::SetCentreOriginInTransform( float centreOrigin[3] )
{

  mat tempMat1;
  mat tempMat2;
  mat tempMat3;
  tempMat1 = eye(4,4);
  tempMat2 = this->matrix;
  double originCorner[3];
  for ( int i=0; i<3; i++)
  {
     originCorner[i] = -(this->dimensions[i]-1)*this->spacing[i]/2;//matrix(i,3)+((this->dimensions[i]-1)*this->spacing[i]/2);
  }
  tempMat1(0,3) = originCorner[0];
  tempMat1(1,3) = originCorner[1];
  tempMat1(2,3) = originCorner[2];
  tempMat3 = MultiplyMatrixAWithB(tempMat1, tempMat2 );

  this->matrix(0,3) = tempMat3(0,3);
  this->matrix(1,3) = tempMat3(1,3);
  this->matrix(2,3) = tempMat3(2,3);

  if ( this->spacingSetCheck && this->dimensionsSetCheck == true )
  {
    SetIGTTransformFromMat();
  }
  else
  {
    std::cerr << "Spacing ans dimensions not set" << std::endl;
  }

  return;

}

void TransformMatrix::SetDirectionInTransform( double dir[9], bool RAS )
{
  this->matrix(0,0) = dir[0]; this->matrix(1,0) = dir[1]; this->matrix(2,0) = dir[2];
  this->matrix(0,1) = dir[3]; this->matrix(1,1) = dir[4]; this->matrix(2,1) = dir[5];
  this->matrix(0,2) = dir[6]; this->matrix(1,2) = dir[7]; this->matrix(2,2) = dir[8];

  //LPS
  if (RAS == false)
  {
      mat tempMat4;
      mat tempMat5;
      tempMat4 = eye(4,4); tempMat4(0,0) = -1; tempMat4(1,1) = -1;
      tempMat5 = MultiplyMatrixAWithB(this->matrix, tempMat4);
      this->matrix = tempMat5;
  }

  if ( this->spacingSetCheck && this->dimensionsSetCheck == true )
  {
    SetIGTTransformFromMat();
  }
  else
  {
    std::cerr << "Spacing ans dimensions not set" << std::endl;
  }

  return;
}

void TransformMatrix::SetLPS()
{

  mat tempmat1;
  tempmat1 = eye(4,4);
  tempmat1(0,0) = -1; tempmat1(1,1) = -1;

  this->MultiplyWith( tempmat1 );

  return;

}

void TransformMatrix::SetDirectionFrom3Angles( double angles[3] )
{

  mat rotZaxis;
  rotZaxis     <<  cos(angles[0]*PI/180) <<      -sin(angles[0]*PI/180) <<    0 <<                      endr
               <<  sin(angles[0]*PI/180) <<      cos(angles[0]*PI/180) <<     0 <<                      endr
               <<  0 <<                          0 <<                         1 <<                      endr;
  mat rotYaxis;
  rotYaxis     <<  cos(angles[1]*PI/180) <<      0 <<                         sin(angles[1]*PI/180) <<  endr
               <<  0 <<                          1 <<                         0 <<                      endr
               <<  -sin(angles[1]*PI/180) <<     0 <<                         cos(angles[1]*PI/180) <<  endr;
  mat rotXaxis;
  rotXaxis     <<  1 <<                   0 <<                                0 <<                      endr
               <<  0 <<                   cos(angles[2]*PI/180) <<           -sin(angles[2]*PI/180) <<  endr
               <<  0 <<                   sin(angles[2]*PI/180) <<            cos(angles[2]*PI/180) <<  endr;
  mat rot;
  rot = rotZaxis * rotYaxis * rotXaxis;

  this->matrix(0,0) = rot(0,0); this->matrix(1,0) = rot(1,0); this->matrix(2,0) = rot(2,0);
  this->matrix(0,1) = rot(0,1); this->matrix(1,1) = rot(1,1); this->matrix(2,1) = rot(2,1);
  this->matrix(0,2) = rot(0,2); this->matrix(1,2) = rot(1,2); this->matrix(2,2) = rot(2,2);


  if ( this->spacingSetCheck && this->dimensionsSetCheck == true )
  {
    SetIGTTransformFromMat();
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

void TransformMatrix::MultiplyWith( mat multiplyMatrix, bool NeedsToSetIGTMatrix )
{

  mat tempMatrix;
  tempMatrix = this->matrix * multiplyMatrix;
  this->matrix = tempMatrix;
  if (NeedsToSetIGTMatrix)
  {
    this->SetIGTTransformFromMat();
  }

  return;

}

void TransformMatrix::ShowMatrix()
{
  std::cerr<< "IGTMatrix:" << std::endl;
  std::cerr << this->IGTMatrix[0][0] << ", "<< this->IGTMatrix[0][1] <<", "<< this->IGTMatrix[0][2] <<", "<< this->IGTMatrix[0][3] <<std::endl;
  std::cerr << this->IGTMatrix[1][0] << ", "<< this->IGTMatrix[1][1] <<", "<< this->IGTMatrix[1][2] <<", "<< this->IGTMatrix[1][3] <<std::endl;
  std::cerr << this->IGTMatrix[2][0] << ", "<< this->IGTMatrix[2][1] <<", "<< this->IGTMatrix[2][2] <<", "<< this->IGTMatrix[2][3] <<std::endl;
  std::cerr << this->IGTMatrix[3][0] << ", "<< this->IGTMatrix[3][1] <<", "<< this->IGTMatrix[3][2] <<", "<< this->IGTMatrix[3][3] <<std::endl;  std::cerr << "Matrix:" << std::endl;
  std::cerr << this->matrix <<std::endl;
  return;

}
