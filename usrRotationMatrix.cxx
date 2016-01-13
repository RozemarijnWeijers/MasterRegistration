#include "usrRotationMatrix.h"

RotationMatrix::RotationMatrix()
{

  this->matrix.eye(4,4);
  this->matrixDouble[0] = 1; this->matrixDouble[1] = 0; this->matrixDouble[2] = 0;
  this->matrixDouble[3] = 0; this->matrixDouble[4] = 1; this->matrixDouble[5] = 0;
  this->matrixDouble[6] = 0; this->matrixDouble[7] = 0; this->matrixDouble[8] = 1;

}

RotationMatrix::~RotationMatrix()
{}

void RotationMatrix::Set3Angels( double angles[3] )
{

  mat matrixZaxis;
  matrixZaxis     <<  cos(angles[0]*PI/180) <<      -sin(angles[0]*PI/180) <<    0 <<                      endr
               <<  sin(angles[0]*PI/180) <<      cos(angles[0]*PI/180) <<     0 <<                      endr
               <<  0 <<                          0 <<                         1 <<                      endr;
  mat matrixYaxis;
  matrixYaxis     <<  cos(angles[1]*PI/180) <<      0 <<                         sin(angles[1]*PI/180) <<  endr
               <<  0 <<                          1 <<                         0 <<                      endr
               <<  -sin(angles[1]*PI/180) <<     0 <<                         cos(angles[1]*PI/180) <<  endr;
  mat matrixXaxis;
  matrixXaxis     <<  1 <<                   0 <<                                0 <<                      endr
               <<  0 <<                   cos(angles[2]*PI/180) <<           -sin(angles[2]*PI/180) <<  endr
               <<  0 <<                   sin(angles[2]*PI/180) <<            cos(angles[2]*PI/180) <<  endr;
  mat matrix;
  this->matrix = matrixZaxis * matrixYaxis * matrixXaxis;

  this->matrixDouble[0] = this->matrix(0,0); this->matrixDouble[1] = this->matrix(1,0); this->matrixDouble[2] = this->matrix(2,0);
  this->matrixDouble[3] = this->matrix(0,1); this->matrixDouble[4] = this->matrix(1,1); this->matrixDouble[5] = this->matrix(2,1);
  this->matrixDouble[6] = this->matrix(0,2); this->matrixDouble[7] = this->matrix(1,2); this->matrixDouble[8] = this->matrix(2,2);

  return;

}

void RotationMatrix::ShowMatrix()
{
  std::cerr<< "matrixDouble:" << std::endl;
  std::cerr << this->matrixDouble[0] << ", "<< this->matrixDouble[1] <<", "<< this->matrixDouble[2] <<std::endl;
  std::cerr << this->matrixDouble[3] << ", "<< this->matrixDouble[4] <<", "<< this->matrixDouble[5] <<std::endl;
  std::cerr << this->matrixDouble[6] << ", "<< this->matrixDouble[7] <<", "<< this->matrixDouble[8] <<std::endl;
  std::cerr << "Matrix:" << std::endl;
  std::cerr << this->matrix <<std::endl;
  return;

}

