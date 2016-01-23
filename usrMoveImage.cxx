#include "usrMoveImage.h"

MoveImage::MoveImage()
{

}

void MoveImage::SetImage( Image* image )
{

  imageToMove = image;

  return;

}

void MoveImage::SetMovement( TransformMatrix* transform )
{

  movementMatrix = transform;

  return;

}

void MoveImage::Move()
{



}
