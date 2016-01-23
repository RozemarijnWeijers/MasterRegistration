#include "usrImage.h"
#include "usrTransformMatrix.h"

class MoveImage
{

  public:

  MoveImage();

  void SetImage( Image* );
  void SetMovement( TransformMatrix* );
  void Move();

  protected:

  Image* imageToMove;
  TransformMatrix* movementMatrix;

};
