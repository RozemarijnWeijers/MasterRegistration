#ifndef IMAGE_H
#define IMAGE_H

#include "itkImage.h"
#include "igtlImageMessage.h"
#include "usrTransformMatrix.h"

typedef  unsigned char PixelType;
typedef  itk::Image< PixelType, 2 > ImageType;

class Image
{

  public:

  Image();
  ~Image();

  int   ConvertIGTtoITKImage();
  int   ConvertITKtoIGTImage();
  void  SetParametersFromIGT();
  void  SetParametersFromITK( double, double, TransformMatrix );

  ImageType::Pointer            imageData;
  igtl::ImageMessage::Pointer   imgMsg;
  TransformMatrix   imageMatrix;   // image origin and orientation matrix

  float             originImage[3];
  float             spacingImage[3];

  //private:

  int               sizeImage[3];


};

#endif //IMAGE_H
