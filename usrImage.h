#ifndef IMAGE_H
#define IMAGE_H

#include "itkImage.h"
#include "igtlImageMessage.h"

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
  void  SetParametersFromITK( double, double );

  ImageType::Pointer            imageData;
  igtl::ImageMessage::Pointer   imgMsg;

  float             originImage[3];
  float             spacingImage[3];

  protected:

  int               sizeImage[3];
  igtl::Matrix4x4   matrixImage;   // image origin and orientation matrix

};

#endif //IMAGE_H
