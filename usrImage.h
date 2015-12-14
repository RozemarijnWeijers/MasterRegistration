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
  void  SetParametersFromITK( double[3], double[3] );//ImageType::PointType, ImageType::SpacingType );

  ImageType::Pointer            imageData;
  igtl::ImageMessage::Pointer   imgMsg;

  protected:

  int               sizeImage[3];
  float             originImage[3];
  float             spacingImage[3];
  igtl::Matrix4x4   matrixImage;   // image origin and orientation matrix

};

#endif //IMAGE_H
