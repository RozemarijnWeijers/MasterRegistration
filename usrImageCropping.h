#include "itkImage.h"
#include "itkResampleImageFilter.h"
#include <itkExtractImageFilter.h>
#include "itkChangeInformationImageFilter.h"
#include "usrImage.h"
#include "usrVolume.h"
#include "itkTileImageFilter.h"
#include <QuickView.h>

const unsigned int InputDimension   = 2;
const unsigned int OutputDimension  = 3;
typedef itk::TileImageFilter< ImageType, VolumeType > Filter2DTo3DType;
typedef itk::ExtractImageFilter< ImageType, ImageType > FilterImageType;

class ImageCropping
{

  public:

  ImageCropping();
  ~ImageCropping();

  void SetImage( Image* );
  void SetCropSizeAndStart( int[2], int[2]);
  void CropImage();
  void Convert2DImageTo3DVolume();
  void SetNumberofImages( int );

  Image                        croppedImage;
  Volume                       croppedVolume;

  private:

  ImageType::IndexType         desiredStart;
  ImageType::SizeType          desiredSize ;
  Image*                       inputImage;
  FilterImageType::Pointer     filter;
  int                          number = 1;

  void SetImageMatrix();

};
