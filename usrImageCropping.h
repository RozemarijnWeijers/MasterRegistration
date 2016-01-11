#include "itkImage.h"
#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkVersorRigid3DTransformOptimizer.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include <itkExtractImageFilter.h>
#include "itkVersorRigid3DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "usrImage.h"
#include "usrTransformMatrix.h"

typedef itk::ExtractImageFilter< ImageType, ImageType > FilterImageType;

class ImageCropping
{

  public:

  ImageCropping();
  ~ImageCropping();

  void SetImage( Image* );
  void SetCropSizeAndStart( int[2], int[2]);
  void CropImage();

  Image                        croppedImage;

  private:

  ImageType::IndexType         desiredStart;
  ImageType::SizeType          desiredSize ;
  Image*                       inputImage;
  FilterImageType::Pointer          filter;

};
