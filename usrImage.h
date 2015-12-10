#include "itkImage.h"
#include "itkImageRegistrationMethod.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkAffineTransform.h"
#include "itkSubtractImageFilter.h"
#include "itkTranslationTransform.h"
#include <itkMatrix.h>
#include <itkNrrdImageIO.h>
#include <itkExtractImageFilter.h>

#include "igtlOSUtil.h"
#include "igtlMessageHeader.h"
#include "igtlTransformMessage.h"
#include "igtlPositionMessage.h"
#include "igtlImageMessage.h"
#include "igtlClientSocket.h"
#include "igtlStatusMessage.h"
#include <igtl_util.h>

typedef  unsigned char PixelType;
typedef  itk::Image< PixelType, 2 > ImageType;

class Images
{

  public:

  Images();
  int   IGTtoITKImage();
  int   ITKtoIGTImage();
  void  SetParametersFromIGT();
  void  SetParametersFromITK( double[3], double[3] );//ImageType::PointType, ImageType::SpacingType );

  ImageType::Pointer            imageData;
  igtl::ImageMessage::Pointer   imgMsg;

  protected:

  int               sizeIm[3];
  float             originIm[3];
  float             spacingIm[3];
  igtl::Matrix4x4   matrixIm;   // image origin and orientation matrix

};
