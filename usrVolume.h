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

#include <vtkNrrdReader.h>
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkMatrix4x4.h"
#include "vtkImageReslice.h"
#include "vtkLookupTable.h"
#include "vtkImageMapToColors.h"

typedef  unsigned char   PixelType;
typedef  itk::Image< PixelType, 3 >  VolumeType;

class Volumes
{

  public:

  Volumes();
  int ITKtoIGTVolume();
  int VTKtoITKVolume();
  void SetParametersFromITK();
  void SetParametersFromVTK( double[3], double[3], int[3] );
  int LoadVolume( char* );

  VolumeType::Pointer volumeData;
  vtkSmartPointer<vtkNrrdReader> VTKreader;
  igtl::ImageMessage::Pointer imgMsg;

  protected:

  int sizeVol[3];
  float originVol[3];
  float spacingVol[3];

};
