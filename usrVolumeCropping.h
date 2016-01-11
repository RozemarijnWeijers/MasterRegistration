#include "itkImage.h"
#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkVersorRigid3DTransformOptimizer.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include <itkExtractImageFilter.h>
#include "itkVersorRigid3DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "usrVolume.h"
#include "usrTransformMatrix.h"

typedef itk::ExtractImageFilter< VolumeType, VolumeType > FilterType;

class VolumeCropping
{

  public:

  VolumeCropping();
  ~VolumeCropping();

  void SetVolume( Volume* );
  void SetCropSizeAndStart( int[3], int[3]);
  void CropVolume();

  Volume                        croppedVolume;

  private:

  VolumeType::IndexType         desiredStart;
  VolumeType::SizeType          desiredSize ;
  Volume*                       inputVolume;
  FilterType::Pointer           filter;

};
