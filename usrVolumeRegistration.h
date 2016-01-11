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

typedef itk::ResampleImageFilter< VolumeType, VolumeType >                    ResampleFilterType3D;
typedef itk::VersorRigid3DTransform< double >                                 TransformType3D;
typedef itk::VersorRigid3DTransformOptimizer                                  OptimizerType3D;
typedef itk::MeanSquaresImageToImageMetric< VolumeType, VolumeType >          MetricType3D;
typedef itk::LinearInterpolateImageFunction< VolumeType, double >            InterpolatorType3D;
typedef itk::ImageRegistrationMethod< VolumeType, VolumeType >                RegistrationType3D;
typedef RegistrationType3D::ParametersType                                    ParametersType3D;

class VolumeRegistration
{

    public:

    VolumeRegistration();
    ~VolumeRegistration();

    void RegisterVolumes();
    void SetFixedVolume( Volume* );
    void SetMovingVolume( Volume* );
    void SetInitialMatrix( TransformMatrix );
    void CreateRegisteredVolume();
    TransformMatrix GetRegistrationMatrix();

    double registrationMatrix[16];
    double initialMatrix[16];
    double metricValue;
    Volume registeredVolume;

    private:

    Volume*  fixedVolume;
    Volume*  movingVolume;
    RegistrationType3D::Pointer   registration;
    TransformType3D::Pointer      transform;
    OptimizerType3D::Pointer      optimizer;
    MetricType3D::Pointer         metric;
    InterpolatorType3D::Pointer   interpolator;
    ResampleFilterType3D::Pointer resampler;

};
