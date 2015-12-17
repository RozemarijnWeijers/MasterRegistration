#include "itkImage.h"
#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkAffineTransform.h"
#include <itkExtractImageFilter.h>
#include "usrImage.h"
#include "usrVolume.h"

typedef itk::AffineTransform< double, 2 > TransformType;
typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
typedef itk::MeanSquaresImageToImageMetric< ImageType, ImageType > MetricType;
typedef itk::LinearInterpolateImageFunction< ImageType, double > InterpolatorType;
typedef itk::ImageRegistrationMethod< ImageType, ImageType >    RegistrationType;
typedef RegistrationType::ParametersType ParametersType;
typedef itk::ResampleImageFilter< ImageType, ImageType > ResampleFilterType;

class ImageRegistration
{

    public:

    ImageRegistration();
    ~ImageRegistration();

    void RegisterImages();
    void SetFixedImage( Image* );
    void SetMovingImage( Image* );
    void SetInitialMatrix( double[9] );
    void CreateRegisteredImage();

    double registrationMatrix[9];
    double initialMatrix[9];
    double metricValue;
    Image  registeredImage;

    private:

    Image*  fixedImage;
    Image*  movingImage;
    RegistrationType::Pointer   registration;
    TransformType::Pointer      transform;
    OptimizerType::Pointer      optimizer;
    MetricType::Pointer         metric;
    InterpolatorType::Pointer   interpolator;
    ResampleFilterType::Pointer resampler;

};
