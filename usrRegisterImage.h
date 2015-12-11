# include "usrImage.h"

//  The transform that will map the fixed image into the moving image.
typedef itk::AffineTransform< double, 2 > TransformType;
//  An optimizer is required to explore the parameter space of the transform in search of optimal values of the metric.
typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
//  The metric will compare how well the two images match each other. Metric types are usually parameterized by the image types as it can be seen in the following type declaration.
typedef itk::MeanSquaresImageToImageMetric< ImageType, ImageType > MetricType;
//  Finally, the type of the interpolator is declared. The interpolator will evaluate the intensities of the moving image at non-grid positions.
typedef itk::LinearInterpolateImageFunction< ImageType, double > InterpolatorType;
//  The registration method type is instantiated using the types of the fixed and moving images. This class is responsible for interconnecting all the components that we have described so far.
typedef itk::ImageRegistrationMethod< ImageType, ImageType >    RegistrationType;

typedef RegistrationType::ParametersType ParametersType;
typedef itk::ResampleImageFilter< ImageType, ImageType > ResampleFilterType;

class RegisterImage
{

    public:

    RegisterImage();
    ~RegisterImage();

    void RegistrationFunction();
    void SetFixedImage(Images*);
    void SetMovingImage(Images*);
    void SetInitialMatrix(double[9]);

    double registrationMatrix[9];
    double initialMatrix[9];


    private:

    int CreateRegisteredImage();

    Images* fixedImage;
    Images* movingImage;
    Images* registeredImage;
    RegistrationType::Pointer registration;
    TransformType::Pointer transform;
    OptimizerType::Pointer optimizer;
    ResampleFilterType::Pointer   resampler;

};
