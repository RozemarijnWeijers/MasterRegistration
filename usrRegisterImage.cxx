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

#include "usrRegisterImage.h"

RegisterImage::RegisterImage()
{

  // Create components registration function
  MetricType::Pointer           metric = MetricType::New();
  InterpolatorType::Pointer     interpolator = InterpolatorType::New();
  transform = TransformType::New();
  optimizer = OptimizerType::New();
  registration = RegistrationType::New();
  //initialMatrix =  // INITIALISEREN!!

  resampler = ResampleFilterType::New();

  // Each component is now connected to the instance of the registration method
  registration->SetMetric( metric );
  registration->SetOptimizer( optimizer );
  registration->SetTransform( transform );
  registration->SetInterpolator( interpolator );

}

RegisterImage::~RegisterImage()
{
}

void RegisterImage::SetFixedImage(Images* imagePointer)
{

  fixedImage = imagePointer;

  return;

}

void RegisterImage::SetMovingImage(Images* imagePointer)
{

  movingImage = imagePointer;

  return;

}

void RegisterImage::SetInitialMatrix(double matrix[9])
{

  registrationMatrix[0] = matrix[0]; registrationMatrix[1] = matrix[1]; registrationMatrix[2] = matrix[2];
  registrationMatrix[3] = matrix[3]; registrationMatrix[4] = matrix[4]; registrationMatrix[5] = matrix[5];
  registrationMatrix[6] = matrix[6]; registrationMatrix[7] = matrix[7]; registrationMatrix[8] = matrix[8];

  return;

}

void RegisterImage::RegistrationFunction()
{

  // Set the registration inputs
  registration->SetFixedImage( fixedImage->imageData );
  registration->SetMovingImage( movingImage->imageData );
  registration->SetFixedImageRegion( fixedImage->imageData->GetLargestPossibleRegion() );

  //  Initialize the transform
  ParametersType initialParameters( transform->GetNumberOfParameters() );

  // USE TRACKER DATA TO SET THE INITIAL PARAMETERS!!!!!!!!!!!!! KLOPT DIT?
  // rotation matrix
  initialParameters[0] = initialMatrix[0];//1.0;
  initialParameters[1] = initialMatrix[1];//0.0;
  initialParameters[2] = initialMatrix[3];//0.0;
  initialParameters[3] = initialMatrix[4];//1.0;
  // translation vector
  initialParameters[4] = initialMatrix[2];//0.0;
  initialParameters[5] = initialMatrix[5];//0.0;

  registration->SetInitialTransformParameters( initialParameters );

  optimizer->SetMaximumStepLength( .2 ); // If this is set too high, you will get a "itk::ERROR: MeanSquaresImageToImageMetric(0xa27ce70): Too many samples map outside moving image buffer: 1818 / 10000" error
  optimizer->SetMinimumStepLength( 0.05 );

  // Set a stopping criterion
  optimizer->SetNumberOfIterations( 100 );

  try
  {
    registration->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return;// EXIT_FAILURE;
  }

  // Create registered version of moving image
  if ( 0 == CreateRegisteredImage() )
  {
    std::cerr << "Registering Image failed" << std::endl;
  }

  //  The result of the registration process is an array of parameters that defines the spatial transformation in an unique way. This final result is obtained using the \code{GetLastTransformParameters()} method.
  ParametersType finalParameters = registration->GetLastTransformParameters();
  registrationMatrix[0] = finalParameters[0];
  registrationMatrix[1] = finalParameters[1];
  registrationMatrix[3] = finalParameters[2];
  registrationMatrix[4] = finalParameters[3];
  registrationMatrix[2] = finalParameters[4];
  registrationMatrix[5] = finalParameters[5];  //KLOPT DIT??

  std::cout << "Final parameters 2D/2D-Registration: " << finalParameters << std::endl;

  //  The value of the image metric corresponding to the last set of parameters can be obtained with the \code{GetValue()} method of the optimizer.
  std::cout << "Metric value: " << optimizer->GetValue() << std::endl;

  return;

}

int RegisterImage::CreateRegisteredImage()
{

  // Set moving image as input
  resampler->SetInput( movingImage->imageData );

  // The Transform produced by the Registration method is passed into the resampling filter
  resampler->SetTransform( registration->GetOutput()->Get() );

  // Specifying parameters of the output image (default pixel value is set to gray in order to highlight the regions that are mapped outside of the moving image)
  resampler->SetSize( fixedImage->imageData->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputOrigin( fixedImage->imageData->GetOrigin() );
  resampler->SetOutputSpacing( fixedImage->imageData->GetSpacing() );
  resampler->SetOutputDirection( fixedImage->imageData->GetDirection() );
  resampler->SetDefaultPixelValue( 0 );
  resampler->Update();

  // Create registered ITKimage
  registeredImage->imageData = resampler->GetOutput();
  //NOG AFMAKEN!!


  return 1;

}
