#include "usrImageRegistration.h"

ImageRegistration::ImageRegistration()
{

  // Create components registration function
  this->metric = MetricType::New();
  this->interpolator = InterpolatorType::New();
  this->transform = TransformType::New();
  this->optimizer = OptimizerType::New();
  this->registration = RegistrationType::New();
  //initialMatrix =  // INITIALISEREN!!

  this->resampler = ResampleFilterType::New();

  // Each component is now connected to the instance of the registration method
  this->registration->SetMetric( this->metric );
  this->registration->SetOptimizer( this->optimizer );
  this->registration->SetTransform( this->transform );
  this->registration->SetInterpolator( this->interpolator );

}

ImageRegistration::~ImageRegistration()
{
}

void ImageRegistration::SetFixedImage(Image* imagePointer)
{

  this->fixedImage = imagePointer;
  //this->registeredImage = imagePointer;

  return;

}

void ImageRegistration::SetMovingImage(Image* imagePointer)
{

  this->movingImage = imagePointer;

  return;

}

void ImageRegistration::SetInitialMatrix(double matrix[9])
{

  this->initialMatrix[0] = matrix[0]; this->initialMatrix[1] = matrix[1]; this->initialMatrix[2] = matrix[2];
  this->initialMatrix[3] = matrix[3]; this->initialMatrix[4] = matrix[4]; this->initialMatrix[5] = matrix[5];
  this->initialMatrix[6] = matrix[6]; this->initialMatrix[7] = matrix[7]; this->initialMatrix[8] = matrix[8];

  return;

}

void ImageRegistration::RegisterImages()
{

  // Set the registration inputs
  this->registration->SetFixedImage( this->fixedImage->imageData );
  this->registration->SetMovingImage( this->movingImage->imageData );
  this->registration->SetFixedImageRegion( this->fixedImage->imageData->GetLargestPossibleRegion() );

  //  Initialize the transform
  ParametersType initialParameters( transform->GetNumberOfParameters() );

  // USE TRACKER DATA TO SET THE INITIAL PARAMETERS!!!!!!!!!!!!! KLOPT DIT?
  // rotation matrix
  initialParameters[0] = this->initialMatrix[0];//1.0;
  initialParameters[1] = this->initialMatrix[1];//0.0;
  initialParameters[2] = this->initialMatrix[3];//0.0;
  initialParameters[3] = this->initialMatrix[4];//1.0;
  // translation vector
  initialParameters[4] = this->initialMatrix[2];//0.0;
  initialParameters[5] = this->initialMatrix[5];//0.0;

  this->registration->SetInitialTransformParameters( initialParameters );

  this->optimizer->SetMaximumStepLength( .2 ); // If this is set too high, you will get a "itk::ERROR: MeanSquaresImageToImageMetric(0xa27ce70): Too many samples map outside moving image buffer: 1818 / 10000" error
  this->optimizer->SetMinimumStepLength( 0.05 );

  // Set a stopping criterion
  this->optimizer->SetNumberOfIterations( 100 );

  try
  {
    this->registration->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return;// EXIT_FAILURE;
  }

  // Create registered version of moving image
  /*if ( 0 == this->CreateRegisteredImage() )
  {
    std::cerr << "Registering Image failed" << std::endl;
  }*/

  //  The result of the registration process is an array of parameters that defines the spatial transformation in an unique way. This final result is obtained using the \code{GetLastTransformParameters()} method.
  ParametersType finalParameters = this->registration->GetLastTransformParameters();
  this->registrationMatrix[0] = finalParameters[0];
  this->registrationMatrix[1] = finalParameters[1];
  this->registrationMatrix[3] = finalParameters[2];
  this->registrationMatrix[4] = finalParameters[3];
  this->registrationMatrix[2] = finalParameters[4];
  this->registrationMatrix[5] = finalParameters[5];

  this->metricValue = this->optimizer->GetValue();

  std::cout << "Final parameters 2D/2D-Registration: " << finalParameters << std::endl;

  //  The value of the image metric corresponding to the last set of parameters can be obtained with the \code{GetValue()} method of the optimizer.
  //std::cout << "Metric value: " << this->optimizer->GetValue() << std::endl;

  return;

}

void ImageRegistration::CreateRegisteredImage()
{

  // Set moving image as input
  this->resampler->SetInput( this->movingImage->imageData );

  // The Transform produced by the Registration method is passed into the resampling filter
  this->resampler->SetTransform( this->registration->GetOutput()->Get() );

  // Specifying parameters of the output image (default pixel value is set to gray in order to highlight the regions that are mapped outside of the moving image)
  this->resampler->SetSize( this->fixedImage->imageData->GetLargestPossibleRegion().GetSize() );
  this->resampler->SetOutputOrigin( this->fixedImage->imageData->GetOrigin() );
  this->resampler->SetOutputSpacing( this->fixedImage->imageData->GetSpacing() );
  this->resampler->SetOutputDirection( this->fixedImage->imageData->GetDirection() );
  this->resampler->SetDefaultPixelValue( 0 );
  this->resampler->Update();

  // Create registered ITKimage'
  this->registeredImage.imageData = this->resampler->GetOutput();
  this->registeredImage.imageData->SetOrigin(this->fixedImage->imageData->GetOrigin());
  this->registeredImage.imageData->SetSpacing(this->fixedImage->imageData->GetSpacing());
  this->registeredImage.SetParametersFromITK( this->fixedImage->originImage[2], this->fixedImage->spacingImage[2] );

  return;

}
