#include "usrVolumeRegistration.h"

VolumeRegistration::VolumeRegistration()
{

  // Create components registration function
  this->metric = MetricType3D::New();
  this->interpolator = InterpolatorType3D::New();
  this->transform = TransformType3D::New();
  this->optimizer = OptimizerType3D::New();
  this->registration = RegistrationType3D::New();
  //initialMatrix =  // INITIALISEREN!!

  this->resampler = ResampleFilterType3D::New();

  // Each component is now connected to the instance of the registration method
  this->registration->SetMetric( this->metric );
  this->registration->SetOptimizer( this->optimizer );
  this->registration->SetTransform( this->transform );
  this->registration->SetInterpolator( this->interpolator );

}

VolumeRegistration::~VolumeRegistration()
{
}

void VolumeRegistration::SetFixedVolume(Volume* VolumePointer)
{

  this->fixedVolume = VolumePointer;

  return;

}

void VolumeRegistration::SetMovingVolume(Volume* VolumePointer)
{

  this->movingVolume = VolumePointer;

  return;

}

void VolumeRegistration::SetInitialMatrix(double matrix[16])
{

  this->initialMatrix[0] = matrix[0]; this->initialMatrix[1] = matrix[1]; this->initialMatrix[2] = matrix[2]; this->initialMatrix[3] = matrix[3];
  this->initialMatrix[4] = matrix[4]; this->initialMatrix[5] = matrix[5]; this->initialMatrix[6] = matrix[6]; this->initialMatrix[7] = matrix[7];
  this->initialMatrix[8] = matrix[8]; this->initialMatrix[9] = matrix[9]; this->initialMatrix[10] = matrix[10]; this->initialMatrix[11] = matrix[11];
  this->initialMatrix[12] = matrix[12]; this->initialMatrix[13] = matrix[13]; this->initialMatrix[14] = matrix[14]; this->initialMatrix[15] = matrix[15];

  return;

}

void VolumeRegistration::RegisterVolumes()
{

  // Set the registration inputs
  this->registration->SetFixedImage( this->fixedVolume->volumeData );
  this->registration->SetMovingImage( this->movingVolume->volumeData );
  this->registration->SetFixedImageRegion( this->fixedVolume->volumeData->GetLargestPossibleRegion() );

  //  Initialize the transform
  typedef itk::CenteredTransformInitializer< TransformType3D, VolumeType, VolumeType >  TransformInitializerType;
  TransformInitializerType::Pointer initializer = TransformInitializerType::New();
  initializer->SetTransform(   transform );
  initializer->SetFixedImage(  this->fixedVolume->volumeData );
  initializer->SetMovingImage( this->movingVolume->volumeData );
  initializer->MomentsOn();
  initializer->InitializeTransform();

  typedef TransformType3D::VersorType  VersorType;
  typedef VersorType::VectorType     VectorType;
  VersorType     rotation;
  VectorType     axis;
  axis[0] = 0.0;
  axis[1] = 0.0;
  axis[2] = 1.0;
  const double angle = 0;
  rotation.Set(  axis, angle  );
  transform->SetRotation( rotation );

  this->registration->SetInitialTransformParameters( transform->GetParameters() );

  typedef OptimizerType3D::ScalesType       OptimizerScalesType;
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );
  const double translationScale = 1.0 / 1000.0;
  optimizerScales[0] = 1.0;
  optimizerScales[1] = 1.0;
  optimizerScales[2] = 1.0;
  optimizerScales[3] = translationScale;
  optimizerScales[4] = translationScale;
  optimizerScales[5] = translationScale;
  optimizer->SetScales( optimizerScales );
  optimizer->SetMaximumStepLength( 0.2000  );
  optimizer->SetMinimumStepLength( 0.0001 );
  optimizer->SetNumberOfIterations( 200 );

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

  // Create registered version of moving Volume
  /*if ( 0 == this->CreateRegisteredVolume() )
  {
    std::cerr << "Registering Volume failed" << std::endl;
  }*/

  //  The result of the registration process is an array of parameters that defines the spatial transformation in an unique way. This final result is obtained using the \code{GetLastTransformParameters()} method.
  OptimizerType3D::ParametersType finalParameters = registration->GetLastTransformParameters();
  const double versorX              = finalParameters[0];
  const double versorY              = finalParameters[1];
  const double versorZ              = finalParameters[2];
  const double finalTranslationX    = finalParameters[3];
  const double finalTranslationY    = finalParameters[4];
  const double finalTranslationZ    = finalParameters[5];
  const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
  const double bestValue = optimizer->GetValue();

  std::cout << std::endl << std::endl;
  std::cout << "Result = " << std::endl;
  std::cout << " versor X      = " << versorX  << std::endl;
  std::cout << " versor Y      = " << versorY  << std::endl;
  std::cout << " versor Z      = " << versorZ  << std::endl;
  std::cout << " Translation X = " << finalTranslationX  << std::endl;
  std::cout << " Translation Y = " << finalTranslationY  << std::endl;
  std::cout << " Translation Z = " << finalTranslationZ  << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Metric value  = " << bestValue          << std::endl;


  //his->metricValue = this->optimizer->GetValue();

  //std::cout << "Final parameters 2D/2D-Registration: " << finalParameters << std::endl;

  //  The value of the Volume metric corresponding to the last set of parameters can be obtained with the \code{GetValue()} method of the optimizer.
  //std::cout << "Metric value: " << this->optimizer->GetValue() << std::endl;

  return;

}

void VolumeRegistration::CreateRegisteredVolume()
{

  // Set moving Volume as input
  this->resampler->SetInput( this->movingVolume->volumeData );

  // The Transform produced by the Registration method is passed into the resampling filter
  this->resampler->SetTransform( this->registration->GetOutput()->Get() );

  // Specifying parameters of the output Volume (default pixel value is set to gray in order to highlight the regions that are mapped outside of the moving Volume)
  this->resampler->SetSize( this->fixedVolume->volumeData->GetLargestPossibleRegion().GetSize() );
  this->resampler->SetOutputOrigin( this->fixedVolume->volumeData->GetOrigin() );
  this->resampler->SetOutputSpacing( this->fixedVolume->volumeData->GetSpacing() );
  this->resampler->SetOutputDirection( this->fixedVolume->volumeData->GetDirection() );
  this->resampler->SetDefaultPixelValue( 0 );
  this->resampler->Update();

  // Create registered ITKVolume'
  this->registeredVolume.volumeData = this->resampler->GetOutput();
  this->registeredVolume.volumeData->SetOrigin(this->fixedVolume->volumeData->GetOrigin());
  this->registeredVolume.volumeData->SetSpacing(this->fixedVolume->volumeData->GetSpacing());
  this->registeredVolume.SetParametersFromITK();

  return;

}
