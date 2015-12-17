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
  ParametersType3D initialParameters( transform->GetNumberOfParameters() );

  // USE TRACKER DATA TO SET THE INITIAL PARAMETERS!!!!!!!!!!!!! KLOPT DIT?
  // rotation matrix
  initialParameters[0] = this->initialMatrix[0];//1.0;
  initialParameters[1] = this->initialMatrix[1];//0.0;
  initialParameters[2] = this->initialMatrix[2];//0.0;
  initialParameters[3] = this->initialMatrix[4];//1.0;
  initialParameters[3] = this->initialMatrix[5];//1.0;
  initialParameters[3] = this->initialMatrix[6];//1.0;
  initialParameters[3] = this->initialMatrix[8];//1.0;
  initialParameters[3] = this->initialMatrix[9];//1.0;
  initialParameters[3] = this->initialMatrix[10];//1.0;
  // translation vector
  initialParameters[4] = this->initialMatrix[3];//0.0;
  initialParameters[5] = this->initialMatrix[7];//0.0;
  initialParameters[3] = this->initialMatrix[11];//1.0;

  this->registration->SetInitialTransformParameters( initialParameters );

  this->optimizer->SetMaximumStepLength( .2 ); // If this is set too high, you will get a "itk::ERROR: MeanSquaresVolumeToVolumeMetric(0xa27ce70): Too many samples map outside moving Volume buffer: 1818 / 10000" error
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

  // Create registered version of moving Volume
  /*if ( 0 == this->CreateRegisteredVolume() )
  {
    std::cerr << "Registering Volume failed" << std::endl;
  }*/

  //  The result of the registration process is an array of parameters that defines the spatial transformation in an unique way. This final result is obtained using the \code{GetLastTransformParameters()} method.
  ParametersType3D finalParameters = this->registration->GetLastTransformParameters();
  this->registrationMatrix[0] = finalParameters[0];
  this->registrationMatrix[1] = finalParameters[1];
  this->registrationMatrix[3] = finalParameters[2];
  this->registrationMatrix[4] = finalParameters[3];
  this->registrationMatrix[2] = finalParameters[4];
  this->registrationMatrix[5] = finalParameters[5];

  this->metricValue = this->optimizer->GetValue();

  std::cout << "Final parameters 2D/2D-Registration: " << finalParameters << std::endl;

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
