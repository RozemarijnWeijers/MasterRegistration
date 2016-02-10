#include "usrVolumeRegistration.h"
 #include <time.h>

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

void VolumeRegistration::SetInitialMatrix( TransformMatrix matrix1 )
{

  this->initialMatrix[0] = matrix1.matrix(0,0); this->initialMatrix[1] = matrix1.matrix(0,1); this->initialMatrix[2] = matrix1.matrix(0,2); this->initialMatrix[3] = matrix1.matrix(0,3);
  this->initialMatrix[4] = matrix1.matrix(1,0); this->initialMatrix[5] = matrix1.matrix(1,1); this->initialMatrix[6] = matrix1.matrix(1,2); this->initialMatrix[7] = matrix1.matrix(1,3);
  this->initialMatrix[8] = matrix1.matrix(2,0); this->initialMatrix[9] = matrix1.matrix(2,1); this->initialMatrix[10] = matrix1.matrix(2,2); this->initialMatrix[11] = matrix1.matrix(2,3);
  this->initialMatrix[12] = matrix1.matrix(3,0); this->initialMatrix[13] = matrix1.matrix(3,1); this->initialMatrix[14] = matrix1.matrix(3,2); this->initialMatrix[15] = matrix1.matrix(3,3);

  return;

}


class CommandIteration : public itk::Command
{

public:
  typedef CommandIteration   Self;
  typedef itk::Command       Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
protected:
  CommandIteration() {};
public:
  typedef  itk::VersorRigid3DTransformOptimizer   OptimizerType;
  typedef const OptimizerType*                    OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >( object );
    if( ! itk::IterationEvent().CheckEvent( &event ) )
    {
      return;
    }
    std::cout << optimizer->GetCurrentIteration() << " : ";
    std::cout << optimizer->GetValue() << " : ";
    std::cout << optimizer->GetCurrentPosition() << std::endl;
  }

};

void VolumeRegistration::RegisterVolumes()
{

  // Set the registration inputs
  this->registration->SetFixedImage( this->fixedVolume->volumeData );
  this->registration->SetMovingImage( this->movingVolume->volumeData );
  this->registration->SetFixedImageRegion( this->fixedVolume->volumeData->GetLargestPossibleRegion() );

  //  Initialize the transform
  /*typedef itk::CenteredTransformInitializer< TransformType3D, VolumeType, VolumeType >  TransformInitializerType;
  TransformInitializerType::Pointer initializer = TransformInitializerType::New();
  initializer->SetTransform(   transform );
  initializer->SetFixedImage(  this->fixedVolume->volumeData );
  initializer->SetMovingImage( this->movingVolume->volumeData );
  initializer->MomentsOn();
  initializer->InitializeTransform();*/

  typedef TransformType3D::VersorType  VersorType;
  typedef VersorType::VectorType     VectorType;
  typedef TransformType3D::TranslationType TranslationType;
  VersorType      rotation;
  VectorType      axis;
  TranslationType translation;
  axis[0] = 0.0;
  axis[1] = 0.0;
  axis[2] = 1.0;
  const double angle = 0;
  rotation.Set(  axis, angle  );
  this->transform->SetRotation( rotation );
  translation[0] = 0;
  translation[1] = 0;
  translation[2] = 0;
  this->transform->SetTranslation(translation);

  this->registration->SetInitialTransformParameters( this->transform->GetParameters() );

  typedef OptimizerType3D::ScalesType OptimizerScalesType;
  OptimizerScalesType optimizerScales( this->transform->GetNumberOfParameters() );
  const double translationScale = 1 / 1000.0;
  optimizerScales[0] = 1.0/1;
  optimizerScales[1] = 1.0/1;
  optimizerScales[2] = 1.0/1;
  optimizerScales[3] = translationScale/1000;
  optimizerScales[4] = translationScale/1000;
  optimizerScales[5] = translationScale/1000;
  optimizer->SetScales( optimizerScales );
  this->optimizer->SetMaximumStepLength( 0.50  );//in mm?
  this->optimizer->SetMinimumStepLength( 0.005 ); //in mm?

  this->metric->SetNumberOfSpatialSamples(50000);

  // Set a stopping criterion
  this->optimizer->SetNumberOfIterations( 100 );

  CommandIteration::Pointer observer = CommandIteration::New();
  this->optimizer->AddObserver( itk::IterationEvent(), observer );


  try
  {
    clock_t time = clock();

    this->registration->Update();
    clock_t time2 = clock();

    printf("\nTime taken: %.4fs\n", (float)(time2 - time)/double(CLOCKS_PER_SEC));
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return;// EXIT_FAILURE;
  }

  // Create registered version of moving Volume

  this->CreateRegisteredVolume();

  //  The result of the registration process is an array of parameters that defines the spatial transformation in an unique way. This final result is obtained using the \code{GetLastTransformParameters()} method.
  OptimizerType3D::ParametersType finalParameters = registration->GetLastTransformParameters();
  const double versorX              = finalParameters[0];
  const double versorY              = finalParameters[1];
  const double versorZ              = finalParameters[2];
  const double finalTranslationX    = finalParameters[3];
  const double finalTranslationY    = finalParameters[4];
  const double finalTranslationZ    = finalParameters[5];
  const unsigned int numberOfIterations = this->optimizer->GetCurrentIteration();
  const double bestValue = this->optimizer->GetValue();

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

  //std::cerr<< this->registration->GetTransform()[0]<< std::endl;

  return;

}

TransformMatrix VolumeRegistration::GetRegistrationMatrix()
{

  TransformMatrix registrationMatrix;
  registrationMatrix.SetTransformFromDouble( this->registrationMatrix );

  return registrationMatrix;

}

void VolumeRegistration::CreateRegisteredVolume()
{

  // Set moving Volume as input


  TransformType3D::Pointer movetransform = TransformType3D::New();
  typedef TransformType3D::VersorType VersorType;
  typedef VersorType::VectorType VectorType;
  typedef TransformType3D::TranslationType TranslationType;
  VersorType      rotation;
  VectorType      axis;
  TranslationType translation;

  translation[0] = 5;
  translation[1] = 0;
  translation[2] = 0;
  movetransform->SetTranslation(translation);

  // The Transform produced by the Registration method is passed into the resampling filter
  this->resampler->SetTransform( this->registration->GetOutput()->Get());
  this->resampler->SetInput( this->movingVolume->volumeData );

  // Specifying parameters of the output Volume (default pixel value is set to gray in order to highlight the regions that are mapped outside of the moving Volume)
  this->resampler->SetOutputOrigin( this->movingVolume->volumeData->GetOrigin() );
  this->resampler->SetOutputSpacing( this->movingVolume->volumeData->GetSpacing() );
  this->resampler->SetSize( this->movingVolume->volumeData->GetLargestPossibleRegion().GetSize() );
  this->resampler->SetOutputDirection( this->movingVolume->volumeData->GetDirection() );
  //this->resampler->SetDefaultPixelValue( 0 );
  this->resampler->Update();

  // Create registered ITKVolume'
  this->registeredVolume.volumeData =resampler->GetOutput();
  //this->registeredVolume.volumeData->SetOrigin(this->movingVolume->volumeData->GetOrigin());
  //this->registeredVolume.volumeData->SetSpacing(this->movingVolume->volumeData->GetSpacing());
  this->registeredVolume.SetParametersFromITK();

  return;

}

/*transform->SetParameters( finalParameters );
  TransformType::MatrixType matrix = transform->GetMatrix();
  TransformType::OffsetType offset = transform->GetOffset();
  std::cout << "Matrix = " << std::endl << matrix << std::endl;
  std::cout << "Offset = " << std::endl << offset << std::endl;*/

