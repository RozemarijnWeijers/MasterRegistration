# include "usrVolumeReslice.h"

VolumeReslice::VolumeReslice()
{

  reslice = vtkSmartPointer<vtkImageReslice>::New();
  resliceAxes = vtkSmartPointer<vtkMatrix4x4>::New();
  vtkImageToImageFilter = VTKImageToImageType::New();

  this->axesSetCheck = false;
  this->volumeSetCheck = false;
  this->originOfResliceSetCheck = false;
  this->resliceDoneCheck = false;

}

VolumeReslice::~VolumeReslice()
{
}

void VolumeReslice::SetResliceAxes( double transformMatrix[9] )
{

  // Set the direction matrix for the reslice function
  this->resliceAxes->DeepCopy( transformMatrix );
  this->axesSetCheck = true;

  return ;

}

void VolumeReslice::SetOriginOfResliceWRTVolume( double origin[3] )
{

  this->resliceOriginWRTVolume[0] = origin[0];
  this->resliceOriginWRTVolume[1] = origin[1];
  this->resliceOriginWRTVolume[2] = origin[2] ;
  this->originOfResliceSetCheck = true;

  return;

}

void VolumeReslice::SetOriginOfReslice() // relative to volume
{

  double    volumeOrigin[3];
  this->volume->VTKReader->Update();
  this->volume->VTKReader->GetOutput()->GetOrigin( volumeOrigin );
  bool RAS = true; //???
  if ( RAS )
  {
    volumeOrigin[0] = -1 * volumeOrigin[0];
    volumeOrigin[1] = -1 * volumeOrigin[1];
  }

  this->resliceOrigin[0] = this->resliceOriginWRTVolume[0] + volumeOrigin[0];
  this->resliceOrigin[1] = this->resliceOriginWRTVolume[1] + volumeOrigin[1];
  this->resliceOrigin[2] = this->resliceOriginWRTVolume[2] + volumeOrigin[2];

  return;

}

void VolumeReslice::SetSpacingOfReslice() // relative to volume
{

  this->volume->VTKReader->Update();
  this->volume->VTKReader->GetOutput()->GetSpacing( resliceSpacing );

  return;

}

void VolumeReslice::SetVolume( Volume* volumeInput )
{

  this->volume = volumeInput;
  this->volumeSetCheck = true;

  return;

}

void VolumeReslice::ResliceVolume()
{

  if( this->axesSetCheck && this->volumeSetCheck && this->originOfResliceSetCheck == true )
  {
    this->SetOriginOfReslice();
    // Set the point through which to slice
    this->resliceAxes->SetElement( 0, 3, this->resliceOrigin[0] );
    this->resliceAxes->SetElement( 1, 3, this->resliceOrigin[1] );
    this->resliceAxes->SetElement( 2, 3, this->resliceOrigin[2] );

    // Extract a slice in the desired orientation through the desired point

    this->reslice->SetInputConnection( this->volume->VTKReader->GetOutputPort() );
    this->reslice->SetOutputDimensionality( 2 );
    this->reslice->SetResliceAxes( this->resliceAxes );
    this->reslice->SetInterpolationModeToLinear();
    this->reslice->Update();
    this->resliceDoneCheck = true;

    return;
  }
  else
  {
    std::cerr << "Not all parameters are set" << std::endl;
  }

  return;

  }

void VolumeReslice::CreateITKReslice()
{

  if( this->resliceDoneCheck == true )
  {
    // Create a greyscale lookup table
    /*vtkSmartPointer<vtkLookupTable> table = vtkSmartPointer<vtkLookupTable>::New();
    table->SetRange(0, 2000); // image intensity range
    table->SetValueRange(0.0, 1.0); // from black to white
    table->SetSaturationRange(0.0, 0.0); // no color saturation
    table->SetRampToLinear();
    table->Build();

    // Map the image through the lookup table
    vtkSmartPointer<vtkImageMapToColors> color = vtkSmartPointer<vtkImageMapToColors>::New();
    color->SetLookupTable(table);
    color->SetInputConnection(reslice->GetOutputPort());*/

    vtkSmartPointer<vtkImageActor> actor = vtkSmartPointer<vtkImageActor>::New();
    actor->GetMapper()->SetInputConnection(this->reslice->GetOutputPort());

    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->AddActor(actor);

    vtkSmartPointer<vtkRenderWindow> window = vtkSmartPointer<vtkRenderWindow>::New();
    window->AddRenderer(renderer);

    vtkSmartPointer<vtkInteractorStyleImage> imageStyle = vtkSmartPointer<vtkInteractorStyleImage>::New();
    //vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    //interactor->SetInteractorStyle(imageStyle);
    //window->SetInteractor(interactor);
    window->Render();
    //interactor->Start();

    // Convert the VTK image to an ITK image for further processing
    this->vtkImageToImageFilter->SetInput( reslice->GetOutput() );
    this->vtkImageToImageFilter->Update();

    //typedef itk::RescaleIntensityImageFilter< ImageType, ImageType >   RescalerType;
    //RescalerType::Pointer rescaler = RescalerType::New();
    //rescaler->SetOutputMinimum(  0  );
    //rescaler->SetOutputMaximum( 255 );
    //rescaler->SetInput(this->vtkImageToImageFilter->GetOutput());

    this->reslicedImage.imageData->Graft( this->vtkImageToImageFilter->GetOutput() );

    QuickView viewer;
    viewer.AddImage( this->reslicedImage.imageData.GetPointer() ); // Need to do this because QuickView can't accept smart pointers
    viewer.Visualize();

    this->SetSpacingOfReslice();
    // Set correct image parameters for the ITK image
    this->reslicedImage.imageData->SetOrigin( this->resliceOrigin );
    this->reslicedImage.imageData->SetSpacing( this->resliceSpacing );
    this->reslicedImage.SetParametersFromITK( this->resliceOrigin[2], this->resliceSpacing[2] ); //Set spacing and origin in imageData
  }
  else
  {
    std::cerr << "Volume hasn't been resliced yet" << std::endl;
  }

  return;

}
