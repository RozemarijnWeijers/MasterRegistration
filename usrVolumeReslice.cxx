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

void VolumeReslice::SetResliceAxesWRTVolume( double transformMatrix[9] )
{

  // Set the direction matrix for the reslice function
  this->transformDirections[0] = transformMatrix[0]; this->transformDirections[1] = transformMatrix[1]; this->transformDirections[2] = transformMatrix[2];
  this->transformDirections[3] = transformMatrix[3]; this->transformDirections[4] = transformMatrix[4]; this->transformDirections[5] = transformMatrix[5];
  this->transformDirections[6] = transformMatrix[6]; this->transformDirections[7] = transformMatrix[7]; this->transformDirections[8] = transformMatrix[8];

  this->axesSetCheck = true;
  return ;

}

void VolumeReslice::SetOriginOfResliceWRTVolume( float origin[3] )
{

  this->resliceOriginWRTVolume[0] = origin[0];
  this->resliceOriginWRTVolume[1] = origin[1];
  this->resliceOriginWRTVolume[2] = origin[2] ;
  this->originOfResliceSetCheck = true;

    return;

}

void VolumeReslice::SetResliceMatrix() // relative to volume
{

  this->resliceMatrixWRTVolume.SetSpacingForIGTMatrix( this->volume->spacingVolume );
  this->resliceMatrixWRTVolume.SetDimensionsForIGTMatrix( this->volume->sizeVolume ); //??
  this->resliceMatrixWRTVolume.SetDirectionInTransform( this->transformDirections );
  this->resliceMatrixWRTVolume.SetOriginInTransform( this->resliceOriginWRTVolume );
  mat tempmat;
  tempmat = this->volume->volumeMatrix.matrix * this->resliceMatrixWRTVolume.matrix;
  this->resliceMatrix.SetSpacingForIGTMatrix( this->volume->spacingVolume );
  this->resliceMatrix.SetDimensionsForIGTMatrix( this->volume->sizeVolume ); //??
  this->resliceMatrix.matrix = tempmat;
  this->resliceMatrix.SetIGTTransformFromMat();

  double transformMatrix1[12];
  transformMatrix1[0] = this->resliceMatrix.matrix(0 ,0); transformMatrix1[1] = this->resliceMatrix.matrix(0 ,1); transformMatrix1[2] = this->resliceMatrix.matrix(0 ,2); transformMatrix1[3] = this->resliceMatrix.matrix(0 ,3);
  transformMatrix1[4] = this->resliceMatrix.matrix(1 ,0); transformMatrix1[5] = this->resliceMatrix.matrix(1 ,1); transformMatrix1[6] = this->resliceMatrix.matrix(1, 2); transformMatrix1[7] = this->resliceMatrix.matrix(1 ,3);
  transformMatrix1[8] = this->resliceMatrix.matrix(2 ,0); transformMatrix1[9] = this->resliceMatrix.matrix(2 ,1); transformMatrix1[10] = this->resliceMatrix.matrix(2 ,2); transformMatrix1[11] = this->resliceMatrix.matrix(2 ,3);
  this->resliceAxes->DeepCopy( transformMatrix1 );

  this->resliceOrigin[0] = this->resliceMatrix.matrix(0,3);// this->resliceMatrix.IGTMatrix[0][3];//
  this->resliceOrigin[1] = this->resliceMatrix.matrix(1,3);// this->resliceMatrix.IGTMatrix[1][3];//
  this->resliceOrigin[2] = this->resliceMatrix.matrix(2,3);// this->resliceMatrix.IGTMatrix[2][3];//
  std::cerr<< this->resliceOrigin[0] <<std::endl;
  /*double    volumeOrigin[3];
  this->volume->VTKReader->Update();
  this->volume->VTKReader->GetOutput()->GetOrigin( volumeOrigin );
  bool RAS = false; //???
  if ( RAS )
  {
    volumeOrigin[0] = -1 * volumeOrigin[0];
    volumeOrigin[1] = -1 * volumeOrigin[1];
  }

  this->resliceOrigin[0] = this->resliceOriginWRTVolume[0] + volumeOrigin[0];
  this->resliceOrigin[1] = this->resliceOriginWRTVolume[1] + volumeOrigin[1];
  this->resliceOrigin[2] = this->resliceOriginWRTVolume[2] + volumeOrigin[2];*/

  return;

}

void VolumeReslice::SetSpacingOfReslice() // relative to volume
{

  this->volume->VTKReader->Update();

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
    this->SetResliceMatrix();
    // Set the point through which to slice
    this->resliceAxes->SetElement( 0, 3, this->resliceOrigin[0] );
    this->resliceAxes->SetElement( 1, 3, this->resliceOrigin[1] );
    this->resliceAxes->SetElement( 2, 3, this->resliceOrigin[2] );
    this->resliceAxes->SetElement( 3, 0, 0 );
    this->resliceAxes->SetElement( 3, 1, 0 );
    this->resliceAxes->SetElement( 3, 2, 0 );
    this->resliceAxes->SetElement( 3, 3, 1 );

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

    // Convert the VTK image to an ITK image for further processing
    vtkSmartPointer<vtkImageCast> ic = vtkSmartPointer<vtkImageCast>::New();
    ic->SetInputConnection(reslice->GetOutputPort());
    ic->SetOutputScalarTypeToUnsignedChar(); //unsigned short to unsigned char DIMENSIONS!!
    ic->Update();
    this->vtkImageToImageFilter->SetInput( ic->GetOutput() );
    this->vtkImageToImageFilter->Update();
    double tempSpac[3];
    this->reslice->GetOutput()->GetSpacing( tempSpac );
    this->resliceSpacing[0] = tempSpac[0];
    this->resliceSpacing[1] = tempSpac[1];
    this->resliceSpacing[2] = tempSpac[2];

    this->reslicedImage.imageData->Graft( this->vtkImageToImageFilter->GetOutput() );

    /*QuickView viewer;
    viewer.AddImage( this->reslicedImage.imageData.GetPointer() ); // Need to do this because QuickView can't accept smart pointers
    viewer.Visualize();*/

    //this->SetSpacingOfReslice();
    this->resliceMatrix.SetSpacingForIGTMatrix( this->resliceSpacing );
    int sizeImage[3] = {1, 1, 1};
    int tempSize[3];
    ic->GetOutput()->GetDimensions( tempSize );
    sizeImage[0] = tempSize[0]; sizeImage[1] = tempSize[1];
    this->resliceMatrix.SetDimensionsForIGTMatrix( sizeImage );
    this->resliceMatrix.matrix = this->resliceMatrix.matrix;
    this->resliceMatrix.SetIGTTransformFromMat();
   // std::cerr<<"reslicematrix"<<std::endl;
    //this->resliceMatrix.ShowMatrix();

    // Set correct image parameters for the ITK image

    this->reslicedImage.SetParametersFromITK(this->resliceOrigin[2], this->resliceSpacing[2], this->resliceMatrix); //Set spacing and origin in imageData

  }
  else
  {
    std::cerr << "ITK Image hasn't been created" << std::endl;
  }
  //std::cerr << resliceOrigin[2] << std::endl;\

  return;

}
