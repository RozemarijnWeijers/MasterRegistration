# include "usrVolume.h"

Volume::Volume()
{

  // Create imageMessage and ITKvolume for this volume
  this->volumeData  = VolumeType::New();
  this->VTKReader = vtkSmartPointer<vtkNrrdReader>::New();
  this->imgMsg = igtl::ImageMessage::New();
  //initialize transformMatrix

}

Volume::~Volume()
{
}

void Volume::SetParametersFromITK()// TransformMatrix ITKimageMatrix )
{

  VolumeType::SizeType      size = this->volumeData->GetLargestPossibleRegion().GetSize();
  VolumeType::PointType     origin = this->volumeData->GetOrigin();
  VolumeType::SpacingType   spacing = this->volumeData->GetSpacing();
  this->sizeVolume[0] = size[0];          this->sizeVolume[1] = size[1];            this->sizeVolume[2] = size[2];
  this->originVolume[0] = origin[0];      this->originVolume[1] = origin[1];        this->originVolume[2] = origin[2];
  this->spacingVolume[0] = spacing[0];    this->spacingVolume[1] = spacing[1];      this->spacingVolume[2] = spacing[2];

  this->volumeMatrix.SetDimensionsForIGTMatrix( this->sizeVolume );
  this->volumeMatrix.SetSpacingForIGTMatrix( this->spacingVolume );

  VolumeType::DirectionType direction = this->volumeData->GetDirection();
  double dir[9];
  dir[0] = direction[0][0]; dir[1] = direction[0][1]; dir[2] = direction[0][2];
  dir[3] = direction[1][0]; dir[4] = direction[1][1]; dir[5] = direction[1][2];
  dir[6] = direction[2][0]; dir[7] = direction[2][1]; dir[8] = direction[2][2];
  this->volumeMatrix.SetDirectionInTransform( dir );
  this->volumeMatrix.SetOriginInTransform( this->originVolume );

  return;

}

void Volume::ConvertITKtoIGTVolume()
{

  // convert the origin (IGT uses the centre as origin instead of the corner)
  /*float originIGT[3];
  for ( int i=0; i<3; i++)
  {
    originIGT[i] = this->originVolume[i]+((this->sizeVolume[i]-1)*this->spacingVolume[i]/2);
  }*/

  int scalarType = igtl::ImageMessage::TYPE_UINT8;

  // Set volume data to image message
  this->imgMsg->SetDimensions( this->sizeVolume );
  this->imgMsg->SetSpacing( this->spacingVolume );
  this->imgMsg->SetScalarType( scalarType );
  //this->imgMsg->SetOrigin( originIGT );
  this->imgMsg->SetMatrix( this->volumeMatrix.IGTMatrix );
  this->imgMsg->AllocateScalars();
  //convert transformMatrix

  // Set copy image data into ITK image
  memcpy( this->imgMsg->GetScalarPointer(), this->volumeData->GetBufferPointer(), this->imgMsg->GetImageSize() );
  // Pack (serialize) and send
  this->imgMsg->Pack();

  return;

}

int Volume::LoadVolume( char* filename )
{

  // Open file with VTKreader and set NRRD file as input (for reslicing)
  this->VTKReader->SetFileName( filename );
  this->VTKReader->Update();

  // Open file ITKreader and set NRRD file as input
  typedef itk::ImageFileReader<VolumeType> FileReaderType;
  FileReaderType::Pointer ITKreader = FileReaderType::New();
  ITKreader->SetFileName(filename);

  try
  {
    // Read volume from file
    ITKreader->Update();
    this->volumeData = ITKreader->GetOutput();
  }
  catch (itk::ExceptionObject &ex)
  {
    std::cerr << ex << std::endl;
    return 0;
  }

  this->SetParametersFromITK();

  std::cerr << "Volume Loaded " << std::endl;

  return 1;

}

void Volume::UpdateVolumeTransform( TransformMatrix updateMatrix )
{

  this->volumeMatrix.matrix = updateMatrix.matrix;
  this->volumeMatrix.SetIGTTransformFromMat( updateMatrix.matrix );
  VolumeType::DirectionType direction ;
  direction[0][0] = updateMatrix.matrix(0,0); direction[1][0] = updateMatrix.matrix(1,0); direction[2][0] = updateMatrix.matrix(2,0);
  direction[1][0] = updateMatrix.matrix(0,1); direction[1][1] = updateMatrix.matrix(1,1); direction[2][1] = updateMatrix.matrix(2,1);
  direction[2][0] = updateMatrix.matrix(0,2); direction[1][2] = updateMatrix.matrix(1,2); direction[2][2] = updateMatrix.matrix(2,2);
  this->volumeData->SetDirection( direction );
  int origin[1];
  origin[0] = updateMatrix.matrix(3,0); origin[1] = updateMatrix.matrix(3,1); origin[2] = updateMatrix.matrix(3,2);
  this->volumeData->SetOrigin( origin );

  return;

}

void Volume::CropVolume( int dStart[3], int dSize[3], Volume* croppedVolume )
{

  // Set parameters for desired image (reslice/ crop)
  VolumeType::IndexType         desiredStart;
  VolumeType::SizeType          desiredSize ;
  desiredStart[0] = dStart[0];  desiredStart[1] = dStart[1];    desiredStart[2] = dStart[2];
  desiredSize[0] = dSize[0];    desiredSize[1] = dSize[1];      desiredSize[2] = dSize[2];

  VolumeType::RegionType        desiredRegion( desiredStart, desiredSize );
  std::cout << "Desired Region: " << desiredRegion << std::endl;

  // Create cropping/reslice
  FilterType::Pointer           filter = FilterType::New();
  filter->SetExtractionRegion( desiredRegion );
  filter->SetInput( this->volumeData );
  filter->SetDirectionCollapseToIdentity(); // This is required.
  filter->Update();

  // Set resliced image to ITK image
  croppedVolume->volumeData = filter->GetOutput();

  // Get and set parameters for reslices image    // spacing (mm/pixel)

  VolumeType::PointType         origin;
  double                        originsliceVolume[3];
  double                        spacingsliceVolume[3];

  float* spacing = this->spacingVolume;
  float* start1 = this->originVolume;
  spacingsliceVolume[0] = spacing[0]; spacingsliceVolume[1] = spacing[1]; spacingsliceVolume[2] = spacing[2];
  originsliceVolume[0] = start1[0] + dStart[0] * spacing[0]; originsliceVolume[1] = start1[1] + dStart[1] * spacing[1];   originsliceVolume[2] = start1[2] + dStart[2] * spacing[2];

  croppedVolume->volumeData->SetSpacing( spacing );
  croppedVolume->volumeData->SetOrigin( originsliceVolume );
  croppedVolume->SetParametersFromITK();

  return;

}
