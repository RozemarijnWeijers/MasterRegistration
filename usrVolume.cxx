# include "usrVolume.h"

Volume::Volume()
{

  // Create imageMessage and ITKvolume for this volume
  this->volumeData  = VolumeType::New();
  this->VTKReader = vtkSmartPointer<vtkNrrdReader>::New();
  this->imgMsg = igtl::ImageMessage::New();

}

Volume::~Volume()
{
}

void Volume::SetParametersFromITK()
{

  VolumeType::SizeType      size = volumeData->GetLargestPossibleRegion().GetSize();
  this->sizeVolume[0] = size[0];          this->sizeVolume[1] = size[1];            this->sizeVolume[2] = size[2];
  VolumeType::PointType     origin = volumeData->GetOrigin();
  this->originVolume[0] = origin[0];      this->originVolume[1] = origin[1];        this->originVolume[2] = origin[2];
  VolumeType::SpacingType   spacing = volumeData->GetSpacing();
  this->spacingVolume[0] = spacing[0];    this->spacingVolume[1] = spacing[1];      this->spacingVolume[2] = spacing[2]; //laatste is 0

  return;

}

void Volume::ConvertITKtoIGTVolume()
{

  // convert the origin (IGT uses the centre as origin instead of the corner)
  float originIGT[3];
  for ( int i=0; i<3; i++)
  {
    originIGT[i] = this->originVolume[i]+((this->sizeVolume[i]-1)*this->spacingVolume[i]/2);
  }

  int scalarType = igtl::ImageMessage::TYPE_UINT8;

  // Set volume data to image message
  this->imgMsg->SetDimensions( this->sizeVolume );
  this->imgMsg->SetSpacing( this->spacingVolume );
  this->imgMsg->SetScalarType( scalarType );
  this->imgMsg->SetOrigin( originIGT );
  this->imgMsg->AllocateScalars();

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

void Volume::UpdateVolumeTranform( double transform[16] )
{
}
