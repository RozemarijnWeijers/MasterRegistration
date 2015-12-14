# include "usrVolume.h"

Volume::Volume()
{

  // Create imageMessage and ITKvolume for this volume
  volumeData  = VolumeType::New();
  VTKReader = vtkSmartPointer<vtkNrrdReader>::New();
  imgMsg = igtl::ImageMessage::New();

}

Volume::~Volume()
{
}

void Volume::SetParametersFromITK()
{

  VolumeType::SizeType      size = volumeData->GetLargestPossibleRegion().GetSize();
  sizeVolume[0] = size[0];          sizeVolume[1] = size[1];            sizeVolume[2] = size[2];
  VolumeType::PointType     origin = volumeData->GetOrigin();
  originVolume[0] = origin[0];      originVolume[1] = origin[1];        originVolume[2] = origin[2];
  VolumeType::SpacingType   spacing = volumeData->GetSpacing();
  spacingVolume[0] = spacing[0];    spacingVolume[1] = spacing[1];      spacingVolume[2] = spacing[2]; //laatste is 0

  return;

}

void Volume::ConvertITKtoIGTVolume()
{

  // convert the origin (IGT uses the centre as origin instead of the corner)
  float originIGT[3];
  for ( int i=0 ; i<3; i++)
  {
    originIGT[i] = originVolume[i]+((sizeVolume[i]-1)*spacingVolume[i]/2);
  }

  int scalarType = igtl::ImageMessage::TYPE_UINT8;

  // Set volume data to image message
  imgMsg->SetDimensions( sizeVolume );
  imgMsg->SetSpacing( spacingVolume );
  imgMsg->SetScalarType( scalarType );
  imgMsg->SetOrigin( originIGT );
  imgMsg->AllocateScalars();

  // Set copy image data into ITK image
  memcpy( imgMsg->GetScalarPointer(), volumeData->GetBufferPointer(), imgMsg->GetImageSize() );
  // Pack (serialize) and send
  imgMsg->Pack();

  return;

}

int Volume::LoadVolume( char* filename )
{

  // Open file with VTKreader and set NRRD file as input (for reslicing)
  VTKReader->SetFileName( filename );
  VTKReader->Update();

  // Open file ITKreader and set NRRD file as input
  typedef itk::ImageFileReader<VolumeType> FileReaderType;
  FileReaderType::Pointer ITKreader = FileReaderType::New();
  ITKreader->SetFileName(filename);

  try
  {
    // Read volume from file
    ITKreader->Update();
    volumeData = ITKreader->GetOutput();
  }
  catch (itk::ExceptionObject &ex)
  {
    std::cerr << ex << std::endl;
    return 0;
  }
  SetParametersFromITK();

  std::cerr << "Volume Loaded " << std::endl;

  return 1;

}
