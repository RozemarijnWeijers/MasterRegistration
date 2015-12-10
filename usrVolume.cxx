# include "usrVolume.h"

Volumes::Volumes()
{

  // Create imageMessage and ITKvolume for this volume
  volumeData  = VolumeType::New();
  VTKreader = vtkSmartPointer<vtkNrrdReader>::New();
  imgMsg = igtl::ImageMessage::New();

}

void Volumes::SetParametersFromITK()
{

  VolumeType::SizeType              size = this->volumeData->GetLargestPossibleRegion().GetSize();
  this->sizeVol[0] = size[0];       this->sizeVol[1] = size[1];        this->sizeVol[2] = size[2];
  VolumeType::PointType             origin = this->volumeData->GetOrigin();
  this->originVol[0] = origin[0];   this->originVol[1] = origin[1];    this->originVol[2] = origin[2];
  VolumeType::SpacingType           spacing = this->volumeData->GetSpacing();
  this->spacingVol[0] = spacing[0]; this->spacingVol[1] = spacing[1];  this->spacingVol[2] = spacing[2]; //laatste is 0

  return;

}

void Volumes::SetParametersFromVTK( double origin[3], double spacing[3], int size[3])
{

  this->sizeVol[0] = size[0];       this->sizeVol[1] = size[1];        this->sizeVol[2] = size[2];
  this->originVol[0] = origin[0];   this->originVol[1] = origin[1];    this->originVol[2] = origin[2];
  this->spacingVol[0] = spacing[0]; this->spacingVol[1] = spacing[1];  this->spacingVol[2] = spacing[2]; //laatste is 0

  return;

}

int Volumes::ITKtoIGTVolume()
{

  // Retrieve the image data from ITK image
  /*VolumeType::SizeType size;
  VolumeType::IndexType start;
  VolumeType::PointType origin;
  VolumeType::SpacingType spacing;

  size = this->volumeData->GetLargestPossibleRegion().GetSize();
  start = this->volumeData->GetLargestPossibleRegion().GetIndex();
  spacing = this->volumeData->GetSpacing();
  origin = this->volumeData->GetOrigin();*/

  // Set volume data to image message
  int   sizeIGT[3];
  int   scalarType;
  float originIGT[3];
  float spacingIGT[3];    // spacing (mm/pixel)
  //int   svsize[3];      // sub-volume size
  //int   svoffset[3];    // sub-volume offset
  sizeIGT[0] = this->sizeVol[0];         sizeIGT[1] = this->sizeVol[1];           sizeIGT[2] = this->sizeVol[2];//sizeIGT[0] = size[0];         sizeIGT[1] = size[1];           sizeIGT[2] = 1;
  spacingIGT[0] = this->spacingVol[0];   spacingIGT[1] = this->spacingVol[1];     spacingIGT[2] = this->spacingVol[2];//spacing[0];   spacingIGT[1] = spacing[1];     spacingIGT[2] = spacing[1]; // third dimension?
  originIGT[0] = this->originVol[0]+((this->sizeVol[0]-1)*this->spacingVol[0]/2);          originIGT[1] = this->originVol[1]+((this->sizeVol[1]-1)*this->spacingVol[1]/2);     originIGT[2] = this->originVol[2]+((this->sizeVol[2]-1)*this->spacingVol[2]/2); // klopt nog niet
  //svsize[0] = sizeIGT[0];     svsize[1] = sizeIGT[1];     svsize[2] = sizeIGT[2];
  //svoffset[0] = start[0];     svoffset[1] = start[1];     svoffset[2] = start[2];
  scalarType = igtl::ImageMessage::TYPE_UINT8;

  this->imgMsg->SetDimensions( sizeIGT );
  this->imgMsg->SetSpacing( spacingIGT );
  this->imgMsg->SetScalarType( scalarType );
  this->imgMsg->SetOrigin( originIGT );
  this->imgMsg->AllocateScalars();

  // Set copy image data into ITK image
  memcpy( this->imgMsg->GetScalarPointer(), volumeData->GetBufferPointer(), this->imgMsg->GetImageSize() );
  // Pack (serialize) and send
  this->imgMsg->Pack();

  return 1;

}

int Volumes::LoadVolume( char* filename )
{

  // Open file VTKreader and set NRRD file as input
  this->VTKreader->SetFileName( filename );
  this->VTKreader->Update();

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

  std::cerr << "VTKVolume Loaded " << std::endl;

  return 1;

}
