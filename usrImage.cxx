# include "usrImage.h"

Images::Images()
{

  // Create imageMessage and ITKimage for this image
  imageData  = ImageType::New();
  imgMsg = igtl::ImageMessage::New();

}

void Images::SetParametersFromIGT()
{

  this->imgMsg->GetDimensions( this->sizeIm );
  this->imgMsg->GetOrigin( this->originIm );
  this->imgMsg->GetSpacing( this->spacingIm );
  this->imgMsg->GetMatrix( this->matrixIm );

  return;

}

void Images::SetParametersFromITK( double origin[3], double spacing[3] )
{

  ImageType::SizeType               size = this->imageData->GetLargestPossibleRegion().GetSize();
  this->sizeIm[0] = size[0];        this->sizeIm[1] = size[1];          this->sizeIm[2] = 1;
  this->originIm[0] = origin[0];    this->originIm[1] = origin[1];      this->originIm[2] = origin[2];
  this->spacingIm[0] = spacing[0];  this->spacingIm[1] = spacing[1];    this->spacingIm[2] = spacing[2]; //laatste is 0

  return;

}

int Images::IGTtoITKImage()
{

  // Retrieve the image data from image message
  //int       size[3];          // image dimension (pixels)
  //float     spacing[3];       // spacing (mm/pixel)
  float     spacingITK[0];    // spacing for 2D image (mm/pixel)
  //int       endian;           // endian (not used)
  //float     origin[3];        // origin ()
  //igtl::Matrix4x4   matrix;   // image origin and orientation matrix

  //endian = this->imgMsg->GetEndian();
  //this->imgMsg->GetDimensions( size );
  //this->imgMsg->GetSpacing( spacing );
  //this->imgMsg->GetMatrix( matrix );
  //this->imgMsg->GetOrigin( origin );

  // Set image data to ITK image (3D information -> 2D)
  ImageType::RegionType     region;
  ImageType::IndexType      start;
  start[0] = this->originIm[0];     start[1] = this->originIm[1];//start[0] = origin[0];     start[1] = origin[1];
  ImageType::SizeType       sizeregion;
  sizeregion[0] = this->sizeIm[0];  sizeregion[1] = sizeIm[1];//sizeregion[0] = size[0];  sizeregion[1] = size[1];
  region.SetSize( sizeregion );
  region.SetIndex( start );
  spacingITK[0] = this->spacingIm[0]; spacingITK[1] = this->spacingIm[1];//spacingITK[0] = spacing[0]; spacingITK[1] = spacing[1];

  this->imageData->SetRegions( region );
  this->imageData->SetSpacing( spacingITK );
  this->imageData->Allocate();

  // Copy image data into ITK image
  memcpy( this->imageData->GetBufferPointer(), this->imgMsg->GetScalarPointer(), this->imgMsg->GetSubVolumeImageSize() );
  this->imageData->Modified();

  return 1;

}

int Images::ITKtoIGTImage()
{

  // Retrieve the image data from ITK image
  ImageType::SizeType       size;
  //ImageType::IndexType      start;
  //ImageType::SpacingType    spacing;
  //ImageType::PointType      origin;

  size = this->imageData->GetLargestPossibleRegion().GetSize();
  //start = this->imageData->GetLargestPossibleRegion().GetIndex();
  //spacing = this->imageData->GetSpacing();
  //origin = this->imageData->GetOrigin();

  // Set image data to image message (2D information -> 3D)
  int                       sizeIGT[3];
  ImageType::RegionType     region;
  float                     originIGT[3];
  float                     spacingIGT[3];  // spacing (mm/pixel)
  int                       scalarType;     // always UINT8
  //int                     svoffset[3];    // sub-volume offset
  //int                     svsize[3];      // sub-volume size
  //sizeIGT[0] = this->sizeIm[0];         sizeIGT[1] = this->sizeIm[1];           sizeIGT[2] = 1;//sizeIGT[0] = size[0];         sizeIGT[1] = size[1];           sizeIGT[2] = 1;
  sizeIGT[0] = size[0];         sizeIGT[1] = size[1];           sizeIGT[2] = 1; //WAAROM??
  spacingIGT[0] = this->spacingIm[0];   spacingIGT[1] = this->spacingIm[1];     spacingIGT[2] = this->spacingIm[2];//spacing[0];   spacingIGT[1] = spacing[1];     spacingIGT[2] = spacing[1]; // third dimension?
  //spacingIGT[0] = spacing[0];   spacingIGT[1] = spacing[1];     spacingIGT[2] = spacing[1]; // third dimension?
  originIGT[0] = this->originIm[0]+((this->sizeIm[0]-1)*this->spacingIm[0]/2);          originIGT[1] = this->originIm[1]+((this->sizeIm[1]-1)*this->spacingIm[1]/2);     originIGT[2] = this->originIm[2];//+((size[2]-1)*spacing[2]/2); // klopt nog niet
  //originIGT[0] = origin[0]+((size[0]-1)*spacing[0]/2);          originIGT[1] = origin[1]+((size[1]-1)*spacing[1]/2);     originIGT[2] = origin[2];//+((size[2]-1)*spacing[2]/2); // klopt nog niet
  //svsize[0] = sizeIGT[0];     svsize[1] = sizeIGT[1];         svsize[2] = sizeIGT[2];
  //svoffset[0] = start[0];     svoffset[1] = start[1];         svoffset[2] = start[2];
  scalarType = igtl::ImageMessage::TYPE_UINT8;

  this->imgMsg->SetDimensions( sizeIGT );
  this->imgMsg->SetSpacing( spacingIGT );
  this->imgMsg->SetScalarType( scalarType );
  this->imgMsg->SetOrigin( originIGT );
  //this->imgMsg->SetSubVolume( sizeIGT, svoffset );
  this->imgMsg->AllocateScalars();

  // Copy image data into ITK image
  memcpy( this->imgMsg->GetScalarPointer(), imageData->GetBufferPointer(), this->imgMsg->GetSubVolumeImageSize() );
  // Pack (serialize) and send
  this->imgMsg->Pack();

  return 1;

}
