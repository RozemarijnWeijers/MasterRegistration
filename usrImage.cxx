#include "usrImage.h"

Image::Image()
{

  // Create imageMessage and ITKimage for this image
  imageData  = ImageType::New();
  imgMsg = igtl::ImageMessage::New();

}

Image::~Image()
{
}

void Image::SetParametersFromIGT()
{

  imgMsg->GetDimensions( sizeImage );
  imgMsg->GetOrigin( originImage );
  imgMsg->GetSpacing( spacingImage );
  imgMsg->GetMatrix( matrixImage );

  return;

}

void Image::SetParametersFromITK( double origin[3], double spacing[3] )
{

  ImageType::SizeType            size = imageData->GetLargestPossibleRegion().GetSize();
  sizeImage[0] = size[0];        sizeImage[1] = size[1];          sizeImage[2] = 1;
  originImage[0] = origin[0];    originImage[1] = origin[1];      originImage[2] = origin[2];  //Waarom deze in arguments????
  spacingImage[0] = spacing[0];  spacingImage[1] = spacing[1];    spacingImage[2] = spacing[2]; //laatste is 0

  return;

}

int Image::ConvertIGTtoITKImage()
{

  // Set image data to ITK image (3D information -> 2D)
  ImageType::RegionType     region;
  ImageType::IndexType      start;
  ImageType::SizeType       sizeregion;

  start[0] = originImage[0];     start[1] = originImage[1];//start[0] = origin[0];     start[1] = origin[1];
  sizeregion[0] = sizeImage[0];  sizeregion[1] = sizeImage[1];//sizeregion[0] = size[0];  sizeregion[1] = size[1];
  region.SetSize( sizeregion );
  region.SetIndex( start );

  imageData->SetRegions( region );
  imageData->SetSpacing( spacingImage );
  imageData->Allocate();

  // Copy image data into ITK image
  memcpy( imageData->GetBufferPointer(), imgMsg->GetScalarPointer(), imgMsg->GetSubVolumeImageSize() );
  imageData->Modified();

  return 1;

}

int Image::ConvertITKtoIGTImage()
{

  // Convert the origin (IGT uses the centre as origin instead of the corner) (2D information -> 3D)
  float originIGT[3];
  for ( int i=0; i<3; i++)
  {
    originIGT[i] = originImage[i]+((sizeImage[i]-1)*spacingImage[i]/2);
  }

  int scalarType = igtl::ImageMessage::TYPE_UINT8;

  imgMsg->SetDimensions( sizeImage );
  imgMsg->SetSpacing( spacingImage );
  imgMsg->SetScalarType( scalarType );
  imgMsg->SetOrigin( originIGT );
  imgMsg->AllocateScalars();

  // Copy image data into ITK image
  memcpy( imgMsg->GetScalarPointer(), imageData->GetBufferPointer(), imgMsg->GetSubVolumeImageSize() );
  // Pack (serialize) and send
  imgMsg->Pack();
  std::cerr<< "reslicing is done" << std::endl;
  return 1;

}

