#include "usrImage.h"

Image::Image()
{

  // Create imageMessage and ITKimage for this image
  this->imageData  = ImageType::New();
  this->imgMsg = igtl::ImageMessage::New();

}

Image::~Image()
{
}

void Image::SetParametersFromIGT()
{

  //is this right?
  this->imgMsg->GetDimensions( this->sizeImage );
  this->imgMsg->GetOrigin( this->originImage );
  this->imgMsg->GetSpacing( this->spacingImage );

  this->imageMatrix.SetDimensionsForIGTMatrix( this->sizeImage );
  this->imageMatrix.SetSpacingForIGTMatrix( this->spacingImage );
  this->imageMatrix.SetTransformFromIGT( this->imgMsg );

  return;

}

void Image::SetParametersFromITK( double origin3th, double spacing3th, TransformMatrix ITKimageMatrix, bool RAS )
{

  ImageType::PointType           origin = this->imageData->GetOrigin();
  ImageType::SpacingType         spacing = this->imageData->GetSpacing();
  ImageType::SizeType            size = this->imageData->GetLargestPossibleRegion().GetSize();
  this->sizeImage[0] = size[0];        this->sizeImage[1] = size[1];          this->sizeImage[2] = 1;
  this->originImage[0] = origin[0];    this->originImage[1] = origin[1];      this->originImage[2] = origin3th;  //Waarom deze in arguments????
  this->spacingImage[0] = spacing[0];  this->spacingImage[1] = spacing[1];    this->spacingImage[2] = spacing3th; //laatste is 0

  this->imageMatrix = ITKimageMatrix;
  this->imageMatrix.SetDimensionsForIGTMatrix( this->sizeImage );
  this->imageMatrix.SetSpacingForIGTMatrix( this->spacingImage );

  if ( RAS == false )
  {
    //this->originImage[0] = -origin[0];    this->originImage[1] = -origin[1];      this->originImage[2] = origin3th;  //Waarom deze in arguments????
    //this->imageMatrix.SetOriginInTransform( this->originImage );
    this->imageMatrix.SetLPS();
  }

  return;

}

int Image::ConvertIGTtoITKImage()
{

  this->SetParametersFromIGT();

  // Set image data to ITK image (3D information -> 2D)
  ImageType::RegionType     region;
  ImageType::IndexType      start;
  ImageType::SizeType       sizeregion;

  start[0] = originImage[0];     start[1] = originImage[1];
  sizeregion[0] = sizeImage[0];  sizeregion[1] = sizeImage[1];
  region.SetSize( sizeregion );
  region.SetIndex( start );

  this->imageData->SetRegions( region );
  this->imageData->SetSpacing( this->spacingImage );
  this->imageData->Allocate();

  // Copy image data into ITK image
  memcpy( this->imageData->GetBufferPointer(), this->imgMsg->GetScalarPointer(), this->imgMsg->GetSubVolumeImageSize() );
  this->imageData->Modified();

  return 1;

}

int Image::ConvertITKtoIGTImage()
{

  // SetParametersFromITK( origin3th, spacing3th );
  // Convert the origin (IGT uses the centre as origin instead of the corner) (2D information -> 3D)
  /*float originIGT[3];
  for ( int i=0; i<3; i++)
  {
    originIGT[i] = this->originImage[i]+((this->sizeImage[i]-1)*this->spacingImage[i]/2);
  }*/

  int scalarType = igtl::ImageMessage::TYPE_UINT8;

  this->imgMsg->SetDimensions( this->sizeImage );
  this->imgMsg->SetSpacing( this->spacingImage );
  this->imgMsg->SetScalarType( scalarType );
  this->imgMsg->SetMatrix( this->imageMatrix.IGTMatrix );
  this->imgMsg->AllocateScalars();

  // Copy image data into ITK image
  memcpy( this->imgMsg->GetScalarPointer(), this->imageData->GetBufferPointer(), this->imgMsg->GetSubVolumeImageSize() );
  // Pack (serialize) and send
  this->imgMsg->Pack();

  return 1;

}

