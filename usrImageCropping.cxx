#include "usrImageCropping.h"

ImageCropping::ImageCropping()
{

  filter = FilterImageType::New();

}

ImageCropping::~ImageCropping()
{

}

void ImageCropping::SetCropSizeAndStart( int dSize[2], int dStart[2] )
{

  this->desiredStart[0] = dStart[0];  this->desiredStart[1] = dStart[1];
  this->desiredSize[0] = dSize[0];    this->desiredSize[1] = dSize[1];

  return;

}

void ImageCropping::SetImage( Image* inputImage )
{

  this->inputImage = inputImage;

  return;

}

void ImageCropping::CropImage()
{

  ImageType::RegionType        desiredRegion( this->desiredStart, this->desiredSize );

  // Create cropping/reslice

  this->filter->SetExtractionRegion( desiredRegion );
  this->filter->SetInput( this->inputImage->imageData );
  this->filter->SetDirectionCollapseToIdentity(); // This is required.
  this->filter->Update();

  // Set resliced image to ITK image
  this->croppedImage.imageData = this->filter->GetOutput();

  // Get and set parameters for reslices image    // spacing (mm/pixel)
  ImageType::PointType         origin;
  float                         originSliceImage[3];
  double                        spacingSliceImage[3];

  float* spacing = this->inputImage->spacingImage;
  float* start1 = this->inputImage->originImage;
  spacingSliceImage[0] = spacing[0]; spacingSliceImage[1] = spacing[1]; spacingSliceImage[2] = spacing[2];
  originSliceImage[0] = start1[0] + desiredStart[0]; originSliceImage[1] = start1[1] + desiredStart[1];   originSliceImage[2] = start1[2] + desiredStart[2];

  this->croppedImage.imageData->SetSpacing( spacing );
  this->croppedImage.imageData->SetOrigin( originSliceImage );
  int tempSize[3];
  /*tempSize[0] = this->desiredSize[0]; tempSize[1] = this->desiredSize[1]; tempSize[2] = this->desiredSize[2];
  this->croppedImage.volumeMatrix.SetDimensionsForIGTMatrix( tempSize );
  this->croppedImage.volumeMatrix.SetSpacingForIGTMatrix( spacing );
  this->croppedImage.volumeMatrix.SetOriginInTransform( originSliceImage );
  this->croppedImage.SetParametersFromITK();//( originSliceImage[2], spacingSliceImage[2], croppedVolume.volumeMatrix );
*/
  return;

}
