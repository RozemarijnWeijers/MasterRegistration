#include "usrVolumeCropping.h"

VolumeCropping::VolumeCropping()
{

  filter = FilterType::New();

}

VolumeCropping::~VolumeCropping()
{

}

void VolumeCropping::SetCropSizeAndStart( int dSize[3], int dStart[3] )
{

  this->desiredStart[0] = dStart[0];  this->desiredStart[1] = dStart[1];    this->desiredStart[2] = dStart[2];
  this->desiredSize[0] = dSize[0];    this->desiredSize[1] = dSize[1];      this->desiredSize[2] = dSize[1];

  return;

}

void VolumeCropping::SetVolume( Volume* inputVolume )
{

  this->inputVolume = inputVolume;

  return;

}

void VolumeCropping::CropVolume()
{

  VolumeType::RegionType        desiredRegion( this->desiredStart, this->desiredSize );

  // Create cropping/reslice

  this->filter->SetExtractionRegion( desiredRegion );
  this->filter->SetInput( this->inputVolume->volumeData );
  this->filter->SetDirectionCollapseToIdentity(); // This is required.
  this->filter->Update();

  // Set resliced image to ITK image
  this->croppedVolume.volumeData = this->filter->GetOutput();

  // Get and set parameters for reslices image    // spacing (mm/pixel)
  VolumeType::PointType         origin;
  float                         originSliceImage[3];
  double                        spacingSliceImage[3];

  float* spacing = this->inputVolume->spacingVolume;
  float* start1 = this->inputVolume->originVolume;
  spacingSliceImage[0] = spacing[0]; spacingSliceImage[1] = spacing[1]; spacingSliceImage[2] = spacing[2];
  originSliceImage[0] = start1[0] + desiredStart[0]; originSliceImage[1] = start1[1] + desiredStart[1];   originSliceImage[2] = start1[2] + desiredStart[2];

  std::cerr<< spacing[0]<< std::endl;
  this->croppedVolume.volumeData->SetSpacing( spacing );
  this->croppedVolume.volumeData->SetOrigin( originSliceImage );
  int tempSize[3];
  tempSize[0] = this->desiredSize[0]; tempSize[1] = this->desiredSize[1]; tempSize[2] = this->desiredSize[2];
  this->croppedVolume.volumeMatrix.SetDimensionsForIGTMatrix( tempSize );
  this->croppedVolume.volumeMatrix.SetSpacingForIGTMatrix( spacing );
  this->croppedVolume.volumeMatrix.SetOriginInTransform( originSliceImage );
  this->croppedVolume.SetParametersFromITK();//( originSliceImage[2], spacingSliceImage[2], croppedVolume.volumeMatrix );

  return;

}
