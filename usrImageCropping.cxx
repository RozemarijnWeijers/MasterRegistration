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
  ImageType::PointType          origin;
  float                         originSliceImage[3];
  double                        spacingSliceImage[3];

  float* spacing = this->inputImage->spacingImage;
  float* start1 = this->inputImage->originImage;
  spacingSliceImage[0] = spacing[0]; spacingSliceImage[1] = spacing[1]; spacingSliceImage[2] = spacing[2];
  originSliceImage[0] = start1[0] + desiredStart[0]; originSliceImage[1] = start1[1] + desiredStart[1];   originSliceImage[2] = start1[2] + desiredStart[2];


  this->croppedImage.imageData->SetSpacing( spacing );
  this->croppedImage.imageData->SetOrigin( originSliceImage );
  this->SetImageMatrix();

  this->croppedImage.SetParametersFromITK( originSliceImage[2], spacingSliceImage[2], croppedImage.imageMatrix );
  std::cerr<< spacing[0] <<std::endl;

  /*QuickView viewer;
  viewer.AddImage( this->croppedImage.imageData.GetPointer() ); // Need to do this because QuickView can't accept smart pointers
  viewer.Visualize();*/

  return;

}

void ImageCropping::SetImageMatrix()
{

  float* spacing = this->inputImage->spacingImage;
  this->croppedImage.imageMatrix.matrix = this->inputImage->imageMatrix.matrix;
  this->croppedImage.imageMatrix.SetDimensionsForIGTMatrix( croppedImage.sizeImage );
  this->croppedImage.imageMatrix.SetSpacingForIGTMatrix( spacing );
  this->croppedImage.imageMatrix.SetIGTTransformFromMat();

  TransformMatrix cropMatrix;
  cropMatrix.SetDimensionsForIGTMatrix( croppedImage.sizeImage );
  cropMatrix.SetSpacingForIGTMatrix( spacing );
  float dStart[3] = { this->desiredStart[0], this->desiredStart[1], 0 };
  cropMatrix.SetOriginInTransform( dStart );

  this->croppedImage.imageMatrix.ShowMatrix();
  this->croppedImage.imageMatrix.MultiplyWith( cropMatrix.matrix );

  std::cerr<<"cropped imageMatrix:"<< std::endl;
  this->croppedImage.imageMatrix.ShowMatrix();

  return;

}

void ImageCropping::Convert2DImageTo3DVolume()
{

  Filter2DTo3DType::Pointer filter2DTo3D = Filter2DTo3DType::New();

  itk::FixedArray< unsigned char, OutputDimension > layout;
  layout[0] = 1;
  layout[1] = 1;
  layout[2] = 0;

  filter2DTo3D->SetLayout( layout );

  ImageType::Pointer input = this->croppedImage.imageData;

  filter2DTo3D->SetInput( 0, input );

  //typedef unsigned char                             PixelType;
  //const PixelType defaultValue = 128;

  //filter->SetDefaultPixelValue( defaultValue );

  return;

}
