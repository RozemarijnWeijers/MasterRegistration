#include "usrImageCropping.h"

ImageCropping::ImageCropping()
{

  filter = FilterImageType::New();

}

ImageCropping::~ImageCropping()
{

}

void ImageCropping::SetNumberofImages( int n )
{

  this->number = n;

  return;

}

void ImageCropping::SetCropSizeAndStart( int dSize[2], float dStart[2] )
{

  this->desiredStart[0] = dStart[0];  this->desiredStart[1] = dStart[1]; // in mm
  this->desiredSize[0] = dSize[0];    this->desiredSize[1] = dSize[1]; // in pixels

  return;

}

void ImageCropping::SetImage( Image* inputImage )
{

  this->inputImage = inputImage;

  return;

}

void ImageCropping::CropImage()
{

  ImageType::IndexType start;
  start[0] = this->desiredStart[0]/this->inputImage->spacingImage[0];
  start[1] = this->desiredStart[1]/this->inputImage->spacingImage[1];
  ImageType::RegionType        desiredRegion(start, this->desiredSize );

  // Create cropping/reslice

  this->filter->SetExtractionRegion( desiredRegion );
  this->filter->SetInput( this->inputImage->imageData );
  this->filter->SetDirectionCollapseToIdentity(); // This is required.
  this->filter->Update();

  // Set resliced image to ITK image
  this->croppedImage.imageData = this->filter->GetOutput();

  // Get and set parameters for reslices image
  ImageType::PointType          origin;
  float                         originSliceImage[3];
  double                        spacingSliceImage[3];

  float* spacing = this->inputImage->spacingImage;
  double start1 = this->croppedImage.imageMatrix.matrix(2,3);
  spacingSliceImage[0] = spacing[0]; spacingSliceImage[1] = spacing[1]; spacingSliceImage[2] = spacing[2];
  this->croppedImage.imageData->SetSpacing( spacing );
  this->SetImageMatrix();

  this->croppedImage.SetParametersFromITK( start1, spacingSliceImage[2], this->croppedImage.imageMatrix );

  /*QuickView viewer;
  viewer.AddImage( this->croppedImage.imageData.GetPointer() ); // Need to do this because QuickView can't accept smart pointers
  viewer.Visualize();*/
  return;

}

void ImageCropping::SetImageMatrix()
{

  float* spacing = this->inputImage->spacingImage;
  this->croppedImage.imageMatrix.matrix = this->inputImage->imageMatrix.matrix;
  int dim[3] = {this->desiredSize[0], this->desiredSize[1], 1 };
  this->croppedImage.imageMatrix.SetDimensionsForIGTMatrix( dim );
  this->croppedImage.imageMatrix.SetSpacingForIGTMatrix( spacing );
  this->croppedImage.imageMatrix.SetIGTTransformFromMat();

  TransformMatrix cropMatrix;
  cropMatrix.SetDimensionsForIGTMatrix( croppedImage.sizeImage );
  cropMatrix.SetSpacingForIGTMatrix( spacing );
  float dStart[3] = { this->desiredStart[0], this->desiredStart[1], 0 };
  cropMatrix.SetOriginInTransform( dStart );

  this->croppedImage.imageMatrix.MultiplyWith( cropMatrix.matrix );

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

  typedef unsigned char PixelType;
  const PixelType defaultValue = 128;

  filter2DTo3D->SetDefaultPixelValue( defaultValue );
  filter2DTo3D->Update();

  typedef itk::ChangeInformationImageFilter< VolumeType > FilterChangeType;
  FilterChangeType::Pointer filterChange = FilterChangeType::New();
  filterChange->SetInput( filter2DTo3D->GetOutput() );

  VolumeType::PointType::VectorType translation;
  translation[0] = this->croppedImage.imageMatrix.matrix(0,3);
  translation[1] = this->croppedImage.imageMatrix.matrix(1,3);
  translation[2] = this->croppedImage.imageMatrix.matrix(2,3);
  VolumeType::PointType origin = filter2DTo3D->GetOutput()->GetOrigin();
  origin = translation;
  filterChange->SetOutputOrigin( origin );
  filterChange->ChangeOriginOn();
  filterChange->Update();

  VolumeType::DirectionType direction;
  direction[0][0] = this->croppedImage.imageMatrix.matrix(0,0);direction[1][0] = this->croppedImage.imageMatrix.matrix(0,1);direction[2][0] = this->croppedImage.imageMatrix.matrix(0,2);
  direction[0][1] = this->croppedImage.imageMatrix.matrix(1,0);direction[1][1] = this->croppedImage.imageMatrix.matrix(1,1);direction[2][1] = this->croppedImage.imageMatrix.matrix(1,2);
  direction[0][2] = this->croppedImage.imageMatrix.matrix(2,0);direction[1][2] = this->croppedImage.imageMatrix.matrix(2,1);direction[2][2] = this->croppedImage.imageMatrix.matrix(2,2);
  filterChange->SetOutputDirection( direction );
  filterChange->ChangeDirectionOn();
  filterChange->Update();

  this->croppedVolume.volumeData = filterChange->GetOutput();

  this->croppedVolume.SetParametersFromITK();

  return;

}
