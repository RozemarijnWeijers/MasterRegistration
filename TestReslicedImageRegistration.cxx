#include "usrClientIGT.h"
#include "usrVolume.h"
#include "usrImage.h"
#include "usrVolumeRegistration.h"
#include "usrTransformMatrix.h"
#include "usrVolumeReslice.h"
#include "usrVolumeCropping.h"
#include "usrImageCropping.h"
#include "usrRotationMatrix.h"

#include "itkResampleImageFilter.h"
#include "itkTranslationTransform.h"
#include <itkCenteredAffineTransform.h>

int main(int argc, char* argv[]) // Why is this one slow? and why does it stop to recognize image messages after the first ca. 20?
{

  if (argc != 12) // check number of arguments
  {
    // If not correct, print usage
    std::cerr << "    <filenameVolume>    : Filename of Volume.nrrd"            << std::endl;
    std::cerr << "    <hostnameReceiver>  : IP or host name"                    << std::endl;
    std::cerr << "    <portReceiver>      : Port # (18944 default)"             << std::endl;
    std::cerr << "    <SliceNumber1>      : number"                             << std::endl;
    std::cerr << "    <Translation1>      : mm "                                << std::endl;
    std::cerr << "    <Translation2>      : mm "                                << std::endl;
    std::cerr << "    <Translation3>      : mm  "                               << std::endl;
    std::cerr << "    <SliceAngle>        : degree"                             << std::endl;
    std::cerr << "    <SliceAngle>        : degree"                             << std::endl;
    std::cerr << "    <SliceAngle>        : degree"                             << std::endl;
    std::cerr << "    <thickness>         : pixels"                             << std::endl;
    exit(0);
  }

  // Load volume data (NRRD file)
  char*    file = argv[1];
  Volume   volume;
  volume.LoadVolume( file );

  // Establish connections with server
  ClientIGT   client1;
  client1.ConnectToServer( argv[2], atoi( argv[3] ) );

  // Send volume to slicer
  volume.ConvertITKtoIGTVolume();
  client1.imgMsg = volume.imgMsg;
  client1.imgMsg->SetDeviceName( "Volume" );
  client1.SendImage();

  // Create images
  Volume*    fixedImage;
  Volume*    movingImage;

  VolumeReslice resliceVolume;
  resliceVolume.SetVolume( &volume );
  Volume CroppedVolumeSlice;

  //create a reslice
  double angles[3] = { 0, 0, 0};
  if (angles[0] == 0 && angles[1] == 0 && angles[2] == 0)
  {
      float slicenumber1 = 0;
      float slicenumber2 = 0;
      float slicenumber3 = volume.sizeVolume[2]/2;//atoi(argv[4]); //your choise

      float resliceOrigin[3];
      resliceOrigin[0] = slicenumber1;
      resliceOrigin[1] = slicenumber2;
      resliceOrigin[2] = slicenumber3;

      RotationMatrix rotMat;
      rotMat.Set3Angels( angles );

      resliceVolume.SetResliceAxesWRTVolume( rotMat.matrixDouble );
      resliceVolume.SetOriginOfResliceWRTVolume( resliceOrigin );
      resliceVolume.ResliceVolume();
      resliceVolume.CreateITKReslice();

      resliceVolume.reslicedImage.ConvertITKtoIGTImage();
      client1.imgMsg = resliceVolume.reslicedImage.imgMsg;
      client1.imgMsg->SetDeviceName( "ReslicedImage" );
      client1.SendImage();

      //Crop resliced image
      ImageCropping cropReslice;
      cropReslice.SetImage( &resliceVolume.reslicedImage );

      int marge = 20;
      int volumeCropSize[3]; // in pixels
      volumeCropSize[0] = resliceVolume.reslicedImage.sizeImage[0]-marge;
      volumeCropSize[1] = resliceVolume.reslicedImage.sizeImage[1]-marge;
      volumeCropSize[2] = 1;

      float volumeCropOrigin[3];
      float dstart[2] = { marge*resliceVolume.reslicedImage.spacingImage[0]/2 , marge*resliceVolume.reslicedImage.spacingImage[1]/2 };//in mmm
      volumeCropOrigin[0] = resliceOrigin[0] + dstart[0];
      volumeCropOrigin[1] = resliceOrigin[1] + dstart[1];
      volumeCropOrigin[2] = resliceOrigin[2];

      cropReslice.SetCropSizeAndStart( volumeCropSize, volumeCropOrigin );
      cropReslice.CropImage();

      cropReslice.Convert2DImageTo3DVolume();
      CroppedVolumeSlice = cropReslice.croppedVolume;

      cropReslice.croppedVolume.ConvertITKtoIGTVolume();
      cropReslice.croppedImage.ConvertITKtoIGTImage();
      client1.imgMsg = cropReslice.croppedVolume.imgMsg;
      client1.imgMsg->SetDeviceName( "CroppedImage" );
      client1.SendImage();
  }

  TransformType3D::Pointer movetransform = TransformType3D::New();
  typedef TransformType3D::VersorType VersorType;
  typedef VersorType::VectorType VectorType;
  typedef TransformType3D::TranslationType TranslationType;
  VersorType      rotation;
  VectorType      axis;
  TranslationType translation;

  double angle;
  double rad2grad = PI/180;
  if (atoi(argv[8]) != 0 )
  {
      axis[0] = 0.0;
      axis[1] = 0.0;
      axis[2] = 1.0;
      angle = atoi(argv[8]) * rad2grad;
      rotation.Set(axis,angle);
      movetransform->SetRotation(rotation);
  }
  if (atoi(argv[9]) != 0 )
  {
      axis[0] = 0.0;
      axis[1] = 1.0;
      axis[2] = 0.0;
      angle = atoi(argv[9]) * rad2grad;
      rotation.Set(axis,angle);
      movetransform->SetRotation(rotation);
  }
  if (atoi(argv[10]) != 0 )
  {
      axis[0] = 1.0;
      axis[1] = 0.0;
      axis[2] = 0.0;
      angle = atoi(argv[10]) * rad2grad;
      rotation.Set(axis,angle);
      movetransform->SetRotation(rotation);
  }

  translation[0] = atoi(argv[5]);
  translation[1] = atoi(argv[6]);
  translation[2] = atoi(argv[7]);
  movetransform->SetTranslation(translation);

  typedef itk::ResampleImageFilter<VolumeType, VolumeType> ResampleVolumeFilterType;
  ResampleVolumeFilterType::Pointer resampleVolumeFilter = ResampleVolumeFilterType::New();
  resampleVolumeFilter->SetTransform(movetransform.GetPointer());
  resampleVolumeFilter->SetInput(volume.volumeData);

  resampleVolumeFilter->SetOutputSpacing( volume.volumeData->GetSpacing() );
  resampleVolumeFilter->SetOutputOrigin( volume.volumeData->GetOrigin() );
  resampleVolumeFilter->SetOutputDirection( volume.volumeData->GetDirection() );

  VolumeType::SizeType   sizeVol = volume.volumeData->GetLargestPossibleRegion().GetSize();
  resampleVolumeFilter->SetSize( sizeVol );

  resampleVolumeFilter->Update();

  Volume Movedvolume;
  Movedvolume.volumeData = resampleVolumeFilter->GetOutput();
  Movedvolume.SetParametersFromITK();

  Movedvolume.ConvertITKtoIGTVolume();
  client1.imgMsg = Movedvolume.imgMsg;
  client1.imgMsg->SetDeviceName( "MovedVolume" );
  client1.SendImage();

  /*MoveTransformType::InputPointType center;
  center[0] = (volume2.originVolume[0] + (volume2.sizeVolume[0] / 2.0)) * volume2.spacingVolume[0];
  center[1] = (volume2.originVolume[1] + (volume2.sizeVolume[1] / 2.0)) * volume2.spacingVolume[1] ;
  center[2] = (volume2.originVolume[2] + (volume2.sizeVolume[2] / 2.0)) * volume2.spacingVolume[2] ;

  movetransform->SetCenter(center);
  const double degreesToRadians = atan(1.0) / 45.0;

  const double angle1 = atoi(argv[8]) * degreesToRadians;
  MoveTransformType::OutputVectorType axis1;
  axis1[0] = 1;
  axis1[1] = 0;
  axis1[2] = 0;
  const double angle2 = atoi(argv[9]) * degreesToRadians;
  MoveTransformType::OutputVectorType axis2;
  axis2[0] = 0;
  axis2[1] = 1;
  axis2[2] = 0;
  const double angle3 = atoi(argv[10]) * degreesToRadians;
  MoveTransformType::OutputVectorType axis3;
  axis3[0] = 0;
  axis3[1] = 0;
  axis3[2] = 1;

  //movetransform->ComputeOffset();
  movetransform->Rotate3D( axis1, angle1, false );
  movetransform->Rotate3D( axis2, angle2, false );
  movetransform->Rotate3D( axis3, angle3, false );*/

  int marges[3] = { 2*20, 2*20, atoi(argv[11]) }; //in pixels
  int volumeCropSize[3] = {volume.sizeVolume[0] - marges[0], volume.sizeVolume[1] - marges[1], marges[2]};
  float volumeCropOrigin[3] = {volume.originVolume[0] + (marges[0]/2), volume.originVolume[1] + (marges[1]/2), volume.sizeVolume[2]/2 - (marges[2]/2)};

  Volume CroppedVolume;
  volume.CropVolume( volumeCropOrigin, volumeCropSize, &CroppedVolume );

  CroppedVolume.ConvertITKtoIGTVolume();
  client1.imgMsg = CroppedVolume.imgMsg;
  client1.imgMsg->SetDeviceName( "CroppedVolume" );
  client1.SendImage();

  fixedImage = &CroppedVolumeSlice;//&CroppedVolume;//&CroppedVolumeSlice;//
  movingImage = &Movedvolume;//&CroppedVolume;

  VolumeRegistration registration;
  registration.SetFixedVolume( fixedImage );
  registration.SetMovingVolume( movingImage );
  registration.RegisterVolumes();

  registration.registeredVolume.ConvertITKtoIGTVolume();
  client1.imgMsg = registration.registeredVolume.imgMsg;
  client1.imgMsg->SetDeviceName( "VolumeRegistered" );
  client1.SendImage();

}
