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

  double angles[3] = { 0, 0, 0};
  if (angles[0] == 0 && angles[1] == 0 && angles[2] == 0)
  {
      float slicenumber1 = 0;
      float slicenumber2 = 0;
      float slicenumber3 = volume.sizeVolume[2]/2;//atoi(argv[4]);

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
  /*if (sliceangle3 == 90 )
  {
      int s[3] = {(volume.sizeVolume[0]-50), 1, (volume.sizeVolume[2]-25)};
      volumeCropSize[0] = s[0]; volumeCropSize[1] = s[1]; volumeCropSize[2] = s[2];
      volumeCropOrigin[0] = volume.originVolume[0] + 20; volumeCropOrigin[1] = atoi(argv[7]); volumeCropOrigin[2] =volume.originVolume[0] + 10;
      //int volumeCropSize[3] = { 160, 40, 65};//cropReslice.croppedVolume.sizeVolume[0]+marge, cropReslice.croppedVolume.sizeVolume[2]+marge, cropReslice.croppedVolume.sizeVolume[1]+marge};
      //float volumeCropOrigin[3] = {20, 50,5};// (resliceOrigin[0]+dstart[0])/volume.spacingVolume[0]-(marge/2), (resliceOrigin[2])/volume.spacingVolume[2]-(marge/2), (resliceOrigin[1]+dstart[1])/volume.spacingVolume[1]-(marge/2)};

      volume.CropVolume( volumeCropOrigin, volumeCropSize, &CroppedVolumeSlice );

      CroppedVolumeSlice.ConvertITKtoIGTVolume();
      client1.imgMsg = CroppedVolumeSlice.imgMsg;
      client1.imgMsg->SetDeviceName( "CroppedImage" );
      client1.SendImage();
  }
  if (sliceangle2 == -90 )
  {
      int s[3] = { 1,volume.sizeVolume[1]-50, volume.sizeVolume[2]-50 };
      volumeCropSize[0] = s[0]; volumeCropSize[1] = s[1]; volumeCropSize[2] = s[2];
      volumeCropOrigin[0] = atoi(argv[7]); volumeCropOrigin[1] = volume.originVolume[1] + 20; volumeCropOrigin[2] =volume.originVolume[0] + 20;

      //int volumeCropSize[3] = { 160, 40, 65};//cropReslice.croppedVolume.sizeVolume[0]+marge, cropReslice.croppedVolume.sizeVolume[2]+marge, cropReslice.croppedVolume.sizeVolume[1]+marge};
      //float volumeCropOrigin[3] = {20, 50,5};// (resliceOrigin[0]+dstart[0])/volume.spacingVolume[0]-(marge/2), (resliceOrigin[2])/volume.spacingVolume[2]-(marge/2), (resliceOrigin[1]+dstart[1])/volume.spacingVolume[1]-(marge/2)};

      volume.CropVolume( volumeCropOrigin, volumeCropSize, &CroppedVolumeSlice );

      CroppedVolumeSlice.ConvertITKtoIGTVolume();
      client1.imgMsg = CroppedVolumeSlice.imgMsg;
      client1.imgMsg->SetDeviceName( "CroppedImage" );
      client1.SendImage();
  }*/


  /*typedef itk::TranslationTransform<double,3> TranslationTransformType;
  TranslationTransformType::Pointer translationtransform = TranslationTransformType::New();
  TranslationTransformType::OutputVectorType tTranslation;
  tTranslation[0] = atoi(argv[4]); // in mm?
  tTranslation[1] = atoi(argv[5]); // in mm?
  tTranslation[2] = atoi(argv[6]); // in mm?
  translationtransform->Translate(tTranslation);*/

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

  /*double initialMatrix[9];

  while ( 1 )
  {
    for ( int i = 0; i < 100; i ++ ) //WHY 100?
    {
      if ( clientIGT1.ReceiveImage() == 0 )
      {
        // Make sure to have two different images
        if ( fixedImage.imgMsg == movingImage.imgMsg )
        {
          std::cerr<<"Fixed and moving are the same"<<std::endl;
        }
        else
        {
        // Register fixed and moving image (2D/2D Registration of two subsequential images), two sequential received images for now, but this wil later be the received image with the reslice of the volume at that position
        initialMatrix[0] = 1; initialMatrix[1] = 0; initialMatrix[2] = 0;
        initialMatrix[3] = 0; initialMatrix[4] = 1; initialMatrix[5] = 0;
        initialMatrix[6] = 0; initialMatrix[7] = 0; initialMatrix[8] = 1;


        //registration.SetInitialMatrix( initialMatrix );
        registration.RegisterImages();
        //registrationRegistrationFunction( &fixedImage, &movingImage, &registeredImage );

        // Create imageMessage for registered image
        //registeredImage = registration.registeredImage;
        //registeredImage.SetParametersFromITK();
        //registeredImage.ITKtoIGTImage();
        //registeredImage.SetParametersFromITK();*/

        // Calculate subtraction of registered and fixed image
        /*typedef itk::SubtractImageFilter <ImageType, ImageType > SubtractImageFilterType;
        SubtractImageFilterType::Pointer subtractFilter = SubtractImageFilterType::New ();
        subtractFilter->SetInput1(fixedImage.imageData);
        subtractFilter->SetInput2(registeredImage.imageData);
        subtractFilter->Update();

        typedef itk::SubtractImageFilter <ImageType, ImageType > SubtractImageFilterType2;
        SubtractImageFilterType2::Pointer subtractFilter2 = SubtractImageFilterType2::New ();
        subtractFilter2->SetInput1(fixedImage.imageData);
        subtractFilter2->SetInput2(movingImage.imageData);
        subtractFilter2->Update();

        Images subtract1;
        subtract1.imageData = subtractFilter->GetOutput();
        subtract1.ITKtoIGTImage(temspacing);
        subtract1.imgMsg->SetDeviceName("FixedMinusRegistered");
        subtract1.imgMsg->Pack();

        Images subtract2;
        subtract2.imageData = subtractFilter2->GetOutput();
        subtract2.ITKtoIGTImage(temspacing);
        subtract2.imgMsg->SetDeviceName("FixedMinusMoving");
        subtract2.imgMsg->Pack();*/

        // Set massages names
        /*clientIGT3.imgMsg = fixedImage.imgMsg;
        clientIGT3.imgMsg->SetDeviceName( "fixedImage" );
        clientIGT3.SendImage();
        //fixedImage.imgMsg->Pack();
        clientIGT3.imgMsg = movingImage.imgMsg;
        clientIGT3.imgMsg->SetDeviceName( "movingImage" );
        clientIGT3.SendImage();
        //movingImage.imgMsg->Pack();
        //registeredImage.imgMsg->SetDeviceName( "registeredImage" );
        //clientIGT3.SendImage( registeredImage.imgMsg );
        //registeredImage.imgMsg->Pack();

        // Send image messages to sliceImager for visualization
        //clientIGT3.socket->Send( fixedImage.imgMsg->GetPackPointer(), fixedImage.imgMsg->GetPackSize() );
        //clientIGT3.socket->Send( movingImage.imgMsg->GetPackPointer(), movingImage.imgMsg->GetPackSize() );
        //clientIGT3.socket->Send( registeredImage.imgMsg->GetPackPointer(), registeredImage.imgMsg->GetPackSize() );
        //clientIGT3.socket->Send(subtract1.imgMsg->GetPackPointer(), subtract1.imgMsg->GetPackSize());
        //clientIGT3.socket->Send(subtract2.imgMsg->GetPackPointer(), subtract2.imgMsg->GetPackSize());

        // Set moving image to fixed image so it can be registered to the next image received
        //Images  temp = movingImage;
        fixedImage = movingImage;
        //fixedImage = temp;
        }
       }
    }
    std::cerr<< "Stopped receiving messages" <<std::endl;
  }

  // Close connection (The example code never reaches this section ...)
  clientIGT1.socket->CloseSocket();
  clientIGT2.socket->CloseSocket();
  clientIGT3.socket->CloseSocket();*/
}
