#include "usrClientIGT.h"
#include "usrVolume.h"
#include "usrImage.h"
#include "usrVolumeRegistration.h"
#include "usrTransformMatrix.h"
#include "usrVolumeReslice.h"
#include "usrVolumeCropping.h"
#include "usrImageCropping.h"
#include "usrRotationMatrix.h"

int main(int argc, char* argv[]) // Why is this one slow? and why does it stop to recognize image messages after the first ca. 20?
{

  if (argc != 10) // check number of arguments
  {
    // If not correct, print usage
    std::cerr << "    <filenameVolume>    : Filename of Volume.nrrd"            << std::endl;
    std::cerr << "    <hostnameReceiver>  : IP or host name"                    << std::endl;
    std::cerr << "    <portReceiver>      : Port # (18944 default)"             << std::endl;
    std::cerr << "    <SliceAngle1>       : Angle (deg)"                        << std::endl;
    std::cerr << "    <SliceAngle2        : Angle (deg)"                        << std::endl;
    std::cerr << "    <SliceAngle3>       : Angle (deg)"                        << std::endl;
    std::cerr << "    <SliceNumber1>      : number"                             << std::endl;
    std::cerr << "    <SliceNumber2>      : number"                             << std::endl;
    std::cerr << "    <SliceNumber3>      : number"                             << std::endl;
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
  client1.imgMsg->SetDeviceName( "volume" );
  client1.SendImage();

  // Create images
  Image    fixedImage;
  Image    movingImage;

  VolumeReslice resliceVolume;
  resliceVolume.SetVolume( &volume );

  // Reslice Volume for simulation of incoming slice

  float sliceangle1 = atoi(argv[4]);
  float sliceangle2 = atoi(argv[5]);
  float sliceangle3 = atoi(argv[6]);
  RotationMatrix rotMat;
  double angles[3] = { sliceangle1, sliceangle2, sliceangle3};
  rotMat.Set3Angels( angles );
  float slicenumber1 = atoi(argv[7]);
  float slicenumber2 = atoi(argv[8]);
  float slicenumber3 = atoi(argv[9]);
  float resliceOrigin[3] = {slicenumber1, slicenumber2, slicenumber3};
  //std::cerr<< "rotmatrix"<< std::endl;
  //rotMat.ShowMatrix();
  resliceVolume.SetResliceAxesWRTVolume( rotMat.matrixDouble );
  resliceVolume.SetOriginOfResliceWRTVolume( resliceOrigin );
  resliceVolume.ResliceVolume();
  resliceVolume.CreateITKReslice();
  //resliceVolume.reslicedImage.imageMatrix.ShowMatrix();

  resliceVolume.reslicedImage.ConvertITKtoIGTImage();
  client1.imgMsg = resliceVolume.reslicedImage.imgMsg;
  client1.imgMsg->SetDeviceName( "ReslicedImage" );
  client1.SendImage();

  //Crop resliced image
  ImageCropping cropReslice;
  double sizeim[2];
  sizeim[0] = resliceVolume.reslicedImage.sizeImage[0];
  sizeim[1] = resliceVolume.reslicedImage.sizeImage[1];
  cropReslice.SetImage( &resliceVolume.reslicedImage );
  int dsize[3] = {sizeim[0]/2, sizeim[1]/2};
  float dstart[3] = {((sizeim[0]/2)-(dsize[0]/2))*resliceVolume.reslicedImage.spacingImage[0], ((sizeim[1]/2)-(dsize[1]/2))*resliceVolume.reslicedImage.spacingImage[1]};
  cropReslice.SetCropSizeAndStart( dsize, dstart );
  cropReslice.CropImage();
  cropReslice.Convert2DImageTo3DVolume();

  cropReslice.croppedVolume.ConvertITKtoIGTVolume();
  cropReslice.croppedImage.ConvertITKtoIGTImage();
  client1.imgMsg = cropReslice.croppedVolume.imgMsg;
  client1.imgMsg->SetDeviceName( "CroppedImage" );
  client1.SendImage();
  // Create a message buffer to receive header and image message
  /*igtl::MessageHeader::Pointer headerMsg = igtl::MessageHeader::New();

  // Allocate a time stamp, WHERE IS THIS USED?????????
  igtl::TimeStamp::Pointer ts = igtl::TimeStamp::New();
  // Receive first image
  ReceiveImage( &clientIGT1, &fixedImage, ts, headerMsg, fixedImage.imgMsg );*/


  /*clientIGT1.ReceiveImage();
  fixedImage.imgMsg = clientIGT1.imgMsg;
  fixedImage.SetParametersFromIGT();
  fixedImage.ConvertIGTtoITKImage();

  ImageRegistration registration;
  double initialMatrix[9];

  while ( 1 )
  {
    for ( int i = 0; i < 100; i ++ ) //WHY 100?
    {
      if ( clientIGT1.ReceiveImage() == 0 )
      {
        movingImage.imgMsg = clientIGT1.imgMsg;
        movingImage.SetParametersFromIGT();
        movingImage.ConvertIGTtoITKImage();

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

        registration.SetFixedImage( &fixedImage );
        registration.SetMovingImage( &movingImage );
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
