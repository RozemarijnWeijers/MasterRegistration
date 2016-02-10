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

  if (argc != 6) // check number of arguments
  {
    // If not correct, print usage
    std::cerr << "    <filenameVolume>    : Filename of Volume.nrrd"            << std::endl;
    std::cerr << "    <hostnameReceiver>  : IP or host name"                    << std::endl;
    std::cerr << "    <portReceiver>      : Port # (18944 default)"             << std::endl;
    std::cerr << "    <hostnameSender1>   : IP or host name"                    << std::endl;
    std::cerr << "    <portSender1>       : Portimage # (18944 default)"        << std::endl;
    exit(0);
  }

  // Load volume data
  char*    file = argv[1];
  Volume   movedVolume;
  movedVolume.LoadVolume( file );


  // Establish connections
  ClientIGT clientIGT1;
  clientIGT1.ConnectToServer( argv[2],atoi( argv[3] ) );
  ClientIGT clientIGT2;
  clientIGT2.ConnectToServer( argv[4],atoi( argv[5] ) );


  // Send ITK volume to Slicer
  movedVolume.ConvertITKtoIGTVolume();
  clientIGT2.imgMsg = movedVolume.imgMsg;
  clientIGT2.imgMsg->SetDeviceName( "MovedVolume" );
  clientIGT2.SendImage();

  // Create images
  Volume*    fixedImage;
  Volume*    movingImage;
  Image*     sliceImage1;
  Image*     sliceImage2;
  Volume     newVolume;
  Image      newImage;

  movingImage = &movedVolume;
  fixedImage = &newVolume;
  sliceImage1 = &newImage;

  clientIGT1.ReceiveImage();

  fixedImage->imgMsg = clientIGT1.imgMsg;
  fixedImage->SetParametersFromIGT();
  fixedImage->ConvertIGTtoITKVolume();

  sliceImage1->imgMsg = clientIGT1.imgMsg;
  sliceImage1->SetParametersFromIGT();
  sliceImage1->ConvertIGTtoITKImage();

  fixedImage->ConvertITKtoIGTVolume();
  clientIGT2.imgMsg = fixedImage->imgMsg;
  clientIGT2.imgMsg->SetDeviceName( "FixedImage" );
  clientIGT2.SendImage();

  VolumeRegistration registration;
  registration.SetFixedVolume( fixedImage );
  registration.SetMovingVolume( movingImage );
  registration.RegisterVolumes();

  registration.registeredVolume.ConvertITKtoIGTVolume();
  clientIGT2.imgMsg = registration.registeredVolume.imgMsg;
  clientIGT2.imgMsg->SetDeviceName( "VolumeRegistered" );
  clientIGT2.SendImage();


  /*ImageRegistration registration;
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
        //registeredImage.SetParametersFromITK();

        // Calculate subtraction of registered and fixed image*/
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

