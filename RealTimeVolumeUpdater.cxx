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

}

