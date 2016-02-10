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

  if (argc != 5) // check number of arguments
  {
    // If not correct, print usage
    std::cerr << "    <filenameVolume>    : Filename of Volume.nrrd"            << std::endl;
    std::cerr << "    <filenameImage>     : Filename of Volume.nrrd"            << std::endl;
    std::cerr << "    <hostnameReceiver>  : IP or host name"                    << std::endl;
    std::cerr << "    <portReceiver>      : Port # (18944 default)"             << std::endl;
    exit(0);
  }

  // Load volume data (NRRD file)
  char*    file = argv[1];
  Volume   movedVolume;
  movedVolume.LoadVolume( file );

  char*     file2 = argv[2];
  Volume    imageSlice;
  imageSlice.LoadVolume( file2 );

  // Establish connections with server
  ClientIGT   client1;
  client1.ConnectToServer( argv[3], atoi( argv[4] ) );

  // Send volume to slicer
  movedVolume.ConvertITKtoIGTVolume();
  client1.imgMsg = movedVolume.imgMsg;
  client1.imgMsg->SetDeviceName( "MovedVolume" );
  client1.SendImage();

  imageSlice.ConvertITKtoIGTVolume();
  client1.imgMsg = imageSlice.imgMsg;
  client1.imgMsg->SetDeviceName( "ImageSlice" );
  client1.SendImage();

  // Create images
  Volume*    fixedImage;
  Volume*    movingImage;

  fixedImage = &imageSlice;
  movingImage = &movedVolume;

  VolumeRegistration registration;
  registration.SetFixedVolume( fixedImage );
  registration.SetMovingVolume( movingImage );
  registration.RegisterVolumes();

  registration.registeredVolume.ConvertITKtoIGTVolume();
  client1.imgMsg = registration.registeredVolume.imgMsg;
  client1.imgMsg->SetDeviceName( "VolumeRegistered" );
  client1.SendImage();

  }

