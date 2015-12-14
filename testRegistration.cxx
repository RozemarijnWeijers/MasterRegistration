#include "usrClientIGT.h"
#include "usrVolume.h"
#include "usrImage.h"
#include "usrImageRegistration.h"
#include "usrVolumeReslice.h"

int cropImageVolume( Volume* volume, int dStart[3], int dsize[2], Image* sliceImage )
{

  typedef itk::ExtractImageFilter< VolumeType, ImageType > FilterType;
  // Set parameters for desired image (reslice/ crop)
  VolumeType::IndexType         desiredStart;
  VolumeType::SizeType          desiredSize ;
  desiredStart[0] = dStart[0];  desiredStart[1] = dStart[1];    desiredStart[2] = dStart[2];
  desiredSize[0] = dsize[0];    desiredSize[1] = dsize[1];      desiredSize[2] = 0;

  VolumeType::RegionType        desiredRegion( desiredStart, desiredSize );
  std::cout << "Desired Region: " << desiredRegion << std::endl;

  // Create cropping/reslice
  FilterType::Pointer           filter = FilterType::New();
  filter->SetExtractionRegion( desiredRegion );
  filter->SetInput( volume->volumeData );
  filter->SetDirectionCollapseToIdentity(); // This is required.
  filter->Update();

  // Set resliced image to ITK image
  sliceImage->imageData = filter->GetOutput();

  // Get and set parameters for reslices image    // spacing (mm/pixel)

  VolumeType::PointType         origin;
  double                        startSliceImage[3];
  double                        spacingSliceImage[3];
  double                        originsliceImage[3];
  float                         originsliceImage1[3];
  VolumeType::SizeType          size1;

  float* spacing = volume->spacingVolume;
  float* start1 = volume->originVolume;
  spacingSliceImage[0] = spacing[0]; spacingSliceImage[1] = spacing[1]; spacingSliceImage[2] = spacing[2];
  startSliceImage[0] = start1[0] + dStart[0]; startSliceImage[1] = start1[1] + dStart[1];   startSliceImage[2] = start1[2] + dStart[2];

  sliceImage->imageData->SetSpacing( spacing );
  sliceImage->imageData->SetOrigin( startSliceImage );
  sliceImage->SetParametersFromITK( startSliceImage, spacingSliceImage );

  return 1;

}

int main(int argc, char* argv[])
{
  // Measure running time
  const clock_t begin_time = clock();

  if (argc != 5) // check number of arguments
  {
    // If not correct, print usage
    std::cerr << "    <filenameVolume>    : Filename of Volume.nrrd"            << std::endl;
    std::cerr << "    <hostnameSender1>   : IP or host name"                    << std::endl;
    std::cerr << "    <portSender1>       : Portimage # (18944 default)"        << std::endl;
    std::cerr << "    <sliceNumber>       : Number between 5 and 20 for now"    << std::endl;
    exit(0);
  }

  // Load volume data (NRRD file)
  char*     file = argv[1];
  Volume   volume;
  volume.LoadVolume( file );

  // Establish connections with server
  ClientIGT   client1;
  client1.ConnectToServer( argv[2], atoi( argv[3] ) );

  // Send volume to slicer
  volume.ConvertITKtoIGTVolume();
  client1.imgMsg = volume.imgMsg;
  client1.SendImage();

  // Create images for registration
  Image    fixedImage;
  Image    movingImage;
  Image    registeredImage;

  // Set parameters for testing resliceImage volume, get (part of) a slice from the volume
  int       sliceNumber = atoi(argv[4]);
  int       translation[2];     translation[0] = 0;                 translation[1] = 0; // Optional
  int       dStart[3];          dStart[0] = translation[0];         dStart[1] = translation[1];        dStart[2] =sliceNumber; // slice number "slicenumber" without the 5 most left pixels (so translated to the left)
  VolumeType::SizeType          size = volume.volumeData->GetLargestPossibleRegion().GetSize();
  int       dSize[2];           dSize[0] = size[0]-translation[0];  dSize[1] = size[1]-translation[1]; // Note the switch in axis for the translation

  // Get (part of) a slice of the volume (by cropping) and set it as fixed image for registration
  cropImageVolume( &volume, dStart, dSize, &fixedImage );

  // Send the fixed image to Slicer
  fixedImage.ConvertITKtoIGTImage();
  client1.imgMsg = volume.imgMsg;
  client1.imgMsg->SetDeviceName( "imageSlice" );
  client1.SendImage();

  // Start registration of the resliced images with other resliced images in the area around it (5 slices before and 5 after) to find the best match
  int               testNumber = 11;
  double            values[testNumber];
  double            bestMatch[2];
  ParametersType    finalParameters;
  double            initialMatrix[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  for ( int i = 0; i < testNumber; i++ )
  {
    // Measure running time for one reslice and registration loop
    const clock_t begin_time_reg = clock();

    // Set parameters for the reslice to match with the fixed image
    dStart[0] = 0;  dStart[1] = 0;  dStart[2] = sliceNumber - ((testNumber-1)/2)+i;

    // Get the new reslice of the volume (by cropping) and set it as the moving image for registration
    cropImageVolume( &volume, dStart, dSize, &movingImage );

    // Register the moving and the fixed image and save th emetric values to find the best match
    ImageRegistration registration;
    registration.SetFixedImage( &fixedImage );
    registration.SetMovingImage( &movingImage );
    registration.SetInitialMatrix( initialMatrix );
    registration.RegisterImages();
    std::cerr<<registration.registrationMatrix[5]<<std::endl;
    //values[i] = RegistrationFunction( &fixedImage, &movingImage, &registeredImage, &finalParameters, 0 );
    /*if ( i > 0 )
    {
      if( values[i] < bestMatch[1] )
      {
        bestMatch[0] = i;
        bestMatch[1] = values[i];
        std::cerr << "Best value: " << bestMatch[1] << std::endl;
      }
    }
    else
    {
      bestMatch[0] = i;
      bestMatch[1] = values[i];
      std::cerr << "Best value: " << bestMatch[1] << std::endl;
    }*/
    std::cerr << "Registration " << i+1 << " is done." << std::endl;
    //std::cout << "Registration time: " << float( clock () - begin_time_reg) / CLOCKS_PER_SEC << std::endl;
  }

  // Give the best mathcing slice and the associated metric value and reistration values
  /*dStart[0] = 0;    dStart[1] = 0;    dStart[2] = sliceNumber - ((testNumber-1)/2) + bestMatch[0];
  resliceImageVolume( &volume, dStart, dSize, &movingImage );
  RegisterImage( &fixedImage, &movingImage, &registeredImage, &finalParameters, 1 ); // the 1 activated the image registration function
  std::cerr<< "Slice number " << sliceNumber << " matches best with slice number " << sliceNumber-((testNumber-1)/2)+bestMatch[0] << " with metric: " << bestMatch[1] << std::endl;
  std::cerr<< "Set translation: " << translation[0] << ", " << translation[1] << " vs. found translation by registration: " << finalParameters[0] << ", " << finalParameters[1] << std::endl;
*/
  // Send the unregistered best matching slice and the registered best slice to Slicer
  /*movingImage.ConvertITKtoIGTImage();
  movingImage.imgMsg->SetDeviceName( "bestMatch" );
  movingImage.imgMsg->Pack();
  client1.socket->Send( movingImage.imgMsg->GetPackPointer(), movingImage.imgMsg->GetPackSize() );
  registeredImage.ConvertITKtoIGTImage();
  registeredImage.imgMsg->SetDeviceName( "registeredImage" );
  registeredImage.imgMsg->Pack();
  client1.socket->Send( registeredImage.imgMsg->GetPackPointer(), registeredImage.imgMsg->GetPackSize() );

  // Close connection
  client1.socket->CloseSocket();
*/
  std::cout << "Total running time: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;

}

