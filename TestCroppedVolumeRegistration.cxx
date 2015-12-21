#include "usrClientIGT.h"
#include "usrVolume.h"
#include "usrImage.h"
#include "usrImageRegistration.h"
#include "usrVolumeReslice.h"
#include "usrVolumeRegistration.h"

void cropVolume( Volume* volume, int dStart[3], int dsize[3], Volume* sliceVolume )
{

  typedef itk::ExtractImageFilter< VolumeType, VolumeType > FilterType;
  // Set parameters for desired image (reslice/ crop)
  VolumeType::IndexType         desiredStart;
  VolumeType::SizeType          desiredSize ;
  desiredStart[0] = dStart[0];  desiredStart[1] = dStart[1];    desiredStart[2] = dStart[2];
  desiredSize[0] = dsize[0];    desiredSize[1] = dsize[1];      desiredSize[2] = dsize[2];

  VolumeType::RegionType        desiredRegion( desiredStart, desiredSize );
  std::cout << "Desired Region: " << desiredRegion << std::endl;

  // Create cropping/reslice
  FilterType::Pointer           filter = FilterType::New();
  filter->SetExtractionRegion( desiredRegion );
  filter->SetInput( volume->volumeData );
  filter->SetDirectionCollapseToIdentity(); // This is required.
  filter->Update();

  // Set resliced image to ITK image
  sliceVolume->volumeData = filter->GetOutput();

  // Get and set parameters for reslices image    // spacing (mm/pixel)

  VolumeType::PointType         origin;
  double                        originsliceVolume[3];
  double                        spacingsliceVolume[3];

  float* spacing = volume->spacingVolume;
  float* start1 = volume->originVolume;
  spacingsliceVolume[0] = spacing[0]; spacingsliceVolume[1] = spacing[1]; spacingsliceVolume[2] = spacing[2];
  originsliceVolume[0] = start1[0] + dStart[0]; originsliceVolume[1] = start1[1] + dStart[1];   originsliceVolume[2] = start1[2] + dStart[2];

  sliceVolume->volumeData->SetSpacing( spacing );
  sliceVolume->volumeData->SetOrigin( originsliceVolume );
  sliceVolume->SetParametersFromITK();

  return;

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
  client1.imgMsg->SetDeviceName( "volume" );
  client1.SendImage();

  // Create images for registration
  Volume    fixedVolume;
  Volume    movingVolume;

  // Set parameters for testing resliceVolume volume, get (part of) a slice from the volume
  int       sliceNumber = atoi(argv[4]);
  int       translation[2];     translation[0] = 0;                 translation[1] = 0; // Optional
  int       dStart[3];          dStart[0] = translation[0];         dStart[1] = translation[1];        dStart[2] =sliceNumber; // slice number "slicenumber" without the 5 most left pixels (so translated to the left)
  VolumeType::SizeType          size = volume.volumeData->GetLargestPossibleRegion().GetSize();
  int       dSize[3];           dSize[0] = size[0]-translation[0];  dSize[1] = size[1]-translation[1];  dSize[2] = 1;// Note the switch in axis for the translation

  // Get (part of) a slice of the volume (by cropping) and set it as fixed image for registration
  cropVolume( &volume, dStart, dSize, &fixedVolume );

  // Send the fixed image to Slicer
  fixedVolume.ConvertITKtoIGTVolume();
  float spac[3];
  fixedVolume.imgMsg->GetSpacing(spac);
  std::cerr<< spac[0]<<std::endl;
  client1.imgMsg = fixedVolume.imgMsg;
  client1.imgMsg->SetDeviceName( "imageSlice" );
  client1.SendImage();

  // Start registration of the resliced images with other resliced images in the area around it (5 slices before and 5 after) to find the best match
  int               testNumber = 11;
  double            values[testNumber];
  double            bestMatch[2];
  ParametersType    finalParameters;
  double            initialMatrix[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};

  VolumeRegistration registration;
  registration.SetFixedVolume( &fixedVolume );

  for ( int i = 0; i < testNumber; i++ )
  {
    // Measure running time for one reslice and registration loop
    const clock_t begin_time_reg = clock();

    // Set parameters for the reslice to match with the fixed image
    dStart[0] = 0;  dStart[1] = 0;  dStart[2] = sliceNumber - ((testNumber-1)/2)+i;
    int       dSize[3];           dSize[0] = size[0]-translation[0];  dSize[1] = size[1]-translation[1];  dSize[2] = 6;

    // Get the new reslice of the volume (by cropping) and set it as the moving image for registration
    cropVolume( &volume, dStart, dSize, &movingVolume );

    // Register the moving and the fixed image and save th emetric values to find the best match
    registration.SetMovingVolume( &movingVolume );
    registration.SetInitialMatrix( initialMatrix );
    registration.RegisterVolumes();
    //registration.CreateRegisteredImage();
    values[i] = registration.metricValue;

    if ( i > 0 )
    {
      if( values[i] < bestMatch[1] )
      {
        bestMatch[0] = i;
        bestMatch[1] = values[i];
        std::cerr << "New best value: " << bestMatch[1] << std::endl;
      }
    }
    else
    {
      bestMatch[0] = i;
      bestMatch[1] = values[i];
      std::cerr << "Best value: " << bestMatch[1] << std::endl;
    }
    std::cerr << "Registration " << i+1 << " is done." << std::endl;
    //std::cout << "Registration time: " << float( clock () - begin_time_reg) / CLOCKS_PER_SEC << std::endl;
  }

  // Give the best mathcing slice and the associated metric value and reistration values
  /*dStart[0] = 0;    dStart[1] = 0;    dStart[2] = sliceNumber - ((testNumber-1)/2) + bestMatch[0];
  cropVolume( &volume, dStart, dSize, &movingVolume );
  registration.SetmovingVolume( &movingVolume );
  registration.SetInitialMatrix( initialMatrix );
  registration.RegisterImages();
  registration.CreateRegisteredImage();
  std::cerr<< "Slice number " << sliceNumber << " matches best with slice number " << sliceNumber-((testNumber-1)/2)+bestMatch[0] << " with metric: " << bestMatch[1] << std::endl;
  //std::cerr<< "Set translation: " << translation[0] << ", " << translation[1] << " vs. found translation by registration: " << finalParameters[0] << ", " << finalParameters[1] << std::endl;

  //Send the unregistered best matching slice and the registered best slice to Slicer
  movingVolume.ConvertITKtoIGTImage();
  client1.imgMsg = movingVolume.imgMsg;
  client1.imgMsg->SetDeviceName( "bestMatch" );
  client1.SendImage();
  registration.registeredImage.ConvertITKtoIGTImage();
  client1.imgMsg = registration.registeredImage.imgMsg;
  client1.imgMsg->SetDeviceName( "registedMatch" );
  client1.SendImage();*/

  // Close connection
  client1.socket->CloseSocket();

  std::cout << "Total running time: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;

}

