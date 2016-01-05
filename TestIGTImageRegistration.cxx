#include "usrClientIGT.h"
#include "usrVolume.h"
#include "usrImage.h"
#include "usrImageRegistration.h"
#include "usrVolumeReslice.h"

/*#include "itkImage.h"
#include "itkImageRegistrationMethod.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkAffineTransform.h"
#include "itkSubtractImageFilter.h"
#include "itkTranslationTransform.h"
#include <itkMatrix.h>
#include <itkNrrdImageIO.h>
#include <itkExtractImageFilter.h>

#include "itkVTKImageToImageFilter.h"

#include <vtkNrrdReader.h>
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkMatrix4x4.h"
#include "vtkImageReslice.h"
#include "vtkLookupTable.h"
#include "vtkImageMapToColors.h"

#include "igtlOSUtil.h"
#include "igtlMessageHeader.h"
#include "igtlTransformMessage.h"
#include "igtlPositionMessage.h"
#include "igtlImageMessage.h"
#include "igtlClientSocket.h"
#include "igtlStatusMessage.h"
#include <igtl_util.h>

typedef  unsigned char   PixelType;
typedef  itk::Image< PixelType, 2 >  ImageType;
typedef  itk::Image< PixelType, 3 >  VolumeType;

//  The transform that will map the fixed image into the moving image.
typedef itk::AffineTransform< double, 2 > TransformType;
//  An optimizer is required to explore the parameter space of the transform in search of optimal values of the metric.
typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
//  The metric will compare how well the two images match each other. Metric types are usually parameterized by the image types as it can be seen in the following type declaration.
typedef itk::MeanSquaresImageToImageMetric< ImageType, ImageType > MetricType;
//  Finally, the type of the interpolator is declared. The interpolator will evaluate the intensities of the moving image at non-grid positions.
typedef itk::LinearInterpolateImageFunction< ImageType, double > InterpolatorType;
//  The registration method type is instantiated using the types of the fixed and moving images. This class is responsible for interconnecting all the components that we have described so far.
typedef itk::ImageRegistrationMethod< ImageType, ImageType >    RegistrationType;

typedef itk::ImageFileReader<VolumeType> FileReaderType;
typedef itk::ExtractImageFilter< VolumeType, ImageType > FilterType;
typedef itk::ResampleImageFilter< ImageType, ImageType > ResampleFilterType;
typedef itk::VTKImageToImageFilter<ImageType> VTKImageToImageType;
typedef RegistrationType::ParametersType ParametersType;*/

/*class Images
{

  public:

  Images();
  int IGTtoITKImage();
  int ITKtoIGTImage();
  void SetParametersFromIGT();
  void SetParametersFromITK( double[3], double[3] );//ImageType::PointType, ImageType::SpacingType );

  ImageType::Pointer imageData;
  igtl::ImageMessage::Pointer imgMsg;

  protected:

  int sizeIm[3];
  float originIm[3];
  float spacingIm[3];
  igtl::Matrix4x4   matrixIm;   // image origin and orientation matrix

};

Images::Images()
{

  // Create imageMessage and ITKimage for this image
  imageData  = ImageType::New();
  imgMsg = igtl::ImageMessage::New();

}

void Images::SetParametersFromIGT()
{

  this->imgMsg->GetDimensions( this->sizeIm );
  this->imgMsg->GetOrigin( this->originIm );
  this->imgMsg->GetSpacing( this->spacingIm );
  this->imgMsg->GetMatrix( this->matrixIm );

  return;

}

void Images::SetParametersFromITK( double origin[3], double spacing[3] )
{

  ImageType::SizeType               size = this->imageData->GetLargestPossibleRegion().GetSize();
  this->sizeIm[0] = size[0];        this->sizeIm[1] = size[1];          this->sizeIm[2] = 1;
  this->originIm[0] = origin[0];    this->originIm[1] = origin[1];      this->originIm[2] = origin[2];
  this->spacingIm[0] = spacing[0];  this->spacingIm[1] = spacing[1];    this->spacingIm[2] = spacing[2]; //laatste is 0

  return;

}

int Images::IGTtoITKImage()
{

  // Retrieve the image data from image message
  //int       size[3];          // image dimension (pixels)
  //float     spacing[3];       // spacing (mm/pixel)
  float     spacingITK[0];    // spacing for 2D image (mm/pixel)
  //int       endian;           // endian (not used)
  //float     origin[3];        // origin ()
  //igtl::Matrix4x4   matrix;   // image origin and orientation matrix

  //endian = this->imgMsg->GetEndian();
  //this->imgMsg->GetDimensions( size );
  //this->imgMsg->GetSpacing( spacing );
  //this->imgMsg->GetMatrix( matrix );
  //this->imgMsg->GetOrigin( origin );

  // Set image data to ITK image (3D information -> 2D)
  ImageType::RegionType     region;
  ImageType::IndexType      start;
  start[0] = this->originIm[0];     start[1] = this->originIm[1];//start[0] = origin[0];     start[1] = origin[1];
  ImageType::SizeType       sizeregion;
  sizeregion[0] = this->sizeIm[0];  sizeregion[1] = sizeIm[1];//sizeregion[0] = size[0];  sizeregion[1] = size[1];
  region.SetSize( sizeregion );
  region.SetIndex( start );
  spacingITK[0] = this->spacingIm[0]; spacingITK[1] = this->spacingIm[1];//spacingITK[0] = spacing[0]; spacingITK[1] = spacing[1];

  this->imageData->SetRegions( region );
  this->imageData->SetSpacing( spacingITK );
  this->imageData->Allocate();

  // Copy image data into ITK image
  memcpy( this->imageData->GetBufferPointer(), this->imgMsg->GetScalarPointer(), this->imgMsg->GetSubVolumeImageSize() );
  this->imageData->Modified();

  return 1;

}

int Images::ITKtoIGTImage()
{

  // Retrieve the image data from ITK image
  ImageType::SizeType       size;
  //ImageType::IndexType      start;
  //ImageType::SpacingType    spacing;
  //ImageType::PointType      origin;

  size = this->imageData->GetLargestPossibleRegion().GetSize();
  //start = this->imageData->GetLargestPossibleRegion().GetIndex();
  //spacing = this->imageData->GetSpacing();
  //origin = this->imageData->GetOrigin();

  // Set image data to image message (2D information -> 3D)
  int                       sizeIGT[3];
  ImageType::RegionType     region;
  float                     originIGT[3];
  float                     spacingIGT[3];  // spacing (mm/pixel)
  int                       scalarType;     // always UINT8
  //int                     svoffset[3];    // sub-volume offset
  //int                     svsize[3];      // sub-volume size
  //sizeIGT[0] = this->sizeIm[0];         sizeIGT[1] = this->sizeIm[1];           sizeIGT[2] = 1;//sizeIGT[0] = size[0];         sizeIGT[1] = size[1];           sizeIGT[2] = 1;
  sizeIGT[0] = size[0];         sizeIGT[1] = size[1];           sizeIGT[2] = 1; //WAAROM??
  spacingIGT[0] = this->spacingIm[0];   spacingIGT[1] = this->spacingIm[1];     spacingIGT[2] = this->spacingIm[2];//spacing[0];   spacingIGT[1] = spacing[1];     spacingIGT[2] = spacing[1]; // third dimension?
  //spacingIGT[0] = spacing[0];   spacingIGT[1] = spacing[1];     spacingIGT[2] = spacing[1]; // third dimension?
  originIGT[0] = this->originIm[0]+((this->sizeIm[0]-1)*this->spacingIm[0]/2);          originIGT[1] = this->originIm[1]+((this->sizeIm[1]-1)*this->spacingIm[1]/2);     originIGT[2] = this->originIm[2];//+((size[2]-1)*spacing[2]/2); // klopt nog niet
  //originIGT[0] = origin[0]+((size[0]-1)*spacing[0]/2);          originIGT[1] = origin[1]+((size[1]-1)*spacing[1]/2);     originIGT[2] = origin[2];//+((size[2]-1)*spacing[2]/2); // klopt nog niet
  //svsize[0] = sizeIGT[0];     svsize[1] = sizeIGT[1];         svsize[2] = sizeIGT[2];
  //svoffset[0] = start[0];     svoffset[1] = start[1];         svoffset[2] = start[2];
  scalarType = igtl::ImageMessage::TYPE_UINT8;

  this->imgMsg->SetDimensions( sizeIGT );
  this->imgMsg->SetSpacing( spacingIGT );
  this->imgMsg->SetScalarType( scalarType );
  this->imgMsg->SetOrigin( originIGT );
  //this->imgMsg->SetSubVolume( sizeIGT, svoffset );
  this->imgMsg->AllocateScalars();

  // Copy image data into ITK image
  memcpy( this->imgMsg->GetScalarPointer(), imageData->GetBufferPointer(), this->imgMsg->GetSubVolumeImageSize() );
  // Pack (serialize) and send
  this->imgMsg->Pack();

  return 1;

}
*/

/*class Volumes
{

  public:

  Volumes();
  int ITKtoIGTVolume();
  int VTKtoITKVolume();
  void SetParametersFromITK();
  void SetParametersFromVTK( double[3], double[3], int[3]);

  VolumeType::Pointer volumeData;
  vtkSmartPointer<vtkNrrdReader> VTKreader;
  igtl::ImageMessage::Pointer imgMsg;

  protected:

  int sizeVol[3];
  float originVol[3];
  float spacingVol[3];

};*/

/*Volumes::Volumes()
{

  // Create imageMessage and ITKvolume for this volume
  volumeData  = VolumeType::New();
  VTKreader = vtkSmartPointer<vtkNrrdReader>::New();
  imgMsg = igtl::ImageMessage::New();

}*/

/*void Volumes::SetParametersFromITK()
{

  VolumeType::SizeType              size = this->volumeData->GetLargestPossibleRegion().GetSize();
  this->sizeVol[0] = size[0];       this->sizeVol[1] = size[1];        this->sizeVol[2] = size[2];
  VolumeType::PointType             origin = this->volumeData->GetOrigin();
  this->originVol[0] = origin[0];   this->originVol[1] = origin[1];    this->originVol[2] = origin[2];
  VolumeType::SpacingType           spacing = this->volumeData->GetSpacing();
  this->spacingVol[0] = spacing[0]; this->spacingVol[1] = spacing[1];  this->spacingVol[2] = spacing[2]; //laatste is 0

  return;

}*/

/*void Volumes::SetParametersFromVTK( double origin[3], double spacing[3], int size[3])
{

  this->sizeVol[0] = size[0];       this->sizeVol[1] = size[1];        this->sizeVol[2] = size[2];
  this->originVol[0] = origin[0];   this->originVol[1] = origin[1];    this->originVol[2] = origin[2];
  this->spacingVol[0] = spacing[0]; this->spacingVol[1] = spacing[1];  this->spacingVol[2] = spacing[2]; //laatste is 0

  return;

}*/

/*int Volumes::ITKtoIGTVolume()
{

  // Retrieve the image data from ITK image
  //VolumeType::SizeType size;
  //VolumeType::IndexType start;
  //VolumeType::PointType origin;
  //VolumeType::SpacingType spacing;

  //size = this->volumeData->GetLargestPossibleRegion().GetSize();
  //start = this->volumeData->GetLargestPossibleRegion().GetIndex();
  //spacing = this->volumeData->GetSpacing();
  //origin = this->volumeData->GetOrigin();

  // Set volume data to image message
  int   sizeIGT[3];
  int   scalarType;
  float originIGT[3];
  float spacingIGT[3];    // spacing (mm/pixel)
  //int   svsize[3];      // sub-volume size
  //int   svoffset[3];    // sub-volume offset
  sizeIGT[0] = this->sizeVol[0];         sizeIGT[1] = this->sizeVol[1];           sizeIGT[2] = this->sizeVol[2];//sizeIGT[0] = size[0];         sizeIGT[1] = size[1];           sizeIGT[2] = 1;
  spacingIGT[0] = this->spacingVol[0];   spacingIGT[1] = this->spacingVol[1];     spacingIGT[2] = this->spacingVol[2];//spacing[0];   spacingIGT[1] = spacing[1];     spacingIGT[2] = spacing[1]; // third dimension?
  originIGT[0] = this->originVol[0]+((this->sizeVol[0]-1)*this->spacingVol[0]/2);          originIGT[1] = this->originVol[1]+((this->sizeVol[1]-1)*this->spacingVol[1]/2);     originIGT[2] = this->originVol[2]+((this->sizeVol[2]-1)*this->spacingVol[2]/2); // klopt nog niet
  //svsize[0] = sizeIGT[0];     svsize[1] = sizeIGT[1];     svsize[2] = sizeIGT[2];
  //svoffset[0] = start[0];     svoffset[1] = start[1];     svoffset[2] = start[2];
  scalarType = igtl::ImageMessage::TYPE_UINT8;

  this->imgMsg->SetDimensions( sizeIGT );
  this->imgMsg->SetSpacing( spacingIGT );
  this->imgMsg->SetScalarType( scalarType );
  this->imgMsg->SetOrigin( originIGT );
  this->imgMsg->AllocateScalars();

  // Set copy image data into ITK image
  memcpy( this->imgMsg->GetScalarPointer(), volumeData->GetBufferPointer(), this->imgMsg->GetImageSize() );
  // Pack (serialize) and send
  this->imgMsg->Pack();

  return 1;

}*/

/*class ClientIGTs
{

  public:

  ClientIGTs( char*, int );

  igtl::ClientIGTSocket::Pointer socket;

};

ClientIGTs::ClientIGTs( char* host, int port )
{

  //Open a socket for th ClientIGT
  socket = igtl::ClientIGTSocket::New();

  //Connect to the server
  int r = socket->ConnectToServer( host, port );

  // Check the connection
  if ((r) != 0)
  {
    std::cerr << "Cannot connect to the server." << std::endl;
    exit(0);
  }

  std::cerr << "ClientIGT is connected to server" << host <<":"<< port<< std::endl;

}*/

/*int ReceiveImage( ClientIGTs* ClientIGT, Images* image, igtl::TimeStamp::Pointer ts, igtl::MessageHeader::Pointer headerMsg, igtl::ImageMessage::Pointer imgMsg )
{

  // Initialize receive buffer
  headerMsg->InitPack();

  // Receive header message
  if (0 == ClientIGT->socket->Receive(headerMsg->GetPackPointer(), headerMsg->GetPackSize()))
  {
    ClientIGT->socket->CloseSocket();
    std::cerr<< " Exit 0 " <<std::endl;
    return 0;
  }

  // Unpack and deserialize the header
  headerMsg->Unpack();

  // Get time stamp
  igtlUint32 sec;
  igtlUint32 nanosec;
  headerMsg->GetTimeStamp(ts);
  ts->GetTimeStamp(&sec, &nanosec);

  if (strcmp(headerMsg->GetDeviceType(), "IMAGE") == 0)
  {
    // Allocate memory for a message buffer to receive transform data
    imgMsg->SetMessageHeader(headerMsg);
    imgMsg->AllocatePack();
    // Receive transform data from the socket
    ClientIGT->socket->Receive(imgMsg->GetPackBodyPointer(), imgMsg->GetPackBodySize());
    // Deserialize the transform data // If you want to skip CRC check, call Unpack() without argument.
    int c = imgMsg->Unpack(1);
    if (c & igtl::MessageHeader::UNPACK_BODY) // if CRC check is OK
    {
      image->SetParametersFromIGT();
      image->IGTtoITKImage();
      std::cerr<<"Image received"<<std::endl;
      return 0;
    }
  }
  else
  {
    ClientIGT->socket->Skip(headerMsg->GetBodySizeToRead(), 0);
  }

  return 1;

}*/

/*void ProjectImage(Images* image, igtl::Matrix4x4 matrixigt, Images* outputImage)
{
  typedef itk::ResampleImageFilter< ImageType, ImageType >    ResampleFilterType; //types?

  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetInput(image->imageData);

    The Transform that is produced as output of the Registration method is
    also passed as input to the resampling filter. Note the use of the
    methods \code{GetOutput()} and \code{Get()}. This combination is needed
    here because the registration method acts as a filter whose output is a
    transform decorated in the form of a \doxygen{DataObject}. For details in
    this construction you may want to read the documentation of the
    \doxygen{DataObjectDecorator}.
  itk::Matrix<double, 4, 4> M;
  M(0,0) = matrixigt[0][0];
  M(0,1) = matrixigt[0][1];
  M(0,2) = matrixigt[0][2];
  M(0,3) = matrixigt[0][3];
  M(1,0) = matrixigt[1][0];
  M(1,1) = matrixigt[1][1];
  M(1,2) = matrixigt[1][2];
  M(1,3) = matrixigt[1][3];
  M(2,0) = matrixigt[2][0];
  M(2,1) = matrixigt[2][1];
  M(2,2) = matrixigt[2][2];
  M(2,3) = matrixigt[2][3];
  M(3,0) = 0;//matrix[3][0];
  M(3,1) = 0;//matrix[3][1];
  M(3,2) = 0;//matrix[3][2];
  M(3,3) = 1;//matrix[3][3];
  resampler->SetTransform(M);

    As described in Section \ref{sec:ResampleImageFilter}, the
    ResampleImageFilter requires additional parameters to be specified, in
    particular, the spacing, origin and size of the output image. The default
    pixel value is also set to a distinct gray level in order to highlight
    the regions that are mapped outside of the moving image.
  resampler->SetSize(image->imageData->GetLargestPossibleRegion().GetSize());
  resampler->SetOutputOrigin(image->imageData->GetOrigin());
  resampler->SetOutputSpacing(image->imageData->GetSpacing());
  resampler->SetOutputDirection(image->imageData->GetDirection());
  resampler->SetDefaultPixelValue(0);
  resampler->Update();
  outputImage->imageData = resampler->GetOutput();

  float position[3];
  float orientation[4];

  // random position
  static float phi = 0.0;
  position[0] = 0;//100.0 * cos(phi);
  position[1] = 0;//100.0 * sin(phi);
  position[2] = 100.0 * cos(phi);;
  phi = phi + 0.1;

  // random orientation
  matrix1[0][0] = 1.0;  matrix1[1][0] = 0.0;  matrix1[2][0] = 0.0;
  matrix1[0][1] = 0.0;  matrix1[1][1] = 1.0;  matrix1[2][1] = 0.0;
  matrix1[0][2] = 0.0;  matrix1[1][2] = 0.0;  matrix1[2][2] = 1.0;
  matrix1[0][3] = 0.0;  matrix1[1][3] = 0.0;  matrix1[2][3] = 0.0;

  igtl::QuaternionToMatrix(orientation, matrix1);

  matrix1[0][3] = position[0];
  matrix1[1][3] = position[1];
  matrix1[2][3] = position[2];
  outputImage->imgMsg->SetMatrix(matrix1);
  return;
}*/

/*int RegisteredImage( Images* movingImagep, Images* secondImagep, Images* registeredImagep, RegistrationType::Pointer registration )
{

  // Use resulting transform from the registration to map the moving image into the moving/fixed image space
  ResampleFilterType::Pointer   resampler = ResampleFilterType::New();

  // Set moving image as input
  resampler->SetInput( movingImagep->imageData );

  // The Transform produced by the Registration method is passed into the resampling filter
  resampler->SetTransform( registration->GetOutput()->Get() );

  // Specifying parameters of the output image (default pixel value is set to gray in order to highlight the regions that are mapped outside of the moving image)
  resampler->SetSize( secondImagep->imageData->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputOrigin( secondImagep->imageData->GetOrigin() );
  resampler->SetOutputSpacing( secondImagep->imageData->GetSpacing() );
  resampler->SetOutputDirection( secondImagep->imageData->GetDirection() );
  resampler->SetDefaultPixelValue( 0 );
  resampler->Update();

  // Create registered ITKimage
  registeredImagep->imageData = resampler->GetOutput();

  // Set orientation of the registered image
  igtl::Matrix4x4           matrix;
  secondImagep->imgMsg->GetMatrix( matrix );
  registeredImagep->imgMsg->SetMatrix( matrix );

  return 1;

}*/

/*int RegistrationFunction( Images* fixedImagep, Images* movingImagep, Images* registeredImagep )
{

  // Create components registration function
  MetricType::Pointer         metric        = MetricType::New();
  TransformType::Pointer      transform     = TransformType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();

  // Each component is now connected to the instance of the registration method
  registration->SetMetric( metric );
  registration->SetOptimizer( optimizer );
  registration->SetTransform( transform );
  registration->SetInterpolator( interpolator );

  // Set the registration inputs
  registration->SetFixedImage( fixedImagep->imageData );
  registration->SetMovingImage( movingImagep->imageData );
  registration->SetFixedImageRegion( fixedImagep->imageData->GetLargestPossibleRegion() );

  //  Initialize the transform
  ParametersType initialParameters( transform->GetNumberOfParameters() );

  // USE TRACKER DATA TO SET THE INITIAL PARAMETERS!!!!!!!!!!!!!
  // rotation matrix
  initialParameters[0] = 1.0;
  initialParameters[1] = 0.0;
  initialParameters[2] = 0.0;
  initialParameters[3] = 1.0;
  // translation vector
  initialParameters[4] = 0.0;
  initialParameters[5] = 0.0;

  registration->SetInitialTransformParameters( initialParameters );

  optimizer->SetMaximumStepLength( .2 ); // If this is set too high, you will get a "itk::ERROR: MeanSquaresImageToImageMetric(0xa27ce70): Too many samples map outside moving image buffer: 1818 / 10000" error
  optimizer->SetMinimumStepLength( 0.05 );

  // Set a stopping criterion
  optimizer->SetNumberOfIterations( 100 );

  try
  {
    registration->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return 0;// EXIT_FAILURE;
  }

  // Create registered version of moving image
  if ( 0 == RegisteredImage( movingImagep, fixedImagep, registeredImagep, registration ) )
  {
    std::cerr << "Registering Image failed" << std::endl;
  }

  //  The result of the registration process is an array of parameters that defines the spatial transformation in an unique way. This final result is obtained using the \code{GetLastTransformParameters()} method.
  ParametersType finalParameters = registration->GetLastTransformParameters();
  std::cout << "Final parameters 2D/2D-Registration: " << finalParameters << std::endl;

  //  The value of the image metric corresponding to the last set of parameters can be obtained with the \code{GetValue()} method of the optimizer.
  std::cout << "Metric value: " << optimizer->GetValue() << std::endl;

  return 1;

}*/

/*int LoadVolumeVTK( char* filename, Volumes* volume )
{

  // Open file reader and set NRRD file as input
  volume->VTKreader->SetFileName( filename );
  volume->VTKreader->Update();

  // Set parameters
  double    spacing[3];
  double    origin[3];
  int       dimensions[3];
  volume->VTKreader->GetOutput()->GetSpacing( spacing );
  volume->VTKreader->GetOutput()->GetOrigin( origin );
  volume->VTKreader->GetOutput()->GetDimensions( dimensions );
  bool RAS = true; //???
  if ( RAS )
  {
    origin[0] = -1 * origin[0];
    origin[1] = -1 * origin[1];
  }
  volume->SetParametersFromVTK( origin, spacing, dimensions );

  std::cerr << "VTKVolume Loaded " << std::endl;

  return 1;

}*/

/*int LoadVolume( char* filename, Volumes* volume )
{

  // Open file reader and set NRRD file as input
  FileReaderType::Pointer reader = FileReaderType::New();
  reader->SetFileName(filename);

  try
  {
    // Read volume from file
    reader->Update();
    volume->volumeData = reader->GetOutput();
    std::cerr << "Volume dimensions: " << volume->volumeData->GetLargestPossibleRegion().GetSize() << std::endl;
  }
  catch (itk::ExceptionObject &ex)
  {
    std::cerr << ex << std::endl;
    return 0;
  }

  // Set parameters
  volume->SetParametersFromITK();

  std::cerr << "Volume Loaded" << std::endl;

  return 1;

}*/

/*int resliceImageVolumeVTK( vtkSmartPointer<vtkNrrdReader> VTKreader, int start[3], double transformmatrix[16], Image* sliceImage ) // Not finished yet
{

  // Get image parameters from vtk volume(nrrd file)
  double    spacing[3];
  double    origin[3];
  int       dimensions[3];

  VTKreader->Update();
  VTKreader->GetOutput()->GetSpacing( spacing );
  VTKreader->GetOutput()->GetOrigin( origin );
  bool RAS = true; //???
  if ( RAS )
  {
    origin[0] = -1 * origin[0];
    origin[1] = -1 * origin[1];
  }
  VTKreader->GetOutput()->GetDimensions( dimensions );

  std::cerr<< "origin: " << origin[0] << "," << origin[1] << "," << origin[2] << std::endl;
  std::cerr<< "spacing: " << spacing[0] << "," << spacing[1] << "," << spacing[2] << std::endl;
  std::cerr<< "dimensions: " << dimensions[0] << "," << dimensions[1] << "," << dimensions[2] << std::endl;

  // Find the centre
  double originslice[3];
  originslice[0] = origin[0] + start[0];//origin[0] + spacing[0] * 0.5 * dimensions[0] + start[0];
  originslice[1] = origin[1] + start[1];//origin[1] + spacing[1] * 0.5 * dimensions[1] + start[1];
  originslice[2] = origin[2] + start[2];//origin[2] + spacing[2] * 0.5 * dimensions[2] + start[2];

  // Set the direction matrix for the reslice function
  vtkSmartPointer<vtkMatrix4x4> resliceAxes = vtkSmartPointer<vtkMatrix4x4>::New();
  resliceAxes->DeepCopy( transformmatrix );

  // Set the point through which to slice
  resliceAxes->SetElement( 0, 3, originslice[0] );
  resliceAxes->SetElement( 1, 3, originslice[1] );
  resliceAxes->SetElement( 2, 3, originslice[2] );

  // Extract a slice in the desired orientation through the desired point
  vtkSmartPointer<vtkImageReslice> reslice = vtkSmartPointer<vtkImageReslice>::New();
  reslice->SetInputConnection( VTKreader->GetOutputPort() );
  reslice->SetOutputDimensionality( 2 );
  reslice->SetResliceAxes( resliceAxes );
  reslice->SetInterpolationModeToLinear();
  reslice->Update();

  // Create a greyscale lookup table
  vtkSmartPointer<vtkLookupTable> table = vtkSmartPointer<vtkLookupTable>::New();
  table->SetRange(0, 2000); // image intensity range
  table->SetValueRange(0.0, 1.0); // from black to white
  table->SetSaturationRange(0.0, 0.0); // no color saturation
  table->SetRampToLinear();
  table->Build();

  // Map the image through the lookup table
  vtkSmartPointer<vtkImageMapToColors> color = vtkSmartPointer<vtkImageMapToColors>::New();
  color->SetLookupTable(table);
  color->SetInputConnection(reslice->GetOutputPort());

  // Convert the VTK image to an ITK image for further processing
  VTKImageToImageType::Pointer vtkImageToImageFilter = VTKImageToImageType::New();
  vtkImageToImageFilter->SetInput( reslice->GetOutput() );
  vtkImageToImageFilter->Update();
  sliceImage->imageData->Graft( vtkImageToImageFilter->GetOutput() );

  // Set correct image parameters for the ITK image
  sliceImage->imageData->SetOrigin( originslice );
  sliceImage->imageData->SetSpacing( spacing );
  sliceImage->SetParametersFromITK(); //Set spacing and origin in imageData
  std::cerr<< originslice[0]<<", "<< originslice[1]<< ", "<< originslice[2]<< std::endl;

  return 1;

}*/

/*int resliceImageVolume( Volume* volume, int dStart[3], int dsize[2], Image* sliceImage )
{

  // Set parameters for desired image (reslice/ crop)
  VolumeType::IndexType         desiredStart;
  VolumeType::SizeType          desiredSize ;
  desiredStart[0] = dStart[0];  desiredStart[1] = dStart[1];    desiredStart[2] = dStart[2];
  desiredSize[0] = dsize[0];    desiredSize[1] = dsize[1];      desiredSize[2] = 0;

  VolumeType::RegionType        desiredRegion( desiredStart, desiredSize );
  std::cout << "Desired Region: " << desiredRegion << std::endl;
  std::cout << "Desired Region2: " << dStart[1] << std::endl;
  std::cout << "Desired Region2: " << desiredStart[1] << std::endl;

  // Create cropping/reslice
  FilterType::Pointer           filter = FilterType::New();
  filter->SetExtractionRegion( desiredRegion );
  filter->SetInput( volume->volumeData );
  filter->SetDirectionCollapseToIdentity(); // This is required.
  filter->Update();

  // Set resliced image to ITK image
  sliceImage->imageData = filter->GetOutput();

  // Get and set parameters for reslices image
  VolumeType::SpacingType       spacing;       // spacing (mm/pixel)
  float                         spacingsliceImage[3];
  VolumeType::IndexType         start1;
  VolumeType::PointType         origin;
  float                         startsliceImage[3];
  float                         originsliceImage[3];
  VolumeType::SizeType          size1;

  start1 = volume->volumeData->GetLargestPossibleRegion().GetIndex();
  spacing = volume->volumeData->GetSpacing();
  origin = volume->volumeData->GetOrigin();
  size1 = volume->volumeData->GetLargestPossibleRegion().GetSize();
  startsliceImage[0] = start1[0]+dStart[0]; startsliceImage[1] = start1[1]+dStart[1];   startsliceImage[2] = start1[2]+dStart[2];
  spacingsliceImage[0] = spacing[0];        spacingsliceImage[1] = spacing[1];          spacingsliceImage[2] = spacing[2];
  originsliceImage[0] = origin[0] + ( (size1[0]-1) * spacing[0]/2 ) + dStart[0];         originsliceImage[1] = origin[1] + ( (size1[1]-1) * spacing[1]/2 ) + dStart[1];       originsliceImage[2] = origin[2] + ( (size1[2]-1) * spacing[2]/2 ) + dStart[2];

  sliceImage->imageData->SetSpacing( spacingsliceImage );
  sliceImage->imageData->SetOrigin( originsliceImage );

  sliceImage->imgMsg->SetOrigin( originsliceImage );
  sliceImage->imgMsg->SetSpacing( spacingsliceImage );

  return 1;

}*/

int main(int argc, char* argv[]) // Why is this one slow? and why does it stop to recognize image messages after the first ca. 20?
{

  if (argc != 8) // check number of arguments
  {
    // If not correct, print usage
    std::cerr << "    <filenameVolume>    : Filename of Volume.nrrd"            << std::endl;
    std::cerr << "    <hostnameReceiver>  : IP or host name"                    << std::endl;
    std::cerr << "    <portReceiver>      : Port # (18944 default)"             << std::endl;
    std::cerr << "    <hostnameSender1>   : IP or host name"                    << std::endl;
    std::cerr << "    <portSender1>       : Portimage # (18944 default)"        << std::endl;
    std::cerr << "    <hostnameSender2>   : IP or host name"                    << std::endl;
    std::cerr << "    <portSender2>       : Portimage # (18944 default) etc.."  << std::endl;
    exit(0);
  }

  // Load volume data
  char* file = argv[1];
  Volume volume;

  // Load ITK volume data
  volume.LoadVolume( file );

  // Establish connections
  ClientIGT clientIGT1;
  clientIGT1.ConnectToServer( argv[2],atoi( argv[3] ) );
  ClientIGT clientIGT2;
  clientIGT2.ConnectToServer( argv[4],atoi( argv[5] ) );
  ClientIGT clientIGT3;
  clientIGT3.ConnectToServer( argv[6],atoi( argv[7] ) );

  // Send ITK volume to Slicer
  volume.ConvertITKtoIGTVolume();
  clientIGT2.imgMsg = volume.imgMsg;
  clientIGT2.SendImage();

  // Create images
  Image    fixedImage;
  Image    movingImage;
  Image    registeredImage;
  Image    sliceImage;

  // Test resliceImage volume, start point relative to volume coordinates (through which to slice)
  /*int       dStart[3];  dStart[0]=0;        dStart[1]=0;    dStart[2]=0;
  // Test matrix for reslicing axial
  static double transformmatrix[16] = {
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1 };

  //resliceImageVolume(&volume, dStart, dSize, &sliceImage);
  resliceImageVolumeVTK( volumeVTK.VTKreader, dStart, transformmatrix, &sliceImage );

  // Send the fixed image to Slicer
  sliceImage.ITKtoIGTImage();
  sliceImage.imgMsg->SetDeviceName( "sliceImage" );
  sliceImage.imgMsg->Pack();
  clientIGT2.socket->Send( sliceImage.imgMsg->GetPackPointer(), sliceImage.imgMsg->GetPackSize() );

  // Create a message buffer to receive header and image message
  igtl::MessageHeader::Pointer headerMsg = igtl::MessageHeader::New();

  // Allocate a time stamp, WHERE IS THIS USED?????????
  igtl::TimeStamp::Pointer ts = igtl::TimeStamp::New();
  // Receive first image
  ReceiveImage( &clientIGT1, &fixedImage, ts, headerMsg, fixedImage.imgMsg );
*/

  clientIGT1.ReceiveImage();
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
        //registeredImage.SetParametersFromITK();

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
        clientIGT3.imgMsg = fixedImage.imgMsg;
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
  clientIGT3.socket->CloseSocket();
}

