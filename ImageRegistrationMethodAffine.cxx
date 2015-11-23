#include "itkCastImageFilter.h"
#include "itkEllipseSpatialObject.h"
#include "itkImage.h"
#include "itkImageRegistrationMethod.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSpatialObjectToImageFilter.h"
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

#include "igtlOSUtil.h"
#include "igtlMessageHeader.h"
#include "igtlTransformMessage.h"
#include "igtlPositionMessage.h"
#include "igtlImageMessage.h"
#include "igtlClientSocket.h"
#include "igtlStatusMessage.h"
#include <igtl_util.h>

const    unsigned int    DimensionImage = 2;
const    unsigned int    DimensionVolume = 3;
typedef  unsigned char   PixelType;
typedef  itk::Image< PixelType, DimensionImage >  ImageType;
typedef  itk::Image< PixelType, DimensionVolume >  VolumeType;

//  The transform that will map the fixed image into the moving image.
typedef itk::AffineTransform< double, DimensionImage > TransformType;
// typedef itk::TranslationTransform< double, Dimension > TransformType;
//  An optimizer is required to explore the parameter space of the transform in search of optimal values of the metric.
typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
//  The metric will compare how well the two images match each other. Metric types are usually parameterized by the image types as it can be seen in the following type declaration.
typedef itk::MeanSquaresImageToImageMetric< ImageType, ImageType > MetricType;
//  Finally, the type of the interpolator is declared. The interpolator will evaluate the intensities of the moving image at non-grid positions.
typedef itk::LinearInterpolateImageFunction< ImageType, double > InterpolatorType;
//  The registration method type is instantiated using the types of the fixed and moving images. This class is responsible for interconnecting all the components that we have described so far.
typedef itk::ImageRegistrationMethod< ImageType, ImageType >    RegistrationType;


class Images
{
  public:

  Images();
  int IGTtoITKImage();
  int ITKtoIGTImage();
  int SetParametersFromIGT();
  int SetParametersFromITK();

  ImageType::Pointer imageData;
  igtl::ImageMessage::Pointer imgMsg;

  protected:

  int sizeIm[3];
  float originIm[3];
  float spacingIm[3];

};

Images::Images()
{

  // Create imageMessage and ITKimage for this image
  imageData  = ImageType::New();
  imgMsg = igtl::ImageMessage::New();

}

int Images::SetParametersFromIGT()
{

  this->imgMsg->GetDimensions(this->sizeIm);
  this->imgMsg->GetOrigin(this->originIm);
  this->imgMsg->GetSpacing(this->spacingIm);

  return 1;

}

int Images::SetParametersFromITK()
{

  ImageType::SizeType size = this->imageData->GetLargestPossibleRegion().GetSize();
  sizeIm[0] = size[0]; sizeIm[1] = size[1]; sizeIm[2] = 1;
  ImageType::PointType origin = this->imageData->GetOrigin();
  originIm[0] = origin[0]; originIm[1] = origin[1]; originIm[2] = origin[2];
  ImageType::SpacingType spacing = this->imageData->GetSpacing();
  spacingIm[0] = spacing[0]; spacingIm[1] = spacing[1]; spacingIm[2] = spacing[2];

  return 1;

}

int Images::IGTtoITKImage()
{

  // Retrieve the image data from image message
  int   size[3];          // image dimension (pixels)
  float spacing[3];       // spacing (mm/pixel)
  float spacingITK[0];    // spacing for 2D image (mm/pixel)
  int   endian;           // endian (not used)
  float origin[3];        // origin ()
  igtl::Matrix4x4 matrix; // image origin and orientation matrix

  endian = this->imgMsg->GetEndian();
  this->imgMsg->GetDimensions(size);
  this->imgMsg->GetSpacing(spacing);
  this->imgMsg->GetMatrix(matrix);  // needed for projection?
  this->imgMsg->GetOrigin(origin);

  // Set image data do ITK image
  ImageType::RegionType region;
  ImageType::IndexType start;
  start[0] = origin[0]; start[1] = origin[1]; // not sur if this is correct
  ImageType::SizeType sizeregion;
  sizeregion[0] = size[0]; sizeregion[1] = size[1];
  region.SetSize(sizeregion);
  region.SetIndex(start);
  spacingITK[0]=spacing[0]; spacingITK[1]=spacing[1];

  this->imageData->SetRegions(region);
  this->imageData->SetSpacing(spacingITK);
  this->imageData->Allocate();

  // Set copy image data into ITK image
  memcpy(this->imageData->GetBufferPointer(), this->imgMsg->GetScalarPointer(), this->imgMsg->GetSubVolumeImageSize());
  this->imageData->Modified();

  return 1;

}

int Images::ITKtoIGTImage()
{

  // Retrieve the image data from ITK image
  ImageType::SizeType size;
  int   sizeIGT[3];
  int   svsize[3];//   = {256, 256, 1};       // sub-volume size
  ImageType::IndexType start;
  int   svoffset[3];// = {0, 0, 0};           // sub-volume offset
  int   scalarType;
  //int   numComponents;
  ImageType::SpacingType spacing;
  float spacingIGT[3];// = {0.0854354, 0.0854352, 0.85353};    // spacing (mm/pixel) nog aanpassen

  size = this->imageData->GetLargestPossibleRegion().GetSize();
  sizeIGT[0] = size[0]; sizeIGT[1] = size[1]; sizeIGT[2] = 1;
  svsize[0] = sizeIGT[0]; svsize[1] = sizeIGT[1]; svsize[2] = sizeIGT[2];
  start = this->imageData->GetLargestPossibleRegion().GetIndex();
  svoffset[0] = start[0]; svoffset[1] = start[1]; svoffset[2] = start[2];
  spacing = this->imageData->GetSpacing();
  spacingIGT[0] = spacing[0]; spacingIGT[1] = spacing[1]; spacingIGT[2] = spacing[1];//spacing[2];!!!!!!!!!!!
  scalarType = igtl::ImageMessage::TYPE_UINT8;

  // Create a new IMAGE type message
  this->imgMsg->SetDimensions(sizeIGT);
  this->imgMsg->SetSpacing(spacingIGT);
  this->imgMsg->SetScalarType(scalarType);
  this->imgMsg->SetSubVolume(sizeIGT, svoffset);
  this->imgMsg->AllocateScalars();

  // Set copy image data into ITK image
  memcpy(this->imgMsg->GetScalarPointer(), imageData->GetBufferPointer(), this->imgMsg->GetSubVolumeImageSize());
  // Pack (serialize) and send
  this->imgMsg->Pack();

  return 1;

}

class Volumes
{

  public:

  Volumes();
  int ITKtoIGTVolume();
  int VTKtoITKVolume();

  VolumeType::Pointer volumeData;
  vtkSmartPointer<vtkNrrdReader> readerVTK;
  igtl::ImageMessage::Pointer imgMsg;

  protected:

  int sizeVol[3];
  float originVol[3];
  float spacingVol[3];

};

Volumes::Volumes()
{

  // Create imageMessage and ITKvolume for this volume
  volumeData  = VolumeType::New();
  readerVTK = vtkSmartPointer<vtkNrrdReader>::New();
  imgMsg = igtl::ImageMessage::New();

}

int Volumes::ITKtoIGTVolume()
{

  // Retrieve the image data from ITK image
  VolumeType::SizeType size;
  int   sizeIGT[3];
  int   svsize[3];//   = {256, 256, 1};       // sub-volume size
  VolumeType::IndexType start;
  int   svoffset[3];// = {0, 0, 0};           // sub-volume offset
  int   scalarType;
  VolumeType::PointType origin;
  float originIGT[3];
  //int   numComponents;
  VolumeType::SpacingType spacing;
  float spacingIGT[3];// = {0.0854354, 0.0854352, 0.85353};    // spacing (mm/pixel) nog aanpassen

  size = this->volumeData->GetLargestPossibleRegion().GetSize();
  sizeIGT[0] = size[0]; sizeIGT[1] = size[1]; sizeIGT[2] = size[2];
  svsize[0] = sizeIGT[0]; svsize[1] = sizeIGT[1]; svsize[2] = sizeIGT[2];
  start = this->volumeData->GetLargestPossibleRegion().GetIndex();
  svoffset[0] = start[0]; svoffset[1] = start[1]; svoffset[2] = start[2];
  spacing = this->volumeData->GetSpacing();
  spacingIGT[0] = spacing[0]; spacingIGT[1] = spacing[1]; spacingIGT[2] = spacing[2];
  origin = this->volumeData->GetOrigin();
  originIGT[0] = origin[0]+((size[0]-1)*spacing[0]/2); originIGT[1] = origin[1]+((size[1]-1)*spacing[1]/2); originIGT[2] = origin[2]+((size[2]-1)*spacing[2]/2);
  scalarType = igtl::ImageMessage::TYPE_UINT8;

  // Create a new IMAGE type message
  this->imgMsg->SetDimensions(sizeIGT);
  this->imgMsg->SetSpacing(spacingIGT);
  this->imgMsg->SetScalarType(scalarType);
  this->imgMsg->SetOrigin(originIGT);
  this->imgMsg->AllocateScalars();

  // Set copy image data into ITK image
  memcpy(this->imgMsg->GetScalarPointer(), volumeData->GetBufferPointer(), this->imgMsg->GetImageSize());
  // Pack (serialize) and send
  this->imgMsg->Pack();

  return 1;

}

int VTKtoITKVolume()
{

  return 1;

}

class Clients
{

  public:

  Clients(char*, int);

  igtl::ClientSocket::Pointer socket;

};

Clients::Clients(char* host, int port)
{

  //Open a socket for th client
  socket = igtl::ClientSocket::New();
  //Connect to the server
  int r = socket->ConnectToServer(host, port);
  // Check the connection
  if ((r) != 0)
  {
    std::cerr << "Cannot connect to the server." << std::endl;
    exit(0);
  }
  std::cerr << "Client is connected to server" << host <<":"<< port<< std::endl;

}

int ReceiveImage(Clients* client, Images* image, igtl::TimeStamp::Pointer ts, igtl::MessageHeader::Pointer headerMsg, igtl::ImageMessage::Pointer imgMsg)
{

  // Initialize receive buffer
  headerMsg->InitPack();

  // Receive header message
  if (0 == client->socket->Receive(headerMsg->GetPackPointer(), headerMsg->GetPackSize()))
  {
    client->socket->CloseSocket();
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
    client->socket->Receive(imgMsg->GetPackBodyPointer(), imgMsg->GetPackBodySize());
    // Deserialize the transform data // If you want to skip CRC check, call Unpack() without argument.
    int c = imgMsg->Unpack(1);
    if (c & igtl::MessageHeader::UNPACK_BODY) // if CRC check is OK
    {
      image->IGTtoITKImage();
      std::cerr<<"Image received"<<std::endl;
      return 0;
    }
  }
  else
  {
    client->socket->Skip(headerMsg->GetBodySizeToRead(), 0);
  }

  return 1;

}

/*=void ProjectImage(Images* image, igtl::Matrix4x4 matrixigt, Images* outputImage)
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

int RegisteredImage(Images* movingImagep, Images* secondImagep, Images* registeredImagep, RegistrationType::Pointer registration)
{
  // Use resulting transform from the registration to map the moving image into the moving/fixed image space
  typedef itk::ResampleImageFilter< ImageType, ImageType >    ResampleFilterType;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  // Set moving image as input
  resampler->SetInput(movingImagep->imageData);

  // The Transform produced by the Registration method is passed into the resampling filter
  // Note: the registration method acts as a filter whose output is a transform decorated in the form of a \doxygen{DataObject}. For details, read the documentation of the \doxygen{DataObjectDecorator}
  resampler->SetTransform(registration->GetOutput()->Get());

  // Specifying parameters of the output image (default pixel value is set to gray in order to highlight the regions that are mapped outside of the moving image)
  resampler->SetSize(secondImagep->imageData->GetLargestPossibleRegion().GetSize());
  resampler->SetOutputOrigin(secondImagep->imageData->GetOrigin());
  resampler->SetOutputSpacing(secondImagep->imageData->GetSpacing());
  resampler->SetOutputDirection(secondImagep->imageData->GetDirection());
  resampler->SetDefaultPixelValue(0);
  resampler->Update();

  // Create registered ITKimage
  registeredImagep->imageData = resampler->GetOutput();

  //Check registered image
  /* if ()
  {
    return 0;
  }/*

  /*float spacing[3];       // spacing (mm/pixel)
  movingImagep->imgMsg->GetSpacing(spacing);
  registeredImagep->imgMsg->SetSpacing(spacing);*/

  // Set orientation of the registered image
  igtl::Matrix4x4 matrix;
  secondImagep->imgMsg->GetMatrix(matrix);
  registeredImagep->imgMsg->SetMatrix(matrix);

  return 1;

}

int RegistrationFunction(Images* fixedImagep, Images* movingImagep, Images* registeredImagep)
{
    // Create components registration function
  MetricType::Pointer         metric        = MetricType::New();
  TransformType::Pointer      transform     = TransformType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();

  // Each component is now connected to the instance of the registration method.
  registration->SetMetric(metric);
  registration->SetOptimizer(optimizer);
  registration->SetTransform(transform);
  registration->SetInterpolator(interpolator);

  /*/ Project moving Image on to fixed Image plane
  igtl::Matrix4x4 T01;
  igtl::Matrix4x4 T02;
  fixedImagep->imgMsg->Getmatrix(T01);
  movingImagep->imgMsg->Getmatrix(T02);
  typedef   itk::Matrix<NumericType,4,4>          MatrixType;
  MatrixType matrix01;
  MatrixType matrix02;

  matrix01(0,0) = T01[0][0];
  matrix01(0,1) = T01[2][0];
  matrix01(0,2) = T01[1][0];
  matrix01(1,0) = T01[0][1];
  matrix01(1,1) = T01[1][1];
  matrix01(1,2) = T01[2][1];
  matrix01(2,0) = T01[0][2];
  matrix01(2,1) = T01[1][2];
  matrix01(2,2) = T01[2][2];
  matrix01(3,0) = T01[0][3];
  matrix01(3,1) = T01[1][3];
  matrix01(3,2) = T01[2][3];
  matrix01(4,4) = 1;

  matrix02(0,0) = T02[0][0];
  matrix02(0,1) = T02[2][0];
  matrix02(0,2) = T02[1][0];
  matrix02(1,0) = T02[0][1];
  matrix02(1,1) = T02[1][1];
  matrix02(1,2) = T02[2][1];
  matrix02(2,0) = T02[0][2];
  matrix02(2,1) = T02[1][2];
  matrix02(2,2) = T02[2][2];
  matrix02(3,0) = T02[0][3];
  matrix02(3,1) = T02[1][3];
  matrix02(3,2) = T02[2][3];
  matrix02(4,4) = 1;*/

  // Set the registration inputs
  registration->SetFixedImage(fixedImagep->imageData);
  registration->SetMovingImage(movingImagep->imageData);
  registration->SetFixedImageRegion(fixedImagep->imageData->GetLargestPossibleRegion());

  //  Initialize the transform
  typedef RegistrationType::ParametersType ParametersType;
  ParametersType initialParameters(transform->GetNumberOfParameters());

  // USE TRACKER DATA TO SET THE INITIAL PARAMETERS!!!!!!!!!!!!!
  // rotation matrix
  initialParameters[0] = 1.0;  // R(0,0)
  initialParameters[1] = 0.0;  // R(0,1)
  initialParameters[2] = 0.0;  // R(1,0)
  initialParameters[3] = 1.0;  // R(1,1)
  // translation vector
  initialParameters[4] = 0.0;
  initialParameters[5] = 0.0;

  registration->SetInitialTransformParameters( initialParameters );

  optimizer->SetMaximumStepLength( .1 ); // If this is set too high, you will get a "itk::ERROR: MeanSquaresImageToImageMetric(0xa27ce70): Too many samples map outside moving image buffer: 1818 / 10000" error
  optimizer->SetMinimumStepLength( 0.01 );

  // Set a stopping criterion
  optimizer->SetNumberOfIterations( 200 );

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

  // Create registered version of moving image. NOTE THAT THE SECOND IMAGE IS NOW ALSO THE MOVING IMAGE!!
  if (0 == RegisteredImage(movingImagep, fixedImagep, registeredImagep, registration))
  {
    std::cerr << "Registering Image failed" << std::endl;
  }

  //  The result of the registration process is an array of parameters that defines the spatial transformation in an unique way. This final result is obtained using the \code{GetLastTransformParameters()} method.
  ParametersType finalParameters = registration->GetLastTransformParameters();
  std::cout << "Final parameters 2D/2D-Registration: " << finalParameters << std::endl;

  //  The value of the image metric corresponding to the last set of parameters can be obtained with the \code{GetValue()} method of the optimizer.
  std::cout << "Metric value  = " << optimizer->GetValue()         << std::endl;

  return 1;

}


int LoadVolumeVTK(char* filename, Volumes* volume)
{

  volume->readerVTK->SetFileName(filename);
  volume->readerVTK->Update();

  /*typedef itk::ImageToVTKImageFilter<VolumeType>  ConnectorType;
  ConnectorType::Pointer connector = ConnectorType::New();
  connector->SetInput(volume->volumeData);*/


  std::cerr << "Volume Loaded" << std::endl;

  return 1;

}

int LoadVolume(char* filename, Volumes* volume)
{

  typedef itk::ImageFileReader<VolumeType> FileReaderType;
  FileReaderType::Pointer reader = FileReaderType::New();
  reader->SetFileName(filename);

  /*std::cerr << "Spacing Volume:" << volume->GetSpacing() << std::endl;
  std::cerr << "Requested Region Volume:" << volume->GetRequestedRegion() << std::endl;
  std::cerr << "Origin Volume:" << volume->GetOrigin() << std::endl;*/
  std::cerr << "VTKVolume Loaded" << std::endl;

  return 1;

}

int resliceImageVolumeVTK(vtkSmartPointer<vtkNrrdReader> readerVTK, int start[3], int size[2], Images* sliceImage)
{

  readerVTK->Update();
  double spacing[3];
  double origin[3];

  static double sagittalElements[16] = {
            0, 0,-1, 0,
            1, 0, 0, 0,
            0,-1, 0, 0,
            0, 0, 0, 1 };

  readerVTK->GetOutput()->GetSpacing(spacing);
  readerVTK->GetOutput()->GetOrigin(origin);

  double center[3];
  center[0] = origin[0] + spacing[0] * 0.5;// * (extent[0] + extent[1]);
  center[1] = origin[1] + spacing[1] * 0.5;// * (extent[2] + extent[3]);
  center[2] = origin[2] + spacing[2] * 0.5;// * (extent[4] + extent[5]);

  vtkSmartPointer<vtkMatrix4x4> resliceAxes = vtkSmartPointer<vtkMatrix4x4>::New();
  resliceAxes->DeepCopy(sagittalElements);
  // Set the point through which to slice
  resliceAxes->SetElement(0, 3, center[0]);
  resliceAxes->SetElement(1, 3, center[1]);
  resliceAxes->SetElement(2, 3, center[2]);

  // Extract a slice in the desired orientation
  vtkSmartPointer<vtkImageReslice> reslice = vtkSmartPointer<vtkImageReslice>::New();
  reslice->SetInputConnection(readerVTK->GetOutputPort());
  reslice->SetOutputDimensionality(2);
  reslice->SetResliceAxes(resliceAxes);
  reslice->SetInterpolationModeToLinear();
  reslice->Update();
  vtkSmartPointer<vtkImageData> slice = vtkSmartPointer<vtkImageData>::New();
  slice= reslice->GetOutput();

  return 1;

}

int resliceImageVolume(Volumes* volume, int start[3], int size[2], Images* sliceImage)
{

  VolumeType::IndexType desiredStart;
  desiredStart[0]=start[0]; desiredStart[1]=start[1]; desiredStart[2]=start[2];
  //desiredStart.Fill(1);

  VolumeType::SizeType desiredSize ;
  desiredSize[0]=size[0]; desiredSize[1]=size[1]; desiredSize[2]=0;
  //desiredSize.Fill(10);
  //esiredSize[2]=0;

  VolumeType::RegionType desiredRegion(desiredStart, desiredSize);

  std::cout << "desiredRegion: " << desiredRegion << std::endl;

  typedef itk::ExtractImageFilter< VolumeType, ImageType > FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetExtractionRegion(desiredRegion);
  filter->SetInput(volume->volumeData);
  #if ITK_VERSION_MAJOR >= 4
  filter->SetDirectionCollapseToIdentity(); // This is required.
  #endif
  filter->Update();

  sliceImage->imageData = filter->GetOutput();
  //output->DisconnectPipeline();
  //output->FillBuffer(2);

  VolumeType::SpacingType spacing;       // spacing (mm/pixel)
  float spacingsliceImage[3];
  VolumeType::IndexType start1;
  VolumeType::PointType origin;
  float startsliceImage[3];
  float originsliceImage[3];
  VolumeType::SizeType size1;

  start1 = volume->volumeData->GetLargestPossibleRegion().GetIndex();
  startsliceImage[0] = start1[0]+start[0]; startsliceImage[1] = start1[1]+start[1]; startsliceImage[2] = start1[2]+start[2];
  spacing = volume->volumeData->GetSpacing();
  spacingsliceImage[0] = spacing[0]; spacingsliceImage[1] = spacing[1]; spacingsliceImage[2] = spacing[2];
  origin = volume->volumeData->GetOrigin();
  size1 = volume->volumeData->GetLargestPossibleRegion().GetSize();
  originsliceImage[0] = origin[0]+((size1[0]-1)*spacing[0]/2)+start[0]; originsliceImage[1] = origin[1]+((size1[1]-1)*spacing[1]/2)+start[1]; originsliceImage[2] = origin[2]+((size1[2]-1)*spacing[2]/2)+start[2];
  //scalarType = igtl::ImageMessage::TYPE_UINT8;

  sliceImage->imageData->SetSpacing(spacingsliceImage);
  sliceImage->imageData->SetOrigin(originsliceImage);
  sliceImage->imgMsg->SetOrigin(originsliceImage);

  sliceImage->imgMsg->SetSpacing(spacingsliceImage);

}

int main(int argc, char* argv[])
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
  Volumes volume;
  LoadVolume(file, &volume);
  //Load VTK volume data
  LoadVolumeVTK(file, &volume);
  //vtkSmartPointer<vtkNrrdReader> readerVTK = vtkSmartPointer<vtkNrrdReader>::New();
  //readerVTK->SetFileName(file);

  // Establish connections
  Clients client1(argv[2],atoi(argv[3]));
  Clients client2(argv[4],atoi(argv[5]));
  Clients client3(argv[6],atoi(argv[7]));

  volume.ITKtoIGTVolume();
  client2.socket->Send(volume.imgMsg->GetPackPointer(), volume.imgMsg->GetPackSize());

  // Create images
  Images fixedImage;
  Images movingImage;
  Images registeredImage;
  Images sliceImage;

  // Test resliceImage volume
  int dStart[3]; dStart[0]=0; dStart[1]=0; dStart[2]=0;
  VolumeType::SizeType size = volume.volumeData->GetLargestPossibleRegion().GetSize();
  int dSize[2]; dSize[0]=size[0]; dSize[1]=size[1];

  //resliceImageVolume(&volume, dStart, dSize, &sliceImage);
  resliceImageVolumeVTK(volume.readerVTK, dStart, dSize, &sliceImage);

  sliceImage.ITKtoIGTImage();
  sliceImage.imgMsg->SetDeviceName("sliceImage");
  sliceImage.imgMsg->Pack();
  client2.socket->Send(sliceImage.imgMsg->GetPackPointer(), sliceImage.imgMsg->GetPackSize());

  // Create a message buffer to receive header and image message
  igtl::MessageHeader::Pointer headerMsg = igtl::MessageHeader::New();

  // Allocate a time stamp, WHERE IS THIS USED?????????
  igtl::TimeStamp::Pointer ts = igtl::TimeStamp::New();

  // Receive first image
  ReceiveImage(&client1, &fixedImage, ts, headerMsg, fixedImage.imgMsg);

  while (1)
  {
    for (int i = 0; i < 100; i ++) //WHY 100????
    {
      if (ReceiveImage(&client1, &movingImage, ts, headerMsg, movingImage.imgMsg) == 0)
      {
        if (fixedImage.imgMsg == movingImage.imgMsg){std::cerr<<"zelfde"<<std::endl;}
        // Register fixed and moving image (2D/2D Registration of two subsequential images)
        RegistrationFunction(&fixedImage, &movingImage, &registeredImage);

        // Create imageMessage for registered image
        registeredImage.ITKtoIGTImage();

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
        fixedImage.imgMsg->SetDeviceName("fixed");
        fixedImage.imgMsg->Pack();
        movingImage.imgMsg->SetDeviceName("moving");
        movingImage.imgMsg->Pack();
        registeredImage.imgMsg->SetDeviceName("registered");
        registeredImage.imgMsg->Pack();

        // Send image messages to sliceImager for visualization
        client3.socket->Send(fixedImage.imgMsg->GetPackPointer(), fixedImage.imgMsg->GetPackSize());
        client3.socket->Send(movingImage.imgMsg->GetPackPointer(), movingImage.imgMsg->GetPackSize());
        //client2.socket->Send(subtract1.imgMsg->GetPackPointer(), subtract1.imgMsg->GetPackSize());
        client3.socket->Send(registeredImage.imgMsg->GetPackPointer(), registeredImage.imgMsg->GetPackSize());
        //client3.socket->Send(subtract2.imgMsg->GetPackPointer(), subtract2.imgMsg->GetPackSize());

        // Set moving image to fixed image so it can be registered to the next image received
        Images temp = movingImage;
        movingImage = fixedImage;
        fixedImage = temp;
       }
    }
    std::cerr<< "Stopped receiving messages" <<std::endl;
  }

  // Close connection (The example code never reaches this section ...)
  client1.socket->CloseSocket();
  client2.socket->CloseSocket();
  client3.socket->CloseSocket();
}

