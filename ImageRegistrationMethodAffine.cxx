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

#include "igtlOSUtil.h"
#include "igtlMessageHeader.h"
#include "igtlTransformMessage.h"
#include "igtlPositionMessage.h"
#include "igtlImageMessage.h"
#include "igtlClientSocket.h"
#include "igtlStatusMessage.h"
#include <igtl_util.h>

const    unsigned int    Dimension = 2;
typedef  unsigned char   PixelType;
typedef  itk::Image< PixelType, Dimension >  ImageType;
//  The transform that will map the fixed image into the moving image.
typedef itk::AffineTransform< double, Dimension > TransformType;
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
 int IGTtoITKImage();//igtl::ImageMessage::Pointer imgMsg, ImageType::Pointer imageData);
 int ITKtoIGTImage();
 ImageType::Pointer imageData;
 igtl::ImageMessage::Pointer imgMsg;
};

Images::Images()
{
  imageData  = ImageType::New();
  imgMsg = igtl::ImageMessage::New();
}

int Images::IGTtoITKImage()
{
  // Retrieve the image data
  int   size[3];          // image dimension
  //float spacing[3];       // spacing (mm/pixel)
  //float spacingITK[0];       // spacing (mm/pixel)
  //int   endian;
  //float origin[3];
  igtl::Matrix4x4 matrix; // Image origin and orientation matrix

  //endian = this->imgMsg->GetEndian();
  this->imgMsg->GetDimensions(size);
  //this->imgMsg->GetSpacing(spacing);
  this->imgMsg->GetMatrix(matrix);
  //this->imgMsg->GetOrigin(origin);

  ImageType::RegionType region;
  ImageType::IndexType start;
  start[0] = 0; start[1] = 0;
  ImageType::SizeType sizeregion;
  sizeregion[0] = size[0]; sizeregion[1] = size[1];
  region.SetSize(sizeregion);
  region.SetIndex(start);
  //spacingITK[0]=spacing[0]; spacingITK[1]=spacing[1];

  this->imageData->SetRegions(region);
  //this->imageData->SetSpacing(spacingITK);
  this->imageData->Allocate();
  //std::cerr<<spacing[0]<< " and " << spacing[1]<< " and " << spacing[2]<< std::endl;
  memcpy(this->imageData->GetBufferPointer(), this->imgMsg->GetScalarPointer(), this->imgMsg->GetSubVolumeImageSize());
  this->imageData->Modified();

  return 1;
}

int Images::ITKtoIGTImage()
{
    // Retrieve the image data
    int   sizeIGT[3];
    float spacingIGT[] = {0.0854354, 0.0854352, 0.85353};    // spacing (mm/pixel) nog aanpassen
    int   svsize[]   = {256, 256, 1};       // sub-volume size
    int   svoffset[] = {0, 0, 0};           // sub-volume offset
    ImageType::SizeType size;
    ImageType::SpacingType spacing;
    int scalarType;
    int numComponents;


    size = this->imageData->GetLargestPossibleRegion().GetSize();
    sizeIGT[0] = size[0]; sizeIGT[1] = size[1]; sizeIGT[2] = 1;
    //spacing = this->imageData->GetSpacing();
    //spacingIGT[0] = spacing[0]; spacingIGT[1] = spacing[1]; spacingIGT[2] = spacing[2];
    scalarType = igtl::ImageMessage::TYPE_UINT8;

    //------------------------------------------------------------
    // Create a new IMAGE type message
    this->imgMsg->SetDimensions(sizeIGT);
    this->imgMsg->SetSpacing(spacingIGT);
    this->imgMsg->SetScalarType(scalarType);
    this->imgMsg->SetSubVolume(sizeIGT, svoffset);
    this->imgMsg->AllocateScalars();
    //this->imgMsg->SetDeviceName("SendimgMsg");

    std::cerr<<this->imgMsg->GetSubVolumeImageSize()<< " and " << size << std::endl;

    /* Get random orientation matrix and set it.
    GetRandomTestMatrix(matrix);*/
    igtl::Matrix4x4 matrix;
    matrix[0][0] = 1.0;  matrix[1][0] = 0.0;  matrix[2][0] = 0.0; matrix[3][0] = 0.0;
    matrix[0][1] = 0.0;  matrix[1][1] = 1.0;  matrix[2][1] = 0.0; matrix[3][1] = 0.0;
    matrix[0][2] = 0.0;  matrix[1][2] = 0.0;  matrix[2][2] = 1.0; matrix[3][2] = 0.0;
    matrix[0][3] = 0.0;  matrix[1][3] = 0.0;  matrix[2][3] = 0.0; matrix[3][3] = 1.0;
    this->imgMsg->SetMatrix(matrix);

    memcpy(this->imgMsg->GetScalarPointer(), imageData->GetBufferPointer(), this->imgMsg->GetSubVolumeImageSize());
    // Pack (serialize) and send
    this->imgMsg->Pack();
    return 1;
}

class Clients
{
  public:
  Clients(char*, int);
  igtl::ClientSocket::Pointer socket;
} ;

Clients::Clients(char* host, int port)
{
  socket = igtl::ClientSocket::New();
  int r = socket->ConnectToServer(host, port);
  if ((r) != 0)
  {
    std::cerr << "Cannot connect to the server." << std::endl;
    exit(0);
  }
}

int ReceiveImage(Clients* client, Images* image, igtl::TimeStamp::Pointer ts, igtl::MessageHeader::Pointer headerMsg, igtl::ImageMessage::Pointer imgMsg)
{
  // Initialize receive buffer
  headerMsg->InitPack();

  if (0 == client->socket->Receive(headerMsg->GetPackPointer(), headerMsg->GetPackSize()))
  {
    client->socket->CloseSocket();
    exit(0);
  }

  // Deserialize the header
  headerMsg->Unpack();

  // Get time stamp
  igtlUint32 sec;
  igtlUint32 nanosec;
  headerMsg->GetTimeStamp(ts);
  ts->GetTimeStamp(&sec, &nanosec);

  if (strcmp(headerMsg->GetDeviceType(), "IMAGE") == 0)
  {
    // Create a message buffer to receive transform data
    imgMsg = igtl::ImageMessage::New();
    imgMsg->SetMessageHeader(headerMsg);
    imgMsg->AllocatePack();
    // Receive transform data from the socket
    client->socket->Receive(imgMsg->GetPackBodyPointer(), imgMsg->GetPackBodySize());
    // Deserialize the transform data // If you want to skip CRC check, call Unpack() without argument.
    int c = imgMsg->Unpack(1);
    if (c & igtl::MessageHeader::UNPACK_BODY) // if CRC check is OK
    {
      image->imgMsg = imgMsg;
      image->IGTtoITKImage();
      //GTtoITKImage(imgMsg, fixedImage);
      return 1;
    }
  }
  else
  {
    client->socket->Skip(headerMsg->GetBodySizeToRead(), 0);
  }
  return 0;
}

void RegisterImage(Images* movingImagep, Images* fixedImagep, Images* registeredImagep, RegistrationType::Pointer registration)
{
  /*  It is common, as the last step of a registration task, to use the
    resulting transform to map the moving image into the fixed image space.
    This is easily done with the \doxygen{ResampleImageFilter}.*/
  typedef itk::ResampleImageFilter< ImageType, ImageType >    ResampleFilterType; //types?

  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetInput(movingImagep->imageData);

  /*  The Transform that is produced as output of the Registration method is
    also passed as input to the resampling filter. Note the use of the
    methods \code{GetOutput()} and \code{Get()}. This combination is needed
    here because the registration method acts as a filter whose output is a
    transform decorated in the form of a \doxygen{DataObject}. For details in
    this construction you may want to read the documentation of the
    \doxygen{DataObjectDecorator}.*/
  resampler->SetTransform(registration->GetOutput()->Get() );

  /*  As described in Section \ref{sec:ResampleImageFilter}, the
    ResampleImageFilter requires additional parameters to be specified, in
    particular, the spacing, origin and size of the output image. The default
    pixel value is also set to a distinct gray level in order to highlight
    the regions that are mapped outside of the moving image.*/
  resampler->SetSize(fixedImagep->imageData->GetLargestPossibleRegion().GetSize());
  resampler->SetOutputOrigin(fixedImagep->imageData->GetOrigin());
  resampler->SetOutputSpacing(fixedImagep->imageData->GetSpacing());
  resampler->SetOutputDirection(fixedImagep->imageData->GetDirection());
  resampler->SetDefaultPixelValue(100);
  resampler->Update();
  registeredImagep->imageData = resampler->GetOutput();
  std::cerr<< registeredImagep->imageData->GetLargestPossibleRegion().GetSize()<<std::endl;

    /*The output of the filter is passed to a writer that will store the
    image in a file. An \doxygen{CastImageFilter} is used to convert the
    pixel type of the resampled image to the final type used by the
    writer. The cast and writer filters are instantiated below.
  //typedef itk::ImageFileWriter< ImageType >  WriterType;
  //typedef unsigned char OutputPixelType;
  //typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  //typedef itk::CastImageFilter< ImageType, ImageType > CastFilterType;

  //WriterType::Pointer      writer =  WriterType::New();
  //CastFilterType::Pointer  caster =  CastFilterType::New();
  //writer->SetFileName("registeredImage.png");

  //caster->SetInput( output );
  //writer->SetInput( output);//caster->GetOutput()   );
  //writer->Update();*/
  return;
}

void RegistrationFunction(Images* fixedImagep, Images* movingImagep, Images* registeredImagep)
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

  // Set the registration inputs
  registration->SetFixedImage(fixedImagep->imageData);
  registration->SetMovingImage(movingImagep->imageData);
  registration->SetFixedImageRegion(fixedImagep->imageData->GetLargestPossibleRegion());

  //  Initialize the transform
  typedef RegistrationType::ParametersType ParametersType;
  ParametersType initialParameters(transform->GetNumberOfParameters());

  // rotation matrix
  initialParameters[0] = 1.0;  // R(0,0)
  initialParameters[1] = 0.0;  // R(0,1)
  initialParameters[2] = 0.0;  // R(1,0)
  initialParameters[3] = 1.0;  // R(1,1)
  // translation vector
  initialParameters[4] = 0.0;
  initialParameters[5] = 0.0;

  registration->SetInitialTransformParameters( initialParameters );

  optimizer->SetMaximumStepLength( .1 ); // If this is set too high, you will get a
  //"itk::ERROR: MeanSquaresImageToImageMetric(0xa27ce70): Too many samples map outside moving image buffer: 1818 / 10000" error

  optimizer->SetMinimumStepLength( 0.01 );

  // Set a stopping criterion
  optimizer->SetNumberOfIterations( 200 );

  // Connect an observer
  //CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  //optimizer->AddObserver( itk::IterationEvent(), observer );
  try
  {
    registration->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return;// EXIT_FAILURE;
  }
  RegisterImage(movingImagep, fixedImagep, registeredImagep, registration);
  //  The result of the registration process is an array of parameters that defines the spatial transformation in an unique way. This final result is obtained using the \code{GetLastTransformParameters()} method.
  ParametersType finalParameters = registration->GetLastTransformParameters();
  std::cout << "Final parameters: " << finalParameters << std::endl;

  //  The value of the image metric corresponding to the last set of parameters can be obtained with the \code{GetValue()} method of the optimizer.
  const double bestValue = optimizer->GetValue();

  // Print out results
  std::cout << "Result = " << std::endl;
  std::cout << " Metric value  = " << bestValue          << std::endl;
  return;
}

int main(int argc, char* argv[])
{

  if (argc != 7) // check number of arguments
  {
    // If not correct, print usage
    std::cerr << "    <hostnameReceiver> : IP or host name"                     << std::endl;
    std::cerr << "    <portReceiver>     : Port # (18944 default)"              << std::endl;
    std::cerr << "    <hostnameSender1>   : IP or host name"                    << std::endl;
    std::cerr << "    <portSender1>       : Portimage # (18944 default)"        << std::endl;
    std::cerr << "    <hostnameSender2>   : IP or host name"                    << std::endl;
    std::cerr << "    <portSender2>       : Portimage # (18944 default) etc.."  << std::endl;
    exit(0);
  }

  // Establish Connection
  Clients client1(argv[1],atoi(argv[2]));
  Clients client2(argv[3],atoi(argv[4]));
  Clients client3(argv[5],atoi(argv[6]));

  // Create two images
  Images fixedImage;
  Images movingImage;
  Images registeredImage;

    // Create a message buffer to receive header
  igtl::MessageHeader::Pointer headerMsg;
  headerMsg = igtl::MessageHeader::New();

  // Allocate a time stamp
  igtl::TimeStamp::Pointer ts;
  ts = igtl::TimeStamp::New();
  igtl::Matrix4x4 matrix;
  igtl::ImageMessage::Pointer imgMsg;

  // Receive first image
  ReceiveImage(&client1, &fixedImage, ts, headerMsg, imgMsg);

  while (1)
  {
    for (int i = 0; i < 100; i ++)
    {
      if (ReceiveImage(&client1, &movingImage, ts, headerMsg, imgMsg) == 1)
      {
        RegistrationFunction(&fixedImage, &movingImage, &registeredImage);
        //registeredImage.imageData = fixedImage.imageData;
        //std::cerr<< registeredImage.imageData->GetLargestPossibleRegion().GetSize()<<std::endl;
        registeredImage.ITKtoIGTImage();

        typedef itk::SubtractImageFilter <ImageType, ImageType > SubtractImageFilterType;
        SubtractImageFilterType::Pointer subtractFilter = SubtractImageFilterType::New ();
        subtractFilter->SetInput1(fixedImage.imageData);
        subtractFilter->SetInput2(registeredImage.imageData);
        subtractFilter->Update();

        typedef itk::SubtractImageFilter <ImageType, ImageType > SubtractImageFilterType2;
        SubtractImageFilterType::Pointer subtractFilter2 = SubtractImageFilterType::New ();
        subtractFilter2->SetInput1(fixedImage.imageData);
        subtractFilter2->SetInput2(movingImage.imageData);
        subtractFilter2->Update();

        Images subtract1;
        subtract1.imageData = subtractFilter->GetOutput();
        subtract1.ITKtoIGTImage();
        subtract1.imgMsg->SetDeviceName("FixedMinusRegistered");
        subtract1.imgMsg->Pack();

        Images subtract2;
        subtract2.imageData = subtractFilter2->GetOutput();
        subtract2.ITKtoIGTImage();
        subtract2.imgMsg->SetDeviceName("FixedMinusMoving");
        subtract2.imgMsg->Pack();

        movingImage.imgMsg->Pack();
        registeredImage.imgMsg->Pack();
        //client2.socket->Send(movingImage.imgMsg->GetPackPointer(), movingImage.imgMsg->GetPackSize());
        client2.socket->Send(subtract1.imgMsg->GetPackPointer(), subtract1.imgMsg->GetPackSize());
        //client3.socket->Send(registeredImage.imgMsg->GetPackPointer(), registeredImage.imgMsg->GetPackSize());
        client3.socket->Send(subtract2.imgMsg->GetPackPointer(), subtract2.imgMsg->GetPackSize());
        Images temp = movingImage;
        movingImage = fixedImage;
        fixedImage = temp;
        //memcpy(fixedImage, movingImage, fixedImage->GetSize());
       }
    }
  }

  // Close connection (The example code never reaches this section ...)
  client1.socket->CloseSocket();
  client2.socket->CloseSocket();
  client3.socket->CloseSocket();
}

