#include "usrClientIGT.h"
#include "usrVolume.h"
#include "usrImage.h"
#include "usrVolumeRegistration.h"
#include "usrTransformMatrix.h"
#include "usrVolumeReslice.h"
#include "usrVolumeCropping.h"
#include "usrImageCropping.h"
#include "usrRotationMatrix.h"


int main(int argc, char* argv[])
{
    if (argc != 4) // check number of arguments
    {
        // If not correct, print usage
        std::cerr << "    <filenameVolume>    : Filename of Volume.nrrd"            << std::endl;
        std::cerr << "    <hostnameReceiver>  : IP or host name"                    << std::endl;
        std::cerr << "    <portReceiver>      : Port # (18944 default)"             << std::endl;
        exit(0);
    }

    char*     file = argv[1];
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

    volume.volumeMatrix.ShowMatrix();

    // Create registration volumes
    Volume    fixedVolume;
    Volume    movingVolume;

    /*TransformMatrix movementMatrix;
    TransformType3D::Pointer transformMovement = TransformType3D::New();
    typedef TransformType3D::VersorType  VersorType;
    typedef VersorType::VectorType     VectorType;
    VersorType     rotation;
    VectorType     axis;
    double angle;
    axis[0] = 0.0;      axis[1] = 0.0;      axis[2] = 1.0;
    typedef TransformType3D::OffsetType OffsetType;
    OffsetType offset;
    rotation.Set(  axis, angle  );
    transformMovement->SetRotation( rotation );*/
    VolumeReslice resliceVolume;
    resliceVolume.SetVolume( &volume );
    //resliceVolume.SetSpacingOfReslice();
    //resliceVolume.SetOriginOfReslice();
    RotationMatrix rotMat;
    double angles[3] = { 0, 0, 0};
    rotMat.Set3Angels( angles );
    float resliceOrigin[3] = {0, 0, 50};
    std::cerr<< "rotmatrix"<< std::endl;
    rotMat.ShowMatrix();
    resliceVolume.SetResliceAxesWRTVolume( rotMat.matrixDouble );
    resliceVolume.SetOriginOfResliceWRTVolume( resliceOrigin );
    resliceVolume.ResliceVolume();
    resliceVolume.CreateITKReslice();
    //resliceVolume.reslicedImage.imageMatrix.ShowMatrix();

    resliceVolume.reslicedImage.ConvertITKtoIGTImage();
    client1.imgMsg = resliceVolume.reslicedImage.imgMsg;
    client1.imgMsg->SetDeviceName( "Volume1" );
    client1.SendImage();

    ImageCropping cropReslice;
    cropReslice.SetImage( &resliceVolume.reslicedImage );
    int dsize[3] = {100, 100};
    int dstart[3] = { 60, 50};
    cropReslice.SetCropSizeAndStart( dsize, dstart );
    cropReslice.CropImage();
    cropReslice.Convert2DImageTo3DVolume();

    cropReslice.croppedVolume.ConvertITKtoIGTVolume();
    cropReslice.croppedImage.ConvertITKtoIGTImage();
    client1.imgMsg = cropReslice.croppedVolume.imgMsg;
    client1.imgMsg->SetDeviceName( "Volume2" );
    client1.SendImage();

    int marge = 0;
    int margez = 4;
    int volumeCropSize[3] = { dsize[0]+marge, dsize[1]+marge, margez };
    int volumeCropOrigin[3] = { (resliceOrigin[0]+dstart[0])/volume.spacingVolume[0]-(marge/2), (resliceOrigin[1]+dstart[1])/volume.spacingVolume[1]-(marge/2), (resliceOrigin[2])/volume.spacingVolume[2]-(margez/2)};
    Volume CroppedVolume;
    volume.CropVolume( volumeCropOrigin, volumeCropSize, &CroppedVolume );

    Volume CroppedVolume2;
    volume.CropVolume( volumeCropOrigin, volumeCropSize, &CroppedVolume2 );

    CroppedVolume2.ConvertITKtoIGTVolume();
    client1.imgMsg = CroppedVolume2.imgMsg;
    client1.imgMsg->SetDeviceName( "Volume3" );
    client1.SendImage();

    // Start registration of the resliced images with other resliced volume of the area around it (5 slices before and 5 after) to find the best match
    typedef TransformType3D::ParametersType ParametersType;
    ParametersType    finalParameters;
    double            initMat[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, marge/2, marge/2, margez/2, 1};
    TransformMatrix   initialMatrix;
    initialMatrix.SetTransformFromDouble( initMat );

    fixedVolume = cropReslice.croppedVolume;//CroppedVolume;
    movingVolume = CroppedVolume;//cropReslice.croppedVolume;
    std::cerr<<movingVolume.volumeData->GetLargestPossibleRegion()<<std::endl;
    std::cerr<<movingVolume.volumeData->GetOrigin()<<std::endl;
    std::cerr<<fixedVolume.volumeData->GetLargestPossibleRegion()<<std::endl;
    std::cerr<<fixedVolume.volumeData->GetOrigin()<<std::endl;
    VolumeRegistration registration;
    registration.SetFixedVolume( &fixedVolume );
    // Register the moving and the fixed image and save th emetric values to find the best match
    registration.SetMovingVolume( &movingVolume );
    registration.SetInitialMatrix( initialMatrix );
    registration.RegisterVolumes();
    //values[i] = registration.metricValue;

    /*for( int i = 0; i < 1; i++)
    {
        // Set movement before taking a slice
        angle = 0 * i;
        axis[0] = 0.0;      axis[1] = 0.0;      axis[2] = 1.0;
        offset[0] = 0.0;    offset[1] = 0.0;    offset[2] = 0.0;
        rotation.Set(  axis, angle  );
        transformMovement->SetRotation( rotation );
        transformMovement->SetOffset( offset );

        // Set parameters for testing resliceVolume volume, get (part of) a slice from the volume
        int       dStart[3];
        dStart[0] = offset[0];         dStart[1] = offset[1];        dStart[2] =offset[2]; // slice number "slicenumber" without the 5 most left pixels (so translated to the left)
        int       dSize[3];
        dSize[0] = 100;  dSize[1] = 100;  dSize[2] = 1;// size?

        // Get (part of) a slice of the volume (by cropping) and set it as fixed image for registration
        volume.CropVolume( dStart, dSize, &fixedVolume );

        // Set the right tracker information for movement

        // Send the fixed image to Slicer
        fixedVolume.ConvertITKtoIGTVolume();
        client1.imgMsg = fixedVolume.imgMsg;
        client1.imgMsg->SetDeviceName( "ImageReslice" );
        client1.SendImage();

        // Start registration of the resliced images with other resliced volume of the area around it (5 slices before and 5 after) to find the best match
        typedef TransformType3D::ParametersType ParametersType;
        ParametersType    finalParameters;
        double            initMat[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
        TransformMatrix   initialMatrix;
        initialMatrix.SetTransformFromDouble( initMat );

        VolumeRegistration registration;
        registration.SetFixedVolume( &fixedVolume );

         // Set parameters for the reslice volume to match with the fixed image
         int     thicknessVolume;
         dStart[0] = 0;  dStart[1] = 0;  dStart[2] = sliceNumber - ((thicknessVolume)/2);
         dSize[0] = size[0]-translation[0];  dSize[1] = size[1]-translation[1];  dSize[2] = thicknessVolume;

         // Get the new reslice of the volume (by cropping) and set it as the moving image for registration
         volume.CropVolume( dStart, dSize, &movingVolume );

         // Register the moving and the fixed image and save th emetric values to find the best match
         registration.SetMovingVolume( &movingVolume );
         registration.SetInitialMatrix( initialMatrix );
         registration.RegisterVolumes();
         values[i] = registration.metricValue;

         // Compare given movement with found registration transformation

         // Update the volume accourding to the found tranformation

         // Send the updated volume to slicer
    }*/


     //const clock_t begin_time_reg = clock();
     //std::cout << "Registration time: " << float( clock () - begin_time_reg) / CLOCKS_PER_SEC << std::endl;

     // Close connection
     client1.socket->CloseSocket();

}
