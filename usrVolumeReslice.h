#include "itkImage.h"
#include "itkVTKImageToImageFilter.h"
#include "QuickView.h"

#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkMatrix4x4.h"
#include "vtkImageReslice.h"
#include "vtkImageCast.h"

#include "usrImage.h"
#include "usrVolume.h"

typedef itk::VTKImageToImageFilter<ImageType> VTKImageToImageType;

class VolumeReslice
{

    public:

    VolumeReslice();
    ~VolumeReslice();

    void ResliceVolume();
    void SetResliceAxesWRTVolume( double [9]);
    void SetOriginOfResliceWRTVolume( float[3] );
    void SetVolume( Volume* );
    void CreateITKReslice();

    Image reslicedImage;
    TransformMatrix resliceMatrix;

    private:

    void SetResliceMatrix();
    void SetSpacingOfReslice();

    vtkSmartPointer<vtkImageReslice> reslice;
    vtkSmartPointer<vtkMatrix4x4> resliceAxes;
    VTKImageToImageType::Pointer vtkImageToImageFilter;
    float resliceOriginWRTVolume[3];
    double resliceOrigin[3];
    float resliceSpacing[3];
    double transformDirections[9];
    Volume* volume;
    TransformMatrix resliceMatrixWRTVolume;

    bool axesSetCheck;
    bool volumeSetCheck;
    bool originOfResliceSetCheck;
    bool resliceDoneCheck;

};
