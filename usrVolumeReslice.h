#include "itkImage.h"
#include "itkVTKImageToImageFilter.h"
#include "QuickView.h"

#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkMatrix4x4.h"
#include "vtkImageReslice.h"
#include "vtkLookupTable.h"
#include "vtkImageMapToColors.h"
#include "vtkImageActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkImageMapper3D.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleImage.h"
#include "itkRescaleIntensityImageFilter.h"

#include "usrImage.h"
#include "usrVolume.h"

typedef itk::VTKImageToImageFilter<ImageType> VTKImageToImageType;

class VolumeReslice
{

    public:

    VolumeReslice();
    ~VolumeReslice();

    void ResliceVolume();
    void SetResliceAxes( double [9]);
    void SetOriginOfResliceWRTVolume( double[3] );
    void SetVolume( Volume* );
    void CreateITKReslice();

    Image reslicedImage;


    private:

    void SetOriginOfReslice();
    void SetSpacingOfReslice();

    vtkSmartPointer<vtkImageReslice> reslice;
    vtkSmartPointer<vtkMatrix4x4> resliceAxes;
    VTKImageToImageType::Pointer vtkImageToImageFilter;
    double resliceOriginWRTVolume[3];
    double resliceOrigin[3];
    double resliceSpacing[3];
    Volume* volume;

    bool axesSetCheck;
    bool volumeSetCheck;
    bool originOfResliceSetCheck;
    bool resliceDoneCheck;

};
