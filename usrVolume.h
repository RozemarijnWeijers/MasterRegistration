#ifndef VOLUME_H
#define VOLUME_H

#include "itkImageFileReader.h"
#include "igtlImageMessage.h"
#include <itkExtractImageFilter.h>
#include <vtkNrrdReader.h>
#include "vtkSmartPointer.h"
#include "usrTransformMatrix.h"


typedef  unsigned char   PixelType;
typedef  itk::Image< PixelType, 3 >  VolumeType;
typedef  itk::ExtractImageFilter< VolumeType, VolumeType > FilterType;

class Volume
{

  public:

  Volume();
  ~Volume();

  void SetParametersFromITK( bool = false );// TransformMatrix );
  void SetParametersFromIGT();
  void ConvertITKtoIGTVolume();
  void ConvertIGTtoITKVolume();
  int LoadVolume( char* );
  void UpdateVolumeTransform( TransformMatrix );
  void CropVolume( float[3], int[3], Volume* );

  VolumeType::Pointer volumeData;
  vtkSmartPointer<vtkNrrdReader> VTKReader;
  igtl::ImageMessage::Pointer imgMsg;
  TransformMatrix volumeMatrix;

  float originVolume[3];
  float spacingVolume[3];
  int sizeVolume[3];

};

#endif //VOLUME_H
