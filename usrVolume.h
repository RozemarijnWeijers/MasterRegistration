#ifndef VOLUME_H
#define VOLUME_H

#include "itkImageFileReader.h"
#include "igtlImageMessage.h"
#include <vtkNrrdReader.h>
#include "vtkSmartPointer.h"


typedef  unsigned char   PixelType;
typedef  itk::Image< PixelType, 3 >  VolumeType;

class Volume
{

  public:

  Volume();
  ~Volume();

  void SetParametersFromITK();
  void ConvertITKtoIGTVolume();
  int LoadVolume( char* );

  VolumeType::Pointer volumeData;
  vtkSmartPointer<vtkNrrdReader> VTKReader;
  igtl::ImageMessage::Pointer imgMsg;

  float originVolume[3];
  float spacingVolume[3];

  private:

  int sizeVolume[3];

};

#endif //VOLUME_H
