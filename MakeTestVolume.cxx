#include "itkImage.h"
#include "itkTranslationTransform.h"
#include "itkImageFileReader.h"
#include "itkNormalizeImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkTileImageFilter.h"
#include <itkNrrdImageIO.h>

#include "QuickView.h"

typedef itk::Image<unsigned char, 2>  ImageType;
typedef itk::Image<unsigned char, 3>  VolumeType;

static void CreateImages();

int main(int argc, char *argv[])
{
  ImageType::Pointer image = ImageType::New();
  CreateImages();

  return EXIT_SUCCESS;
}

void CreateImages()
{

  int number = 100;
  ImageType::Pointer image[number];
  int begin1 = 10;
  int begin2 = 20;
  int end1 = 190;
  int end2 = 180;
  for ( int i = 0; i < number; i++ )
  {
    image[i] = ImageType::New();
    ImageType::IndexType start;
    start.Fill(0);

    ImageType::SizeType size;
    size.Fill(200);

    ImageType::RegionType region(start, size);
    image[i]->SetRegions(region);
    image[i]->Allocate();
    image[i]->FillBuffer(0);

    // Make a square
    for(unsigned int r = begin1; r < end1; r++)
    {
      for(unsigned int c = begin2; c < end2; c++)
      {
        ImageType::IndexType pixelIndex;
        pixelIndex[0] = r;
        pixelIndex[1] = c;

        image[i]->SetPixel(pixelIndex, 255);
      }
    }
    begin1++;
    begin2++;
    end1--;
    end2--;
  }

  typedef itk::TileImageFilter< ImageType, VolumeType > Filter2DTo3DType;
  Filter2DTo3DType::Pointer filter2DTo3D = Filter2DTo3DType::New();
  itk::FixedArray< unsigned char, 3 > layout;
  layout[0] = 1;
  layout[1] = 1;
  layout[2] = 0;
  filter2DTo3D->SetLayout( layout );

  unsigned int inputImageNumber = 0;
  ImageType::Pointer inputImageTile;

  for (int i = 0; i < number ; i++)
  {
    inputImageTile = image[i];
    filter2DTo3D->SetInput( i, inputImageTile );
  }
  typedef unsigned char PixelType;
  const PixelType defaultValue = 128;

  filter2DTo3D->SetDefaultPixelValue( defaultValue );
  filter2DTo3D->Update();

  itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();

  typedef itk::ImageFileWriter< VolumeType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter2DTo3D->GetOutput() );
  writer->SetImageIO(io);
  writer->SetFileName( "volume.nrrd" );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return ;
    }

  return;

}
