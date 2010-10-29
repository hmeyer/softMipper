#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <vector>
#include <numeric>
#include <string>
#include <itkCastImageFilter.h>
#include "itkbasics.h"
#include "toString.h"
#include "itkVariableProjectImageFilter.h"
#include "itkFilterFunctions.h"
#include "itkJoinSeriesImageFilter.h"
#include "itkExtractImageFilter.h"
#include <boost/concept_check.hpp>
#include <boost/lexical_cast.hpp>

using std::string;
using namespace itkBasic;

typedef itk::VariableProjectImageFilter< DicomInputImageType, DicomInputImageType, itk::MaximumProjector > MipFilterType;
typedef itk::VariableProjectImageFilter< DicomInputImageType, DicomInputImageType, itk::softMipProjector > softMipFilterType;
typedef itk::VariableProjectImageFilter< DicomInputImageType, DicomInputImageType, itk::AverageProjector > AverageFilterType;

const char PATHSEPARATOR = '/';

template<class TImagePointerType, template < typename, typename > class Projector>
TImagePointerType CreateProjectedSlices(TImagePointerType image, int projectionDimension, double slice_thinkness, double slice_increment,
				     Projector<typename TImagePointerType::ObjectType::PixelType, typename TImagePointerType::ObjectType::PixelType> proj = 
					Projector<typename TImagePointerType::ObjectType::PixelType, typename TImagePointerType::ObjectType::PixelType>());

const string projectionMIP("mip");
const string projectionSoftMIP("softmip");

void parseProjectionCommand( const std::string &cmd, std::string &projectionType, float &softMipStrength, int &dimension, float &thickness, float &increment );


int main( int argc, char* argv[] )
{
	if( argc < 4 )	{
		std::cerr << "Usage: " << argv[0] << " InDicomDirectory OutDicomDirectory [(mip|softmip:strength):dimension:thickness:increment]" << std::endl;
		return EXIT_FAILURE;
	}

	SeriesReader sreader( argv[1] );
	sreader.readSeriesData( 2 );
	string outDir = argv[2];
	FilterImageType::Pointer image;
	{
		ReaderType::Pointer imageReader = ReaderType::New();
		FileNamesContainer fc;
		sreader.getSeriesFileNames(0, fc);
		DicomInputImageType::Pointer image = getDicomSerie( fc, imageReader, 1 );
		int outDicomIndex = 0;

		for(int k = 3; k < argc; k++) {
		  string projectionType;
		  float softMipStrength;
		  int dimension;
		  float thickness;
		  float increment;
		  try {
		    parseProjectionCommand( argv[k], projectionType, softMipStrength, dimension, thickness, increment );
		    DicomInputImageType::Pointer resultSlices;
		    if (projectionType == projectionMIP) {
		      MipFilterType::ProjectorType mipProj;
		      resultSlices = CreateProjectedSlices( image, dimension, thickness, increment, mipProj );
		    } else if (projectionType == projectionSoftMIP) {
		      softMipFilterType::ProjectorType softMipProj( softMipStrength );
		      resultSlices = CreateProjectedSlices( image, dimension, thickness, increment, softMipProj );
		    }
		    replaceMetaData( imageReader->GetMetaDataDictionaryArray(), "0008|0008", "DERIVED\\SECONDARY\\AXIAL");
		    writeDicomSeries( resultSlices, outDir + PATHSEPARATOR + "IM%06d", imageReader, outDicomIndex);
		  
		  } catch( itk::ExceptionObject & exp ) {
		    std::cerr << "Exception caught !" << std::endl;
		    std::cerr << exp << std::endl;
		  }
		}
	}
	return EXIT_SUCCESS;
}

void parseProjectionCommand( const std::string &cmd, std::string &projectionType, float &softMipStrength, int &dimension, float &thickness, float &increment ) {
  const unsigned int TS = 10;
  char t[TS];
  std::istringstream cmdstream( cmd );
  cmdstream.exceptions( std::istringstream::failbit );
  cmdstream.getline( t, TS, ':');
  if (!projectionMIP.compare(t)) {
    projectionType = projectionMIP;
  } else if (!projectionSoftMIP.compare(t)) {
    projectionType = projectionSoftMIP;
    cmdstream.getline( t, TS, ':');
    softMipStrength = boost::lexical_cast< float >(t);
  } else {
    throw std::runtime_error(std::string("unknown projection:") + t);
  }
  cmdstream.getline( t, TS, ':');
  dimension = boost::lexical_cast< int >(t);

  cmdstream.getline( t, TS, ':');
  thickness = boost::lexical_cast< float >(t);

  cmdstream.getline( t, TS, ':');
  increment = boost::lexical_cast< float >(t);
}

template<class TImagePointerType, template < typename, typename > class Projector>
TImagePointerType CreateProjectedSlices(TImagePointerType image, int projectionDimension, double slice_thinkness, double slice_increment,
				     Projector<typename TImagePointerType::ObjectType::PixelType, typename TImagePointerType::ObjectType::PixelType> proj = 
					Projector<typename TImagePointerType::ObjectType::PixelType, typename TImagePointerType::ObjectType::PixelType>()) {
  typedef typename TImagePointerType::ObjectType TImageType;
  typedef itk::Image< typename TImageType::PixelType, 2> Image2DType;
  typedef itk::JoinSeriesImageFilter< Image2DType, TImageType > JoinSeriesFilterType;
  typedef itk::ExtractImageFilter< TImageType, Image2DType > Extract2DFilterType;
  typedef itk::ExtractImageFilter< TImageType, TImageType > Extract3DFilterType;
  typedef typename TImageType::SizeType::SizeValueType ImageSizeValueType;

  typename JoinSeriesFilterType::Pointer joiner = JoinSeriesFilterType::New();

  typename TImageType::RegionType region = image->GetLargestPossibleRegion();
  typename TImageType::RegionType sliceRegion = region;
  sliceRegion.SetSize(projectionDimension, 0 );
  double projectionSpacing = image->GetSpacing()[projectionDimension];
  double projectionsPos = 0.0;
  int maxSize = region.GetSize(projectionDimension);

  region.SetSize(projectionDimension, std::min( region.GetSize(projectionDimension), 
	  static_cast<ImageSizeValueType>(slice_thinkness / projectionSpacing )) );

  typename TImageType::DirectionType projectionDirection;
  while (static_cast<int>(projectionsPos / projectionSpacing) < maxSize) {
    region.SetIndex( projectionDimension, projectionsPos / projectionSpacing );
    region.SetSize(projectionDimension, std::min( region.GetSize(projectionDimension), 
	    static_cast<ImageSizeValueType>(maxSize - region.GetIndex(projectionDimension) )) );
    typename Extract3DFilterType::Pointer extractor = Extract3DFilterType::New();
    extractor->SetInput( image );
    extractor->SetExtractionRegion( region );
    extractor->Update();

    TImagePointerType extractImage = extractor->GetOutput();
    TImagePointerType projectionResult = ImageProjector(extractImage, proj, projectionDimension);
    projectionDirection = projectionResult->GetDirection();

    typename Extract2DFilterType::Pointer sliceExtractor = Extract2DFilterType::New();
    sliceExtractor->SetExtractionRegion( sliceRegion );
    sliceExtractor->SetInput( projectionResult );
    sliceExtractor->Update();
    joiner->PushBackInput( sliceExtractor->GetOutput() );
    projectionsPos += slice_increment;
  }
  joiner->Update();
  TImagePointerType result;
  result = joiner->GetOutput();
  result->SetOrigin( image->GetOrigin() );
  result->SetDirection( projectionDirection );
  typename TImageType::SpacingType resultSpacing = result->GetSpacing();
  resultSpacing[2] = slice_increment;
  result->SetSpacing( resultSpacing );
  return result;
}
