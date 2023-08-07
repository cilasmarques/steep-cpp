#include "tiffs.h"

PixelReader read_line_tiff(TIFF *tif, tdata_t tif_line, int line)
{
  uint32 width_band;
  uint16 sample_bands;

  TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &sample_bands);
  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width_band);

  tdata_t tif_buff = _TIFFmalloc(TIFFScanlineSize(tif));
  unsigned short tif_byte_size = TIFFScanlineSize(tif) / width_band;
  PixelReader tif_reader = PixelReader(sample_bands, tif_byte_size, tif_buff);

  if (TIFFReadScanline(tif, tif_buff, line) < 0)
  {
    cerr << "Read problem" << endl;
    exit(3);
  }

  return tif_reader;
};

/**
 * @brief  Verifies if a TIFF was open correctly.
 * @param  tif: TIFF to be verified
 * @throws Throw an error with exit code 1 if the TIFF isn't open.
 */
void check_open_tiff(TIFF *tif)
{
  if (!tif)
  {
    cerr << "Open tiff problem" << endl;
    exit(1);
  }
};

/**
 * @brief  Configures a TIFF based on a second TIFF.
 * @param  new_tif: TIFF to be configured.
 * @param  base_tif: TIFF used to provide the configurations.
 */
void setup(TIFF *new_tif, TIFF *base_tif)
{
  uint32 image_width, image_length;

  TIFFGetField(base_tif, TIFFTAG_IMAGEWIDTH, &image_width);
  TIFFGetField(base_tif, TIFFTAG_IMAGELENGTH, &image_length);

  TIFFSetField(new_tif, TIFFTAG_IMAGEWIDTH, image_width);
  TIFFSetField(new_tif, TIFFTAG_IMAGELENGTH, image_length);
  TIFFSetField(new_tif, TIFFTAG_BITSPERSAMPLE, 64);
  TIFFSetField(new_tif, TIFFTAG_SAMPLEFORMAT, 3);
  TIFFSetField(new_tif, TIFFTAG_COMPRESSION, 1);
  TIFFSetField(new_tif, TIFFTAG_PHOTOMETRIC, 1);
  TIFFSetField(new_tif, TIFFTAG_SAMPLESPERPIXEL, 1);
  TIFFSetField(new_tif, TIFFTAG_ROWSPERSTRIP, 1);
  TIFFSetField(new_tif, TIFFTAG_RESOLUTIONUNIT, 1);
  TIFFSetField(new_tif, TIFFTAG_XRESOLUTION, 1);
  TIFFSetField(new_tif, TIFFTAG_YRESOLUTION, 1);
  TIFFSetField(new_tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
};