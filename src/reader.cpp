#include "reader.h"

/**
 * @brief  Empty constructor.
 */
Reader::Reader()
{
  sampleFormat = 0;
  byteSize = 0;
  buffer = NULL;
};

/**
 * @brief  Constructor
 * @param  _sampleFormat: Sample format.
 * @param  _byteSize: Byte size.
 * @param  _buffer: Buffer.
 */
Reader::Reader(uint16 _sampleFormat, uint8 _byteSize, tdata_t _buffer)
{
  sampleFormat = _sampleFormat;
  byteSize = _byteSize;
  buffer = _buffer;
};

/**
 * @brief  Read the value of pixel in a TIFF.
 * @param  column: Number of column of the pixel.
 * @retval Value of the pixel.
 */
double Reader::read_tiff_pixel(uint32 column)
{
  double ret = 0;
  switch (sampleFormat)
  {
  case 1:
  {
    uint64 value = 0;
    memcpy(&value, buffer + (column * byteSize), byteSize);
    ret = value;
  }
  break;
  case 2:
  {
    int64 value = 0;
    memcpy(&value, buffer + (column * byteSize), byteSize);
    ret = value;
  }
  break;
  case 3:
    switch (byteSize)
    {
    case 4:
    {
      float value = 0;
      memcpy(&value, buffer + (column * byteSize), byteSize);
      ret = value;
    }
    break;
    case 8:
    {
      double value = 0;
      memcpy(&value, buffer + (column * byteSize), byteSize);
      ret = value;
    }
    break;
    case 16:
    {
      long double value = 0;
      memcpy(&value, buffer + (column * byteSize), byteSize);
      ret = value;
    }
    break;
    default:
      cerr << "Unsupported operation!" << endl;
      exit(7);
    }
    break;
  default:
    cerr << "Unsupported operation!" << endl;
    exit(7);
  }
  return ret;
};

void Reader::check_open_tiff(TIFF *tif)
{
  if (!tif)
  {
    cerr << "Open tiff problem" << endl;
    exit(1);
  }
};
