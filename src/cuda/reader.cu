#include "reader.h"

Reader::Reader()
{
  sampleFormat = 0;
  byteSize = 0;
  buffer = NULL;
};

Reader::Reader(uint16_t _sampleFormat, uint8_t _byteSize, tdata_t _buffer)
{
  sampleFormat = _sampleFormat;
  byteSize = _byteSize;
  buffer = _buffer;
};

double Reader::read_tiff_pixel(uint32_t column)
{
  double ret = 0;
  switch (sampleFormat)
  {
  case 1:
  {
    uint64_t value = 0;
    memcpy(&value, static_cast<unsigned char *>(buffer) + (column * byteSize), byteSize);
    ret = value;
  }
  break;
  case 2:
  {
    int64_t value = 0;
    memcpy(&value, static_cast<unsigned char *>(buffer) + (column * byteSize), byteSize);
    ret = value;
  }
  break;
  case 3:
    switch (byteSize)
    {
    case 4:
    {
      float value = 0;
      memcpy(&value, static_cast<unsigned char *>(buffer) + (column * byteSize), byteSize);
      ret = value;
    }
    break;
    case 8:
    {
      double value = 0;
      memcpy(&value, static_cast<unsigned char *>(buffer) + (column * byteSize), byteSize);
      ret = value;
    }
    break;
    case 16:
    {
      long double value = 0;
      memcpy(&value, static_cast<unsigned char *>(buffer) + (column * byteSize), byteSize);
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
