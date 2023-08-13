#pragma once

#include "constants.h"

/**
 * @brief  Auxiliary struct to read data from TIFFs.
 */
struct Reader{
	uint16 sampleFormat;
	uint8 byteSize;
	tdata_t buffer;

	/**
	 * @brief  Empty constructor.
	 */
	Reader();

	/**
	 * @brief  Constructor
	 * @param  _sampleFormat: Sample format.
	 * @param  _byteSize: Byte size.
	 * @param  _buffer: Buffer.
	 */
	Reader(uint16 _sampleFormat, uint8 _byteSize, tdata_t _buffer);

	/**
	 * @brief  Read the value of pixel in a TIFF.
	 * @param  column: Number of column of the pixel.
	 * @retval Value of the pixel.
	 */
	double read_tiff_pixel(uint32 column);

  void check_open_tiff(TIFF *tif);
};
