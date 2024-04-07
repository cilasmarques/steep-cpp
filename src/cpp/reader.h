#pragma once

#include "constants.h"

/**
 * @brief  Auxiliary struct to read data from TIFFs.
 */
struct Reader{
	uint16_t sampleFormat;
	uint8_t byteSize;
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
	Reader(uint16_t _sampleFormat, uint8_t _byteSize, tdata_t _buffer);

	/**
	 * @brief  Read the value of pixel in a TIFF.
	 * @param  column: Number of column of the pixel.
	 * @retval Value of the pixel.
	 */
	double read_tiff_pixel(uint32_t column);

  void check_open_tiff(TIFF *tif);
};
