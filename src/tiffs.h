#include "types.h"
#include "pixel_reader.h"

PixelReader read_line_tiff(TIFF *tif, tdata_t tif_line, int line);

void check_open_tiff(TIFF *tif);

void setup(TIFF *new_tif, TIFF *base_tif);
