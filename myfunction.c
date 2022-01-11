#include <stdbool.h>

typedef struct {
   unsigned char red;
   unsigned char green;
   unsigned char blue;
//   unsigned char c1;
//   unsigned short s1, s2;
//   unsigned int i1;
} pixel;

typedef struct {
    int red;
    int green;
    int blue;
    //int num;
} pixel_sum;

#define MIN(a, b) (a < b ? a : b)
#define MAX(a, b) (a > b ? a : b)
#define CALC_INDEX(i, j, n) (i * n + j)


/* Compute min and max of two integers, respectively */
int min(int a, int b) { return (a < b ? a : b); }
int max(int a, int b) { return (a > b ? a : b); }

int calcIndex(int i, int j, int n) {
	return ((i)*(n)+(j));
}

/*
 * initialize_pixel_sum - Initializes all fields of sum to 0
 */
void initialize_pixel_sum(pixel_sum *sum) {
	sum->red = sum->green = sum->blue = 0;
	// sum->num = 0;
	return;
}

/*
 * assign_sum_to_pixel - Truncates pixel's new value to match the range [0,255]
 */
static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum, int kernelScale) {

	// divide by kernel's weight
	sum.red = sum.red / kernelScale;
	sum.green = sum.green / kernelScale;
	sum.blue = sum.blue / kernelScale;

	// truncate each pixel's color values to match the range [0,255]
	current_pixel->red = (unsigned char) (min(max(sum.red, 0), 255));
	current_pixel->green = (unsigned char) (min(max(sum.green, 0), 255));
	current_pixel->blue = (unsigned char) (min(max(sum.blue, 0), 255));
	return;
}

/*
* sum_pixels_by_weight - Sums pixel values, scaled by given weight
*/
static void sum_pixels_by_weight(pixel_sum *sum, pixel p, int weight) {
	sum->red += ((int) p.red) * weight;
	sum->green += ((int) p.green) * weight;
	sum->blue += ((int) p.blue) * weight;
	// sum->num++;
	return;
}

/*
 *  Applies kernel for pixel at (i,j)
 */
static pixel applyKernel(int dim, int i, int j, pixel *src, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

	int ii, jj;
	int currRow, currCol;
	pixel_sum sum;
	pixel current_pixel;
	int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
	int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
	int min_row, min_col, max_row, max_col;
	pixel loop_pixel;

    sum.red = sum.green = sum.blue = 0;

	for(ii = max(i-1, 0); ii <= min(i+1, dim-1); ii++) {
		for(jj = max(j-1, 0); jj <= min(j+1, dim-1); jj++) {

			int kRow, kCol;

			// compute row index in kernel
			if (ii < i) {
				kRow = 0;
			} else if (ii > i) {
				kRow = 2;
			} else {
				kRow = 1;
			}

			// compute column index in kernel
			if (jj < j) {
				kCol = 0;
			} else if (jj > j) {
				kCol = 2;
			} else {
				kCol = 1;
			}

			// apply kernel on pixel at [ii,jj]
			sum_pixels_by_weight(&sum, src[calcIndex(ii, jj, dim)], kernel[kRow][kCol]);
		}
	}

	if (filter) {
		// find min and max coordinates
		for(ii = max(i-1, 0); ii <= min(i+1, dim-1); ii++) {
			for(jj = max(j-1, 0); jj <= min(j+1, dim-1); jj++) {
				// check if smaller than min or higher than max and update
				loop_pixel = src[calcIndex(ii, jj, dim)];
				if ((((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue)) <= min_intensity) {
					min_intensity = (((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue));
					min_row = ii;
					min_col = jj;
				}
				if ((((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue)) > max_intensity) {
					max_intensity = (((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue));
					max_row = ii;
					max_col = jj;
				}
			}
		}
		// filter out min and max
		sum_pixels_by_weight(&sum, src[calcIndex(min_row, min_col, dim)], -1);
		sum_pixels_by_weight(&sum, src[calcIndex(max_row, max_col, dim)], -1);
	}

	// assign kernel's result to pixel at [i,j]
	assign_sum_to_pixel(&current_pixel, sum, kernelScale);
	return current_pixel;
}


void myfunction(Image *image, char* srcImgpName, char* blurRsltImgName, char* sharpRsltImgName, char* filteredBlurRsltImgName, char* filteredSharpRsltImgName, char flag) {

    // these variables are very common, so we want them as registers for quick access.
    register int i, j;          // used in the loops.
    register int position;      // used to "move around" each pixel to sum.
    register int redSum, greenSum, blueSum; // used to sum. it avoids a lot of memory calls of the struct.
    register pixel curPixel;    // the current pixel to calculate.
    register pixel p;           // the current neighbour of curPixel.

    /*
    * [1, 1, 1]
    * [1, 1, 1]
    * [1, 1, 1]
    */
	int blurKernel[3][3] = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

	/*
	* [-1, -1, -1]
	* [-1, 9, -1]
	* [-1, -1, -1]
	*/
	int sharpKernel[3][3] = {{-1,-1,-1},{-1,9,-1},{-1,-1,-1}};
    int size = m * m * sizeof(pixel); // calculating only once
    pixel* pixelsImg = malloc(size);
    pixel* backupOrg = malloc(size);

	if (flag == '1') {
//                  --------------- blur image ---------------
//      we do doConvolution's job here instead of calling to function. this applies to all inside-functions used.

        // just like charsToPixels but faster
        memcpy(pixelsImg, image->data, size);
        // just like copyPixels but faster
        memcpy(backupOrg, pixelsImg, size);

        // This is the same as smooth func. Instead of using for loops, I use while, and the counting is going back,
        // which is faster. The loops sum the pixels 1,1 to m-1,m-1, so it isn't necessary to check the bounds.
        // We go "around" every pixel to sum it and its neighbours, it is done always 9 times so there are no loops inside.
        i = m - 1;
        while (--i) {
            j = m - 1;
            while (--j) {

                // calculate current [0][0]
                position = (i - 1) * m;
                position += j - 1;

                // backupOrg[0][0]
                p = backupOrg[position];
                redSum = p.red;
                greenSum = p.green;
                blueSum = p.blue;

                // backupOrg[0][1]
                p = backupOrg[++position];
                redSum += p.red;
                greenSum += p.green;
                blueSum += p.blue;

                // backupOrg[0][2]
                p = backupOrg[++position];
                redSum += p.red;
                greenSum += p.green;
                blueSum += p.blue;

                // backupOrg[1][0]
                position += (m - 2);
                p = backupOrg[position];
                redSum += p.red;
                greenSum += p.green;
                blueSum += p.blue;

                // backupOrg[1][1]
                p = backupOrg[++position];
                redSum += p.red;
                greenSum += p.green;
                blueSum += p.blue;

                // backupOrg[1][2]
                p = backupOrg[++position];
                redSum += p.red;
                greenSum += p.green;
                blueSum += p.blue;

                // backupOrg[2][0]
                position += (m - 2);
                p = backupOrg[position];
                redSum += p.red;
                greenSum += p.green;
                blueSum += p.blue;

                // backupOrg[2][1]
                p = backupOrg[++position];
                redSum += p.red;
                greenSum += p.green;
                blueSum += p.blue;

                // backupOrg[2][2]
                p = backupOrg[++position];
                redSum += p.red;
                greenSum += p.green;
                blueSum += p.blue;

                // The kernel scale is 9, so there's a division for each sum.
                curPixel.red = (unsigned char)(MIN(MAX(redSum / 9, 0), 255));
                curPixel.green = (unsigned char)(MIN(MAX(greenSum / 9, 0), 255));
                curPixel.blue = (unsigned char)(MIN(MAX(blueSum / 9, 0), 255));

                pixelsImg[CALC_INDEX(i, j, m)] = curPixel;
            }
        }

        // just like pixelsToChars but faster
        memcpy(image->data, pixelsImg, size);

		// write result image to file
		writeBMP(image, srcImgpName, blurRsltImgName);

//              --------------- sharpen the resulting image ---------------

//      we do doConvolution's job here instead of calling to function. this applies to all inside-functions used.

        // just like charsToPixels but faster
        memcpy(pixelsImg, image->data, size);
        // just like copyPixels but faster
        memcpy(backupOrg, pixelsImg, size);

        // just like the smooth function from blur, but with the sharp matrix.
        i = m - 1;
        while (--i) {
            j = m - 1;
            while (--j) {

                // calculate current [0][0]
                position = (i - 1) * m;
                position += j - 1;

                // backupOrg[0][0]
                p = backupOrg[position];
                redSum = - p.red;
                greenSum = - p.green;
                blueSum = - p.blue;

                // backupOrg[0][1]
                p = backupOrg[++position];
                redSum -= p.red;
                greenSum -= p.green;
                blueSum -= p.blue;

                // backupOrg[0][2]
                p = backupOrg[++position];
                redSum -= p.red;
                greenSum -= p.green;
                blueSum -= p.blue;

                // backupOrg[1][0]
                position += (m - 2);
                p = backupOrg[position];
                redSum -= p.red;
                greenSum -= p.green;
                blueSum -= p.blue;

                // backupOrg[1][1]
                p = backupOrg[++position];
                redSum += 9 * p.red;
                greenSum += 9 * p.green;
                blueSum += 9 * p.blue;

                // backupOrg[1][2]
                p = backupOrg[++position];
                redSum -= p.red;
                greenSum -= p.green;
                blueSum -= p.blue;

                // backupOrg[2][0]
                position += (m - 2);
                p = backupOrg[position];
                redSum -= p.red;
                greenSum -= p.green;
                blueSum -= p.blue;

                // backupOrg[2][1]
                p = backupOrg[++position];
                redSum -= p.red;
                greenSum -= p.green;
                blueSum -= p.blue;

                // backupOrg[2][2]
                p = backupOrg[++position];
                redSum -= p.red;
                greenSum -= p.green;
                blueSum -= p.blue;

                // Here the scale is 1 so no need to divide.
                curPixel.red = (unsigned char)(MIN(MAX(redSum, 0), 255));
                curPixel.green = (unsigned char)(MIN(MAX(greenSum, 0), 255));
                curPixel.blue = (unsigned char)(MIN(MAX(blueSum, 0), 255));

                pixelsImg[CALC_INDEX(i, j, m)] = curPixel;
            }
        }

        // just like pixelsToChar but faster
        memcpy(image->data, pixelsImg, size);

        // write result image to file
		writeBMP(image, srcImgpName, sharpRsltImgName);

	} else {
		// apply extermum filtered kernel to blur image
		//doConvolution(image, 3, blurKernel, 7, true);
        // memcpy is faster than charsToPixels
        memcpy(pixelsImg, image->data, size);
        // memcpy is faster than copyPixels
        memcpy(backupOrg, pixelsImg, size);
// V

        //smooth(m, backupOrg, pixelsImg, 3, blurKernel, 9, false);
        int i, j;
        for (i = 1; i < m - 1; i++) {
            for (j =  1; j < m - 1; j++) {
                pixelsImg[CALC_INDEX(i, j, m)] = applyKernel(m, i, j, backupOrg, 3, blurKernel, 7, true); //                       change applyKernel
            }
        }

        //pixelsToChars(pixelsImg, image);
        // memcpy is faster than pixelsToChar
        memcpy(image->data, pixelsImg, size);
		// write result image to file
		writeBMP(image, srcImgpName, filteredBlurRsltImgName);

		// sharpen the resulting image

		//doConvolution(image, 3, sharpKernel, 1, false);
// memcpy is faster than charsToPixels
        memcpy(pixelsImg, image->data, size);

        // memcpy is faster than copyPixels
        memcpy(backupOrg, pixelsImg, size);

        //smooth(m, backupOrg, pixelsImg, 3, sharpKernel, 1, false);
        for (i = 1; i < m - 1; i++) {
            for (j =  1; j < m - 1; j++) {
                pixelsImg[CALC_INDEX(i, j, m)] = applyKernel(m, i, j, backupOrg, 3, sharpKernel, 1, false); //                       change applyKernel
            }
        }

        //pixelsToChars(pixelsImg, image);
        // memcpy is faster than pixelsToChar
        memcpy(image->data, pixelsImg, size);
		// write result image to file
		writeBMP(image, srcImgpName, filteredSharpRsltImgName);
	}

    free(pixelsImg);
    free(backupOrg);
}

