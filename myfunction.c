// 208642884 Guy Adani
typedef struct {
   unsigned char red;
   unsigned char green;
   unsigned char blue;
} pixel;

#define MIN(a, b) (a < b ? a : b)
#define MAX(a, b) (a > b ? a : b)
#define CALC_INDEX(i, j, n) (i * n + j)

void myfunction(Image *image, char* srcImgpName, char* blurRsltImgName, char* sharpRsltImgName, char* filteredBlurRsltImgName, char* filteredSharpRsltImgName, char flag) {

    // these variables are very common, so we want them as registers for quick access.
    register int i;          // used in the loops.
    register int position;      // used to "move around" each pixel to sum.
    register int redSum, greenSum, blueSum; // used to sum.
    register pixel p;           // the current pixel/neighbour in the loops.
    register int mMinus1 = m - 1;
    register int index;         // save CALC_INDEX at the end of every loop
    register int im = m * m;    // save multiply every loop

    int size = im * sizeof(pixel); // calculating only once
    pixel* pixelsImg = malloc(size);
    pixel* backupOrg = malloc(size);

	if (flag == '1') { //                  @--------------- Filter OFF ---------------@

//                  --------------- blur image ---------------
//      we do doConvolution's job here instead of calling to function. this applies to all inside-functions used.

        // just like charsToPixels but faster
        memcpy(pixelsImg, image->data, size);
        // just like copyPixels but faster
        memcpy(backupOrg, pixelsImg, size);

        // This is the same as smooth func. Instead of using for loops, I use while, and the counting is going back,
        // which is faster. The loops sum the pixels 1,1 to m-1,m-1, so it isn't necessary to check the bounds.
        // We go "around" every pixel to sum it and its neighbours, it is done always 9 times so there are no loops inside.

        i = mMinus1;
        im = im - m - m;                // to start the loop right
        while (--i) {
            register int j = mMinus1;
            im -= m;
            while (--j) {

                // calculate current [0][0]
                position = im;
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
                position += mMinus1 - 1;
                p = backupOrg[position];
                redSum += p.red;
                greenSum += p.green;
                blueSum += p.blue;

                // backupOrg[1][1]
                index = ++position;                         // saving this index for the end
                p = backupOrg[index];
                redSum += p.red;
                greenSum += p.green;
                blueSum += p.blue;

                // backupOrg[1][2]
                p = backupOrg[++position];
                redSum += p.red;
                greenSum += p.green;
                blueSum += p.blue;

                // backupOrg[2][0]
                position += mMinus1 - 1;
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
                p.red = (unsigned char)(MIN(MAX(redSum / 9, 0), 255));
                p.green = (unsigned char)(MIN(MAX(greenSum / 9, 0), 255));
                p.blue = (unsigned char)(MIN(MAX(blueSum / 9, 0), 255));

                pixelsImg[index] = p;
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
        i = mMinus1;
        im = i * m - m;
        while (--i) {
            register int j = mMinus1;
            im -= m;
            while (--j) {

                // calculate current [0][0]
                position = im;
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
                position += mMinus1 - 1;
                p = backupOrg[position];
                redSum -= p.red;
                greenSum -= p.green;
                blueSum -= p.blue;

                // backupOrg[1][1]
                index = ++position;                         // saving this index for the end
                p = backupOrg[index];
                redSum += (p.red << 3) + p.red;             // using shift and
                greenSum += (p.green << 3) + p.green;       // add instead of
                blueSum += (p.blue << 3) + p.blue;          // multiply

                // backupOrg[1][2]
                p = backupOrg[++position];
                redSum -= p.red;
                greenSum -= p.green;
                blueSum -= p.blue;

                // backupOrg[2][0]
                position += mMinus1 - 1;
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
                p.red = (unsigned char)(MIN(MAX(redSum, 0), 255));
                p.green = (unsigned char)(MIN(MAX(greenSum, 0), 255));
                p.blue = (unsigned char)(MIN(MAX(blueSum, 0), 255));

                pixelsImg[index] = p;
            }
        }

        // just like pixelsToChar but faster
        memcpy(image->data, pixelsImg, size);

        // write result image to file
		writeBMP(image, srcImgpName, sharpRsltImgName);

	} else { //                        @--------------- Filter ON ---------------@

        // This section is very similar to the "filter off" section, except it apply the
        // extermum filtered kernel by finding the max and min intensity.
        // In the original applyKernel function it happens after blurring, but here it is
        // already known that the filter is on, so we can do it faster at the same time.

        // These variables are very common in this section, but haven't been used in the previous.
        register int maxInt, minInt; // used to find the max/min intensity.
        register int rgbSum;         // used to sum the RGB values for every pixel.

//                  --------------- blur & filter image ---------------

        // just like charsToPixels but faster
        memcpy(pixelsImg, image->data, size);
        // just like copyPixels but faster
        memcpy(backupOrg, pixelsImg, size);

        // As mentioned, these loops are similar to the normal blur from above, but also finding the
        // min/max intensity. I noticed the same sum of RGB happens 4 times every loop (worst case) so
        // here it is calculated once for every pixel, using the register variable.


        int minRow, minCol, maxRow, maxCol;     // used to find the min/max position

        i = mMinus1;
        im = i * m - m;
        while (--i) {               // i is the row
            register int j = mMinus1;
            im -= m;
            while (--j) {           // j is the column
                // Note: although the order is 0,0 -> 0,1 -> ... 2,2, i/j are the row/col number of 1,1 - current pixel.

                maxInt = -1; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
                minInt = 766;  // arbitrary value that is lower than minimum possible intensity, which is 0

                // calculate current [0][0]
                position = im;
                position += j - 1;

                // backupOrg[0][0]
                p = backupOrg[position];
                redSum = p.red;
                greenSum = p.green;
                blueSum = p.blue;

                // check for min/max intensity
                rgbSum = ((int)p.red) + ((int)p.green) + ((int)p.blue);
                if (rgbSum <= minInt) {
                    minInt = rgbSum;
                    minRow = i - 1;
                    minCol = j - 1;
                }
                if (rgbSum > maxInt) {
                    maxInt = rgbSum;
                    maxRow = i - 1;
                    maxCol = j - 1;
                }

                // backupOrg[0][1]
                p = backupOrg[++position];
                redSum += p.red;
                greenSum += p.green;
                blueSum += p.blue;

                // check for min/max intensity
                rgbSum = ((int)p.red) + ((int)p.green) + ((int)p.blue);
                if (rgbSum <= minInt) {
                    minInt = rgbSum;
                    minRow = i - 1;
                    minCol = j;
                } else if (rgbSum > maxInt) {
                    maxInt = rgbSum;
                    maxRow = i - 1;
                    maxCol = j;
                }

                // backupOrg[0][2]
                p = backupOrg[++position];
                redSum += p.red;
                greenSum += p.green;
                blueSum += p.blue;

                // check for min/max intensity
                rgbSum = ((int)p.red) + ((int)p.green) + ((int)p.blue);
                if (rgbSum <= minInt) {
                    minInt = rgbSum;
                    minRow = i - 1;
                    minCol = j + 1;
                } else if (rgbSum > maxInt) {
                    maxInt = rgbSum;
                    maxRow = i - 1;
                    maxCol = j + 1;
                }

                // backupOrg[1][0]
                position += mMinus1 - 1;
                p = backupOrg[position];
                redSum += p.red;
                greenSum += p.green;
                blueSum += p.blue;

                // check for min/max intensity
                rgbSum = ((int)p.red) + ((int)p.green) + ((int)p.blue);
                if (rgbSum <= minInt) {
                    minInt = rgbSum;
                    minRow = i;
                    minCol = j - 1;
                } else if (rgbSum > maxInt) {
                    maxInt = rgbSum;
                    maxRow = i;
                    maxCol = j - 1;
                }

                // backupOrg[1][1]
                index = ++position;                         // saving this index for the end
                p = backupOrg[index];
                redSum += p.red;
                greenSum += p.green;
                blueSum += p.blue;

                // check for min/max intensity
                rgbSum = ((int)p.red) + ((int)p.green) + ((int)p.blue);
                if (rgbSum <= minInt) {
                    minInt = rgbSum;
                    minRow = i;
                    minCol = j;
                } else if (rgbSum > maxInt) {
                    maxInt = rgbSum;
                    maxRow = i;
                    maxCol = j;
                }

                // backupOrg[1][2]
                p = backupOrg[++position];
                redSum += p.red;
                greenSum += p.green;
                blueSum += p.blue;

                // check for min/max intensity
                rgbSum = ((int)p.red) + ((int)p.green) + ((int)p.blue);
                if (rgbSum <= minInt) {
                    minInt = rgbSum;
                    minRow = i;
                    minCol = j + 1;
                } else if (rgbSum > maxInt) {
                    maxInt = rgbSum;
                    maxRow = i;
                    maxCol = j + 1;
                }

                // backupOrg[2][0]
                position += mMinus1 - 1;
                p = backupOrg[position];
                redSum += p.red;
                greenSum += p.green;
                blueSum += p.blue;

                // check for min/max intensity
                rgbSum = ((int)p.red) + ((int)p.green) + ((int)p.blue);
                if (rgbSum <= minInt) {
                    minInt = rgbSum;
                    minRow = i + 1;
                    minCol = j - 1;
                } else if (rgbSum > maxInt) {
                    maxInt = rgbSum;
                    maxRow = i + 1;
                    maxCol = j - 1;
                }

                // backupOrg[2][1]
                p = backupOrg[++position];
                redSum += p.red;
                greenSum += p.green;
                blueSum += p.blue;

                // check for min/max intensity
                rgbSum = ((int)p.red) + ((int)p.green) + ((int)p.blue);
                if (rgbSum <= minInt) {
                    minInt = rgbSum;
                    minRow = i + 1;
                    minCol = j;
                } else if (rgbSum > maxInt) {
                    maxInt = rgbSum;
                    maxRow = i + 1;
                    maxCol = j;
                }

                // backupOrg[2][2]
                p = backupOrg[++position];
                redSum += p.red;
                greenSum += p.green;
                blueSum += p.blue;

                // check for min/max intensity
                rgbSum = ((int)p.red) + ((int)p.green) + ((int)p.blue);
                if (rgbSum <= minInt) {
                    minInt = rgbSum;
                    minRow = i + 1;
                    minCol = j + 1;
                } else if (rgbSum > maxInt) {
                    maxInt = rgbSum;
                    maxRow = i + 1;
                    maxCol = j + 1;
                }

                // filter the min/max
                p = backupOrg[CALC_INDEX(minRow, minCol, m)];
                redSum -= ((int)p.red);
                greenSum -= ((int)p.green);
                blueSum -= ((int)p.blue);

                p = backupOrg[CALC_INDEX(maxRow, maxCol, m)];
                redSum -= ((int)p.red);
                greenSum -= ((int)p.green);
                blueSum -= ((int)p.blue);

                // The kernel scale is 7, so there's a division for each sum.
                p.red = (unsigned char)(MIN(MAX(redSum / 7, 0), 255));
                p.green = (unsigned char)(MIN(MAX(greenSum / 7, 0), 255));
                p.blue = (unsigned char)(MIN(MAX(blueSum / 7, 0), 255));

                pixelsImg[index] = p;
            }
        }

        // just like pixelsToChar but faster
        memcpy(image->data, pixelsImg, size);

		// write result image to file
		writeBMP(image, srcImgpName, filteredBlurRsltImgName);

//              --------------- sharpen the resulting image ---------------

//      Exactly the same as the "filter off" sharp!

        // just like charsToPixels but faster
        memcpy(pixelsImg, image->data, size);
        // just like copyPixels but faster
        memcpy(backupOrg, pixelsImg, size);

        // just like the smooth function from blur, but with the sharp matrix.
        i = mMinus1;
        im = i * m - m;
        while (--i) {
            register int j = mMinus1;
            im -= m;
            while (--j) {

                // calculate current [0][0]
                position = im;
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
                position += mMinus1 - 1;
                p = backupOrg[position];
                redSum -= p.red;
                greenSum -= p.green;
                blueSum -= p.blue;

                // backupOrg[1][1]
                index = ++position;                         // saving this index for the end
                p = backupOrg[index];
                redSum += (p.red << 3) + p.red;             // using shift and
                greenSum += (p.green << 3) + p.green;       // add instead of
                blueSum += (p.blue << 3) + p.blue;          // multiply

                // backupOrg[1][2]
                p = backupOrg[++position];
                redSum -= p.red;
                greenSum -= p.green;
                blueSum -= p.blue;

                // backupOrg[2][0]
                position += mMinus1 - 1;
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
                p.red = (unsigned char)(MIN(MAX(redSum, 0), 255));
                p.green = (unsigned char)(MIN(MAX(greenSum, 0), 255));
                p.blue = (unsigned char)(MIN(MAX(blueSum, 0), 255));

                pixelsImg[index] = p;
            }
        }

        // just like pixelsToChar but faster
        memcpy(image->data, pixelsImg, size);

		// write result image to file
		writeBMP(image, srcImgpName, filteredSharpRsltImgName);
	}

    free(pixelsImg);
    free(backupOrg);
}

