#include <stdio.h>
#include <stdlib.h>

int main() {

    // Global Variables
    unsigned int height, width, maxVal;

    // File Pointer
    FILE *fp;
    FILE *fp_out;
    
    const char *fileName="monalisa.ascii.pgm";   // Change file name as required
    const char *outfileName = "monalisa-reverse.pgm";
            
    fp = fopen(fileName, "r");

    // READING THE IMAGE FILE
    if (fp == NULL) {
        printf("File not found!");
        return -1;
    }
    else {
        char magic_type[3];
        fscanf(fp, "%s\n", magic_type);
        printf("Image MAGIC Type: %s\n", magic_type);

        if ( magic_type[0] == 'P' & magic_type[1] == '2' ) {
            // Check the File type MAGIC
            char ch[255]; // For storing comments
            
            while ((ungetc(getc(fp), fp)) == '#') {
                fgets(ch, 255, fp);
            }

            fscanf(fp, "%d %d", &height, &width);

            printf("Height: %d pixels\n", height);
            printf("Width: %d pixels\n", width);

            // Allocate memory for Image 2D array
            unsigned int imgArray[height][width];

            // Image Type (Binary or Grayscale)
            fscanf(fp, "%d", &maxVal);
            if (maxVal == 255) {
                printf("Grayscale Image\n");
            } else if (maxVal == 1) {
                printf("Binary Image\n");
            }

            // Image Data
            unsigned int imgData;
            for (int i = 0; i < height; i++) {
                for (int j = 0; j < width; j++) {
                    fscanf(fp, "%d", &imgData);
                    imgArray[i][j] = imgData;
                }
            }

            
            // Invert and write in a file
            fp_out = fopen(outfileName, "w");

            fprintf(fp_out, "P2\n%d %d\n%d\n", height, width, maxVal);

            for (int i = 0; i < height; i++) {
                for (int j = 0; j < width; j++) {
                    fprintf(fp_out, "%d ", maxVal - imgArray[i][j]);
                }
                fprintf(fp_out, "\n");
            }

            
        } else {
            printf("Not Supported File type");
        }
    }

    fclose(fp); // Close the file and free the pointer
    fclose(fp_out); // Close the file and free the pointer

    return 0;
}