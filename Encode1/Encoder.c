
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <windows.h>
#include <math.h>

BITMAPFILEHEADER    Header1;
BITMAPINFOHEADER    Header2;




#define mm  4           /* RS code over GF(2**4) - change to suit:(x^0 x^1 x^2 x^3 x^4 -> 2^4 =x^4 +x+1=11001)*/
#define nn  15          /* nn=2**mm -1   length of codeword  */
#define tt  3           /* number of errors that can be corrected */
#define kk  9           /* kk = nn-2*tt  */

int pp[mm + 1] = { 1, 1, 0, 0, 1 }; /* specify irreducible polynomial coeffts  פולינום הגנרי 4*/
int alpha_to[nn + 1], index_of[nn + 1], gg[nn - kk + 1];
int recd[nn] ,data[kk + 2 + kk], bb[nn - kk], currentData[kk];
int saveIndexPading = 0, lenght = 0, orginal[kk + 2 + kk], afterChange[nn + nn], t;
int numOfAllChunks = 0;


void generate_gf()
/* generate GF(2**mm) from the irreducible polynomial p(X) in pp[0]..pp[mm]
   lookup tables:  index->polynomial form   alpha_to[] contains j=alpha**i;
                   polynomial form -> index form  index_of[j=alpha**i] = i
   alpha=2 is the primitive element of GF(2**mm)
*/
{
    register int i, mask;

    mask = 1;
    alpha_to[mm] = 0;
    for (i = 0; i < mm; i++)
    {
        alpha_to[i] = mask;
        index_of[alpha_to[i]] = i;
        if (pp[i] != 0)
            alpha_to[mm] ^= mask;
        mask <<= 1;
    }
    index_of[alpha_to[mm]] = mm;
    mask >>= 1;
    for (i = mm + 1; i < nn; i++)
    {
        if (alpha_to[i - 1] >= mask)
            alpha_to[i] = alpha_to[mm] ^ ((alpha_to[i - 1] ^ mask) << 1);
        else alpha_to[i] = alpha_to[i - 1] << 1;
        index_of[alpha_to[i]] = i;
    }
    index_of[0] = -1;
}

void gen_poly()
/* Obtain the generator polynomial of the tt-error correcting, length
  nn=(2**mm -1) Reed Solomon code  from the product of (X+alpha**i), i=1..2*tt
*/
{
    register int i, j;

    gg[0] = 2;    /* primitive element alpha = 2  for GF(2**mm)  */
    gg[1] = 1;    /* g(x) = (X+alpha) initially */
    for (i = 2; i <= nn - kk; i++)
    {
        gg[i] = 1;
        for (j = i - 1; j > 0; j--)
            if (gg[j] != 0)  gg[j] = gg[j - 1] ^ alpha_to[(index_of[gg[j]] + i) % nn];
            else gg[j] = gg[j - 1];
        gg[0] = alpha_to[(index_of[gg[0]] + i) % nn];     /* gg[0] can never be zero */
    }
    /* convert gg[] to index form for quicker encoding */
    for (i = 0; i <= nn - kk; i++)  gg[i] = index_of[gg[i]];
}

void encode_rs()
/* take the string of symbols in data[i], i=0..(k-1) and encode systematically
   to produce 2*tt parity symbols in bb[0]..bb[2*tt-1]
   data[] is input and bb[] is output in polynomial form.
   Encoding is done by using a feedback shift register with appropriate
   connections specified by the elements of gg[], which was generated above.
   Codeword is   c(X) = data(X)*X**(nn-kk)+ b(X)          */
{
    register int i, j;
    int feedback;

    for (i = 0; i < nn - kk; i++)   bb[i] = 0;
    for (i = kk - 1; i >= 0; i--)
    {
        feedback = index_of[currentData[i] ^ bb[nn - kk - 1]];
        if (feedback != -1)
        {
            for (j = nn - kk - 1; j > 0; j--)
                if (gg[j] != -1)
                    bb[j] = bb[j - 1] ^ alpha_to[(gg[j] + feedback) % nn];
                else
                    bb[j] = bb[j - 1];
            bb[0] = alpha_to[(gg[0] + feedback) % nn];
        }
        else
        {
            for (j = nn - kk - 1; j > 0; j--)
                bb[j] = bb[j - 1];
            bb[0] = 0;
        };
    };
};

void  writeOrginalData() {
    FILE* fp;
    fopen_s(&fp, "OrginalData.bin", "wb");
    for (int i = 0; i <= lenght; i++) {
        fwrite(&data[i], sizeof(data[i]), 1, fp);
    }
    fclose(fp);
}

void readOrginalData() {
    FILE* fp;
    fopen_s(&fp, "OrginalData.bin", "rb");
    for (size_t i = 0; i <= sizeof(orginal) / sizeof(int); i++)
    {
        fread(&orginal[i], sizeof(int), 1, fp);
        printf("%d ", orginal[i]);
    }
    fclose(fp);
}

void noiseActivation()
{
    recd[nn / 2] = 4;
    recd[0] = 6;
    recd[3]=1;
}

void  writeIndexPading() {
    FILE* fp;
    fopen_s(&fp, "IndexPading.bin", "wb");
    fwrite(&saveIndexPading, sizeof(int), 1, fp);
    fclose(fp);
}

void saveDataAndRedundancy() {
    int i = 0;

    for (i = 0; i < nn - kk; i++)  recd[i] = bb[i];
    for (i = 0; i < kk; i++) recd[i + nn - kk] = currentData[i];
    printf("\n rect after encoder:\n");
    for (i = 0; i < nn; i++) {
        /* recd[i] = index_of[recd[i]];*/
        printf("%d  ", recd[i]);
    }
    printf("\n");
}

void writeDataAfterChange() {

    FILE* fw;
    fopen_s(&fw, "afterChange.bin", "ab");
    /* if (fw == NULL)
         printf("\n open file erro\n");*/
    printf("\n write Data AfterChange: \n");
    for (int i = 0; i < nn; i++) {
        fwrite(&recd[i], sizeof(int), 1, fw);
        printf("%d ", recd[i]);
    }
    fclose(fw);
}

void pading(int IndexStartPading) {

    for (int i = IndexStartPading; i < kk; i++)
    {

        currentData[i] = 0;
    }
}

void writenNumOfAllChunks() {
    FILE* fp;
    fopen_s(&fp, "numOfAllChunks.bin", "wb");


    fwrite(&numOfAllChunks, sizeof(int), 1, fp);
    fclose(fp);
}

void writeArrGenerate_gf() {
    FILE* fp;
    printf("\nWe will write these arr ​​to a binary file for Decode \n");

    fopen_s(&fp, "ArrGenerate.bin", "wb");
    printf("\n index_of: \n");
    for (int i = 0; i < nn + 1; i++) {
        fwrite(&index_of[i], sizeof(int), 1, fp);
        printf("%d ", index_of[i]);
    }
    printf("\n alpha_to: \n");
    for (int i = 0; i < nn + 1; i++) {
        fwrite(&alpha_to[i], sizeof(int), 1, fp);
        printf("%d ", alpha_to[i]);
    }
    printf("\n The generic polynomial: \n");
    for (int i = 0; i < nn - kk + 1; i++) {
        fwrite(&gg[i], sizeof(int), 1, fp);
        printf("%d ", gg[i]);
    }
    fclose(fp);
}

void deleteBinFile() {

    if (remove("OrginalData.bin") == 0) {
        printf("The file is deleted successfully.");
    }
    else {
        printf("The file is not deleted.");
    }

    if (remove("IndexPading.bin") == 0) {
        printf("The file is deleted successfully.");
    }
    else {
        printf("The file is not deleted.");
    }

    if (remove("afterChange.bin") == 0) {
        printf("The file is deleted successfully.");
    }
    else {
        printf("The file is not deleted.");
    }

}

main()
{
    int i=0,j = 0, k = 0;
    int saveIndexPading = 0;
#pragma region image
    /* FILE* myPalette;
      uint8_t* str = (uint8_t*)malloc(50 * sizeof(uint8_t));


      FILE* BMP;
      fopen_s(&BMP, "C:/Users/pc/Pictures/bird.bmp", "rb");

      fread(&Header1, sizeof(Header1), 1, BMP);

      printf("Type:%hd  and Type in  %x\n", Header1.bfType, Header1.bfType);
      printf("bfSize %d\n", Header1.bfSize);
      printf("bfOffBits %d\n", Header1.bfOffBits);

      fread(&Header2, sizeof(Header2), 1, BMP);
      int img_size = (Header2.biWidth * Header2.biHeight);
      printf("img size =  %d\n", img_size);
      uint8_t* data = (uint8_t*)malloc(img_size * sizeof(uint8_t));
      fread(&myPalette, 4, 2, BMP);

      int x, y;
      uint8_t pixel;
      printf("the array is: \n");
      for (int i = 0; i < img_size; i++) {
          fread(&pixel, sizeof(pixel), 1, BMP);
          data[i] = pixel;
          printf("{%d} ", data[i]);
      }*/
      ////////////////////////////////////////////////////////////////////////////////////
#pragma endregion  
   // data={ 8,6,8,1,2,4,8,9,9,9,7,7,7,1,2,3,4,5,6,7,8 };
#pragma region data arr
    data[0] = 8;
    data[1] = 6;
    data[2] = 8;
    data[3] = 1;
    data[4] = 2;
    data[5] = 4;
    data[6] = 8;
    data[7] = 9;
    data[8] = 9;
    data[9] = 9;
    data[10] = 7;
    data[11] = 7;
    data[12] = 7;
    data[13] = 1;
    data[14] = 2;
    data[15] = 3;
    data[16] = 4;
    data[17] = 5;
    data[18] = 6;
    data[19] = 7;
    data[20] = 8;
#pragma endregion
   
    deleteBinFile();
    generate_gf();
    
    

    printf("Look-up tables for GF(2**%2d)\n", mm);
    printf("  i   alpha_to[i]  index_of[i]\n");
    for (int i = 0; i <= nn; i++)
        printf("%3d      %3d          %3d\n", i, alpha_to[i], index_of[i]);
    printf("\n\n");
    
    gen_poly();
  
    lenght = sizeof(data) / sizeof(int);
    printf("\nlenght Arr Bmp = %d \n",lenght);
    int numOfChunks = lenght / 9;
    int carry = lenght % 9;
    printf("\ncarry: %d   \n", carry);
    printf("\nnumOfChunks: %d   \n", numOfChunks);
    writeOrginalData();
    printf("\n Orginal data read from file :\n");
    readOrginalData();
    while (numOfChunks > 0)
    {
        for (j = 0; j < 9; j++) {
            currentData[j] = data[k];
            k++;
        }
        printf("\n The current data is for a size 9 shank:\n");
        for (int i = 0; i < 9; i++)
        {
            printf("%d ", currentData[i]);
        }

        encode_rs();

        printf("\nThe redundancy:\n");
        for (int i = 0; i < 6; i++)
        {
            printf("%d ", bb[i]);
        }

        saveDataAndRedundancy();
        numOfChunks--;

        noiseActivation();
        printf("\n The current data after the noise:\n");
        for (int i = 0; i < 15; i++)
        {
            printf("%d ", recd[i]);
        }
        numOfAllChunks++;
        writeDataAfterChange();
    }
    
    if (carry != 0)
    {

        while (k <= lenght) {
            currentData[saveIndexPading] = data[k];
            k++;
            saveIndexPading++;
        }

     
        writeIndexPading(saveIndexPading);
        pading(saveIndexPading);

        printf("\n The current data after pading\n");
        for (size_t i = 0; i < kk; i++)
        {
            printf("%d ", currentData[i]);
        }

        encode_rs();

        printf("\nThe redundancy:\n");
        for (int i = 0; i < 6; i++)
        {
            printf("%d ", bb[i]);
        }
        saveDataAndRedundancy();

        noiseActivation();
        printf("\n The current data after the noise::\n");
        for (int i = 0; i < 15; i++)
        {
            printf("%d ", recd[i]);
        }
        numOfAllChunks++;
        writeDataAfterChange();
      
    }
   /* in order to perform Decod*/
    writenNumOfAllChunks();
    writeArrGenerate_gf();

}
