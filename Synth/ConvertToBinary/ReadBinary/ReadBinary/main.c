//
//  main.c
//  ReadBinary
//
//  Created by 방정호 on 3/15/17.
//  Copyright © 2017 Bangtoven. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>

typedef struct Model {
    double frequency;
    uint32_t partials;
    uint32_t length;
    uint64_t empty[6]; // for later use...
} Model;

int main(int argc, const char * argv[]) {
    // insert code here...
    FILE *file = fopen("/Users/bangtoven/Desktop/Bassoon.ff.C2B2-2-70.msm","rb");
    if (!file) {
        printf("Unable to open file!");
        return 1;
    }
    
    Model model;
    fread(&model, sizeof(Model), 1, file);
    printf("fundamental freq: %lf\n", model.frequency);
    printf("size: %d %d\n", model.partials, model.length);
    
    int cells = model.partials * model.length;
    double *mat = malloc(cells * sizeof(double));
    fread(mat, sizeof(double), cells, file);
    printf("For t=0, \n");
    for(int i=0; i<model.partials; i++) {
        printf("%lf ", mat[i]);
    }
    printf("\n");
    
    fclose(file);
    
    
    return 0;
}
