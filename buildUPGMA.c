/**
 * Author: Alemdar Salmoor
 * Date: December 14, 2019
**/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

struct TreeNode
{
    double height;
    char * name;
    int namelen;
    int numOfGenes;
};


int NuctoDec(char c){
    if (c == 'A'){ return 0; }
    else if(c == 'C'){ return 1; }
    else if(c == 'G'){ return 2; }
    else { return 3; }
}

int max(int A, int B){
    if(A > B){
        return A;
    }
    else
    {
        return B;
    }
}
void findIndexOfMin( double ( * tMatrix)[25], int nodeSize, int * a,int *b){
    int min = INT_MAX;

    int first = 0;
    int second = first + 1;
    while (first < nodeSize)
    {
        second = first + 1;
        while (second < nodeSize)
        {
            if(tMatrix[first][second] <= min){
                min = tMatrix[first][second];
                (*a) = first; (*b) = second;
            }

            second++;
        }
        first++;
        
    }
    
}

void computeAGlobal(char ** seq, int * seqSize, int (* sMatrix)[4], int gapopen, int gapext, char *** seqOutReverse, int * seqOutLen, int * score){

    int ** V = malloc((seqSize[0]+1) * sizeof(int *));
    char ** P = malloc((seqSize[0]+1) * sizeof(char *));
    int ** G = malloc((seqSize[0]+1) * sizeof(int *));


    char letter;

    *seqOutReverse = malloc(2 * sizeof(char *));
    int seqOutCapacity = max(seqSize[0], seqSize[1]);
    int seqOutSize = 0;
    (*seqOutReverse)[0] = malloc(seqOutCapacity * sizeof(char));
    (*seqOutReverse)[1] = malloc(seqOutCapacity * sizeof(char));

    //sup grid
    int i = 0, j;
    int ** F = malloc((seqSize[0] + 1) * sizeof(int *));
    int ** E = malloc((seqSize[0] + 1) * sizeof(int *));


    //initialize
    i = 0;
    while (i < seqSize[0] + 1)
    {
        V[i] = malloc((seqSize[1]+1) * sizeof(int));
        P[i] = malloc((seqSize[1]+1) * sizeof(char));
        E[i] = malloc((seqSize[1]+1) * sizeof(int));
        F[i] = malloc((seqSize[1]+1) * sizeof(int));
        G[i] = malloc((seqSize[1]+1) * sizeof(int));
        

        j = 0;
        while (j < seqSize[1] + 1)
        {
            if(i == 0 && j == 0){
                
                E[i][j] = INT_MIN/2;
                F[i][j] = INT_MIN/2;
                G[i][j] = INT_MIN/2;
                V[i][j] = 0;
                P[i][j] = 't';
            }

            if(i == 0 && j != 0){
                
                E[i][j] = 0;
                F[i][j] = gapopen + j * gapext;
                G[i][j] = 0;
                V[i][j] = gapopen + j * gapext;
                P[i][j] = 'l';

            }

            if(j == 0 && i != 0){

                E[i][j] = gapopen + i * gapext;
                F[i][j] = 0;
                G[i][j] = 0;
                V[i][j] = gapopen + i * gapext;
                P[i][j] = 'u';

            }
            
            j++;
        }
        
        i++;
    }

    //fill in
    i = 1;
    while (i < seqSize[0] + 1)
    {

        j = 1;
        while (j < seqSize[1] + 1)
        {

            E[i][j] = max(E[i][j-1] + gapext, (V[i][j-1] + gapopen + gapext));
            F[i][j] = max(F[i-1][j] + gapext, (V[i-1][j] + gapopen + gapext));
            G[i][j] = V[i-1][j-1] + sMatrix[NuctoDec(seq[0][i-1])][NuctoDec(seq[1][j-1])];

        
            if(F[i][j] >= G[i][j] && F[i][j] >= E[i][j]){ V[i][j] = F[i][j]; P[i][j] = 'u';}
            else if(E[i][j] >= G[i][j] && E[i][j] >= F[i][j]){ V[i][j] = E[i][j]; P[i][j] = 'l'; }
            else{ V[i][j] = G[i][j]; P[i][j] = 'd'; }

            j++;
        }
        
        i++;
    }

    //find
    letter = P[seqSize[0]][seqSize[1]];

    i = seqSize[0];
    j = seqSize[1];

    while (letter != 't')
    {
        if(seqOutSize >= seqOutCapacity)
        {
            (*seqOutReverse)[0] = realloc((*seqOutReverse)[0], seqOutCapacity * 2);
            (*seqOutReverse)[1] = realloc((*seqOutReverse)[1], seqOutCapacity * 2);
            seqOutCapacity = seqOutCapacity * 2;
        }

        if(letter == 'd'){
            (*seqOutReverse)[0][seqOutSize] = seq[0][i-1];
            (*seqOutReverse)[1][seqOutSize] = seq[1][j-1];

            i--;
            j--;
        }
        else if(letter == 'l')
        {
            (*seqOutReverse)[0][seqOutSize] = '-';
            (*seqOutReverse)[1][seqOutSize] = seq[1][j-1];
            j--;
        }
        else
        {
            (*seqOutReverse)[0][seqOutSize] = seq[0][i-1];
            (*seqOutReverse)[1][seqOutSize] = '-';
            i--;
        }

        seqOutSize++;

        letter = P[i][j];
        
    }
    

    *score = V[seqSize[0]][seqSize[1]];
    *seqOutLen = seqOutSize;

    i = 0;
    while (i < seqSize[0]+1)
    {
        free(P[i]);
        free(V[i]);
        free(E[i]);
        free(F[i]);
        free(G[i]);
        i++;
    }
    
    free(V);
    free(P);
    free(E);
    free(F);
    free(G);

}


int main(int argc, char** argv){

    char * seq = argv[2];
    int match = atoi(argv[4]);
    int mismatch = atoi(argv[6]);
    int gapopen = atoi(argv[8]);
    int gapext = atoi(argv[10]);
    char * out = argv[12];

    //Create the sMatrix
    int sMatrix[4][4];

    int i = 0, j;

    while (i < 4)
    {
        j = 0;
        while (j < 4)
        {
            if(i == j){ sMatrix[i][j] = match; }
            else{ sMatrix[i][j] = mismatch; }
            
            j++;
            
        }
        i++;
    }


    char letter;

    char sequences[25][500];
    int seqSize[25];
    memset(seqSize, 0, 25 * sizeof(int));

    char * seqName[25];
    int seqNameSize[25];
    int seqNameCapacity[25];
    memset(seqNameSize, 0, 25 * sizeof(int));
    memset(seqNameCapacity, 1, 25 * sizeof(int));

    FILE * finput = fopen(seq, "r");

    int seqCount = 0;
    int seqNameCount = 0;

    

    //input the sequences
    while (fscanf(finput ,"%c", &letter) > 0)
    {
        if(letter == '>'){
            name:
            seqName[seqNameCount] = malloc(sizeof(char));
            while (fscanf(finput, "%c", &letter) > 0){
                if (letter == '\n'){
                    seqNameCount++;
                    goto seque;
                }
                if(seqNameSize[seqNameCount] >= seqNameCapacity[seqNameCount])
                {
                    seqName[seqNameCount] = realloc(seqName[seqNameCount], seqNameCapacity[seqNameCount] * 2);
                    seqNameCapacity[seqNameCount] = seqNameCapacity[seqNameCount] * 2;
                }
                seqName[seqNameCount][seqNameSize[seqNameCount]] = letter;
                seqNameSize[seqNameCount]++;
            }
            
        }
        else{
            seque:
            while (fscanf(finput, "%c", &letter) > 0)
            {
                if (letter == '>'){
                    seqCount++;
                    goto name;
                }

                if(letter != '\n'){
                    sequences[seqCount][seqSize[seqCount]] = letter;
                    seqSize[seqCount]++;
                }
                

            }
            
        }

    }

    

    //determine the number of sequences

    int cnt = 0;
    int numOfSeqs = 0;

    while (cnt < 25)
    {
        
        if(seqSize[cnt] == 0){
            break;
        }
        cnt++;
        numOfSeqs = cnt;
    }


    seqNameCount = 0;
    while (seqNameCount < numOfSeqs)
    {
        if(seqNameSize[seqNameCount] >= seqNameCapacity[seqNameCount])
        {
            seqName[seqNameCount] = realloc(seqName[seqNameCount], seqNameCapacity[seqNameCount] * 2);
            seqNameCapacity[seqNameCount] = seqNameCapacity[seqNameCount] * 2;
        }
        seqName[seqNameCount][seqNameSize[seqNameCount]] = '\0';
        seqNameSize[seqNameCount]++;

        seqNameCount++;
    }


    
    

    char ** seqOutReverse;
    int seqOutLen;
    int score;

    double tMatrix[25][25];



    int first = 0; int second;

    char * mSeq[2];
    int mSeqSize[2];
    

    cnt = 0;
    int cnt2;
    int distance;

    while (first < numOfSeqs)
    {
        second = first + 1;
        while (second < numOfSeqs)
        {
            mSeq[0] = malloc(seqSize[first] * sizeof(char));
            strncpy(mSeq[0], sequences[first], seqSize[first]);
            mSeq[1] = malloc(seqSize[second] * sizeof(char));
            strncpy(mSeq[1], sequences[second], seqSize[second]);

            mSeqSize[0] = seqSize[first];
            mSeqSize[1] = seqSize[second];

            //Compute the match
            computeAGlobal(mSeq, mSeqSize, sMatrix, gapopen, gapext, &seqOutReverse, &seqOutLen, &score);
            
            cnt2 = seqOutLen - 1;
            distance = 0;
            while (cnt2 >= 0)
            {

                if(seqOutReverse[0][cnt2] != seqOutReverse[1][cnt2]){
                    distance++;
                }

                cnt2--;
            }

            tMatrix[first][second] = (double) distance;
            //printf("Distance between %s and %s: %d\n", seqName[first], seqName[second], distance);
            

            //printf("\n");
            free(seqOutReverse[0]);
            free(seqOutReverse[1]);
            free(mSeq[0]);
            free(mSeq[1]);
            second++;
        }

        first++;
        
    }

    
    

    struct TreeNode nodes[25];

    int nodeSize = numOfSeqs;
    //printf("nodeSize:%d\n", nodeSize);
    
    cnt = 0;
    while (cnt < nodeSize)
    {
        struct TreeNode newNode;
        newNode.height = 0;
        newNode.namelen = seqNameSize[cnt];
        newNode.name = malloc(newNode.namelen * sizeof(char));
        strncpy(newNode.name, seqName[cnt], seqNameSize[cnt]);
        newNode.numOfGenes = 1;
        nodes[cnt] = newNode;

        cnt++;
    }

    

    int a, b;
    double leftHeight;
    double rightHeight;
    double numerator;
    double denominator;
    double firstSum;
    double secondSum;
    int rowIndex;
    int columnIndex;

    int index;

    

    while (nodeSize > 1)
    {

        findIndexOfMin(tMatrix, nodeSize, &a, &b);

        struct TreeNode mergeNode;
        mergeNode.height = tMatrix[a][b]/2.0;

        leftHeight = mergeNode.height - nodes[a].height;
        rightHeight = mergeNode.height - nodes[b].height;
        mergeNode.namelen = 20 + nodes[a].namelen + nodes[b].namelen;
        mergeNode.name = malloc(mergeNode.namelen * sizeof(char));
        sprintf(mergeNode.name, "(%s:%.2lf, %s:%.2lf)%c", nodes[a].name, leftHeight, nodes[b].name, rightHeight, '\0');

        mergeNode.numOfGenes = nodes[a].numOfGenes + nodes[b].numOfGenes;

        

        cnt = 0;
        while (cnt < nodeSize)
        {
            cnt2 = cnt + 1;
            while (cnt2 < nodeSize)
            {
                if((cnt == a || cnt2 == a) && ((cnt != b && cnt2 != b))){
                    if(a == cnt){
                        firstSum = tMatrix[a][cnt2] * (double) nodes[a].numOfGenes;
                        if (b < cnt2)
                        {
                            secondSum = tMatrix[b][cnt2]* (double) nodes[b].numOfGenes;
                        }
                        else
                        {
                            secondSum = tMatrix[cnt2][b]* (double) nodes[b].numOfGenes;
                        }
                        
                    }
                    else{
                        firstSum = tMatrix[cnt][a]* (double) nodes[a].numOfGenes;
                        if (b < cnt)
                        {
                            secondSum = tMatrix[b][cnt]* (double) nodes[b].numOfGenes;
                        }
                        else
                        {
                            secondSum = tMatrix[cnt][b]* (double) nodes[b].numOfGenes;
                        }
                    }

                    numerator = firstSum + secondSum;
                    denominator = nodes[a].numOfGenes + nodes[b].numOfGenes;

                    

                    tMatrix[cnt][cnt2] = numerator/denominator;
                    
                }
                cnt2++;
                
            }
            

            cnt++;
            
        }


        rowIndex = b;
        while (rowIndex < nodeSize - 1)
        {
            columnIndex = 0;
            while (columnIndex < nodeSize)
            {
                tMatrix[rowIndex][columnIndex] = tMatrix[rowIndex + 1][columnIndex];
                columnIndex++;
            }
             

            rowIndex++;
        }

        columnIndex = b;
        while (columnIndex < nodeSize - 1)
        {
            rowIndex = 0;
            while (rowIndex < nodeSize)
            {
                tMatrix[rowIndex][columnIndex] = tMatrix[rowIndex][columnIndex+1];
                rowIndex++;
            }

            columnIndex++;
            
        }
        
        

        index = b;
        while (index < nodeSize - 1)
        {
            nodes[index].height = nodes[index+1].height;
            nodes[index].name = nodes[index+1].name;
            nodes[index].namelen = nodes[index+1].namelen;
            nodes[index].numOfGenes = nodes[index+1].numOfGenes;

            index++;
        }

        nodes[a].height = mergeNode.height;
        nodes[a].name = mergeNode.name;
        nodes[a].namelen = mergeNode.namelen;
        nodes[a].numOfGenes = mergeNode.numOfGenes;

        nodeSize--;
        
    }


    //Out
    FILE *ouf = fopen(out, "w+");
    
    fprintf(ouf, "%s;", nodes[0].name);

    fclose(ouf);


    return 0;    

}