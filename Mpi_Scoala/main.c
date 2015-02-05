//
//  main.c
//  MPIDEMO
//
//  Created by Vlad Berila on 09/10/14.
//  Copyright (c) 2014 ___VLADBERILA___. All rights reserved.
//

#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>

//Unde poate fi
typedef struct Pozitie_struct
{
    //N,E,S,V
    char directie;
    int i,j;
    int iAnterior, jAnterior;
} Pozitie;
MPI_Datatype mpi_pozitie;

typedef struct Vedere_struct
{
    //T,C,O,R
    char fata, spate, stanga, dreapta;
    //Doar pentru a ne ajuta sa nu intram in cicluri
    int i,j;
} Vedere;
MPI_Datatype mpi_vedere;

typedef struct Profesor_struct
{
    int i,j;
    char directie;
} Profesor;

Profesor profesor;
char Harta[100][100];
int hartaPatternMutari[100][100];
int N,M;

// orientarea profului e publica si globala doar ca un move helper pt el. sclavii nu stiu de ea

/*
 fata = 0
 spate = 1
 dreapta = 2
 stanga = 3
 */
int dirIProf[4], dirJProf[4];
// nr de procese
int nNumOfProcs, rank;


Pozitie pozitiiPosibile[10000];
int nPozitiiPosibile = 0;

void master();
void slave();
void computeOrientationVectors(char orientare, int dirI[4], int dirJ[4]);
void init();
int checkMatch( Vedere vedere,int i, int j,char orientare);
int checkMatchCuDeplasare(Vedere, Pozitie, char directie);
void sendToSlaveToCompute(int, Pozitie[], int, char directie, Vedere);
void receiveFromSlaves( Pozitie[], int*);
void stopSlaves();

void createMPIStruct()
{
    //Pozitie
    const int nitems=5;
    int          blocklengths[5] = {1,1,1,1,1};
    MPI_Datatype types[5] = {MPI_CHAR, MPI_INT, MPI_INT,MPI_INT, MPI_INT};
    MPI_Aint     offsets[5];
    
    offsets[0] = offsetof(Pozitie, directie);
    offsets[1] = offsetof(Pozitie, i);
    offsets[2] = offsetof(Pozitie, j);
    offsets[3] = offsetof(Pozitie, iAnterior);
    offsets[4] = offsetof(Pozitie, jAnterior);
    
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_pozitie);
    MPI_Type_commit(&mpi_pozitie);
    
    //Vedere
    const int nitems2=6;
    int          blocklengths2[6] = {1,1,1,1,1,1};
    MPI_Datatype types2[6] = {MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_INT, MPI_INT};
    MPI_Aint     offsets2[6];
    
    offsets2[0] = offsetof(Vedere, fata);
    offsets2[1] = offsetof(Vedere, spate);
    offsets2[2] = offsetof(Vedere, stanga);
    offsets2[3] = offsetof(Vedere, dreapta);
    offsets2[4] = offsetof(Vedere, i);
    offsets2[5] = offsetof(Vedere, j);
    
    MPI_Type_create_struct(nitems2, blocklengths2, offsets2, types2, &mpi_vedere);
    MPI_Type_commit(&mpi_vedere);
}

int main(int argc, char * argv[])
{
    
    // insert code here...
    MPI_Init(&argc , &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &nNumOfProcs);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    createMPIStruct();
    
    if(rank == 0)
        init();
    
    MPI_Bcast(&Harta, 100 * 100, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if(rank == 0)
        master();
    else
        slave();
    
    MPI_Type_free(&mpi_pozitie);
    MPI_Type_free(&mpi_vedere);
    MPI_Finalize();
    
    return 0;
}

void init()
{
    // citire din fisier
    // TO DO
    FILE *fp = fopen("/Users/vladberila/Documents/Dev/Mpi_Scoala/Mpi_Scoala/harta.txt", "r");
    if(fp == NULL)
        exit(7);
    fscanf(fp, "%d %d", &N, &M);
    
    char c;
    
    for(int i = 0; i <= N + 1; i++)
        Harta[i][0] = Harta[i][M + 1] = 'T';
    
    for(int i = 0; i <= M + 1; i++)
        Harta[0][i] = Harta[N + 1][i] = 'T';
    
    for(int i = 1 ; i <= N; i++)
        for(int j = 1; j <= M; j++)
        {
            fscanf(fp, "%c", &c);
            if(c != ' ' && c!= '\n')
                Harta[i][j] = c;
            else
                j--;
//            printf("%c",Harta[i][j]);
        }
    
    fscanf(fp, "%d %d", &profesor.i, &profesor.j);
    fscanf(fp,"%c", &c);
    fscanf(fp,"%c", &profesor.directie);
    
    
    fclose(fp);
    
    computeOrientationVectors(profesor.directie, dirIProf, dirJProf);
}

void computeOrientationVectors(char orientare, int dirI[4], int dirJ[4])
{
    switch (orientare)
    {
        case 'N': dirI[0] = -1, dirI[1] = 1, dirI[2] = dirI[3] = 0, dirJ[0] = 0, dirJ[1] = 0, dirJ[2] = 1, dirJ[3] = -1;
            break;
        case 'S': dirI[0] = 1, dirI[1] = -1, dirI[2] = dirI[3] = 0, dirJ[0] = 0, dirJ[1] = 0, dirJ[2] = -1, dirJ[3] = 1;
            break;
        case 'E': dirI[0] = dirI[1] = 0, dirI[2] = 1, dirI[3] = -1, dirJ[0] = 1, dirJ[1] = -1, dirJ[2] = 0, dirJ[3] = 0;
            break;
        case 'V': dirI[0] = dirI[1] = 0, dirI[2] = -1, dirI[3] = 1, dirJ[0] = -1, dirJ[1] = 1, dirJ[2] = 0, dirJ[3] = 0;
            break;
    }
}

Vedere getVedere()
{
    Vedere vedere;
    int i = profesor.i;
    int j = profesor.j;
    
    vedere.fata = Harta[i + dirIProf[0]][j + dirJProf[0]];
    vedere.spate = Harta[i + dirIProf[1]][j + dirJProf[1]];
    vedere.dreapta = Harta[i + dirIProf[2]][j + dirJProf[2]];
    vedere.stanga = Harta[i + dirIProf[3]][j + dirJProf[3]];
    vedere.i = i;
    vedere.j = j;
    
    return vedere;
}

void MoveProfessor(char directie)
{
    switch (directie)
    {
        case 'f' : profesor.i += dirIProf[0], profesor.j += dirJProf[0]; break;
        case 'b' : profesor.i += dirIProf[1], profesor.j += dirJProf[1]; break;
        case 'r' : profesor.i += dirIProf[2], profesor.j += dirJProf[2]; break;
        case 'l' : profesor.i += dirIProf[3], profesor.j += dirJProf[3]; break;
    }
    
//    printf("Profesor mutat\n");
}

//directie = f,b,l,r
//forward, back, left, right
int deplaseaza(char directie, Vedere vedereCurenta, Pozitie pozitiiPosibile[])
{
    //Nu il deplasam daca da de copac sau prapastie
    if( (directie == 'f' && (vedereCurenta.fata == 'T' || vedereCurenta.fata == 'C')) ||
       (directie == 'b' && (vedereCurenta.spate == 'T' || vedereCurenta.spate == 'C')) ||
       (directie == 'l' && (vedereCurenta.stanga == 'T' || vedereCurenta.stanga == 'C')) ||
       (directie == 'r' && (vedereCurenta.dreapta == 'T' || vedereCurenta.dreapta == 'C'))
       )
        return 0;
    
    //Il mutam pe profesor
    MoveProfessor(directie);
    
    //Updatam ce vede
    Vedere vedereNoua = getVedere();
    
    //Updatam harta daca putem merge in noua pozitie
    if(hartaPatternMutari[vedereNoua.i][vedereNoua.j] == 1)
        return 0;
    hartaPatternMutari[vedereNoua.i][vedereNoua.j] = 1;
    
    
    ///
    // Trimitem la sclavi pozitiile, directia si noua vedere
    ///
    
    int nNrCalculeSclav = nPozitiiPosibile >= nNumOfProcs ? nPozitiiPosibile / (nNumOfProcs - 1) + 1 : nNumOfProcs;
    if(nPozitiiPosibile % ( nNumOfProcs - 1 ) == 0)
        nNrCalculeSclav--;
    
    int nNrPozitiiToCompute = 0;
    for(int rank = 1; rank < nNumOfProcs; ++rank)
    {
        Pozitie pozitiiPtSclav[1000];
        nNrPozitiiToCompute = 0;
        
        for(int j = 0; j < nNrCalculeSclav && j + nNrCalculeSclav * (rank - 1) < nPozitiiPosibile; ++j)
        {
            pozitiiPtSclav[nNrPozitiiToCompute++] = pozitiiPosibile[j + nNrCalculeSclav * (rank - 1)];
        }
        
        
        //Trimitem la sclav
        printf("Sclavul %d sa proceseze %d pozitii!\n",rank,nNrPozitiiToCompute);
        
        sendToSlaveToCompute(rank, pozitiiPtSclav, nNrPozitiiToCompute, directie, vedereNoua);
    }
    
    ///
    // Ne returneaza pozitiile care se potrivesc si facem o noua lista
    ///
    Pozitie pozitiiPosibileUpdatate[1000];
    int nNrPozPosibileUpdatate = 0;
    
    receiveFromSlaves(pozitiiPosibileUpdatate, &nNrPozPosibileUpdatate);
    
    //Daca noua list are 1 element, am gasit solutia
    
    if(nNrPozPosibileUpdatate == 1)
    {
        stopSlaves();
        Pozitie p = pozitiiPosibileUpdatate[0];
        printf("%d %d %c",p.i,p.j,p.directie);
        return 1;
    }
    
    //Incercam mutare
    if(deplaseaza('f',vedereNoua,pozitiiPosibileUpdatate) == 0) if(deplaseaza('b',vedereNoua,pozitiiPosibileUpdatate) == 0) if(deplaseaza('l',vedereNoua,pozitiiPosibileUpdatate) == 0) if(deplaseaza('r',vedereNoua,pozitiiPosibileUpdatate) == 0);
    
    //Resetam Patternul, o poate lua si pe aici acum
    if(directie == 'f') MoveProfessor('b');
    if(directie == 'b') MoveProfessor('f');
    if(directie == 'l') MoveProfessor('r');
    if(directie == 'r') MoveProfessor('l');
    hartaPatternMutari[vedereNoua.i][vedereNoua.j] = 0;
    
    return 0;
    
}

int checkMatch(Vedere vedere,int i, int j,char orientare)
{
    if(Harta[i][j] != 'O' && Harta[i][j] != 'R')
        return 0;
    
    Vedere vedereCurenta = vedere;
    int dirI[4], dirJ[4];
    
    computeOrientationVectors(orientare, dirI, dirJ);
    
    if(vedereCurenta.fata == Harta[i + dirI[0]][j + dirJ[0]]
       && vedereCurenta.dreapta == Harta[i + dirI[2]][j + dirJ[2]]
       && vedereCurenta.spate == Harta[i + dirI[1]][j + dirJ[1]]
       && vedereCurenta.stanga == Harta[i + dirI[3]][j + dirJ[3]])
        return 1;
    
    
    return 0;
}

int checkMatchCuDeplasare(Vedere vedere, Pozitie pozitie, char directie)
{
    int dirI[4], dirJ[4];
    computeOrientationVectors(pozitie.directie, dirI, dirJ);
    switch (directie)
    {
        case 'f' : pozitie.i += dirI[0], pozitie.j += dirJ[0]; break;
        case 'b' : pozitie.i += dirI[1], pozitie.j += dirJ[1]; break;
        case 'r' : pozitie.i += dirI[2], pozitie.j += dirJ[2]; break;
        case 'l' : pozitie.i += dirI[3], pozitie.j += dirJ[3]; break;
    }
    return checkMatch(vedere, pozitie.i, pozitie.j, pozitie.directie);
}

void master()
{
    Vedere vedereCurenta = getVedere();
    
//    printf("%d %d %c",profesor.i,profesor.j,profesor.directie);
    
    ///
    // Gasim pozitiile initiale posibile
    ///
    for(int i = 1; i <= N; i++)
        for(int j = 1; j <= M; j++)
        {
            //Directia N
            if(checkMatch(vedereCurenta,i,j,'N'))
            {
                Pozitie pozitie;
                pozitie.directie = 'N';
                pozitie.i = i;
                pozitie.j = j;
                pozitie.iAnterior = -1;
                pozitie.jAnterior = -1;
                pozitiiPosibile[nPozitiiPosibile++] = pozitie;
            }
            
            //Directia S
            if(checkMatch(vedereCurenta,i,j,'S'))
            {
                Pozitie pozitie;
                pozitie.directie = 'S';
                pozitie.i = i;
                pozitie.j = j;
                pozitie.iAnterior = -1;
                pozitie.jAnterior = -1;
                pozitiiPosibile[nPozitiiPosibile++] = pozitie;
            }
            
            //Directia E
            if(checkMatch(vedereCurenta,i,j,'E'))
            {
                Pozitie pozitie;
                pozitie.directie = 'E';
                pozitie.i = i;
                pozitie.j = j;
                pozitie.iAnterior = -1;
                pozitie.jAnterior = -1;
                pozitiiPosibile[nPozitiiPosibile++] = pozitie;
            }
            
            //Directia V
            if(checkMatch(vedereCurenta,i,j,'V'))
            {
                Pozitie pozitie;
                pozitie.directie = 'V';
                pozitie.i = i;
                pozitie.j = j;
                pozitie.iAnterior = -1;
                pozitie.jAnterior = -1;
                pozitiiPosibile[nPozitiiPosibile++] = pozitie;
            }
            
        }
    
    //printf("\n%d\n",nPozitiiPosibile);
    
    for(int i=0;i<nPozitiiPosibile;i++)
    {
        Pozitie p = pozitiiPosibile[i];
        printf("%d %d %c\n",p.i,p.j,p.directie);
    }
    
    ///
    // Il plimbam pe profesor prin padure
    ///
    hartaPatternMutari[vedereCurenta.i][vedereCurenta.j] = 1;
    if(deplaseaza('f',vedereCurenta,pozitiiPosibile) == 0)
        if(deplaseaza('b',vedereCurenta,pozitiiPosibile) == 0)
            if(deplaseaza('l',vedereCurenta,pozitiiPosibile) == 0)
                if(deplaseaza('r',vedereCurenta,pozitiiPosibile) == 0);
    
    
    
    
    
    
}

void stopSlaves()
{
    for(int rank = 1; rank < nNumOfProcs; ++rank)
    {
        int n = -321;
        MPI_Send(&n, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
    }
}

void sendToSlaveToCompute( int rank, Pozitie pozitii[], int nNrPozitiiToCompute, char directie, Vedere vedere)
{
    printf("Rank: %d\n",rank);
    for(int i=0;i<nNrPozitiiToCompute;i++)
    {
        Pozitie p = pozitii[i];
        printf("M: %d %d %c\n",p.i,p.j,p.directie);
    }
    
    printf("Vdere trimisa: %c %c %c %c %d %d\n",vedere.fata,vedere.spate,vedere.stanga,vedere.dreapta,vedere.i,vedere.j);
    
    MPI_Send(&nNrPozitiiToCompute, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
  //    MPI_Send(&pozitii, nNrPozitiiToCompute, mpi_pozitie, rank, 1, MPI_COMM_WORLD);
    int I[1000],J[1000];
    char C[1000];
    for(int i=0; i < nNrPozitiiToCompute; i++)
    {
        I[i] = pozitii[i].i;
        J[i] = pozitii[i].j;
        C[i] = pozitii[i].directie;
    }
    MPI_Send(&I, nNrPozitiiToCompute, MPI_INT, rank, 10, MPI_COMM_WORLD);
    MPI_Send(&J, nNrPozitiiToCompute, MPI_INT, rank, 11, MPI_COMM_WORLD);
    MPI_Send(&C, nNrPozitiiToCompute, MPI_CHAR, rank, 12, MPI_COMM_WORLD);

//    MPI_Send(&vedere, 1, mpi_vedere, rank, 2, MPI_COMM_WORLD);
    MPI_Send(&(vedere.fata), 1, MPI_CHAR, rank, 20, MPI_COMM_WORLD);
    MPI_Send(&(vedere.spate), 1, MPI_CHAR, rank, 21, MPI_COMM_WORLD);
    MPI_Send(&(vedere.stanga), 1, MPI_CHAR, rank, 22, MPI_COMM_WORLD);
    MPI_Send(&(vedere.dreapta), 1, MPI_CHAR, rank, 23, MPI_COMM_WORLD);
    MPI_Send(&(vedere.i), 1, MPI_INT, rank, 24, MPI_COMM_WORLD);
    MPI_Send(&(vedere.j), 1, MPI_INT, rank, 25, MPI_COMM_WORLD);
    
    MPI_Send(&directie, 1, MPI_CHAR, rank, 3, MPI_COMM_WORLD);
}

void receiveFromSlaves(Pozitie pozitii[], int* nPozitii)
{
    (*nPozitii) = 0;
    Pozitie pozitiiDePrimit[1000];
    int nPozitiiDePrimit = 0;
    
    MPI_Status st;
    
    for(int rank = 1; rank < nNumOfProcs; ++rank)
    {
        nPozitiiDePrimit = 0;
        MPI_Recv(&nPozitiiDePrimit, 1, MPI_INT, rank, 0, MPI_COMM_WORLD, &st);
        if(st.MPI_ERROR != MPI_SUCCESS)
            printf("Primu in receiveFromSlaves");
        MPI_Recv(&pozitiiDePrimit, nPozitiiDePrimit, mpi_pozitie, 0, 1, MPI_COMM_WORLD, &st);
        if(st.MPI_ERROR != MPI_SUCCESS)
            printf("Al doilea in receiveFromSlaves");
        
        for(int j = 0; j < nPozitiiDePrimit; ++j)
            pozitii[(*nPozitii)++] = pozitiiDePrimit[j];
    }
}

void slave()
{
    Pozitie pozitiiDeProcesat[1000];
    Vedere vedere;
    Pozitie pozitiiDeReturnat[1000];
    int nPozitii, nPozitiiDeReturnat = 0;
    char directie;
    
    MPI_Status st;
    
//    printf("Sclav\n");
    while(1)
    {
//        printf("Intru in while\n");
        nPozitiiDeReturnat = 0;
        MPI_Recv(&nPozitii, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &st);
        
        printf("S %d: Nr pozitii de procesat: %d\n", rank, nPozitii);
        
        if(nPozitii == -321)
            break;
        
//        if(st.MPI_ERROR != MPI_SUCCESS)
//            printf("Primu in slaves");
        
        //MPI_Recv(&pozitiiDeProcesat, nPozitii, mpi_pozitie, 0, 1, MPI_COMM_WORLD, &st);
        int I[1000],J[1000];
        char C[1000];
        MPI_Recv(&I, nPozitii, MPI_INT, 0, 10, MPI_COMM_WORLD, &st);
        MPI_Recv(&J, nPozitii, MPI_INT, 0, 11, MPI_COMM_WORLD, &st);
        MPI_Recv(&C, nPozitii, MPI_CHAR, 0, 12, MPI_COMM_WORLD, &st);
        for(int i=0; i < nPozitii; i++)
        {
            pozitiiDeProcesat[i].i = I[i];
            pozitiiDeProcesat[i].j = J[i];
            pozitiiDeProcesat[i].directie = C[i];
        }

//        if(st.MPI_ERROR != MPI_SUCCESS)
//            printf("Al doilea in slaves");

        //MPI_Recv(&vedere, 1, mpi_vedere, 0, 2, MPI_COMM_WORLD, &st);
        MPI_Recv(&(vedere.fata), 1, MPI_CHAR, 0, 20, MPI_COMM_WORLD, &st);
        MPI_Recv(&(vedere.spate), 1, MPI_CHAR, 0, 21, MPI_COMM_WORLD, &st);
        MPI_Recv(&(vedere.stanga), 1, MPI_CHAR, 0, 22, MPI_COMM_WORLD, &st);
        MPI_Recv(&(vedere.dreapta), 1, MPI_CHAR, 0, 23, MPI_COMM_WORLD, &st);
        MPI_Recv(&(vedere.i), 1, MPI_INT, 0, 24, MPI_COMM_WORLD, &st);
        MPI_Recv(&(vedere.j), 1, MPI_INT, 0, 25, MPI_COMM_WORLD, &st);

        MPI_Recv(&directie, 1, MPI_CHAR, 0, 3, MPI_COMM_WORLD, &st);
        
//        if(st.MPI_ERROR != MPI_SUCCESS)
//            printf("Al treilea in slaves");
        
        printf("S %d: Vdere primita: %c %c %c %c %d %d\n",rank,vedere.fata,vedere.spate,vedere.stanga,vedere.dreapta,vedere.i,vedere.j);
    
        for(int i=0; i < nPozitii; i++)
        {
            Pozitie p = pozitiiDeProcesat[i];
            printf("S %d: %d %d %c\n",rank,p.i,p.j,p.directie);
        
            if(checkMatchCuDeplasare(vedere, p, directie))
                pozitiiDeReturnat[nPozitiiDeReturnat++] = p;
        }
        
        printf("S %d: Pozitii inca bune: %d\n",rank,nPozitiiDeReturnat);
    
        MPI_Send(&nPozitiiDeReturnat, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&pozitiiDeReturnat, nPozitiiDeReturnat, mpi_pozitie, 0, 1, MPI_COMM_WORLD);
    }
}