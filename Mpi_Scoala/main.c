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
void sendToSlaveToCompute(int, Pozitie[], Vedere, char directie);
void receiveFromSlaves( Pozitie[], int*);

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
    
    /*
    if(rank == 0)
        master();
    else
        slave();
     */

    MPI_Type_free(&mpi_pozitie);
    MPI_Type_free(&mpi_vedere);
    MPI_Finalize();
    
    return 0;
}

void init()
{
    // citire din fisier
    // TO DO
    FILE *fp = fopen("harta.txt", "r");
    if(fp == NULL)
        exit(7);
    fscanf(fp, "%d %d", &N, &M);
    
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
}

//directie = f,b,l,r
//forward, back, left, right
int deplaseaza(char directie, Vedere vedereCurenta, Pozitie pozitiiPosibile[])
{
    if(hartaPatternMutari[vedereCurenta.i][vedereCurenta.j] == 1)
        return 0;
    //Nu il deplasam daca da de copac sau prapastie
    if( (directie == 'f' && (vedereCurenta.fata == 'T' || vedereCurenta.fata == 'C')) ||
       (directie == 'b' && (vedereCurenta.spate == 'T' || vedereCurenta.spate == 'C')) ||
       (directie == 'l' && (vedereCurenta.stanga == 'T' || vedereCurenta.stanga == 'C')) ||
       (directie == 'r' && (vedereCurenta.dreapta == 'T' || vedereCurenta.dreapta == 'C'))
       )
        return 0;
    
    //Il mutam pe profesor
    MoveProfessor(directie);
    
    hartaPatternMutari[vedereCurenta.i][vedereCurenta.j] = 1;
    
    //Updatam ce vede
    Vedere vedereNoua = getVedere();
    
    ///
    // Trimitem la sclavi pozitiile, directia si noua vedere
    ///

    int nNrCalculeSclav = nPozitiiPosibile >= nNumOfProcs ? nPozitiiPosibile / (nNumOfProcs - 1) + 1 : nNumOfProcs;
    if(nPozitiiPosibile % ( nNumOfProcs - 1 ) == 0)
        nNrCalculeSclav--;
    
    for(int rank = 1; rank < nNumOfProcs; ++rank)
    {
        Pozitie pozitiiPtSclav[1000];
        
        for(int j = 0; j < nNrCalculeSclav && j + nNrCalculeSclav * (rank - 1) < nPozitiiPosibile; ++j)
        {
            pozitiiPtSclav[j] = pozitiiPosibile[j + nNrCalculeSclav * (rank - 1)];
        }
        
        sendToSlaveToCompute(rank, pozitiiPtSclav, vedereNoua, directie);
    }
    
    ///
    // Ne returneaza pozitiile care se potrivesc si facem o noua lista
    ///
    Pozitie pozitiiPosibileUpdatate[1000];
    int nNrPozPosibileUpdatate = 0;
    
    receiveFromSlaves(pozitiiPosibileUpdatate, &nNrPozPosibileUpdatate);
    
    //Daca noua list are 1 element, am gasit solutia

    if(nNrPozPosibileUpdatate == 1)
        return 1;
    
    //Incercam mutare
    if(deplaseaza('f',vedereNoua,pozitiiPosibileUpdatate) == 0) if(deplaseaza('b',vedereNoua,pozitiiPosibileUpdatate) == 0) if(deplaseaza('l',vedereNoua,pozitiiPosibileUpdatate) == 0) if(deplaseaza('r',vedereNoua,pozitiiPosibileUpdatate) == 0);
    
    //Resetam Patternul, o poate lua si pe aici acum
    if(directie == 'f') MoveProfessor('b');
    if(directie == 'b') MoveProfessor('f');
    if(directie == 'l') MoveProfessor('r');
    if(directie == 'r') MoveProfessor('l');
    hartaPatternMutari[vedereCurenta.i][vedereCurenta.j] = 0;
    
    return 0;
    
}

int checkMatch(Vedere vedere,int i, int j,char orientare)
{
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

void master()
{
    Vedere vedereCurenta = getVedere();
    
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
    
    ///
    // Il plimbam pe profesor prin padure
    ///
    hartaPatternMutari[vedereCurenta.i][vedereCurenta.j] = 1;
    if(deplaseaza('f',vedereCurenta,pozitiiPosibile) == 0)
        if(deplaseaza('b',vedereCurenta,pozitiiPosibile) == 0)
            if(deplaseaza('l',vedereCurenta,pozitiiPosibile) == 0)
                if(deplaseaza('r',vedereCurenta,pozitiiPosibile) == 0);
    
    
    
    
    
    
}

void sendToSlaveToCompute( int rank, Pozitie pozitii[], Vedere vedere, char directie)
{
    
}

void receiveFromSlaves(Pozitie pozitii[], int* nPozitii)
{
    

void slave()
{
    
}}