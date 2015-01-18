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

//Unde poate fi
struct Pozitie
{
    //N,E,S,V
    char directie;
    int i,j;
    int iAnterior, jAnterior;
};

struct Vedere
{
    //T,C,O,R
    char fata, spate, stanga, dreapta;
    //Doar pentru a ne ajuta sa nu intram in cicluri
    int i,j;
};

struct Profesor
{
    int i,j;
    char directie;
} profesor;

int Harta[100][100];
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


struct Pozitie pozitiiPosibile[10000];
int nPozitiiPosibile = 0;

void master();
void slave();
void computeOrientationVectors(char orientare, int dirI[4], int dirJ[4]);
void init();
int checkMatch(struct Vedere vedere,int i, int j,char orientare);

int main(int argc, char * argv[])
{
    
    // nr de procese
    int nNumOfProcs, rank, numarPePrc;
    
    // insert code here...
    MPI_Init(&argc , &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &nNumOfProcs);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if(rank == 0)
        master();
    else
        slave();
    
    MPI_Finalize();
    
    return 0;
}

void init()
{
    // citire din fisier
    // TO DO
    
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

struct Vedere getVedere()
{
    struct Vedere vedere;
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
int deplaseaza(char directie, struct Vedere vedereCurenta, struct Pozitie pozitiiPosibile[])
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
    struct Vedere vedereNoua = getVedere();
    
    ///
    // Trimitem la sclavi pozitiile, directia si noua vedere
    ///
    //TO DO
    
    ///
    // Ne returneaza pozitiile care se potrivesc si facem o noua lista
    ///
    struct Pozitie pozitiiPosibileUpdatate[1000];
    
    //Daca noua list are 1 element, am gasit solutia
    /*
        return 1;
     */
    
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

int checkMatch(struct Vedere vedere,int i, int j,char orientare)
{
    struct Vedere vedereCurenta = vedere;
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
    struct Vedere vedereCurenta = getVedere();
    
    ///
    // Gasim pozitiile initiale posibile
    ///
    for(int i = 1; i <= N; i++)
        for(int j = 1; j <= M; j++)
        {
            char dD = 'X';
            //Directia N
            if(checkMatch(vedereCurenta,i,j,'N'))
                dD = 'N';
            
            //Directia S
            if(checkMatch(vedereCurenta,i,j,'S'))
                dD = 'S';
            
            //Directia E
            if(checkMatch(vedereCurenta,i,j,'E'))
                dD = 'E';

            //Directia V
            if(checkMatch(vedereCurenta,i,j,'V'))
                dD = 'V';

            if(dD != 'X')
            {
                struct Pozitie pozitie;
                pozitie.directie = dD;
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