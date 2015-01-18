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
};

struct Profesor
{
    int i,j;
    char directie;
} profesor;

int Harta[100][100];
int N,M;

struct Pozitie pozitiiPosibile[10000];
int nPozitiiPosibile = 0;

void master();
void slave();

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

struct Vedere getVedere()
{
    struct Vedere vedere;
    int i = profesor.i;
    int j = profesor.j;
    
    if(profesor.directie == 'N')
    {
        vedere.fata = Harta[i-1][j];
        vedere.dreapta = Harta[i][j+1];
        vedere.spate = Harta[i+1][j];
        vedere.stanga = Harta[i][j-1];
    }
    
    if(profesor.directie == 'S')
    {
        vedere.fata = Harta[i+1][j];
        vedere.dreapta = Harta[i][j-1];
        vedere.spate = Harta[i-1][j];
        vedere.stanga = Harta[i][j+1];
    }
    
    if(profesor.directie == 'E')
    {
        vedere.fata = Harta[i][j+1];
        vedere.dreapta = Harta[i+1][j];
        vedere.spate = Harta[i][j-1];
        vedere.stanga = Harta[i-1][j];
    }
    
    if(profesor.directie == 'V')
    {
        vedere.fata = Harta[i][j-1];
        vedere.dreapta = Harta[i-1][j];
        vedere.spate = Harta[i][j+1];
        vedere.stanga = Harta[i+1][j];
    }
    
    return vedere;
}

//directie = f,b,l,r
//forward, back, left, right
int deplaseaza(char directie, struct Pozitie pozitiiPosibile[])
{
    struct Vedere vedereCurenta = getVedere();
    
    //Nu il deplasam daca da de copac sau prapastie
    if( (directie == 'f' && (vedereCurenta.fata == 'T' || vedereCurenta.fata == 'C')) ||
       (directie == 'b' && (vedereCurenta.spate == 'T' || vedereCurenta.spate == 'C')) ||
       (directie == 'l' && (vedereCurenta.stanga == 'T' || vedereCurenta.stanga == 'C')) ||
       (directie == 'r' && (vedereCurenta.dreapta == 'T' || vedereCurenta.dreapta == 'C'))
       )
        return 0;
    
    //Il mutam pe profesor
    //TO DO
    
    //Updatam ce vede
    struct Vedere vedereNoua = getVedere();
    
    ///
    // Trimitem la sclavi pozitiile, directia si noua vedere
    ///
    //TO DO
    
    // Ne returneaza pozitiile care se potrivesc si refacem lista
    //
    // Obs! Nu trebuie sa retinem pozitiile eliminate in eventualitatea intoarcerii in recursivitate
    // deoarece pozitiile eliminate nu pot fi solutii
    ///
    
    //Incercam mutare
    if(deplaseaza('f') == 0) if(deplaseaza('b') == 0) if(deplaseaza('l') == 0) if(deplaseaza('r') == 0);
    
    //Daca ne blocam ne intoarcem
    //TO DO
    
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
            //Directia N
            if(vedereCurenta.fata == Harta[i-1][j] && vedereCurenta.dreapta == Harta[i][j+1] && vedereCurenta.spate == Harta[i+1][j] && vedereCurenta.stanga == Harta[i][j-1])
            {
                struct Pozitie pozitie;
                pozitie.directie = 'N';
                pozitie.i = i;
                pozitie.j = j;
                pozitie.iAnterior = -1;
                pozitie.jAnterior = -1;
                pozitiiPosibile[nPozitiiPosibile++] = pozitie;
            }
            
            //Directia S
            if(vedereCurenta.fata == Harta[i+1][j] && vedereCurenta.dreapta == Harta[i][j-1] && vedereCurenta.spate == Harta[i-1][j] && vedereCurenta.stanga == Harta[i][j+1])
            {
                struct Pozitie pozitie;
                pozitie.directie = 'S';
                pozitie.i = i;
                pozitie.j = j;
                pozitie.iAnterior = -1;
                pozitie.jAnterior = -1;
                pozitiiPosibile[nPozitiiPosibile++] = pozitie;
            }
            
            //Directia E
            if(vedereCurenta.fata == Harta[i][j+1] && vedereCurenta.dreapta == Harta[i+1][j] && vedereCurenta.spate == Harta[i][j-1] && vedereCurenta.stanga == Harta[i-1][j])
            {
                struct Pozitie pozitie;
                pozitie.directie = 'E';
                pozitie.i = i;
                pozitie.j = j;
                pozitie.iAnterior = -1;
                pozitie.jAnterior = -1;
                pozitiiPosibile[nPozitiiPosibile++] = pozitie;
            }
            
            //Directia V
            if(vedereCurenta.fata == Harta[i][j-1] && vedereCurenta.dreapta == Harta[i-1][j] && vedereCurenta.spate == Harta[i][j+1] && vedereCurenta.stanga == Harta[i+1][j])
            {
                struct Pozitie pozitie;
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
    if(deplaseaza('f',pozitiiPosibile) == 0)
        if(deplaseaza('b',pozitiiPosibile) == 0)
            if(deplaseaza('l',pozitiiPosibile) == 0)
                if(deplaseaza('r',pozitiiPosibile) == 0);
    
    
    
    
    
    
}