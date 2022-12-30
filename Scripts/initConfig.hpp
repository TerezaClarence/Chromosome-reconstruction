/*
    File: initConfig.hpp
    Function: Initialize chromosomes
    Model: chromoCell
    Created: 7 December, 2017 (XF)
*/

#ifndef INITCONFIG_HPP_INCLUDED
#define INITCONFIG_HPP_INCLUDED

#include <stdio.h>
//#include <math.h>
#include <iostream>
#include <fstream>

#include <cstring>
#include <string>
#include <sstream>
#include <iomanip>

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cmath>

#include<random>
#include<chrono>

#include<omp.h>

using namespace std;

typedef struct {double x,y,z;} dVec;
typedef struct {int i,j,k;} iVec;
typedef struct {int i,j;} iPair;

// system settings
const double DT                 = 1E-3;     // s; 1E-3 used for RECONSTRUCTION
const double T                  = 1E5;      // s; 1E5 used for RECONSTRUCTION
const int N_ITER_TOTAL          = T/DT;
const int CNT_START_INTERACT    = 0.0*N_ITER_TOTAL;
const int FREQ_WRITE            = N_ITER_TOTAL/1000;
const int FREQ_PRINT            = N_ITER_TOTAL/10000;
const int NUM_THREADS           = 8;   // ????why 8?
const bool OMP_ON               = false; // ???what it means?
// constants
const double PI                 = 3.1415926535897;
const double AVOG               = 6.02214E+23;        // 1/mol, Avogadro number
const double NM_PER_BP          = 0.34;               // nm/bp; old value used in RECONSTRUCTION: 4.6875E-4 (where does it come from?)
const double RAD_NUCLEUS        = 300;               // nm; 6000 nm diameter average for mammalian nucleus;710 used for budding yeast; 1000(?) for S. pombe ; generally 1700-2000(??) for animal cell nuclues
/* ------------ Parameters Reconstruct ------------- */
const string SIM_TYPE           = "RECONSTRUCT";   // ???explain the syntax?

// structure
const double GENOMIC_RESOLUTION = 150000;     // bp (base pair ) "RECONSTRUCT" 3000 for budding data; 2000 for fission data;
//const double RAD_BEAD           = GENOMIC_RESOLUTION*NM_PER_BP;
const double RAD_BEAD           = 0.01;       // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!need some estimation!
const double MOLW_PER_BP_DNA    = 650;      // g/mol or dalton
const double MOLW_PER_NUCLEO    = 2.1E+5;   // g/mol or dalton
const int    N_NUCLEO_PER_BP    = 146;      // 146 base pairs wrap around nucleosome    ???constant with double precision = const double, constant no precision = const int??
const double MASS_BEAD          = (MOLW_PER_BP_DNA+MOLW_PER_NUCLEO/N_NUCLEO_PER_BP)/AVOG*GENOMIC_RESOLUTION *1E-3;    // kg, ~1E-23 (this is pure DNA. how to account for nucleosome?)

// dynamics
const double CUTOFF_NORMMAT     = 75.0;      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!something between 0.001 and 0.005 
const int    GENOMIC_LEAST_SEP_INTRA  = 0;  // make it 0 when nearest interactions are accounted for according to HiC map
const double K_NEAR             = 5E-23;    // nN/nm
const double K_INTRA            = 1.0*K_NEAR;
const double K_INTER            = 0.04*K_NEAR; 
const double K_REP		= 0.01*K_NEAR;

// data input
const string SUFFIX_normMatrix = ".normMatrix";
const string SUFFIX_numberBead = ".numberBead";
/* -------------------------------------------------- */
// ===================== Classes =========================

/* -------------- Node ---------------- */
class Node
{
    int id, hostChromoId;
    dVec coord; // nm
    dVec veloc; // nm/s
    dVec accel, accelPrev; // nm/s/s   ????? difference between dVec and accelPrev????

public:
    Node ();
    Node (int id, int hostChromoId);
    Node (int id, int hostChromoId, dVec coord, dVec veloc, dVec accel, dVec accelPrev);
    ~Node ();

    // setters  ???? what setters do?
    void set_id(int id) {this->id=id;}
    void set_hostChromoId(int hostChromoId) {this->hostChromoId=hostChromoId;}
    void set_coord(dVec coord) {this->coord=coord;}
    void set_veloc(dVec veloc) {this->veloc=veloc;}
    void set_accel(dVec accel) {this->accel=accel;}
    void set_accelPrev(dVec accelPrev) {this->accelPrev=accelPrev;}

    // getters  ???what getters do?
    int get_id() const {return this->id;}
    int get_hostChromoId() const {return this->hostChromoId;}
    dVec get_coord() const {return this->coord;}
    dVec get_veloc() const {return this->veloc;}
    dVec get_accel() const {return this->accel;}
    dVec get_accelPrev() const {return this->accelPrev;}
};

/* -------------- Chromosomes ---------------- */
class Chromosome
{
    int id;
    vector<int> vecBeadIds; // may change to vector of pointer to Node objects in vecNode
    vector<int> vecGeneIds; // may change to vector of iPairs {start, end}
    //iPair lociTelomere;
    //iPair lociCentromere;
public:
    Chromosome ();
    Chromosome (int id, vector<int> vecBeadIds);
    Chromosome (int id, vector<int> vecBeadIds, vector<int> vecGeneIds);
    ~Chromosome();

    // setters
    void set_id(int id) {this->id = id;}
    void set_vecBeadIds(vector<int> vecBeadIds) {this->vecBeadIds = vecBeadIds;}
    void set_vecGeneIds(vector<int> vecGeneIds) {this->vecGeneIds = vecGeneIds;}

    // getters
    int get_id() const {return this->id;}
    vector<int> get_vecBeadIds() const {return this->vecBeadIds;}
    vector<int> get_vecGeneIds() const {return this->vecGeneIds;}
};

/* -------------- VoxMap ----------------- */           // ??? what does VoxMap do?
class VoxMap
{
    iVec dim;
    unordered_map<int, vector<int>> mapBeadIds;         // contains Bead ids in VoxBox; key is collapsed from 3d indexing
    unordered_map<int, vector<int>> mapBeadIdsHalo;     // contains Bead ids in halo of VoxBox;

public:
    VoxMap ();
    VoxMap (iVec dim, unordered_map<int, vector<int>> mapBeadIds, unordered_map<int, vector<int>> mapBeadIdsHalo);
    ~VoxMap ();

    // setters
    void set_dim(iVec dim) {this->dim = dim;}
    void set_mapBeadIds(unordered_map<int, vector<int>> mapBeadIds) {this->mapBeadIds = mapBeadIds;}
    void set_mapBeadIdsHalo(unordered_map<int, vector<int>> mapBeadIdsHalo) {this->mapBeadIdsHalo = mapBeadIdsHalo;}

    // getters
    iVec get_dim() const {return this->dim;}
    unordered_map<int, vector<int>> get_mapBeadIds() const {return this->mapBeadIds;}
    unordered_map<int, vector<int>> get_mapBeadIdsHalo() const {return this->mapBeadIdsHalo;}
};
// ===================== Functions =========================                   
// --- following functions are for RECONSTRUCTION objective ---
void initChromosome(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<iPair> & vecBeadPairInteract, \
                    vector<double> & vecBeadPairInteractFreq, int argc, char ** argv);
void initDynamics(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<int> & vecBeadIdsCentromere, int argc, char ** argv, string initType, vector<double> & randNum01);
// ------------------------------------------------------------
void printChromoInfo(vector<Chromosome> & vecChromosome, vector<Node> & vecNode);
void writeChromoInfo(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<int> & vecCondensinSite, vector<int> & vecBeadIdsCentromere, double t_now);

#endif // INITCONFIG_HPP_INCLUDED

