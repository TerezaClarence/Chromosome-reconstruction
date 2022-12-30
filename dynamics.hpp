/*
    File: dynamics.hpp
    Function: Simulate dynamics of chromosomes and nuclear envelope
    Model: chromoCell
    Created: 11 December, 2017 (XF)
*/

#ifndef DYNAMICS_HPP_INCLUDED
#define DYNAMICS_HPP_INCLUDED

#include "initConfig.hpp"

// ===================== Functions =========================
void oneIter(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<iPair> vecBeadPairInteract, vector<double> & vecBeadPairInteractFreq, vector<int> & vecBeadIdsCentromere, 
             unordered_map<string, double> & energy, VoxMap & voxMap, bool FLAG_INTERACT_ON);   // if using VoxMap
// STEP 1
void updateCoord(vector<Node> & vecNode, vector<int> & vecBeadIdsCentromere);
// STEP 2
// --- following functions are for RECONSTRUCTION objective ---
void updateAccel(vector<Chromosome> & vecChromosome, vector<Node> & vecNode,
                 vector<iPair> vecBeadPairInteract, vector<double> & vecBeadPairInteractFreq,
                 VoxMap & voxMap, unordered_map<string, double> & energy, bool FLAG_INTERACT_ON);
// STEP 3
void updateVeloc(vector<Node> & vecNode, unordered_map<string, double> & energy);
void updateVelocLinearDamp(vector<Node> & vecNode, unordered_map<string, double> & energy);

double dVecDist(dVec & p1, dVec & p2);

// VoxMap algorithm
void updateVoxMap(VoxMap & voxMap, vector<Node> & vecNode);

#endif // DYNAMICS_HPP_INCLUDED
