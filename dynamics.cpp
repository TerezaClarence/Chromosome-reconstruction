/*
    File: dynamics.cpp
    Function: Simulate dynamics of chromosomes and nuclear envelope
    Model: chromoCell
    Created: 11 December, 2017 (XF)
*/

#include "dynamics.hpp"

// ===================== Functions =========================
void oneIter(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<iPair> vecBeadPairInteract, vector<double> & vecBeadPairInteractFreq, vector<int> & vecBeadIdsCentromere,
             unordered_map<string, double> & energy, VoxMap & voxMap, bool FLAG_INTERACT_ON)
{
    /* --- Apply VoxMap algorithm to update the voxMap for detecting collision
    STEP 1 : assign each bead into VoxMap by doing division and modulo
    --- */
    updateVoxMap(voxMap, vecNode);
    //exit(123);

    /* --- Apply verlet integration (https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet)
    STEP 1 : x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt*dt
    STEP 2 : a(t+dt) updated according to x(t+dt)
    STEP 3 : v(t+dt) = v(t) + 0.5* ( a(t)+a(t+dt) ) * dt
    --- */

    /* -- STEP 1 -- */
    updateCoord(vecNode, vecBeadIdsCentromere);
    //updateCoord_OMP(vecNode);

    /* -- STEP 2 -- */
    updateAccel(vecChromosome, vecNode, vecBeadPairInteract, vecBeadPairInteractFreq, voxMap, energy, FLAG_INTERACT_ON);
    
    /* -- STEP 3 -- */
    //updateVeloc(vecNode, energy);
    updateVelocLinearDamp(vecNode, energy);
}
// voxmap
void updateVoxMap(VoxMap & voxMap, vector<Node> & vecNode)
{
    /* indexing i = nz*(dim.i*dim.j) + ny*(dim.k) + nx */
    iVec dim = voxMap.get_dim();
    int nMax = dim.i*dim.j*dim.k;
    double rMax = RAD_NUCLEUS;
    double aVoxX = rMax*2./dim.i;
    double aVoxY = rMax*2./dim.j;
    double aVoxZ = rMax*2./dim.k;

    unordered_map<int, vector<int>> mapBeadIds;
    unordered_map<int, vector<int>> mapBeadIdsHalo;
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it ++)
    {
        double xx, yy, zz;
        int nx, ny, nz;
        dVec coord = it->get_coord();
        xx = (coord.x+rMax) / aVoxX;
        yy = (coord.y+rMax) / aVoxY;
        zz = (coord.z+rMax) / aVoxZ;

        if (xx < 0 || yy < 0 || zz < 0) // warn negative values
        {
            cout << "Please offset coordinates so that x, y, z values are positive !" << endl;
            exit(234);
        }

        nx = (int) floor( xx );
        ny = (int) floor( yy );
        nz = (int) floor( zz );
        int n = nz*(dim.i*dim.j) + ny*dim.i + nx;
        int bId = it->get_id();
        mapBeadIds[n].push_back(bId);

        // OFFSET in space to examine if this Node is within Halo space of neighboring 26 boxes
        int nx2, ny2, nz2;
        double aOffsetX = 2.*RAD_BEAD/aVoxX, aOffsetY = 2.*RAD_BEAD/aVoxY, aOffsetZ = 2.*RAD_BEAD/aVoxZ;

        // TYPE 1 : 6 faces
        if (true)
        {
            nx2 = (int) floor( xx+aOffsetX );         // plus x
            if (nx2 == nx+1 && n+1 < nMax)
                mapBeadIdsHalo[n+1].push_back(bId);
            nx2 = (int) floor( xx-aOffsetX );         // minus x
            if (nx2 == nx-1 && n-1 >= 0)
                mapBeadIdsHalo[n-1].push_back(bId);
            // ---------------------------------------------------
            ny2 = (int) floor( yy+aOffsetY );         // plus y
            if (ny2 == ny+dim.i && n+dim.i < nMax)
                mapBeadIdsHalo[n+dim.i].push_back(bId);
            ny2 = (int) floor( yy-aOffsetY );         // minus y
            if (ny2 == ny-dim.i && n-dim.i >= 0)
                mapBeadIdsHalo[n-dim.i].push_back(bId);
            // ---------------------------------------------------
            nz2 = (int) floor( zz+aOffsetZ );         // plus y
            if (nz2 == nz+dim.i*dim.j && n+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // minus y
            if (nz2 == nz-dim.i*dim.j && n-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n-dim.i*dim.j].push_back(bId);
        }
        // TYPE 2 : 12 edges
        if (true)
        {
            nx2 = (int) floor( xx+aOffsetX );         // plus x, plus y
            ny2 = (int) floor( yy+aOffsetY );
            if (nx2 == nx+1 && ny2 == ny+1 && n+1+dim.i < nMax)
                mapBeadIdsHalo[n+1+dim.i].push_back(bId);
            ny2 = (int) floor( yy-aOffsetY );         // plus x, minus y
            if (nx2 == nx+1 && ny2 == ny-1 && n+1-dim.i >= 0)
                mapBeadIdsHalo[n+1-dim.i].push_back(bId);
            nx2 = (int) floor( xx-aOffsetX );         // minus x, plus y
            ny2 = (int) floor( yy+aOffsetY );
            if (nx2 == nx-1 && ny2 == ny+1 && n-1+dim.i < nMax)
                mapBeadIdsHalo[n-1+dim.i].push_back(bId);
            ny2 = (int) floor( yy-aOffsetY );         // minus x, minus y
            if (nx2 == nx-1 && ny2 == ny-1 && n-1-dim.i >= 0)
                mapBeadIdsHalo[n-1-dim.i].push_back(bId);
            // ----------------------------------------------------------
            nx2 = (int) floor( xx+aOffsetX );         // plus x, plus z
            nz2 = (int) floor( zz+aOffsetZ );
            if (nx2 == nx+1 && nz2 == nz+1 && n+1+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n+1+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // plus x, minus z
            if (nx2 == nx+1 && nz2 == nz-1 && n+1-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n+1-dim.i*dim.j].push_back(bId);
            nx2 = (int) floor( xx-aOffsetX );         // minus x, plus z
            nz2 = (int) floor( zz+aOffsetZ );
            if (nx2 == nx-1 && nz2 == nz+1 && n-1+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n-1+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // minus x, minus z
            if (nx2 == nx-1 && nz2 == nz-1 && n-1-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n-1-dim.i*dim.j].push_back(bId);
            // ----------------------------------------------------------
            ny2 = (int) floor( yy+aOffsetY );         // plus y, plus z
            nz2 = (int) floor( zz+aOffsetZ );
            if (ny2 == ny+1 && nz2 == nz+1 && n+dim.i+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n+dim.i+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // plus y, minus z
            if (ny2 == ny+1 && nz2 == nz-1 && n+dim.i-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n+dim.i-dim.i*dim.j].push_back(bId);
            ny2 = (int) floor( yy-aOffsetY );         // minus y, plus z
            nz2 = (int) floor( zz+aOffsetZ );
            if (ny2 == ny-1 && nz2 == nz+1 && n-dim.i+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n-dim.i+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // minus y, minus z
            if (ny2 == ny-1 && nz2 == nz-1 && n-dim.i-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n-dim.i-dim.i*dim.j].push_back(bId);
        }
        // TYPE 3 : 8 corners
        if (true)
        {
            nx2 = (int) floor( xx+aOffsetX );         // plus x, plus y, plus z
            ny2 = (int) floor( yy+aOffsetY );
            nz2 = (int) floor( zz+aOffsetZ );
            if (nx2 == nx+1 && ny2 == ny+1 && nz2 == nz+1 && n+1+dim.i+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n+1+dim.i+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // plus x, plus y, minus z
            if (nx2 == nx+1 && ny2 == ny+1 && nz2 == nz-1 && n+1+dim.i-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n+1+dim.i-dim.i*dim.j].push_back(bId);
            ny2 = (int) floor( yy-aOffsetY );         // plus x, minus y, plus z
            nz2 = (int) floor( zz+aOffsetZ );
            if (nx2 == nx+1 && ny2 == ny-1 && nz2 == nz+1 && n+1-dim.i+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n+1-dim.i+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // plus x, minus y, minus z
            if (nx2 == nx+1 && ny2 == ny-1 && nz2 == nz-1 && n+1-dim.i-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n+1-dim.i-dim.i*dim.j].push_back(bId);
            // ----------------------------------------------------------
            nx2 = (int) floor( xx-aOffsetX );         // minus x, plus y, plus z
            ny2 = (int) floor( yy+aOffsetY );
            nz2 = (int) floor( zz+aOffsetZ );
            if (nx2 == nx-1 && ny2 == ny+1 && nz2 == nz+1 && n-1+dim.i+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n-1+dim.i+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // minus x, plus y, minus z
            if (nx2 == nx-1 && ny2 == ny+1 && nz2 == nz-1 && n-1+dim.i-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n-1+dim.i-dim.i*dim.j].push_back(bId);
            ny2 = (int) floor( yy-aOffsetY );         // minus x, minus y, plus z
            nz2 = (int) floor( zz+aOffsetZ );
            if (nx2 == nx-1 && ny2 == ny-1 && nz2 == nz+1 && n-1-dim.i+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n-1-dim.i+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // minus x, minus y, minus z
            if (nx2 == nx-1 && ny2 == ny-1 && nz2 == nz-1 && n-1-dim.i-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n-1-dim.i-dim.i*dim.j].push_back(bId);
        }
    }
    voxMap.set_mapBeadIds(mapBeadIds);
    voxMap.set_mapBeadIdsHalo(mapBeadIdsHalo);

    // printing
    bool flagPrint = false;
    if (flagPrint)
    {
        for (unordered_map<int, vector<int>>::iterator it = mapBeadIds.begin(); it != mapBeadIds.end(); it ++)
        {
            int boxId = it->first;
            vector<int> vecBeadIds = it->second;
            cout << "boxId = " << boxId << " contains beadIds: ";
            for (vector<int>::iterator it2 = vecBeadIds.begin(); it2 != vecBeadIds.end(); it2 ++)
                cout << *it2 << "\t";

            unordered_map<int, vector<int>>::const_iterator got = mapBeadIdsHalo.find(boxId);
            if ( got != mapBeadIdsHalo.end() )
            {
                cout << "with Halo beadIds: ";
                for (vector<int>::iterator it3 = mapBeadIdsHalo[boxId].begin(); it3 != mapBeadIdsHalo[boxId].end(); it3 ++)
                    cout << *it3 << "\t";
            }

            cout << "\n" << endl;
        }
    }
}
// dynamics
void updateCoord(vector<Node> & vecNode, vector<int> & vecBeadIdsCentromere)
{
	
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
    {
        dVec coord = it->get_coord(), veloc = it->get_veloc(), accel = it->get_accel();
        dVec coord_new;
        coord_new = { coord.x + veloc.x*DT + 0.5*accel.x*DT*DT,
                      coord.y + veloc.y*DT + 0.5*accel.y*DT*DT,
                      coord.z + veloc.z*DT + 0.5*accel.z*DT*DT };

        if (SIM_TYPE == "RECONSTRUCT")
        {
            //it->set_coord(coord_new);

            // check confinement by rigid nucleus
            double dist2centerSq = coord_new.x*coord_new.x + coord_new.y*coord_new.y + coord_new.z*coord_new.z;
            dVec center_virtual = {0, 0, -3*RAD_NUCLEUS};
            double dist2centerSq_virtual = (coord_new.x-center_virtual.x)*(coord_new.x-center_virtual.x) \
										 + (coord_new.y-center_virtual.y)*(coord_new.y-center_virtual.y) \
										 + (coord_new.z-center_virtual.z)*(coord_new.z-center_virtual.z);
            int beadId = it - vecNode.begin();
            if (dist2centerSq < RAD_NUCLEUS*RAD_NUCLEUS)
            {
				//if (  find(vecBeadIdsCentromere.begin(), vecBeadIdsCentromere.end(), beadId) != vecBeadIdsCentromere.end())	// this is centromere
					//continue;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//no centromere attached by commenting the above if statement
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				if (dist2centerSq_virtual < RAD_NUCLEUS*RAD_NUCLEUS)	// close to nucleolus
					continue;
				else
					it->set_coord(coord_new);	
			}

        }

        if (SIM_TYPE == "SIMULATION")
        {
            //it->set_coord(coord_new);

            // check confinement by rigid nucleus
            double dist2centerSq = coord_new.x*coord_new.x + coord_new.y*coord_new.y + coord_new.z*coord_new.z;
            dVec center_virtual = {0, 0, -3*RAD_NUCLEUS};
            double dist2centerSq_virtual = (coord_new.x-center_virtual.x)*(coord_new.x-center_virtual.x) \
										 + (coord_new.y-center_virtual.y)*(coord_new.y-center_virtual.y) \
										 + (coord_new.z-center_virtual.z)*(coord_new.z-center_virtual.z);
            int beadId = it - vecNode.begin();
            if (dist2centerSq < RAD_NUCLEUS*RAD_NUCLEUS)
            {
				it->set_coord(coord_new);	
			}

        }
    }
}
// --- following functions are for RECONSTRUCTION objective ---
void updateAccel(vector<Chromosome> & vecChromosome, vector<Node> & vecNode,
                 vector<iPair> vecBeadPairInteract, vector<double> & vecBeadPairInteractFreq,
                 VoxMap & voxMap, unordered_map<string, double> & energy, bool FLAG_INTERACT_ON)
{
    // (0) clear accel values in Nodes
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
    {
        it->set_accelPrev(it->get_accel()); // save previous accel for verlet integration
        it->set_accel({0, 0, 0});
    }

    // (1) beads connected at nearby genomic site
    vector<int> vecBeadIdLocked;
    if (true)   // interaction between beads close in genomic position; note that the Hi-C map includes the nearest interactions!
    {
        for (vector<Chromosome>::iterator it = vecChromosome.begin(); it != vecChromosome.end(); it++)
        {
            vector<int> vecBeadIds = it->get_vecBeadIds();  // temporary
            // bi-node interaction (harmonic spring)
            for (int i = 0; i < vecBeadIds.size()-1; i ++)
            {
                int bId1 = vecBeadIds[i], bId2 = vecBeadIds[i+1];   // consecutive ids of Bead
                if (i == 0)
                    vecBeadIdLocked.push_back(bId1);
                dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
                dVec accel1 = vecNode[bId1].get_accel(), accel2 = vecNode[bId2].get_accel(), accel1_new, accel2_new;
                double L0 = 2.*RAD_BEAD, L = dVecDist(coord1, coord2), L_inv;   // nm
                double k = K_NEAR;  // nN/nm     (N = kg*m/s/s;    nN = kg*nm/s/s)
                double m_inv = 1./MASS_BEAD;    // 1/kg, ~ 1E+23
                dVec unitVec12;
                if (L > 0)
                    L_inv = 1./L;   // 1/nm
                else if (L == 0)
                    L_inv = 10.;
                unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
                double a = k*(L-L0)*m_inv;  //  nm /s/s
                //cout << "a = " << a << endl;
                double ax = a * unitVec12.x, ay = a * unitVec12.y, az = a * unitVec12.z;

                // !!! Assume that 1st bead of each chromosome is locked in position !!!
                //if (i != 0)
                if (true)
                {
                    accel1_new = {  accel1.x + ax, accel1.y + ay, accel1.z + az };
                    vecNode[bId1].set_accel(accel1_new);
                }
                accel2_new = {  accel2.x - ax, accel2.y - ay, accel2.z - az };
                vecNode[bId2].set_accel(accel2_new);

                double ePot12 = 0.5*k*(L-L0)*(L-L0);
                energy["potential"] += ePot12;
                energy["total"]     += ePot12;
            }
            // tri-node interaction (confine angle)
            if (false)
            {
                for (int i = 0; i < vecBeadIds.size()-2; i ++)
                {
                    int bId1 = vecBeadIds[i], bId2 = vecBeadIds[i+2];   // consecutive ids of Bead
                    if (i == 0)
                        vecBeadIdLocked.push_back(bId1);
                    dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
                    dVec accel1 = vecNode[bId1].get_accel(), accel2 = vecNode[bId2].get_accel(), accel1_new, accel2_new;
                    double angle = PI*2./3.;
                    double L0 = 2.*RAD_BEAD*sin(angle), L = dVecDist(coord1, coord2), L_inv;   // nm
                    double k = K_NEAR;  // nN/nm     (N = kg*m/s/s;    nN = kg*nm/s/s)
                    double m_inv = 1./MASS_BEAD;    // 1/kg, ~ 1E+23
                    dVec unitVec12;
                    if (L > 0)
                        L_inv = 1./L;   // 1/nm
                    else if (L == 0)
                        L_inv = 10.;
                    unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
                    double a = k*(L-L0)*m_inv;  //  nm /s/s
                    //cout << "a = " << a << endl;
                    double ax = a * unitVec12.x, ay = a * unitVec12.y, az = a * unitVec12.z;

                    // !!! Assume that 1st bead of each chromosome is locked in position !!!
                    //if (i != 0)
                    if (true)
                    {
                        accel1_new = {  accel1.x + ax, accel1.y + ay, accel1.z + az };
                        vecNode[bId1].set_accel(accel1_new);
                    }
                    accel2_new = {  accel2.x - ax, accel2.y - ay, accel2.z - az };
                    vecNode[bId2].set_accel(accel2_new);

                    double ePot12 = 0.5*k*(L-L0)*(L-L0);
                    energy["potential"] += ePot12;
                    energy["total"]     += ePot12;
                }
            }
        }

    }

    // (2) beads interacting according to Hi-C map
    if (FLAG_INTERACT_ON)   // interaction between beads close in physical position but not in genomic position
    {
        double freqMaxHiC_inv = 1./ *max_element(vecBeadPairInteractFreq.begin(), vecBeadPairInteractFreq.end());
        for (vector<iPair>::iterator it = vecBeadPairInteract.begin(); it != vecBeadPairInteract.end(); it++)
        {
            int pos = it-vecBeadPairInteract.begin();
            int bId1 = (*it).i, bId2 = (*it).j;
            int hostChromoId1 = vecNode[bId1].get_hostChromoId();
            int hostChromoId2 = vecNode[bId2].get_hostChromoId();
            
            //---------- this mean we will have only INTRA-chromosomal contacts --------------
            if (hostChromoId1 == hostChromoId2) // this limit the following calculations to intra-chromosomal interactions
            //-----------------------------------------------------
            if (true)
            {
                //cout << "beadId: " << bId1 << ", " << bId2 << endl;
                //cout << "chromoID: " << hostChromoId1 << ", " << hostChromoId2;
                //cout << " ... " << endl;
                dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
                dVec accel1 = vecNode[bId1].get_accel(), accel2 = vecNode[bId2].get_accel(), accel1_new, accel2_new;
                double L0 = 2.*RAD_BEAD, L = dVecDist(coord1, coord2), L_inv;   // nm

                double k = K_INTRA;  // nN/nm     (N = kg*m/s/s;    nN = kg*nm/s/s)
                if (hostChromoId1 == hostChromoId2)
					k = K_INTRA;
                if (hostChromoId1 != hostChromoId2)
                    k = K_INTER;
                k = k*freqMaxHiC_inv*vecBeadPairInteractFreq[pos];  // linearly relate k to the frequency in HiC map

                double m_inv = 1./MASS_BEAD;    // 1/kg, ~ 1E+23
                dVec unitVec12;
                if (L > 0)
                    L_inv = 1./L;   // 1/nm
                else if (L == 0)
                    L_inv = 10.;
                unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
                double a = k*(L-L0)*m_inv;  //  nm /s/s
                //cout << "a = " << a << endl;
                double ax = a * unitVec12.x, ay = a * unitVec12.y, az = a * unitVec12.z;

                //auto findLocked = find(vecBeadIdLocked.begin(), vecBeadIdLocked.end(), bId1);   // temporary
                //if (findLocked == vecBeadIdLocked.end())
                if (true)
                {
                    accel1_new = {  accel1.x + ax, accel1.y + ay, accel1.z + az };
                    vecNode[bId1].set_accel(accel1_new);
                }
                //auto findLocked2 = find(vecBeadIdLocked.begin(), vecBeadIdLocked.end(), bId1);  // temporary
                //if (findLocked2 == vecBeadIdLocked.end())
                if (true)
                {
                    accel2_new = {  accel2.x - ax, accel2.y - ay, accel2.z - az };
                    vecNode[bId2].set_accel(accel2_new);
                }

                double ePot12 = 0.5*k*(L-L0)*(L-L0);
                energy["potential"] += ePot12;
                energy["total"]     += ePot12;
            }
        }
        //exit(123);
    }

    // (3) beads volume exclusion / collision detection
    if (true)   // interaction between beads that collide (physically too close)
    {
        iVec dim = voxMap.get_dim();
        int nMax = dim.i*dim.j*dim.k;
        unordered_map<int, vector<int>> mapBeadIds = voxMap.get_mapBeadIds(); // use reference to save memory
        unordered_map<int, vector<int>> mapBeadIdsHalo = voxMap.get_mapBeadIdsHalo();
        for (int boxId = 0; boxId < nMax; boxId ++)
        {
            vector<int> vecBeadIds0, vecBeadIds1;
            unordered_map<int, vector<int>>::const_iterator got = mapBeadIds.find(boxId);
            if ( got != mapBeadIds.end() )  // this box has Bead inside
            {
                for (vector<int>::iterator it = mapBeadIds[boxId].begin(); it != mapBeadIds[boxId].end(); it ++)
                    vecBeadIds0.push_back(*it);
            }
            unordered_map<int, vector<int>>::const_iterator got2 = mapBeadIdsHalo.find(boxId);
            if ( got2 != mapBeadIdsHalo.end() )  // this box has Bead inside
            {
                for (vector<int>::iterator it = mapBeadIdsHalo[boxId].begin(); it != mapBeadIdsHalo[boxId].end(); it ++)
                    vecBeadIds1.push_back(*it);
            }

            // STEP 1 : repulsion between Beads within the Box
            int nBeadInBox = vecBeadIds0.size();
            if (nBeadInBox > 1)
            {
                for (int l = 0; l < nBeadInBox-1; l++)
                {
                    for (int m = l+1; m < nBeadInBox; m++)
                    {
                        int bId1 = vecBeadIds0[l], bId2 = vecBeadIds0[m];
                        dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
                        double L0 = 2.*RAD_BEAD, L = dVecDist(coord1, coord2);   // nm
                        if (L < L0) // repulsion only when overlapping of Beads occurs
                        {
                            double L_inv;
                            double m_inv = 1./MASS_BEAD;    // 1/kg, ~ 1E+23
                            dVec unitVec12;
                            if (L > 0)
                                L_inv = 1./L;   // 1/nm
                            else if (L == 0)
                                L_inv = 10.;
                            double k = K_INTRA;
                            dVec accel1 = vecNode[bId1].get_accel(), accel2 = vecNode[bId2].get_accel(), accel1_new, accel2_new;
                            unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
                            double a = k*(L-L0)*m_inv;  //  nm /s/s
                            double ax = a * unitVec12.x, ay = a * unitVec12.y, az = a * unitVec12.z;

                            accel1_new = {  accel1.x + ax, accel1.y + ay, accel1.z + az };
                            vecNode[bId1].set_accel(accel1_new);

                            accel2_new = {  accel2.x - ax, accel2.y - ay, accel2.z - az };
                            vecNode[bId2].set_accel(accel2_new);
                        }
                    }
                }
            }
            // STEP 2 : repulsion between one Bead within the Box and another Bead within Halo
            int nBeadInHalo = vecBeadIds1.size();
            if (nBeadInBox > 0 && nBeadInHalo > 0)
            {
                for (int l = 0; l < nBeadInBox-1; l++)
                {
                    for (int m = 0; m < nBeadInHalo; m++)
                    {
                        int bId1 = vecBeadIds0[l], bId2 = vecBeadIds1[m];
                        dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
                        double L0 = 2.*RAD_BEAD, L = dVecDist(coord1, coord2);   // nm
                        if (L < L0) // repulsion only when overlapping of Beads occurs
                        {
                            double L_inv;
                            double m_inv = 1./MASS_BEAD;    // 1/kg, ~ 1E+23
                            dVec unitVec12;
                            if (L > 0)
                                L_inv = 1./L;   // 1/nm
                            else if (L == 0)
                                L_inv = 10.;
                            double k = K_INTRA;
                            dVec accel1 = vecNode[bId1].get_accel(), accel2 = vecNode[bId2].get_accel(), accel1_new, accel2_new;
                            unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
                            double a = k*(L-L0)*m_inv;  //  nm /s/s
                            double ax = a * unitVec12.x, ay = a * unitVec12.y, az = a * unitVec12.z;

                            accel1_new = {  accel1.x + ax, accel1.y + ay, accel1.z + az };
                            vecNode[bId1].set_accel(accel1_new);
                        }
                    }
                }
            }
        }
    }


}
void updateVeloc(vector<Node> & vecNode, unordered_map<string, double> & energy)
{
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
    {
        dVec veloc = it->get_veloc(), accel = it->get_accel(), accelPrev = it->get_accelPrev();
        dVec veloc_new;
        veloc_new = { veloc.x + 0.5*(accel.x+accelPrev.x)*DT,
                      veloc.y + 0.5*(accel.y+accelPrev.y)*DT,
                      veloc.z + 0.5*(accel.z+accelPrev.z)*DT };
        it->set_veloc(veloc_new);

        // summ up kinetic energy
        double speedSq = veloc.x*veloc.x + veloc.y*veloc.y + veloc.z*veloc.z;
        //double speedSq = veloc_new.x*veloc_new.x + veloc_new.y*veloc_new.y + veloc_new.z*veloc_new.z;
        //double eKin = 0.5*MASS_BEAD*speedSq;
        double eKin = 0.5*speedSq;  // kinetic energy per unit mass
        energy["kinetic"] += eKin;
        //cout << "kinetic energy: add = " << eKin << " , sum = " << energy["kinetic"] << endl;
        energy["total"]   += eKin;
    }
}
void updateVelocLinearDamp(vector<Node> & vecNode, unordered_map<string, double> & energy)
{
    // refer to modified verlet method ( equation (34) in cmotion.pdf )
    // STEP 3 : v(t+dt) = ( v(t)*( 1 - 0.5*dt*gamma/m )   0.5* ( a(t)+a(t+dt) ) * dt ) / (1 + 0.5*dt*gamma/m)
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
    {
        dVec veloc = it->get_veloc(), accel = it->get_accel(), accelPrev = it->get_accelPrev();
        dVec veloc_new;
        double GAMMA = 1E-20;
        double factor = 0.5*DT*GAMMA/MASS_BEAD;
        veloc_new = { ((1-factor)*veloc.x + 0.5*(accel.x+accelPrev.x)*DT)/(1+factor),
                      ((1-factor)*veloc.y + 0.5*(accel.y+accelPrev.y)*DT)/(1+factor),
                      ((1-factor)*veloc.z + 0.5*(accel.z+accelPrev.z)*DT)/(1+factor) };
        it->set_veloc(veloc_new);

        double speedSq = veloc.x*veloc.x + veloc.y*veloc.y + veloc.z*veloc.z;
        //double eKin = 0.5*MASS_BEAD*speedSq;
        double eKin = 0.5*speedSq;
        energy["kinetic"] += eKin;
        energy["total"]   += eKin;
    }
}

double dVecDist(dVec & p1, dVec & p2)
{
    return sqrt( (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y) + (p1.z-p2.z)*(p1.z-p2.z) );
}
