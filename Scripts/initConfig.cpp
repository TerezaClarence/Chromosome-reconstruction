/*
    File: initConfig.cpp
    Function: Initialize chromosomes
    Model: chromoCell
    Created: 7 December, 2017 (XF)
*/

#include "initConfig.hpp"

// ===================== Classes =========================

/* -------------- Node ---------------- */
Node::Node() {}
Node::Node(int id, int hostChromoId)
{
    this->id = id;
    this->hostChromoId = hostChromoId;
}
Node::Node(int id, int hostChromoId, dVec coord, dVec veloc, dVec accel, dVec accelPrev)
{
    this->id = id;
    this->hostChromoId = hostChromoId;
    this->coord = coord;
    this->veloc = veloc;
    this->accel = accel;
    this->accelPrev = accelPrev;
}
Node::~Node() {}

/* -------------- Chromosome ---------------- */
Chromosome::Chromosome() {}
Chromosome::Chromosome(int id, vector<int> vecBeadIds)
{
    this->id = id;
    this->vecBeadIds = vecBeadIds;
}
Chromosome::Chromosome(int id, vector<int> vecBeadIds, vector<int> vecGeneIds)
{
    this->id = id;
    this->vecBeadIds = vecBeadIds;
    this->vecGeneIds = vecGeneIds;
}
Chromosome::~Chromosome() {}

/* -------------- VoxMap ---------------- */
VoxMap::VoxMap() {}
VoxMap::VoxMap(iVec dim, unordered_map<int, vector<int>> mapBeadIds, unordered_map<int, vector<int>> mapBeadIdsHalo)
{
    this->dim = dim;
    this->mapBeadIds = mapBeadIds;
    this->mapBeadIdsHalo = mapBeadIdsHalo;
}
VoxMap::~VoxMap() {}

// ===================== Functions =========================
// --- following functions are for RECONSTRUCTION objective ---
void initChromosome(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, \
                    vector<iPair> & vecBeadPairInteract, vector<double> & vecBeadPairInteractFreq, int argc, char ** argv)
{
    double **matChromoContFreq = NULL;
    string iFilename_normMatrix, iFilename_numberBead;
    if (argc == 1)
    {
        cout << "\nPlease input filename (without suffix) as argument! \n";
        exit(0);
    }
    if (argc == 2)
    {
        iFilename_normMatrix = argv[1]+SUFFIX_normMatrix;
        iFilename_numberBead = argv[1]+SUFFIX_numberBead;

        cout << ">> reading file " << iFilename_numberBead << " ... "<< endl;
        ifstream file_numberBead (iFilename_numberBead);
        string line_numberBead;
        if (file_numberBead.is_open())
        {
            getline (file_numberBead, line_numberBead); // skip the first line
            while (getline (file_numberBead, line_numberBead))
            {
                stringstream ss (line_numberBead);
                vector<string> items;
                string buf;
                while (ss >> buf)
                    items.push_back(buf);
                //cout << line_numberBead << endl;
                //cout << items[0] << ", " << items[1] << endl;   // items[1] is chromosome id, items[2] is number of beads

                /* ---------- create Node and Chromosome ----------- */
                int chromoId, nBead;
                chromoId = stoi(items[0]);
                nBead    = stoi(items[1]);
                //cout << chromoId << ", " << nBead << endl;

                int currBeadId = vecNode.size();
                vector<int> vecBeadIds;
                for (int i = 0; i < nBead; i ++)
                {
                    Node bead0 = Node(currBeadId, chromoId);
                    vecNode.push_back(bead0);
                    vecBeadIds.push_back(currBeadId);
                    currBeadId ++;
                }
                Chromosome chromo0 = Chromosome(chromoId, vecBeadIds);
                vecChromosome.push_back(chromo0);
            }
            file_numberBead.close();
            cout << "    done!" << endl;
        }
        else
        {
            cout << "Unable to open file " << iFilename_numberBead << endl;
            exit(1);
        }

        // allocate memory to contact frequency map
        const int nBeadTotal = vecNode.size();
       /* cout << ">> allocating memory to temporary matChromoContFreq with dimension " << nBeadTotal << " by " << nBeadTotal << endl;
        matChromoContFreq = new double* [nBeadTotal];
        for (int i = 0; i<nBeadTotal; i++)
            matChromoContFreq[i] = new double[nBeadTotal];
        for (int j = 0; j<nBeadTotal; j++)
        {
            for (int k = 0; k<nBeadTotal; k++)
                matChromoContFreq[j][k] = 0;
        }
*/
        cout << "    done!" << endl;

        cout << ">> reading matrix file " << iFilename_normMatrix << " ... "<< endl;
        ifstream file_normMatrix (iFilename_normMatrix);
        string line_normMatrix;
        if (file_normMatrix.is_open())
        {
            //getline (file_normMatrix, line_normMatrix); // skip the first line
            int currRow = 0;
            double freqMax = 0;
            int cntOverCutoff = 0;
            while(getline (file_normMatrix, line_normMatrix))
            {
                stringstream ss (line_normMatrix);
                vector<string> items;
                string buf;
                while (ss >> buf)
                    items.push_back(buf);
                /* --------- initiate contact frequency matrix --------- */
               
					double freq = stod (items[2]);			//double freq = stod (*it);
					int currRow = stod (items[0]);		//int currCol = it-items.begin();
                    int currCol = stod(items[1]);
     //               matChromoContFreq[currRow][currCol] = freq;
                    if (freq > freqMax)
                        freqMax = freq;
                    if (freq != 0 && freq >= CUTOFF_NORMMAT && currCol > currRow + GENOMIC_LEAST_SEP_INTRA)
                    {
                        cntOverCutoff ++;
                        cout << "Bead id pairs with interaction: " << currRow << " " << currCol << endl;
                        vecBeadPairInteract.push_back({currRow, currCol});
                        vecBeadPairInteractFreq.push_back(freq);
                    }
                currRow ++;
            }
            cout << "freqMax = " << freqMax << "; " << cntOverCutoff << " over cutoff (" << CUTOFF_NORMMAT << ")" << endl;
            //exit(88);
            file_normMatrix.close();

            cout << "    done!" << endl;
        }
        else
        {
            cout << "Unable to open file " << iFilename_normMatrix << endl;
            exit(1);
        }
    }
}
void initDynamics(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<int> & vecBeadIdsCentromere, int argc, char ** argv,
                    string initType, vector<double> & randNum01)
{
    int k = 0;  // position in randNum01
    cout << ">> initializing bead coordinates with *" << initType << "* method ... " << endl;
    if (initType == "reconstruct")
    {

    }
    else if (initType == "random")
    {
        double rInit = 0.75*RAD_NUCLEUS;
        double r, theta, phi;
        for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it ++)
        {
            // ref1: https://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability
            // ref2: http://corysimon.github.io/articles/uniformdistn-on-sphere/
            dVec coord;
            int beadId = it - vecNode.begin();
            if (  find(vecBeadIdsCentromere.begin(), vecBeadIdsCentromere.end(), beadId) != vecBeadIdsCentromere.end())	// this is centromere
            {
				coord = {0, 0, RAD_NUCLEUS};
			}
            else
            {
				r     = rInit * cbrt(randNum01[k++]);
				theta = acos(2.* randNum01[k++] - 1);
				phi   = 2.* PI * randNum01[k++];	
				coord = {r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta)};
			}
            it->set_coord(coord);
        }
    }
    cout << "    done!" << endl;

    cout << ">> initializing bead velocities & accelerations ... " << endl;
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it ++)
    {
        //double s = 10.; // nm/s
        //dVec veloc = {s*randNum01[k++], s*randNum01[k++], s*randNum01[k++]};
        dVec veloc = {0., 0., 0.};
        dVec accel = {0., 0., 0.}, accelPrev = {0., 0., 0.};
        it->set_veloc(veloc);
        it->set_accel(accel);
        it->set_accelPrev(accelPrev);
    }
    cout << "    done!" << endl;
}
// ------------------------------------------------------------
void printChromoInfo(vector<Chromosome> & vecChromosome, vector<Node> & vecNode)
{
    cout << "Chromosome_id\tBead_id\tSpeed_avg(nm/s)\tSpeed_max(nm/s)" << endl;
    for (vector<Chromosome>::iterator it_c = vecChromosome.begin(); it_c != vecChromosome.end(); it_c++)
    {
        vector<int> vecBeadIds = it_c->get_vecBeadIds();
        int chromoId = it_c->get_id();

        int beadIdM;
        dVec coordM, velocM, accelM;
        double speedSqM = -1.;
        double speedSqAvg = 0;
        double speedAvg = 0;
        for (vector<int>::iterator it_n = vecBeadIds.begin(); it_n != vecBeadIds.end(); it_n++)
        {
            int beadId = *it_n;
            dVec veloc = vecNode[beadId].get_veloc();
            double speedSq = veloc.x*veloc.x+veloc.y*veloc.y+veloc.z*veloc.z;
            //speedSqAvg += speedSq;
            speedAvg   += sqrt(speedSq);
            if ( speedSq > speedSqM )
            {
                speedSqM = speedSq;
                beadIdM = beadId;
                coordM = vecNode[beadId].get_coord();
                velocM = veloc;
                accelM = vecNode[beadId].get_accel();
            }
        }

        //cout << chromoId << "\t" << beadIdM << "\t" << sqrt(speedSqAvg/vecBeadIds.size()) << "\t" << sqrt(speedSqM) << endl;
        cout << chromoId << "\t" << beadIdM << "\t" << speedAvg/vecBeadIds.size() << "\t" << sqrt(speedSqM) << endl;
        /*
        cout << "Bead id = " << beadIdM << ": c = {" << coordM.x << ", " << coordM.y << ", " << coordM.z << "}; "
                                        <<  " v = {" << velocM.x << ", " << velocM.y << ", " << velocM.z << "}; "
                                        <<  " a = {" << accelM.x << ", " << accelM.y << ", " << accelM.z << "}; "
                                        << endl;
        */
    }
}

void writeChromoInfo(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<int> & vecCondensinSite, vector<int> & vecBeadIdsCentromere, double t_now)
{
    cout << ">> writing bead information at t_now = " << t_now << " ... " << endl;

    // write in x, y, z format (maybe with color information)
    cout << ">> >> coordinate & color ..." << endl;
    ofstream beadFile;
    string beadFileName;
    stringstream ss;
    ss << fixed << setprecision(1) << t_now;
    beadFileName = "_bead_coord_color_" + ss.str() + ".txt";
    if (t_now == 0)
        beadFile.open(beadFileName, ios::out | ios::trunc); beadFile.close();

    // write in PDB format
    // ....... format ........ see https://www.ichemlabs.com/166
    // ATOM   {atom_id as global bead id} CA  ASP {chain_id as A, B, C, ...}  {residue_id as bead_id in Chromosome} {x} {y} {z} 1.00 0.00
    // TER  {atom_id_prev + 1}  {empty} ASP {chain_id}  {residue_id}
    bool flagWritePDB = true;
    if (flagWritePDB)   // has an error about large residue number on chain
    {
        cout << ">> >> PDB file ..." << endl;
        ofstream beadFilePDB;
        string beadFileNamePDB;
        stringstream ssPDB;
        ssPDB << fixed << setprecision(1) << t_now;
        beadFileNamePDB = "_chromoPDB_" + ss.str() + ".pdb";
        if (t_now == 0)
            beadFilePDB.open(beadFileNamePDB, ios::out | ios::trunc); beadFilePDB.close();

        beadFile.open(beadFileName, ios::app | ios::binary);
        beadFilePDB.open(beadFileNamePDB, ios::app | ios::binary);
        int nChromo = vecChromosome.size();
        unordered_map<int, string> int2letter;
        int2letter[1] = "A"; int2letter[2] = "B"; int2letter[3] = "C"; int2letter[4] = "D"; int2letter[5] = "E"; int2letter[6] = "F";
        int2letter[7] = "G"; int2letter[8] = "H"; int2letter[9] = "I"; int2letter[10] = "J"; int2letter[11] = "K"; int2letter[12] = "L";
        int2letter[13] = "M"; int2letter[14] = "N"; int2letter[15] = "O"; int2letter[16] = "P"; int2letter[17] = "Q"; int2letter[18] = "R";
        for (vector<Chromosome>::iterator it_c = vecChromosome.begin(); it_c != vecChromosome.end(); it_c++)
        {
            vector<int> vecBeadIds = it_c->get_vecBeadIds();
            int chromoId = it_c->get_id();
            double cx = (float)chromoId/nChromo, cz = 0, cy = 1-cx;

            int chromoSegLength = 9999; // this is used for visualization purpose because the length is too long for nucleosome level simulation
            int segId = 0;
            for (vector<int>::iterator it_n = vecBeadIds.begin(); it_n != vecBeadIds.end(); it_n++)
            {
                int beadId = *it_n;
                bool flagIsCondensinSite = false, flagIsCentromere = false;
                if (SIM_TYPE == "SIMULATION" && find(vecCondensinSite.begin(), vecCondensinSite.end(), beadId) != vecCondensinSite.end())
                    flagIsCondensinSite = true;
                if (SIM_TYPE == "RECONSTRUCT" && find(vecBeadIdsCentromere.begin(), vecBeadIdsCentromere.end(), beadId) != vecBeadIdsCentromere.end())
					flagIsCentromere = false;

                dVec coord = vecNode[beadId].get_coord();
                double x = coord.x, y = coord.y, z = coord.z;

                beadFile <<  x << "\t" <<  y << "\t" <<  z << "\t"
                         << cx << "\t" << cy << "\t" << cz << endl;

                stringstream ss0, ss1;
                stringstream ssx, ssy, ssz;
                ss0 << setw(5) << beadId+chromoId;
                ss1 << setw(4) << 1+(it_n-vecBeadIds.begin());
                ssx << fixed << setw(8) << setprecision(3) << x;
                ssy << fixed << setw(8) << setprecision(3) << y;
                ssz << fixed << setw(8) << setprecision(3) << z;
                if (SIM_TYPE == "RECONSTRUCT")
                {
					if (flagIsCentromere == true)	// this is centromere
						beadFilePDB << "ATOM  " << ss0.str() << "  CA" << "  " << "CYS" << " " << int2letter[chromoId] << ss1.str() << "    "
								<< ssx.str() << ssy.str() << ssz.str() << "  1.00" << "  0.00" << endl;	
					else
						beadFilePDB << "ATOM  " << ss0.str() << "  CA" << "  " << "ASP" << " " << int2letter[chromoId] << ss1.str() << "    "
								<< ssx.str() << ssy.str() << ssz.str() << "  1.00" << "  0.00" << endl;	
				}
           }
            stringstream ss2;
            ss2 << setw(5) << vecBeadIds[vecBeadIds.size()-1]+1+chromoId;
            beadFilePDB << "TER   " << ss2.str() << "    " << "  " << "ASP" << " " << int2letter[chromoId] << endl;
        }
        beadFile.close();
        beadFilePDB.close();
    }

    // write pairwise distance
    bool flagWriteDistMat = false;
    if (flagWriteDistMat)
    {
        cout << ">> >> bead pairwise distance ..." << endl;
        ofstream beadDistFile;
        string beadDistFileName;
        stringstream ss3;
        ss3 << fixed << setprecision(1) << t_now;
        beadDistFileName = "_bead_dist_matrix_" + ss3.str() + ".txt";
        if (t_now == 0)
            beadDistFile.open(beadDistFileName, ios::out | ios::trunc); beadDistFile.close();
        beadDistFile.open(beadDistFileName, ios::app | ios::binary);

        int nBeadTotal = vecNode.size();
        double **matBeadDist = NULL;
        matBeadDist = new double* [nBeadTotal];
        for (int i = 0; i<nBeadTotal; i++)
            matBeadDist[i] = new double[nBeadTotal];
        for (int j = 0; j<nBeadTotal; j++)
        {
            for (int k = 0; k<nBeadTotal; k++)
                matBeadDist[j][k] = 0;
        }
        cout << "    ____calculating____" << endl;
        int cnt = 0;
        for (int i = 0; i < nBeadTotal; i ++)
        {
            for (int j = i+1; j < nBeadTotal; j ++)
            {
                dVec coord1 = vecNode[i].get_coord();
                dVec coord2 = vecNode[j].get_coord();
                matBeadDist[i][j] = sqrt( (coord1.x-coord2.x)*(coord1.x-coord2.x) + (coord1.y-coord2.y)*(coord1.y-coord2.y) + (coord1.z-coord2.z)*(coord1.z-coord2.z) );
                matBeadDist[j][i] = matBeadDist[i][j];
                cnt ++;
                if (cnt % ((int) nBeadTotal*nBeadTotal/40) == 0)
                    cout << "     ... " << (int) (cnt*200./nBeadTotal/nBeadTotal) << "\% finished ..." << endl;
            }
        }
        cout << "    ____writing____" << endl;
        int cnt2 = 0;
        for (int i = 0; i < nBeadTotal; i ++)
        {
            for (int j = 0; j < nBeadTotal; j ++)
            {
                beadDistFile <<  matBeadDist[i][j] << "\t";
                cnt2 ++;
                if (cnt2 % ((int) nBeadTotal*nBeadTotal/40) == 0)
                    cout << "     ... " << (int) (cnt2*200./nBeadTotal/nBeadTotal) << "\% finished ..." << endl;
            }
            beadDistFile << endl;
        }
        beadDistFile.close();
    }

    cout << "    done!" << endl;
}
