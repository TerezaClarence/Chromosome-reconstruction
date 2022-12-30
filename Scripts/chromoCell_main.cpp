#include "initConfig.hpp"
#include "dynamics.hpp"

int main(int argc, char** argv)
{
    // Set up random number generators
    int seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 generator (seed);
    uniform_real_distribution<double> uniform01(0.0, 1.0);

    // initialize variables
    unordered_map<string, double> energy;
    energy["total"] = 0.;
    energy["kinetic"] = 0.;
    energy["potential"] = 0.;

    // initialize voxMap -- if using VoxMap method
    //const int nSubdiv = 3;
    //const int dimX = pow(2, nSubdiv);
    const int dimX = 10;
    const iVec dim = {dimX, dimX, dimX};
    unordered_map<int, vector<int>> mapBeadIds;
    unordered_map<int, vector<int>> mapBeadIdsHalo;
    VoxMap voxMap = VoxMap(dim, mapBeadIds, mapBeadIdsHalo);

    // initialize vectors
    vector<Chromosome> vecChromosome;           // collect all Chromosome objects
    vector<Node>       vecNode;                 // collect all Node objects (Beads, Vertics)
    vector<int>        vecCondensinSite;        // collect all condensin sites
    //double **matChromoContFreq;                 // contact frequency matrix for intra- and inter-chromosome interactions
    vector<iPair> vecBeadPairInteract;          //
    vector<double> vecBeadPairInteractFreq;

    // initialize chromosomes (beads)
    if (SIM_TYPE == "RECONSTRUCT")
        initChromosome(vecChromosome, vecNode, vecBeadPairInteract, vecBeadPairInteractFreq, argc, argv);  // pass the reference to vectors
	
	// initialize centromeres
    // ... CEN1: 3753687 - 3789421; chr1 offset 567
    // ... CEN2: 1602264 - 1644747; chr2 offset 903
    // ... CEN3: 1070904 - 1137003; chr3 offset 442
    int cen1_start = 3753687, cen1_end = 3789421;
    int cen2_start = 1602264, cen2_end = 1644747;
    int cen3_start = 1070904, cen3_end = 1137003;
    int offset_bp_chr1 = 567;
    int offset_bp_chr2 = 903;
    int offset_bp_chr3 = 442;
    unordered_map<int, int> mapCentromere;	// chromoId -> bead rank in chromo (location_centromere)
    vector<int> vecBeadIdsCentromere;       // contains bead ids that are centromere
    mapCentromere[0] = (cen1_start-offset_bp_chr1 + cen1_end-offset_bp_chr1) / 2 / GENOMIC_RESOLUTION;
    mapCentromere[1] = (cen2_start-offset_bp_chr2 + cen2_end-offset_bp_chr2) / 2 / GENOMIC_RESOLUTION;
    mapCentromere[2] = (cen3_start-offset_bp_chr3 + cen3_end-offset_bp_chr3) / 2 / GENOMIC_RESOLUTION;
    
    
    for (vector<Chromosome>::iterator it = vecChromosome.begin(); it != vecChromosome.end(); it++)
    {
		vector<int> vecBeadIds = it->get_vecBeadIds();
		int chromoId = it - vecChromosome.begin();
		int beadIdCen = vecBeadIds[ mapCentromere[chromoId] ];
		cout << mapCentromere[chromoId] << " (" << vecBeadIds.size() << ") " << endl;
		vecBeadIdsCentromere.push_back(beadIdCen);
		cout << "beadId = " << beadIdCen << " is centromere on chromoId = " << chromoId << endl;
	}
	
    // initialize bead coordinate
    vector<double> randNum01;
    for (int i = 0; i < vecNode.size()*6; i ++)
        randNum01.push_back(uniform01(generator));  // used for coord and veloc initiation
    if (SIM_TYPE == "RECONSTRUCT")
    {
        string initType = "random";
        initDynamics(vecChromosome, vecNode, vecBeadIdsCentromere, argc, argv, initType, randNum01);
    }
    
    
    
    /* ---------- simulation starts here ----------- */

    double t_now = 0;
    long nIter = 0;
    while (t_now <= T)
    {
        // single iteration of simulation
        energy["total"] = 0; energy["kinetic"] = 0; energy["potential"] = 0;
        bool FLAG_INTERACT_ON = false;
        if (nIter >= CNT_START_INTERACT)
            FLAG_INTERACT_ON = true;

        if (SIM_TYPE == "RECONSTRUCT")
            oneIter(vecChromosome, vecNode, vecBeadPairInteract, vecBeadPairInteractFreq, vecBeadIdsCentromere, energy, voxMap, FLAG_INTERACT_ON); // if using voxMap
        
		// print information
        if (nIter % FREQ_PRINT == 0)
        {
            cout << "\n --------- t_now = " << t_now << " ---------" << endl;
            cout << " --------- Energy (x1e-18 J/kg ) = " << energy["total"] << endl;
            cout << " ---------        (potential   ) = " << energy["potential"] << endl;
            cout << " ---------        (kinetic     ) = " << energy["kinetic"] << endl;
            printChromoInfo(vecChromosome, vecNode);
        }
        
		if (energy["total"] < -100 && t_now > 1.0 )
		{
			cout << "Converged! Total energy is " << energy["total"] << endl;
			break;	
		}

        // write coordinate to file
        if (nIter % FREQ_WRITE == 0)
        {
            writeChromoInfo(vecChromosome, vecNode, vecCondensinSite, vecBeadIdsCentromere, t_now);
        }

        t_now += DT;
        nIter ++;

    }

    return 0;
}

