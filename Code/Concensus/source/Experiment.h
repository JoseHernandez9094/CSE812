#ifndef HP_EXPERIMENT_H
#define HP_EXPERIMENT_H

#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <deque>
#include <string>
#include <utility>
#include <fstream>


#include "hp_config.h"
#include "../../Empirical/source/tools/Random.h"
#include "../../Empirical/source/tools/random_utils.h"
#include "../../Empirical/source/hardware/AvidaGP.h"
#include "../../Empirical/source/Evolve/World.h"

///mnt/scratch/herna383/CSE812

using std::cout; using std::endl; using std::unordered_map; using std::vector;
using std::unordered_set; using std::find; using std::deque; using std::pair;
using std::make_pair; using std::ofstream; using std::string;

constexpr size_t SYS_POS = 0;
constexpr size_t SYS_VOTE = 1;
constexpr size_t SYS_ID = 2;
constexpr size_t MAX_MESS = 50;

class Experiment
{
    struct System;
    struct Agent;

    // World of Systems
    using world_t = emp::World<Agent>; 
    // Type of representation
    using hardware_t = emp::AvidaGP;
    // Instruction object for hardware
    using inst_t = hardware_t::inst_t;
    // Instruction Library
    using inst_lib_t = hardware_t::inst_lib_t;
    // Hardware Genome
    using program_t = hardware_t::Genome;
    // Random number generator pointer
    using rng_t = emp::Ptr<emp::Random>;
    // World Pointer
    using wrd_t = emp::Ptr<world_t>; 

    struct Agent
    {
        //Agents score
        double score = 0;
        double time = -1;
        program_t mGenome;

        Agent(const program_t & p) : mGenome(p) {;}
        program_t & GetGenome() {return mGenome;}

        void SetGenome(const program_t & p) {mGenome = p;}
    };

    struct System
    {
        // System score
        double mScore = 0;
        // Systems program
        program_t mGenome;        
        // Nodes in the System
        int N;
        // Maximum Score
        double MAX_SCORE;
        // System with GP system
        vector<emp::Ptr<hardware_t>> system;
        // System Messaging Queue
        vector<deque<int>> message;
        // Adjacency List
        unordered_map<int, vector<int>> adj_list;

        // Genome getter
        program_t & GetGenome() {return mGenome;}
        // Genome setter
        void SetGenome(const program_t p) 
        {
            mGenome = p; 
            for(int i = 0; i < N; ++i)
                system[i]->SetGenome(mGenome);
        }
        // Set Size of vectors
        // dimenstionality of the map
        void SetSize(int n) 
        {
            N = n * n; 
            MAX_SCORE = N + N;
            message.resize(N);

            for(int i = 0; i < N; ++i)
            {
                auto hw = emp::NewPtr<hardware_t>(mGenome);
                system.push_back(hw);
                system[i]->SetTrait(SYS_POS, i);
                system[i]->SetTrait(SYS_VOTE, -1);
            } 
        }
        // Generate Random UIDS for System
        void SetUID(rng_t rng)
        {
            unordered_set<int> ids;
            while(ids.size() != N)
            {
                ids.insert(rng->GetUInt(1000000000,2000000000));
            }
            int i = 0;
            for(auto & id : ids)
            {
                system[i]->SetTrait(SYS_ID, id);
                ++i;
            }
        }
        // Set Adjacency List
        void SetAdj(unordered_map<int, vector<int>> & adj){adj_list = adj;}
        // Calculate score and store it!
        double GetScore()
        {
            double cnt = 0.0;

            vector<int> votes;
            vector<int> UIDS;

            for(int i = 0; i < system.size(); ++i)
            {
                votes.push_back(system[i]->GetTrait(SYS_VOTE));
                UIDS.push_back(system[i]->GetTrait(SYS_ID));
            }

            // Find the vote with the most support!
            // <id, cnt> 
            unordered_map<int, int> score;
            //Go through all the votes the system decided on
            for(auto & v : votes)
            {
                //Check if vote is in the UIDS
                if(find(UIDS.begin(), UIDS.end(), v) != UIDS.end())
                {
                    if(score.find(v) != score.end())
                        ++score[v];
                    else
                        score[v] = 1;
                }
            }

            double low = 0; 
            for(auto & p : score)
            {
                if(p.second > low)
                    low = p.second;
            }
            cnt += low;

            // Count all the legal votes!
            for(auto &v : votes)
            {
                if(find(UIDS.begin(), UIDS.end(), v) != UIDS.end())
                    ++cnt;
            }

            

            if(cnt > 0.0)
                mScore = cnt;

            return mScore;
        }
        // Allow One instruction to be executed
        // Return if consensus is formed
        bool Run(vector<size_t> order)
        {
            for(auto i : order)
            {
                system[i]->SingleProcess();
            }

            if(GetScore() == MAX_SCORE)
            {
                cout << "*******MAX" <<  MAX_SCORE << endl;
                cout << "*******Score" << mScore << endl;
                return true;
            }

            return false;
        }
        // Reset the System
        void Reset()
        {
            for(auto &m : message)
                m.clear();
            mScore = 0;
            for(int i = 0; i < N; ++i)
            {
                system[i]->Reset();
                system[i]->SetTrait(SYS_POS, i);
                system[i]->SetTrait(SYS_VOTE, -1);
            } 
        }
        // Print internal traits
        void PrintTraits()
        {
            for(int i = 0; i < system.size(); ++i)
            {
                if (message[i].size() > 0)
                    cout << "( ID: " << system[i]->GetTrait(SYS_POS) << ", Vote: " << system[i]->GetTrait(SYS_VOTE) << ", UID: " << system[i]->GetTrait(SYS_ID) << ", Message: " << message[i].front()  << ")" << ", " <<endl;
                else
                    cout << "( ID: " << system[i]->GetTrait(SYS_POS) << ", Vote: " << system[i]->GetTrait(SYS_VOTE) << ", UID: " << system[i]->GetTrait(SYS_ID) << ", Message: " << -1  << ")" << ", " <<endl;
                
                system[i]->PrintState();
                cout << "Messages: ";
                for(auto m : message[i])
                    cout << m << ", ";
                cout << endl;
                cout << endl;

            }
        }
        // Print genomes!
        void PrintGenomes()
        {
            for(int i = 0; i < system.size(); ++i)
            {
                cout << i << ": " << endl;
                system[i]->PrintGenome();
                cout << endl;
                cout << "+++++++++++++++++++++++++++++++" << endl;
                cout <<endl;
            }
        }
    };

  public:
    Experiment(const HPConfig & config) : 
    DIM(config.DIM()), BEST_TIME(config.NUM_ITER()), DIR(config.DIR()), NUM_ITER(config.NUM_ITER()), 
    POP_SIZE(config.POP_SIZE()), NUM_GENS(config.NUM_GENS()), RNG_SEED(config.RNG_SEED()), 
    TOURN_SIZE(config.TOURN_SIZE()), SNAP_SHOT(config.SNAP_SHOT()), INST_MUT_RATE(config.INST_MUT_RATE()), 
    ARG_MUT_RATE(config.ARG_MUT_RATE()), DEL_MUT_RATE(config.DEL_MUT_RATE()), INS_MUT_RATE(config.INS_MUT_RATE()), 
    MAX_SIZE_G(config.MAX_SIZE_G()), MIN_SIZE_G(config.MIN_SIZE_G())
    {
        mRng = emp::NewPtr<emp::Random>(RNG_SEED);
        mWorld = emp::NewPtr<world_t>(*mRng, "System_World");
        mSys = emp::NewPtr<System>();
        BEST_TIME++;

        for(size_t i = 0; i < (DIM*DIM); ++i)
            mOrder.push_back(i);

        for(int i = 0; i < POP_SIZE; ++i)
            mPopIdx.push_back(i);
    }

    ~Experiment()
    {
        mRng.Delete();
        mWorld.Delete();
        for(auto p : mSys->system)
            p.Delete();

        mSys.Delete();
    }

    /* EXPERIMENT FUNCTIONS */

    // Initialized the Experiment!
    void Run();
    // Evaluate each organism
    pair<size_t, double> Evaluate(size_t gen);
    // Create troidal graph
    void Toroidal(int N, unordered_map<int, vector<int>> &adj);
    // Configure all requirements 
    void Config_All();
    // Configure instructions
    void Config_Inst();
    // Configure world
    void Config_Wrd();
    // Configure hardware
    void Config_HW();
    // Configure system
    void Config_SYS();
    // Get id inst
    static void Inst_GID(hardware_t & hw, const inst_t & inst);
    


    /* TESTING FUNCTIONS */

    // Test 1: Load a genome into a distributed system
    void test1();
    // Test 2: Check id scoring for system
    void test2();
    // Test 3: Execute an actual system!
    void test3();
    // Test 4: Evovle the organism!
    void test4();


    /* EXPERIMENT VARIABLES */

    // Training Set Adjacency List
    unordered_map<int, vector<int>> train_adj;
    // Random Number Generator Pointer
    rng_t mRng;
    // World Pointer
    wrd_t mWorld;
    // World index holder
    vector<int> mPopIdx;
    // Node order of executing
    vector<size_t> mOrder;
    // Instruction library pointer
    inst_lib_t inst_lib;
    // Dimensionality of grid
    size_t DIM;
    // Ponter of the system
    emp::Ptr<System> mSys;
    // Best overall time
    size_t BEST_TIME;
    // Best overall score
    double BEST_SCORE = 0.0;
    // Data Directory
    string DIR; 
    // Writer to csv
    ofstream CSV;
    ofstream fs;



    /* HARDWARE SETTING VARIABLES  - CONFIG*/

    // Location where the UID 
    size_t NUM_ITER; 
    // Population size
    size_t POP_SIZE;
    // Number of Generations
    size_t NUM_GENS;
    // Random Seed
    int RNG_SEED;
    // Tournament size
    size_t TOURN_SIZE;
    // Snapshot time
    size_t SNAP_SHOT;


    /* MUTATION SPECIFIC PARAMATERS - CONFIG */

    // Instruction mutation rate
    double INST_MUT_RATE;
    // Argument mutation rate
    double ARG_MUT_RATE;
    // Deletion instruction rate
    double DEL_MUT_RATE;
    // Insertion instruction rate
    double INS_MUT_RATE;
    // Maximum genome size
    size_t MAX_SIZE_G;
    // Minimum genome size
    size_t MIN_SIZE_G;
};

// Create NxN Toroidal Graph
void Experiment::Toroidal(int N, unordered_map<int, vector<int>> &adj)
{
    vector<vector<int>> temp(N);
    for(auto &v : temp)
    {
        v.resize(N);
    }

    //Fill vectors with label numbers
    int cnt = 0;
    for(int i = 0; i < N; ++i)
    {
        for(int j = 0; j < N; ++j)
        {
            temp[i][j] = cnt;
            ++cnt;
        }
    }

    //Create adjacency lists
    for(int i = 0; i < N; ++i)
    {
        for(int j = 0; j < N; ++j)
        {
            vector<int> t;
            // Get neighboors!
            t.push_back(temp[(N+i+1)%N][j]);
            t.push_back(temp[(N+i-1)%N][j]);
            t.push_back(temp[i][(N+j+1)%N]);
            t.push_back(temp[i][(N+j-1)%N]);   

            //Erase duplicates
            t.erase(unique(t.begin(), t.end()), t.end());

            adj[temp[i][j]] = t;           
        }
    }
}

// Evaluate each organism
pair<size_t, double> Experiment::Evaluate(size_t gen)
{
    size_t best_time = NUM_ITER;
    double best_score = 0.0;
    size_t best_org = 0;

    for(size_t org_id : mPopIdx)
    {
        // Get the genome
        auto & org = mWorld->GetOrg(org_id);
        // Set the genome in the system
        mSys->SetGenome(org.GetGenome());
        // Randomize the UID
        mSys->SetUID(mRng);
    
        for(size_t i = 0; i < NUM_ITER; ++i)
        {
            emp::Shuffle(*mRng, mOrder);

            // If consensus is reached
            if(mSys->Run(mOrder))
            {
                cout << "CONSENSUS!" << endl;
                hardware_t hw(org.GetGenome());
                hw.PrintGenome();
                if(i < best_time)
                    best_time = i;

                org.time = i; 
                org.score = mSys->mScore + (NUM_ITER - i);

                if(org.score > best_score)
                {
                    best_score = org.score;
                    best_org = org_id;
                }

                if(org.score > BEST_SCORE)
                    BEST_SCORE = org.score;

                if(org.time < BEST_TIME && org.time != -1)
                    BEST_TIME = org.time;

                break;
            }
            else
                org.score = mSys->mScore;

            if(org.score > best_score)
            {
                best_score = org.score;
                best_org = org_id;
            }
            if(org.score > BEST_SCORE)
                BEST_SCORE = org.score;

            if(org.time < BEST_TIME && org.time != -1)
                BEST_TIME = org.time;
        }
        mSys->Reset();
    }

    if(gen % SNAP_SHOT == 0)
    {
        hardware_t p(mWorld->GetOrg(best_org).GetGenome());
        cout << endl;
        p.PrintGenome();
        cout << endl;
    }

    return make_pair(best_time, best_score);
}

// Run the experiment!
void Experiment::Run()
{    
    // Configure all of the configs
    Experiment::Config_All();

    string dir = DIR + "/data.csv";
    CSV.open(dir);
    CSV << "GEN, BEST_TIME, BEST_SCORE" << endl;
    for(size_t i = 0; i < NUM_GENS; ++i)
    {
        cout << "Gen " << i << ": " << endl;
        auto best = Experiment::Evaluate(i);
        cout << "best_curr_time: " << best.first << ", best_curr_score: " << best.second \
        << ", best_all_time: " << BEST_TIME << ", best_all_score: " << BEST_SCORE << endl;

        CSV << i << ", " << BEST_TIME << ", " << BEST_SCORE << endl;

        // Do selection and mutations
        emp::TournamentSelect(*mWorld, TOURN_SIZE, mWorld->GetSize());
        mWorld->DoMutations();
    }    
}

// Configure all requirements 
void Experiment::Config_All()
{
    Experiment::Config_Inst();
    Experiment::Config_Wrd();
    Experiment::Config_HW();
    Experiment::Config_SYS();
}

// Configure world
void Experiment::Config_Wrd()
{
    mWorld->Reset();
    mWorld->SetPopStruct_Mixed(true);
    mWorld->SetFitFun([](Agent & a) {return a.score;});
    mWorld->SetMutFun([this](Agent & s, emp::Random & r)
    {
        hardware_t org(inst_lib);
        org.SetGenome(s.GetGenome());

        // Go through the whole genome
        for(size_t i = 0; i < org.GetSize(); ++i)
        {
            // Mutate the instruction?
            if(r.P(INST_MUT_RATE))
                org.RandomizeInst(i, r);

            // Mutate the instruction arguments?
            for(size_t j = 0; j <  emp::AvidaGP::base_t::INST_ARGS; ++j)
            {
                if(r.P(ARG_MUT_RATE))
                    org.genome.sequence[i].args[j] = r.GetUInt(org.CPU_SIZE);
            }

            // Insert an instruction?
            if(r.P(INS_MUT_RATE))
            {
                if(org.GetSize() < MAX_SIZE_G)
                {
                    org.genome.sequence.insert(org.genome.sequence.begin() + (int)i, emp::AvidaGP::genome_t::sequence_t::value_type());
                    org.RandomizeInst(i, r);
                }
            }

            // Delete an instruction?
            if(r.P(DEL_MUT_RATE))
            {
                if(org.GetSize() > MIN_SIZE_G)
                {
                    org.genome.sequence.erase(org.genome.sequence.begin() + (int)i);
                }
            }
        }
        s.SetGenome(org.GetGenome());
        return 0;
    });
}

// Configure hardware
void Experiment::Config_HW()
{
    for(size_t i = 0; i < POP_SIZE; ++i)
    {
        program_t p(inst_lib);
        Agent s(p);
        hardware_t org(s.GetGenome());
        int upper = mRng->GetUInt(MIN_SIZE_G, MAX_SIZE_G);
        // Add different amount and random instructions in a system genome
        for(size_t j = 0; j < upper; ++j)
        {
            auto inst = org.GetRandomInst(*mRng);
            org.PushInst(inst);
        }

        s.SetGenome(org.GetGenome());
        mWorld->Inject(s.GetGenome(), 1);
    }
}

// Configure hardware
void Experiment::Config_SYS()
{
    // Create toroidal graph
    Experiment::Toroidal(DIM, train_adj);
    // Load the adj list
    mSys->SetAdj(train_adj);
    mSys->SetSize(DIM);
}

// Configure instructions
void Experiment::Config_Inst()
{
    // Get default AvidaGP instructions
    inst_lib = emp::AvidaGP::inst_lib_t::DefaultInstLib();

    // Add instruction to send information
    inst_lib.AddInst("Broadcast", [this](hardware_t &hw, const inst_t &inst)
    {
        auto neigh = this->mSys->adj_list;

        for(int n : neigh[hw.GetTrait(SYS_POS)])
        {
            if(this->mSys->message[n].size() < MAX_MESS)
                this->mSys->message[n].push_back(hw.regs[inst.args[0]]);            
        }
    }, 
    1, "Local Memory Arg1 => All neighbors");

    // Add instruction to load hw id
    inst_lib.AddInst("GetID", Inst_GID, 1, "TRAIT[ID] => Local Memory Arg1");

    // Get the oldest message of queue and store it in registers
    inst_lib.AddInst("OldMessage", [this](hardware_t &hw, const inst_t &inst)
    {
        //auto messages = this->mSys->message[hw.GetTrait(SYS_POS)];
        if(this->mSys->message[hw.GetTrait(SYS_POS)].size() > 0)
        {
            hw.regs[inst.args[0]] = this->mSys->message[hw.GetTrait(SYS_POS)].front();
            this->mSys->message[hw.GetTrait(SYS_POS)].pop_front();
        }
        else
        {
            hw.regs[inst.args[0]] = -1;
        }
    }, 
    1, "MESSAGE.first => Local Memory Arg1");

    // Get the newest message of queue asd store it in registers
    inst_lib.AddInst("RecentMessage", [this](hardware_t &hw, const inst_t &inst)
    {
        //auto messages = this->mSys->message[hw.GetTrait(SYS_POS)];
        if(this->mSys->message[hw.GetTrait(SYS_POS)].size() > 0)
        {
            hw.regs[inst.args[0]] = this->mSys->message[hw.GetTrait(SYS_POS)].back();
            this->mSys->message[hw.GetTrait(SYS_POS)].pop_back();
        }
        else
        {
            hw.regs[inst.args[0]] = -1;
        }
    }, 
    1, "MESSAGE.last => Local Memory Arg1");

    // Store the vote 
    inst_lib.AddInst("SetVote", [](hardware_t &hw, const inst_t &inst)
    {
        hw.SetTrait(SYS_VOTE ,hw.regs[inst.args[0]]);
    }, 
    1, " Local Memory Arg 1 => VOTE");
    // Store the vote 
    inst_lib.AddInst("GetVote", [](hardware_t &hw, const inst_t &inst)
    {
        hw.regs[inst.args[0]] = hw.GetTrait(SYS_VOTE);
    }, 
    1, " Local Memory Arg 1 => VOTE");
}

// Get id inst
void Experiment::Inst_GID(hardware_t & hw, const inst_t & inst)
{
    hw.regs[inst.args[0]] = hw.GetTrait(SYS_ID);
}


/* TESTING FUNCTION */

//test 1: Load a genome into a distributed system
void Experiment::test1()
{
    // Create random program
    program_t p(inst_lib);
    hardware_t org(p);

    // Load random instructions
    for(int i = 0; i < 20; ++i)
    {
        auto inst = org.GetRandomInst(*mRng);
        org.PushInst(inst);
    }
    cout << "SEED: " << mRng->GetSeed() << endl;
    cout << "Genome: " << endl;
    org.PrintGenome();
    cout << endl;

    // load it into the system
    mSys->SetGenome(org.GetGenome());
    mSys->SetSize(DIM * DIM);

    auto system = mSys->system;
    cout << "SYSTEM N: " << mSys->N << endl;
    cout << endl;
    for(int i = 0; i < mSys->N; ++i)
    {
        cout << "***i: " << i << endl;
        system[i]->PrintGenome();
        cout << "---------------------------------\n" << endl;
    }


}

//test 2: Test system for scoring
void Experiment::test2()
{
    // Create random program
    program_t p(inst_lib);
    hardware_t org(p);

    // Load random instructions
    for(int i = 0; i < 20; ++i)
    {
        auto inst = org.GetRandomInst(*mRng);
        org.PushInst(inst);
    }
    cout << "SEED: " << mRng->GetSeed() << endl;

    // load it into the system
    mSys->SetSize(DIM * DIM);
    mSys->SetUID(mRng);

    vector<int> id;
    for(auto hw : mSys->system)
        id.push_back(hw->GetTrait(SYS_ID));
    
    mSys->PrintTraits();
    cout << "Score: " << mSys->GetScore() << endl;
    cout << endl;

    // Add legal votes
    mSys->system[0]->SetTrait(SYS_VOTE, id[0]);
    mSys->PrintTraits();
    cout << "Score: " << mSys->GetScore() << endl;
    cout << endl;

    mSys->system[0]->SetTrait(SYS_VOTE, id[0]);
    mSys->system[1]->SetTrait(SYS_VOTE, id[0]);
    mSys->PrintTraits();
    cout << "Score: " << mSys->GetScore() << endl;
    cout << endl;

    mSys->system[0]->SetTrait(SYS_VOTE, id[0]);
    mSys->system[1]->SetTrait(SYS_VOTE, id[0]);
    mSys->system[2]->SetTrait(SYS_VOTE, id[2]);
    mSys->PrintTraits();
    cout << "Score: " << mSys->GetScore() << endl;
    cout << endl;

    mSys->system[0]->SetTrait(SYS_VOTE, id[0]);
    mSys->system[1]->SetTrait(SYS_VOTE, id[0]);
    mSys->system[2]->SetTrait(SYS_VOTE, id[2]);
    mSys->system[3]->SetTrait(SYS_VOTE, id[2]);
    mSys->PrintTraits();
    cout << "Score: " << mSys->GetScore() << endl;
    cout << endl;

    mSys->system[0]->SetTrait(SYS_VOTE, id[0]);
    mSys->system[1]->SetTrait(SYS_VOTE, id[0]);
    mSys->system[2]->SetTrait(SYS_VOTE, id[0]);
    mSys->system[3]->SetTrait(SYS_VOTE, id[2]);
    mSys->PrintTraits();
    cout << "Score: " << mSys->GetScore() << endl;
    cout << endl;
}

//test 3: Test system for scoring
void Experiment::test3()
{
    int N = 4;
    Experiment::Toroidal(N, train_adj);
    // Load the adj list
    mSys->SetAdj(train_adj);

    // Create random program
    program_t p(inst_lib);
    hardware_t org(p);

    org.SetTrait(SYS_ID, 17);
    org.SetTrait(SYS_VOTE, -1);
    org.SetTrait(SYS_POS, 0);

    // Store ID in R[1]
    org.PushInst("GetID", 0);
    // Vote for ID in R[1]
    org.PushInst("SetVote", 0);
    // Send everyone your vote in R[1]
    org.PushInst("Broadcast", 0);
    // While loooooop
    org.PushInst("While", 1, 1);
    // Grab most recent message and put in R[2]
    org.PushInst("OldMessage", 1);
    // See if VOTE < RECENT MESSAGE
    org.PushInst("GetVote", 0);
    org.PushInst("TestLess", 0, 1, 3);
    org.PushInst("If", 3, 2);
    //Set vote to new higher vote
    org.PushInst("SetVote", 1);
    org.PushInst("Broadcast", 1);
    org.PushInst("If", 1);
    cout << endl;

    mSys->SetGenome(org.GetGenome());
    mSys->SetSize(N);
    mSys->SetUID(mRng);

    vector<size_t> order;
    for(int i = 0; i < N*N; ++i)
        order.push_back(i);

    for(int i = 0; i < 100; ++i)
    {
        cout << "I: " << i << endl;
        if(mSys->Run(order))
        {
            cout << "*******Reached Consensus*********" << endl;
            break;
        }
        else
        {
            cout << "Max: " << mSys->mScore << endl;
        }
    }
}

//test 3: Evolve an organism!
void Experiment::test4()
{
    int N = 3;
    // Create random program
    program_t p(inst_lib);
    hardware_t org(p);

    org.PushRandom(*mRng, (size_t) mRng->GetUInt(10,20));
    org.PrintGenome();

    mWorld->Inject(org.GetGenome());
    mWorld->DoMutations();
    org.SetGenome(mWorld->GetGenomeAt(0));
    cout << endl;
    cout << "__________________________________" << endl;
    cout << endl;
    org.PrintGenome();

    

    mSys->SetGenome(org.GetGenome());
    mSys->SetSize(N);
    mSys->SetUID(mRng);
}

#endif