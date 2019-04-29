#ifndef HP_EXPERIMENT_H
#define HP_EXPERIMENT_H

#include <iostream>

#include "hp_config.h"
#include "Graph.h"
#include "../../Empirical/source/tools/Random.h"
#include "../../Empirical/source/tools/random_utils.h"
#include "../../Empirical/source/hardware/EventDrivenGP.h"
#include "../../Empirical/source/Evolve/World.h"
#include "../../Alex/signalgp-adventures/adventures/utility_belt/source/utilities.h"

/* CONSTEXPR FOR HARDWARE */

//Number of bits in the tag per hardware
constexpr size_t TAG_WIDTH = 16;
constexpr double VALUE = .0000001;

const std::string NOP_PATH = "genome1.txt";
const std::string BASIC_PATH = "genome.txt";
const std::string ADVANCE_PATH = "genome2.txt";

/* NEW TYPE DECLARATIONS FOR HARDWARE*/

//Type for a actural hardware
using hardware_t = emp::EventDrivenGP_AW<TAG_WIDTH>;
//Type for the hardware genome
using program_t = hardware_t::Program;
//Type for the hardwares state 
using state_t = hardware_t::State;
//Instruction object for hardware
using inst_t = hardware_t::inst_t;
//Instruction library for hardware
using inst_lib_t = hardware_t::inst_lib_t;
//Type for Events for each hardware
using event_t = hardware_t::event_t;
//Event Libarary for hardware
using event_lib_t = hardware_t::event_lib_t;
//Function type for the hardware
using function_t = hardware_t::Function;
//Memory type for hardware
using memory_t = hardware_t::memory_t;
//Mutator
using mutant_t = toolbelt::SignalGPMutator<hardware_t>;

/* NEW TYPE DECLARATIONS FOR SIMPLICITY*/

using coor_t = std::pair<size_t, size_t>;


class Experiment
{
  struct Agent;

  //World of Agents
  using world_t = emp::World<Agent>;  

  struct Agent
  {
    //Agents score
    double mScore = 0;

    program_t mGenome;

    Agent(const program_t & p) : mGenome(p) {;}

    program_t & GetGenome() {return mGenome;}
  };

  public:
    Experiment(const HPConfig & config) :
    POP_SIZE(config.POP_SIZE()), NUM_GENS(config.NUM_GENS()),
    RNG_SEED(config.RNG_SEED()), EVAL_SIZE(config.EVAL_SIZE()),
    TOURN_SIZE(config.TOURN_SIZE()), GRA_DIM(config.GRA_DIM()),
    GRA_TYPE(config.GRA_TYPE()), SNAP_SHOT(config.SNAP_SHOT()),
    NUM_ITER(config.NUM_ITER()), MIN_FUN_CNT(config.MIN_FUN_CNT()),
    MAX_FUN_CNT(config.MAX_FUN_CNT()), MIN_FUN_LEN(config.MIN_FUN_LEN()), 
    MAX_FUN_LEN(config.MAX_FUN_LEN()), MAX_TOT_LEN(config.MAX_TOT_LEN())    
    {
      mRng = emp::NewPtr<emp::Random>(RNG_SEED);
      inst_lib = emp::NewPtr<inst_lib_t>();
      event_lib = emp::NewPtr<event_lib_t>();
      mGraph = emp::NewPtr<Graph>(config, mRng);
      mWorld = emp::NewPtr<world_t>(*mRng, "World");
      mMutant = emp::NewPtr<mutant_t>(MIN_FUN_CNT, MAX_FUN_CNT, MIN_FUN_LEN, MAX_FUN_LEN, MAX_TOT_LEN);
      THEORY_MAX = NUM_ITER * GRA_DIM * GRA_DIM + 1;
      SetVote(config.VOTE()); SetUID(config.UID());
      SetPOSX(config.POSX()); SetPOSY(config.POSY());
    }

    ~Experiment()
    {
      mGraph.Delete();
      mWorld.Delete();
      mRng.Delete();
      mMutant.Delete();
      inst_lib.Delete();
      event_lib.Delete();
    }

    /* FUNCTIONS DEDICATED TO THE EXPERIMENT */

    //Run the experiment
    void Run();
    
    //Confiugre all the neccesary things
    void Config_All();

    //Evalute each agent for 
    size_t Evaluation_step();

    //Selection
    void Selection_step();

    //Update
    void Update_step();

    //Return an genome full of nops
    program_t Genome_NOP();

    //Return basic program
    program_t Genome_BASIC();

    //Return advance program
    program_t Genome_ADVANCE();


    /* FUNCTIONS DEDICATED TO THE CONFIGURATIONS, INSTRUCTIONS, EVENTS*/

    //Will make the instruction library
    void Config_Inst();

    //Will set the output memory of a hardware to all of its neighboors
    static void Inst_BroadcastMail(hardware_t & hw, const inst_t & inst);

    //Will send the vote to all of a hardwares neighboors
    static void Inst_BroadcastVote(hardware_t & hw, const inst_t & inst);

    //Will load the UID of a hardware into its working buffer
    static void Inst_GetUID(hardware_t & hw, const inst_t & inst);
    
    //Will get the vote of a hardware
    static void Inst_GetVote(hardware_t & hw, const inst_t & inst);

    //Will set the vote of the hardware
    static void Inst_SetVote(hardware_t & hw, const inst_t & inst);

    //Will make the event library
    void Config_Events();

    //Will actually do the event
    void Dispatch_Broadcast(hardware_t & hw, const event_t & e);

    //Will do the event to send vote out
    void Dispatch_BroadcastVote(hardware_t & hw, const event_t & e);

    //Will spawn a core for the event
    static void Handle_Broadcast(hardware_t & hw, const event_t & e);

    //Will create the hardware and add it to the world
    void Config_HW(program_t p);

    //Will configure the world
    void Config_World();


    /* FUNCTIONS DEDICATED TO SETTERS */
    static void SetVote(size_t x) {VOTE = x;}
    static void SetUID(size_t x) {UID = x;}
    static void SetPOSX(size_t x) {POSX = x;}
    static void SetPOSY(size_t x) {POSY = x;}


    /* FUNCTIONS DEDICATED TO TEST GRAPH */
    
    //Will test if the graph is actually GRA_DIM X GRA_DIM
    void GraphTest1();

    //Will test the adjacency list
    void GraphTest2();

    //Will test setting up initial traits
    void GraphTest3();

    //Will test to see if the FinalVotes vector is correct
    void GraphTest4();    

    //Will test to se if scoring functions are correct
    void GraphTest5();

    //Will test loading GENOME
    void GraphTest6();

    //Will test step by step process of the hardware
    void GraphTest7();

  private:

    /* EXPERIMENT SPECIFIC PARAMATERS */

    //Size of the population we are evaluating.
    size_t POP_SIZE;
    //Number of genearations evaluating during experiments
    size_t NUM_GENS;
    //Random Number Seed Number
    size_t RNG_SEED;
    //Number of bad guys a good guy will face per Run
    size_t EVAL_SIZE;
    //Number of organisms competing during tournament do_selection_step
    size_t TOURN_SIZE;
    //Dimension of the graph
    size_t GRA_DIM;
    //Pointer for random number generator
    emp::Ptr<emp::Random> mRng;
    //Pointer for Graph
    emp::Ptr<Graph> mGraph;
    //Graph tyep
    size_t GRA_TYPE;
    //vector to hold positions of hardware
    std::vector<size_t> mHwPos;
    //World to hold organisms
    emp::Ptr<world_t> mWorld;
    //Mutator
    emp::Ptr<mutant_t> mMutant;
    //Snapshot index
    size_t SNAP_SHOT;
    //All the iterations
    size_t NUM_ITER;
    //Count of how many time broadcast vote is called
    double mCount = 0;
    //Theoretical max
    size_t THEORY_MAX;

    /* HARDWARE SPECIFIC PARAMATERS */

    //Position of UID within hw trait vector
    static size_t UID;
    //Position of VOTE within hw trait vector
    static size_t VOTE;
    //Position of X coordinate within hw trait vector
    static size_t POSX;
    //Position of Y coordinate within hw trait vector
    static size_t POSY;
    //Instruction Library for hardware
    emp::Ptr<inst_lib_t> inst_lib;
    //Event Library for hardware
    emp::Ptr<event_lib_t> event_lib;

    /* MUTATION SPECIFIC PARAMATERS */
    size_t MIN_FUN_CNT;
    size_t MAX_FUN_CNT;
    size_t MIN_FUN_LEN;
    size_t MAX_FUN_LEN;
    size_t MAX_TOT_LEN;
};

size_t Experiment::UID = 0;
size_t Experiment::VOTE = 0;
size_t Experiment::POSX = 0;
size_t Experiment::POSY = 0;

/* FUNCTIONS DEDICATED TO THE EXPERIMENT */

//Run the experiment
void Experiment::Run()
{
  std::cout << "SETTING UP CONFIGS!" << std::endl;
  Experiment::Config_All();
  std::cout << "CONFIGS SET!" << std::endl;

  for(size_t i = 0; i < NUM_GENS; ++i)
  {
    std::cout << "GEN: " << i;
    size_t best_org = Experiment::Evaluation_step();
    if((i%SNAP_SHOT) == 0)
    {
      std::cout << std::endl;
      Agent & agent = mWorld->GetOrg(best_org);
      agent.mGenome.PrintProgramFull();
      std::cout << std::endl;
    }
    Experiment::Selection_step();
    Experiment::Update_step();
  }
}

//Confiugre all the neccesary things
void Experiment::Config_All()
{
  Experiment::Config_Inst();
  Experiment::Config_Events();
  Experiment::Config_World();

  program_t pro = Experiment::Genome_NOP();
  //program_t pro = Experiment::Genome_ADVANCE();  
  //program_t pro = Experiment::Genome_BASIC();

  Experiment::Config_HW(pro);


  std::cout << "CREATING THE GRAPH!" << std::endl;
  mGraph->CreateGraph(GRA_DIM, GRA_TYPE, inst_lib, event_lib);
  mGraph->ConfigureTraits();
  mGraph->CreateAdjList(GRA_TYPE, GRA_DIM);
  std::cout << "GRAPH CREATED!" << std::endl;
}

//Evalute each agent for 
size_t Experiment::Evaluation_step()
{
  double best = -999;
  size_t best_org = 0;
  for(size_t i = 0; i < POP_SIZE; ++i)
  {
    mCount = 0;
    mGraph->Reset();
    Agent & agent = mWorld->GetOrg(i);
    program_t & pro = agent.GetGenome();
    mGraph->SetGenome(pro);
    double score = mGraph->RunGraph();

    if(mCount > 10)
    {
      score += 1;
    }

    else
    {
      score += (mCount * VALUE);
    }    
    agent.mScore = score;

    if(score > best)
    {
      best_org = i;
      best = score;
      //mGraph->PrintTraits();
    }    
  }

  std::cout << " Best Score: " << best  << " THEORY_MAX: " << THEORY_MAX << " SUCESS%: " << (best / THEORY_MAX) << std::endl;
  return best_org;
}

//Selection
void Experiment::Selection_step()
{
  emp::TournamentSelect(*mWorld, TOURN_SIZE, POP_SIZE);
}

//Update
void Experiment::Update_step()
{
  mWorld->Update();
  mWorld->DoMutations();
}

//Return an genome full of nops
program_t Experiment::Genome_NOP()
{
  program_t ancestor_prog(inst_lib);
  std::ifstream ancestor_fstream(NOP_PATH);

  if (!ancestor_fstream.is_open())
  {
    std::cout << "Failed to open ancestor program file(" << NOP_PATH << "). Exiting..." << std::endl;
    exit(-1);
  }

  ancestor_prog.Load(ancestor_fstream);

  return ancestor_prog;
}

//Return basic program
program_t Experiment::Genome_BASIC()
{
  program_t ancestor_prog(inst_lib);
  std::ifstream ancestor_fstream(BASIC_PATH);

  if (!ancestor_fstream.is_open())
  {
    std::cout << "Failed to open ancestor program file(" << BASIC_PATH << "). Exiting..." << std::endl;
    exit(-1);
  }

  ancestor_prog.Load(ancestor_fstream);

  return ancestor_prog;
}

//Return advance program
program_t Experiment::Genome_ADVANCE()
{
  program_t ancestor_prog(inst_lib);
  std::ifstream ancestor_fstream(ADVANCE_PATH);

  if (!ancestor_fstream.is_open())
  {
    std::cout << "Failed to open ancestor program file(" << ADVANCE_PATH << "). Exiting..." << std::endl;
    exit(-1);
  }

  ancestor_prog.Load(ancestor_fstream);

  return ancestor_prog;
}


/* FUNCTIONS DEDICATED TO THE CONFIGURATIONS, INSTRUCTIONS, EVENTS */

//Will make the instruction library
void Experiment::Config_Inst()
{

  //Experiment specific instructions
  inst_lib->AddInst("GetUID", Inst_GetUID, 1, "UID => Local Memory Arg1");
  inst_lib->AddInst("Broadcast", Inst_BroadcastMail, 1, "Output Memory => hw.neighboors");
  //inst_lib->AddInst("BroadcastVote", Inst_BroadcastVote, 1, "Output Memory[Arg1] = Vote => hw.neighboors");
  inst_lib->AddInst("GetVote", Inst_GetVote, 1, "Vote => Local Memory Arg1");
  inst_lib->AddInst("SetVote", Inst_SetVote, 1, "Local Memory Arg1 => Hw.Trait[VOTE]");

  // - Setup the instruction set. -
  // Standard instructions:

  inst_lib->AddInst("Inc", hardware_t::Inst_Inc, 1, "Increment value in local memory Arg1");
  inst_lib->AddInst("Dec", hardware_t::Inst_Dec, 1, "Decrement value in local memory Arg1");
  inst_lib->AddInst("Not", hardware_t::Inst_Not, 1, "Logically toggle value in local memory Arg1");
  inst_lib->AddInst("Add", hardware_t::Inst_Add, 3, "Local memory: Arg3 = Arg1 + Arg2");
  inst_lib->AddInst("Sub", hardware_t::Inst_Sub, 3, "Local memory: Arg3 = Arg1 - Arg2");
  inst_lib->AddInst("Mult", hardware_t::Inst_Mult, 3, "Local memory: Arg3 = Arg1 * Arg2");
  inst_lib->AddInst("Div", hardware_t::Inst_Div, 3, "Local memory: Arg3 = Arg1 / Arg2");
  inst_lib->AddInst("Mod", hardware_t::Inst_Mod, 3, "Local memory: Arg3 = Arg1 % Arg2");
  inst_lib->AddInst("TestEqu", hardware_t::Inst_TestEqu, 3, "Local memory: Arg3 = (Arg1 == Arg2)");
  inst_lib->AddInst("TestNEqu", hardware_t::Inst_TestNEqu, 3, "Local memory: Arg3 = (Arg1 != Arg2)");
  inst_lib->AddInst("TestLess", hardware_t::Inst_TestLess, 3, "Local memory: Arg3 = (Arg1 < Arg2)");
  inst_lib->AddInst("If", hardware_t::Inst_If, 1, "Local memory: If Arg1 != 0, proceed; else, skip block.", emp::ScopeType::BASIC, 0, {"block_def"});
  inst_lib->AddInst("While", hardware_t::Inst_While, 1, "Local memory: If Arg1 != 0, loop; else, skip block.", emp::ScopeType::BASIC, 0, {"block_def"});
  inst_lib->AddInst("Countdown", hardware_t::Inst_Countdown, 1, "Local memory: Countdown Arg1 to zero.", emp::ScopeType::BASIC, 0, {"block_def"});
  inst_lib->AddInst("Close", hardware_t::Inst_Close, 0, "Close current block if there is a block to close.", emp::ScopeType::BASIC, 0, {"block_close"});
  inst_lib->AddInst("Break", hardware_t::Inst_Break, 0, "Break out of current block.");
  inst_lib->AddInst("Call", hardware_t::Inst_Call, 0, "Call function that best matches call affinity.", emp::ScopeType::BASIC, 0, {"affinity"});
  inst_lib->AddInst("Return", hardware_t::Inst_Return, 0, "Return from current function if possible.");
  inst_lib->AddInst("SetMem", hardware_t::Inst_SetMem, 2, "Local memory: Arg1 = numerical value of Arg2");
  inst_lib->AddInst("CopyMem", hardware_t::Inst_CopyMem, 2, "Local memory: Arg1 = Arg2");
  inst_lib->AddInst("SwapMem", hardware_t::Inst_SwapMem, 2, "Local memory: Swap values of Arg1 and Arg2.");
  inst_lib->AddInst("Input", hardware_t::Inst_Input, 2, "Input memory Arg1 => Local memory Arg2.");
  inst_lib->AddInst("Output", hardware_t::Inst_Output, 2, "Local memory Arg1 => Output memory Arg2.");
  inst_lib->AddInst("Commit", hardware_t::Inst_Commit, 2, "Local memory Arg1 => Shared memory Arg2.");
  inst_lib->AddInst("Pull", hardware_t::Inst_Pull, 2, "Shared memory Arg1 => Local memory Arg2.");
  inst_lib->AddInst("Nop", hardware_t::Inst_Nop, 0, "No operation.");

  
}

//Will set the output memory of a hardware to all of its neighboors
void Experiment::Inst_BroadcastMail(hardware_t & hw, const inst_t & inst)
{
  state_t & state = hw.GetCurState();
  hw.TriggerEvent("BroadcastMail", inst.affinity, state.output_mem);
}

//Will set the output memory of a hardware to all of its neighboors
void Experiment::Inst_BroadcastVote(hardware_t & hw, const inst_t & inst)
{
  double vote = hw.GetTrait(VOTE);
  state_t & state = hw.GetCurState();
  state.SetOutput(inst.args[0], vote);

  hw.TriggerEvent("BroadcastVote", inst.affinity, state.output_mem);
}

//Will load the UID of a hardware into its working buffer
void Experiment::Inst_GetUID(hardware_t & hw, const inst_t & inst)
{
  state_t & state = hw.GetCurState();
  state.SetLocal(inst.args[0], hw.GetTrait(UID));
}

//Will get the vote of a hardware
void Experiment::Inst_GetVote(hardware_t & hw, const inst_t & inst)
{
  state_t & state = hw.GetCurState();
  state.SetLocal(inst.args[0], hw.GetTrait(VOTE));
}

//Will set the vote of the hardware
void Experiment::Inst_SetVote(hardware_t & hw, const inst_t & inst)
{
  state_t & state = hw.GetCurState();
  double vote = state.GetLocal(inst.args[0]);
  hw.SetTrait(VOTE, vote);
}

//Will make the event library
void Experiment::Config_Events()
{
  event_lib->AddEvent("BroadcastMail", Experiment::Handle_Broadcast, "Send output memory to all neighbors.");
  event_lib->RegisterDispatchFun("BroadcastMail", [this](hardware_t & hw, const event_t & e)
  {
    // state_t & state = hw.GetCurState();
    // auto output = state.output_mem;
    // auto graph = this->mGraph;

    // for(auto x : output)
    // {
    //   if(graph->Find(x.second))
    //   {
    //     this->mCount += 1;
    //     break;
    //   }
    // }

    this->Dispatch_Broadcast(hw, e);
  });


  event_lib->AddEvent("BroadcastVote", Experiment::Handle_Broadcast, "Send output memory to all neighbors.");
  event_lib->RegisterDispatchFun("BroadcastVote", [this](hardware_t & hw, const event_t & e)
  {
    auto graph = this->mGraph;
    double vote = hw.GetTrait(VOTE);

    if(graph->Find(vote))
    {
      this->mCount+= 1;
    }


    this->Dispatch_Broadcast(hw, e);
  });
}

//Will spawn a core for the event
void Experiment::Handle_Broadcast(hardware_t & hw, const event_t & e)
{
  hw.SpawnCore(e.affinity, hw.GetMinBindThresh(), e.msg);
}

//Will actually do the event
void Experiment::Dispatch_Broadcast(hardware_t & hw, const event_t & e)
{
  coor_t p = std::make_pair(hw.GetTrait(POSX), hw.GetTrait(POSY));
  auto team = mGraph->GetNodeNeig(p.first, p.second);

  for(auto pair : team)
  {
    auto node = mGraph->GetNode(pair.first, pair.second);
    node->mHW->QueueEvent(e);
  }
}

//Will create the hardware
void Experiment::Config_HW(program_t p)
{
  mWorld->Inject(p, POP_SIZE);
}

//Will configure the world
void Experiment::Config_World()
{
  mWorld->Reset();
  mWorld->SetPopStruct_Mixed(true);
  mWorld->SetFitFun([this](Agent & agent) {return agent.mScore;});
  mWorld->SetMutFun([this](Agent & agent, emp::Random & rnd)
  {
    program_t & p = agent.mGenome;
    return this->mMutant->ApplyMutations(p, rnd);
  });
}


/* FUNCTIONS DEDICATED TO TEST GRAPH */

//Will test if the graph is actually GRA_DIM X GRA_DIM
//Pass
void Experiment::GraphTest1()
{
  mGraph->CreateGraph(GRA_DIM, GRA_TYPE, inst_lib, event_lib);

  std::cout << "The graph should be " << GRA_DIM << "x" << GRA_DIM << std::endl;
  auto g = mGraph->GetGraph();
  if(g.size() != GRA_DIM)
    std::cout << "WRONG" << std::endl;

  for(auto col : g)
  {
    if(col.size() != GRA_DIM)
      std::cout << "WRONG" << std::endl;
  }
  std::cout << "Done with GraphTest1!" << std::endl;
}

//Will test the adjacency list
//Tested by hand.
//Pass
void Experiment::GraphTest2()
{
  mGraph->CreateGraph(GRA_DIM, GRA_TYPE, inst_lib, event_lib);
  mGraph->CreateAdjList(GRA_TYPE, GRA_DIM);

  auto g = mGraph->GetGraph();

  for(size_t i = 0; i < g.size(); ++i)
  {
    for(size_t j = 0; j < g[i].size(); ++j)
    {
      std::cout << "(" << i << "," << j << "): ";
      
      for(auto p : g[i][j]->mFriends)
      {
        std::cout << "(" << p.first << "," << p.second << "), ";
      }
      std::cout << std::endl;
    }
  }
}

//Will test to see if we can randomly assign a UID
//Tested by comparison of mRanNum Vector comparison
//Pass
void Experiment::GraphTest3()
{
  mGraph->CreateGraph(GRA_DIM, GRA_TYPE, inst_lib, event_lib);

  std::cout << "RUN 1" << std::endl;
  mGraph->ConfigureTraits();
  mGraph->PrintRanUID();
  mGraph->PrintTraits();

  std::cout << "RUN 2" << std::endl;
  mGraph->ConfigureTraits();
  mGraph->PrintRanUID();
  mGraph->PrintTraits();
}

//Will test to see if FinalVotes vector is correct
//Pass
void Experiment::GraphTest4()
{
  //Creates graph
  mGraph->CreateGraph(GRA_DIM, GRA_TYPE, inst_lib, event_lib);

  std::cout << "RUN1-SCORE=dim*dim\n" << std::endl;
  //Sets up traits
  mGraph->ConfigureTraits();
  auto uid = mGraph->GetRandNums();

  //Set all votes to first thing in the random number vec
  mGraph->SetVotes(uid[0]);

  //Create the final votes map score keepr
  mGraph->MakeFinalVotes();
  auto votes = mGraph->GetFinVotes();
  mGraph->PrintRanUID();
  mGraph->PrintFinVotes();
  mGraph->PrintTraits();

  std::cout << "RUN2-SCORE=dim*dim-1\n" << std::endl;

  mGraph->ConfigureTraits();
  auto uid2 = mGraph->GetRandNums();
  mGraph->SetVotes(uid2[0]);
  mGraph->SetVote(0,0,0);
  mGraph->MakeFinalVotes();
  auto votes2 = mGraph->GetFinVotes();
  mGraph->PrintRanUID();
  mGraph->PrintFinVotes(); 
  mGraph->PrintTraits(); 
}

//Will test to se if scoring functions are correct
void Experiment::GraphTest5()
{
  mGraph->CreateGraph(GRA_DIM, GRA_TYPE, inst_lib, event_lib);
  mGraph->ConfigureTraits();
  auto uid = mGraph->GetRandNums();

  std::cout << "RUN1: All hw should be in consnesus\n" << std::endl;
  mGraph->PrintRanUID();
  mGraph->SetVotes(uid[0]);
  mGraph->MakeFinalVotes();
  mGraph->PrintVotes();
  mGraph->PrintFinVotes();
  std::cout << "Consensus: " << (GRA_DIM*GRA_DIM) << " == " << mGraph->Consensus() << std::endl;
  std::cout << "Legal Votes: " << (GRA_DIM*GRA_DIM) << " == " << mGraph->LegalVotes() << std::endl;
  std::cout << "Largest Legal Votes: " << (GRA_DIM*GRA_DIM) << " == " << mGraph->LargestLegalVotes() << std::endl;
  std::cout << std::endl;
  
  mGraph->ConfigureTraits();
  uid = mGraph->GetRandNums();
  std::cout << "RUN2: All hw should be in consnesus except one ILLEGAL VOTE\n" << std::endl;
  mGraph->PrintRanUID();
  mGraph->SetVotes(uid[0]);
  mGraph->SetVote(0,0,0);
  mGraph->MakeFinalVotes();
  mGraph->PrintVotes();
  mGraph->PrintFinVotes();
  std::cout << "Consensus: " << 0 << " == " << mGraph->Consensus() << std::endl;
  std::cout << "Legal Votes: " << (GRA_DIM*GRA_DIM) - 1 << " == " << mGraph->LegalVotes() << std::endl;
  std::cout << "Largest Legal Votes: " << (GRA_DIM*GRA_DIM) -1 << " == " << mGraph->LargestLegalVotes() << std::endl;
  std::cout << std::endl;

  mGraph->ConfigureTraits();
  uid = mGraph->GetRandNums();
  std::cout << "RUN3: All hw should be in consnesus except one\n" << std::endl;
  mGraph->PrintRanUID();
  mGraph->SetVotes(uid[0]);
  mGraph->SetVote(0,0,uid[1]);
  mGraph->MakeFinalVotes();
  mGraph->PrintVotes();
  mGraph->PrintFinVotes();
  std::cout << "Consensus: " << 0 << " == " << mGraph->Consensus() << std::endl;
  std::cout << "Legal Votes: " << (GRA_DIM*GRA_DIM) << " == " << mGraph->LegalVotes() << std::endl;
  std::cout << "Largest Legal Votes: " << (GRA_DIM*GRA_DIM) -1 << " == " << mGraph->LargestLegalVotes() << std::endl;
  std::cout << std::endl;

  mGraph->ConfigureTraits();
  uid = mGraph->GetRandNums();
  std::cout << "RUN4: All hw should be in consnesus except one and one ILLEGAL VOTE\n" << std::endl;
  mGraph->PrintRanUID();
  mGraph->SetVotes(uid[0]);
  mGraph->SetVote(0,0,uid[1]);
  mGraph->SetVote(1,1,0);
  mGraph->MakeFinalVotes();
  mGraph->PrintVotes();
  mGraph->PrintFinVotes();
  std::cout << "Consensus: " << 0 << " == " << mGraph->Consensus() << std::endl;
  std::cout << "Legal Votes: " << (GRA_DIM*GRA_DIM) -1 << " == " << mGraph->LegalVotes() << std::endl;
  std::cout << "Largest Legal Votes: " << (GRA_DIM*GRA_DIM) - 2 << " == " << mGraph->LargestLegalVotes() << std::endl;
  std::cout << std::endl;
}

//Will test loading GENOME
void Experiment::GraphTest6()
{
  Experiment::Config_Inst();
  Experiment::Config_Events();
  Experiment::Config_World();

  mGraph->CreateGraph(GRA_DIM, GRA_TYPE, inst_lib, event_lib);
  mGraph->ConfigureTraits();
  mGraph->CreateAdjList();

  program_t p = Experiment::Genome_ADVANCE();
  p.PrintProgramFull();
  mGraph->SetGenome(p);
  std::cout << "Advance Score: " << mGraph->RunGraph() << std::endl;
  std::cout << std::endl;

  std::cout << "Before Reset" << std::endl;
  std::cout << "mRandomNums.size(): " << mGraph->GetRandNums().size() << std::endl;
  std::cout << "mFinalVotes.size(): " << mGraph->GetFinVotes().size() << std::endl;
  mGraph->PrintRanUID();
  mGraph->PrintFinVotes();
  mGraph->PrintGenomes();
  std::cout << std::endl;

  mGraph->Reset();
  std::cout << "After Reset" << std::endl;
  std::cout << "mRandomNums.size(): " << mGraph->GetRandNums().size() << std::endl;
  std::cout << "mFinalVotes.size(): " << mGraph->GetFinVotes().size() << std::endl;
  mGraph->PrintRanUID();
  mGraph->PrintFinVotes();
  
  std::cout << std::endl;

  program_t p1 = Experiment::Genome_NOP();
  p1.PrintProgramFull();
  mGraph->SetGenome(p1);
  mGraph->PrintGenomes();
  std::cout << "NOP Score: " << mGraph->RunGraph() << std::endl;
}

//Will test step by step process of the hardware
void Experiment::GraphTest7()
{
  Experiment::Config_Inst();
  Experiment::Config_Events();
  Experiment::Config_World();

  mGraph->CreateGraph(GRA_DIM, GRA_TYPE, inst_lib, event_lib);  
  mGraph->ConfigureTraits();
  mGraph->CreateAdjList(GRA_TYPE, GRA_DIM);
  mGraph->Reset();

  program_t p = Experiment::Genome_BASIC();
  p.PrintProgramFull();
  mGraph->SetGenome(p);

  mGraph->PrintFriends();

  //std::cout << "Advance Score: " << mGraph->RunGraph(15) << std::endl;
  std::cout << std::endl;
}

#endif