#ifndef HP_GRAPH_H
#define HP_GRAPH_H

#include <iostream>
#include <set>

#include "hp_config.h"
#include "../../Empirical/source/tools/Random.h"
#include "../../Empirical/source/tools/random_utils.h"
#include "../../Empirical/source/hardware/EventDrivenGP.h"
#include "../../Empirical/source/base/assert.h"

struct Node;

/* CONSTEXPR FOR HARDWARE */

//Number of bits in the tag per hardware
constexpr size_t TAG_WIDTH_ = 16; 

/* NEW TYPE DECLARATIONS FOR HARDWARE*/

//Type for a actural hardware
using hardware_t = emp::EventDrivenGP_AW<TAG_WIDTH_>;
//Type for the hardware genome
using program_t = hardware_t::Program;
//Instruction object for hardware
using ins_t = hardware_t::inst_t;
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

/* NEW TYPE DECLARATIONS FOR SIMPLICITY*/
using graph_t = std::vector<std::vector<emp::Ptr< Node >>>;
using coor_t = std::pair<size_t, size_t>;
using map_t = std::unordered_map<double, double>;
using nodes_t = std::vector<emp::Ptr<Node>>;
using randnum_t = std::vector<size_t>;

/* STRUCT USED TO IMITATE A NODE WITHIN THE SYSTEM */
struct Node 
{
  //Hardware within the node
  emp::Ptr<hardware_t> mHW;
  //All of the nodes that this node has access too
  std::set<coor_t> mFriends;

  Node(emp::Ptr<inst_lib_t> ilib, emp::Ptr<event_lib_t> elib, emp::Ptr<emp::Random> rng)
  {
    mHW = emp::NewPtr<hardware_t>(ilib, elib, rng);
  }
};

class Graph
{
  public:
    Graph(const HPConfig & config, emp::Ptr<emp::Random> rng = nullptr) :
    GRA_DIM(config.GRA_DIM()), NUM_ITER(config.NUM_ITER()), 
    NUM_FRI(config.NUM_FRI()), NUM_ENE(config.NUM_ENE()),
    mRng(rng), RNG_SEED(config.RNG_SEED()), MIN_BIN_THSH(config.MIN_BIN_THSH()),
    UID(config.UID()), VOTE(config.VOTE()), POSX(config.POSX()), 
    POSY(config.POSY()), MAX_BND(config.MAX_BND()), MIN_BND(config.MIN_BND()),
    MAX_CORES(config.MAX_CORES())
    {;}

    //Delete all pointers in the class
    ~Graph()
    {
      for(size_t i = 0; i < mGraph.size(); ++i)
      {
        for(size_t j = 0; j < mGraph[i].size(); ++j)
        {
          mGraph[i][j]->mHW.Delete();
          mGraph[i][j].Delete();
        }
      }
      
      mNodes.clear();
    }

    /* FUNCTIONS DEDICATED TO THE STRUCTURE */

    //Give the graph NUM_ITER single processes to figure it out
    double RunGraph(size_t iter = -1);

    //Function will create a general graph structure and set the x and y position per hardware
    void CreateGraph(size_t dim = 2, size_t type = 0, emp::Ptr<inst_lib_t> ilib = nullptr, emp::Ptr<event_lib_t> elib  = nullptr);

    //Will create adjacency list for each node
    void CreateAdjList(size_t type = 0, size_t dim = 2);

    /* FUNCTIONS DEDICATED TO RUNNING EXPERIMENT */

    //Will see if element exits in mRanNums vector
    bool Find(size_t me);

    //Will configure traits for hardware UID
    void ConfigureTraits();

    //Will create a dictinary with legal votes that the hardware has come up with
    void MakeFinalVotes();

    //Will return all legal votes
    double LegalVotes();

    //Will return the largest vote
    double LargestLegalVotes();

    //Will return if the system is in consensus or not
    double Consensus();

    //Will reset the graph to rerun with different program
    void Reset();
    

    /* FUNCTIONS DEDICATED TO PRINTING OUT CRAP */

    //Print the traits
    void PrintTraits();

    //Print the random UID
    void PrintRanUID();

    //Print the FinalVotes vector
    void PrintFinVotes();

    //Print Vote Trait
    void PrintVotes();

    //Print each nodes genome
    void PrintGenomes();

    //Print the mSchedule
    void PrintSchedule();

    //Print all the friends
    void PrintFriends();


    /* FUNCTIONS DEDICATED TO BE GETTERS */

    //Reuturn the graph in matrix form
    graph_t GetGraph() const {return mGraph;}

    //Return the graph in list form
    nodes_t GetNodes() const {return mNodes;} 

    //Return Final Votes at time i called
    map_t GetFinVotes() const {return mFinalVotes;}

    //Return Random Number UIDs
    randnum_t GetRandNums() const {return mRandomNums;}

    //Return a friends list of a node
    std::set<coor_t> GetNodeNeig(size_t x, size_t y) const {return mGraph[x][y]->mFriends; }

    //Will return a node
    Node* GetNode(size_t x, size_t y) {return mGraph[x][y];}


    /* FUNCTIONS DEDICATED TO BE Setters */

    //Function will set votes for all nodes in graph
    void SetVotes(size_t x);

    //Function will set vote for one node
    void SetVote(size_t x, size_t y, size_t z);

    //Load the dna into all the hardware
    void SetGenome(program_t & pro);
    

    /* FUNCTIONS DEDICATED TO CLEAN UP CRAP */
    void CleanUp();

  private:
    /* GRAPH ITSELF AND ASSISANTS */

    //Holder of the graph object
    graph_t mGraph;
    //Holder of random numbers generated
    randnum_t mRandomNums;
    //Will hold the same pointers as mGraph, but only to set traits
    nodes_t mNodes;
    //Will hold the final votes at time called upon
    map_t mFinalVotes;
    //The scheduler
    std::vector<coor_t> mSchedule;

    /* GRAPH SPECIFIC PARAMATERS */

    //Dimension of the structure
    size_t GRA_DIM;
    //Number of iterations
    size_t NUM_ITER;
    //Number of Friends
    size_t NUM_FRI;
    //Number of enemies
    size_t NUM_ENE;
    //Random number generator
    emp::Ptr<emp::Random> mRng;
    //Random Number seed
    size_t RNG_SEED;
    //Minimum threshold
    double MIN_BIN_THSH;

    /* HARDWARE SPECIFIC PARAMATERS */

    //Position of UID within hw trait vector
    size_t UID;
    //Position of VOTE within hw trait vector
    size_t VOTE;
    //Position of X coordinate within hw trait vector
    size_t POSX;
    //Position of Y coordinate within hw trait vector
    size_t POSY;
    //Uperbound on random numbers
    size_t MAX_BND;
    //Lowerbound on ranodm numbers
    size_t MIN_BND;
    //Max Cores
    size_t MAX_CORES;
};

/* FUNCTIONS DEDICATED TO THE STRUCTURE */

//Function will create a general graph structure and 
//set the x and y position per hardware and wrong vote
//and will spawn a core (0, memory_t(), false)
void Graph::CreateGraph(size_t dim, size_t type, emp::Ptr<inst_lib_t> ilib, emp::Ptr<event_lib_t> elib)
{
  //0 => Toroidal Graph
  if(type == 0)
  {
    mGraph.resize(dim);
    
    for(size_t i = 0; i < mGraph.size(); ++i)
    {
      for(size_t j = 0; j < dim; ++j)
      {
        mSchedule.emplace_back(std::make_pair(i,j));

        emp::Ptr<Node> n =  emp::NewPtr<Node>(ilib, elib, mRng);
        n->mHW->SetMinBindThresh(MIN_BIN_THSH);
        n->mHW->SetTrait(POSX, i);
        n->mHW->SetTrait(POSY, j);
        n->mHW->SetMaxCores(100);
        mNodes.push_back(n);
        mGraph[i].push_back(n);
      }
    }

    if(mNodes.size() != (GRA_DIM * GRA_DIM))
    {
      std::cout << "WE NOT THE RIGHT SIZE" << std::endl;
    }
  }
}

//Will create adjacency list for each node
void Graph::CreateAdjList(size_t type, size_t dim)
{
  //0 => Toroidal Graph
  if(type == 0)
  {
    for(size_t i = 0; i < mGraph.size(); ++i)
    {
      for(size_t j = 0; j < mGraph[i].size(); ++j)
      {
        //Right
        mGraph[i][j]->mFriends.emplace(std::make_pair((((i+dim)+1) % dim), j));
        //Left
        mGraph[i][j]->mFriends.emplace(std::make_pair((((i+dim)-1) % dim), j));
        //Up
        mGraph[i][j]->mFriends.emplace(std::make_pair(i, (((j+dim)+1) % dim)));
        //Down
        mGraph[i][j]->mFriends.emplace(std::make_pair(i, (((j+dim)-1) % dim)));
      }
    }
  }
}

/* FUNCTIONS DEDICATED TO RUNNING EXPERIMENT */

//Give the graph NUM_ITER single processes to figure it out
double Graph::RunGraph(size_t iter)
{
  if(iter == -1)
    iter = NUM_ITER;

  double score = 0.0;

  for(size_t i = 0; i < iter; ++i)
  {
    emp::Shuffle(*mRng, mSchedule);

    for(auto p : mSchedule)
    {
      mGraph[p.first][p.second]->mHW->SingleProcess();
    }
    Graph::MakeFinalVotes();
    score += Graph::Consensus();
  }

  Graph::MakeFinalVotes();
  score += Graph::LegalVotes();
  score += Graph::LargestLegalVotes();  
  return score;
}

//Will see if element exits in vector
//True if inside mRandomNums
//False if not
bool Graph::Find(size_t me)
{
  for(auto num : mRandomNums)
  {
    if(num == me)
      return true;
  }

  return false;
}

//Will configure traits for hardware UID
void Graph::ConfigureTraits()
{
  mRandomNums.clear();

  while(mRandomNums.size() != (GRA_DIM * GRA_DIM))
  {
    size_t num = mRng->GetUInt(MIN_BND, MAX_BND);

    if(Graph::Find(num))
    {
      continue;
    }

    else
    {
      mRandomNums.push_back(num);
    }
  }
  
  if(mRandomNums.size() != mNodes.size())
  {
    std::cout << "Graph::ConfigureTraits() vector sizes not equal!" << std::endl;
    std::cout << "RandomNums: " << mRandomNums.size() << " mNodes: " << mNodes.size() << std::endl;
    exit(0);
  }

  for(size_t i = 0; i < mRandomNums.size(); ++i)
  {
    mNodes[i]->mHW->SetTrait(UID, mRandomNums[i]);
    mNodes[i]->mHW->SetTrait(VOTE, -999);
  }
}

//Will create a dictinary with legal votes that the hardware has come up with
void Graph::MakeFinalVotes()
{
  mFinalVotes.clear();

  for(size_t i = 0; i < mNodes.size(); ++i)
  {
    double vote = mNodes[i]->mHW->GetTrait(VOTE);

    if(Graph::Find(vote))
    {
      auto iter = mFinalVotes.find(vote);

      //Not found
      if(iter == mFinalVotes.end() )
      {
        mFinalVotes[vote] = 1.0;
      }

      else
      {
        mFinalVotes[vote] += 1.0;
      }
    }
  }
}

//Will return all legal votes by using mFinalVotes
double Graph::LegalVotes()
{
  double score = 0.0;

  for(auto p : mFinalVotes)
  {
    score += p.second;
  }

  return score;
}

//Will return the largest vote by using mFinalVotes
double Graph::LargestLegalVotes()
{
  double score = 0.0;

  for(auto p : mFinalVotes)
  {
    if(score < p.second)
    {
      score = p.second;
    }
  }

  return score;
}

//Will return if the system is in consensus or not by using mFinalVotes
double Graph::Consensus()
{
  double size = 0.0;

  for(auto p : mFinalVotes)
  {
    if(size < p.second)
    {
      size = p.second;
    }
  }

  if(size == (GRA_DIM * GRA_DIM))
  {
    return (double) (GRA_DIM * GRA_DIM);
  }

  else
  {
    return 0;
  }
}

//Will reset the graph to rerun with different program
void Graph::Reset()
{
  mRandomNums.clear();
  mFinalVotes.clear();
  Graph::ConfigureTraits();

  for(size_t i = 0; i < mNodes.size(); ++i)
  {
    mNodes[i]->mHW->ResetHardware();
    mNodes[i]->mHW->SpawnCore(0, memory_t(), true);
  }
}


/* FUNCTIONS DEDICATED TO BE Setters */

//Function to test scoring functions
void Graph::SetVotes(size_t x)
{
  for(size_t i = 0; i < mNodes.size(); ++i)
  {
    mNodes[i]->mHW->SetTrait(VOTE, x);
  }
}

//Function will set vote for one node
void Graph::SetVote(size_t x, size_t y, size_t z)
{
  mGraph[x][y]->mHW->SetTrait(VOTE, z);
}

//Load the dna into all the hardware
void Graph::SetGenome(program_t & pro)
{
  for(size_t i = 0; i < mGraph.size(); ++i)
  {
    for(size_t j = 0; j < mGraph.size(); ++j)
    {
      mGraph[i][j]->mHW->SetProgram(pro);
    }
  }
}


/* FUNCTIONS DEDICATED TO PRINTING OUT CRAP */

//Print out the triats
void Graph::PrintTraits()
{
  for(size_t i = 0; i < mGraph.size(); ++i)
  {
    for(size_t j = 0; j < mGraph.size(); ++j)
    {
      std::cout << "(" << i << "," << j << "): " << std::endl;
      std::cout << "UID: " << mGraph[i][j]->mHW->GetTrait(UID) << std::endl;
      std::cout << "VOTE: " << mGraph[i][j]->mHW->GetTrait(VOTE) << std::endl;
      std::cout << "POSX: " << mGraph[i][j]->mHW->GetTrait(POSX) << std::endl;
      std::cout << "POSY: " << mGraph[i][j]->mHW->GetTrait(POSY) << std::endl;
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

//Print the random UID
void Graph::PrintRanUID()
{
  std::cout << "Random UID: ";
  for(auto x : mRandomNums)
    std::cout << x << " ";
  std::cout << "\n" << std::endl;
}

//Print the FinalVotes vector
void Graph::PrintFinVotes()
{
  for(auto p : mFinalVotes)
  {
    std::cout << "Vote: " << p.first << " Count: " << p.second << std::endl;
  }
  std::cout << std::endl;
}

//Print Vote Trait
void Graph::PrintVotes()
{
  std::cout << "Hardware Votes: " << std::endl;
  for(size_t i = 0; i < mGraph.size(); ++i)
  {
    for(size_t j = 0; j < mGraph[i].size(); ++j)
    {
      auto node = mGraph[i][j]->mHW;
      auto s = node->GetCurState();
      std::cout << "(" << i << "," << j << "): " << node->GetTrait(VOTE) << " FP: " << s.func_ptr << " IP: " << s.inst_ptr << std::endl;
    }
  }
  std::cout << std::endl;
}

//Print Vote Trait
void Graph::PrintGenomes()
{
  std::cout << "Hardware Votes: " << std::endl;
  for(size_t i = 0; i < mGraph.size(); ++i)
  {
    for(size_t j = 0; j < mGraph[i].size(); ++j)
    {
      std::cout << "(" << i << "," << j << "): " << std::endl;
      mGraph[i][j]->mHW->PrintProgramFull();
    }
  }
  std::cout << std::endl;
}

//Print the mSchedule
void Graph::PrintSchedule()
{
  std::cout << "mSchedule.size(): " << mSchedule.size() << std::endl;
  for(auto p : mSchedule)
  {
    std::cout << "(" << p.first << ", " << p.second <<  ")" << std::endl;
  }
}

//Print all friends
void Graph::PrintFriends()
{
  for(size_t i = 0; i < mGraph.size(); ++i)
  {
    for(size_t j = 0; j < mGraph[i].size(); ++j)
    {
      std::cout << "(" << i << "," << j << "): "; 
      for(auto p : mGraph[i][j]->mFriends)
      {
        std::cout << "(" << p.first << ", " << p.second << "), ";
      }
      std::cout << std::endl;
    }
  }
}
#endif