#include <set>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <unordered_map>
#include "Map.h"
#include <cmath>
#include <algorithm>
#include "MR1PermutationPDB.h"
#include <math.h>

template<class Env, class State>
class TASS;

enum kAnchorSelection
{
    Temporal,
    Closest,
    Random
};

template<class State>
struct OpenClosedData
{
public:

    double g = 0;
    State parent;
    int dir;
    int loc = -1;
};

template<class State>
struct OpenData
{
public:
    State state;
	bool valid;
	uint64_t hash;
	OpenData()
	{

	}
	OpenData(State _state, bool _valid, uint64_t _hash)
	{
		state = _state;
		valid = _valid;
		hash = _hash;
	}
};

template <class Env, class State>
class TASSFrontier
{
public:
    Env *env;
	bool isIdle;
	unsigned int seed;
	std::vector<State> candidates;
	std::vector<State> otherCandidates;
	std::vector<State> children;
	std::vector<OpenData<State>> open;
	std::vector<State> firstPriority;
	//std::unordered_map<uint64_t, double> gValues;
	//std::unordered_map<uint64_t, State> parents;
	//std::unordered_map<uint64_t, int> loc;
    TASS<Env, State>* tass;
    int dir;
    
	std::vector<int> occurrences;
	std::vector<int> _occurrences;
	int numOfExp = 0;
	State rendezvous;
	TASSFrontier *other;
    int sampleCount;
	State anchor;
	State start;
    Heuristic<State> *h;
	TASSFrontier(Env *_env, State _start, Heuristic<State> *h, int _sampleCount, TASS<Env, State>* _tass, int _dir);
	~TASSFrontier(){}
	bool DoSingleSearchStep();
	void GetPath(std::vector<State> &path);
	void SetSeed(unsigned int seed);
	std::vector<State> GetPath(State node, bool forward);
	void ExtractPath(std::vector<State> &path);
	bool validSolution;
    double anchorH, anchorG;

    kAnchorSelection anchorSelection;

	//double comps = 0;



	double AdditiveDH(Heuristic<State> *heu, State s1, State s2)
	{
		double res = 0;
	    for (int i = 0; i < heu->heuristics.size(); i++)
	    {
	    	res += fabs(heu->heuristics[i]->HCost(s1, s1) - heu->heuristics[i]->HCost(s2, s2));
	    }
		return res;
	}

	double MaxDH(State s1, State s2)
	{
		return max(AdditiveDH(h, s1, s2), AdditiveDH(other->h, s1, s2));
	}

	double HCost(State s1, State s2)
	{
		//return abs(h->HCost(s1, s1) - h->HCost(s2, s2)) + abs(other->h->HCost(s1, s1) - other->h->HCost(s2, s2));

		return h->HCost(s1, s2); 
		//return MaxDH(s1, s2);
	}

};


template <class Env,  class State>
TASSFrontier<Env, State>::TASSFrontier(Env *_env, State _start, Heuristic<State> *_h, int _sampleCount, TASS<Env, State>* _tass, int _dir)
{
    dir = _dir;
    tass = _tass;
    sampleCount = _sampleCount;
	env = _env;
	open.resize(0);
	open.reserve(10000000);
	h = _h;
	start = _start;
    //openClosed.clear();
    //openClosed.reserve(10000000);
	auto shash = env->GetStateHash(start);
	open.push_back(OpenData<State>(start, true, shash));
    OpenClosedData<State> data;
    data.g = 0;
    data.loc = 0;
    data.parent = start;
    data.dir = dir;
    tass->openClosed.insert({shash, data});
	anchor = start;
	numOfExp = 0;
	validSolution = true;
}


template <class Env,  class State>
void TASSFrontier<Env, State>::ExtractPath(std::vector<State> &path)
{
	auto current = rendezvous;
	while (true)
	{
		path.push_back(current);
		if (tass->openClosed.find(env->GetStateHash(current)) != tass->openClosed.end())
			current = tass->openClosed[env->GetStateHash(current)].parent;
		else
			break;
	}
	path.push_back(start);
}




template <class Env,  class State>
bool TASSFrontier<Env, State>::DoSingleSearchStep()
{
	if (open.size() == 0)
	{
		//std::cout << "******* No Solution Found! ********" << std::endl;
		validSolution = false;
		return true;
	}
		
	double minDist = 99999999;
	auto data = open.back();
	State bestCandidate = data.state;
	uint64_t bestCandidateHash = data.hash;
	int _samples = sampleCount;
	int index = open.size() - 1;
    int bestIndex = index;
    while (index >= 0 && _samples > 0)
    {
        if (!open[index].valid)
        {
            open.pop_back();
            index--;
            continue;
        }
        State c = open[index].state;
        auto hash = open[index].hash;
		auto dist = HCost(c, other->anchor);
		if (dist < minDist || dist == minDist && tass->openClosed[hash].g > tass->openClosed[bestCandidateHash].g)
		{
			minDist = dist;
			bestCandidate = c;
			bestCandidateHash = hash;
            bestIndex = index;
		}
        index--;
		_samples--;
    }
	open[bestIndex].valid = false;
	tass->openClosed[bestCandidateHash].loc = -1;

	numOfExp++;
	children.resize(0);
	env->GetSuccessors(bestCandidate, children);
	for (State neighbor : children)
	{
		//occurrences[int(HCost(neighbor, other->anchor))]++;
		auto nhash = env->GetStateHash(neighbor);
		double g = tass->openClosed[bestCandidateHash].g + env->GCost(bestCandidate, neighbor);
        auto ent = tass->openClosed.find(nhash);
		if (ent != tass->openClosed.end())
		{
            if (ent->second.dir != dir)
            {
                rendezvous = bestCandidate;
		        other->rendezvous = neighbor;
		        return true;
            }
            else
            {
                if (g < ent->second.g)
			    {
			    	ent->second.g = g;
			    	ent->second.parent = bestCandidate;
			    }
			    else if (ent->second.loc != -1)
			    {
                    open.push_back(open[ent->second.loc]);
                    open[ent->second.loc].valid = false;
			        ent->second.loc = open.size() - 1;
			    }
            }
		}
		else
		{
			open.push_back(OpenData<State>(neighbor, true, nhash));
            OpenClosedData<State> data;
            data.g = g;
            data.loc = open.size() - 1;
            data.parent = bestCandidate;
            data.dir = dir;
            tass->openClosed[nhash] = data;
		}
	}

    switch (anchorSelection)
    {
        case Temporal:
            anchor = bestCandidate;
            break;
        case Closest:
		{
			anchorH = HCost(anchor, other->start);
			auto hh = HCost(bestCandidate, other->start);
            if (hh < anchorH)
            {
                anchor = bestCandidate;
                anchorH = hh;
                anchorG = tass->openClosed[bestCandidateHash].g;
            }
            else if (hh == anchorH)
            {
                if (tass->openClosed[bestCandidateHash].g < anchorG)
                {
                    anchor = bestCandidate;
                    anchorG = tass->openClosed[bestCandidateHash].g;
                }
            }
            break;
		}
        case Random:
		{
			int random_index = rand_r(&seed) % open.size();
            anchor = open[random_index].state;
            break;
		}
		default:
			break;
    }

    
	//anchor = bestCandidate;
    /*
	if (other->openClosed.find(bestCandidateHash) != other->openClosed.end())// || bestCandidate == other->start)
	{
		rendezvous = bestCandidate;
		other->rendezvous = rendezvous;
		return true;
	}
    */
	return false;
}


template <class Env,  class State>
void TASSFrontier<Env, State>::SetSeed(unsigned int seed)
{
	this->seed = seed;
}

template <class Env, class State>
class TASS
{
private:
	int turn = 0;
public:
    std::unordered_map<uint64_t, OpenClosedData<State>> openClosed;
	int episode = -1;
    double pathRatio;
	TASSFrontier<Env, State>* ff;
	TASSFrontier<Env, State>* bf;
	TASS();
	TASS(Env *_env, State _start, State _goal, Heuristic<State> *hf, Heuristic<State> *hb, int _sampleCount);
	~TASS(){}
	void Init(Env *_env, State _start, State _goal, Heuristic<State> *hf, Heuristic<State> *hb, int _sampleCount)
	{
		ff = new TASSFrontier<Env, State>(_env, _start, hf, _sampleCount, this);
		bf = new TASSFrontier<Env, State>(_env, _goal, hb, _sampleCount, this);
		ff->other = bf;
		bf->other = ff;
		ff.anchorH = ff.HCost(_start, _goal);
		bf.anchorH = bf.HCost(_goal, _start);
	}
	std::function<void(State, State, Heuristic<State>*&)> heuristicUpdate;
	void SetHeuristicUpdate(std::function<void(State, State, Heuristic<State>*&)> p)
	{
		heuristicUpdate = p;
	}
	int GetPathLength()
	{

	}
	int GetNodesExpanded()
	{
		return ff->numOfExp + bf->numOfExp;
	}
	bool DoSingleSearchStep()
	{
		if (turn == 0)
		{
			turn = 1;
			return ff->DoSingleSearchStep();
		}
		else
		{
			turn = 0;
			return bf->DoSingleSearchStep();
		}
	}
	void ExtractPath(std::vector<State> &path)
	{
		auto current = bf->rendezvous;
		std::vector<State> front, back;
		while (true)
		{
			back.push_back(current);
			if (openClosed.find(bf->env->GetStateHash(current))->second.parent != current)
				current = openClosed[bf->env->GetStateHash(current)].parent;
			else
				break;
		}
		back.push_back(bf->start);
		current = ff->rendezvous;
		while (true)
		{
			front.push_back(current);
			if (openClosed.find(ff->env->GetStateHash(current))->second.parent != current)
				current = openClosed[ff->env->GetStateHash(current)].parent;
			else
				break;
		}
		front.push_back(ff->start);
		path.resize(0);
        pathRatio = (double)front.size() / (double)(front.size() + back.size());
		for (int i = front.size() - 1; i >= 0; i--)
			path.push_back(front[i]);
		for (int i = 0; i < back.size(); i++)
			path.push_back(back[i]);
	}
	void GetPath(std::vector<State> &path)
	{
		path.resize(0);
		int count = 0;
		while (true)
		{
			if (episode > 0 && count % episode == 0)
			{
				//Heuristic<State> *h1 = new Heuristic<State>();
				heuristicUpdate(ff->anchor, bf->anchor, bf->h);
				//bf->h = h1;
				//std::cout << "TEST: " << bf->h->HCost(ff->anchor, bf->anchor) << std::endl;

				//Heuristic<State> *h2 = new Heuristic<State>();
				//heuristicUpdate(bf->anchor, ff->h);
				ff->h = bf->h;
				//std::cout << "TEST: " << bf->h->HCost(ff->anchor, bf->anchor) << std::endl;
				std::cout << "episode(" << episode << "): " << count / episode << std::endl;
				std::cout << ff->HCost(ff->anchor, bf->start) << " " << ff->anchorH << std::endl;
			}
			count++;
			if (DoSingleSearchStep())
			{
				break;
			}
		}
		if (!ff->validSolution || !bf->validSolution)
			return;
		ExtractPath(path);
		openClosed.clear();
		ff->open.clear();
		bf->open.clear();
		delete ff;
		delete bf;
	}
	void SetSeed(unsigned int seed)
	{
		this.seed = seed;
	}

    void SetAnchorSelection(kAnchorSelection selection)
    {
        ff->anchorSelection = selection;
        bf->anchorSelection = selection;
    }
};


template <class Env, class State>
TASS<Env, State>::TASS(Env *_env, State _start, State _goal, Heuristic<State> *hf, Heuristic<State> *hb, int _sampleCount)
{
    openClosed.clear();
    openClosed.reserve(100000);
	ff = new TASSFrontier<Env, State>(_env, _start, hf, _sampleCount, this, 1);
	bf = new TASSFrontier<Env, State>(_env, _goal, hb, _sampleCount, this, 2);
	ff->other = bf;
	bf->other = ff;
}