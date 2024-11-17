/*
 *  hog2
 *
 *  Created by Nathan Sturtevant on 12/13/22.
 *
 * This file is part of HOG2. See https://github.com/nathansttt/hog2 for licensing information.
 *
 */

#include "Common.h"
#include "Driver.h"
#include "Map2DEnvironment.h"
#include "Racetrack.h"
#include "TemplateAStar.h"
#include "ScenarioLoader.h"
#include "FPUtil.h"
#include "Graphics.h"
#include "SVGUtil.h"
#include "DSDWAStar.h"
#include "MapGenerators.h"
#include "DynamicPotentialSearch.h"
#include <ctime>
#include <cmath>
#include "GridHeuristics.h"
#include "MNPuzzle.h"
#include "STPInstances.h"

int stepsPerFrame = 1;
float bound = 4;
int problemNumber = 0;
float testScale = 1.0;
void GetNextWeightRange(float &minWeight, float &maxWeight, point3d currPoint, float nextSlope);
float GetPriority(float h, float g);
float ChooseWeightForTargetPriority(point3d point, float priority, float minWeight, float maxWeight, point3d last, float &K);
DSDWAStar<xyLoc, tDirection, MapEnvironment> dsd;
TemplateAStar<xyLoc, tDirection, MapEnvironment> tas;
DynamicPotentialSearch<xyLoc, tDirection, MapEnvironment> dps;
std::vector<xyLoc> solution;
std::vector<xyLoc> theList;
Map * m = 0;
MapEnvironment * me = 0;
ScenarioLoader * sl = 0;
xyLoc start, goal, swampedloc;
TemplateAStar<RacetrackState, RacetrackMove, Racetrack> tas_track;
DSDWAStar<RacetrackState, RacetrackMove, Racetrack> dsd_track;
DynamicPotentialSearch<RacetrackState, RacetrackMove, Racetrack> dps_track;
std::vector<RacetrackState> path;
Racetrack *r = 0;
RacetrackState from;
RacetrackState end;
int numExtendedGoals = 300;
int exper=8; //Exper 8 is to make DW domains. Use any other number for other domains.
float terrain_size=30, terrain_width=10, terrain_height=10;
double Tcosts[4], rdm, hardness[4];
bool showPlane = false;
bool searchRunning = false;
bool saveSVG = false;
bool useDH = false;
bool useLookUpTable = false;
bool limitScenarios = false;
int numLimitedScenarios = 1000, lowerLimit=50, upperLimit=2500, numScenario=598;
int randomIndex;
xyLoc xyLocRandomState;
MNPuzzleState<4, 4> mnpRandomState;
GridEmbedding *dh;

// 1=DSD random room map, 2=DSD random map, 3=DSD random maze, 4=DSD designed map,
// 5=DSD Load map, 6=DSD Load RaceTrack, 7=DPS Load map, 8=DPS Load RaceTrack
int mapcmd = 2;

//std::string mapload = "mazes/maze512-32-6";
 std::string mapload = "dao/ost003d";
// std::string mapload = "dao/orz000d";
// std::string mapload = "dao/den520d";
// std::string mapload = "dao/arena";
// std::string mapload = "da2/ca_caverns1";
// std::string mapload = "da2/lt_undercityserialkiller";
// std::string mapload = "da2/lt_foundry_n";

//#include "Plot2D.h"
struct DSDdata {
    float slope;
    float weight;
    float K;
    point3d crossPoint; // cached for simplicity
};
std::vector<DSDdata> data;

struct DSDdata_v2 {
    float weight;
    float K;
    point3d crossPoint; // cached for simplicity
};
std::unordered_map<double, DSDdata_v2> look_up_table;
double last_in_lookup;

std::vector<DSDdata_v2> LookUpVector;

int main(int argc, char* argv[])
{
    setvbuf(stdout, NULL, _IONBF, 0);
    InstallHandlers();
    RunHOGGUI(argc, argv, 1500, 1500);
    return 0;
}

/**
 * Allows you to install any keyboard handlers needed for program interaction.
 */
void InstallHandlers()
{
    InstallKeyboardHandler(MyDisplayHandler, "Reset lines", "Reset incremenetal lines", kAnyModifier, 'r');
    InstallKeyboardHandler(MyDisplayHandler, "Show plane", "Show gradient of plane heights", kAnyModifier, 'p');
    InstallKeyboardHandler(MyDisplayHandler, "Faster", "Speed up search animation", kAnyModifier, ']');
    InstallKeyboardHandler(MyDisplayHandler, "Slower", "Slow down search animation", kAnyModifier, '[');
    InstallKeyboardHandler(MyDisplayHandler, "Policy", "Increment policy", kAnyModifier, '}');
    InstallKeyboardHandler(MyDisplayHandler, "Policy", "Decrement policy", kAnyModifier, '{');
    InstallKeyboardHandler(MyDisplayHandler, "Problem", "Increment problem", kAnyModifier, '+');
    InstallKeyboardHandler(MyDisplayHandler, "Problem", "Decrement problem", kAnyModifier, '-');
    InstallKeyboardHandler(MyDisplayHandler, "Swamped Problems", "Increment swamped problem", kAnyModifier, 's');
    InstallKeyboardHandler(MyDisplayHandler, "Bound", "Increment bound", kAnyModifier, 'w');
    InstallKeyboardHandler(MyDisplayHandler, "toggle_lookUp_Table", "toggle lookUp Table", kAnyModifier, 't');

    InstallCommandLineHandler(MyCLHandler, "-stpDSD", "-stpDSD <problem> <alg> <weight> <puzzleW>", "Test STP <problem> <algorithm> <weight> <puzzleW>");
    InstallCommandLineHandler(MyCLHandler, "-stpBaseLines", "-stpBaseLines <problem> <alg> <weight> <puzzleW>", "Test STP <problem> <algorithm> <weight> <puzzleW>");
    InstallCommandLineHandler(MyCLHandler, "-gridDSD", "-gridDSD <map> <scenario> <alg> <weight> <TerrainSize> <SwampHardness> <Experiment>", "Test grid <map> on <scenario> with <algorithm> <weight> <TerrainSize> <SwampHardness> and <Experiment>");
    InstallCommandLineHandler(MyCLHandler, "-gridDPS", "-gridDPS <map> <scenario> <weight> <TerrainSize> <SwampHardness> <Experiment>", "Test grid <map> on <scenario> with <weight> <TerrainSize> <SwampHardness> and <Experiment>");
    InstallCommandLineHandler(MyCLHandler, "-gridBaseLines", "-gridBaseLines <map> <scenario> <alg> <weight> <TerrainSize> <SwampHardness> <Experiment>", "Test grid <map> on <scenario> with <algorithm> <weight> <TerrainSize> <SwampHardness> and <Experiment>");
    InstallCommandLineHandler(MyCLHandler, "-rtBaseLines", "-rtBaseLines <map> <scenario> <alg> <weight> <numExtendedGoals>", "Test grid <map> on <scenario> with <algorithm> <weight> <numExtendedGoals>");
    InstallCommandLineHandler(MyCLHandler, "-rtDSD", "-rtDSD <map> <scenario> <alg> <weight> <numExtendedGoals>", "Test grid <map> on <scenario> with <algorithm> <weight> <numExtendedGoals>");
    InstallCommandLineHandler(MyCLHandler, "-exp0BaseLines", "-exp0BaseLines <alg> <weight> <TerrainSize> <mapType> <SwampHardness>", "Test grid with <algorithm> <weight> <TerrainSize> <mapType> <SwampHardness>");
    InstallCommandLineHandler(MyCLHandler, "-exp0DSD", "-exp0DSD <alg> <weight> <TerrainSize> <mapType> <SwampHardness>", "Test grid with <algorithm> <weight> <TerrainSize> <mapType> <SwampHardness>");
    InstallCommandLineHandler(MyCLHandler, "-exp0DPS", "-exp0DPS <weight> <TerrainSize> <mapType> <SwampHardness>", "Test grid <weight> <TerrainSize> <mapType> <SwampHardness>");
    InstallCommandLineHandler(MyCLHandler, "-stpAstar", "-stpAstar <problem> <alg> <weight>", "Test STP <problem> <algorithm> <weight>");
    InstallCommandLineHandler(MyCLHandler, "-stpDPS", "-stpDPS problem weight puzzleW", "Test STP <problem> <weight> <puzzleW>");
    InstallCommandLineHandler(MyCLHandler, "-rtDPS", "-rtDPS <map> <scenario> <weight> <numExtendedGoals>", "Test grid <map> on <scenario> with <weight> <numExtendedGoals>");
    InstallCommandLineHandler(MyCLHandler, "-timeDSWA", "-timeDSWA stp problem weight", "Test STP <problem> <weight>");
    InstallCommandLineHandler(MyCLHandler, "-map", "-map <map> <scenario> alg weight", "Test grid <map> on <scenario> with <algorithm> <weight>");
    
    InstallWindowHandler(MyWindowHandler);
    
    InstallMouseClickHandler(MyClickHandler);
}

void MyWindowHandler(unsigned long windowID, tWindowEventType eType){
    if (eType == kWindowDestroyed)
    {
        printf("Window %ld destroyed\n", windowID);
        RemoveFrameHandler(MyFrameHandler, windowID, 0);
    }
    else if (eType == kWindowCreated){
        printf("Window %ld created\n", windowID);
        InstallFrameHandler(MyFrameHandler, windowID, 0);
        ReinitViewports(windowID, {-1, -1, 0, 1}, kScaleToSquare);
        AddViewport(windowID, {0, -1, 1, 1}, kScaleToSquare);
        srandom(20221228);

        //load map or generate a random one.
        if(mapcmd==1){
            m = new Map(200, 200);
            BuildRandomRoomMap(m, 30);
            // default 8-connected with ROOT_TWO edge costs
            me = new MapEnvironment(m);
        }
        else if(mapcmd==2){
            m = new Map(200, 200);
            MakeRandomMap(m, 10);
            // default 8-connected with ROOT_TWO edge costs
            me = new MapEnvironment(m);
        }
        else if(mapcmd==3){
            m = new Map(200, 200);
            MakeMaze(m, 10);
            // default 8-connected with ROOT_TWO edge costs
            me = new MapEnvironment(m);
        }
        else if(mapcmd==4){
            m = new Map(200, 200);
             MakeDesignedMap(m, 30, 6);
            // default 8-connected with ROOT_TWO edge costs
            me = new MapEnvironment(m);
        }
        else if(mapcmd>=5){
            std::string tmpscenload = "./Users/mohammadrezahami/Documents/University/Project/hog2-PDB-refactor/DynamicPhiFunction/scenarios/"+ mapload + ".map.scen";
            std::string tmpmapload = "./Users/mohammadrezahami/Documents/University/Project/hog2-PDB-refactor/DynamicPhiFunction/maps/"+ mapload + ".map";
            sl = new ScenarioLoader (tmpscenload.c_str());
            m = new Map(tmpmapload.c_str());
            me = new MapEnvironment(m);
            r = new Racetrack(m);
        }

        //load the scenario, or initiate it randomly.
        if(mapcmd <=4){
            start = {1,99};
            goal = {198, 99};

            me->SetDiagonalCost(1.41);
            // me->SetDiagonalCost(1.5);

            dsd.policy = kWA;
            me->SetInputWeight(bound);

            swampedloc = {1,1};

            //Set the cost of each terrain type randomly.
            for(int i=0; i<4; i++){
                //[0]=kSwamp, [1]=kWater,[2]=kGrass, [3]=kTrees
                string type;
                if(i==0) type="Swamp";
                if(i==1) type="Water";
                if(i==2) type="Grass";
                if(i==3) type="Trees";

                // Define the Hardness
                // 1. Random Hardness
                // rdm = random()%101;
                // hardness[i] = rdm/101+1;
                //Or 2. Hardcoded Hardness
                hardness[0]=2; hardness[1]=1.45; hardness[2]=1.25; hardness[3]=1.95;
                
                // Use Hardness to define Cost
                // 1. Hardness with respect to input w
                Tcosts[i] = hardness[i]*(me->GetInputWeight())-(hardness[i]-1);
                // std::cout<<"Cost "<<type<<"="<<Tcosts[i]<<" ("<<hardness[i]<<"*"<<me->GetInputWeight()<<"-"<<(hardness[i]-1)<<")"<<std::endl;
                //Or 2. Hardness as the exact Cost
                // Tcosts[i] = hardness[i];
                //Or 3. The exact Cost hardcoded
                // Tcosts[i] = 4.5;
                std::cout<<"Cost "<<type<<"="<<Tcosts[i]<<std::endl;
            }
            me->SetTerrainCost(Tcosts);

            dsd.SetWeight(bound);

            if(useLookUpTable){
                // look_up_table.clear();
                // dsd.InitializeSearch_v2(me, start, goal, solution);
                LookUpVector.clear();
                dsd.InitializeSearch_v3(me, start, goal, solution);

            }
            else{
                data.resize(0);
                dsd.InitializeSearch(me, start, goal, solution);
            }
            
            if(useDH){
                dh = new GridEmbedding(me, 10, kLINF);
                for (int x = 0; x < 10; x++)
                    dh->AddDimension(kDifferential, kFurthest);
                dsd.SetHeuristic(dh);
            }
        }
        else if(mapcmd == 5 || mapcmd == 7){
            //Grid Map
            Experiment exp = sl->GetNthExperiment(numScenario);
            start.x = exp.GetStartX();
            start.y = exp.GetStartY();
            goal.x = exp.GetGoalX();
            goal.y = exp.GetGoalY();

            me->SetDiagonalCost(1.41);
            // me->SetDiagonalCost(1.5);

            me->SetInputWeight(bound);

            swampedloc = {1,1};

            //Set the cost of each terrain type randomly.
            for(int i=0; i<4; i++){
                //[0]=kSwamp, [1]=kWater,[2]=kGrass, [3]=kTrees
                string type;
                if(i==0) type="Swamp";
                if(i==1) type="Water";
                if(i==2) type="Grass";
                if(i==3) type="Trees";

                // Define the Hardness
                // 1. Random Hardness
                // rdm = random()%101;
                // hardness[i] = rdm/101+1;
                //Or 2. Hardcoded Hardness
                hardness[0]=1.5; hardness[1]=1.45; hardness[2]=1.25; hardness[3]=1.95;
                
                // Use Hardness to define Cost
                // 1. Hardness with respect to input w
                Tcosts[i] = hardness[i]*(me->GetInputWeight())-(hardness[i]-1);
                // std::cout<<"Cost "<<type<<"="<<Tcosts[i]<<" ("<<hardness[i]<<"*"<<me->GetInputWeight()<<"-"<<(hardness[i]-1)<<")"<<std::endl;
                //Or 2. Hardness as the exact Cost
                // Tcosts[i] = hardness[i];
                //Or 3. The exact Cost hardcoded
                // Tcosts[i] = 4.5;
                std::cout<<"Cost "<<type<<"="<<Tcosts[i]<<std::endl;
            }
            me->SetTerrainCost(Tcosts);

            if(useDH){
                dh = new GridEmbedding(me, 10, kLINF);
                for (int x = 0; x < 10; x++)
                    dh->AddDimension(kDifferential, kFurthest);
                dsd.SetHeuristic(dh);
            }
            
            if(mapcmd == 5){
                dsd.policy = kWA;
                dsd.SetWeight(bound);
                dsd.InitializeSearch(me, start, goal, solution);
            }
            else if(mapcmd == 7){
                dps.SetOptimalityBound(bound);
                dps.InitializeSearch(me, start, goal, solution);
            }
        
        }
        else if(mapcmd == 6 || mapcmd == 8){
            //RaceTrack
            Experiment exp = sl->GetNthExperiment(numScenario);
            std::cout<<"Path Length is "<<exp.GetDistance()<<" (scen "<<numScenario<<"/"<<sl->GetNumExperiments()<<")\n";
            if(numExtendedGoals) numExtendedGoals = exp.GetDistance()/3;

            //Set the start and goal states, weight and policy of the search.
            from.xLoc = exp.GetStartX();
            from.yLoc = exp.GetStartY();
            from.xVelocity = 0;
            from.yVelocity = 0;
            end.xLoc = exp.GetGoalX();
            end.yLoc = exp.GetGoalY();
            end.xVelocity = 0;
            end.yVelocity = 0;
            
            //Set the start and goal TerrainType.
            m->SetTerrainType(from.xLoc, from.yLoc, kStartTerrain);
            
            if(numExtendedGoals == 0){//Not extend the goal.
                m->SetTerrainType(end.xLoc, end.yLoc, kEndTerrain);
            }
            else{//Extend the goal
                tas.SetPhi([=](double h,double g){return g;});
                goal.x = end.xLoc;
                goal.y = end.yLoc;
                tas.ExtendGoal(me, goal, theList, numExtendedGoals);
                for(int tNode=0; tNode<theList.size(); tNode++){
                    m->SetTerrainType(theList[tNode].x, theList[tNode].y, kEndTerrain);
                }
            }
            
            r->UpdateMap(m);

            if(mapcmd == 6){
                dsd_track.policy = kWA;
                dsd_track.SetWeight(bound);
                dsd_track.InitializeSearch(r, from, end, path);
            }
            else if(mapcmd == 8){
                dps_track.SetOptimalityBound(bound);
                dps_track.InitializeSearch(r, from, end, path);
            }
        }

        searchRunning = true;

    }

}

point3d zero(0, 0);
point3d origin(-1, 1);

point3d HOGToLocal(point3d p)
{
    return point3d((p.x+1)*bound/2.0f , (p.y-1)*bound/-2.0);
}

point3d LocalToHOG(point3d p)
{
    return point3d(2.0f*(p.x/bound-0.5), -2.0*(p.y/bound-0.5));
}

// 1. resize window to 1024x768
void MyFrameHandler(unsigned long windowID, unsigned int viewport, void *)
{
    Graphics::Display &display = getCurrentContext()->display;
    if(mapcmd <= 5 || mapcmd==7){
        if (viewport == 0)
        {
            if (me)
                me->Draw(display);
            if (searchRunning)
            {
                for (int x = 0; x < stepsPerFrame; x++)
                {
                    if (solution.size() == 0)
                    {
                        if(mapcmd <= 5 && useLookUpTable){
                            // if (dsd.DoSingleSearchStep_v2(solution)){
                            //     std::cout << "Node Expansions: " << dsd.GetNodesExpanded() << "\n";
                            // }
                            if (dsd.DoSingleSearchStep_v3(solution)){
                                std::cout << "Node Expansions: " << dsd.GetNodesExpanded() << "\n";
                            }
                        }
                        if(mapcmd <= 5 && !useLookUpTable){
                            if (dsd.DoSingleSearchStep(solution)){
                                std::cout << "Node Expansions: " << dsd.GetNodesExpanded() << "\n";
                            }
                        }
                        else if(mapcmd == 7){
                            if (dps.DoSingleSearchStep(solution)){
                            std::cout << "Node Expansions: " << dps.GetNodesExpanded() << "\n";
                            }
                        }
                    }
                }
                if(mapcmd <= 5){
                    dsd.Draw(display);
                }
                else if(mapcmd == 7){
                    dps.Draw(display);
                }

            }
        }
        if (viewport == 1)
        {
            if (searchRunning)
            {
                // if(mapcmd <= 5){
                //     dsd.DrawPriorityGraph(display);
                // }
                if(mapcmd <= 5 && !useLookUpTable){
                    dsd.DrawPriorityGraph(display);
                }
                if(mapcmd <= 5 && useLookUpTable){
                    // dsd.DrawPriorityGraph_v2(display);
                    dsd.DrawPriorityGraph_v3(display);
                }
                else if(mapcmd == 7){
                    //Does not apply here
                    // dps.DrawPriorityGraph(display);
                }
            }
        }
        if (viewport == 3) // TODO: Add mode for exploration
        {
            if (showPlane)
            {
                //        float divisor = GetPriority(bound, bound);
                float divisor = GetPriority(1, 0);
                if (isinf(divisor))
                    divisor = 4;
                const float delta = 0.015;//0.025;
                for (float x = 0; x < bound; x+=delta)
                {
                    for (float y = 0; y < bound; y+=delta)
                    {
                        display.FillSquare(LocalToHOG({x, y}), delta/2, rgbColor::hsl(fmodf(GetPriority(x, y)/divisor, 1.0f), 1.0, 0.5));
                    }
                }
            }
            
            //    point3d lastCrossPoint(1/bound, 0);
            float priority = 1;//HOGToLocal({-1, -1}).y;//2.0f/bound;
            
            // draw bounding line
            point3d bl1(priority, 0), bl2(0, priority*bound); // main suboptimality line
            point3d bl3(priority*bound, 0);//(0+2, priority*bound-2); //
            point3d bl4(priority*bound/(2*bound-1), 0);//priority*bound-(2*bound-1));
            point3d bl2a(0, priority);
            point3d bl2c(priority-priority*bound/(2*bound-1), priority*bound);
            // main priority line
            display.DrawLine(LocalToHOG(bl1), LocalToHOG(bl2), 1/100.0f, Colors::yellow);
            display.DrawLine(LocalToHOG(bl1*testScale), LocalToHOG(bl2*testScale), 1/100.0f, Colors::darkyellow);
            // 45° upper bound line
            display.DrawLine(LocalToHOG(bl2), LocalToHOG(bl3), 1/100.0f, Colors::darkgray);
            // 2w-1 upper bound line
            display.DrawLine(LocalToHOG(bl1), LocalToHOG(bl2c), 1/100.0f, Colors::darkgray);
            
            // 45° lower bound line
            display.DrawLine(LocalToHOG(bl2a), LocalToHOG(bl1), 1/100.0f, Colors::lightgray);
            // 2w-1 lower bound line
            display.DrawLine(LocalToHOG(bl2), LocalToHOG(bl4), 1/100.0f, Colors::lightgray);
            
            // Draw actual priority line across
            for (int x = 0; x < data.size(); x++)
            {
                point3d value = origin;
                
                // y = slope * x // x=1 -> y = slope; y=1 -> x = 1/slope;
                if (data[x].slope < 1)
                {
                    value.x += 2;
                    value.y -= 2*data[x].slope;
                }
                else {
                    value.x += 2/data[x].slope;
                    value.y -= 2;
                }
                display.DrawLine(origin, value, 1/200.0f, Colors::blue);
                
                point3d crossPoint1, crossPoint2;
                crossPoint1.x = priority/(data[x].K*(data[x].slope+data[x].weight));
                crossPoint1.y = crossPoint1.x*data[x].slope;
                float lastSlope = ((x==0)?(0):(data[x-1].slope));
                crossPoint2.x = priority/(data[x].K*(lastSlope+data[x].weight));
                crossPoint2.y = crossPoint2.x*lastSlope;
                display.DrawLine(LocalToHOG(crossPoint1), LocalToHOG(crossPoint2), 1/100.0f, Colors::red);
            }
            for (int x = 0; x < data.size(); x++)
                display.FillCircle(LocalToHOG(data[x].crossPoint), 0.01, Colors::darkgreen);
            
            display.DrawLine(origin, {1, 1}, 1./100.0f, Colors::white);
            display.DrawLine(origin, {-1, -1}, 1./100.0f, Colors::white);
        }
    }
    else if(mapcmd == 6 || mapcmd == 8){
        //Racetrack
        if (viewport == 0)
        {
            if (r){
                r->Draw(display);
            }
            if (searchRunning)
            {
                for (int x = 0; x < stepsPerFrame; x++)
                {
                    if (path.size() == 0)
                    {
                        if(mapcmd == 6){
                            if (dsd_track.DoSingleSearchStep(path)){
                            std::cout << "Node Expansions: " << dsd_track.GetNodesExpanded() << "\n";
                            searchRunning = false;
                            }
                        }
                        else if(mapcmd == 8){
                            if (dps_track.DoSingleSearchStep(path)){
                            std::cout << "Node Expansions: " << dps_track.GetNodesExpanded() << "\n";
                            searchRunning = false;
                            }
                        }
                    }
                }
            }
            if(mapcmd == 6){
                dsd_track.Draw(display);
            }
            else if(mapcmd == 8){
                dps_track.Draw(display);
            }
        }
        if (viewport == 1){
            if (searchRunning)
            {
                if(mapcmd == 6){
                    dsd_track.DrawPriorityGraph(display);
                }
                else if(mapcmd == 8){
                    //Does not apply here
                    // dps_track.DrawPriorityGraph(display);
                }
            }
        }
        if (viewport == 3) // TODO: Add mode for exploration
        {
            if (showPlane)
            {
                //        float divisor = GetPriority(bound, bound);
                float divisor = GetPriority(1, 0);
                if (isinf(divisor))
                    divisor = 4;
                const float delta = 0.015;//0.025;
                for (float x = 0; x < bound; x+=delta)
                {
                    for (float y = 0; y < bound; y+=delta)
                    {
                        display.FillSquare(LocalToHOG({x, y}), delta/2, rgbColor::hsl(fmodf(GetPriority(x, y)/divisor, 1.0f), 1.0, 0.5));
                    }
                }
            }
            
            //    point3d lastCrossPoint(1/bound, 0);
            float priority = 1;//HOGToLocal({-1, -1}).y;//2.0f/bound;
            
            // draw bounding line
            point3d bl1(priority, 0), bl2(0, priority*bound); // main suboptimality line
            point3d bl3(priority*bound, 0);//(0+2, priority*bound-2); //
            point3d bl4(priority*bound/(2*bound-1), 0);//priority*bound-(2*bound-1));
            point3d bl2a(0, priority);
            point3d bl2c(priority-priority*bound/(2*bound-1), priority*bound);
            // main priority line
            display.DrawLine(LocalToHOG(bl1), LocalToHOG(bl2), 1/100.0f, Colors::yellow);
            display.DrawLine(LocalToHOG(bl1*testScale), LocalToHOG(bl2*testScale), 1/100.0f, Colors::darkyellow);
            // 45° upper bound line
            display.DrawLine(LocalToHOG(bl2), LocalToHOG(bl3), 1/100.0f, Colors::darkgray);
            // 2w-1 upper bound line
            display.DrawLine(LocalToHOG(bl1), LocalToHOG(bl2c), 1/100.0f, Colors::darkgray);
            
            // 45° lower bound line
            display.DrawLine(LocalToHOG(bl2a), LocalToHOG(bl1), 1/100.0f, Colors::lightgray);
            // 2w-1 lower bound line
            display.DrawLine(LocalToHOG(bl2), LocalToHOG(bl4), 1/100.0f, Colors::lightgray);
            
            // Draw actual priority line across
            for (int x = 0; x < data.size(); x++)
            {
                point3d value = origin;
                
                // y = slope * x // x=1 -> y = slope; y=1 -> x = 1/slope;
                if (data[x].slope < 1)
                {
                    value.x += 2;
                    value.y -= 2*data[x].slope;
                }
                else {
                    value.x += 2/data[x].slope;
                    value.y -= 2;
                }
                display.DrawLine(origin, value, 1/200.0f, Colors::blue);
                
                point3d crossPoint1, crossPoint2;
                crossPoint1.x = priority/(data[x].K*(data[x].slope+data[x].weight));
                crossPoint1.y = crossPoint1.x*data[x].slope;
                float lastSlope = ((x==0)?(0):(data[x-1].slope));
                crossPoint2.x = priority/(data[x].K*(lastSlope+data[x].weight));
                crossPoint2.y = crossPoint2.x*lastSlope;
                display.DrawLine(LocalToHOG(crossPoint1), LocalToHOG(crossPoint2), 1/100.0f, Colors::red);
            }
            for (int x = 0; x < data.size(); x++)
                display.FillCircle(LocalToHOG(data[x].crossPoint), 0.01, Colors::darkgreen);
            
            display.DrawLine(origin, {1, 1}, 1./100.0f, Colors::white);
            display.DrawLine(origin, {-1, -1}, 1./100.0f, Colors::white);
        }
    }
}

int MyCLHandler(char *argument[], int maxNumArgs)
{
    if (strcmp(argument[0], "-stpBaseLines") == 0)
    {
        assert(maxNumArgs >= 5);
        
        TemplateAStar<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>> tas_mnp;
        std::vector<MNPuzzleState<4, 4>> path;
        MNPuzzle<4, 4> mnp;
        MNPuzzleState<4, 4> start = STP::GetKorfInstance(atoi(argument[1]));
        MNPuzzleState<4, 4> goal;
        //random.randint(10, int(data[10])-10)
        std::vector<int> hardcodedNumbers = {23, 12, 25, 20, 31, 14, 14, 10, 13, 13, 14, 15, 23, 20, 21, 12, 13, 28, 23, 17, 14, 22, 12, 22, 12, 12, 12, 22, 28, 25, 15, 13, 27, 21, 10, 19, 20, 15, 21, 22, 15, 12, 36, 20, 15, 19, 23, 12, 12, 10, 12, 25, 22, 11, 14, 13, 12, 23, 15, 28, 12, 13, 17, 10, 11, 16, 11, 19, 21, 17, 11, 17, 14, 24, 16, 24, 15, 23, 13, 15, 28, 22, 20, 13, 11, 12, 21, 20, 19, 22, 16, 14, 15, 15, 21, 17, 22, 12, 27, 13 };

        double proveBound = atof(argument[3]);
        tas_mnp.SetWeight(proveBound);

        // printf("Solving STP Korf instance [%d of %d] using DSD weight %f\n", atoi(argument[1])+1, 100, atof(argument[3]));

        // Order of calling these functions matters.
        mnp.SetInputWeight(atof(argument[3])); //0

        //If DW mode
        if(atoi(argument[4])==5){

            //Assign a heuristic value, so all states within that range of heuristic value
            //from the start state are going to be weighted.

            mnp.SetMiddleState(start);
            mnp.SetTerrainSize(hardcodedNumbers[atoi(argument[1])]);
            // std::cout<<"SetTerrainSize is: "<<hardcodedNumbers[atoi(argument[1])]<<"\n";
        }
        
        mnp.SetPuzzleWeight(atoi(argument[4])); //1
        // case 0:UnitWeight, case 1:SquareRoot, case 2:Squared,
        // case 3:UnitPlusFrac, case 4:SquarePlusOneRoot, case 5:DW

        mnp.SetMaxMinTileCost(start); //2

        if(atoi(argument[2]) == 0){ //WA*
            tas_mnp.SetPhi([=](double h,double g){return g+proveBound*h;});
        }
        else if(atoi(argument[2]) == 1){ //PWXD
            tas_mnp.SetPhi([=](double h,double g){return (h>g)?(g+h):(g/proveBound+h*(2*proveBound-1)/proveBound);});
        }
        else if(atoi(argument[2]) == 2){ //PWXU
            tas_mnp.SetPhi([=](double h,double g){return (h*(2*proveBound-1)>g)?(g/(2*proveBound-1)+h):(1/proveBound*(g+h));});
        }
        else if(atoi(argument[2]) == 3){ //XDP
            tas_mnp.SetPhi([=](double h,double g){return (g+(2*proveBound-1)*h+sqrt((g-h)*(g-h)+4*proveBound*g*h))/(2*proveBound);});
        }
        else if(atoi(argument[2]) == 4){ //XUP
            tas_mnp.SetPhi([=](double h,double g){return (g+h+sqrt((g+h)*(g+h)+4*proveBound*(proveBound-1)*h*h))/(2*proveBound);});
        }
        else if(atoi(argument[3]) == 5){ //A*
            tas.SetPhi([=](double h,double g){return g+h;});
        }
        
        // tas_mnp.InitializeSearch(&mnp, start, goal, path);
        // tas_mnp.GetPath(&mnp, start, goal, path);
        // printf("STP %d ALG %d weight %1.2f Nodes %llu path %lu\n", atoi(argument[1]), atoi(argument[2]), atof(argument[3]), tas_mnp.GetNodesExpanded(), path.size());

        tas_mnp.InitializeSearch(&mnp, start, goal, path);
        clock_t start_time, end_time;
        start_time = clock();

        float avg_runtime_per_node = tas_mnp.GetPath_v2(&mnp, start, goal, path);
        avg_runtime_per_node *= pow(10, 9);

        end_time = clock();
        float total_runningTime = (float) (end_time - start_time) / CLOCKS_PER_SEC;
        total_runningTime *= pow(10, 9);
            
        printf("Time - STP %d ALG %d weight %1.2f Nodes %llu total_runningTime %1.3f avg_runtime_per_node %1.3f \n", atoi(argument[1]), atoi(argument[2]), atof(argument[3]), tas_mnp.GetNodesExpanded(), total_runningTime, avg_runtime_per_node);
            
        exit(0);
    }
    else if (strcmp(argument[0], "-stpDSD") == 0)
    {
        assert(maxNumArgs >= 5);
        
        DSDWAStar<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>> dsd_mnp;
        std::vector<MNPuzzleState<4, 4>> path;
        MNPuzzle<4, 4> mnp;
        MNPuzzleState<4, 4> start = STP::GetKorfInstance(atoi(argument[1]));
        MNPuzzleState<4, 4> goal;
        //random.randint(10, int(data[10])-10)
        std::vector<int> hardcodedNumbers = {23, 12, 25, 20, 31, 14, 14, 10, 13, 13, 14, 15, 23, 20, 21, 12, 13, 28, 23, 17, 14, 22, 12, 22, 12, 12, 12, 22, 28, 25, 15, 13, 27, 21, 10, 19, 20, 15, 21, 22, 15, 12, 36, 20, 15, 19, 23, 12, 12, 10, 12, 25, 22, 11, 14, 13, 12, 23, 15, 28, 12, 13, 17, 10, 11, 16, 11, 19, 21, 17, 11, 17, 14, 24, 16, 24, 15, 23, 13, 15, 28, 22, 20, 13, 11, 12, 21, 20, 19, 22, 16, 14, 15, 15, 21, 17, 22, 12, 27, 13 };

        dsd_mnp.InitializeSearch(&mnp, start, goal, path);

        dsd_mnp.policy = (tExpansionPriority)atoi(argument[2]);
        dsd_mnp.SetWeight(atof(argument[3]));
        
        // printf("Solving STP Korf instance [%d of %d] using DSD weight %f\n", atoi(argument[1])+1, 100, atof(argument[3]));

        // Order of calling these functions matters.
        mnp.SetInputWeight(atof(argument[3])); //0

        // std::cout<<"Heuristic start to goal in alg "<<argument[2]<<" weight "<<argument[3]<<" is "<<mnp.HCost(start, goal)<<"\n";

        //if DW mode
        if(atoi(argument[4])==5){

            //Assign a heuristic value, so all states within that range of heuristic value
            //from the start state are going to be weighted.

            mnp.SetMiddleState(start);
            mnp.SetTerrainSize(hardcodedNumbers[atoi(argument[1])]);
            // std::cout<<"SetTerrainSize is: "<<hardcodedNumbers[atoi(argument[1])]<<"\n";
        }

        mnp.SetPuzzleWeight(atoi(argument[4])); //1
        // case 0:UnitWeight, case 1:SquareRoot, case 2:Squared,
        // case 3:UnitPlusFrac, case 4:SquarePlusOneRoot, case 5:DW

        mnp.SetMaxMinTileCost(start); //2

        // dsd_mnp.InitializeSearch(&mnp, start, goal, path);
        // dsd_mnp.GetPath(&mnp, start, goal, path);
        // printf("STP %d ALG %d weight %1.2f Nodes %llu path %lu\n", atoi(argument[1]), atoi(argument[2]), atof(argument[3]), dsd_mnp.GetNodesExpanded(), path.size());

        dsd_mnp.InitializeSearch_v3(&mnp, start, goal, path);
        clock_t start_time, end_time;
        start_time = clock();

        float avg_runtime_per_node = dsd_mnp.GetPath_v3(&mnp, start, goal, path);
        avg_runtime_per_node *= pow(10, 9);

        end_time = clock();
        float total_runningTime = (float) (end_time - start_time) / CLOCKS_PER_SEC;
        total_runningTime *= pow(10, 9);
            
        printf("Time - STP %d ALG %d weight %1.2f Nodes %llu total_runningTime %1.3f avg_runtime_per_node %1.3f \n", atoi(argument[1]), atoi(argument[2]), atof(argument[3]), dsd_mnp.GetNodesExpanded(), total_runningTime, avg_runtime_per_node);
        
        
        //Save the SVG of the isoloines plots
        if(saveSVG && (dsd_mnp.policy==5)){
            Graphics::Display d;
            dsd_mnp.DrawPriorityGraph(d);
            std::string s = "stp="+ string(argument[1])+"_w="+string(argument[3])+".svg";
            MakeSVG(d, s.c_str(), 750,750,0);
        }
            
        exit(0);
    }
    else if (strcmp(argument[0], "-stpDPS") == 0)
    {
        assert(maxNumArgs >= 4);
        
        DynamicPotentialSearch<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>> DPS_mnp;
        std::vector<MNPuzzleState<4, 4>> path;
        MNPuzzle<4, 4> mnp;
        MNPuzzleState<4, 4> start = STP::GetKorfInstance(atoi(argument[1]));
        MNPuzzleState<4, 4> goal;
        //random.randint(10, int(data[10])-10)
        std::vector<int> hardcodedNumbers = {23, 12, 25, 20, 31, 14, 14, 10, 13, 13, 14, 15, 23, 20, 21, 12, 13, 28, 23, 17, 14, 22, 12, 22, 12, 12, 12, 22, 28, 25, 15, 13, 27, 21, 10, 19, 20, 15, 21, 22, 15, 12, 36, 20, 15, 19, 23, 12, 12, 10, 12, 25, 22, 11, 14, 13, 12, 23, 15, 28, 12, 13, 17, 10, 11, 16, 11, 19, 21, 17, 11, 17, 14, 24, 16, 24, 15, 23, 13, 15, 28, 22, 20, 13, 11, 12, 21, 20, 19, 22, 16, 14, 15, 15, 21, 17, 22, 12, 27, 13 };

        DPS_mnp.SetOptimalityBound(atof(argument[2]));
        // printf("Solving STP Korf instance [%d of %d] using DSD weight %f\n", atoi(argument[1])+1, 100, atof(argument[2]));

        //if DW mode
        if(atoi(argument[3])==5){
            //Assign a heuristic value, so all states within that range of heuristic value
            //from the start state are going to be weighted.

            mnp.SetMiddleState(start);
            mnp.SetTerrainSize(hardcodedNumbers[atoi(argument[1])]);
            // std::cout<<"SetTerrainSize is: "<<hardcodedNumbers[atoi(argument[1])]<<"\n";
        }

        mnp.SetPuzzleWeight(atoi(argument[3])); //1
        // case 0:UnitWeight, case 1:SquareRoot, case 2:Squared,
        // case 3:UnitPlusFrac, case 4:SquarePlusOneRoot, case 5:DW

        // DPS_mnp.InitializeSearch(&mnp, start, goal, path);
        // DPS_mnp.GetPath(&mnp, start, goal, path);
        // printf("STP %d ALG %d weight %1.2f Nodes %llu path %lu\n", atoi(argument[1]), 6, atof(argument[2]), DPS_mnp.GetNodesExpanded(), path.size());
        
        DPS_mnp.InitializeSearch(&mnp, start, goal, path);
        clock_t start_time, end_time;
        start_time = clock();

        float avg_runtime_per_node = DPS_mnp.GetPath_v2(&mnp, start, goal, path);
        avg_runtime_per_node *= pow(10, 9);

        end_time = clock();
        float total_runningTime = (float) (end_time - start_time) / CLOCKS_PER_SEC;
        total_runningTime *= pow(10, 9);
            
        printf("Time - STP %d ALG %d weight %1.2f Nodes %llu total_runningTime %1.3f avg_runtime_per_node %1.3f \n", atoi(argument[1]), 6, atoi(argument[2]), DPS_mnp.GetNodesExpanded(), total_runningTime, avg_runtime_per_node);
        

        exit(0);
    }
    else if (strcmp(argument[0], "-gridBaseLines") == 0){
        // Runs all five baselines.
        assert(maxNumArgs >= 6);
        me = new MapEnvironment(new Map(argument[1]));
        exper = atoi(argument[7]);

        //Set the heuristic
        if(useDH){
            dh = new GridEmbedding(me, 10, kLINF);
            for (int x = 0; x < 10; x++)
                dh->AddDimension(kDifferential, kFurthest);
            tas.SetHeuristic(dh);
        }

        // me->SetDiagonalCost(1.5);
        me->SetDiagonalCost(1.41);
        swampedloc = {1,1};

        // ScenarioLoader sl(argument[2]);
        sl = new ScenarioLoader(argument[2]);
        int approvedScenaios = 0;
        
        for (int x = 0; x < sl->GetNumExperiments(); x++)
        {
            // Experiment exp = sl.GetNthExperiment(x);
            // if(exp.GetDistance()<30) continue;
            if(approvedScenaios >= numLimitedScenarios && limitScenarios) break;
            Experiment exp = sl->GetNthExperiment(x);
            if(limitScenarios && (exp.GetDistance()<lowerLimit || exp.GetDistance()>upperLimit)) continue;
            else approvedScenaios ++;
            
            //Reset the last problem's swamped states.
            if(exper==8){
                if(swampedloc.x == 0){
                    for(int i=0; i < me->GetMap()->GetMapWidth(); i++)
                        for(int j=swampedloc.y-int(terrain_height/2); j<=swampedloc.y+int(terrain_height/2); j++)
                            if(me->GetMap()->GetTerrainType(i, j) == kSwamp)
                                me->GetMap()->SetTerrainType(i, j, kGround);
                }
                else if(swampedloc.y == 0){
                    for(int j=0; j < me->GetMap()->GetMapHeight(); j++)
                        for(int i=swampedloc.x-int(terrain_width/2); i<=swampedloc.x+int(terrain_width/2); i++)
                            if(me->GetMap()->GetTerrainType(i, j) == kSwamp)
                                me->GetMap()->SetTerrainType(i, j, kGround);
                }
            }

            //Set the start and goal states, weight and policy of the search.
            start.x = exp.GetStartX();
            start.y = exp.GetStartY();
            goal.x = exp.GetGoalX();
            goal.y = exp.GetGoalY();
            double proveBound = atof(argument[4]);
            tas.SetWeight(proveBound);
            me->SetInputWeight(atof(argument[4]));
            
            //Set the cost of each terrain type randomly.
            for(int i=0; i<4; i++){
                //[0]=kSwamp, [1]=kWater, [2]=kGrass, [3]=kTrees

                // Define the Hardness
                // 1. Random Hardness
                // rdm = random()%101;
                // hardness[i] = rdm/101+1;
                //Or 2. Hardcoded Hardness
                // hardness[0]=1.1; hardness[1]=1.45; hardness[2]=1.25; hardness[3]=1.95;
                //Or 3. Get Hardness from Input Argument
                hardness[0]=atof(argument[6]);/*Followings are "Don't Care":*/hardness[1]=100.0; hardness[2]=100.0; hardness[3]=100.0;

                // Use Hardness to define Cost
                // 1. Hardness with respect to input w
                // Tcosts[i] = hardness[i]*(me->GetInputWeight())-(hardness[i]-1);
                //Or 2. Hardness as the exact Cost
                Tcosts[i] = hardness[i];
                //Or 3. The exact Cost hardcoded
                // Tcosts[i] = 4.5;

                string type;
                if(i==0) type="Swamp";
                if(i==1) type="Water";
                if(i==2) type="Grass";
                if(i==3) type="Trees";
                // std::cout<<"Cost "<<type<<"="<<Tcosts[i]<<" ("<<hardness[i]<<"*"<<me->GetInputWeight()<<"-"<<(hardness[i]-1)<<")"<<std::endl;
            }
            me->SetTerrainCost(Tcosts);
            
            //Change one or more squares of ground states to swamp type.
            if(exper==8){
                //Set the size of the swamp area using the terrain_size argument.
                terrain_width = max(terrain_size/100*abs(goal.x-start.x), 10);
                terrain_height = max(terrain_size/100*abs(goal.y-start.y), 10);

                if(abs(goal.x-start.x) > abs(goal.y-start.y) && abs(goal.x-start.x)!=0){
                    // swampedloc.x = (goal.x+start.x)/2;
                    swampedloc.x = random()%(abs(goal.x-start.x))+min(goal.x,start.x);
                    for(int j=0; j < me->GetMap()->GetMapHeight(); j++)
                        for(int i=swampedloc.x-int(terrain_width/2); i<=swampedloc.x+int(terrain_width/2); i++)
                            if(me->GetMap()->GetTerrainType(i, j) == kGround)
                                me->GetMap()->SetTerrainType(i, j, kSwamp);
                    swampedloc.y = 0;
                }
                else if(abs(goal.y-start.y) > abs(goal.x-start.x) && abs(goal.y-start.y)!=0){
                    // swampedloc.y = (goal.y+start.y)/2;
                    swampedloc.y = random()%(abs(goal.y-start.y))+min(goal.y,start.y);
                    for(int i=0; i < me->GetMap()->GetMapWidth(); i++)
                        for(int j=swampedloc.y-int(terrain_height/2); j<=swampedloc.y+int(terrain_height/2); j++)
                            if(me->GetMap()->GetTerrainType(i, j) == kGround)
                                me->GetMap()->SetTerrainType(i, j, kSwamp);
                    swampedloc.x = 0;
                }
                else if(abs(goal.x-start.x) > abs(goal.y-start.y)){
                    swampedloc.x = (goal.x+start.x)/2;
                    // swampedloc.x = random()%(abs(goal.x-start.x))+min(goal.x,start.x);
                    for(int j=0; j < me->GetMap()->GetMapHeight(); j++)
                        for(int i=swampedloc.x-int(terrain_width/2); i<=swampedloc.x+int(terrain_width/2); i++)
                            if(me->GetMap()->GetTerrainType(i, j) == kGround)
                                me->GetMap()->SetTerrainType(i, j, kSwamp);
                    swampedloc.y = 0;
                }
                else{
                    swampedloc.y = (goal.y+start.y)/2;
                    // swampedloc.y = random()%(abs(goal.y-start.y))+min(goal.y,start.y);
                    for(int i=0; i < me->GetMap()->GetMapWidth(); i++)
                        for(int j=swampedloc.y-int(terrain_height/2); j<=swampedloc.y+int(terrain_height/2); j++)
                            if(me->GetMap()->GetTerrainType(i, j) == kGround)
                                me->GetMap()->SetTerrainType(i, j, kSwamp);
                    swampedloc.x = 0;
                }
            }

            //Set the Baseline Policy.
            if(atoi(argument[3]) == 0){ //WA*
            tas.SetPhi([=](double h,double g){return g+proveBound*h;});
            }
            else if(atoi(argument[3]) == 1){ //PWXD
            tas.SetPhi([=](double h,double g){return (h>g)?(g+h):(g/proveBound+h*(2*proveBound-1)/proveBound);});
            }
            else if(atoi(argument[3]) == 2){ //PWXU
            tas.SetPhi([=](double h,double g){return (h*(2*proveBound-1)>g)?(g/(2*proveBound-1)+h):(1/proveBound*(g+h));});
            }
            else if(atoi(argument[3]) == 3){ //XDP
            tas.SetPhi([=](double h,double g){return (g+(2*proveBound-1)*h+sqrt((g-h)*(g-h)+4*proveBound*g*h))/(2*proveBound);});
            }
            else if(atoi(argument[3]) == 4){ //XUP
            tas.SetPhi([=](double h,double g){return (g+h+sqrt((g+h)*(g+h)+4*proveBound*(proveBound-1)*h*h))/(2*proveBound);});
            }
            else if(atoi(argument[3]) == 5){ //A*
            tas.SetPhi([=](double h,double g){return g+h;});
            }
            
            // tas.InitializeSearch(me, start, goal, solution);
            // tas.GetPath(me, start, goal, solution);
            // printf("MAP %s #%d %1.2f ALG %d weight %1.2f Nodes %llu path %f\n", argument[1], x, exp.GetDistance(), atoi(argument[3]), atof(argument[4]), tas.GetNodesExpanded(), me->GetPathLength(solution));

            tas.InitializeSearch(me, start, goal, solution);

            float avg_runtime_per_node = tas.GetPath_v2(me, start, goal, solution);
            avg_runtime_per_node *= pow(10, 9);
  
            printf("Time - MAP %s #%d %1.2f ALG %d weight %1.2f Nodes %llu avg_runtime_per_node %1.3f\n", argument[1], x, exp.GetDistance(), atoi(argument[3]), atof(argument[4]), tas.GetNodesExpanded(), avg_runtime_per_node);

        }
        exit(0);
    }
    else if (strcmp(argument[0], "-gridDSD") == 0){
        // Arguments: -DSMAP mapAddress $scen $alg $weight $TerrainSize $SwampHardness

        //Thie is Prior Knowledge Map.
        //The knowledge we have prior to the search about the domain, is it's a dynamic weighted grid map.
        //Some Terrain Types exist that makes action costs different in different parts of the map.
        assert(maxNumArgs >= 8);
        me = new MapEnvironment(new Map(argument[1]));
        exper = atoi(argument[7]);

        //Set the heuristic
        if(useDH){
            dh = new GridEmbedding(me, 10, kLINF);
            for (int x = 0; x < 10; x++)
                dh->AddDimension(kDifferential, kFurthest);
            dsd.SetHeuristic(dh);
        }

        // me->SetDiagonalCost(1.5);
        me->SetDiagonalCost(1.41);
        swampedloc = {1,1};

        // ScenarioLoader sl(argument[2]);
        sl = new ScenarioLoader(argument[2]);
        int approvedScenaios = 0;
        
        for (int x = 0; x < sl->GetNumExperiments(); x++)
        {
            // Experiment exp = sl.GetNthExperiment(x);
            // if(exp.GetDistance()<30) continue;
            if(approvedScenaios >= numLimitedScenarios && limitScenarios) break;
            Experiment exp = sl->GetNthExperiment(x);
            if(limitScenarios && (exp.GetDistance()<lowerLimit || exp.GetDistance()>upperLimit)) continue;
            else approvedScenaios ++;

            //Reset the last problem's swamped states.
            if(exper==8){
                if(swampedloc.x == 0){
                    for(int i=0; i < me->GetMap()->GetMapWidth(); i++)
                        for(int j=swampedloc.y-int(terrain_height/2); j<=swampedloc.y+int(terrain_height/2); j++)
                            if(me->GetMap()->GetTerrainType(i, j) == kSwamp)
                                me->GetMap()->SetTerrainType(i, j, kGround);
                }
                else if(swampedloc.y == 0){
                    for(int j=0; j < me->GetMap()->GetMapHeight(); j++)
                        for(int i=swampedloc.x-int(terrain_width/2); i<=swampedloc.x+int(terrain_width/2); i++)
                            if(me->GetMap()->GetTerrainType(i, j) == kSwamp)
                                me->GetMap()->SetTerrainType(i, j, kGround);
                }
            }

            //Set the start and goal states, weight and policy of the search.
            start.x = exp.GetStartX();
            start.y = exp.GetStartY();
            goal.x = exp.GetGoalX();
            goal.y = exp.GetGoalY();
            dsd.policy = (tExpansionPriority)atoi(argument[3]);
            dsd.SetWeight(atof(argument[4]));
            me->SetInputWeight(atof(argument[4]));

            //Set the cost of each terrain type randomly.
            for(int i=0; i<4; i++){
                //[0]=kSwamp, [1]=kWater, [2]=kGrass, [3]=kTrees

                //Define the Hardness
                //1. Random Hardness
                // rdm = random()%101;
                // hardness[i] = rdm/101+1;
                //Or 2. Hardcoded Hardness
                // hardness[0]=1.1; hardness[1]=1.45; hardness[2]=1.25; hardness[3]=1.95;
                //Or 3. Get Hardness from Input Argument
                hardness[0]=atof(argument[6]);/*Followings are "Don't Care":*/hardness[1]=100.0; hardness[2]=100.0; hardness[3]=100.0;

                //Use Hardness to define Cost
                // 1. Hardness with respect to input w
                // Tcosts[i] = hardness[i]*(me->GetInputWeight())-(hardness[i]-1);
                //Or 2. Hardness as the exact Cost
                Tcosts[i] = hardness[i];
                //Or 3. The exact Cost hardcoded
                // Tcosts[i] = 4.5;

                string type;
                if(i==0) type="Swamp";
                if(i==1) type="Water";
                if(i==2) type="Grass";
                if(i==3) type="Trees";
                // std::cout<<"Cost "<<type<<"="<<Tcosts[i]<<" ("<<hardness[i]<<"*"<<me->GetInputWeight()<<"-"<<(hardness[i]-1)<<")"<<std::endl;
            }
            me->SetTerrainCost(Tcosts);
            
            //Change one or more squares of ground states to swamp type.
            if(exper==8){
                //Set the size of the swamp area using the terrain_size argument.
                terrain_width = max(terrain_size/100*abs(goal.x-start.x), 10);
                terrain_height = max(terrain_size/100*abs(goal.y-start.y), 10);

                if(abs(goal.x-start.x) > abs(goal.y-start.y) && abs(goal.x-start.x)!=0){
                    // swampedloc.x = (goal.x+start.x)/2;
                    swampedloc.x = random()%(abs(goal.x-start.x))+min(goal.x,start.x);
                    for(int j=0; j < me->GetMap()->GetMapHeight(); j++)
                        for(int i=swampedloc.x-int(terrain_width/2); i<=swampedloc.x+int(terrain_width/2); i++)
                            if(me->GetMap()->GetTerrainType(i, j) == kGround)
                                me->GetMap()->SetTerrainType(i, j, kSwamp);
                    swampedloc.y = 0;
                }
                else if(abs(goal.y-start.y) > abs(goal.x-start.x) && abs(goal.y-start.y)!=0){
                    // swampedloc.y = (goal.y+start.y)/2;
                    swampedloc.y = random()%(abs(goal.y-start.y))+min(goal.y,start.y);
                    for(int i=0; i < me->GetMap()->GetMapWidth(); i++)
                        for(int j=swampedloc.y-int(terrain_height/2); j<=swampedloc.y+int(terrain_height/2); j++)
                            if(me->GetMap()->GetTerrainType(i, j) == kGround)
                                me->GetMap()->SetTerrainType(i, j, kSwamp);
                    swampedloc.x = 0;
                }
                else if(abs(goal.x-start.x) > abs(goal.y-start.y)){
                    swampedloc.x = (goal.x+start.x)/2;
                    // swampedloc.x = random()%(abs(goal.x-start.x))+min(goal.x,start.x);
                    for(int j=0; j < me->GetMap()->GetMapHeight(); j++)
                        for(int i=swampedloc.x-int(terrain_width/2); i<=swampedloc.x+int(terrain_width/2); i++)
                            if(me->GetMap()->GetTerrainType(i, j) == kGround)
                                me->GetMap()->SetTerrainType(i, j, kSwamp);
                    swampedloc.y = 0;
                }
                else{
                    swampedloc.y = (goal.y+start.y)/2;
                    // swampedloc.y = random()%(abs(goal.y-start.y))+min(goal.y,start.y);
                    for(int i=0; i < me->GetMap()->GetMapWidth(); i++)
                        for(int j=swampedloc.y-int(terrain_height/2); j<=swampedloc.y+int(terrain_height/2); j++)
                            if(me->GetMap()->GetTerrainType(i, j) == kGround)
                                me->GetMap()->SetTerrainType(i, j, kSwamp);
                    swampedloc.x = 0;
                }
            }
            
            // dsd.InitializeSearch(me, start, goal, solution);
            // dsd.GetPath(me, start, goal, solution);
            // printf("MAP %s #%d %1.2f ALG %d weight %1.2f Nodes %llu path %f\n", argument[1], x, exp.GetDistance(), atoi(argument[3]), atof(argument[4]), dsd.GetNodesExpanded(), me->GetPathLength(solution));

            dsd.InitializeSearch_v3(me, start, goal, solution);
            
            float avg_runtime_per_node = dsd.GetPath_v3(me, start, goal, solution);
            avg_runtime_per_node *= pow(10, 9);

            printf("Time - MAP %s #%d %1.2f ALG %d weight %1.2f Nodes %llu avg_runtime_per_node %1.3f\n", argument[1], x, exp.GetDistance(), atoi(argument[3]), atof(argument[4]), dsd.GetNodesExpanded(), avg_runtime_per_node);

        }
        exit(0);
    }
    else if (strcmp(argument[0], "-gridDPS") == 0){
        assert(maxNumArgs >= 7);
        me = new MapEnvironment(new Map(argument[1]));
        exper = atoi(argument[6]);

        //Set the heuristic
        if(useDH){
            dh = new GridEmbedding(me, 10, kLINF);
            for (int x = 0; x < 10; x++)
                dh->AddDimension(kDifferential, kFurthest);
            dsd.SetHeuristic(dh);
        }

        // me->SetDiagonalCost(1.5);
        me->SetDiagonalCost(1.41);
        swampedloc = {1,1};

        // ScenarioLoader sl(argument[2]);
        sl = new ScenarioLoader(argument[2]);
        int approvedScenaios = 0;
        
        for (int x = 0; x < sl->GetNumExperiments(); x++)
        {
            // Experiment exp = sl.GetNthExperiment(x);
            // if(exp.GetDistance()<30) continue;
            if(approvedScenaios >= numLimitedScenarios && limitScenarios) break;
            Experiment exp = sl->GetNthExperiment(x);
            if(limitScenarios && (exp.GetDistance()<lowerLimit || exp.GetDistance()>upperLimit)) continue;
            else approvedScenaios ++;

            //Reset the last problem's swamped states.
            if(exper==8){
                if(swampedloc.x == 0){
                    for(int i=0; i < me->GetMap()->GetMapWidth(); i++)
                        for(int j=swampedloc.y-int(terrain_height/2); j<=swampedloc.y+int(terrain_height/2); j++)
                            if(me->GetMap()->GetTerrainType(i, j) == kSwamp)
                                me->GetMap()->SetTerrainType(i, j, kGround);
                }
                else if(swampedloc.y == 0){
                    for(int j=0; j < me->GetMap()->GetMapHeight(); j++)
                        for(int i=swampedloc.x-int(terrain_width/2); i<=swampedloc.x+int(terrain_width/2); i++)
                            if(me->GetMap()->GetTerrainType(i, j) == kSwamp)
                                me->GetMap()->SetTerrainType(i, j, kGround);
                }
            }

            //Set the start and goal states, weight and policy of the search.
            start.x = exp.GetStartX();
            start.y = exp.GetStartY();
            goal.x = exp.GetGoalX();
            goal.y = exp.GetGoalY();
                
            dps.SetOptimalityBound(atof(argument[3]));
            me->SetInputWeight(atof(argument[3]));

            //Set the cost of each terrain type randomly.
            for(int i=0; i<4; i++){
                //[0]=kSwamp, [1]=kWater, [2]=kGrass, [3]=kTrees

                //Define the Hardness
                //1. Random Hardness
                // rdm = random()%101;
                // hardness[i] = rdm/101+1;
                //Or 2. Hardcoded Hardness
                // hardness[0]=1.1; hardness[1]=1.45; hardness[2]=1.25; hardness[3]=1.95;
                //Or 3. Get Hardness from Input Argument
                hardness[0]=atof(argument[5]);/*Followings are "Don't Care":*/hardness[1]=100.0; hardness[2]=100.0; hardness[3]=100.0;

                //Use Hardness to define Cost
                // 1. Hardness with respect to input w
                // Tcosts[i] = hardness[i]*(me->GetInputWeight())-(hardness[i]-1);
                //Or 2. Hardness as the exact Cost
                Tcosts[i] = hardness[i];
                //Or 3. The exact Cost hardcoded
                // Tcosts[i] = 4.5;

                string type;
                if(i==0) type="Swamp";
                if(i==1) type="Water";
                if(i==2) type="Grass";
                if(i==3) type="Trees";
                // std::cout<<"Cost "<<type<<"="<<Tcosts[i]<<" ("<<hardness[i]<<"*"<<me->GetInputWeight()<<"-"<<(hardness[i]-1)<<")"<<std::endl;
            }
            me->SetTerrainCost(Tcosts);
            
            //Change one or more squares of ground states to swamp type.
            if(exper==8){
                //Set the size of the swamp area using the terrain_size argument.
                terrain_width = max(terrain_size/100*abs(goal.x-start.x), 10);
                terrain_height = max(terrain_size/100*abs(goal.y-start.y), 10);

                if(abs(goal.x-start.x) > abs(goal.y-start.y) && abs(goal.x-start.x)!=0){
                    // swampedloc.x = (goal.x+start.x)/2;
                    swampedloc.x = random()%(abs(goal.x-start.x))+min(goal.x,start.x);
                    for(int j=0; j < me->GetMap()->GetMapHeight(); j++)
                        for(int i=swampedloc.x-int(terrain_width/2); i<=swampedloc.x+int(terrain_width/2); i++)
                            if(me->GetMap()->GetTerrainType(i, j) == kGround)
                                me->GetMap()->SetTerrainType(i, j, kSwamp);
                    swampedloc.y = 0;
                }
                else if(abs(goal.y-start.y) > abs(goal.x-start.x) && abs(goal.y-start.y)!=0){
                    // swampedloc.y = (goal.y+start.y)/2;
                    swampedloc.y = random()%(abs(goal.y-start.y))+min(goal.y,start.y);
                    for(int i=0; i < me->GetMap()->GetMapWidth(); i++)
                        for(int j=swampedloc.y-int(terrain_height/2); j<=swampedloc.y+int(terrain_height/2); j++)
                            if(me->GetMap()->GetTerrainType(i, j) == kGround)
                                me->GetMap()->SetTerrainType(i, j, kSwamp);
                    swampedloc.x = 0;
                }
                else if(abs(goal.x-start.x) > abs(goal.y-start.y)){
                    swampedloc.x = (goal.x+start.x)/2;
                    // swampedloc.x = random()%(abs(goal.x-start.x))+min(goal.x,start.x);
                    for(int j=0; j < me->GetMap()->GetMapHeight(); j++)
                        for(int i=swampedloc.x-int(terrain_width/2); i<=swampedloc.x+int(terrain_width/2); i++)
                            if(me->GetMap()->GetTerrainType(i, j) == kGround)
                                me->GetMap()->SetTerrainType(i, j, kSwamp);
                    swampedloc.y = 0;
                }
                else{
                    swampedloc.y = (goal.y+start.y)/2;
                    // swampedloc.y = random()%(abs(goal.y-start.y))+min(goal.y,start.y);
                    for(int i=0; i < me->GetMap()->GetMapWidth(); i++)
                        for(int j=swampedloc.y-int(terrain_height/2); j<=swampedloc.y+int(terrain_height/2); j++)
                            if(me->GetMap()->GetTerrainType(i, j) == kGround)
                                me->GetMap()->SetTerrainType(i, j, kSwamp);
                    swampedloc.x = 0;
                }
            }
            
            // dps.InitializeSearch(me, start, goal, solution);
            // dps.GetPath(me, start, goal, solution);
            // printf("MAP %s #%d %1.2f ALG %d weight %1.2f Nodes %llu path %f\n", argument[1], x, exp.GetDistance(), 6, atof(argument[3]), dps.GetNodesExpanded(), me->GetPathLength(solution));

            dps.InitializeSearch(me, start, goal, solution);

            float avg_runtime_per_node = dps.GetPath_v2(me, start, goal, solution);
            avg_runtime_per_node *= pow(10, 9);
                
            printf("Time - MAP %s #%d %1.2f ALG %d weight %1.2f Nodes %llu avg_runtime_per_node %1.3f\n", argument[1], x, exp.GetDistance(), 6, atof(argument[3]), dps.GetNodesExpanded(), avg_runtime_per_node);

        }
        exit(0);
    }
    else if (strcmp(argument[0], "-rtBaseLines") == 0){
        // Runs all five baselines.
        assert(maxNumArgs >= 6);
        Map * m = new Map(argument[1]);
        r = new Racetrack(m);
        me = new MapEnvironment(m);

        // ScenarioLoader sl(argument[2]);
        sl = new ScenarioLoader(argument[2]);
        int approvedScenaios = 0, x=0;
        while ((approvedScenaios < numLimitedScenarios) && (x < sl->GetNumExperiments()))
        {
            Experiment exp = sl->GetNthExperiment(x);
            x=x+(sl->GetNumExperiments()/1000)+1;
            if((exp.GetDistance()<lowerLimit) || (exp.GetDistance()>upperLimit)) continue;
            else approvedScenaios ++;

            //Reset the previous start and goal
            m->SetTerrainType(from.xLoc, from.yLoc, kGround);
            if(numExtendedGoals == 0){//Not extended the goal.
                m->SetTerrainType(end.xLoc, end.yLoc, kGround);
            }
            else{//Extended the goal
                for(int tNode=0; tNode<theList.size(); tNode++){
                    m->SetTerrainType(theList[tNode].x, theList[tNode].y, kGround);
                }
            }

            //Set the start and goal states, weight and policy of the search.
            from.xLoc = exp.GetStartX();
            from.yLoc = exp.GetStartY();
            from.xVelocity = 0;
            from.yVelocity = 0;
            end.xLoc = exp.GetGoalX();
            end.yLoc = exp.GetGoalY();
            end.xVelocity = 0;
            end.yVelocity = 0;

            //Set the start and goal TerrainType.
            numExtendedGoals = atoi(argument[5]);
            if(numExtendedGoals) numExtendedGoals = exp.GetDistance()/3;
            m->SetTerrainType(from.xLoc, from.yLoc, kStartTerrain);
            if(numExtendedGoals == 0){//Not extend the goal.
                m->SetTerrainType(end.xLoc, end.yLoc, kEndTerrain);
            }
            else{//Extend the goal
                tas.SetPhi([=](double h,double g){return g;});
                goal.x = end.xLoc;
                goal.y = end.yLoc;
                tas.ExtendGoal(me, goal, theList, numExtendedGoals);
                for(int tNode=0; tNode<theList.size(); tNode++){
                    m->SetTerrainType(theList[tNode].x, theList[tNode].y, kEndTerrain);
                }
            }
            
            r->UpdateMap(m);

            double proveBound = atof(argument[4]);
            tas_track.SetWeight(proveBound);

            //Set the Baseline Policy.
            if(atoi(argument[3]) == 0){ //WA*
            tas_track.SetPhi([=](double h,double g){return g+proveBound*h;});
            }
            else if(atoi(argument[3]) == 1){ //PWXD
            tas_track.SetPhi([=](double h,double g){return (h>g)?(g+h):(g/proveBound+h*(2*proveBound-1)/proveBound);});
            }
            else if(atoi(argument[3]) == 2){ //PWXU
            tas_track.SetPhi([=](double h,double g){return (h*(2*proveBound-1)>g)?(g/(2*proveBound-1)+h):(1/proveBound*(g+h));});
            }
            else if(atoi(argument[3]) == 3){ //XDP
            tas_track.SetPhi([=](double h,double g){return (g+(2*proveBound-1)*h+sqrt((g-h)*(g-h)+4*proveBound*g*h))/(2*proveBound);});
            }
            else if(atoi(argument[3]) == 4){ //XUP
            tas_track.SetPhi([=](double h,double g){return (g+h+sqrt((g+h)*(g+h)+4*proveBound*(proveBound-1)*h*h))/(2*proveBound);});
            }

            // tas_track.InitializeSearch(r, from, end, path);
            // tas_track.GetPath(r, from, end, path);
            // printf("MAP %s #%d %1.2f ALG %d weight %1.2f Nodes %llu \n", argument[1], x, exp.GetDistance(), atoi(argument[3]), atof(argument[4]), tas_track.GetNodesExpanded());

            tas_track.InitializeSearch(r, from, end, path);
            clock_t start_time, end_time;
            start_time = clock();

            float avg_runtime_per_node = tas_track.GetPath_v2(r, from, end, path);
            avg_runtime_per_node *= pow(10, 9);

            end_time = clock();
            float total_runningTime = (float) (end_time - start_time) / CLOCKS_PER_SEC;
            total_runningTime *= pow(10, 9);
                
            printf("Time - MAP %s #%d %1.2f ALG %d weight %1.2f Nodes %llu total_runningTime %1.3f avg_runtime_per_node %1.3f\n", argument[1], x, exp.GetDistance(), atoi(argument[3]), atof(argument[4]), tas_track.GetNodesExpanded(), total_runningTime, avg_runtime_per_node);


        }
        exit(0);
    }
    else if (strcmp(argument[0], "-rtDSD") == 0){
        assert(maxNumArgs >= 6);
        Map * m = new Map(argument[1]);
        r = new Racetrack(m);
        me = new MapEnvironment(m);

        // ScenarioLoader sl(argument[2]);
        sl = new ScenarioLoader(argument[2]);
        int approvedScenaios = 0, x=0;
        // numLimitedScenarios = max(sl->GetNumExperiments(), numLimitedScenarios);
        
        while ((approvedScenaios < numLimitedScenarios) && (x < sl->GetNumExperiments()))
        {
            Experiment exp = sl->GetNthExperiment(x);
            x=x+(sl->GetNumExperiments()/1000)+1;
            if((exp.GetDistance()<lowerLimit) || (exp.GetDistance()>upperLimit)) continue;
            else approvedScenaios ++;

            //Reset the previous start and goal
            m->SetTerrainType(from.xLoc, from.yLoc, kGround);
            if(numExtendedGoals == 0){//Not extended the goal.
                m->SetTerrainType(end.xLoc, end.yLoc, kGround);
            }
            else{//Extended the goal
                for(int tNode=0; tNode<theList.size(); tNode++){
                    m->SetTerrainType(theList[tNode].x, theList[tNode].y, kGround);
                }
            }

            //Set the start and goal states, weight and policy of the search.
            from.xLoc = exp.GetStartX();
            from.yLoc = exp.GetStartY();
            from.xVelocity = 0;
            from.yVelocity = 0;
            end.xLoc = exp.GetGoalX();
            end.yLoc = exp.GetGoalY();
            end.xVelocity = 0;
            end.yVelocity = 0;
            dsd_track.policy = (tExpansionPriority)atoi(argument[3]);
            dsd_track.SetWeight(atof(argument[4]));
            
            //Set the start and goal TerrainType.
            numExtendedGoals = atoi(argument[5]);
            if(numExtendedGoals) numExtendedGoals = exp.GetDistance()/3;
            m->SetTerrainType(from.xLoc, from.yLoc, kStartTerrain);
            if(numExtendedGoals == 0){//Not extend the goal.
                m->SetTerrainType(end.xLoc, end.yLoc, kEndTerrain);
            }
            else{//Extend the goal
                tas.SetPhi([=](double h,double g){return g;});
                goal.x = end.xLoc;
                goal.y = end.yLoc;
                tas.ExtendGoal(me, goal, theList, numExtendedGoals);
                for(int tNode=0; tNode<theList.size(); tNode++){
                    m->SetTerrainType(theList[tNode].x, theList[tNode].y, kEndTerrain);
                }
            }

            r->UpdateMap(m);

            // dsd_track.InitializeSearch(r, from, end, path);
            // dsd_track.GetPath(r, from, end, path);
            // printf("MAP %s #%d %1.2f ALG %d weight %1.2f Nodes %llu \n", argument[1], x, exp.GetDistance(), atoi(argument[3]), atof(argument[4]), dsd_track.GetNodesExpanded());

            dsd_track.InitializeSearch_v3(r, from, end, path);
            clock_t start_time, end_time;
            start_time = clock();

            float avg_runtime_per_node = dsd_track.GetPath_v3(r, from, end, path);
            avg_runtime_per_node *= pow(10, 9);

            end_time = clock();
            float total_runningTime = (float) (end_time - start_time) / CLOCKS_PER_SEC;
            total_runningTime *= pow(10, 9);
                
            printf("Time - MAP %s #%d %1.2f ALG %d weight %1.2f Nodes %llu total_runningTime %1.3f avg_runtime_per_node %1.3f\n", argument[1], x, exp.GetDistance(), atoi(argument[3]), atof(argument[4]), dsd_track.GetNodesExpanded(), total_runningTime, avg_runtime_per_node);


        }
        exit(0);
    }
    else if (strcmp(argument[0], "-rtDPS") == 0)
    {
        assert(maxNumArgs >= 5);
        Map * m = new Map(argument[1]);
        r = new Racetrack(m);
        me = new MapEnvironment(m);

        // ScenarioLoader sl(argument[2]);
        sl = new ScenarioLoader(argument[2]);
        int approvedScenaios = 0, x=0;

        while ((approvedScenaios < numLimitedScenarios) && (x < sl->GetNumExperiments()))
        {
            Experiment exp = sl->GetNthExperiment(x);
            x=x+(sl->GetNumExperiments()/1000)+1;
            if((exp.GetDistance()<lowerLimit) || (exp.GetDistance()>upperLimit)) continue;
            else approvedScenaios ++;

            //Reset the previous start and goal
            m->SetTerrainType(from.xLoc, from.yLoc, kGround);
            if(numExtendedGoals == 0){//Not extended the goal.
                m->SetTerrainType(end.xLoc, end.yLoc, kGround);
            }
            else{//Extended the goal
                for(int tNode=0; tNode<theList.size(); tNode++){
                    m->SetTerrainType(theList[tNode].x, theList[tNode].y, kGround);
                }
            }

            //Set the start and goal states, weight and policy of the search.
            from.xLoc = exp.GetStartX();
            from.yLoc = exp.GetStartY();
            from.xVelocity = 0;
            from.yVelocity = 0;
            end.xLoc = exp.GetGoalX();
            end.yLoc = exp.GetGoalY();
            end.xVelocity = 0;
            end.yVelocity = 0;

            //Set the start and goal TerrainType.
            numExtendedGoals = atoi(argument[4]);
            if(numExtendedGoals) numExtendedGoals = exp.GetDistance()/3;
            m->SetTerrainType(from.xLoc, from.yLoc, kStartTerrain);
            if(numExtendedGoals == 0){//Not extend the goal.
                m->SetTerrainType(end.xLoc, end.yLoc, kEndTerrain);
            }
            else{//Extend the goal
                tas.SetPhi([=](double h,double g){return g;});
                goal.x = end.xLoc;
                goal.y = end.yLoc;
                tas.ExtendGoal(me, goal, theList, numExtendedGoals);
                for(int tNode=0; tNode<theList.size(); tNode++){
                    m->SetTerrainType(theList[tNode].x, theList[tNode].y, kEndTerrain);
                }
            }
            
            r->UpdateMap(m);

            double proveBound = atof(argument[3]);
            dps_track.SetOptimalityBound(proveBound);

            // dps_track.InitializeSearch(r, from, end, path);
            // dps_track.GetPath(r, from, end, path);
            // printf("MAP %s #%d %1.2f ALG %d weight %1.2f Nodes %llu \n", argument[1], x, exp.GetDistance(), 6, atof(argument[3]), dps_track.GetNodesExpanded());

            dps_track.InitializeSearch(r, from, end, path);
            clock_t start_time, end_time;
            start_time = clock();

            float avg_runtime_per_node = dps_track.GetPath_v2(r, from, end, path);
            avg_runtime_per_node *= pow(10, 9);

            end_time = clock();
            float total_runningTime = (float) (end_time - start_time) / CLOCKS_PER_SEC;
            total_runningTime *= pow(10, 9);
                
            printf("Time - MAP %s #%d %1.2f ALG %d weight %1.2f Nodes %llu total_runningTime %1.3f avg_runtime_per_node %1.3f\n", argument[1], x, exp.GetDistance(), 6, atof(argument[3]), dps_track.GetNodesExpanded(), total_runningTime, avg_runtime_per_node);


        }
        
        exit(0);
    }
    else if (strcmp(argument[0], "-exp0BaseLines") == 0)
    {
        Map *m = new Map(200, 200);
        MakeDesignedMap(m, atoi(argument[3]), atoi(argument[4]));
        // MakeDesignedMap(m, 30, 6);

        me = new MapEnvironment(m);
        // me->SetDiagonalCost(1.5);
        me->SetDiagonalCost(1.41);

        assert(maxNumArgs >= 5);
        //Set the cost of each terrain type randomly.
        for(int i=0; i<4; i++){
            //[0]=kSwamp, [1]=kWater, [2]=kGrass, [3]=kTrees

            // Define the Hardness
            // 1. Random Hardness
            // rdm = random()%101;
            // hardness[i] = rdm/101+1;
            //Or 2. Hardcoded Hardness
            // hardness[0]=2; hardness[1]=1.45; hardness[2]=1.25; hardness[3]=1.95;
            //Or 3. Get Hardness from Input Argument
            hardness[0]=atof(argument[5]);/*Followings are "Don't Care":*/hardness[1]=100.0; hardness[2]=100.0; hardness[3]=100.0;

            // Use Hardness to define Cost
            // 1. Hardness with respect to input w
            // Tcosts[i] = hardness[i]*(me->GetInputWeight())-(hardness[i]-1);
            //Or 2. Hardness as the exact Cost
            Tcosts[i] = hardness[i];
            //Or 3. The exact Cost hardcoded
            // Tcosts[i] = 4.5;

            string type;
            if(i==0) type="Swamp";
            if(i==1) type="Water";
            if(i==2) type="Grass";
            if(i==3) type="Trees";
            // std::cout<<"Cost "<<type<<"="<<Tcosts[i]<<" ("<<hardness[i]<<"*"<<me->GetInputWeight()<<"-"<<(hardness[i]-1)<<")"<<std::endl;
        }
        me->SetTerrainCost(Tcosts);

        double proveBound = atof(argument[2]);
        tas.SetWeight(proveBound);
        me->SetInputWeight(atof(argument[2]));

        start = {1,99};
        goal = {198, 99};

        //Set the Baseline Policy.
        if(atoi(argument[1]) == 0){ //WA*
        tas.SetPhi([=](double h,double g){return g+proveBound*h;});
        }
        else if(atoi(argument[1]) == 1){ //PWXD
        tas.SetPhi([=](double h,double g){return (h>g)?(g+h):(g/proveBound+h*(2*proveBound-1)/proveBound);});
        }
        else if(atoi(argument[1]) == 2){ //PWXU
        tas.SetPhi([=](double h,double g){return (h*(2*proveBound-1)>g)?(g/(2*proveBound-1)+h):(1/proveBound*(g+h));});
        }
        else if(atoi(argument[1]) == 3){ //XDP
        tas.SetPhi([=](double h,double g){return (g+(2*proveBound-1)*h+sqrt((g-h)*(g-h)+4*proveBound*g*h))/(2*proveBound);});
        }
        else if(atoi(argument[1]) == 4){ //XUP
        tas.SetPhi([=](double h,double g){return (g+h+sqrt((g+h)*(g+h)+4*proveBound*(proveBound-1)*h*h))/(2*proveBound);});
        }
        else if(atoi(argument[1]) == 5){ //A*
        tas.SetPhi([=](double h,double g){return g+h;});
        }
        
        tas.InitializeSearch(me, start, goal, solution);

        tas.GetPath(me, start, goal, solution);
        printf("MAP %s #%d %1.2f ALG %d weight %1.2f Nodes %llu path %f\n", "Exp0", 0, 197, atoi(argument[1]), atof(argument[2]), tas.GetNodesExpanded(), me->GetPathLength(solution));
        exit(0);
    }
    else if (strcmp(argument[0], "-exp0DSD") == 0){
        Map *m = new Map(200, 200);
        MakeDesignedMap(m, atoi(argument[3]), atoi(argument[4]));
        // MakeDesignedMap(m, 30, 6);

        me = new MapEnvironment(m);
        // me->SetDiagonalCost(1.5);
        me->SetDiagonalCost(1.41);

        assert(maxNumArgs >= 5);
        //Set the cost of each terrain type randomly.
        for(int i=0; i<4; i++){
            //[0]=kSwamp, [1]=kWater, [2]=kGrass, [3]=kTrees

            // Define the Hardness
            // 1. Random Hardness
            // rdm = random()%101;
            // hardness[i] = rdm/101+1;
            //Or 2. Hardcoded Hardness
            // hardness[0]=2; hardness[1]=1.45; hardness[2]=1.25; hardness[3]=1.95;
            //Or 3. Get Hardness from Input Argument
            hardness[0]=atof(argument[5]);/*Followings are "Don't Care":*/hardness[1]=100.0; hardness[2]=100.0; hardness[3]=100.0;

            // Use Hardness to define Cost
            // 1. Hardness with respect to input w
            // Tcosts[i] = hardness[i]*(me->GetInputWeight())-(hardness[i]-1);
            //Or 2. Hardness as the exact Cost
            Tcosts[i] = hardness[i];
            //Or 3. The exact Cost hardcoded
            // Tcosts[i] = 4.5;

            string type;
            if(i==0) type="Swamp";
            if(i==1) type="Water";
            if(i==2) type="Grass";
            if(i==3) type="Trees";
            // std::cout<<"Cost "<<type<<"="<<Tcosts[i]<<" ("<<hardness[i]<<"*"<<me->GetInputWeight()<<"-"<<(hardness[i]-1)<<")"<<std::endl;
        }
        me->SetTerrainCost(Tcosts);

        dsd.policy = (tExpansionPriority)atoi(argument[1]);
        dsd.SetWeight(atof(argument[2]));
        me->SetInputWeight(atof(argument[2]));

        start = {1,99};
        goal = {198, 99};

        dsd.InitializeSearch(me, start, goal, solution);
        dsd.GetPath(me, start, goal, solution);
        printf("MAP %s #%d %1.2f ALG %d weight %1.2f Nodes %llu path %f\n", "Exp0", 0, 197, atoi(argument[1]), atof(argument[2]), dsd.GetNodesExpanded(), me->GetPathLength(solution));
        
        exit(0);
    }
    else if (strcmp(argument[0], "-exp0DPS") == 0){
        Map *m = new Map(200, 200);
        MakeDesignedMap(m, atoi(argument[2]), atoi(argument[3]));
        // MakeDesignedMap(m, 30, 6);

        me = new MapEnvironment(m);
        // me->SetDiagonalCost(1.5);
        me->SetDiagonalCost(1.41);

        assert(maxNumArgs >= 5);
        //Set the cost of each terrain type randomly.
        for(int i=0; i<4; i++){
            //[0]=kSwamp, [1]=kWater, [2]=kGrass, [3]=kTrees

            // Define the Hardness
            // 1. Random Hardness
            // rdm = random()%101;
            // hardness[i] = rdm/101+1;
            //Or 2. Hardcoded Hardness
            // hardness[0]=2; hardness[1]=1.45; hardness[2]=1.25; hardness[3]=1.95;
            //Or 3. Get Hardness from Input Argument
            hardness[0]=atof(argument[4]);/*Followings are "Don't Care":*/hardness[1]=100.0; hardness[2]=100.0; hardness[3]=100.0;

            // Use Hardness to define Cost
            // 1. Hardness with respect to input w
            // Tcosts[i] = hardness[i]*(me->GetInputWeight())-(hardness[i]-1);
            //Or 2. Hardness as the exact Cost
            Tcosts[i] = hardness[i];
            //Or 3. The exact Cost hardcoded
            // Tcosts[i] = 4.5;

            string type;
            if(i==0) type="Swamp";
            if(i==1) type="Water";
            if(i==2) type="Grass";
            if(i==3) type="Trees";
            // std::cout<<"Cost "<<type<<"="<<Tcosts[i]<<" ("<<hardness[i]<<"*"<<me->GetInputWeight()<<"-"<<(hardness[i]-1)<<")"<<std::endl;
        }
        me->SetTerrainCost(Tcosts);

        start = {1,99};
        goal = {198, 99};

        dps.SetOptimalityBound(atof(argument[1]));
        me->SetInputWeight(atof(argument[1]));

        dps.InitializeSearch(me, start, goal, solution);
            
        dps.GetPath(me, start, goal, solution);
        printf("MAP %s #%d %1.2f ALG %d weight %1.2f Nodes %llu path %f\n", "Exp0", 0, 197, 6, atof(argument[1]), dps.GetNodesExpanded(), me->GetPathLength(solution));
        
        exit(0);
    }
    else if (strcmp(argument[0], "-timeDSWA") == 0)
    {
        assert(maxNumArgs >= 4);
        
        DSDWAStar<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>> dsd_mnp;
        std::vector<MNPuzzleState<4, 4>> path;
        MNPuzzle<4, 4> mnp;
        MNPuzzleState<4, 4> start = STP::GetKorfInstance(atoi(argument[1]));
        MNPuzzleState<4, 4> goal;
        dsd_mnp.policy = (tExpansionPriority)atoi(argument[2]);
        dsd_mnp.SetWeight(atof(argument[3]));
        printf("Solving STP Korf instance [%d of %d] using DSD weight %f\n", atoi(argument[1])+1, 100, atof(argument[3]));
        
        clock_t start_time, end_time;
        start_time = clock();

        // dsd_mnp.GetPath(&mnp, start, goal, path, true);
        dsd_mnp.GetPath(&mnp, start, goal, path);

        end_time = clock();
        float runningTime = (float) (end_time - start_time) / CLOCKS_PER_SEC;

        printf("STP %d ALG %d weight %1.2f nodePERsec %f path %lu\n", atoi(argument[1]), atoi(argument[2]), atof(argument[3]), dsd_mnp.GetNodesExpanded()/runningTime, path.size());
        exit(0);
    }
    else if (strcmp(argument[0], "-stpAstar") == 0)
    {
        assert(maxNumArgs >= 4);
        
        TemplateAStar<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>> dsd_mnp;
        std::vector<MNPuzzleState<4, 4>> path;
        MNPuzzle<4, 4> mnp;
        MNPuzzleState<4, 4> start = STP::GetKorfInstance(atoi(argument[1]));
        MNPuzzleState<4, 4> goal;
        // dsd_mnp.policy = (tExpansionPriority)atoi(argument[2]);
        dsd_mnp.SetWeight(atof(argument[3]));
        printf("Solving STP Korf instance [%d of %d] using DSD weight %f\n", atoi(argument[1])+1, 100, atof(argument[3]));

        clock_t start_time, end_time;
        start_time = clock();

        // dsd_mnp.GetPath(&mnp, start, goal, path, true);
        dsd_mnp.GetPath(&mnp, start, goal, path);

        end_time = clock();
        float runningTime = (float) (end_time - start_time) / CLOCKS_PER_SEC;

        printf("STP %d ALG %d weight %1.2f nodePERsec %f path %lu\n", atoi(argument[1]), atoi(argument[2]), atof(argument[3]), dsd_mnp.GetNodesExpanded()/runningTime, path.size());
        exit(0);
    }
    else if (strcmp(argument[0], "-map") == 0)
    {
        assert(maxNumArgs >= 5);
        me = new MapEnvironment(new Map(argument[1]));
        // ScenarioLoader sl(argument[2]);
        sl = new ScenarioLoader(argument[2]);
        // std::cout<<"number of experiments is "<<sl.GetNumExperiments()<<std::endl;
        for (int x = 0; x < sl->GetNumExperiments(); x++)
        {
            Experiment exp = sl->GetNthExperiment(x);
            start.x = exp.GetStartX();
            start.y = exp.GetStartY();
            goal.x = exp.GetGoalX();
            goal.y = exp.GetGoalY();
            dsd.policy = (tExpansionPriority)atoi(argument[3]);
            dsd.SetWeight(atof(argument[4]));
            // dsd.GetPath(me, start, goal, solution, true);
            dsd.GetPath(me, start, goal, solution);
            printf("MAP %s #%d %1.2f ALG %d weight %1.2f Nodes %llu path %f\n", argument[1], x, exp.GetDistance(), atoi(argument[3]), atof(argument[4]), dsd.GetNodesExpanded(), me->GetPathLength(solution));
        }
        exit(0);
    }
    
    return 0;
}

void MyDisplayHandler(unsigned long windowID, tKeyboardModifier mod, char key)
{
    switch (key)
    {
        case 's':
            {
                //Reset the last problem's swamped states.
                if(exper==8 && (mapcmd<=5 || mapcmd==7)){
                    if(swampedloc.x == 0){
                        for(int i=0; i < me->GetMap()->GetMapWidth(); i++)
                            for(int j=swampedloc.y-int(terrain_height/2); j<=swampedloc.y+int(terrain_height/2); j++)
                                if(me->GetMap()->GetTerrainType(i, j) == kSwamp)
                                    me->GetMap()->SetTerrainType(i, j, kGround);
                    }
                    else if(swampedloc.y == 0){
                        for(int j=0; j < me->GetMap()->GetMapHeight(); j++)
                            for(int i=swampedloc.x-int(terrain_width/2); i<=swampedloc.x+int(terrain_width/2); i++)
                                if(me->GetMap()->GetTerrainType(i, j) == kSwamp)
                                    me->GetMap()->SetTerrainType(i, j, kGround);
                    }
                }
                
                //Set the start and goal states of the search.
                //load the scenario, or initiate it randomly
                if(mapcmd <= 4){
                    do {
                        start.x = random()%me->GetMap()->GetMapWidth();
                        start.y = random()%me->GetMap()->GetMapHeight();
                    } while (me->GetMap()->GetTerrainType(start.x, start.y) != kGround);
                    do {
                        goal.x = random()%me->GetMap()->GetMapWidth();
                        goal.y = random()%me->GetMap()->GetMapHeight();
                    } while (me->GetMap()->GetTerrainType(goal.x, goal.y) != kGround);
                }
                else if(mapcmd==5 || mapcmd==7){
                    numScenario += 1;
                    Experiment exp = sl->GetNthExperiment(numScenario);
                    start.x = exp.GetStartX();
                    start.y = exp.GetStartY();
                    goal.x = exp.GetGoalX();
                    goal.y = exp.GetGoalY();
                }

                //Change one or more squares of ground states to swamp type.
                if(exper==8 && (mapcmd<=5 || mapcmd==7)){
                    //Set the size of the swamp area using the terrain_size argument.
                    terrain_width = max(terrain_size/100*abs(goal.x-start.x), 10);
                    terrain_height = max(terrain_size/100*abs(goal.y-start.y), 10);

                    if(abs(goal.x-start.x) > abs(goal.y-start.y) && abs(goal.x-start.x)!=0){
                        // swampedloc.x = (goal.x+start.x)/2;
                        swampedloc.x = random()%(abs(goal.x-start.x))+min(goal.x,start.x);
                        for(int j=0; j < me->GetMap()->GetMapHeight(); j++)
                            for(int i=swampedloc.x-int(terrain_width/2); i<=swampedloc.x+int(terrain_width/2); i++)
                                if(me->GetMap()->GetTerrainType(i, j) == kGround)
                                    me->GetMap()->SetTerrainType(i, j, kSwamp);
                        swampedloc.y = 0;
                        // std::cout<<"columns\n";
                    
                    }
                    else if(abs(goal.y-start.y) > abs(goal.x-start.x) && abs(goal.y-start.y)!=0){
                        // swampedloc.y = (goal.y+start.y)/2;
                        swampedloc.y = random()%(abs(goal.y-start.y))+min(goal.y,start.y);
                        for(int i=0; i < me->GetMap()->GetMapWidth(); i++)
                            for(int j=swampedloc.y-int(terrain_height/2); j<=swampedloc.y+int(terrain_height/2); j++)
                                if(me->GetMap()->GetTerrainType(i, j) == kGround)
                                    me->GetMap()->SetTerrainType(i, j, kSwamp);
                        swampedloc.x = 0;
                        // std::cout<<"rows\n";
                    }
                    else if(abs(goal.x-start.x) > abs(goal.y-start.y)){
                        swampedloc.x = (goal.x+start.x)/2;
                        // swampedloc.x = random()%(abs(goal.x-start.x))+min(goal.x,start.x);
                        for(int j=0; j < me->GetMap()->GetMapHeight(); j++)
                            for(int i=swampedloc.x-int(terrain_width/2); i<=swampedloc.x+int(terrain_width/2); i++)
                                if(me->GetMap()->GetTerrainType(i, j) == kGround)
                                    me->GetMap()->SetTerrainType(i, j, kSwamp);
                        swampedloc.y = 0;
                    }
                    else{
                        swampedloc.y = (goal.y+start.y)/2;
                        // swampedloc.y = random()%(abs(goal.y-start.y))+min(goal.y,start.y);
                        for(int i=0; i < me->GetMap()->GetMapWidth(); i++)
                            for(int j=swampedloc.y-int(terrain_height/2); j<=swampedloc.y+int(terrain_height/2); j++)
                                if(me->GetMap()->GetTerrainType(i, j) == kGround)
                                    me->GetMap()->SetTerrainType(i, j, kSwamp);
                        swampedloc.x = 0;
                    }
                }
                
                if(mapcmd<=5){
                    problemNumber +=1;
                    printf("==============\n");
                    printf("Problem: %d\n", problemNumber);
                    printf("Policy: %d\n", dsd.policy);
                    printf("Bound: %f\n", bound);
                    for(int i=0; i<4; i++){
                        //[0]=kSwamp, [1]=kWater ,[2]=kGrass, [3]=kTrees
                        string type;
                        if(i==0) type="Swamp";
                        if(i==1) type="Water";
                        if(i==2) type="Grass";
                        if(i==3) type="Trees";
                        std::cout<<"Cost "<<type<<"="<<Tcosts[i]<<" ("<<hardness[i]<<"*"<<me->GetInputWeight()<<"-"<<(hardness[i]-1)<<")"<<std::endl;
                    }
                    if(useLookUpTable){
                        // look_up_table.clear();
                        // dsd.InitializeSearch_v2(me, start, goal, solution);
                        LookUpVector.clear();
                        dsd.InitializeSearch_v3(me, start, goal, solution);
                        
                    }
                    else{
                        data.resize(0);
                        dsd.InitializeSearch(me, start, goal, solution);
                    }
                }
                else if(mapcmd == 6){
                    //Does not apply to Racetrack
                    if(useLookUpTable){
                        look_up_table.clear();
                        dsd.InitializeSearch_v2(me, start, goal, solution);
                    }
                    else{
                        data.resize(0);
                        dsd.InitializeSearch(me, start, goal, solution);
                    }
                }
                else if(mapcmd == 7){
                    problemNumber +=1;
                    printf("==============\n");
                    printf("Problem: %d\n", problemNumber);
                    printf("DPS Algorithm\n");
                    printf("Bound: %f\n", bound);
                    for(int i=0; i<4; i++){
                        //[0]=kSwamp, [1]=kWater ,[2]=kGrass, [3]=kTrees
                        string type;
                        if(i==0) type="Swamp";
                        if(i==1) type="Water";
                        if(i==2) type="Grass";
                        if(i==3) type="Trees";
                        std::cout<<"Cost "<<type<<"="<<Tcosts[i]<<" ("<<hardness[i]<<"*"<<me->GetInputWeight()<<"-"<<(hardness[i]-1)<<")"<<std::endl;
                    }
                    data.resize(0);
                    dps.InitializeSearch(me, start, goal, solution);
                }
                else if(mapcmd == 8){
                    //Does not apply to Racetrack
                    data.resize(0);
                    dps_track.InitializeSearch(r, from, end, path);
                }
                
                searchRunning = true;
                break;
            }
        case 'r':
            {
                if(mapcmd <= 5){
                    printf("==============\n");
                    printf("Problem: %d\n", problemNumber);
                    printf("Policy: %d\n", dsd.policy);
                    printf("Bound: %f\n", bound);
                    for(int i=0; i<4; i++){
                        //[0]=kSwamp, [1]=kWater ,[2]=kGrass, [3]=kTrees
                        string type;
                        if(i==0) type="Swamp";
                        if(i==1) type="Water";
                        if(i==2) type="Grass";
                        if(i==3) type="Trees";
                        std::cout<<"Cost "<<type<<"="<<Tcosts[i]<<" ("<<hardness[i]<<"*"<<me->GetInputWeight()<<"-"<<(hardness[i]-1)<<")"<<std::endl;
                    }
                    if(useLookUpTable){
                        // look_up_table.clear();
                        // dsd.InitializeSearch_v2(me, start, goal, solution);
                        LookUpVector.clear();
                        dsd.InitializeSearch_v3(me, start, goal, solution);
                    }
                    else{
                        data.resize(0);
                        dsd.InitializeSearch(me, start, goal, solution);
                    }
                }
                else if(mapcmd == 6){
                    //RaceTrack
                    printf("==============\n");
                    printf("Problem: %d\n", problemNumber);
                    printf("Policy: %d\n", dsd_track.policy);
                    printf("Bound: %f\n", bound);
                    data.resize(0);
                    dsd_track.InitializeSearch(r, from, end, path);
                }
                else if(mapcmd == 7){
                    printf("==============\n");
                    printf("Problem: %d\n", problemNumber);
                    printf("DPS Algorithm\n");
                    printf("Bound: %f\n", bound);
                    data.resize(0);
                    dps.InitializeSearch(me, start, goal, solution);
                }
                else if(mapcmd == 8){
                    //RaceTrack
                    printf("==============\n");
                    printf("Problem: %d\n", problemNumber);
                    printf("DPS Algorithm\n");
                    printf("Bound: %f\n", bound);
                    data.resize(0);
                    dps_track.InitializeSearch(r, from, end, path);
                }
                
                searchRunning = true;
                break;
            }
        case 'w':
            {
                if(bound>=9) bound=1.25;
                else bound = (bound-1)*2+1;

                MyWindowHandler(windowID, kWindowDestroyed);
                MyWindowHandler(windowID, kWindowCreated);

                if(mapcmd <= 5){
                    //Set the input weight to the new bound
                    me->SetInputWeight(bound);
                    dsd.SetWeight(bound);

                    //Set the cost of each terrain type randomly.
                    for(int i=0; i<4; i++){
                        //[0]=kSwamp, [1]=kWater ,[2]=kGrass, [3]=kTrees

                        // 1. Random Hardness
                        // rdm = random()%101;
                        // hardness[i] = rdm/101+1;
                        //Or 2. Hardcoded Hardness
                        hardness[0]=1.5; hardness[1]=1.45; hardness[2]=1.25; hardness[3]=1.95;

                        Tcosts[i] = hardness[i]*(me->GetInputWeight())-(hardness[i]-1);
                        string type;
                        if(i==0) type="Swamp";
                        if(i==1) type="Water";
                        if(i==2) type="Grass";
                        if(i==3) type="Trees";
                        std::cout<<"Cost "<<type<<"="<<Tcosts[i]<<" ("<<hardness[i]<<"*"<<me->GetInputWeight()<<"-"<<(hardness[i]-1)<<")"<<std::endl;
                    }
                    me->SetTerrainCost(Tcosts);

                    if(mapcmd<=4) problemNumber = 0;

                    printf("Problem: %d\n", problemNumber);
                    printf("Policy: %d\n", dsd.policy);
                    printf("Bound: %f\n", bound);
                    // data.resize(0);
                    // dsd.InitializeSearch(me, start, goal, solution);
                    look_up_table.clear();
                    dsd.InitializeSearch_v2(me, start, goal, solution);
                }
                else if(mapcmd == 6){
                    //RaceTrack

                    //Set the input weight to the new bound
                    dsd_track.SetWeight(bound);

                    printf("Problem: %d\n", problemNumber);
                    printf("Policy: %d\n", dsd_track.policy);
                    printf("Bound: %f\n", bound);
                    data.resize(0);
                    dsd_track.InitializeSearch(r, from, end, path);
                }
                else if(mapcmd == 7){
                    //Set the input weight to the new bound
                    dps.SetOptimalityBound(bound);

                    printf("==============\n");
                    printf("Problem: %d\n", problemNumber);
                    printf("DPS Algorithm\n");
                    printf("Bound: %f\n", bound);
                    data.resize(0);
                    dps.InitializeSearch(me, start, goal, solution);
                }
                else if(mapcmd == 8){
                    //RaceTrack

                    //Set the input weight to the new bound
                    dps_track.SetOptimalityBound(bound);

                    printf("Problem: %d\n", problemNumber);
                    printf("DPS Algorithm\n");
                    printf("Bound: %f\n", bound);
                    data.resize(0);
                    dps_track.InitializeSearch(r, from, end, path);
                }
                
                searchRunning = true;
                break;
            }
        case 'p': showPlane = !showPlane; break;
        case '[': stepsPerFrame = std::max(stepsPerFrame/2, 1); break;
        case ']': stepsPerFrame = stepsPerFrame*2; break;
        case '}':
            {
                if(mapcmd <= 5){
                    dsd.policy = (tExpansionPriority)((dsd.policy+1)%kDSDPolicyCount);
                    printf("Problem: %d\n", problemNumber);
                    printf("Policy: %d\n", dsd.policy);
                    printf("Bound: %f\n", bound);
                    if(useLookUpTable){
                        // look_up_table.clear();
                        // dsd.InitializeSearch_v2(me, start, goal, solution);
                        LookUpVector.clear();
                        dsd.InitializeSearch_v3(me, start, goal, solution);
                    }
                    else{
                        data.resize(0);
                        dsd.InitializeSearch(me, start, goal, solution);
                    }
                }
                else if(mapcmd == 6){
                    //RaceTrack
                    dsd_track.policy = (tExpansionPriority)((dsd_track.policy+1)%kDSDPolicyCount);
                    printf("Problem: %d\n", problemNumber);
                    printf("Policy: %d\n", dsd_track.policy);
                    printf("Bound: %f\n", bound);
                    data.resize(0);
                    dsd_track.InitializeSearch(r, from, end, path);
                }
                else if(mapcmd == 7){
                    //Does not apply here
                    printf("Problem: %d\n", problemNumber);
                    printf("DPS Algorithm\n");
                    printf("Bound: %f\n", bound);
                    data.resize(0);
                    dps.InitializeSearch(me, start, goal, solution);
                }
                else if(mapcmd == 8){
                    //Does not apply here
                    printf("Problem: %d\n", problemNumber);
                    printf("DPS Algorithm\n");
                    printf("Bound: %f\n", bound);
                    data.resize(0);
                    dps_track.InitializeSearch(r, from, end, path);
                }
                
                searchRunning = true;
                break;
            }
        case '{':
            {
                if(mapcmd <= 5){
                    dsd.policy = (tExpansionPriority)((dsd.policy-1)%kDSDPolicyCount);
                    printf("Problem: %d\n", problemNumber);
                    printf("Policy: %d\n", dsd.policy);
                    printf("Bound: %f\n", bound);
                    if(useLookUpTable){
                        // look_up_table.clear();
                        // dsd.InitializeSearch_v2(me, start, goal, solution);
                        LookUpVector.clear();
                        dsd.InitializeSearch_v3(me, start, goal, solution);
                    }
                    else{
                        data.resize(0);
                        dsd.InitializeSearch(me, start, goal, solution);
                    }
                    
                }
                else if(mapcmd == 6){
                    //RaceTrack
                    dsd_track.policy = (tExpansionPriority)((dsd_track.policy-1)%kDSDPolicyCount);
                    printf("Problem: %d\n", problemNumber);
                    printf("Policy: %d\n", dsd_track.policy);
                    printf("Bound: %f\n", bound);
                    data.resize(0);
                    dsd_track.InitializeSearch(r, from, end, path);
                }
                else if(mapcmd == 7){
                    //Does not apply here
                    printf("Problem: %d\n", problemNumber);
                    printf("DPS Algorithm\n");
                    printf("Bound: %f\n", bound);
                    data.resize(0);
                    dps.InitializeSearch(me, start, goal, solution);
                }
                else if(mapcmd == 8){
                    //Does not apply here
                    printf("Problem: %d\n", problemNumber);
                    printf("DPS Algorithm\n");
                    printf("Bound: %f\n", bound);
                    data.resize(0);
                    dps_track.InitializeSearch(r, from, end, path);
                }
                
                searchRunning = true;
                break;
            }
        case '+':
        {
            //load the scenario, or initiate it randomly
            if(mapcmd <= 4){
                do {
                    start.x = random()%me->GetMap()->GetMapWidth();
                    start.y = random()%me->GetMap()->GetMapHeight();
                } while (me->GetMap()->GetTerrainType(start.x, start.y) != kGround);
                do {
                    goal.x = random()%me->GetMap()->GetMapWidth();
                    goal.y = random()%me->GetMap()->GetMapHeight();
                } while (me->GetMap()->GetTerrainType(goal.x, goal.y) != kGround);
            }
            else if(mapcmd == 5 || mapcmd == 7){
                numScenario = (numScenario + sl->GetNumExperiments()/20+1) % sl->GetNumExperiments();
                Experiment exp = sl->GetNthExperiment(numScenario);
                printf("==============\n");
                std::cout<<"Path Length is "<<exp.GetDistance()<<" (scen "<<numScenario<<"/"<<sl->GetNumExperiments()<<")\n";
                start.x = exp.GetStartX();
                start.y = exp.GetStartY();
                goal.x = exp.GetGoalX();
                goal.y = exp.GetGoalY();
            }
            else if(mapcmd == 6 || mapcmd == 8){
                //Racetrack
                numScenario = (numScenario + sl->GetNumExperiments()/20+1) % sl->GetNumExperiments();
                // numScenario = numScenario + +1 % sl->GetNumExperiments();
                Experiment exp = sl->GetNthExperiment(numScenario);
                std::cout<<"Path Length is "<<exp.GetDistance()<<" (scen "<<numScenario<<"/"<<sl->GetNumExperiments()<<")\n";
                if(numExtendedGoals) numExtendedGoals = exp.GetDistance()/3;

                //Reset the previous start and goal
                m->SetTerrainType(from.xLoc, from.yLoc, kGround);
                if(numExtendedGoals == 0){//Not extended the goal.
                    m->SetTerrainType(end.xLoc, end.yLoc, kGround);
                }
                else{//Extended the goal
                    for(int tNode=0; tNode<theList.size(); tNode++){
                        m->SetTerrainType(theList[tNode].x, theList[tNode].y, kGround);
                    }
                }

                //Set the start and goal states, weight and policy of the search.
                from.xLoc = exp.GetStartX();
                from.yLoc = exp.GetStartY();
                from.xVelocity = 0;
                from.yVelocity = 0;
                end.xLoc = exp.GetGoalX();
                end.yLoc = exp.GetGoalY();
                end.xVelocity = 0;
                end.yVelocity = 0;
                
                //Set the start and goal TerrainType.
                m->SetTerrainType(from.xLoc, from.yLoc, kStartTerrain);
                if(numExtendedGoals == 0){//Not extend the goal.
                    m->SetTerrainType(end.xLoc, end.yLoc, kEndTerrain);
                }
                else{//Extend the goal
                    tas.SetPhi([=](double h,double g){return g;});
                    goal.x = end.xLoc;
                    goal.y = end.yLoc;
                    tas.ExtendGoal(me, goal, theList, numExtendedGoals);
                    for(int tNode=0; tNode<theList.size(); tNode++){
                        m->SetTerrainType(theList[tNode].x, theList[tNode].y, kEndTerrain);
                    }
                }
                
                r->UpdateMap(m);
            }

            if(mapcmd <= 5){
                problemNumber +=1;
                printf("Problem: %d\n", problemNumber);
                printf("Policy: %d\n", dsd.policy);
                printf("Bound: %f\n", bound);
                if(useLookUpTable){
                    // look_up_table.clear();
                    // dsd.InitializeSearch_v2(me, start, goal, solution);
                    LookUpVector.clear();
                    dsd.InitializeSearch_v3(me, start, goal, solution);
                }
                else{
                    data.resize(0);
                    dsd.InitializeSearch(me, start, goal, solution);
                }
            }
            else if(mapcmd == 6){
                //RaceTrack
                problemNumber += 1;
                printf("Problem: %d\n", problemNumber);
                printf("Policy: %d\n", dsd_track.policy);
                printf("Bound: %f\n", bound);
                data.resize(0);
                dsd_track.InitializeSearch(r, from, end, path);
            }
            else if(mapcmd == 7){
                problemNumber += 1;
                printf("Problem: %d\n", problemNumber);
                printf("DPS Algorithm\n");
                printf("Bound: %f\n", bound);
                data.resize(0);
                dps.InitializeSearch(me, start, goal, solution);
            }
            else if(mapcmd == 8){
                problemNumber += 1;
                printf("Problem: %d\n", problemNumber);
                printf("DPS Algorithm\n");
                printf("Bound: %f\n", bound);
                data.resize(0);
                dps_track.InitializeSearch(r, from, end, path);
            }

            searchRunning = true;
            break;
        }
        case '-':
        {
            //load the scenario, or initiate it randomly
            if(mapcmd <= 4){
                do {
                    start.x = random()%me->GetMap()->GetMapWidth();
                    start.y = random()%me->GetMap()->GetMapHeight();
                } while (me->GetMap()->GetTerrainType(start.x, start.y) != kGround);
                do {
                    goal.x = random()%me->GetMap()->GetMapWidth();
                    goal.y = random()%me->GetMap()->GetMapHeight();
                } while (me->GetMap()->GetTerrainType(goal.x, goal.y) != kGround);
            }
            else if(mapcmd == 5 || mapcmd == 7){
                numScenario = (numScenario - sl->GetNumExperiments()/20-1+sl->GetNumExperiments())  % sl->GetNumExperiments();
                Experiment exp = sl->GetNthExperiment(numScenario);
                printf("==============\n");
                std::cout<<"Path Length is "<<exp.GetDistance()<<" (scen "<<numScenario<<"/"<<sl->GetNumExperiments()<<")\n";
                start.x = exp.GetStartX();
                start.y = exp.GetStartY();
                goal.x = exp.GetGoalX();
                goal.y = exp.GetGoalY();
            }
            else if(mapcmd == 6 || mapcmd == 8){
                //Racetrack
                numScenario = (numScenario - sl->GetNumExperiments()/20-1+sl->GetNumExperiments()) % sl->GetNumExperiments();
                // numScenario -= 1;
                Experiment exp = sl->GetNthExperiment(numScenario);
                printf("==============\n");
                std::cout<<"Path Length is "<<exp.GetDistance()<<" (scen "<<numScenario<<"/"<<sl->GetNumExperiments()<<")\n";
                if(numExtendedGoals) numExtendedGoals = exp.GetDistance()/3;

                //Reset the previous start and goal
                m->SetTerrainType(from.xLoc, from.yLoc, kGround);
                if(numExtendedGoals == 0){//Not extended the goal.
                    m->SetTerrainType(end.xLoc, end.yLoc, kGround);
                }
                else{//Extended the goal
                    for(int tNode=0; tNode<theList.size(); tNode++){
                        m->SetTerrainType(theList[tNode].x, theList[tNode].y, kGround);
                    }
                }

                //Set the start and goal states, weight and policy of the search.
                from.xLoc = exp.GetStartX();
                from.yLoc = exp.GetStartY();
                from.xVelocity = 0;
                from.yVelocity = 0;
                end.xLoc = exp.GetGoalX();
                end.yLoc = exp.GetGoalY();
                end.xVelocity = 0;
                end.yVelocity = 0;
                
                //Set the start and goal TerrainType.
                m->SetTerrainType(from.xLoc, from.yLoc, kStartTerrain);
                if(numExtendedGoals == 0){//Not extend the goal.
                    m->SetTerrainType(end.xLoc, end.yLoc, kEndTerrain);
                }
                else{//Extend the goal
                    tas.SetPhi([=](double h,double g){return g;});
                    goal.x = end.xLoc;
                    goal.y = end.yLoc;
                    tas.ExtendGoal(me, goal, theList, numExtendedGoals);
                    for(int tNode=0; tNode<theList.size(); tNode++){
                        m->SetTerrainType(theList[tNode].x, theList[tNode].y, kEndTerrain);
                    }
                }
                
                r->UpdateMap(m);
            }

            if(mapcmd <= 5){
                problemNumber -=1;
                printf("Problem: %d\n", problemNumber);
                printf("Policy: %d\n", dsd.policy);
                printf("Bound: %f\n", bound);
                if(useLookUpTable){
                    // look_up_table.clear();
                    // dsd.InitializeSearch_v2(me, start, goal, solution);
                    LookUpVector.clear();
                    dsd.InitializeSearch_v3(me, start, goal, solution);
                }
                else{
                    data.resize(0);
                    dsd.InitializeSearch(me, start, goal, solution);
                }
            }
            else if(mapcmd == 6){
                //RaceTrack
                problemNumber -= 1;
                printf("Problem: %d\n", problemNumber);
                printf("Policy: %d\n", dsd_track.policy);
                printf("Bound: %f\n", bound);
                data.resize(0);
                dsd_track.InitializeSearch(r, from, end, path);
            }
            else if(mapcmd == 7){
                problemNumber -= 1;
                printf("Problem: %d\n", problemNumber);
                printf("DPS Algorithm\n");
                printf("Bound: %f\n", bound);
                data.resize(0);
                dps.InitializeSearch(me, start, goal, solution);
            }
            else if(mapcmd == 8){
                problemNumber -= 1;
                printf("Problem: %d\n", problemNumber);
                printf("DPS Algorithm\n");
                printf("Bound: %f\n", bound);
                data.resize(0);
                dps_track.InitializeSearch(r, from, end, path);
            }
            
            searchRunning = true;
            break;
        }
        case 't':
        {
            useLookUpTable = !useLookUpTable;

            if(mapcmd <= 5){
                printf("==============\n");
                printf("Problem: %d\n", problemNumber);
                printf("Policy: %d\n", dsd.policy);
                printf("Bound: %f\n", bound);
                for(int i=0; i<4; i++){
                    //[0]=kSwamp, [1]=kWater ,[2]=kGrass, [3]=kTrees
                    string type;
                    if(i==0) type="Swamp";
                    if(i==1) type="Water";
                    if(i==2) type="Grass";
                    if(i==3) type="Trees";
                    std::cout<<"Cost "<<type<<"="<<Tcosts[i]<<" ("<<hardness[i]<<"*"<<me->GetInputWeight()<<"-"<<(hardness[i]-1)<<")"<<std::endl;
                }
                if(useLookUpTable){
                    // look_up_table.clear();
                    // dsd.InitializeSearch_v2(me, start, goal, solution);
                    LookUpVector.clear();
                    dsd.InitializeSearch_v3(me, start, goal, solution);
                }
                else{
                    data.resize(0);
                    dsd.InitializeSearch(me, start, goal, solution);
                }
            }
            else if(mapcmd == 6){
                //RaceTrack
                printf("==============\n");
                printf("Problem: %d\n", problemNumber);
                printf("Policy: %d\n", dsd_track.policy);
                printf("Bound: %f\n", bound);
                data.resize(0);
                dsd_track.InitializeSearch(r, from, end, path);
            }
            else if(mapcmd == 7){
                printf("==============\n");
                printf("Problem: %d\n", problemNumber);
                printf("DPS Algorithm\n");
                printf("Bound: %f\n", bound);
                data.resize(0);
                dps.InitializeSearch(me, start, goal, solution);
            }
            else if(mapcmd == 8){
                //RaceTrack
                printf("==============\n");
                printf("Problem: %d\n", problemNumber);
                printf("DPS Algorithm\n");
                printf("Bound: %f\n", bound);
                data.resize(0);
                dps_track.InitializeSearch(r, from, end, path);
            }
            
            searchRunning = true;
            break;
        }
        default:
            break;
    
    }
}

void SetNextPriority(float h, float g, float target) // const Graphics::point &loc
{
    float slope = g/h;
    if (h > 0 && g > 0)
    {
        if (data.size() == 0 || data.back().slope < slope)
        {
            //std::cout << "Virtual hit of " << loc << " slope " << slope << "\n";
            float minWeight, maxWeight;
            if (data.size() > 0)
                GetNextWeightRange(minWeight, maxWeight, data.back().crossPoint, slope);
            else
                GetNextWeightRange(minWeight, maxWeight, {1, 0}, slope);
                        
            float K;
            // K (g + [w] * h) = 1 at previous point
            point3d last;
            if (data.size() == 0)
            {
                last = point3d(1, 0);
            }
            else {
                last = data.back().crossPoint;
            }
            // returns nextWeight and K
            float nextWeight = ChooseWeightForTargetPriority({h, g}, target, minWeight, maxWeight, last, K);
            
            // our cross point of next slope
            point3d crossPoint1;
            crossPoint1.x = 1.0f/(K*(slope+nextWeight));
            crossPoint1.y = crossPoint1.x*slope;
            
            data.push_back({slope, nextWeight, K, crossPoint1});
            
//            std::cout << "Cross Priorities: ";
//            for (const auto &i : data)
//            {
//                std::cout << i.crossPoint << ": ";
//                std::cout << GetPriority(i.crossPoint.x, i.crossPoint.y) << " ";
//            }
//            std::cout << "\n";
        }
    }
}

bool MyClickHandler(unsigned long windowID, int viewport, int, int, point3d loc, tButtonType button, tMouseEventType mType)
{
//    return false;
//    static point3d startLoc;
//    loc.x = (loc.x+1)/2;
//    loc.y = 1-(loc.y+1)/2;
//    loc *= 2;
    loc = HOGToLocal(loc);
    if (mType == kMouseDown)
    {
        switch (button)
        {
            case kRightButton: //printf("Right button\n"); break;
            {
                std::cout << "Priority of " << loc << " is " << GetPriority(loc.x, loc.y) << "\n";
            }
                break;
            case kLeftButton: //printf("Left button\n");
                
                if (viewport == 3) // TODO: add interactive mode back in later
                    SetNextPriority(loc.x, loc.y, testScale); // h and g
                break;
            case kMiddleButton: printf("Middle button\n"); break;
            case kNoButton: break;
        }
    }
    if ((button == kMiddleButton) && (mType == kMouseDown))
    {}
    if (button == kRightButton)
    {}
    return false;
}

float GetPriority(float h, float g)
{
    if (data.size() == 0)
        return INFINITY;
    float slope = g/h;
    if (fgreater(slope, data.back().slope))
        return INFINITY;
    // range includes low but not high
//    int low = 0, high = data.size();
    // dumb/slow but correct
    for (int x = 0; x < data.size(); x++)
        if (flesseq(slope, data[x].slope))
            return data[x].K*(g + data[x].weight * h);
//    while (true)
//    {
//        // f = K * ( g + w_i * h )
//        if (low >= high-1)
//        {
//            return data[low].K*(g + data[low].weight * h);
//        }
//        int mid = (low+high)/2;
//        if (data[mid].slope > slope)
//        {
//            high = mid+1;
//        }
//        else {
//            low = mid+1;
//        }
//    }
    return INFINITY;
//    return data[low].K*(g + data[low].weight * h);
}

/**
 * Given the slope of the next bounding line, give the possbile range of weights that can be used in the priority function
 *
 * \param minWeight (returned) The minimum weight that can be used without going under the lower limit
 * \param maxWeight (returned) The maximum weight that can be used without going over the upper limit
 * \param currPoint The point on the previous bounding line with priorirty 1.0
 * \param nextSlope The slope of the next bounding line
 **/
void GetNextWeightRange(float &minWeight, float &maxWeight, point3d currPoint, float nextSlope)
{
    // defaults
    minWeight = 1;
    maxWeight = 2*bound-1;

    // 0. next slope line is y = slope*x
    // 1. cannot go over the [y = -x + w] line
    // slope*x = -x + w; slope*x + x = w; x = w/(slope+1)
    // y/slope = w-y; y = slope * w - slope *y; y*(1+slope) = slope*w; y = slope*w/(1+slope)
    point3d upperPoint(bound/(nextSlope+1),
                       nextSlope*bound/(1+nextSlope));
//    std::cout << "Upper point: " << upperPoint << "\n";
    // 2. cannot go under the [y = -(2w-1)x + w] line
    // slope*x = -(2w-1)x + w; x*(slope + (2w-1)) = w
    // y/slope = (w-y)/(2w-1)
    // y (2w-1) = slope * w - slope*y
    // y (slope + 2w-1) = slope*w
    point3d lowerPoint(bound/(nextSlope+2*bound-1),
                       nextSlope*bound / (nextSlope + 2*bound-1));
    // get (negative) slopes to upper and lower points
    minWeight = std::max(minWeight, (lowerPoint.y-currPoint.y)/(currPoint.x-lowerPoint.x));
    if (upperPoint.x < currPoint.x)
        maxWeight = std::min(maxWeight, (upperPoint.y-currPoint.y)/(currPoint.x-upperPoint.x));
//    printf("Weight needs to be [%f, %f]\n", minWeight, maxWeight);
}

/*
 * Given location, and a range of priorities, try to give the point the exact desired priority
 * Returned priority will always be in the minWeight/maxWeight range
 */
float ChooseWeightForTargetPriority(point3d loc, float priority, float minWeight, float maxWeight, point3d last, float &K)
{
    float weight;
    // New: Last point is the priority 1 crossing point.
    // 1. Find the point on the previous line with priority {priority}
    point3d projectedPoint = last*priority;
    // 2. FYI: loc is the point on the new line that we want to have the desired priority
    // 3. Find the slope between the last point and our new point.
    //    Which is delta y / delta x
    if (flesseq(loc.y, projectedPoint.y))
    {
        printf("Ill defined case (new y < old y); defaulting to min\n");
        K = 1/(last.y+minWeight*last.x);
        return minWeight;
    }
    if (fgreatereq(loc.x, projectedPoint.x))
    {
        printf("Ill defined case (new x > old x); defaulting to max\n");
        K = 1/(last.y+maxWeight*last.x);
        return maxWeight;
    }
    
    // Then try to extend that point to the new point giving it the desired priority
    weight = (loc.y-projectedPoint.y)/(projectedPoint.x-loc.x);
    // and bound by min/max weight
    weight = std::max(std::min(weight, maxWeight), minWeight);
    K = 1/(last.y+weight*last.x);
    return weight;

/*
    // Old approach. (Keeping to verify should be same for priority=1)
    // Only difference is in checking for ill defined cases
    // Requires that loc be above and to the left of last
    if (loc.x > last.x)
    {
        K = 1/(last.y+maxWeight*last.x);
        return maxWeight;
    }
    if (loc.y < last.y)
    {
        K = 1/(last.y+minWeight*last.x);
        return minWeight;
    }
    // Old: Find the crossing point of two equations of the form:
    // f = K * (y +  w*x);
    // K = 1/(last.y+nextWeight*last.x);
    // priority = K*(loc.y + weight * loc.x);
    // priority = (loc.y + weight * loc.x)/(last.y+weight*last.x);
    // priority*last.y+priority*weight*last.x = loc.y + weight * loc.x;
    // priority*weight*last.x-weight * loc.x = loc.y - priority*last.y
    // weight*(priority*last.x-loc.x) = loc.y - priority*last.y
    weight = (loc.y - priority*last.y)/(priority*last.x-loc.x);
    weight = std::max(std::min(weight, maxWeight), minWeight);
    K = 1/(last.y+weight*last.x);
    return weight;
*/
 }
