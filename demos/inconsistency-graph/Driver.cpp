/*
 *  $Id: sample.cpp
 *  hog2
 *
 *  Created by Nathan Sturtevant on 5/31/05.
 *  Modified by Nathan Sturtevant on 02/29/20.
 *
 * This file is part of HOG2. See https://github.com/nathansttt/hog2 for licensing information.
 *
 */

#include "Common.h"
#include "Driver.h"
#include "GraphEnvironment.h"
#include "TemplateAStar.h"
#include "TextOverlay.h"
#include <string>
#include "SVGUtil.h"
#include "GraphInconsistencyInstances.h"
#include "IncrementalBGS.h"

enum mode {
	kAddNodes,
	kAddEdges,
	kMoveNodes,
	kFindPath,
	kFindPathIbex
};

TemplateAStar<graphState, graphMove, GraphEnvironment> astar;
IncrementalBGS<graphState, graphMove> ibex;

std::vector<graphState> path;
const int kHeuristic = GraphSearchConstants::kTemporaryLabel;
std::string MyToString(double val);

mode m = kFindPath;

bool recording = false;
bool running = false;
bool paused = true;
bool BPMX = false;
double edgeCost = 1.0;
double weight = 1.0;
int ticksPerNode = 1;
int currentTick = 0;

int numExampleNodes = 5;
char graphType = 'b';
Graph *g = 0;
GraphEnvironment *ge;
graphState from=-1, to=-1;

TextOverlay te(35);

void SaveGraph(const char *file);
void LoadGraph(const char *file);
void ShowSearchInfo();
static char name[2] = "a";

double distance(unsigned long n1, unsigned long n2);

class GraphDistHeuristic : public Heuristic<graphState> {
public:
	double HCost(const graphState &a, const graphState &b) const
	{
		if (graphType == 'w')
			return g->GetNode(a)->GetLabelF(kHeuristic);
		return g->GetNode(a)->GetLabelL(kHeuristic);
	}
};

GraphDistHeuristic h;

int main(int argc, char* argv[])
{
	InstallHandlers();
	RunHOGGUI(argc, argv, 1600, 800);
	return 0;
}

/**
 * Allows you to install any keyboard handlers needed for program interaction.
 */
void InstallHandlers()
{
	InstallKeyboardHandler(MyDisplayHandler, "Record", "Record a movie", kAnyModifier, 'r');
	InstallKeyboardHandler(MyDisplayHandler, "Pause Simulation", "Pause simulation execution.", kNoModifier, 'p');
	InstallKeyboardHandler(MyDisplayHandler, "Step Simulation", "If the simulation is paused, step forward .1 sec.", kAnyModifier, 'o');
	InstallKeyboardHandler(MyDisplayHandler, "Scale+", "Larger scale", kAnyModifier, '}');
	InstallKeyboardHandler(MyDisplayHandler, "Scale-", "Smaller scale", kAnyModifier, '{');
	InstallKeyboardHandler(MyDisplayHandler, "Larger", "Larger examples", kAnyModifier, ']');
	InstallKeyboardHandler(MyDisplayHandler, "Smaller", "Smaller examples", kAnyModifier, '[');
	InstallKeyboardHandler(MyDisplayHandler, "Clear", "Clear graph", kAnyModifier, '|');
	InstallKeyboardHandler(MyDisplayHandler, "BPMX", "Toggle BPMX", kAnyModifier, 'x');
	InstallKeyboardHandler(MyDisplayHandler, "IBEX", "Run IBEX", kAnyModifier, 'i');
//	InstallKeyboardHandler(MyDisplayHandler, "Help", "Draw help", kAnyModifier, '?');
//	InstallKeyboardHandler(MyDisplayHandler, "Weight", "Toggle Dijkstra & A*", kAnyModifier, 'w');
//	InstallKeyboardHandler(MyDisplayHandler, "Save", "Save current graph", kAnyModifier, 's');
//	InstallKeyboardHandler(MyDisplayHandler, "Load", "Load last saved graph", kAnyModifier, 'l');
	InstallKeyboardHandler(BuildGraphFromPuzzle, "Incons.", "Build Inconsistency Graph", kAnyModifier, 'a', 'd');
	InstallKeyboardHandler(BuildGraphFromPuzzle, "WA*", "Build WA* graph", kAnyModifier, 'w');

	//InstallCommandLineHandler(MyCLHandler, "-map", "-map filename", "Selects the default map to be loaded.");
	
	InstallWindowHandler(MyWindowHandler);

	InstallMouseClickHandler(MyClickHandler);
}

void MyWindowHandler(unsigned long windowID, tWindowEventType eType)
{
	if (eType == kWindowDestroyed)
	{
		printf("Window %ld destroyed\n", windowID);
		RemoveFrameHandler(MyFrameHandler, windowID, 0);
	}
	else if (eType == kWindowCreated)
	{
		ReinitViewports(windowID, {-1, -1, 0.0, 1}, kScaleToSquare);
		AddViewport(windowID, {0, -1, 1, 1}, kScaleToSquare);

		printf("Window %ld created\n", windowID);
		glClearColor(0.99, 0.99, 0.99, 1.0);
		InstallFrameHandler(MyFrameHandler, windowID, 0);
		SetNumPorts(windowID, 2);
		g = new Graph();
		ge = new GraphEnvironment(g);
		ge->SetDirected(true);
		ge->SetDrawEdgeCosts(true);
		ge->SetDrawNodeLabels(true);
		ge->SetIntegerEdgeCosts(true);
		astar.SetReopenNodes(true);
		astar.SetDirected(true);
		astar.SetWeight(1.0);
		astar.SetHeuristic(&h);
		BuildGraphFromPuzzle(windowID, kNoModifier, graphType);
//		te.AddLine("A* with Inconsistent Heuristics");
//		te.AddLine("Current mode: add nodes (click to add node)");
//		te.AddLine("Press [ or ] to change modes. '?' for help.");
	}
}

int frameCnt = 0;

void MyFrameHandler(unsigned long windowID, unsigned int viewport, void *)
{
	Graphics::Display &display = getCurrentContext()->display;

	if (viewport == 0)
	{
		display.FillRect({-1, -1, 1, 1}, Colors::white);

		if (ge == 0 || g == 0)
			return;
		ge->SetColor(Colors::lightgray);
		ge->Draw(display);
				
		ge->SetColor(Colors::white);
		for (int x = 0; x < g->GetNumNodes(); x++)
		{
			ge->Draw(display, x);
		}
		if (running)
		{
			if (!paused)
			{
				currentTick++;
				if (0 == currentTick%ticksPerNode)
					MyDisplayHandler(windowID, kNoModifier, 'o');
			}
			if (m == kFindPath)
				astar.Draw(display);
			else if (m == kFindPathIbex)
			{
				ibex.Draw(display);
				std::string tmp = "f-limit: "+ibex.fEquation;
				//tmp += ibex.GetCurrentFLimit()!=DBL_MAX?to_string_with_precision(ibex.GetCurrentFLimit(), 2):"∞";
				display.DrawText(tmp.c_str(), {-1+0.05, -0.9}, Colors::black, 0.05, Graphics::textAlignLeft);

				uint64_t small, large;
				ibex.GetNodeInterval(small, large);
				tmp = "nodes: ["+std::to_string(small) + ","+ (large==ibex.infiniteWorkBound?"∞":std::to_string(large)) +"]";
				display.DrawText(tmp.c_str(), {-1+0.05, -0.84}, Colors::black, 0.05, Graphics::textAlignLeft);


				tmp = "expand: "+std::to_string(ibex.GetIterationNodesExpanded());
				display.DrawText(tmp.c_str(), {-1+0.05, -0.78}, Colors::black, 0.05, Graphics::textAlignLeft);
			}
		}
		
		if (path.size() > 0)
		{
			ge->SetColor(0, 1, 0);
			for (int x = 1; x < path.size(); x++)
			{
				ge->DrawLine(display, path[x-1], path[x], 10);
			}
		}

		for (int x = 0; x < g->GetNumNodes(); x++)
		{
			ge->SetColor(Colors::darkpurple);
			if (graphType == 'w')
				ge->DrawStateLabel(display, x, MyToString(g->GetNode(x)->GetLabelF(kHeuristic)).c_str());
			else
				ge->DrawStateLabel(display, x, MyToString(g->GetNode(x)->GetLabelL(kHeuristic)).c_str());
		}
	}
	if (viewport == 1)
	{
		te.Draw(display);
	}
	
	if (recording && viewport == GetNumPorts(windowID)-1)
	{
		char fname[255];
		sprintf(fname, "/Users/nathanst/Movies/tmp/inconsistency-%d%d%d%d.svg",
				(frameCnt/1000)%10, (frameCnt/100)%10, (frameCnt/10)%10, frameCnt%10);
		
		MakeSVG(GetContext(windowID)->display, fname, 1200, 1200, 0, "", false);

//		SaveScreenshot(windowID, fname);
		printf("Saved %s\n", fname);
		frameCnt++;
//		if (path.size() == 0)
//		{
//			MyDisplayHandler(windowID, kNoModifier, 'o');
//		}
//		else {
//			recording = false;
//		}
	}
	return;
	
}

int MyCLHandler(char *argument[], int maxNumArgs)
{
	if (maxNumArgs <= 1)
		return 0;
	strncpy(gDefaultMap, argument[1], 1024);
	return 2;
}

void MyDisplayHandler(unsigned long windowID, tKeyboardModifier mod, char key)
{
	switch (key)
	{
		case '{':
		{
			ge->SetNodeScale(ge->GetNodeScale()*0.9);
			break;
		}
		case '}':
		{
			ge->SetNodeScale(ge->GetNodeScale()*1.1);
			break;
		}
		case ']':
		{
			if (numExampleNodes < 15)
			{
				numExampleNodes++;
				BuildGraphFromPuzzle(windowID, kNoModifier, graphType);
				if (numExampleNodes > 8)
				{
					ge->SetNodeScale(2.5);
					ge->SetDrawEdgeCosts(false);
					ge->SetDrawNodeLabels(false);
				}
				else {
					ge->SetNodeScale(1.75);
					ge->SetDrawEdgeCosts(true);
					ge->SetDrawNodeLabels(true);
				}
			}
		}
			break;
		case '[':
		{
			if (numExampleNodes > 4)
			{
				numExampleNodes--;
				BuildGraphFromPuzzle(windowID, kNoModifier, graphType);
			}
			if (numExampleNodes <= 8)
			{
				ge->SetNodeScale(1.75);
				ge->SetDrawEdgeCosts(true);
				ge->SetDrawNodeLabels(true);
			}
			else {
				ge->SetNodeScale(2.5);
				ge->SetDrawEdgeCosts(false);
				ge->SetDrawNodeLabels(false);
			}
		}
			break;
		case '|':
		{
			BuildGraphFromPuzzle(windowID, kNoModifier, graphType);
//			name[0] = 'a';
//			g->Reset();
//			//te.AddLine("Current mode: add nodes");
//			m = kFindPath;
//			path.resize(0);
//			running = false;
		}
			break;
		case 'i':
			if (m == kFindPathIbex)
			{
				m = kFindPath;
			}
			else if (m == kFindPath)
			{
				m = kFindPathIbex;
			}
			ShowSearchInfo();
			break;
//		case 'w':
//			if (weight > 0.5)
//				weight = 0.0;
//			else
//				weight = 1.0;
//			m = kFindPath;
//			astar.SetWeight(weight);
//			astar.SetUseBPMX(BPMX?1:0);
//			astar.InitializeSearch(ge, astar.start, astar.goal, path);
//			ShowSearchInfo();
//
//			running = true;
//			break;
		case 'r':
			recording = !recording;
			break;
//		case '0':
		case 'x':
			BPMX = !BPMX;
			if (running)
			{
				ShowSearchInfo();
			}
			break;
		case 'p':
			paused = !paused;
			break;
		case 'o':
		{
			if (running && path.size() == 0)
			{
				if (m == kFindPath)
				{
					if (astar.GetNumOpenItems() > 0)
						printf("%1.2f\n", astar.GetOpenItem(0).g+astar.GetOpenItem(0).h);
					if (astar.DoSingleSearchStep(path))
						paused = true;
				}
				else if (m == kFindPathIbex)
				{
					if (ibex.DoSingleSearchStep(path))
						paused = true;
				}
				ShowSearchInfo();
			}
		}
			break;
		case '?':
		{
//			te.AddLine("Help:");
//			te.AddLine("-----");
//			te.AddLine("Press '[' and ']' to switch between modes.");
//			te.AddLine("Add nodes: click to add a node to the graph");
//			te.AddLine("Add edges: drag between nodes to add an edge");
//			te.AddLine("           Press 1-9 to change edge weight");
//			te.AddLine("Move nodes: drag to move node locations");
//			te.AddLine("Find Path: Drag to find path between nodes");
//			te.AddLine("           'o' to step pathfinding forward");
//			te.AddLine("           Red nodes: closed list");
//			te.AddLine("           Green nodes: open list");
//			te.AddLine("           Yellow node: next on open list");
//			te.AddLine("           Pink node: goal state");
		}
			break;
		case 's':
//			SaveGraph("save.graph");
		{
			Graphics::Display tmp;
			//tmp.FillRect({-1, -1, 1, 1}, Colors::black);
			ge->SetColor(Colors::black);
			ge->Draw(tmp);
			for (int x = 0; x < g->GetNumNodes(); x++)
			{
				ge->SetColor(Colors::black);
				ge->Draw(tmp, x);
				ge->SetColor(Colors::red);
				ge->DrawStateLabel(tmp, x, MyToString(g->GetNode(x)->GetLabelL(kHeuristic)).c_str());
				//ge->GLLabelState(x, MyToString(g->GetNode(x)->GetLabelL(kHeuristic)).c_str());
			}

			MakeSVG(tmp, "/Users/nathanst/Desktop/inc.svg", 400, 400, 0, "", false);
		}
			break;
		case 'l':
//			LoadGraph("save.graph");
			break;
		default:
			break;
	}
	
}

void BuildGraphFromPuzzle(unsigned long windowID, tKeyboardModifier mod, char key)
{
	if (key == 'w') // weighted inconsistency
	{
		graphType = key;
		delete g;
		delete ge;
//		double weight = 3;
		weight = 3;

		g = GraphInconsistencyExamples::GetWeightedInconsistency(weight, numExampleNodes);
		ge = new GraphEnvironment(g);
		ge->SetDirected(false);
		ge->SetDrawEdgeCosts(true);
		ge->SetDrawNodeLabels(true);
		ge->SetIntegerEdgeCosts(false);
		from = 0;
		to = g->GetNumNodes()-1;
		astar.SetWeight(3);
		astar.SetUseBPMX(0);
		astar.InitializeSearch(ge, from, to, path);
		m = kFindPath;
		ShowSearchInfo();
		running = true;
		paused = true;
		path.clear();
		if (numExampleNodes > 8)
		{
			ge->SetNodeScale(2.5);
			ge->SetDrawEdgeCosts(false);
			ge->SetDrawNodeLabels(false);
		}
		else {
			ge->SetNodeScale(1.75);
			ge->SetDrawEdgeCosts(true);
			ge->SetDrawNodeLabels(true);
		}
	}
	if (key == 'a')
	{
		graphType = key;
		delete g;
		delete ge;
		
		g = GraphInconsistencyExamples::GetPolyGraph(numExampleNodes);
		ge = new GraphEnvironment(g);
		ge->SetDirected(true);
		ge->SetDrawEdgeCosts(true);
		ge->SetDrawNodeLabels(true);
		ge->SetIntegerEdgeCosts(true);
	
		from = 0;
		to = g->GetNumNodes()-1;
		weight = 1.0;
		astar.SetWeight(weight);
		astar.SetUseBPMX(BPMX?1:0);
		astar.InitializeSearch(ge, from, to, path);
		ibex.InitializeSearch(ge, from, to, &h, path);
		ShowSearchInfo();
		running = true;
		paused = true;
		path.clear();
		if (numExampleNodes > 8)
		{
			ge->SetNodeScale(2.5);
			ge->SetDrawEdgeCosts(false);
			ge->SetDrawNodeLabels(false);
		}
		else {
			ge->SetNodeScale(1.75);
			ge->SetDrawEdgeCosts(true);
			ge->SetDrawNodeLabels(true);
		}
	}
	if (key == 'b' || key == 'c')
	{
		graphType = key;
		delete g;
		delete ge;
		
		if (key == 'b')
			g = GraphInconsistencyExamples::GetExpoGraphA(numExampleNodes);
		else
			g = GraphInconsistencyExamples::GetExpoGraphB(numExampleNodes);
		ge = new GraphEnvironment(g);
		ge->SetDirected(true);
		ge->SetDrawEdgeCosts(true);
		ge->SetDrawNodeLabels(true);
		ge->SetIntegerEdgeCosts(true);
		to = 0;
		from = g->GetNumNodes()-1;
		weight = 1.0;
		astar.SetWeight(weight);
		astar.SetUseBPMX(BPMX?1:0);
		astar.InitializeSearch(ge, from, to, path);
		ibex.InitializeSearch(ge, from, to, &h, path);
		ShowSearchInfo();
		running = true;
		paused = true;
		path.clear();
		if (numExampleNodes > 8)
		{
			ge->SetNodeScale(2.5);
			ge->SetDrawEdgeCosts(false);
			ge->SetDrawNodeLabels(false);
		}
		else {
			ge->SetNodeScale(1.75);
			ge->SetDrawEdgeCosts(true);
			ge->SetDrawNodeLabels(true);
		}
	}
}
#include <iomanip>
#include <sstream>

std::string MyToString(double val)
{
	std::stringstream ss;
	ss << std::fixed << std::setprecision(0) << val;
	return ss.str();
}

void ShowSearchInfo()
{
	const int colWidth = 24;
	std::string s;
	te.Clear();
	if (m == kFindPathIbex)
	{
		s = "--> IBEX (BGS)";
		s += " Searching from ";
	}
	else {
		s = "--> A*";
		if (BPMX)
			s += "(BPMX)";
		else
			s += "(No BPMX)";
		s += " Searching from ";
	}

	s +=g->GetNode(astar.start)->GetName();
	s +=" to ";
	s += g->GetNode(astar.goal)->GetName();
	s += " <--";
	//printf("Start: %d, goal %d. total (%d)\n", astar.start, astar.goal, g->GetNumNodes());

	te.AddLine(s.c_str());
	
	//te.AddLine("Press 'o' to advance search.");
	if (numExampleNodes <= 8)
	{
		for (int x = 0; x < g->GetNumNodes(); x++)
		{
			double gcost;
			s = g->GetNode(x)->GetName();
			switch (astar.GetStateLocation(x))
			{
				case kClosedList:
				{
					s += ": Closed  (g: ";
					astar.GetClosedListGCost(x, gcost);
					s += MyToString(gcost);
					s += ", h: ";
					s += MyToString(h.HCost(x, astar.goal));
					s += ")";
				}
					break;
				case kOpenList:
				{
					s += ": Open    (g: ";
					astar.GetOpenListGCost(x, gcost);
					s += MyToString(gcost);
					s += ", h: ";
					s += MyToString(h.HCost(x, astar.goal));
					s += ")";
				}
					break;
					
				case kNotFound:
					s += ": Ungenerated (h: ";
					s += MyToString(h.HCost(x, astar.goal));
					s += ")";
					break;
			}
			
			te.AddLine(s.c_str());
		}
		
		te.AddLine("");
		te.AddLine("Open List:");
		size_t length = strlen(te.GetLastLine());
		for (int x = length; x < colWidth; x++)
			te.AppendToLine(" ");
		te.AppendToLine("Closed List:");
		
		std::vector<std::string> open, closed;
		
		for (int x = 0; x < astar.GetNumOpenItems(); x++)
		{
			auto item = astar.GetOpenItem(x);
			s = g->GetNode(item.data)->GetName();
			s += ": ";
			s += MyToString(item.g+item.h);
			s += "=";
			s += MyToString(item.g);
			s += "+";
			s += MyToString(item.h);
			s += " p: ";
			s += g->GetNode(astar.GetItem(item.parentID).data)->GetName();
			//te.AddLine(s.c_str());
			open.push_back(s);
		}
		
		for (int x = 0; x < astar.GetNumItems(); x++)
		{
			auto item = astar.GetItem(x);
			if (item.where == kClosedList)
			{
				s = g->GetNode(item.data)->GetName();
				s += ": ";
				s += MyToString(item.g+item.h);
				s += "=";
				s += MyToString(item.g);
				s += "+";
				s += MyToString(item.h);
				s += " p: ";
				s += g->GetNode(astar.GetItem(item.parentID).data)->GetName();
				//te.AddLine(s.c_str());
				closed.push_back(s);
			}
		}
		for (size_t x = 0; x < open.size(); x++)
		{
			te.AddLine(open[x].c_str());
			for (size_t t = open[x].length(); t < colWidth; t++)
				te.AppendToLine(" ");
			if (x < closed.size())
				te.AppendToLine(closed[x].c_str());
		}
		for (size_t x = open.size(); x < closed.size(); x++)
		{
			te.AddLine("                        ");
			te.AppendToLine(closed[x].c_str());
		}

	}
}

double distsquared(unsigned long node, point3d loc)
{
	double dx = g->GetNode(node)->GetLabelF(GraphSearchConstants::kXCoordinate);
	double dy = g->GetNode(node)->GetLabelF(GraphSearchConstants::kYCoordinate);

	return (dx-loc.x)*(dx-loc.x) + (dy-loc.y)*(dy-loc.y);
}

double distance(unsigned long n1, unsigned long n2)
{
	double dx1 = g->GetNode(n1)->GetLabelF(GraphSearchConstants::kXCoordinate);
	double dy1 = g->GetNode(n1)->GetLabelF(GraphSearchConstants::kYCoordinate);

	double dx2 = g->GetNode(n2)->GetLabelF(GraphSearchConstants::kXCoordinate);
	double dy2 = g->GetNode(n2)->GetLabelF(GraphSearchConstants::kYCoordinate);

	return sqrt((dx1-dx2)*(dx1-dx2)+(dy1-dy2)*(dy1-dy2));
}

node *FindClosestNode(Graph *gr, point3d loc)
{
	if (gr->GetNumNodes() == 0)
		return 0;
	unsigned long best = 0;
	double dist = distsquared(0, loc);
	for (unsigned long x = 1; x < gr->GetNumNodes(); x++)
	{
		if (fless(distsquared(x, loc), dist))
		{
			dist = distsquared(x, loc);
			best = x;
		}
	}
	return gr->GetNode(best);
}


bool MyClickHandler(unsigned long , int windowX, int windowY, point3d loc, tButtonType button, tMouseEventType mType)
{
	if (mType == kMouseDown)
	{
		switch (button)
		{
			case kRightButton: printf("Right button\n"); break;
			case kLeftButton: printf("Left button\n"); break;
			case kMiddleButton: printf("Middle button\n"); break;
		}
	}
	if (button != kLeftButton)
		return false;
	switch (mType)
	{
		case kMouseDown:
		{
			printf("Hit (%f, %f, %f)\n", loc.x, loc.y, loc.z);
			if (m == kAddNodes)
			{
				if (loc.x > 1 || loc.x < -1 || loc.y > 1 || loc.y < -1)
					return false;
				node *n = new node(name);
				name[0]++;
				g->AddNode(n);
				n->SetLabelF(GraphSearchConstants::kXCoordinate, loc.x);
				n->SetLabelF(GraphSearchConstants::kYCoordinate, loc.y);
				n->SetLabelF(GraphSearchConstants::kZCoordinate, 0);
				printf("Added node %d to graph\n", g->GetNumNodes());
			}
			if (m == kAddEdges || m == kFindPath)
			{
				if (g->GetNumNodes() > 0)
				{
					from = to = FindClosestNode(g, loc)->GetNum();
				}
			}
			if (m == kMoveNodes)
			{
				from = to = FindClosestNode(g, loc)->GetNum();

				if (loc.x > 1) loc.x = 1;
				if (loc.x < -1) loc.x = -1;
				if (loc.y > 1) loc.y = 1;
				if (loc.y < -1) loc.y = -1;
				g->GetNode(from)->SetLabelF(GraphSearchConstants::kXCoordinate, loc.x);
				g->GetNode(from)->SetLabelF(GraphSearchConstants::kYCoordinate, loc.y);
				g->GetNode(from)->SetLabelF(GraphSearchConstants::kZCoordinate, 0);
				edge_iterator i = g->GetNode(from)->getEdgeIter();
				for (edge *e = g->GetNode(from)->edgeIterNext(i); e; e = g->GetNode(from)->edgeIterNext(i))
					e->setWeight(distance(e->getFrom(), e->getTo()));
			}
			return true;
		}
		case kMouseDrag:
		{
			if (m == kAddEdges || m == kFindPath)
			{
				if (g->GetNumNodes() > 0)
					to = FindClosestNode(g, loc)->GetNum();
			}
			if (m == kMoveNodes)
			{
				if (loc.x > 1) loc.x = 1;
				if (loc.x < -1) loc.x = -1;
				if (loc.y > 1) loc.y = 1;
				if (loc.y < -1) loc.y = -1;
				g->GetNode(from)->SetLabelF(GraphSearchConstants::kXCoordinate, loc.x);
				g->GetNode(from)->SetLabelF(GraphSearchConstants::kYCoordinate, loc.y);
				g->GetNode(from)->SetLabelF(GraphSearchConstants::kZCoordinate, 0);
				edge_iterator i = g->GetNode(from)->getEdgeIter();
				for (edge *e = g->GetNode(from)->edgeIterNext(i); e; e = g->GetNode(from)->edgeIterNext(i))
					e->setWeight(distance(e->getFrom(), e->getTo()));
			}
			return true;
		}
		case kMouseUp:
		{
			printf("UnHit at (%f, %f, %f)\n", loc.x, loc.y, loc.z);
			if (m == kAddEdges)
			{
				to = FindClosestNode(g, loc)->GetNum();
				if (from != to)
				{
					edge *e;
					if ((e = g->FindEdge(from, to)) != 0)
					{
						//e->setWeight(distance(from, to));
					}
					else {
						g->AddEdge(new edge(from, to, distance(from, to)));
					}
				}
			}
			if (m == kFindPath)
			{
				if (g->GetNumNodes() > 0)
					to = FindClosestNode(g, loc)->GetNum();
				if (from != to)
				{
					weight = 1.0;
					astar.SetWeight(weight);
					astar.SetUseBPMX(BPMX?1:0);
					astar.InitializeSearch(ge, from, to, path);
					ShowSearchInfo();
					
					running = true;
				}
			}
			if (m == kFindPathIbex)
			{
				if (g->GetNumNodes() > 0)
					to = FindClosestNode(g, loc)->GetNum();
				if (from != to)
				{
					ibex.InitializeSearch(ge, from, to, &h, path);
					ShowSearchInfo();
					running = true;
				}
			}

			if (m == kMoveNodes)
			{
				if (loc.x > 1) loc.x = 1;
				if (loc.x < -1) loc.x = -1;
				if (loc.y > 1) loc.y = 1;
				if (loc.y < -1) loc.y = -1;

				g->GetNode(from)->SetLabelF(GraphSearchConstants::kXCoordinate, loc.x);
				g->GetNode(from)->SetLabelF(GraphSearchConstants::kYCoordinate, loc.y);
				g->GetNode(from)->SetLabelF(GraphSearchConstants::kZCoordinate, 0);
				edge_iterator i = g->GetNode(from)->getEdgeIter();
				for (edge *e = g->GetNode(from)->edgeIterNext(i); e; e = g->GetNode(from)->edgeIterNext(i))
					e->setWeight(distance(e->getFrom(), e->getTo()));
			}
			from = to = -1;
			return true;
		}
	}
	return false;
}

void SaveGraph(const char *file)
{
	FILE *f = fopen(file, "w+");
	if (f)
	{
		fprintf(f, "%d %d\n", g->GetNumNodes(), g->GetNumEdges());
		for (int x = 0; x < g->GetNumNodes(); x++)
		{
			fprintf(f, "%d %f %f %f %s\n", x,
					g->GetNode(x)->GetLabelF(GraphSearchConstants::kXCoordinate),
					g->GetNode(x)->GetLabelF(GraphSearchConstants::kYCoordinate),
					g->GetNode(x)->GetLabelF(GraphSearchConstants::kZCoordinate),
					g->GetNode(x)->GetName()
					);
		}
		edge_iterator ei = g->getEdgeIter();
		while (1)
		{
			edge *e = (edge *)g->edgeIterNext(ei);
			if (e)
				fprintf(f, "%d %d %f\n", e->getFrom(), e->getTo(), e->GetWeight());
			else
				break;
		}
		fclose(f);
	}
}

void LoadGraph(const char *file)
{
	g->Reset();
	char nodeName[255];
	FILE *f = fopen(file, "r");
	if (f)
	{
		int numNodes, numEdges;
		fscanf(f, "%d %d\n", &numNodes, &numEdges);
		for (int n = 0; n < numNodes; n++)
		{
			int which;
			float x, y, z;
			fscanf(f, "%d %f %f %f %s\n", &which, &x, &y, &z, nodeName);
			assert(which == n);
			node *next = new node(nodeName);
			next->SetLabelF(GraphSearchConstants::kXCoordinate, x);
			next->SetLabelF(GraphSearchConstants::kYCoordinate, y);
			next->SetLabelF(GraphSearchConstants::kZCoordinate, z);
			g->AddNode(next);
		}
		for (int e = 0; e < numEdges; e++)
		{
			int from, to;
			float weight;
			fscanf(f, "%d %d %f", &from, &to, &weight);
			g->AddEdge(new edge(from, to, weight));
		}
		fclose(f);
	}
}
