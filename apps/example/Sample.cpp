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
#include "Sample.h"

int main(int argc, char* argv[])
{
	InstallHandlers();
	RunHOGGUI(argc, argv);
	return 0;
}

/**
 * Allows you to install any keyboard handlers needed for program interaction.
 */
void InstallHandlers()
{
	InstallKeyboardHandler(MyDisplayHandler, "Numbers", "Will get events for numbers 0-9", kAnyModifier, '0', '9');
	
	// Command-line argument, argument with parameters, description
	InstallCommandLineHandler(MyCLHandler, "-test", "-test <id>", "Sample command-line handler.");
	
	// Gets window events (created/destroyed)
	InstallWindowHandler(MyWindowHandler);
	
	// Gets mouse events - handles all windows
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
		printf("Window %ld created\n", windowID);
		InstallFrameHandler(MyFrameHandler, windowID, 0);
		// Basic port layout
		// 1 port with (0,0) in the middle. Top left is (-1, -1), bottom right is (1, 1)
		SetNumPorts(windowID, 1);
		// Default text line at top of window is hidden
		setTextBufferVisibility(false);
	}
}

void MyFrameHandler(unsigned long windowID, unsigned int viewport, void *)
{
	// Get reference to display
	Graphics::Display &d = GetContext(windowID)->display;
	
	// Do any drawing here
	d.FillRect({{-1, -1}, {1, 1}}, Colors::darkblue);
}

int MyCLHandler(char *argument[], int maxNumArgs)
{
	if (strcmp(argument[0], "-test" ) == 0 )
	{
		printf("Running in test mode.");
		if (maxNumArgs <= 1)
			return 1;
		printf("Got argument '%s'\n", argument[1]);
		return 2;
	}
	return 0;
}

void MyDisplayHandler(unsigned long windowID, tKeyboardModifier mod, char key)
{
	switch (key)
	{
		case 'a':
			printf("Hit 'a'\n");
			break;
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
			printf("Hit %c\n", key);
			break;
		default:
			break;
	}
}

bool MyClickHandler(unsigned long windowID, int, int, point3d loc, tButtonType button, tMouseEventType mType)
{
	if (button == kLeftButton)
	{
		switch (mType)
		{
			case kMouseDown:
				printf("Mouse down\n");
				break;
			case kMouseDrag:
				printf("Mouse drag\n");
				break;
			case kMouseUp:
				printf("Mouse up\n");
				break;
			kMouseMove: // You have to request this event explicitly if you want to receive it
				break;
		}
		return true; // handled
	}
	return false; // not handled
}
