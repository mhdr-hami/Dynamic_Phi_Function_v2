/*
 *  $Id: sample.h
 *  hog2
 *
 *  Created by Nathan Sturtevant on 5/31/05.
 *  Modified by Nathan Sturtevant on 02/29/20.
 *
 * This file is part of HOG2. See https://github.com/nathansttt/hog2 for licensing information.
 *
 */

void MyWindowHandler(unsigned long windowID, tWindowEventType eType);
void MyFrameHandler(unsigned long windowID, unsigned int viewport, void *data);
void MyDisplayHandler(unsigned long windowID, tKeyboardModifier, char key);
void MyPancakeHandler(unsigned long windowID, tKeyboardModifier, char key);
void MyFlingHandler(unsigned long windowID, tKeyboardModifier, char key);
//void MyPathfindingKeyHandler(unsigned long windowID, tKeyboardModifier, char key);
//void MyRandomUnitKeyHandler(unsigned long windowID, tKeyboardModifier, char key);

void MyPDBKeyHandler(unsigned long windowID, tKeyboardModifier, char key);


int MyCLHandler(char *argument[], int maxNumArgs);
bool MyClickHandler(unsigned long windowID, int x, int y, point3d loc, tButtonType, tMouseEventType);
void InstallHandlers();
