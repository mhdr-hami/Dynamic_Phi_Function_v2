#-----------------------------------------------------------------------------
# GNU Makefile for static libraries: project dependent part
#
# $Id: Makefile.prj.inc,v 1.2 2006/11/01 23:34:14 nathanst Exp $
# $Source: /usr/cvsroot/project_hog/build/gmake/graph/Makefile.prj.inc,v $
#-----------------------------------------------------------------------------

NAME = envutil
DBG_NAME = lib$(NAME).a
REL_NAME = lib$(NAME).a

ROOT = ../../..
VPATH = $(ROOT)

DBG_OBJDIR = $(ROOT)/objs_html/$(NAME)/debug
REL_OBJDIR = $(ROOT)/objs_html/$(NAME)/release
DBG_BINDIR = $(ROOT)/html/debug
REL_BINDIR = $(ROOT)/html/release

PROJ_CXXFLAGS = -I$(ROOT)/graph -I$(ROOT)/environments -I$(ROOT)/utils -I$(ROOT)/abstraction -I$(ROOT)/simulation -I$(ROOT)/graphalgorithms  -I$(ROOT)/generic -I$(ROOT)/algorithms -I$(ROOT)/search -I$(ROOT)/gui
PROJ_DBG_CXXFLAGS = $(PROJ_CXXFLAGS)
PROJ_REL_CXXFLAGS = $(PROJ_CXXFLAGS)

PROJ_DBG_LNFLAGS =
PROJ_REL_LNFLAGS =

PROJ_DBG_LIB =
PROJ_REL_LIB =

default : all

SRC_CPP = \
	envutil/VoxelTriangleExtraction.cpp \
	envutil/STPInstances.cpp \
	envutil/PancakeInstances.cpp \
	envutil/GraphInconsistencyInstances.cpp \
