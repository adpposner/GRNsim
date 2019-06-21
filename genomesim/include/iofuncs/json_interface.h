//json_interface.h - Interface for JSON I/O for sim
/*
    
    Copyright (C) 2018, Russell Posner

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
*/
#ifndef JSON_INTERFACE_H__
#define JSON_INTERFACE_H__
#include "cjson.h"
#include "../globals.h"


struct GenesList;

void writeGenesList_JSON(struct GenesList * g,const char * filename, SkipDisabled skipEnabled);
struct GenesList * readGenesList_JSON(const char * filename);

#endif