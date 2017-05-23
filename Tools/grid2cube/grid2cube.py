#!/usr/bin/python
# -*- coding: utf-8 -*-

#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#                                                                      *
# Copyright (C) 2013, Ignacio Fdez. Galvan                             *
#               2008,2009, Neil Martinsen-Burrell (FortranFile)        *
#***********************************************************************
# ==============================================================================
# Python script for converting Molcas "grid" files into Gaussian "cube" format
#
# For reading binary grid files, it uses part of the FortranFile library
# (see below)
#
# Last modified: 2013 August 23
#            by: Ignacio Fdez. Galván
# ==============================================================================

import sys, re, os

def usage():
  print ""
  print "Convert a Molcas grid file into a Gaussian cube file."
  print "The input file can be ASCII or binary (non-packed),"
  print "the generated output file will be ASCII."
  print ""
  print "USAGE:"
  print "  {0} input_file output_file".format( os.path.basename(__file__) )
  print ""
  print "If there are several grids in the input file,"
  print "the user will be asked which one to convert."
  print ""

# Define symbol list, to get atomic numbers
symbol = [                                              "X",\
  "H",  "HE", "LI", "BE", "B",  "C",  "N",  "O",  "F",  "NE",
  "NA", "MG", "AL", "SI", "P",  "S",  "CL", "AR", "K",  "CA",
  "SC", "TI", "V",  "CR", "MN", "FE", "CO", "NI", "CU", "ZN",
  "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y",  "ZR",
  "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN",
  "SB", "TE", "I",  "XE", "CS", "BA", "LA", "CE", "PR", "ND",
  "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB",
  "LU", "HF", "TA", "W",  "RE", "OS", "IR", "PT", "AU", "HG",
  "TL", "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC", "TH",
  "PA", "U",  "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM",
  "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS",
  "RG", "CN", "UUT", "FL", "UUP", "LV", "UUS", "UUO"
]

# First argument is a grid input file
try:
  grid_input = sys.argv[1]
except IndexError:
  usage()
  sys.exit("** Missing input grid file")

# Second argument is a cube output file
try:
  cube_output = sys.argv[2]
except IndexError:
  usage()
  sys.exit("** Missing output cube file")

#===============================================================================
# THE FOLLOWING IS A MINIMIZED VERSION OF FortranFile
#    ( http://scipy-cookbook.readthedocs.io/items/FortranIO_FortranFile.html )
#
# Copyright 2008, 2009 Neil Martinsen-Burrell
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import struct, numpy

class FortranFile(file):

  def _get_header_length(self):
    return struct.calcsize(self._header_prec)

  _header_length = property(fget=_get_header_length)

  def _set_header_prec(self, prec):
    if prec in 'hilq':
      self._header_prec = prec
    else:
      raise ValueError('Cannot set header precision')

  def _get_header_prec(self):
    return self._header_prec

  HEADER_PREC = property(fset=_set_header_prec,
                         fget=_get_header_prec,
                         doc="Possible header precisions are 'h', 'i', 'l', 'q'"
                        )

  def __init__(self, fname, endian='@', header_prec='i', *args, **kwargs):
    file.__init__(self, fname, *args, **kwargs)
    self.ENDIAN = endian
    self.HEADER_PREC = header_prec

  def _read_exactly(self, num_bytes):
    data = ''
    while True:
      l = len(data)
      if l == num_bytes:
        return data
      else:
        read_data = self.read(num_bytes - l)
      if read_data == '':
        raise IOError('Could not read enough data.'
                      '  Wanted %d bytes, got %d.' % (num_bytes, l))
      data += read_data

  def _read_check(self):
    return struct.unpack(self.ENDIAN+self.HEADER_PREC,
                         self._read_exactly(self._header_length)
                        )[0]

  def readRecord(self):
    l = self._read_check()
    data_str = self._read_exactly(l)
    check_size = self._read_check()
    if check_size != l:
      raise IOError('Error reading record from data file')
    return data_str

  _real_precisions = 'df'

  def readReals(self, prec='f'):
    _numpy_precisions = {'d': numpy.float64,
                         'f': numpy.float32
                        }
    if prec not in self._real_precisions:
      raise ValueError('Not an appropriate precision')
    data_str = self.readRecord()
    num = len(data_str)/struct.calcsize(prec)
    numbers =struct.unpack(self.ENDIAN+str(num)+prec,data_str)
    return numpy.array(numbers, dtype=_numpy_precisions[prec])

# HERE ENDS THE FortranFile PART
#===============================================================================

#===============================================================================
# READ INPUT FILE
#===============================================================================

# Detect ASCII or binary file
binary = False
file_grid = open(grid_input, "r")
if (file_grid.next()[0] != "0"):
  file_grid.close()
  file_grid = FortranFile(grid_input)
  binary = True

# Read only the header (until the first "Title=" is found)
grid = []
while True:
  # Get next line/record in binary or ASCII mode
  if (binary):
    line = file_grid.readRecord()
  else:
    line = file_grid.next()

  # Number of atoms, the geometry follows
  match = re.match("Natom=\s*(\d+)", line)
  if (match):
    natoms = int(match.group(1))
    atom = []
    for i in range(natoms):
      if (binary):
        tmp = file_grid.readRecord().split()
      else:
        tmp = file_grid.next().split()
      atom.append( { "name": tmp[0], "x": float(tmp[1]), "y": float(tmp[2]), "z": float(tmp[3]) } )
  # Number of data per block (per grid)
  match = re.match("Block_Size=(.*)", line)
  if (match):
    block_size = int(match.group(1))
  # Number of grid divisions per dimension
  match = re.match("Net=(.*)", line)
  if (match):
    tmp = match.group(1).split()
    net = map(int, tmp)
    # Number of points is one more in each dimension
    npt = map(lambda x: x+1, net)
  # Origin of the grid (smallest corner)
  match = re.match("Origin=(.*)", line)
  if (match):
    tmp = match.group(1).split()
    origin = map(float, tmp)
  # The three edges
  match = re.match("Axis_1=(.*)", line)
  if (match):
    tmp = match.group(1).split()
    axis_1 = map(float, tmp)
  match = re.match("Axis_2=(.*)", line)
  if (match):
    tmp = match.group(1).split()
    axis_2 = map(float, tmp)
  match = re.match("Axis_3=(.*)", line)
  if (match):
    tmp = match.group(1).split()
    axis_3 = map(float, tmp)
  # Name of all grids present in the file
  match = re.match("GridName=(.*)", line)
  if (match):
    grid.append(match.group(1).strip())
  # The first grid data starts, stop reading header
  match = re.match("Title=(.*)", line)
  if (match):
    break

#---------------------------------------
# Select a grid from the input file
#---------------------------------------

# Print a menu showing the grids found in the file
if (len(grid) == 1):
  ng = 1
if (len(grid) < 1):
  sys.exit("** No grids found in input file")
if (len(grid) > 1):
  print "Grids in input file:"
  for i in range(len(grid)):
    print "{0:3}: {1}".format( i+1, grid[i] )
  # Read the user selection
  try:
    ng = int(raw_input("Which one to convert? "))
  except ValueError:
    sys.exit("** Sorry, I don't understand")
  if (ng < 1 or ng > len(grid)):
    sys.exit("** Sorry, no such grid")

# Total number of grid points to read
np = npt[0]*npt[1]*npt[2]
data = []

# The data for each grid is written across several blocks separated by "Title="
i = 0
j = 0
file_grid.seek(0)
if (binary):
  line = file_grid.readRecord()
else:
  line = file_grid.next()
while True:
  # When a new block starts, step the grid count by 1
  # (and start again when the last grid is reached)
  match = re.match("Title=(.*)", line)
  if (match):
    i = (i % len(grid))+1 
  # Only read the data if this is the grid we want
  if (i == ng):
    if (binary):
      # Assume we are reading double precision
      line = file_grid.readReals(prec='d')
      data.extend(line)
      j += len(line)
    else:
      # Each block is at most of block_size length (the last one is shorter)
      for k in range(min(np-j, block_size)): 
        data.append( float(file_grid.next().split()[0]) )
      j += k+1
  if (j >= np): break
  if (binary):
    line = file_grid.readRecord()
  else:
    line = file_grid.next()

file_grid.close()

#===============================================================================
# WRITE OUTPUT FILE
#===============================================================================

# Set the delta values in each dimension
dx = map( lambda x: x/float(net[0]), axis_1 )
dy = map( lambda x: x/float(net[1]), axis_2 )
dz = map( lambda x: x/float(net[2]), axis_3 )
# Assign an atomic number for each atom
# (from symbol, which is assumed to be the first 1 or 2 letters)
for i in atom:
  match = re.match("([A-Za-z]{1,2})", i["name"])
  i["number"] = symbol.index(match.group(1).upper())

file_cube = open(cube_output, "w")

# Write the header
file_cube.write("File converted from MOLCAS grid format\n")
file_cube.write("Title = {0}\n".format( grid[ng-1] ))
# Grid properties
file_cube.write("{0:5}{1:12.6f}{2:12.6f}{3:12.6f}\n".format( natoms, *origin ))
file_cube.write("{0:5}{1:12.6f}{2:12.6f}{3:12.6f}\n".format( npt[0], *dx ))
file_cube.write("{0:5}{1:12.6f}{2:12.6f}{3:12.6f}\n".format( npt[1], *dy ))
file_cube.write("{0:5}{1:12.6f}{2:12.6f}{3:12.6f}\n".format( npt[2], *dz ))
# Molecular geometry
for i in atom:
  file_cube.write("{0:5}{1:12.6f}{2:12.6f}{3:12.6f}{4:12.6f}\n".format( i["number"], 0.0, i["x"], i["y"], i["z"] ))

# Write the grid data, nd values in each line, new line every npt[2] elements
k = 0
nd = 6 # default 6
for i in range(npt[0]):
  for j in range(npt[1]):
    # kk is the final element of this block (i,j pair)
    kk = k+npt[2]
    # end is the final relative element of this (previous) line
    end = 0
    while (k+end < kk):
      # next element to print is k+end
      k += end
      # final element of this line is either kk or k+nd
      end = min(kk-k, nd) 
      # write elements from k to k+end
      file_cube.write("".join([ "{0:13.5E}".format(num) for num in data[k:k+end] ] )+"\n")
    # update for the final line of the block
    k += end

file_cube.close()
