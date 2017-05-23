************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Thomas Dresselhaus                                     *
************************************************************************
      subroutine embPotInit(preparingOutput)
************************************************************************
*                                                                      *
* Object: routine to read in an embedding potential on a grid from a   *
*         file.                                                        *
*                                                                      *
* Called from: Seward                                                  *
*              DrvMO                                                   *
*              EmbPotOutput                                            *
*                                                                      *
* Calling    : GetMem                                                  *
*                                                                      *
*     Author: Thomas Dresselhaus                                       *
*                                                                      *
************************************************************************

      Implicit None

****** Includes
#include "WrkSpc.fh"

      ! Holds the data which is read in in this subroutine
#include "embpotdata.fh"

****** Variables

      ! Unit of input file (the embedding potential)
      Integer iunit, isFreeUnit

      ! Index
      Integer i

      ! Switch to toggle whether only information relevant for the
      ! output needs to be read
      Logical preparingOutput

******
      embDebug=.false.

      ! Open the file
      iunit = isFreeUnit(1)
      if (preparingOutput.and.outGridPathGiven) then
       call molcas_open(iunit, outGridPath)
      else
       call molcas_open(iunit, embPotPath)
      end if

      ! TODO MAKE THIS READING PROCEDURE SAFE!!!

      ! Read in header of the file (just a line with one integer)
      read(iunit,*) nEmbGridPoints

      ! Allocate memory for the grid points, potential and weights
      Call GetMem('embG','ALLO','REAL',posEmbGridCoord,nEmbGridPoints*3)
      Call GetMem('embP','ALLO','REAL',posEmbPotVal,nEmbGridPoints)
      Call GetMem('embW','ALLO','REAL',posEmbWeight,nEmbGridPoints)

      ! Read in data
      do i=1, nEmbGridPoints
       if (preparingOutput.and.outGridPathGiven) then
        read(iunit, *) Work(posEmbGridCoord+i*3-3),
     &                 Work(posEmbGridCoord+i*3-2),
     &                 Work(posEmbGridCoord+i*3-1)
       else
        read(iunit, *) Work(posEmbGridCoord+i*3-3),
     &                 Work(posEmbGridCoord+i*3-2),
     &                 Work(posEmbGridCoord+i*3-1),
     &                 Work(posEmbWeight+i-1),
     &                 Work(posEmbPotVal+i-1)
       end if
      end do

      close(iunit)

      if (embDebug) then
        write(6,*) "---------------------------------------------------"
        write(6,*) "---------------------------------------------------"
        write(6,*) "Potential has been read in. Coords:"
        do i=1, nEmbGridPoints
         if (mod(i,587).eq.0) then
          write(6,*) i, Work(posEmbGridCoord+i*3-3),
     &               Work(posEmbGridCoord+i*3-2),
     &               Work(posEmbGridCoord+i*3-1)
         end if
        end do
        write(6,*) "---------------------------------------------------"
        write(6,*) "Potential value, weight"
        do i=1, nEmbGridPoints
         if (mod(i,587).eq.0) then
          write(6,*) i, Work(posEmbPotVal+i-1),
     &               Work(posEmbWeight+i-1)
         end if
        end do
        write(6,*) "---------------------------------------------------"
        write(6,*) "---------------------------------------------------"
      end if

      return
      end
