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
* Copyright (C) 2008, Igor Schapiro                                    *
************************************************************************
C
C *********************************************************************
C *                                                                   *
C *  Writes out the forces and energies for Gromacs                   *
C *                                                                   *
C * 18/01/2008                                                        *
C * Igor Schapiro                                                     *
C *                                                                   *
C *********************************************************************

C   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8

      SUBROUTINE GROM(irc)
#include "warnings.fh"
#include "Molcas.fh"
#include "prgm.fh"
#include "stdalloc.fh"
      PARAMETER   (ROUTINE='GROM')
#include "MD.fh"
#include "WrkSpc.fh"
      External    IsFreeUnit
      INTEGER     natom,i,j,irc,file,IsFreeUnit
      CHARACTER   filname*80
      REAL*8, ALLOCATABLE ::     xyz(:),force(:)
      CHARACTER, ALLOCATABLE ::  atom(:)*2
*
      IF(IPRINT.EQ.INSANE) WRITE(6,*)' Entering ',ROUTINE
      CALL QENTER(ROUTINE)

      WRITE(6,*)'**** Writes out Forces and Energies for Gromacs ****'
*
      CALL DxRdNAtomStnd(natom)
      CALL mma_allocate(atom,natom)
      CALL mma_allocate(xyz,natom*3)
      CALL mma_allocate(force,natom*3)
*
C
C     Read atom, their coordinates and forces
C
      CALL DxRdStnd(natom,atom,xyz,force)
C
C     Write the energies and forces to file
C
      file=IsFreeUnit(81)
      filname='MOL2GROM'
      Call Molcas_Open(file,filname)
      WRITE(file,*) natom
      DO i=1, natom
         WRITE(file,'(3D20.10)') (force((i-1)*3+j),j=1,3)
      ENDDO
      CLOSE(file)
*
      CALL mma_deallocate(atom)
      CALL mma_deallocate(xyz)
      CALL mma_deallocate(force)
*
      irc=_RC_ALL_IS_WELL_
      CALL qExit(ROUTINE)
      RETURN
*
      END
