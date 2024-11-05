!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
Module Centra
# include "Molcas.fh"
! The paramaters defined in Molcas should be private
Private MaxBfn,MaxBfn_Aux, MxAO, mxAtom, mxroot, mxNemoAtom, Mxdbsc, lCache, mxact, mxina, mxbas, mxOrb, &
        mxSym, mxGAS, LENIN, LENIN1, LENIN2, LENIN3, LENIN4, LENIN5, LENIN6, LENIN8

! Note: MXATOM to be taken from Molcas.fh
INTEGER               NGROUP,IGROUP(8),NATOMS,COOR(3,MXATOM)
! Atom labels, 4 bytes each.
CHARACTER(LEN=LENIN) ATLBL(MXATOM)
End Module Centra
