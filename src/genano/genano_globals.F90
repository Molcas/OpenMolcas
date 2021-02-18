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
!----------------------------------------------------------------------*
!                                                                      *
! This module contains various parameters and variables for GenAno.    *
!                                                                      *
! o MxLqn  - The highest l-q.n. that the program can handle.           *
! x MxSets - The highest number of orbitals sets that can be used.     *
! x MxBasX - The highest number of basis functions in any of the sets. *
! x MxS    - The maximum number of s-primitives.                       *
! x MxP    - The maximum number of p-primitives.                       *
! x MxD    - The maximum number of d-primitives.                       *
! x MxF    - The maximum number of f-primitives.                       *
! x MxG    - The maximum number of g-primitives.                       *
! x MxH    - The maximum number of h-primitives.                       *
! x MxI    - The maximum number of i-primitives.                       *
! x MxK    - The maximum number of k-primitives.                       *
! x MxDsym - The space allocated for l-q.n. blocked density matrices.  *
!                                                                      *
! o nPrim  - Number of primitives per l quantum number.                *
! o nCore  - Number of core orbitals per l quantum number.             *
!            Note: currently not used.                                 *
! o nSets  - Number of orbital sets to include in averaging.           *
! o kSet   - Set number counter.                                       *
! o nSym   - Number of symmetries blocks of current set.               *
! o nBas   - Number of basis functions per symmetry block of           *
!            current set.                                              *
! x nDsym  - Size of density matrix in l-q.n. blocked form stored as   *
!            lower triangles.                                          *
! o iProj  - Switch for performing projection.                         *
! o iSymBk - Pointer to start of l-q.n. blocks of density matrix.      *
! o kRfSet - Reference set from which the metric (overlap matrix) is   *
!            read.                                                     *
!                                                                      *
! o wSet   - Weight of the individual set in the averaging. Not        *
!            normalized if read from input.                            *
! o pDsym  - Partial l-q.n. blocked density matrix.                    *
! o tDsym  - Total l-q.n. blocked density matrix.                      *
! o thr    - Threshold for reporting NO's.                             *
! o Ssym   - Metric (overlap matrix) from reference set.               *
! o Cmo    - Natural orbitals.                                         *
! o Occ    - Occupation numbers from vector file.                      *
! o Eps    - Orbital energies                                          *
! o wc0    - Weight coefficient 0                                      *
! o wc1    - Weight coefficient 1                                      *
!                                                                      *
! o Center - Atomic label for which to generate ano's.                 *
! o BasName- Atomic labels and basis function label from one-electron  *
!            file.                                                     *
! o Title  - First title line from input.                              *
!                                                                      *
! o rowise - Switch for printing ano's as rows instead of columns.     *
! o lftdeg - Switch for lifting occupation number degeneracy           *
! o rydgen - Switch for generating rydberg orbitals                    *
!                                                                      *
! o symlab - Symmetry label information for basis functions generated  *
!            by seward.                                                *
!----------------------------------------------------------------------*

module Genano_globals

use Definitions, only: wp, iwp

implicit none
private

#include "Molcas.fh"
!--- MkS is set to 90 to circumvent a logical error in the code!
!#define MxS 90
!#define MxP 25
!#define MxD 20
!#define MxF 15
!#define MxG 12
!#define MxH 10
!#define MxI  8
!#define MxK  8
! MxDsym =   MxS*(MxS+1)/2 +  3*MxP*(MxP+1)/2 +  5*MxD*(MxD+1)/2 +  7*MxF*(MxF+1)/2 + &
!          9*MxG*(MxG+1)/2 + 11*MxH*(MxH+1)/2 + 13*MxI*(MxI+1)/2 + 15*MxK*(MxK+1)/2
integer(kind=iwp), parameter :: MxLqn = 7
integer(kind=iwp) :: nPrim(0:MxLqn), nCore(0:MxLqn), nSets, kSet, nSym, nBas(1:MxSym), nDsym, iProj, iSymBk((MxLqn+1)**2), kRfSet, &
                     isUHF
real(kind=wp) :: thr, wc0, wc1
real(kind=wp), allocatable :: wSet(:), pDsym(:), tDsym(:), Ssym(:), Cmo(:), Occ(:), Cmo2(:), Occ2(:), Eps(:)
logical(kind=iwp) :: rowise, lftdeg, rydgen
character(len=LenIn) :: Center
character(len=LenIn+8), allocatable :: BasName(:)
character(len=80) :: Title
character(len=8) :: symlab((MxLqn+1)**2)

public :: MxLqn, nPrim, nCore, nSets, kSet, nSym, nBas, nDsym, iProj, iSymBk, kRfSet, isUHF, wSet, pDsym, tDsym, thr, Ssym, Cmo, &
          Occ, Cmo2, Occ2, Eps, wc0, wc1, rowise, lftdeg, rydgen, LenIn, Center, BasName, Title, symlab

end module Genano_globals
