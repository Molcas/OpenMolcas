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
! module file for module lnklst (linked list stuff)
!
Module LnkLst
use MxDM, only: MxIter
Private
Public :: Debug_LnkLst, lLList, nLList, MAXnodes, NodSiz
Public :: SCF_V
Integer, Parameter :: NodSiz=6
Integer, Parameter :: MAXnodes=MxIter*5
Logical :: Debug_LnkLst
Integer :: lLList=0
Integer :: nLList(MAXnodes,0:NodSiz-1)

Type Vector
     Real*8, Allocatable :: A(:)
End Type Vector

Type (Vector) :: SCF_V(Maxnodes)

!----------------------------------------------------------------------*
!     Define linked lists for storage of vectors of subsequent iters   *
!----------------------------------------------------------------------*
!     LLGrad - linked list of gradients                                *
!     LLdGrd - linked list of gradient diffs                           *
!     LLDelt - linked list of Delta vectors                            *
!     LLy    - linked list of y-vectors                                *
!     LLx    - linked list of x orbital rotation parameter vectors     *
!              (only for QNR/DIIS combination)                         *
!----------------------------------------------------------------------*
!
Integer, Public ::  LLGrad,LLdGrd,LLDelt,LLy,LLx
Logical, Public ::  Init_LLs=.False.
End Module LnkLst

