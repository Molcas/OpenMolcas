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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_AllocateAuxBasVector(Label,ip_V)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: Allocate an auxiliary basis vector,
C              such as the vector V(J)=sum_uv C(uv,J)*D(uv) used for
C              computing the Coulomb contribution to the Fock matrix).
C              On exit, ip_V points to pointers in iWork such that
C
C              ip=iWork(ip_V-1+iAtom)
C
C              is the pointer (in array Work) to the vector block for
C              atom iAtom. Two-center functions are placed at the end
C              such that
C
C              ip=iWork(ip_V-1+nAtom+iAtomPair)
C
C              is the pointer (in array Work) to the vector block for
C              atom pair iAtomPair.
C
C              Deallocation is done by calling subroutine
C              LDF_DeallocateBlockVector(Label,ip_V).
C
C              Label is used to label the allocation in GetMem.
C              (Any three-character string is accepted.)
C
      Implicit None
      Character*3 Label
      Integer ip_V
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nAtom, LDF_nBasAux_Atom
      External LDF_nAtom, LDF_nBasAux_Atom

      Integer l
      Integer ip0
      Integer iAtom, nAtom
      Integer ip
      Integer iAtomPair

      Character*8 myLabel

      Integer i, j
      Integer AP_2CFunctions
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      ! Get number of atoms
      nAtom=LDF_nAtom()

      ! Allocate pointer array
      Write(myLabel,'(A3,A5)') Label,'Blk_P'
      l=nAtom+NumberOfAtomPairs
      Call GetMem(myLabel,'Allo','Inte',ip_V,l)

      ! Set pointers relative to 0 (i.e. offsets)
      l=0
      ip0=ip_V-1
      Do iAtom=1,nAtom
         iWork(ip0+iAtom)=l
         l=l+LDF_nBasAux_Atom(iAtom)
      End Do
      ip0=ip0+nAtom
      Do iAtomPair=1,NumberOfAtomPairs
         iWork(ip0+iAtomPair)=l
         l=l+AP_2CFunctions(1,iAtomPair)
      End Do

      ! Allocate matrix blocks
      Write(myLabel,'(A3,A5)') Label,'Block'
      Call GetMem(myLabel,'Allo','Real',ip,l)

      ! Shift pointers so that they become relative to allocated memory
      ip0=ip_V-1
      Do i=1,nAtom+NumberOfAtomPairs
         iWork(ip0+i)=iWork(ip0+i)+ip
      End Do

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_DeallocateAuxBasVector(Label,ip_V)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: Deallocate aux bas vector (see LDF_AllocateAuxBasVector).
C
      Implicit None
      Character*3 Label
      Integer ip_V
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nAtom, LDF_nBasAux_Atom
      External LDF_nAtom, LDF_nBasAux_Atom

      Integer l
      Integer iAtom, nAtom
      Integer ip
      Integer iAtomPair

      Character*8 myLabel

      Integer i, j
      Integer AP_2CFunctions
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      ! Get number of atoms
      nAtom=LDF_nAtom()

      ! Get total length of block matrix
      l=0
      Do iAtom=1,nAtom
         l=l+LDF_nBasAux_Atom(iAtom)
      End Do
      Do iAtomPair=1,NumberOfAtomPairs
         l=l+AP_2CFunctions(1,iAtomPair)
      End Do

      ! Deallocate matrix blocks
      Write(myLabel,'(A3,A5)') Label,'Block'
      ip=iWork(ip_V)
      Call GetMem(myLabel,'Free','Real',ip,l)

      ! Deallocate pointer array
      Write(myLabel,'(A3,A5)') Label,'Blk_P'
      l=nAtom+NumberOfAtomPairs
      Call GetMem(myLabel,'Free','Inte',ip_V,l)

      End
