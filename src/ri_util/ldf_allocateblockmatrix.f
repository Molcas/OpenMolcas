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
      Subroutine LDF_AllocateBlockMatrix(Label,ip_Blocks)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: Allocate atom pair blocks of a quadratic matrix (e.g.
C              density matrix, Fock matrix). On exit, ip_Blocks points
C              to pointers in iWork such that
C
C              ip=iWork(ip_Blocks-1+iAtomPair)
C
C              is the pointer (in array Work) to the matrix block for
C              atom pair iAtomPair.
C
C              Deallocation is done by calling subroutine
C              LDF_DeallocateBlockMatrix(Label,ip_Blocks).
C
C              Label is used to label the allocation in GetMem.
C              (For example, Label could be 'Den' or 'Fck', but
C              any three-character string is accepted).
C
      Implicit None
      Character*3 Label
      Integer ip_Blocks
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom

      Integer l
      Integer ip0
      Integer iAtomPair
      Integer iAtom, jAtom
      Integer ip

      Character*8 myLabel

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      ! Allocate pointer array
      Write(myLabel,'(A3,A5)') Label,'Blk_P'
      l=NumberOfAtomPairs
      Call GetMem(myLabel,'Allo','Inte',ip_Blocks,l)
      ip0=ip_Blocks-1

      ! Set pointers relative to 0 (i.e. offsets)
      l=0
      Do iAtomPair=1,NumberOfAtomPairs
         iAtom=AP_Atoms(1,iAtomPair)
         jAtom=AP_Atoms(2,iAtomPair)
         iWork(ip0+iAtomPair)=l
         l=l+LDF_nBas_Atom(iAtom)*LDF_nBas_Atom(jAtom)
      End Do

      ! Allocate matrix blocks
      Write(myLabel,'(A3,A5)') Label,'Block'
      Call GetMem(myLabel,'Allo','Real',ip,l)

      ! Shift pointers so that they become relative to allocated memory
      Do iAtomPair=1,NumberOfAtomPairs
         iWork(ip0+iAtomPair)=iWork(ip0+iAtomPair)+ip
      End Do

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_DeallocateBlockMatrix(Label,ip_Blocks)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: Deallocate block matrix (see LDF_AllocateBlockMatrix).
C
      Implicit None
      Character*3 Label
      Integer ip_Blocks
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom

      Integer l
      Integer iAtomPair
      Integer iAtom, jAtom
      Integer ip

      Character*8 myLabel

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      ! Get total length of block matrix
      l=0
      Do iAtomPair=1,NumberOfAtomPairs
         iAtom=AP_Atoms(1,iAtomPair)
         jAtom=AP_Atoms(2,iAtomPair)
         l=l+LDF_nBas_Atom(iAtom)*LDF_nBas_Atom(jAtom)
      End Do

      ! Deallocate matrix blocks
      Write(myLabel,'(A3,A5)') Label,'Block'
      ip=iWork(ip_Blocks)
      Call GetMem(myLabel,'Free','Real',ip,l)

      ! Deallocate pointer array
      Write(myLabel,'(A3,A5)') Label,'Blk_P'
      l=NumberOfAtomPairs
      Call GetMem(myLabel,'Free','Inte',ip_Blocks,l)

      End
