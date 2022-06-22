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
      Subroutine LDF_AllocateBlockVector(Label,ip_Blocks)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: Allocate atom pair blocks of an auxiliary basis vector,
C              such as the vector V(J)=sum_uv C(uv,J)*D(uv) used for
C              computing the Coulomb contribution to the Fock matrix).
C              On exit, ip_Blocks points to pointers in iWork such that
C
C              ip=iWork(ip_Blocks-1+iAtomPair)
C
C              is the pointer (in array Work) to the vector block for
C              atom pair iAtomPair.
C
C              Deallocation is done by calling subroutine
C              LDF_DeallocateBlockVector(Label,ip_Blocks).
C
C              Label is used to label the allocation in GetMem.
C              (Any three-character string is accepted.)
C
      Implicit None
      Character*3 Label
      Integer ip_Blocks
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nBasAux_Pair
      External LDF_nBasAux_Pair

      Integer l
      Integer ip0
      Integer iAtomPair
      Integer ip

      Character*8 myLabel

      ! Allocate pointer array
      Write(myLabel,'(A3,A5)') Label,'Blk_P'
      l=NumberOfAtomPairs
      Call GetMem(myLabel,'Allo','Inte',ip_Blocks,l)
      ip0=ip_Blocks-1

      ! Set pointers relative to 0 (i.e. offsets)
      l=0
      Do iAtomPair=1,NumberOfAtomPairs
         iWork(ip0+iAtomPair)=l
         l=l+LDF_nBasAux_Pair(iAtomPair)
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
      Subroutine LDF_DeallocateBlockVector(Label,ip_Blocks)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: Deallocate block vector (see LDF_AllocateBlockVector).
C
      Implicit None
      Character*3 Label
      Integer ip_Blocks
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nBasAux_Pair
      External LDF_nBasAux_Pair

      Integer l
      Integer iAtomPair
      Integer ip

      Character*8 myLabel

      ! Get total length of block matrix
      l=0
      Do iAtomPair=1,NumberOfAtomPairs
         l=l+LDF_nBasAux_Pair(iAtomPair)
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
