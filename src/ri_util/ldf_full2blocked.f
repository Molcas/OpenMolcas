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
      Subroutine LDF_Full2Blocked(FullMatrix,Packed,ip_Blocks)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: reorder FullMatrix in atom pair blocks.
C              FullMatrix may be stored triangularly (Packed=.True.)
C              or as a quadratic array (Packed=.False.).
C              iWork(ip_Blocks-1+iAtomPair) must point to the
C              location in Work where the block of atom pair iAtomPair
C              is to be stored. Note that the blocks are always stored
C              rectangularly.
C
      Implicit None
      Real*8  FullMatrix(*)
      Logical Packed
      Integer ip_Blocks
#include "WrkSpc.fh"
#include "localdf_bas.fh"

#if defined (_DEBUG_)
      Character*16 SecNam
      Parameter (SecNam='LDF_Full2Blocked')

      Logical  LDF_AtomPairInfoIsSet
      External LDF_AtomPairInfoIsSet
#endif

      Integer  LDF_nAtomPair
      External LDF_nAtomPair

      Integer iAtomPair
      Integer ip_SB, l_SB
      Integer n, ip0, iS
      Integer ip

      Integer i
      Integer nBasSh
      nBasSh(i)=iWork(ip_nBasSh-1+i)

#if defined (_DEBUG_)
      ! Check that atom pair info is set
      If (.not.LDF_AtomPairInfoIsSet()) Then
         Call WarningMessage(2,SecNam//': atom pair info not set!')
         Call LDF_Quit(1)
      End If
#endif

      ! Allocate and set offset array to shell blocks
      l_SB=nShell_Valence
      Call GetMem('SB','Allo','Inte',ip_SB,l_SB)
      n=0
      ip0=ip_SB-1
      Do iS=1,nShell_Valence
         iWork(ip0+iS)=n
         n=n+nBasSh(iS)
      End Do
#if defined (_DEBUG_)
      If (n.ne.nBas_Valence) Then
         Call WarningMessage(2,SecNam//': nBas_Valence problem!')
         Call LDF_Quit(1)
      End If
#endif

      ! Block-by-block sorting
      Do iAtomPair=1,LDF_nAtomPair()
         ip=iWork(ip_Blocks-1+iAtomPair)
         Call LDF_Full2Block(iAtomPair,FullMatrix,Packed,iWork(ip_SB),
     &                       Work(ip))
      End Do

      ! Deallocation
      Call GetMem('SB','Free','Inte',ip_SB,l_SB)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Full2Block(iAtomPair,FullMatrix,Packed,iSB,Block)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: extract atom pair block from FullMatrix.
C
      Implicit None
      Integer iAtomPair
      Real*8  FullMatrix(*)
      Logical Packed
      Integer iSB(*)
      Real*8  Block(*)

      If (Packed) Then
         Call LDF_Full2Block_Packed(iAtomPair,FullMatrix,iSB,Block)
      Else
         Call LDF_Full2Block_Quadratic(iAtomPair,FullMatrix,iSB,Block)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Full2Block_Packed(iAtomPair,FullMatrix,iSB,Block)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: extract atom pair block from packed FullMatrix.
C
      Implicit None
      Integer iAtomPair
      Real*8  FullMatrix(*)
      Integer iSB(*)
      Real*8  Block(*)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nShell_Atom, LDF_lShell_Atom

      Integer iAtom, jAtom
      Integer nShell_iAtom, nShell_jAtom
      Integer iS, jS
      Integer ipi, ipj
      Integer iShell, jShell
      Integer iOffi, iOffj
      Integer j_
      Integer ij0
      Integer iOffB

      Integer i, j
      Integer nBasSh
      Integer AP_Atoms, iTri
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

      ! Get atoms of atom pair
      iAtom=AP_Atoms(1,iAtomPair)
      jAtom=AP_Atoms(2,iAtomPair)

      ! Get number of shells on atoms
      nShell_iAtom=LDF_nShell_Atom(iAtom)
      nShell_jAtom=LDF_nShell_Atom(jAtom)

      ! Get pointers to atomic shell lists
      ipi=LDF_lShell_Atom(iAtom)-1
      ipj=LDF_lShell_Atom(jAtom)-1

      ! Extract block
      iOffB=0
      Do jS=1,nShell_jAtom
         jShell=iWork(ipj+jS)
         iOffj=iSB(jShell)
         Do iS=1,nShell_iAtom
            iShell=iWork(ipi+iS)
            iOffi=iSB(iShell)
            Do j=1,nBasSh(jShell)
               j_=iOffj+j
               ij0=iOffB+nBasSh(iShell)*(j-1)
               Do i=1,nBasSh(iShell)
                  Block(ij0+i)=FullMatrix(iTri(iOffi+i,j_))
               End Do
            End Do
            iOffB=iOffB+nBasSh(iShell)*nBasSh(jShell)
         End Do
      End Do

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Full2Block_Quadratic(iAtomPair,FullMatrix,iSB,
     &                                    Block)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: extract atom pair block from quadratic FullMatrix.
C
      Implicit None
      Integer iAtomPair
      Real*8  FullMatrix(*)
      Integer iSB(*)
      Real*8  Block(*)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nShell_Atom, LDF_lShell_Atom

      Integer iAtom, jAtom
      Integer nShell_iAtom, nShell_jAtom
      Integer iS, jS
      Integer ipi, ipj
      Integer iShell, jShell
      Integer iOffi, iOffj
      Integer j_
      Integer ij0, ij0_
      Integer iOffB

      Integer i, j
      Integer nBasSh
      Integer AP_Atoms
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      ! Get atoms of atom pair
      iAtom=AP_Atoms(1,iAtomPair)
      jAtom=AP_Atoms(2,iAtomPair)

      ! Get number of shells on atoms
      nShell_iAtom=LDF_nShell_Atom(iAtom)
      nShell_jAtom=LDF_nShell_Atom(jAtom)

      ! Get pointers to atomic shell lists
      ipi=LDF_lShell_Atom(iAtom)-1
      ipj=LDF_lShell_Atom(jAtom)-1

      ! Extract block
      iOffB=0
      Do jS=1,nShell_jAtom
         jShell=iWork(ipj+jS)
         iOffj=iSB(jShell)
         Do iS=1,nShell_iAtom
            iShell=iWork(ipi+iS)
            iOffi=iSB(iShell)
            Do j=1,nBasSh(jShell)
               j_=iOffj+j
               ij0=iOffB+nBasSh(iShell)*(j-1)
               ij0_=nBas_Valence*(j_-1)+iOffi
               Do i=1,nBasSh(iShell)
                  Block(ij0+i)=FullMatrix(ij0_+i)
               End Do
            End Do
            iOffB=iOffB+nBasSh(iShell)*nBasSh(jShell)
         End Do
      End Do

      End
