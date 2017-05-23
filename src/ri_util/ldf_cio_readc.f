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
      Subroutine LDF_CIO_ReadC(iAtomPair,C,l_C)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: Read LDF fitting coefficients for atom pair iAtomPair.
C              The coefficients are read from memory (buffer) if
C              possible (though the caller need not worry about that).
C              Columns of C corresponding to linearly dependent
C              auxiliary (one-center) functions are excluded, i.e.
C              C contains only linearly independent columns.
C              For reading with linearly dependent columns,
C              use LDF_CIO_ReadC_WithLinDep.
C
      Implicit None
      Integer iAtomPair
      Integer l_C
      Real*8  C(l_C)
#include "WrkSpc.fh"
#include "ldf_cio.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nBas_Atom, LDF_nBasAux_Pair
      External LDF_nBas_Atom, LDF_nBasAux_Pair

      Integer iAddr
      Integer iAtom, jAtom
      Integer l

#if defined (_DEBUG_)
      Real*8 Byte
      Character*2 Unt
#endif

      Integer i, j
      Integer AP_DiskC, AP_Atoms
      AP_DiskC(i)=iWork(ip_AP_DiskC-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      If (Lu_LDFC.lt.1) Then
         Call WarningMessage(2,'LDF_CIO_ReadC: Lu_LDFC<1')
         Call LDF_Quit(1)
      End If

      iAtom=AP_Atoms(1,iAtomPair)
      jAtom=AP_Atoms(2,iAtomPair)
      l=LDF_nBas_Atom(iAtom)*LDF_nBas_Atom(jAtom)
     & *LDF_nBasAux_Pair(iAtomPair)
      If (l.gt.l_C) Then
         Call WarningMessage(2,
     &                    'LDF_CIO_ReadC: insufficient array dimension')
         Call LDF_Quit(1)
      End If

      If (iAtomPair.gt.LastAtomPair) Then ! read from disk
         iAddr=AP_DiskC(iAtomPair)
         Call dDAFile(Lu_LDFC,2,C,l,iAddr)
#if defined (_DEBUG_)
         Call Cho_Word2Byte(l,8,Byte,Unt)
         Write(6,'(A,I9,A,I9,A,F7.2,A,A)')
     &   'LDF_CIO_ReadC: atom pair=',iAtomPair,' Dim=',l,' (',Byte,Unt,
     &   ') coefficients read from disk'
         Call xFlush(6)
#endif
      Else ! copy from buffer
         iAddr=iWork(ip_LDFC_Blocks-1+iAtomPair)
         Call dCopy_(l,Work(iAddr),1,C,1)
#if defined (_DEBUG_)
         Call Cho_Word2Byte(l,8,Byte,Unt)
         Write(6,'(A,I9,A,I9,A,F7.2,A,A)')
     &   'LDF_CIO_ReadC: atom pair=',iAtomPair,' Dim=',l,' (',Byte,Unt,
     &   ') coefficients copied from buffer'
         Call xFlush(6)
#endif
      End If

      End
