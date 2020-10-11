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
* Copyright (C) 2012, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_SortCoefficients(StorageMode,AB,n,M,C)
C
C     Thomas Bondo Pedersen, September 2012.
C
C     Sort coefficients according to storage mode.
C       StorageMode=0: shell-pair blocks [nothing done]
C                  =1: canonical
C                  anything else is an error.
C
C     Can of course be used for 3-index integrals as well.
C
      Implicit None
      Integer StorageMode
      Integer AB
      Integer n, M
      Real*8  C(n,M)
#include "ldf_atom_pair_info.fh"
#include "WrkSpc.fh"

      Character*20 SecNam
      Parameter (SecNam='LDF_SortCoefficients')
      Character*53 text ! length=len(SecNam)+33

      Integer  LDF_nShell_Atom
      External LDF_nShell_Atom
#if defined (_DEBUGPRINT_)
      Logical  LDF_AtomPairInfoIsSet
      External LDF_AtomPairInfoIsSet
      Integer  LDF_AtomPair_DiagDim, LDF_nBasAux_Pair
      Integer  LDF_nBasAux_Pair_wLD
      External LDF_AtomPair_DiagDim, LDF_nBasAux_Pair
      External LDF_nBasAux_Pair_wLD
#endif

      Logical OffsetDefined

      Integer A, B
      Integer nSA, nSB
      Integer ip_iOff, l_iOff
      Integer ip_Scr, l_Scr
      Integer K

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      ! Check StorageMode
      If (StorageMode.eq.0) Then
         Return
      Else If (StorageMode.ne.1) Then
         Write(text,'(A,A,I4,A)') SecNam,': StorageMode',StorageMode,
     &                          ' not implemented'
         Call WarningMessage(2,text)
         Call LDF_Quit(1)
      End If

#if defined (_DEBUGPRINT_)
      ! Check
      If (.not.LDF_AtomPairInfoIsSet()) Then
         Call WarningMessage(2,SecNam//': atom pair info not set')
         Call LDF_Quit(1)
      End If
      If (AB.lt.1 .or. AB.gt.NumberOfAtomPairs) Then
         Call WarningMessage(2,SecNam//': AB out of bounds')
         Call LDF_Quit(1)
      End If
      If (n.ne.LDF_AtomPair_DiagDim(AB)) Then
         Call WarningMessage(2,
     &                     SecNam//': C array not properly dimensioned')
         Call LDF_Quit(1)
      End If
      If (M.ne.LDF_nBasAux_Pair(AB) .and. M.ne.LDF_nBasAux_Pair_wLD(AB))
     & Then
         Call WarningMessage(0,
     &              SecNam//': C array may not be properly dimensioned')
      End If
#endif

      ! Get atoms of atom pair AB
      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)

      ! Allocate shell-pair offset array
      nSA=LDF_nShell_Atom(A)
      nSB=LDF_nShell_Atom(B)
      l_iOff=nSA*nSB
      Call GetMem('SrtOff','Allo','Inte',ip_iOff,l_iOff)
      OffsetDefined=.False.

      ! Allocate scratch array
      l_Scr=n
      Call GetMem('SrtScr','Allo','Real',ip_Scr,l_Scr)

      ! Sort
      Do K=1,M
         Call LDF_SortCanonical(AB,l_Scr,Work(ip_Scr),
     &                          OffsetDefined,
     &                          nSA,nSB,iWork(ip_iOff),
     &                          n,C(1,K))
      End Do

      ! Deallocations
      Call GetMem('SrtScr','Free','Real',ip_Scr,l_Scr)
      Call GetMem('SrtOff','Free','Inte',ip_iOff,l_iOff)

      End
