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
      Subroutine LDF_PrintBlockVector(VectorName,ip_Blocks)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: Print block vector.
C
      Implicit None
      Character*(*) VectorName
      Integer ip_Blocks
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nBasAux_Pair
      External LDF_nBasAux_Pair

      real*8 ddot_
      external ddot_

      Integer l_max
      Parameter (l_max=80)
      Character*(l_max) myName

      Integer l
      Integer iAtomPair
      Integer iAtom, jAtom
      Integer M
      Integer ip

      Real*8 tnorm, bnorm

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      l=min(len(VectorName),l_max)
      If (l.gt.0) Then
         Write(myName,'(A)') VectorName(1:l)
      Else
         Write(myName,'(A)') '<unknown> '
      End If

      tnorm=0.0d0
      Do iAtomPair=1,NumberOfAtomPairs
         iAtom=AP_Atoms(1,iAtomPair)
         jAtom=AP_Atoms(2,iAtomPair)
         M=LDF_nBasAux_Pair(iAtomPair)
         ip=iWork(ip_Blocks-1+iAtomPair)
         bnorm=dDot_(M,Work(ip),1,Work(ip),1)
         tnorm=tnorm+bnorm
         Write(6,'(/,A,A,I9,A,I9,1X,I9,A)')
     &   myName(1:l),' block',iAtomPair,' (Atoms:',iAtom,jAtom,')'
         Write(6,'(A,I9,A,1P,D15.6)')
     &   'Dimension:',M,'    Norm:',sqrt(bnorm)
         Call Cho_Output(Work(ip),1,1,1,M,1,M,1,6)
      End Do
      Write(6,'(/,A,A,1P,D15.6)')
     & myName(1:l),' total norm:',sqrt(tnorm)
      Call xFlush(6)

      End
