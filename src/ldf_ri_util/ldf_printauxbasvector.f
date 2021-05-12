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
      Subroutine LDF_PrintAuxBasVector(VectorName,ip_V)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: Print aux bas vector.
C
      Implicit None
      Character*(*) VectorName
      Integer ip_V
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nAtom, LDF_nBasAux_Atom
      External LDF_nAtom, LDF_nBasAux_Atom

      real*8 ddot_
      external ddot_

      Integer l_max
      Parameter (l_max=80)
      Character*(l_max) myName

      Integer l
      Integer nAtom
      Integer iAtomPair
      Integer iAtom
      Integer M
      Integer ip

      Real*8 tnorm, bnorm

      Integer i, j
      Integer AP_2CFunctions
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      l=min(len(VectorName),l_max)
      If (l.gt.0) Then
         Write(myName,'(A)') VectorName(1:l)
      Else
         Write(myName,'(A)') '<unknown> '
      End If

      nAtom=LDF_nAtom()
      tnorm=0.0d0
      Do iAtom=1,nAtom
         M=LDF_nBasAux_Atom(iAtom)
         ip=iWork(ip_V-1+iAtom)
         bnorm=dDot_(M,Work(ip),1,Work(ip),1)
         tnorm=tnorm+bnorm
         Write(6,'(/,A,A,I9)')
     &   myName(1:l),' aux bas block for atom ',iAtom
         Write(6,'(A,I9,A,1P,D15.6)')
     &   'Dimension:',M,'    Norm:',sqrt(bnorm)
         Call Cho_Output(Work(ip),1,1,1,M,1,M,1,6)
      End Do
      Do iAtomPair=1,NumberOfAtomPairs
         M=AP_2CFunctions(1,iAtomPair)
         If (M.gt.0) Then
            ip=iWork(ip_V-1+nAtom+iAtomPair)
            bnorm=dDot_(M,Work(ip),1,Work(ip),1)
            tnorm=tnorm+bnorm
            Write(6,'(/,A,A,I9)')
     &      myName(1:l),' aux bas block for atom pair',iAtomPair
            Write(6,'(A,I9,A,1P,D15.6)')
     &      'Dimension:',M,'    Norm:',sqrt(bnorm)
            Call Cho_Output(Work(ip),1,1,1,M,1,M,1,6)
         End If
      End Do
      Write(6,'(/,A,A,1P,D15.6)')
     & myName(1:l),' total norm:',sqrt(tnorm)
      Call xFlush(6)

      End
