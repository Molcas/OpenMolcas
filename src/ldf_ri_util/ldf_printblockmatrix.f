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
      Subroutine LDF_PrintBlockMatrix(MatrixName,ip_Blocks)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: Print block matrix.
C
      Implicit None
      Character*(*) MatrixName
      Integer ip_Blocks
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom

      real*8 ddot_
      external ddot_

      Integer l_max
      Parameter (l_max=80)
      Character*(l_max) myName

      Real*8 tnorm, bnorm

      Integer l
      Integer iAtomPair
      Integer iAtom, jAtom
      Integer nBas_i, nBas_j
      Integer ip
      Integer nShell_i, nShell_j
      Integer ipi, ipj
      Integer iS, jS
      Integer iShell, jShell

      Integer i, j
      Integer nBasSh
      Integer AP_Atoms
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      l=min(len(MatrixName),l_max)
      If (l.gt.0) Then
         Write(myName,'(A)') MatrixName(1:l)
      Else
         Write(myName,'(A)') '<unknown> '
      End If

      Call Cho_Head(myName(1:l)//' (blocked)','=',120,6)

      tnorm=0.0d0
      Do iAtomPair=1,NumberOfAtomPairs
         ip=iWork(ip_Blocks-1+iAtomPair)
         iAtom=AP_Atoms(1,iAtomPair)
         jAtom=AP_Atoms(2,iAtomPair)
         nBas_i=LDF_nBas_Atom(iAtom)
         nBas_j=LDF_nBas_Atom(jAtom)
         bnorm=dDot_(nBas_i*nBas_j,Work(ip),1,Work(ip),1)
         tnorm=tnorm+bnorm
         Write(6,'(/,A,A,I9,A,I9,1X,I9,A)')
     &   myName(1:l),' block',iAtomPair,' (Atoms:',iAtom,jAtom,')'
         Write(6,'(A,I9,A,I9,A,1P,D15.6)')
     &   'Block dimension:',nBas_i,' by',nBas_j,'    Block norm:',
     &   sqrt(bnorm)
         nShell_i=LDF_nShell_Atom(iAtom)
         nShell_j=LDF_nShell_Atom(jAtom)
         ipi=LDF_lShell_Atom(iAtom)-1
         ipj=LDF_lShell_Atom(jAtom)-1
         Do jS=1,nShell_j
            jShell=iWork(ipj+jS)
            Do iS=1,nShell_i
               iShell=iWork(ipi+iS)
               Write(6,'(/,A,A,I9,A,I9,1X,I9,A)')
     &         myName(1:l),' block',iAtomPair,' (Atoms:',iAtom,jAtom,')'
               Write(6,'(A,I9,1X,I9,A,I9,1X,I9,A,I9)')
     &         'Shells:',iS,jS,' (',iShell,jShell,')   Location:',ip
               Write(6,'(A,I9,A,I9,A,1P,D15.6)')
     &         'Dimension:',nBasSh(iShell),' by',nBasSh(jShell),
     &         '    Norm:',
     &         sqrt(ddot_(nBasSh(iShell)*nBasSh(jShell),Work(ip),1,
     &                                                 Work(ip),1))
               Call Cho_Output(Work(ip),1,nBasSh(iShell),
     &                         1,nBasSh(jShell),nBasSh(iShell),
     &                         nBasSh(jShell),1,6)
               ip=ip+nBasSh(iShell)*nBasSh(jShell)
            End Do
         End Do
      End Do
      Write(6,'(/,A,A,1P,D15.6)')
     & myName(1:l),' total norm:',sqrt(tnorm)
      Call xFlush(6)

      End
