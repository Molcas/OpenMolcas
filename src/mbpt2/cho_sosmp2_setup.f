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
* Copyright (C) 2007, Francesco Aquilante                              *
************************************************************************
      SubRoutine Cho_SOSmp2_Setup(irc)
C
C     Francesco Aquilante   May 2007.
C
C     Purpose: setup of SOS-MP2 program.
C
#include "implicit.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"
************************************************************************
      MulD2h(i,j)=iEor(i-1,j-1) + 1
************************************************************************

      irc = 0

C     Setup index arrays and counters.
C     --------------------------------

      If (DecoMP2 .and. ThrMP2.le.0.0D0) Then
         Call Get_dScalar('Cholesky Threshold',ThrMP2)
      End If

      Call ChoMP2_GetInf(nOrb,nOcc,nFro,nDel,nVir)
      iOcc(1) = 0
      iVir(1) = 0
      nOccT = nOcc(1)
      nVirT = nVir(1)
      Do iSym = 2,nSym
         iOcc(iSym) = nOccT
         iVir(iSym) = nVirT
         nOccT = nOccT + nOcc(iSym)
         nVirT = nVirT + nVir(iSym)
      End Do

      Do iSym = 1,nSym
         nT1am(iSym) = 0
         Do iSymi = 1,nSym
            iSyma = MulD2h(iSymi,iSym)
            iT1am(iSyma,iSymi) = nT1am(iSym)
            nT1am(iSym) = nT1am(iSym)
     &                  + nVir(iSyma)*nOcc(iSymi)
         End Do
      End Do

      Do iSym = 1,nSym
         nT1AOT(iSym) = 0
         Do iSymAl = 1,nSym
            iSymi = MulD2h(iSymAl,iSym)
            iT1AOT(iSymi,iSymAl) = nT1AOT(iSym)
            nT1AOT(iSym) = nT1AOT(iSym)
     &                   + nOcc(iSymi)*nBas(iSymAl)
         End Do
      End Do

      Do iSym = 1,nSym
         nAOVir(iSym) = 0
         Do iSyma = 1,nSym
            iSymAl = MulD2h(iSyma,iSym)
            iAOVir(iSymAl,iSyma) = nAOVir(iSym)
            nAOVir(iSym) = nAOVir(iSym)
     &                   + nBas(iSymAl)*nVir(iSyma)
         End Do
      End Do

      If (ChoAlg .eq. 2) Then
         Do iSym = 1,nSym
            nMatab(iSym) = 0
            Do iSymb = 1,nSym
               iSyma = MulD2h(iSymb,iSym)
               iMatab(iSyma,iSymb) = nMatab(iSym)
               nMatab(iSym) = nMatab(iSym) + nVir(iSyma)*nVir(iSymb)
            End Do
         End Do
      Else
         Call Cho_iZero(nMatab,8)
         Call Cho_iZero(iMatab,64)
      End If

C     If batching over occuped orbitals is forced by user, then
C        turn it Off !
C     -----------------------------------------------------------------

      ForceBatch = .false.

      nBatch = 1

C     Initialize file units.
C     ----------------------

      Do iSym = 1,nSym
         Do iTyp = 1,nTypF
            Call ChoMP2_OpenF(0,iTyp,iSym)
         End Do
      End Do

      End

************************************************************************
      SubRoutine Cho_SOSmp2_Setup_Prt(irc)
C
C     Francesco Aquilante  May 2007
C
C     Purpose: print setup for SOS-MP2.
C
#include "implicit.fh"
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"


      irc = 0

      Call Cho_Head('Cholesky SOS-MP2 Setup','=',80,6)
      Write(6,*)

      If (nBatch .gt. 1) Then
         Write(6,'(A,I6,A,I6,A)')
     &   'The list of',nOccT,' occupied orbitals has been split in',
     &   nBatch,' batches:'
         Write(6,*)'Batching is not allowed in SOS-MP2 : I stop here! '
         Call Abend()
      Else If (nBatch .eq. 1) Then
         Write(6,'(A,I6,A)')
     &   'The list of',nOccT,' occupied orbitals is not split:'
      Else
         Write(6,*) 'Oops, #batches over occupied orbitals ',
     &              'is non-positive: ',nBatch
         irc = -101
         Return
      End If

      Write(6,'(//,A)')
     & 'The following tasks will be performed:'
      Write(6,'(A)')
     & ' * AO-to-MO transformation of original Cholesky vectors.'
      If (DecoMP2) Then
         Write(6,'(A)')
     &   ' * Cholesky decomposition of M=(ai|bj)^2 matrix.'
      End If
      Write(6,*)
     & ' * Calculation of SOS-MP2 correlation energy.'

      Call xFlush(6)

      End
