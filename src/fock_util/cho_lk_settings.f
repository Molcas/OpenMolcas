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
* Copyright (C) 2013, Thomas Bondo Pedersen                            *
************************************************************************
C#define _DEBUGPRINT_
      Real*8 Function Cho_LK_ScreeningThreshold(delta)
C
C     Thomas Bondo Pedersen, May 2013.
C
C     Return the basic LK screening threshold. Input is a damping
C     parameter, f.ex. the max. Fock matrix element in the previous
C     iteration. Use delta<0.0d0 or delta>1.0d0 to avoid damping.
C
C     The Cholesky environment must have been set up prior to calling
C     this function [by calling Cho_X_Init(..)]
C
      Implicit None
      Real*8 delta
#include "cholesky.fh"
      Real*8 thr0
      Parameter (thr0=1.0d-6)
      Real*8 thr
      thr=min(ThrCom,thr0)
      If (delta.ge.0.0d0 .and. delta.le.1.0d0) Then
         thr=thr*delta
      End If
      Cho_LK_ScreeningThreshold=max(thr,1.0d-15)
#if defined (_DEBUGPRINT_)
      Write(6,'(1P,4(A,D15.6))')
     & 'ThrCom=',ThrCom,' thr0=',thr0,' delta=',delta,
     & ' Cho_LK_ScreeningThreshold=',Cho_LK_ScreeningThreshold
#endif
      End
      Integer Function Cho_LK_MaxVecPerBatch()
C
C     Thomas Bondo Pedersen, May 2013.
C
C     Return the max number of vectors to be handled in core in the LK
C     algorithm.
C
C     Implemented as a function to make all LK flavors (SCF, RASSCF,
C     RASSI) behave the same way. The value 25 has been chosen to mimic
C     the approximate number of vectors in each reduced set of full CD
C     with the one-step algorithm. This should make the LK algorithm
C     (nearly) independent of available memory.
C
      Implicit None
      Cho_LK_MaxVecPerBatch=25
#if defined (_DEBUGPRINT_)
      Write(6,'(A,I10)')
     & 'Cho_LK_MaxVecPerBatch=',Cho_LK_MaxVecPerBatch
#endif
      End
