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
* Copyright (C) 2019, Roland Lindh                                     *
************************************************************************
      Subroutine USOTRANS(USOR,USOI,NSS,
     &                    EigVec,MSTATE,
     &                    VSOR,VSOI)
      use rassi_global_arrays, only: JBNUM
      IMPLICIT Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Integer NSS, MSTATE
      Real*8 USOR(NSS,NSS), USOI(NSS,NSS), EigVec(MSTATE,MSTATE)
      Real*8 VSOR(NSS,NSS), VSOI(NSS,NSS)
      Integer, Allocatable:: MAPST(:)
      REAL*8 tmp_R, tmp_I
      Integer ISTATE, JOB, MPLET, MSPROJ, ISS, JSS, JSS_, KSS, KSS_
*                                                                      *
************************************************************************
*                                                                      *
*     Before we start we need to backtransform the coefficients of the
*     SO states from the basis of the SF states which diagonalize the
*     SF Hamiltonian to the basis of the original SF states. This since
*     all transition moments, whether or retrived from disk or
*     recomputed, are in the basis of the original SF states.
*
C Mapping from spin states to spin-free state:
      Call mma_allocate(MAPST,nSS,Label='MAPST')
      ISS=0
      DO ISTATE=1,MSTATE
         JOB=JBNUM(ISTATE)
         MPLET=MLTPLT(JOB)
         DO MSPROJ=-MPLET+1,MPLET-1,2
            ISS=ISS+1
            MAPST(ISS)=ISTATE
       END DO
      END DO
*
*     Let us transform the coefficients in USOR and USOI
*
      Do iSS = 1, nSS
         Do JSS = 1, nSS
            jSS_=MAPST(JSS)
            tmp_R=0.0D0
            tmp_I=0.0D0
            Do kSS = 1, nSS
               kSS_=MAPST(kSS)
               tmp_R=tmp_R + USOR(kSS,iSS)*EigVec(kSS_,jSS_)
               tmp_I=tmp_I + USOI(kSS,iSS)*EigVec(kSS_,jSS_)
            End Do
            VSOR(JSS,ISS)=tmp_R
            VSOI(JSS,ISS)=tmp_I
         End Do
      End Do
      Call mma_deallocate(MAPST)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End Subroutine USOTRANS
