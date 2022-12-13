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
* Copyright (C) 2022, Jie J. Bao                                       *
************************************************************************
******************************************************************
* history:                                                       *
* Jie J. Bao, on Apr. 11, 2022, created this file.               *
******************************************************************


      Subroutine CMSNewton(R,GDorbit,GDstate,Dgorbit,Dgstate,nGD)
      use CMS, only:CMSNotConverged,CMSThres,NeedMoreStep,
     &              nPosHess,LargestQaaGrad,NCMSScale
      use stdalloc, only : mma_allocate, mma_deallocate
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.h"
      INTEGER nGD
      Real*8 R(lRoots**2),
     &       GDorbit(nGD),GDstate(nGD),
     &       Dgorbit(nGD),Dgstate(nGD)
      Real*8,DIMENSION(:),Allocatable::X,Hess,Grad,EigVal,deltaR,DDg,
     &                                 XScr,GScr,ScrDiag,
     &                                 RCopy,GDCopy,DgCopy

      Real*8,DIMENSION(:,:),Allocatable::RotMat
      INTEGER iStep,nDDg,lRoots2,NAC2,
     &        nSPair,nSPair2,nScr
      Real*8 Qnew,Qold
      Logical Saved

*     preparation
      lRoots2=lRoots**2
      NAC2=NAC**2
      nDDg=lRoots2**2
      nSPair=(lRoots-1)*lRoots/2
      nSPair2=nSPair**2
      CMSThres=CMSThreshold
      CALL mma_allocate(DDg    ,nDDg     )
      Call mma_allocate(X      ,nSPair   )
      Call mma_allocate(XScr   ,nSPair   )
      Call mma_allocate(GScr   ,nSPair   )
      Call mma_allocate(Grad   ,nSPair   )
      Call mma_allocate(Hess   ,nSPair2  )
      Call mma_allocate(EigVal ,nSPair   )
      CALL mma_allocate(DeltaR ,lRoots2  )
      CALL mma_allocate(GDCopy ,nGD      )
      CALL mma_allocate(DgCopy ,nGD      )
      CALL mma_allocate(RCopy  ,lRoots2  )
      CALL mma_allocate(RotMat ,lRoots,lRoots)
*     Step 0
      iStep=0
      Qold=0.0d0
*     Note that the following six lines appear as a group
      CALL RotGD(GDstate,R,nGD,lRoots,NAC2)
      CALL RotGD(Dgstate,R,nGD,lRoots,NAC2)
      CALL TransposeMat(Dgorbit,Dgstate,nGD,lRoots2,NAC2)
      CALL TransposeMat(GDorbit,GDstate,nGD,lRoots2,NAC2)
      CALL CalcDDg(DDg,GDorbit,Dgorbit,nDDg,nGD,lRoots2,NAC2)
      CALL CalcQaa(Qnew,DDg,lRoots,nDDg)
      nPosHess=0
      LargestQaaGrad=0.0d0
      Qold=Qnew
      CALL PrintCMSIter(iStep,Qnew,Qold,R,lRoots)
      CALL CalcGradCMS(Grad,DDg,nDDg,lRoots,nSPair)
      CALL CalcHessCMS(Hess,DDg,nDDg,lRoots,nSPair)
      CALL GetDiagScr(nScr,Hess,EigVal,nSPair)
      CALL mma_allocate(ScrDiag,nScr)

*     Starting iteration
      DO WHILE(CMSNotConverged)
       iStep=iStep+1
       Qold=Qnew
       IF(iStep.gt.iCMSIterMax) THEN
         write(6,'(4X,A)')'NOT CONVERGED AFTER MAX NUMBER OF CYCLES'
         Exit
       END IF

       CALL DCopy_(lRoots2,R,1,RCopy,1)
       CALL DCopy_(nGD,GDState,1,GDCopy,1)
       CALL DCopy_(nGD,DgState,1,DgCopy,1)

       CALL CalcNewX(X,Hess,Grad,nSPair,
     &               XScr,GScr,EigVal,ScrDiag,nScr)
       CALL UpDateRotMat(R,DeltaR,X,lRoots,nSPair)

       CALL RotGD(GDstate,DeltaR,nGD,lRoots,NAC2)
       CALL RotGD(Dgstate,DeltaR,nGD,lRoots,NAC2)
       CALL TransposeMat(Dgorbit,Dgstate,nGD,lRoots2,NAC2)
       CALL TransposeMat(GDorbit,GDstate,nGD,lRoots2,NAC2)
       CALL CalcDDg(DDg,GDorbit,Dgorbit,nDDg,nGD,lRoots2,NAC2)
       CALL CalcQaa(Qnew,DDg,lRoots,nDDg)

       NCMSScale=0
       Saved=.true.
       IF((Qold-Qnew).gt.CMSThreshold) THEN
        If(iStep.gt.ICMSIterMin) Then
*        When Onew is less than Qold, scale the rotation matrix
         CALL CMSScaleX(X,R,DeltaR,Qnew,Qold,
     &                  RCopy,GDCopy,DgCopy,
     &                  GDstate,GDOrbit,Dgstate,DgOrbit,DDg,
     &                  nSPair,lRoots2,nGD,NAC2,nDDg,Saved)
        End If
       END IF
       CALL PrintCMSIter(iStep,Qnew,Qold,R,lRoots)
       CALL AntiOneDFoil(RotMat,R,lRoots,lRoots)
       CALL PrintMat('ROT_VEC','CMS-PDFT temp',
     &                RotMat,lroots,lroots,7,13,'T')

       IF(.not. Saved) THEN
        CMSNotConverged=.true.
*        Exit
       END IF
*      sanity check
       IF(abs(Qnew-Qold).lt.CMSThreshold) THEN
        CMSNotConverged=.false.
        If(NeedMoreStep)         CMSNotConverged=.true.
        If(iStep.lt.iCMSIterMin) CMSNotConverged=.true.
        If(NCMSScale.gt.0)       CMSNotConverged=.true.
       END IF
       IF(CMSNotConverged) THEN
        CALL CalcGradCMS(Grad,DDg,nDDg,lRoots,nSPair)
        CALL CalcHessCMS(Hess,DDg,nDDg,lRoots,nSPair)
       ELSE
        write(6,'(4X,A)')'CONVERGENCE REACHED'
       END IF
      END DO

      CALL mma_deallocate(DDg    )
      Call mma_deallocate(X      )
      Call mma_deallocate(XScr   )
      Call mma_deallocate(GScr   )
      Call mma_deallocate(Grad   )
      Call mma_deallocate(Hess   )
      Call mma_deallocate(EigVal )
      CALL mma_deallocate(DeltaR )
      CALL mma_deallocate(ScrDiag)
      CALL mma_deallocate(GDCopy )
      CALL mma_deallocate(DgCopy )
      CALL mma_deallocate(RCopy  )
      CALL mma_deallocate(RotMat )
      RETURN
      End Subroutine



