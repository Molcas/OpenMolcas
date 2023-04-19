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
      Subroutine CMSOpt(TUVX)
* ****************************************************************
* history:                                                       *
* Jie J. Bao, on Apr. 07, 2022, created this file.               *
* ****************************************************************
      use stdalloc, only : mma_allocate, mma_deallocate
      use CMS, only: CMSNotConverged,RGD
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.h"

      Real*8,DIMENSION(NACPR2)::TUVX
      Real*8,DIMENSION(:),Allocatable::Gtuvx,R,
     &                                 GDstate,GDorbit,
     &                                 Dgstate,Dgorbit
      Real*8,DIMENSION(:,:),Allocatable::RotMat

      INTEGER nTUVX,nGD,lRoots2,NAC2

      CHARACTER(len=16)::VecName
******************************************************************
*     some notes on the arrays:                                  *
*     Gtuvx : two-electron integral, g_tuvx                      *
*     GD    : "generalized 1-e density matrix"                   *
*              GD^KL: transition density matrix from L to K      *
*              GD^KK: density matrix for state K                 *
*     Dg    :  sum_{vx}{GD^KL_vx * g_tuvx}                       *
*     In GDorbit and Dgorbit, the leading index is orbital index;*
*     In GDstate and Dgstate, the leading index is state index.  *
*                                                                *
*     DDg   :  sum_{tuvx}{GD^KL_tu * GD^MN_vx * g_tuvx}          *
*      namely, sum_{tu}{GD^KL_tu * Dg^MN_tu}                     *
******************************************************************

      NAC2=NAC**2
      nTUVX=NAC2**2
      lRoots2=lRoots**2
      nGD=lRoots2*NAC2

      CMSNotConverged=.true.

******Memory Allocation
      CALL mma_allocate(R      ,lRoots2)
      CALL mma_allocate(GDstate,nGD    )
      CALL mma_allocate(Dgstate,nGD    )
      CALL mma_allocate(GDorbit,nGD    )
      CALL mma_allocate(Dgorbit,nGD    )
      CALL mma_allocate(Gtuvx  ,nTUVX  )
      CALL mma_allocate(RGD    ,lRoots2)
      CALL mma_allocate(RotMat ,lRoots,lRoots)

******Calculate generalized density mtrix
      CALL UnzipTUVX(TUVX,Gtuvx,nTUVX)

C      write(6,*) 'Gtuvx matrix'
C      CALL RecPrt(' ',' ',Gtuvx,NAC2,NAC2)

      CALL CalcGD(GDorbit,nGD)
C      write(6,*) 'GD matrix orbital-leading'
C      CALL RecPrt(' ',' ',GDorbit,NAC2,lRoots2)
      CALL CalcDg(Dgorbit,GDorbit,Gtuvx,nGD,nTUVX,NAC,lRoots)
C      write(6,*) 'Dg matrix orbital-leading'
C      CALL RecPrt(' ',' ',Dgorbit,NAC2,lRoots2)

      CALL mma_deallocate(Gtuvx  )

      CALL TransposeMat(Dgstate,Dgorbit,nGD,NAC2,lRoots2)
      CALL TransposeMat(GDstate,GDorbit,nGD,NAC2,lRoots2)

******Load initial rotation matrix
      CALL InitRotMat(RotMat,lRoots,
     &                trim(CMSStartMat),len_trim(CMSStartMat))

      CALL OneDFoil(R,RotMat,lRoots,lRoots)

******Print header of CMS iterations
      CALL CMSHeader(trim(CMSStartMat),len_trim(CMSStartMat))

******Start CMS Optimization
      CMSNotConverged=.true.
      CALL CMSNewton(R,GDorbit,GDstate,Dgorbit,Dgstate,nGD)

******Print end of CMS intermediate-state optimization
      CALL CMSTail()

******Save rotation matrix
      CALL AntiOneDFoil(RotMat,R,lRoots,lRoots)
      VecName='CMS-PDFT'
      CALL PrintMat('ROT_VEC',VecName,RotMat,lroots,lroots,7,16,'T')

******releasing memory
      CALL mma_deallocate(R      )
      CALL mma_deallocate(GDstate)
      CALL mma_deallocate(Dgstate)
      CALL mma_deallocate(GDorbit)
      CALL mma_deallocate(Dgorbit)
      CALL mma_deallocate(RGD    )
      CALL mma_deallocate(RotMat )

******check convergence
      IF(CMSNotConverged) THEN
       Call WarningMessage(2,'CMS Intermediate States Not Converged')
       Call Quit(_RC_NOT_CONVERGED_)
      END IF

      RETURN
      End Subroutine

