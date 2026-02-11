!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
!*****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 11, 2022, created this file.               *
!*****************************************************************

      Subroutine CalcGD(GD,nGD)
      use stdalloc, only: mma_allocate, mma_deallocate
      use lucia_data, only: DStmp, Dtmp
      use Lucia_Interface, only: Lucia_Util
      use rasscf_global, only: lRoots, NAC, iADR15
      use general_data, only: JOBIPH,NCONF
      Implicit None

#include "warnings.h"
      INTEGER nGD
      Real*8 GD(nGD)

      INTEGER CIDisk1,CIDisk2
      INTEGER p,q,ipq,iqp,NAC2,IOffNIJ1,IOffNIJ2, jRoot, kRoot
      Real*8, Allocatable:: SDtmp(:), TmpD(:)
      Real*8, Allocatable, Target:: VecL(:), VecR(:)

      NAC2=NAC**2
      Call mma_allocate(VecL,NConf,Label='VecL')
      Call mma_allocate(VecR,NConf,Label='VecR')
      Call mma_allocate(TmpD,NAC**2,Label='TmpD')
      Call mma_allocate(SDtmp,NAC**2,Label='SDtmp')
      SDtmp(:)=DStmp(:)
      TmpD(:)=DTmp(:)
      CIDisk1=IADR15(4)
      Do jRoot=1,lRoots
       Call DDafile(JOBIPH,2,VecL,nConf,CIDisk1)
       CIDisk2=IADR15(4)
       Do kRoot=1,jRoot-1
        Call DDafile(JOBIPH,2,VecR,nConf,CIDisk2)
        Call Lucia_Util('Densi',                                        &
     &                   CI_Vector=VecL(:),                             &
     &                   RVec=VecR(:))
        IOffNIJ1=(lRoots*(jRoot-1)+kRoot-1)*NAC2
        IOffNIJ2=(lRoots*(kRoot-1)+jRoot-1)*NAC2
!        write(6,*)'GD matrix',jRoot,kRoot
!        CALL RecPrt(' ',' ',Dtmp,NAC,NAC)
        Call DCopy_(NAC2,Dtmp,1,GD(IOffNIJ1+1),1)
         dO q=1,NAC
          do p=1,NAC
          ipq=(q-1)*NAC+p
          iqp=(p-1)*NAC+q
          GD(IOffNIJ2+iqp)=Dtmp(ipq)
!          GDMat(NIJ2,q,p)=Dtmp(q+(p-1)*NAC)
          end do
         eND dO
       End Do
       kRoot=jRoot
       Call DDafile(JOBIPH,2,VecR,nConf,CIDisk2)
       Call Lucia_Util('Densi',                                         &
     &                 CI_Vector=VecL(:),                               &
     &                 RVec=VecR(:))
       IOffNIJ1=(lRoots+1)*(jRoot-1)*NAC2
!       write(6,*)'GD matrix',jRoot,kRoot
!       CALL RecPrt(' ',' ',Dtmp,NAC,NAC)
       Call DCopy_(NAC2,Dtmp,1,GD(IOffNIJ1+1),1)
      End DO
      DStmp(:)=SDtmp(:)
      Dtmp(:)=TmpD(:)
      Call mma_deallocate(SDtmp)
      Call mma_deallocate(TmpD)
      Call mma_deallocate(VecL)
      Call mma_deallocate(VecR)
      END Subroutine


