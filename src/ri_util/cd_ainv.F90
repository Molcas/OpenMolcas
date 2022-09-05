!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine CD_AInv(A,n,AInV,Thr_CD)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
      Real*8 A(n,n), AInv(n,n)
      Real*8, Allocatable :: ADiag(:), QVec(:,:)
      Integer, Allocatable :: iADiag(:)
#ifdef _ACCURACY_
      Real*8, Allocatable :: Tmp(:,:), Tmp2(:,:)
#endif
!
      Call mma_allocate(ADiag,n,Label='ADiag')
      Call mma_allocate(iADiag,n,Label='iADiag')
!
      iSeed=77
      Lu_A=IsFreeUnit(iSeed)
      Call DaName_MF_WA(Lu_A,'AMat09')
!
      iDisk=0
      Call dDaFile(Lu_A,1,A,n**2,iDisk)
!
!     Call RecPrt('A',' ',A,n,n)
!
      iSeed=iSeed+1
      Lu_Q=IsFreeUnit(iSeed)
      Call DaName_MF_WA(Lu_Q,'QMat09')
!
      call dcopy_(n,A,n+1,ADiag,1)
!
      Call CD_AInv_(n,m,ADiag,iADiag,Lu_A,Lu_Q,Thr_CD)
!
      Call mma_deallocate(ADiag)
      Call mma_deallocate(iADiag)
!
      Call mma_allocate(QVec,n,m,Label='QVec')
!
      iDisk=0
      Call dDaFile(Lu_Q,2,QVec,n*m,iDisk)
!
!     Call RecPrt('QVec','(6G20.10)',QVec,n,m)
      Call DGEMM_('N','T',n,n,m,                                        &
     &            One,QVec,n,                                           &
     &                  QVec,n,                                         &
     &            Zerp,AInv,n)
!     Call RecPrt('AInv',' ',AInv,n,n)
      Call DaEras(Lu_Q)
      Call mma_deallocate(QVec)
!                                                                      *
!***********************************************************************
!                                                                      *
!     Check the accuracy I-AA^1
!
#ifdef _ACCURACY_
      Call mma_allocate(Tmp,n,n,Label='Tmp')
!---
      Tmp(:,:)=Zero
!     I
      call dcopy_(n,One,0,Tmp,n+1)
!     I-AA^-1
      Call DGEMM_('N','N',n,n,n,                                        &
     &           -One,A,n,                                              &
     &                  AInv,n,                                         &
     &            One,Tmp,n)
      Call RecPrt('I-AA^-1','(6G20.12)',Tmp,n,n)
!
      Call DGEMM_('N','N',n,n,n,                                        &
     &            One,A,n,                                              &
     &                  AInv,n,                                         &
     &            Zero,Tmp,n)

      Call mma_allocate(Tmp2,n,n,Label='Tmp2')
      Tmp2(:,:)=Zero
      call dcopy_(n,One,0,Tmp2,n+1)
      Call DGEMM_('N','N',n,n,n,                                        &
     &           -One,Tmp,n,                                            &
     &                Tmp,n,                                            &
     &            One,Tmp2,n)
      Call RecPrt('I-AA^-1AA^-1','(6G20.12)',Tmp2,n,n)
!---
      Call mma_deallocate(Tmp2)
      Call mma_deallocate(Tmp)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
      Return
      End
