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
! Copyright (C) 1996, Niclas Forsberg                                  *
!***********************************************************************
!!-----------------------------------------------------------------------!
!!
      Subroutine PotEnergy(A,nMat,iCre,iAnn,energy,grad,Hess,           &
     &       D3,D4,max_term,W,max_ord,nosc,nOscOld)
!!
!       Use TabMod
!!
!!  Purpose:
!!    Calculate matrix elements of potential energy terms.
!!
!!  Written by:
!!    Niclas Forsberg,
!!    Dept. of Theoretical Chemistry, Lund University, 1996.
!!
      Implicit Real*8 ( a-h,o-z )
#include "dims.fh"
      Real*8 A  (0:mdim1,0:ndim1)
      Integer nMat  (0:ndim1,ndim2)
      Integer iAnn  (0:ndim1,ndim2)
      Integer iCre  (0:ndim1,ndim2)
      Real*8 rdx (4)
      Real*8 grad  (noscold)
      Real*8 Hess  (noscold,noscold)
      Real*8 D3  (noscold,noscold,noscold)
      Real*8 D4   (noscold,noscold,noscold,noscold)
      Real*8 W    (noscold,nosc)
#include "WrkSpc.fh"
!!
!!---- Zeroth order term.
      call dcopy_(max_ord+1,[Energy],0,A,max_Ord+2)
      rdx(1) = 1.0d0
      rdx(2) = 1.0d0
      rdx(3) = 1.0d0
      rdx(4) = 1.0d0
      Call GetMem('Temp','Allo','Real',ipTemp,nOscOld**4)
!!
!!---- First order terms.
      If ( max_term.gt.0 ) Then
      Call GetMem('grad_2','Allo','Real',ipgrad_2,nOsc)
      Call DGEMM_('T','N',                                              &
     &            1,nOsc,nOscOld,                                       &
     &            1.0d0,grad,nOscOld,                                   &
     &            W,nOscOld,                                            &
     &            0.0d0,Work(ipgrad_2),1)
      Call Mul1(nMat,A,icre,iann,Work(ipgrad_2),max_ord,nosc,rdx)
      Call GetMem('grad_2','Free','Real',ipgrad_2,nOsc)
      End If
!!
!!---- Second order terms.
      If ( max_term.gt.1 ) Then
      Call GetMem('Hess_2','Allo','Real',ipHess_2,nOsc**2)
      Call DGEMM_('T','N',                                              &
     &            nOscOld,nOsc,nOscOld,                                 &
     &            1.0d0,Hess,nOscOld,                                   &
     &            W,nOscOld,                                            &
     &            0.0d0,Work(ipTemp),nOscOld)
      Call DGEMM_('T','N',                                              &
     &            nOsc,nOsc,nOscOld,                                    &
     &            1.0d0,Work(ipTemp),nOscOld,                           &
     &            W,nOscOld,                                            &
     &            0.0d0,Work(ipHess_2),nOsc)
      Call Mul2(nMat,A,icre,iann,Work(ipHess_2),max_ord,nosc,rdx)
      Call GetMem('Hess_2','Free','Real',ipHess_2,nOsc**2)
      End If
!!
!!---- Third order terms.
      If ( max_term.gt.2 ) Then
      Call GetMem('D3_2','Allo','Real',ipD3_2,nOsc**3)
      Call DGEMM_('T','N',                                              &
     &            nOscOld**2,nOsc,nOscOld,                              &
     &            1.0d0,D3,nOscOld,                                     &
     &            W,nOscOld,                                            &
     &            0.0d0,Work(ipTemp),nOscOld**2)
      Call DGEMM_('T','N',                                              &
     &            nOsc*nOscOld,nOsc,nOscOld,                            &
     &            1.0d0,Work(ipTemp),nOscOld,                           &
     &            W,nOscOld,                                            &
     &            0.0d0,Work(ipD3_2),nOsc*nOscOld)
      Call DGEMM_('T','N',                                              &
     &            nOsc**2,nOsc,nOscOld,                                 &
     &            1.0d0,Work(ipD3_2),nOscOld,                           &
     &            W,nOscOld,                                            &
     &            0.0d0,Work(ipTemp),nOsc**2)
      call dcopy_(nOsc**3,Work(ipTemp),1,Work(ipD3_2),1)
      Call Mul3(nMat,A,icre,iann,Work(ipD3_2),max_ord,nosc,rdx)
      Call GetMem('D3_2','Free','Real',ipD3_2,nOsc**3)
      End If
!!
!!---- Fourth order terms.
      If ( max_term.gt.3 ) Then
      Call GetMem('D4_2','Allo','Real',ipD4_2,nOsc**4)

      Call DGEMM_('T','N',                                              &
     &            nOscOld**3,nOsc,nOscOld,                              &
     &            1.0d0,D4,nOscOld,                                     &
     &            W,nOscOld,                                            &
     &            0.0d0,Work(ipTemp),nOscOld**3)
      Call DGEMM_('T','N',                                              &
     &            nOsc*nOscOld**2,nOsc,nOscOld,                         &
     &            1.0d0,Work(ipTemp),nOscOld,                           &
     &            W,nOscOld,                                            &
     &            0.0d0,Work(ipD4_2),nOsc*nOscOld**2)
      Call DGEMM_('T','N',                                              &
     &            nOsc**2*nOscOld,nOsc,nOscOld,                         &
     &            1.0d0,Work(ipD4_2),nOscOld,                           &
     &            W,nOscOld,                                            &
     &            0.0d0,Work(ipTemp),nOsc**2*nOscOld)
      Call DGEMM_('T','N',                                              &
     &            nOsc**3,nOsc,nOscOld,                                 &
     &            1.0d0,Work(ipTemp),nOscOld,                           &
     &            W,nOscOld,                                            &
     &            0.0d0,Work(ipD4_2),nOsc**3)
      Call Mul4(nMat,A,icre,iann,Work(ipD4_2),max_ord,nosc,rdx)
      Call GetMem('D4_2','Free','Real',ipD4_2,nOsc**4)
      End If
!!
      Call GetMem('Temp','Free','Real',ipTemp,nOscOld**4)
!!
      End
