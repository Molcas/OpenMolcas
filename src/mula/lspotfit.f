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
      Subroutine LSPotFit(r01,energy1,grad1,Hess1,D3_1,D4_1,            &
     &     r02,energy2,grad2,Hess2,D3_2,D4_2,                           &
     &     r00,energy0,r_min,FitCoef,mMat,stand_dev,max_err,            &
     &     use_weight,max_term,pot,nosc,numcoef)
!!
!!  Purpose:
!!    Perform a least squares fit of the potentiag at two different
!!    centra, r01 and r02.
!!
!!  Uses:
!!    Constants
!!    Linalg
!!    FCMod
!!
!!  Written by:
!!    Niclas Forsberg,
!!    Dept. of Theoretical Chemistry, Lund University, 1996.
!!
!       Use Linalg
!       Use FCMod
!       Use TabMod
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
      Parameter ( mxdeg = 6)
      Real*8 r01(nosc),r02(nosc),r00(nosc),r_min(nosc)
      Real*8 grad1(nosc),grad2(nosc)
      Real*8 Hess1(nosc,nosc),Hess2(nosc,nosc)
      Dimension D3_1(nosc,nosc,nosc),D3_2 (nosc,nosc,nosc)
      Dimension D4_1(nosc,nosc,nosc,nosc),D4_2 (nosc,nosc,nosc,nosc)
      Real*8 FitCoef (numcoef,1)
      Integer mMat (0:numcoef-1,nosc)
      Real*8   stand_dev,max_err

      Logical   use_weight,pot
      Integer  nTabDim
#include "WrkSpc.fh"
!!
!!---- Initialize.
      nvar = nosc
      mTabDim = numcoef-1
      Call TabDim_drv(max_term,nOsc,nTabDim)
      max_mOrd = nTabDim-1
      Call TabDim_drv(max_term-1,nOsc,nTabDim)
      max_mInc = nTabDim-1
      n_mDec=mTabDim+1
      Call GetMem('mDec','Allo','Inte',ipmDec,n_mDec*nvar)
      Call GetMem('mInc','Allo','Inte',ipmInc,n_mDec*nvar)
      l_mMat=nOsc
      Call MakeTab(max_term,max_mOrd,max_mInc,mMat,iWork(ipmInc),       &
     &  iWork(ipmDec),                                                  &
     &  l_mMat )
      nterm = max_mOrd+1
!!
!!---- Set up weight matrix and right hand side.
      nDim = 2*nterm
      !nDim = nVar*nterm !!!!
      If ( pot ) nDim = nDim+1
      Call GetMem('rhs','Allo','Real',iprhs,nDim)

      Call GetMem('weight','Allo','Real',ipweight,nDim*nDim)

      call dcopy_(ndim**2,[0.0d0],0,Work(ipweight),1)
      n = 1
      Work(ipweight+n+nDim*(n-1)-1) = 1.0d4
      Work(ipweight+n+nterm+nDim*(n+nterm-1)-1) = 1.0d4
      Work(iprhs+n-1) = energy1
      Work(iprhs+n+nterm-1) = energy2
      n = n+1
      If ( max_term.gt.0 ) Then
      Do i = 1,nvar
      Work(ipweight+n+nDim*(n-1)-1) = 1.0d3
      Work(ipweight+n+nterm+nDim*(n+nterm-1)-1) = 1.0d3
      Work(iprhs+n-1) = grad1(i)
      Work(iprhs+n+nterm-1) = grad2(i)
      n = n+1
      End Do
      If ( max_term.gt.1 ) Then
      Do i = 1,nvar
      Do j = i,nvar
      Work(ipweight+n+nDim*(n-1)-1) = 1.0d2
      Work(ipweight+n+nterm+nDim*(n+nterm-1)-1) = 1.0d2
      Work(iprhs+n-1) = Hess1(i,j)
      Work(iprhs+n+nterm-1) = Hess2(i,j)
      n = n+1
      End Do
      End Do
      If ( max_term.gt.2 ) Then
      Do i = 1,nvar
      Do j = i,nvar
      Do k = j,nvar
      Work(ipweight+n+nDim*(n-1)-1) = 1.0d1
      Work(ipweight+n+nterm+nDim*(n+nterm-1)-1) =                       &
     &                                                         1.0d1
      Work(iprhs+n-1) = D3_1(i,j,k)
      Work(iprhs+n+nterm-1) = D3_2(i,j,k)
      n = n+1
      End Do
      End Do
      End Do
      If ( max_term.gt.3 ) Then
      Do i = 1,nvar
      Do j = i,nvar
      Do k = j,nvar
      Do l = k,nvar
      Work(ipweight+n+nDim*(n-1)-1) = 1.0d0
      Work(ipweight+n+nterm+                                            &
     &                               nDim*(n+nterm-1)-1) = 1.0d0
      Work(iprhs+n-1) = D4_1(i,j,k,l)
      Work(iprhs+n+nterm-1) = D4_2(i,j,k,l)
      n = n+1
      End Do
      End Do
      End Do
      End Do
      End If
      End If
      End If
      End If
      If ( pot ) Then
      Work(iprhs+nDim-1) = energy0+20000.0d0/HarToRcm
      Work(ipweight+nDim+nDim*(nDim-1)-1) = 1.0d4
      End If
!!
!!----
      l_vpow=mxdeg+1
      Call GetMem('vpow','Allo','Real',ipvpow,l_vpow*nvar)
      Call GetMem('Tmat','Allo','Real',ipTmat,nDim*nterm)
      call dcopy_(nDim*nterm,[0.0d0],0,work(ipTmat),1)
!       Tmat = 0.0d0
      Call GetMem('x','Allo','Real',ipx,nvar)

      nrow = 1
      Do m = 1,2
      mrow = 1
      If ( m.eq.1 ) Then
      do iv=1,nvar
      Work(ipx+iv-1) = r01(iv)-r00(iv)
      enddo
      Else
      do iv=1,nvar
      Work(ipx+iv-1) = r02(iv)-r00(iv)
      enddo
!             x = r02-r00
      End If
!!
!!---- Calculate powers of individual variable values.
      Do ivar = 1,nvar
      pow = 1.0d0
      work(ipvpow+1+l_vpow*(ivar-1)-1) = 1.0d0
      Do i = 1,mxdeg
      pow = pow*Work(ipx+ivar-1)
      Work(ipvpow+i+1+l_vpow*(ivar-1)-1) = pow
      End Do
      End Do
!!
!!---- Calculate value of each polynomial term at this point.
      Do iterm = 1,nterm
      ip = mMat(iterm-1,1)
      t = Work(ipvpow+ip)
      Do ivar = 2,nvar
      ip = mMat(iterm-1,ivar)
      t = t*Work(ipvpow+ip+1+l_vpow*(ivar-1)-1)
      End Do
      Work(ipTmat+nrow+nDim*(iterm-1)-1) = t
      End Do
      mrow = mrow+1
      nrow = nrow+1
!!
!!--- First derivatives.
      If ( max_term.gt.0 ) Then
      Do ivar = 1,nvar
      Do iterm = 2,nterm
      irow =                                                            &
     &             iWork(ipmDec+mrow+n_mDec*(ivar-1)-1)+1+((m-1)*       &
     &             (nDim/2))
      jvar = nvar
      Do While (( mMat(iterm-1,jvar).eq.0 ).and.                        &
     &             ( jvar.gt.1 ))
      jvar = jvar-1
      End Do
      jterm = iWork(ipmDec+iterm+n_mDec*(jvar-1)-1)+1
      Work(ipTmat+nrow+nDim*(iterm-1)-1) =                              &
     &             work(ipx+jvar-1)*                                    &
     &             Work(ipTmat+nrow+nDim*(jterm-1)-1)
      If ( ivar.eq.jvar ) Then
      Work(ipTmat+nrow+nDim*(iterm-1)-1) =                              &
     &                Work(ipTmat+nrow+nDim*(iterm-1)-1)+               &
     &                Work(ipTmat+irow+nDim*(jterm-1)-1)
      End If
      End Do
      mrow = mrow+1
      nrow = nrow+1
      End Do
      End If
!!
!!--- Second derivatives.
      If ( max_term.gt.1 ) Then
      Do ivar = 1,nvar
      Do jvar = ivar,nvar
      Do iterm = 2,nterm
      irow = iWork(ipmDec+mrow+n_mDec*(ivar-1)-1)+1+                    &
     &                       ((m-1)*(nDim/2))
      jrow = iWork(ipmDec+mrow+n_mDec*(jvar-1)-1)+1+                    &
     &                       ((m-1)*(nDim/2))
      kvar = nvar
      Do While (( mMat(iterm-1,kvar).eq.0 ).and.                        &
     &                ( kvar.gt.1 ))
      kvar = kvar-1
      End Do
      jterm = iWork(ipmDec+iterm+n_mDec*(kvar-1)-1)+1
      Work(ipTmat+nrow+nDim*(iterm-1)-1) =                              &
     &                Work(ipx+kvar-1)*                                 &
     &                Work(ipTmat+nrow+nDim*(jterm-1)-1)
      If ( ivar.eq.kvar ) Then
      Work(ipTmat+nrow+nDim*(iterm-1)-1) =                              &
     &                   Work(ipTmat+nrow+nDim*(iterm-1)-1)+            &
     &                   Work(ipTmat+irow+nDim*(jterm-1)-1)
      End If
      If ( jvar.eq.kvar ) Then
      Work(ipTmat+nrow+nDim*(iterm-1)-1) =                              &
     &                   Work(ipTmat+nrow+nDim*(iterm-1)-1)+            &
     &                   Work(ipTmat+jrow+nDim*(jterm-1)-1)
      End If
      End Do
      mrow = mrow+1
      nrow = nrow+1
      End Do
      End Do
      End If
!!
!!--- Third derivatives.
      If ( max_term.gt.2 ) Then
      Do ivar = 1,nvar
      Do jvar = ivar,nvar
      Do kvar = jvar,nvar
      Do iterm = 2,nterm
      irow = iWork(ipmDec+mrow+                                         &
     &                   n_mDec*(ivar-1)-1)+1+                          &
     &                          ((m-1)*(nDim/2))
      jrow = iWork(ipmDec+mrow+n_mDec*(jvar-1)-1)+1+                    &
     &                          ((m-1)*(nDim/2))
      krow = iWork(ipmDec+mrow+n_mDec*(kvar-1)-1)+1+                    &
     &                          ((m-1)*(nDim/2))
      lvar = nvar
      Do While (( mMat(iterm-1,lvar).eq.0 ).and.                        &
     &                   ( lvar.gt.1 ))
      lvar = lvar-1
      End Do
      jterm = iWork(ipmDec+iterm+n_mDec*(lvar-1)-1)+1
      Work(ipTmat+nrow+nDim*(iterm-1)-1) =                              &
     &                   Work(ipx+lvar-1)*                              &
     &                   Work(ipTmat+nrow+nDim*(jterm-1)-1)
      If ( ivar.eq.lvar ) Then
      Work(ipTmat+nrow+nDim*(iterm-1)-1) =                              &
     &                      Work(ipTmat+nrow+nDim*(iterm-1)-1)+         &
     &                      Work(ipTmat+irow+nDim*(jterm-1)-1)
      End If
      If ( jvar.eq.lvar ) Then
      Work(ipTmat+nrow+nDim*(iterm-1)-1) =                              &
     &                      Work(ipTmat+nrow+nDim*(iterm-1)-1)+         &
     &                      Work(ipTmat+jrow+nDim*(jterm-1)-1)
      End If
      If ( kvar.eq.lvar ) Then
      Work(ipTmat+nrow+nDim*(iterm-1)-1) =                              &
     &                      Work(ipTmat+nrow+nDim*(iterm-1)-1)+         &
     &                      Work(ipTmat+krow+nDim*(jterm-1)-1)
      End If
      End Do
      mrow = mrow+1
      nrow = nrow+1
      End Do
      End Do
      End Do
      End If
!!
!!--- Fourth derivatives.
      If ( max_term.gt.3 ) Then
      Do ivar = 1,nvar
      Do jvar = ivar,nvar
      Do kvar = jvar,nvar
      Do lvar = kvar,nvar
      Do iterm = 2,nterm
      irow =                                                            &
     &  iWork(ipmDec+mrow+n_mDec*(ivar-1)-1)+1+((m-1)*(nDim/2))
      jrow =                                                            &
     &  iWork(ipmDec+mrow+n_mDec*(jvar-1)-1)+1+((m-1)*(nDim/2))
      krow =                                                            &
     &  iWork(ipmDec+mrow+n_mDec*(kvar-1)-1)+1+((m-1)*(nDim/2))
      lrow =                                                            &
     &  iWork(ipmDec+mrow+n_mDec*(lvar-1)-1)+1+((m-1)*(nDim/2))
      mvar = nvar
      Do While                                                          &
     &                      (( mMat(iterm-1,mvar).eq.0 ).and.           &
     &                      ( mvar.gt.1 ))
      mvar = mvar-1
      End Do
      jterm =                                                           &
     &                      iWork(ipmDec+iterm+n_mDec*(mvar-1)-1)+1
      Work(ipTmat+nrow+nDim*(iterm-1)-1) =                              &
     &                      Work(ipx+mvar-1)*                           &
     &                   Work(ipTmat+nrow+nDim*(jterm-1)-1)
      If ( ivar.eq.mvar ) Then
      Work(ipTmat+nrow+nDim*(iterm-1)-1) =                              &
     &                         Work(ipTmat+nrow+nDim*(iterm-1)-1)+      &
     &                         Work(ipTmat+irow+ndim*(jterm-1)-1)
      End If
      If ( jvar.eq.mvar ) Then
      work(ipTmat+nrow+nDim*(iterm-1)-1) =                              &
     &                         Work(ipTmat+nrow+nDim*(iterm-1)-1)+      &
     &                         Work(ipTmat+jrow+nDim*(jterm-1)-1)
      End If
      If ( kvar.eq.mvar ) Then
      Work(ipTmat+nrow+nDim*(iterm-1)-1) =                              &
     &                         Work(ipTmat+nrow+nDim*(iterm-1)-1)+      &
     &                         Work(ipTmat+krow+nDim*(jterm-1)-1)
      End If
      If ( lvar.eq.mvar ) Then
      Work(ipTmat+nrow+nDim*(iterm-1)-1) =                              &
     &                         Work(ipTmat+nrow+nDim*(iterm-1)-1)+      &
     &                         Work(ipTmat+lrow+nDim*(jterm-1)-1)
      End If
      End Do
      mrow = mrow+1
      nrow = nrow+1
      End Do
      End Do
      End Do
      End Do
      End If
      End Do
!!
      If ( pot ) Then
      Work(ipx) = r_min(1)-r00(1)
      Work(ipx+1) = r_min(2)-r00(2)
      Work(ipx+2) = 2.0d0*rpi-r_min(3)-r00(3)
!!
!!---- Calculate powers of individual variable values.
      Do ivar = 1,nvar
      pow = 1.0d0
      Work(ipvpow+1+l_vpow*(ivar-1)-1) = 1.0d0
      Do i = 1,mxdeg
      pow = pow*Work(ipx+ivar-1)
      Work(ipvpow+i+1+l_vpow*(ivar-1)-1) = pow
      End Do
      End Do
!!
!!---- Calculate value of each polynomial term at this point.
      Do iterm = 1,nterm
      ip = mMat(iterm-1,1)
      t = Work(ipvpow+ip)
      Do ivar = 2,nvar
      ip = mMat(iterm-1,ivar)
      t = t*Work(ipvpow+ip+1+l_vpow*(ivar-1)-1)
      End Do
      Work(ipTmat+nrow+nDim*(iterm-1)-1) = t
      End Do
      End If
!!
!!---- Calculate equation matrix, T(t)*weight*T, and T(t)*weight*rhs.
      Call GetMem('Temp','Allo','Real',ipTemp,nterm*nDim)
      Call GetMem('Equmat','Allo','Real',ipEqumat,nterm*nterm)

      Call DGEMM_('T','N',                                              &
     &            nterm,nDim,nDim,                                      &
     &            1.0d0,Work(ipTmat),nDim,                              &
     &            Work(ipweight),nDim,                                  &
     &            0.0d0,Work(ipTemp),nterm)
      Call DGEMM_('N','N',                                              &
     &            nterm,nterm,nDim,                                     &
     &            1.0d0,Work(ipTemp),nterm,                             &
     &            Work(ipTmat),nDim,                                    &
     &            0.0d0,Work(ipEqumat),nterm)
!!
      n = 2
      If ( max_term.gt.0 ) Then
      Do i = 1,nvar
      n = n+1
      End Do
      If ( max_term.gt.1 ) Then
      Do i = 1,nvar
      Do j = i,nvar
      n = n+1
      End Do
      End Do
      If ( max_term.gt.2 ) Then
      Do i = 1,nvar
      Do j = i,nvar
      Do k = j,nvar
      Work(ipEqumat+n+nterm*(n-1)-1) =                                  &
     &                   Work(ipEqumat+n+nterm*(n-1)-1)+1.0d-1
      n = n+1
      End Do
      End Do
      End Do
      If ( max_term.gt.3 ) Then
      Do i = 1,nvar
      Do j = i,nvar
      Do k = j,nvar
      Do l = k,nvar
      Work(ipEqumat+n+nterm*(n-1)-1) =                                  &
     &                         Work(ipEqumat+n+nterm*(n-1)-1)+1.0d-1
      n = n+1
      End Do
      End Do
      End Do
      End Do
      End If
      End If
      End If
      End If
!!
      Call DGEMM_('N','N',                                              &
     &            nterm,1,nDim,                                         &
     &            1.0d0,Work(ipTemp),nterm,                             &
     &            Work(iprhs),nDim,                                     &
     &            0.0d0,FitCoef,nterm)
!!
!!---- Solve the resulting equation system.
      Call Dool_MULA(Work(ipEqumat),nterm,nterm,FitCoef,nterm,nterm,det)
!!
      Call GetMem('weight','Free','Real',ipweight,nDim*nDim)
      Call GetMem('mDec','Free','Inte',ipmDec,n_mDec*nvar)
      Call GetMem('mInc','Free','Inte',ipmInc,n_mDec*nvar)
      Call GetMem('vpow','Allo','Real',ipvpow,l_vpow*nvar)
      Call GetMem('x','Free','Real',ipx,nvar)
      Call GetMem('Tmat','Free','Real',ipTmat,nDim*nterm)
      Call GetMem('rhs','Free','Real',iprhs,nDim)
      Call GetMem('Equmat','Free','Real',ipEqumat,nterm*nterm)
      Call GetMem('Temp','Free','Real',ipTemp,nterm*nDim)
!!
! Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(stand_dev)
         Call Unused_real(max_err)
         Call Unused_logical(use_weight)
      End If
      End
