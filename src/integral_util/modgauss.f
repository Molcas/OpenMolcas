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
      Subroutine ModGauss(Z,A,Xi,w)
      use Constants, only: rBohr
      Implicit None
      Real*8 Z, Xi, w
      Integer A

      Real*8 Facts(2,0:12), Errors(0:12), g(2), H(2,2), HInv(2,2),
     &       Step(2)
      Data Facts/ 0.0D0, 0.0D0,
     &            1.0D0, 0.0D0,
     &           -1.0D0, 0.0D0,
     &            2.0D0, 0.0D0,
     &           -2.0D0, 0.0D0,
     &            0.0D0, 1.0D0,
     &            0.0D0,-1.0D0,
     &            0.0D0, 2.0D0,
     &            0.0D0,-2.0D0,
     &            1.0D0, 1.0D0,
     &           -1.0D0, 1.0D0,
     &            1.0D0,-1.0D0,
     &           -1.0D0,-1.0D0/
      Real*8 A3, RMS, T, r_90, Thr, X, Delta_W, Delta_R, W0, R0, Det,
     &       e, f, r, x1, x2
      Integer MaxIter, Iter, i, iNeg
!                                                                      *
!***********************************************************************
!                                                                      *
      f(x,w)=(1.0D0+w*x**2)*exp(-x**2)
      R(w)=Sqrt(2.0D0*RMS**2*(3.0D0*w+2.0D0)/(3.0D0*(2.0D0+5.0D0*w)))
      E(x1,x2,w) = (f(x1,w)-0.9D0)**2 + (f(x2,w)-0.1D0)**2
!                                                                      *
!***********************************************************************
!                                                                      *
!     Write (*,*) 'A=',A
      A3    = 1.0d0*DBLE(A)**(1.0d0/3.0d0)
!     DA Eq. 51
      RMS=0.836d0*A3+0.570d0                    ! fm
      RMS=RMS*1.0D-15/rBohr  ! bohr
!     Write (6,*) 'RMS:',RMS
      w=0.0D0
      Xi=1.0D0/R(w)**2
      If (A.le.9) Then
!        Write (*,*) 'Use the Gaussian model!'
         Return
      End If
!     Write (*,*) 'Xi:',Xi
!     Xi=1.5D0/RMS**2
!     Write (*,*) 'Xi:',Xi
!
      T=2.30D0                              ! fm
      T=T*1.0D-15/rBohr  ! bohr
!
!     Start seeds
!
      w=0.5D0
      r_90=RMS/2.0D0
!                                                                      *
!***********************************************************************
!                                                                      *
!     Start of iterations
!
      Thr=1.0D-7
      MaxIter=100
      Iter=0
!
 777  Continue
!
      Iter=Iter+1
!     Write (6,*)
!     Write (6,*) 'Iteration:', Iter
!     Write (6,*) 'w.r_90:', w,r_90
      Delta_w=0.0001D0*w
      Delta_r=0.0001D0*r_90
!                                                                      *
!***********************************************************************
!                                                                      *
!     The error is a function of r_90 and w
!
!     RMS_s=RMS
!     RMS=RMS*1.0D-15/rBohr  ! bohr
!     Write (6,*) 'RMS:',RMS
!     Xi = 1.0D0/R(w)**2
!     Write (*,*) 'Xi=',Xi
!     RMS=RMS_s
      Do i = 0, 12
         w0 = w    + Facts(1,i) * Delta_w
         r0 = r_90 + Facts(2,i) * Delta_r
!
         Errors(i)= E(r0/R(w0),(r0+T)/R(w0),w0)
!
!        If (i.eq.0) Then
!           Write (*,*) 'r_90,f(x_90)=',r0,f(r0    /R(w0),w0)
!           Write (*,*) 'r_10,f(x_10)=',r0+T,f((r0+T)/R(w0),w0)
!           Write (*,*) 'w,R(w)      =',w0,R(w0)
!        End If
!
      End Do
!     Write (*,*) 'Error=',Errors(0)
!     Call RecPrt('Errors',' ',Errors,13,1)
      g(1)=(Errors(1)-Errors(2))/(2.0D0*Delta_w)
      H(1,1)=(Errors(3)+Errors(4)-2.0D0*Errors(0))/(2.0D0*Delta_w)**2
      g(2)=(Errors(5)-Errors(6))/(2.0D0*Delta_r)
      H(2,2)=(Errors(7)+Errors(8)-2.0D0*Errors(0))/(2.0D0*Delta_r)**2
      H(1,2)=(Errors(9)+Errors(12)-Errors(10)-Errors(11))
     &      / ((2.0D0*Delta_w)*(2.0D0*Delta_r))
      H(2,1)=H(1,2)
!     Call RecPrt('gradient',' ',g,1,2)
!     Call RecPrt('Hessian ',' ',H,2,2)
!
      Call DiagMtrx_x(H,2,iNeg)
!     Write (6,*) 'iNeg=',iNeg
!     Call RecPrt('Hessian ',' ',H,2,2)
      Call MInv(H,HInv,Det,2)
!     Call RecPrt('HInv ',' ',HInv,2,2)
      Step(1)=HInv(1,1)*g(1)+HInv(1,2)*g(2)
      Step(2)=HInv(2,1)*g(1)+HInv(2,2)*g(2)
!     Call RecPrt('Step',' ',Step,1,2)
      Step(1)=Sign(Min(Abs(Step(1)),0.1D0*w),Step(1))
      Step(2)=Sign(Min(Abs(Step(2)),0.1D0*r_90),Step(2))
!     Call RecPrt('Step',' ',Step,1,2)
      w    = w    - Step(1)
      r_90 = r_90 - Step(2)

      If (Iter.lt.MaxIter.and.Errors(0).gt.Thr) Go To 777
      w0=w
      r0=r_90
!     Write (*,*)
!     Write (*,*) 'Iterations:',Iter
!     Write (*,*) 'Error=',Errors(0)
!     Write (*,*) 'r_90,f(x_90)=',r0,f(r0    /R(w0),w0)
!     Write (*,*) 'r_10,f(x_10)=',r0+T,f((r0+T)/R(w0),w0)
!     Write (*,*) 'w,R(w)      =',w0,R(w0)
!     Write (*,*)
!                                                                      *
!***********************************************************************
!                                                                      *
      Xi=1.0D0/R(w)**2
      w=w*Xi
!     Write (*,*) 'Xi:',Xi
!     Write (*,*) ' w:', w
!                                                                      *
!***********************************************************************
!                                                                      *
      Return
! Avoid unused argument warnings
      If (.False.) Call Unused_real(Z)
      End Subroutine ModGauss

      Subroutine DiagMtrx_x(H,nH,iNeg)
      use Constants, only: Zero, One
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer nH, iNeg
      Real*8 H(nH,nH)

      Real*8, Allocatable :: EVal(:), EVec(:,:), Diag(:,:), HU(:,:)
      Real*8 SumHii, Temp
      Integer i, j, ij, ii
!
!     Lu=6
!
      Call mma_allocate(EVal,nH*(nH+1)/2,label='EVal')
      Call mma_allocate(EVec,nH,nH,label='EVec')
!
!---- Copy elements for H
!
      SumHii=Zero
      Do i = 1, nH
         Do j = 1, i
            ij = i*(i-1)/2 + j
            EVal(ij)=H(i,j)
         End Do
         SumHii=SumHii+H(i,i)
      End Do
!     Write (Lu,*) ' SumHii=',SumHii
!
!---- Set up a unit matrix
!
      call dcopy_(nH*nH,[Zero],0,EVec,1)
      call dcopy_(nH,[One],0,EVec,nH+1)
!
!---- Compute eigenvalues and eigenvectors
!
      Call Jacob (EVal,EVec,nH,nH)
      Call Jacord(EVal,EVec,nH,nH)
!
!---- Print out the result
!
      iNeg=0
      Do i = 1, nH
         ii = i*(i+1)/2
         If (EVal(ii).lt.Zero) iNeg=iNeg+1
      End Do
!
      Call mma_allocate(Diag,nH,nH,label='Diag')
      Call mma_allocate(HU,nH,nH,label='HU')
!
      call dcopy_(nH*nH,[Zero],0,Diag,1)
      Do i = 1, nH
         ii=i*(i+1)/2
         temp = EVal(ii)
!        Write (Lu,'(A,G10.4)') 'Hii=',temp
         Diag(i,i)=Max(Abs(temp),1.0D-15)
      End Do
!
      Call DGEMM_('N','N',
     &            nH,nH,nH,
     &            1.0d0,EVec,nH,
     &                  Diag,nH,
     &            0.0d0,HU,nH)
      Call DGEMM_('N','T',
     &            nH,nH,nH,
     &            1.0d0,HU,nH,
     &                  EVec,nH,
     &            0.0d0,H,nH)
!
      Call mma_deallocate(HU)
      Call mma_deallocate(Diag)
      Call mma_deallocate(EVec)
      Call mma_deallocate(EVal)
!
      Return
      End Subroutine DiagMtrx_x
