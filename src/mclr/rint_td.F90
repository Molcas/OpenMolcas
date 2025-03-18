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
! Copyright (C) Jonna Stalring                                         *
!***********************************************************************
      SubRoutine RInt_td(ekappa,mkappa,isym)
!
!
!******************************************************
!            [2]                                      *
! Modifies  E    *kappa by subtracting Omega*S*kappa  *
! to make timedep PT possible.                        *
!                                                     *
!******************************************************
!
!     mkappa       is d/dx(kappa) in matrix form (derivative of orb rot mat w r t PT)
!     ekappa       is E*d/dx(kappa) in matrix form, REPLACED BY
!                  the entire timedep expression
!     isym         Rinttd is called once for each symmtry
!
! Local variables
!
!     dens The density matrix
!     wDKt   Omega*(density matrix)*(kappa transposed)
!     wKtD   As above but different order
!
!
      use Arrays, only: G1t
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: Zero, Two
      use MCLR_Data, only: nDens2,ipCM,ipMat,nA
      use input_mclr, only: Omega,nSym,nAsh,nBas,nIsh
      Implicit None
      Real*8 ekappa(ndens2),mkappa(ndens2)
      Integer iSym
!
      Real*8, Allocatable:: Dens(:), wDKt(:), wKtD(:)
      Integer lDens,iS,ip3,Inc,iB,jB,ip,iA,jA,ip2,jS,IncX,Incy,Length

      integer i,j,itri
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
!
!---------------------------------------------------------
!
!  Allocate memory
!
      lDens = 0
      Do iS=1,nSym
         lDens = lDens + nBas(iS)**2
      End Do
      Call mma_allocate(Dens,lDens,Label='Dens')
      Call mma_allocate(wDKt,ndens2,Label='wDKt')
      Call mma_allocate(wKtD,ndens2,Label='wKtD')
!
!---------------------------------------------------------
!***************************************
! Construct the density matrix
!***************************************
!
! For HF
!
      Dens(:)=Zero
      ip3 = 1
      Do iS=1,nSym
          inc = nBas(iS) + 1
          call dcopy_(nIsh(iS),[2.0d0],0,Dens(ip3),inc)
          ip3 = ip3 + nBas(iS)*nBas(iS)
      End Do
!
!         do is=1,nsym
!            size=nbas(is)*nbas(is)      ! The size of the sym i and j martix
!            do k=-1,(size-1)             ! Loop over the entire submatrix
!               residual=DMOD((k+1),(nbas(is)+1))
!               if ((residual.eq.zero).and.(k.lt.(nbas(is)*nIsh
!     &      (is))).and.(nIsh(is).ne.0)) then
!                     Dens(1+ipCM(is)+k)=Two
!               else
!                     Dens(1+ipCM(is)+k)=Zero
!               end if
!            end do
!          end do
!
! For a CASSCF wavefunc. From Anders subrut r2elint
! Add the active active dens
!
          Do iS=1,nSym
             Do iB=1,nAsh(iS)
                Do jB=1,nAsh(iS)
                   ip=ipCM(iS)+ib+nIsh(is)+(jB+nIsh(is)-1)*nBas(is)-1
                   iA=nA(is)+ib
                   jA=nA(is)+jb
                   ip2=itri(iA,jA)
                   Dens(ip)=G1t(ip2)
                End Do
             End Do
          End Do
!*******************************************
! Multiply D and mkappa wDKt and wKtD
!*******************************************
!      Call RecPrt('dens ',' ',dens,ldens,1)
!      Call RecPrt('ekappa ',' ',ekappa,ndens2,1)
      do is=1,nsym
           js=iEOr(is-1,isym-1)+1
!      wDKt
           If ( (nBas(iS).gt.0).and.(nBas(jS).gt.0) ) Then
             call DGEMM_('n','n',nbas(is),nbas(js),nbas(is),            &
     &                   Two*Omega,Dens(ipCM(is)),nbas(is),             &
     &                             mkappa(ipmat(is,js)),nbas(is),       &
     &                   Zero,wDKt(ipmat(is,js)),nbas(is))
!      wKtD
             call DGEMM_('n','n',nbas(is),nbas(js),nbas(js),            &
     &                   Two*Omega,mkappa(ipmat(is,js)),nbas(is),       &
     &                             Dens(ipCM(js)),nbas(js),             &
     &                   Zero,wKtD(ipmat(is,js)),nbas(is))

!*****************************************************
!            Replace ekappa ekappa=ekappa-wDKt+wKtD
!*****************************************************
             incx=1
             incy=incx
             length=nbas(is)*nbas(js)
             call daxpy_(length,1.0d0,wDKt(ipmat(is,js)),               &
     &            incx,ekappa(ipmat(is,js)),incy)
             call daxpy_(length,-1.0d0,wKtD(ipmat(is,js)),              &
     &            incx,ekappa(ipmat(is,js)),incy)
          End If
      end do
!      Call RecPrt('wDKt ',' ',wDKt,ndens2,1)
!      Call RecPrt('wKtD ',' ',wKtD,ndens2,1)
!      Call RecPrt('ekappa ',' ',ekappa,ndens2,1)
!
!  Free memory
!
       Call mma_deallocate(Dens)
       Call mma_deallocate(wDKt)
       Call mma_deallocate(wKtD)
      end SubRoutine RInt_td
