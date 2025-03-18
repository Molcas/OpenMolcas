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
! Copyright (C) 2000, Jonna Stalring                                   *
!***********************************************************************
      Subroutine calcerr(kappa,iestate)
!
!-------------------------------------------------
! Jonna 000411
!
! Calculates the derivative d E_j
!                           ------ = 2\sum_pq Fpq kappa_pq
!                            d w
!
! which estimates the error in the SA.
!---------------------------------------------------
!
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: One, Two
      use MCLR_Data, only: ipMat, nDens2
      use MCLR_Data, only: ISTATE
      use input_mclr, only: nSym,nBas,ntAsh
      Implicit None
      Real*8 kappa(*)
      Integer ieState
#include "SysDef.fh"
      Real*8, Allocatable:: G1q(:), G1r(:), G2r(:),                     &
     &                      T(:), Q(:)
      Integer i, j, itri
      Integer iS, nG1
      Real*8 dejdw
      Real*8, External:: DDot_
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
!
      ng1=itri(ntash,ntash)
!
      Call mma_allocate(G1Q,ng1,Label='G1Q')
      Call mma_allocate(G1R,ntash**2,Label='G1R')
      Call mma_allocate(G2R,ntash**4,Label='G2R')
      Call mma_allocate(T,ndens2,Label='T')
      Call mma_allocate(Q,ndens2,Label='Q')
!
! Get densities of iestate in triangular and rectangular form.
!
      Call rddj(G1R,G1Q,G2R,iestate)
!
! Get Fock matrix ipT
!
      Call FockGen(One,G1R,G2R,T,Q,1)
!
! Print energy
!
!     Call prnte(G1q,T)
!
!     Call Recprt('KAPP',' ',kappa,nDens2,1)
!     Call Recprt('FOCK',' ',T,nDens2,1)
!
      dejdw = 0.0d0
      Do iS=1,nsym
         dejdw = dejdw + ddot_(nBas(is)**2,T(ipmat(is,is)),1,           &
     &                                 kappa(ipmat(is,is)),1)
      End do
      dejdw = Two*dejdw
!     is=1
!     write(*,*)'T(ipmat(is,is))',T(ipmat(is,is))
!     write(*,*)'kappa(ipmat(is,is)+1)',kappa(ipmat(is,is)+1)
!
!     d E_j / d w_i, without constraint
!     Write(6,200) iestate,istate,dejdw
!
      If (iestate.eq.istate) Then
!        d E_i / d w_i, with the constraint sum(w_i)=1,
!        multiplied by 1-w_i... which happens to be equal to
!        d E_i / d w_i, without constraint
!        change sign, because this is the *error*
         Write(6,100) iestate,-dejdw
      End If
100   Format(' **********'/'                 ',                         &
     &       ' Estimated error in the energy of state ',I5,': ',ES12.5/ &
     &       ' **********')
!200   Format(' Derivative of the energy of state ',I5,
!     &       ' with the weight ',I5,': ',ES12.5)
!
      Call mma_deallocate(Q)
      Call mma_deallocate(T)
      Call mma_deallocate(G2R)
      Call mma_deallocate(G1R)
      Call mma_deallocate(G1Q)
!
      End Subroutine calcerr
