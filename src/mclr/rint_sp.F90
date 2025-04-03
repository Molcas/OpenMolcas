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
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************

subroutine RInt_SP(rkappa,rmos,rmoa,Focki,Sigma)
!                      ^   ~
! Constructs  F  = <0|[Q  ,H]|0>
!              pq       pq

use MCLR_Data, only: FAMO_SpinP, FAMO_SpinM, SFock, G2mm, G2mp, G2pp, Fp, Fm, G1p, G1m
use MCLR_Data, only: nDensC, nDens2, nMBA
use MCLR_Data, only: rBetaS, rBetaA
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half

implicit none
real*8 rkappa(nDensC), rMOs(*), rmoa(*), Focki(ndens2), Sigma(nDensC)
real*8, allocatable :: MT1(:), MT2(:), MT3(:), Scr(:)

! D,FA used in oit of FA
! p11 -> (~~|  )
! p12 -> (  |~~)
! sign1 k**t=sign K
! sign2=Sign K * Lambda

call mma_allocate(MT1,nmba,Label='MT1')
call mma_allocate(MT2,nmba,Label='MT2')
call mma_allocate(MT3,nmba,Label='MT3')
call mma_allocate(Scr,ndensc,Label='Scr')

call Oit_sp(rkappa,Scr,-1,-One,G2mp,One,Fm,G1m,FAMO_Spinm,MT1,MT2,Focki)
!sigma(:) = rbetaA*Half*SCR(:)
sigma(:) = SCR(:)
call Recprt(' ',' ',SCR,ndensc,1)

! kappa_S

call Oit_sp(rkappa,Scr,1,-One,G2mp,One,Fm,G1m,FAMO_Spinm,MT1,MT2,Focki)
!sigma(:) = sigma(:)-rbetaA*Half*SCR(:)
sigma(:) = sigma(:)-SCR(:)
call Recprt(' ',' ',SCR,ndensc,1)

! alpha_S

if (rbetas /= Zero) then
  call Oit_sp(rkappa,Scr,-1,One,G2pp,One,Fp,G1p,FAMO_Spinp,MT1,MT2,Focki)
  !sigma(:) = sigma(:)+rbetaS*Half*SCR(:)
  sigma(:) = sigma(:)+SCR(:)
  call Recprt(' ',' ',SCR,ndensc,1)
end if
call Oit_sp(rkappa,Scr,-1,One,G2pp,One,G2mm,G1p,Famo_spinp,MT1,MT2,Focki)
!sigma(:) = sigma(:)+ralphas*SCR(:)
call Recprt(' ',' ',SCR,ndensc,1)

call AddGrad_sp(rKappa,Scr,SFock,1,One,One)
call Recprt(' ',' ',SCR,ndensc,1)
sigma(:) = sigma(:)+rbetaA*Half*SCR(:)

MT3(:) = MT1(:)+MT2(:)
call PickMO_MCLR(MT3,rmos,1)

MT3(:) = MT1(:)-MT2(:)
call PickMO_MCLR(MT3,rmoa,1)

call mma_deallocate(Scr)
call mma_deallocate(MT3)
call mma_deallocate(MT2)
call mma_deallocate(MT1)

end subroutine RInt_SP
