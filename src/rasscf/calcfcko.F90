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
! Copyright (C) 2020, Jie J. Bao                                       *
!***********************************************************************

subroutine CalcFckO(CMO,FI,FA,FckO)

use rasscf_global, only: NAC
use general_data, only: NTOT1, NTOT2, NSYM, NASH, NBAS, NFRO, NISH
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
!*****Input
real*8, dimension(NTOT1) :: FI, FA
real*8, dimension(NTOT2) :: CMO
!*****Output
real*8 :: FckO(NAC,NAC)
!*****Auxiliary quantities
integer NB, NA, NI, IOff1, IOff2, IOff3
integer IBas, JBas, ISym, IOrb, JOrb
real*8, allocatable :: FIAAO(:,:), Scr(:,:), FckOt(:,:)

FckO(:,:) = 0.0d0

IOff1 = 0
IOff2 = 1
IOff3 = 0
do ISym=1,NSym
  NB = NBas(ISym)
  NA = NASH(ISym)
  NI = NISH(ISym)+NFro(ISym)
  if (NASH(ISym) > 0) then
    call mma_allocate(FIAAO,nB,nB,Label='FIAAO')
    call mma_allocate(Scr,nB,nA,Label='Scr')
    call mma_allocate(FckOt,NA,NA,Label='FckOt')
    FckOt(:,:) = 0.0d0
    !write(6,*)'Print FI Matrix'
    !do IBas=1,NB
    !  write(6,*) (FI(IOff1+(IBas-1)*IBas/2+JBas),JBas=1,IBas)
    !end do
    !write(6,*)'Print FA Matrix'
    !do IBas=1,NB
    ! write(6,*) ( FA(IOff1+(IBas-1)*IBas/2+JBas),JBas=1,IBas)
    !end do
    !write(6,*) 'Active CMO mat for sym',ISym
    !call RecPrt(' ',' ',CMO(IOff2+NI*NB),NB,NA)
    do IBas=1,NB
      do JBas=1,IBas
        FIAAO(jBas,iBas) = FI(IOff1+(IBas-1)*IBas/2+JBas)+FA(IOff1+(IBas-1)*IBas/2+JBas)
        FIAAO(iBas,jBas) = FIAAO(jBas,iBas)
      end do
    end do
    !write(6,*) 'FIA mat for sym',ISym
    !call RecPrt(' ',' ',FIAAO,NB,NB)
    call DGEMM_('n','n',NB,NA,NB,1.0d0,FIAAO,NB,CMO(IOff2+NI*NB),NB,0.0d0,Scr,NB)
    call DGEMM_('t','n',NA,NA,NB,1.0d0,CMO(IOff2+NI*NB),NB,Scr,NB,0.0d0,FckOt,NA)
    !write(6,*) 'FckO mat for sym',ISym
    !call RecPrt(' ',' ',FckOt,NA,NA)
    do IOrb=1,NA
      do JOrb=1,NA
        FckO(IOrb+IOff3,JOrb+IOff3) = FckOt(jOrb,iOrb)
      end do
    end do
    call mma_deallocate(FIAAO)
    call mma_deallocate(Scr)
    call mma_deallocate(FckOt)
  end if
  IOff1 = IOff1+NB*(NB+1)/2
  IOff2 = IOff2+NB**2
  IOff3 = IOff3+NA
end do

end subroutine CalcFckO
