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

use Index_Functions, only: iTri, nTri_Elem
use rasscf_global, only: NAC
use general_data, only: NASH, NBAS, NFRO, NISH, NSYM, NTOT1, NTOT2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: CMO(NTOT2), FI(NTOT1), FA(NTOT1)
real(kind=wp), intent(out) :: FckO(NAC,NAC)
integer(kind=iwp) :: IBas, IOff1, IOff2, IOff3, IOrb, ISym, JBas, NA, NB, NI
real(kind=wp), allocatable :: FIAAO(:,:), FckOt(:,:), Scr(:,:)

FckO(:,:) = Zero

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
    FckOt(:,:) = Zero
    !write(u6,*)'Print FI Matrix'
    !do IBas=1,NB
    !  write(u6,*) (FI(IOff1+iTri(IBas,JBas)),JBas=1,IBas)
    !end do
    !write(u6,*)'Print FA Matrix'
    !do IBas=1,NB
    ! write(u6,*) ( FA(IOff1+iTri(IBas,JBas)),JBas=1,IBas)
    !end do
    !write(u6,*) 'Active CMO mat for sym',ISym
    !call RecPrt(' ',' ',CMO(IOff2+NI*NB),NB,NA)
    do IBas=1,NB
      do JBas=1,IBas
        FIAAO(jBas,iBas) = FI(IOff1+iTri(IBas,JBas))+FA(IOff1+iTri(IBas,JBas))
        FIAAO(iBas,jBas) = FIAAO(jBas,iBas)
      end do
    end do
    !write(u6,*) 'FIA mat for sym',ISym
    !call RecPrt(' ',' ',FIAAO,NB,NB)
    call DGEMM_('n','n',NB,NA,NB,One,FIAAO,NB,CMO(IOff2+NI*NB),NB,Zero,Scr,NB)
    call DGEMM_('t','n',NA,NA,NB,One,CMO(IOff2+NI*NB),NB,Scr,NB,Zero,FckOt,NA)
    !write(u6,*) 'FckO mat for sym',ISym
    !call RecPrt(' ',' ',FckOt,NA,NA)
    do IOrb=1,NA
      FckO(IOrb+IOff3,IOff3+1:IOff3+NA) = FckOt(:,iOrb)
    end do
    call mma_deallocate(FIAAO)
    call mma_deallocate(Scr)
    call mma_deallocate(FckOt)
  end if
  IOff1 = IOff1+nTri_Elem(NB)
  IOff2 = IOff2+NB**2
  IOff3 = IOff3+NA
end do

end subroutine CalcFckO
