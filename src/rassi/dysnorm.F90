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
! Copyright (C) 2020, Bruno Tenorio                                    *
!***********************************************************************
!  SUBROUTINE DYSNORM
!  PURPOSE: CALCULATE CORRECTED DYSON NORMS FOR CI EXPANSIONS IN A
!  BIORTH. BASIS
!***********************************************************************

subroutine DYSNORM(CMOA,DYSCMO,NORM)

use Symmetry_Info, only: nIrrep
use Index_Functions, only: nTri_Elem
use rassi_data, only: NBASF, NCMO, NOSH, NOSHT
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: CMOA(*), DYSCMO(*)
real(kind=wp), intent(out) :: NORM
integer(kind=iwp) :: IC, iOff1, iOff2, iOpt, iRc, ist, ista, istacc(8), istao(8), istc, istca, istcb, istcc, istcmo(8), isy1, &
                     iSyLbl, isym, nb, nb1, nbast1, nbast2, ndys, no1, nscr
real(kind=wp) :: NORMSCR
character(len=8) :: Label
real(kind=wp), allocatable :: Dysab(:), Dysab2(:), IAO(:), SAO(:), Scr(:), Scr2(:)
real(kind=wp), external :: DDOT_

!============================================================
nbast1 = sum(nTri_Elem(NBASF(1:nIrrep)))
nbast2 = sum(NBASF(1:nIrrep)**2)

call mma_allocate(DYSAB,NOSHT)
DYSAB(:) = DYSCMO(1:NOSHT)

call mma_allocate(SAO,NBAST1)
call mma_allocate(IAO,NBAST2)
IAO = Zero

iRc = -1
iOpt = 6
IC = 1
iSyLbl = 1
Label = 'Mltpl  0'
call RdOne(iRc,iOpt,Label,IC,SAO,iSyLbl)
iOff1 = 0
iOff2 = 0
do iSym=1,nIrrep
  nb = nBasf(iSym)
  if (nb > 0) call Square(SAO(1+iOff1),IAO(iOff2+1),1,nb,nb)
  iOff1 = iOff1+nTri_Elem(nb)
  iOff2 = iOff2+nb**2
end do
call mma_deallocate(SAO)

!============================================================

IST = 1
ISTA = 1
ISTCC = 1
do ISY1=1,nIrrep
  ISTCMO(ISY1) = IST
  ISTAO(ISY1) = ISTA
  ISTACC(ISY1) = ISTCC
  NO1 = NOSH(ISY1)
  IST = IST+NO1*NBASF(ISY1)
  ISTA = ISTA+(NBASF(ISY1)*NBASF(ISY1))
  ISTCC = ISTCC+NO1**2
end do

norm = Zero
NSCR = NCMO
NDYS = 1
do ISY1=1,nIrrep
  ISTCB = ISTCMO(ISY1)
  ISTCA = ISTAO(ISY1)
  ISTC = ISTACC(ISY1)
  NO1 = NOSH(ISY1)
  NB1 = NBASF(ISY1)
  if (NB1*NO1 == 0) cycle

  call mma_allocate(scr,nscr)
  call mma_allocate(scr2,nscr)

  call DGEMM_('N','N',NB1,NO1,NB1,One,IAO(ISTCA),NB1,CMOA(ISTCB),NB1,Zero,Scr(ISTCB),NB1)

  call DGEMM_('T','N',NO1,NO1,NB1,One,CMOA(ISTCB),NB1,Scr(ISTCB),NB1,Zero,Scr2(ISTC),NO1)

  ! Src2 is the M matrix
  ! Proceed to compute norm=Dab*(Dab*M)
  ! Where Dab*M=Dysab2

  call mma_allocate(DYSAB2,NO1)
  call DGEMV_('N',NO1,NO1,One,Scr2(ISTC),NO1,DYSAB(NDYS),1,Zero,DYSAB2,1)

  normscr = DDOT_(NO1,DYSAB(NDYS),1,DYSAB2,1)
  norm = norm+normscr
  call mma_deallocate(DYSAB2)

  NDYS = NDYS+NO1
  call mma_deallocate(Scr)
  call mma_deallocate(Scr2)
end do

call mma_deallocate(IAO)
call mma_deallocate(DYSAB)

end subroutine DYSNORM
