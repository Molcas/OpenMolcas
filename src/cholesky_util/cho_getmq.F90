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

subroutine Cho_GetMQ(MQ,l_MQ,List_QShp,nQShp)

use Cholesky, only: iiBstR, iiBstRSh, IndRed, IndRSh, iQuAB, LuSel, nnBstR, nnBstRSh, nnShl, nQual, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: l_MQ, nQShp, List_QShp(nQShp)
real(kind=wp), intent(out) :: MQ(l_MQ)
integer(kind=iwp) :: iAB, iabSh, iAdr, iL_Shp, iL_ShpG, ipfr, ipS, ipto, iQ, iQoff, isAB, iShAB, iShABG, iShp, iShpAdr, jQ, jSym, &
                     Lint, Lread
integer(kind=iwp), allocatable :: kOff_Shp(:)
real(kind=wp), allocatable :: Scr(:)
integer(kind=iwp), parameter :: iOpt = 2
integer(kind=iwp), external :: Cho_F2SP, Cho_P_LocalSP

if (sum(nQual(1:nSym)) < 1) return ! this test makes sense for parallel run

call mma_allocate(kOff_Shp,nnShl,Label='kOff_Shp')

iQoff = 0
do jSym=1,nSym

  if (nQual(jSym) < 1) cycle ! next symmetry

  Lint = 0
  do iShp=1,nQShp  ! set only the needed offsets
    iL_ShpG = List_QShp(iShp) ! Shell pair
    iL_Shp = Cho_P_LocalSP(iL_ShpG) ! local shell pair
    kOff_Shp(iL_Shp) = Lint
    Lint = Lint+nnBstRSh(jSym,iL_Shp,2)
  end do

  call mma_allocate(Scr,Lint,Label='Scr')

  ! Read the integrals
  !-------------------
  do jQ=1,nQual(jSym)

    iAdr = nnBstr(jSym,2)*(jQ-1)

    do iShp=1,nQShp

      iL_ShpG = List_QShp(iShp)
      iL_Shp = Cho_P_LocalSP(iL_ShpG) ! local shell pair
      ipS = 1+kOff_Shp(iL_Shp)
      iShpAdr = iAdr+iiBstRSh(jSym,iL_Shp,2)
      Lread = nnBstRSh(Jsym,iL_Shp,2)
      call dDaFile(LuSel(JSym),iOpt,Scr(ipS),Lread,iShpAdr)

    end do

    ! Extract the matrix of qualified integrals
    !------------------------------------------
    do iQ=1,nQual(jSym)

      iAB = iQuAB(iQ,jSym)  ! addr curr red set
      isAB = iAB-iiBstR(jSym,2)  ! symm. reduction
      iShABG = IndRsh(IndRed(iAB,2)) ! glob. SP it belongs to
      iShAB = Cho_P_LocalSP(Cho_F2SP(iShABG)) ! local SP
      iabSh = isAB-iiBstRSh(jSym,iShAB,2) ! addr within Shp
      ipfr = kOff_Shp(iShAB)+iabSh
      ipto = iQoff+nQual(jSym)*(jQ-1)+iQ
      MQ(ipto) = Scr(ipfr)

    end do

  end do

  iQoff = iQoff+nQual(jSym)**2

  call mma_deallocate(Scr)

end do

call mma_deallocate(kOff_Shp)

return

end subroutine Cho_GetMQ
