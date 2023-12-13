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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************
!*************************************************************
! Initializes info needed by the Cholesky integral generator
!
! F. Aquilante
!*************************************************************

subroutine INIT_GETINT(RC)

use Index_Functions, only: nTri_Elem
use GetInt_mod, only: LuCVec, mNeed, nBas, nPQ, nRS, NumCho, nVec, pq1, Vec2
use RICD_Info, only: Do_DCCD
use TwoDat, only: rcTwo
use stdalloc, only: mma_allocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: RC
integer(kind=iwp) :: LWork, nSym

rc = 0
call get_iscalar('nSym',nSym)
call get_iarray('nBas',nBas,nSym)
call INIT_NumCV(NumCho,nSym)

if (Do_DCCD) then
  if (NumCho(1) < 1) then
    write(u6,*) 'Init_GetInt: NumCho(1) < 1'
    call Abend()
  end if

  nPQ = nTri_Elem(nBas(1))
  nRS = nPQ
  ! Memory management
  mNeed = 2*nPQ

  if (mNeed <= 0) then
    ! ***QUIT*** bad initialization
    write(u6,*) 'Gen_Int: bad initialization'
    rc = rcTwo%RD11
    call Abend()
  end if

  ! Set up the batch procedure
  ! --------------------------
  call mma_maxDBLE(LWORK)
  LWORK = LWORK-LWORK/10

  nVec = min(LWORK/mNeed,NumCho(1))

  if (nVec <= 0) then
    ! ***QUIT*** insufficient memory
    write(u6,*) 'Gen_Int: Insufficient memory for batch'
    write(u6,*) 'LWORK= ',LWORK
    write(u6,*) 'mNeed= ',mNeed
    write(u6,*) 'NumCho= ',NumCho(1)
    rc = rcTwo%RD05
    call Abend()
  end if

  ! Allocate memory for reading the vectors and do the transposition
  call mma_allocate(Vec2,Npq,nVec,label='MemC2')

end if

LuCVec(1) = -1
LuCVec(2) = -1

pq1 = 0

return

end subroutine INIT_GETINT
