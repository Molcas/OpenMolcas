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

subroutine CHORAS_DRV(nSym,nBas,nOcc,W_DSQ,W_DLT,W_FLT,ExFac,FSQ,W_CMO)

use Data_Structures, only: Allocate_DSBA, Deallocate_DSBA, DSBA_Type, Integer_Pointer
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nSym, nBas(8)
integer(kind=iwp), target :: nOcc(nSym)
real(kind=wp) :: W_DSQ(*), W_DLT(*), W_FLT(*), ExFac, W_CMO(*)
type(DSBA_Type) :: FSQ
#include "chounit.fh"
#include "choras.fh"
integer(kind=iwp), parameter :: MaxDs = 1
integer(kind=iwp) :: i, iUHF, ja, loff1, MinMem(8), nDen, NumV, rc
real(kind=wp) :: FactC(MaxDs), FactX(MaxDs), Thr, Ymax
logical(kind=iwp) :: DoCoulomb(MaxDs), DoExchange(MaxDs)
type(DSBA_Type) :: CMO, DDec, DLT, DSQ, FLT, MSQ(MaxDs), Vec
type(Integer_Pointer) :: pNocc(1)
integer(kind=iwp), allocatable, target :: nVec(:)
!                                                                      *
!***********************************************************************
!                                                                      *

rc = 0

Lunit(:) = -1

nDen = 1
DoCoulomb(1) = .true.
DoExchange(1) = ExFac /= Zero
FactC(1) = One
FactX(1) = Half*ExFac ! ExFac used for hybrid functionals

call Allocate_DSBA(CMO,nBas,nBas,nSym,Ref=W_CMO)

call Allocate_DSBA(DLT,nBas,nBas,nSym,aCase='TRI',Ref=W_DLT)
call Allocate_DSBA(FLT,nBas,nBas,nSym,aCase='TRI',Ref=W_FLT)

call Allocate_DSBA(DSQ,nBas,nBas,nSym,Ref=W_DSQ)

iUHF = 0

if (DECO) then ! use decomposed density
  ! ==============  Alternative A: Use decomposed density matrix =====
  call mma_allocate(nVec,nSym,Label='nVec')

  ! Allocate vectors representing decomposed density matrix:
  call Allocate_DSBA(Vec,nBas,nBas,nSym)

  ! ------------------------------------------------------------------
  call Allocate_DSBA(Ddec,nBas,nBas,nSym)
  DDec%A0(:) = DSQ%A0(:)
  do i=1,nSym
    ! Loop over symmetries
    if (nBas(i) > 0) then
      Ymax = Zero
      do ja=1,nBas(i)
        Ymax = max(Ymax,DDec%SB(i)%A2(ja,ja))
      end do
      Thr = 1.0e-13_wp*Ymax
      ! Call for decomposition:
      call CD_InCore(DDec%SB(i)%A2,nBas(i),Vec%SB(i)%A2,nBas(i),NumV,Thr,rc)
      if (rc /= 0) call Error(rc)
      nVec(i) = NumV

      if (NumV /= nOcc(i)) then
        write(u6,*) 'Warning! The number of occupied from the decomposition of the density matrix is ',numV,' in symm. ',i
        write(u6,*) 'Expected value = ',nOcc(i)
        write(u6,*) 'Max diagonal of the density in symm. ',i,' is equal to ',Ymax
      end if

    else
      nVec(i) = 0
    end if
  ! End of loop over symmetries
  end do
  call Deallocate_DSBA(DDec)
  ! ------------------------------------------------------------------

  pNocc(1)%I1(1:) => nVec(1:)

  call Allocate_DSBA(MSQ(1),nBas,nBas,nSym,Ref=Vec%A0)

  ! ========End of  Alternative A: Use decomposed density matrix =====
else

  pNocc(1)%I1(1:) => nOcc(1:)

  call Allocate_DSBA(MSQ(1),nBas,nBas,nSym,Ref=CMO%A0)

end if

call CHOSCF_MEM(nSym,nBas,iUHF,DoExchange,pNocc,ALGO,REORD,MinMem,loff1)

! Here follows a long if nest with six combinations:
! ALGO is 1 ,  REORD is .true. or .false., or
! ALGO is 2,  REORD is .true. or .false., DECO is .true.or .false.
if (ALGO == 1) then
  if (REORD) then
    ! (ALGO == 1) .and. REORD:
    call CHO_FOCKTWO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,FactC,FactX,[DLT],[DSQ],[FLT],[FSQ],pNocc,MinMem)
    if (rc /= 0) call Error(rc)

  else
    ! (ALGO == 1) .and. (.not. REORD):
    call CHO_FOCKTWO_RED(rc,nBas,nDen,DoCoulomb,DoExchange,FactC,FactX,[DLT],[DSQ],[FLT],[FSQ],pNocc,MinMem)
    if (rc /= 0) call Error(rc)
  end if

else if (ALGO == 2) then
  if (DECO) then ! use decomposed density

    FactX(1) = Half*ExFac ! vectors are scaled by construction
    if (REORD) then
      ! (ALGO == 2) .and. DECO .and. REORD:
      call CHO_FTWO_MO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,lOff1,FactC,FactX,[DLT],[DSQ],[FLT],[FSQ],MinMem,MSQ,pNocc)
      if (rc /= 0) call Error(rc)

    else
      ! (ALGO == 2) .and. DECO .and. (.not. REORD):
      call CHO_FMO_red(rc,nDen,DoCoulomb,DoExchange,lOff1,FactC,FactX,[DLT],[DSQ],[FLT],[FSQ],MinMem,MSQ,pNocc)
      if (rc /= 0) call Error(rc)
    end if

  else
    if (REORD) then
      ! (ALGO == 2) .and. (.not. DECO) .and. REORD:
      FactX(1) = One*ExFac ! because MOs coeff. are not scaled
      call CHO_FTWO_MO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,lOff1,FactC,FactX,[DLT],[DSQ],[FLT],[FSQ],MinMem,MSQ,pNocc)
      if (rc /= 0) call Error(rc)
    else
      ! (ALGO == 2) .and. (.not. DECO) .and. (.not. REORD):
      FactX(1) = One*ExFac ! because MOs coeff. are not scaled
      call CHO_FMO_red(rc,nDen,DoCoulomb,DoExchange,lOff1,FactC,FactX,[DLT],[DSQ],[FLT],[FSQ],MinMem,MSQ,pNocc)
      if (rc /= 0) call Error(rc)
    end if
  end if

else
  rc = 99
  write(u6,*) 'Illegal Input. Specified Cholesky Algorithm= ',ALGO
  call QUIT(rc)
end if

call CHO_SUM(rc,nSym,nBas,iUHF,DoExchange,[FLT],[FSQ])

if (rc /= 0) call Error(rc)

pNocc(1)%I1 => null()
call Deallocate_DSBA(MSQ(1))
if (DECO) then
  call Deallocate_DSBA(Vec)
  call mma_deallocate(nVec)
end if
call Deallocate_DSBA(DSQ)
call Deallocate_DSBA(DLT)
call Deallocate_DSBA(FLT)
call Deallocate_DSBA(CMO)

return

contains

subroutine Error(rc)

  integer(kind=iwp), intent(in) :: rc

  write(u6,*) 'CHORAS_DRV. Non-zero return code. rc= ',rc
  call QUIT(rc)

end subroutine Error

end subroutine CHORAS_DRV
