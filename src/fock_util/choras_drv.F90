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

use Fock_util_global, only: ALGO, Deco, Lunit, REORD
use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type, Integer_Pointer
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(8)
integer(kind=iwp), target, intent(in) :: nOcc(nSym)
real(kind=wp), intent(in) :: W_DSQ(*), W_DLT(*), ExFac, W_CMO(*)
real(kind=wp), intent(inout) :: W_FLT(*)
type(DSBA_Type), intent(inout) :: FSQ(1)
integer(kind=iwp), parameter :: MaxDs = 1
integer(kind=iwp) :: i, nD, ja, loff1, MinMem(8), nDen, NumV, rc
real(kind=wp) :: FactC(MaxDs), FactX(MaxDs), Thr, Ymax
logical(kind=iwp) :: DoCoulomb(MaxDs), DoExchange(MaxDs)
type(DSBA_Type) :: DDec, DLT(1), DSQ(1), FLT(1), MSQ(MaxDs), Vec
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

call Allocate_DT(DLT(1),nBas,nBas,nSym,aCase='TRI',Ref=W_DLT)
call Allocate_DT(FLT(1),nBas,nBas,nSym,aCase='TRI',Ref=W_FLT)

call Allocate_DT(DSQ(1),nBas,nBas,nSym,Ref=W_DSQ)

nD = 1

if (DECO) then ! use decomposed density
  ! ==============  Alternative A: Use decomposed density matrix =====
  call mma_allocate(nVec,nSym,Label='nVec')

  ! Allocate vectors representing decomposed density matrix:
  call Allocate_DT(Vec,nBas,nBas,nSym)

  ! ------------------------------------------------------------------
  call Allocate_DT(Ddec,nBas,nBas,nSym)
  DDec%A0(:) = DSQ(1)%A0(:)
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
  call Deallocate_DT(DDec)
  ! ------------------------------------------------------------------

  pNocc(1)%I1(1:) => nVec(1:)

  call Allocate_DT(MSQ(1),nBas,nBas,nSym,Ref=Vec%A0)

  ! ========End of  Alternative A: Use decomposed density matrix =====
else

  pNocc(1)%I1(1:) => nOcc(1:)

  call Allocate_DT(MSQ(1),nBas,nBas,nSym,Ref=W_CMO)

end if

call CHOSCF_MEM(nSym,nBas,nD,DoExchange,pNocc,ALGO,REORD,MinMem,loff1)

! Here follows a long if nest with six combinations:
! ALGO is 1,  REORD is .true. or .false., or
! ALGO is 2,  REORD is .true. or .false., DECO is .true.or .false.
if (ALGO == 1) then
  if (REORD) then
    ! (ALGO == 1) .and. REORD:
    call CHO_FOCKTWO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,FactC,FactX,DLT,DSQ,FLT,FSQ,pNocc,MinMem)
    if (rc /= 0) call Error(rc)

  else
    ! (ALGO == 1) .and. (.not. REORD):
    call CHO_FOCKTWO_RED(rc,nBas,nDen,DoCoulomb,DoExchange,FactC,FactX,DLT,DSQ,FLT,FSQ,pNocc,MinMem)
    if (rc /= 0) call Error(rc)
  end if

else if (ALGO == 2) then
  if (DECO) then ! use decomposed density

    FactX(1) = Half*ExFac ! vectors are scaled by construction
    if (REORD) then
      ! (ALGO == 2) .and. DECO .and. REORD:
      call CHO_FTWO_MO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,lOff1,FactC,FactX,DLT,DSQ,FLT,FSQ,MinMem,MSQ,pNocc)
      if (rc /= 0) call Error(rc)

    else
      ! (ALGO == 2) .and. DECO .and. (.not. REORD):
      call CHO_FMO_red(rc,nDen,DoCoulomb,DoExchange,lOff1,FactC,FactX,DLT,DSQ,FLT,FSQ,MinMem,MSQ,pNocc)
      if (rc /= 0) call Error(rc)
    end if

  else
    if (REORD) then
      ! (ALGO == 2) .and. (.not. DECO) .and. REORD:
      FactX(1) = One*ExFac ! because MOs coeff. are not scaled
      call CHO_FTWO_MO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,lOff1,FactC,FactX,DLT,DSQ,FLT,FSQ,MinMem,MSQ,pNocc)
      if (rc /= 0) call Error(rc)
    else
      ! (ALGO == 2) .and. (.not. DECO) .and. (.not. REORD):
      FactX(1) = One*ExFac ! because MOs coeff. are not scaled
      call CHO_FMO_red(rc,nDen,DoCoulomb,DoExchange,lOff1,FactC,FactX,DLT,DSQ,FLT,FSQ,MinMem,MSQ,pNocc)
      if (rc /= 0) call Error(rc)
    end if
  end if

else
  rc = 99
  write(u6,*) 'Illegal Input. Specified Cholesky Algorithm= ',ALGO
  call QUIT(rc)
end if

call CHO_SUM(rc,nSym,nBas,nD,DoExchange,FLT,FSQ)

if (rc /= 0) call Error(rc)

pNocc(1)%I1 => null()
call Deallocate_DT(MSQ(1))
if (DECO) then
  call Deallocate_DT(Vec)
  call mma_deallocate(nVec)
end if
call Deallocate_DT(DSQ(1))
call Deallocate_DT(DLT(1))
call Deallocate_DT(FLT(1))

return

contains

subroutine Error(rc)

  integer(kind=iwp), intent(in) :: rc

  write(u6,*) 'CHORAS_DRV. Non-zero return code. rc= ',rc
  call QUIT(rc)

end subroutine Error

end subroutine CHORAS_DRV
