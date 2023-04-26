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

subroutine FockTwo_Drv(nSym,nBas,nAux,Keep,DLT,DSQ,FLT,nFLT,ExFac,nBMX)

use Fock_util_interface, only: CHOras_drv
use Fock_util_global, only: ALGO
use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(8), nAux(8), Keep(8), nFLT, nBMX
real(kind=wp), intent(in) :: DLT(*), DSQ(*), ExFac
real(kind=wp), intent(inout) :: FLT(nFLT)
integer(kind=iwp) :: LBUF
real(kind=wp) :: CMO_DUMMY(1)
logical(kind=iwp) :: DoCholesky, GenInt
type(DSBA_Type) :: WFSQ(1)
real(kind=wp), allocatable :: Temp(:), W1(:), W2(:)

!                                                                      *
!***********************************************************************
!                                                                      *
! nAux is the number of occupied orbitals
GenInt = .false.
DoCholesky = .false.
if (ALGO == 0) GenInt = .true. ! use GenInt to regenerate integrals

call DecideOnCholesky(DoCholesky)

call Allocate_DT(WFSQ(1),nBas,nBas,nSym)
WFSQ(1)%A0(:) = Zero

if ((.not. DoCholesky) .or. GenInt) call mma_allocate(W2,NBMX**2,Label='W2')

call mma_allocate(Temp,nFlt,Label='Temp')
Temp(:) = Zero

call mma_maxDBLE(LBUF)

!                                                                      *
!***********************************************************************
!                                                                      *
! Standard building of the Fock matrix from two-el integrals
! OR
! Building of the Fock matrix regenerating the two-el integrals on the fly from CD/RI vectors

if ((.not. DoCholesky) .or. GenInt) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (.not. DoCholesky) then
    call mma_allocate(W1,LBUF,Label='W1')
  else
    LBUF = max(LBUF-LBUF/10,0) ! save some space for Gen_Int
    call mma_allocate(W1,LBUF,Label='W1')
  endif

  if (LBUF < 1+NBMX**2) then
    write(u6,*) ' FockTwo_Drv Error: Too little memory remains for the call to FOCKTWO.'
    write(u6,*) ' Largest allocatable array size LBUF=',LBUF
    write(u6,*) ' Max nr of bf in any symmetry,  NBMX=',NBMX
    write(u6,*) ' Required minimum size     1+NBMX**2=',1+NBMX**2
    write(u6,*) '    (All in Real*8-size words)'
    call ABEND()
  end if

  call FOCKTWO(nSym,nBas,nAux,Keep,DLT,DSQ,Temp,nFlt,WFSQ(1)%A0,LBUF,W1,W2,ExFac)

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !   Building of the Fock matrix directly from Cholesky vectors

else if (DoCholesky .and. (.not. GenInt)) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! CMO_DUMMY is required call argument of choras_drv:
  ! (Not used, see logical flags in choras_drv)
  ! subroutine choras_drv(nSym,nBas,nOcc,DSQ,DLT,FLT,ExFac,WFSQ,CMO)
  call CHOras_drv(nSym,nBas,nAux,DSQ,DLT,Temp,ExFac,WFSQ,CMO_DUMMY)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call DaXpY_(nFlt,One,Temp,1,FLT,1)

call mma_deallocate(Temp)
if (allocated(W1)) call mma_deallocate(W1)
if (allocated(W2)) call mma_deallocate(W2)
call Deallocate_DT(WFSQ(1))

return

end subroutine FockTwo_Drv
