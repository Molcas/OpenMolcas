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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

!#define _DEBUGPRINT_
subroutine OptClc_X(CInter,nCI,nD,Array,mOV,Ind,MxOptm,kOptim,kOV,LL,DD)

use LnkLst, only: GetNod, GetVec, iVPtr
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nCI, nD, mOV, MxOptm, Ind(MxOptm), kOptim, kOV(2), LL
real(kind=wp), intent(in) :: CInter(nCI,nD)
real(kind=wp), intent(out) :: Array(mOV)
real(kind=wp), optional, intent(out) :: DD
integer(kind=iwp) :: i, iD, iEnd, inode, iSt, ivec
real(kind=wp), allocatable :: Aux(:)

! QNR/DIIS case: compute extrapolated Gradient grd'(n),
! extrapolated Orb Rot Param x'(n), and from this, the
! new, extrapolated displacement direction delta(n)
call mma_allocate(Aux,mOV,Label='Aux')
Aux(:) = Zero

! get last array(n) from LL

call GetVec(Ind(kOptim),LL,inode,Array,mOV)

iEnd = 0
do iD=1,nD
  iSt = iEnd+1
  iEnd = iEnd+kOV(iD)
  Array(iSt:iEnd) = CInter(kOptim,iD)*Array(iSt:iEnd)
end do
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'Initial scaled entities.'
call NrmClc(Array,mOV,'OptClc_X','Array')
write(u6,*)
#endif

do i=1,kOptim-1
  ivec = Ind(i)

  ! get proper gradient from LList.
  call GetNod(ivec,LL,inode)
  if (inode == 0) then
    write(u6,*) 'DIIS: no entry found in LList!'
    call Abend()
  end if
  call iVPtr(Aux,mOV,inode)
  iEnd = 0
  do iD=1,nD
    iSt = iEnd+1
    iEnd = iEnd+kOV(iD)
    Array(iSt:iEnd) = Array(iSt:iEnd)+CInter(i,iD)*Aux(iSt:iEnd)
  end do

end do
#ifdef _DEBUGPRINT_
write(u6,*)
call NrmClc(Array,mOV,'OptClc_X','Array')
write(u6,*)
#endif
if (present(DD)) then
  DD = Zero
  iEnd = 0
  do iD=1,nD
    iSt = iEnd+1
    iEnd = iEnd+kOV(iD)
    DD = DD+sum(Array(iSt:iEnd)**2)
  end do
  DD = sqrt(DD)
end if

call mma_deallocate(Aux)

return

end subroutine OptClc_X
