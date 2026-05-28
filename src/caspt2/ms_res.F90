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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine MS_Res(MODE,IST,JST,Scal)
! Compute the derivative of E^PT2 with respct to the T amplitude

use Index_Functions, only: nTri3_Elem
use EQSOLV, only: IVECC, IVECC2, IVECW
use caspt2_global, only: IDTCEX, LUCIEX
use caspt2_module, only: ISCF, MXCI, NASHT, NCONF, NSTATE, STSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: MODE, IST, JST
real(kind=wp), intent(in) :: Scal
integer(kind=iwp) :: I, IDCI, NTG1, NTG2, NTG3
real(kind=wp) :: DUMMY(1), DVALUE, OVL
real(kind=wp), allocatable :: CI1(:), CI2(:), TG1(:), TG2(:), TG3(:)

! We evaluate the effective Hamiltonian matrix element in two steps.

NTG1 = NASHT**2
NTG2 = NASHT**4
NTG3 = nTri3_Elem(NTG1)
! Note: Need proper allocation even if unused, sinced allocated
! arrays are in subroutine parameter lists of MKTG3, HCOUP.
NTG1 = max(1,NTG1)
NTG2 = max(1,NTG2)
NTG3 = max(1,NTG3)
call mma_allocate(TG1,NTG1,Label='TG1')
call mma_allocate(TG2,NTG2,Label='TG2')
call mma_allocate(TG3,NTG3,Label='TG3')
TG1(:) = Zero
TG2(:) = Zero
TG3(:) = Zero

call mma_allocate(CI1,MXCI,Label='MCCI1')
call mma_allocate(CI2,MXCI,Label='MCCI2')
if (ISCF == 0) then
  ! Read root vectors nr. IST and JST from LUCI.
  do I=1,NSTATE
    IDCI = IDTCEX(I)
    if (I == IST) then
      call DDAFILE(LUCIEX,2,CI1,NCONF,IDCI)
      if (I == JST) CI2(1:NCONF) = CI1(1:NCONF)
    else if (I == JST) then
      call DDAFILE(LUCIEX,2,CI2,NCONF,IDCI)
    else
      call DDAFILE(LUCIEX,0,DUMMY,NCONF,IDCI)
    end if
  end do
end if

call MKTG3(STSYM,STSYM,CI1,CI2,OVL,TG1,TG2,NTG3,TG3)
call mma_deallocate(CI1)
call mma_deallocate(CI2)

!! Do similar to RHS_STRANS. Multiply the solution vector (T) with
!! the overlap-like term constructe with transition density
!! matrices. The output is IVECC (MODE=1) or IVECC2 (MODE=2)
if (MODE == 1) then
  call MS_STRANS(IVECW,IVECC,NASHT,NTG3,OVL,TG1,TG2,TG3,DVALUE,SCAL)
else if (MODE == 2) then
  call MS_STRANS(IVECC,IVECC2,NASHT,NTG3,OVL,TG1,TG2,TG3,DVALUE,SCAL)
else if (MODE == 3) then
  call MS_STRANS(IVECW,IVECC,NASHT,NTG3,OVL,TG1,TG2,TG3,DVALUE,SCAL)
end if

call mma_deallocate(TG1)
call mma_deallocate(TG2)
call mma_deallocate(TG3)

end subroutine MS_Res
