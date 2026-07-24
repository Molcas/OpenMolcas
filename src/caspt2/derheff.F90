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

subroutine DerHEff(nConf,nRoots,nState,CLag,VECROT)

use Index_Functions, only: nTri3_Elem
use caspt2_global, only: IDTCEX, LUCIEX
use EQSOLV, only: IVECC, IVECW
use general_data, only: STSYM
use caspt2_module, only: ISCF, JSTATE, MXCI, NASHT
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nConf, nRoots, nState
real(kind=wp), intent(inout) :: CLag(nConf,nRoots)
real(kind=wp), intent(in) :: VECROT(nState)
integer(kind=iwp) :: I, IDCI, IST, JST, NTG1, NTG2, NTG3
real(kind=wp) :: DUMMY(1), OVL
real(kind=wp), allocatable :: CI1(:), CI2(:), CI3(:), DTG1(:), DTG2(:), DTG3(:)

!return

! We evaluate the effective Hamiltonian matrix element in two steps.

NTG1 = NASHT**2
NTG2 = NASHT**4
NTG3 = nTri3_Elem(NTG1)
! Note: Need proper allocation even if unused, sinced allocated
! arrays are in subroutine parameter lists of MKTG3, HCOUP.
NTG1 = max(1,NTG1)
NTG2 = max(1,NTG2)
NTG3 = max(1,NTG3)
call mma_allocate(DTG1,NTG1,Label='DTG1')
call mma_allocate(DTG2,NTG2,Label='DTG2')
call mma_allocate(DTG3,NTG3,Label='DTG3')
DTG1(:) = Zero
DTG2(:) = Zero
DTG3(:) = Zero

!! OVL will contain the derivative contribution?
!! It should be ignored
OVL = Zero
call DerHeffX(IVECW,IVECC,NASHT,NTG3,OVL,DTG1,DTG2,DTG3)

call mma_allocate(CI1,MXCI,Label='MCCI1')
call mma_allocate(CI2,MXCI,Label='MCCI2')
call mma_allocate(CI3,MXCI,Label='MCCI3')

if (ISCF == 0) then
  JST = jState
  CI1(1:NCONF) = Zero
  do I=1,NSTATE
    IDCI = IDTCEX(I)
    IST = I
    if (IST == JST) then
      call DDAFILE(LUCIEX,2,CI2,NCONF,IDCI)
    else if (abs(VECROT(IST)) <= 1.0e-12_wp) then
      call DDAFILE(LUCIEX,0,DUMMY,NCONF,IDCI)
    else
      call DDAFILE(LUCIEX,2,CI3,NCONF,IDCI)
      CI1(1:NCONF) = CI1(1:NCONF)+VECROT(IST)*CI3(1:NCONF)
    end if
  end do

  CI3(1:NCONF) = Zero
  call DERTG3(.true.,STSYM,STSYM,NCONF,NASHT,CI1,CI2,OVL,DTG1,DTG2,NTG3,DTG3,CI3,CLag(1,JST))

  do I=1,NSTATE
    IST = I
    if (IST == JST) then
      cycle
    else if (abs(VECROT(IST)) <= 1.0e-12_wp) then
      cycle
    else
      CLag(1:NCONF,IST) = CLag(1:NCONF,IST)+VECROT(IST)*CI3(1:NCONF)
    end if
  end do
end if

call mma_deallocate(CI1)
call mma_deallocate(CI2)
call mma_deallocate(CI3)

call mma_deallocate(DTG1)
call mma_deallocate(DTG2)
call mma_deallocate(DTG3)

end subroutine DerHEff
