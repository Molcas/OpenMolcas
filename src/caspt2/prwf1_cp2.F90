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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine PRWF1_CP2(ISYCI,CI,mCI,THR)

use CI_Lib_Info, only: CI_Lib_Space
use CI_Lib_States, only: CI_Walk_Next
use caspt2_module, only: PRSD
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ISYCI, mCI
real(kind=wp), intent(in) :: CI(mCI), THR
integer(kind=iwp) :: ICONF, IMS, ISY, it(4), itc(4), K, LENCSF, LEV
real(kind=wp) :: COEF
character(len=256) :: LINE
integer(kind=iwp), allocatable :: ICS(:), LEX(:)
character, parameter :: CODE(0:3) = ['0','u','d','2']
integer(kind=iwp), parameter :: istate = 1

! The CSFs come from the CI interface in the order of the CI vector, with
! their split-graph position (midvertex, upper symmetry, upper and lower
! walk) and step vector.

call mma_allocate(ICS,CI_Lib_Space(iState)%nLev,Label='ICS')

! SVC: set up a CSF string length as LENCSF
LINE = ' '
LENCSF = 0
ISY = 0
do LEV=1,CI_Lib_Space(iState)%nLev
  if (ISY /= CI_Lib_Space(istate)%ISm(LEV)) then
    ISY = CI_Lib_Space(istate)%ISm(LEV)
    LENCSF = LENCSF+1
  end if
  LENCSF = LENCSF+1
end do
LENCSF = min(LENCSF,256)
LENCSF = max(LENCSF,10)

! Size of occup/spin coupling part of line:
write(u6,*) ' Occupation of active orbitals, and spin coupling'
write(u6,*) ' of open shells. (u,d: Spin up or down).'
write(u6,*) ' SGUGA info is (Midvert:IsyUp:UpperWalk/LowerWalk)'
LINE(1:10) = 'Occupation'
write(u6,'(2X,A10,2X,A16,2X,A,2(2X,A13))') 'Conf','SGUGA info      ',LINE(1:LENCSF),'Coefficient','Weight'

! SVC2010:
! allocate scratch memory for determinant expansion
if (PRSD) call mma_allocate(LEX,CI_Lib_Space(iState)%nLev,LABEL='LEX')

LINE = ''

it(:) = 1
do ICONF=1,mCI
  if (it(1) == 0) exit
  itc(:) = it(:)
  call CI_Walk_Next(istate,ISYCI,it,ICS)
  COEF = CI(ICONF)
  ! SKIP OR PRINT IT OUT?
  if (abs(COEF) < THR) cycle
  ! PRINT IT!
  K = 0
  ISY = 0
  do LEV=1,CI_Lib_Space(iState)%nLev
    if (ISY /= CI_Lib_Space(istate)%ISm(LEV)) then
      ISY = CI_Lib_Space(istate)%ISm(LEV)
      K = K+1
      LINE(K:K) = ' '
    end if
    K = K+1
    LINE(K:K) = CODE(ICS(LEV))
  end do
  write(u6,'(2X,I10,2X,"(",I2,":",I1,":",I4,"/",I4,")",2X,A,2(2X,F13.6))') ICONF,itc(1),itc(2),itc(4),itc(3),LINE(1:LENCSF), &
                                                                          COEF,COEF**2
  ! SVC2010 experimental: add determinant expansion
  if (PRSD) then
    ! Specify projected spin in half integer units
    ! Default: use maximum spin projection
    IMS = CI_Lib_Space(iState)%iSpin-1
    write(u6,*)
    call EXPCSF(ICS,CI_Lib_Space(iState)%nLev,IMS,LEX,coef,0)
    write(u6,*)
  end if
end do

call mma_deallocate(ICS)

! SVC2010: free scratch for determinant expansion
if (PRSD) call mma_deallocate(LEX)
write(u6,*)

end subroutine PRWF1_CP2
