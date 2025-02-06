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

subroutine SPSPCLS_GAS(NOCTPA,NOCTPB,IOCA,IOCB,NELFGP,MXPNGAS,NGAS,ISPSPCLS,ICLS,NCLS,IPRNT)
! Obtain mapping a-supergroup X b-supergroup => class
!
! =====
! Input
! =====
!
! NOCTPA : Number of alpha types
! NOCTPB : Number of beta types
!
! IOCA(IGAS,ISTR) occupation of AS IGAS for alpha string type ISTR
! IOCB(IGAS,ISTR) occupation of AS IGAS for beta  string type ISTR
!
! MXPNGAS : Largest allowed number of gas spaces
! NGAS    : Actual number of gas spaces
!
! ======
! Output
! ======
!
! ISPSPCLS(IATP,IBTP) => Class of this block of determinants
!                        =0 indicates unallowed(class less) combination

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: NOCTPA, NOCTPB, MXPNGAS, IOCA(MXPNGAS,NOCTPA), IOCB(MXPNGAS,NOCTPB), NELFGP(*), NGAS, &
                     ISPSPCLS(NOCTPA,NOCTPB), NCLS, ICLS(NGAS,NCLS), IPRNT
integer(kind=iwp) :: IAMOKAY, IATP, IBTP, IEL, IGAS, IICLS, KCLS, NTEST

NTEST = 0
NTEST = max(NTEST,IPRNT)
if (NTEST >= 10) then
  write(u6,*) ' ISPSPCLS_GAS entered'
  write(u6,*) ' ===================='
  write(u6,*)
  write(u6,*) ' IOCA and IOCB'
  call IWRTMA(IOCA,NGAS,NOCTPA,MXPNGAS,NGAS)
  call IWRTMA(IOCB,NGAS,NOCTPB,MXPNGAS,NGAS)
  write(u6,*)
  write(u6,*) ' ICLS'
  call IWRTMA(ICLS,NGAS,NCLS,NGAS,NCLS)
end if

do IATP=1,NOCTPA
  do IBTP=1,NOCTPB
    IICLS = 0
    do KCLS=1,NCLS
      IAMOKAY = 1
      do IGAS=1,NGAS
        IEL = NELFGP(IOCA(IGAS,IATP))+NELFGP(IOCB(IGAS,IBTP))
        if (IEL /= ICLS(IGAS,KCLS)) IAMOKAY = 0
      end do
      if (IAMOKAY == 1) IICLS = KCLS
    end do
    ISPSPCLS(IATP,IBTP) = IICLS
  end do
end do

if (NTEST >= 10) then
  write(u6,*)
  write(u6,*) ' Matrix giving classes for alpha-beta supergroups'
  write(u6,*)
  call IWRTMA(ISPSPCLS,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
end if

end subroutine SPSPCLS_GAS
