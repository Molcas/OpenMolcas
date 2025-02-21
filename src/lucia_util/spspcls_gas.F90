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

!#define _DEBUGPRINT_
subroutine SPSPCLS_GAS(NOCTPA,NOCTPB,IOCA,IOCB,NELFGP,MXPNGAS,NGAS,ISPSPCLS,ICLS,NCLS)
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

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NOCTPA, NOCTPB, MXPNGAS, IOCA(MXPNGAS,NOCTPA), IOCB(MXPNGAS,NOCTPB), NELFGP(*), NGAS, NCLS, &
                                 ICLS(NGAS,NCLS)
integer(kind=iwp), intent(out) :: ISPSPCLS(NOCTPA,NOCTPB)
integer(kind=iwp) :: IAMOKAY, IATP, IBTP, IEL, IGAS, IICLS, KCLS

#ifdef _DEBUGPRINT_
write(u6,*) ' ISPSPCLS_GAS entered'
write(u6,*) ' ===================='
write(u6,*)
write(u6,*) ' IOCA and IOCB'
call IWRTMA(IOCA,NGAS,NOCTPA,MXPNGAS,NGAS)
call IWRTMA(IOCB,NGAS,NOCTPB,MXPNGAS,NGAS)
write(u6,*)
write(u6,*) ' ICLS'
call IWRTMA(ICLS,NGAS,NCLS,NGAS,NCLS)
#endif

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

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' Matrix giving classes for alpha-beta supergroups'
write(u6,*)
call IWRTMA(ISPSPCLS,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
#endif

end subroutine SPSPCLS_GAS
