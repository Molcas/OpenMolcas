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
subroutine IAIBCM_GAS(LCMBSPC,ICMBSPC,MNMXOC,NOCTPA,NOCTPB,IOCA,IOCB,NELFTP,MXPNGAS,NGAS,IOCOC)
! Allowed combinations of alpha and beta types, GAS version
!
! =====
! Input
! =====
!
! LCMBSPC : Number of GAS spaces included in this expnasion
! ICMBSPC : Gas spaces included in this expansion
!
! MXMNOC(IGAS,1,IGASSPC) : Min accumulated occ for AS 1-IGAS for space IGASSPC
! MXMNOC(IGAS,2,IGASSPC) : Max accumulated occ for AS 1-IGAS for space IGASSPC
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
! IOCOC(IATP,IBTP) == 1 =>     allowed combination
! IOCOC(IATP,IBTP) == 0 => not allowed combination

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: LCMBSPC, ICMBSPC(LCMBSPC), MXPNGAS, MNMXOC(MXPNGAS,2,*), NOCTPA, NOCTPB, IOCA(MXPNGAS,NOCTPA), &
                                 IOCB(MXPNGAS,NOCTPB), NELFTP(*), NGAS
integer(kind=iwp), intent(out) :: IOCOC(NOCTPA,NOCTPB)
integer(kind=iwp) :: IAMOKAY, IATP, IBTP, IEL, IGAS, INC, JCMBSPC, JJCMBSPC
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: II

write(u6,*) ' IAICBM_GAS entered'
write(u6,*) ' ==================='
write(u6,*)
write(u6,*) ' Number of GAS spaces included ',LCMBSPC
write(u6,*) ' GAS spaces included ',(ICMBSPC(II),II=1,LCMBSPC)
write(u6,*)
write(u6,*) ' IOCA and IOCB'
call IWRTMA(IOCA,NGAS,NOCTPA,MXPNGAS,NGAS)
call IWRTMA(IOCB,NGAS,NOCTPB,MXPNGAS,NGAS)
#endif

IOCOC(:,:) = 0
do IATP=1,NOCTPA
  do IBTP=1,NOCTPB
    ! is this combination allowed in any of the GAS spaces included
    INC = 0
    do JJCMBSPC=1,LCMBSPC
      JCMBSPC = ICMBSPC(JJCMBSPC)
      IEL = 0
      IAMOKAY = 1
      do IGAS=1,NGAS
        IEL = IEL+NELFTP(IOCA(IGAS,IATP))+NELFTP(IOCB(IGAS,IBTP))
        if ((IEL < MNMXOC(IGAS,1,JCMBSPC)) .or. (IEL > MNMXOC(IGAS,2,JCMBSPC))) IAMOKAY = 0
      end do
      if (IAMOKAY == 1) INC = 1
    end do

    !if (I_RE_MS2_SPACE /= 0) then
    !  ! Spin projection after space I_RE_MS2_SPACE :
    !  MS2_INTERM = 0
    !  do IGAS=1,I_RE_MS2_SPACE
    !    MS2_INTERM = MS2_INTERM+NELFTP(IOCA(IGAS,IATP))-NELFTP(IOCB(IGAS,IBTP))
    !  end do
    !  if (MS2_INTERM /= I_RE_MS2_VALUE) INC = 0
    !end if

    ! Congratulations, you are allowed
    if (INC == 1) IOCOC(IATP,IBTP) = 1
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' Matrix giving allowed combinations of types'
write(u6,*)
call IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
#endif

end subroutine IAIBCM_GAS
