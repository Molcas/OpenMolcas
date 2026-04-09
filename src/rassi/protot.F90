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
! Copyright (C) 1999, Per Ake Malmqvist                                *
!***********************************************************************

subroutine PROTOT(NPORB,NPSDSZ,IPSDMS,NPCSFSZ,IPCSFCP,PCSFTOSD)
! Expand csf's in terms of determinants by the Grabenstetter method
!  (I.J.Q.C. 10, P142 (1976))
! Recoded by PAM 1999, after Jeppe Olsen.
!
! Input :
!         NPORB    : NUMBER OF OPEN ORBITALS
!         IPSDMS  : OCCUPATION OF PROTO-SDs
!         NPSDSZ   : NUMBER OF PROTOSD's
!         NPCSFSZ  : NUMBER OF PROTOCSF's
!         IPCSFCP  : SPIN COUPLINGS IN PROTO-CSFs
! Output :
!         PCSFTOSD :  NPSDSZ X NPCSFSZ MATRIX
!                GIVING EXPANSION FROM P-SD'S TO P-CSF'S

use rassi_aux, only: ipglob
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NPORB, NPSDSZ, IPSDMS(NPORB,NPSDSZ), NPCSFSZ, IPCSFCP(NPORB,NPCSFSZ)
real(kind=wp), intent(inout) :: PCSFTOSD(NPSDSZ,NPCSFSZ)
integer(kind=iwp) :: IC, IM, INDSMM, INDSPM, IOPEN, JCSF, JDET
real(kind=wp) :: COEF1, COEF2
integer(kind=iwp), parameter :: ASPIN = 1, BSPIN = 0, DWNCPL = 0, UPCPL = 1

do JCSF=1,NPCSFSZ
  if (IPGLOB >= 5) write(u6,*) ' ....Output for P-CSF ',JCSF
  do JDET=1,NPSDSZ
    ! EXPANSION COEFFICIENT OF DETERMINANT JDET FOR P-CSF JCSF
    COEF1 = One
    COEF2 = One
    INDSMM = 0
    INDSPM = 0
    do IOPEN=1,NPORB
      IC = 0
      if (IPCSFCP(IOPEN,JCSF) == UPCPL) IC = 1
      IM = 0
      if (IPSDMS(IOPEN,JDET) == ASPIN) IM = 1
      if (IC == 0) then
        if (IM == 0) then
          INDSPM = INDSPM-1
          COEF1 = COEF1*sqrt(real(INDSPM+1,kind=wp))
          ! If COEF1 has gone down to 0 exactly.
          if (INDSPM+1 == 0) exit
        else
          INDSMM = INDSMM-1
          COEF1 = -COEF1*sqrt(real(INDSMM+1,kind=wp))
          ! If COEF1 has gone down to 0 exactly.
          if (INDSMM+1 == 0) exit
        end if
        COEF2 = COEF2*sqrt(real(INDSPM+INDSMM+2,kind=wp))
      else
        if (IM == 0) then
          INDSMM = INDSMM+1
          COEF1 = COEF1*sqrt(real(INDSMM,kind=wp))
        else
          INDSPM = INDSPM+1
          COEF1 = COEF1*sqrt(real(INDSPM,kind=wp))
        end if
        COEF2 = COEF2*sqrt(real(INDSPM+INDSMM,kind=wp))
      end if
    end do

    PCSFTOSD(JDET,JCSF) = COEF1/COEF2

  end do
end do

end subroutine PROTOT
