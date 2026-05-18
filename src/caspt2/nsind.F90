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

subroutine NSIND(INS,ISYM,ICASE,IP,IQ,IR)

use Symmetry_Info, only: Mul
use definitions, only: iwp, u6
use SUPERINDEX, only: MIGEJ, MIGTJ, MAGTB, MAGEB
use caspt2_module, only: NIES, IINAIS, NIGEJES, NSES, NSYM, NIGTJES, IEXTIS, NISH, NSSH, NIGEJ, NIGTJ, NAGEBES, NAGTBES, NAGEB, &
                         NAGTB

implicit none
integer(kind=iwp), intent(in) :: INS, ISYM, ICASE
integer(kind=iwp), intent(out) :: IP, IQ, IR
integer(kind=iwp) IA, IAABS, IAB, IABABS, IAIJ, IBABS, II, IIA, IIAB, IIABS, IIJ, IIJABS, IJABS, ISYMA, ISYMAB, ISYMI, ISYMIJ, NA, &
                  NAB, NAIJ, NI, NIA, NIAB, NIJ

select case (ICASE)
  case (1)
    IIABS = INS+NIES(ISYM)
    IP = IINAIS(IIABS)
    IQ = 0
    IR = 0
  case (2,12)
    IIJABS = INS+NIGEJES(ISYM)
    IIABS = MIGEJ(1,IIJABS)
    IJABS = MIGEJ(2,IIJABS)
    IP = IINAIS(IIABS)
    IQ = IINAIS(IJABS)
    IR = 0
  case (3,13)
    IIJABS = INS+NIGTJES(ISYM)
    IIABS = MIGTJ(1,IIJABS)
    IJABS = MIGTJ(2,IIJABS)
    IP = IINAIS(IIABS)
    IQ = IINAIS(IJABS)
    IR = 0
  case (4)
    IAABS = INS+NSES(ISYM)
    IP = IEXTIS(IAABS)
    IQ = 0
    IR = 0
    return
  case (5)
    IIA = INS
    do ISYMA=1,NSYM
      ISYMI = Mul(ISYMA,ISYM)
      NIA = NISH(ISYMI)*NSSH(ISYMA)
      if (IIA <= NIA) then
        IA = 1+(IIA-1)/NISH(ISYMI)
        II = IIA-NISH(ISYMI)*(IA-1)
        IAABS = IA+NSES(ISYMA)
        IIABS = II+NIES(ISYMI)
        IP = IINAIS(IIABS)
        IQ = IEXTIS(IAABS)
        IR = 0
        return
      end if
      IIA = IIA-NIA
    end do
    write(u6,*) 'NSIND AIVX: Impossible situation.'
    call ABEND()
  case (6,7)
    IAIJ = INS
    NIJ = 0 ! dummy initialize
    do ISYMA=1,NSYM
      ISYMIJ = Mul(ISYMA,ISYM)
      if (ICASE == 6) NIJ = NIGEJ(ISYMIJ)
      if (ICASE == 7) NIJ = NIGTJ(ISYMIJ)
      NA = NSSH(ISYMA)
      NAIJ = NA*NIJ
      if (IAIJ <= NAIJ) then
        IIJ = 1+(IAIJ-1)/NA
        IA = IAIJ-NA*(IIJ-1)
        IAABS = IA+NSES(ISYMA)
        if (ICASE == 6) then
          IIJABS = IIJ+NIGEJES(ISYMIJ)
          IIABS = MIGEJ(1,IIJABS)
          IJABS = MIGEJ(2,IIJABS)
        else
          IIJABS = IIJ+NIGTJES(ISYMIJ)
          IIABS = MIGTJ(1,IIJABS)
          IJABS = MIGTJ(2,IIJABS)
        end if
        IP = IEXTIS(IAABS)
        IQ = IINAIS(IIABS)
        IR = IINAIS(IJABS)
        return
      end if
      IAIJ = IAIJ-NAIJ
    end do
    write(u6,*) 'NSIND VJAI: Impossible situation.'
    call ABEND()
  case (8)
    IABABS = INS+NAGEBES(ISYM)
    IAABS = MAGEB(1,IABABS)
    IBABS = MAGEB(2,IABABS)
    IP = IEXTIS(IAABS)
    IQ = IEXTIS(IBABS)
    IR = 0
  case (9)
    IABABS = INS+NAGTBES(ISYM)
    IAABS = MAGTB(1,IABABS)
    IBABS = MAGTB(2,IABABS)
    IP = IEXTIS(IAABS)
    IQ = IEXTIS(IBABS)
    IR = 0
  case (10,11)
    IIAB = INS
    NAB = 0 ! dummy initialize
    do ISYMI=1,NSYM
      ISYMAB = Mul(ISYMI,ISYM)
      if (ICASE == 10) NAB = NAGEB(ISYMAB)
      if (ICASE == 11) NAB = NAGTB(ISYMAB)
      NI = NISH(ISYMI)
      NIAB = NI*NAB
      if (IIAB <= NIAB) then
        IAB = 1+(IIAB-1)/NI
        II = IIAB-NI*(IAB-1)
        IIABS = II+NIES(ISYMI)
        if (ICASE == 10) then
          IABABS = IAB+NAGEBES(ISYMAB)
          IAABS = MAGEB(1,IABABS)
          IBABS = MAGEB(2,IABABS)
        else
          IABABS = IAB+NAGTBES(ISYMAB)
          IAABS = MAGTB(1,IABABS)
          IBABS = MAGTB(2,IABABS)
        end if
        IP = IINAIS(IIABS)
        IQ = IEXTIS(IAABS)
        IR = IEXTIS(IBABS)
        return
      end if
      IIAB = IIAB-NIAB
    end do
  case default
    write(u6,*) 'NSIND BJAT: Impossible situation.'
    call ABEND()
end select

end subroutine NSIND
