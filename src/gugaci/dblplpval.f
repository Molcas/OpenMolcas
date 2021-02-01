************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
C dbl space partial loop values
      SUBROUTINE Value_of_PL_IN_DBL()
#include "drt_h.fh"
      DATA DZERO,DONE,DHALF/0.D0,1.D0,0.5D0/
      DATA DSQ2,VSQ2/1.414213562373095D0, 0.7071067811865D0 /
        DB=JB_SYS
        IB=mod(JB_SYS,2)
        FG=DONE
        IF(IB.EQ.1) FG=-DONE
!------------------- diag ---------------------------------------
!  W0=+(-)FG*DSQ2  or  W0=+(-)FG/DSQ2
!W1=-(B_SYS+3)/(2*B_SYS+2)
      W1_CS21=-(DB+3)/(2*DB+2)    ! W1  for  D&R&L(2)D^R^L(1)
      W1_CS12=-(DB-1)/(2*DB+2)   ! W1  for  D&R&L(1)D^R^L(2)
      W1_CT11=0
      W1_CT22=0
      PD = FG*sqrt((DB+3)/(2*DB+2))
      PT =-FG*sqrt((DB+4)/(2*DB+4))
      PS1= FG*sqrt(DB/(2*DB+4))                   ! (2)D&R&L(1)
      PS4=-FG*(DB+3)*sqrt((DB)/(2*DB+4))/(DB+1)      !D&R&L(2)C"(1)

      IF(JB_SYS.EQ.0) GOTO 100
      PDD=-FG*sqrt((DB-1)/(2*DB+2))
      PS2=-FG*sqrt((DB+2)/(2*DB))                ! (1)D&R&L(2)
      PS3= FG*(DB-1)*sqrt((DB+2)/(2*DB))/(DB+1)      !D&R&L(1)C"(2)

      IF(JB_SYS.EQ.1) GOTO 100
      PTT= FG*sqrt((DB-2)/(2*DB))

!------------------- off diag ---------------------------------------
100     FGSQ2 =FG*DSQ2
        FGVSQ2=FG*VSQ2

!SS(1-2)  Ar(02)-Bl(31)-
!SS(1-4)  Ar(23)-Bl(10)-
!SS(1-5)  (22)-Ar(13)-Bl(31)-
!SS(1-9)  Ar(23)-C'(11)-Bl(32)-
!SS(1-10) Ar(23)-C'(12)-Bl(31)-
!SS(1-14) Ar(23)-Bl(32)-C"(11)-
!SS(1-15) (22)-Drl(11)-
!SS(1-17) Drl(22)-C"(11)-
!SS(1-20) Drl(33)-C"(00)-
!SS(1-20) Drl(33)-C"(22)-C"(11)-
!SS(1-20) (22)Drl(33)-C"(11)-
!SS(1-20) (22)(11)Drl(33)-
        W0_SS(2) = FG*sqrt((DB+2)/(2*DB+2))
        W1_SS(2) =-FG*sqrt(DB/(2*DB+2))
        W0_SS(4) = FG*sqrt((DB+2)/(2*DB+2))
        W1_SS(4) = FG*sqrt(DB/(2*DB+2))
        W0_SS(5) = FGVSQ2
        W1_SS(5) =-FG*sqrt(DB/(2*DB+4))
        W0_SS(9) =-FGVSQ2*sqrt(DB*(DB+2))/(DB+1)
        W1_SS(9) =-FGVSQ2*(DB+2)/(DB+1)
        W0_SS(10) =-FGVSQ2/(DB+1)
        W1_SS(10) = FG*sqrt(DB/(2*DB+4))/(DB+1)
        W0_SS(14) = FGVSQ2
        W1_SS(14) = FG*(DB+3)*sqrt(DB/(2*DB+4))/(DB+1)
        W0_SS(15) = FGVSQ2
        W1_SS(15) = FG*sqrt(DB/(2*DB+4))
        W0_SS(17) = FGVSQ2
        W1_SS(17) =-FG*(DB+3)*sqrt(DB/(2*DB+4))/(DB+1)      !??????
        W0_SS(20) = FGSQ2
        W1_SS(20) = DZERO

!ST(2-1) Ar(02)-Bl(32)-
!ST(2-2) (22)Ar(13)-Bl(32)-
!ST(2-4) Ar(23)-C'(12)-Bl(32)-
!ST(2-4) Ar(23)-Bl(32)-C'(12)-
!ST(2-5) (22)Drl(12)-
!ST(2-6) Drl(22)-C"(12)-
        W1_ST(1)= DONE
        W1_ST(2)= sqrt((DB+1)/(DB+2))
        W1_ST(4)=-1/sqrt((DB+1)*(DB+2))
        W1_ST(5)=-sqrt((DB+1)/(DB+2))
        W1_ST(6)=1/sqrt((DB+1)*(DB+2))

!TS(3-1) Ar(23)-Bl(20)-
!TS(3-2) (22)Ar(23)-Bl(31)-
!TS(3-2) Ar(23)-C'(22)-Bl(31)-
!TS(3-4) Ar(23)-Bl(32)-C"(21)-
        W1_TS(1)= sqrt((DB+3)/(DB+1))
        W1_TS(2)=-sqrt((DB+3)/(DB+2))
        W1_TS(4)= sqrt((DB+3)/(DB+2))/(DB+1)

!SD(6-1) A&r(02)-
!SD(6-2) C(22)A&(13)-
!SD(6-4) A&r(23)C'(12)-
!SD(6-5) A&r(23)B&r(13)B^r(32)
!SD(6-8) A&r(23)B&l(32)B^l(13)
!SD(6-9) D&r&r(03)B^r(32)
!SD(6-11) D&r&l(22)B^l(13)
!SD(6-12) D&r&l(33)B^l(02)
!SD(6-14) (22)D&r&l(33)B^l(13)
!SD(6-14) D&r&l(33)C"(22)B^l(13)
!SD(6-16) D&r&l(33)B^l(23)C'(12)
        W0_SD(1)= DONE
        W1_SD(1)= DONE
        W0_SD(2)= sqrt((DB+1)/(DB+2))
        W1_SD(2)= W0_SD(2)
        W0_SD(4)=-1/sqrt((DB+1)*(DB+2))
        W1_SD(4)= W0_SD(4)
        W0_SD(5)= sqrt((DB+2)/(DB+1))/2
        W1_SD(5)= DB/(2*sqrt((DB+1)*(DB+2)))
        W0_SD(8)=-sqrt((DB+1)/(DB+2))/2
        W1_SD(8)=-(DB+3)/(2*sqrt((DB+1)*(DB+2)))
        W0_SD(9)=-DONE
        W1_SD(9)= DZERO
        W0_SD(11)=-sqrt((DB+1)/(DB+2))/2
        W1_SD(11)= (DB+3)/(2*sqrt((DB+1)*(DB+2)))
        W0_SD(12)=-DONE
        W1_SD(12)= DZERO
        W0_SD(14)=-sqrt((DB+1)/(DB+2))
        W1_SD(14)= DZERO
        W0_SD(16)= 1/sqrt((DB+1)*(DB+2))
        W1_SD(16)= DZERO

!DS(7-1) Ar(23)-DRl(30)-
!DS(7-3) Ar(23)-Bl(32)-BR(31)-
        W0_DS(1)=-sqrt((DB+2)/(DB+1))
        W1_DS(1)= DZERO
        W0_DS(3)= DHALF
        W1_DS(3)= (DB+3)/(2*DB+2)

!SV(10-2) Ar(23)-Br(13)-
!SV(10-3) Drr(03)-
        W0_SV(2)= sqrt((DB+2)/(2*DB+2))
        W1_SV(2)=-sqrt(DB/(2*DB+2))
        W0_SV(3)=-DSQ2
        W1_SV(3)= DZERO

!TT(11-1) (22)Ar(23)-Bl(32)-
!TT(11-1) Ar(23)-C'(22)-Bl(32)-
!TT(11-1) Ar(23)-Bl(32)-C"(22)-
!TT(11-2) (22)Drl(22)-
!TT(11-2) Drl(22)-C"(22)-
!TT(11-3) (22)Drl(33)-
!TT(11-3) Drl(33)-C"(22)-
!TT(11-3) Drl(33)-C"(22)-C"(22)-
        W0_TT(1)= FGVSQ2
        W1_TT(1)= FG*sqrt((DB+4)/(2*DB+4))
        W0_TT(2)= FGVSQ2
        W1_TT(2)=-FG*sqrt((DB+4)/(2*DB+4))
        W0_TT(3)= FGSQ2
        W1_TT(3)= DZERO

!TD(13-1) (22)A&(23)
!TD(13-1) A&(23)C'(22)
!TD(13-2) A&(23)B&r(23)B^r(32)
!TD(13-3) A&(23)B&l(32)B^l(23)
!TD(13-4) D&r&l(22)B^l(23)
!TD(13-5) (22)D&&l(33)B^l(23)
!TD(13-5) D&rl(33)C"(22)B^l(23)
!TD(13-5) D&rl(33)B^l(23)C'(22)
        W0_TD(1)= FG*sqrt((DB+3)/(DB+2))
        W0_TD(2)= DZERO
        W1_TD(2)= W0_TD(1)
        W0_TD(3)=-W0_TD(1)/2
        W1_TD(3)= W0_TD(1)/2
        W0_TD(4)= W0_TD(3)
        W1_TD(4)= W0_TD(3)
        W0_TD(5)=-W0_TD(1)
        W1_TD(5)= DZERO

!DT(14) Ar(23)-Bl(32)-BR(32)-
        W0_DT=-FG/2
        W1_DT= FG/2

!TV(17) Ar(23)-Br(23)-
        W1_TV= DZERO
        W1_TV=-FG*sqrt((DB+3)/(DB+1))

!DD(19-1) Ar(23)-Bl(32)-
!DD(19-2) Drl(22)-
!DD(19-3) (22)Drl(33)-
!DD(19-3) Drl(33)-C"(22)-
        W0_DD(1)=-FGVSQ2
        W1_DD(1)=-FG*sqrt((DB+3)/(2*DB+2))
        W0_DD(2)=-FGVSQ2
        W1_DD(2)= FG*sqrt((DB+3)/(2*DB+2))
        W0_DD(3)=-FGSQ2
        W1_DD(3)= DZERO

!DV(23-1) Ar(23)-
!DV(23-2) Drl(33)-BL(23)-
        W0_DV(1)=-FG*sqrt((DB+2)/(DB+1))
        W1_DV(1)= W0_DV(1)
        W0_DV(2)= FG*sqrt((DB+2)/(DB+1))
        W1_DV(2)= DZERO

!VV(25) Drl(33)-
        W0_VV   = FGSQ2

C===============================================================
        IF(JB_SYS.EQ.0) RETURN

!SS(1-1)  Ar(01)-Bl(32)-
!SS(1-3)  Ar(13)-Bl(20)-
!SS(1-6)  (11)-Ar(23)-Bl(32)-
!SS(1-7)  Ar(13)-C'(21)-Bl(32)-
!SS(1-8)  Ar(13)-C'(22)-Bl(31)-
!SS(1-11) Ar(13)-Bl(31)-C"(22)-
!SS(1-12) Ar(13)-Bl(32)-C"(21)-
!SS(1-13) Ar(23)-Bl(31)-C"(12)-
!SS(1-16) (11)-Drl(22)-
!SS(1-18) Drl(11)-C"(22)-
!SS(1-19) Drl(12)-C"(21)-
!SS(1-20) Drl(33)-C"(11)-C"(22)-
!SS(1-20) (11)Drl(33)-C"(22)-
!SS(1-20) (11)(22)Drl(33)-
        W0_SS(1) = FGVSQ2*sqrt(DB/(DB+1))
        W1_SS(1) = FGVSQ2*sqrt((DB+2)/(DB+1))
        W0_SS(3) = FGVSQ2*sqrt(DB/(DB+1))
        W1_SS(3) =-FGVSQ2*sqrt((DB+2)/(DB+1))
        W0_SS(6) = FGVSQ2
        W1_SS(6) = FGVSQ2*sqrt((DB+2)/DB)
        W0_SS(7) = FGVSQ2/(DB+1)
        W1_SS(7) = FGVSQ2*sqrt((DB+2)/DB)/(DB+1)
        W0_SS(8) =-FGVSQ2*sqrt(DB*(DB+2))/(DB+1)
        W1_SS(8) = FGVSQ2*DB/(DB+1)
        W0_SS(11) = FGVSQ2
        W1_SS(11) =-FG*sqrt((DB+2)/(2*DB))*(DB-1)/(DB+1)
        W0_SS(12) = DZERO
        W1_SS(12) =-FGSQ2/(DB+1)
        W0_SS(13) = DZERO
        W1_SS(13) =-FGSQ2/(DB+1)
        W0_SS(16) = FGVSQ2
        W1_SS(16) =-FGVSQ2*sqrt((DB+2)/DB)
        W0_SS(18) = FGVSQ2
        W1_SS(18) = FGVSQ2*(DB-1)*sqrt((DB+2)/DB)/(DB+1)
        W0_SS(19) = DZERO
        W1_SS(19) = FGSQ2/(DB+1)

!ST(2-3) Ar(13)-C'(22)-Bl(32)-
!ST(2-3) Ar(13)-Bl(32)-C'(22)-
!ST(2-7) Drl(12)-C"(22)-
        W1_ST(3)=-sqrt(DB/(DB+1))
        W1_ST(7)= sqrt(DB/(DB+1))

!TS(3-3) Ar(23)-Bl(31)-C"(22)-
        W1_TS(3)= sqrt(DB*(DB+3))/(DB+1)

!SD(6-3) A&r(13)C'(22)-
!SD(6-6) A&r(13)B&r(23)B^r(32)
!SD(6-7) A&r(13)B&l(32)B^l(23)
!SD(6-10) D&r&l(12)B^l(23)
!SD(6-15) D&r&l(33)B^l(13)C'(22)
        W0_SD(3)=-sqrt(DB/(DB+1))
        W1_SD(3)= W0_SD(3)
        W0_SD(6)= sqrt(DB/(DB+1))/2
        W1_SD(6)=-W0_SD(6)
        W0_SD(7)= DZERO
        W1_SD(7)=-sqrt(DB/(DB+1))
        W0_SD(10)= DZERO
        W1_SD(10)= sqrt(DB/(DB+1))
        W0_SD(15)= sqrt(DB/(DB+1))
        W1_SD(15)= DZERO

!DS(7-2) Ar(23)-Bl(31)-BR(32)-
        W0_DS(2)= DZERO
        W1_DS(2)= sqrt(DB*(DB+2))/(DB+1)

!SD1(8-1)    Ar(01)-
!SD1(8-2)    Ar(23)-
!SD1(8-3)    Ar(13)-C'(21)-
!SD1(8-4)    Ar(23)-C'(11)-
!SD1(8-5)    Ar(13)-Br(23)-BR(31)-
!SD1(8-6)    Ar(23)-Br(13)-BR(31)-
!SD1(8-7)    Ar(13)-Bl(31)-BL(23)-
!SD1(8-8)    Ar(23)-Bl(31)-BL(13)-
!SD1(8-9)    Drr(03)-BR(31)-
!SD1(8-9)    Drl(33)-BL(01)-
!SD1(8-10) Drl(11)-BL(23)-
!SD1(8-11) (11)Drl(33)-BL(23)-
!SD1(8-11) Drl(33)-C"(11)-BL(23)-
!SD1(8-12) Drl(33)-BL(13)-C'(21)-
!SD1(8-13) Drl(33)-BL(23)-C'(11)-
        W0_SD1(1) = FG
        W1_SD1(1) = FG
        W0_SD1(2) = FG*sqrt((DB+1)/DB)
        W1_SD1(2) = W0_SD1(2)
        W0_SD1(3) = FG/sqrt(DB*(DB+1))
        W1_SD1(3) = W0_SD1(3)
        W0_SD1(4) =-FG*sqrt((DB+2)/(DB+1))
        W1_SD1(4) = W0_SD1(4)
        W0_SD1(5) = FG*sqrt(DB/(DB+1))/2
        W1_SD1(5) = FG*(DB+2)/(2*sqrt(DB*(DB+1)))
        W0_SD1(6) = FG*sqrt((DB+2)/(DB+1))/2
        W1_SD1(6) =-W0_SD1(6)
        W0_SD1(7) =-FG*sqrt((DB+1)/DB)/2
        W1_SD1(7) =-FG*(DB-1)/(2*sqrt(DB*(DB+1)))
        W0_SD1(8) = DZERO
        W1_SD1(8) =-FG*sqrt((DB+2)/(DB+1))
        W0_SD1(9) =-FG
        W1_SD1(9) = DZERO
        W0_SD1(10) =-FG*sqrt((DB+1)/DB)/2
        W1_SD1(10) = FG*(DB-1)/(2*sqrt(DB*(DB+1)))
        W0_SD1(11) =-FG*sqrt((DB+1)/DB)
        W1_SD1(11) = DZERO
        W0_SD1(12) =-FG/sqrt(DB*(DB+1))
        W1_SD1(12) = DZERO
        W0_SD1(13) = FG*sqrt((DB+2)/(DB+1))
        W1_SD1(13) = DZERO

!D1S(9-1)    Ar(13)-DlR(30)-
!D1S(9-2)    Ar(13)-Bl(31)-BR(32)-
!D1S(9-3)    Ar(13)-Bl(32)-BR(31)-
!D1S(9-4)    Drl(12)-Br(31)-
        W0_D1S(1) = FG*sqrt(DB/(DB+1))
        W1_D1S(1) = DZERO
        W0_D1S(2) =-FG/2
        W1_D1S(2) =-FG*(DB-1.D0)/(2*DB+2)
        W0_D1S(3) = DZERO
        W1_D1S(3) =-FG*sqrt(DB*(DB+2))/(DB+1)
        W0_D1S(4) = DZERO
        W1_D1S(4) = FG*sqrt(DB*(DB+2))/(DB+1)

!SV(10-1) Ar(13)-Br(23)-
        W0_SV(1)= sqrt(DB/(2*DB+2))
        W1_SV(1)= sqrt((DB+2)/(2*DB+2))

!D1D1(20-1) Ar(13)-BL(31)-
!D1D1(20-1) Drl(11)-
!D1D1(20-1) Drl(33)-
!D1D1(20-1) Drl(33)-C"(11)-
        W0_D1D1(1) =-FGVSQ2
        W1_D1D1(1) = FGVSQ2*sqrt((DB-1)/(DB+1))
        W0_D1D1(2) = W0_D1D1(1)
        W1_D1D1(2) =-W1_D1D1(1)
        W0_D1D1(3) =-FGSQ2
        W1_D1D1(3) = DZERO

!DD1(21)Ar(23)-Bl(31)-
        W1_DD1     =-sqrt((DB+2)/(DB+1))

!D1D(22-1)   Ar(13)-Bl(32)-
!D1D(22-2)   Drl(12)-
        W0_D1D(1) = DZERO
        W1_D1D(1) = sqrt(DB/(DB+1))
        W0_D1D(2) = DZERO
        W1_D1D(2) =-W1_D1D(1)
!D1V(24-1)  Ar(13)-
!D1V(24-2)  Drl(33)-BL(13)-

        W0_D1V(1) = sqrt(DB/(DB+1))
        W1_D1V(1) = W0_D1V(1)
        W0_D1V(2) =-W0_D1V(1)
        W1_D1V(2) = DZERO

C=======================================================================
        IF(JB_SYS.EQ.1) RETURN

!ST1(4-1) Ar(01)-Bl(31)-
!ST1(4-2) Ar(23)-Bl(31)-
!ST1(4-3) Ar(13)-C'(21)-Bl(31)-
!ST1(4-3) Ar(13)-Bl(31)-C"(21)-
!ST1(4-4) Ar(23)-C'(11)-Bl(31)-
!ST1(4-4) Ar(23)-Bl(31)-C"(11)-
        W1_ST1(1) =-1.d0
        W1_ST1(2) =-sqrt((DB+1)/DB)
        W1_ST1(3) =-1/sqrt(DB*(DB+1))
        W1_ST1(4) = sqrt((DB+2)/(DB+1))

!T1S(5-1)   Ar(13)-Bl(10)-
!T1S(5-2)   Ar(13)-Bl(32)-
!T1S(5-2)   Ar(13)-C'(11)-Bl(32)-
!T1S(5-3)   Ar(13)-Bl(31)-C"(12)-
!T1S(5-4)   Ar(13)-Bl(32)-C"(11)-
!T1S(5-5)   Drl(12)-
!T1S(5-6)   Drl(12)-C"(12)-
!T1S(5-7)   Drl(12)-C"(11)-
        W1_T1S(1) =-sqrt((DB-1)/(DB+1))
        W1_T1S(2) = sqrt((DB-1)/DB)
        W1_T1S(3) = sqrt((DB-1)/DB)/(DB+1)
        W1_T1S(4) =-sqrt((DB-1)*(DB+2))/(DB+1)
        W1_T1S(5) =-sqrt((DB-1)/DB)
        W1_T1S(6) =-sqrt((DB-1)/DB)/(DB+1)
        W1_T1S(7) = sqrt((DB-1)*(DB+2))/(DB+1)

!T1T1(12-1)  Ar(13)-Bl(31)-
!T1T1(12-1)  Ar(13)-C'(11)-Bl(31)-
!T1T1(12-1)  Ar(13)-Bl(31)-C"(11)-
!T1T1(12-2)  Drl(11)-
!T1T1(12-2)  Drl(11)-C"(11)-
!T1T1(12-3)  Drl(33)-
!T1T1(12-3)  Drl(33)-C"(11)-
!T1T1(12-3)  Drl(33)-C"(11)-C"(11)-
        W0_T1T1(1) = FGVSQ2
        W1_T1T1(1) =-FGVSQ2*sqrt((DB-2)/DB)
        W0_T1T1(2) = FGVSQ2
        W1_T1T1(2) =-W1_T1T1(1)
        W0_T1T1(3) = FGSQ2
        W1_T1T1(3) = DZERO

!T1D1(15-1)  Ar(13)-
!T1D1(15-1)  Ar(13)-C'(11)-
!T1D1(15-2)  Ar(13)-Br(13)-BR(31)-
!T1D1(15-3)  Ar(13)-Bl(31)-Bl(13)-
!T1D1(15-4)  Drl(11)-BL(13)-
!T1D1(15-5)  Drl(33)-BL(13)-
!T1D1(15-5)  Drl(33)-C"(11)-BL(13)-
!T1D1(15-5)  Drl(33)-BL(13)-C'(11)-
        W0_T1D1(1) = sqrt((DB-1)/DB)
        W1_T1D1(1) = W0_T1D1(1)
        W0_T1D1(2) = DZERO
        W1_T1D1(2) = W0_T1D1(1)
        W0_T1D1(3) =-W0_T1D1(1)/2
        W1_T1D1(3) = W0_T1D1(1)/2
        W0_T1D1(4) =-W0_T1D1(1)/2
        W1_T1D1(4) = W0_T1D1(4)
        W0_T1D1(5) =-W0_T1D1(1)
        W1_T1D1(5) = DZERO

!D1T1(16) Ar(13)-Bl(31)-BR(31)-
        W0_D1T1=0.5D0
        W1_D1T1=-0.5D0

!T1V(18) Ar(13)-Br(13)-
        W0_T1V= DZERO
        W1_T1V=   FG*sqrt((DB-1)/(DB+1))

      END

      SUBROUTINE SS2_EXT(LRI,LRJ,NK)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      LMI=LSM_INN(LRI)
      LMJ=LSM_INN(LRJ)
      LMIJ=MUL_TAB(LMI,LMJ)
      NK=0
      IF(JML.NE.1.OR.LMIJ.NE.JMR) RETURN
!SS(1-2)  Ar(02)-Bl(31)-
      NK=1
      W0SS2=W0_SS(2)
      W1SS2=W1_SS(2)
      NI=MOD(LRJ-LRI,2)
      IF(NI.EQ.0) THEN
        W0SS2=-W0SS2
        W1SS2=-W1SS2
      ENDIF
      DO MPL=1,MTYPE
        VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS2
        VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS2
      ENDDO
      IWDL=JUST(LRI,LRI)
      IWDR=JUST(LRI,LRJ)
      DO MPL=1,MHLP
         IWAL=LPNEW_LWEI(MPL)
         IWAR=LPNEW_RWEI(MPL)
         LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
         LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      ENDDO
      RETURN
      END

      SUBROUTINE SS3_EXT(LRI,LRJ)
!SS(1-3)  Ar(13)-Bl(20)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      W0SS3=W0_SS(3)
      W1SS3=W1_SS(3)
      NI=MOD(LRJ-LRI,2)
      IF(NI.EQ.0) THEN
        W0SS3=-W0SS3
        W1SS3=-W1SS3
      ENDIF
      DO MPL=1,MTYPE
        VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS3
        VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS3
      ENDDO
      IWDL=JUST(LRJ,LRI)
      IWDR=JUST(LRJ,LRJ)
      DO MPL=1,MHLP
         IWAL=LPNEW_LWEI(MPL)
         IWAR=LPNEW_RWEI(MPL)
         LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
         LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      ENDDO
      RETURN
      END

      SUBROUTINE SS4_EXT(LRI,LRJ,NK)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      LMI=LSM_INN(LRI)
      LMJ=LSM_INN(LRJ)
      LMIJ=MUL_TAB(LMI,LMJ)
      NK=0
      IF(JMR.NE.1.OR.LMIJ.NE.JML) RETURN
      NK=1
!SS(1-4)  Ar(23)-Bl(10)-
      W0SS4=W0_SS(4)
      W1SS4=W1_SS(4)
      NI=MOD(LRJ-LRI,2)
      IF(NI.EQ.0) THEN
        W0SS4=-W0SS4
        W1SS4=-W1SS4
      ENDIF
      DO MPL=1,MTYPE
        VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS4
        VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS4
      ENDDO
      IWDL=JUST(LRI,LRJ)
      IWDR=JUST(LRJ,LRJ)
      DO MPL=1,MHLP
         IWAL=LPNEW_LWEI(MPL)
         IWAR=LPNEW_RWEI(MPL)
         LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
         LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      ENDDO
      RETURN
      END

      SUBROUTINE SS5_EXT(LRI,LRJ,NK)
!SS(1-5)  (22)-Ar(13)-Bl(31)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      DIMENSION IWDL(MAX_INNORB),IWDR(MAX_INNORB)
      NK=0
      LMI=LSM_INN(LRI)
      LMJ=LSM_INN(LRJ)
      DO LRK=NORB_FRZ+1,LRI-1
        LMK=LSM_INN(LRK)
        LMKI=MUL_TAB(LMK,LMI)
        LMKJ=MUL_TAB(LMK,LMJ)
        IF(LMKI.NE.JML.OR.LMKJ.NE.JMR) CYCLE
        NK=NK+1
        IWDL(NK)=JUST(LRK,LRI)
        IWDR(NK)=JUST(LRK,LRJ)
      ENDDO
      IF(NK.EQ.0) RETURN

      W0SS5=W0_SS(5)
      W1SS5=W1_SS(5)
      NI=MOD(LRJ-LRI,2)
      IF(NI.EQ.0) THEN
        W0SS5=-W0SS5
        W1SS5=-W1SS5
      ENDIF
      DO MPL=1,MTYPE
        VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS5
        VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS5
      ENDDO
      NPL=0
      DO MPL=1,MHLP
        IWAL=LPNEW_LWEI(MPL)
        IWAR=LPNEW_RWEI(MPL)
        DO K=1,NK
          NPL=NPL+1
         LP_LWEI(NPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL(K))
          LP_RWEI(NPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR(K))
        ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE SS10_EXT(LRI,LRJ,NK)
!SS(1-10) Ar(23)-C'(12)-Bl(31)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      DIMENSION IWDL(MAX_INNORB),IWDR(MAX_INNORB)
      NK=0
      LMI=LSM_INN(LRI)
      LMJ=LSM_INN(LRJ)
      DO LRK=LRI+1,LRJ-1
        LMK=LSM_INN(LRK)
        LMKI=MUL_TAB(LMK,LMI)
        LMKJ=MUL_TAB(LMK,LMJ)
        IF(LMKI.NE.JML.OR.LMKJ.NE.JMR) CYCLE
        NK=NK+1
        IWDL(NK)=JUST(LRI,LRK)
        IWDR(NK)=JUST(LRK,LRJ)
      ENDDO
      IF(NK.EQ.0) RETURN

      W0SS10=-W0_SS(10)
      W1SS10=-W1_SS(10)
      NI=MOD(LRJ-LRI,2)
      IF(NI.EQ.0) THEN
        W0SS10=-W0SS10
        W1SS10=-W1SS10
      ENDIF
      DO MPL=1,MTYPE
        VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS10
        VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS10
      ENDDO
      NPL=0
      DO MPL=1,MHLP
        IWAL=LPNEW_LWEI(MPL)
        IWAR=LPNEW_RWEI(MPL)
        DO K=1,NK
          NPL=NPL+1
          LP_LWEI(NPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL(K))
          LP_RWEI(NPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR(K))
        ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE SS14_EXT(LRI,LRJ,NK)
!SS(1-14) Ar(23)-Bl(32)-C"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      DIMENSION IWDL(MAX_INNORB),IWDR(MAX_INNORB)
      NK=0
      LMI=LSM_INN(LRI)
      LMJ=LSM_INN(LRJ)
      DO LRK=LRJ+1,NORB_DZ
        LMK=LSM_INN(LRK)
        LMKI=MUL_TAB(LMK,LMI)
        LMKJ=MUL_TAB(LMK,LMJ)
        IF(LMKI.NE.JML.OR.LMKJ.NE.JMR) CYCLE
        NK=NK+1
        IWDL(NK)=JUST(LRI,LRK)
        IWDR(NK)=JUST(LRJ,LRK)
      ENDDO
      IF(NK.EQ.0) RETURN

      W0SS14=W0_SS(14)
      W1SS14=W1_SS(14)
      NI=MOD(LRJ-LRI,2)
      IF(NI.EQ.0) THEN
        W0SS14=-W0SS14
        W1SS14=-W1SS14
      ENDIF
      DO MPL=1,MTYPE
        VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0SS14
        VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1SS14
      ENDDO
      NPL=0
      DO MPL=1,MHLP
        IWAL=LPNEW_LWEI(MPL)
        IWAR=LPNEW_RWEI(MPL)
        DO K=1,NK
          NPL=NPL+1
          LP_LWEI(NPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL(K))
          LP_RWEI(NPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR(K))
        ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE TT1_EXT(LRI,LRJ,NK,IGF)
!SS(1-14) Ar(23)-Bl(32)-C"(11)-
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      DIMENSION IWDL(MAX_INNORB),IWDR(MAX_INNORB)
      NK=0
      LMI=LSM_INN(LRI)
      LMJ=LSM_INN(LRJ)
      IF(IGF.EQ.-1) GOTO 200
!TT(11-1) (22)Ar(23)-Bl(32)-
      DO LRK=NORB_FRZ+1,LRI-1
        LMK=LSM_INN(LRK)
        LMKI=MUL_TAB(LMK,LMI)
        LMKJ=MUL_TAB(LMK,LMJ)
        IF(LMKI.NE.JML.OR.LMKJ.NE.JMR) CYCLE
        NK=NK+1
        IWDL(NK)=JUST(LRK,LRI)
        IWDR(NK)=JUST(LRK,LRJ)
      ENDDO
!TT(11-1) Ar(23)-Bl(32)-C"(22)-    ACT -C"-
      DO LRK=LRJ+1,NORB_DZ
        LMK=LSM_INN(LRK)
        LMKI=MUL_TAB(LMK,LMI)
        LMKJ=MUL_TAB(LMK,LMJ)
        IF(LMKI.NE.JML.OR.LMKJ.NE.JMR) CYCLE
        NK=NK+1
        IWDL(NK)=JUST(LRI,LRK)     !
        IWDR(NK)=JUST(LRJ,LRK)     !
      ENDDO
      W0TT1=W0_TT(1)
      W1TT1=W1_TT(1)
      GOTO 300
!TT(11-1) Ar(23)-C'(22)-Bl(32)-
200   DO LRK=LRI+1,LRJ-1
        LMK=LSM_INN(LRK)
        LMKI=MUL_TAB(LMK,LMI)
        LMKJ=MUL_TAB(LMK,LMJ)
        IF(LMKI.NE.JML.OR.LMKJ.NE.JMR) CYCLE
        NK=NK+1
        IWDL(NK)=JUST(LRI,LRK)     !
        IWDR(NK)=JUST(LRK,LRJ)     !
      ENDDO
      W0TT1=-W0_TT(1)
      W1TT1=-W1_TT(1)

300   IF(NK.EQ.0) RETURN

      NI=MOD(LRJ-LRI,2)
      IF(NI.EQ.0) THEN
        W0TT1=-W0TT1
        W1TT1=-W1TT1
      ENDIF
      DO MPL=1,MTYPE
        VPLP_W0(MPL)=VPLPNEW_W0(MPL)*W0TT1
        VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TT1
      ENDDO
      NPL=0
      DO MPL=1,MHLP
        IWAL=LPNEW_LWEI(MPL)
        IWAR=LPNEW_RWEI(MPL)
        DO K=1,NK
          NPL=NPL+1
         LP_LWEI(NPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL(K))
          LP_RWEI(NPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR(K))
        ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE TS1_EXT(LRI,LRJ,NK)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      LMI=LSM_INN(LRI)
      LMJ=LSM_INN(LRJ)
      LMIJ=MUL_TAB(LMI,LMJ)
      NK=0
      IF(JMR.NE.1.OR.LMIJ.NE.JML) RETURN
      NK=1
      W1TS1=W1_TS(1)
      NI=MOD(LRJ-LRI,2)
      IF(NI.EQ.0) W1TS1=-W1TS1
!TS(3-1) Ar(23)-Bl(20)-
      DO MPL=1,MTYPE
        VPLP_W0(MPL)=0.D0
        VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TS1
      ENDDO
      IWDL=JUST(LRI,LRJ)      !
      IWDR=JUST(LRJ,LRJ)
      DO MPL=1,MHLP
        IWAL=LPNEW_LWEI(MPL)
        IWAR=LPNEW_RWEI(MPL)
        LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      ENDDO
      RETURN
      END

      SUBROUTINE TS2_EXT(LRI,LRJ,NK,IGF)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      DIMENSION IWDL(MAX_INNORB),IWDR(MAX_INNORB)
      NK=0
      LMI=LSM_INN(LRI)
      LMJ=LSM_INN(LRJ)
      IF(IGF.EQ.-1) GOTO 200
!TS(3-2) (22)Ar(23)-Bl(31)-
      DO LRK=NORB_FRZ+1,LRI-1
        LMK=LSM_INN(LRK)
        LMKI=MUL_TAB(LMK,LMI)
        LMKJ=MUL_TAB(LMK,LMJ)
        IF(LMKI.NE.JML.OR.LMKJ.NE.JMR) CYCLE
        NK=NK+1
        IWDL(NK)=JUST(LRK,LRI)      !
        IWDR(NK)=JUST(LRK,LRJ)
      ENDDO
      W1TS2=W1_TS(2)
      GOTO 300
!TS(3-2) Ar(23)-C'(22)-Bl(31)-   ACT -C"-
200   DO LRK=LRI+1,LRJ-1
        LMK=LSM_INN(LRK)
        LMKI=MUL_TAB(LMK,LMI)
        LMKJ=MUL_TAB(LMK,LMJ)
        IF(LMKI.NE.JML.OR.LMKJ.NE.JMR) CYCLE
        NK=NK+1
        IWDL(NK)=JUST(LRI,LRK)        !
        IWDR(NK)=JUST(LRK,LRJ)
      ENDDO
      W1TS2=-W1_TS(2)
300   IF(NK.EQ.0) RETURN
      NI=MOD(LRJ-LRI,2)
      IF(NI.EQ.0) THEN
        W1TS2=-W1TS2
      ENDIF
      DO MPL=1,MTYPE
        VPLP_W0(MPL)=0.D0
        VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TS2
      ENDDO
      NPL=0
      DO MPL=1,MHLP
        IWAL=LPNEW_LWEI(MPL)
        IWAR=LPNEW_RWEI(MPL)
        DO K=1,NK
          NPL=NPL+1
         LP_LWEI(NPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL(K))
          LP_RWEI(NPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR(K))
        ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE TS4_EXT(LRI,LRJ,NK)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      DIMENSION IWDL(MAX_INNORB),IWDR(MAX_INNORB)
      LMI=LSM_INN(LRI)
      LMJ=LSM_INN(LRJ)
      NK=0
!TS(3-4) Ar(23)-Bl(32)-C"(21)-
      DO LRK=LRJ+1,NORB_DZ
        LMK=LSM_INN(LRK)
        LMKI=MUL_TAB(LMK,LMI)
        LMKJ=MUL_TAB(LMK,LMJ)
        IF(LMKI.NE.JML.OR.LMKJ.NE.JMR) CYCLE
        NK=NK+1
        IWDL(NK)=JUST(LRI,LRK)          !
        IWDR(NK)=JUST(LRJ,LRK)
      ENDDO
      IF(NK.EQ.0) RETURN
      W1TS4=W1_TS(4)
      NI=MOD(LRJ-LRI,2)
      IF(NI.EQ.0) THEN
        W1TS4=-W1TS4
      ENDIF
      DO MPL=1,MTYPE
        VPLP_W0(MPL)=0.D0
        VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1TS4
      ENDDO
      NPL=0
      DO MPL=1,MHLP
        IWAL=LPNEW_LWEI(MPL)
        IWAR=LPNEW_RWEI(MPL)
        DO K=1,NK
          NPL=NPL+1
         LP_LWEI(NPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL(K))
          LP_RWEI(NPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR(K))
        ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE ST1_EXT(LRI,LRJ,NK)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      LMI=LSM_INN(LRI)
      LMJ=LSM_INN(LRJ)
      LMIJ=MUL_TAB(LMI,LMJ)
      NK=0
      IF(JML.NE.1.OR.LMIJ.NE.JMR) RETURN
      NK=1
      W1ST1=W1_ST(1)
      NI=MOD(LRJ-LRI,2)
      IF(NI.EQ.0) W1ST1=-W1ST1
!ST(2-1) Ar(02)-Bl(32)-
      IWDL=JUST(LRI,LRI)
      IWDR=JUST(LRI,LRJ)     !
      DO MPL=1,MHLP
        IWAL=LPNEW_LWEI(MPL)
        IWAR=LPNEW_RWEI(MPL)
        LP_LWEI(MPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL)
        LP_RWEI(MPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR)
      ENDDO
      DO MPL=1,MTYPE
        VPLP_W0(MPL)=0.D0
        VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1ST1
      ENDDO
      RETURN
      END

      SUBROUTINE ST2_EXT(LRI,LRJ,NK)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      DIMENSION IWDL(MAX_INNORB),IWDR(MAX_INNORB)
      LMI=LSM_INN(LRI)
      LMJ=LSM_INN(LRJ)
      NK=0
!ST(2-2) (22)Ar(13)-Bl(32)-
      DO LRK=NORB_FRZ+1,LRI-1
        LMK=LSM_INN(LRK)
        LMKI=MUL_TAB(LMK,LMI)
        LMKJ=MUL_TAB(LMK,LMJ)
        IF(LMKI.NE.JML.OR.LMKJ.NE.JMR) CYCLE
        NK=NK+1
        IWDL(NK)=JUST(LRK,LRI)
        IWDR(NK)=JUST(LRK,LRJ)     !
      ENDDO
      IF(NK.EQ.0) RETURN

      W1ST2=W1_ST(2)
      NI=MOD(LRJ-LRI,2)
      IF(NI.EQ.0) THEN
        W1ST2=-W1ST2
      ENDIF
      DO MPL=1,MTYPE
        VPLP_W0(MPL)=0.D0
        VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1ST2
      ENDDO
      NPL=0
      DO MPL=1,MHLP
        IWAL=LPNEW_LWEI(MPL)
        IWAR=LPNEW_RWEI(MPL)
        DO K=1,NK
          NPL=NPL+1
         LP_LWEI(NPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL(K))
          LP_RWEI(NPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR(K))
        ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE ST4_EXT(LRI,LRJ,NK,IGF)
#include "drt_h.fh"
#include "intsort_h.fh"
#include "pl_structure_h.fh"
#include "lpextmode_h.fh"
      DIMENSION IWDL(MAX_INNORB),IWDR(MAX_INNORB)
      NK=0
      LMI=LSM_INN(LRI)
      LMJ=LSM_INN(LRJ)
      IF(IGF.EQ.-1) GOTO 200
!ST(2-4) Ar(23)-Bl(32)-C'(12)-
      DO LRK=LRJ+1,NORB_DZ
        LMK=LSM_INN(LRK)
        LMKI=MUL_TAB(LMK,LMI)
        LMKJ=MUL_TAB(LMK,LMJ)
        IF(LMKI.NE.JML.OR.LMKJ.NE.JMR) CYCLE
        NK=NK+1
        IWDL(NK)=JUST(LRI,LRK)      !
        IWDR(NK)=JUST(LRJ,LRK)
      ENDDO
      W1ST4=W1_ST(4)
      GOTO 300

!ST(2-4) Ar(23)-C'(12)-Bl(32)-
200   DO LRK=LRI+1,LRJ-1
        LMK=LSM_INN(LRK)
        LMKI=MUL_TAB(LMK,LMI)
        LMKJ=MUL_TAB(LMK,LMJ)
        IF(LMKI.NE.JML.OR.LMKJ.NE.JMR) CYCLE
        NK=NK+1
        IWDL(NK)=JUST(LRI,LRK)        !
        IWDR(NK)=JUST(LRK,LRJ)
      ENDDO
      W1ST4=-W1_ST(4)

300   IF(NK.EQ.0) RETURN
      NI=MOD(LRJ-LRI,2)
      IF(NI.EQ.0) THEN
        W1ST4=-W1ST4
      ENDIF
      DO MPL=1,MTYPE
        VPLP_W0(MPL)=0.D0
        VPLP_W1(MPL)=VPLPNEW_W1(MPL)*W1ST4
      ENDDO
      NPL=0
      DO MPL=1,MHLP
        IWAL=LPNEW_LWEI(MPL)
        IWAR=LPNEW_RWEI(MPL)
        DO K=1,NK
          NPL=NPL+1
         LP_LWEI(NPL)=IWALK_AD(JPADL,IpaeL,IWAL,IWDL(K))
          LP_RWEI(NPL)=IWALK_AD(JPAD,Ipae,IWAR,IWDR(K))
        ENDDO
      ENDDO
      RETURN
      END
