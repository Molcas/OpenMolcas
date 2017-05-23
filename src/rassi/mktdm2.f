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
      SUBROUTINE MKTDM2(LSYM1,MSPROJ1,IFSBTAB1,
     &                  LSYM2,MSPROJ2,IFSBTAB2,ISSTAB,
     &                  MAPORB,DET1,DET2,NTDM2,TDM2)
      IMPLICIT NONE
      INTEGER IFSBTAB1(*),IFSBTAB2(*),ISSTAB(*),MAPORB(*),NTDM2
      REAL*8 DET1(*),DET2(*)
      REAL*8 TDM2(NTDM2)
      INTEGER IORB,ITABS,IUABS,JORB,LORBTB
      INTEGER NASHT,NASORB
      REAL*8 SGNJL,SGNIK
      REAL*8 GVAL,GAAAA,GABBA,GBAAB,GBBBB,GABAB,GBABA
      INTEGER LSYM1,MSPROJ1,LSYM2,MSPROJ2,ISYOP,MS2OP
      INTEGER IAAAA,IABAB,IABBA,IAKA,IAKB,IBAAB,IBABA,IBBBB,IBIA
      INTEGER IBKA,IBKB,IJ,IJIJ,IORBA,IORBB,ITU,ITUVX
      INTEGER IVABS,IVX,IXABS,JALA,JALB,JBJA,JBLA,JBLB
      INTEGER JORBA,JORBB,KORB,KORBA,KORBB,LORB,LORBA,LORBB
      INTEGER LSPD2,NASGEM,NSPD2
#include "symmul.fh"
#include "WrkSpc.fh"

C Given two CI expansions, using a biorthonormal set of SD''s,
C calculate the spin-summed 2-particle transition density matrix
C in the biorthonormal active orbital basis.

      LORBTB=ISSTAB(3)
C Pick out nr of active orbitals from orbital table:
      NASORB=IWORK(LORBTB+3)
      NASHT=NASORB/2
      NASGEM=(NASORB*(NASORB-1))/2
      NSPD2=NASGEM**2
      CALL GETMEM('SPD2','Allo','Real',LSPD2,NSPD2)
      CALL DCOPY_(NSPD2,0.0D0,0,WORK(LSPD2),1)
      ISYOP=MUL(LSYM1,LSYM2)
      MS2OP=MSPROJ1-MSPROJ2
      CALL SPIND2(ISYOP,MS2OP,IWORK(LORBTB),ISSTAB,
     &            IFSBTAB1,IFSBTAB2,
     &            DET1,DET2,WORK(LSPD2))
      SGNJL=1.0D0 ! dummy initialize
      SGNIK=1.0D0 ! dummy initialize
      IAKA=0      ! dummy initialize
      IAKB=0      ! dummy initialize
      IBIA=0      ! dummy initialize
      IBKA=0      ! dummy initialize
      IBKB=0      ! dummy initialize
      JALA=0      ! dummy initialize
      JALB=0      ! dummy initialize
      JBLA=0      ! dummy initialize
      JBLB=0      ! dummy initialize
      JBJA=0      ! dummy initialize
      DO JORB=1,NASHT
       JORBA=2*JORB-1
       JORBB=2*JORB
       IUABS=MAPORB(JORBA)
       DO IORB=1,NASHT
        IORBA=2*IORB-1
        IORBB=2*IORB
        ITABS=MAPORB(IORBA)
        ITU=ITABS+NASHT*(IUABS-1)
        DO LORB=1,NASHT
         LORBA=2*LORB-1
         LORBB=2*LORB
         IXABS=MAPORB(LORBA)
      IF(JORB.GT.LORB) THEN
        SGNJL=1.0D0
        JALA=((JORBA-1)*(JORBA-2))/2+LORBA
        JALB=((JORBA-1)*(JORBA-2))/2+LORBB
        JBLA=((JORBB-1)*(JORBB-2))/2+LORBA
        JBLB=((JORBB-1)*(JORBB-2))/2+LORBB
      ELSE IF(JORB.EQ.LORB) THEN
        JBJA=((JORBB-1)*(JORBB-2))/2+JORBA
      ELSE
        SGNJL=-1.0D0
        JALA=((LORBA-1)*(LORBA-2))/2+JORBA
        JALB=((LORBB-1)*(LORBB-2))/2+JORBA
        JBLA=((LORBA-1)*(LORBA-2))/2+JORBB
        JBLB=((LORBB-1)*(LORBB-2))/2+JORBB
      END IF
         DO KORB=1,NASHT
          KORBA=2*KORB-1
          KORBB=2*KORB
          IVABS=MAPORB(KORBA)
          IVX=IVABS+NASHT*(IXABS-1)
          IF(ITU.LT.IVX) GOTO 123
      IF(IORB.GT.KORB) THEN
        SGNIK=1.0D0
        IAKA=((IORBA-1)*(IORBA-2))/2+KORBA
        IAKB=((IORBA-1)*(IORBA-2))/2+KORBB
        IBKA=((IORBB-1)*(IORBB-2))/2+KORBA
        IBKB=((IORBB-1)*(IORBB-2))/2+KORBB
      ELSE IF(IORB.EQ.KORB) THEN
        IBIA=((IORBB-1)*(IORBB-2))/2+IORBA
      ELSE
        SGNIK=-1.0D0
        IAKA=((KORBA-1)*(KORBA-2))/2+IORBA
        IAKB=((KORBB-1)*(KORBB-2))/2+IORBA
        IBKA=((KORBA-1)*(KORBA-2))/2+IORBB
        IBKB=((KORBB-1)*(KORBB-2))/2+IORBB
      END IF
      IF(IORB.NE.KORB) THEN
       IF(JORB.NE.LORB) THEN
        IAAAA=IAKA+NASGEM*(JALA-1)
        IABBA=IAKB+NASGEM*(JALB-1)
        IBAAB=IBKA+NASGEM*(JBLA-1)
        IBBBB=IBKB+NASGEM*(JBLB-1)
        GAAAA=WORK(LSPD2-1+IAAAA)
        GABBA=WORK(LSPD2-1+IABBA)
        GBAAB=WORK(LSPD2-1+IBAAB)
        GBBBB=WORK(LSPD2-1+IBBBB)
        GVAL=SGNIK*SGNJL*(GAAAA+GABBA+GBAAB+GBBBB)
       ELSE
        IABAB=IAKB+NASGEM*(JBJA-1)
        IBAAB=IBKA+NASGEM*(JBJA-1)
        GABAB=WORK(LSPD2-1+IABAB)
        GBAAB=WORK(LSPD2-1+IBAAB)
        GVAL=SGNIK*(-GABAB+GBAAB)
       END IF
      ELSE
       IF(JORB.NE.LORB) THEN
        IBABA=IBIA+NASGEM*(JALB-1)
        IBAAB=IBIA+NASGEM*(JBLA-1)
        GBABA=WORK(LSPD2-1+IBABA)
        GBAAB=WORK(LSPD2-1+IBAAB)
        GVAL=SGNJL*(-GBABA+GBAAB)
       ELSE
        IBAAB=IBIA+NASGEM*(JBJA-1)
        GBAAB=WORK(LSPD2-1+IBAAB)
        GVAL=2.0D0*GBAAB
       END IF
      END IF
C Position determined by active orbital index in external order:
          ITUVX=(ITU*(ITU-1))/2+IVX
          TDM2(ITUVX)=GVAL
 123      CONTINUE
         END DO
        END DO
       END DO
      END DO
      CALL GETMEM('SPD2','Free','Real',LSPD2,NSPD2)
C DIAGONAL ELEMENTS HALF-SIZED (This is for proper contraction with TUVX):
      IJIJ=0
      DO IJ=1,NASHT**2
        IJIJ=IJIJ+IJ
        TDM2(IJIJ)=0.5D0*TDM2(IJIJ)
      END DO
      RETURN
      END
