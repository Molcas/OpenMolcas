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
      SUBROUTINE CNF2TXT(IFORM,NORB,NCLS,NOPN,ICONF,LENGTH,TEXT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ICONF(*)
      CHARACTER*(*) TEXT
      CHARACTER*1 SEP

      MXWR=LEN(TEXT)
      NWR=1
      TEXT(1:1)='('
      IWORD = 0 ! dummy initialize
      IF(IFORM.EQ.1 .OR. IFORM.EQ.3) THEN
        NOCC=NCLS+NOPN
        IF(NCLS.EQ.0) THEN
          NWR=2
          TEXT(2:2)=';'
        END IF
        IF(IFORM.EQ.1) THEN
          DO IOCC=1,NOCC
            IORB=ICONF(IOCC)
            SEP=','
            IF(IOCC.EQ.NCLS) SEP=';'
            IF(IORB.LT.10) THEN
              NWR=MIN(MXWR,NWR+2)
              WRITE(TEXT(NWR-1:NWR),'(I1,A1)')IORB,SEP
            ELSE IF(IORB.LT.100) THEN
              NWR=MIN(MXWR,NWR+3)
              WRITE(TEXT(NWR-2:NWR),'(I2,A1)')IORB,SEP
            ELSE
              NWR=MIN(MXWR,NWR+4)
              WRITE(TEXT(NWR-3:NWR),'(I3,A1)')IORB,SEP
            END IF
          END DO
        ELSE
          DO IOCC=1,NOCC
            IW=(3+IOCC)/4
            IREST=(3+IOCC)-4*IW
            IF(IREST.EQ.0) THEN
              IWORD=ICONF(IW)
            END IF
            IORB=MOD(IWORD,256)
            IWORD=IWORD/256
            SEP=','
            IF(IOCC.EQ.NCLS) SEP=';'
            IF(IORB.LT.10) THEN
              NWR=MIN(MXWR,NWR+2)
              WRITE(TEXT(NWR-1:NWR),'(I1,A1)')IORB,SEP
            ELSE IF(IORB.LT.100) THEN
              NWR=MIN(MXWR,NWR+3)
              WRITE(TEXT(NWR-2:NWR),'(I2,A1)')IORB,SEP
            ELSE
              NWR=MIN(MXWR,NWR+4)
              WRITE(TEXT(NWR-3:NWR),'(I3,A1)')IORB,SEP
            END IF
          END DO
        END IF
      ELSE IF(IFORM.EQ.2 .OR. IFORM.EQ.4) THEN
        IF(IFORM.EQ.2) THEN
          DO IORB=1,NORB
            IOCC=ICONF(IORB)
            NWR=MIN(MXWR,NWR+1)
            WRITE(TEXT(NWR:NWR),'(I1)')IOCC
          END DO
        ELSE
          DO IORB=1,NORB
            IW=(IORB+14)/15
            IREST=IORB+14-15*IW
            IF(IREST.EQ.0) THEN
              IWORD=ICONF(IW)
            END IF
            IOCC=MOD(IWORD,4)
            IWORD=IWORD/4
            NWR=MIN(MXWR,NWR+1)
            WRITE(TEXT(NWR:NWR),'(I1)')IOCC
          END DO
        END IF
      END IF
      TEXT(NWR:NWR)=')'
      LENGTH=NWR
      RETURN
      END
