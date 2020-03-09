************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1998, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1998  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE COMMWEW(IVEC,JVEC,DCOM)
      USE SUPERINDEX
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "stdalloc.fh"
      DIMENSION DCOM(NASHT,NASHT)
      REAL*8, ALLOCATABLE :: CBLK(:), TBLK(:), SMAT(:)

C This subroutine is one of the components needed to compute the active/active
C transition density matrix elements for the two first-order vectors IVEC and
C JVEC.
C Adds, into the matrix DCOM, a correction obtained by commutation relations.
C Present assumption: The two vectors nr. IVEC and JVEC, stored on LUSOLV,
C are both in contravariant representation. Possibly, IVEC equals JVEC.

      CALL QENTER('COMMWEW')
      DO ICASE=1,11
        DO ISYM=1,NSYM
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          NCBLK=NAS*NIS
CTEST       WRITE(*,*)' COMMWEW ISYM,ICASE:',ISYM,ICASE
CTEST       WRITE(*,*)'                NAS:',NAS
CTEST       WRITE(*,*)'                NIS:',NIS
CTEST       WRITE(*,*)'              NCBLK:',NCBLK
          IF(NCBLK.EQ.0) GOTO 200
C Allocate CBLK, TBLK
          CALL MMA_ALLOCATE(CBLK,NCBLK)
          CALL MMA_ALLOCATE(TBLK,NCBLK)
C First, read in the ICASE, ISYM block of coefficients from JVEC into CBLK:
C Note carefully, this is not a mistake: vector JVEC into CBLK it is!
          CALL RDBLKC(ISYM,ICASE,JVEC,CBLK)
C Allocate overlap matrix:
          NS=(NAS*(NAS+1))/2
          CALL MMA_ALLOCATE(SMAT,NS)
          IDS=IDSMAT(ISYM,ICASE)
C          CALL DAFILE(LUSBT,2,SMAT,RtoI*NS,IDS)
          CALL DDAFILE(LUSBT,2,SMAT,NS,IDS)
C Compute TBLK as the covariant representation of vector JVEC, by multiplying
C with the overlap matrix. Then get rid of the overlap matrix.
          CALL DCOPY_(NCBLK,[0.0D0],0,TBLK,1)
          CALL TRIMUL(NAS,NIS,1.0D00,SMAT,CBLK,NAS,TBLK,NAS)
          CALL MMA_DEALLOCATE(SMAT)
C Finally, if IVEC not equals JVEC, read in the contravariant block of vector
C IVEC into CBLK:
          IF(IVEC.NE.JVEC) THEN
            CALL RDBLKC(ISYM,ICASE,IVEC,CBLK)
          END IF
C Finally, branch to the appropriate code section:
          GOTO(1,2,3,4,5,6,7,8,9,10,11) ICASE

C Case 1 code section:
   1  CONTINUE
      K000=NTUVES(ISYM)

      DO ISYMX=1,NSYM
        NAX=NASH(ISYMX)
        DO IX=1,NAX
          IXABS=NAES(ISYMX)+IX
          DO IY=1,NAX
            IYABS=NAES(ISYMX)+IY

            SUM=0.0D0
            ISYMTU=MUL(ISYMX,ISYM)
            DO ITABS=1,NASHT
              ISYMT=IASYM(ITABS)
              ISYMU=MUL(ISYMT,ISYMTU)
              DO IU=1,NASH(ISYMU)
                IUABS=NAES(ISYMU)+IU

                IYTU=KTUV(IYABS,ITABS,IUABS)-K000
                IXTU=KTUV(IXABS,ITABS,IUABS)-K000
                ITYU=KTUV(ITABS,IYABS,IUABS)-K000
                ITXU=KTUV(ITABS,IXABS,IUABS)-K000
                ITUY=KTUV(ITABS,IUABS,IYABS)-K000
                ITUX=KTUV(ITABS,IUABS,IXABS)-K000

                DO IIS=1,NIS
                  SUM=SUM+CBLK(IXTU+NAS*(IIS-1))
     &                                *TBLK(IYTU+NAS*(IIS-1))
                  SUM=SUM+CBLK(ITXU+NAS*(IIS-1))
     &                                *TBLK(ITYU+NAS*(IIS-1))
                  SUM=SUM-CBLK(ITUY+NAS*(IIS-1))
     &                                *TBLK(ITUX+NAS*(IIS-1))
                END DO

              END DO
            END DO
            DCOM(IXABS,IYABS)=DCOM(IXABS,IYABS)+SUM

          END DO
        END DO
      END DO
      GOTO 100

C Case 2 code section:
   2  CONTINUE
      DO ISYMX=1,NSYM
        NAX=NASH(ISYMX)
        DO IX=1,NAX
          IXABS=NAES(ISYMX)+IX
          DO IY=1,NAX
            IYABS=NAES(ISYMX)+IY

            SUM=0.0D0
            ISYMT=MUL(ISYMX,ISYM)
            DO IT=1,NASH(ISYMT)
              ITABS=NAES(ISYMT)+IT
              IF(ITABS.GE.IXABS) THEN
                IXT=KTGEU(ITABS,IXABS)-NTGEUES(ISYM)
              ELSE
                IXT=KTGEU(IXABS,ITABS)-NTGEUES(ISYM)
              END IF
              IF(ITABS.GE.IYABS) THEN
                IYT=KTGEU(ITABS,IYABS)-NTGEUES(ISYM)
              ELSE
                IYT=KTGEU(IYABS,ITABS)-NTGEUES(ISYM)
              END IF
              PARTSUM=0.0d0
              DO IIS=1,NIS
                PARTSUM=PARTSUM+CBLK(IXT+NAS*(IIS-1))
     &                               *TBLK(IYT+NAS*(IIS-1))
              END DO
              IF(ITABS.EQ.IXABS) PARTSUM=2.0D0*PARTSUM
              SUM=SUM+PARTSUM
            END DO
            DCOM(IXABS,IYABS)=DCOM(IXABS,IYABS)+SUM

          END DO
        END DO
      END DO
      GOTO 100

C Case 3 code section:
   3  CONTINUE
      DO ISYMX=1,NSYM
        NAX=NASH(ISYMX)
        DO IX=1,NAX
          IXABS=NAES(ISYMX)+IX
          DO IY=1,NAX
            IYABS=NAES(ISYMX)+IY

            SUM=0.0D0
            ISYMT=MUL(ISYMX,ISYM)
            DO IT=1,NASH(ISYMT)
              ITABS=NAES(ISYMT)+IT
              IF(ITABS.EQ.IXABS) GOTO 390
              IF(ITABS.EQ.IYABS) GOTO 390
              IF(ITABS.GT.IXABS) THEN
                IXT=KTGTU(ITABS,IXABS)-NTGTUES(ISYM)
                SGN=1.0d0
              ELSE
                IXT=KTGTU(IXABS,ITABS)-NTGTUES(ISYM)
                SGN=-1.0d0
              END IF
              IF(ITABS.GT.IYABS) THEN
                IYT=KTGTU(ITABS,IYABS)-NTGTUES(ISYM)
              ELSE
                IYT=KTGTU(IYABS,ITABS)-NTGTUES(ISYM)
                SGN=-SGN
              END IF
              PARTSUM=0.0d0
              DO IIS=1,NIS
                PARTSUM=PARTSUM+CBLK(IXT+NAS*(IIS-1))
     &                                *TBLK(IYT+NAS*(IIS-1))
              END DO
              SUM=SUM+SGN*PARTSUM

 390          CONTINUE
            END DO
            DCOM(IXABS,IYABS)=DCOM(IXABS,IYABS)+SUM

          END DO
        END DO
      END DO
      GOTO 100

C Case 4 code section:
   4  CONTINUE
      K000=NTUVES(ISYM)

      DO ISYMX=1,NSYM
        NAX=NASH(ISYMX)
        DO IX=1,NAX
          IXABS=NAES(ISYMX)+IX
          DO IY=1,NAX
            IYABS=NAES(ISYMX)+IY

            SUM=0.0D0
            ISYMTU=MUL(ISYMX,ISYM)
            DO ITABS=1,NASHT
              ISYMT=IASYM(ITABS)
              ISYMU=MUL(ISYMT,ISYMTU)
              DO IU=1,NASH(ISYMU)
                IUABS=NAES(ISYMU)+IU

                IYTU=KTUV(IYABS,ITABS,IUABS)-K000
                IXTU=KTUV(IXABS,ITABS,IUABS)-K000
                ITYU=KTUV(ITABS,IYABS,IUABS)-K000
                ITXU=KTUV(ITABS,IXABS,IUABS)-K000
                ITUY=KTUV(ITABS,IUABS,IYABS)-K000
                ITUX=KTUV(ITABS,IUABS,IXABS)-K000

                DO IIS=1,NIS
                  SUM=SUM-CBLK(IYTU+NAS*(IIS-1))
     &                         *TBLK(IXTU+NAS*(IIS-1))
                  SUM=SUM+CBLK(ITXU+NAS*(IIS-1))
     &                         *TBLK(ITYU+NAS*(IIS-1))
                  SUM=SUM-CBLK(ITUY+NAS*(IIS-1))
     &                         *TBLK(ITUX+NAS*(IIS-1))
                END DO

              END DO
            END DO
            DCOM(IXABS,IYABS)=DCOM(IXABS,IYABS)+SUM

          END DO
        END DO
      END DO
      GOTO 100

C Case 5 code section:
   5  CONTINUE
      NAS1=NAS/2

      DO ISYMX=1,NSYM
        NAX=NASH(ISYMX)
        DO IX=1,NAX
          IXABS=NAES(ISYMX)+IX
          DO IY=1,NAX
            IYABS=NAES(ISYMX)+IY

            SUM=0.0D0
            ISYMT=MUL(ISYMX,ISYM)
            DO IT=1,NASH(ISYMT)
              ITABS=NAES(ISYMT)+IT
              IXT1=KTU(IXABS,ITABS)-NTUES(ISYM)
              IYT1=KTU(IYABS,ITABS)-NTUES(ISYM)
              ITX1=KTU(ITABS,IXABS)-NTUES(ISYM)
              ITY1=KTU(ITABS,IYABS)-NTUES(ISYM)
              IXT2=IXT1+NAS1
              IYT2=IYT1+NAS1
              ITX2=ITX1+NAS1
              ITY2=ITY1+NAS1
              DO IIS=1,NIS
                SUM=SUM+CBLK(IXT1+NAS*(IIS-1))
     &                      *TBLK(IYT1+NAS*(IIS-1))
                SUM=SUM-CBLK(ITY1+NAS*(IIS-1))
     &                      *TBLK(ITX1+NAS*(IIS-1))
                SUM=SUM+CBLK(IXT2+NAS*(IIS-1))
     &                      *TBLK(IYT2+NAS*(IIS-1))
                SUM=SUM-CBLK(ITY2+NAS*(IIS-1))
     &                      *TBLK(ITX2+NAS*(IIS-1))
              END DO
            END DO
            DCOM(IXABS,IYABS)=DCOM(IXABS,IYABS)+SUM

          END DO
        END DO
      END DO
      GOTO 100

C Case 6 code section:
   6  CONTINUE
      NAX=NASH(ISYM)
      DO IX=1,NAX
        IXABS=NAES(ISYM)+IX
        DO IY=1,NAX
          IYABS=NAES(ISYM)+IY

          SUM=0.0D0
          DO IIS=1,NIS
            SUM=SUM+CBLK(IX+NAS*(IIS-1))
     &                      *TBLK(IY+NAS*(IIS-1))
          END DO
          DCOM(IXABS,IYABS)=DCOM(IXABS,IYABS)+SUM

        END DO
      END DO
      GOTO 100

C Case 7 code section:
   7  CONTINUE
      NAX=NASH(ISYM)
      DO IX=1,NAX
        IXABS=NAES(ISYM)+IX
        DO IY=1,NAX
          IYABS=NAES(ISYM)+IY

          SUM=0.0D0
          DO IIS=1,NIS
            SUM=SUM+CBLK(IX+NAS*(IIS-1))
     &                   *TBLK(IY+NAS*(IIS-1))
          END DO
          DCOM(IXABS,IYABS)=DCOM(IXABS,IYABS)+SUM

        END DO
      END DO
      GOTO 100

C Case 8 code section:
   8  CONTINUE
      DO ISYMX=1,NSYM
        NAX=NASH(ISYMX)
        DO IX=1,NAX
          IXABS=NAES(ISYMX)+IX
          DO IY=1,NAX
            IYABS=NAES(ISYMX)+IY

            SUM=0.0D0
            ISYMT=MUL(ISYMX,ISYM)
            DO IT=1,NASH(ISYMT)
              ITABS=NAES(ISYMT)+IT
              IF(ITABS.GE.IXABS) THEN
                IXT=KTGEU(ITABS,IXABS)-NTGEUES(ISYM)
              ELSE
                IXT=KTGEU(IXABS,ITABS)-NTGEUES(ISYM)
              END IF
              IF(ITABS.GE.IYABS) THEN
                IYT=KTGEU(ITABS,IYABS)-NTGEUES(ISYM)
              ELSE
                IYT=KTGEU(IYABS,ITABS)-NTGEUES(ISYM)
              END IF
              PARTSUM=0.0d0
              DO IIS=1,NIS
                PARTSUM=PARTSUM-CBLK(IYT+NAS*(IIS-1))
     &                       *TBLK(IXT+NAS*(IIS-1))
              END DO
              IF(ITABS.EQ.IYABS) PARTSUM=2.0D0*PARTSUM
              SUM=SUM+PARTSUM

            END DO
            DCOM(IXABS,IYABS)=DCOM(IXABS,IYABS)+SUM

          END DO
        END DO
      END DO
      GOTO 100

C Case 9 code section:
   9  CONTINUE
      DO ISYMX=1,NSYM
        NAX=NASH(ISYMX)
        DO IX=1,NAX
          IXABS=NAES(ISYMX)+IX
          DO IY=1,NAX
            IYABS=NAES(ISYMX)+IY

            SUM=0.0D0
            ISYMT=MUL(ISYMX,ISYM)
            DO IT=1,NASH(ISYMT)
              ITABS=NAES(ISYMT)+IT
              IF(ITABS.EQ.IXABS) GOTO 990
              IF(ITABS.EQ.IYABS) GOTO 990
              IF(ITABS.GT.IXABS) THEN
                IXT=KTGTU(ITABS,IXABS)-NTGTUES(ISYM)
                SGN=1.0d0
              ELSE
                IXT=KTGTU(IXABS,ITABS)-NTGTUES(ISYM)
                SGN=-1.0d0
              END IF
              IF(ITABS.GT.IYABS) THEN
                IYT=KTGTU(ITABS,IYABS)-NTGTUES(ISYM)
              ELSE
                IYT=KTGTU(IYABS,ITABS)-NTGTUES(ISYM)
                SGN=-SGN
              END IF
              PARTSUM=0.0d0
              DO IIS=1,NIS
                PARTSUM=PARTSUM-CBLK(IYT+NAS*(IIS-1))
     &                               *TBLK(IXT+NAS*(IIS-1))
              END DO
              SUM=SUM+SGN*PARTSUM

 990          CONTINUE
            END DO
            DCOM(IXABS,IYABS)=DCOM(IXABS,IYABS)+SUM

          END DO
        END DO
      END DO
      GOTO 100

C Case 10 code section:
  10  CONTINUE
      NAX=NASH(ISYM)
      DO IX=1,NAX
        IXABS=NAES(ISYM)+IX
        DO IY=1,NAX
          IYABS=NAES(ISYM)+IY

          SUM=0.0D0
          DO IIS=1,NIS
            SUM=SUM-CBLK(IY+NAS*(IIS-1))
     &                   *TBLK(IX+NAS*(IIS-1))
          END DO
          DCOM(IXABS,IYABS)=DCOM(IXABS,IYABS)+SUM

        END DO
      END DO
      GOTO 100

C Case 11 code section:
  11  CONTINUE
      NAX=NASH(ISYM)
      DO IX=1,NAX
        IXABS=NAES(ISYM)+IX
        DO IY=1,NAX
          IYABS=NAES(ISYM)+IY

          SUM=0.0D0
          DO IIS=1,NIS
            SUM=SUM-CBLK(IY+NAS*(IIS-1))
     &                   *TBLK(IX+NAS*(IIS-1))
          END DO
          DCOM(IXABS,IYABS)=DCOM(IXABS,IYABS)+SUM

        END DO
      END DO
      GOTO 100


 100  CONTINUE
      CALL MMA_DEALLOCATE(CBLK)
      CALL MMA_DEALLOCATE(TBLK)

 200  CONTINUE

C Here ends the loops over ISYM and ICASE.
        END DO
      END DO

      CALL QEXIT('COMMWEW')
      RETURN
      END
