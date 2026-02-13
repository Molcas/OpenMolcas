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
      SUBROUTINE PRSCTAB(SCTAB,TRANS)
      use definitions, only: iwp, wp, u6
      IMPLICIT NONE
      INTEGER(kind=iwp), intent(in):: SCTAB(*)
      REAL(kind=wp), intent(in):: TRANS(*)

      INTEGER(kind=iwp) NSIZE,ITYPE,MLTPL,MS2,MINOP,MAXOP,LTRANS,NTRANS,
     &                  N,IOPEN,NCP,NBLK,IBLK,NCPL,ND,KSPCPL,KSPDET
      INTEGER(kind=iwp), External:: ngene

      WRITE(u6,*)
      WRITE(u6,*)'------------------------------------------'
      WRITE(u6,*)' Spin Coupling Table printout'
      WRITE(u6,*)'------------------------------------------'
      NSIZE=SCTAB(1)
      ITYPE=SCTAB(2)
      MLTPL=SCTAB(3)
      MS2  =SCTAB(4)
      MINOP=SCTAB(5)
      MAXOP=SCTAB(6)
      LTRANS=SCTAB(7)
      NTRANS=SCTAB(8)
      WRITE(u6,'(1x,A,I16)')' Table size       :',NSIZE
      WRITE(u6,'(1x,A,I16)')' Table type ID    :',ITYPE
      WRITE(u6,'(1x,A,I16)')' Spin multiplicity:',MLTPL
      WRITE(u6,'(1x,A,I16)')' Spin projection  :',MS2
      WRITE(u6,'(1x,A,I16)')' Open shells; min :',MINOP
      WRITE(u6,'(1x,A,I16)')' Open shells; max :',MAXOP
      WRITE(u6,'(1x,A,I16)')' Transf data; addr:',LTRANS
      WRITE(u6,'(1x,A,I16)')' Transf data; wrds:',NTRANS
C Number of (non-trivial) values of IOPEN:
      N=0
      DO IOPEN=MINOP,MAXOP
        NCP=NGENE(IOPEN,MLTPL)
        IF(NCP.GT.0) N=N+1
      END DO
      IF(N.EQ.0) THEN
        WRITE(u6,*)
        WRITE(u6,*)' There is no such spin-coupling scheme.'
        WRITE(u6,*)
      ELSE
       WRITE(u6,'(1x,A,I9)')'   Nr of schemes  :',N
       NBLK=MAXOP-MINOP+1
       DO IBLK=1,NBLK
        IOPEN=SCTAB(9+(IBLK-1)*6)
        NCPL=SCTAB(10+(IBLK-1)*6)
        IF(NCPL.NE.0) THEN
         ND=SCTAB(11+(IBLK-1)*6)
         KSPCPL=SCTAB(12+(IBLK-1)*6)
         KSPDET=SCTAB(13+(IBLK-1)*6)
         LTRANS=SCTAB(14+(IBLK-1)*6)
         WRITE(u6,*)'------------------------------------------'
         WRITE(u6,'(1x,A,I16)')' Nr of open shells  :',IOPEN
         WRITE(u6,'(1x,A,I16)')' Nr of proto-CSF    :',NCPL
         WRITE(u6,'(1x,A,I16)')' Nr of proto-SD     :',ND
         WRITE(u6,'(1x,A,I16)')' Addr of proto-CSF  :',KSPCPL
         WRITE(u6,'(1x,A,I16)')' Addr of proto-SD   :',KSPDET
         WRITE(u6,'(1x,A,I16)')' Addr of transf matr:',LTRANS
         WRITE(u6,*)' proto-CSF''s:'
         CALL PRPCSF(IOPEN,NCPL,SCTAB(KSPCPL))
         WRITE(u6,*)' proto-SD''s:'
         CALL PRPDET(IOPEN,ND,SCTAB(KSPDET))
         WRITE(u6,*)' Transformation matrix:'
         CALL PRPTRA(ND,NCPL,TRANS)
        END IF
       END DO
      END IF
      WRITE(u6,*)'------------------------------------------'
      END SUBROUTINE PRSCTAB

      SUBROUTINE PRPCSF(IOPEN,NCPL,ICOUP)
      use definitions, only: iwp, u6
      IMPLICIT NONE
      Integer(kind=iwp), intent(in):: IOPEN, NCPL
      Integer(kind=iwp), intent(in):: ICOUP(IOPEN,NCPL)

      CHARACTER(LEN=1) :: CPLSMB(0:1)=['d','u']
      CHARACTER(LEN=24) FORM
      Integer(kind=iwp) N,I,J
      IF(IOPEN.LT.0 .OR.NCPL.LT.0) THEN
        Call WarningMessage(2,'Program bug: Erroneous call to PRPCSF.')
        WRITE(u6,*)'PRPCSF error: Wrong arguments.'
        WRITE(u6,*)'PRPCSF: IOPEN=',IOPEN
        WRITE(u6,*)'PRPCSF: NCPL =',NCPL
        CALL ABEND()
      END IF
      IF(IOPEN.EQ.0 .OR.NCPL.EQ.0) THEN
        Call WarningMessage(1,'Program bug? Strange call to PRPCSF.')
        WRITE(u6,*)'PRPCSF warning: Strange arguments.'
        WRITE(u6,*)'PRPCSF: IOPEN=',IOPEN
        WRITE(u6,*)'PRPCSF: NCPL =',NCPL
      ELSE
       N=80/(7+IOPEN)
       WRITE(FORM,'(A1,I2,A,I2,A4)') '(',N,'(1X,I5,1X,',IOPEN,'A1))'
       WRITE(u6,FORM)(I,(CPLSMB(ICOUP(J,I)),J=1,IOPEN),I=1,NCPL)
      END IF
      END SUBROUTINE PRPCSF

      SUBROUTINE PRPTRA(ND,NCPL,TRA)
      use definitions, only: iwp, wp, u6
      IMPLICIT NONE
      integer(kind=iwp), intent(in):: ND, NCPL
      Real(kind=wp), intent(in):: TRA(ND,NCPL)

      integer(kind=iwp) ISTA,IEND,ID,ICPL
      IF(ND.LT.0 .OR.NCPL.LT.0) THEN
        Call WarningMessage(2,'Program bug: Erroneous call to PRPTRA.')
        WRITE(u6,*)'PRPTRA error: Wrong arguments.'
        WRITE(u6,*)'PRPTRA: ND,NCPPL=',ND,NCPL
        CALL ABEND()
      END IF
      IF(ND.EQ.0 .OR.NCPL.EQ.0) THEN
        Call WarningMessage(1,'Program bug? Strange call to PRPCSF.')
        WRITE(u6,*)'PRPTRA warning: Strange arguments.'
        WRITE(u6,*)'PRPTRA: ND,NCPPL=',ND,NCPL
      ELSE
        DO ISTA=1,NCPL,5
         IEND=MIN(NCPL,ISTA+4)
         WRITE(u6,*)
         WRITE(u6,'(8x,5(I8,8X))')(ICPL,ICPL=ISTA,IEND)
         DO ID=1,ND
           WRITE(6,'(1x,5F16.8)')(TRA(ID,ICPL),ICPL=ISTA,IEND)
         END DO
        END DO
      END IF
      END SUBROUTINE PRPTRA

      SUBROUTINE PRPDET(IOPEN,ND,IDET)
      use definitions, only: iwp, u6
      IMPLICIT NONE
      integer(kind=iwp), intent(in):: IOPEN,ND
      Integer(kind=iwp), intent(in):: IDET(IOPEN,ND)

      CHARACTER(LEN=1) :: SPNSMB(0:1)=['b','a']
      CHARACTER(LEN=24) FORM
      Integer(kind=iwp) N,I,J
      IF(IOPEN.LT.0 .OR.ND.LT.0) THEN
        Call WarningMessage(2,'Program bug: Erroneous call to PRPDET.')
        WRITE(u6,*)'PRPDET error: Wrong arguments.'
        WRITE(u6,*)'PRPDET: IOPEN=',IOPEN
        WRITE(u6,*)'PRPDET: ND =',ND
        CALL ABEND()
      END IF
      IF(IOPEN.EQ.0 .OR.ND.EQ.0) THEN
        Call WarningMessage(1,'Program bug? Strange call to PRPDET.')
        WRITE(u6,*)'PRPDET warning: Strange arguments.'
        WRITE(u6,*)'PRPDET: IOPEN=',IOPEN
        WRITE(u6,*)'PRPDET: ND =',ND
      ELSE
       N=80/(7+IOPEN)
       WRITE(FORM,'(A1,I2,A,I2,A4)') '(',N,'(1X,I5,1X,',IOPEN,'A1))'
       WRITE(u6,FORM)(I,(SPNSMB(IDET(J,I)),J=1,IOPEN),I=1,ND)
      END IF
      END SUBROUTINE PRPDET
