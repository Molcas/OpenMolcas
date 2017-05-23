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
      SUBROUTINE PRWVF(IORBTAB,ISSTAB,IFSBTAB,PRTHR,CI)
      IMPLICIT NONE
      REAL*8 CI(*),PRTHR
      INTEGER IORBTAB(*),ISSTAB(*),IFSBTAB(*)
C     INTEGER MORSBITS,IBLKDET,NPRINTED
      INTEGER          IBLKDET,NPRINTED
C     INTEGER NASORB,NASPRT,KSPART,KOINFO,KSSTTB,KSBSMRS
      INTEGER        NASPRT,KSPART,       KSSTTB,KSBSMRS
      INTEGER NFSB,KSTARR,LSBSET,NSSTP,ISUM,ISST,NSBS
      INTEGER IFSB,KPOS,ISPART,NBLKDET
      INTEGER ISSTARR(50),NDIMARR(50)
      INTEGER ICREL,ICI,INDX,ISORB,N,I,ISBS,IMORS
      CHARACTER*80 DETTXT
#include "WrkSpc.fh"

C The orbital table:
      NASPRT= IORBTAB(9)
      KSPART= IORBTAB(10)
C The substring table:
      NSSTP   =ISSTAB(7)
      KSSTTB=15
      KSBSMRS=ISSTAB(11)
C The FS blocks of the SGM wave function:
      NFSB=IFSBTAB(3)
      KSTARR=8
C Make an array with nr of earlier substrings for each
C substring type:
      CALL GETMEM('NSBSET','Allo','Inte',LSBSET,NSSTP)
      ISUM=0
      DO ISST=1,NSSTP
        IWORK(LSBSET-1+ISST)=ISUM
        NSBS=ISSTAB(KSSTTB+5*(ISST-1))
        ISUM=ISUM+NSBS
      END DO
C Loop over FS blocks of the SGM wave function
      NPRINTED=0
      DO IFSB=1,NFSB
        KPOS=KSTARR+(NASPRT+2)*(IFSB-1)
        DO ISPART=1,NASPRT
          ISSTARR(ISPART)=IFSBTAB(KPOS-1+ISPART)
        END DO
        NBLKDET =IFSBTAB(KPOS+NASPRT  )
        IBLKDET =IFSBTAB(KPOS+NASPRT+1)
C Dimension of each substring type:
        DO ISPART=1,NASPRT
          ISST=ISSTARR(ISPART)
          NSBS=ISSTAB(KSSTTB+5*(ISST-1))
          NDIMARR(ISPART)=NSBS
        END DO
        DO ICREL=1,NBLKDET
         ICI=IBLKDET-1+ICREL
         IF(ABS(CI(ICI)).GE.PRTHR) THEN
C Get occupation array in the form of string DETTXT:
           INDX=ICREL-1
           ISORB=0
           DO ISPART=1,NASPRT
            N=NDIMARR(ISPART)
            I=MOD(INDX,N)
            INDX=INDX-N*I
            ISST=ISSTARR(ISPART)
            ISBS=IWORK(LSBSET-1+ISST)+I+1
            IMORS=ISSTAB(KSBSMRS+2*(ISBS-1))
            N=IORBTAB(KSPART-1+ISPART)
            CALL MORSWRITE(IMORS,DETTXT(ISORB+1:ISORB+N))
            ISORB=ISORB+N
           END DO
           WRITE(6,'(1x,a,5x,f16.8)') DETTXT(1:ISORB),CI(ICI)
           NPRINTED=NPRINTED+1
         END IF
        END DO
C End of loop over FS blocks
      END DO
      IF(NPRINTED.EQ.0) WRITE(6,*)' (PRWVF: Nothing worth printing)'
      CALL GETMEM('NSBSET','free','Inte',LSBSET,NSSTP)
      RETURN
      END
