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
* Copyright (C) Francesco Aquilante                                    *
************************************************************************
      SUBROUTINE CHO_SUM(rc,nSym,nBas,iUHF,DoExchange,
     &                  ipFLT,ipFSQ)

*****************************************************************
*  Author : F. Aquilante
*
*  Purpose:
*           Accumulates the Coulomb and Exchange contribution
*           to the frozen AO-Fock matrices for alpha and beta
*           spin as defined in the calling routine
******************************************************************

      Implicit Real*8 (a-h,o-z)
      Integer   rc,nBas(8),nSym,iUHF
      Integer   ISTSQ(8),ISTLT(8)
      Integer   ipFLT(3),ipFSQ(3)
      Logical DoExchange(3),Debug


#include "WrkSpc.fh"

**************************************************
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
**************************************************

#ifdef _DEBUG_
      Debug=.true.
#else
      Debug=.false.
#endif

       if(iUHF.eq.1)then
               nDen=3
       else
               nDen=1
       endif


c ISTSQ: Offsets to full symmetry block in DSQ,FSQ
c ISTLT: Offsets to packed LT symmetry blocks in DLT,FLT
        ISTSQ(1)=0
        ISTLT(1)=0
      DO ISYM=2,NSYM
        NB=NBAS(ISYM-1)
        NB2=NB*NB
        NB3=(NB2+NB)/2
        ISTSQ(ISYM)=ISTSQ(ISYM-1)+NB2
        ISTLT(ISYM)=ISTLT(ISYM-1)+NB3
      END DO


c Accumulate the contributions and Square the final matrix
c FLT is in lower triangular storage
c FSQ is in squared storage
c
c the lower triangular part of FSQ is added to FLT
c
      IF(nDen.eq.1) THEN

      DO ISYM=1,NSYM
       NB = NBAS(ISYM)
       IF (NB.gt.0) THEN
       koff1 = ipFLT(1) - 1
       koff2 = ipFSQ(1) - 1
       If (DoExchange(1)) Then
        K1 = koff1 + ISTLT(ISYM)
        K2 = koff2 + ISTSQ(ISYM)
        DO IB=1,NB
          DO JB=IB,NB
              Work(K1+iTri(JB,IB)) = Work(K1+iTri(JB,IB)) + Work(K2+JB)
          END DO
          K2 = K2 + NB
        END DO
       End If
      CALL SQUARE(Work(ipFLT(1)+ISTLT(ISYM)),Work(ipFSQ(1)+ISTSQ(ISYM))
     &            ,1,NB,NB)
       ENDIF
      END DO


      ELSE  ! nDen=3

      DO ISYM=1,NSYM
       NB = NBAS(ISYM)
       IF (NB.gt.0) THEN
       koff0= ipFLT(1) - 1
       koff1= ipFLT(2) - 1
       koff2= ipFSQ(2) - 1
       koff3= ipFSQ(3) - 1
       If (DoExchange(2)) Then
        K0 = koff0 + ISTLT(ISYM)
        K1 = koff1 + ISTLT(ISYM)
        K2 = koff2 + ISTSQ(ISYM)
        K3 = koff3 + ISTSQ(ISYM)
        DO IB=1,NB
          DO JB=IB,NB
            Work(K0+iTri(JB,IB)) = Work(K0+iTri(JB,IB)) + Work(K2+JB)
            Work(K1+iTri(JB,IB)) = Work(K1+iTri(JB,IB)) + Work(K3+JB)
          END DO
          K2 = K2 + NB
          K3 = K3 + NB
        END DO
       End If

       CALL SQUARE(Work(ipFLT(1)+ISTLT(ISYM)),Work(ipFSQ(2)+ISTSQ(ISYM))
     &            ,1,NB,NB)
       CALL SQUARE(Work(ipFLT(2)+ISTLT(ISYM)),Work(ipFSQ(3)+ISTSQ(ISYM))
     &            ,1,NB,NB)

       ENDIF
      END DO

      ENDIF  ! nDen=3


c Print the Fock-matrix
#ifdef _DEBUG_
      WRITE(6,'(6X,A)')'TEST PRINT FROM CHO_SUM.'
      WRITE(6,'(6X,A)')'FROZEN FOCK MATRIX IN AO BASIS.'

      if (nDen.gt.1) then

      do jDen=1,2
      if(jDen.eq.1)WRITE(6,'(6X,A)')'SPIN ALPHA'
      if(jDen.eq.2)WRITE(6,'(6X,A)')'SPIN BETA'
      icount=0
      DO ISYM=1,NSYM
       ILT=ipFLT(jDen)+icount
        NB=NBAS(ISYM)
        IF ( NB.GT.0 ) THEN
          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
          CALL TRIPRT(' ',' ',Work(ILT),NB)
          icount=icount+NB*(NB+1)/2
        END IF
      END DO
      end do

      else ! nDen=1

      icount=0
      DO ISYM=1,NSYM
      ILT=ipFLT(1)+icount
        NB=NBAS(ISYM)
        IF ( NB.GT.0 ) THEN
          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
          CALL TRIPRT(' ',' ',Work(ILT),NB)
          icount=icount+NB*(NB+1)/2
        END IF
      END DO

      endif

#endif

      rc=0

      Return
      END

**************************************************************
