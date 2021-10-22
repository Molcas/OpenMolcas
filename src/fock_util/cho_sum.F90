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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************
      SUBROUTINE CHO_SUM(rc,nSym,nBas,iUHF,DoExchange,FLT,FSQ)

!****************************************************************
!  Author : F. Aquilante
!
!  Purpose:
!           Accumulates the Coulomb and Exchange contribution
!           to the frozen AO-Fock matrices for alpha and beta
!           spin as defined in the calling routine
!*****************************************************************
      use Data_Structures, only: DSBA_Type
      Implicit Real*8 (a-h,o-z)
      Integer   rc,nBas(8),nSym,iUHF

      Type (DSBA_Type) FLT(*), FSQ(*)
      Logical DoExchange(*)

!*************************************************
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
!*************************************************

      if (iUHF.eq.1)then
         nDen=3
      else
         nDen=1
      endif

! Accumulate the contributions and Square the final matrix
! FLT is in lower triangular storage
! FSQ is in squared storage
!
! the lower triangular part of FSQ is added to FLT
!
      IF(nDen.eq.1) THEN

      DO ISYM=1,NSYM
       NB = NBAS(ISYM)
       IF (NB.gt.0) THEN
       If (DoExchange(1)) Then
        DO IB=1,NB
          DO JB=IB,NB
             IJB=iTri(JB,IB)
             FLT(1)%SB(ISYM)%A1(IJB)= FLT(1)%SB(ISYM)%A1(IJB)           &
     &                              + FSQ(1)%SB(ISYM)%A2(JB,IB)
          END DO
        END DO
       End If
       CALL SQUARE(FLT(1)%SB(ISYM)%A1,FSQ(1)%SB(ISYM)%A2,1,NB,NB)
       ENDIF
      END DO


      ELSE  ! nDen=3

      DO ISYM=1,NSYM
       NB = NBAS(ISYM)
       IF (NB.gt.0) THEN
       If (DoExchange(2)) Then
        DO IB=1,NB
          DO JB=IB,NB
            IJB=iTri(JB,IB)
            FLT(1)%SB(ISYM)%A1(IJB) = FLT(1)%SB(ISYM)%A1(IJB)           &
     &                              + FSQ(2)%SB(ISYM)%A2(JB,IB)
            FLT(2)%SB(ISYM)%A1(IJB) = FLT(2)%SB(ISYM)%A1(IJB)           &
     &                              + FSQ(3)%SB(ISYM)%A2(JB,IB)
          END DO
        END DO
       End If

       CALL SQUARE(FLT(1)%SB(ISYM)%A1,FSQ(2)%SB(ISYM)%A2,1,NB,NB)
       CALL SQUARE(FLT(2)%SB(ISYM)%A1,FSQ(3)%SB(ISYM)%A2,1,NB,NB)

       ENDIF
      END DO

      ENDIF  ! nDen=3


! Print the Fock-matrix
#ifdef _DEBUGPRINT_
      WRITE(6,'(6X,A)')'TEST PRINT FROM CHO_SUM.'
      WRITE(6,'(6X,A)')'FROZEN FOCK MATRIX IN AO BASIS.'

      if (nDen.gt.1) then

      do jDen=1,2
      if(jDen.eq.1)WRITE(6,'(6X,A)')'SPIN ALPHA'
      if(jDen.eq.2)WRITE(6,'(6X,A)')'SPIN BETA'
      DO ISYM=1,NSYM
        NB=NBAS(ISYM)
        IF ( NB.GT.0 ) THEN
          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
          CALL TRIPRT(' ',' ',FLT(jDen)%SB(ISYM)%A1,NB)
        END IF
      END DO
      end do

      else ! nDen=1

      DO ISYM=1,NSYM
        NB=NBAS(ISYM)
        IF ( NB.GT.0 ) THEN
          WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
          CALL TRIPRT(' ',' ',FLT(1)%SB(ISYM)%A1,NB)
        END IF
      END DO

      endif

#endif

      rc=0

      Return
      END

!*************************************************************
