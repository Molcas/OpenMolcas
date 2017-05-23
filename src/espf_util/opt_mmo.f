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
#ifdef _GROMACS_
      SUBROUTINE Opt_MMO(nAtIn,Coord,nAtOut,CoordMMO,nAtGMX,AT,ipGMS)

      IMPLICIT NONE

#include "espf.fh"
#include "opt_mmo.fh"
#include "stdalloc.fh"

      INTEGER, INTENT(IN) :: ipGMS,nAtGMX,nAtIn,nAtOut
      INTEGER, DIMENSION(nAtGMX), INTENT(IN) :: AT
      REAL*8, DIMENSION(3,nAtIn), INTENT(IN) :: Coord
      REAL*8, DIMENSION(3,nAtOut), INTENT(INOUT) :: CoordMMO

      INTEGER :: i,iAtIn,iAtOut,iOk,iPL,j,MMIter
      REAL*8 :: EnergyGMX,MaxF,OldEn,Step
      REAL*8, PARAMETER :: TinyStep = 1.0D-50*AuToNm
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: CoordGMX,FieldGMX,ForceGMX
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: GradMMO,NewCoord,OldCoord
      REAL*8, DIMENSION(:), ALLOCATABLE :: PotGMX
      CHARACTER(LEN=256) :: Message

      INTEGER, EXTERNAL :: iPL_espf,mmslave_calc_energy
      REAL*8, EXTERNAL :: ddot_

      CALL qEnter('Opt_MMO')

      iPL = iPL_espf()

      IF (nAtOut==0) THEN
         Message = 'Opt_MMO: subroutine called with nAtOut=0'
         CALL WarningMessage(2,Message)
         CALL Abend()
      END IF

! Set up arrays for use with Gromacs
      CALL mma_allocate(CoordGMX,3,nAtGMX)
      CALL mma_allocate(ForceGMX,3,nAtGMX)
      CALL mma_allocate(FieldGMX,3,nAtGMX)
      CALL mma_allocate(PotGMX,nAtGMX)
      iAtIn = 1
      iAtOut = 1
      DO i = 1,nAtGMX
         IF (AT(i)==QM.OR.AT(i)==MMI) THEN
            CALL dcopy_(3,Coord(1,iAtIn),1,CoordGMX(1,i),1)
            iAtIn = iAtIn+1
         ELSE
            CALL dcopy_(3,CoordMMO(1,iAtOut),1,CoordGMX(1,i),1)
            iAtOut = iAtOut+1
         END IF
      END DO
      CALL dscal_(3*nAtGMX,AuToNm,CoordGMX,1)

! Set up arrays with MMO coordinates and gradient
      CALL mma_allocate(NewCoord,3,nAtOut)
      CALL mma_allocate(OldCoord,3,nAtOut)
      CALL mma_allocate(GradMMO,3,nAtOut)
      CALL dcopy_(3*nAtOut,CoordMMO,1,NewCoord,1)
      CALL dscal_(3*nAtOut,AuToNm,NewCoord,1)

      IF (iPL>=2) THEN
100      FORMAT(I6,3(F12.6,1X),2X,A)
         CALL CollapseOutput(1,'Gromacs microiterations')
         WRITE(6,*)
         WRITE(6,*) 'Initial coordinates (angstrom)'
         WRITE(6,*) '------------------------------'
         DO i = 1,nAtGMX
            IF (AT(i)==QM) THEN
               WRITE(6,100) i,(CoordGMX(j,i)*NmToAng,j=1,3),'QM'
            ELSE IF (AT(i)==MMI) THEN
               WRITE(6,100) i,(CoordGMX(j,i)*NmToAng,j=1,3),'MMI'
            ELSE
               WRITE(6,100) i,(CoordGMX(j,i)*NmToAng,j=1,3),'MMO'
            END IF
         END DO
         WRITE(6,*)
      END IF

! MMO optimization cycle, performed in Gromacs units
      Step = 0.1D0*AuToNm
      OldEn = HUGE(OldEn)
      MaxF = HUGE(MaxF)
      MMIter = 0
      DO WHILE ((MMIter<MMIterMax).AND.(MaxF>ConvF).AND.(Step>TinyStep))
         MMIter = MMIter+1
         ! Get gradient from Gromacs
         iOk = mmslave_calc_energy(%val(ipGMS),CoordGMX,ForceGMX,
     &                             FieldGMX,PotGMX,EnergyGMX)
         IF (iOk/=1) THEN
            Message = 'Opt_MMO: mmslave_calc_energy is not ok'
            CALL WarningMessage(2,Message)
            CALL Abend()
         END IF
         iAtOut = 1
         DO i = 1,nAtGMX
            IF (AT(i)==MMO) THEN
               GradMMO(1,iAtOut) = -ForceGMX(1,i)
               GradMMO(2,iAtOut) = -ForceGMX(2,i)
               GradMMO(3,iAtOut) = -ForceGMX(3,i)
               iAtOut = iAtOut+1
            END IF
         END DO
#ifdef _DEBUG_
200      FORMAT(I5,6ES12.4)
         IF (MMIter==1) THEN
            WRITE(6,*)
     &         'Iter       E          |F|         Fmax       Step'
            WRITE(6,*)
     &         '                     (Atomic units)'
            WRITE(6,*)
     &         '-----------------------------------------------------'
         END IF
         IF ((MOD(MMIter,10)==0).OR.(MMIter==1)) THEN
            WRITE(6,200) MMIter,EnergyGMX/AuToKjPerMol,
     &         SQRT(DDot_(3*nAtOut,GradMMO,1,GradMMO,1))/AuToKjPerMolNm,
     &         MaxF/AuToKjPerMolNm,Step/AuToNm
         END IF
#endif
         ! Steepest descent with adaptive step a la Gromacs
         IF ((EnergyGMX<OldEn).OR.(MMIter==1)) THEN
            Step = 1.2D0*Step
            MaxF = Zero
            DO iAtOut = 1,nAtOut
               DO j = 1,3
                  MaxF = MAX(MaxF,ABS(GradMMO(j,iAtOut)))
               END DO
            END DO
            OldEn = EnergyGMX
            CALL dcopy_(3*nAtOut,NewCoord,1,OldCoord,1)
            CALL daxpy_(3*nAtOut,-Step/MaxF,GradMMO,1,NewCoord,1)
         ELSE
            Step = 0.2D0*Step
            CALL dcopy_(3*nAtOut,OldCoord,1,NewCoord,1)
            CALL daxpy_(3*nAtOut,-Step/MaxF,GradMMO,1,NewCoord,1)
         END IF
         ! Update coordinates for Gromacs
         iAtOut = 1
         DO i = 1,nAtGMX
            IF (AT(i)==MMO) THEN
               CoordGMX(1,i) = NewCoord(1,iAtOut)
               CoordGMX(2,i) = NewCoord(2,iAtOut)
               CoordGMX(3,i) = NewCoord(3,iAtOut)
               iAtOut = iAtOut+1
            END IF
         END DO
      END DO

! Issue warning if calculation stopped for other reason than small F
      IF (MaxF>ConvF) THEN
         IF (MMIter==MMIterMax) THEN
            Message = 'Maximum number of microiterations reached'
            CALL WarningMessage(1,Message)
         ELSE IF (Step<=TinyStep) THEN
            Message = 'Microiterations stopped due to zero step size'
            CALL WarningMessage(1,Message)
         END IF
      END IF

! Undo the last move to recover the "converged" coordinates
      IF (MMIter>0) THEN
         CALL dcopy_(3*nAtOut,OldCoord,1,NewCoord,1)
         iAtOut = 1
         DO i = 1,nAtGMX
            IF (AT(i)==MMO) THEN
               CoordGMX(1,i) = NewCoord(1,iAtOut)
               CoordGMX(2,i) = NewCoord(2,iAtOut)
               CoordGMX(3,i) = NewCoord(3,iAtOut)
               iAtOut = iAtOut+1
            END IF
         END DO
      END IF

      IF (iPL>=2) THEN
         WRITE(6,*)
300      FORMAT(' Performed: ',I6,' MM iterations')
         WRITE(6,300) MMIter
400      FORMAT(' Max. force: ',F12.6,' (convergence at: ',F12.6,')')
         WRITE(6,400) MaxF/AuToKjPerMolNm,ConvF/AuToKjPerMolNm
         CALL dscal_(3*nAtGMX,-One/AuToKjPerMolNm,ForceGMX,1)
         EnergyGMX = EnergyGMX/AuToKjPerMol
#ifdef _DEBUG_
500      FORMAT(3(F12.6,1X))
         WRITE(6,*)
         WRITE(6,*) 'Properties'
         WRITE(6,*) '----------'
         WRITE(6,*) 'Energy: ',EnergyGMX
         DO i = 1,nAtGMX
            IF (AT(i)==MMO) THEN
               WRITE(6,500) (ForceGMX(j,i),j=1,3)
            ELSE
               WRITE(6,500) Zero,Zero,Zero
            END IF
         END DO
#endif
         WRITE(6,*)
      END IF

! Put optimized MMO coordinates on runfile
      CALL dcopy_(3*nAtOut,NewCoord,1,CoordMMO,1)
      CALL dscal_(3*nAtOut,One/AuToNm,CoordMMO,1)
      CALL Put_dArray('MMO Coords',CoordMMO,3*nAtOut)

      IF (iPL>=2) THEN
         IF (MMIter>1) THEN
            WRITE(6,*) 'Final coordinates (angstrom)'
            WRITE(6,*) '----------------------------'
            DO i = 1,nAtGMX
               IF (AT(i)==QM) THEN
                  WRITE(6,100) i,(CoordGMX(j,i)*NmToAng,j=1,3),'QM'
               ELSE IF (AT(i)==MMI) THEN
                  WRITE(6,100) i,(CoordGMX(j,i)*NmToAng,j=1,3),'MMI'
               ELSE
                  WRITE(6,100) i,(CoordGMX(j,i)*NmToAng,j=1,3),'MMO'
               END IF
            END DO
         END IF
         CALL CollapseOutput(0,'Gromacs microiterations')
         WRITE(6,*)
      END IF

! Clean up
      CALL mma_deallocate(CoordGMX)
      CALL mma_deallocate(ForceGMX)
      CALL mma_deallocate(FieldGMX)
      CALL mma_deallocate(PotGMX)
      CALL mma_deallocate(NewCoord)
      CALL mma_deallocate(OldCoord)
      CALL mma_deallocate(GradMMO)

      CALL qExit('Opt_MMO')

      RETURN
      END
#elif defined (NAGFOR)
! Some compilers do not like empty files
      SUBROUTINE empty_Opt_MMO()
      END
#endif
