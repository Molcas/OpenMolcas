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
C   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8

      SUBROUTINE Dynamix(iReturn)
      USE Isotopes
      IMPLICIT REAL*8 (a-h,o-z)
#include "Molcas.fh"
#include "warnings.fh"
#include "MD.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "dyn.fh"
#include "constants2.fh"
      INTEGER AixRm
      EXTERNAL IsFreeUnit,AixRm,IsStructure
      PARAMETER   (nTasks=3)
      PARAMETER  (nh=6)
      CHARACTER   StdIn*16, caption*15, ENV*8
      REAL*8      time, mean, kb
      REAL*8      Epot,Ekin,Etot0
      REAL*8      NHC(nh)
      INTEGER     Task(nTasks),natom,IsFreeUnit,irc,Itr,MxItr
      INTEGER     iRlxRoot,nRoots,i
      LOGICAL     Found,lHop
      INTEGER     VelVer, VV_First, VV_Second, Gromacs, VV_Dump
      PARAMETER   (au_time = CONST_AU_TIME_IN_SI_*1.0D15)
      PARAMETER   (kb = CONST_BOLTZMANN_/
     &             (CONV_AU_TO_KJ_*1.0D3))
      PARAMETER  (VelVer=1,VV_First=2,VV_Second=3,Gromacs=4,VV_Dump=5)
      PARAMETER  (iQ1=1,iQ2=2,iX1=3,iX2=4,iVx1=5,iVx2=6)
      CHARACTER, ALLOCATABLE :: atom(:)*2
      REAL*8, ALLOCATABLE ::    Mass(:),vel(:),pcoo(:,:)
      INTEGER Iso

*
      CALL QEnter('Dynamix')
      iReturn=99
*
C
C     Initialize Dynamix and set default values
C
#ifdef _DEBUG_
      WRITE(6,*)' Dynamix calls Init_Dynamix.'
#endif
      Call Init_Dynamix
#ifdef _DEBUG_
      WRITE(6,*)' Dynamix back from Init_Dynamix.'
#endif
C

C     Read the input
C
#ifdef _HDF5_
      call cre_dyn
#endif
#ifdef _DEBUG_
      WRITE(6,*)' Dynamix calls Readin_Dynamix.'
#endif
      CALL Readin_Dynamix(Task,nTasks,mTasks)
#ifdef _DEBUG_
      WRITE(6,*)' Dynamix back from Readin_Dynamix.'
#endif
C
C     Check if this is an initial run of Dynamix
C
      CALL Qpg_dScalar('MD_Time',Found)
C
#ifdef _HDF5_
      if (.not.found .and. lH5Restart) then
         call restart_dynamix(file_h5res)
         found = .true.
      endif
#endif

C     Generate or read velocities if this is an initial run
C
      IF (.NOT.Found) THEN
C     Check if the RESTART keyword was used.
         IF (RESTART.EQ.0.0D0) THEN
            time=0.000D0
         ELSE
            time=RESTART
            WRITE(6,'(5X,A,T55,F9.2,A)') 'MD restart time = ',
     &                                    RESTART, ' a.u.'
         END IF

         CALL Put_dScalar('MD_Time',time)
         CALL Get_nAtoms_Full(natom)

         CALL mma_allocate(atom,natom)
         CALL mma_allocate(Mass,natom)
         CALL mma_allocate(vel,natom*3)
         CALL mma_allocate(pcoo,POUT,natom*3)

         CALL Get_nAtoms_All(matom)
         CALL Get_Mass_All(Mass,matom)

C Initialize Thermostat Variables

         IF (THERMO.eq.2) THEN
            Freq = 1.D0/(2.2D1/au_time)
            Q1 = 3.D0*dble(natom)*TEMP*Kb/(Freq*Freq)
            Q2 = TEMP*Kb/(Freq*Freq)

            NHC(iQ1) = Q1
            NHC(iQ2) = Q2
            NHC(iX1) = 0.D0
            NHC(iX2) = 0.D0
            NHC(iVx1) = 0.D0
            NHC(iVx2) = 0.D0

            CALL Put_NHC(NHC,nh)
#ifdef _HDF5_
            call mh5_put_dset(dyn_nh,NHC)
#endif


         END IF

C Check if nuclear coordinates to project out from the dynamics
         IF (POUT.eq.0) THEN
            WRITE(6,'(5X,A,T55)') 'Dynamics in full dimensionality.'
         ELSE
            WRITE(6,'(5X,A,T55)') 'Dynamics in reduced dimensionality.'
            CALL DxRdOut(pcoo,POUT,natom)
C Save on RUNFILE
            CALL Put_dArray('Proj_Coord',pcoo,POUT*natom*3)
         ENDIF


         IF (VELO.eq.1) THEN
            CALL DxRdVel(vel,natom)
            WRITE(6,'(5X,A,T55)')
     &      'The initial velocities (bohr/au) are read in.'
         ELSEIF (VELO.eq.2) THEN
            CALL DxRdVel(vel,natom)
            DO i=1, natom
               IF (i.gt.matom) THEN
                  CALL LeftAd(atom(i))
                  Iso=0
                  CALL Isotope(Iso,atom(i),Mass(i))
               END IF
C-------------------------------------------

               DO j=1, 3
                  vel(3*(i-1)+j)=vel(3*(i-1)+j)/SQRT(Mass(i))
               END DO
            END DO
            WRITE(6,'(5X,A,T55)')
     &     'The initial mass weighted velocities (bohr/au) are read in.'

C Maxwell-Boltzmann distribution
         ELSEIF (VELO.eq.3) THEN
            nFlag=0
            val=0.d0
            buffer=0.D0
            CALL Get_Name_Full(atom)

C   . |  1    .    2    .    3    .    4    .    5    .    6    .    7 |  .    8

            WRITE(6,'(5X,A,T55)')
     & 'The initial velocities (bohr/au) are taken '
            WRITE(6,'(5X,A,f9.2,A)')
     & 'from a Boltzmann distribution at', TEMP, ' kelvin'

            CALL getSeed(iseed)

            DO i=1, natom
               IF (i.gt.matom) THEN
                  CALL LeftAd(atom(i))
                  Iso=0
                  CALL Isotope(Iso,atom(i),Mass(i))
               END IF
C-------------------------------------------

               arg=TEMP*Kb/Mass(i)
               Sigma=SQRT(arg)
               mean = 0.D0
               DO j=1, 3
                  CALL RandomGauss(mean,Sigma,iseed,nflag,buffer,Val)
                  vel(3*(i-1)+j)= Val

C                  WRITE(6,'(5X,A,T55,D16.8)') 'Vel = ', Val

               END DO
            END DO
         ELSE
            DO i=1, 3*natom
               vel(i)=0.000000000000D0
            END DO
            WRITE(6,'(5X,A,T55)')
     &      'The initial velocities are set to zero.'
         END IF
         CALL Get_Name_Full(atom)
         caption='Velocities'
         CALL DxPtTableWithoutMassForce(caption,time,natom,
     &        atom,vel)

C     Calculate the kinetic energy
         IF (VELO.gt.0) THEN
            Ekin=0.000000000000D0
            CALL Get_Name_Full(atom)
            DO i=1, natom
               IF (i.GT.matom) THEN
                  CALL LeftAd(atom(i))
                  Iso=0
                  CALL Isotope(Iso,atom(i),Mass(i))
               END IF
C-------------------------------------------
               DO j=1, 3
                  Ekin=Ekin+(5.0D-01)*Mass(i)*(vel(3*(i-1)+j)**2)
               END DO
            END DO
         ELSE
            Ekin=0.000000000000D0
         END IF
         WRITE(6,'(5X,A,6X,D19.12,1X,A)') 'Kinetic energy',Ekin,'a.u.'
C     Save the velocities on RUNFILE
         CALL Put_Velocity(vel,3*natom)
C     Save the total energy on RUNFILE if the total energy should be conserved.
         CALL Get_dScalar('Last Energy',Epot)
         Etot0 = Epot + Ekin
         CALL Put_dScalar('MD_Etot0',Etot0)
         CALL Put_dScalar('MD_Etot',Etot0)
#ifdef _HDF5_
         call mh5_put_dset(dyn_vel,vel)
         call mh5_put_dset(dyn_etot0,Etot0)
         call mh5_put_dset(dyn_etot,Etot0)
#endif
         CALL DxEnergies(time,Epot,Ekin,Etot0)
         WRITE(6,'(5X,A,8X,D19.12,1X,A)') 'Total Energy',Etot0,'a.u.'
         CALL mma_deallocate(atom)
         CALL mma_deallocate(Mass)
         CALL mma_deallocate(vel)
         CALL mma_deallocate(pcoo)
      END IF

C
C     Execute the tasks
C
      DO iTask = 1, mTasks

         IF (Task(iTask).eq.VelVer) THEN

            IF (Found) THEN

#ifdef _DEBUG_
      WRITE(6,*)' Dynamix calls VelVer_Second.'
#endif
               CALL VelVer_Second(irc)
#ifdef _DEBUG_
      WRITE(6,*)' Dynamix back from VelVer_Second.'
#endif
C
C     Check for Hopping?
C
               lHop=.FALSE.
               CALL qpg_iScalar('MaxHops',lHop)
               IF (lHop) THEN
C
C     Read the roots
C
                  CALL Get_iScalar('Number of roots',nRoots)
                  CALL Get_iScalar('Relax CASSCF root',iRlxRoot)
C
C     Run RASSI
C
                  LuInput=11
                  LuInput=IsFreeUnit(LuInput)
                  Call StdIn_Name(StdIn)
                  Call Molcas_Open(LuInput,StdIn)
                  Write (LuInput,'(A)')
     &                  '>export DYN_OLD_TRAP=$MOLCAS_TRAP'
                  Write (LuInput,'(A)') '>export MOLCAS_TRAP=ON'
                  Write (LuInput,'(A)') ' &RASSI &End'
                  Write (LuInput,'(A)') ' NR OF JOBIPHS'
                  Write (LuInput,*) ' 1 ',nRoots
                  Write (LuInput,*) (i,i=1,nRoots)
*                  Write (LuInput,'(X,I1,1X,I1)') inxtState,iRlxRoot
                  Write (LuInput,'(A)') ' HOP'
                  Write (LuInput,'(A)') 'End of Input'
                  Write (LuInput,'(A)') ' &Dynamix &End'
                  Write (LuInput,'(A)') ' VV_First'
                  Write (LuInput,'(A)') ' DT'
                  Write (LuInput,*) DT
                  Write (LuInput,'(A)') 'THERMO'
                  Write (LuInput,*) THERMO
                  Write (LuInput,'(A)') 'VELO'
                  Write (LuInput,*) VELO
                  Write (LuInput,'(A)') 'OUT'
                  Write (LuInput,*) POUT
                  Write (LuInput,'(A)') 'End of Input'
                  Write (LuInput,'(A)')
     &                  '>export MOLCAS_TRAP=$DYN_OLD_TRAP'
                  Close(LuInput)
                  Call Finish(_RC_INVOKED_OTHER_MODULE_)
               ELSE
#ifdef _DEBUG_
      WRITE(6,*)' Dynamix calls VelVer_First.'
#endif
                  CALL VelVer_First(irc)
               END IF

#ifdef _DEBUG_
      WRITE(6,*)' Dynamix back from VelVer_First.'
#endif
            ELSE
#ifdef _DEBUG_
      WRITE(6,*)' Dynamix calls VelVer_First.'
#endif

               CALL VelVer_First(irc)
#ifdef _DEBUG_
      WRITE(6,*)' Dynamix back from VelVer_First.'
#endif
            END IF
         ELSE IF (Task(iTask).eq.VV_First) THEN
             CALL VelVer_First(irc)

         ELSE IF (Task(iTask).eq.VV_Second) THEN
             CALL VelVer_Second(irc)

         ELSE IF (Task(iTask).eq.Gromacs) THEN
            CALL GROM(irc)

         ELSE IF (Task(iTask).eq.VV_Dump) THEN
            CALL VelVer_Dump(irc)

         ELSE
            WRITE(6,*) 'Illegal task'
            CALL QTrace()
            CALL Abend()
         END IF
      END DO

C
C-----Remove the GRADS file
C
      Call f_Inquire('GRADS',Found)
      If (Found) Then
         If (AixRm('GRADS').ne.0) Call Abend()
      End If

#ifdef _HDF5_
      call mh5_close_file(dyn_fileid)
#endif

C
C-----If running in a DoWhile loop, we turn a successful
C     return code into "continue loop", except on the
C     last iteration
C
      If ((IsStructure().eq.1).and.(irc.eq.0)) Then
         MxItr=0
         Call GetEnvf('MOLCAS_MAXITER', ENV)
         If (ENV.ne.' ') Then
            Read (ENV,*) MxItr
         End If
         Itr=1
         Call GetEnvf('MOLCAS_ITER', ENV)
         If (ENV.ne.' ') Then
            Read (ENV,*) Itr
         End If
         If (Itr.lt.MxItr) Then
            iReturn=_RC_CONTINUE_LOOP_
         Else
            iReturn=irc
         End If
      Else
         iReturn=irc
      End If
      CALL QExit('Dynamix')
      RETURN
*
      END
