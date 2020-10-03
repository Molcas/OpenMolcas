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
      SUBROUTINE Tully(CIBigArray, NSTATE, NCI)
      use Tully_variables
      implicit none
#include "surfacehop.fh"
      INTEGER      values(8)
      character*8  date
      character*10 time
      character*5  zone
      real*8       Random_Molcas
      EXTERNAL     Random_Molcas
      logical    HOPPED, normalTully, found, lmaxHop, lnhop
      integer    NSTATE, NCI, maxhop, nhop
      real*8     DT, LO, EKIN, TAU(NSTATE)
      real*8     CIBigArray(NCI*NSTATE), CIBigArrayP(NCI*NSTATE)
      real*8     CIBigArrayPP(NCI*NSTATE), Etot, ediffcheck
      real*8     Dmatrix(NSTATE,NSTATE), sp(NSTATE,NSTATE)
      real*8     D32matrix(NSTATE,NSTATE), D12matrix(NSTATE,NSTATE)
      real*8     ExtrSlope(NSTATE,NSTATE), ExtrInter(NSTATE,NSTATE)
      real*8     VenergySlope(NSTATE), tempVector(NSTATE)
      real*8     tempVector2(NSTATE), VenergyInter(NSTATE)
      real*8     V(NSTATE,NSTATE), Bmatrix(NSTATE,NSTATE)
      real*8     Gprobab(NSTATE), Popul(NSTATE)
      real*8     VenergyP(NSTATE), Venergy(NSTATE), temp
      real*8     SumProb, scalarprod, prod,populOS
      integer    k,l,j,i,ii,jjj
      integer    rightOrder(NSTATE), decVec(NSTATE)
      integer    nstatesq, nciquery, stateRi, temproot, nsatom
      integer    ISTATE2, iseed, irlxroot
      complex*16 Amatrix(NSTATE,NSTATE), AmatrixDT(NSTATE,NSTATE)
      complex*16 ArelaxPrev


      CIBigArrayP(:)=0.0D0
      CIBigArrayPP(:)=0.0D0
      V(:,:) = 0.0D0
      TAU(:)=0.0D0
      write(6,*) ''
      write(6,*) '------------------------------------------'
      write(6,*) '            TULLY ALGORITHM'
      write(6,*) '------------------------------------------'
      write(6,*) ''

      call get_darray('Last energies',Venergy,NSTATE)
      call Get_iScalar('Relax CASSCF root',iRlxRoot)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C             Watch if dt exists.  then if Amatrix exist               C
C        if exist, it gets it, if not, it creates a new one            C
C                                                                      C

      call qpg_dscalar('Timestep',Found)

#ifdef _HDF5_
      if (.not.found .and. lH5Restart) then
         call restart_surfacehop
         found = .true.
      endif
#endif

      if (Found) then
         call Get_dScalar('Timestep',DT)
      end if

      call Qpg_zArray('AmatrixV', Found, nStateSq)
      write(6,*) 'Did the density matrix exists? ', Found
      IF (.NOT.Found) THEN
        do i=1, NSTATE
           do j=1, NSTATE
              Amatrix(i,j)=(0.0, 0.0)
           end do
        end do
        Amatrix(iRlxRoot,iRlxRoot)=(1.0, 0.0)

        call Put_zArray('AmatrixV',Amatrix,NSTATE*NSTATE)
      ELSE
        call get_zarray('AmatrixV',Amatrix,NSTATE*NSTATE)
      END IF

      call Qpg_dArray('AllCIP',Found,nCiQuery)
      write(6,*) 'Did the Pre-coefficients array exists? ', Found
      if (.NOT.Found) then
         call put_darray('AllCIP', CIBigArray, NCI*NSTATE)
         call put_dArray('VenergyP', Venergy, NSTATE)
         do i=1,NSTATE
            Popul(i)=Real(Amatrix(i,i))
         end do
         write(6,*) 'Gnuplot:', (Popul(j), j=1,NSTATE,1),
     &      (Venergy(j), j=1,NSTATE,1 ), Venergy(iRlxRoot)
         write(6,*) 'Cannot do deltas at first step, see you later! '
         RETURN
      else
         call Get_dArray('AllCIP', CIBigArrayP,NCI*NSTATE)
         call Get_dArray('VenergyP', VenergyP, NSTATE)
      end if

C now check for the CI coefficients at Pre-Pre-Step (PP)

      call Qpg_dArray('AllCIPP',Found,nCiQuery)
      write(6,*) 'Did the Pre-Pre-coefficients array exists? ', Found
      if (.NOT.Found) then
         normalTully=.true.
         call put_darray('AllCIPP', CIBigArrayP, NCI*NSTATE)
         call put_darray('AllCIP', CIBigArray, NCI*NSTATE)
         write(6,*) 'At second step we will use normal Tully Algorithm:'
      else
         normalTully=.false.
         call Get_dArray('AllCIPP', CIBigArrayPP,NCI*NSTATE)
      end if


      call Get_dScalar('MD_Etot',Etot)
      call Get_iScalar('Unique atoms',nsAtom)

      write(6,*) 'Density Matrix elements (i=1..#states):'
      do i=1, NSTATE
         write(6,*)( Amatrix(i,j), j=1,NSTATE,1 )
      end do


C Timestep:                  DT
C Total Energy               Etot
C Coefficients:              CIBigArray(i)      length = NCI*NSTATE
C Prev step coefficients:    CIBigArrayP(i)     length = NCI*NSTATE
C Prev-Prev coefficients:    CIBigArrayPP(i)    length = NCI*NSTATE
C V energy:                  Venergy(i)         length = NSTATE
C V Prev step energy:        VenergyP(i)        length = NSTATE
C MatriX A:                  Amatrix(i,j)       length = NSTATE*NSTATE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                        C
C                         Sign Corrector                                 C
C                                                                        C

C so first of all I create 2 temp vectors that store the absolute value
C of the energy difference

      do i=1, NSTATE
            tempVector(i) = ABS(Venergy(i)-Venergy(irlxRoot))
            tempVector2(i) = ABS(Venergy(i)-Venergy(irlxRoot))
      end do

C then I sort one of them, (relaxroot becomes first, it's zero)

      do j = 1,(NSTATE-1)
        do k = (j+1), NSTATE
           if (tempVector(j) .gt. tempVector(k)) then
              temp=tempVector(j)
              tempVector(j)=tempVector(k)
              tempVector(k)=temp
           end if
        end do
      end do

C I get the right order I need to process roots

      do i=1,NSTATE
         do j=1,NSTATE
            if (tempVector(i) .eq. tempVector2(j)) then
             rightOrder(i)= j
            end if
         end do
      end do
C ii counter on CURRENT STEP
      do ii=1,NSTATE
       do i=1,NSTATE
        scalarprod = 0.0
         do j=1,NCI
        scalarprod=scalarprod+
     &      CIBigArray(NCI*(ii-1)+j)*CIBigArrayP(NCI*(i-1)+j)
         end do
        sp(ii,i)=scalarprod
       end do
        decVec(ii)=1
      end do

      do ii=1,NSTATE
         stateRi=rightOrder(ii)
         prod=0.0
         jjj=0
         do i=1,NSTATE
            if (decVec(i) .eq. 1) then
               if (ABS(sp(stateRi,i)) .gt. ABS(prod)) then
                  prod = sp(stateRi,i)
                  jjj = i
               end if
            end if
         end do
         decVec(jjj)=0
         if (prod .lt. 0) then
             do k=1,NCI
          CIBigArray(NCI*(stateRi-1)+k)=-CIBigArray(NCI*(stateRi-1)+k)
             end do
         end if
      end do
C                                                                       C
C                         end of sign corrector                         C
C                                                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C create D matrix normal tully

      IF (normalTully) THEN

      write (6,*) ''
      write (6,*) 'Executing Normal Tully !!'
      write (6,*) ''

      do i=1, NSTATE
         do j=1, NSTATE
            if (i.ne.j) then
              Dmatrix(i,j) = 0.0
              do ii=1, NCI
                 Dmatrix(i,j) = Dmatrix(i,j) +
     &    CIBigArray(NCI*(i-1)+ii)*CIBigArrayP(NCI*(j-1)+ii)
              end do
              Dmatrix(i,j) = -Dmatrix(i,j)/DT
            else
              Dmatrix(i,i)=0.0
            end if
         end do
      end do
      normalTully=.false.

      do i=1, NSTATE
         do j=1, NSTATE
           ExtrInter(i,j)= Dmatrix(i,j)
           ExtrSlope(i,j)= 0.0
         end do
           VenergyInter(i)=VenergyP(i)
           VenergySlope(i)=(Venergy(i)-VenergyP(i))/DT
      end do

      ELSE

C Create D matrix according to Hammes-Schiffer-Tully (interpolating extrapolating)

C      D32matrix

      do i=1, NSTATE
         do j=1, NSTATE
            if (i.ne.j) then
            D32matrix(i,j) = 0.0
            do ii=1, NCI
                D32matrix(i,j) = D32matrix(i,j) +
     &            CIBigArrayPP(NCI*(i-1)+ii)*CIBigArrayP(NCI*(j-1)+ii)
                D32matrix(i,j) = D32matrix(i,j) -
     &            CIBigArrayP(NCI*(i-1)+ii)*CIBigArrayPP(NCI*(j-1)+ii)
            end do
            D32matrix(i,j) = D32matrix(i,j)/(2*DT)
            else
            D32matrix(i,i)=0.0
            end if
         end do
      end do

C      D12matrix

      do i=1, NSTATE
         do j=1, NSTATE
            if (i.ne.j) then
            D12matrix(i,j) = 0.0
            do ii=1, NCI
                D12matrix(i,j) = D12matrix(i,j) +
     &              CIBigArrayP(NCI*(i-1)+ii)*CIBigArray(NCI*(j-1)+ii)
                D12matrix(i,j) = D12matrix(i,j) -
     &              CIBigArray(NCI*(i-1)+ii)*CIBigArrayP(NCI*(j-1)+ii)
            end do
            D12matrix(i,j) = D12matrix(i,j)/(2*DT)
            else
            D12matrix(i,i)=0.0
            end if
         end do
      end do

C definition of Y intercept (ExtrInter) and slope (ExtrSlope) for EXTRapolation line

      do i=1, NSTATE
         do j=1, NSTATE
           ExtrSlope(i,j) = (1/DT) * (D12matrix(i,j) - D32matrix(i,j))
           ExtrInter(i,j) = D12matrix(i,j)
         end do
           VenergyInter(i)=VenergyP(i)
           VenergySlope(i)=(Venergy(i)-VenergyP(i))/DT
      end do

      END IF

C   UNCOMMENT to print coefficients !!!
C     just few coefficients
      write(6,*) 'WaveFunctionsCoefficients are: ',NCI,'*',NSTATE
      write(6,*) '       This step        Previous step:        PP:'

      write(6,*) CIBigArray(1), CIBigArrayP(1), CIBigArrayPP(1)
      write(6,*) CIBigArray(2), CIBigArrayP(2), CIBigArrayPP(2)
      write(6,*) CIBigArray(3), CIBigArrayP(3), CIBigArrayPP(3)
      write(6,*) CIBigArray(4), CIBigArrayP(4), CIBigArrayPP(4)

C    All coefficients
C      do i=1,NCI*NSTATE
C      write(6,*) CIBigArray(i), CIBigArrayP(i), CIBigArrayPP(i)
C      END DO
C   UNCOMMENT to print coefficients !!!

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                        GET RANDOM number LO                          C
C                                                                      C


C     Check if Random number is fixed by user
      if (fixedrandL) then
        write(6,*) 'From input Random is Fixed at:', FixedRand
        LO=FixedRand/dble(NSUBSTEPS)
        write(6,*) 'Normalized by Substep number:', LO
      else
C       check if seed number is read from input / RunFile
        if (iseedL) then
          call qpg_iscalar('Seed',Found)
          if (Found) then
            call get_iscalar('Seed',iseed)
            write(6,*) 'Seed number read from the RunFile: ', iseed
          else
            iseed = InitSeed ! initial seed read from input
            write(6,*) 'Seed number read from input file: ', iseed
          end if
C       or generate a new random seed number
        else
          call date_and_time( date, time, zone, values )
C         Just milliseconds multiplied by seconds
          iseed = ((values(7)+1)*values(8)+1)
        endif
        LO = Random_Molcas(iseed)
        call put_iscalar('Seed',iseed)
        if (LO.lt.RandThreshold) then
            LO=RandThreshold
        end if
      end if

C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C                             INTEGRATOR                               C
C                   NSUBSTEPS is the substep number                    C
C                                                                      C

      write(6,*) 'executing the integration:'
      write(6,*) 'Time Step is:    ', DT
      write(6,*) 'Substeps are: ',NSUBSTEPS

      temproot=iRlxRoot

      DO ii=1,NSUBSTEPS
      if (fixedrandL) then
         continue
      else
        LO = Random_Molcas(iseed)
        call put_iscalar('Seed',iseed)
        if (LO.lt.RandThreshold) then
           LO=RandThreshold
        end if
      end if

      do i=1,NSTATE
        do j=1,NSTATE
          Dmatrix(i,j) = ExtrInter(i,j) + ExtrSlope(i,j) *
     &                   (((ii-1)*DT/NSUBSTEPS) - DT/2)
          if(i.ne.j) then
             V(i,j) = 0.0
          else
             V(i,j) = VenergyInter(i) +
     &                VenergySlope(i)*(ii-1)*DT/NSUBSTEPS
          end if
        end do
      end do

      do i=1,NSTATE
       do j=1, NSTATE
         AmatrixDT(i,j) = (0.0, 0.0)
         do l=1,NSTATE
               AmatrixDT(i,j)=AmatrixDT(i,j) + Amatrix(i,l)*Dmatrix(l,j)
     &                                       - Amatrix(l,j)*Dmatrix(i,l)
               AmatrixDT(i,j) = AmatrixDT(i,j) +
     &               (0.0,1.0)*Amatrix(i,l)*V(l,j)-
     &               (0.0,1.0)*Amatrix(l,j)*V(i,l)
         end do

       end do
      end do

      do i=1,NSTATE
         do j=1, NSTATE
            Amatrix(i,j) = Amatrix(i,j) + AmatrixDT(i,j)*DT/NSUBSTEPS
         end do
      end do

      do i=1,NSTATE
        do j=1,NSTATE
         if(i.eq.j) then
!         B(i,i) not used
          Bmatrix(i,j)=2*aimag(CONJG(Amatrix(i,j))*V(i,j))
         else
          Bmatrix(i,j)=-2*real(CONJG(Amatrix(i,j))*Dmatrix(i,j))
         endif
        end do
      end do

      do i=1,NSTATE
        if (i.ne.temproot) then
          Gprobab(i) = Bmatrix(i,temproot)*DT/
     &    (real(Amatrix(temproot,temproot))*NSUBSTEPS)
          if (Gprobab(i).lt.0.0) then
             Gprobab(i)=0.0
          end if
         else
          Gprobab(i) = 0.0
        end if
      end do

      SumProb=0.0
      do i=1,NSTATE
       if (i.ne.temproot) then
            SumProb=SumProb + Gprobab(i)
        if (LO.le.SumProb) then
           write(6,*)'Following Tully, should hop from',temproot,'to',i
           if (V(i,i).lt.Etot) then
               write(6,*)'this root has an energy lower than the total',
     &               Etot,' (thus, permitted)'
               Ediffcheck=(V(i,i)-V(temproot,temproot))
               write(6,*)'Ediffcheck is:',Ediffcheck
               Ediffcheck=ABS(Ediffcheck)
               if (Ediffcheck.lt.Ethreshold) then
                  write(6,*)'lower than the threshold:',Ethreshold
                  write(6,*)'temproot set to:', i
                  temproot=i
                  GOTO 9000
               end if
               GOTO 9000
           end if
        GOTO 9000
        end if
       end if
      end do
 9000 continue

C Persico-Granucci all in a THEN branch

      IF (decoherence) THEN
      EKIN=Etot-V(temproot,temproot)
      If (EKIN.le.0.0D0) then
         write(6,*) 'WARNING! Negative Kinetic Energy. Ekin= ',
     &               EKIN,' a.u.'
         write(6,*) 'Kinetic energy rescaled to 10 e-5.'
         EKIN=0.00001D0
      End If
      do i=1,NSTATE
       if(i.ne.temproot) then
        TAU(i)=abs(1.0D0/(V(temproot,temproot)-V(i,i)))*
     &   (1.0D0+DECO/EKIN)
       end if
      end do

      do i=1,NSTATE
       do j=1,NSTATE
        if(i.ne.temproot.and.j.ne.temproot) then
         Amatrix(i,j)=Amatrix(i,j)*exp(-(DT/NSUBSTEPS)/TAU(i))*
     &                 exp(-(DT/NSUBSTEPS)/TAU(j))
        end if
       end do
      end do

      populOS=0.0D0
      do i=1,NSTATE
       if(i.ne.temproot) then
        populOS=populOS+real(Amatrix(i,i))
       end if
      end do

      ArelaxPrev=Amatrix(temproot,temproot)
      Amatrix(temproot,temproot)=1.-populOS

      do i=1,NSTATE
       do j=1,NSTATE
        if(i.ne.temproot.and.j.eq.temproot) then
         Amatrix(i,j)=Amatrix(i,j)*exp(-(DT/NSUBSTEPS)/TAU(i))*
     &     sqrt(Amatrix(temproot,temproot)/ArelaxPrev)
        elseif(i.eq.temproot.and.j.ne.temproot) then
         Amatrix(i,j)=Amatrix(i,j)*exp(-(DT/NSUBSTEPS)/TAU(j))*
     &     sqrt(Amatrix(temproot,temproot)/ArelaxPrev)
        endif
       end do
      end do
      END IF

      do i=1,NSTATE
          Popul(i)=real(Amatrix(i,i))
      end do

      if (tullySubVerb) then
        write(6,*) 'Substep:', ii
        write(6,*) 'Temproot:', temproot
        write(6,*) 'MATRIX D:'
        do i=1, NSTATE
          write (6,*) (Dmatrix(i,j), j=1,NSTATE,1 )
        end do
        write(6,*) 'Density Matrix elements (i=1..#states):'
        do i=1, NSTATE
         write(6,*)( Amatrix(i,j), j=1,NSTATE,1 )
        end do
        write(6,*) 'B Matrix:'
        do i=1, NSTATE
         write(6,*)( Bmatrix(i,j), j=1,NSTATE,1 )
        end do
        write(6,*) 'Probabilities:'
        write(6,*)( Gprobab(j), j=1,NSTATE,1 )
        write(6,*) 'Random Number is: ', LO
        write(6,*) 'Populations Energies:', ( Popul(j), j=1,NSTATE,1 ),
     &           (V(j,j), j=1,NSTATE,1 ), V(temproot,temproot)
        write(6,*) ''
      end if


      END DO

C                                                                      C
C                        END OF INTEGRATOR                             C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      write(6,*) 'Gnuplot:', (Popul(j), j=1,NSTATE,1),
     &      (Venergy(j), j=1,NSTATE,1 ), Venergy(temproot)
C      write(6,*) 'Gnuplot:', ( Popul(j), j=1,NSTATE,1 ),
C     &           (V(j,j), j=1,NSTATE,1 ), V(temproot,temproot)

      if (temproot.eq.iRlxRoot) then
        HOPPED=.false.
      else
C       Is the keyword MAXHOP set?
        CALL Qpg_iScalar('MaxHopsTully',lmaxHop)
        IF (lmaxHop) THEN
           CALL get_iScalar('MaxHopsTully',maxHop)
           CALL qpg_iScalar('Number of Hops',lnHop)
           if (.not.lnHop) then
              nHop=0
           else
              call get_iScalar('Number of Hops',nHop)
           end if
           write (6,*) 'User set a max of', maxHop
           write (6,*) 'nhop is:', nHop
           if (maxHop.le.nHop) then
              write (6,*) 'This surface HOP is not allowed'
              HOPPED=.false.
           else
              nHop = nHop + 1
              call put_iScalar('Number of Hops',nHop)
              HOPPED=.true.
           end if
        ELSE
C In case temproot is different then iRlxRoot, and no lmaxHop has been
C set
           HOPPED=.true.
        END IF
      end if

C give the "HOP TO:" state name ISTATE2 !!!!
C start the HOPPING procedure
      if (HOPPED) then
      ISTATE2=temproot
      write(6,*)'HOP ALLOWED'
      write(6,'(6X,120A1)')'+',('-',i=1,118),'+'
      write(6,'(6X,A,118X,A)')'|','|'
      write(6,'(6X,A,2(47X,A))')'|','A HOP event is detected!',
     &           '|'
      write(6,'(6X,A,118X,A)')'|','|'
      write(6,'(6X,A,44X,2(A,I3,4X),40X,A)')'|','From state:',
     &           iRlxRoot,'To state:',ISTATE2,'|'
      write(6,'(6X,A,118X,A)')'|','|'
      write(6,'(6X,120A1,//)')'+',('-',i=1,118),'+'

      call Put_iScalar('NumGradRoot',ISTATE2)
      call Put_iScalar('Relax CASSCF root',ISTATE2)
      call Put_iScalar('Relax Original root',ISTATE2)
      call Put_dScalar('Last energy',Venergy(ISTATE2))
      end if
C scale velocities
C
C        CALL get_dArray('Velocities',vel,nsAtom*3)
C
C        write(6,*) 'Velocities before Hop:'
C        do i=1, nsAtom
C          write(6,*) vel(i*3-2), vel(i*3-1), vel(i*3)
C        end do
C        EKIN=Etot-Venergy(iRlxRoot)
C        EKIN_target=Etot-Venergy(temproot)
C        scalfac = sqrt(Ekin_target / Ekin)
C        write(6,*) Etot, Venergy(iRlxRoot), Venergy(temproot), EKIN,
C     &             EKIN_target, scalfac
C        do i=1, nsAtom
C           do j=1,3
C              vel(3*(i-1)+j)=scalfac*vel(3*(i-1)+j)
C           end do
C        end do
C        write(6,*) 'Velocities after Hop:'
C        do i=1, nsAtom
C          write(6,*) vel(i*3-2), vel(i*3-1), vel(i*3)
C        end do
C
C        CALL put_dArray('Velocities',vel,nsAtom*3)
C

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C                             SAVING                                   C
C                                                                      C

      call put_lscalar('hopped',HOPPED)
      call Put_dArray('VenergyP', Venergy, NSTATE)
      call Put_dArray('AllCIPP', CIBigArrayP, NCI*NSTATE)
      call Put_dArray('AllCIP', CIBigArray, NCI*NSTATE)
      call put_zarray('AmatrixV',Amatrix,NSTATE*NSTATE)

C                                                                      C
C                           END SAVING                                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      RETURN

      END

