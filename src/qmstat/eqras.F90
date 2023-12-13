!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine EqRas(iQ_Atoms,nAtomsCC,Coord,nBas,nBas_C)

use qmstat_global, only: AvElcPot, CAFieldG, CBFieldG, CFexp, ChaNuc, CordIm, Cordst, DelFi, DelR, DelX, Diel, DispDamp, dLJrep, &
                         FieldDamp, Forcek, HmatState, iExtr_Eig, iExtra, iLuSaIn, iLuSaUt, info_atom, Inter, iPrint, iRead, &
                         iSeed, lExtr, lSlater, nAtom, nCent, nMacro, nMicro, nPart, nPol, nState, nTemp, OldGeo, outxyzRAS, &
                         PertNElcInt, Pres, Qmeq, QmProd, ParallelT, RasCha, RasDip, RasQua, rStart, SaFilIn, SaFilUt, SimEx, &
                         SURF, Temp, xyzMyI, xyzMyQ, xyzQuQ
use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Three, Four, Ten, Half, Pi, Angstrom, atmToau, auTokJ, deg2rad, KBoltzmann
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iQ_Atoms, nAtomsCC, nBas(1), nBas_C(1)
real(kind=wp), intent(in) :: Coord(3,iQ_Atoms)
integer(kind=iwp) :: i, i9, iAcc, iAt, iCi, iCNum, iCStart, iDisk, iDiskSa, iDum(1), iHowMSampIN, iHowMSampUT, iLuExtr, iMacro, &
                     iMicro, IndMa, Inte, iProdMax, iSnurr, it1h, it1m, it2h, it2m, it3h, it3m, it4h, it4m, iTriBasQ, iTriState, &
                     jMacro, nBaseC, nBaseQ, nClas, NCountField, nVarv
real(kind=wp) :: AddRep, Ax, Ay, Az, BetaBol, Cpu1, Cpu2, Cpu3, Cpu4, Cpu5, dele, DiFac, Dum, E2Die, E_Nuc_Part, E_Nuc_Rubbet, &
                 Edisp, EEDisp, EHam, Elene, EnCLAS, Energy, Eold, Esav, Etot, Exrep, Expe, Expran, ExDie, Gam, Gmma, &
                 PertElcInt(1), Pr, Ract, Rold, RRRnVarv, s90um, Sum1, t1s, t2s, t3s, t4s, Tim1, Tim2, Tim3, timeCLAS, timeEL, &
                 timeEX, timeMC
logical(kind=iwp) :: CalledBefore, DidWeAccept, Exists, Haveri, InCutOff, Loop, SampleThis, Skip
character(len=200) :: Head
character(len=4) :: Labjhr
real(kind=wp), allocatable :: AOSum(:), BoMaH(:), BoMaO(:), Dist(:,:,:), DistIm(:,:,:,:), DT(:,:), Eint(:,:), Eint_Nuc(:), &
                              ExpCento(:,:), ExpVal(:,:), FI(:,:), Fil(:,:,:,:), FP(:,:), GP(:,:), Poli(:,:), Smat(:), &
                              SmatPure(:), STC(:,:), SumElcPot(:,:), Vmat(:)
real(kind=wp), parameter :: BoltzK = 1.0e-3_wp*KBoltzmann/auTokJ, &
                            ExLim = Ten !Over how long distance the exchange rep. is computed, the solv-solv.
integer(kind=iwp), external :: IsFreeUnit
real(kind=wp), external :: Random_Molcas
#include "warnings.h"
!****Jose** Interaction with Slater type to consider Penetration
!           Eint_Nuc
!****JoseMEP**New variables for the MEP calculation
!             SumElcPot, PertElcInt, Labjhr

! Enter eqras.

! Numbers, initializations, conversions.

Ract = rStart         !Initial radie
delX = delX/Angstrom  !angstrom-->bohr
delFi = delFi*deg2rad !degree-->radian
delR = delR/Angstrom  !angstrom-->bohr
CalledBefore = .false.
SampleThis = .true.
if (ParallelT) nMacro = nTemp*nMacro

! If we have vacuum, then no volume-pressure work is done, nor do we
! have any cavitation free-energy.

if ((abs(Diel-One) > 1.0e-4_wp) .and. (SURF > Zero)) then
  Gmma = Pres*Four/Three*Pi*atmToau
  Gam = 5.973455e-4_wp/74.0_wp*SURF !What are these numbers? Is 74.0 supposed to be the surface tension of water in mN/m?
else
  Gmma = Zero
  Gam = Zero
end if

! More numbers.

if (Temp <= Zero) then
  BetaBol = 1.0e23_wp
else
  BetaBol = One/(Temp*BoltzK)
end if
DiFac = -(Diel-One)/(Diel+One)
Expran = Zero
iHowMSampIN = 0
iHowMSampUT = 0
nBaseC = nBas_C(1)
nBaseQ = nBas(1)
iTriBasQ = nTri_Elem(nBaseQ)
iTriState = nTri_Elem(nState)
iCi = nTri_Elem(iQ_Atoms)
timeCLAS = Zero
timeEX = Zero
timeEL = Zero
timeMC = Zero

! If some damping has been requested, prepare it here and print.

if ((DispDamp .or. FieldDamp) .and. (iPrint >= 8)) then
  write(u6,*)
  write(u6,*) '-----Various damping data.'
  write(u6,*)
end if
call mma_allocate(BoMaH,iQ_Atoms,label='BoMaH')
call mma_allocate(BoMaO,iQ_Atoms,label='BoMaO')

! Construct the Born-Mayer parameters, a la Brdarski-Karlstrom
if (DispDamp) call BornMayerBK(iQ_Atoms,BoMaH,BoMaO)

! Damping of field.

if (FieldDamp .and. (iPrint >= 8)) then
  write(u6,*) '  Damping the field between Qm-region and solvent.'
  write(u6,*) '  E_damp=E_0*(1-exp(alpha*distance))^N'
  write(u6,*) '  alpha(QM-oxygen)   =',CAFieldG
  write(u6,*) '  alpha(QM-hydrogen) =',CBFieldG
  write(u6,*) '  N                  =',CFexp
end if

! Check what type of simulation to run, and generate some output of utmost beauty.

call mma_allocate(SumElcPot,iCi,10,label='SumElcPot')

if (Qmeq .and. (iRead /= 9)) then
  call NiceOutPut('EIQ')
else if (QmProd .and. (iRead /= 9)) then
  call NiceOutPut('PIQ')
  call DaName(iLuSaUt,SaFilUt)
  iDisk = 0 !Put some dummy on the sampfile so we have space for the real number later.
  iDum(1) = iHowMSampUT
  call iDaFile(iLuSaUt,1,iDum,1,iDisk)
  ! Below we make a check for extreme cases. Our algorithm to
  ! select sampling configurations sets this limit.
  iProdMax = (2**30-1)*2
  if ((nMacro*nMicro) >= iProdMax) then
    write(u6,*)
    write(u6,*) 'WARNING! Too large numbers for nMacro and nMicro to run production!'
    write(u6,*) '   Their product must not be greater than 2**31!'
    write(u6,*) '   If you wish to make such large samples, you can run several samplings and collect several sampfiles.'
    call Quit(_RC_INTERNAL_ERROR_)
  end if
else if (iRead == 9) then
  ! If we read from sampfile: open the
  ! sampfile and read how many configurations and open extract file.
  call NiceOutPut('SSS')
  call DaName(iLuSaIn,SaFilIn)
  iDiskSa = 0
  call iDaFile(iLuSaIn,2,iDum,1,iDiskSa)
  iHowMSampIN = iDum(1)
  iLuExtr = IsFreeUnit(54)
  call OpnFl(SimEx,iLuExtr,Exists)
  write(iLuExtr,*) 'Extract-File'
  write(iLuExtr,*)
  ! And put some words in the output
  write(u6,*)
  write(u6,*) '   Total number of sampled configurations:',iHowMSampIN
  write(u6,*) '   Reading from the file ',SaFilIn
  write(u6,*) '   Summarizing data put on ',SimEx
  !*****JoseMEP
  ! If we perform MEP calculation, first we make some zeros and allocate some memory.
  if (lExtr(8)) then
    call mma_allocate(AvElcPot,iCi,10,label='AvElcPot')
    SumElcPot(:,:) = Zero
    AvElcPot(:,:) = Zero
    NCountField = 0

    PertNElcInt(:) = Zero
    call mma_allocate(AOSum,iTriBasQ,label='SumOvlAOQ')
    AOSum(:) = Zero
  end if
  !*********
else
  write(u6,*)
  write(u6,*) 'An invalid number of iRead detected.'
  call Quit(_RC_INTERNAL_ERROR_)
end if
if (.not. allocated(AOSum)) call mma_allocate(AOSum,0,label='SumOvlAOQ')
!----------------------------------------------------------------------*
! If we have input file, then read from it.                            *
!----------------------------------------------------------------------*
call mma_allocate(Fil,nPol*nPart,3,iCi,10,label='Fil')
call mma_allocate(Eint,iCi,10,label='Eint')
call mma_allocate(Eint_Nuc,iQ_Atoms,label='Eint_Nuc')
call mma_allocate(Poli,iCi,10,label='Poli')
call mma_allocate(Smat,iTriState,label='Smat')
call mma_allocate(Vmat,iTriState,label='Vmat')
call mma_allocate(SmatPure,iTriState,label='SmatPure')
call mma_allocate(STC,nState,nState,label='Coeff')
call mma_allocate(ExpVal,4,nState,label='ExpVals')
call mma_allocate(ExpCento,4,nState,label='ExpCento')
if (.not. allocated(CordIm)) call mma_allocate(CordIm,3,nPart*nCent,label='CordIm')
iCStart = (((iQ_Atoms-1)/nAtom)+1)*nCent+1
iCNum = (iCStart-1)/nCent
i9 = 0 !i9 is active if iRead == 9 and we are collecting configurations from the sampfile.
outer: do
  Loop = .false.
  i9 = i9+1
  if ((iRead <= 8) .and. (iRead >= 6)) then
    call Get8(Ract,Dum)
  else if (iRead == 9) then
    call Get9(Ract,Coord,info_atom,iQ_Atoms,iDiskSa)
  else
    if (iExtra > 0) call NyPart(iExtra,nPart,Cordst,rStart,nCent,iSeed)
    if (iPrint >= 10) then
      write(Head,*) 'Coordinates of the initial distribution.'
      call Cooout(Head,Cordst,nPart,nCent)
    end if
  end if

  ! Give a startvalue for the Total energy. The effect is that we
  ! always accept the first microstep.

  Etot = 1.0e10_wp

  ! Some numbers.

  nClas = nPart-iCNum
  IndMa = nPart*nPol

  ! Put QM-molecule in its place.

  if ((iRead == 8) .or. (iRead == 0)) then
    call PlaceIt(Coord,iQ_Atoms,iCNum)
  else if (iRead == 6) then
    call PlaceIt9(Coord,Cordst,info_atom,iQ_Atoms)
    if (iPrint >= 10) then
      write(Head,*) 'CM-centred coordinates after substitution.'
      call Cooout(Head,Cordst,nPart,nCent)
    end if
  end if

  !--------------------------------------------------------------------*
  !                                                                    *
  !------------------------- START SIMULATION -------------------------*
  !                                                                    *
  !--------------------------------------------------------------------*

  iSnurr = 0 !How many steps taken totally.

  ! The Macrosteps.

  do iMacro=1,nMacro
    Esav = Zero
    if (iRead == 9) then
      iAcc = 1
    else
      iAcc = 0
    end if

    ! If we are running parallel tempering, then...
    if (ParallelT) call ParaRoot(Ract,BetaBol,Etot,CalledBefore,SampleThis)

    ! The Microsteps.

    Skip = .false.
    do iMicro=1,nMicro
      call Timing(Cpu1,Tim1,Tim2,Tim3)
      Eold = Etot
      iSnurr = iSnurr+1

      ! Generate new configuration, both solvent and QM-region.

      call GeoGen(Ract,Rold,iCNum,iQ_Atoms)

      ! Compute Solvent-solvent interaction.

      call mma_allocate(FP,3,nPol*nPart,label='FP')
      call mma_allocate(GP,3,nPol*nPart,label='GP')
      call mma_allocate(DT,3,nPol*nPart,label='DT')
      call mma_allocate(FI,3,nPol*nPart,label='FI')
      call mma_allocate(Dist,nCent,nCent,nTri_Elem(nClas-1),label='Dist')
      call mma_allocate(DistIm,nCent,nClas,nCent,nClas,label='DistIm')
      call ClasClas(iCNum,nClas,FP,GP,DT,FI,Dist,DistIm,Elene,Edisp,Exrep,E2Die,ExDie)
      call QMPosition(EHam,Cordst,Coord(:,1),Forcek,dLJrep,Ract,iQ_Atoms)
      call Timing(Cpu2,Tim1,Tim2,Tim3)
      timeCLAS = timeCLAS+(Cpu2-Cpu1)
      !----------------------------------------------------------------*
      ! Work a bit with the quantum part.                              *
      !----------------------------------------------------------------*
      xyzMyQ(:) = Zero !Dipoles for the QM-part, see polink.
      xyzMyI(:) = Zero
      Fil(:,:,:,:) = Zero
      Eint(:,:) = Zero
      Eint_Nuc(:) = Zero

      ! Compute the exchange operator.

      call ExRas(iCStart,nBaseQ,nBaseC,iQ_Atoms,nAtomsCC,Ax,Ay,Az,iTriState,Smat,SmatPure,InCutOff,AOSum)
      call Timing(Cpu3,Tim1,Tim2,Tim3)
      timeEX = timeEX+(Cpu3-Cpu2)

      ! Electrostatics commencing.

      ! Compute various gradients of 1/r.

      if (lSlater) then
        call OneOverR_Sl(Fil,Ax,Ay,Az,BoMaH,BoMaO,EEDisp,iCNum,Eint,iQ_Atoms,outxyzRAS,Eint_Nuc)
      else
        call OneOverR(Fil,Ax,Ay,Az,BoMaH,BoMaO,EEDisp,iCNum,Eint,iQ_Atoms,outxyzRAS)
      end if

      ! Couple the point-charges in the solvent to the QM-region.

      call HelState(Eint,nState,iCi,RasCha,RasDip,RasQua,Vmat)

      ! Let QM-region and solvent polarize.
      call PolRas(Dist,DistIM,DT,FI,FP,Fil,iCStart,iTriState,VMat,Smat,DiFac,Ract,iCNum,Energy,nVarv,STC,Haveri,iQ_Atoms,ExpVal, &
                  Poli)

      ! Energy from QM-nuclei interacting with solvent field.

      if (lSlater) then
        do i=1,iQ_Atoms
          Energy = Energy-Eint_Nuc(i)*ChaNuc(i)
        end do
      else
        do i=1,iQ_Atoms
          Energy = Energy-Eint(i,1)*ChaNuc(i)
        end do
      end if

      ! Some additional boost of short-range repulsion.

      call BoostRep(AddRep,SmatPure,STC,nState,InCutOff)

      ! Sum-up what we will call QM-region energy.

      Energy = Energy-EEdisp+AddRep

      !----------------------------------------------------------------*
      ! Final induction and reaction field energies.                   *
      !----------------------------------------------------------------*
      call ReaInd(GP,DT,DistIm,iCNum,IndMa,nClas,Sum1,s90um)
      call Timing(Cpu4,Tim1,Tim2,Tim3)
      timeEL = timeEL+(Cpu4-Cpu3)
      !----------------------------------------------------------------*
      ! Construct the final energy.                                    *
      !----------------------------------------------------------------*
      EnCLAS = Elene+EHam-Edisp+Exrep+E2Die+ExDie
      Etot = EnCLAS-Half*S90um-Sum1+Gmma*Ract**3+Energy+Gam*Ract**2
      Dele = Etot-Eold
      !----------------------------------------------------------------*
      ! Printing and various if requested.                             *
      !----------------------------------------------------------------*
      if (iPrint >= 10) then
        if (Haveri) Etot = 999999.0_wp
        write(u6,*)
        write(u6,*) '   ----Microstep',iMicro
        write(u6,*) '         Number of iterations:',nVarv
        write(u6,*) '         Total energy:',Etot
        write(u6,*) '         Of which is'
        write(u6,*) '            Pairwise solvent-solvent interaction:',EnCLAS
        write(u6,*) '              Solvent-solvent Electrostatic',Elene
        write(u6,*) '              Harmonic Spring:',EHam
        write(u6,*) '              Solvent-solvent Dispersion:',-Edisp
        write(u6,*) '              Solvent-solvent Exchange:',Exrep
        write(u6,*) '            Energy of induced dipoles in field from explicit solvent:',-Sum1
        write(u6,*) '            Energy of solvent charge distribution in reaction field:',-Half*S90um
        write(u6,*) '            Solvent E-interaction with image:',E2Die
        write(u6,*) '            Solvent Repulsion with boundary:',ExDie
        write(u6,*) '            Solvent-Solute dispersion:',EEdisp
        write(u6,*) '            Energy of QM-region:',Energy
        write(u6,*) '              Higher order overlap exchange pair-term:',AddRep
        write(u6,*) '            Surface tension term:',Gam*Ract**2
        write(u6,*) '            Volume-pressure term:',Gmma*Ract**3
        write(u6,*) '         Previous accepted energy:',Eold
        write(u6,*) '         Difference:',Dele
        write(u6,*) '         Total dipole in QM-region:(',-xyzMyQ(1),',',-xyzMyQ(2),',',-xyzMyQ(3),')'
        write(u6,*) '         Radie:',Ract
        if (Haveri) then
          write(u6,*) '    WARNING! SOME OF THE NUMBERS ABOVE HAVE NO MEANING SINCE THE POLARIZATION DID NOT CONVERGE!!!'
          Skip = .true.
          exit
        end if
      end if

      ! If we are collecting stuff from a sampfile, now is the time to put
      ! data on the extract file. If center-specific expectation values
      ! are requested, call Allen.

      if (iRead == 9) then
        if (lExtr(6)) then
          E_Nuc_Rubbet = Zero
          if (lSlater) then
            do iAt=1,iQ_Atoms
              E_Nuc_Rubbet = E_Nuc_Rubbet-(Eint_Nuc(iAt)+Poli(iAt,1))*ChaNuc(iAt)
            end do
          else
            do iAt=1,iQ_Atoms
              E_Nuc_Rubbet = E_Nuc_Rubbet-(Eint(iAt,1)+Poli(iAt,1))*ChaNuc(iAt)
            end do
          end if
        end if
        if (lExtr(7)) call AllenGinsberg('RASSI',Eint,Poli,ChaNuc,RasCha,RasDip,RasQua,STC,nState,lExtr(4),iExtr_Eig,iQ_Atoms, &
                                         ExpCento,E_Nuc_Part,lSlater,Eint_Nuc)

        call Extract(iLuExtr,i9,Etot,xyzMyQ,HMatState,STC,nState,xyzQuQ,ExpVal,ExpCento,E_Nuc_Rubbet,E_Nuc_Part)
        !***JoseMEP**********
        ! If MEP option. Add electr. potential, field, etc. for all solvent config.
        if (lExtr(8)) then
          Labjhr = 'Add '
          call AverMEP(Labjhr,Eint,Poli,iCi,SumElcPot,NCountField,PertElcInt,1,1,[1],[1],1)
          NCountField = NCountField+1
        end if
        !********
      else
        ! Resume the MC-wrap up.

        Dele = Dele*BetaBol
        DidWeAccept = .true.
        if (Dele < Zero) then
          iAcc = iAcc+1
        else
          Expe = exp(-Dele)
          Expran = Random_Molcas(iSeed)
          if (iPrint >= 10) then
            write(u6,*) '         Positive energy change!'
            write(u6,*) '         Boltzmann weight:',Expe
            write(u6,*) '         Random number:',ExpRan
          end if
          if (Expe >= ExpRan) then
            iAcc = iAcc+1
          else
            Etot = Eold
            Ract = Rold
            Cordst(:,:) = OldGeo
            DidWeAccept = .false.
            if (iPrint >= 10) write(u6,*) '         Not accepted!'
          end if
        end if
        if (DidWeAccept) Esav = Esav+Etot
        !--------------------------------------------------------------*
        ! If this is a production run, then put stuff on the sampfile. *
        !--------------------------------------------------------------*
        if (QmProd .and. (iRead /= 9)) then
          if (SampleThis) then
            if (Inter /= 0) then
              Inte = (iSnurr/Inter)*Inter
              if (Inte == ((iMacro-1)*nMicro+iMicro)) call Put9(Etot,Ract,iHowMSampUT,Gmma,Gam,Esav,iDisk)
            end if
          end if
        end if
      end if

      ! Free memory.

      call mma_deallocate(Dist)
      call mma_deallocate(DistIm)
      call mma_deallocate(FP)
      call mma_deallocate(GP)
      call mma_deallocate(DT)
      call mma_deallocate(FI)
      call Timing(Cpu5,Tim1,Tim2,Tim3)
      timeMC = timeMC+(Cpu5-Cpu4)

      !----------------------------------------------------------------*
      ! End of Microstep.                                              *
      !----------------------------------------------------------------*
    end do
    !------------------------------------------------------------------*
    ! Have we collected all sampled configurations? If no, go up again.*
    ! If yes, then close some files and take a little jump downwards to*
    ! the END!!!                                                       *
    !------------------------------------------------------------------*
    !Jose***************************************
    ! This point is also used to perform the Average of the Potential,
    ! Field and Field gradients to obtain and average Electrostatic
    ! perturbation. The Non-Electr. perturbation is also obtained here
    ! it will be added directly to the One-electron file.
    !******************************************

    if (.not. Skip) then
      if ((i9 < iHowMSampIN) .and. (iRead == 9)) then
        cycle outer
      else if ((i9 >= iHowMSampIN) .and. (iRead == 9)) then
        !*****JoseMEP**
        ! If MEP option. Obtain the mean Potential, Field and Field Gradients.
        ! It is also obtained the average of the Non-Electrostatic perturbation
        if (lExtr(8)) then
          Labjhr = 'Aver'
          call AverMEP(Labjhr,Eint,Poli,iCi,SumElcPot,NCountField,PertElcInt,1,1,[1],[1],1)

          PertNElcInt(:) = PertNElcInt+AOSum/real(NCountField,kind=wp)
          call mma_deallocate(AOSum)
        end if
        !********

        call DaClos(iLuSaIn)
        close(iLuExtr)
        exit outer
      end if
    end if
    !------------------------------------------------------------------*
    ! Write to startfile.                                              *
    !------------------------------------------------------------------*
    ESav = Esav/real(iAcc,kind=wp)
    call Put8(Ract,Etot,Gmma,Gam,ESav)
    if (Haveri) call Quit(_RC_NOT_CONVERGED_)
    !------------------------------------------------------------------*
    ! Print some things here at the end of the macrostep.              *
    !------------------------------------------------------------------*
    if (.not. ParallelT) then
      jMacro = iMacro
    else
      jMacro = 1+(iMacro-1)/nTemp
    end if
    write(u6,'(A,i4)') '---Macrostep ',jMacro
    write(u6,'(A,i4)') '     Number of microsteps:',nMicro
    Pr = 100.0_wp*(real(iAcc,kind=wp)/real(nMicro,kind=wp))
    write(u6,'(A,i4,A,f5.1,A)') '     Number of acceptances:',iAcc,'(',Pr,'%)'
    write(u6,'(A,f12.4)') '     Radie (a.u.):',Ract
    write(u6,'(A,f16.8)') '     Average Energy (a.u.) in Macrostep:',Esav
    write(u6,'(A,3(f12.4))') '     Total dipole in QM-region last microstep (a.u.):',-xyzMyQ(1),-xyzMyQ(2),-xyzMyQ(3)
    write(u6,*)
    call xFlush(u6)
    !------------------------------------------------------------------*
    ! End of Macrostep.                                                *
    !------------------------------------------------------------------*
  end do
  !--------------------------------------------------------------------*
  ! Put some things on info-file. Used to make tests.                  *
  !--------------------------------------------------------------------*
  call Add_Info('Total Energy',[Etot],1,6)
  call Add_Info('Induction of system',[Sum1],1,6)
  call Add_Info('React. field int.',[s90um],1,6)
  call Add_Info('Solv-Solu Disp.',[EEdisp],1,6)
  call Add_Info('QM-region Energy',[Energy],1,6)
  call Add_Info('QM-region dipole',xyzMyQ,3,5)
  RRRnVarv = real(nVarv,kind=wp)
  call Add_Info('Pol.Iterations',[RRRnVarv],1,8)
  if (.not. Loop) exit outer
end do outer
call mma_deallocate(BoMaH)
call mma_deallocate(BoMaO)
call mma_deallocate(Fil)
call mma_deallocate(Eint)
call mma_deallocate(Eint_Nuc)
call mma_deallocate(Poli)
call mma_deallocate(Smat)
call mma_deallocate(Vmat)
call mma_deallocate(SmatPure)
call mma_deallocate(STC)
call mma_deallocate(ExpVal)
call mma_deallocate(ExpCento)
call mma_deallocate(SumElcPot)
if (allocated(AOSum)) call mma_deallocate(AOSum)
!----------------------------------------------------------------------*
! Close some external files.                                           *
!----------------------------------------------------------------------*
if (QmProd .and. (iRead /= 9)) then
  iDisk = 0
  iDum(1) = iHowMSampUT
  call iDaFile(iLuSaUt,1,iDum,1,iDisk)
  call DaClos(iLuSaUt)
end if
!----------------------------------------------------------------------*
! The End... be happy!                                                 *
!----------------------------------------------------------------------*

write(u6,*)
write(u6,*)
write(u6,*) ' Time statistics. (hour:minute:second)'
write(u6,*) ' -----------------------------------------------------'
it1h = int(timeCLAS)/3600
it1m = int(timeCLAS-it1h*3600)/60
t1s = timeCLAS-it1h*3600-it1m*60
write(u6,9) '   Time spent on pair-wise solvent-solvent interactions: ',it1h,':',it1m,':',t1s
it2h = int(timeEX)/3600
it2m = int(timeEX-it2h*3600)/60
t2s = timeEX-it2h*3600-it2m*60
write(u6,9) '   Time spent on solvent-solute overlap calculations: ',it2h,':',it2m,':',t2s
it3h = int(timeEL)/3600
it3m = int(timeEL-it3h*3600)/60
t3s = timeEL-it3h*3600-it3m*60
write(u6,9) '   Time spent on solvent and solute electrostatic interaction: ',it3h,':',it3m,':',t3s
it4h = int(timeMC)/3600
it4m = int(timeMC-it4h*3600)/60
t4s = timeMC-it4h*3600-it4m*60
write(u6,9) '   Time spent on the Metropolis-Monte Carlo decision: ',it4h,':',it4m,':',t4s

return

9 format(A,I4,A,I3,A,F5.2)

end subroutine EqRas
