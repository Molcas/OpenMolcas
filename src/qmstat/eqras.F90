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

subroutine EqRas(iQ_Atoms,nAtomsCC,Coord,nBas,nBas_C,nCnC_C,iBigForDeAll,nSizeBig,ip_UNCLE_MOE,nB)

use Index_Functions, only: nTri_Elem
use Constants, only: Zero, One, Three, Four, Ten, Half, Pi, Angstrom, atmToau, auTokJ, deg2rad, KBoltzmann
use Definitions, only: wp, iwp, u6

implicit none
#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"
#include "qmcom.fh"
#include "qm2.fh"
#include "WrkSpc.fh"
#include "warnings.h"
integer(kind=iwp) :: iQ_Atoms, nAtomsCC, nBas(1), nBas_C(1), nCnC_C(MxBasC), iBigForDeAll, nSizeBig, ip_UNCLE_MOE, nB
real(kind=wp) :: Coord(MxAt*3)
integer(kind=iwp) :: i, i9, iAcc, iAt, iCi, iCNum, iCStart, iDisk, iDiskSa, iDist, iDistIm, iDT(3), iDum(1), iFi(3), &
                     iFil(MxQCen,10), iFP(3), iGP(3), iHowMSampIN, iHowMSampUT, ijhr, iLuExtr, iMacro, iMicro, indma, Inte, & !IFG
                     ip_ExpCento, ip_ExpVal, ipAOSum, iProdMax, iSnurr, iSTC, it1h, it1m, it2h, it2m, it3h, it3m, it4h, it4m, &
                     iTriBasQ, iTriMaxBasQ, iTriState, j, jjhr, jMacro, k, nBaseC, nBaseQ, nClas, NCountField, ncParm, ncpart, &
                     nSize, nSizeIm, nVarv
real(kind=wp) :: AddRep, AverFact, Ax, Ay, Az, BetaBol, BoMaH(MxAt), BoMaO(MxAt), Cpu1, Cpu2, Cpu3, Cpu4, Cpu5, dele, DiFac, Dum, & !IFG
                 E2Die, E_Nuc_Part, E_Nuc_Rubbet, Edisp, EEDisp, EHam, Eint(MxQCen,10), Eint_Nuc(MxAt), Elene, EnCLAS, Energy, & !IFG
                 Eold, Esav, Etot, Exrep, Expe, Expran, ExDie, Gam, Gmma, PertElcInt(nTri_Elem(MxBas)), Poli(MxQCen,10), Pr, Ract, & !IFG
                 Rold, RRRnVarv, s90um, Sum1, Smat(MxStOT), SmatPure(MxStOT), SumElcPot(MxQCen,10), t1s, t2s, t3s, t4s, Tim1, & !IFG
                 Tim2, Tim3, timeCLAS, timeEL, timeEX, timeMC, Vmat(MxStOT) !IFG
logical(kind=iwp) :: CalledBefore, DidWeAccept, Exists, Haveri, InCutOff, Loop, SampleThis, Skip
character(len=200) :: Head
character(len=20) :: MemLaaaa, Memlaaab, Memlaabe, Memlabel, MemQFal
character(len=4) :: Labjhr
character(len=2) :: ChCo, ChCo2
real(kind=wp), parameter :: BoltzK = 1.0e-3_wp*KBoltzmann/auTokJ, &
                            ExLim = Ten !Over how long distance the exchange rep. is computed, the solv-solv.
integer(kind=iwp), external :: IsFreeUnit
real(kind=wp), external :: Ranf
!****Jose** Interaction with Slater type to consider Penetration
!           Eint_Nuc
!****JoseMEP**New variables for the MEP calculation
!             SumElcPot, PertElcInt, Labjhr

! Enter eqras.

! Numbers, initializations, conversions.

Ract = Rstart         !Initial radie
delX = delX/Angstrom  !angstrom-->bohr
delFi = delFi*deg2rad !degree-->radian
delR = delR/Angstrom  !angstrom-->bohr
CalledBefore = .false.
SampleThis = .true.
iBigForDeAll = iBigT
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
iTriBasQ = nBaseQ*(nBaseQ+1)/2
iTriState = nState*(nState+1)/2
timeCLAS = Zero
timeEX = Zero
timeEL = Zero
timeMC = Zero

! If some damping has been requested, prepare it here and print.

if ((Dispdamp .or. FieldDamp) .and. (iPrint >= 8)) then
  write(u6,*)
  write(u6,*) '-----Various damping data.'
  write(u6,*)
end if
if (Dispdamp) then

  ! Construct the Born-Mayer parameters, a la Brdarski-Karlstrom

  call BornMayerBK(iQ_Atoms,BoMaH,BoMaO)
end if

! Damping of field.

if (Fielddamp) then
  if (iPrint >= 8) then
    write(u6,*) '  Damping the field between Qm-region and solvent.'
    write(u6,*) '  E_damp=E_0*(1-exp(alpha*distance))^N'
    write(u6,*) '  alpha(QM-oxygen)   =',CAFieldG
    write(u6,*) '  alpha(QM-hydrogen) =',CBFieldG
    write(u6,*) '  N                  =',CFexp
  end if
end if

! Check what type of simulation to run, and generate some output of
! utmost beauty.

if (Qmeq .and. (iRead /= 9)) then
  call NiceOutPut('EIQ',Gam,Gmma,BetaBol)
else if (QmProd .and. (iRead /= 9)) then
  call NiceOutPut('PIQ',Gam,Gmma,BetaBol)
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
  call NiceOutPut('SSS',Gam,Gmma,BetaBol)
  call DaName(iLuSaIn,SaFilIn)
  iDiskSa = 0
  call iDaFile(iLuSaIn,2,iDum,1,iDiskSa)
  iHowMSampIN = iDum(1)
  iLuExtr = 54
  iLuExtr = IsFreeUnit(iLuExtr)
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
    do ijhr=1,MxQCen
      do jjhr=1,10
        SumElcPot(ijhr,jjhr) = Zero
        AvElcPot(ijhr,jjhr) = Zero
      end do
    end do
    NCountField = 0

    iTriMaxBasQ = MxBas*(MxBas+1)/2
    call dcopy_(iTriMaxBasQ,[Zero],0,PertNElcInt,1)
    call GetMem('SumOvlAOQ','Allo','Real',ipAOSum,iTriBasQ)
    call dcopy_(iTriBasQ,[Zero],0,Work(ipAOSum),1)
  end if
  !*********
else
  write(u6,*)
  write(u6,*) 'An invalid number of iRead detected.'
  call Quit(_RC_INTERNAL_ERROR_)
end if
!----------------------------------------------------------------------*
! If we have input file, then read from it.                            *
!----------------------------------------------------------------------*
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
    if (iExtra > 0) then
      call NyPart(iExtra,nPart,Cordst,Rstart,nCent,iSeed)
    end if
    if (iPrint >= 10) then
      write(Head,*) 'Coordinates of the initial distribution.'
      call Cooout(Head,Cordst,nPart,nCent)
    end if
  end if

  ! Give a startvalue for the Total energy. The effect is that we
  ! always accept the first microstep.

  Etot = 1.0e10_wp

  ! Some numbers.

  ncpart = Ncent*nPart
  ncParm = ncPart-(nCent*icNum)
  nClas = nPart-iCNum
  indma = npart*npol
  iCi = (iQ_Atoms*(iQ_Atoms+1))/2

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

      call ClasClas(iCNum,iCStart,ncParm,Coord,iFP,iGP,iDT,iFI,iDist,iDistIm,Elene,Edisp,Exrep,E2Die,ExDie)
      call QMPosition(EHam,Cordst,Coord,Forcek,dLJrep,Ract,iQ_Atoms)
      call Timing(Cpu2,Tim1,Tim2,Tim3)
      timeCLAS = timeCLAS+(Cpu2-Cpu1)
      !----------------------------------------------------------------*
      ! Work a bit with the quantum part.                              *
      !----------------------------------------------------------------*
      do i=1,3
        xyzMyQ(i) = Zero !Dipoles for the QM-part, see polink.
        xyzMyI(i) = Zero
      end do
      nSize = 3*nPol*nPart
      do i=1,iCi !Allocate memory for the field on the QM-mol. iCi: number of quantum molecule sites.
        do j=1,10
          write(MemQFal,'(A,i2.2,i2.2)') 'Falt',i,j
          call GetMem(MemQFal,'Allo','Real',iFil(i,j),nSize)
        end do
      end do
      do i=1,iCi
        do j=1,10 !Charges (1),Dipoles(3),Quadrupoles(6)
          do k=1,nPart*nPol !Classical polarisation sites including quantum molecule.
            Work(iFil(i,j)-1+k) = Zero
            Work(iFil(i,j)-1+k+nPart*nPol) = Zero
            Work(iFil(i,j)-1+k+2*nPart*nPol) = Zero
          end do
          Eint(i,j) = Zero
        end do
        if (i <= MxAt) Eint_Nuc(i) = Zero
      end do

      ! Compute the exchange operator.

      call ExRas(iCStart,nBaseQ,nBaseC,nCnC_C,iQ_Atoms,nAtomsCC,Ax,Ay,Az,iTriState,Smat,SmatPure,InCutOff,ipAOSum)
      call Timing(Cpu3,Tim1,Tim2,Tim3)
      timeEX = timeEX+(Cpu3-Cpu2)

      ! Electrostatics commencing.

      ! Compute various gradients of 1/r.

      if (lSlater) then
        call OneOverR_Sl(iFil,Ax,Ay,Az,BoMaH,BoMaO,EEDisp,iCNum,Eint,iQ_Atoms,outxyzRAS,Eint_Nuc)
      else
        call OneOverR(iFil,Ax,Ay,Az,BoMaH,BoMaO,EEDisp,iCNum,Eint,iQ_Atoms,outxyzRAS)
      end if

      ! Couple the point-charges in the solvent to the QM-region.

      call HelState(Eint,nState,iCi,RasCha,RasDip,RasQua,Vmat,iPrint)

      ! Let QM-region and solvent polarize.
      call PolRas(iDist,iDistIM,iDT,iFI,iFP,iFil,iCStart,iTriState,VMat,Smat,DiFac,Ract,iCNum,Energy,nVarv,iSTC,Haveri,iQ_Atoms, &
                  ip_ExpVal,Poli)

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

      call BoostRep(AddRep,SmatPure,iSTC,nState,InCutOff)

      ! Sum-up what we will call QM-region energy.

      Energy = Energy-EEdisp+AddRep

      !----------------------------------------------------------------*
      ! Final induction and reaction field energies.                   *
      !----------------------------------------------------------------*
      call ReaInd(iGP,iDT,iDistIm,iCNum,IndMa,NcParm,Sum1,s90um)
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
        if (Haveri) Etot = 999999
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
        if (lExtr(7)) call AllenGinsberg('RASSI',Eint,Poli,ChaNuc,RasCha,RasDip,RasQua,MxStOT,iSTC,nState,iExtr_Atm,lExtr(4), &
                                         iExtr_Eig,iQ_Atoms,ip_ExpCento,E_Nuc_Part,lSlater,Eint_Nuc)

        call Extract(iLuExtr,i9,Etot,xyzMyQ,HMatState,iSTC,iDt,nState,HMatSOld,xyzQuQ,ip_ExpVal,ip_ExpCento,E_Nuc_Rubbet,E_Nuc_Part)
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
          Expran = ranf(iseed)
          if (iPrint >= 10) then
            write(u6,*) '         Positive energy change!'
            write(u6,*) '         Boltzmann weight:',Expe
            write(u6,*) '         Random number:',ExpRan
          end if
          iAcc = iAcc+1
          if (Expe < ExpRan) then
            call Oldge(iAcc,Etot,Eold,Ract,Rold)
            DidWeAccept = .false.
            if (iPrint >= 10) then
              write(u6,*) '         Not accepted!'
            end if
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
              if (Inte == ((iMacro-1)*nMicro+iMicro)) then
                call Put9(Etot,Ract,iDT,iHowMSampUT,Gmma,Gam,Esav,iDisk)
              end if
            end if
          end if
        end if
      end if

      ! Free memory.

      nSize = (nClas*(nClas-1)/2)*(nCent**2)
      nSizeIm = (nClas*nCent)**2
      call GetMem('DistMat','Free','Real',iDist,nSize)
      call GetMem('DistMatIm','Free','Real',iDistIm,nSizeIm)
      do i=1,3
        write(ChCo,'(I1.1)') i
        write(MemLabel,*) 'FP'//ChCo
        write(MemLaabe,*) 'GP'//ChCo
        write(MemLaaab,*) 'DT'//ChCo
        write(MemLaaaa,*) 'FI'//ChCo
        call GetMem(MemLabel,'Free','Real',iFP(i),IndMa)
        call GetMem(MemLaabe,'Free','Real',iGP(i),IndMa)
        call GetMem(MemLaaab,'Free','Real',iDT(i),IndMa)
        call GetMem(MemLaaaa,'Free','Real',iFI(i),IndMa)
      end do
      nSize = 3*nPol*nPart
      do i=1,iCi
        do j=1,10
          write(ChCo,'(I2.2)') i
          write(ChCo2,'(I2.2)') j
          write(MemQFal,*) 'Falt'//ChCo//ChCo2
          call GetMem(MemQFal,'Free','Real',iFil(i,j),nSize)
        end do
      end do
      call GetMem('Coeff','Free','Real',iSTC,nState**2)
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

          AverFact = One/real(NCountField,kind=wp)
          call DaxPy_(iTriBasQ,AverFact,Work(ipAOSum),1,PertNElcInt,1)
          call GetMem('SumOvlAOQ','Free','Real',ipAOSum,iTriBasQ)
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
if (MoAveRed) then
  nB = nRedMO
  ip_UNCLE_MOE = ipAvRed
else
  nB = nBaseQ
end if
nSizeBig = nState*(nState+1)*nB*(nB+1)/4

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
