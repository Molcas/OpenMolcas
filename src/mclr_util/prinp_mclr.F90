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

subroutine PrInp_MCLR(iPL)
!***********************************************************************
!                                                                      *
!     Echo input                                                       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use MCLR_Data, only: ChDisp, DspVec, IRLXROOT, ISNAC, ISTATE, lDisp, NACSTATES, nexp_max, NSSA, SA, SwLbl, XISPSM
use input_mclr, only: AtLbl, ChIrr, Coor, Eps, ERASSCF, ESCF, Header1I, iMCPD, iMethod, iPT2, iRoot, iSpin, mTit, nActEl, nAsh, &
                      nAtoms, nBas, nCSF, nDel, nDisp, nElec3, NewCho, nFro, nHole1, nIsh, nIter, nOrb, nRoots, nRS1, nRS2, nRS3, &
                      nSkip, nSym, ntAsh, ntBas, ntIsh, nTPert, Perturbation, PotNuc, PT2, SpinPol, State_Sym, State_Sym, &
                      StepType, TitleIn, TwoStep, Weight
use PCM_grad, only: RFPERT
use Constants, only: Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iPL
integer(kind=iwp) :: i, iAT, iDisp, iSym, j, jDisp, Left, lLine, nLine, nTSsh
logical(kind=iwp) :: RICD
character(len=100) :: Line
character(len=8) :: Fmt1, Fmt2
character, parameter :: XYZ(3) = ['X','Y','Z']

!----------------------------------------------------------------------*
!     Initialize blank and header lines                                *
!----------------------------------------------------------------------*

lLine = len(Line)
!lPaper = 132
!left = (lPaper-lLine)/2
left = 5
write(Fmt1,'(A,I3.3,A)') '(',left,'X,A)'
write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
!----------------------------------------------------------------------*
call DecideOnCholesky(RICD)
!----------------------------------------------------------------------*
!     Print the project title                                          *
!----------------------------------------------------------------------*
!if (mTit > 0) then
if (iPL >= 3) write(u6,*)
nLine = mTit+5
do i=1,nLine
  Line = ''
  if ((i == 1) .or. (i == nLine)) Line = repeat('*',lLine)
  if (i == 3) Line = 'Project:'
  if ((i >= 4) .and. (i <= nLine-2)) write(Line,'(18A4)') (TitleIN((i-4)*18+j),j=1,18)
  if (iPL >= 3) then
    call Center_Text(Line)
    write(u6,Fmt1) '*'//Line//'*'
  end if
end do
if (iPL >= 3) write(u6,*)
!end if
!----------------------------------------------------------------------*
!     Print file identifiers                                           *
!----------------------------------------------------------------------*
if (iPL >= 3) then
  write(u6,*)
  write(u6,Fmt1) 'Header of the ONEINT file:'
  write(u6,Fmt1) '--------------------------'
  write(Line,Fmt1) Header1I(1)
  write(u6,'(A)') trim(Line)
  write(Line,Fmt1) Header1I(2)
  write(u6,'(A)') trim(Line)
  write(u6,*)
  !--------------------------------------------------------------------*
  !     Print cartesian coordinates of the system                      *
  !--------------------------------------------------------------------*
  write(u6,*)
  write(u6,Fmt1) 'Cartesian coordinates:'
  write(u6,Fmt1) '----------------------'
  write(u6,*)
  write(u6,Fmt1) '----------------------------------------------'
  write(u6,Fmt1) ' No.    Label       X         Y         Z'
  write(u6,Fmt1) '----------------------------------------------'
  do iAt=1,nAtoms
    write(u6,Fmt2//'I3,5X,A6,2X,3F10.5)') iAt,AtLbl(iAt),Coor(1,iAt),Coor(2,iAt),Coor(3,iAt)
  end do
  write(u6,Fmt1) '----------------------------------------------'
  write(u6,Fmt2//'A,F20.10)') 'Nuclear repulsion energy =',PotNuc
end if
!----------------------------------------------------------------------*
!     Print orbital and wavefunction specifications                    *
!----------------------------------------------------------------------*
if (iMethod == 2) then
  ntIsh = sum(nIsh(1:nSym))
  ntAsh = sum(nAsh(1:nSym))
  ntBas = sum(nBas(1:nSym))
  ntSsh = ntBas-ntIsh-ntAsh-sum(nFro(1:nSym)+nDel(1:nSym))
  if (iPL >= 2) then
    write(u6,*)
    Line = ''
    write(Line(left-2:),'(A)') 'Wave function specifications:'
    call CollapseOutput(1,Line)
    write(u6,Fmt1) '-----------------------------'
    write(u6,*)
    write(u6,Fmt2//'A,T47,I6)') 'Number of closed shell electrons',2*ntIsh
    write(u6,Fmt2//'A,T47,I6)') 'Number of electrons in active shells',nActEl
    write(u6,Fmt2//'A,T47,I6)') 'Max number of holes in RAS1 space',nHole1
    write(u6,Fmt2//'A,T47,I6)') 'Max number of electrons in RAS3 space',nElec3
    write(u6,Fmt2//'A,T47,I6)') 'Number of inactive orbitals',ntIsh
    write(u6,Fmt2//'A,T47,I6)') 'Number of active orbitals',ntAsh
    write(u6,Fmt2//'A,T47,I6)') 'Number of secondary orbitals',ntSsh
    write(u6,Fmt2//'A,T47,F6.1)') 'Spin quantum number',real(iSpin-1,kind=wp)*Half
    write(u6,Fmt2//'A,T47,I6)') 'State symmetry',State_Sym
    write(u6,Fmt2//'A,T47,I6)') 'Number of CI roots',nroots
    write(u6,Fmt2//'A,(T47,10I6))') 'States considered',(iroot(i),i=1,nroots)
    write(u6,Fmt2//'A,(T47,10F6.3))') 'Weights',(weight(i),i=1,nroots)
    write(u6,*)
    write(u6,Fmt2//'A,T47,8I6)') 'Symmetry species',(i,i=1,nSym)
    write(u6,Fmt2//'A,T47,8I6)') 'Skipped sym. species',(nSkip(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T47,8I6)') 'Frozen orbitals',(nFro(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T47,8I6)') 'Inactive orbitals',(nIsh(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T47,8I6)') 'Active orbitals',(nAsh(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T47,8I6)') 'RAS1 orbitals',(nRs1(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T47,8I6)') 'RAS2 orbitals',(nRs2(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T47,8I6)') 'RAS3 orbitals',(nRs3(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T47,8I6)') 'Deleted orbitals',(nDel(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T47,8I6)') 'Number of basis functions',(nBas(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T47,8I6)') 'Number of orbitals',(nOrb(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T47,8I6)') 'Number of configurations',(ncsf(isym),isym=1,nsym)

    write(u6,Fmt2//'A,T47,8I6)') 'Number of combinations',(nint(xispsm(isym,1)),isym=1,nsym)

    if (iPt2 == 0) then
      write(u6,Fmt1) 'Natural orbitals are used in the last CI'
    else
      write(u6,Fmt1) 'Pseudo canonical orbitals are used in the last CI'
    end if

    write(u6,Fmt2//'A,T33,F20.10)') 'RASSCF state energy =',ERASSCF(istate)
    write(u6,Fmt2//'A,T47,I6)') 'Size of explicit Hamiltonian in PCG:',nExp_Max
    call CollapseOutput(0,'Wave function specifications:')
  end if
else
  if (iPL >= 2) then
    write(u6,*)
    Line = ''
    write(Line(left-2:),'(A)') 'Wave function specifications:'
    call CollapseOutput(1,Line)
    write(u6,Fmt1) '-----------------------------'
    write(u6,Fmt2//'A,T50,A)') 'Wavefunction type:','SCF'
    write(u6,Fmt2//'A,T45,I6)') 'Number of irreducible symmetry groups',nSym
    write(u6,Fmt2//'A,T45,8I6)') 'Number of basis functions',(nBas(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T45,8I6)') 'Number of occupied orbitals',(nish(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T31,F20.10)') 'SCF energy =',ESCF
    call CollapseOutput(0,'Wave function specifications:')
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPL >= 2) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  write(u6,*)
  write(u6,Fmt2//'A,T42,ES11.4)') 'Convergence threshold=',Eps
  write(u6,Fmt2//'A,T45,I8)') 'Max number of iterations in PCG:',nIter
  if (RICD) then
    if (NewCho) then
      write(u6,Fmt2//'A)') 'Using the Cho-Fock Algorithm'
    else
      write(u6,Fmt2//'A)') 'Using the Cho-MO Algorithm'
    end if
  end if

  if (SPINPOL) then
    write(u6,Fmt1) 'Calculating spin polarization'
  else if (SA .or. iMCPD .or. PT2) then
    if (PT2) write(u6,Fmt2//'A)') 'Calculating Lagrangian multipliers for CASPT2'
    if (isNAC) then
      write(u6,Fmt2//'A,I3,"/",I3)') 'Lagrangian multipliers are calculated for states no. ',NACstates(1),NACstates(2)
      if ((NSSA(1) /= NACstates(1)) .or. (NSSA(2) /= NACstates(2))) &
        write(u6,Fmt2//'39X,A,I3,"/",I3,A)') '(SA roots no. ',NSSA(1),NSSA(2),')'
    else
      write(u6,Fmt2//'A,I3)') 'Lagrangian multipliers are calculated for state no. ',irlxroot
      if (istate /= irlxroot) write(u6,Fmt2//'39X,A,I3,A)') '(SA root no. ',istate,')'
    end if
    if (TwoStep) then
      if (StepType == 'RUN1') write(u6,Fmt1) 'TwoStep activated. Run 1 (preparation).'
      if (StepType == 'RUN2') write(u6,Fmt1) 'TwoStep activated. Run 2 (final run).'
    end if
  else
    if (ndisp /= 0) then
      write(u6,*)
      Line = ''
      write(Line(left-2:),'(A)') 'Perturbation specifications:'
      call CollapseOutput(1,Line)
      write(u6,Fmt1) '----------------------------'
      write(u6,*)
      write(u6,Fmt2//'A,T49,8I4)') 'Number of perturbations in each symmetry',(ldisp(iSym),iSym=1,nSym)
      write(u6,Fmt2//'A,T52,A)') 'Type of perturbation:',Perturbation
      call CollapseOutput(0,'Perturbation specifications:')
      write(u6,*)
      Line = ''
      write(Line(left-2:),'(A)') 'Perturbations:'
      call CollapseOutput(1,Line)
      write(u6,Fmt1) '--------------'
      write(u6,*)
      write(u6,Fmt1) '-------------------------------------'
      write(u6,Fmt1) ' No.    Symmetry    Center Direction'
      write(u6,Fmt1) '-------------------------------------'
      jDisp = 0
      do iSym=1,nSym
        do iDisp=1,lDisp(iSym)
          jDisp = jDisp+1
          if (btest(ntpert(jdisp),4)) then
            write(u6,Fmt2//'I3,T16,A3,T29,A)') jDisp,chIrr(isym),ChDisp(dspvec(jDisp))
          else
            write(u6,Fmt2//'I3,T16,A3,T29,A8,A,A)') jDisp,chIrr(isym),swlbl(jDisp),' ',XYZ(dspvec(jDisp))
          end if
        end do
      end do
      write(u6,Fmt1) '-------------------------------------'
      call CollapseOutput(0,'Perturbations:')
    end if
  end if
  write(u6,*)
  !--------------------------------------------------------------------*
  !     Print reference state information                              *
  !--------------------------------------------------------------------*
  if (.not. SA) then
    write(u6,*)
    if (iMethod == 2) then
      write(u6,Fmt2//'A,I3)') 'Linear response function is computed for root no. = ',irlxroot
    else
      write(u6,Fmt2//'A,I3)') 'Linear response function is computed for Restricted Hartree-Fock wavefunction'
    end if
  end if

  if (RFpert) then
    write(u6,*)
    write(u6,Fmt1) 'Reaction field specifications:'
    write(u6,Fmt1) '------------------------------'
    write(u6,*)
    write(u6,'(6X,A)') 'The Reaction field is added as a perturbation and has been determined in a previous calculation'
    write(u6,*)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
write(u6,*)
!                                                                      *
!***********************************************************************
!                                                                      *
if (isNAC .and. (nSym > 1)) then
  call WarningMessage(2,'NAC is not supported with symmetry')
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine PrInp_MCLR
