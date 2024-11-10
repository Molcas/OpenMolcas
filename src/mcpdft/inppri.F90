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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!               2024, Matthew R. Hennefarth                            *
!***********************************************************************
subroutine InpPri_m()
  use stdalloc,only:mma_allocate,mma_deallocate
  use OneDat,only:sNoOri
  use constants,only:zero,two
  Use Functionals,only:Print_Info
  Use KSDFT_Info,only:CoefR,CoefX
  use printlevel,only:silent,terse,usual,verbose
  use mcpdft_output,only:iPrLoc
  use Fock_util_global,only:docholesky
  use rctfld_module,only:lRF
  use mcpdft_input,only:mcpdft_options
  use definitions,only:iwp,wp,u6
  use rasscf_global,only:NAC,NFR,NIN,NONEQ,NROOTS,NSEC,Tot_Charge,tot_el_charge, &
                          tot_nuc_charge,header
  use general_data,only:nfro,nish,ndel,nbas,nash,nrs1,nrs2,nrs3,ispin,nactel,nconf,nelec3,nhole1,nsym,ntot1,stsym,nssh
  implicit none

  Character(len=8) :: Fmt1,Fmt2,Label
  Character(len=120) :: Line
  Character(len=3),dimension(8) :: lIrrep

  integer(kind=iwp) :: i,icharge,icomp,iOpt,iPrLev
  integer(kind=iwp) :: iRc,iSyLbl,iSym,left
  integer(kind=iwp) :: lPaper

  real(kind=wp),allocatable :: Tmp0(:)

  IPRLEV = IPRLOC(1)

  ! This should not be done in this function
  Call Put_dScalar('DFT exch coeff',CoefX)
  Call Put_dScalar('DFT corr coeff',CoefR)

  if(mcpdft_options%extparam) then
    call CheckFuncParam(mcpdft_options%extparamfile)
  endif

  IF(IPRLEV == silent) then
    return
  endif

  ! Define the paper width
  lPaper = 132

  left = (lPaper-len(line))/2
  Write(Fmt1,'(A,I3.3,A)') '(',left,'X,A)'
  Write(Fmt2,'(A,I3.3,A)') '(',left,'X,'

  if(iPrLev >= verbose) then
    ! Print the ONEINT file identifier
    write(u6,*)
    write(u6,Fmt1) 'Header of the ONEINT file:'
    write(u6,Fmt1) '--------------------------'
    Write(Line,'(36A2)')(Header(i),i=1,36)
    write(u6,Fmt1) trim(adjustl(Line))
    Write(Line,'(36A2)')(Header(i),i=37,72)
    write(u6,Fmt1) trim(adjustl(Line))
    write(u6,*)

    ! Print cartesian coordinates of the system
    Call PrCoor

  endif

  if(iPrLev >= usual) then
    ! Print orbital and wavefunction specifications
    write(u6,*)
    Line = ' '
    Write(Line(left-2:),'(A)') 'Wave function specifications:'
    Call CollapseOutput(1,Line)
    write(u6,Fmt1) '-----------------------------'
    write(u6,*)
    If(NFR > 0) then
      write(u6,Fmt2//'A,T45,I6)') 'Number of frozen shell electrons',2*NFR
    endif
    write(u6,Fmt2//'A,T45,I6)') 'Number of closed shell electrons',2*NIN
    write(u6,Fmt2//'A,T45,I6)') 'Number of electrons in active shells',NACTEL
    write(u6,Fmt2//'A,T45,I6)') 'Max number of holes in RAS1 space',NHOLE1
    write(u6,Fmt2//'A,T45,I6)') 'Max nr of electrons in RAS3 space',NELEC3

    If(NFR > 0) then
      write(u6,Fmt2//'A,T45,I6)') 'Number of frozen orbitals',NFR
    endif
    write(u6,Fmt2//'A,T45,I6)') 'Number of inactive orbitals',NIN
    write(u6,Fmt2//'A,T45,I6)') 'Number of active orbitals',NAC
    write(u6,Fmt2//'A,T45,I6)') 'Number of secondary orbitals',NSEC
    write(u6,Fmt2//'A,T45,F6.1)') 'Spin quantum number',(DBLE(ISPIN-1))/two
    write(u6,Fmt2//'A,T45,I6)') 'State symmetry',STSYM
    write(u6,fmt2//'A,T40,I11)') 'Number of CSFs',NCONF
    write(u6,Fmt2//'A,T45,I6)') 'Number of RASSCF root(s) available',nroots
    Call CollapseOutput(0,'Wave function specifications:')

    Call Get_cArray('Irreps',lIrrep,24)
    Do iSym = 1,nSym
      lIrrep(iSym) = adjustr(lIrrep(iSym))
    EndDo

    write(u6,*)
    Line = ' '
    Write(Line(left-2:),'(A)') 'Orbital specifications:'
    Call CollapseOutput(1,Line)
    write(u6,Fmt1) '-----------------------'
    write(u6,*)
    write(u6,Fmt2//'A,T47,8I4)') 'Symmetry species',(iSym,iSym=1,nSym)
    write(u6,Fmt2//'A,T47,8(1X,A))') '                ',(lIrrep(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T47,8I4)') 'Frozen orbitals',(nFro(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T47,8I4)') 'Inactive orbitals',(nIsh(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T47,8I4)') 'Active orbitals',(nAsh(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T47,8I4)') 'RAS1 orbitals',(nRs1(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T47,8I4)') 'RAS2 orbitals',(nRs2(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T47,8I4)') 'RAS3 orbitals',(nRs3(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T47,8I4)') 'Secondary orbitals',(nSsh(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T47,8I4)') 'Deleted orbitals',(nDel(iSym),iSym=1,nSym)
    write(u6,Fmt2//'A,T47,8I4)') 'Number of basis functions',(nBas(iSym),iSym=1,nSym)
    Call CollapseOutput(0,'Orbital specifications:')

  endif

  ! Reaction Field Specification
  If(lRF) then
    iRc = -1
    iOpt = ibset(0,sNoOri)
    iComp = 1
    iSyLbl = 1
    Label = 'Mltpl  0'

    call mma_allocate(Tmp0,nTot1+4,Label="Ovrlp")
    Call RdOne(iRc,iOpt,Label,iComp,Tmp0,iSyLbl)
    If(iRc /= 0) then
      write(u6,*) 'InpPri: iRc from Call RdOne not 0'
      write(u6,*) 'Label = ',Label
      write(u6,*) 'iRc = ',iRc
      Call Abend
    Endif
    tot_nuc_charge = Tmp0(nTot1+4)
    call mma_deallocate(Tmp0)
    Tot_El_Charge = Zero
    Do iSym = 1,nSym
      Tot_El_Charge = Tot_El_Charge-two*DBLE(nFro(iSym)+nIsh(iSym))
    EndDo
    Tot_El_Charge = Tot_El_Charge-DBLE(nActEl)
    Tot_Charge = Tot_Nuc_Charge+Tot_El_Charge
    iCharge = Int(Tot_Charge)
    Call PrRF(.False.,NonEq,iCharge,2)
  EndIf

  if(iPrLev >= terse) then
    write(u6,*)
    Line = ' '
    Write(Line(left-2:),'(A)') 'MCPDFT specifications:'
    Call CollapseOutput(1,Line)
    write(u6,Fmt1) '----------------------'
    write(u6,*)
    if(DoCholesky) then
      write(u6,Fmt2//'A,T50,A)') 'Cholesky decomposition','On'
    else
      write(u6,Fmt2//'A,T50,A)') 'Cholesky decomposition','Off'
    endif
    if(mcpdft_options%mspdft) then
      write(u6,Fmt2//'A,T50,A)') 'Type of calculation','MS-PDFT'
    else
      write(u6,Fmt2//'A,T50,A)') 'Type of calculation','MC-PDFT'
    endif
    write(u6,Fmt2//'A,T50,A)') 'On-Top Functional',trim(mcpdft_options%otfnal%otxc)
    write(u6,Fmt2//'A,T45,F9.2)') 'Exchange scaling factor',CoefX
    write(u6,Fmt2//'A,T45,F9.2)') 'Correlation scaling factor',CoefR
    write(u6,Fmt2//'A,T45,F9.2)') 'Wave function energy weight',mcpdft_options%otfnal%lambda
    if(mcpdft_options%wjob) then
      write(u6,Fmt1) 'Final energies (and CI vectors) will be written to wave function file'
    endif
    If(mcpdft_options%grad) then
      write(u6,Fmt1) 'On-top potentials are computed'
      if(mcpdft_options%nac) then
        write(u6,fmt2//'A,T45,I6,1X,"/",1X,I6)') 'MSPDFT states for NAC',mcpdft_options%nac_states(1),mcpdft_options%nac_states(2)
      else
        write(u6,Fmt2//'A,T45,I6)') 'Root chosen for geometry opt.',mcpdft_options%rlxroot
      endif
    endif
    Call CollapseOutput(0,'MCPDFT specifications:')

    ! Print out grid information
    Call Funi_Print()
  endif

  if(iPrLev >= usual) then
    ! Print our DFT functional specifications
    write(u6,*)
    Line = ' '
    Write(Line(left-2:),'(A)') 'DFT functional specifications:'
    call CollapseOutput(1,Line)
    write(u6,Fmt1) '------------------------------'
    Call libxc_version()
    Call Print_Info()
    call CollapseOutput(0,'DFT functional specifications:')
  endif

endsubroutine
