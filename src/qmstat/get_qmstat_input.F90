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

! Process input to QMSTAT. All input variables are stored in
! qminp.fh which in turn are initialized in qmstat_init.
subroutine Get_Qmstat_Input(iQ_Atoms)

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"
#include "warnings.h"
character*180 Key
character*20 Kword
character*180 Get_Ln
character VecsQue*3
dimension CoTEMP1(3), CoTEMP2(3), CoTEMP3(3), CoTEMP4(3), CoTEMP5(3)
dimension SlFacTemp(6)
external Get_Ln, iClast
logical YesNo(20), Changed

! Say what is done and set all YesNo to false; their purpose is to
! keep track on compulsory keywords and certain keyword combinations.

!write(6,*)
!write(6,*)'Input processed...'
do i=1,20
  YesNo(i) = .false.
end do

! Use some nice routines to collect input.

LuRd = 79
LuRd = IsFreeUnit(LuRd)
call SpoolInp(LuRd)
rewind(LuRd)
call RdNlst(LuRd,'QMSTAT')

! The turning-point in this do-while loop.

1000 continue

! Use Get_Ln to read the lines; it takes care of commented lines.

Key = Get_Ln(LuRd)
Kword = trim(Key)
call UpCase(Kword)

! The keywords and their labels.

if (Kword(1:4) == 'TITL') Go To 101
if (Kword(1:4) == 'SIMU') Go To 103
if (Kword(1:4) == 'THRE') Go To 104
if (Kword(1:4) == 'STEP') Go To 105
if (Kword(1:4) == 'RUN ') Go To 106
if (Kword(1:4) == 'PRIN') Go To 107
if (Kword(1:4) == 'EXTE') Go To 108
if (Kword(1:4) == 'EDIT') Go To 109
if (Kword(1:4) == 'CONF') Go To 110
if (Kword(1:4) == 'QMSU') Go To 111
if (Kword(1:4) == 'SOLV') Go To 112
if (Kword(1:4) == 'RASS') Go To 113
if (Kword(1:4) == 'SCFS') Go To 114
if (Kword(1:4) == 'SING') Go To 115
if (Kword(1:4) == 'ANAL') Go To 116
if (Kword(1:4) == 'EXTR') Go To 117
if (Kword(1:4) == 'END ') Go To 99999

! This code is only reached if an illegal keyword in the
! first tier is encountered.

iChrct = len(Kword)
Last = iCLast(Kword,iChrct)
write(6,*)
write(6,*) 'ERROR!'
write(6,'(1X,A,A)') Kword(1:Last),' is not a valid keyword!'
call Quit(_RC_INPUT_ERROR_)

! <<<TITLe>>>   Read title.

101 continue
Key = Get_Ln(LuRd)
Joblab = trim(Key)
ATitle = .true.
Go To 1000

! <<<SIMUlation parameters>>>   Read a variety of simulation parameters.

103 continue
Key = Get_Ln(LuRd)
Kword = trim(Key)
call UpCase(Kword)
if (Kword(1:4) == 'RADI') then
  ! <<<RADIe>>>
  Key = Get_Ln(LuRd)
  call Get_F1(1,Rstart)
  Go To 103
else if (Kword(1:4) == 'PERM') then
  ! <<<PERMitivity>>>
  Key = Get_Ln(LuRd)
  call Get_F1(1,Diel)
  Go To 103
else if (Kword(1:4) == 'TEMP') then
  ! <<<TEMPerature>>>
  Key = Get_Ln(LuRd)
  call Get_F1(1,Temp)
  Go To 103
else if (Kword(1:4) == 'PRES') then
  ! <<<PRESsure>>>
  Key = Get_Ln(LuRd)
  call Get_F1(1,Pres)
  Go To 103
else if (Kword(1:4) == 'SURF') then
  ! <<<SURFace>>>
  Key = Get_Ln(LuRd)
  call Get_F1(1,Surf)
  Go To 103
else if (Kword(1:4) == 'TRAN') then
  ! <<<TRANslation>>>
  Key = Get_Ln(LuRd)
  call Get_F1(1,DelX)
  Go To 103
else if (Kword(1:4) == 'ROTA') then
  ! <<<ROTAtion>>>
  Key = Get_Ln(LuRd)
  call Get_F1(1,DelFi)
  Go To 103
else if (Kword(1:4) == 'CAVI') then
  ! <<<CAVIty>>>
  Key = Get_Ln(LuRd)
  call Get_F1(1,DelR)
  Go To 103
else if (Kword(1:4) == 'FORC') then
  ! <<<FORCe>>>
  Key = Get_Ln(LuRd)
  call Get_F1(1,Forcek)
  Go To 103
else if (Kword(1:4) == 'BREP') then
  ! <<<BREPulsion>>>
  Key = Get_Ln(LuRd)
  call Get_F1(1,dLJRep)
  Go To 103
else if (Kword(1:4) == 'SEED') then
  ! <<<SEED>>>
  Key = Get_Ln(LuRd)
  call Get_I1(1,iSeed)
  Go To 103
else if (Kword(1:4) == 'PARA') then
  ! <<<PARAlleltemp>>>
  ParallelT = .true.
  Key = Get_Ln(LuRd)
  call Get_I1(1,nTemp)
  Key = Get_Ln(LuRd)
  call Get_I(1,nStFilT,nTemp)
  Key = Get_Ln(LuRd)
  call Get_F(1,ParaTemps,nTemp)
  Go To 103
else if (Kword(1:4) == 'END ') then
  ! <<<END simulation parameters>>>
  Go To 1000
end if
! Here we come if something gets wrong above
write(6,*)
write(6,*) ' Unrecognized keyword in the SIMUlation parameter section:',Kword(1:4)
call Quit(_RC_INPUT_ERROR_)

! <<<THREshold>>>    Get the polarization thresholds.

104 continue
Key = Get_Ln(LuRd)
call Get_F1(1,Pollim)
call Get_F1(2,Enelim)
call Get_I1(3,itMax)
Go To 1000

! <<<STEPs>>>   Specify how many macro- and microsteps.

105 continue
Key = Get_Ln(LuRd)
call Get_I1(1,NMacro)
call Get_I1(2,NMicro)
Go To 1000

! <<<RUN >>>   What type of simulation are we to run?

106 continue
Key = Get_Ln(LuRd)
Kword = trim(Key)
call UpCase(Kword)
if (Kword(1:4) == 'ANAL') Anal = .true.
if (Kword(1:4) == 'QMEQ') Qmeq = .true.
if (Kword(1:4) == 'QMPR') then
  QmProd = .true.
  read(LuRd,*) Inter,iNrUt
  write(SaFilUt(6:6),'(i1.1)') iNrUt
  iLuSaUt = 32+iNrUt
end if
if (Kword(1:2) == 'SM') then
  write(6,*)
  write(6,*) 'No classical simulations are available.'
  call Quit(_RC_INPUT_ERROR_)
end if
YesNo(8) = .true.
Go To 1000

! <<<PRINt>>>   Specify print-level.

107 continue
Key = Get_Ln(LuRd)
call Get_I1(1,iPrint)
Go To 1000

! <<<EXTErnal>>>   External one-electron perturbation
!                  should be added on the hamiltonian.

108 continue
AddExt = .true.
Key = Get_Ln(LuRd)
call Get_I1(1,nExtAddOns)
if (nExtAddOns > MxExtAddOn) then
  write(6,*)
  write(6,*) 'Too many external perturbations asked for.'
  call Quit(_RC_INPUT_ERROR_)
end if
do i=1,nExtAddOns
  read(LuRd,*) ScalExt(i),ExtLabel(i),iCompExt(i)
end do
Go To 1000

! <<<EDITstartfile>>> Section for editing and displaying stuff on given startfile.

109 continue
EdSt = .true.
Key = Get_Ln(LuRd)
Kword = trim(Key)
call UpCase(Kword)
if (Kword(1:4) == 'DELE') then
  ! <<<DELEte>>>   Delete solvent molecules.
  DelOrAdd(1) = .true.
  Key = Get_Ln(LuRd)
  call Get_I1(1,NrStarti)
  call Get_I1(2,NrStartu)
  Key = Get_Ln(LuRd)
  call Get_I1(1,nDel)
  Go To 109
end if
if (Kword(1:4) == 'ADD ') then
  ! <<<ADD >>>    Add solvent molecules.
  DelOrAdd(2) = .true.
  Key = Get_Ln(LuRd)
  call Get_I1(1,NrStarti)
  call Get_I1(2,NrStartu)
  Key = Get_Ln(LuRd)
  call Get_I1(1,nAdd)
  Go To 109
end if
if (Kword(1:4) == 'QMDE') then
  ! <<<QMDElete>>>   Substitute all slots with non-water coordinates with water coordinates.
  DelOrAdd(3) = .true.
  Key = Get_Ln(LuRd)
  call Get_I1(1,NrStarti)
  call Get_I1(2,NrStartu)
  Go To 109
end if
if (Kword(1:4) == 'DUMP') then
  ! <<<DUMP coordinates>>>  Dump coordinates in a way suitable for graphical display.
  DelOrAdd(4) = .true.
  Key = Get_Ln(LuRd)
  call UpCase(Key)
  cDumpForm = Key(1:4)
  Key = Get_Ln(LuRd)
  call Get_I1(1,NrStarti)
  Go To 109
end if
if (Kword(1:4) == 'END ') then
  ! <<<END editstartfile>>>
  write(StFilIn(6:6),'(i1.1)') NrStarti
  write(StFilUt(6:6),'(i1.1)') NrStartu
  Go To 1000
end if
! Here we come if something gets wrong above
write(6,*)
write(6,*) ' Unrecognized keyword in the EDITstartfile section:',Kword(1:4)
call Quit(_RC_INPUT_ERROR_)

! <<<CONFiguration>>>    Where is the initial configuration to be
!                        obtained. Also, if we wish to edit the startfile.

110 continue
Key = Get_Ln(LuRd)
Kword = trim(Key)
call UpCase(Kword)
if (Kword(1:4) == 'ADD ') then
  ! <<<ADD >>>  How many to add at random.
  Key = Get_Ln(LuRd)
  call Get_I1(1,iExtra)
  if (iExtra > MxPut) then
    write(6,*)
    write(6,*) 'The present limit of explicit solvent molecules is',MxPut,'.'
    call Quit(_RC_INPUT_ERROR_)
  end if
  Go To 110
else if (Kword(1:4) == 'FILE') then
  ! <<<FILE>>>  Read configuration/s from a file and put them there.
  Key = Get_Ln(LuRd)
  Kword = trim(Key)
  call UpCase(Kword)
  if (Kword(1:4) == 'STAR') then
    ! <<<STARtfile>>>  Read from startfile.
    Key = Get_Ln(LuRd)
    Kword = trim(Key)
    call UpCase(Kword)
    if (Kword(1:4) == 'SCRA') then
      ! <<<SCRAtch>>> Just put QM as given on RUNFILE.
      Key = Get_Ln(LuRd)
      call Get_I1(1,iNrIn)
      call Get_I1(2,iNrUt)
      iRead = 8
    else if (Kword(1:4) == 'COPY') then
      ! <<<COPY>>> Collect place of QM from startfile. WARNING!
      !            You must use consistent startfile and RUNFILE!
      Key = Get_Ln(LuRd)
      call Get_I1(1,iNrIn)
      call Get_I1(2,iNrUt)
      iRead = 7
    else if (Kword(1:2) == 'CM') then
      ! <<<CM  >>> Put QM in CM of QM-place on startfile.
      Key = Get_Ln(LuRd)
      call Get_I1(1,iNrIn)
      call Get_I1(2,iNrUt)
      iRead = 6
    else
      write(6,*)
      write(6,*) 'Illegal StartFile option.'
      call Quit(_RC_INPUT_ERROR_)
    end if
  else if (Kword(1:4) == 'SAMP') then
    ! <<<SAMPfile>>>  Read configurations from sampfile and collect
    !                 the extracted information in iNrExtr.
    Key = Get_Ln(LuRd)
    call Get_I1(1,iNrIn)
    call Get_I1(2,iNrExtr)
    iRead = 9
    write(SimEx(6:6),'(i1.1)') iNrExtr
    YesNo(9) = .true.
  else
    ! CRASH-BOOM-BANG!
    write(6,*)
    write(6,*) ' Error in CONFiguration section, FILE subsection.'
    call Quit(_RC_INPUT_ERROR_)
  end if
  write(StFilIn(6:6),'(i1.1)') iNrIn
  write(StFilUt(6:6),'(i1.1)') iNrUt
  write(SaFilIn(6:6),'(i1.1)') iNrIn
  iLuStIn = 8+iNrIn
  iLuStUt = 16+iNrUt
  iLuSaIn = 24+iNrIn
  Go To 110
else if (Kword(1:4) == 'INPU') then
  ! <<<INPUt>>>  Signify that the first configuration will be given
  !              explicitly in input. The coordinates are then given
  !              in the solvent section.
  Key = Get_Ln(LuRd)
  call Get_I1(1,iNrUt)
  iLuStUt = 16+iNrUt
  YesNo(5) = .true. !User will give coords in input.
  Go To 110
else if (Kword(1:4) == 'END ') then
  YesNo(2) = .true. !Signify that source of starting conf. is specified.
  Go To 1000
end if
! Here we come if something gets wrong above
write(6,*)
write(6,*) ' Unrecognized keyword in the CONFiguration section:',Kword(1:4)
call Quit(_RC_INPUT_ERROR_)

! <<<QMSUrrounding>>>  Give parameters for the QM-Stat.Mech. interaction.

111 continue
Key = Get_Ln(LuRd)
Kword = trim(Key)
call UpCase(Kword)
if (Kword(1:4) == 'DPAR') then
  ! <<<DPARameters>>>  Dispersion
  do i=1,iQ_Atoms
    Key = Get_Ln(LuRd)
    call Get_F(1,Udisp(i,1),1)
    call Get_F(2,Udisp(i,2),1)
  end do
  Go To 111
else if (Kword(1:4) == 'ELEC') then
  ! <<<ELECtrostatic>>> Electrostatic Slater Numbers
1131 continue
  Key = Get_Ln(LuRd)
  Kword = trim(Key)
  call UpCase(Kword)
  if (Kword(1:4) == 'THRE') then
    ! <<<THREsholds>>>  First is the Cutoff (distance Quantum Site-
    !                   Classical molecule) to evaluate Penetration
    !                   effects. Second, difference between two Slater
    !                   exponents to not be consider the same value.
    Key = Get_Ln(LuRd)
    call Get_F1(1,Cut_Elc)
    call Get_F1(2,DifSlExp)
    Go To 1131
  end if
  if (Kword(1:4) == 'NOPE') then
    ! <<<NOPEnetration>>>  Electrostatic Penetration Not Computed
    lSlater = .false.
    Go To 1131
  end if
  if (Kword(1:4) == 'QUAD') then
    ! <<QUADrupoles>>>  Electrostatic Penetration Computed in quadrupoles.
    lQuad = .true.
    Go To 1131
  end if
  if (Kword(1:4) == 'END ') then
    ! <<<END Electrostatic>>>
    Go To 111
  end if
  write(6,*)
  write(6,*) ' Error in QMSUrrounding section, ELECtrostatic subsection.'
  call Quit(_RC_INPUT_ERROR_)
else if (Kword(1:4) == 'XPAR') then
  ! <<<XPARameters>>>  Exchange repulsion
1112 continue
  Key = Get_Ln(LuRd)
  Kword = trim(Key)
  call UpCase(Kword)
  if (Kword(1:4) == 'S2  ') then
    ! <<<S2  >>>  The S2 parameter
    Key = Get_Ln(LuRd)
    call Get_F1(1,Exrep2)
    Go To 1112
  end if
  if (Kword(1:4) == 'S4  ') then
    ! <<<S4  >>>  The S4 parameter
    Key = Get_Ln(LuRd)
    call Get_F1(1,Exrep4)
    Go To 1112
  end if
  if (Kword(1:4) == 'S6  ') then
    ! <<<S6  >>>  The S6 parameter
    Key = Get_Ln(LuRd)
    call Get_F1(1,Exrep6)
    Go To 1112
  end if
  if (Kword(1:4) == 'S10 ') then
    ! <<<S10 >>>  The S10 parameter
    Key = Get_Ln(LuRd)
    call Get_F1(1,Exrep10)
    Go To 1112
  end if
  if (Kword(1:4) == 'CUTO') then
    ! <<<CUTOff>>> The cut-off radii for repulsion. The first is
    !              outer radius that says EX=0 if R > Cut_Ex1, while
    !              the second is a EX=infinity if R < Cut_Ex2.
    Key = Get_Ln(LuRd)
    call Get_F1(1,Cut_Ex1)
    call Get_F1(2,Cut_Ex2)
    Go To 1112
  end if
  if (Kword(1:4) == 'END ') then
    ! <<<END xparameters>>>
    Go To 111
  end if
  write(6,*)
  write(6,*) ' Error in QMSUrrounding section, XPARameters subsection.'
  call Quit(_RC_INPUT_ERROR_)
else if (Kword(1:4) == 'DAMP') then
  ! <<<DAMPing>>>  Damping parameters for Qm-Surrounding interaction.
1119 continue
  Key = Get_Ln(LuRd)
  Kword = trim(Key)
  call UpCase(Kword)
  if (Kword(1:4) == 'DISP') then
    ! <<<DISPersion>>>  Dispersion damping parameters. This part
    !                   should have a AUTO keyword which collects
    !                   default parameters from the MpProp file.
    Dispdamp = .true.
    Key = Get_Ln(LuRd)
    ! Damping numbers for solvent
    call Get_F(1,CharDi(1),1)
    call Get_F(2,QuaDi(1,1),1)
    call Get_F(3,QuaDi(2,1),1)
    call Get_F(4,QuaDi(3,1),1)
    Key = Get_Ln(LuRd)
    call Get_F(1,CharDi(2),1)
    call Get_F(2,QuaDi(1,2),1)
    call Get_F(3,QuaDi(2,2),1)
    call Get_F(4,QuaDi(3,2),1)
    ! Damping numbers for solute
    do i=1,iQ_Atoms
      Key = Get_Ln(LuRd)
      call Get_F(1,CharDiQ(i),1)
      call Get_F(2,QuaDiQ(1,i),1)
      call Get_F(3,QuaDiQ(2,i),1)
      call Get_F(4,QuaDiQ(3,i),1)
    end do
    Go To 1119
  else if (Kword(1:4) == 'FIEL') then
    ! <<<FIELd>>> Parameters for damping electric field.
    Fielddamp = .true.
    Key = Get_Ln(LuRd)
    call Get_F1(1,CAFieldG)
    call Get_F1(2,CBFieldG)
    call Get_F1(3,CFExp)
    Go To 1119
  else if (Kword(1:4) == 'END ') then
    ! <<<END damping>>>
    Go To 111
  end if
  write(6,*)
  write(6,*) ' Error in QMSUrrounding section, DAMPing subsection.'
  call Quit(_RC_INPUT_ERROR_)
else if (Kword(1:4) == 'END ') then
  ! <<<END qmsurrounding>>>
  Go To 1000
end if
! And here we only go if unrecognized keyword is encountered.
write(6,*)
write(6,*) ' Unrecognized keyword in the QMSUrrounding section:',Kword(1:4)
call Quit(_RC_INPUT_ERROR_)

! <<<SOLVent>>>  Specify stuff about the solvent. Usually, these
!                parameters should not be altered.

112 continue
Key = Get_Ln(LuRd)
Kword = trim(Key)
call UpCase(Kword)
if (Kword(1:4) == 'EXCH') then
  ! <<<EXCHange>>>  Exchange repulsion parameters to solvent-solvent.
  do i=1,nAtom
    do j=1,i
      Key = Get_Ln(LuRd)
      call Get_F(1,Sexrep(i,j),1)
      call Get_F(2,Sexre1(i,j),1)
      call Get_F(3,Sexre2(i,j),1)
      Sexrep(j,i) = Sexrep(i,j)
      Sexre1(j,i) = Sexre1(i,j)
      Sexre2(j,i) = Sexre2(i,j)
    end do
  end do
  Go To 112
end if
if (Kword(1:4) == 'DISP') then
  ! <<<DISPersion>>>  Dispersion parameters to solvent-solvent.
  do i=1,nPol
    do j=1,i
      Key = Get_Ln(LuRd)
      call Get_F(1,Disp(i,j),1)
      Disp(j,i) = Disp(i,j)
    end do
  end do
  Go To 112
end if
if (Kword(1:4) == 'COOR') then
  ! <<<COORdinates>>>  Explicitly given coordinates of solvent
  !                    molecules. Also need number of particles.
  YesNo(6) = .true. !Signify that user gives coordinates.
  Key = Get_Ln(LuRd)
  call Get_I1(1,nPart)
  kaunt = 0
  do i=1,nPart
    do j=1,nAtom
      kaunt = kaunt+1
      Key = Get_Ln(LuRd)
      call Get_F(1,Cordst(kaunt,1),1)
      call Get_F(2,Cordst(kaunt,2),1)
      call Get_F(3,Cordst(kaunt,3),1)
    end do
    do kk=1,3
      CoTEMP1(kk) = Cordst(kaunt-2,kk)
      CoTEMP2(kk) = Cordst(kaunt-1,kk)
      CoTEMP3(kk) = Cordst(kaunt-0,kk)
    end do
    call OffAtom(CoTEMP1,CoTEMP2,CoTEMP3,CoTEMP4,CoTEMP5)
    kaunt = kaunt+1
    do kk=1,3
      Cordst(kaunt,kk) = CoTEMP4(kk)
    end do
    kaunt = kaunt+1
    do kk=1,3
      Cordst(kaunt,kk) = CoTEMP5(kk)
    end do
  end do
  Go To 112
end if
if (Kword(1:4) == 'CAVR') then
  ! <<<CAVRepulsion>>>  Repulsion parameters between solvent and cavity boundary.
  Key = Get_Ln(LuRd)
  call Get_F1(1,Exdtal)
  call Get_F1(2,Exdt1)
  Go To 112
end if
if (Kword(1:4) == 'OCOR') then
  ! <<<OCORbitals>>>  Occupied Orbitals for the solvent molecule
  Key = Get_Ln(LuRd)
  call Get_I(1,iOrb(2),1)
  Go To 112
end if
if (Kword(1:4) == 'ATCE') then
  ! <<<ATCEchpol>>>  Number of atoms, centers, charges and polarizabilities.
  ! Jose Slater Sites
  Key = Get_Ln(LuRd)
  call Get_I1(1,nAtom)
  call Get_I1(2,nCent)
  call Get_I1(3,nCha)
  call Get_I1(4,nPol)
  call Get_I1(5,nSlSiteC)
  Go To 112
end if
if (Kword(1:4) == 'CHAR') then
  ! <<<CHARge>>>  Magnitude of the charges.
  Key = Get_Ln(LuRd)
  call Get_F(1,Qsta,nCha)
  Go To 112
end if
if (Kword(1:4) == 'POLA') then
  ! <<<POLArizability>>>  Magnitude of polarizabilities.
  Key = Get_Ln(LuRd)
  call Get_F(1,Pol,nPol)
  Go To 112
end if
!Jose+++++++++++++
if (Kword(1:4) == 'SLAT') then
  ! <<<SLATer>>> Magnitude of Slater PreFactors and Exponents.
  Key = Get_Ln(LuRd)
  call Get_I1(1,lMltSlC)
  if (lMltSlC > 1) then
    write(6,*)
    write(6,*) 'Too high order of multipole in classical system'
    write(6,*) '              Higher order is 1'
    call Quit(_RC_INPUT_ERROR_)
  end if
  do i=1,nSlSiteC
    do j=0,lMltSlC
      nS = j*(j+1)*(j+2)/6
      nT = (j+1)*(j+2)*(j+3)/6
      Key = Get_Ln(LuRd)
      call Get_F1(1,SlExpTemp)
      SlExpC(j+1,i) = SlExpTemp
      njhr = nT-nS
      Key = Get_Ln(LuRd)
      call Get_F(1,SlFacTemp,njhr)
      njhr = 1
      do k=nS+1,nT
        SlFactC(k,i) = SlFacTemp(njhr)
        njhr = njhr+1
      end do
    end do
    Key = Get_Ln(LuRd)
    call Get_F(1,SlPC(i),1)
  end do
  Go To 112
end if

if (Kword(1:4) == 'END ') then
  ! <<<END solvent>>>
  Go To 1000
end if
! And the bla bla bla if something gets wrong.
write(6,*)
write(6,*) ' Unrecognized keyword in the SOLVent section:',Kword(1:4)
call Quit(_RC_INPUT_ERROR_)

! <<<RASSisection>>>  Give some numbers specific for the handling
!                     of the RASSI-construction of the wave-func.

113 continue
QmType = 'RASSI'
Key = Get_Ln(LuRd)
Kword = trim(Key)
call UpCase(Kword)
if (Kword(1:4) == 'JOBF') then
  ! <<<JOBFiles>>>  How many jobfiles and how many states in them.
  Key = Get_Ln(LuRd)
  call Get_I1(1,NrFiles)
  Key = Get_Ln(LuRd)
  call Get_I(1,NrStates,NrFiles)
  Go To 113
end if
if (Kword(1:4) == 'EQST') then
  ! <<<EQSTate>>> Which state is to be equilibrated.
  Key = Get_Ln(LuRd)
  call Get_I1(1,nEqState)
  Go To 113
end if
if (Kword(1:4) == 'MORE') then
  ! <<<MOREduce>>> Work in reduced MO-basis.
  Key = Get_Ln(LuRd)
  call Get_F1(1,ThrsRedOcc)
  MoAveRed = .true.
  Go To 113
end if
if (Kword(1:4) == 'CONT') then
  ! <<<CONTract>>> Contract the RASSI state basis.
  Key = Get_Ln(LuRd)
  call Get_F1(1,ThrsCont)
  ContrStateB = .true.
  Go To 113
end if
if (Kword(1:4) == 'LEVE') then
  ! <<<LEVElshift>>> Introduce levelshift of RASSI states.
  Key = Get_Ln(LuRd)
  call Get_I1(1,nLvlShift)
  Key = Get_Ln(LuRd)
  call Get_I(1,iLvlShift,nLvlShift)
  Key = Get_Ln(LuRd)
  call Get_F(1,dLvlShift,nLvlShift)
  ! Just a little sorting.
7485 continue
  Changed = .false.
  do i=1,nLvlShift-1
    if (iLvlShift(i) > iLvlShift(i+1)) then
      iTemp = iLvlShift(i)
      iLvlShift(i) = iLvlShift(i+1)
      iLvlShift(i+1) = iTemp
      dTemp = dLvlShift(i)
      dLvlShift(i) = dLvlShift(i+1)
      dLvlShift(i+1) = dTemp
      Changed = .true.
    end if
  end do
  if (Changed) goto 7485
  Go To 113
end if
if (Kword(1:4) == 'CISE') then
  ! <<<CISElect>>> Use overlap criterion in choosing state.
  lCiSelect = .true.
  Key = Get_Ln(LuRd)
  call Get_I1(1,nCIRef)
  Key = Get_Ln(LuRd)
  call Get_I(1,iCIInd,nCIRef)
  Key = Get_Ln(LuRd)
  call Get_F(1,dCIRef,nCIRef)
  goto 113
end if
if (Kword(1:4) == 'END ') then
  ! <<<END rassisection>>>
  YesNo(3) = .true.  !Rassi section has been visited.
  Go To 1000
end if
! HOW COULD IT GET WRONG HERE?
write(6,*)
write(6,*) ' Unrecognized keyword in the RASSisection section:',Kword(1:4)
call Quit(_RC_INPUT_ERROR_)

! <<<SCFSection>>>   Numbers for a SCF-QmStat run.

114 continue
QmType = 'SCF'
Key = Get_Ln(LuRd)
Kword = trim(Key)
call UpCase(Kword)
if (Kword(1:4) == 'ORBI') then
  ! <<<ORBItals>>>  Specifiy the reduced orbital space in which the problem is solved.
  Key = Get_Ln(LuRd)
  call Get_I(1,iOrb(1),1)
  call Get_I1(2,iOcc1)
  if (iOrb(1) > MxOrb) then
    write(6,*)
    write(6,*) 'The parameter MxOrb is set too low, or your total number of orbitals too high.'
    call Quit(_RC_INPUT_ERROR_)
  end if
  Go To 114
end if
if (Kword(1:4) == 'MP2D') then
  ! <<<MP2Denscorr>>>
  Mp2DensCorr = .true.
  Go To 114
end if
if (Kword(1:4) == 'END ') then
  ! <<<END scfsection>>>
  YesNo(4) = .true.  !Scf section has been visited.
  Go To 1000
end if
! ETWAS FALSCH!
write(6,*)
write(6,*) ' Unrecognized keyword in the SCFSection:',Kword(1:4)
call Quit(_RC_INPUT_ERROR_)

! <<<SINGle-point>>>   Signify that a set of single point calculations are to be done.

115 continue
SingPoint = .true.
YesNo(7) = .true.
Go To 1000

! <<<ANALyze section>>> Give details what analysis of the sampfile coordinates that is to be done.

116 continue
Go To 1000

! <<<EXTRact section>>> Give details what QM and QM/MM analysis
!                       that is to be done from the sampfile coordinates.

117 continue
Key = Get_Ln(LuRd)
Kword = trim(Key)
call UpCase(Kword)
if (Kword(1:4) == 'TOTA') then
  ! <<<TOTAl energy>>>
  lExtr(1) = .true.
  Go To 117
end if
if (Kword(1:4) == 'DIPO') then
  ! <<<DIPOle>>>
  lExtr(2) = .true.
  Go To 117
end if
if (Kword(1:4) == 'QUAD') then
  ! <<<QUADrupole>>>
  lExtr(3) = .true.
  Go To 117
end if
if (Kword(1:4) == 'EIGE') then
  ! <<<EIGEn things>>>
  lExtr(4) = .true.
  Key = Get_Ln(LuRd)
  call Get_I1(1,iExtr_Eig)
  call Get_S(2,VecsQue,1)
  call UpCase(VecsQue)
  if (VecsQue(1:3) == 'YES') lExtr(5) = .true.
  Go To 117
end if
if (Kword(1:4) == 'EXPE') then
  ! <<<EXPEctation values>>>
  lExtr(6) = .true.
  Go To 117
end if
if (Kword(1:4) == 'ELOC') then
  ! <<<ELOCal>>>
  lExtr(7) = .true.
  Key = Get_Ln(LuRd)
  call Get_I1(1,NExtr_Atm)
  Key = Get_Ln(LuRd)
  call Get_I(1,iExtr_Atm,NExtr_Atm)
  Go To 117
end if
!*****JoseMEP****************
if (Kword(1:4) == 'MESP') then
  ! <<<MESP>>>
  ! The Main Electrostatic potential, field and field gradients will
  ! be obtained in order to produce perturbation integrals that will
  ! be used to optimize the intramolecular geometry of the QM system.
  lExtr(8) = .true.
  Go To 117
end if
!*****************************
if (Kword(1:4) == 'END ') then
  ! <<<END extract section>>>
  YesNo(10) = .true.
  goto 1000
end if
! ETWAS FALSCH!
write(6,*)
write(6,*) ' Unrecognized keyword in the EXTRact section:',Kword(1:4)
call Quit(_RC_INPUT_ERROR_)

! Exit

99999 continue

! Check if mandatory input was included and that no blatant
! inconsistencies exist. Not fool-proof, fools!

call MandatoryInp(YesNo)

! Good bye.

return

end subroutine Get_Qmstat_Input
