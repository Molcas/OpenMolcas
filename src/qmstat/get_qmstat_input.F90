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
! qmstat_global which in turn are initialized in qmstat_init.
subroutine Get_Qmstat_Input(iQ_Atoms)

use qmstat_global, only: AddExt, Anal, ATitle, CAFieldG, CBFieldG, cDumpForm, CFExp, CharDi, CharDiQ, ContrStateB, Cordst, &
                         Cut_Elc, Cut_Ex1, Cut_Ex2, dCIRef, DelFi, DelOrAdd, DelR, DelX, Diel, DifSlExp, Disp, DispDamp, dLJRep, &
                         dLvlShift, EdSt, Enelim, ExtLabel, Exdt1, Exdtal, Exrep10, Exrep2, Exrep4, Exrep6, FieldDamp, Forcek, &
                         iCIInd, iCompExt, iExtr_Atm, iExtr_Eig, iExtra, iLuSaIn, iLuSaUt, iLuStIn, iLuStUt, iLvlShift, iNrIn, &
                         iNrUt, Inter, iOcc1, iOrb, iPrint, iRead, iSeed, itMax, Joblab, lCiSelect, lExtr, lMltSlC, lQuad, &
                         lSlater, MoAveRed, Mp2DensCorr, MxPut, nAdd, nAtom, nCent, nCha, nCIRef, nDel, nEqState, nExtAddOns, &
                         nLvlShift, nMacro, nMicro, nPart, nPol, NrFiles, NrStarti, NrStartu, NrStates, nSlSiteC, nStFilT, nTemp, &
                         ParallelT, ParaTemps, Pol, Pollim, Pres, Qmeq, QmProd, QmType, Qsta, QuaDi, QuaDiQ, rStart, SaFilIn, &
                         SaFilUt, ScalExt, Sexre1, Sexre2, Sexrep, SimEx, SingPoint, SlExpC, SlFactC, SlPC, StFilIn, StFilUt, &
                         Surf, Temp, ThrsCont, ThrsRedOcc, Udisp
use Index_Functions, only: nTri3_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iQ_Atoms
#include "warnings.h"
integer(kind=iwp) :: i, iChrct, iNrExtr, iTemp, j, kaunt, Last, LuRd, NExtr_Atm, njhr, nS, nT
real(kind=wp) :: CoTEMP1(3), CoTEMP2(3), CoTEMP3(3), CoTEMP4(3), CoTEMP5(3), dTemp, SlExpTemp, SlFacTemp(6)
logical(kind=iwp) :: Changed, YesNo(20)
character(len=180) :: Key
character(len=20) :: Kword
character(len=3) :: VecsQue
integer(kind=iwp), external :: iClast, IsFreeUnit
character(len=180), external :: Get_Ln
real(kind=wp), allocatable :: Tmp(:), Tmp2(:,:)

! Say what is done and set all YesNo to false; their purpose is to
! keep track on compulsory keywords and certain keyword combinations.

!write(u6,*)
!write(u6,*)'Input processed...'
YesNo(:) = .false.

! Use some nice routines to collect input.

LuRd = IsFreeUnit(79)
call SpoolInp(LuRd)
rewind(LuRd)
call RdNlst(LuRd,'QMSTAT')

! The turning-point in this do-while loop.

do

  ! Use Get_Ln to read the lines; it takes care of commented lines.

  Key = Get_Ln(LuRd)
  Kword = trim(Key)
  call UpCase(Kword)

  ! The keywords and their labels.

  select case (Kword(1:4))
    case default
      ! This code is only reached if an illegal keyword in the
      ! first tier is encountered.

      iChrct = len(Kword)
      Last = iCLast(Kword,iChrct)
      write(u6,*)
      write(u6,*) 'ERROR!'
      write(u6,'(1X,A,A)') Kword(1:Last),' is not a valid keyword!'
      call Quit(_RC_INPUT_ERROR_)

    case ('TITL')
      ! <<<TITLe>>>   Read title.

      Key = Get_Ln(LuRd)
      Joblab = trim(Key)
      ATitle = .true.

    case ('SIMU')
      ! <<<SIMUlation parameters>>>   Read a variety of simulation parameters.

      do
        Key = Get_Ln(LuRd)
        Kword = trim(Key)
        call UpCase(Kword)
        select case (KWord(1:4))
          case default
            ! Here we come if something gets wrong
            write(u6,*)
            write(u6,*) ' Unrecognized keyword in the SIMUlation parameter section:',Kword(1:4)
            call Quit(_RC_INPUT_ERROR_)
          case ('RADI')
            ! <<<RADIe>>>
            Key = Get_Ln(LuRd)
            call Get_F1(1,rStart)
          case ('PERM')
            ! <<<PERMitivity>>>
            Key = Get_Ln(LuRd)
            call Get_F1(1,Diel)
          case ('TEMP')
            ! <<<TEMPerature>>>
            Key = Get_Ln(LuRd)
            call Get_F1(1,Temp)
          case ('PRES')
            ! <<<PRESsure>>>
            Key = Get_Ln(LuRd)
            call Get_F1(1,Pres)
          case ('SURF')
            ! <<<SURFace>>>
            Key = Get_Ln(LuRd)
            call Get_F1(1,Surf)
          case ('TRAN')
            ! <<<TRANslation>>>
            Key = Get_Ln(LuRd)
            call Get_F1(1,DelX)
          case ('ROTA')
            ! <<<ROTAtion>>>
            Key = Get_Ln(LuRd)
            call Get_F1(1,DelFi)
          case ('CAVI')
            ! <<<CAVIty>>>
            Key = Get_Ln(LuRd)
            call Get_F1(1,DelR)
          case ('FORC')
            ! <<<FORCe>>>
            Key = Get_Ln(LuRd)
            call Get_F1(1,Forcek)
          case ('BREP')
            ! <<<BREPulsion>>>
            Key = Get_Ln(LuRd)
            call Get_F1(1,dLJRep)
          case ('SEED')
            ! <<<SEED>>>
            Key = Get_Ln(LuRd)
            call Get_I1(1,iSeed)
          case ('PARA')
            ! <<<PARAlleltemp>>>
            ParallelT = .true.
            Key = Get_Ln(LuRd)
            call Get_I1(1,nTemp)
            call mma_allocate(nStFilT,nTemp,label='nStFilT')
            call mma_allocate(Paratemps,nTemp,label='Paratemps')
            Key = Get_Ln(LuRd)
            call Get_I(1,nStFilT,nTemp)
            Key = Get_Ln(LuRd)
            call Get_F(1,ParaTemps,nTemp)
          case ('END ')
            ! <<<END simulation parameters>>>
            exit
        end select
      end do

    case ('THRE')
      ! <<<THREshold>>>    Get the polarization thresholds.

      Key = Get_Ln(LuRd)
      call Get_F1(1,Pollim)
      call Get_F1(2,Enelim)
      call Get_I1(3,itMax)

    case ('STEP')
      ! <<<STEPs>>>   Specify how many macro- and microsteps.

      Key = Get_Ln(LuRd)
      call Get_I1(1,nMacro)
      call Get_I1(2,nMicro)

    case ('RUN ')
      ! <<<RUN >>>   What type of simulation are we to run?

      Key = Get_Ln(LuRd)
      Kword = trim(Key)
      call UpCase(Kword)
      select case (Kword(1:4))
        case ('ANAL')
          Anal = .true.
        case ('QMEQ')
          Qmeq = .true.
        case ('QMPR')
          QmProd = .true.
          read(LuRd,*) Inter,iNrUt
          write(SaFilUt(6:6),'(i1.1)') iNrUt
          iLuSaUt = 32+iNrUt
        case default
          if (Kword(1:2) == 'SM') then
            write(u6,*)
            write(u6,*) 'No classical simulations are available.'
            call Quit(_RC_INPUT_ERROR_)
          end if
      end select
      YesNo(8) = .true.

    case ('PRIN')
      ! <<<PRINt>>>   Specify print-level.

      Key = Get_Ln(LuRd)
      call Get_I1(1,iPrint)

    case ('EXTE')
      ! <<<EXTErnal>>>   External one-electron perturbation
      !                  should be added on the hamiltonian.

      AddExt = .true.
      Key = Get_Ln(LuRd)
      call Get_I1(1,nExtAddOns)
      call mma_allocate(ScalExt,nExtAddOns,label='ScalExt')
      call mma_allocate(ExtLabel,nExtAddOns,label='ExtLabel')
      call mma_allocate(iCompExt,nExtAddOns,label='iCompExt')
      do i=1,nExtAddOns
        read(LuRd,*) ScalExt(i),ExtLabel(i),iCompExt(i)
      end do

    case ('EDIT')
      ! <<<EDITstartfile>>> Section for editing and displaying stuff on given startfile.

      do
        EdSt = .true.
        Key = Get_Ln(LuRd)
        Kword = trim(Key)
        call UpCase(Kword)
        select case (Kword(1:4))
          case default
            ! Here we come if something gets wrong
            write(u6,*)
            write(u6,*) ' Unrecognized keyword in the EDITstartfile section:',Kword(1:4)
            call Quit(_RC_INPUT_ERROR_)
          case ('DELE')
            ! <<<DELEte>>>   Delete solvent molecules.
            DelOrAdd(1) = .true.
            Key = Get_Ln(LuRd)
            call Get_I1(1,NrStarti)
            call Get_I1(2,NrStartu)
            Key = Get_Ln(LuRd)
            call Get_I1(1,nDel)
          case ('ADD ')
            ! <<<ADD >>>    Add solvent molecules.
            DelOrAdd(2) = .true.
            Key = Get_Ln(LuRd)
            call Get_I1(1,NrStarti)
            call Get_I1(2,NrStartu)
            Key = Get_Ln(LuRd)
            call Get_I1(1,nAdd)
          case ('QMDE')
            ! <<<QMDElete>>>   Substitute all slots with non-water coordinates with water coordinates.
            DelOrAdd(3) = .true.
            Key = Get_Ln(LuRd)
            call Get_I1(1,NrStarti)
            call Get_I1(2,NrStartu)
          case ('DUMP')
            ! <<<DUMP coordinates>>>  Dump coordinates in a way suitable for graphical display.
            DelOrAdd(4) = .true.
            Key = Get_Ln(LuRd)
            call UpCase(Key)
            cDumpForm = Key(1:4)
            Key = Get_Ln(LuRd)
            call Get_I1(1,NrStarti)
          case ('END ')
            ! <<<END editstartfile>>>
            write(StFilIn(6:6),'(i1.1)') NrStarti
            write(StFilUt(6:6),'(i1.1)') NrStartu
            exit
        end select
      end do

    case ('CONF')
      ! <<<CONFiguration>>>    Where is the initial configuration to be
      !                        obtained. Also, if we wish to edit the startfile.

      do
        Key = Get_Ln(LuRd)
        Kword = trim(Key)
        call UpCase(Kword)
        select case (Kword(1:4))
          case default
            ! Here we come if something gets wrong
            write(u6,*)
            write(u6,*) ' Unrecognized keyword in the CONFiguration section:',Kword(1:4)
            call Quit(_RC_INPUT_ERROR_)
          case ('ADD ')
            ! <<<ADD >>>  How many to add at random.
            Key = Get_Ln(LuRd)
            call Get_I1(1,iExtra)
            if (iExtra > MxPut) then
              write(u6,*)
              write(u6,*) 'The present limit of explicit solvent molecules is',MxPut,'.'
              call Quit(_RC_INPUT_ERROR_)
            end if
          case ('FILE')
            ! <<<FILE>>>  Read configuration/s from a file and put them there.
            Key = Get_Ln(LuRd)
            Kword = trim(Key)
            call UpCase(Kword)
            select case (Kword(1:4))
              case default
                ! CRASH-BOOM-BANG!
                write(u6,*)
                write(u6,*) ' Error in CONFiguration section, FILE subsection.'
                call Quit(_RC_INPUT_ERROR_)
              case ('STAR')
                ! <<<STARtfile>>>  Read from startfile.
                Key = Get_Ln(LuRd)
                Kword = trim(Key)
                call UpCase(Kword)
                select case (KWord(1:4))
                  case default
                    write(u6,*)
                    write(u6,*) 'Illegal StartFile option.'
                    call Quit(_RC_INPUT_ERROR_)
                  case ('SCRA')
                    ! <<<SCRAtch>>> Just put QM as given on RUNFILE.
                    Key = Get_Ln(LuRd)
                    call Get_I1(1,iNrIn)
                    call Get_I1(2,iNrUt)
                    iRead = 8
                  case ('COPY')
                    ! <<<COPY>>> Collect place of QM from startfile. WARNING!
                    !            You must use consistent startfile and RUNFILE!
                    Key = Get_Ln(LuRd)
                    call Get_I1(1,iNrIn)
                    call Get_I1(2,iNrUt)
                    iRead = 7
                  case ('CM  ')
                    ! <<<CM  >>> Put QM in CM of QM-place on startfile.
                    Key = Get_Ln(LuRd)
                    call Get_I1(1,iNrIn)
                    call Get_I1(2,iNrUt)
                    iRead = 6
                end select
              case ('SAMP')
                ! <<<SAMPfile>>>  Read configurations from sampfile and collect
                !                 the extracted information in iNrExtr.
                Key = Get_Ln(LuRd)
                call Get_I1(1,iNrIn)
                call Get_I1(2,iNrExtr)
                iRead = 9
                write(SimEx(6:6),'(i1.1)') iNrExtr
                YesNo(9) = .true.
            end select
            write(StFilIn(6:6),'(i1.1)') iNrIn
            write(StFilUt(6:6),'(i1.1)') iNrUt
            write(SaFilIn(6:6),'(i1.1)') iNrIn
            iLuStIn = 8+iNrIn
            iLuStUt = 16+iNrUt
            iLuSaIn = 24+iNrIn
          case ('INPU')
            ! <<<INPUt>>>  Signify that the first configuration will be given
            !              explicitly in input. The coordinates are then given
            !              in the solvent section.
            Key = Get_Ln(LuRd)
            call Get_I1(1,iNrUt)
            iLuStUt = 16+iNrUt
            YesNo(5) = .true. !User will give coords in input.
          case ('END ')
            YesNo(2) = .true. !Signify that source of starting conf. is specified.
            exit
        end select
      end do

    case ('QMSU')
      ! <<<QMSUrrounding>>>  Give parameters for the QM-Stat.Mech. interaction.

      do
        Key = Get_Ln(LuRd)
        Kword = trim(Key)
        call UpCase(Kword)
        select case (Kword(1:4))
          case default
            ! Here we only go if unrecognized keyword is encountered.
            write(u6,*)
            write(u6,*) ' Unrecognized keyword in the QMSUrrounding section:',Kword(1:4)
            call Quit(_RC_INPUT_ERROR_)
          case ('DPAR')
            ! <<<DPARameters>>>  Dispersion
            do i=1,iQ_Atoms
              Key = Get_Ln(LuRd)
              call Get_F(1,Udisp(1,i),1)
              call Get_F(2,Udisp(2,i),1)
            end do
          case ('ELEC')
            ! <<<ELECtrostatic>>> Electrostatic Slater Numbers
            do
              Key = Get_Ln(LuRd)
              Kword = trim(Key)
              call UpCase(Kword)
              select case (Kword(1:4))
                case default
                  write(u6,*)
                  write(u6,*) ' Error in QMSUrrounding section, ELECtrostatic subsection.'
                  call Quit(_RC_INPUT_ERROR_)
                case ('THRE')
                  ! <<<THREsholds>>>  First is the Cutoff (distance Quantum Site-
                  !                   Classical molecule) to evaluate Penetration
                  !                   effects. Second, difference between two Slater
                  !                   exponents to not be consider the same value.
                  Key = Get_Ln(LuRd)
                  call Get_F1(1,Cut_Elc)
                  call Get_F1(2,DifSlExp)
                case ('NOPE')
                  ! <<<NOPEnetration>>>  Electrostatic Penetration Not Computed
                  lSlater = .false.
                case ('QUAD')
                  ! <<QUADrupoles>>>  Electrostatic Penetration Computed in quadrupoles.
                  lQuad = .true.
                case ('END ')
                  ! <<<END Electrostatic>>>
                  exit
              end select
            end do
          case ('XPAR')
            ! <<<XPARameters>>>  Exchange repulsion
            do
              Key = Get_Ln(LuRd)
              Kword = trim(Key)
              call UpCase(Kword)
              select case (Kword(1:4))
                case default
                  write(u6,*)
                  write(u6,*) ' Error in QMSUrrounding section, XPARameters subsection.'
                  call Quit(_RC_INPUT_ERROR_)
                case ('S2  ')
                  ! <<<S2  >>>  The S2 parameter
                  Key = Get_Ln(LuRd)
                  call Get_F1(1,Exrep2)
                case ('S4  ')
                  ! <<<S4  >>>  The S4 parameter
                  Key = Get_Ln(LuRd)
                  call Get_F1(1,Exrep4)
                case ('S6  ')
                  ! <<<S6  >>>  The S6 parameter
                  Key = Get_Ln(LuRd)
                  call Get_F1(1,Exrep6)
                case ('S10 ')
                  ! <<<S10 >>>  The S10 parameter
                  Key = Get_Ln(LuRd)
                  call Get_F1(1,Exrep10)
                case ('CUTO')
                  ! <<<CUTOff>>> The cut-off radii for repulsion. The first is
                  !              outer radius that says EX=0 if R > Cut_Ex1, while
                  !              the second is a EX=infinity if R < Cut_Ex2.
                  Key = Get_Ln(LuRd)
                  call Get_F1(1,Cut_Ex1)
                  call Get_F1(2,Cut_Ex2)
                case ('END ')
                  ! <<<END xparameters>>>
                  exit
              end select
            end do
          case ('DAMP')
            ! <<<DAMPing>>>  Damping parameters for Qm-Surrounding interaction.
            do
              Key = Get_Ln(LuRd)
              Kword = trim(Key)
              call UpCase(Kword)
              select case (Kword(1:4))
                case default
                  write(u6,*)
                  write(u6,*) ' Error in QMSUrrounding section, DAMPing subsection.'
                  call Quit(_RC_INPUT_ERROR_)
                case ('DISP')
                  ! <<<DISPersion>>>  Dispersion damping parameters. This part
                  !                   should have a AUTO keyword which collects
                  !                   default parameters from the MpProp file.
                  DispDamp = .true.
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
                  call mma_allocate(CharDiQ,iQ_Atoms,label='CharDiQ')
                  call mma_allocate(QuaDiQ,3,iQ_Atoms,label='QuaDiQ')
                  do i=1,iQ_Atoms
                    Key = Get_Ln(LuRd)
                    call Get_F(1,CharDiQ(i),1)
                    call Get_F(2,QuaDiQ(1,i),1)
                    call Get_F(3,QuaDiQ(2,i),1)
                    call Get_F(4,QuaDiQ(3,i),1)
                  end do
                case ('FIEL')
                  ! <<<FIELd>>> Parameters for damping electric field.
                  FieldDamp = .true.
                  Key = Get_Ln(LuRd)
                  call Get_F1(1,CAFieldG)
                  call Get_F1(2,CBFieldG)
                  call Get_F1(3,CFExp)
                case ('END ')
                  ! <<<END damping>>>
                  exit
              end select
            end do
          case ('END ')
            ! <<<END qmsurrounding>>>
            exit
        end select
      end do

    case ('SOLV')
      ! <<<SOLVent>>>  Specify stuff about the solvent. Usually, these
      !                parameters should not be altered.

      do
        Key = Get_Ln(LuRd)
        Kword = trim(Key)
        call UpCase(Kword)
        select case (KWord(1:4))
          case default
            ! The bla bla bla if something gets wrong.
            write(u6,*)
            write(u6,*) ' Unrecognized keyword in the SOLVent section:',Kword(1:4)
            call Quit(_RC_INPUT_ERROR_)
          case ('EXCH')
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
          case ('DISP')
            ! <<<DISPersion>>>  Dispersion parameters to solvent-solvent.
            do i=1,nPol
              do j=1,i
                Key = Get_Ln(LuRd)
                call Get_F(1,Disp(i,j),1)
                Disp(j,i) = Disp(i,j)
              end do
            end do
          case ('COOR')
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
                call Get_F(1,Cordst(1,kaunt),1)
                call Get_F(2,Cordst(2,kaunt),1)
                call Get_F(3,Cordst(3,kaunt),1)
              end do
              CoTEMP1(:) = Cordst(:,kaunt-2)
              CoTEMP2(:) = Cordst(:,kaunt-1)
              CoTEMP3(:) = Cordst(:,kaunt-0)
              call OffAtom(CoTEMP1,CoTEMP2,CoTEMP3,CoTEMP4,CoTEMP5)
              kaunt = kaunt+1
              Cordst(:,kaunt) = CoTEMP4(:)
              kaunt = kaunt+1
              Cordst(:,kaunt) = CoTEMP5(:)
            end do
          case ('CAVR')
            ! <<<CAVRepulsion>>>  Repulsion parameters between solvent and cavity boundary.
            Key = Get_Ln(LuRd)
            call Get_F1(1,Exdtal)
            call Get_F1(2,Exdt1)
          case ('OCOR')
            ! <<<OCORbitals>>>  Occupied Orbitals for the solvent molecule
            Key = Get_Ln(LuRd)
            call Get_I(1,iOrb(2),1)
          case ('ATCE')
            ! <<<ATCEchpol>>>  Number of atoms, centers, charges and polarizabilities.
            ! Jose Slater Sites
            Key = Get_Ln(LuRd)
            call Get_I1(1,nAtom)
            call Get_I1(2,nCent)
            call Get_I1(3,nCha)
            call Get_I1(4,nPol)
            call Get_I1(5,nSlSiteC)
            ! Reallocate Qsta
            call mma_allocate(Tmp,max(size(Qsta),nCha),label='Tmp')
            Tmp(1:size(Qsta)) = Qsta
            call mma_deallocate(Qsta)
            call move_alloc(Tmp,Qsta)
            ! Reallocate Pol
            call mma_allocate(Tmp,max(size(Pol),nPol),label='Tmp')
            Tmp(1:size(Pol)) = Pol
            call mma_deallocate(Pol)
            call move_alloc(Tmp,Pol)
            ! Reallocate Disp
            call mma_allocate(Tmp2,max(size(Disp,1),nPol),max(size(Disp,2),nPol),label='Tmp')
            Tmp2(1:size(Disp,1),1:size(Disp,2)) = Disp
            call mma_deallocate(Disp)
            call move_alloc(Tmp2,Disp)
            ! Reallocate SlExpC
            call mma_allocate(Tmp2,4,max(size(SlExpC,2),nSlSiteC),label='Tmp')
            Tmp2(:,1:size(SlExpC,2)) = SlExpC
            call mma_deallocate(SlExpC)
            call move_alloc(Tmp2,SlExpC)
            ! Reallocate SlFactC
            call mma_allocate(Tmp2,4,max(size(SlFactC,2),nSlSiteC),label='Tmp')
            Tmp2(:,1:size(SlFactC,2)) = SlFactC
            call mma_deallocate(SlFactC)
            call move_alloc(Tmp2,SlFactC)
            ! Reallocate SlPC
            call mma_allocate(Tmp,max(size(SlPC),nSlSiteC),label='Tmp')
            Tmp(1:size(SlPC)) = SlPC
            call mma_deallocate(SlPC)
            call move_alloc(Tmp,SlPC)
            ! Reallocate Sexrep
            call mma_allocate(Tmp2,max(size(Sexrep,1),nAtom),max(size(Sexrep,2),nAtom),label='Tmp')
            Tmp2(1:size(Sexrep,1),1:size(Sexrep,2)) = Sexrep
            call mma_deallocate(Sexrep)
            call move_alloc(Tmp2,Sexrep)
            ! Reallocate Sexre1
            call mma_allocate(Tmp2,max(size(Sexre1,1),nAtom),max(size(Sexre1,2),nAtom),label='Tmp')
            Tmp2(1:size(Sexre1,1),1:size(Sexre1,2)) = Sexre1
            call mma_deallocate(Sexre1)
            call move_alloc(Tmp2,Sexre1)
            ! Reallocate Sexre2
            call mma_allocate(Tmp2,max(size(Sexre2,1),nAtom),max(size(Sexre2,2),nAtom),label='Tmp')
            Tmp2(1:size(Sexre2,1),1:size(Sexre2,2)) = Sexre2
            call mma_deallocate(Sexre2)
            call move_alloc(Tmp2,Sexre2)
          case ('CHAR')
            ! <<<CHARge>>>  Magnitude of the charges.
            Key = Get_Ln(LuRd)
            call Get_F(1,Qsta,nCha)
          case ('POLA')
            ! <<<POLArizability>>>  Magnitude of polarizabilities.
            Key = Get_Ln(LuRd)
            call Get_F(1,Pol,nPol)
          case ('SLAT')
            !Jose+++++++++++++
            ! <<<SLATer>>> Magnitude of Slater PreFactors and Exponents.
            Key = Get_Ln(LuRd)
            call Get_I1(1,lMltSlC)
            if (lMltSlC > 1) then
              write(u6,*)
              write(u6,*) 'Too high order of multipole in classical system'
              write(u6,*) '              Highest order is 1'
              call Quit(_RC_INPUT_ERROR_)
            end if
            do i=1,nSlSiteC
              do j=0,lMltSlC
                nS = nTri3_Elem(j)
                nT = nTri3_Elem(j+1)
                Key = Get_Ln(LuRd)
                call Get_F1(1,SlExpTemp)
                SlExpC(j+1,i) = SlExpTemp
                njhr = nT-nS
                Key = Get_Ln(LuRd)
                call Get_F(1,SlFacTemp,njhr)
                SlFactC(nS+1:nT,i) = SlFacTemp(1:nT-nS)
              end do
              Key = Get_Ln(LuRd)
              call Get_F(1,SlPC(i),1)
            end do
            !+++++++++++++++++
          case ('END ')
            ! <<<END solvent>>>
            exit
        end select
      end do

    case ('RASS')
      ! <<<RASSisection>>>  Give some numbers specific for the handling
      !                     of the RASSI-construction of the wave-func.

      QmType = 'RASSI'
      do
        Key = Get_Ln(LuRd)
        Kword = trim(Key)
        call UpCase(Kword)
        select case (Kword(1:4))
          case default
            ! HOW COULD IT GET WRONG HERE?
            write(u6,*)
            write(u6,*) ' Unrecognized keyword in the RASSisection section:',Kword(1:4)
            call Quit(_RC_INPUT_ERROR_)
          case ('JOBF')
            ! <<<JOBFiles>>>  How many jobfiles and how many states in them.
            Key = Get_Ln(LuRd)
            call Get_I1(1,NrFiles)
            call mma_allocate(NrStates,NrFiles,label='NrStates')
            Key = Get_Ln(LuRd)
            call Get_I(1,NrStates,NrFiles)
          case ('EQST')
            ! <<<EQSTate>>> Which state is to be equilibrated.
            Key = Get_Ln(LuRd)
            call Get_I1(1,nEqState)
          case ('MORE')
            ! <<<MOREduce>>> Work in reduced MO-basis.
            Key = Get_Ln(LuRd)
            call Get_F1(1,ThrsRedOcc)
            MoAveRed = .true.
          case ('CONT')
            ! <<<CONTract>>> Contract the RASSI state basis.
            Key = Get_Ln(LuRd)
            call Get_F1(1,ThrsCont)
            ContrStateB = .true.
          case ('LEVE')
            ! <<<LEVElshift>>> Introduce levelshift of RASSI states.
            Key = Get_Ln(LuRd)
            call Get_I1(1,nLvlShift)
            call mma_allocate(iLvlShift,nLvlShift,label='iLvlShift')
            call mma_allocate(dLvlShift,nLvlShift,label='dLvlShift')
            Key = Get_Ln(LuRd)
            call Get_I(1,iLvlShift,nLvlShift)
            Key = Get_Ln(LuRd)
            call Get_F(1,dLvlShift,nLvlShift)
            ! Just a little sorting.
            do
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
              if (.not. Changed) exit
            end do
          case ('CISE')
            ! <<<CISElect>>> Use overlap criterion in choosing state.
            lCiSelect = .true.
            call mma_allocate(iCIInd,nCIRef,label='iCIInd')
            call mma_allocate(dCIRef,nCIRef,label='dCIRef')
            Key = Get_Ln(LuRd)
            call Get_I1(1,nCIRef)
            Key = Get_Ln(LuRd)
            call Get_I(1,iCIInd,nCIRef)
            Key = Get_Ln(LuRd)
            call Get_F(1,dCIRef,nCIRef)
          case ('END ')
            ! <<<END rassisection>>>
            YesNo(3) = .true.  !Rassi section has been visited.
            exit
        end select
      end do

    case ('SCFS')
      ! <<<SCFSection>>>   Numbers for a SCF-QmStat run.

      QmType = 'SCF'
      do
        Key = Get_Ln(LuRd)
        Kword = trim(Key)
        call UpCase(Kword)
        select case (Kword(1:4))
          case default
            ! ETWAS FALSCH!
            write(u6,*)
            write(u6,*) ' Unrecognized keyword in the SCFSection:',Kword(1:4)
            call Quit(_RC_INPUT_ERROR_)
          case ('ORBI')
            ! <<<ORBItals>>>  Specifiy the reduced orbital space in which the problem is solved.
            Key = Get_Ln(LuRd)
            call Get_I(1,iOrb(1),1)
            call Get_I1(2,iOcc1)
          case ('MP2D')
            ! <<<MP2Denscorr>>>
            Mp2DensCorr = .true.
          case ('END ')
            ! <<<END scfsection>>>
            YesNo(4) = .true.  !Scf section has been visited.
            exit
        end select
      end do

    case ('SING')
      ! <<<SINGle-point>>>   Signify that a set of single point calculations are to be done.

      SingPoint = .true.
      YesNo(7) = .true.

    case ('ANAL')
      ! <<<ANALyze section>>> Give details what analysis of the sampfile coordinates that is to be done.

    case ('EXTR')
      ! <<<EXTRact section>>> Give details what QM and QM/MM analysis
      !                       that is to be done from the sampfile coordinates.

      do
        Key = Get_Ln(LuRd)
        Kword = trim(Key)
        call UpCase(Kword)
        select case (Kword(1:4))
          case default
            ! ETWAS FALSCH!
            write(u6,*)
            write(u6,*) ' Unrecognized keyword in the EXTRact section:',Kword(1:4)
            call Quit(_RC_INPUT_ERROR_)
          case ('TOTA')
            ! <<<TOTAl energy>>>
            lExtr(1) = .true.
          case ('DIPO')
            ! <<<DIPOle>>>
            lExtr(2) = .true.
          case ('QUAD')
            ! <<<QUADrupole>>>
            lExtr(3) = .true.
          case ('EIGE')
            ! <<<EIGEn things>>>
            lExtr(4) = .true.
            Key = Get_Ln(LuRd)
            call Get_I1(1,iExtr_Eig)
            call Get_S(2,VecsQue,1)
            call UpCase(VecsQue)
            if (VecsQue(1:3) == 'YES') lExtr(5) = .true.
          case ('EXPE')
            ! <<<EXPEctation values>>>
            lExtr(6) = .true.
          case ('ELOC')
            ! <<<ELOCal>>>
            lExtr(7) = .true.
            Key = Get_Ln(LuRd)
            call Get_I1(1,NExtr_Atm)
            call mma_deallocate(iExtr_Atm)
            call mma_allocate(iExtr_Atm,NExtr_Atm,label='iExtr_Atm')
            Key = Get_Ln(LuRd)
            call Get_I(1,iExtr_Atm,NExtr_Atm)
          case ('MESP')
            !*****JoseMEP****************
            ! <<<MESP>>>
            ! The Main Electrostatic potential, field and field gradients will
            ! be obtained in order to produce perturbation integrals that will
            ! be used to optimize the intramolecular geometry of the QM system.
            lExtr(8) = .true.
            !*****************************
          case ('END ')
            ! <<<END extract section>>>
            YesNo(10) = .true.
        end select
      end do

    case ('END ')
      ! <<<END >>>

      exit

  end select
end do

! Check if mandatory input was included and that no blatant
! inconsistencies exist. Not fool-proof, fools!

call MandatoryInp(YesNo)

! Good bye.

return

end subroutine Get_Qmstat_Input
