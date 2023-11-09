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
! Copyright (C) Yannick Carissan                                       *
!               Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine Readinp_localisation()
! Author: Y. Carissan [heavily modified by T.B. Pedersen].

use Localisation_globals, only: AnaAtom, AnaDomain, Analysis, AnaNrm, AnaPAO, AnaPAO_Save, ChoStart, DoCNOs, DoDomain, EvalER, &
                                iWave, LocCanOrb, LocModel, LocNatOrb, LocPAO, LuSpool, Maximisation, MxConstr, nActa, NamAct, &
                                nConstr, nFro, NMxIter, nOccInp, nOrb, nOrb2Loc, nSym, nVirInp, Order, PrintMOs, Silent, Skip, &
                                Test_Localisation, ThrDomain, ThrGrad, ThrPairDomain, ThrRot, Thrs, ThrSel, Timing, Wave
#ifdef _DEBUGPRINT
use Localisation_globals, only: nBas
#endif
use stdalloc, only: mma_allocate
use Constants, only: Ten
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, iPL, istatus, iSym, j, LocOrb
character(len=180) :: Key, Line
logical(kind=iwp) :: Thrs_UsrDef, LocModel_UsrDef, nFro_UsrDef, nOrb2Loc_UsrDef, Freeze
integer(kind=iwp), parameter :: Occupied = 0, Virtual = 1, AllOrb = 2
real(kind=wp), parameter :: ThrsDef = 1.0e-6_wp, ThrRotDef = 1.0e-10_wp, ThrGradDef = 1.0e-2_wp
character(len=*), parameter :: SecNam = 'Readinp_localisation'
integer(kind=iwp), external :: iPrintLevel, isFreeUnit
character(len=180), external :: Get_Ln

LuSpool = 17
LuSpool = isFreeUnit(LuSpool)
call SpoolInp(LuSpool)

! Locate "start of input"
rewind(LuSpool)
call RdNLst(LuSpool,'LOCALISATION')

! Get print level

iPL = iPrintLevel(-1)

! Default Parameters

do iSym=1,nSym
  nOrb2Loc(iSym) = 0
  nFro(iSym) = 0
  nConstr(iSym) = 0
end do
Skip = .false.
LocOrb = Occupied
!LocVir = .false.
Thrs_UsrDef = .false.
nOrb2Loc_UsrDef = .false.
nFro_UsrDef = .false.
Freeze = .false.
Maximisation = .true.
ChoStart = .false.
if (iPL < 3) then
  Silent = .true.
else
  Silent = .false.
end if
LocModel = 1  ! Pipek-Mezey localisation
if (nSym > 1) LocModel = 3  ! Cholesky localisation
LocModel_UsrDef = .false.
Test_Localisation = .false.
NMxIter = 300
Thrs = ThrsDef
ThrRot = ThrRotDef
ThrGrad = ThrGradDef
Analysis = .false.
AnaAtom = nSym == 1
AnaNrm = 'Fro'
PrintMOs = .true.
Timing = .true.
EvalER = .false.
Order = .false.
LocPAO = .false.
AnaPAO = .false.
AnaPAO_Save = AnaPAO
DoDomain = .false.
AnaDomain = .false.
ThrDomain(1) = 0.9_wp
ThrDomain(2) = 2.0e-2_wp
ThrPairDomain(1) = 1.0e-10_wp
ThrPairDomain(2) = Ten
ThrPairDomain(3) = 15.0_wp
LocNatOrb = .false.
LocCanOrb = .false.
Wave = .false.
iWave = 0
DoCNOs = .false.

! End Default Parameters

do
  Key = Get_Ln(LuSpool)
  Line = Key
  call UpCase(Line)
  select case (Line(1:4))
    case ('NORB')
      ! NORBitals

      Line = Get_Ln(LuSpool)
      call Get_I(1,nOrb2Loc,nSym)
      nOrb2Loc_UsrDef = .true.

    case ('NFRO')
      ! NFROzen

      Line = Get_Ln(LuSpool)
      call Get_I(1,nFro,nSym)
      nFro_UsrDef = .true.
      Freeze = .false.
      !if (LocVir) then
      if (LocOrb == Virtual) then
        write(u6,*)
        write(u6,*) 'WARNING!!!'
        write(u6,*) 'You have chosen to freeze some orbitals AND asked to localize the virtual space. This implies that'
        write(u6,*) ' some virtual orbitals will be kept frozen, thus they will be left unchanged. Is that OK ?'
        write(u6,*)
      end if

    case ('FREE')
      ! FREEze core orbitals (as defined by seward)

      if (.not. nFro_UsrDef) then
        call Get_iArray('Non valence orbitals',nFro,nSym)
        Freeze = .true.
        !if (LocVir) then
        if (LocOrb == Virtual) then
          write(u6,*)
          write(u6,*) 'WARNING!!!'
          write(u6,*) 'You have chosen to freeze some orbitals AND asked to localize the virtual space. This implies that'
          write(u6,*) ' some virtual orbitals will be kept frozen, thus they will be left unchanged. Is that OK ?'
          write(u6,*)
        end if
      end if

    case ('MAXI')
      ! MAXImisation

      Maximisation = .true.

    case ('MINI')
      ! MINImisation

      Maximisation = .false.

    case ('NITE','ITER')
      ! NITErations or ITERations

      Line = Get_Ln(LuSpool)
      call Get_I1(1,NMxIter)

    case ('THRE')
      ! THREshold

      Line = Get_Ln(LuSpool)
      call Get_F1(1,Thrs)
      Thrs_UsrDef = .true.

    case ('THRG')
      ! THRGrad

      Line = Get_Ln(LuSpool)
      call Get_F1(1,ThrGrad)

    case ('CHOS')
      ! CHOStart (use Cholesky orbitals as start guess for PM/Boys/ER)

      ChoStart = .true.

    case ('THRR')
      ! THRRotation

      Line = Get_Ln(LuSpool)
      call Get_F1(1,ThrRot)

    case ('PIPE','PM  ')
      ! PIPEk-Mezey or PM

      LocModel = 1
      LocModel_UsrDef = .true.

    case ('BOYS')
      ! BOYS

      LocModel = 2
      LocModel_UsrDef = .true.

    case ('CHOL')
      ! CHOLesky

      LocModel = 3
      LocModel_UsrDef = .true.

    case ('EDMI','ER  ')
      ! EDMIston-Ruedenberg or ER

      LocModel = 4
      LocModel_UsrDef = .true.

    case ('SILE')
      ! SILEnt mode

      Silent = .true.

    case ('TEST')
      ! TEST localisation (orthonormality, density, etc.)

      Test_Localisation = .true.

    case ('ANAL')
      ! ANALysis: generate bitmaps + histograms for density, original, and
      !           local MOs. Analysis "per atom" is default, but does not work
      !           with symmetry.

      Analysis = .true.
      AnaAtom = nSym == 1

    case ('ANAA')
      ! ANAAtom: generate bitmaps + histograms for density, original, and
      !          local MOs. Analysis "per atom".

      Analysis = .true.
      AnaAtom = .true.

    case ('ANAS')
      ! ANAShell: generate bitmaps + histograms for density, original, and
      !           local MOs. Analysis "per shell".

      Analysis = .true.
      AnaAtom = .false.

    case ('MAX ')
      ! MAX norm: use max element as norm in analysis.

      AnaNrm = 'Max'

    case ('FROB')
      ! FROBenius norm: use Frobenius norm in analysis.

      AnaNrm = 'Fro'

    case ('NOMO')
      ! NOMO: do not print MOs.

      PrintMOs = .false.

    case ('TIME')
      ! TIME: time the localisation procedure explicitly

      Timing = .true.

    case ('NOTI')
      ! NOTI: do not time the localisation procedure explicitly

      Timing = .false.

    case ('ERFU')
      ! ERFUnctional: evaluate Edmiston-Ruedenberg functional with original
      !               and localised orbitals.

      EvalER = .true.

    case ('ORDE')
      ! ORDEr: order localised orbitals according to Cholesky orbital
      !        ordering, defined according to the overlap between the two sets
      !        orbitals.

      Order = .true.

    case ('VIRT')
      ! VIRTual: localise virtual orbitals

      LocOrb = Virtual
      !LocVir = .true.
      if (nFro_UsrDef .or. Freeze) then
        write(u6,*)
        write(u6,*) 'WARNING!!!'
        write(u6,*) 'You have chosen to freeze some orbitals AND asked to localize the virtual space. This implies that'
        write(u6,*) ' some virtual orbitals will be kept frozen, thus they will be left unchanged. Is that OK ?'
        write(u6,*)
      end if

    case ('OCCU')
      ! OCCUpied: localise occupied orbitals

      LocOrb = Occupied
      !LocVir = .false.

    case ('ALL ')
      ! ALL: localise all orbitals

      LocOrb = AllOrb

    case ('PAO ')
      ! PAO : compute projected AOs that span the virtual space
      !       (= all - occupied - frozen) using Cholesky decomposition to
      !       remove linear dependence.

      LocPAO = .true.
      LocModel = 3
      LocModel_UsrDef = .true.

    case ('ANAP')
      ! ANAP: Special analysis for Cholesky PAOs before orthonormalization.

      AnaPAO = .true.

    case ('DOMA')
      ! DOMAin: set up orbital domains (Pulay-style) and pair domains

      DoDomain = .true.

    case ('THRD')
      ! THRDomain: thresholds for setting up orbital domains.
      !            First value is the minimum sum of gross atomic Mulliken
      !            charge for each orbital.
      !            Second value is the threshold used for the completeness
      !            check of Boughton and Pulay.

      DoDomain = .true.
      Line = Get_Ln(LuSpool)
      call Get_F(1,ThrDomain,2)

    case ('THRP')
      ! THRPairdomain: thresholds for setting up pair domains.
      !                3 values needed: Rs, Rw, Rd (in bohr).
      !                Let R be the minimum distance between any two atoms in
      !                a pair domain. Then,
      !                     R <= Rs: strong pair
      !                Rs < R <= Rw: weak pair
      !                Rw < R <= Rd: distant pair
      !                Rd < R <= Rd: very distant pair

      DoDomain = .true.
      Line = Get_Ln(LuSpool)
      call Get_F(1,ThrPairDomain,3)

    case ('ANAD')
      ! ANADomain: analysis of orbital domains

      AnaDomain = .true.
      DoDomain = .true.

    case ('SKIP')
      ! SKIP: skip localisation completely; i.e. only perform analysis etc.

      Skip = .true.

    case ('LOCN','LOCC')
      ! LOCN: localized natural orbitals
      ! LOCC: localized canonical orbitals

      if (Line(1:4) == 'LOCN') then
        LocNatOrb = .true.
      else if (Line(1:4) == 'LOCC') then
        LocCanOrb = .true.
      end if
      read(LuSpool,*) nActa,ThrSel
      ! Read atom names
      do
        read(LuSpool,'(A)',iostat=istatus) Line
        if (istatus /= 0) call Error()
        if ((Line(1:1) /= '*') .and. (Line /= ' ')) exit
      end do
      call UpCase(Line)
      call mma_allocate(NamAct,nActa,label='NamAct')
      do i=1,nActa
        if (Line == ' ') call Error()
        Line = adjustl(Line)
        j = index(Line,' ')
        NamAct(i) = Line(1:j-1)
        Line(1:j-1) = ' '
      end do

    case ('WAVE')
      ! WAVE: wavelet transform of the MO basis

      read(LuSpool,*) iWave
      if ((iWave /= 0) .and. (iWave /= 1)) then
        write(u6,*) ' WARNING!!!'
        write(u6,*) ' Incorrect bit-switch input parameter for WAVE.'
        write(u6,*) ' I will continue with the default value.'
        iWave = 0
      end if
      Skip = .false.
      Wave = .true.
      LocModel = 0

    case ('CONS')
      ! CONS: constrained natural orbital analysis

      Line = Get_Ln(LuSpool)
      call Get_I(1,nConstr,nSym)
      MxConstr = 0
      do i=1,nSym
        MxConstr = max(MxConstr,nConstr(i))
      end do
      if (MxConstr > 16) then
        write(u6,*) ' ERROR  !!!'
        write(u6,*) ' Cannot handle more than 16 constraints/symm !!'
        write(u6,*) ' Increase size of indxC in Get_CNOs and recompile.'
      end if
      DoCNOs = .true.
      Skip = .false.
      PrintMOs = .false.
      LocModel = 0

    case ('FILE')
      ! FILE: filename with input orbitals

      Line = Get_Ln(LuSpool)
      ! Filename is read in localisation.F90

    case ('END ')
      exit

    case default
      write(u6,*) 'Unidentified key word  : ',Key
      write(u6,*) 'Internal representation: ',Line(1:4)
      call FindErrorLine()
      call Quit_OnUserError()
  end select
end do

! END of Input

! ==============
! Postprocessing
! ==============

! Only Cholesky localisation can be run with symmetry.
! Reset if the user explicitly requested a localisation
! model different from Cholesky.
! ------------------------------------------------------------------

if ((nSym > 1) .and. LocModel_UsrDef .and. (LocModel /= 3)) then
  write(u6,*)
  write(u6,*) 'WARNING!!!'
  write(u6,*) 'The localisation model you have suggested in combination with symmetry will not work.'
  write(u6,*) 'In this case the program defaults to the Cholesky localisation scheme!'
  write(u6,*)
  LocModel = 3
end if

! Special settings for PAO: LocModel must be 3 (i.e. Cholesky), else
! we cancel PAO. For PAO, special analysis may be activated by the
! DEBUGPRINT case.
! ------------------------------------------------------------------

LocPAO = LocPAO .and. (LocModel == 3)
if (LocPAO) then
# ifdef _DEBUGPRINT_
  AnaPAO = .true.
# endif
else
  AnaPAO = .false.
end if

! Set orbitals to localize.
! If the user specified which orbitals to localise, frozen orbitals
! is not our problem (the user must set those as well or use the
! default).
! If the user did not specify which orbitals to localise, we have
! several cases (remembering that virtual localisation is default
! for the PAO method):
! For virtual localisation (specified with keyword VIRTual):
!    user-defined frozen orbitals are regarded as frozen virtual,
!    else we freeze all occupied and localise the rest. Note that
!    keyword FREEze is ignored for virtual localisation!
! For occupied localisation (default for all methods but PAO):
!    user-defined frozen orbitals (possibly by FREEze keyword)
!    are subtracted from the set of occupied orbitals, else we
!    localise all occupied.
! -----------------------------------------------------------------

if (nOrb2Loc_UsrDef) then
  LocOrb = Occupied
  !LocVir = .false.
else
  !if ((LocOrb == Virtual) .or. LocPAO) then
  !  LocOrb = Virtual
  !end if
  !LocVir = LocVir .or. LocPAO ! virtual is default for PAO
  if (LocOrb == Virtual) then ! virtual localisation
    !if (LocVir) then ! virtual localisation
    if (nFro_UsrDef) then
      do iSym=1,nSym
        nFro(iSym) = nOccInp(iSym)+nFro(iSym)
        nOrb2Loc(iSym) = nOrb(iSym)-nFro(iSym)
      end do
    else
      do iSym=1,nSym
        nFro(iSym) = nOccInp(iSym)
        nOrb2Loc(iSym) = nVirInp(iSym)
      end do
    end if
  !else ! occupied localisation
  else if (LocOrb == Occupied) then ! occupied localisation
    if (nFro_UsrDef .or. Freeze) then
      do iSym=1,nSym
        nOrb2Loc(iSym) = nOccInp(iSym)-nFro(iSym)
      end do
    else
      do iSym=1,nSym
        nOrb2Loc(iSym) = nOccInp(iSym)
      end do
    end if
  else ! AllOrb
    if (nFro_UsrDef .or. Freeze) then
      do iSym=1,nSym
        nOrb2Loc(iSym) = nOccInp(iSym)+nVirInp(iSym)-nFro(iSym)
      end do
    else
      do iSym=1,nSym
        nOrb2Loc(iSym) = nOccInp(iSym)+nVirInp(iSym)
      end do
    end if
  end if
end if

#ifdef _DEBUGPRINT
write(u6,'(/,A,A)') SecNam,': orbital definitions:'
write(u6,'(A,8I9)') 'nBas    : ',(nBas(iSym),iSym=1,nSym)
write(u6,'(A,8I9)') 'nOrb    : ',(nOrb(iSym),iSym=1,nSym)
write(u6,'(A,8I9)') 'nOccInp : ',(nOccInp(iSym),iSym=1,nSym)
write(u6,'(A,8I9)') 'nVirInp : ',(nVirInp(iSym),iSym=1,nSym)
write(u6,'(A,8I9)') 'nFro    : ',(nFro(iSym),iSym=1,nSym)
write(u6,'(A,8I9,/)') 'nOrb2Loc: ',(nOrb2Loc(iSym),iSym=1,nSym)
#endif

! If Cholesky, reset default threshold (unless user defined).
! -----------------------------------------------------------

if ((LocModel == 3) .and. (.not. Thrs_UsrDef)) then
  Thrs = 1.0e-8_wp
end if

! No need to order Cholesky MOs.
! ------------------------------

Order = Order .and. (LocModel /= 3)

! No need to evaluate ER function for Edmiston-Ruedenberg.
! --------------------------------------------------------

EvalER = EvalER .and. (LocModel /= 4)

! Turn on localisation test if debug or analysis was specified.
! -------------------------------------------------------------

Test_Localisation = Test_Localisation .or. Analysis

contains

subroutine Error()
  write(u6,*) ' READIN: Premature end of file when reading selected'
  write(u6,*) ' atoms in keyword LOCN'
  call Abend()
end subroutine Error

end subroutine Readinp_localisation
