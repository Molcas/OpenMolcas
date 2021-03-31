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

implicit real*8(a-h,o-z)
#include "Molcas.fh"
#include "inflocal.fh"
#include "debug.fh"

!TBP  Namelist /LOCALISATION/ dummy

character*20 SecNam
parameter(SecNam='Readinp_localisation')

character*180 Key, Line, Blank
character*180 Get_Ln
external Get_Ln

integer iPrintLevel
external iPrintLevel

logical Thrs_UsrDef, LocModel_UsrDef, Freeze
parameter(ThrsDef=1.0d-6,ThrRotDef=1.0d-10)
parameter(ThrGradDef=1.0d-2)

LuSpool = 17
LuSpool = isFreeUnit(LuSpool)
call SpoolInp(LuSpool)

Debug = .false.
! Locate "start of input"
rewind(LuSpool)
call RdNLst(LuSpool,'LOCALISATION')
Blank = ' '

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
if (iPL >= 4) then
  Debug = .true.
else
  Debug = .false.
end if
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
ThrDomain(1) = 9.0d-1
ThrDomain(2) = 2.0d-2
ThrPairDomain(1) = 1.0d-10
ThrPairDomain(2) = 1.0d1
ThrPairDomain(3) = 1.5d1
LocNatOrb = .false.
LocCanOrb = .false.
Wave = .false.
iWave = 0
DoCNOs = .false.

! End Default Parameters

999 continue
Key = Get_Ln(LuSpool)
Line = Key
call UpCase(Line)
if (Line(1:4) == 'DEBU') Go To 1000
if (Line(1:4) == 'NORB') Go To 2000
if (Line(1:4) == 'NFRO') Go To 2100
if (Line(1:4) == 'FREE') Go To 2200
if (Line(1:4) == 'MAXI') Go To 3000
if (Line(1:4) == 'MINI') Go To 3100
if (Line(1:4) == 'NITE') Go To 4000
if (Line(1:4) == 'ITER') Go To 4000
if (Line(1:4) == 'THRE') Go To 5000
if (Line(1:4) == 'THRG') Go To 5001
if (Line(1:4) == 'CHOS') Go To 5100
if (Line(1:4) == 'THRR') Go To 6000
if (Line(1:4) == 'PIPE') Go To 7000
if (Line(1:4) == 'PM  ') Go To 7000
if (Line(1:4) == 'BOYS') Go To 7100
if (Line(1:4) == 'CHOL') Go To 7200
if (Line(1:4) == 'EDMI') Go To 7300
if (Line(1:4) == 'ER  ') Go To 7300
if (Line(1:4) == 'SILE') Go To 8000
if (Line(1:4) == 'TEST') Go To 9000
if (Line(1:4) == 'ANAL') Go To 10000
if (Line(1:4) == 'ANAA') Go To 10100
if (Line(1:4) == 'ANAS') Go To 10200
if (Line(1:4) == 'MAX ') Go To 10300
if (Line(1:4) == 'FROB') Go To 10310
if (Line(1:4) == 'NOMO') Go To 11000
if (Line(1:4) == 'TIME') Go To 12000
if (Line(1:4) == 'NOTI') Go To 12100
if (Line(1:4) == 'ERFU') Go To 13000
if (Line(1:4) == 'ORDE') Go To 14000
if (Line(1:4) == 'VIRT') Go To 15000
if (Line(1:4) == 'OCCU') Go To 15100
if (Line(1:4) == 'ALL ') Go To 15200
if (Line(1:4) == 'PAO ') Go To 16000
if (Line(1:4) == 'ANAP') Go To 16100
if (Line(1:4) == 'DOMA') Go To 17000
if (Line(1:4) == 'THRD') Go To 17100
if (Line(1:4) == 'THRP') Go To 17200
if (Line(1:4) == 'ANAD') Go To 17300
if (Line(1:4) == 'SKIP') Go To 18000
if (Line(1:4) == 'LOCN') Go To 18100
if (Line(1:4) == 'LOCC') Go To 18200
if (Line(1:4) == 'WAVE') Go To 18300
if (Line(1:4) == 'CONS') Go To 18400
if (Line(1:4) == 'FILE') Go To 18500
if (Line(1:4) == 'END ') Go To 99999
write(6,*) 'Unidentified key word  : ',Key
write(6,*) 'Internal representation: ',Line(1:4)
call FindErrorLine()
call Quit_OnUserError()

! DEBUg

1000 continue
Debug = .true.
Go To 999

! NORBitals

2000 continue
Line = Get_Ln(LuSpool)
call Get_I(1,nOrb2Loc,nSym)
nOrb2Loc_UsrDef = .true.
Go To 999

! NFROzen

2100 continue
Line = Get_Ln(LuSpool)
call Get_I(1,nFro,nSym)
nFro_UsrDef = .true.
Freeze = .false.
!if (LocVir) then
if (LocOrb == Virtual) then
  write(6,*)
  write(6,*) 'WARNING!!!'
  write(6,*) 'You have chosen to freeze some orbitals AND asked to localize the virtual space. This implies that'
  write(6,*) ' some virtual orbitals will be kept frozen, thus they will be left unchanged. Is that OK ?'
  write(6,*)
end if
Go To 999

! FREEze core orbitals (as defined by seward)

2200 continue
if (.not. nFro_UsrDef) then
  call Get_iArray('Non valence orbitals',nFro,nSym)
  Freeze = .true.
  !if (LocVir) then
  if (LocOrb == Virtual) then
    write(6,*)
    write(6,*) 'WARNING!!!'
    write(6,*) 'You have chosen to freeze some orbitals AND asked to localize the virtual space. This implies that'
    write(6,*) ' some virtual orbitals will be kept frozen, thus they will be left unchanged. Is that OK ?'
    write(6,*)
  end if
end if
Go To 999

! MAXImisation

3000 continue
Maximisation = .true.
Go To 999

! MINImisation

3100 continue
Maximisation = .false.
Go To 999

! NITErations or ITERations

4000 continue
Line = Get_Ln(LuSpool)
call Get_I1(1,NMxIter)
Go To 999

! THREshold

5000 continue
Line = Get_Ln(LuSpool)
call Get_F1(1,Thrs)
Thrs_UsrDef = .true.
Go To 999

! THRGrad

5001 continue
Line = Get_Ln(LuSpool)
call Get_F1(1,ThrGrad)
Go To 999

! CHOStart (use Cholesky orbitals as start guess for PM/Boys/ER)

5100 continue
ChoStart = .true.
Go To 999

! THRRotation

6000 continue
Line = Get_Ln(LuSpool)
call Get_F1(1,ThrRot)
Go To 999

! PIPEk-Mezey or PM

7000 continue
LocModel = 1
LocModel_UsrDef = .true.
Go To 999

! BOYS

7100 continue
LocModel = 2
LocModel_UsrDef = .true.
Go To 999

! CHOLesky

7200 continue
LocModel = 3
LocModel_UsrDef = .true.
Go To 999

! EDMIston-Ruedenberg or ER

7300 continue
LocModel = 4
LocModel_UsrDef = .true.
Go To 999

! SILEnt mode

8000 continue
Silent = .true.
Go To 999

! TEST localisation (orthonormality, density, etc.)

9000 continue
Test_Localisation = .true.
Go To 999

! ANALysis: generate bitmaps + histograms for density, original, and
!           local MOs. Analysis "per atom" is default, but does not work
!           with symmetry.

10000 continue
Analysis = .true.
AnaAtom = nSym == 1
Go To 999

! ANAAtom: generate bitmaps + histograms for density, original, and
!          local MOs. Analysis "per atom".

10100 continue
Analysis = .true.
AnaAtom = .true.
Go To 999

! ANAShell: generate bitmaps + histograms for density, original, and
!           local MOs. Analysis "per shell".

10200 continue
Analysis = .true.
AnaAtom = .false.
Go To 999

! MAX norm: use max element as norm in analysis.

10300 continue
AnaNrm = 'Max'
Go To 999

! FROBenius norm: use Frobenius norm in analysis.

10310 continue
AnaNrm = 'Fro'
Go To 999

! NOMO: do not print MOs.

11000 continue
PrintMOs = .false.
Go To 999

! TIME: time the localisation procedure explicitly

12000 continue
Timing = .true.
Go To 999

! NOTI: do not time the localisation procedure explicitly

12100 continue
Timing = .false.
Go To 999

! ERFUnctional: evaluate Edmiston-Ruedenberg functional with original
!               and localised orbitals.

13000 continue
EvalER = .true.
Go To 999

! ORDEr: order localised orbitals according to Cholesky orbital
!        ordering, defined according to the overlap between the two sets
!        orbitals.

14000 continue
Order = .true.
Go To 999

! VIRTual: localise virtual orbitals

15000 continue
LocOrb = Virtual
!LocVir = .true.
if (nFro_UsrDef .or. Freeze) then
  write(6,*)
  write(6,*) 'WARNING!!!'
  write(6,*) 'You have chosen to freeze some orbitals AND asked to localize the virtual space. This implies that'
  write(6,*) ' some virtual orbitals will be kept frozen, thus they will be left unchanged. Is that OK ?'
  write(6,*)
end if
Go To 999

! OCCUpied: localise occupied orbitals

15100 continue
LocOrb = Occupied
!LocVir = .false.
Go To 999

! ALL: localise all orbitals

15200 continue
LocOrb = All
Go To 999

! PAO : compute projected AOs that span the virtual space
!       (= all - occupied - frozen) using Cholesky decomposition to
!       remove linear dependence.

16000 continue
LocPAO = .true.
LocModel = 3
LocModel_UsrDef = .true.
Go To 999

! ANAP: Special analysis for Cholesky PAOs before orthonormalization.

16100 continue
AnaPAO = .true.
Go To 999

! DOMAin: set up orbital domains (Pulay-style) and pair domains

17000 continue
DoDomain = .true.
Go To 999

! THRDomain: thresholds for setting up orbital domains.
!            First value is the minimum sum of gross atomic Mulliken
!            charge for each orbital.
!            Second value is the threshold used for the completeness
!            check of Boughton and Pulay.

17100 continue
DoDomain = .true.
Line = Get_Ln(LuSpool)
call Get_F(1,ThrDomain,2)
Go To 999

! THRPairdomain: thresholds for setting up pair domains.
!                3 values needed: Rs, Rw, Rd (in bohr).
!                Let R be the minimum distance between any two atoms in
!                a pair domain. Then,
!                     R <= Rs: strong pair
!                Rs < R <= Rw: weak pair
!                Rw < R <= Rd: distant pair
!                Rd < R <= Rd: very distant pair

17200 continue
DoDomain = .true.
Line = Get_Ln(LuSpool)
call Get_F(1,ThrPairDomain,3)
Go To 999

! ANADomain: analysis of orbital domains

17300 continue
AnaDomain = .true.
DoDomain = .true.
Go To 999

! SKIP: skip localisation completely; i.e. only perform analysis etc.

18000 continue
Skip = .true.
Go To 999

! LOCN: localized natural orbitals

18100 continue
LocNatOrb = .true.
read(LuSpool,*) nActa,ThrSel
! Read atom names
18101 read(LuSpool,'(A)',end=9940) Line
if (Line(1:1) == '*') goto 18101
if (Line == Blank) goto 18101
call UpCase(Line)
do i=1,nActa
  call LeftAd(Line)
  if (Line == Blank) goto 9940
  j = index(Line,' ')
  NamAct(i) = Line(1:j-1)
  Line(1:j-1) = Blank(1:j-1)
end do
Go To 999

! LOCC: localized canonical orbitals

18200 continue
LocCanOrb = .true.
read(LuSpool,*) nActa,ThrSel
! Read atom names
read(LuSpool,'(A)',end=9940) Line
if (Line(1:1) == '*') goto 18101
if (Line == Blank) goto 18101
call UpCase(Line)
do i=1,nActa
  call LeftAd(Line)
  if (Line == Blank) goto 9940
  j = index(Line,' ')
  NamAct(i) = Line(1:j-1)
  Line(1:j-1) = Blank(1:j-1)
end do
Go To 999

! WAVE: wavelet transform of the MO basis

18300 continue
read(LuSpool,*) iWave
if ((iWave /= 0) .and. (iWave /= 1)) then
  write(6,*) ' WARNING!!!'
  write(6,*) ' Incorrect bit-switch input parameter for WAVE.'
  write(6,*) ' I will continue with the default value.'
  iWave = 0
end if
Skip = .false.
Wave = .true.
LocModel = 0
Go To 999

! CONS: constrained natural orbital analysis

18400 continue
Line = Get_Ln(LuSpool)
call Get_I(1,nConstr,nSym)
MxConstr = 0
do i=1,nSym
  MxConstr = max(MxConstr,nConstr(i))
end do
if (MxConstr > 16) then
  write(6,*) ' ERROR  !!!'
  write(6,*) ' Cannot handle more than 16 constraints/symm !!'
  write(6,*) ' Increase size of indxC in Get_CNOs and recompile.'
end if
DoCNOs = .true.
Skip = .false.
PrintMOs = .false.
LocModel = 0
Go To 999

! FILE: filename with input orbitals

18500 continue
Line = Get_Ln(LuSpool)
! Filename is read in localisation.f
Go To 999

! END of Input

99999 continue

! ==============
! Postprocessing
! ==============

! Only Cholesky localisation can be run with symmetry.
! Reset if the user explicitly requested a localisation
! model different from Cholesky.
! ------------------------------------------------------------------

if ((nSym > 1) .and. LocModel_UsrDef .and. (LocModel /= 3)) then
  write(6,*)
  write(6,*) 'WARNING!!!'
  write(6,*) 'The localisation model you have suggested in combination with symmetry will not work.'
  write(6,*) 'In this case the program defaults to the Cholesky localisation scheme!'
  write(6,*)
  LocModel = 3
end if

! Special settings for PAO: LocModel must be 3 (i.e. Cholesky), else
! we cancel PAO. For PAO, special analysis may be activated by the
! DEBUg keyword.
! ------------------------------------------------------------------

LocPAO = LocPAO .and. (LocModel == 3)
if (LocPAO) then
  AnaPAO = AnaPAO .or. Debug
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
  else ! All
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

if (Debug) then
  write(6,'(/,A,A)') SecNam,': orbital definitions:'
  write(6,'(A,8I9)') 'nBas    : ',(nBas(iSym),iSym=1,nSym)
  write(6,'(A,8I9)') 'nOrb    : ',(nOrb(iSym),iSym=1,nSym)
  write(6,'(A,8I9)') 'nOccInp : ',(nOccInp(iSym),iSym=1,nSym)
  write(6,'(A,8I9)') 'nVirInp : ',(nVirInp(iSym),iSym=1,nSym)
  write(6,'(A,8I9)') 'nFro    : ',(nFro(iSym),iSym=1,nSym)
  write(6,'(A,8I9,/)') 'nOrb2Loc: ',(nOrb2Loc(iSym),iSym=1,nSym)
end if

! If Cholesky, reset default threshold (unless user defined).
! -----------------------------------------------------------

if ((LocModel == 3) .and. (.not. Thrs_UsrDef)) then
  Thrs = 1.0d-8
end if

! No need to order Cholesky MOs.
! ------------------------------

Order = Order .and. (LocModel /= 3)

! No need to evaluate ER function for Edmiston-Ruedenberg.
! --------------------------------------------------------

EvalER = EvalER .and. (LocModel /= 4)

! Turn on localisation test if debug or analysis was specified.
! -------------------------------------------------------------

Test_Localisation = Test_Localisation .or. Debug .or. Analysis

Go To 9950

9940 continue
write(6,*) ' READIN: Premature end of file when reading selected'
write(6,*) ' atoms in keyword LOCN'
call Abend()

9950 continue

end subroutine Readinp_localisation
