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
!                                                                      *
!***********************************************************************

subroutine SCF_Init()
!***********************************************************************
!                                                                      *
!     purpose: set up parameters that has to be predefined in SCF      *
!                                                                      *
!***********************************************************************

use InfSCF, only: AccCon, AddFragments, C1DIIS, Damping, DelThr, DIIS, DIISTh, DltnTH, DoFMM, DThr, E1V, E2V, EKin, EmConv, &
                  Energy, EneV, EThr, ExFac, FckAuf, FThr, iCoCo, iDisk, IDKeep, InVec, iPrint, iPrOrb, Iter, IterSO, jPrint, &
                  jVOut, kIVO, kOptim, KSDFT, lPaper, LstVec, MapDns, MiniDn, MxIter, NamFld, nAtoms, nBas, nCore, nD, nDel, &
                  nDens, nDisc, nFro, nFrz, nIter, nIterP, nOcc, NoProp, nOrb, nSym, One_Grid, OnlyProp, pmTime, PotNuc, PreSch, &
                  qNRTh, QudThr, RFPert, RGEK, RSRFO, rTemp, Scrmbl, TemFac, Thize, TimFld, TStop, Two_Thresholds, WarnCFG, &
                  WarnPOCC, WarnSlow, WrOutD
use NDDO, only: twoel_NDDO
use hfc_logical, only: UHF_HFC
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iFMM, iPrintLevel, nData
logical(kind=iwp) :: Found
logical(kind=iwp), external :: Reduce_Prt

!----------------------------------------------------------------------*
!     Initialize global variables                                      *
!----------------------------------------------------------------------*
call Peek_iScalar('nSym',nSym)
nAtoms = 0
nBas(:) = 0
nOrb(:) = 0
nOcc(:,:) = 0
if (nD == 1) then
  call Put_iArray('nIsh',nOcc(1,1),nSym)
else
  call Put_iArray('nIsh',nOcc(1,1),nSym)
  call Put_iArray('nIsh beta',nOcc(1,2),nSym)
end if
nFro(:) = 0
call Put_iArray('nFro',nFro,nSym)
nFrz(:) = 0
nDel(:) = 0
call qpg_iarray('nDel',Found,ndata)
if (.not. Found) call Put_iArray('nDel',nDel,nSym)
Thize = 1.0e-6_wp
EThr = 1.0e-9_wp
DThr = 1.0e-4_wp
call Qpg_dScalar('S delete thr',Found)
if (Found) then
  call Get_dScalar('S delete thr',DelThr)
else
  DelThr = 1.0e-5_wp
  call Put_dScalar('S delete thr',DelThr)
end if
DiisTh = 0.15_wp
QNRTh = 0.075_wp
!DltNTh = 0.2e-4_wp
DltNTh = 0.1e-2_wp
FThr = 1.5e-4_wp
QudThr = 1.0e-5_wp
Energy = Zero
EKin = Zero
PotNuc = Zero
EneV = Zero
E1V = Zero
E2V = Zero
OnlyProp = .false.
NoProp = .false.
FckAuf = .true.

! Order: SCF(0), Guessorb(1), Lumorb(2), Core(4)
! New order: SCF(0), Lumorb(2), Guessorb(1), Core(4)

InVec = 1
LstVec(:) = -1
LstVec(1) = 0
LstVec(2) = 2
LstVec(3) = 1
LstVec(4) = 4
nIter(0) = MxIter
nIter(1) = MxIter
nIterP = 1
iter = 1
iterso = 0

iPrint = iPrintLevel(-1)
jPrint = iPrint
if (Reduce_Prt() .and. (iPrint < 3)) jPrint = 0

iPrOrb = 1
kIvo = 0
nD = 1
UHF_HFC = .false.
iCoCo = 0
jVOut = 2
Diis = .true.
Damping = .true.
One_Grid = .false.
Two_Thresholds = .true.
iDKeep = 4
lPaper = 132
nDens = 0
kOptim = 1
AccCon = '         '
nDisc = 2000
nCore = 512
iDisk(:,1) = -1
MapDns(:) = 0
PreSch = .false.
MiniDn = .true.
WrOutD = .false.
c1Diis = .false.
RSRFO = .false.
RGEK = .false.
Scrmbl = .false.
RFpert = .false.
PmTime = .false.
EmConv = .false.
AddFragments = .false.
TimFld(:) = Zero
NamFld(1) = '1) Input processing                         :'
NamFld(2) = '2) Wave function optimization               :'
NamFld(3) = '     Line Search Iterations    (QNR steps)  :'
NamFld(4) = '  a ) calculation of the density            :'
NamFld(5) = '  b ) contraction with integrals            :'
NamFld(6) = '  c ) acceleration of convergence           :'
NamFld(7) = '        recursive BFGS         (QNR steps)  :'
NamFld(8) = '  d ) solution to Roothaan-Hall equations   :'
NamFld(9) = '  e ) rotate MOs C with U       (QNR steps) :'
NamFld(10) = '        U=exp(kap)                          :'
NamFld(11) = '  f ) transformation to new reference       :'
NamFld(12) = '  g ) s-GEK/RVO microiterations             :'
NamFld(13) = '  h ) transf. Fck Mat. with C  (QNR steps)  :'
NamFld(14) = '  i ) other calculations                    :'
NamFld(15) = '3) Final processing (generation of outputs) :'
NamFld(16) = '   T O T A L                                :'
twoel_NDDO = .false.
RTemp = 0.04_wp
TemFac = 0.02_wp
TStop = Zero
KSDFT = 'SCF '
ExFac = One
WarnCfg = .false.
WarnPocc = .false.
WarnSlow = .false.

call Get_iScalar('FMM',iFMM)
DoFMM = iFMM == 1

! Initialize energy due to external potential on the run file. This
! to make sure that it can be read unconditionally and is zero if
! not redefined by DrvXV
call Poke_dScalar('KSDFT energy',Zero)
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

end subroutine SCF_Init
