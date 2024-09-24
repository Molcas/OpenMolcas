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

use InfSO, only: qNRTh, DltnTH, Energy, IterSO
use InfSCF, only: nAtoms, nOcc, nD, Thize, EThr, DThr, DelThr, DIISTh, FThr, QudThr, E1, E2, EKin, OnlyProp, NoProp, FckAuf, &
                  InVec, nIterP, Iter, iPrint, jPrint, iPrOrb, kIVO, iCoCo, jVOut, DIIS, Damping, One_Grid, Two_Thresholds, &
                  IDKeep, lPaper, nDens, kOptim, AccCon, nDisc, nCore, kDisk, PotNuc, EneV, E1V, E2V, iDisk, EkInv, MapDns, &
                  PreSch, MiniDn, WrOutD, C1DIIS, RSRFO, RGEK, Scrmbl, RFPert, pmTime, EmConv, AddFragments, rTemp, TemFac, TStop, &
                  KSDFT, ExFac, WarnCFG, WarnPOCC, WarnSlow, DoFMM, nIter, NamFld, LstVec, nBas, nDel, nFro, nOrb, nSym, TimFld, &
                  nFrz
use Constants, only: Zero, One
use MxDM, only: MxIter
use NDDO, only: twoel_NDDO
use Constants, only: Zero, One

implicit none
#include "hfc_logical.fh"
integer iFMM, iPrintLevel, nData
logical Found, Reduce_Prt
external Reduce_Prt

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
Thize = 1.0d-6
EThr = 1.0d-9
DThr = 1.0d-4
call Qpg_dScalar('S delete thr',Found)
if (Found) then
  call Get_dScalar('S delete thr',DelThr)
else
  DelThr = 1.0d-5
  call Put_dScalar('S delete thr',DelThr)
end if
DiisTh = 0.15d+00
QNRTh = 0.075d+00
!DltNTh = 0.2d-4
DltNTh = 0.1d-2
FThr = 1.5d-4
QudThr = 1.0d-5
Energy = Zero
E1 = Zero
E2 = Zero
EKin = Zero
PotNuc = Zero
EneV = Zero
E1V = Zero
E2V = Zero
EKinV = Zero
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
kDisk(:) = -1
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
NamFld(9) = '  d'') rotate MOs C with U      (QNR steps)  :'
NamFld(10) = '        U=exp(kap)                          :'
NamFld(11) = '          via Taylor expansion (sin/cos)    :'
NamFld(12) = '          via transformation to Schur basis :'
NamFld(13) = '  e'') transf. Fck Mat. with C  (QNR steps)  :'
NamFld(14) = '  f ) other calculations                    :'
NamFld(15) = '3) Final processing (generation of outputs) :'
NamFld(16) = '   T O T A L                                :'
twoel_NDDO = .false.
RTemp = 0.04d0
TemFac = 0.02d0
TStop = 0.00d0
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
