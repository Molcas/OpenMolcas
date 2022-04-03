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

! Subroutine with purpose to initialize and set defaults for the input section.
subroutine Qmstat_Init()

use qmstat_global, only: Anal, ATitle, ChargedQM, Cordst, Cut_Elc, Cut_Ex1, Cut_Ex2, DelFi, DelOrAdd, DelR, DelX, Diel, DifSlExp, &
                         Disp, DispDamp, dLJRep, EdSt, EigV, EneLim, Exdt1, Exdtal, Exrep10, Exrep2, Exrep4, Exrep6, FieldDamp, &
                         FieldNuc, Forcek, iExtr_Atm, iExtra, iLuSaIn, iLuSaUt, iLuStIn, iLuStUt, iNrIn, iNrUt, iOrb, iPrint, &
                         iRead, iSeed, itMax, lCiSelect, lExtr, lMltSlC, lQuad, lSlater, MoAveRed, Mp2DensCorr, MxPut, nAtom, &
                         nCent, nCha, nEqState, nLvlShift, nMacro, nMicro, nPart, nPol, nSlSiteC, ParallelT, Pol, PolLim, Pres, &
                         Qmeq, QmProd, Qsta, RassiM, rStart, SaFilIn, SaFilUt, Sexre1, Sexre2, Sexrep, SimEx, SlExpC, SlFactC, &
                         SlPC, StFilIn, StFilUt, Surf, Temp
use stdalloc, only: mma_allocate
use Constants, only: Zero, One, Six, Ten, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: i, j

!IO_stuff
StFilIn = 'STFIL0'
SaFilIn = 'SAFIL0'
StFilUt = 'STFIL0'
SaFilUt = 'SAFIL0'
SimEx = 'EXTRA0'
!Jose: File for the optimization procedure
FieldNuc = 'AVENUC'
RassiM = 'RASSIM'
EigV = 'EIGV'
iNrIn = -1
iNrUt = 0
iLuStIn = 8+iNrIn
iLuStUt = 16+iNrUt
iLuSaIn = 24+iNrIn
iLuSaUt = 32+iNrUt
iRead = 0
! Defaults
nEqState = 1
Cut_Ex1 = Ten
Cut_Ex2 = Zero
DelX = Zero
DelFi = Zero
DelR = Zero
Temp = 300.0_wp
call GetSeed(ISEED)
iPrint = 1
NMACRO = 1
NMICRO = 1
RSTART = 80.0_wp
NPART = 0
nAtom = 3
NCENT = 5
NPOL = 3
NCHA = 4
!Jose.Slater Penetration
nSlSiteC = 5
lMltSlC = 0
!****
nLvlShift = 0
call mma_allocate(iExtr_Atm,0,label='iExtr_Atm')
call mma_allocate(Qsta,NCHA,label='Qsta')
QSTA(1) = 0.5836_wp
QSTA(2) = 0.5836_wp
QSTA(3) = -0.5836_wp
QSTA(4) = -0.5836_wp
call mma_allocate(Pol,NPOL,label='Pol')
POL(1) = 5.932_wp
POL(2) = 0.641_wp
POL(3) = 0.641_wp
!Jose.Slater Penetration
Cut_Elc = Six
DifSlExp = 1.0e-3_wp

call mma_allocate(SlFactC,4,nSlSiteC,label='SlFactC')
SlFactC(1,1) = -Half
SlFactC(1,2) = -0.4164_wp
SlFactC(1,3) = -0.4164_wp
SlFactC(1,4) = -0.5836_wp
SlFactC(1,5) = -0.5836_wp
SlFactC(2:4,:) = Zero
call mma_allocate(SlExpC,2,nSlSiteC,label='SlExpC')
SlExpC(1,1) = 2.5552_wp
SlExpC(1,2) = 2.6085_wp
SlExpC(1,3) = 2.6085_wp
SlExpC(1,4) = 2.5552_wp
SlExpC(1,5) = 2.5552_wp
SlExpC(2,:) = Zero
call mma_allocate(SlPC,nSlSiteC,label='SlPC')
SlPC(1) = Half
SlPC(2) = One
SlPC(3) = One
SlPC(4) = Zero
SlPC(5) = Zero
!******************************
Exrep2 = Zero
Exrep4 = Zero
Exrep6 = Zero
Exrep10 = Zero
call mma_allocate(sExRep,nAtom,nAtom,label='sExRep')
call mma_allocate(sExRe1,nAtom,nAtom,label='sExRe1')
call mma_allocate(sExRe2,nAtom,nAtom,label='sExRe2')
sExRep(1,1) = 2.092338_wp
sExRe1(1,1) = 158.998_wp
sExRe2(1,1) = 4.66009e10_wp
sExRep(2,1) = 2.112447_wp
sExRe1(2,1) = 8.31922_wp
sExRe2(2,1) = 97560.62_wp
sExRep(2,2) = 1.075803_wp
sExRe1(2,2) = 0.06521_wp
sExRe2(2,2) = 1121941276.0_wp
sExRep(3,1) = 2.112447_wp
sExRe1(3,1) = 8.31922_wp
sExRe2(3,1) = 97560.62_wp
sExRep(3,2) = 1.075803_wp
sExRe1(3,2) = 0.06521_wp
sExRe2(3,2) = 1121941276.0_wp
sExRep(3,3) = 1.075803_wp
sExRe1(3,3) = 0.06521_wp
sExRe2(3,3) = 1121941276.0_wp
do I=1,NATOM
  do J=1,I
    SEXREP(J,I) = SEXREP(I,J)
    SEXRE1(J,I) = SEXRE1(I,J)
    SEXRE2(J,I) = SEXRE2(I,J)
  end do
end do
call mma_allocate(Disp,NPOL,NPOL,label='Disp')
Disp(1,1) = 11.338_wp
Disp(2,1) = 3.38283_wp
Disp(2,2) = 0.627068_wp
Disp(3,1) = 3.38283_wp
Disp(3,2) = 0.627068_wp
Disp(3,3) = 0.627068_wp
do I=1,NPOL
  do J=1,I
    DISP(J,I) = DISP(I,J)
  end do
end do
call mma_allocate(Cordst,3,MxPut*nCent,label='Cordst')
CORDST(1,1) = Zero
CORDST(1,2) = Zero
CORDST(1,3) = Zero
CORDST(1,4) = 0.3126_wp
CORDST(1,5) = -0.3126_wp
CORDST(2,1) = Zero
CORDST(2,2) = 1.43_wp
CORDST(2,3) = -1.43_wp
CORDST(2,4) = Zero
CORDST(2,5) = Zero
CORDST(3,1) = 0.3_wp
CORDST(3,2) = -0.807_wp
CORDST(3,3) = -0.807_wp
CORDST(3,4) = -0.1191_wp
CORDST(3,5) = -0.1191_wp
ForceK = 1.0e-3_wp
dLJrep = Zero
Pres = One
PolLim = 1.0e-4_wp
EneLim = 1.0e-7_wp
itMax = 30
Exdtal = 30.0_wp
Exdt1 = 0.06_wp
Surf = 30.0_wp
iOrb(2) = 5
Diel = 80.0_wp
iExtra = 0
Qmeq = .false.
FieldDamp = .false.
DispDamp = .false.
QmProd = .false.
ChargedQM = .false.
ATitle = .false.
Anal = .false.
ParallelT = .false.
Mp2DensCorr = .false.
MoAveRed = .false.
lCiSelect = .false.
EdSt = .false.
DelOrAdd(:) = .false.
lExtr(:) = .false.
lSlater = .true.
lQuad = .false.

return

end subroutine Qmstat_Init
