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

use Constants, only: Zero, One, Six, Ten, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: i, j
#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"

!IO_stuff
StFilIn = 'STFIL0'
SAFilIn = 'SAFIL0'
StFilUt = 'STFIL0'
SAFilUt = 'SAFIL0'
BlockIn = 'BLOCIN'
BlockUt = 'BLOCUT'
SimEx = 'EXTRA0'
!Jose: File for the optimization procedure
FieldNuc = 'AVENUC'
do i=1,MxJobs
  write(JbName(i),'(A,i3.3)') 'JOB',i
end do
RassiM = 'RASSIM'
GammaO = 'GAMORB'
EigV = 'EIGV'
AddOns(1) = 'ADDON1'
AddOns(2) = 'ADDON2'
AddOns(3) = 'ADDON3'
iNrIn = -1
iNrUt = 0
iLuStIn = 8+iNrIn
iLuStUt = 16+iNrUt
iLuSaIn = 24+iNrIn
iLuSaUt = 32+iNrUt
iLuBlockIn = 3
iLuBlockUt = 4
iRead = 0
! Defaults
nEqState = 1
Cut_Ex1 = Ten
Cut_Ex2 = Zero
DelX = Zero
DelFi = Zero
DelR = Zero
Temp = 300.0_wp
ISEED = 791204
IPrint = 1
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
do i=1,MxAt
  iExtr_Atm(i) = -1
end do
QSTA(1) = 0.5836_wp
QSTA(2) = 0.5836_wp
QSTA(3) = -0.5836_wp
QSTA(4) = -0.5836_wp
POL(1) = 5.932_wp
POL(2) = 0.641_wp
POL(3) = 0.641_wp
!Jose.Slater Penetration
Cut_Elc = Six
DifSlExp = 1.0e-3_wp

SlFactC(1,1) = -Half
SlFactC(1,2) = -0.4164_wp
SlFactC(1,3) = -0.4164_wp
SlFactC(1,4) = -0.5836_wp
SlFactC(1,5) = -0.5836_wp
do i=1,5
  do j=2,4
    SlFactC(j,i) = Zero
  end do
end do
SlExpC(1,1) = 2.5552_wp
SlExpC(1,2) = 2.6085_wp
SlExpC(1,3) = 2.6085_wp
SlExpC(1,4) = 2.5552_wp
SlExpC(1,5) = 2.5552_wp
do i=1,5
  SlExpC(2,i) = Zero
end do
SlPC(1) = Half
SlPC(2) = One
SlPC(3) = One
SlPC(4) = Zero
SlPC(5) = 0.0d0
!******************************
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
Disp(1,1) = 11.338_wp
Disp(2,1) = 3.38283_wp
Disp(2,2) = 0.627068_wp
Disp(3,1) = 3.38283_wp
Disp(3,2) = 0.627068_wp
Disp(3,3) = 0.627068_wp
do I=1,NPOL
  do J=1,I
    DISP(J,I) = DISP(I,J)
    SEXREP(J,I) = SEXREP(I,J)
    SEXRE1(J,I) = SEXRE1(I,J)
    SEXRE2(J,I) = SEXRE2(I,J)
  end do
end do
CORDST(1,1) = Zero
CORDST(2,1) = Zero
CORDST(3,1) = Zero
CORDST(4,1) = 0.3126_wp
CORDST(5,1) = -0.3126_wp
CORDST(1,2) = Zero
CORDST(2,2) = 1.43_wp
CORDST(3,2) = -1.43_wp
CORDST(4,2) = Zero
CORDST(5,2) = Zero
CORDST(1,3) = 0.3_wp
CORDST(2,3) = -0.807_wp
CORDST(3,3) = -0.807_wp
CORDST(4,3) = -0.1191_wp
CORDST(5,3) = -0.1191_wp
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
Smeq = .false.
Qmeq = .false.
Fielddamp = .false.
Dispdamp = .false.
Smprod = .false.
QmProd = .false.
ChargedQM = .false.
ATitle = .false.
Anal = .false.
ParallelT = .false.
Mp2DensCorr = .false.
MoAveRed = .false.
lCiSelect = .false.
EdSt = .false.
!JoseMEP***** The dimension was increased from 8 to 12
do i=1,12
  DelOrAdd(i) = .false.
  lExtr(i) = .false.
  lAnal(i) = .false.
end do
lSlater = .true.
lQuad = .false.

return

end subroutine Qmstat_Init
