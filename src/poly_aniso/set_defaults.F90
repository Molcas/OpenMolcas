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

subroutine set_defaults(nneq,nTempMagn,nDir,nDirZee,nMult,neq,nexch,nK,mG,ncut,nP,AngPoints,nBlock,encut_definition,iopt,iPrint, &
                        dltT0,dltH0,zJ,tmin,tmax,hmin,hmax,XField,thrs,TempMagn,cryst,coord,encut_rate,gtens_input,D_fact, &
                        EoverD_fact,riso,decompose_exchange,AnisoLines1,AnisoLines3,AnisoLines9,DM_exchange,Dipol,KE, &
                        JITO_exchange,fitCHI,fitM,Do_structure_abc,doplot,compute_g_tensors,tinput,compute_susceptibility, &
                        compute_torque,compute_barrier,compute_magnetization,hinput,compute_Mdir_vector,zeeman_energy,m_paranoid, &
                        m_accurate,smagn,itype)
! definition of the cluster
!  nneq, neq, gtens_input, D_fact, EoverD_fact, riso, itype
! definition of exchange interaction
!  nexch
! definition of g and D tensors
!  nMult, compute_g_tensors
! Zeeman energy and M vector
!  nDir, nDirZee
! definition of the exchange:
!  AnisoLines1, AnisoLines3, AnisoLines9
! options used in connection with Dipol-Dipol interaction
!  Dipol
! options used in connection with KE
!  KE
! option for Dzialoshinsky-Morya antisymmetric exchange
!  DM_exchange
! option for ITO exchange
!  JITO_exchange
! options for automatic fitting of parameters:
!  fitCHI : not used so far
!  fitM   : not used so far
! definition of data for susceptibility
!  tinput, compute_susceptibility, tmin, tmax, dltT0
! options related to XT_MoverH
!  Xfield
! definition of data for magnetization:
!  nTempMagn, TempMagn, iopt, dltH0, thrs, hmin, hmax, hinput, compute_magnetization, compute_Mdir_vector, zeeman_energy,
!  m_paranoid, m_accurate, smagn
! options used to set up nM and EM
!  nK, mG     : encut_definition=1;
!  ncut       : encut_definition=2;
!  encut_rate : encut_definition=3;
!  encut_definition
! decompose exchange
!  decompose_exchange
! magnetization torque
!  compute_torque, nP, AngPoints
! mean field parameter
!  zJ
! definition of the crystal axes:
!  cryst : a, b, c, alpha, beta, gamma
!  Do_structure_abc
! Cartesian coordinates of the main metal site, or center
!  coord
! definitions for blocking barrier
!  nBlock, compute_barrier
! default print level
!  iPrint, DoPlot

use Constants, only: Zero, Two, Ten, gElectron
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nneq, nTempMagn, nDir, nDirZee, nMult
integer(kind=iwp), intent(out) :: neq(nneq), nexch(nneq), nK, mG, ncut, nP, AngPoints, nBlock, encut_definition, iopt, iPrint
real(kind=wp), intent(out) :: dltT0, dltH0, zJ, tmin, tmax, hmin, hmax, Xfield, thrs, TempMagn(nTempMagn), cryst(6), coord(3), &
                              encut_rate, gtens_input(3,nneq), D_fact(nneq), EoverD_fact(nneq), riso(nneq,3,3)
logical(kind=iwp), intent(out) :: decompose_exchange, AnisoLines1, AnisoLines3, AnisoLines9, DM_exchange, Dipol, KE, &
                                  JITO_exchange, fitCHI, fitM, Do_structure_abc, DoPlot, compute_g_tensors, tinput, &
                                  compute_susceptibility, compute_torque, compute_barrier, compute_magnetization, hinput, &
                                  compute_Mdir_vector, zeeman_energy, m_paranoid, m_accurate, smagn
character, intent(out) :: itype(nneq)
#include "mgrid.fh"
integer(kind=iwp) :: i, j, l

!------------------------------------------------------------------
! at this point, the follwing variables have been already assigned
! their values
!integer :: nneq, exch, nLoc, nCenter, nT, nH, nTempMagn, nDir, nDirZee, nMult, nPair
!
!ntotcenter = 1
!local_states = 1
!nDirZee = 0
!nCenter = 1
!exch = 1
!nT = 301
! --  in this subroutine we define other default values of other variables
! --  these values will be (partly) overwritten by the subsequent
!     "readin_poly" subroutine

!-----------------------------------------------------------------------
! variables related to the cluster:
if (nneq > 0) then
  do i=1,nneq
    neq(i) = 1
    nexch(i) = 1
    itype(i) = ' '
    D_fact(i) = Zero
    EoverD_fact(i) = Zero
    call unitmat(riso(i,:,:),3)
    do l=1,3
      gtens_input(l,i) = -gElectron
    end do
  end do
end if
!-----------------------------------------------------------------------
! magnetic exchange
Dipol = .false.
AnisoLines1 = .false.
AnisoLines3 = .false.
AnisoLines9 = .false.
KE = .false.
DM_exchange = .false.
decompose_exchange = .false.
JITO_exchange = .false.
!-----------------------------------------------------------------------
! g and D tensors
compute_g_tensors = nMult > 0
!-----------------------------------------------------------------------
! Magnetization vector and Zeeman energy spliting
compute_Mdir_vector = nDir > 0
zeeman_energy = nDirZee > 0
!-----------------------------------------------------------------------
! powder magnetization:
compute_magnetization = .false.
m_accurate = .true.
m_paranoid = .true.
smagn = .false.
THRS = 1.0e-10_wp ! threshold for convergence of average spin when zJ /= 0.0
iopt = 1
dltH0 = 0.001_wp ! the non-zero field point
NK = 50
MG = 50
HMIN = Zero
HMAX = Ten
encut_definition = 2
encut_rate = 1.0_wp
if (nTempMagn > 0) TempMagn(1) = Two
ncut = 0
!-----------------------------------------------------------------------
! magnetic susceptibility
compute_susceptibility = .false.
TMIN = Zero
TMAX = 300.0_wp
DLTT0 = 0.001_wp
Xfield = Zero
!-----------------------------------------------------------------------
! magnetic torque
compute_torque = .false.
nP = 0
AngPoints = 0
!-----------------------------------------------------------------------
! perform comparison with experiment
TINPUT = .false.
HINPUT = .false.
!-----------------------------------------------------------------------
! relate to crystallographic axes
Do_structure_abc = .false.
do i=1,6
  Cryst(i) = Zero
end do
do i=1,3
  coord(i) = Zero
end do
!-----------------------------------------------------------------------
! mean field paraemeter
ZJ = Zero
!-----------------------------------------------------------------------
! automatic fitting:
fitCHI = .false.
fitM = .false.
!-----------------------------------------------------------------------
! print level
iPrint = 2
!-----------------------------------------------------------------------
! blocking barrier
nBlock = 0
compute_barrier = .false.
!ccc
DoPlot = .false.
!------------------------------------------------------------------
! variables in "mgrid.fh"
do i=1,32
  do j=1,3
    get_nP(j,i) = 0
  end do
end do
nsymm = 1
ngrid = 15
get_nP(1,1) = 5
get_nP(1,2) = 9
get_nP(1,3) = 17
get_nP(1,4) = 25
get_nP(1,5) = 29
get_nP(1,6) = 45
get_nP(1,7) = 49
get_nP(1,8) = 61
get_nP(1,9) = 77
get_nP(1,10) = 93
get_nP(1,11) = 105
get_nP(1,12) = 125
get_nP(1,13) = 141
get_nP(1,14) = 161
get_nP(1,15) = 185
get_nP(1,16) = 229
get_nP(1,17) = 309
get_nP(1,18) = 401
get_nP(1,19) = 505
get_nP(1,20) = 621
get_nP(1,21) = 749
get_nP(1,22) = 889
get_nP(1,23) = 1041
get_nP(1,24) = 1205
get_nP(1,25) = 1381
get_nP(1,26) = 1569
get_nP(1,27) = 1769
get_nP(1,28) = 1981
get_nP(1,29) = 2205
get_nP(1,30) = 2441
get_nP(1,31) = 2689
get_nP(1,32) = 2949
get_nP(2,1) = 4
get_nP(2,2) = 6
get_nP(2,3) = 11
get_nP(2,4) = 16
get_nP(2,5) = 17
get_nP(2,6) = 27
get_nP(2,7) = 28
get_nP(2,8) = 34
get_nP(2,9) = 41
get_nP(2,10) = 51
get_nP(2,11) = 57
get_nP(2,12) = 68
get_nP(2,13) = 75
get_nP(2,14) = 86
get_nP(2,15) = 98
get_nP(2,16) = 121
get_nP(2,17) = 162
get_nP(2,18) = 209
get_nP(2,19) = 262
get_nP(2,20) = 321
get_nP(2,21) = 386
get_nP(2,22) = 457
get_nP(2,23) = 534
get_nP(2,24) = 617
get_nP(2,25) = 706
get_nP(2,26) = 801
get_nP(2,27) = 902
get_nP(2,28) = 1009
get_nP(2,29) = 1122
get_nP(2,30) = 1241
get_nP(2,31) = 1366
get_nP(2,32) = 1497
get_nP(3,1) = 3
get_nP(3,2) = 4
get_nP(3,3) = 7
get_nP(3,4) = 10
get_nP(3,5) = 10
get_nP(3,6) = 16
get_nP(3,7) = 16
get_nP(3,8) = 19
get_nP(3,9) = 22
get_nP(3,10) = 28
get_nP(3,11) = 31
get_nP(3,12) = 37
get_nP(3,13) = 40
get_nP(3,14) = 46
get_nP(3,15) = 52
get_nP(3,16) = 64
get_nP(3,17) = 85
get_nP(3,18) = 109
get_nP(3,19) = 136
get_nP(3,20) = 166
get_nP(3,21) = 199
get_nP(3,22) = 235
get_nP(3,23) = 274
get_nP(3,24) = 316
get_nP(3,25) = 361
get_nP(3,26) = 409
get_nP(3,27) = 460
get_nP(3,28) = 514
get_nP(3,29) = 571
get_nP(3,30) = 631
get_nP(3,31) = 694
get_nP(3,32) = 760

return

end subroutine set_defaults
