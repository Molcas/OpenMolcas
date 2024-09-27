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
                        nsymm,ngrid,dltT0,dltH0,zJ,tmin,tmax,hmin,hmax,XField,thrs,TempMagn,cryst,coord,encut_rate,gtens_input, &
                        D_fact,EoverD_fact,riso,decompose_exchange,AnisoLines1,AnisoLines3,AnisoLines9,DM_exchange,Dipol,KE, &
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
integer(kind=iwp), intent(out) :: neq(nneq), nexch(nneq), nK, mG, ncut, nP, AngPoints, nBlock, encut_definition, iopt, iPrint, &
                                  nsymm, ngrid
real(kind=wp), intent(out) :: dltT0, dltH0, zJ, tmin, tmax, hmin, hmax, Xfield, thrs, TempMagn(nTempMagn), cryst(6), coord(3), &
                              encut_rate, gtens_input(3,nneq), D_fact(nneq), EoverD_fact(nneq), riso(3,3,nneq)
logical(kind=iwp), intent(out) :: decompose_exchange, AnisoLines1, AnisoLines3, AnisoLines9, DM_exchange, Dipol, KE, &
                                  JITO_exchange, fitCHI, fitM, Do_structure_abc, DoPlot, compute_g_tensors, tinput, &
                                  compute_susceptibility, compute_torque, compute_barrier, compute_magnetization, hinput, &
                                  compute_Mdir_vector, zeeman_energy, m_paranoid, m_accurate, smagn
character, intent(out) :: itype(nneq)
integer(kind=iwp) :: i

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
  neq(:) = 1
  nexch(:) = 1
  itype(:) = ' '
  D_fact(:) = Zero
  EoverD_fact(:) = Zero
  gtens_input(:,:) = -gElectron
  do i=1,nneq
    call unitmat(riso(:,:,i),3)
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
Cryst(:) = Zero
coord(:) = Zero
!-----------------------------------------------------------------------
! mean field parameter
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
!-----------------------------------------------------------------------
! angular grid
nsymm = 1
ngrid = 15

return

end subroutine set_defaults
