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
!***********************************************************************
!                                                                      *
!-----Data for reaction field calculations.                            *
!                                                                      *
!     lMax: highest angular momentum for multipole expansion           *
!     Eps  : dielectric constant                                       *
!     EpsInf: refraction index                                         *
!     rds: radius of cavity                                            *
!     latato: Gitter type                                              *
!     polsi: site polarizability                                       *
!     dipsi: site dipole moment                                        *
!     radlat: maximum extension of the lattice                         *
!     scala,scalb,scalc: lengths of the cell dimensions                *
!     scaaa: overall scaling                                           *
!     gatom: atoms in the lattice.                                     *
!     diedel: controls deletion of gitter polarizabilities             *
!     tK: Boltzmann factor                                             *
!     clim: convergence threshold for solver                           *
!     afac: controls equation solver (valid values 0.1-0.97)           *
!     nexpo: exponent of the potential with which the QM system kills  *
!            the gitter polarizabilities                               *
!     prefac: scaling of the polarization energies                     *
!     fmax: the square of the largest field in any grid point          *
!                                                                      *
!     PCM: whether to use the PCM solvent model                        *
!     Conductor: whether to use the Conductor-PCM model                *
!     Solvent: the name of the solvent we are using                    *
!             (allowed solvents in datasol.f)                          *
!     ISlPar: 100 integers to pass quickly PCM information             *
!            (defaulted and explained in pcmdef.f)                     *
!     RSlPar: 100 reals to pass quickly PCM information                *
!            (defaulted and explained in pcmdef.f)                     *
!     MxA: maximum number of atoms for the building of PCM cavity      *
!     NSinit: initial number of atomic spheres                         *
!     NS: actual number of spheres (initial+smoothing)                 *
!     nTs: number of surface tesserae                                  *
!     NOrdInp: number of atoms where explicit spheres are centered     *
!     RadInp: radius of spheres explicitly given in the input          *
!                                                                      *
!***********************************************************************

module Rctfld_Module

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

integer, parameter :: MxA = 1000, MxPar = 100
integer(kind=iwp) :: iCharge_ref, ISlPar(MxPar), latato, lMax = -1, maxa, maxb, maxc, nCavxyz, nexpo, nGrid, nGrid_Eff, &
                     nGridAverage, nGridSeed, NOrdInp(MxA), nPCM_Info, NS, NSinit, nSparse = 0, nTs
real(kind=wp) :: afac, clim = Zero, Cordsi(3,4) = Zero, dampIter = Zero, diedel, dipCutoff = Zero, dipsi = Zero, &
                 distSparse = Zero, Eps, Eps_User, EpsInf, EpsInf_User, fmax, gatom, polsi = Zero, prefac, RadInp(MxA), &
                 radlat = Zero, rds, rotAlpha = Zero, rotBeta = Zero, rotGamma = Zero, rsca, RSlPar(MxPar), RSolv, scaaa, scal14, &
                 scala = Zero, scalb, scalc, tK = Zero, VMol
logical(kind=iwp) :: Conductor, DoDeriv, Done_Lattice, lAmberPol, lDamping, lDiprestart, lFirstIter, LGridAverage, lLangevin, lRF, &
                     lRFCav, LSparse, NonEq_ref, PCM
character(len=32) :: Solvent
real(kind=wp), allocatable :: MM(:,:)

public :: aFac, cLim, Conductor, Cordsi, DampIter, DieDel, DipCutOff, Dipsi, DistSparse, DoDeriv, Done_Lattice, Eps, Eps_USER, &
          EpsInf, EpsInf_USER, fMax, gAtom, iCharge_Ref, ISlPar, lAmberPol, LatAto, lDamping, lDipRestart, lFirstIter, &
          lGridAverage, lLangevin, lMax, lRF, lRFCav, lSparse, MaxA, MaxB, MaxC, MM, MXA, nCavxyz, nExpO, nGrid, nGrid_Eff, &
          nGridAverage, nGridSeed, NonEQ_Ref, nOrdInp, nPCM_Info, nS, nSInit, nSparse, nTS, PCM, PCM_Info_Dmp, PCM_Info_Get, &
          Polsi, PreFac, RadInp, RadLat, RDS, RotAlpha, RotBeta, RotGamma, rSca, RSlPar, RSolv, ScaAA, Scal14, Scala, ScalB, &
          Scalc, Solvent, TK, vMol

contains

subroutine PCM_Info_Dmp()

  integer(kind=iwp) :: i, l_cRF, l_iRF, l_bRF, l_rRF
  integer(kind=iwp), allocatable :: bRF(:), iRF(:)
  real(kind=wp), allocatable :: rRF(:)
  character(len=:), allocatable :: cRF

  l_iRF = 17+size(ISlPar)+size(NOrdInp)
  call mma_allocate(iRF,l_iRF,label='iRF')
  i = 1
  iRF(i) = lMax
  i = i+1
  iRF(i) = latato
  i = i+1
  iRF(i) = nexpo
  i = i+1
  iRF(i) = maxa
  i = i+1
  iRF(i) = maxb
  i = i+1
  iRF(i) = maxc
  i = i+1
  iRF(i) = nCavxyz
  i = i+1
  iRF(i) = nGrid
  i = i+1
  iRF(i) = nGrid_Eff
  i = i+1
  iRF(i) = nSparse
  i = i+1
  iRF(i:i+size(ISlPar)-1) = ISlPar(:)
  i = i+size(ISlPar)
  iRF(i) = NSinit
  i = i+1
  iRF(i) = NS
  i = i+1
  iRF(i) = nTs
  i = i+1
  iRF(i:i+size(NOrdInp)-1) = NOrdInp(:)
  i = i+size(NOrdInp)
  iRF(i) = nPCM_Info
  i = i+1
  iRF(i) = iCharge_ref
  i = i+1
  iRF(i) = nGridAverage
  i = i+1
  iRF(i) = nGridSeed
  call Put_iArray('RFiInfo',iRF,l_iRF)
  call mma_deallocate(iRF)

  l_rRF = 29+size(Cordsi)+size(RSlPar)+size(RadInp)
  call mma_allocate(rRF,l_rRF,label='rRF')
  i = 1
  rRF(i) = EpsInf_User
  i = i+1
  rRF(i) = Eps_User
  i = i+1
  rRF(i) = rds
  i = i+1
  rRF(i:i+size(Cordsi)-1) = pack(Cordsi,.true.)
  i = i+size(Cordsi)
  rRF(i) = polsi
  i = i+1
  rRF(i) = dipsi
  i = i+1
  rRF(i) = radlat
  i = i+1
  rRF(i) = scala
  i = i+1
  rRF(i) = scalb
  i = i+1
  rRF(i) = scalc
  i = i+1
  rRF(i) = scaaa
  i = i+1
  rRF(i) = gatom
  i = i+1
  rRF(i) = diedel
  i = i+1
  rRF(i) = tK
  i = i+1
  rRF(i) = rotAlpha
  i = i+1
  rRF(i) = rotBeta
  i = i+1
  rRF(i) = rotGamma
  i = i+1
  rRF(i) = distSparse
  i = i+1
  rRF(i) = clim
  i = i+1
  rRF(i) = afac
  i = i+1
  rRF(i) = prefac
  i = i+1
  rRF(i) = fmax
  i = i+1
  rRF(i) = rsca
  i = i+1
  rRF(i) = Eps
  i = i+1
  rRF(i) = EpsInf
  i = i+1
  rRF(i) = RSolv
  i = i+1
  rRF(i) = VMol
  i = i+1
  rRF(i:i+size(RadInp)-1) = RadInp(:)
  i = i+size(RadInp)
  rRF(i:i+size(RSlPar)-1) = RSlPar(:)
  i = i+size(RSlPar)
  rRF(i) = dampIter
  i = i+1
  rRF(i) = dipCutoff
  i = i+1
  rRF(i) = scal14
  call Put_dArray('RFrInfo',rRF,l_rRF)
  call mma_deallocate(rRF)

  l_bRF = 14
  call mma_allocate(bRF,l_bRF,label='bRF')
  i = 1
  bRF(i) = merge(1,0,lRF)
  i = i+1
  bRF(i) = merge(1,0,lLangevin)
  i = i+1
  bRF(i) = merge(1,0,PCM)
  i = i+1
  bRF(i) = merge(1,0,Conductor)
  i = i+1
  bRF(i) = merge(1,0,NonEq_ref)
  i = i+1
  bRF(i) = merge(1,0,DoDeriv)
  i = i+1
  bRF(i) = merge(1,0,lRFCav)
  i = i+1
  bRF(i) = merge(1,0,LSparse)
  i = i+1
  bRF(i) = merge(1,0,LGridAverage)
  i = i+1
  bRF(i) = merge(1,0,lDamping)
  i = i+1
  bRF(i) = merge(1,0,lAmberPol)
  i = i+1
  bRF(i) = merge(1,0,Done_Lattice)
  i = i+1
  bRF(i) = merge(1,0,lFirstIter)
  i = i+1
  bRF(i) = merge(1,0,lDiprestart)
  call Put_iArray('RFlInfo',bRF,l_bRF)
  call mma_deallocate(bRF)

  l_cRF = len(Solvent)
  call mma_allocate(cRF,l_cRF,label='cRF')
  i = 1
  cRF(i:i+len(Solvent)-1) = Solvent
  call Put_cArray('RFcInfo',cRF,l_cRF)
  call mma_deallocate(cRF)

end subroutine PCM_Info_Dmp

subroutine PCM_Info_Get()

  integer(kind=iwp) :: i, l_cRF, l_iRF, l_bRF, l_rRF
  integer(kind=iwp), allocatable :: bRF(:), iRF(:)
  real(kind=wp), allocatable :: rRF(:)
  character(len=:), allocatable :: cRF

  l_iRF = 17+size(ISlPar)+size(NOrdInp)
  call mma_allocate(iRF,l_iRF,label='iRF')
  call Get_iArray('RFiInfo',iRF,l_iRF)
  i = 1
  lMax = iRF(i)
  i = i+1
  latato = iRF(i)
  i = i+1
  nexpo = iRF(i)
  i = i+1
  maxa = iRF(i)
  i = i+1
  maxb = iRF(i)
  i = i+1
  maxc = iRF(i)
  i = i+1
  nCavxyz = iRF(i)
  i = i+1
  nGrid = iRF(i)
  i = i+1
  nGrid_Eff = iRF(i)
  i = i+1
  nSparse = iRF(i)
  i = i+1
  ISlPar(:) = iRF(i:i+size(ISlPar)-1)
  i = i+size(ISlPar)
  NSInit = iRF(i)
  i = i+1
  NS = iRF(i)
  i = i+1
  nTs = iRF(i)
  i = i+1
  NOrdInp = iRF(i:i+size(NOrdInp)-1)
  i = i+size(NOrdInp)
  nPCM_Info = iRF(i)
  i = i+1
  iCharge_ref = iRF(i)
  i = i+1
  nGridAverage = iRF(i)
  i = i+1
  nGridSeed = iRF(i)
  call mma_deallocate(iRF)

  l_rRF = 29+size(Cordsi)+size(RSlPar)+size(RadInp)
  call mma_allocate(rRF,l_rRF,label='rRF')
  call Get_dArray('RFrInfo',rRF,l_rRF)
  i = 1
  EpsInf_User = rRF(i)
  i = i+1
  Eps_User = rRF(i)
  i = i+1
  rds = rRF(i)
  i = i+1
  Cordsi(:,:) = reshape(rRF(i:i+size(Cordsi)-1), shape(Cordsi))
  i = i+size(Cordsi)
  polsi = rRF(i)
  i = i+1
  dipsi = rRF(i)
  i = i+1
  radlat = rRF(i)
  i = i+1
  scala = rRF(i)
  i = i+1
  scalb = rRF(i)
  i = i+1
  scalc = rRF(i)
  i = i+1
  scaaa = rRF(i)
  i = i+1
  gatom = rRF(i)
  i = i+1
  diedel = rRF(i)
  i = i+1
  tK = rRF(i)
  i = i+1
  rotAlpha = rRF(i)
  i = i+1
  rotBeta = rRF(i)
  i = i+1
  rotGamma = rRF(i)
  i = i+1
  distSparse = rRF(i)
  i = i+1
  clim = rRF(i)
  i = i+1
  afac = rRF(i)
  i = i+1
  prefac = rRF(i)
  i = i+1
  fmax = rRF(i)
  i = i+1
  rsca = rRF(i)
  i = i+1
  Eps = rRF(i)
  i = i+1
  EpsInf = rRF(i)
  i = i+1
  RSolv = rRF(i)
  i = i+1
  VMol = rRF(i)
  i = i+1
  RadInp(:) = rRF(i:i+size(RadInp)-1)
  i = i+size(RadInp)
  RSlPar(:) = rRF(i:i+size(RSlPar)-1)
  i = i+size(RSlPar)
  dampIter = rRF(i)
  i = i+1
  dipCutoff = rRF(i)
  i = i+1
  scal14 = rRF(i)
  call mma_deallocate(rRF)


  l_bRF = 14
  call mma_allocate(bRF,l_bRF,label='bRF')
  call Get_iArray('RFlInfo',bRF,l_bRF)
  i = 1
  lRF = bRF(i) > 0
  i = i+1
  lLangevin = bRF(i) > 0
  i = i+1
  PCM = bRF(i) > 0
  i = i+1
  Conductor = bRF(i) > 0
  i = i+1
  NonEq_ref = bRF(i) > 0
  i = i+1
  DoDeriv = bRF(i) > 0
  i = i+1
  lRFCav = bRF(i) > 0
  i = i+1
  LSparse = bRF(i) > 0
  i = i+1
  LGridAverage = bRF(i) > 0
  i = i+1
  lDamping = bRF(i) > 0
  i = i+1
  lAmberPol = bRF(i) > 0
  i = i+1
  Done_Lattice = bRF(i) > 0
  i = i+1
  lFirstIter = bRF(i) > 0
  i = i+1
  lDiprestart = bRF(i) > 0
  call mma_deallocate(bRF)

  l_cRF = len(Solvent)
  call mma_allocate(cRF,l_cRF,label='cRF')
  call Get_cArray('RFcInfo',cRF,l_cRF)
  i = 1
  Solvent = cRF(i:i+len(Solvent)-1)
  call mma_deallocate(cRF)

end subroutine PCM_Info_Get

end module Rctfld_Module
