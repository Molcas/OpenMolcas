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
!     diedel: controlls deletion of gitter polarizabilities            *
!     tK: Boltzman factor                                              *
!     clim: convergence threshold for solver                           *
!     afac: controlls equation solver (valid values 0.1-0.97)          *
!     nexpo: exponent of the potential with which the QM system kills  *
!            the gitter polatizabilities                               *
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
!     NOrdInp: number of atom where explicit spheres are centered      *
!     RadInp: radius of spheres explicitly given in the input          *
!                                                                      *
!***********************************************************************

module Rctfld_Module

use Definitions, only: wp, iwp

implicit none
private

integer, parameter :: MxA = 1000, MxPar = 100
integer(kind=iwp) :: iCharge_ref, ISlPar(MxPar), latato, lMax, maxa, maxb, maxc, nCavxyz, nexpo, nGrid, nGrid_Eff, nGridAverage, &
                     nGridSeed, NOrdInp(MxA), nPCM_Info, NS, NSinit, nSparse, nTs
real(kind=wp) :: afac, clim, Cordsi(3,4), dampIter, diedel, dipCutoff, dipsi, distSparse, Eps, Eps_User, EpsInf, EpsInf_User, &
                 fmax, gatom, polsi, prefac, RadInp(MxA), radlat, rds, rotAlpha, rotBeta, rotGamma, rsca, RSlPar(MxPar), RSolv, &
                 scaaa, scal14, scala, scalb, scalc, tK, VMol
logical(kind=iwp) :: Conductor, DoDeriv, Done_Lattice, lAmberPol, lDamping, lDiprestart, lFirstIter, LGridAverage, lLangevin, lRF, &
                     lRFCav, LSparse, NonEq_ref, PCM
character(len=32) :: Solvent
real(kind=wp), allocatable :: MM(:,:)

integer(kind=iwp) :: cRFEnd, cRFStrt, iRFEnd, iRFStrt, lRFEnd, lRFStrt
real(kind=wp) :: rRFEnd, rRFStrt

common/iRct/iRFStrt,lMax,latato,nexpo,maxa,maxb,maxc,nCavxyz,nGrid,nGrid_Eff,nSparse,ISlPar,NSinit,NS,nTs,NOrdInp,nPCM_Info, &
            iCharge_ref,nGridAverage,nGridSeed,iRFEnd
common/rRct/rRFStrt,EpsInf_User,Eps_User,rds,Cordsi,polsi,dipsi,radlat,scala,scalb,scalc,scaaa,gatom,diedel,tK,rotAlpha,rotBeta, &
            rotGamma,distSparse,clim,afac,prefac,fmax,rsca,Eps,EpsInf,RSolv,VMol,RadInp,RSlPar,dampIter,dipCutoff,scal14,rRFEnd
common/lRct/lRFStrt,lRF,lLangevin,PCM,Conductor,NonEq_ref,DoDeriv,lRFCav,LSparse,LGridAverage,lDamping,lAmberPol,Done_Lattice, &
            lFirstIter,lDiprestart,lRFEnd
common/cRct/cRFStrt,Solvent,cRFEnd

public :: aFac, cLim, Conductor, Cordsi, CRFEnd, CRFStrt, DampIter, DieDel, DipCutOff, Dipsi, DistSparse, DoDeriv, Done_Lattice, &
          Eps, Eps_USER, EpsInf, EpsInf_USER, fMax, gAtom, iCharge_Ref, iRFEnd, iRFStrt, ISlPar, lAmberPol, LatAto, lDamping, &
          lDipRestart, lFirstIter, lGridAverage, lLangevin, lMax, lRF, lRFCav, lRFEnd, lRFStrt, lSparse, MaxA, MaxB, MaxC, MM, &
          MXA, nCavxyz, nExpO, nGrid, nGrid_Eff, nGridAverage, nGridSeed, NonEQ_Ref, nOrdInp, nPCM_Info, nS, nSInit, nSparse, nTS, &
          PCM, Polsi, PreFac, RadInp, RadLat, RDS, RotAlpha, RotBeta, RotGamma, rRFEnd, rRFStrt, rSca, RSlPar, RSolv, ScaAA, &
          Scal14, Scala, ScalB, Scalc, Solvent, TK, vMol

end module Rctfld_Module
