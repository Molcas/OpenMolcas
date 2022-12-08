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
! Copyright (C) 1992,2000, Roland Lindh                                *
!***********************************************************************

subroutine InpRct(LuSpool)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!                                                                      *
!             Modified for Langevin polarizabilities, March 2000 (RL)  *
!***********************************************************************

use Constants, only: Zero, One, Two, Three, Four, Ten, Half, Pi, deg2rad, auTokJ, kBoltzmann
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LuSpool
integer(kind=iwp) :: i, I_Sph, i_sph_inp, iChrct, ii, iPrint, iRout, istatus, ITypRad, jRout, Last, n
real(kind=wp) :: aArea, epscm, poltot, r_min_Sphere, Radius, tal, Temp, val
character(len=180) :: KWord, Key
integer(kind=iwp), external :: iCLast, nToken, NumSolv
real(kind=wp), external :: Anal_Gitt
character(len=180), external :: Get_Ln
#include "print.fh"
#include "rctfld.fh"
#include "covradt_data.fh"

iRout = 1
iPrint = nPrint(iRout)
!                                                                      *
!***********************************************************************
!                                                                      *
! Initialize data for dielectric medium

Eps = One
EpsInf = One
Eps_User = -One
EpsInf_User = Zero
rds = Zero
lMax = -1
lRF = .false.
lRFCav = .false.
RF_Basis = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
! Initialize data for PCM

PCM = .false.
i_sph_inp = 0
! Default PCM parameters
call PCMDef(ISlPar,RSlPar,iPrint)
!                                                                      *
!***********************************************************************
!                                                                      *
! Setting up default values for parameters used in Langevin
! polarizabilities

! latato: Gitter type (1 or 4)
latato = 1
Cordsi(:,latato) = Half

! polsi: site polarizability
polsi = Zero

! dipsi: site dipole moment
dipsi = Zero

! radlat: maximum extension of the lattic
radlat = Zero

! scala,scalb,scalc: lengths of the cell dimensions
scala = Zero ! should be set by input
scalb = Zero
scalc = Zero

! scaaa: overall scaling
scaaa = One

! rotation of grid
rotAlpha = Zero
rotBeta = Zero
rotGamma = Zero

! rsca: scaling of radii
rsca = One

! Controls wheather a sparse grid should be used outside
! a distance distSparse from ALL atoms and XF points
! nSparse = scala(sparse) / scala(normal)
LSparse = .false.
nSparse = 1
distSparse = Zero

! Minimum distance for handling dipole-dipole interactions
dipCutoff = Zero

! Damping of dipole-dipole interactions
lDamping = .true.

! Use a dipole-dipole cutoff of Amber type, i.e. same exclusions
! as for the static field
lAmberPol = .false.

! Controls scaling of contributions from negative entries in
! the exclusion list, typically 1-4 interactions
scal14 = One

! Controls wheather an average of different grids is done
LGridAverage = .false.

! diedel: controls deletion of gitter polarizabilities
diedel = 0.01_wp

! tK: Boltzman factor
tK = 0.001_wp

! clim: convergence threshold for solver
clim = 1.0e-15_wp

! afac: controls equation solver (valid values 0.3-0.97)
afac = Half

! iterDamp: controls update of dipoles in equation solver
!newDip = oldDip*iterDamp + estimatedDip*(1-iterDamp)
dampIter = 0.4_wp
lDiprestart = .false.

! nexpo: exponent of the potential with which the QM system kills
!        the gitter polatizabilities
nexpo = 12

! prefac: scaling of the polarization energies
prefac = One

lLangevin = .false.

! default solvent
Solvent = 'WATER'
ISlPar(15) = NumSolv(Solvent)
!                                                                      *
!***********************************************************************
!                                                                      *
iPrint = 5
!                                                                      *
!***********************************************************************
!                                                                      *
! Process the input

do
  read(LuSpool,'(A72)',iostat=istatus) Key
  if (istatus < 0) call Error(1)
  if (istatus > 0) call Error(2)
  KWord = Key
  call UpCase(KWord)
  if ((KWord(1:1) == '*') .or. (KWord == '')) cycle
  select case (KWord(1:4))
    case ('REAC')
      !                                                                *
      !***** REAC ******************************************************
      !                                                                *
      ! Read reaction field parameters
      KWord = Get_Ln(LuSpool)
      call Get_F1(1,Eps)
      call Get_F1(2,rds)
      call Get_I1(3,lmax)
      if (nToken(KWord) > 3) call Get_F1(4,EpsInf)
      lRF = .true.
      lRFCav = .true.
      write(KWord,'(A,F10.5,A,F10.5,A,I4)') 'eps=',Eps,' radius=',rds,' higest moment=',lMax
    case ('PRIN')
      !                                                                *
      !***** PRIN ******************************************************
      !                                                                *
      ! Print level
      KWord = Get_Ln(LuSpool)
      call Get_I1(1,n)
      do i=1,n
        KWord = Get_Ln(LuSpool)
        call Get_I1(1,jRout)
        call Get_I1(2,iPrint)
        nPrint(jRout) = iPrint
      end do
    case ('LANG')
      !                                                                *
      !***** LANG ******************************************************
      !                                                                *
      KWord = Get_Ln(LuSpool)
      call Get_F1(1,scala)
      call Get_F1(2,scalb)
      call Get_F1(3,scalc)
      lLangevin = .true.
      lRF = .true.
    case ('GITT')
      !                                                                *
      !***** GITT ******************************************************
      !                                                                *
      KWord = Get_Ln(LuSpool)
      call Get_I1(1,latato)
      do i=1,latato
        KWord = Get_Ln(LuSpool)
        call Get_F(1,Cordsi(1,i),3)
      end do
    case ('POLA')
      !                                                                *
      !***** POLA ******************************************************
      !                                                                *
      KWord = Get_Ln(LuSpool)
      call Get_F1(1,polsi)
    case ('DIPO')
      !                                                                *
      !***** DIPO ******************************************************
      !                                                                *
      KWord = Get_Ln(LuSpool)
      call Get_F1(1,dipsi)
    case ('OSCA')
      !                                                                *
      !***** OSCA ******************************************************
      !                                                                *
      KWord = Get_Ln(LuSpool)
      call Get_F1(1,scaaa)
    case ('DIED')
      !                                                                *
      !***** DIED ******************************************************
      !                                                                *
      KWord = Get_Ln(LuSpool)
      call Get_F1(1,diedel)
    case ('MXLX')
      !                                                                *
      !***** MXLX ******************************************************
      !                                                                *
      KWord = Get_Ln(LuSpool)
      call Get_F1(1,radlat)
    case ('AFAC')
      !                                                                *
      !***** AFAC ******************************************************
      !                                                                *
      KWord = Get_Ln(LuSpool)
      call Get_F1(1,afac)
      if (afac >= One) then
        call WarningMessage(2,'InpRct: afac invalid value!;        afac >= 1.0 !')
        call Quit_OnUserError()
      end if
    case ('RFBA')
      !                                                                *
      !***** RFBA ******************************************************
      !                                                                *
      !RF_Basis = .true.
    case ('PCM-')
      !                                                                *
      !***** PCM- ******************************************************
      !                                                                *
      PCM = .true.
      lRF = .true.
      lLangevin = .false.
    case ('SOLV')
      !                                                                *
      !***** SOLV ******************************************************
      !                                                                *
      Solvent = trim(Get_Ln(LuSpool))
      ISlPar(15) = NumSolv(Solvent)
    case ('DIEL')
      !                                                                *
      !***** DIEL ******************************************************
      !                                                                *
      KWord = Get_Ln(LuSpool)
      call Get_F1(1,Eps_User)
      if (nToken(KWord) > 1) call Get_F1(2,EpsInf_User)
    case ('COND')
      !                                                                *
      !***** COND ******************************************************
      !                                                                *
      Conductor = .true.
      ISlPar(16) = 2
    case ('AARE')
      !                                                                *
      !***** AARE ******************************************************
      !                                                                *
      KWord = Get_Ln(LuSpool)
      call Get_F1(1,aArea)
      RSlPar(7) = aArea
    case ('R-MI')
      !                                                                *
      !***** RMIN ******************************************************
      !                                                                *
      KWord = Get_Ln(LuSpool)
      call Get_F1(1,r_min_Sphere)
      RSlPar(3) = r_min_Sphere
    case ('PAUL')
      !                                                                *
      !***** PAUL ******************************************************
      !                                                                *
      ITypRad = 2
      ISlPar(9) = ITypRad
    case ('SPHE')
      !                                                                *
      !***** SPHE ******************************************************
      !                                                                *
      KWord = Get_Ln(LuSpool)
      call Get_I1(1,I_Sph)
      call Get_F1(2,Radius)
      ITypRad = 3
      ISlPar(9) = ITypRad
      i_sph_inp = i_sph_inp+1
      if (i_sph_inp > MxA) then
        call WarningMessage(2,'InpRct: i_sph_inp > MxA')
        call Abend()
      end if
      NOrdInp(i_sph_inp) = I_Sph
      RadInp(i_sph_inp) = Radius
      ISlPar(14) = i_sph_inp
    case ('TEMP')
      !                                                                *
      !***** TEMP ******************************************************
      !                                                                *
      ! Temperature for Langevin model
      KWord = Get_Ln(LuSpool)
      call Get_F1(1,Temp)
      tK = kBoltzmann/auTokJ*1.0e-3_wp*Temp
    case ('RSCA')
      !                                                                *
      !***** RSCA ******************************************************
      !                                                                *
      ! Simultaneous scaling of all radii
      KWord = Get_Ln(LuSpool)
      call Get_F1(1,rsca)
    case ('ROTA')
      !                                                                *
      !***** ROTA ******************************************************
      !                                                                *
      ! Rotation of the grid with Euler angles alpha, beta, gamma (degrees)
      KWord = Get_Ln(LuSpool)
      call Get_F1(1,rotAlpha)
      call Get_F1(2,rotBeta)
      call Get_F1(3,rotGamma)
      rotAlpha = rotAlpha*deg2rad
      rotBeta = rotBeta*deg2rad
      rotGamma = rotGamma*deg2rad
    case ('SPAR')
      !                                                                *
      !***** SPAR ******************************************************
      !                                                                *
      ! Use sparse grid outside a distance distSparse from all atoms
      KWord = Get_Ln(LuSpool)
      LSparse = .true.
      call Get_I1(1,nSparse)
      call Get_F1(2,distSparse)
    case ('IDAM')
      !                                                                *
      !***** IDAM ******************************************************
      !                                                                *
      ! change dampIter
      KWord = Get_Ln(LuSpool)
      call Get_F1(1,dampIter)
    case ('AVER')
      !                                                                *
      !***** AVER ******************************************************
      !                                                                *
      ! requests average among a number of random grids in the Langevin
      ! routine. Note: Average is taken for each SCF iteration and only
      ! the last grid is used for updating the density
      ! nGridSeed is RNG seed, 0=no seed, -1=system timer
      KWord = Get_Ln(LuSpool)
      LGridAverage = .true.
      call Get_I1(1,nGridAverage)
      call Get_I1(2,nGridSeed)
    case ('CONV')
      !                                                                *
      !***** CONV ******************************************************
      !                                                                *
      ! Change convergence threshold
      KWord = Get_Ln(LuSpool)
      call Get_F1(1,clim)
    case ('SKIP')
      !                                                                *
      !***** SKIP ******************************************************
      !                                                                *
      ! Set minimal distance for dipole-dipole interactions
      KWord = Get_Ln(LuSpool)
      call Get_F1(1,dipCutoff)
    case ('NODA')
      !                                                                *
      !***** NODA ******************************************************
      !                                                                *
      ! Turn off damping of dipole-dipole interactions
      lDamping = .false.
    case ('LRAD')
      !                                                                *
      !***** LRAD ******************************************************
      !                                                                *
      ! Input non-standard Langevin radii
      do
        Kword = Get_Ln(LuSpool)
        call UpCase(KWord)
        if (KWord(1:1) == '*') cycle
        if (KWord == '') cycle
        if (KWord(1:3) == 'END') exit
        read(KWord,*,iostat=istatus) ii,val
        if (istatus > 0) call Error(2)
        CovRadT_(ii) = val
        !write(u6,*) 'COVR',ii,val
      end do
    case ('AMBE')
      !                                                                *
      !***** AMBE ******************************************************
      !                                                                *
      ! Use dipole-dipole cutoff of AMBER type
      lAmberPol = .true.
    case ('SC14')
      !                                                                *
      !***** SC14 ******************************************************
      !                                                                *
      ! Scaling of negative entries in exclusion list
      KWord = Get_Ln(LuSpool)
      call Get_F1(1,scal14)
    case ('DRES')
      !                                                                *
      !***** DRES ******************************************************
      !                                                                *
      ! Restart dipoles from scratch in each QM iteration
      ! This sometimes gives better convergence.
      lDiprestart = .true.
    case ('END ')
      !                                                                *
      !***** END  ******************************************************
      !                                                                *
      ! End of input
      exit
    case default
      iChrct = len(KWord)
      Last = iCLast(KWord,iChrct)
      write(u6,'(1X,A,A)') KWord(1:Last),' is not a keyword!'
      call WarningMessage(2,'InpRct: Error in keyword.')
      call Quit_OnUserError()
  end select
end do
!                                                                      *
!***********************************************************************
!                                                                      *
if (lLangevin) then

  polsi = polsi*scaaa**Three
  scala = scala*scaaa
  scalb = scalb*scaaa
  scalc = scalc*scaaa
  dipsi = dipsi*scaaa**(Three/Two)

  ! Analyze the lattice

  gatom = Anal_Gitt(cordsi,latato)

  ! Scale the grid.

  do i=1,latato
    cordsi(1,i) = cordsi(1,i)*scala
    cordsi(2,i) = cordsi(2,i)*scalb
    cordsi(3,i) = cordsi(3,i)*scalc
  end do

  if (lRFCav) then
    if (radlat == Zero) radlat = rds-One/Ten
    if (radlat > rds) then
      call WarningMessage(2,'InpRct: radlat > rds')
      call Abend()
    end if
  else
    call WarningMessage(1,'Running Langevin without reaction field')
  end if

  tK = One/tK
  poltot = gatom*(polsi+dipsi**2*tK/Three)
  tal = poltot*Four*Pi/(Three*scala*scalb*scalc)
  epscm = (One+Two*tal)/(One-tal)
  if (Eps < One) Eps = epscm

  tk5 = Half*tk

end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

contains

subroutine Error(code)
  integer(kind=iwp), intent(in) :: code
  select case (code)
    case (1)
      call WarningMessage(2,'InpRct: Premature end of input file.')
    case (2)
      call WarningMessage(2,'InpRct: Error while reading input file.')
  end select
  call Quit_OnUserError()
end subroutine Error

end subroutine InpRct
