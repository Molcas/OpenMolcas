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

use Constants, only: Zero, One, Two, Three, Four, Ten, Half, Pi, auTokJ, kBoltzmann
use Definitions, only: wp, iwp, u6, r8

implicit none
integer(kind=iwp), intent(in) :: LuSpool
integer(kind=iwp) :: i, I_Sph, i_sph_inp, iChrct, ii, iPrint, iRout, ITypRad, jRout,Last, n
real(kind=wp) :: aArea, epscm, poltot, r_min_Sphere, Radius, tal, Temp, val
character(len=180) :: KWord, Key
integer(kind=iwp), external :: iCLast, nToken, NumSolv
real(kind=r8), external :: Anal_Gitt
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
Cordsi(1,latato) = Half
Cordsi(2,latato) = Half
Cordsi(3,latato) = Half

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
! a distance distSparce from ALL atoms and XF points
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
!                                                                      *
!***********************************************************************
!                                                                      *
iPrint = 5
!                                                                      *
!***********************************************************************
!                                                                      *
! Process the input

998 continue
read(LuSpool,'(A72)',end=977,ERR=988) Key
KWord = Key
call UpCase(KWord)
if (KWord(1:1) == '*') Go To 998
if (KWord == '') Go To 998
if (KWord(1:4) == 'REAC') Go To 900
if (KWord(1:4) == 'PRIN') Go To 910
if (KWord(1:4) == 'LANG') Go To 911
if (KWord(1:4) == 'GITT') Go To 912
if (KWord(1:4) == 'POLA') Go To 913
if (KWord(1:4) == 'DIPO') Go To 914
if (KWord(1:4) == 'OSCA') Go To 915
if (KWord(1:4) == 'DIED') Go To 916
if (KWord(1:4) == 'MXLX') Go To 917
if (KWord(1:4) == 'AFAC') Go To 918
if (KWord(1:4) == 'RFBA') Go To 919
if (KWord(1:4) == 'PCM-') Go To 920
if (KWord(1:4) == 'SOLV') Go To 921
if (KWord(1:4) == 'DIEL') Go To 922
if (KWord(1:4) == 'COND') Go To 923
if (KWord(1:4) == 'AARE') Go To 925
if (KWord(1:4) == 'R-MI') Go To 926
if (KWord(1:4) == 'PAUL') Go To 927
if (KWord(1:4) == 'SPHE') Go To 928
if (KWord(1:4) == 'TEMP') Go To 929
if (KWord(1:4) == 'RSCA') Go To 930
if (KWord(1:4) == 'ROTA') Go To 931
if (KWord(1:4) == 'SPAR') Go To 932
if (KWord(1:4) == 'IDAM') Go To 933
if (KWord(1:4) == 'AVER') Go To 934
if (KWord(1:4) == 'CONV') Go To 935
if (KWord(1:4) == 'SKIP') Go To 936
if (KWord(1:4) == 'NODA') Go To 937
if (KWord(1:4) == 'LRAD') Go To 938
if (KWord(1:4) == 'AMBE') Go To 939
if (KWord(1:4) == 'SC14') Go To 940
if (KWord(1:4) == 'DRES') Go To 941
if (KWord(1:4) == 'END ') Go To 997
iChrct = len(KWord)
Last = iCLast(KWord,iChrct)
write(u6,'(1X,A,A)') KWord(1:Last),' is not a keyword!'
call WarningMessage(2,'InpRct: Error in keyword.')
call Quit_OnUserError()
!                                                                      *
!***** REAC ************************************************************
!                                                                      *
! Read reaction field parameters

900 continue
KWord = Get_Ln(LuSpool)
call Get_F1(1,Eps)
call Get_F1(2,rds)
call Get_I1(3,lmax)
if (nToken(KWord) > 3) call Get_F1(4,EpsInf)

lRF = .true.
lRFCav = .true.
write(KWord,'(A,F10.5,A,F10.5,A,I4)') 'eps=',Eps,' radius=',rds,' higest moment=',lMax
Go To 998
!                                                                      *
!***** PRIN ************************************************************
!                                                                      *
! Print level

910 continue
KWord = Get_Ln(LuSpool)
call Get_I1(1,n)
do i=1,n
  KWord = Get_Ln(LuSpool)
  call Get_I1(1,jRout)
  call Get_I1(2,iPrint)
  nPrint(jRout) = iPrint
end do
Go To 998
!                                                                      *
!***** LANG ************************************************************
!                                                                      *
911 continue
KWord = Get_Ln(LuSpool)
call Get_F1(1,scala)
call Get_F1(2,scalb)
call Get_F1(3,scalc)
lLangevin = .true.
lRF = .true.
Go To 998
!                                                                      *
!***** GITT ************************************************************
!                                                                      *
912 continue
KWord = Get_Ln(LuSpool)
call Get_I1(1,latato)
do i=1,latato
  KWord = Get_Ln(LuSpool)
  call Get_F(1,Cordsi(1,i),3)
end do
Go To 998
!                                                                      *
!***** POLA ************************************************************
!                                                                      *
913 continue
KWord = Get_Ln(LuSpool)
call Get_F1(1,polsi)
Go To 998
!                                                                      *
!***** DIPO ************************************************************
!                                                                      *
914 continue
KWord = Get_Ln(LuSpool)
call Get_F1(1,dipsi)
Go To 998
!                                                                      *
!***** OSCA ************************************************************
!                                                                      *
915 continue
KWord = Get_Ln(LuSpool)
call Get_F1(1,scaaa)
Go To 998
!                                                                      *
!***** DIED ************************************************************
!                                                                      *
916 continue
KWord = Get_Ln(LuSpool)
call Get_F1(1,diedel)
Go To 998
!                                                                      *
!***** MXLX ************************************************************
!                                                                      *
917 continue
KWord = Get_Ln(LuSpool)
call Get_F1(1,radlat)
Go To 998
!                                                                      *
!***** AFAC ************************************************************
!                                                                      *
918 continue
KWord = Get_Ln(LuSpool)
call Get_F1(1,afac)
if (afac >= One) then
  call WarningMessage(2,'InpRct: afac invalid value!;        afac >= 1.0 !')
  call Quit_OnUserError()
end if
Go To 998
!                                                                      *
!***** RFBA ************************************************************
!                                                                      *
919 continue
!RF_Basis = .true.
Go To 998
!                                                                      *
!***** PCM- ************************************************************
!                                                                      *
920 continue
PCM = .true.
lRF = .true.
lLangevin = .false.
Go To 998
!                                                                      *
!***** SOLV ************************************************************
!                                                                      *
921 continue
Solvent = trim(Get_Ln(LuSpool))
ISlPar(15) = NumSolv(Solvent)
Go To 998
!                                                                      *
!***** DIEL ************************************************************
!                                                                      *
922 continue
KWord = Get_Ln(LuSpool)
call Get_F1(1,Eps_User)
if (nToken(KWord) > 1) call Get_F1(2,EpsInf_User)
Go To 998
!                                                                      *
!***** COND ************************************************************
!                                                                      *
923 continue
Conductor = .true.
ISlPar(16) = 2
Go To 998
!                                                                      *
!***** AARE ************************************************************
!                                                                      *
925 continue
KWord = Get_Ln(LuSpool)
call Get_F1(1,aArea)
RSlPar(7) = aArea
Go To 998
!                                                                      *
!***** RMIN ************************************************************
!                                                                      *
926 continue
KWord = Get_Ln(LuSpool)
call Get_F1(1,r_min_Sphere)
RSlPar(3) = r_min_Sphere
Go To 998
!                                                                      *
!***** PAUL ************************************************************
!                                                                      *
927 continue
ITypRad = 2
ISlPar(9) = ITypRad
Go To 998
!                                                                      *
!***** SPHE ************************************************************
!                                                                      *
928 continue
KWord = Get_Ln(LuSpool)
call Get_I1(1,I_Sph)
call Get_F1(2,Radius)
ITypRad = 3
ISlPar(9) = ITypRad
i_sph_inp = i_sph_inp+1
NOrdInp(i_sph_inp) = I_Sph
RadInp(i_sph_inp) = Radius
ISlPar(14) = i_sph_inp
Go To 998
!                                                                      *
!***** TEMP ************************************************************
!                                                                      *
! Temperature for Langevin model
929 continue
KWord = Get_Ln(LuSpool)
call Get_F1(1,Temp)
tK = kBoltzmann/auTokJ*1.0e-3_wp*Temp
Go To 998
!                                                                      *
!***** RSCA ************************************************************
!                                                                      *
! Simultaneous scaling of all radii
930 continue
KWord = Get_Ln(LuSpool)
call Get_F1(1,rsca)
Go To 998
!                                                                      *
!***** ROTA ************************************************************
!                                                                      *
! Rotation of the grid with Euler angles alpha, beta, gamma (degrees)
931 continue
KWord = Get_Ln(LuSpool)
call Get_F1(1,rotAlpha)
call Get_F1(2,rotBeta)
call Get_F1(3,rotGamma)

rotAlpha = rotAlpha/180.0_wp*Pi
rotBeta = rotBeta/180.0_wp*Pi
rotGamma = rotGamma/180.0_wp*Pi
Go To 998
!                                                                      *
!***** SPAR ************************************************************
!                                                                      *
! Use sparse grid outside a distance distSparse from all atoms
932 continue
KWord = Get_Ln(LuSpool)
LSparse = .true.
call Get_I1(1,nSparse)
call Get_F1(2,distSparse)
Go To 998
!                                                                      *
!***** IDAM ************************************************************
!                                                                      *
! change dampIter
933 continue
KWord = Get_Ln(LuSpool)
call Get_F1(1,dampIter)
Go To 998
!                                                                      *
!***** AVER ************************************************************
!                                                                      *
! requests average among a number of random grids in the Langevin
! routine. Note: Average is taken for each SCF iteration and only
! the last grid is used for updating the density
! nGridSeed is RNG seed, 0=no seed, -1=system timer
934 continue
KWord = Get_Ln(LuSpool)
LGridAverage = .true.
call Get_I1(1,nGridAverage)
call Get_I1(2,nGridSeed)
Go To 998
!                                                                      *
!***** CONV ************************************************************
!                                                                      *
! Change convergence threshold
935 continue
KWord = Get_Ln(LuSpool)
call Get_F1(1,clim)
Go To 998
!                                                                      *
!***** SKIP ************************************************************
!                                                                      *
! Set minimal distance for dipole-dipole interactions
936 continue
KWord = Get_Ln(LuSpool)
call Get_F1(1,dipCutoff)
Go To 998
!                                                                      *
!***** NODA ************************************************************
!                                                                      *
! Turn off damping of dipole-dipole interactions
937 continue
lDamping = .false.
Go To 998
!                                                                      *
!***** LRAD ************************************************************
!                                                                      *
! Input non-standard Langevin radii
938 continue
Kword = Get_Ln(LuSpool)
call UpCase(KWord)
if (KWord(1:1) == '*') Go To 938
if (KWord == '') Go To 938
if (KWord(1:3) == 'END') Go To 998
read(KWord,*,err=988) ii,val
CovRadT_(ii) = val
!write(u6,*) 'COVR',ii,val
Go To 938
!                                                                      *
!***** AMBE ************************************************************
!                                                                      *
! Use dipole-dipole cutoff of AMBER type
939 continue
lAmberPol = .true.
Go To 998
!                                                                      *
!***** SC14 ************************************************************
!                                                                      *
! Scaling of negative entries in exclusion list
940 continue
KWord = Get_Ln(LuSpool)
call Get_F1(1,scal14)
Go To 998
!                                                                      *
!***** DRES ************************************************************
!                                                                      *
! Restart dipoles from scratch in each QM iteration
! This sometimes gives better convergence.
941 continue
lDiprestart = .true.
Go To 998
!                                                                      *
!***** END  ************************************************************
!                                                                      *
! End of input

997 continue
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
!                                                                      *
!***********************************************************************
!                                                                      *
! Error handling

977 continue
call WarningMessage(2,'InpRct: Premature end of input file.')
call Quit_OnUserError()
988 continue
call WarningMessage(2,'InpRct: Error while reading input file.')
call Quit_OnUserError()

end subroutine InpRct
