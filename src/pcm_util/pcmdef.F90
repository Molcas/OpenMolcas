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

subroutine PCMDef(ISlPar,RSlPar,iPrint)
! Set PCM defaults.

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: ISlPar(100)
real(kind=wp), intent(out) :: RSlPar(100)
integer(kind=iwp), intent(in) :: iPrint
integer(kind=iwp) :: I

! Initialize the integer array.

do I=1,100
  ISlPar(I) = 0
end do

ISlPar(1) = 0       ! SCRF flag
ISlPar(2) = 0       ! derivative level
ISlPar(3) = 0       ! 0-3=invqsep/invqtot/iterqsep/iterqtot
ISlPar(4) = 0       ! geomview
ISlPar(5) = 2       ! escaped charge compensation type
ISlPar(6) = 1       ! cavitation energy
ISlPar(7) = 1       ! repulsion energy
ISlPar(8) = 1       ! dispersion energy
ISlPar(9) = 1       ! type of atomic radii (UATM/Pauling/input)
ISlPar(10) = 2      ! type of polyhedra
ISlPar(11) = -400   ! initial number of tesserae
ISlPar(12) = 1      ! level of division in PolyGen
ISlPar(13) = 0      ! derivatives style (default/pcm/cond)
ISlPar(14) = 0      ! number of explicit spheres
ISlPar(15) = 1      ! solvent id number
ISlPar(16) = 1      ! PCM / Conductor
ISlPar(17) = 1      ! [unused]
ISlPar(18) = 1      ! [unused]
ISlPar(19) = 0      ! [unused]
ISlPar(20) = 0      ! non equilibrium solvation
ISlPar(21) = 0      ! fixed cavity gradients
ISlPar(22) = 0      ! fixed cavity 2nd derivative
ISlPar(23) = 0      ! number of spheres
ISlPar(24) = 3      ! atom types in solvent
ISlPar(25) = 0      ! [unused]
ISlPar(26) = 1      ! [unused]
ISlPar(27) = 0      ! explicit solute charges
ISlPar(28) = 0      ! gives PCM help
ISlPar(29) = 0      ! fits charges to electrostatic potential
ISlPar(30) = 0      ! [unused]
ISlPar(31) = 0      ! [unused]
ISlPar(32) = 0      ! [unused]
ISlPar(33) = IPrint ! specific printing level for solvent routines
ISlPar(34) = 0      ! number of spheres
ISlPar(35) = 0      ! number of tesserae
ISlPar(36) = 0      ! [unused]
ISlPar(37) = 0      ! read in radius and alpha for atom/element
ISlPar(38) = 0      ! read in element parameters for UA method
ISlPar(39) = 0      ! local field factor switch
ISlPar(40) = 0      ! [unused]
ISlPar(41) = 0      ! vacuum SCF flag
ISlPar(42) = 0      ! number of solute atoms

! Initialize the real array.

do I=1,100
  RSlPar(I) = Zero
end do

RSlPar(1) = Zero        ! XCosmo
RSlPar(2) = 40.0_wp     ! Omega
RSlPar(3) = 0.2_wp      ! Ret, minimum radius of added sphere
RSlPar(4) = 0.7_wp      ! Fro
RSlPar(5) = 1.0e-4_wp   ! DR
RSlPar(6) = 1.0e-6_wp   ! accuracy for PCM
RSlPar(7) = 0.4_wp      ! average tesserae area
RSlPar(8) = 78.39_wp    ! dielectric constant
RSlPar(9) = 1.2_wp      ! Alpha scaling of radii
RSlPar(10) = 78.39_wp   ! Exx
RSlPar(11) = Zero       ! Exy
RSlPar(12) = Zero       ! Exz
RSlPar(13) = Zero       ! Eyy
RSlPar(14) = Zero       ! Eyz
RSlPar(15) = Zero       ! Ezz
RSlPar(16) = 298.15_wp  ! absolute temperature
RSlPar(17) = 1.776_wp   ! limiting dielectric constant
RSlPar(18) = -0.3562_wp ! derivative of eps wrt T
RSlPar(19) = 1.385_wp   ! RSolv
RSlPar(20) = 2.57e-4_wp ! TCE
RSlPar(21) = 78.39_wp   ! dielectric constant for shell  1
RSlPar(22) = One        ! dielectric constant for shell  2
RSlPar(23) = One        ! dielectric constant for shell  3
RSlPar(24) = One        ! dielectric constant for shell  4
RSlPar(25) = One        ! dielectric constant for shell  5
RSlPar(26) = One        ! dielectric constant for shell  6
RSlPar(27) = One        ! dielectric constant for shell  7
RSlPar(28) = One        ! dielectric constant for shell  8
RSlPar(29) = One        ! dielectric constant for shell  9
RSlPar(30) = One        ! dielectric constant for shell 10
RSlPar(31) = 1.2_wp     ! scaling factor for neutrals
RSlPar(32) = One        ! scaling factor for acidic hydrogens
RSlPar(33) = 1.1_wp     ! scaling factor for charged
RSlPar(34) = 18.07_wp   ! solvent molecular volume
RSlPar(35) = 71.81_wp   ! Sten
RSlPar(36) = 0.650_wp   ! DSten
RSlPar(37) = 1.277_wp   ! CMF
RSlPar(38) = 3.348e-2_wp! solvent density
RSlPar(39) = One        ! DK2
RSlPar(40) = One        ! DAlp
RSlPar(41) = One        ! DetEps
RSlPar(42) = 30.0_wp    ! hard cutoff for iterative PCM
RSlPar(43) = 15.0_wp    ! multiplier cutoff for iterative PCM
RSlPar(44) = 0.8_wp     ! overlap factor for GePol
RSlPar(45) = 0.5_wp     ! min radius for added spheres by GePol
RSlPar(46) = Zero       ! GCavP
RSlPar(47) = Zero       ! SCavP
RSlPar(48) = Zero       ! HCavP
RSlPar(49) = Zero       ! GDisp
RSlPar(50) = Zero       ! GRep
RSlPar(51) = Zero       ! Total Cavity surface (SCav)
RSlPar(52) = Zero       ! Total Cavity volume (VCav)

return

end subroutine PCMDef
