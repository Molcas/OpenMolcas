************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine PCMDef(ISlPar,RSlPar,iPrint)
      Implicit Real*8 (A-H,O-Z)
*
*     Set PCM defaults.
*
      Integer ISlPar(100)
      Real*8 RSlPar(100)
*
*     Initialize the integer array.
*
      Do 10 I = 1, 100
        ISlPar(I) = 0
  10  Continue
*
      ISlPar(  1) =  0      ! SCRF flag
      ISlPar(  2) =  0      ! derivative level
      ISlPar(  3) =  0      ! 0-3=invqsep/invqtot/iterqsep/iterqtot
      ISlPar(  4) =  0      ! geomview
      ISlPar(  5) =  2      ! escaped charge compensation type
      ISlPar(  6) =  1      ! cavitation energy
      ISlPar(  7) =  1      ! repulsion energy
      ISlPar(  8) =  1      ! dispersion energy
      ISlPar(  9) =  1      ! type of atomic radii (UATM/Pauling/input)
      ISlPar( 10) =  2      ! type of polyhedra
      ISlPar( 11) = -400    ! initial number of tesserae
      ISlPar( 12) =  1      ! level of division in PolyGen
      ISlPar( 13) =  0      ! derivatives style (default/pcm/cond)
      ISlPar( 14) =  0      ! number of explicit spheres
      ISlPar( 15) =  1      ! solvent id number
      ISlPar( 16) =  1      ! PCM / Conductor
      ISlPar( 17) =  1      ! [unused]
      ISlPar( 18) =  1      ! [unused]
      ISlPar( 19) =  0      ! [unused]
      ISlPar( 20) =  0      ! non equilibrium solvation
      ISlPar( 21) =  0      ! fixed cavity gradients
      ISlPar( 22) =  0      ! fixed cavity 2nd derivative
      ISlPar( 23) =  0      ! number of spheres
      ISlPar( 24) =  3      ! atom types in solvent
      ISlPar( 25) =  0      ! [unused]
      ISlPar( 26) =  1      ! [unused]
      ISlPar( 27) =  0      ! explicit solute charges
      ISlPar( 28) =  0      ! gives PCM help
      ISlPar( 29) =  0      ! fits charges to electrostatic potential
      ISlPar( 30) =  0      ! [unused]
      ISlPar( 31) =  0      ! [unused]
      ISlPar( 32) =  0      ! [unused]
      ISlPar( 33) =  IPrint ! specific printing level for
                            !            solvent routines
      ISlPar( 34) =  0      ! number of spheres
      ISlPar( 35) =  0      ! number of tesserae
      ISlPar( 36) =  0      ! [unused]
      ISlPar( 37) =  0      ! read in radius and alpha for atom/element
      ISlPar( 38) =  0      ! read in element parameters for UA method
      ISlPar( 39) =  0      ! local field factor switch
      ISlPar( 40) =  0      ! [unused]
      ISlPar( 41) =  0      ! vacuum SCF flag
      ISlPar( 42) =  0      ! number of solute atoms

**     Initialize the real array.
*
      Do 20 I = 1, 100
        RSlPar(I) = 0.0d0
  20  Continue
*
      RSlPar(  1) =   0.0d0     ! XCosmo
      RSlPar(  2) =   4.0d+1    ! Omega
      RSlPar(  3) =   2.0d-1    ! Ret, minimum radius of added sphere
      RSlPar(  4) =   7.0d-1    ! Fro
      RSlPar(  5) =   1.0d-4    ! DR
      RSlPar(  6) =   1.0d-6    ! accuracy for PCM
      RSlPar(  7) =   4.0d-1    ! average tesserae area
      RSlPar(  8) =  78.39d0    ! dielectric constant
      RSlPar(  9) =   1.2d0     ! Alpha scaling of radii
      RSlPar( 10) =  78.39d0    ! Exx
      RSlPar( 11) =   0.0d0     ! Exy
      RSlPar( 12) =   0.0d0     ! Exz
      RSlPar( 13) =   0.0d0     ! Eyy
      RSlPar( 14) =   0.0d0     ! Eyz
      RSlPar( 15) =   0.0d0     ! Ezz
      RSlPar( 16) = 298.15d0    ! absolute temperature
      RSlPar( 17) =   1.776d0   ! limiting dielectric constant
      RSlPar( 18) =  -0.3562d0  ! derivative of eps wrt T
      RSlPar( 19) =   1.385d0   ! RSolv
      RSlPar( 20) =   2.57d-4   ! TCE
      RSlPar( 21) =  78.39d0    ! dielectric constant for shell  1
      RSlPar( 22) =   1.0d0     ! dielectric constant for shell  2
      RSlPar( 23) =   1.0d0     ! dielectric constant for shell  3
      RSlPar( 24) =   1.0d0     ! dielectric constant for shell  4
      RSlPar( 25) =   1.0d0     ! dielectric constant for shell  5
      RSlPar( 26) =   1.0d0     ! dielectric constant for shell  6
      RSlPar( 27) =   1.0d0     ! dielectric constant for shell  7
      RSlPar( 28) =   1.0d0     ! dielectric constant for shell  8
      RSlPar( 29) =   1.0d0     ! dielectric constant for shell  9
      RSlPar( 30) =   1.0d0     ! dielectric constant for shell 10
      RSlPar( 31) =   1.2d0     ! scaling factor for neutrals
      RSlPar( 32) =   1.0d0     ! scaling factor for acidic hydrogens
      RSlPar( 33) =   1.1d0     ! scaling factor for charged
      RSlPar( 34) =  18.07d0    ! solvent molecular volume
      RSlPar( 35) =  71.81d0    ! Sten
      RSlPar( 36) =   0.650d0   ! DSten
      RSlPar( 37) =   1.277d0   ! CMF
      RSlPar( 38) =   3.348d-2  ! solvent density
      RSlPar( 39) =   0.0d0     ! DK2
      RSlPar( 40) =   0.0d0     ! DAlp
      RSlPar( 41) =   0.0d0     ! DetEps
      RSlPar( 42) =  30.0d0     ! hard cutoff for iterative PCM
      RSlPar( 43) =  15.0d0     ! multiplier cutoff for iterative PCM
      RSlPar( 44) =   0.8d0     ! overlap factor for GePol
      RSlPar( 45) =   0.5d0     ! min radius for added spheres by GePol
      RSlPar( 46) =   0.0d0     ! GCavP
      RSlPar( 47) =   0.0d0     ! SCavP
      RSlPar( 48) =   0.0d0     ! HCavP
      RSlPar( 49) =   0.0d0     ! GDisp
      RSlPar( 50) =   0.0d0     ! GRep
      RSlPar( 51) =   0.0d0     ! Total Cavity surface (SCav)
      RSlPar( 52) =   0.0d0     ! Total Cavity volume (VCav)
*
      Return
      End
