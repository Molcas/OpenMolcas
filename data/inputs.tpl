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

@SCF/DFT single point
*
* SCF/DFT single energy
*
&GATEWAY
&SEWARD
&SCF
&GRID_IT

@SCF/DFT geometry optimization
*
* SCF/DFT geometry optimization
*
&GATEWAY
>> DoWhile
&SEWARD
&SCF
&SLAPAF
>> EndDo
&GRID_IT

@SCF/DFT geo. opt. and freq. calculation
*
* SCF/DFT geometry optimization and frequency calculation
*
&GATEWAY
>> DoWhile
&SEWARD
&SCF
&SLAPAF
>> EndDo
&MCKINLEY
&GRID_IT

@MP2 single point
*
* MP2 single energy
*
&GATEWAY
&SEWARD
&SCF
&MBPT2
&GRID_IT

@MP2 geometry optimization
*
* MP2 geometry optimization
*
&GATEWAY
>> DoWhile
&SEWARD
&SCF
&MBPT2
&SLAPAF
>> EndDo
&GRID_IT

@MP2 geo. opt. and freq. calculation
*
* MP2 geometry optimization and frequency calculation
*
&GATEWAY
>> DoWhile
&SEWARD
&SCF
&MBPT2
&SLAPAF
>> EndDo
&MCKINLEY
&GRID_IT

@CASSCF/RASSCF single point
*
* CASSCF/RASSCF single energy
*
&GATEWAY
&SEWARD
&SCF
&RASSCF
  LumOrb
&GRID_IT

@CASSCF/RASSCF geometry optimization
*
* CASSCF/RASSCF geometry optimization
*
&GATEWAY
>> DoWhile
&SEWARD
&RASSCF
  FileOrb = _Start_Orbitals
&SLAPAF
>> EndDo
&GRID_IT

@CASSCF/RASSCF geo. opt. and freq. calculation
*
* CASSCF/RASSCF geometry optimization and frequency calculation in RASSCF level
*
&GATEWAY
>> DoWhile
&SEWARD
&RASSCF
  FileOrb = _Start_Orbitals
&SLAPAF
>> EndDo
&MCKINLEY

@CASSCF/RASSCF IRC from a transition state
*
* CASSCF/RASSCF intrinsic reaction path from a transition state
*
&GATEWAY
>>DoWhile
&SEWARD
&RASSCF
  FileOrb = _Start_Orbitals
&SLAPAF
  IRC
  nIRC = 5
  ReactionVector = 3*n_real_numbers
>>EndDo

@CASPT2/MS-CASPT2 single point
*
* CASPT2/MS-CASPT2 single energy
*
&GATEWAY
&SEWARD
&SCF
&RASSCF
  LumOrb
&CASPT2
&GRID_IT

@CASPT2/MS-CASPT2 geometry optimization
*
* CASPT2/MS-CASPT2 geometry optimization
*
&GATEWAY
>> DoWhile
&SEWARD
&RASSCF
  FileOrb = _Start_Orbitals
&CASPT2
&SLAPAF
>> EndDo

@CASPT2/MS-CASPT2 geo. opt. and freq. calculation
*
* CASPT2/MS-CASPT2 geometry optimization and frequency calculation
*
&GATEWAY
>> DoWhile
&SEWARD
&RASSCF
  FileOrb = _Start_Orbitals
&CASPT2
&SLAPAF
>> EndDo
&MCKINLEY

@Localized Molecular Orbitals for CASSCF/RASSCF
*
* Generating Starting Orbitals for CASSCF/RASSCF in ANO-RCC-MB
*
&GATEWAY
  Basis = ANO-RCC-MB
  Group = Full
  RICD
&SEWARD
&SCF
&LOCALISATION
  FileOrb = $Project.ScfOrb
  Occupied
&LOCALISATION
  FileOrb = $Project.LocOrb
  Virtual
&GRID_IT
  FileOrb = $Project.LocOrb
  All
>> Unix molcas gv $Project.grid
&RASSCF
  FileOrb = $Project.GvOrb
&GRID_IT
  FileOrb=$Project.RasOrb
  All
>> Unix molcas gv $Project.grid

@Localized Atomic Orbitals for CASSCF/RASSCF
*
* Generating Starting Orbitals for CASSCF/RASSCF
*
&GATEWAY
  Basis = ANO-RCC-MB
  Group = Full
  RICD
&LOCALISATION
  FileOrb = $Project.GssOrb
  All
&GRID_IT
  FileOrb = $Project.LocOrb
  All
>> Unix molcas gv $Project.grid
&RASSCF
  FileOrb = $Project.GvOrb
&GRID_IT
  FileOrb = $Project.RasOrb
  All
>> Unix molcas gv $Project.grid

@SCF Orbitals for CASSCF/RASSCF
*
* Generating Starting Orbitals for CASSCF/RASSCF in ANO-RCC-MB
*
&GATEWAY
  Basis = ANO-RCC-MB
  Group = Full
  RICD
&SEWARD
&SCF
&GRID_IT
  All
  FileOrb = $Project.ScfOrb
>> UnIX molcas gv $Project.grid
&RASSCF
  FileOrb = $Project.GvOrb
&GRID_IT
  FileOrb = $Project.RasOrb
  All
>> Unix molcas gv $Project.grid

@Expand Orbitals from one basis set to another
*
* Generate a set of orbitals based on orbitals from a smaller basis set
*
&GATEWAY
  RICD
>> Copy $Project.RunFile RUNFIL1
&GATEWAY
  RICD
>> Copy $Project.RunFile RUNFIL2
>> Copy External_Orbital_File $Project.RasOrb_External
&EXPBAS
  FileOrb = $Project.RasOrb_External
&SEWARD
&RASSCF
  FileOrb = $Project.ExpOrb

@Localized SCF Orbitals for CASSCF/RASSCF & Expand
*
* Generating Starting Orbitals for CASSCF/RASSCF
*
&GATEWAY
  Basis = ANO-RCC-MB
  Group = NoSymm
  RICD
&SEWARD
&SCF
&LOCALISATION
  FileOrb = $Project.ScfOrb
  Occupied
&LOCALISATION
  FileOrb = $Project.LocOrb
  Virtual
&GRID_IT
  FileOrb = $Project.LocOrb
  Name = localized
  All
>> Unix molcas gv $Project.localized.grid
&RASSCF
  FileOrb = $Project.localized.GvOrb
>> Copy $Project.RunFile RUNFIL1
&GATEWAY
  Basis = ANO-RCC-VDZP
  Group = NoSymm
  RICD
>> Copy $Project.RunFile RUNFIL2
&EXPBAS
  FileOrb = $Project.RasOrb
&SEWARD
&RASSCF
  FileOrb = $Project.ExpOrb
&GRID_IT
  File = $Project.RasOrb
  All

@MEP with surface hopping
*
* MEP with surface hopping
*
&GATEWAY
>> DoWhile
&SEWARD
&RASSCF
  FileOrb = _Start_Orbitals
&RASSI
  HOP
&Slapaf
  MEP-Search
  nMEP = 10
  MaxStep = 0.1
>> EndDo
&GRID_IT
  All

@Computation of Spectrum
*
* Compute a spectrum
*
&GATEWAY
&SEWARD
&RASSCF
  FileOrb = _Start_Orbitals
&CASPT2
>> Copy $Project.JobMix JOB001
&RASSI
  EJob

