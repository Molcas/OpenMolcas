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
  Basis (XYZ) = ANO-RCC-MB
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
>> shell pegamoid.py $Project.guessorb.h5 $Project.LocOrb
&RASSCF
  FileOrb = $Project.LocOrb
>> shell pegamoid.py $Project.guessorb.h5 $Project.RasOrb

@Localized Atomic Orbitals for CASSCF/RASSCF
*
* Generating Starting Orbitals for CASSCF/RASSCF
*
&GATEWAY
  Basis (XYZ) = ANO-RCC-MB
  Group = Full
  RICD
&LOCALISATION
  FileOrb = $Project.GssOrb
  All
>> shell pegamoid.py $Project.guessorb.h5 $Project.LocOrb
&RASSCF
  FileOrb = $Project.LocOrb
>> shell pegamoid.py $Project.guessorb.h5 $Project.RasOrb

@SCF Orbitals for CASSCF/RASSCF
*
* Generating SCF Starting Orbitals for CASSCF/RASSCF in ANO-RCC-MB
*
&GATEWAY
  Basis (XYZ) = ANO-RCC-MB
  Group = Full
  RICD
&SEWARD
&SCF
>> shell pegamoid.py $Project.guessorb.h5 $Project.ScfOrb
&RASSCF
  FileOrb = $Project.ScfOrb
>> shell pegamoid.py $Project.guessorb.h5 $Project.RasOrb

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
  Basis (XYZ) = ANO-RCC-MB
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
>> shell pegamoid.py $Project.guessorb.h5 $Project.LocOrb
&RASSCF
  FileOrb = $Project.LocOrb
>> Copy $Project.RunFile RUNFIL1
&GATEWAY
  Basis (XYZ) = ANO-RCC-VDZP
  Group = NoSymm
  RICD
>> Copy $Project.RunFile RUNFIL2
&EXPBAS
  FileOrb = $Project.RasOrb
&SEWARD
&RASSCF
  FileOrb = $Project.ExpOrb

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

