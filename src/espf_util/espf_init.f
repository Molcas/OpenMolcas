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
      Subroutine espf_init (natom,nAtQM,ipCord,ipIsMM,ipExt)
      Implicit Real*8 (A-H,O-Z)
*
*     ESPF initialization:
*       natom: number of atoms
*       ipCord: pointer to the atom coordinates
*       ipIsMM: pointer to the "logical" MM or QM array
*       ipExt: pointer to the external potential array
*
#include "espf.fh"
*
      Call QEnter('espf_init')
*
      Call Get_iScalar('Unique atoms',natom)
      Call GetMem('AtomCoord','Allo','Real',ipCord,3*natom)
      Call Get_dArray('Unique Coordinates',Work(ipCord),3*natom)
      Call MMCount(natom,nAtMM,ipIsMM)
      nAtQM = natom - nAtMM
      Call GetMem('ExtPot','ALLO','REAL',ipExt,natom*MxExtPotComp)
      call dcopy_(MxExtPotComp*natom,[Zero],0,Work(ipExt),1)
*
      Call QExit('espf_init')
      iReturn=0
      Return
      End
