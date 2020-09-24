************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1999, Anders Bernhardsson                              *
*               1999, Coen de Graaf                                    *
************************************************************************
      Subroutine Dens_IF_SCF(COEFF,CMO,Mode)
*
*     A small stupid interface for creating the alpha and beta
*     occupation numbers and corresponding molecular orbitals,
*     used in MOLDEN.
*
*     EAW 990118
*
*     For SCF it gets even more 'stupid', just reading the MO coefficients and
*     dumping it into a large matrix of dimension nTot x nTot
*       (Coen de Graaf, Oct. 99)
*
      Implicit Real*8 (a-h,o-z)

#include "SysDef.fh"
#include "real.fh"
      Real*8 COEFF(*), CMO(*)
      Character Mode*1
      Integer nBas(0:7)
*
*     COEN WANTED IT AS A BLOCKED MATRIX, SO HERE THEY COME...
*
#ifdef _DEBUG_
#endif
      Call Get_iScalar('nSym',nIrrep)
      Call Get_iArray('nBas',nBas,nIrrep)
      nTot=0
      nTot2=0
      Do iS=0,nIrrep-1
         nTot=nTot+nBas(iS)
         nTot2=nTot2+nBas(iS)**2
      End Do
*
      ip1=1
      ip2=1
      If (Mode.eq.'F') call dcopy_(nTot**2,[Zero],0,COEFF,1)
C     Call RecPrt('COEFF',' ',Coeff,nTot,nTot)
      Do iS=0,nIrrep-1
         If (Mode.eq.'B') call dcopy_(nbas(is)**2,[Zero],0,CMO(ip1),1)
C        Call RecPrt('CMO',' ',CMO(ip1),nbas(is),nbas(is))
         Do i=1,nbas(is)
            If (Mode.eq.'F') call dcopy_(nbas(is),CMO(ip1),  1,
     &                                           COEFF(ip2),1)
            If (Mode.eq.'B') call dcopy_(nbas(is),COEFF(ip2),1,
     &                                           CMO(ip1),  1)
            ip1=ip1+nbas(is)
            ip2=ip2+nTot
         End Do
C        Call RecPrt('CMO',' ',CMO(ip1-nbas(is)**2),nbas(is),nbas(is))
         ip2=ip2+nbas(is)
      End Do
C     Call RecPrt('COEFF',' ',Coeff,nTot,nTot)
*
#ifdef _DEBUG_
#endif
      Return
      End
