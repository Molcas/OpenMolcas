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
      Subroutine Init_SCF()
      use SCF_Arrays, only: Dens, TwoHam, Vxc
      use InfSCF, only: iUHF, MapDns, nBT, nDens, Two_Thresholds
      use MxDM, only: MxKeep
      use RICD_Info, only: Do_DCCD
      use Constants, only: Zero
      Implicit None
      Integer i, iZero, nActEl, nD
      Integer nAsh(8)
*
      nD = 1
      If (iUHF.eq.1) nD=2
*
*---- Clear Dens and TwoHam matrices
      Dens  (:,:,:)=Zero
      TwoHam(:,:,:)=Zero
      Vxc   (:,:,:)=Zero
*
*---- Set number of active shells on the RUNFILE to zero
*
      nAsh(:)=0
      Call Peek_iScalar('nSym',i)
* PAM Jan 2007 -- deactivated, improper. Fixed in nqutil in
* another way (query rather than get from runfile)
*     Call Put_iArray('nAsh',nAsh,i)
      NACTEL = 0
      Call Put_iScalar('nActel',NACTEL)
*
      Call IniLLs()
*     clear MapDns ...
      MapDns(:)=0
*
      Two_Thresholds=.NOT.Do_DCCD
*
      Return
      End
