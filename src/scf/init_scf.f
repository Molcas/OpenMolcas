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
      use SCF_Arrays
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
      Integer nAsh(8)
*
      nD = 1
      If (iUHF.eq.1) nD=2
*
*---- Clear Dens and TwoHam matrices
      Call FZero(Dens  ,nBT*nD*nDens)
      Call FZero(TwoHam,nBT*nD*nDens)
      Call FZero(Vxc   ,nBT*nD*nDens)
*
*---- Set number of active shells on the RUNFILE to zero
*
      Call ICopy(8,0,0,nAsh,1)
      Call Peek_iScalar('nSym',i)
* PAM Jan 2007 -- deactivated, improper. Fixed in nqutil in
* another way (query rather than get from runfile)
*      Call Put_iArray('nAsh',nAsh,i)
      NACTEL = 0
      Call Put_iScalar('nActel',NACTEL)
*
      Call IniLLs
*     clear MapDns ...
      iZero=0
      Call ICopy(MxKeep,iZero,0,MapDns,1)
*
      Return
      End
