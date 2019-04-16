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
* Copyright (C) Jonas Bostrom                                          *
************************************************************************
      Subroutine Compute_A_jk_mp2(iSO,jVec,kVec,Ajk,fac_ij,fac_kl,
     &                        nVec,iOpt)
**************************************************************************
*     Author: J Bostrom
*
*     Purpose: Loading A-matrix for mp2 from disk
*
**************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "exterm.fh"
#include "chomp2g_alaska.fh"
#include "WrkSpc.fh"
      Real*8  Ajk,Fac_ij, Fac_kl, dum(1)

      Character*16 SECNAM
      Parameter (SECNAM = 'Compute_A_jk_mp2')

      Ajk = 0.0d0
      If(imp2prpt .eq. 2) Then
         lTot = 1
         iAdrA = nVec*(kVec-1) + jVec
         Call dDaFile(LuAVector(iOpt),2,dum,lTot,iAdrA)
         Ajk_mp2=dum(1)
         Ajk = Ajk + (Ajk_mp2*Fac_kl*Fac_ij)
      End If

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iSO)
      End
