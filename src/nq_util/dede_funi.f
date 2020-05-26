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
      Subroutine DeDe_Funi(Dens,nDens,nr_of_Densities,mDens,ipDq)
      use k2_arrays
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "setup.fh"
      Real*8 Dens(nDens,nr_of_Densities)
      Logical Special_NoSym, DFT_Storage
*
      nIndij=nShlls*(nShlls+1)/2
      nField=2+nr_of_Densities
*
*
      Call mma_allocate(ipOffD,nField,nIndij,label='ipOffD')
      Call GetMem('DeDe2','Allo','Real',ipDeDe,nDeDe_DFT+MaxDe*MaxDCR)
      ipD00 = ipDeDe + nDeDe_DFT
      Call FZero(Work(ipD00),MaxDe*MaxDCR)
*
      Special_NoSym=.False.
      DFT_Storage=.True.
      Call DeDe(Dens,nDens,nr_of_Densities,ipOffD,nIndij,ipDeDe,
     &          ipD00,MaxDe,mDeDe,mIndij,Special_NoSym,DFT_Storage,
     &          Work,1)
      If (mDeDe.ne.nDeDe_DFT) Then
         Call WarningMessage(2,'DeDe_Funi: mDeDe.ne.nDeDe_DFT')
         Write (6,*) ' mDeDe =', mDeDe,' nDeDe_DFT =', nDeDe_DFT
         Call Abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      ipDq=ip_of_Work(Dens(1,1))
      mDens=nDens
*
      Return
      End
