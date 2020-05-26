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
      Subroutine DeDe_FAIEMP( Dens, TwoHam, nDens, mDens, ipDq, ipFq)
      use k2_arrays
      Implicit None
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "setup.fh"
      Integer  nDens, mDens, nr_of_Densities, ipDq, ipFq
      Real*8   Dens(nDens), TwoHam(nDens)
      Logical  Special_NoSym, DFT_Storage
* local variables
      Integer  i, ij, mDeDe, nIndij, nField, mIndij
      Integer  ip_of_work
      External ip_of_work

#ifdef _DEBUG_
      Call qEnter('DeDe_FAIEMP')
#endif
      nr_of_Densities=1  ! Hardwired option
*
      nIndij=nShlls*(nShlls+1)/2
      nField=2+nr_of_Densities
      Call mma_allocate(ipOffD,nField,nIndij,label='ipOffD')
      Call GetMem('DeDe2','Allo','Real',ipDeDe,nDeDe+MaxDe*MaxDCR)
      ipD00=ipDeDe+nDeDe
      call dcopy_(MaxDe*MaxDCR,[Zero],0,Work(ipD00),1)
*
      Special_NoSym=.True.
      DFT_Storage=.False.
      Call DeDe(Dens,nDens,nr_of_Densities,ipOffD,nIndij,ipDeDe,
     &          ipD00,MaxDe,mDeDe,mIndij,Special_NoSym,DFT_Storage,
     &          Work,1)
      If (mDeDe.ne.nDeDe) Then
c         Write (6,*) ' mDeDe =', mDeDe,' nDeDe =', nDeDe
         Call ErrTra
         Call Abend
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     In case of no symmetry do a temporary square copy of the
*     density matrix.
*
*---- Change the folded density to pure triangular form,
*     i.e. off-diagonal elements are divided by two.
*
      If (nIrrep.eq.1) Then
         Call DScal_(nDens,Half,Dens,1)
         ij=0
         Do i = 1, nBas(0)
            ij = ij + i
            Dens(ij)=Two*Dens(ij)
         End Do
         mDens=nbas(0)*nbas(0)
         Call GetMem('DENQ','Allo','Real',ipDq,mDens)
         Call GetMem('FMAQ','Allo','Real',ipFq,mDens)
         Call Square(Dens,Work(ipDq),1,nbas(0),nbas(0))
         Call fzero(work(ipFq),mDens)
      Else
         ipDq=ip_of_Work(Dens(1))
         ipFq=ip_of_Work(TwoHam(1))
         mDens=nDens
      End If
#ifdef _DEBUG_
      Call qExit('DeDe_FAIEMP')
#endif
*
      Return
      End
