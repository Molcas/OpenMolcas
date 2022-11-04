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
* Copyright (C) 1995, Martin Schuetz                                   *
************************************************************************
      SubRoutine UpdFck(nIter_)
************************************************************************
*                                                                      *
*     purpose: Update Fock matrix from actual OneHam & TwoHam matrices *
*                                                                      *
*     called from: WfCtl                                               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.Schuetz                                                        *
*     University of Lund, Sweden, 1995                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Use SCF_Arrays, only: OneHam, TwoHam, Vxc, FockAO
      Implicit None
#include "real.fh"
#include "mxdm.fh"
*
      Integer nIter_

      Integer nD,nDT,NumDT,iD,i2Hm
*
      NumDT=Size(TwoHam,3)
      nDT  =Size(FockAO,1)
      nD   =Size(FockAO,2)

      i2Hm=NumDT
      If (nIter_.eq.1) i2Hm=1
#ifdef _DEBUGPRINT_
      Write (6,*) 'i2Hm=',i2Hm
      Call NrmClc(OneHam,nDT,'UpdFck','OneHam')
      Call NrmClc(TwoHam(1,1,i2Hm),nDT*nD,'UpdFck','T in i2Hm')
      Call NrmClc(FockAO          ,(nDT)*nD,'UpdFck','FockAO')
      Call NrmClc(Vxc(1,1,i2Hm  ),nDT*nD,'UpdFck','V in i2Hm  ')
#endif
      Do iD = 1, nD
*
*        F = h + G(D), where D is the inter-/extra-polated
*        densitity matrix. This part is exact.
*
*        G(D) = Sum_i c_i G(D_i)
*
         Call DZAXPY(nDT,One,TwoHam(1,iD,i2Hm),1,OneHam,1,
     &               FockAO(1,iD),1)
*
*        Here we add the contribution to the Fock matrix from
*        the external field (the none linear or bi-linear parts).
*        Note that Vxc(D_a+D_b)=/= Vxc(D_a) + Vxc(D_b). However,
*        using a first order Taylor expansion we can demonstate that
*        this is correct to first order. Hence, we pick up the
*        inter-/extra-polated external potential here as generated
*        by OptClc.f.
*
*        Proof:
*
*        D = Sum_i c_i * D_i; Sum_i c_i = 1
*
*        First order Taylor expansion of Vxc around D_0 (arbitrary)
*
*        Vxc(D) = Vxc(D_0) + dVxc(D_0)/dD (D-D_0)
*
*        Vxc(D) = Vxc(D_0) + Sum_i C_i dVxc(D_0)/dD (D_i_D_0)
*
*        Vxc(D) = Vxc(D_0) + Sum_i C_i (Vxc(D_i) - Vxc(D_0)
*
*        Vxc(D) = Sum_i c_i Vxc(D_i), which is produced by OptClc.f
*
         Call DaXpY_(nDT,One,Vxc(1,iD,i2Hm  ),1,FockAO(1,iD),1)
*
#ifdef _DEBUGPRINT_
         write(6,*) 'Fock'
         write(6,'(5f12.6)') (FockAO(ivv,iD),ivv=1,nDT)
         Call NrmClc(FockAO(1,iD),nDT,'Fock  ','UpdFck ')
         Call NrmClc(TwoHam(1,iD,i2Hm),nDT,'TwoHam','UpdFck ')
         Call NrmClc(Vxc(1,iD,i2Hm),nDT,'Vxc','UpdFck ')
#endif
*
      End Do

      End
