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
      SubRoutine UpdFck(OneHam,TwoHam,Vxc,nDT,NumDT,Fock,niter_,nD)
************************************************************************
*                                                                      *
*     purpose: Update Fock matrix from actual OneHam & TwoHam matrices *
*                                                                      *
*     input:                                                           *
*       OneHam  : one-electron hamiltonian of length nDT               *
*       TwoHam  : a few last two-electron hamiltonians (nDT,NumDT)     *
*                                                                      *
*     output:                                                          *
*       Fock    : Fock matrix of length nFO (= OneHam+Optimal TwoHam)  *
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
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
*
      Integer nDT,NumDT
      Real*8 OneHam(nDT),TwoHam(nDT,nD,NumDT),Vxc(nDT,nD,NumDT),
     &       Fock(nDT,nD)
*
*define _DEBUG_
#ifdef _DEBUG_
#endif
*
      i2Hm=NumDT
      If (nIter_.eq.1) i2Hm=1
#ifdef _DEBUG_
      Write (6,*) 'i2Hm=',i2Hm
      Call NrmClc(OneHam,nDT,'UpdFck','OneHam')
      Call NrmClc(TwoHam(1,1,i2Hm),nDT*nD,'UpdFck','T in i2Hm')
      Call NrmClc(Fock            ,(nDT)*nD,'UpdFck','Fock')
      Call NrmClc(Vxc(1,1,i2Hm  ),nDT*nD,'UpdFck','V in i2Hm  ')
#endif
      Do iD = 1, nD
*
*        F = h + G(D), where D is the inter-/extra-polated
*        densitity matrix. This part is exact.
*
*        G(D) = Sum_i c_i G(D_i)
*
         Call DZAXPY(nDT,One,TwoHam(1,iD,i2Hm),1,OneHam,1,Fock(1,iD),1)
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
         Call DaXpY_(nDT,One,Vxc(1,iD,i2Hm  ),1,Fock(1,iD),1)
*
#ifdef _DEBUG_
         write(6,*) 'Fock'
         write(6,'(5f12.6)') (Fock(ivv,iD),ivv=1,nDT)
         Call NrmClc(Fock(1,iD),nDT,'Fock  ','UpdFck ')
         Call NrmClc(TwoHam(1,iD,i2Hm),nDT,'TwoHam','UpdFck ')
         Call NrmClc(Vxc(1,iD,i2Hm),nDT,'Vxc','UpdFck ')
#endif
*
      End Do
#ifdef _DEBUG_
#endif
      Return
      End
