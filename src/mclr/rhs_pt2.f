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
* Copyright (C) 1998, Anders Bernhardsson                              *
************************************************************************
#define NOCODE
#ifdef NOCODE
      Subroutine RHS_PT2(rkappa,iprci)
      Implicit Real*8(a-h,o-z)
      Real*8 rKappa(*)
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(rkappa)
         Call Unused_integer(iprci)
      End If
#else
      Subroutine RHS_PT2(rkappa,iprci)
      use Arrays, only: CMO, Int1
      use Str_Info, only: DTOC, CNSM
      use ipPage, only: W
      Implicit Real*8(a-h,o-z)
      Real*8 rKappa(*)
*
*    Calculate and reads in from CASPT2 the RHS needed
*    for the calculation of the Lagrangian multipliers
*    for CASPT2.
*
*    From the same Fock matrix the effective Fock matrix
*    needed for the renormalization term is  calculated.

*    Here is the structure, it is not debugged and one
*    needs to check the detail, but all "input" is there.
*
#include "Pointers.fh"

#include "Input.fh"
#include "stdalloc.fh"
#include "glbbas_mclr.fh"
#include "Files_mclr.fh"
#ifndef _DEBUGED_
#include "detdim.fh"
#endif
      Real*8, Allocatable:: DCAS(:), TempK(:), TempCI(:), TempCI2(:),
     &                      T2(:), FAO1(:), FAO2(:), FMO1(:), FMO2(:),
     &                      DP(:), DP2(:), DCAS2(:), K1(:), K2(:)
      Real*8 rDum(1)
      Half=0.5d0
*
*     Read in a and b part of effective gradient from CASPT2
*                                                                     *
***********************************************************************
*                                                                     *
      Call mma_allocate(TempK,ndens2,Label='TempK')
      If (imethod.eq.2)Then
         i=0
         Call mma_allocate(TempCI,nconf1,Label='TempCI')
         Call mma_allocate(TempCI2,nconf1,Label='TempCI2')
         Call dDaFile(LuPT2,2,TempK,ndens2,i)
#ifdef _DEBUGED_
*
*    CASPT2 mode
*
         Call dDaFile(LuPT2,2,TempCI,nconf1,i)
#else
*
* lucia mode
*
         n=nint(xispsm(State_SYM,1))
         Call dDaFile(LuPT2,2,TempCI,n,i)
#endif
*
*---- Transform from split GUGA to symmetric group
*
#ifdef _DEBUGED_
*
*    CASPT2 mode (split graph)
*
         Call Gugactl_MCLR(TempCI(1),1)
#else
*
* lucia mode (Symmetric group)
*
      iprdia=0
      Call INCSFSD(STATE_SYM,STATE_SYM,.false.)
      CALL CSDTVC_MCLR(TempCI2,TempCI,2,
     &                 TOC,CNSM(1)%ICTS,
     &                 State_SYM,1,IPRDIA)
#endif
      Call mma_deallocate(TempCI2)
*                                                                     *
***********************************************************************
*                                                                     *
      Else
*                                                                     *
***********************************************************************
*                                                                     *
*     MP2
*
        Call mma_allocate(TempCI,1,Label='TempCI')
        i=0
        Call dDaFile(LuPT2,2,TempK,ndens2,i)
*       Call ThreeP(Kappa)
*       MP2
*                                                                     *
***********************************************************************
*                                                                     *
      End if
* ----
*
      Call mma_allocate(T2,ndens2,Label='T2')
      Call mma_allocate(FAO1,ndens2,Label='FAO1')
      Call mma_allocate(FAO2,ndens2,Label='FAO2')
      Call mma_allocate(FMO1,ndens2,Label='FMO1')
      Call mma_allocate(FMO2,ndens2,Label='FMO2')
      Call mma_allocate(DP,ndens2,Label='DP')
      Call mma_allocate(DP2,ndens2,Label='DP2')
      Call mma_allocate(DCAS2,ndens2,Label='DCAS2')
      Call mma_allocate(K1,ndens2,Label='K1')
      Call mma_allocate(K2,ndens2,Label='K2')
*
*---  Read in necessary densities.
*
      Call Qpg_dArray('D1ao',Found,nDens)
      If (Found .and. nDens/=0) Then
         Call mma_allocate(DCAS,nDens,Label='DCAS*)
      Else
         Write (6,*) 'RHS_PT2: Density not found'
         Call Abend()
      End If
      Call Get_D1ao(DCAS,nDens)
      irc=-1
      iopt=0
      Call RdRlx(irc,iopt,'D1PT22',DP)
      If (irc.ne.0) Goto 100
      irc=-1
      iopt=0
      Call RdRlx(irc,iopt,'OVLP',rovlp)
      If (irc.ne.0) Goto 100
*
*--- Squared density
*
      Call UnFold_MCLR(DP,DP2)
*
      Call UnFold_MCLR(DCAS,DCAS2)
*
*======================================================================*
*
*
*---  Make Fock matrixes
*
      ExFac=1.0D0
*
*  1) P(PT2)
*
      nFlt=0
      nBMX=0
      Do iSym = 1, nSym
         nFlt=nFlt+nBas(iSym)*(nBas(iSym)+1)/2
         nBMX=Max(nBMX,nBas(iSym))
      End Do
*
      FAO1(:)=0.0D0
      Call FockTwo_Drv(nSym,nBas,nBas,nSkip,
     &                 DP,DP2),FAO1,nFlt,
     &                 ExFac,nDens2,nBMX)
      Call AO2MO(FAO1,FMO1)
*
*  2) P(CAS)
*
      FAO2(:)=0.0D0
      Call FockTwo_Drv(nSym,nBas,nBas,nSkip,
     &                 DAS,DCAS2,FAO2,nFlt,
     &                 ExFac,nDens2,nBMX)
      Call AO2MO(FAO2,FMO2)
*
*-- Add one particle hamiltonian here ???
*
  ??? Call DaXpY_(ndens2,-rovlp,Int1,1,FMO1,1)
*
*======================================================================*
*
*---- CI Vector
*
*     <i|Sigma> = <i|F|0> - <0|F|0><i|0>+CI_a+CI_b
*
      Call CISigma(0,State_sym,State_sym,FMO1,nDens2,rdum,1,rdum,ipci,
     &             iprci,.False.)
      irc=ipin(iprci)
      irc=ipin(ipci)
      rE=ddot_(nconf1,W(iprci)%Vec,1,W(ipci)%Vec,1)
      Call Daxpy_(nconf1,1.0d0,TempCI,1,W(iprci)%Vec,1)
      Call Daxpy_(nconf1,-rE,W(ipCI)%Vec,1,W(iprci)%Vec,1)
*==============================================================================*
*
*     D^CAS*P(PT2+CAS)
*
      Call DFock(DCAS2,FMO1,K1)
*
*     D^(PT2+CAS)*P(CAS)
*
      Call DFock(DP,FMO2,K2)
*
*==============================================================================*
*
*     Add all orbital terms together
*
      Call Daxpy_(nDens2,1.0d0,K2,1,K1,1)
      call daxpy_(ndens2,1.0d0,TempK,1,K1,1)
*
*==============================================================================*
*
*  OKIDOKI Two things are needed, first of all the symmetric fockmatrix for the
*  connection term:
*
*==============================================================================*
*
*---  Calculate efficent Fock matrix in AO basis (contravariant,folded)
*     and write to disk. Woba
*
      Do iS=1,nSym
        If (nbas(is).ne.0)
     *  Call DGEadd(K1(ipMat(is,is)),nBas(is),'N',
     *              K1(ipMat(is,is)),nBas(is),'T',
     *              K2(ipMat(is,is)),nBas(is),
     *              nBas(is),nBas(is))
      End Do
      Do iS=1,nSym
        If (nBas(is).ne.0) Then
           Call DGEMM_('N','N',
     &                 nBas(iS),nBas(iS),nBas(iS),
     &                 1.0d0,CMO(ipCM(iS)),nBas(iS),
     &                 K2(ipMat(is,is)),nBas(iS),
     &                 0.0d0,T2,nBas(iS))
           Call DGEMM_('N','T',
     &                 nBas(is),nBas(iS),nBAs(iS),
     &                 1.0d0,T2,nBas(iS),
     &                 CMO(ipCM(iS)),nBas(iS),
     &                 0.0d0,FAO1(ipMat(iS,iS)),nBas(is))
        End If
      End Do
      Call Fold2(nsym,nbas,FAO1,FAO2)
      Call Put_Fock_Occ(FAO2,nDens2)
*
*     And then I want the unsymmetric part to the MCLR
*
*
*---  Keep gradient terms woba!!!!!!!
*
      Do iS=1,nSym
        If (nbas(is).ne.0)
     *  Call DGESUB(K1(ipMat(is,is)),nBas(is),'N',
     *              K1(ipMat(is,is)),nBas(is),'T',
     *              rKappa(ipMat(is,is)),nBas(is),
     *              nBas(is),nBas(is))
      End Do

*
*    OK Thats all folks!!
*
      Call mma_deallocate(T2)
      Call mma_deallocate(FAO1)
      Call mma_deallocate(FAO2)
      Call mma_deallocate(FMO1)
      Call mma_deallocate(FMO2)
      Call mma_deallocate(DP)
      Call mma_deallocate(DP2)
      Call mma_deallocate(DCAS)
      Call mma_deallocate(DCAS2)
      Call mma_deallocate(K1)
      Call mma_deallocate(K2)
      Call mma_deallocate(TempK)
      Call mma_deallocate(TempCI)
*
*............
      return
 100  Call SysHalt('rhs_pt2')
      end
      Subroutine DFock(DAO,FockMO,Fock)
      use Arrays, only: CMO
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"

#include "Input.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Real*8 DAO(*),FockMO(*),Fock(*)
      Real*8, Allocatable:: Temp1(:), Temp2(:)

      Call mma_allocate(Temp1,ndens2,Label='Temp1')
      Call mma_allocate(Temp2,ndens2,Label='Temp2')
      Do iS=1,nSym
        If (nBas(is).ne.0) Then
           Call DGEMM_('T','N',
     &                 nBas(iS),nBas(iS),nBas(iS),
     &                 1.0d0,CMO(ipCM(iS)),nBas(iS),
     &                 DAO(ipCM(is)),nBas(iS),
     &                 0.0d0,Temp2,nBas(iS))
           Call DGEMM_('N','N',
     &                 nBas(is),nBas(iS),nBAs(iS),
     &                 1.0d0,Temp2,nBas(iS),
     &                 CMO(ipCM(iS)),nBas(iS),
     &                 0.0d0,Temp1,nBas(is))
           Call DGEMM_('N','N',
     &                 nBas(is),nBas(iS),nBAs(iS),
     &                 1.0d0,Temp1,nBas(iS),
     &                 FockMO(ipCM(iS)),nBas(iS),
     &                 0.0d0,Fock(ipCM(is)),nBas(is))
        End If
      End Do
      Call mma_deallocate(Temp2)
      Call mma_deallocate(Temp1)
      Return
      End
      Subroutine AO2MO(FAO ,FMO)
      use Arrays, only: CMO
      Implicit Real*8 (a-h,o-z)
      Real*8 FAO(*),FMO(*)

#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Real*8, Allocatable:: Temp1(:), Temp2(:)

      Call mma_allocate(Temp1,ndens2,Label='Temp1')
      Call mma_allocate(Temp2,ndens2,Label='Temp2')
      ip=1
      Do iS=1,nSym
        If (nBas(is).ne.0) Then
           Call Square(FAO(ip),Temp1,1,nBas(is),nBas(is))
           Call DGEMM_('T','N',
     &                 nBas(iS),nBas(iS),nBas(iS),
     &                 1.0d0,CMO(ipCM(iS)),nBas(iS),
     &                       Temp1,nBas(iS),
     &                 0.0d0,Temp2,nBas(iS))
           Call DGEMM_('N','N',
     &                 nBas(is),nBas(iS),nBAs(iS),
     &                 1.0d0,Temp2,nBas(iS),
     &                       CMO(ipCM(iS)),nBas(iS),
     &                 0.0d0,FMO(ipMat(iS,iS)),nBas(is))
           ip=ip+nBas(is)*(nBas(iS)+1)/2
        End If
      End Do
      Call mma_deallocate(Temp2)
      Call mma_deallocate(Temp1)
#endif
      Return
      End
