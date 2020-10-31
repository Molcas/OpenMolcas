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
       SubRoutine CISigma_sa(iispin,iCsym,iSSym,Int1,int2s,
     &                       Int2a,ipCI1,ipCI2)
       use ipPage, only: W
       Implicit Real*8(a-h,o-z)
*
#include "Pointers.fh"

#include "Input.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "genop.fh"
#include "glbbas_mclr.fh"
#include "lbbas1.fh"
#include "cands.fh"
#include "detdim.fh"
#include "cstate_mclr.fh"

#include "cicisp_mclr.fh"
       integer kic(2),opout,nbb(8)
       Real*8 Int1(*), Int2s(*), Int2a(*)
       Real*8, Allocatable:: CIDET(:)
*
*      Interface Anders to Jeppe
*      This interface initiates Jeppes common block
*      and will make it easier to use Anders modifications
*      of the CI routines
*
*      OK first tell Jeppe where the integrals are.

       ipInt1 = ip_of_Work(Int1(1))
       ipInt2s= ip_of_Work(Int2s(1))
       ipInt2a= ip_of_Work(Int2a(1))
       !> yma: notice the nconf1 if DMRG
*
       If (nconf1.eq.0) return
*
*      One electron integrals
*
       KAIN1=ipInt1
*
*      Two electron integrals
*      symmetric in perticle one and two
       if (ipint1.eq.ipfimo) then
        Call icopy(nsym,nbas,1,nbb,1)
       Else
         Call icopy(nsym,norb,1,nbb,1)
       End if
*
       KINT2=ipint2s
*
*      Two electron integrals
*      anti symmetric in perticle one and two
*
       KINT2a=ipint2a
*
       irefsm=iCSym
*
*      Do we have any twoelectron integrals?
*
       If (ipInt2s.eq.0) Then
         i12=1
       Else
         i12=2
       End If
*
*      Symmetry of Sigma vector
*
       iSSM=iSSym
       kic(2)=2
       if (issm.eq.State_sym) kic(2)=1
*
*      Symmetry of CI vector
*
       iCSM=iCSym
       kic(1)=2
       if (icsm.eq.State_sym) kic(1)=1
*
*      Symmetry properties of operator
*
       ndet=nint(max(xispsm(iSSym,1),xispsm(iCSym,1)))
       ndet=Max(ndet,ncsf(icsym),ncsf(issym))

       If (ndet.eq.0) Return
       iOP=iEOr(iCSM-1,iSSm-1)+1
       If (iOp.eq.1) Then
         Call iCopy(nSym,ipCM,1,iWork(kapin1),1)
       Else
         Do iS=1,nSym
          jS=iEor(iS-1,iOp-1)+1
          iWork(kapin1+is-1)=ipMat(is,jS)
         End Do
       End If
*
*      Triplet/Singlet operator
*
       ist=iispin+1
       square=.false.
*
       If (.not.page) Then
          Call mma_allocate(CIDET,nDet,Label='CIDET')
          irc=ipin(ipCI1)
          irc=ipin(ipci2)
          Do i=0,nroots-1
             call dcopy_(nCSF(iCSM),W(ipCI1)%Vec(1+i*ncsf(icsm)),1,
     &                   CIDET,1)
             Call SigmaVec(CIDET,W(ipci2)%Vec(1+i*ncsf(issm)),
     &                     kic)
             Call DSCAL_(nCSF(iCSM),weight(i+1),
     &                   W(ipci2)%Vec(1+i*ncsf(issm)),1)
          End Do
          Call mma_deallocate(CIDET)
       Else
          irc=ipnout(ipci2)
          irc=ipin1(ipCI1,ndet)
          irc=ipin(ipci2)
          Call SigmaVec(W(ipCI1)%Vec,W(ipci2)%Vec,kic)
          irc=opout(ipci1)
       End If
*
       Return
       End
