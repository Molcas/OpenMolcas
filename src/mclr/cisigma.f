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
       SubRoutine CISigma(iispin,iCsym,iSSym,Int1,nInt1,Int2s,nInt2s,
     &                    Int2a,nInt2a,ipCI1,ipCI2, Have_2_el)
       use ipPage, only: W
       use Arrays, only: FIMO, KAIN1, KINT2, KINT2A
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
       Real*8, Target:: Int1(nInt1), Int2s(nInt2s), Int2a(nInt2a)
       Logical Have_2_el
       integer kic(2),opout,nbb(8)
       Real*8, Allocatable:: CIDET(:)
*
*      Interface Anders to Jeppe
*      This interface initiates Jeppes common block
*      and will make it easier to use Anders modifications
*      of the CI routines
*
*      OK first tell Jeppe where the integrals are.
*
       !> yma: notice the nconf1 if DMRG

       If (nconf1.eq.0) return
*
*      One electron integrals
*
       KAIN1=>Int1
*
*      Two electron integrals
*      symmetric in perticle one and two
*
       if (ip_of_Work(Int1(1)).eq.ip_of_Work(FIMO(1))) then
        Call icopy(nsym,nbas,1,nbb,1)
       Else
         Call icopy(nsym,norb,1,nbb,1)
       End if
*
       KINT2=>Int2s
*
*      Two electron integrals
*      anti symmetric in perticle one and two
*
       KINT2A=>Int2a
*
       irefsm=iCSym
*
*      Do we have any twoelectron integrals?
*
       If (Have_2_el) Then
         i12=2
       Else
         i12=1
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
#ifdef _MS_
       irc=ipin(ipCI1)
       irc=ipin(ipci2)
       Do i=0,nroots-1
          call dcopy_(nCSF(iCSM),W(ipCI1)%Vec(1+i*ncsf(icsm)),1,
     &                        CIDET,1)
          Call SigmaVec(CIDET,
     &               W(ipci2)%Vec(1+i*ncsf(issm)),kic)
       End Do
#else
       call dcopy_(nCSF(iCSM),W(ipCI1)%Vec,1,CIDET,1)

       Call SigmaVec(CIDET,W(ipci2)%Vec,kic)

#endif
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
