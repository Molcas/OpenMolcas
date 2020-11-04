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
       SubRoutine CISigma_td(iispin,iCsym,iSSym,Int1,nInt1,Int2s,nInt2s,
     &                       Int2a,nInt2a,ipCI1,ipCI2,NT, Have_2_el )
       use ipPage, only: W
       use Arrays, only: KAIN1, KINT2, KINT2A, TI1, TI2, pInt1
       Implicit Real*8(a-h,o-z)
c
c For the timeindep case ipS1 and ipS2 will be half as long
c Avoid sigmavec calls. 95% of the time in mclr is spent in sigmavec
c
*
#include "Pointers.fh"
#include "Input.fh"
#include "stdalloc.fh"
#include "genop.fh"
#include "cands.fh"
#include "detdim.fh"
#include "cstate_mclr.fh"
#include "cicisp_mclr.fh"
       Character NT
       integer kic(2),opout
       Logical Have_2_el
       Real*8, Target:: Int1(nInt1), Int2s(nInt2s), Int2a(nInt2a)
       Real*8, Allocatable:: CIDET(:)

       itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*
*      Interface Anders to Jeppe
*      This interface initiates Jeppes common block
*      and will make it easier to use Anders modifications
*      of the CI routines
*
*      OK first tell Jeppe where the integrals are.
*
       If (nconf1.eq.0) return
*
*      One electron integrals
*
       KAIN1=>Int1
*
*      Two electron integrals
*      symmetric in perticle one and two
*
*
       KINT2 => Int2s
       KINT2a=> Int2a
*
*      Two electron integrals
*      anti symmetric in perticle one and two
*
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
C
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
         Call iCopy(nSym,ipCM,1,pInt1,1)
       Else
         Do iS=1,nSym
          jS=iEor(iS-1,iOp-1)+1
          pInt1(is)=ipMat(is,jS)
         End Do
       End If
*
*      Triplet/Singlet operator
*
       ist=iispin+1
       square=.false.
*                                                                      *
************************************************************************
*                                                                      *
       If (TIMEDEP) Then
*                                                                      *
************************************************************************
*                                                                      *
          If (NT.eq.'T')  square=.true.  ! The operator is not sym


          If (page) Then
             Write(6,*) 'Page not implemented for Timedependent'
     &              //' perturbations'
             Call Abend()
          End If

c         CIDET is here because sigmavec will destroy the first
C         input vector.
          Call mma_allocate(CIDET,nDet,Label='CIDET')
          irc=ipin(ipCI1)
          call dcopy_(nCSF(iCSM),W(ipCI1)%Vec,1,CIDET,1)

          irc=ipin(ipci2)
          Call SigmaVec(CIDET,W(ipci2)%Vec,kic)
C
          If (NT.eq.'N') Then
             Call mma_deallocate(CIDET)
             Return
          End If
c
          If (NT.eq.'S') Then

C.......... Symmetric operator, no transpose of integrals needed!
            irc=ipin(ipCI1)
            call dcopy_(nCSF(iCSM),W(ipCI1)%Vec(1+nConf1),1,
     &                  CIDET,1)
C
            irc=ipin(ipci2)
            Call SigmaVec(CIDET,W(ipci2)%Vec(1+nconf1),kic)

          Else  ! NT.ne.'S'

C.......... The operator is not sym --> transpose integrals! NT.ne.S
            irc=ipin(ipCI1)
            call dcopy_(nCSF(iCSM),W(ipCI1)%Vec,1,CIDET,1)

            Call mma_allocate(TI1,ndens2,Label='TI1')
            Call mma_allocate(TI2,ntash**4,Label='TI2')

            Do i=1,ntash
              Do j=1,ntash
               ij=i+ntash*(j-1)
               ji=j+ntash*(i-1)
               Do k=1,ntash
                Do l=1,ntash
                 kl=k+ntash*(l-1)
                 lk=l+ntash*(k-1)
                 If (ij.ge.kl) Then
                  ijkl=itri(ij,kl)
                  jilk=itri(ji,lk)
                  TI2(jilk)=int2s(ijkl)
                 End if
                End Do
               End Do
              End Do
            End Do

            Do is=1,nSym
             js=ieor(ieor(icsym-1,issym-1),is-1)+1
             If (nbas(js)*nbas(is).ne.0)
     &       Call DGETMO(Int1(ipmat(is,js)),nbas(is),
     &                 nbas(is),nbas(js),TI1(ipmat(js,is)),
     &                 nbas(js))
            End Do

            KAIN1=>TI1
            KINT2=>TI2

            irc=ipin(ipci2)
            Call SigmaVec(CIDET,W(ipci2)%Vec(1+nconf1),kic)

            KAIN1=>Null()
            KINT2=>Null()
            Call mma_deallocate(TI1)
            Call mma_deallocate(TI2)

         End If  ! End the transpose of integrals.
*
         Call mma_deallocate(CIDET)
*                                                                      *
************************************************************************
*                                                                      *
       Else   ! If not timedep
*                                                                      *
************************************************************************
*                                                                      *

          If (.not.page) Then
             Call mma_allocate(CIDET,nDet,Label='CIDET')
             irc=ipin(ipCI1)
             call dcopy_(nCSF(iCSM),W(ipCI1)%Vec,1,CIDET,1)
             irc=ipin(ipci2)
             Call SigmaVec(CIDET,W(ipci2)%Vec,kic)
             Call mma_deallocate(CIDET)
          Else
             irc=ipnout(ipci2)
             irc=ipin1(ipCI1,ndet)
             irc=ipin(ipci2)
             Call SigmaVec(W(ipCI1)%Vec,W(ipci2)%Vec,kic)
             irc=opout(ipci1)
          End If
*                                                                      *
************************************************************************
*                                                                      *
       End If
*                                                                      *
************************************************************************
*                                                                      *
*
       Return
       End
