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
       SubRoutine CISigma_td(iispin,iCsym,iSSym,Int1,Int2s,
     &                       Int2a,ipCI1,ipCI2,NT)
       Implicit Real*8(a-h,o-z)
c
c For the timeindep case ipS1 and ipS2 will be half as long
c Avoid sigmavec calls. 95% of the time in mclr is spent in sigmavec
c
*
#include "Pointers.fh"
#include "Input.fh"
#include "WrkSpc.fh"
#include "genop.fh"
#include "glbbas_mclr.fh"
#include "lbbas1.fh"
#include "cands.fh"
#include "detdim.fh"
#include "cstate_mclr.fh"
#include "cicisp_mclr.fh"
       Character NT
       integer kic(2),opout
       Real*8 Int1(*), Int2s(*), Int2a(*)
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
       ipInt1 =ip_of_Work(Int1(1))
       ipInt2s=ip_of_Work(Int2s(1))
       ipInt2a=ip_of_Work(Int2a(1))
*
*      One electron integrals
*
       KAIN1=ipInt1
*
*      Two electron integrals
*      symmetric in perticle one and two
*
*
       KINT2= ipInt2s
       KINT2a=ipint2a
*
*      Two electron integrals
*      anti symmetric in perticle one and two
*
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
       If (TIMEDEP) Then
        If (NT.eq.'T')  square=.true.  ! The operator is not sym
        If (.not.page) Then
c ipcidet is here because sigmavec will destroy the first input vector.
         Call GetMem('CIDET_TD','ALLO','REAL',ipCIDET,nDet)
             call dcopy_(nCSF(iCSM),Work(ipin(ipCI1)),1,Work(ipCIDET),1)
C
           Call SigmaVec(Work(ipCIDET),Work(ipin(ipci2)),kic)
C
         If (NT.eq.'N') Then
            Call GetMem('CIDET_TD','FREE','REAL',ipCIDET,nDet)
            Return
         End If
c
         If (NT.eq.'S') Then
C......... Symmetric operator, no transpose of integrals needed!
           call dcopy_(nCSF(iCSM),Work(ipin(ipCI1)+nConf1),1,
     &     Work(ipCIDET),1)
C
         Else
C......... The operator is not sym --> transpose integrals! NT.ne.S
           call dcopy_(nCSF(iCSM),Work(ipin(ipCI1)),1,Work(ipCIDET),1)
           Call GetMem('TEMPINT1','ALLO','REAL',ipTI1,ndens2)
           Call GetMem('TEMPINT2','ALLO','REAL',ipTI2,ntash**4)
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
                 Work(ipTI2+jilk-1)=int2s(ijkl)
                End if
               End Do
              End Do
             End Do
           End Do
           Do is=1,nSym
            js=ieor(ieor(icsym-1,issym-1),is-1)+1
            If (nbas(js)*nbas(is).ne.0)
     &      Call DGETMO(Int1(ipmat(is,js)),nbas(is),
     &                nbas(is),nbas(js),Work(ipTI1+ipmat(js,is)-1),
     &                nbas(js))
           End Do
           kain1=ipTI1
           KINT2=ipTI2
         End If  ! End the transpose of integrals.
*
*
C
            Call SigmaVec(Work(ipCIDET),Work(ipin(ipci2)+nconf1),kic)
C
         Call GetMem('CIDET_TD','FREE','REAL',ipCIDET,nDet)
         If (NT.ne.'S') Then
           Call GetMem('TEMPINT1','FREE','REAL',ipTI1,ndens2)
           Call GetMem('TEMPINT2','FREE','REAL',ipTI2,ntash**4)
         End If
        Else   ! If page
         Write(6,*) 'Page not implemented for Timedependent'
     &            //' perturbations'
         Call Abend()
        End If
       Else   ! If not timedep
        If (.not.page) Then
         Call GetMem('CIDET_TD','ALLO','REAL',ipCIDET,nDet)
         call dcopy_(nCSF(iCSM),Work(ipin(ipCI1)),1,Work(ipCIDET),1)
         Call SigmaVec(Work(ipCIDET),Work(ipin(ipci2)),kic)
         Call GetMem('CIDET_TD','FREE','REAL',ipCIDET,ndet)
        Else
         irc=ipnout(ipci2)
         Call SigmaVec(Work(ipin1(ipCI1,ndet)),Work(ipin(ipci2)),kic)
         irc=opout(ipci1)
        End If
       End If
*
       Return
       End
