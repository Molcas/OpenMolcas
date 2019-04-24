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
       SubRoutine CISigma_td(iispin,iCsym,iSSym,ipInt1,ipint2s,
     &                    ipint2a,ipCI1,ipCI2,NT)
       Implicit Real*8(a-h,o-z)
c
c For the timeindep case ipS1 and ipS2 will be halv as long
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
       KAIN1=ipInt1
*
*      Two electron integrals
*      symmetric in perticle one and two
*
*
C
*       Write(*,*)'ipint1,ipint2s, ipint2a',ipint1,ipint2s, ipint2a
*       Call RecPrt('int1 ',' ',Work(ipint1),ndens2,1)
*       Call RecPrt('int2s ',' ',Work(ipint2s),ndens2,1)
C
       KINT2=ipint2s
       KINT2a=ipint2a
C      If (nt.eq.'T') Then
       If (.false.) Then
       call dcopy_(n2dens,[0.0d0],0,Work(KINT2),1)
       Write(6,*) 'INTEGRAL CHANGE'
       Read(*,*) i,j,k,l
       ij=i+(j-1)*ntash
       kl=k+(l-1)*ntash
       ijkl=itri(ij,kl)
       Work(KINT2+ijkl-1)=1.0D0
       Write(6,*) 'dada',ijkl,Work(KINT2+ijkl-1)
       End If
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
       If (ipint1.eq.0) Then
          Call GetMem('TEMP1INT','ALLO','REAL',ipT1Int,ndens2)
          call dcopy_(ndens2,[0.0d0],0,Work(ipT1Int),1)
c One el int pointer
          KAIN1=ipT1Int
       End If
C
*       Call Getmem('cisigma1','CHECK','REAL',idum,idum)
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
*       Call Getmem('cisigma2','CHECK','REAL',idum,idum)
*
       ist=iispin+1
       square=.false.
*
       If (TIMEDEP) Then
        If (NT.eq.'T')  square=.true.  ! The operator is not sym
        If (.not.page) Then
c ipcidet is here because sigmavec will destroy the first input vector.
         Call GetMem('CIDET_TD','ALLO','REAL',ipCIDET,nDet)
*         If (NT.eq.'S') Then
*             call dcopy_(nCSF(iCSM),Work(ipin(ipCI1)+nConf1),
*     &                    1,Work(ipCIDET),1)
*         Else
*
*       Call Getmem('cisigma3','CHECK','REAL',idum,idum)
*
             call dcopy_(nCSF(iCSM),Work(ipin(ipCI1)),1,Work(ipCIDET),1)
*         End If
C        Do i=1,ntash
C        Do j=1,ntash
C        Do k=1,ntash
C        Do l=1,ntash
C        Write(*,'(I1,I1,I1,I1,F12.5)') i,j,k,l,
C    &    Work(kint2-1+itri(i+(j-1)*ntash,k+(l-1)*ntash))
C        End Do
C        End Do
C        End Do
C        End Do
C         Call RecPrt(' ',' ',Work(ipCIDET),nconf1,1)
*
C
*         If (NT.eq.'T') Then
*           Call SigmaVec(Work(ipCIDET),Work(ipin(ipci2)+nConf1),kic)
*         Else
*           Call RecPrt('ipcidet in cisigma',' ',Work(ipcidet),nconf1,1)
C
           Call SigmaVec(Work(ipCIDET),Work(ipin(ipci2)),kic)
*         End If
C
*          Call Getmem('cisigma','CHECK','REAL',idum,idum)
C
*         call dscal_(nConf1,-1.0d0,Work(ipin(ipci2)),1)
C
*        Call RecPrt('ipci2 in cisigma',' ',Work(ipin(ipci2)),nconf1,1)
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
*           call dcopy_(nCSF(iCSM),Work(ipin(ipCI1)),1,
*     &     Work(ipCIDET),1)
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
                 Work(ipTI2+jilk-1)=Work(ipint2s+ijkl-1)
                End if
               End Do
              End Do
             End Do
           End Do
           Do is=1,nSym
            js=ieor(ieor(icsym-1,issym-1),is-1)+1
            If (nbas(js)*nbas(is).ne.0)
     &      Call DGETMO(Work(ipint1+ipmat(is,js)-1),nbas(is),
     &                nbas(is),nbas(js),Work(ipTI1+ipmat(js,is)-1),
     &                nbas(js))
           End Do
           kain1=ipTI1
           KINT2=ipTI2
         End If  ! End the transpose of integrals.
*
C        Do i=1,ntash
C        Do j=1,ntash
C        Do k=1,ntash
C        Do l=1,ntash
C        Write(*,'(I1,I1,I1,I1,F12.5)') i,j,k,l,
C    &    Work(kint2-1+itri(i+(j-1)*ntash,k+(l-1)*ntash))
C        End Do
C        End Do
C        End Do
C        End Do
C         Call RecPrt('ipcidet in cisig',' ',Work(ipcidet),nconf1,1)
*
C
*         If (NT.eq.'T') Then
*            Call SigmaVec(Work(ipCIDET),Work(ipin(ipci2)),kic)
*         Else
            Call SigmaVec(Work(ipCIDET),Work(ipin(ipci2)+nconf1),kic)
*         End If
C
C        call dscal_(nConf1,-1.0d0,Work(ipin(ipci2)+nconf1),1)
C         Call RecPrt('ipci2 in cisig',' ',Work(ipin(ipci2)),2*nconf1,1)
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
