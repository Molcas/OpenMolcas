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
      SubRoutine Prec_td(pre2,DigPrec,isym)
*
*     pre2      Preconditioner from Prec
*     DigPrec Output - Diagonal of prec2
*     isym      Symmetry of PT
*
      Implicit Real*8 (a-h,o-z)
#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
      Real*8 nonzero
      Real*8 DigPrec(*),pre2(*)
      Logical jump
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*
*---------------------------------------------
* Construct one el density in MO ipdens
*---------------------------------------------
*
      nBasTot = 0
      Do iS=1,nSym
         nBasTot = nBasTot + nBas(iS)*nBas(iS)
      End Do
      Call GetMem('densmat','Allo','Real',ipDens,nBasTot)
      call dcopy_(nBasTot,0.0d0,0,Work(ipDens),1)
*
      ip3 = 0
      Do iS=1,nSym
          inc = nBas(iS)+1
          call dcopy_(nIsh(iS),2.0d0,0,Work(ipDens+ip3),inc)
          ip3 = ip3 + nBas(iS)*nBas(iS)
      End Do
*
* For a CASSCF wavefunc. From Anders subrut r2elint
* Add the active active dens
*
      Do iS=1,nSym
          Do iB=1,nAsh(iS)
              Do jB=1,nAsh(iS)
                   ip=ipCM(iS)+ib+nIsh(is)+(jB+nIsh(is)-1)*nBas(is)-1
                   iA=nA(is)+ib
                   jA=nA(is)+jb
                   ip2=itri(iA,jA)
                   Work(ipDens+ip-1)=Work(ipG1t+ip2-1)
              End Do
           End Do
      End Do
*
C
*      Call RECPRT('dens',' ',Work(ipDens),nBasTot,1)
*      Write(*,*)'Diagnonal elements in D'
*      Do iS=1,nSym
*         Do k=0,nBas(iS)-1
*            Write(*,*) Work(ipDens + ipCM(iS)-1 + k*(nBas(iS)+1))
*         End Do
*      End Do
*      Stop
C
*
*-------------------------------------------------------------------
* Construct the diagonal approximation to the orbital prec, ipPrecTd
*-------------------------------------------------------------------
*
      Call GetMem('prectd','Allo','Real',ipPreTd,nDens2)
      call dcopy_(nDens2,0.0d0,0,Work(ipPreTd),1)
      ip1 = 0
      ip2 = 1
      ipsave = 0
      Do iS=1,nSym
         jS=iEOr(iS-1,iSym-1)+1
         nD = nBas(jS) - nIsh(jS)
         Do k=1,nIsh(iS)
            ip1 = ip1 + nIsh(jS)
            Do l=1, nD
               Work(ipPreTd + ip1) = Pre2(ip2)
               ip1 = ip1 + 1
               ip2 = ip2 + 1
            End Do
         End Do
         nD = nBas(jS) - nAsh(jS)
         Do k=1,nAsh(iS)
            jump = .true.
            Do l=1, nD
               If (l.gt.nIsh(jS).and.jump) Then
                   ip1 = ip1 + nAsh(jS)
                   jump = .false.
               End If
               Work(ipPreTd + ip1) = Pre2(ip2)
               ip1 = ip1 + 1
               ip2 = ip2 + 1
            End Do
            If ((nBas(jS) - nAsh(jS) - nIsh(jS)).eq.0) Then
               ip1 = ip1 + nAsh(jS)
            End If
         End Do
*         Call RECPRT('PreTd1',' ',Work(ipPreTd + ipsave),
*     &              nBas(jS),nBas(iS))
         ip1 = ip1 + (nBas(iS)-nIsh(iS)-nAsh(iS))*nBas(jS)
         ipsave = ip1
      End Do
*      Call RECPRT('PreTd',' ',Work(ipPreTd),nDens2,1)
*
*-----------------------------------------------------------
* Symmetrize ipPreTd
*-----------------------------------------------------------
      Call GetMem('temptd','Allo','Real',ipTempTd,nDens2)
      call dcopy_(nDens2,0.0d0,0,Work(ipTempTd),1)
*
      Do iS=1,nSym
         jS=iEOr(iS-1,iSym-1)+1
         call dcopy_(nDens2,0.0d0,0,Work(ipTempTd),1)
         Call Trans(Work(ipPreTd+ipMat(jS,iS)-1),nBas(iS),
     &               nBas(jS),Work(ipTempTd))
         nD = nBas(iS)*nBas(jS)
         Do i=0, nD-1
            nonzero = Work(ipPreTd + ipMat(iS,jS) -1 +i)
            If (nonzero.ne.0.0d0) Then
                Work(ipTempTd +i) =
     &              Work(ipPreTd+ipMat(iS,jS)-1+i)
            End If
         End Do
         Call Trans(Work(ipTempTd),nBas(jS),nBas(iS),
     &              Work(ipPreTd+ipMat(jS,iS)-1))
C
      End Do
*
*------------------------------------------------------------------
* Add the density part ipPreTd_at = ipPreTd_at - omega(D_aa - D_tt)
*------------------------------------------------------------------
      i = 0
      Do iS=1,nSym
         jS=iEOr(iS-1,iSym-1)+1
         nD = nBas(iS)*nBas(jS)
         j = 0
         l = 0
         Do k=0, nD-1
            If (l.eq.nBas(jS)) l = 0
            If (k.eq.(j+1)*nBas(jS)) j = j +1
C
                Work(ipPreTd +i) = Work(ipPreTd +i) +
     &     2.0d0*omega*(Work(ipDens + ipCM(iS)-1 + j*(nBas(iS)+1)) +
     &     Work(ipDens + ipCM(jS)-1 + l*(nBas(jS)+1)))
C
            l = l + 1
            i=i+1
         End Do
C
      End Do
*
*-----------------------------------------------------------------------------
* Symmetry transpose ipPreTd - To get the same order fo sym as in b_x and
* as required by compress.
*-----------------------------------------------------------------------------
*
      call dcopy_(nDens2,0.0d0,0,Work(ipTempTd),1)
*
      Do iS=1,nSym
         jS=iEOr(iS-1,iSym-1)+1
         nD = nBas(iS)*nBas(jS)
         Do k=0, nD-1
            Work(ipTempTd + ipmat(iS,jS)-1 +k)=
     &            Work(ipPreTd +ipmat(jS,iS)-1 +k)
         End Do
      End Do
*
C
      Call Compress(Work(ipTempTd),DigPrec,isym)
C
      Call GetMem('densmat','Free','Real',ipDens,nBasTot)
      Call GetMem('prectd','Free','Real',ipPreTd,nDens2)
      Call GetMem('temptd','Free','Real',ipTempTd,nDens2)
*
      Return
      End
