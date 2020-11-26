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
      use Arrays, only: G1t
*
*     pre2      Preconditioner from Prec
*     DigPrec Output - Diagonal of prec2
*     isym      Symmetry of PT
*
      Implicit Real*8 (a-h,o-z)
#include "Input.fh"
#include "Pointers.fh"
#include "stdalloc.fh"
      Real*8 nonzero
      Real*8 DigPrec(*),pre2(*)
      Logical jump
      Real*8, Allocatable:: Dens(:), PreTd(:), TempTd(:)
*                                                                      *
************************************************************************
*                                                                      *
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
*
*---------------------------------------------
* Construct one el density in MO Dens
*---------------------------------------------
*
      nBasTot = 0
      Do iS=1,nSym
         nBasTot = nBasTot + nBas(iS)*nBas(iS)
      End Do
      Call mma_allocate(Dens,nBasTot,Label='Dens')
      Dens(:)=0.0D0
*
      ip3 = 1
      Do iS=1,nSym
          inc = nBas(iS)+1
          call dcopy_(nIsh(iS),[2.0d0],0,Dens(ip3),inc)
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
                   Dens(ip)=G1t(ip2)
              End Do
           End Do
      End Do
*
C
*      Call RECPRT('Dens',' ',Dens,nBasTot,1)
*      Write(*,*)'Diagnonal elements in D'
*      Do iS=1,nSym
*         Do k=0,nBas(iS)-1
*            Write(*,*) Dens(ipCM(iS) + k*(nBas(iS)+1))
*         End Do
*      End Do
*      Stop
C
*
*-------------------------------------------------------------------
* Construct the diagonal approximation to the orbital prec, PreTd
*-------------------------------------------------------------------
*
      Call mma_allocate(PreTd,nDens2,Label='PreTd')
      PreTd(:)=0.0d0
      ip1 = 1
      ip2 = 1
      ipsave = 0
      Do iS=1,nSym
         jS=iEOr(iS-1,iSym-1)+1
         nD = nBas(jS) - nIsh(jS)
         Do k=1,nIsh(iS)
            ip1 = ip1 + nIsh(jS)
            Do l=1, nD
               PreTd(ip1) = Pre2(ip2)
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
               PreTd(ip1) = Pre2(ip2)
               ip1 = ip1 + 1
               ip2 = ip2 + 1
            End Do
            If ((nBas(jS) - nAsh(jS) - nIsh(jS)).eq.0) Then
               ip1 = ip1 + nAsh(jS)
            End If
         End Do
*         Call RECPRT('PreTd',' ',PreTd(1 + ipsave),
*     &              nBas(jS),nBas(iS))
         ip1 = ip1 + (nBas(iS)-nIsh(iS)-nAsh(iS))*nBas(jS)
         ipsave = ip1
      End Do
*      Call RECPRT('PreTd',' ',PreTd,nDens2,1)
*
*-----------------------------------------------------------
* Symmetrize PreTd
*-----------------------------------------------------------
      Call mma_allocate(TempTd,nDens2,Label='TempTd')
*
      Do iS=1,nSym
         jS=iEOr(iS-1,iSym-1)+1
         TempTd(:)=0.0d0
         Call Trans(PreTd(ipMat(jS,iS)),nBas(iS),
     &               nBas(jS),TempTd)
         nD = nBas(iS)*nBas(jS)
         Do i=0, nD-1
            nonzero = PreTd(ipMat(iS,jS) +i)
            If (nonzero.ne.0.0d0) Then
                TempTd(1+i) = PreTd(ipMat(iS,jS)+i)
            End If
         End Do
         Call Trans(TempTd,nBas(jS),nBas(iS),
     &              PreTd(ipMat(jS,iS)))
C
      End Do
*
*------------------------------------------------------------------
* Add the density part PreTd_at = PreTd_at - omega(D_aa - D_tt)
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
            i=i+1
            PreTd(i) = PreTd(i) + 2.0D0*Omega*
     &               (  Dens(ipCM(iS) + j*(nBas(iS)+1)) +
     &                  Dens(ipCM(jS) + l*(nBas(jS)+1)) )
C
            l = l + 1
         End Do
C
      End Do
*
*-----------------------------------------------------------------------------
* Symmetry transpose PreTd - To get the same order fo sym as in b_x and
* as required by compress.
*-----------------------------------------------------------------------------
*
      TempTd(:)=0.0d0
*
      Do iS=1,nSym
         jS=iEOr(iS-1,iSym-1)+1
         nD = nBas(iS)*nBas(jS)
         Do k=0, nD-1
            TempTd(ipmat(iS,jS)+k)=PreTd(ipmat(jS,iS)+k)
         End Do
      End Do
*
      Call Compress(TempTd,DigPrec,isym)

      Call mma_deallocate(TempTd)
      Call mma_deallocate(PreTd)
      Call mma_deallocate(Dens)
*
      Return
      End
