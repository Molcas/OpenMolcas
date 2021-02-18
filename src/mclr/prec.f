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
* Copyright (C) 1996,1997, Anders Bernhardsson                         *
************************************************************************
      SubRoutine Prec(rpre,idsym)
************************************************************************
*
*  idsym, symmetry of orbital hessian of interest
*  CMtx preconditioner
*
*
* The orbital hessian is dominated of elements that couples
*
* kappa  -> kappa      where i is occupied and p,q is general.
*      ip        iq
*
* we therefore approximate the hessian with thoose diagonal
* terms in the preconditioner
*
*  Anders Bernhardsson 96
*
*     active; active,general is needed for rasscf calculation
*     and is not coded yet (ugly bastard) (970109, AB )
************************************************************************
      use Arrays, only: FAMO, FIMO, F0SQMO
      Implicit Real*8(a-h,o-z)
#include "Pointers.fh"
#include "Input.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "machine.fh"
      Real*8 rpre(*)
*
      Call Prec_internal(rpre)
*
*     This is to allow type punning without an explicit interface
      Contains
      Subroutine Prec_internal(rpre)
      Use Iso_C_Binding
      Real*8, Target :: rpre(*)
      Integer, Pointer :: ipre(:)
      Real*8, Allocatable:: JA(:), KA(:), Scr(:),
     &                      Temp1(:,:), Temp2(:), Temp3(:)
      nmm=0
      nmmm=0
      Do iS=1,nSym
         nMM=Max(nMM,nAsh(is)+nIsh(iS))
         nMMM=Max(nmmM,nBas(is))
      End Do
      n2=nMMM**2
      nmmm=((nmmm-1)/nRec+1)*nRec
      nmm=nmm*nMMM
      nmm=nmm**2

*
      Call mma_allocate(JA,n2,Label='JA')
      Call mma_allocate(KA,n2,Label='KA')
      Call mma_allocate(Scr,n2,Label='Scr')
*
      ip=1
      iAdr=0
      iAdr2=0
      Do iS=1,nSym
         jS=iEOr(is-1,iDSym-1)+1
         nD=nOrb(js)-nIsh(jS)
         ni=nBas(js)**2
         sign=1.0d0
         Call mma_allocate(Temp2,ni,Label='Temp2')
         Call mma_allocate(Temp3,ni,Label='Temp3')
         Call mma_MaxDBLE(nTemp)
         nTemp=Min(nmm,nTemp/2)
         Call mma_allocate(Temp1,nTemp,2,Label='Temp1')
         If (nd.eq.0) Goto 100
         Do iB=1,nIsh(iS)
            call dcopy_(nD**2,[0.0d0],0,Temp3,1)
            ibb=nOrb(is)*(ib-1)+ib-2
*
**          Cholesky code
*
            If (newCho) Then
              Call preci_cho(ib,is,jS,nD,Temp3,
     &                       nOrb(is),nOrb(js),
     &                       FIMO(1+ipCM(is)+ibb),
     &                       FAMO(1+ipCM(is)+ibb),
     &                       F0sqMO(1+ipCM(is)+ibb),
     &                       FIMO(ipCM(js)),
     &                       FAMO(ipCM(js)),
     &                       F0sqMO(ipCM(js)),sign,
     &                       JA,KA,Scr,n2,
     &                       iAdr) ! OK

            Else
            If (iMethod.eq.2) Then
*                                                                      *
************************************************************************
*                                                                      *
*              G
*               iaib
*
               If (nash(js).gt.0)
     &            Call Preciaa(ib,is,js,nd,Temp3,
     &                         nOrb(is),nOrb(js),
     &                         FIMO(1+ipCM(is)+ibb),
     &                         FAMO(1+ipCM(is)+ibb),
     &                         F0sqMO(1+ipCM(is)+ibb),
     &                         FIMO(ipCM(js)),
     &                         FAMO(ipCM(js)),
     &                         F0sqMO(ipCM(js)),sign,
     &                         JA,KA,Scr,n2) ! OK
*                                                                      *
************************************************************************
*                                                                      *
*              G
*               ipia
*
               If ((nOrb(js)-nish(js)-nash(js))*nash(js).gt.0)
     &            Call Preciba(ib,is,js,nd,Temp3,nOrb(js),
     &                         FIMO(ipCM(js)),
     &                         FAMO(ipCM(js)),
     &                         F0sqMO(ipCM(js)),sign,
     &                         JA,KA,Scr,n2) ! OK
*
            End If
*                                                                      *
************************************************************************
*                                                                      *
*           G
*            ipiq
*
            If ((nOrb(js)-nish(js)-nash(js)) .gt.0)
     &         Call Precibb(ib,is,js,nd,Temp3,
     &                      nbas(js),norb(js),
     &                      Temp1(:,1),Temp1(:,2),Temp2,
     &                      FiMo(1+ipCM(is)+ibb),
     &                      FAMO(1+ipcm(is)+ibb),  ! OK
     &                      FiMo(ipCM(js)),
     &                      FAMO(ipcm(js)),sign)  ! OK
            EndIf ! newCho
*                                                                      *
************************************************************************
*                                                                      *
*           Factorize G:
*
*               T
*           G=LL
*
*            write(6,*) 'Preconditioner i =',iB
*            Do i=1,min(nd,10)
*             write(6,'(10F12.8)') (Temp3(1+(j-1)*(2*nd-j+2)/2+i-j),
*     &                             j=1,i)
*            End Do

            Call SQM(Temp3,rpre(ip),nd)

!            write(*,*)" ====== rpre ====== "
!            do i=1,nd*nd
!              write(*,*)i,"rpre",rpre(ip+i-1)
!            end do

            irc=0
            call c_f_pointer(c_loc(rpre(ip+nd**2)),ipre,[nd])
            call dgetrf_(nd,nd,rpre(ip),nd,ipre,irc)
            nullify(ipre)
            If (irc.ne.0) then
               Write(6,*) 'Error in DGETRF called from prec'
               Call Abend
            End If
            ip=ip+nD*(nd+1)
*
         End Do

 100     Continue
*                                                                      *
************************************************************************
*                                                                      *
         Do iB=1,nAsh(iS)
            ibb=nOrb(is)*(nish(is)+ib-1)+nish(is)+ib-2
            If (ib.le.nRs1(iS)+nRs2(is)+nRs3(is)) iR=3
            If (ib.le.nRs1(iS)+nRs2(is)) iR=2
            If (ib.le.nRs1(iS)) iR=1
            If (ir.eq.1) nD=nOrb(js)-nRs1(js)
            If (ir.eq.2) nD=nOrb(js)-nRs2(js)
            If (ir.eq.3) nD=nOrb(js)-nRs3(js)
            If (nd.eq.0) Goto 110
            call dcopy_(nD**2,[0.0d0],0,Temp3,1)
*
**  New Cholesky code
*
            If (newCho) Then
               Call Preca_cho(ib,is,js,nd,ir,Temp3,
     &                        nOrb(is),nOrb(js),
     &                        FIMO(1+ipCM(is)+ibb),
     &                        FAMO(1+ipCM(is)+ibb),
     &                        F0SqMO(1+ipCM(is)+ibb),
     &                        FIMO(ipCM(js)),
     &                        FAMO(ipCM(js)),
     &                        F0SqMO(ipCM(js)),sign,
     &                        JA,KA,Scr,n2,
     &                        iAdr2)
            Else
            If (nish(js).gt.0)
     &         Call Precaii(ib,is,js,nd,ir,Temp3,
     &                      nOrb(is),nOrb(js),
     &                      FIMO(1+ipCM(is)+ibb),
     &                      FAMO(1+ipCM(is)+ibb),
     &                      F0SqMO(1+ipCM(is)+ibb),
     &                      FIMO(ipCM(js)),
     &                      FAMO(ipCM(js)),
     &                      F0SqMO(ipCM(js)),sign,
     &                      JA,KA,Scr,n2) ! OK
*           Call Precaai(ib,nd,ir,rpre(ip))
*           Call Precaaa(ib,nd,ir,rpre(ip))
            If (nish(js)*nOrb(js).gt.0)
     &         Call Precabi(ib,is,js,ir,nd,Temp3,nOrb(js),
     &                      FIMO(ipCM(js)),
     &                      FAMO(ipCM(js)),
     &                      F0SQMO(ipCM(js)),sign,
     &                      JA,KA,Temp1(:,2),n2) !+/-?
*           Call Precaba(ib,nd,ir,rpre(ip))
            If (nOrb(js).gt.0)
     &            Call Precabb_2(ib,is,js,nd,nbas(js),nOrb(js),
     &                           Temp3,
     &                           Temp1(:,1),ntemp,Temp1(:,2),
     &                           Temp2,
     &                           F0SQMO(1+ipCM(is)+ibb),
     &                           FiMo(ipCM(js)),
     &                           FAMO(ipcm(js)) ,
     &                           F0SQMO(ipCM(js)),sign)
*
            EndIf ! newCho

            Call SQM(Temp3,rpre(ip),nD)
            irc=0
            call c_f_pointer(c_loc(rpre(ip+nd**2)),ipre,[nd])
            call dgetrf_(nd,nd,rpre(ip),nd,ipre,irc)
            nullify(ipre)
            If (irc.ne.0) then
               Write(6,*) 'Error in DGETRF called from prec'
               Call Abend
            End If
            ip=ip+nD*(nd+1)
         End Do
110      Continue
         Call mma_deallocate(Temp1)
         Call mma_deallocate(Temp2)
         Call mma_deallocate(Temp3)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(Scr)
      Call mma_deallocate(KA)
      Call mma_deallocate(JA)

*                                                                      *
************************************************************************
*                                                                      *
      Return
      End Subroutine Prec_internal
*
      End
