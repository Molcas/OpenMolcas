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
      SubRoutine Prec_dig(rpre,idsym)
************************************************************************
*                                                                      *
*     idsym, symmetry of orbital hessian of interest                   *
*     CMtx preconditioner                                              *
*                                                                      *
*     The orbital hessian is dominated of elements that couples        *
*                                                                      *
*     kappa  -> kappa      where i is occupied and p,q is general.     *
*          ip        iq                                                *
*                                                                      *
*     we therefore approximate the hessian with thoose diagonal        *
*     terms in the preconditioner                                      *
*                                                                      *
*     Anders Bernhardsson 96                                           *
*                                                                      *
*     active; active,general is needed for rasscf calculation          *
*     and is not coded yet (ugly bastard) (970109, AB )                *
************************************************************************
      use Arrays, only: FAMO, FIMO, F0SQMO
      Implicit Real*8(a-h,o-z)
#include "Pointers.fh"
#include "Input.fh"
#include "stdalloc.fh"
#include "machine.fh"
      Real*8 rpre(*)
      Real*8, Allocatable:: JInt(:), KInt(:), Scr(:)
      Real*8, Allocatable:: Temp1(:,:), Temp2(:), Temp3(:), Temp4(:)
*                                                                      *
************************************************************************
*                                                                      *
*
      Call Prec_dig_internal(rpre)
*
*     This is to allow type punning without an explicit interface
      Contains
      Subroutine Prec_dig_internal(rpre)
      Use Iso_C_Binding
      Real*8, Target :: rpre(*)
      Integer, Pointer :: ipre(:)
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
      Call mma_allocate(JInt,n2,Label='JInt')
      Call mma_allocate(KInt,n2,Label='KInt')
      Call mma_allocate(Scr,n2,Label='Scr')
*
      ip=1
      sign=1.0d0
      Do iS=1,nSym
         jS=iEOr(is-1,iDSym-1)+1
         nD=nBas(js)-nIsh(jS)
         ni=nBas(js)**2
         Call mma_allocate(Temp2,ni,Label='Temp2')
         Call mma_allocate(Temp3,ni,Label='Temp3')
         Call mma_allocate(Temp4,ni,Label='Temp4')
         Temp4(:)=0.0d0
         call mma_MaxDBLE(nTemp)
         nTemp=Min(nmm,nTemp/2)
         Call mma_allocate(Temp1,nTemp,2,Label='Temp1')
         If (nD.eq.0) Goto 100
*
         Do iB=1,nIsh(iS)
            Temp3(1:nD**2)=0.0D0
            ibb=nBas(is)*(ib-1)+ib-2
*
            If (iMethod.eq.iCASSCF) Then
*
*              G
*               iaib
*
               If (nash(js).gt.0)
     &         Call Preciaa(ib,is,js,nd,Temp3,
     &                      nbas(is),nbas(js),
     &                      FIMO(1+ipCM(is)+ibb),
     &                      FAMO(1+ipCM(is)+ibb),
     &                      F0sqMO(1+ipCM(is)+ibb),
     &                      FIMO(ipCM(js)),
     &                      FAMO(ipCM(js)),
     &                      F0sqMO(ipCM(js)),sign,
     &                      JInt,KInt,Scr,n2) ! OK
*
*              G
*               ipia
*
               If ((nbas(js)-nish(js)-nash(js))*nash(js).gt.0)
     &         Call Preciba(ib,is,js,nd,Temp3,nbas(js),
     y                      FIMO(ipCM(js)),
     &                      FAMO(ipCM(js)),
     &                      F0sqMO(ipCM(js)),sign,
     &                      JInt,KInt,Scr,n2) ! OK
            End If
*
*           G
*            ipiq
*
            If ((nbas(js)-nish(js)-nash(js)) .gt.0)
     &      Call Precibb_td(ib,is,js,nd,Temp3,nBas(js),
     &                   Temp1(:,1),Temp1(:,2),Temp2,
     &                   FiMo(1+ipCM(is)+ibb),
     &                   FAMO(1+ipcm(is)+ibb),  ! OK
     &                   FiMo(ipCM(js)),
     &                   FAMO(ipcm(js)),sign)  ! OK
*
*           Factorize G:
*
*               T
*           G=LL
*
            If (.not.timedep) Then
               Call SQM(Temp3,rpre(ip),nd)
#ifdef RS6K
               Call DGEF(rPre(ip),nD,nD,rpre(ip+nD**2))
#else
               irc=0
               call c_f_pointer(c_loc(rpre(ip+nd**2)),ipre,[nd])
               call dgetrf_(nd,nd,rpre(ip),nd,ipre,irc)
               nullify(ipre)
               If (irc.ne.0) Then
                  Write(6,*) 'Error in DGETRF called from prec_dig'
                  Call Abend
               End If
#endif
            Else
               Call SQM(Temp3,Temp4,nD)
               Call SortOutDiagonal2(Temp4,rpre(ip),nd)
            End if
            If (TimeDep) then
               ip=ip+nD
            Else
               ip=ip+nD*(nd+1)
            End If
*
         End Do   ! iB, inactive
 100     Continue
*
         Temp4(1:ni)=0.0d0
         Do iB=1,nAsh(iS)
            ibb=nBas(is)*(nish(is)+ib-1)+nish(is)+ib-2
            If (ib.le.nRs1(iS)+nRs2(is)+nRs3(is)) iR=3
            If (ib.le.nRs1(iS)+nRs2(is)) iR=2
            If (ib.le.nRs1(iS)) iR=1
            If (ir.eq.1) nD=nBas(js)-nRs1(js)
            If (ir.eq.2) nD=nBas(js)-nRs2(js)
            If (ir.eq.3) nD=nBas(js)-nRs3(js)
            If (nd.eq.0) Goto 110
            Temp3(1:nD**2)=0.0d0
            If (nish(js).gt.0)
     &         Call Precaii(ib,is,js,nd,ir,Temp3,
     &                      nbas(is),nbas(js),
     &                      FIMO(1+ipCM(is)+ibb),
     &                      FAMO(1+ipCM(is)+ibb),
     &                      F0SqMO(1+ipCM(is)+ibb),
     &                      FIMO(ipCM(js)),
     &                      FAMO(ipCM(js)),
     &                      F0SqMO(ipCM(js)),sign,
     &                      JInt,KInt,Scr,n2) ! OK
*           Call Precaai(ib,nd,ir,rpre(ip))
*           Call Precaaa(ib,nd,ir,rpre(ip))
            If (nish(js)*nBas(js).gt.0)
     &         Call Precabi(ib,is,js,ir,nd,Temp3,nBas(js),
     &                      FIMO(ipCM(js)),
     &                      FAMO(ipCM(js)),
     &                      F0SQMO(ipCM(js)),sign,
     &                      JInt,KInt,Scr,n2) !+/-?

*           Call Precaba(ib,nd,ir,rpre(ip))
            If (nBas(js).gt.0)
     &         Call Precabb(ib,is,js,nd,nbas(js),Temp3,
     &                      Temp1(:,1),ntemp,Temp1(:,2),
     &                      Temp2,
     &                      F0SQMO(1+ipCM(is)+ibb),
     &                      FiMo(ipCM(js)),
     &                      FAMO(ipcm(js)) ,
     &                      F0SQMO(ipCM(js)),sign)
            If (.not.timedep) then
               Call SQM(Temp3,rpre(ip),nD)
#ifdef RS6K
               Call DGEF(rPre(ip),nD,nd,rpre(ip+nd**2))
#else
               irc=0
               call c_f_pointer(c_loc(rpre(ip+nd**2)),ipre,[nd])
               call dgetrf_(nd,nd,rpre(ip),nd,ipre,irc)
               nullify(ipre)
               If (irc.ne.0) then
                  Write(6,*) 'Error in DGETRF called from prec_dig'
                  Call Abend
               End If
#endif
            Else
*              From Triang mat
               Call SQM(Temp3,Temp4,nD)
               Call SortOutDiagonal2(Temp4,rpre(ip),nd)
            End If
            If (timedep) Then
               ip=ip+nd
            Else
               ip=ip+nD*(nd+1)
            End if
         End Do ! iB
110      Continue
         Call mma_deallocate(Temp1)
         Call mma_deallocate(Temp2)
         Call mma_deallocate(Temp3)
         Call mma_deallocate(Temp4)
      End Do ! End loop over symmetries
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(Scr)
      Call mma_deallocate(KInt)
      Call mma_deallocate(JInt)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End Subroutine Prec_dig_internal
*
      End
*===============================================================
      subroutine SortOutDiagonal(Matrix,diagonal,nb)
*
*    Copy the diagonal elements from Matrix to the vector
*    Diagonal
*
      Implicit Real*8(a-h,o-z)
      Real*8 Matrix(*),diagonal(*)
      call triPrt(' ',' ',matrix,nb)
      ipM=0
      do i=1,nb
       ipM=ipM+i
       Diagonal(i)=Matrix(ipM)
      end do
*      call recPrt(' ',' ',diagonal,nb,1)
      return
      end
*===============================================================
      subroutine SortOutDiagonal2(Matrix,diagonal,nb)
*
*    Copy the diagonal elements from Matrix to the vector
*    Diagonal
*
      Implicit Real*8(a-h,o-z)
      Real*8 Matrix(*),diagonal(*)
      ipM=1
      Do i=1,nb
       Diagonal(i)=Matrix(ipM)
       ipM=ipM+nb+1
      end do
      Return
      End
