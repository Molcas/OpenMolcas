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
      Implicit Real*8(a-h,o-z)
#include "Pointers.fh"
#include "Input.fh"
#include "WrkSpc.fh"
#include "machine.fh"
      Real*8 rpre(*)
*                                                                      *
************************************************************************
*                                                                      *
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
      Call GetMem('JInt','Allo','Real',ipJ,n2)
      Call GetMem('KInt','Allo','Real',ipK,n2)
      Call GetMem('Scr ','Allo','Real',ipS,n2)
*
      ip=1
      sign=1.0d0
      Do iS=1,nSym
         jS=iEOr(is-1,iDSym-1)+1
         nD=nBas(js)-nIsh(jS)
         ni=nBas(js)**2
         Call GetMem('2TEMP2','ALLO','REAL',ipTemp2,ni)
         Call GetMem('3TEMP3','ALLO','REAL',ipTemp3,ni)
         Call GetMem('4TEMP4','ALLO','REAL',ipTemp4,ni)
         call dcopy_(ni,0.0d0,0,Work(ipTemp4),1)
         Call GetMem('1Temp1','MAX','Real',ipT1,nTemp)
         nTemp=Min(nmm,nTemp/2)
         Call GetMem('1Temp1','ALLO','Real',ipTemp1,2*nTemp)
         ipScr=ipTemp1+nTemp
         If (nD.eq.0) Goto 100
*
         Do iB=1,nIsh(iS)
            call dcopy_(nD**2,0.0d0,0,Work(ipTemp3),1)
            ibb=nBas(is)*(ib-1)+ib-2
*
            If (iMethod.eq.iCASSCF) Then
*
*              G
*               iaib
*
               If (nash(js).gt.0)
     &         Call Preciaa(ib,is,js,nd,Work(ipTemp3),
     &                      nbas(is),nbas(js),
     &                      Work(ipFIMO+ipCM(is)+ibb),
     &                      Work(ipFAMO+ipCM(is)+ibb),
     &                      Work(ipF0sqMO+ipCM(is)+ibb),
     &                      Work(ipFIMO+ipCM(js)-1),
     &                      Work(ipFAMO+ipCM(js)-1),
     &                      Work(ipF0sqMO+ipCM(js)-1),sign,
     &                      Work(ipJ),Work(ipK),Work(ipS),n2) ! OK
*
*              G
*               ipia
*
               If ((nbas(js)-nish(js)-nash(js))*nash(js).gt.0)
     &         Call Preciba(ib,is,js,nd,Work(ipTemp3),nbas(js),
     y                      Work(ipFIMO+ipCM(js)-1),
     &                      Work(ipFAMO+ipCM(js)-1),
     &                      Work(ipF0sqMO+ipCM(js)-1),sign,
     &                      Work(ipJ),Work(ipK),Work(ipS),n2) ! OK
            End If
*
*           G
*            ipiq
*
            If ((nbas(js)-nish(js)-nash(js)) .gt.0)
     &      Call Precibb_td(ib,is,js,nd,Work(ipTemp3),nBas(js),
     &                   Work(ipTemp1),Work(ipScr),Work(ipTemp2),
     &                   Work(ipFiMo+ipCM(is)+ibb),
     &                   Work(ipFAMO+ipcm(is)+ibb),  ! OK
     &                   Work(ipFiMo+ipCM(js)-1),
     &                   Work(ipFAMO+ipcm(js)-1),sign)  ! OK
*
*           Factorize G:
*
*               T
*           G=LL
*
            If (.not.timedep) Then
               Call SQM(Work(ipTemp3),rpre(ip),nd)
#ifdef RS6K
               Call DGEF(rPre(ip),nD,nD,rpre(ip+nD**2))
#else
               irc=0
               call dgetrf_(nd,nd,rpre(ip),nd,rpre(ip+nd**2),irc)
               If (irc.ne.0) Then
                  Write(6,*) 'Error in DGETRF called from prec_dig'
                  Call Abend
               End If
#endif
            Else
               Call SQM(Work(ipTemp3),Work(ipTemp4),nD)
               Call SortOutDiagonal2(Work(ipTemp4),rpre(ip),nd)
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
         call dcopy_(ni,0.0d0,0,Work(ipTemp4),1)
         Do iB=1,nAsh(iS)
            ibb=nBas(is)*(nish(is)+ib-1)+nish(is)+ib-2
            If (ib.le.nRs1(iS)+nRs2(is)+nRs3(is)) iR=3
            If (ib.le.nRs1(iS)+nRs2(is)) iR=2
            If (ib.le.nRs1(iS)) iR=1
            If (ir.eq.1) nD=nBas(js)-nRs1(js)
            If (ir.eq.2) nD=nBas(js)-nRs2(js)
            If (ir.eq.3) nD=nBas(js)-nRs3(js)
            If (nd.eq.0) Goto 110
            call dcopy_(nD**2,0.0d0,0,Work(ipTemp3),1)
            If (nish(js).gt.0)
     &         Call Precaii(ib,is,js,nd,ir,Work(ipTemp3),
     &                      nbas(is),nbas(js),
     &                      Work(ipFIMO+ipCM(is)+ibb),
     &                      Work(ipFAMO+ipCM(is)+ibb),
     &                      Work(ipF0SqMO+ipCM(is)+ibb),
     &                      Work(ipFIMO+ipCM(js)-1),
     &                      Work(ipFAMO+ipCM(js)-1),
     &                      Work(ipF0SqMO+ipCM(js)-1),sign,
     &                      Work(ipJ),Work(ipK),Work(ipS),n2) ! OK
*           Call Precaai(ib,nd,ir,rpre(ip))
*           Call Precaaa(ib,nd,ir,rpre(ip))
            If (nish(js)*nBas(js).gt.0)
     &         Call Precabi(ib,is,js,ir,nd,Work(ipTemp3),nBas(js),
     &                      Work(ipFIMO+ipCM(js)-1),
     &                      Work(ipFAMO+ipCM(js)-1),
     &                      Work(ipF0SQMO+ipCM(js)-1),sign,
     &                      Work(ipJ),Work(ipK),Work(ipS),n2) !+/-?

*           Call Precaba(ib,nd,ir,rpre(ip))
            If (nBas(js).gt.0)
     &         Call Precabb(ib,is,js,nd,nbas(js),Work(ipTemp3),
     &                      Work(ipTemp1),ntemp,Work(ipScr),
     &                      Work(ipTemp2),
     &                      Work(ipF0SQMO+ipCM(is)+ibb),
     &                      Work(ipFiMo+ipCM(js)-1),
     &                      Work(ipFAMO+ipcm(js)-1) ,
     &                      Work(ipF0SQMO+ipCM(js)-1),sign)
            If (.not.timedep) then
               Call SQM(Work(ipTemp3),rpre(ip),nD)
#ifdef RS6K
               Call DGEF(rPre(ip),nD,nd,rpre(ip+nd**2))
#else
               irc=0
               call dgetrf_(nd,nd,rpre(ip),nd,rpre(ip+nd**2),irc)
               If (irc.ne.0) then
                  Write(6,*) 'Error in DGETRF called from prec_dig'
                  Call Abend
               End If
#endif
            Else
*              From Triang mat
               Call SQM(Work(ipTemp3),Work(ipTemp4),nD)
               Call SortOutDiagonal2(Work(ipTemp4),rpre(ip),nd)
            End If
            If (timedep) Then
               ip=ip+nd
            Else
               ip=ip+nD*(nd+1)
            End if
         End Do ! iB
110      Continue
         Call GetMem('1TEMP1','FREE','REAL',ipTemp1,nTemp)
         Call GetMem('2TEMP2','FREE','REAL',ipTemp2,nBas(js)**2)
         Call GetMem('3TEMP3','FREE','REAL',ipTemp3,nBas(js)**2)
         Call GetMem('4TEMP4','FREE','REAL',ipTemp4,ni)
      End Do ! End loop over symmetries
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('Scr ','Free','Real',ipS,n2)
      Call GetMem('KInt','Free','Real',ipK,n2)
      Call GetMem('JInt','Free','Real',ipJ,n2)
*                                                                      *
************************************************************************
*                                                                      *
      Return
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
