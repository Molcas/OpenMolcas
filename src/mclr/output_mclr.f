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
* Copyright (C) 1996, Anders Bernhardsson                              *
************************************************************************
       SubRoutine OutPut_Mclr(iKapDisp,isigdisp,iCiDisp,
     &                        iCiSigDisp,iRHSDisp,iRHSCIDisp,
     &                        converged)
********************************************************************
*                                                                  *
* Contracts the response coefficient to the hessian                *
*                                                                  *
* Input                                                            *
*       iKapDisp : Disk locations of solutions to respons equation *
*       iSigDisp : Disk locations of RHS                           *
*       iCIDisp  : Disk locations of CI Soulutions to response     *
*       iCISigDisp : Disk locations of RHS                         *
*       nHess    : Length of hessian                               *
*                                                                  *
* Output to disk                                                   *
*                                                                  *
*       RespHess : Hessian etc                                     *
*       Hess     : Hessian etc                                     *
*                                                                  *
* Author: Anders Bernhardsson, 1996                                *
*         Theoretical Chemistry, University of Lund                *
********************************************************************
       Implicit Real*8 (a-h,o-z)
#include "detdim.fh"

#include "Input.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "disp_mclr.fh"
#include "cicisp_mclr.fh"
#include "WrkSpc.fh"
       Character*8 Label
#ifdef _DEBUG_
       Character*20 Label2
#endif
       Integer Pstate_sym,ldisp2(8),ielec(3)
       Integer iKapDisp(nDisp),isigdisp(nDisp),
     &         iCiDisp(nDisp),iCiSigDisp(nDisp),
     &         iRHSDisp(nDisp),iRHSCiDisp(nDisp)
       Logical elec,converged(8),CI
       Real*8 Pola(6)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
*
       Call QEnter('Output')
#ifdef _DEBUG_
       debug=.True.
#else
       debug=.false.
#endif
       nHss=0
       Do iS=1,nSym
        nHss=nHss+lDisp(is)*(lDisp(is)+1)/2
       End Do
       nhess=nDIsp*(nDisp+1)/2
       Call GetMem('RESPH','ALLO','REAL',ipRHss,nHss)
       call dcopy_(nHss,[0.0d0],0,Work(ipRHss),1)
*
*-------------------------------------------------------------------*
*
* Ok construct hessian
*
*-------------------------------------------------------------------*
*
       mSym=0
       kSym=0
       idum=1
       idisp=0
       Do 100 iSym=1,nSym
*
*         Calculate length of the density, Fock and Kappa matrix etc
*         notice that this matrixes not necessary are symmetric.
*         Store pointers.
*
*         Input:
*         iSym : Symmetry of perturbation
*
*         Output: Commonblocks (Pointers.fh)
*
          Call Setup_MCLR(iSym)
          PState_SYM=iEor(State_Sym-1,iSym-1)+1
          nconfM=Max(ncsf(PState_Sym),nint(xispsm(Pstate_Sym,1)))
          nconf1=ncsf(PState_Sym)
          CI=.false.
          If (iMethod.eq.2.and.nconf1.gt.0) CI=.true.
          If (CI.and.nconf1.eq.1.and.isym.eq.1) CI=.false.
*
*         Allocate areas for scratch and state variables
*
          Call GetMem('kappa1','Allo','Real',ipKap1,nDensC)
          Call GetMem('kappa2','Allo','Real',ipKap2,nDensC)
          Call GetMem('skappa','Allo','Real',ipsKap,nDensC)
          Call Getmem('rkappa1','ALLO','Real',iprkap1,nDensC)
          Call Getmem('rkappa2','ALLO','Real',iprkap2,nDensC)
*
*
          If (CI) Then
             ipcip1=ipget(nconfM)
             ipcip2=ipget(nconfM)
             ipsp=ipget(nconfM)
             iprp2=ipget(nconfM)
             iprp1=ipget(nconfM)
          End If
*
*                                    [2]
*         Calculate the diagonal of E    and store in core/disc
*
*
*
        Do 110 jDisp=1,lDisp(iSym)
          iDisp=iDisp+1
          jspin=0
          If (iAnd(nTPert(idisp),1).eq.1) jSpin=1
          If (jspin.eq.0) Then
             nconf1=ncsf(Pstate_sym)
          Else
             nconf1=nint(xispsm(Pstate_Sym,1))
          End If
*
*---------------------------------------------------------------*
*                                                               *
*    Read response from disk                                    *
*                                                               *
*---------------------------------------------------------------*
*
          iDisk=iKapDisp(iDisp)
*
*-------- If disk address =/= -1 arrays on file.
*
          If (iDisk.ne.-1) Then
              Len=nDensC
              Call dDaFile(LuTemp,2,Work(ipKap1),Len,iDisk)
              iDisk=iSigDisp(iDisp)
              Call dDaFile(LuTemp,2,Work(ipSKap),Len,iDisk)
              iDisk=iRHSDisp(iDisp)
              Call dDaFile(LuTemp,2,Work(iprKap1),Len,iDisk)
              Do i=0,ndensC-1
                 Work(ipSKap+i)=-Work(ipSKap+i)-Work(iprKap1+i)
              End Do
C
*             Call Recprt('ORB-RHS',' ',Work(iprKap1),nDensC,1)
*             Write(*,*)'ddot orb-resp',
*     &            ddot_(ndensC,Work(ipKap1),1,Work(ipKap1),1)
*             Write(*,*)'ddot orb-sigma',
*     &            ddot_(ndensC,Work(ipSKap),1,Work(ipSKap),1)
*             Write(*,*)'ddot orb-rhs',
*     &            ddot_(ndensC,Work(iprKap1),1,Work(iprKap1),1)
C
*
             Call GADSum(Work(ipKap1),Len)
             Call GADSum(Work(ipSKap),Len)
             Call GADSum(Work(iprKap1),Len)

             If (CI) Then
                ilen=nconf1
                idis=iCIDisp(iDisp)
                Call dDaFile(LuTemp,2,Work(ipin(ipCIp1)),iLen,iDis)
                idis=iCISigDisp(idisp)
                Call dDaFile(LuTemp,2,Work(ipin(ipSp)),iLen,iDis)
                idis=iRHSCIDisp(idisp)
                Call dDaFile(LuTemp,2,Work(ipin(iprp1)),iLen,iDis)
                ii=ipin(ipSp)
                jj=ipin(iprp1)
                Do i=0,nConf1-1
                   Work(ii+i)= -Work(ii+i)-Work(jj+i)
                End Do
C
*               Write(*,*)'ddot ci-resp',
*     &               ddot_(nConf1,Work(ipin(ipcip1)),1,
*     &                      Work(ipin(ipcip1)),1)
*               Write(*,*)'ddot ci-sigma',
*     &               ddot_(nConf1,Work(ipin(ipSp)),1,
*     &                      Work(ipin(ipSp)),1)
*               Write(*,*)'ddot ci-rhs',
*     &               ddot_(nConf1,Work(ipin(iprp1)),1,
*     &                      Work(ipin(iprp1)),1)
C
                Call GADSum(Work(ipin(ipCIp1)),iLen)
                Call GADSum(Work(ipin(ipSp  )),iLen)
                Call GADSum(Work(ipin(ipRp1 )),iLen)
             End If
*
          Else
*
                            Len = nDensC
             Call FZero(Work(ipKap1),Len)
             Call GADSum(Work(ipKap1),Len)
             Call FZero(Work(ipSKap),Len)
             Call GADSum(Work(ipSKap),Len)
             Call FZero(Work(iprKap1),Len)
             Call GADSum(Work(iprKap1),Len)
*
             If (CI) Then
*
                ilen=nconf1
                Call FZero(Work(ipin(ipCIp1)),iLen)
                Call GADSum(Work(ipin(ipCIp1)),iLen)
                Call FZero(Work(ipin(ipSp  )),iLen)
                Call GADSum(Work(ipin(ipSp  )),iLen)
                Call FZero(Work(ipin(ipRp1 )),iLen)
                Call GADSum(Work(ipin(ipRp1 )),iLen)
*
             End If
*
          End If
*
************************************************************************
*
           Do 120 kDisp=1, jdisp
*
*                    (x) (2) (y)   (x) (y)    (y) (x)
*            E    = k   E   k   + F   k   +  F   k
*             Resp
*
               kspin=0
               If (iAnd(nTPert(kdisp+ksym),1).eq.1) kSpin=1
               If (kspin.eq.0) Then
                 nconf1=ncsf(PState_Sym)
               Else
                 nConf1=nint(xispsm(Pstate_Sym,1))
               End If
               If (.not.lCalc(kDisp+ksym)) Goto 120
*
C
*                  Write(*,*)'kDisp+kSym',kDisp+kSym
*                  Write(*,*)'iKapDisp(kdisp+ksym)',iKapDisp(kdisp+ksym)
C
               iDisk=iKapDisp(kDisp+kSym)
               If (iDisk.ne.-1) Then
                  Len=nDensC
                  Call dDaFile(LuTemp,2,Work(ipKap2),Len,iDisk)
                  iDisk=iRHSDisp(kDisp+kSym)
                  Call dDaFile(LuTemp,2,Work(iprKap2),Len,iDisk)

                  Call GASync()
                  Call GADSum(Work(ipKap2 ),Len)
                  Call GADSum(Work(iprKap2),Len)

                  If (CI) Then
                     ilen=nconf1
                     idis=iCIDisp(kDisp+ksym)
                     Call dDaFile(LuTemp,2,Work(ipin(ipCIp2)),iLen,iDis)
                     idis=iRHSCIDisp(kdisp+ksym)
                     Call dDaFile(LuTemp,2,Work(ipin(iprp2)),iLen,iDis)
                     rTempc1=DDot_(nConf1,Work(ipin(ipCIp2)),1,
     &                                   Work(ipin(ipsp)),1)

                          Call GASync()
                     Call GADSum(Work(ipin(ipCIp2)),iLen)
                     Call GADSum(Work(ipin(iprp2 )),iLen)

                  Else
                     rtempc1=0.0d0
                  End If
*
               Else
*
                  Call GASync()
                  Len=nDensC
                  Call FZero(Work(ipKap2 ),Len)
                  Call GADSum(Work(ipKap2 ),Len)
                  Call FZero(Work(iprKap2),Len)
                  Call GADSum(Work(iprKap2),Len)
                  If (CI) Then
                     ilen=nconf1
                          Call GASync()
                     Call FZero(Work(ipin(ipCIp2)),iLen)
                     Call GADSum(Work(ipin(ipCIp2)),iLen)
                     Call FZero(Work(ipin(iprp2 )),iLen)
                     Call GADSum(Work(ipin(iprp2 )),iLen)
                     rTempc1=DDot_(nConf1,Work(ipin(ipCIp2)),1,
     &                                   Work(ipin(ipsp)),1)
                  Else
                     rtempc1=0.0d0
                  End If
*
               End If
*
               rTempk1=DDot_(nDensC,Work(ipKap2),1,Work(ipSKap),1)
*
               Fact=1.0d0
               If (kdisp.eq.jdisp) Fact=2.0d0
               rTempk2=Fact*
     &                DDot_(nDensC,Work(ipKap1),1,Work(iprKap2),1)
               If (kdisp.ne.jdisp) Then
               rtempk3=1.0d0*DDot_(nDensC,Work(iprKap1),1,
     &                  Work(ipKap2),1)
               Else
                rTempk3=0.0d0
               End if
               If (CI) Then
                 Fact=1.0d0
                 If (kdisp.eq.jdisp) Fact=2.0d0
                  rTempc2=Fact*
     &                DDot_(nConf1,Work(ipin(ipCip1)),1,
     &                            Work(ipin(iprp2)),1)
                 If (kdisp.ne.jdisp) Then
                  rtempc3=1.0d0*
     &                DDot_(nConf1,Work(ipin(iprp1)),1,
     &                            Work(ipin(ipCIp2)),1)
                 Else
                  rTempc3=0.0d0
                 End if
                Else
                  rtempc2=0.0d0
                  rtempc3=0.0d0
               End if
C
*              Write(*,*) kdisp,jdisp
*              Write(*,*) rTempk1,rtempk2,rtempk3
*              Write(*,*) rtempc1,rtempc2,rtempc3
C
               Maxi=Max(kDisp,jDisp)
               Mini=Min(kDisp,jDisp)
               index=mSym+Maxi*(Maxi-1)/2+Mini
*
               Work(ipRhss-1+Index)=Work(ipRhss-1+Index)+
     &             rTempk1+rtempk2+rtempk3+
     &             rtempc1+rtempc2+rtempc3
*
 120       Continue

**********************************************************************
*
 110     Continue
       kSym=kSym+lDisp(iSym)
       mSym=mSym+lDisp(iSym)*(lDisp(iSym)+1)/2
*
*    Free areas for scratch and state variables
*
          Call Getmem('rkappa2','FREE','Real',iprkap2,nDensC)
          Call Getmem('rkappa1','FREE','Real',iprkap1,nDensC)
          Call GetMem('skappa','FREE','Real',ipsKap,nDensC)
          Call GetMem('kappa2','FREE','Real',ipKap2,nDensC)
          Call GetMem('kappa1','FREE','Real',ipKap1,nDensC)
          If (CI)  irc=ipclose(ipcip1)
 100  Continue
      Call GetMem('HESS','ALLO','REAL',ipHess,nHss)
      Call GetMem('HESS','ALLO','REAL',ipHess2,nHss)
      Call GetMem('Temp','ALLO','REAL',ipTemp,nHss)
      Call FZero(Work(ipTemp),nHss)
      Call GetMem('Temp','ALLO','REAL',ipELEC,3*ndisp)
C     Call GetMem('Temp','ALLO','REAL',ipELEC2,3*ndisp)
      Call GetMem('Temp','ALLO','REAL',ipEG  ,3*ndisp)
      Call GetMem('Temp','ALLO','REAL',ipELOUT,3*ndisp)
      irc=ipclose(-1)
*
*-------------------------------------------------------------------*
*
*     OK now when we have out Hessian, what should  we do with it!
*
*-------------------------------------------------------------------*
*
*
*     If a basis set is dependent on perturbation add terms
*     constructed in mckinley.
*
      call dcopy_(6,[0.0d0],0,pola,1)
      idum=1
      iopt=128
      irc=3*ndisp
      Label='DOTELGR'
      Call drdMCk(irc,iopt,LaBeL,idum,Work(ipEG),idum)
      elec=.true.
      if (irc.ne.0) elec=.false.
                Call GADsum(Work(iphss), nhss)
      call dcopy_(nhss,Work(iphss),1,Work(iphess2),1)
#ifdef _DEBUG_
      If (debug) Then
         ip=ipHess2
         Do iSym=1,nSym
           Write(label2,'(A,I2)') 'CHessian symmetry',iSym
            If (lDisp(iSym).ne.0)
     &         Call TriPrt(label2,' ',Work(ip),lDisp(iSym))
            ip=ip+ldisp(isym)*(1+ldisp(isym))/2
         End Do
      End If
*
      If (debug) Then
       Call MMSORT2(Work(ipHESS2),Work(ipELEC),pola,ielec)
       Call Recprt('CONN',' ',Work(ipElec),3*nDisp,1)
      End If
*
C
c       Write(*,*)'I am here 1'
       Call Recprt('ipRhss','(5G20.10) ',Work(ipRHss),nhss,1)
       Call Recprt('iphss','(5G20.10) ',Work(ipHss),nhss,1)
#endif
C
      Call DaXpY_(mSym,1.0d0,Work(ipRHss),1,Work(ipHess2),1)
*
#ifdef _DEBUG_
      If (debug) Then
       Call MMSORT2(Work(ipRHSS),Work(ipELEC),pola,ielec)
       Call Recprt('RESP',' ',Work(ipElec),3*nDisp,1)
      End If
#endif
*
      Call MMSORT2(Work(ipHESS2),Work(ipELEC),pola,ielec)
*
#ifdef _DEBUG_
      If (debug) Then
       Call Recprt('R+C',' ',Work(ipElec),3*nDisp,1)
       ip=ipHess2
       Do iSym=1,nSym
        Write(label2,'(A,I2)') 'Hessian symmetry',iSym
        If (lDisp(iSym).ne.0)
     &   Call TriPrt(label2,' ',Work(ip),lDisp(iSym))
        ip=ip+ldisp(isym)*(1+ldisp(isym))/2
       End Do
      End If
#endif
*
      Call mmSort(Work(ipHess2),Work(ipHess),ldisp2)
*
#ifdef _DEBUG_
C
cvv       Write(*,*)'I am here'
c       Call Recprt('iphess2',' ',Work(ipHess2),nhss,1)
c       Call HssPrt_MCLR(iwork(ipdegdisp),Work(ipHess2),ldisp2)
#endif
      If (McKinley) Then
*
         iRC=-1
         iOpt=0
         Label='StatHess'
         Call dRdMck(iRC,iOpt,Label,idum,Work(ipTemp),idum)
         If (iRC.ne.0) Then
            Write (6,*) 'OutPut: Error reading MCKINT'
            Write (6,'(A,A)') 'Label=',Label
            Call QTrace
            Call Abend()
         End If
*
#ifdef _DEBUG_
         If (debug) Then
            ip=ipTemp
            Do iSym=1,nSym
               Write(label2,'(a,i2)') 'SHessian symmetry',iSym
               If (lDisp2(iSym).ne.0)
     &            Call TriPrt(label2,' ',Work(ip),lDisp2(iSym))
               ip=ip+ldisp2(isym)*(1+ldisp2(isym))/2
            End Do
         End If
#endif
         Call DaXpY_(mSym,1.0d0,Work(ipTemp),1,Work(ipHess),1)
      End If
#ifdef _DEBUG_
      If (debug) Then
        ip=ipHess
        Do iSym=1,nSym
          Write(label2,'(a,i2)') 'Hessian symmetry',iSym
          If (lDisp2(iSym).ne.0)
     &    Call TriPrt(label2,' ',Work(ip),lDisp2(iSym))
          ip=ip+ldisp2(isym)*(1+ldisp2(isym))/2
        End Do
      End If
#endif
*
*
      If (McKinley) Then
         iRC=-1
         iOpt=0
         Label='Hess    '
         Call dWrMck(iRC,iOpt,Label,iDum,Work(ipHess),iDum)
         If (iRC.ne.0) Then
            Write (6,*) 'OutPut: Error writing to MCKINT'
            Write (6,'(A,A)') 'Label=',Label
            Call QTrace
            Call Abend()
         End If
         Call Put_iScalar('No of Internal coordinates',ldisp2(1))
         Call Put_AnalHess(Work(ipHess),ldisp2(1)*(ldisp2(1)+1)/2)
      End If
*
      If (.true.) Then
       iRC=-1
       iOpt=0
       Call GetMem('NRDISP','ALLO','INTE',ipnrdisp,ndisp)
       Call RdMck(irc,iopt,'NRCTDISP',idum,iWork(ipnrDisp),idum)
       iRC=-1
       iOpt=0
       Call GetMem('DEGDISP','ALLO','INTE',ipdegdisp,ndisp)
       Label='DEGDISP'
       Call RdMck(irc,iopt,Label,idum,iWork(ipDegDisp),idum)
       If (iRC.ne.0) Then
          Write (6,*) 'OutPut: Error reading RELAX'
          Write (6,'(A,A)') 'Label=',Label
          Call QTrace
          Call Abend()
       End If
C
*       If (debug)
*        Call HssPrt_MCLR(iwork(ipdegdisp),Work(ipHess),ldisp2)
*       Call Recprt('iphess',' ',Work(ipHess),nhss,1)
C
       call daxpy_(3*ndisp,-1.0d0,Work(ipEG),1,Work(ipELEC),1)
#ifdef _DEBUG_
       If (debug.and.elec)
     &  Call Recprt('ELEC-ST',' ',Work(ipEG),3*nDisp,1)
       If (debug.and.elec)
     &  Call Recprt('ELEC-TOT',' ',Work(ipElec),3*nDisp,1)
#endif
*
       Lu_10=10
       Lu_10=IsFreeUnit(Lu_10)
       call molcas_open(lu_10,'UNSYM')
c       Open(unit=Lu_10, file='UNSYM')
*
       If (Mckinley) Then
           Call FreqAnal(iwork(ipdegdisp),iWork(ipnrdisp),work(ipHess),
     &               converged,Work(ipELEC),ielec,Work(ipelout),
     &               ldisp2,Lu_10)
           Call Niclas(work(ipHess),coor,Lu_10)
       End If
       Write(6,*)
       Write(6,*)
       Write(6,*)'************************************'
       Write(6,*)'*                                  *'
       Write(6,*)'*       Polarizabilities           *'
       Write(6,*)'*                                  *'
       Write(6,*)'************************************'
       Write(6,*)
       Write(6,*)
       Call Add_Info('POLARIZABILITIES',Pola,6,2)
*
*      Go from energy derivative to polarizability, there is a difference
*      in the sign in the definition.
*
       Call DScal_(6,-1.0D0,Pola,1)
*
       Call TriPrt(' ',' ',Pola,3)
       close(Lu_10)
       Call GetMem('NRDISP','FREE','INTE',ipnrdisp,ndisp)
       Call GetMem('DEGDISP','FREE','INTE',ipdegdisp,ndisp)
      End If
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*

      Call GetMem('Temp','FREE','REAL',ipTemp,nHss)
      Call GetMem('HESS','FREE','REAL',ipHess,nHss)
      Call GetMem('HESS','FREE','REAL',ipHess2,nHss)
      Call GetMem('RESPH','FREE','REAL',ipRHss,nHss)
      Call GetMem('Temp','FREE','REAL',ipELEC,3*ndisp)
      Call GetMem('Temp','Free','REAL',ipEG  ,3*ndisp)
      Call GetMem('Temp','FREE','REAL',ipELOUT,3*ndisp)
*
      Call QExit('Output')
      Return
      End
