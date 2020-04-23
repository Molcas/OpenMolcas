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
       SubRoutine OutPut_td(iKapDisp,isigdisp,iCiDisp,
     &                      iCiSigDisp,iRHSDisp,iRHSCIDisp,
     &                      converged)
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
       Character*20 Label2
       Integer Pstate_sym,ldisp2(8),ielec(3)
       Integer iKapDisp(nDisp),isigdisp(nDisp),
     &         iCiDisp(nDisp),iCiSigDisp(nDisp),
     &         iRHSDisp(nDisp),iRHSCiDisp(nDisp)
       Logical elec,converged(8),CI
       Real*8 Pola(6)
*
       debug=.false.
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
*
*          Calculate length of the density, Fock and Kappa matrix etc
*          notice that this matrixes not necessary are symmetric.
*          Store pointers.
*
*          Input:
*                 iSym : Symmetry of perturbation
*
*          Output: Commonblocks (Pointers.fh)
*
          Call Setup_MCLR(iSym)
          PState_SYM=iEor(State_Sym-1,iSym-1)+1
          nconfM=Max(ncsf(PState_Sym),nint(xispsm(Pstate_Sym,1)))
          nconf1=ncsf(PState_Sym)
          if (TimeDep) nconf1=nconf1*2
          if (TimeDep) nconfM=nconfM*2
          CI=.false.
          If (iMethod.eq.2.and.nconf1.gt.0) CI=.true.
          If (CI.and.nconf1.eq.1.and.isym.eq.1) CI=.false.

*
*    Allocate areas for scratch and state variables
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
          if (jspin.eq.0) then
           nconf1=ncsf(Pstate_sym)
          Else
           nconf1=nint(xispsm(Pstate_Sym,1))
          end if
          nCI=nconf1
          if (Timedep) nCI=nCI*2
          If (.not.lCalc(iDisp)) Goto 110


*
          iDisk=iKapDisp(iDisp)
          Len=nDensC
          Call dDaFile(LuTemp,2,Work(ipKap1),Len,iDisk)
          iDisk=iSigDisp(iDisp)
          Call dDaFile(LuTemp,2,Work(ipSKap),Len,iDisk)
          iDisk=iRHSDisp(iDisp)
          Call dDaFile(LuTemp,2,Work(iprKap1),Len,iDisk)
          Do i=0,ndensC-1
           Work(ipSKap+i)=-Work(ipSKap+i)-Work(iprKap1+i)
          End Do

          If (CI) Then
            ilen=nCI
            idis=iCIDisp(iDisp)
            Call dDaFile(LuTemp,2,Work(ipin(ipCIp1)),iLen,iDis)
            idis=iCISigDisp(idisp)
            Call dDaFile(LuTemp,2,Work(ipin(ipSp)),iLen,iDis)
            idis=iRHSCIDisp(idisp)
            Call dDaFile(LuTemp,2,Work(ipin(iprp1)),iLen,iDis)
            ii=ipin(ipsp)
            jj=ipin(iprp1)
            Do i=0,nCI -1
             Work(ii+i)=
     &             -Work(ii+i)-Work(jj+i)
            End Do
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
               If (iAnd(nTPert(kdisp+ksym),1).eq.1)
     &           kSpin=1
               if (kspin.eq.0) Then
                 nconf1=ncsf(PState_Sym)
               else
                 nConf1=nint(xispsm(Pstate_Sym,1))
               end if
               nCI=nconf1
               if (Timedep) nCI=nCI*2
               If (.not.lCalc(kDisp+ksym)) Goto 120
               iDisk=iKapDisp(kDisp+kSym)
               Len=nDensC
               Call dDaFile(LuTemp,2,Work(ipKap2),Len,iDisk)
*
*               Call Recprt('ipkap2',' ',Work(ipkap2),nDensC,1)
*               Call Recprt('ipSkap',' ',Work(ipSkap),nDensC,1)
              rTempk1=-1.0d0*DDot_(nDensC,Work(ipKap2),1,Work(ipSKap),1)
*
               if (CI) Then
                ilen=nCI
                idis=iCIDisp(kDisp+ksym)
                Call dDaFile(LuTemp,2,Work(ipin(ipCIp2)),iLen,iDis)
*
                rTempc1=DDot_(nCI,Work(ipin(ipCIp2)),1,
     &                              Work(ipin(ipsp)),1)
               Else
                rtempc1=0.0d0
               End If
*
               iDisk=iRHSDisp(kDisp+kSym)
               Len=nDensC
               Call dDaFile(LuTemp,2,Work(iprKap2),Len,iDisk)
               If (CI) Then
                 ilen=nCI
                 idis=iRHSCIDisp(kdisp+ksym)
                 Call dDaFile(LuTemp,2,Work(ipin(iprp2)),iLen,iDis)
               End If
*
               Fact=1.0d0
               If (kdisp.eq.jdisp) Fact=2.0d0
*               Call Recprt('ipkap1',' ',Work(ipkap1),nDensC,1)
*               Call Recprt('iprkap2',' ',Work(iprkap2),nDensC,1)
               rTempk2=-1.0d0*Fact*
     &                DDot_(nDensC,Work(ipKap1),1,Work(iprKap2),1)
               If (kdisp.ne.jdisp) Then
               rtempk3=-1.0d0*DDot_(nDensC,Work(iprKap1),1,
     &                  Work(ipKap2),1)
               Else
                rTempk3=0.0d0
               End if
               If (CI) Then
                 Fact=1.0d0
                 If (kdisp.eq.jdisp) Fact=2.0d0
                  rTempc2=Fact*
     &                DDot_(nCI   ,Work(ipin(ipCip1)),1,
     &                            Work(ipin(iprp2)),1)
                 If (kdisp.ne.jdisp) Then
                  rtempc3=1.0d0*
     &                DDot_(nCI,Work(ipin(iprp1)),1,
     &                            Work(ipin(ipCIp2)),1)
                 Else
                  rTempc3=0.0d0
                 End if
                Else
                  rtempc2=0.0d0
                  rtempc3=0.0d0
               End if
*
*              Write(*,*) kdisp,jdisp
*              Write(*,*)'rtempk', rTempk1,rtempk2,rtempk3
*              Write(*,*)'rtempc', rtempc1,rtempc2,rtempc3
*
               Maxi=Max(kDisp,jDisp)
               Mini=Min(kDisp,jDisp)
               index=mSym+Maxi*(Maxi-1)/2+Mini
*
               Work(ipRhss-1+Index)=Work(ipRhss-1+Index)+
     &             0.5d0*(rTempk1+rtempk2+rtempk3+
     &             rtempc1+rtempc2+rtempc3)
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
*
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
      Call GetMem('ELEC','ALLO','REAL',ipELEC,3*ndisp)
      Call GetMem('ELEC2','ALLO','REAL',ipELEC2,3*ndisp)
      Call GetMem('EG  ','ALLO','REAL',ipEG  ,3*ndisp)
      Call GetMem('ELOUT','ALLO','REAL',ipELOUT,3*ndisp)
      If (iMethod.eq.2)  irc=ipclose(-1)
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
      elec=.False.
      If (Mckinley) Then
         idum=1
         iopt=128
         irc=3*ndisp
         Label='DOTELGR'
         Call drdMCk(irc,iopt,LaBeL,idum,Work(ipEG),idum)
         elec=.true.
         if (irc.ne.0) elec=.false.
      End If
      call dcopy_(nhss,Work(iphss),1,Work(iphess2),1)
      If (debug) Then
         ip=ipHess2
         Do iSym=1,nSym
          Write(label2,'(A,I2)') 'CHessian symmetry',iSym
          If (lDisp(iSym).ne.0)
     &    Call TriPrt(label2,' ',Work(ip),lDisp(iSym))
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
*      Call recprt('iprhss',' ',Work(ipRHss),nHss,1)
*      Call recprt('iphess',' ',Work(ipHess),nHss,1)
C
      Call DaXpY_(mSym,1.0d0,Work(ipRHss),1,Work(ipHess2),1)
*
      If (debug) Then
       Call MMSORT2(Work(ipRHSS),Work(ipELEC),pola,ielec)
       Call Recprt('RESP',' ',Work(ipElec),3*nDisp,1)
      End If
*
      Call MMSORT2(Work(ipHESS2),Work(ipELEC),pola,ielec)
*
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
      Call mmSort(Work(ipHess2),Work(ipHess),ldisp2)

      If (McKinley) Then
*
        iRC=-1
        iOpt=0
        Label='StatHess'
        Call dRdMck(iRC,iOpt,Label,idum,Work(ipTemp),idum)
        If (iRC.ne.0) Then
           Write (6,*)
           Write (6,*) ' *** Error in subroutine OUTPUT_TD ***'
           Write (6,*) ' Reading from MCKINT file failed '
           Write (6,*)
        End If
*
        If (debug) Then
         ip=ipTemp
         Do iSym=1,nSym
          Write(label2,'(a,i2)') 'SHessian symmetry',iSym
          If (lDisp2(iSym).ne.0)
     &    Call TriPrt(label2,' ',Work(ip),lDisp2(iSym))
         ip=ip+ldisp2(isym)*(1+ldisp2(isym))/2
         End Do
        End If
        Call DaXpY_(mSym,1.0d0,Work(ipTemp),1,Work(ipHess),1)
      End If
      If (debug) Then
        ip=ipHess
        Do iSym=1,nSym
          Write(label2,'(a,i2)') 'Hessian symmetry',iSym
          If (lDisp2(iSym).ne.0)
     &    Call TriPrt(label2,' ',Work(ip),lDisp2(iSym))
          ip=ip+ldisp2(isym)*(1+ldisp2(isym))/2
        End Do
      End If
*
*
      If (McKinley) Then
       iRC=-1
       iOpt=0
       Label='Hess    '
       Call dWrMck(iRC,iOpt,Label,iDum,Work(ipHess),iDum)
       If (iRC.ne.0) Then
         Write (6,*)
         Write (6,*) ' *** Error in subroutine OUTPUT_TD ***'
         Write (6,*) ' Writing',Label,' to MCKINT file failed '
         Write (6,*)
       End If
*
       Call Put_iScalar('No of Internal coordinates',ldisp2(1))
       Call Put_AnalHess(Work(ipHess),ldisp2(1)*(ldisp2(1)+1)/2)
*
      End If
      If (.true.) Then
       iRC=-1
       iOpt=0
       Call GetMem('NRDISP','ALLO','INTE',ipnrdisp,ndisp)
       Call RdMck(irc,iopt,'NRCTDISP',idum,iWork(ipnrDisp),idum)
       iRC=-1
       iOpt=0
       Call GetMem('DEGDISP','ALLO','INTE',ipdegdisp,ndisp)
       Label='DEGDISP '
       Call RdMck(irc,iopt,Label,idum,iWork(ipDegDisp),idum)
       If (iRC.ne.0) Then
         Write (6,*)
         Write (6,*) ' *** Error in subroutine OUTPUT_TD ***'
         Write (6,*) 'Reading ',Label,' from MCKINT file failed '
         Write (6,*)
       End If
       If (debug)
     &  Call HssPrt_MCLR(iwork(ipdegdisp),Work(ipHess),ldisp2)
       call daxpy_(3*ndisp,-1.0d0,Work(ipEG),1,Work(ipELEC),1)
       If (debug.and.elec)
     &  Call Recprt('ELEC-ST',' ',Work(ipEG),3*nDisp,1)
       If (debug.and.elec)
     &  Call Recprt('ELEC-TOT',' ',Work(ipElec),3*nDisp,1)
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
*          Call Niclas(work(ipHess),coor)
       End If
       Write(6,*)
       Write(6,*)
       Write(6,*)'************************************'
       Write(6,*)'*                                  *'
       Write(6,*)'*       Time Dependent             *'
       Write(6,*)'*       Polarizabilities           *'
       Write(6,*)'*                                  *'
       Write(6,*)'************************************'
       Write(6,*)
       Write(6,*)
       Call Add_Info('TimeDep_Pol',Pola,6,2)
*
*      Go from energy derivative to polarizability, there is a difference
*      in the sign in the definition.
*
       Call DScal_(6,-1.0D0,Pola,1)
*
       Call TriPrt(' ',' ',Pola,3)
       close(Lu_10)
       Call GetMem('NRDISP','FREE','INTE',ipnrdisp,ndisp)
       Call GetMem('NRDISP','FREE','INTE',ipdegdisp,ndisp)
      End If
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
*      Call Getmem('output','CHECK','REAL',idum,idum)
      Call GetMem('Temp','FREE','REAL',ipTemp,nHss)
      Call GetMem('HESS','FREE','REAL',ipHess,nHss)
      Call GetMem('HESS','FREE','REAL',ipHess2,nHss)
      Call GetMem('RESPH','FREE','REAL',ipRHss,nHss)
      Call GetMem('EG  ','Free','REAL',ipEG  ,3*ndisp)
      Call GetMem('ELEC2','Free','REAL',ipELEC2,3*ndisp)
      Call GetMem('Temp','FREE','REAL',ipELEC,3*ndisp)
      Call GetMem('Temp','FREE','REAL',ipELOUT,3*ndisp)
*
      Return
      End
