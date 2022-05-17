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
       SubRoutine OutPut_MCLR(iKapDisp,isigdisp,iCiDisp,
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
      Use Arrays, only: Hss
      use ipPage, only: W
      Implicit Real*8 (a-h,o-z)
#include "detdim.fh"
#include "Input.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "disp_mclr.fh"
#include "cicisp_mclr.fh"
#include "stdalloc.fh"
      Character(LEN=8) Label
#ifdef _DEBUGPRINT_
      Character(LEN=20) Label2
#endif
      Integer Pstate_sym,ldisp2(8),ielec(3)
      Integer iKapDisp(nDisp),isigdisp(nDisp),
     &        iCiDisp(nDisp),iCiSigDisp(nDisp),
     &        iRHSDisp(nDisp),iRHSCiDisp(nDisp)
      Logical elec_On,converged(8),CI
      Real*8 Pola(6)
      Real*8, Allocatable:: RHss(:)
      Real*8, Allocatable:: Kap1(:), Kap2(:), sKap(:),
     &                     rKap1(:),rKap2(:)
      Real*8, Allocatable:: Hess(:), Hess2(:), Temp(:), ELEC(:),
     &                      EG(:), ELOUT(:)
      Integer, Allocatable:: NrDisp(:), DegDisp(:)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*
#ifdef _DEBUGPRINT_
       debug=.True.
#else
       debug=.false.
#endif
      nHss=SIZE(Hss)
      nhess=nDisp*(nDisp+1)/2
      Call mma_allocate(RHss,nHss,Label='RHss')
      RHss(:)=0.0d0
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
          Call mma_allocate(Kap1,nDensC,Label='Kap1')
          Call mma_allocate(Kap2,nDensC,Label='Kap2')
          Call mma_allocate(sKap,nDensC,Label='sKap')
          Call mma_allocate(rkap1,nDensC,Label='rKap1')
          Call mma_allocate(rkap2,nDensC,Label='rKap2')
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
              Call dDaFile(LuTemp,2,Kap1,Len,iDisk)
              iDisk=iSigDisp(iDisp)
              Call dDaFile(LuTemp,2,SKap,Len,iDisk)
              iDisk=iRHSDisp(iDisp)
              Call dDaFile(LuTemp,2,rKap1,Len,iDisk)
              Do i=1,ndensC
                 SKap(i)=-SKap(i)-rKap1(i)
              End Do
C
*             Call Recprt('ORB-RHS',' ',rKap1,nDensC,1)
*             Write(*,*)'ddot orb-resp',
*     &            ddot_(ndensC,Kap1,1,Kap1,1)
*             Write(*,*)'ddot orb-sigma',
*     &            ddot_(ndensC,SKap,1,SKap,1)
*             Write(*,*)'ddot orb-rhs',
*     &            ddot_(ndensC,rKap1,1,rKap1,1)
C
*
             Call GADSum(Kap1,Len)
             Call GADSum(SKap,Len)
             Call GADSum(rKap1,Len)

             If (CI) Then
                ilen=nconf1
                idis=iCIDisp(iDisp)
                irc=ipin(ipCIp1)
                Call dDaFile(LuTemp,2,W(ipCIp1)%Vec,iLen,iDis)
                idis=iCISigDisp(idisp)
                irc=ipin(ipSp)
                Call dDaFile(LuTemp,2,W(ipSp)%Vec,iLen,iDis)
                idis=iRHSCIDisp(idisp)
                irc=ipin(iprp1)
                Call dDaFile(LuTemp,2,W(iprp1)%Vec,iLen,iDis)
                ii=ipin(ipSp)
                jj=ipin(iprp1)
                Do i=1,nConf1
                   W(ipSp)%Vec(i)= -W(ipSp)%Vec(i)-W(iprp1)%Vec(i)
                End Do
C
*               Write(*,*)'ddot ci-resp',
*     &               ddot_(nConf1,W(ipcip1)%Vec,1,
*     &                            W(ipcip1)%Vec,1)
*               Write(*,*)'ddot ci-sigma',
*     &               ddot_(nConf1,W(ipSp)%Vec,1,
*     &                            W(ipSp)%Vec,1)
*               Write(*,*)'ddot ci-rhs',
*     &               ddot_(nConf1,W(iprp1)%Vec,1,
*     &                            W(iprp1)%Vec,1)
C
                Call GADSum(W(ipCIp1)%Vec,iLen)
                Call GADSum(W(ipSp  )%Vec,iLen)
                Call GADSum(W(ipRp1 )%Vec,iLen)
             End If
*
          Else
*
                            Len = nDensC
             Kap1(1:Len)=0.0D0
             Call GADSum(Kap1,Len)
             sKap(1:Len)=0.0D0
             Call GADSum(SKap,Len)
             rKap1(1:Len)=0.0D0
             Call GADSum(rKap1,Len)
*
             If (CI) Then
*
                ilen=nconf1
                irc=ipin(ipCIp1)
                Call FZero(W(ipCIp1)%Vec,iLen)
                Call GADSum(W(ipCIp1)%Vec,iLen)
                irc=ipin(ipSp  )
                Call FZero(W(ipSp  )%Vec,iLen)
                Call GADSum(W(ipSp  )%Vec,iLen)
                irc=ipin(ipRp1 )
                Call FZero(W(ipRp1 )%Vec,iLen)
                Call GADSum(W(ipRp1 )%Vec,iLen)
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
                  Call dDaFile(LuTemp,2,Kap2,Len,iDisk)
                  iDisk=iRHSDisp(kDisp+kSym)
                  Call dDaFile(LuTemp,2,rKap2,Len,iDisk)

                  Call GASync()
                  Call GADSum(Kap2,Len)
                  Call GADSum(rKap2,Len)

                  If (CI) Then
                     ilen=nconf1
                     idis=iCIDisp(kDisp+ksym)
                     irc=ipin(ipCIp2)
                     Call dDaFile(LuTemp,2,W(ipCIp2)%Vec,iLen,iDis)
                     idis=iRHSCIDisp(kdisp+ksym)
                     irc=ipin(iprp2)
                     Call dDaFile(LuTemp,2,W(iprp2)%Vec,iLen,iDis)
                     irc=ipin(ipsp)
                     rTempc1=DDot_(nConf1,W(ipCIp2)%Vec,1,W(ipsp)%Vec,1)

                     Call GASync() ! <----------------- NOTE!
                     Call GADSum(W(ipCIp2)%Vec,iLen)
                     Call GADSum(W(iprp2 )%Vec,iLen)

                  Else
                     rtempc1=0.0d0
                  End If
*
               Else
*
                  Call GASync()
                  Len=nDensC
                  Call FZero(Kap2,Len)
                  Call GADSum(Kap2,Len)
                  Call FZero(rKap2,Len)
                  Call GADSum(rKap2,Len)
                  If (CI) Then
                     ilen=nconf1
                     Call GASync()   ! <----------------- NOTE!
                     irc=ipin(ipCIp2)
                     Call FZero(W(ipCIp2)%Vec,iLen)
                     Call GADSum(W(ipCIp2)%Vec,iLen)
                     irc=ipin(iprp2 )
                     Call FZero(W(iprp2 )%Vec,iLen)
                     Call GADSum(W(iprp2 )%Vec,iLen)
                     irc=ipin(ipsp)
                     rTempc1=DDot_(nConf1,W(ipCIp2)%Vec,1,W(ipsp)%Vec,1)
                  Else
                     rtempc1=0.0d0
                  End If
*
               End If
*
               rTempk1=DDot_(nDensC,Kap2,1,SKap,1)
*
               Fact=1.0d0
               If (kdisp.eq.jdisp) Fact=2.0d0
               rTempk2=Fact*DDot_(nDensC,Kap1,1,rKap2,1)
               If (kdisp.ne.jdisp) Then
               rtempk3=1.0d0*DDot_(nDensC,rKap1,1,Kap2,1)
               Else
                rTempk3=0.0d0
               End if
               If (CI) Then
                 Fact=1.0d0
                 If (kdisp.eq.jdisp) Fact=2.0d0
                  irc=ipin(ipCip1)
                  irc=ipin(iprp2)
                  rTempc2=Fact*
     &                DDot_(nConf1,W(ipCip1)%Vec,1,W(iprp2)%Vec,1)
                 If (kdisp.ne.jdisp) Then
                  irc=ipin(iprp1)
                  irc=ipin(ipCIp2)
                  rtempc3=1.0d0*
     &                DDot_(nConf1,W(iprp1)%Vec,1,W(ipCIp2)%Vec,1)
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
               Rhss(Index)=Rhss(Index)+
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
          Call mma_deallocate(rKap2)
          Call mma_deallocate(rKap1)
          Call mma_deallocate(sKap)
          Call mma_deallocate(Kap2)
          Call mma_deallocate(Kap1)
          If (CI)  irc=ipclose(ipcip1)
 100  Continue

      Call mma_allocate(Hess,nHss,Label='Hess')
      Call mma_allocate(Hess2,nHss,Label='Hess2')
      Call mma_allocate(Temp,nHss,Label='Temp')
      Temp(:)=0.0D0
      Call mma_allocate(ELEC,3*ndisp,Label='ELEC')
      Call mma_allocate(EG  ,3*ndisp,Label='EG')
      Call mma_allocate(ELOUT,3*ndisp,Label='ELOUT')
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
      Call drdMCk(irc,iopt,LaBeL,idum,EG,idum)
      elec_On=.true.
      if (irc.ne.0) elec_On=.false.
                Call GADsum(Hss,nHss)
      call dcopy_(nHss,Hss,1,Hess2,1)
#ifdef _DEBUGPRINT_
      If (debug) Then
         ip=1
         Do iSym=1,nSym
           Write(label2,'(A,I2)') 'CHessian symmetry',iSym
            If (lDisp(iSym).ne.0)
     &         Call TriPrt(label2,' ',Hess2(ip),lDisp(iSym))
            ip=ip+ldisp(isym)*(1+ldisp(isym))/2
         End Do
      End If
*
      If (debug) Then
       Call MMSORT2(HESS2,ELEC,pola,ielec)
       Call Recprt('CONN',' ',Elec,3*nDisp,1)
      End If
*
C
c       Write(*,*)'I am here 1'
       Call Recprt('Rhss','(5G20.10) ',RHss,nhss,1)
       Call Recprt('Hss','(5G20.10) ',Hss,nHss,1)
#endif
C
      Call DaXpY_(mSym,1.0d0,RHss,1,Hess2,1)
*
#ifdef _DEBUGPRINT_
      If (debug) Then
       Call MMSORT2(RHSS,ELEC,pola,ielec)
       Call Recprt('RESP',' ',Elec,3*nDisp,1)
      End If
#endif
*
      Call MMSORT2(HESS2,ELEC,pola,ielec)
*
#ifdef _DEBUGPRINT_
      If (debug) Then
       Call Recprt('R+C',' ',Elec,3*nDisp,1)
       ip=1
       Do iSym=1,nSym
        Write(label2,'(A,I2)') 'Hessian symmetry',iSym
        If (lDisp(iSym).ne.0)
     &   Call TriPrt(label2,' ',Hess2(ip),lDisp(iSym))
        ip=ip+ldisp(isym)*(1+ldisp(isym))/2
       End Do
      End If
#endif
*
      Call mmSort(Hess2,Hess,ldisp2)
*
      If (McKinley) Then
*
         iRC=-1
         iOpt=0
         Label='StatHess'
         Call dRdMck(iRC,iOpt,Label,idum,Temp,idum)
         If (iRC.ne.0) Then
            Write (6,*) 'OutPut: Error reading MCKINT'
            Write (6,'(A,A)') 'Label=',Label
            Call Abend()
         End If
*
#ifdef _DEBUGPRINT_
         If (debug) Then
            ip=1
            Do iSym=1,nSym
               Write(label2,'(a,i2)') 'SHessian symmetry',iSym
               If (lDisp2(iSym).ne.0)
     &            Call TriPrt(label2,' ',Temp(ip),lDisp2(iSym))
               ip=ip+ldisp2(isym)*(1+ldisp2(isym))/2
            End Do
         End If
#endif
         Call DaXpY_(mSym,1.0d0,Temp,1,Hess,1)
      End If
#ifdef _DEBUGPRINT_
      If (debug) Then
        ip=1
        Do iSym=1,nSym
          Write(label2,'(a,i2)') 'Hessian symmetry',iSym
          If (lDisp2(iSym).ne.0)
     &    Call TriPrt(label2,' ',Hess(ip),lDisp2(iSym))
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
         Call dWrMck(iRC,iOpt,Label,iDum,Hess,iDum)
         If (iRC.ne.0) Then
            Write (6,*) 'OutPut: Error writing to MCKINT'
            Write (6,'(A,A)') 'Label=',Label
            Call Abend()
         End If
         Call Put_iScalar('No of Internal coordinates',ldisp2(1))
         Call Put_AnalHess(Hess,ldisp2(1)*(ldisp2(1)+1)/2)
      End If
*
      If (.true.) Then
       iRC=-1
       iOpt=0
       Call mma_allocate(NrDisp,ndisp,Label='NrDisp')
       Call RdMck(irc,iopt,'NRCTDISP',idum,NrDisp,idum)
       iRC=-1
       iOpt=0
       Call mma_allocate(DegDisp,ndisp,Label='DegDisp')
       Label='DegDisp'
       Call RdMck(irc,iopt,Label,idum,DegDisp,idum)
       If (iRC.ne.0) Then
          Write (6,*) 'OutPut: Error reading RELAX'
          Write (6,'(A,A)') 'Label=',Label
          Call Abend()
       End If
C
*       If (debug)
*        Call HssPrt_MCLR(DegDisp,Hess,ldisp2)
*       Call Recprt('hess',' ',Hess,nhss,1)
C
       call daxpy_(3*ndisp,-1.0d0,EG,1,ELEC,1)
#ifdef _DEBUGPRINT_
       If (debug.and.elec_On)
     &  Call Recprt('ELEC-ST',' ',EG,3*nDisp,1)
       If (debug.and.elec_On)
     &  Call Recprt('ELEC-TOT',' ',Elec,3*nDisp,1)
#endif
*
       Lu_10=10
       Lu_10=IsFreeUnit(Lu_10)
       call molcas_open(lu_10,'UNSYM')
c       Open(unit=Lu_10, file='UNSYM')
*
       If (Mckinley) Then
           Call FreqAnal(DegDisp,NrDisp,Hess,
     &               converged,ELEC,ielec,elout,
     &               ldisp2,Lu_10)
           Call Niclas(Hess,coor,Lu_10)
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
       Call mma_deallocate(NrDisp)
       Call mma_deallocate(DegDisp)
      End If
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Call mma_deallocate(ELOUT)
      Call mma_deallocate(EG)
      Call mma_deallocate(ELEC)
      Call mma_deallocate(Temp)
      Call mma_deallocate(Hess2)
      Call mma_deallocate(Hess)
      Call mma_deallocate(RHss)
*
      Return
      End
