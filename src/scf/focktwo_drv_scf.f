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
      Subroutine FockTwo_Drv_scf(nSym,nBas,nAux,Keep,
     &                       DLT,DSQ,FLT,nFLT,
     &                       ExFac,nBSQT,nBMX,iUHF,DLT_ab,
     &                       DSQ_ab,FLT_ab,nOcc,nOcc_ab,iDummy_run)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
      Integer nSym,nBas(8), nAux(8), Keep(8)
      Integer nOcc(nSym),nOcc_ab(nSym)
      Logical DoCholesky,GenInt,DoLDF
      Integer ALGO,nSCREEn
      Logical REORD,DECO,timings
      Real*8 DLT(*),DSQ(*),FLT(nFLT),dmpk,dFKmat
      Real*8 DLT_ab(*),DSQ_ab(*),FLT_ab(*)
      Character*50 CFmt

      Common /CHOSCF / REORD,DECO,dmpk,dFKmat,ALGO,nSCREEN
      Common /CHOTIME / timings
      Logical Do_OFemb, KEonly, OFE_first
      COMMON  / OFembed_L / Do_OFemb,KEonly,OFE_first
      COMMON  / OFembed_R / Rep_EN,Func_AB,Func_A,Func_B,Energy_NAD,
     &                      V_Nuc_Ab,V_Nuc_BA,V_emb
      COMMON  / OFembed_I / ipFMaux, ip_NDSD, l_NDSD
*

      GenInt=.false.
      DoCholesky=.false.
      if(ALGO.eq.0) GenInt=.true. !use GenInt to regenerate integrals
      Call DecideOnCholesky(DoCholesky)
      Call DecideOnLocalDF(DoLDF)
      If (DoLDF) GenInt=.false. ! regenerate ints not impl. for LDF

c      write(6,*)'*************************'
c      write(6,*)'ONLY COULOMB CONTRIBUTION'
c      write(6,*)'*************************'
c      exFac=0.d0
c      write(6,*)'ExFac= ',ExFac
*
      If (Do_OFemb) Then ! Coul. potential from subsys B
         nFM=1
         If (iUHF.eq.1) nFM=2
         If (OFE_first) Call GetMem('FMaux','Allo','Real',ipFMaux,nFlt)
         Call Coul_DMB(OFE_first,nFM,Rep_EN,Work(ipFMaux),
     &                 DLT,DLT_ab,nFlt)
         OFE_first=.false.
      End If
*
      Call GetMem('LWFSQ','Allo','Real',LWFSQ,NBSQT)
C zeroing the elements
      call dcopy_(NBSQT,Zero,0,Work(LWFSQ),1)

      if((.not.DoCholesky).or.(GenInt)) then
      Call GetMem('LW2','Allo','Real',LW2,NBMX*NBMX)
      end if
*
* nFlt is the total dimension of the LT fock matrix
      Call Getmem('tempFLT','Allo','Real',ipTemp,nFlt)
      Call FZero(Work(ipTemp),nFlt)
*
      IF(iUHF.eq.1) THEN
        Call GetMem('LWFSQ_ab','Allo','Real',LWFSQ_ab,NBSQT)
        call dcopy_(NBSQT,Zero,0,Work(LWFSQ_ab),1)
        Call Getmem('FLT_ab','Allo','Real',ipTemp_ab,nFlt)
        Call FZero(Work(ipTemp_ab),nFlt)
*
        if((.not.DoCholesky).or.(GenInt)) then
          Call GetMem('LW2_ab','Allo','Real',LW2_ab,NBMX*NBMX)
        endif
*
      ENDIF
*
*
      Call GetMem('LW1','MAX','Real',LW1,LBUF)
*
* Standard building of the Fock matrix from Two-el integrals
*
      Call CWTIME(TotCPU1,TotWALL1)

      IF (.not.DoCholesky) THEN
         Call GetMem('LW1','Allo','Real',LW1,LBUF)
*
       If (LBUF.LT.NBMX**2) Then
         WRITE(6,*)'FockTwo_Drv_SCF Error: Too little memory remains'
     &     //' for the call to FOCKTWO_SCF.'
         WRITE(6,*)' Largest allocatable array size LBUF=',LBUF
         WRITE(6,*)' Max nr of bf in any symmetry,  NBMX=',NBMX
         WRITE(6,*)' Required minimum size       NBMX**2=',NBMX**2
         WRITE(6,*)'    (All in Real*8-size words)'
         Call QTRACE()
         Call  ABEND()
       End If
*
       If (iUHF.eq.1) Then

         Call FOCKTWO_scf(nSym,nBas,nAux,Keep,
     &             DLT,DSQ,Work(ipTemp),nFlt,
     &             Work(LWFSQ),LBUF,Work(LW1),Work(LW2),ExFac,iUHF,
     &             DLT_ab,DSQ_ab,Work(ipTemp_ab),Work(LWFSQ_ab))

       Else  ! RHF calculation

         Call FOCKTWO_scf(nSym,nBas,nAux,Keep,
     &             DLT,DSQ,Work(ipTemp),nFlt,
     &             Work(LWFSQ),LBUF,Work(LW1),Work(LW2),ExFac,iUHF,
     &             Work(ip_Dummy),Work(ip_Dummy),Work(ip_Dummy),
     &             Work(ip_Dummy))

       EndIf

      ENDIF
*
* Building of the Fock matrix regenerating the integrals on the fly
*
      IF ((DoCholesky).and.(GenInt)) THEN ! save some space for GenInt
         LBUF = MAX(LBUF-LBUF/10,0)
         Call GetMem('LW1','Allo','Real',LW1,LBUF)
*
       If (LBUF.LT.NBMX**2) Then
         WRITE(6,*)' FockTwo_Drv Error: Too little memory remains for'
     &     //' the call to FOCKTWO_SCF.'
         WRITE(6,*)' Largest allocatable array size LBUF=',LBUF
         WRITE(6,*)' Max nr of bf in any symmetry,  NBMX=',NBMX
         WRITE(6,*)' Required minimum size       NBMX**2=',NBMX**2
         WRITE(6,*)'    (All in Real*8-size words)'
         Call QTRACE()
         Call  ABEND()
       End If
*
       If (iUHF.eq.1) Then

         Call FOCKTWO_scf(nSym,nBas,nAux,Keep,
     &             DLT,DSQ,Work(ipTemp),nFlt,
     &             Work(LWFSQ),LBUF,Work(LW1),Work(LW2),ExFac,iUHF,
     &             DLT_ab,DSQ_ab,Work(ipTemp_ab),Work(LWFSQ_ab))

       Else  ! RHF calculation

         Call FOCKTWO_scf(nSym,nBas,nAux,Keep,
     &             DLT,DSQ,Work(ipTemp),nFlt,
     &             Work(LWFSQ),LBUF,Work(LW1),Work(LW2),ExFac,iUHF,
     &             Work(ip_Dummy),Work(ip_Dummy),Work(ip_Dummy),
     &             Work(ip_Dummy))

       EndIf

      ENDIF
*
      Call CWTIME(TotCPU2,TotWALL2)
      TOTCPU  = TotCPU2 - TotCPU1
      TOTWALL = TotWALL2 - TotWALL1

*---- Timings information for conventional or Cholesky with ALGO=0
      IF ( .not.DoCholesky .or. (DoCholesky.and.GenInt) ) THEN
      if(timings)then
      CFmt='(2x,A)'
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      if(DoCholesky)then
      Write(6,CFmt)'---    Cholesky SCF - Integral regeneration   ---'
      else
      Write(6,CFmt)'-----------     Conventional SCF     ------------'
      endif
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,CFmt)'Fock matrix construction        CPU       WALL   '
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
         Write(6,'(2x,A26,2f10.2)')'TOTAL                            '
     &                           //'         ',TOTCPU,TOTWALL
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)
      endif
      ENDIF
*
*
* Building of the Fock matrix directly from Cholesky vectors
*

      IF (DoCholesky .and. .not.GenInt.and.iDummy_run.eq.1) THEN
       Write(6,*) '*** Warning: missing feature in Cholesky code'
       Write(6,*) 'Use the results with extra care!'
      endif

      IF (DoCholesky .and. .not.GenInt.and.iDummy_run.eq.0) THEN
*
        If (iUHF.eq.1) Then

           CALL CHOscf_drv(iUHF,nSym,nBas,DSQ,DLT,DSQ_ab,DLT_ab,
     &                 Work(ipTemp),Work(ipTemp_ab),nFLT,ExFac,
     &                 LWFSQ,LWFSQ_ab,nOcc,nOcc_ab)
        Else
*
           CALL CHOscf_drv(iUHF,nSym,nBas,DSQ,DLT,
     &                 Work(ip_Dummy),Work(ip_Dummy),
     &                 Work(ipTemp),Work(ip_Dummy),nFLT,ExFac,
     &                 LWFSQ,ip_Dummy,nOcc,iWork(ip_iDummy))
        EndIf

      ENDIF
*
      Call DaXpY_(nFlt,One,Work(ipTemp),1,FLT,1)
      if(iUHF.eq.1) then
        Call DaXpY_(nFlt,One,Work(ipTemp_ab),1,FLT_ab,1)
      endif
*
      Call GetMem('tempFLT','Free','Real',ipTemp,nFlt)
      if(iUHF.eq.1) then
       Call GetMem('FLT_ab','Free','Real',ipTemp_ab,nFlt)
      endif
*
      If (Do_OFemb) Then ! add FM from subsystem B
        Call DaXpY_(nFlt,One,Work(ipFMaux),1,FLT,1)
        If (iUHF.eq.1) Call DaXpY_(nFlt,One,Work(ipFMaux),1,FLT_ab,1)
      EndIf
*
      IF ((.not.DoCholesky).or.(GenInt)) THEN
      Call GetMem('LW1','Free','Real',LW1,LBUF)
      Call GetMem('LW2','Free','Real',LW2,NBMX*NBMX)
      END IF

      Call GetMem('LWFSQ','Free','Real',LWFSQ,NBSQT)
      if(iUHF.eq.1) then
         IF ((.not.DoCholesky).or.(GenInt)) THEN
            Call GetMem('LW2_ab','Free','Real',LW2_ab,NBMX*NBMX)
         END IF
       Call GetMem('LWFSQ_ab','Free','Real',LWFSQ_ab,NBSQT)
      endif
*
      Return
      End
