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
* Copyright (C) 2021, Yoshio Nishimoto                                 *
************************************************************************
      Subroutine OLagVVVO(iSym,DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,T2AO,
     *                    DIA,DI,FIFA,FIMO,A_PT2,NumCho)
      USE iSD_data
      USE CHOVEC_IO
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "caspt2_grad.fh"
#include "nsd.fh"
C
      DIMENSION DPT2AO(*),DPT2CAO(*),FPT2AO(*),FPT2CAO(*),T2AO(*)
      Dimension DIA(*),DI(*),FIFA(*),FIMO(*),A_PT2(*)
C
      DIMENSION nBasX(8),KEEP(8)
      logical dorys
      Allocatable :: T_hbf(:,:,:,:),iOffAO(:)
      Character*4096 RealName
      Logical DoCholesky,is_error,Square
C
C     ----- (VV|VO)
C
C     Compute L_{pq} = (pa|jb) * T_{qj}^{ab}, in particular for
C     (p,q) \in  (virt, inact+act). This operation involves
C     (VV|VO) integrals which are stored in neither disk nor memory.
C     The back-transformation is applied to occupied orbital derivatives
C     of two-electron integrals that have two virtual indices (F, G, H).
C     In principle, the algorithm is to avoid (VV|VO) integrals,
C     i.e. U_{pq} for (p,q) = (vir, inact+act), but can also be applied
C     to U_{pq} for (p,q) = (all, inact+act).
C
      Call GetMem('WRK1','Allo','Real',ipWRK1,nBasT*nBasT)
      Call GetMem('WRK2','Allo','Real',ipWRK2,nBasT*nBasT)
      Call GetOrd(IRC,Square,nSymX,nBasX,KEEP)
C
      ! nTot1 = nBasT*(nBasT+1)/2
      nTot2 = nBast*nBasT
      ipvLag = ipWRK1
      Call DCopy_(nBasT*nBasT,[0.0D+00],0,Work(ipvLag),1)
C
      Call DecideOncholesky(DoCholesky)
      If (DoCholesky) Then
        !! No need to save CMOPT2. Just save A_PT2 and B_PT2.
        !! First, save A_PT2 in LuCMOPT2
        Call PrgmTranslate('CMOPT2',RealName,lRealName)
        If (IFMSCOUP.and.jState.ne.1) Then
    !       Call MOLCAS_Open_Ext2(LuCMOPT2,RealName(1:lRealName),
    !  &                          'DIRECT','UNFORMATTED',
    !  &                          iost,.FALSE.,
    !  &                          1,'OLD',is_error)
          Call GetMem('WRK3','ALLO','Real',ipWRK3,NumCho*NumCho)
    !       Read (LuCMOPT2) Work(ipWRK3:ipWRK3+NumCho**2-1)

          ! read A_PT2 from LUAPT2
          id = 0
          call ddafile(LUAPT2, 2, Work(ipWRK3), NumCho**2, id)

          Call DaXpY_(NumCho*NumCho,1.0D+00,Work(ipWRK3),1,A_PT2,1)
          Call GetMem('WRK3','FREE','Real',ipWRK3,NumCho*NumCho)
          ! REWIND LuCMOPT2
        Else
    !       Call MOLCAS_Open_Ext2(LuCMOPT2,RealName(1:lRealName),
    !  &                         'DIRECT','UNFORMATTED',
    !  &                          iost,.FALSE.,
    !  &                          1,'REPLACE',is_error)
        End If
C       write(6,*) "write...",numcho,lucmopt2

        ! For SS-CASPT2 I should write A_PT2 on disk only
        ! for the correct iRlxRoot
        if (jState.eq.iRlxRoot .or. nStLag.gt.1) then
          ! Write (LuCMOPT2) A_PT2(1:NumCho**2)

          ! write A_PT2 in LUAPT2
          id = 0
          call ddafile(LUAPT2,1,A_PT2,NumCho**2,id)

        end if
        ! Close (LuCMOPT2)

        !! Prepare for saving B_PT2. B_PT2 is saved in VVVOX2
        Call PrgmTranslate('GAMMA',RealName,lRealName)

        ! open LuGamma here, as we will write on it in vvvo_drv
        If (IFMSCOUP.and.jState.ne.1) Then
          Call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),
     &                         'DIRECT','UNFORMATTED',
     &                          iost,.TRUE.,
     &                          nBas(iSym)**2*8,'OLD',is_error)
        Else
          ! open only for the correct root
          if (jState.eq.iRlxRoot .or. nStLag.gt.1) then
            Call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),
     &                           'DIRECT','UNFORMATTED',
     &                            iost,.TRUE.,
     &                            nBas(iSym)**2*8,'REPLACE',is_error)
          end if
        End If
      End If
      !! 2) Compute ERI (mu rho | nu sigma)
      !! 3) Quarter-transformation of ERI
      !!    (mu rho | j sigma) = sum_{j} C_{nu j} (mu rho | nu sigma)
      !! 4) Contract with AO amplitude
      !!    L_{mu i} = sum_{j,rho sigma} T_{ij}^{mu nu}*(mu rho|j sigma)
      !! the third argument is used wrongly, but should be no problem

      !! D_{mu nu} := DPT2AO
      !! FPT2AO = ((mu nu|rho sigma)-(mu rho|nu sigma)/4
      !!          -(mu sigma|nu rho)/4)*D_{rho sigma}
      isymi = 1
      isymj = 1
      isyma = 1
      isymb = 1
C     nocc = nfro(1)+nish(1)+nash(1)
      nocc = nish(1)+nash(1)
      Call VVVO_Drv(nSym,nBas,nIsh,nFro,KEEP,
     *              iSym,iSymI,iSymA,iSymJ,iSymB,
     *              T2AO,Work(ipvLag),
     *              nOcc,nBasT,nTot2,nBMX,
     *              Work(LCMOPT2+nBasT*nFro(iSymA)),
     *              DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,
     *              DIA,DI,FIFA,FIMO)
C     write(6,*) "fpt2cao"
C     call sqprt(fpt2cao,12)
C
      !! Save the half transformed integrals on disk.
      !! It will be used in drvg1.f etc for gradient.
C     write(6,*) "1"
C     write(6,*) "lugamma = ", lugamma
C     !! MP2: DaName_MF_MA, MRCI: DaName_MF
C     Call DaName_MF(LuGamma,'GAMMA')
C     iDisk = 0
C     write(6,*) "2"
C     write(6,*) "nbast = ", nbast
C     CALL DaName_MF_WA(LuGamma,'GAMMA')
C     Call Molcas_BinaryOpen_Vanilla(LuGamma,'GAMMA')
C     write(6,*) "lugamma mod = ", lugamma
C     write(6,*) "aa"
C     Call dDaFile(LuGamma,1,Work(LCMOPT2),nBasT*nBasT,iDisk)
      If (DoCholesky) Then
        if (jState.eq.iRlxRoot .or. nStLag.gt.1) then
          Close (LuGamma)
        end if
      Else
        !! This is only for conventional calculations!!
        Call PrgmTranslate('CMOPT2',RealName,lRealName)
        Call MOLCAS_Open_Ext2(LuCMOPT2,RealName(1:lRealName),
     &                       'DIRECT','UNFORMATTED',
     &                        iost,.FALSE.,
     &                        1,'REPLACE',is_error)
        !! First, CMOPT2 has to be saved. The MO coefficient matrix in
        !! grvg1.f may be different from CMOPT2.
C       Write (LuCMOPT2) (Work(LCMOPT2+i-1),i=1,nBasT*nBasT)
        Do i = 1, nBasT*nBasT
          Write (LuCMOPT2) Work(LCMOPT2+i-1)
C         Write (LuCMOPT2) Work(LCMO+i-1)
        End Do
        write (LuCMOPT2) (nIsh(1)+nAsh(1))
        write (LuCMOPT2) (nIsh(2)+nAsh(2))
        write (LuCMOPT2) (nIsh(3)+nAsh(3))
        write (LuCMOPT2) (nIsh(4)+nAsh(4))
        write (LuCMOPT2) (nIsh(5)+nAsh(5))
        write (LuCMOPT2) (nIsh(6)+nAsh(6))
        write (LuCMOPT2) (nIsh(7)+nAsh(7))
        write (LuCMOPT2) (nIsh(8)+nAsh(8))
        write (LuCMOPT2) (nFro(1))
        write (LuCMOPT2) (nFro(2))
        write (LuCMOPT2) (nFro(3))
        write (LuCMOPT2) (nFro(4))
        write (LuCMOPT2) (nFro(5))
        write (LuCMOPT2) (nFro(6))
        write (LuCMOPT2) (nFro(7))
        write (LuCMOPT2) (nFro(8))
        !! Number of state-specific density matrix
        !! It is needed, because the separable contribution to the
        !! electron repulsion integral is different
        nSSDM = 0
        If (.not.IFSADREF) Then
          If (nState.eq.1) Then
            nSSDM = 0
          Else If (IFDW .or. IFRMS) Then
            !! For (X)DW, use nState density matrix
            nSSDM = nState
          Else If (IFXMS) Then
            !! For XMS, use SA density matrix
            nSSDM = 0
          Else If (IFMSCOUP) Then
            !! For MS, use nState density matrix
            nSSDM = nState
          Else
            !! Otherwise, use one SS density matrix
            nSSDM = 1
          End If
        End If
        write (LuCMOPT2) nSSDM
C
        !! Write the number of occupied orbitals
C       write (LuCMOPT2) (dble(nIsh(iSym)+nAsh(iSym)),iSym=1,8)
C       write (LuCMOPT2) (dble(nFro(iSym)),iSym=1,8)
C       do isym = 1, 8
C       write (LuCMOPT2) dble(nIsh(iSym)+nAsh(iSym))
C       end do
C       do isym = 1, 8
C       write (LuCMOPT2) dble(nFro(iSym))
C       end do
        Close (LuCMOPT2)
C       write(6,*) "mo saved"
C       call sqprt(Work(LCMOPT2),nbast)
C
C       write(6,*) "going to open LuGamma"
        Call PrgmTranslate('GAMMA',RealName,lRealName)
        if (ifmscoup.and.jstate.ne.1) then
          Call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),
     &                         'DIRECT','UNFORMATTED',
     &                          iost,.TRUE.,
     &                          nOcc*nOcc*8,'OLD',is_error)
        else
          Call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),
     &                         'DIRECT','UNFORMATTED',
     &                          iost,.TRUE.,
     &                          nOcc*nOcc*8,'REPLACE',is_error)
        end if
        if (is_error) then
         write (6,*) "Something is wrong in opening LuGamma in olagvvvo"
          call abend
        end if
C       Call dDaFile(LuGamma,1,Work(LCMOPT2),nBasT*nBasT,iDisk)
C       write(6,*) "Idisk = ", idisk
C       Call dDaFile(LuGamma,1,1.0d+00,1,iDisk)
C       write(6,*) "Idisk = ", idisk
C       write(6,*) "3"
        !  Setup for shell. Why do I have to call IniSew damn here?
        !  The number of shells should be able to be referred globally.
        nDiff=1
        DoRys=.True.
        Call IniSew(DoRys,nDiff)
        Call Nr_Shells(nSkal)
        Call Setup_iSD()
C       write(6,*) "4"
C       write(6,*) "nskal = ", nskal
        !! see Include/info.fh
C       call bbbb
        Allocate (iOffAO(nSkal+1))
        MaxShlAO = 0
        iOffAO(1) = 0
        Do iSh = 1, nSkal
          nBasI = iSD(2,iSh)*iSD(3,iSh)
          If (nBasI.gt.MaxShlAO) MaxShlAO = nBasI
          iOffAO(iSh+1) = iOffAO(iSh)+nBasI
        End Do
        ! nMax = MaxShlAO*MaxShlAO*nOcc*nOcc
C       write(6,*) "maxshlao = ", maxshlao

        Allocate (T_hbf(nOcc,nOcc,MaxShlAO,MaxShlAO))
C       write(6,*) "nocc = ", nocc
        Do iSh = 1, nSkal
          !! iSD(2,iSh): number of AOs of the shell
          !! iSD(3,iSh): number of cont. func. of the shell
          nBasI = iSD(2,iSh)*iSD(3,iSh)
C         write(6,*) isd(2,ish),isd(3,ish),nbasi
C         write(6,*) "ioffao = ", ioffao(ish)
          Do jSh = 1, nSkal
            nBasJ = iSD(2,jSh)*iSD(3,jSh)
            Do iBas0 = 1, nBasI
              iBas = iOffAO(iSh) + iBas0
              Do jBas0 = 1, nBasJ
                jBas = iOffAO(jSh) + jBas0
                Do iOcc = 1, nOcc
                  Do jOcc = 1, nOcc
                    loc1 = jOcc-1 + (jBas-1)*nOcc
     *                   + (iOcc-1)*nOcc*nBasT
     *                   + (iBas-1)*nOcc*nBasT*nOcc
                    loc2 = iOcc-1 + (iBas-1)*nOcc
     *                   + (jOcc-1)*nOcc*nBasT
     *                   + (jBas-1)*nOcc*nBasT*nOcc
                    T_hbf(iOcc,jOcc,iBas0,jBas0)
     *                = T2AO(loc1+1)+T2AO(loc2+1)
C       if (nocc.eq.10.and.(iocc.le.2.or.jocc.le.2)) then
C                   T_hbf(iOcc,jOcc,iBas0,jBas0)=0.0d+00
C       end if
                  End Do
                End Do
                iRec = iBas+nBasT*(jBas-1)
        if (ifmscoup.and.jstate.ne.1) then
          read (lugamma,rec=irec) (work(ipwrk2+i-1),i=1,nocc*nocc)
          call daxpy_(nocc*nocc,1.0d+00,work(ipwrk2),1,
     *                                t_hbf(1,1,ibas0,jbas0),1)
        end if
                Write (LuGamma,Rec=iRec)
     *            (T_hbf(i,1,iBas0,jBas0),i=1,nOcc*nOcc)
              End Do
            End Do
C           write(6,*) "nbasi,nbasj = ", nbasi,nbasj
C           Call dDaFile(LuGamma,1,T_hbf,nOcc*nOcc*nBasI*nBasJ,iDisk)
          End Do
        End Do
C       call abend
C       Call DaClos(LuGamma)
        Close (LuGamma)
C       write(6,*) "end"
        Deallocate (iOffAO)
        Deallocate (T_hbf)
C       write(6,*) "end1"
C       call abend
C       Call DaClos(LuGamma)
        Call Free_iSD()
C       write(6,*) "end2"
        call clssew
      End If
C
      !! 5) L_{ai} = sum_{mu} C_{mu a} L_{mu i}
C     CALL DGEMM_('T','N',nSsh(iSym),nOcc,nBasT,
C    *            1.0D+00,Work(LCMOPT2+nBasT*nOcc),nBasT,
C    *                    Work(ipvLag),nBasT,
C    *            1.0D+00,Work(ipOLAG+nOCC),nOrb(iSymA))
C     write(6,*) "olag before vvvo"
C     call sqprt(work(ipolag),nbast)
      nOrbA = nFro(iSymA)+nIsh(iSymA)+nAsh(iSymA)+nSsh(iSymA)
      If (DoCholesky) nOcc = nOrbA-nFro(iSymA)
C     write(6,*) "vLag"
C     call sqprt(work(ipvLag),nbast)
      CALL DGEMM_('T','N',nOrbA,nOcc,nBasT,
     *            1.0D+00,Work(LCMOPT2),nBasT,
     *                    Work(ipvLag),nBasT,
     *            1.0D+00,Work(ipOLAG+nOrbA*nFro(iSymA)),nOrbA)
C     write(6,*) "olag after vvvo"
C     call sqprt(work(ipolag),nbast)
C
      Call GetMem('WRK1','Free','Real',ipWRK1,nBasT*nBasT)
      Call GetMem('WRK2','Free','Real',ipWRK2,nBasT*nBasT)
C
      END SUBROUTINE OLagVVVO

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
      !! ftwo_drv.f
      SUBROUTINE VVVO_Drv(nSym,nBas,nAsh,nFro,nSkipX,
     *                    iSym,iSymI,iSymJ,iSymK,iSymL,
     &                    T2AO,vLag,nOcc,nBasT,nTot2,
     &                    nBMX,CMO,DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,
     *                    DIA,DI,FIFA,FIMO)


      Implicit Real*8 (A-H,O-Z)
#include "rasdim.fh"
#include "WrkSpc.fh"
#include "real.fh"

      Integer nBas(8), nAsh(8), nSkipX(8), nfro(8)
      Dimension CMO(*),T2AO(*),vLag(*)
      Dimension DPT2AO(*),DPT2CAO(*),FPT2AO(*),FPT2CAO(*)
      Dimension DIA(*),DI(*),FIFA(*),FIMO(*)
      Logical DoCholesky,REORD,DECO
      Integer ALGO
      COMMON /CHORAS/ REORD,DECO,ALGO


      Call DecideOncholesky(DoCholesky)

C     IF (DoCholesky.and.ALGO.eq.2)THEN
*
* Building of the Fock matrix directly from Cholesky
* vectors
*
C        Call GetMem('LWFSQ','Allo','Real',LWFSQ,nTot2)
C        call dcopy_(nTot2,[Zero],0,Work(LWFSQ),1)

C        Call Allocate_Work(ipTemp,nTot1)
C        Call FZero(Work(ipTemp),nTot1)
C
C        CALL CHORAS_DRV(nSym,nBas,nAsh,D1A,DI,Work(ipTemp),
C    &                   ExFac,LWFSQ,CMO)

C        Call DaXpY_(nTot1,One,Work(ipTemp),1,FA,1)
*
C        Call Free_Work(ipTemp)
C        Call GetMem('LWFSQ','Free','Real',LWFSQ,nTot2)

C     ELSE

*
* Standard building of the Fock matrix from Two-el integrals
*
C        write(6,*) "calling drv2"
         Call VVVO_Drv2(nSym,nBas,nAsh,nFro,nSkipX,
     *                  iSym,iSymI,iSymJ,iSymK,iSymL,
     &                  T2AO,vLag,nOcc,nBasT,
     &                  nTot2,nBMX,CMO,DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,
     *                  DIA,DI,FIFA,FIMO)

C     ENDIF


      RETURN
C
      END SUBROUTINE VVVO_Drv
C
C-----------------------------------------------------------------------
C
      !! focktwo_drv.f
      Subroutine VVVO_Drv2(nSym,nBas,nAux,nFro,Keep,
     *                     iSym,iSymI,iSymJ,iSymK,iSymL,
     &                     T2AO,vLag,nOcc,nBasT,
     &                     nBSQT,nBMX,CMO,DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,
     *                     DIA,DI,FIFA,FIMO)
C
      USE CHOVEC_IO
C
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 T2AO(*),vLag(*),CMO(*)
      Dimension DPT2AO(*),DPT2CAO(*),FPT2AO(*),FPT2CAO(*)
      DImension DIA(*),DI(*),FIFA(*),FIMO(*)
      Integer nBas(8), nAux(8), Keep(8), nfro(8)
      Logical DoCholesky,GenInt
      Integer ALGO
      Logical REORD,DECO
C     Real*8 CMO_DUMMY(1)

      Common /CHORAS / REORD,DECO,ALGO
*
* nAux is the number of occupied orbitals
      GenInt=.false.
      DoCholesky=.false.
      if(ALGO.eq.0) GenInt=.true. !use GenInt to regenerate integrals

      Call DecideOnCholesky(DoCholesky)

      Call GetMem('LWFSQ','Allo','Real',LWFSQ,NBSQT)
      call dcopy_(NBSQT,[Zero],0,Work(LWFSQ),1)

C     if((.not.DoCholesky).or.(GenInt)) then
      Call GetMem('LW2','Allo','Real',LW2,NBMX*NBMX)
C     end if
*
C     Call Allocate_Work(ipTemp,nFlt)
C     Call FZero(Work(ipTemp),nFlt)
*
C     write(6,*) "lbuf = ", lbuf
      Call GetMem('LW1','MAX','Real',LW1,LBUF)
      If (DoCholesky) LBUF = NBMX*NBMX+1
C     write(6,*) "lbuf = ", lbuf
*
* Standard building of the Fock matrix from Two-el integrals
*

C     IF (.not.DoCholesky) THEN
         Call GetMem('LW1','Allo','Real',LW1,LBUF)

      If (LBUF.LT.1+NBMX**2) Then
         WRITE(6,*)' FockTwo_Drv Error: Too little memory remains for'
     &     //' the call to FOCKTWO.'
         WRITE(6,*)' Largest allocatable array size LBUF=',LBUF
         WRITE(6,*)' Max nr of bf in any symmetry,  NBMX=',NBMX
         WRITE(6,*)' Required minimum size     1+NBMX**2=',1+NBMX**2
         WRITE(6,*)'    (All in Real*8-size words)'
         Call  ABEND()
      End If
*
C        write(6,*) "calling vvvox"
      Call GetMem('LWRK','Allo','Real',LWRK,nBasT*nBasT)
      IF (DoCholesky) Then
C       write(6,*) "calling vvvox2"
        Call VVVOX2(nAux,Keep,iSym,iSymI,iSymJ,iSymK,iSymL,
     *              vLag,CMO,Work(LWRK),
     *              DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,
     *              FIFA,FIMO)
C       write(6,*) "finished vvvox2"
      Else
        Call VVVOX(nSym,nBas,nFro,Keep,
     *             iSymI,iSymJ,iSymK,iSymL,
     *             T2AO,vLag,CMO,nOcc,nBasT,
     *             LBUF,Work(LW1),Work(LW2),Work(LWRK),
     *             DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,
     *             DIA,DI,FIFA,FIMO)
      End If
C     write(6,*) "aaa"
      !! vLag must be transposed
      !! In VVVOX(2) subroutines, vLag(p,mu) is constructed.
      call dcopy_(nbast*nbast,vlag,1,Work(LWRK),1)
      do i = 1, nbast
        do j = 1, nbast
          vlag(i+nbast*(j-1)) = Work(LWRK+j-1+nbast*(i-1))
        end do
      end do
      Call GetMem('LWRK','Free','Real',LWRK,nBasT*nBasT)

C     ENDIF

*
* Building of the Fock matrix directly from Cholesky vectors
*
      IF (DoCholesky .and. .not.GenInt) THEN

*
* CMO_DUMMY is required call argument of choras_drv:
* (Not used, see logical flags in choras_drv)
*      SUBROUTINE CHORAS_DRV(nSym,nBas,nOcc,DSQ,DLT,FLT,
*     &                      ExFac,LWFSQ,CMO)
C      CALL CHOras_drv(nSym,nBas,nAux,DSQ,DLT,
C    &                 Work(ipTemp),ExFac,LWFSQ,CMO_DUMMY)
*
      ENDIF
*


C     Call DaXpY_(nFlt,One,Work(ipTemp),1,FLT,1)
*
C     Call Free_Work(ipTemp)

C     if(.not.DoCholesky)then
      Call GetMem('LW1','Free','Real',LW1,LBUF)
      Call GetMem('LW2','Free','Real',LW2,NBMX*NBMX)
C     endif

      Call GetMem('LWFSQ','Free','Real',LWFSQ,NBSQT)
*
      Return
C
      End SUBROUTINE VVVO_Drv2
C
C-----------------------------------------------------------------------
C
      !! focktwo.f
      SUBROUTINE VVVOX(NSYM,NBAS,NFRO,KEEP,
     *                 iSymI,iSymJ,iSymK,iSymL,
     &                 T2AO,vLag,CMO,nOcc,nBasT,LBUF,X1,X2,WRK,
     *                 DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,
     *                 DIA,DI,FIFA,FIMO)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      Real*8 T2AO(nOcc,nBasT,nOcc,nBasT),vLag(nBasT,*),CMO(nBasT,*),
     *       X1(*),X2(*),WRK(*)
      Dimension DPT2AO(*),DPT2CAO(*),FPT2AO(*),FPT2CAO(*)
      Dimension DIA(*),DI(*),FIFA(*),FIMO(*)
      Integer ISTLT(8),ISTSQ(8),KEEP(8),NBAS(8),NFRO(8)
C
C     write(6,*) "start vvvox"
C     MUL(I,J)=IEOR(I-1,J-1)+1
      ISTSQ(1)=0
      ISTLT(1)=0
      Do iSym = 2, nSym
        nB  = nBas(iSym-1)
        nB2 = nB*nB
        nB3 = (nB2+nB)/2
        ISTSQ(iSym) = ISTSQ(iSym-1) + nB2
        ISTLT(iSym) = ISTLT(iSym-1) + nB3
      End Do
C     write(6,*) "a"
      nFroT = 0
      Do iSym = 1, nSym
        nFroT = nFroT + nFro(iSym)
      End Do
C
      nBasI  = nBas(iSymI)
      KEEPI  = KEEP(iSymI)
      ! nAuxI  = nAux(iSymI)
      nBasJ  = nBas(iSymJ)
      KEEPJ  = KEEP(iSymJ)
      ! nAuxJ  = nAux(iSymJ)
      iSymIJ = 1+iEor(iSymI-1,iSymJ-1)
      nBasIJ = nBasI*nBasJ
      If (iSymI.EQ.iSymJ) nBasIJ = (nBasI*(nBasI+1))/2
      If (nBasIJ.eq.0) Return
C     write(6,*) "b"

      nBasK  = nBas(iSymK)
      KEEPK  = KEEP(iSymK)
      ! nAuxK  = nAux(iSymK)
      iSMax  = iSymK
      If (iSymK.EQ.iSymI) iSMax = iSymJ
      iSymL  = 1+iEor(iSymIJ-1,iSymK-1)
      IF (iSymL.GT.iSMax) Return !! should not
      nBasL  = nBas(iSymL)
      KEEPL  = KEEP(iSymL)
      ! nAuxL  = nAux(iSymL)
      nBasKL = nBasK*nBasL
      IF (iSymK.EQ.iSymL) nBasKL = (nBasK*(nBasK+1))/2
      If (nBasKL.eq.0) Return
C     write(6,*) "c"
C
      ! INTEGRAL BLOCK EXCLUDED BY SETTING KEEP PARAMETERS?
      IF (KEEPI+KEEPJ+KEEPK+KEEPL.NE.0) Return
C     write(6,*) "d"
      !! This will not work when the number of the inactive orbital is 0
C     IF (nAuxI+nAuxJ+nAuxK+nAuxL.EQ.0) Return ! frozen orbitals
C     write(6,*) "e"
C     write(6,*) "nbasij = ", nbasij, 6*13
C     write(6,*) "keep=",keepi,keepj,keepk,keepl
C     write(6,*) "CMO"
C     call sqprt(cmo,nbast)
C
      !! (ij|kl)
C     write(6,*) "doing actual calculation"
      If (iSymI.EQ.iSymJ .AND. iSymI.EQ.iSymK) Then
C       write(6,*) "aa"
C       write(6,*) "nocc,nbast = ", nocc,nbast
        IOPT=1
        LPQ=0
        IPQ=0
        NPQ=0
C       write(6,*) "nbasi = ", nbasi
        DO IP = 1, nBasI
          DO JQ = 1, IP
C           write(6,*) "ip,jq = ", ip,jq
            IPQ=IPQ+1
            LPQ=LPQ+1
            IF (IPQ.GT.NPQ) THEN
              CALL RDORD(IRC,IOPT,iSymI,iSymJ,iSymK,iSymL,X1,LBUF,NPQ)
              IF(IRC.GT.1) GOTO 999
                IOPT=2
                IPQ=1
            ENDIF
C           write(6,*) "do"
            ISX=(IPQ-1)*nBasKL+1
            ISF=ISTLT(iSymI)+LPQ
            ISD=ISTLT(iSymI)+1
            !! coulomb
C           TEMP=DDOT_(KLB,X1(ISX),1,DLT(ISD),1)
C           FLT(ISF)=FLT(ISF)+TEMP
            !! exchange
            !! f_{ik} = (ij|kl)*D_{jl}
            !! read (ij|kl) for given ij as (kl)
            !! first  dgemv: (kl)*D_{l,i} -> f_{ki}
            !! second dgemv: (kl)*D_{l,i} -> f_{ki}

            !! For (mu nu), read (mu nu|rho sigma) as I(rho sigma)
            !! (mu j|rho sigma) = sum_{nu} C_{nu j} (mu nu|rho sigma)

            !! L_{mu i} = T_{ij}^{mu nu} * (mu rho | j sigma)
            !!          = T_{ij}^{rho sigma} * C_{nu j}
            !!                 * (mu rho | nu sigma)

            !! IP = mu, JQ = rho
            !! X2 = (nu sigma) --> X2' = C_{nu j}*(nu sigma) = (j sigma)
            !! T_{ij}^{rho sigma} * (j sigma)
            !!   -> U_{i rho} for (mu, rho) pairs

            CALL SQUARE (X1(ISX),X2(1),1,nBasK,nBasL)
C     write(6,*) "ip,jq= ",ip,jq
C     write(6,*) "integral"
C     call sqprt(x2,nbask)
C           If (DoCholesky) Then
C           Else
            !! (mu(ip) rho(jq) | nu sigma) -> (mu(ip) rho(jq) | j sigma)
            call dgemm_('T','N',nOcc,nBasT,nBasT,
     *                  1.0D+00,CMO,nBasT,X2,nBasT,
     *                  0.0D+00,WRK,nOcc)
C           call dgemm_('T','T',nOcc,nBasT,nBasT,
C    *                  1.0D+00,CMO,nBasT,X2,nBasT,
C    *                  1.0D+00,WRK,nOcc)
C           write(6,*) "dgemm finished"
            !! wrk(j,sigma) for the given mu(ip), mu(jq)
            !! T2AO(j,sigma,i,rho) = T_{ij}^{rho sigma}
            !! rather than L_{mu i}, L_{i mu} is computed
            !! L_{i mu} = wrk(j,sigma)*(T2AO(j,sigma,i,rho)
            call dgemv_('t',nOcc*nBasT,nOcc,
     *                  1.0d+00,t2ao(1,1,1,jq),nOcc*nBasT,wrk,1,
     *                  1.0d+00,vlag(1,ip),1)
            if (ip.ne.jq) then
              call dgemv_('t',nOcc*nBasT,nOcc,
     *                    1.0d+00,t2ao(1,1,1,ip),nOcc*nBasT,wrk,1,
     *                    1.0d+00,vlag(1,jq),1)
C           call dgemm_('T','T',nOcc,nBasT,nBasT,
C    *                  1.0D+00,CMO,nBasT,X2,nBasT,
C    *                  1.0D+00,WRK,nOcc)
C           call dgemv_('n',nOcc*nBasT,nOcc,
C    *                  1.0d+00,t2ao(1,ip,1,1),nOcc*nBasT,wrk,1,
C    *                  1.0d+00,vlag(1,jq),1)
            end if
C           End If
C
            !! DPT2AO -> FPT2AO transformation
            !! FPT2 = G(DPT2)
            Call FDGTRF(DPT2AO,FPT2AO)
            Call FDGTRF(DPT2CAO,FPT2CAO)
            If (nFroT.ne.0) Then
              Call FDGTRF(DIA,FIFA)
              Call FDGTRF(DI,FIMO)
            End If
            !! Coulomb
!           Val = DDot_(nBasK*nBasL,X2,1,DPT2AO(ISTSQ(iSymI)+1),1)
!           iSF = ISTSQ(iSYmI) + iP + nBasI*(jQ-1)
!           FPT2AO(iSF) = FPT2AO(iSF) + Val
!           !! Exchange
!           iSF = ISTSQ(iSymI) + (jQ-1)*nBasJ+1
!           iSD = ISTSQ(iSymI) + (iP-1)*nBasI+1
!           CALL DGEMV_('N',nBasK,nBasL,
!    *                  -0.5D+00,X2,nBasK,DPT2AO(iSD),1,
!    *                  1.0D+00,FPT2AO(iSF),1)
!           If (iP.ne.jQ) Then
!             iSF = ISTSQ(iSymI) + jQ + nBasJ*(iP-1)
!             FPT2AO(iSF) = FPT2AO(iSF) + Val
!             iSF = ISTSQ(iSymI) + (iP-1)*nBasI+1
!             iSD = ISTSQ(iSymI) + (jQ-1)*nBasJ+1
!             CALL DGEMV_('N',nBasK,nBasL,
!    *                    -0.5D+00,X2,nBasK,DPT2AO(iSD),1,
!    *                    1.0D+00,FPT2AO(iSF),1)
!           End If
!
!           !! DPT2CAO -> FPT2CAO transformation
!           !! Coulomb
!           Val = DDot_(nBasK*nBasL,X2,1,DPT2CAO(ISTSQ(iSymI)+1),1)
!           iSF = ISTSQ(iSYmI) + iP + nBasI*(jQ-1)
!           FPT2CAO(iSF) = FPT2CAO(iSF) + Val
!           !! Exchange
!           iSF = ISTSQ(iSymI) + (jQ-1)*nBasJ+1
!           iSD = ISTSQ(iSymI) + (iP-1)*nBasI+1
!           CALL DGEMV_('N',nBasK,nBasL,
!    *                  -0.5D+00,X2,nBasK,DPT2CAO(iSD),1,
!    *                  1.0D+00,FPT2CAO(iSF),1)
!           If (iP.ne.jQ) Then
!             iSF = ISTSQ(iSymI) + jQ + nBasJ*(iP-1)
!             FPT2CAO(iSF) = FPT2CAO(iSF) + Val
!             iSF = ISTSQ(iSymI) + (iP-1)*nBasI+1
!             iSD = ISTSQ(iSymI) + (jQ-1)*nBasJ+1
!             CALL DGEMV_('N',nBasK,nBasL,
!    *                    -0.5D+00,X2,nBasK,DPT2CAO(iSD),1,
!    *                    1.0D+00,FPT2CAO(iSF),1)
!           End If
          End Do
        End Do
      End If

!     ELSE IF ( IS.EQ.JS .AND. IS.NE.KS ) THEN
!       ! CASE 2: Integrals are of symmetry type (II/JJ)
!       ! Coulomb terms need to be accumulated only
!     ELSE IF ( IS.EQ.KS .AND. JS.EQ.LS ) THEN
!       ! CASE 3: Integrals are of symmetry type (IJ/IJ)
!       ! Exchange terms need to be accumulated only
!         IOPT=1
!         LPQ=0
!         IPQ=0
!         NPQ=0
!         DO IP=1,IB
!           DO JQ=1,JB
!             IPQ=IPQ+1
!             LPQ=LPQ+1
!             IF ( IPQ.GT.NPQ ) THEN
!               CALL RDORD(IRC,IOPT,IS,JS,KS,LS,X1,LBUF,NPQ)
!               IF(IRC.GT.1) GOTO 999
!               IOPT=2
!               IPQ=1
!             ENDIF
!             ISX=(IPQ-1)*KLB+1
!             IF ( NFI.NE.0 ) THEN
!               ISD=ISTSQ(IS)+(IP-1)*IB+1
!               ISF=ISTSQ(JS)+(JQ-1)*JB+1
!               CALL DGEMV_('N',LB,KB,(-0.5D0*ExFac),X1(ISX),LB,
!    &                       DSQ(ISD),1,1.0D0,FSQ(ISF),1)
!             ENDIF
!             IF ( NFJ.NE.0 ) THEN
!               ISD=ISTSQ(JS)+(JQ-1)*JB+1
!               ISF=ISTSQ(IS)+(IP-1)*IB+1
!               CALL DGEMV_('T',LB,KB,(-0.5D0*ExFac),X1(ISX),LB,
!    &                       DSQ(ISD),1,1.0D0,FSQ(ISF),1)
!             ENDIF
!           End Do
!         End Do
!     End If

      RETURN
 999  CONTINUE
      WRITE(6,*)' Error return code IRC=',IRC
      WRITE(6,*)' from RDORD call, in FTWOI.'
      CALL Abend
C
      Contains
C
      Subroutine FDGTRF(DD,FF)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension DD(*),FF(*)
C
      !! Coulomb
      Val = DDot_(nBasK*nBasL,X2,1,DD(ISTSQ(iSymI)+1),1)
      iSF = ISTSQ(iSYmI) + iP + nBasI*(jQ-1)
      FF(iSF) = FF(iSF) + Val
C
      !! Exchange
      iSF = ISTSQ(iSymI) + (jQ-1)*nBasJ+1
      iSD = ISTSQ(iSymI) + (iP-1)*nBasI+1
      CALL DGEMV_('N',nBasK,nBasL,
     *           -0.5D+00,X2,nBasK,DD(iSD),1,
     *            1.0D+00,FF(iSF),1)
      If (iP.ne.jQ) Then
        iSF = ISTSQ(iSymI) + jQ + nBasJ*(iP-1)
        FF(iSF) = FF(iSF) + Val
        iSF = ISTSQ(iSymI) + (iP-1)*nBasI+1
        iSD = ISTSQ(iSymI) + (jQ-1)*nBasJ+1
        CALL DGEMV_('N',nBasK,nBasL,
     *             -0.5D+00,X2,nBasK,DD(iSD),1,
     *              1.0D+00,FF(iSF),1)
      End If
C
      End Subroutine FDGTRF
C
      End Subroutine VVVOX
C
C-----------------------------------------------------------------------
C
      !! mabe taken from tracho3.f
      SUBROUTINE VVVOX2(nAux,KEEP,iSym,iSymI,iSymJ,iSymK,iSymL,
     *                  vLag,CMO,WRK,
     *                  DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,
     *                  FIFA,FIMO)

      use ChoVec_io
      use ChoSwp, only: InfVec
      use ChoArr, only: nDimRS

      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "warnings.h"
#include "caspt2.fh"
#include "eqsolv.fh"
#include "chocaspt2.fh"
#include "choglob.fh"
#include "WrkSpc.fh"
#include "caspt2_grad.fh"

      Real*8 vLag(nBasT,*),CMO(nBasT,*),WRK(nBasT,nBasT)
      Dimension DPT2AO(*),DPT2CAO(*),FPT2AO(*),FPT2CAO(*)
      Dimension FIFA(*),FIMO(*)
      Integer ISTLT(8),ISTSQ(8),nAux(8),KEEP(8),ipWRK(8)

      Integer iSkip(8)

!     INFVEC(I,J,K)=IWORK(ip_INFVEC-1+MAXVEC*N2*(K-1)+MAXVEC*(J-1)+I)

      Do jSym = 1, nSym
        iSkip(jSym) = 1
      End Do

      ISTSQ(1)=0
      ISTLT(1)=0
      nAuxT = 0
      Do jSym = 2, nSym
        nB  = nBas(jSym-1)
        nB2 = nB*nB
        nB3 = (nB2+nB)/2
        ISTSQ(jSym) = ISTSQ(jSym-1) + nB2
        ISTLT(jSym) = ISTLT(jSym-1) + nB3
        nAuxT = nAuxT + nAux(jSym)
      End Do

C     write(6,*) "sym=",isymi,isymj,isymk,isyml
      nBasI  = nBas(iSymI)
      KEEPI  = KEEP(iSymI)
      ! nAuxI  = nAux(iSymI)
      nBasJ  = nBas(iSymJ)
      KEEPJ  = KEEP(iSymJ)
      ! nAuxJ  = nAux(iSymJ)
      iSymIJ = 1+iEor(iSymI-1,iSymJ-1)
      nBasIJ = nBasI*nBasJ
      If (iSymI.EQ.iSymJ) nBasIJ = (nBasI*(nBasI+1))/2
C     write(6,*) "nbasij = ", nbasij
      If (nBasIJ.eq.0) Return

      nBasK  = nBas(iSymK)
      KEEPK  = KEEP(iSymK)
      ! nAuxK  = nAux(iSymK)
      iSMax  = iSymK
      If (iSymK.EQ.iSymI) iSMax = iSymJ
      iSymL  = 1+iEor(iSymIJ-1,iSymK-1)
C     write(6,*) "isyml,ismax = ", isyml,ismax
      IF (iSymL.GT.iSMax) Return !! should not
      nBasL  = nBas(iSymL)
      KEEPL  = KEEP(iSymL)
      ! nAuxL  = nAux(iSymL)
      nBasKL = nBasK*nBasL
      IF (iSymK.EQ.iSymL) nBasKL = (nBasK*(nBasK+1))/2
C     write(6,*) "nbaskl = ", nbaskl
      If (nBasKL.eq.0) Return
C
C     write(6,*) "keep=",keepi,keepj,keepk,keepl
      ! INTEGRAL BLOCK EXCLUDED BY SETTING KEEP PARAMETERS?
      IF (KEEPI+KEEPJ+KEEPK+KEEPL.NE.0) Return
C     write(6,*) "nAux=",nAuxi,nAuxj,nAuxk,nAuxl
      !! This will not work when the number of the inactive orbital is 0
C     IF (nAuxI+nAuxJ+nAuxK+nAuxL.EQ.0) Return ! frozen orbitals
C
      jSym = iSymJ
      ! kSym = iSymK
      ! lSym = iSymL
      nIshI = nIsh(iSym)
      ! nIshJ = nIsh(jSym)
      ! nIshK = nIsh(kSym)
      ! nIshL = nIsh(lSym)
      nAshI = nAsh(iSym)
      ! nAshJ = nAsh(jSym)
      ! nAshK = nAsh(kSym)
      ! nAshL = nAsh(lSym)
      nSshI = nSsh(iSym)
      ! nSshJ = nSsh(jSym)
      ! nSshK = nSsh(kSym)
      ! nSshL = nSsh(lSym)
      nOrbI = nIshI+nAshI+nSshI
C
C     write(6,*) "nchspc = ", nchspc
      CALL GETMEM('CHSPC','ALLO','REAL',IP_CHSPC,NCHSPC)
      CALL GETMEM('HTSPC','ALLO','REAL',IP_HTSPC,NHTSPC)
      CALL GETMEM('HTVEC','ALLO','REAL',ipHTVec,nBasT*nBasT)
      CALL GETMEM('WRK  ','ALLO','REAL',ipWRK(iSym),nBasT*nBasT)
C
      IBATCH_TOT=NBTCHES(iSym)

      IF(NUMCHO_PT2(iSym).EQ.0) Return

      ! ipnt=ip_InfVec+MaxVec_PT2*(1+InfVec_N2_PT2*(iSym-1))
      ! JRED1=iWork(ipnt)
      ! JRED2=iWork(ipnt-1+NumCho_PT2(iSym))
      JRED1=InfVec(1,2,jSym)
      JRED2=InfVec(NumCho_PT2(jSym),2,jSym)
C     write(6,*) "jred1,jred2 = ", jred1,jred2

* Loop over JRED
      DO JRED=JRED1,JRED2

        CALL Cho_X_nVecRS(JRED,iSym,JSTART,NVECS_RED)
        IF(NVECS_RED.EQ.0) Cycle

        ILOC=3
        CALL CHO_X_SETRED(IRC,ILOC,JRED)
* For a reduced set, the structure is known, including
* the mapping between reduced index and basis set pairs.
* The reduced set is divided into suitable batches.
* First vector is JSTART. Nr of vectors in r.s. is NVECS_RED.
        ! JEND=JSTART+NVECS_RED-1

* Determine batch length for this reduced set.
* Make sure to use the same formula as in the creation of disk
* address tables, etc, above:
        NBATCH=1+(NVECS_RED-1)/MXNVC

* Loop over IBATCH
        JV1=JSTART
        DO IBATCH=1,NBATCH
C         write(6,*) "ibatch,nbatch = ", ibatch,nbatch
          IBATCH_TOT=IBATCH_TOT+1

          JNUM=NVLOC_CHOBATCH(IBATCH_TOT)
          JV2=JV1+JNUM-1
C
C         ----- Construct orbital Lagrangian -----
C
          !! Work(ip_CHSPC) :: (mu nu|P)
          !! Work(ip_HTSPC) :: ( q nu|P)
          !! Bra and Ket    :: (ia|P) = T_{ij}^{ab}*(jb|P)
          !! In 1), HT_{i mu,P} = C_{mu a}*(ia|P)
          !!        (ia|P) read from disk
          !! In 3), L_{i mu} = HT_{i nu,P} * (mu nu|P)
          !! Then, L_{pi} = C_{mu p} * L_{i mu}^T
          !! Transpose is done in VVVO_Drv2,
          !! and C_{mu p} is in OLagVVVO
C
          !! 1) Half back-transformation of Bra and Ket density
          !! Read the 3c-2e pseudo-density (in MO), and half transform
          CALL VVVOTRA_RI(CMO,WORK(IP_CHSPC),WORK(IP_HTSPC),
     *                    JNUM,IBATCH_TOT,IBATCH_TOT,nOrbI)
C
          !! 2) read AO Cholesky vectors,
          !!    then, (strange) reduced form -> squared AO (mu nu|iVec)
          JREDC=JRED
* Read a batch of reduced vectors
          CALL CHO_VECRD(WORK(IP_CHSPC),NCHSPC,JV1,JV2,iSym,
     &                            NUMV,JREDC,MUSED)
          IF(NUMV.ne.JNUM) THEN
            write(6,*)' Rats! CHO_VECRD was called, assuming it to'
            write(6,*)' read JNUM vectors. Instead it returned NUMV'
            write(6,*)' vectors: JNUM, NUMV=',JNUM,NUMV
            write(6,*)' Back to the drawing board?'
            CALL QUIT(_RC_INTERNAL_ERROR_)
          END IF
          IF(JREDC.NE.JRED) THEN
            write(6,*)' Rats! It was assumed that the Cholesky vectors'
            write(6,*)' in HALFTRNSF all belonged to a given reduced'
            write(6,*)' set, but they don''t!'
            write(6,*)' JRED, JREDC:',JRED,JREDC
            write(6,*)' Back to the drawing board?'
            write(6,*)' Let the program continue and see what happens.'
          END IF
C
          !! (strange) reduced form -> squared AO (mu nu|iVec)
          !! is it possible to avoid this transformation?
      ! choptr.fh
          Call R2FIP(Work(ip_CHSPC),Work(ipWRK(iSym)),ipWRK(iSym),NUMV,
     *               size(nDimRS),infVec,nDimRS,
     *               nBasT,nSym,iSym,iSkip,irc,JREDC)
C
C           ----- Fock-like transformations (if needed) -----
C
          If (nFroT.eq.0) Then
            Do iVec = 1, NUMV
              Call FDGTRF(Work(ip_CHSPC+nBasT**2*(iVec-1)),
     *                    DPT2AO,FPT2AO)
              Call FDGTRF(Work(ip_CHSPC+nBasT**2*(iVec-1)),
     *                    DPT2CAO,FPT2CAO)
            End Do
          End If
C
          !! 3) Contract with Cholesky vectors
          Call DGemm_('N','T',nOrbI,nBasI,nBasI*JNUM,
     *                2.0D+00,Work(ip_HTSPC),nOrbI,
     *                        Work(ip_CHSPC),nBasI,
     *                1.0D+00,vLag,nBasI)
C
          !! 4) Construct the 3c-2e pseudo-density in AO
          !! D_{p nu} -> D_{mu nu}
          !! i.e., construct B_PT2, used in ALASKA
          Call DGemm_('N','N',nBasI,nBasI*JNUM,nOrbI,
     *                1.0D+00,CMO(1,1),nBasI,
     *                        Work(ip_HTSPC),nOrbI,
     *                0.0D+00,Work(ip_CHSPC),nBasI)
C
          !! 5) Save the 3c-2e pseudo-density in the disk
          !! it may be replaced with ddafile
          Do iVec = 1, NUMV
            If (IFMSCOUP.and.jState.ne.1) Then
              Read (LuGamma,Rec=iVec+JV1-1)
     *          (Work(ipHTVec+i-1),i=1,nBasI**2)
              Call DaXpY_(nBasI**2,1.0D+00,
     *                    Work(ip_CHSPC+nBasI**2*(iVec-1)),1,
     *                    Work(ipHTVec),1)
              Write (LuGamma,Rec=iVec+JV1-1)
     *          Work(ipHTVec:ipHTVec+nBasI**2-1)
            Else
             if (jState.eq.iRlxRoot .or. nStLag.gt.1) then
              Write (LuGamma,Rec=iVec+JV1-1)
     *        Work(ip_CHSPC+nBasI**2*(iVec-1):ip_CHSPC+nBasI**2*iVec-1)
             end if
            End If
          End Do
C
          JV1=JV1+JNUM
        End Do
      End Do
C
      CALL GETMEM('CHSPC','FREE','REAL',IP_CHSPC,NCHSPC)
      CALL GETMEM('HTSPC','FREE','REAL',IP_HTSPC,NHTSPC)
      CALL GETMEM('HTVEC','FREE','REAL',ipHTVec,nBasT*nBasT)
      CALL GETMEM('WRK  ','FREE','REAL',ipWRK(iSym),nBasT*nBasT)
C
      !! Have to (?) symmetrize Fock-transformed matrices
      If (nFroT.eq.0) Then
      Do i = 1, nBasI
        Do j = 1, i-1
          tmp = (FPT2AO(i+nBasI*(j-1))+FPT2AO(j+nBasI*(i-1)))*0.5d+00
          FPT2AO(i+nBasI*(j-1)) = Tmp
          FPT2AO(j+nBasI*(i-1)) = Tmp
          tmp = (FPT2CAO(i+nBasI*(j-1))+FPT2CAO(j+nBasI*(i-1)))*0.5d+00
          FPT2CAO(i+nBasI*(j-1)) = Tmp
          FPT2CAO(j+nBasI*(i-1)) = Tmp
          If (nFroT.ne.0) Then
            tmp = (FIFA(i+nBasI*(j-1))+FIFA(j+nBasI*(i-1)))*0.5d+00
            FIFA(i+nBasI*(j-1)) = Tmp
            FIFA(j+nBasI*(i-1)) = Tmp
            tmp = (FIMO(i+nBasI*(j-1))+FIMO(j+nBasI*(i-1)))*0.5d+00
            FIMO(i+nBasI*(j-1)) = Tmp
            FIMO(j+nBasI*(i-1)) = Tmp
          End If
        End Do
      End Do
      End If
C
      Return
C
      Contains
C
      Subroutine FDGTRF(ChoVec,DD,FF)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension ChoVec(*),DD(*),FF(*)
C
      !! Coulomb
      Scal = DDot_(nBasI**2,DD,1,ChoVec,1)
      Call DaXpY_(nBasI**2,Scal,ChoVec,1,FF,1)
C
      !! Exchange
      Call DGEMM_('T','N',nBasI,nBasI,nBasI,
     *            1.0D+00,ChoVec,nBasI,DD,nBasI,
     *            0.0D+00,WRK,nBasI)
      Call DGEMM_('T','T',nBasI,nBasI,nBasI,
     *           -0.5D+00,ChoVec,nBasI,WRK,nBasI,
     *            1.0D+00,FF,nBasI)
C
      End Subroutine FDGTRF
C
      Subroutine VVVOTRA_RI(CMO,CHSPC,HTSPC,NVEC,IBSTA,IBEND,nOrbI)
C
      Implicit Real*8(A-H,O-Z)
C
      Integer Active, Inactive, Virtual
      Parameter (Inactive=1, Active=2, Virtual=3)
C
      !! CHSPC is used as a temporary array
      Dimension CMO(nBasI,nOrbI),CHSPC(*),HTSPC(nOrbI,nBasT,*)
C
      !! BraAI
      Call Cholesky_Vectors(2,Inactive,Active,JSYM,CHSPC,nBra,
     *                      IBSTA,IBEND)
      IPQ = nAshI*nIshI
      Do jVec = 1, NVEC
        ! a. AI -> mu I
        Call DGemm_('T','T',nIshI,nBasI,nAshI,
     *              1.0D+00,CHSPC(1+IPQ*(jVec-1)),nAshI,
     *                      CMO(1,1+nIshI),nBasI,
     *              0.0D+00,HTSPC(1,1,jVec),nOrbI)
        ! a. AI -> A mu
        Call DGemm_('N','T',nAshI,nBasI,nIshI,
     *              1.0D+00,CHSPC(1+IPQ*(jVec-1)),nAshI,
     *                      CMO(1,1),nBasI,
     *              0.0D+00,HTSPC(1+nIshI,1,jVec),nOrbI)
      End Do
C
      !! BraSI
      Call Cholesky_Vectors(2,Inactive,Virtual,JSYM,CHSPC,nBra,
     *                      IBSTA,IBEND)
      IPQ = nIshI*nSshI
      Do jVec = 1, NVEC
        ! b. SI -> mu I
        Call DGemm_('T','T',nIshI,nBasI,nSshI,
     *              1.0D+00,CHSPC(1+IPQ*(jVec-1)),nSshI,
     *                      CMO(1,1+nIshI+nAshI),nBasI,
     *              1.0D+00,HTSPC(1,1,jVec),nOrbI)
        ! a. SI -> S mu
        Call DGemm_('N','T',nSshI,nBasI,nIshI,
     *              1.0D+00,CHSPC(1+IPQ*(jVec-1)),nSshI,
     *                      CMO(1,1),nBasI,
     *              0.0D+00,HTSPC(1+nIshI+nAshI,1,jVec),nOrbI)
      End Do
C
      !! BraSA
      Call Cholesky_Vectors(2,Active,Virtual,JSYM,CHSPC,nBra,
     *                      IBSTA,IBEND)
      IPQ = nAshI*nSshI
      Do jVec = 1, NVEC
        ! d. SA -> mu A
        Call DGemm_('T','T',nAshI,nBasI,nSshI,
     *              1.0D+00,CHSPC(1+IPQ*(jVec-1)),nSshI,
     *                      CMO(1,1+nIshI+nAshI),nBasI,
     *              1.0D+00,HTSPC(1+nIshI,1,jVec),nOrbI)
        ! b. SA -> S mu
        Call DGemm_('N','T',nSshI,nBasI,nAshI,
     *              1.0D+00,CHSPC(1+IPQ*(jVec-1)),nSshI,
     *                      CMO(1,1+nIshI),nBasI,
     *              1.0D+00,HTSPC(1+nIshI+nAshI,1,jVec),nOrbI)
      End Do
C
      !! BraAA
      Call Cholesky_Vectors(2,Active,Active,JSYM,CHSPC,nBra,
     *                      IBSTA,IBEND)
      IPQ = nAshI*nAshI
      Do jVec = 1, NVEC
        ! b. AA -> mu A
        Call DGemm_('T','T',nAshI,nBasI,nAshI,
     *              1.0D+00,CHSPC(1+IPQ*(jVec-1)),nAshI,
     *                      CMO(1,1+nIshI),nBasI,
     *              1.0D+00,HTSPC(1+nIshI,1,jVec),nOrbI)
        ! c. AA -> A mu
        Call DGemm_('N','T',nAshI,nBasI,nAshI,
     *              1.0D+00,CHSPC(1+IPQ*(jVec-1)),nAshI,
     *                      CMO(1,1+nIshI),nBasI,
     *              1.0D+00,HTSPC(1+nIshI,1,jVec),nOrbI)
      End Do
C
      End Subroutine VVVOTRA_RI
C
      End Subroutine VVVOX2
C
C-----------------------------------------------------------------------
C
      subroutine getritrfinfo(nnbstr_,maxvec_,n2_)
      implicit real*8(a-h,o-z)
#include "cholesky.fh"
      dimension nnbstr_(8,3)
      maxvec_ = maxvec
      n2_     = infvec_n2
      nnbstr_ = nnbstr
      end subroutine getritrfinfo
C
C-----------------------------------------------------------------------
C
      Subroutine R2FIP(CHSPC,WRK,ipWRK,NUMV,l_NDIMRS,INFVEC,
     *                 nDimRS,nBasT,nSym0,iSym,iSkip,irc,JREDC)
C
      Implicit Real*8 (A-H,O-Z)
C
#include "cholesky.fh"
C
      Dimension CHSPC(nBasT**2,*),WRK(*),ipWRK(*)
      Dimension INFVEC(MAXVEC,INFVEC_N2,*),nDimRS(nSym0,*),iSkip(8)
C
C     Transform the reduced form to the full form in place
C
      kloc = 0
      Do iVec = 1, NUMV
        If (l_NDIMRS.LT.1) Then
          lscr  = NNBSTR(iSym,3)
        Else
          JREDL = INFVEC(iVec,2,iSym)
          lscr  = nDimRS(iSym,JREDL) !! JRED?
        End If
        kloc = kloc + lscr
      End Do
C
      ipVecL = 1 + kloc !! lscr*(JNUM-1)
      jloc = NUMV
      Do iVec = NUMV, 1, -1
        If (l_NDIMRS.LT.1) Then
          lscr  = NNBSTR(iSym,3)
        Else
          JREDL = INFVEC(iVec,2,iSym)
          lscr  = nDimRS(iSym,JREDL) !! JRED?
        End If
        ipVecL = ipVecL - lscr
        Call DCopy_(nBasT**2,[0.0D+00],0,WRK,1)
        Call Cho_ReOrdr(irc,CHSPC(ipVecL,1),lscr,1,
     *                  1,1,1,iSym,JREDC,2,ipWRK,
     *                  iSkip)
        Call DCopy_(nBasT**2,WRK,1,CHSPC(1,jloc),1)
        jloc = jloc-1
      End Do
C
      Return
C
      End Subroutine R2FIP
