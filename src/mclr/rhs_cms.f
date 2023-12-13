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
* Copyright (C) 2021, Jie J. Bao                                       *
************************************************************************
* ****************************************************************
* history:                                                       *
* Jie J. Bao, on Aug. 06, 2020, created this file.               *
* ****************************************************************
      subroutine RHS_CMS(Fock,CICSF)
      use stdalloc, only : mma_allocate, mma_deallocate
*#include "stdalloc.fh"
#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "real.fh"
#include "sa.fh"

******Input
******Output
      Real*8,DIMENSION(nDens2+6)::Fock
      Real*8,DIMENSION(nconf1*nroots)::CICSF
******Auxiliary Quantities
      Real*8,DIMENSION((nRoots+1)*nRoots/2,nnA,nnA)::GDMat
      Real*8,DIMENSION((nRoots+1)*nRoots/2,(nRoots+1)*nRoots/2)::W
      INTEGER,DIMENSION(ntBas,ntAsh,ntAsh,ntAsh)::IndPUVX
      INTEGER,DIMENSION(ntAsh,ntAsh,ntAsh,ntAsh)::IndTUVX

      Real*8,DIMENSION(:),Allocatable::PUVX,R,H,AXkzx,AXPzx,AXX,
     &bk,bP,bX,FMO1t,FMO2t,zX
******R:     rotation matrix from inter states to final states
*      Real*8,DIMENSION(:,:),Allocatable::GDMat
*      INTEGER,DIMENSION(:,:,:,:),Allocatable::IndPUVX,IndTUVX
      INTEGER NPUVX,NTri

******MEMORY ALLOCATION
      CALL mma_allocate(AXkzx,nDens2)
      CALL mma_allocate(AXPzx,NConf1*nRoots)
      CALL mma_allocate(AXX,((nRoots-1)*nRoots/2)**2)
      CALL mma_allocate(R,nRoots**2)
      CALL mma_allocate(H,nRoots**2)
      CALL mma_allocate(bk,nDens2)
      CALL mma_allocate(bP,nConf1*nRoots)
      CALL mma_allocate(bX,(nRoots-1)*nRoots/2)
      CALL mma_allocate(zX,(nRoots-1)*nRoots/2)
*      CALL mma_allocate(GDMat,(nRoots+1)*nRoots/2,nnA,nnA)
*      CALL mma_allocate(IndPUVX,ntBas,ntAsh,ntAsh,ntAsh)
*      CALL mma_allocate(IndTUVX,ntAsh,ntAsh,ntAsh,ntAsh)
      Call Get_PUVXLen(NPUVX)
      CALL mma_allocate(PUVX,NPUVX)
      CALL Get_Ntri(nTri)
      CALL mma_allocate(FMO1t,nRoots*nTri)
      NACPAR=(nnA+1)*nnA/2
      NAcPr2=(nacpar+1)*nacpar/2
      CALL mma_allocate(FMO2t,nRoots*nacpr2)
******MAIN COURSE
******First, read results printed in MCPDFT module
      CALL CMSRdMat(H,nRoots,nRoots,'ROT_HAM',7)
      CALL Get_DArray('MS_FINAL_ROT    ',R,nRoots**2)
      Call Read_PUVX(PUVX,NPUVX)
      CALL Get_Two_Ind(IndPUVX,IndTUVX)
      CALL GetPDFTFocks(FMO1t,FMO2t,nTri)
******Calculate six additional terms in CMS Lagrangian equaiton
      CALL CMSRHSGDMat(GDMat)
      CALL CalcW(W,GDMAt,PUVX,NPUVX,IndTUVX)

      CALL CalcAXX(AXX,W)

      CALL CalcbXbP(bX,bP,FMO1t,FMO2t,R,H,nTri)

      CALL SolveforzX(zX,AXX,bX)

      CALL CalcAXkzx(AXkzx,GDMat,PUVX,NPUVX,IndPUVX,zx)

      CALL CalcAXPzx(AXPzx,GDMat,PUVX,NPUVX,IndTUVX,W,zx)

      CALL Calcbk(bk,R,nTri,GDMat,zX)

      CALL SolveforRHS(Fock,CICSF,AXkzx,AXPzx,bk,bP)

******MEMORY DEALLOCATION
      CALL mma_deallocate(AXkzx)
      CALL mma_deallocate(AXPzx)
      CALL mma_deallocate(AXX)
      CALL mma_deallocate(R)
      CALL mma_deallocate(H)
      CALL mma_deallocate(bk)
      CALL mma_deallocate(bP)
      CALL mma_deallocate(bX)
      CALL mma_deallocate(zX)
      CALL mma_deallocate(FMO1t)
      CALL mma_deallocate(FMO2t)
      CALL mma_deallocate(PUVX)
      RETURN
      end subroutine
******************************************************
      Subroutine CMSRdMat(Mat,NRow,NCol,FileName,NameLen)

      INTEGER NRow,NCol,NameLen
      Real*8,DIMENSION(NRow*NCol)::Mat
      CHARACTER(Len=NameLen)::FileName
      INTEGER I,J,LU
      External IsFreeUnit

      LU=233
      LU=IsFreeUnit(LU)
      CALL Molcas_Open(LU,FileName)
      DO I=1,NRow
       read(LU,*)(Mat((I-1)*NCol+J),J=1,NCol)
      END DO
      CLOSE(LU)

      RETURN
      END Subroutine
******************************************************


******************************************************
      Subroutine Get_PUVXLen(NPUVX)
******Rewritten from mcpdft/alloc.f
      INTEGER NPUVX
      INTEGER iSp,iSq,iSr,iSs,nAq,iSpq,iSpqr,nAr,nAs,nOp
#include "stdalloc.fh"
#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "real.fh"
#include "sa.fh"
#include "crun_mclr.fh"

      NPUVX=0
      DO iSp=1,nSym
       nOp=NORB(iSp)
       Do iSq=1,nSym
        nAq=NASH(iSq)
        iSpq=IEOR(iSp-1,iSq-1)
        dO iSr=1,nSym
         iSpqr=IEOR(iSpq,iSr-1)+1
         nAr=NASH(iSr)
         do iSs=1,iSr
          IF(iSpqr.ne.iSs) GO TO 11
          nAs=NASH(iSs)
          nRS=nAr*nAs
          IF(iSs.eq.iSr) nRS=(nAr+nAr**2)/2
          NPUVX=NPUVX+nOp*nAq*nRS
11        CONTINUE
         end do
        eND dO
       End Do
      END DO

      RETURN
      End Subroutine
******************************************************

******************************************************
      Subroutine Read_PUVX(PUVX,NPUVX)
      INTEGER NPUVX
      Real*8,DIMENSION(NPUVX)::PUVX
      CALL Get_Darray('TwoEIntegral    ',PUVX,nPUVX)
      RETURN
      End Subroutine
******************************************************

      Subroutine Get_Two_Ind(Ind_PUVX,IndTUVX)
************************************************************************
*     Readapted from src/fock_util/get_tuvx.f
*     Return to an index in the PUVX array given
*     four MO indices.
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "Input.fh"
******Output
      INTEGER,DIMENSION(ntBas,ntAsh,ntAsh,ntAsh)::Ind_PUVX
      INTEGER,DIMENSION(ntAsh,ntAsh,ntAsh,ntAsh)::IndTUVX
******Auxiliaries
      Integer,DIMENSION(nSym):: off_Ash,off_PUVX,off_Orb

      iTri(i) = (i*i-i)/2

*      generate offsets

***** Initialization
      DO lOrb=1, ntAsh
       Do KOrb=1, ntAsh
        do JOrb=1, ntAsh
         do iOrb=1, ntAsh
          IndTUVX(iOrb,jOrb,kOrb,lOrb)=0
          Ind_PUVX(iOrb,jOrb,kOrb,lOrb)=0
         end do
         do iOrb=ntAsh+1,ntBas
          Ind_PUVX(iOrb,jOrb,kOrb,lOrb)=0
         end do
        end do
       End Do
      END DO

      iStack = 0
      Do iSym = 1,nSym
        off_Orb(iSym) = iStack
        iStack = iStack + nOrb(iSym)
      End Do

      iStack = 0
      Do iSym = 1,nSym
        off_Ash(iSym) = iStack
        iStack = iStack + nAsh(iSym)
      End Do

      iStack = 0
      Do iSym = 1,nSym
        off_PUVX(iSym) = iStack
        iOrb = nOrb(iSym)
        Do jSym = 1,nSym
          jAsh = nAsh(jSym)
          ijSym = 1 + ieor(iSym-1,jSym-1)
          Do kSym = 1,nSym
            kAsh = nAsh(kSym)
            Do lSym = 1,kSym
              lAsh = nAsh(lSym)
              klSym = 1 + ieor(kSym-1,lSym-1)
              If ( ijSym.eq.klSym) then
                kl_Orb_pairs = kAsh*lAsh
                If ( kSym.eq.lSym ) kl_Orb_pairs = (kAsh*kAsh+kAsh)/2
                iStack = iStack + iOrb*jAsh*kl_Orb_pairs
              End If
            End Do
          End Do
        End Do
      End Do

*     select integrals with all 4 indices active

      Do iSym = 1,nSym
        iOrb = nOrb(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        iPUVX = off_PUVX(iSym)
        Do jSym = 1,nSym
          jAsh = nAsh(jSym)
          ijSym = 1 + ieor(iSym-1,jSym-1)
          Do kSym = 1,nSym
            kAsh = nAsh(kSym)
            lSym = 1 + ieor(ijSym-1,kSym-1)
            lAsh = nAsh(lSym)

            If ( lSym.le.kSym .and.
     &           iAsh*jAsh*kAsh*lAsh.ne.0 ) then
              Do iV = 1,kAsh
                lMax = lAsh
                If ( kSym.eq.lSym ) lMax = iV
                Do iX = 1,lMax
                  Do iU = 1,jAsh
                    Do iP = 1,iOrb
                      iT = iP - iIsh
                      iPUVX=iPUVX+1
      io=iP+Off_Orb(Isym)
      jo=iU+Off_Ash(Jsym)
      ko=iV+Off_Ash(ksym)
      lo=iX+Off_Ash(lsym)
      Ind_PUVX(io,jo,ko,lo)=iPUVX
      Ind_PUVX(io,jo,lo,ko)=iPUVX
                      If ( iT.gt.0 .and. iT.le.iAsh ) then
                        iiT = iT + off_Ash(iSym)
                        iiU = iU + off_Ash(jSym)
                        If ( iiU.gt.iiT ) then
                          iiT = iU + off_Ash(jSym)
                          iiU = iT + off_Ash(iSym)
                        End If
*
                        iTU = iiU + iTri(iiT)
                        iiV = iV + off_Ash(kSym)
                        iiX = iX + off_Ash(lSym)
                        If ( iiX.gt.iiV ) then
                          iiV = iX + off_Ash(lSym)
                          iiX = iV + off_Ash(kSym)
                        End If
                        iVX = iiX + iTri(iiV)
                        If ( iVX.gt.iTU ) then
                          iTemp = iTU
                          iTU = iVX
                          iVX = iTemp
                        End If
                        IndTUVX(iiT,iiU,iiV,iiX)=iPUVX
                        IndTUVX(iiT,iiU,iiX,iiV)=iPUVX
                        IndTUVX(iiU,iiT,iiV,iiX)=iPUVX
                        IndTUVX(iiU,iiT,iiX,iiV)=iPUVX
                        IndTUVX(iiV,iiX,iiT,iiU)=iPUVX
                        IndTUVX(iiV,iiX,iiU,iiT)=iPUVX
                        IndTUVX(iiX,iiV,iiT,iiU)=iPUVX
                        IndTUVX(iiX,iiV,iiU,iiT)=iPUVX
                      End If
                    End Do
                  End Do
                End Do
              End Do
            End If
          End Do
        End Do
      End Do
      Return
      End
******************************************************


******************************************************
      Subroutine CMSRHSGDMat(GDMat)
      use ipPage, only: W
#include "Pointers.fh"
#include "Input.fh"
#include "Files_mclr.fh"
#include "stdalloc.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
*      Input
*      Output
       Real*8,DIMENSION(nRoots*(nRoots+1)/2,nnA,nnA)::GDMat
*      Auxiliary quantities
       Real*8,DIMENSION(:),Allocatable::GDArray
       Real*8,DIMENSION(n2dens)::rdum
       INTEGER I,J,IOrb,JOrb,NIJ
*      I:state index, "excited state" of state J when I .ne. J
*      IOrb:row index,    orbital index for state I
*      JOrb:column index, orbital index for state J
       INTEGER nConfL,nConfR,iL,iR
       Real*8,Allocatable::CIL(:),CIR(:)
*      (D^IJ)_pq = <I|E_pq|J>, setting I>=J
*       <I|E_pq|J>=<J|E_qp|I>
       Call mma_allocate(GDArray,n1dens)
       iL=state_sym
       iR=state_sym
       nConfR=Max(ncsf(iR),nint(xispsm(iR,1)))
       nConfL=Max(ncsf(iL),nint(xispsm(iL,1)))
       Call mma_allocate(CIR,nConfR)
       Call mma_allocate(CIL,nConfL)
       DO I=1,nRoots
        Call CSF2SD(W(ipCI)%Vec(1+(I-1)*ncsf(iL)),CIL,iL)
        Do J=1,I
        Call CSF2SD(W(ipCI)%Vec(1+(J-1)*ncsf(iR)),CIR,iR)
         Call Densi2(1,GDArray,rdum,CIL,CIR,0,0,0,n1dens,n2dens)
         NIJ=I*(I-1)/2+J
         do IOrb=1,nnA
          do JOrb=1,nnA
           GDMat(NIJ,IOrb,JOrb)=GDArray((JOrb-1)*nnA+IOrb)
          end do
         end do
        End Do
       END DO
       Call mma_deallocate(GDArray)
       Call mma_deallocate(CIL)
       Call mma_deallocate(CIR)
       RETURN
       END Subroutine

******************************************************

      subroutine CalcW(W,GDMat,PUVX,NPUVX,IndTUVX)
#include "Input.fh"
#include "Pointers.fh"

******Output
      Real*8,DIMENSION((nRoots+1)*nRoots/2,(nRoots+1)*nRoots/2)::W
******Input
      Integer NPUVX
      Real*8,DIMENSION((nRoots+1)*nRoots/2,nnA,nnA)::GDMat
      Real*8,DIMENSION(NPUVX)::PUVX
      INTEGER,DIMENSION(ntAsh,ntAsh,ntAsh,ntAsh)::IndTUVX
******Auxiliary Quantities
      INTEGER K,L,M,N,IKL,IMN,it,iu,iv,ix

      DO K=1,nRoots
       DO L=1,K
       IKL=(K-1)*K/2+L
       Do M=1,nRoots
        Do N=1,M
         IMN=(M-1)*M/2+N
         W(IKL,IMN)=0.0d0
         do it=1,nnA
          do iu=1,nnA
           do iv=1,nnA
            do ix=1,nnA
             IF(IndTUVX(it,iu,iv,ix).ne.0) THEN
            W(IKL,IMN)=W(IKL,IMN)+GDMat(IKL,it,iu)*GDMat(IMN,iv,ix)*
     &       PUVX(IndTUVX(it,iu,iv,ix))
             END IF
            end do
           end do
          end do
         end do
        End Do
       End Do
       END DO
      END DO

      RETURN
      End Subroutine

      subroutine CalcAXX(AXX,W)
#include "Input.fh"
#include "Pointers.fh"
******Input
      Real*8,DIMENSION((nRoots+1)*nRoots/2,(nRoots+1)*nRoots/2)::W
******Output
      Real*8,DIMENSION(((nRoots-1)*nRoots/2)**2)::AXX

******Auxiliary Quantities
      INTEGER K,L,M,N,IKL,IMN,IKL2,IMN2,IKK,ILL,IMM,INN,IC,nRTri
      Real*8  VKLMN,VLKNM,VKLNM,VLKMN

      nRTri=(nRoots-1)*nRoots/2
      DO K=1,nRoots
      DO L=1,K-1
       IKL=(K-1)*K/2+L
       IKK=(K+1)*K/2
       ILL=(L+1)*L/2
       IKL2=(K-2)*(K-1)/2+L
        Do M=1,nRoots
        Do N=1,M-1
         IMN=(M-1)*M/2+N
         IMM=(M+1)*M/2
         INN=(N+1)*N/2
         IMN2=(M-2)*(M-1)/2+N
         VKLMN=0.0d0
         VLKNM=0.0d0
         VLKMN=0.0d0
         VKLNM=0.0d0
         IF(L.eq.M) THEN
          If(N.lt.K) Then
           IC=(K-1)*K/2+N
          Else
           IC=(N-1)*N/2+K
          End If
          VKLMN=W(IC,IKK)+W(IC,INN)-2.0d0*W(IC,ILL)-4.0d0*W(IKL,IMN)
         END IF
         IF(K.eq.N) THEN
          If(M.lt.L) Then
           IC=(L-1)*L/2+M
          Else
           IC=(M-1)*M/2+L
          End If
          VLKNM=W(IC,ILL)+W(IC,IMM)-2.0d0*W(IC,IKK)-4.0d0*W(IKL,IMN)
         END IF
         IF(K.eq.M) THEN
          If(N.lt.L) Then
           IC=(L-1)*L/2+N
          Else
           IC=(N-1)*N/2+L
          End If
          VLKMN=W(IC,ILL)+W(IC,INN)-2.0d0*W(IC,IKK)-4.0d0*W(IKL,IMN)
         END IF
         IF(L.eq.N) THEN
          If(M.lt.K) Then
           IC=(K-1)*K/2+M
          Else
           IC=(M-1)*M/2+K
          End If
          VKLNM=W(IC,IKK)+W(IC,IMM)-2.0d0*W(IC,ILL)-4.0d0*W(IKL,IMN)
         END IF
         AXX((IKL2-1)*nRTri+IMN2)=VKLMN+VLKNM-VKLNM-VLKMN
        End Do
        End Do
       END DO
       END DO
      RETURN
      END SUBROUTINE
******************************************************
******************************************************
      Subroutine Get_Ntri(nTri)
#include "Input.fh"

      INTEGER nTri,kSym
      nTri=0
      DO kSym=1,nSym
       nTri=nTri+nBas(kSym)*(nBas(kSym)+1)/2
      END DO
      RETURN
      END Subroutine
******************************************************

******************************************************
      Subroutine GetPDFTFocks(FMO1t,FMO2t,nTri)
#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "sa.fh"

      INTEGER nTri
      Real*8,DIMENSION(nRoots*nTri)::FMO1t
      Real*8,DIMENSION(nRoots*NACPR2)::FMO2t
      CALL Get_DArray('F1MS            ',FMO1t,nRoots*nTri  )
      CALL Get_DArray('F2MS            ',FMO2t,nRoots*NACPR2)
      RETURN
      end subroutine
******************************************************
