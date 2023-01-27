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
      SUBROUTINE DERFG3(CI,CLAG,DG1,DG2,DG3,DF1,DF2,DF3,idxG3,
     *                  DEPSA,G1,G2)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
      use caspt2_output, only:iPrGlb,verbose,debug
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "SysDef.fh"
#include "WrkSpc.fh"
#include "pt2_guga.fh"

      LOGICAL RSV_TSK

      REAL*8, INTENT(IN) :: CI(MXCI)
      INTEGER*1 idxG3(6,*)
      REAL*8 DEPSA(NLEV,NLEV)
      INTEGER, PARAMETER :: I1=KIND(idxG3)

      REAL*8 DG1(NLEV,NLEV),DG2(NLEV,NLEV,NLEV,NLEV),
     *       DG3(*),DF1(NLEV,NLEV),DF2(NLEV,NLEV,NLEV,NLEV),DF3(*),
     *       CLAG(NCONF)
      REAL*8 G1(NLEV,NLEV),G2(NLEV,NLEV,NLEV,NLEV)

      INTEGER I,J,IDX,JDX
      INTEGER IB,IBMN,IBMX,IBUF,NB,NBTOT,IBUF1
      INTEGER IP1,IP2,IP3,IP1MN,IP1MX,IP1I,IP1STA,IP1END,IP3MX,IQ1
      INTEGER IG3,IG3OFF,IG3BK
      INTEGER ISTU,ISVX,ISYZ
      INTEGER IT,IU,IV,IX,IY,IZ
      INTEGER ITLEV,IULEV,IVLEV,IXLEV,IYLEV,IZLEV,IXLEV0
      INTEGER LBUF1,LBUF2,LBUFD,LBUFT,LBUF3,LBUF4,LBUFX,
     *        LDTU,LDYZ,LDAB
      INTEGER NBUF1,NBUF2,NBUFD,NBUFT,NBUF3,NBUF4,NBUFX,
     *        NDTU,NDYZ,NDAB
      INTEGER LIBUF1,LIP1STA,LIP1END,LOFFSET,IOFFSET
      INTEGER ISSG1,ISSG2,ISP1
      INTEGER ITASK,ISUBTASK,ID,NTASKS,NSUBTASKS,
     &        LTASK_LIST,MXTASK,MYTASK,MYBUFFER
      INTEGER NSGM1,NSGM2
      INTEGER NTRI1,NTRI2
      INTEGER L,LTO,LFROM
      INTEGER MEMMAX, MEMMAX_SAFE
      INTEGER NLEV2
      INTEGER LDUM
      INTEGER NCI,ICSF

      REAL*8, EXTERNAL :: DDOT_

      ! translation tables for levels i,j to and from pair indices idx
      INTEGER IJ2IDX(MXLEV,MXLEV)
      INTEGER IDX2IJ(2,MXLEV**2)
      INTEGER ICNJ(MXLEV**2)
      INTEGER IP1_BUF(MXLEV**2)

      ! result buffer, maximum size is the largest possible ip1 range,
      ! which is set to nbuf1 later, i.e. a maximum of nlev2 <= mxlev**2
C     REAL*8 BUFR(MXLEV**2)
C
C     INTEGER LFCDer1,LFCDer2
      INTEGER ialev,iblev
      REAL*8  SCAL,ScalG,ScalF
C     REAL*8 tmp,tmp2
C
C Put in zeroes. Recognize special cases:
      IF(nlev.EQ.0) GOTO 999

      IF(NACTEL.EQ.0) GOTO 999

      NCI=NCSF(STSYM)
* This should not happen, but...
      IF(NCI.EQ.0) GOTO 999

C Here, for regular CAS or RAS cases.

C Special pair index idx2ij allows true RAS cases to be handled:
      nlev2=nlev**2
      ntri1=(nlev2-nlev)/2
      ntri2=(nlev2+nlev)/2
      idx=0
      do i=1,nlev-1
        do j=i+1,nlev
          idx=idx+1
          ij2idx(i,j)=idx
          idx2ij(1,idx)=i
          idx2ij(2,idx)=j
          jdx=nlev2+1-idx
          ij2idx(j,i)=jdx
          idx2ij(1,jdx)=j
          idx2ij(2,jdx)=i
        end do
      end do
      do i=1,nlev
        idx=ntri1+i
        ij2idx(i,i)=idx
        idx2ij(1,idx)=i
        idx2ij(2,idx)=i
      end do
      do idx=1,nlev2
        i=idx2ij(1,idx)
        j=idx2ij(2,idx)
        jdx=ij2idx(j,i)
        icnj(idx)=jdx
      end do
C
* Correction to G3: It is now <0| E_tu E_vx E_yz |0>
* Similar for F3 values.
* These corrections are already done in CLagDXA_FG3 and CLagDXC_FG3
C
      !! Some preparations
      !! For DF1
      !! F_{tu}
      !! = D_{tu,vw}*f_{vw}
      !! = <0|E_{tu,vw}|0>*f_{vw}
      !! = <0|E_{tu}E_{vw} - \delta_{uv}E_{tw}|0>*f_{vw}
      !! = <0|E_{tu}E_{vw}|0>*f_{vw} - delta_{uv}<0|E_{tw}|0>*f_{vw}
      !! = <0|E_{tu}E_{vv}|0>*f_{vv} - <0|E_{tu}|0>*f_{uu}

      !!  <0|E_{tu}E_{yz}|0>
      !!= <0|t+ u y+ z|0>
      !!= <0|t+ (delta(uy) - y+ u) z|0>
      !!= delta(uy) <0|t+ z|0> - <0|t+ y+ u z|0>
      !! G2(y,t,u,z) = EX2(t,u,y,z) - delta(uy) G1(t,z)
      !! -> G2(u,t,u,z) = EX2(t,u,u,z) - G1(t,z)

      !!   r+ i s+ j - del(is) r+ j
      !! = r+ (i s+ - del(is)) j
      !! = r+ (-s+ i) j
      !! = -r+ s+ i j
      !! = D_{rs,ij}
C     write(6,*) "EPSA"
C     do i = 1, 5
C       write(6,'(i3,f20.10)') i,epsa(i)
C     end do

C     !! F2_{tuvw}
C     !! = <0|E_{tu,vw,xy}|0>*f_{xy}
C     !! = <0|t+ v+ x+ y w u|0> * f_{xy}
C     !!   t+ v+ x+ y w u
C     !! = t+ x+ v+ w y u
C     !! = t+ x+ (del(vw) - w v+) y u
C     !! = del(vw) t+ x+ y u - t+ x+ w v+ y u
C     !! = del(vw)*F1(tu) - t+ (del(xw)-w x+)*(del(vy)-y v+) u
C     !! = del(vw)*F1(tu) - del(xw)*del(vy)*t+ u
C     !!   + del(xw)* t+ y v+ u + del(vy)*t+ w x+ u - t+ w x+ y v+ u
C     !! = del(vw)*F1(tu) - del(xw)*del(vy)*G1_{tu}*f_{xy}
C     !!   + del(xw)*<0|E_{ty}E_{vu}|0>*f_{xy} + del(vy)*<0|E_{tw}E_{xu}|0>*f_{xy}
C     !!   - <0|E_{tw}E_{xy}E_{vu}|0>*f_{xy}
C     !! (if canonical)
C     !! = F1(tu) - G1{tu}*f_{wv} + <0|E_{tw}E_{vu}|0*f_{ww}
C     !!   + <0|E_{tw}E_{vu}|0>*f_{vv} - <0|E_{tw}tE_{xy}E_{vu}|0>
C     !!
C     !!   E_{tu}E_{vw}E_{yz}*f_{vw}
C     !! = t+ u v+ w y+ z * f_{vw}
C     !! = t+ (del(uv) - v+ u)*(del(wy) - y+ w) z * f_{vw}
C     !! = del(uv)*del(wy)*t+ z f_{vw}
C     !!   - del(uv)* t+ y+ w z * f_{vw}
C     !!   - del(wy)* t+ v+ u z * f_{vw}
C     !!   + t+ v+ u y+ w z * f_{vw}
C     !! = del(uv)*del(wy)*G1(tz)*f_{vw}
C     !!   - del(uv)*E_{ty,zw}*f_{vw}
C     !!   - del(wy)*E_{tv,zu}*f_{vw}
C     !!   + t+ v+ (del(uy)-y+ u) w z * f_{vw}
C     !! = del(uv)*del(wy)*G1(tz)*f_{vw}
C     !!   - del(uv)*E_{ty,zw}*f_{vw}
C     !!   - del(wy)*E_{tv,zu}*f_{vw}
C     !!   + del(uy) t+ v+ w z * f_{vw}
C     !!   - t+ v+ y+ u w z * f_{vw}
C     !! = del(uv)*del(wy)*G1(tz)*f_{vw}
C     !!   - del(uv)*E_{ty,zw}*f_{vw}
C     !!   - del(wy)*E_{tv,zu}*f_{vw}
C     !!   + del(uy)*F1_{tz}
C     !!   + t+ y+ v+ w z u * f_{vw}
C     !! = del(uv)*del(wy)*G1(tz)*f_{vw}
C     !!   - del(uv)*E_{ty,zw}*f_{vw}
C     !!   - del(wy)*E_{tv,zu}*f_{vw}
C     !!   + del(uy)*F1_{tz}
C     !!   + F2_{tuyz}

C     !! F2_{tuvw}
C     !! = <0|E_{tu,vw,xy}|0>*f_{xy}
C     !! = <0|t+ v+ x+ y w u|0> * f_{xy}
C     !!   t+ v+ x+ y w u
C     !! =-t+ v+ x+ y u w
C     !! =-t+ x+ v+ u y w * f_{xy}
C     !! =-t+ x+ (del(vu) - u v+) y w * f_{xy}
C     !! =-del(vu) E_{tx,wy}*f_{xy} + t+ x+ u v+ y w * f_{xy}
C     !! =-del(vu) E_{tx,wy}*f_{xy}
C     !!  + t+ (del(ux) - u x+) (del(vy) - y v+) w * f_{xy}
C     !! =-del(vu) E_{tx,wy}*f_{xy}
C     !!  + del(ux)del(vy) G(t,w)*f_{xy} - del(ux) E_{ty}E_{vw}*f_{xy}
C     !!  - del(vy) E_{tu}E_{xw}*f_{xy} + E_{tu}E_{xy}E_{vw}*f_{xy}
C     !! = -del(vu) E_{tx,wy}*f_{xy} + del(ux)del(vy) G1(t,w)*f_{xy}
C     !!   -del(ux) E_{ty}E_{vw}*f_{xy}

C     !! = G1(tz)*f_{uy}
C     !!   + E_{ty,wz}*f_{uw} (in the code, E(it,iw(iu),iy,iz))*e(u))
C     !!   + E_{tv,uz}*f_{vy} (E(it,iu,iv(iy),iz)*e(y))
C     !!   + del(uy)*F1_{tz}
C     !!   + F2_{tuyz}
* Correction to F2: It is now = <0| E_tu H0Diag E_yz |0>
C     !! Etu H0Diag Eyz
C     !! = t+ u w+ w y+ z fww
C     !! = del(uw) t+ w y+ z fww - t+ w+ u w y+ z fww
C     !! = t+ u y+ z fuu - del(wy) t+ w+ u z fww + t+ w+ u y+ w z fww
C     !! = t+ u y+ z fuu - t+ y+ u z fyy + del(uy) t+ w+ w z fww  - t+ w+ y+ u w z fww
C     !! = t+ u y+ z fuu - t+ y+ u z fyy + del(uy) Ftz - t+ y+ w+ w u z fww
C     !! = del(uy) Gtz fuu - Gty,zu fuu - Gty,zu fyy + del(uy) Ftz - Fty,zu
C     !! = del(uy) (Ftz+Gtz fuu) + Gty,uz (fuu+fyy) + Fty,uz
C     !!-> del(uy) (F1(tz)+G1(tz) fuu) + G2(tuyz) (fuu+fyy) + F2(tuyz)
* Correction to F2: Some values not computed follow from symmetry
      do ip1=1,nlev2-1
       itlev=idx2ij(1,ip1)
       iulev=idx2ij(2,ip1)
       it=L2ACT(itlev)
       iu=L2ACT(iulev)
       do ip3=ip1+1,nlev2
        iylev=idx2ij(1,ip3)
        izlev=idx2ij(2,ip3)
        iy=L2ACT(iylev)
        iz=L2ACT(izlev)
        SCAL=DF2(iy,iz,it,iu)+DF2(it,iu,iy,iz)
        DF2(it,iu,iy,iz)=0.0D+00
        DF2(iy,iz,it,iu)=Scal
       end do
      end do
C-SVC20100310: took some spurious mirroring of F2 values out
C-of the loops and put them here, after the parallel section has
C-finished, so that GAdSUM works correctly.
      do ip1=ntri2+1,nlev2
       itlev=idx2ij(1,ip1)
       iulev=idx2ij(2,ip1)
       it=L2ACT(itlev)
       iu=L2ACT(iulev)
       do ip3=ntri1+1,ip1
        iylev=idx2ij(1,ip3)
        izlev=idx2ij(2,ip3)
        iy=L2ACT(iylev)
        iz=L2ACT(izlev)
        Scal=DF2(iz,iy,iu,it)+DF2(it,iu,iy,iz)
        DF2(it,iu,iy,iz)=0.0d+00
        DF2(iz,iy,iu,it)=Scal
       end do
      end do
* Correction to F2: It is now = <0| E_tu H0Diag E_yz |0>
       do iz=1,nlev
        do iy=1,nlev
         do iu=1,nlev
          do it=1,nlev
           DG2(it,iu,iy,iz) = DG2(it,iu,iy,iz)
     *      - DF2(it,iu,iy,iz)*(EPSA(iu)+EPSA(iy))
           do ix=1,nlev
            DEPSA(iu,ix) = DEPSA(iu,ix)
     *        - DF2(it,iu,iy,iz)*G2(it,ix,iy,iz)
            DEPSA(ix,iy) = DEPSA(ix,iy)
     *        - DF2(it,iu,iy,iz)*G2(it,iu,ix,iz)
           end do
          end do
         end do
        end do
       end do
       do iz=1,nlev
        do iu=1,nlev
         do it=1,nlev
          DG1(it,iz) = DG1(it,iz) - EPSA(iu)*DF2(it,iu,iu,iz)
          DF1(it,iz) = DF1(it,iz) - DF2(it,iu,iu,iz)
          do iy = 1, nlev
            DEPSA(iu,iy) = DEPSA(iu,iy) - G1(it,iz)*DF2(it,iu,iy,iz)
          end do
         end do
        end do
       end do
C
C
C
* Correction to G2: Some values not computed follow from symmetry
      do ip1=1,nlev2-1
       itlev=idx2ij(1,ip1)
       iulev=idx2ij(2,ip1)
       it=L2ACT(itlev)
       iu=L2ACT(iulev)
       do ip3=ip1+1,nlev2
        iylev=idx2ij(1,ip3)
        izlev=idx2ij(2,ip3)
        iy=L2ACT(iylev)
        iz=L2ACT(izlev)
        SCAL=DG2(iy,iz,it,iu)+DG2(it,iu,iy,iz)
        DG2(it,iu,iy,iz)=0.0D+00
        DG2(iy,iz,it,iu)=Scal
       end do
      end do
C-SVC20100310: took some spurious mirroring of G2 values out
C-of the loops and put them here, after the parallel section has
C-finished, so that GAdSUM works correctly.
      do ip1=ntri2+1,nlev2
       itlev=idx2ij(1,ip1)
       iulev=idx2ij(2,ip1)
       it=L2ACT(itlev)
       iu=L2ACT(iulev)
       do ip3=ntri1+1,ip1
        iylev=idx2ij(1,ip3)
        izlev=idx2ij(2,ip3)
        iy=L2ACT(iylev)
        iz=L2ACT(izlev)
        Scal=DG2(iz,iy,iu,it)+DG2(it,iu,iy,iz)
        DG2(it,iu,iy,iz)=0.0d+00
        DG2(iz,iy,iu,it)=Scal
       end do
      end do
* Correction to G2: It is now = <0| E_tu E_yz |0>
      do iu=1,nlev
       do iz=1,nlev
        do it=1,nlev
         DG1(it,iz) = DG1(it,iz) - DG2(it,iu,iu,iz)
        end do
       end do
      end do
* Additional correction terms for F1
      Do iT = 1, NLEV
        Do iU = 1, NLEV
          !! With the improved algorithm, symmetrizing the DF1
          !! contribution to DF1 here is somehow required
          DG1(iT,iU) = DG1(iT,iU)
     *      - (DF1(iT,iU)+DF1(iU,iT))*EPSA(iU)*0.5d+00
        End Do
      End Do
      Do iT = 1, NLEV
        Do iU = 1, NLEV
          Do iV = 1, NLEV
            Do iX = 1, nLEV
              DEPSA(iV,iX) = DEPSA(iV,iX) + G2(iT,iU,iV,iX)*DF1(iT,iU)
            End Do
          End DO
        End DO
      End Do
C
      Call CLagSym(nLev,DG1,DG2,DF1,DF2,2)
C
* Dummy values necessary for fooling syntax checkers:
      ldum=1
      ! ndum=1
      call getmem('memmx','max','real',ldum,memmax)

* Use *almost* all remaining memory:
      memmax_safe=int(dble(memmax)*0.95D0)

* Buffers to compute CI expansion vectors into:
* <Psi0|E_ip1 | E_ip2 E_ip3|Psi0>
* buf1: bra buffer with E_ip1 excitations of Psi0
*       holds multiple CI vectors (allowed by memory)
* buf2: ket buffer for an E_ip3 excitation of Psi0
* buft: ket buffer for an E_ip2 excitation of E_ip3|Psi0>
* bufd: diagonal matrix elements to compute the F matrix
      nbuf1=max(1,min(nlev2,(memmax_safe-(6+nlev)*mxci)/mxci/3))
      nbuf2= 1
      nbuft= 1
      nbufd= 1
C
      ndtu =max(1,min(nlev2,(memmax_safe-(6+nlev)*mxci)/mxci/3))
      ndyz = 1
      ndab =max(1,min(nlev2,(memmax_safe-(6+nlev)*mxci)/mxci/3))
      nbuf3= 1
      nbuf4= 1
      nbufx= nlev
      CALL GETMEM('BUF1','ALLO','REAL',LBUF1,NBUF1*MXCI)
      CALL GETMEM('BUF2','ALLO','REAL',LBUF2,NBUF2*MXCI)
      CALL GETMEM('BUFT','ALLO','REAL',LBUFT,NBUFT*MXCI)
      CALL GETMEM('BUFD','ALLO','REAL',LBUFD,NBUFD*MXCI)
C
      CALL GETMEM('DTU ','ALLO','REAL',LDTU ,NDTU *MXCI)
      CALL GETMEM('DYZ ','ALLO','REAL',LDYZ ,NDYZ *MXCI)
      CALL GETMEM('DAB ','ALLO','REAL',LDAB ,NDAB *MXCI)
      CALL GETMEM('BUF3','ALLO','REAL',LBUF3,NBUF3*MXCI)
      CALL GETMEM('BUF4','ALLO','REAL',LBUF4,NBUF4*MXCI)
C
      CALL GETMEM('BUFX','ALLO','REAL',LBUFX,NBUFX*MXCI)

C-SVC20100301: calculate maximum number of tasks possible
      MXTASK=(NTRI2-1)/NBUF1+1+(NTRI1-1)/NBUF1+1
      CALL GETMEM ('TASKLIST','ALLO','INTE',lTask_List,4*mxTask)
      lip1sta=lTask_List
      lip1end=lTask_List+mxTask
      libuf1=lTask_List+2*mxTask
      lOffSet=lTask_List+3*mxTask

      IF(iPrGlb.GE.VERBOSE) THEN
        WRITE(6,*)
        WRITE(6,'(2X,A)') 'Constructing derivatives of G3/F3'
        WRITE(6,'(2X,A,F16.9,A)') ' memory avail: ',
     &    (memmax*RtoB)/1.0D9, ' GB'
        WRITE(6,'(2X,A,F16.9,A)') ' memory used:  ',
     &    (((3*nbuf1+6+nlev)*MXCI)*RtoB)/1.0D9, ' GB'
        call xFlush(6)
      ENDIF
C     CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
C     CPUT =CPTF10-CPTF0
C     WALLT=TIOTF10-TIOTF0
C     write(6,*) "PREP    : CPU/WALL TIME=", cput,wallt

      iG3OFF=0
* A *very* long loop over the symmetry of Sgm1 = E_ut Psi as segmentation.
* This also allows precomputing the Hamiltonian (H0) diagonal elements.
      DO issg1=1,nsym
       isp1=mul(issg1,STSYM)
       nsgm1=ncsf(issg1)
       !! Work(LBufD) = \sum_t <I|E_{tt}|I>*f_{tt}
       CALL H0DIAG_CASPT2(ISSG1,WORK(LBUFD),IWORK(LNOW),IWORK(LIOW))

C-SVC20100301: calculate number of larger tasks for this symmetry, this
C-is basically the number of buffers we fill with sigma1 vectors.
      iTask=1
      ibuf1=0
      DO ip1=1,nlev2
        itlev=idx2ij(1,ip1)
        iulev=idx2ij(2,ip1)
        istu=mul(ism(itlev),ism(iulev))
        IF (istu.EQ.isp1) THEN
          ibuf1=ibuf1+1
          ip1_buf(ibuf1)=ip1
          IF (ibuf1.EQ.1) iwork(lip1sta+iTask-1)=ip1
        ENDIF
        IF (ibuf1.EQ.nbuf1.OR.(ibuf1.GT.0.AND.
     &         (ip1.EQ.ntri2.OR.ip1.EQ.nlev2))) THEN
            iwork(lip1end+iTask-1)=ip1_buf(ibuf1)
            iwork(libuf1+iTask-1)=ibuf1
            iTask=iTask+1
            ibuf1=0
        ENDIF
      ENDDO
      nTasks=iTask
C     write(6,*) "nTasks = ", nTasks
      IF (ibuf1.EQ.0) nTasks=nTasks-1
C-SVC20100309: calculate number of inner loop iteration tasks.
      iOffSet=0
      DO iTask=1,nTasks
        iWork(lOffSet+iTask-1)=iOffSet
        ip1sta=iwork(lip1sta+iTask-1)
        ip1end=iwork(lip1end+iTask-1)
        ip3mx=ntri2
        if(ip1end.le.ntri2) ip3mx=ip1end
        if(ip1sta.gt.ntri2) ip3mx=ntri1
C       write(6,*) "iTask = ", iTask
C       write(6,*) "start,end=",ip1sta,ip1end
C       write(6,*) "ip3mx = ", ip3mx
C-SVC20100309: Currently -we are going to limit this to the ip3-loop and
C-leave the ip2-loop intact.  This was based on the large overhead which
C-was observed for a very large number of small tasks.
C       iOffSet=iOffSet+ip3mx*ntri2-((ip3mx**2-ip3mx)/2)
        iOffSet=iOffSet+ip3mx
      ENDDO
      nSubTasks=iOffSet
C     write(6,*) "nSubTasks = ", nSubTasks

      IF(iPrGlb.GE.VERBOSE) THEN
        WRITE(6,'(2X,A,I3,A,I6)') 'Sym: ',issg1,', #Tasks: ',nSubTasks
        call xFlush(6)
      ENDIF

      IF(iPrGlb.GE.DEBUG) THEN
        IF (nSubTasks .GT. 0) THEN
          WRITE(6,'("DEBUG> ",A8,1X,A12,1X,A4,1X,A9)')
C-position 12345678901234567890
     &    "--------",
     &    "------------",
     &    "----",
     &    "---------"
          WRITE(6,'("DEBUG> ",A8,1X,A12,1X,A4,1X,A9)')
C-position 12345678901234567890
     &    "task ID ",
     &    " ip1 range  ",
     &    "ip3 ",
     &    "#elements"
          WRITE(6,'("DEBUG> ",A8,1X,A12,1X,A4,1X,A9)')
C-position 12345678901234567890
     &    "--------",
     &    "------------",
     &    "----",
     &    "---------"
          call xFlush(6)
        END IF
      END IF

C-SVC20100301: initialize the series of subtasks
      Call Init_Tsk(ID, nSubTasks)

      myBuffer=0

      !! loop start
 500  CONTINUE
C-SVC20100908: first check: can I actually do any task?
      IF ((NG3-iG3OFF).LT.nbuf1*ntri2) GOTO 501
C-SVC20100831: initialize counter for offset into G3
C-SVC20100302: BEGIN SEPARATE TASK EXECUTION
C     write(6,*) rsv_tsk(id,isubtask)
      If (.NOT.Rsv_Tsk(ID,iSubTask)) GOTO 501

      myTask=nTasks
      DO iTask=1,nTasks
        iBuf=iSubTask-iWork(lOffSet+iTask-1)
        IF (iBuf.LE.0) THEN
          myTask=iTask-1
          goto 666
        ENDIF
      ENDDO
666   continue
      iTask=myTask

      iOffSet=iWork(lOffSet+iTask-1)

C-SVC20100310: one task handles a range of ip1 values
C-that are in the buffer and one ip3 value, for which
C-a loop over ip2 values is then executed.
      ip1sta=iWork(lip1sta+iTask-1)
      ip1end=iWork(lip1end+iTask-1)
      ip3=iSubTask-iOffSet

C-SVC20100301: fill the buffer with sigma vectors if they
C-have not been computed yet, else just get the number of
C-sigma vectors in the buffer.
C     write(6,*) "myBuffer,iTask = ", myBuffer,iTask
      IF (myBuffer.NE.iTask) THEN
        ibuf1=0
        do ip1i=ip1sta,ip1end
         itlev=idx2ij(1,ip1i)
         iulev=idx2ij(2,ip1i)
         istu=mul(ism(itlev),ism(iulev))
         it=L2ACT(itlev)
         iu=L2ACT(iulev)
         if(istu.eq.isp1) then
          ibuf1=ibuf1+1
          ip1_buf(ibuf1)=ip1i
          lto=lbuf1+mxci*(ibuf1-1)
          call dcopy_(nsgm1,[0.0D0],0,work(lto),1)
          CALL SIGMA1_CP2(IULEV,ITLEV,1.0D00,STSYM,CI,WORK(LTO),
     &     IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &     IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &     WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
         end if
        end do
        myBuffer=iTask
        Call DCopy_(MXCI*ibuf1,[0.0D+00],0,Work(LDTU),1)
        Call DCopy_(MXCI*ibuf1,[0.0D+00],0,Work(LDAB),1)
      ELSE
        ibuf1=iWork(libuf1+iTask-1)
      ENDIF
C-SVC20100301: necessary batch of sigma vectors is now in the buffer

      ! The ip1 buffer could be the same on different processes
      ! so only compute the G1 contribution when ip3 is 1, as
      ! this will only be one task per buffer.
      if (issg1.eq.STSYM.AND.ip3.eq.1) then
        !! lbuf1 = <Psi0|E_ip1|I>
        !! <0|E_{tu}I> = <I|E_{ut}|0>
C       write(6,*) "ib loop"
        do ib=1,ibuf1
          idx=ip1_buf(ib)
          itlev=idx2ij(1,idx)
          iulev=idx2ij(2,idx)
          it=L2ACT(itlev)
          iu=L2ACT(iulev)
C         write(6,'("itlev,iulev,it,iu = ",4i3)') itlev,iulev,it,iu
          lto=lbuf1+mxci*(ib-1)
C         write(6,'(5f20.10)') (work(lto+i-1),i=1,nsgm1)
C         write(6,'(f20.10)') ddot_(nsgm1,work(lto),1,ci,1)
          !! DG1 contribution
          SCAL = DG1(iT,iU) + DG1(iT,iU)
          Call DaXpY_(nsgm1,SCAL,Work(lto),1,CLag,1)
C
          !! left derivative of DF1
          ibuf = ldtu + mxci*(ib-1)
          Do icsf = 1, nsgm1
            Work(ibuf+icsf-1) = Work(ibuf+icsf-1)
     *        + DF1(it,iu)*Work(LBUFD+icsf-1)*CI(icsf)
          End Do
          !! right derivative of DF1
          Do icsf = 1, nsgm1
            CLag(icsf) = CLag(icsf)
     *        + DF1(it,iu)*Work(lto+icsf-1)*Work(LBUFD+icsf-1)
          End Do

C         G1(it,iu)=DDOT_(nsgm1,ci,1,work(lto),1)
C         IF(IFF.ne.0) then
C           F1sum=0.0D0
C           do i=1,nsgm1
C             F1sum=F1sum+CI(i)*work(lto-1+i)*work(lbufd-1+i)
C           end do
C           F1(it,iu)=F1sum-EPSA(iu)*G1(it,iu)
C         end if
        end do
      end if

C     ip3mx=ntri2
C     if(ip1end.le.ntri2) ip3mx=ip1end
C     if(ip1sta.gt.ntri2) ip3mx=ntri1
C-SVC20100309: loop over ip3, ip2
C     do ip3=1,ip3mx

C-SVC20100309: PAM's magic formula
*     iCnt=iSubTask-iOffSet
*     ip3=int(dble(ntri2)+1.5D0 -
*    &     sqrt((dble(ntri2)+0.5d0)**2-2*iCnt+0.000001D0))
*     ip2=iCnt-((ip3-1)*ntri2-((ip3-1)*(ip3-2))/2 )+ip3-1

C-SVC20100309: use simpler procedure by keeping inner ip2-loop intact

C     write(6,*) "ip3 = ", ip3
C     CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
      iq1=icnj(ip3)
* The indices corresponding to pair index p3:
      iylev=idx2ij(1,ip3)
      izlev=idx2ij(2,ip3)
      isyz=mul(ism(iylev),ism(izlev))
      issg2=mul(isyz,STSYM)
      nsgm2=ncsf(issg2)
      iy=L2ACT(iylev)
      iz=L2ACT(izlev)
      lto=lbuf2
      call dcopy_(nsgm2,[0.0D0],0,work(lto),1)
      CALL SIGMA1_CP2(IYLEV,IZLEV,1.0D00,STSYM,CI,WORK(LTO),
     &     IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &     IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &     WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
      Call Dcopy_(nsgm1,[0.0D+00],0,Work(LDYZ),1)
      if(issg2.eq.issg1) then
        call dcopy_(nsgm2,[0.0D0],0,work(lbuf3),1)
        do ib=1,ibuf1
          idx=ip1_buf(ib)
          itlev=idx2ij(1,idx)
          iulev=idx2ij(2,idx)
          it=L2ACT(itlev)
          iu=L2ACT(iulev)
C
          ScalG=DG2(iT,iU,iY,iZ)
          ScalF=DF2(iT,iU,iY,iZ)
          If (ScalG.eq.0.0d+00.and.ScalF.eq.0.0D+00) Cycle
C
          !! left derivative
          Do icsf = 1, nsgm1
            Work(LBUFT+icsf-1)
     *        = ScalG*Work(LTO+icsf-1)
     *        + ScalF*Work(LTO+icsf-1)*Work(LBUFD+icsf-1)
          End Do
          ibuf = ldtu + mxci*(ib-1)
          Call DaXpY_(nsgm1,1.0d+00,Work(LBUFT),1,Work(ibuf),1)
C
          !! right derivative
          call dcopy_(mxci,work(lbuf1+mxci*(ib-1)),1,work(lbuft),1)
          Do icsf = 1, nsgm1
            Work(LBUF3+icsf-1) = Work(LBUF3+icsf-1)
     *        + ScalG*Work(LBUFT+icsf-1)
     *        + ScalF*Work(LBUFT+icsf-1)*Work(LBUFD+icsf-1)
          End Do
C
          !! For DEPSA
          ibuf = ldab + mxci*(ib-1)
          Call DaXpY_(nsgm1,ScalF,Work(LTO),1,Work(IBUF),1)
        end do
        !! Save for Eyz
        Call DaXpY_(nsgm2,1.0d+00,Work(LBUF3),1,Work(LDYZ),1)
      end if
      nbtot=0
C
C
C
      !! Prepare for DEPSA for the -epsa(iv) term with square
      lfrom=lbuf2
      iG3bk = iG3OFF
      Do ixlev0 = 1, nlev
        lfrom=lbuf2
        Do ivlev = 1, nlev
          L = LBUFX + MXCI*(ivlev-1)
          Call DCopy_(nsgm1,[0.0D0],0,Work(L),1)
          CALL SIGMA1_CP2(IVLEV,IXLEV0,1.0D+0,STSYM,Work(LFROM),Work(L),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        End Do
        iG3OFF = iG3bk
      do ip2=ip3,ntri2
        ivlev=idx2ij(1,ip2)
        ixlev=idx2ij(2,ip2)
        isvx=mul(ism(ivlev),ism(ixlev))
        iv=L2ACT(ivlev)
        ix=L2ACT(ixlev)
        if(isvx.ne.mul(issg1,issg2)) goto 99
C       lfrom=lbuf2
        !! <I|EvxEyz|0>
        lto=lbuft
        If (IXLEV.EQ.IXLEV0) THEN
        L = LBUFX + MXCI*(IVLEV-1)
        Call DCopy_(nsgm1,Work(L),1,Work(LBUFT),1)
        END IF
*-----------
* Max and min values of index p1:
        ip1mx=ntri2
        if(ip3.le.ntri1) then
          ip1mx=nlev2
          if(ip2.gt.ntri1) ip1mx=iq1
        end if
        ip1mn=max(ip2,ip1sta)
        ip1mx=min(ip1mx,ip1end)
* The corresponding locations in the Sgm1 buffer:
        ibmn=999999
        ibmx=-999999
        do ib=ibuf1,1,-1
          ip1=ip1_buf(ib)
          if(ip1.ge.ip1mn)ibmn=ib
        end do
        do ib=1,ibuf1
          ip1=ip1_buf(ib)
          if(ip1.le.ip1mx)ibmx=ib
        end do
        nb=ibmx-ibmn+1
        if(nb.le.0) goto 99
        if (ixlev.ne.ixlev0) then
          iG3OFF = iG3OFF + nb
          cycle
        end if

C
C       ----- left derivative
C
        do icsf = 1, nsgm1
          ! BUF3 = (<I|Ett|I>-EPSA(V))*<I|EvxEyz|0> = <I|fEvxEyz|0>
          work(lbuf3-1+icsf)
     *      = (work(lbufd-1+icsf)-epsa(iv))*work(lbuft-1+icsf)
        end do
        do ib=1,nb
          iG3=iG3OFF+ib
          idx=ip1_buf(ibmn-1+ib)
C
C         ibuf = ldtu + mxci*(idx-ip1sta)
          ibuf = ldtu + mxci*(ibmn+ib-2)
          !! <I|EvxEyz|0>*Dtuvxyz
          Call DaXpY_(nsgm1,DG3(iG3),Work(LBUFT),1,Work(ibuf),1)
          !! <I|fEvxEyz|0>*Ftuvxyz
          Call DaXpY_(nsgm1,DF3(iG3),Work(LBUF3),1,Work(ibuf),1)
          !! DEPSA of the Work(LBUFD) term
C         ibuf = ldab + mxci*(idx-ip1sta)
          ibuf = ldab + mxci*(ibmn+ib-2)
          Call DaXpY_(nsgm1,DF3(iG3),Work(LBUFT),1,Work(ibuf),1)
        end do
C
C       ----- right derivative
C
        Call DCopy_(nsgm1,[0.0D0],0,Work(LBUF3),1)
        Call DCopy_(nsgm1,[0.0D0],0,Work(LBUF4),1)
        !! right derivative (1): <0|Etu|I>*Dtuvxyz and <0|Etu|I>*Ftuvxyz
        !! <0|EtuEvxEyz|I> -> <I|EzyExvEut|0>
        do ib=1,nb
          iG3=iG3OFF+ib
          idx=ip1_buf(ibmn-1+ib)
C
C         ibuf = lbuf1+mxci*(idx-ip1sta)
          ibuf = lbuf1+mxci*(ibmn+ib-2)
          !! BUF3 = <0|Etu|I>*Dtuvxyz
          Call DaXpY_(nsgm1,DG3(iG3),Work(IBUF),1,Work(LBUF3),1)
          !! BUFC = <0|Etu|I>*Ftuvxyz
          Call DaXpY_(nsgm1,DF3(iG3),Work(IBUF),1,Work(LBUF4),1)
        end do
C
        !! DEPSA of the -EPSA(iv) term
        Call DGEMV_('T',nsgm1,NLEV,
     *             -1.0d+00,Work(lbufx),mxci,
C    *             -1.0d+00,Work(lbufx+mxci*nlev*(ixlev-1)),mxci,
     *                      Work(lbuf4),1,
     *              1.0d+00,DEPSA(1,IVLEV),1)
C
        !! Scale the DF3 contribution with the diagonal Fock
        !! and add to the DG3 contribution
        do icsf = 1, nsgm1
          work(lbuf3-1+icsf) = work(lbuf3-1+icsf)
     *      + work(lbuf4+icsf-1)*(work(lbufd+icsf-1)-epsa(iv))
        end do
        !! right derivative (2): <0|EtuEvx|I>*Dtuvxyz
       CALL SIGMA1_CP2(IXLEV,IVLEV,1.0D+00,STSYM,WORK(LBUF3),WORK(LDYZ),
     &      IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &      IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &      WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C
        iG3OFF=iG3OFF+nb
        nbtot=nbtot+nb
 99     continue
      end do !! end of ip2 loop
      End Do !! end of ixlev0 loop
C
      !! Complete the right derivative contribution:
      !! <0|EtuEyz|I> and <0|EtuEvxEyz|I>
      CALL SIGMA1_CP2(IZLEV,IYLEV,1.0D+00,STSYM,WORK(LDYZ),CLAG,
     &     IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &     IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &     WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C
      IF(iPrGlb.GE.DEBUG) THEN
        WRITE(6,'("DEBUG> ",I8,1X,"[",I4,"..",I4,"]",1X,I4,1X,I9)')
     &    iSubTask, ip1sta, ip1end, ip3, nbtot
        call xFlush(6)
      END IF

CSVC: The master node now continues to only handle task scheduling,
C     needed to achieve better load balancing. So it exits from the task
C     list.  It has to do it here since each process gets at least one
C     task.
#if defined (_MOLCAS_MPP_) && !defined (_GA_)
      IF (IS_REAL_PAR().AND.KING().AND.(NPROCS.GT.1)) GOTO 501
#endif
C
      !! Complete the left derivative and DEPSA contribution
      If ((ip1end.le.ntri2.and.ip3.eq.ip1end).or.
     *    (ip1sta.gt.ntri2.and.ip3+ntri2.eq.nlev2)) Then
        do ib=1,ibuf1
          idx=ip1_buf(ib)
          itlev=idx2ij(1,idx)
          iulev=idx2ij(2,idx)
          lto=ldtu+mxci*(ib-1)
          !! left derivative
          CALL SIGMA1_CP2(ITLEV,IULEV,1.0D00,STSYM,WORK(LTO),CLAG,
     &     IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &     IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &     WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          !! the rest is DEPSA contribution
          IBUF = LDAB + MXCI*(ib-1)
          Do IALEV = 1, NLEV
            Do IBLEV = 1, NLEV
              Call DCopy_(nsgm1,[0.0D0],0,Work(LBUF2),1)
       CALL SIGMA1_CP2(IALEV,IBLEV,1.0D+00,STSYM,Work(IBUF),Work(LBUF2),
     &          IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &          IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &          WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
              DEPSA(IALEV,IBLEV) = DEPSA(IALEV,IBLEV)
     *          + DDot_(nsgm1,Work(LBUF1+MXCI*(IB-1)),1,Work(LBUF2),1)
            End Do
          End Do
        end do
      End If

C-SVC20100301: end of the task
      GOTO 500

 501  CONTINUE


C-SVC20100302: no more tasks, wait here for the others, then proceed
C with next symmetry
      CALL Free_Tsk(ID)

      IF(iPrGlb.GE.DEBUG) THEN
        IF (nSubTasks .GT. 0) THEN
          WRITE(6,'("DEBUG> ",A8,1X,A12,1X,A4,1X,A9)')
C-position 12345678901234567890
     &    "--------",
     &    "------------",
     &    "----",
     &    "---------"
        END IF
      END IF

* End of sectioning loop over symmetry of Sgm1 wave functions.
      END DO
C
      CALL GETMEM ('TASKLIST','FREE','INTE',lTask_List,4*mxTask)
      ! free CI buffers
      CALL GETMEM('BUF1','FREE','REAL',LBUF1,NBUF1*MXCI)
      CALL GETMEM('BUF2','FREE','REAL',LBUF2,NBUF2*MXCI)
      CALL GETMEM('BUFT','FREE','REAL',LBUFT,NBUFT*MXCI)
      CALL GETMEM('BUFD','FREE','REAL',LBUFD,NBUFD*MXCI)
C
      CALL GETMEM('DTU ','FREE','REAL',LDTU ,NDTU *MXCI)
      CALL GETMEM('DYZ ','FREE','REAL',LDYZ ,NDYZ *MXCI)
      CALL GETMEM('DAB ','FREE','REAL',LDAB ,NDAB *MXCI)
C
      CALL GETMEM('BUF3','FREE','REAL',LBUF3,NBUF3*MXCI)
      CALL GETMEM('BUF4','FREE','REAL',LBUF4,NBUF4*MXCI)
      CALL GETMEM('BUFX','FREE','REAL',LBUFX,NBUFX*MXCI)
C
 999  continue
      RETURN
      END
C
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
* Copyright (C) 2006, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 2006  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE DERSPE(DF1,DF2,DF3,idxG3,DEPSA,G1,G2,G3)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DF1(NASHT,NASHT),DF2(NASHT,NASHT,NASHT,NASHT),DF3(*)
      DIMENSION G1(NASHT,NASHT),G2(NASHT,NASHT,NASHT,NASHT),G3(*)
      DIMENSION DEPSA(NASHT,NASHT)
      INTEGER*1 idxG3(6,*)
      INTEGER I1
      PARAMETER (I1=KIND(idxG3))
C SPECIAL-CASE ROUTINE. DELIVERS G AND F MATRICES FOR A HIGH-SPIN
C OR CLOSED-SHELL SCF CASE.
#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"

      LOGICAL RSV_TSK

      ESUM=0.0D0
      DESUM=0.0D+00
      DO I=1,NLEV
        ESUM=ESUM+ETA(I)
      END DO
C ISCF=1 for closed-shell, =2 for hispin
      OCC=2.0D0
      IF(ISCF.EQ.2) OCC=1.0D0

      IF(NACTEL.EQ.1) THEN
        NG3=0
        GOTO 111
      END IF
C
      IF(NACTEL.EQ.2) THEN
        NG3=0
        GOTO 222
      END IF
C
C
C
      write(6,*) "I have not implemented for non-standard Psi0",
     *", when A and C subspaces contribute to the energy, in particular"
      write(6,*)"I cannot debug, because I do not know when it happens"
C     call abend
C
      NLEV2=NLEV**2
      NLEV4=NLEV**4

      iG3=0
      nTask=NLEV4
C SVC20100908 initialize the series of tasks
      Call Init_Tsk(ID, nTask)

 500  CONTINUE
#ifdef _MOLCAS_MPP_
      IF ((NG3-iG3).LT.NLEV2) GOTO 501
#endif
 502  IF (.NOT.Rsv_Tsk(ID,iTask)) GOTO 501

      IND1=MOD(iTask-1,NLEV2)+1
      IND2=((iTask-IND1)/(NLEV2))+1
      IF(IND2.GT.IND1) GOTO 502

      IT1=MOD(IND1-1,NASHT)+1
      IU1=(IND1-IT1)/NASHT+1
      LU1=LEVEL(IU1)
      IT2=MOD(IND2-1,NASHT)+1
      IU2=(IND2-IT2)/NASHT+1
      LU2=LEVEL(IU2)

      DO IT3=1,NLEV
       DO IU3=1,NLEV
        IND3=IT3+NASHT*(IU3-1)
        IF(IND3.GT.IND2) GOTO 198
        LU3=LEVEL(IU3)
C       VAL=G1(IT1,IU1)*G1(IT2,IU2)*G1(IT3,IU3)

C Here VAL is the value <PSI1|E(IT1,IU1)E(IT2,IU2)E(IT3,IU3)|PSI2>
C Add here the necessary Kronecker deltas times 2-body matrix
C elements and lower, so we get a true normal-ordered density matrix
C element.

C <PSI1|E(T1,U1,T2,U2,T3,U3)|PSI2>
C = <PSI1|E(T1,U1)E(T2,U2)E(T3,U3)|PSI2>
C -D(T3,U2)*(G2(T1,U1,T2,U3)+D(T2,U1)*G1(T1,U3))
C -D(T2,U1)*G2(T1,U2,T3,U3)
C -D(T3,U1)*G2(T2,U2,T1,U3)

C       IF(IT3.EQ.IU2) THEN
C         VAL=VAL-G2(IT1,IU1,IT2,IU3)
C         IF(IT2.EQ.IU1) THEN
C           VAL=VAL-G1(IT1,IU3)
C         END IF
C       END IF
C       IF(IT2.EQ.IU1) THEN
C         VAL=VAL-G2(IT1,IU2,IT3,IU3)
C       END IF
C       IF(IT3.EQ.IU1) THEN
C         VAL=VAL-G2(IT2,IU2,IT1,IU3)
C       END IF

C VAL is now =<PSI1|E(IT1,IU1,IT2,IU2,IT3,IU3)|PSI2>
        iG3=iG3+1
        idxG3(1,iG3)=INT(iT1,I1)
        idxG3(2,iG3)=INT(iU1,I1)
        idxG3(3,iG3)=INT(iT2,I1)
        idxG3(4,iG3)=INT(iU2,I1)
        idxG3(5,iG3)=INT(iT3,I1)
        idxG3(6,iG3)=INT(iU3,I1)
C       G3(iG3)=VAL
C       F3(iG3)=(ESUM*OCC-ETA(LU1)-ETA(LU2)-ETA(LU3))*VAL
        DESUM = DESUM + OCC*G3(iG3)*DF3(iG3)
        DEPSA(LU1,LU1) = DEPSA(LU1,LU1) - OCC*G3(iG3)*DF3(iG3)
        DEPSA(LU2,LU2) = DEPSA(LU2,LU2) - OCC*G3(iG3)*DF3(iG3)
        DEPSA(LU3,LU3) = DEPSA(LU3,LU3) - OCC*G3(iG3)*DF3(iG3)

 198    CONTINUE
        END DO
      END DO

CSVC: The master node now continues to only handle task scheduling,
C     needed to achieve better load balancing. So it exits from the task
C     list.  It has to do it here since each process gets at least one
C     task.
#if defined (_MOLCAS_MPP_) && !defined (_GA_)
      IF (IS_REAL_PAR().AND.KING().AND.(NPROCS.GT.1)) GOTO 501
#endif

      GO TO 500
 501  CONTINUE

C SVC2010: no more tasks, wait here for the others.
      CALL Free_Tsk(ID)

      NG3=iG3
C
C
C
 222  CONTINUE
      DO IT=1,NASHT
       LT=LEVEL(IT)
       DO IU=1,NASHT
        LU=LEVEL(IU)
C       G2(IT,IT,IU,IU)=G1(IT,IT)*G1(IU,IU)
C       IF(IU.EQ.IT) THEN
C        G2(IT,IT,IU,IU)=G2(IT,IT,IU,IU)-G1(IT,IU)
C       ELSE
C        G2(IT,IU,IU,IT)=-G1(IT,IT)
C       END IF
C       F2(IT,IT,IU,IU)=(ESUM*OCC-ETA(LT)-ETA(LU))*G2(IT,IT,IU,IU)
C       F2(IT,IU,IU,IT)=(ESUM*OCC-ETA(LT)-ETA(LU))*G2(IT,IU,IU,IT)
        DESUM = DESUM + OCC*G2(IT,IT,IU,IU)*DF2(IT,IT,IU,IU)
        DESUM = DESUM + OCC*G2(IT,IU,IU,IT)*DF2(IT,IU,IU,IT)
        DO IV=1,NASHT
          LV=LEVEL(IV)
          DEPSA(LT,LV)=DEPSA(LT,LV)-OCC*G2(IT,IT,IU,IU)*DF2(IT,IV,IU,IU)
          DEPSA(LU,LV)=DEPSA(LU,LV)-OCC*G2(IT,IT,IU,IU)*DF2(IT,IT,IU,IV)
          DEPSA(LT,LV)=DEPSA(LT,LV)-OCC*G2(IT,IU,IU,IT)*DF2(IT,IU,IU,IV)
          DEPSA(LU,LV)=DEPSA(LU,LV)-OCC*G2(IT,IU,IU,IT)*DF2(IT,IU,IV,IT)
        END DO
       END DO
      END DO
C
 111  CONTINUE
      DO IT=1,NASHT
C       G1(IT,IT)=OCC
        LT=LEVEL(IT)
C       F1(IT,IT)=(ESUM*OCC-ETA(LT))*G1(IT,IT)
        DESUM = DESUM + OCC*G1(IT,IT)*DF1(IT,IT)
        Do IU=1, NASHT
          LU=LEVEL(IU)
          DEPSA(LT,LU) = DEPSA(LT,LU) - OCC*G1(IT,IT)*DF1(IT,IU)
        End Do
      END DO

      RETURN
      END
