      SUBROUTINE DERFG3(IFF,CI,CLAG,DG1,DG2,DG3,DF1,DF2,DF3,idxG3,
     *                  DEPSA,G1,G2,G3,F1,F2,F3,TRF)
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "SysDef.fh"
#include "WrkSpc.fh"
#include "pt2_guga.fh"

#include "para_info.fh"
      LOGICAL RSV_TSK

      INTEGER, INTENT(IN) :: IFF
      REAL*8, INTENT(IN) :: CI(MXCI)
C     REAL*8, INTENT(OUT) :: DG1(NLEV,NLEV),DG2(NLEV,NLEV,NLEV,NLEV)
C     REAL*8, INTENT(OUT) :: DF1(NLEV,NLEV),DF2(NLEV,NLEV,NLEV,NLEV)
C     REAL*8, INTENT(OUT) :: DG3(*), DF3(*)
      INTEGER*1 idxG3(6,*)
      REAL*8 DEPSA(NLEV,NLEV)
      INTEGER, PARAMETER :: I1=KIND(idxG3)

      REAL*8 DG1(NLEV,NLEV),DG2(NLEV,NLEV,NLEV,NLEV),
     *       DG3(*),DF1(NLEV,NLEV),DF2(NLEV,NLEV,NLEV,NLEV),DF3(*),
     *       CLAG(NCONF)
      REAL*8 G1(NLEV,NLEV),G2(NLEV,NLEV,NLEV,NLEV),G3(*),
     *       F1(NLEV,NLEV),F2(NLEV,NLEV,NLEV,NLEV),F3(*)
      REAL*8 TRF(NLEV,NLEV)
      REAL*8 F1SUM,F2SUM,DF3tmp
C     real*8 df2bk(nlev,nlev,nlev,nlev)

      INTEGER I,J,IDX,JDX
      INTEGER IB,IBMN,IBMX,IBUF,NB,NBTOT,IBUF1
      INTEGER IP1,IP2,IP3,IP1MN,IP1MX,IP1I,IP1STA,IP1END,IP3MX,IQ1
      INTEGER IG3,IG3OFF
      INTEGER ISTU,ISVX,ISYZ
      INTEGER IT,IU,IV,IX,IY,IZ,iA
      INTEGER ITLEV,IULEV,IVLEV,IXLEV,IYLEV,IZLEV
      INTEGER LBUF1,LBUF2,LBUFD,LBUFT,LBUF3,LBUF4,LBUF5,LBUF6,LBUFE,
     *        LBUFA
      INTEGER NBUF1,NBUF2,NBUFD,NBUFT,NBUF3,NBUF4,NBUF5,NBUF6,NBUFE,
     *        NBUFA
      INTEGER LIBUF1,LIP1STA,LIP1END,LOFFSET,IOFFSET
      INTEGER ISSG1,ISSG2,ISP1
      INTEGER ITASK,ISUBTASK,ID,NTASKS,NSUBTASKS,
     &        LTASK_LIST,MXTASK,MYTASK,MYBUFFER
      INTEGER NSGM1,NSGM2
      INTEGER NTRI1,NTRI2
      INTEGER L1,LTO,LFROM,LTO2
      INTEGER MEMMAX, MEMMAX_SAFE
      INTEGER NLEV2
      INTEGER LDUM,NDUM
      INTEGER NCI,ICSF

      REAL*8, EXTERNAL :: DDOT_,DNRM2_

      ! translation tables for levels i,j to and from pair indices idx
      INTEGER IJ2IDX(MXLEV,MXLEV)
      INTEGER IDX2IJ(2,MXLEV**2)
      INTEGER ICNJ(MXLEV**2)
      INTEGER IP1_BUF(MXLEV**2)

      ! result buffer, maximum size is the largest possible ip1 range,
      ! which is set to nbuf1 later, i.e. a maximum of nlev2 <= mxlev**2
      REAL*8 BUFR(MXLEV**2)
C
      INTEGER LFCDer1,LFCDer2
      INTEGER iEPS,jEPS,ialev,iblev,nadd
      REAL*8  SCAL,ScalGYZL,ScalGZYL,ScalFYZL,ScalFZYL,
     *             ScalGYZR,ScalGZYR,ScalFYZR,ScalFZYR
      REAL*8 tmp,tmp2
      LOGICAL RAS,RASSPE
      integer iw,itu,iyz
      logical dodebug
C
      !!
C     Real*8 CPTF0,CPTF10,CPE,TIOTF0,TIOTF10,TIOE,CPTF,CPUT,WALLT
C
      dodebug = .false.
      if (.not.dodebug) ng3 = ng3*6
      RAS = NRAS1T+NRAS3T.GT.0
C     do it=1,nlev
C     do iu=1,nlev
C     do iy=1,nlev
C     do iz=1,nlev
C     write(6,'(4i3,f20.10)') it,iu,iy,iz,g2(it,iy,iu,iz)
C     end do
C     end do
C     end do
C     end do
C     It seems G2 is not a standard normal-ordered RDM.
C       G2(standard)(a,b,c,d) = G2(MOLCAS)(a,c,b,d)
C     That is, in MOLCAS, G2(a,b,c,d) = <0|a+ c+ d b|0>
C     The indexing is, in a sense,  similar to E_{ab}E_{cd}:
C       E_{ab}E_{cd} = del(bc)*E_{ad} + E_{ac,bd}
C     In mkfg3.f, the del(bc)*E_{ad} term is subtracted, but
C     the index was not permutated: E_{ab}E_{cd} -> E_{ac,bd}
C     (see after "Correction to G2" in mkfg3.f).
C     So, in the standard index,
C       F1(iT,iU) = \sum_iV G2(iT,iV,iU,iV) * EPSA(iV)
C                                 ~~~~~
C     but in MOLCAS implmenetation, it has to be written (in the code):
C       F1(iT,iU) = \sum_iV G2(iT,iU,iV,iV) * EPSA(iV)
C                                 ~~~~~
C     Indexing of E_{ab}E_{cd} is normal.
C
C     G3 (and F3) is also similar.
C     G3(standard)(t,u,v,x,y,z) = G3(MOLCAS)(t,y,v,x,u,z)
C                               = G3(MOLCAS)(t,v,y,u,x,z)
C     Also, in the standard normal-ordered RDM,
C     EtuEvxEyz = del(uv)del(xy)G(tz)
C               + del(uv)G(ty,xz)
C               + del(xy)G(tv,uz)
C               + del(uy)G(vt,xz)
C               + G(tyv,xuz) or G(tvy,uxz)
C
C     In the mkfg3.f loop,
C     EtuEvxEyzww*e(w) is constructed, and then the following
C     corrections are applied:
C     EtuEvxEyz,ww
C     = t+ u v+ x y+ w+ w z f_{ww}
C     = del(ww) t+ u v+ x y+ z * f_{ww} - t+ u v+ x y+ w w+ z * f_{ww}
C     = F_{tvy,uxz}
C      + del(uv)del(xy)F_{tz} + del(uv)F_{ty,xz} + del(xy)F_{tv,uz} + del(uy)F_{tv,zx}
C      + del(uv)del(xw)G_{ty,wz}f_{ww} + del(xy)del(uw)G_{tv,wz}f_{ww} + del(uy)del(xw)F_{tv,zw}f_{ww}
C      + del(xw)G_{tvy,uwz}f_{uw} + del(uw)G_{tvy,wxz}f_{ww}
C     = (in canonical EPS) F_{tvy,uxz}
C      + del(uv)del(xy)F_{tz}
C      + del(uv)F_{ty,xz}    + del(xy)F_{tv,uz}    + del(uy)F_{tv,zx}
C      + del(uv)G_{ty,xz}e_x + del(xy)G_{tv,uz}e_u + del(uy)F_{tv,zx}e_x
C      + G_{tvy,uxz}(e_x+e_u)
C
C     For derivatives, we need an expression of non-canonical EPS:
C     EtuEvxEyz,ab*fab (a,b=active)
C     = t+ u v+ x y+ a+ b z * f_{ab}
C     = del(ab) t+ u v+ x y+ z - t+ u v+ x y+ b a+ z * f_{ab}
C     = del(uv)del(xy)F_{tz}
C      + del(uv)F_{ty,xz}        + del(xy)F_{tv,uz}        + del(uy)F_{tv,zx}
C      +(del(uv)del(xa)G_{ty,bz} + del(xy)del(ua)G_{tv,bz} + del(uy)del(xa)G_{tv,zb})*f_{ab}
C      +(del(xa)G_{tvy,ubz} + del(ua)G_{tvy,bxz})*f_{ab} + F_{tvy,uxz}
C
C
C
C     In the mkfg3.f loop,
C     EtuEwv,xwEyz*fww is constructed, and the following corrections are applied:
C     EtuEwv,xwEyz*fww
C     = t+ u w+ v+ w x y+ z f_{ww}
C     = del(vw) t+ u w+ x y+ z fww - t+ u w+ w v+ x y+ z fww
C     = t+ u v+ x y+ z fvv - t+ u w+ w v+ x y+ z fww
C     In non-canonical:
C     EtuEav,xbEyz * fab
C     = t+ u a+ v+ b x y+ z fab
C     = del(ua) t+ v+ b x y+ z fab - t+ a+ u v+ b x y+ z fab
C     = del(ua) del(xy) t+ v+ b z fab - del(ua) t+ v+ b y+ x z fab
C      -del(uv) t+ a+ b x y+ z fab + t+ a+ v+ u b x y+ z fab
C     = del(ua) del(xy) Gtv,zb fab - del(ua) del(by) t+ v+ x z fab + del(ua) t+ v+ y+ b x z fab
C      -del(uv) del(xy) t+ a+ b z fab + del(uv) t+ a+ b y+ x z fab
C      +del(xy) t+ a+ v+ u b z fab - t+ a+ v+ u b y+ x z fab
C     = del(ua) del(xy) Gtv,zb fab - del(ua) del(by) Gtv,zx fab + del(ua) Gtvy,zxb fab
C      -del(uv) del(xy) Ftz + del(uv) del(by) t+ a+ x z fab - del(uv) t+ a+ y+ b x z fab
C      +del(xy) Ftv,zu - del(by) t+ a+ v+ u x z fab + t+ a+ v+ u y+ b x z fab
C     = del(ua) del(xy) Gtv,zb fab - del(ua) del(by) Gtv,zx fab + del(ua) Gtvy,zxb fab
C      -del(uv) del(xy) Ftz + del(uv) del(by) Gta,zx fab + del(uv) Fty,zx
C      +del(xy) Ftv,zu - del(by) Gtav,zxu fab + del(uy) t+ a+ v+ b x z fab - t+ a+ v+ y+ u b x z fab
C     = del(ua) del(xy) Gtv,zb fab - del(ua) del(by) Gtv,zx fab + del(ua) Gtvy,zxb fab
C      -del(uv) del(xy) Ftz + del(uv) del(by) Gta,zx fab + del(uv) Fty,zx
C      +del(xy) Ftv,zu - del(by) Gtav,zxu fab - del(uy) Ftv,zx + Ftvy,zxu
C     =-del(ua) del(xy) Gtv,bz fab - del(ua) del(by) Gtv,zx fab - del(ua) Gtvy,bxz fab
C      -del(uv) del(xy) Ftz - del(uv) del(by) Gta,xz fab - del(uv) Fty,xz
C      -del(xy) Ftv,uz - del(by) Gtva,uxz fab - del(uy) Fvt,xz - Ftvy,uxz

C     =-del(ua) del(xy) G2(tbvz) fab - del(ua) del(by) G2(vxtz) fab + del(ua) G3(tuvxyz) fab
C      -del(uv) del(xy) F1(tz) - del(uv) del(by) G2(txaz) fab - del(uv) F2(txyz)
C      -del(xy) F2(tuvz) - del(by) G3(tuvxaz) fab - del(uy) F2(vxtz) - F3(tuvxyz)
C
C
C     --- Memo for RASPT2 Derivative ---
C
C     With RAS, EtuEvx is somehow wrong if both are true
C      - v>x and v,x are in a different RAS space
C      - t<u and t,u are in a different RAS space
C     Seeing the construction RDM2, at least v<=x has to be.
C
C     For instance, E13E41 can be evaluated -> OK for left
C     For right, E14E31 is needed, because <0|E13E41|I> = <I|E14E31|0>,
C     and this cannot be computed in a usual way. However,
C     <I|E14E31|0> = <I|E31E14|0>, so use this relation for RAS.
C
C     This applies to DG2 and DF2. DF1 has no problem, even though this
C     term is computed similarly. It is because the indices of the
C     latter excitation operator are the same.
C
      CALL QENTER('DERFG3')
C     call dcopy_(nlev**4,[0.0d+00],0,df2,1)
C     CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
C
C Put in zeroes. Recognize special cases:
      IF(nlev.EQ.0) GOTO 999

      IF(NACTEL.EQ.0) GOTO 999

      NCI=NCSF(LSYM)
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

      !! F2_{tuvw}
      !! = <0|E_{tu,vw,xy}|0>*f_{xy}
      !! = <0|t+ v+ x+ y w u|0> * f_{xy}
      !!   t+ v+ x+ y w u
      !! = t+ x+ v+ w y u
      !! = t+ x+ (del(vw) - w v+) y u
      !! = del(vw) t+ x+ y u - t+ x+ w v+ y u
      !! = del(vw)*F1(tu) - t+ (del(xw)-w x+)*(del(vy)-y v+) u
      !! = del(vw)*F1(tu) - del(xw)*del(vy)*t+ u
      !!   + del(xw)* t+ y v+ u + del(vy)*t+ w x+ u - t+ w x+ y v+ u
      !! = del(vw)*F1(tu) - del(xw)*del(vy)*G1_{tu}*f_{xy}
      !!   + del(xw)*<0|E_{ty}E_{vu}|0>*f_{xy} + del(vy)*<0|E_{tw}E_{xu}|0>*f_{xy}
      !!   - <0|E_{tw}E_{xy}E_{vu}|0>*f_{xy}
      !! (if canonical)
      !! = F1(tu) - G1{tu}*f_{wv} + <0|E_{tw}E_{vu}|0*f_{ww}
      !!   + <0|E_{tw}E_{vu}|0>*f_{vv} - <0|E_{tw}tE_{xy}E_{vu}|0>
      !!
      !!   E_{tu}E_{vw}E_{yz}*f_{vw}
      !! = t+ u v+ w y+ z * f_{vw}
      !! = t+ (del(uv) - v+ u)*(del(wy) - y+ w) z * f_{vw}
      !! = del(uv)*del(wy)*t+ z f_{vw}
      !!   - del(uv)* t+ y+ w z * f_{vw}
      !!   - del(wy)* t+ v+ u z * f_{vw}
      !!   + t+ v+ u y+ w z * f_{vw}
      !! = del(uv)*del(wy)*G1(tz)*f_{vw}
      !!   - del(uv)*E_{ty,zw}*f_{vw}
      !!   - del(wy)*E_{tv,zu}*f_{vw}
      !!   + t+ v+ (del(uy)-y+ u) w z * f_{vw}
      !! = del(uv)*del(wy)*G1(tz)*f_{vw}
      !!   - del(uv)*E_{ty,zw}*f_{vw}
      !!   - del(wy)*E_{tv,zu}*f_{vw}
      !!   + del(uy) t+ v+ w z * f_{vw}
      !!   - t+ v+ y+ u w z * f_{vw}
      !! = del(uv)*del(wy)*G1(tz)*f_{vw}
      !!   - del(uv)*E_{ty,zw}*f_{vw}
      !!   - del(wy)*E_{tv,zu}*f_{vw}
      !!   + del(uy)*F1_{tz}
      !!   + t+ y+ v+ w z u * f_{vw}
      !! = del(uv)*del(wy)*G1(tz)*f_{vw}
      !!   - del(uv)*E_{ty,zw}*f_{vw}
      !!   - del(wy)*E_{tv,zu}*f_{vw}
      !!   + del(uy)*F1_{tz}
      !!   + F2_{tuyz}

      !! F2_{tuvw}
      !! = <0|E_{tu,vw,xy}|0>*f_{xy}
      !! = <0|t+ v+ x+ y w u|0> * f_{xy}
      !!   t+ v+ x+ y w u
      !! =-t+ v+ x+ y u w
      !! =-t+ x+ v+ u y w * f_{xy}
      !! =-t+ x+ (del(vu) - u v+) y w * f_{xy}
      !! =-del(vu) E_{tx,wy}*f_{xy} + t+ x+ u v+ y w * f_{xy}
      !! =-del(vu) E_{tx,wy}*f_{xy}
      !!  + t+ (del(ux) - u x+) (del(vy) - y v+) w * f_{xy}
      !! =-del(vu) E_{tx,wy}*f_{xy}
      !!  + del(ux)del(vy) G(t,w)*f_{xy} - del(ux) E_{ty}E_{vw}*f_{xy}
      !!  - del(vy) E_{tu}E_{xw}*f_{xy} + E_{tu}E_{xy}E_{vw}*f_{xy}
      !! = -del(vu) E_{tx,wy}*f_{xy} + del(ux)del(vy) G1(t,w)*f_{xy}
      !!   -del(ux) E_{ty}E_{vw}*f_{xy}

      !! = G1(tz)*f_{uy}
      !!   + E_{ty,wz}*f_{uw} (in the code, E(it,iw(iu),iy,iz))*e(u))
      !!   + E_{tv,uz}*f_{vy} (E(it,iu,iv(iy),iz)*e(y))
      !!   + del(uy)*F1_{tz}
      !!   + F2_{tuyz}
* Correction to F2: It is now = <0| E_tu H0Diag E_yz |0>
      !! Etu H0Diag Eyz
      !! = t+ u w+ w y+ z fww
      !! = del(uw) t+ w y+ z fww - t+ w+ u w y+ z fww
      !! = t+ u y+ z fuu - del(wy) t+ w+ u z fww + t+ w+ u y+ w z fww
      !! = t+ u y+ z fuu - t+ y+ u z fyy + del(uy) t+ w+ w z fww  - t+ w+ y+ u w z fww
      !! = t+ u y+ z fuu - t+ y+ u z fyy + del(uy) Ftz - t+ y+ w+ w u z fww
      !! = del(uy) Gtz fuu - Gty,zu fuu - Gty,zu fyy + del(uy) Ftz - Fty,zu
      !! = del(uy) (Ftz+Gtz fuu) + Gty,uz (fuu+fyy) + Fty,uz
      !!-> del(uy) (F1(tz)+G1(tz) fuu) + G2(tuyz) (fuu+fyy) + F2(tuyz)
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
         If (iy.ne.it.and.iz.ne.iu) Then
           DF2(iy,iz,it,iu)=DF2(iy,iz,it,iu)+DF2(it,iu,iy,iz)
           DF2(it,iu,iy,iz)=0.0d+00
         End If
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
         If (iz.ne.it.and.iy.ne.iu) Then
          DF2(iz,iy,iu,it)=DF2(iz,iy,iu,it)+DF2(it,iu,iy,iz)
          DF2(it,iu,iy,iz)=0.0d+00
         End If
        end do
       end do
C
       do iz=1,nlev
        do iy=1,nlev
         do iu=1,nlev
          do it=1,nlev
           if (dodebug) then
           DG2(it,iu,iy,iz) = DG2(it,iu,iy,iz)
     *      - DF2(it,iu,iy,iz)*(EPSA(iu)+EPSA(iz))
           else
        if (ras.and.it.lt.iu.and.iy.gt.iz) then
           !! why is this if-clause needed?
           DG2(iy,iz,it,iu) = DG2(iy,iz,it,iu)
     *      - DF2(it,iu,iy,iz)*(EPSA(iz)+EPSA(it))
        else
           DG2(it,iu,iy,iz) = DG2(it,iu,iy,iz)
     *      - DF2(it,iu,iy,iz)*(EPSA(iu)+EPSA(iy))
        end if
           end if
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
* Correction to G2: It is now = <0| E_tu E_yz |0>
      do it=1,nlev
       do iu=1,nlev
        itu = it+nasht*(iu-1)
        do iy=1,nlev
         do iz=1,nlev
          iyz = iy+nlev*(iz-1)
C         if (iyz.gt.itu) then
C           g2(it,iu,iv,ix) = g2(it,iu,iv,ix) + g2(iv,ix,it,iu)
C           g2(iv,ix,it,iu) = 0.0d+00
C           f2(it,iu,iv,ix) = f2(it,iu,iv,ix) + f2(iv,ix,it,iu)
C           f2(iv,ix,it,iu) = 0.0d+00
C         end if
         end do
        end do
       end do
      end do
C     call dcopy_(nlev**4,dg2,1,g2,1)
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
        If (iy.ne.it.and.iz.ne.iu) Then
        DG2(iy,iz,it,iu)=DG2(iy,iz,it,iu)+DG2(it,iu,iy,iz)
        DG2(it,iu,iy,iz)=0.0D+00
        End If
       end do
      end do
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
        If (iz.ne.it.and.iy.ne.iu) Then
        DG2(iz,iy,iu,it)=DG2(iz,iy,iu,it)+DG2(it,iu,iy,iz)
        DG2(it,iu,iy,iz)=0.0d+00
        End If
       end do
      end do
      do iu=1,nlev
       do iz=1,nlev
        do it=1,nlev
         DG1(it,iz) = DG1(it,iz) - DG2(it,iu,iu,iz)
        end do
       end do
      end do
* Additional terms for F1
      Do iT = 1, NLEV
        Do iU = 1, NLEV
          DG1(iT,iU) = DG1(iT,iU) - DF1(iT,iU)*EPSA(iU)
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
      Call CLagSym(nLev,DG1,DG2,DF1,DF2,1)
C
* Dummy values necessary for fooling syntax checkers:
      ldum=1
      ndum=1
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
      nbuf2= 1
      nbuf3= 1
      nbuf4= 1
      nbuf5= 1
      nbuf6= 1
      nbuft= 1
      nbufd= 1
      nbufe= 1
      nbuf1=max(1,min(nlev2,(memmax_safe-8*mxci)/mxci/2))
      nbufa=max(1,min(nlev2,(memmax_safe-8*mxci)/mxci/2))
      if (nbufa.ne.nlev2) then
        write(6,*) "nlev2.ne.nbufa..."
        write(6,*) "need more memory?"
        call abend()
      end if
      CALL GETMEM('BUF1','ALLO','REAL',LBUF1,NBUF1*MXCI)
      CALL GETMEM('BUF2','ALLO','REAL',LBUF2,NBUF2*MXCI)
      CALL GETMEM('BUF3','ALLO','REAL',LBUF3,NBUF3*MXCI)
      CALL GETMEM('BUF4','ALLO','REAL',LBUF4,NBUF4*MXCI)
      CALL GETMEM('BUF5','ALLO','REAL',LBUF5,NBUF5*MXCI)
      CALL GETMEM('BUF6','ALLO','REAL',LBUF6,NBUF6*MXCI)
      CALL GETMEM('BUFT','ALLO','REAL',LBUFT,NBUFT*MXCI)
      CALL GETMEM('BUFD','ALLO','REAL',LBUFD,NBUFD*MXCI)
      CALL GETMEM('BUFE','ALLO','REAL',LBUFE,NBUFE*MXCI)
      CALL GETMEM('BUFA','ALLO','REAL',LBUFA,NBUFA*MXCI)

      !! For the derivative of the Fock-contracted portion
C     CALL GETMEM('FCDER1','ALLO','REAL',LFCDer1,NBUFD*MXCI)
C     Call DCopy_(NBUFD*MXCI,[0.0D+00],0,Work(LFCDer1),1)
C     CALL GETMEM('FCDER2','ALLO','REAL',LFCDer2,NLEV*NLEV)
C     Call DCopy_(NLEV*NLEV ,[0.0D+00],0,Work(LFCDer2),1)

C-SVC20100301: calculate maximum number of tasks possible
      MXTASK=(NTRI2-1)/NBUF1+1+(NTRI1-1)/NBUF1+1
      CALL GETMEM ('TASKLIST','ALLO','INTE',lTask_List,4*mxTask)
      lip1sta=lTask_List
      lip1end=lTask_List+mxTask
      libuf1=lTask_List+2*mxTask
      lOffSet=lTask_List+3*mxTask

      IF(iPrGlb.GE.VERBOSE) THEN
        WRITE(6,*)
        WRITE(6,'(2X,A)') 'Constructing G3/F3'
        WRITE(6,'(2X,A,F16.9,A)') ' memory avail: ',
     &    (memmax*RtoB)/1.0D9, ' GB'
        WRITE(6,'(2X,A,F16.9,A)') ' memory used:  ',
     &    (((nbuf1+3)*MXCI)*RtoB)/1.0D9, ' GB'
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
       isp1=mul(issg1,lsym)
       nsgm1=ncsf(issg1)
       !! EASUM = sum_I Work(LBUFD)*CI(I)*CI(I)
       !! (LHS)
       !!   \sum_{tt} <0|E_{tt}|0>*f_{tt}
       !! = \sum_{IJ}\sum_{tt} CI(I)*CI(J)*<I|E_{tt}|J>*f_{tt}
       !! = \sum_{I} \sum_{tt} CI(I)*CI(I)*<I|E_{tt}|I>*f_{tt}
       !! It's because E_{tt} is valid only when I=J
       !!
       !! Work(LBufD) = \sum_t <I|E_{tt}|I>*f_{tt}
       CALL H0DIAG_CASPT2(ISSG1,WORK(LBUFD),IWORK(LNOW),IWORK(LIOW))
       Call DCopy_(MXCI*NBUFA,[0.0D+00],0,Work(LBUFA),1)
C      tmp = 0.0d+00
C      tmp2=0.0d+00
C      do i = 1, nsgm1
C        tmp = tmp + work(lbufd+i-1)*ci(i)
C        tmp2= tmp2+ work(lbufd+i-1)*ci(i)*ci(i)
C       end do
C     write(6,*) "before"
C      do i = 1, nsgm1
C        write(6,'(i3,f20.10)') i,work(lbufd+i-1)
C      end do

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


      if (.false.) then
        !! asdfasdf
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
          CALL SIGMA1_CP2(1,2,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C     do i = 1, nsgm1
C     write(6,'(i3,f20.10)') i,work(lbuf1+i-1)
C     end do
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
          CALL SIGMA1_CP2(2,1,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C     do i = 1, nsgm1
C     write(6,'(i3,f20.10)') i,work(lbuf1+i-1)
C     end do
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
          CALL SIGMA1_CP2(1,2,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVR),IWORK(LMVL))
C     do i = 1, nsgm1
C     write(6,'(i3,f20.10)') i,work(lbuf1+i-1)
C     end do
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
          CALL SIGMA1_CP2(2,1,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVR),IWORK(LMVL))
C     do i = 1, nsgm1
C     write(6,'(i3,f20.10)') i,work(lbuf1+i-1)
C     end do
        write(6,*) "Etu"
        do it = 1, nlev
        do iu = 1, nlev
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
          CALL SIGMA1_CP2(IT,IU,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          scal=ddot_(nsgm1,ci,1,work(lbuf1),1)
          write(6,'(2i3,f20.10)') it,iu,scal
        end do
        end do
        call sqprt(g1,nlev)
        write(6,*) "EtuEvx"
        do ix = 1, nlev
        do iy = 1, nlev
          do it = 1, nlev
          do iu = 1, nlev
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
          CALL SIGMA1_CP2(IX,IY,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
            CALL SIGMA1_CP2(IT,IU,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
     &           IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &           IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &           WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          scal=ddot_(nsgm1,ci,1,work(lbuf2),1)
          tmp=scal
          if (iu.eq.ix) then
            tmp=tmp-g1(it,iy)
          end if
          tmp2=tmp-g2(it,iu,ix,iy)
          if (it.lt.iu.and.ix.gt.iy) then
            !!  t+ u x+ y - del(ux) t+ y
            ! = x+ y t+ u - del(yt) x+ u
            !! left derivative of t+ u x+ y :: <I|(tu)(xy)|CI1>
            !! right derivative of t+ u x+ y ::
            !!   <CI1|(tu)(xy)|I> = <I|(yx)(ut)|CI1>
            !   = <I|(ut)|(yx)|CI1>
            !! Similar permutation applies to F2 too
            call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
            CALL SIGMA1_CP2(IT,IU,1.0D00,LSYM,CI,Work(LBUF1),
     &           IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &           IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &           WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
            call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
            CALL SIGMA1_CP2(IX,IY,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
     &           IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &           IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &           WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
            scal=ddot_(nsgm1,ci,1,work(lbuf2),1)
            tmp=scal
            if (iy.eq.it) then
              tmp=tmp-g1(ix,iu)
            end if
            tmp2=tmp-g2(it,iu,ix,iy)
          end if
      if (abs(tmp2).ge.1.0d-10)
     *    write(6,'(4i3,4f20.10)') it,iu,ix,iy,scal,tmp,
     *      g2(it,iu,ix,iy),tmp2
          end do
          end do
        end do
        end do
C
        write(6,*) "EtuEvvEyz*ev"
        do ix = 1, nlev
        do iy = 1, nlev
          do it = 1, nlev
          do iu = 1, nlev
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
          CALL SIGMA1_CP2(IX,IY,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          do icsf = 1, nsgm1
            work(lbuf1+icsf-1) = work(lbuf1+icsf-1)*work(lbufd+icsf-1)
          end do
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
            CALL SIGMA1_CP2(IT,IU,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
     &           IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &           IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &           WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          scal=ddot_(nsgm1,ci,1,work(lbuf2),1)
          tmp=scal
          tmp = tmp - (epsa(iu)+epsa(ix))*g2(it,iu,ix,iy)
          if (iu.eq.ix) then
            tmp=tmp-f1(it,iy)-epsa(iu)*g1(it,iy)
          end if
          tmp2=tmp-f2(it,iu,ix,iy)
          if (it.lt.iu.and.ix.gt.iy) then
            call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
            CALL SIGMA1_CP2(IT,IU,1.0D00,LSYM,CI,Work(LBUF1),
     &           IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &           IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &           WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
            do icsf = 1, nsgm1
              work(lbuf1+icsf-1) = work(lbuf1+icsf-1)*work(lbufd+icsf-1)
            end do
            call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
            CALL SIGMA1_CP2(IX,IY,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
     &           IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &           IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &           WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
            scal=ddot_(nsgm1,ci,1,work(lbuf2),1)
            tmp=scal
            tmp = tmp - (epsa(iy)+epsa(it))*g2(ix,iy,it,iu)
            if (iy.eq.it) then
              tmp=tmp-f1(ix,iu)-epsa(iy)*g1(ix,iu)
            end if
            tmp2=tmp-f2(it,iu,ix,iy)
          end if
      if (abs(tmp2).ge.1.0d-10)
     *    write(6,'(4i3,4f20.10)') it,iu,ix,iy,scal,tmp,
     *      f2(it,iu,ix,iy),tmp2
          end do
          end do
        end do
        end do
C
        !! G3
        write(6,*) "g3 loop"
        do ig3 = 1, ng3
        iT=idxG3(1,iG3)
        iU=idxG3(2,iG3)
        iV=idxG3(3,iG3)
        iX=idxG3(4,iG3)
        iY=idxG3(5,iG3)
        iZ=idxG3(6,iG3)
        if (it.eq.0) exit
        call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
        CALL SIGMA1_CP2(IY,IZ,1.0D+00,LSYM,CI,Work(LBUF1),
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
        CALL SIGMA1_CP2(IV,IX,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        call dcopy_(nsgm1,[0.0d+00],0,work(lbuf3),1)
        CALL SIGMA1_CP2(IT,IU,1.0d+00,LSYM,Work(LBUF2),Work(LBUF3),
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))

        call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
        CALL SIGMA1_CP2(IU,IT,1.0D+00,LSYM,CI,Work(LBUF1),
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
        CALL SIGMA1_CP2(IX,IV,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        call dcopy_(nsgm1,[0.0d+00],0,work(lbuf3),1)
        CALL SIGMA1_CP2(IZ,IY,1.0d+00,LSYM,Work(LBUF2),Work(LBUF3),
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        scal=ddot_(nsgm1,ci,1,work(lbuf3),1)
        tmp = scal
       if(iY.eq.iX) then
        tmp=tmp-G2(iT,iU,iV,iZ)
        if(iv.eq.iu) then
         tmp=tmp-G1(iT,iZ)
        end if
       end if
       if(iV.eq.iU) then
         tmp=tmp-G2(iT,iX,iY,iZ)
       end if
       if(iY.eq.iU) then
         tmp=tmp-G2(iV,iX,iT,iZ)
       end if
            tmp2=tmp-g3(ig3)
      if (abs(tmp2).ge.1.0d-10)
     *  write(6,'(6i3,3f20.10)') it,iu,iv,ix,iy,iz,tmp,
     * g3(ig3),tmp2
        end do
C       call sqprt(g2,nlev**2)
C       write(6,*) "mktg3"
C       call mktg33(ci)
        call abend
        do it = 1, nlev
        do iu = 1, nlev
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
          CALL SIGMA1_CP2(IT,IU,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          do ix = 1, nlev
          do iy = 1, nlev
            scal = dg2(ix,iy,it,iu)+dg2(iu,it,iy,ix)
            CALL SIGMA1_CP2(IX,IY,SCAL,LSYM,Work(LBUF1),CLAG,
     &           IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &           IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &           WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          end do
          end do
        end do
        end do
        call dcopy_(nlev**4,[0.0d+00],0,dg2,1)
      end if

C-SVC20100301: initialize the series of subtasks
      Call Init_Tsk(ID, nSubTasks)

      myBuffer=0
      if (dodebug) goto 1234
      if (.false.) Then
        do it = 1, nlev
        do iu = 1, nlev
        do ix = 1, nlev
        do iy = 1, nlev
          scal = dg2(it,iu,ix,iy)
          !! left derivative (EtuExy)
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
          if (it.lt.iu.and.ix.gt.iy) then
          CALL SIGMA1_CP2(IT,IU,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          CALL SIGMA1_CP2(IX,IY,SCAL,LSYM,Work(LBUF1),CLAG,
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          else
          CALL SIGMA1_CP2(IX,IY,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          CALL SIGMA1_CP2(IT,IU,SCAL,LSYM,Work(LBUF1),CLAG,
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          end if
          !! right derivative (EyxEut)
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
          if (iy.lt.ix.and.iu.gt.it) then
          CALL SIGMA1_CP2(IY,IX,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          CALL SIGMA1_CP2(IU,IT,SCAL,LSYM,Work(LBUF1),CLAG,
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          else
          CALL SIGMA1_CP2(IU,IT,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          CALL SIGMA1_CP2(IY,IX,SCAL,LSYM,Work(LBUF1),CLAG,
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          end if
C
C
C
          scal = df2(it,iu,ix,iy)
          !! left derivative (EtuExy)
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
          if (it.lt.iu.and.ix.gt.iy) then
          CALL SIGMA1_CP2(IT,IU,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          do icsf = 1, nsgm1
            work(lbuf1+icsf-1) = work(lbuf1+icsf-1)*work(lbufd+icsf-1)
          end do
          CALL SIGMA1_CP2(IX,IY,SCAL,LSYM,Work(LBUF1),CLAG,
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          else
          CALL SIGMA1_CP2(IX,IY,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          do icsf = 1, nsgm1
            work(lbuf1+icsf-1) = work(lbuf1+icsf-1)*work(lbufd+icsf-1)
          end do
          CALL SIGMA1_CP2(IT,IU,SCAL,LSYM,Work(LBUF1),CLAG,
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          end if
          !! right derivative (EyxEut)
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
          if (iy.lt.ix.and.iu.gt.it) then
          CALL SIGMA1_CP2(IY,IX,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          do icsf = 1, nsgm1
            work(lbuf1+icsf-1) = work(lbuf1+icsf-1)*work(lbufd+icsf-1)
          end do
          CALL SIGMA1_CP2(IU,IT,SCAL,LSYM,Work(LBUF1),CLAG,
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          else
          CALL SIGMA1_CP2(IU,IT,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          do icsf = 1, nsgm1
            work(lbuf1+icsf-1) = work(lbuf1+icsf-1)*work(lbufd+icsf-1)
          end do
          CALL SIGMA1_CP2(IY,IX,SCAL,LSYM,Work(LBUF1),CLAG,
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          end if
          ! depsa
          do iv = 1, nlev
          do iz = 1, nlev
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf3),1)
          if (it.lt.iu.and.ix.gt.iy) then
          CALL SIGMA1_CP2(IT,IU,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          CALL SIGMA1_CP2(IV,IZ,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          CALL SIGMA1_CP2(IX,IY,1.0D+00,LSYM,Work(LBUF2),Work(LBUF3),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          else
          CALL SIGMA1_CP2(IX,IY,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          CALL SIGMA1_CP2(IV,IZ,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          CALL SIGMA1_CP2(IT,IU,1.0D+00,LSYM,Work(LBUF2),Work(LBUF3),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          end if
          depsa(iv,iz) = depsa(iv,iz)
     *      + scal*ddot_(nsgm1,ci,1,work(lbuf3),1)
          end do
          end do
        end do
        end do
        end do
        end do
        call dcopy_(nlev**4,[0.0d+00],0,dg2,1)
        call dcopy_(nlev**4,[0.0d+00],0,df2,1)
      end if

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
        !! ip1sta=1, ip1end=nact*(nact+1)/2
        !! 1 2
        !! 1 3
        !! ...
        !! 1 1
        !! 2 2
        !! ...
        !! Work(LTO) :: TO(nconf,nact*(nact+1)/2) array
C       write(6,*) "ip1sta,ip1end=",ip1sta,ip1end
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
C        write(6,'(6i3)') ip1i,itlev,iulev,it,iu,ibuf1
          call dcopy_(nsgm1,[0.0D0],0,work(lto),1)
      !! According to mktg3.f, it seems that SIGMA1_CP2 computes
      !!   TO(I,t-u) = \sum_J CI(J) * <J|E_{tu}|I>
      !!             = <0|E_{tu}|I>
      !! it seems that SIGMA1_CP2 computes
      !!   <I|E_{tu}|0>
      !! 1.0D+00 is the scaling factor.
      !! Work(LTO) has to be initialized
      !! In the following, SIGMA1_CP2, "iULEV,ITLEV",
      !! so it is <I|E_{ut}|0> = <0|E_{tu}|I>
      !! Note that <0|E_{tu}|I> \= <0|E_{ut}|I>

      ! <0|E_{tu}|I>*CI(I) = <0|E_{ut}|I>*CI(I)

      !! The right derivative can be computed easily
      !! Suppose D_{tu}*<0|E_{tu}|I>*dCI(I)/da = CLAGR(I)*dCI(I)/da
      !! where CLAGR(I) = D_{tu}*<0|E_{tu}|I>
      !! It is just a DAXPY operation
      !! The left derivative is not straightforward.
      !! CLAGL(I) = D_{tu}*<I|E_{tu}|0>
      !!          = D_{tu}*<I|E_{tu}|J>*CI(J)
      !!          = D_{tu}*<J|E_{ut}|I>*CI(J)
      !!          = D_{ut}*<J|E_{tu}|I>*CI(J)
      !! In the subroutine SGM is alternated with CLAG,
      !! and CLAG(I) = D_{ut} * CI(J)*<J|E_{tu}|I>
          CALL SIGMA1_CP2(IULEV,ITLEV,1.0D00,LSYM,CI,WORK(LTO),
     &     IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &     IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &     WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C         CALL SIGMA1_CP2ld(IULEV,ITLEV,1.0D00,LSYM,CI,WORK(LTO),
C    &     IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &     IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &     WORK(LVTAB),IWORK(LMVL),IWORK(LMVR),dg1(it,iu),clag)
         end if
        end do
        myBuffer=iTask
      ELSE
        ibuf1=iWork(libuf1+iTask-1)
      ENDIF
C-SVC20100301: necessary batch of sigma vectors is now in the buffer

      ! The ip1 buffer could be the same on different processes
      ! so only compute the G1 contribution when ip3 is 1, as
      ! this will only be one task per buffer.
      if (issg1.eq.lsym.AND.ip3.eq.1) then
      ! lbuf1 = <Psi0|E_ip1|I> ?
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
          SCAL = DG1(iT,iU) + DG1(iT,iU)
          Call DaXpY_(nsgm1,SCAL,Work(lto),1,CLag,1) ! asdf
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
      issg2=mul(isyz,lsym)
      nsgm2=ncsf(issg2)
      iy=L2ACT(iylev)
      iz=L2ACT(izlev)
      lto=lbuf2
      call dcopy_(nsgm2,[0.0D0],0,work(lto),1)
      CALL SIGMA1_CP2(IYLEV,IZLEV,1.0D00,LSYM,CI,WORK(LTO),
     &     IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &     IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &     WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
      call dcopy_(nsgm2,[0.0D0],0,work(lbuf3),1)
      If (IYLEV.NE.IZLEV) Then
        Call GETSGM2(IZLEV,IYLEV,LSYM,CI,Work(LBUF4))
        call dcopy_(nsgm2,[0.0D0],0,work(lbuf5),1)
      End If
      Call DCopy_(nsgm1,[0.0D+00],0,Work(LBUFE),1)
      if(issg2.eq.issg1) then
        do ib=1,ibuf1
          idx=ip1_buf(ib)
          itlev=idx2ij(1,idx)
          iulev=idx2ij(2,idx)
          it=L2ACT(itlev)
          iu=L2ACT(iulev)
C         write(6,'("0 ",4i1)') it,iu,iy,iz
          !! E_{tu}E_{yz} = E_{ut}E_{zy} = E_{yz}E_{tu} = E_{zy}E_{ut}
          !! However, this does not apply to the configuration derivative.
C
          !! CLAG(I) = <I|E_{tu}E_{yz}|0> + <0|E_{tu}E_{yz}|I>
          !! left derivative
          !! <I|E_{tu}E_{yz}|0>*D_{tu,yz}
          !! = <I|E_{tu}|J>*<J|E_{yz}|0>*D_{tu,yz}
          !! right derivative
          !! <0|E_{tu}E_{yz}|I>*D_{tu,yz}
          !! = <0|E_{tu}|J>*<J|E_{yz}|I>*D_{tu,yz}
          !! = <J|E_{ut}|0>*<I|E_{zy}|J>*D_{tu,yz}
          !! = <I|E_{zy}|J>*<J|E_{ut}|0>*D_{tu,yz}
C
          RASSPE = RAS.and.iT.LT.iU.and.iY.NE.iZ
C
          !! Set some parameters for left derivative
          !! EtuEyz
          ScalGYZL=DG2(iT,iU,iY,iZ)
          If (iY.eq.iZ) ScalGYZL=ScalGYZL+DF1(iT,iU)*EPSA(iY)
          Call DCopy_(nsgm1,[0.0D+00],0,Work(LBUFT),1)
          If (ScalGYZL.ne.0.0D+00) Then
            DG2(iT,iU,iY,iZ) = 0.0D+00
            If (.not.RASSPE.or.iY.LE.iZ) Then ! .not.RAS.or.iT.GE.iU.or.iY.LE.iZ) Then
              Call DaXpY_(nsgm1,ScalGYZL,Work(LTO),1,Work(LBUFT),1)
            End If
          End If
          ScalFYZL=DF2(iT,iU,iY,iZ)
          If (ScalFYZL.ne.0.0D+00) Then
            DF2(iT,iU,iY,iZ) = 0.0D+00
          End If
          !! EtuEzy
          ScalGZYL=DG2(iT,iU,iZ,iY)
          If (ScalGZYL.ne.0.0D+00) Then ! .and.IYLEV.NE.IZLEV) Then
            DG2(iT,iU,iZ,iY) = 0.0D+00
            If (.not.RASSPE.or.iZ.LE.iY) Then ! .or.iT.GE.iU.or.iZ.LE.iY) Then
              Call DaXpY_(nsgm1,ScalGZYL,Work(LBUF4),1,Work(LBUFT),1)
            End If
          End If
          ScalFZYL=DF2(iT,iU,iZ,iY)
          DF2(iT,iU,iZ,iY) = 0.0D+00
          !!
          If (ScalFYZL.ne.0.0D+00.or.ScalFZYL.ne.0.0D+00.or.RASSPE) Then
            If (RASSPE) Then
              Call GETSGM2(ITLEV,IULEV,LSYM,CI,Work(LBUFE))
              If (iY.GT.iZ) Then
                Do icsf = 1, nsgm1
                  Work(LBUFT+icsf-1) = Work(LBUFT+icsf-1)
     *              + ScalFZYL*Work(LBUF4+icsf-1)*Work(LBUFD+icsf-1)
                  Work(LBUFE+icsf-1) =
     *              + ScalGYZL*Work(LBUFE+icsf-1)
     *              + ScalFYZL*Work(LBUFE+icsf-1)*Work(LBUFD+icsf-1)
                End Do
              Else If (iZ.GT.iY) Then
                Do icsf = 1, nsgm1
                  Work(LBUFT+icsf-1) = Work(LBUFT+icsf-1)
     *              + ScalFYZL*Work(LTO+icsf-1)*Work(LBUFD+icsf-1)
                  Work(LBUFE+icsf-1) =
     *              + ScalGZYL*Work(LBUFE+icsf-1)
     *              + ScalFZYL*Work(LBUFE+icsf-1)*Work(LBUFD+icsf-1)
                End Do
              End If
            Else
              Do icsf = 1, nsgm1
                Work(LBUFT+icsf-1) = Work(LBUFT+icsf-1)
     *            + (ScalFYZL*Work(LTO+icsf-1)
     *              +ScalFZYL*Work(LBUF4+icsf-1))
     *            *Work(LBUFD+icsf-1)
              End Do
            End If
          End If
C
          If (iT.ne.iU) Then
            ScalGYZR = DG2(iU,iT,iY,iZ)
            If (iY.eq.iZ) ScalGYZR = ScalGYZR+DF1(iU,iT)*EPSA(iY)
            If (ScalGYZR.ne.0.0D+00) THen
              DG2(iU,iT,iY,iZ) = 0.0D+00
            End If
            SCALFYZR=DF2(iU,iT,iY,iZ)
            If (SCALFYZR.ne.0.0D+00) Then
              DF2(iU,iT,iY,iZ) = 0.0D+00
            End If
            SCALGZYR=DG2(iU,iT,iZ,iY)
            If (SCALGZYR.ne.0.0D+00) Then
              DG2(iU,iT,iZ,iY) = 0.0D+00
            End If
            SCALFZYR=DF2(iU,iT,iZ,iY)
            If (SCALFZYR.ne.0.0D+00) Then
              DF2(iU,iT,iZ,iY) = 0.0D+00
            End If
          Else
            ScalGYZR = 0.0D+00
            ScalFYZR = 0.0D+00
            ScalGZYR = 0.0D+00
            ScalFZYR = 0.0D+00
          End If
C
          !! Complete the left derivative
          CALL SIGMA1_CP2(ITLEV,IULEV,1.0D+00,LSYM,Work(LBUFT),CLAG,
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          If (RASSPE) Then
            If (iY.GT.iZ) Then
            CALL SIGMA1_CP2(IYLEV,IZLEV,1.0D+00,LSYM,Work(LBUFE),CLAG,
     &           IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &           IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &           WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
            Else If (iZ.GT.iY) Then
            CALL SIGMA1_CP2(IZLEV,IYLEV,1.0D+00,LSYM,Work(LBUFE),CLAG,
     &           IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &           IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &           WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
            End If
          End If
C
          !! right derivative of EtuEyz and EtuEzy
          !! First, <0|EtuEyz|I> = <I|EzyEut|0>
          !! For RAS (special): EzyEut = EutEzy
          RASSPE = RAS.and.iZ.NE.iY.and.iU.GT.iT
          If (ScalGYZL.NE.0.0D+00.OR.ScalGZYL.NE.0.0D+00.OR.
     *        ScalFYZL.NE.0.0D+00.OR.ScalFZYL.NE.0.0D+00) Then
            Call GETSGM2(IULEV,ITLEV,LSYM,CI,Work(LBUFT))
            If (.not.RASSPE.or.iZ.GE.iY) Then
              Call DaXpY_(nsgm1,ScalGYZL,WORK(LBUFT),1,Work(LBUF3),1)
              Do icsf = 1, nsgm1
                Work(LBUF3+icsf-1) = Work(LBUF3+icsf-1)
     *            + ScalFYZL*Work(LBUFT+icsf-1)*Work(LBUFD+icsf-1)
              End Do
            End If
            If (.not.RASSPE.or.iY.GE.iZ) Then
              Call DaXpY_(nsgm1,ScalGZYL,WORK(LBUFT),1,Work(LBUF5),1)
              Do icsf = 1, nsgm1
                Work(LBUF5+icsf-1) = Work(LBUF5+icsf-1)
     *            + ScalFZYL*Work(LBUFT+icsf-1)*Work(LBUFD+icsf-1)
              End Do
            End If
            If (RASSPE) Then
              If (iZ.LT.iY) Then
                Call GETSGM2(IZLEV,IYLEV,LSYM,CI,Work(LBUFT))
                Do icsf = 1, nsgm1
                  Work(LBUFT+icsf-1) =
     *              + ScalGYZL*Work(LBUFT+icsf-1)
     *              + ScalFYZL*Work(LBUFT+icsf-1)*Work(LBUFD+icsf-1)
                End Do
              Else If (iY.LT.iZ) Then
                Call GETSGM2(IYLEV,IZLEV,LSYM,CI,Work(LBUFT))
                Do icsf = 1, nsgm1
                  Work(LBUFT+icsf-1) =
     *              + ScalGZYL*Work(LBUFT+icsf-1)
     *              + ScalFZYL*Work(LBUFT+icsf-1)*Work(LBUFD+icsf-1)
                End Do
              End If
              CALL SIGMA1_CP2(IULEV,ITLEV,1.0D+00,LSYM,Work(LBUFT),CLAG,
     &             IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &             IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &             WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
            End If
C           !! For DEPSA
C           Call DCopy_(nsgm1,Work(LBUFT),1,Work(LBUF6),1)
          Else
C           Call dcopy_(nsgm1,[0.0d+00],0,Work(LBUF6),1)
          End If
C
          !! For EutEyz and EutEzy cases
          !! First, left derivative
          !! For special RAS, EutEyz -> EyzEut
          RASSPE = RAS.and.iU.LT.iT.and.iY.NE.iZ
          If (iT.ne.iU) Then
            Call DCopy_(nsgm1,[0.0D+00],0,Work(LBUFT),1)
            If (ScalGYZR.ne.0.0D+00) Then
              If (.not.RASSPE.or.iY.LE.iZ) Then
                Call DaXpY_(nsgm1,ScalGYZR,Work(LTO),1,Work(LBUFT),1)
              End If
            End If
            If (SCALGZYR.ne.0.0D+00) Then
              If (.not.RASSPE.or.iZ.LE.iY) Then
                Call DaXpY_(nsgm1,ScalGZYR,Work(LBUF4),1,Work(LBUFT),1)
              End If
            End If
            If (RASSPE) Then
              Call GETSGM2(IULEV,ITLEV,LSYM,CI,Work(LBUFE))
              If (iY.GT.iZ) Then
                Do icsf = 1, nsgm1
                 Work(LBUFE+icsf-1) =
     *             + ScalGYZR*Work(LBUFE+icsf-1)
     *             + ScalFYZR*Work(LBUFE+icsf-1)*Work(LBUFD+icsf-1)
                End Do
              CALL SIGMA1_CP2(IYLEV,IZLEV,1.0D+00,LSYM,Work(LBUFE),CLAG,
     &             IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &             IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &             WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
              Else If (iZ.GT.iY) Then
                Do icsf = 1, nsgm1
                 Work(LBUFE+icsf-1) =
     *             + ScalGZYR*Work(LBUFE+icsf-1)
     *             + ScalFZYR*Work(LBUFE+icsf-1)*Work(LBUFD+icsf-1)
                End Do
              CALL SIGMA1_CP2(IZLEV,IYLEV,1.0D+00,LSYM,Work(LBUFE),CLAG,
     &             IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &             IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &             WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
              End If
            Else
              If (ScalFYZR.ne.0.0D+00.or.ScalFZYR.ne.0.0D+00) Then
                Do icsf = 1, nsgm1
                 Work(LBUFT+icsf-1) = Work(LBUFT+icsf-1)
     *           + (ScalFYZR*Work(LTO+icsf-1)
     *             +ScalFZYR*Work(LBUF4+icsf-1))
     *             *Work(LBUFD+icsf-1)
                End Do
              End If
            End If
            CALL SIGMA1_CP2(IULEV,ITLEV,1.0D+00,LSYM,Work(LBUFT),CLAG,
     &           IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &           IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &           WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
            !! right derivative of EutEyz and EutEzy
            !! <0|EutEyz|I> = <I|EzyEtu|0>
            !! For special RAS, EzyEtu = EtuEzy
            If (ScalGYZR.NE.0.0D+00.OR.ScalGZYR.NE.0.0D+00.OR.
     *          ScalFYZR.NE.0.0D+00.OR.ScalFZYR.NE.0.0D+00) Then
              RASSPE = RAS.and.iZ.NE.iY.and.iT.GT.iU
              Call GETSGM2(ITLEV,IULEV,LSYM,CI,Work(LBUFT))
              If (.not.RASSPE.or.iZ.GE.iY) Then
                Call DaXpY_(nsgm1,ScalGYZR,WORK(LBUFT),1,Work(LBUF3),1)
                Do icsf = 1, nsgm1
                  Work(LBUF3+icsf-1) = Work(LBUF3+icsf-1)
     *              + ScalFYZR*Work(LBUFT+icsf-1)*Work(LBUFD+icsf-1)
                End Do
              End If
              If (.not.RASSPE.or.iY.GE.iZ) Then
                Call DaXpY_(nsgm1,ScalGZYR,WORK(LBUFT),1,Work(LBUF5),1)
                Do icsf = 1, nsgm1
                  Work(LBUF5+icsf-1) = Work(LBUF5+icsf-1)
     *              + ScalFZYR*Work(LBUFT+icsf-1)*Work(LBUFD+icsf-1)
                End Do
              End If
              If (RASSPE) Then
                If (iZ.LT.iY) Then
                  Call GETSGM2(IZLEV,IYLEV,LSYM,CI,Work(LBUFT))
                  Do icsf = 1, nsgm1
                    Work(LBUFT+icsf-1) =
     *                + ScalGYZR*Work(LBUFT+icsf-1)
     *                + ScalFYZR*Work(LBUFT+icsf-1)*Work(LBUFD+icsf-1)
                  End Do
                Else If (iY.LT.iZ) Then
                  Call GETSGM2(IYLEV,IZLEV,LSYM,CI,Work(LBUFT))
                  Do icsf = 1, nsgm1
                    Work(LBUFT+icsf-1) =
     *                + ScalGZYR*Work(LBUFT+icsf-1)
     *                + ScalFZYR*Work(LBUFT+icsf-1)*Work(LBUFD+icsf-1)
                  End Do
                End If
              CALL SIGMA1_CP2(ITLEV,IULEV,1.0D00,LSYM,WORK(LBUFT),CLAG,
     &             IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &             IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &             WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
              End If
            End If
          End If
C
          !! For DEPSA
          IBUF = LBUFA + MXCI*(ITLEV-1+NLEV*(IULEV-1))
          Call DaXpY_(nsgm1,ScalFYZL,Work(LTO  ),1,Work(IBUF),1)
          Call DaXpY_(nsgm1,ScalFZYL,Work(LBUF4),1,Work(IBUF),1)
          IBUF = LBUFA + MXCI*(IULEV-1+NLEV*(ITLEV-1))
          Call DaXpY_(nsgm1,ScalFYZR,Work(LTO  ),1,Work(IBUF),1)
          Call DaXpY_(nsgm1,ScalFZYR,Work(LBUF4),1,Work(IBUF),1)
C
C         G2(it,iu,iy,iz)=DDOT_(nsgm1,work(lto),1,
C    &         work(lbuf1+mxci*(ib-1)),1)
C         IF(IFF.ne.0) THEN
C           F2sum=0.0D0
C           do i=1,nsgm1
C             F2sum=F2sum+work(lto-1+i)*work(lbufd-1+i)*
C    &             work(lbuf1-1+i+mxci*(ib-1))
C           end do
C           F2(it,iu,iy,iz)=F2sum
C         END IF
        end do
        CALL SIGMA1_CP2(IZLEV,IYLEV,1.0D00,LSYM,WORK(LBUF3),CLAG,
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        If (IYLEV.NE.IZLEV) Then
          CALL SIGMA1_CP2(IYLEV,IZLEV,1.0D00,LSYM,WORK(LBUF5),CLAG,
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        End If
      end if
      nbtot=0
C     CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
C     CPUT =CPTF10-CPTF0
C     WALLT=TIOTF10-TIOTF0
C     write(6,*) "FG2     : CPU/WALL TIME=", cput,wallt
C     CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
C
C
C
      !! Work(LBUF4) is the temporary storage for <I|ExvEut|0>*Dtuvxyz
      call dcopy_(nsgm1,[0.0D0],0,work(lbuf4),1)
      do ip2=ip3,ntri2
        ivlev=idx2ij(1,ip2)
        ixlev=idx2ij(2,ip2)
        isvx=mul(ism(ivlev),ism(ixlev))
        iv=L2ACT(ivlev)
        ix=L2ACT(ixlev)
        if(isvx.ne.mul(issg1,issg2)) goto 99
        lfrom=lbuf2
        lto=lbuft
        call dcopy_(nsgm1,[0.0D0],0,work(lto),1)
        !! <0|E_{vx}E_{yz}|I>
        CALL SIGMA1_CP2(IVLEV,IXLEV,1.0D00,ISSG2,WORK(LFROM),WORK(LTO),
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
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

*-----------
* Contract the Sgm1 wave functions with the Tau wave function.
        l1=lbuf1+mxci*(ibmn-1)
C       call DGEMV_ ('T',nsgm1,nb,1.0D0,work(l1),mxci,
C    &       work(lbuft),1,0.0D0,bufr,1)
* and distribute this result into G3:
C       call dcopy_(nb,bufr,1,G3(iG3OFF+1),1)
* and copy the active indices into idxG3:
        !! Work(LBUF3) is the temporary storage for <I|Eut|0>*Dtuvxyz
        Call DCopy_(nsgm1,[0.0D0],0,Work(LBUF3),1)
        do ib=1,nb
          iG3=iG3OFF+ib
          idx=ip1_buf(ibmn-1+ib)
          itlev=idx2ij(1,idx)
          iulev=idx2ij(2,idx)
C         iT=l2act(itlev)
C         iU=l2act(iulev)
C         idxG3(1,iG3)=int(iT,I1)
C         idxG3(2,iG3)=int(iU,I1)
C         idxG3(3,iG3)=int(iV,I1)
C         idxG3(4,iG3)=int(iX,I1)
C         idxG3(5,iG3)=int(iY,I1)
C         idxG3(6,iG3)=int(iZ,I1)
C
          !! left derivative
          !! <I|EtuEvxEyz|0>
          Do icsf=1,nsgm1
            work(lbuf5-1+icsf)=
     &        DG3(iG3)*work(lbuft-1+icsf)
     &      + DF3(iG3)*(work(lbufd-1+icsf)-epsa(iv))*work(lbuft-1+icsf)
          End Do
          CALL SIGMA1_CP2(ITLEV,IULEV,1.0D+00,LSYM,WORK(LBUF5),CLAG,
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          !! right derivative (1): just construct <I|Eut|0>*Dtuvxyz
          !! <0|EtuEvxEyz|I> -> <I|EzyExvEut|0>
          If (idx.ge.ip1sta.and.idx.le.ip1end) Then
            ibuf = lbuf1+mxci*(idx-ip1sta)
            Call DCopy_(nsgm1,Work(ibuf),1,Work(LBUF5),1)
          Else
            Call DCopy_(nsgm1,[0.0D0],0,Work(LBUF5),1)
            CALL SIGMA1_CP2(IULEV,ITLEV,1.0D+00,LSYM,CI,WORK(LBUF5),
     &           IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &           IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &           WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          End If
          Do icsf=1,nsgm1
            Work(lbuf3-1+icsf)=work(lbuf3-1+icsf)
     &       + DG3(iG3)*work(LBUF5-1+icsf)
     &       + DF3(iG3)*(work(lbufd-1+icsf)-epsa(iv))*work(lbuf5-1+icsf)
          End Do
          !! DEPSA for the Work(LBUFD) term
          ibuf = lbufa + mxci*(itlev-1+nlev*(iulev-1))
          Call DaXpY_(nsgm1,DF3(iG3),Work(LBUFT),1,Work(ibuf),1)
        end do
        !! right derivative (2): construct <I|ExvEut|0>*Dtuvxyz
       CALL SIGMA1_CP2(IXLEV,IVLEV,1.0D+00,LSYM,WORK(LBUF3),WORK(LBUF4),
     &      IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &      IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &      WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C
        !! DEPSA for the -epsa(iv) term
        do ialev = 1, nlev
            Call DCopy_(nsgm1,[0.0D0],0,Work(LBUFE),1)
       CALL SIGMA1_CP2(IALEV,IXLEV,1.0D+00,LSYM,Work(LFROM),Work(LBUFE),
     &      IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &      IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &      WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          do ib = 1, nb
            iG3=iG3OFF+ib
            idx=ip1_buf(ibmn-1+ib)
            itlev=idx2ij(1,idx)
            iulev=idx2ij(2,idx)
            if (idx.ge.ip1sta.and.idx.le.ip1end) then
              ibuf = lbuf1+mxci*(idx-ip1sta)
            else
              ibuf = lbuf5
              Call DCopy_(nsgm1,[0.0D0],0,Work(LBUF5),1)
              CALL SIGMA1_CP2(IULEV,ITLEV,1.0D+00,LSYM,CI,WORK(LBUF5),
     &             IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &             IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &             WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
            end if
            DEPSA(IALEV,IVLEV) = DEPSA(IALEV,IVLEV)
     *        - DDot_(nsgm1,Work(ibuf),1,Work(LBUFE),1)*DF3(iG3)
          end do
        end do
C       IF(IFF.ne.0) THEN
* Elementwise multiplication of Tau with H0 diagonal - EPSA(IV):
C         do icsf=1,nsgm1
C           work(lbuft-1+icsf)=
C    &           (work(lbufd-1+icsf)-epsa(iv))*work(lbuft-1+icsf)
C         end do
* so Tau is now = Sum(eps(w)*E_vxww) Psi. Contract and distribute:
C         call DGEMV_ ('T',nsgm1,nb,1.0D0,work(l1),mxci,
C    &         work(lbuft),1,0.0D0,bufr,1)
C         call dcopy_(nb,bufr,1,F3(iG3OFF+1),1)
C       END IF
        iG3OFF=iG3OFF+nb
        nbtot=nbtot+nb
 99     continue
      end do
      !! right derivative (3): finally, <I|EzyExvEut|0>*Dtuvxyz
      CALL SIGMA1_CP2(IZLEV,IYLEV,1.0D+00,LSYM,WORK(LBUF4),CLAG,
     &     IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &     IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &     WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C     CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
C     CPUT =CPTF10-CPTF0
C     WALLT=TIOTF10-TIOTF0
C     write(6,*) "FG3     : CPU/WALL TIME=", cput,wallt
C
      IF(iPrGlb.GE.DEBUG) THEN
        WRITE(6,'("DEBUG> ",I8,1X,"[",I4,"..",I4,"]",1X,I4,1X,I9)')
     &    iSubTask, ip1sta, ip1end, ip3, nbtot
        call xFlush(6)
      END IF

      !! Finally, right derivative contributions of the Fock-weighted
      !! part
      !! It is the right place?
C     CALL H0DIAG_CASPT2_der(ISSG1,WORK(LFCDer),IWORK(LNOW),IWORK(LIOW),
C    *                       CLag)

CSVC: The master node now continues to only handle task scheduling,
C     needed to achieve better load balancing. So it exits from the task
C     list.  It has to do it here since each process gets at least one
C     task.
#if defined (_MOLCAS_MPP_) && !defined (_GA_)
      IF (IS_REAL_PAR().AND.KING().AND.(NPROCS.GT.1)) GOTO 501
#endif

C-SVC20100301: end of the task
      GOTO 500

 501  CONTINUE

      !! DEPSA
C     CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
      Do ITLEV = 1, NLEV
        DO IULEV = 1, NLEV
          Call DCopy_(nsgm1,[0.0D0],0,Work(LBUF5),1)
          CALL SIGMA1_CP2(IULEV,ITLEV,1.0D+00,LSYM,CI,WORK(LBUF5),
     &      IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &      IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &      WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          IBUF = LBUFA + MXCI*(ITLEV-1+NLEV*(IULEV-1))
          Do IALEV = 1, NLEV
            Do IBLEV = 1, NLEV
              Call DCopy_(nsgm1,[0.0D0],0,Work(LBUFE),1)
        CALL SIGMA1_CP2(IALEV,IBLEV,1.0D+00,LSYM,Work(IBUF),Work(LBUFE),
     &          IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &          IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &          WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
              DEPSA(IALEV,IBLEV) = DEPSA(IALEV,IBLEV)
     *          + DDot_(nsgm1,Work(LBUF5),1,Work(LBUFE),1)
            End Do
          End Do
          Call DCopy_(nsgm1,Work(LBUF5),1,Work(IBUF),1)
        End Do
      End Do
C     CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
C     CPUT =CPTF10-CPTF0
C     WALLT=TIOTF10-TIOTF0
C     write(6,*) "DEPSA(2): CPU/WALL TIME=", cput,wallt


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
C-SVC20100831: set correct number of elements in new G3
      NG3=iG3OFF
C
C     ------------------------------------------------------------------
C     --------------------DEBUG PURPOSE BEGIN---------------------------
C     ------------------------------------------------------------------
C     may be symmetrizations are mandatory here?
C      call dcopy_(nlev**4,df2bk,1,df2,1)
C
C     write(6,*) "dg2 after"
C     do it = 1, nlev
C     do iu = 1, nlev
C     do iy = 1, nlev
C     do iz = 1, nlev
C     write(6,'(4i3,f20.10)') it,iu,iy,iz,dg2(it,iu,iy,iz)
C     end do
C     end do
C     end do
C     end do
C       do it = 1, nlev
C       do iu = 1, nlev
C         call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
C         CALL SIGMA1_CP2(IT,IU,1.0d+00,LSYM,CI,Work(lbuf1),
C    &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C      scal=ddot_(nsgm1,work(lbuf1),1,clag,1)
C      write(6,'(2i3,f20.10)') it,iu,scal
C       end do
C       end do
C     write(6,*) "clag=",clag(1)
1234  continue
      if (.false..or.dodebug) then
        !! RDM1
        do it = 1, nlev
        do iu = 1, nlev
          scal = dg1(it,iu)+dg1(iu,it)
          CALL SIGMA1_CP2(IT,IU,scal,LSYM,CI,CLAG,
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        end do
        end do
      write(6,*) "clag=",clag(1)

        !! DF1
        !! right derivative
        do it = 1, nlev
        do iu = 1, nlev
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
          CALL SIGMA1_CP2(IT,IU,1.0d+00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          do iv = 1, nlev
            scal = df1(it,iu)*epsa(iv)
            CALL SIGMA1_CP2(IV,IV,SCAL,LSYM,Work(LBUF1),CLAG,
     &           IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &           IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &           WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          end do
        end do
        end do

        !! left derivative
        do iv = 1, nlev
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
          CALL SIGMA1_CP2(IV,IV,1.0d+00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          do it = 1, nlev
          do iu = 1, nlev
            scal = df1(it,iu)*epsa(iv)
            CALL SIGMA1_CP2(IU,IT,SCAL,LSYM,Work(LBUF1),CLAG,
     &           IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &           IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &           WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          end do
          end do
        end do
      write(6,*) "clag=",clag(1)

        !! RDM2

        !! right derivative
        !! <0|E_{tu}E_{xy}|I> * D_{tuxy}
        !! = CI(J)*<J|E_{tu}|K>*<K|E_{xy}|I> * D_{tuxy}
        do it = 1, nlev
        do iu = 1, nlev
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
          CALL SIGMA1_CP2(IT,IU,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          do ix = 1, nlev
          do iy = 1, nlev
            scal = dg2(it,iu,ix,iy)
            CALL SIGMA1_CP2(IX,IY,SCAL,LSYM,Work(LBUF1),CLAG,
     &           IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &           IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &           WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          end do
          end do
        end do
        end do
        !! left derivative
        !! <I|E_{tu}E_{xy}|0> * D_{tuxy}
        !! = <I|E_{tu}|K>*<K|E_{xy}|J>*CI(J) * D_{tuxy}
        !! = <I|E_{tu}|K>*<J|E_{yx}|K>*CI(J) * D_{tuxy}
        !! = <K|E_{ut}|I>*<0|E_{yx}|K> * D_{tuxy}
        do ix = 1, nlev
        do iy = 1, nlev
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
          CALL SIGMA1_CP2(IY,IX,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          do it = 1, nlev
          do iu = 1, nlev
            scal = dg2(it,iu,ix,iy)
            CALL SIGMA1_CP2(IU,IT,SCAL,LSYM,Work(LBUF1),CLAG,
     &           IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &           IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &           WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          end do
          end do
        end do
        end do
      write(6,*) "clag=",clag(1)

        !! DF2

        !! EPSA derivative and right derivative
        do it = 1, nlev
        do iu = 1, nlev
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
          CALL SIGMA1_CP2(IT,IU,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          do iv = 1, nlev
          do ix = 1, nlev
            call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
            CALL SIGMA1_CP2(IV,IX,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
     &           IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &           IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &           WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
            do iy = 1, nlev
            do iz = 1, nlev
              call dcopy_(nsgm1,[0.0d+00],0,work(lbuf3),1)
            !! for EPSA derivative
           CALL SIGMA1_CP2(IY,IZ,1.0D+00,LSYM,Work(LBUF2),Work(Lbuf3),
     &          IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &          IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &          WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
              scal = ddot_(nsgm1,ci,1,Work(lbuf3),1)
              depsa(iv,ix) = depsa(iv,ix) + scal*df2(it,iu,iy,iz)
            !! for right derivative
            if (iv.eq.ix) then
            scal = df2(it,iu,iy,iz)*epsa(iv)
            CALL SIGMA1_CP2(IY,IZ,SCAL,LSYM,Work(LBUF2),CLAG,
     &           IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &           IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &           WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
            end if
            end do
            end do
          end do
          end do
        end do
        end do
        !! left derivative
        !!   <I|E_{tu}E_{vx}E_{yz}|0>*f_{vx}
        !! = <I|E_{tu}|J>*<J|E_{vx}|K>*<K|E_{yz}|L>*CI(L)*f_{vx} * DF2(tuyz)
        !! = <J|E_{ut}|I>*<K|E_{xv}|J>*<L|E_{zy}|K>*CI(L)*f_{vx} * DF2(tuyz)
        !! = CI(L)*<L|E_{zy}|K>*<K|E_{xv}|J>*<J|E_{ut}|I>*f_{vx} * DF2(tuyz)
        do iy = 1, nlev
        do iz = 1, nlev
          call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
          CALL SIGMA1_CP2(IZ,IY,1.0D00,LSYM,CI,Work(LBUF1),
     &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          do iv = 1, nlev
            call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
            CALL SIGMA1_CP2(IV,IV,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
     &           IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &           IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &           WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
            do it = 1, nlev
            do iu = 1, nlev
              scal = df2(it,iu,iy,iz)*epsa(iv)
           CALL SIGMA1_CP2(IU,IT,scal,LSYM,Work(LBUF2),CLAG,
     &          IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &          IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &          WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
            end do
            end do
          end do
        end do
        end do
      write(6,*) "clag=",clag(1)
C
      !! DG3
      !! <0|EtuEvxEyz|0>*Dtu,vx,yz
      !! -> <I|EtuEvxEyz|0>*Dtu,vx,yz
      !!  + <0|EtuEvxEyz|I>*Dtu,vx,yz

      ! CJ *<J| t+ u v+ x y+ z |I>
      !=CJ *<I|z+ y x+ v u+ t|J>
      tmp = 0.0d+00
      do i = 1, nlev
        tmp = tmp + epsa(i)
      end do
      write(6,*) "ng3 = ", ng3
      Do iG3 = 1, NG3
C     write(6,'(i4,2f20.10)') ig3,df3(ig3),dg3(ig3)
        iT=idxG3(1,iG3)
        iU=idxG3(2,iG3)
        iV=idxG3(3,iG3)
        iX=idxG3(4,iG3)
        iY=idxG3(5,iG3)
        iZ=idxG3(6,iG3)
C       write(6,'(6i3)') it,iu,iv,ix,iy,iz
        if (it.eq.0) exit

C       iT=idxG3(1,iG3)
C       iV=idxG3(2,iG3)
C       iY=idxG3(3,iG3)
C       iU=idxG3(4,iG3)
C       iX=idxG3(5,iG3)
C       iZ=idxG3(6,iG3)
        call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
        CALL SIGMA1_CP2(IY,IZ,1.0D+00,LSYM,CI,Work(LBUF1),
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
        CALL SIGMA1_CP2(IV,IX,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        scal = dg3(ig3)
        scal = scal - df3(ig3)*epsa(iv)
        CALL SIGMA1_CP2(IT,IU,SCAL,LSYM,Work(LBUF2),CLAG,
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))

        call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
        CALL SIGMA1_CP2(IU,IT,1.0D+00,LSYM,CI,Work(LBUF1),
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
        CALL SIGMA1_CP2(IX,IV,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        scal = dg3(ig3)
        scal = scal - df3(ig3)*epsa(iv)
        CALL SIGMA1_CP2(IZ,IY,SCAL,LSYM,Work(LBUF2),CLAG,
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))

        !! depsa
        do iw = 1, nlev
        call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
        CALL SIGMA1_CP2(IY,IZ,1.0D+00,LSYM,CI,Work(LBUF1),
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
        CALL SIGMA1_CP2(Iw,IX,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
        CALL SIGMA1_CP2(IT,IU,1.0D+00,LSYM,work(lbuf2),Work(LBUF1),
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
         depsa(iv,iw) = depsa(iv,iw)
     *     - df3(ig3)*ddot_(nsgm1,ci,1,work(lbuf1),1)
        end do



      !! EtuEvxEyz,ab*fab ?
      !! = t+ u v+ x y+ a+ b z * fab
      !! = del(ab) t+ u v+ x y+ z - t+ u v+ x y+ b a+ z

      do iw = 1, nlev
        call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
        CALL SIGMA1_CP2(IY,IZ,1.0D+00,LSYM,CI,Work(LBUF1),
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
        CALL SIGMA1_CP2(IV,IX,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
        CALL SIGMA1_CP2(IW,IW,1.0D+00,LSYM,work(lbuf2),Work(LBUF1),
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        scal = df3(ig3)*epsa(iw)
        CALL SIGMA1_CP2(IT,IU,SCAL,LSYM,Work(LBUF1),CLAG,
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))

        call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
        CALL SIGMA1_CP2(IU,IT,1.0D+00,LSYM,CI,Work(LBUF1),
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
        CALL SIGMA1_CP2(IW,IW,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
        CALL SIGMA1_CP2(IX,IV,1.0D+00,LSYM,Work(LBUF2),Work(LBUF1),
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        scal = df3(ig3)*epsa(iw)
        CALL SIGMA1_CP2(IZ,IY,SCAL,LSYM,Work(LBUF1),CLAG,
     &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
      end do
        !! depsa
      ! do ia = 1, nlev
      ! do ib = 1, nlev
      !   call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
      !   CALL SIGMA1_CP2(IY,IZ,1.0D+00,LSYM,CI,Work(LBUF1),
     &!        IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &!        IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &!        WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
      !   call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
      !   CALL SIGMA1_CP2(IV,IX,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
     &!        IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &!        IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &!        WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
      !   call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
      !   CALL SIGMA1_CP2(IA,IB,1.0D+00,LSYM,work(lbuf2),Work(LBUF1),
     &!        IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &!        IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &!        WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
      !   call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
      !   CALL SIGMA1_CP2(IT,IU,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
     &!        IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &!        IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &!        WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
      !  depsa(ia,ib) = depsa(ia,ib)
     *!    + df3(ig3)*ddot_(nsgm1,ci,1,work(lbuf2),1)
      ! end do
      ! end do
      End Do
      write(6,*) "asdf"
C       do it = 1, nlev
C       do iu = 1, nlev
C         call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
C         CALL SIGMA1_CP2(IT,IU,1.0D00,LSYM,CI,Work(LBUF1),
C    &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C         call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
C         do iv = 1, nlev
C           CALL SIGMA1_CP2(IV,IV,epsa(iv),LSYM,Work(LBUF1),Work(LBUF2),
C    &           IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &           IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &           WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C         end do
C         write(6,*) it,iu,ddot_(nsgm1,work(lbuf2),1,ci,1)
C       end do
C       end do
C
        !! if all indices are different, only the following symmetry holds true:
        !! <0| E_tu H0Diag E_yz |0> = <0| E_zy H0Diag E_ut |0>
        !! For F2: F2(tuyz) = F2(utzy) = F2(yztu) = F2(zyut)
        !! Also, E_{tu}E_{yz} = E_{ut}E_{zy} = E_{yz}E_{tu} = E_{zy}E_{ut}
      ! do iy = 1, nlev
      ! do iz = 1, nlev
      !   call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
      !   call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
      !   CALL SIGMA1_CP2(IY,IZ,SCAL,LSYM,CI,Work(LBUF1),
     &!        IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &!        IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &!        WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
      !   CALL SIGMA1_CP2(IZ,IY,SCAL,LSYM,CI,Work(LBUF2),
     &!        IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &!        IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &!        WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
      !   do i = 1, nsgm1
      !     work(lbuf1+i-1) = work(lbuf1+i-1)*work(lbufd+i-1)
      !     work(lbuf2+i-1) = work(lbuf2+i-1)*work(lbufd+i-1)
      !   end do
      !   do it = 1, nlev
      !   do iu = 1, nlev
      !     !! for right derivative
      ! !   scal = df2(it,iu,iy,iz)
      ! !   write(6,'(4i3,f20.10)') it,iu,iy,iz,scal
      ! !   CALL SIGMA1_CP2(IT,IU,SCAL,LSYM,Work(LBUF2),CLAG,
     &! !        IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &! !        IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &! !        WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C     ! !   scal = df2(iz,iy,iu,it)
      ! !   CALL SIGMA1_CP2(IU,IT,SCAL,LSYM,Work(LBUF1),CLAG,
     &! !        IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &! !        IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &! !        WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))

C     !   call dcopy_(nsgm1,[0.0d+00],0,work(lfcder2),1)
C     !   CALL SIGMA1_CP2(IY,IZ,1.0D+00,LSYM,Work(LBUF1),work(lfcder2),
C    &!        IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &!        IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &!        WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C     !   tmp = ddot_(nsgm1,work(lfcder2),1,ci,1)
C     !   call dcopy_(nsgm1,[0.0d+00],0,work(lfcder2),1)
C     !   CALL SIGMA1_CP2(IY,IZ,1.0D+00,LSYM,Work(LBUF2),work(lfcder2),
C    &!        IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &!        IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &!        WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C     !   tmp2= ddot_(nsgm1,work(lfcder2),1,ci,1)
C     !     write(6,'(4i3,2f20.10)') it,iu,iy,iz,tmp,tmp2
      !   end do
      !   end do
      ! end do
      ! end do




        !! DF3
C     Do iG3 = 1, NG3
C       iT=idxG3(1,iG3)
C       iU=idxG3(2,iG3)
C       iV=idxG3(3,iG3)
C       iX=idxG3(4,iG3)
C       iY=idxG3(5,iG3)
C       iZ=idxG3(6,iG3)
C       call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
C       CALL SIGMA1_CP2(IY,IZ,1.0D+00,LSYM,CI,Work(LBUF1),
C    &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C       call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
C       CALL SIGMA1_CP2(IV,IX,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
C    &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C       scal = dg3(ig3)
C       CALL SIGMA1_CP2(IT,IU,SCAL,LSYM,Work(LBUF2),CLAG,
C    &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))

C       call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
C       CALL SIGMA1_CP2(IU,IT,1.0D+00,LSYM,CI,Work(LBUF1),
C    &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C       call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
C       CALL SIGMA1_CP2(IX,IV,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
C    &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C       scal = dg3(ig3)
C       CALL SIGMA1_CP2(IZ,IY,SCAL,LSYM,Work(LBUF2),CLAG,
C    &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C     End Do


        !! new in 2020/08/30
        !! left derivative
      ! do i = 1, nsgm1
      !   work(lbuf1+i-1) = work(lbufd+i-1)*ci(i)
      ! end do
      ! do iy = 1, nlev
      ! do ix = 1, nlev
      !   call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
      !   CALL SIGMA1_CP2(IX,IY,1.0d+00,LSYM,work(lbuf1),Work(LBUF2),
     &!        IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &!        IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &!        WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
      !   do it = 1, nlev
      !   do iu = 1, nlev
      !     scal = df2(it,iu,ix,iy)
C     !     CALL SIGMA1_CP2(IT,IU,SCAL,LSYM,work(lbuf2),CLAG,
C    &!          IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &!          IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &!          WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
      !   end do
      !   end do
      ! end do
      ! end do
      ! !! right derivative
      ! do it = 1, nlev
      ! do iu = 1, nlev
      !   call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
      !   CALL SIGMA1_CP2(IU,IT,1.0d+00,LSYM,CI,Work(LBUF2),
     &!        IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &!        IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &!        WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
      !   do ix = 1, nlev
      !   do iy = 1, nlev
      !     call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
      !     CALL SIGMA1_CP2(IY,IX,1.0d+00,LSYM,work(lbuf1),work(lbuf2),
     &!          IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &!          IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &!          WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
      !     do iv = 1, nlev
      !       scal = df2(it,iu,ix,iy)*epsa(iv)
C     !       CALL SIGMA1_CP2(IV,IV,SCAL,LSYM,work(lbuf2),CLAG,
C    &!            IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &!            IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &!            WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
      !     end do
      !   end do
      !   end do
      ! end do
      ! end do


      !! test for f3
      !! Etu*(Eww*fww - fvv)*Evx*Eyz
      !! = t+ u ( w+ w fww - fvv) v+ x y+ z
      !! = t+ u ( w+ w fww - del(wv) fwv) v+ x y+ z
      !! = t+ u w+ v+ w x y+ z * f_{ww}
      !! = EtuEwv,xwEyz
C     call dcopy_(ng3,[0.0d+00],0,f3,1)
C     tmp = 0.0d+00
C     do i = 1, nlev
C       tmp = tmp + epsa(i)
C     end do
C     do ig3 = 1, ng3
C       iT=idxG3(1,iG3)
C       iU=idxG3(2,iG3)
C       iV=idxG3(3,iG3)
C       iX=idxG3(4,iG3)
C       iY=idxG3(5,iG3)
C       iZ=idxG3(6,iG3)

C       !! del(ww) t+ u v+ x y+ z * f_{ww}
C       call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
C       CALL SIGMA1_CP2(IY,IZ,1.0D+00,LSYM,CI,Work(LBUF1),
C    &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C       call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
C       CALL SIGMA1_CP2(IV,IX,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
C    &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C       call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
C       CALL SIGMA1_CP2(IT,IU,1.0d+00,LSYM,Work(LBUF2),Work(LBUF1),
C    &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C       f3(ig3) = f3(ig3) - epsa(iv)*ddot_(nsgm1,work(lbuf1),1,ci,1)

C       !! - t+ u v+ x y+ w w+ z * f_{ww}
C       do iw = 1, nlev
C       call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
C       CALL SIGMA1_CP2(IY,IZ,1.0D+00,LSYM,CI,Work(LBUF1),
C    &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C       call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
C       CALL SIGMA1_CP2(IV,IX,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
C    &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C       call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
C       CALL SIGMA1_CP2(IW,IW,1.0d+00,LSYM,Work(LBUF2),Work(LBUF1),
C    &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C       call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
C       CALL SIGMA1_CP2(IT,IU,1.0d+00,LSYM,Work(LBUF1),Work(LBUF2),
C    &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C       f3(ig3) = f3(ig3) + epsa(iw)*ddot_(nsgm1,work(lbuf2),1,ci,1)
C       end do
C       !! write(6,'(i5,6i2,f20.10)') ig3,it,iu,iv,ix,iy,iz,f3(ig3)
C     end do



      !!   t+ u v+ v y+ z * f_{vv}
      !! = kkkk

      !! DF2(t,u,x,y)
      !! = t+ u+ v+ w y x * f_{vw}
      !! = t+ u+ v+ y x w * f_{vw}
      !! = t+ u+ (del(vy) - y v+) x w * f_{vw}
      !! = del(vy) t+ u+ x w * f_{vw}
      !!   - t+ u+ y v+ x w * f_{vw}
      !! = del(vy) E_{tu,wx} * f_{vw}
      !!   - t+ (del(uy) - y u+)*(del(vx) - x v+) w * f_{vw}
      !! = del(vy) E_{tu,wx} * f_{vw}
      !!   - del(uy)*del(vx)*G1(t,w) * f_{vw}
      !!   + del(uy) t+ x v+ w * f_{vw}
      !!   + del(vx) t+ y u+ w * f_{vw}
      !!   - t+ y u+ x v+ w * f_{vw}
      !! = del(vy) E_{tu,wx} * f_{vw}
      !!   - del(uy)*del(vx)*G1(t,w) * f_{vw}
      !!   + del(uy) E_{tx}*E_{vw} * f_{vw}
      !!   + del(vx) E_{ty}*E_{uw} * f_{vw}
      !!   - E_{ty}E_{ux}E_{vw}*f_{vw}
      !! = del(vy) E_{tu,wx} * f_{vw}
      !!   - del(uy)*del(vx)*G1(t,w) * f_{vw}
      !!   + del(uy) (F1(t,x) + del(xv)*G1(t,w)*f_{vw})
      !!   + del(vx) (del(yu) G1(t,w) - E(tu,wy)) * f_{vw}
      !!   - E_{ty}E_{ux}E_{vw}*f_{vw}
      !! ->  E_{tu,yx}*e(y)
      !!    -E_{tu,xy}*e(x)
      !!    + del(uy)*F1(t,x)
      !!    - del(uy)*G1(t,x)*e(x)
      !!    + del(uy)*G1(t,x)*e(x)
      !!    + del(uy)*G1(t,x)*e(x)
      !!    - E_{ty}E_{ux}E_{vw}*f_{vw}

      !! DF2(t,u,x,y)
      !! = t+ u+ v+ w y x * f_{vw}
      !! = t+ u+ v+ y x w * f_{vw}
      !! =-t+ u+ v+ x y w * f_{vw}
      !! =-t+ u+ (del(vx) - x v+) y w * f_{vw}
      !! =-del(vx) t+ u+ y w * f_{vw}
      !!  + t+ u+ x v+ y w * f_{vw}
      !! =-del(vx) E_{tu,wy}*f_{vw}
      !!  + t+ (del(ux) - x u+)*(del(vy) - y v+) w * f_{vw}
      !! =-del(vx) E_{tu,wy}*f_{vw}
      !!  + del(ux)*del(vy) G1(t,w)*f_{vw}
      !!  - del(ux) t+ y v+ w * f_{vw}
      !!  - del(vy) t+ x u+ w * f_{vw}
      !!  + t+ x u+ y v+ w * f_{vw}
      !! =-del(vx) E_{tu,wy}*f_{vw}
      !!  + del(ux)*del(vy) G1(t,w)*f_{vw}
      !!  - del(ux) E_{ty}E_{vw}*f_{vw}
      !!  - del(vy) E_{tx}E_{uw}*f_{vw}
      !!  + E_{tx}E_{uy}E_{vw}*f_{vw}
      !! =-del(vx) E_{tu,wy}*f_{vw}
      !!  + del(ux)*del(vy) G1(t,w)*f_{vw}
      !!  - del(ux)*(F1(ty)+del(vy)*G1(tw)*f_{vw})
      !!  - del(vy)*(del(xu)*G1(t,w)+E_{tu,xw})*f_{vw}
      !!  + E_{tx}E_{uy}E_{vw}*f_{vw}
      !! =- E_{tu,xy}*e_(x)
      !!  + del(ux)*G1(t,y)*e(y)
      !!  + del(ux)*G1(t,w)*e(y)
      !!  - del(ux)*F1(t,y)
      !!  - del(ux)*G1(t,y)*e(y)
      !!  - del(ux)*G1(t,y)*e(y)
      !!  - E_{tu,xy}*e(y)
      !!  + E_{tx}E_{uy}E_{vw}*f_{vw}

      !! DF2(t,u,x,y)
      !! = -E_{tu,xy}*(e(x)+e(y))
      !!   -del(ux)(F1(t,y)+G1(t,y)*eps(y))
      !!   +E_{tx}E_{uy}E_{vw}*f_{vw}

      !! DF1(t,u)
      !! = t+ v+ w u * f_{vw}
      !! = -t+ v+ u w * f_{vw}
      !! = -del(vu)*G1(t,w)*f_{vw} + E_{tu}E_{vw}*f_{vw}
      !! E_{tu}E_{xy} = t+ u x+ y
      !! = t+ (del(ux) - x+ u) y
      !! = del(ux) G1(t,y) - E_{tx,yu}
      !! = del(ux) G1(t,y) + E_{tx,uy}
C     call dcopy_(nlev**4,[0.0d+00],0,df2,1)
C     call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
C     call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
C     call dcopy_(nsgm1,work(lbufd),1,work(lbuf1),1)
C     do i = 1, nsgm1
C       work(Lbuf1+i-1) = work(lbuf1+i-1)*ci(i)
C     end do
C     do iu = 1, nlev
C     do iy = 1, nlev
C       call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
C       CALL SIGMA1_CP2(IU,IY,1.0D+00,LSYM,Work(LBUF1),Work(LBUF2),
C    &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C       do it = 1, nlev
C       do ix = 1, nlev
C         call dcopy_(nsgm1,[0.0d+00],0,work(lbufd),1)
C         CALL SIGMA1_CP2(IT,IX,1.0D+00,LSYM,Work(LBUF2),Work(lbufd),
C    &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C         scal = +ddot_(nsgm1,ci,1,work(lbufd),1)
C         scal = scal - g2(it,ix,iu,iy)*epsa(iy)
C    *                - g2(it,ix,iu,iy)*epsa(ix)
C         if (iu.eq.ix) scal = scal - f1(it,iy)
C         if (iu.eq.ix) scal = scal - g1(it,iy)*epsa(iy)
C         df2(it,ix,iu,iy) = scal
C       end do
C       end do
C     end do
C     end do

      !! g2(tu,xy)
      !! = t+ u+ y x
      !! = del(uy) G1(t,x) - E_{ty}E_{ux}
      !! = -t+ u+ x y
      !! = -del(ux)*G1(t,y) + E_{tx}E_{uy}

      !! g2(tu,xy)
      !!\= t+ u+ x y
      !! = del(ux) G1(ty) - E_{tx}E_{uy}
C     do iu = 1, nlev
C     do ix = 1, nlev
C       call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
C       CALL SIGMA1_CP2(IU,IX,1.0D+00,LSYM,CI,Work(LBUF1),
C    &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C       do it = 1, nlev
C       do iy = 1, nlev
C         call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
C         CALL SIGMA1_CP2(IT,IY,1.0D+00,LSYM,Work(LBUF1),Work(lbuf2),
C    &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C         scal = -ddot_(nsgm1,ci,1,work(lbuf2),1)
C         if (iu.eq.iy) scal = scal + g1(it,ix)
C         dg2(it,iu,ix,iy) = scal
C       end do
C       end do
C     end do
C     end do

      do it = 1, nlev
      do iu = 1, nlev
      do iv = 1, nlev
      do ix = 1, nlev
C       tmp =-df2(it,iu,iv,ix)
C       tmp2=  f2(it,iv,iu,ix)
        tmp = df2(it,iu,iv,ix)
        tmp2=  f2(it,iu,iv,ix)
C       tmp = dg2(it,iu,iv,ix)
C       tmp2=  g2(it,iv,iu,ix)
C     write(6,'(4i3,3f20.10)') it,iu,iv,ix,tmp,tmp2,abs(tmp-tmp2)
      end do
      end do
      end do
      end do
C     do it = 1, nlev
C     do iu = 1, nlev-1
C       depsa(it,iu) = (depsa(it,iu) + depsa(iu,it))*0.5d+00
C       depsa(iu,it) = (depsa(it,iu) + depsa(iu,it))*0.5d+00
C     end do
C     end do


        !! <CLAG|Etu|CI>
      ! call dcopy_(nlev*nlev,[0.0d+00],0,g1,1)
      ! do it = 1, nlev
      ! do iu = 1, nlev
      !   call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
      !   CALL SIGMA1_CP2(IT,IU,1.0D+00,LSYM,CI,Work(LBUF1),
     &!        IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &!        IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &!        WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
      !   scal = ddot_(nsgm1,clag,1,work(lbuf1),1)
      !   g1(it,iu) = g1(it,iu) + scal
      !   g1(iu,it) = g1(iu,it) + scal
      ! end do
      ! end do
      ! write(6,*) "<CLAG|Etu|CI>+<CI|Etu|CLAG>"
      ! call sqprt(g1,nlev)
C
      ! scal = ddot_(nsgm1,ci,1,clag,1)
      ! call daxpy_(nsgm1,-scal,ci,1,clag,1)
      ! !! <CLAG|Etu|CI>
      ! call dcopy_(nlev*nlev,[0.0d+00],0,g1,1)
      ! do it = 1, nlev
      ! do iu = 1, nlev
      !   call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
      !   CALL SIGMA1_CP2(IT,IU,1.0D+00,LSYM,CI,Work(LBUF1),
     &!        IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &!        IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &!        WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
      !   scal = ddot_(nsgm1,clag,1,work(lbuf1),1)
      !   g1(it,iu) = g1(it,iu) + scal
      !   g1(iu,it) = g1(iu,it) + scal
      ! end do
      ! end do
      ! write(6,*) "<CLAG|Etu|CI>+<CI|Etu|CLAG> after projection"
      ! call sqprt(g1,nlev)
      end if
C     write(6,*) "EtuEyz"
C     do iy = 1, nlev
C     do iz = 1, nlev
C       call dcopy_(nsgm1,[0.0d+00],0,work(lbuf1),1)
C       CALL SIGMA1_CP2(IY,IZ,1.0D+00,LSYM,CI,Work(LBUF1),
C    &       IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &       IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &       WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C       do it = 1, nlev
C       do iu = 1, nlev
C         call dcopy_(nsgm1,[0.0d+00],0,work(lbuf2),1)
C         CALL SIGMA1_CP2(IT,IU,1.0D+00,LSYM,Work(LBUF1),Work(lbuf2),
C    &         IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
C    &         IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
C    &         WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
C         scal = ddot_(nsgm1,ci,1,work(lbuf2),1)
C         write(6,'(4i3,f20.10)') it,iu,iy,iz,scal
C       end do
C       end do
C     end do
C     end do
C
C     ------------------------------------------------------------------
C     --------------------DEBUG PURPOSE END-----------------------------
C     ------------------------------------------------------------------
C
C     do it = 1, nlev
C     do iu = 1, nlev
C     do ix = 1, nlev
C        G2(it,ix,ix,iu)=G2(it,ix,ix,iu)+G1(it,iu)
C     end do
C     end do
C     end do
C     scal = ddot_(625,dg2,1,g2,1)
C     write(6,*) "scal = ", scal
C     scal = ddot_(nsgm1,clag,1,ci,1)
C     write(6,*) "scal = ", scal
C        write(6,*) "dg1"
C        call sqprt(dg1,5)
C        call dcopy_(25,[0.0d+00],0,dg1,1)
C       write(6,*) "after"
C       do i = 1, nsgm1
C         write(6,'(i3,f20.10)') i,work(lbufd+i-1)
C       end do

      CALL GETMEM ('TASKLIST','FREE','INTE',lTask_List,4*mxTask)
      ! free CI buffers
      CALL GETMEM('BUF1','FREE','REAL',LBUF1,NBUF1*MXCI)
      CALL GETMEM('BUF2','FREE','REAL',LBUF2,NBUF2*MXCI)
      CALL GETMEM('BUF3','FREE','REAL',LBUF3,NBUF3*MXCI)
      CALL GETMEM('BUF4','FREE','REAL',LBUF4,NBUF4*MXCI)
      CALL GETMEM('BUF5','FREE','REAL',LBUF5,NBUF5*MXCI)
      CALL GETMEM('BUF6','FREE','REAL',LBUF6,NBUF6*MXCI)
      CALL GETMEM('BUFT','FREE','REAL',LBUFT,NBUFT*MXCI)
      CALL GETMEM('BUFD','FREE','REAL',LBUFD,NBUFD*MXCI)
      CALL GETMEM('BUFE','FREE','REAL',LBUFE,NBUFE*MXCI)
      CALL GETMEM('BUFA','FREE','REAL',LBUFA,NBUFA*MXCI)
C
C     CALL GETMEM('FCDER1','FREE','REAL',LFCDer1,NBUFD*MXCI)
C     CALL GETMEM('FCDER2','FREE','REAL',LFCDer2,NLEV*NLEV)

C-SVC20100302: Synchronized add into the densitry matrices
C  only for the G1 and G2 replicate arrays
C     CALL GADSUM(G1,NG1)
C     CALL GADSUM(G2,NG2)

C     CALL GADSUM(F1,NG1)
C     CALL GADSUM(F2,NG2)

* Correction to G2: It is now = <0| E_tu E_yz |0>
C     do iu=1,nlev
C      do iz=1,nlev
C       do it=1,nlev
C        G2(it,iu,iu,iz)=G2(it,iu,iu,iz)-G1(it,iz)
C       end do
C      end do
C     end do
C-SVC20100310: took some spurious mirroring of G2 values out
C-of the loops and put them here, after the parallel section has
C-finished, so that GAdSUM works correctly.
C     do ip1=ntri2+1,nlev2
C      itlev=idx2ij(1,ip1)
C      iulev=idx2ij(2,ip1)
C      it=L2ACT(itlev)
C      iu=L2ACT(iulev)
C      do ip3=ntri1+1,ip1
C       iylev=idx2ij(1,ip3)
C       izlev=idx2ij(2,ip3)
C       iy=L2ACT(iylev)
C       iz=L2ACT(izlev)
C       G2(it,iu,iy,iz)=G2(iz,iy,iu,it)
C      end do
C     end do
* Correction to G2: Some values not computed follow from symmetry
C     do ip1=1,nlev2-1
C      itlev=idx2ij(1,ip1)
C      iulev=idx2ij(2,ip1)
C      it=L2ACT(itlev)
C      iu=L2ACT(iulev)
C      do ip3=ip1+1,nlev2
C       iylev=idx2ij(1,ip3)
C       izlev=idx2ij(2,ip3)
C       iy=L2ACT(iylev)
C       iz=L2ACT(izlev)
C       G2(it,iu,iy,iz)=G2(iy,iz,it,iu)
C      end do
C     end do
C     IF(IFF.ne.0) THEN
* Correction to F2: It is now = <0| E_tu H0Diag E_yz |0>
C      do iz=1,nlev
C       do iy=1,nlev
C        do iu=1,nlev
C         do it=1,nlev
C          F2(it,iu,iy,iz)=F2(it,iu,iy,iz)-
C    &           (EPSA(iu)+EPSA(iy))*G2(it,iu,iy,iz)
C         end do
C        end do
C       end do
C      end do
C      do iz=1,nlev
C       do iu=1,nlev
C        do it=1,nlev
C         F2(it,iu,iu,iz)=F2(it,iu,iu,iz)-
C    &          (F1(it,iz)+EPSA(iu)*G1(it,iz))
C        end do
C       end do
C      end do
C-SVC20100310: took some spurious mirroring of F2 values out
C-of the loops and put them here, after the parallel section has
C-finished, so that GAdSUM works correctly.
C      do ip1=ntri2+1,nlev2
C       itlev=idx2ij(1,ip1)
C       iulev=idx2ij(2,ip1)
C       it=L2ACT(itlev)
C       iu=L2ACT(iulev)
C       do ip3=ntri1+1,ip1
C        iylev=idx2ij(1,ip3)
C        izlev=idx2ij(2,ip3)
C        iy=L2ACT(iylev)
C        iz=L2ACT(izlev)
C        F2(it,iu,iy,iz)=F2(iz,iy,iu,it)
C       end do
C      end do
* Correction to F2: Some values not computed follow from symmetry
C      do ip1=1,nlev2-1
C       itlev=idx2ij(1,ip1)
C       iulev=idx2ij(2,ip1)
C       it=L2ACT(itlev)
C       iu=L2ACT(iulev)
C       do ip3=ip1+1,nlev2
C        iylev=idx2ij(1,ip3)
C        izlev=idx2ij(2,ip3)
C        iy=L2ACT(iylev)
C        iz=L2ACT(izlev)
C        F2(it,iu,iy,iz)=F2(iy,iz,it,iu)
C       end do
C      end do
C     END IF

* Correction to G3: It is now <0| E_tu E_vx E_yz |0>
* Similar for F3 values.
C     DO iG3=1,NG3
C      iT=idxG3(1,iG3)
C      iU=idxG3(2,iG3)
C      iV=idxG3(3,iG3)
C      iX=idxG3(4,iG3)
C      iY=idxG3(5,iG3)
C      iZ=idxG3(6,iG3)
* Correction: From <0| E_tu E_vx E_yz |0>, form <0| E_tuvxyz |0>
C      if(iY.eq.iX) then
C       G3(iG3)=G3(iG3)-G2(iT,iU,iV,iZ)
C       IF(IFF.ne.0) F3(iG3)=
C    &           F3(iG3)-(F2(iT,iU,iV,iZ)+EPSA(iu)*G2(iT,iU,iV,iZ))
C       if(iv.eq.iu) then
C        G3(iG3)=G3(iG3)-G1(iT,iZ)
C        IF(IFF.ne.0) F3(iG3)=F3(iG3)-F1(iT,iZ)
C       end if
C      end if
C      if(iV.eq.iU) then
C        G3(iG3)=G3(iG3)-G2(iT,iX,iY,iZ)
C        IF(IFF.ne.0) F3(iG3)=F3(iG3)-
C    &            (F2(iT,iX,iY,iZ)+EPSA(iY)*G2(iT,iX,iY,iZ))
C      end if
C      if(iY.eq.iU) then
C        G3(iG3)=G3(iG3)-G2(iV,iX,iT,iZ)
C        IF(IFF.ne.0) F3(iG3)=
C    &            F3(iG3)-(F2(iV,iX,iT,iZ)+EPSA(iU)*G2(iV,iX,iT,iZ))
C      end if
C      IF(IFF.ne.0) F3(iG3)=F3(iG3)-(EPSA(iU)+EPSA(iY))*G3(iG3)
C     END DO

 999  continue
      CALL QEXIT('DERFG3')
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      Subroutine DerFG3_Trf(NLEV,TRF,DG1,DF1)
C
      Implicit Real*8 (A-H,O-Z)
C
#include "WrkSpc.fh"
C
      Dimension TRF(NLEV,NLEV),DG1(NLEV,NLEV),DF1(NLEV,NLEV)
C
      CALL GETMEM('WRK1','ALLO','REAL',LWRK1,NLEV*NLEV)
C
      !! DG1
      Call DGemm_('t','N',NLEV,NLEV,NLEV,
     *            1.0D+00,TRF,NLEV,DG1,NLEV,
     *            0.0D+00,Work(LWRK1),NLEV)
      Call DGemm_('N','n',NLEV,NLEV,NLEV,
     *            1.0D+00,Work(LWRK1),NLEV,TRF,NLEV,
     *            0.0D+00,DG1,NLEV)
C
      !! DF1
      Call DGemm_('t','N',NLEV,NLEV,NLEV,
     *            1.0D+00,TRF,NLEV,DF1,NLEV,
     *            0.0D+00,Work(LWRK1),NLEV)
      Call DGemm_('N','n',NLEV,NLEV,NLEV,
     *            1.0D+00,Work(LWRK1),NLEV,TRF,NLEV,
     *            0.0D+00,DF1,NLEV)

      write(6,*) "after transformation"
      write(6,*) "DG1"
      call sqprt(dg1,nlev)
      write(6,*) "DF1"
      call sqprt(df1,nlev)
C
      CALL GETMEM('WRK1','FREE','REAL',LWRK1,NLEV*NLEV)
C
      Return
C
      End Subroutine DerFG3_Trf
C
C-----------------------------------------------------------------------
C
      SUBROUTINE H0DIAG_CASPT2_der(ISYCI,DIAG,NOW,IOW,CLAG)
      IMPLICIT REAL*8 (A-H,O-Z)
C INPUT ARRAYS:

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
      DIMENSION DIAG(MXCI),NOW(2,NSYM,NMIDV),IOW(2,NSYM,NMIDV),
     *          CLAG(MXCI)

C PURPOSE: FORM AN ARRAY OF DIAGONAL HAMILTONIAN MATRIX ELEMENTS
C FOR THE SPECIFIED TOTAL SYMMETRY ISYCI

      CALL QENTER('H0DIAGDER')
      IEMU=1
      DO MV=1,NMIDV
        DO ISYUP=1,NSYM
          NUP=NOW(1,ISYUP,MV)
          IF(NUP.EQ.0) GOTO 30
          ISYDWN=MUL(ISYUP,ISYCI)
          NDWN=NOW(2,ISYDWN,MV)
          IF(NDWN.EQ.0) GOTO 30
          ICS=LICASE+IOW(1,ISYUP,MV)
          JCS=LICASE+IOW(2,ISYDWN,MV)
          NC=NUP*NDWN
          CALL DIELMV_der(IWORK(ICS),IWORK(JCS),NUP,NDWN,
     *                    DIAG(IEMU),CLAG(IEMU))
          IEMU=IEMU+NC
  30      CONTINUE
        END DO
      END DO
      CALL QEXIT('H0DIAGDER')
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DIELMV_der(ICASE,JCASE,NUP,NDWN,EMU,CLAG)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ICASE(*),JCASE(*)
      DIMENSION EMU(NUP,NDWN),CLAG(NUP,NDWN)
#include "pt2_guga.fh"

      DO I=1,NUP
        II=NIPWLK*(I-1)
        SUM=0.0D0
        DO LV1=MIDLEV+1,NLEV,15
          II=II+1
          IC=ICASE(II)
          DO LEV=LV1,MIN(LV1+14,NLEV)
            IC1=IC/4
            ISTEP=IC-4*IC1
            IOC=(ISTEP+1)/2
            SUM=SUM+DBLE(IOC)*ETA(LEV)
            IC=IC1
          END DO
        END DO
        DO J=1,NDWN
C         EMU(I,J)=EMU(I,J)+SUM
          CLAG(I,J) = CLAG(I,J) + SUM*EMU(I,J)
        END DO
      END DO
C THEN THE LOWER HALF:
      DO I=1,NDWN
        II=NIPWLK*(I-1)
        SUM=0.0D0
        DO LV1=1,MIDLEV,15
        II=II+1
        IC=JCASE(II)
        DO LEV=LV1,MIN(LV1+14,MIDLEV)
          IC1=IC/4
          ISTEP=IC-4*IC1
          IOC=(ISTEP+1)/2
          SUM=SUM+DBLE(IOC)*ETA(LEV)
          IC=IC1
          END DO
        END DO
        DO J=1,NUP
C         EMU(J,I)=EMU(J,I)+SUM
          CLAG(I,J) = CLAG(I,J) + SUM*EMU(I,J)
        END DO
      END DO

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SIGMA1_CP2ld(IP,IQ,CPQ,ISYCI,CI,SGM,NOCSF,IOCSF,NOW,
     &                 IOW,NOCP,IOCP,ICOUP,VTAB,MVL,MVR,dcont,clag)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NOCSF(NSYM,NMIDV,NSYM),IOCSF(NSYM,NMIDV,NSYM)
      DIMENSION NOW(2,NSYM,NMIDV),IOW(2,NSYM,NMIDV)
      DIMENSION NOCP(MXEO,NSYM,NMIDV),IOCP(MXEO,NSYM,NMIDV)
      DIMENSION VTAB(NVTAB)
      INTEGER ICOUP
      DIMENSION ICOUP(3,NICOUP)
      DIMENSION MVL(NMIDV,2),MVR(NMIDV,2)
      DIMENSION CI(*),SGM(*)

      DIMENSION CLAG(*)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
      INTRINSIC MOD
*****************************************************************
*  GIVEN ACTIVE LEVEL INDICES IP AND IQ, AND INPUT CI ARRAYS
*  CI AND SGM, THIS ROUTINE ADDS TO SGM THE RESULT OF ACTING ON
*  CI WITH THE NUMBER CPQ TIMES THE EXCITATION OPERATOR E(IP,IQ).
*  THE ADDITIONAL ENTRIES IN THE PARAMETER LIST ARE TABLES THAT
*  WERE PREPARED BY GINIT AND ITS SUBROUTINES.
*****************************************************************

      CALL QENTER('SIGMA1')

      IF(ABS(CPQ).LT.1.0D-12) GOTO 999
      iLUW=0
C SYMMETRY OF ORBITALS:
      ISYP=ISM(IP)
      ISYQ=ISM(IQ)
      ISYPQ=MUL(ISYP,ISYQ)
C SYMMETRY OF SIGMA ARRAY:
      ISYSGM=MUL(ISYPQ,ISYCI)

      IF(IP.GT.IQ) THEN
C THEN THIS IS AN EXCITING OPERATOR.
        IF(IP.LE.MIDLEV) GOTO 1600
        IF(IQ.GT.MIDLEV) GOTO 1700
        GOTO 1800
      ELSE IF(IQ.GT.IP) THEN
C THEN THIS IS A DEEXCITING OPERATOR.
        IF(IQ.LE.MIDLEV) GOTO 1300
        IF(IP.GT.MIDLEV) GOTO 1400
        GOTO 1500
      END IF
      IF(IP.GT.MIDLEV) GOTO 1200

C SPECIAL CASE: WEIGHT OPERATOR, IP=IQ.
C IP=IQ < MIDLEV.
      DO 100 MVSGM=1,NMIDV
        DO 101 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 101
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          NDWNSG=NOW(2,ISYDSG,MVSGM)
          IOLW=IOW(2,ISYDSG,MVSGM)
          IPSHFT=2*(IP-1)
          LLW=LICASE+IOLW-NIPWLK+IPSHFT/30
cVV next statement is a cheat to forbid optimization
          if(iLUW.eq.666) iLUW=iLUW+LLW
cVV
          IPSHFT=MOD(IPSHFT,30)
          IPPOW=2**IPSHFT
          DO 102 J=1,NDWNSG
            JC=IWORK(LLW+J*NIPWLK)
CPAM96           ICS=IAND(ISHFT(JC,-IPSHFT),3)
            ICS=MOD(JC/IPPOW,4)
            IF(ICS.EQ.0) GOTO 102
            X=CPQ*DBLE((1+ICS)/2)
            JSTA=ISGSTA+NUPSG*(J-1)
            CALL DAXPY_(NUPSG,X,CI(JSTA),1,SGM(JSTA),1)
            CALL DAXPY_(NUPSG,X*DCONT,CI(JSTA),1,CLAG(JSTA),1)
 102      CONTINUE
 101    CONTINUE
 100  CONTINUE

      GOTO 999

 1200 CONTINUE
C SPECIAL CASE: WEIGHT OPERATOR, IP=IQ.
C IP=IQ>MIDLEV
      DO 200 MVSGM=1,NMIDV
        DO 201 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 201
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          NDWNSG=NOW(2,ISYDSG,MVSGM)
          IOUW=IOW(1,ISYUSG,MVSGM)
          IPSHFT=2*(IP-1-MIDLEV)
          LUW=LICASE+IOUW-NIPWLK+IPSHFT/30
cVV next statement is a cheat to forbid optimization
          if(iLUW.eq.666) iLUW=iLUW+LUW
cVV
          IPSHFT=MOD(IPSHFT,30)
          IPPOW=2**IPSHFT
          DO 202 I=1,NUPSG
            IC=IWORK(LUW+I*NIPWLK)
CPAM96            ICS=IAND(ISHFT(IC,-IPSHFT),3)
            ICS=MOD(IC/IPPOW,4)
            IF(ICS.EQ.0) GOTO 202
            X=CPQ*DBLE((1+ICS)/2)
            ISTA=ISGSTA-1+I
            CALL DAXPY_(NDWNSG,X,CI(ISTA),NUPSG,SGM(ISTA),NUPSG)
            CALL DAXPY_(NDWNSG,X*DCONT,CI(ISTA),NUPSG,CLAG(ISTA),NUPSG)
 202      CONTINUE
 201    CONTINUE
 200  CONTINUE

      GOTO 999

 1300 CONTINUE
C DEEXCITING OPERATOR, IP<IQ<=MIDLEV.
      DO 300 MVSGM=1,NMIDV
        DO 301 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 301
          ISYDSG=MUL(ISYUSG,ISYSGM)
          ISYDC=MUL(ISYPQ,ISYDSG)
          NDWNC=NOW(2,ISYDC,MVSGM)
          IF(NDWNC.EQ.0) GOTO 301
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          IOC=IOCSF(ISYUSG,MVSGM,ISYCI)
          INDEO=2*NLEV+(IQ*(IQ-1))/2+IP
          NCP=NOCP(INDEO,ISYDSG,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYDSG,MVSGM)
C CASE IS: LOWER HALF, DEEXCITE:
            CALL DEX1_CP2 (CPQ,NUPSG,CI(IOC+1),SGM(ISGSTA),
     &             NCP,ICOUP(1,LICP),VTAB)
            CALL DEX1_CP2 (CPQ*DCONT,NUPSG,
     *             CI(IOC+1),SGM(ISGSTA),
     &             NCP,ICOUP(1,LICP),VTAB)
          END IF
 301    CONTINUE
 300  CONTINUE

      GOTO 999

 1400 CONTINUE
C DEEXCITING OPERATOR, MIDLEV<IP<IQ
      DO 400 MVSGM=1,NMIDV
        DO 401 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 401
          ISYUC=MUL(ISYPQ,ISYUSG)
          NUPC=NOW(1,ISYUC,MVSGM)
          IF (NUPC.EQ.0) GOTO 401
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          NDWNSG=NOW(2,ISYDSG,MVSGM)
          IOC=IOCSF(ISYUC,MVSGM,ISYCI)
          INDEO=2*NLEV+(IQ*(IQ-1))/2+IP
          NCP=NOCP(INDEO,ISYUSG,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYUSG,MVSGM)
C CASE IS: UPPER HALF, DEEXCITE:
            CALL DEX2_CP2 (CPQ,NDWNSG,NUPC,CI(IOC+1),NUPSG,SGM(ISGSTA),
     &             NCP,ICOUP(1,LICP),VTAB)
          END IF
 401    CONTINUE
 400  CONTINUE

      GOTO 999

 1500 CONTINUE
C DEEXCITING CASE, IP<=MIDLEV<IQ.
C ALLOCATE TEMPORARY WORK AREA:
      DO 500 MVSGM=1,NMIDV
        MV4=MVR(MVSGM,1)
        MV5=MVR(MVSGM,2)
        DO 501 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 501
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          ISYUC=MUL(ISYQ,ISYUSG)
          ISYDC=MUL(ISYP,ISYDSG)
          IF(MV4.EQ.0) GOTO 499
          NUPC=NOW(1,ISYUC,MV4)
          IF(NUPC.EQ.0) GOTO 499
          NDWNC=NOW(2,ISYDC,MV4)
          IF(NDWNC.EQ.0) GOTO 499
          INDEO=IQ
          NCP=NOCP(INDEO,ISYUSG,MVSGM)
          IF(NCP.EQ.0) GOTO 499
          NTMP=NUPSG*NDWNC
          CALL FZERO(WORK(LSGTMP),NTMP)
          LICP=1+IOCP(INDEO,ISYUSG,MVSGM)
          IOC=IOCSF(ISYUC,MV4,ISYCI)
C CASE IS: UPPER HALF, DEEXCITE:
          CALL DEX2_CP2 (CPQ,NDWNC,NUPC,CI(IOC+1),NUPSG,WORK(LSGTMP),
     &             NCP,ICOUP(1,LICP),VTAB)
          INDEO=IP
          NCP=NOCP(INDEO,ISYDSG,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYDSG,MVSGM)
C CASE IS: LOWER HALF, DEEXCITE:
            X=1.0D00
            CALL DEX1_CP2 (X,NUPSG,WORK(LSGTMP),SGM(ISGSTA),
     &               NCP,ICOUP(1,LICP),VTAB)
          END IF
 499      CONTINUE
          IF(MV5.EQ.0) GOTO 501
          NUPC=NOW(1,ISYUC,MV5)
          IF(NUPC.EQ.0) GOTO 501
          NDWNC=NOW(2,ISYDC,MV5)
          IF(NDWNC.EQ.0) GOTO 501
          INDEO=NLEV+IQ
          NCP=NOCP(INDEO,ISYUSG,MVSGM)
          IF(NCP.EQ.0) GOTO 501
          NTMP=NUPSG*NDWNC
          CALL FZERO(WORK(LSGTMP),NTMP)
          LICP=1+IOCP(INDEO,ISYUSG,MVSGM)
          IOC=IOCSF(ISYUC,MV5,ISYCI)
C CASE IS: UPPER HALF, DEEXCITE:
          CALL DEX2_CP2 (CPQ,NDWNC,NUPC,CI(IOC+1),NUPSG,WORK(LSGTMP),
     &             NCP,ICOUP(1,LICP),VTAB)
          INDEO=NLEV+IP
          NCP=NOCP(INDEO,ISYDSG,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYDSG,MVSGM)
C CASE IS: LOWER HALF, DEEXCITE:
            X=1.0D00
            CALL DEX1_CP2 (X,NUPSG,WORK(LSGTMP),SGM(ISGSTA),
     &               NCP,ICOUP(1,LICP),VTAB)
          END IF
 501    CONTINUE
 500  CONTINUE

      GOTO 999

 1600 CONTINUE
C EXCITING CASE, IQ<IP<=MIDLEV.
      DO 600 MVSGM=1,NMIDV
        DO 601 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 601
          ISYDSG=MUL(ISYUSG,ISYSGM)
          ISYDC=MUL(ISYPQ,ISYDSG)
          NDWNC=NOW(2,ISYDC,MVSGM)
          IF(NDWNC.EQ.0) GOTO 601
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          IOC=IOCSF(ISYUSG,MVSGM,ISYCI)
          INDEO=2*NLEV+(IP*(IP-1))/2+IQ
          NCP=NOCP(INDEO,ISYDC,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYDC,MVSGM)
C CASE IS: LOWER HALF, EXCITE:
            CALL EXC1_CP2 (CPQ,NUPSG,CI(IOC+1),SGM(ISGSTA),
     &             NCP,ICOUP(1,LICP),VTAB)
          END IF
 601    CONTINUE
 600  CONTINUE

      GOTO 999

 1700 CONTINUE
C EXCITING CASE, MIDLEV<IQ<IP
      DO 700 MVSGM=1,NMIDV
        DO 701 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 701
          ISYUC=MUL(ISYPQ,ISYUSG)
          NUPC=NOW(1,ISYUC,MVSGM)
          IF (NUPC.EQ.0) GOTO 701
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          NDWNSG=NOW(2,ISYDSG,MVSGM)
          IOC=IOCSF(ISYUC,MVSGM,ISYCI)
          INDEO=2*NLEV+(IP*(IP-1))/2+IQ
          NCP=NOCP(INDEO,ISYUC,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYUC,MVSGM)
C CASE IS: UPPER HALF, EXCITE:
            CALL EXC2_CP2 (CPQ,NDWNSG,NUPC,CI(IOC+1),NUPSG,SGM(ISGSTA),
     &             NCP,ICOUP(1,LICP),VTAB)
          END IF
 701    CONTINUE
 700  CONTINUE

      GOTO 999

 1800 CONTINUE
C EXCITING CASE, IQ<=MIDLEV<IP
C ALLOCATE TEMPORARY WORK AREA:
      DO 800 MVSGM=1,NMIDV
        MV1=MVL(MVSGM,2)
        MV2=MVL(MVSGM,1)
        IF((MV1.EQ.0).AND.(MV2.EQ.0)) GOTO 800
        DO 801 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 801
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          ISYUC=MUL(ISYP,ISYUSG)
          ISYDC=MUL(ISYQ,ISYDSG)
          IF(MV2.EQ.0) GOTO 799
          NUPC=NOW(1,ISYUC,MV2)
          IF(NUPC.EQ.0) GOTO 799
          NDWNC=NOW(2,ISYDC,MV2)
          IF(NDWNC.EQ.0) GOTO 799
          INDEO=IP
          NCP=NOCP(INDEO,ISYUC,MV2)
          IF(NCP.EQ.0) GOTO 799
          NTMP=NUPSG*NDWNC
          CALL FZERO(WORK(LSGTMP),NTMP)
          LICP=1+IOCP(INDEO,ISYUC,MV2)
          IOC=IOCSF(ISYUC,MV2,ISYCI)
C CASE IS: UPPER HALF, EXCITE:
          CALL EXC2_CP2 (CPQ,NDWNC,NUPC,CI(IOC+1),NUPSG,WORK(LSGTMP),
     &             NCP,ICOUP(1,LICP),VTAB)
          INDEO=IQ
          NCP=NOCP(INDEO,ISYDC,MV2)
          IF(NCP.EQ.0) GOTO 799
          LICP=1+IOCP(INDEO,ISYDC,MV2)
C CASE IS: LOWER HALF, EXCITE:
          X=1.0D00
          CALL EXC1_CP2 (X,NUPSG,WORK(LSGTMP),SGM(ISGSTA),
     &               NCP,ICOUP(1,LICP),VTAB)
  799     CONTINUE
          IF(MV1.EQ.0) GOTO 801
          NUPC=NOW(1,ISYUC,MV1)
          IF(NUPC.EQ.0) GOTO 801
          NDWNC=NOW(2,ISYDC,MV1)
          IF(NDWNC.EQ.0) GOTO 801
          INDEO=NLEV+IP
          NCP=NOCP(INDEO,ISYUC,MV1)
          IF(NCP.EQ.0) GOTO 801
          NTMP=NUPSG*NDWNC
          CALL FZERO(WORK(LSGTMP),NTMP)
          LICP=1+IOCP(INDEO,ISYUC,MV1)
          IOC=IOCSF(ISYUC,MV1,ISYCI)
C CASE IS: UPPER HALF, EXCITE:
          CALL EXC2_CP2 (CPQ,NDWNC,NUPC,CI(IOC+1),NUPSG,WORK(LSGTMP),
     &             NCP,ICOUP(1,LICP),VTAB)
          INDEO=NLEV+IQ
          NCP=NOCP(INDEO,ISYDC,MV1)
          IF(NCP.EQ.0) GOTO 801
          LICP=1+IOCP(INDEO,ISYDC,MV1)
C CASE IS: LOWER HALF, EXCITE:
          X=1.0D00
          CALL EXC1_CP2 (X,NUPSG,WORK(LSGTMP),SGM(ISGSTA),
     &               NCP,ICOUP(1,LICP),VTAB)
 801  CONTINUE
 800  CONTINUE

 999  CONTINUE
      CALL QEXIT('SIGMA1')
      RETURN
      END
C
      Subroutine DFG3Cor(nLev,iT,iU,iV,iX,iY,iZ,
     *                   DF1,DF2,DF3,DG1,DG2,DG3)
C
      Use SUPERINDEX
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
C
      !! DF3 and DG3 is scalar
      Dimension DF1(nLev,nLev),DF2(nLev,nLev,nLev,nLev),
     *          DG1(nLev,nLev),DG2(nLev,nLev,nLev,nLev)
C
      iST=IASYM(iT)
      iSU=IASYM(iU)
      iSV=IASYM(iV)
      iSX=IASYM(iX)
      iSY=IASYM(iY)
      iSZ=IASYM(iZ)
      iTUVs=MUL(IST,MUL(ISU,ISV))
      iXYZs=MUL(ISX,MUL(ISY,ISZ))
      If (iTUVs.ne.iXYZs) Return
C
      iTU = iT+nLev*(iU-1)
      iVX = iV+nLev*(iX-1)
      iYZ = iY+nLev*(iZ-1)
C
      SYM = MUL(IASYM(iT),MUL(IASYM(iU),IASYM(iV)))
      If (jSYM.eq.iSYM) THEN
        Call DFG3CorX(iT,iU,iV,iX,iY,iZ)
      End If
C
      If (iTU.eq.iVX.and.iVX.eq.iYZ) go to 300
      If (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 200
C
      jSYM = MUL(IASYM(iV),MUL(IASYM(iX),IASYM(iT)))
      If (jSYM.EQ.iSYM) THEN
        Call DFG3CorX(iV,iX,iT,iU,iY,iZ)
      ENDIF
      jSYM=MUL(IASYM(iT),MUL(IASYM(iZ),IASYM(iY)))
      IF (jSYM.EQ.iSYM) THEN
        Call DFG3CorX(iT,iZ,iY,iX,iV,iU)
      ENDIF
      jSYM=MUL(IASYM(iY),MUL(IASYM(iU),IASYM(iV)))
      IF (jSYM.EQ.iSYM) THEN
        Call DFG3CorX(iY,iU,iV,iZ,iT,iX)
      ENDIF
C
 200  CONTINUE
C
      jSYM = MUL(IASYM(iV),MUL(IASYM(iZ),IASYM(iY)))
      IF (jSYM.EQ.iSYM) THEN
        Call DFG3CorX(iV,iZ,iY,iU,iT,iX)
      ENDIF
      jSYM=MUL(IASYM(iY),MUL(IASYM(iX),IASYM(iV)))
      IF (jSYM.EQ.iSYM) THEN
        Call DFG3CorX(iY,iX,iT,iZ,iV,iU)
      ENDIF
C
 300  CONTINUE
C
      If (iV.eq.iU.and.iT.eq.iX.and.iY.eq.iZ) Return
      If (iV.eq.iU.and.iT.eq.iZ.and.iX.eq.iY) Return
      If (iX.eq.iT.and.iV.eq.iZ.and.iU.eq.iY) Return
      If (iZ.eq.iY.and.iT.eq.iU.and.iX.eq.iV) Return
C
      jSYM = MUL(IASYM(iX),MUL(IASYM(iV),IASYM(iU)))
      IF (jSYM.EQ.iSYM) THEN
        Call DFG3CorX(iX,iV,iU,iT,iZ,iY)
      ENDIF
C
      If (iTU.eq.iVX.and.iVX.eq.iYZ) Return
      If (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 400
C
      jSYM=MUL(IASYM(iU),MUL(IASYM(iT),IASYM(iX)))
      IF (jSYM.EQ.iSYM) THEN
        Call DFG3CorX(iU,iT,iX,iV,iZ,iY)
      ENDIF
      jSYM=MUL(IASYM(iX),MUL(IASYM(iY),IASYM(iZ)))
      IF (jSYM.EQ.iSYM) THEN
        Call DFG3CorX(iX,iY,iZ,iT,iU,iV)
      ENDIF
      jSYM=MUL(IASYM(iZ),MUL(IASYM(iV),IASYM(iU)))
      IF (jSYM.EQ.iSYM) THEN
        Call DFG3CorX(iZ,iV,iU,iY,iX,iT)
      ENDIF
C
 400   CONTINUE
C
      jSYM=MUL(IASYM(iU),MUL(IASYM(iY),IASYM(iZ)))
      IF (jSYM.EQ.iSYM) THEN
        Call DFG3CorX(iU,iY,iZ,iV,iX,iT)
      ENDIF
      jSYM=MUL(IASYM(iZ),MUL(IASYM(iT),IASYM(iX)))
      IF (jSYM.EQ.iSYM) THEN
        Call DFG3CorX(iZ,iT,iX,iY,iU,iV)
      ENDIF
C
      Return
C
      Contains

      Subroutine DFG3CorX(iA,iB,iC,iD,iE,iF)
C
      Implicit Real*8 (A-H,O-Z)
C
      If (iB.eq.iC .and. iD.eq.iE) Then
        DG1(iA,iF) = DG1(iA,iF) - DG3
      End If
C
      If (iB.eq.iC) Then
        DG2(iA,iE,iD,iF) = DG2(iA,iE,iD,iF) - DG3
      End If
C
      If (iD.eq.iE) Then
        DG2(iA,iC,iB,iF) = DG2(iA,iC,iB,iF) - DG3
      End If
C
      If (iB.eq.iE) Then
        DG2(iC,iA,iD,iF) = DG2(iC,iA,iD,iF) - DG3
      End If
C
      Return
C
      End Subroutine DFG3CorX
C
      End Subroutine DFG3Cor
      subroutine mktg33(ci)
      implicit real*8 (a-h,o-z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
      dimension ci(*)
      NTG1=NASHT**2
      NTG2=NASHT**4
      NTG3=(NTG1*(NTG1+1)*(NTG1+2))/6
      NTG1=MAX(1,NTG1)
      NTG2=MAX(1,NTG2)
      NTG3=MAX(1,NTG3)
      CALL GETMEM('TG1','ALLO','REAL',LTG1,NTG1)
      CALL GETMEM('TG2','ALLO','REAL',LTG2,NTG2)
      CALL GETMEM('TG3','ALLO','REAL',LTG3,NTG3)
      WORK(LTG1)=0.0D0
      WORK(LTG2)=0.0D0
      WORK(LTG3)=0.0D0
      CALL MKTG3(1,1,CI,CI,OVL,
     &           WORK(LTG1),WORK(LTG2),NTG3,WORK(LTG3))
      CALL GETMEM('TG1','FREE','REAL',LTG1,NTG1)
      CALL GETMEM('TG2','FREE','REAL',LTG2,NTG2)
      CALL GETMEM('TG3','FREE','REAL',LTG3,NTG3)
      end subroutine
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
      SUBROUTINE DERSPE(DG1,DDG2,DG3,DF1,DF2,DF3,idxG3,
     *                  DEPSA,G1,G2,G3,F1,F2,F3)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION G1(NASHT,NASHT),G2(NASHT,NASHT,NASHT,NASHT),G3(*)
      DIMENSION F1(NASHT,NASHT),F2(NASHT,NASHT,NASHT,NASHT),F3(*)
      DIMENSION DG1(NASHT,NASHT),DG2(NASHT,NASHT,NASHT,NASHT),DG3(*)
      DIMENSION DF1(NASHT,NASHT),DF2(NASHT,NASHT,NASHT,NASHT),DF3(*)
      DIMENSION DEPSA(NASHT,NASHT)
      INTEGER*1 idxG3(6,*)
      INTEGER I1
      PARAMETER (I1=KIND(idxG3))
C SPECIAL-CASE ROUTINE. DELIVERS G AND F MATRICES FOR A HIGH-SPIN
C OR CLOSED-SHELL SCF CASE.
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"

#include "para_info.fh"
      LOGICAL RSV_TSK

      CALL QENTER('SPECIAL')

C     CALL DCOPY_(NG1,[0.0D0],0,G1,1)
C     CALL DCOPY_(NG2,[0.0D0],0,G2,1)
C     CALL DCOPY_(NG3,[0.0D0],0,G3,1)
C     CALL DCOPY_(NG1,[0.0D0],0,F1,1)
C     CALL DCOPY_(NG2,[0.0D0],0,F2,1)
C     CALL DCOPY_(NG3,[0.0D0],0,F3,1)

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

      CALL QEXIT('SPECIAL')
      RETURN
      END
