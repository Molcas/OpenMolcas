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
      Real*8 Function VExch (ZP,NP,ZQ,NQ,LA,nProj,iCoShll)
************************************************************************
*                                                                      *
*     VExch calculates the atomic integral                             *
*     <zp,np,la | Sum(core) Kc |zq,nq,la>                              *
*                                                                      *
************************************************************************
      Use Basis_Info
      Implicit Real*8 (A-H,O-Z)
C...  auxiliar constant pool:       ready only up to g-valence/g-core
      PARAMETER (lp1=5,lp12=lp1*lp1,lp13=(lp1*lp1+lp1)/2)
      COMMON/CONST/RCA(lp1,lp13),DFAC(lp12),KOSUU(lp13),NYU(lp1,lp13)
*
*     Define statement function
      VF(N,X)=DFAC(N)/sqrt(X)**(N+1)
*
      If (nProj.gt.4) Then
        Write (6,*) 'VExch: nProj',nProj
        Write (6,*) 'Abend: Implementation ready only up to g-core.'
        Write (6,*) '       Update common block /CONST/.'
        Call Abend
      Endif
      If (NP.gt.5.or.NQ.gt.5) Then
        Write (6,*) 'VExch: NP,NQ',NP,NQ
        Write (6,*) 'Abend: Implementation ready only up to g-valence.'
        Write (6,*) '       Update common block /CONST/.'
        Call Abend
      Endif
*
*     PI=2.D0*ACOS(0.D0)
*     PIPPI=(0.5D0/PI)**0.75D0
      L1=LA
      VPQ=VF(2*NP,ZP)*VF(2*NQ,ZQ)
*     PSMT=0.D0
*ls start adding
      Vexch=0.D0
*ls end adding
*     loop over angular momentum
      DO 55 ISIM=1,nProj+1
        iCoSh=iCoShll+ISIM-1
*       EDUM=0.5D0*DBLE(ISIM)+0.25D0
        NR=ISIM
        NS=ISIM
        L2=ISIM
        DL2=DBLE(2*(L2-1)+1)
        LMT=((L1-1)*L1)/2+L2
        IF(L1.LT.L2) LMT=((L2-1)*L2)/2+L1
        KOMAX=KOSUU(LMT)
*       loop over core orbitals of a given angular momentum
        DO 54 ICORB=1,Shells(iCoSh)%nBasis
          OrbPS=0d0
          DO 50 INU=1,KOMAX
            NUT=NYU(INU,LMT)
            NU=NUT
            RCAT=RCA(INU,LMT)*DL2
C                BEGIN KSM.
C           COMPUTATION OF K SUPER MATRIX OVER CGTO
C           NU IS UPPER INDEX NYU IN ROOTHAAN'S NOTATION
C
            IT1=NP+NR-NU-1
            IT2=NQ+NS+NU
            IT3=NQ+NS-NU-1
            IT4=NP+NR+NU
            SUMA=0.D0
            nExpon=Shells(iCoSh)%nExp
            DO K=1,nExpon
              ZR=Shells(iCoSh)%Exp(K)
              CR=Shells(iCoSh)%Cff_c(K,ICORB,2)
CLS           coeff of unnormalized non-diagonal cartesian GTF
* to be used in ecpaimp       CR=CR/(PIPPI*(4.D0*ZR)**EDUM)
              VR=VF(2*NR,ZR)
              RTT1=0.5D0*(ZP+ZR)
              DO L=1,nExpon
                ZS=Shells(iCoSh)%Exp(L)
                CS=Shells(iCoSh)%Cff_c(L,ICORB,2)
* to be used in ecpaimp          CS=CS/(PIPPI*(4.D0*ZS)**EDUM)
                VS=VF(2*NS,ZS)
                RTT2=0.5D0*(ZQ+ZS)
                RTT3=RTT1/RTT2
                RTT4=1.D0/RTT3
                CALL AUXC((IT1+1)/2,IT2,RTT3,RTT5)
                CALL AUXC((IT3+1)/2,IT4,RTT4,RTT6)
                RTT7=VF(IT1,RTT1)*VF(IT2,RTT2)*RTT5
     .              +VF(IT3,RTT2)*VF(IT4,RTT1)*RTT6
                SUMA=SUMA+RTT7*CR*CS/sqrt(VPQ*VR*VS)
              END DO
            END DO
C                END KSM.
*ls c       PSMT=PSMT+RCAT*0.797884561D0*SUMA
*ls start adding
            OrbPS=OrbPS+RCAT*0.797884561D0*SUMA
*ls end adding
  50      CONTINUE
*ls start adding
          FOcc=Shells(iCoSh)%Occ(iCOrb)
          Vexch=Vexch+2d0*OrbPS * FOcc
*ls start adding
 54     CONTINUE
*       end of loop over core orbitals of a given angular momentum
55    CONTINUE
*     end of loop over angular momentum
*ls c VEXCH= 2.D0*PSMT
C
      RETURN
      END
