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
      Subroutine CalcAMt (iOpt,LUQRP,MPLbl,
     &                    lMax,iSRShll,ipAkl,ip_Occ,
     &                    nProj,iCoShll,
     &                    ipExp,ipCff,nExp,nBasis,MxShll,
     &                    rcharge,DInf,nDInf)
C
C...       calculates the non-diagonal spectral representation
C          A matrix for an atom. Note that its signs is shuch
C          that the spectral representation must be ADDED to the
C          one-electron hamiltonian.
C
C     Internal matrices fixed the maximum number of primitives
C     per symmetry to 'maxprim'
C
      Implicit Real*8 (A-H,O-Z)
      External Agin, Ovlmp, Vexch, Vqr
#include "relmp.fh"
      Real*8 DInf(nDInf)
      Character*20 MPLbl
      Integer ipExp(MxShll), ipCff(MxShll), ipAkl(MxShll),
     &           nExp(MxShll), nBasis(MxShll),  ip_Occ(MxShll)

C...  working variables (change this)
      Parameter (maxprim=40)
      Parameter (mx100=100)
      Parameter (Mxlpq=(mx100*(mx100+1)/2))
      Parameter (Nrel=7*Mxlpq+5*mx100*mx100+5*mx100)
      Real*8 COREK(maxprim,maxprim,2),OVL(maxprim,maxprim),
     &       rel(Nrel),
     &       srel(maxprim*(maxprim+1)/2),
     &       trel(maxprim*(maxprim+1)/2),
     &       urel(maxprim*(maxprim+1)/2),
     &       tnrel(maxprim*(maxprim+1)/2),
     &       unrel(maxprim*(maxprim+1)/2),
     &       hcorr(maxprim*(maxprim+1)/2)
      Real*8 VEXTT(maxprim*(maxprim+1)/2),PVPT(maxprim*(maxprim+1)/2),
     &       EVN1(maxprim,maxprim),EVN2(maxprim,maxprim),
     &       RE1R(maxprim,maxprim),
     &       AUXI(maxprim,maxprim),
     &       W1W1(maxprim,maxprim),W1E0W1(maxprim,maxprim)
      Data iExch/1/,iMVPot/2/,iDWPot/4/,iNPPot/8/
*
      Call qEnter('CalcAMt')
      PI=2D0*ACOS(0D0)
      iprint=0
      Call Agin
      call dcopy_(2*MaxPrim**2,0.0D0,0,Corek,1)
c
c     calculate relativistic integrals if needed
      lpq=0
      do i=1,lmax+1
         nnexp=nexp(isrshll+i-1)
         lpq=lpq+nnexp*(nnexp+1)/2
      enddo
      if(4*lpq.gt.Nrel) then
         write(6,*) ' problem. Nrel and lpq are', nrel,lpq
         write(6,*) ' The dimension of rel must somehow be increased.'
         Call Abend
      endif
         if(iprint.ge.10) then
            lpq=0
            write(6,*) ' basis:', (nexp(isrshll+i-1),i=1,lmax+1)
            nbias=ipExp(isrshll)-1
            do i=1,lmax+1
              nnexp=nexp(isrshll+i-1)
              write(6,*) ' number of exponents', nnexp
              write(6,*) ' exponents, symmetry', i
              nbias=ipexp(isrshll+i-1)-1
              do j=1,nnexp
                 nbias=nbias+1
                 write(6,*) DInf(nbias)
              enddo
            lpq=lpq+nnexp*(nnexp+1)/2
            enddo
      endif
      Do 1000 lP1=1,lMax+1
        iSRSh=iSRShll+lP1-1
        nP=nExp(iSRSh)
        If (np.gt.maxprim) Then
          Write (6,*) 'CalcAMt: np.gt.maxprim',np,maxprim
          Write (6,*) 'Abend: Increase MaxPrim !'
          Call Abend
        End If
        N=lP1
        LAM=lP1
        If (nP.le.0) Go To 1000
C
        call dcopy_(nRel,0.0D0,0,Rel,1)
        If (iAnd(iOpt,iMVPot).ne.0 .and.
     &      iAnd(iOpt,iDWPot).ne.0 ) Then
C...      Mass-velocity and/or Darwin potentials
          Call Vqr(LUQRP,MPLbl,lP1,DInf(ipExp(iSRSh)),nP,rel)
        Else If ((iAnd(iOpt,iMVPot).ne.0 .and.
     &            iAnd(iOpt,iDWPot).eq.0 ) .or.
     &           (iAnd(iOpt,iMVPot).eq.0 .and.
     &            iAnd(iOpt,iDWPot).ne.0 )) Then
          Write (6,*) 'Mass-Velocity and Darwin potentials must be'
          Write (6,*) 'active simultaneosly.'
          Call Abend
        Endif
C
      If (iAnd(iOpt,iNPPot).ne.0) then   ! Zero
         call oeisg(rel,srel,trel,urel,DInf(ipExp(iSRSh)),
     &   rCharge,mx100,lp1,nExp(iSRSh),unrel,tnrel,hcorr,iprint,
     &   VEXTT,PVPT,EVN1,EVN2,RE1R,AUXI,W1W1,W1E0W1)
         nmat=(nExp(iSRSh)*(nExp(iSRSh)+1))/2
         if(iprint.ge.10) then
            write(6,*) ' relativistic integrals'
            write(6,12) (hcorr(i),i=1,nmat)
 12      format(4d19.12)
         endif
         call dcopy_(nRel,0.0D0,0,Rel,1)
      endif
C
C...    Overlap and, if neccesary, exchange.
        IJ=0
        DO 101 I=1,NP
          ZI=DInf(ipExp(iSRSh)+I-1)
          DO 102 J=1,I
            IJ=IJ+1
            ZJ=DInf(ipExp(iSRSh)+J-1)
            COREK(I,J,1)=rel(ij)
            If (iAnd(iOpt,iNPPot).ne.0) Then
               corek(i,j,2)=hcorr(ij)
            End If
            If (iAnd(iOpt,iExch).ne.0) Then
C...          minus exchange potential
              AuxLs=VExch(ZI,N,ZJ,N,LAM,
     &                    ipExp,ipCff,nExp,nBasis,MxShll,
     &                    nProj,iCoShll, ip_Occ,DInf,nDInf)
              COREK(I,J,1)=COREK(I,J,1)-AuxLs
            ENDIF
            OVL(I,J)=OVLMP(N,ZI,N,ZJ)
            OVL(J,I)=OVL(I,J)
            COREK(J,I,1)=COREK(I,J,1)
            COREK(J,I,2)=COREK(I,J,2)
102       CONTINUE
101     CONTINUE
C
        CALL MATINV (OVL,rel,NP,0,maxprim)
C
        IJAM0=ipAkl(iSRSh)-1
        PreFac = sqrt((2D0/PI)**3) * 4D0**(lP1-1)
*
        Do iq = 1, 2
        DO 201 I=1,NP
          ZI=DInf(ipExp(iSRSh)+I-1)
          DO 204 L=1,NP
          rel(L)=0.D0
            DO 203 K=1,NP
              rel(L)=rel(L)+OVL(I,K)*COREK(K,L,iq)
203         CONTINUE
204       CONTINUE
          DO 202 J=1,I
            ZJ=DInf(ipExp(iSRSh)+J-1)
            ADUM=0.D0
            DO 214 L=1,NP
              ADUM=ADUM+rel(L)*OVL(L,J)
214         CONTINUE
*           in MOLCAS3:X1
*           multiply by the (radial) normalization constants
*           of the primitives i and j, so that the spectral
*           representation coeffients correspond to the
*           (radially) unnormalized primitives.
            ADUM=ADUM*PreFac*sqrt(sqrt( (ZI*ZJ)**(3+2*(lP1-1)) ))
            DInf(IJAM0+(J-1)*NP+I)=ADUM
            DInf(IJAM0+(I-1)*NP+J)=ADUM
202       CONTINUE
201     CONTINUE
        IJAM0=IJAM0+NP**2
        End Do
1000  CONTINUE
      Call qExit('CalcAMt')
      RETURN
      END
