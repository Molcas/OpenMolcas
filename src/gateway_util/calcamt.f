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
      Subroutine CalcAMt (iOpt,LUQRP,MPLbl,lMax,iSRShll,
     &                    nProj,iCoShll,rcharge)
************************************************************************
*                                                                      *
*...       calculates the non-diagonal spectral representation         *
*          A matrix for an atom. Note that its signs is such           *
*          that the spectral representation must be ADDED to the       *
*          one-electron hamiltonian.                                   *
*                                                                      *
*     Internal matrices fixed the maximum number of primitives         *
*     per symmetry to 'maxprim'                                        *
*                                                                      *
************************************************************************
      Use Basis_Info
      Implicit Real*8 (A-H,O-Z)
      External Agin, Ovlmp, Vexch, Vqr
#include "relmp.fh"
#include "stdalloc.fh"
      Character*20 MPLbl

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
      PI=2D0*ACOS(0D0)
      iprint=0
      Call Agin
      call dcopy_(2*MaxPrim**2,[0.0D0],0,Corek,1)
c
c     calculate relativistic integrals if needed
      lpq=0
      do i=1,lmax+1
         nnexp=Shells(isrshll+i-1)%nExp
         lpq=lpq+nnexp*(nnexp+1)/2
      enddo
      if(4*lpq.gt.Nrel) then
         write(6,*) ' problem. Nrel and lpq are', nrel,lpq
         write(6,*) ' The dimension of rel must somehow be increased.'
         Call Abend
      endif
#ifdef _DEBUGPRINT_
      write(6,*) ' basis:', (Shells(isrshll+i-1)%nExp,i=1,lmax+1)
      do i=1,lmax+1
         nnExp=Shells(isrshll+i-1)%nExp
         write(6,*) ' number of exponents', nnExp
         write(6,*) ' exponents, symmetry', i
         do j=1,nnExp
            write(6,*) Shells(iSRShll+i-1)%Exp(j)
         enddo
      enddo
#endif
      Do 1000 lP1=1,lMax+1
        iSRSh=iSRShll+lP1-1
        nP=Shells(iSRSh)%nExp
        If (np.gt.maxprim) Then
          Write (6,*) 'CalcAMt: np.gt.maxprim',np,maxprim
          Write (6,*) 'Abend: Increase MaxPrim !'
          Call Abend
        End If
        N=lP1
        LAM=lP1
        If (nP.le.0) Go To 1000
C
        call dcopy_(nRel,[0.0D0],0,Rel,1)
        If (iAnd(iOpt,iMVPot).ne.0 .and.
     &      iAnd(iOpt,iDWPot).ne.0 ) Then
C...      Mass-velocity and/or Darwin potentials
          Call Vqr(LUQRP,MPLbl,lP1,Shells(iSRSh)%Exp,nP,rel)
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
         call oeisg(rel,srel,trel,urel,Shells(iSRSh)%Exp,
     &   rCharge,mx100,lp1,Shells(iSRSh)%nExp,unrel,tnrel,hcorr,iprint,
     &   VEXTT,PVPT,EVN1,EVN2,RE1R,AUXI,W1W1,W1E0W1)
         nmat=(Shells(iSRSh)%nExp*(Shells(iSRSh)%nExp+1))/2
         if(iprint.ge.10) then
            write(6,*) ' relativistic integrals'
            write(6,12) (hcorr(i),i=1,nmat)
 12      format(4d19.12)
         endif
         call dcopy_(nRel,[0.0D0],0,Rel,1)
      endif
C
C...    Overlap and, if neccesary, exchange.
        IJ=0
        DO 101 I=1,NP
          ZI=Shells(iSRSh)%Exp(I)
          DO 102 J=1,I
            IJ=IJ+1
            ZJ=Shells(iSRSh)%Exp(J)
            COREK(I,J,1)=rel(ij)
            If (iAnd(iOpt,iNPPot).ne.0) Then
               corek(i,j,2)=hcorr(ij)
            End If
            If (iAnd(iOpt,iExch).ne.0) Then
C...          minus exchange potential
              AuxLs=VExch(ZI,N,ZJ,N,LAM,nProj,iCoShll)
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
        PreFac = sqrt((2D0/PI)**3) * 4D0**(lP1-1)
*
        Call mma_Allocate(Shells(ISRSh)%Akl,NP,NP,2,Label='Akl')
        Shells(ISRSh)%nAkl=NP
        DO iq = 1, 2
           DO I=1,NP
              ZI=Shells(iSRSh)%Exp(I)
              DO L=1,NP
                 rel(L)=0.D0
                 DO K=1,NP
                    rel(L)=rel(L)+OVL(I,K)*COREK(K,L,iq)
                 END DO
              END DO
              DO J=1,I
                 ZJ=Shells(iSRSh)%Exp(J)
                 ADUM=0.D0
                 DO L=1,NP
                    ADUM=ADUM+rel(L)*OVL(L,J)
                 END DO
*                in MOLCAS3:X1
*                multiply by the (radial) normalization constants
*                of the primitives i and j, so that the spectral
*                representation coeffients correspond to the
*                (radially) unnormalized primitives.
                 ADUM=ADUM*PreFac*sqrt(sqrt( (ZI*ZJ)**(3+2*(lP1-1)) ))
                 Shells(iSRSh)%Akl(I,J,iq)=ADUM
                 Shells(iSRSh)%Akl(J,I,iq)=ADUM
              END DO
           END DO
           IJAM0=IJAM0+NP**2
        END DO
*
1000  CONTINUE
      RETURN
      END
