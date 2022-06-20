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
      SubRoutine Precaaa(iC,is,js,nd,ir,rOut,nbai,nbaj,
     &                   focki,fock,sign,
     &                   Scr,nScr,ActInt)
************************************************************************
*                                                                      *
*                                        [2]                           *
*   Calculates the diagonal submatrix of E    that couple              *
*                                                                      *
*   kappa                with   kappa                for a             *
*        kactive,kactive            kactive,kactive                    *
*                                                                      *
*   single active index.                                               *
*   Used for preconditioner.                                           *
*                                                                      *
*   See Olsen,Yeager, Joergensen:                                      *
*    "Optimization and characterization of an MCSCF state"             *
*                                                                      *
*     Eq. C.12e                                                        *
*                                                                      *
*   Called by prec                                                     *
*                                                                      *
*   ib,is       :       active index for the submatrix                 *
*   js          :       symmetry of inactive,inactive                  *
*   rOut        :       Submatrix                                      *
*                                                                      *
************************************************************************
      use Arrays, only: G1t, G2t
      Implicit Real*8(a-h,o-z)
#include "Input.fh"
#include "Pointers.fh"
      Real*8 Fock(nbaj,nbaj),Focki(nbaj,nbaj)
      Real*8 rout(nd*(nd+1)/2), Scr(nScr)
      Real*8 ActInt(ntAsh,ntAsh,ntAsh,ntAsh)
*                                                                      *
************************************************************************
*                                                                      *
      iTri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
      iTri1(i,j)=nTri-itri(nd-Min(i,j)+1,nd-Min(i,j)+1)
     &          +Max(i,j)-Min(i,j)+1
*                                                                      *
************************************************************************
*                                                                      *
      nTri=itri(nd,nd)
      iCC = iC+nA(iS)
      ! iiC = iC+nIsh(iS)
      iA  = iC
      iAA = iCC
*
      i=itri(iAA,iCC)
      !! Construct for all active orbitals first
      Call DCopy_(ntAsh*(ntAsh+1)/2,[0.0d+00],0,Scr,1)
      Scr(i) = 1.0d+00
*                                                                      *
************************************************************************
*                                                                      *
      Do jB=1,nAsh(jS) !! index B
        jBB=jB+nA(jS)
        jjB=jB+nIsh(jS)
        Do jD=1,jB    !! index D
          jDD=jD+nA(jS)
          jjD=jD+nIsh(jS)
C         i=itri1(jjB,jjD)
          i=itri(jBB,jDD)
          Do kS=1,nSym
C           Call Coul(kS,kS,jS,jS,jB,jD,A_J,Scr)
C           Call Exch(kS,jS,kS,jS,jB,jD,A_K,Scr)
            Do jE=1,nAsh(ks) !! index E
               jEE=jE+nA(ks)
               Do jF=1,nAsh(kS) !! index F
                  jFF=jF+nA(ks)
C
                  !! first term
                  aecf=ActInt(iAA,jEE,iCC,jFF)
                  bedf=ActInt(jBB,jEE,jDD,jFF)
                  becf=ActInt(jBB,jEE,iCC,jFF)
                  aedf=ActInt(iAA,jEE,jDD,jFF)
                  rDbedf=G2t(itri(itri(jBB,jEE),itri(jDD,jFF)))
                  rDaecf=G2t(itri(itri(iAA,jEE),itri(iCC,jFF)))
                  rDaedf=G2t(itri(itri(iAA,jEE),itri(jDD,jFF)))
                  rDbecf=G2t(itri(itri(jBB,jEE),itri(iCC,jFF)))
                  Scr(i) = Scr(i) + 4.0d+00*(aecf*rDbedf+bedf*rDaecf
     *                                        -becf*rDaedf-aedf*rDbecf)
     *                     *sign
C
                  !! second term
                  acef=ActInt(iAA,iCC,jEE,jFF)
                  bdef=ActInt(jBB,jDD,jEE,jFF)
                  bcef=ActInt(jBB,iCC,jEE,jFF)
                  adef=ActInt(iAA,jDD,jEE,jFF)
                  rDbdef=G2t(itri(itri(jBB,jDD),itri(jEE,jFF)))
                  rDacef=G2t(itri(itri(iAA,iCC),itri(jEE,jFF)))
                  rDadef=G2t(itri(itri(iAA,jDD),itri(jEE,jFF)))
                  rDbcef=G2t(itri(itri(jBB,iCC),itri(jEE,jFF)))
                  Scr(i) = Scr(i) + 2.0d+00*(acef*rDbdef+bdef*rDacef
     *                                        -bcef*rDadef-adef*rDbcef)
     *                     *sign
               End Do
            End Do
          End Do
        End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      !! Work(ipG1) -> \rho_{ab}^{(1)}
      i=0 ! dummy initialize
*
      !! the remaining third and fourth terms
      Do jB=1,nAsh(jS)
         jBB=jB+nA(jS)
         jjB=jB+nIsh(jS)
         Do jD=1,jB
            jDD=jD+nA(jS)
            jjD=jD+nIsh(jS)
*
C           i=itri1(jjB,jjD)
            i=itri(jBB,jDD)
*
            rDbd = G1t(itri(jBB,jDD))
            rDac = G1t(itri(iAA,iCC))
            rDad = G1t(itri(iAA,jDD))
            rDbc = G1t(itri(jBB,iCC))
C
            !! third term
            Scr(i) = Scr(i) + sign*2.0d+00*
     *         (rDbd*Focki(iA+nIsh(iS),iC+nIsh(iS))
     *         +rDac*Focki(jB+nIsh(jS),jD+nIsh(jS))
     *         -rDad*Focki(jB+nIsh(jS),iC+nIsh(iS))
     *         -rDbc*Focki(iA+nIsh(iS),jD+nIsh(jS)))
            !! fourth term
            if (iA.eq.jD) Scr(i) = Scr(i)
     *        + sign*2.0d+00*Fock(iC+nIsh(iS),jB+nIsh(jS))
            if (jB.eq.iC) Scr(i) = Scr(i)
     *        + sign*2.0d+00*Fock(jD+nIsh(jS),iA+nIsh(iS))
            if (jB.eq.jD) Scr(i) = Scr(i)
     *        - sign*2.0d+00*Fock(iC+nIsh(iS),iA+nIsh(iS))
            if (iA.eq.iC) Scr(i) = Scr(i)
     *        - sign*2.0d+00*Fock(jD+nIsh(jS),jB+nIsh(jS))
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      !! Then, Scr -> rOut
      !! Scr is nAsh*(nAsh+1)/2 dimension
      !! rOut is iOrb-nRASx space, where x is the space iC belongs to.
C     nseq = 0
C     do i = 1, 5
C     do j = 1, 5
C     nseq =nseq + 1
C     nseq = itri(i,j)
C     a_j(i+5*(j-1)) = scr(nseq)
C     a_j(j+5*(i-1)) = scr(nseq)
C     end do
C     end do
C     call sqprt(a_j,5)
C     write (*,*) "ir = ", ir
      If (iR.eq.1) Then
        Do jB = nRs1(jS)+1, nAsh(jS)
          jBB=jB+nA(jS)
          jjB=jB-nRs1(jS)+nIsh(jS)
          Do jD = nRs1(jS)+1, jB
            jDD=jD+nA(jS)
            jjD=jD-nRs1(jS)+nIsh(jS)
            i=itri(jBB,jDD)
            j=itri1(jjB,jjD)
            rOut(j) = rOut(j) + Scr(i)
          End Do
        End Do
      Else If (iR.eq.2) Then
        Do jB = 1, nRs1(jS)+nRs3(jS)
          jBB=jB+nA(jS)
          If (jB.gt.nRs1(jS)) jBB = jBB+nRs2(jS)
          jjB=jB+nIsh(jS)
          Do jD = 1, jB
            jDD=jD+nA(jS)
            If (jD.gt.nRs1(jS)) jDD = jDD+nRs2(jS)
            jjD=jD+nIsh(jS)
            i=itri(jBB,jDD)
            j=itri1(jjB,jjD)
            rOut(j) = rOut(j) + Scr(i)
          End Do
        End Do
      Else If (iR.eq.3) Then
        Do jB = 1, nRs1(jS)+nRs2(jS)
          jBB=jB+nA(jS)
          jjB=jB+nIsh(jS)
          Do jD = 1, jB
            jDD=jD+nA(jS)
            jjD=jD+nIsh(jS)
            i=itri(jBB,jDD)
            j=itri1(jjB,jjD)
            rOut(j) = rOut(j) + Scr(i)
          End Do
        End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
C        Call Unused_integer(ir)
         Call Unused_integer(nbai)
         Call Unused_real_array(fock)
      End If
      End
C
C
C
      Subroutine Precaaa_Pre(ActInt,A_J,Scr)
C
      Implicit Real*8 (A-H,O-Z)
C
#include "Input.fh"
C#include "WrkSpc.fh"
#include "Pointers.fh"
C
      Real*8 ActInt(ntAsh,ntAsh,ntAsh,ntAsh),A_J(*),Scr(*)
C
      Do iSym = 1, nSym
        Do iA = 1, nAsh(iSym)
          iAabs = iA + nA(iSym)
          iAtot = iA + nIsh(iSym)
          Do iB = 1, nAsh(iSym) ! iA
            iBabs = iB + nA(iSym)
            iBtot = iB + nIsh(iSym)
            Do jSym = 1, iSym
              Call Coul(iSym,iSym,jSym,jSym,iAtot,iBtot,A_J,Scr)
              Do iC = 1, nAsh(jSym)
                iCabs = iC + nA(jSym)
                iCtot = iC + nIsh(jSym)
                Do iD = 1, nAsh(jSym) ! iC
                  iDabs = iD + nA(jSym)
                  iDtot = iD + nIsh(jSym)
                  Val = A_J(iCtot+nBas(jSym)*(iDtot-1))
                  ActInt(iAabs,iBabs,iCabs,iDabs) = Val
C                 ActInt(iBabs,iAabs,iCabs,iDabs) = Val
C                 ActInt(iAabs,iBabs,iDabs,iCabs) = Val
C                 ActInt(iBabs,iAabs,iDabs,iCabs) = Val
C                 ActInt(iCabs,iDabs,iAabs,iBabs) = Val
C                 ActInt(iCabs,iDabs,iBabs,iAabs) = Val
C                 ActInt(iDabs,iCabs,iAabs,iBabs) = Val
C                 ActInt(iDabs,iCabs,iBabs,iAabs) = Val
                End Do
              End Do
            End Do
          End Do
        End Do
      End Do
C
      End Subroutine Precaaa_Pre
