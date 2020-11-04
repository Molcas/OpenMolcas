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
* Copyright (C) 2014, Mickael G. Delcey                                *
************************************************************************
      SubRoutine Preca_cho(iB,is,js,nd,ir,rOut,nbai,nbaj,
     &                   fockii,fockai,fockti,
     &                   focki,focka,fock,sign,
     &                   A_J,A_K,Scr,nScr,iAdr)
*     &                   A_J,A_K,Scr,nScr,iAdr2)
************************************************************************
*                                                                      *
*     This routine remplaces precaii, precabi and precabb              *
*     in case the new Cholesky alrgorithm is used,                     *
*     that is if only (ii|ab) and (ia|ib) integrals were computed      *
*                                                                      *
*     The code should be slightly more efficient as the integrals      *
*     are only read once and not for each distinct element             *
*                                                                      *
*     Written by M.G. Delcey, november 2014                            *
*                                                                      *
************************************************************************
      use Arrays, only: G1t
      Implicit Real*8(a-h,o-z)
#include "Input.fh"
#include "WrkSpc.fh"
#include "Pointers.fh"
      Real*8 Fock(nbaj,nbaj),FockA(nBaj,nBaj),Focki(nbaj,nbaj)
      Real*8 rout(nd*(nd+1)/2), A_J(nScr), A_K(nScr), Scr(nScr)
*                                                                      *
************************************************************************
*                                                                      *
      iTri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
      iTri1(i,j)=nTri-itri(nd-Min(i,j)+1,nd-Min(i,j)+1)
     &          +Max(i,j)-Min(i,j)+1
*                                                                      *
************************************************************************
*                                                                      *
      nO=nAsh(js)+nIsh(js)
      nTri=itri(nd,nd)
      iBB=ib+nA(is)
      iiB=ib+nish(is)
      jVert=nOrb(js)-nIsh(js)-nAsh(js)
      itAA=itri(iBB,iBB)
      i2=nD-jVert+1
*                                                                      *
************************************************************************
*    Integral contribution                                             *
*
*      iAdr=iAdr2 ! do not update iAdr2!
      iAdr=0
      Do ksym=1,nsym
*        ksym=iEor(jsym,)
!!!!!!!!!!!!!
        Do iv=1,nAsh(kSym)
          jCC=iv+nA(kSym)
          Do iu=1,iv
            jDD=iu+nA(kSym)
            Do iJK=1,2 ! 1=coulomb, 2=exchange
              factor=2.0d0
              factor2=1.0d0
              If (iJK.eq.2) factor=1.0d0 ! exchange
              If (iJK.eq.2) factor2=2.0d0 ! exchange
              Do lsym=1,nsym
                nl=nBas(lsym)
                If (nl.le.0) Go to 20
                Call DDAFILE(LuChoInt(2),2,A_J,nl**2,iAdr)

                If (lsym.eq.js) Then
                  If (iJK.eq.1) Then
                    rDens2=2.0d0*sign*
     &                Work(ipg2-1+itri(itri(jCC,jDD),itri(iBB,iBB)))
                  Else
                    rDens2=2.0d0*sign*
     &                work(ipg2-1+itri(itri(iBB,jDD),itri(jCC,iBB)))
                  EndIf
*
                  rDensaii=factor*rDens2
                  rDensabi=-factor2*rDens2
                  rDensabb=factor2*rDens2
                  rDensaiil=0
                  rDensaiiu=0
                  rDensabil=0
                  rDensabiu=0
*
                  If (iBB.eq.jDD) Then
                    rDens1=-2.0d0*G1t(itri(jCC,iBB))
                    If (jCC.eq.iBB) rdensaii=rdensaii-2.0d0*factor
                    rDensaiil=rDensaiil-factor*rDens1
                    rDensabil=rDensabil+sign*rDens1
                  ElseIf ((iBB.eq.jCC).and.(iJK.eq.2)) Then
                    rDens1=-2.0d0*G1t(itri(jDD,iBB))
                    rDensaiiu=rDensaiiu-factor*rDens1
                    rDensabiu=rDensabiu+sign*rDens1
                  EndIf
                  If ((iBB.eq.jCC).and.(iJK.eq.2)) Then
                    rDensaiil=rDensaiil-factor*
     &                       14.0d0*sign*G1t(itri(jDD,iBB))
                    rDensabil=rDensabil+
     &                       8.0d0*sign*G1t(itri(jDD,iBB))
                    If (jDD.eq.iBB) rdensaii=rdensaii+factor*14.0d0*sign
                  ElseIf ((iBB.eq.jDD).and.(iJK.eq.2)) Then
                    rDensaiiu=rDensaiiu-factor*
     &                       14.0d0*sign*G1t(itri(jCC,iBB))
                    rDensabiu=rDensabiu+
     &                       8.0d0*sign*G1t(itri(jCC,iBB))
                  EndIf
                  rDensaiiu=rDensaiiu+rDensaii
                  rDensaii=rDensaiil+rDensaii
                  rDensabiu=rDensabiu+rDensabi
                  rDensabi=rDensabil+rDensabi
                  If ((iJK.eq.1).and.(jCC.gt.jDD)) Then
                      rDensaii=rDensaii*2.0d0
                      rDensabi=rDensabi*2.0d0
                      rDensabb=rDensabb*2.0d0
                  EndIf

                  Do ii=1,nIsh(js)
*
** aii
*
                    ni=nIsh(jS)-ii+1
                    ip=iTri1(ii,ii)
                    Call DaXpY_(ni,rDensaii,
     &                       A_J((ii-1)*nl+ii),1,
     &                       rout(ip),1)
                    If ((iJK.eq.2).and.(jCC.gt.jDD)) Then
                      Call DaXpY_(ni,rDensaiiu,
     &                       A_J((ii-1)*nl+ii),nl,
     &                       rout(ip),1)
                    EndIf
*
** abi
*
                    ip=itri1(ii,nd-jVert+1)
                    Call DaXpY_(jvert,rDensabi,
     &                       A_J((ii-1)*nl+no+1),1,
     &                       rout(ip),1)
                    If ((iJK.eq.2).and.(jCC.gt.jDD)) Then
                      Call DaXpY_(jvert,rDensabiu,
     &                         A_J(no*nl+ii),nl,
     &                         rout(ip),1)
                    EndIf
                  End Do
*
** abb
*
                  Do ii=1,jvert
                    ni=jvert-ii+1
                    ip=itri1(nIsh(jS)+ii,nIsh(jS)+ii)
                    Call DaXpY_(ni,rDensabb,
     &                       A_J((nO+ii-1)*nl+no+ii),1,
     &                       rout(ip),1)
                    If ((iJK.eq.2).and.(jCC.gt.jDD)) Then
                      Call DaXpY_(ni,rDensabb,
     &                         A_J((nO+ii-1)*nl+no+ii),nl,
     &                         rout(ip),1)
                    EndIf
                  End Do
                EndIf
 20             Continue
*
              End Do
            End Do
          End Do
        End Do
      End Do
*                                                                      *
************************************************************************
*    Fock matrix contribution                                          *
*
      rFock = sign*2.0d0*Fockii + sign*2.0d0*Fockai - sign*Fockti
      rdens=sign*2.0d0*G1t(itri(ibb ,ibb))
*
**    aii
*
      Do jA=1,nIsh(jS)
         Do jB=1,jA
            i=itri1(ja,jb)
            rout(i) = rout(i) - sign*4.0d0*( Focka(jA,jB)
     &                                      +Focki(jA,jB) )
     &                        + rdens*Focki(ja,jb)
         End Do
         rout(i) = rout(i) + 2.0d0*rfock
*
**       abi
*
         ip=itri1(jA,nd-jVert+1)
         Fact=(2.0d0-2.0d0*G1t(itAA))
         Call DaxPy_(jVert,Sign*Fact,FockI(nO+1,jA),1,rOut(ip),1)
         Fact=2.0d0
         Call DaxPy_(jVert,Sign*Fact,FockA(nO+1,jA),1,rOut(ip),1)
      End Do
*
**       abb
*
      ip=iTri1(i2,i2)
      rF=sign*Fockti
      Do iI=nAsh(js)+nIsh(js)+1,nBas(js)
         rOut(ip)=rout(ip)-2.0d0*rF+rDens*FockI(iI,ii)
         ip=ip+1
         Do iJ=iI+1,Nbas(js)
            rOut(ip)=rout(ip)+rDens*FockI(iI,iJ)
            ip=ip+1
         End Do
      End Do
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_integer(ir)
        Call Unused_integer(nbai)
        Call Unused_real_array(fock)
        Call Unused_real_array(A_K)
        Call Unused_real_array(Scr)
      End If
      end
