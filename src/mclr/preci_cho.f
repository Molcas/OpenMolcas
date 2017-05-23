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
      SubRoutine Preci_cho(iB,iS,jS,nd,rOut,nbai,nbaj,
     &                   fockii,fockai,fockti,
     &                   focki,focka,fock,sign,
     &                   A_J,A_K,Scr,nScr,iAdr)
************************************************************************
*                                                                      *
*     This routine remplaces preciaa, preciba and precibb              *
*     in case the new Cholesky alrgorithm is used,                     *
*     that is if only (ii|ab) and (ia|ib) integrals were computed      *
*                                                                      *
*     The code should be slightly more efficient as the integrals      *
*     are only read once and not for each distinct element             *
*                                                                      *
*     Written by M.G. Delcey, November 2014                            *
*                                                                      *
************************************************************************
      Implicit Real*8(a-h,o-z)
#include "Input.fh"
#include "Pointers.fh"
#include "standard_iounits.fh"
#include "WrkSpc.fh"
      Real*8 focki(nbaj,nbaj),fock(nbaj,nbaj),focka(nbaj,nbaj),
     &       rout(*), A_J(nScr), A_K(nScr), Scr(nScr)
*                                                                      *
************************************************************************
*                                                                      *
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
      itri1(i,j)=nTri-itri(nd-Min(i,j)+1,nd-Min(i,j)+1)
     &          +Max(i,j)-Min(i,j)+1
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
*                                                                      *
************************************************************************
*                                                                      *
      nTri=itri(nd,nd)
      nO=nAsh(jS)+nIsh(jS)
      jVert=nOrb(jS)-nAsh(jS)-nIsh(jS)
      nvirt=nOrb(jS)-nIsh(jS)
      i1=nD-jVert+1
*                                                                      *
************************************************************************
*     Integral contribution                                            *
*
      Do iJK=1,2 ! 1=coulomb, 2=exchange
        factor=-1.0d0
        If (ijK.eq.2) factor=3.0d0
*
**    Read integrals
*
        Do isym=1,nsym
          nvirt2=nBas(isym)-nIsh(isym)
          call DDAFILE(LuChoInt(1),2,A_J,nvirt2**2,iAdr)

*
**        iaa
*
          Do jC=1,nAsh(isym)
            jjC=JC+nA(iSym)
            Do jD=1,nAsh(isym)
              jjD=JD+nA(iSym)
              ip1=(jD-1)*nvirt2+jC
              aabb=A_J(ip1) ! abab if exchange
              Do jA=1,nAsh(jS)
                jjA=jA+nA(jS)
                Do jB=1,jA
                  jjB=jB+nA(jS)
                  i=itri1(jA,jB)
                  If (iJK.eq.1) Then
                    rDens1=2.0d0*sign*work(ipg2-1+
     &                          (itri(itri(jjC,jjD),itri(jjB,jjA))))
                  Else
                   rDens1=4.0d0*sign*work(ipg2-1+
     &                          (itri(itri(jjB,jjD),itri(jjC,jjA))))
                  EndIf
                  rout(i)=rout(i)+rDens1*aabb
                End Do
              End Do
            End Do
          End Do
*
          If (iSym.eq.jS) Then
            Do jA=1,nAsh(jS)
              jAA=jA+nA(jS)
              Do jB=1,jA
                jBB=jB+nA(jS)
                i=itri1(jA,jB)
                Do jC=1,nAsh(jS)
                  jCC=jC+nA(jS)
                  iBC=(jC-1)*nvirt+jB
                  BCbb=A_J(iBC) ! BbCb if exchange
                  iAC=(jC-1)*nvirt+jA
                  ACbb=A_J(iAC) ! AbCb if exchange
                  rDens1=-sign*Work(ipg1-1+itri(jAA,jCC))
                  rDens2=-sign*Work(ipg1-1+itri(jBB,jCC))
                  If (jAA.eq.jCC) rDens1=rdens1+sign
                  If (jBB.eq.jCC) rDens2=rdens2+sign
                  rout(i)=rout(i)
     &                   +2.0d0*rdens1*factor*BCbb
     &                   +2.0d0*rdens2*factor*ACbb

                End Do
              End Do
            End Do
*
**          iba
*
            Do jA=1,nAsh(jS)
              ip=itri1(ja,nd-jVert+1)
              Do jB=1,nAsh(jS)
                rDens=-sign*Work(ipg1-1+iTri(jA+nA(jS),jB+nA(jS)))
                If (jA.eq.jB) rDens=rdens+sign*2.0d0
*
                ivB=(jB-1)*nvirt + nAsh(jS) + 1
                Call DaXpY_(jVert,2.0d0*factor*rDens,
     &                     A_J(ivB),1,
     &                     rOut(ip),1)
              End Do
            End Do
*
**    ibb
*
            i=itri1(i1,i1)
            Do kB=nAsh(jS),nvirt-1
              nlB=nvirt-kb
              ilB=kB+1+nvirt*kb
              call daxpy_(nlB,sign*4.0d0*factor,A_J(ilB),
     &                   nvirt,rout(i),1)
              i=i+nlB
            End Do
          EndIf
        End Do
      End Do
*                                                                      *
************************************************************************
*    Fock matrix contribution                                          *
*
      rFock=sign*(Fockii+FockAi)
*MGD tmp
*      Go to 10
      Do jA=1,nAsh(jS)
*
**       iaa
*
         jAA=jA+nA(jS)
         jjA=jA+nIsh(js)
         Do jB=1,JA
            jBB=jB+nA(jS)
            jjB=jB+nIsh(js)
            i=itri1(jA,jB)
            rDens=Work(ipG1+itri(jbb,jAA)-1)
            rout(i)=rout(i)+Sign*(2.0d0*rdens*Fockii + ! (ib,ib)+
     &              2.0d0*(2.0d0*Focki(jjA,jjB)+
     &              2.0d0*FockA(jjA,jjB)-Fock(jjB,jjA)))
         End Do
         rout(i) = rout(i) - 4.0d0*rFock
*
**       iba
*
         ip=iTri1(ja,nAsh(js)+1)
         Call DaXpY_(jVert,sign*4.0d0,Focki(nO+1,ja+nIsh(js)),1,
     &               rout(ip),1)
         Call DaXpY_(jVert,sign*4.0d0,FockA(nO+1,ja+nIsh(js)),1,
     &               rout(ip),1)
         Call DaXpY_(jVert,-sign,     Fock (nO+1,ja+nIsh(js)),1,
     &               rout(ip),1)
      End Do
*
**       ibb
*
      i=itri1(i1,i1)-1
      Do kB=nIsh(jS)+nAsh(jS),nOrb(jS)-1
         rOut(i+1)=rout(i+1)-4.0d0*rFock
         Do lB=kb,nOrb(JS)-1
            i=i+1
            rOut(i)=rout(i)+
     &              sign*4.0d0*Focki(kb+1,lb+1)+
     &              sign*4.0d0*Focka(kb+1,lb+1)
         End Do
      End Do
* 10   Continue
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_integer(iB)
        Call Unused_integer(iS)
        Call Unused_integer(nbai)
        Call Unused_real(fockti)
        Call Unused_real_array(A_K)
        Call Unused_real_array(Scr)
      End If
      end
