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
* Copyright (C) 1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine ClrBuf(idcrr,idcrs,idcrt,ngr,
     &                  istb,jstb,kstb,lstb,
     &                  Shijij,iAnga,iCmp,iCmpa,
     &                  iShll,iShell,jShell,
     &                  iBasi,jBasj,kBask,lBasl,
     &                  Dij1,Dij2,mDij,nDij,
     &                  Dkl1,Dkl2,mDkl,nDkl,
     &                  Dik1,Dik2,mDik,nDik,
     &                  Dil1,Dil2,mDil,nDil,
     &                  Djk1,Djk2,mDjk,nDjk,
     &                  Djl1,Djl2,mDjl,nDjl,
     &                  Final,nFinal,
     &                  FckTmp,nFT,Scrtch1,nS1,Scrtch2,nS2,
     &                  Temp,nTemp,
     &                  TwoHam,nTwo,IndGrd,Index,iAO,iAOst,
     &                  iuvwx,IfG,n8,ltri,moip,nAcO,
     &                  rmoin,nmoin,ntemptot,Buffer,c,nop,din,dan,
     &                  new_fock)
************************************************************************
*                                                                      *
*       Called from: Twoel                                             *
*       takes care of the integrals                                    *
*       integrals -> fckmatrix,MO                                      *
*       in the near feature a disk based version                       *
*                                                                      *
*       Calling:   CntrDens : Gets the indexes for d1                  *
*                  MkFck : Add up the integrals on the Fock Matrix     *
*                                                                      *
*       Author: Anders Bernhardsson, Theoretical Chemistry,            *
*               University of Lund, Sweden, June '95                   *
************************************************************************
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "pso.fh"
#include "itmax.fh"
#include "info.fh"
#include "disp.fh"
#include "disp2.fh"
#include "buffer.fh"
#include "cputime.fh"
*
      Integer iAnga(4), iShll(4),iShell(4),jShell(4),
     &        jOp(6), iCmp(4),icmpa(4),
     &        nop(4),Index(3,4),iuvwx(4),
     &        moip(0:nIrrep-1),
     &        IndGrd(3,4,0:nIrrep-1),iAO(4),iAOst(4)
      Real*8 Dij1(mDij,nDij),Dkl1(mDkl,nDkl),
     &       Dik1(mDik,nDik),Dil1(mDil,nDil),
     &       Djk1(mDjk,nDjk),Djl1(mDjl,nDjl),
     &       Dij2(mDij,nDij),Dkl2(mDkl,nDkl),
     &       Dik2(mDik,nDik),Dil2(mDil,nDil),
     &       Djk2(mDjk,nDjk),Djl2(mDjl,nDjl),
     &       FckTmp(nFT),Scrtch1(nS1),Temp(nTemp),
     &       Scrtch2(nS2),TwoHam(nTwo),Final(nFinal),
     &       rmoin(nMOIN),c(12),buffer(*),din(*),dan(*)
      Logical   Shijij,n8,IfG(4),pert(0:7),ltri,new_fock
*
      nijkl=iBasi*jBasj*kBask*lBasl
      nabcd=iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4)
      Call Timing(dum1,Time,dum2,dum3)
*
      ExFac=One
*
      If (ltri) Then
      if (.not.new_fock) then
*----------------------------------------------------------------*
*
*       Get the size of the work area that should be contracted
*
*----------------------------------------------------------------*
        If (jShell(1).ge.jShell(2)) Then
         ij1 = iBasi
         ij2 = jBasj
         ij3 = iCmp(1)
         ij4 = iCmp(2)
        Else
         ij1 = jBasj
         ij2 = iBasi
         ij3 = iCmp(2)
         ij4 = iCmp(1)
        End If
        If (jShell(3).ge.jShell(4)) Then
         kl1 = kBask
         kl2 = lBasl
         kl3 = iCmp(3)
         kl4 = iCmp(4)
        Else
         kl1 = lBasl
         kl2 = kBask
         kl3 = iCmp(4)
         kl4 = iCmp(3)
        End If
        If (jShell(1).ge.jShell(3)) Then
         ik1 = iBasi
         ik2 = kBask
         ik3 = iCmp(1)
         ik4 = iCmp(3)
        Else
         ik1 = kBask
         ik2 = iBasi
         ik3 = iCmp(3)
         ik4 = iCmp(1)
        End If
        If (jShell(1).ge.jShell(4)) Then
         il1 = iBasi
         il2 = lBasl
         il3 = iCmp(1)
         il4 = iCmp(4)
        Else
         il1 = lBasl
         il2 = iBasi
         il3 = iCmp(4)
         il4 = iCmp(1)
        End If
        If (jShell(2).ge.jShell(3)) Then
         jk1 = jBasj
         jk2 = kBask
         jk3 = iCmp(2)
         jk4 = iCmp(3)
        Else
         jk1 = kBask
         jk2 = jBasj
         jk3 = iCmp(3)
         jk4 = iCmp(2)
        End If
        If (jShell(2).ge.jShell(4)) Then
         jl1 = jBasj
         jl2 = lBasl
         jl3 = iCmp(2)
         jl4 = iCmp(4)
        Else
         jl1 = lBasl
         jl2 = jBasj
         jl3 = iCmp(4)
         jl4 = iCmp(2)
        End If
*
*    Here we go
*
        nao=iBasi*jBasj*kBask*lBasl*
     &    iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4)
        Call CtlDns(iDCRR,iDCRS,iDCRT,jOp)
         end if
*
*             Add this contribution to the (Inactive) Fock Matrix
*
*                Out from  CD : jOp
*
        Do iCent=1,4
         Do iCar=1,3
            Call lCopy(8,[.false.],0,pert,1)
*
*           Too which irreps does this derivative contribute ?
*
            Do iIrrep=0,nIrrep-1
               If (indgrd(iCar,iCent,iirrep).ne.0)
     &          pert(iIrrep)=.true.
            End Do
*
           If (Index(iCar,iCent).gt.0) Then
            iGr=Index(iCar,iCent)-1
            ipFin=1+iGr*nijkl*nabcd
            If (.not.new_fock) Then
            Call MkFck(iAnga,iCmpa,iCmp,
     &                 Shijij,
     &                 iShll,iShell,
     &                 iBasi,jBasj,kBask,lBasl,
     &                 iAO,iAOst,nop,jop,
     &                 Dij1,mDij,nDij,ij1,ij2,ij3,ij4,
     &                 Dkl1,mDkl,nDkl,kl1,kl2,kl3,kl4,
     &                 Dik1,mDik,nDik,ik1,ik2,ik3,ik4,
     &                 Dil1,mDil,nDil,il1,il2,il3,il4,
     &                 Djk1,mDjk,nDjk,jk1,jk2,jk3,jk4,
     &                 Djl1,mDjl,nDjl,jl1,jl2,jl3,jl4,
     &                 Final(ipFin),
     &                 nAO,TwoHam,nTwo,
     &                 Scrtch1,nS1,Scrtch2,nS2,
     &                 iDCRR,iDCRS,iDCRT,
     &                 FckTmp,nFT,pert,iuvwx(iCent),
     &                 iCent,iCar,indgrd,ipdisp)
            If (nMethod.eq.RASSCF)
     &      Call MkFck(iAnga,iCmpa,iCmp,
     &                 Shijij,
     &                 iShll,iShell,
     &                 iBasi,jBasj,kBask,lBasl,
     &                 iAO,iAOst,nop,jop,
     &                 Dij2,mDij,nDij,ij1,ij2,ij3,ij4,
     &                 Dkl2,mDkl,nDkl,kl1,kl2,kl3,kl4,
     &                 Dik2,mDik,nDik,ik1,ik2,ik3,ik4,
     &                 Dil2,mDil,nDil,il1,il2,il3,il4,
     &                 Djk2,mDjk,nDjk,jk1,jk2,jk3,jk4,
     &                 Djl2,mDjl,nDjl,jl1,jl2,jl3,jl4,
     &                 Final(ipFin),
     &                 nAO,TwoHam,nTwo,
     &                 Scrtch1,nS1,Scrtch2,nS2,
     &                 iDCRR,iDCRS,iDCRT,
     &                 FckTmp,nFT,pert,iuvwx(iCent),
     &                 iCent,iCar,indgrd,ipdisp2)
        else
         ip=ipDisp(abs(indgrd(iCar,iCent,0)))
         Call FckAcc_NoSym(iAnga,
     &                  iCmpa(1),iCmpa(2),iCmpa(3),iCmpa(4), Shijij,
     &                  iShll, iShell, nijkl,
     &                  Final(ipFin),TwoHam(ip),dan,ndens,
     &                  iAO,iAOst,iBasi,jBasj,kBask,lBasl,ExFac)
         If (nMethod.eq.RASSCF) Then
         ip=ipDisp2(abs(indgrd(iCar,iCent,0)))
         Call FckAcc_NoSym(iAnga,
     &                  iCmpa(1),iCmpa(2),iCmpa(3),iCmpa(4), Shijij,
     &                  iShll, iShell, nijkl,
     &                  Final(ipFin),TwoHam(ip),din,nDens,
     &                  iAO,iAOst,iBasi,jBasj,kBask,lBasl,ExFac)
        end if
        end if
*
*
          Else If (Index(iCar,iCent).lt.0) Then
             call dcopy_(nabcd*nijkl,[Zero],0,Temp,1)
             Do iCnt=1,4
                iGr=Index(iCar,iCnt)
                If (iGr.gt.0)Then
                   ipFin=1+(iGr-1)*nijkl*nabcd
                   Do ii=1,nabcd*nijkl
                      Temp(ii)=Temp(ii)-Final(ipFin-1+ii)
                   End Do
                End If
             End Do
*
           if (.not.new_fock) Then
           Call MkFck(iAnga,iCmpa,iCmp,
     &                Shijij,
     &                iShll,iShell,
     &                iBasi,jBasj,kBask,lBasl,
     &                iAO,iAOst,nop,jop,
     &                Dij1,mDij,nDij,ij1,ij2,ij3,ij4,
     &                Dkl1,mDkl,nDkl,kl1,kl2,kl3,kl4,
     &                Dik1,mDik,nDik,ik1,ik2,ik3,ik4,
     &                Dil1,mDil,nDil,il1,il2,il3,il4,
     &                Djk1,mDjk,nDjk,jk1,jk2,jk3,jk4,
     &                Djl1,mDjl,nDjl,jl1,jl2,jl3,jl4,
     &                Temp,
     &                nAO,TwoHam,nTwo,
     &                Scrtch1,nS1,Scrtch2,nS2,
     &                iDCRR,iDCRS,iDCRT,
     &                FckTmp,nFT,pert,iuvwx(iCent),
     &                icent,iCar,indgrd,ipdisp)
           If (nMethod.eq.RASSCF)
     &     Call MkFck(iAnga,iCmpa,iCmp,
     &                Shijij,
     &                iShll,iShell,
     &                iBasi,jBasj,kBask,lBasl,
     &                iAO,iAOst,nop,jop,
     &                Dij2,mDij,nDij,ij1,ij2,ij3,ij4,
     &                Dkl2,mDkl,nDkl,kl1,kl2,kl3,kl4,
     &                Dik2,mDik,nDik,ik1,ik2,ik3,ik4,
     &                Dil2,mDil,nDil,il1,il2,il3,il4,
     &                Djk2,mDjk,nDjk,jk1,jk2,jk3,jk4,
     &                Djl2,mDjl,nDjl,jl1,jl2,jl3,jl4,
     &                Temp,
     &                nAO,TwoHam,nTwo,
     &                Scrtch1,nS1,Scrtch2,nS2,
     &                iDCRR,iDCRS,iDCRT,
     &                FckTmp,nFT,pert,iuvwx(iCent),
     &                icent,iCar,indgrd,ipdisp2)
*
         else
         ip=ipDisp(abs(indgrd(iCar,iCent,0)))
         Call FckAcc_NoSym(iAnga,
     &                  iCmpa(1),iCmpa(2),iCmpa(3),iCmpa(4), Shijij,
     &                  iShll, iShell, nijkl,
     &                  Temp,TwoHam(ip),dan,nDens,
     &                  iAO,iAOst,iBasi,jBasj,kBask,lBasl,ExFac)
            If (nMethod.eq.RASSCF) Then
         ip=ipDisp2(abs(indgrd(iCar,iCent,0)))
         Call FckAcc_NoSym(iAnga,
     &                  iCmpa(1),iCmpa(2),iCmpa(3),iCmpa(4), Shijij,
     &                  iShll, iShell, nijkl,
     &                  Temp,TwoHam(ip),din,nDens,
     &                  iAO,iAOst,iBasi,jBasj,kBask,lBasl,ExFac)
        end if
        end if
*
          End If
         End Do
        End Do
        Call Timing(dum1,Time,dum2,dum3)
        CPUStat(nFckAck)=CPUStat(nFckAck)+Time
      End If
*
      If (n8.and.nmethod.eq.RASSCF)
     &Call MakeMO(Final,Scrtch1,nTempTot,nFinal,
     &                  TwoHam,nTwo,
     &                  iCmp,iCmpa,
     &                  iBasi,jbasj,kbask,lbasl,
     &                  nGr,index,
     &                  moip,naco,nop,indgrd,
     &                  ishll,ishell,rmoin,nMOIN,
     &                  iuvwx,iao,iaost,Buffer,ianga,c)
*
*
      Call Timing(dum1,Time,dum2,dum3)
      CPUStat(nMOTrans)=CPUStat(nMOTrans)+Time
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(istb)
         Call Unused_integer(jstb)
         Call Unused_integer(kstb)
         Call Unused_integer(lstb)
         Call Unused_logical_array(IfG)
      End If
      End
