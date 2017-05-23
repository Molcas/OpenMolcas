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
      Subroutine EditStart
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"
#include "WrkSpc.fh"
#include "warnings.fh"

      Parameter (ThrdSpread=1.0d0)
      Dimension iC(3),iC2(3)
      Dimension Coord(MxCen*MxPut,3),Coo(MxCen,3),CooRef(MxCen,3)
      Character Filstart*6,FilSlut*6,Head*200
      Logical Exist,ValidOrNot

*
*-- Inquire if file exists, and if so open it.
*
      Write(FilStart,'(A5,i1.1)')'STFIL',NrStarti
      Call f_Inquire(FilStart,Exist)
      If(.not.Exist) then
        Write(6,*)
        Write(6,*)'The input startfile given in the EDITstartfile'
     &//' section was not found.'
        Call Quit(_RC_IO_ERROR_READ_)
      Endif
      iLu=73
      Call DaName(iLu,Filstart)

*
*-- Read header and coordinates.
*
      iDisk=0
      Call WrRdSim(iLu,2,iDisk,iTcSim,64,Etot,Ract,nPart,Gamold
     &            ,GaOld,Esub)
      iDisk=iTcSim(1)
      Do 5,l=1,3
        Call GetMem('Coordinates','Allo','Real',iC(l),nPart*nCent)
        Call dDaFile(iLu,2,Work(iC(l)),nPart*nCent,iDisk)
5     Continue
      Call DaClos(iLu)

*
*-- Now take different paths depending of what user have requested in
*   the input.
*

*
*-- If deleting solvent molecules.
*
      If(DelOrAdd(1)) then
*
*----Find the solvent molecules fartherst away from origo and
*    delete them.
*
        Do 10, i=1,nDel
          rMax=0.0d0
          indMax=0
          Do 20, j=1,nPart
            r=0.0d0
            Do 30, k=1,3
              r=r+Work(iC(k)+(j-1)*nCent)**2
30          Continue
            If(r.gt.rMax) then
              rMax=r
              indMax=j
            Endif
20        Continue
          Do 40, j=indMax,nPart-i
            Do 50, l=1,3
              Do 60, ll=1,nCent
                Work(iC(l)+(j-1)*nCent+ll-1)=Work(iC(l)+j*nCent+ll-1)
60            Continue
50          Continue
40        Continue
10      Continue
*
*----Print the new set of coordinates to a startfile.
*
        iLu=74
        Write(FilSlut,'(A5,i1.1)')'STFIL',NrStartu
        Call DaName(iLu,FilSlut)
        iDisk=0
        Call WrRdSim(iLu,1,iDisk,iTcSim,64,Etot,Ract,nPart-nDel,Gamold
     &              ,Gaold,Esub)
        iTcSim(1)=iDisk
        Do 70, l=1,3
          Call dDaFile(iLu,1,Work(iC(l)),(nPart-nDel)*nCent,iDisk)
          iTcSim(1+l)=iDisk
70      Continue
        iDisk=0
        Call WrRdSim(iLu,1,iDisk,iTcSim,64,Etot,Ract,nPart-nDel,Gamold
     &              ,Gaold,Esub)
        Call DaClos(iLu)
*
*----If user wants, print print print.
*
        If(iPrint.ge.10) then
          Do 1001, i=1,(nPart-nDel)*nCent
            Coord(i,1)=Work(iC(1)+i-1)
            Coord(i,2)=Work(iC(2)+i-1)
            Coord(i,3)=Work(iC(3)+i-1)
1001      Continue
          Write(Head,*)'Final coordinates'
          Call Cooout(Head,Coord,nPart-nDel,nCent)
        Endif
      Endif

*
*-- If adding solvent molecules.
*
      If(DelOrAdd(2)) then
        Do 301, i=1,nPart*nCent
          Coord(i,1)=Work(iC(1)+i-1)
          Coord(i,2)=Work(iC(2)+i-1)
          Coord(i,3)=Work(iC(3)+i-1)
301     Continue
        If(nAdd.ne.0) then
*---- Just an ugly trick for using nypart. It requires that the first
*     slot contains the solvent coordinates, so we, temporarily, put
*     them there.
          Do 3010, i=1,nCent
            Do 3011, j=1,3
              Coord(i,j)=Cordst(i,j)
3011        Continue
3010      Continue
*---- Introduce the new particles. nPart is redefined.
          Call NyPart(nAdd,nPart,Coord,rStart,nCent,iSeed)
*---- The ugly trick is reversed, and the first slot is retained.
          Do 3012, i=1,nCent
            Do 3013, j=1,3
              Coord(i,j)=Work(iC(j)+i-1)
3013        Continue
3012      Continue
        Endif
*
*----Then dump new coordinates on designated startfile.
*
        Do 3001, k=1,3
          Call GetMem('NewCoo','Allo','Real',iC2(k),nPart*nCent)
3001    Continue
        Do 302, i=1,nPart*nCent
          Work(iC2(1)+i-1)=Coord(i,1)
          Work(iC2(2)+i-1)=Coord(i,2)
          Work(iC2(3)+i-1)=Coord(i,3)
302     Continue
        iLu=74
        Write(FilSlut,'(A5,i1.1)')'STFIL',NrStartu
        Call DaName(iLu,FilSlut)
        iDisk=0
        Call WrRdSim(iLu,1,iDisk,iTcSim,64,Etot,RStart,nPart,Gamold
     &              ,Gaold,Esub)
        iTcSim(1)=iDisk
        Do 270, l=1,3
          Call dDaFile(iLu,1,Work(iC2(l)),nPart*nCent,iDisk)
          iTcSim(1+l)=iDisk
270     Continue
        iDisk=0
        Call WrRdSim(iLu,1,iDisk,iTcSim,64,Etot,RStart,nPart,Gamold
     &              ,Gaold,Esub)
        Call DaClos(iLu)
        If(iPrint.ge.10) then
          Write(Head,*)'Final coordinates'
          Call Cooout(Head,Coord,nPart,nCent)
        Endif
        Do 3002, k=1,3
          Call GetMem('NewCoo','Free','Real',iC2(k),nPart*nCent)
3002    Continue
      Endif

*
*-- If requested, substitute all particles that are not of valid water
*   geometry for, you guessed it, valid water molecules.
*
      If(DelOrAdd(3)) then
        nRemoved=0
        Do 451, iPart=1,nPart
          ind=nCent*(iPart-1)
          Do 452, iCent=1,nCent
            Coo(iCent,1)=Work(iC(1)+ind+iCent-1)
            Coo(iCent,2)=Work(iC(2)+ind+iCent-1)
            Coo(iCent,3)=Work(iC(3)+ind+iCent-1)
            CooRef(iCent,1)=Cordst(iCent,1)
            CooRef(iCent,2)=Cordst(iCent,2)
            CooRef(iCent,3)=Cordst(iCent,3)
452       Continue
          Call IsItValid(Coo,CooRef,ValidOrNot)
          If(.not.ValidOrNot) then
            dCMx=0.0d0
            dCMy=0.0d0
            dCMz=0.0d0
            Do 453, iCent=1,nCent
              dCMx=dCMx+Work(iC(1)+ind+iCent-1)
              dCMy=dCMy+Work(iC(2)+ind+iCent-1)
              dCMz=dCMz+Work(iC(3)+ind+iCent-1)
453         Continue
            dCMx=dCMx*(1.0d0/dble(nCent))
            dCMy=dCMy*(1.0d0/dble(nCent))
            dCMz=dCMz*(1.0d0/dble(nCent))
*------ Check if the points are spread out, otherwise just delete.
            dSpread=0.0d0
            Do 455, iCent=1,nCent
              dSpread=dSpread+(Work(iC(1)+ind+iCent-1)-dCMx)**2
              dSpread=dSpread+(Work(iC(2)+ind+iCent-1)-dCMy)**2
              dSpread=dSpread+(Work(iC(3)+ind+iCent-1)-dCMz)**2
455         Continue
            dSpread=dSpread*(1.0d0/dble(nCent))
            If(dSpread.lt.ThrdSpread) then
              nRemoved=nRemoved+1
              Do 456, jP=iPart,nPart-1
                jnd1=(jP-1)*nCent
                jnd2=(jP)*nCent
                Do 457, iCent=1,nCent
                  Work(iC(1)+jnd1+iCent-1)=Work(iC(1)+jnd2+iCent-1)
                  Work(iC(2)+jnd1+iCent-1)=Work(iC(2)+jnd2+iCent-1)
                  Work(iC(3)+jnd1+iCent-1)=Work(iC(3)+jnd2+iCent-1)
457             Continue
456           Continue
            Else
              Do 454, iCent=1,nCent
                Work(iC(1)+ind+iCent-1)=dCMx+CooRef(iCent,1)
                Work(iC(2)+ind+iCent-1)=dCMy+CooRef(iCent,2)
                Work(iC(3)+ind+iCent-1)=dCMz+CooRef(iCent,3)
454           Continue
            Endif
          Endif
451     Continue
        nPart=nPart-nRemoved
*
*----Print the new set of coordinates to a startfile.
*
        iLu=74
        Write(FilSlut,'(A5,i1.1)')'STFIL',NrStartu
        Call DaName(iLu,FilSlut)
        iDisk=0
        Call WrRdSim(iLu,1,iDisk,iTcSim,64,Etot,Ract,nPart,Gamold
     &              ,Gaold,Esub)
        iTcSim(1)=iDisk
        Do 71, l=1,3
          Call dDaFile(iLu,1,Work(iC(l)),nPart*nCent,iDisk)
          iTcSim(1+l)=iDisk
71      Continue
        iDisk=0
        Call WrRdSim(iLu,1,iDisk,iTcSim,64,Etot,Ract,nPart,Gamold
     &              ,Gaold,Esub)
        Call DaClos(iLu)
*
*----If user wants, print print print.
*
        If(iPrint.ge.10) then
          Do 1002, i=1,nPart*nCent
            Coord(i,1)=Work(iC(1)+i-1)
            Coord(i,2)=Work(iC(2)+i-1)
            Coord(i,3)=Work(iC(3)+i-1)
1002      Continue
          Write(Head,*)'Final coordinates'
          Call Cooout(Head,Coord,nPart,nCent)
        Endif
      Endif

*
*-- If the user want to, print the coordinates in some format suitable
*   for graphical representation.
*
      If(DelOrAdd(4)) then
        If(cDumpForm(1:4).eq.'MOLD') then
          Do 444, iCent=1,nCent
            CooRef(iCent,1)=Cordst(iCent,1)
            CooRef(iCent,2)=Cordst(iCent,2)
            CooRef(iCent,3)=Cordst(iCent,3)
444       Continue
          Call MoldenDump(iC,CooRef,nPart,nAtom,nCent)
        Endif
      Endif

*
*-- Deallocate.
*
      Do 501,l=1,3
        Call GetMem('Coordinates','Free','Real',iC(l),nPart*nCent)
501   Continue

*
*-- This routine ends now!
*
      Return
      End
