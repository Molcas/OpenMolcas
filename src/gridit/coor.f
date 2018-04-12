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

      Subroutine MyCoor_nosupport(iRun,isAuto,
     &      Ox,Oy,Oz, Rx,Ry,Rz, iGx,iGy,iGz, iMagic,iCustOrig)
************************************************************************
*                                                                      *
*   Read Coordinates and calculate a cub for grid                      *
*   isAuto=1 - real job, else only print                               *
*   Origin(3) fix cub in space                                         *
*   Rx,Ry,Rz - size of cub                                             *
*   iGx,iGy,iGz - net                                                  *
*   iMagic = magic guess for net                                       *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
      Include 'WrkSpc.fh'
      Include 'grid.nosupport.fh'
      Include 'real.fh'
      Include 'periodic_table.fh'
c      Include 'wdw.fh'
      Dimension iOper(8)
      Dimension RotVec(3)
      Character*(LENIN) Byte4
      Character*4 tt
      Character*80 Line
*----------------------------------------------------------------------*
*     Prologue                                                         *
*----------------------------------------------------------------------*
#ifdef _DEBUG_
      Call qEnter('Coor')
#endif
*----------------------------------------------------------------------*
*     Read no.of symm. species                                         *
*----------------------------------------------------------------------*
      Call get_iScalar('nSym',nSym)
*----------------------------------------------------------------------*
*     Read symm. oper per symm. species                                *
*----------------------------------------------------------------------*
      Call Get_iArray('Symmetry operations',iOper,nSym)
*----------------------------------------------------------------------*
*     Read no. of unique atoms in the system                           *
*----------------------------------------------------------------------*
      Call Get_iScalar('Unique atoms',nAtoms)
*----------------------------------------------------------------------*
*     Read atom labels                                                 *
*----------------------------------------------------------------------*
      lw2=1
      Call Get_cArray('Unique Atom Names',AtomLbl,LENIN*nAtoms)
*----------------------------------------------------------------------*
*     Read coordinates of atoms                                        *
*----------------------------------------------------------------------*
      Call GetMem('Coor','ALLO','REAL',ipCoor,3*nSym*nAtoms)
      Call Get_dArray('Unique Coordinates',Work(ipCoor),3*nAtoms)
*----------------------------------------------------------------------*
*     Read nuclear repulsion energy                                    *
*----------------------------------------------------------------------*
*     Call Get_PotNuc(PotNuc)
c      Call Get_dScalar('PotNuc',PotNuc)
*----------------------------------------------------------------------*
*     Apply the symmetry operations                                    *
*----------------------------------------------------------------------*
      nOper=0
      If ( nSym.eq.2 ) nOper=1
      If ( nSym.eq.4 ) nOper=2
      If ( nSym.eq.8 ) nOper=3
      nCenter=nAtoms
      Do i=1,nOper
        jOper=i+1
        If ( i.eq.3 ) jOper=5
        RotVec(1)=One
        If ( IAND(iOper(jOper),1).eq.1 ) RotVec(1)=-One
        RotVec(2)=One
        If ( IAND(iOper(jOper),2).eq.2 ) RotVec(2)=-One
        RotVec(3)=One
        If ( IAND(iOper(jOper),4).eq.4 ) RotVec(3)=-One
        newAt=0
        mCenter=nCenter
        Do iAt=0,mCenter-1
          Xold=Work(ipCoor+iAt*3+0)
          Yold=Work(ipCoor+iAt*3+1)
          Zold=Work(ipCoor+iAt*3+2)
          Byte4=AtomLbl(lw2+iAt)
          Xnew=RotVec(1)*Xold
          Ynew=RotVec(2)*Yold
          Znew=RotVec(3)*Zold
          Do jAt=0,nCenter-1
             If (Byte4.eq.AtomLbl(Lw2+jAt)) Then
                Xold2=Work(ipCoor+jAt*3+0)
                Yold2=Work(ipCoor+jAt*3+1)
                Zold2=Work(ipCoor+jAt*3+2)

          If ( Xnew.eq.Xold2.and.Ynew.eq.Yold2.and.Znew.eq.Zold2)
     &                   goto 999
             Endif
            Enddo
            Work(ipCoor+nCenter*3+0)=Xnew
            Work(ipCoor+nCenter*3+1)=Ynew
            Work(ipCoor+nCenter*3+2)=Znew
            AtomLbl(lw2+nCenter)=Byte4
            nCenter=nCenter+1
            newAt=newAt+1
 999    Continue
        End Do
c        nCenter=nCenter+newAt
      End Do
* Check : are any strange names?
      NoOrig=0
      Do iAt=0,nCenter-1
       write(line,'(a)') AtomLbl(lw2+iAt)
       if(index( line,'Ori').ne.0) NoOrig=NoOrig+1
      enddo
*
      if(isBinary.eq.0) then
      if(iGauss.eq.0.and.isLine.eq.0) then
      write(LuVal,'(a,i8)') 'Natom= ',nCenter-NoOrig
      if(isUHF.eq.1) write(LuVal_ab,'(a,i8)') 'Natom= ',nCenter-NoOrig
      Do iAt=0,nCenter-1
       write(line,'(a)') AtomLbl(lw2+iAt)

       if(index( line,'Ori').eq.0) then
        Write(LuVal,'(A,2X,3F10.5)') AtomLbl(lw2+iAt),
     &  Work(ipCoor+3*iAt),Work(ipCoor+3*iAt+1),Work(ipCoor+3*iAt+2)
        if(isUHF.eq.1) Write(LuVal_ab,'(A,2X,3F10.5)') AtomLbl(lw2+iAt),
     &  Work(ipCoor+3*iAt),Work(ipCoor+3*iAt+1),Work(ipCoor+3*iAt+2)
       endif
      enddo
      else
       nGatoms=nCenter-NoOrig
       call molcas_open(37,'coord.gtmp')
c       open(37,file='coord.gtmp')
       do iip=1,103
       call upcase(PTab(iip))
       enddo
      Do iAt=0,nCenter-1
       write(Byte4,'(a)') AtomLbl(lw2+iAt)
c Well, let's draw He in a worst case
       ip=2
       call upcase(Byte4)
        do ii=1,4
        if(index ('0123456789',Byte4(ii:ii)).ne.0) Byte4(ii:ii)=' '
        enddo
        if(Byte4(2:2).eq.' ') then
        Byte4(2:2)=Byte4(1:1)
        Byte4(1:1)=' '
        endif
       do iip=1,103
       if(Byte4.eq.PTab(iip)) then
         ip=iip
       endif
       enddo
       if(index( Byte4,'Ori').eq.0)
     & write(37,'(I5,4F12.6)')
     &  ip,0.0,
     &  Work(ipCoor+3*iAt+2),Work(ipCoor+3*iAt+1),Work(ipCoor+3*iAt)
      enddo
      close(37)
      endif
      endif
      if(isBinary.eq.1) then
      write(Line,'(a,i8)') 'Natom= ',nCenter-NoOrig
      write(LuVal) Line(1:15)
      if(isUHF.eq.1) write(LuVal_ab) Line(1:15)
      Do iAt=0,nCenter-1
       write(line,'(a)') AtomLbl(lw2+iAt)
       if(index( line,'Ori').eq.0)
     & Write(line,'(A,2X,3F10.5)')
     &  AtomLbl(lw2+iAt),
     &  Work(ipCoor+3*iAt),Work(ipCoor+3*iAt+1),Work(ipCoor+3*iAt+2)
      write(LuVal) Line(1:40)
      if(isUHF.eq.1) write(LuVal_ab) Line(1:40)
      enddo
      endif
      if(iCustOrig.eq.1) goto 500

      if(isAuto.eq.1) Then
*----------------------------------------------------------------------*
*     Find Cub parameters                                              *
*                 Ox->RxMin, Rx->RxMax
*----------------------------------------------------------------------*
      Ox=9999.9D0
      Oy=9999.9D0
      Oz=9999.9D0
      Rx=-9999.9D0
      Ry=-9999.9D0
      Rz=-9999.9D0

      Do iAt=0,nCenter-1
       rrx=Work(ipCoor+3*iAt)
       rry=Work(ipCoor+3*iAt+1)
       rrz=Work(ipCoor+3*iAt+2)
       if(rrx.lt.Ox) Ox=rrx
       if(rrx.gt.Rx) Rx=rrx
       if(rry.lt.Oy) Oy=rry
       if(rry.gt.Ry) Ry=rry
       if(rrz.lt.Oz) Oz=rrz
       if(rrz.gt.Rz) Rz=rrz
      End Do
       Rx=Rx-Ox
       Ry=Ry-Oy
       Rz=Rz-Oz
*
* and now, expand this cub to place all atoms inside
*
      Ox=dble(int(Ox-TheGap))
      Oy=dble(int(Oy-TheGap))
      Oz=dble(int(Oz-TheGap))
      Rx=dble(int(Rx+2.0d0*TheGap))
      Ry=dble(int(Ry+2.0d0*TheGap))
      Rz=dble(int(Rz+2.0d0*TheGap))
* finish of iAuto.
      endif
*----------------------------------------------------------------------*
*     Calculate corrected coords
*----------------------------------------------------------------------*
*
* make a stupid Patch: Cerius2 works well only with even nets!
*
       Rx=dble(int(Rx)/2*2)
       Ry=dble(int(Ry)/2*2)
       Rz=dble(int(Rz)/2*2)

       if(iMagic.gt.0) then
           iGx=int(abs(Rx))*iMagic
           iGy=int(abs(Ry))*iMagic
           iGz=int(abs(Rz))*iMagic
        endif
*
* make a stupid Patch: Cerius2 works well only with even nets!
*
      if(iGauss.eq.0) then
      iGx=(iGx+1)/2*2
      iGy=(iGy+1)/2*2
      iGz=(iGz+1)/2*2
      else
      Rx=max(Rx,Ry,Rz)
      Ry=Rx
      Rz=Rx
      iGx=max(iGx,iGy,iGz)
      iGy=iGx
      iGz=iGx
      Ox=min(Ox,Oy,Oz)
      Oy=Ox
      Oz=Ox
      endif
      mynCenter=nCenter

*----------------------------------------------------------------------*
*     Print coordinates of the system                                  *
*----------------------------------------------------------------------*
c      call bXML('Coord')
c      call iXML('nCoord',nCenter)
c      Do iAt=0,nCenter-1
c        Call cXML('Atom',AtomLbl(lw2+iAt))
c        call daXML('Atom coord',Work(ipCoor+3*iAt),3)
c      End Do
c      call eXML('Coord')
      isWDW=1
      call getmem('WDW','ALLO','REAL',ipWdW,nCenter)
      Do iAt=0,nCenter-1
        Byte4=AtomLbl(lw2+iAt)
       call upcase(Byte4)
        do ii=1,4
        if(index ('0123456789',Byte4(ii:ii)).ne.0) Byte4(ii:ii)=' '
        enddo
        if(Byte4(2:2).eq.' ') then
        Byte4(2:2)=Byte4(1:1)
        Byte4(1:1)=' '
        endif
        do j=1,103
          tt=PTab(j)
          call upcase(tt)
          if(Byte4 .eq. tt) Work(ipWdW+iAt)=0
        enddo
c      print *,'here',Byte4,Work(ipWdW+iAt)
      enddo


500   if(iRun.eq.1.and.levelprint.ge.3) then
      Write(6,*)
      Write(6,'(6X,A)')'Cartesian coordinates:'
      Write(6,'(6X,A)')'-----------------------------------------'
      Write(6,'(6X,A)')'No.  Label     X         Y         Z     '
      Write(6,'(6X,A)')'-----------------------------------------'
      Do iAt=0,nCenter-1
        Write(6,'(4X,I4,3X,A,2X,3F10.5)')
     &  iAt+1,AtomLbl(lw2+iAt),
     &  Work(ipCoor+3*iAt),Work(ipCoor+3*iAt+1),Work(ipCoor+3*iAt+2)
      End Do
      Write(6,'(6X,A)')'-----------------------------------------'
      Write(6,'(6X,A,3F12.6)')'Grid Origin      = ',Ox,Oy,Oz
      Write(6,'(6X,A,3F12.6)')'Grid Axis Length = ',Rx,Ry,Rz
      Write(6,'(6X,A)')'-----------------------------------------'
      Write(6,*)
      Write(6,*)
      endif
c       if (isAtom.eq.0)
c     &   Call GetMem('Coor','FREE','REAL',ipCoor,3*nSym*nAtoms)
*
*----------------------------------------------------------------------*
*     Normal exit                                                      *
*----------------------------------------------------------------------*
#ifdef _DEBUG_
      Call qExit('Coor')
#endif
      Return
      End
*----------------------------------------------------------------------*
