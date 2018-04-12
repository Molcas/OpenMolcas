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

      SubRoutine DrvMO_nosupport(iRun,INPORB)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from: RhoCmp                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              MOEval                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             Valera Veryazov, Dept. Theoretical Chemistry             *
*                                                                      *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
      Include 'real.fh'
      Include 'WrkSpc.fh'
#include <SysDef.fh>
      Include 'grid.nosupport.fh'
      Logical Debug, long_prt
      Logical is_error
      Character str*80, Crypt*7
      Character INPORB*(*)
c      Character bb
c      character namepack*3
      Character*80 myTitle
      Character*64 status
      Real*8 pp(3)
      Integer nTypes(7)
c      Character *10 Pack, oPack
c      Dimension xLimits(4)
c      Integer iYDelta(3)
c      Real tmp(3)
      Data Crypt/'fi123sd'/
*
*---- Set size of batches
*
* vv: to print in Gaussian cube format it's better
*     if nInc/6*6=nInc
*
* ag: symbolic constant -- to allow to be array dimension
* ag: FIXME: temporary hack
*
      parameter (nIncPack=18*1024)

      Character cMoBlock(nIncPack*2)

* ag: actual nInc is variable
      nInc=nIncPack

*                                                                      *
************************************************************************
*
*...  Prologue
#ifdef _DEBUG_
      Debug=.true.
#else
      Debug=.false.
#endif
      Call qEnter('DrvMO')
      isEner=1
      iRout = 2
c      iPack=1
       ipCutOff=ip_iDummy

      iPrint = nPrint(iRout)
      If (Debug) iPrint=99
      dNorm=0
      ddNorm=0
      if(iRun.eq.1.and.levelprint.ge.3) then
       Write (6,*)
       Write (6,'(A,8I5)') 'Irreps  : ',(i,i=1,nIrrep)
       Write (6,'(A,8I5)') 'Basis   : ',(nBas(i),i=0,nIrrep-1)
       Write (6,'(A,3I5)') 'Grid Net: ',(iGridNpt(i),i=1,3)
       Write (6,*)
      endif
*
*...  Compute the size of the densities
*
      nCMO=0
      nMOs=0
      Do iIrrep = 0, nIrrep - 1
         nCMO  = nCMO + nBas(iIrrep)**2
         nMOs  = nMOs + nBas(iIrrep)
      End Do
*
      Call GetMem('CMO','Allo','Real',ipCMO,nCMO)
      Call GetMem('Ener','Allo','Real',ipE,nMOs)
      Call GetMem('Occu','Allo','Real',ipOcc,nMOs)
      Call GetMem('Occ2','Allo','Real',ipOoo,nMOs)
      if(isVirt.eq.1) then
        Call GetMem('ddNo','Allo','Real',ipdd,nMOs*nMOs)
        Do i=1,nMOs*nMOs
         Work(ipdd+i-1)=0.0d0
        Enddo
      endif
      Call GetMem('ITyp','Allo','INTE',ipType,nMOs)
      Call GetMem('Vol','Allo','Real',ipVol,nMOs)
      Call GetMem('Sort','Allo','INTE',ipSort,nMOs)
      Call GetMem('Nzer','Allo','INTE',ipNZ,nMOs*2)
      Call GetMem('NRef','Allo','INTE',ipGRef,nMOs)
      Call GetMem('DoIt','Allo','INTE',ipDoIt,nMOs)
      Call GetMem('Pab','Allo','Real',ipPab,nMOs)
        Do i=0, nMOs-1
         iWork(ipGRef+i)=-1
        enddo
*
*  Read information from INPORB file
*
       LuOrb=isFreeUnit(46)
c      Call Chk_Rdvec('INPORB',LuOrb,isUHF)
c       print *,'inporb=',inporb
      if(isUHF.eq.0) then
      Call RdVec(INPORB,LuOrb,'COE',nIrrep,NBAS,NBAS,
     &  Work(ipCMO),Work(ipOcc),Work(ipE),iDummy,
     &  myTitle,0,iErr)
c construct Pab
       if(NoOrb.eq.1) then
c        print *,'nCMO,nMOs', nCMO,nMOs
        Call makePab_nosupport(Work(ipCMO),Work(ipOcc),
     *  Work(ipPab),nMOs,nMOs,
     &  nIrrep, nBas    )
c       print *,'Pab=', (Work(ipPab+i),i=0,nMOs-1)
       endif
        if(iErr.eq.1) then
          do j=0, nMOs-1
          Work(ipE+j)=Zero
          enddo
        endif
      Call RdVec(INPORB,LuOrb,'I',nIrrep,NBAS,NBAS,
     &  Dummy,Dummy,Dummy,iWork(ipType),
     &  myTitle,0,iErr)
        if(iErr.eq.1) then
          do j=0, nMOs-1
          iWork(ipType+j)=Zero
          enddo
        endif
      else  ! UHF case
c allocate memory for extra arrays.
      Call GetMem('CMO_ab','Allo','Real',ipCMO_ab,nCMO)
      Call GetMem('Ener_ab','Allo','Real',ipE_ab,nMOs)
      Call GetMem('Occu_ab','Allo','Real',ipOcc_ab,nMOs)
      Call GetMem('Sort_ab','Allo','INTE',ipSort_ab,nMOs)
      Call GetMem('NRef_ab','Allo','INTE',ipGRef_ab,nMOs)
      Call GetMem('DoIt_ab','Allo','INTE',ipDoIt_ab,nMOs)

      Call RdVec_(INPORB,LuOrb,'COE',1,nIrrep,NBAS,NBAS,
     &  Work(ipCMO),Work(ipCMO_ab),Work(ipOcc),Work(ipOcc_ab),
     &  Work(ipE),Work(ipE_ab),iDummy,
     &  myTitle,0,iErr,iWFtype)
c it can be only after SCF, so we do not need TypeIndex info
          do j=0, nMOs-1
          iWork(ipType+j)=Zero
          enddo

      endif
      do j=1,7
      nTypes(j)=0
      enddo
      do j=0, nMOs-1
      jj=iWork(ipType+j)
      if(jj.gt.0) nTypes(jj)=nTypes(jj)+1
      enddo
*
*  Calculate net.
*
       if(iRun.eq.1) then
        Write(6,'(A)') '   Input vectors read from INPORB'
        Write(6,'(A,A)') '   Orbital file label: ',myTitle
       endif
      do j=1,3
       iGridNpt(j)=iGridNpt(j)+1
      enddo
      if(iGauss.eq.1) then
        nInc=iGridNpt(1)**3
        if(nInc.gt.30000) then
       call WarningMessage(1,'The grid is too dense for CUBE option;'//
     & 'Trying to reset it!;;'//
     & 'If you need a better quality for the picture '//
     & 'DO NOT use CUBE keyword !')
        write(6,*)
        iGridNpt(1)=25
        iGridNpt(2)=25
        iGridNpt(3)=25
        nInc=iGridNpt(1)**3

        endif
      endif

      nCoor=iGridNpt(1)*iGridNpt(2)*iGridNpt(3)
      iiCoord=nCoor
*
*---- if isCutOff is in used - recalculate nCoor
*
       if(isCutOff.eq.1) then
        if(isUserGrid.eq.1) then
            write(6,*) 'Not implemented'
            call Quit_OnUserError
        endif
c       print *,'Natom=',nAtoms*nSym
c         do ii=1,nAtoms*nSym
c          print *,'coo',Work(ipCoor+ii*3-3), Work(ipCoor+ii*3-3+1),
c     &   Work(ipCoor+ii*3-3+2)
c         enddo
       call GetMem('CUTFL','ALLO','INTE',ipCutOff,nCoor)
       ie1=max(iGridNpt(1)-1,1)
       ie2=max(iGridNpt(2)-1,1)
       ie3=max(iGridNpt(3)-1,1)
c       print *,'vv', ie1, ie2, ie3
       iiCoord=0
       iiiCoord=0
       do i1=0,ie1
       do i2=0,ie2
       do i3=0,ie3
         pp(1)=GridOrigin(1)+
     &     GridAxis1(1)*i1/ie1+
     &     GridAxis2(1)*i2/ie2+
     &     GridAxis3(1)*i3/ie3
         pp(2)=GridOrigin(2)+
     &     GridAxis1(2)*i1/ie1+
     &     GridAxis2(2)*i2/ie2+
     &     GridAxis3(2)*i3/ie3
         pp(3)=GridOrigin(3)+
     &     GridAxis1(3)*i1/ie1+
     &     GridAxis2(3)*i2/ie2+
     &     GridAxis3(3)*i3/ie3
c       write(6,'(3f8.4)') pp
         ishow=0
         do ii=1,nAtoms*nSym
           iaia=0
           if(abs(Work(ipCoor+ii*3-3)  -pp(1)).lt.CutOff) iaia=iaia+1
           if(abs(Work(ipCoor+ii*3-3+1)-pp(2)).lt.CutOff) iaia=iaia+1
           if(abs(Work(ipCoor+ii*3-3+2)-pp(3)).lt.CutOff) iaia=iaia+1
c        write(6,'(7f8.4)') Work(ipCoor+ii*3-3),pp(1),
c     & Work(ipCoor+ii*3-3+1),pp(2),
c     &  Work(ipCoor+ii*3-3+2),pp(3),CutOff
            if(iaia.eq.3) ishow=1
          enddo
         if(ishow.eq.1) then
            iiCoord=iiCoord+1
            iWork(ipCutOff+iiiCoord)=1
          else
            iWork(ipCutOff+iiiCoord)=0
          endif
          iiiCoord=iiiCoord+1
       enddo
       enddo
       enddo
       write(6,*) nCoor-iiCoord, ' points are eliminated'
c       print *,'old=',nCoor,' New=', iiCoord
c       nCoor=iiCoord
       endif

      if(isTheOne.eq.1) nCoor=int(OneCoor(7)+0.3)
      if(isAtom.eq.1) nCoor=nAtoms
      if(isUserGrid.eq.1) nCoor=nGridPoints
      if(iRun.eq.1.and.levelprint.ge.1) then
      Write (6,*) ' Number of grid points in file:  ',nCoor
c      call iXML('nPoints',nCoor)
      endif

*************************************************************************
* And now we had to choose orbitals to draw.
*
* Sometime we had to make an automatic guess....
*
*************************************************************************

       Call PickOrb_nosupport(ipNz,ipSort,ipGref,ipSort_ab,
     &  ipGref_ab,ipVol,ipE,ipOcc,ipE_ab,ipOcc_ab,
     &  nShowMOs,nShowMOs_ab,isener,nMOs,myTitle,ipType)
*                                                                      *
*---- Start run over sets of grid points
*
555     continue
c       if(isBinary.eq.1) then
       nShowMOs=nShowMOs+isDensity+isSphere+isColor
c       else
c       nShowMOs=1
c       endif
      if(isUHF.eq.1) nShowMOs_ab=nShowMOs_ab+isDensity+isSphere+isColor
      nShowMOs2=nShowMOs+nShowMOs_ab
      if(iRun.eq.1.and.levelprint.ge.1) then
       Write (6,*)
       Write (6,*) ' Total number of MOs               :',nMOs
c       call iXML('nMOs',nMOs)
       Write (6,*) ' Number MOs for grid               :',nShowMOs2
       Write (6,*) ' Batches processed in increments of:',nInc
       Write (6,*)
      endif
      iPrintCount=0
      Call GetMem('MOValue','Allo','Real',ipMO,nInc*nMOs)
      Call GetMem('DOValue','Allo','Real',ipOut,nInc)
      if (imoPack .ne. 0) then
        Call GetMem('PackedBlock','Allo','INTE',ipPBlock,nInc)
      else
        Call GetMem('PackedBlock','Allo','INTE',ipPBlock,1)
        iWork(ipPBlock)=0
      endif
*
*.... Allocate memory for the some grid points
*
      Call GetMem('Coor','Allo','Real',ipC,nInc*3)
      iSphrDist =ip_Dummy
      iSphrColor=ip_Dummy

      if(isSphere.eq.1) then
      Call GetMem('SpDi','Allo','Real',iSphrDist,nInc)
      endif
      if(isColor.eq.1) then
      Call GetMem('SpCo','Allo','Real',iSphrColor,nInc)
      endif

*
*   check grids to calculate
*

      Do i=0, nMOs-1
        iWork(ipDoIt+i)=0
        if(isUHF.eq.1) iWork(ipDoIt_ab+i)=0
      enddo
c
      if(NoOrb.eq.1) goto 550
      Do i=1,nShowMOs-isDensity-isSphere-isColor
        iWork(ipDoIt+iWork(ipGRef+i-1)-1)=1
      enddo
      if(isUHF.eq.1) then
      Do i=1,nShowMOs_ab-isDensity-isSphere-isColor
        iWork(ipDoIt_ab+iWork(ipGRef_ab+i-1)-1)=1
      enddo
      endif
      ifpartial=1-isTotal
      if(isTotal.eq.1) then
        Do i=0, nMOs-1
          if(abs(Work(ipOcc+i)).gt.Zero) then
                     iWork(ipDoIt+i)=1
          endif
          if(isUHF.eq.1) then
          if(abs(Work(ipOcc_ab+i)).gt.Zero) then
                     iWork(ipDoIt_ab+i)=1
          endif
          endif
        enddo
      else
        Do i=0, nMOs-1
          if(Work(ipOcc+i).gt.Zero.and.iWork(ipDoIt+i).ne.1) then
                     ifpartial=1
          endif
          if(isUHF.eq.1) then
          if(Work(ipOcc_ab+i).gt.Zero.and.iWork(ipDoIt_ab+i).ne.1) then
                     ifpartial=1
          endif
          endif
        enddo
      endif
550   Continue
*
*   If plotting VB orbitals : count active orbitals and make sure
*   all are included in DOIT :
*
* nothing to do with UHF

      if(isVB.eq.1)then
        nActOrb=0
        dActEl=0d0
        Do i=0, nMOs-1
          if(Work(ipOcc+i).gt.0d0.and.Work(ipOcc+i).lt.2d0) then
            nActOrb=nActOrb+1
            dActEl=dActEl+Work(ipOcc+i)
            iWork(ipDoIt+i)=1
          endif
        enddo
        nActEl=nint(dActEl)
        vbOcc=DBLE(nActEl)/DBLE(nActOrb)
      endif
*
*   Print a header.
*
c
c AG: attempt to generalize binary & ascii
c
      if(iGauss.eq.1) then
      if(isUHF.eq.1) then
        write(6,*) 'Cube files and UHF is not implemented'
        call Quit_OnUserError()
      endif
        call PrintCubeHeader(nShowMOs)
        goto 3838
      endif
      call PrintHeader_nosupport(nMOs,nShowMOs,nShowMOs_ab,nCoor,nInc,
     & iiCoord,nTypes )
      Call PrintTitles_nosupport(LuVal,nShowMOs,isDensity,nMOs,
     &  iWork(ipGRef), isEner, isVB, Work(ipOcc), iWork(ipType), Crypt,
     &  iWork(ipNZ), Work(ipE), VBocc, ifpartial,isLine,isSphere,
     &  isColor)
      if(isUHF.eq.1) then
      Call PrintTitles_nosupport(LuVal_ab,nShowMOs_ab,isDensity,nMOs,
     &  iWork(ipGRef_ab), isEner, isVB, Work(ipOcc_ab),
     &  iWork(ipType), Crypt,
     &  iWork(ipNZ), Work(ipE_ab), VBocc, ifpartial,isLine,isSphere,
     &  isColor)
      endif
*                                                                      *
3838    continue

      If (isMULL.eq.1) Then
         long_prt=isLONGPRT.eq.1
         Call charge_grid_it(nSym,nBas,Work(ipCMO),nCMO,Work(ipOcc),
     &                       iWork(ipDoIt),long_prt)
         If (isUHF.eq.1) Then
            Write (6,*)
            Write (6,'(A)')       '      ----------------------'
            Write (6,'(A)')       '       BETA spinorbitals    '
            Write (6,'(A)')       '      ----------------------'
            Write (6,*)
            Call charge_grid_it(nSym,nBas,Work(ipCMO_ab),nCMO,
     &                          Work(ipOcc_ab),
     &                          iWork(ipDoIt_ab),long_prt)
         EndIf
      EndIf
*                                                                      *
      If(iRun.eq.1.and.levelprint.ge.1) then
       Call CollapseOutput(1, ' list of grids ')
      endif
*                                                                      *
************************************************************************
*                                                                      *
*....... Get VB transformation matrix
*
      if(isVB.eq.1)then
        Call GetMem('VBtmp','Allo','Real',ipVBtmp,nActOrb*nActOrb)
        call getvb2mo2_cvb(work(ipVBtmp),nActOrb,nActEl)
        Call GetMem('VBmat','Allo','Real',ipVBmat,nMOs*nActOrb)
        call fzero(work(ipVBmat),nMOs*nActOrb)
        iActOrb=0
        do iMOs=1,nMOs
          if(Work(ipOcc-1+iMOs).gt.0d0
     >       .and.
     >       Work(ipOcc-1+iMOs).lt.2d0) then
            iActOrb=iActOrb+1
            do jAct=1,nActOrb
              Work(ipVBmat-1+iMOs+(jAct-1)*nMOs)=
     *         Work(ipVBtmp-1+iActOrb+(jAct-1)*nActOrb)
            enddo
          endif
        enddo
        Call GetMem('VBtmp','Free','Real',ipVBtmp,nActOrb*nActOrb)
      else
        Call GetMem('VBmat','Allo','Real',ipVBmat,1)
        Work(ipVBmat)=0.0d0
      endif
*
*---- Loop over grid points in batches of nInc
*
      ive3=max(iGridNpt(3)-1,1)
      ive2=max(iGridNpt(2)-1,1)
      ive1=max(iGridNpt(1)-1,1)
       iv3=0
       iv2=0
       iv1=0
c      if(isAtom.eq.1) goto 6000

      iiiCoord=0

c      if(isCutOff.eq.1) nCoor=iiCoord

      iShiftCut=0
      Do iSec = 1, nCoor, nInc
      mCoor=Min(nInc,nCoor-iSec+1)
      write(status,'(a,i8,a,i8)') ' batch ',iSec,' out of ',nCoor/nInc
      call StatusLine('grid_it: ',status)
* Generate next portion of points..
      if (isTheOne.eq.0) then
c general case: we have a CUBIC box.
       ipPO=0
667    continue
       iiiCoord=iiiCoord+1
        gv3=1.0d+0*iv3/ive3
        gv2=1.0d+0*iv2/ive2
        gv1=1.0d+0*iv1/ive1
        if(isUserGrid.eq.0) then
          if(isCutOff.eq.0) then
            Work(ipC+ipPO*3)=GridOrigin(1)+
     *        GridAxis1(1)*gv1+GridAxis2(1)*gv2+GridAxis3(1)*gv3
            Work(ipC+ipPO*3+1)=GridOrigin(2)+
     *        GridAxis1(2)*gv1+GridAxis2(2)*gv2+GridAxis3(2)*gv3
            Work(ipC+ipPO*3+2)=GridOrigin(3)+
     *        GridAxis1(3)*gv1+GridAxis2(3)*gv2+GridAxis3(3)*gv3
          else
c using ipCutOff
            Work(ipC+ipPO*3)=40
            Work(ipC+ipPO*3+1)=40
            Work(ipC+ipPO*3+2)=40
             if(iWork(ipCutOff+iiiCoord-1).eq.1) then
              Work(ipC+ipPO*3)=GridOrigin(1)+
     *           GridAxis1(1)*gv1+GridAxis2(1)*gv2+GridAxis3(1)*gv3
              Work(ipC+ipPO*3+1)=GridOrigin(2)+
     *           GridAxis1(2)*gv1+GridAxis2(2)*gv2+GridAxis3(2)*gv3
              Work(ipC+ipPO*3+2)=GridOrigin(3)+
     *           GridAxis1(3)*gv1+GridAxis2(3)*gv2+GridAxis3(3)*gv3
              endif
          endif
        else
          Work(ipC+ipPO*3+1-1)=Work(ipGrid+(iSec+ipPO-1)*3+1-1)
          Work(ipC+ipPO*3+2-1)=Work(ipGrid+(iSec+ipPO-1)*3+2-1)
          Work(ipC+ipPO*3+3-1)=Work(ipGrid+(iSec+ipPO-1)*3+3-1)
        endif
        iv3=iv3+1
        if (iv3.gt.ive3) then
           iv3=0
           iv2=iv2+1
        endif
        if (iv2.gt.ive2) then
           iv2=0
           iv1=iv1+1
        endif
        ipPO=ipPO+1
c           print *,'ipo',ipP0
            if(ipPO.le.mCoor-1) goto 667
C end of CUBIC box
      else
C coords are given in specific points
* coords for DEBUG mode
       if(isLine.eq.0) then
         do i=0,nCoor-1
           Work(ipC+i*3)=OneCoor(1)+OneCoor(4)*i
           Work(ipC+1+i*3)=OneCoor(2)+OneCoor(5)*i
           Work(ipC+2+i*3)=OneCoor(3)+OneCoor(6)*i
         enddo
       else
c LINE keyword
         do i=0,nCoor-1
            Work(ipC+i*3)=
     &        OneCoor(1)+(OneCoor(4)-OneCoor(1))*i/(nCoor-1)
            Work(ipC+1+i*3)=
     &        OneCoor(2)+(OneCoor(5)-OneCoor(2))*i/(nCoor-1)
            Work(ipC+2+i*3)=
     &        OneCoor(3)+(OneCoor(6)-OneCoor(3))*i/(nCoor-1)
         enddo
       endif
      endif
C end of coordinates.
CVV: FIXME; separate color and sphere
      if(isSphere.eq.1.and.isColor.eq.1) then
c        Call Sphr_Grid(Work(ipCoor),mCoor,Work(ipC),Work(iSphrDist),
c     &       Work(iSphrColor))
      endif
        if(NoOrb.eq.1) then
        Call MOEval(Work(ipMO),nMOs,mCoor,Work(ipC),Work(ipPab),nMOs,
     &              nCoor,iWork(ipDoIt),Debug)

        else

        Call MOEval(Work(ipMO),nMOs,mCoor,Work(ipC),Work(ipCMO),nCMO,
     &              nCoor,iWork(ipDoIt),Debug)

        endif
*
*....... Write out values
*
       nLine=min(nShowMOs,20)
       nSLine=nLine*mCoor
       if(isLine.eq.0) then
         nSLine=1
         nLine=1
        endif
      Call GetMem('Line','Allo','Real',ipLine,nSLine)
       if(levelprint.lt.2) iPrintCount=100
CVV BUG Update ipcutOFF
      Call DumpM2Msi_nosupport(iRun,Luval,nShowMOs,isDensity,nMOs,
     &  iWork(ipGRef), isVb, Work(ipOcc), Work(ipMO),
     &  Work(ipOut), mCoor,
     &  Work(ipVBmat), iGauss, nInc, imoPack, iWork(ipPBlock),
     &  cMoBlock,nBytesPackedVal, dnorm, Crypt, VbOcc,
     &  isTheOne,isLine,isBinary, isEner,iWork(ipType),
     &  iWork(ipNZ),Work(ipE),Work(ipLine),nLine,Work(ipC),
     &  iPrintCount,isDebug,isCutOff,iWork(ipCutOff+iShiftCut),
     &  isSphere,Work(iSphrDist),isColor,Work(iSphrColor))
        if(isUHF.eq.1) then
cVV:
        Call MOEval(Work(ipMO),nMOs,mCoor,Work(ipC),Work(ipCMO_ab),nCMO,
     &              nCoor,iWork(ipDoIt_ab),Debug)

*
*....... Write out values
*
      Call DumpM2Msi_nosupport(iRun,Luval_ab,nShowMOs_ab,isDensity,nMOs,
     &  iWork(ipGRef_ab), isVb, Work(ipOcc_ab), Work(ipMO),
     &  Work(ipOut), mCoor,
     &  Work(ipVBmat), iGauss, nInc, imoPack, iWork(ipPBlock),
     &  cMoBlock,nBytesPackedVal, dnorm, Crypt, VbOcc,
     &  isTheOne,isLine,isBinary, isEner,iWork(ipType),
     &  iWork(ipNZ),Work(ipE_ab),Work(ipLine),nLine,Work(ipC),
     &  iPrintCount,isDebug,isCutOff,iWork(ipCutOff+iShiftCut),
     &  isSphere,Work(iSphrDist),isColor,Work(iSphrColor))
        endif
       Call GetMem('Line','Free','Real',ipLine,nSLine)
        iShiftCut=iShiftCut+mCoor
       if(isVirt.eq.1) then
         do iiMO=1,nMOs
         do jjMO=1,nMOs
         if(Work(ipOcc+iiMO-1).gt.1.1
     &      .and. Work(ipOcc+jjMO-1).lt.0.9) then
c  here if this is a pair Occ-Virt

          do i=1,nMOs
          Work(ipOoo+i-1)=0
          if(i.eq.iiMO) Work(ipOoo+i-1) = 2.0d0-Virt
          if(i.eq.jjMO) Work(ipOoo+i-1) = Virt
          enddo

          call outmo_nosupport(0,2,Work(ipMO),Work(ipOoo),Work(ipOut),
     >               mCoor,nMOs)
          dd=0
          do j=1, mCoor
            dd=dd+Work(ipOut+j-1)
          enddo
c         ddNorm=ddNorm+dd
          call save_ddNorm(dd,iiMO,jjMO,Work(ipdd),nMOs)
         endif
         enddo
         enddo
        endif
       Enddo


      if(iRun.eq.1.and.levelprint.ge.1) then
       Call CollapseOutput(0, ' list of grids ')
      endif
      Write (6,*)
* Check norms
*
      if(isAtom.ne.1) then

        det3=GridAxis1(1)*(GridAxis2(2)*GridAxis3(3)-
     -                     GridAxis2(3)*GridAxis3(2))-
     -        GridAxis1(2)*(GridAxis2(1)*GridAxis3(3)-
     -                     GridAxis2(3)*GridAxis3(1))+
     +        GridAxis1(3)*(GridAxis2(1)*GridAxis3(2)-
     -                     GridAxis2(2)*GridAxis3(1))
        det3=det3/( (iGridNpt(1)-1)*(iGridNpt(2)-1)*(iGridNpt(3)-1))
        dNorm=dNorm*det3

c        print *,'dNorm=',dNorm

       if(isVirt.eq.1) then
        call print_ddNorm(nMOs,Work(ipdd),det3)

c       print *,'ddNorm=',ddNorm*det3

        Call GetMem('ddNo','FREE','Real',ipdd,nMOs*nMOs)
       endif

        Call Add_Info('GRIDIT_NORM',dNorm,1,6)
      endif
*                                                                      *
************************************************************************
*                                                                      *
**      Call DAClos(LuVal)
      if(isLine.eq.0) then
      write(str,'(9i8)') nIrrep, (nBas(i),i=0,nIrrep-1)
      if(isBinary.eq.0) Then
        write(LuVal,'(a)') str
      else
        write(LuVal) str
      endif
c Well, to avoid rewritting of Cerius2 we use old INPORB format
c temporary!
      LuOrb=isFreeUnit(46)
      call molcas_open_ext2(luorb,INPORB,'sequential','formatted',
     & iostat,.false.,irecl,'old',is_error)
c      Open(Unit=LuOrb, File=INPORB, status='old')
c 4001  read(LuOrb,'(a)',err=5000,end=5000) str
c      if(str.ne.'#ORB') goto 4001
5001  read(LuOrb,'(a)',err=5000,end=5000) str
c      if(str(1:1).eq.'#') goto 5001
      iLen=len(str)
      do 5050 i=iLen,1,-1
        if(str(i:i).ne.' ') goto 5060
5050  continue
5060  if(isBinary.eq.0) Then
        write(LuVal,'(a)') str(1:i)
      else
        write(LuVal) str(1:i)
      endif
      goto 5001
5000  close(unit=LuOrb)
      endif ! isLine
      close(unit=LuVal)
      if(isUHF.eq.1) close(unit=LuVal_ab)
      if(isTheOne.eq.1) then
        MM=mCoor-1
        if(MM.gt.10) MM=10
        Call Add_Info('GRIDIT_ONE',Work(ipOut),MM,6)
      endif

*                                                                      *
************************************************************************
*                                                                      *
*...  Epilogue, end
*
c6000  continue

      if(isAtom.eq.1)  then
        mCoor=nCoor

       Call MOEval(Work(ipMO),nMOs,mCoor,Work(ipCoor),Work(ipCMO),
     &             nCMO,nCoor,iWork(ipDoIt),Debug)
       call outmo_nosupport(0,2,Work(ipMO),Work(ipOcc),Work(ipOut),
     &     nCoor,nMOs)
          write(6,'(60a1)') ('*',i=1,60)
      if(ifpartial.eq.0) Then
       write(6,'(a5,3a10,a20)') 'Atom','x','y','z','Density'
      else
       write(6,'(a5,3a10,a20)') 'Atom','x','y','z','Density (partial)'
      endif
          do i=0, nAtoms-1
           write(6,'(a5,3f10.3,e20.10)') AtomLbl(i+1),Work(ipCoor+i*3),
     &   Work(ipCoor+i*3+1),Work(ipCoor+i*3+2),Work(ipOut+i)
          enddo
        Call Add_Info('GRIDIT_ATOM',Work(ipOut),nAtoms,6)


       Call GetMem('Coor','FREE','REAL',ipCoor,3*nSym*nAtoms)

      endif

      Call GetMem('MOValue','Free','Real',ipMO,nInc*nMOs)
      Call GetMem('DOValue','Free','Real',ipOut,nInc)
      if(isSphere.eq.1) then
      Call GetMem('SpDi','Free','Real',iSphrDist,nInc)
      endif
      if(isColor.eq.1) then
      Call GetMem('SpCo','Free','Real',iSphrColor,nInc)
      endif
      if(isWDW.eq.1) call getmem('WDW','FREE','REAL',ipWdW,nCenter)
      Call GetMem('Coor','Free','Real',ipC,nInc*3)
      Call GetMem('CMO','Free','Real',ipCMO,nCMO)
      Call GetMem('iTyp','Free','INTE',ipType,nMOs)
      Call GetMem('Ener','Free','Real',ipE,nMOs)
      Call GetMem('Occu','Free','Real',ipOcc,nMOs)
      Call GetMem('Occ2','Free','Real',ipOoo,nMOs)
      Call GetMem('Sort','Free','Inte',ipSort,nMOs)
      Call GetMem('Nzer','Free','Inte',ipNZ,nMOs*2)
      Call GetMem('Vol','Free','real',ipVol,nMOs)
      Call GetMem('NRef','Free','Inte',ipGRef,nMOs)
      Call GetMem('DoIt','Free','Inte',ipDoIt,nMOs)
      Call GetMem('Pab','Free','Real',ipPab,nMOs)
      if(isUHF.eq.1) then
      Call GetMem('CMO_ab','Free','Real',ipCMO_ab,nCMO)
      Call GetMem('Ener_ab','Free','Real',ipE_ab,nMOs)
      Call GetMem('Occu_ab','Free','Real',ipOcc_ab,nMOs)
      Call GetMem('Sort_ab','Free','INTE',ipSort_ab,nMOs)
      Call GetMem('NRef_ab','Free','INTE',ipGRef_ab,nMOs)
      Call GetMem('DoIt_ab','Free','INTE',ipDoIt_ab,nMOs)
      endif
      if(isVB.eq.1)then
      Call GetMem('VBmat','Free','Real',ipVBmat,nMOs*nActOrb)
      else
      Call GetMem('VBmat','Free','Real',ipVBmat,1)
      endif
      if(isUserGrid.eq.1)
     *  Call GetMem('Grid','FREE','REAL',ipGrid,nGridPoints*3)
      if (imoPack .ne. 0) then
        Call GetMem('PackedBlock','Free','INTE',ipPBlock,nInc)
      else
        Call GetMem('PackedBlock','Free','INTE',ipPBlock,1)
      endif
       if (isAtom.eq.0)
     &   Call GetMem('Coor','FREE','REAL',ipCoor,3*nSym*nAtoms)

       if(isCutOff.eq.1)
     &      Call GetMem('CUTFL','Free','Inte',ipCutOff,nCoor)
*
************************************************************************
*                                                                      *
*     Call GetMem('DrvMO_E','Check','Real',iDum,iDum)

      Call qExit('DrvMO')
      Return
2000  write(6,*) 'Error during open InpOrb file'
      Call Quit_OnUserError()
1000  write(6,*) 'Error during read from InpOrb file'
      Call Quit_OnUserError()
3000  Write (6,*) 'Section not found in INPORB file'
      Call Quit_OnUserError()
      End

      Subroutine save_ddNorm(ddNorm,ii,jj,Wdd,nMOs)
      Implicit Real*8 (A-H,O-Z)
      Dimension Wdd(nMOs,nMOs)
      Wdd(ii,jj)=Wdd(ii,jj)+ddNorm
      return
      end

      subroutine print_ddNorm(nMOs,Wdd,det3)
      Implicit Real*8 (A-H,O-Z)
      Dimension Wdd(nMOs,nMOs)
        do i=1,nMOs
        write (6,'(20F8.3)') (Wdd(i,j)*det3,j=1,nMOs)
        enddo
      return
      end
