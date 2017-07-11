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
* Copyright (C) 1995, Roland Lindh                                     *
*               2000, Valera Veryazov                                  *
************************************************************************
      SubRoutine DrvMO(iRun,INPORB)
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
*     The ag comments hits on that Alexander Gaenko has modified the   *
*     code too.                                                        *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "grid.fh"
#ifdef _FDE_
#include "embpotdata.fh"
#endif
#include "para_info.fh"

      Logical Debug, long_prt
      Logical is_error
      Character str*128, Crypt*7
      Character INPORB*(*)
c      Character bb
c      character namepack*3
      Character*80 myTitle
      character*128 line
      Character*64 status
      Real*8 pp(3)
      Real*8 dInteg(3),maxes(3,3)
      Integer nTypes(7)
      integer nSym
c      Character *10 Pack, oPack
c      Dimension xLimits(4)
c      Integer iYDelta(3)
c      Real tmp(3)
      Data Crypt/'fi123sd'/
*
*---- Set size of batches
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

      if(isUHF.ne.0.and.(isDerivative.ne.0.or.isCurDens.ne.0)) then
        write(6,*) "ERROR - Current density or derivatives not
     &              implemented for UHF!!!"
        call ABEND()
      end if
      if(NoOrb.ne.0.and.(isDerivative.ne.0.or.isCurDens.ne.0)) then
        write(6,*) "ERROR - Current density or derivatives not
     &              implemented with the NOORB keyword!!!"
        call ABEND()
      end if
      if(isAtom.ne.0.and.(isDerivative.ne.0.or.isCurDens.ne.0)) then
        write(6,*) "ERROR - Current density or derivatives not
     &              implemented with the ATOM keyword!!!"
        call ABEND()
      end if


      if(isCurDens.gt.0) then
        if(iuseMaxes.gt.0) then
c read the magnetic axes from a file (if exists)
          LUMAXES=71
          LUMAXES=IsFreeUnit(LUMAXES)
          CALL Molcas_Open(LUMAXES,'SODIAG.MAXES')
          READ(LUMAXES,*) MAXES
          CLOSE(LUMAXES)
        else
          CALL DCOPY_(9,0.0d0,0,maxes,1)
          CALL DCOPY_(9,1.0d0,0,maxes,4)
        end if

        write(6,*) "Magnetic axes:"
        write(6,*) (maxes(1,k),k=1,3)
        write(6,*) (maxes(2,k),k=1,3)
        write(6,*) (maxes(3,k),k=1,3)

        Call GetMem('CMOR','Allo','Real',ipCMOR,nCMO)
        Call GetMem('CMOI','Allo','Real',ipCMOI,nCMO)
        dInteg(1) = 0.0d0
        dInteg(2) = 0.0d0
        dInteg(3) = 0.0d0
      else
        Call GetMem('CMO','Allo','Real',ipCMO,nCMO)
      end if

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


      if(isUHF.eq.0) then

      if(isCurDens.gt.0) then
        Call RdVec(INPORB,LuOrb,'COE',nIrrep,NBAS,NBAS,
     &    Work(ipCMOR),Work(ipOcc),Work(ipE),iDummy,
     &    myTitle,0,iErr)
        Call RdVec(INPORB(1:LEN_TRIM(INPORB))//'.I',LuOrb,'COE',nIrrep,
     &    NBAS,NBAS,Work(ipCMOI),Work(ipOcc),Work(ipE),iDummy,
     &    myTitle,0,iErr)
c        write(6,*) "From INPORB"
c        write(6,*) Work(ipCMOR:ipCMOR-1+nCMO)
c        write(6,*) "From INPORB.I"
c        write(6,*) Work(ipCMOI:ipCMOI-1+nCMO)
      else
        Call RdVec(INPORB,LuOrb,'COE',nIrrep,NBAS,NBAS,
     &    Work(ipCMO),Work(ipOcc),Work(ipE),iDummy,
     &    myTitle,0,iErr)
c        write(6,*) "From INPORB"
c        write(6,*) Work(ipCMO:ipCMO-1+nCMO)
      end if


c construct Pab
       if(NoOrb.eq.1) then
c        print *,'nCMO,nMOs', nCMO,nMOs
        Call makePab(Work(ipCMO),Work(ipOcc),Work(ipPab),nMOs,nMOs,
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
          iWork(ipType+j)=0
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
          iWork(ipType+j)=0
          enddo

      endif
      do j=1,7
      nTypes(j)=0
      enddo
      do j=0, nMOs-1
      jj=iWork(ipType+j)
      if(jj.gt.0) nTypes(jj)=nTypes(jj)+1
      enddo

      close(LuOrb)
*
*  Calculate net.
*
       if(iRun.eq.1) then
        Write(6,'(A)') '   Input vectors read from INPORB'
        Write(6,'(A,A)') '   Orbital file label: ',
     &                   myTitle(:mylen(myTitle))
       endif
      do j=1,3
       iGridNpt(j)=iGridNpt(j)+1
      enddo

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

       Call PickOrb(ipNz,ipSort,ipGref,ipSort_ab,
     &  ipGref_ab,ipVol,ipE,ipOcc,ipE_ab,ipOcc_ab,
     &  nShowMOs,nShowMOs_ab,isener,nMOs,myTitle,ipType)
*                                                                      *
*---- Start run over sets of grid points
*
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

      if(isCurDens.gt.0) then
        Call GetMem('MOValDR','Allo','Real',ipMODR,4*nInc*nMOs)
        Call GetMem('MOValDI','Allo','Real',ipMODI,4*nInc*nMOs)
        Call GetMem('DOValue','Allo','Real',ipOutR,4*nInc)
        Call GetMem('DOValue','Allo','Real',ipOutI,4*nInc)
      else
        Call GetMem('MOValue','Allo','Real',ipMO,nInc*nMOs)
        Call GetMem('DOValue','Allo','Real',ipOut,nInc)
        if(isXField.eq.1)
     &     Call GetMem('DOValXF','Allo','Real',ipOutXF,nInc)
      end if

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
      if(isXField.eq.1) Call GetMem('Wght','Allo','Real',ipWt,nInc)
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
c  FIXME!
c
      iCRSIZE=1
      NBYTES=10
      NINLINE=10

      call PrintHeader(nMOs,nShowMOs,nShowMOs_ab,nCoor,nInc,
     & iiCoord,nTypes,iCRSIZE,NBYTES,NINLINE,nBlocks )
      Call PrintTitles(LuVal,nShowMOs,isDensity,nMOs,
     &  iWork(ipGRef), isEner, isVB, Work(ipOcc), iWork(ipType), Crypt,
     &  iWork(ipNZ), Work(ipE), VBocc, ifpartial,isLine,isSphere,
     &  isColor, ISLUSCUS,ncoor,nBlocks,nInc)
      if(isUHF.eq.1) then
      Call PrintTitles(LuVal_ab,nShowMOs_ab,isDensity,nMOs,
     &  iWork(ipGRef_ab), isEner, isVB, Work(ipOcc_ab),
     &  iWork(ipType), Crypt,
     &  iWork(ipNZ), Work(ipE_ab), VBocc, ifpartial,isLine,isSphere,
     &  isColor, ISLUSCUS,ncoor,nBlocks,nInc)
      endif
*                                                                      *

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
      dtot=0.d0
c      if(isCutOff.eq.1) nCoor=iiCoord
      if(isXField.eq.1) then
         LuXFieldFile=18
         LuXFieldFile=IsFreeUnit(LuXFieldFile)
         Call Molcas_Open(LuXFieldFile,'CHARGEFILE')
           if (debug) write(6,*) 'XFIELD:  dtot = ', dtot
                      write(6,*) 'XFIELD:  dtot = ', dtot
         Call Get_iScalar('Unique atoms',nAtoms)
         Call Allocate_Work(ipCU2,3*nAtoms)
         Call Allocate_Work(ipCMu2, nAtoms)
         Call Get_dArray('Unique Coordinates',Work(ipCU2),3*nAtoms)
         Call Get_dArray('Nuclear charge',    Work(ipCMu2), nAtoms)
c         write(LuXFieldFile,'(I20)') nAtoms+nCoor
c ensure that the nuclear charges are added only on the master node,
c when running in parallel
      IF (MyRank.eq.0) THEN
         do i=1, nAtoms
            write(LuXFieldFile,'(3ES24.14,2x,ES24.14,A)')
     &      Work(ipCU2+3*i-3), Work(ipCU2+3*i-2),Work(ipCU2+3*i-1),
     &      Work(ipCMu2+i-1),' 0.0 0.0 0.0'
         enddo
      END IF
         Call Free_Work(ipCMu2)
         Call Free_Work(ipCU2)
      endif

cccccccccccccc  main loop starts here  ccccccccccccccccccccccccccccccccc
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
      if(isXField.eq.1) call dcopy_(mCoor,0.d0,0,Work(ipOutXF),1)

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
          !make a local copy of the weights of the corresponding grid points:
          if(isXField.eq.1) Work(ipWt+ipPO)=Work(ipWeight+iSec+ipPO-1)
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
        Call Sphr_Grid(Work(ipCoor),mCoor,Work(ipC),Work(iSphrDist),
     &       Work(iSphrColor))
      endif

        if(NoOrb.eq.1) then
        nDrv=0
        Call MOEval(Work(ipMO),nMOs,mCoor,Work(ipC),Work(ipPab),nMOs,
     &              nCoor,iWork(ipDoIt),nDrv,1,Debug)
        else if(isCurDens.gt.0) then
           Call MOEvalDel(Work(ipMODR),
     &                    nMOs,mCoor,Work(ipC),Work(ipCMOR),
     &                    nCMO,nCoor,iWork(ipDoIt),Debug)
           Call MOEvalDel(Work(ipMODI),
     &                    nMOs,mCoor,Work(ipC),Work(ipCMOI),
     &                    nCMO,nCoor,iWork(ipDoIt),Debug)
        else if(isDerivative.gt.0) then
           Call MOEvalDer(Work(ipMO),isDerivative,
     &                    nMOs,mCoor,Work(ipC),Work(ipCMO),
     &                    nCMO,nCoor,iWork(ipDoIt),Debug)
        else
           nDrv=0
           Call MOEval(Work(ipMO),nMOs,mCoor,Work(ipC),Work(ipCMO),nCMO,
     &                 nCoor,iWork(ipDoIt),nDrv,1,Debug)
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

      if(isCurDens.gt.0) then
#ifdef _UNUSED_
        Call DumpM2MsiCD(iRun,Luval,LID,nShowMOs,isDensity,nMOs,
     &    iWork(ipGRef), isVb, Work(ipOcc), Work(ipMODR), Work(ipMODI),
     &    Work(ipOutR),Work(ipOutI), mCoor,
     &    Work(ipVBmat), iGauss, nInc, imoPack, iWork(ipPBlock),
     &    cMoBlock,nBytesPackedVal, dnorm, Crypt, VbOcc,
     &    isTheOne,isLine,isBinary, isEner,iWork(ipType),
     &    iWork(ipNZ),Work(ipE),Work(ipLine),nLine,Work(ipC),
     &    iPrintCount,isDebug,isCutOff,iWork(ipCutOff+iShiftCut),
     &    isSphere,Work(iSphrDist),isColor,Work(iSphrColor),
     &    isCurDens,isRxJ,dInteg,maxes,isLuscus,NBYTES,NINLINE)
#endif
      else
        Call DumpM2Msi(iRun,Luval,LID,nShowMOs,isDensity,nMOs,
     &    iWork(ipGRef), isVb, Work(ipOcc), Work(ipMO),
     &    Work(ipOut), mCoor,
     &    Work(ipVBmat), iGauss, nInc, imoPack, iWork(ipPBlock),
     &    cMoBlock,nBytesPackedVal, dnorm, Crypt, VbOcc,
     &    isTheOne,isLine,isBinary, isEner,iWork(ipType),
     &    iWork(ipNZ),Work(ipE),Work(ipLine),nLine,Work(ipC),
     &    iPrintCount,isDebug,isCutOff,iWork(ipCutOff+iShiftCut),
     &    isSphere,Work(iSphrDist),isColor,Work(iSphrColor),
     &    isLuscus,NBYTES,NINLINE)
c
c      if(isXField.eq.1) call dcopy_(mCoor,Work(ipOut),1,Work(ipOutXF),1)
          if(isXField.eq.1) then
            do i=1,mCoor
            Work(ipOutXF+i-1)=Work(ipOut+i-1)
            enddo
          endif
        if(isUHF.eq.1) then
cVV:
          nDrv=0
          Call MOEval(Work(ipMO),nMOs,mCoor,Work(ipC),Work(ipCMO_ab),
     &                nCMO,nCoor,iWork(ipDoIt_ab),nDrv,1,Debug)

*....... Write out values
*
          Call DumpM2Msi(iRun,Luval_ab,LID,nShowMOs_ab,isDensity,nMOs,
     &      iWork(ipGRef_ab), isVb, Work(ipOcc_ab), Work(ipMO),
     &      Work(ipOut), mCoor,
     &      Work(ipVBmat), iGauss, nInc, imoPack, iWork(ipPBlock),
     &      cMoBlock,nBytesPackedVal, dnorm, Crypt, VbOcc,
     &      isTheOne,isLine,isBinary, isEner,iWork(ipType),
     &      iWork(ipNZ),Work(ipE_ab),Work(ipLine),nLine,Work(ipC),
     &      iPrintCount,isDebug,isCutOff,iWork(ipCutOff+iShiftCut),
     &      isSphere,Work(iSphrDist),isColor,Work(iSphrColor),
     &      isLuscus,NBYTES,NINLINE)

          if(isXField.eq.1) then
            do i=1,mCoor
            Work(ipOutXF+i-1)=Work(ipOutXF+i-1)+Work(ipOut+i-1)
            enddo
          endif
        endif ! end if(isUHF.eq.1)

        if(isXField.eq.1) then
            do i=1,mCoor
              write(LuXFieldFile,'(3ES24.14,2x,ES24.14,A)')
     &        Work(ipC+3*i-3),Work(ipC+3*i-2),Work(ipC+3*i-1), ! cartesian coordinate of the grid point
     &       -Work(ipOutXF+i-1)*Work(ipWt +i-1),  ! electronic charge on it
     &        ' 0.0 0.0 0.0'
              dtot=dtot+Work(ipOutXF+i-1)*Work(ipWt +i-1)  ! accumulate total density, for debug
            enddo
            call xFlush(LuXFieldFile)
        endif

      endif !end if(isCurDens.gt.0)
c

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

          call outmo(0,2,Work(ipMO),Work(ipOoo),Work(ipOut),
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
       Enddo !iSec
cccccccccccccc  main loop ends here  ccccccccccccccccccccccccccccccccc


      if(iRun.eq.1.and.levelprint.ge.1) then
       Call CollapseOutput(0, ' list of grids ')
      endif
      Write (6,*)
* Check norms
*
      if(isXField.eq.1) then
      call xFlush(LuXFieldFile)
      close(LuXFieldFile)
      endif

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
        IF (ISLUSCUS .EQ. 1) THEN
            LINE=' '
            CALL PRINTLINE(LUVAL,LINE,1,0)
            LINE=' </DENSITY>'
            CALL PRINTLINE(LUVAL,LINE,11,0)
            LINE=' <BASIS>'
            CALL PRINTLINE(LUVAL,LINE,8,0)
            LINE=' </BASIS>'
            CALL PRINTLINE(LUVAL,LINE,9,0)
            LINE=' <INPORB>'
            CALL PRINTLINE(LUVAL,LINE,9,0)
          endif


      if(isLine.eq.0) then
      write(str,'(9i8)') nIrrep, (nBas(i),i=0,nIrrep-1)
      CALL PrintLine(luval, STR, 72,0)
C      if(isBinary.eq.0) Then
C        write(LuVal,'(a)') str
C      else
C        write(LuVal) str
C      endif
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
 5060 CONTINUE
      CALL PrintLine(luval, STR, I,0)
C5060  if(isBinary.eq.0) Then
C        write(LuVal,'(a)') str(1:i)
C      else
C        write(LuVal) str(1:i)
C      endif
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

#ifdef _FDE_
      ! Embedding module, output data needed for polarization of the
      ! environment
      ! Thomas Dresselhaus
      if (embPot.and.
     &    (embWriteDens.or.embWriteGrad.or.embWriteHess)) then
       if (embWriteHess) then
         mAO = 10
       else if (embWriteGrad) then
         mAO = 4
       else
         mAO = 1
       end if
       Call EmbPotInit(.true.)
       ! Calculate density on the grid
       Call GetMem('MyMOValue','ALLO','REAL',ipMyMOval,
     &             mAO*nMOs*nEmbGridPoints)
       do i=0, (nMOs-1)
         iWork(ipDoIt+i)=1
       end do
       if (embWriteDens) then
        Call GetMem('Density','ALLO','REAL',ipDensGrid,nEmbGridPoints)
        nDrv=0
       end if
       if (embWriteGrad) then
        Call GetMem('Gradient','ALLO','REAL',ipGradGrid,
     &              3*nEmbGridPoints)
        nDrv=1
       end if
       if (embWriteHess) then
        Call GetMem('Hessian','ALLO','REAL',ipHessGrid,6*nEmbGridPoints)
        nDrv=2
       end if
       ! Get the values of the MOs at the grid points to use them later
       ! to calculate the density
       Call MOEval(Work(ipMyMOval),nMOs,nEmbGridPoints,
     &             Work(posEmbGridCoord),
     &             Work(ipCMO),nCMO,nCoor,iWork(ipDoIt),nDrv,mAO,Debug)
       sumDensEmb=0.0d0
       if (embWriteDens) then
        iUnitEmb = isFreeUnit(1)
        call molcas_open(iUnitEmb, embOutDensPath)
       end if
       if (embWriteGrad) then
         iUnitEmbGrad = isFreeUnit(1)
         call molcas_open(iUnitEmbGrad, embOutGradPath)
       end if
        if (embWriteHess) then
         iUnitEmbHess = isFreeUnit(1)
         call molcas_open(iUnitEmbHess, embOutHessPath)
       end if
       do i=0, (nEmbGridPoints-1)
        if (embWriteDens) then
         Work(ipDensGrid+i)=0.0d0
        end if
        if (embWriteGrad) then
         Work(ipGradGrid+i*3)=0.0d0
         Work(ipGradGrid+i*3+1)=0.0d0
         Work(ipGradGrid+i*3+2)=0.0d0
        end if
        if (embWriteHess) then
         Work(ipHessGrid+i*6)=0.0d0
         Work(ipHessGrid+i*6+1)=0.0d0
         Work(ipHessGrid+i*6+2)=0.0d0
         Work(ipHessGrid+i*6+3)=0.0d0
         Work(ipHessGrid+i*6+4)=0.0d0
         Work(ipHessGrid+i*6+5)=0.0d0
        end if
        ! Sum up values for all MOs to get the total density/grad
        do j=0, (nMOs-1)
         if (embWriteDens) then
          Work(ipDensGrid+i) = Work(ipDensGrid+i)+
     &                        Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO)*
     &                        Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO)*
     &                        Work(ipOcc+j)
         end if
         if (embWriteGrad) then
          ! grad(rho(r)) = 2 * sum_i(phi_i(r) * grad(phi_i(r)))
          Work(ipGradGrid+i*3) = Work(ipGradGrid+i*3)+
     &                       Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO)*
     &                       Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+1)*
     &                       Work(ipOcc+j)*2
          Work(ipGradGrid+i*3+1) = Work(ipGradGrid+i*3+1)+
     &                       Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO)*
     &                       Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+2)*
     &                       Work(ipOcc+j)*2
          Work(ipGradGrid+i*3+2) = Work(ipGradGrid+i*3+2)+
     &                       Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO)*
     &                       Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+3)*
     &                       Work(ipOcc+j)*2
         end if
         if (embWriteHess) then
          ! Hess(rho(r)) = 2 * sum_i(phi_i(r) * Hess(phi_i(r)) +
          !                         grad(phi_i(r)) * (grad(phi_i(r)))*)
          Work(ipHessGrid+i*6) = Work(ipHessGrid+i*6)+
     &                    (Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO)*
     &                     Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+4) +
     &                    (Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+1)*
     &                     Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+1)))*
     &                     Work(ipOcc+j)*2
          Work(ipHessGrid+i*6+1) = Work(ipHessGrid+i*6+1)+
     &                    (Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO)*
     &                     Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+5) +
     &                    (Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+1)*
     &                     Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+2)))*
     &                     Work(ipOcc+j)*2
          Work(ipHessGrid+i*6+2) = Work(ipHessGrid+i*6+2)+
     &                    (Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO)*
     &                     Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+6) +
     &                    (Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+1)*
     &                     Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+3)))*
     &                     Work(ipOcc+j)*2
          Work(ipHessGrid+i*6+3) = Work(ipHessGrid+i*6+3)+
     &                    (Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO)*
     &                     Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+7) +
     &                    (Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+2)*
     &                     Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+2)))*
     &                     Work(ipOcc+j)*2
          Work(ipHessGrid+i*6+4) = Work(ipHessGrid+i*6+4)+
     &                    (Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO)*
     &                     Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+8) +
     &                    (Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+2)*
     &                     Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+3)))*
     &                     Work(ipOcc+j)*2
          Work(ipHessGrid+i*6+5) = Work(ipHessGrid+i*6+5)+
     &                    (Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO)*
     &                     Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+9) +
     &                    (Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+3)*
     &                     Work(ipMyMOval+(i+j*nEmbGridPoints)*mAO+3)))*
     &                     Work(ipOcc+j)*2
         end if
        end do
        ! Write to file
        if (embWriteDens) then
          write(iUnitEmb,'(e18.10,2X)') Work(ipDensGrid+i)
        end if
        if (embWriteGrad) then
          write(iUnitEmbGrad,'(3(e18.10,2X))') Work(ipGradGrid+i*3),
     &             Work(ipGradGrid+i*3+1), Work(ipGradGrid+i*3+2)
        end if
        if (embWriteHess) then
          write(iUnitEmbHess,'(6(e18.10,2X))') Work(ipHessGrid+i*6),
     &             Work(ipHessGrid+i*6+1), Work(ipHessGrid+i*6+2),
     &             Work(ipHessGrid+i*6+3), Work(ipHessGrid+i*6+4),
     &             Work(ipHessGrid+i*6+5)
        end if
        sumDensEmb=sumDensEmb+Work(posEmbWeight+i)*Work(ipDensGrid+i)
       end do
       if (embWriteDens) close(iUnitEmb)
       if (embWriteGrad) close(iUnitEmbGrad)
       if (embWriteHess) close(iUnitEmbHess)
       write(6,*) "EMBPOT: Integrated density on the grid=", sumDensEmb
       ! Tidy up
       Call GetMem('MOValue','FREE','REAL',ipMyMOval,
     &             mAO*nMOs*nEmbGridPoints)
       if (embWriteDens) then
        Call GetMem('Density','FREE','REAL',ipDensGrid,nEmbGridPoints)
       end if
       if (embWriteGrad) then
        Call GetMem('Gradient','FREE','REAL',ipGradGrid,
     &              3*nEmbGridPoints)
       end if
       if (embWriteHess) then
        Call GetMem('Hessian','FREE','REAL',ipHessGrid,6*nEmbGridPoints)
       end if
       Call embPotFreeMem
      end if
#endif
*                                                                      *
************************************************************************
*                                                                      *
*...  Epilogue, end
*
c6000  continue

      if(isAtom.eq.1)  then
        mCoor=nCoor
        nDrv=0
        Call MOEval(Work(ipMO),nMOs,mCoor,Work(ipCoor),Work(ipCMO),
     &             nCMO,nCoor,iWork(ipDoIt),nDrv,1,Debug)
        call outmo(0,2,Work(ipMO),Work(ipOcc),Work(ipOut),nCoor,nMOs)
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

      IF (ISLUSCUS .EQ. 1) THEN
        CALL PRTLUSENDGRID()
      END IF

      if(isSphere.eq.1) then
      Call GetMem('SpDi','Free','Real',iSphrDist,nInc)
      endif
      if(isColor.eq.1) then
      Call GetMem('SpCo','Free','Real',iSphrColor,nInc)
      endif
      if(isWDW.eq.1) call getmem('WDW','FREE','REAL',ipWdW,nCenter)
      Call GetMem('Coor','Free','Real',ipC,nInc*3)

      if(isCurDens.gt.0) then
        Call GetMem('CMOR','Free','Real',ipCMOR,nCMO)
        Call GetMem('CMOI','Free','Real',ipCMOI,nCMO)
        Call GetMem('MOValDR','Free','Real',ipMODR,4*nInc*nMOs)
        Call GetMem('MOValDI','Free','Real',ipMODI,4*nInc*nMOs)
        Call GetMem('DOValue','Free','Real',ipOutR,4*nInc)
        Call GetMem('DOValue','Free','Real',ipOutI,4*nInc)

        write(6,*) "CURRENT DENSITY INTEGRAL:"
        boxSize=GridAxis1(1)/(1.0d0*(iGridNpt(1)-1))
     &         *GridAxis2(2)/(1.0d0*(iGridNpt(2)-1))
     &         *GridAxis3(3)/(1.0d0*(iGridNpt(3)-1))
        write(6,*) "GridAxis1,2,3"
        write(6,*) GridAxis1,GridAxis2,GridAxis3
        write(6,*) "iGridNpt"
        write(6,*) iGridNpt
        write(6,*) "Box Size: ",boxSize

        write(6,*) "X: ",dInteg(1)*boxSize
        write(6,*) "Y: ",dInteg(2)*boxSize
        write(6,*) "Z: ",dInteg(3)*boxSize

        call Add_info("CURDINT",dInteg, 3, 5)

      else
        Call GetMem('CMO','Free','Real',ipCMO,nCMO)
        Call GetMem('MOValue','Free','Real',ipMO,nInc*nMOs)
        Call GetMem('DOValue','Free','Real',ipOut,nInc)
      end if

        Call GetMem('iTyp','Free','INTE',ipType,nMOs)
        Call GetMem('Ener','Free','Real',ipE,nMOs)
        Call GetMem('Occu','Free','Real',ipOcc,nMOs)
        Call GetMem('Occ2','Free','Real',ipOoo,nMOs)
        Call GetMem('Sort','Free','Inte',ipSort,nMOs)
        Call GetMem('Nzer','Free','Inte',ipNZ,nMOs*2)
        Call GetMem('Vol' ,'Free','real',ipVol,nMOs)
        Call GetMem('NRef','Free','Inte',ipGRef,nMOs)
        Call GetMem('DoIt','Free','Inte',ipDoIt,nMOs)
        Call GetMem('Pab' ,'Free','Real',ipPab,nMOs)
      if(isUHF.eq.1) then
        Call GetMem('CMO_ab' ,'Free','Real',ipCMO_ab,nCMO)
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
      if(isUserGrid.eq.1.and.isXfield.eq.1)
     *  Call GetMem('Weight','FREE','REAL',ipWeight,nGridPoints)
      if (imoPack .ne. 0) then
        Call GetMem('PackedBlock','Free','INTE',ipPBlock,nInc)
      else
        Call GetMem('PackedBlock','Free','INTE',ipPBlock,1)
      endif
      if (isAtom.eq.0.and.isXField.eq.0)
     &  Call GetMem('Coor','FREE','REAL',ipCoor,3*nSym*nAtoms)

      if(isCutOff.eq.1)
     &  Call GetMem('CUTFL','Free','Inte',ipCutOff,nCoor)
      if(isXField.eq.1) then
        Call GetMem('Wght','Free','Real',ipWt,nInc)
        Call GetMem('DOValXF','Free','Real',ipOutXF,nInc)
      endif
*
************************************************************************
*                                                                      *
*     Call GetMem('DrvMO_E','Check','Real',iDum,iDum)

      Call qExit('DrvMO')
      Return
c2000  write(6,*) 'Error during open InpOrb file'
c      Call Quit_OnUserError()
c1000  write(6,*) 'Error during read from InpOrb file'
c      Call Quit_OnUserError()
c3000  Write (6,*) 'Section not found in INPORB file'
c      Call Quit_OnUserError()
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
