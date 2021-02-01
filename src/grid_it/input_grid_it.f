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
      SubRoutine Input_Grid_It(iRun,INPORB,iReturn)
************************************************************************
* Adapted from SAGIT to work with OpenMolcas (October 2020)            *
************************************************************************
*                                                                      *
* Object: input module for grid code                                   *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "grid.fh"

      Character Key*80
      Character INPORB*(*)
      Integer iTemp(3)
      Dimension dTemp(4)
      Character AllKeys*265
      Character SelectStr*120
      Character FileStr*256, FileIn*256
      Character MULLprt*80
*
      AllKeys=
     *'PRIN BINA ASCI NPOI DENS '//
     *'SPAR ORBI REGI ONE  TITL '//
     *'GAP  END  NODE TOTA NAME '//
     *'VB   ALL  ATOM CUBE GRID '//
     *'PACK PKLI PKBI NOOR LINE '//
     *'ORAN ERAN DEBU CUTO NOPA '//
     *'GORI SELE NOSO FILE SPHR '//
     *'COLO VIRT MULL SUBB XDER '//
     *'YDER ZDER GDER CURD CRXJ '//
     *'UMAX NOLU XFIE LUS1 LUS2 '//
     *'PLUS MINU XFMI'
*
c      Do 108 i = 1, nRout
c         nPrint(i) = 5
c 108  Continue
      isBinary=3
      isAuto=1
      isNet=0
      isReadNet=0
      isTheOne=0
      Title1=' '
      TheName=' '
      TheGap=4.0d0
      nMOmin=0
      nGrid=5
      nReq=-1
      isDensity=1
      isDerivative=0
      isCurDens=0
      isRxJ=0
      iuseMaxes=0
      isAuMO=-1
      isAtom=0
      isTotal=0
c      isVB=0
      isAll=0
      isUserGrid=0
      iGauss=0
      NoOrb=0
      iMaxUp=7
      iMaxDown=7
      isLine=0
      itRange=1
      isDebug=0
      isCutOff=0
      CutOff=2.5D0
      iCustOrig=0
      NoSort=0
      isFileOrb=0
      isSphere=0
      isColor=0
      isVirt=0
      isMULL=0
      isLONGPRT=0
      isWDW=0
      iSubBlock=0
      isLuscus=1
c      isLusMath=0
c      aLusMath=-1
!not preparing the GRIDCHARGE file as external source for XFIELD input
      isXField=0
      XFminCh=0.d0
c Default values for packing
      imoPack=0
      isBinPack=0
      xLeft=0.005d0
      xRight=0.7d0
c   (really, half range:)
      iyRange=128
      nBytesPackedVal=1
      iMinYLeft=4
      xLoErr=0.10D0
      xHiErr=0.25D0
      INPORB='INPORB'
      if(iRun.eq.0) then
c make defaults for a fast run via call from other module
        isNet=1  ! set sparse
c        isBinary=0 ! temporary set Ascii output
c       iMaxUp=1
c       iMaxDown=5
c        imoPack=0 ! packed grids not working with 64bit (?)
c        isCutOff=1
c        goto 500
      endif

*
*     KeyWord directed input
*
      InUnit=5
*      Function MyGetKey
*     * (InUnit, What, IValue, RValue, SValue, N, IArray, RArray)

      Call RdNLst(InUnit,'GRID_IT')
 998  if(MyGetKey(InUnit,'S',iD,rD,Key,iD,[iD],[rD]).ne.0) goto 997
      iKey=index(AllKeys,Key(1:4))
      if(iKey.eq.0.or.(iKey-1)/5*5.ne.(iKey-1)) then
      write(6,'(a,a)') 'Unrecognized keyword in input file:', Key(1:4)
        Call Quit_OnUserError()
      endif
      iKey=(iKey-1)/5+1
*
      if(iKey.eq.1) then
* PRIN
      if(MyGetKey(InUnit,'I',n,rD,Key,iD,[iD],[rD]).ne.0) goto 666
      Do j = 1, n
      if(MyGetKey(InUnit,'A',iD,rD,Key,2,iTemp,[rD]).ne.0) goto 666
c      print *,'debug'
c         nPrint(iTemp(1))=iTemp(2)
      enddo
      endif
      if(iKey.eq.2) then
* BINARY = default
        isBinary=1
      endif
      if(iKey.eq.3) then
* ASCII = for debug
c       Write(6,*) ' Keyword ASCII is obsolete'
c       Write(6,*) ' It can be used only for debugging purpose'
c       Write(6,*) ' Note that .lus files produced with this option '
c       Write(6,*) '      can not be visualised'
        isBinary=0
      endif
      if(iKey.eq.4) then
* NPOI
      if(MyGetKey(InUnit,'A',iD,rD,Key,3,iGridNpt,[rD]).ne.0) goto 666
      isNet=-1
      isReadNet=isReadNet+1
      endif
      if(iKey.eq.5) then
* DENSE - dense grid network..
        isNet=2
        isReadNet=isReadNet+1
      endif
      if(iKey.eq.6) then
* SPARSE - rare grid network..
        isNet=1
        isReadNet=isReadNet+1
      endif
      if(iKey.eq.7) then
* ORBI Orbitals
      if(nReq.gt.0) then
        write(6,*) 'ORBI keyword can not be used together with SELEct'
        call Quit_OnUserError()
      endif
      if(MyGetKey(InUnit,'I',nReq,rD,Key,iD,[iD],[rD]).ne.0) goto 666

      if(nReq.gt.MAXGRID) then
      write(6,'(a,i5,a,i5)')
     *      'Too many requested orbitals ',nReq,'>',MAXGRID
      Call Quit_OnUserError()
        endif
       read(inUnit,*, err=666, end=666) (iReq(i),i=1,nReq*2)
cc      if(MyGetKey(InUnit,'A',iD,rD,Key,nReq*2,iReq,[rD]).ne.0) goto 666
      isAuMO=0
      endif
      if(iKey.eq.8) then
* REGION
      if(MyGetKey(InUnit,'D',iD,rD,Key,2,[iD],Region).ne.0) goto 666
      itRange=1
      isAuMO=1
      write(6,*) ' *** Warning keyword REGION is obsolete'
      write(6,*) ' ***         assumimg Energy range '
      endif
      if(iKey.eq.9) then
* ONE - debug option
      if(MyGetKey(InUnit,'D',iD,rD,Key,7,[iD],OneCoor).ne.0) goto 666
      isTheOne=1
      isBinary=0
      endif
      if(iKey.eq.10) then
* TITLE
* NOTE: Title can be only ONE line here!!!
      if(MyGetKey(InUnit,'S',iD,rD,Title1,iD,[iD],[rD]).ne.0) goto 666
      endif
      if(iKey.eq.11) then
* GAP
      if(MyGetKey(InUnit,'R',iD,TheGap,Key,iD,[iD],[rD]).ne.0) goto 666
      endif

      if(iKey.eq.12) then
* END
      goto 997
      endif
      if(iKey.eq.13) then
* NODENSITY
      isDensity=0
      endif
      if(iKey.eq.14) then
* TOTAL
      isTotal=1
      endif
      if(iKey.eq.15) then
* NAME
      read(InUnit,'(a)') TheName
c      if(MyGetKey(InUnit,'S',iD,rD,TheName,iD,[iD],[rD]).ne.0) goto 666
c unfortunately MyGetKey uppercases strings!
      endif
      if(iKey.eq.16) then
* VB
c      isVB=1
      call Quit_OnUserError()
      endif
      if(iKey.eq.17) then
* All
      isAll=1
      endif

      if(iKey.eq.18) then
* Atom
      isAtom=1
      isNet=-1
      iGridNpt(1)=0
      iGridNpt(2)=0
      iGridNpt(3)=0
      endif
      if(iKey.eq.19) then
* CUBE
      iGauss=1
      isBinary=0
      write(6,*) 'Cube option is moved to grid2cube'
      call Quit_OnUserError()
      endif
* Grid
      if(iKey.eq.20) then
      isUserGrid=1
      isBinary=0
      isNet=-1
      iGridNpt(1)=0
      iGridNpt(2)=0
      iGridNpt(3)=0
      if(MyGetKey(InUnit,'I',nGridPoints,rD,Key,iD,[iD],[rD]).ne.0)
     *      goto 666
        Call GetMem('Grid','ALLO','REAL',ipGrid,nGridPoints*3)
       Read(InUnit,*,Err=666, end=666)
     *      (Work(ipGrid+i-1),i=1,nGridPoints*3)
      endif
* Pack
      if (iKey .eq. 21) then
c        imoPack=1
      endif
* PkLims
      if (iKey .eq. 22) then
        if (MyGetKey(InUnit,'D',iD,rD,Key,4,[iD],dTemp).ne.0) goto 666
        xLeft=dTemp(1)
        xRight=dTemp(2)
        xLoErr=dTemp(3)
        xHiErr=dTemp(4)
      endif
* PkBits
      if (iKey .eq. 23) then
        if (MyGetKey(InUnit,'I',ibits,rD,Key,iD,[iD],[rD]).ne.0)
     *      goto 666
        if (ibits .eq. 16) then
          iyRange=32768
          nBytesPackedVal=2
        endif
      endif
* NoOrbitals
      if (iKey .eq. 24) then
        NoOrb=1
      endif
* LINE - density on line
      if (iKey .eq. 25) then
      if(MyGetKey(InUnit,'D',iD,rD,Key,7,[iD],OneCoor).ne.0) goto 666
      isTheOne=1
      isTotal=1
      isBinary=0
      isLine=1
      endif
      if(iKey.eq.26) then
* ORANGE
      if(MyGetKey(InUnit,'D',iD,rD,Key,2,[iD],Region).ne.0) goto 666
      itRange=0
      isAuMO=1
      NoSort=1
      endif
      if(iKey.eq.27) then
* ERANGE
      if(MyGetKey(InUnit,'D',iD,rD,Key,2,[iD],Region).ne.0) goto 666
      itRange=1
      isAuMO=1
      endif
* DEBUG
      if(iKey.eq.28) then
      isBinary=0
      isDebug=1
      endif
      if(iKey.eq.29) then
* CUTOFF
      if(MyGetKey(InUnit,'R',iD,CutOff,Key,iD,[iD],[rD]).ne.0) goto 666
      isCutOff=1
      endif
      if(iKey.eq.30) then
* NOPACK
c      imoPack=0
      endif
      if(iKey.eq.31) then
* GORI
      iCustOrig=1
      if(MyGetKey(InUnit,'D',iD,rD,Key,3,[iD],GridOrigin).ne.0) goto 666
      if(MyGetKey(InUnit,'D',iD,rD,Key,3,[iD],GridAxis1).ne.0) goto 666
      if(MyGetKey(InUnit,'D',iD,rD,Key,3,[iD],GridAxis2).ne.0) goto 666
      if(MyGetKey(InUnit,'D',iD,rD,Key,3,[iD],GridAxis3).ne.0) goto 666
      endif
      if(iKey.eq.32) then
* SELEct
      if(MyGetKey(InUnit,'S',iD,rD,SelectStr,iD,[iD],[rD]).ne.0) then
         goto 666
      endif
      if(nReq.gt.0) then
      write(6,*) 'SELEct keyword can not be used together with ORBItals'
        call Quit_OnUserError()
      endif
      call gridExpandSelect(SelectStr)
       isAuMO=0
      endif
      if(iKey.eq.33) then
* NOSOrt
       NoSort=1
      endif
      if(iKey.eq.34) then
* FILE
      read(InUnit,'(A)') FileIn
      isFileOrb=1
      call fileorb(FileIn,FileStr)
      write(6,*) 'INPORB file: ',FileStr(:mylen(FileStr))
      endif
      if(iKey.eq.35) then
* SPHR
      isSphere=1
      endif
      if(iKey.eq.36) then
* COLOr
      isColor=1
      endif
      if(iKey.eq.37) then
* VIRT
      isVirt=1
      if(MyGetKey(InUnit,'R',iD,Virt,Key,iD,[iD],[rD]).ne.0) goto 666
      endif
      if(iKey.eq.38) then
* MULLiken charges per MO
      isMULL=1
      read(InUnit,'(A)') MULLPRT
      Call upCASE(MULLPRT)
      Call LeftAd(MULLPRT)
      If (MULLPRT(1:4).eq.'LONG') isLONGPRT=1
      endif
* SUBBLOCK
      if(iKey.eq.39) then
      if(MyGetKey(InUnit,'D',iD,rD,Key,3,[iD],SubBlock).ne.0) goto 666
      if(MyGetKey(InUnit,'R',iD,rSubBlock,Key,iD,[iD],[rD]).ne.0)
     *      goto 666
      iSubBlock=1
      endif
* XDER,YDER,ZDER,GDER
      if(iKey.eq.40) then
         isDerivative=1
         call Quit_OnUserError()
         endif
      if(iKey.eq.41) then
         isDerivative=2
         call Quit_OnUserError()
         endif
      if(iKey.eq.42) then
         isDerivative=3
         call Quit_OnUserError()
         endif
      if(iKey.eq.43) then
         isDerivative=4
         call Quit_OnUserError()
         endif

* CURD (current density)
      if(iKey.eq.44) then
         isCurDens=1
         call Quit_OnUserError()
         endif

* CRXJ (current density, rxj)
      if(iKey.eq.45) then
        isCurDens=1
        isRxJ=1
        call Quit_OnUserError()
      endif
* UMAX (use magnetic axes)
      if(iKey.eq.46) iuseMaxes=1
* NOLUSCUS
      if(iKey.eq.47) then
        isLuscus=0
        isBinary=0
      endif
* XFIEld - ask Grid_It to compute electronic density on a DFT integration grid
      if(iKey.eq.48) then
        isXField=1
        isReadNet=isReadNet+1 !make the grid definition exclusive
      endif
* LUS1
      if(iKey.eq.49) then
      write(6,*) 'Not implemented'
c      isLusMath=1
c      read(InUnit,'(a)') LUS1
      endif
* LUS2
      if(iKey.eq.50) then
      write(6,*) 'Not implemented'
c      isLusMath=1
c      read(InUnit,'(a)') LUS2
      endif
* PLUS
      if(iKey.eq.51) then
      write(6,*) 'Not implemented'
c      aLusMath=1
      endif
* MINUS
      if(iKey.eq.52) then
      write(6,*) 'Not implemented'
c      aLusMath=-1
      endif
* XFMI xfield minimum charge of each grid point to be stored
      if(iKey.eq.53) then
         if(MyGetKey(InUnit,'R',iD,XFminCh,Key,iD,[iD],[rD]).ne.0)
     *     goto 666
      endif
      goto 998

 666  write(6,'(a,a,a)')
     *     'Error during reading ', Key(1:20), 'section in input file'
      Call Quit_OnUserError()

************************************************************************
*                                                                      *
*                       End of input section.                          *
*                                                                      *
************************************************************************
 997  Continue
c      if(isLusMath.eq.1) return
      if(isLuscus.eq.1.and.isBinary.eq.0) then
        write(6,*) 'ASCII keyword is set, but NoLUSCUS is not'
        write(6,*) 'calling abend as the best option available'
        call Quit_OnUserError()
      endif
      close(InUnit)
      IF (isLuscus .EQ. 1) THEN
        IF (isLine .NE. 0) THEN
          WRITE(6,*) 'LUSCUS and LINE options are not compatible'
          Call Quit_OnUserError()
        END IF
      END IF
*
      if(isReadNet.GT.1) write(6,'(a)')
     *    'Warning: Double definition of GRID net'

*
*  Well, there is something to do!
*
       if(isFileOrb.eq.1) then
       INPORB=FileStr
       endif
       Call OpenGrid(INPORB)
*
*  try to generate grid position automatically
*

          magicValue=0
          if(isNet.ge.0) then
              magicValue=GridNormal
              if(isNet.eq.1) magicValue=GridSparse
              if(isNet.eq.2) magicValue=GridDense
          endif
          if(iCustOrig.eq.1) then
            if(isNet.ne.-1) then
              write(6,*) 'GORI can be used only with NPOI'
              call Quit_OnUserError()
            endif
          endif
          Call MyCoor(isAuto,
     &       GridOrigin(1),GridOrigin(2),GridOrigin(3),
     &       GridAxis1(1),GridAxis2(2),GridAxis3(3),
     &       iGridNpt(1),iGridNpt(2),iGridNpt(3), magicValue,
     &       iCustOrig)

          if(iSubBlock.eq.1) then
          do i=1,3
          GridOrigin(i)=SubBlock(i)-rSubBlock
          enddo
          GridAxis1(1)=rSubBlock*2
          GridAxis2(2)=rSubBlock*2
          GridAxis3(3)=rSubBlock*2
          iGridNpt(1)=40
          iGridNpt(2)=40
          iGridNpt(3)=40
          endif


      Return
* Avoid unused argument warnings
      If (.False.) Call Unused_integer(iReturn)
      End
      Subroutine gridExpandSelect(SelectStr)
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "grid.fh"
      Character SelectStr*120, tmp*120
       ifirst=0
       ilen=0
       tmp=' '
       do i=1,120
         if(ifirst.eq.0.and.SelectStr(i:i).ne.' ') then
           ilen=ilen+1
           tmp(ilen:ilen)=SelectStr(i:i)
           ifirst=1
           goto 10
         endif
         if(ifirst.eq.1.and.SelectStr(i:i).ne.' ') then
           ilen=ilen+1
           tmp(ilen:ilen)=SelectStr(i:i)
         endif
         if(ifirst.eq.1.and.SelectStr(i:i).eq.' ') then
           ilen=ilen+1
           tmp(ilen:ilen)=SelectStr(i:i)
           ifirst=0
         endif
10     continue
       enddo
      if(ilen.lt.2) then
         write(6,*) 'SELEct section is incomplete'
         call Quit_OnUserError()
      endif
c      print *,'current',tmp
      nReq=0
1     istart=1
      iend=index(tmp(istart:),' ')
       ibr=index(tmp(istart:iend),':')
       if(ibr.eq.0) then
         write(6,*) 'Wrong format in SELEct section'
         write(6,*) 'Expecting : sign in >',tmp(istart:iend),'<'
         call Quit_OnUserError()
       endif
c       print *,'v01 >',tmp(istart:istart+ibr-2),'<'
       read(tmp(istart:istart+ibr-2),*,err=20,end=20) isymm
       ibrm=index(tmp(istart+ibr+1:iend),'-')
       if(ibrm.eq.0) then
c the only number
c       print *,'v02 >',tmp(istart+ibr:iend),'<'
         read(tmp(istart+ibr:iend),*,err=20,end=20) iibeg
         iReq(2*nReq+1)=isymm
         iReq(2*nReq+2)=iibeg
         nReq=nReq+1
         if(nReq.gt.MAXGRID) goto 30
       else
c         print *,'v03 >',tmp(istart+ibr:istart+ibr+ibrm-1),'<'
         read(tmp(istart+ibr:istart+ibr+ibrm-1),*,err=20,end=20) iibeg
c        print *,'v04 >',tmp(istart+ibr+ibrm+1:iend),'<'
         read(tmp(istart+ibr+ibrm+1:iend),*,err=20,end=20) iiend
         if(iiend.lt.iibeg) then
          write(6,*) 'Wrong data in SELEct section'
          call Quit_OnUserError()
         endif
        do i=iibeg,iiend
         iReq(2*nReq+1)=isymm
         iReq(2*nReq+2)=i
         nReq=nReq+1
         if(nReq.gt.MAXGRID) goto 30
        enddo
       endif
c      print *,'current',tmp(iend+1:)
      tmp=tmp(iend+1:)
      if(tmp.ne.' ') goto 1
      Return
20    write(6,*) 'Error in analyzing SELECT section'
      call Quit_OnUserError()
30    write(6,*) 'Too many Grids requested'
      call Quit_OnUserError()
      end

      integer function LineNr(Lu)
      implicit none
      integer Lu, ios
      real*8 foo
      LineNr=0
      do while(.true.)
      read(Lu,*,iostat=ios) foo, foo, foo, foo
        if( ios > 0 ) then
           call Abend() !'problem somewhere'
        else if( ios < 0 ) then ! end of file is reached
           return
        else
           LineNr = LineNr + 1
        end if
      end do
      return
#ifdef _WARNING_WORKAROUND_
      if (.false.) call Unused_real(foo)
#endif
      end function LineNr
