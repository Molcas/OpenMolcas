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

      SubRoutine Input_Grid_It(iRun,INPORB)
************************************************************************
*                                                                      *
* Object: input module for grid code                                   *
*                                                                      *
* Called from: grid                                                    *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*              MyCoor                                                  *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             September '91                                            *
*                                                                      *
*             Modified to complement GetInf, January '92.              *
*                                                                      *
* Copyright 1991 R. Lindh, Dept. of Theor. Chem. Univ. of Lund, Sweden *
* Copyright 1992 R. Lindh, Dept. of Theor. Chem. Univ. of Lund, Sweden *
*                                                                      *
* Patch 1999, Valera Veryazov                                          *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
      Include 'real.fh'
      Include 'WrkSpc.fh'
      Include 'grid.nosupport.fh'
      Character Key*80
      Character INPORB*(*)
      Integer iTemp(3)
      Dimension dTemp(4)
      Character AllKeys*190
      Character SelectStr*120
      Character FileStr*256, FileIn*256
      Character MULLprt*80
      Character Slash
      Slash='/'
*
      AllKeys=
     *'PRIN BINA ASCI NPOI DENS '//
     *'SPAR ORBI REGI ONE  TITL '//
     *'GAP  END  NODE TOTA NAME '//
     *'VB   ALL  ATOM CUBE GRID '//
     *'PACK PKLI PKBI NOOR LINE '//
     *'ORAN ERAN DEBU CUTO NOPA '//
     *'GORI SELE NOSO FILE SPHR '//
     *'COLO VIRT MULL '
*
      iRout = 99
      iPrint = nPrint(iRout)
#ifdef _DEBUG_
      Call qEnter('Input')
#endif
      Do 108 i = 1, nRout
         nPrint(i) = 5
 108  Continue
      isBinary=1
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
      isAuMO=-1
      isAtom=0
      isTotal=0
      isVB=0
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
      if(iRun.eq.0) then
c make defaults for a fast run via call from other module
        isNet=1  ! set sparse
c        isBinary=0 ! temporary set Ascii output
c       iMaxUp=1
c       iMaxDown=5
        imoPack=1
c        isCutOff=1
        goto 500
      endif

*
*     KeyWord directed input
*
      InUnit=5

      Call RdNLst(InUnit,'GRID_IT')
 998  if(MyGetKey_nosupport(InUnit,'S',iD,rD,Key,iD,iD,rD).ne.0) 
     * goto 997
      iKey=index(AllKeys,Key(1:4))

      if(iKey.eq.0.or.(iKey-1)/5*5.ne.(iKey-1)) then
      write(6,'(a,a)') 'Unrecognized keyword in input file:', Key(1:4)
        Call Quit_OnUserError()
      endif
      iKey=(iKey-1)/5+1
*
      if(iKey.eq.1) then
* PRIN
      if(mygetkey_nosupport(InUnit,'I',n,rD,Key,iD,iD,rD).ne.0) 
     * goto 666
      Do j = 1, n
      if(mygetkey_nosupport(InUnit,'A',iD,rD,Key,2,iTemp,rD).ne.0) 
     *goto 666
         nPrint(iTemp(1))=iTemp(2)
      enddo
      endif
      if(iKey.eq.2) then
* BINARY = default
        isBinary=1
      endif
      if(iKey.eq.3) then
* ASCII = for debug
c        Write(6,*) ' Keyword ASCII is obsolete'
        isBinary=0
      endif
      if(iKey.eq.4) then
* NPOI
      if(mygetkey_nosupport(InUnit,'A',iD,rD,Key,3,iGridNpt,rD).ne.0) 
     * goto 666
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
      if(mygetkey_nosupport(InUnit,'I',nReq,rD,Key,iD,iD,rD).ne.0) 
     * goto 666

      if(nReq.gt.MAXGRID) then
      write(6,'(a,i5,a,i5)')
     *      'Too many requested orbitals ',nReq,'>',MAXGRID
      Call Quit_OnUserError()
        endif
       read(inUnit,*, err=666, end=666) (iReq(i),i=1,nReq*2)
cc      if(mygetkey_nosupport(InUnit,'A',iD,rD,Key,nReq*2,iReq,rD).ne.0) goto 666
      isAuMO=0
      endif
      if(iKey.eq.8) then
* REGION
      if(mygetkey_nosupport(InUnit,'D',iD,rD,Key,2,iD,Region).ne.0) 
     * goto 666
      itRange=1
      isAuMO=1
      write(6,*) ' *** Warning keyword REGION is obsolete'
      write(6,*) ' ***         assumimg Energy range '
      endif
      if(iKey.eq.9) then
* ONE - debug option
      if(mygetkey_nosupport(InUnit,'D',iD,rD,Key,7,iD,OneCoor).ne.0) 
     * goto 666
      isTheOne=1
      isBinary=0
      endif
      if(iKey.eq.10) then
* TITLE
* NOTE: Title can be only ONE line here!!!
      if(mygetkey_nosupport(InUnit,'S',iD,rD,Title1,iD,iD,rD).ne.0) 
     * goto 666
      endif
      if(iKey.eq.11) then
* GAP
      if(mygetkey_nosupport(InUnit,'R',iD,TheGap,Key,iD,iD,rD).ne.0) 
     * goto 666
      endif

      if(iKey.eq.12) then
      goto 997
* END
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
c      if(mygetkey_nosupport(InUnit,'S',iD,rD,TheName,iD,iD,rD).ne.0) goto 666
c unfortunately mygetkey_nosupport uppercases strings!
      endif
      if(iKey.eq.16) then
* VB
      isVB=1
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
      endif
* Grid
      if(iKey.eq.20) then
      isUserGrid=1
      isBinary=0
      isNet=-1
      iGridNpt(1)=0
      iGridNpt(2)=0
      iGridNpt(3)=0
      if(mygetkey_nosupport(InUnit,'I',nGridPoints,rD,Key,1,iD,rD).ne.0) 
     * goto 666
        Call GetMem('Grid','ALLO','REAL',ipGrid,nGridPoints*3)
       Read(InUnit,*,Err=666, end=666)
     *      (Work(ipGrid+i-1),i=1,nGridPoints*3)
      endif
* Pack
      if (iKey .eq. 21) then
        imoPack=1
      endif
* PkLims
      if (iKey .eq. 22) then
        if (mygetkey_nosupport(InUnit,'D',iD,rD,Key,4,iD,
     *      dTemp).ne.0) goto 666
        xLeft=dTemp(1)
        xRight=dTemp(2)
        xLoErr=dTemp(3)
        xHiErr=dTemp(4)
      endif
* PkBits
      if (iKey .eq. 23) then
        if (mygetkey_nosupport(InUnit,'I',ibits,rD,Key,1,iD,rD).ne.0) 
     * goto 666
        if (ibits .eq. 16) then
          iyRange=32768
          nBytesPackedVal=2
        endif
      endif
* Pack
      if (iKey .eq. 24) then
        NoOrb=1
      endif
* LINE - density on line
      if (iKey .eq. 25) then
      if(mygetkey_nosupport(InUnit,'D',iD,rD,Key,7,iD,OneCoor).ne.0) 
     * goto 666
      isTheOne=1
      isTotal=1
      isBinary=0
      isLine=1
      endif
      if(iKey.eq.26) then
* ORANGE
      if(mygetkey_nosupport(InUnit,'D',iD,rD,Key,2,iD,Region).ne.0) 
     * goto 666
      itRange=0
      isAuMO=1
      NoSort=1
      endif
      if(iKey.eq.27) then
* ERANGE
      if(mygetkey_nosupport(InUnit,'D',iD,rD,Key,2,iD,Region).ne.0) 
     * goto 666
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
      if(mygetkey_nosupport(InUnit,'R',iD,CutOff,Key,1,iD,rD).ne.0) 
     * goto 666
      isCutOff=1
      endif
      if(iKey.eq.30) then
* NOPACK
      imoPack=0
      endif
      if(iKey.eq.31) then
* GORI
      iCustOrig=1
      if(mygetkey_nosupport(InUnit,'D',iD,rD,Key,3,iD,GridOrigin).ne.0) 
     * goto 666
      if(mygetkey_nosupport(InUnit,'D',iD,rD,Key,3,iD,GridAxis1).ne.0) 
     *goto 666
      if(mygetkey_nosupport(InUnit,'D',iD,rD,Key,3,iD,GridAxis2).ne.0) 
     *goto 666
      if(mygetkey_nosupport(InUnit,'D',iD,rD,Key,3,iD,GridAxis3).ne.0) 
     *goto 666
      endif
      if(iKey.eq.32) then
* SELEct
      if(mygetkey_nosupport(InUnit,'S',iD,rD,SelectStr,iD,iD,rD).ne.0) 
     * then
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
c      call fileorb(FileIn,FileStr)
      write(6,*) 'INPORB file: ',FileStr
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
      if(mygetkey_nosupport(InUnit,'R',iD,Virt,Key,1,iD,rD).ne.0) 
     * goto 666
      endif
      if(iKey.eq.38) then
* MULLiken charges per MO
      isMULL=1
      read(InUnit,'(A)') MULLPRT
      Call upCASE(MULLPRT)
      Call LeftAd(MULLPRT)
      If (MULLPRT(1:4).eq.'LONG') isLONGPRT=1
      endif

      goto 998
 666  write(6,'(a,a,a)')
     *     'Error during reading ', Key(1:20), 'section in input file'
      Call Quit_OnUserError()

*****************************************************************************
*                                                                           *
*                          End of input section.                            *
*                                                                           *
*****************************************************************************
 997  Continue
      close(InUnit)
      if(imoPack.eq.1.and.isBinary.eq.0) then
       Write(6,*) 'PACK and ASCII options are not compatible'
       Call Quit_OnUserError()
      endif
      if(imoPack.eq.1.and.isSphere.eq.1) then
       Write(6,*) 'PACK and SPHERE options are not compatible'
       Call Quit_OnUserError()
      endif
      if(imoPack.eq.1.and.isColor.eq.1) then
       Write(6,*) 'PACK and COLOR options are not compatible'
       Call Quit_OnUserError()
      endif
c if Pack/Nopack was not set - set it to Pack if not dense
      if(imoPack.eq.-1) then
       if(isBinary.eq.0) then
         imoPack=0
       else
         if(isNet.eq.2) then
            imoPack=0
         else
            imoPack=1
         endif
       endif
      endif
      if(imoPack.eq.-1) imoPack=0
*
      if(isReadNet.GT.1) write(6,'(a)')
     *    'Warning: Double definition of GRID net'

*
*  Well, there is something to do!
*
500    Continue
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
           call abend
         endif
      endif
      Call MyCoor_nosupport(iRun,isAuto,
     &      GridOrigin(1),GridOrigin(2),GridOrigin(3),
     &      GridAxis1(1),GridAxis2(2),GridAxis3(3),
     &      iGridNpt(1),iGridNpt(2),iGridNpt(3), magicValue,
     &      iCustOrig)

#ifdef _DEBUG_
      Call qExit('Input')
#endif
      Return
      End
      Subroutine gridExpandSelect(SelectStr)
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
      Include 'real.fh'
      include 'grid.nosupport.fh'
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
