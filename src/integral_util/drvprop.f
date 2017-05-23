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
* Copyright (C) 1991, Roland Lindh                                     *
*               2001, Hans-Joachim Werner                              *
************************************************************************
      SubRoutine Drvprop(Opname,Ccoor,opmol,opnuc,ncmp,idirect,isyop,
     >                   ptchrg,ngrid,iaddpot)
************************************************************************
*                                                                      *
* Object: driver for computation of one-electron property matrices     *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              OneEl                                                   *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January '91                                              *
*                                                                      *
*     Modified for Properties only by HJW AUg 2001                     *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
c     External Get_ProgName
      Character*(*) Opname
      External MltInt, KnEInt, MVeInt, VeInt,  D1Int,  NAInt,  EFInt,
     &         OAMInt, OMQInt, DMSInt, Potint,
     &         AMPInt
      External MltMem, KnEMem, MVeMem, VeMem,  D1Mem,  NAMem,  EFMem,
     &         OAMMem, OMQMem, DMSMem,
     &         AMPMem
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "wldata.fh"
#include "property_label.fh"
      Character*8 Label
c     Character*100 ProgName, Get_ProgName
      Real*8 Ccoor(3),opmol(*),opnuc(*),ptchrg(*)
      Real*8 Dum(1)
      Real*8, Allocatable :: Chrg(:), Centr(:,:)
cnf
      Logical Do_ESPF
cnf
*                                                                      *
************************************************************************
*                                                                      *
c     ProgName=Get_ProgName()
c     Call UpCase(ProgName)
*                                                                      *
************************************************************************
*                                                                      *
*
      PLabel=' '
      If(opname.eq.'POT'.or.iaddpot.lt.0) then
        call inisewm('mltpl',0)
      else
        call inisewm('mltpl',3)
      end if
*
      isymxy=1
      isymxz=1
      isymyz=1
      isyxyz=1
      Call qEnter('DrvProp')
      ipad=1
cnf
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)
      ntdg = 0
      Do iIrrep = 0, nIrrep - 1
         ntdg = ntdg + nBas(iIrrep)*(nBas(iIrrep)+1)/2
      End Do
      Call DecideOnESPF(Do_ESPF)
cnf
c
      Call mma_allocate(Centr,3,mCentr)
      Call mma_allocate(Chrg,mCentr)
      ndc = 0
      nc = 1
      Do jCnttp = 1, nCnttp
         mCnt = nCntr(jCnttp)
         If (AuxCnttp(jCnttp)) mCnt = 0
         jxyz = ipCntr(jCnttp)
         Do jCnt = 1, mCnt
            ndc = ndc + 1
            x1 = Work(jxyz)
            y1 = Work(jxyz+1)
            z1 = Work(jxyz+2)
            Do i = 0, nIrrep/nStab(ndc) - 1
               iFacx=iPhase(1,iCoset(i,0,ndc))
               iFacy=iPhase(2,iCoset(i,0,ndc))
               iFacz=iPhase(3,iCoset(i,0,ndc))
               Centr(1,nc) = x1*DBLE(iFacx)
               Centr(2,nc) = y1*DBLE(iFacy)
               Centr(3,nc) = z1*DBLE(iFacz)
               nchr=iAtmNr(jCnttp)
               Chrg(nc) = DBLE(nchr)
               nc = nc + 1
            End Do
            jxyz = jxyz + 3
         End Do
      End Do
      nc = nc-1
c
      idone=0

************************************************************************
*                                                                      *
*     Multipole Moments starting with the overlap. If SEWARD is run in *
*     the property mode we will skip the overlap integrals.            *
*                                                                      *
************************************************************************

      rHrmt=One
      imltpl=-1
      if(Opname.eq.'OV') imltpl=0
      if(Opname.eq.'DM') imltpl=1
      if(Opname.eq.'SM') imltpl=2
      if(Opname.eq.'TM') imltpl=3
      if(imltpl.ge.0) then
         idone=1
         Write (Label,'(A,I2)') 'Mltpl ', iMltpl
         nComp = (iMltpl+1)*(iMltpl+2)/2
         Call GetMem('ip    ','ALLO','INTE',ip1,nComp)
         Call GetMem('lOper ','ALLO','INTE',ip2,nComp)
         Call GetMem('kOper ','ALLO','INTE',ip3,nComp)
         Call GetMem('Nuc   ','ALLO','REAL',ipNuc,nComp)
         Call GetMem('Ccoor' ,'ALLO','REAL',ipCc,nComp*3)
         iComp=0
         Do ix = iMltpl, 0, -1
            If (Mod(ix,2).eq.0) Then
               iSymX=1
            Else
               ixyz=1
               iSymX=2**IrrFnc(ixyz)
               If (Ccoor(1).ne.Zero) iSymX = iOr(iSymX,1)
            End If
            Do iy = iMltpl-ix, 0, -1
               If (Mod(iy,2).eq.0) Then
                  iSymY=1
               Else
                  ixyz=2
                  iSymY=2**IrrFnc(ixyz)
                  If (Ccoor(2).ne.Zero) iSymY = iOr(iSymY,1)
               End If
               iz = iMltpl-ix-iy
               If (Mod(iz,2).eq.0) Then
                  iSymZ=1
               Else
                  ixyz=4
                  iSymZ=2**IrrFnc(ixyz)
                  If (Ccoor(3).ne.Zero) iSymZ = iOr(iSymZ,1)
               End If
               iChO = Mod(ix,2)*iChBas(2)
     &              + Mod(iy,2)*iChBas(3)
     &              + Mod(iz,2)*iChBas(4)
*
               iWork(ip2+iComp) = MltLbl(iSymX,MltLbl(iSymY,iSymZ,
     &                            nIrrep),nIrrep)
               iWork(ip3+iComp) = iChO
               call dcopy_(3,Ccoor,1,Work(ipCc+iComp*3),1)
               iComp = iComp + 1
            End Do
         End Do
*
         Call MltNuc(Work(ipCc),Chrg,Centr,mCentr,
     &               Work(ipNuc),iMltpl,nComp)
         Call OneEl(MltInt,MltMem,Label,iWork(ip1),iWork(ip2),nComp,
     &              Work(ipCc),iMltpl,Work(ipNuc),rHrmt,iWork(ip3),
     &              opmol,ipad,opnuc,iopadr,idirect,isyop,
     &              Dum,1,0)
*
         Call GetMem('Ccoor' ,'FREE','REAL',ipCc,nComp*3)
         Call GetMem('Nuc   ','FREE','REAL',ipNuc,nComp)
         Call GetMem('kOper ','FREE','INTE',ip3,nComp)
         Call GetMem('lOper ','FREE','INTE',ip2,nComp)
         Call GetMem('ip    ','FREE','INTE',ip1,nComp)
       End If
*
************************************************************************
*                                                                      *
*     Kinetic energy, nuclear attraction                               *
*                                                                      *
*     Mass-velocity and One-electron Darwin contact term integrals.    *
*                                                                      *
************************************************************************
      rHrmt=One
      nComp=1
      nOrdOp = 0
      Call GetMem('ip    ','ALLO','INTE',ip1,nComp)
      Call GetMem('lOper ','ALLO','INTE',ip2,nComp)
      Call GetMem('kOper ','ALLO','INTE',ip3,nComp)
      Call GetMem('Ccoor ','ALLO','REAL',ipCc,nComp*3)
      call dcopy_(3,Zero,0,Work(ipCc),1)
      iWork(ip2) = 1
      iWork(ip3) = iChBas(1)
*
      If(Opname.eq.'EKIN'.or.Opname.eq.'KINETIC') then
         idone=1
         Label='Kinetic '
         Call OneEl(KnEInt,KnEMem,Label,iWork(ip1),iWork(ip2),nComp,
     &              Work(ipCC),nOrdOp,Zero,rHrmt,iWork(ip3),
     &              opmol,ipad,opnuc,iopadr,idirect,isyop,
     &              Dum,1,0)
      End If
*
      If (Opname.eq.'ATTR') then
         idone=1
         Label='Attract '
         Call OneEl(NAInt,NAMem,Label,iWork(ip1),iWork(ip2),nComp,
     &              Work(ipCC),nOrdOp,PotNuc,rHrmt,iWork(ip3),
     &              opmol,ipad,opnuc,iopadr,idirect,isyop,
     &              Dum,1,0)
      End If
*
      If (Opname.eq.'POT') then
         idone=1
         Label='Pot '
         If (iaddpot.le.0.and..not.Do_ESPF) Then
cnf         If (iaddpot.le.0.and.ProgName.ne.'ALASKA') Then
            Call GetMem('Nuc ','ALLO','REAL',ipNuc,ngrid)
            Call Pot_nuc(CCoor,work(ipnuc),ngrid)
         Else
           ipnuc=ip_Dummy
         End if
         If (iaddpot.lt.0) Then
            If (iaddpot.eq.-1) Then
               Call Get_D1ao_Var(ipdens,Length)
            Else
               Call Get_D1ao(ipdens,Length)
            End If
            call Drv1_Pot(work(ipdens),CCoor,ptchrg,ngrid,1,0)
            Call GetMem('DENS','FREE','REAL',ipdens,ntdg)
            If (.not.Do_ESPF) Then
cnf         If (ProgName.ne.'ALASKA') Then
               Call AddVec(ptchrg,ptchrg,work(ipnuc),ngrid)
               Call fMove(work(ipnuc),opnuc,ngrid)
            End If
         Else
           iWork(ip2) = 2**nirrep-1
           isyoper=1
           Call OneEl(PotInt,NAMem,Label,iWork(ip1),iWork(ip2),ncmp,
     &                Ccoor,nOrdOp,work(ipnuc),rHrmt,iWork(ip3),
     &                opmol,ipad,opnuc,iopadr,idirect,isyoper,
     &                ptchrg,ngrid,iaddpot)
            If (iaddpot.eq.0.and..not.Do_ESPF)
c           If (iaddpot.eq.0.and.ProgName.ne.'ALASKA')
     &         opnuc(1)=work(ipnuc)
         End If
         If (iaddpot.le.0.and..not.Do_ESPF)
cnf      If (iaddpot.le.0.and.ProgName.ne.'ALASKA')
     &      Call GetMem('Nuc ','FREE','REAL',ipNuc,ngrid)
      End If
*
      If (Opname.eq.'REL'.or.Opname.eq.'MASSV') then
*
         idone=1
         Label='MassVel '
         Call OneEl(MVeInt,MVeMem,Label,iWork(ip1),iWork(ip2),nComp,
     &              Work(ipCc),nOrdOp,Zero,rHrmt,iWork(ip3),
     &              opmol,ipad,opnuc,iopadr,idirect,isyop,
     &              Dum,1,0)
      End If
*
      If (Opname.eq.'REL'.or.Opname.eq.'DARW') then
         idone=1
         Label='Darwin  '
         Call OneEl(D1Int,D1Mem,Label,iWork(ip1),iWork(ip2),nComp,
     &              Work(ipCc),nOrdOp,Zero,rHrmt,iWork(ip3),
     &              opmol,ipad,opnuc,iopadr,idirect,isyop,
     &              Dum,1,0)
      End If
*
      Call GetMem('Ccoor ','FREE','REAL',ipCc,nComp*3)
      Call GetMem('kOper ','FREE','INTE',ip3,nComp)
      Call GetMem('lOper ','FREE','INTE',ip2,nComp)
      Call GetMem('ip    ','FREE','INTE',ip1,nComp)
************************************************************************
*                                                                      *
*     Velocity integrals.                                              *
*                                                                      *
************************************************************************
      if(Opname.eq.'VELO') then
         idone=1
         rHrmt=-One
         nOrdOp = 1
         Label='Velocity'
         nComp = 3
         Call GetMem('ip    ','ALLO','INTE',ip1,nComp)
         Call GetMem('lOper ','ALLO','INTE',ip2,nComp)
         Call GetMem('kOper ','ALLO','INTE',ip3,nComp)
         Call GetMem('Nuc   ','ALLO','REAL',ipNuc,nComp)
         Call GetMem('Ccoor ','ALLO','REAL',ipCc,nComp*3)
         call dcopy_(3*nComp,Zero,0,Work(ipCc),1)
         ixyz=1
         iWork(ip2  ) = 2**IrrFnc(ixyz)
         iWork(ip3  ) = iChBas(2)
         ixyz=2
         iWork(ip2+1) = 2**IrrFnc(ixyz)
         iWork(ip3+1) = iChBas(3)
         ixyz=4
         iWork(ip2+2) = 2**IrrFnc(ixyz)
         iWork(ip3+2) = iChBas(4)
*
         call dcopy_(3,Zero,0,Work(ipNuc),1)
         Call OneEl(VeInt,VeMem,Label,iWork(ip1),iWork(ip2),nComp,
     &              Work(ipCc),nOrdOp,Work(ipNuc),rHrmt,iWork(ip3),
     &              opmol,ipad,opnuc,iopadr,idirect,isyop,
     &              Dum,1,0)
*
         Call GetMem('Ccoor ','FREE','REAL',ipCc,nComp*3)
         Call GetMem('Nuc   ','FREE','REAL',ipNuc,nComp)
         Call GetMem('kOper ','FREE','INTE',ip3,nComp)
         Call GetMem('lOper ','FREE','INTE',ip2,nComp)
         Call GetMem('ip    ','FREE','INTE',ip1,nComp)
      End If
************************************************************************
*                                                                      *
*     Electric field integrals.                                        *
*                                                                      *
************************************************************************
      If(Opname.eq.'EF') then
        idone=1
        iEF=1
        rHrmt=One
        nOrdOp = 1
        nComp = 3
        if(iaddpot.lt.0) then
         Call GetMem('DENS','ALLO','REAL',ipdens,ntdg)
         call Drv1_Pot(work(ipdens),CCoor,ptchrg,ngrid,ncomp,nOrdOp)
         Call GetMem('DENS','FREE','REAL',ipdens,ntdg)
         call addvec(ptchrg,ptchrg,opnuc,ngrid*ncomp)
        else
         ixyz=1
         iSymX = 2**IrrFnc(ixyz)
         ixyz=2
         iSymY = 2**IrrFnc(ixyz)
         ixyz=4
         iSymZ = 2**IrrFnc(ixyz)
         ixyz=3
         iSymXY = 2**IrrFnc(ixyz)
         ixyz=5
         iSymXZ = 2**IrrFnc(ixyz)
         ixyz=6
         iSymYZ = 2**IrrFnc(ixyz)
         ixyz=7
         iSyXYZ = 2**IrrFnc(ixyz)
         Write (Label,'(A,I5)') 'EF ',iEF
         Call GetMem('ip    ','ALLO','INTE',ip1,nComp)
         Call GetMem('lOper ','ALLO','INTE',ip2,nComp)
         Call GetMem('kOper ','ALLO','INTE',ip3,nComp)
         Call GetMem('Nuc   ','ALLO','REAL',ipNuc,nComp)
         Call GetMem('Ccoor ','ALLO','REAL',ipCc,nComp*3)
         iSymC = 1
         If (Ccoor(1).ne.Zero) iSymC = iOr(iSymC,iSymX)
         If (Ccoor(2).ne.Zero) iSymC = iOr(iSymC,iSymY)
         If (Ccoor(3).ne.Zero) iSymC = iOr(iSymC,iSymZ)
         If (Ccoor(1).ne.Zero .and. Ccoor(2).ne.Zero)
     &      iSymC = iOr(iSymC,iSymXY)
         If (Ccoor(1).ne.Zero .and. Ccoor(3).ne.Zero)
     &      iSymC = iOr(iSymC,iSymXZ)
         If (Ccoor(2).ne.Zero .and. Ccoor(3).ne.Zero)
     &      iSymC = iOr(iSymC,iSymYZ)
         If (Ccoor(1).ne.Zero .and. Ccoor(2).ne.Zero .and.
     &       Ccoor(3).ne.Zero) iSymC = iOr(iSymC,iSyXYZ)
*
         iComp=0
         Do ix = nOrdOp, 0, -1
            Do iy = nOrdOp-ix, 0, -1
               iComp = iComp + 1
               iz = nOrdOp-ix-iy
               ixyz=0
               If (Mod(ix,2).ne.0) ixyz=iOr(ixyz,1)
               If (Mod(iy,2).ne.0) ixyz=iOr(ixyz,2)
               If (Mod(iz,2).ne.0) ixyz=iOr(ixyz,4)
               iSym = 2**IrrFnc(ixyz)
               If (Ccoor(iComp).ne.Zero ) iSym = iOr(iSym,1)
               iWork(ip2+(iComp-1)) = MltLbl(iSymC,iSym,nIrrep)
               iWork(ip3+(iComp-1)) = iChBas(iComp+1)
               call dcopy_(3,Ccoor,1,Work(ipCc+(iComp-1)*3),1)
            End Do
         End Do
*
         Call EFNuc(Work(ipCc),Chrg,Centr,mCentr,
     &              Work(ipNuc),nOrdOp)
         Call OneEl(EFInt,EFMem,Label,iWork(ip1),iWork(ip2),nComp,
     &              Work(ipCc),nOrdOp,Work(ipNuc),rHrmt,iWork(ip3),
     &              opmol,ipad,opnuc,iopadr,idirect,isyop,
     &              Dum,1,0)
*
         Call GetMem('Ccoor ','FREE','REAL',ipCc,nComp*3)
         Call GetMem('Nuc   ','FREE','REAL',ipNuc,nComp)
         Call GetMem('kOper ','FREE','INTE',ip3,nComp)
         Call GetMem('lOper ','FREE','INTE',ip2,nComp)
         Call GetMem('ip    ','FREE','INTE',ip1,nComp)
       End If
      End If
************************************************************************
*                                                                      *
*     Electric field gradient integrals.                               *
*                                                                      *
************************************************************************
      If(Opname.eq.'FG'.or.Opname.eq.'EFG') then
        idone=1
        iEFg=1
        rHrmt=One
        nComp = 6
        nOrdOp = 2
        if(iaddpot.lt.0) then
         Call GetMem('DENS','ALLO','REAL',ipdens,ntdg)
         call Drv1_Pot(work(ipdens),CCoor,ptchrg,ngrid,ncomp,nOrdOp)
         Call GetMem('DENS','FREE','REAL',ipdens,ntdg)
         call addvec(ptchrg,ptchrg,opnuc,ngrid*ncomp)
        else
         ixyz=1
         iSymX = 2**IrrFnc(ixyz)
         ixyz=2
         iSymY = 2**IrrFnc(ixyz)
         ixyz=4
         iSymZ = 2**IrrFnc(ixyz)
         ixyz=3
         iSymXY = 2**IrrFnc(ixyz)
         ixyz=5
         iSymXZ = 2**IrrFnc(ixyz)
         ixyz=6
         iSymYZ = 2**IrrFnc(ixyz)
         ixyz=7
         iSyXYZ = 2**IrrFnc(ixyz)
         Write (Label,'(A,I5)') 'EFg',iEFg
         Call GetMem('ip    ','ALLO','INTE',ip1,nComp)
         Call GetMem('lOper ','ALLO','INTE',ip2,nComp)
         Call GetMem('kOper ','ALLO','INTE',ip3,nComp)
         Call GetMem('Nuc   ','ALLO','REAL',ipNuc,nComp)
         Call GetMem('Ccoor ','ALLO','REAL',ipCc,nComp*3)
         iSymC = 1
         If (Ccoor(1).ne.Zero) iSymC = iOr(iSymC,iSymX)
         If (Ccoor(2).ne.Zero) iSymC = iOr(iSymC,iSymY)
         If (Ccoor(3).ne.Zero) iSymC = iOr(iSymC,iSymZ)
         If (Ccoor(1).ne.Zero .and. Ccoor(2).ne.Zero)
     &      iSymC = iOr(iSymC,iSymXY)
         If (Ccoor(1).ne.Zero .and. Ccoor(3).ne.Zero)
     &      iSymC = iOr(iSymC,iSymXZ)
         If (Ccoor(2).ne.Zero .and. Ccoor(3).ne.Zero)
     &      iSymC = iOr(iSymC,iSymYZ)
         If (Ccoor(1).ne.Zero .and. Ccoor(2).ne.Zero .and.
     &       Ccoor(3).ne.Zero) iSymC = iOr(iSymC,iSyXYZ)
*
         ijComp = 0
         iComp=0
         Do 2010 ix = 1, 0, -1
         Do 2010 iy = 1-ix, 0, -1
            iComp = iComp + 1
            iz = 1-ix-iy
            ixyz=0
            If (Mod(ix,2).ne.0) ixyz=iOr(ixyz,1)
            If (Mod(iy,2).ne.0) ixyz=iOr(ixyz,2)
            If (Mod(iz,2).ne.0) ixyz=iOr(ixyz,4)
            iSym=2**IrrFnc(ixyz)

            If (Ccoor(iComp).ne.Zero .and. iChBas(iComp+1).ne.0 )
     &         iSym = iOr(iSym,1)
            jComp=iComp-1
            Do 2011 jx = ix, 0, -1
            jyMax=1-jx
            If (jx.eq.ix) jyMax=iy
            Do 2011 jy = jyMax, 0, -1
               call dcopy_(3,Ccoor,1,Work(ipCc+ijComp*3),1)
               ijComp = ijComp + 1
               jComp=jComp+1
               jz = 1-jx-jy
               jxyz=0
               If (Mod(jx,2).ne.0) jxyz=iOr(jxyz,1)
               If (Mod(jy,2).ne.0) jxyz=iOr(jxyz,2)
               If (Mod(jz,2).ne.0) jxyz=iOr(jxyz,4)
               jSym=2**IrrFnc(jxyz)
               If (Ccoor(jComp).ne.Zero .and. iChBas(jComp+1).ne.0 )
     &            jSym = iOr(jSym,1)
               iWork(ip2+(ijComp-1)) = MltLbl(iSymC,MltLbl(iSym,jSym,
     &                                 nIrrep),nIrrep)
               iWork(ip3+(ijComp-1)) = Mod(ix+jx,2)*iChBas(2)
     &                               + Mod(iy+jy,2)*iChBas(3)
     &                               + Mod(iz+jz,2)*iChBas(4)
*              Write (*,*) iSymC,iSym,jSym
 2011       Continue
 2010    Continue
*
         Call EFNuc(Work(ipCc),Chrg,Centr,mCentr,
     &               Work(ipNuc),nOrdOp)
         Call OneEl(EFInt,EFMem,Label,iWork(ip1),iWork(ip2),nComp,
     &              Work(ipCc),nOrdOp,Work(ipNuc),rHrmt,iWork(ip3),
     &              opmol,ipad,opnuc,iopadr,idirect,isyop,
     &              Dum,1,0)
*
         Call GetMem('Ccoor ','FREE','REAL',ipCc,nComp*3)
         Call GetMem('Nuc   ','FREE','REAL',ipNuc,nComp)
         Call GetMem('kOper ','FREE','INTE',ip3,nComp)
         Call GetMem('lOper ','FREE','INTE',ip2,nComp)
         Call GetMem('ip    ','FREE','INTE',ip1,nComp)
       End If
      End If
************************************************************************
*                                                                      *
*     Orbital angular momentum integrals.                              *
*                                                                      *
************************************************************************
      If (Opname.eq.'LOP'.or.Opname.eq.'ANGM') Then
         idone=1
         rHrmt=-One
         Label='AngMom  '
         nComp = 3
         nOrdOp = 1
         Call GetMem('ip    ','ALLO','INTE',ip1,nComp)
         Call GetMem('lOper ','ALLO','INTE',ip2,nComp)
         Call GetMem('kOper ','ALLO','INTE',ip3,nComp)
         Call GetMem('Nuc   ','ALLO','REAL',ipNuc,nComp)
         Call GetMem('Ccoor ','ALLO','REAL',ipCc,nComp*3)
         call dcopy_(3,Ccoor,1,Work(ipCc  ),1)
         call dcopy_(3,Ccoor,1,Work(ipCc+3),1)
         call dcopy_(3,Ccoor,1,Work(ipCc+6),1)
         ixyz=1
         iSymX = 2**IrrFnc(ixyz)
         ixyz=2
         iSymY = 2**IrrFnc(ixyz)
         ixyz=4
         iSymZ = 2**IrrFnc(ixyz)
         iSymCx = iSymX
         If (Ccoor(1).ne.Zero) iSymCx = iOr(iSymCx,1)
         iSymCy = iSymY
         If (Ccoor(2).ne.Zero) iSymCy = iOr(iSymCy,1)
         iSymCz = iSymZ
         If (Ccoor(3).ne.Zero) iSymCz = iOr(iSymCz,1)
*
         iSymLx = iOr(MltLbl(iSymCy,iSymZ,nIrrep),
     &                MltLbl(iSymCz,iSymY,nIrrep))
         iChOx = iChBas(3) + iChBas(4)
         iWork(ip2  ) = iSymLx
         iWork(ip3  ) = iChOx
         iSymLy = iOr(MltLbl(iSymCz,iSymX,nIrrep),
     &                MltLbl(iSymCx,iSymZ,nIrrep))
         iChOy = iChBas(4) + iChBas(2)
         iWork(ip2+1) = iSymLy
         iWork(ip3+1) = iChOy
         iSymLz = iOr(MltLbl(iSymCx,iSymY,nIrrep),
     &                MltLbl(iSymCy,iSymX,nIrrep))
         iChOz = iChBas(2) + iChBas(3)
         iWork(ip2+2) = iSymLz
         iWork(ip3+2) = iChOz
*
         call dcopy_(nComp,Zero,0,Work(ipNuc),1)
         Call OneEl(OAMInt,OAMMem,Label,iWork(ip1),iWork(ip2),nComp,
     &              Work(ipCc),nOrdOp,Work(ipNuc),rHrmt,iWork(ip3),
     &              opmol,ipad,opnuc,iopadr,idirect,isyop,
     &              Dum,1,0)
*
         Call GetMem('Ccoor ','FREE','REAL',ipCc,nComp*3)
         Call GetMem('Nuc   ','FREE','REAL',ipNuc,nComp)
         Call GetMem('kOper ','FREE','INTE',ip3,nComp)
         Call GetMem('lOper ','FREE','INTE',ip2,nComp)
         Call GetMem('ip    ','FREE','INTE',ip1,nComp)
      End If
************************************************************************
*                                                                      *
*     Orbital Magnetic Quadrupole integrals.                           *
*                                                                      *
************************************************************************
      If (Opname.eq.'OMQ') Then
         idone=1
         rHrmt=-One
         Label='OMQ     '
         nComp = 9
         nOrdOp = 3
         Call GetMem('ip    ','ALLO','INTE',ip1,nComp)
         Call GetMem('lOper ','ALLO','INTE',ip2,nComp)
         Call GetMem('kOper ','ALLO','INTE',ip3,nComp)
         Call GetMem('Nuc   ','ALLO','REAL',ipNuc,nComp)
         Call GetMem('Ccoor ','ALLO','REAL',ipCc,nComp*3)
*
         Call dcopy_(nComp,Work(ipOMQ),1,Work(ipCc  ),1) ! Change from 3 to ncomp?
         Call dcopy_(nComp,Work(ipOMQ),1,Work(ipCc+ncomp),1) !
         Call dcopy_(nComp,Work(ipOMQ),1,Work(ipCc+2*ncomp),1) !
         Call dcopy_(nComp,Work(ipOMQ),1,Ccoor,1) ! Should then not all be copied?
*
         ixyz=1
         iSymX = 2**IrrFnc(ixyz)
         ixyz=2
         iSymY = 2**IrrFnc(ixyz)
         ixyz=4
         iSymZ = 2**IrrFnc(ixyz)
         iSymCx = iSymX
         If (Ccoor(1).ne.Zero) iSymCx = iOr(iSymCx,1)
         iSymCy = iSymY
         If (Ccoor(2).ne.Zero) iSymCy = iOr(iSymCy,1)
         iSymCz = iSymZ
         If (Ccoor(3).ne.Zero) iSymCz = iOr(iSymCz,1)
*
         iSymLx = iOr(MltLbl(iSymCy,iSymZ,nIrrep),
     &                MltLbl(iSymCz,iSymY,nIrrep))
         iSymLy = iOr(MltLbl(iSymCz,iSymX,nIrrep),
     &                MltLbl(iSymCx,iSymZ,nIrrep))
         iSymLz = iOr(MltLbl(iSymCx,iSymY,nIrrep),
     &                MltLbl(iSymCy,iSymX,nIrrep))
*
* Calculates M_ij = r_j*L_i + L_j*r_i = r_j*L_i + r_i*L_j - i hbar E_ijk r_k
* Since the -i hbar is included outside we do
* M_ij = r_j*L_i + L_j*r_i = r_j*L_i + r_i*L_j + E_ijk r_k
*
* Mxx
         iSymxLx = MltLbl(iSymCx,iSymLx,nIrrep)
         iChOxx = iChBas(15)
         iWork(ip2  ) = iSymxLx
         iWork(ip3  ) = iChOxx
* Mxy
         iSymxLy = iSymCz
         iChOxy = iChBas(4)
         iWork(ip2+1) = iSymxLy
         iWork(ip3+1) = iChOxy
* Mxz
         iSymxLz = iSymCy
         iChOxz = iChBas(3)
         iWork(ip2+2) = iSymxLz
         iWork(ip3+2) = iChOxz
* Myx
         iSymyLx = iSymCz
         iChOyx = iChBas(4)
         iWork(ip2+3) = iSymyLx
         iWork(ip3+3) = iChOyx
* Myy
         iSymyLy = MltLbl(iSymCy,iSymLy,nIrrep)
         iChOyy = iChBas(15)
         iWork(ip2+4) = iSymyLy
         iWork(ip3+4) = iChOyy
* Myz
         iSymyLz = iSymCx
         iChOyz = iChBas(4)
         iWork(ip2+5) = iSymyLz
         iWork(ip3+5) = iChOyz
* Mzx
         iSymzLx = iSymCy
         iChOzx = iChBas(3)
         iWork(ip2+6) = iSymzLx
         iWork(ip3+6) = iChOzx
* Mzy
         iSymzLy = iSymCx
         iChOzy = iChBas(4)
         iWork(ip2+7) = iSymzLy
         iWork(ip3+7) = iChOzy
* Mzz
         iSymzLz = MltLbl(iSymCz,iSymLz,nIrrep)
         iChOzz = iChBas(15)
         iWork(ip2+8) = iSymzLz
         iWork(ip3+8) = iChOzz
*
         Call dcopy_(nComp,Zero,0,Work(ipNuc),1)
         Call OneEl(OMQInt,OMQMem,Label,iWork(ip1),iWork(ip2),nComp,
     &              Work(ipCc),nOrdOp,Work(ipNuc),rHrmt,iWork(ip3),
     &              opmol,ipad,opnuc,iopadr,idirect,isyop,
     &              Dum,1,0)
*
         Call GetMem('Ccoor ','FREE','REAL',ipCc,nComp*3)
         Call GetMem('Nuc   ','FREE','REAL',ipNuc,nComp)
         Call GetMem('kOper ','FREE','INTE',ip3,nComp)
         Call GetMem('lOper ','FREE','INTE',ip2,nComp)
         Call GetMem('ip    ','FREE','INTE',ip1,nComp)
      End If
************************************************************************
*                                                                      *
*     Diamagnetic shielding integrals.                                 *
*                                                                      *
************************************************************************
      If(Opname.eq.'DMS') then
         idone=1
         iDMS=1
         rHrmt=One
         Write (Label,'(A,I2,I2)') 'DMS ',1,iDMS
         nComp = 9
         Call GetMem('ip    ','ALLO','INTE',ip1,nComp)
         Call GetMem('lOper ','ALLO','INTE',ip2,nComp)
         Call GetMem('kOper ','ALLO','INTE',ip3,nComp)
         Call GetMem('Nuc   ','ALLO','REAL',ipNuc,nComp)
         Call GetMem('Ccoor ','ALLO','REAL',ipCc,nComp*3)
         iSymC = 1
         If (Ccoor(1).ne.Zero) iSymC = iOr(iSymC,iSymX)
         If (Ccoor(2).ne.Zero) iSymC = iOr(iSymC,iSymY)
         If (Ccoor(3).ne.Zero) iSymC = iOr(iSymC,iSymZ)
         If (Ccoor(1).ne.Zero .and. Ccoor(2).ne.Zero)
     &      iSymC = iOr(iSymC,iSymXY)
         If (Ccoor(1).ne.Zero .and. Ccoor(3).ne.Zero)
     &      iSymC = iOr(iSymC,iSymXZ)
         If (Ccoor(2).ne.Zero .and. Ccoor(3).ne.Zero)
     &      iSymC = iOr(iSymC,iSymYZ)
         If (Ccoor(1).ne.Zero .and. Ccoor(2).ne.Zero .and.
     &       Ccoor(3).ne.Zero) iSymC = iOr(iSymC,iSyXYZ)
*
         iComp = 0
         iC = 0
         Do 1510 ix = 1, 0, -1
         Do 1510 iy = 1-ix, 0, -1
            iz=1-ix-iy
            iC = iC + 1
            iChO1 = iChBas(iC+1)
            ixyz=0
            If (Mod(ix,2).ne.0) ixyz=iOr(ixyz,1)
            If (Mod(iy,2).ne.0) ixyz=iOr(ixyz,2)
            If (Mod(iz,2).ne.0) ixyz=iOr(ixyz,4)
            iSym = 2**IrrFnc(ixyz)
            If (Ccoor(iC).ne.Zero) iSym = iOr(iSym,1)
            iD = 0
            Do 1511 jx = 1, 0, -1
            Do 1511 jy = 1-jx, 0, -1
               jz=1-jx-jy
               iD = iD + 1
               iChO2 = iChBas(iD+1)
               jxyz=0
               If (Mod(jx,2).ne.0) jxyz=iOr(jxyz,1)
               If (Mod(jy,2).ne.0) jxyz=iOr(jxyz,2)
               If (Mod(jz,2).ne.0) jxyz=iOr(jxyz,4)
               iSymD = 2**IrrFnc(jxyz)
               If (Dxyz(iD).ne.Zero) iSymD = iOr(iSymD,1)
               If (iC.eq.iD) Then
                  i2 = iD + 1
                  If (i2.gt.3) i2 = i2 - 3
                  i3 = iD + 2
                  If (i3.gt.3) i3 = i3 - 3
                  iChO = iAnd(iChBas(i2+1),iChBas(i3+1))
               Else
                  iChO = iOr(iChO1,iChO2)
               End If
               iWork(ip2+iComp) = MltLbl(iSymD,MltLbl(iSym,iSymC,
     &                            nIrrep),nIrrep)
               iWork(ip3+iComp) = iChO
               call dcopy_(3,Ccoor,1,Work(ipCc+iComp*3),1)
               iComp = iComp + 1
 1511       Continue
 1510    Continue
         call dcopy_(3,Dxyz,1,Work(ipCc+3),1)
*
         call dcopy_(nComp,Zero,0,Work(ipNuc),1)
         Call OneEl(DMSInt,DMSMem,Label,iWork(ip1),iWork(ip2),nComp,
     &              Work(ipCc),1,Work(ipNuc),rHrmt,iWork(ip3),
     &              opmol,ipad,opnuc,iopadr,idirect,isyop,
     &              Dum,1,0)
*
         Call GetMem('Ccoor ','FREE','REAL',ipCc,nComp*3)
         Call GetMem('Nuc   ','FREE','REAL',ipNuc,nComp)
         Call GetMem('kOper ','FREE','INTE',ip3,nComp)
         Call GetMem('lOper ','FREE','INTE',ip2,nComp)
         Call GetMem('ip    ','FREE','INTE',ip1,nComp)
      End If
************************************************************************
*                                                                      *
*     Angular momentum products                                        *
*                                                                      *
************************************************************************
* Hermitized products of angular momentum integrals
* Component(1) is Lx*Lx
* Component(2) is (Lx*Ly+Ly*Lx)/2, etc.
* Coded P-A Malmqvist, Garching, Nov 1996
      If(Opname.eq.'LOP2'.or.Opname.eq.'PAM') then
         idone=1
         rHrmt=-One
         Label='AMProd  '
         nComp = 6
         nOrdOp = 2
         Call GetMem('ip    ','ALLO','INTE',ip1,nComp)
         Call GetMem('lOper ','ALLO','INTE',ip2,nComp)
         Call GetMem('kOper ','ALLO','INTE',ip3,nComp)
         Call GetMem('Nuc   ','ALLO','REAL',ipNuc,nComp)
         Call GetMem('Ccoor ','ALLO','REAL',ipCc,nComp*3)
         call dcopy_(nComp,Ccoor(1),0,Work(ipCc  ),3)
         call dcopy_(nComp,Ccoor(2),0,Work(ipCc+1),3)
         call dcopy_(nComp,Ccoor(3),0,Work(ipCc+2),3)
C Symmetry labels iSymX  for operator d/dx, etc.
C Symmetry labels iSymLx for operator Lx, etc.
C Characters iChOx for operator Lx, etc.
         ixyz=1
         iSymX = 2**IrrFnc(ixyz)
         ixyz=2
         iSymY = 2**IrrFnc(ixyz)
         ixyz=4
         iSymZ = 2**IrrFnc(ixyz)
         iSymCx = iSymX
         If (Ccoor(1).ne.Zero) iSymCx = iOr(iSymCx,1)
         iSymCy = iSymY
         If (Ccoor(2).ne.Zero) iSymCy = iOr(iSymCy,1)
         iSymCz = iSymZ
         If (Ccoor(3).ne.Zero) iSymCz = iOr(iSymCz,1)
         iSymLx = iOr(MltLbl(iSymCy,iSymZ,nIrrep),
     &                MltLbl(iSymCz,iSymY,nIrrep))
         iChOx = iChBas(3) + iChBas(4)
         iSymLy = iOr(MltLbl(iSymCz,iSymX,nIrrep),
     &                MltLbl(iSymCx,iSymZ,nIrrep))
         iChOy = iChBas(4) + iChBas(2)
         iSymLz = iOr(MltLbl(iSymCx,iSymY,nIrrep),
     &                MltLbl(iSymCy,iSymX,nIrrep))
         iChOz = iChBas(2) + iChBas(3)

C Symmetry labels and characters of products. Let G be the full
C  molecular point group, and Gsub=subgroup of G=stabilizer of
C gauge origin. The totally symmetric irrep of Gsub can be
C decomposed into irreps of G.
C Then symmetry label=packed array of bits, one for each irrep
C of G. The bit is set, if that irrep is included in the
C decomposition of the totally symmetric irrep of Gsub.
         iwork(ip2  )=MltLbl(iSymLx,iSymLx,nIrrep)
         iwork(ip2+1)=MltLbl(iSymLx,iSymLy,nIrrep)
         iwork(ip2+2)=MltLbl(iSymLx,iSymLz,nIrrep)
         iwork(ip2+3)=MltLbl(iSymLy,iSymLy,nIrrep)
         iwork(ip2+4)=MltLbl(iSymLy,iSymLz,nIrrep)
         iwork(ip2+5)=MltLbl(iSymLz,iSymLz,nIrrep)
         iwork(ip3  )=0
         iwork(ip3+1)=iEOr(iChOx,iChOy)
         iwork(ip3+2)=iEOr(iChOx,iChOz)
         iwork(ip3+3)=0
         iwork(ip3+4)=iEOr(iChOy,iChOz)
         iwork(ip3+5)=0
*
         call dcopy_(nComp,Zero,0,Work(ipNuc),1)
         Call OneEl(AMPInt,AMPMem,Label,iWork(ip1),iWork(ip2),nComp,
     &              Work(ipCc),nOrdOp,Work(ipNuc),rHrmt,iWork(ip3),
     &              opmol,ipad,opnuc,iopadr,idirect,isyop,
     &              Dum,1,0)
*
         Call GetMem('Ccoor ','FREE','REAL',ipCc,nComp*3)
         Call GetMem('Nuc   ','FREE','REAL',ipNuc,nComp)
         Call GetMem('kOper ','FREE','INTE',ip3,nComp)
         Call GetMem('lOper ','FREE','INTE',ip2,nComp)
         Call GetMem('ip    ','FREE','INTE',ip1,nComp)
      End If
*
      If(idone.eq.0) then
        Call WarningMessage(2,'Unknown operator in Prpdrv!')
        Write (6,*) 'Opname=',Opname
        Call Abend()
      end if
*
      Call mma_deallocate(Centr)
      Call mma_deallocate(Chrg)
      Call Free_iSD()
      Call qExit('DrvProp')
      Return
      End
