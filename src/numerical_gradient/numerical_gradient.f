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
      Subroutine Numerical_Gradient(ireturn)
#ifndef _HAVE_EXTRA_
      Use Prgm
#endif
      Use Para_Info, Only: MyRank, nProcs, Set_Do_Parallel
#if defined (_MOLCAS_MPP_) && !defined(_GA_)
      Use Para_Info, Only: King
#endif
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "standard_iounits.fh"
#include "WrkSpc.fh"
#include "timtra.fh"
#include "real.fh"
#include "warnings.fh"
#include "constants2.fh"
#include "stdalloc.fh"
      Real*8 Energy_Ref
      Integer iOper(0:7), jStab(0:7), iCoSet(0:7,0:7),
     &        iDispXYZ(3)
      Character*8 Method
      Character AtomLbl(MxAtom)*(LENIN)
      Character Namei*(LENIN)
      Character*10 ESPFKey
      Character*180 Line,Get_Ln
      External Get_Ln
      External Rsv_Tsk
#if defined (_MOLCAS_MPP_) && !defined(_GA_)
      Character*80  SSTMNGR
      Integer  SSTMODE
      Logical  Rsv_Tsk_Even
      External Rsv_Tsk_Even
#endif
      Logical DispX, DispY, DispZ,Rsv_Tsk, Is_Roots_Set, Found,
     &        External_Coor_List, Do_ESPF, StandAlone, Exist, DoTinker,
     &        NMCart, DynExtPot, DoDirect, KeepOld, Reduce_Prt
      External Reduce_Prt
      Real*8 FX(3), rDum(1)
      Real*8, Allocatable, Dimension(:,:) :: EnergyArray, GradArray,
     &                                       OldGrads
      Real*8, Allocatable:: Grad(:), GNew(:), MMGrd(:,:)
      Integer, Allocatable:: IsMM(:)
      Real*8, Allocatable:: BMtrx(:,:), TMtrx(:,:), Coor(:,:),
     &                      Energies_Ref(:), XYZ(:,:), All(:,:),
     &                      Disp(:), Deg(:,:), Mltp(:), C(:,:),
     &                      Tmp2(:), Tmp(:,:)
      Integer rc, Read_Grad
      External Read_Grad
      Parameter (ToHartree = CONV_CAL_TO_J_ /
     &           CONV_AU_TO_KJ_PER_MOLE_)
*                                                                      *
************************************************************************
*                                                                      *
      Interface
      Subroutine RunTinker(nAtom,Cord,ipMltp,IsMM,MltOrd,DynExtPot,
     &                     iQMChg,nAtMM,StandAlone,DoDirect)
      Integer, Intent(In):: nAtom
      Real*8, Intent(In):: Cord(3,nAtom)
      Integer, Intent(In):: ipMltp
      Integer, Intent(In):: IsMM(nAtom)
      Integer, Intent(In):: MltOrd
      Logical, Intent(InOut):: DynExtPot
      Integer, Intent(In):: iQMChg
      Integer, Intent(InOut):: nAtMM
      Logical, Intent(In):: StandAlone
      Logical, Intent(In):: DoDirect
      End Subroutine RunTinker
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
*     Get the print level.
*
      iPL_Save=iPrintLevel(-1)
      iPL=iPL_Save

      If (Reduce_Prt().and.iPL.lt.3) iPL=0
*                                                                      *
************************************************************************
*                                                                      *
      ireturn=_RC_ALL_IS_WELL_
*
*     Get information regarding the last method used
*
      Call Get_cArray('Relax Method',Method,8)
      Call DecideOnESPF(Do_ESPF)
      Is_Roots_Set = .False.
      Call Qpg_iScalar('Number of roots',Is_Roots_Set)
      If (Is_Roots_Set) Then
         Call Get_iScalar('Number of roots',nRoots)
         If (nRoots.eq.1) Then
            iRoot=1
         Else
            Call qpg_iScalar('NumGradRoot',Found)
            If (Found) Then
               Call Get_iScalar('NumGradRoot',iRoot)
            Else
               iRoot=1
            End If
         End If
      Else
         nRoots = 1
         iRoot  = 1
      End If
      Call mma_allocate(Energies_Ref,nRoots,Label='Energies_Ref')
C     Print *,'Is_Roots_Set, nRoots, iRoot = ',Is_Roots_Set,nRoots,iRoot
      If (iRoot .gt. nRoots) Then
         Write(LuWr,*)
         Write(LuWr,*) '****************** ERROR ******************'
         Write(LuWr,*) 'It was selected to run numerical gradient'
         Write(LuWr,*) 'using energies for root/state number ',iRoot
         Write(LuWr,*) 'but only ',nRoots,' exists.'
         Write(LuWr,*) '*******************************************'
         Call Abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (nRoots .gt. 1) Then
         Call Get_dArray('Last energies',Energies_Ref,nRoots)
      Else
         Call Get_dScalar('Last energy',Energies_Ref(iRoot))
      End If
      Call Get_iScalar('Unique atoms',nAtoms)
      Call mma_allocate(Coor,3,nAtoms,Label='Coor')
      Call Get_cArray('Unique Atom Names',AtomLbl,LENIN*nAtoms)
      Call Get_dArray('Unique Coordinates',Coor,3*nAtoms)
      If (iPL_Save.ge.3)
     &     Call RecPrt('Original coordinates',' ',Coor,3,nAtoms)
*                                                                      *
************************************************************************
*                                                                      *
*     If this is a QM/MM calculation, the gradient on the MM atoms
*     is analytical (using the ESPF method for the QM/MM electrostatics)
*     Accordingly, no need to loop over the MM atoms
*
      DoTinker = .False.
      DoDirect = .False.
      ipMltp = ip_Dummy
      Call F_Inquire('ESPF.DATA',Exist)
      If (Exist) Then
         IPotFl = IsFreeUnit(15)
         Call Molcas_Open(IPotFl,'ESPF.DATA')
         Line = ' '
         Do While (Index(Line,'ENDOFESPF ') .eq. 0)
            Line = Get_Ln(IPotFl)
            If (Index(Line,'TINKER ') .ne. 0) Then
              DoTinker = .True.
            Else If (Index(Line,'DIRECT ') .ne. 0) Then
              DoDirect = .True.
            Else If (Index(Line,'MLTORD ') .ne. 0) Then
              Call Get_I1(2,MltOrd)
              ibla = 0
              Do ii = 0, MltOrd
                 ibla = ibla + (ii+2)*(ii+1)/2
              End Do
              MltOrd = ibla
            Else If (Index(Line,'MULTIPOLE ') .ne. 0) Then
              Call Get_I1(2,nMult)
              Call mma_allocate(Mltp,nMult,Label='Mltp')
              ipMltp = ip_of_Work(Mltp)
              Do iMlt = 1, nMult, MltOrd
                Line = Get_Ln(IPotFl)
                Call Get_I1(1,iAt)
                Call Get_F(2,Mltp(iMlt),MltOrd)
              End Do
            End If
         End Do
         Close (IPotFl)
      End If
*
      nAtMM = 0
      If (DoTinker) Then
         Call mma_allocate(IsMM,nAtoms,Label='IsMM')
         Call MMCount(nAtoms,nAtMM,IsMM)
         If (nAtMM .gt. 0) Then
            iQMChg = 0
            StandAlone = .False.
            Call RunTinker(nAtoms,Coor,ipMltp,IsMM,MltOrd,
     &               DynExtPot,iQMchg,iBlabla,StandAlone,DoDirect)
            Call mma_allocate(MMGrd,3,nAtoms,Label='MMGrd')
            MMGrd(:,:)=Zero
            ITkQMMM=IsFreeUnit(15)
            Call Molcas_Open(ITkQMMM,'QMMM')
            Line = ' '
            Do While (Index(Line,'TheEnd') .eq. 0)
               Line = Get_Ln(ITkQMMM)
               If (Index(Line,'MMGradient') .ne. 0) Then
                  Call Get_I1(2,iAtom)
                  Call Get_F(3,FX,3)
                  If (IsMM(iAtom) .eq. 1) MMGrd(:,iAtom)=FX(:)
               End If
            End Do
            Call DScal_(3*nAtoms,Angstrom*ToHartree,MMGrd,1)
            Close(ITkQMMM)
            If (iPL_Save .ge. 3) Call RecPrt('MM Grad:',' ',MMGrd,3,
     &                                       nAtoms)
         End If
         If (Allocated(Mltp)) Call mma_deallocate(Mltp)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Pick up rDelta from the runfile
      Call Get_dScalar('Numerical Gradient rDelta',rDelta)
*                                                                      *
************************************************************************
*                                                                      *
*     Check if there is a coordinate list from Slapaf on the run file.
*     Read the B-matrix and T-matrix. To be used in case of
*     differentiation in internal coordinates according to Slapaf.
*
      If (nAtMM .eq. 0) Call GenCxCTL(iRC,NMCart,rDelta)
*
      Call qpg_dArray('CList',Found,nCList)
      If (Found) Then
         External_Coor_List=.True.
         Call Get_iScalar('No of Internal Coordinates',mInt)
         Call Get_iScalar('nLambda',nLambda)
         nBMtrx=(3*nAtoms)*mInt
         Call mma_allocate(BMtrx,3*nAtoms,mInt,Label='BMtrx')
         Call mma_allocate(TMtrx,mInt,mInt,Label='TMtrx')
         Call Get_dArray('BMtrx',BMtrx,SIZE(BMtrx))
         Call Get_dArray('T-Matrix',TMtrx,mInt**2)
      Else
         NMCart=.FALSE.
         External_Coor_List=.False.
         mInt=3*nAtoms
         nLambda=0
         nBMtrx=(3*nAtoms)**2
         Call mma_allocate(BMtrx,3*nAtoms,3*nAtoms,Label='BMtrx')
         Call mma_allocate(TMtrx,mInt,mInt,Label='TMtrx')
         BMtrx(:,:)=Zero
         Do i = 1, 3*nAtoms
            BMtrx(i,i)=One
         End Do
         TMtrx(:,:)=Zero
         Do i = 1, mInt
            TMtrx(i,i)=One
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (Method(5:7) .eq. 'SCF'    .OR.
     &    Method(1:6) .eq. 'KS-DFT' .OR.
     &    Method(1:6) .eq. 'CASSCF' .OR.
     &    Method(1:6) .eq. 'RASSCF' .OR.
     &    Method(1:6) .eq. 'GASSCF' .OR.
     &    Method(1:6) .eq. 'CASPT2' .OR.
     &    Method(1:5) .eq. 'MBPT2'  .OR.
     &    Method(1:5) .eq. 'CCSDT'  .OR.
     &    Method(1:4) .eq. 'CHCC'   .OR.
     &    Method(1:6) .eq. 'MCPDFT' .OR.
     &    Method(1:4) .eq. 'CHT3'   .OR.
     &    Method(1:8) .eq. 'EXTERNAL') Then
         If (iPL_Save.ge.3) Then
            Write (LuWr,*)
            Write (LuWr,'(A,A,A)')
     &    ' Numerical_Gradient: Original ',Method,' Energies:'
            Write (LuWr,'(G21.14)') Energies_Ref(:)
            Write (LuWr,*)
         End If
      Else
         Write (LuWr,'(A,A,A)') 'Numerical gradient for ',Method,
     &                    ' is not implemented yet.'
         Call Abend()
      End If
*
      nDisp2 = 2*3*nAtoms
      Call mma_Allocate(EnergyArray,nRoots,nDisp2)
      Call FZero(EnergyArray,nRoots*nDisp2)
*                                                                      *
************************************************************************
*                                                                      *
*     Pick up symmetry information
*
      Call Get_iScalar('nSym',nSym)
      nIrrep=nSym
      Call Get_iArray('Symmetry operations',iOper,nSym)
      MaxDCR = nIrrep
*                                                                      *
************************************************************************
*                                                                      *
*     Set up displacement vector
*
      Call Get_nAtoms_All(nAll)
      Call mma_allocate(All,3,nAll,Label='All')
      Call Get_Coord_All(All,nAll)
*                                                                      *
************************************************************************
*                                                                      *
      If (External_Coor_List) Then
*
*        Externally define displacement list
*
         nDisp=mInt
         Call mma_allocate(Disp,mInt,Label='Disp')
         Call Get_dArray('DList',Disp,mInt)
*        Call RecPrt('Dlist',' ',Disp,1,mInt)
*
      Else
*
*        Cartesian displacement list
*
         nDisp=3*nAtoms
         Call mma_allocate(Disp,nDisp,Label='Disp')
         Call FZero(Disp,nDisp)
*
         Do i = 1, nAtoms
*
*           Find the stabilizer of this center
*
            iChxyz=iChAtm(Coor(:,i))
            Call Stblz(iChxyz,nStab,jStab,MaxDCR,iCoSet)
*
            Call IZero(iDispXYZ,3)
            Do j = 0, nStab-1
               If (iAnd(jStab(j),1).ne.0) Then
                  iDispXYZ(1)=iDispXYZ(1)-1
               Else
                  iDispXYZ(1)=iDispXYZ(1)+1
               End If
               If (iAnd(jStab(j),2).ne.0) Then
                  iDispXYZ(2)=iDispXYZ(2)-1
               Else
                  iDispXYZ(2)=iDispXYZ(2)+1
               End If
               If (iAnd(jStab(j),4).ne.0) Then
                  iDispXYZ(3)=iDispXYZ(3)-1
               Else
                  iDispXYZ(3)=iDispXYZ(3)+1
               End If
            End Do
*
*           If this is a MM atom, do not make displacements
*
            If (DoTinker .and. IsMM(i) .eq. 1)
     &        Call IZero(iDispXYZ,3)
            DispX=iDispXYZ(1).ne.0
            DispY=iDispXYZ(2).ne.0
            DispZ=iDispXYZ(3).ne.0
            x0= Coor(1,i)
            y0= Coor(2,i)
            z0= Coor(3,i)
*
*           Find the shortest distance to another atom!
*
            rMax=1.0D19
            Do j = 1, nAll
               x = All(1,j)
               y = All(2,j)
               z = All(3,j)
               rTest=Sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
               If (rTest.eq.Zero) rTest=1.0D19
               rMax=Min(rMax,rTest)
            End Do
            If (DispX) Disp((i-1)*3+1)=rDelta*rMax
            If (DispY) Disp((i-1)*3+2)=rDelta*rMax
            If (DispZ) Disp((i-1)*3+3)=rDelta*rMax
         End Do
*        Call RecPrt('Disp',' ',Disp,1,nDisp)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(Deg,3,nAtoms,Label='Deg')
      Deg(:,:)=Zero
      Do i = 1, nAtoms
         rDeg=DBLE(iDeg(Coor(:,i)))
         Deg(1,i)=rDeg
         Deg(2,i)=rDeg
         Deg(3,i)=rDeg
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(XYZ,3*nAtoms,2*nDisp,Label='XYZ')
      XYZ(:,:)=Zero
      If (External_Coor_List) Then
         mDisp=nDisp*2
         Call Get_dArray('CList',XYZ,3*nAtoms*mDisp)
      Else
         mDisp=0
         Do iDisp = 1, nDisp2
            iCoor=(iDisp+1)/2
            If (Disp(icoor).eq.Zero) Go To 101
            mDisp=mDisp+1
*
*--------   Modify the geometry
*
            jpXYZ = ipXYZ + (mDisp-1)*3*nAtoms
            call dcopy_(3*nAtoms,Coor,1,XYZ(:,mDisp),1)
            Sign=One
            If (Mod(iDisp,2).eq.0) Sign=-One
            XYZ(iCoor,mDisp) = XYZ(iCoor,mDisp)
     &                          + Sign*Disp(icoor)
 101        Continue
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Save the "new geometry" field from the RunFile, if any
*
      Call qpg_dArray('GeoNew',Found,nGNew)
      If (.not.Found) nGNew=0
      If (nGNew.gt.0) Then
        Call mma_allocate(GNew,nGNew)
        Call Get_dArray('GeoNew',GNew,nGNew)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Save global print level
*
      iPL_Save=iPrintLevel(-1)
*     iPL_Base=0
*     If (iPL_Save.ge.3) iPl_Base=iPL_Save
*
#ifdef _DEBUGPRINT_
      Call RecPrt('BMtrx',' ',BMtrx,3*nAtoms,mInt)
      Call RecPrt('TMtrx',' ',TMtrx,mInt,mInt)
      Call RecPrt('Degeneracy vector',' ',Deg,3,nAtoms)
      Call RecPrt('Coordinate List',' ',XYZ,3*nAtoms,mDisp)
#endif

*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Loop over displacements
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      If (iPL.ge.2) Then
         Write (LuWr,*) 'Root to use: ',iRoot
         Write (LuWr,*) 'Number of internal degrees            ', nDisp
         Write (LuWr,*) 'Number of constraints                 ',
     &                  nLambda
         Write (LuWr,*) 'Number of displacements               ',
     &                   nDisp*2
         If (nAtMM.ne.0 .and. DoTinker .and. .not.DoDirect) Then
            Write(LuWr,*) 'Number of MM degrees (analytical grad)',
     &                    nAtMM*3
         End If
         Write (LuWr,*) 'Effective number of displacements     ',
     &                   2*(nDisp-nLambda-3*nAtMM)
         Write (LuWr,*) 'Relative displacements                ',
     &                   rDelta
         Write (LuWr,*)
      End If
*
*     This printout is disabled since it means different things for
*     externally defined coordinates or not.
*
      If (iPL.ge.3) Then
#ifdef _HIDE_
         Write (LuWr,'(1x,A)')
     &                 '---------------------------------------------'
         Write (LuWr,'(1x,A)')
     &                 '               X           Y           Z     '
         Write (LuWr,'(1x,A)')
     &                 '---------------------------------------------'
         Do iAtom = 1, nAtoms
            TempX = Disp(3*(iAtom-1)+1)
            TempY = Disp(3*(iAtom-1)+2)
            TempZ = Disp(3*(iAtom-1)+3)
            Namei = AtomLbl(iAtom)
            Write (LuWr,'(2X,A,3X,3F12.6)') Namei, TempX, TempY, TempZ
         End Do
#endif
         Write (LuWr,'(1x,A)')
     &              '---------------------------------------------'
         Write (LuWr,*)
         Write (LuWr,*) 'Here we go ...'
         Write (LuWr,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Change output unit
*
      LuWr_save=LuWr
      If (MyRank.ne.0) Then
         LuWr=55
         LuWr=isFreeUnit(LuWr)
         call molcas_open(luwr,'Temp_OutPut')
      End If
*
* FM 16/4/2013
* ESPF charges are set to zero so that no microiterations will be
* performed during numerical gradient
* IFG: swap files, as now NG uses a subdirectory
*
      iSave=15
      If (Do_ESPF) Then
         iSave = IsFreeUnit(iSave)
         Call Molcas_Open(iSave,'ESPF.SAV')
         iData = IsFreeUnit(iSave)
         Call Molcas_Open(iData,'ESPF.DATA')
96       Line = Get_Ln(iData)
         ESPFKey = Line(1:10)
         If (ESPFKey.eq.'ENDOFESPF ') Then
            Write(iSave,'(A132)') Line
            Goto 98
         End If
         If (ESPFKey.ne.'MULTIPOLE ') Then
            Write(iSave,'(A132)') Line
         Else
            Call Get_I1(2,nMult)
            Do iMlt = 1, nMult
               Line = Get_Ln(iData)
            End Do
         End If
         Goto 96
98       Close (iSave)
         Close (iData)
      EndIf
* FM End
*                                                                      *
************************************************************************
*                                                                      *
* Parallel loop over the nDisp displacements.
* Reserve task on global task list and get task range in return.
* Function will be false if no more tasks to execute.
#if !defined(_GA_) && defined(_MOLCAS_MPP_)
      Call getenvf('MOLCAS_SSTMNGR',SSTMNGR)
      If (SSTMNGR(1:1).eq.'Y') Then
          SSTMODE=1
      Else
          SSTMODE=0
      End If
      If (SSTMODE.eq.1) Then
         Call Init_Tsk(id_Tsk,mDisp-2*nLambda)
      Else
         Call Init_Tsk_Even(id_Tsk,mDisp-2*nLambda)
      End If
#else
      Call Init_Tsk(id_Tsk,mDisp-2*nLambda)
#endif
      Call mma_allocate(C,3,nAtoms,Label='C')
   10 Continue
#if defined (_MOLCAS_MPP_) && !defined(_GA_)
         If (SSTMODE.eq.1) Then
           If (.Not. Rsv_Tsk(id_Tsk,iDisp)) Goto 11
         Else
           If (.Not. Rsv_Tsk_Even(id_Tsk,iDisp)) Goto 11
         End If
#else
         If (.Not. Rsv_Tsk(id_Tsk,iDisp)) Goto 11
#endif
*
*        Offset for the constraints
*
         iDisp = iDisp + 2*nLambda
*
*------- Get the displaced geometry
*
         call dcopy_(3*nAtoms,xyz(:,iDisp),1,C,1)
*        Call RecPrt('C',' ',C,3*nAtoms,1)
         Call Put_Coord_New(C,nAtoms)
*
*        Compute integrals
*
*        jPL     =iPrintLevel(Max(iPL_Base,0)) ! Silent
         Call Set_Do_Parallel(.False.)
*
*        Switch to a new directory
*        WARNING WARNING WARNING WARNING WARNING
*          ugly hack, do not try this at home
         Call SubWorkDir()
*        WARNING WARNING WARNING WARNING WARNING
*
         Call StartLight('seward')
         Call init_run_use()
         Call init_ppu(.True.)
         Call Disable_Spool()
         Call Seward(ireturn)
         Call ReClose()
         If (iReturn .ne. 0) Then
            Write(LuWr,*) 'Numerical_Gradient failed ...'
            Write(LuWr,*)'Seward returned with return code, rc = ',
     &                    iReturn
            Write(LuWr,*) 'for the perturbation iDisp = ',iDisp
            Call Abend()
         End If
*
*        Compute the ESPF stuff
*
         If (Do_ESPF) Then
            Call StartLight('espf')
            Call init_run_use()
            Call Disable_Spool()
            StandAlone=.True.
            Call ESPF(ireturn,StandAlone)
            Call ReClose()
            If (iReturn .ne. 0) Then
               Write(LuWr,*) 'Numerical_Gradient failed ...'
               Write(LuWr,*)'ESPF returned with return code, rc = ',
     &                       iReturn
               Write(LuWr,*) 'for the perturbation iDisp = ',iDisp
               Call Abend()
            End If
         End If
*
*        Compute the wave function
*
         If (Method(5:7) .eq. 'SCF' .OR.
     &       Method(1:6) .eq. 'KS-DFT' .OR.
     &       Method(1:5) .eq. 'MBPT2' .OR.
     &       Method(1:4) .eq. 'CHCC' .OR.
     &       Method(1:4) .eq. 'CHT3') Then
            Call StartLight('scf')
            Call init_run_use()
            Call Disable_Spool()
            Call xml_open('module',' ',' ',0,'scf')
            Call SCF(iReturn)
            Call xml_close('module')
            Call ReClose()
            If (iReturn .ne. 0) Then
               Write(LuWr,*) 'Numerical_Gradient failed ...'
               Write(LuWr,*) 'SCF returned with return code, rc = ',
     &                     iReturn
               Write(LuWr,*) 'for the perturbation iDisp = ',iDisp
               Call Abend()
            End If
         Else If (Method(1:6) .eq. 'RASSCF' .OR.
     &            Method(1:6) .eq. 'GASSCF' .OR.
     &            Method(1:6) .eq. 'CASSCF' .OR.
     &            Method(1:6) .eq. 'MCPDFT' .OR.
     &            Method(1:6) .eq. 'CASPT2' .OR.
     &            Method(1:5) .eq. 'CCSDT') Then
            Call StartLight('rasscf')
            Call init_run_use()
            Call Disable_Spool()
            Call RASSCF(ireturn)
            Call ReClose()
            If (iReturn .ne. 0) Then
               Write(LuWr,*) 'Numerical_Gradient failed ...'
               Write(LuWr,*) 'RASSCF returned with return code, rc = ',
     &                     iReturn
               Write(LuWr,*) 'for the perturbation iDisp = ',iDisp
               Call Abend()
            End If
         Else If (Method(1:8) .eq. 'EXTERNAL') Then
            Call StartLight('false')
            Call init_run_use()
            Call Disable_Spool()
            Call False_program(ireturn)
            Call ReClose()
            If (iReturn .ne. 0) Then
               Write(LuWr,*) 'Numerical_Gradient failed ...'
               Write(LuWr,*) 'FALSE returned with return code, rc = ',
     &                     iReturn
               Write(LuWr,*) 'for the perturbation iDisp = ',iDisp
               Call Abend()
            End If
         End If
*
         If (Method(1:5) .eq. 'MBPT2') Then
            Call StartLight('mbpt2')
            Call init_run_use()
            Call Disable_Spool()
            Call MP2_Driver(ireturn)
            Call ReClose()
            If (iReturn .ne. 0) Then
               Write(LuWr,*) 'Numerical_Gradient failed ...'
               Write(LuWr,*) 'MBPT2 returned with return code, rc = ',
     &                     iReturn
               Write(LuWr,*) 'for the perturbation iDisp = ',iDisp
               Call Abend()
            End If
         End If
*
         If (Method(1:5) .eq. 'CCSDT') Then
            Call StartLight('motra')
            Call init_run_use()
            Call Disable_Spool()
            Call Motra(ireturn)
            Call ReClose()
            If (iReturn .ne. 0) Then
               Write(LuWr,*) 'Numerical_Gradient failed ...'
               Write(LuWr,*) 'Motra returned with return code, rc = ',
     &                     iReturn
               Write(LuWr,*) 'for the perturbation iDisp = ',iDisp
               Call Abend()
            End If
*
            Call StartLight('ccsdt')
            Call init_run_use()
            Call Disable_Spool()
            Call CCSDT(ireturn)
            Call ReClose()
            If (iReturn .ne. 0) Then
               Write(LuWr,*) 'Numerical_Gradient failed ...'
               Write(LuWr,*) 'CCSDT returned with return code, rc = ',
     &                     iReturn
               Write(LuWr,*) 'for the perturbation iDisp = ',iDisp
               Call Abend()
            End If
         End If
*
         If (Method(1:4) .eq. 'CHCC' .OR.
     &       Method(1:4) .eq. 'CHT3') Then
            Call StartLight('chcc')
            Call init_run_use()
            Call Disable_Spool()
            Call CHCC(ireturn)
            Call ReClose()
            If (iReturn .ne. 0) Then
               Write(LuWr,*) 'Numerical_Gradient failed ...'
               Write(LuWr,*) 'CHCC returned with return code, rc = ',
     &                     iReturn
               Write(LuWr,*) 'for the perturbation iDisp = ',iDisp
               Call Abend()
            End If
         End If
*
         If (Method(1:4) .eq. 'CHT3') Then
            Call StartLight('cht3')
            Call init_run_use()
            Call Disable_Spool()
            Call CHT3(ireturn)
            Call ReClose()
            If (iReturn .ne. 0) Then
               Write(LuWr,*) 'Numerical_Gradient failed ...'
               Write(LuWr,*) 'CHT3 returned with return code, rc = ',
     &                     iReturn
               Write(LuWr,*) 'for the perturbation iDisp = ',iDisp
               Call Abend()
            End If
         End If
*
         If (Method(1:6) .eq. 'CASPT2') Then
            Call StartLight('caspt2')
            Call init_run_use()
            Call Disable_Spool()
            Call CASPT2(ireturn)
            Call ReClose()
            If (iReturn .ne. 0) Then
               Write(LuWr,*) 'Numerical_Gradient failed ...'
               Write(LuWr,*) 'CASPT2 returned with return code, rc = ',
     &                     iReturn
               Write(LuWr,*) 'for the perturbation iDisp = ',iDisp
               Call Abend()
            End If
         End If

         If (Method(1:6) .eq. 'MCPDFT') Then
            Call StartLight('mcpdft')
            Call init_run_use()
            Call Disable_Spool()
            Call MCPDFT(ireturn)
            Call ReClose()
            If (iReturn .ne. 0) Then
               Write(LuWr,*) 'Numerical_Gradient failed ...'
               Write(LuWr,*) 'MCPDFT returned with return code, rc = ',
     &                     iReturn
               Write(LuWr,*) 'for the perturbation iDisp = ',iDisp
               Call Abend()
            End If
         End If
*
*        Call Get_Energy(EnergyArray(iRoot,iDisp))
         If (nRoots .gt. 1) Then
            Call Get_dArray('Last energies',EnergyArray(1,iDisp),nRoots)
         Else
            Call Get_dScalar('Last energy',EnergyArray(iRoot,iDisp))
         End If
*
*        Restore directory and prgm database
*
         Call ParentWorkDir()
         Call prgmfree()
         Call prgminit('numerical_gradient')
*
         Call Set_Do_Parallel(.True.)

         If (iPL.ge.2) Then
            Write (LuWr,200)'   * Point #',iDisp-2*nLambda,' of ',
     &                                2*(nDisp-nLambda-3*nAtMM),' done.'
            Write (LuWr,200)'    (Perturbation ',iDisp,')'
         End If
 200     Format(A,I4,A,I4,A)
*                                                                      *
************************************************************************
*                                                                      *
#if defined (_MOLCAS_MPP_) && !defined(_GA_)
         If(King().and.nProcs.gt.1.and.SSTMODE.eq.1) Go To 11
#endif
         Go To 10
  11  Continue
C_MPP End Do
#if !defined(_GA_) && defined(_MOLCAS_MPP_)
      If (SSTMODE.eq.1) Then
         Call Free_Tsk(id_Tsk)
      Else
         Call Free_Tsk_Even(id_Tsk)
      End If
#else
      Call Free_Tsk(id_Tsk)
#endif
      Call GADSum(EnergyArray,nRoots*mDisp)
      Call mma_deallocate(C)
      Call mma_deallocate(XYZ)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     End of Loop                                                      *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Flush the output from the other nodes.
*
      Do iRank = 1, nProcs-1
         Call GASync()
         If (iRank.eq.MyRank) Then
            ReWind(LuWr)
 777        Read(LuWr,'(A)',END=778) Line
            Write (LuWr_Save,*) Line
            Go To 777
 778        Continue
         End If
         Call GASync()
      End Do
      If (MyRank.ne.0) Then
         Close(LuWr)
         LuWr=LuWr_save
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Restore the "new geometry" field to the RunFile, if any
*
      If (nGNew.eq.0) Then
        Call Put_Coord_New([Zero],0)
      Else
        Call Put_Coord_New(GNew,nGNew/3)
        Call mma_deallocate(GNew)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If(nProcs.ge.2. and. iPL.ge.2) Then
         Write (LuWr,*)
         Write (LuWr,*) ' Points were printed only by master node'
         Write (LuWr,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Read old gradient(s) and convert to internal coordinates
*
      Call mma_Allocate(OldGrads,3*nAtoms,nRoots)
      Call FZero(OldGrads,3*nAtoms*nRoots)
      Call Get_lScalar('Keep old gradient',KeepOld)
      If (KeepOld) Then
        Call Query_Grads(Found,i,j)
        If (Found.and.(i.ge.nRoots).and.(j.eq.3*nAtoms)) Then
*         If there is a valid GRADS file, read each gradient
          Do iR=1,nRoots
            rc=Read_Grad(OldGrads(1,iR),3*nAtoms,iR,0,0)
            If (rc.le.0) Then
              Write(LuWr,*)
              Write(LuWr,'(2X,A,I4,A)') 'No gradient found for root ',
     &                                  iR,', using 0'
              Write(LuWr,*)
            End If
          End Do
        Else
*         If the GRADS file does not exist or has the wrong sizes,
*         read a single gradient from the runfile
          Call qpg_dArray('GRAD',Found,nGrad)
          If (Found) Then
            Call Get_dArray('GRAD',OldGrads(1,1),nGrad)
            Do iR=2,nRoots
              Call dCopy_(nGrad,OldGrads(1,1),1,OldGrads(1,iR),1)
            End Do
            If (nRoots.gt.1) Then
              Write(LuWr,*)
              Write(LuWr,'(2X,A)') 'Using stored gradient for all roots'
              Write(LuWr,*)
            End If
          Else
            Write(LuWr,*)
            Write(LuWr,'(2X,A)') 'No gradient found, using 0'
            Write(LuWr,*)
          End If
        End If
      End If
      Call mma_allocate(Tmp2,nDisp,Label='Tmp2')
      Do iR=1,nRoots
        Call Eq_Solver('N',3*nAtoms,nDisp,1,BMtrx,.True.,
     &                 rDum(1),OldGrads(1,iR),Tmp2)
        Call FZero(OldGrads(1,iR),3*nAtoms)
        Call Eq_Solver('N',nDisp,nDisp,1,TMtrx,.True.,
     &                 rDum(1),Tmp2,OldGrads(1,iR))
      End Do
      Call mma_deallocate(Tmp2)
*                                                                      *
************************************************************************
*                                                                      *
*     jPL=iPrintLevel(iPL_Save)
      If (iPL_Save.ge.3)
     &    Call RecPrt('Energies','(8G16.10)',EnergyArray,nRoots,mDisp)
      Call mma_Allocate(GradArray,nDisp,nRoots)
      Call FZero(GradArray,nDisp*nRoots)
      Call mma_Allocate(Grad,nRoots)
*
      iDisp = nLambda
      Do i = 1, nDisp
*
         If (Disp(i).eq.Zero) then
            If (iPL_Save .ge. 3) Then
               If (KeepOld) Then
                  Write(6,*) 'gradient set to old value'
               Else
                  Write(6,*) 'gradient set to zero'
               End If
            End If
            Do iR = 1, nRoots
               Grad(iR) = OldGrads(i,iR)
            End Do
         Else
            iDisp = iDisp + 1
            Do iR = 1, nRoots
               iEp = iDisp*2-1
               Eplus=EnergyArray(iR,iEp)
               iEm = iDisp*2
               EMinus=EnergyArray(iR,iEm)
               Dsp=Disp(i)
               Grad(iR) = (Eplus-EMinus)/(Two*Dsp)
*
*              If the gradient is not close to zero check that it is
*              consistent. For CASPT2/CASSCF sometimes the active space
*              breaks down and the computed gradient is rubbish. If we
*              are lucky this happens only in one of the displacement
*              directions. In that case compute the gradient with the
*              one-point equation. The one with the lowest gradient is
*              the one which is most likely to be correct.
*
               If (Abs(Grad(iR)).gt.1.0D-1) Then
                  Energy_Ref=Energies_Ref(iR)
                  Grada= (Eplus-Energy_Ref)/Dsp
                  Gradb= (Energy_Ref-EMinus)/Dsp
                  If (Abs(Grad(iR)/Grada).gt.1.5D0 .or.
     &                Abs(Grad(iR)/Gradb).gt.1.5D0 ) Then
                     If (Abs(Grada).le.Abs(Gradb)) Then
                        Grad(iR)=Grada
                     Else
                        Grad(iR)=Gradb
                     End If
                  End If
               End If
            End Do
         End If
*
         call dCopy_(nRoots,Grad,1,GradArray(i,1),nDisp)
*
      End Do
*     Call RecPrt('Grads (old)',' ',OldGrads,3*nAtoms,nRoots)
*     Call RecPrt('BMtrx',' ',BMtrx,3*nAtoms,nDisp)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the gradient in Cartesian coordinates
*
      Call mma_allocate(Tmp2,nDisp,Label='Tmp2')
      Call mma_allocate(Tmp,3,nAtoms,Label='Tmp')
      Do iR = 1, nRoots
         Tmp2(:)=Zero
         Tmp(:,:)=Zero
*
*        Transform the gradient in PCO basis to the internal coordinate
*        format.
*
         Call DGEMM_('N','N',nDisp,1,nDisp,
     &               1.0D0,TMtrx,nDisp,
     &                     GradArray(1,iR),nDisp,
     &               0.0D0,Tmp2,nDisp)
*
*        Transform internal coordinates to Cartesian.
*
         Call DGEMM_('N','N',3*nAtoms,1,nDisp,
     &               1.0d0,BMtrx,3*nAtoms,
     &                     Tmp2,nDisp,
     &               0.0d0,Tmp,3*nAtoms)
*        Call RecPrt('Tmp',' ',Tmp,3,nAtoms)
*
*        Modify with degeneracy factors.
*
         If (.NOT.NMCart) Then
            Do iAtom = 1, nAtom
               Do ixyz = 1, 3
                  Tmp(ixyz,iAtom) = Tmp(ixyz,iAtom)/Deg(ixyz,iAtom)
               End Do
            End Do
         End If
*
*        Add the MM contribution for MM atoms
*
         If (nAtMM .ne. 0) Then
            call daxpy_(3*nAtoms,One,MMGrd,1,Tmp,1)
         End If
*
*        Apply Morokuma's scheme if needed
*
         Call F_Inquire('QMMM',Exist)
         If (Exist .and. DoTinker) Then
            Call LA_Morok(nAtoms,Tmp,1)
         End If
*
         If (iR.eq.iRoot) Call Put_Grad(Tmp,3*nAtoms)
         Call Add_Info('Grad',Tmp,3*nAtoms,6)
         Call Store_grad(Tmp,3*nAtoms,iR,0,0)
*
         If (iPL.ge.2) Then
            Write (LuWr,*)
            Write (LuWr,'(1x,A,I5)') 'Numerical gradient, root ',iR
            Write (LuWr,'(2x,A)')
     &                  '---------------------------------------------'
            Write (LuWr,'(2x,A)')
     &                  '               X           Y           Z'
            Write (LuWr,'(2x,A)')
     &                  '---------------------------------------------'
            Do iAtom = 1, nAtoms
               TempX = Tmp(1,iAtom)
               TempY = Tmp(2,iAtom)
               TempZ = Tmp(3,iAtom)
               Namei = AtomLbl(iAtom)
               Write (LuWr,'(2X,A,3X,3F12.6)')
     &               Namei, TempX, TempY, TempZ
            End Do
            Write (LuWr,'(2x,A)')
     &                  '---------------------------------------------'
         End If
      End Do
      Call mma_deallocate(Tmp2)
      Call mma_deallocate(Tmp)
*                                                                      *
************************************************************************
*                                                                      *
*     Deallocate
*
      Call mma_deallocate(Deg)
      Call mma_deallocate(Disp)
      Call mma_deallocate(All)
      Call mma_Deallocate(TMtrx)
      Call mma_Deallocate(BMtrx)
      Call mma_Deallocate(EnergyArray)
      Call mma_Deallocate(GradArray)
      Call mma_Deallocate(OldGrads)
      Call mma_Deallocate(Grad)
      Call mma_deallocate(Coor)
      Call mma_deallocate(Energies_Ref)
      If (DoTinker) Call mma_deallocate(IsMM)
      If (nAtMM.gt.0) Call mma_deallocate(MMGrd)
*
*     Since Numerical_Gradient itself cleans up after the modules
*     then it's necessary to force nfld_tim and nfld_stat to zero,
*     or finish will scream.
*
      nfld_tim  = 0
      nfld_stat = 0
*
*     Restore iRlxRoot if changed as set by the RASSCF module.
*
      If (Method(1:6) .eq. 'CASSCF' .OR.
     &    Method(1:6) .eq. 'RASSCF' ) Then
         Call Get_iScalar('Relax CASSCF root',irlxroot1)
         Call Get_iScalar('Relax Original root',irlxroot2)
         If (iRlxRoot1.ne.iRlxRoot2) Then
            Call Put_iScalar('Relax CASSCF root',irlxroot2)
            Call Put_iScalar('NumGradRoot',irlxroot2)
         End If
      End If

*
      Return
      End
