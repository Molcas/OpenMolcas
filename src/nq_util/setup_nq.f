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
* Copyright (C) 1999, Roland Lindh                                     *
************************************************************************
      Subroutine Setup_NQ(Maps2p,nShell,nSym,nNQ,Do_Grad,On_Top,nD,
     &                    Pck_Old,PMode_old,R_Min,nR_Min)
************************************************************************
*                                                                      *
* Object: to set up information for calculation of integrals via a     *
*         numerical quadrature.                                        *
* Warning: The exponents of each shell are reordered diffuse to compact*
*                                                                      *
* Called from: Drvnq                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              Quit                                                    *
*              Nr_Shells                                               *
*              GetMem                                                  *
*              DSwap                                                   *
*              GauLeg                                                  *
*              RecPrt                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh,                                            *
*             Dept of Chemical Physics,                                *
*             University of Lund, Sweden                               *
*             August 1999                                              *
************************************************************************
      use Real_Spherical
      use iSD_data
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "status.fh"
#include "nq_info.fh"
#include "nsd.fh"
#include "setup.fh"
#include "grid_on_disk.fh"
#include "print.fh"
      Real*8 Coor(3), C(3)
      Logical EQ, Do_Grad, On_Top, PMode_Old
      Real*8 Alpha(2),rm(2), R_Min(0:nR_Min)
      Integer Maps2p(nShell,0:nSym-1)
      Integer iDCRR(0:7)
      Dimension Dummy(1)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement Functions
*
#include "nq_structure.fh"
      declare_ip_coor
      declare_ip_a_high
      declare_ip_a_low
      declare_ip_r_rs
      declare_ip_l_max
      declare_ip_r_quad
      declare_ip_angular
      declare_ip_atom_nr
      declare_ip_dodx
      nElem(i)=(i+1)*(i+2)/2
*                                                                      *
************************************************************************
*                                                                      *
      Call ICopy(nShell*nSym,[-99999999],0,Maps2p,1)
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
c     Write(6,*) '********** Setup_NQ ***********'
C     Call QEnter('Setup_NQ')
      ntotgp=0
*                                                                      *
************************************************************************
*                                                                      *
*-----Check if NQ environment has been activated
*
      If (NQ_Status.ne.Active.and.NQ_Status.ne.Inactive) Then
         Call WarningMessage(2,'Setup_NQ: NQ_Status not initialized')
         Call Quit_OnUserError()
      End If
      If (NQ_Status.eq.Active) Return
      NQ_Status=Active
*                                                                      *
************************************************************************
*                                                                      *
*---- Get the coordinates to the centers of all Voronoi polyhedra
*
*     Note that this will be all centers with valence basis sets on
*     them. Hence this will also include any pseudo centers!
*
      Call Allocate_Work(ipTempC,3*nShell*nSym)
      nAtoms = 0
      If (nShell.gt.nskal_iSD) Then
         Write (6,*) 'nShell.gt.nSkal_iSD'
         Write (6,*) 'nShell=',nShell
         Write (6,*) 'nSkal_iSD=',nSkal_iSD
         Call AbEnd()
      End If
      Do iShell = 1, nShell
         iCnttp=iSD(13,iShell)
         iCnt  =iSD(14,iShell)
         Coor(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
         Call Process_Coor(Coor,Work(ipTempC),nAtoms,nSym,iOper)
      End Do
      Call Allocate_Work(ipCoor,3*nAtoms)
      call dcopy_(3*nAtoms,Work(ipTempC),1,Work(ipCoor),1)
C     Call RecPrt('Coor',' ',Work(ipCoor),3,nAtoms)
      Call Free_Work(ipTempC)
*                                                                      *
************************************************************************
*                                                                      *
*---- Get the symmetry unique coordinates
*
      nNQ=nAtoms
      Call GetMem('nq_centers','Allo','Real',ipNQ,nNQ*l_NQ)
      ipTmp1=ipCoor
      Do iNQ = 1, nNQ
         call dcopy_(3,Work(ipTmp1),1,Work(ip_Coor(iNQ)),1)
         ipTmp1=ipTmp1+3
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*-----Pick up the requested accuracy.
*
*     Dr=-Log10(Thr)
*                                                                      *
************************************************************************
*                                                                      *
*-----Loop over each unique center, and find the highest and lowest    *
*     Gaussians exponents associated with this center. The later       *
*     information will be used to design the radial grid associated    *
*     with this center.                                                *
*                                                                      *
*-----Initialize the exponents to extreme values.
      Do iNQ=1,nNQ
         Work(ip_lMax(iNQ))=DBLE(-1)
         Work(ip_A_high(iNQ))=-1.0D99
         Work(ip_A_low(iNQ) )= 1.0D99
      End Do
*
      iAngMax=0
      NbrMxBas=0
      Do iShell=1,nShell
         iAng  =iSD(1,iShell)
         nCntrc=iSD(3,iShell) !Get the # of contracted functions
                              !for iShell
         mExp  =iSD(5,iShell) ! Get the number of exponents of ishell
         NbrMxBas=Max(NbrMxbas,nCntrc)
         iAngMax=Max(iAngMax,iAng)
      End Do
*
*-----Loop over the shells
*
      nMaxExp=0
      nAOMax=0
      Do iShell = 1, nShell
*
*------- Get the Atom number
         iANr=dbsc(iSD(13,iShell))%AtmNr
*
         iShll=iSD(0,iShell)   ! Get the angular momentum of ishell
         iAng=iSD(1,iShell)   ! Get the angular momentum of ishell
         iCmp=iSD(2,iShell)   ! Get the # of angular components
         nCntrc=iSD(3,iShell) !Get the # of contracted functions
                              !for iShell
         mExp=iSD(5,iShell)   ! Get the number of exponents of ishell
         nMaxExp=Max(nMaxExp,mExp)
         nAOMax=Max(nAOMax,iCmp*nCntrc)
*
************************************************************************
*                                                                      *
*     Order the exponents diffuse to compact for the active shell      *
*                                                                      *
************************************************************************
         Call OrdExpD2C(mExp,Shells(iShll)%Exp,nCntrc,
     &                       Shells(iShll)%pCff)
*
*-----Get the extreme exponents for the active shell.
         A_low =Shells(iShll)%Exp(1)
         A_high=Shells(iShll)%Exp(mExp)
*
         iCnttp=iSD(13,iShell)
         iCnt  =iSD(14,iShell)
         C(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
         Do iIrrep = 0, nIrrep-1
            Coor(1) = C(1)*DBLE(iPhase(1,iOper(iIrrep)))
            Coor(2) = C(2)*DBLE(iPhase(2,iOper(iIrrep)))
            Coor(3) = C(3)*DBLE(iPhase(3,iOper(iIrrep)))
         Do iNQ=1,nNQ
            jNQ=ip_Coor(iNQ)
*
            If (EQ(Work(jNQ  ),Coor)) Then
*
               Work(ip_Atom_Nr(iNQ))=DBLE(iANR)
*
*------------- Assign the BS radius to the center
               Work(ip_R_RS(iNQ))=Bragg_Slater(iANr)
*
*------------- What is the maximum angular momentum for the active center ?
               Work(ip_lMax(iNQ))=Max(Work(ip_lMax(iNQ)),DBLE(iAng))
*
*------------- Get the extreme exponents for the atom
               Work(ip_A_high(iNQ))=Max(Work(ip_A_high(iNQ)),A_High)
               Work(ip_A_low(iNQ) )=Min(Work(ip_A_low(iNQ) ),A_low)
*
               Maps2p(iShell,iIrrep)=iNQ
               Go To 100
            End If
         End Do                 ! iNQ
         Call WarningMessage(2,
     &              'Didn''t find a center associated with the shell!')
         Call Abend()
 100     Continue
         End Do                 ! iIrrep
*
      End Do                    !iShell
*
************************************************************************
*                                                                      *
*               END OF THE LOOP OVER THE SHELLS.                       *
*     Now we have the number of unique center and their associated     *
*     exponents and maximum angular momentum. And for each shell the   *
*     exponents are ordered diffuse to compact.                        *
*                                                                      *
************************************************************************
*                                                                      *
*     Loop over all the atoms to create the radial quadrature          *
*                                                                      *
************************************************************************
*
*-----Allocate memory to store the number of effective radial points for
*     each center and the radius of this center.
      Call GetMem('NumRadEff','Allo','Inte',ip_nR_eff,nNQ)
*
      iNQ_MBC=0
      iReset=0
      Threshold_tmp=Zero
      nR_tmp=0
      Do iNQ=1,nNQ
*--------Get the extreme exponents for the atom
         Alpha(1)=Work(ip_A_low(iNQ) )
         Alpha(2)=Work(ip_A_high(iNQ))
*
*--------Get the coordinates of the atom
         call dcopy_(3,Work(ip_Coor(iNQ)),1,Coor,1)
*
*        For a special center we can increase the accuracy.
*
         jNQ=ip_Coor(iNQ)
         If (MBC.ne.' ') Then
            Do iS = 1, nShell
               iCnttp=iSD(13,iS)
               iCnt  =iSD(14,iS)
               C(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
               If ( EQ(Work(jNQ),C) ) Then
                  mdci=iSD(10,iS)
                  Write (6,*) 'LblCnt(mdci)=',LblCnt(mdci)
                  If (LblCnt(mdci).eq.MBC) Then
                     nR_tmp=nR
                     nR=INT(DBLE(nR)*2.0D0)
                     Threshold_tmp=Threshold
                     Threshold=Threshold*1.0D-6
*
                     iReset=1
                     iNQ_MBC=iNQ
                     Go To 1771
                  End If
               End If
            End Do
         End If
 1771    Continue
*
*        Max angular momentum for the atom -> rm(1)
*        Max Relative Error -> rm(2)
         rm(1)=Work(ip_lMax(iNQ))
         rm(2)=Threshold
*
         Call GenVoronoi(Coor,iWork(ip_nR_eff),nNQ,Alpha,rm,iNQ)
*
         If (iReset.eq.1) Then
            nR=nR_tmp
            Threshold=Threshold_tmp
            iReset=0
         End If
*
      End Do  ! iNQ
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the principal axis system and optionally to
*     compute derivatives of the principal axis. Needed in order to
*     compute the gradient of the rotationally invariant DFT energy.
*
      Call GetMem('O','Allo','Real',ip_O,3*3)
      Call GetMem('dOdx','Allo','Real',ipdOdx,3*3*nNQ*3)
      Call FZero(Work(ipdOdx),3*3*nNQ*3)
      Call GetMem('ZA','Allo','Real',ip_ZA,nNQ)
      Call Allocate_Work(ip_Crd,3*nNQ)
*
*     Collect coordinates and charges of the nuclei
*
      iOff = ip_Crd
      Do iNQ = 1, nNQ
         Work(ip_ZA+(iNQ-1))=Work(ip_Atom_Nr(iNQ))
         call dcopy_(3,Work(ip_Coor(iNQ)),1,Work(iOff),1)
         iOff = iOff + 3
      End Do
*
#ifdef _DEBUG_
      If (Do_Grad) Then
      delta=1.0D-8
      Call Allocate_Work(ip_Dbg,9)
      Do iNQ = 1, nNQ
         Do iCar = 1, 3
*
            Write (6,*) 'iCar,iNQ=',iCar,iNQ
*
            Call RotGrd(Work(ip_Crd),Work(ip_ZA),Work(ip_O),
     &                  Dummy,Dummy,nNQ,.False.,.False.)
            Call RecPrt('O(original)',' ',Work(ip_O),3,3)
*
            i3 = (iNQ-1)*3+iCar-1+ip_Crd
            temp=Work(i3)
*
            Work(i3)=temp+delta
            Call RotGrd(Work(ip_Crd),Work(ip_ZA),Work(ip_O),
     &                  Dummy,Dummy,nNQ,.False.,.False.)
            call dcopy_(9,Work(ip_O),1,Work(ip_Dbg),1)
            Call RecPrt('O',' ',Work(ip_O),3,3)
*
            Work(i3)=temp-delta
            Call RotGrd(Work(ip_Crd),Work(ip_ZA),Work(ip_O),
     &                  Dummy,Dummy,nNQ,.False.,.False.)
            Call RecPrt('O',' ',Work(ip_O),3,3)
*
            Work(i3)=temp
*
            Do i = 0, 8
               temp = (Work(ip_Dbg+i) - Work(ip_O+i))/(2.0D0*delta)
               Work(ip_O+i) = temp
            End Do
            Call RecPrt('dOdx(numerical)',' ',Work(ip_O),3,3)
            call dcopy_(9,Work(ip_O),1,Work(ip_dOdx(iNQ,iCar)),1)
*
         End Do
      End Do
      Call Free_Work(ip_Dbg)
      End If
#endif
      Call RotGrd(Work(ip_Crd),Work(ip_ZA),Work(ip_O),Work(ipdOdx),
     &            Dummy,nNQ,Do_Grad,.False.)
*
*     Distribute derivative of the principle axis system
*
      If (Do_Grad) Then
      iOff = ipdOdx
      Do iCar = 1, 3
         Do iNQ = 1, nNQ
#ifdef _DEBUG_
         Call RecPrt('dOdx','(3G20.10)',Work(iOff),3,3)
         Call RecPrt('dOdx(numerical)','(3G20.10)',
     &               Work(ip_dOdx(iNQ,iCar)),3,3)
#endif
            call dcopy_(9,Work(iOff),1,Work(ip_dOdx(iNQ,iCar)),1)
            iOff = iOff + 9
         End Do
      End Do
      End If
#ifdef _DEBUG_
*
*     Check translational invariance
*
      Call Allocate_Work(ip_debug,9)
      Do iCar = 1, 3
         Call FZero(Work(ip_debug),9)
         Do iNQ = 1, nNQ
            Call DaXpY_(9,1.0D0,Work(ip_dOdx(iNQ,iCar)),1,
     &                         Work(ip_debug),1)
         End Do
         Call RecPrt('Debug',' ',Work(ip_debug),3,3)
      End Do
      Call Free_Work(ip_debug)
*
#endif
*

      Call Free_Work(ipdOdx)
      Call Free_Work(ip_Crd)
      Call Free_Work(ip_ZA)
*                                                                      *
************************************************************************
*                                                                      *
      If (Rotational_Invariance.eq.Off) Then
         Call FZero(Work(ip_O),9)
         call dcopy_(3,[One],0,Work(ip_O),4)
         Do iNQ = 1, nNQ
            Do iCar = 1, 3
               Call FZero(Work(ip_dOdx(iNQ,iCar)),9)
            End Do
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Generate the angular grid
*
      Call Angular_grid()
*
      Crowding_tmp=Zero
      Do iNQ = 1, nNQ
*
*        Assign the angular grid to be used with each radial grid point
*
         nR_Eff=iWork(ip_nR_eff-1+iNQ)
         Call GetMem('ip_Angular','Allo','Inte',ip_A,nR_Eff)
         ip_iA=ip_of_iWork_d(Work(ip_Angular(iNQ)))
         iWork(ip_iA)=ip_A
         Call ICopy(nR_Eff,[nAngularGrids],0,iWork(ip_A),1)
*
*        Prune the angular grid
*
         If (Angular_Prunning.eq.On) Then
*
*
*---------- Find the R_min values of each angular shell
*
            lAng=Int(Work(ip_lMax(iNQ)))
            Do iAng = 0, lAng
               R_Min(iAng)=Zero
               ValExp=-One
               iSet=0
               Do iShell=1,nShell
                  iShll =iSD(0,iShell)
                  iAng_ =iSD(1,iShell)
                  NrExp =iSD(5,iShell)
*                 Write (6,*) 'iAng_,iAng=',iAng_,iAng
                  If (iAng_.eq.iAng.and.NrExp.ge.1) Then
                     Do iSym = 0, nSym-1
                        iNQ_=Maps2p(iShell,iSym)
*                       Write (6,*) 'iNQ_,iNQ=',iNQ_,iNQ
                        If (iNQ_.eq.iNQ) Then
                           ValExp=Shells(iShll)%Exp(NrExp)
                           iSet=1
                        End If
                     End Do
                  End If
               End Do
               If (ValExp.lt.Zero.and.iSet.eq.1) Then
                  Call WarningMessage(2,'ValExp.lt.Zero')
                  Call Abend()
               End If
               If (iSet.eq.1) Then
                  R_Min(iAng)=Eval_RMin(ValExp,iAng,Threshold)
                  If (iAng.eq.0) R_Min(iAng)=Zero
               End If
            End Do
*
            ip_iRx=ip_of_iWork_d(Work(ip_R_Quad(iNQ)))
            ip_Rx=iWork(ip_iRx)
            R_BS = Work(ip_R_RS(iNQ))
*
            If (iNQ.eq.iNQ_MBC) Then
               Crowding_tmp=Crowding
               Crowding=One + (Crowding-One)*0.25D0
               iReset=1
            End If
*
            Call Angular_Prune(Work(ip_Rx),nR_Eff,iWork(ip_A),Crowding,
     &                         Fade,R_BS,L_Quad,R_Min,lAng,
     &                         nAngularGrids,Info_Ang,LMax_NQ)
*
            If (iReset.eq.1) Then
               Crowding=Crowding_tmp
               iReset=0
            End If
*
         End if
*
      End Do
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      Write (6,*)
      Write (6,'(A)') ' =================================='
      Write (6,'(A)') ' =        Grid information        ='
      Write (6,'(A)') ' =================================='
      Write (6,'(A)') ' Legend             '
      Write (6,'(A)') ' ----------------------------------'
      Write (6,'(A)') ' ANr: element number'
      Write (6,'(A)') ' nR : number of radial grid points'
      Write (6,'(A)') ' iNQ: grid index'
      Write (6,'(A)') ' ----------------------------------'
      Write (6,*)
      Write (6,'(A)') ' iNQ ANr  nR'
      Do iNQ=1,nNQ
         iANr=Int(Work(ip_Atom_Nr(iNQ)))
         kR=iWork(ip_nR_Eff+iNQ-1)
         Write (6,'(3I4)') iNQ, iANr, kR
      End Do
      Write (6,*)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*-----Determine the spatial extension of the molecular system
*
      Box_Size=Two            ! Angstrom
      !Box_Size=1.0d0/Two            ! Angstrom
      Block_size=Box_Size
      x_min=1.0D99
      y_min=1.0D99
      z_min=1.0D99
      x_max=-1.0D99
      y_max=-1.0D99
      z_max=-1.0D99
      Do iAt = 0, nAtoms-1
         x_min=Min(x_min,Work(ipCoor+iAt*3+0))
         y_min=Min(y_min,Work(ipCoor+iAt*3+1))
         z_min=Min(z_min,Work(ipCoor+iAt*3+2))
         x_max=Max(x_max,Work(ipCoor+iAt*3+0))
         y_max=Max(y_max,Work(ipCoor+iAt*3+1))
         z_max=Max(z_max,Work(ipCoor+iAt*3+2))
      End Do
*
*---- Add half an box size around the whole molecule
*
      x_min=x_min-Box_Size/Two
      y_min=y_min-Box_Size/Two
      z_min=z_min-Box_Size/Two
      x_max=x_max+Box_Size/Two
      y_max=y_max+Box_Size/Two
      z_max=z_max+Box_Size/Two
*
*---- At least one finite box. Adjust to an even number of boxes.
*
      nx=Int((x_max-x_min+Box_Size)/Box_Size)
      nx=2*((nx+1)/2)
      ny=Int((y_max-y_min+Box_Size)/Box_Size)
      ny=2*((ny+1)/2)
      nz=Int((z_max-z_min+Box_Size)/Box_Size)
      nz=2*((nz+1)/2)
*
*---- Adjust extremal values to fit exactly with the
*     box size.
*
      dx=(DBLE(nx)*Box_Size-(x_max-x_min))/Two
      dy=(DBLE(ny)*Box_Size-(y_max-y_min))/Two
      dz=(DBLE(nz)*Box_Size-(z_max-z_min))/Two
*
      x_min=x_min-dx
      y_min=y_min-dy
      z_min=z_min-dz
      x_max=x_max+dx
      y_max=y_max+dy
      z_max=z_max+dz
*
*---- Add the infinite edge boxes
*
      nx=nx+2
      ny=ny+2
      nz=nz+2
#ifdef _DEBUG_
      Write (6,*) 'x_min=',x_min,dx
      Write (6,*) 'y_min=',y_min,dy
      Write (6,*) 'z_min=',z_min,dz
      Write (6,*) 'x_max=',x_max
      Write (6,*) 'y_max=',y_max
      Write (6,*) 'z_max=',z_max
      Write (6,*) 'nx,ny,nz=',nx,ny,nz
      Write (6,*) 'Total number of blocks=',nx*ny*nz
#endif
      number_of_subblocks=nx*ny*nz
*                                                                      *
************************************************************************
*                                                                      *
*     nFOrd: the order of the functional. nFOrd-1 is the number of times
*            the basis functions has to be differentiated to compute the
*            energy contribution.
*     mTmp: seems to be redundant?
*     mRad: number of different radial functions associated with a
*           basis function. This number depends on the type of
*           functional and the number of times the basis function has
*           to be differentiated in order to produce the values of the
*           parameters which the functional depends on (rho, grad rho,
*           and nabla rho).
*     nScr: Used to assemble integrals, (rho, grad rho, nabla rho)
*           (1,4, or 5). This is not needed in gradient calculations!
*     mAO: number of elements a basis function generates upon
*          differentiation (1,4,10,20, etc.)
*
      If (Functional_type.eq.LDA_type) Then
         nFOrd=1
         mAO=(nFOrd*(nFOrd+1)*(nFOrd+2))/6
         if(do_grad) mAO=4!AMS - GRADIENTS?
         If (.Not.Do_Grad) Then
            mTmp=1
            mRad=nFOrd
            nScr=nD*nAOMax
         Else
            mTmp=7
            mRad=nFOrd+1
            nScr=0
         End If
*
      Else If (Functional_type.eq.GGA_type) Then
         nFOrd=2
         mAO=(nFOrd*(nFOrd+1)*(nFOrd+2))/6
         if(do_grad) mAO=10
         If (.Not.Do_Grad) Then
            mTmp=7
            mRad=nFOrd
            nScr=nD*4*nAOMax
         Else
            mTmp=28
            mRad=nFOrd+1
            nScr=0
         End If
      Else If (Functional_type.eq.CASDFT_type) Then
*        I need to discuss this with Sergey!
*        nFOrd=3 !?
*        mAO=(nFOrd*(nFOrd+1)*(nFOrd+2))/6
         nFOrd=2
         mAO=10
         If (.Not.Do_Grad) Then
            mTmp=7
            mRad=nFOrd
            nScr=nD*4*nAOMax
         Else
            mTmp=28
            mRad=nFOrd+1
            nScr=0
         End If
*
      Else If (Functional_type.eq.meta_GGA_type1) Then
         nFOrd=2
         mAO=(nFOrd*(nFOrd+1)*(nFOrd+2))/6
         If (.Not.Do_Grad) Then
            mTmp=7
            mRad=nFOrd
            nScr=nD*4*nAOMax
         Else
            mTmp=28
            mRad=NFOrd+1
            nScr=0
         End If
*
      Else If (Functional_type.eq.meta_GGA_type2) Then
         nFOrd=3
         mAO=(nFOrd*(nFOrd+1)*(nFOrd+2))/6
         If (.Not.Do_Grad) Then
            mTmp=7
            mRad=nFOrd
            nScr=nD*5*nAOMax
         Else
            mTmp=28
            mRad=NFOrd+1
            nScr=0
         End If
*
*     Else If (Functional_type.eq.Other_type) Then
      Else
         mTmp=0 ! Dummy initialize
         mRad=0 ! Dummy initialize
         mAO=0  ! Dummy initialize
         nScr=0 ! Dummy initialize
         Call WarningMessage(2,'Functional_type.eq.Other_type')
         Call Abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate scratch for processing AO's on the grid
*
      nTmp=nGridMax*Max(mTmp,nScr)
      Call GetMem('Tmp','Allo','Real',ipTmp,nTmp+nShell*nIrrep)

      ipTabAOMax=ipTmp+nTmp
*
      nMem=0
      nSO =0
      Do ish = 1, nShell
         iAng  = iSD( 1,iSh)
         iCmp  = iSD( 2,iSh)
         iBas  = iSD( 3,iSh)
         iPrim = iSD( 5,iSh)
*
         nxyz    = nGridMax*3*(iAng+mRad)
         nDrv     = mRad - 1
         nForm    = 0
         Do iDrv  = 0, nDrv
            nForm = nForm + nElem(iDrv)
         End Do
         nTerm    = 2**nDrv
         nAngular = 5*nForm*nTerm
         nRad    = iPrim*nGridMax*mRad
         nRadial = iBas*nGridMax*mRad
         If (On_Top) Then
            mdci  = iSD(10,iSh)
            kAO=iCmp*iBas*nGridMax
            nSO=kAO*nSym/nStab(mdci)*mAO
         End If
c         nMem=Max(nMem,nxyz+nAngular+nRad+nRadial+nSO)
         nMem=Max(nMem,nxyz+nAngular+nRad+nRadial,2*nSO)
      End Do
*
      Call GetMem('nMem','Allo','Real',ipMem,nMem)
*                                                                      *
************************************************************************
*                                                                      *
*     Access the file with Grid points and weights.
*
*---- Open the file.
      Lu_Grid=88
      Call DaName_MF_WA(Lu_Grid,'NQGRID')
*
      If (iGrid_Set.eq.Not_Specified) iGrid_Set=Final
*
*---- Read the status flag.
      iDisk_Grid=0
      Call iDaFile(Lu_Grid,2,G_S,5,iDisk_Grid)
*
      Grid_Status=G_S(iGrid_Set)
      If (Old_Functional_Type.ne.Functional_Type) Then
         G_S(Final)=Regenerate
         G_S(Intermediate)=Regenerate
         Grid_Status=Regenerate
      End If
      iDisk_Grid=iDisk_Set(iGrid_Set)
*
*---- Allocate memory for the master TOC.
      Call GetMem('GridInfo','Allo','Inte',ip_GridInfo,
     &            2*number_of_subblocks)
*
*---- Retrieve the TOC or regenerate it.
*
*     The table contains two data items per subblock.
*     1) disk address and 2) number of batches.
*
      If (Grid_Status.eq.Regenerate) Then
C        Write (6,*) 'Grid_Status.eq.Regenerate'
         Grid_Status=Regenerate
         Call ICopy(2*number_of_subblocks,[0],0,iWork(ip_GridInfo),1)
         Call iDaFile(Lu_Grid,1,iWork(ip_GridInfo),
     &                2*number_of_subblocks,iDisk_Grid)
         Old_Functional_Type=Functional_Type
      Else If (Grid_Status.eq.Use_Old) Then
C        Write (6,*) 'Grid_Status.eq.Use_Old'
         Call iDaFile(Lu_Grid,2,iWork(ip_GridInfo),
     &                2*number_of_subblocks,iDisk_Grid)
      Else
         Call WarningMessage(2,'Illegal Grid Status!')
         Call Abend()
      End If
*
      Call ParmPkR8(Pck_Old,PMode_old)
      Call IniPkR8(T_X,.True.)
*                                                                      *
************************************************************************
*                                                                      *
*     Setup some symmetry stuff outside the loop
*
      ndc = 0
      Do iSh = 1, nShell
         ndc = Max(ndc,iSD(10,iSh))
      End Do
      ndc2 = ndc**2
      Call GetMem('Fact','Allo','Real',ip_Fact,ndc2)
      Do mdci = 1, ndc
         nDegi=nIrrep/nStab(mdci)
         Do mdcj = 1, ndc
            nDegj=nIrrep/nStab(mdcj)
*
            Call DCR(LmbdR,iOper,nIrrep,jStab(0,mdci),
     &               nStab(mdci),jStab(0,mdcj),
     &               nStab(mdcj),iDCRR,nDCRR)
*
            iuv = nStab(mdci)*nStab(mdcj)
            If (MolWgh.eq.1) Then
               Fact = DBLE(nIrrep) / DBLE(LmbdR)
            Else If (MolWgh.eq.0) Then
               Fact = DBLE(iuv) / DBLE(nIrrep * LmbdR)
            Else
               Fact = Sqrt(DBLE(iuv))/ DBLE(LmbdR)
            End If
            Fact=Fact*DBLE(nDCRR)/DBLE(nDegi*nDegj)
*
*---------- Save: Fact
*
            ij = (mdcj-1)*ndc + mdci
            Work(ip_Fact+ij-1) = Fact
*
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
C     Call QExit('Setup_NQ')
      Return
      End
