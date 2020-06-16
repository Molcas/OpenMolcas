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
      Subroutine MkGrid(natom,ipCord,ipGrd,nGrdPt,iRMax,DeltaR,
     &                  Forces,ipIsMM,iGrdTyp,ipDGrd,nAtQM)
      use PCM_arrays
      Implicit Real*8 (A-H,O-Z)
*
#include "espf.fh"
*
#include "rctfld.fh"
#include "stdalloc.fh"
      Logical Forces,Process,Dirty
*
      Call QEnter('mkgrid')
      iPL = iPL_espf()
*
      iPrint = 5
      If (iPL.ge.3) iPrint = 50
      If (iPL.ge.4) iPrint = 99
      nDiff = 0
      Call GetMem('Atomic Numbers','Allo','Inte',ipAN,natom)
      Call GetMem('Get_Atoms','Allo','Real',ipChrg,natom)
      Call Get_dArray('Nuclear charge',Work(ipChrg),natom)
      Do iatom = 0, natom-1
         iWork(ipAN+iatom) = Int(Work(ipChrg+iatom))
      End Do
      Call GetMem('Get_Atoms','Free','Real',ipChrg,natom)
      nGrdPt_old = nGrdPt
*
*     PNT grid (it uses Angstroms !!!)
*
      If (Abs(iGrdTyp).eq.1) Then
         DeltaR = DeltaR * Angstrom
         Process = (iGrdTyp.eq.1)
         Call DScal_(3*natom,Angstrom,Work(ipCord),1)
         Call PNT(6,natom,Work(ipCord),iRMax,DeltaR,iWork(ipAN),
     &            nGrdPt,Work(ipGrd),iWork(ipIsMM),Process)
         Call DScal_(3*natom,One/Angstrom,Work(ipCord),1)
         DeltaR = DeltaR / Angstrom
         If (nGrdPt.le.0) Then
            Write(6,'(A)') ' Error in espf/mkgrid: nGrdPt < 0 !!!'
            Call Quit_OnUserError()
         End If
*
*        Printing the PNT point coordinates
*
         If(Process .and. .not. DoDeriv) Then
            If (iPL.ge.4) Then
               Write(6,'(A,I5,A)') ' PNT grid (in Angstrom) '
               Do iPt = 1, nGrdPt
                  iCur = 3*(iPt-1)
                  Write(6,'(A4,3F15.6)') ' X  ',Work(ipGrd+iCur),
     &                                          Work(ipGrd+iCur+1),
     &                                          Work(ipGrd+iCur+2)
               End Do
            End If
            Call DScal_(3*nGrdPt,One/Angstrom,Work(ipGrd),1)
         End If
c
c     GEPOL grid (made of iRMax surfaces)
c
      Else
         iXPolType = 0
         PCM = .True.
         DoDeriv = Forces
         nDer = nAtQM*3
         Do J = 0, iRMax-1
            If (iPL.ge.3) Write(6,'(A13,I1)') ' GEPOL shell ',J+1
            Call GetMem('LcCoor','Allo','Real',ip_LcCoor,3*natom)
            Call GetMem('LcANr','Allo','Inte',ip_LcANr,natom)
            nPCM_info = 0
            Call PCM_Cavity(iPrint,0,natom,Angstrom,Work(ipCord),
     &                      iWork(ipAN),iWork(ipIsMM),Work(ip_LcCoor),
     &                      iWork(ip_LcANr),J)
            Call GetMem('LcANr','Free','Inte',ip_LcANr,natom)
            Call GetMem('LcCoor','Free','Real',ip_LcCoor,3*natom)
            If(J.eq.0) Then
              nTmp = 0
              nGrdPt = nTs
              Call Allocate_Work(ipTmp,3*nGrdPt)
            Else
              nTmp = nGrdPt
              Call Allocate_Work(ipBla,3*nTmp)
              call dcopy_(3*nTmp,Work(ipTmp),1,Work(ipBla),1)
              nGrdPt = nTs + nGrdPt
              Call Free_Work(ipTmp)
              Call Allocate_Work(ipTmp,3*nGrdPt)
              call dcopy_(3*nTmp,Work(ipBla),1,Work(ipTmp),1)
              Call Free_Work(ipBla)
            End If
            If (DoDeriv) Call GetMem('ESPF_DGrid','Allo','Real',ipDGrd,
     &                               3*nGrdPt*nDer)
            Do I = 1, nTs
               call dcopy_(3,Work(ip_Tess + 4*(I-1)),1,
     &                      Work(ipTmp   + 3*nTmp+3*(I-1)),1)
            End Do
            If (DoDeriv) call dcopy_(3*nGrdPt*nDer,DPnt,1,
     &                                            Work(ipDGrd ),1)
            Call GetMem('PCMSph','Free','Real',ip_Sph,nPCM_info)
            If (DoDeriv) Then
               LcNAtm = ISlPar(42)
               NDeg = 3*LcNAtm
               Call mma_deallocate(dTes)
               Call mma_deallocate(dPnt)
               Call mma_deallocate(dRad)
               Call GetMem('DerCentr','Free','Real',ip_DCntr,3*nS*NDeg)
               Call GetMem('PCM-Q','Free','Real',ip_Q,2*nTs)
            End If
         End Do
         Call GetMem('ESPF_Grid','Allo','Real',ipGrd,3*nGrdPt)
         call dcopy_(3*nGrdPt,Work(ipTmp),1,Work(ipGrd),1)
         Call Free_Work(ipTmp)
*
*        Cleaning the GEPOL grid:
*        all grid points must be distant by 1 bohr at least
*
10       Dirty = .False.
         Call Allocate_iWork(ipKeep,nGrdPt)
         Do iPnt = 0, nGrdPt-1
            iWork(ipKeep+iPnt) = 1
         End Do
         Do iPnt = 0, nGrdPt-2
            If (iWork(ipKeep+iPnt) .eq. 0) Goto 11
            Do jPnt = iPnt+1, nGrdPt-1
               X = Work(ipGrd+3*jPnt  ) - Work(ipGrd+3*iPnt  )
               Y = Work(ipGrd+3*jPnt+1) - Work(ipGrd+3*iPnt+1)
               Z = Work(ipGrd+3*jPnt+2) - Work(ipGrd+3*iPnt+2)
               R = Sqrt(X*X+Y*Y+Z*Z)
               If (R .lt. One) iWork(ipKeep+jPnt) = 0
            End Do
11          Continue
         End Do
         New_nGrdPt = 0
         Do iPnt = 0, nGrdPt-1
            If(iWork(ipKeep+iPnt) .eq. 1) New_nGrdPt = New_nGrdPt + 1
         End Do
         Dirty = New_nGrdPt.lt.nGrdPt
         If(Dirty) Then
            Call Allocate_Work(ipTmpG,3*nGrdPt)
            call dcopy_(3*nGrdPt,Work(ipGrd),1,Work(ipTmpG),1)
            Call GetMem('ESPF_Grid','Free','Real',ipGrd,3*nGrdPt)
            Call GetMem('ESPF_Grid','Allo','Real',ipGrd,3*New_nGrdPt)
            If (DoDeriv) Then
               Call Allocate_Work(ipTmpDG,3*nGrdPt*NDer)
               call dcopy_(3*nGrdPt*NDer,Work(ipDGrd),1,Work(ipTmpDG),1)
               Call GetMem('ESPF_DGrid','Free','Real',ipDGrd,
     &                  3*nGrdPt*NDer)
               Call GetMem('ESPF_DGrid','Allo','Real',ipDGrd,
     &                  3*New_nGrdPt*NDer)
            End If
            ibla = -1
            Do iPnt = 0, nGrdPt-1
               If (iWork(ipKeep+iPnt) .eq. 1) Then
                  ibla = ibla + 1
                  call dcopy_(3,Work(ipTmpG+3*iPnt),1,
     &                         Work(ipGrd +3*ibla),1)
                  If (DoDeriv)
     &               call dcopy_(9*nAtQM,Work(ipTmpDG+iPnt),nGrdPt,
     &                                  Work(ipDGrd +ibla),New_nGrdPt)
               End If
            End Do
            If (DoDeriv) Call Free_Work(ipTmpDG)
            Call Free_Work(ipTmpG)
            nGrdPt = New_nGrdPt
         End If
         Call Free_iWork(ipKeep)
         If(Dirty) Goto 10
*
*        Printing the GEPOL point coordinates
*
         If (.not.DoDeriv .and. iPL.ge.4) Then
            Call DScal_(3*nGrdPt,Angstrom,Work(ipGrd),1)
            Write(6,'(A)') 'PCM grid (in Angstroms):'
            Do iPnt = 0, nGrdPt-1
               iCur = 3*iPnt
               Write(6,'(A4,3F15.6)') ' X  ',Work(ipGrd+iCur  ),
     &                                       Work(ipGrd+iCur+1),
     &                                       Work(ipGrd+iCur+2)
            End Do
            Call DScal_(3*nGrdPt,One/Angstrom,Work(ipGrd),1)
         End If
         PCM = .False.
      End If
      If(nGrdPt_old.ne.0 .and. nGrdPt.ne.nGrdPt_old) Then
        Write(6,'(A,2i10)') 'MkGrid: inconsistency in nGrdPt:',
     &                      nGrdPt_old,nGrdPt
        Call Abend()
      End If
      Call GetMem('Atomic Numbers','Free','Inte',ipAN,natom)
*
      Call QExit('mkgrid')
      Return
      End
