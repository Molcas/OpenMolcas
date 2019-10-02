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
* Copyright (C) 1990,1996, Markus P. Fuelscher                         *
************************************************************************
      SUBROUTINE EXPLH2(DIAG,ONEINT,TUVX,ISEL,EXPLE,EXPLV)
************************************************************************
*                                                                      *
*     Compute and diagonalize the explicit Hamiltonian in a            *
*     subspace smaller or identical to nSel. nSel may be chosen        *
*     somewhat smaller in order to avoid selecting only one of         *
*     several degenerate diagonal matrix elements.                     *
*                                                                      *
*     calling arguments:                                               *
*     DIAG    : array of Real*8                                        *
*               diagonal Hamiltonian                                   *
*     ONEINT  : array of Real*8                                        *
*               one-electron integrals                                 *
*     TUVX    : array of Real*8                                        *
*               two-electron integrals                                 *
*     ISEL    : array of integer                                       *
*               index array                                            *
*     EXPLE   : array of Real*8                                        *
*               eigenvalues of the explicit Hamiltonian                *
*     EXPLV   : array of Real*8                                        *
*               eigenvectors of the explicit Hamiltonian               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher and J. Olsen                                      *
*     University of Lund, Sweden, 1990                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*     - updated for integral direct and reaction field calculations    *
*       M.P. Fuelscher, University of Lund, Sweden, 1996               *
*                                                                      *
************************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "ciinfo.fh"
#include "spinfo.fh"
#include "csfbas.fh"
#include "strnum.fh"
#include "WrkSpc.fh"
#include "timers.fh"
#include "output_ras.fh"
      PARAMETER (ROUTINE='EXPLH2  ')
      DIMENSION DIAG(*)
      DIMENSION ONEINT(*)
      DIMENSION TUVX(*)
      DIMENSION EXPLE(*)
      DIMENSION EXPLV(*)
      DIMENSION ISEL(*)
*
      CALL QENTER('EXPLH')
      Call Timing(Omega_1,Swatch,Swatch,Swatch)
      IPRLEV=IPRLOC(3)

      ECORE=0.0D0
      MXXSEL=NSEL
      NHEX=NSEL*(NSEL+1)/2
*
* ALLOCATE LOCAL MEMORY
*
      CALL GETMEM('IPCNF','ALLO','INTE',LW1,NCNASM(LSYM))
      CALL GETMEM('HONE','ALLO','REAL',LOCONE,NAC**2)
      CALL GETMEM('EXHAM','ALLO','REAL',LEXHAM,NHEX)
*
* EXPAND ONE-INTS FROM TRIANGULAR PACKING TO FULL STORAGE MODE
*
      CALL TRIEXP(ONEINT,Work(LOCONE),NAC)
*
* Load the diagonal approximation of the CI Hamiltonian
*
      Call Load_H_diag(nConf,DIAG,LuDavid)
*
* CONSTRUCT THE EXPLICIT HAMILTONIAN
*
      IPRINT=0
      IF(IPRLEV.EQ.INSANE) IPRINT=40
      CALL GETMEM('IREOTS','ALLO','INTEGER',IREOTS,NAC)
      CALL GETMEM('EXHSCR','MAX','REAL',LW2,MXXWS)
      CALL GETMEM('EXHSCR','ALLO','REAL',LW2,MXXWS)
      CALL GET_IREOTS(IWORK(IREOTS),NAC)
      CALL PHPCSF(Work(LEXHAM),ISEL,iWork(LW1),
     &            MXXSEL,Work(KDTOC),
     &            iWork(KDFTP),iWork(KICONF(1)),
     &            LSYM,Work(LOCONE),ECORE,NAC,
     &            Work(LW2),NCNASM(LSYM),
     &            (NAEL+NBEL),NAEL,NBEL,
     &            NSEL,NPCNF,DIAG,TUVX,IPRINT,ExFac,IWORK(IREOTS))
c      IF ( IPRLEV.EQ.INSANE) then
        Call Square(Work(LEXHAM),EXPLV,1,NSEL,NSEL)
        CALL RECPRT('Square Explicit Hamiltonian',' ',EXPLV,NSEL,NSEL)
c      END IF
      CALL GETMEM('EXHSCR','FREE','REAL',LW2,MXXWS)
      CALL GETMEM('IREOTS','FREE','INTEGER',IREOTS,NAC)
      CALL GETMEM('HONE','FREE','REAL',LOCONE,NAC**2)
      CALL GETMEM('IPCNF','FREE','INTE',LW1,NCNASM(LSYM))
*
* DIAGONALIZE THE EXPLICIT HAMILTONIAN.
*
*     If ( nSel.eq.nConf ) then
      If ( .true. ) then
        CALL DCOPY_(MXXSEL*MXXSEL,[0.0D0],0,EXPLV,1)
        DO 10 I=1,NSEL
          II=I+NSEL*(I-1)
          EXPLV(II)=1.0D00
10      CONTINUE
*       CALL Jacob(Work(LEXHAM),EXPLV,NSEL,NSEL)
#ifdef _DEBUG_
        CALL NIdiag(Work(LEXHAM),EXPLV,NSEL,NSEL,0)
#else
        CALL NIdiag_new(Work(LEXHAM),EXPLV,NSEL,NSEL,0)
#endif
        CALL JACORD(Work(LEXHAM),EXPLV,NSEL,NSEL)
        DO 15 I=1,NSEL
          EXPLE(I)=Work(LEXHAM-1+I*(I+1)/2)
15      CONTINUE
      Else
        Call GetMem('ExHscr','Allo','Real',lwscr,nSel)
        Call Square(Work(LEXHAM),EXPLV,1,NSEL,NSEL)
        Call Eigen_Molcas(NSEL,EXPLV,EXPLE,Work(lwscr))
        Call GetMem('ExHscr','Free','Real',lwscr,nSel)
      End If
      CALL GETMEM('EXHAM','FREE','REAL',LEXHAM,NHEX)
      IF( IPRLEV.GE.INSANE)
     &  CALL IVCPRT
     &    ('Configurations included in the explicit Hamiltonian',' ',
     &      ISEL,NSEL)
      IF( IPRLEV.GE.INSANE)
     &  CALL DVCPRT('Eigenvalues of the explicit Hamiltonian',' ',
     &               EXPLE,NSEL)
      IF( IPRLEV.GE.INSANE)
     &  CALL RECPRT('Eigenvectors of the explicit Hamiltonian',' ',
     &               EXPLV,NSEL,NSEL)
*
      Call Timing(Omega_2,Swatch,Swatch,Swatch)
      Omega_2 = Omega_2 - Omega_1
      Omega_3 = Omega_3 + Omega_2
      CALL QEXIT('EXPLH')

      RETURN
      END
