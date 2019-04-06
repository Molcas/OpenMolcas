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
* Copyright (C) Anders Bernhardsson                                    *
************************************************************************
      Subroutine PutRlx(D,DS,P,DAO,C)
      Implicit Real*8 (a-h,o-z)
#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
#include "output_ras.fh"
#include "WrkSpc.fh"
      Parameter (ROUTINE='PUTRLX  ')
      Real*8 D(*),DS(*),P(*),DAO(*),C(*)
      Dimension rdum(1)

      IPRLEV=IPRLOC(3)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF
*
*     Read in corresponding density matrixes
*

      jDisk = IADR15(3)
      NZ=NTOT2
      Do i=1,iRlxRoot-1
        Call DDaFile(JOBIPH,0,rdum,NACPAR,jDisk)
        Call DDaFile(JOBIPH,0,rdum,NACPAR,jDisk)
        Call DDaFile(JOBIPH,0,rdum,NACPR2,jDisk)
        Call DDaFile(JOBIPH,0,rdum,NACPR2,jDisk)
      End Do
      Call DDaFile(JOBIPH,2,D,NACPAR,jDisk)
      Call DDaFile(JOBIPH,2,DS,NACPAR,jDisk)
      Call DDaFile(JOBIPH,2,P,NACPR2,jDisk)
      Call DDaFile(JOBIPH,0,rdum,NACPR2,jDisk)
*
*     Construc D-ACTIVE AND D-INACTIVE IN AO BASIS
*
      Call GetMem('TEMP','ALLO','REAL',ipDA,NZ)
      Call GetMem('TEMP','ALLO','REAL',ipDI,NZ)
      Call GetMem('TEMP','ALLO','REAL',ipDS,NZ)
*
      Call Get_D1I_RASSCF(C,Work(ipDI))
*
      Call GetMem('TEMP','ALLO','REAL',ipD,NACPAR)
*
      CALL DCOPY_(NACPAR,DS,1,Work(ipD),1)
      Call DBLOCK(Work(ipD))
      CALL Get_D1A_RASSCF(C,Work(ipD),Work(ipDS))
*
      CALL DCOPY_(NACPAR,D,1,Work(ipD),1)
      Call DBLOCK(Work(ipD))
      CALL Get_D1A_RASSCF(C,Work(ipd),Work(ipDA))
*
*     Construct the Fock matrix used for the connection term.
*
      NFSIZE=MAX(NACPAR,NTOT4)
      Call GetMem('TEMP','ALLO','REAL',ipF,NFSIZE)
      Call GetMem('TEMP','ALLO','REAL',ipB,NZ)
      Call GetMem('TEMP','ALLO','REAL',ipQ,NZ)
      Call GetMem('TEMP','ALLO','REAL',ipFA,NZ)
      Call GetMem('TEMP','ALLO','REAL',ipFI,NZ)
      Call GetMem('TEMP','ALLO','REAL',ipPUVX,NFINT)
      Call GetMem('TEMP','ALLO','REAL',iptuvx,nacpr2)
      IPR=0
      IF(IPRLEV.EQ.3) IPR=1
      IF(IPRLEV.EQ.4) IPR=5
      IF(IPRLEV.EQ.5) IPR=10
      Call FZero(Work(ipPUVX),NFINT)
      Call TraCtl2(C,Work(ipPUVX),Work(ipTUVX),
     &             Work(ipDI),Work(ipFI),
     &             Work(ipDA),Work(ipFA),ipr,lsquare,ExFac)
      CALL SGFCIN(C,WORK(ipF),Work(ipFI),Work(ipDI),Work(ipDA),
     &                                              Work(ipDS))
      call dcopy_(ntot4,[0.0d0],0,Work(ipF),1)
      call dcopy_(ntot4,[0.0d0],0,Work(ipB),1)
*
*     Prevent FMAT from changing Active fock matrix
*
      iTmp=newfock
      newFock=-99999
*
      Call Fmat(C,Work(ipPUVX),Work(ipD),Work(ipDA),Work(ipFI),
     &          Work(ipFA))
      newFock=itmp
      IFINAL=1
      rTmp=   CBLBM
      itmp= IBLBM
      jtmp= JBLBM
      istmp= ISYMBB

      IF(ISTORP(NSYM+1).GT.0) THEN
          CALL GETMEM('ISTRP','ALLO','REAL',ipP,ISTORP(NSYM+1))
          CALL PMAT_RASSCF(P,WORK(ipP))
      Else
          ipP = ip_Dummy
      END IF
*
      Call FOCK(Work(ipF),Work(ipB),Work(ipFI),
     &          Work(ipFA),Work(ipd),WORK(ipP),Work(ipQ),Work(ipPUVX),
     &          IFINAL,C)
      CBLBM=rtmp
      IBLBM=itmp
      JBLBM=jtmp
      ISYMBB=istmp
      IF(ISTORP(NSYM+1).GT.0)
     &    CALL GETMEM('ISTRP','FREE','REAL',ipP,ISTORP(NSYM+1))

      RlxGrd=DNRM2_(NSXS,Work(ipB),1)
*
      Call GetMem('TEMP','FREE','REAL',ipD,NACPAR)
      Call GetMem('TEMP','FREE','REAL',iptuvx,idum)
      Call GetMem('TEMP','FREE','REAL',ippuvx,idum)
      Call GetMem('TEMP','FREE','REAL',ipFI,idum)
      Call GetMem('TEMP','FREE','REAL',ipFA,idum)
      Call GetMem('TEMP','FREE','REAL',ipQ,idum)
      Call GetMem('TEMP','FREE','REAL',ipB,idum)
      Call GetMem('TEMP','FREE','REAL',ipF,NFSIZE)
*
* Add up one electron densities
*
      call daxpy_(nZ,1.0d0,Work(ipDA),1,Work(ipDI),1)
      Call Fold(nSym,nBas,Work(ipDI),DAO)

      Call GetMem('TEMP','FREE','REAL',ipDS,NZ)
      Call GetMem('TEMP','FREE','REAL',ipDA,idum)
      Call GetMem('TEMP','FREE','REAL',ipDI,idum)

      return
      end
