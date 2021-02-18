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
* Copyright (C) 1995, Niclas Forsberg                                  *
************************************************************************
C!-----------------------------------------------------------------------!
C!
      Subroutine RotTranRem(Sinv,S,Mass,AtCoord,NumOfAt,NumInt)
C!
C!  Purpose:
C!    Project total translation and total rotation out of inverted
C!    S matrix.
C!
C!  Input:
C!    S        : Real*8 three dimensional array.
C!    Mass     : Real*8 array - masses of the atoms.
C!    AtCoord  : Real*8 two dimensional array - coordinates
C!               of the atoms.
C!
C!  Output:
C!    Sinv     : Real*8 three dimensional array - Inverted
C!               S matrix with rotation and translation projected out.
C!
C!  Calls:
C!    Daxpy  (ESSL)
C!    Dcopy  (ESSL)
C!    Ddot_  (ESSL)
C!
C!  Uses:
C!    LinAlg
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1995.
C!
c       Use Linalg
c       Implicit None
#include "Constants_mula.fh"
      Integer NumOfAt,NumInt,nFree
      Real*8 AtCoord (3,NumofAt)
      Real*8 S (3, NumOfAt, NumInt)
      Real*8 Sinv (3, NumOfAt, NumInt)
      Real*8 Mass (NumOfAt)
      Real*8 Det,X
      Integer iAtom,n,k,i,j
#include "WrkSpc.fh"
C!
C!---- Initialize.
      nFree = 6
      n=3*NumOfAt
      Call GetMem('Amat','Allo','Real',ipAmat,3*NumOfAt*nFree)
      lAmat=n
c       Amat = 0.0d0
      call dcopy_(3*NumOfAt*nFree,[0.0d0],0,Work(ipAmat),1)
C!
C!---- Invert S matrix.
      Call GetMem('Temp2','Allo','Real',ipTemp2,NumInt*Numint)

      Call DGEMM_('T','N',
     &            NumInt,NumInt,3*NumOfAt,
     &            1.0d0,S,3*NumOfAt,
     &            S,3*NumOfAt,
     &            0.0d0,Work(ipTemp2),NumInt)
C! Invert, by solving eq Temp2*X=Temp1. Solution computed in-place.
C! Thus, Temp1=Unit matrix(in) and contains solution (out).
      Call GetMem('Temp1','Allo','Real',ipTemp1,NumInt*Numint)

      call dcopy_(NumInt**2,[0.0d0],0,Work(ipTemp1),1)
      call dcopy_(NumInt,[1.0d0],0,Work(ipTemp1),NumInt+1)
      Call Dool_MULA(Work(ipTemp2),NumInt,NumInt,Work(ipTemp1),
     &  NumInt,NumInt,det)
C!PAM01 Replacement for Dool_MULA, if superstable solution is wanted:
C! Eps=1.0D-8
C! Call SymSolve(Temp2,Temp1,Temp3,Temp4,Eps)
C! call dcopy_(NumInt**2,Temp3,1,Temp1,1)
C!PAM01 End of replacement code.

      Call DGEMM_('N','N',
     &            3*NumOfAt,NumInt,NumInt,
     &            1.0d0,S,3*NumOfAt,
     &            Work(ipTemp1),NumInt,
     &            0.0d0,Sinv,3*NumOfAt)
      Call GetMem('Temp1','Free','Real',ipTemp1,NumInt*Numint)
      Call GetMem('Temp2','Free','Real',ipTemp2,NumInt*Numint)
C!
C!---- Pure translation.
      k = 1
      Do iAtom = 1,NumOfAt
      Work(ipAmat+k-1) = 1.0d0
      Work(ipAmat+k+lAmat) = 1.0d0
      Work(ipAmat+k+1+lAmat*2) = 1.0d0
      k = k+3
      End Do
C!
C!---- Pure rotation.
      k = 1
      Do iAtom = 1,NumOfAt
      Work(ipAmat+k  +lAmat*3-1) =-AtCoord(2,iAtom)
      Work(ipAmat+k+1+lAmat*3-1) = AtCoord(1,iAtom)
      Work(ipAmat+k+2+lAmat*3-1) = 0.0d0
      Work(ipAmat+k + lAmat*4-1) = AtCoord(3,iAtom)
      Work(ipAmat+k+1+lAmat*4-1) = 0.0d0
      Work(ipAmat+k+2+lAmat*4-1) =-AtCoord(1,iAtom)
      Work(ipAmat+k  +lAmat*5-1) = 0.0d0
      Work(ipAmat+k+1+lAmat*5-1) =-AtCoord(3,iAtom)
      Work(ipAmat+k+2+lAmat*5-1) = AtCoord(2,iAtom)
      k = k+3
      End Do
C!
C!---- Scale with the mass of the atom.
      n=3*NumOfAt
      Call GetMem('AmatMass','Allo','Real',
     &  ipAmatMass,3*NumOfAt*nFree)
c       AmatMass = Amat
      call dcopy_(lAmat*nFree,Work(ipAmat),1,Work(ipAmatMass),1)
* PAM04: Replace the following code section...
*       Do i = 1,nFree
*          Acol => AmatMass(:,i)
*          k = 1
*          Do j = 1,NumOfAt
*             jMass = (k+2)/3
*             Acol(k  ) = uToAu*Mass(jMass)*Acol(k  )
*             Acol(k+1) = uToAu*Mass(jMass)*Acol(k+1)
*             Acol(k+2) = uToAu*Mass(jMass)*Acol(k+2)
*             k = k+3
*          End Do
*       End Do
* PAM04: ... with the following...
      Do i=1,nFree
      Do j=1,NumOfAt
      X=uToAu*Mass(j)
      Work(ipAmatMass+1+3*(j-1)+lAmat*(i-1)-1)=
     &    X*Work(ipAmatMass+1+3*(j-1)+lAmat*(i-1)-1)
      Work(ipAmatMass+2+3*(j-1)+lAmat*(i-1)-1)=
     &    X*Work(ipAmatMass+2+3*(j-1)+lAmat*(i-1)-1)
      Work(ipAmatMass+3+3*(j-1)+lAmat*(i-1)-1)=
     &    X*Work(ipAmatMass+3+3*(j-1)+lAmat*(i-1)-1)
      End Do
      End Do
* PAM04: ... until here.
C!
C!---- Project rotation and translation out of S matrix.
      Call GetMem('Temp1','Allo','Real',ipTemp1,nFree*nFree)
      Call DGEMM_('T','N',
     &            nFree,nFree,3*NumOfAt,
     &            1.0d0,Work(ipAmatMass),3*NumOfAt,
     &            Work(ipAmat),3*NumOfAt,
     &            0.0d0,Work(ipTemp1),nFree)
      Call GetMem('Ainv','Allo','Real',ipAinv,nFree*nFree)
      call dcopy_(nFree**2,[0.0d0],0,Work(ipAinv),1)
      call dcopy_(nFree,[1.0d0],0,Work(ipAinv),nFree+1)
      Call Dool_MULA(Work(ipTemp1),nFree,nFree,Work(ipAinv),nFree,
     &  nFree,det)
      Call GetMem('Temp1','Free','Real',ipTemp1,nFree*nFree)
C!
      Call GetMem('Temp2','Allo','Real',ipTemp2,nFree*Numint)
      Call DGEMM_('T','N',
     &            nFree,NumInt,3*NumOfAt,
     &            1.0d0,Work(ipAmatMass),3*NumOfAt,
     &            Sinv,3*NumOfAt,
     &            0.0d0,Work(ipTemp2),nFree)
      Call GetMem('AmatMass','Free','Real',
     &  ipAmatMass,3*NumOfAt*nFree)
C!
      Call GetMem('Temp3','Allo','Real',ipTemp3,nFree*Numint)

      Call DGEMM_('N','N',
     &            nFree,NumInt,nFree,
     &            1.0d0,Work(ipAinv),nFree,
     &            Work(ipTemp2),nFree,
     &            0.0d0,Work(ipTemp3),nFree)
      Call GetMem('Ainv','Free','Real',ipAinv,nFree*nFree)
      Call GetMem('Temp2','Free','Real',ipTemp2,nFree*Numint)
C!
      Call GetMem('Stemp','Allo','Real',ipStemp,3*NumOfAt*Numint)
      Call DGEMM_('N','N',
     &            3*NumOfAt,NumInt,nFree,
     &            1.0d0,Work(ipAmat),3*NumOfAt,
     &            Work(ipTemp3),nFree,
     &            0.0d0,work(ipStemp),3*NumOfAt)
      Call GetMem('Amat','Free','Real',ipAmat,3*NumOfAt*nFree)
      Call GetMem('Temp3','Free','Real',ipTemp3,nFree*Numint)

C!
c       Sinv = Sinv-Stemp
      call daxpy_(3*NumOfAt*Numint,-1.0d0,Work(ipStemp),1,Sinv,1)
C!
      Call GetMem('Stemp','Free','Real',ipStemp,3*NumOfAt*Numint)
C!
      End
