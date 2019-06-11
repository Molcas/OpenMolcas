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
      Subroutine hdir( nDir, nDirZee, dirX, dirY, dirZ, dir_weight,
     &                 nP, nsymm, ngrid, nDirTot, dHX, dHY, dHZ, dHW)
c      this routine generates the directions of the applied magnetic
c      field according to Lebedev-Laikov grids using the the given parameters (nsymm, ngrid)

      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer                      :: nP,nDirTot,nDir,nDirZee,i,j
      Integer                      :: nsymm,ngrid

      Real(kind=wp) :: dirX(nDir), dirY(nDir), dirZ(nDir),
     &                 dir_weight(nDirZee,3)
      Real(kind=wp) :: dHX(nDirTot), dHY(nDirTot), dHZ(nDirTot),
     &                 dHW(nDirTot)
      Real(kind=wp) :: X(nP), Y(nP), Z(nP), W(nP)

      Call qEnter('hdir')

      If ((nDirTot-nDir-nDirZee-nP).ne.0) Then
        Write(6,'(A   )') 'the number of directions of applied '//
     &                    'magnetic field is not consistent:'
        Write(6,'(A,i4)') 'nDir    = ', nDir
        Write(6,'(A,i4)') 'nDirZee = ', nDirZee
        Write(6,'(A,i4)') 'nP      = ', nP
        Write(6,'(A,i4)') 'nDirTot = ', nDirTot
        Write(6,'(A,i4)') 'The rule is :'
        Write(6,'(A   )') 'nDir + nDirZee + nP = nDirTot'
        Call xFlush(6)
        Call abend()
      End If
c intialization
      Call dcopy_(nDirTot, [0.0_wp], 0, dHX(1), 1)
      Call dcopy_(nDirTot, [0.0_wp], 0, dHY(1), 1)
      Call dcopy_(nDirTot, [0.0_wp], 0, dHZ(1), 1)
      Call dcopy_(nDirTot, [0.0_wp], 0, dHW(1), 1)
      Call dcopy_(     nP, [0.0_wp], 0,   X(1), 1)
      Call dcopy_(     nP, [0.0_wp], 0,   Y(1), 1)
      Call dcopy_(     nP, [0.0_wp], 0,   Z(1), 1)
      Call dcopy_(     nP, [0.0_wp], 0,   W(1), 1)

c      If ( nDir.gt.0) then
c        Call DCOPY_(   nDir, dirX, 1, dHX(1), 1)
c        Call DCOPY_(   nDir, dirY, 1, dHY(1), 1)
c        Call DCOPY_(   nDir, dirZ, 1, dHZ(1), 1)
c      End If
c
c      If ( nDirZee.gt.0) then
c        Call DCOPY_(nDirZee, dir_weight(1:3,1), 1, dHX(1+nDir), 1)
c        Call DCOPY_(nDirZee, dir_weight(1:3,2), 1, dHY(1+nDir), 1)
c        Call DCOPY_(nDirZee, dir_weight(1:3,3), 1, dHZ(1+nDir), 1)
c      End If

      Do i=1,nDir
        dHX(i) = dirX(i)
        dHY(i) = dirY(i)
        dHZ(i) = dirZ(i)
      End Do

      Do i=1,nDirZee
        dHX(i+nDir) = dir_weight(i,1)
        dHY(i+nDir) = dir_weight(i,2)
        dHZ(i+nDir) = dir_weight(i,3)
      End Do


      Call Lebedev_Laikov(nSymm,nGrid,nP,X,Y,Z,W)

      Do i=1, nP
              j=i+nDir+nDirZee
        dHX(j) = X(i)
        dHY(j) = Y(i)
        dHZ(j) = Z(i)
        dHW(j) = W(i)
      End Do

c      Call DCOPY_(     nP, X(1), 1, dHX(1+nDir+nDirZee), 1)
c      Call DCOPY_(     nP, Y(1), 1, dHY(1+nDir+nDirZee), 1)
c      Call DCOPY_(     nP, Z(1), 1, dHZ(1+nDir+nDirZee), 1)
c      Call DCOPY_(     nP, W(1), 1, dHW(1+nDir+nDirZee), 1)
      Call qExit('hdir')
      Return
      End
