
!
!
 SUBROUTINE MultLogMat(n,vpl,L1,D1,Nb1,L2,D2,Nb2,L3,D3,Nb3)
! **********************************************************************
! *                                                                    *
! * Multiplies transition matrices in the log and log(1-t) formalisms  *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       11/02/2016                                             *
! * Version:    1.2                                                    *
! *                                                                    *
! **********************************************************************
!
   USE ParamsRXN, ONLY: ln2
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: vpl
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n,n) :: L1, L2
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n,n) :: D1, D2
   LOGICAL, INTENT(IN), DIMENSION(n,n) :: Nb1, Nb2
   DOUBLE PRECISION, INTENT(OUT), DIMENSION(n,n) :: L3
   DOUBLE PRECISION, INTENT(OUT), DIMENSION(n,n) :: D3
   LOGICAL, INTENT(OUT), DIMENSION(n,n) :: Nb3
!
! -------------------------------------------------------------------
!
   INTEGER :: i, j, k
   DOUBLE PRECISION :: pl, ltemp, pl2, pplus, pminus
   DOUBLE PRECISION :: h
   INTEGER :: nind, nplus, nminus
   DOUBLE PRECISION, DIMENSION(n) :: vlnow, vplus, vminus
   DOUBLE PRECISION, DIMENSION(2) :: vUjo
   DOUBLE PRECISION, DIMENSION(3*n) :: vtr
   LOGICAL :: btemp, bl
   LOGICAL, DIMENSION(3*n) :: vbtr
   LOGICAL, DIMENSION(2) :: vbUjo
   LOGICAL, DIMENSION(n,n) :: Nbted
   PARAMETER (h=-3.0d0)
!
! -------------------------------------------------------------------
!
!
!  L3 = L1 * L2
!
   cycCalcRowL3: DO j=1,n
      cycCalcColL3: DO i=1,n
         cycCorrL12: DO k=1,n
            vlnow(k) = L1(i,k)+L2(k,j)
         ENDDO cycCorrL12
         CALL LogSumExp(n,vlnow,ltemp)
         L3(i,j) = ltemp
      ENDDO cycCalcColL3
   ENDDO cycCalcRowL3
!
!  symmetrise L3
!
   cycDoRows: DO i=1,n-1
      cycDoCols: DO j=i+1,n
         vUjo(1) = L3(i,j)
         vUjo(2) = L3(j,i)
         CALL LogSumExp(2,vUjo,ltemp)
         vUjo(1) = vpl(i)
         vUjo(2) = vpl(j)
         CALL LogSumExp(2,vUjo,pl)
         L3(i,j) = vpl(i) + ltemp - pl
         L3(j,i) = vpl(j) + ltemp - pl
      ENDDO cycDoCols
   ENDDO cycDoRows
!
!  normalise L3
!
   cycL3norm: DO j=1,n
      cycGetL3Trms: DO i=1,n
         vlnow(i) = L3(i,j)
      ENDDO cycGetL3Trms
      CALL LogSumExp(n,vlnow,ltemp)
      cycNormL3: DO i=1,n
         L3(i,j) = L3(i,j) - ltemp
      ENDDO cycNormL3
   ENDDO cycL3norm
!
!  calculate D3 AND Nb3
!
   cycCalcRowD3: DO i=1,n
      cycCalcColD3: DO j=1,n
         IF (L3(i,j)-vpl(i).GT.h) THEN
            cycFillvtr: DO k=1,n
               vtr((k-1)*3+1)  = vpl(k) + D1(i,k)
               vtr((k-1)*3+2)  = vpl(k)           + D2(k,j)
               vtr(k*3)        = vpl(k) + D1(i,k) + D2(k,j)
               vbtr((k-1)*3+1) = Nb1(i,k)
               vbtr((k-1)*3+2) = Nb2(k,j)
               vbtr(k*3)       = Nb1(i,k).EQV.Nb2(k,j)
            ENDDO cycFillvtr
            CALL LogSumDiff(3*n,vtr,vbtr,ltemp,btemp)
            Nb3(i,j) = btemp
            D3(i,j) = ltemp
         ELSE
            IF (L3(i,j).GT.vpl(i)) THEN
               Nb3(i,j) = .TRUE.
!               D3(i,j) = LOG(EXP(L3(i,j)-vpl(i))-1)
               CALL LogDiffExp(L3(i,j)-vpl(i),0.0d0,pl)
               D3(i,j) = pl
            ELSE
               Nb3(i,j) = .FALSE.
!               D3(i,j) = LOG(1-EXP(L3(i,j)-vpl(i)))
               CALL LogDiffExp(0.0d0,L3(i,j)-vpl(i),pl)
               D3(i,j) = pl
            ENDIF
         ENDIF
      ENDDO cycCalcColD3
   ENDDO cycCalcRowD3
!
!  symmetrise D3
!
   cycSymmD3Row: DO i=1,n
      cycSymmD3Col: DO j=1,n
         IF (D3(i,j).GT.D3(j,i)) THEN
            Nb3(j,i) = Nb3(i,j)
         ELSE
            Nb3(i,j) = Nb3(j,i)
         ENDIF
         vUjo(1) = D3(i,j)
         vUjo(2) = D3(j,i)
         CALL LogSumExp(2,vUjo,pl)
         D3(i,j) = pl - ln2
         D3(j,i) = D3(i,j)
      ENDDO cycSymmD3Col
   ENDDO cycSymmD3Row
!
!  normalise d3
!
   cycNormCols: DO j=1,n
      nplus = 0
      nminus = 0
      cycGetTrms: DO i=1,n
         vtr(i)  = D3(i,j)+vpl(i)
!         vtr(i)  = D3(i,j) ... this was the mistake
         vbtr(i) = Nb3(i,j)
         IF (Nb3(i,j)) THEN
            nplus = nplus+1
            vplus(nplus) = vtr(i)
         ELSE
            nminus = nminus+1
            vminus(nminus) = vtr(i)
         ENDIF
      ENDDO cycGetTrms
      CALL LogSumDiff(n,vtr,vbtr,pl,bl)
      IF (nplus.GT.0) THEN
         CALL LogSumExp(nplus,vplus,pplus)
      ELSE
         pplus = -99e9
         WRITE(*,*) "warning: no positive terms"
      ENDIF
      IF (nminus.GT.0) THEN
         CALL LogSumExp(nminus,vminus,pminus)
      ELSE
         pplus = -99e9
         WRITE(*,*) "warning: no positive terms"
      ENDIF
      cycNormD3: DO i=1,n
         
!      WRITE(*,*) "D3 col before norm", j, "sum row", pl, "smerom", bl
!      cycDoTrms: DO i=1,n
!         vUjo(1) = MIN(pl-1.0e-2,D3(i,j)-1.0e-2)
!         vUjo(1) = MIN(pl+vpl(i),D3(i,j)-1e-5)
!         vUjo(2) = D3(i,j)
!         vbUjo(1) = .NOT.bl
!         vbUjo(2) = Nb3(i,j)
!         CALL LogSumDiff(2,vujo,vbujo,ltemp,btemp)
!         D3(i,j) = ltemp
!         vlnow(i) = ltemp
         IF (Nb3(i,j).EQV.bl) THEN
            D3(i,j) = D3(i,j) - ABS(pplus-pminus)
         ENDIF
      ENDDO cycNormD3
      cycCheckTrms: DO i=1,n
         vtr(i)  = vlnow(i)+vpl(i)
         vbtr(i) = Nb3(i,j)
      ENDDO cycCheckTrms
      CALL LogSumDiff(n,vtr,vbtr,pl2,bl)
!
!  i have no idea how this normalisation can make things worse
!    correct D3 only if normalisation improves it
!
      IF (ABS(pl).GT.ABS(pl2)) THEN
         DO i=1,n
            D3(i,j) = vlnow(i)
         ENDDO
      ENDIF
   ENDDO cycNormCols
!
   RETURN
!
! -------------------------------------------------------------------
!
 END SUBROUTINE MultLogMat
!
