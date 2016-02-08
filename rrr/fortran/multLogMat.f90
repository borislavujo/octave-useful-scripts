
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
! * Date:       03/02/2016                                             *
! * Version:    1.1                                                    *
! *                                                                    *
! **********************************************************************
!
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
   DOUBLE PRECISION :: pl, ltemp
   DOUBLE PRECISION :: h
   INTEGER :: nind
   DOUBLE PRECISION, DIMENSION(n) :: vlnow
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
!
!  normalise L3
!
      cycCalcSumColL3: DO i=1,n
         vlnow = L3(i,j)
      ENDDO cycCalcSumColL3
      CALL LogSumExp(n,vlnow,ltemp)
      cycNormColL3: DO i=1,n
         L3(i,j) = L3(i,j) - ltemp
      ENDDO cycNormColL3
   ENDDO cycCalcRowL3
!
!  calculate D3 AND Nb3
!
   cycCalcRowD3: DO i=1,n
      cycCalcColD3: DO j=1,n
         IF (L3(i,j)-vpl(i).GT.h) THEN
            Nbted(i,j) = .TRUE.
            cycFillvtr: DO k=1,n
               vtr((k-1)*3+1)  = vpl(k) + D1(i,k)
               vtr((k-1)*3+2)  = vpl(k)           + D2(k,j)
               vtr(k*3)        = vpl(k) + D1(i,k) + D2(k,j)
               vbtr((k-1)*3+1) = Nb1(i,k)
               vbtr((k-1)*3+2) = Nb2(k,j)
               vbtr(k*3)       = Nb1(i,k).EQV.Nb2(k,j)
!(Nb1(i,k).AND.Nb2(k,j)).OR. ((.NOT.Nb1(i,k)).AND.(.NOT.Nb2(k,j)))
            ENDDO cycFillvtr
            CALL LogSumDiff(3*n,vtr,vbtr,ltemp,btemp)
            Nb3(i,j) = btemp
            D3(i,j) = ltemp
         ELSE
            Nbted = .FALSE.
            IF (L3(i,j).GT.vpl(i)) THEN
               Nb3(i,j) = .TRUE.
               D3(i,j) = LOG(EXP(L3(i,j)-vpl(i))-1)
            ELSE
               Nb3(i,j) = .FALSE.
               D3(i,j) = LOG(1-EXP(L3(i,j)-vpl(i)))
            ENDIF
         ENDIF
      ENDDO cycCalcColD3
   ENDDO cycCalcRowD3
   DO i=1,3
      WRITE(*,'(A10,3F12.7)') "beforenorm", D3(i,:)
   ENDDO
!
!  normalisation
!
   cycNormCols: DO j=1,n
      cycGetTrms: DO i=1,n
         vtr(i)  = D3(i,j)+vpl(i)
!         vtr(i)  = D3(i,j) ... this was the mistake
         vbtr(i) = Nb3(i,j)
      ENDDO cycGetTrms
      CALL LogSumDiff(n,vtr,vbtr,pl,bl)
!      WRITE(*,*) "D3 col before norm", j, "sum row", pl, "smerom", bl
!      WRITE(*,*) "vtr", vtr(1:n), "vbtr", vbtr(1:n)
      cycDoTrms: DO i=1,n
         vUjo(1) = MIN(pl,D3(i,j)-LOG(2.0d0))
!         vUjo(1) = MIN(pl+vpl(i),D3(i,j)-1e-5)
         vUjo(2) = D3(i,j)
         vbUjo(1) = .NOT.bl
         vbUjo(2) = Nb3(i,j)
         CALL LogSumDiff(2,vujo,vbujo,ltemp,btemp)
         D3(i,j) = ltemp
      ENDDO cycDoTrms
      cycCheckTrms: DO i=1,n
         vtr(i)  = D3(i,j)+vpl(i)
         vbtr(i) = Nb3(i,j)
      ENDDO cycCheckTrms
      CALL LogSumDiff(n,vtr,vbtr,pl,bl)
!      WRITE(*,*) "D3 column", j, "sum row", pl, "smerom", bl
!      WRITE(*,*) "vtr", vtr(1:n), "vbtr", vbtr(1:n)
      cycGetL3Trms: DO i=1,n
         vlnow(i) = L3(i,j)
      ENDDO cycGetL3Trms
      CALL LogSumExp(n,vlnow,ltemp)
      cycNormL3: DO i=1,n
         L3(i,j) = L3(i,j) - ltemp
      ENDDO cycNormL3
   ENDDO cycNormCols
!
   RETURN
!
! -------------------------------------------------------------------
!
 END SUBROUTINE MultLogMat
!
