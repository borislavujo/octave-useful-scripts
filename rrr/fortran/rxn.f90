!
!
 SUBROUTINE RXN(n,Kl0,vpl0,vba,lkab,lkba)
! **********************************************************************
! *                                                                    *
! * Calculates log rate constants between 2 states                     *
! *   using enhanced exponential lumping                               *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       03/02/2016                                             *
! * Version:    1.1                                                    *
! *                                                                    *
! **********************************************************************
!
   USE ParamsRXN, ONLY: howFine, lstartdt, thresh
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: vpl0
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n,n) :: Kl0
   LOGICAL, INTENT(IN), DIMENSION(n) :: vba
   DOUBLE PRECISION, INTENT(OUT) :: lkab, lkba
!
! -------------------------------------------------------------------
!
   INTEGER :: i, j, k
   DOUBLE PRECISION :: pl, ltemp, xl, doldd, dmax, pl1, pl2
   DOUBLE PRECISION :: lds, ldt, ltnow, ltau
   INTEGER :: nFine, na
   DOUBLE PRECISION, DIMENSION(n,n) :: Kl, L1, D1, Ls, Ds, Lp, Dp, &
        Dold, Ltem, Dtem
   LOGICAL, DIMENSION(n,n) :: Nb1, Nbs, Nbp, Nbtem
   DOUBLE PRECISION, DIMENSION(n) :: vlnow, vpla, vpl, vtemp
   DOUBLE PRECISION, DIMENSION(2) :: vUjo
   DOUBLE PRECISION, DIMENSION(3*n) :: vtr
   LOGICAL :: btemp, bl
   LOGICAL, DIMENSION(3*n) :: vbtr
   DOUBLE PRECISION, DIMENSION(9999) :: vXdata1, vXdata2 ! max ratio of timescales = 2^(9999/(2^3))
   INTEGER :: sizeX, nind
!   PARAMETER (howFine=3)
!   PARAMETER (lstartdt = -6d0)
!   PARAMETER (thresh = -6d0)
!
! -------------------------------------------------------------------
!
   WRITE(*,*) "Subroutine RXN"
!
!  make sure populations sum to 1 and rate mat satisfies detailed bal
!
   CALL LogSumExp(n,vpl0,pl)
   cycNormvpl: DO i=1,n
      vpl(i) = vpl0(i) - pl
   ENDDO cycNormvpl
   Kl = Kl0
   CALL SymmetriseRateMat(n,Kl,vpl)
!   WRITE(*,*) "Symmetric Kl"
!   DO i=1,3
!      WRITE(*,*) Kl(i,1:3)
!   ENDDO
!
!  calculate quilibrium populations
!
   WRITE(*,*) "Calc equil populs"
   na = 0
   cycCountA: DO i=1,n
      IF (vba(i)) THEN 
         na = na + 1 
         vpla(na) = vpl(i)
      ENDIF
   ENDDO cycCountA
   CALL LogSumExp(na,vpla,pl1)
   CALL LogDiffExp(0.0d0,pl1,pl2)
!
!  calculate optimal minimum (log) time step lds
!
   nFine = 2**howFine
   cycDiagRates: DO j=1,n
      cycRows: DO i=1,n
         vlnow(i) = Kl(i,j)
      ENDDO cycRows
      CALL LogSumExp(n,vlnow,pl)
      vtemp(j) = pl
   ENDDO cycDiagRates
!   WRITE(*,*) "vmaxes", vtemp
   lds = -MAXVAL(vtemp) + lstartdt - REAL(howFine)*LOG(2.0d0) ! - LOG(2.0d0)
   WRITE(*,'(A8,F12.7)') "lds", lds
!
!  produce log transition matrix
!
   WRITE(*,*) "Make the log trans mat"
   cycTMCol: DO j=1,n
      cycTMRow: DO i=1,n
         Ls(i,j) = Kl(i,j) + lds
         vlnow(i) = Ls(i,j)
      ENDDO cycTMRow
      vlnow(j) = vlnow(n)
      CALL LogSumExp(n-1,vlnow,ltemp)
      CALL LogDiffExp(0.0d0,ltemp,pl)
!      vUjo(1) = pl
!      vUjo(2) = Ls(j,j)
!      CALL LogSumExp(2,vUjo,ltemp)
!      WRITE(*,*) "ltemp", ltemp, "pl", pl
      Ls(j,j) = pl
   ENDDO cycTMCol
!
!  produce D and N matrices for log(1-t) formalism
!
   WRITE(*,*) "produce D and N matrices"
   cycDoRows: DO i=1,n
      cycDoCols: DO j=1,n
         IF (Ls(i,j).GT.vpl(i)) THEN
            Nbs(i,j) = .TRUE.
            CALL LogDiffExp(0.0d0,vpl(i)-Ls(i,j),ltemp)
            Ds(i,j) = ltemp + Ls(i,j) - vpl(i)
         ELSE
            Nbs(i,j) = .FALSE.
            CALL LogDiffExp(0.0d0,Ls(i,j)-vpl(i),ltemp)
            Ds(i,j) = ltemp
         ENDIF
      ENDDO cycDoCols
   ENDDO cycDoRows
!   WRITE(*,*) "Ls, Ds, Nbs"
!   DO j=1,3
!      WRITE(*,'(6F12.7,3L3)') Ls(j,:), Ds(j,:), Nbs(j,:)
!   ENDDO
!
!  log transition matrices for dt = nFine*exp(lds)
!
   WRITE(*,*) "get the matrices for nFine*lds"
   L1 = Ls
   D1 = Ds
   Nb1 = Nbs
   cycDoubling: DO i=1,howFine
!      WRITE(*,*) "which fine", i, "out of", howFine
!      CALL LPopul(n,L1,D1,Nb1,vpl,vba,xl)
      CALL MultLogMat(n,vpl,L1,D1,Nb1,L1,D1,Nb1,Ltem,Dtem,Nbtem)
      L1 = Ltem
      D1 = Dtem
      Nb1 = Nbtem
!      WRITE(*,*) "L1, D1, Nb1"
!      DO j=1,3
!         WRITE(*,'(6F12.7,3L3)') L1(j,:), D1(j,:), Nb1(j,:)
!      ENDDO
   ENDDO cycDoubling
   xl = 0
   vXdata1(1) = -99e9
   vXdata2(1) = xl
   ldt = lds + REAL(howFine)*LOG(2.0d0)
   doldd = 1.0d1
   dmax = 1.0d1
   nind = 1
!
!  main cycle
!
   WRITE(*,*) "Main cycle"
   cycMain: DO WHILE( ((xl.GT.-5).OR.(doldd.GT.-3)) .AND.(dmax.GT.thresh))
      Dold = D1
!      WRITE(*,*) "L(1,1)", L1(1,1)
!      WRITE(*,*) "L1, D1, Nb1", nind
!      DO j=1,3
!         WRITE(*,'(6F12.7,3L3)') L1(j,:), D1(j,:), Nb1(j,:)
!      ENDDO
      CALL LPopul(n,L1,D1,Nb1,vpl,vba,xl)
      vXdata1(nind) = ldt
      vXdata2(nind) = xl
!      WRITE(*,*) nind, ldt, xl
      nind = nind + 1
      Lp = L1
      Dp = D1
      Nbp = Nb1
      ltnow = ldt
      cycLinProp: DO i=1,nFine-1
         CALL MultLogMat(n,vpl,Ls,Ds,Nbs,Lp,Dp,Nbp,Ltem,Dtem,Nbtem)
         Lp = Ltem
         Dp = Dtem
         Nbp = Nbtem
!         WRITE(*,*) "Lp, Dp, Nbp", nind
!         DO j=1,3
!            WRITE(*,'(6F12.7,3L3)') Lp(j,:), Dp(j,:), Nbp(j,:)
!         ENDDO
         CALL LPopul(n,Lp,Dp,Nbp,vpl,vba,xl)
         vUjo(1) = ltnow
         vUjo(2) = lds
         CALL LogSumExp(2,vUjo,ltemp)
         ltnow = ltemp
         vXdata1(nind) = ltnow
         vXdata2(nind) = xl
         WRITE(*,'(3A10,I6,2F12.7)') "nind", "ltnow", "xl", nind, ltnow, xl
         nind = nind + 1
      ENDDO cycLinProp
      CALL MultLogMat(n,vpl,L1,D1,Nb1,L1,D1,Nb1,Ltem,Dtem,Nbtem)
      L1 = Ltem
      D1 = Dtem
      Nb1 = Nbtem
      CALL MultLogMat(n,vpl,Ls,Ds,Nbs,Ls,Ds,Nbs,Ltem,Dtem,Nbtem)
      Ls = Ltem
      Ds = Dtem
      Nbs = Nbtem
      ldt = ldt + LOG(2.0d0)
      lds = lds + LOG(2.0d0)
      doldd = LOG(SUM(ABS(D1-Dold)))
      dmax = MAXVAL(D1)
   ENDDO cycMain
   CALL TrapzLog(nind-1,vXdata1,vXdata2,ltau)
!   WRITE(*,*) "vypoctov", nind
!   WRITE(*,*) "ltau", ltau
   lkab = pl1 - ltau
   WRITE(*,*) "lkab", lkab
   lkba = pl2 - ltau
!
   RETURN
!
! -------------------------------------------------------------------
!
 END SUBROUTINE RXN
!
!
!
 SUBROUTINE SymmetriseRateMat(n,Kl,vpl)
! **********************************************************************
! *                                                                    *
! * Ensures the rate matrix satisfies detailed balance                 *
! *   Kl must have -99e9 for zero rates and -1e99 at its diagonal      *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       02/02/2016                                             *
! * Version:    1.0                                                    *
! *                                                                    *
! **********************************************************************
!
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: vpl
   DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n,n) :: Kl
!
! -------------------------------------------------------------------
!
   INTEGER :: i, j
   DOUBLE PRECISION :: recTau, pl
   DOUBLE PRECISION, DIMENSION(2) :: vUjo
   LOGICAL, DIMENSION(3*n) :: vbtr
!
! -------------------------------------------------------------------
!
!   cycDoDiags: DO i=1,n
!      Kl(i,i) = -1e99
!   ENDDO cycDoDiags
   cycDoRows: DO i=1,n-1
      cycDoCols: DO j=i+1,n
         vUjo(1) = Kl(i,j)
         vUjo(2) = Kl(j,i)
         CALL LogSumExp(2,vUjo,recTau)
         vUjo(1) = vpl(i)
         vUjo(2) = vpl(j)
         CALL LogSumExp(2,vUjo,pl)
         Kl(i,j) = vpl(i) + recTau - pl
         Kl(j,i) = vpl(j) + recTau - pl
      ENDDO cycDoCols
   ENDDO cycDoRows
!
   RETURN
!
! -------------------------------------------------------------------
!
 END SUBROUTINE SymmetriseRateMat
!
!
!
 SUBROUTINE TrapzLog(n,vX1,vX2,ltau)
! **********************************************************************
! *                                                                    *
! * Integrates a vector in log formalism usind trapezoidal ruel        *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       02/02/2016                                             *
! * Version:    1.0                                                    *
! *                                                                    *
! **********************************************************************
!
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: vX1, vX2
   DOUBLE PRECISION, INTENT(OUT) :: ltau
!
! -------------------------------------------------------------------
!
   INTEGER :: i
   DOUBLE PRECISION :: ldt, pl
   DOUBLE PRECISION, DIMENSION(2) :: vUjo
!
! -------------------------------------------------------------------
!
!   WRITE(*,*) "X"
   ltau = -9.9d99
   cycTrapez: DO i=1,n-1
!      WRITE(*,'(4F12.7)') vX1(i), vX2(i), vX1(i+1), vX2(i+1)
      CALL LogDiffExp(vX1(i+1),vX1(i),ldt)
      vUjo(1) = vX2(i)
      vUjo(2) = vX2(i+1)
      CALL LogSumExp(2,vUjo,pl)
      pl = pl - LOG(2.0d0)
      vUjo(1) = ltau
      vUjo(2) = ldt + pl
      CALL LogSumExp(2,vUjo,ltau)
!      WRITE(*,*) "ltau", ltau
  ENDDO cycTrapez
 !
   RETURN
!
! -------------------------------------------------------------------
!
 END SUBROUTINE TrapzLog
!
