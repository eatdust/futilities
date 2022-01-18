INTERFACE
  SUBROUTINE fcn1(m, n, x, fvec, iflag)
    use precision, only : fdp
    IMPLICIT NONE
    INTEGER, PARAMETER         :: lmdp = fdp
    INTEGER, INTENT(IN)        :: m, n
    REAL (lmdp), INTENT(IN)      :: x(:)
    REAL (lmdp), INTENT(IN OUT)  :: fvec(:)
    INTEGER, INTENT(IN OUT)    :: iflag
  END SUBROUTINE fcn1



  SUBROUTINE fcn2(m, n, x, fvec, fjac, iflag)
    use precision, only : fdp
    IMPLICIT NONE
    INTEGER, PARAMETER         :: dp = fdp
    INTEGER, INTENT(IN)        :: m, n
    REAL (dp), INTENT(IN)      :: x(:)
    REAL (dp), INTENT(IN OUT)  :: fvec(:)
    REAL (dp), INTENT(OUT)     :: fjac(:,:)
    INTEGER, INTENT(IN OUT)    :: iflag
  END SUBROUTINE fcn2

END INTERFACE
