PROGRAM GaussSeidel

	implicit none

	integer, parameter::n=10
	integer i, j, k, stat, imax
	double precision a(n,n), e(n,n), b(n), x(n), x1(n), c(n), r(n), es, error(n), f(n), lambda !real kind 8

	!tol=1.0d-10

!Definicion Matriz A

	OPEN(UNIT=10,FILE='mat_diag.txt',IOSTAT=stat,ACTION='READ')
	IF(stat/=0) THEN
		WRITE(*,*) 'Opening unit 10 failed with iostat',stat, '.'
	END IF

	DO i=1,n
		READ (10,*) (a(i,j), j=1,n)
			do j=1 , n
			e(i,j)=a(i,j)
			end do
	END DO

	READ(10,*)
	DO i=1,n
		READ(10,*) b(i)
        f(i)=b(i)
	END DO

	CLOSE(UNIT=10, IOSTAT=stat)
	IF (stat /=0) THEN
		WRITE(*,*) 'closing unit 10 failed with',stat,'.'
	END IF

!Presentacion matriz
	WRITE(*,*) 'La matriz A='
	do i=1,n
		write(*,01) (a(i,j),j=1,n)
	end do

	WRITE(*,*) 'La matriz B='
	do i=1,n
		write(*,01) b(i)
    end do

	

	call gaussito(a,b,n,x,imax,es,lambda)
    
    WRITE(*,*) 'x= '
    DO i=1,n
            WRITE(*,01) x(i)
    END DO
    
    01 format(1x,1p3e15.7)
   
END PROGRAM GaussSeidel

!-------------------------------------------------------------------------------------------------------------------------------------------

subroutine gaussito(a,b,n,x,imax,es,lambda)

    integer n,i,j,k,m,iter, imax
	double precision a(n,n), b(n), x(n), es, lambda, dummy, suma, centinela, old,s

    es=0.00001
    lambda=0.05
    imax=1000000
    
    !do i=1,n 
     !   s=0
      !  do j=1,n
       !     if(i/=j)then 
        !        s=s+abs(a(i,j))
         !   end if
    !    end do
     !   if(abs(a(i,i)).lt.s) then
      !      write(*,*) 'these equations are not diagonal'
       !     stop
        !end if
    !end do
    
    
    DO i=1, n
        dummy=a(i,i)
        DO j=1, n
            a(i,j)=a(i,j)/dummy
        END DO
        b(i)=b(i)/dummy
    END DO
    
    DO i=1,n
        sum=b(i)
        DO j=1,n
            IF (i/=j) Then
                sum=sum-a(i,j)*x(j)
            END IF
        END DO
        x(i)=sum
    END DO
    
    iter=1
    centinela=0
    DO k=1,imax
        centinela=1
        DO i=1, n
            old=x(i)
            sum=b(i)
            DO j=1, n
                IF (i/=j) Then
                    sum=sum-a(i,j)*x(j)
                END IF
            END DO
            x(i)=lambda*sum+(1.0d0-lambda)*old
            IF ((centinela == 1) .AND. x(i)/=0) then
                ea=abs((x(i)-old)/x(i))*100
                IF (ea > es) Then
                    centinela = 0
                END IF
            END IF
        END DO
        iter=iter+1
        IF (centinela==1) EXIT
    END DO
    write(*,*) iter
                
end subroutine
