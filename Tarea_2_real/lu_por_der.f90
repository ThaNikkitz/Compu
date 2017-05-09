	program metodos_directos

	integer, parameter :: np = 10, n = 10
	integer indx(n), stat
!	double precision a(np, np), y(np, np), b(n), d, a1(np, np), num, c(np,np)
	real(kind=8) :: a(np, np), y(np, np), b(n), d, a1(np, np), num, c(np,np)
!	Definición Matriz A

    OPEN(UNIT=11,FILE='viento_der.txt',IOSTAT=stat,ACTION='READ')
    IF (stat /= 0) THEN
        WRITE(*,*) 'Opening Unit 10 Failed with iostat ', stat, '.'
    END IF
    
    DO i = 1,n
        READ(11,*) (a(i,j), j = 1,n)
    END DO
    READ(11,*) 
    DO i = 1,n
        READ(11,*) b(i)
    END DO
    
    CLOSE(UNIT=11,IOSTAT=stat)
    IF (stat /= 0) THEN
        WRITE(*,*) 'Closing Unit 10 Failed with iostat ', stat, '.'
    END IF

	a1 = a

01  format(1x,1p6e15.7)

!	Presentación Matriz

!    write(*,*) 'La matriz A='
    do i=1,n
!    write(*,01) (a(i,j), j = 1,n)
    end do

!    write(*,*) 'El vector b='
    do i=1,n
    write(*,01) b(i)
    end do
	

	call lu_pivoteado(a, n, np, indx, d)
	
	call luevaluacion(a, n, np, indx, b)

	call lu_inversa(a, n, np, indx, y)

	call numero_condicion(a1, y, n, np, num)

	write(*,*) 'Las fuerzas que soporta cada nodo son:'
	do i = 1, n
	write(*,*) 'F',i , '=', b(i)
	end do

	write(*,*)
	write(*,*)'El número de Condición de A es'
	write(*,01) num

	end


!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------




	subroutine lu_pivoteado(a, n, np, indx, d)
	
	integer, parameter :: NMAX = 500
	integer i, j, k, imax, n, np, indx(n)
	double precision aamax, dum, suma, vv(NMAX), d, a(np,np)
	double precision, parameter :: min = 1.0d-20

!	Dada una matriz a(1:n,1:n), con dimensiones físicas np x np,
!	esta rutina reemplaza esta por la descomposición LU con 
!	permutaciones por fila en a. La matriz a y el entero n son 
!	entradas. Además a es salida conteniendo en sí a la descomposicion
!   LU omitiendo la la diagonal de 1 de la matriz L; indx(1:n) es un
!	vector de salida que registra las permutaciones por filas efectuadas
!	por el pivoteo parcial; d es una salida +1 o -1 dependiendo si el
!	número de filas intercambiado	es par o impar respectivamente.
!	El arreglo vv registra el escalamiento implicito en cada fila.

	d = 1.0d0

!	Loop sobre las filas para obtener la información del escalamiento 
!	implicito.

	do i = 1, n
	  aamax = 0.0d0
	  do j = 1, n
	    if ( abs(a(i, j)) > aamax ) aamax = abs( a(i, j) )
	  enddo
	  if ( aamax == 0.0d0 ) stop 'Matriz singular, caso 1' 
	  vv(i) = 1.0d0 / aamax
	enddo

!	Ahora se aplicará un loop para el método de Crout (o reducción de Gauss).

	do j = 1, n
	  do i = 1, j - 1
	    suma = a(i, j)
		do k = 1, i - 1
		  suma = suma - a(i, k) * a(k, j)
		enddo
		a(i, j) = suma
	  enddo
	  aamax = 0.0d0

!	Inicio de búsqueda de elemento pivote.

	  do i = j, n
	    suma = a(i, j)
		do k = 1, j - 1
		  suma = suma - a(i, k) * a(k, j)
		enddo
		a(i, j) = suma
		dum = vv(i) * abs(suma)
		if ( dum >= aamax ) then
		  imax = i
		  aamax = dum
		endif
	  enddo
	  if ( j /= imax ) then
	    do k = 1, n
		  dum = a(imax, k)
		  a(imax, k) = a(j, k)
		  a(j, k) = dum
		enddo
		d = - 1.0d0 * d
		vv(imax) =vv(j)
	  endif
	  indx(j) = imax
	  if ( a(j, j) == 0.0d0 ) a(j, j) = min
	  if ( j /= n ) then
	    dum = 1.0d0 / a(j, j)
		do i = j + 1, n
		  a(i, j) = a(i, j) * dum
		enddo
	  endif
	enddo
	
	end subroutine


!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------



	subroutine luevaluacion(a, n, np, indx, b)

	integer i, j, ii, ll, n, np, indx(n)
	double precision a(np, np), b(n), suma

!	Solver para un set de n ecuaciones lineales en la forma Ax = b, donde
!	A es la matriz correspondiente a la descomposición LU determinada en 
!	la subrutina lu_pivoteada. Aquí b(1:n) es una entrada correspondiente
!	al vector columna b del lado derecho del sistema y sale como el vector
!	columna x del sistema (esta entrada es la única que se modifica, por 
!	la que el resto puede ocuparse nuevamente con otro vector b de entrada).
!	El arreglo indx es ingresado como un vector que contiene las
!	permutaciones realizadas en la subrutina lu_pivoteada.

	ii = 0.0d0
	do i = 1, n
	  ll = indx(i)
	  suma = b(ll)
	  b(ll) = b(i)
	  if ( ii /= 0.0d0 ) then
	    do j = ii, i - 1
		  suma = suma - a(i, j) * b(j)
		enddo
	  elseif ( suma /= 0.0d0 ) then
	    ii = i
	  endif
	  b(i) = suma
	enddo

	do i = n, 1, -1
	  suma = b(i)
	  do j = i + 1, n
	    suma = suma - a(i, j) * b(j)
	  enddo
	  b(i) = suma / a(i, i)
	enddo

	end subroutine



!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------



	subroutine lu_inversa(a, n, np, indx, y)

	implicit none
	integer n, np, i, j, indx(np)
	double precision a(np, np), y(np, np), d

	do i = 1, n
	  do j = 1, n
	    y(i, j) = 0.0d0
	  enddo
	  y(i, i) = 1.0d0
	enddo
!	write(*,*) y	
	do j = 1, n
	  call luevaluacion(a, n, np, indx, y(1, j))
	enddo

	end subroutine



!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------



	subroutine numero_condicion(a, y, n, np, num)

!	Esta subrutina esta destinada al cálculo del número de condición de 
!	una matriz a(1:n,1:n) con dimensión física np x np. Aquí a e y son 
!	la matriz A y su inversa A^(-1) respectiva. num es la única salida 
!	que contiene el número de condición calculado mediante la norma 
!	infinito de matrices.

	integer n, np
	double precision a(np,np), y(np,np), num, suma, d, d1
	write(*,*)
	write(*,*) 'La matriz invertida es:'
	write(*,*) y
	write(*,*)

	OPEN(UNIT=8,FILE='mat_inv_p4b.txt', ACTION='write')
		
	do i = 1, n 
		write(8,*) (y(i,j), j = 1,n)
	end do

	Close(UNIT = 8)

	suma = 0.0d0
	d = 0.0d0
	d1 = 0.0d0

!	Cálculo tamaño de A utilizando norma infinito.

	do i = 1, n
	  suma =0.0d0
	  do j = 1, n
	    suma = suma + abs(a(i,j))
	  enddo
	  if ( abs(suma) >= d ) then
	    d = suma
	  endif
	enddo

!	Cálculo tamaño de A^(-1) utilizando norma infinito.

	do i = 1, n
	  suma =0.0d0
	  do j = 1, n
	    suma = suma + abs(y(i,j))
	  enddo
	  if ( abs(suma) >= d1 ) then
	    d1 = suma
	  endif
	enddo

!	Cálculo número de Condición

	num = d * d1

	end subroutine
