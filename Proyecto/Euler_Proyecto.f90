Program Euler_Proyecto
implicit none

integer :: i, j, kk, ptiempo
integer, parameter :: m = 128, n = 64
double precision :: dt, delta, Temp(n,m), Tempk(n,m), Calores(n,m), tiempo, alfa, k, Cp, ro
k = 100.0d0
ro = 8862.0d0
 Cp = 421.0d0

alfa = k/(ro*Cp)
delta = 2.0d0/m

write(*,*) 'Ingrese tamaño del paso de tiempo, número de pasos'
read(*,*) dt, tiempo

ptiempo = tiempo/dt
Tempk = 0.0d0
Temp = Tempk

Do kk = 1, ptiempo
	Do i = 1, n
		Do j = 2, m-1
			If (i == 1 .and. j <= m/2) Then !Mitad izquierda de la primera fila
				Temp(i,j) = Tempk(i,j) + alfa*(2*Tempk(i+1,j) + Tempk(i,j+1) + Tempk(i,j-1) + 200000*delta/k - 4*Tempk(i,j))*dt/(delta**2)

			Elseif (i == 1 .and. j > m/2) Then !Mitad derecha de la primera fila
				Temp(i,j) = Tempk(i,j) + alfa*(2*Tempk(i+1,j) + Tempk(i,j+1) + Tempk(i,j-1) - 4*Tempk(i,j))*dt/(delta**2)

			Elseif (i == n .and. j <= m/2) Then !Mitad izquierda de la última fila
				Temp(i,j) = Tempk(i,j) + alfa*(200000*delta/k + Tempk(i,j+1) + Tempk(i,j-1) + 2*Tempk(i-1,j) - 4*Tempk(i,j))*dt/(delta**2)

			Elseif (i == n .and. j > m/2) Then !Mitad derecha de la última fila
				Temp(i,j) = Tempk(i,j) + alfa*(Tempk(i,j+1) + Tempk(i,j-1) + 2*Tempk(i-1,j) - 4*Tempk(i,j))*dt/(delta**2)

			Else !Nodos internos
				Temp(i,j) = Tempk(i,j) + alfa*(Tempk(i,j+1) + Tempk(i+1,j) + Tempk(i-1,j) + Tempk(i,j-1) - 4*Tempk(i,j))*dt/(delta**2)
			End If
		End Do
		Temp(:,1) = 0.0d0
		Temp(:,m) = 0.0d0
	End Do
	Tempk = Temp
End Do

!write(*,*) Temp(1,1)

Open(Unit = 10, file = 'Blah.txt')
Do i = 1, n
	write(10,*) Temp(i,:)
End Do
 Close(Unit = 10)

!write(*,*) Temp

End Program
