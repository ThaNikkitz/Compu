Program Euler_Proyecto
implicit none

integer :: i, j, kk, ptiempo
integer, parameter :: m = 128, n = 64
double precision :: dt, delta, Temp(n,m), Tempk(n,m), Calores(n,m), tiempo, alfa, k, Cp, ro
k = 100.0d0
ro = 8862.0d0
Cp = 421.0d0

alfa = k/(ro*Cp)

write(*,*) 'Ingrese tamaño del paso de tiempo, número de pasos'
read(*,*) dt, tiempo

ptiempo = tiempo/dt
Tempk = 0.0d0
Temp = Tempk

Do kk = 1, ptiempo
	Do i = 1, n
		Do j = 1, m
			Tempk(1,:) = 0.0d0
			Tempk(n,:) = 0.0d0
			If (j == 1 .and. i < n/2) Then !Mitad izquierda de la primera fila
				Temp(i,j) = Tempk(i,j) + alfa*(2*Tempk(i,j+1) + Tempk(i+1,j) + Tempk(i-1,j) + 2000*delta/k - 4*Tempk(i,j))*dt/(delta**2)
			Elseif (j == 1 .and. i > n/2) Then !Mitad derecha de la primera fila
				Temp(i,j) = Tempk(i,j) + alfa*(2*Tempk(i,j+1) + Tempk(i+1,j) + Tempk(i,j-1) - 4*Tempk(i,j))*dt/(delta**2)
			Elseif (j == m .and. i < n/2) Then !Mitad izquierda de la última fila
				Temp(i,j) = Tempk(i,j) + alfa*(2000*delta/k + Tempk(i+1,j) + Tempk(i-1,j) + 2*Tempk(i,j-1) - 4*Tempk(i,j))*dt/(delta**2)
			Elseif (j == m .and. i > n/2) Then !Mitad derecha de la última fila
				Temp(i,j) = Tempk(i,j) + alfa*(Tempk(i+1,j) + Tempk(i-1,j) + 2*Tempk(i,j-1) - 4*Tempk(i,j))*dt/(delta**2)
			Else !Nodos internos
				Temp(i,j) = Tempk(i,j) + alfa*(Tempk(i,j+1) + Tempk(i+1,j) + Tempk(i-1,j) + Tempk(i,j-1) - 4*Tempk(i,j))*dt/(delta**2)
			End If
		End Do
	End Do
	Tempk = Temp
End Do

write(*,*) Temp

End Program
