function Pos_sensores=Calcula_Pos_sensores(Area,Pixel, sigma_sensor_x,sigma_sensor_y,Nceldas)
%Calcula la posicion de los sensores cuando siguen una variable aleatoria
%gaussiana
Pos_lin_sensores_x=(1:fix(Area(1)/Pixel(1)))*Pixel(1)-Pixel(1)/2;
Pos_lin_sensores_y=(1:fix(Area(2)/Pixel(2)))*Pixel(2)-Pixel(2)/2;
[PosicionSensoresx1,PosicionSensoresy1]=meshgrid(Pos_lin_sensores_x,Pos_lin_sensores_y);

Pos_sensores=[PosicionSensoresx1(:)'+sigma_sensor_x*randn(1,Nceldas);PosicionSensoresy1(:)'+sigma_sensor_y*randn(1,Nceldas)];