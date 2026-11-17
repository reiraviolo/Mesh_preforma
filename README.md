## IDEA PRINCIPAL:
* Generar mallas pre-establecidas para componentes de cañerias ( tee, codos, reducciones y caños) y luego combinarlos para formar una malla de elementos 2D (placa) de una geometría compuesta a evaluar en un programa FEM ( no alcanzado en este proyecto)

## Partes de tubería a pre-mallar

**Sistema de Piping
  * Caño  - variables: OD, Longitud
  * Codo -  variables: OD, H, Radio
  * Tee  - variables: OD_1 (Directo), OD_2 (derivación), L_1 (ingreso), L_2 (derivación), L_3 (salida)
  * Reducción - variables: OD_1 (entrada), OD_2 ( salida), Longitud. (perfil pre-formado)


##  Variables comunes para poder acoplar las mallas
* OD = Díametro de tubería 
* N = Número de divisiones en la circunferencia

## Primera etapa: Mallar cada componente de manera individual
* Mallado de piezas individuales.

## Segunda etapa: Armar un tramo de cañeria combinando varios  de estos accesorios
* Generación de aplicación para integración de componentes.
