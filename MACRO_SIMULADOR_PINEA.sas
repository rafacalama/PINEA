/* 31/08/22. Macro prueba Marta Ezquerro*/
/* Código para para la Unidad Natural: Viana, asumimos índice de sitio=17. Puede cambiarse para otras unidades*/
/* Cosas que puedo cambiar:
- Bases de datos iniciales
- Datos climáticos
- Periodo de simulación
- Año inciial de simulaciój (en este caso fijo 2010)
- Esquema de claras y corta final
- Turno de corta
- Límites mortalidad por autoaclareo
- Diámetro mínimo para volumen de sierra*/

/* Puedo simular distintas parcelas de 1 ha a la vez*/

/* Ojo! Todas las bases de datos necesarias para la simulación, también disponibles, se incluirían en una librería (en mi caso llamada FTERRA) */


%macro viana_BIOMASA; 

/* Cargo una base de datos inicial, que debe incluir al menos un codigo de parcela, un listado de los árboles y un diámetro normal*/
/* Hay que asignarle una edad de rodal, un índice de sitio, y conocer la Unidad Natura*/
/* La aplicación simula rodales - parcelas de 1 hectárea*/

/* Mi base de datos inicial simula 5 parcelas de 1 ha, cada una de una clase de edad (20-40-60-80-100)*/
/* Las condiciones iniciales son simuladas para una SI=17. Pueden cambiarse por otras */

data datos_ini; set FTERRA.viana_INI;run; 
data datos_ini;set datos_ini;

aux=1;year=2010; parcela=cledad;t=t-1;si=17; run;

data datos_ini; set datos_ini; drop cledad;run;



/* lo primero le asignamos el valor qq del ef aleatorio arbol a cada observacion*/

data datos_ini; set datos_ini;
qq = rannor(1)*sqrt(0.0088);
run;

/* Importo de los datos climáticos simulados a partir de ADAPTECCA. 
Eligo el escenario que quiero HIST es histórico, 4.5 el RCp4.5 y 8.5 el RCP 8.5 (el peor)*/

/* Importante. Los datos climaticos deben cubrir hasta el último año de simulación*/

/* Seleccionamos el escenario climático*/

data clima; set fterra.escenarios_clima;

if escenario = "HIST";
/*if escenario = "4_5";
if escenario = "8_5";*/

run;


/* generamos la base de datos de regenerado. Cada vez que el rodal termine cortas y se regenere lo que se regenera es esto*/

data reg_ini;set fterra.reg_ini;

;
parcela=1; aux=1; SI=17; pon=1; year=2010; t=20; run;

data reg_ini; set reg_ini; 
qq = rannor(1)*sqrt(0.0088);d=pon*d; run; 

/* Importo las cabeceras de las tablas que luego voy a obtener*/
data tabla; set fterra.tabla;run;
data cortatot; set fterra.corta;run; 
data autoaclatot; set fterra.autoaclareo;run; 


/*A continuación programo la macro estado, que es la que me simula el estado de la parcela en cada año de simulacion*/

%macro estado;

dm log 'clear';

/*Genero una variable year, que se corresponde en cada ciclo con el año de simulación*/

data datos; set datos; year=&year;  run;

/*Junto la base de datos correspondiente a ese año con los datos climáticos*/

data datos; merge datos clima; by year; run;

data datos;set datos; if d^='.'; run;


/* Si queremos hacer una simulación asumiendo todos los años un clima medio incluiría el siguiente código*/

/*data datos; set datos; 
NVSTD=0;	mrstd=0;	spstd=0;	jl3std=0;	jl7std=0;	helstd=0;
run;*/


/* calcular la densidad, a partir del conteo del numero de observaciones */

proc means data=datos n noprint; var arbol; output out=salida n=dens;RUN;
data salida; set salida; keep dens; run;
data salida; set salida; aux=1; run;

data datos; merge datos salida; by aux;run;

/* Calculamos los percentiles 90 (=diam dominante)y 10 de la distribucion de diametros */
proc means data=datos noprint p90 p10 ; var d; output out=perc p90=d90 p10=d10 ; run:
data perc; set perc; p9010=d90-d10;run;
data perc; set perc;keep d90 d10 p9010;run;
data perc; set perc;aux=1; run;
data datos;merge datos perc; by aux;run;
/* calculo de área basimétrica: BA, diam medio cuadrático: dg, índice Reineke: SDI, cociente entre d y dg: ddg, altura dominante ho */
data datos; set datos; gi=3.14159*(d**2)/40000;
proc means noprint data=datos sum; var gi; output out= ab1 sum=ab ;run;
data ab1; set ab1;aux=1;run;
data ab1; set ab1; keep ab aux;run;
data datos; merge datos ab1;by aux;run;
data datos; set datos;
dg=sqrt(40000*(ab/dens)/3.14159);
SDI=dens*(dg/25)**1.605;
ddg=d/dg;
ho=exp(4.1437+((log(si)-4.1437)*((t/100)**-0.3935)));
dn=d90;
hn=(ho);
ldens=log(dens);
logd=log(d);
logab=log(ab);

/* Cálculo de las alturas total de cada uno de los árboles*/

n1=5.5862+(-0.4563*ldens);
m11=n1*((1/d)-(1/dn));
k1= m11+ ((1/(hn-1.3))**0.5) ;
cho1= k1**(-2);
h=1.3+(cho1) ;


/* Calculo de la produccion anual de piña de cada árbol m1*/
/* UNH! y UNH2 son específicos de la Unidad Natural Viana. Para otras unidades consultar Calama et al. 2016. Forest Systems https://revistas.inia.es/index.php/fs/article/view/9671 */

UNH1=-1.1506; UNH2=-0.3565;

logist1=0.3997+0.7697*nvstd+0.6209*mrstd+0.3389*spstd+0.3956*jl3std-0.3237*helstd+unh1-58.0943*(1/dg)-0.0020*sdi+2.3315*ddg;

abund1 = -0.8654 + 0.4772*nvstd+0.0992*mrstd+0.2780*spstd+0.2496*jl7std-25.4095*(1/dg)-0.0008*sdi+1.0931*ddg + 0.0322*d;

m1 = (exp(logist1)/(1+exp(logist1)))*exp(0.127)*exp(abund1);

run;

/* Creo una variable que es producción de piña aprovechable, que es m1 si vale más de 1 kg, y 0 en caso contrario.
Estoya asumiendo que pies con menos de 1 kg de piña no se cosechan. Este límite podría ser cambiado por el usuario*/

data datos; set datos; if m1<1 then m2=0;else m2=m1;run;

/* Cálculo de las dimensiones de copa de cada pie*/

data datos; set datos;

hbc = h*exp((-0.1254*(d/h))+(-11.07/t)+(-2.9504*(d/h)*(1/t)));
cw = 0.813 -0.202*hbc+0.169*d;
run;

/* Modulo transicion = calculo del incremento en diámetro (id) que el pie va a experimentar este año.
Empleamos modelo de BAI de Calama et al. 2019*/

data id;set datos;

logbai = qq + 1.6319+0.07028*d-0.0004*(d*d)-0.07557*Ho-0.2937*logab+0.06172*si+0.2009*ppstd-0.1186*tstd;
bai=exp(logbai);
 
run;
data id; set id;
id1= -d + sqrt((d*d)+4*(bai/3.141592));
data id;set id; if id1<0 then id1=0;run;
data id; set id;if d>22.5 then zlim=1; else zlim=0;
run;

/* Ecuaciones de biomasa. Ruiz- Peinado et al. 2012*/
data id;set id;
bf=0.0224*(d**1.923)*(h**1.0193);
br7=zlim*(0.247*((d-22.5)**2));
br27=0.0525*(d**2);
br2=21.927+0.0707*(d**2)-2.827*h;
braiz=0.117*(d**2);
ba= bf+br7+br27+br2;
bt=ba+braiz;

run;

/* Calculo del volumen y función perfil. Aquí entra la macro cubica, que está progrmada más abajo*/


%cubica;

/*data idvol;set idvol;
if h<2 then vmad=0; if hpred2<2 then  vtroza=0;
if hpred2<2 then vsierra=0; if hpred2<2 then  vtrit=0;run;*/


/* Genero la linea de la tabla de valores resumenes de salida correspondiente al año y parcela*/

/* Aquí los valores medios: ab, dg, h, t, dens...*/
 
proc means data=idvol noprint mean;var parcela year  dg d90 ho h t ab dens SDI vmad m1;
output out= medias mean=parcela year dg ddom ho hm t ab dens SDI vmad pina;run;
/* aquílos que son suma: volúmenes, biomasa, piña*/


proc means data=idvol noprint sum; var m1 vtroza vmad vsierra vtrit m2 bf br7 br27 br2 ba bt braiz;

output out= sumas sum= pina_ha vtroza Vtotal Vsierra Vtrit pina_aprov bf br7 br27 br2 ba bt braiz; run;


data tabla&year; merge medias sumas;   run;

/*Añade la fila del año a la del año anterior*/

data tabla; set tabla tabla&year;run;

data datos1; set idvol;
run;

/* Módulo de paso. Mortalidad por autoaclareo y definición de claras y cortas*/

/*Genero un número aleatorio para seleccionar los pies que mueren de forma natural*/
/* Asumo que si SDI > 500 se muere un 2% de la masa. Puedo cambiarlo*/

data datos1; set datos1; 
autac=ranuni(0);
data autoaclareo; set datos1; if sdi>=500 and autac>0.98; run; 
proc means data=autoaclareo noprint sum; var vtroza vmad vsierra vtrit bf br7 br27 br2 ba bt braiz;

output out= autoaclareomedia sum= vtroza Vtotal Vsierra Vtrit bf br7 br27 br2 ba bt braiz; run;

proc means data=autoaclareo noprint mean; var parcela year;  
output out= autoaclamedia mean=parcela year ;run; 

data autoaclareo&year; merge autoaclamedia autoaclareomedia; run;

data autoaclatot; set autoaclatot autoaclareo&year; run; 

data datos1; set datos1; if sdi>=500 and autac>0.98 then delete; 


/* aquí metemos las definición de las claras y cortas a realizar, salvo la corta final. 
Las claras y cortas de regeneración quedan definidas por una edad (t) y la densidad máxima a esa edad.
En el caso de abajo plantemos una clara a los 30, quedando 250 pies/ha, otra a los 45, quedando 150 pies/ha, una corta a los 100 donde quedan
75 pies /ha y otra a los 110 con 25 pies/ha*/
/* Si quieres simular más cortas añadir más líneas, siempre en orden de edad creciente y densidad decreciente*/
/*La corta final queda definida en la linea donde se invoca al data corta, definiendo que si t=turno se corta todo lo que queda*/



proc sort data=datos1; by descending d; 
data datos1; set datos1;
ord=_n_;

data datos2; set datos1; if t>=30 and ord>250 then marcado=1; else marcado=0; run;
data datos2; set datos2; if t>=45 and ord>150 then marcado=marcado+1; else marcado=marcado;run;
data datos2; set datos2; if t>=100 and ord>75 then marcado=marcado+1; else marcado=marcado;run; 
data datos2; set datos2; if t>=110 and ord>25 then marcado=marcado+1; else marcado=marcado;run; 

data datos1; set datos2; if marcado=0;run;
data corta; set datos2; if marcado>0 or t>=120;
proc means data=corta noprint mean; var parcela year;  
output out= cortamedias mean=parcela year ;run;

/*Genero la fila resumen con lo que se ha cortado ese año y parcela*/

proc means data=corta noprint sum; var vtroza vmad vsierra vtrit bf br7 br27 br2 ba bt braiz;

output out= cortas sum= vtroza Vtotal Vsierra Vtrit bf br7 br27 br2 ba bt braiz; run;

data corta&year; merge cortamedias cortas; run; 

data cortatot; set cortatot corta&year;run;

 

/* A partir de aquí genero la masa que va a entrar en simulación el año que viene, una vez quitado lo aclarado o muerto*/
/* En el caso de que haya legado a la corta final, para el año siguiente se simula un nuevo regenerado con 20 años de edad*/

/* Simulacion del crecimiento */
data datos1;set datos1;
d=d+id1;
t=t+1;
year=year+1;
run;

data datos1;set datos1;

keep parcela year arbol d t si qq aux ;

run;

proc means data=datos1 noprint mean; var t;
output out=cortafinal mean=t_par;run;

data cortafinal;	set cortafinal; aux=1; 


run;

data reg; set reg_ini; parcela=&parcela;run;
data datos1; set datos1 reg; run;

data datos1; merge datos1 cortafinal; by aux; run;

data datos1; set datos1;
if t_par>=121 then do;
if t>=121 then delete; 
end; 

else do;
if t=20 then delete;
end;
run;

data datos; set datos1;  run;

%mend estado;

/*Fin de la macro estado*/



/* macro CUBICA para incluir en la macro estado*/

/* Esta macro se basa en la ecuación de perfil de Calama et al. 2006.
El programa secciona la altura total del árbol en rodajas de 20 cm, calcula el diámetro de cada sección, y la cubica*/

%macro cubica;

proc sort data=id; by arbol;run;
data volumen; set id; keep arbol t d h dens ldens dg ho ab;run;

data volumen; set volumen;
h1=h;
di_dg1=d-dg;run;
data volumen; set volumen;
troza=1+round(h/0.20);run;


data volumen; set volumen;
do j=1 to troza;
output;

end;
run;
data volumen; set volumen;
h=j*0.20;run;
data volumen; set volumen; if h<h1;run;

/* La ecuación de perfil viene en dm (altura) y mm(diámetros)*/
data volumen; set volumen;
d1=d*10;
dg1=dg*10;
h1=h1*10;
h=h*10;
di_dg1=di_dg1*10;

run;
data volumen; set volumen;
htt=h1**3.5;
x0= (h1-h)/(h1-13);
x1= ((h1**1.5)-(h**1.5))*(h-13)/(h1**1.5);
x2= ((h1-h)**4)*(h-13)/(htt);
run;
data volumen; set volumen;
if h<=13 then

dSEC= (d1*x0)+ ((7.4085-(0.7306*ldens)+(0.0065*di_dg1))*x1)+(-0.4176*x2);

if h>13 then

dSEC= (d1*x0)+ ((0.4034+(0.05907*ho)+(0.0033*di_dg1))*x1)+((0.1417-(0.0093*ab))*x2);
run;
data volumen; set volumen;
sec=3.141516*(dSEC**2)/4000000;
arbol1=lag1(arbol);
sec1=lag1(sec);
secm=(sec+sec1)/2;run;
data volumen; set volumen;
if sec>sec1 then secm='.';
run;
data volumen; set volumen;
vtroza=secm*0.2;
run;

/* En esta línea defino el diámetro mínimo de sierra y el maderable en punta delgada. En este caso 300 mm y 70 mm.Podría cambiarlo*/

data volumen; set volumen;
if dSEC>=70 then vmad=vtroza; else vmad=0;
if dSEC>=300 and d1>=300 then vsierra=vtroza; else vsierra=0;
vtrit=vmad-vsierra;
run;
proc means data =volumen mean sum noprint; var vtroza vmad vsierra vtrit; by arbol;
output out=volarbol sum= vtroza vmad vsierra vtrit;run;

data idvol; merge id volarbol;
by arbol;run;

%mend cubica;

/* Aquí acaba la macro cubica*/


/* EMpieza la macro simula. Aquí voy a indicar los parámetros de simulación*/


%macro simula;

/* Año inicial y final de simulación*/

%do year = 2010 %to 2049;

%estado;

%end;

%mend simula;

/*
data datos; set datos_ini; if parcela=2;

%simula; */

/* Aquí indico que parcela (clase de edad) quiero simular*/

%macro bucle;
%do parcela = 1 %to 2;
data datos; set datos_ini; if parcela=&parcela;run;
%simula;
%end;
%mend bucle;

%bucle;

/* Genero las tablas de cortas, autoaclareo y tabla producción*/

data fterra.cortas_viana_prueba; set cortatot;run;
data fterra.tabla_viana_prueba; set tabla;run;
data fterra.autoac_viana_prueba; set autoaclatot;run;

PROC EXPORT  DATA=work.cortatot
OUTFILE='/home/u138706/RAFA/FTRRA/corta_PRUEBA_17'
DBMS=xlsx replace;  run;

PROC EXPORT DATA=work.tabla
OUTFILE='/home/u138706/RAFA/FTRRA/tabla_PRUEBA_17' 
DBMS=xlsx replace;  run;

PROC EXPORT DATA=work.autoaclatot
OUTFILE='/home/u138706/RAFA/FTRRA/autoac_PRUEBA_17' 
DBMS=xlsx replace;  run;

%mend viana_BIOMASA;


/* Si corro todo lo anterior carga el programa. Luego solo tengo que correr el comando siguiente para que se ejecute*/


%viana_BIOMASA;
