;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;get the 1D model results 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (var.ne."QTD" .and. var.ne."THVD") then
filename="mminput.nc"
f=addfile(pathin+"/"+casename+"_01_01_01_01_01_01/"+filename,"r")
z=(f->z)
zi=(f->zi)
time=f->time
vref=f->$var$
dims=dimsizes(vref)
if (dimsizes(dims).gt.1) then
nz=dims(1)
else
nz=1
end if


else


filename="rmse_thvd_"+casename
rmses = asciiread(pathin+"/"+casename+"_01_01_01_01_01_01/"+filename,(/1,9/),"float")
if (var.eq."QTD") then
  vref=rmses(0,8)
  var_ind=8
else
  vref=rmses(0,7)
  var_ind=7
end if
nz=1
time=0.0
dims=dimsizes(vref)

end if


nump=6

; read parameter coordinates
pc=new((/nump,3/),float)
do i=0,nump-1
  pc(i,:)=readAsciiTable (pathin+"param_"+i,1,"float",1)
end do

do i=0,nump-1
 pc(i,:)=pc(i,:)-pc(i,1)
end do


;allocate a and B 
ac=new((/dimsizes(time),nz,nump/),float)
bc=new((/dimsizes(time),nz,nump,nump/),float)
ac=-999.
bc=-999.
ac@_FillValue=-999.
bc@_FillValue=-999.

; time loop
do nt=0,dimsizes(time)-1

if (nz.gt.1) then

ind_min=nz

if (num(ismissing(vref(nt,:))) .gt. 0 ) then 
ind_miss=ind(ismissing(vref(nt,:)))
ind_n = min(ind_miss)
delete(ind_miss)
if (ind_n.lt.ind_min) then
  ind_min=ind_n
end if
end if


end if ; nz .gt. 1

; do parameter loop to get ai and Bii
do i=0,nump-1

; read min, max profiles for parameter i
rnamemin=casename
rnamemax=casename
do j=0,nump-1
  rnamemin=rnamemin+"_"
  rnamemax=rnamemax+"_"
  if (j.eq.i) then
  rnamemin=rnamemin+"00"
  rnamemax=rnamemax+"02"
  else
  rnamemin=rnamemin+"01"
  rnamemax=rnamemax+"01"
  end if
end do

if (var.ne."QTD" .and. var.ne."THVD") then

f=addfile(pathin+"/"+rnamemin+"/"+filename,"r")
vmin=f->$var$
f=addfile(pathin+"/"+rnamemax+"/"+filename,"r")
vmax=f->$var$

else

rmses = asciiread(pathin+"/"+rnamemin+"/"+filename,(/1,9/),"float")
vmin=rmses(0,var_ind)

rmses = asciiread(pathin+"/"+rnamemax+"/"+filename,(/1,9/),"float")
  vmax=rmses(0,var_ind)

end if

if (nz.gt.1) then

if (num(ismissing(vmin(nt,:))) .gt. 0 ) then 
ind_miss=ind(ismissing(vmin(nt,:)))
ind_n = min(ind_miss)
delete(ind_miss)
if (ind_n.lt.ind_min) then 
  ind_min=ind_n
end if
end if

if (num(ismissing(vmax(nt,:))) .gt. 0 ) then 
ind_miss=ind(ismissing(vmax(nt,:)))
ind_n = min(ind_miss)
delete(ind_miss)
if (ind_n.lt.ind_min) then 
  ind_min=ind_n
end if
end if


; do level loop
do l=0,ind_min-1

; compute Bii
bc(nt,l,i,i)=(vmin(nt,l)-vref(nt,l)-pc(i,0)/pc(i,2)*(vmax(nt,l)-vref(nt,l)))/(pc(i,0)*(pc(i,0)-pc(i,2)))

; compute aii
ac(nt,l,i)=(vmax(nt,l)-vref(nt,l)-pc(i,2)^2 * bc(nt,l,i,i) )/pc(i,2)

end do ; level loop

else   ;nz=1

; compute Bii
bc(nt,0,i,i)=(vmin(nt)-vref(nt)-pc(i,0)/pc(i,2)*(vmax(nt)-vref(nt)))/(pc(i,0)*(pc(i,0)-pc(i,2)))

; compute aii
ac(nt,0,i)=(vmax(nt)-vref(nt)-pc(i,2)^2 * bc(nt,0,i,i) )/pc(i,2)


end if

end do ; parameter loop

; do pair-wise loop to get interaction terms Bij
do i=0,nump-1
do j=0,nump-1

; only one side of the matrix needed
if (j.gt.i) then

; load 4 profiles for this parameter pair
namevec=(/"01","01","01","01","01","01"/)
n00=namevec
n00(j)="00"
n00(i)="00"
n01=namevec
n01(j)="00"
n01(i)="02"
n10=namevec
n10(j)="02"
n10(i)="00"
n11=namevec
n11(j)="02"
n11(i)="02"

if (var.ne."QTD" .and. var.ne."THVD") then

f=addfile(pathin+"/"+casename+"_"+n00(0)+"_"+n00(1)+"_"+n00(2)+"_"+n00(3)+"_"+n00(4)+"_"+n00(5)+"/"+filename,"r")
v00=f->$var$
f=addfile(pathin+"/"+casename+"_"+n01(0)+"_"+n01(1)+"_"+n01(2)+"_"+n01(3)+"_"+n01(4)+"_"+n01(5)+"/"+filename,"r")
v01=f->$var$
f=addfile(pathin+"/"+casename+"_"+n10(0)+"_"+n10(1)+"_"+n10(2)+"_"+n10(3)+"_"+n10(4)+"_"+n10(5)+"/"+filename,"r")
v10=f->$var$
f=addfile(pathin+"/"+casename+"_"+n11(0)+"_"+n11(1)+"_"+n11(2)+"_"+n11(3)+"_"+n11(4)+"_"+n11(5)+"/"+filename,"r")
v11=f->$var$

else

rmses = asciiread(pathin+"/"+casename+"_"+n00(0)+"_"+n00(1)+"_"+n00(2)+"_"+n00(3)+"_"+n00(4)+"_"+n00(5)+"/"+filename,(/1,9/),"float")
v00=rmses(0,var_ind)
rmses = asciiread(pathin+"/"+casename+"_"+n01(0)+"_"+n01(1)+"_"+n01(2)+"_"+n01(3)+"_"+n01(4)+"_"+n01(5)+"/"+filename,(/1,9/),"float")
v01=rmses(0,var_ind)
rmses = asciiread(pathin+"/"+casename+"_"+n10(0)+"_"+n10(1)+"_"+n10(2)+"_"+n10(3)+"_"+n10(4)+"_"+n10(5)+"/"+filename,(/1,9/),"float")
v10=rmses(0,var_ind)
rmses = asciiread(pathin+"/"+casename+"_"+n11(0)+"_"+n11(1)+"_"+n11(2)+"_"+n11(3)+"_"+n11(4)+"_"+n11(5)+"/"+filename,(/1,9/),"float")
v11=rmses(0,var_ind)

end if

if (.not.isvar("x")) then
   x=new(4,float)
end if
x(0) = pc(j,0) * pc(i,0)
x(1) = pc(j,0) * pc(i,2)
x(2) = pc(j,2) * pc(i,0)
x(3) = pc(j,2) * pc(i,2)

if (.not.isvar("y")) then
   y=new(4,float)
end if

if (nz.gt.1) then

if (num(ismissing(v00(nt,:))) .gt. 0 ) then 
ind_miss=ind(ismissing(v00(nt,:)))
ind_n = min(ind_miss)
delete(ind_miss)
if (ind_n.lt.ind_min) then
ind_min=ind_n
end if
end if
if (num(ismissing(v01(nt,:))) .gt. 0 ) then 
ind_miss=ind(ismissing(v01(nt,:)))
ind_n = min(ind_miss)
delete(ind_miss)
if (ind_n.lt.ind_min) then
ind_min=ind_n
end if
end if
if (num(ismissing(v10(nt,:))) .gt. 0 ) then 
ind_miss=ind(ismissing(v10(nt,:)))
ind_n = min(ind_miss)
delete(ind_miss)
if (ind_n.lt.ind_min) then
ind_min=ind_n
end if
end if
if (num(ismissing(v11(nt,:))) .gt. 0 ) then 
ind_miss=ind(ismissing(v11(nt,:)))
ind_n = min(ind_miss)
delete(ind_miss)
if (ind_n.lt.ind_min) then
ind_min=ind_n
end if
end if
;
; do level loop
do l=0,ind_min-1

; compute bij from linear regression using 4 corner points of the i,j plane
y(0)=v00(nt,l)-vref(nt,l) - ac(nt,l,j)*pc(j,0) - ac(nt,l,i)*pc(i,0) - bc(nt,l,j,j)*pc(j,0)^2 - bc(nt,l,i,i)*pc(i,0)^2
y(1)=v01(nt,l)-vref(nt,l) - ac(nt,l,j)*pc(j,0) - ac(nt,l,i)*pc(i,2) - bc(nt,l,j,j)*pc(j,0)^2 - bc(nt,l,i,i)*pc(i,2)^2
y(2)=v10(nt,l)-vref(nt,l) - ac(nt,l,j)*pc(j,2) - ac(nt,l,i)*pc(i,0) - bc(nt,l,j,j)*pc(j,2)^2 - bc(nt,l,i,i)*pc(i,0)^2
y(3)=v11(nt,l)-vref(nt,l) - ac(nt,l,j)*pc(j,2) - ac(nt,l,i)*pc(i,2) - bc(nt,l,j,j)*pc(j,2)^2 - bc(nt,l,i,i)*pc(i,2)^2

y=y/2.

;linear regression
bc(nt,l,j,i) = regline(x,y)
bc(nt,l,i,j) = bc(nt,l,j,i) 

end do ; end level loop

else ; nz=1

; compute bij from linear regression using 4 corner points of the i,j plane
y(0)=v00(nt)-vref(nt) - ac(nt,0,j)*pc(j,0) - ac(nt,0,i)*pc(i,0) - bc(nt,0,j,j)*pc(j,0)^2 - bc(nt,0,i,i)*pc(i,0)^2
y(1)=v01(nt)-vref(nt) - ac(nt,0,j)*pc(j,0) - ac(nt,0,i)*pc(i,2) - bc(nt,0,j,j)*pc(j,0)^2 - bc(nt,0,i,i)*pc(i,2)^2
y(2)=v10(nt)-vref(nt) - ac(nt,0,j)*pc(j,2) - ac(nt,0,i)*pc(i,0) - bc(nt,0,j,j)*pc(j,2)^2 - bc(nt,0,i,i)*pc(i,0)^2
y(3)=v11(nt)-vref(nt) - ac(nt,0,j)*pc(j,2) - ac(nt,0,i)*pc(i,2) - bc(nt,0,j,j)*pc(j,2)^2 - bc(nt,0,i,i)*pc(i,2)^2

y=y/2.

;linear regression
bc(nt,0,j,i) = regline(x,y)
bc(nt,0,i,j) = bc(nt,0,j,i)

end if

end if


end do
end do

if (nz.gt.1) then
if (ind_min.lt.nz) then
  ac(nt,ind_min:nz-1,:)   = -999.
  bc(nt,ind_min:nz-1,:,:) = -999.
end if
end if


end do ; time loop

; write matrix B and a out as netcdf for the chosen case
             
 system("/bin/rm -f "+pathin+"/"+casename+"_"+var+"_mmcoeff_rmse.nc" )   ; remove any pre-existing file
       ncdf = addfile(pathin+"/"+casename+"_"+var+"_mmcoeff_rmse.nc" ,"c")  ; open output netCDF file

    ;===================================================================
    ; create global attributes of the file (optional)
    ;===================================================================
       fAtt               = True            ; assign file attributes
       fAtt@title         = "Metamodel coefficients a,b for variable "+var
       fAtt@source_file   =  " "
       fAtt@Conventions   = "None"
       fAtt@creation_date = systemfunc ("date")
       fileattdef( ncdf, fAtt )            ; copy file attributes

    ;===================================================================
    ; make time an UNLIMITED dimension; recommended  for most applications
    ;===================================================================
       filedimdef(ncdf,"time",nt,False) 
       if (nz.gt.1) then
         filedimdef(ncdf,vref!1,nz,False) 
       else
         filedimdef(ncdf,"level",1,False) 
       end if
       filedimdef(ncdf,"p1",nump,False) 
       filedimdef(ncdf,"p2",nump,False) 

       ac!0="time"
       ac&time=time
       bc!0="time"
       bc&time=time
       if (nz.gt.1) then
         ac!1=vref!1
         bc!1=vref!1
         if (vref!1.eq."z") then
           ac&z=z
           bc&z=z
         end if
         if (vref!1.eq."zi") then
           ac&zi=zi
           bc&zi=zi
         end if
       else
         ac!1="level"
         ac&level=0
         bc!1="level"
         bc&level=0
       end if
       ac!2="p1"
       ac&p1=ispan(1,nump,1)
       bc!2="p1"
       bc&p1=ispan(1,nump,1)
       bc!3="p2"
       bc&p2=ispan(1,nump,1)
       ac@long_name = "linear coefficients"           ; assign attributes
       ac@_FillValue = -999.           ; assign attributes
       bc@long_name = "quadratic coefficients"           ; assign attributes
       bc@_FillValue = -999.           ; assign attributes

       ncdf->a  = ac
       ncdf->b  = bc
