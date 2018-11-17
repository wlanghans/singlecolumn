; get parameter setting for reference run
pref=(/2.,2.,1000.,1.,550.,1./)
nparam=dimsizes(pref)
pmin=(/0.,0.,0.,0.,100.,0./)
pmax=(/4.,4.,2000.,2.,1000.,2./)
p=new(nparam,float)
p1_out=new(nruns,float)
p2_out=new(nruns,float)
p3_out=new(nruns,float)
p4_out=new(nruns,float)
p5_out=new(nruns,float)
p6_out=new(nruns,float)

casenames=(/"BOMEX","CPBL2"/)
rmse_out=new((/dimsizes(casenames),nruns/),float)

random_setallseed(36484749+mpiid, 9494848+mpiid)
unf   = random_uniform(0., 1., (/nparam,nruns/))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; get the benchmark profile
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
f=addfile(benchdir+"/benchmark_SAM_BOMEXhighlarge_qt.nc","r")
zmean=f->height1
zflx=f->height2
nzmean = dimsizes(zmean)
nzflx = dimsizes(zflx)

qt =  f->QT
thl = f->THL
thvd=f->THVD
qtd=f->QTD

wthl=f->THLflx
wqt=f->QTflx


do case=0,dimsizes(casenames)-1

casename=casenames(case)
f=addfile(benchdir+"/benchmark_SAM_BOMEXhighlarge_qt.nc","r")

if (casename.eq."BOMEX") then

qc  = f->QC
variables=(/"THL","QT","WTHL","WQT","THVD","QTD","QC"/)
;variables=(/"THL","QT","WTHL","WQT","QC"/)
ntime=1

else

ntime=2
thgrad=new(2,float)
do i =0,1
tmp=thl(1+i,:)
zlow=1200.+i*200.
zhigh=zlow+400.
zlow=700.
tmp=where(zmean.gt.zhigh.or.zmean.lt.zlow,-999.,tmp)
tmp@_FillValue=-999.
thgrad(i)=regline(zmean,tmp)
delete(tmp)
end do

variables=(/"THL","QT","WTHL","WQT","THVD","QTD","THGRAD"/)
;variables=(/"THL","QT","WTHL","WQT","THGRAD"/)

end if

numvar=dimsizes(variables)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;get the metamodel results 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


do i0=0,nruns-1 ; loop over parameter space

p(0)=pmin(0) + unf(0,i0) * (pmax(0)-pmin(0))
p(1)=pmin(1) + unf(1,i0) * (pmax(1)-pmin(1))
p(2)=pmin(2) + unf(2,i0) * (pmax(2)-pmin(2))
p(3)=pmin(3) + unf(3,i0) * (pmax(3)-pmin(3))
p(4)=pmin(4) + unf(4,i0) * (pmax(4)-pmin(4))
p(5)=pmin(5) + unf(5,i0) * (pmax(5)-pmin(5))

p=p-pref

rmse = new((/ntime,numvar/),float)
rmse=0.0
rmse_tot=0.

; do loop over variables
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
do ivar=0,dimsizes(variables)-1

if (ivar.eq.4 .or. ivar .eq. 5) then
filename="rmse_thvd_"+casename
rmses = asciiread("/global/scratch/langhans/metamodel_mm3/"+"/"+casename+"_01_01_01_01_01_01/"+filename,(/1,9/),"float")
if (ivar.eq.5) then
  vref=rmses(0,8)
  var_ind=8
else
  vref=rmses(0,7)
  var_ind=7
end if

filein=casename + "_" + variables(ivar) + "_mmcoeff_rmse.nc"

else

;get reference profile in the center of paramater space
filename="mminput.nc"
fref=addfile("/global/scratch/langhans/metamodel_mm3/"+"/"+casename+"_01_01_01_01_01_01/"+filename,"r")
;load reference profile
vref=fref->$variables(ivar)$

filein=casename + "_" + variables(ivar) + "_mmcoeff_rmse.nc"
end if




f=addfile("/global/scratch/langhans/metamodel_mm3/"+"/"+filein,"r")

a=f->a
b=f->b
time=f->time

if (casename.eq."BOMEX") then
  rind=0
else
  rind=1
end if

kind=0

;do time loop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
do nt=0,ntime-1

 zlow=0.

 if (variables(ivar).eq."THL") then
    varb=thl(rind+nt,:)
    rnorm=0.298500898501127
    if (casename.eq."CPBL2") then
       zlow=10.
       rnorm=0.142648591929376
    end if
 end if
 if (variables(ivar).eq."QT") then
    varb=qt(rind+nt,:)
    rnorm = 4.709051607813674e-4
    if (casename.eq."CPBL2") then
       zlow=10.
       rnorm = 1.561694507888805e-004
    end if
 end if
 if (variables(ivar).eq."WTHL") then
    varb=wthl(rind+nt,:)
    kind=1
    if (casename.eq."BOMEX") then
      rnorm=8.29152631298648 
    else
      rnorm=5.55608093415176
    end if
 end if
 if (variables(ivar).eq."WQT") then
    varb=wqt(rind+nt,:)
    kind=1
    if (casename.eq."BOMEX") then
      rnorm=37.5370361230699
    else
      rnorm=18.0438175345439
    end if
 end if
 if (variables(ivar).eq."THVD") then
    varb=thvd(rind+nt,:)
    rnorm=0.110751867145755
    if (casename.eq."CPBL2") then
       zlow=10.
       rnorm=0.119297832160030
    end if
 end if
 if (variables(ivar).eq."QTD") then
    varb=qtd(rind+nt,:)
    rnorm=4.754294124718259e-4
    if (casename.eq."CPBL2") then
       zlow=10.
       rnorm=1.220929376408715e-4
    end if
 end if
 if (variables(ivar).eq."THGRAD") then
    varb=thgrad(nt)
    rnorm=9.579204282494366e-005
 end if
 if (variables(ivar).eq."QC") then
    varb=qc(0,:)
    rnorm=2.291507888805414e-6
 end if

 if (variables(ivar).ne."THGRAD".and.ivar.ne.4.and.ivar.ne.5) then

   fkmax=1000000.


 ikmax=floattointeger(floor(fkmax))

 ncount = 0.
 do k=0,nzmean-1

  if (zmean(k).gt.zlow.and.zmean(k).lt.2900..and.k.le.ikmax) then

  var = vref(nt,k+kind)
  do i=0,nparam-1
    var = var + a(nt,k+kind,i)*p(i) + b(nt,k+kind,i,i)*p(i)^2
  end do
  do i=0,nparam-1
  do j=0,nparam-1
    if (j.ne.i) then
      var = var + b(nt,k+kind,j,i)*p(i)*p(j)
    end if
  end do
  end do


  rmse(nt,ivar)=rmse(nt,ivar)+(var-varb(k))^2
  ncount=ncount+1.

  end if ;

 end do ; level loop

 if (ncount.eq.0.0) then
   rmse(nt,ivar)=0.
 else
   rmse(nt,ivar)=sqrt(rmse(nt,ivar)/ncount)
 end if

 
 else   ; variables is thgrad or rmse_qtd or rmse_thvd

 if (ivar.ne.4.and.ivar.ne.5) then
 var = vref(nt)
 do i=0,nparam-1
   var = var + a(nt,0,i)*p(i) + b(nt,0,i,i)*p(i)^2
 end do
 do i=0,nparam-1
 do j=0,nparam-1
   if (j.ne.i) then
     var = var + b(nt,0,j,i)*p(i)*p(j)
   end if
 end do
 end do

 else

 var = vref
 do i=0,nparam-1
   var = var + a(0,0,i)*p(i) + b(0,0,i,i)*p(i)^2
 end do
 do i=0,nparam-1
 do j=0,nparam-1
   if (j.ne.i) then
     var = var + b(0,0,j,i)*p(i)*p(j)
   end if
 end do
 end do

 end if

 if (ivar.eq.4.or.ivar.eq.5) then
   rmse(nt,ivar)=var
 else

   rmse(nt,ivar)=rmse(nt,ivar)+abs(var-varb)

 end if

 end if
   

end do


delete(varb)
delete(vref)
delete(a)
delete(b)
delete(time)

if(isvar("var")) then
  delete(var)
end if

rmse_tot=rmse_tot + 1./7. * avg(rmse(:,ivar))/rnorm

; end loop over variables
end do


;write to output
rmse_out(case,i0) = rmse_tot
if (case.eq.0.) then
p1_out(i0)=p(0)
p2_out(i0)=p(1)
p3_out(i0)=p(2)
p4_out(i0)=p(3)
p5_out(i0)=p(4)
p6_out(i0)=p(5)
end if
delete(rmse)
delete(rmse_tot)

end do ; loop over parameter space


end do ; loop over two cases

delete(unf)

;write output
mat=new((/nparam+1,nruns/),float)
mat(0,:) = 0.5*(rmse_out(0,:)+rmse_out(1,:))
mat(1,:) = p1_out
mat(2,:) = p2_out
mat(3,:) = p3_out
mat(4,:) = p4_out
mat(5,:) = p5_out
mat(6,:) = p6_out

print("Writing output to file " + fileout)
asciiwrite (fileout , mat)
