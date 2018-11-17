; get parameter setting for reference run
pref=(/2.,2.,1000.,1.,550.,1./)
nparam=dimsizes(pref)

;also read in index for each parameter range
pindex=new(nparam,integer)

rname=stringtochar(runname)
lenc=sizeof(rname)-1
tmp=tostring(rname(lenc-2:lenc-1))
pindex(5)=tointeger(tmp)
tmp=tostring(rname(lenc-5:lenc-4))
pindex(4)=tointeger(tmp)
tmp=tostring(rname(lenc-8:lenc-7))
pindex(3)=tointeger(tmp)
tmp=tostring(rname(lenc-11:lenc-10))
pindex(2)=tointeger(tmp)
tmp=tostring(rname(lenc-14:lenc-13))
pindex(1)=tointeger(tmp)
tmp=tostring(rname(lenc-17:lenc-16))
pindex(0)=tointeger(tmp)
delete(rname)
delete(tmp)

; then read in parameter vectors
p=new(nparam,float)
do n=0,nparam-1
  pall=readAsciiTable ( pathin+"/"+"param_"+n,1,"float",1)
  p(n)=pall(pindex(n),0)
  delete(pall)
end do

;compute perturbation parameter vector
p=p-pref

rname=runname


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

;get reference profile in the center of paramater space
filename="mminput.nc"
fref=addfile(pathin+"/"+casename+"_01_01_01_01_01_01/"+filename,"r")

rmse = new((/ntime,numvar/),float)
rmse=0.0
rmse_tot=0.

; do loop over variables
do ivar=0,dimsizes(variables)-1

;load reference profile
vref=fref->$variables(ivar)$


filein=casename + "_" + variables(ivar) + "_mmcoeff.nc"


if (ivar.eq.4 .or. ivar .eq. 5) then
  f=addfile("/global/scratch/langhans/metamodel_n"+"/"+filein,"r")
else
  f=addfile(pathin+"/"+filein,"r")
end if

a=f->a
b=f->b
time=f->time
ntime=dimsizes(time)

if (ivar.eq.4 .or. ivar .eq. 5) then
  filein=casename + "_" + "PLKMAX" + "_mmcoeff.nc"
  f=addfile(pathin+"/"+filein,"r")
  akmax=f->a
  bkmax=f->b
  vrefkmax=fref->PLKMAX
end if



if (casename.eq."BOMEX") then
  rind=0
else
  rind=1
end if

kind=0

;do time loop
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

 if (variables(ivar).ne."THGRAD") then

 if (ivar.eq.4.or.ivar.eq.5) then
   fkmax = vrefkmax(nt)
   do i=0,nparam-1
     fkmax = fkmax + akmax(nt,0,i)*p(i) + bkmax(nt,0,i,i)*p(i)^2
   end do
   do i=0,nparam-1
   do j=0,nparam-1
     if (j.ne.i) then
       fkmax = fkmax + bkmax(nt,0,j,i)*p(i)*p(j)
     end if
   end do
   end do

   if (casename.eq."BOMEX".and.fkmax.lt.2.0) then
      fkmax=2.0
   end if
   if (casename.eq."CPBL2".and.fkmax.lt.5.0) then
      fkmax=5.0
   end if

   if (fkmax.gt.42..and.casename.eq."BOMEX") then
     fkmax=42.
   end if

   if (casename.ne."BOMEX") then
    if (nt.eq.0.and.fkmax.gt.102) then
      fkmax=102.
    end if
    if (nt.eq.1.and.fkmax.gt.111) then
      fkmax=111.
    end if
   end if


 else
   fkmax=1000000.

 end if

 print(fkmax)

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

 
 else   ; variables is thgrad

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

 rmse(nt,ivar)=rmse(nt,ivar)+abs(var-varb)

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

if (casename.eq."BOMEX") then

system("echo "+ sprintf("%13.8f",rmse_tot) + " "+sprintf("%13.8f",avg(rmse(:,0))) + " "+sprintf("%13.8f",-999.) + " "+ sprintf("%13.8f",avg(rmse(:,1))) \
  +" "+sprintf("%13.8f",avg(rmse(:,6)))+" "+sprintf("%13.8f",avg(rmse(:,2)))+" "+sprintf("%13.8f",avg(rmse(:,3))) \ 
 +" "+sprintf("%13.8f",avg(rmse(:,4)))+" "+sprintf("%13.8f",avg(rmse(:,5)))+" > " +pathin+"/"+rname+"/mmrmse_thvd_"+casename)

else

system("echo "+ sprintf("%13.8f",rmse_tot) + " "+sprintf("%13.8f",avg(rmse(:,0))) + " "+sprintf("%13.8f",avg(rmse(:,6))) + " "+ sprintf("%13.8f",avg(rmse(:,1))) \
  +" "+sprintf("%13.8f",-999.)+" "+sprintf("%13.8f",avg(rmse(:,2)))+" "+sprintf("%13.8f",avg(rmse(:,3))) \
 +" "+sprintf("%13.8f",avg(rmse(:,4)))+" "+sprintf("%13.8f",avg(rmse(:,5)))+" > " +pathin+"/"+rname+"/mmrmse_thvd_"+casename)

end if
