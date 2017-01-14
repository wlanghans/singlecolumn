;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; get the benchmark profile
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
f=addfile(benchdir+"/benchmark.nc","r")
zmean=f->height1
zflx=f->height2
nzmean = dimsizes(zmean)
nzflx = dimsizes(zflx)

qt =  f->QT
thl = f->THL
if (casename.eq."BOMEX") then
qc  = f->QC
wthl=f->THLflx
wqt=f->QTflx
end if


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;get the 1D model results 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
f=addfile(pathin+"/"+filein,"r")

tavg=-1.
if (casename.eq."CPBL2".or.casename.eq."CPBL4") then
tavg=0.
t_ind=(/3.,5./)*3600.
end if
if (casename.eq."BOMEX") then
tavg=1800.
t_ind=(/2.5/)*3600.
end if

if (tavg.eq.-1.) then
  print("ERROR: Case not specified in NCL postproc!")
  exit
end if

z=(f->z)
zi=(f->z2)
time=f->time
nzm=dimsizes(z)

thl1D = new((/dimsizes(t_ind),nzm/),float)
qt1D = new((/dimsizes(t_ind),nzm/),float)

if (casename.eq."BOMEX") then
qc1D = new(nzm,float)
wthl1D = new(nzm+1,float)
wqt1D = new(nzm+1,float)
rhoi1D= new(nzm+1,float)
end if

do t=0,dimsizes(t_ind)-1
    tmin = ind(abs(time-t_ind(t)+tavg).eq.min(abs(time-t_ind(t)+tavg)))
    tmax = ind(abs(time-t_ind(t)-tavg).eq.min(abs(time-t_ind(t)-tavg)))
    nt  =tmax-tmin+1
    avgind = ispan(tmin,tmax,1)
    thl1D(t,:)=(dim_avg_n(f->thl(tmin:tmax,:),0))
    qt1D(t,:)=(dim_avg_n(f->qt(tmin:tmax,:),0))
  
    if (casename.eq."BOMEX") then
    rhoi1D(1:nzm-1)  = ((0.5*(dim_avg_n(f->rho(avgind,0:nzm-2),0) + dim_avg_n(f->rho(avgind,1:nzm-1),0))))
    rhoi1D(0) = 2.*rhoi1D(1) - rhoi1D(2)
    rhoi1D(nzm)= 2.*rhoi1D(nzm-1) - rhoi1D(nzm-2)
    qc1D(:)=(dim_avg_n(f->qcl(avgind,:),0))
    wthl1D(:)=(dim_avg_n(f->tflx(avgind,:),0) )
    wqt1D(:)=(dim_avg_n(f->qtflx(avgind,:),0)) 
    end if
end do

rmse_thl = new(dimsizes(t_ind),float)
rmse_qt  = new(dimsizes(t_ind),float)
rmse_thl = 0.
rmse_qt  = 0.

; convert fluxes to W/m2 if BOMEX case
if (casename.eq."BOMEX") then
  wthl(0:nzflx-1)=1005. * wthl(0:nzflx-1) * rhoi1D(1:nzflx)
  wqt(0:nzflx-1)=2.5e06 * wqt(0:nzflx-1) * rhoi1D(1:nzflx)
  rmse_qc   = 0.
  rmse_wthl = 0.
  rmse_wqt  = 0.
end if

do n=0,dimsizes(t_ind) - 1 

ncount = 0
do k=0,nzmean-1

  if (casename.eq."BOMEX") then
    if (zmean(k).lt.2900.) then
      ncount=ncount+1
      rmse_thl(n) = rmse_thl(n) + (thl(n,k)-thl1D(n,k+1))^2. 
      rmse_qt(n)  = rmse_qt(n)  + (qt(n,k)-qt1D(n,k+1))^2. 
    end if
  end if
  if (casename.eq."CPBL2") then
    rmse_thl(n) = rmse_thl(n) + (thl(n+1,k)-thl1D(n,k+1))^2. 
    rmse_qt(n)  = rmse_qt(n)  + (qt(n+1,k)-qt1D(n,k+1))^2. 
  end if
  if (casename.eq."CPBL4") then
    rmse_thl(n) = rmse_thl(n) + (thl(n+3,k)-thl1D(n,k+1))^2. 
    rmse_qt(n)  = rmse_qt(n)  + (qt(n+3,k)-qt1D(n,k+1))^2. 
    if (n.eq.0) then
      ;print(zmean(k)+" "+thl(n+3,k) + " "+thl1D(n,k+1)+" "+(thl(n+3,k)-thl1D(n,k+1))^2.)
    end if
  end if

end do
end do

if (casename.eq."BOMEX") then
ncount=0
do k=0,nzmean-1

  if (zmean(k).lt.2900.) then
  ncount=ncount+1
  rmse_wthl = rmse_wthl + (wthl(k)-wthl1D(k+1))^2. 
  rmse_wqt  = rmse_wqt  + (wqt(k)-wqt1D(k+1))^2. 
  rmse_qc   = rmse_qc   + (qc(k)-qc1D(k+1))^2.
  end if

end do
  rmse_wthl=sqrt(rmse_wthl/ncount)
  rmse_wqt=sqrt(rmse_wqt/ncount)
  rmse_qt=sqrt(rmse_qt/ncount)
  rmse_qc=sqrt(rmse_qc/ncount)
  rmse_thl(:) = sqrt(rmse_thl(:)/ncount)
  rmse_qt(:)  = sqrt(rmse_qt(:)/ncount)
else
  rmse_thl(:) = sqrt(rmse_thl(:)/nzmean)
  rmse_qt(:)  = sqrt(rmse_qt(:)/nzmean)
end if

if (casename.eq."BOMEX") then
  thl_0 = 300. 
  qt_0  = 14.e-3
  qc_0  = 3.e-6
  wthl_0=8.
  wqt_0 =150.
  rmse= 1./6.*(rmse_thl(0)/thl_0 + rmse_qt(0)/qt_0) + \
  2./9.*(rmse_wqt/wqt_0 + rmse_wthl/wthl_0+rmse_qc/qc_0)
else
  thl_0 = 300. 
  qt_0  = 5.e-3
  rmse = 0.25 * (sum(rmse_thl)/thl_0 + sum(rmse_qt)/qt_0) 
end if


; write case rmse's to run directory
if (casename.eq."CPBL2".or.casename.eq."CPBL4") then
  rmse_wqt =-999.
  rmse_wthl=-999.
  rmse_qc  =-999.
  ;print("RMSE THL="+sprintf("%13.8f",avg(rmse_thl(:))))
end if

system("echo "+ sprintf("%13.8f",rmse) + " "+sprintf("%13.8f",avg(rmse_thl(:))) + " "+ sprintf("%13.8f",avg(rmse_qt(:))) \
+" "+sprintf("%13.8f",rmse_qc)+" "+sprintf("%13.8f",rmse_wthl)+" "+sprintf("%13.8f",rmse_wqt) +" > " +pathin+"/rmse_"+casename)

