load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

bnd=(/0.1,0.2,0.3,0.4,0.5/)

if (casein.eq."CPBL2".or.casein.eq."CPBL4") then
  plt_times=3600.*(/3.,5./)
  tavg=600.
else
  plt_times=3600.*(/2.5/)
  tavg=1800.
end if

colors=(/"blue"/)
dashes=(/0,1,2,3,4/)
thickness=(/3.,1.5,1.5,1.5,1.5/)
plot=new(20,graphic)

wks   = gsn_open_wks ("eps","pltflx_"+casein)                   ; open workstation

gsn_define_colormap(wks,"BlGrYeOrReVi200")

f=addfile(filein,"r")

rd=287.
cp=1005.

z=f->z
zi=f->z2
nzm=dimsizes(z)
time=f->time
qt1D=new((/dimsizes(plt_times),dimsizes(z)/),float)
presi1D=new((/dimsizes(plt_times),dimsizes(zi)/),float)
thl1D=new((/dimsizes(plt_times),dimsizes(z)/),float)
thr1D=new((/dimsizes(plt_times),dimsizes(z)/),float)
thvd1D=new((/dimsizes(plt_times),dimsizes(bnd),dimsizes(zi)/),float)
thvd1D@_FillValue=-999.
thvd1D=0.
thgrad1D=new((/dimsizes(plt_times)/),float)

wqt1D=new((/dimsizes(plt_times),dimsizes(zi)/),float)
wthl1D=new((/dimsizes(plt_times),dimsizes(zi)/),float)
wqt1D_mf=new((/dimsizes(plt_times),dimsizes(zi)/),float)
wthl1D_mf=new((/dimsizes(plt_times),dimsizes(zi)/),float)
wqt1D_ed=new((/dimsizes(plt_times),dimsizes(zi)/),float)
wthl1D_ed=new((/dimsizes(plt_times),dimsizes(zi)/),float)
rhoi1D =new(dimsizes(zi),float)
BU=new((/3,dimsizes(plt_times),dimsizes(z)/),float)
WU=new((/3,dimsizes(plt_times),dimsizes(zi)/),float)

if (casein.eq."BOMEX") then
  qc1D=new((/dimsizes(plt_times),dimsizes(z)/),float)
end if
qtd1D=new((/dimsizes(plt_times),dimsizes(bnd),dimsizes(zi)/),float)
qtd1D@_FillValue=-999.
qtd1D=0.

n1=39
n2=24
n3=4
t_ind=new(dimsizes(plt_times),integer)
do i=0,dimsizes(t_ind)-1
  index_vec = ind(abs(plt_times(i)-time).eq.min(abs(plt_times(i)-time)))
  t_ind(i) = index_vec(0)
  is_vec=ind(abs(plt_times(i)-tavg-time).eq.min(abs(plt_times(i)+tavg-time)))
  ie_vec=ind(abs(plt_times(i)-time+tavg).eq.min(abs(plt_times(i)-time+tavg)))
  is=is_vec(0)
  ie=ie_vec(0)
  delete(index_vec)
  qt1D(i,:)  = dim_avg_n(f->qt(is:ie,:),0)
  presi1D(i,:) = dim_avg_n(f->pi(is:ie,:),0) 
  thl1D(i,:) = dim_avg_n(f->thl(is:ie,:),0)
  thr1D(i,:) = dim_avg_n(f->thr(is:ie,:),0)
  rhoi1D(1:nzm-1)  = ((0.5*(dim_avg_n(f->rho(is:ie,0:nzm-2),0) + dim_avg_n(f->rho(is:ie,1:nzm-1),0))))
  rhoi1D(0) = 2.*rhoi1D(1) - rhoi1D(2)
  rhoi1D(nzm)= 2.*rhoi1D(nzm-1) - rhoi1D(nzm-2)
  wqt1D(i,:)  = dim_avg_n(f->qtflx(is:ie,:),0)
  wthl1D(i,:)  = dim_avg_n(f->tflx(is:ie,:),0)
  wqt1D_mf(i,:)  = 2.5e06* dim_avg_n(f->qtflx_mf(is:ie,:),0) * rhoi1D
  wthl1D_mf(i,:)  = 1005. * dim_avg_n(f->tflx_mf(is:ie,:),0) * rhoi1D
  wqt1D_ed(i,:)  = 2.5e06 * dim_avg_n(f->qtflx_ed(is:ie,:),0) * rhoi1D
  wthl1D_ed(i,:)  = 1005. * dim_avg_n(f->tflx_ed(is:ie,:),0) * rhoi1D
  BU(0,i,:)  = dim_avg_n(f->B(is:ie,n1,:),0)
  WU(0,i,:)  = dim_avg_n(f->upw(is:ie,n1,:),0)
  BU(1,i,:)  = dim_avg_n(f->B(is:ie,n2,:),0)
  WU(1,i,:)  = dim_avg_n(f->upw(is:ie,n2,:),0)
  BU(2,i,:)  = dim_avg_n(f->B(is:ie,n3,:),0)
  WU(2,i,:)  = dim_avg_n(f->upw(is:ie,n3,:),0)
  if (casein.eq."BOMEX") then
    qc1D(i,:)  = dim_avg_n(f->qcl(is:ie,:),0)
  end if

  uptabs=f->upt(is:ie,:,:)
  upqcl=f->upqcl(is:ie,:,:)
  upw=f->upw(is:ie,:,:)
  upqt=f->upqt(is:ie,:,:)
  fraca=f->upa(is:ie,:,:)
  fractot=dim_sum_n(fraca,1)

  nt=ie-is+1

    do k=0,nzm-1
    navg=0
    do indi=0,nt-1
      if (fractot(indi,k).ge.0.01) then
         fracsmall=0.
         fracsmallold=0.
         np2=0
         lfrac=True
         do while (np2.le.39)
          fracsmall=fracsmallold+fraca(indi,np2,0)
          do nbnd=0,dimsizes(bnd)-1
            if (bnd(nbnd).ge.fracsmall) then
              wgt = fraca(indi,np2,k)
            else
              wgt = bnd(nbnd)-fracsmallold
              if (wgt.lt.0.0.or.fraca(indi,np2,k).eq.0.0) then
                 wgt = 0.0
              end if
            end if
          thvd1D(i,nbnd,k) = thvd1D(i,nbnd,k) + rhoi1D(k) *1005.* upw(indi,np2,k) * ((1000./presi1D(i,k))^(rd/cp)*(uptabs(indi,np2,k)-2.5e06/cp*upqcl(indi,np2,k))  -thl1D(i,k)) * wgt
          qtd1D(i,nbnd,k) = qtd1D(i,nbnd,k) +  rhoi1D(k) *2.5e06* upw(indi,np2,k) * (upqt(indi,np2,k)-qt1D(i,k)) * wgt
          end do
          fracsmallold=fracsmall
          np2=np2+1
         end do
         navg=navg+1
      end if
    end do
    if (navg.gt.0) then
       thvd1D(i,:,k) = thvd1D(i,:,k) / int2flt(navg)
       qtd1D(i,:,k) = qtd1D(i,:,k) / int2flt(navg)
    else
       thvd1D(i,:,k) = -999.
       qtd1D(i,:,k) = -999.
    end if
    end do

    
 
    delete(is_vec)
    delete(ie_vec)

end do

thvd1D(:,:,0)=-999.
qtd1D(:,:,0)=-999.


; get benchmark LES data

f=addfile("/global/home/users/langhans/singlecolumn/src_ensemble/benchmark_SAM_BOMEXhighlarge_qt.nc","r")

if (casein.eq."CPBL2") then
   Bt_ind = (/1,2/)
end if
if (casein.eq."CPBL4") then
   Bt_ind = (/3,4/)
end if
if (casein.eq."BOMEX")
   Bt_ind = (/0/)
end if

zles=f->height1
zflx=f->height2
nzles=dimsizes(zles)
nzflx=dimsizes(zflx)
thl=new((/dimsizes(Bt_ind),nzles/),float)
qt=new((/dimsizes(Bt_ind),nzles/),float)
wqt=new((/dimsizes(Bt_ind),nzflx/),float)
wthl=new((/dimsizes(Bt_ind),nzflx/),float)
thvd=new((/dimsizes(Bt_ind),nzles/),float)
thgrad=new((/dimsizes(Bt_ind)/),float)


thl=f->THL(Bt_ind,:)
qt =f->QT (Bt_ind,:)
wthl=f->THLflx(Bt_ind,:)
wqt=f->QTflx(Bt_ind,:)
thvd=f->THVD(Bt_ind,:)
qtd=new((/dimsizes(Bt_ind),nzles/),float)
qtd=f->QTD(Bt_ind,:)


 res                   = True                       ; plot mods desired
 res@gsnFrame = False
 res@gsnDraw  = False
 res@vpHeightF = 0.5
; res@tiMainString      = "Potential temperature"             ; add title
 str=floattointeger((201.-36.)/dimsizes(bnd)) + 1
 res@xyLineColors = ispan(36,201,str)
 res@xyLineThicknessF = 3.
 res@xyDashPattern = 0
 res@vpWidthF=0.6
 res@trYMinF=00.
 if (casein.eq."BOMEX") then
   res@trYMaxF=2000.
 else
   res@trYMaxF=3000.
 end if
 res@tiYAxisString = "Height (m)"

 res@pmLegendDisplayMode    = "Never"            ; turn on legend
 res@pmLegendSide           = "Bottom"               ; Change location of
 res@pmLegendParallelPosF   = .7             ; move units right
 res@pmLegendWidthF         = 0.08                ; Change width and
 res@pmLegendHeightF        = 0.18                ; height of legend.
 res@lgPerimOn              = False               ; turn off box around
 res@lgLabelFontHeightF     = .015                 ; label font height


 res@tmXBLabelsOn         = True            ; do not draw bottom labels
 res@tiXAxisOn            = True

 res@trXMinF=-20.
 res@trXMaxF=80.

 res@tiXAxisString = "<w !F8!q!F!!B!l!N!> (W m!S!-2!N!)"
 res@gsnXRefLine = 0.0
 plot(0)  = gsn_csm_xy (wks,thvd1D(0,:,:),zi,res) ; create plot
 res@xyLineColor = "black"
 res@xyMonoLineColor = True
 plot(1)  = gsn_csm_xy (wks,wthl1D_mf(0,:),zflx,res) ; create plot
 res@xyDashPattern = 1
 plot(3)  = gsn_csm_xy (wks,wthl1D_ed(0,:),zflx,res) ; create plot
 ;overlay(plot(0),plot(1))
 ;overlay(plot(0),plot(3))

 res@trXMinF=300.
 res@trXMaxF=303.

 res@tiXAxisString = "!F8!q!F!> (K)"
 res@xyDashPattern = 0
 plot(2)  = gsn_csm_xy (wks,thl1D(0,:),z,res) ; create plot
 print(thl1D(0,:) +" "+ z(:))

 resP                    = True                 ; modify the panel plot
  resP@gsnFrame           = False                ; don't advance panel plot
  resP@gsnPanelLabelBar   = False                 ; add common colorbar
 resP@gsnPanelMainString = casein     ; new way of setting main title (as of NCL V6.4.0)
;  if (casein.eq."BOMEX") then
;  gsn_panel(wks,plot((/0,2,4,6,8,14,15,13,17/)),(/3,4/),resP)           
;  gsn_panel(wks,plot((/0,6,12,13/)),(/1,4/),resP)           
;  else
;  gsn_panel(wks,plot((/0,2,6,8,14,15,13,17/)),(/2,4/),resP)           
;  end if
  gsn_panel(wks,plot((/0,2/)),(/1,2/),resP)           




