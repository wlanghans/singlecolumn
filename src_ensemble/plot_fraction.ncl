; ***********************************************
; xy_1.ncl
;
; Concepts illustrated:
;   - Drawing a black-and-white XY plot
;
;************************************************
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;************************************************
begin
;************************************************
; read in data
;************************************************
data1 = asciiread("data/dycor_fractionbetterthanuntuned.txt",-1,"float")
data2 = asciiread("data/stoch_fractionbetterthanuntuned.txt",-1,"float")
neval=dimsizes(data1)
data=new((/2,neval/),float)
data(0,:)=data1
data(1,:)=data2
x_vec=fspan(1,neval,neval)
;************************************************
; plotting parameters
;************************************************
 wks   = gsn_open_wks ("eps","fracbetterthanuntuned")               ; send graphics to PNG file

 res                  = True                     ; plot mods desired
 res@xyLineThicknessF = 4.0
 res@xyLineColors = (/"black","red"/)
 res@xyDashPatterns = 0
 res@trXMaxF = 100
 res@tiYAxisString    = "Fraction better than default"
 res@tiXAxisString    = "Number of evaluations"

 res@pmLegendDisplayMode    = "Always"            ; turn on legend
  res@pmLegendSide           = "Bottom"               ; Change location of
  res@pmLegendParallelPosF   = .8                  ; move units right
  res@pmLegendOrthogonalPosF = -0.5                ; move units down
  res@pmLegendWidthF         = 0.15                ; Change width and
  res@pmLegendHeightF        = 0.18                ; height of legend.
  res@lgPerimOn              = False               ; turn off box around
  res@lgLabelFontHeightF     = .03                 ; label font height
  res@xyExplicitLegendLabels = (/"Dycor","Stoch"/)         ; create explicit labels


 plot  = gsn_csm_xy (wks,x_vec,data,res) ; create plot
end

