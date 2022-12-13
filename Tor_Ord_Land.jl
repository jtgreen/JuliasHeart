using DifferentialEquations, Plots
using Sundials
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

#function output=model_ToRORd_Land(t,X,flag_ode, cellType, ICaL_Multiplier, ...
function model_ToRORd_Land(u,p,t,flag_ode, cellType, ICaL_Multiplier, INa_Multiplier, Ito_Multiplier, INaL_Multiplier, 
    IKr_Multiplier, IKs_Multiplier, IK1_Multiplier, IKb_Multiplier,INaCa_Multiplier,
    INaK_Multiplier, INab_Multiplier, ICab_Multiplier, IpCa_Multiplier, ICaCl_Multiplier, IClb_Multiplier, 
    Jrel_Multiplier,Jup_Multiplier, nao,cao,ko,ICaL_fractionSS,INaCa_fractionSS, stimAmp, stimDur, vcParameters, 
    apClamp, extraParams)

    celltype = cellType # %endo = 0, epi = 1, mid = 2

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %give names to the state vector values
    v=u[1]
    nai=u[2]
    nass=u[3]
    ki=u[4]
    kss=u[5]
    cai=u[6]
    cass=u[7]
    cansr=u[8]
    cajsr=u[9]
    m=u[10]
    hp=u[11]
    h=u[12]
    j=u[13]

    jp=u[14]
    mL=u[15]
    hL=u[16]
    hLp=u[17]
    a=u[18]
    iF=u[19]
    iS=u[20]
    ap=u[21]
    iFp=u[22]
    iSp=u[23]
    
    # % ical
    d=u[24]
    ff=u[25]
    fs=u[26]
    fcaf=u[27]
    fcas=u[28]
    jca=u[29]
    nca=u[30]
    nca_i=u[31]
    ffp=u[32]
    fcafp=u[33]
    # % end ical
    
    xs1=u[34]
    xs2=u[35]
    Jrel_np=u[36]
    CaMKt=u[37]
    
    # % new MM ICaL states
    ikr_c0 = u[38]
    ikr_c1 = u[39]
    ikr_c2 = u[40]
    ikr_o = u[41]
    ikr_i = u[42]
    Jrel_p = u[43]

    cli = 24   # % Intracellular Cl  [mM]
    clo = 150  # % Extracellular Cl  [mM]
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    XS=max(0,u[44])
    XW=max(0,u[45])
    Ca_TRPN=max(0,u[46])
    TmBlocked=u[47]
    ZETAS=u[48]
    ZETAW=u[49]
    Istim=u[50]


    
    # %% Land-Niederer model
    # %==========================================================================
    # % input coded here:
    mode="intact"
    lambda_=1
    lambda_rate=0

    # %EC parameters
    perm50=0.35
    TRPN_n=2
    koff=0.1
    dr=0.25
    wfrac=0.5
    TOT_A=25
    ktm_unblock= 0.021 
    beta_1=-2.4
    beta_0=2.3
    gamma=0.0085
    gamma_wu=0.615
    phi=2.23  

    if mode=="skinned"
        nperm=2.2
        ca50=2.5
        Tref=40.5
        nu=1
        mu=1
    else
        nperm=2.036 
        ca50=0.805
        Tref=120
        nu=7
        mu=3
    end
    
    k_ws=0.004
    k_uw=0.026

    lambda_min=0.87
    lambda_max=1.2

    # du = zeros((1,6))

    k_ws=k_ws*mu
    k_uw=k_uw*nu

    cdw=phi*k_uw*(1-dr)*(1-wfrac)/((1-dr)*wfrac)
    cds=phi*k_ws*(1-dr)*wfrac/dr
    k_wu=k_uw*(1/wfrac-1)-k_ws
    k_su=k_ws*(1/dr-1)*wfrac 
    A=(0.25*TOT_A)/((1-dr)*wfrac+dr)*(dr/0.25)

    # %XB model
    lambda0=min(lambda_max, lambda_)                                         # JTG array or not?
    Lfac=max(0, 1+beta_0*(lambda0+min(lambda_min,lambda0)-(1+lambda_min)))

    XU=(1-TmBlocked)-XW-XS # % unattached available xb = all - tm blocked - already prepowerstroke - already post-poststroke - no overlap
    xb_ws=k_ws*XW
    xb_uw=k_uw*XU 
    xb_wu=k_wu*XW
    xb_su=k_su*XS

    # gamma_rate=gamma*max((ZETAS>0).*ZETAS,(ZETAS<-1).*(-ZETAS-1)) 
    gamma_rate=gamma*max((ZETAS>0)*ZETAS,(ZETAS<-1)*(-ZETAS-1))
    xb_su_gamma=gamma_rate*XS
    gamma_rate_w=gamma_wu*abs(ZETAW) # % weak xbs don't like being strained
    xb_wu_gamma=gamma_rate_w*XW     

    XS=max(0,u[44])
    XW=max(0,u[45])
    Ca_TRPN=max(0,u[46])
    TmBlocked=u[47]
    ZETAS=u[48]
    ZETAW=u[49]


    #du[1]=xb_ws-xb_su-xb_su_gamma
    dXS=xb_ws-xb_su-xb_su_gamma
    #du[2]=xb_uw-xb_wu-xb_ws-xb_wu_gamma
    dXW=xb_uw-xb_wu-xb_ws-xb_wu_gamma

    ca50=ca50+beta_1*min(0.2,lambda_-1)
    #du[3]=koff*(((cai*1000)/ca50)^(TRPN_n*(1-Ca_TRPN)-Ca_TRPN)) # % untouched
    dCa_TRPN=koff*( (((cai*1000)/ca50)^TRPN_n) * (1-Ca_TRPN)-Ca_TRPN)

    XSSS=dr*0.5
    XWSS=(1-dr)*wfrac*0.5
    ktm_block=(ktm_unblock*(perm50^nperm)*0.5/(0.5-XSSS-XWSS))
    #du[4]=ktm_block*min(100,(Ca_TRPN^(-(nperm/2))))*XU-ktm_unblock*(Ca_TRPN^(nperm/2))*TmBlocked
    dTmBlock=ktm_block*min(100,(Ca_TRPN^(-(nperm/2)) ) )*XU-ktm_unblock*(Ca_TRPN^(nperm/2))*TmBlocked

    # %velocity dependence -- assumes distortion resets on W->S
    #du[5]=A*lambda_rate-cds*ZETAS # % - gamma_rate * ZETAS
    dZETAS=A*lambda_rate-cds*ZETAS # % - gamma_rate * ZETAS
    #du[6]=A*lambda_rate-cdw*ZETAW # % - gamma_rate_w * ZETAW
    dZETAW=A*lambda_rate-cdw*ZETAW # % - gamma_rate_w * ZETAW

    # % Active Force
    Ta=Lfac*(Tref/dr)*((ZETAS+1)*XS+(ZETAW)*XW)
    
    # %========================================================================

    # %physical constants
    R=8314.0
    T=310.0
    F=96485.0

    # %cell geometry
    L=0.01
    rad=0.0011
    vcell=1000*3.14*rad*rad*L
    Ageo=2*3.14*rad*rad+2*3.14*rad*L
    Acap=2*Ageo
    vmyo=0.68*vcell
    vnsr=0.0552*vcell
    vjsr=0.0048*vcell
    vss=0.02*vcell

    # %CaMK constants
    KmCaMK=0.15

    aCaMK=0.05
    bCaMK=0.00068
    CaMKo=0.05
    KmCaM=0.0015

    # %update CaMK
    CaMKb=CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass)
    CaMKa=CaMKb+CaMKt
    dCaMKt=aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # %reversal potentials
    ENa=(R*T/F)*log(nao/nai)
    EK=(R*T/F)*log(ko/ki)
    PKNa=0.01833
    EKs=(R*T/F)*log((ko+PKNa*nao)/(ki+PKNa*nai))

    # %convenient shorthand calculations
    vffrt=v*F*F/(R*T)
    vfrt=v*F/(R*T)
    frt=F/(R*T)




    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fINap=(1.0/(1.0+KmCaMK/CaMKa))
    fINaLp=(1.0/(1.0+KmCaMK/CaMKa))
    fItop=(1.0/(1.0+KmCaMK/CaMKa))
    fICaLp=(1.0/(1.0+KmCaMK/CaMKa))

    # %% INa
    INa, dm, dh, dhp, dj, djp = getINa_Grandi(v, m, h, hp, j, jp, fINap, ENa, INa_Multiplier)


    # %% INaL
    INaL,dmL,dhL,dhLp = getINaL_ORd2011(v, mL, hL, hLp, fINaLp, ENa, celltype, INaL_Multiplier)

    # %% ITo
    Ito,da,diF,diS,dap,diFp, diSp = getITo_ORd2011(v, a, iF, iS, ap, iFp, iSp, fItop, EK, celltype, Ito_Multiplier)


    # %% ICaL
    ICaL_ss,ICaNa_ss,ICaK_ss,ICaL_i,ICaNa_i,ICaK_i,dd,dff,dfs,dfcaf,dfcas,djca,dnca,dnca_i, dffp,dfcafp, PhiCaL_ss, PhiCaL_i, gammaCaoMyo, gammaCaiMyo = getICaL_ORd2011_jt(v, d,ff,fs,fcaf,fcas,jca,nca,nca_i,ffp,fcafp,
        fICaLp, cai, cass, cao, nai, nass, nao, ki, kss, ko, cli, clo, celltype, ICaL_fractionSS, ICaL_Multiplier)

    ICaL = ICaL_ss + ICaL_i
    ICaNa = ICaNa_ss + ICaNa_i
    ICaK = ICaK_ss + ICaK_i
    ICaL_tot = ICaL + ICaNa + ICaK

    # %% IKr
    IKr, dt_ikr_c0, dt_ikr_c1, dt_ikr_c2, dt_ikr_o, dt_ikr_i  = getIKr_ORd2011_MM(v,ikr_c0,ikr_c1, ikr_c2, ikr_o, ikr_i,
        ko, EK, celltype, IKr_Multiplier)

    # %% IKs
    IKs,dxs1, dxs2 = getIKs_ORd2011(v,xs1, xs2, cai,  EKs,  celltype, IKs_Multiplier)

    # %% IK1
    IK1 = getIK1_CRLP(v, ko, EK, celltype, IK1_Multiplier)

    # %% INaCa
    INaCa_i, INaCa_ss = getINaCa_ORd2011(v,F,R,T, nass, nai, nao, cass, cai, cao, celltype, INaCa_Multiplier, INaCa_fractionSS)

    # %% INaK
    INaK = getINaK_ORd2011(v, F, R, T, nai, nao, ki, ko, celltype, INaK_Multiplier)

    # %% Minor/background currents
    # %calculate IKb
    xkb=1.0/(1.0+exp(-(v-10.8968)/(23.9871)))
    GKb=0.0189*IKb_Multiplier
    if celltype==1
        GKb=GKb*0.6
    end
    IKb=GKb*xkb*(v-EK)

    # %calculate INab
    PNab=1.9239e-09*INab_Multiplier
    INab=PNab*vffrt*(nai*exp(vfrt)-nao)/(exp(vfrt)-1.0)

    # %calculate ICab
    PCab=5.9194e-08*ICab_Multiplier
    # % 
    ICab=PCab*4.0*vffrt*(gammaCaiMyo*cai*exp(2.0*vfrt)-gammaCaoMyo*cao)/(exp(2.0*vfrt)-1.0)

    # %calculate IpCa
    GpCa=5e-04*IpCa_Multiplier
    IpCa=GpCa*cai/(0.0005+cai)

    # %% Chloride
    # % I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current

    ecl = (R*T/F)*log(cli/clo);            # % [mV]

    Fjunc = 1
    Fsl = 1-Fjunc # % fraction in SS and in myoplasm - as per literature, I(Ca)Cl is in junctional subspace

    Fsl = 1-Fjunc # % fraction in SS and in myoplasm
    GClCa = ICaCl_Multiplier * 0.2843 # % [mS/uF]
    GClB = IClb_Multiplier * 1.98e-3 # % [mS/uF] %
    KdClCa = 0.1 # % [mM]

    I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/cass)*(v-ecl)
    I_ClCa_sl = Fsl*GClCa/(1+KdClCa/cai)*(v-ecl)

    I_ClCa = I_ClCa_junc+I_ClCa_sl
    I_Clbk = GClB*(v-ecl)

    # %% Calcium handling
    # %calculate ryanodione receptor calcium induced calcium release from the jsr
    fJrelp=(1.0/(1.0+KmCaMK/CaMKa))

    # %% Jrel
    Jrel, dJrel_np, dJrel_p = getJrel_ORd2011(Jrel_np, Jrel_p, ICaL_ss,cass, cajsr, fJrelp, celltype, Jrel_Multiplier)


    fJupp=(1.0/(1.0+KmCaMK/CaMKa))
    Jup, Jleak = getJup_ORd2011(cai, cansr, fJupp, celltype, Jup_Multiplier)

    # %calculate tranlocation flux
    Jtr=(cansr-cajsr)/60

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # %% calculate the stimulus current, Istim
    # amp=stimAmp
    # duration=stimDur
    # if t<=duration
    #     Istim=amp
    # else
    #     Istim=0.0
    # end
    
    # stimAmp=amp
    # if stimAmp != 0 
    #     Istim=stimAmp
    # else
    #     Istim=0.0
    # end

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # %update the membrane voltage

    # dv=-(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+ + I_ClCa+I_Clbk + Istim);
    # Is it missing a chanel or is this a typo?
    dv=-(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+I_ClCa+I_Clbk+Istim)


    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # %calculate diffusion fluxes
    JdiffNa=(nass-nai)/2.0
    JdiffK=(kss-ki)/2.0
    Jdiff=(cass-cai)/0.2

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # %calcium buffer constants 
    cmdnmax=0.05
    if celltype==1
        cmdnmax=cmdnmax*1.3
    end

    kmcmdn=0.00238
    trpnmax=0.07
    kmtrpn=0.0005
    BSRmax=0.047
    KmBSR = 0.00087
    BSLmax=1.124
    KmBSL = 0.0087
    csqnmax=10.0
    kmcsqn=0.8

    # %update intracellular concentrations, using buffers for cai, cass, cajsr
    dnai=-(ICaNa_i+INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo
    dnass=-(ICaNa_ss+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa

    dki=-(ICaK_i+Ito+IKr+IKs+IK1+IKb+Istim-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo
    dkss=-(ICaK_ss)*Acap/(F*vss)-JdiffK

    # % Bcai=1.0/(1.0+cmdnmax*kmcmdn/(kmcmdn+cai)^2.0+trpnmax*kmtrpn/(kmtrpn+cai)^2.0)
    # % dcai=Bcai*(-(ICaL_i + IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo);

    Bcai=1.0/(1.0+cmdnmax*kmcmdn/(kmcmdn+cai)^2.0)

    # JTG: dydt looks awkward here
    dcai=Bcai*(-(ICaL_i+IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo- dCa_TRPN *trpnmax) #dCa_TRPN was du[3]

    Bcass=1.0/(1.0+BSRmax*KmBSR/(KmBSR+cass)^2.0+BSLmax*KmBSL/(KmBSL+cass)^2.0)
    dcass=Bcass*(-(ICaL_ss-2.0*INaCa_ss)*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff)

    dcansr=Jup-Jtr*vjsr/vnsr

    Bcajsr=1.0/(1.0+csqnmax*kmcsqn/(kmcsqn+cajsr)^2.0)
    dcajsr=Bcajsr*(Jtr-Jrel)

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # %output the state vector when ode_flag==1, and the calculated currents and fluxes otherwise
    # JTG: converted a high-risk operator here, an ' at end of if/output; 
    # it flips the matrix around dimensions and takes the complex conjugate of the matrix
    # https://stackoverflow.com/questions/39885495/what-is-the-meaning-of-single-quote-in-matlab-and-how-to-change-it-to-python#:~:text=It%20flips%20the%20matrix%20around,a'%20in%20Python%20is%20
    if flag_ode==1

        # atleast_2d(a).T.conj()

        # output=[dv, dnai, dnass, dki, dkss, dcai, dcass, dcansr, dcajsr, dm, dhp, dh, dj, djp, dmL, dhL, dhLp, da, diF, diS, dap, diFp, diSp,
        #     dd, dff, dfs, dfcaf, dfcas, djca, dnca, dnca_i, dffp, dfcafp, dxs1, dxs2, dJrel_np, dCaMKt,
        #     dt_ikr_c0, dt_ikr_c1, dt_ikr_c2, dt_ikr_o, dt_ikr_i, dJrel_p, dydt(1), dydt(2), dydt(3), dydt(4), dydt(5), dydt(6)]'

        # return atleast_2d([dv, dnai, dnass, dki, dkss, dcai, dcass, dcansr, dcajsr, dm, dhp, dh, dj, djp, dmL, dhL, dhLp, da, diF, diS, dap, diFp, diSp,
        #     dd, dff, dfs, dfcaf, dfcas, djca, dnca, dnca_i, dffp, dfcafp, dxs1, dxs2, dJrel_np, dCaMKt,
        #     dt_ikr_c0, dt_ikr_c1, dt_ikr_c2, dt_ikr_o, dt_ikr_i, dJrel_p, dydt[1], dydt[2], dydt[3], dydt[4], dydt[5], dydt[6]]).T.conj()[0]

        du=zeros(50)
        du[1]=dv
        du[2]=dnai
        du[3]=dnass
        du[4]=dki
        du[5]=dkss
        du[6]=dcai
        du[7]=dcass
        du[8]=dcansr
        du[9]=dcajsr
        du[10]=dm
        du[11]=dhp
        du[12]=dh
        du[13]=dj
        du[14]=djp
        du[15]=dmL
        du[16]=dhL
        du[17]=dhLp
        du[18]=da
        du[19]=diF
        du[20]=diS
        du[21]=dap
        du[22]=diFp
        du[23]=diSp
        du[24]=dd
        du[25]=dff
        du[26]=dfs
        du[27]=dfcaf
        du[28]=dfcas
        du[29]=djca
        du[30]=dnca
        du[31]=dnca_i
        du[32]=dffp
        du[33]=dfcafp
        du[34]=dxs1
        du[35]=dxs2
        du[36]=dJrel_np
        du[37]=dCaMKt
        du[38]=dt_ikr_c0
        du[39]=dt_ikr_c1
        du[40]=dt_ikr_c2
        du[41]=dt_ikr_o
        du[42]=dt_ikr_i
        du[43]=dJrel_p
        du[44]=dXS
        du[45]=dXW
        du[46]=dCa_TRPN
        du[47]=dTmBlock
        du[48]=dZETAS
        du[49]=dZETAW
        #du[50]=u[50]
 
        return du
        # return [dv, dnai, dnass, dki, dkss, dcai, dcass, dcansr, dcajsr, dm, dhp, dh, dj, djp, dmL, dhL, dhLp, da, diF, diS, dap, diFp, diSp,
        #     dd, dff, dfs, dfcaf, dfcas, djca, dnca, dnca_i, dffp, dfcafp, dxs1, dxs2, dJrel_np, dCaMKt,
        #     dt_ikr_c0, dt_ikr_c1, dt_ikr_c2, dt_ikr_o, dt_ikr_i, dJrel_p, dXS, dXW, dCa_TRPN, dTmBlock, dZETAS, dZETAW]
    else
        # JTG: Confused, some were comma sepparated some were not!
        return [INa, INaL, Ito, ICaL, IKr, IKs, IK1, INaCa_i, INaCa_ss, INaK, IKb, INab, ICab, IpCa, Jdiff, JdiffNa, JdiffK, Jup, Jleak, Jtr, Jrel, CaMKa, Istim, fINap, 
            fINaLp, fICaLp, fJrelp, fJupp, cajsr, cansr, PhiCaL_ss, v, ICaL_i, I_ClCa, I_Clbk, ICaL_tot, Ta]
    end
end

# %% INa formulations
# function [INa, dm, dh, dhp, dj, djp] = getINa_Grandi(v, m, h, hp, j, jp, fINap, ENa, INa_Multiplier)
function getINa_Grandi(v, m, h, hp, j, jp, fINap, ENa, INa_Multiplier)
    # % The Grandi implementation updated with INa phosphorylation.
    # %% m gate
    mss = 1 / ((1 + exp( -(56.86 + v) / 9.03 ))^2)
    taum = 0.1292 * exp(-((v+45.79)/15.54)^2) + 0.06487 * exp(-((v-4.823)/51.12)^2)
    dm = (mss - m) / taum

    # %% h gate
    # JTG: This looks... hard to convert
    ah = (v >= -40) * (0) + (v < -40) * (0.057 * exp( -(v + 80) / 6.8 ))
    bh = (v >= -40) * (0.77 / (0.13*(1 + exp( -(v + 10.66) / 11.1 )))) + (v < -40) * ((2.7 * exp( 0.079 * v) + 3.1*10^5 * exp(0.3485 * v)))
    tauh = 1 / (ah + bh)
    hss = 1 / ((1 + exp( (v + 71.55)/7.43 ))^2)
    dh = (hss - h) / tauh
    
    # %% j gate
    aj = (v >= -40) * (0) +(v < -40) * (((-2.5428 * 10^4*exp(0.2444*v) - 6.948*10^-6 * exp(-0.04391*v)) * (v + 37.78)) / (1 + exp( 0.311 * (v + 79.23) )))
    bj = (v >= -40) * ((0.6 * exp( 0.057 * v)) / (1 + exp( -0.1 * (v + 32) ))) + (v < -40) * ((0.02424 * exp( -0.01052 * v )) / (1 + exp( -0.1378 * (v + 40.14) )))
    tauj = 1 / (aj + bj)
    jss = 1 / ((1 + exp( (v + 71.55)/7.43 ))^2)
    dj = (jss - j) / tauj

    # %% h phosphorylated
    hssp = 1 / ((1 + exp( (v + 71.55 + 6)/7.43 ))^2)
    dhp = (hssp - hp) / tauh
    
    # %% j phosphorylated
    taujp = 1.46 * tauj
    djp = (jss - jp) / taujp

    GNa = 11.7802
    INa=INa_Multiplier * GNa*(v-ENa)*m^3.0*((1.0-fINap)*h*j+fINap*hp*jp)

    return INa, dm, dh, dhp, dj, djp
end

# %% INaL
# function [INaL,dmL,dhL,dhLp] = getINaL_ORd2011(v, mL, hL, hLp, fINaLp, ENa, celltype, INaL_Multiplier)
function getINaL_ORd2011(v, mL, hL, hLp, fINaLp, ENa, celltype, INaL_Multiplier)

    # %calculate INaL
    mLss=1.0/(1.0+exp((-(v+42.85))/5.264))
    tm = 0.1292 * exp(-((v+45.79)/15.54)^2) + 0.06487 * exp(-((v-4.823)/51.12)^2)
    tmL=tm
    dmL=(mLss-mL)/tmL
    hLss=1.0/(1.0+exp((v+87.61)/7.488))
    thL=200.0
    dhL=(hLss-hL)/thL
    hLssp=1.0/(1.0+exp((v+93.81)/7.488))
    thLp=3.0*thL
    dhLp=(hLssp-hLp)/thLp
    GNaL=0.0279 * INaL_Multiplier
    if celltype==1
        GNaL=GNaL*0.6
    end

    INaL=GNaL*(v-ENa)*mL*((1.0-fINaLp)*hL+fINaLp*hLp)
    
    return INaL,dmL,dhL,dhLp
end

# %% ITo
# function [Ito,da,diF,diS,dap,diFp, diSp] = getITo_ORd2011(v, a, iF, iS, ap, iFp, iSp, fItop, EK, celltype, Ito_Multiplier)
function getITo_ORd2011(v, a, iF, iS, ap, iFp, iSp, fItop, EK, celltype, Ito_Multiplier)

    # %calculate Ito
    ass=1.0/(1.0+exp((-(v-14.34))/14.82))
    ta=1.0515/(1.0/(1.2089*(1.0+exp(-(v-18.4099)/29.3814)))+3.5/(1.0+exp((v+100.0)/29.3814)))
    da=(ass-a)/ta
    iss=1.0/(1.0+exp((v+43.94)/5.711))
    if celltype==1
        delta_epi=1.0-(0.95/(1.0+exp((v+70.0)/5.0)))
    else
        delta_epi=1.0
    end
    
    tiF=4.562+1/(0.3933*exp((-(v+100.0))/100.0)+0.08004*exp((v+50.0)/16.59))
    tiS=23.62+1/(0.001416*exp((-(v+96.52))/59.05)+1.780e-8*exp((v+114.1)/8.079))
    tiF=tiF*delta_epi
    tiS=tiS*delta_epi
    AiF=1.0/(1.0+exp((v-213.6)/151.2))
    AiS=1.0-AiF
    diF=(iss-iF)/tiF
    diS=(iss-iS)/tiS
    i=AiF*iF+AiS*iS
    assp=1.0/(1.0+exp((-(v-24.34))/14.82))
    dap=(assp-ap)/ta
    dti_develop=1.354+1.0e-4/(exp((v-167.4)/15.89)+exp(-(v-12.23)/0.2154))
    dti_recover=1.0-0.5/(1.0+exp((v+70.0)/20.0))
    tiFp=dti_develop*dti_recover*tiF
    tiSp=dti_develop*dti_recover*tiS
    diFp=(iss-iFp)/tiFp
    diSp=(iss-iSp)/tiSp
    ip=AiF*iFp+AiS*iSp
    Gto=0.16 * Ito_Multiplier
    if celltype==1
        Gto=Gto*2.0
    elseif celltype==2
        Gto=Gto*2.0
    end

    Ito=Gto*(v-EK)*((1.0-fItop)*a*i+fItop*ap*ip)
    
    return Ito,da,diF,diS,dap,diFp,diSp
end

# % a variant updated by jakub, using a changed activation curve
# % it computes both ICaL in subspace and myoplasm (_i)
# function [ICaL_ss,ICaNa_ss,ICaK_ss,ICaL_i,ICaNa_i,ICaK_i,dd,dff,dfs,dfcaf,dfcas,...
#     djca,dnca, dnca_i, dffp,dfcafp, PhiCaL_ss, PhiCaL_i, gammaCaoMyo, gammaCaiMyo] = getICaL_ORd2011_jt(v, d,ff,fs,fcaf,fcas,jca,nca, nca_i,ffp,fcafp,...
#     fICaLp, cai, cass, cao, nai, nass, nao, ki,kss,ko, cli, clo, celltype, ICaL_fractionSS, ICaL_PCaMultiplier)
function getICaL_ORd2011_jt(v, d,ff,fs,fcaf,fcas,jca,nca, nca_i, ffp, fcafp, fICaLp, cai, cass, cao, nai, nass, nao, ki, kss,ko, cli, clo, celltype, ICaL_fractionSS, ICaL_PCaMultiplier)


    # %physical constants
    R=8314.0
    T=310.0
    F=96485.0
    vffrt=v*F*F/(R*T)
    vfrt=v*F/(R*T)

    # %calculate ICaL, ICaNa, ICaK

    dss=1.0763*exp(-1.0070*exp(-0.0829*(v))) # % magyar
    if(v >31.4978) # % activation cannot be greater than 1
        dss = 1
    end

    td= 0.6+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)))

    dd=(dss-d)/td
    fss=1.0/(1.0+exp((v+19.58)/3.696))
    tff=7.0+1.0/(0.0045*exp(-(v+20.0)/10.0)+0.0045*exp((v+20.0)/10.0))
    tfs=1000.0+1.0/(0.000035*exp(-(v+5.0)/4.0)+0.000035*exp((v+5.0)/6.0))
    Aff=0.6
    Afs=1.0-Aff
    dff=(fss-ff)/tff
    dfs=(fss-fs)/tfs
    f=Aff*ff+Afs*fs
    fcass=fss
    tfcaf=7.0+1.0/(0.04*exp(-(v-4.0)/7.0)+0.04*exp((v-4.0)/7.0))
    tfcas=100.0+1.0/(0.00012*exp(-v/3.0)+0.00012*exp(v/7.0))

    Afcaf=0.3+0.6/(1.0+exp((v-10.0)/10.0))

    Afcas=1.0-Afcaf
    dfcaf=(fcass-fcaf)/tfcaf
    dfcas=(fcass-fcas)/tfcas
    fca=Afcaf*fcaf+Afcas*fcas

    tjca = 75
    jcass = 1.0/(1.0+exp((v+18.08)/(2.7916)))
    djca=(jcass-jca)/tjca
    tffp=2.5*tff
    dffp=(fss-ffp)/tffp
    fp=Aff*ffp+Afs*fs
    tfcafp=2.5*tfcaf
    dfcafp=(fcass-fcafp)/tfcafp
    fcap=Afcaf*fcafp+Afcas*fcas

    # %% SS nca
    Kmn=0.002
    k2n=500.0
    km2n=jca*1
    anca=1.0/(k2n/km2n+(1.0+Kmn/cass)^4.0)
    dnca=anca*k2n-nca*km2n

    # %% myoplasmic nca
    anca_i = 1.0/(k2n/km2n+(1.0+Kmn/cai)^4.0)
    dnca_i = anca_i*k2n-nca_i*km2n

    # %% SS driving force
    clo = 150
    cli = 2
    Io = 0.5*(nao + ko + clo + 4*cao)/1000 # % ionic strength outside. /1000 is for things being in micromolar
    Ii = 0.5*(nass + kss + cli + 4*cass)/1000 # % ionic strength outside. /1000 is for things being in micromolar
    
    # % The ionic strength is too high for basic DebHuc. We'll use Davies
    dielConstant = 74 # % water at 37
    temp = 310 # % body temp in kelvins.
    constA = 1.82*10^6*(dielConstant*temp)^(-1.5)

    gamma_cai = exp(-constA * 4 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii))
    gamma_cao = exp(-constA * 4 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io))
    gamma_nai = exp(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii))
    gamma_nao = exp(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io))
    gamma_ki = exp(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii))
    gamma_kao = exp(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io))


    PhiCaL_ss =  4.0*vffrt*(gamma_cai*cass*exp(2.0*vfrt)-gamma_cao*cao)/(exp(2.0*vfrt)-1.0)
    PhiCaNa_ss =  1.0*vffrt*(gamma_nai*nass*exp(1.0*vfrt)-gamma_nao*nao)/(exp(1.0*vfrt)-1.0)
    PhiCaK_ss =  1.0*vffrt*(gamma_ki*kss*exp(1.0*vfrt)-gamma_kao*ko)/(exp(1.0*vfrt)-1.0)

    # %% Myo driving force
    Io = 0.5*(nao + ko + clo + 4*cao)/1000 # % ionic strength outside. /1000 is for things being in micromolar
    Ii = 0.5*(nai + ki + cli + 4*cai)/1000 # % ionic strength outside. /1000 is for things being in micromolar
    
    # % The ionic strength is too high for basic DebHuc. We'll use Davies
    dielConstant = 74 # % water at 37�.
    temp = 310 # % body temp in kelvins.
    constA = 1.82*10^6*(dielConstant*temp)^(-1.5)

    gamma_cai = exp(-constA * 4 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii))
    gamma_cao = exp(-constA * 4 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io))
    gamma_nai = exp(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii))
    gamma_nao = exp(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io))
    gamma_ki = exp(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii))
    gamma_kao = exp(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io))

    gammaCaoMyo = gamma_cao
    gammaCaiMyo = gamma_cai

    PhiCaL_i = 4.0*vffrt*(gamma_cai*cai*exp(2.0*vfrt)-gamma_cao*cao)/(exp(2.0*vfrt)-1.0)
    PhiCaNa_i = 1.0*vffrt*(gamma_nai*nai*exp(1.0*vfrt)-gamma_nao*nao)/(exp(1.0*vfrt)-1.0)
    PhiCaK_i = 1.0*vffrt*(gamma_ki*ki*exp(1.0*vfrt)-gamma_kao*ko)/(exp(1.0*vfrt)-1.0)
    
    #%% The rest
    PCa=8.3757e-05 * ICaL_PCaMultiplier

    if celltype==1
        PCa=PCa*1.2
    elseif celltype==2
        PCa=PCa*1.8 # %2;
    end

    PCap=1.1*PCa
    PCaNa=0.00125*PCa
    PCaK=3.574e-4*PCa
    PCaNap=0.00125*PCap
    PCaKp=3.574e-4*PCap

    ICaL_ss=(1.0-fICaLp)*PCa*PhiCaL_ss*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCap*PhiCaL_ss*d*(fp*(1.0-nca)+jca*fcap*nca)
    ICaNa_ss=(1.0-fICaLp)*PCaNa*PhiCaNa_ss*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaNap*PhiCaNa_ss*d*(fp*(1.0-nca)+jca*fcap*nca)
    ICaK_ss=(1.0-fICaLp)*PCaK*PhiCaK_ss*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaKp*PhiCaK_ss*d*(fp*(1.0-nca)+jca*fcap*nca)

    ICaL_i=(1.0-fICaLp)*PCa*PhiCaL_i*d*(f*(1.0-nca_i)+jca*fca*nca_i)+fICaLp*PCap*PhiCaL_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i)
    ICaNa_i=(1.0-fICaLp)*PCaNa*PhiCaNa_i*d*(f*(1.0-nca_i)+jca*fca*nca_i)+fICaLp*PCaNap*PhiCaNa_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i)
    ICaK_i=(1.0-fICaLp)*PCaK*PhiCaK_i*d*(f*(1.0-nca_i)+jca*fca*nca_i)+fICaLp*PCaKp*PhiCaK_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i)


    # % And we weight ICaL (in ss) and ICaL_i
    ICaL_i = ICaL_i * (1-ICaL_fractionSS)
    ICaNa_i = ICaNa_i * (1-ICaL_fractionSS)
    ICaK_i = ICaK_i * (1-ICaL_fractionSS)
    ICaL_ss = ICaL_ss * ICaL_fractionSS
    ICaNa_ss = ICaNa_ss * ICaL_fractionSS
    ICaK_ss = ICaK_ss * ICaL_fractionSS

    return ICaL_ss,ICaNa_ss,ICaK_ss,ICaL_i,ICaNa_i,ICaK_i,dd,dff,dfs,dfcaf,dfcas, djca,dnca, dnca_i, dffp,dfcafp, PhiCaL_ss, PhiCaL_i, gammaCaoMyo, gammaCaiMyo
end
# % Variant based on Lu-Vandenberg
# function [IKr, dc0, dc1, dc2, do, di ] = getIKr_ORd2011_MM(V,c0,c1, c2, o, i,...
#     ko, EK, celltype, IKr_Multiplier)

function getIKr_ORd2011_MM(V,c0,c1, c2, o, i, ko, EK, celltype, IKr_Multiplier)

    # %physical constants
    R=8314.0
    T=310.0
    F=96485.0

    # % Extracting state vector
    # % c3 = y(1);
    # % c2 = y(2);
    # % c1 = y(3);
    # % o = y(4);
    # % i = y(5);
    b = 0 # % no channels blocked in via the mechanism of specific MM states
    vfrt = V*F/(R*T)

    # % transition rates
    # % from c0 to c1 in l-v model,
    alpha = 0.1161 * exp(0.2990 * vfrt)
    # % from c1 to c0 in l-v/
    beta =  0.2442 * exp(-1.604 * vfrt)

    # % from c1 to c2 in l-v/
    alpha1 = 1.25 * 0.1235 
    # % from c2 to c1 in l-v/
    beta1 =  0.1911

    # % from c2 to o/           c1 to o
    alpha2 =0.0578 * exp(0.9710 * vfrt)
    # % from o to c2/
    beta2 = 0.349e-3* exp(-1.062 * vfrt)

    # % from o to i
    alphai = 0.2533 * exp(0.5953 * vfrt)
    # % from i to o
    betai = 1.25* 0.0522 * exp(-0.8209 * vfrt)

    # % from c2 to i (from c1 in orig)
    alphac2ToI = 0.52e-4 * exp(1.525 * vfrt)
    # % from i to c2
    # % betaItoC2 = 0.85e-8 * exp(-1.842 * vfrt); %
    betaItoC2 = (beta2 * betai * alphac2ToI)/(alpha2 * alphai)
    # % transitions themselves
    # % for reason of backward compatibility of naming of an older version of a
    # % MM IKr, c3 in code is c0 in article diagram, c2 is c1, c1 is c2.

    dc0 = c1 * beta - c0 * alpha # % delta for c0
    dc1 = c0 * alpha + c2*beta1 - c1*(beta+alpha1) # % c1
    dc2 = c1 * alpha1 + o*beta2 + i*betaItoC2 - c2 * (beta1 + alpha2 + alphac2ToI) # % subtraction is into c2, to o, to i. % c2
    _do = c2 * alpha2 + i*betai - o*(beta2+alphai)
    di = c2*alphac2ToI + o*alphai - i*(betaItoC2 + betai)

    GKr = 0.0321 * sqrt(ko/5) * IKr_Multiplier # % 1st element compensates for change to ko (sqrt(5/5.4)* 0.0362)
    if celltype==1
        GKr=GKr*1.3
    elseif celltype==2
        GKr=GKr*0.8
    end

    IKr = GKr * o  * (V-EK)
    
    return IKr, dc0, dc1, dc2, _do, di
end

# function [IKs,dxs1, dxs2] = getIKs_ORd2011(v,xs1, xs2, cai, EKs, celltype, IKs_Multiplier)
function getIKs_ORd2011(v,xs1, xs2, cai, EKs, celltype, IKs_Multiplier)
    # %calculate IKs
    xs1ss=1.0/(1.0+exp((-(v+11.60))/8.932))
    txs1=817.3+1.0/(2.326e-4*exp((v+48.28)/17.80)+0.001292*exp((-(v+210.0))/230.0))
    dxs1=(xs1ss-xs1)/txs1
    xs2ss=xs1ss
    txs2=1.0/(0.01*exp((v-50.0)/20.0)+0.0193*exp((-(v+66.54))/31.0))
    dxs2=(xs2ss-xs2)/txs2
    KsCa=1.0+0.6/(1.0+(3.8e-5/cai)^1.4)
    GKs= 0.0011 * IKs_Multiplier
    
    if celltype==1
        GKs=GKs*1.4
    end
    
    IKs=GKs*KsCa*xs1*xs2*(v-EKs)
    
    return IKs,dxs1, dxs2
end

# function [IK1] = getIK1_CRLP(v,  ko , EK, celltype, IK1_Multiplier)
function getIK1_CRLP(v,  ko , EK, celltype, IK1_Multiplier)
    # % IK1
    aK1 = 4.094/(1+exp(0.1217*(v-EK-49.934)))
    bK1 = (15.72*exp(0.0674*(v-EK-3.257))+exp(0.0618*(v-EK-594.31)))/(1+exp(-0.1629*(v-EK+14.207)))
    K1ss = aK1/(aK1+bK1)

    GK1=IK1_Multiplier*0.6992; # %0.7266 # %* sqrt(5/5.4))
    
    if celltype==1
        GK1=GK1*1.2
    elseif celltype==2
        GK1=GK1*1.3
    end
    
    IK1=GK1*sqrt(ko/5)*K1ss*(v-EK)
    
    return IK1
end

# function [ INaCa_i, INaCa_ss] = getINaCa_ORd2011(v,F,R,T, nass, nai, nao, cass, cai, cao, celltype, INaCa_Multiplier, INaCa_fractionSS)
function getINaCa_ORd2011(v,F,R,T, nass, nai, nao, cass, cai, cao, celltype, INaCa_Multiplier, INaCa_fractionSS)
    
    zca = 2.0
    kna1=15.0
    kna2=5.0
    kna3=88.12
    kasymm=12.5
    wna=6.0e4
    wca=6.0e4
    wnaca=5.0e3
    kcaon=1.5e6
    kcaoff=5.0e3
    qna=0.5224
    qca=0.1670
    hca=exp((qca*v*F)/(R*T))
    hna=exp((qna*v*F)/(R*T))
    h1=1+nai/kna3*(1+hna)
    h2=(nai*hna)/(kna3*h1)
    h3=1.0/h1
    h4=1.0+nai/kna1*(1+nai/kna2)
    h5=nai*nai/(h4*kna1*kna2)
    h6=1.0/h4
    h7=1.0+nao/kna3*(1.0+1.0/hna)
    h8=nao/(kna3*hna*h7)
    h9=1.0/h7
    h10=kasymm+1.0+nao/kna1*(1.0+nao/kna2)
    h11=nao*nao/(h10*kna1*kna2)
    h12=1.0/h10
    k1=h12*cao*kcaon
    k2=kcaoff
    k3p=h9*wca
    k3pp=h8*wnaca
    k3=k3p+k3pp
    k4p=h3*wca/hca
    k4pp=h2*wnaca
    k4=k4p+k4pp
    k5=kcaoff
    k6=h6*cai*kcaon
    k7=h5*h2*wna
    k8=h8*h11*wna
    x1=k2*k4*(k7+k6)+k5*k7*(k2+k3)
    x2=k1*k7*(k4+k5)+k4*k6*(k1+k8)
    x3=k1*k3*(k7+k6)+k8*k6*(k2+k3)
    x4=k2*k8*(k4+k5)+k3*k5*(k1+k8)
    E1=x1/(x1+x2+x3+x4)
    E2=x2/(x1+x2+x3+x4)
    E3=x3/(x1+x2+x3+x4)
    E4=x4/(x1+x2+x3+x4)
    KmCaAct=150.0e-6
    allo=1.0/(1.0+(KmCaAct/cai)^2.0)
    zna=1.0
    JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp
    JncxCa=E2*k2-E1*k1
    Gncx= 0.0034* INaCa_Multiplier
    
    if celltype==1
        Gncx=Gncx*1.1
    elseif celltype==2
        Gncx=Gncx*1.4
    end
    
    INaCa_i=(1-INaCa_fractionSS)*Gncx*allo*(zna*JncxNa+zca*JncxCa)

    # %calculate INaCa_ss
    h1=1+nass/kna3*(1+hna)
    h2=(nass*hna)/(kna3*h1)
    h3=1.0/h1
    h4=1.0+nass/kna1*(1+nass/kna2)
    h5=nass*nass/(h4*kna1*kna2)
    h6=1.0/h4
    h7=1.0+nao/kna3*(1.0+1.0/hna)
    h8=nao/(kna3*hna*h7)
    h9=1.0/h7
    h10=kasymm+1.0+nao/kna1*(1+nao/kna2)
    h11=nao*nao/(h10*kna1*kna2)
    h12=1.0/h10
    k1=h12*cao*kcaon
    k2=kcaoff
    k3p=h9*wca
    k3pp=h8*wnaca
    k3=k3p+k3pp
    k4p=h3*wca/hca
    k4pp=h2*wnaca
    k4=k4p+k4pp
    k5=kcaoff
    k6=h6*cass*kcaon
    k7=h5*h2*wna
    k8=h8*h11*wna
    x1=k2*k4*(k7+k6)+k5*k7*(k2+k3)
    x2=k1*k7*(k4+k5)+k4*k6*(k1+k8)
    x3=k1*k3*(k7+k6)+k8*k6*(k2+k3)
    x4=k2*k8*(k4+k5)+k3*k5*(k1+k8)
    E1=x1/(x1+x2+x3+x4)
    E2=x2/(x1+x2+x3+x4)
    E3=x3/(x1+x2+x3+x4)
    E4=x4/(x1+x2+x3+x4)
    KmCaAct=150.0e-6 
    allo=1.0/(1.0+(KmCaAct/cass)^2.0)
    JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp
    JncxCa=E2*k2-E1*k1
    INaCa_ss=INaCa_fractionSS*Gncx*allo*(zna*JncxNa+zca*JncxCa)
    
    return INaCa_i, INaCa_ss
end

# function INaK = getINaK_ORd2011(v, F, R, T, nai, nao, ki, ko, celltype, INaK_Multiplier)
function getINaK_ORd2011(v, F, R, T, nai, nao, ki, ko, celltype, INaK_Multiplier)

    # %calculate INaK
    zna=1.0
    k1p=949.5
    k1m=182.4
    k2p=687.2
    k2m=39.4
    k3p=1899.0
    k3m=79300.0
    k4p=639.0
    k4m=40.0
    Knai0=9.073
    Knao0=27.78
    delta=-0.1550
    Knai=Knai0*exp((delta*v*F)/(3.0*R*T))
    Knao=Knao0*exp(((1.0-delta)*v*F)/(3.0*R*T))
    Kki=0.5
    Kko=0.3582
    MgADP=0.05
    MgATP=9.8
    Kmgatp=1.698e-7
    H=1.0e-7
    eP=4.2
    Khp=1.698e-7
    Knap=224.0
    Kxkur=292.0
    P=eP/(1.0+H/Khp+nai/Knap+ki/Kxkur)
    a1=(k1p*(nai/Knai)^3.0)/((1.0+nai/Knai)^3.0+(1.0+ki/Kki)^2.0-1.0)
    b1=k1m*MgADP
    a2=k2p
    b2=(k2m*(nao/Knao)^3.0)/((1.0+nao/Knao)^3.0+(1.0+ko/Kko)^2.0-1.0)
    a3=(k3p*(ko/Kko)^2.0)/((1.0+nao/Knao)^3.0+(1.0+ko/Kko)^2.0-1.0)
    b3=(k3m*P*H)/(1.0+MgATP/Kmgatp)
    a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp)
    b4=(k4m*(ki/Kki)^2.0)/((1.0+nai/Knai)^3.0+(1.0+ki/Kki)^2.0-1.0)
    x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2
    x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4
    x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1
    x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1
    E1=x1/(x1+x2+x3+x4)
    E2=x2/(x1+x2+x3+x4)
    E3=x3/(x1+x2+x3+x4)
    E4=x4/(x1+x2+x3+x4)
    zk=1.0
    JnakNa=3.0*(E1*a3-E2*b3)
    JnakK=2.0*(E4*b1-E3*a1)
    Pnak= 15.4509 * INaK_Multiplier

    if celltype==1
        Pnak=Pnak*0.9
    elseif celltype==2
        Pnak=Pnak*0.7
    end
    
    INaK=Pnak*(zna*JnakNa+zk*JnakK)

    return INaK
end

# %% Jrel
# function [Jrel, dJrelnp, dJrelp] = getJrel_ORd2011(Jrelnp, Jrelp, ICaL, cass, cajsr, fJrelp, celltype, Jrel_Multiplier)
function getJrel_ORd2011(Jrelnp, Jrelp, ICaL, cass, cajsr, fJrelp, celltype, Jrel_Multiplier)
    jsrMidpoint = 1.7

    bt=4.75
    a_rel=0.5*bt
    Jrel_inf=a_rel*(-ICaL)/(1.0+(jsrMidpoint/cajsr)^8.0)
    if celltype==2
        Jrel_inf=Jrel_inf*1.7
    end
    
    tau_rel=bt/(1.0+0.0123/cajsr)

    if tau_rel<0.001
        tau_rel=0.001
    end

    dJrelnp=(Jrel_inf-Jrelnp)/tau_rel
    btp=1.25*bt
    a_relp=0.5*btp
    Jrel_infp=a_relp*(-ICaL)/(1.0+(jsrMidpoint/cajsr)^8.0)

    if celltype==2
        Jrel_infp=Jrel_infp*1.7
    end
    
    tau_relp=btp/(1.0+0.0123/cajsr)

    if tau_relp<0.001
        tau_relp=0.001
    end

    dJrelp=(Jrel_infp-Jrelp)/tau_relp

    Jrel=Jrel_Multiplier * 1.5378 * ((1.0-fJrelp)*Jrelnp+fJrelp*Jrelp)

    return Jrel, dJrelnp, dJrelp
end

# %% Jup
# function [Jup, Jleak] = getJup_ORd2011(cai, cansr, fJupp, celltype, Jup_Multiplier)
function getJup_ORd2011(cai, cansr, fJupp, celltype, Jup_Multiplier)
    # %calculate serca pump, ca uptake flux
    # % camkFactor = 2.4;
    # % gjup = 0.00696;
    # % Jupnp=Jup_Multiplier * gjup*cai/(cai+0.001);
    # % Jupp=Jup_Multiplier * camkFactor*gjup*cai/(cai + 8.2500e-04);
    # % if celltype==1
    # %     Jupnp=Jupnp*1.3;
    # %     Jupp=Jupp*1.3;
    # % end
    # % 
    # % 
    # % Jleak=Jup_Multiplier * 0.00629 * cansr/15.0;
    # % Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;

    # %calculate serca pump, ca uptake flux
    Jupnp=Jup_Multiplier * 0.005425*cai/(cai+0.00092)
    Jupp=Jup_Multiplier * 2.75*0.005425*cai/(cai+0.00092-0.00017)
    if celltype==1
        Jupnp=Jupnp*1.3
        Jupp=Jupp*1.3
    end

    Jleak=Jup_Multiplier* 0.0048825*cansr/15.0
    Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak

    return Jup, Jleak
end

function getStartingState(varargin)
    # Added a 0 starting state for Istim to end of each: JTG
    if varargin == "m_endo"
            X0 = [-88.6369922306458,11.8973412949238,11.8976610470850,141.234464714982,141.234423402713,7.26747296460659e-05,6.33786975780735e-05,1.53265306371970,1.53394579180493,0.000828007761976018,0.666527193684116,0.826020806005678,0.826055985895856,0.825850881115628,0.000166868626513013,0.522830604669169,0.285969584294187,0.000959137028030184,0.999601150012565,0.593401639836100,0.000488696137242056,0.999601147267179,0.654668660159696,9.50007519781516e-32,0.999999992317577,0.939258048397962,0.999999992317557,0.999898379647465,0.999978251560040,0.000444816183420527,0.000755072490632667,0.999999992318446,0.999999992318445,0.242404683449520,0.000179537726989804,-6.88308558109975e-25,0.0111749845355653,0.998036620213316,0.000858801779013532,0.000709744678350176,0.000381261722195702,1.35711566929992e-05,2.30252452954649e-23,0.000156194131630688,0.000235128869057516,0.00807763107171867,0.999373435811656,0.0,0.0,0.0]
    elseif varargin == "m_epi"
            X0 = [-89.0462806262884,12.7218980311997,12.7222039977392,142.248960281735,142.248911688304,6.54105789316085e-05,5.68443136844764e-05,1.80911728399381,1.80970235621251,0.000758182108180449,0.679839847935577,0.834150231581688,0.834188252920967,0.834081731522592,0.000154387698861246,0.538295069820379,0.302769394159465,0.000933035060391086,0.999628705730844,0.999626204093615,0.000475390662092180,0.999628705544664,0.999628513430851,1.74213411952898e-37,0.999999993122906,0.947952168523141,0.999999993122889,0.999932686646139,0.999982915381882,0.000291544679470133,0.000502604507932921,0.999999993124187,0.999999993123756,0.228815500940270,0.000171497784228012,-1.13118992668881e-26,0.0129505221481656,0.998194356754674,0.000834232097912889,0.000683865770895308,0.000277878501096440,9.66775862738005e-06,8.16930403133409e-24,0.000125999575445634,0.000189952183128226,0.00655149435193622,0.999493972060333,0.0,0.0,0.0]
    elseif varargin == "m_mid"
            X0 = [-89.5379994049964,14.9292004720038,14.9296673679334,144.844718688810,144.844658476157,7.50228807455408e-05,6.10763598140135e-05,1.79043480744558,1.79484249993962,0.000681936485046493,0.695380653101535,0.843488797335149,0.843520761455969,0.843226224403045,0.000140621109700401,0.545314876174586,0.292496735833565,0.000902612655601118,0.999659345906191,0.563119679366890,0.000459883274920751,0.999659343029625,0.623696443871387,-1.31418873360667e-33,0.999999993979673,0.920408593154793,0.999999993979652,0.999761950174748,0.999962530196306,0.000385359469667100,0.000853529194511867,0.999999993978835,0.999999993980401,0.266415111925392,0.000162310655612839,1.20976169203982e-24,0.0178243652102213,0.997971986641796,0.000805399061926759,0.000678179976274546,0.000526536308931670,1.78956481798154e-05,7.05916237956270e-23,0.000167065448972034,0.000250679417924283,0.00860262481216925,0.999331445205816,0.0,0.0,0.0]
    end
    return X0
end

function modelRunner( X0, options, parameters, beats, ignoreFirst)
    # % Parameters are set here
    # % defaults which may be overwritten

    haskey(parameters, "cellType")==true ? cellType = parameters["cellType"] : cellType = 0

    # % if the number of simulated beat is to be printed out.
    haskey(parameters, "verbose")==true ? verbose = parameters["verbose"] : verbose = true
    haskey(parameters, "nao")==true ? nao = parameters["nao"] : nao = 140    
    haskey(parameters, "cao")==true ? cao = parameters["cao"] : cao = 1.8
    haskey(parameters, "ko")==true ? ko = parameters["ko"] : ko = 5
    haskey(parameters, "ICaL_fractionSS")==true ? ICaL_fractionSS = parameters["ICaL_fractionSS"] : ICaL_fractionSS = 0.8
    haskey(parameters, "INaCa_fractionSS")==true ? INaCa_fractionSS = parameters["INaCa_fractionSS"] : INaCa_fractionSS = 0.35 
    haskey(parameters, "INa_Multiplier")==true ? INa_Multiplier = parameters["INa_Multiplier"] : INa_Multiplier = 1.0 
    haskey(parameters, "ICaL_Multiplier")==true ? ICaL_Multiplier = parameters["ICaL_Multiplier"] : ICaL_Multiplier = 1.0
    haskey(parameters, "Ito_Multiplier")==true ? Ito_Multiplier = parameters["Ito_Multiplier"] : Ito_Multiplier = 1.0
    haskey(parameters, "INaL_Multiplier")==true ? INaL_Multiplier = parameters["INaL_Multiplier"] : INaL_Multiplier = 1.0
    haskey(parameters, "IKr_Multiplier")==true ? IKr_Multiplier = parameters["IKr_Multiplier"] : IKr_Multiplier = 1.0
    haskey(parameters, "IKs_Multiplier")==true ? IKs_Multiplier = parameters["IKs_Multiplier"] : IKs_Multiplier = 1.0
    haskey(parameters, "IK1_Multiplier")==true ? IK1_Multiplier = parameters["IK1_Multiplier"] : IK1_Multiplier = 1.0
    haskey(parameters, "IKb_Multiplier")==true ? IKb_Multiplier = parameters["IKb_Multiplier"] : IKb_Multiplier = 1.0  
    haskey(parameters, "INaCa_Multiplier")==true ? INaCa_Multiplier = parameters["INaCa_Multiplier"] : INaCa_Multiplier = 1.0
    haskey(parameters, "INaK_Multiplier")==true ? INaK_Multiplier = parameters["INaK_Multiplier"] : INaK_Multiplier = 1.0
    haskey(parameters, "INab_Multiplier")==true ? INab_Multiplier = parameters["INab_Multiplier"] : INab_Multiplier = 1.0
    haskey(parameters, "ICab_Multiplier")==true ? ICab_Multiplier = parameters["ICab_Multiplier"] : ICab_Multiplier = 1.0
    haskey(parameters, "IpCa_Multiplier")==true ? IpCa_Multiplier = parameters["IpCa_Multiplier"] : IpCa_Multiplier = 1.0
    haskey(parameters, "ICaCl_Multiplier")==true ? ICaCl_Multiplier = parameters["ICaCl_Multiplier"] : ICaCl_Multiplier = 1.0
    haskey(parameters, "IClb_Multiplier")==true ? IClb_Multiplier = parameters["IClb_Multiplier"] : IClb_Multiplier = 1.0
    haskey(parameters, "Jrel_Multiplier")==true ? Jrel_Multiplier = parameters["Jrel_Multiplier"] : Jrel_Multiplier = 1.0
    haskey(parameters, "Jup_Multiplier")==true ? Jup_Multiplier = parameters["Jup_Multiplier"] : Jup_Multiplier = 1.0

    extraParams = []
    # if 'extraParams' in parameters: extraParams = parameters['extraParams']
    vcParameters = []
    # if 'vcParameters' in parameters: vcParameters = parameters['vcParameters']
    apClamp = []
    # if 'apClamp' in parameters: apClamp = parameters['apClamp']

    haskey(parameters, "stimAmp")==true ? stimAmp = parameters["stimAmp"] : stimAmp = -53.0
    haskey(parameters, "stimDur")==true ? stimDur = parameters["stimDur"] : stimDur = 1.0

    #if 'model' not in parameters: parameters['model'] = "@model11"

    haskey(parameters, "bcl")==true ? bcl = parameters["bcl"] : bcl = 500

    X = []
    flag_ode = 1
    sol = []
    for n in range(1,beats)

        prob = ODEProblem((u,p,t) -> model_ToRORd_Land(u,p,t, flag_ode, cellType, ICaL_Multiplier,
                 INa_Multiplier, Ito_Multiplier, INaL_Multiplier, IKr_Multiplier, IKs_Multiplier, IK1_Multiplier, IKb_Multiplier,INaCa_Multiplier,
                 INaK_Multiplier, INab_Multiplier, ICab_Multiplier, IpCa_Multiplier, ICaCl_Multiplier, IClb_Multiplier, Jrel_Multiplier,Jup_Multiplier,nao,cao,ko,ICaL_fractionSS,INaCa_fractionSS, stimAmp, stimDur, vcParameters, apClamp, extraParams), 
                 X0, (0.0, bcl))

        #sol = solve(prob, alg_hints=[:stiff], reltol=1e-9, abstol=1e-9, progress = true,
        #progress_steps = 1, ) #dense=true);
        if haskey(parameters, "callback")
            # This would break if tstops wasn't passed
            # Make this extensible to other callbacks
            sol = solve(prob, CVODE_BDF(), callback=parameters["callback"], reltol=1e-8, abstol=1e-8, tstops=parameters["tstops"]) #dtmax=0.5)
        else
            sol = solve(prob, CVODE_BDF(), reltol=1e-8, abstol=1e-8) #dtmax=0.5)
        end

        # if verbose:
        #     print('Beat = ', str(n))
        
        # Email Elissa about this
        # if ((exist('ode15sTimed')) == 2) # % if timed version provided, it is preferred
        #     [time{n}, X{n}]=ode15sTimed(parameters.model,[0 CL],X0,options,1,  cellType, ICaL_Multiplier, ...
        #         INa_Multiplier, Ito_Multiplier, INaL_Multiplier, IKr_Multiplier, IKs_Multiplier, IK1_Multiplier, IKb_Multiplier,INaCa_Multiplier,...
        #         INaK_Multiplier, INab_Multiplier, ICab_Multiplier, IpCa_Multiplier, ICaCl_Multiplier, IClb_Multiplier, Jrel_Multiplier,Jup_Multiplier,nao,cao,ko,ICaL_fractionSS,INaCa_fractionSS, stimAmp, stimDur, vcParameters, apClamp, extraParams);
            
        # X.append(
        #     odeint(parameters['model'], X0, t, (flag_ode, cellType, ICaL_Multiplier,
        #         INa_Multiplier, Ito_Multiplier, INaL_Multiplier, IKr_Multiplier, IKs_Multiplier, IK1_Multiplier, IKb_Multiplier,INaCa_Multiplier,
        #         INaK_Multiplier, INab_Multiplier, ICab_Multiplier, IpCa_Multiplier, ICaCl_Multiplier, IClb_Multiplier, Jrel_Multiplier,Jup_Multiplier,nao,cao,ko,ICaL_fractionSS,INaCa_fractionSS, stimAmp, stimDur, vcParameters, apClamp, extraParams))
        #     )
        
        # X0=X{n}(size(X{n},1),:);
        # JTG: Interpretation -- From the last beat, grab last row

        # %n %output beat number to the screen to monitor runtime progress

    end
    return sol
end