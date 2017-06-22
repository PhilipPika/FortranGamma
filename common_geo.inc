        implicit real*8(a-h,o-z) 
        parameter (nsolid=13,ndiss=19,ncomp=32,nreac=37,nx=401)
        common/kinetics/sw01,sw02,sw03,sw04,sw05,sw06,sw07,sw08,sw09,sw1
     +0,sw11,sw12,sw13,sw14,sw15,sw16,sw17,sw18,sw19,xbiot,inita,nu,kmo2
     +,kmno3,kmmno2,kmfeoh3,kmso4,knit,kmnox,kfemno2,kfeo2,kh2so2,kh2smn
     +o2,kh2sfeoh3,kaom,kaomo2,kfeso2,kmnco3precip,KsMnCO3,kfesprecip,kf
     +esdiss,KsFeS,kfeco3precip,KsFeCO3,kpyr,kfess0,kcaldiss,kapa,po4_eq
     +,kdis,h2sstar,ch4eq,kgas,kdi,ksnh4,kspo4,ksfp,ksfe,kf1,kb1,kf2,kb2
     +,kf3,kb3,kf4,kb4,kf5,kb5,kf6,kb6,kf7,kb7,kf8,kb8,kspcal,keq1,keq2,
     +keq3,keq4
        real*8 sw01,sw02,sw03,sw04,sw05,sw06,sw07,sw08,sw09,sw10,sw11,sw
     +12,sw13,sw14,sw15,sw16,sw17,sw18,sw19,xbiot,inita,nu,kmo2,kmno3,km
     +mno2,kmfeoh3,kmso4,knit,kmnox,kfemno2,kfeo2,kh2so2,kh2smno2,kh2sfe
     +oh3,kaom,kaomo2,kfeso2,kmnco3precip,KsMnCO3,kfesprecip,kfesdiss,Ks
     +FeS,kfeco3precip,KsFeCO3,kpyr,kfess0,kcaldiss,kapa,po4_eq,kdis,h2s
     +star,ch4eq,kgas,kdi,ksnh4,kspo4,ksfp,ksfe,kf1,kb1,kf2,kb2,kf3,kb3,
     +kf4,kb4,kf5,kb5,kf6,kb6,kf7,kb7,kf8,kb8,kspcal,keq1,keq2,keq3,keq4
        common/physics/al,q0,w0,Db0,por0,area0,t_celsius,salin,delt,dept
     +hmax,alph0,inann,xirr,porX,attpor,kvalue,fvalue
        real*8 al,q0,w0,Db0,por0,area0,t_celsius,salin,delt,depthmax,alp
     +h0,inann(7,18),xirr,porX,attpor,kvalue(1,14),fvalue(1,14)
        common/physics2/iq,iw,iDb,ipor,igrid,iarea,ic,endt
        integer iq,iw,iDb,ipor,igrid,iarea,ic,endt
